"""Exact local second order of the FLEX general path (design
2026-09-08, rev 10): the on-site interaction compiler, the off-site
density-slot arrays, the sparse factors and the batched kernel.

This module carries the on-site compiler, the off-site density-slot
arrays (:func:`build_offsite`), the sparse factor pack
(:func:`build_factors`) and the batched kernel (:func:`accumulate_batch`,
plus the test helper :func:`dense_w2`). The compiler turns the reader's
on-site interaction rows into the
spin-resolved antisymmetrised tensor ``Gamma`` of the normal-ordered Hamiltonian

    H_int = 1/2 sum_{pqrs} V_{pq,rs} c^dag_p c^dag_q c_s c_r ,
    Gamma[p, q, r, s] = V[p, q, r, s] - V[p, q, s, r],

with ``p = (orbital a, spin s)`` collapsed to the generalised index
``s * norb + a`` (the spin-block order of the Hartree-Fock kernel's
``(2*norb, 2*norb)`` accumulator). The per-type records below are the
normal-ordered form of the SAME file operators the Hartree-Fock kernel
implements (:mod:`hwave.solver.hartree_fock`), which is the AUTHORITY on
their conventions: Hund carries a minus sign, Exchange contributes both
spin directions per row, PairHop is recorded in one orientation (its
reverse hop is the transposed row, which the reversal closure conjugates),
PairLift contributes a record and its Hermitian conjugate per row, and a
bond whose transposed row repeats the same operator is half-weighted per
ordered entry (see :data:`_MIRRORED_TYPES`). The whole set is pinned
against ``accumulate_hf`` by
``tests/test_second_order_compiler.py::TestContracts``.
"""
import itertools
import logging
from dataclasses import dataclass as _dataclass

import numpy as np

from . import backend as _bk
from .hartree_fock import NonFiniteError

logger = logging.getLogger(__name__)

TWO_ORBITAL_TYPES = ("CoulombInter", "Hund", "Ising", "Exchange", "PairHop", "PairLift")
OFFSITE_DENSITY_TYPES = ("CoulombInter", "Hund", "Ising")
UP, DN = 0, 1


class DegenerateRowError(ValueError):
    """An on-site same-orbital row of a two-orbital type (spec D7)."""


_D7_REMEDY = {
    "CoulombInter": "v n_a n_a = 2 v n_{a up} n_{a dn} + v n_a: declare CoulombIntra with 2 v "
                    "(the one-body part v n_a is a level shift: add it to the Transfer file or "
                    "absorb it in the chemical potential)",
    "Hund": "-v (n_{a up}^2 + n_{a dn}^2) = -v n_a: a one-body level shift with no two-body "
            "content -- remove the row or add the shift to the Transfer file",
    "Ising": "v (n_{a up} - n_{a dn})^2 = v n_a - 2 v n_{a up} n_{a dn}: declare CoulombIntra "
             "with -2 v and add the level shift v n_a to the Transfer file (or absorb it in the "
             "chemical potential)",
    "Exchange": "v c+_{a up} c_{a up} c+_{a dn} c_{a dn} = v n_{a up} n_{a dn} (both directions: "
                "2 v): declare CoulombIntra with 2 v",
    "PairHop": "the row and its Hermitian-closed partner give (v + v*) n_{a up} n_{a dn}: declare "
               "CoulombIntra with 2 Re v (an imaginary part has no physical content here)",
    "PairLift": "identically zero (c+_{a up} twice): remove the row",
}


def gen_index(a, s, norb):
    """Generalised on-site index (orbital a, spin s) -> s*norb + a."""
    return s * norb + a


def _orbital_pair(key):
    """The orbital pair ``(a, b)`` of one table key, or ``None`` if the key
    is not an on-site row. Accepts the reader's ``((rx, ry, rz), (a, b))``
    keys and the already-reduced ``(a, b)`` keys."""
    first, second = key
    if isinstance(first, (tuple, list)):
        if tuple(int(x) for x in first) != (0, 0, 0):
            return None
        return (int(second[0]), int(second[1]))
    return (int(first), int(second))


def _raw_onsite_rows(onsite_tbl):
    """``{type: {(a, b): complex}}`` -- the declared on-site rows summed per
    ``(a, b)``, without the reversal closure."""
    out = {}
    for itype, tbl in onsite_tbl.items():
        raw = {}
        for key, v in tbl.items():
            ab = _orbital_pair(key)
            if ab is None:
                continue
            raw[ab] = raw.get(ab, 0.0) + complex(v)
        out[itype] = raw
    return out


def close_onsite_rows(onsite_tbl):
    """The kernel's reversal closure at R = 0: v_closed(a, b) = (v_ab + conj(v_ba)) / 2
    per type (``_reverse_closed`` of hwave.solver.hartree_fock restricted to the
    on-site cell). Returns {type: {(a, b): complex}}.

    Like the kernel's array-valued closure, the result carries BOTH ordered
    entries of every declared bond: a lone row ``v_ab`` closes to
    ``v_ab / 2`` at ``(a, b)`` and ``conj(v_ab) / 2`` at ``(b, a)``."""
    out = {}
    for itype, raw in _raw_onsite_rows(onsite_tbl).items():
        keys = set(raw) | set((q, p) for (p, q) in raw)
        out[itype] = {(a, b): 0.5 * (raw.get((a, b), 0.0)
                                     + np.conj(raw.get((b, a), 0.0)))
                      for (a, b) in keys}
    return out


# Types whose transposed row (b, a) carries the SAME operator as (a, b).
# The reversal closure always emits both ordered entries of a bond
# (:func:`close_onsite_rows`), so for these types each ordered entry is
# half of the bond: the declared matrix means H = 1/2 sum_{a b} V_ab ...
# This is the Hartree-Fock kernel's convention -- its Hartree and Fock
# einsums contract the closed table over the ordered pair (a, b) once,
# so a bond declared in both directions with V couples with V, not 2 V.
# PairHop is NOT in this set: its (b, a) row is the reverse pair hop
# (the Hermitian conjugate), a different operator, and the kernel gives
# it full weight. CoulombIntra has no transposed partner at all.
_MIRRORED_TYPES = ("CoulombInter", "Hund", "Ising", "Exchange", "PairLift")


def _records(itype, a, b, v, norb):
    """Normal-ordered records (p, q, s, r, x) of the monomial(s) of one
    closed row (spec 2.2 table): x c+_p c+_q c_s c_r."""
    g = lambda o, sp: gen_index(o, sp, norb)
    w = 0.5 * v if itype in _MIRRORED_TYPES else v
    if itype == "CoulombIntra":
        return [(g(a, UP), g(a, DN), g(a, DN), g(a, UP), w)]
    if itype == "CoulombInter":
        return [(g(a, s1), g(b, s2), g(b, s2), g(a, s1), w)
                for s1, s2 in itertools.product((UP, DN), repeat=2)]
    if itype == "Hund":
        return [(g(a, s1), g(b, s1), g(b, s1), g(a, s1), -w) for s1 in (UP, DN)]
    if itype == "Ising":
        return [(g(a, s1), g(b, s2), g(b, s2), g(a, s1), w if s1 == s2 else -w)
                for s1, s2 in itertools.product((UP, DN), repeat=2)]
    if itype == "Exchange":
        return [(g(a, s1), g(b, 1 - s1), g(b, s1), g(a, 1 - s1), -w) for s1 in (UP, DN)]
    if itype == "PairHop":
        # w c+_{a up} c+_{a dn} c_{b dn} c_{b up}: written on the DIRECT
        # placement (r, s) = (b up, b dn) so that V alone -- not only the
        # antisymmetrised Gamma -- carries the kernel's PairHop term.
        return [(g(a, UP), g(a, DN), g(b, DN), g(b, UP), w)]
    if itype == "PairLift":
        # -w c+_{a up} c+_{b up} c_{a dn} c_{b dn} and its Hermitian
        # conjugate, both written on the DIRECT placement (see PairHop).
        return [(g(a, UP), g(b, UP), g(b, DN), g(a, DN), w),
                (g(b, DN), g(a, DN), g(a, UP), g(b, UP), w)]
    raise ValueError("second_order: unknown interaction type {!r}".format(itype))


def _check_norb(norb):
    """A usable orbital count. Checked ONCE at the entry of
    :func:`compile_onsite_v`, not per declared type: an EMPTY table has no
    type to iterate, so a per-type check would let ``norb = 0`` through and
    build a degenerate ``(0, 0, 0, 0)`` tensor instead of saying what is
    wrong."""
    if isinstance(norb, bool) or not isinstance(norb, (int, np.integer)) or norb < 1:
        raise ValueError(
            "flex_second_order = \"local\": norb must be a positive integer, got {!r}"
            .format(norb))


def _check_rows(itype, tbl, norb):
    """Row-level validation of one type's on-site table: orbital indices
    inside the model, and finite coefficients.

    Both are refusals a programmatic caller can trigger by bypassing the
    reader. Without them a negative orbital index wraps silently (numpy
    indexes from the end), a large one raises an ``IndexError`` that names
    no row, and a non-finite coupling propagates into ``Gamma`` and then
    past the Hermiticity guard below -- ``dev > herm_tol * scale`` is FALSE
    when ``dev`` is NaN, so a NaN table would compile without complaint."""
    for (a, b), v in tbl.items():
        if not (0 <= a < norb and 0 <= b < norb):
            raise ValueError(
                "flex_second_order = \"local\": the on-site {} row (orbitals {}, {}) is "
                "outside the {} orbitals of the model".format(itype, a + 1, b + 1, norb))
        if not np.isfinite(complex(v).real) or not np.isfinite(complex(v).imag):
            raise ValueError(
                "flex_second_order = \"local\": the on-site {} row (orbitals {}, {}) has a "
                "non-finite coefficient {!r}".format(itype, a + 1, b + 1, v))


def _check_hermitian_closure(raw_rows, herm_tol=1e-12):
    """Refuse a table whose transposed row is not the Hermitian partner of
    its row, NAMING the offending type and row.

    Every accepted type's transposed entry carries the adjoint of the row's
    operator -- the same operator for the mirrored types
    (:data:`_MIRRORED_TYPES`), the reverse pair hop for ``PairHop`` -- so a
    Hermitian Hamiltonian means ``v_{ba} = conj(v_{ab})`` per declared pair
    of ordered rows, which is what the k-space readers have enforced since
    #93. The reversal closure (:func:`close_onsite_rows`) SYMMETRISES
    whatever it is given, so a table that violates this is not rejected by
    the closure: it is quietly replaced by its Hermitian part (a complex
    same-operator row by its real part, say), and only the tensor-wide guard
    in :func:`compile_onsite` can notice -- and then only when the
    replacement happens to leave a non-Hermitian ``Gamma``, with a
    diagnostic that names no row.

    The check is on the RAW rows and applies only where BOTH ordered rows of
    a pair are declared: a lone row is legal input, and the closure supplies
    its conjugate partner (see :func:`close_onsite_rows`). ``a == b`` rows
    are their own partner, so for them the condition is that the coefficient
    is real."""
    for itype, tbl in raw_rows.items():
        for (a, b), v in tbl.items():
            if (b, a) not in tbl:
                continue                    # a lone row: the closure conjugates it
            w = tbl[(b, a)]
            dev = abs(complex(v) - np.conj(complex(w)))
            scale = max(1.0, abs(complex(v)), abs(complex(w)))
            if dev > herm_tol * scale:
                raise ValueError(
                    "flex_second_order = \"local\": the on-site interaction is not "
                    "Hermitian-closed: the {} row (orbitals {}, {}) is {} but its "
                    "transposed row is {}, and a Hermitian table needs the complex "
                    "conjugate {} there".format(itype, a + 1, b + 1, complex(v),
                                                complex(w), np.conj(complex(v))))


def compile_onsite_v(onsite_tbl, norb, *, closed=False):
    """The non-antisymmetrised ``V`` ``(M, M, M, M)``, ``M = 2 norb``, from the
    on-site table ``{type: {((0,0,0),(a,b)): v}}`` (with ``closed=True`` the rows
    are taken as already reversal-closed -- tests only). Refuses degenerate
    ``a == b`` rows of the two-orbital types (spec D7) and unknown types."""
    _check_norb(norb)
    raw = _raw_onsite_rows(onsite_tbl)
    # validate the DECLARED rows, before the closure touches them: the
    # closure averages a row with its transpose, and a non-finite
    # coefficient would first surface there as a numpy warning about an
    # arithmetic the caller never asked for.
    for itype, tbl in raw.items():
        _check_rows(itype, tbl, norb)
    rows = raw if closed else close_onsite_rows(onsite_tbl)
    M = 2 * norb
    V = np.zeros((M, M, M, M), dtype=np.complex128)
    for itype, tbl in rows.items():
        if itype == "CoulombIntra":
            pass
        elif itype in TWO_ORBITAL_TYPES:
            for (a, b) in tbl:
                if a == b:
                    raise DegenerateRowError(
                        "flex_second_order = \"local\": the on-site same-orbital {} row (orbital "
                        "{}) is not a two-body term: {}".format(itype, a + 1, _D7_REMEDY[itype]))
        else:
            raise ValueError("second_order: unknown interaction type {!r}".format(itype))
        for (a, b), v in tbl.items():
            if itype == "CoulombIntra" and a != b:
                continue          # the kernel drops these too (discarded report)
            for (p, q, s_, r, x) in _records(itype, a, b, v, norb):
                V[p, q, r, s_] += x
                V[q, p, s_, r] += x
    return V


def compile_onsite(onsite_tbl, norb, *, herm_tol=1e-12, closed=False):
    """Gamma_on (M, M, M, M), M = 2 norb, from the on-site table
    {type: {((0,0,0),(a,b)): v}} (or, with closed=True, rows already closed --
    tests only). ``Gamma[p, q, r, s] = V[p, q, r, s] - V[p, q, s, r]``.
    Refuses degenerate a = b rows of the two-orbital types (D7) and a
    non-Hermitian result."""
    V = compile_onsite_v(onsite_tbl, norb, closed=closed)
    # Per-row first: it names the type and the row. The tensor-wide guard
    # below stays as the last resort -- it is the one that would catch a
    # convention error inside the record table itself, which no row check
    # can see.
    _check_hermitian_closure(_raw_onsite_rows(onsite_tbl), herm_tol)
    Gamma = V - V.transpose(0, 1, 3, 2)
    herm = Gamma.conj().transpose(2, 3, 0, 1)
    dev = float(np.max(np.abs(Gamma - herm))) if Gamma.size else 0.0
    scale = max(1.0, float(np.max(np.abs(Gamma)))) if Gamma.size else 1.0
    if dev > herm_tol * scale:
        raise ValueError(
            "flex_second_order = \"local\": the compiled on-site interaction is not Hermitian "
            "(relative deviation {:.3e}); the interaction table is not Hermitian-closed"
            .format(dev / scale))
    return Gamma


def hf_first_order(tensor, rho):
    """TEST HELPER (host/numpy only, not used by any solve): the first-order
    (mean-field) self-energy of a rank-4 on-site tensor.

    The production first order is ``hartree_fock.accumulate_hf``; this is the
    same functional written as one contraction, so the compiler's ``Gamma``
    can be confronted with it (spec 2.2 contract (i)). It materialises the
    dense ``(2 norb)^4`` tensor and runs ``np.einsum``, so it is host-side by
    construction.

    Contraction:

        Sigma1[p, r] = sum_{q s} tensor[p, q, r, s] rho[q, s],
        rho[q, s] = <c^dag_q c_s>.

    Which contribution comes out is chosen by WHICH tensor the caller passes:
    ``Gamma`` from :func:`compile_onsite` gives Hartree + Fock (it carries the
    exchange placement), the non-antisymmetrised ``V`` from
    :func:`compile_onsite_v` gives the Hartree (direct) part alone. The
    contraction itself is the same, so there is no flag to get wrong.
    """
    return np.einsum('pqrs,qs->pr', np.asarray(tensor), rho)


def density_slots(Gamma, norb):
    """d[s, s', a, b] = V[(a s),(b s'),(a s),(b s')] read from Gamma
    (for (a s) != (b s') the direct entry is Gamma[p, q, p, q])."""
    d = np.zeros((2, 2, norb, norb), dtype=np.complex128)
    for s1, s2, a, b in itertools.product((UP, DN), (UP, DN), range(norb), range(norb)):
        p, q = gen_index(a, s1, norb), gen_index(b, s2, norb)
        if p != q:
            d[s1, s2, a, b] = Gamma[p, q, p, q]
    return d


_SPIN_WEIGHT = {   # w[itype][(s, s')]: the density-slot spin structure of spec 2.4
    "CoulombInter": lambda s1, s2: 1.0,
    "Hund": lambda s1, s2: -1.0 if s1 == s2 else 0.0,
    "Ising": lambda s1, s2: 1.0 if s1 == s2 else -1.0,
}


def build_offsite(split, lattice, norb):
    """vpair[s, s', q, a, b] = v^{pair, s s'}_{(aa),(bb)}(q) = v^{file, s s'}_{ba}(q):
    the off-site density-slot arrays in the ED-validated pair-space
    orientation of ``hwave.sc._build_interaction_k`` (transpose=True), with
    the reader's reversal closure (``symmetrise_k``). None without off-site
    density types."""
    from hwave.sc import _build_interaction_k
    from hwave.solver.declarations import symmetrise_k
    types = [t for t in OFFSITE_DENSITY_TYPES if split.offsite_tbl.get(t)]
    if not types:
        return None
    nx, ny, nz = (int(x) for x in lattice.shape)
    nvol = nx * ny * nz
    kx = np.linspace(0, 2.0 * np.pi, nx, endpoint=False)
    ky = np.linspace(0, 2.0 * np.pi, ny, endpoint=False)
    kz = np.linspace(0, 2.0 * np.pi, nz, endpoint=False)
    inter_k = _build_interaction_k(kx, ky, kz, {t: split.offsite_tbl[t] for t in types}, norb)
    inter_k = symmetrise_k(inter_k)
    vpair = np.zeros((2, 2, nvol, norb, norb), dtype=np.complex128)
    for t in types:
        M = np.asarray(inter_k[t]).reshape(norb, norb, nvol).transpose(2, 0, 1)   # (nvol, x, y)
        # _build_interaction_k(transpose=True): M[q, x, y] = v^{file}_{yx}(q) = v^{pair}_{(xx),(yy)}
        for s1, s2 in itertools.product((UP, DN), repeat=2):
            w = _SPIN_WEIGHT[t](s1, s2)
            if w != 0.0:
                vpair[s1, s2] += w * M
    return vpair


def factor_bytes(norb, nvol, offsite):
    nd = norb * norb
    return 16 * (8 * nd * nd + (4 * nvol * norb * norb if offsite else 0))


@_dataclass(frozen=True)
class SecondOrderFactors:
    norb: int
    nd: int
    nvol: int
    triples: tuple          # (sigma_s, sigma_r, sigma_q) of the nonzero on-site pairs
    A_on: tuple             # (nd, nd) complex128 each, read-only
    B_on: tuple
    vpair: object           # (2, 2, nvol, norb, norb) or None
    nbytes: int


def _readonly(a):
    a = np.ascontiguousarray(a)
    a.setflags(write=False)
    return a


def build_factors(split, lattice, norb):
    """Spec 2.5: the sparse factors of the local kernel, from the locality
    split (on-site rows -> Gamma_on -> A/B per nonzero spin triple; off-site
    density rows -> vpair)."""
    nd = norb * norb
    nvol = int(np.prod([int(x) for x in lattice.shape]))
    Gamma = compile_onsite(split.onsite_tbl, norb)
    triples, As, Bs = [], [], []
    for ss, sr, sq in itertools.product((UP, DN), repeat=3):
        A = np.zeros((nd, nd), dtype=np.complex128)
        B = np.zeros((nd, nd), dtype=np.complex128)
        for c, a, r, q in itertools.product(range(norb), repeat=4):
            A[c * norb + a, r * norb + q] = Gamma[gen_index(a, UP, norb), gen_index(q, sq, norb),
                                                  gen_index(r, sr, norb), gen_index(c, ss, norb)]
            B[r * norb + q, c * norb + a] = Gamma[gen_index(r, sr, norb), gen_index(c, ss, norb),
                                                  gen_index(a, UP, norb), gen_index(q, sq, norb)]
        if np.abs(A).max() == 0.0 and np.abs(B).max() == 0.0:
            continue
        triples.append((ss, sr, sq))
        As.append(_readonly(A))
        Bs.append(_readonly(B))
    vpair = build_offsite(split, lattice, norb)
    if vpair is not None:
        vpair = _readonly(vpair)
    nbytes = 16 * (2 * len(triples) * nd * nd + (4 * nvol * norb * norb if vpair is not None else 0))
    return SecondOrderFactors(norb=norb, nd=nd, nvol=nvol, triples=tuple(triples), A_on=tuple(As),
                              B_on=tuple(Bs), vpair=vpair, nbytes=nbytes)


def _dens_slice(norb):
    """The ``(a, a)`` density pair slots as a BASIC slice.

    Slot ``(a, a)`` sits at pair index ``a * norb + a = a * (norb + 1)``, so
    the density slots of the ``nd = norb^2`` pair axis are the arithmetic
    sequence ``0, norb + 1, ..., (norb - 1) (norb + 1)``: a regular stride.
    Selecting them with this slice is BASIC indexing, so every gather is a
    view and every ``out[..., slots] +=`` is a genuine in-place update --
    no copy of the gathered block and no scatter temporary, which is what
    keeps :func:`accumulate_batch` inside its documented budget.
    """
    return slice(None, None, norb + 1)


def _carve(T, shape):
    """A C-contiguous view of ``shape`` carved out of the front of the
    C-contiguous scratch buffer ``T``.

    The off-site terms need ``(nb, nvol, nd, norb)``-shaped and smaller
    output arrays, and the natural sub-block of ``T`` -- ``T[:, :, :, :norb]``
    and friends -- is STRIDED. numpy's ``matmul`` accepts a strided ``out=``,
    but that is not portable (CuPy's need not), so the scratch is taken from
    the buffer's flat memory instead: same bytes, no allocation, contiguous.
    """
    n = 1
    for d in shape:
        n *= int(d)
    return T.reshape(-1)[:n].reshape(shape)


def _byte_extent(arr):
    """``(lo, hi)`` -- the byte offsets, relative to the array's own base
    pointer, that ``arr`` actually spans, computed from its shape and
    strides. ``lo`` is <= 0 for an array with a negative stride.

    ``nbytes`` is NOT this: it counts the bytes the array's ELEMENTS occupy,
    which for any strided view is smaller than the region the view reaches
    across. ``a[:, ::2]`` has half the ``nbytes`` of its base but spans
    nearly all of it, so an overlap test built on ``[ptr, ptr + nbytes)``
    would call such a view disjoint from a buffer that sits in the gap it
    strides over -- a false NEGATIVE, which for an aliasing guard is the
    dangerous direction.

    An empty array spans nothing and is reported as ``(0, 0)``.
    """
    shape = tuple(int(x) for x in arr.shape)
    if any(n == 0 for n in shape):
        return 0, 0
    lo = hi = 0
    for n, stride in zip(shape, tuple(int(x) for x in arr.strides)):
        step = (n - 1) * stride
        if step < 0:
            lo += step
        else:
            hi += step
    return lo, hi + int(arr.itemsize)


def _shares_storage(a, b):
    """``True`` when ``a`` and ``b`` MAY share memory, ``False`` when they
    are known not to, ``None`` when the array module cannot say.

    numpy answers exactly (``np.shares_memory``). For a device module the
    fallback is an EXTENT overlap test on the absolute byte ranges
    ``[ptr + lo, ptr + hi)`` of :func:`_byte_extent`, using the base pointer
    the module exposes (CuPy: ``arr.data.ptr``, which already includes a
    view's own offset). The rule is deliberately CONSERVATIVE in one
    direction: two arrays whose extents overlap need not share a single
    element (interleaved strides do not), so this can report ``True`` for a
    disjoint pair -- a refusal the caller can work around by lending
    contiguous buffers -- but it cannot report ``False`` for a pair that
    does overlap, which is the answer that would corrupt a result. A module
    that exposes no base pointer returns ``None`` and the caller documents
    the precondition instead of guessing.
    """
    if a is b:
        return True
    if isinstance(a, np.ndarray) and isinstance(b, np.ndarray):
        return bool(np.shares_memory(a, b))
    pa = getattr(getattr(a, "data", None), "ptr", None)
    pb = getattr(getattr(b, "data", None), "ptr", None)
    if pa is None or pb is None:
        return None
    lo_a, hi_a = _byte_extent(a)
    lo_b, hi_b = _byte_extent(b)
    if hi_a == lo_a or hi_b == lo_b:
        return False                      # one of them spans nothing
    return not (pa + hi_a <= pb + lo_b or pb + hi_b <= pa + lo_a)


def _validate_batch(out_b, chibar_b, factors, work):
    """Shape and aliasing preconditions of :func:`accumulate_batch`.

    Every one of these is a silent-corruption route rather than a crash:
    a ``chibar_b`` whose pair axes do not match the factors contracts the
    wrong slots (or raises deep inside a ``matmul`` with a message naming
    no argument), a mismatched ``out_b`` broadcasts instead of accumulating,
    and a lent buffer that overlaps an operand is overwritten mid-product.
    They are checked BEFORE any work is done, including when ``triples`` is
    empty and no product runs at all."""
    nd = factors.nd
    if nd != factors.norb ** 2:
        raise ValueError(
            "accumulate_batch: the factor pack is inconsistent, nd = {} but norb^2 = {}"
            .format(nd, factors.norb ** 2))
    # The volume axis is constrained only when the pack carries an off-site
    # vertex: the on-site factors are (nd, nd) matrices broadcast over q, so
    # an on-site-only pack is usable at ANY nvol (the bond store's synthetic
    # fixtures rely on that). vpair is indexed by q and must match.
    nvol = factors.nvol if factors.vpair is not None else None
    for name, arr in (("chibar_b", chibar_b), ("out_b", out_b)):
        shape = tuple(int(x) for x in arr.shape)
        ok = (len(shape) == 4 and shape[2] == nd and shape[3] == nd
              and (nvol is None or shape[1] == nvol))
        if not ok:
            raise ValueError(
                "accumulate_batch: {} has shape {}, expected (nb, {}, {}, {}) from the "
                "factor pack (nvol, nd, nd)".format(
                    name, shape, factors.nvol if nvol is not None else "nvol", nd, nd))
    if tuple(out_b.shape) != tuple(chibar_b.shape):
        raise ValueError(
            "accumulate_batch: out_b has shape {} but chibar_b has {}; the two must "
            "cover the same frequency batch".format(tuple(out_b.shape), tuple(chibar_b.shape)))
    # The kernel ACCUMULATES into out_b while still reading chibar_b: the
    # on-site sandwich reads the bubble again on every spin triple, and the
    # off-site terms read its density slots after out_b has been written.
    # Sharing storage between the two is the same silent corruption as a
    # lent buffer overlapping an operand, so it is refused the same way.
    if _shares_storage(out_b, chibar_b):
        raise ValueError(
            "accumulate_batch: out_b shares memory with chibar_b; the kernel accumulates "
            "into out_b while it is still reading chibar_b")
    if work is None:
        return
    T1, T2 = work
    if T1 is T2:
        raise ValueError(
            "accumulate_batch: work[0] and work[1] are the same array; the kernel needs "
            "two independent scratch buffers")
    for name, T in (("work[0]", T1), ("work[1]", T2)):
        if tuple(T.shape) != tuple(chibar_b.shape):
            raise ValueError(
                "accumulate_batch: {} has shape {}, expected chibar_b's {}"
                .format(name, tuple(T.shape), tuple(chibar_b.shape)))
        for oname, operand in (("chibar_b", chibar_b), ("out_b", out_b)):
            if _shares_storage(T, operand):
                raise ValueError(
                    "accumulate_batch: {} shares memory with {}; the scratch buffers are "
                    "overwritten while both operands are still being read"
                    .format(name, oname))
    if _shares_storage(T1, T2):
        raise ValueError(
            "accumulate_batch: work[0] and work[1] share memory; the kernel needs two "
            "independent scratch buffers")


def accumulate_batch(out_b, chibar_b, l0, factors, work=None):
    """Add W2 (spec 2.3 + 2.4) for the bosonic frequencies [l0, l0 + nb)
    into out_b in place. Two (nb, nvol, nd, nd) temporaries, reused.

    ``work`` optionally supplies those two temporaries as a pair
    ``(T1, T2)`` of arrays shaped exactly like ``chibar_b``, on the same
    array module. A caller that already holds scratch of that shape --
    the general path's batch loop does (:meth:`hwave.solver.flex.FLEX.
    _calc_veff_general`) -- passes it here so the whole assembly peaks at
    TWO ``(nb, nvol, nd, nd)`` temporaries rather than four, which is the
    ``T_bytes = 2 nb nvol nd^2 * 16`` budget of spec 2.5. With ``None``
    the buffers are allocated per call, as before. The buffers are
    scratch: their contents on entry are irrelevant and are overwritten.

    Preconditions are validated before any work is done, including when the
    factor pack has no nonzero spin triple at all (:func:`_validate_batch`):
    the pair axes of ``chibar_b``/``out_b`` against the factor pack, the two
    arrays against each other, and the lent buffers against BOTH operands
    for aliasing (numpy exactly; on a device module the pointer ranges where
    the module exposes them, and otherwise not at all -- see
    :func:`_shares_storage`).

    The lent buffers must be C-contiguous, which is how both production
    callers allocate them: the off-site terms need sub-blocks of them as
    ``matmul`` output (:func:`_carve`), and a strided ``out=`` is not
    portable across array modules (numpy accepts it, CuPy need not). A
    non-contiguous pair of the right shape is not an error -- the kernel
    then allocates its own two contiguous buffers and ignores the loan --
    while a mis-shaped one still raises.

    Allocation contract (pinned by
    ``tests/test_second_order_kernel.py::TestKernel::
    test_allocation_peak_with_lent_buffers``): with C-contiguous ``work``
    supplied (a non-contiguous loan is discarded and two buffers allocated,
    see above) the kernel allocates NOTHING that scales with
    ``nb * nvol * nd^2``. Every
    product is written with ``matmul(..., out=)`` into a CONTIGUOUS view of
    the lent buffers, and every density-slot gather and scatter goes through
    :func:`_dens_slice`, i.e. basic indexing -- a view, not a fancy-index
    copy. The one allocation left is the boolean mask of the finiteness
    checkpoint below, ``nb * nvol * nd^2`` bytes (a sixteenth of one
    complex batch array, ``T_bytes / 32``); on top of that only array-module
    scratch that does not scale with the batch (matmul's own buffering for
    strided operands, a few KiB).
    """
    # FIRST: nothing below may index an axis the arrays need not have. A
    # rank-0/1/2 input used to raise IndexError from the unpacking on the
    # next line, naming no argument, before the guard could say what was
    # wrong.
    _validate_batch(out_b, chibar_b, factors, work)
    xp = _bk.array_module_of(chibar_b)
    norb = factors.norb
    nb, nvol, nd = chibar_b.shape[0], chibar_b.shape[1], chibar_b.shape[2]
    lent = work is not None
    if lent:
        T1, T2 = work
        # the off-site terms carve contiguous sub-blocks out of these
        lent = bool(T1.flags.c_contiguous and T2.flags.c_contiguous)
    if not lent:
        # xp.empty (not empty_like): C-contiguous whatever chibar_b's layout
        T1 = xp.empty(tuple(chibar_b.shape), dtype=chibar_b.dtype)
        T2 = xp.empty(tuple(chibar_b.shape), dtype=chibar_b.dtype)
    ds = _dens_slice(norb)
    vp = None
    if factors.vpair is not None:
        # A_v is the crossed entry Gamma_{(a up)(q sigma),(q sigma)(a up)} =
        # -v^{up sigma}_{a q}(q-vector), i.e. the orbital of the EXTERNAL leg
        # indexes the ROW. :func:`build_offsite` stores the ring's pair-space
        # placement v^{pair}_{(aa),(bb)} = v^{file}_{ba} -- the row index is the
        # orbital of the row's DISPLACED cell. Once a file row (r, a, b, v) is
        # read as v n_{j,a} n_{j+r,b} (the DOCUMENTED orientation, which the
        # mean-field kernel also implements since spec 2026-09-16, #193), those
        # two are the same matrix, so the kernel uses the stored array as it is
        # -- no transpose here. On a bond with v_{ab}(R) != v_{ba}(R) the
        # reversed reading is a different matrix (the two differ by q -> -q);
        # the one below is what the independent real-space oracle requires on
        # every pair of coupling types (tests/test_second_order_oracle.py,
        # gate G2 (a)).
        vp = xp.asarray(factors.vpair)                          # (2, 2, nvol, norb, norb)
    # on-site: 1/2 sum_triples A chibar B, and -- for the triples the mixed
    # term selects -- the two on-site/off-site cross terms of the SAME
    # triple, which reuse the ``A chibar`` product in T1 instead of
    # recomputing it.
    for (ss, sr, sq), A, B in zip(factors.triples, factors.A_on, factors.B_on):
        Bb = xp.asarray(B)[None, None]
        xp.matmul(xp.asarray(A)[None, None], chibar_b, out=T1)
        xp.matmul(T1, Bb, out=T2)
        T2 *= 0.5
        out_b += T2
        if vp is None or ss != UP or sr != sq:
            continue
        sig = sq
        # A_on chibar B_v : (A chibar)[:, :, :, dens] @ (-vpair[sig, up]).
        # T2 is free again (its content is already in out_b); T1 still holds
        # A chibar, whose density-slot columns are the left operand.
        P = _carve(T2, (nb, nvol, nd, norb))
        xp.matmul(T1[:, :, :, ds], vp[sig, UP][None], out=P)
        out_b[:, :, :, ds] -= P
        # A_v chibar B_on : (-vpair[up, sig]) @ chibar[:, :, dens, :] @ B.
        # The first product lands in T2 (P has been consumed), the second in
        # T1 (A chibar is no longer needed in this iteration).
        Q = _carve(T2, (nb, nvol, norb, nd))
        xp.matmul(vp[UP, sig][None], chibar_b[:, :, ds, :], out=Q)
        R = _carve(T1, (nb, nvol, norb, nd))
        xp.matmul(Q, Bb, out=R)
        out_b[:, :, ds, :] -= R
    if vp is not None:
        # A_v chibar B_v : sum_sig (-vpair[up,sig]) chibar[dens,dens] (-vpair[sig,up]).
        # blk is a view and does not depend on sig.
        blk = chibar_b[:, :, ds, ds]                           # (nb, nvol, norb, norb)
        X = _carve(T1, (nb, nvol, norb, norb))
        Y = _carve(T2, (nb, nvol, norb, norb))
        for sig in (UP, DN):
            xp.matmul(vp[UP, sig][None], blk, out=X)
            xp.matmul(X, vp[sig, UP][None], out=Y)
            out_b[:, :, ds, ds] += Y
    if not bool(xp.all(xp.isfinite(out_b))):
        raise NonFiniteError("non-finite second-order kernel W2 in the frequency batch [{}, {})"
                             .format(l0, l0 + nb))
    return out_b


def dense_w2(chibar, factors):
    """TEST HELPER (host/numpy only, not used by any solve): ``W2``
    materialised for the WHOLE frequency axis at once.

    Production never does this -- :meth:`hwave.solver.flex.FLEX.
    _calc_veff_general` accumulates in frequency batches precisely so that no
    full-size ``W2`` is ever allocated -- so this helper is for tests small
    enough to hold one. ``np.asarray`` pulls its argument to the host, which
    is what makes it numpy-only."""
    out = np.zeros_like(np.asarray(chibar))
    accumulate_batch(out, np.asarray(chibar), 0, factors)
    return out
