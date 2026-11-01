"""Test-side oracles for the spin-orbital bond-resolved vertex.

Everything here is an independent, brute-force evaluation of the label-form
formulas below on the production grids, and an exact diagonalization built
from the interaction tables. The module never imports
``hwave.solver.so_bond`` or ``hwave.solver.so_records``: the maps are
duck-typed on a geometry object carrying ``delta_r`` (B, 3), ``reverse`` (B,)
and ``spatial_shape`` (Nx, Ny, Nz), every convention choice is an explicit
flag tuple, and the operators come from the table-driven builder
:func:`build_terms`, so a production compiler error is never common-mode.

Grids (the production ones)
---------------------------
* Momenta are in C order over ``spatial_shape``: grid index
  ``(i, j, l) -> (i * Ny + j) * Nz + l``, ``k = 2 pi (i/Nx, j/Ny, l/Nz)``.
* Fermionic ``omega_n = (2n + 1 - nmat) pi T``; bosonic
  ``nu_l = (2l - nmat) pi T``, zero at ``l = nmat // 2`` (``nmat`` even).

Label and code indices
----------------------
* A label ``(d; X, Y)`` names ``rho^{(d)}_{XY}(q) = sum_k e^{+ik.d}
  c+_{X,k+q} c_{Y,k}``. The label displacements are the geometry's own
  ``delta_r`` rows, so the label index of ``(delta_r[m]; X, Y)`` is
  ``m * nd**2 + X * nd + Y`` (:func:`bond_label_index`).
* A code index ``(m; x, y)`` has the same flat position
  ``m * nd**2 + x * nd + y``; it is the slot of the production bond bubble
  (``bubble.bond_bubble_static``: ``I = m * npair + (a * p + c)``).

Index-convention boundary
-------------------------
* Monomials (:func:`expand_tables`), labels, :func:`file_gamma_label`, the
  Lindhard functions and every code/label map use the SPIN-BLOCK mode index
  ``X = s * norb_phys + a`` (``s = 0`` up, ``1`` down).
* :func:`build_terms`, :func:`build_h1`, :func:`apply_terms`, :class:`JWED`
  and :func:`rho_matrix` use the INTERLEAVED exact-diagonalization mode
  ``nd * j + 2 a + s`` (``tests.ed_oracle_util.EDFixture.mode``), ``j`` the
  C-order cell index. The one conversion is :func:`ed_mode`.
* Transfer tables are keyed by the interleaved in-cell index ``A = 2 a + s``
  (the spin-orbital file format); interaction tables by physical orbitals.
"""

import functools
import itertools
from typing import NamedTuple

import numpy as np
import scipy.sparse as sp

# ---------------------------------------------------------------------------
# Shared convention fixtures: the 5x3x1 grid with on-site, x, y and diagonal bonds
# ---------------------------------------------------------------------------

CONVENTION_SHAPE = (5, 3, 1)
CONVENTION_DR = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0],
                          [1, 1, 0], [-1, -1, 0]], dtype=np.int64)
CONVENTION_REV = np.array([0, 2, 1, 4, 3, 6, 5], dtype=np.int64)
CONVENTION_ND = 2
CONVENTION_BETA = 2.0
CONVENTION_NMAT_STATIC = 8
CONVENTION_NMAT_DYNAMIC = 6
CONVENTION_GREEN_SEEDS = (2701, 2702)
CONVENTION_SAMPLE_SEED = 2707
CONVENTION_SAMPLE_SIZE = 256


def random_green(seed, nmat, nvol, nd):
    """A random NON-factorized Green function ``(nmat, nvol, nd, nd)``.

    Every entry, the spin-off-diagonal ones included, is an independent
    complex normal number, so no symmetry of a physical G can hide a
    convention error. Production entry points take the leading block axis
    ``green[None]`` (``(1, nmat, nvol, nd, nd)``)."""
    rng = np.random.default_rng(seed)
    shape = (nmat, nvol, nd, nd)
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


# ---------------------------------------------------------------------------
# Momentum-grid helpers
# ---------------------------------------------------------------------------

def mesh_grid(spatial_shape):
    """Integer grid indices ``(nvol, 3)`` in C order."""
    shape = tuple(int(x) for x in spatial_shape)
    return np.array(list(itertools.product(*(range(n) for n in shape))),
                    dtype=np.int64).reshape(-1, 3)


def mesh_index(vectors, spatial_shape):
    """Flat C-order index of integer grid vectors (wrapped onto the mesh)."""
    shape = np.array(spatial_shape, dtype=np.int64)
    v = np.mod(np.asarray(vectors, dtype=np.int64), shape)
    return (v[..., 0] * shape[1] + v[..., 1]) * shape[2] + v[..., 2]


def mesh_momenta(spatial_shape):
    """Momenta ``k = 2 pi (i/Nx, j/Ny, l/Nz)`` ``(nvol, 3)`` in C order."""
    return 2.0 * np.pi * mesh_grid(spatial_shape) / np.array(spatial_shape, float)


def mesh_negate(spatial_shape):
    """``neg[q]`` is the index of ``-q``."""
    return mesh_index(-mesh_grid(spatial_shape), spatial_shape)


def mesh_add(spatial_shape):
    """``add[q, k]`` is the index of ``k + q``."""
    g = mesh_grid(spatial_shape)
    return mesh_index(g[:, None, :] + g[None, :, :], spatial_shape)


# ---------------------------------------------------------------------------
# Label form
# ---------------------------------------------------------------------------

def bond_label_index(geom, nd):
    """``{(d, X, Y): m * nd**2 + X * nd + Y}`` with ``d = tuple(delta_r[m])``.

    Refuses a geometry with a repeated displacement (the label would be
    ambiguous)."""
    out = {}
    for m, d in enumerate(np.asarray(geom.delta_r, dtype=np.int64)):
        d = tuple(int(x) for x in d)
        for x in range(nd):
            for y in range(nd):
                key = (d, x, y)
                if key in out:
                    raise ValueError("bond_label_index: displacement {} appears twice"
                                     .format(d))
                out[key] = m * nd * nd + x * nd + y
    return out


def lindhard_discrete(green, beta, geom):
    """The section-3.0 label bubble on the production grid, by direct sum.

    ``chi[l, q, (d; X, Y), (d'; Z, W)] = -(T/N) sum_{k, n} e^{+ik.d}
    e^{+i(k+q).d'} G_YZ(k, n) G_WX(k+q, n')`` with the fermionic index of
    the second leg wrapped circularly, ``n' = (n - l + nmat/2) mod nmat``
    (the shift by ``nu_l``), and no high-frequency tail.

    Parameters
    ----------
    green : ndarray (nmat, nvol, nd, nd)
        No block axis.
    beta : float
    geom : object with ``delta_r`` and ``spatial_shape``
        The labels are the rows of ``delta_r``.

    Returns
    -------
    ndarray complex128 (nmat, nvol, NL, NL), ``NL = B * nd**2``.
    """
    green = np.asarray(green, dtype=np.complex128)
    nmat, nvol, nd, _ = green.shape
    shape = tuple(int(x) for x in geom.spatial_shape)
    if int(np.prod(shape)) != nvol:
        raise ValueError("lindhard_discrete: green has nvol={} but the mesh {} has {}"
                         .format(nvol, shape, int(np.prod(shape))))
    D = np.asarray(geom.delta_r, dtype=np.int64)
    B = D.shape[0]
    NL = B * nd * nd
    kv = mesh_momenta(shape)
    pk = np.exp(1j * (kv @ D.T.astype(float)))            # (k, a)
    add = mesh_add(shape)                                  # (q, k) -> k + q
    first = pk[None, :, :, None, None] * green[:, :, None, :, :]      # (n, k, a, y, z)
    pkq = pk[add]                                          # (q, k, b)
    out = np.empty((nmat, nvol, NL, NL), dtype=np.complex128)
    for l in range(nmat):
        shifted = green[(np.arange(nmat) - l + nmat // 2) % nmat]    # (n, k, w, x)
        second = pkq[None, :, :, :, None, None] * shifted[:, add][:, :, :, None, :, :]
        t = np.einsum("nkayz,nqkbwx->qaxybzw", first, second, optimize=True)
        out[l] = (-(1.0 / beta) / nvol) * t.reshape(nvol, NL, NL)
    return out


# ---------------------------------------------------------------------------
# Candidate label -> code maps
# ---------------------------------------------------------------------------

def _check_flags(flags):
    flags = tuple(int(f) for f in flags)
    if len(flags) != 9:
        raise ValueError("flags must be (lambda_r, lambda_c, s_q, s_nu, rowrev, colrev, t, "
                         "p_r, p_c), got {}".format(flags))
    lam_r, lam_c, s_q, s_nu, rr, cr, t, p_r, p_c = flags
    if any(v not in (1, -1) for v in (lam_r, lam_c, s_q, s_nu)) \
            or any(v not in (0, 1) for v in (rr, cr, t)) \
            or any(v not in (-1, 0, 1) for v in (p_r, p_c)):
        raise ValueError("flags out of range: {}".format(flags))
    return flags


def candidate_map(chi_label, geom, nd, flags, *, entries=None):
    """One candidate label -> code map of the bond bubble.

    With code row ``r = (m; x, y)`` and column ``c = (m'; z, w)``,
    ``L_r(r) = (lambda_r delta_r[m]; (x, y), or (y, x) if rowrev)``,
    ``L_c(c) = (lambda_c delta_r[m']; (z, w), or (w, z) if colrev)``:

    * ``t = 0``: ``code[r, c](q, nu) = e^{iq.(p_r dr_m + p_c dr_m')}
      chi(s_q q, s_nu nu)[L_r(r), L_c(c)]``;
    * ``t = 1``: the same with ``chi(s_q q, s_nu nu)[L_c(c), L_r(r)]``.

    The phase uses the code momentum ``q``; ``-nu_l`` is the slice
    ``(nmat - l) mod nmat``.

    Parameters
    ----------
    chi_label : ndarray (nmat, nvol, NL, NL)
        From :func:`lindhard_discrete` on the same geometry.
    flags : 9-tuple ``(lambda_r, lambda_c, s_q, s_nu, rowrev, colrev, t, p_r, p_c)``
    entries : tuple of four int arrays ``(l, q, r, c)``, optional
        Evaluate only these code entries (returns a 1-D array); default:
        the full ``(nmat, nvol, NC, NC)`` array.
    """
    lam_r, lam_c, s_q, s_nu, rr, cr, t, p_r, p_c = _check_flags(flags)
    chi_label = np.asarray(chi_label)
    nmat, nvol = chi_label.shape[:2]
    shape = tuple(int(x) for x in geom.spatial_shape)
    D = np.asarray(geom.delta_r, dtype=np.int64)
    B = D.shape[0]
    npair = nd * nd
    lab = bond_label_index(geom, nd)

    def slot_labels(lam, rev):
        out = np.empty(B * npair, dtype=np.int64)
        for m in range(B):
            d = tuple(int(v) for v in lam * D[m])
            if (d, 0, 0) not in lab:
                raise ValueError("candidate_map: displacement {} is not a label of the "
                                 "geometry".format(d))
            for x in range(nd):
                for y in range(nd):
                    a, b = (y, x) if rev else (x, y)
                    out[m * npair + x * nd + y] = lab[(d, a, b)]
        return out

    lr = slot_labels(lam_r, rr)
    lc = slot_labels(lam_c, cr)
    qsel = np.arange(nvol) if s_q > 0 else mesh_negate(shape)
    lsel = np.arange(nmat) if s_nu > 0 else (nmat - np.arange(nmat)) % nmat
    qd = mesh_momenta(shape) @ D.T.astype(float)            # (q, m) = q . dr_m
    ph_r = np.repeat(np.exp(1j * p_r * qd), npair, axis=1)  # (q, r)
    ph_c = np.repeat(np.exp(1j * p_c * qd), npair, axis=1)  # (q, c)

    if entries is None:
        L = lsel[:, None, None, None]
        Q = qsel[None, :, None, None]
        R = lr[None, None, :, None]
        C = lc[None, None, None, :]
        vals = chi_label[L, Q, R, C] if t == 0 else chi_label[L, Q, C, R]
        return vals * ph_r[None, :, :, None] * ph_c[None, :, None, :]
    l, q, r, c = (np.asarray(e, dtype=np.int64) for e in entries)
    if t == 0:
        vals = chi_label[lsel[l], qsel[q], lr[r], lc[c]]
    else:
        vals = chi_label[lsel[l], qsel[q], lc[c], lr[r]]
    return vals * ph_r[q, r] * ph_c[q, c]


def bubble_label_to_code(chi_label, geom, nd, flags):
    """The full label -> code map of the bubble for a ``t = 0`` flag tuple
    (:func:`candidate_map` without the row/column transposition)."""
    flags = _check_flags(flags)
    if flags[6] != 0:
        raise ValueError("bubble_label_to_code: the transposition flag t must be 0, got {}"
                         .format(flags[6]))
    return candidate_map(chi_label, geom, nd, flags)


# ---------------------------------------------------------------------------
# Shared interaction tables and small helpers
# ---------------------------------------------------------------------------

#: every interaction type of the spin-orbital bond vertex
TYPES = ("CoulombIntra", "CoulombInter", "Hund", "Ising", "Exchange", "PairLift", "PairHop")
#: the types that may sit on an off-site bond
TYPES_OFFSITE = TYPES[1:]
_CANONICAL = {t.lower(): t for t in TYPES}
#: the coefficient of each type in the shared fixtures (PairHop complex);
#: next-nearest-neighbour rows carry ``_NNN`` times the nearest-neighbour one
BASE = {"CoulombIntra": 0.83, "CoulombInter": 0.37, "Hund": 0.53, "Ising": 0.71,
         "Exchange": 0.29, "PairLift": 0.61, "PairHop": 0.43 + 0.27j}
_NNN = 0.47

#: the Bravais cell of every 1-D fixture (cell coordinates); read-only
CELL = np.eye(3)
CELL.flags.writeable = False


def rel(a, b):
    """``max|a - b| / max|b|``, the max-norm relative deviation."""
    a, b = np.asarray(a), np.asarray(b)
    return float(np.abs(a - b).max() / np.abs(b).max())


def _vec3(R):
    """An integer displacement as a 3-tuple (an int ``r`` is ``(r, 0, 0)``)."""
    if isinstance(R, (int, np.integer)):
        return (int(R), 0, 0)
    out = tuple(int(x) for x in R)
    if len(out) != 3:
        raise ValueError("displacement must have three components, got {!r}".format(R))
    return out


def _neg(R):
    return tuple(-x for x in R)


def closed_row(R, a, b, v):
    """The row ``(R, a, b) = v`` and its partner ``(-R, b, a) = conj(v)``
    (0-based orbitals; an int ``R`` is ``(R, 0, 0)``). Hermitian-closed for
    every type; for a real ``v`` the partner is ``v``."""
    R = _vec3(R)
    return {(R, (a, b)): v, (_neg(R), (b, a)): np.conj(v)}


#: an off-site CoulombInter table (L = 5, norb_phys = 2), mirror-asymmetric
T1 = {**closed_row(1, 0, 1, 0.37), **closed_row(-1, 0, 1, 0.22)}
#: a Hund table on two shells (L = 5, norb_phys = 2)
T2 = {**closed_row(1, 0, 1, 0.53), **closed_row(2, 1, 0, 0.25)}
#: an on-site CoulombIntra table (norb_phys = 2)
INTRA = {((0, 0, 0), (0, 0)): 0.83, ((0, 0, 0), (1, 1)): 0.66}


def merge_tables(*tables):
    """Interaction tables ``{type: {key: v}}`` merged type by type into a new
    dict; an entry given in two tables is refused (``ValueError``)."""
    out = {}
    for tabs in tables:
        for t, rows in tabs.items():
            dst = out.setdefault(t, {})
            common = set(dst) & set(rows)
            if common:
                raise ValueError("merge_tables: {} entries {} given twice"
                                 .format(t, sorted(common)))
            dst.update(rows)
    return out


def _compiler_tables():
    """L = 5, norb_phys = 2: every type on R = +-1, +-2 with mirror-asymmetric
    inter-orbital rows (``X_01(R) != X_01(-R)``, ``X_01 != X_10``), a
    same-orbital off-site row, an on-site inter-orbital row, complex PairHop."""
    out = {"CoulombIntra": dict(INTRA)}
    for t in TYPES_OFFSITE:
        v = BASE[t]
        tbl = {}
        for row in (closed_row(1, 0, 1, v), closed_row(-1, 0, 1, 0.6 * v),
                    closed_row(2, 0, 1, _NNN * v), closed_row(-2, 0, 1, _NNN * 0.6 * v),
                    closed_row(1, 1, 1, 0.3 * v), closed_row(0, 0, 1, 0.8 * v)):
            tbl.update(row)
        out[t] = tbl
    return out


#: the compiler fixture (20 modes). Shared by several test modules: copy
#: before mutating.
COMPILER_TABLES = _compiler_tables()


# ---------------------------------------------------------------------------
# Table-driven operator builder
# ---------------------------------------------------------------------------

class Monomial(NamedTuple):
    """``v c+_1 c_2 c+_3 c_4``: spin-block modes ``X = s * norb_phys + a`` and
    integer cells (3-tuples), operator 1 in cell 0."""
    v: complex
    modes: tuple
    cells: tuple


def _canonical_tables(tables, where):
    out = {}
    for key, tbl in tables.items():
        name = _CANONICAL.get(str(key).lower())
        if name is None:
            raise ValueError("{}: unknown interaction type {!r}".format(where, key))
        if name in out:
            raise ValueError("{}: interaction type {!r} given twice".format(where, name))
        out[name] = {(_vec3(R), (int(ab[0]), int(ab[1]))): complex(v)
                     for (R, ab), v in dict(tbl).items()}
    return out


def _symmetrised(name, tbl):
    """``T~[(R,a,b)] = (T[(R,a,b)] + P(T[(-R,b,a)])) / 2`` on the literal
    keys and their partners (no wrapping), a missing entry counting 0;
    ``P = conj`` for PairHop, the identity otherwise."""
    p = np.conj if name == "PairHop" else (lambda z: z)
    keys = set(tbl) | {(_neg(R), (b, a)) for (R, (a, b)) in tbl}
    return {(R, (a, b)): 0.5 * (tbl.get((R, (a, b)), 0.0) + p(tbl.get((_neg(R), (b, a)), 0.0)))
            for (R, (a, b)) in sorted(keys)}


def expand_tables(tables, norb_phys):
    """The documented operator of every entry of every interaction table.

    Each table is first read symmetrised (``T~``, see :func:`_symmetrised`);
    then EVERY entry ``(R, a, b)`` of ``T~`` (both members of a bond) is
    expanded with the operator of the interaction-file documentation,
    including its 1/2 and its "+ h.c." as written, ``a`` in cell 0, ``b`` in
    cell ``R``, ``sigma = +1`` (up), ``-1`` (down):

    ============  ===========================================================
    CoulombIntra  ``U c+_{a up} c_{a up} c+_{a dn} c_{a dn}`` (R = 0, a = b only)
    CoulombInter  ``V/2 sum_{s,s'} n_{a s}[0] n_{b s'}[R]``
    Hund          ``-J/2 sum_s n_{a s}[0] n_{b s}[R]``
    Ising         ``J/2 sum_{s,s'} sigma_s sigma_s' n_{a s}[0] n_{b s'}[R]``
    Exchange      ``J c+_{a up}[0] c_{b up}[R] c+_{b dn}[R] c_{a dn}[0]``
    PairHop       ``J/2 c+_{a up}[0] c_{b up}[R] c+_{a dn}[0] c_{b dn}[R]
                  + J*/2 c+_{b dn}[R] c_{a dn}[0] c+_{b up}[R] c_{a up}[0]``
    PairLift      ``J/2 c+_{a up}[0] c_{a dn}[0] c+_{b up}[R] c_{b dn}[R]
                  + J*/2 c+_{b dn}[R] c_{b up}[R] c+_{a dn}[0] c_{a up}[0]``
    ============  ===========================================================

    (``n_{a s}[c]`` stands for ``c+_{a s}[c] c_{a s}[c]``; every monomial is
    the four-operator string ``c+ c c+ c`` in the order written, not normal
    ordered.) Every monomial is translated
    so that operator 1 sits in cell 0. Type keys are case-insensitive.

    Returns
    -------
    list of :class:`Monomial` (spin-block modes), in type order then
    ascending entry key.
    """
    norb = int(norb_phys)
    tabs = _canonical_tables(tables, "expand_tables")
    o = (0, 0, 0)

    def X(s, a):
        return s * norb + a
    out = []
    for name in TYPES:
        if name not in tabs:
            continue
        for (R, (a, b)), w in _symmetrised(name, tabs[name]).items():
            mR = _neg(R)
            if name == "CoulombIntra":
                if R != o or a != b:
                    raise ValueError("expand_tables: CoulombIntra entry R={} orbitals ({}, {}) "
                                     "is not on-site and same-orbital".format(R, a, b))
                out.append(Monomial(w, (X(0, a), X(0, a), X(1, a), X(1, a)), (o, o, o, o)))
            elif name in ("CoulombInter", "Ising"):
                for s, ss in ((0, 1.0), (1, -1.0)):
                    for s2, ss2 in ((0, 1.0), (1, -1.0)):
                        c = 0.5 * w * (ss * ss2 if name == "Ising" else 1.0)
                        out.append(Monomial(c, (X(s, a), X(s, a), X(s2, b), X(s2, b)),
                                            (o, o, R, R)))
            elif name == "Hund":
                for s in (0, 1):
                    out.append(Monomial(-0.5 * w, (X(s, a), X(s, a), X(s, b), X(s, b)),
                                        (o, o, R, R)))
            elif name == "Exchange":
                out.append(Monomial(w, (X(0, a), X(0, b), X(1, b), X(1, a)), (o, R, R, o)))
            elif name == "PairHop":
                out.append(Monomial(0.5 * w, (X(0, a), X(0, b), X(1, a), X(1, b)), (o, R, o, R)))
                out.append(Monomial(0.5 * np.conj(w), (X(1, b), X(1, a), X(0, b), X(0, a)),
                                    (o, mR, o, mR)))
            else:                                               # PairLift
                out.append(Monomial(0.5 * w, (X(0, a), X(1, a), X(0, b), X(1, b)), (o, o, R, R)))
                out.append(Monomial(0.5 * np.conj(w), (X(1, b), X(0, b), X(1, a), X(0, a)),
                                    (o, o, mR, mR)))
    return [Monomial(complex(m.v), m.modes, m.cells) for m in out]


def _shape(L):
    """``L`` (int, a chain) or a 3-tuple, as the mesh ``(Nx, Ny, Nz)``."""
    if isinstance(L, (int, np.integer)):
        return (int(L), 1, 1)
    shape = tuple(int(x) for x in L)
    if len(shape) != 3:
        raise ValueError("the mesh must be an int or three ints, got {!r}".format(L))
    return shape


def ed_mode(j, X, norb_phys):
    """The interleaved exact-diagonalization mode of spin-block index ``X`` in
    cell ``j`` (C-order index): ``nd * j + 2 a + s`` with
    ``X = s * norb_phys + a``. Works elementwise on integer arrays."""
    s, a = np.divmod(X, norb_phys)
    out = 2 * norb_phys * np.asarray(j) + 2 * a + s
    return int(out) if np.ndim(out) == 0 else out


def translate_monomials(modes, cells, L, norb_phys):
    """The interleaved exact-diagonalization modes ``(nvol, n)`` of one
    operator string with spin-block ``modes`` ``(n,)`` at integer ``cells``
    ``(n, 3)``, translated over every cell of the periodic mesh ``L`` (int or
    3-tuple; row ``j`` = translation by the C-order cell ``j``)."""
    shape = _shape(L)
    grid = mesh_grid(shape)
    idx = mesh_index(grid[:, None, :] + np.asarray(cells, dtype=np.int64)[None], shape)
    return ed_mode(idx, np.asarray(modes, dtype=np.int64)[None, :], norb_phys)


def build_terms(tables, L, norb_phys):
    """Quartic terms ``(p, q, r, s, v)`` = ``v c+_p c_q c+_r c_s`` of the
    tables, every monomial of :func:`expand_tables` translated over every
    cell of the periodic mesh ``L`` (int or 3-tuple), interleaved modes."""
    out = []
    for m in expand_tables(tables, norb_phys):
        modes = translate_monomials(m.modes, m.cells, L, norb_phys)               # (nvol, 4)
        out.extend((int(p), int(q), int(r), int(s), m.v) for p, q, r, s in modes)
    return out


def _sb_of_interleaved(A, norb_phys):
    a, s = divmod(int(A), 2)
    return s * norb_phys + a


def build_h1(transfer, L, norb_phys, *, zeeman=0.0):
    """The one-body matrix ``(nmode, nmode)`` of a spin-orbital Transfer table.

    ``transfer`` maps ``(R, (A, B))`` (interleaved in-cell indices
    ``A = 2 a + s``) to ``t``, the coefficient of ``c+_{A,j} c_{B,j+R}``:
    ``h1[mode(j, A), mode(j+R, B)] += t``. ``zeeman`` adds ``+zeeman`` on
    every spin-up mode and ``-zeeman`` on every spin-down one. Refuses a
    non-Hermitian result."""
    shape = _shape(L)
    nd = 2 * int(norb_phys)
    grid = mesh_grid(shape)
    nvol = len(grid)
    h = np.zeros((nd * nvol, nd * nvol), dtype=np.complex128)
    rows0 = nd * np.arange(nvol)
    for (R, (A, B)), t in dict(transfer).items():
        jj = mesh_index(grid + np.array(_vec3(R)), shape)
        h[rows0 + int(A), nd * jj + int(B)] += t
    if zeeman:
        diag = np.arange(nd * nvol)
        h[diag, diag] += np.where(diag % 2 == 0, zeeman, -zeeman)
    if not np.allclose(h, h.conj().T, rtol=0.0, atol=1e-12 * max(1.0, np.abs(h).max())):
        raise ValueError("build_h1: the transfer table is not Hermitian-closed")
    return h


def band_k(transfer, spatial_shape, norb_phys, *, zeeman=0.0):
    """``eps_XY(k) = sum_R e^{+ik.R} t_XY(R)`` ``(nvol, nd, nd)`` in the
    spin-block index, from the same table as :func:`build_h1` (and the same
    Zeeman term)."""
    shape = _shape(spatial_shape)
    nd = 2 * int(norb_phys)
    kv = mesh_momenta(shape)
    hk = np.zeros((len(kv), nd, nd), dtype=np.complex128)
    for (R, (A, B)), t in dict(transfer).items():
        X, Y = _sb_of_interleaved(A, norb_phys), _sb_of_interleaved(B, norb_phys)
        hk[:, X, Y] += t * np.exp(1j * (kv @ np.array(_vec3(R), dtype=float)))
    if zeeman:
        z = np.where(np.arange(nd) < norb_phys, zeeman, -zeeman)
        hk[:, np.arange(nd), np.arange(nd)] += z
    return hk


# ---------------------------------------------------------------------------
# Fock space and term action (bit manipulation)
# ---------------------------------------------------------------------------

_POP8 = np.unpackbits(np.arange(256, dtype=np.uint8)[:, None], axis=1).sum(axis=1).astype(np.int64)


def _popcount(x, nbits=64):
    """Set bits of ``x`` among its lowest ``nbits`` bits (x must be zero above)."""
    out = np.zeros(x.shape, dtype=np.int64)
    for shift in range(0, nbits, 8):
        out += _POP8[(x >> shift) & 0xFF]
    return out


@functools.lru_cache(maxsize=None)
def fock_states(nmode, n):
    """Sorted occupation bitstrings (int64; bit ``m`` = mode ``m``) with
    ``n`` particles among ``nmode`` modes. Read-only."""
    if nmode > 62:
        raise ValueError("fock_states: at most 62 modes, got {}".format(nmode))
    st = np.array(sorted(sum(1 << m for m in c) for c in itertools.combinations(range(nmode), n)),
                  dtype=np.int64)
    st.flags.writeable = False
    return st


def _act(states, ops):
    """Apply ``ops`` (``(mode, create)``, rightmost operator first) to every
    basis state: ``(valid, new_states, sign)`` with the Jordan-Wigner sign
    ``(-1)^(occupied modes below)`` of every step."""
    st = np.array(states, dtype=np.int64)
    sign = np.ones(st.shape, dtype=np.int64)
    valid = np.ones(st.shape, dtype=bool)
    for m, create in ops:
        bit = (st >> m) & 1
        valid &= (bit == 0) if create else (bit == 1)
        sign *= 1 - 2 * (_popcount(st & ((1 << m) - 1), m) & 1)
        st ^= (1 << m)
    return valid, st, sign


def _ops_of(modes):
    if len(modes) == 2:
        p, q = modes
        return ((q, False), (p, True))
    p, q, r, s = modes
    return ((s, False), (r, True), (q, False), (p, True))


def _term_coo(terms, states):
    """``(rows, cols, data)`` of a term list on a fixed-N basis."""
    rows, cols, data = [], [], []
    n = len(states)
    for term in terms:
        *modes, v = term
        if v == 0:
            continue
        valid, new, sign = _act(states, _ops_of(tuple(int(m) for m in modes)))
        if not valid.any():
            continue
        idx = np.searchsorted(states, new[valid])
        if (idx >= n).any() or (states[np.minimum(idx, n - 1)] != new[valid]).any():
            raise ValueError("apply_terms: a term does not conserve the particle number")
        rows.append(idx)
        cols.append(np.nonzero(valid)[0])
        data.append(v * sign[valid])
    if not rows:
        return (np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0, complex))
    return np.concatenate(rows), np.concatenate(cols), np.concatenate(data).astype(complex)


def apply_terms(terms, states, vecs):
    """``H vecs`` for a term list (quartic ``(p, q, r, s, v)`` or one-body
    ``(p, q, v)``, interleaved modes) on vectors in the fixed-N basis
    ``states`` (sorted bitstrings, :func:`fock_states`); ``vecs`` is
    ``(len(states),)`` or ``(len(states), k)``."""
    states = np.asarray(states, dtype=np.int64)
    vecs = np.asarray(vecs, dtype=np.complex128)
    flat = vecs.ndim == 1
    v2 = vecs[:, None] if flat else vecs
    rows, cols, data = _term_coo(terms, states)
    out = np.zeros_like(v2)
    np.add.at(out, rows, data[:, None] * v2[cols])
    return out[:, 0] if flat else out


def _terms_block(terms, states):
    rows, cols, data = _term_coo(terms, states)
    n = len(states)
    return sp.csr_matrix((data, (rows, cols)), shape=(n, n))


@functools.lru_cache(maxsize=None)
def _onebody_coo(nmode, n):
    """``(rows, cols, pair, sign)``: ``c+_p c_q`` (``pair = p * nmode + q``)
    maps basis state ``cols`` to ``sign`` times state ``rows``."""
    states = fock_states(nmode, n)
    out = [[], [], [], []]
    for p in range(nmode):
        for q in range(nmode):
            valid, new, sign = _act(states, ((q, False), (p, True)))
            if not valid.any():
                continue
            out[0].append(np.searchsorted(states, new[valid]))
            out[1].append(np.nonzero(valid)[0])
            out[2].append(np.full(int(valid.sum()), p * nmode + q))
            out[3].append(sign[valid])
    if not out[0]:
        return tuple(np.zeros(0, np.int64) for _ in range(4))
    return tuple(np.concatenate(x) for x in out)


def _onebody_dense(Ms, nmode, n):
    """The fixed-N blocks ``(k, dim, dim)`` of ``sum_pq M_pq c+_p c_q`` for a
    stack of matrices ``Ms`` ``(k, nmode, nmode)`` (dense)."""
    rows, cols, pair, sign = _onebody_coo(nmode, n)
    dim = len(fock_states(nmode, n))
    out = np.zeros((len(Ms), dim * dim), dtype=np.complex128)
    flat = rows * dim + cols
    vals = Ms.reshape(len(Ms), -1)[:, pair] * sign
    for s in range(len(Ms)):
        out[s] = np.bincount(flat, weights=vals[s].real, minlength=dim * dim) \
            + 1j * np.bincount(flat, weights=vals[s].imag, minlength=dim * dim)
    return out.reshape(len(Ms), dim, dim)


def _onebody_block(M, nmode, n):
    """The fixed-N block of ``sum_pq M_pq c+_p c_q`` (sparse)."""
    rows, cols, pair, sign = _onebody_coo(nmode, n)
    data = np.asarray(M, dtype=np.complex128).reshape(-1)[pair] * sign
    nz = data != 0
    dim = len(fock_states(nmode, n))
    return sp.csr_matrix((data[nz], (rows[nz], cols[nz])), shape=(dim, dim))


# ---------------------------------------------------------------------------
# Exact diagonalization (Jordan-Wigner, blocked by total N)
# ---------------------------------------------------------------------------

def _lehmann_kernel(E, p, beta):
    """``K[m, n] = (p_m - p_n) / (E_n - E_m)``, and ``beta p_m`` where
    ``|E_n - E_m| <= 1e-10`` (the static Lehmann kernel of one block)."""
    dE = E[None, :] - E[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(np.abs(dE) > 1e-10, (p[:, None] - p[None, :]) / dE, beta * p[:, None])


class JWED:
    """Grand-canonical exact diagonalization of ``h1 + terms - mu N``.

    Sparse Jordan-Wigner blocks, one per total particle number, each
    diagonalized densely; the Boltzmann weights are normalized over the
    whole Fock space.

    Parameters
    ----------
    h1 : ndarray (nmode, nmode)
        One-body matrix, interleaved modes (:func:`build_h1`).
    terms : list
        Quartic terms ``(p, q, r, s, v)`` (:func:`build_terms`).
    mu, beta : float
    """

    def __init__(self, h1, terms, *, mu, beta):
        h1 = np.asarray(h1, dtype=np.complex128)
        nmode = h1.shape[0]
        hams = [_onebody_block(h1, nmode, n) + _terms_block(terms, fock_states(nmode, n))
                for n in range(nmode + 1)]
        self._setup(nmode, hams, mu, beta)

    @classmethod
    def from_blocks(cls, nmode, hams, *, mu, beta):
        """From the number-conserving blocks ``hams[N]`` (without ``-mu N``)."""
        self = cls.__new__(cls)
        self._setup(nmode, hams, mu, beta)
        return self

    def _setup(self, nmode, hams, mu, beta):
        self.nmode, self.mu, self.beta = int(nmode), float(mu), float(beta)
        eig = []
        for n, H in enumerate(hams):
            Hd = H.toarray() if sp.issparse(H) else np.array(H, dtype=np.complex128)
            if not np.allclose(Hd, Hd.conj().T, rtol=0.0,
                               atol=1e-12 * max(1.0, np.abs(Hd).max())):
                raise ValueError("JWED: the Hamiltonian block N={} is not Hermitian".format(n))
            E, V = np.linalg.eigh(Hd - self.mu * n * np.eye(len(Hd)))
            eig.append((n, E, V))
        e0 = min(E.min() for _, E, _ in eig)
        z = sum(np.exp(-self.beta * (E - e0)).sum() for _, E, _ in eig)
        self.blocks = [(n, E - e0, V, np.exp(-self.beta * (E - e0)) / z) for n, E, V in eig]

    def static_corr(self, A, B):
        """``int_0^beta <A(tau) B> dtau - beta <A><B>`` for one-body operators
        ``A = sum A_pq c+_p c_q``, ``B`` likewise (``(nmode, nmode)`` matrices)."""
        return complex(self.static_corr_many([(A, B)])[0])

    def static_corr_many(self, pairs):
        """:meth:`static_corr` of every ``(A, B)`` in ``pairs`` (one kernel and
        one batched pair of dense products per block)."""
        As = np.array([a for a, _ in pairs], dtype=np.complex128)
        Bs = np.array([b for _, b in pairs], dtype=np.complex128)
        total = np.zeros(len(pairs), dtype=np.complex128)
        ea = np.zeros(len(pairs), dtype=np.complex128)
        eb = np.zeros(len(pairs), dtype=np.complex128)
        for n, E, V, p in self.blocks:
            K = _lehmann_kernel(E, p, self.beta)
            # bounded batches: at most ~2e6 complex entries per stacked block
            step = max(1, 2_000_000 // (len(E) * len(E)))
            for lo in range(0, len(pairs), step):
                sl = slice(lo, lo + step)
                Ae = V.conj().T @ _onebody_dense(As[sl], self.nmode, n) @ V
                Be = V.conj().T @ _onebody_dense(Bs[sl], self.nmode, n) @ V
                total[sl] += np.einsum("smn,snm,mn->s", Ae, Be, K)
                ea[sl] += np.einsum("smm,m->s", Ae, p)
                eb[sl] += np.einsum("smm,m->s", Be, p)
        return total - self.beta * ea * eb

    def static_corr_matrix(self, As, Bs):
        """``out[i, j]`` = :meth:`static_corr` of ``(As[i], Bs[j])``, each
        operator rotated into the eigenbasis once."""
        As = np.asarray(As, dtype=np.complex128)
        Bs = np.asarray(Bs, dtype=np.complex128)
        total = np.zeros((len(As), len(Bs)), dtype=np.complex128)
        ea = np.zeros(len(As), dtype=np.complex128)
        eb = np.zeros(len(Bs), dtype=np.complex128)
        for n, E, V, p in self.blocks:
            K = _lehmann_kernel(E, p, self.beta)
            Ae = V.conj().T @ _onebody_dense(As, self.nmode, n) @ V
            Be = V.conj().T @ _onebody_dense(Bs, self.nmode, n) @ V
            total += np.einsum("imn,jnm->ij", Ae * K, Be, optimize=True)
            ea += np.einsum("smm,m->s", Ae, p)
            eb += np.einsum("smm,m->s", Be, p)
        return total - self.beta * ea[:, None] * eb[None, :]

    def density(self):
        """``n[a, b] = <c+_a c_b>`` ``(nmode, nmode)``."""
        out = np.zeros(self.nmode * self.nmode, dtype=np.complex128)
        for n, E, V, p in self.blocks:
            rows, cols, pair, sign = _onebody_coo(self.nmode, n)
            rho = (V * p) @ V.conj().T
            np.add.at(out, pair, sign * rho[cols, rows])
        return out.reshape(self.nmode, self.nmode)


# ---------------------------------------------------------------------------
# One-body (free) references and the Hartree-Fock shift
# ---------------------------------------------------------------------------

def _fermi(x, beta):
    return 0.5 * (1.0 - np.tanh(0.5 * beta * x))


def _pair_kernel(xa, xb, beta):
    """``(f(xb) - f(xa)) / (xa - xb)``, and ``-f'(xa)`` where
    ``|xa - xb| < 1e-10``."""
    fa, fb = _fermi(xa, beta), _fermi(xb, beta)
    d = xa - xb
    deg = np.abs(d) < 1e-10
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(deg, beta * fa * (1.0 - fa), (fb - fa) / np.where(deg, 1.0, d))


def free_static_corr(h, mu, beta, A, B):
    """The free (one-body ``h``) value of :meth:`JWED.static_corr`."""
    e, U = np.linalg.eigh(h)
    xi = e - mu
    F = _pair_kernel(xi[:, None], xi[None, :], beta)
    At = U.conj().T @ A @ U
    Bt = U.conj().T @ B @ U
    return complex(np.sum(At * Bt.T * F))


def free_density(h, mu, beta):
    """``n[a, b] = <c+_a c_b>`` of the one-body ``h``."""
    e, U = np.linalg.eigh(h)
    return ((U * _fermi(e - mu, beta)) @ U.conj().T).T


def hf_onebody(terms, n):
    """The Hartree-Fock one-body matrix ``dE/dn`` of a quartic term list at
    the density ``n[a, b] = <c+_a c_b>``, with
    ``E = sum v [n_pq n_rs + n_ps (delta_qr - n_rq)]`` (Wick's theorem for
    ``v c+_p c_q c+_r c_s``)."""
    n = np.asarray(n)
    h = np.zeros(n.shape, dtype=np.complex128)
    for (p, q, r, s, v) in terms:
        h[p, q] += v * n[r, s]
        h[r, s] += v * n[p, q]
        h[p, s] += v * ((1.0 if q == r else 0.0) - n[r, q])
        h[r, q] -= v * n[p, s]
    return h


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def freeze_table(table):
    """A table (mapping, or iterable of ``(key, value)``) as a sorted tuple of
    ``(R, (a, b), complex v)``; refuses a key given twice."""
    items = table.items() if hasattr(table, "items") else table
    out = {}
    for (R, ab), v in items:
        key = (_vec3(R), (int(ab[0]), int(ab[1])))
        if key in out:
            raise ValueError("freeze_table: entry R={} orbitals ({}, {}) given twice"
                             .format(key[0], key[1][0], key[1][1]))
        out[key] = complex(v)
    return tuple((R, ab, v) for (R, ab), v in sorted(out.items()))


def freeze_tables(tables):
    """Interaction tables as a sorted tuple of ``(type, frozen table)``,
    canonical type names."""
    tabs = _canonical_tables(dict(tables), "freeze_tables")
    return tuple((t, freeze_table(tabs[t])) for t in sorted(tabs))


class OracleGeometry:
    """A duck-typed bond geometry: ``delta_r`` (B, 3), ``reverse`` (B,),
    ``spatial_shape`` (read-only arrays)."""
    __slots__ = ("delta_r", "reverse", "spatial_shape")

    def __init__(self, delta_r, reverse, spatial_shape):
        dr = np.array(delta_r, dtype=np.int64).reshape(-1, 3)
        rev = np.array(reverse, dtype=np.int64)
        dr.flags.writeable = False
        rev.flags.writeable = False
        self.delta_r, self.reverse, self.spatial_shape = dr, rev, tuple(spatial_shape)

    @property
    def B(self):
        return int(self.delta_r.shape[0])

    n_channels = B

    @property
    def nvol(self):
        return int(np.prod(self.spatial_shape))


class FixtureSpec(NamedTuple):
    """An exact-diagonalization fixture, frozen and hashable by content.

    ``L`` is an int (a chain) or a 3-tuple mesh; ``transfer`` and ``tables``
    are frozen (:func:`freeze_table`, :func:`freeze_tables`); ``zeeman`` is
    ``+zeeman`` on spin up, ``-zeeman`` on spin down."""
    L: object
    norb_phys: int
    beta: float
    mu: float
    transfer: tuple
    tables: tuple
    zeeman: float = 0.0

    @classmethod
    def make(cls, L, norb_phys, beta, mu, transfer, tables, *, zeeman=0.0):
        return cls(L if isinstance(L, (int, np.integer)) else tuple(int(x) for x in L),
                   int(norb_phys), float(beta), float(mu), freeze_table(transfer),
                   freeze_tables(tables), float(zeeman))

    @property
    def spatial_shape(self):
        return _shape(self.L)

    @property
    def nvol(self):
        return int(np.prod(self.spatial_shape))

    @property
    def nd(self):
        return 2 * self.norb_phys

    @property
    def transfer_dict(self):
        return {(R, ab): v for R, ab, v in self.transfer}

    @property
    def tables_dict(self):
        return {t: {(R, ab): v for R, ab, v in tbl} for t, tbl in self.tables}

    @property
    def h1(self):
        """The one-body matrix (interleaved modes)."""
        return build_h1(self.transfer_dict, self.L, self.norb_phys, zeeman=self.zeeman)

    @property
    def band(self):
        """``(h1_k, mu, beta)``, the arguments of :func:`lindhard_static`."""
        return (band_k(self.transfer_dict, self.spatial_shape, self.norb_phys,
                       zeeman=self.zeeman), self.mu, self.beta)

    @property
    def geom(self):
        return geom_of(self)


@functools.lru_cache(maxsize=None)
def geom_of(fixture):
    """The label geometry of a fixture: ``0`` and every declared
    displacement (interaction and transfer) with its negative, ordered by
    length then lexicographically (channel 0 is on-site), one representative
    per class modulo the mesh; ``reverse[m]`` is the channel of ``-delta_r[m]``."""
    shape = fixture.spatial_shape
    ds = {(0, 0, 0)}
    for _, tbl in fixture.tables:
        for R, _, _ in tbl:
            ds.update({R, _neg(R)})
    for R, _, _ in fixture.transfer:
        ds.update({R, _neg(R)})
    keep, seen = [], {}
    for d in sorted(ds, key=lambda d: (sum(x * x for x in d), d)):
        key = int(mesh_index(np.array(d), shape))
        if key not in seen:
            seen[key] = len(keep)
            keep.append(d)
    rev = [seen[int(mesh_index(-np.array(d), shape))] for d in keep]
    return OracleGeometry(keep, rev, shape)


def label_list(geom, nd):
    """Labels ``(d, X, Y)`` in flat index order ``m * nd**2 + X * nd + Y``."""
    return [(tuple(int(x) for x in d), X, Y) for d in np.asarray(geom.delta_r)
            for X in range(nd) for Y in range(nd)]


def _fixture_i_transfer(*, spin_flip=True):
    """Fixture (i) band (norb_phys = 1): the on-site Hermitian spin matrix
    and a random complex spin matrix (seed 7, times 1/2) at R = +1 with its
    partner at R = -1; ``spin_flip=False`` drops every spin-off-diagonal
    entry."""
    rng = np.random.default_rng(7)
    t1 = 0.5 * (rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))
    h0 = np.array([[0.15, 0.07 - 0.05j], [0.07 + 0.05j, -0.1]])
    tr = {}
    for s in (0, 1):
        for s2 in (0, 1):
            if not spin_flip and s != s2:
                continue
            tr[((0, 0, 0), (s, s2))] = h0[s, s2]
            tr[((1, 0, 0), (s, s2))] = t1[s, s2]
            tr[((-1, 0, 0), (s2, s))] = np.conj(t1[s, s2])
    return tr


def _fixture_i_tables(*, pairlift=True):
    """Fixture (i) interactions: CoulombIntra on site, every other type on
    R = +-1 and (times 0.47) R = +-2, Hermitian-closed."""
    out = {"CoulombIntra": {((0, 0, 0), (0, 0)): BASE["CoulombIntra"]}}
    for t in TYPES_OFFSITE:
        if t == "PairLift" and not pairlift:
            continue
        out[t] = {**closed_row(1, 0, 0, BASE[t]), **closed_row(2, 0, 0, _NNN * BASE[t])}
    return out


#: Transfer rows ``(r, a, b)``: ``t c+_{j,a,s} c_{j+r,b,s}`` of the two-orbital
#: chain (the values of ``tests/test_so_pairing_ed.py``)
TRANSFER_II = {(0, 0, 0): -0.4, (0, 1, 1): 0.6, (0, 0, 1): -0.3, (0, 1, 0): -0.3,
               (1, 0, 1): -1.0, (-1, 1, 0): -1.0, (1, 0, 0): -0.5, (-1, 0, 0): -0.5}
#: spin-flip hopping ``(r, a, b)``: ``v c+_{j,a,up} c_{j+r,b,dn} + h.c.``
SOC_II = {(1, 0, 1): 0.7 + 0.5j, (1, 0, 0): -0.4 + 0.6j, (0, 0, 1): 0.6j,
          (1, 1, 1): 0.3 - 0.6j}


def _fixture_ii_transfer():
    tr = {}

    def add(key, v):
        tr[key] = tr.get(key, 0.0) + v
    for (r, a, b), t in TRANSFER_II.items():
        for s in (0, 1):
            add(((r, 0, 0), (2 * a + s, 2 * b + s)), t)
    for (r, a, b), v in SOC_II.items():
        add(((r, 0, 0), (2 * a, 2 * b + 1)), v)
        add(((-r, 0, 0), (2 * b + 1, 2 * a)), np.conj(v))
    return tr


def _fixture_ii_tables():
    """Fixture (ii) interactions (norb_phys = 2): CoulombIntra on site, every
    other type on the NN inter-orbital rows ``X_01(+1) = v``,
    ``X_01(-1) = 0.6 v`` with their partners (complex PairHop)."""
    out = {"CoulombIntra": {((0, 0, 0), (0, 0)): BASE["CoulombIntra"],
                            ((0, 0, 0), (1, 1)): 0.8 * BASE["CoulombIntra"]}}
    for t in TYPES_OFFSITE:
        out[t] = {**closed_row(1, 0, 1, BASE[t]), **closed_row(-1, 0, 1, 0.6 * BASE[t])}
    return out


#: (i): L = 5, norb_phys = 1 (10 modes), every type on NN and NNN bonds,
#: complex spin-flip hopping
FIXTURE_I = FixtureSpec.make(5, 1, 2.0, 0.3, _fixture_i_transfer(), _fixture_i_tables())
#: (i) with a Zeeman field
FIXTURE_I_ZEEMAN = FIXTURE_I._replace(zeeman=0.15)
#: (i) with S_z conserved: no spin-flip hopping, no PairLift, Zeeman on
FIXTURE_I_SZ = FixtureSpec.make(5, 1, 2.0, 0.3, _fixture_i_transfer(spin_flip=False),
                                _fixture_i_tables(pairlift=False), zeeman=0.15)
#: (ii): L = 3, norb_phys = 2 (12 modes, largest sector 924)
FIXTURE_II = FixtureSpec.make(3, 2, 2.0, 0.8, _fixture_ii_transfer(), _fixture_ii_tables())

FIXTURES = {"i": FIXTURE_I, "ii": FIXTURE_II, "i_zeeman": FIXTURE_I_ZEEMAN,
            "i_sz": FIXTURE_I_SZ}
#: the central-difference step of each fixture's first-order vertex part
FIXTURE_X = {"i": 2e-4, "ii": 2e-4, "i_zeeman": 2e-4, "i_sz": 2e-4}
#: Richardson extrapolation is used where the measured truncation
#: ``|d(x) - d_richardson| / max|d_richardson|`` exceeds this guard (half the
#: 1e-6 first-order gate tolerance)
TRUNCATION_GUARD = 5e-7
#: whether each fixture's measured truncation exceeds the guard. Measured
#: with :func:`truncation` at the pinned step, q in {1, 2}: (i) 6.1e-8,
#: (i)-Zeeman 5.7e-8, (i)-Sz 5.0e-8 (fast sketches), (ii) 2.6e-8 (full
#: sketches)
NEEDS_RICHARDSON = {"i": False, "ii": False, "i_zeeman": False, "i_sz": False}


def fixture_name(fixture):
    """The name of a pinned fixture, or None."""
    return next((k for k, f in FIXTURES.items() if f == fixture), None)


# ---------------------------------------------------------------------------
# Bond bilinears, sketches and the exact first-order vertex part
# ---------------------------------------------------------------------------

def rho_matrix(label, q, L, norb_phys):
    """The one-body matrix of ``rho^{(d)}_{XY}(q) = sum_j e^{+iq.j}
    c+_{X,j} c_{Y,j+d}``, label ``(d, X, Y)`` (spin-block ``X``, ``Y``),
    ``q`` the C-order mesh index; interleaved modes."""
    shape = _shape(L)
    grid = mesh_grid(shape)
    d, X, Y = label
    qv = mesh_momenta(shape)[q]
    nmode = 2 * norb_phys * len(grid)
    M = np.zeros((nmode, nmode), dtype=np.complex128)
    jd = mesh_index(grid + np.array(_vec3(d)), shape)
    M[ed_mode(np.arange(len(grid)), X, norb_phys), ed_mode(jd, Y, norb_phys)] = \
        np.exp(1j * (grid @ qv))
    return M


@functools.lru_cache(maxsize=64)
def _rho_stack(shape, norb_phys, labels, q):
    out = np.array([rho_matrix(lab, q, shape, norb_phys) for lab in labels])
    out.flags.writeable = False
    return out


def _sketch_ops(u, v, q, labels, L, norb_phys):
    """``A_u = sum_I u_I rho_I(q)`` and ``B_v = sum_J v_J rho_J(-q)``."""
    shape = _shape(L)
    labels = tuple(labels)
    mq = int(mesh_negate(shape)[q])
    A = np.tensordot(np.asarray(u), _rho_stack(shape, norb_phys, labels, int(q)), 1)
    B = np.tensordot(np.asarray(v), _rho_stack(shape, norb_phys, labels, mq), 1)
    return A, B


def sketch_value(ed, u, v, q, labels, L, norb_phys):
    """``u^T chi(q) v`` of the exact static bond correlator
    ``chi_IJ(q) = (1/N) int <rho_I(q, tau) rho_J(-q)>_conn`` (no
    conjugation of ``u`` or ``v``)."""
    A, B = _sketch_ops(u, v, q, labels, L, norb_phys)
    return ed.static_corr(A, B) / int(np.prod(_shape(L)))


def label_sectors(fixture):
    """The S_z change ``sigma = s_Y - s_X`` (``s`` = 0 up, 1 down) of every
    label ``rho_XY`` of ``geom_of(fixture)``, in flat label order."""
    norb = fixture.norb_phys
    return np.array([Y // norb - X // norb for (_, X, Y) in label_list(fixture.geom, fixture.nd)])


def has_pairlift(fixture):
    """Whether the fixture has a nonzero PairLift entry (PairLift changes
    S_z by 2, so it couples label sector ``sigma`` with ``sigma``)."""
    return any(abs(v) > 0 for t, tbl in fixture.tables if t == "PairLift" for _, _, v in tbl)


def sketch_set(fixture, kind, seed=270):
    """Sketch pairs ``(u, v)`` over the labels of ``geom_of(fixture)``.

    * 8 random complex pairs over every label;
    * 3 sector-targeted pairs: ``u`` on the labels with S_z change
      ``sigma = s_Y - s_X`` of ``rho_XY`` (``s`` = 0 up, 1 down) for
      ``sigma = 0, +1, -1``, ``v`` on ``-sigma``; ``v`` on ``+sigma``
      instead when the fixture has a nonzero PairLift entry (PairLift
      changes S_z by 2, so it couples sector ``sigma`` with ``sigma``);
    * ``kind="full"`` adds: ``u`` on bond ``+R`` with ``v`` on ``-R``,
      ``u`` on ``-R`` with ``v`` on ``+R`` (``+R`` the first positive
      off-site channel), and ``u`` on the mode pair ``(A, C)`` with ``v`` on
      ``(C, A)``, ``A = 0``, ``C = nd - 1``.
    """
    if kind not in ("fast", "full"):
        raise ValueError("sketch_set: kind must be 'fast' or 'full', got {!r}".format(kind))
    geom, nd = fixture.geom, fixture.nd
    labels = label_list(geom, nd)
    NL = len(labels)
    rng = np.random.default_rng(seed)

    def rand(mask):
        z = rng.normal(size=NL) + 1j * rng.normal(size=NL)
        return np.where(mask, z, 0.0)
    every = np.ones(NL, dtype=bool)
    out = [(rand(every), rand(every)) for _ in range(8)]
    sector = label_sectors(fixture)
    lift = has_pairlift(fixture)
    for s in (0, +1, -1):
        out.append((rand(sector == s), rand(sector == (s if lift else -s))))
    if kind == "full":
        dr = np.asarray(geom.delta_r)
        mp = next(m for m in range(len(dr)) if tuple(dr[m]) > (0, 0, 0))
        mm = int(geom.reverse[mp])
        chan = np.repeat(np.arange(len(dr)), nd * nd)
        out.append((rand(chan == mp), rand(chan == mm)))
        out.append((rand(chan == mm), rand(chan == mp)))
        pairs = np.array([(X, Y) for (_, X, Y) in labels])
        A, C = 0, nd - 1
        out.append((rand((pairs[:, 0] == A) & (pairs[:, 1] == C)),
                    rand((pairs[:, 0] == C) & (pairs[:, 1] == A))))
    return out


#: the label-sector pairs ``(sigma_u, sigma_v)`` of each bond-correlator
#: family: the longitudinal ``B`` against its conjugate leg, and the
#: transverse ``S+`` against its conjugate leg (``S-``) and the reverse
FAMILIES = {"longitudinal": ((0, 0),), "transverse": ((+1, -1), (-1, +1))}


def family_sketches(fixture, family, n_random, seed=271):
    """``n_random`` random complex sketch pairs ``(u, v)`` within one
    correlator family of :data:`FAMILIES`: ``u`` on the labels of sector
    ``sigma_u`` and ``v`` on ``sigma_v`` (:func:`label_sectors`), cycling
    through the family's sector pairs."""
    if family not in FAMILIES:
        raise ValueError("family_sketches: family must be one of {}, got {!r}"
                         .format(sorted(FAMILIES), family))
    sector = label_sectors(fixture)
    rng = np.random.default_rng((seed, sorted(FAMILIES).index(family)))
    pairs = FAMILIES[family]

    def rand(s):
        z = rng.normal(size=len(sector)) + 1j * rng.normal(size=len(sector))
        return np.where(sector == s, z, 0.0)
    return [tuple(rand(s) for s in pairs[i % len(pairs)]) for i in range(int(n_random))]


def bond_observable(fixture, q, R, a, b, s_a, s_b):
    """The label vectors ``(u, v)`` of the bond bilinear
    ``O(q) = sum_j e^{+iq.j} c+_{j+R, a s_a} c_{j, b s_b}`` and of its
    conjugate leg ``O(q)^dagger``, as sketch operands
    (``u^T rho(q)``, ``v^T rho(-q)``):

    ``O(q) = e^{-iq.R} rho^{(-R)}_{XY}(q)`` and
    ``O(q)^dagger = rho^{(R)}_{YX}(-q)``, ``X = s_a norb_phys + a``,
    ``Y = s_b norb_phys + b``; ``R`` is matched to the label displacement
    equal to it modulo the mesh. ``s_a = s_b`` is the longitudinal
    ``B_{R,a,b,sigma}``; ``(s_a, s_b) = (0, 1)`` the transverse
    ``S+_{R,a,b}``. So ``u^T chi(q) v`` is the static correlator
    ``(1/N) int <O(q, tau) O(q)^dagger>_conn``."""
    geom, nd, norb = fixture.geom, fixture.nd, fixture.norb_phys
    shape = fixture.spatial_shape
    lab = bond_label_index(geom, nd)
    rep = {int(mesh_index(np.asarray(d), shape)): tuple(int(x) for x in d) for d in geom.delta_r}
    R = np.array(_vec3(R))
    dm, dp = rep[int(mesh_index(-R, shape))], rep[int(mesh_index(R, shape))]
    X, Y = s_a * norb + a, s_b * norb + b
    u = np.zeros(len(lab), dtype=np.complex128)
    v = np.zeros(len(lab), dtype=np.complex128)
    u[lab[(dm, X, Y)]] = np.exp(-1j * (mesh_momenta(shape)[q] @ R))
    v[lab[(dp, Y, X)]] = 1.0
    return u, v


def sketch_matrix(ed, us, vs, q, labels, L, norb_phys):
    """``out[i, j]`` = :func:`sketch_value` of ``(us[i], vs[j])``: every
    operator is built once (:meth:`JWED.static_corr_matrix`)."""
    shape = _shape(L)
    labels = tuple(labels)
    mq = int(mesh_negate(shape)[q])
    As = np.tensordot(np.asarray(us), _rho_stack(shape, norb_phys, labels, int(q)), 1)
    Bs = np.tensordot(np.asarray(vs), _rho_stack(shape, norb_phys, labels, mq), 1)
    return ed.static_corr_matrix(As, Bs) / int(np.prod(shape))


@functools.lru_cache(maxsize=None)
def _fixture_parts(fixture):
    """``(h1, terms, h1 blocks, interaction blocks, Hartree-Fock one-body)``."""
    h1 = fixture.h1
    terms = build_terms(fixture.tables_dict, fixture.L, fixture.norb_phys)
    nmode = h1.shape[0]
    h1b = tuple(_onebody_block(h1, nmode, n) for n in range(nmode + 1))
    hib = tuple(_terms_block(terms, fock_states(nmode, n)) for n in range(nmode + 1))
    hhf = hf_onebody(terms, free_density(h1, fixture.mu, fixture.beta))
    return h1, terms, h1b, hib, hhf


# One entry holds the eigenpairs of every particle-number block,
# sum_N C(nmode, N)^2 = C(2 nmode, nmode) complex entries: about 3 MB at 10
# modes (fixture i) and 43 MB at 12 modes (fixture ii), so 16 entries stay
# below 0.7 GB.
@functools.lru_cache(maxsize=16)
def ed_cached(fixture, x):
    """The :class:`JWED` of ``h1 + x H_int`` of a fixture, cached in this
    process by the full fixture content and ``x`` (never on disk)."""
    h1, _, h1b, hib, _ = _fixture_parts(fixture)
    hams = [a + float(x) * b for a, b in zip(h1b, hib)]
    return JWED.from_blocks(h1.shape[0], hams, mu=fixture.mu, beta=fixture.beta)


def _step_of(fixture):
    name = fixture_name(fixture)
    if name is None:
        raise ValueError("vertex_part_sketches: no pinned step for this fixture; pass x=")
    return FIXTURE_X[name]


def _central_difference(fixture, sketches, qs, x):
    h1, _, _, _, hhf = _fixture_parts(fixture)
    labels = label_list(fixture.geom, fixture.nd)
    nvol = fixture.nvol
    ops = [[_sketch_ops(u, v, q, labels, fixture.L, fixture.norb_phys) for (u, v) in sketches]
           for q in qs]

    def full(y):
        ed = ed_cached(fixture, y)
        return np.array([ed.static_corr_many(row) for row in ops]) / nvol

    def hf(y):
        h = h1 + y * hhf
        return np.array([[free_static_corr(h, fixture.mu, fixture.beta, A, B) for (A, B) in row]
                         for row in ops]) / nvol
    return (full(x) - full(-x) - hf(x) + hf(-x)) / (2.0 * x)


@functools.lru_cache(maxsize=64)
def _central_difference_kind(fixture, kind, qs, x):
    """:func:`_central_difference` on ``sketch_set(fixture, kind)``, cached in
    this process by content (the default seed). Read-only."""
    out = _central_difference(fixture, sketch_set(fixture, kind), qs, x)
    out.flags.writeable = False
    return out


def vertex_part_sketches(fixture, *, sketches, qs, richardson=False, x=None):
    """The first-order vertex part of the exact bond correlators, sketched.

    ``d(x) = [S(+x) - S(-x)] / 2x - [S_HF(+x) - S_HF(-x)] / 2x`` where
    ``S(y)`` is :func:`sketch_value` of the exact diagonalization of
    ``h1 + y H_int`` and ``S_HF(y)`` the same correlator of the free
    one-body ``h1 + y Sigma_HF[n0]`` (``n0`` the free density; evaluated by
    the exact one-body formula, no diagonalization of the Fock space).
    ``richardson=True`` returns ``(4 d(x/2) - d(x)) / 3``.

    Parameters
    ----------
    sketches : list of ``(u, v)``, or a :func:`sketch_set` kind (then the
        result is cached in this process by the fixture content)
    qs : C-order mesh indices
    x : float, optional
        The step; default the fixture's pinned :data:`FIXTURE_X` entry.

    Returns
    -------
    ndarray complex ``(len(qs), len(sketches))``.
    """
    x = _step_of(fixture) if x is None else float(x)
    qs = tuple(int(q) for q in qs)
    if isinstance(sketches, str):
        def cd(step):
            return _central_difference_kind(fixture, sketches, qs, step).copy()
    else:
        def cd(step):
            return _central_difference(fixture, sketches, qs, step)
    d = cd(x)
    if richardson:
        return (4.0 * cd(0.5 * x) - d) / 3.0
    return d


def truncation(fixture, *, qs=(1, 2), sketches="fast", x=None):
    """The measured central-difference truncation
    ``max|d(x) - d_richardson| / max|d_richardson|`` of
    :func:`vertex_part_sketches`."""
    d = vertex_part_sketches(fixture, sketches=sketches, qs=qs, x=x)
    r = vertex_part_sketches(fixture, sketches=sketches, qs=qs, x=x, richardson=True)
    return rel(d, r)


# ---------------------------------------------------------------------------
# Label-form Lindhard function, filing rule and the code-slot map
# ---------------------------------------------------------------------------

def lindhard_static(h1_k, mu, beta, geom):
    """The static label bubble from the band, analytically:

    ``chi[q, (d;X,Y), (d';Z,W)] = (1/N) sum_k e^{+ik.d} e^{+i(k+q).d'}
    sum_{a,b} U_Ya(k) U*_Za(k) U_Wb(k+q) U*_Xb(k+q) F(xi_a(k), xi_b(k+q))``

    with ``F(xa, xb) = (f(xb) - f(xa)) / (xa - xb)`` and ``-f'(xa)`` where
    ``|xa - xb| < 1e-10`` (the Matsubara sum of
    ``-T sum_n G_YZ(k, i w_n) G_WX(k+q, i w_n)``), ``xi = eps - mu``.

    Parameters
    ----------
    h1_k : ndarray (nvol, nd, nd), spin-block index
    geom : labels from ``delta_r``; mesh from ``spatial_shape``

    Returns
    -------
    ndarray complex128 (nvol, NL, NL), ``NL = B * nd**2``.
    """
    h1_k = np.asarray(h1_k, dtype=np.complex128)
    nvol, nd, _ = h1_k.shape
    shape = tuple(int(x) for x in geom.spatial_shape)
    if int(np.prod(shape)) != nvol:
        raise ValueError("lindhard_static: h1_k has nvol={} but the mesh {} has {}"
                         .format(nvol, shape, int(np.prod(shape))))
    D = np.asarray(geom.delta_r, dtype=np.int64)
    B = D.shape[0]
    NL = B * nd * nd
    e, U = np.linalg.eigh(h1_k)
    xi = e - mu
    pk = np.exp(1j * (mesh_momenta(shape) @ D.T.astype(float)))       # (k, m)
    add = mesh_add(shape)
    out = np.empty((nvol, NL, NL), dtype=np.complex128)
    for q in range(nvol):
        kq = add[q]
        F = _pair_kernel(xi[:, :, None], xi[kq][:, None, :], beta)     # (k, a, b)
        T = np.einsum("kya,kza,kab,kwb,kxb->kyzwx", U, U.conj(), F, U[kq], U[kq].conj(),
                      optimize=True)
        chi = np.einsum("km,kn,kyzwx->mxynzw", pk, pk[kq], T, optimize=True) / nvol
        out[q] = chi.reshape(NL, NL)
    return out


def file_gamma_label(monomials, geom, nd, *, q_sign=1, direct=True, crossed=True,
                     crossed_sign=-1.0):
    """The label-form interaction ``Gamma(q)`` ``(nvol, NL, NL)`` of a list of
    monomials. Either part filed alone gives ``H_int = (1/2N) sum_q sum_IJ
    Gamma_IJ(q) rho_I(q) rho_J(-q)``; the default direct plus crossed filing
    gives it with ``1/4N`` (one-body remainders are not part of Gamma).

    A monomial ``v c+_1(0) c_2(d1) c+_3(t) c_4(t+d2)`` is filed as

    * direct: ``Gamma[(d1;1,2),(d2;3,4)](q) += v e^{+iq.t}`` and
      ``Gamma[(d2;3,4),(d1;1,2)](q) += v e^{-iq.t}``;
    * crossed: ``Gamma[(t+d2;1,4),(d1-t;3,2)](q) += -v e^{+iq.t}`` and
      ``Gamma[(d1-t;3,2),(t+d2;1,4)](q) += -v e^{-iq.t}``,

    from ``c+_1 c_2 c+_3 c_4 = -c+_1 c_4 c+_3 c_2 + delta_23 c+_1 c_4 +
    delta_34 c+_1 c_2``. The keyword arguments select wrong readings
    (``q_sign=-1``: filed at ``-q``; ``direct``/``crossed``: drop a part;
    ``crossed_sign``), used only to show that a check can tell them apart.
    Refuses a displacement that is not a label of ``geom``.
    """
    lab = bond_label_index(geom, nd)
    shape = tuple(int(x) for x in geom.spatial_shape)
    nvol = int(np.prod(shape))
    NL = len(lab)
    qv = q_sign * mesh_momenta(shape)
    G = np.zeros((nvol, NL, NL), dtype=np.complex128)

    def idx(d, x, y):
        key = (tuple(int(c) for c in d), int(x), int(y))
        if key not in lab:
            raise ValueError("file_gamma_label: displacement {} is not a label of the geometry"
                             .format(key[0]))
        return lab[key]
    for m in monomials:
        c = np.array(m.cells, dtype=np.int64)
        if c[0].any():
            raise ValueError("file_gamma_label: operator 1 must sit in cell 0, got {}"
                             .format(tuple(c[0])))
        d1, t, d2 = c[1], c[2], c[3] - c[2]
        ph = np.exp(1j * (qv @ t.astype(float)))
        x1, x2, x3, x4 = m.modes
        if direct:
            I, J = idx(d1, x1, x2), idx(d2, x3, x4)
            G[:, I, J] += m.v * ph
            G[:, J, I] += m.v * ph.conj()
        if crossed:
            I, J = idx(t + d2, x1, x4), idx(d1 - t, x3, x2)
            G[:, I, J] += crossed_sign * m.v * ph
            G[:, J, I] += crossed_sign * m.v * ph.conj()
    return G


def first_order_prediction(chibar, gamma_label, spatial_shape, sketches, qs):
    """``u^T [-chibar(q) Gamma(-q) chibar(q)] v`` ``(len(qs), len(sketches))``
    (label form)."""
    neg = mesh_negate(tuple(spatial_shape))
    out = np.empty((len(qs), len(sketches)), dtype=np.complex128)
    for i, q in enumerate(qs):
        P = -chibar[q] @ gamma_label[neg[q]] @ chibar[q]
        for s, (u, v) in enumerate(sketches):
            out[i, s] = u @ P @ v
    return out


def _gamma_slot_maps(gamma, geom, nd, slot_sign, where):
    sr, sc = (int(s) for s in slot_sign)
    if sr not in (1, -1) or sc not in (1, -1):
        raise ValueError("{}: slot_sign must be a pair of +1/-1, got {!r}".format(where, slot_sign))
    D = np.asarray(geom.delta_r, dtype=np.int64)
    rev = np.asarray(geom.reverse, dtype=np.int64)
    shape = tuple(int(x) for x in geom.spatial_shape)
    B, npair = D.shape[0], nd * nd
    NL = B * npair
    if np.shape(gamma) != (int(np.prod(shape)), NL, NL):
        raise ValueError("{}: expected shape {}, got {}".format(
            where, (int(np.prod(shape)), NL, NL), np.shape(gamma)))
    lr = np.empty(NL, dtype=np.int64)
    lc = np.empty(NL, dtype=np.int64)
    for m in range(B):
        mr = m if sr == 1 else rev[m]                 # row label displacement sr dr_m
        mc = m if sc == -1 else rev[m]                # column label displacement -sc dr_m
        for x in range(nd):
            for y in range(nd):
                lr[m * npair + x * nd + y] = mr * npair + x * nd + y
                lc[m * npair + x * nd + y] = mc * npair + y * nd + x
    phase = np.repeat(np.exp(-1j * sc * (mesh_momenta(shape) @ D.T.astype(float))), npair, axis=1)
    return lr, lc, phase


def gamma_label_to_code(gamma_label, geom, nd, slot_sign):
    """The code-slot interaction from the label form:

    ``Gamma_B(q)[(m;x,y),(m';z,w)] = e^{-iq.sc dr_m'}
    Gamma(q)[(sr dr_m; x,y), (-sc dr_m'; w,z)]``, ``(sr, sc) = slot_sign``
    (``H_int = (1/4N) sum Gamma_B O_r O_c^dagger`` for the antisymmetrized
    vertex, with row slot
    ``rho^{(sr dr_m)}_{xy}`` and column slot ``rho^{(sc dr_m')}_{zw}``).
    ``geom`` supplies ``delta_r``, ``reverse`` and the mesh."""
    lr, lc, phase = _gamma_slot_maps(gamma_label, geom, nd, slot_sign, "gamma_label_to_code")
    g = np.asarray(gamma_label)
    return phase[:, None, :] * g[:, lr[:, None], lc[None, :]]


def gamma_code_to_label(gamma_code, geom, nd, slot_sign):
    """The inverse of :func:`gamma_label_to_code`."""
    lr, lc, phase = _gamma_slot_maps(gamma_code, geom, nd, slot_sign, "gamma_code_to_label")
    out = np.empty(np.shape(gamma_code), dtype=np.complex128)
    out[:, lr[:, None], lc[None, :]] = np.asarray(gamma_code) * phase.conj()[:, None, :]
    return out


# ---------------------------------------------------------------------------
# Matsubara Lindhard function, free Green function and the label RPA ring
# ---------------------------------------------------------------------------

def matsubara_pair_kernel(xa, xb, beta, nu):
    """``-T sum_n 1 / ((i w_n - xa)(i w_n - i nu - xb))`` for a bosonic
    ``nu``: ``(f(xb) - f(xa)) / (xa - xb - i nu)``; at ``nu = 0`` it is
    :func:`_pair_kernel` (with its ``-f'`` limit)."""
    if nu == 0:
        return _pair_kernel(xa, xb, beta).astype(np.complex128)
    return (_fermi(xb, beta) - _fermi(xa, beta)) / (xa - xb - 1j * nu)


def lindhard_matsubara(h1_k, mu, beta, geom, ls):
    """The label bubble of the band at the bosonic frequencies
    ``nu_m = 2 pi m T`` for every integer ``m`` of ``ls``, analytically:

    ``chi[i, q, (d;X,Y), (d';Z,W)] = (1/N) sum_k e^{+ik.d} e^{+i(k+q).d'}
    sum_{a,b} U_Ya(k) U*_Za(k) U_Wb(k+q) U*_Xb(k+q) F_nu(xi_a(k), xi_b(k+q))``

    with ``F_nu`` = :func:`matsubara_pair_kernel` at ``nu = nu_{ls[i]}``,
    the Matsubara sum of ``-T sum_n G_YZ(k, i w_n) G_WX(k+q, i w_n - i nu)``
    (the frequency convention of :func:`lindhard_discrete`). ``ls = (0,)``
    is :func:`lindhard_static`.

    Returns
    -------
    ndarray complex128 (len(ls), nvol, NL, NL), ``NL = B * nd**2``.
    """
    h1_k = np.asarray(h1_k, dtype=np.complex128)
    nvol, nd, _ = h1_k.shape
    shape = tuple(int(x) for x in geom.spatial_shape)
    if int(np.prod(shape)) != nvol:
        raise ValueError("lindhard_matsubara: h1_k has nvol={} but the mesh {} has {}"
                         .format(nvol, shape, int(np.prod(shape))))
    D = np.asarray(geom.delta_r, dtype=np.int64)
    NL = D.shape[0] * nd * nd
    e, U = np.linalg.eigh(h1_k)
    xi = e - mu
    pk = np.exp(1j * (mesh_momenta(shape) @ D.T.astype(float)))       # (k, m)
    add = mesh_add(shape)
    ls = tuple(int(m) for m in ls)
    out = np.empty((len(ls), nvol, NL, NL), dtype=np.complex128)
    for i, m in enumerate(ls):
        nu = 2.0 * np.pi * m / beta
        for q in range(nvol):
            kq = add[q]
            F = matsubara_pair_kernel(xi[:, :, None], xi[kq][:, None, :], beta, nu)
            T = np.einsum("kya,kza,kab,kwb,kxb->kyzwx", U, U.conj(), F, U[kq], U[kq].conj(),
                          optimize=True)
            chi = np.einsum("km,kn,kyzwx->mxynzw", pk, pk[kq], T, optimize=True) / nvol
            out[i, q] = chi.reshape(NL, NL)
    return out


def free_green(h1_k, mu, beta, nmat):
    """``G(k, i w_n) = (i w_n + mu - h(k))^{-1}`` ``(nmat, nvol, nd, nd)`` on
    ``w_n = (2n + 1 - nmat) pi T`` (spin-block index)."""
    h1_k = np.asarray(h1_k, dtype=np.complex128)
    nvol, nd, _ = h1_k.shape
    w = (2 * np.arange(nmat) + 1 - nmat) * np.pi / beta
    eye = np.eye(nd)
    A = (1j * w + mu)[:, None, None, None] * eye - h1_k[None]
    return np.linalg.inv(A)


def ring_label(chi_label, gamma_label, spatial_shape):
    """The label RPA ring ``chi(q) = [1 + chibar(q) Gamma(-q)]^{-1}
    chibar(q)``, ``Gamma(-q)`` read at the grid index ``(-q) mod mesh``
    (:func:`mesh_negate` of ``spatial_shape``).

    ``chi_label`` is ``(nvol, NL, NL)`` or ``(nl, nvol, NL, NL)`` (every
    frequency slice dressed by the same ``Gamma``); ``gamma_label`` is
    ``(nvol, NL, NL)``. Returns an array of the shape of ``chi_label``."""
    chi = np.asarray(chi_label, dtype=np.complex128)
    gamma = np.asarray(gamma_label, dtype=np.complex128)
    static = chi.ndim == 3
    if static:
        chi = chi[None]
    nvol, NL = chi.shape[1], chi.shape[2]
    if gamma.shape != (nvol, NL, NL) or int(np.prod(spatial_shape)) != nvol:
        raise ValueError("ring_label: chibar {} and Gamma {} do not share the mesh {}".format(
            chi.shape, gamma.shape, tuple(spatial_shape)))
    neg = mesh_negate(tuple(spatial_shape))
    out = np.linalg.solve(np.eye(NL) + chi @ gamma[neg][None], chi)
    return out[0] if static else out


#: the bonds of the 5x3 grid fixture: x, y and the diagonal x + y
GRID_BONDS = ((1, 0, 0), (0, 1, 0), (1, 1, 0))


def _grid_transfer():
    """norb_phys = 1 on the 5x3 grid: an on-site Hermitian spin matrix and a
    random complex spin matrix (seed 11, times 1/2) on each bond of
    :data:`GRID_BONDS`, with its partner at ``-R``."""
    rng = np.random.default_rng(11)
    h0 = np.array([[0.1, 0.05 + 0.08j], [0.05 - 0.08j, -0.12]])
    tr = {}
    for s in (0, 1):
        for s2 in (0, 1):
            tr[((0, 0, 0), (s, s2))] = h0[s, s2]
    for R in GRID_BONDS:
        t = 0.5 * (rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))
        for s in (0, 1):
            for s2 in (0, 1):
                tr[(R, (s, s2))] = t[s, s2]
                tr[(_neg(R), (s2, s))] = np.conj(t[s, s2])
    return tr


def _grid_tables():
    """CoulombIntra on site; CoulombInter, Exchange and complex PairHop on
    every bond of :data:`GRID_BONDS` (the diagonal at 0.6 times the axis
    value), Hermitian-closed."""
    out = {"CoulombIntra": {((0, 0, 0), (0, 0)): BASE["CoulombIntra"]}}
    for t in ("CoulombInter", "Exchange", "PairHop"):
        tbl = {}
        for R, f in zip(GRID_BONDS, (1.0, 0.8, 0.6)):
            tbl.update(closed_row(R, 0, 0, f * BASE[t]))
        out[t] = tbl
    return out


#: a two-dimensional 5x3 mesh (norb_phys = 1) with complex spin-flip hopping
#: and interactions on the x, y and diagonal bonds (no exact diagonalization)
FIXTURE_GRID_5X3 = FixtureSpec.make((5, 3, 1), 1, 2.0, 0.2, _grid_transfer(), _grid_tables())
