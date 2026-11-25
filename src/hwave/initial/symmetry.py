"""Spin-rotation and time-reversal checks, and the maps that deduplicate
candidate patterns (design section 5, step 6).

Two candidates are merged only under a map that has been VERIFIED on the
model: a check of the solver's own Hartree-Fock map, not a reading of the
interaction tables. No map is assumed.

* :func:`su2_check`: is the solver's mean-field map equivariant under
  global spin rotations, and does ``H0`` commute with the total spin?
* :func:`tr_check`: is the model time-reversal invariant (``H0`` and the
  Hartree-Fock map)?
* :func:`find_map`: the verified map, if any, that takes one pattern to
  another (a lattice translation of the cell first, then time reversal,
  then a global spin rotation).
* :func:`sign_map_verified`: does time reversal or a spin rotation take a
  pattern ``O`` to ``-O``?

**Conventions.** The checks evaluate ``Sigma'[G] = hf_potential(G, ...)``
(``initial_facade``) on the STORED Green function ``G[r, s, a, t, b] =
<c^dag_{sa}(r) c_{tb}(0)>`` with the run's Fock flag and blocks. The
physics is written for ``rho = G^T``, on which a spin rotation ``U`` acts as
``rho -> U rho U^dag`` and the potential as ``Sigma -> U Sigma U^dag``. In
the stored form this reads ``G -> U* G U^T`` (:func:`rotate_green`), and the
generator identity checked is

    Sigma'[i [S_mu^T, e_G]] = -i [S_mu, Sigma'[e_G]],   S_mu = sigma^mu / 2 (x) 1.

Time reversal ``Theta = (i sigma^y) K`` acts on ``G`` as ``G -> Y conj(G)
Y^T`` and on a one-body operator in real space as ``O(R) -> Y conj(O(R))
Y^T`` (``Y = i sigma^y`` on the spin index; ``O(k) -> Y conj(O(-k)) Y^T``).

Patterns are given as the Hamiltonian-convention ``O(k)`` ``(nvol, nd,
nd)`` of the folded cell, in the solver's index order (``realise``
returns them so); the problem passed alongside is the problem of that
cell. A folded orbital is ``aa = a + norb_orig * ir`` with the in-cell index
``ir = ix + Bx*(iy + By*iz)`` (spin-block ``s*norb + aa`` in normal mode,
interleaved inside ``a`` in spin-orbital mode).
"""
import itertools
import logging
from dataclasses import dataclass
from typing import Literal, Optional, Tuple

import numpy as np

from hwave.solver import initial_facade

__all__ = [
    "SU2_TOL", "TR_TOL", "MAP_TOL", "MAX_HF_APPLICATIONS", "SZ_WARNING",
    "Su2Result", "TrResult", "su2_check", "tr_check", "find_map",
    "sign_map_verified", "table_probes", "hermitian_green", "rotate_green",
    "spin_rotation", "time_reverse_green",
]

logger = logging.getLogger(__name__)

#: the relative tolerance of the spin-rotation check (design section 5 step 6)
SU2_TOL = 1e-12

#: the relative tolerance of the time-reversal check
TR_TOL = 1e-12

#: the entrywise relative tolerance of a deduplication map
MAP_TOL = 1e-10

#: above this many Hartree-Fock applications (``4 N_probe``) the generator
#: check is replaced by the random regression (design section 4.4)
MAX_HF_APPLICATIONS = 20000

#: the warning of a failed S_z check in 2Sz mode (design section 5 step 6)
SZ_WARNING = "the interaction tables break spin conservation about z"

_STATUS_SO = "not applicable (spin-orbital)"
_STATUS_2SZ = "not applicable (2Sz fixed)"

_SIGMA = (np.array([[0, 1], [1, 0]], dtype=np.complex128),
          np.array([[0, -1j], [1j, 0]], dtype=np.complex128),
          np.array([[1, 0], [0, -1]], dtype=np.complex128))

#: ``i sigma^y``
_ISY = np.array([[0, 1], [-1, 0]], dtype=np.complex128)


@dataclass(frozen=True)
class Su2Result:
    """The spin-rotation check.

    ``status``: ``"passed"``, ``"failed"``, ``"not applicable (2Sz fixed)"``
    or ``"not applicable (spin-orbital)"``. ``max_residual``: the largest
    relative residual (Hartree-Fock map and ``H0`` commutator), ``None``
    when not applicable. ``method``: ``"generator"`` (the complete
    table-support probes) or ``"random"`` (3 generic G x 3 rotations), or
    ``None``. ``sz_residual``: in 2Sz mode the S_z-only residual, else
    ``None``. ``residuals``: the generator residual of the Hartree-Fock map
    per axis ``(x, y, z)`` (``None`` for ``"random"`` or when not checked;
    in 2Sz mode only z is set). ``h0_residual``: the relative
    ``max_{mu,R} ||[S_mu, H0(R)]||``. ``warnings``: messages for the
    report.

    Every residual is relative: divided by ``max(1, ||Sigma'||)`` (the
    largest Frobenius norm of the potential over the probes) or by
    ``max(1, ||H0||)``; the check passes when each is ``<= SU2_TOL``.
    """
    status: Literal["passed", "failed", "not applicable (2Sz fixed)",
                    "not applicable (spin-orbital)"]
    max_residual: Optional[float]
    method: Optional[Literal["generator", "random"]]
    sz_residual: Optional[float]
    residuals: Optional[Tuple[Optional[float], Optional[float], Optional[float]]] = None
    h0_residual: Optional[float] = None
    warnings: Tuple[str, ...] = ()


@dataclass(frozen=True)
class TrResult:
    """The time-reversal check: ``h0_residual`` of ``Theta H0(k) Theta^-1 =
    H0(-k)`` and ``hf_residual`` of ``Sigma'[Theta G Theta^-1] = Theta
    Sigma'[G] Theta^-1`` on 3 generic G, both relative; ``invariant`` when
    both are ``<= TR_TOL`` and ``2Sz`` is free or 0."""
    invariant: bool
    h0_residual: float
    hf_residual: float


# --------------------------------------------------------------------------
# helpers on the stored layout
# --------------------------------------------------------------------------

def _neg_r_index(shape):
    """The flat index of ``-r`` for each flat ``r`` (C order over ``shape``)."""
    nx, ny, nz = shape
    ix, iy, iz = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    return (((-iz) % nz) + nz * (((-iy) % ny) + ny * ((-ix) % nx))).ravel()


def hermitian_green(green, shape):
    """The Hermitian part of a stored Green function: ``(G(r) + G(-r)^dag)
    / 2`` in the ``(nd, nd)`` view, so that ``G(r)_ab = conj(G(-r)_ba)``."""
    g = np.asarray(green, dtype=np.complex128)
    nvol = g.shape[0]
    nd = g.shape[1] * g.shape[2]
    m = g.reshape(nvol, nd, nd)
    partner = np.conj(np.swapaxes(m[_neg_r_index(shape)], -1, -2))
    return ((m + partner) / 2).reshape(g.shape)


def spin_rotation(axis, theta, norb):
    """``U = exp(-i theta n.S) (x) 1_norb`` in the spin-block order."""
    n = np.asarray(axis, dtype=float)
    n = n / np.linalg.norm(n)
    ns = sum(c * s for c, s in zip(n, _SIGMA))
    u = np.cos(theta / 2) * np.eye(2) - 1j * np.sin(theta / 2) * ns
    return np.kron(u, np.eye(norb))


def rotate_green(green, u):
    """The stored Green function of the rotated state: ``U* G U^T`` per
    ``r`` (``rho -> U rho U^dag`` for ``rho = G^T``)."""
    g = np.asarray(green)
    nvol = g.shape[0]
    nd = u.shape[0]
    m = g.reshape(nvol, nd, nd)
    return (np.conj(u) @ m @ u.T).reshape(g.shape)


def _theta_matrix(problem):
    """``Y = i sigma^y`` on the spin index, in the solver's index order."""
    if problem.enable_spin_orbital:
        return np.kron(np.eye(problem.nd // 2), _ISY)
    return np.kron(_ISY, np.eye(problem.norb))


def time_reverse_green(green, problem):
    """The stored Green function of the time-reversed state, ``Y conj(G)
    Y^T`` per ``r``."""
    g = np.asarray(green)
    y = _theta_matrix(problem)
    nd = y.shape[0]
    m = g.reshape(g.shape[0], nd, nd)
    return (y @ np.conj(m) @ y.T).reshape(g.shape)


def _sigma_map(problem):
    tables = initial_facade.hf_tables(problem)

    def sigma(green):
        return initial_facade.hf_potential(green, tables, include_fock=problem.iflag_fock,
                                           blocks=problem.blocks)
    return sigma


def _green_shape(problem):
    return (problem.nvol, problem.ns, problem.norb, problem.ns, problem.norb)


# --------------------------------------------------------------------------
# probes
# --------------------------------------------------------------------------

def _table_slots(problem):
    """The ``(r, i, j)`` orbital slots of the stored G (``G[r, ., i, ., j]``)
    the two-body tables read, closed under the Hermitian partner
    ``(-r, j, i)``; sorted."""
    tables = problem.tables
    shape = problem.shape
    nvol = problem.nvol
    neg = _neg_r_index(shape)
    slots = set()
    for name in sorted(tables.inter_table):
        tab = tables.inter_table[name]
        if tab is None:
            continue
        n = tab.shape[-1]
        jab = np.asarray(tab).reshape(nvol, n, n)
        for r, a, b in zip(*np.nonzero(jab)):
            r, a, b = int(r), int(a), int(b)
            # the Fock and PairHop terms read G[r, ., b, ., a]; the Hartree
            # term the on-site densities of b (and, by the reverse-closed
            # table, of a)
            for s in ((r, b, a), (0, a, a), (0, b, b)):
                slots.add(s)
                slots.add((int(neg[s[0]]), s[2], s[1]))
    return sorted(slots)


def table_probes(problem):
    """The complete Hermitian probe basis of the table support (design
    section 5 step 6): for each slot ``(r, i, j)`` read by a two-body table
    and each spin pair ``(s, t)``, the stored Green functions ``E + E^dag``
    and ``i (E - E^dag)``, ``E`` the unit at ``G[r, s, i, t, j]`` and
    ``E^dag`` its partner at ``G[-r, t, j, s, i]`` (zero probes dropped,
    each pair once). Normal mode only. Returns a list of arrays ``(nvol, 2,
    norb, 2, norb)``."""
    shape = problem.shape
    neg = _neg_r_index(shape)
    gshape = _green_shape(problem)
    probes = []
    seen = set()
    for (r, i, j) in _table_slots(problem):
        for s, t in itertools.product(range(problem.ns), repeat=2):
            idx = (r, s, i, t, j)
            partner = (int(neg[r]), t, j, s, i)
            key = min(idx, partner)
            if key in seen:
                continue
            seen.add(key)
            e = np.zeros(gshape, dtype=np.complex128)
            e[idx] += 1.0
            e[partner] += 1.0
            probes.append(e)
            if idx != partner:
                e = np.zeros(gshape, dtype=np.complex128)
                e[idx] += 1j
                e[partner] -= 1j
                probes.append(e)
    return probes


# --------------------------------------------------------------------------
# SU(2)
# --------------------------------------------------------------------------

def _h0_commutator(problem, axes):
    h0 = np.asarray(problem.ham_trans_r)
    nd = problem.nd
    mats = h0.reshape(-1, nd, nd)
    scale = max(1.0, float(max(np.linalg.norm(m) for m in mats)))
    res = 0.0
    for mu in axes:
        S = np.kron(_SIGMA[mu] / 2, np.eye(problem.norb))
        res = max(res, float(np.max(np.linalg.norm(S @ mats - mats @ S, axis=(1, 2)))))
    return res / scale


def _generator_residuals(problem, probes, axes):
    sigma = _sigma_map(problem)
    nd = problem.nd
    gshape = _green_shape(problem)
    ops = {mu: np.kron(_SIGMA[mu] / 2, np.eye(problem.norb)) for mu in axes}
    res = {mu: 0.0 for mu in axes}
    norm = 0.0
    for e in probes:
        s0 = sigma(e)
        norm = max(norm, float(np.linalg.norm(s0)))
        e2 = e.reshape(problem.nvol, nd, nd)
        for mu in axes:
            S = ops[mu]
            St = S.T
            lhs = sigma((1j * (St @ e2 - e2 @ St)).reshape(gshape))
            rhs = -1j * (S @ s0 - s0 @ S)
            res[mu] = max(res[mu], float(np.linalg.norm(lhs - rhs)))
    scale = max(1.0, norm)
    return {mu: r / scale for mu, r in res.items()}


def _random_residual(problem, axes, rng_seed):
    """3 generic G x 3 rotations (about z only when ``axes == (2,)``)."""
    sigma = _sigma_map(problem)
    rng = np.random.default_rng(rng_seed)
    gshape = _green_shape(problem)
    res = 0.0
    norm = 0.0
    for _ in range(3):
        g = hermitian_green(rng.normal(size=gshape) + 1j * rng.normal(size=gshape),
                            problem.shape)
        s0 = sigma(g)
        norm = max(norm, float(np.linalg.norm(s0)))
        for _ in range(3):
            n = (0.0, 0.0, 1.0) if axes == (2,) else rng.normal(size=3)
            theta = rng.uniform(0.1, np.pi)
            u = spin_rotation(n, theta, problem.norb)
            lhs = sigma(rotate_green(g, u))
            rhs = u @ s0 @ u.conj().T
            res = max(res, float(np.linalg.norm(lhs - rhs)))
    return res / max(1.0, norm)


def su2_check(problem, *, rng_seed=364, max_hf_applications=MAX_HF_APPLICATIONS):
    """The spin-rotation equivariance check of ``problem`` (design section 5
    step 6).

    Normal mode with free S_z: the generator identity ``Sigma'[i[S_mu^T,
    e]] = -i[S_mu, Sigma'[e]]`` on every probe of :func:`table_probes`,
    mu = x, y, z, and ``[S_mu, H0(R)] = 0`` for every R. When ``4 N_probe >
    max_hf_applications`` the probes are replaced by 3 generic G x 3
    random rotations (``default_rng(rng_seed)``), ``method = "random"``.

    2Sz mode: only S_z is checked; the status is ``"not applicable (2Sz
    fixed)"`` and the residual goes to ``sz_residual``; a failure adds the
    warning :data:`SZ_WARNING` (logged and in ``warnings``) and changes
    nothing else. Spin-orbital mode: ``"not applicable (spin-orbital)"``.
    """
    if problem.enable_spin_orbital:
        return Su2Result(status=_STATUS_SO, max_residual=None, method=None,
                         sz_residual=None)
    axes = (0, 1, 2) if problem.sz_free else (2,)
    h0 = _h0_commutator(problem, axes)
    probes = table_probes(problem)
    if 4 * len(probes) > max_hf_applications:
        method = "random"
        hf = _random_residual(problem, axes, rng_seed)
        per_axis = None
    else:
        method = "generator"
        r = _generator_residuals(problem, probes, axes)
        hf = max(r.values()) if r else 0.0
        per_axis = tuple(r.get(mu) for mu in range(3))
    worst = max(hf, h0)
    if not problem.sz_free:
        warnings = ()
        if worst > SU2_TOL:
            logger.warning(SZ_WARNING)
            warnings = (SZ_WARNING,)
        return Su2Result(status=_STATUS_2SZ, max_residual=None, method=method,
                         sz_residual=worst, residuals=per_axis, h0_residual=h0,
                         warnings=warnings)
    status = "passed" if (hf <= SU2_TOL and h0 <= SU2_TOL) else "failed"
    return Su2Result(status=status, max_residual=worst, method=method, sz_residual=None,
                     residuals=per_axis, h0_residual=h0)


# --------------------------------------------------------------------------
# time reversal
# --------------------------------------------------------------------------

def _reverse_k(a, shape):
    """``a(-k)`` for ``a`` ``(nvol, ...)`` on the C-ordered k-grid."""
    return np.asarray(a)[_neg_r_index(shape)]


def tr_check(problem, *, rng_seed=364):
    """The time-reversal check of ``problem``: ``Theta H0(k) Theta^-1 =
    H0(-k)`` (in real space ``Y conj(H0(R)) Y^T = H0(R)`` for every R) and
    ``Sigma'[Theta G Theta^-1] = Theta Sigma'[G] Theta^-1`` on 3 generic G
    (``default_rng(rng_seed)``), each to ``TR_TOL`` relative. A fixed
    ``2Sz != 0`` is never invariant (time reversal maps the sector to
    ``-2Sz``)."""
    y = _theta_matrix(problem)
    nd = problem.nd
    h0 = np.asarray(problem.ham_trans_r).reshape(-1, nd, nd)
    scale = max(1.0, float(max(np.linalg.norm(m) for m in h0)))
    h0_res = float(np.max(np.linalg.norm(y @ np.conj(h0) @ y.T - h0, axis=(1, 2)))) / scale

    sigma = _sigma_map(problem)
    rng = np.random.default_rng(rng_seed)
    gshape = _green_shape(problem)
    res = 0.0
    norm = 0.0
    for _ in range(3):
        g = hermitian_green(rng.normal(size=gshape) + 1j * rng.normal(size=gshape),
                            problem.shape)
        s0 = sigma(g)
        norm = max(norm, float(np.linalg.norm(s0)))
        lhs = sigma(time_reverse_green(g, problem))
        rhs = y @ np.conj(_reverse_k(s0, problem.shape)) @ y.T
        res = max(res, float(np.linalg.norm(lhs - rhs)))
    hf_res = res / max(1.0, norm)
    sz_ok = problem.two_sz is None or int(problem.two_sz) == 0
    invariant = bool(sz_ok and h0_res <= TR_TOL and hf_res <= TR_TOL)
    return TrResult(invariant=invariant, h0_residual=h0_res, hf_residual=hf_res)


# --------------------------------------------------------------------------
# maps between patterns
# --------------------------------------------------------------------------

def _to_r(o_k, problem):
    """``O(R)`` ``(nx, ny, nz, nd, nd)`` of a Hamiltonian-convention
    ``O(k)`` (the inverse of ``O(k) = sum_R e^{ikR} O(R)``)."""
    nx, ny, nz = problem.shape
    nd = problem.nd
    o = np.asarray(o_k, dtype=np.complex128)
    if o.shape != (problem.nvol, nd, nd):
        raise ValueError("a pattern of this problem has shape {}, got {}".format(
            (problem.nvol, nd, nd), o.shape))
    return np.fft.fftn(o.reshape(nx, ny, nz, nd, nd), axes=(0, 1, 2), norm="forward")


def _subshape(problem):
    return tuple(int(c) // int(n) for c, n in zip(problem.cellshape, problem.shape))


def _translations(problem):
    """The in-cell translations ``tau`` of the folded cell (identity first)."""
    bx, by, bz = _subshape(problem)
    return list(itertools.product(range(bx), range(by), range(bz)))


def _translate(o_r, problem, tau):
    """``O`` translated by the primitive-cell vector ``tau``.

    In the Hamiltonian convention ``O(R)[aa, bb]`` couples folded orbital
    ``aa`` of supercell 0 (in-cell offset ``r_i``) with ``bb`` of supercell
    ``R`` (offset ``r_j``). Moving both by ``tau`` sends ``r -> (r + tau)
    mod B`` with the carries ``c``, so ``O'(R + c_j - c_i)[aa', bb'] =
    O(R)[aa, bb]``."""
    B = np.array(_subshape(problem))
    shape = np.array(problem.shape)
    nd = problem.nd
    if problem.enable_spin_orbital:
        norb_cell, nspin = problem.norb_orig, 1
    else:
        norb_cell, nspin = problem.norb_orig, 2
    norb = problem.norb
    idx = np.arange(nd)
    if nspin == 2:
        s, aa = idx // norb, idx % norb
    else:
        s, aa = np.zeros(nd, dtype=int), idx
    a, ir = aa % norb_cell, aa // norb_cell
    off = np.stack([ir % B[0], (ir // B[0]) % B[1], ir // (B[0] * B[1])], axis=1)
    moved = off + np.array(tau)
    carry = moved // B
    moved = moved % B
    ir2 = moved[:, 0] + B[0] * (moved[:, 1] + B[1] * moved[:, 2])
    aa2 = a + norb_cell * ir2
    new = s * norb + aa2 if nspin == 2 else aa2
    # d[i, j] = c_j - c_i; the map (R, i, j) -> (R + d[i, j], new[i], new[j])
    # is a bijection, so one scatter writes every entry once
    d = carry[None, :, :] - carry[:, None, :]
    dst = [(np.arange(n).reshape([-1 if ax == k else 1 for k in range(3)] + [1, 1])
            + d[None, None, None, :, :, ax]) % n for ax, n in enumerate(shape)]
    out = np.zeros_like(o_r)
    out[dst[0], dst[1], dst[2], new[:, None], new[None, :]] = o_r
    return out


def _time_reverse_r(o_r, problem):
    y = _theta_matrix(problem)
    return y @ np.conj(o_r) @ y.T


def _spin_components(o_r, norb):
    """``c^mu = tr_spin(sigma^mu O)`` (mu = 0, x, y, z), each ``(..., norb,
    norb)``, normal mode."""
    uu = o_r[..., :norb, :norb]
    ud = o_r[..., :norb, norb:]
    du = o_r[..., norb:, :norb]
    dd = o_r[..., norb:, norb:]
    return uu + dd, ud + du, 1j * (ud - du), uu - dd


def _from_components(c0, cv, norb):
    cx, cy, cz = cv
    out = np.empty(c0.shape[:-2] + (2 * norb, 2 * norb), dtype=np.complex128)
    out[..., :norb, :norb] = (c0 + cz) / 2
    out[..., norb:, norb:] = (c0 - cz) / 2
    out[..., :norb, norb:] = (cx - 1j * cy) / 2
    out[..., norb:, :norb] = (cx + 1j * cy) / 2
    return out


def _rotated_fit(o1, o2, norb):
    """``U O1 U^dag`` for the proper rotation that best maps the spin part of
    ``o1`` onto that of ``o2`` (Kabsch); the caller checks the fit."""
    c1 = _spin_components(o1, norb)
    c2 = _spin_components(o2, norb)
    m1 = np.array([np.concatenate([c.ravel().real, c.ravel().imag]) for c in c1[1:]])
    m2 = np.array([np.concatenate([c.ravel().real, c.ravel().imag]) for c in c2[1:]])
    u, _, vt = np.linalg.svd(m1 @ m2.T)
    d = np.sign(np.linalg.det(vt.T @ u.T)) or 1.0
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    cv = [sum(rot[i, j] * c1[1 + j] for j in range(3)) for i in range(3)]
    return _from_components(c1[0], cv, norb)


def _close(a, b, tol):
    return float(np.max(np.abs(a - b))) <= tol


def _map_tol(o1_r, o2_r):
    return MAP_TOL * max(float(np.linalg.norm(o1_r)), float(np.linalg.norm(o2_r)))


def _su2_ok(problem, su2):
    return (not problem.enable_spin_orbital) and su2 is not None and su2.status == "passed"


def _tr_ok(problem, tr):
    return tr is not None and bool(tr.invariant) and (
        problem.two_sz is None or int(problem.two_sz) == 0)


def find_map(O1, O2, *, problem, su2,
             tr) -> Optional[Literal["translation", "time_reversal", "spin_rotation"]]:
    """The verified map taking pattern ``O1`` to ``O2``, or ``None``.

    Tried in this order, each entrywise to ``MAP_TOL * max(||O1||, ||O2||)``
    on ``O(R)``: (a) ``"translation"``, a lattice translation of the cell
    (the identity included); (b) ``"time_reversal"``, time reversal followed
    by a translation, only when ``tr.invariant`` (and ``2Sz`` free or 0),
    also combined with a global spin rotation when (c) is available; (c)
    ``"spin_rotation"``, a global spin rotation followed by a translation,
    only when ``su2.status == "passed"`` (never assumed).

    ``O1``, ``O2``: Hamiltonian-convention ``O(k)`` ``(nvol, nd, nd)`` of
    ``problem``'s folded cell. Raises ``ValueError`` for another shape.
    """
    o1 = _to_r(O1, problem)
    o2 = _to_r(O2, problem)
    tol = _map_tol(o1, o2)
    taus = _translations(problem)
    for tau in taus:
        if _close(_translate(o1, problem, tau), o2, tol):
            return "translation"
    tr_on = _tr_ok(problem, tr)
    su2_on = _su2_ok(problem, su2)
    if tr_on:
        o1t = _time_reverse_r(o1, problem)
        for tau in taus:
            if _close(_translate(o1t, problem, tau), o2, tol):
                return "time_reversal"
    if su2_on:
        norb = problem.norb
        for tau in taus:
            moved = _translate(o1, problem, tau)
            if _close(_rotated_fit(moved, o2, norb), o2, tol):
                return "spin_rotation"
        if tr_on:
            for tau in taus:
                moved = _translate(_time_reverse_r(o1, problem), problem, tau)
                if _close(_rotated_fit(moved, o2, norb), o2, tol):
                    return "time_reversal"
    return None


def sign_map_verified(O, *, problem, su2, tr) -> bool:
    """Whether a verified map (b) time reversal or (c) a global spin
    rotation takes ``O`` to ``-O`` (no translation: at q not equivalent to
    -q a translation is a phase shift, design section 5 step 6)."""
    o = _to_r(O, problem)
    tol = _map_tol(o, o)
    if _tr_ok(problem, tr) and _close(_time_reverse_r(o, problem), -o, tol):
        return True
    if _su2_ok(problem, su2) and _close(_rotated_fit(o, -o, problem.norb), -o, tol):
        return True
    return False
