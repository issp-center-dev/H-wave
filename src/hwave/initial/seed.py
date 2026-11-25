"""Field seeding (design section 7, owner decision O1).

A candidate with the pattern ``O`` (the Hamiltonian-convention ``O(k)`` of
the candidate cell, :func:`hwave.initial.patterns.realise`) is seeded by
the density of ``H_base - h O``:

* ``G(h) = green_function(H_base - h O, T_run, blocks, fill="uniform")``,
  the stored Green function, with the solver's own blocks and mu-groups,
  so ``Ncond`` and ``2Sz`` are exact; ``rho(h) = G(h)^T``. At T = 0 the
  levels within ``shell_tol = 1e-10 * max(1, W0, h ||O||_op)`` of the Fermi
  level of a group form its shell and are filled uniformly (O7).
* ``m(h) = Re tr[O rho(h)] = Re sum_ab O_ab G_ab(h)`` per primitive cell
  (summed over k and bond slots); ``Re tr[O G]`` would be ``<O^T>``,
  which reverses every imaginary part.
* ``m(h)`` is non-decreasing (the free energy is concave in h). Its exact
  limit ``m_inf`` is the Ky Fan sum: per mu-group the ``n_group`` largest
  eigenvalues of ``O(k)`` over the group's blocks and every k. The
  occupations inside an O eigenspace that the cut splits stay selected by
  ``H_base``; only their total has a limit.
* ``h`` is found on the change ``Delta m(h) = m(h) - m(0)``: target
  ``Delta m* = amplitude * (m_inf - m(0))``; doubling from ``W_floor =
  max(W(H_base), 1e-3)`` (W the all-band span) up to ``2^40 W_floor``, then
  bisection to ``|Delta m - Delta m*| <= 1e-3 (m_inf - m(0))`` or 60 steps
  (then the smallest h found with ``Delta m >= Delta m*``).

Before any seeding O is checked against the solver blocks (a component
outside them is refused, naming it), then against the constraint modes
(the electron number of a mu-group cannot respond).

``T_hint`` never enters: ``W0`` only scales the shell tolerance.
"""
import dataclasses
import math
import numbers
from dataclasses import dataclass
from typing import Literal, Mapping, Optional, Tuple

import numpy as np

from hwave.solver import initial_facade

__all__ = [
    "SEED_BASES", "BLOCK_RTOL", "NO_RESPONSE_TOL", "W_FLOOR_MIN", "SHELL_RTOL",
    "BISECTION_RTOL", "BISECTION_STEPS", "CAP_DOUBLINGS", "CAP_REASON",
    "DEGENERATE_REASON", "CONSTRAINT_REASON", "BLOCK_REASON", "SeedResult",
    "block_check", "m_value", "m_inf", "constraint_mode", "w_floor", "op_norm",
    "density_at", "moments_of", "seed_pattern",
]

#: the recorded base Hamiltonian of a seed
SEED_BASES = ("para_hf", "bare")
#: a pattern component outside the solver blocks with more than this
#: fraction of ``||O||_F`` is refused
BLOCK_RTOL = 1e-10
#: ``m_inf - m(0)`` below this: "no response (degenerate with para)"
NO_RESPONSE_TOL = 1e-6
#: the floor of ``W_floor`` (input energy units)
W_FLOOR_MIN = 1e-3
#: the T = 0 shell tolerance, relative to ``max(1, W0, h ||O||_op)``
SHELL_RTOL = 1e-10
#: the bisection stops at ``|Delta m - Delta m*| <=`` this times
#: ``m_inf - m(0)``
BISECTION_RTOL = 1e-3
#: the bisection stops after this many steps
BISECTION_STEPS = 60
#: the numerical guard: at most ``2^CAP_DOUBLINGS * W_floor``
CAP_DOUBLINGS = 40
_CAP_DOUBLINGS = CAP_DOUBLINGS
#: a constraint mode is ``O`` within this fraction of its norm of a
#: combination of the mu-group number operators
CONSTRAINT_RTOL = 1e-10

CAP_REASON = "target not reached at the cap"
DEGENERATE_REASON = "no response (degenerate with para)"
CONSTRAINT_REASON = "no response (constraint)"
BLOCK_REASON = "refused: block check"


@dataclass(frozen=True, eq=False)
class SeedResult:
    """The seed of one row (design section 7; the null rules of its table).

    ``tried``, ``not_tried_reason`` (``None`` when tried). ``G``: the
    stored Green function ``(nvol, ns, norb, ns, norb)`` written as the
    initial state (tried rows only). ``seed_base`` (``"para_hf"`` or
    ``"bare"``). ``shell``: ``{"h0": [...], "h_used": [...]}``, one record
    per mu-group, ``{"group", "open_shell", "shell_degeneracy",
    "shell_fill"}`` (``"h_used"`` is ``None`` when no seed was made).
    ``moments_h0``: ``{"spin": [[x, y, z] per site], "charge": [...]}`` of
    ``rho(0)``, per site of the folded cell (site group ``g`` of in-cell
    offset ``ir`` is site ``g + n_groups * ir``). ``m0``, ``m_inf``,
    ``m_target`` (``m0 + Delta m*``), ``m_used``, ``h_hi`` (the upper end of
    the doubling bracket), ``h_used``, ``amplitude_requested``,
    ``amplitude_used`` (``Delta m(h_used) / (m_inf - m0)``) and
    ``cap_reason`` (only at the cap).

    Set fields per state: seeded: all (``cap_reason`` only at the cap);
    para: ``G``, ``seed_base``, ``shell``, ``moments_h0``, ``h_used = 0``;
    no response (degenerate with para): ``seed_base``, ``shell``,
    ``moments_h0``, ``m0``, ``m_inf``; no response (constraint):
    ``seed_base``, ``m0``; refused (block check) and seeding not performed
    (:meth:`not_performed`): none.
    """
    tried: bool
    not_tried_reason: Optional[str] = None
    G: Optional[np.ndarray] = None
    seed_base: Optional[str] = None
    shell: Optional[Mapping] = None
    moments_h0: Optional[Mapping] = None
    m0: Optional[float] = None
    m_inf: Optional[float] = None
    m_target: Optional[float] = None
    m_used: Optional[float] = None
    h_hi: Optional[float] = None
    h_used: Optional[float] = None
    amplitude_requested: Optional[float] = None
    amplitude_used: Optional[float] = None
    cap_reason: Optional[str] = None

    @classmethod
    def not_performed(cls, reason: str) -> "SeedResult":
        """A row whose seeding was not performed (quota, unrepresentable,
        duplicate, a null ``Gamma x``): every seed field null."""
        return cls(tried=False, not_tried_reason=str(reason))


# --------------------------------------------------------------------------
# checks on a pattern
# --------------------------------------------------------------------------

def _block_of(blocks, nd):
    out = np.full(nd, -1)
    for k, blk in enumerate(blocks.block_info):
        out[np.asarray(blk, dtype=int)] = k
    return out


def _as_pattern(O_k, nvol=None, nd=None):
    o = np.asarray(O_k, dtype=np.complex128)
    if o.ndim != 3 or o.shape[1] != o.shape[2]:
        raise ValueError(f"a pattern O(k) must be (nvol, nd, nd), got shape {o.shape}")
    if nvol is not None and o.shape != (nvol, nd, nd):
        raise ValueError(f"the pattern O(k) of this cell has shape {(nvol, nd, nd)}, "
                         f"got {o.shape}")
    return o


def block_check(O_k, blocks, *, ns=None, norb=None) -> list:
    """The components of ``O_k`` outside the solver ``blocks``
    (:class:`SolverBlocks`) that hold more than :data:`BLOCK_RTOL` of
    ``||O||_F``, one message each; ``[]`` when O passes. With ``ns = 2``
    and ``norb`` (normal mode, spin-block index ``s*norb + a``) the
    out-of-block part is split into transverse spin (spin-off-diagonal)
    and spin-diagonal entries, so the message names the channel (e.g.
    transverse spin in 2Sz mode)."""
    o = _as_pattern(O_k)
    nd = o.shape[-1]
    total = float(np.linalg.norm(o))
    if total == 0.0:
        return []
    block_of = _block_of(blocks, nd)
    cross = block_of[:, None] != block_of[None, :]
    if ns == 2 and norb is not None:
        spin = np.arange(nd) // int(norb)
        off = spin[:, None] != spin[None, :]
        parts = (("transverse spin (spin-off-diagonal) entries", cross & off),
                 ("spin-diagonal entries", cross & ~off))
    else:
        parts = (("entries", cross),)
    msgs = []
    for name, mask in parts:
        norm = float(np.linalg.norm(o[:, mask]))
        if norm > BLOCK_RTOL * total:
            big = np.max(np.abs(o), axis=0) * mask > BLOCK_RTOL * total
            pairs = [(int(i), int(j)) for i, j in zip(*np.nonzero(big))][:4]
            msgs.append("{} outside the solver blocks ({:.3g} of the pattern norm; "
                        "index pairs {})".format(name, norm / total, pairs))
    return msgs


def constraint_mode(O_k, facts) -> bool:
    """Whether ``O_k`` is a constraint mode of the cell of ``facts``: a
    combination of the mu-group number operators (the uniform charge of a
    group; with ``2Sz`` fixed, one group per spin, so the uniform S^z too).
    Such a pattern cannot move the seed (each group's electron number is
    fixed)."""
    p = facts.problem
    o = _as_pattern(O_k, p.nvol, p.nd)
    total = float(np.linalg.norm(o))
    if total == 0.0:
        return False
    group_of = np.full(p.nd, -1)
    for blk, g in zip(p.blocks.block_info, p.blocks.block_to_group):
        group_of[np.asarray(blk, dtype=int)] = g
    fit = np.zeros(p.nd, dtype=np.complex128)
    diag = np.einsum("kii->i", o) / p.nvol
    for g in set(int(x) for x in group_of):
        sel = group_of == g
        fit[sel] = np.mean(diag[sel])
    residual = o - np.diag(fit)[None, :, :]
    return float(np.linalg.norm(residual)) <= CONSTRAINT_RTOL * total


# --------------------------------------------------------------------------
# m, m_inf and the scales
# --------------------------------------------------------------------------

def _bvol(problem):
    return int(np.prod(problem.subshape))


def m_value(O_k, G, *, problem) -> float:
    """``m = Re sum_ab O_ab G_ab`` summed over k, per primitive cell
    (``Re tr[O rho]``, rho = G^T): the expectation of ``sum O_ab c^dag_a
    c_b`` per primitive cell. ``G`` is the stored Green function ``(nvol,
    ns, norb, ns, norb)`` of ``problem``'s cell, brought to k by the
    solver's ``to_k`` (``G(k) = sum_r e^{ikr} G(r)``, the inverse of the
    fill)."""
    p = problem
    o = _as_pattern(O_k, p.nvol, p.nd)
    g = np.asarray(G)
    if g.shape != (p.nvol, p.ns, p.norb, p.ns, p.norb):
        raise ValueError(f"G of this cell has shape {(p.nvol, p.ns, p.norb, p.ns, p.norb)}, "
                         f"got {g.shape}")
    gk = initial_facade.to_k(g.reshape(tuple(p.shape) + (p.nd, p.nd)))
    return float(np.real(np.einsum("kab,kab->", o, gk))) / (p.nvol * _bvol(p))


def m_inf(O_k, blocks, *, bvol) -> float:
    """The exact limit of ``m(h)`` as ``h -> infinity`` at fixed electron
    numbers: per mu-group the sum of the ``n_group`` largest eigenvalues
    of ``O(k)`` over the group's blocks and every k (the Ky Fan maximum,
    unchanged by a degeneracy at the cut), summed over the groups, per
    primitive cell (``nvol * bvol`` cells). O must lie within the blocks
    (:func:`block_check`)."""
    o = _as_pattern(O_k)
    nvol = o.shape[0]
    by_group = {}
    for blk, g in zip(blocks.block_info, blocks.block_to_group):
        idx = np.asarray(blk, dtype=int)
        ev = np.linalg.eigvalsh(o[:, idx[:, None], idx[None, :]]).ravel()
        by_group.setdefault(int(g), []).append(ev)
    total = 0.0
    for g, n in enumerate(blocks.group_nconds):
        n = int(round(float(n)))
        if n <= 0 or g not in by_group:
            continue
        ev = np.sort(np.concatenate(by_group[g]))
        total += float(np.sum(ev[len(ev) - n:]))
    return total / (nvol * int(bvol))


def w_floor(base_ham_k) -> float:
    """``W_floor = max(W(H_base), 1e-3)``, W the all-band span (max - min
    eigenvalue over every band and k)."""
    ev = np.linalg.eigvalsh(np.asarray(base_ham_k))
    return max(float(ev.max() - ev.min()), W_FLOOR_MIN)


def op_norm(O_k) -> float:
    """``||O||_op``: the largest ``|eigenvalue|`` of ``O(k)`` over k."""
    o = _as_pattern(O_k)
    return float(np.max(np.abs(np.linalg.eigvalsh(o)), initial=0.0))


# --------------------------------------------------------------------------
# the density
# --------------------------------------------------------------------------

def density_at(facts, base_ham_k, O_k, h, *, T_run, W0, o_norm=None):
    """``(G(h), FillInfo)``: the stored Green function of ``H_base - h O``
    on the cell of ``facts.problem`` with the solver's blocks, mu rule
    and Fermi cut, uniform shell fill at T = 0 with ``shell_tol = 1e-10 *
    max(1, W0, h ||O||_op)`` (``o_norm`` may pass ``||O||_op``)."""
    p = facts.problem
    base = np.asarray(base_ham_k)
    h = float(h)
    if O_k is None or h == 0.0:
        ham, scale = base, 0.0
    else:
        o = _as_pattern(O_k, p.nvol, p.nd)
        ham = base - h * o
        scale = h * (op_norm(o) if o_norm is None else float(o_norm))
    tol = SHELL_RTOL * max(1.0, float(W0), scale)
    return initial_facade.green_function(
        ham, T_run, p.blocks, shape=p.shape, ns=p.ns, norb=p.norb, fill="uniform",
        shell_tol=tol, ene_cutoff=p.ene_cutoff, return_info=True)


def _sites(facts):
    """The site of each folded physical orbital ``j = a + norb_phys * ir``."""
    p = facts.problem
    nphys = int(facts.norb_phys)
    site_of = np.empty(nphys, dtype=int)
    groups = facts.site_groups.groups
    for g, members in enumerate(groups):
        for a in members:
            site_of[a] = g
    j = np.arange(nphys * _bvol(p))
    return site_of[j % nphys] + len(groups) * (j // nphys)


def moments_of(G, facts) -> dict:
    """``{"spin": [[x, y, z], ...], "charge": [...]}`` per site of the folded
    cell (:class:`SeedResult`), from PR 1a's
    :func:`hwave.results.order_parameters_from_density` fed the stored
    Green slot ``rho[o, s, t] = G[0, s, o, t, o]`` (normal mode), or
    :func:`hwave.results.so_density_from_green` in spin-orbital mode."""
    from hwave.results import order_parameters_from_density, so_density_from_green
    p = facts.problem
    g = np.asarray(G)
    if p.enable_spin_orbital:
        rho = so_density_from_green(g[0, 0, :, 0, :], int(facts.norb_phys), _bvol(p))
    else:
        rho = np.einsum("sata->ast", g[0])
    op = order_parameters_from_density(rho, site_of_orbital=_sites(facts), n_cell=None)
    return {"spin": np.asarray(op.spin).tolist(), "charge": np.asarray(op.charge).tolist()}


def _shell(info):
    return [dataclasses.asdict(r) for r in info.shell]


def _read_only(a):
    a = np.array(a, copy=True)
    a.setflags(write=False)
    return a


# --------------------------------------------------------------------------
# the seed
# --------------------------------------------------------------------------

def _check_amplitude(amplitude):
    if isinstance(amplitude, bool) or not isinstance(amplitude, numbers.Real) \
            or not math.isfinite(amplitude) or not 0.0 < float(amplitude) < 1.0:
        raise ValueError(
            "amplitude must be in (0, 1), got {!r}: the target is amplitude * "
            "(m_inf - m(0)), and amplitude = 1 needs h = infinity at T > 0".format(amplitude))
    return float(amplitude)


def seed_pattern(facts, base_ham_k, O_k, *, T_run, W0, amplitude=0.5,
                 seed_base: Literal["para_hf", "bare"]) -> SeedResult:
    """Seed the pattern ``O_k`` (``None``: the para row) on the cell of
    ``facts.problem`` (design section 7).

    Parameters
    ----------
    facts : :class:`~hwave.initial.facts.ModelFacts`
        Loaded on the candidate cell (its problem's blocks, mu-groups,
        ``Ncond``/``2Sz`` and Fermi cut are used).
    base_ham_k : array ``(nvol, nd, nd)``
        ``H_base`` in the solver's index order: the para-HF Hamiltonian
        folded to the cell (``seed_base = "para_hf"``) or the bare
        ``T(k)`` (``"bare"``).
    O_k : array ``(nvol, nd, nd)`` or None
        The Hermitian, Hamiltonian-convention pattern
        (:func:`hwave.initial.patterns.realise`).
    T_run : float
        The temperature of the run (the input T).
    W0 : float
        The active-window span of the T_hint helper (scales the T = 0
        shell tolerance only).
    amplitude : float in (0, 1)
        The fraction of ``m_inf - m(0)`` to reach.
    seed_base : ``"para_hf"`` or ``"bare"``

    Returns
    -------
    SeedResult
    """
    amplitude = _check_amplitude(amplitude)
    if seed_base not in SEED_BASES:
        raise ValueError(f"seed_base must be one of {SEED_BASES}, got {seed_base!r}")
    p = facts.problem
    base = np.asarray(base_ham_k)
    if base.shape != (p.nvol, p.nd, p.nd):
        raise ValueError(f"base_ham_k of this cell has shape {(p.nvol, p.nd, p.nd)}, "
                         f"got {base.shape}")

    if O_k is None:                                         # para
        G0, info0 = density_at(facts, base, None, 0.0, T_run=T_run, W0=W0)
        shell0 = _shell(info0)
        return SeedResult(tried=True, G=_read_only(G0), seed_base=seed_base,
                          shell={"h0": shell0, "h_used": shell0},
                          moments_h0=moments_of(G0, facts), h_used=0.0)

    o = _as_pattern(O_k, p.nvol, p.nd)
    norm = float(np.linalg.norm(o))
    if norm == 0.0:
        raise ValueError("seed_pattern: the pattern is zero (the para row passes O_k=None)")
    if float(np.linalg.norm(o - np.conj(np.swapaxes(o, -1, -2)))) > 1e-10 * norm:
        raise ValueError("seed_pattern: the pattern O(k) is not Hermitian")

    msgs = block_check(o, p.blocks, ns=p.ns, norb=p.norb)
    if msgs:
        return SeedResult(tried=False,
                          not_tried_reason="{}: {}".format(BLOCK_REASON, "; ".join(msgs)))

    G0, info0 = density_at(facts, base, o, 0.0, T_run=T_run, W0=W0)
    m0 = m_value(o, G0, problem=p)
    if constraint_mode(o, facts):
        return SeedResult(tried=False, not_tried_reason=CONSTRAINT_REASON,
                          seed_base=seed_base, m0=m0)
    shell0 = _shell(info0)
    moments = moments_of(G0, facts)
    minf = m_inf(o, p.blocks, bvol=_bvol(p))
    dm_sat = minf - m0
    if dm_sat < NO_RESPONSE_TOL:
        return SeedResult(tried=False, not_tried_reason=DEGENERATE_REASON,
                          seed_base=seed_base, shell={"h0": shell0, "h_used": None},
                          moments_h0=moments, m0=m0, m_inf=minf)

    target = amplitude * dm_sat
    tol = BISECTION_RTOL * dm_sat
    o_norm = op_norm(o)

    def evaluate(h):
        G, info = density_at(facts, base, o, h, T_run=T_run, W0=W0, o_norm=o_norm)
        return m_value(o, G, problem=p) - m0, G, info

    h = w_floor(base)
    h_lo = 0.0
    found = None
    for j in range(_CAP_DOUBLINGS + 1):
        d, G, info = evaluate(h)
        if d >= target:
            found = (h, d, G, info)
            break
        if j == _CAP_DOUBLINGS:
            break
        h_lo = h
        h *= 2.0
    cap_reason = None
    if found is None:                                       # the numerical guard
        h_hi = h_used = h
        d_used, G_used, info_used = d, G, info
        cap_reason = CAP_REASON
    else:
        h_hi = found[0]
        h_used, d_used, G_used, info_used = found
        if abs(d_used - target) > tol:
            a, b = h_lo, h_hi
            for _ in range(BISECTION_STEPS):
                mid = 0.5 * (a + b)
                d, G, info = evaluate(mid)
                if abs(d - target) <= tol:
                    h_used, d_used, G_used, info_used = mid, d, G, info
                    break
                if d >= target:
                    b = mid
                    h_used, d_used, G_used, info_used = mid, d, G, info
                else:
                    a = mid
    return SeedResult(
        tried=True, G=_read_only(G_used), seed_base=seed_base,
        shell={"h0": shell0, "h_used": _shell(info_used)}, moments_h0=moments,
        m0=m0, m_inf=minf, m_target=m0 + target, m_used=m0 + d_used, h_hi=h_hi,
        h_used=h_used, amplitude_requested=amplitude, amplitude_used=d_used / dm_sat,
        cap_reason=cap_reason)
