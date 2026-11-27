"""Instability hints: the hint temperature ``T_hint`` and the paramagnetic
Hartree-Fock reference.

:func:`t_hint_bare` sets the temperature at which the static para hints are
evaluated, from the bare one-body Hamiltonian ``H0 = T(k)`` and the k mesh
only (never from a result):

* a T > 0 input: ``T_hint = T`` (``t_hint_source = "input"``);
* a T = 0 input: ``T_hint = max(W0/(2*L_min), 1e-4*W0)`` with ``W0 =
  max(W_act, 1e-3)`` (``"rule"``); ``L_min`` is the smallest ``CellShape``
  axis longer than 2 (an axis of length 1 or 2 is structural: a ladder leg,
  a bilayer, a slab). With no such axis ``T_hint = W0/8`` and a warning;
* ``t_user`` (``--t-hint``) overrides both (``"user"``).

``W_act`` is the active-window span. Per mu-group g (one group with a free
Sz; the up and down sectors with a fixed ``2Sz``), on that group's
eigenvalues only:

* bands are energy ranks of the original (primitive) cell on the
  ``CellShape`` mesh: band n is the n-th eigenvalue of the group at each k
  (all its indices merged, sorted ascending). The folded cell of a
  ``SubShape`` has the same spectrum but more ranks, so the primitive cell
  keeps the result independent of ``SubShape``;
* ``E_F,g`` is the midpoint of the ``n_g``-th and ``(n_g+1)``-th sorted
  eigenvalues over the mesh;
* the Fermi bands are the bands whose range contains ``E_F,g``; with none
  (a gap at ``E_F,g``) the highest band below and the lowest above;
  ``D_g`` is the span of the Fermi bands' eigenvalues;
* ``W_act,g`` is the span of every band with an eigenvalue in
  ``[E_F,g - D_g, E_F,g + D_g]``.

An empty or full group has no window (its entries are ``None``);
``W_act`` is the largest ``W_act,g``, or the all-band span ``W_all`` when
every group is empty or full. Every quantity is shift-invariant: ``H0 ->
H0 + C`` moves each ``E_F,g`` by C and nothing else.

The band-range and window membership tests are inclusive with a
tolerance of ``1e-12 * W_all,g`` (the span of the group's eigenvalues; no
absolute energy enters, so it is shift-invariant), so a tie at an edge
does not depend on the rounding of the absolute energy origin. The
bands come from :func:`hwave.solver.initial_facade.primitive_h0_k`, the
solver's own ``H(k)`` builder.

Warnings: the wide-window warning is evaluated per mu-group, ``W_all,g >
3*W_act,g`` (strictly greater) for a group with a window, so the splitting
between the sectors of a fixed ``2Sz`` never triggers it; it is recorded
for every ``t_hint_source``, ``input`` and ``user`` included, because ``W0``
is also the scale of the field seeding; the no-mesh-axis warning only when
the ``W0/8`` cap is the ``T_hint`` used.

``T_hint`` affects the hint eigenvalues and eigenvectors, hence which hint
rows fill the instability quota; never the q set, the cells, or the
reserved, baseline, user and literature rows. The UHFk candidate runs stay
at the input T.

**The para reference.** :func:`para_hf` runs UHFk on the primitive cell at
``T_hint`` from the free density (:class:`ParaResult`): its Hamiltonian
``H_para = H0 + Sigma'[G_para]`` gives the bands of the hints, and, folded
to a candidate cell (:func:`fold_para`), the base Hamiltonian of the seeds
(``seed_base = "para_hf"``). The reference is accepted only when the SCF
converged and the state keeps the symmetries of the model (residual at most
:data:`PARA_RESIDUAL_TOL`): the spin rotations about every axis that ``H0``
and the Hartree-Fock map both keep (a Zeeman field along z keeps only the
rotations about z, and the z polarisation it induces is recorded, not
checked); with a fixed 2Sz the transverse spin and the electron count of
each sector; in spin-orbital mode time reversal when the model has it.
Otherwise the hints are skipped with the reason and the seeds use the bare
``H0``. :func:`para_from_green` rebuilds a reference from a saved para
Green function; :func:`hint_record` writes the ``hint`` record of the
coverage report.
"""
import dataclasses
import json
import logging
import math
import numbers
import os
import shutil
import tempfile
from dataclasses import dataclass
from types import MappingProxyType
from typing import Literal, Mapping, Optional, Tuple

import numpy as np

from hwave.solver import initial_facade

logger = logging.getLogger("qlms").getChild("initial").getChild("hints")

__all__ = ["THint", "t_hint_bare", "W_FLOOR_ACT", "T_HINT_REL_FLOOR",
           "STRUCTURAL_AXIS_MAX", "WIDE_WINDOW_FACTOR", "group_bands",
           "ParaResult", "para_hf", "para_from_green", "fold_para", "hint_record", "model_key",
           "PARA_RESIDUAL_TOL", "PARA_CHECK_TOL", "PARA_SYMMETRY_TOL"]

#: the floor of ``W0`` in input energy units (a single flat level)
W_FLOOR_ACT = 1e-3
#: ``T_hint >= T_HINT_REL_FLOOR * W0``
T_HINT_REL_FLOOR = 1e-4
#: a ``CellShape`` axis of at most this length is structural, not a k mesh
STRUCTURAL_AXIS_MAX = 2
#: the wide-window warning fires when ``W_all,g > WIDE_WINDOW_FACTOR * W_act,g``
WIDE_WINDOW_FACTOR = 3.0

#: the relative tolerance (of the group span ``W_all,g``) of the band-range
#: and window membership tests, so that a tie at a window edge does not depend
#: on the rounding of the absolute energy origin
_MEMBER_RTOL = 1e-12

RULE_MESH = "T_hint = max(W0/(2*L_min), 1e-4*W0)"
RULE_CAP = "T_hint = W0/8 (no k-mesh axis longer than 2)"
RULE_INPUT = "T_hint = T (the input temperature)"
RULE_USER = "T_hint = t_user (given with --t-hint)"

WIDE_WARNING = ("wide Wannier window: T_hint uses the active-window span W_act = {}, "
                "not the all-band span W_all = {}")
NO_MESH_WARNING = ("no k-mesh axis longer than 2: T_hint = W0/8 is a fixed cap, not a "
                   "mesh-derived scale; give --t-hint")

_EFFECT = (
    "The static hints are evaluated at T_hint = {:.6g}: each transition enters chi0 "
    "with weight (f(e) - f(e'))/(e' - e) <= 1/(4*T_hint), so structure within about "
    "T_hint of E_F or from energy differences below T_hint (small gaps, nesting "
    "mismatch, split van Hove peaks) is suppressed: lambda drops, peaks broaden, q* "
    "can shift and close rows can swap; lambda < 1 at T_hint is normal. T_hint "
    "affects lambda and the hint eigenvectors, hence which hint rows fill the "
    "instability quota; never the q set, the cells, or the reserved, baseline, user "
    "and literature rows. The UHFk candidate runs stay at the input T = {:.6g}.")


@dataclass(frozen=True)
class THint:
    """The hint temperature and its inputs.

    ``T_hint``; ``rule`` (the formula used, a sentence); ``W0 = max(W_act,
    1e-3)``, ``W_act`` (the active-window span), ``W_all`` (the all-band
    span of ``H0``). Per mu-group, in the solver's group order: ``E_F``,
    ``D``, ``W_act_groups`` (``None`` for an empty or full group).
    ``W_source`` (``"active_window"`` or ``"floor"``), ``L_min`` (the
    smallest ``CellShape`` axis longer than 2, ``None`` when there is
    none), ``t_hint_rule`` (``"mesh"`` or ``"no_mesh_axis_cap"``; ``None``
    when ``t_hint_source`` is ``"input"`` or ``"user"``), ``t_hint_source``
    (``"input"``, ``"rule"`` or ``"user"``), ``effect`` (the sentence read
    aloud) and ``warnings``.
    """
    T_hint: float
    rule: str
    W0: float
    W_act: float
    W_all: float
    E_F: Tuple[Optional[float], ...]
    D: Tuple[Optional[float], ...]
    W_act_groups: Tuple[Optional[float], ...]
    W_source: Literal["active_window", "floor"]
    L_min: Optional[int]
    t_hint_rule: Optional[Literal["mesh", "no_mesh_axis_cap"]]
    t_hint_source: Literal["input", "rule", "user"]
    effect: str
    warnings: Tuple[str, ...]


def _primitive_groups(problem, nd_orig):
    """The mu-group of each primitive index: the group of its image in the
    first primitive cell of the folded cell (folded orbital ``a +
    norb_orig * ir`` at ``ir = 0``)."""
    blocks = problem.blocks
    group_of = np.empty(problem.nd, dtype=int)
    for blk, g in zip(blocks.block_info, blocks.block_to_group):
        group_of[np.asarray(blk)] = g
    norb, norb_orig = problem.norb, problem.norb_orig
    if problem.enable_spin_orbital:
        folded = np.arange(nd_orig)
    else:
        prim = np.arange(nd_orig)
        folded = (prim // norb_orig) * norb + prim % norb_orig
    groups = group_of[folded]
    # every folded index of a group is an image of a primitive index of it
    subvol = int(np.prod(problem.subshape))
    for g in range(len(blocks.group_nconds)):
        if np.count_nonzero(group_of == g) != subvol * np.count_nonzero(groups == g):
            raise ValueError("t_hint: the mu-groups of the folded cell are not "
                             "images of primitive-cell groups")
    return groups


def group_bands(problem):
    """The bare bands of each mu-group of ``problem``.

    The bands are those of the original (primitive) cell on the
    ``CellShape`` mesh, so they do not depend on ``SubShape`` (the folded
    cell has the same spectrum but ``prod(SubShape)`` times as many energy
    ranks). Returns one ``(nvol_prim, nb_g)`` array per group (the
    solver's group order): at each k the eigenvalues of ``H0(k)`` on the
    group's indices, sorted ascending (band n is the n-th energy rank).
    No two groups are coupled by ``H0`` (they are unions of the solver's
    blocks).
    """
    h0 = initial_facade.primitive_h0_k(problem)
    groups = _primitive_groups(problem, h0.shape[1])
    out = []
    for g in range(len(problem.blocks.group_nconds)):
        idx = np.nonzero(groups == g)[0]
        sub = h0[:, idx[:, None], idx[None, :]]
        out.append(np.sort(np.linalg.eigvalsh(sub), axis=1))
    return out


def _group_window(eps, n):
    """``(E_F, D, W_act)`` of one group with bands ``eps`` ``(nvol, nb)``
    and ``n`` electrons; ``None`` for an empty or full group."""
    flat = np.sort(eps, axis=None)
    if n <= 0 or n >= flat.size:
        return None
    tol = _MEMBER_RTOL * float(flat[-1] - flat[0])
    e_f = 0.5 * (float(flat[n - 1]) + float(flat[n]))
    lo = eps.min(axis=0)
    hi = eps.max(axis=0)
    fermi = (lo <= e_f + tol) & (hi >= e_f - tol)
    if not fermi.any():
        below = np.nonzero(hi < e_f)[0]
        above = np.nonzero(lo > e_f)[0]
        if below.size:
            fermi[below.max()] = True
        if above.size:
            fermi[above.min()] = True
    d = float(hi[fermi].max() - lo[fermi].min())
    in_window = ((eps >= e_f - d - tol) & (eps <= e_f + d + tol)).any(axis=0)
    in_window |= fermi
    w_act = float(hi[in_window].max() - lo[in_window].min())
    return e_f, d, w_act


def _l_min(cellshape):
    axes = [int(n) for n in cellshape if int(n) > STRUCTURAL_AXIS_MAX]
    return min(axes) if axes else None


def _check_t_user(t_user):
    if (isinstance(t_user, bool) or not isinstance(t_user, numbers.Real)
            or not math.isfinite(float(t_user)) or float(t_user) <= 0.0):
        raise ValueError("t_hint: t_user must be a finite number > 0, got {!r}"
                         .format(t_user))
    return float(t_user)


def t_hint_bare(problem, *, t_user: Optional[float] = None) -> THint:
    """The hint temperature of ``problem`` (a ``UHFkProblem``) from its bare
    bands and k mesh (the rules of the module docstring).

    ``t_user``, when given, is the ``T_hint`` used (also for a T > 0
    input; the runs stay at the input T); it must be finite and > 0
    (``ValueError``). Otherwise ``T_hint = T`` for a T > 0 input and the
    mesh rule for a T = 0 input. ``W0``, ``W_act``, the per-group windows
    and ``L_min`` are computed in every case.
    """
    if t_user is not None:
        t_user = _check_t_user(t_user)
    bands = group_bands(problem)
    w_all = float(max(e.max() for e in bands) - min(e.min() for e in bands))
    e_fs, ds, w_groups, wide = [], [], [], []
    for eps, n in zip(bands, problem.blocks.group_nconds):
        res = _group_window(eps, int(n))
        if res is None:
            e_fs.append(None); ds.append(None); w_groups.append(None)
            continue
        e_fs.append(res[0]); ds.append(res[1]); w_groups.append(res[2])
        w_all_g = float(eps.max() - eps.min())
        if w_all_g > WIDE_WINDOW_FACTOR * res[2]:
            text = WIDE_WARNING.format("{:.6g}".format(res[2]), "{:.6g}".format(w_all_g))
            if text not in wide:
                wide.append(text)
    active = [w for w in w_groups if w is not None]
    w_act = max(active) if active else w_all
    w0 = max(w_act, W_FLOOR_ACT)
    w_source = "active_window" if w_act >= W_FLOOR_ACT else "floor"
    l_min = _l_min(problem.cellshape)

    warnings = list(wide)
    t_hint_rule = None
    if t_user is not None:
        t_hint, rule, source = t_user, RULE_USER, "user"
    elif problem.T > 0:
        t_hint, rule, source = float(problem.T), RULE_INPUT, "input"
    elif l_min is not None:
        t_hint = max(w0 / (2 * l_min), T_HINT_REL_FLOOR * w0)
        rule, source, t_hint_rule = RULE_MESH, "rule", "mesh"
    else:
        t_hint = w0 / 8
        rule, source, t_hint_rule = RULE_CAP, "rule", "no_mesh_axis_cap"
        warnings.append(NO_MESH_WARNING)

    return THint(
        T_hint=float(t_hint), rule=rule, W0=float(w0), W_act=float(w_act),
        W_all=w_all, E_F=tuple(e_fs), D=tuple(ds), W_act_groups=tuple(w_groups),
        W_source=w_source, L_min=l_min, t_hint_rule=t_hint_rule,
        t_hint_source=source, effect=_EFFECT.format(t_hint, float(problem.T)),
        warnings=tuple(warnings))


# --------------------------------------------------------------------------
# the paramagnetic Hartree-Fock reference
# --------------------------------------------------------------------------

#: the para state is accepted when its symmetry residual is at most this
PARA_RESIDUAL_TOL = 1e-8
#: the two para starts of the para check disagree when their free energies
#: per primitive cell differ by more than this
PARA_CHECK_TOL = 1e-8
#: a symmetry of ``H0`` or of the Hartree-Fock map holds to this (relative)
PARA_SYMMETRY_TOL = 1e-12

PARA_CHECK_WARNING = (
    "para check: the paramagnetic Hartree-Fock runs from the free density and from "
    "the uniform orbital occupation reach free energies per cell {:.12g} and {:.12g} "
    "(difference {:.3g} > {:.0e}); the lower ({}) is kept: the other start is a "
    "metastable (orbital-ordered) state")
PARA_CHECK_FAILED_WARNING = (
    "para check: the run from the uniform orbital occupation did not give a "
    "paramagnetic state ({}); the run from the free density is kept")

_AXES = ("x", "y", "z")
_SIGMA_P = (np.array([[0, 1], [1, 0]], dtype=np.complex128),
            np.array([[0, -1j], [1j, 0]], dtype=np.complex128),
            np.array([[1, 0], [0, -1]], dtype=np.complex128))


# eq=False: the fields hold arrays
@dataclass(frozen=True, eq=False)
class ParaResult:
    """The paramagnetic Hartree-Fock reference of a model.

    ``status``: ``"converged"`` (the SCF converged and the symmetry residual
    is at most :data:`PARA_RESIDUAL_TOL`), ``"not_converged"`` (the SCF did
    not converge, or its state breaks a symmetry the reference must keep)
    or ``"skipped"`` (the step did not run: no interacting orbital, or the
    solver refused). Only a converged reference is used; otherwise the
    hints are skipped with ``skip_reason`` and the seeds use the bare
    ``H0``.

    ``G_para``: the stored Green function of the primitive cell
    (``SubShape = [1, 1, 1]``). ``H_para_k``: ``H0(k) + Sigma'[G_para]``
    ``(prod(CellShape), nd_orig, nd_orig)`` (``H0`` from
    :func:`initial_facade.primitive_h0_k`, ``Sigma'`` the solver's
    block-projected potential with its Fock flag). ``bands``:
    :func:`initial_facade.bands_of` of ``H_para_k`` at ``T``. ``T``: the
    temperature of the reference (``T_hint``). ``free_energy``: the
    solver's ``FreeEnergy_Total`` (whole lattice); ``free_energy_per_cell``
    per primitive cell. ``residual``: the symmetry residual, the largest
    of ``residual_parts`` (keys among ``"spin"``, ``"transverse spin"``,
    ``"sector counts"``, ``"time reversal"``); ``checked``: what it covers
    (a spin axis, or time reversal, is included only when ``H0`` and the
    Hartree-Fock map both keep it). ``polarisation``: the polarisation the
    reference is allowed to have, recorded (``"spin"``: the on-site spin
    density ``[x, y, z]`` of each physical orbital; ``"sector_counts"``:
    the electrons of each mu-group; ``"time_reversal_odd"``: the largest
    time-reversal-odd entry, spin-orbital mode; ``"allowed"``: the spin
    axes along which a polarisation is allowed by the symmetries of the
    model, or ``"sectors"`` with a fixed 2Sz). ``iterations``,
    ``scf_residual``: the SCF record; ``eps``: the ``EPS`` exponent of the
    run (the input's, raised to at least
    :data:`initial_facade.PARA_EPS_MIN`). ``free_energy_check``: the free
    energy per cell of the para-check start (``None`` without the check).
    ``warnings``: messages for the report. ``input_sha256``: of the model
    the reference belongs to; ``model``: its ``cellshape``, ``norb_orig``,
    ``ncond``, ``two_sz``, ``iflag_fock`` and ``spin_orbital`` (the TOML
    settings the input files do not fix).
    """
    status: Literal["converged", "not_converged", "skipped"]
    G_para: Optional[np.ndarray]
    H_para_k: Optional[np.ndarray]
    bands: Optional["initial_facade.Bands"]
    free_energy: Optional[float]
    residual: Optional[float]
    polarisation: Optional[Mapping]
    skip_reason: Optional[str]
    warnings: Tuple[str, ...]
    T: float
    free_energy_per_cell: Optional[float] = None
    residual_parts: Optional[Mapping] = None
    checked: Tuple[str, ...] = ()
    iterations: Optional[int] = None
    scf_residual: Optional[float] = None
    free_energy_check: Optional[float] = None
    input_sha256: Optional[Mapping] = None
    model: Optional[Mapping] = None
    eps: Optional[int] = None


def model_key(facts):
    """The settings of the model of ``facts`` a para reference must share
    besides its input files (:attr:`ParaResult.model`)."""
    p = facts.problem
    return MappingProxyType({
        "cellshape": tuple(int(n) for n in p.cellshape), "norb_orig": int(p.norb_orig),
        "ncond": p.ncond, "two_sz": p.two_sz, "iflag_fock": bool(p.iflag_fock),
        "spin_orbital": bool(p.enable_spin_orbital)})


def _primitive_facts(facts):
    from hwave.initial.facts import model_facts
    if tuple(facts.problem.subshape) == (1, 1, 1):
        return facts
    return model_facts(facts.toml_path, subshape=(1, 1, 1))


def _on_site_spin(G, p):
    """``[x, y, z]`` of the on-site spin density of each physical orbital,
    from the stored G (``S^mu = 1/2 sum_st sigma^mu_st <c^dag_s c_t>``)."""
    g0 = np.asarray(G)[0].reshape(p.nd, p.nd)
    out = []
    nphys = p.norb_phys
    for a in range(nphys):
        if p.enable_spin_orbital:
            idx = [initial_facade.interleaved_index(a, s) for s in range(2)]
        else:
            idx = [s * p.norb + a for s in range(2)]
        blk = g0[np.ix_(idx, idx)]
        out.append([0.5 * float(np.real(np.sum(sig * blk))) for sig in _SIGMA_P])
    return np.array(out)


def _group_counts(G, p):
    g0 = np.asarray(G)[0].reshape(p.nd, p.nd)
    counts = []
    for g in range(len(p.blocks.group_nconds)):
        idx = np.concatenate([np.asarray(b) for b, gg in
                              zip(p.blocks.block_info, p.blocks.block_to_group) if gg == g])
        counts.append(float(p.nvol * np.sum(np.real(g0[idx, idx]))))
    return counts


def _symmetric_axes(p):
    """The spin axes that ``H0`` and the Hartree-Fock map both keep
    (normal mode, free Sz)."""
    from hwave.initial import symmetry
    h0 = symmetry.h0_spin_commutators(p)
    hf = symmetry.hf_axis_residuals(p)
    return tuple(mu for mu in range(3)
                 if h0[mu] <= PARA_SYMMETRY_TOL and hf[mu] <= PARA_SYMMETRY_TOL)


def _residual(G, p):
    """``(residual_parts, checked, polarisation)`` of the stored G of the
    primitive-cell problem ``p``."""
    from hwave.initial import symmetry
    spin = _on_site_spin(G, p)
    parts, checked = {}, []
    pol = {"spin": spin.tolist(), "sector_counts": None, "time_reversal_odd": None,
           "allowed": None}
    if p.enable_spin_orbital:
        tr = symmetry.tr_check(p)
        odd = float(np.max(np.abs(np.asarray(G) - symmetry.time_reverse_green(G, p)))) / 2
        pol["time_reversal_odd"] = odd
        if tr.invariant:
            parts["time reversal"] = odd
            checked.append("time reversal")
            pol["allowed"] = []
        else:
            pol["allowed"] = ["time-reversal odd (the model breaks time reversal)"]
        return parts, tuple(checked), pol
    if not p.sz_free:
        counts = _group_counts(G, p)
        pol["sector_counts"] = counts
        pol["allowed"] = ["sectors"]
        parts["transverse spin"] = float(np.max(np.hypot(spin[:, 0], spin[:, 1])))
        parts["sector counts"] = max(
            abs(c - n) / max(1.0, float(n)) for c, n in zip(counts, p.blocks.group_nconds))
        return parts, ("transverse spin", "sector counts"), pol
    axes = _symmetric_axes(p)
    if len(axes) >= 2:
        allowed = ()
    elif len(axes) == 1:
        allowed = axes
    else:
        allowed = (0, 1, 2)
    forbidden = [mu for mu in range(3) if mu not in allowed]
    pol["allowed"] = [_AXES[mu] for mu in allowed]
    if forbidden:
        parts["spin"] = float(np.max(np.linalg.norm(spin[:, forbidden], axis=1)))
        checked.extend("S_" + _AXES[mu] for mu in forbidden)
    return parts, tuple(checked), pol


def _start_free(facts1, thint):
    """The free density of the primitive cell at ``T_hint``: the bare-H0
    fill with the solver's mu-groups."""
    from hwave.initial.seed import seed_pattern
    p = facts1.problem
    return np.asarray(seed_pattern(facts1, initial_facade.primitive_h0_k(p), None,
                                   T_run=thint.T_hint, W0=thint.W0, seed_base="bare").G)


def _start_uniform(p):
    """The uniform orbital occupation: on the r = 0 diagonal of every index
    of mu-group g the value ``n_g / (nvol * dim_g)``, zero elsewhere."""
    g = np.zeros((p.nvol, p.nd, p.nd), dtype=np.complex128)
    for grp, n in enumerate(p.blocks.group_nconds):
        idx = np.concatenate([np.asarray(b) for b, gg in
                              zip(p.blocks.block_info, p.blocks.block_to_group) if gg == grp])
        g[0, idx, idx] = float(n) / (p.nvol * idx.size)
    return g.reshape(p.nvol, p.ns, p.norb, p.ns, p.norb)


def para_from_green(facts, thint, G_para, *, free_energy, converged, iterations=None,
                    scf_residual=None, eps=None, warnings=()) -> ParaResult:
    """The :class:`ParaResult` of a para Green function ``G_para`` of the
    primitive cell (from :func:`para_hf`, or read back from the green file
    of a paramagnetic run of the input at ``thint.T_hint``).

    ``H_para_k = H0 + Sigma'[G_para]``, the bands, the symmetry residual
    and the polarisation are rebuilt from ``G_para``; ``free_energy`` is
    that run's ``FreeEnergy_Total``, ``converged`` its SCF status and
    ``eps`` its ``EPS`` exponent."""
    facts1 = _primitive_facts(facts)
    p = facts1.problem
    G = np.array(G_para, dtype=np.complex128)
    want = (p.nvol, p.ns, p.norb, p.ns, p.norb)
    if G.shape != want:
        raise ValueError("para_from_green: G_para must be the primitive-cell Green function "
                         "{}, got {}".format(want, G.shape))
    G.flags.writeable = False
    H = initial_facade.primitive_h0_k(p) + initial_facade.hf_potential(
        G, initial_facade.hf_tables(p), include_fock=p.iflag_fock, blocks=p.blocks)
    H.flags.writeable = False
    T = float(thint.T_hint)
    bands = initial_facade.bands_of(H, p.blocks, T, shape=p.shape, ns=p.ns, norb=p.norb,
                                    ene_cutoff=p.ene_cutoff)
    parts, checked, pol = _residual(G, p)
    residual = max(parts.values()) if parts else 0.0
    ncell = int(np.prod(p.cellshape))
    f_cell = float(free_energy) / ncell
    if not converged:
        status = "not_converged"
        reason = ("the paramagnetic Hartree-Fock run did not converge at T_hint = {:.6g} "
                  "({} iterations, residual {})".format(
                      T, iterations, "unknown" if scf_residual is None
                      else "{:.3g}".format(scf_residual)))
    elif residual > PARA_RESIDUAL_TOL:
        status = "not_converged"
        worst = max(parts, key=parts.get)
        reason = ("the paramagnetic Hartree-Fock state breaks a symmetry of the model: "
                  "{} residual {:.3g} > {:.0e}".format(worst, residual, PARA_RESIDUAL_TOL))
    else:
        status, reason = "converged", None
    return ParaResult(
        status=status, G_para=G, H_para_k=H, bands=bands, free_energy=float(free_energy),
        residual=float(residual), polarisation=MappingProxyType(pol), skip_reason=reason,
        warnings=tuple(warnings), T=T, free_energy_per_cell=f_cell,
        residual_parts=MappingProxyType(dict(parts)), checked=checked,
        iterations=iterations, scf_residual=scf_residual,
        input_sha256=MappingProxyType(dict(facts.input_sha256)), model=model_key(facts),
        eps=None if eps is None else int(eps))


def _skipped(facts, thint, reason, warnings=()):
    return ParaResult(status="skipped", G_para=None, H_para_k=None, bands=None,
                      free_energy=None, residual=None, polarisation=None,
                      skip_reason=reason, warnings=tuple(warnings),
                      T=float(thint.T_hint),
                      input_sha256=MappingProxyType(dict(facts.input_sha256)),
                      model=model_key(facts))


def _one_run(facts1, thint, start, work_dir):
    run = initial_facade.run_para(facts1.toml_path, T=thint.T_hint, start_G=start,
                                  work_dir=work_dir)
    res = para_from_green(facts1, thint, run.G, free_energy=run.free_energy,
                          converged=run.converged, iterations=run.iterations,
                          scf_residual=run.residual, eps=run.eps)
    return run, res


def _record(res, selected, work_dir):
    out = {
        "status": res.status, "skip_reason": res.skip_reason, "T": res.T,
        "free_energy": res.free_energy, "free_energy_per_cell": res.free_energy_per_cell,
        "residual": res.residual,
        "residual_parts": None if res.residual_parts is None else dict(res.residual_parts),
        "checked": list(res.checked), "iterations": res.iterations,
        "scf_residual": res.scf_residual, "eps": res.eps,
        "free_energy_check": res.free_energy_check,
        "polarisation": None if res.polarisation is None else dict(res.polarisation),
        "selected_start": selected, "warnings": list(res.warnings),
        "input_sha256": dict(res.input_sha256 or {}),
        "model": None if res.model is None else
        {k: (list(v) if isinstance(v, tuple) else v) for k, v in res.model.items()},
    }
    with open(os.path.join(work_dir, "para_hf.json"), "w") as f:
        json.dump(out, f, indent=1, sort_keys=True, allow_nan=True)
        f.write("\n")


def para_hf(facts, thint, *, para_check=False, work_dir=None) -> ParaResult:
    """The paramagnetic Hartree-Fock reference of the model of ``facts``.

    UHFk runs on the primitive cell (``SubShape = [1, 1, 1]``) at
    ``thint.T_hint`` with the input's ``Ncond``, ``2Sz``, Fock flag,
    tables and mu-groups, through :func:`initial_facade.run_para`. It
    starts from the free density (the bare-``H0`` fill at ``T_hint``; with
    a fixed 2Sz the free density of each sector), which the Hartree-Fock
    map keeps paramagnetic. The result is accepted when the SCF converged
    and the symmetry residual (:class:`ParaResult`) is at most
    :data:`PARA_RESIDUAL_TOL`.

    ``para_check`` runs a second start, the uniform orbital occupation, and
    keeps the lower free energy; a difference per primitive cell above
    :data:`PARA_CHECK_TOL` is a warning (a metastable orbital-ordered
    state passes the residual).

    ``work_dir``: the directory of the runs (``main/``, and ``check/`` with
    the check), kept; the selected run's Green function is copied to
    ``work_dir/green.npz`` and the result recorded in
    ``work_dir/para_hf.json``. Without it the runs use a temporary
    directory that is removed.

    A model without an interacting orbital, or a run the solver refuses,
    gives ``status = "skipped"`` with the reason (no exception).
    """
    if facts.n_c == 0:
        return _skipped(facts, thint, "no interacting orbital (n_c = 0)")
    facts1 = _primitive_facts(facts)
    tmp = None
    if work_dir is None:
        tmp = tempfile.mkdtemp(prefix="hwave-para-")
        root = tmp
    else:
        root = os.path.abspath(os.fspath(work_dir))
        os.makedirs(root, exist_ok=True)
    try:
        try:
            run, res = _one_run(facts1, thint, _start_free(facts1, thint),
                                os.path.join(root, "main"))
        except (ValueError, initial_facade.MuNotConvergedError) as exc:
            return _skipped(facts, thint, "the paramagnetic Hartree-Fock run failed: "
                            "{}".format(exc))
        selected = "main"
        if para_check:
            res, run, selected = _check(facts1, thint, res, run, root)
        res = dataclasses.replace(res, input_sha256=MappingProxyType(dict(facts.input_sha256)))
        if work_dir is not None:
            shutil.copyfile(os.path.join(run.work_dir, "green.npz"),
                            os.path.join(root, "green.npz"))
            _record(res, selected, root)
        return res
    finally:
        if tmp is not None:
            shutil.rmtree(tmp, ignore_errors=True)


def _check(facts1, thint, res, run, root):
    p = facts1.problem
    try:
        run2, res2 = _one_run(facts1, thint, _start_uniform(p), os.path.join(root, "check"))
    except (ValueError, initial_facade.MuNotConvergedError) as exc:
        msg = PARA_CHECK_FAILED_WARNING.format(exc)
        logger.warning(msg)
        return dataclasses.replace(res, warnings=res.warnings + (msg,)), run, "main"
    if res2.status != "converged":
        msg = PARA_CHECK_FAILED_WARNING.format(res2.skip_reason)
        logger.warning(msg)
        return (dataclasses.replace(res, warnings=res.warnings + (msg,),
                                    free_energy_check=res2.free_energy_per_cell),
                run, "main")
    if res.status != "converged":
        # the check start reached a valid reference where the free start did not
        return dataclasses.replace(res2, free_energy_check=res2.free_energy_per_cell), \
            run2, "check"
    f1, f2 = res.free_energy_per_cell, res2.free_energy_per_cell
    keep, krun, sel = (res2, run2, "check") if f2 < f1 else (res, run, "main")
    warnings = keep.warnings
    if abs(f1 - f2) > PARA_CHECK_TOL:
        msg = PARA_CHECK_WARNING.format(
            f1, f2, abs(f1 - f2), PARA_CHECK_TOL,
            "uniform orbital occupation" if sel == "check" else "free density")
        logger.warning(msg)
        warnings = warnings + (msg,)
    return dataclasses.replace(keep, warnings=warnings, free_energy_check=f2), krun, sel


def fold_para(H_para_k, facts, subshape):
    """``H_para_k`` (the primitive-cell para Hamiltonian of
    :class:`ParaResult`) folded onto the ``SubShape`` ``subshape`` cell of
    the model of ``facts``: ``(nvol_B, nd_B, nd_B)`` in the solver's folded
    index order (:func:`initial_facade.fold_hamiltonian_k`), the base
    Hamiltonian of the seeds on that cell."""
    p = facts.problem
    return initial_facade.fold_hamiltonian_k(
        H_para_k, cellshape=facts.cellshape, subshape=tuple(int(b) for b in subshape),
        norb_orig=p.norb_orig, ns=p.ns)


def hint_record(thint, para=None, *, skip_reason=None):
    """The ``hint`` record of ``coverage.json`` for the hint temperature
    ``thint`` and the para reference ``para`` (``None``: no para step).

    The record is skipped with ``skip_reason`` when it is given, with
    ``para.skip_reason`` when the para reference is not converged, and when
    there is neither a reason nor a para reference ``skip_reason`` must be
    given. A skipped record has every numeric para field null (its
    ``para_status`` is kept); the vertex fields (``gamma_self_adjoint``,
    ``method``, ``q_status``, ``q_set``) are left empty for the eigen-solve
    to fill."""
    if skip_reason is None and para is None:
        raise ValueError("hint_record: a record without a para reference needs skip_reason")
    if skip_reason is None and para.status != "converged":
        skip_reason = para.skip_reason
    skipped = skip_reason is not None
    out = {
        "status": "skipped" if skipped else "computed", "skip_reason": skip_reason,
        "T_hint": thint.T_hint, "rule": thint.rule, "W0": thint.W0, "W_act": thint.W_act,
        "W_all": thint.W_all, "E_F": list(thint.E_F), "D": list(thint.D),
        "W_act_groups": list(thint.W_act_groups), "W_source": thint.W_source,
        "q_status": {}, "L_min": thint.L_min, "t_hint_rule": thint.t_hint_rule,
        "t_hint_source": thint.t_hint_source, "effect": thint.effect,
        "para_status": None if para is None else para.status,
        "para_free_energy": None, "para_residual": None, "para_polarisation": None,
        "para_eps": None,
        "gamma_self_adjoint": None, "method": None, "q_set": [],
    }
    if not skipped:
        out["para_free_energy"] = para.free_energy_per_cell
        out["para_residual"] = para.residual
        out["para_polarisation"] = _plain(dict(para.polarisation))
        out["para_eps"] = para.eps
    return out


def _plain(obj):
    if isinstance(obj, Mapping):
        return {k: _plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_plain(v) for v in obj]
    if isinstance(obj, np.generic):
        return obj.item()
    return obj
