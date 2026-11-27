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

**The vertex.** :func:`pair_space` lists the slots the solver's
Hartree-Fock map reads and writes (:class:`PairSpace`, with its label
coordinates), :func:`vertex` and :func:`vertex_sparse` give the Jacobian
``Gamma(q) = dSigma'/drho`` of that map on them (dense for the small tier,
CSR for the medium tier), and :func:`vertex_self_adjoint` checks that it is
Hermitian. :func:`lindhard` and :func:`lindhard_apply` give the static
Lindhard response ``chi0(q)`` of the para bands on the same pair space and
sectors (dense for the small tier, a matrix-free product for the medium
tier; one kernel, :func:`hwave.solver.initial_facade.static_lindhard`).
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
           "ParaResult", "para_hf", "para_from_green", "fold_para", "hint_record", "model_key", "input_digest",
           "PARA_RESIDUAL_TOL", "PARA_CHECK_TOL", "PARA_SYMMETRY_TOL",
           "PairSpace", "pair_space", "vertex", "vertex_sparse", "vertex_self_adjoint",
           "SECTOR_FULL", "SECTOR_CHARGE", "SECTOR_SZ", "VERTEX_SELF_ADJOINT_TOL"]

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
PARA_CHECK_SWITCHED_WARNING = (
    "para check: the run from the free density did not give a paramagnetic state "
    "({}); the run from the uniform orbital occupation is kept")
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
    :data:`initial_facade.PARA_EPS_MIN`), ``iteration_max`` its
    ``IterationMax`` (the input's). ``free_energy_check``: the free
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
    iteration_max: Optional[int] = None


def input_digest(input_sha256):
    """The SHA-256 of the canonical JSON (sorted keys, no whitespace) of an
    ``input_sha256`` mapping: one string naming the input files of a model
    (the para record and the seed provenance carry it)."""
    import hashlib
    text = json.dumps(dict(input_sha256), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


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
                    scf_residual=None, eps=None, iteration_max=None,
                    warnings=()) -> ParaResult:
    """The :class:`ParaResult` of a para Green function ``G_para`` of the
    primitive cell (from :func:`para_hf`, or read back from the green file
    of a paramagnetic run of the input at ``thint.T_hint``).

    ``H_para_k = H0 + Sigma'[G_para]``, the bands, the symmetry residual
    and the polarisation are rebuilt from ``G_para``; ``free_energy`` is
    that run's ``FreeEnergy_Total``, ``converged`` its SCF status and
    ``eps`` and ``iteration_max`` its ``EPS`` exponent and
    ``IterationMax``."""
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
                  "({} iterations, residual {}; EPS = {}, IterationMax = {}: the para "
                  "run converges to EPS >= 12, so it may need a larger IterationMax "
                  "than the input's)".format(
                      T, iterations, "unknown" if scf_residual is None
                      else "{:.3g}".format(scf_residual),
                      "unknown" if eps is None else eps,
                      "unknown" if iteration_max is None else iteration_max))
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
        eps=None if eps is None else int(eps),
        iteration_max=None if iteration_max is None else int(iteration_max))


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
                          scf_residual=run.residual, eps=run.eps,
                          iteration_max=run.iteration_max)
    return run, res


def _record(res, selected, work_dir):
    out = {
        "status": res.status, "skip_reason": res.skip_reason, "T": res.T,
        "free_energy": res.free_energy, "free_energy_per_cell": res.free_energy_per_cell,
        "residual": res.residual,
        "residual_parts": None if res.residual_parts is None else dict(res.residual_parts),
        "checked": list(res.checked), "iterations": res.iterations,
        "scf_residual": res.scf_residual, "eps": res.eps,
        "iteration_max": res.iteration_max,
        "free_energy_check": res.free_energy_check,
        "polarisation": None if res.polarisation is None else dict(res.polarisation),
        "selected_start": selected, "warnings": list(res.warnings),
        "input_sha256": dict(res.input_sha256 or {}),
        "input_digest": None if res.input_sha256 is None else input_digest(res.input_sha256),
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
    ``work_dir/para_hf.json`` (also for a skipped step, then without a
    green file and with ``selected_start = null``). Without it the runs use a temporary
    directory that is removed.

    A model without an interacting orbital, or a run the solver refuses,
    gives ``status = "skipped"`` with the reason (no exception).
    """
    root = None if work_dir is None else os.path.abspath(os.fspath(work_dir))
    if root is not None:
        os.makedirs(root, exist_ok=True)
    if facts.n_c == 0:
        res = _skipped(facts, thint, "no interacting orbital (n_c = 0)")
        if root is not None:
            _record(res, None, root)
        return res
    facts1 = _primitive_facts(facts)
    tmp = None
    if root is None:
        tmp = tempfile.mkdtemp(prefix="hwave-para-")
    base = root if tmp is None else tmp
    try:
        try:
            run, res = _one_run(facts1, thint, _start_free(facts1, thint),
                                os.path.join(base, "main"))
        except (ValueError, initial_facade.MuNotConvergedError) as exc:
            res = _skipped(facts, thint, "the paramagnetic Hartree-Fock run failed: "
                           "{}".format(exc))
            if root is not None:
                _record(res, None, root)
            return res
        selected = "main"
        if para_check:
            res, run, selected = _check(facts1, thint, res, run, base)
        res = dataclasses.replace(res, input_sha256=MappingProxyType(dict(facts.input_sha256)))
        if root is not None:
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
        msg = PARA_CHECK_SWITCHED_WARNING.format(res.skip_reason)
        logger.warning(msg)
        return dataclasses.replace(res2, warnings=res2.warnings + (msg,),
                                   free_energy_check=res2.free_energy_per_cell), \
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


#: the vertex fields of a computed hint record
_VERTEX_FIELDS = ("gamma_self_adjoint", "method", "q_status", "q_set")


def hint_record(thint, para=None, *, skip_reason=None, vertex=None):
    """The ``hint`` record of ``coverage.json`` for the hint temperature
    ``thint`` and the para reference ``para`` (``None``: no para step).

    The record is skipped with ``skip_reason`` when it is given, or with
    ``para.skip_reason`` when the para reference is not converged; a record
    without a para reference needs ``skip_reason``. A skipped record has
    every numeric para field null (its ``para_status`` is kept) and empty
    vertex fields. A computed record (a converged para and no skip reason)
    needs ``vertex``, the mapping of the eigen-solve's
    ``gamma_self_adjoint``, ``method``, ``q_status`` and ``q_set``
    (``ValueError`` otherwise)."""
    if skip_reason is None and para is None:
        raise ValueError("hint_record: a record without a para reference needs skip_reason")
    if skip_reason is None and para.status != "converged":
        skip_reason = para.skip_reason
    skipped = skip_reason is not None
    if not skipped and vertex is None:
        raise ValueError("hint_record: a computed record needs the vertex fields "
                         "({})".format(", ".join(_VERTEX_FIELDS)))
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
        missing = [k for k in _VERTEX_FIELDS if k not in vertex]
        if missing:
            raise ValueError("hint_record: vertex lacks {}".format(", ".join(missing)))
        for k in _VERTEX_FIELDS:
            out[k] = _plain(vertex[k])
    return out


def _plain(obj):
    if isinstance(obj, Mapping):
        return {k: _plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_plain(v) for v in obj]
    if isinstance(obj, np.generic):
        return obj.item()
    return obj


# ---------------------------------------------------------------------------
# The Hartree-Fock vertex on the pair space
# ---------------------------------------------------------------------------

#: the one sector of a pair space without a passed spin-rotation check
SECTOR_FULL = "full"
#: the charge and S^z sectors of a pair space with a passed spin-rotation check
SECTOR_CHARGE = "charge"
SECTOR_SZ = "S^z"
#: the relative tolerance of :func:`vertex_self_adjoint`
VERTEX_SELF_ADJOINT_TOL = 1e-12
#: the dense vertex is built in batches of columns of at most this many
#: slot-grid entries in total
_VERTEX_BATCH_ENTRIES = 2 ** 21
#: the Hartree-loop types of the solver's kernel, in its order (PairHop last)
_HF_LOOP = ("CoulombIntra", "CoulombInter", "Hund", "Ising", "PairLift", "Exchange")


# eq=False: the fields hold arrays
@dataclass(frozen=True, eq=False)
class PairSpace:
    """The pair space ``P`` of the Hartree-Fock vertex of a model.

    **Coordinates.** A label ``m = (a, s, b, t, r)`` is the matrix slot of a
    lattice one-body operator ``sum_R O_m c^dag_{as}(R + r) c_{bt}(R)``
    (``a``, ``b`` physical orbitals of the primitive cell, 0-based; ``s``,
    ``t`` the spins, 0 = up, explicit in every mode, the spin-orbital index
    ``2*a + s`` included; ``r`` the bond in primitive cells). Its Green
    function slot is ``G_m = <c^dag_{as}(r) c_{bt}(0)> = G[r, s, a, t, b]``
    of the stored layout, and ``<O> = sum_m O_m G_m`` per cell. A density
    vector ``x`` on ``P`` holds rho-components (``rho = G^T``: ``x_m =
    G_mbar = conj(G_m)`` for a Hermitian state, ``mbar = (b, t, a, s, -r)``);
    the vertex output and the conjugate field ``O = Gamma x`` hold
    potential components ``Sigma'_m``, the coefficient of the operator at
    ``m``. A q-modulated vector carries ``e^{iq.R}`` at the anchor ``R``,
    the cell of the annihilation operator. With these coordinates the
    vertex is Hermitian in the plain product ``sum_m conj(u_m) v_m``.

    ``labels``: ``(|P|, 7)`` int64 ``(a, s, b, t, rx, ry, rz)``, sorted by
    the slot ``(r, s, a, t, b)`` of the stored layout (``r`` as the FFT
    index of the k-mesh); ``r`` is the representative with ``-N/2 < r_i <=
    N/2``. ``P`` holds every slot the solver's Hartree-Fock map reads or
    writes (the Hartree terms always; the Fock terms and the second
    PairHop contraction with the model's Fock flag), restricted to the
    interacting orbitals and to slots inside one solver block.

    ``sectors``: ``("full",)``, or ``("charge", "S^z")`` when the
    spin-rotation check passed (normal mode, free S_z): the vertex is then
    block diagonal in the basis ``(e_up,up +- e_dn,dn)/sqrt 2`` of each
    orbital slot and the spin-flip labels, and by the spin rotations the
    ``S^x`` and ``S^y`` blocks have the nonzero spectrum of the ``S^z``
    block, so only charge and ``S^z`` are kept. ``sector_basis``: per sector the
    ``(|P|, d)`` CSR matrix of orthonormal columns (``None`` for the full
    space); ``sector_labels``: per sector ``(d, 6)`` int64 ``(a, mu, b, rx,
    ry, rz)`` with ``mu = 0`` (charge) or ``3`` (S^z), or ``None`` for the
    full space (its labels are ``labels``).

    ``shape``: the k-mesh (``CellShape``), ``norb_phys``, ``include_fock``
    (the model's Fock flag), ``spin_orbital``, ``su2_status`` (the status
    of the spin-rotation check used) and ``sector_note`` (why a passed
    check still gave one sector, else ``None``).
    """
    labels: np.ndarray
    sectors: Tuple[str, ...]
    sector_basis: Tuple[Optional[object], ...]
    sector_labels: Tuple[Optional[np.ndarray], ...]
    shape: Tuple[int, int, int]
    norb_phys: int
    include_fock: bool
    spin_orbital: bool
    su2_status: str
    sector_note: Optional[str] = None

    @property
    def size(self):
        """``|P|``."""
        return int(self.labels.shape[0])

    def sector_dim(self, sector):
        """The dimension of ``sector``."""
        k = _sector_index(self, sector)
        basis = self.sector_basis[k]
        return self.size if basis is None else int(basis.shape[1])


def _sector_index(P, sector):
    if sector is None:
        return None
    if sector not in P.sectors:
        raise ValueError("pair space: unknown sector {!r}; the sectors are {}"
                         .format(sector, ", ".join(P.sectors)))
    return P.sectors.index(sector)


def _grid_dims(problem):
    """``(nvol, 2, n, 2, n)``: the slot grid ``(r, s, a, t, b)`` of the
    primitive cell (``n`` physical orbitals)."""
    n = int(problem.tables.norb_phys)
    return (int(problem.nvol), 2, n, 2, n)


def _neg_r(shape):
    nx, ny, nz = shape
    ix, iy, iz = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    return (((-iz) % nz) + nz * (((-iy) % ny) + ny * ((-ix) % nx))).ravel()


def _r_vectors(shape):
    """The centred bond vector of each flat r index, ``(nvol, 3)``."""
    nx, ny, nz = shape
    idx = np.stack(np.unravel_index(np.arange(nx * ny * nz), shape), axis=1)
    n = np.array(shape)
    return np.where(idx <= n // 2, idx, idx - n).astype(np.int64)


def _phases(shape, q):
    """``e^{2 pi i q.r}`` for every flat r index (q on the mesh of
    ``shape``, reduced exactly before the exponential)."""
    from hwave.initial.cells import to_qfrac
    qf = to_qfrac(q)
    m = []
    for i, (x, n) in enumerate(zip(qf, shape)):
        y = x * int(n)
        if y.denominator != 1:
            raise ValueError("vertex: q = {} is not on the k-mesh {} (component {} "
                             "times {} is not an integer)".format(
                                 ",".join(str(v) for v in qf), list(shape), i, n))
        m.append(int(y))
    L = int(np.lcm.reduce([int(n) for n in shape]))
    idx = np.stack(np.unravel_index(np.arange(int(np.prod(shape))), shape), axis=1)
    num = sum(m[i] * idx[:, i] * (L // int(shape[i])) for i in range(3)) % L
    return np.exp(2j * np.pi * num / L)


def _nd_index(problem, a, s):
    """The solver index of physical orbital ``a``, spin ``s`` (arrays)."""
    if problem.enable_spin_orbital:
        return initial_facade.interleaved_index(a, s)
    return s * problem.norb + a


def _vertex_terms(problem):
    """The terms of the vertex, in the order of the solver's kernel (types
    in its order; per type the Hartree, then the Fock terms; per term the
    table entries in C order, then the spin entries): ``(out, inp, coef,
    rid, kind)``, the slot-grid flat indices of the output and the input
    label, the coefficient without its phase, the flat r of the table entry
    and the phase kind (``+1``: ``e^{+iq.r}``, ``0``: none, ``-1``:
    ``e^{-iq.r}``). Every term has a nonzero coefficient."""
    tables = problem.tables
    dims = _grid_dims(problem)
    nvol, n = dims[0], dims[2]
    neg = _neg_r(problem.shape)
    fock = bool(problem.iflag_fock)
    parts = []

    def add(r, J, w, out, inp, kind):
        # entries x spin entries, entry-major
        ne, ns_ = r.size, w.size
        if ne == 0 or ns_ == 0:
            return
        parts.append((np.ravel_multi_index(out, dims).reshape(ne * ns_),
                      np.ravel_multi_index(inp, dims).reshape(ne * ns_),
                      (J[:, None] * w[None, :]).reshape(ne * ns_),
                      np.repeat(r, ns_), np.full(ne * ns_, kind, dtype=np.int8)))

    def entries(name):
        tab = tables.inter_table.get(name)
        if tab is None:
            return None
        jab = np.asarray(tab).reshape(nvol, n, n)
        r, a, b = np.nonzero(jab)
        spin = np.asarray(tables.spin_table[name])
        sp = np.argwhere(spin != 0)
        w = spin[tuple(sp.T)].astype(np.complex128)
        return r, a, b, jab[r, a, b].astype(np.complex128), sp, w

    for name in _HF_LOOP:
        e = entries(name)
        if e is None:
            continue
        r, a, b, J, sp, w = e
        R, A, B = r[:, None], a[:, None], b[:, None]
        Z = np.zeros_like(R)
        p0, p1, p2, p3 = (sp[:, i][None, :] for i in range(4))
        # Hartree: Sigma'^{st}_{aa}(0) += J_ab(r) spin[s,u,v,t] rho(b v, b u; 0),
        # b sits at +r from a: e^{+iq.r}
        add(r, J, w, (Z, p0, A, p3, A), (Z, p2, B, p1, B), 1)
        if fock:
            # Fock: Sigma'(a s, b t; -r) -= J_ab(r) spin[s,u,t,v] rho(a v, b u; -r)
            add(r, -J, w, (neg[R], p0, A, p2, B), (neg[R], p3, A, p1, B), 0)
    e = entries("PairHop")
    if e is not None:
        r, a, b, J, sp, w = e
        R, A, B = r[:, None], a[:, None], b[:, None]
        p0, p1, p2, p3 = (sp[:, i][None, :] for i in range(4))
        # the conjugated read conj(G[r, v,b, u,a]) is the linear read
        # rho(b v, a u; +r), anchored r cells away from the output: e^{-iq.r}
        add(r, J, w, (neg[R], p0, A, p3, B), (R, p2, B, p1, A), -1)
        if fock:
            add(r, -J, w, (neg[R], p0, A, p2, B), (R, p3, B, p1, A), -1)
    if not parts:
        empty = np.zeros(0, dtype=np.int64)
        return (empty, empty, np.zeros(0, dtype=np.complex128), empty,
                np.zeros(0, dtype=np.int8))
    return tuple(np.concatenate([p[i] for p in parts]) for i in range(5))


def _allowed_slots(problem, facts1, flat):
    """Whether each slot-grid flat index is inside one solver block and on
    interacting orbitals."""
    dims = _grid_dims(problem)
    _r, s, a, t, b = np.unravel_index(flat, dims)
    block_of = np.full(problem.nd, -1, dtype=np.int64)
    for k, blk in enumerate(problem.blocks.block_info):
        block_of[np.asarray(blk, dtype=np.int64)] = k
    ic = np.zeros(dims[2], dtype=bool)
    ic[list(facts1.interacting_orbitals)] = True
    i = np.asarray(_nd_index(problem, a, s), dtype=np.int64)
    j = np.asarray(_nd_index(problem, b, t), dtype=np.int64)
    return (block_of[i] == block_of[j]) & ic[a] & ic[b]


def _kept_terms(problem, facts1):
    out, inp, coef, rid, kind = _vertex_terms(problem)
    keep = _allowed_slots(problem, facts1, out) & _allowed_slots(problem, facts1, inp)
    return out[keep], inp[keep], coef[keep], rid[keep], kind[keep]


def _labels_of(problem, flat):
    dims = _grid_dims(problem)
    r, s, a, t, b = np.unravel_index(flat, dims)
    rv = _r_vectors(problem.shape)[r]
    return np.column_stack([a, s, b, t, rv]).astype(np.int64)


def _flat_of(P, labels):
    """The slot-grid flat index of each label row (r taken modulo the mesh)."""
    labels = np.asarray(labels, dtype=np.int64)
    n = int(P.norb_phys)
    nvol = int(np.prod(P.shape))
    rr = tuple(labels[:, 4 + i] % int(P.shape[i]) for i in range(3))
    rid = np.ravel_multi_index(rr, P.shape)
    return np.ravel_multi_index((rid, labels[:, 1], labels[:, 0], labels[:, 3], labels[:, 2]),
                                (nvol, 2, n, 2, n))


def _su2_sectors(labels):
    """``(bases, sector_labels)`` of the charge and S^z sectors, or ``None``
    when an orbital slot carries only one spin-diagonal label."""
    import scipy.sparse
    diag = labels[:, 1] == labels[:, 3]
    pos = {}
    for k in np.nonzero(diag)[0]:
        key = (int(labels[k, 0]), int(labels[k, 2])) + tuple(int(v) for v in labels[k, 4:])
        pos.setdefault(key, {})[int(labels[k, 1])] = int(k)
    slots = sorted(pos, key=lambda key: min(pos[key].values()))
    if any(len(pos[key]) != 2 for key in slots):
        return None
    d = len(slots)
    rows = np.array([[pos[key][0], pos[key][1]] for key in slots], dtype=np.int64).reshape(d, 2)
    cols = np.repeat(np.arange(d), 2)
    h = 1.0 / math.sqrt(2.0)
    bases, slabels = [], []
    for mu, sign in ((0, 1.0), (3, -1.0)):
        data = np.tile(np.array([h, sign * h], dtype=np.complex128), d)
        bases.append(scipy.sparse.csr_array((data, (rows.reshape(-1), cols)),
                                            shape=(labels.shape[0], d)))
        sl = np.array([[key[0], mu, key[1], key[2], key[3], key[4]] for key in slots],
                      dtype=np.int64).reshape(d, 6)
        sl.flags.writeable = False
        slabels.append(sl)
    return tuple(bases), tuple(slabels)


def pair_space(facts, *, su2=None) -> PairSpace:
    """The pair space of the Hartree-Fock vertex of the model of ``facts``
    (:class:`PairSpace`), on its primitive cell.

    ``P`` is the support of the solver's Hartree-Fock map in the label
    coordinates of :class:`PairSpace`: every slot a Hartree, Fock or
    PairHop term of a nonzero table entry reads or writes, with the model's
    Fock flag (without Fock: the on-site orbital-diagonal Hartree slots
    and the PairHop slots), on the interacting orbitals, and inside one
    solver block (so the 2Sz block mask and the blocks that separate
    orbitals without a connecting Transfer remove what the solver
    removes). The tables are reverse-closed, so the slots read and the
    slots written are the same set.

    ``su2``: the :class:`~hwave.initial.symmetry.Su2Result` of the model
    (computed on the primitive cell when ``None``). When its status is
    ``"passed"`` the pair space carries the charge and ``S^z`` sectors. A
    given ``"passed"`` result is re-checked cheaply on this model (normal
    mode, free S_z, ``H0`` and the map with this Fock flag, 3 generic G x 3
    rotations); when it does not hold, the full check is run again.
    """
    from hwave.initial import symmetry
    facts1 = _primitive_facts(facts)
    p = facts1.problem
    out, inp, _coef, _rid, _kind = _kept_terms(p, facts1)
    flat = np.unique(np.concatenate([out, inp]))
    labels = _labels_of(p, flat)
    labels.flags.writeable = False
    if su2 is None:
        su2 = symmetry.su2_check(p)
    elif su2.status == "passed" and not _su2_holds(p):
        logger.warning("pair_space: the given spin-rotation result says passed, but it "
                       "does not hold for this model (its Sz mode, spin-orbital mode or "
                       "Fock setting); the check is run again")
        su2 = symmetry.su2_check(p)
    status = su2.status
    sectors, bases, slabels, note = (SECTOR_FULL,), (None,), (None,), None
    if status == "passed":
        split = _su2_sectors(labels)
        if split is None:
            note = ("the spin-rotation check passed, but an orbital slot of the pair "
                    "space has only one spin-diagonal label; one sector is kept")
            logger.warning(note)
        else:
            sectors, bases, slabels = (SECTOR_CHARGE, SECTOR_SZ), split[0], split[1]
    return PairSpace(labels=labels, sectors=sectors, sector_basis=bases,
                     sector_labels=slabels, shape=tuple(int(v) for v in p.shape),
                     norb_phys=int(p.tables.norb_phys), include_fock=bool(p.iflag_fock),
                     spin_orbital=bool(p.enable_spin_orbital), su2_status=status,
                     sector_note=note)


def _su2_holds(p):
    """A cheap re-check of a passed spin-rotation result on ``p``: normal
    mode with a free S_z, ``[S, H0] = 0`` and the Hartree-Fock map (with
    ``p``'s Fock flag and blocks) equivariant under 3 generic G x 3 random
    rotations."""
    from hwave.initial import symmetry
    if p.enable_spin_orbital or not p.sz_free:
        return False
    if symmetry._h0_commutator(p, (0, 1, 2)) > symmetry.SU2_TOL:
        return False
    return symmetry._random_residual(p, (0, 1, 2), 364) <= symmetry.SU2_TOL


def _check_labels(P, problem, facts1):
    """``ValueError`` unless ``P`` holds exactly the slots of this model's
    Hartree-Fock map (a pair space of another model on the same mesh)."""
    out, inp, _coef, _rid, _kind = _kept_terms(problem, facts1)
    want = np.unique(np.concatenate([out, inp]))
    have = np.sort(_flat_of(P, P.labels))
    if want.shape != have.shape or not np.array_equal(want, have):
        raise ValueError("vertex: the pair space does not hold exactly the slots of the "
                         "Hartree-Fock map of this model ({} labels, the model has {}); "
                         "build it with pair_space(facts) of this model".format(
                             have.size, want.size))


def _check_pair_space(P, problem):
    if not isinstance(P, PairSpace):
        raise ValueError("vertex: P must be a PairSpace, got {!r}".format(type(P).__name__))
    if (tuple(P.shape) != tuple(int(v) for v in problem.shape)
            or P.norb_phys != int(problem.tables.norb_phys)
            or P.include_fock != bool(problem.iflag_fock)
            or P.spin_orbital != bool(problem.enable_spin_orbital)):
        raise ValueError("vertex: the pair space belongs to another model (mesh {}, {} "
                         "physical orbitals); build it with pair_space(facts) of this "
                         "model".format(list(P.shape), P.norb_phys))


def _apply_vertex(problem, ph, x):
    """The vertex on a batch of slot-grid vectors ``x`` ``(c, nvol, 2, n, 2,
    n)`` (rho-components; ``ph`` the phases ``e^{iq.r}``): the solver's
    kernel written on the slot grid, with the conjugated PairHop read
    replaced by its linear read. No block projection (the caller restricts
    to the pair space)."""
    tables = problem.tables
    nvol, _, n = _grid_dims(problem)[:3]
    neg = _neg_r(problem.shape)
    fock = bool(problem.iflag_fock)
    out = np.zeros_like(x)
    diag = np.arange(n)
    for name in _HF_LOOP:
        tab = tables.inter_table.get(name)
        if tab is None:
            continue
        J = np.asarray(tab).reshape(nvol, n, n)
        spin = np.asarray(tables.spin_table[name])
        xb = x[:, 0][:, :, diag, :, diag]            # (n, c, v, u): rho(b v, b u; 0)
        hh = np.einsum("rab,r,suvt,bcvu->csta", J, ph, spin, xb, optimize=True)
        out[:, 0][:, :, diag, :, diag] += np.moveaxis(hh, 3, 0)
        if fock:
            f = np.einsum("rab,sutv,crvaub->crsatb", J, spin, x[:, neg], optimize=True)
            out -= f[:, neg]
    tab = tables.inter_table.get("PairHop")
    if tab is not None:
        J = np.asarray(tab).reshape(nvol, n, n)
        spin = np.asarray(tables.spin_table["PairHop"])
        phc = np.conj(ph)
        h = np.einsum("rab,r,suvt,crvbua->crsatb", J, phc, spin, x, optimize=True)
        if fock:
            h = h - np.einsum("rab,r,sutv,crvbua->crsatb", J, phc, spin, x, optimize=True)
        out += h[:, neg]
    return out


def vertex(facts, q, P, *, sector=None) -> np.ndarray:
    """The Hartree-Fock vertex ``Gamma(q) = dSigma'/drho`` on the pair space
    ``P`` (:func:`pair_space` of the same model), dense ``(|P|, |P|)``
    complex128 (the small tier).

    ``Gamma`` is the complex-linear extension of the solver's real-linear
    map (``hf_potential`` with the model's Fock flag and blocks) in the
    coordinates of :class:`PairSpace`: the PairHop read ``conj(G[r, v, b,
    u, a])`` is replaced by the linear read ``G[-r, u, a, v, b]`` (equal on
    Hermitian G). It does not depend on the para state: the map is linear.
    ``q`` (three rationals, reduced lattice units) must be on the k-mesh of
    the primitive cell (``ValueError`` otherwise). Facts of a folded cell
    are reloaded on the primitive cell at every call, so a caller with
    many q passes primitive-cell facts. ``sector``: one of ``P.sectors``; the block
    ``U^dag Gamma U`` of that sector's basis (``None``: the full space).
    """
    facts1 = _primitive_facts(facts)
    p = facts1.problem
    _check_pair_space(P, p)
    _check_labels(P, p, facts1)
    k = _sector_index(P, sector)
    ph = _phases(p.shape, q)
    dims = _grid_dims(p)
    gsize = int(np.prod(dims))
    flat = _flat_of(P, P.labels)
    npair = flat.size
    gamma = np.zeros((npair, npair), dtype=np.complex128)
    batch = max(1, _VERTEX_BATCH_ENTRIES // gsize)
    for c0 in range(0, npair, batch):
        cols = np.arange(c0, min(npair, c0 + batch))
        x = np.zeros((cols.size, gsize), dtype=np.complex128)
        x[np.arange(cols.size), flat[cols]] = 1.0
        y = _apply_vertex(p, ph, x.reshape((cols.size,) + dims)).reshape(cols.size, gsize)
        gamma[:, cols] = y[:, flat].T
    if k is None or P.sector_basis[k] is None:
        return gamma
    u = P.sector_basis[k].toarray()
    return u.conj().T @ gamma @ u


def vertex_sparse(facts, q, P, *, sector=None):
    """The vertex of :func:`vertex` as a ``scipy.sparse.csr_array``
    (complex128), assembled entry by entry from the interaction tables
    without a dense ``|P|^2`` array (the medium tier; its ``nnz`` feeds the
    cost model).

    The entries are added in the order of the solver's kernel (types,
    Hartree then Fock, table entries, spin entries), duplicates are summed
    (``sum_duplicates``) and the indices sorted; an entry whose terms
    cancel stays stored, so ``nnz`` does not depend on q. ``sector`` as in
    :func:`vertex` (``U^dag Gamma U`` in sparse arithmetic)."""
    import scipy.sparse
    facts1 = _primitive_facts(facts)
    p = facts1.problem
    _check_pair_space(P, p)
    _check_labels(P, p, facts1)
    k = _sector_index(P, sector)
    ph = _phases(p.shape, q)
    out, inp, coef, rid, kind = _kept_terms(p, facts1)
    flat = _flat_of(P, P.labels)
    order = np.argsort(flat, kind="stable")
    sflat = flat[order]
    row = np.searchsorted(sflat, out)
    col = np.searchsorted(sflat, inp)
    npair = flat.size
    ok = (row < npair) & (col < npair)
    ok[ok] &= (sflat[row[ok]] == out[ok]) & (sflat[col[ok]] == inp[ok])
    if not np.all(ok):
        raise ValueError("vertex_sparse: the pair space does not hold every slot of the "
                         "Hartree-Fock map of this model; build it with pair_space(facts)")
    phase = np.where(kind == 1, ph[rid], np.where(kind == -1, np.conj(ph[rid]), 1.0))
    data = coef * phase
    g = scipy.sparse.coo_array((data, (order[row], order[col])), shape=(npair, npair))
    g = scipy.sparse.csr_array(g)
    g.sum_duplicates()
    g.sort_indices()
    if k is None or P.sector_basis[k] is None:
        return g
    u = P.sector_basis[k]
    g = scipy.sparse.csr_array(u.conj().T @ g @ u)
    g.sum_duplicates()
    g.sort_indices()
    return g


def vertex_self_adjoint(Gamma) -> Tuple[bool, float]:
    """``(ok, residual)``: ``residual = ||Gamma - Gamma^dag||_F / max(1,
    ||Gamma||_F)`` and ``ok`` when it is at most
    :data:`VERTEX_SELF_ADJOINT_TOL` (dense or sparse ``Gamma``). A failure
    is a convention or implementation defect, not a property of a model."""
    import scipy.sparse
    if scipy.sparse.issparse(Gamma):
        g = scipy.sparse.csr_array(Gamma)
        diff = g - g.conj().T
        dn = float(np.sqrt(np.sum(np.abs(diff.data) ** 2))) if diff.nnz else 0.0
        gn = float(np.sqrt(np.sum(np.abs(g.data) ** 2))) if g.nnz else 0.0
    else:
        g = np.asarray(Gamma)
        dn = float(np.linalg.norm(g - g.conj().T))
        gn = float(np.linalg.norm(g))
    res = dn / max(1.0, gn)
    return res <= VERTEX_SELF_ADJOINT_TOL, res


def _check_bands(bands, P, who):
    """``ValueError`` unless ``bands`` live on the mesh and index space of
    the pair space ``P`` (the primitive cell, its mode). The block structure
    is checked by the facade: a label of ``P`` joining two blocks of the
    bands (bands of a fixed 2Sz against a pair space without it) is
    refused there."""
    if not isinstance(bands, initial_facade.Bands):
        raise ValueError("{}: bands must be a Bands (initial_facade.bands_of), got {!r}"
                         .format(who, type(bands).__name__))
    so = int(bands.ns) == 1
    nphys = int(bands.norb) // 2 if so else int(bands.norb)
    if (tuple(int(n) for n in bands.kgrid) != tuple(P.shape) or so != P.spin_orbital
            or nphys != P.norb_phys):
        raise ValueError("{}: the bands (mesh {}, {} physical orbitals, spin-orbital {}) "
                         "do not belong to the pair space (mesh {}, {} physical orbitals, "
                         "spin-orbital {}); use the bands of the primitive cell of the same "
                         "model".format(who, list(bands.kgrid), nphys, so, list(P.shape),
                                        P.norb_phys, P.spin_orbital))


def lindhard(bands, q, P, *, sector=None) -> np.ndarray:
    """The static Lindhard response ``chi0(q)`` on the pair space ``P``
    (:func:`initial_facade.static_lindhard` on ``P.labels``), dense,
    Hermitian and positive semidefinite, in the coordinates of
    :class:`PairSpace` (rho-components in and out: ``delta x = chi0
    (-h)``). ``bands``: the para bands of the primitive cell
    (:attr:`ParaResult.bands`). ``sector``: one of ``P.sectors``; the block
    ``U^dag chi0 U`` of that sector's basis (``None``: the full space)."""
    _check_bands(bands, P, "lindhard")
    k = _sector_index(P, sector)
    chi = initial_facade.static_lindhard(bands, q, P.labels)
    if k is None or P.sector_basis[k] is None:
        return chi
    u = P.sector_basis[k].toarray()
    return u.conj().T @ chi @ u


def lindhard_apply(bands, q, P, v, *, sector=None) -> np.ndarray:
    """``chi0(q) v`` of :func:`lindhard` without forming ``chi0``
    (:func:`initial_facade.apply_lindhard`); ``v`` ``(d,)`` or ``(d, c)``,
    ``d`` the dimension of ``sector`` (``|P|`` for ``None``)."""
    _check_bands(bands, P, "lindhard_apply")
    k = _sector_index(P, sector)
    if k is None or P.sector_basis[k] is None:
        return initial_facade.apply_lindhard(bands, q, P.labels, v)
    u = P.sector_basis[k]
    v = np.asarray(v)
    if v.ndim not in (1, 2) or v.shape[0] != u.shape[1]:
        raise ValueError("lindhard_apply: v must be (d,) or (d, c) with d = {} (sector "
                         "{}), got {}".format(u.shape[1], sector, v.shape))
    y = initial_facade.apply_lindhard(bands, q, P.labels, u @ v.astype(np.complex128))
    return u.conj().T @ y
