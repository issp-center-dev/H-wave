"""Instability hints: the hint temperature ``T_hint`` (design section 6.4).

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
"""
import math
import numbers
from dataclasses import dataclass
from typing import Literal, Optional, Tuple

import numpy as np

from hwave.solver import initial_facade

__all__ = ["THint", "t_hint_bare", "W_FLOOR_ACT", "T_HINT_REL_FLOOR",
           "STRUCTURAL_AXIS_MAX", "WIDE_WINDOW_FACTOR", "group_bands"]

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
    """The hint temperature and its inputs (design section 6.4).

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
    bands and k mesh (design section 6.4).

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
