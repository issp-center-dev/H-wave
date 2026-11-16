"""The spin-orbital bond-resolved FLEX (issue #270): one map of the FLEX
self-consistency on the bond channels (:func:`so_bond_flex_map`) and the
admission rows of the FLEX SCF beyond the bond plan's own
(:func:`flex_so_bond_rows`, for :func:`hwave.solver.so_bond.plan_so_bond`
with the consumer ``flex_scf``).
"""
import contextlib
import math
from typing import NamedTuple, Optional

import numpy as np

from . import flex_bond, so_bond

__all__ = ["SOBondFlexMap", "so_bond_flex_map", "flex_so_bond_rows", "SOBondFlexAdmission",
           "admit_so_bond_flex", "CAP_KEY", "MAX_SHELLS_KEY", "FREQ_BATCH_CANDIDATES"]

_ENTRY = 16                     # bytes per complex128 entry
#: the per-map Green-function digest hashes at most this many bytes per copy
_DIGEST_COPY_CAP = 64 * 2 ** 20


def flex_so_bond_rows(*, nmat, nvol, nd, depth, mixing, n_types, split_seed=False):
    """The FLEX SCF rows of a spin-orbital bond plan: ``(resident_rows,
    extra_rows)`` for :func:`hwave.solver.so_bond.plan_so_bond`
    (``resident_rows=`` and ``extra_rows=``).

    With ``G = 16 nmat nvol nd^2`` (one spin-orbital Green function) and
    ``C = 16 nmat nvol nd^4`` (one channel-0 susceptibility block over the
    frequency grid):

    * resident: the persistent rows of :func:`hwave.solver.flex_bond.
      flex_scf_rows` with ``nmodes = nd`` (``state``, ``seed_envelope``,
      ``anderson_history``, ``anderson_work``, ``eigenpairs_hf``,
      ``flex_arrays``, ``hf_tables``) and ``collapses = 2 C`` (the chi0q
      and chiq buffers of the last map);
    * extra: its phase rows (``green_mu``, ``density_hf``, ``convergence``,
      ``mixing``; each row named by its phase), ``("convergence",
      "green_digest", 2 min(G, 64 MiB))`` (the per-map Green-function
      digest), ``("transport", "transport", 4 G + 4 C)`` and
      ``("post_scf", "post_scf", G + C)``.

    ``split_seed`` is False on the spin-orbital path (its seeds are
    refused), so ``seed_envelope`` is 0 by default. Every byte count is an
    ``int``."""
    nmat, nvol, nd = int(nmat), int(nvol), int(nd)
    persistent, phase = flex_bond.flex_scf_rows(nmat=nmat, nvol=nvol, nmodes=nd, depth=depth,
                                                mixing=mixing, n_types=n_types,
                                                split_seed=split_seed)
    G = _ENTRY * nmat * nvol * nd * nd
    C = _ENTRY * nmat * nvol * nd ** 4
    resident = tuple((k, int(v)) for k, v in persistent.items()) + (("collapses", int(2 * C)),)
    extra = tuple((k, k, int(v)) for k, v in phase.items()) + (
        ("convergence", "green_digest", int(2 * min(G, _DIGEST_COPY_CAP))),
        ("transport", "transport", int(4 * G + 4 * C)),
        ("post_scf", "post_scf", int(G + C)))
    return resident, extra


# ---------------------------------------------------------------------------
# One map of the spin-orbital bond FLEX
# ---------------------------------------------------------------------------

class SOBondFlexMap(NamedTuple):
    """The result of :func:`so_bond_flex_map`.

    Fields
    ------
    sigma_fluct : ndarray (1, nmat, nvol, nd, nd), complex128
        The fluctuation self-energy of the map.
    chi0_c0, chi_c0 : ndarray (nmat, nvol, nd**2, nd**2), complex128
        The channel-0 blocks of the bare bubble ``chibar_B`` and of its
        dressing ``chi_B`` on every bosonic slice (the caller's ``chi0_out``
        / ``chi_out`` when given).
    cond_min : float
        The smallest conditioning score of the dressing guard over the map;
        ``nan`` when the guard is off (``cond_tol = None``).
    sigma_min_sector, sigma_max_sector : str or None
        S_z-packed layout: the sectors of the smallest and largest singular
        value at the worst ``(l, q)`` of the map; ``None`` dense or without
        the guard.
    pole : tuple
        :func:`hwave.solver.so_bond.static_pole_indicator` of the ``nu =
        0`` slice of ``chibar_B``: ``(entry, per_sector, worst)``.
    """

    sigma_fluct: np.ndarray
    chi0_c0: np.ndarray
    chi_c0: np.ndarray
    cond_min: float
    sigma_min_sector: Optional[str]
    sigma_max_sector: Optional[str]
    pole: tuple


def _null_trace(name):
    return contextlib.nullcontext()


def _check_buffer(name, buf, shape):
    if buf is None:
        return np.empty(shape, dtype=np.complex128)
    if not (isinstance(buf, np.ndarray) and buf.dtype == np.complex128 and buf.shape == shape):
        raise ValueError("so_bond_flex_map: {} must be a complex128 array of shape {}, got {}"
                         .format(name, shape, so_bond._describe(buf)))
    return buf


def so_bond_flex_map(store, gamma_b, green_kw, green_scf, green0_tail, beta, geom, *,
                     plan, cond_tol, workers=1, trace=None, chi0_out=None, chi_out=None):
    """One map of the spin-orbital bond FLEX: the fluctuation self-energy of
    the Green function ``green_kw``.

    * bubble (phase ``bubble``): ``chibar_B`` of ``green_scf`` with its tail
      ``green0_tail`` on every bosonic slice into the store slot ``chibar``
      (:func:`hwave.solver.so_bond.so_bond_bubble`);
    * per frequency batch ``[l0, l1)`` of ``plan.nb`` slices (phase
      ``dressing``): the channel-0 block of ``chibar_B`` into ``chi0_c0``,
      the static pole indicator of the ``nu = 0`` slice (``l = nmat // 2``),
      the dressing ``chi_B = [1 + chibar_B Gamma_B]^{-1} chibar_B``
      (:func:`hwave.solver.so_bond.so_bond_dress` with ``cond_tol``) and its
      channel-0 block into ``chi_c0``; then (phase ``w_conversion``) ``W =
      P_B[Gamma_B (chi_B - chibar_B / 2) Gamma_B]``
      (:func:`hwave.solver.so_bond.so_bond_w_batch`) into the store slot
      ``W``, which may be an alias of ``chibar`` (the batch ``W`` is complete
      before it overwrites the batch of ``chibar`` it was built from);
    * transport (phase ``transport``):
      ``flex_bond.calc_self_energy_bond(store, green_kw, beta, geom,
      geom.spatial_shape, nd, workers)``.

    Parameters
    ----------
    store : the block store of ``plan`` (:func:`hwave.solver.so_bond.open_store`)
    gamma_b : the vertex ``Gamma_B`` in the plan's layout
        (:func:`hwave.solver.so_bond.so_bond_vertex`).
    green_kw, green_scf : ndarray (1, nmat, nvol, nd, nd)
        The Green function of the transport and the (tail-subtracted) one of
        the bubble; ``green0_tail`` the tail of ``green_scf`` (same shape) or
        ``None``.
    beta : float
    geom : hwave.solver.so_bond.BondSlotGeometry
        The slot geometry of the bond channels.
    plan : SOBondPlan
        A plan with the consumer ``flex_scf``.
    cond_tol : float or None
        The conditioning guard of the dressing (``None`` disables it).
    workers : int
        FFT workers of the bubble's and the transport's FFTs.
    trace : callable, optional
        ``trace(phase)`` returns a context manager entered around each phase
        (e.g. a tracemalloc scope); ``None`` is a null context.
    chi0_out, chi_out : ndarray, optional
        Caller-owned complex128 ``(nmat, nvol, nd**2, nd**2)`` buffers,
        overwritten in place and returned as ``chi0_c0`` / ``chi_c0``;
        ``None`` allocates them.

    Every module function is looked up at call time. Returns
    :class:`SOBondFlexMap`."""
    if not isinstance(plan, so_bond.SOBondPlan):
        raise ValueError("so_bond_flex_map: plan must be a SOBondPlan (plan_so_bond), got {}"
                         .format(type(plan).__name__))
    if "flex_scf" not in plan.consumers:
        raise ValueError("so_bond_flex_map: the plan needs the consumer flex_scf (closed "
                         "consumers: {})".format(", ".join(plan.consumers)))
    nmat, nvol, nd, nb = plan.nmat, plan.nvol, plan.nd, plan.nb
    gshape = (1, nmat, nvol, nd, nd)
    for name, g in (("green_kw", green_kw), ("green_scf", green_scf),
                    ("green0_tail", green0_tail)):
        if g is None and name == "green0_tail":
            continue
        if tuple(np.shape(g)) != gshape:
            raise ValueError("so_bond_flex_map: {} must have shape {}, got {}".format(
                name, gshape, tuple(np.shape(g))))
    P = nd * nd
    cshape = (nmat, nvol, P, P)
    chi0_c0 = _check_buffer("chi0_out", chi0_out, cshape)
    chi_c0 = _check_buffer("chi_out", chi_out, cshape)
    if np.shares_memory(chi0_c0, chi_c0):
        raise ValueError("so_bond_flex_map: chi0_out and chi_out share memory")
    traced = _null_trace if trace is None else trace
    shape = tuple(geom.spatial_shape)
    dense = plan.layout == "dense"

    with traced("bubble"):
        so_bond.so_bond_bubble(green_scf, green0_tail, beta, geom, shape, plan=plan, freq="all",
                               sink=store)

    cond_min, s_min, s_max = math.inf, None, None
    pole = None
    l_static = nmat // 2
    for l0 in range(0, nmat, nb):
        l1 = min(l0 + nb, nmat)
        diag = {}
        with traced("dressing"):
            X = store.get_freq_batch("chibar", l0, l1)
            chi0_c0[l0:l1] = so_bond.channel0_block(X, nd=nd)
            if l0 <= l_static < l1:
                i = l_static - l0
                static = X[i] if dense else so_bond.PackedArrays(
                    X.sector_map, tuple(b[i] for b in X.blocks))
                pole = so_bond.static_pole_indicator(static, gamma_b, shape)
                del static
            Y = so_bond.so_bond_dress(X, gamma_b, plan=plan, l0=l0, spatial_shape=shape,
                                      cond_tol=cond_tol, diagnostics=diag)
            chi_c0[l0:l1] = so_bond.channel0_block(Y, nd=nd)
        with traced("w_conversion"):
            Wb = so_bond.so_bond_w_batch(X, Y, gamma_b, plan=plan)
            del Y
            store.put_freq_batch("W", l0, l1, Wb)
            del Wb, X
        c = diag.get("cond_min")
        if c is not None and c < cond_min:
            cond_min = float(c)
            s_min, s_max = diag.get("sigma_min_sector"), diag.get("sigma_max_sector")
    if cond_tol is None or not math.isfinite(cond_min):
        cond_min, s_min, s_max = float("nan"), None, None

    with traced("transport"):
        sigma = flex_bond.calc_self_energy_bond(store, green_kw, beta, geom, shape, nd, workers)
    return SOBondFlexMap(sigma, chi0_c0, chi_c0, cond_min, s_min, s_max, pole)


# ---------------------------------------------------------------------------
# The admission of the spin-orbital bond FLEX
# ---------------------------------------------------------------------------

#: the settings named by the admission's texts
CAP_KEY = "[mode.param] longitudinal_bond_memory_cap_gb"
MAX_SHELLS_KEY = "[mode.param] longitudinal_bond_max_shells"
#: the frequency batches tried, largest first, when longitudinal_bond_freq_batch is unset
FREQ_BATCH_CANDIDATES = (16, 8, 4, 2, 1)
_B1_INFO = ("spin-orbital bond FLEX: no off-site row; the bond machinery runs with the "
            "on-site channel only (B = 1).")
_NB_INFO = ("spin-orbital bond FLEX: frequency batch nb = {nb} (the largest of 16, 8, 4, 2, "
            "1 not above Nmat that fits the cap)")
_REMEDY_AUTO = ("The frequency batch was selected automatically; the smallest batch (nb = 1) "
                "was tried.")
_REMEDY_SET = ("longitudinal_bond_freq_batch = {nb} is set; a smaller value, or removing it so "
               "that the run chooses, lowers the dressing rows.")
_REMEDY_SET_1 = "longitudinal_bond_freq_batch = 1 is set (the smallest batch)."


@contextlib.contextmanager
def _quiet_plan():
    """The candidate plans of the batch selection are not logged (their
    reason line would name an unlimited cap); the admitting plan is."""
    drop = lambda record: False     # noqa: E731
    so_bond.logger.addFilter(drop)
    try:
        yield
    finally:
        so_bond.logger.removeFilter(drop)


class SOBondFlexAdmission(NamedTuple):
    """The admission of a spin-orbital bond FLEX run
    (:func:`admit_so_bond_flex`), decided before any bond array exists.

    Fields
    ------
    topo, geom
        The bond topology (:func:`hwave.solver.so_bond.so_bond_topology`)
        and its slot geometry.
    layout, reason : str
        ``"dense"`` or ``"sz_packed"`` and the selection's reason line.
    sector_map
        The S_z sector map (packed layout; ``None`` dense).
    plan : SOBondPlan
        The admitting plan (consumer ``flex_scf``), the one memory authority
        of the run.
    caps : SOBondCaps
        The resolved host cap.
    bond_inputs : BondInputs
        The bond provenance of the topology.
    nb : int
        The frequency batch of the plan.
    sector_dyson : bool
        Whether the Dyson step inverts per S_z sector (the band and every
        interaction record conserve S_z)."""
    topo: object
    geom: object
    layout: str
    reason: str
    sector_map: object
    plan: object
    caps: object
    bond_inputs: object
    nb: int
    sector_dyson: bool


def admit_so_bond_flex(*, compiled, transfer, extern, cell, shape, nmat, max_shells,
                       layout_request, cap_raw, freq_batch, depth, mixing, n_types, log):
    """The :class:`SOBondFlexAdmission` of a spin-orbital bond FLEX run.

    ``compiled`` is the run's compiled interaction, ``transfer`` / ``extern``
    its pre-fold Transfer and Extern tables, ``cell`` the lattice vectors
    (rows; the Geometry ``rvec``), ``shape`` the k mesh, ``nmat`` the
    Matsubara count, ``max_shells`` the shell cutoff (``None``: every
    shell), ``layout_request`` the ``so_bond_layout`` value (``None``:
    ``"auto"``), ``cap_raw`` the raw memory cap in GiB (``None``: ``0.8 x``
    available, :func:`hwave.solver.so_bond.resolve_host_cap`),
    ``freq_batch`` the set frequency batch (``None``: the largest of
    :data:`FREQ_BATCH_CANDIDATES` not above ``nmat`` whose peak fits the
    cap), ``depth`` / ``mixing`` / ``n_types`` the SCF rows' parameters
    (:func:`flex_so_bond_rows`).

    Logs (on ``log``) the B = 1 notice, the layout's reason line, the cap
    and the auto-selected batch. An over-cap run is refused with
    :func:`hwave.solver.so_bond.plan_so_bond`'s text, which names how the
    batch was chosen; no other exception of a candidate plan is caught."""
    nd = 2 * int(compiled.norb_phys)
    topo, bond_inputs, _ = so_bond.bond_run(compiled, cell, max_shells)
    shape = tuple(int(x) for x in shape)
    geom = so_bond.slot_geometry(topo, shape)
    B = int(geom.B)
    if B == 1:
        log.info(_B1_INFO)
    layout, reason = so_bond.select_layout(
        compiled, transfer, override="auto" if layout_request is None else layout_request)
    log.info(reason)
    sector_map = so_bond.sz_sector_map(geom, nd) if layout == "sz_packed" else None
    caps = so_bond.resolve_host_cap(cap_raw, key=CAP_KEY, log=log)
    nmat = int(nmat)
    nvol = int(np.prod(shape))
    resident, extra = flex_so_bond_rows(nmat=nmat, nvol=nvol, nd=nd, depth=depth,
                                        mixing=mixing, n_types=n_types)
    common = dict(nmat=nmat, nvol=nvol, B=B, nd=nd, consumers=("flex_scf",), layout=layout,
                  layout_reason=reason, cap_key=CAP_KEY,
                  shell_rows=so_bond.shell_rows(compiled, topo, cell),
                  max_shells_key=MAX_SHELLS_KEY, envelope=so_bond.envelope_text(layout),
                  resident_rows=resident, extra_rows=extra)
    if freq_batch is None:
        candidates = [n for n in FREQ_BATCH_CANDIDATES if n <= nmat]
        nb = candidates[-1]
        unlimited = so_bond.SOBondCaps(None, None, None)
        for n in candidates:
            with _quiet_plan():
                trial = so_bond.plan_so_bond(caps=unlimited, nb=n, **common)
            if caps.host_bytes is None or trial.peak_host <= caps.host_bytes:
                nb = n
                break
        note = _REMEDY_AUTO
    else:
        nb = int(freq_batch)
        note = _REMEDY_SET_1 if nb == 1 else _REMEDY_SET.format(nb=nb)
    plan = so_bond.plan_so_bond(caps=caps, nb=nb, remedy_note=note, **common)
    if freq_batch is None:
        log.info(_NB_INFO.format(nb=plan.nb))
    sector_dyson = bool(so_bond.band_conserves_sz(transfer, extern, compiled.norb_phys).conserving
                        and so_bond.records_conserve_sz(compiled) is None)
    return SOBondFlexAdmission(topo, geom, layout, reason, sector_map, plan, caps, bond_inputs,
                               int(plan.nb), sector_dyson)
