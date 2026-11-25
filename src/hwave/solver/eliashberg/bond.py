"""The static bond-resolved pairing path of hwave_sc: topology, vertices, the bond kernel operator, the bond Green function and the diagnostics record.

Moved verbatim from hwave.sc. hwave.sc re-exports only its public API
(calc_eliashberg, main, reject_spin_orbital_mode and the public
constants); import the other names here from this module.
"""
from __future__ import annotations

import logging
from typing import NamedTuple
import numpy as np
from scipy.sparse.linalg import LinearOperator

from hwave.solver.rpa import TAIL_ENDPOINT_CONVENTION
from hwave.solver import bond_channels
from hwave.solver import bubble
from hwave.solver import green as green_mod
from hwave.solver.so_pairing import _SO_WHAT
from hwave.solver.eliashberg.config import _validate_bond_green_nmat, _SO_BOND_CAP_KEY
from hwave.solver.eliashberg.band import _determine_mu
from hwave.solver.eliashberg.kernel import _green_sc_to_canonical, _green_canonical_to_sc, _BondGreen
from hwave.solver.eliashberg.eigen import _odd_harmonic_basis, _project_gap_parity, _SOAntisymmetricBasis
from hwave.solver.eliashberg.io import _refuse_so_green_file, _load_green_npz, _so_green_fp, _so_borrowed_green, _so_bond_green_member_bytes

logger = logging.getLogger("hwave_sc")


# Degeneracy tolerance of the opt-in bond diagnostics' eigenvalue clustering
# ([eliashberg] bond_diagnostics; spec S7.7 default).
_BOND_DIAGNOSTICS_DEG_TOL = 1.0e-3


# Number of (N_q, ND, ND)-shaped complex128 arrays that the bond-channel
# resource preflight (below) assumes are alive simultaneously. The high-water
# mark of ``_build_bond_operator`` is the moment a dressed channel is added
# to the pairing vertex: chi_bar, S, C, the pairing vertex Gamma, the dressed
# channel chi, the sandwich product V chi V and its coefficient-scaled copy
# (``PairVertexAccumulator``), plus one transient (the dressing denominator,
# the matrix-product intermediate V chi, or the kernel's transformed blocks).
# Measured at 7.1 arrays (``tests/test_sc_bond.py::TestPreflightOperatorBudget``).
# Named so the preflight's byte estimate and this list can never silently
# drift apart.
_BOND_N_Q_ARRAYS = 8


def _bond_memory_estimate(norb, bond_set, Nx, Ny, Nz, nmat, tail_on):
    """Byte budget for the bond path, broken down by buffer family.

    Split out of :func:`_bond_resource_preflight` so the accounting can be
    asserted directly against a MEASURED peak (``tests/test_sc_bond.py``):
    a preflight that undercounts is worse than no preflight, because a run it
    waves through then OOMs anyway.

    Everything is complex128 (16 B) end to end on the bond path (the
    unified-bubble-kernel spec's "Memory contract": ``build_green`` produces
    complex128, the external-npz loader promotes to complex128 at load, so
    there is one itemsize and no input-vs-compute dtype split). ``unit =
    N_q * nmat * 16`` is one ``norb``-pair-free frequency-grid buffer, so
    ``S_in = norb**2 * unit`` is one FULL-size Green-function-shaped array
    (``(norb, norb, Nx, Ny, Nz, nmat)`` worth of complex128).

    ``q_bytes``
        ``_BOND_N_Q_ARRAYS`` q-resolved ``(N_q, ND, ND)`` arrays alive
        simultaneously in the operator assembly (see the list beside that
        constant). ``bond_bubble_static``'s output ``chi_bar`` is the first of
        these, so it is NOT counted again in ``bubble_bytes``.
    ``bubble_bytes``
        ``bond_bubble``'s simultaneous working set at its high-water mark, the
        Green carrier excluded (that is ``carrier_bytes``):
        ``BOND_BUBBLE_N4_BUFFERS`` ``norb**4``-sized plus
        ``BOND_BUBBLE_N2_BUFFERS`` ``norb**2``-sized frequency-grid buffers.
        The buffer-by-buffer list, and the ``del``/in-place discipline in
        ``bond_bubble`` that holds the real peak down to it, are documented
        beside those two constants in ``hwave.solver.bond_channels``; they are
        imported from there so the two sides cannot silently desync.
    ``carrier_bytes``
        The Green carrier's (:class:`_BondGreen`) contribution, replacing
        the pre-tail-correction ``green_bytes`` (one bare ``green_kw``
        array). Derived from TWO phases the carrier passes through (spec
        "Memory contract"), each expressed in units of ``S_in``:

        * carrier CONSTRUCTION (:func:`_build_bond_green`, transient):
          ``4 * S_in`` when ``tail_on`` (``build_green``'s canonical
          ``full_kw`` + canonical ``deflated_kw`` + ``green0_tail`` all
          alive at once, PLUS ``full_sc`` during the layout conversion,
          before the canonical ``full_kw`` is released) and ``2 * S_in``
          tail-off (``full``/``deflated`` alias one array; the sc-layout
          conversion adds the second). The external-npz branch peaks at
          ``2 * S_in`` (source + complex128-promoted copy) before the
          source is released -- covered by the same tail-off bound.

          Pass B review, round 1, measured ~5-7x ``S_in`` here on a SMALL
          fixture and (wrongly) concluded the ``4x`` bound undercounted;
          round 2 found the true cause and fixed it at the source instead
          of inflating this estimate: ``build_green``'s locals ``Vg``
          (full-size) and ``g_deflated`` were never released and stayed
          alive through the ``green0_tail``/``full_kw`` assembly (see
          ``green.build_green``'s ``del Vg, g_deflated`` -- load-bearing
          for this bound), and round 1's tiny fixture additionally paid a
          FIXED small-buffer overhead (``V_conj_t``/``VVt``/frequency-grid
          arrays, not ``S_in``-scaled) that is disproportionate relative to
          a small ``S_in`` but negligible at realistic problem sizes. With
          both effects accounted for, the measured peak is ``4 * S_in``
          plus a small FIXED residual (measured 888-952 bytes across
          several fixture sizes -- NOT ``S_in``-scaled, so it vanishes in
          relative terms as ``S_in`` grows), pinned by
          ``tests.test_sc_bond.TestPreflightMemoryBudget.
          test_bond_green_construction_peak_within_carrier_budget`` (which
          carries a small, honestly-documented absolute slack for exactly
          that residual in its own comparison -- not in this estimate).
        * carrier SETTLED / bubble phase: ``n_green * S_in`` with
          ``n_green = 3`` (``full_sc`` + canonical ``deflated_kw`` +
          ``green0_tail``) when ``tail_on`` else ``2``, PLUS one extra
          ``S_in``-sized endpoint-branch buffer while ``tail_on`` (the
          bubble's per-pair tail-endpoint correction) -- ``4 * S_in``
          tail-on, ``2 * S_in`` tail-off, both equal to their respective
          construction-phase bound.

        Both phases top out at ``4 * S_in`` (tail-on) / ``2 * S_in``
        (tail-off), so ONE closed-form bound covers the whole carrier
        lifetime: ``carrier_bytes = (4 if tail_on else 2) * norb**2 *
        unit``.
    ``chi_bar_bytes``
        one ``(N_q, ND, ND)`` array; reported separately (it is already inside
        ``q_bytes``) because it is the one q-resolved array ``bond_bubble``
        itself allocates, so ``bubble_bytes + chi_bar_bytes`` is the budget a
        measurement of ``bond_bubble`` alone must fit into.
    ``peak``
        ``q_bytes + bubble_bytes + carrier_bytes`` -- the SUM combine rule;
        what the cap is applied to.

    Parameters
    ----------
    tail_on : bool
        ``coeff_tail_applied != 0`` (NEVER the requested value -- an
        external Green function always resolves ``tail_on = False``
        regardless of what was requested, since ``coeff_tail_applied`` is
        forced to ``0.0`` on that branch).
    """
    nd = norb * norb
    B = int(bond_set.n_channels)
    ND = nd * B
    Nq = int(Nx) * int(Ny) * int(Nz)
    itemsize = 16  # complex128
    unit = Nq * int(nmat) * itemsize

    chi_bar_bytes = Nq * ND * ND * itemsize
    q_bytes = _BOND_N_Q_ARRAYS * chi_bar_bytes
    bubble_bytes = (bond_channels.BOND_BUBBLE_N4_BUFFERS * norb ** 4
                    + bond_channels.BOND_BUBBLE_N2_BUFFERS * norb ** 2) * unit
    carrier_bytes = (4 if tail_on else 2) * (norb ** 2) * unit
    return {"nd": nd, "B": B, "ND": ND, "Nq": Nq,
            "chi_bar_bytes": chi_bar_bytes,
            "q_bytes": q_bytes,
            "bubble_bytes": bubble_bytes,
            "carrier_bytes": carrier_bytes,
            "peak": q_bytes + bubble_bytes + carrier_bytes}


def _bond_resource_preflight(norb, bond_set, Nx, Ny, Nz, nmat, cap_gb, *,
                             tail_on):
    """Estimate the peak memory of the bond path and refuse to exceed the cap.

    Spec S3.2 "Resource guard": only here are ``nd``, ``N_q``, ``B`` and the
    dtype all known, so this is where the ``ND = nd*B`` blow-up is caught --
    never a silent runaway allocation. Called BEFORE the Green carrier is
    built (review fix I-2, preserved through the unified-bubble-kernel
    switch: previously the bond-channel dispatch branch built the Green
    function first, so a large-Nmat run could OOM before the cap was
    ever consulted); the estimate includes the carrier's own allocation so
    the preflight is not blind to it.

    ``tail_on`` : bool
        Whether the bond path's tail correction will be applied
        (``coeff_tail_applied != 0``); known at this call site even though
        the carrier itself is built afterwards, because an external
        ``bond_green`` always resolves ``coeff_tail_applied = 0`` (see
        :func:`_build_bond_green`) regardless of what was requested.

    See :func:`_bond_memory_estimate` for the buffer-by-buffer accounting.
    """
    est = _bond_memory_estimate(norb, bond_set, Nx, Ny, Nz, nmat, tail_on)
    B, ND, Nq = est["B"], est["ND"], est["Nq"]
    q_bytes = est["q_bytes"]
    bubble_bytes = est["bubble_bytes"]
    carrier_bytes = est["carrier_bytes"]
    peak = est["peak"]

    logger.info(
        "Bond-channel preflight: B = %d channels, ND = nd*B = %d, "
        "N_q = %d, tail_on = %s, estimated peak memory = %.3f GB "
        "(cap %.3f GB)", B, ND, Nq, tail_on, peak / 1.0e9, cap_gb)

    if peak > cap_gb * 1.0e9:
        raise ValueError(
            "[eliashberg] bond_channels: estimated peak memory {:.3f} GB "
            "exceeds bond_memory_cap_gb = {:.3f} GB. The cost is driven by "
            "B = {} bond channels (Delta r = {}) giving ND = nd*B = {} on "
            "N_q = {} q-points ({:.3f} GB of q-resolved ND x ND arrays) plus "
            "{:.3f} GB of bond-bubble work buffers plus {:.3f} GB for the "
            "Green carrier (tail_on={}). Use fewer bond channels (declare "
            "fewer off-site shells, or cut declared-zero shells with "
            "bond_max_shells), reduce the k-grid/Nmat, or raise "
            "bond_memory_cap_gb.".format(
                peak / 1.0e9, cap_gb, B, list(bond_set.delta_r), ND, Nq,
                q_bytes / 1.0e9, bubble_bytes / 1.0e9,
                carrier_bytes / 1.0e9, tail_on))
    return peak


def _bond_locality_split(interactions):
    """The pre-fold locality split of the static bond path's interactions.

    ``hwave_sc`` refuses sublattice folding (``SubShape`` must be
    ``[1, 1, 1]``), so the table read from the files is already the pre-fold
    one and :func:`hwave.solver.offsite.split_locality` runs on it directly.
    No ``orbital_site`` labels are passed: the bond channels treat the
    exchange of every bond exactly, so -- as on the FLEX bond gate -- only
    ``R = 0`` rows are on-site.
    """
    from types import SimpleNamespace
    from hwave.solver.offsite import split_locality
    return split_locality(SimpleNamespace(param_ham=interactions),
                          SimpleNamespace(subshape=(1, 1, 1)))


def _resolve_bond_topology(interactions, geom_info, norb, bond_max_shells):
    """Resolve the bond topology of the static bond path.

    The same resolver as the FLEX bond gate and the dynamic bond pairing
    (``bond_channels.resolve_bond_topology``) over every off-site type the
    bond channels carry (``bond_channels.longitudinal_bond_types``:
    CoulombInter, Hund and Ising, plus Exchange and PairHop when declared
    off-site). Channel 0 is the on-site ``R = 0`` point; the off-site shells
    follow by bond length (measured with the lattice vectors of the
    geometry file), then lexicographically. ``bond_max_shells`` is the
    resolver's ``max_shells``: it keeps that many off-site shells and
    refuses a cutoff that would drop a declared nonzero coefficient.

    Returns
    -------
    (topo, view, split)
        The :class:`~hwave.solver.bond_channels.BondTopology`, its
        :class:`~hwave.solver.bond_channels.BondSetView` and the locality
        split (:func:`_bond_locality_split`) the channel-0 block is built
        from.
    """
    split = _bond_locality_split(interactions)
    if "PairLift" in split.offsite_types:
        logger.info(
            "[eliashberg] bond_channels=true: off-site PairLift declarations "
            "carry no longitudinal (spin/charge) content and are not part of "
            "the bond channels.")
    topo = bond_channels.resolve_bond_topology(
        split.whole_tbl, geom_info["rvec"], norb, max_shells=bond_max_shells,
        active_types=bond_channels.longitudinal_bond_types(split.offsite_types))
    return topo, bond_channels.BondSetView(topo), split


def _build_bond_vertices(topo, split, norb, spatial_shape):
    """The enlarged spin/charge vertices ``S``, ``C`` on the bond basis.

    Built exactly as the FLEX bond gate builds them: the channel-0 block is
    the locality-split pair-space matrix
    (:func:`hwave.solver.offsite.sc_matrices_from_split`; the off-site
    Hartree ``V(q)`` of CoulombInter / Hund / Ising sits in its density
    slots, the Fock crossing of an ``R = 0`` row in its cross slots), and
    every off-site row of every type the topology carries enters through
    ``bond_channels.build_sc_bond_channel`` (the bond-diagonal Fock
    crossing; the PairHop bond entries).

    Returns
    -------
    S, C : ndarray, shape (nvol, ND, ND), complex
        ``ND = B * norb**2``, bond-major index ``m * norb**2 + a * norb + b``.
    """
    from hwave.solver.offsite import sc_matrices_from_split
    Nx, Ny, Nz = (int(x) for x in spatial_shape)
    nvol, nd = Nx * Ny * Nz, norb * norb
    S0, C0 = sc_matrices_from_split(
        split, bond_channels._LONGITUDINAL_ACTIVE_TYPES, norb, Nx, Ny, Nz)
    S0 = S0.reshape(nvol, nd, nd)
    C0 = C0.reshape(nvol, nd, nd)
    types = tuple(topo.coeffs)
    S = bond_channels.build_sc_bond_channel(
        topo, S0, "S", types=types, spatial_shape=(Nx, Ny, Nz))
    C = bond_channels.build_sc_bond_channel(
        topo, C0, "C", types=types, spatial_shape=(Nx, Ny, Nz))
    return S, C


def _bond_kernel_parts(S, C, dressed, green_kw, view, norb, spatial_shape,
                       beta, pairing_type, *, topo, g2_tail=False):
    """The static bond Eliashberg operator and its two parts, built from the
    bond pairing pieces of ``hwave.solver.eliashberg.bond_dynamic``.

    Parameters
    ----------
    S, C : ndarray, shape (nvol, ND, ND) or (Nx, Ny, Nz, ND, ND)
        The enlarged bare spin/charge vertices (:func:`_build_bond_vertices`).
    dressed : iterable of (channel, chi)
        ``("spin", chi_s)`` and ``("charge", chi_c)``, each ``chi`` of the
        shape of ``S``: the dressed susceptibilities on the bond basis. Each
        is added to the pairing vertex as it arrives, so a generator keeps
        one channel resident at a time.
    green_kw : ndarray, shape (norb, norb, Nx, Ny, Nz, nmat)
        The full Green function (``sc.py`` layout) of the pair bubble.
    view : BondSetView
        The bond topology.
    pairing_type : {"singlet", "triplet"}
    topo : bond_channels.BondTopology or None
        The topology ``S`` and ``C`` were built on, for the singlet bare
        pairing of its off-site Exchange and PairHop rows
        (``eliashberg_bond.instantaneous_vertex``); ``None`` states that
        there are none.

    The fluctuation vertex ``Gamma = c_s S chi_s S + c_c C chi_c C``
    (``PairVertexAccumulator``) and the bare vertex (``instantaneous_vertex``:
    the channel-0 block, plus the singlet pairing of the off-site Exchange
    and PairHop rows of ``topo``) are applied by ``BondPairKernel`` on a
    one-point frequency axis, where its Matsubara transforms are the
    identity: one matvec is ``-(1/N) sum e^{-ik.R_b} Gamma_ab(k - k')
    e^{-ik'.R_a} X(k')`` with the static pair bubble ``X = G2 phi``,
    ``G2 = T sum_n G(k, iw_n) G(-k, -iw_n)``.

    The instantaneous part acts on the gap projected onto the parity sector
    of the pairing channel (``_project_gap_parity``): the bare Cooper vertex
    acts on the pair amplitude, which is even (singlet) or odd (triplet)
    under ``Delta_ab(k) -> Delta_ba(-k)``. On that sector the projection is
    the identity; on the opposite sector it keeps the bare vertex from
    acting on a pair amplitude the channel cannot carry. The fluctuation
    part is applied unprojected, as on the non-bond static path.

    Returns
    -------
    A, A_fl, A_pp : scipy.sparse.linalg.LinearOperator
        The full kernel ``K = K^fl + K^pp``, its fluctuation part and its
        instantaneous (bare) part, acting on a gap ``(norb, norb, Nx, Ny,
        Nz)`` raveled.
    vec_size : int
        ``norb**2 * Nx * Ny * Nz``.
    """
    from types import SimpleNamespace
    from hwave.solver.eliashberg import bond_dynamic as _eb
    from hwave.solver import flex_bond as _fb

    Nx, Ny, Nz = (int(x) for x in spatial_shape)
    shape = (Nx, Ny, Nz)
    nvol, nd = Nx * Ny * Nz, norb * norb
    S = np.asarray(S)
    C = np.asarray(C)
    ND = S.shape[-1]
    S = S.reshape(nvol, ND, ND)
    C = C.reshape(nvol, ND, ND)
    if pairing_type not in _eb.COEFF:
        raise ValueError(
            "unknown pairing_type '{}'. Use 'singlet' or 'triplet'.".format(
                pairing_type))

    with _fb.BondDeviceContext(np, S, C) as dev:
        acc = _eb.PairVertexAccumulator(
            dev, pairing_types=(pairing_type,), nb=1, nmat=1, nvol=nvol,
            nd=nd, spatial_shape=shape)
        for channel, chi in dressed:
            chi = np.asarray(chi).reshape(1, nvol, ND, ND)
            acc.add_channel(channel, _eb.ArrayBlockSource({"chi": chi}, nd),
                            "chi")
            del chi
        vertex = acc.finish()[pairing_type]
        del acc
    V_inst = _eb.instantaneous_vertex(S, C, nd, pairing_type, shape, topo=topo)

    # the static pair bubble on a one-point frequency axis
    G2 = bond_channels._g2_from_green(green_kw, beta, tail=g2_tail)[..., None]
    # The static path has its own resource preflight (bond_memory_cap_gb),
    # run before the Green function was built; the kernel's admission is
    # therefore not capped a second time.
    common = dict(xp=np, spatial_shape=shape, norb=norb, beta=beta, nfreq=1,
                  residency="host", host_cap=np.inf)
    kernel_fl = _eb.BondPairKernel(vertex, G2, view, **common)
    # The bare vertex lives on channel 0 alone, so its kernel needs only the
    # on-site channel and a zero fluctuation block.
    view0 = SimpleNamespace(n_channels=1, delta_r=((0, 0, 0),), reverse=(0,))
    zero = _eb.PairVertexUniform(
        source=_eb.ArrayBlockSource(
            {"zero": np.zeros((1, nvol, nd, nd), dtype=np.complex128)}, nd),
        slot="zero", B=1, nd=nd, nmat=1, nvol=nvol)
    kernel_pp = _eb.BondPairKernel(zero, G2, view0, V_inst=V_inst, **common)
    del G2

    gap_shape = (norb, norb, Nx, Ny, Nz)
    vec_size = nd * nvol

    def _fl(v):
        return kernel_fl.matvec(v)

    def _pp(v):
        phi = np.asarray(v).reshape(gap_shape)
        return kernel_pp.matvec(_project_gap_parity(phi, pairing_type).ravel())

    def _full(v):
        return _fl(v) + _pp(v)

    A, A_fl, A_pp = (LinearOperator((vec_size, vec_size), matvec=f,
                                    dtype=complex)
                     for f in (_full, _fl, _pp))
    return A, A_fl, A_pp, vec_size


#: what a static bond run past the RPA pole means for its result (issue #345)
_BOND_POLE_CONSEQUENCE = (
    "The pairing eigenvalue lambda of this run is built on that susceptibility, so it "
    "describes no physical state; it is still reported (eigen_converged keeps its meaning "
    "of eigen-solver convergence).")


def _build_bond_operator(topo, split, green_carrier, norb, kx_array, ky_array,
                         kz_array, beta, pairing_type, *, bond_max_shells,
                         bond_memory_cap_gb, precondition_opts=None,
                         green_source=None, g2_tail=False,
                         bond_coeff_tail_source="default"):
    """Build the static bond-resolved Eliashberg operator.

    The static path is the zero-bosonic-frequency member of the bond-resolved
    pairing family the FLEX bond gate and the dynamic bond pairing share, and
    it is built from the same modules:

    1. the bond bubble ``chibar(q)``: ``bubble.bond_bubble_static`` -- the
       ``Omega = 0`` slice of the per-pair pipeline that
       ``flex_bond.assemble_bubble`` streams over every bosonic frequency;
    2. the vertices ``S``, ``C`` (:func:`_build_bond_vertices`);
    3. the RPA dressing ``chi = [1 -/+ chibar W]^{-1} chibar``:
       ``bond_channels.dress_batch`` on a one-slice batch, one channel at a
       time (its conditioning guard names the channel and the q-point of a
       singular denominator);
    4. the fluctuation vertex ``Gamma = c_s S chi_s S + c_c C chi_c C``:
       ``eliashberg_bond.PairVertexAccumulator``;
    5. the bare (instantaneous) vertex: ``eliashberg_bond.instantaneous_
       vertex`` on the channel-0 block;
    6. the kernel: ``eliashberg_bond.BondPairKernel`` on a one-point
       frequency axis (:func:`_bond_kernel_parts`).

    ``topo``/``split`` come from :func:`_resolve_bond_topology` and must have
    been pre-flighted (``_bond_resource_preflight``) by the caller BEFORE
    ``green_carrier`` was built (review fix I-2: the preflight runs before the
    potentially large Green-function allocation).

    ``green_carrier`` : :class:`_BondGreen`
        The bubble reads ONLY ``deflated_kw``/``green0_tail``; the pair
        bubble ``G2`` and the pair weight read ``full_sc``.

    ``precondition_opts`` (from ``_read_bond_config``) is forwarded verbatim
    to ``bond_channels.check_hermitian_preconditions``.

    ``green_source`` is provenance only: the path the Green function was
    loaded from (``[eliashberg] bond_green``), or ``None`` when the caller
    built the bare Green function from the transfer Hamiltonian.

    ``bond_max_shells`` and ``bond_memory_cap_gb`` are likewise PROVENANCE
    ONLY here: the cutoff and the preflight already happened in the caller.

    Returns
    -------
    (A, vec_size) : tuple
        The ``scipy.sparse.linalg.LinearOperator`` kernel and its vector size,
        in the same convention as ``_make_kernel_operator``.
    provenance : dict
        Additive output-provenance record (spec S5).
    attribution : dict
        Context for the ``lambda = lambda^pp + lambda^fl`` decomposition: the
        two kernel parts (``pp``/``fl``), the full kernel and the pair weight
        (``weight``), consumed by ``bond_channels.attribute_lambda``.
    """
    Nx, Ny, Nz = len(kx_array), len(ky_array), len(kz_array)
    shape = (Nx, Ny, Nz)
    view = bond_channels.BondSetView(topo)
    green_kw = green_carrier.full_sc

    logger.info("Computing the bond-resolved bubble chi-bar...")
    chi_bar = bubble.bond_bubble_static(
        green_carrier.deflated_kw, green_carrier.green0_tail, beta, view,
        spatial_shape=shape)

    S, C = _build_bond_vertices(topo, split, norb, shape)
    # the static RPA pole (issue #345): read-only, before the dressing
    pole = bond_channels.check_rpa_pole(
        {ch: bond_channels.static_pole_indicator(chi_bar, W, ch, shape)
         for ch, W in (("spin", S), ("charge", C))},
        where="eliashberg bond_channels (static)", consequence=_BOND_POLE_CONSEQUENCE,
        log=logger)

    def _dressed(chi_bar, S, C):
        # one dressed channel at a time, on a one-slice batch at the zero
        # bosonic frequency (l0 = nmat = 0, so a refusal names index 0)
        for channel, W in (("spin", S), ("charge", C)):
            chi, _ = bond_channels.dress_batch(
                chi_bar[None], W, channel, l0=0, nmat=0, spatial_shape=shape)
            yield channel, chi[0]
            del chi

    logger.info("Building the bond-resolved Eliashberg operator "
                "(pairing_type=%s)...", pairing_type)
    A, A_fl, A_pp, vec_size = _bond_kernel_parts(
        S, C, _dressed(chi_bar, S, C), green_kw, view, norb, shape, beta,
        pairing_type, topo=topo, g2_tail=g2_tail)
    del chi_bar, S, C

    # v1 runtime preconditions of the Hermitian static path (spec S4.5): the
    # pair weight w = GG must be real (Hermitian) and >= 0 and the symmetrized
    # kernel K~ = -sqrt(w) Gamma sqrt(w) Hermitian, else the reported lambda is
    # not a real eigenvalue/Rayleigh quotient. Violations RAISE -- there is no
    # non-Hermitian biorthogonal fallback in v1.
    weight = bond_channels.pair_weight(green_kw, beta, g2_tail=g2_tail)
    diagnostics = bond_channels.check_hermitian_preconditions(
        A, weight, label="eliashberg bond_channels",
        parts={"pp": A_pp, "fl": A_fl},
        **(precondition_opts or {}))
    logger.info(
        "Bond kernel Hermiticity preconditions OK (%s check): "
        "||K~ - K~^dag||_F = %.3e (relative %.3e; the same residual on "
        "M = W K is %.3e), min eig w = %.3e",
        diagnostics["method"],
        diagnostics["kernel_hermiticity_residual_ktilde"],
        diagnostics["kernel_hermiticity_relative_ktilde"],
        diagnostics["kernel_hermiticity_residual"],
        diagnostics["weight_min_eigenvalue"])

    attribution = {
        "full": A,
        "pp": A_pp,
        "fl": A_fl,
        "weight": weight,
        "diagnostics": diagnostics,
    }

    provenance = {
        "bond_channels": True,
        "bond_n_channels": int(view.n_channels),
        "bond_delta_r": [list(dr) for dr in view.delta_r],
        "bond_max_shells": ("all" if bond_max_shells is None
                            else int(bond_max_shells)),
        "bond_memory_cap_gb": float(bond_memory_cap_gb),
        "g2_tail": bool(g2_tail),
        "bond_coeff_tail_requested": float(green_carrier.coeff_tail_requested),
        "bond_coeff_tail_applied": float(green_carrier.coeff_tail_applied),
        "bond_coeff_tail_source": bond_coeff_tail_source,
        "bond_endpoint_convention": TAIL_ENDPOINT_CONVENTION,
        "bond_green_source": green_carrier.source,
        "external_full_green_assumed": (green_carrier.source == "external_npz"),
        "approximation": (
            "static RPA-ladder bond dressing on "
            + ("an EXTERNAL Green function supplied via [eliashberg] "
               "bond_green (its own approximation level -- e.g. FLEX-dressed "
               "and self-consistent -- is that of the file, not of this run); "
               "the bond chi-bar is still built here from that green, so "
               "this is NOT a conserving FLEX result"
               if green_source else
               "the BARE Green function built from the transfer Hamiltonian "
               "(supply [eliashberg] bond_green to feed an external/FLEX "
               "green instead); NOT a conserving FLEX result")
            + " -- absolute lambda is not comparable to FLEX/dynamic values"),
    }
    if green_source:
        provenance["bond_green"] = str(green_source)
    if view.n_channels == 1:
        provenance["collapsed_to_pure_hubbard"] = True
    provenance["kernel_hermiticity_residual"] = "{:.3e}".format(
        diagnostics["kernel_hermiticity_residual"])
    provenance["kernel_hermiticity_residual_ktilde"] = "{:.3e}".format(
        diagnostics["kernel_hermiticity_residual_ktilde"])
    provenance["kernel_hermiticity_relative_ktilde"] = "{:.3e}".format(
        diagnostics["kernel_hermiticity_relative_ktilde"])
    provenance["kernel_hermiticity_method"] = diagnostics["method"]
    provenance["pair_weight_min_eigenvalue"] = "{:.6e}".format(
        diagnostics["weight_min_eigenvalue"])
    # written only past the pole: a stable run's eigenvalue.dat is unchanged
    provenance.update(pole.text_members())

    return (A, vec_size), provenance, attribution


def _build_bond_green(eigenvalues, eigenvectors, mu, beta, nmat,
                      coeff_tail_requested, bond_green_path,
                      spin_orbital_mode=False):
    """Build the bond path's :class:`_BondGreen` carrier.

    Internal branch (``bond_green_path is None``): flattens sc.py's
    ``(Nx, Ny, Nz, p)`` / ``(Nx, Ny, Nz, p, p)`` eigen-decomposition to the
    kernel's flattened-only ``(nvol, p)`` / ``(nvol, p, p)`` contract (C
    order -- the SAME flattening :func:`_green_sc_to_canonical` uses, so
    k-point ordering is consistent by construction), calls
    ``hwave.solver.green.build_green`` with the CALLER's already-computed
    ``mu`` (no mu search is ever re-run here), and converts the canonical
    ``full_kw`` to ``full_sc`` ONCE -- the canonical array is released
    (``del``) immediately after, so the two full-size layouts coexist only
    for the duration of that one conversion. ``coeff_tail_applied ==
    coeff_tail_requested``, ``source = "internal"``.

    External branch (``bond_green_path`` set): loads the npz via
    :func:`_load_green_npz` (promoted to complex128 at load if needed),
    ``full_sc = loaded``, ``deflated_kw = _green_sc_to_canonical(full_sc)``,
    ``green0_tail = None``, ``coeff_tail_applied = 0.0``, ``source =
    "external_npz"``. A nonzero ``coeff_tail_requested`` logs a WARNING
    (the tail correction is not applicable to an externally supplied Green
    function, whose high-frequency asymptote is not known here) rather
    than raising -- recorded as requested-but-not-applied, not an error.

    Parameters
    ----------
    eigenvalues : ndarray, shape (Nx, Ny, Nz, p)
        Only used for its VALUES on the internal branch; on the external
        branch only its SHAPE is read (to validate the loaded file's
        ``norb``/``Nx``/``Ny``/``Nz`` against the model).
    eigenvectors : ndarray, shape (Nx, Ny, Nz, p, p)
        As above.
    mu : float
        The caller's already-determined chemical potential.
    beta : float
        Inverse temperature.
    nmat : int
        Target Matsubara-frequency count for an internal build; unused on
        the external branch (the file's own frequency count wins -- the
        caller is responsible for reconciling the two, exactly as the
        pre-carrier code did).
    coeff_tail_requested : float
        The resolved ``[eliashberg] bond_coeff_tail`` value.
    bond_green_path : str or None
        ``[eliashberg] bond_green``; ``None`` selects the internal build.
    spin_orbital_mode : bool
        ``False`` (the non-spin-orbital bond path): an external file is
        checked for spin-orbital markers first and refused (the #271
        preflight does not run on the bond path). ``True`` skips the check.

    Returns
    -------
    _BondGreen
    """
    if bond_green_path is not None:
        if not spin_orbital_mode:
            _refuse_so_green_file(bond_green_path)
        Nx, Ny, Nz, norb = eigenvalues.shape
        green_sc, nfreq = _load_green_npz(
            bond_green_path, norb, Nx, Ny, Nz, label="bond_green")
        _validate_bond_green_nmat(nfreq, bond_green_path)
        if green_sc.dtype != np.complex128:
            green_sc = green_sc.astype(np.complex128)
        if coeff_tail_requested != 0.0:
            logger.warning(
                "[eliashberg] bond_coeff_tail = %s requested but "
                "bond_green = '%s' supplies an EXTERNAL Green function: "
                "the tail correction needs the analytic high-frequency "
                "asymptote of an internally built bare Green function, "
                "which is not available for a file. NOT applying it here "
                "(recorded as requested but not applied).",
                coeff_tail_requested, bond_green_path)
        return _BondGreen(
            full_sc=green_sc,
            deflated_kw=_green_sc_to_canonical(green_sc),
            green0_tail=None,
            coeff_tail_requested=float(coeff_tail_requested),
            coeff_tail_applied=0.0,
            source="external_npz")

    Nx, Ny, Nz, p = eigenvalues.shape
    nvol = Nx * Ny * Nz
    ev_flat = eigenvalues.reshape(nvol, p)
    evec_flat = eigenvectors.reshape(nvol, p, p)
    full_kw, deflated_kw, green0_tail = green_mod.build_green(
        ev_flat, evec_flat, mu, beta, nmat, coeff_tail_requested)
    full_sc = _green_canonical_to_sc(full_kw, (Nx, Ny, Nz))
    del full_kw
    return _BondGreen(
        full_sc=full_sc,
        deflated_kw=deflated_kw,
        green0_tail=green0_tail,
        coeff_tail_requested=float(coeff_tail_requested),
        coeff_tail_applied=float(coeff_tail_requested),
        source="internal")


class _SOBondGreen(NamedTuple):
    """The Green function of a static spin-orbital bond run.

    ``carrier`` is the :class:`_BondGreen` of the run (``full_sc`` ``(nd,
    nd, Nx, Ny, Nz, nmat)`` in spin-block order, on the centered fermionic
    grid ``w_n = (2n + 1 - nmat) pi / beta``): built on the canonical band
    (``source = "internal"``), or read from ``[eliashberg] bond_green``
    (``source = "external_npz"``, no tail correction). ``fp_green`` is the
    digest of the Green function the kernel uses (:func:`_so_green_fp`;
    for a borrowed one projected onto the S_z sectors, the projected one);
    ``conserving`` whether the band conserves S_z (decided by the Transfer
    table alone on both paths; then the band's G has exactly zero
    S_z-off-sector entries); ``eigenvalues`` the band ``(Nx, Ny, Nz, nd)``,
    ascending per k (``None`` for a borrowed G: no band is diagonalized);
    ``mu`` the chemical potential (a borrowed G: the file's ``state_mu``).

    A borrowed G also records: ``state``, the producer's recorded state of
    that Green function (the gap archive records it);
    ``fp_green_source``, the file's digest when the G was projected (else
    ``None``); ``leak``, its S_z leak (``None`` when the input does not
    conserve S_z and no leak was measured); ``scf_converged``, the file's
    ``scf_converged`` (``None`` when it records none). On the internal path
    these are ``None`` and ``source`` is ``"internal"``."""
    carrier: object
    fp_green: str
    conserving: bool
    eigenvalues: object
    mu: float
    state: object = None
    fp_green_source: object = None
    leak: object = None
    source: str = "internal"
    scf_converged: object = None


#: INFO (``%s``: path, producer label, leak text, ``"yes"`` / ``"no"``)
_SO_BOND_GREEN_ACCEPTED_INFO = (
    "[eliashberg] bond_green: '%s' (written by %s) passed the archive checks; S_z leak %s; "
    "projected: %s.")
_SO_BOND_GREEN_LEAK_UNCHECKED = "not checked (the input does not conserve S_z)"
#: WARNING (``%s`` path)
_SO_BOND_GREEN_UNCONVERGED = (
    "[eliashberg] bond_green '%s' records scf_converged = false: the FLEX run did not "
    "converge; its Green function is used as stored.")
#: WARNING (``%r`` the configured bond_coeff_tail, ``%s`` path)
_SO_BOND_TAIL_NOT_APPLIED = (
    "[eliashberg] bond_coeff_tail = %r is not applied with bond_green = '%s': the tail "
    "correction is implemented for the band Green function only (bond_coeff_tail_applied = "
    "0), so the bond bubble carries the Matsubara truncation error, which decreases only "
    "about as 1/Nmat; check the eigenvalue against Nmat.")
#: INFO
_SO_BOND_TAIL_DEFAULT_INFO = (
    "[eliashberg] bond_coeff_tail: the default 1.0 applies to the band Green function only; "
    "the Green function read from bond_green is used without a tail correction.")


def _so_conserving_input(so_run, norb_phys):
    """Whether the whole input of ``so_run`` conserves S_z (its band AND its
    interaction records; :func:`hwave.solver.so_bond.input_conserves_sz`):
    the rule that decides whether a borrowed Green function is checked for
    S_z leakage and projected. The ``so_bond_layout`` override never
    decides it."""
    from hwave.solver import so_bond
    return bool(so_bond.input_conserves_sz(
        so_run.band.transfer, so_run.band.extern, so_run.context.compiled,
        norb_phys).conserving)


def _so_bond_green(epsilon_k, so_run, *, norb_phys, filling, beta, nmat,
                   coeff_tail, bond_green=None, coeff_tail_source="default"):
    """The :class:`_SOBondGreen` of a static spin-orbital bond run.

    ``epsilon_k`` is the run's spin-block band (:func:`_so_epsilon_k`) and
    ``so_run`` its :class:`_SORun`. Whether the band conserves S_z is decided
    by the Transfer table alone
    (:func:`hwave.solver.so_bond.band_conserves_sz`; never by the layout of
    the bond objects, so a forced dense layout gives the same G).

    Without ``bond_green`` the Green function is built from the band: the
    band is diagonalized per spin sector when it conserves S_z
    (:func:`hwave.solver.so_bond.canonical_band_eigen`), mu holds
    ``filling`` electrons per orbital per spin, counted over the ``nd = 2
    norb_phys`` spin-orbitals, and ``coeff_tail`` (the resolved
    ``[eliashberg] bond_coeff_tail``, :func:`_resolve_bond_coeff_tail`)
    changes only the deflated Green function of the bubble, never
    ``full_sc``.

    With ``bond_green`` (a path) the Green function is read from that file
    (:func:`_so_borrowed_green`; the band still decides the checks and,
    with the interaction, whether the G is projected onto the S_z sectors,
    :func:`_so_conserving_input`); mu and the state are the file's, and no
    tail correction is applied (``coeff_tail_applied = 0``): the default
    ``bond_coeff_tail`` (``coeff_tail_source = "default"``) is noted at
    INFO, a configured non-zero one gets a WARNING. Once the file is
    accepted, one INFO line gives its leak and whether it was projected,
    and a file that records ``scf_converged = false`` gets one WARNING.

    The carrier's arrays are returned read-only."""
    from hwave.solver import so_bond
    norb_phys = int(norb_phys)
    if coeff_tail is None:
        raise ValueError("_so_bond_green: coeff_tail must be the resolved "
                         "[eliashberg] bond_coeff_tail, got None")
    conserving = bool(so_bond.band_conserves_sz(
        so_run.band.transfer, so_run.band.extern, norb_phys).conserving)
    if bond_green is not None:
        shape = tuple(int(x) for x in np.shape(epsilon_k)[2:5])
        conserving_input = _so_conserving_input(so_run, norb_phys)
        borrowed = _so_borrowed_green(bond_green, so_run, norb_phys=norb_phys, nmat=nmat,
                                      shape=shape, conserving_input=conserving_input)
        full_sc = _green_canonical_to_sc(borrowed.used, shape)
        full_sc.setflags(write=False)
        carrier = _BondGreen(full_sc=full_sc, deflated_kw=borrowed.used, green0_tail=None,
                             coeff_tail_requested=float(coeff_tail), coeff_tail_applied=0.0,
                             source="external_npz")
        logger.info(_SO_BOND_GREEN_ACCEPTED_INFO, bond_green, borrowed.label,
                    _SO_BOND_GREEN_LEAK_UNCHECKED if borrowed.leak is None
                    else "%.3e" % borrowed.leak, "yes" if borrowed.projected else "no")
        if borrowed.scf_converged is False:
            logger.warning(_SO_BOND_GREEN_UNCONVERGED, bond_green)
        if coeff_tail_source == "default":
            logger.info(_SO_BOND_TAIL_DEFAULT_INFO)
        elif float(coeff_tail) != 0.0:
            logger.warning(_SO_BOND_TAIL_NOT_APPLIED, float(coeff_tail), bond_green)
        return _SOBondGreen(
            carrier, borrowed.fp_used, conserving, None, float(borrowed.state.mu),
            state=borrowed.state,
            fp_green_source=borrowed.fp_file if borrowed.projected else None,
            leak=borrowed.leak, source="external_npz", scf_converged=borrowed.scf_converged)
    eigenvalues, eigenvectors = so_bond.canonical_band_eigen(
        epsilon_k, norb_phys, conserving=conserving)
    mu = _determine_mu(eigenvalues, beta, filling, 2 * norb_phys)
    carrier = _build_bond_green(eigenvalues, eigenvectors, mu, beta, nmat,
                                coeff_tail, None, spin_orbital_mode=True)
    del eigenvectors
    # fp_green is taken once here and travels with the carrier to the end of
    # the solve, so the carrier's arrays are read-only from now on
    for a in (carrier.full_sc, carrier.deflated_kw, carrier.green0_tail):
        if a is not None:
            a.setflags(write=False)
    return _SOBondGreen(carrier, _so_green_fp(carrier.full_sc), conserving,
                        eigenvalues, float(mu))


#: ``{method}`` = the configured ``[eliashberg] eigenvalue_method`` (repr)
_SO_BOND_METHOD_REFUSAL = (
    _SO_WHAT + ": [eliashberg] eigenvalue_method = {method} is not available on the "
    "spin-orbital bond path; only \"arnoldi\" is. Set eigenvalue_method = \"arnoldi\" "
    "(or remove the key).")


class _SOBondAdmission(NamedTuple):
    """The admission of a static spin-orbital bond run, decided before the
    vertex, the bubble or the kernel is built: the resolved memory ``caps``
    and their one-line ``cap_summary``, the slot geometry ``geom``, the
    ``layout`` (``"dense"`` or ``"sz_packed"``) and its ``layout_reason``
    line, the S_z ``sector_map`` (packed layout; ``None`` dense) and the
    admitting ``plan`` (:func:`hwave.solver.so_bond.plan_so_bond`), the one
    memory authority of the run."""
    caps: object
    cap_summary: str
    geom: object
    layout: str
    layout_reason: str
    sector_map: object
    plan: object


class _SOBondStatic(NamedTuple):
    """The static spin-orbital bond pairing operator.

    Dense layout: ``operator`` is the kernel on the flat gap ``(nd, nd, Nx,
    Ny, Nz)`` (a ``LinearOperator``), ``basis`` its fermion-antisymmetric
    basis and ``kernel`` ``None``. S_z-packed layout: ``kernel`` is the
    :class:`hwave.solver.so_bond.PackedBondPairKernel` (one operator per pair
    S_z sector, :meth:`sector_op`) and ``operator``, ``basis`` are ``None``.
    ``r_V`` and ``r_G`` are the S_z-conservation residuals of the pairing
    vertex and of the Green function (``r_V = 0`` by construction in the
    packed layout). ``pole`` is the static RPA pole check of the run
    (:func:`_so_bond_static_pole`; ``None`` when not made)."""
    layout: str
    layout_reason: str
    geom: object
    plan: object
    operator: object
    basis: object
    kernel: object
    r_V: float
    r_G: float
    pole: object = None


def _admit_so_bond_static(so_run, *, cfg, cell, nmat, nd, shape, solver_mode,
                          num_eigenvalues, max_iter, spectral_shift, eigenvalue_method, log):
    """The :class:`_SOBondAdmission` of a static spin-orbital bond run.

    The host cap is resolved from ``[eliashberg] bond_memory_cap_gb``
    (``cfg.cap_raw``), the layout is selected from the compiled interaction
    and the Transfer table of ``so_run`` (``cfg.layout_request``; its reason
    line is logged), an ``eigenvalue_method`` other than ``"arnoldi"`` is
    refused in either layout, and the plan is built with the rows of every
    phase of the static path (``spectral_shift`` is the run's raw
    ``[eliashberg] spectral_shift``: the eigensolver rows keep the outputs it
    implies; a Green function read from ``bond_green`` adds its resident read
    row, sized from the zip records and npy headers of the file's members,
    :func:`_so_bond_green_member_bytes`, without reading a member).
    An over-cap run is refused here, before any bond array is allocated."""
    from hwave.solver import so_bond
    info = {}
    caps = so_bond.resolve_host_cap(cfg.cap_raw, key=_SO_BOND_CAP_KEY, log=log,
                                    info_out=info)
    cap_summary = info["summary"]
    topo = so_run.bond.topo
    compiled = so_run.context.compiled
    geom = so_bond.slot_geometry(topo, shape)
    layout, reason = so_bond.select_layout(compiled, so_run.band.transfer,
                                           override=cfg.layout_request)
    log.info(reason)
    if eigenvalue_method != "arnoldi":
        raise ValueError(_SO_BOND_METHOD_REFUSAL.format(method=repr(eigenvalue_method)))
    nvol = int(np.prod(shape))
    plan = so_bond.plan_so_bond(
        nmat=nmat, nvol=nvol, B=geom.B, nd=nd, consumers=("pairing_static",), caps=caps,
        layout=layout, layout_reason=reason, cap_key=_SO_BOND_CAP_KEY,
        shell_rows=so_bond.shell_rows(compiled, topo, cell),
        max_shells_key="[eliashberg] bond_max_shells",
        envelope=so_bond.envelope_text(layout),
        extra_rows=so_bond.static_so_bond_rows(
            nmat=nmat, nvol=nvol, B=geom.B, nd=nd, layout=layout, solver_mode=solver_mode,
            num_eigenvalues=num_eigenvalues, max_iter=max_iter, spatial_shape=shape,
            spectral_shift=spectral_shift),
        resident_rows=(so_bond.static_so_bond_resident_rows() if cfg.green is None
                       else so_bond.static_so_bond_resident_rows(
                           green_source="file",
                           member_bytes=_so_bond_green_member_bytes(cfg.green))))
    sector_map = so_bond.sz_sector_map(geom, nd) if layout == "sz_packed" else None
    return _SOBondAdmission(caps, cap_summary, geom, layout, reason, sector_map, plan)


def _so_bond_static_pole(X, gamma_b, shape):
    """The static RPA pole check of a spin-orbital bond run (issue #345):
    the indicator ``max_q eig[-chibar_B Gamma_B](q, 0)`` of the spinful
    channel, over the dense block or over the S_z sectors of the packed
    layout (whose spectra are together the dense block's), logged and
    returned as the :class:`~hwave.solver.bond_channels.RpaPoleReport`."""
    from hwave.solver import so_bond
    entry, per, worst = so_bond.static_pole_indicator(X, gamma_b, shape)
    if per is not None:
        logger.info("spin-orbital bond (static, S_z-packed): the spinful pole indicator "
                    "per S_z sector: %s", "; ".join("{} {:.6f}".format(k, v[0]) for k, v in per.items()))
        where = ("eliashberg bond_channels (spin-orbital, static, S_z-packed{})".format(
            "" if worst is None else "; S_z sector " + worst))
    else:
        where = "eliashberg bond_channels (spin-orbital, static)"
    return bond_channels.check_rpa_pole({"spinful": entry}, where=where,
                                        consequence=_BOND_POLE_CONSEQUENCE, log=logger)


def _so_bond_static_vertex(so_run, admission, so_green, *, beta, nmat, pole_out=None):
    """``(Gamma_B, V)``: the bond vertex of the run's compiled interaction and
    the static pairing vertex ``V = V_bare + V_fluct``, ``V_bare = 0.25
    (Gamma_B + Gamma_B^P)``, ``V_fluct = -Gamma_B chi_B Gamma_B`` with the
    dressed ``chi_B = [1 + chibar_B Gamma_B]^{-1} chibar_B`` of the static
    bubble of the carrier's deflated Green function (a batch of one at the
    grid index ``l0 = nmat // 2``, the zero bosonic frequency). ``V`` is
    formed in place in ``V_bare``; the temporaries are released before
    return. ``pole_out`` (a dict, optional) receives ``"report"``, the RPA
    pole check of the static bubble (:func:`_so_bond_static_pole`); it is
    made in either case."""
    from hwave.solver import so_bond
    plan, geom = admission.plan, admission.geom
    shape = tuple(geom.spatial_shape)
    nvol = geom.nvol
    carrier = so_green.carrier
    gamma_b = so_bond.so_bond_vertex(so_run.context.compiled, so_run.bond.topo, shape,
                                     plan=plan)
    if admission.layout == "sz_packed":
        sm = admission.sector_map
        sink = so_bond.PackedArrays.zeros(sm, (nvol,))
    else:
        sink = np.zeros((nvol, plan.ND, plan.ND), dtype=np.complex128)
    X = so_bond.so_bond_bubble(carrier.deflated_kw, carrier.green0_tail, beta, geom, shape,
                               plan=plan, freq="static", sink=sink)
    del sink
    report = _so_bond_static_pole(X, gamma_b, shape)
    if pole_out is not None:
        pole_out["report"] = report
    if admission.layout == "sz_packed":
        Xb = so_bond.PackedArrays(sm, tuple(b[None] for b in X.blocks))
        cb = so_bond.so_bond_dress(Xb, gamma_b, plan=plan, l0=nmat // 2, spatial_shape=shape)
        chi_b = so_bond.PackedArrays(sm, tuple(b[0] for b in cb.blocks))
        v_bare, v_fluct = so_bond.so_bond_pairing_vertex_parts(gamma_b, chi_b, geom, plan=plan)
        for b, f in zip(v_bare.blocks, v_fluct.blocks):
            b += f
        del b, f, v_fluct, chi_b, cb, Xb, X
    else:
        chi_b = so_bond.so_bond_dress(X[None], gamma_b, plan=plan, l0=nmat // 2,
                                      spatial_shape=shape)[0]
        v_bare, v_fluct = so_bond.so_bond_pairing_vertex_parts(gamma_b, chi_b, geom, plan=plan)
        v_bare += v_fluct
        del v_fluct, chi_b, X
    V = v_bare
    return gamma_b, V


def _so_bond_static_g2(so_green, *, beta, g2_tail):
    """The static pair propagator ``G2 (nd, nd, nd, nd, Nx, Ny, Nz, 1)`` of
    the carrier's full Green function (``1/beta`` included; the Matsubara
    tail correction when ``g2_tail``)."""
    return bond_channels._g2_from_green(so_green.carrier.full_sc, beta, tail=g2_tail)[..., None]


def _so_bond_static_kernel(so_run, admission, V, G2, so_green, *, beta, norb_phys):
    """The :class:`_SOBondStatic` of the pairing vertex ``V`` and the pair
    propagator ``G2``: the kernel ``out_AD(k) = -(1/N) sum_{alpha, beta}
    e^{-ik.R_beta} sum_{k'} V_{alpha beta; ABCD}(k - k') e^{-ik'.R_alpha}
    F_BC(k')`` on the uniform static grid (no instantaneous vertex), checked
    against the plan's kernel rows before it is built. The kernels' own
    memory checks are off (the plan admitted the run)."""
    from hwave.solver.eliashberg import bond_dynamic as _eb
    from hwave.solver import so_bond, so_pairing
    plan, geom = admission.plan, admission.geom
    shape = tuple(geom.spatial_shape)
    nvol, B = geom.nvol, geom.B
    nd = 2 * int(norb_phys)
    plan.require("so_bond_static_kernel", nvol=nvol, B=B, nd=nd, ND=B * nd * nd,
                 layout=admission.layout, n_gap=nd * nd * nvol)
    r_G = so_bond.green_leak(so_green.carrier.full_sc, norb_phys, axes=(0, 1))
    if admission.layout == "sz_packed":
        sm = admission.sector_map
        kernel = so_bond.PackedBondPairKernel(
            so_bond.PackedArrays(sm, tuple(b[None] for b in V.blocks)), G2, geom, sm,
            plan=plan, spatial_shape=shape, beta=beta, nfreq=1, host_cap=None)
        return _SOBondStatic(admission.layout, admission.layout_reason, geom, plan, None,
                             None, kernel, 0.0, float(r_G))
    vertex = _eb.PairVertexUniform(_eb.ArrayBlockSource({"V": V[None]}, nd * nd), "V", B,
                                   nd * nd, 1, nvol)
    kernel = _eb.BondPairKernel(vertex, G2, geom, xp=np, spatial_shape=shape, norb=nd,
                                beta=beta, nfreq=1, V_inst=None, residency="host",
                                host_cap=np.inf)
    n = nd * nd * nvol
    operator = LinearOperator((n, n), matvec=kernel.matvec, dtype=complex)
    basis = _SOAntisymmetricBasis((nd, nd) + shape)
    r_V = so_pairing.bond_sz_conservation_residual(V, nd=nd, B=B)
    return _SOBondStatic(admission.layout, admission.layout_reason, geom, plan, operator,
                         basis, None, float(r_V), float(r_G))


def _build_so_bond_operator(so_run, admission, so_green, *, beta, nmat, g2_tail, norb_phys):
    """The :class:`_SOBondStatic` of an admitted static spin-orbital bond run:
    :func:`_so_bond_static_vertex`, :func:`_so_bond_static_g2` and
    :func:`_so_bond_static_kernel` in order (each phase releases its
    temporaries before the next starts)."""
    pole = {}
    gamma_b, V = _so_bond_static_vertex(so_run, admission, so_green, beta=beta, nmat=nmat,
                                        pole_out=pole)
    # Gamma_B is released here; the admission rows keep counting it as resident
    # to the end of the solve (conservative)
    del gamma_b
    G2 = _so_bond_static_g2(so_green, beta=beta, g2_tail=g2_tail)
    return _so_bond_static_kernel(so_run, admission, V, G2, so_green, beta=beta,
                                  norb_phys=norb_phys)._replace(pole=pole["report"])


def _bond_diagnostics_record(provenance, gap, eigenvalues, weight, norb,
                             kx_array, ky_array, kz_array,
                             deg_tol=_BOND_DIAGNOSTICS_DEG_TOL,
                             sector_eigenvalues=None, leader_available=True):
    """Opt-in character analysis of the leading bond state (``[eliashberg]
    bond_diagnostics``; spec S7.7).

    Writes into the ADDITIVE provenance record (``#`` comment lines of
    ``eigenvalue.dat``, see ``_save_results``) two things a single solve can
    honestly report about its leading state:

    1. ``bond_diagnostics_harmonics`` -- the odd-harmonic decomposition of the
       returned gap on ``_odd_harmonic_basis`` in the ``sqrt(GG)`` metric
       (``bond_channels.harmonic_decomposition``). Each number is the fraction
       of that harmonic captured by the state, in ``[0, 1]``; it is the
       V-sweep "character" quantity of spec S7.7 evaluated at ONE point.
    2. ``bond_diagnostics_eigenvalue_clusters`` -- the near-degenerate grouping
       of the computed eigenvalues (``bond_channels.cluster_eigenvalues``), so
       a user can see whether the leading state is a doublet (the odd E
       doublet of the square lattice) or a singlet branch. Only present when a
       spectrum was actually computed (``solver_mode`` includes
       ``"eigenvalue"``); the iteration path returns one vector and no
       spectrum, and a one-eigenvalue "cluster list" would be a fake.

    The ``lambda^pp`` / ``lambda^fl`` attribution sits in the same block and is
    written by the caller regardless of this flag.

    **Eigenvector orientation.** ``eigenvalues``/``gap`` come from
    ``_solve_eigenvalue``, which already returns eigenvectors as one gap array
    PER LEADING INDEX (it transposes ``scipy.sparse.linalg.eigs``'s COLUMNS);
    the row stack this function hands to ``harmonic_decomposition`` therefore
    holds one eigenvector per row, which is the convention the
    ``bond_channels`` subspace helpers (and ``track_subspace``) consume.

    Multi-point tracking across a ``V`` sweep -- following ONE invariant
    subspace through the spectrum rather than characterizing one solve -- is
    ``bond_channels.track_subspace``; a single ``calc_eliashberg`` run has only
    one point, so that helper is the documented API for sweep drivers instead.
    """
    provenance["bond_diagnostics"] = True

    basis = _odd_harmonic_basis(norb, kx_array, ky_array, kz_array)
    if basis and gap is not None:
        vec = np.asarray(gap).ravel()[np.newaxis, :]
        harmonics = bond_channels.harmonic_decomposition(vec, basis, weight)
        provenance["bond_diagnostics_harmonics"] = ", ".join(
            "{}={:.6f}".format(name, harmonics[name])
            for name in sorted(harmonics))
        provenance["bond_diagnostics_harmonics_note"] = (
            "odd-harmonic decomposition of the LEADING (returned) gap on the "
            "basis {sin kx, sin ky, sin kx (cos kx - cos ky), "
            "sin ky (cos ky - cos kx)} in the sqrt(GG) pair-weight metric; "
            "each value is the fraction of that harmonic captured by the "
            "state, in [0, 1] (spec S7.7). p_x/p_y/f_x/f_y span the SAME odd "
            "2D representation on a square lattice, so this is a character "
            "report, not a point-group separation")
        logger.info("bond diagnostics -- odd-harmonic content of the leading "
                    "state: %s", provenance["bond_diagnostics_harmonics"])
    elif not basis:
        provenance["bond_diagnostics_harmonics"] = (
            "n/a (every odd harmonic vanishes on this k-grid)")

    if not leader_available:
        # no reported leader (NaN): there is no state to describe
        provenance["bond_diagnostics_leading_cluster"] = "n/a"
        provenance["bond_diagnostics_leading_cluster_size"] = "n/a"
        provenance["bond_diagnostics_leading_cluster_source"] = (
            "n/a (no leading eigenvalue was obtained)")
    elif sector_eigenvalues is not None and len(sector_eigenvalues) > 0:
        # the reported leader came from the channel's parity-sector solve:
        # its degeneracy is read off THAT spectrum (index 0 is the leader),
        # not off the global listing, whose rows may describe other states
        clusters = bond_channels.cluster_eigenvalues(
            np.asarray(sector_eigenvalues), deg_tol=deg_tol)
        leading = next(c for c in clusters if 0 in c)
        if eigenvalues is not None and len(eigenvalues) > 0:
            # the table's own groups, for reference (indices into its rows;
            # row 0 is the leader, whose state the global rows may list
            # again, so this grouping is not a degeneracy count)
            provenance["bond_diagnostics_eigenvalue_clusters"] = str(
                bond_channels.cluster_eigenvalues(np.asarray(eigenvalues),
                                                  deg_tol=deg_tol))
        provenance["bond_diagnostics_sector_eigenvalue_clusters"] = str(clusters)
        provenance["bond_diagnostics_leading_cluster"] = str(leading)
        provenance["bond_diagnostics_leading_cluster_size"] = len(leading)
        provenance["bond_diagnostics_leading_cluster_source"] = (
            "channel parity-sector spectrum")
        provenance["bond_diagnostics_deg_tol"] = "{:.3e}".format(deg_tol)
        provenance["bond_diagnostics_clusters_note"] = (
            "near-degenerate groups (single-linkage on |lambda_i - lambda_j| "
            "<= deg_tol) of the eigenvalues of the requested channel's "
            "parity-sector solve, which produced the reported leader (its "
            "index 0); bond_diagnostics_leading_cluster is the group "
            "containing it, i.e. the degeneracy of the state this run "
            "returns; size 2 is the odd E doublet")
        logger.info("bond diagnostics -- sector eigenvalue clusters (deg_tol "
                    "%.1e): %s; the returned state's cluster is %s (size %d)",
                    deg_tol, clusters, leading, len(leading))
    elif eigenvalues is not None and len(eigenvalues) > 0:
        clusters = bond_channels.cluster_eigenvalues(np.asarray(eigenvalues),
                                                     deg_tol=deg_tol)
        # The cluster of the RETURNED eigenpair, i.e. of row 0 -- NOT
        # clusters[0]. cluster_eigenvalues orders by descending Re lambda over
        # every computed eigenpair, while row 0 is the eigenpair the run
        # actually reports (the channel's leader is row 0 of the table); an
        # opposite-parity mode with a larger Re lambda sits below it but is
        # still clustered here. Reporting clusters[0]'s size as "the leading
        # state's degeneracy" would then describe a different state.
        leading = next(c for c in clusters if 0 in c)
        provenance["bond_diagnostics_eigenvalue_clusters"] = str(clusters)
        provenance["bond_diagnostics_leading_cluster"] = str(leading)
        provenance["bond_diagnostics_leading_cluster_size"] = len(leading)
        provenance["bond_diagnostics_leading_cluster_source"] = (
            "eigen-analysis table")
        provenance["bond_diagnostics_deg_tol"] = "{:.3e}".format(deg_tol)
        provenance["bond_diagnostics_clusters_note"] = (
            "near-degenerate groups of the COMPUTED eigenvalues (indices into "
            "the rows below, single-linkage on |lambda_i - lambda_j| <= "
            "deg_tol, ordered by descending Re lambda over ALL computed "
            "eigenpairs -- including any opposite-parity mode below row 0). "
            "bond_diagnostics_leading_cluster is the group containing ROW 0, "
            "i.e. the degeneracy of the state this run actually returns; size "
            "2 is the odd E doublet")
        logger.info("bond diagnostics -- eigenvalue clusters (deg_tol %.1e): "
                    "%s; the returned state's cluster is %s (size %d)",
                    deg_tol, clusters, leading, len(leading))
    return provenance
