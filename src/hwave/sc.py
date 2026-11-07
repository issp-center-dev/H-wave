"""Linearized Eliashberg equation solver for superconducting instability analysis.

This module implements a post-processing tool that uses the bare susceptibility
chi0q from H-wave's RPA solver to solve the linearized Eliashberg equation
and analyze superconducting instabilities.

The tool reads chi0q.npz output from RPA calculations along with interaction
parameters, reconstructs the Green's function, computes RPA vertices,
and solves the Eliashberg equation via self-consistent iteration or
eigenvalue analysis.
"""

from __future__ import annotations

import os
import sys
import logging
import zipfile
from dataclasses import dataclass
from typing import NamedTuple

import numpy as np

from hwave.solver.vertex_table import sc_coefficients
from hwave.solver.kgrid import reverse_fft_axes
from hwave.solver.declarations import symmetrise_k
from hwave.solver import npy_header as _npy_header
from numpy.fft import fftn, ifftn
from scipy.optimize import brentq
from scipy.sparse.linalg import LinearOperator, eigs, bicgstab, gmres, lgmres

import hwave
import hwave.qlmsio.wan90 as wan90
from hwave.solver.rpa import validate_chi0q_index_convention, TAIL_ENDPOINT_CONVENTION
from hwave.solver.ir_axis import is_ir_native, ir_native_meta
from hwave.solver import backend
from hwave.solver import bond_channels
from hwave.solver import bubble
from hwave.solver import green as green_mod


# Moved to hwave.solver.eliashberg (dev/sc-split-plan.md). Re-exported here so
# that `hwave.sc.<name>` keeps resolving for the tests and the solver modules
# until they import from the new modules directly.
from hwave.solver.eliashberg.vertex import (  # noqa: F401,E402
    _build_interaction_k,
    _symmetrise_interactions_k,
    _accumulate_coeff,
    SiteSplitInteractionK,
    _offsite_rows,
    _nonbond_site,
    _build_interaction_k_for_sites,
    _build_sc_matrices_all_q,
    _build_sc_matrices,
    _declarations_partner_closed,
    _compute_vertices,
    _compute_vertices_simple,
    _compute_vertices_general,
    _off_density_sc_weight,
    _build_vertex_sc_matrices,
    _reject_reduced_flex_unsupported,
    _warn_reduced_flex_missing_components,
    _triplet_bare_vertex,
    _compute_vertices_flex,
    _offsite_mirror_asymmetry,
    _decode_version_stamp,
    placement_sensitive_types,
    check_placement_stamp,
    PLACEMENT_VERSION,
    _PLACEMENT_TYPES,
    _REDUCED_FLEX_PARTIAL,
    _REDUCED_FLEX_REJECTED,
    _NONBOND_FOCK_TYPES,
)
from hwave.solver.eliashberg.config import (  # noqa: F401,E402
    _bond_bool_option,
    _bond_float_option,
    _bond_int_option,
    _eliashberg_frequency,
    _validate_dynamic_prereqs,
    _warn_if_static_ignores_channel_flags,
    _resolve_bond_coeff_tail,
    _read_bond_config,
    _validate_bond_prereqs,
    _reject_ir_native,
    _check_nonbond_offsite,
    _check_spin_block_discarded,
    _validate_bond_green_nmat,
    _finite_positive_float64,
    _coerce_run_beta,
    _coerce_config_bool,
    _coerce_g2_tail,
    _resolve_spin_orbital_flag,
    _warn_if_g2_tail_outside_asymptotic_regime,
    _warn_if_g2_indefinite,
    _resolve_init_gap,
    _validate_spectral_shift,
    _G2_CHECK_MAX_WORK,
    _G2_CHECK_MAX_BYTES,
    _G2_TAIL_EDGE_DEV_THRESHOLD,
    _SPIN_DISCARD_ROUNDOFF_RATIO,
    _NONBOND_REFUSED_OFFSITE,
    _BOND_IMAG_TOL,
    _BOND_OPTION_KEYS,
    _BOND_MEMORY_CAP_GB,
    _BOND_BOOL_TRUE,
    _BOND_BOOL_FALSE,
)
from hwave.solver.eliashberg.band import (  # noqa: F401,E402
    _MeanFieldBand,
    _one_body_band_sources,
    _warn_one_body_band_ignored,
    _spin_dependence,
    _refuse_spin_dependent,
    _AssembledH0Ready,
    _checked_band,
    _check_green_init_spin_free,
    _build_hamiltonian_k,
    _calc_eigenvalues,
    _determine_mu,
    SPIN_FREE_ATOL,
    _EPS_K_HAM_ARRAYS,
    HERMITIAN_RTOL,
    HERMITIAN_ATOL,
    _FLEX_BAND_REASON,
    _ONE_BODY_FILE_KEYS,
)
from hwave.solver.eliashberg.kernel import (  # noqa: F401,E402
    _green_sc_to_canonical,
    _green_canonical_to_sc,
    _BondGreen,
    _calc_g2,
    _eliashberg_kernel_fft,
    _make_kernel_operator,
    _G2_BETA_MAX,
)

from hwave.solver.eliashberg.eigen import (  # noqa: F401,E402
    _NonFiniteEigenvalues,
    _is_arpack_zero_start_error,
    _null_flag,
    _operator_is_exactly_zero,
    _initialize_gap,
    _odd_harmonic_basis,
    _solve_iteration,
    _arpack_start,
    _order_by_seed_overlap,
    _order_eigenpairs,
    _reverse_k_and_orbital,
    _project_gap_parity,
    _is_gap_parity,
    _reorder_eigenpairs_by_parity,
    _parity_sector_basis,
    _listing_with_leader,
    _channel_eigenpairs,
    _shift_from_eigenvalues,
    _normalize_spectral_shift,
    _resolve_spectral_shift,
    _ShiftProbeFailed,
    _estimate_spectral_radius,
    _resolve_iteration_shift,
    _signed_rayleigh,
    _validate_shifted_eigenvalue,
    _stall_advice,
    _solve_leading,
    _solve_eigenvalue,
    _eigs_shift_invert,
    _solve_subspace_iteration,
    _solve_shifted_bicg,
    _SOAntisymmetricBasis,
    _SEED_EIGENVECTOR_TOL,
    _DENSE_FALLBACK_MAX,
    _TWO_CYCLE_WINDOW,
    _COLLAPSE_REL,
    _NULL_PROBE_COUNT,
    _DIFF_ROUNDOFF,
    _SEED_ANNIHILATED_REL,
    _RATE_CLIP,
    _AUTO_SHIFT_MARGIN,
    _SEED_ADMIX_SEED,
    _RATE_WINDOW,
    _RAYLEIGH_RESIDUAL_TOL_FACTOR,
    _RAYLEIGH_RESIDUAL_ATOL,
    _AUTO_SHIFT_PROBE_STEPS,
    _PARITY_GUARD_PROBES,
    _ARPACK_START_SEED,
    _ODD_HARMONICS,
    _NULL_PROBE_SEED,
)

logger = logging.getLogger("hwave_sc")


# Default [eliashberg] solver_mode of hwave_sc (static, dynamic, spin-orbital
# and the bond file entry) and of the in-process FLEX bond pairing step. Both
# are the ARPACK eigenvalue solver: it searches the whole spectrum for the
# largest real part, while the power iteration ("iteration", opt-in) follows
# one symmetry sector and can only report what it can certify.
DEFAULT_SOLVER_MODE = "eigenvalue"
IN_PROCESS_DEFAULT_SOLVER_MODE = "eigenvalue"

# Exit status of the hwave_sc command when the eigen-solver did not converge
# (the reported leading eigenvalue is NaN; the outputs are still written).
EXIT_EIGEN_NOT_CONVERGED = 3

# Default Matsubara grid size (mode.param.Nmat). Single definition so the
# legacy-file fallback in _static_freq_position and the main solver loop can
# never silently disagree on the grid.
_DEFAULT_NMAT = 1024


# Degeneracy tolerance of the opt-in bond diagnostics' eigenvalue clustering
# ([eliashberg] bond_diagnostics; spec S7.7 default).
_BOND_DIAGNOSTICS_DEG_TOL = 1.0e-3


# ---------------------------------------------------------------------------
# Eliashberg frequency mode dispatch
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Bond-resolved interaction channels (opt-in; spec S3-S5)
# ---------------------------------------------------------------------------


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
    bond pairing pieces of ``hwave.solver.eliashberg_bond``.

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
    from hwave.solver import eliashberg_bond as _eb
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

    return (A, vec_size), provenance, attribution


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _static_freq_position(freq_index, nfreq, config_nmat, file_name,
                          file_nmat=None):
    """Locate the zero bosonic frequency along the chi0q frequency axis.

    RPA writes ``freq_index`` (the ORIGINAL Matsubara indices, zero frequency
    at original index Nmat//2) into chi0q.npz; the matsubara_frequency option
    can restrict the stored axis, so the static slice is generally NOT at
    nfreq//2 of the stored array.  Newer files also record ``nmat`` (the full
    grid size of the producing run), which resolves the position without any
    assumption about the hwave_sc configuration.

    Parameters
    ----------
    freq_index : array-like or None
        The freq_index metadata from the npz file (None for legacy files).
    nfreq : int
        Length of the stored frequency axis.
    config_nmat : int
        Nmat from mode.param of the hwave_sc run (fallback only).
    file_name : str
        For error messages.
    file_nmat : int or None
        The nmat metadata from the npz file, if present.  When present it
        is the authoritative full-grid reference everywhere in this
        function -- including the two metadata-poor fallbacks, where a
        stored axis equal to the producer's own nmat is provably a full
        grid even if the consuming run's Nmat differs.

    Returns
    -------
    int or None
        Position of the zero bosonic frequency along the stored axis, or
        None when the file carries no usable metadata -- the caller must
        then slice the CENTER of the frequency axis it actually uses (only
        the caller can identify that axis reliably, e.g. for the 6D
        reference format).
    """
    # The full-grid reference: the PRODUCER's own recorded grid size when the
    # file carries one (authoritative provenance -- a stored axis equal to it
    # is provably a full grid whatever the consuming run's Nmat is), else the
    # configured Nmat.
    ref_nmat = int(file_nmat) if file_nmat is not None else config_nmat
    ref_desc = ("the file's recorded nmat = {}".format(ref_nmat)
                if file_nmat is not None
                else "mode.param.Nmat = {}".format(ref_nmat))
    # the same reference phrased as a requirement, for the message that
    # states what nfreq would have to equal
    ref_req = ("the file's recorded nmat ({})".format(ref_nmat)
               if file_nmat is not None
               else "mode.param.Nmat ({})".format(ref_nmat))
    if file_nmat is not None:
        ref_advice = (
            "The file records nmat = {}, so its stored axis is a "
            "restriction of that grid: regenerate it holding the full "
            "grid, or with a matsubara_frequency range covering nmat//2 "
            "(which records freq_index).".format(ref_nmat))
    else:
        ref_advice = (
            "If the file holds a full grid, set mode.param.Nmat = {}; "
            "otherwise regenerate it with a newer version (which records "
            "nmat in the file).".format(nfreq))

    if freq_index is None:
        # Pre-provenance file: with no metadata at all, the stored axis can
        # only be centered safely when the choice of slice is unambiguous.
        # A single-frequency axis (nfreq == 1) is such a case: index 0 is the
        # ONLY slice, so centering cannot pick a different (finite) frequency
        # -- this is the common static-only chi0q reduced to its zero
        # bosonic-frequency component.  A full grid of the reference size
        # (nfreq == ref_nmat) is the other.  Any other length could be a
        # matsubara_frequency restriction of some run's grid, so nfreq//2
        # might point at a finite frequency -- refuse to guess.
        if nfreq == 1:
            # Ruling (#186 review): a static-only chi0q holds exactly the
            # zero bosonic frequency by construction, so keep accepting the
            # singleton -- but state the assumption instead of reporting a
            # generic centering.
            logger.warning(
                "chi0q file '{}' has no freq_index metadata and a single "
                "stored frequency slice; assuming it is the zero bosonic "
                "frequency (a static-only file) -- regenerate the file "
                "with a newer version to record freq_index/nmat if it is "
                "not.".format(file_name))
            return None
        if nfreq == ref_nmat:
            logger.warning(
                "chi0q file '{}' has no freq_index metadata; using the "
                "center of the stored frequency axis as the static slice."
                .format(file_name))
            return None
        raise ValueError(
            "chi0q file '{}' has no freq_index metadata and its stored "
            "frequency axis length {} does not match {}, so the "
            "zero-frequency position cannot be determined: centering the "
            "axis could pick a finite frequency of a restricted grid. {}"
            .format(file_name, nfreq, ref_desc, ref_advice))

    freq_index = np.asarray(freq_index).ravel()
    if freq_index.size == 0:
        raise ValueError(
            "chi0q file '{}' records no frequencies (matsubara_frequency="
            "\"none\"?); it cannot provide the static susceptibility."
            .format(file_name))
    if freq_index.size != nfreq:
        # Legacy FLEX files store the FULL grid but a restricted freq_index
        # (FLEX never applied the matsubara_frequency filter): the DATA axis
        # is authoritative, so fall back to the pre-metadata behavior -- but
        # ONLY when the data axis is provably a full grid of the reference
        # size (nfreq == ref_nmat, the file's own recorded nmat when it has
        # one, else the configured Nmat) AND freq_index is a 0-based subset
        # of it (set(freq_index) <= set(range(nfreq))).  Otherwise the stored
        # axis could be a restriction of some other grid and centering could
        # pick a finite frequency, so refuse to guess.
        if nfreq == ref_nmat and set(freq_index.tolist()) <= set(
                range(nfreq)):
            logger.warning(
                "chi0q file '{}': freq_index length {} does not match the "
                "frequency axis length {}; ignoring the metadata and using "
                "the center of the stored frequency axis as the static "
                "slice.".format(file_name, freq_index.size, nfreq))
            return None
        raise ValueError(
            "chi0q file '{}': freq_index length {} does not match the "
            "frequency axis length {}, and the axis is not provably a full "
            "grid (nfreq must equal {}, with a 0-based freq_index), so the "
            "zero-frequency position cannot be determined: centering the "
            "axis could pick a finite frequency of a restricted grid. {}"
            .format(file_name, freq_index.size, nfreq, ref_req, ref_advice))

    def _find(nmat_orig):
        # the zero bosonic frequency has ORIGINAL index nmat_orig//2
        pos = np.where(freq_index == nmat_orig // 2)[0]
        if pos.size == 0:
            raise ValueError(
                "chi0q file '{}' does not contain the zero bosonic frequency "
                "(original index Nmat//2 = {}; stored freq_index = {}..{}). "
                "Regenerate chi0q with a matsubara_frequency range covering "
                "Nmat//2.".format(file_name, nmat_orig // 2,
                                  freq_index[0], freq_index[-1]))
        return int(pos[0])

    if file_nmat is not None:
        # authoritative: the producing run recorded its own grid size
        return _find(int(file_nmat))

    if np.array_equal(freq_index, np.arange(nfreq)):
        # Legacy file (no nmat metadata) with a contiguous 0-based
        # freq_index: EITHER a full grid of its own (zero at nfreq//2) OR a
        # matsubara_frequency=[0,K] restriction of SOME run's grid (config
        # or larger).  Resolve only when the readings coincide; silently
        # picking either would return a finite frequency whenever the other
        # reading is the true one.
        if nfreq == config_nmat:
            return nfreq // 2
        raise ValueError(
            "chi0q file '{}' is ambiguous: freq_index = 0..{} with "
            "mode.param.Nmat = {} could be either a full {}-point grid "
            "(zero frequency at {}) or a [0,{}] restriction of a run with "
            "a different Nmat (zero frequency elsewhere or absent). If the "
            "file holds a full grid, set mode.param.Nmat = {}; otherwise "
            "regenerate it with a newer version (which records nmat in "
            "the file)."
            .format(file_name, nfreq - 1, config_nmat, nfreq, nfreq // 2,
                    nfreq - 1, nfreq))

    # Legacy restricted range: fall back to the configured Nmat, which must
    # be the producing run's grid size (the standard workflow shares one
    # TOML).  Indices reaching past the config grid disprove that reading,
    # so refuse to guess (looking up config_nmat//2 could land on a finite
    # frequency of the true, larger grid).
    if int(freq_index.max()) >= config_nmat:
        raise ValueError(
            "chi0q file '{}': freq_index reaches {} >= mode.param.Nmat = "
            "{}, so the file was produced by a run with a different Nmat "
            "and the zero-frequency position cannot be determined. Set "
            "mode.param.Nmat to the producing run's value, or regenerate "
            "the file with a newer version (which records nmat)."
            .format(file_name, int(freq_index.max()), config_nmat))
    return _find(config_nmat)


def _read_freq_meta(data):
    """Decode the frequency metadata of a chi0q/chiq npz file.

    Returns
    -------
    (freq_index or None, nmat or None)
        Single definition shared by all loaders so chi0q and chiq_s/chiq_c
        can never interpret the same file's frequency axis differently.
    """
    freq_index = data["freq_index"] if "freq_index" in data else None
    file_nmat = int(data["nmat"]) if "nmat" in data else None
    return freq_index, file_nmat


def _require_so_chi0q_marker(data, file_name, hint=None, producer="RPA"):
    """Refuse a chi0q that does not record ``spin_orbital = True``.

    Every RPA run writes ``index_convention = "spin_block"``, so that tag does
    not identify a spin-orbital file, and the shape cannot either: a
    non-spin-orbital chi0q over 2N orbitals has the shape of a spin-orbital
    one over N physical orbitals. The RPA therefore writes an explicit
    ``spin_orbital = True`` marker, ONLY in a spin-orbital run (a
    non-spin-orbital archive carries no such key); spin-orbital mode accepts
    only a 0-d boolean ``True``. Checked on the open npz before the chi0q
    array itself is read. Also used for the chi0q companion of a
    spin-orbital FLEX run (``producer = "FLEX"``)."""
    if hint is None:
        hint = ("Regenerate the file with a spin-orbital RPA run "
                "([mode] enable_spin_orbital = true, calc_scheme = \"general\") "
                "that computes chi0q itself or reads a chi0q_init marked "
                "spin_orbital = True.")
    if "spin_orbital" not in data:
        raise ValueError(
            "chi0q file '{}' has no spin_orbital marker, so it cannot be told "
            "from a non-spin-orbital chi0q of the same shape. {}".format(
                file_name, hint))
    val = np.asarray(data["spin_orbital"])
    if val.ndim != 0 or val.dtype != np.bool_:
        raise ValueError(
            "chi0q file '{}': malformed spin_orbital marker ({!r}); expected a "
            "boolean. {}".format(file_name, data["spin_orbital"], hint))
    if not bool(val[()]):
        raise ValueError(
            "chi0q file '{}' records spin_orbital = False: it was not written "
            "by a spin-orbital {} run and cannot be used with [mode] "
            "enable_spin_orbital = true. {}".format(file_name, producer, hint))


def _refuse_so_archive_in_non_so_mode(data, file_name,
                                      check_index_convention=False):
    """Refuse a spin-orbital susceptibility archive in a non-spin-orbital run.

    A spin-orbital archive over N physical orbitals has nd = 2N and so the
    shape of a non-spin-orbital general archive over 2N orbitals; read as the
    latter it would be solved with the wrong index interpretation and no
    error. Spin-orbital RPA/FLEX producers write the scalar bool
    ``spin_orbital = True`` (and a FLEX chiq the tag ``chi_convention =
    "spin_orbital"``); non-spin-orbital archives carry neither, so they load
    unchanged. A marker that is not a scalar bool is refused too (the
    spin-orbital side does the same). Checked on the open npz, before any
    array is read.

    ``check_index_convention`` (Green-function files only) also refuses
    ``index_convention == "spin_block"``: a spin-orbital FLEX green.npz can
    carry that tag alone, whereas the tag sits on EVERY chi0q/chiq archive
    (non-spin-orbital runs included), so the susceptibility readers leave
    it off."""
    hint = ("Enable [mode] enable_spin_orbital = true, or regenerate the "
            "archive with a non-spin-orbital run.")
    if "spin_orbital" in data:
        # a spin-orbital FLEX green.npz carries the marker too; it is named
        # by its payload member (every other archive keeps this wording)
        noun = "Green-function file" if "green" in data else "susceptibility file"
        val = np.asarray(data["spin_orbital"])
        if val.ndim != 0 or val.dtype != np.bool_:
            raise ValueError(
                "{} '{}': malformed spin_orbital marker "
                "({!r}); expected a boolean. {}".format(
                    noun, file_name, data["spin_orbital"], hint))
        if bool(val[()]):
            raise ValueError(
                "{} '{}' was written by a spin-orbital run "
                "(spin_orbital = True) and cannot be used without [mode] "
                "enable_spin_orbital = true: its shape would be misread as a "
                "non-spin-orbital file over twice as many orbitals. {}".format(
                    noun, file_name, hint))
    if ("chi_convention" in data
            and str(data["chi_convention"]) == _SO_FLEX_CHI_CONVENTION):
        raise ValueError(
            "FLEX file '{}' is tagged chi_convention = {!r} (a spin-orbital "
            "FLEX output) and cannot be used without [mode] "
            "enable_spin_orbital = true. {}".format(
                file_name, _SO_FLEX_CHI_CONVENTION, hint))
    from hwave.solver import spin_orbital as _spin_orbital
    if (check_index_convention and "index_convention" in data
            and str(data["index_convention"]) == _spin_orbital.INDEX_CONVENTION):
        raise ValueError(
            "Green-function file '{}' is tagged index_convention = {!r} (a "
            "spin-orbital run's output, spin_orbital = True) and cannot be "
            "used without [mode] enable_spin_orbital = true: its orbital "
            "axis is the spin-block spin-orbital index, not a physical "
            "orbital. {}".format(file_name, _spin_orbital.INDEX_CONVENTION,
                                 hint))


def _refuse_so_green_file(green_path, label="bond_green"):
    """Open a Green-function npz lazily and refuse a spin-orbital one.

    O(1): only the scalar marker keys are read, never the ``green`` array.
    A missing / unreadable file is left to :func:`_load_green_npz`, whose
    error is the more specific one."""
    if not os.path.isfile(green_path):
        return
    try:
        data = np.load(green_path)
    except (OSError, ValueError, EOFError, zipfile.BadZipFile,
            AttributeError):
        return
    try:
        _refuse_so_archive_in_non_so_mode(data, green_path,
                                          check_index_convention=True)
    finally:
        close = getattr(data, "close", None)
        if close is not None:
            close()


def _preflight_refuse_so_archives(input_dict, chi0q_mode):
    """Metadata-only form of :func:`_refuse_so_archive_in_non_so_mode`.

    Run by the public entries before the Hamiltonian is diagonalized or any
    Green function / interaction table is built, so a refusal costs nothing.
    Opens each archive that the non-spin-orbital run will read (the chi0q for
    ``chi0q_mode = "load"``, the FLEX chi_s/chi_c for ``"flex"``) and reads
    only the scalar markers. A missing file is left to the loader, whose
    error is the more specific one. The loader checks stay as defence in
    depth."""
    if chi0q_mode == "load":
        output_info = input_dict["file"]["output"]
        paths = [os.path.join(output_info["path_to_output"],
                              output_info.get("chi0q", "chi0q") + ".npz")]
    elif chi0q_mode == "flex":
        paths = list(_resolve_flex_paths(input_dict)[:2])
    else:
        return
    for path in paths:
        if not os.path.isfile(path):
            continue
        try:
            data = np.load(path)
        except (OSError, ValueError, EOFError, zipfile.BadZipFile,
                AttributeError):
            continue    # unreadable / not an npz: the loader reports it
        if not hasattr(data, "files"):
            continue    # a bare .npy under an .npz name: the loader reports it
        with data:
            _refuse_so_archive_in_non_so_mode(data, path)


def _load_chi0q(input_dict, norb=None, so_run=None):
    """Load chi0q from NPZ file produced by H-wave RPA solver.

    Parameters
    ----------
    input_dict : dict
        Parsed TOML input dictionary.
    so_run : _SORun, optional
        Spin-orbital mode (#88): the archive's provenance is
        checked against this run (:func:`_so_check_archive`, the load-mode
        nmat rule) after the marker checks and before the chi0q array is
        read, and the result is returned as a third item.

    Returns
    -------
    chi0q : ndarray
        Bare susceptibility array.
    static_index : int
        Position of the zero bosonic frequency along the frequency axis
        (determined from the freq_index metadata).
    so_info : SOArchiveInfo
        Only with ``so_run``.
    """
    output_info = input_dict["file"]["output"]
    path_to_output = output_info["path_to_output"]
    chi0q_name = output_info.get("chi0q", "chi0q")
    file_name = os.path.join(path_to_output, chi0q_name + ".npz")

    logger.info("Loading chi0q from {}".format(file_name))
    data = np.load(file_name)
    _reject_ir_native(
        data, file_name,
        "the static Eliashberg solver requires uniform-grid files. Re-run "
        "FLEX with [mode.param] write_densified = true, or switch to "
        "frequency = \"dynamic\" with [eliashberg] matsubara_basis = "
        "\"ir\".")
    enable_spin_orbital = _resolve_spin_orbital_flag(input_dict)
    validate_chi0q_index_convention(data, enable_spin_orbital, file_name)
    so_info = None
    if enable_spin_orbital:
        _require_so_chi0q_marker(data, file_name)
        if so_run is not None:
            so_info = _so_check_archive(data, file_name, so_run, "load", "chi0q",
                                        _SO_RPA_PRODUCER)
    else:
        _refuse_so_archive_in_non_so_mode(data, file_name)
    chi0q = data["chi0q"]
    logger.info("chi0q shape: {}".format(chi0q.shape))

    # coeff_tail provenance check (issue #80): the Matsubara tail correction
    # changes chi0q at O(1) (lambda differed by 50% in the reproduction), so
    # a config whose effective coeff_tail differs from the file's would
    # recompute DIFFERENT physics under chi0q_mode="calc". Warn so load-vs-
    # calc comparisons are not silently inconsistent. Files without the key
    # (older builds) load silently as before.
    if "coeff_tail" in data:
        # type-strict provenance normalization, mirroring
        # rpa._validate_chi0q_provenance (round-8 review)
        import numbers
        _val = np.asarray(data["coeff_tail"])
        _item = _val[()] if _val.ndim == 0 else None
        if (_item is None or isinstance(_item, (bool, np.bool_))
                or not isinstance(_item, numbers.Real)
                or not np.isfinite(float(_item))):
            raise ValueError(
                "chi0q file '{}': malformed coeff_tail ({!r})".format(
                    file_name, data["coeff_tail"]))
        file_tail = float(_item)
        config_tail = float(input_dict.get("mode", {}).get("param", {})
                            .get("coeff_tail", 0.0))
        if file_tail != config_tail:
            logger.warning(
                "chi0q file '{}' was produced with coeff_tail = {} but this "
                "config's effective value is coeff_tail = {}; results are "
                "NOT comparable with a chi0q_mode=\"calc\" run under this "
                "config. Set [mode.param] coeff_tail = {} to match the "
                "file.".format(file_name, file_tail, config_tail, file_tail))
        # Endpoint-convention gate (issue #134): a nonzero-tail chi0q
        # produced before the branch-mean endpoint fix carries a pre-fix
        # O(1/Nmat) error indistinguishable from the array itself; refuse
        # it rather than feed it to the Eliashberg solver. Zero-tail files
        # and legacy files without the coeff_tail key are exempt (same
        # policy as the RPA chi0q_init reader).
        if file_tail != 0.0:
            from hwave.solver.rpa import TAIL_ENDPOINT_CONVENTION
            te = None
            if "tail_endpoint" in data:
                val = np.asarray(data["tail_endpoint"])
                te = str(val[()]) if val.ndim == 0 else repr(val)
            if te is None:
                raise ValueError(
                    "chi0q file '{}' records coeff_tail = {} but no "
                    "tail_endpoint marker: it was produced before the "
                    "equal-time endpoint fix (issue #134) and its "
                    "tail-corrected values carry the pre-fix O(1/Nmat) "
                    "endpoint error. Recompute the bubble with this "
                    "version.".format(file_name, file_tail))
            if te != TAIL_ENDPOINT_CONVENTION:
                raise ValueError(
                    "chi0q file '{}' carries unrecognized tail_endpoint = "
                    "{!r} (this build implements {!r}); refusing to use a "
                    "bubble whose endpoint treatment is unknown.".format(
                        file_name, te, TAIL_ENDPOINT_CONVENTION))

    # Unconditional marker validation (round-10 review): a PRESENT
    # marker must be checked even when the layout/grid cannot be
    # identified (absent CellShape, unknown shape).
    from hwave.solver.rpa import check_momentum_marker
    check_momentum_marker(data, file_name)
    # Fourier-sign provenance gate (issue #133): chi0q is q-labeled; a
    # pre-#133 file carries flipped labels. Legacy files are accepted only
    # when the payload is elementwise q-even (see the validator). The q
    # axes are chosen STRUCTURALLY from the accepted layouts (round-2
    # review: a size search probed the frequency axis when nfreq == nvol,
    # and skipped the already-expanded reference layouts entirely):
    #   4D/6D raw  (nfreq, nvol, ...)            -> axis 1
    #   5D/7D spin-diag (2, nfreq, nvol, ...)    -> axis 2
    #   6D ref  (no, no, Nx, Ny, Nz, nfreq)      -> axes (2, 3, 4)
    #   8D ref  (no, no, no, no, Nx, Ny, Nz, nfreq) -> axes (4, 5, 6)
    from hwave.solver.rpa import validate_momentum_convention
    _cs = input_dict.get("mode", {}).get("param", {}).get("CellShape")
    if _cs is None:
        # no lattice configured (unit-level callers): the grid, and with
        # it the momentum axes, cannot be identified -- production
        # configs always carry CellShape, so no production bypass
        _grid = None
    else:
        _cs = list(_cs)
        while len(_cs) < 3:
            _cs.append(1)
        _grid = tuple(int(x) for x in _cs)
    _nvol = int(np.prod(_grid)) if _grid is not None else -1
    # With the orbital count known (the production entry always passes
    # it), validate the EXACT supported layouts up front (round-6 review:
    # without norb a malformed 8D file passed both marked and unmarked;
    # with norb the (8,8,2,2,2,2)-style raw/ref ambiguity also vanishes,
    # since norb cannot be two values at once).
    _qax = None
    _resolved = False
    if norb is not None and _grid is not None:
        nv, no = _nvol, int(norb)
        # Full trailing shapes for every layout (round-7 review: partial
        # slices accepted malformed spin-diag arrays), and the layout
        # matched HERE directly selects the q axes -- re-running the
        # norb-blind structural routing below raised 'matches BOTH' for
        # shapes this gate had already disambiguated.
        if chi0q.ndim == 4 and tuple(chi0q.shape[1:]) == (nv, no, no):
            _qax = 1
        elif (chi0q.ndim == 6
                and tuple(chi0q.shape[1:]) == (nv, no, no, no, no)):
            _qax = 1
        elif (chi0q.ndim == 6
                and tuple(chi0q.shape[:5]) == (no, no) + _grid):
            _qax = (2, 3, 4)
        elif (chi0q.ndim == 8
                and tuple(chi0q.shape[:7]) == (no,) * 4 + _grid):
            _qax = (4, 5, 6)
        elif (chi0q.ndim == 5 and chi0q.shape[0] == 2
                and tuple(chi0q.shape[2:]) == (nv, no, no)):
            _qax = 2
        elif (chi0q.ndim == 7 and chi0q.shape[0] == 2
                and tuple(chi0q.shape[2:]) == (nv, no, no, no, no)):
            _qax = 2
        else:
            raise ValueError(
                "chi0q file '{}': shape {} matches no supported layout "
                "for norb = {} and CellShape {}; the file is malformed "
                "or from a different system. Regenerate it with the "
                "current version.".format(file_name, chi0q.shape, no,
                                          list(_grid)))
        _resolved = True
    if _resolved or _grid is None:
        pass
    elif chi0q.ndim in (5, 7) and chi0q.shape[0] == 2:
        _qax = 2
    elif chi0q.ndim == 8:
        _qax = (4, 5, 6)
    elif chi0q.ndim == 6:
        # raw (nfreq, nvol, norb x4) vs ref (norb, norb, Nx, Ny, Nz,
        # nfreq): decide by the FULL structure of both patterns
        # (round-3 review: testing shape[1] != nvol misrouted a ref file
        # with norb == nvol onto the orbital axis). Truly ambiguous
        # shapes fail closed unless the marker decides.
        _is_ref = (chi0q.shape[0] == chi0q.shape[1]
                   and tuple(chi0q.shape[2:5]) == _grid)
        _is_raw = (chi0q.shape[1] == _nvol
                   and len(set(chi0q.shape[2:6])) == 1)
        if _is_ref and _is_raw:
            if "momentum_convention" in getattr(data, "files", []):
                _qax = (2, 3, 4)   # any axis: the validator accepts on tag
            else:
                raise ValueError(
                    "chi0q file '{}': shape {} matches BOTH the raw and "
                    "the reference 6D layout and carries no "
                    "momentum_convention marker (issue #133) -- the "
                    "momentum axes cannot be identified safely. "
                    "Regenerate the file with the current version, which "
                    "stamps the marker.".format(file_name, chi0q.shape))
        elif _is_ref:
            _qax = (2, 3, 4)
        elif _is_raw:
            _qax = 1
        else:
            # neither complete pattern matches: fail closed HERE,
            # independent of provenance (round-5 review: routing through
            # the validator let a matching marker return early, and the
            # downstream converter then silently RESHAPED the unknown
            # layout, reinterpreting orbital axes as data -- a marker can
            # establish the Fourier sign, never the array layout)
            raise ValueError(
                "chi0q file '{}': shape {} matches neither the raw "
                "(nfreq, nvol, norb x4) nor the reference "
                "(norb, norb, Nx, Ny, Nz, nfreq) 6D layout for CellShape "
                "{}; the file is malformed or from a different lattice. "
                "Regenerate it with the current version.".format(
                    file_name, chi0q.shape, list(_grid)))
    elif chi0q.ndim == 4:
        _qax = 1
    if _qax is not None:
        validate_momentum_convention(data, file_name, chi0q, _qax, _grid)

    freq_index, file_nmat = _read_freq_meta(data)
    # Identify the frequency axis from the array LAYOUT, never from the
    # freq_index length (a restricted freq_index can coincidentally match
    # an orbital axis).  4D raw: axis 0.  8D ref: last axis.  5D/7D
    # spin-diagonal (leading 2): axis 1.  6D is either
    # raw (nmat, nvol, norb^4; the last four axes are equal) or ref
    # (norb, norb, Nx, Ny, Nz, nmat; the first two axes are equal) --
    # disambiguate structurally, with the freq_index length only as the
    # tiebreaker for degenerate shapes.  The length is computed even for
    # metadata-less files: _static_freq_position now needs it to decide
    # whether the stored axis IS the configured full grid (nfreq ==
    # config_nmat, the only case where centering is unambiguous).
    nfi = np.asarray(freq_index).size if freq_index is not None else None
    if _qax == 1:
        # raw layout (nfreq, nvol, norb...): frequency axis is axis 0
        nfreq = chi0q.shape[0]
    elif _qax == 2:
        # spin-diagonal block-leading layout (2, nfreq, nvol, norb...)
        nfreq = chi0q.shape[1]
    elif _qax in ((2, 3, 4), (4, 5, 6)):
        # reference layout (norb..., Nx, Ny, Nz, nfreq): axis is last
        nfreq = chi0q.shape[-1]
    elif chi0q.ndim in (5, 7) and chi0q.shape[0] == 2:
        # spin-diagonal block-leading layout, identifiable from the array
        # structure alone: the routing above leaves _qax unset when the run
        # carries no CellShape, and the last axis is an ORBITAL axis here,
        # so reading it as the frequency length would feed the full-grid
        # check the wrong number.
        nfreq = chi0q.shape[1]
    elif chi0q.ndim == 4:
        nfreq = chi0q.shape[0]
    elif chi0q.ndim == 8:
        nfreq = chi0q.shape[-1]
    elif chi0q.ndim == 6:
        raw_like = len(set(chi0q.shape[2:])) == 1
        ref_like = chi0q.shape[0] == chi0q.shape[1]
        if raw_like and not ref_like:
            nfreq = chi0q.shape[0]
        elif ref_like and not raw_like:
            nfreq = chi0q.shape[-1]
        elif nfi is not None and nfi == chi0q.shape[0]:
            nfreq = chi0q.shape[0]
        else:
            nfreq = chi0q.shape[-1]
    else:
        # layout not resolved by any of the above: keep the historical
        # last-axis reading, but say so -- for an unsupported shape that
        # axis is a guess, and the static-slice guard downstream judges a
        # length this function could not identify.
        logger.debug(
            "chi0q file '{}': shape {} matches no known layout; taking the "
            "last axis as the frequency axis".format(file_name,
                                                     chi0q.shape))
        nfreq = chi0q.shape[-1]
    config_nmat = input_dict.get("mode", {}).get("param", {}).get("Nmat",
                                                                _DEFAULT_NMAT)
    static_index = _static_freq_position(freq_index, nfreq, config_nmat,
                                         file_name, file_nmat=file_nmat)
    if static_index is not None and static_index != nfreq // 2:
        logger.info("static (zero bosonic frequency) slice at index {} "
                    "of {}".format(static_index, nfreq))
    if so_run is not None:
        return chi0q, static_index, so_info
    return chi0q, static_index


def _auto_chi0q_tensor_scheme(files):
    """hwave_sc's ``chi0q_tensor = "auto"`` rule, behaviour-frozen at its
    1.0.x semantics via ``scheme.CAPABILITIES[...].sc_legacy_forcing``
    (True for CoulombInter/Hund/Ising/Exchange/PairHop/Coulomb; False for
    CoulombIntra/PairLift). Case-insensitive; an unknown key is IGNORED
    (fail-closed discovery is deliberately NOT adopted here -- recorded
    follow-up). Pinned by tests/test_sc_auto_scheme_freeze.py."""
    from hwave.solver.scheme import CAPABILITIES, canonical_name
    for k in (files or {}):
        canon = canonical_name(k)
        if canon is not None and CAPABILITIES[canon].sc_legacy_forcing:
            return "general"
    return "reduced"


def _calc_chi0q_internal(input_dict, chi0q_tensor="auto",
                         precomputed_mu=None, prepared=None,
                         static_only=False, so_context=None):
    """Compute chi0q internally using H-wave's RPA module.

    Instead of loading chi0q from a file, this function creates an RPA solver
    instance and computes chi0q from scratch using the Transfer and interaction
    parameters specified in the TOML config.

    Parameters
    ----------
    input_dict : dict
        Parsed TOML input dictionary.
    chi0q_tensor : str
        Index structure of chi0q. Options:
        - "reduced": 2-index (nmat, nvol, norb, norb). Exact for single-orbital
          and for multi-orbital with CoulombIntra only (no inter-orbital
          interactions).
        - "general": 4-index (nmat, nvol, norb, norb, norb, norb). Required
          when inter-orbital interactions (CoulombInter, Hund, Exchange) are
          present, as their S/C matrices couple different orbital-pair indices.
        - "auto": Use "general" if CoulombInter/Hund/Exchange are present,
          otherwise "reduced".
    precomputed_mu : float, optional
        If provided, skip chemical potential determination and use this value
        directly. This avoids redundant eigenvalue decomposition and mu search
        when the caller has already computed mu (e.g. calc_eliashberg).
    prepared : _MeanFieldBand, optional
        A band already assembled by :func:`_mean_field_band` (issue #273).
        Its RPA solver -- input read once, ``_calc_epsilon_k`` already run --
        is reused, so chi0q is built on exactly the H0(k) the caller used
        for mu and the kernel Green function. ``chi0q_tensor`` is then
        ignored (the solver was constructed with the matching scheme).
    static_only : bool, optional
        If true, compute only the ``nmat // 2`` bosonic slice (the zero
        bosonic frequency the static Eliashberg solve reads) and return it as
        a one-frequency axis (``nmat`` -> 1 below; read it with static index
        0). It equals the full-grid slice to round-off; memory drops from
        ``O(nmat * nvol * nd**4)`` to ``O(nvol * nd**4)`` (issue #272).
        Applies to the ``prepared`` band path too. Dynamic callers must leave
        it false: they need every bosonic frequency.
    so_context : spin_orbital.SOEntryContext, optional
        The entry context of a spin-orbital hwave_sc run: the nested RPA
        reuses it (no second compile of the interaction, no second
        ``orbital_site`` WARNING).

    Returns
    -------
    chi0q : ndarray
        Bare susceptibility array. Shape depends on chi0q_tensor mode:
        - reduced: (nmat, nvol, norb, norb)
        - general: (nmat, nvol, norb, norb, norb, norb)
    """
    if prepared is not None:
        solver = prepared.solver
        logger.info("Computing chi0q internally on the mean-field one-body "
                    "band ({}; calc_scheme={})...".format(
                        ", ".join(prepared.sources), solver.calc_scheme))
        return _chi0q_from_prepared_solver(solver, precomputed_mu,
                                           static_only=static_only)

    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.rpa as sol_rpa

    info_mode = input_dict.get("mode", {})
    info_file = input_dict.get("file", {"input": {}, "output": {}})
    info_inputfile = info_file.get("input", {})
    info_inputfile["path_to_input"] = info_inputfile.get("path_to_input", "")
    info_log = input_dict.get("log", {})
    info_log["print_level"] = info_log.get("print_level", 1)
    info_log["print_step"] = info_log.get("print_step", 1)

    # Determine calc_scheme from chi0q_tensor option
    calc_scheme = _chi0q_tensor_scheme(chi0q_tensor,
                                       info_inputfile.get("interaction", {}))

    info_mode_rpa = dict(info_mode)
    info_mode_rpa["calc_scheme"] = calc_scheme
    # hwave_sc has already reported an undeclared orbital_site (#234)
    info_mode_rpa["_orbital_site_reported"] = True
    # Forward the RESOLVED boolean, not the raw value: RPA mixes truthiness
    # and == True checks on this flag, so a string "false" would diverge
    # internally (round-1 review of #83's guard).
    info_mode_rpa["enable_spin_orbital"] = _resolve_spin_orbital_flag(
        input_dict)
    if so_context is not None:
        from hwave.solver import spin_orbital
        spin_orbital.context_for_nested(info_mode_rpa, so_context)

    logger.info("Computing chi0q internally (calc_scheme={})...".format(calc_scheme))

    # Read input files via QLMSkInput
    read_io = read_input_k.QLMSkInput(info_inputfile)
    ham_info = read_io.get_param("ham")

    # Create RPA solver (sets up Lattice, Interaction, parameters)
    solver = sol_rpa.RPA(ham_info, info_log, info_mode_rpa)
    if info_mode_rpa["enable_spin_orbital"]:
        # the spin-orbital Eliashberg vertex needs the full spinful chi0
        # (nd = 2 * norb spin-orbitals) even when the transfer happens to be
        # spin-diagonal, as the spin-orbital FLEX does (#83)
        solver._keep_spinful = True

    # Get green_info (may contain chi0q_init, trans_mod, green_init)
    green_info = read_io.get_param("green")
    green_info.update(solver.read_init(info_inputfile))

    # Compute chi0q: eigenvalues -> mu -> Green's function -> chi0q
    beta = 1.0 / solver.T

    solver._calc_epsilon_k(green_info)

    if precomputed_mu is not None:
        # Reuse mu from caller to skip redundant bisection
        mu = precomputed_mu
        logger.info("Using precomputed mu = {}".format(mu))
    elif solver.calc_mu:
        if solver.spin_mode == "spin-free":
            Ncond = solver.Ncond / 2
        else:
            Ncond = solver.Ncond
        dist, mu = solver._find_mu(Ncond, solver.T)
    else:
        mu = solver.mu_value

    green0, green0_tail = solver._calc_green(beta, mu)

    chi0q = solver._calc_chi0q(green0, green0_tail, beta,
                               static_only=static_only)
    del green0, green0_tail

    # For spin-free mode, chi0q has shape (1, nmat, nvol, ...)
    # Remove the block index
    if solver.spin_mode in ["spin-free", "spinful"]:
        assert chi0q.shape[0] == 1
        chi0q = chi0q[0]
    else:
        # spin-diag: shape (2, nmat, nvol, ...)
        assert chi0q.shape[0] == 2
        pass

    logger.info("chi0q computed internally, shape: {}".format(chi0q.shape))
    return chi0q


def _chi0q_tensor_scheme(chi0q_tensor, interaction_files):
    """RPA ``calc_scheme`` for an ``[eliashberg] chi0q_tensor`` value:
    ``"auto"`` follows :func:`_auto_chi0q_tensor_scheme`, ``"general"`` is
    general, anything else is reduced (the historical rule)."""
    if chi0q_tensor == "auto":
        return _auto_chi0q_tensor_scheme(interaction_files)
    if chi0q_tensor == "general":
        return "general"
    return "reduced"


def _chi0q_from_prepared_solver(solver, precomputed_mu, static_only=False):
    """chi0q of an RPA solver whose ``_calc_epsilon_k`` has already run.

    Same arithmetic as the tail of :func:`_calc_chi0q_internal`; only the
    spin-free layout reaches here (:func:`_mean_field_band` refuses the
    others), so the block axis is always length 1.
    """
    beta = 1.0 / solver.T
    if precomputed_mu is not None:
        mu = precomputed_mu
        logger.info("Using precomputed mu = {}".format(mu))
    elif solver.calc_mu:
        dist, mu = solver._find_mu(solver.Ncond / 2, solver.T)
    else:
        mu = solver.mu_value
    green0, green0_tail = solver._calc_green(beta, mu)
    chi0q = solver._calc_chi0q(green0, green0_tail, beta,
                               static_only=static_only)
    del green0, green0_tail
    assert chi0q.shape[0] == 1
    chi0q = chi0q[0]
    logger.info("chi0q computed internally, shape: {}".format(chi0q.shape))
    return chi0q


def _read_spin_free_tol(eli_param):
    """``[eliashberg] spin_free_tol`` (case-insensitive): finite float >= 0,
    default :data:`SPIN_FREE_ATOL`."""
    import numbers
    from requests.structures import CaseInsensitiveDict
    raw = CaseInsensitiveDict(eli_param or {}).get("spin_free_tol")
    if raw is None:
        return SPIN_FREE_ATOL
    # type-strict (the _resolve_bond_coeff_tail pattern): a boolean would
    # become 1.0 and a quoted number would be coerced -- both silently
    # change the tolerance, so both are refused
    if isinstance(raw, (bool, np.bool_)) or not isinstance(raw, numbers.Real):
        raise ValueError(
            "[eliashberg] spin_free_tol must be a real number, got "
            "{!r}.".format(raw))
    tol = float(raw)
    if not np.isfinite(tol) or tol < 0.0:
        raise ValueError(
            "[eliashberg] spin_free_tol must be a finite number >= 0; got "
            "{!r}.".format(raw))
    return float(tol)


def _load_band_into_solver(solver, up):
    """Make ``solver`` diagonalize exactly ``up`` (the checked band) for
    chi0q: it is handed to RPA's own ``_calc_epsilon_k`` as an explicitly
    spin-free trans_mod (both spin blocks equal to ``up``, no mixing), with
    the Extern term switched off for that call because it is already inside
    ``up``. RPA's H0_k must then equal ``up`` exactly."""
    nvol, norb, _ = up.shape
    tm = np.zeros((nvol, 2 * norb, 2 * norb), dtype=np.complex128)
    tm[:, :norb, :norb] = up
    tm[:, norb:, norb:] = up
    ham = solver.ham_info
    saved_extern = ham.ham_extern_q
    ham.ham_extern_q = None
    try:
        solver._calc_epsilon_k({"trans_mod": tm})
    finally:
        ham.ham_extern_q = saved_extern
    H0 = np.asarray(solver.H0_k)
    if (solver.spin_mode != "spin-free" or H0.shape != (1,) + up.shape
            or not np.array_equal(H0[0], up)):
        raise RuntimeError(
            "hwave_sc: the RPA solver did not reproduce the checked "
            "one-body band (spin_mode={}, shape {}); refusing to compute "
            "chi0q on a different band.".format(solver.spin_mode, H0.shape))


def _mean_field_band(input_dict, sources, calc_scheme, spatial_shape, norb,
                     spin_free_tol=SPIN_FREE_ATOL, allow_average=True,
                     allow_hermitize=True):
    """Assemble the one-body band through the RPA machinery (issue #273).

    The band is built exactly as the internal chi0q builds it --
    ``QLMSkInput`` + ``RPA.read_init`` + ``RPA._calc_epsilon_k`` -- so that
    trans_mod (T + Sigma_HF written by UHFk), T + Sigma_HF[green_init] (UHFk's
    Hartree + Fock mean field of G, issue #278) and the Extern field define
    the SAME H0(k) for chi0q, the chemical potential
    and the kernel Green function.

    Only a spin-independent band is accepted: the singlet/triplet kernel
    carries no spin index, so a spin-diagonal or spin-mixing mean field
    cannot be represented. The band is assembled ONCE by RPA
    (:func:`_checked_band`: finiteness, then the absolute tolerance
    :data:`SPIN_FREE_ATOL`), after the input files are read by the RPA
    readers and before H0 is diagonalized -- so before any chi0q, pair
    bubble or kernel is built, and before a chi0q file is loaded. The
    checked band (its up-up block, Extern included) is the band used for
    epsilon(k), mu, the kernel Green function and -- through
    :func:`_load_band_into_solver`, which verifies RPA reproduces it
    exactly -- chi0q.
    """
    from requests.structures import CaseInsensitiveDict
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.rpa as sol_rpa

    fin = (input_dict.get("file", {}) or {}).get("input", {}) or {}
    for key in _ONE_BODY_FILE_KEYS:
        spellings = [k for k in fin if str(k).lower() == key]
        if len(spellings) > 1:
            raise ValueError(
                "[file.input] {} is given more than once with different "
                "spellings ({}); keys are case-insensitive.".format(
                    key, ", ".join(repr(k) for k in spellings)))
    ci = CaseInsensitiveDict(fin)
    # Only the keys the band needs: chi0q_init / UHF initial files are not
    # read here (an RPA read_init would otherwise load a chi0q_init array).
    # The lower-case spelling is what RPA.read_init looks up.
    info_inputfile = {"path_to_input": ci.get("path_to_input", ""),
                      "interaction": ci.get("interaction", {})}
    for key in _ONE_BODY_FILE_KEYS:
        if key in ci:
            info_inputfile[key] = ci[key]

    info_log = dict(input_dict.get("log", {}) or {})
    info_log.setdefault("print_level", 1)
    info_log.setdefault("print_step", 1)
    import copy
    # deep copy: the RPA constructor normalizes [mode.param] in place, and
    # the caller's dict has to stay as the user wrote it
    info_mode_rpa = copy.deepcopy(dict(input_dict.get("mode", {}) or {}))
    info_mode_rpa["calc_scheme"] = calc_scheme
    info_mode_rpa["_orbital_site_reported"] = True
    info_mode_rpa["enable_spin_orbital"] = _resolve_spin_orbital_flag(
        input_dict)

    logger.info("One-body band from {} (RPA H0 assembly; used for chi0q, "
                "mu and the kernel Green function)".format(", ".join(sources)))
    read_io = read_input_k.QLMSkInput(info_inputfile)
    solver = sol_rpa.RPA(read_io.get_param("ham"), info_log, info_mode_rpa)
    green_info = read_io.get_param("green")
    green_info.update(solver.read_init(info_inputfile))
    if "green_init" in green_info and "trans_mod" not in green_info:
        # G-level refusal, before the band is assembled
        _check_green_init_spin_free(green_info["green_init"], solver.norb,
                                    spin_free_tol)

    # Assemble the band ONCE, check it (spin independence, finiteness),
    # before any diagonalization and before any chi0q / pair-bubble /
    # kernel work; the checked band is the band that is used.
    up = _checked_band(solver, green_info, sources, spin_free_tol,
                       allow_average=allow_average,
                       allow_hermitize=allow_hermitize)

    Nx, Ny, Nz = spatial_shape
    if up.shape != (Nx * Ny * Nz, norb, norb):
        raise ValueError(
            "hwave_sc: the one-body band assembled from {} has shape {}, "
            "expected {} for the {}x{}x{} grid with norb={}.".format(
                " / ".join(sources), up.shape,
                (Nx * Ny * Nz, norb, norb), Nx, Ny, Nz, norb))
    # chi0q is computed by this solver on exactly the checked band
    _load_band_into_solver(solver, up)
    H0 = up[np.newaxis]
    # RPA (1, nvol, a, b) with nvol in C order over (Nx, Ny, Nz) -- the same
    # flattening as _green_sc_to_canonical -- to sc.py's (a, b, Nx, Ny, Nz).
    epsilon_k = np.ascontiguousarray(
        H0[0].reshape(Nx, Ny, Nz, norb, norb).transpose(3, 4, 0, 1, 2),
        dtype=np.complex128)
    bare_q = getattr(solver.ham_info, "ham_trans_q", None)
    same_as_bare = (bare_q is not None
                    and np.shape(bare_q) == H0[0].shape
                    and np.array_equal(H0[0], bare_q))
    return _MeanFieldBand(epsilon_k=epsilon_k, solver=solver,
                          sources=tuple(sources), same_as_bare=same_as_bare)


def _read_interaction_files(input_dict):
    """Read Transfer and interaction files using wan90 module.

    Parameters
    ----------
    input_dict : dict
        Parsed TOML input dictionary.

    Returns
    -------
    geom_info : dict
        Geometry information (norb, rvec, center).
    hr : dict
        Transfer (hopping) parameters.
    interactions : dict
        Dictionary of interaction parameters keyed by type name.
    """
    # Case-insensitive on purpose (round-4 review): the FLEX reader stores
    # these keys in a CaseInsensitiveDict, so a run configured with e.g.
    # 'coulombinter' PRODUCES a susceptibility with that interaction while
    # an exact-case lookup here would silently read an EMPTY interaction
    # set -- the compatibility gate and the pairing vertex would both be
    # built from nothing.
    from requests.structures import CaseInsensitiveDict
    files = CaseInsensitiveDict(input_dict["file"]["input"]["interaction"])
    path_to_input = files.get("path_to_input", "")

    geom_file = os.path.join(path_to_input, files["Geometry"])
    geom_info = wan90.read_geom(geom_file)
    logger.info("norb = {}".format(geom_info["norb"]))

    transfer_file = os.path.join(path_to_input, files["Transfer"])
    hr = wan90.read_w90(transfer_file)

    interaction_types = ["CoulombIntra", "CoulombInter", "Hund", "Exchange",
                        "Ising", "PairLift", "PairHop"]
    interactions = {}
    from hwave.solver.declarations import (validate_hermitian_closure,
                                           validate_no_onsite_same_orbital_inter)

    for itype in interaction_types:
        if itype in files:
            f = os.path.join(path_to_input, files[itype])
            logger.info("Reading {} from {}".format(itype, f))
            tbl = wan90.read_w90(f)
            # issue #93: fail fast on non-Hermitian-closed declarations,
            # with the same rule and tolerance as the k-space reader
            validate_hermitian_closure(itype, tbl, source=f)
            validate_no_onsite_same_orbital_inter(itype, tbl, source=f)
            interactions[itype] = tbl

    # combined 'Coulomb' input: same decomposition as UHFk/RPA
    # (wan90.split_coulomb; r=0 diagonal -> CoulombIntra, rest -> CoulombInter)
    if "Coulomb" in files:
        if "CoulombIntra" in interactions or "CoulombInter" in interactions:
            raise ValueError(
                "Coulomb cannot be specified together with "
                "CoulombIntra or CoulombInter")
        f = os.path.join(path_to_input, files["Coulomb"])
        logger.info("Reading Coulomb from {}".format(f))
        _coulomb_tbl = wan90.read_w90(f)
        validate_hermitian_closure("Coulomb", _coulomb_tbl, source=f)
        coulomb_intra, coulomb_inter = wan90.split_coulomb(_coulomb_tbl)
        if coulomb_intra:
            interactions["CoulombIntra"] = coulomb_intra
        if coulomb_inter:
            interactions["CoulombInter"] = coulomb_inter

    return geom_info, hr, interactions


# ---------------------------------------------------------------------------
# k-space construction
# ---------------------------------------------------------------------------


def _read_so_extern(input_dict, norb_phys):
    """The spin-orbital ``Extern`` table, read and checked (issue #88).

    In spin-orbital mode the Extern file holds ``H_ab(R)`` in PHYSICAL
    orbital indices ``1..norb_phys``; the operator is ``coeff_extern *
    sum H_ab(R) sigma^z_{ss'} c+_{as} c_{bs'}``. The header (line 2) must
    declare ``norb_phys`` and every index must be in range
    (:func:`hwave.solver.spin_orbital.check_extern_indices`); a malformed
    header is refused by :func:`hwave.qlmsio.wan90.read_w90_header_dim`.
    Run in every ``chi0q_mode`` (calc, load, flex), before any heavy work.

    Returns
    -------
    (extern_table, file_name) : tuple
        ``(None, None)`` when ``[file.input.interaction]`` names no Extern
        file (the key is matched case-insensitively).
    """
    from requests.structures import CaseInsensitiveDict
    from hwave.solver import spin_orbital
    files = CaseInsensitiveDict(input_dict["file"]["input"]["interaction"])
    if "Extern" not in files:
        return None, None
    file_name = os.path.join(files.get("path_to_input", ""), files["Extern"])
    dim = wan90.read_w90_header_dim(file_name)
    table = wan90.read_w90(file_name)
    spin_orbital.check_extern_indices(table, norb_phys, file_name,
                                      declared_norb=dim)
    logger.info("Reading Extern from {} ({} physical orbitals)".format(
        file_name, norb_phys))
    return table, file_name


def _so_coeff_extern(input_dict):
    """``[mode.param] coeff_extern`` (case-insensitive, default 0.0): a
    finite real number."""
    from requests.structures import CaseInsensitiveDict
    from hwave.solver import spin_orbital
    raw = CaseInsensitiveDict(input_dict["mode"].get("param", {}) or {}).get(
        "coeff_extern", 0.0)
    return spin_orbital.validate_coeff_extern(raw)


class _SORun(NamedTuple):
    """This spin-orbital hwave_sc run as the archive checks see it (#88):
    its band (:func:`hwave.solver.fingerprint.band_inputs`), its state
    (``beta``, ``nmat``, ``calc_mu = True``, ``target`` in electrons per cell;
    ``mu``, ``density`` and ``mu_tol`` are NaN, in memory only, never
    compared and never serialized) and the entry context
    (:class:`hwave.solver.spin_orbital.SOEntryContext`: the interaction
    compiled once, the ``orbital_site`` WARNING state) that the nested RPA
    runs reuse."""
    band: object
    state: object
    context: object = None


def _so_run_inputs(input_dict, hr, extern_table, norb_phys, beta, nmat,
                   n_filling, context=None):
    """The :class:`_SORun` of a spin-orbital hwave_sc run, built from its own
    input before any archive is read. hwave_sc always searches mu
    (``filling`` is mandatory), so the target is ``filling * nd`` electrons
    per cell over all ``nd = 2 * norb_phys`` spin-orbitals. The band is the
    Transfer table as read (interleaved), the physical-orbital Extern table
    (or None), and the metadata of the spin-orbital RPA/FLEX producers
    (``coeff_extern``, the index convention, the orbital counts, CellShape
    and the resolved SubShape). ``context`` is stored as
    :attr:`_SORun.context`."""
    from requests.structures import CaseInsensitiveDict
    from hwave.solver import fingerprint, spin_orbital
    # case-insensitive, like every [mode.param] key of the producers
    param = CaseInsensitiveDict(input_dict["mode"].get("param", {}) or {})
    cs = param["CellShape"]
    cs = list(cs) if isinstance(cs, (list, tuple)) else [cs]
    cs = cs + [1] * (3 - len(cs))
    ss = param.get("SubShape", cs)
    ss = list(ss) if isinstance(ss, (list, tuple)) else [ss]
    ss = ss + [1] * (3 - len(ss))
    meta = fingerprint.BandMeta(
        coeff_extern=_so_coeff_extern(input_dict), enable_spin_orbital=True,
        index_convention=spin_orbital.INDEX_CONVENTION,
        norb=2 * int(norb_phys), norb_phys=int(norb_phys),
        cell_shape=tuple(int(x) for x in cs), sub_shape=tuple(int(x) for x in ss))
    nan = float("nan")
    state = fingerprint.SOState(float(beta), int(nmat), nan, nan, True,
                                float(n_filling) * 2 * int(norb_phys), nan)
    return _SORun(fingerprint.band_inputs(hr, extern_table, meta), state, context)


#: the producers a spin-orbital hwave_sc run names when an archive it reads
#: must be regenerated: the RPA run of a chi0q load, the FLEX run of a flex
#: entry (static and dynamic)
_SO_RPA_PRODUCER = "the spin-orbital RPA run (hwave RPA, enable_spin_orbital = true)"
_SO_FLEX_PRODUCER = "the spin-orbital FLEX run (hwave FLEX, enable_spin_orbital = true)"


def _so_check_archive(data, path, so_run, nmat_rule, kind, producer, **kw):
    """:func:`hwave.solver.fingerprint.check_so_archive` of the ``kind``
    archive (``"chi0q"``, ``"chiq"``, ``"green"`` or ``"sigma"``) against
    ``so_run``: its band, state and the interaction fingerprints of its
    entry context; this run consumes on-site (``"so_onsite"``) artifacts.
    ``producer`` (:data:`_SO_RPA_PRODUCER` or :data:`_SO_FLEX_PRODUCER`)
    is the run that writes the archive."""
    from hwave.solver import fingerprint
    return fingerprint.check_so_archive(
        data, path, kind=kind, run_band=so_run.band, run_state=so_run.state,
        nmat_rule=nmat_rule, run_interaction=so_run.context.interaction,
        need="so_onsite", expected_producer=producer, **kw)


def _so_own_state(eigenvalues, beta, mu, nmat, target):
    """The state of a kernel built on this run's own band and mu: the masked
    Fermi count of the eigenvalues at ``mu`` (``rpa._masked_fermi_delta_n``,
    the count ``_determine_mu`` roots) per cell, the ``calc_mu = True``
    branch."""
    from hwave.solver import fingerprint
    from hwave.solver import rpa as _rpa
    w = np.asarray(backend.to_host(eigenvalues))
    nvol = int(np.prod(w.shape[:3]))
    count = _rpa._masked_fermi_delta_n(w, 1.0 / beta, mu, 0.0, 100.0)[0]
    return fingerprint.so_state(beta=beta, nmat=nmat, mu=mu,
                                density=float(count) / nvol, calc_mu=True,
                                target=target)


def _so_flex_consumed_state(green_sc, epsilon_k, sigma_static, mu, beta,
                            target):
    """The state of a borrowed spin-orbital FLEX Green function as consumed:
    the closure density FLEX itself uses
    (``hartree_fock.equal_time_density`` with the reference
    ``heff_eigenpairs(H0_k, sigma_static)``, flex.py
    ``_phase_b_density_and_hf``) per cell, at the mu the kernel applied.
    ``green_sc`` is ``(nd, nd, Nx, Ny, Nz, nmat)``, ``epsilon_k`` this run's
    spin-block band ``(nd, nd, Nx, Ny, Nz)`` and ``sigma_static`` the FLEX
    run's static self-energy ``(nvol, nd, nd)``."""
    from hwave.solver import fingerprint, hartree_fock as _hf
    eps = np.asarray(backend.to_host(epsilon_k))
    nd = eps.shape[0]
    shape = tuple(int(x) for x in eps.shape[2:5])
    nvol = int(np.prod(shape))
    H0 = np.moveaxis(eps.reshape(nd, nd, nvol), -1, 0)
    g = _green_sc_to_canonical(np.asarray(backend.to_host(green_sc)))
    dens = _hf.equal_time_density(g, _hf.heff_eigenpairs(H0, sigma_static), mu,
                                  beta, shape, sym_tol=None)
    return fingerprint.so_state(beta=beta, nmat=int(g.shape[1]), mu=mu,
                                density=float(np.real(dens.n_per_spin)) / nvol,
                                calc_mu=True, target=target)


def _so_writer_members(stage, label, so_run, fp_green, state):
    """The fingerprint members of a hwave_sc gap archive: this run's band
    and interaction fingerprints, and the Green function and state its
    kernel used (never the producer's); ``fp_green`` is
    :func:`_so_green_fp` of that Green function."""
    from hwave.solver import fingerprint
    return fingerprint.so_provenance_members(
        producer_stage=stage, producer_label=label,
        transfer=so_run.band.transfer, extern=so_run.band.extern,
        meta=so_run.band.meta, fp_green=fp_green, state=state,
        interaction=so_run.context.interaction)


def _so_green_fp(green_kw):
    """``fp_green`` of a Green function ``(nd, nd, Nx, Ny, Nz, nmat)`` as a
    hwave_sc gap archive records it."""
    from hwave.solver import fingerprint
    return fingerprint.green_digest(green_kw, "sc_o_o_xyz_f")


def _so_epsilon_k(input_dict, kx_array, ky_array, kz_array, hr, extern_table,
                  norb_phys):
    """The spin-orbital one-body band ``(nd, nd, Nx, Ny, Nz)`` in spin-block
    order (``nd = 2 * norb_phys``), as the spin-orbital RPA assembles it::

        interleaved_to_spinblock(T(k)) + coeff_extern * H1(k) (x) sigma^z

    ``T(k)`` from the interleaved Transfer table ``hr``, ``H1(k)`` from the
    physical-orbital Extern table (``None``: no field). The Zeeman sign is
    that of the Extern definition: ``+h`` shifts spin up by ``+h`` and spin
    down by ``-h``. A Zeeman field written into the Transfer file as well
    (Route 1) adds to it, with a WARNING."""
    from hwave.solver import spin_orbital
    nd = 2 * norb_phys
    epsilon_k = spin_orbital.interleaved_to_spinblock(
        _build_hamiltonian_k(kx_array, ky_array, kz_array, hr, nd),
        norb_phys, axes=(0, 1))
    if extern_table is None:
        return epsilon_k
    coeff = _so_coeff_extern(input_dict)
    if coeff == 0.0:
        return epsilon_k
    up, dn = 0, 1
    for a in range(norb_phys):
        t_up = hr.get(((0, 0, 0), (spin_orbital.interleaved_index(a, up),
                                   spin_orbital.interleaved_index(a, up))), 0.0)
        t_dn = hr.get(((0, 0, 0), (spin_orbital.interleaved_index(a, dn),
                                   spin_orbital.interleaved_index(a, dn))), 0.0)
        if t_up != t_dn:
            logger.warning(
                "Extern adds to spin-dependent on-site terms already present "
                "in the Transfer file (orbital {}); the two fields add".format(
                    a + 1))
            break
    h1 = coeff * _build_hamiltonian_k(kx_array, ky_array, kz_array,
                                      extern_table, norb_phys)
    if np.any(h1 != 0):
        logger.info("Zeeman sign: +h raises spin-up (favours spin-down "
                    "occupation)")
    epsilon_k[:norb_phys, :norb_phys] += h1
    epsilon_k[norb_phys:, norb_phys:] -= h1
    return epsilon_k


# ---------------------------------------------------------------------------
# Green's function
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# RPA vertices
# ---------------------------------------------------------------------------


def _resolve_flex_dir(input_dict):
    """The directory the FLEX outputs are read from.

    ``[file.input] path_to_flex_output``, falling back to ``[file.output]
    path_to_output`` and finally to ``"output"``. Factored out so every
    consumer of the FLEX outputs (the chi/green loaders here and the bond
    archive of the dynamic bond path) resolves the SAME directory.
    """
    file_input = input_dict.get("file", {}).get("input", {})
    return file_input.get("path_to_flex_output",
                          input_dict.get("file", {}).get("output", {}).get(
                              "path_to_output", "output"))


def _resolve_flex_paths(input_dict):
    """Resolve the FLEX output directory and chi_s/chi_c/green file paths.

    Shared by `_load_flex_susceptibilities_full` (which reads the files) and
    `_load_flex_susceptibilities` (which re-derives the static frequency
    index from the same files' metadata), so the two never disagree about
    which files are in play.
    """
    flex_dir = _resolve_flex_dir(input_dict)

    eli_param = input_dict.get("eliashberg", {})
    chi_s_file = eli_param.get("flex_chi_s", "chiq_s.npz")
    chi_c_file = eli_param.get("flex_chi_c", "chiq_c.npz")
    green_file = eli_param.get("flex_green", "green.npz")

    chi_s_path = os.path.join(flex_dir, chi_s_file)
    chi_c_path = os.path.join(flex_dir, chi_c_file)
    green_path = os.path.join(flex_dir, green_file)
    return chi_s_path, chi_c_path, green_path


def _load_flex_susceptibilities_full(input_dict, norb, Nx, Ny, Nz,
                                     allow_ir=False, interactions=None):
    """Load FLEX-computed susceptibilities from NPZ files (full frequency axis).

    Parameters
    ----------
    input_dict : dict
        Parsed TOML configuration.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        Grid dimensions.

    Returns
    -------
    chis_w : ndarray
        Spin susceptibility, shape (Nx, Ny, Nz, nd, nd, nmat) with nd = norb^2
        (after the spin-orbital block-extraction/expansion, if needed) and
        the full bosonic Matsubara axis as the last dimension.
    chic_w : ndarray
        Charge susceptibility, same shape/convention as `chis_w`.
    green_w : ndarray or None
        Dressed Green's function if available, shape
        (norb, norb, Nx, Ny, Nz, nmat) with the full fermionic axis.
    chi_convention : str
        Provenance/S-C-convention tag of chis_w/chic_w: "myo" (general
        full-vertex FLEX) or "kuroki" (reduced FLEX / legacy files). It selects
        which S/C interaction matrices _compute_vertices_flex pairs the
        susceptibilities with (one implementation since the #113
        adjudication; the tag is a layout/provenance discriminator)
        and which orbital-pair shape family the file uses (orbital-pair
        nd=norb^2 for "myo", spin-orbital nd=norb*ns for "kuroki"). It is NOT an
        index-layout flag: chis_w/chic_w are returned in the public RPA
        [a,c,b,d] orbital-pair order for both (issue #78).
    ir_meta : dict or None
        Only when ``allow_ir=True``: ``None`` for uniform files, else the
        per-file IR node metadata ``{"chis":..., "chic":..., "green":...}``
        ("green" absent when there is no green file). Mixed encodings raise.
    """
    if not allow_ir:
        (chi_s_raw, chi_c_raw, chi_convention,
         legacy_tags) = _read_flex_chi_raw(
            input_dict, interactions=interactions)
        ir_meta = None
        green_w = _load_flex_green(input_dict, norb, Nx, Ny, Nz,
                                   interactions=interactions)
    else:
        (chi_s_raw, chi_c_raw, chi_convention, ir_meta,
         legacy_tags) = _read_flex_chi_raw(
            input_dict, allow_ir=True, interactions=interactions)
        green_w, green_meta = _load_flex_green(input_dict, norb, Nx, Ny, Nz,
                                               allow_ir=True,
                                               interactions=interactions)
        if green_w is not None and (green_meta is None) != (ir_meta is None):
            raise ValueError(
                "all dynamic-Eliashberg inputs must share one encoding; "
                "re-run FLEX so chiq_s/chiq_c/green are all densified or "
                "all IR-native (chi: {}, green: {}).".format(
                    "IR-native" if ir_meta is not None else "densified",
                    "IR-native" if green_meta is not None else "densified"))
        if ir_meta is not None and green_meta is not None:
            ir_meta = dict(ir_meta, green=green_meta)

    # Expand the FULL frequency axis (the static slice is selected by the
    # caller). The frequency axis is moved from leading to trailing position.
    legacy_s, legacy_c = legacy_tags
    if _accept_up_block_only(input_dict):
        legacy_s = legacy_c = "config_override"
    _check_spin_block_discarded(chi_s_raw, norb, chi_convention, "chi_s",
                                legacy_s)
    _check_spin_block_discarded(chi_c_raw, norb, chi_convention, "chi_c",
                                legacy_c)
    chis_w = np.moveaxis(
        _expand_flex_chi(chi_s_raw, norb, Nx, Ny, Nz, chi_convention), 0, -1)
    chic_w = np.moveaxis(
        _expand_flex_chi(chi_c_raw, norb, Nx, Ny, Nz, chi_convention), 0, -1)

    if not allow_ir:
        return chis_w, chic_w, green_w, chi_convention
    return chis_w, chic_w, green_w, chi_convention, ir_meta


_STATIC_IR_HINT = (
    "the static Eliashberg solver requires uniform-grid files. Re-run FLEX "
    "with [mode.param] write_densified = true, or switch to frequency = "
    "\"dynamic\" with [eliashberg] matsubara_basis = \"ir\".")


def _accept_up_block_only(input_dict):
    """Whether ``[eliashberg] accept_up_block_only`` asserts the layout.

    The file tag is the better signal, but files written before this check
    existed cannot carry one -- so without a configuration route the escape
    hatch would be unreachable for exactly the files that need it. The user
    takes responsibility here rather than the loader guessing from values it
    cannot disambiguate.
    """
    eli = input_dict.get("eliashberg", {})
    value = eli.get("accept_up_block_only", False)
    if isinstance(value, bool):
        return value
    # Not bool(value): a TOML typo like the string "false" is truthy in Python,
    # and silently enabling an override that relaxes a correctness guard is the
    # worst possible reading of a malformed value.
    raise ValueError(
        "[eliashberg] accept_up_block_only must be a boolean (true/false), "
        "got {!r}.".format(value))


def _legacy_up_block_only(data):
    """Whether an npz declares the legacy single-populated-block layout.

    Opt-in only: the values cannot distinguish "the other blocks were never
    filled" from "this run is fully polarized", so the file has to say so.
    """
    if "chi_spin_blocks" not in data:
        return False
    return str(data["chi_spin_blocks"]).strip().lower() == "up_only"


def _onsite_transpose_asymmetry(table):
    """Max |X[0,(a,b)] - X[0,(b,a)]| over the on-site entries of a raw
    declaration table. Off-site entries are ignored on purpose: their
    same-operator partner lives at -R (and possibly under a wrapped
    canonical key), and the general FLEX path rejects off-site two-body
    terms anyway -- the reachable orientation-sensitive case (#101) is
    an asymmetric ON-SITE inter-orbital entry."""
    worst = 0.0
    for (irvec, orbvec), v in table.items():
        if tuple(irvec) != (0, 0, 0):
            continue
        a, b = orbvec
        if a == b:
            continue
        partner = table.get((tuple(irvec), (b, a)), 0.0)
        diff = abs(complex(v) - complex(partner))
        if not np.isfinite(diff):
            # fail CLOSED: max() drops NaN (nan > x is False), so a
            # non-finite declaration would otherwise read as symmetric
            return float("inf")
        worst = max(worst, diff)
    return worst


def _onsite_hermitian_pair_mismatch(table):
    """Max |X[0,(a,b)] - conj(X[0,(b,a)])| over on-site entries.

    PairHop's two declarations of one pair are HERMITIAN partners, so its
    orientation never depended on the #99 transpose -- but the version-2
    stamp introduced the conjugated-mean symmetrisation of those
    declarations, and a pair that is NOT Hermitian-closed was read
    differently by older builds (round-4 review measured |dS| = 0.15 in
    the myo matrices for a real asymmetric declaration). Fails closed on
    non-finite values like the transpose helper."""
    worst = 0.0
    for (irvec, orbvec), v in table.items():
        if tuple(irvec) != (0, 0, 0):
            continue
        a, b = orbvec
        if a == b:
            continue
        partner = table.get((tuple(irvec), (b, a)), 0.0)
        diff = abs(complex(v) - np.conj(complex(partner)))
        if not np.isfinite(diff):
            return float("inf")
        worst = max(worst, diff)
    return worst


def _read_flex_chi_raw(input_dict, allow_ir=False, interactions=None):
    """Read the raw FLEX chi_s / chi_c NPZ arrays and their orbital convention.

    Returns ``(chi_s_raw, chi_c_raw, chi_convention, legacy_tags)`` in the
    H-wave layout. ``legacy_tags`` is the pair of ``chi_spin_blocks`` flags read
    from the SAME open handles as the arrays: reopening the files afterwards
    would let a tagged replacement authorize the untagged -- possibly polarized
    -- array already in memory. No values from the replacement would be used,
    but the wrong ones would be accepted in silence.
    ``(nfreq, nvol, nd, nd)`` -- no reshape/expansion, so callers that only
    need one static frequency can slice before expanding.

    The tuple is ``(chi_s_raw, chi_c_raw, chi_convention, legacy_tags)``; with
    ``allow_ir=True`` it is ``(chi_s_raw, chi_c_raw, chi_convention, ir_meta,
    legacy_tags)``. In that case the return value
    gains a fourth element ``ir_meta``: ``None`` for uniform files, or
    ``{"chis": meta, "chic": meta}`` for IR-native ones (each file carries
    its own node set; the caller refits both independently). The arity
    changes ONLY under ``allow_ir=True``, so every existing call site keeps
    its unpacking. Mixed encodings (one native, one densified) are rejected
    -- refitting one of the pair would silently pair susceptibilities of
    different provenance."""
    chi_s_path, chi_c_path, _ = _resolve_flex_paths(input_dict)

    logger.info("Loading FLEX chi_s from: {}".format(chi_s_path))
    data_s = np.load(chi_s_path)
    _non_so = not _resolve_spin_orbital_flag(input_dict)
    if _non_so:
        _refuse_so_archive_in_non_so_mode(data_s, chi_s_path)
    if not allow_ir:
        _reject_ir_native(data_s, chi_s_path, _STATIC_IR_HINT)
    chi_s_raw = data_s["chiq_s"] if "chiq_s" in data_s else data_s["chiq"]
    from hwave.solver.rpa import (validate_momentum_convention,
                                  check_momentum_marker)
    # unconditional (round-10): a present marker is validated even when
    # the layout-dependent evenness gate below cannot run
    check_momentum_marker(data_s, chi_s_path)
    # Fourier-sign provenance (issue #133): the H-wave chi layout --
    # uniform AND IR-native alike -- has the flattened q volume on axis 1
    # (round-4 review: IR-native files were wrongly exempted; no IR
    # consumer gates them).
    _cs = list(input_dict.get("mode", {}).get("param", {}).get(
        "CellShape", [1, 1, 1]))
    while len(_cs) < 3:
        _cs.append(1)
    _grid = tuple(int(x) for x in _cs)
    if chi_s_raw.ndim >= 2 and chi_s_raw.shape[1] == int(np.prod(_grid)):
        validate_momentum_convention(data_s, chi_s_path, chi_s_raw, 1,
                                     _grid)
    # Orbital convention tag (general FLEX writes "myo", reduced "kuroki").
    # The tag and the general-path MYO consumption ship together, so any
    # untagged file from a released build is necessarily a reduced/Kuroki
    # output -> default to "kuroki". (A pre-tag general output could only exist
    # as a transient artifact of an unreleased dev build; the s/c-agreement
    # check below still guards against accidentally mixing conventions.)
    # This default is legacy-only: an IR-native file missing the tag is
    # rejected below instead of silently defaulting (an IR-native norb=2 MYO
    # file would otherwise be mis-read as spin-orbital "kuroki").
    chi_convention_present_s = "chi_convention" in data_s
    chi_convention = (str(data_s["chi_convention"])
                      if chi_convention_present_s else "kuroki")

    logger.info("Loading FLEX chi_c from: {}".format(chi_c_path))
    data_c = np.load(chi_c_path)
    if _non_so:
        _refuse_so_archive_in_non_so_mode(data_c, chi_c_path)
    if not allow_ir:
        _reject_ir_native(data_c, chi_c_path, _STATIC_IR_HINT)
    chi_c_raw = data_c["chiq_c"] if "chiq_c" in data_c else data_c["chiq"]
    check_momentum_marker(data_c, chi_c_path)
    if chi_c_raw.ndim >= 2 and chi_c_raw.shape[1] == int(np.prod(_grid)):
        validate_momentum_convention(data_c, chi_c_path, chi_c_raw, 1,
                                     _grid)
    # The spin and charge files must share one convention; combining e.g. an MYO
    # chi_s with a Kuroki chi_c would build a meaningless pairing vertex.
    chi_convention_present_c = "chi_convention" in data_c
    chi_convention_c = (str(data_c["chi_convention"])
                        if chi_convention_present_c else "kuroki")
    if chi_convention_c != chi_convention:
        raise ValueError(
            "FLEX chi_s and chi_c have different conventions ('{}' vs '{}'); "
            "they must come from the same run. Check flex_chi_s/flex_chi_c.".format(
                chi_convention, chi_convention_c))

    # S/C vertex-content versioning (#113). The Hund / Exchange / Ising vertex
    # entries were corrected against exact diagonalization; susceptibilities
    # computed with the OLD entries must not be paired with the corrected
    # matrices when the interaction set contains an affected type -- the
    # kernel would silently mix two different interactions. U/V-only inputs
    # are provably unchanged, so legacy files stay usable there.
    if interactions is not None:
        # activation requires actual CONTENT, not key presence: an explicitly
        # configured but empty Hund/Exchange/Ising file contributes nothing
        affected = [t for t in ("Hund", "Exchange", "Ising")
                    if interactions.get(t)]
        # #101: a SECOND independent reason a no-version file is unusable.
        # PR #99 changed the orbital orientation _build_interaction_k stores
        # (which the general/myo S/C matrices consume), and version-2 files
        # are necessarily post-#99 (the version tag came later), so for a
        # myo file the version requirement doubles as the orientation
        # marker. Orientation only matters when an interaction is not
        # invariant under the orbital transpose. PairHop is excluded (its
        # declaration partner is the HERMITIAN entry, a different pairing);
        # PairLift is excluded because its particle-hole S/C contribution
        # is exactly zero on both the producer and consumer sides -- an
        # asymmetric PairLift cannot change either the stored chi or the
        # vertex, so rejecting on it would be a pure false positive.
        orientation = []
        resym = []
        if str(chi_convention).lower() == "myo":
            orientation = [
                t for t in ("CoulombInter", "Hund", "Ising", "Exchange")
                if interactions.get(t)
                and _onsite_transpose_asymmetry(interactions[t]) > 0.0]
            # PairHop is NOT exempt wholesale: its orientation never
            # changed, but the conjugated-mean reading of its Hermitian-
            # partner declarations arrived with the version stamp, so a
            # non-Hermitian-closed declaration has unverifiable
            # symmetrisation semantics in an unversioned file.
            if (interactions.get("PairHop")
                    and _onsite_hermitian_pair_mismatch(
                        interactions["PairHop"]) > 0.0):
                resym = ["PairHop"]
        # #232: version 3 puts an off-site row (R, a, b) on the bond of a
        # Transfer row (R, a, b); every earlier file was computed with it on
        # the opposite bond. That only matters for a table that is not
        # symmetric under R -> -R.
        # PairLift is exempt for the same reason as in the orientation rule
        # above: its particle-hole content is exactly zero on both sides.
        placement = placement_sensitive_types(
            interactions,
            grid=(_grid if "CellShape" in input_dict.get("mode", {}).get(
                "param", {}) else None))
        if affected or orientation or resym or placement:
            versions = {}
            for data, path in ((data_s, chi_s_path), (data_c, chi_c_path)):
                if "sc_vertex_version" not in data:
                    reasons = []
                    if placement:
                        reasons.append(
                            "{} declare(s) off-site rows that are not "
                            "symmetric under R -> -R, which a file this old "
                            "placed on the bond opposite their Transfer rows "
                            "(#232)".format(", ".join(placement)))
                    if affected:
                        reasons.append(
                            "the interaction set contains {} whose vertex "
                            "content changed in the #113 corrections".format(
                                ", ".join(affected)))
                    if resym:
                        reasons.append(
                            "PairHop declares on-site entries that are "
                            "not Hermitian partners, whose symmetrised "
                            "reading cannot be verified for a file this "
                            "old: the conjugated-mean reading of PairHop "
                            "declarations arrived with the version stamp "
                            "(#101)")
                    if orientation:
                        reasons.append(
                            "{} declare(s) an asymmetric on-site "
                            "inter-orbital coupling, whose orientation and "
                            "symmetrisation semantics cannot be verified "
                            "for a file this old: the interaction "
                            "orientation changed in #99 and the "
                            "declaration symmetrisation arrived with the "
                            "version-2 stamp, so an unversioned file may "
                            "pair either or both differently from the "
                            "current vertex (#101)".format(
                                ", ".join(orientation)))
                    raise ValueError(
                        "FLEX susceptibility file '{}' predates the current "
                        "vertex conventions (no sc_vertex_version field), "
                        "and {}. Pairing the old susceptibilities with the "
                        "current vertices would silently mix two different "
                        "interactions -- regenerate the susceptibilities "
                        "with the current code.".format(
                            path, "; and ".join(reasons)))
                # strict decode: the tag must be a single integral scalar.
                # A plain int() would silently truncate (2.9 -> 2, accepted)
                # and non-finite values would escape as OverflowError.
                shown = "<unreadable>"
                try:
                    raw = data["sc_vertex_version"]
                    shown = repr(raw)
                    arr = np.asarray(raw)
                    if arr.dtype.kind not in "iufc":    # no bool, str or object stamps
                        raise ValueError("not numeric")
                    if arr.size != 1:
                        raise ValueError("not a scalar")
                    val = complex(arr.reshape(())[()])
                    if val.imag != 0.0 or not (
                            np.isfinite(val.real)
                            and float(val.real).is_integer()):
                        raise ValueError("not an integer")
                    versions[path] = int(val.real)
                except (TypeError, ValueError, OverflowError):
                    raise ValueError(
                        "FLEX susceptibility file '{}' carries a malformed "
                        "sc_vertex_version field ({}).".format(path, shown))
            if len(set(versions.values())) != 1:
                raise ValueError(
                    "FLEX chi_s and chi_c carry different sc_vertex_version "
                    "values ({}); they must come from the same run.".format(
                        versions))
            ver = next(iter(versions.values()))
            if ver not in (2, 3):
                raise ValueError(
                    "FLEX susceptibility files carry sc_vertex_version = {} "
                    "but this code supports version 2 or version 3 (the "
                    "#113 exact-diagonalization vertex content; version 3 "
                    "also places off-site rows on the bond of their Transfer "
                    "rows, #232); "
                    "regenerate the susceptibilities with the current "
                    "code.".format(ver))
            if ver == 2 and placement:
                raise ValueError(
                    "FLEX susceptibility files carry sc_vertex_version = 2, "
                    "and {} declare(s) off-site rows that are not symmetric "
                    "under R -> -R. Version-2 runs placed such a row on the "
                    "bond opposite its Transfer row (#232), so pairing those "
                    "susceptibilities with the current vertices would "
                    "silently mix two different interactions -- regenerate "
                    "the susceptibilities with the current code.".format(
                        ", ".join(placement)))
    # Self-describing index-order marker (issue #78 follow-up). Current files
    # always store [a,c,b,d]; the marker exists so a pre-#78 dev output (the
    # general path stored MYO-transposed arrays under the SAME "myo" tag,
    # indistinguishable by tag alone) is rejected instead of silently building
    # a transposed pairing vertex, and so any future layout change fails fast.
    # Marker-less "kuroki"/untagged files stay readable: the reduced-path
    # layout never changed.
    for data, path in ((data_s, chi_s_path), (data_c, chi_c_path)):
        if "chi_orbital_layout" in data:
            layout = str(data["chi_orbital_layout"])
            if layout != "acbd":
                raise ValueError(
                    "FLEX chi file '{}' declares chi_orbital_layout='{}' but "
                    "this reader only supports 'acbd'.".format(path, layout))
        elif chi_convention == "myo":
            raise ValueError(
                "FLEX chi file '{}' is tagged chi_convention='myo' but lacks "
                "the chi_orbital_layout marker: it was produced by a pre-fix "
                "development build that stored the arrays orbital-pair-"
                "TRANSPOSED (issue #78) and would yield a wrong pairing "
                "vertex. Re-run FLEX with the current build to regenerate "
                "it.".format(path))
    if not allow_ir:
        return (chi_s_raw, chi_c_raw, chi_convention,
                (_legacy_up_block_only(data_s),
                 _legacy_up_block_only(data_c)))

    native_s, native_c = is_ir_native(data_s), is_ir_native(data_c)
    if native_s != native_c:
        raise ValueError(
            "all dynamic-Eliashberg inputs must share one encoding; re-run "
            "FLEX so chiq_s/chiq_c/green are all densified or all IR-native "
            "(chiq_s: {}, chiq_c: {}).".format(
                "IR-native" if native_s else "densified",
                "IR-native" if native_c else "densified"))
    if native_s and not (chi_convention_present_s and chi_convention_present_c):
        raise ValueError(
            "FLEX chi_s/chi_c are IR-native but missing 'chi_convention'; "
            "the file appears to be an MYO/IR-native FLEX output and its "
            "orbital layout (spin-orbital 'kuroki' vs orbital-pair 'myo') "
            "cannot be safely defaulted -- re-run FLEX with a build that "
            "tags chi_convention, or re-tag the npz explicitly.")
    ir_meta = ({"chis": ir_native_meta(data_s),
                "chic": ir_native_meta(data_c)} if native_s else None)
    return (chi_s_raw, chi_c_raw, chi_convention, ir_meta,
            (_legacy_up_block_only(data_s),
             _legacy_up_block_only(data_c)))


def _expand_flex_chi(chi_raw, norb, Nx, Ny, Nz, convention):
    """Reshape H-wave chi ``(nfreq, nvol, nd, nd)`` to
    ``(nfreq, Nx, Ny, Nz, nd, nd)`` in the ``nd = norb^2`` Eliashberg space,
    resolving the spin-orbital-vs-orbital-pair layout from ``convention``.

    The two FLEX conventions have DIFFERENT physical layouts that shape alone
    cannot always tell apart:

    - ``"myo"`` (general full-vertex FLEX) is already in orbital-pair space
      ``nd_chi = norb^2``; passed through unchanged.
    - ``"kuroki"`` (reduced FLEX) is in spin-orbital reduced space
      ``nd_chi = norb*ns`` (spin-block ordered ``s*norb + a``), and its matrix
      index is a DENSITY PAIR: ``X[a, b]`` is ``chi_{(a,a),(b,b)}``.  The
      spin-up block ``[:norb, :norb]`` is extracted and embedded at the
      density-pair positions ``[a*norb + a, b*norb + b]`` of the
      ``norb^2 x norb^2`` orbital-pair space, with every off-density component
      left exactly zero (see the body for why any other placement is wrong).

    Whether the spin-orbital block must be extracted is decided by ``nd_chi``:
    ``nd_chi == norb*ns`` means spin-orbital (extract), ``nd_chi == norb^2``
    means orbital-pair (pass through). These two collide ONLY for ``norb == 2``
    (``norb^2 == norb*ns == 4``); there the shape is ambiguous and the layout is
    resolved from ``convention`` (``"kuroki"`` -> spin-orbital extract,
    ``"myo"`` -> orbital-pair). The previous shape-only heuristic silently
    treated a norb=2 kuroki spin-orbital chi as orbital-pair, skipping the
    extraction and building a wrong pairing vertex.

    For ``norb != 2`` the shape is unambiguous, but ``convention`` is still
    REQUIRED to agree with the shape-inferred layout (and to be a recognized
    value): the caller forwards ``convention`` unchanged to
    ``_compute_vertices_flex``, which uses it (not the shape) to pick the
    MYO vs Kuroki S/C matrices, so a shape/tag mismatch would silently build
    the wrong pairing vertex just as in the norb=2 case.

    The mapping is elementwise in frequency, so it may be applied to a single
    static slice or the full axis identically.
    """
    ns = 2
    nd = norb * norb          # orbital-pair dimension
    nd_so = norb * ns         # spin-orbital reduced dimension
    nfreq = chi_raw.shape[0]
    # Derive nd_chi from the trailing AXES, not from sqrt of the flattened
    # element count: the latter accepts any rectangle whose product happens to
    # be a perfect square -- (nfreq, nvol, 2, 8) has 16 trailing elements and
    # was silently reshaped into a 4x4 susceptibility.
    #
    # Two layouts arrive here: the reduced/kuroki one is (nfreq, nvol, nd, nd),
    # the general/myo one is (nfreq, nvol, norb, norb, norb, norb) -- an
    # orbital-pair matrix stored with its pair indices unflattened.
    tail = chi_raw.shape[2:]
    if len(tail) == 2 and tail[0] == tail[1]:
        nd_chi = tail[0]
    elif len(tail) == 4 and len(set(tail)) == 1:
        nd_chi = tail[0] * tail[1]
    else:
        raise ValueError(
            "FLEX chi must be (nfreq, nvol, nd, nd) or "
            "(nfreq, nvol, norb, norb, norb, norb); got shape {}.".format(
                chi_raw.shape))
    chi_full = chi_raw.reshape(nfreq, Nx, Ny, Nz, nd_chi, nd_chi)

    if nd_chi == nd and nd_chi == nd_so:
        # ambiguous (norb == 2): the convention tag is the only disambiguator,
        # and choosing the wrong branch silently corrupts the pairing vertex, so
        # require an explicitly known tag rather than defaulting on mismatch.
        if convention == "kuroki":
            is_spin_orbital = True
        elif convention == "myo":
            is_spin_orbital = False
        else:
            raise ValueError(
                "norb=2 FLEX chi has the shape-ambiguous dimension nd_chi={} "
                "(norb^2 == norb*ns); a known chi_convention ('myo' or "
                "'kuroki') is required to resolve the layout, got '{}'.".format(
                    nd_chi, convention))
    elif nd_chi == nd_so and nd_chi != nd:
        # Unambiguous spin-orbital shape: the convention tag must still agree,
        # otherwise a file shaped spin-orbital but mistagged "myo" (or with an
        # unrecognized tag) would be extracted correctly here yet the wrong
        # tag would later select the MYO S/C matrices downstream, silently
        # building the wrong pairing vertex.
        if convention not in ("myo", "kuroki"):
            raise ValueError(
                "FLEX chi has unrecognized chi_convention='{}' (expected "
                "'myo' or 'kuroki').".format(convention))
        if convention != "kuroki":
            raise ValueError(
                "FLEX chi dimension nd_chi={} (norb={}) is unambiguously "
                "spin-orbital (nd_so=norb*ns={}, nd=norb^2={}) but is tagged "
                "chi_convention='{}'; expected 'kuroki'.".format(
                    nd_chi, norb, nd_so, nd, convention))
        is_spin_orbital = True
    elif nd_chi == nd and nd_chi != nd_so:
        # Unambiguous orbital-pair shape: same agreement check, mirrored.
        if convention not in ("myo", "kuroki"):
            raise ValueError(
                "FLEX chi has unrecognized chi_convention='{}' (expected "
                "'myo' or 'kuroki').".format(convention))
        if convention != "myo":
            raise ValueError(
                "FLEX chi dimension nd_chi={} (norb={}) is unambiguously "
                "orbital-pair (nd=norb^2={}, nd_so=norb*ns={}) but is tagged "
                "chi_convention='{}'; expected 'myo'.".format(
                    nd_chi, norb, nd, nd_so, convention))
        is_spin_orbital = False
    else:
        hint = ""
        if nd_chi == norb and norb % 2 == 0:
            # The spin-orbital case. FLEX writes chi in its reduced spin-orbital
            # space, nd = norb_phys * ns; with geom.dat's norb being the
            # spin-orbital count (= 2 * norb_phys) that lands on nd_chi == norb.
            # Say so -- the bare dimensions look like a malformed file.
            hint = (" nd_chi == norb, which is what a "
                    "spin-orbital FLEX run (mode.enable_spin_orbital = true) "
                    "produces: FLEX writes chi with the PHYSICAL orbital count "
                    "while hwave_sc reads norb from geom.dat, where "
                    "spin-orbital mode stores the spin-orbital count "
                    "(= 2 x physical). The Eliashberg step does not support "
                    "spin-orbital / spin-mixing models -- its pairing vertex is "
                    "paramagnetic -- so this combination is not usable rather "
                    "than merely misconfigured.")
        raise ValueError(
            "FLEX chi dimension nd_chi={} matches neither the orbital-pair "
            "size norb^2={} nor the spin-orbital size norb*ns={} (norb={})."
            "{}".format(nd_chi, nd, nd_so, norb, hint))

    if not is_spin_orbital:
        return chi_full

    # Spin-orbital reduced (spin-block ordered s*norb+a): extract the spin-up
    # block and embed it at the DENSITY-PAIR positions of the orbital-pair
    # space.
    #
    # The reduced scheme keeps only the density-density diagonal of the
    # susceptibility: the stored X[a, b] IS chi_{(a,a),(b,b)} (its companion
    # interaction reduction in FLEX._inflate_chi0q_and_ham is
    # einsum('ksasatbtb->ksatb', ...), the density-density diagonal of the
    # vertex).  With the orbital-pair flat index (l1,l2) -> l1*norb + l2 used by
    # the S/C builders, the one faithful embedding is therefore
    #
    #     out[(a,a), (b,b)] = X[a, b],   every other component zero.
    #
    # The historical placement out[(l1,l2), (l3,l2)] = X[l1, l3] (a delta_{l2,l4}
    # "spectator" scatter, i.e. kron(X, I_norb)) instead read the density-pair
    # index as an ordinary orbital index.  That dropped the genuine
    # inter-orbital density coupling chi_{(0,0),(1,1)} -- which S @ chi @ S does
    # reference -- and scattered X into pair positions the reduced scheme never
    # computed, so chi0q_mode="flex" on a reduced run disagreed with the
    # equivalent "load" run for identical Sigma=0 physics.
    #
    # Off-density rows/columns stay exactly zero: a reduced run carries no
    # information about pair indices (a,b) with a != b, and fabricating them
    # from the density block is what caused the mismatch.  See
    # _compute_vertices_flex, which warns when the interaction actually needs
    # those missing components.
    chi_orb = chi_full[:, :, :, :, :norb, :norb]
    out = np.zeros((nfreq, Nx, Ny, Nz, nd, nd), dtype=complex)
    dens = np.arange(norb) * norb + np.arange(norb)   # flat index of (a,a)
    out[..., dens[:, None], dens[None, :]] = chi_orb
    return out


def _load_green_npz(green_path, norb, Nx, Ny, Nz, label="Green"):
    """Load a Green-function npz written in the H-wave layout, path-explicit.

    Shared by ``_load_flex_green`` (which resolves the path from the FLEX
    output directory) and by the ``[eliashberg] bond_green`` key (an
    externally supplied Green function, spec Goal). The file must hold a
    ``green`` array of shape ``(nblock, nfreq, nvol, norb, norb)`` -- what
    ``flex.save_results`` writes -- with ``nvol = Nx*Ny*Nz`` and the same
    ``norb`` as the model; anything else raises rather than being reshaped
    into silent nonsense.

    Returns ``(green, nfreq)`` with ``green`` in the sc.py layout
    ``(norb, norb, Nx, Ny, Nz, nfreq)``.
    """
    if not os.path.exists(green_path):
        raise FileNotFoundError(
            "{} file not found: {}".format(label, green_path))
    logger.info("Loading {} from: {}".format(label, green_path))
    data_g = np.load(green_path)
    _reject_ir_native(data_g, green_path, _STATIC_IR_HINT)
    if "green" not in data_g.files:
        raise ValueError(
            "{} file {} has no 'green' array (keys: {}); the expected format "
            "is the H-wave green.npz written by the FLEX/UHF solvers."
            .format(label, green_path, list(data_g.files)))
    green_raw = data_g["green"]
    if green_raw.ndim != 5:
        raise ValueError(
            "{} file {}: 'green' must have shape (nblock, nfreq, nvol, norb, "
            "norb), got {}".format(label, green_path, green_raw.shape))
    nblock_g, nfreq, nvol_g, norb1, norb2 = green_raw.shape
    if nblock_g != 1:
        raise ValueError(
            "{} file {}: 'green' carries {} spin blocks, but this loader "
            "unconditionally consumes block 0 -- the bond/FLEX Green "
            "consumers here are single-block. A file with nblock=2 (e.g. a "
            "legitimate spin-diagonal H-wave Green with separate G_up/G_down "
            "blocks) would otherwise be silently halved to its up block "
            "while being recorded as a trusted full Green. Re-run the "
            "producing calculation paramagnetically (nblock=1), or extract "
            "and save the single block you intend to use."
            .format(label, green_path, nblock_g))
    if (nvol_g != Nx * Ny * Nz) or (norb1 != norb) or (norb2 != norb):
        raise ValueError(
            "{} file {} does not match the model: it holds nvol={}, norb={}x{} "
            "while this run has nvol={} ({}x{}x{}) and norb={}."
            .format(label, green_path, nvol_g, norb1, norb2,
                    Nx * Ny * Nz, Nx, Ny, Nz, norb))
    green = green_raw[0].reshape(
        nfreq, Nx, Ny, Nz, norb, norb).transpose(4, 5, 1, 2, 3, 0).copy()
    return green, nfreq


class _UnsupportedNpyHeaderVersion(Exception):
    """Internal to :func:`_peek_green_npz_nfreq`: raised when a ``.npy``
    header version is one neither numpy's public ``read_array_header_X_Y``
    wrappers nor its internal version-generic dispatcher can parse.

    This is deliberately NOT a plain ``ValueError`` catch target inside
    :func:`_peek_green_npz_nfreq`'s own except clause: a genuinely unsupported
    version must surface as a clear error (review fix), not be swallowed into
    a silent ``None`` that reopens the double-preflight path the header peek
    exists to close.
    """


def _read_npy_header_shape(fh):
    """Read only the ``shape`` from a ``.npy`` header, for any format version
    ``np.load`` itself accepts (1.0, 2.0, 3.0, and any future version numpy
    grows a dedicated ``read_array_header_X_Y`` wrapper for).

    Dispatches on the version tuple from ``np.lib.format.read_magic`` to the
    matching public ``read_array_header_X_Y`` function when one exists (1.0,
    2.0). Format 3.0 (and any future version numpy's public API has not grown
    a dedicated wrapper for) has no such function, so this falls back to
    ``numpy.lib.format``'s own internal version-generic dispatcher -- exactly
    what ``np.load`` uses under the hood to parse ANY version it accepts, so
    the fallback keeps every such version working here too instead of
    silently giving up.

    Raises :class:`_UnsupportedNpyHeaderVersion` if the version is genuinely
    unknown to :mod:`hwave.solver.npy_header`.
    """
    try:
        return _npy_header.read_npy_header_shape(fh)
    except _npy_header.UnsupportedNpyHeaderVersion as exc:
        raise _UnsupportedNpyHeaderVersion(str(exc))


def _peek_green_npz_nfreq(green_path, label="Green"):
    """Return the ``nfreq`` axis of a Green npz WITHOUT materialising it.

    The bond path's resource preflight must run with the frequency count the
    run will really use -- and for an external ``bond_green`` that is the
    FILE's, not the configured ``Nmat`` (review fix: the preflight used to run
    with the configured value first, so an over-large ``Nmat`` could refuse a
    run that fits). Reading the npz member's array header alone keeps the
    preflight strictly ahead of the (possibly large) allocation.

    Supports every ``.npy`` header version ``np.load`` accepts (see
    :func:`_read_npy_header_shape`); a genuinely unsupported version raises
    rather than silently returning ``None``, since falling back to the
    configured Nmat for the first preflight is exactly the double-preflight
    behaviour this function exists to eliminate.

    Returns ``None`` when the frequency count cannot be determined from the
    header alone for an ordinary reason (missing file, no ``green`` member,
    unexpected rank); the caller then falls back to :func:`_load_green_npz`,
    which raises the informative error.
    """
    try:
        if not zipfile.is_zipfile(green_path):
            return None
        with zipfile.ZipFile(green_path) as zf:
            member = None
            for name in zf.namelist():
                if name in ("green.npy", "green"):
                    member = name
                    break
            if member is None:
                return None
            with zf.open(member) as fh:
                shape = _read_npy_header_shape(fh)
    except _UnsupportedNpyHeaderVersion as exc:
        raise ValueError(
            "{} file {}: {}".format(label, green_path, exc)) from exc
    except (OSError, ValueError):
        # unreadable/corrupt: let the real loader produce the diagnostic
        return None
    if len(shape) != 5:
        return None
    logger.info("%s %s carries %d Matsubara frequencies (read from the npz "
                "header)", label, green_path, int(shape[1]))
    return int(shape[1])


def _load_flex_green(input_dict, norb, Nx, Ny, Nz, allow_ir=False,
                     interactions=None):
    """Load the FLEX dressed Green's function, or ``None`` if absent.

    Returns shape ``(norb, norb, Nx, Ny, Nz, nfreq)`` (full fermionic axis;
    ``nfreq`` = Nmat for uniform files, the node count for IR-native ones).
    With ``allow_ir=True`` returns ``(green, ir_meta)`` instead (``ir_meta``
    is ``None`` for uniform files) -- arity changes only under that flag.
    """
    _, _, green_path = _resolve_flex_paths(input_dict)
    if not os.path.exists(green_path):
        return (None, None) if allow_ir else None
    logger.info("Loading FLEX dressed Green from: {}".format(green_path))
    data_g = np.load(green_path)
    if interactions is not None:
        _grid_cs = input_dict.get("mode", {}).get("param", {}).get("CellShape")
        check_placement_stamp(
            data_g, green_path, interactions,
            grid=(tuple(int(x) for x in list(_grid_cs) + [1] * (3 - len(_grid_cs)))
                  if _grid_cs is not None else None),
            what="FLEX Green file")
    if not allow_ir:
        _reject_ir_native(data_g, green_path, _STATIC_IR_HINT)
    green_raw = data_g["green"]
    # H-wave format: (nblock, nfreq, nvol, norb, norb)
    if green_raw.ndim != 5:
        raise ValueError(
            "green file '{}': expected a rank-5 array (nblock, nfreq, "
            "nvol, norb, norb), got shape {}.".format(
                green_path, green_raw.shape))
    nblock, nmat_g, nvol, norb1, norb2 = green_raw.shape
    # Validate every dimension BEFORE any use: an element-count-compatible
    # but misshapen file would otherwise be silently reinterpreted by the
    # reshape below (measured: (1, 3, 1, 4, 4) on a 4-site 2-orbital
    # system was accepted and scrambled), an empty frequency axis would
    # contract to an all-zero kernel and return a bogus eigenvalue, and a
    # zero-block file died on an incidental IndexError.
    expected_nvol = Nx * Ny * Nz
    if nblock < 1 or nmat_g < 1:
        raise ValueError(
            "green file '{}': needs at least one spin block and one "
            "frequency, got shape {}.".format(green_path, green_raw.shape))
    if nvol != expected_nvol or norb1 != norb or norb2 != norb:
        raise ValueError(
            "green file '{}': shape {} does not match this run "
            "(nvol = Nx*Ny*Nz = {}, norb = {}).".format(
                green_path, green_raw.shape, expected_nvol, norb))
    # A spin-diag FLEX run writes TWO blocks, G_up and G_down. Taking block 0
    # would discard the down-spin propagator exactly as the reduced chi loader
    # used to discard the down-spin susceptibility -- and the pair bubble built
    # from it feeds the Eliashberg kernel directly. The susceptibility guard
    # normally rejects such a run first, but relying on that is fragile: it is
    # a different file, and a run whose chi happens to be redundant while its
    # Green functions are not would slip straight through. (The dynamic loader
    # even reads the Green function BEFORE the susceptibility check.) So the
    # same policy applies here: discarding the other blocks is lossless
    # exactly when they are redundant; when they are not, the pair bubble
    # cannot represent the stored physics and the loader RAISES unless the
    # user takes responsibility via [eliashberg] accept_up_block_only.
    # Non-finite content is rejected up front (below); the block
    # comparison then runs on finite data with plain component equality.
    # Non-finite content is rejected at the file boundary, as the
    # susceptibility loader already does: a NaN/Inf dressed Green function
    # reaches the pair bubble G2 directly, and the iteration solver can
    # turn the resulting non-finite norm into a SAVED zero eigenvalue with
    # converged=False (round-7 review). Per-frequency scan, bounded
    # temporaries.
    for b in range(nblock):
        for f in range(nmat_g):
            if not np.all(np.isfinite(green_raw[b, f])):
                raise ValueError(
                    "green file '{}' contains non-finite values (block "
                    "{}, frequency index {}); a NaN/Inf dressed Green "
                    "function cannot feed the pair bubble. Re-run the "
                    "producing FLEX calculation.".format(green_path, b, f))

    if nblock > 1:
        def _differs(b):
            # Per-frequency, per-component comparison (finiteness is
            # already established above, so plain equality suffices);
            # frequency-sliced scanning bounds the comparison temporaries
            # to one (nvol, norb, norb) slice. .real/.imag are views.
            for f in range(green_raw.shape[1]):
                gb = green_raw[b, f]
                g0 = green_raw[0, f]
                if not (np.array_equal(gb.real, g0.real)
                        and np.array_equal(gb.imag, g0.imag)):
                    return True
            return False

        for b in range(1, nblock):
            if _differs(b):
                if _accept_up_block_only(input_dict):
                    logger.warning(
                        "green file '%s' carries %d spin blocks with "
                        "DIFFERING content; [eliashberg] "
                        "accept_up_block_only is set, so the up block "
                        "alone feeds the pair bubble and the discarded "
                        "spin content is the user's responsibility.",
                        green_path, nblock)
                    break
                raise ValueError(
                    "green file '{}' carries {} spin blocks whose contents "
                    "differ (block {} != block 0). The Eliashberg pair "
                    "bubble is built from the up block only, so a "
                    "spin-resolved dressed Green function cannot be "
                    "represented -- the discarded block is real physics, "
                    "not redundancy. Re-run FLEX paramagnetically, or set "
                    "[eliashberg] accept_up_block_only = true to take "
                    "responsibility for the approximation.".format(
                        green_path, nblock, b))
    # Grid/temperature provenance (round-3 review). The consumer rebuilds
    # the Matsubara frequencies from ITS OWN beta and assumes the file holds
    # the full centered fermionic grid, so a Green function produced at a
    # different temperature -- or a densified/restricted frequency axis --
    # would be consumed silently with a wrong grid, and the tail correction
    # would then apply a wrong analytic complement. Uniform files written
    # since this change carry beta; older files get a warning, not an error.
    if not is_ir_native(data_g):
        # An odd uniform axis cannot be a centered fermionic grid (it
        # contains wn = 0); reject at the file boundary, before any vertex
        # or pair-bubble construction (round-6 review).
        if nmat_g % 2 != 0:
            raise ValueError(
                "green file '{}': the frequency axis has {} points, which "
                "cannot be a centered fermionic Matsubara grid (even count "
                "required). Re-run the producing FLEX calculation.".format(
                    green_path, nmat_g))
        if "freq_index" in data_g:
            fi = np.asarray(data_g["freq_index"]).ravel()
            if fi.size != nmat_g or not np.array_equal(
                    fi, np.arange(nmat_g)):
                raise ValueError(
                    "green file '{}': freq_index does not describe the "
                    "full centered frequency grid of the stored axis "
                    "({} entries for {} frequencies); the pair bubble "
                    "needs the full fermionic grid. Re-run the producing "
                    "FLEX calculation with the full grid.".format(
                        green_path, fi.size, nmat_g))
        run_T = input_dict.get("mode", {}).get("param", {}).get("T")
        run_beta = None if run_T is None else _coerce_run_beta(run_T)
        if "beta" in data_g:
            # One REAL-NUMERIC finite positive scalar, validated on the
            # NARROWED float64 (see _finite_positive_float64): NaN/Inf make
            # every tolerance comparison False (silently "matching"), an
            # empty array would die on an incidental IndexError, a vector
            # would silently use its first element, float() would silently
            # discard a complex value's imaginary part, booleans would read
            # as 0/1, a string would raise an incidental TypeError, and a
            # wider-than-binary64 longdouble could pass elementwise checks
            # yet narrow to inf.
            file_beta = _finite_positive_float64(data_g["beta"])
            if file_beta is None:
                raise ValueError(
                    "green file '{}': beta metadata must be a single "
                    "finite positive real scalar, got {!r}. Regenerate "
                    "the file.".format(green_path,
                                       np.asarray(data_g["beta"])))
            # Explicitly symmetric relative tolerance (max of both
            # magnitudes), no absolute floor: an absolute term would wave
            # through large relative mismatches at small beta (high
            # temperature), where the grid differs the most.
            if run_beta is not None and abs(file_beta - run_beta) \
                    > 1.0e-8 * max(abs(file_beta), abs(run_beta)):
                raise ValueError(
                    "green file '{}' was produced at beta = {} but this "
                    "run uses beta = {}; the Matsubara grid (and the tail "
                    "correction built on it) would be inconsistent. Use "
                    "the producing run's temperature or regenerate the "
                    "file.".format(green_path, file_beta, run_beta))
        elif run_beta is not None:
            logger.warning(
                "green file '%s' carries no beta metadata (written before "
                "this field existed); the run's beta = %g is assumed to "
                "match the producing FLEX run. Verify the temperatures "
                "agree -- a mismatch silently corrupts the pair bubble.",
                green_path, run_beta)
    from hwave.solver.rpa import validate_momentum_convention
    # green layout (nblock, nfreq_or_nodes, nvol, norb, norb): q axis 2,
    # for uniform AND IR-native files alike (round-4 review)
    validate_momentum_convention(data_g, green_path, green_raw, 2,
                                 (Nx, Ny, Nz))
    # Convert to sc.py format: (norb, norb, Nx, Ny, Nz, nfreq)
    green = green_raw[0].reshape(
        nmat_g, Nx, Ny, Nz, norb, norb
    ).transpose(4, 5, 1, 2, 3, 0).copy()
    if not allow_ir:
        return green
    meta = ir_native_meta(data_g) if is_ir_native(data_g) else None
    return green, meta


def _load_flex_susceptibilities(input_dict, norb, Nx, Ny, Nz,
                                interactions=None):
    """Load FLEX-computed susceptibilities at the static (zero bosonic
    frequency) limit from NPZ files.

    Parameters
    ----------
    input_dict : dict
        Parsed TOML configuration.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        Grid dimensions.

    Returns
    -------
    chis : ndarray
        Spin susceptibility at static limit, shape (Nx, Ny, Nz, nd, nd).
    chic : ndarray
        Charge susceptibility at static limit, shape (Nx, Ny, Nz, nd, nd).
    green_dressed : ndarray or None
        Dressed Green's function if available, shape (norb, norb, Nx, Ny, Nz, nmat).
    chi_convention : str
        Orbital convention of chis/chic: "myo" (general full-vertex FLEX) or
        "kuroki" (reduced FLEX / legacy files). Pass this to
        _compute_vertices_flex so the matching S/C matrices are used.
    """
    # Read the raw chi (H-wave layout, frequency = axis 0) WITHOUT expanding,
    # so the static slice is taken before the spin-orbital expansion -- the full
    # loader would otherwise allocate the whole Nmat-long expanded array only to
    # keep one frequency (a memory regression proportional to Nmat).
    (chi_s_raw, chi_c_raw, chi_convention,
     legacy_tags) = _read_flex_chi_raw(
        input_dict, interactions=interactions)

    # The zero bosonic frequency is located via the freq_index/nmat metadata
    # (RPA chiq files can carry a restricted matsubara_frequency axis whose
    # center is NOT the static limit); FLEX files always hold the full grid.
    # Re-derive the same paths the raw reader used so the metadata read here
    # is guaranteed to describe the same files.
    chi_s_path, chi_c_path, _ = _resolve_flex_paths(input_dict)
    config_nmat = input_dict.get("mode", {}).get("param", {}).get(
        "Nmat", _DEFAULT_NMAT)

    def _static_center(path, nfreq):
        data = np.load(path)
        fi, fn = _read_freq_meta(data)
        pos = _static_freq_position(fi, nfreq, config_nmat, path,
                                    file_nmat=fn)
        # no usable metadata: center of the actual data axis (axis 0 in the
        # H-wave chiq layout)
        return nfreq // 2 if pos is None else pos

    center_s = _static_center(chi_s_path, chi_s_raw.shape[0])
    center_c = _static_center(chi_c_path, chi_c_raw.shape[0])

    # Slice the static frequency FIRST, then expand only that single slice.
    # Check the WHOLE stored axis, not just the slice consumed here.
    #
    # Both choices have a failure mode and neither case is common: validating
    # only the static slice lets through a run that is redundant at omega=0 and
    # polarized elsewhere, which then returns a paramagnetic eigenvalue in
    # silence; validating everything rejects a file whose unused frequencies are
    # malformed. The first failure is silent and the second is loud and
    # diagnosable, so the axis is checked in full. Paramagnetism is a property
    # of the producing run rather than of one frequency, which is the same
    # reason.
    legacy_s, legacy_c = legacy_tags
    if _accept_up_block_only(input_dict):
        legacy_s = legacy_c = "config_override"
    _check_spin_block_discarded(chi_s_raw, norb, chi_convention, "chi_s",
                                legacy_s)
    _check_spin_block_discarded(chi_c_raw, norb, chi_convention, "chi_c",
                                legacy_c)
    chis = _expand_flex_chi(chi_s_raw[center_s:center_s + 1],
                            norb, Nx, Ny, Nz, chi_convention)[0]
    chic = _expand_flex_chi(chi_c_raw[center_c:center_c + 1],
                            norb, Nx, Ny, Nz, chi_convention)[0]

    green_dressed = _load_flex_green(input_dict, norb, Nx, Ny, Nz,
                                     interactions=interactions)
    return chis, chic, green_dressed, chi_convention


# ---------------------------------------------------------------------------
# G2 and Eliashberg kernel
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Gap function initialization
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# Solvers
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Rayleigh validation of a SHIFTED power iteration
#
# ``||(K + sigma I) v|| - sigma`` is a signed eigenvalue of K only once the
# iteration has converged to a mode whose SHIFTED eigenvalue is positive real.
# Neither an explicit sigma (which the user may set too small) nor
# ``"auto"`` (a finite-step power estimate of the spectral radius, which
# approaches it from below) guarantees that, and the NON-converged exit
# returns the same quantity. Example: ``K = diag(-3, -8, -5)`` with
# ``sigma = 1`` drives the iterate to the -8 mode, whose shifted eigenvalue is
# -7, so the raw quantity is ``|-7| - 1 = +6`` -- not an eigenvalue of K at
# all. The signed Rayleigh quotient below (one extra matvec on the UNSHIFTED
# operator) plus its residual decide whether the run may call its number an
# eigenvalue; nothing here touches the unshifted (default) path.
# ---------------------------------------------------------------------------


def _shifted_eigenvalue_note(solver_mode, spectral_shift, converged, n_iter,
                             iter_info, kernel_label="kernel"):
    """The ``eigenvalue.dat`` comment describing a power iteration's number.

    Split out of ``calc_eliashberg`` so the static and dynamic entry points --
    and the tests -- share one wording, and so the "this is NOT an eigenvalue"
    case can never be silently dropped from the file. ``converged`` is the
    ``eigen_converged`` verdict (``_solve_leading``'s ``info["converged"]``):
    when it is False the reported value is NaN and the note says why.
    ``spectral_shift`` None means the unshifted iteration
    (``spectral_shift = "none"``).
    """
    info = iter_info or {}
    resid = info.get("rayleigh_residual")
    resid_txt = ("{:.3e}".format(resid) if isinstance(resid, float)
                 and np.isfinite(resid) else str(resid))
    shift_txt = ("spectral_shift={!r}".format(spectral_shift)
                 if spectral_shift is not None else "no spectral_shift")
    if not converged:
        if info.get("two_cycle"):
            why = ("the iterate settled into a sign-alternating two-cycle "
                   "after {} steps -- the dominant (largest-magnitude) "
                   "eigenvalue of the iterated operator is NEGATIVE"
                   .format(n_iter))
            cyc = info.get("cycle_eigenvalue")
            if cyc is not None:
                why += " ({:.8e} on the {})".format(cyc, kernel_label)
        elif info.get("diff_converged"):
            why = ("the iterate difference converged after {} steps but the "
                   "Rayleigh check left a residual ||K v - lambda v||/||v|| = "
                   "{} (tolerance {})".format(
                       n_iter, resid_txt, info.get("rayleigh_residual_tol")))
        else:
            why = "the iteration did NOT converge after {} steps".format(n_iter)
        return (
            "solver_mode='{}' with {}: {}. eigen_converged = false: NO "
            "eigenvalue is reported (the value below is NaN, NOT an "
            "eigenvalue) and the gap written with it is the last iterate, NOT "
            "an eigenvector. Use solver_mode=\"eigenvalue\", or "
            "spectral_shift=\"auto\" (the default for the iteration), or "
            "raise max_iter.".format(solver_mode, shift_txt, why))
    if spectral_shift is None:
        return (
            "solver_mode='{}' with no spectral_shift: the value below is the "
            "SIGNED Rayleigh eigenvalue <v|K|v>/<v|v> of the converged "
            "iterate (residual {}); the iteration converged after {} steps. "
            "The unshifted iteration follows the largest-MAGNITUDE "
            "eigenvalue.".format(solver_mode, resid_txt, n_iter))
    note = (
        "solver_mode='{}' with spectral_shift={!r}: the power iteration "
        "ran on the shifted {} K + sigma*I and the value below is the "
        "SIGNED eigenvalue of K (the shift has been subtracted back), "
        "VALIDATED by the Rayleigh quotient <v|K|v>/<v|v> with residual "
        "||K v - lambda v||/||v|| = {}; the iteration converged after {} "
        "steps.".format(solver_mode, spectral_shift, kernel_label, resid_txt,
                        n_iter))
    if info.get("shift_sufficient") is False:
        note += (" WARNING: spectral_shift appears INSUFFICIENT -- the "
                 "mode found has lambda + sigma < 0, so it is the "
                 "largest-MAGNITUDE eigenvalue of K + sigma*I, not the "
                 "algebraically largest (physical) one; re-run with a "
                 "larger spectral_shift or \"auto\".")
    return note


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _save_results(output_dir, sigma, eigenvalue, eigenvalues_eig, kx_array, ky_array, kz_array,
                  gap_file="gap.dat", eigenvalue_file="eigenvalue.dat",
                  eigenvalue_match=None, provenance=None,
                  eigenvalue_note=None, eigenvalue_weights=None,
                  gap_comments=None, eigen_converged=None,
                  iteration_converged=None, analysis_note=None,
                  leading_eigenvalue=None):
    """Save gap function and eigenvalue results to files.

    Parameters
    ----------
    output_dir : str
        Output directory path.
    sigma : ndarray or None
        Gap function from iteration, shape (norb, norb, Nx, Ny, Nz).
    eigenvalue : float or None
        Leading eigenvalue from iteration.
    eigenvalues_eig : ndarray or None
        Eigenvalues from eigenvalue analysis.
    kx_array, ky_array, kz_array : ndarray
        k-point arrays.
    gap_file : str
        Output filename for gap function.
    eigenvalue_file : str
        Output filename for eigenvalues.
    provenance : dict, optional
        Extra approximation-level metadata (currently the bond-channel record,
        spec S5). Written as ``# key = value`` COMMENT lines at the top of the
        eigenvalue file, so this is a strictly ADDITIVE schema change: existing
        readers (``np.loadtxt``, ``tsweep.parse_leading_eig``) skip ``#`` lines
        and every existing key/column is untouched. When None (the default) the
        file is byte-for-byte what it always was.
    eigenvalue_note : str, optional
        Free text describing WHAT the single ``eigenvalue`` number is, written
        as ``#`` comment line(s) between the ``# Iteration eigenvalue`` header
        and the value itself. Used to state that an iteration-mode value is an
        unsigned iterate norm (and whether it converged), so it can never be
        silently confused with the signed ``lambda_rayleigh`` sitting in the
        provenance block above it. Additive in the same sense as
        ``provenance``: comment lines only, no numeric row changes.
    eigenvalue_weights : ndarray, optional
        Spin-orbital mode only: ``(n_eig, 4)`` spin weights
        ``(psi, dx, dy, dz)`` of each eigenvector
        (:func:`hwave.solver.so_pairing.so_spin_decomposition`), appended as
        four columns after the ``match`` column. None (the default) leaves the
        eigenvalue rows unchanged.
    gap_comments : list of str, optional
        Spin-orbital mode only: ``#`` comment lines written before the column
        header of the gap file. None (the default) leaves the gap file
        unchanged.
    eigen_converged : bool, optional
        Written as ``# eigen_converged = <true|false>`` -- right after
        ``# Iteration eigenvalue`` when the iteration ran, else at the top.
        ``false``: the power iteration did not converge, the iteration value
        is NaN and the gap file holds the last iterate, not an eigenvector.
        With ``iteration_converged`` given (``solver_mode = "both"``) it
        refers to the eigenvalue analysis and is written right before
        ``# Eigenvalue analysis``.
    iteration_converged : bool, optional
        ``solver_mode = "both"`` only: the power iteration's own verdict,
        written as ``# iteration_converged = <true|false>`` after
        ``# Iteration eigenvalue``.
    analysis_note : str, optional
        Where the eigenvalue analysis' reported leader comes from (the
        channel's parity-sector solve), written as ``#`` lines right after
        the ``# eigen_converged`` line.
    leading_eigenvalue : complex, optional
        The eigenvalue analysis' reported leader (NaN when none), written as
        the machine-readable line ``# leading_eigenvalue = <Re> <Im>`` right
        before ``# Eigenvalue analysis``; it equals row 0 of the table.
    """
    os.makedirs(output_dir, exist_ok=True)

    if sigma is not None:
        norb = sigma.shape[0]
        Nx, Ny, Nz = len(kx_array), len(ky_array), len(kz_array)
        filepath = os.path.join(output_dir, gap_file)
        logger.info("Saving gap function to {}".format(filepath))

        with open(filepath, "w") as fw:
            for line in (gap_comments or ()):
                fw.write("# {}\n".format(line))
            # Header
            header_parts = ["# kx", "ky", "kz"]
            for i in range(norb):
                for j in range(norb):
                    header_parts.append("Re(sigma_{}{})".format(i, j))
                    header_parts.append("Im(sigma_{}{})".format(i, j))
            fw.write(" ".join(header_parts) + "\n")

            for ix in range(Nx):
                kx = kx_array[ix]
                if kx > np.pi:
                    kx -= 2.0 * np.pi
                for iy in range(Ny):
                    ky = ky_array[iy]
                    if ky > np.pi:
                        ky -= 2.0 * np.pi
                    for iz in range(Nz):
                        kz = kz_array[iz]
                        if kz > np.pi:
                            kz -= 2.0 * np.pi
                        parts = ["{:.8f}".format(kx),
                                 "{:.8f}".format(ky),
                                 "{:.8f}".format(kz)]
                        for i in range(norb):
                            for j in range(norb):
                                val = sigma[i, j, ix, iy, iz]
                                parts.append("{:.8e}".format(val.real))
                                parts.append("{:.8e}".format(val.imag))
                        fw.write(" ".join(parts) + "\n")

    if eigenvalue is not None or eigenvalues_eig is not None:
        filepath = os.path.join(output_dir, eigenvalue_file)
        logger.info("Saving eigenvalues to {}".format(filepath))
        with open(filepath, "w") as fw:
            if provenance:
                for key in sorted(provenance):
                    fw.write("# {} = {}\n".format(key, provenance[key]))
            flag_line = (None if eigen_converged is None else
                         "# eigen_converged = {}\n".format(
                             str(bool(eigen_converged)).lower()))
            if eigenvalue is not None:
                fw.write("# Iteration eigenvalue\n")
                if iteration_converged is not None:
                    fw.write("# iteration_converged = {}\n".format(
                        str(bool(iteration_converged)).lower()))
                elif flag_line:
                    fw.write(flag_line)
                    if analysis_note:
                        fw.write("# {}\n".format(analysis_note))
                if eigenvalue_note:
                    for line in str(eigenvalue_note).splitlines():
                        fw.write("# {}\n".format(line))
                fw.write("{:.8e}\n".format(eigenvalue))
            if eigenvalues_eig is not None:
                if flag_line and (eigenvalue is None
                                  or iteration_converged is not None):
                    fw.write(flag_line)
                    if analysis_note:
                        fw.write("# {}\n".format(analysis_note))
                if leading_eigenvalue is not None:
                    lv = complex(leading_eigenvalue)
                    fw.write("# leading_eigenvalue = {:.8e} {:.8e}\n".format(
                        lv.real, lv.imag))
                fw.write("# Eigenvalue analysis\n")
                # The optional trailing `match` column is 1 when the eigenvector
                # lies (dominantly) in the channel's parity sector (even for
                # singlet, odd for triplet) and 0 when it lies in the opposite
                # sector. The leading index/Re/Im/|ev| columns are unchanged, so
                # position-based parsers (e.g. reading column 1 as Re) are
                # unaffected; only the column count grows from 4 to 5.
                if eigenvalue_weights is not None:
                    fw.write("# index  Re(eigenvalue)  Im(eigenvalue)  "
                             "|eigenvalue|  match(1=fermion-antisymmetric)  "
                             "w_psi  w_dx  w_dy  w_dz\n")
                elif eigenvalue_match is not None:
                    fw.write("# index  Re(eigenvalue)  Im(eigenvalue)  "
                             "|eigenvalue|  match(1=channel-parity)\n")
                else:
                    fw.write("# index  Re(eigenvalue)  Im(eigenvalue)  "
                             "|eigenvalue|\n")
                for i, ev in enumerate(eigenvalues_eig):
                    row = "{:4d} {:15.8e} {:15.8e} {:15.8e}".format(
                        i, ev.real, ev.imag, abs(ev))
                    if eigenvalue_match is not None:
                        row += " {:d}".format(int(bool(eigenvalue_match[i])))
                    if eigenvalue_weights is not None:
                        row += "".join(" {:15.8e}".format(w)
                                       for w in eigenvalue_weights[i])
                    fw.write(row + "\n")


# ---------------------------------------------------------------------------
# chi0q format conversion
# ---------------------------------------------------------------------------

def _convert_chi0q_to_ref_format(chi0q, norb, Nx, Ny, Nz):
    """Convert chi0q from H-wave format to reference code format.

    Supports both 2-index (reduced) and 4-index (general) chi0q:

    2-index:
        H-wave: (nmat, nvol, norb, norb) -> Ref: (norb, norb, Nx, Ny, Nz, nmat)
    4-index:
        H-wave: (nmat, nvol, norb, norb, norb, norb)
        -> Ref: (norb, norb, norb, norb, Nx, Ny, Nz, nmat)

    Parameters
    ----------
    chi0q : ndarray
        chi0q in H-wave format.
    norb, Nx, Ny, Nz : int
        System parameters.

    Returns
    -------
    chi0q_ref : ndarray
        chi0q in reference format:
        - 2-index: (norb, norb, Nx, Ny, Nz, nmat)
        - 4-index: (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
    """
    if chi0q.ndim == 4:
        # 2-index: (nmat, nvol, norb, norb) -> (norb, norb, Nx, Ny, Nz, nmat)
        nf = chi0q.shape[0]
        chi0q_3d = chi0q.reshape(nf, Nx, Ny, Nz, norb, norb)
        chi0q_ref = chi0q_3d.transpose(4, 5, 1, 2, 3, 0)
    elif chi0q.ndim == 6:
        if chi0q.shape[0] == norb and chi0q.shape[1] == norb:
            # Already in ref format (norb, norb, Nx, Ny, Nz, nmat)
            chi0q_ref = chi0q
        else:
            # 4-index H-wave: (nmat, nvol, norb, norb, norb, norb)
            # -> (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
            nf = chi0q.shape[0]
            chi0q_3d = chi0q.reshape(nf, Nx, Ny, Nz, norb, norb, norb, norb)
            chi0q_ref = chi0q_3d.transpose(4, 5, 6, 7, 1, 2, 3, 0)
    elif chi0q.ndim == 8:
        # Already in ref format (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
        chi0q_ref = chi0q
    else:
        raise ValueError(
            "Unexpected chi0q shape: {}. Expected 4D, 6D, or 8D.".format(chi0q.shape)
        )
    return chi0q_ref


# ---------------------------------------------------------------------------
# Main calculation
# ---------------------------------------------------------------------------

def reject_spin_orbital_mode(input_dict):
    """Raise when [mode] enable_spin_orbital is set, for the Eliashberg
    entries that do not support it: the dynamic bond path
    (``eliashberg_bond_io.solve_dynamic_bond``) and the temperature sweep.
    The static entry (:func:`calc_eliashberg`, issue #83) and the on-site
    dynamic solver (``frequency = "dynamic"``, issue #284) support
    spin-orbital models and refuse only their unsupported combinations
    (:func:`hwave.solver.so_pairing.check_so_supported`).

    Until issue #83 a spin-orbital config RAN through these entries and
    printed eigenvalues built on four internal inconsistencies (norb never
    halved to the physical count; the internally computed chi0q in
    interleaved order against spin-block consumers; interaction files read
    at the spin-orbital dimension; paramagnetic orbital-space S/C rules
    applied to spin-orbital indices). A silent wrong result, not an
    approximation.

    Resolution is strict (shared _resolve_spin_orbital_flag): recognized
    false forms proceed as non-SO, unrecognized values raise -- a user who
    misspells "true" must not silently get a non-SO run that then
    misreads spin-orbital-formatted inputs.
    """
    if _resolve_spin_orbital_flag(input_dict):
        raise ValueError(
            "[mode] enable_spin_orbital = true is not supported by this "
            "Eliashberg entry (the dynamic bond path, bond_channels = true, "
            "and the temperature sweep): it would produce silently wrong "
            "eigenvalues (issue #83). hwave_sc supports spin-orbital models "
            "with frequency = \"static\" (chi0q_mode = \"calc\", \"load\" "
            "or \"flex\") and with frequency = \"dynamic\" (chi0q_mode = "
            "\"calc\" or \"flex\", on-site interactions).")


def _so_rpa_vertex(input_dict, so_context=None):
    """The antisymmetrized spinful vertex ``Gamma`` ``(nvol, nd^2, nd^2)``
    (:func:`hwave.solver.rpa.spinful_bubble_vertex`) of a spin-orbital
    config, built from the interaction objects of the spin-orbital RPA
    (``ham_inter_q`` and ``spinful_exchange()``, the tensors the RPA and FLEX
    SO resum; interaction files in physical orbital indices). With
    ``so_context`` (the entry context of the hwave_sc run) the RPA reuses it:
    no second compile of the interaction, no second ``orbital_site``
    WARNING."""
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.rpa as sol_rpa

    info_inputfile = dict(input_dict.get("file", {}).get("input", {}))
    info_inputfile["path_to_input"] = info_inputfile.get("path_to_input", "")
    info_mode_rpa = dict(input_dict.get("mode", {}))
    info_mode_rpa["calc_scheme"] = "general"
    info_mode_rpa["_orbital_site_reported"] = True
    info_mode_rpa["enable_spin_orbital"] = True
    if so_context is not None:
        from hwave.solver import spin_orbital
        spin_orbital.context_for_nested(info_mode_rpa, so_context)
    read_io = read_input_k.QLMSkInput(info_inputfile)
    solver = sol_rpa.RPA(read_io.get_param("ham"),
                         {"print_level": 1, "print_step": 1}, info_mode_rpa)
    nvol, nd = solver.lattice.nvol, 2 * solver.norb
    ham = np.asarray(backend.to_host(solver.ham_info.ham_inter_q)).reshape(
        nvol, nd, nd, nd, nd)
    return sol_rpa.spinful_bubble_vertex(ham, solver.ham_info.spinful_exchange())


def _so_gap_npz_name(gap_file):
    """The spin-orbital gap archive: the gap file name with its extension
    replaced by ``.npz`` (``gap.dat`` -> ``gap.npz``); a gap file that is
    itself named ``*.npz`` gets ``.npz`` appended instead, so the two never
    collide."""
    name = os.path.splitext(gap_file)[0] + ".npz"
    return gap_file + ".npz" if name == gap_file else name


#: producer labels of the hwave_sc spin-orbital gap archives (#88)
_SO_STATIC_LABEL = "hwave_sc (spin-orbital, static)"
_SO_DYNAMIC_LABEL = "hwave_sc (spin-orbital, dynamic)"


def _so_flex_writer_members(stage, label, so_run, borrowed, green_sc, epsilon_k,
                            beta):
    """The writer members of a run whose kernel used a borrowed spin-orbital
    FLEX Green function: its digest as consumed (the digest of the archive
    check), the mu the kernel applied (the validated archive value) and the
    closure density of that G (:func:`_so_flex_consumed_state`)."""
    state = _so_flex_consumed_state(
        green_sc, epsilon_k, borrowed.sigma_static, borrowed.info.producer_state.mu,
        beta, so_run.state.target)
    return _so_writer_members(stage, label, so_run, borrowed.fp_green, state)


def _so_prov_line(v):
    """The ``# key = value`` text of a scalar provenance member."""
    a = np.asarray(v)
    if a.dtype.kind == "b":
        return "true" if bool(a) else "false"
    if a.dtype.kind == "f":
        return repr(float(a))
    if a.dtype.kind in "iu":
        return str(int(a))
    return str(a[()])


def _so_provenance_lines(members):
    """The scalar members of ``members`` as ``{key: text}`` (the arrays stay
    in the npz only)."""
    return {k: _so_prov_line(v) for k, v in (members or {}).items()
            if np.ndim(v) == 0}


def _so_output_record(output_dir, gap_file, norb_phys, gap, eigenvalues,
                      eigenvectors, match, kx_array, ky_array, kz_array,
                      sz_residual=None, provenance_members=None,
                      eigen_converged=None):
    """The spin-orbital output schema (issue #83).

    ``eigen_converged`` (bool, optional) is stored as the archive key of
    that name: False when the gap is the last iterate of a power iteration
    that did not converge (not an eigenvector).

    Writes the gap archive (:func:`_so_gap_npz_name`, default ``gap.npz``)
    and returns ``(provenance, weights, gap_comments)`` for
    :func:`_save_results`:

    * ``provenance`` -- the ``# key = value`` lines of ``eigenvalue.dat``
      (``spin_orbital = true``, ``index_convention = spin_block``, the gap
      layout, the spin decomposition and the weights of the written gap);
    * ``weights`` -- ``(n_eig, 4)`` spin weights ``(psi, dx, dy, dz)`` of the
      eigenvectors (:func:`hwave.solver.so_pairing.so_spin_decomposition`),
      or None without an eigenvalue analysis;
    * ``gap_comments`` -- the comment lines of ``gap.dat``.

    The archive holds ``gap`` ``(nd, nd, Nx, Ny, Nz)`` (the gap of
    ``gap.dat``), ``index_convention``, ``spin_orbital``, ``norb_phys``,
    ``kx``/``ky``/``kz``, ``spin_components`` and ``gap_spin_weights``; with
    an eigenvalue analysis also ``eigenvalues``, ``eigenvectors``
    ``(n_eig, nd, nd, Nx, Ny, Nz)``, ``eigenvector_spin_weights`` and
    ``eigenvalue_match``.

    With ``provenance_members`` (#88: the fingerprint members of
    :func:`hwave.solver.fingerprint.so_provenance_members` for this run's band
    and the Green function and state its kernel used) the archive gets every
    member and ``provenance`` gets the scalar ones.

    With ``sz_residual`` (the S_z-conservation residual of the run,
    :func:`hwave.solver.so_pairing.sz_conservation_residual`) the S_z report
    of :func:`hwave.solver.so_pairing.sz_report` is logged and its members
    (``sz_weight_upup``, ``sz_weight_dndn``, ``sz_weight_singlet0``,
    ``sz_weight_triplet0`` ``(n_ev,)``, ``sz_conservation_residual``,
    ``sz_conserved``) are added to the archive; ``n_ev`` counts the
    eigenvectors, or is 1 (the written gap) without an eigenvalue analysis."""
    from hwave.solver import so_pairing, spin_orbital
    comps = so_pairing.SPIN_COMPONENTS
    nd = 2 * int(norb_phys)
    layout = ("Delta_AB(k), A = spin * {} + orb (spin 0 = up, spin 1 = down); "
              "fermion-antisymmetric, Delta_AB(k) = -Delta_BA(-k)"
              .format(norb_phys))
    npz_name = _so_gap_npz_name(gap_file)
    prov = {
        "spin_orbital": "true",
        "index_convention": spin_orbital.INDEX_CONVENTION,
        "gap_layout": layout,
        "gap_npz": npz_name,
        "match_column": "1 = fermion-antisymmetric eigenvector",
        "spin_decomposition": (
            "Delta_ab(k) = (psi_ab(k) + d_ab(k) . sigma) i sigma_y per orbital "
            "pair (a, b); w_X = sum_{k,a,b} |X_ab(k)|^2 / sum_{k,a,b} "
            "(|psi_ab(k)|^2 + |d_ab(k)|^2), X = psi, dx, dy, dz"),
        "weight_columns": "w_psi w_dx w_dy w_dz (after the match column)",
    }
    record = {"index_convention": np.array(spin_orbital.INDEX_CONVENTION),
              "spin_orbital": np.array(True),
              "norb_phys": np.array(int(norb_phys)),
              "kx": np.asarray(kx_array), "ky": np.asarray(ky_array),
              "kz": np.asarray(kz_array),
              "spin_components": np.array(comps)}
    if gap is not None:
        gap = np.asarray(backend.to_host(gap))
        w = so_pairing.so_spin_decomposition(gap, norb_phys)
        prov["gap_spin_weights"] = " ".join(
            "{}={:.8e}".format(c, w[c]) for c in comps)
        record["gap"] = gap.reshape((nd, nd) + gap.shape[2:])
        record["gap_spin_weights"] = np.array([w[c] for c in comps])
    if eigen_converged is not None:
        record["eigen_converged"] = np.array(bool(eigen_converged))
    weights = None
    if eigenvalues is not None:
        vecs = np.stack([np.asarray(backend.to_host(v)) for v in eigenvectors])
        weights = np.array([[so_pairing.so_spin_decomposition(v, norb_phys)[c]
                             for c in comps] for v in vecs])
        record["eigenvalues"] = np.asarray(eigenvalues)
        record["eigenvectors"] = vecs
        record["eigenvector_spin_weights"] = weights
        if match is not None:
            record["eigenvalue_match"] = np.asarray(match, dtype=bool)
    if sz_residual is not None:
        if eigenvalues is not None:
            sz_gaps = vecs
        elif gap is not None:
            sz_gaps = [gap]
        else:
            sz_gaps = []
        if len(sz_gaps):
            record.update(so_pairing.sz_report(
                [so_pairing.so_sz_decomposition(v, norb_phys) for v in sz_gaps],
                sz_residual, logger))
    if provenance_members:
        # #88: the fingerprint members of this run (its band, the G and state its
        # kernel used): the npz holds all of them, eigenvalue.dat the scalars
        record.update(provenance_members)
        prov.update(_so_provenance_lines(provenance_members))
    gap_comments = [
        "spin_orbital = true",
        "index_convention = {}".format(spin_orbital.INDEX_CONVENTION),
        "gap_layout = {}; sigma_ij below is Delta_ij".format(layout),
    ]
    if gap is not None:
        os.makedirs(output_dir, exist_ok=True)
        path = os.path.join(output_dir, npz_name)
        logger.info("Saving the spin-orbital gap to {}".format(path))
        np.savez(path, **record)
    return prov, weights, gap_comments


def _so_static_vertex_from_chi(gam, chi_static, nd, spatial_shape):
    """The static spin-orbital pairing vertex ``Vs_q``
    ``(nd, nd, nd, nd, Nx, Ny, Nz)`` from the spinful vertex ``gam`` and a
    static spinful susceptibility ``chi_static`` (both ``(nvol, nd^2,
    nd^2)``, bubble ``[a, c, b, d]`` layout, spin-block order):
    :func:`hwave.solver.so_pairing.so_pairing_vertex`. Pure; the susceptibility
    may come from the RPA (:func:`_so_static_vertex`) or from FLEX."""
    from hwave.solver.so_pairing import so_pairing_vertex
    return so_pairing_vertex(gam, chi_static, nd, spatial_shape)


def _so_static_chi0(chi0q_raw, nd, Nx, Ny, Nz, static_index):
    """The static spinful bubble ``(nvol, nd^2, nd^2)`` (bubble ``[a, c, b,
    d]`` layout) of a general spinful chi0q -- raw ``(nmat, nvol, nd, nd, nd,
    nd)`` or reference ``(nd, nd, nd, nd, Nx, Ny, Nz, nmat)`` -- at the zero
    bosonic frequency ``static_index`` (``nmat // 2`` when ``None``). Any other
    layout (reduced 2-index, spin-free / spin-diagonal) is refused."""
    nvol = Nx * Ny * Nz
    chi0q_raw = np.asarray(chi0q_raw)
    if chi0q_raw.ndim == 6 and chi0q_raw.shape[1:] == (nvol,) + (nd,) * 4:
        nmat = chi0q_raw.shape[0]
        si = nmat // 2 if static_index is None else static_index
        chi0 = chi0q_raw[si].reshape(nvol, nd * nd, nd * nd)
    elif chi0q_raw.ndim == 8 and chi0q_raw.shape[:7] == (nd,) * 4 + (Nx, Ny, Nz):
        nmat = chi0q_raw.shape[-1]
        si = nmat // 2 if static_index is None else static_index
        chi0 = chi0q_raw[..., si].reshape(nd * nd, nd * nd, nvol).transpose(2, 0, 1)
    else:
        raise ValueError(
            "[mode] enable_spin_orbital = true: the chi0q has shape {}, not the "
            "general (4-index) spinful layout over nd = {} spin-orbitals "
            "(nmat, {}, {nd}, {nd}, {nd}, {nd}). A reduced (2-index) or a "
            "spin-free / spin-diagonal chi0q cannot carry the spin-flip "
            "particle-hole pairs; use chi0q_mode = \"calc\", or load the "
            "chi0q of a spin-orbital RPA run with calc_scheme = \"general\" and a "
            "spin-mixing transfer.".format(chi0q_raw.shape, nd, nvol, nd=nd))
    return chi0


def _so_static_vertex(input_dict, chi0q_raw, norb_phys, Nx, Ny, Nz,
                      static_index, so_context=None):
    """The static spin-orbital pairing vertex of the RPA paths
    (``chi0q_mode = "calc"`` / ``"load"``): ``Gamma`` from the spin-orbital RPA
    objects (:func:`_so_rpa_vertex`), ``chi = [1 + chi0 Gamma]^-1 chi0`` at the
    zero bosonic frequency ``static_index`` (``nmat // 2`` when ``None``) of
    the general spinful ``chi0q_raw`` -- raw ``(nmat, nvol, nd, nd, nd, nd)``
    or reference ``(nd, nd, nd, nd, Nx, Ny, Nz, nmat)`` -- then
    :func:`_so_static_vertex_from_chi`. ``so_context`` is passed to
    :func:`_so_rpa_vertex`."""
    nd = 2 * int(norb_phys)
    chi0 = _so_static_chi0(chi0q_raw, nd, Nx, Ny, Nz, static_index)
    chi0 = np.ascontiguousarray(chi0, dtype=complex)
    gam = _so_rpa_vertex(input_dict, so_context=so_context)
    eye = np.eye(nd * nd, dtype=complex)
    chi = np.linalg.solve(eye[None] + chi0 @ gam, chi0)
    return _so_static_vertex_from_chi(gam, chi, nd, (Nx, Ny, Nz))


#: ``chi_convention`` of the single antisymmetrized spinful channel a
#: spin-orbital FLEX run stores as ``chiq`` (flex.py save_results).
_SO_FLEX_CHI_CONVENTION = "spin_orbital"

#: Backward-error bound of the Gamma consistency check of
#: :func:`_check_so_flex_vertex` (round-off of the FLEX solve is ~1e-15; an
#: interaction that differs from the FLEX run's shows up at O(dU/U)).
_SO_FLEX_GAMMA_TOL = 1.0e-8


def _so_flex_static_slice(data, arr, path, nd, nvol, config_nmat, what):
    """The zero-bosonic-frequency slice ``(nvol, nd^2, nd^2)`` of a raw
    spinful FLEX array ``arr`` ``(nmat, nvol, nd, nd, nd, nd)`` of the open
    npz ``data`` read from ``path``, located through the file's
    ``freq_index`` / ``nmat`` metadata."""
    if arr.ndim != 6 or arr.shape[1:] != (nvol,) + (nd,) * 4:
        raise ValueError(
            "FLEX {} in '{}' has shape {}, not the general spinful layout "
            "(nmat, {}, {nd}, {nd}, {nd}, {nd}) over nd = {nd} spin-orbitals of "
            "this run. Re-run FLEX on the same model and grid.".format(
                what, path, arr.shape, nvol, nd=nd))
    fi, fn = _read_freq_meta(data)
    pos = _static_freq_position(fi, arr.shape[0], config_nmat, path,
                                file_nmat=fn)
    si = arr.shape[0] // 2 if pos is None else pos
    out = np.ascontiguousarray(arr[si].reshape(nvol, nd * nd, nd * nd),
                               dtype=complex)
    if not np.all(np.isfinite(out)):
        raise ValueError(
            "FLEX {} in '{}' contains non-finite values at the zero bosonic "
            "frequency; re-run the producing FLEX calculation.".format(what, path))
    return out


#: IR-native refusal hint of the spin-orbital dynamic reader (issue #284)
_SO_DYNAMIC_IR_HINT = (
    "the spin-orbital dynamic Eliashberg solver reads uniform-grid files only "
    "(the IR basis is not supported in spin-orbital mode yet). Re-run FLEX "
    "with [mode.param] write_densified = true.")


def _so_flex_full_axis(data, arr, path, nd, nvol, config_nmat, what):
    """The whole bosonic axis ``(nmat, nvol, nd^2, nd^2)`` of a raw spinful
    FLEX array ``arr`` ``(nmat, nvol, nd, nd, nd, nd)`` for the dynamic
    solver (issue #284): the stored axis must be the full uniform grid of the
    configured ``Nmat`` (``freq_index``, when recorded, ``0 .. Nmat - 1``;
    ``nmat``, when recorded, ``Nmat``), with bosonic index ``m`` at
    ``nu_m = (2 m - Nmat) pi T``."""
    if arr.ndim != 6 or arr.shape[1:] != (nvol,) + (nd,) * 4:
        raise ValueError(
            "FLEX {} in '{}' has shape {}, not the general spinful layout "
            "(nmat, {}, {nd}, {nd}, {nd}, {nd}) over nd = {nd} spin-orbitals of "
            "this run. Re-run FLEX on the same model and grid.".format(
                what, path, arr.shape, nvol, nd=nd))
    fi, fn = _read_freq_meta(data)
    nmat = int(config_nmat)
    full = (arr.shape[0] == nmat and (fn is None or int(fn) == nmat)
            and (fi is None or np.array_equal(np.asarray(fi).ravel(),
                                              np.arange(nmat))))
    if not full:
        raise ValueError(
            "FLEX {} in '{}' does not hold the full bosonic grid of [mode.param] "
            "Nmat = {} (stored axis {}, recorded nmat {}): frequency = "
            "\"dynamic\" needs the susceptibility at every bosonic frequency of "
            "the run's grid. Re-run FLEX with the same Nmat.".format(
                what, path, nmat, arr.shape[0], fn))
    out = np.ascontiguousarray(arr.reshape(nmat, nvol, nd * nd, nd * nd),
                               dtype=complex)
    for m in range(nmat):
        if not np.all(np.isfinite(out[m])):
            raise ValueError(
                "FLEX {} in '{}' contains non-finite values at bosonic index {}; "
                "re-run the producing FLEX calculation.".format(what, path, m))
    return out


class _SOFlexBorrowed(NamedTuple):
    """What a spin-orbital FLEX consumer borrowed (#88): ``info``, the
    :class:`hwave.solver.fingerprint.SOArchiveInfo` of the chiq archive;
    ``fp_green``, the digest of the Green function as consumed (computed once,
    in the verification against the archive); ``sigma_static``, the FLEX
    run's static self-energy ``(nvol, nd, nd)`` (for the closure density of
    that Green function)."""
    info: object
    fp_green: object
    sigma_static: object


def _so_flex_sigma_path(eli, flex_dir):
    """The sigma archive of a spin-orbital FLEX run: ``[eliashberg]
    flex_sigma`` (default ``"sigma.npz"``) in ``flex_dir``."""
    return os.path.join(flex_dir, eli.get("flex_sigma", "sigma.npz"))


def _load_flex_so_sigma(eli, flex_dir, nd, nvol, so_run, fp_green, green_path,
                        hint):
    """The static self-energy ``(nvol, nd, nd)`` of the spin-orbital FLEX run
    whose Green function digest is ``fp_green``, from ``[eliashberg]
    flex_sigma`` (default ``"sigma.npz"``): the archive is validated
    (:func:`_so_check_archive`) and must record the same ``fp_green``."""
    from hwave.solver import spin_orbital
    sigma_path = _so_flex_sigma_path(eli, flex_dir)
    if not os.path.exists(sigma_path):
        raise ValueError(
            "FLEX self-energy archive '{}' not found. {}, and also its sigma: "
            "hwave_sc records the closure density of the borrowed Green "
            "function, which needs the static self-energy. Set [file.output] "
            "sigma in the FLEX run and [eliashberg] flex_sigma to that "
            "file.".format(sigma_path, hint))
    with np.load(sigma_path) as ds:
        sinfo = _so_check_archive(ds, sigma_path, so_run, "exact", "sigma",
                                  _SO_FLEX_PRODUCER)
        if sinfo.producer_fp_green != fp_green:
            raise ValueError(
                "{} and {} come from different runs: '{}' records fp_green {} "
                "but the Green function '{}' has digest {}. Use the sigma and "
                "green archives of one spin-orbital FLEX run.".format(
                    os.path.basename(sigma_path), os.path.basename(green_path),
                    sigma_path, (sinfo.producer_fp_green or "none")[:12],
                    green_path, fp_green[:12]))
        conv = str(ds["sigma_convention"]) if "sigma_convention" in ds else None
        tag = str(ds["index_convention"]) if "index_convention" in ds else None
        if conv != "split" or tag != spin_orbital.INDEX_CONVENTION:
            raise ValueError(
                "FLEX self-energy archive '{}' records sigma_convention = {!r} "
                "and index_convention = {!r}; expected 'split' and {!r} (a "
                "spin-orbital Hartree-Fock FLEX run writes them).".format(
                    sigma_path, conv, tag, spin_orbital.INDEX_CONVENTION))
        st = np.asarray(ds["sigma_static"]) if "sigma_static" in ds else None
        want = (1, 1, nvol, nd, nd)
        # complex128 in either byte order (the archive check accepts both);
        # returned in native order
        if (st is None or st.dtype.kind != "c" or st.dtype.itemsize != 16
                or st.shape != want):
            raise ValueError(
                "FLEX self-energy archive '{}': sigma_static must be complex128 "
                "of shape {}, got {}. Regenerate it with a spin-orbital "
                "Hartree-Fock FLEX run.".format(
                    sigma_path, want,
                    None if st is None else "{} {}".format(st.dtype, st.shape)))
    return np.ascontiguousarray(st[0, 0], dtype=np.complex128)


def _load_flex_so(input_dict, nd, Nx, Ny, Nz, interactions=None, full=False,
                  so_run=None):
    """The outputs of a spin-orbital FLEX run (``calc_scheme = "general"``,
    ``[mode] enable_spin_orbital = true``) for ``chi0q_mode = "flex"``.

    ``full = True`` (the dynamic solver, issue #284) returns ``chi`` and
    ``chi0`` on the whole bosonic axis, ``(nmat, nvol, nd^2, nd^2)``
    (:func:`_so_flex_full_axis`), and requires the Green function on the
    configured Nmat; the tag checks are the same.

    Returns ``(chi, green_kw, chi0, chi_path, borrowed)``:

    * ``chi`` ``(nvol, nd^2, nd^2)``: the stored ``chiq`` -- the single
      antisymmetrized spinful channel ``[1 + chi0 Gamma]^-1 chi0`` in the
      bubble ``[a, c, b, d]`` pair layout over spin-block spin-orbitals
      (``Flex._flex_compute_veff_spinful``) -- at the zero bosonic frequency;
    * ``green_kw`` ``(nd, nd, Nx, Ny, Nz, nmat)``: the stored dressed Green
      function (spin-block order, as stored; no reorder);
    * ``chi0`` ``(nvol, nd^2, nd^2)``: the stored static bubble of the same
      map (for :func:`_check_so_flex_vertex`). The chi0q archive is REQUIRED:
      ``chiq`` carries no fingerprint of the vertex it was resummed with, so
      without it interaction files that differ from the FLEX run's would give
      silently wrong eigenvalues;
    * ``chi_path``: the chiq archive read;
    * ``borrowed`` (:class:`_SOFlexBorrowed`), or None without ``so_run``.

    With ``so_run`` (:class:`_SORun`, #88) the provenance of the chiq and
    chi0q archives is checked against this run (:func:`_so_check_archive`,
    exact nmat) after their tag checks, and so are the green and sigma
    archives of the same run, all before any payload is read; the chiq,
    chi0q and green archives are checked again with the digest of the
    loaded Green function, and the sigma archive is then read
    (:func:`_load_flex_so_sigma`, which binds it to the same digest).

    The archive is identified by its TAGS, never by its shape (a one-orbital
    spin-orbital file has nd = 2, the shape of a two-orbital non-spin-orbital
    one): ``chi_convention = "spin_orbital"`` and ``index_convention`` equal
    to :data:`hwave.solver.spin_orbital.INDEX_CONVENTION`. A non-spin-orbital
    ("myo"), reduced ("kuroki") or untagged archive is refused before any
    heavy work, and so is a Green function without the spin-block tag and a
    chi0q without the ``spin_orbital = True`` marker (written only by a
    spin-orbital run)."""
    from requests.structures import CaseInsensitiveDict
    from hwave.solver import spin_orbital
    from hwave.solver.rpa import (validate_momentum_convention,
                                  check_momentum_marker)
    eli = CaseInsensitiveDict(input_dict.get("eliashberg", {}) or {})
    flex_dir = _resolve_flex_dir(input_dict)
    chi_path = os.path.join(flex_dir, eli.get("flex_chiq", "chiq.npz"))
    chi0_path = os.path.join(flex_dir, eli.get("flex_chi0q", "chi0q.npz"))
    # the Green file name is read case-insensitively like the keys above; the
    # shared loader below is handed a view of the config with that name
    green_name = eli.get("flex_green", "green.npz")
    green_path = os.path.join(flex_dir, green_name)
    green_input = dict(input_dict)
    green_input["eliashberg"] = dict(
        {k: v for k, v in (input_dict.get("eliashberg", {}) or {}).items()
         if str(k).lower() != "flex_green"}, flex_green=green_name)
    nvol = Nx * Ny * Nz
    config_nmat = input_dict.get("mode", {}).get("param", {}).get(
        "Nmat", _DEFAULT_NMAT)
    hint = ("[mode] enable_spin_orbital = true with chi0q_mode = \"flex\" reads "
            "the outputs of a spin-orbital FLEX run (calc_scheme = \"general\", "
            "[mode] enable_spin_orbital = true, [file.output] chiq, chi0q, green "
            "and sigma)")

    if not os.path.exists(chi_path):
        raise ValueError(
            "FLEX susceptibility archive '{}' not found. {}; set [eliashberg] "
            "flex_chiq to the file the FLEX run wrote as [file.output] "
            "chiq.".format(chi_path, hint))
    if not os.path.exists(chi0_path):
        # the stored chiq carries no fingerprint of the vertex it was resummed
        # with; the run's chi0q is what pins Gamma (_check_so_flex_vertex)
        raise ValueError(
            "FLEX chi0q archive '{}' not found. {}, and also its chi0q: the "
            "check that this run's interaction files reproduce the FLEX run's "
            "vertex needs it. Set [file.output] chi0q in the FLEX run and "
            "[eliashberg] flex_chi0q to that file.".format(chi0_path, hint))
    ir_hint = _SO_DYNAMIC_IR_HINT if full else _STATIC_IR_HINT
    logger.info("Loading spin-orbital FLEX chiq from: {}".format(chi_path))
    with np.load(chi_path) as data:
        _reject_ir_native(data, chi_path, ir_hint)
        if "chi_convention" not in data:
            raise ValueError(
                "FLEX archive '{}' has no chi_convention tag: an untagged file "
                "cannot be told from a non-spin-orbital one by its shape (one "
                "spin-orbital orbital has the shape of two orbitals). {}.".format(
                    chi_path, hint))
        conv = str(data["chi_convention"])
        if conv != _SO_FLEX_CHI_CONVENTION:
            what = {"myo": "a non-spin-orbital general FLEX output (separate "
                           "spin and charge channels)",
                    "kuroki": "a reduced FLEX output (separate spin and "
                              "charge channels)"}.get(
                        conv.lower(), "not a spin-orbital FLEX output")
            raise ValueError(
                "FLEX archive '{}' carries chi_convention = {!r}: {}, which "
                "cannot carry the spin-flip pairs of a spin-orbital model. "
                "{}.".format(chi_path, conv, what, hint))
        if "index_convention" not in data:
            raise ValueError(
                "FLEX archive '{}' has no index_convention tag; the spin-orbital "
                "order of its axes cannot be verified. {}.".format(chi_path, hint))
        if str(data["index_convention"]) != spin_orbital.INDEX_CONVENTION:
            raise ValueError(
                "FLEX archive '{}' carries index_convention = {!r}, expected {!r}. "
                "{}.".format(chi_path, str(data["index_convention"]),
                             spin_orbital.INDEX_CONVENTION, hint))
        if "chiq" not in data:
            raise ValueError(
                "FLEX archive '{}' has no 'chiq' array (the spinful channel). "
                "{}.".format(chi_path, hint))
        if ("chi_orbital_layout" in data
                and str(data["chi_orbital_layout"]) != "acbd"):
            raise ValueError(
                "FLEX archive '{}' declares chi_orbital_layout = {!r} but this "
                "reader only supports 'acbd'.".format(
                    chi_path, str(data["chi_orbital_layout"])))
        if "scf_converged" in data and not bool(np.asarray(data["scf_converged"])):
            logger.warning(
                "FLEX archive '%s' records scf_converged = false: the FLEX run "
                "did not converge, and its chiq and green may not belong to one "
                "self-consistent state.", chi_path)
        check_momentum_marker(data, chi_path)
        # the chi0q companion must be a spin-orbital FLEX output: its marker
        # is checked after the lightweight chiq metadata above, but before
        # any array (chiq or chi0q) is materialized and before the Green
        # function or any vertex is built (_check_so_flex_vertex still
        # follows: the marker does not establish equal interactions)
        with np.load(chi0_path) as d0:
            _require_so_chi0q_marker(
                d0, chi0_path, producer="FLEX",
                hint="{}; and its chi0q ([file.output] chi0q).".format(hint))
        if so_run is not None:
            _so_check_archive(data, chi_path, so_run, "exact", "chiq",
                              _SO_FLEX_PRODUCER)
            # the chi0q companion carries the same members (FLEX writes one
            # build into every archive); checked before any array is read
            with np.load(chi0_path) as d0:
                _so_check_archive(d0, chi0_path, so_run, "exact", "chi0q",
                                  _SO_FLEX_PRODUCER)
            # the green and sigma archives of the same run against this run
            # (structure, band, interaction, state), also before any payload
            # is read; their Green digest is compared once the Green function
            # is loaded (a missing file gets its own message below)
            sigma_path = _so_flex_sigma_path(eli, flex_dir)
            for path, kind in ((green_path, "green"), (sigma_path, "sigma")):
                if os.path.exists(path):
                    with np.load(path) as other:
                        _so_check_archive(other, path, so_run, "exact", kind,
                                          _SO_FLEX_PRODUCER)
        chiq = data["chiq"]
        if chiq.ndim >= 2 and chiq.shape[1] == nvol:
            validate_momentum_convention(data, chi_path, chiq, 1, (Nx, Ny, Nz))
        if full:
            chi = _so_flex_full_axis(data, chiq, chi_path, nd, nvol, config_nmat,
                                     "chiq")
        else:
            chi = _so_flex_static_slice(data, chiq, chi_path, nd, nvol,
                                        config_nmat, "chiq")
        del chiq

    # the Green function: required (the stored chi belongs to the dressed G)
    if not os.path.exists(green_path):
        raise ValueError(
            "FLEX Green archive '{}' not found. {}; set [eliashberg] flex_green "
            "to the file the FLEX run wrote as [file.output] green.".format(
                green_path, hint))
    with np.load(green_path) as gdata:
        gtag = (str(gdata["index_convention"]) if "index_convention" in gdata
                else None)
    if gtag != spin_orbital.INDEX_CONVENTION:
        raise ValueError(
            "FLEX Green archive '{}' carries index_convention = {!r}, expected "
            "{!r} (a spin-orbital FLEX run tags it). {}.".format(
                green_path, gtag, spin_orbital.INDEX_CONVENTION, hint))
    green_kw = _load_flex_green(green_input, nd, Nx, Ny, Nz,
                                interactions=interactions)
    if green_kw is None or _resolve_flex_paths(green_input)[2] != green_path:
        raise ValueError(
            "internal error: the FLEX Green archive '{}' resolved to a different "
            "path in the shared loader.".format(green_path))
    if full and green_kw.shape[-1] != int(config_nmat):
        raise ValueError(
            "FLEX Green archive '{}' holds {} fermionic frequencies, but "
            "[mode.param] Nmat = {}: frequency = \"dynamic\" needs the Green "
            "function on the run's grid. Re-run FLEX with the same Nmat.".format(
                green_path, green_kw.shape[-1], int(config_nmat)))

    with np.load(chi0_path) as d0:
        _reject_ir_native(d0, chi0_path, ir_hint)
        validate_chi0q_index_convention(d0, True, chi0_path)
        check_momentum_marker(d0, chi0_path)
        if full:
            chi0 = _so_flex_full_axis(d0, d0["chi0q"], chi0_path, nd, nvol,
                                      config_nmat, "chi0q")
        else:
            chi0 = _so_flex_static_slice(d0, d0["chi0q"], chi0_path, nd, nvol,
                                         config_nmat, "chi0q")
    borrowed = None
    if so_run is not None:
        consumed_fp = _so_green_fp(green_kw)
        with np.load(chi_path) as data:
            so_info = _so_check_archive(data, chi_path, so_run, "exact", "chiq",
                                        _SO_FLEX_PRODUCER, green_fp=consumed_fp)
        with np.load(chi0_path) as d0:
            _so_check_archive(d0, chi0_path, so_run, "exact", "chi0q",
                              _SO_FLEX_PRODUCER, green_fp=consumed_fp)
        with np.load(green_path) as gd:
            _so_check_archive(gd, green_path, so_run, "exact", "green",
                              _SO_FLEX_PRODUCER, green_fp=consumed_fp)
        sigma_static = _load_flex_so_sigma(eli, flex_dir, nd, nvol, so_run,
                                           consumed_fp, green_path, hint)
        borrowed = _SOFlexBorrowed(so_info, consumed_fp, sigma_static)
    return chi, green_kw, chi0, chi_path, borrowed


def _so_flex_backward_error(gam, chi, chi0):
    """The backward error of :func:`_check_so_flex_vertex`: the maximum over
    the q blocks of ``|(1 + chi0 gam) chi - chi0| / (|chi0| + |chi0 gam|
    |chi|)`` (Frobenius per block; a zero-scale block counts 0 if its residual
    vanishes and infinity otherwise; NaN propagates). ``chi`` and ``chi0`` are
    ``(..., nd^2, nd^2)`` with any leading axes."""
    c0g = chi0 @ gam
    res = np.linalg.norm(chi + c0g @ chi - chi0, axis=(-2, -1))
    scale = (np.linalg.norm(chi0, axis=(-2, -1))
             + np.linalg.norm(c0g, axis=(-2, -1))
             * np.linalg.norm(chi, axis=(-2, -1)))
    pos = scale > 0
    per_q = np.where(res == 0, 0.0, np.inf)
    per_q = np.where(pos, res / np.where(pos, scale, 1.0), per_q)
    # NaN in any block propagates (and is refused by the callers)
    per_q = np.where(np.isnan(res) | np.isnan(scale), np.nan, per_q)
    return float(np.max(per_q)) if per_q.size else 0.0


def _check_so_flex_vertex(gam, chi, chi0, chi_path):
    """Refuse when ``gam`` (built here from the interaction files) is not the
    vertex the FLEX run resummed: the stored ``chi`` must solve
    ``(1 + chi0 gam) chi = chi0`` with the stored ``chi0`` of the same map.
    Measured as a backward error, ``|(1 + chi0 gam) chi - chi0| /
    (|chi0| + |chi0 gam| |chi|)`` (Frobenius) PER q block, and the maximum
    over q is compared with the tolerance: the bound does not depend on how
    close the run is to an instability, nor on the number of q points (a
    global norm would dilute one bad block by ``1/sqrt(nvol)``). A block of
    zero scale (``chi0_q = 0`` and ``chi0_q gam chi_q = 0``) has error 0 if
    its residual vanishes and infinity otherwise."""
    err = _so_flex_backward_error(gam, chi, chi0)
    logger.info("spin-orbital FLEX: vertex consistency backward error %.3e", err)
    if not err <= _SO_FLEX_GAMMA_TOL:
        raise ValueError(
            "the interaction files of this run do not reproduce the vertex of "
            "the FLEX run that wrote '{}' (backward error {:.3e} > {:.1e} of "
            "(1 + chi0 Gamma) chi = chi0 at the zero bosonic frequency). Use the "
            "same interaction files and [mode.param] as the FLEX run.".format(
                chi_path, err, _SO_FLEX_GAMMA_TOL))
    return err


#: Bound on the relative antisymmetry residual ``||v - P v|| / ||v||`` of a
#: spin-orbital gap (:func:`_so_antisymmetry_residual`). The solvers return
#: gaps that are antisymmetric by construction (residual exactly 0, see
#: :class:`_SOAntisymmetricBasis`); anything above round-off is refused.
_SO_ANTISYM_TOL = 1.0e-10


class SectorOp(NamedTuple):
    """One pair-S_z sector of the spin-orbital pairing kernel.

    ``operator`` acts on the sector vector of ``basis`` (a
    :meth:`_SOAntisymmetricBasis.from_subset` basis); ``has_adjoint = True``
    promises a working ``operator.rmatvec``, and a solver then reduces it
    with ``basis.reduce(operator, adjoint=True)``."""
    pp_sector: int
    operator: object
    basis: object
    has_adjoint: bool


def _so_antisymmetry_residual(gap):
    """``||v - P v|| / ||v||`` of a spin-orbital gap ``(nd, nd, Nx, Ny, Nz)``
    with ``P`` the projection onto fermion-antisymmetric gaps; ``inf`` for a
    zero or non-finite gap. Unlike :func:`_is_gap_parity` (a dominant-sector
    test used by the non-spin-orbital channels) this is a round-off-scale
    criterion: compare with :data:`_SO_ANTISYM_TOL`."""
    gap = np.asarray(gap)
    n = float(np.linalg.norm(gap))
    if n == 0.0 or not np.isfinite(n):
        return float("inf")
    sym = 0.5 * (gap + _reverse_k_and_orbital(gap))
    return float(np.linalg.norm(sym)) / n


def _so_is_antisymmetric(gap):
    return _so_antisymmetry_residual(gap) <= _SO_ANTISYM_TOL


def _so_solve_iteration(A, basis, sigma_init, max_iter, alpha, tol,
                        spectral_shift=None, info_out=None):
    """The power iteration of the spin-orbital kernel on the antisymmetric
    coordinates of ``basis`` (the shared ``_solve_leading`` driver on
    ``E^dagger K E``); returns ``(gap, eigenvalue, converged, n_iter)`` with
    the gap mapped back by ``E``."""
    red = basis.reduce(A)
    c0 = basis.to_reduced(np.asarray(sigma_init).ravel())
    eigenvalue, c, info = _solve_leading(
        lambda: (red, basis.dim), basis.dim, "iteration",
        max_iter=max_iter, convergence_tol=tol, init_vec=c0, alpha=alpha,
        project_fn=None, spectral_shift=spectral_shift)
    if info_out is not None:
        info_out.update(info)
    gap = basis.to_full(c).reshape(basis.shape)
    return gap, eigenvalue, info["converged"], info["n_iter"]


def _so_dense_reduced_eigen(red, dim, num_eigenvalues, seed_vec=None):
    """The dense diagonalization of the reduced operator ``red`` (``dim x
    dim``): ``min(max(1, num_eigenvalues), dim)`` eigenpairs selected by the
    largest real part, ordered by descending real part (by the overlap with
    ``seed_vec`` when given). Returns ``(vals, vecs, M)`` with the
    eigenvectors as columns and ``M`` the dense reduced matrix. Logs
    nothing (the dense branch of :func:`_so_solve_eigenvalue`)."""
    dense = np.asarray(red.matmat(np.eye(dim, dtype=complex)))
    vals, vecs = np.linalg.eig(dense)
    n_keep = min(max(1, int(num_eigenvalues)), dim)
    # the dense spectrum is complete, so select by the largest REAL part
    # (the physical leading eigenvalue) whether or not a spectral_shift
    # was given; a largest-magnitude cut could drop it behind repulsive
    # modes
    sel = np.argsort(-vals.real, kind="stable")[:n_keep]
    vals, vecs = _order_eigenpairs(vals[sel], vecs[:, sel])
    if seed_vec is not None:
        vals, vecs = _order_by_seed_overlap(vals, vecs, seed_vec)
    return vals, vecs, dense


def _so_solve_eigenvalue(A, basis, num_eigenvalues, method, sigma_shift=None,
                         spectral_shift=None, seed_vec=None, map_back=True,
                         info_out=None):
    """The eigenvalue analysis of the spin-orbital kernel on the
    antisymmetric coordinates of ``basis``; returns ``(vals, gaps)`` with
    the eigenvectors mapped back by ``E``, shape ``(n_ev,) + shape``,
    ordered by descending real part.

    ``seed_vec`` (reduced coordinates, optional): an eigenvector-continuation
    seed, used as the ARPACK start vector and to order the eigenpairs by
    their overlap with it (``_solve_leading``; on the dense branch the
    selected eigenpairs are ordered the same way). ``map_back = False``
    returns the reduced eigenvectors as columns ``(dim, n_ev)`` instead of
    the gaps (the dynamic solver maps them back one at a time: each gap is
    ``nd^2 Nk nmat`` long). The defaults are the static solver's call.

    ARPACK needs ``k < dim - 1``; when ``num_eigenvalues >= dim - 1`` the
    reduced operator (then small) is diagonalized densely instead and
    ``min(num_eigenvalues, dim)`` eigenpairs are returned, selected by the
    largest real part (the dense spectrum is complete), so no physical
    eigenpair is dropped at the full-spectrum boundary. Only ``method = "arnoldi"`` is supported.

    ``info_out`` (dict, optional) receives ``selection``: the criterion that
    produced the eigenpairs (``"dense-LR"`` on the dense branch, else
    ``_solve_leading``'s ``"LM"`` / ``"LR"`` / ``"LR_retry"``)."""
    if method != "arnoldi":
        raise ValueError(
            "[mode] enable_spin_orbital = true supports only "
            "eigenvalue_method = \"arnoldi\", got {!r}.".format(method))
    spectral_shift = _normalize_spectral_shift(spectral_shift)
    _validate_spectral_shift(spectral_shift, method)
    red = basis.reduce(A)
    dim = basis.dim
    if num_eigenvalues >= dim - 1:
        if sigma_shift is not None:
            # the same advisory the ARPACK path (_solve_leading) gives
            logger.warning(
                "[eliashberg] sigma_shift is ignored for eigenvalue_method="
                "'arnoldi' (it targets the shift-invert methods); use "
                "spectral_shift to bias the arnoldi selection.")
        vals, vecs, _ = _so_dense_reduced_eigen(red, dim, num_eigenvalues,
                                                seed_vec)
        n_keep = len(vals)
        logger.info("spin-orbital: num_eigenvalues = %d >= dim - 1 = %d of the "
                    "antisymmetric space; dense diagonalization, %d eigenpairs.",
                    num_eigenvalues, dim - 1, n_keep)
        selection = "dense-LR"
    else:
        _, _, eig = _solve_leading(
            lambda: (red, dim), dim, method,
            num_eigenvalues=num_eigenvalues, sigma_shift=sigma_shift,
            spectral_shift=spectral_shift, seed_vec=seed_vec)
        vals, vecs = eig["eigenvalues"], eig["eigenvectors"]
        selection = eig.get("selection")
    if info_out is not None:
        info_out["selection"] = selection
    if not map_back:
        return vals, vecs
    gaps = np.array([basis.to_full(vecs[:, i]).reshape(basis.shape)
                     for i in range(len(vals))])
    return vals, gaps


# --- the per-sector solvers of the S_z-packed spin-orbital kernel -------------------
#
# With S_z conservation the pairing kernel is block-diagonal in the pair-S_z
# sector p = S_z(B) + S_z(C) of a gap component Delta_BC (p in +1, 0, -1). Each
# sector is a SectorOp: an operator on the sector vector and the subset basis
# of its antisymmetric coordinates. The solvers below work sector by sector and
# merge the results; the dense layout keeps _so_solve_eigenvalue /
# _so_solve_iteration, and both layouts meet in the one report
# (_so_bond_report).

#: pair-S_z sectors in their tie-break order
_SO_PP_ORDER = (1, 0, -1)
#: per-pair residual bound ||K v - lambda v|| <= tol * n_K * ||v||
_SO_SECTOR_RESIDUAL_TOL = 1.0e-10
#: eigenvalues closer than this (relative, floor 1) form one tie cluster
_SO_SECTOR_TIE_TOL = 1.0e-12
#: a sector part of the seed at most this fraction of the seed (2-norm) is
#: round-off: the sector is solved unseeded
_SO_SEED_PART_RTOL = 1.0e-8
#: the seeded union warns when no selected pair overlaps the seed this much
_SO_SEED_OVERLAP_WARN = 0.5
#: steps and relative-change stop of the operator-norm power estimate
_SO_NORM_STEPS = 30
_SO_NORM_RTOL = 1.0e-3
#: pp_sector of a result that is not sector-resolved (dense layout, or none)
_SO_PP_DENSE = 2


def _pp_label(p):
    return "+1" if p == 1 else str(int(p))


class SectorLeader(NamedTuple):
    """The eigen solve of one pair-S_z sector.

    ``value``: the largest real part among the returned eigenvalues (``0``
    for a null sector, NaN when the sector is empty or did not converge).
    ``selection``: ``"LM"``, ``"LR_retry"``, ``"LR"`` (a user
    ``spectral_shift``), ``"dense-LR"``, ``"null-operator"``,
    ``"not-converged"``, ``"residual-failed"`` or ``"empty"``.
    ``residual``: the largest ``||K v - lambda v|| / (n_K ||v||)`` over the
    returned pairs (certified when ``<= 1e-10``). ``norm_estimate`` /
    ``norm_kind``: the operator-norm scale of ``n_K`` -- ``"exact"`` (the
    2-norm of the dense matrix), ``"adjoint"`` (power iteration on ``K^dagger
    K``) or ``"forward"`` (power iteration on ``K``, a lower bound).
    ``n_returned``: the number of pairs the sector solve returned."""
    pp_sector: int
    dim: int
    value: complex
    selection: str
    certified: bool
    residual: float
    norm_estimate: float
    norm_kind: str
    n_returned: int


@dataclass(frozen=True)
class SectorSolve:
    """The union of the sector eigen solves.

    ``vals (n,)``, ``gaps (n,) + gap_shape`` and ``pp_sector (n,) int8``:
    the certified pairs of the sectors that did not fail, at most ``k`` of
    them, ordered by descending real part (ties: see
    :func:`_so_solve_sectors`; by the seed overlap in a seeded solve).
    ``leaders``: one :class:`SectorLeader` per given sector, in the order
    ``(+1, 0, -1)``. ``failed_sectors``: the
    sectors whose solve did not converge or returned a pair that failed the
    residual check. ``converged``: no failed sector and at least one pair;
    ``leading`` is then ``vals[0]``, else NaN. ``selection``: the per-sector
    selections, e.g. ``"+1:dense-LR;0:LR;-1:null-operator"``."""
    vals: np.ndarray
    gaps: np.ndarray
    pp_sector: np.ndarray
    leaders: tuple
    failed_sectors: tuple
    converged: bool
    leading: complex
    selection: str


class IterationLeader(NamedTuple):
    """The power iteration of one pair-S_z sector.

    ===============  ===================  ==========================  =========
    status           value                gap                         n_iter
    ===============  ===================  ==========================  =========
    ``"empty"``      NaN                  ``None``                    0
    ``"null"``       0.0                  the first lifted unit       0
                                          vector ``e_0``
    ``"converged"``  the signed Rayleigh  the converged iterate       the steps
                     eigenvalue
    ``"failed"``     NaN                  the last iterate            the steps
    ===============  ===================  ==========================  =========

    A one-dimensional sector is solved exactly (``"converged"``, its single
    matrix element, ``n_iter = 0``; ``"failed"`` when that element is not
    finite). A run that reports convergence with a non-finite value is
    ``"failed"`` too. ``gap`` is lifted to the full gap shape."""
    pp_sector: int
    status: str
    value: object
    gap: object
    n_iter: int


@dataclass(frozen=True)
class SectorIteration:
    """The reduction of the sector power iterations.

    Without a failed sector the leader is the largest real part over the
    non-empty sectors (null sectors at 0; ties as in
    :func:`_so_solve_sectors`) and ``value``, ``gap``, ``pp_sector`` are its.
    With a failed sector ``value`` is NaN, ``converged`` false, and ``gap`` /
    ``pp_sector`` are the last iterate of the first failed sector in the
    order ``(+1, 0, -1)``. With every sector empty ``value`` is NaN, ``gap``
    ``None`` and ``pp_sector`` 2. ``n_iter`` is the largest sector step
    count; ``n_iter_per_sector`` lists them in sector order (0 for empty,
    null and exactly solved sectors). ``selection`` is ``"null-operator"``
    when the leader is a null sector, else ``"iteration"``."""
    value: object
    gap: object
    pp_sector: int
    leaders: tuple
    failed_sectors: tuple
    converged: bool
    n_iter: int
    n_iter_per_sector: tuple
    selection: str


class EigenLeg(NamedTuple):
    """The eigenvalue leg of a spin-orbital run, in either layout."""
    vals: np.ndarray
    gaps: np.ndarray
    pp_sector: np.ndarray
    converged: bool
    leading: complex
    selection: str
    failed_sectors: tuple

    @classmethod
    def from_dense(cls, vals, gaps, *, selection):
        """The result of :func:`_so_solve_eigenvalue` (the arrays are passed
        through, not copied); converged when the first eigenvalue is finite."""
        converged = bool(len(vals) > 0 and np.isfinite(vals[0]))
        leading = complex(vals[0]) if len(vals) else complex(np.nan)
        return cls(vals, gaps, np.full(len(vals), _SO_PP_DENSE, dtype=np.int8),
                   converged, leading, selection, ())

    @classmethod
    def from_sectors(cls, solve):
        return cls(solve.vals, solve.gaps, solve.pp_sector, solve.converged,
                   solve.leading, solve.selection, solve.failed_sectors)


class IterationLeg(NamedTuple):
    """The power-iteration leg of a spin-orbital run, in either layout."""
    value: object
    gap: object
    pp_sector: int
    converged: bool
    n_iter: int
    failed_sectors: tuple
    selection: str

    @classmethod
    def from_dense(cls, gap, value, converged, n_iter, info):
        """The result of :func:`_so_solve_iteration` (``info`` is its
        ``info_out`` record, accepted for the call's symmetry)."""
        return cls(value, gap, _SO_PP_DENSE, bool(converged), int(n_iter), (),
                   "iteration")

    @classmethod
    def from_sectors(cls, iteration):
        return cls(iteration.value, iteration.gap, iteration.pp_sector,
                   iteration.converged, iteration.n_iter,
                   iteration.failed_sectors, iteration.selection)


class SectorReport(NamedTuple):
    """What a spin-orbital run reports, one schema for both layouts.

    ``reported``: the reported eigenvalue (NaN when not certified).
    ``gap`` (``gap_shape``) with ``gap_status``: ``"certified"`` (the gap of
    the certified leader; exactly when ``reported`` is finite),
    ``"uncertified"`` (a real candidate -- a listed pair or the last iterate
    -- written while the result is not certified) or ``"none"`` (no candidate
    exists; a zero array). ``pp_sector``: the sector of the written gap (2
    when not sector-resolved or none). ``eigen_converged``;
    ``iteration_converged`` (only in ``"both"``); ``selection``; the listing
    ``vals``, ``gaps``, ``pp_sectors`` (empty in ``"iteration"``);
    ``failed_sectors`` (of the authoritative leg) and
    ``iteration_failed_sectors`` (the iteration leg's, whenever it ran)."""
    reported: float
    gap: np.ndarray
    gap_status: str
    pp_sector: int
    eigen_converged: bool
    iteration_converged: object
    selection: str
    vals: np.ndarray
    gaps: np.ndarray
    pp_sectors: np.ndarray
    failed_sectors: tuple
    iteration_failed_sectors: object


def _so_canonical_sector_ops(sector_ops, solver):
    """``sector_ops`` validated (each pp_sector +1, 0 or -1, none twice) and
    put in the order ``(+1, 0, -1)``."""
    seen = {}
    for op in sector_ops:
        p = op.pp_sector
        if (isinstance(p, (bool, np.bool_))
                or not isinstance(p, (int, np.integer)) or p not in _SO_PP_ORDER):
            raise ValueError("{}: pp_sector must be +1, 0 or -1, got {}"
                             .format(solver, p))
        if int(p) in seen:
            raise ValueError("{}: pp_sector {} appears twice".format(solver, p))
        seen[int(p)] = op
    return tuple(seen[p] for p in _SO_PP_ORDER if p in seen)


def _so_check_sector_op(op, gap_shape, solver):
    basis = op.basis
    if tuple(int(x) for x in basis.shape) != tuple(gap_shape):
        raise ValueError("{}: the basis of sector {} lives on a gap of shape {}, "
                         "not {}".format(solver, _pp_label(op.pp_sector),
                                         tuple(basis.shape), tuple(gap_shape)))
    if tuple(op.operator.shape) != (basis.n, basis.n):
        raise ValueError("{}: the operator of sector {} has shape {}, its basis "
                         "acts on vectors of length {}".format(
                             solver, _pp_label(op.pp_sector),
                             tuple(op.operator.shape), basis.n))


def _so_lifted_unit_gaps(basis, count, gap_shape):
    """The first ``count`` reduced unit vectors lifted to full gaps of
    ``gap_shape``, each normalized to unit 2-norm."""
    out = np.zeros((count,) + tuple(gap_shape), dtype=complex)
    for i in range(count):
        e = np.zeros(basis.dim, dtype=complex)
        e[i] = 1.0
        g = basis.scatter(basis.to_full(e))
        g /= np.linalg.norm(g)
        out[i] = g.reshape(gap_shape)
    return out


def _so_sector_norm_estimate(red, dim, adjoint):
    """Fixed-seed power estimate of ``||red||_2``.

    ``adjoint``: ``x_0`` from ``default_rng(0)``, normalized; ``y = K x_j``,
    ``z = K^dagger y``, ``mu_j = ||z||``, ``x_{j+1} = z / ||z||``; stop after
    30 steps or when ``|mu_j - mu_{j-1}| < 1e-3 mu_j``; the estimate is
    ``sqrt(mu_last)``. Without an adjoint: ``r_j = ||K x_j||`` over the
    normalized steps ``x_{j+1} = K x_j / ||K x_j||`` with the same stop, the
    estimate ``max_j r_j`` (a lower bound). A zero image at step 0 redraws
    ``x_0`` once; a second zero gives 0. Returns ``(estimate, steps,
    relative change, converged)``."""
    rng = np.random.default_rng(0)

    def draw():
        x = rng.standard_normal(dim) + 1j * rng.standard_normal(dim)
        return x / np.linalg.norm(x)

    def image(x):
        y = np.asarray(red.matvec(x)).ravel()
        if adjoint:
            y = np.asarray(red.rmatvec(y)).ravel()
        return y

    x = draw()
    z = image(x)
    if not np.any(z):
        x = draw()
        z = image(x)
        if not np.any(z):
            return 0.0, 1, 0.0, True
    mu_prev = None
    best = 0.0
    change = float("inf")
    steps = 0
    converged = False
    for j in range(_SO_NORM_STEPS):
        if j > 0:
            z = image(x)
        mu = float(np.linalg.norm(z))
        steps = j + 1
        best = max(best, mu)
        if mu == 0.0:
            break
        if mu_prev is not None:
            change = abs(mu - mu_prev) / mu
            if change < _SO_NORM_RTOL:
                converged = True
                break
        mu_prev = mu
        x = z / mu
    if adjoint:
        est = float(np.sqrt(mu))
    else:
        est = float(best)
    return est, steps, change, converged


def _so_union(cands, num_eigenvalues):
    """Order the candidate pairs ``(value, rank, in-sector index, gap,
    pp_sector)``: tie clusters by single linkage (``|l_i - l_j| <= 1e-12
    max(1, |l_i|, |l_j|)``) ordered by descending largest real part, then
    descending largest imaginary part, then the smallest sector rank of a
    member; inside a cluster by sector rank, then in-sector index. Returns
    the first ``num_eigenvalues``."""
    n = len(cands)
    parent = list(range(n))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    vals = [complex(c[0]) for c in cands]
    for i in range(n):
        for j in range(i + 1, n):
            tol = _SO_SECTOR_TIE_TOL * max(1.0, abs(vals[i]), abs(vals[j]))
            if abs(vals[i] - vals[j]) <= tol:
                parent[root(i)] = root(j)
    clusters = {}
    for i in range(n):
        clusters.setdefault(root(i), []).append(i)
    groups = sorted(
        clusters.values(),
        key=lambda m: (-max(vals[i].real for i in m),
                       -max(vals[i].imag for i in m),
                       min(cands[i][1] for i in m)))
    order = [i for m in groups
             for i in sorted(m, key=lambda i: (cands[i][1], cands[i][2]))]
    return [cands[i] for i in order[:num_eigenvalues]]


def _so_order_by_seed(kept, seed_unit):
    """Order the kept union pairs ``(value, rank, index, gap, pp_sector)``
    by descending normalized overlap ``|<seed, gap>| / ||gap||`` with the
    unit full seed, then descending real part, then sector rank (the
    sector form of :func:`_order_by_seed_overlap`). Logs a WARNING when
    the best overlap is below 0.5: the seed's branch is then not among the
    selected pairs."""
    ov = [abs(np.vdot(seed_unit, np.ravel(c[3])))
          / (float(np.linalg.norm(c[3])) + 1.0e-300) for c in kept]
    order = sorted(range(len(kept)),
                   key=lambda i: (-ov[i], -complex(kept[i][0]).real, kept[i][1]))
    best = ov[order[0]]
    if best < _SO_SEED_OVERLAP_WARN:
        logger.warning("spin-orbital: the largest overlap of a selected eigenvector with "
                       "the seed is %.3f (< %.1f); the seed's branch is not among the "
                       "selected eigenpairs.", best, _SO_SEED_OVERLAP_WARN)
    return [kept[i] for i in order]


def _so_solve_sectors(sector_ops, num_eigenvalues, method, *, gap_shape,
                      spectral_shift=None, seed_vec=None):
    """The eigen solve of the spin-orbital kernel sector by sector.

    ``sector_ops``: :class:`SectorOp` per pair-S_z sector (any subset of
    ``+1, 0, -1``, each at most once; taken in that order). Per sector, with
    ``k = num_eigenvalues`` and ``dim`` its antisymmetric dimension:

    * ``dim = 0``: ``"empty"``, no pairs, certified;
    * an exactly zero reduced operator (:func:`_operator_is_exactly_zero`):
      ``"null-operator"``, ``min(k, dim)`` eigenvalues 0 with the first
      reduced unit vectors as gaps (lifted, unit 2-norm), certified;
    * ``dim <= max(2k + 1, 20)``: dense diagonalization (``"dense-LR"``);
    * else :func:`_so_solve_eigenvalue` on the sector (ARPACK, ``"LM"`` /
      ``"LR_retry"`` / ``"LR"``).

    The seed ``seed_vec`` (a full gap, optional) is split per sector into
    reduced coordinates; a part with ``||seed_p||_2 <= 1e-8 ||seed||_2``
    means no seed for that sector (an INFO line when it is not exactly
    zero). With a nonzero seed the kept pairs are then ordered by their
    normalized overlap with the seed, then real part, then sector rank, as
    the dense layout orders a seeded solve (:func:`_so_order_by_seed`); the
    selection of the ``k`` pairs is unchanged. Every returned pair must
    satisfy ``||K v - lambda v|| <= 1e-10 n_K ||v||`` with ``n_K = max(1,
    |lambda|, est)``, ``est`` the operator-norm estimate
    (:func:`_so_sector_norm_estimate`; exact on the dense branch), else the
    whole sector is ``"residual-failed"``. A ``"not-converged"`` or
    ``"residual-failed"`` sector is failed: its pairs are not listed and the
    result is not converged (``leading = NaN``). The certified pairs are
    merged by :func:`_so_union` and the first ``k`` kept, their gaps placed
    in full gaps of ``gap_shape``. Only ``method = "arnoldi"`` is
    supported."""
    solver = "_so_solve_sectors"
    if method != "arnoldi":
        raise ValueError(
            "{}: the per-sector solve supports eigenvalue_method = \"arnoldi\" "
            "only, got {}".format(solver, method))
    ops = _so_canonical_sector_ops(sector_ops, solver)
    spectral_shift = _normalize_spectral_shift(spectral_shift)
    _validate_spectral_shift(spectral_shift, method)
    if (isinstance(num_eigenvalues, (bool, np.bool_))
            or not isinstance(num_eigenvalues, (int, np.integer))
            or num_eigenvalues < 1):
        raise ValueError("{}: num_eigenvalues must be a positive integer, got {}"
                         .format(solver, num_eigenvalues))
    k = int(num_eigenvalues)
    gap_shape = tuple(int(x) for x in gap_shape)
    for op in ops:
        _so_check_sector_op(op, gap_shape, solver)
    seed_flat = None if seed_vec is None else np.asarray(seed_vec).ravel()
    if seed_flat is not None and seed_flat.size != int(np.prod(gap_shape)):
        raise ValueError("{}: seed_vec has {} entries, the gap {} has {}".format(
            solver, seed_flat.size, gap_shape, int(np.prod(gap_shape))))
    if seed_flat is not None and not np.all(np.isfinite(seed_flat)):
        raise ValueError("{}: seed_vec has non-finite entries".format(solver))
    seed_norm = 0.0 if seed_flat is None else float(np.linalg.norm(seed_flat))
    if seed_norm == 0.0:
        seed_flat = None                    # a zero seed is no seed

    leaders = []
    cands = []          # (value, rank, index, gap, pp_sector)
    failed = []
    for op in ops:
        p = int(op.pp_sector)
        rank = _SO_PP_ORDER.index(p)
        basis = op.basis
        dim = int(basis.dim)
        if dim == 0:
            leaders.append(SectorLeader(p, 0, complex(np.nan), "empty", True,
                                        0.0, 0.0, "exact", 0))
            continue
        red = basis.reduce(op.operator, adjoint=bool(op.has_adjoint))
        n_ret = min(k, dim)
        if _operator_is_exactly_zero(red, dim):
            gaps = _so_lifted_unit_gaps(basis, n_ret, gap_shape)
            for i in range(n_ret):
                cands.append((0j, rank, i, gaps[i], p))
            leaders.append(SectorLeader(p, dim, 0j, "null-operator", True, 0.0,
                                        0.0, "exact", n_ret))
            continue
        seed = None
        if seed_flat is not None:
            part = basis.to_reduced(seed_flat[basis.indices]
                                    if basis.indices is not None else seed_flat)
            part_norm = float(np.linalg.norm(part))
            if part_norm > _SO_SEED_PART_RTOL * seed_norm:
                seed = part
            elif part_norm > 0.0:
                logger.info("spin-orbital sector %s: the seed part (norm %.3e, %.3e of "
                            "the seed) is round-off; the sector is solved unseeded.",
                            _pp_label(p), part_norm, part_norm / seed_norm)
        M = None
        if dim <= max(2 * k + 1, 20):
            try:
                vals, vecs, M = _so_dense_reduced_eigen(red, dim, k, seed)
                selection = "dense-LR"
            except np.linalg.LinAlgError as exc:
                logger.warning("spin-orbital sector %s: the dense "
                               "diagonalization failed (%s).", _pp_label(p), exc)
                vals = np.full(n_ret, np.nan + 0j)
                vecs = np.zeros((dim, n_ret), dtype=complex)
                selection = "not-converged"
        else:
            info = {}
            vals, vecs = _so_solve_eigenvalue(
                op.operator, basis, k, method, spectral_shift=spectral_shift,
                seed_vec=seed, map_back=False, info_out=info)
            selection = info.get("selection")
            if selection == "not-converged" or not np.all(np.isfinite(vals)):
                selection = "not-converged"
        vals = np.asarray(vals)
        vecs = np.asarray(vecs)
        kind = ("exact" if dim <= max(2 * k + 1, 20)
                else "adjoint" if op.has_adjoint else "forward")
        if selection == "not-converged":
            logger.warning("spin-orbital sector %s: the eigen solve did not "
                           "converge; the sector is failed.", _pp_label(p))
            leaders.append(SectorLeader(p, dim, complex(np.nan), selection,
                                        False, float("inf"), float("nan"),
                                        kind, len(vals)))
            failed.append(p)
            continue
        if kind == "exact":
            est = float(np.linalg.norm(M, 2))
            images = M @ vecs
        else:
            est, steps, change, conv = _so_sector_norm_estimate(
                red, dim, op.has_adjoint)
            images = np.asarray(red.matmat(vecs))
        scale = max(1.0, float(np.max(np.abs(vals))), est)
        if kind == "forward" or (kind == "adjoint" and not conv):
            logger.info(
                "spin-orbital sector %s: operator-norm estimate %s after %d "
                "steps (relative change %.2e); residual scale %.6e",
                _pp_label(p), kind, steps, change, scale)
        resid = 0.0
        for i in range(len(vals)):
            v = vecs[:, i]
            n_k = max(1.0, abs(vals[i]), est)
            r = (float(np.linalg.norm(images[:, i] - vals[i] * v))
                 / (n_k * float(np.linalg.norm(v))))
            resid = max(resid, r) if np.isfinite(r) else float("inf")
        certified = bool(resid <= _SO_SECTOR_RESIDUAL_TOL)
        lead_value = complex(vals[np.argmax(vals.real)]) if len(vals) else complex(np.nan)
        if not certified:
            logger.warning(
                "spin-orbital sector %s: a returned eigenpair fails the "
                "residual check (%.3e > %.0e relative to the scale %.6e); "
                "the sector is failed.", _pp_label(p), resid,
                _SO_SECTOR_RESIDUAL_TOL, scale)
            leaders.append(SectorLeader(p, dim, lead_value, "residual-failed",
                                        False, resid, est, kind, len(vals)))
            failed.append(p)
            continue
        leaders.append(SectorLeader(p, dim, lead_value, selection, True, resid,
                                    est, kind, len(vals)))
        for i in range(len(vals)):
            g = basis.scatter(basis.to_full(vecs[:, i])).reshape(gap_shape)
            cands.append((complex(vals[i]), rank, i, g, p))

    if all(lead.dim == 0 for lead in leaders):
        logger.warning("spin-orbital: the antisymmetric gap space is empty on "
                       "this grid; no eigenvalue is reported")
    kept = _so_union(cands, k)
    if seed_flat is not None and kept:
        kept = _so_order_by_seed(kept, seed_flat / seed_norm)
    vals = np.array([c[0] for c in kept], dtype=complex)
    gaps = (np.array([c[3] for c in kept], dtype=complex) if kept
            else np.zeros((0,) + gap_shape, dtype=complex))
    pp = np.array([c[4] for c in kept], dtype=np.int8)
    converged = bool(not failed and len(kept) > 0)
    leading = complex(vals[0]) if converged else complex(np.nan)
    selection = ";".join("{}:{}".format(_pp_label(lead.pp_sector), lead.selection)
                         for lead in leaders)
    return SectorSolve(vals, gaps, pp, tuple(leaders), tuple(failed), converged,
                       leading, selection)


def _so_iterate_sectors(sector_ops, sigma_init, *, gap_shape, max_iter, alpha,
                        tol, spectral_shift=None):
    """The power iteration of the spin-orbital kernel sector by sector.

    Per sector (order ``+1, 0, -1``): ``dim = 0`` is skipped (``"empty"``);
    an exactly zero reduced operator is ``"null"`` (``lambda = 0``, the first
    lifted unit vector, not iterated); ``dim = 1`` is solved exactly;
    otherwise ``_solve_leading(..., "iteration")`` runs on the reduced sector
    operator from the sector's part of ``sigma_init`` (a zero part starts
    from the complex ``default_rng(0)`` vector of the reduced space), and a
    sector whose run did not converge is failed. The reduction is
    :class:`SectorIteration`'s."""
    solver = "_so_iterate_sectors"
    ops = _so_canonical_sector_ops(sector_ops, solver)
    spectral_shift = _normalize_spectral_shift(spectral_shift)
    _validate_spectral_shift(spectral_shift, "iteration")
    gap_shape = tuple(int(x) for x in gap_shape)
    for op in ops:
        _so_check_sector_op(op, gap_shape, solver)
    init_flat = np.asarray(sigma_init).ravel()
    if init_flat.size != int(np.prod(gap_shape)):
        raise ValueError("{}: sigma_init has {} entries, the gap {} has {}".format(
            solver, init_flat.size, gap_shape, int(np.prod(gap_shape))))

    leaders = []
    for op in ops:
        p = int(op.pp_sector)
        basis = op.basis
        dim = int(basis.dim)
        if dim == 0:
            leaders.append(IterationLeader(p, "empty", float("nan"), None, 0))
            continue
        red = basis.reduce(op.operator)
        if _operator_is_exactly_zero(red, dim):
            gap = _so_lifted_unit_gaps(basis, 1, gap_shape)[0]
            leaders.append(IterationLeader(p, "null", 0.0, gap, 0))
            continue
        if dim == 1:
            # a one-dimensional sector: its single element is the eigenvalue
            lam = complex(np.asarray(red.matvec(np.ones(1, dtype=complex))).ravel()[0])
            gap = _so_lifted_unit_gaps(basis, 1, gap_shape)[0]
            if np.isfinite(lam):
                leaders.append(IterationLeader(p, "converged", lam, gap, 0))
            else:
                logger.warning("spin-orbital sector %s: the matrix element of the "
                               "one-dimensional sector is not finite (%s); the sector "
                               "is failed.", _pp_label(p), lam)
                leaders.append(IterationLeader(p, "failed", float("nan"), gap, 0))
            continue
        part = basis.to_reduced(init_flat[basis.indices]
                                if basis.indices is not None else init_flat)
        if not np.any(part != 0):
            rng = np.random.default_rng(0)
            part = rng.standard_normal(dim) + 1j * rng.standard_normal(dim)
        value, vec, info = _solve_leading(
            lambda: (red, dim), dim, "iteration", max_iter=max_iter,
            convergence_tol=tol, init_vec=part, alpha=alpha,
            spectral_shift=spectral_shift)
        gap = basis.scatter(basis.to_full(np.asarray(vec).ravel())).reshape(gap_shape)
        n_iter = int(info.get("n_iter", 0))
        if info.get("converged") and np.isfinite(value):
            leaders.append(IterationLeader(p, "converged", value, gap, n_iter))
        else:
            leaders.append(IterationLeader(p, "failed", float("nan"), gap, n_iter))

    leaders = tuple(leaders)
    per_sector = tuple(lead.n_iter for lead in leaders)
    n_iter = max(per_sector) if per_sector else 0
    failed = tuple(lead.pp_sector for lead in leaders if lead.status == "failed")
    if failed:
        first = next(lead for lead in leaders if lead.status == "failed")
        return SectorIteration(float("nan"), first.gap, first.pp_sector, leaders,
                               failed, False, n_iter, per_sector, "iteration")
    live = [lead for lead in leaders if lead.status != "empty"]
    if not live:
        # the empty-space WARNING is the eigen leg's (one per run)
        return SectorIteration(float("nan"), None, _SO_PP_DENSE, leaders, (),
                               False, n_iter, per_sector, "iteration")
    best = _so_union([(complex(lead.value), _SO_PP_ORDER.index(lead.pp_sector), 0,
                       None, lead.pp_sector) for lead in live], 1)[0]
    lead = next(x for x in live if x.pp_sector == best[4])
    selection = "null-operator" if lead.status == "null" else "iteration"
    return SectorIteration(lead.value, lead.gap, lead.pp_sector, leaders, (),
                           True, n_iter, per_sector, selection)


def _so_bond_report(solver_mode, eigen, iteration, *, gap_shape):
    """The one report of a spin-orbital run (:class:`SectorReport`) from its
    legs (:class:`EigenLeg` for ``"eigenvalue"`` and ``"both"``,
    :class:`IterationLeg` for ``"iteration"`` and ``"both"``; the other
    ``None``).

    * ``"eigenvalue"``: the eigen leg's leader (``reported = Re leading``,
      NaN if not converged); the gap of the leader, else of the first listed
      pair (``"uncertified"``), else zeros (``"none"``).
    * ``"iteration"``: the iteration value and gap (``eigen_converged`` is the
      iteration's convergence, which needs a finite value); zeros and
      ``"none"`` when it has no gap.
    * ``"both"``: the eigen leg decides ``reported``, ``eigen_converged``,
      the listing and ``failed_sectors``; the gap is the eigen leader's when
      the eigen leg converged, else the iteration's (``"uncertified"``),
      else zeros; ``iteration_converged`` is the iteration's convergence
      (again with a finite value).

    A listing whose first value is NaN (the failed-solve sentinel of the
    dense layout) holds no candidate: in ``"eigenvalue"`` mode its gap is
    zeros with ``"none"``.

    ``gap_status`` is ``"certified"`` exactly when ``reported`` is finite."""
    if solver_mode not in ("eigenvalue", "iteration", "both"):
        raise ValueError("_so_bond_report: solver_mode must be \"eigenvalue\", "
                         "\"iteration\" or \"both\", got {!r}".format(solver_mode))
    if solver_mode in ("eigenvalue", "both") and eigen is None:
        raise ValueError("_so_bond_report: solver_mode {!r} needs the eigenvalue "
                         "leg".format(solver_mode))
    if solver_mode in ("iteration", "both") and iteration is None:
        raise ValueError("_so_bond_report: solver_mode {!r} needs the iteration "
                         "leg".format(solver_mode))
    gap_shape = tuple(int(x) for x in gap_shape)

    def zeros():
        return np.zeros(gap_shape, dtype=complex)

    def shaped(g):
        g = np.asarray(g)
        return g if g.shape == gap_shape else g.reshape(gap_shape)

    nan = float("nan")
    iteration_converged = None
    iteration_failed = None
    def iteration_ok():
        # convergence of the iteration leg needs a finite value
        return bool(iteration.converged) and bool(np.isfinite(iteration.value))

    if solver_mode == "iteration":
        converged = iteration_ok()
        reported = float(np.real(iteration.value)) if converged else nan
        if iteration.gap is None:
            gap, pp = zeros(), _SO_PP_DENSE
        else:
            gap, pp = shaped(iteration.gap), int(iteration.pp_sector)
        selection = iteration.selection
        vals = np.zeros(0, dtype=complex)
        gaps = np.zeros((0,) + gap_shape, dtype=complex)
        pps = np.zeros(0, dtype=np.int8)
        failed = tuple(iteration.failed_sectors)
        iteration_failed = failed
        has_candidate = iteration.gap is not None
    else:
        converged = bool(eigen.converged)
        reported = float(np.real(eigen.leading)) if converged else nan
        if not np.isfinite(reported):
            converged = False
        n = len(eigen.vals)
        if converged:
            gap, pp = shaped(eigen.gaps[0]), int(eigen.pp_sector[0])
            has_candidate = True
        elif solver_mode == "both":
            has_candidate = iteration.gap is not None
            gap, pp = ((shaped(iteration.gap), int(iteration.pp_sector))
                       if has_candidate else (zeros(), _SO_PP_DENSE))
        else:
            # a NaN-valued listing (the failed-solve sentinel of the dense
            # layout) holds placeholders, not candidates
            has_candidate = n > 0 and bool(np.isfinite(eigen.vals[0]))
            gap, pp = ((shaped(eigen.gaps[0]), int(eigen.pp_sector[0]))
                       if has_candidate else (zeros(), _SO_PP_DENSE))
        if solver_mode == "both":
            iteration_converged = iteration_ok()
            iteration_failed = tuple(iteration.failed_sectors)
        selection = eigen.selection
        vals, gaps, pps = eigen.vals, eigen.gaps, eigen.pp_sector
        failed = tuple(eigen.failed_sectors)
    if np.isfinite(reported):
        status = "certified"
    elif has_candidate:
        status = "uncertified"
    else:
        status = "none"
    return SectorReport(reported, gap, status, pp,
                        bool(converged), iteration_converged, selection, vals,
                        gaps, pps, failed, iteration_failed)


def _so_check_antisymmetric_eigenpairs(vals, match):
    """Defensive check of the spin-orbital eigenpairs (``match`` from
    :func:`_so_is_antisymmetric`): refuse when the LEADING eigenvector --
    the gap written under the fermion-antisymmetric header -- is not
    antisymmetric, and warn about any further pair that is not (match
    column 0). The solvers work on antisymmetric coordinates, so this fires
    only if that construction is broken."""
    if np.size(vals) and not np.all(np.isfinite(np.asarray(vals))):
        # the eigenvalue solver could not determine the largest real part
        # (NaN sentinel, _solve_leading): the vectors are placeholders, not
        # eigenvectors, and the run is reported with eigen_converged = false
        logger.warning(
            "spin-orbital: no eigenvalue was determined (NaN); the written gap "
            "is a placeholder, not an eigenvector.")
        return
    match = np.asarray(match, dtype=bool)
    if match.size and not match[0]:
        raise ValueError(
            "[mode] enable_spin_orbital = true: the leading eigenvector is not "
            "fermion-antisymmetric (Delta_AB(k) = -Delta_BA(-k)) to round-off; "
            "it is Pauli-forbidden and is not written (eigenvalues {})."
            .format(np.asarray(vals)))
    bad = np.where(~match)[0]
    if bad.size:
        logger.warning(
            "spin-orbital eigenpairs %s are not fermion-antisymmetric to "
            "round-off (eigenvalues %s); they are written with match = 0.",
            bad.tolist(), np.asarray(vals)[bad])


def calc_eliashberg(input_dict):
    """Main calculation orchestration for linearized Eliashberg equation.

    Supports ``[eliashberg] gpu = true`` for both ``frequency = "static"`` and
    ``"dynamic"``: the kernel matvec/matmat runs on the GPU (CuPy) while the
    host eigensolvers (ARPACK/power/subspace/shift-invert) consume host arrays.
    Falls back to CPU with a warning when CuPy/CUDA is unusable, unless
    ``gpu_required = true`` (then it fails fast).

    Parameters
    ----------
    input_dict : dict
        Parsed TOML configuration dictionary.

    Returns
    -------
    float or None
        The reported leading eigenvalue lambda: the eigenvalue analysis' for
        ``solver_mode = "eigenvalue"`` (the default) and ``"both"``, the power
        iteration's for ``"iteration"``; ``NaN`` when it could not be
        determined (``eigen_converged = false`` in the outputs; the
        ``hwave_sc`` command then exits with status 3). The dynamic and
        spin-orbital entries return the same quantity. (Earlier versions
        returned None from the static path.)
    """
    # --- Parse parameters ---
    # [mode] enable_spin_orbital (issue #83): the static entry solves the
    # spin-orbital gap equation; every combination it does not support is
    # refused here, before any file is read (check_so_config), and after the
    # interaction files are read (check_so_supported, below).
    spin_orbital_mode = _resolve_spin_orbital_flag(input_dict)
    if spin_orbital_mode:
        from hwave.solver import so_pairing
        so_pairing.check_so_config(input_dict)
        if so_pairing.so_frequency(input_dict) == "dynamic":
            # the spin-orbital dynamic solver (issue #284); [eliashberg] is
            # read case-insensitively there, as the checks above read it
            from hwave.solver import eliashberg_dynamic_so
            return eliashberg_dynamic_so.solve_dynamic_so(input_dict)
    mode_param = input_dict["mode"]["param"]
    T = mode_param["T"]
    beta = _coerce_run_beta(T)
    cell_shape = mode_param["CellShape"]
    # Resolve SubShape by the PACKAGE convention (documented default:
    # CellShape, i.e. the whole cell as one supercell) and guard the
    # RESOLVED value: a guard keyed on the explicit key alone let an
    # omitted SubShape default to a fully folded configuration and reach
    # the file reader with the very mismatch the guard exists to stop
    # (round-7 review).
    if isinstance(cell_shape, (list, tuple)):
        _cs = list(cell_shape)
    else:
        _cs = [cell_shape]
    while len(_cs) < 3:
        _cs.append(1)
    _ss = list(mode_param.get("SubShape", _cs))
    while len(_ss) < 3:
        _ss.append(1)
    if _ss != [1, 1, 1]:
        # supported nowhere in this module: the geometry and
        # interactions are consumed UNFOLDED here, so a folded
        # susceptibility mismatches the expected orbital count and an
        # off-site bond would fold onto an on-site supercell entry
        # (round-4 review); failing late produced an unhelpful shape
        # error instead of this actionable one
        raise ValueError(
            "SubShape (sublattice folding) is not supported by the "
            "Eliashberg module: fold the model into the unit cell "
            "yourself, or set SubShape = [1, 1, 1] explicitly. (Note: "
            "omitting SubShape defaults it to CellShape -- the whole "
            "cell as one supercell -- per the package convention, so it "
            "must be set explicitly here.)")
    sub_shape = _ss
    nmat = mode_param.get("Nmat", _DEFAULT_NMAT)

    # Filling
    if "filling" in mode_param:
        n_filling = mode_param["filling"]
    elif "Ncond" in mode_param:
        raise ValueError("Ncond not yet supported in eliashberg. Use filling.")
    else:
        raise ValueError("filling must be specified.")

    # Lattice dimensions
    if isinstance(cell_shape, list):
        while len(cell_shape) < 3:
            cell_shape.append(1)
    Lx, Ly, Lz = cell_shape
    if isinstance(sub_shape, list):
        while len(sub_shape) < 3:
            sub_shape.append(1)
    Bx, By, Bz = sub_shape
    Nx, Ny, Nz = Lx // Bx, Ly // By, Lz // Bz
    nvol = Nx * Ny * Nz

    # Eliashberg parameters
    eli_param = input_dict.get("eliashberg", {})
    if spin_orbital_mode:
        # the spin-orbital checks read [eliashberg] case-insensitively, so the
        # run must too (a key the checks see must be the key the run uses)
        from requests.structures import CaseInsensitiveDict
        eli_param = CaseInsensitiveDict(eli_param or {})

    # Dispatch to dynamic Eliashberg if requested
    if _eliashberg_frequency(input_dict) == "dynamic":
        _validate_dynamic_prereqs(input_dict)
        _warn_one_body_band_ignored(input_dict, 'frequency = "dynamic"')
        if _bond_bool_option(eli_param, "bond_channels", False):
            # bond-resolved dynamic pairing kernel, fed by the FLEX run's
            # bond archive (the on-site solver below reads chiq_s/chiq_c)
            from hwave.solver import eliashberg_bond_io
            return eliashberg_bond_io.solve_dynamic_bond(input_dict)
        from hwave.solver import eliashberg_dynamic
        return eliashberg_dynamic.solve_dynamic(input_dict)

    # Static path from here on: the channel-decomposition flags only affect the
    # dynamic vertex, so warn rather than silently ignore them.
    _warn_if_static_ignores_channel_flags(eli_param)

    # GPU backend for the static kernel (the dynamic branch returned above and
    # handles its own gpu flag). get_backend falls back to numpy with a warning
    # when CuPy/CUDA is unusable, unless gpu_required=true (then it raises).
    use_gpu = backend.as_bool(eli_param.get("gpu", False))
    gpu_required = backend.as_bool(eli_param.get("gpu_required", False))
    xp, gpu_active = backend.get_backend(use_gpu, logger, required=gpu_required)

    solver_mode = eli_param.get("solver_mode", DEFAULT_SOLVER_MODE)
    max_iter = eli_param.get("max_iter", 1000)
    alpha = eli_param.get("alpha", 0.5)
    tol = eli_param.get("convergence_tol", 1.0e-5)
    num_eigenvalues = eli_param.get("num_eigenvalues", 10)
    eigenvalue_method = eli_param.get("eigenvalue_method", "arnoldi")
    pairing_type = eli_param.get("pairing_type", "singlet")
    # Default the initial gap to the channel's parity (even for singlet, odd
    # for triplet) when the user did not specify one. For a centrosymmetric
    # system the kernel preserves parity in exact arithmetic, but numerical
    # noise leaks the wrong parity; _solve_iteration projects every iterate back
    # onto the channel's sector (pairing_type below), so the seed stays in -- and
    # converges within -- the physical sector. (For a non-centrosymmetric kernel
    # _solve_iteration detects [A,P] != 0 and disables the projection.)
    init_gap_mode = _resolve_init_gap(eli_param.get("init_gap"), pairing_type)
    chi0q_mode = eli_param.get("chi0q_mode", "load")
    chi0q_tensor = eli_param.get("chi0q_tensor", "auto")
    # Bond-resolved interaction channels (opt-in, spec S5). The options are
    # ignored with a warning when the flag is off; the guards themselves run
    # right after the interaction files are read (below), before any
    # chi0q/FLEX data is touched.
    (use_bond_channels, bond_green, bond_max_shells, bond_memory_cap_gb,
     bond_precondition_opts, bond_diagnostics) = _read_bond_config(eli_param)
    from requests.structures import CaseInsensitiveDict
    g2_tail = _coerce_g2_tail(
        CaseInsensitiveDict(eli_param).get("g2_tail", True))
    gap_file = eli_param.get("output_gap", "gap.dat")
    eigenvalue_file = eli_param.get("output_eigenvalue", "eigenvalue.dat")

    output_dir = input_dict["file"]["output"]["path_to_output"]
    if spin_orbital_mode:
        so_npz = _so_gap_npz_name(gap_file)
        if os.path.normpath(so_npz) == os.path.normpath(eigenvalue_file):
            raise ValueError(
                "[mode] enable_spin_orbital = true (hwave_sc, issue #83): the "
                "spin-orbital gap archive '{}' (derived from [eliashberg] "
                "output_gap = {!r}) would overwrite [eliashberg] "
                "output_eigenvalue = {!r}. Change output_gap or "
                "output_eigenvalue.".format(so_npz, gap_file, eigenvalue_file))

    logger.info("=== Linearized Eliashberg Equation Solver ===")
    logger.info("T = {}, beta = {:.4f}".format(T, beta))
    logger.info("Cell: {}x{}x{}, Sub: {}x{}x{}, Grid: {}x{}x{}".format(
        Lx, Ly, Lz, Bx, By, Bz, Nx, Ny, Nz))
    logger.info("Nmat = {}, filling = {}".format(nmat, n_filling))
    logger.info("solver_mode = {}, chi0q_mode = {}".format(solver_mode, chi0q_mode))

    # --- Step 2: Read input files ---
    geom_info, hr, interactions = _read_interaction_files(input_dict)
    norb = geom_info["norb"]
    if spin_orbital_mode:
        # Geometry norb and the Transfer file count SPIN-ORBITALS (interleaved
        # 2*orb+spin); the interaction files use PHYSICAL orbital indices. The
        # gap, the Green function and the vertex run over nd = norb spin-
        # orbitals in spin-block order (spin_orbital.INDEX_CONVENTION).
        from hwave.solver import spin_orbital
        norb_phys = spin_orbital.physical_norb(
            norb, True, source="the Geometry file ({})".format(
                input_dict["file"]["input"]["interaction"].get("Geometry", "")))
        # every R = 0 row is treated exactly; [mode.param] orbital_site is
        # never parsed (one WARNING when given, from the entry context)
        so_pairing.check_so_supported(input_dict, interactions, norb_phys=norb_phys)
        # the entry context: the interaction compiled once for the whole
        # run (the nested RPA runs reuse it)
        so_context = spin_orbital.build_so_entry_context(
            input_dict["mode"].get("param"), interactions, norb_phys, logger)
        # the Zeeman field sigma^z H_ab(R) of an Extern file (issue #88),
        # read and checked in every chi0q_mode before any heavy work
        so_extern_table, _so_extern_file = _read_so_extern(input_dict, norb_phys)
        logger.info("enable_spin_orbital: %d physical orbitals, nd = %d "
                    "spin-orbitals (spin-block order)", norb_phys, norb)
        # this run's band and state, for the archive checks (#88)
        so_run = _so_run_inputs(input_dict, hr, so_extern_table, norb_phys,
                                beta, nmat, n_filling, context=so_context)

    if use_bond_channels:
        # Every guard raises before any susceptibility/Green data is built, so
        # a misconfigured bond run fails fast and never emits numbers from a
        # path it is not validated for (spec S5).
        _validate_bond_prereqs(chi0q_mode, norb, interactions)
        logger.info(
            "bond_channels=true: using the bond-resolved (Z+1)x(Z+1) "
            "interaction channels for the static Eliashberg vertex.")
        logger.warning(
            "[eliashberg] bond_channels=true is a STATIC RPA-ladder bond "
            "dressing on %s. It demonstrates the qualitative V-dependence of "
            "the pairing eigenvalue; it is NOT a conserving FLEX result "
            "(the bond chi-bar is built here at RPA-ladder level), so "
            "absolute lambda values are not comparable with FLEX/dynamic "
            "references.",
            ("the EXTERNAL Green function '{}' ([eliashberg] bond_green)"
             .format(bond_green) if bond_green else
             "the Green function of the mean-field one-body band ({})"
             .format(" / ".join(_one_body_band_sources(input_dict)))
             if _one_body_band_sources(input_dict) else
             "the BARE Green function built from the transfer Hamiltonian "
             "(set [eliashberg] bond_green to feed an external/FLEX green)"))
        logger.warning(
            "[eliashberg] bond_channels=true builds its own bond-resolved "
            "bubble directly from the Green function (the m=m'=0 block IS the "
            "ordinary chi0q), so chi0q_mode='%s' and chi0q_tensor are not "
            "used: no chi0q file is read or computed on this path.",
            chi0q_mode)
        if use_gpu:
            logger.warning(
                "[eliashberg] gpu=true is ignored for bond_channels=true: the "
                "bond-resolved kernel is CPU-only in v1 (GPU support is a "
                "deferred follow-up).")

    if not spin_orbital_mode and not use_bond_channels:
        # a spin-orbital archive is refused before any heavy work (#271)
        _preflight_refuse_so_archives(input_dict, chi0q_mode)

    # --- One-body band (issue #273) ---
    # trans_mod / green_init / an active Extern replace the bare Transfer
    # band. The internal chi0q has always been built on that band; mu and
    # every kernel Green function must come from the SAME band, or the run
    # mixes two Hamiltonians. Built (and a spin-dependent band refused)
    # before any chi0q is computed or loaded. The FLEX path takes chi and G
    # from the FLEX run, whose outputs already include its one-body band.
    band = None
    spin_free_tol = _read_spin_free_tol(eli_param)
    if spin_orbital_mode:
        # spin-orbital mode refuses trans_mod / green_init; an Extern field
        # enters its band directly (_so_epsilon_k, Step 4), as it enters the
        # spin-orbital RPA and FLEX bands, so nothing here applies to it
        pass
    elif chi0q_mode == "flex":
        _warn_one_body_band_ignored(input_dict, "chi0q_mode = 'flex'")
    elif use_bond_channels and bond_green is not None:
        # the external Green function replaces the band entirely: the bond
        # path builds its bubble and kernel from bond_green
        _warn_one_body_band_ignored(
            input_dict, "bond_channels with an external bond_green",
            "the Green function of the pairing kernel and of the bond bubble "
            "is read from bond_green = '{}'".format(bond_green))
    else:
        band_sources = _one_body_band_sources(input_dict)
        if band_sources:
            if chi0q_mode == "calc" and not use_bond_channels:
                # the same solver computes chi0q below
                band_scheme = _chi0q_tensor_scheme(
                    chi0q_tensor,
                    input_dict["file"]["input"].get("interaction", {}))
            else:
                # band only; 'general' keeps the RPA reader from issuing its
                # reduced-scheme advisory for a chi0q that is never built
                band_scheme = "general"
            # spin averaging and hermitization need hwave_sc to build chi0q
            # itself (calc, or the bond path's own bubble); a loaded chi0q
            # was built by an external RPA run on the raw band
            own_chi0q = (chi0q_mode == "calc") or use_bond_channels
            band = _mean_field_band(input_dict, band_sources, band_scheme,
                                    (Nx, Ny, Nz), norb,
                                    spin_free_tol=spin_free_tol,
                                    allow_average=own_chi0q,
                                    allow_hermitize=own_chi0q)
            if chi0q_mode != "calc" or use_bond_channels:
                band = band._replace(solver=None)

    # --- Step 3: Setup k-mesh ---
    kx_array = np.linspace(0, 2.0 * np.pi, Nx, endpoint=False)
    ky_array = np.linspace(0, 2.0 * np.pi, Ny, endpoint=False)
    kz_array = np.linspace(0, 2.0 * np.pi, Nz, endpoint=False)

    # --- Step 4: Build Hamiltonian and diagonalize ---
    if spin_orbital_mode:
        # the Transfer file is interleaved; everything downstream is
        # spin-block. Extern adds coeff_extern * sigma^z H_ab(k), exactly as
        # in the spin-orbital RPA band (spin-orbital mode refuses trans_mod /
        # green_init, so it never takes the mean-field band below)
        logger.info("Building Hamiltonian in k-space...")
        epsilon_k = _so_epsilon_k(input_dict, kx_array, ky_array, kz_array,
                                  hr, so_extern_table, norb_phys)
    elif band is None or band.same_as_bare:
        if band is not None:
            logger.info("The one-body band from {} equals the Transfer band "
                        "exactly".format(", ".join(band.sources)))
        logger.info("Building Hamiltonian in k-space...")
        epsilon_k = _build_hamiltonian_k(kx_array, ky_array, kz_array, hr, norb)
    else:
        logger.info("Using the mean-field one-body band ({}) as the "
                    "Hamiltonian in k-space".format(", ".join(band.sources)))
        epsilon_k = band.epsilon_k

    logger.info("Diagonalizing Hamiltonian...")
    eigenvalues, eigenvectors = _calc_eigenvalues(epsilon_k)

    # --- Step 5: Determine chemical potential ---
    logger.info("Determining chemical potential...")
    mu = _determine_mu(eigenvalues, beta, n_filling, norb)
    logger.info("mu = {:.6f}".format(mu))

    # --- Step 7: Build interaction in k-space ---
    logger.info("Building interactions in k-space...")
    # only the non-bond paths read the k-space table; the bond path builds its
    # vertices from the declarations themselves
    inter_k = None
    if not use_bond_channels and not spin_orbital_mode:
        _check_nonbond_offsite(interactions, "hwave_sc (bond_channels = false)",
                               _nonbond_site(input_dict, norb))
        inter_k = _build_interaction_k_for_sites(kx_array, ky_array, kz_array, interactions,
                                                 norb, input_dict)

    # Pre-built Eliashberg operator + provenance record; only the bond path
    # sets them (every other path builds its operator from Vs_q/G2 as before).
    bond_operator = None
    bond_provenance = None
    bond_attribution = None

    if spin_orbital_mode:
        # --- Spin-orbital mode (issue #83): one spinful channel ---
        if chi0q_mode == "flex":
            # the stored chi and G of a spin-orbital FLEX run (tags checked
            # before any heavy work); Gamma from the interaction files, checked
            # against the run's through its (required) chi0q archive
            chi_so, green_kw, chi0_so, chi_path, so_borrowed = _load_flex_so(
                input_dict, norb, Nx, Ny, Nz, interactions=interactions,
                so_run=so_run)
            logger.info("Using the spin-orbital FLEX chiq and dressed Green's "
                        "function")
            so_members = _so_flex_writer_members(
                "sc_so_static", _SO_STATIC_LABEL, so_run, so_borrowed, green_kw,
                epsilon_k, beta)
            gam = _so_rpa_vertex(input_dict, so_context=so_run.context)
            _check_so_flex_vertex(gam, chi_so, chi0_so, chi_path)
            del chi0_so
            Vs_q = _so_static_vertex_from_chi(gam, chi_so, norb, (Nx, Ny, Nz))
            del chi_so, gam
        else:
            if chi0q_mode == "calc":
                # only the zero bosonic frequency is read (issue #272)
                chi0q_raw = _calc_chi0q_internal(input_dict, chi0q_tensor="general",
                                                 precomputed_mu=mu,
                                                 static_only=True,
                                                 so_context=so_run.context)
                static_index = 0
            else:
                from hwave.solver.fingerprint import SOArchiveError
                try:
                    chi0q_raw, static_index, _so_info = _load_chi0q(
                        input_dict, norb=norb, so_run=so_run)
                except SOArchiveError:
                    raise       # a provenance refusal is not a layout error
                except ValueError as err:
                    raise ValueError(
                        "{} ([mode] enable_spin_orbital = true with chi0q_mode = "
                        "\"load\" needs the general (4-index) spinful chi0q of a "
                        "spin-orbital RPA run over nd = {} spin-orbitals; a model the "
                        "RPA reduces to spin-free or spin-diagonal form writes a "
                        "physical-orbital chi0q instead -- use chi0q_mode = \"calc\" "
                        "for it.)".format(err, norb)) from err
                # refuse a reduced / non-spinful layout before any further work
                # (_so_static_vertex re-checks it)
                _so_static_chi0(chi0q_raw, norb, Nx, Ny, Nz, static_index)
            logger.info("Calculating Green's function G(k, iwn) (spin-orbital)...")
            green_kw = _build_bond_green(
                eigenvalues, eigenvectors, mu, beta, nmat, 0.0, None).full_sc
            # the kernel's own G and state (in load mode only chi0 is borrowed)
            so_members = _so_writer_members(
                "sc_so_static", _SO_STATIC_LABEL, so_run, _so_green_fp(green_kw),
                _so_own_state(eigenvalues, beta, mu, nmat, so_run.state.target))
            logger.info("Computing the spin-orbital pairing vertex...")
            Vs_q = _so_static_vertex(input_dict, chi0q_raw, norb_phys,
                                     Nx, Ny, Nz, static_index,
                                     so_context=so_run.context)
            del chi0q_raw
        logger.info("Spin-orbital vertex shape: {}".format(Vs_q.shape))

    elif chi0q_mode == "flex":
        # --- FLEX mode: use pre-computed dressed susceptibilities ---
        logger.info("FLEX mode: loading dressed susceptibilities and Green's function")

        chis, chic, green_dressed, chi_convention = _load_flex_susceptibilities(
            input_dict, norb, Nx, Ny, Nz, interactions=interactions)
        logger.info("FLEX susceptibility convention: {}".format(chi_convention))

        # Use dressed Green's function if available, otherwise use bare
        if green_dressed is not None:
            green_kw = green_dressed
            logger.info("Using FLEX dressed Green's function")
        else:
            logger.info("Calculating bare Green's function G(k, iwn)...")
            # Only ``full_sc`` is read on this path (no bond bubble here),
            # and ``full_sc`` is numerically independent of the tail
            # coefficient (build_green's full = deflated + tail add-back
            # recovers the plain bare Green function regardless of
            # ``coeff_tail``), so 0.0 is passed to skip the extra
            # tail-buffer construction rather than because it is the
            # physically meaningful choice here. At tail_off,
            # _build_bond_green transiently holds the canonical build_green
            # output and its sc-layout conversion (``full_sc``) at once
            # before releasing the canonical array -- unlike the bond
            # branch, no ``_bond_resource_preflight`` covers this transient
            # here (this path has no memory cap of its own).
            green_kw = _build_bond_green(
                eigenvalues, eigenvectors, mu, beta, nmat, 0.0, None).full_sc

        # Compute pairing vertex from FLEX susceptibilities
        logger.info("Computing FLEX vertices (pairing_type={})...".format(pairing_type))
        # Reject BEFORE the S/C build: the pair is O(Nq * norb^4) and an
        # unsupported calculation must not allocate heavily on its way to
        # the validation error (round-10 review).
        _reject_reduced_flex_unsupported(inter_k, chi_convention)
        # One S/C build for the run: the pair feeds the missing-component
        # diagnostic AND the vertex contraction (round-9 review; previously
        # each built its own full-grid pair).
        sc_mats = _build_vertex_sc_matrices(chi_convention, inter_k,
                                            norb, Nx, Ny, Nz)
        _warn_reduced_flex_missing_components(
            inter_k, norb, Nx, Ny, Nz, chi_convention,
            sc_matrices=(sc_mats if str(chi_convention).lower() == "kuroki"
                         else None))
        Vs_q = _compute_vertices_flex(chis, chic, inter_k, norb, Nx, Ny, Nz,
                                      pairing_type=pairing_type,
                                      convention=chi_convention,
                                      sc_matrices=sc_mats)
        logger.info("FLEX vertex shape: {}".format(Vs_q.shape))

    elif use_bond_channels:
        # --- Bond-resolved mode: the enlarged (orbital-pair)x(bond) ladder ---
        # This path builds its own bubble (bond_bubble_static) directly from the
        # Green function, so no chi0q file/tensor is read: the m=m'=0 block of
        # chi-bar IS the ordinary chi0q, and the m!=0 blocks (which no chi0q
        # file carries) are what the bond channels add.
        #
        # Resolve the bond topology and run the resource preflight (S3.2)
        # BEFORE building the Green carrier (review fix I-2: the carrier's
        # own bytes are folded into the preflight's estimate, and a
        # large-Nmat run must be refused here, not after the (possibly
        # large) Green-function allocation has already happened --
        # unified-bubble-kernel spec/test_sc_bond.py:382-399 pin this
        # ordering across the carrier switch too).
        bond_topo, bond_set, bond_split = _resolve_bond_topology(
            interactions, geom_info, norb, bond_max_shells)
        logger.info("Bond channels: B = %d, Delta r = %s, types = %s",
                    bond_set.n_channels, list(bond_set.delta_r),
                    list(bond_topo.coeffs))
        bond_channels.warn_pairing_approximations(
            bond_topo, "[eliashberg] bond_channels=true", log=logger)

        bond_coeff_tail_requested, bond_coeff_tail_source = (
            _resolve_bond_coeff_tail(eli_param))

        # The frequency grid the bond bubble will really be built on: for an
        # external bond_green that is the FILE's Nmat, which therefore has to
        # be known BEFORE the cap is applied -- otherwise an over-large
        # configured Nmat could refuse a run that fits comfortably (the file
        # is never reached). The header peek reads the npz member's shape
        # only, so the preflight still precedes every large allocation.
        nmat_bond = nmat
        if bond_green is not None:
            _refuse_so_green_file(bond_green)
            nmat_peek = _peek_green_npz_nfreq(bond_green, label="bond_green")
            if nmat_peek is not None:
                # Validate the EFFECTIVE (file) Nmat before the preflight: an
                # odd or non-positive grid is a silent-wrong-result, and it
                # must not even be used to size the cap check.
                nmat_bond = _validate_bond_green_nmat(nmat_peek, bond_green)
                if nmat_peek != nmat:
                    logger.warning(
                        "[eliashberg] bond_green %s carries %d Matsubara "
                        "frequencies while [mode.param] Nmat = %d; the FILE "
                        "wins (the Green function defines the frequency grid "
                        "the bond bubble is built on). The resource preflight "
                        "uses the file's Nmat.",
                        bond_green, nmat_peek, nmat)
        # tail_on = coeff_tail_applied != 0: known before the carrier is
        # built because an external bond_green ALWAYS resolves
        # coeff_tail_applied = 0.0 regardless of what was requested (see
        # _build_bond_green).
        tail_on = (bond_green is None) and (bond_coeff_tail_requested != 0.0)
        _bond_resource_preflight(norb, bond_set, Nx, Ny, Nz, nmat_bond,
                                 bond_memory_cap_gb, tail_on=tail_on)

        if bond_green is not None:
            logger.info(
                "Loading external Green's function from bond_green=%s...",
                bond_green)
        else:
            logger.info("Calculating Green's function G(k, iwn)...")
        # Externally supplied Green function (spec Goal): the milestone
        # green is FLEX-dressed and self-consistent, which the bare
        # transfer-Hamiltonian green is not. The bond path never reads a
        # chi file, so this is the only external input it takes. Built
        # AFTER the preflight above (review fix I-2 / the ordering pin).
        green_carrier = _build_bond_green(
            eigenvalues, eigenvectors, mu, beta, nmat,
            bond_coeff_tail_requested, bond_green,
            spin_orbital_mode=False)
        green_kw = green_carrier.full_sc
        if bond_green is not None:
            nmat_green = green_kw.shape[-1]
            if nmat_green != nmat_bond:
                # The header peek could not resolve the count (unusual npz
                # layout); the preflight above therefore ran with the
                # configured Nmat, so re-run it with the authoritative one.
                logger.warning(
                    "[eliashberg] bond_green %s carries %d Matsubara "
                    "frequencies while [mode.param] Nmat = %d; the FILE wins "
                    "(the Green function defines the frequency grid the bond "
                    "bubble is built on). Re-running the resource preflight "
                    "with the file's Nmat.", bond_green, nmat_green, nmat)
                _bond_resource_preflight(norb, bond_set, Nx, Ny, Nz,
                                         nmat_green, bond_memory_cap_gb,
                                         tail_on=tail_on)
        logger.info("g2_tail = %s (Matsubara tail correction for the pair "
                    "bubble); bond_coeff_tail = %s (source=%s, applied=%s); "
                    "Green frequency axis = %d points", g2_tail,
                    bond_coeff_tail_requested, bond_coeff_tail_source,
                    green_carrier.coeff_tail_applied, green_kw.shape[-1])
        if g2_tail:
            _warn_if_g2_tail_outside_asymptotic_regime(green_kw, beta)
        bond_operator, bond_provenance, bond_attribution = _build_bond_operator(
            bond_topo, bond_split, green_carrier, norb,
            kx_array, ky_array, kz_array, beta, pairing_type,
            bond_max_shells=bond_max_shells,
            bond_memory_cap_gb=bond_memory_cap_gb,
            precondition_opts=bond_precondition_opts,
            green_source=bond_green,
            g2_tail=g2_tail,
            bond_coeff_tail_source=bond_coeff_tail_source)
        Vs_q = None

    else:
        # --- Standard RPA mode ---

        # Step 6: Calculate bare Green's function
        logger.info("Calculating Green's function G(k, iwn)...")
        # Same rationale as the FLEX-fallback branch above: only full_sc is
        # consumed here, and it is independent of coeff_tail. Also same
        # transient-memory caveat: _build_bond_green briefly holds the
        # canonical build_green output alongside its sc-layout conversion
        # before releasing the former, and no resource preflight covers
        # this Standard-RPA path.
        green_kw = _build_bond_green(
            eigenvalues, eigenvectors, mu, beta, nmat, 0.0, None).full_sc

        # Step 1: Load or compute chi0q
        if chi0q_mode == "calc":
            # the static solve reads only the zero bosonic frequency, so only
            # that slice is computed (one-frequency axis, issue #272)
            chi0q_raw = _calc_chi0q_internal(input_dict, chi0q_tensor=chi0q_tensor,
                                                precomputed_mu=mu,
                                                prepared=band,
                                                static_only=True)
            band = None
            static_index = 0
        else:
            chi0q_raw, static_index = _load_chi0q(input_dict, norb=norb)

        # Step 8: Convert chi0q format; the frequency axis is last in the
        # reference format, so read its length after the conversion (a 6D
        # input can be either raw H-wave (nmat first) or ref (nmat last))
        chi0q = _convert_chi0q_to_ref_format(chi0q_raw, norb, Nx, Ny, Nz)
        nmat_chi0q = chi0q.shape[-1]
        logger.info("chi0q converted to shape: {}".format(chi0q.shape))

        # Step 9: Compute RPA vertices
        logger.info("Computing RPA vertices (pairing_type={})...".format(pairing_type))
        vertex_result = _compute_vertices(
            chi0q, inter_k, norb, Nx, Ny, Nz, nmat_chi0q,
            pairing_type=pairing_type, static_index=static_index,
            declarations_closed=_declarations_partner_closed(
                interactions, Nx, Ny, Nz, norb))
        if isinstance(vertex_result, tuple):
            Pc_q, Ps_q = vertex_result
            Vs_q = Pc_q + Ps_q
            logger.info("Simple mode: Pc + Ps vertex, shape {}".format(Vs_q.shape))
        else:
            Vs_q = vertex_result
            logger.info("General mode: 4-index V^s vertex, shape {}".format(Vs_q.shape))

    # --- Step 10: Compute G2 ---
    # The bond operator carries its own pair bubble (the same construction as
    # _calc_g2, ``bond_channels._g2_from_green`` inside
    # _build_bond_operator), so G2 is not
    # needed -- and must not be built -- on that path.
    if use_bond_channels:
        G2 = None
    else:
        logger.info("Computing G2...")
        # New config keys are read case-insensitively (the config layers
        # disagree on case handling; see the PR #128 sweep) and through the
        # strict parser: a misspelled value must fail, not silently flip the
        # physics.
        logger.info("g2_tail = %s (Matsubara tail correction for the pair "
                    "bubble); Green frequency axis = %d points", g2_tail,
                    green_kw.shape[-1])
        if g2_tail:
            _warn_if_g2_tail_outside_asymptotic_regime(green_kw, beta)
        G2 = _calc_g2(green_kw, beta, tail=g2_tail)
        _warn_if_g2_indefinite(G2, norb, g2_tail)

    # --- Step 11: Initialize gap function ---
    if spin_orbital_mode:
        sigma_init = so_pairing.so_initial_gap(
            init_gap_mode, norb_phys, kx_array, ky_array, kz_array)
        sigma_init = _project_gap_parity(sigma_init, "triplet")
        if np.linalg.norm(sigma_init) == 0.0:
            raise ValueError(
                "[eliashberg] init_gap = {!r} has no fermion-antisymmetric "
                "component on this k grid; choose another init_gap.".format(
                    init_gap_mode))
        sigma_init = sigma_init / np.linalg.norm(sigma_init)
    else:
        sigma_init = _initialize_gap(init_gap_mode, norb, kx_array, ky_array, kz_array)

    # S_z-conservation residual of the spin-orbital vertex and Green function
    # (issue #88), once per run, on the host before any device transfer
    so_sz_residual = None
    if spin_orbital_mode:
        so_sz_residual = so_pairing.sz_conservation_residual(
            V=Vs_q, green=green_kw, norb_phys=norb_phys)

    # GPU: park the two large invariants (pairing vertex and pair bubble) on the
    # device once; each matvec then only moves the gap vector across PCIe. The
    # solver entry points below pass Vs_q/G2 straight to _make_kernel_operator,
    # which derives its backend from them, so no other change is needed here.
    if gpu_active and not use_bond_channels:
        # Resident device tensors are more than the two inputs: the operator
        # also holds V_r (~Vs_q) and G2_pre (~G2), and in general mode
        # G2_gemm (~G2) + V_r_gemm (~Vs_q) -- i.e. up to ~3x the inputs. The
        # per-matmat block workspace (gap-sized x num columns) adds on top; a
        # subspace/eigenvalue run with a wide block can still exceed this. It is
        # only a pre-warning -- CuPy raises a clear OutOfMemoryError regardless.
        est_bytes = 3 * (Vs_q.nbytes + G2.nbytes)
        backend.warn_if_device_memory_short(
            est_bytes, logger, label="the static Eliashberg kernel")
        logger.info("GPU backend active (CuPy): moving the pairing vertex and "
                    "G2 to the device (%.2f GB).",
                    (Vs_q.nbytes + G2.nbytes) / 1e9)
        Vs_q = xp.asarray(Vs_q)
        G2 = xp.asarray(G2)

    # --- Step 12: Solve ---
    sigma_result = None
    eigenvalue_iter = None
    eigenvalues_eig = None
    eigenvectors_eig = None
    eigenvalue_match = None
    # Live only on the iteration path; kept in scope so the output artifact can
    # record WHICH quantity the written number is and whether it converged
    # (review fix C1).
    converged = None
    n_iter = None

    # Spectral shift of the POWER ITERATION (None unless requested). When set,
    # the iteration runs on K + sigma*I and its reported eigenvalue is the
    # signed eigenvalue of K, which changes how the result must be labelled
    # below; keep it in scope for that.
    # An absent spectral_shift means "auto" for the power iteration
    # (_resolve_spectral_shift); "none" keeps it unshifted.
    iteration_spectral_shift = (
        _resolve_spectral_shift(eli_param.get("spectral_shift"), "iteration")
        if solver_mode in ("iteration", "both") else None)
    eigen_spectral_shift = _resolve_spectral_shift(
        eli_param.get("spectral_shift"), "eigenvalue")

    # Rayleigh-validation record of a SHIFTED power iteration (empty on the
    # unshifted default path); decides whether the number written to
    # eigenvalue.dat may be called an eigenvalue at all.
    iteration_info = {}

    if spin_orbital_mode:
        # both solvers below work on the antisymmetric coordinates of the
        # gap (E^dagger K E): the Pauli-forbidden sector is not in the space
        so_basis = _SOAntisymmetricBasis((norb, norb, Nx, Ny, Nz))
        A_so, _ = _make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)
        kernel_operator = None
        # the gap is one spin-orbital matrix: the iteration's channel parity
        # projection/leakage probe does not apply (P is in the operator), and
        # the eigenpairs are tagged by fermion antisymmetry
        pairing_type = None
        parity_channel = "triplet"
    else:
        kernel_operator = bond_operator
        parity_channel = pairing_type

    if solver_mode in ("iteration", "both"):
        logger.info("=== Self-consistent iteration ===")
        if spin_orbital_mode:
            sigma_result, eigenvalue_iter, converged, n_iter = (
                _so_solve_iteration(
                    A_so, so_basis, sigma_init, max_iter, alpha, tol,
                    spectral_shift=iteration_spectral_shift,
                    info_out=iteration_info))
        else:
            sigma_result, eigenvalue_iter, converged, n_iter = _solve_iteration(
                green_kw, Vs_q, G2, sigma_init, norb,
                max_iter=max_iter, alpha=alpha, tol=tol,
                pairing_type=pairing_type,
                operator=kernel_operator,
                spectral_shift=iteration_spectral_shift,
                info_out=iteration_info,
            )
        logger.info("Iteration result: eigenvalue = {:.6f}, converged = {}, n_iter = {}".format(
            eigenvalue_iter, converged, n_iter))
        if not converged:
            logger.warning(
                "[eliashberg] the power iteration did not converge (%s): the "
                "iteration eigenvalue is written as NaN with eigen_converged "
                "= false, and the gap file holds the last iterate, NOT an "
                "eigenvector. Use solver_mode = \"eigenvalue\" or "
                "spectral_shift = \"auto\".",
                "sign-alternating two-cycle: negative dominant eigenvalue"
                if iteration_info.get("two_cycle")
                else "after {} steps".format(n_iter))
        if spin_orbital_mode and not _so_is_antisymmetric(sigma_result):
            raise ValueError(
                "[mode] enable_spin_orbital = true: the power iteration "
                "returned a gap that is not fermion-antisymmetric "
                "(Delta_AB(k) = -Delta_BA(-k)) to round-off; it is not "
                "written.")

    eigen_leg_ok = True
    eigen_leading = None
    channel_lead = None
    channel_lead_gap = None
    analysis_note = None
    sector_eigenvalues = None
    channel_found_flag = False
    if solver_mode in ("eigenvalue", "both"):
        logger.info("=== Eigenvalue analysis ===")
        if spin_orbital_mode:
            eigenvalues_eig, eigenvectors_eig = _so_solve_eigenvalue(
                A_so, so_basis, num_eigenvalues, eigenvalue_method,
                sigma_shift=eli_param.get("sigma_shift"),
                spectral_shift=eigen_spectral_shift)
        else:
            eigenvalues_eig, eigenvectors_eig = _solve_eigenvalue(
                Vs_q, G2, norb, Nx, Ny, Nz,
                num_eigenvalues=num_eigenvalues,
                method=eigenvalue_method,
                sigma_shift=eli_param.get("sigma_shift"),
                spectral_shift=eigen_spectral_shift,
                operator=kernel_operator,
            )
        # The kernel preserves parity; promote the eigenpairs whose gap has the
        # requested channel parity (singlet even / triplet odd). On the
        # non-spin-orbital path with a parity-commuting kernel the reported
        # leader comes from a solve on the channel's parity sector and is
        # prepended to the table as row 0 (_channel_eigenpairs,
        # _listing_with_leader); a non-commuting kernel reports the retained
        # eigenpair of (dominantly) the requested parity, or NaN.
        if spin_orbital_mode:
            eigenvalues_eig, eigenvectors_eig = _reorder_eigenpairs_by_parity(
                eigenvalues_eig, eigenvectors_eig, parity_channel)
        else:
            op_static = (kernel_operator if kernel_operator is not None
                         else _make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz))
            (eigenvalues_eig, eigenvectors_eig, channel_lead,
             channel_lead_gap, channel_found, analysis_note,
             sector_eigenvalues) = \
                _channel_eigenpairs(eigenvalues_eig, eigenvectors_eig,
                                    op_static, (norb, norb, Nx, Ny, Nz),
                                    parity_channel, num_eigenvalues,
                                    spectral_shift=eigen_spectral_shift)
            channel_found_flag = bool(channel_found
                                      and analysis_note is not None)
            if not channel_found:
                eigen_leg_ok = False
            elif analysis_note is not None:
                # a failed global solve (its NaN sentinel, placeholder
                # vectors) lists nothing; the sector solve's leader stands
                listed = np.isfinite(np.asarray(eigenvalues_eig))
                if not np.all(listed):
                    logger.warning(
                        "The global eigen-analysis set could not be "
                        "computed; the table lists only the '%s' leader of "
                        "the parity-sector solve.", parity_channel)
                    eigenvalues_eig = np.asarray(eigenvalues_eig)[listed]
                    eigenvectors_eig = np.asarray(eigenvectors_eig)[listed]
                eigenvalues_eig, eigenvectors_eig = _listing_with_leader(
                    eigenvalues_eig, eigenvectors_eig, channel_lead,
                    channel_lead_gap)
        # Flag whether each eigenvector lies (dominantly) in the channel's
        # parity sector (even for singlet, odd for triplet). _is_gap_parity is a
        # threshold test, so this is "dominantly the channel parity" rather than
        # an exact label; for a centrosymmetric model parity is an exact kernel
        # symmetry and the opposite-parity modes are unphysical for the channel
        # (e.g. an even mode in the triplet kernel is Pauli-forbidden).
        if spin_orbital_mode:
            # round-off criterion, not the dominant-sector test
            eigenvalue_match = np.array(
                [_so_is_antisymmetric(g) for g in eigenvectors_eig])
            _so_check_antisymmetric_eigenpairs(eigenvalues_eig,
                                               eigenvalue_match)
        else:
            eigenvalue_match = np.array(
                [_is_gap_parity(g, parity_channel) for g in eigenvectors_eig])
        logger.info("Leading eigenvalues:")
        for i, ev in enumerate(eigenvalues_eig):
            tag = "" if eigenvalue_match[i] else "  [opposite-parity sector]"
            logger.info("  {:3d}: {:.6f} (|ev| = {:.6f}){}".format(
                i, ev.real, abs(ev), tag))

        # the reported leader: the channel solve on the non-spin-orbital
        # path (_channel_eigenpairs), the leading listed eigenpair otherwise
        if channel_lead is None:
            channel_lead = eigenvalues_eig[0] if len(eigenvalues_eig) \
                else complex(np.nan)
            channel_lead_gap = eigenvectors_eig[0] if len(eigenvectors_eig) \
                else None
        eigen_leg_ok = bool(eigen_leg_ok and np.isfinite(channel_lead))
        eigen_leading = (float(np.real(channel_lead)) if eigen_leg_ok
                         else float("nan"))
        lead_gap = (channel_lead_gap if channel_lead_gap is not None
                    else eigenvectors_eig[0])
        # Use the leading eigenvector as gap if no iteration result; in
        # "both" the reported value is the eigenvalue analysis', so the gap
        # written with it is its eigenvector too (when it was obtained)
        if sigma_result is None:
            sigma_result = lead_gap
            eigenvalue_iter = eigen_leading
        elif eigen_leg_ok:
            sigma_result = lead_gap

    # --- Step 12b: bare-pp vs fluctuation attribution of lambda (S4.5) ---
    # Real Rayleigh quotients of the two kernel parts on the symmetrized
    # kernel, evaluated at the converged gap. The preconditions that make them
    # real were already enforced when the operator was built.
    if bond_attribution is not None and sigma_result is not None:
        attr = bond_channels.attribute_lambda(
            sigma_result, bond_attribution["weight"],
            bond_attribution["pp"], bond_attribution["fl"],
            op_full=bond_attribution["full"])
        logger.info(
            "lambda attribution (spec S4.5): lambda = %.8f = lambda^pp %.8f "
            "+ lambda^fl %.8f (sum residual %.2e)",
            attr["lambda"], attr["lambda_pp"], attr["lambda_fl"],
            attr["sum_residual"])
        if not attr["imag_within_tol"]:
            logger.warning(
                "lambda attribution: the Rayleigh quotients carry a non-"
                "negligible imaginary part (Im lambda = %.3e, pp %.3e, fl "
                "%.3e). The reported real parts are still the projections "
                "onto the Hermitian path, but check the gap's convergence.",
                attr["imag"], attr["imag_pp"], attr["imag_fl"])
        bond_provenance["lambda_rayleigh"] = "{:.8e}".format(attr["lambda"])
        bond_provenance["lambda_pp"] = "{:.8e}".format(attr["lambda_pp"])
        bond_provenance["lambda_fl"] = "{:.8e}".format(attr["lambda_fl"])
        bond_provenance["lambda_attribution"] = (
            "lambda = lambda_pp + lambda_fl, real Rayleigh quotients of the "
            "instantaneous (bare particle-particle) and fluctuation parts on "
            "the symmetrized kernel K~ = -sqrt(w) Gamma sqrt(w) (spec S4.5); "
            "evaluated at the RETURNED gap, so lambda_rayleigh equals the "
            "solver's eigenvalue only when that gap is a converged "
            "eigenvector (and, unlike the iteration-mode norm, it carries "
            "the sign of the eigenvalue)")
        # Which solver produced the gap the Rayleigh quotient was evaluated at,
        # and -- on the iteration path -- whether that gap converged. Without
        # this the file carries two numbers (lambda_rayleigh and the iteration
        # value) that a reader cannot tell apart (review fix C1a).
        # In "both" the written gap is the eigenvalue analysis' eigenvector
        # when that leg succeeded, so the attribution refers to it; the
        # iteration's semantics apply only when the iteration gap is written.
        gap_from_eigen = (solver_mode == "both" and eigen_leg_ok)
        bond_provenance["lambda_rayleigh_solver_mode"] = solver_mode
        if solver_mode == "both":
            bond_provenance["lambda_rayleigh_gap_source"] = (
                "eigenvalue analysis" if gap_from_eigen else "iteration")
        bond_provenance["lambda_rayleigh_converged"] = (
            str(bool(eigen_leg_ok)) if gap_from_eigen
            else str(bool(converged)) if solver_mode in ("iteration", "both")
            else "n/a")

        # The reported value (the iteration's SIGNED Rayleigh eigenvalue, or
        # the eigenvalue analysis' leader in "both") must agree with the
        # symmetrized-kernel Rayleigh quotient at the same gap.
        reported = eigen_leading if gap_from_eigen else eigenvalue_iter
        if (solver_mode in ("iteration", "both")
                and reported is not None
                and np.isfinite(reported)):
            lam_signed = attr["lambda"]
            if abs(lam_signed - reported) > 1.0e-6 * max(
                    1.0, abs(lam_signed)):
                logger.warning(
                    "[eliashberg] solver_mode='%s' reports %.8e but the "
                    "symmetrized-kernel Rayleigh quotient at the same gap is "
                    "lambda_rayleigh = %.8e; check the convergence "
                    "(convergence_tol) or use solver_mode=\"eigenvalue\".",
                    solver_mode, reported, lam_signed)

    # --- Step 12c: opt-in bond diagnostics (S7.7 character analysis) ---
    # Additive comment lines only; off by default, so a run without the flag is
    # byte-identical to what it was before the diagnostics existed.
    if (bond_diagnostics and bond_provenance is not None
            and bond_attribution is not None):
        _bond_diagnostics_record(
            bond_provenance, sigma_result, eigenvalues_eig,
            bond_attribution["weight"], norb,
            kx_array, ky_array, kz_array,
            sector_eigenvalues=(sector_eigenvalues if channel_found_flag
                                else None),
            leader_available=(eigen_leg_ok if eigenvalues_eig is not None
                              else True))

    # --- Step 13: Save results ---
    # Label the iteration number in the file itself, so it can never be
    # silently confused with the signed lambda_rayleigh sitting next to it
    # (review fix C1c).
    #
    # The note is written when it carries information:
    #   * a spectral_shift (the iteration default "auto") -- the value is the
    #     SIGNED eigenvalue of K, validated by the Rayleigh residual;
    #   * a non-converged iteration -- the value is NaN, and the note says
    #     why (cap, sign-alternating two-cycle, failed residual check); or
    #   * bond_channels -- eigenvalue.dat then also carries lambda_rayleigh,
    #     and the two numbers must be told apart.
    # The "# eigen_converged = <true|false>" line is written for every run
    # (_save_results).
    eigenvalue_note = None
    eigen_converged = None
    # solver_mode = "both": eigen_converged, the returned value and the exit
    # status follow the eigenvalue analysis; the iteration's own convergence
    # is recorded as iteration_converged
    iteration_converged = (bool(converged) if solver_mode == "both"
                           and eigenvalue_iter is not None else None)
    if solver_mode in ("iteration", "both") and eigenvalue_iter is not None:
        eigen_converged = bool(converged)
        if iteration_spectral_shift is not None or not eigen_converged:
            # what the number is: the SIGNED (validated) eigenvalue of K, or
            # NaN with the reason the iteration did not converge
            eigenvalue_note = _shifted_eigenvalue_note(
                solver_mode, iteration_spectral_shift, eigen_converged,
                n_iter, iteration_info)
        elif use_bond_channels:
            eigenvalue_note = (
                "solver_mode='{}': the value below is the SIGNED Rayleigh "
                "eigenvalue of the converged power iterate ({} steps)."
                .format(solver_mode, n_iter))
        if use_bond_channels and bond_provenance \
                and "lambda_rayleigh" in bond_provenance:
            eigenvalue_note += (
                " lambda_rayleigh = {} above is the symmetrized-kernel "
                "Rayleigh quotient at the same gap.".format(
                    bond_provenance["lambda_rayleigh"]))
    if solver_mode in ("eigenvalue", "both") and eigenvalues_eig is not None:
        # NaN when the largest-real-part solve failed on an operator too large
        # for the dense fallback (_solve_leading), or when no eigenpair of the
        # requested channel could be obtained (_channel_eigenpairs)
        eigen_converged = bool(eigen_leg_ok)

    so_weights = None
    so_gap_comments = None
    if spin_orbital_mode:
        so_out = _so_output_record(
            output_dir, gap_file, norb_phys, sigma_result, eigenvalues_eig,
            eigenvectors_eig,
            eigenvalue_match, kx_array, ky_array, kz_array,
            sz_residual=so_sz_residual, provenance_members=so_members,
            eigen_converged=eigen_converged)
        bond_provenance, so_weights, so_gap_comments = so_out

    _save_results(
        output_dir, sigma_result, eigenvalue_iter, eigenvalues_eig,
        kx_array, ky_array, kz_array,
        gap_file=gap_file, eigenvalue_file=eigenvalue_file,
        eigenvalue_match=eigenvalue_match,
        provenance=bond_provenance,
        eigenvalue_note=eigenvalue_note,
        eigenvalue_weights=so_weights,
        gap_comments=so_gap_comments,
        eigen_converged=eigen_converged,
        iteration_converged=iteration_converged,
        analysis_note=analysis_note,
        leading_eigenvalue=(
            (complex(channel_lead) if eigen_leg_ok
             else complex(np.nan, np.nan))
            if eigenvalues_eig is not None and channel_lead is not None
            else None),
    )

    logger.info("=== Done ===")
    # the reported leading eigenvalue -- the eigenvalue analysis' for
    # solver_mode "eigenvalue" and "both", the iteration's for "iteration";
    # NaN when it could not be determined. main() turns NaN into exit 3.
    if eigen_leading is not None:
        return eigen_leading
    if eigenvalue_iter is None:
        return None
    return float(np.real(eigenvalue_iter))


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    """Command-line interface for the Eliashberg equation solver."""
    import tomli
    import argparse

    parser = argparse.ArgumentParser(
        description="Linearized Eliashberg equation solver (post-tool for H-wave RPA)",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("input", type=str, help="input TOML file (same as RPA input + [eliashberg] section)")
    parser.add_argument("-q", "--quiet", action="store_true", help="suppress output")
    parser.add_argument(
        "-v", "--version", action="version",
        version="hwave_sc v{}".format(hwave.__version__),
    )
    args = parser.parse_args()

    # Setup logging
    log_level = logging.WARNING if args.quiet else logging.INFO
    logging.basicConfig(level=log_level, format="%(name)s: %(message)s")

    file_toml = args.input
    if not os.path.exists(file_toml):
        logger.error("Input file does not exist: {}".format(file_toml))
        sys.exit(1)

    with open(file_toml, "rb") as f:
        input_dict = tomli.load(f)

    lam = calc_eliashberg(input_dict)
    if isinstance(lam, float) and np.isnan(lam):
        # the outputs are written (lambda = NaN, eigen_converged = false);
        # a script or job chain must still see that no eigenvalue exists
        logger.error(
            "The Eliashberg eigen-solver did not converge: no leading "
            "eigenvalue was obtained (written as NaN with eigen_converged = "
            "false). Use [eliashberg] solver_mode = \"eigenvalue\" or "
            "spectral_shift = \"auto\". Exit status %d.",
            EXIT_EIGEN_NOT_CONVERGED)
        sys.exit(EXIT_EIGEN_NOT_CONVERGED)


if __name__ == "__main__":
    main()
