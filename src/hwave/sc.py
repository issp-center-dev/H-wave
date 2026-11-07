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

from hwave.solver.eliashberg.io import (  # noqa: F401,E402
    _static_freq_position,
    _read_freq_meta,
    _require_so_chi0q_marker,
    _refuse_so_archive_in_non_so_mode,
    _refuse_so_green_file,
    _preflight_refuse_so_archives,
    _load_chi0q,
    _auto_chi0q_tensor_scheme,
    _calc_chi0q_internal,
    _chi0q_tensor_scheme,
    _chi0q_from_prepared_solver,
    _read_spin_free_tol,
    _load_band_into_solver,
    _mean_field_band,
    _read_interaction_files,
    _read_so_extern,
    _so_check_archive,
    _so_green_fp,
    _resolve_flex_dir,
    _resolve_flex_paths,
    _load_flex_susceptibilities_full,
    _accept_up_block_only,
    _legacy_up_block_only,
    _onsite_transpose_asymmetry,
    _onsite_hermitian_pair_mismatch,
    _read_flex_chi_raw,
    _expand_flex_chi,
    _load_green_npz,
    _UnsupportedNpyHeaderVersion,
    _read_npy_header_shape,
    _peek_green_npz_nfreq,
    _load_flex_green,
    _load_flex_susceptibilities,
    _convert_chi0q_to_ref_format,
    _so_flex_static_slice,
    _so_flex_full_axis,
    _SOFlexBorrowed,
    _so_flex_sigma_path,
    _load_flex_so_sigma,
    _load_flex_so,
    _SO_DYNAMIC_IR_HINT,
    _STATIC_IR_HINT,
    _SO_FLEX_CHI_CONVENTION,
    _SO_FLEX_PRODUCER,
    _DEFAULT_NMAT,
    _SO_RPA_PRODUCER,
)
from hwave.solver.eliashberg.output import (  # noqa: F401,E402
    _shifted_eigenvalue_note,
    _save_results,
)

from hwave.solver.eliashberg.bond import (  # noqa: F401,E402
    _bond_memory_estimate,
    _bond_resource_preflight,
    _bond_locality_split,
    _resolve_bond_topology,
    _build_bond_vertices,
    _bond_kernel_parts,
    _build_bond_operator,
    _build_bond_green,
    _bond_diagnostics_record,
    _BOND_DIAGNOSTICS_DEG_TOL,
    _BOND_N_Q_ARRAYS,
)
from hwave.solver.eliashberg.so_solve import (  # noqa: F401,E402
    _so_coeff_extern,
    _SORun,
    _so_run_inputs,
    _so_own_state,
    _so_flex_consumed_state,
    _so_writer_members,
    _so_epsilon_k,
    _so_rpa_vertex,
    _so_gap_npz_name,
    _so_flex_writer_members,
    _so_prov_line,
    _so_provenance_lines,
    _so_output_record,
    _so_static_vertex_from_chi,
    _so_static_chi0,
    _so_static_vertex,
    _so_flex_backward_error,
    _check_so_flex_vertex,
    SectorOp,
    _so_antisymmetry_residual,
    _so_is_antisymmetric,
    _so_solve_iteration,
    _so_dense_reduced_eigen,
    _so_solve_eigenvalue,
    _pp_label,
    SectorLeader,
    SectorSolve,
    IterationLeader,
    SectorIteration,
    EigenLeg,
    IterationLeg,
    SectorReport,
    _so_canonical_sector_ops,
    _so_check_sector_op,
    _so_lifted_unit_gaps,
    _so_sector_norm_estimate,
    _so_union,
    _so_order_by_seed,
    _so_solve_sectors,
    _so_iterate_sectors,
    _so_bond_report,
    _so_check_antisymmetric_eigenpairs,
    _SO_PP_DENSE,
    _SO_PP_ORDER,
    _SO_SECTOR_RESIDUAL_TOL,
    _SO_SEED_PART_RTOL,
    _SO_SEED_OVERLAP_WARN,
    _SO_SECTOR_TIE_TOL,
    _SO_NORM_STEPS,
    _SO_NORM_RTOL,
    _SO_ANTISYM_TOL,
    _SO_FLEX_GAMMA_TOL,
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


# ---------------------------------------------------------------------------
# Eliashberg frequency mode dispatch
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# k-space construction
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Green's function
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# RPA vertices
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# G2 and Eliashberg kernel
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Gap function initialization
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# chi0q format conversion
# ---------------------------------------------------------------------------


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


#: producer labels of the hwave_sc spin-orbital gap archives (#88)
_SO_STATIC_LABEL = "hwave_sc (spin-orbital, static)"
_SO_DYNAMIC_LABEL = "hwave_sc (spin-orbital, dynamic)"


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
