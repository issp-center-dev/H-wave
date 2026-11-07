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
from typing import Any, NamedTuple

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

from hwave.solver.eliashberg.static import (  # noqa: F401,E402
    reject_spin_orbital_mode,
    EliashbergRun,
    _read_mode_params,
    _read_static_options,
    _read_model,
    _check_bond_and_archives,
    _prepare_band,
    _build_one_body,
    _vertex_spin_orbital,
    _vertex_flex,
    _vertex_bond,
    _vertex_rpa,
    _build_pairing_vertex,
    _build_g2_and_seed,
    _prepare_solvers,
    _iteration_leg,
    _eigenvalue_leg,
    _attribute_lambda,
    _write_results,
    calc_eliashberg,
    main,
    EXIT_EIGEN_NOT_CONVERGED,
    _SO_STATIC_LABEL,
    DEFAULT_SOLVER_MODE,
    _SO_DYNAMIC_LABEL,
    IN_PROCESS_DEFAULT_SOLVER_MODE,
)

logger = logging.getLogger("hwave_sc")


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


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    main()
