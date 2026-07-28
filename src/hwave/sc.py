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

import numpy as np
from numpy.fft import fftn, ifftn
from scipy.optimize import bisect
from scipy.sparse.linalg import LinearOperator, eigs, bicgstab, gmres, lgmres

import hwave
import hwave.qlmsio.wan90 as wan90
from hwave.solver.rpa import validate_chi0q_index_convention
from hwave.solver.ir_axis import is_ir_native, ir_native_meta
from hwave.solver import backend
from hwave.solver import bond_channels

logger = logging.getLogger("hwave_sc")

# Number of independent random probes used by _solve_iteration to test whether
# the Eliashberg kernel commutes with parity before enabling gap-parity
# projection. More than one guards against a single probe accidentally lying
# near the nullspace of the commutator and missing a non-centrosymmetric kernel.
_PARITY_GUARD_PROBES = 3

# spectral_shift="auto" on the POWER-ITERATION path (see
# _resolve_iteration_shift): number of probe matvecs per probe vector used to
# estimate the spectral radius rho(K), and the safety factor applied to it.
# sigma = _AUTO_SHIFT_MARGIN * rho_est + 1e-6 sits just above the estimated
# radius, which is the smallest shift that makes every eigenvalue of K + sigma*I
# have a positive real part (so the power iteration converges to the
# algebraically largest eigenvalue of K). The margin absorbs the fact that a
# power-iteration estimate approaches rho FROM BELOW; it is deliberately small
# because a larger shift compresses the relative eigenvalue gaps and therefore
# slows the iteration down.
_AUTO_SHIFT_PROBE_STEPS = 15
_AUTO_SHIFT_MARGIN = 1.05

# Default Matsubara grid size (mode.param.Nmat). Single definition so the
# legacy-file fallback in _static_freq_position and the main solver loop can
# never silently disagree on the grid.
_DEFAULT_NMAT = 1024

# Default peak-memory ceiling (GB) of the bond-channel preflight
# ([eliashberg] bond_memory_cap_gb; spec S3.2/S5).
_BOND_MEMORY_CAP_GB = 8.0

# Tolerance on |Im V_ab(R)| above which a CoulombInter entry counts as complex
# and the v1 bond path refuses to run (spec S5 guard).
_BOND_IMAG_TOL = 1.0e-12

# Every [eliashberg] key that only means something when bond_channels=true;
# set while the flag is off they are ignored with a warning (spec S5).
_BOND_OPTION_KEYS = (
    "bond_diagnostics",
    "bond_green",
    "bond_max_shells",
    "bond_memory_cap_gb",
    "bond_precondition_atol",
    "bond_precondition_rtol",
    "bond_precondition_dense_limit",
)

# Degeneracy tolerance of the opt-in bond diagnostics' eigenvalue clustering
# ([eliashberg] bond_diagnostics; spec S7.7 default).
_BOND_DIAGNOSTICS_DEG_TOL = 1.0e-3

# Boolean spellings accepted for the bond-channel [eliashberg] keys -- exactly
# the ones backend.as_bool documents. ANY other string/type is REFUSED (see
# _bond_bool_option): as_bool's plain-truthiness fallback would read a typo
# like "ture" as false, silently disabling the requested feature AND the
# safety guards that only run with it.
_BOND_BOOL_TRUE = ("true", "yes", "on", "1")
_BOND_BOOL_FALSE = ("false", "no", "off", "0")


def _bond_bool_option(eli_param, key, default=False):
    """Read a bond-channel boolean option, refusing unrecognised spellings.

    A key present with the value ``None`` counts as unset (TOML has no null,
    so that can only come from a programmatic caller meaning "not
    specified").
    """
    val = eli_param.get(key)
    if val is None:
        return default
    if isinstance(val, bool):
        return val
    if isinstance(val, str):
        spelling = val.strip().lower()
        if spelling in _BOND_BOOL_TRUE:
            return True
        if spelling in _BOND_BOOL_FALSE:
            return False
    raise ValueError(
        "[eliashberg] {} must be a boolean: true/false, or one of the strings "
        "{} / {}; got {!r}. An unrecognised value is REFUSED rather than read "
        "as false -- a typo would otherwise silently disable both the "
        "requested feature and the guards that come with it.".format(
            key, "/".join(_BOND_BOOL_TRUE), "/".join(_BOND_BOOL_FALSE), val))


def _bond_float_option(eli_param, key, default=None):
    """Read a bond-channel numeric option, naming the key on a bad value.

    A bare ``float(...)`` raises a message that never says WHICH
    ``[eliashberg]`` key was wrong; booleans are rejected outright (``True``
    would silently become ``1.0``). ``None`` counts as unset.
    """
    val = eli_param.get(key)
    if val is None:
        return default
    if isinstance(val, bool):
        raise ValueError(
            "[eliashberg] {} must be a number, not a boolean; got {!r}."
            .format(key, val))
    try:
        return float(val)
    except (TypeError, ValueError):
        raise ValueError(
            "[eliashberg] {} must be a number; got {!r}.".format(key, val)
        ) from None


def _bond_int_option(eli_param, key, default=None):
    """Read a bond-channel integer option WITHOUT silently truncating it.

    ``int(1.5)`` used to turn a malformed value into a different, valid
    configuration (``bond_max_shells = 1``), and TOML booleans were accepted
    as 0/1. Both change the MEANING of a bad config, so they are refused here,
    at the TOML boundary, with the key named. ``None`` counts as unset.
    """
    val = eli_param.get(key)
    if val is None:
        return default
    if isinstance(val, bool):
        raise ValueError(
            "[eliashberg] {} must be an integer, not a boolean; got {!r}."
            .format(key, val))
    fval = _bond_float_option(eli_param, key)
    if not np.isfinite(fval) or fval != np.floor(fval):
        raise ValueError(
            "[eliashberg] {} must be an integer; got {!r} (a non-integral or "
            "non-finite value is refused rather than truncated).".format(
                key, val))
    return int(fval)


# ---------------------------------------------------------------------------
# Eliashberg frequency mode dispatch
# ---------------------------------------------------------------------------

def _eliashberg_frequency(input_dict):
    """Validate and return the eliashberg frequency mode.

    Parameters
    ----------
    input_dict : dict
        Parsed TOML configuration dictionary.

    Returns
    -------
    str
        Either "static" (default) or "dynamic".

    Raises
    ------
    ValueError
        If the frequency mode is not in the allowed set.
    """
    freq = input_dict.get("eliashberg", {}).get("frequency", "static")
    if freq not in ("static", "dynamic"):
        raise ValueError(
            "eliashberg.frequency must be 'static' or 'dynamic', got '{}'"
            .format(freq))
    return freq


def _warn_dynamic_bond_ignores_chi(eli_param):
    """Warn that ``chi0q_mode`` is inert in the DYNAMIC bond branch.

    State-machine row 3 of the dynamic-bond spec (S5.2): with
    ``bond_channels=true`` the frequency-resolved chi-bar is built from the
    Green function inside the bond machinery, so no chi file is read or
    computed -- ``chi0q_mode`` (whatever it says, including the ``'flex'``
    the scalar dynamic path requires) is meaningless there. Warn-and-ignore
    rather than silently accepting a setting that has no effect.
    """
    mode = eli_param.get("chi0q_mode")
    if mode is not None:
        logger.warning(
            "[eliashberg] chi0q_mode='%s' is IGNORED with bond_channels=true "
            "and frequency='dynamic': the bond kernel builds its own "
            "frequency-resolved chi-bar from the Green function "
            "([eliashberg] bond_green, or the bare H0 green), so NO chi file "
            "(chiq_s/chiq_c/chi0q) is read or computed on this path.", mode)


def _validate_dynamic_prereqs(input_dict):
    """Validate prerequisites for dynamic Eliashberg calculation.

    Parameters
    ----------
    input_dict : dict
        Parsed TOML configuration dictionary.

    Raises
    ------
    ValueError
        If chi0q_mode is not "flex" on the scalar dynamic path, or if Nmat is
        odd.

    Notes
    -----
    ``bond_channels=true`` is now SUPPORTED on the dynamic path (dynamic
    bond-channel spec S3.4/S5.2). That branch bypasses the chi files
    entirely, so the ``chi0q_mode='flex'`` requirement -- which exists only
    because the scalar dynamic vertex INGESTS FLEX chiq_s/chiq_c -- does not
    apply to it; ``chi0q_mode`` is instead warned-and-ignored inside
    ``solve_dynamic`` (:func:`_warn_dynamic_bond_ignores_chi`). The remaining
    bond guards (real V, CoulombInter present, norb = 1) need the interaction
    files, so they run in ``solve_dynamic`` right after they are read
    (:func:`_validate_bond_interactions`) -- the same place the static path
    runs them.
    """
    eli = input_dict.get("eliashberg", {})
    # Goes through the strict validator, so a typo ("ture") is refused here
    # rather than quietly reading as false and disabling both the feature and
    # its guards.
    use_bond = _bond_bool_option(eli, "bond_channels", False)
    if not use_bond and eli.get("chi0q_mode") != "flex":
        raise ValueError(
            "eliashberg.frequency='dynamic' requires chi0q_mode='flex' "
            "(full-frequency chiq_s/chiq_c and a dressed green are only "
            "produced by the FLEX path); got chi0q_mode='{}'"
            .format(eli.get("chi0q_mode")))
    nmat = int(input_dict["mode"]["param"].get("Nmat", 1024))
    if nmat % 2 != 0:
        raise ValueError(
            "eliashberg.frequency='dynamic' requires an even Nmat "
            "(centered Matsubara grid); got Nmat={}".format(nmat))


def _warn_if_static_ignores_channel_flags(eli_param):
    """Warn that the channel-decomposition flags are inert on the static path.

    ``zero_chi_s``/``zero_chi_c`` split the DYNAMIC pairing vertex into its
    spin-/charge-fluctuation and instantaneous-bare pieces (see
    ``eliashberg_dynamic``). The static solver builds its kernel from the RPA/
    FLEX static susceptibility and never reads these flags, so a user who sets
    one with ``frequency='static'`` gets no decomposition. Warn loudly instead
    of silently ignoring the request.
    """
    from hwave.solver.backend import as_bool

    ignored = [name for name in ("zero_chi_s", "zero_chi_c")
               if as_bool(eli_param.get(name, False))]
    if ignored:
        logger.warning(
            "[eliashberg] %s set but frequency='static'; the channel-"
            "decomposition flags only apply to frequency='dynamic' and are "
            "ignored here.", ", ".join(ignored))


# ---------------------------------------------------------------------------
# Bond-resolved interaction channels (opt-in; spec S3-S5)
# ---------------------------------------------------------------------------

def _read_bond_config(eli_param):
    """Read the ``[eliashberg]`` bond-channel options (spec S5).

    Returns ``(use_bond_channels, bond_green, bond_max_shells,
    bond_memory_cap_gb, precondition_opts, bond_diagnostics)``; ``bond_green``
    is the path to an externally supplied Green function (``None`` = build the
    bare one from the transfer Hamiltonian), ``precondition_opts`` a kwargs
    dict for ``bond_channels.check_hermitian_preconditions`` holding only the
    keys the user actually set (so the function's own defaults stay
    authoritative), and ``bond_diagnostics`` the opt-in character analysis of
    the leading state (default ``False``).
    Every bond option set while ``bond_channels`` is off is IGNORED WITH A
    WARNING (they have no meaning on the default path) and never parsed, so a
    stale option cannot fail an otherwise valid run.

    Every value is validated centrally (``_bond_bool_option`` /
    ``_bond_float_option`` / ``_bond_int_option``): an unrecognised boolean
    spelling, a non-integral integer or an unparsable number raises naming the
    offending ``[eliashberg]`` key instead of being coerced into a DIFFERENT
    valid configuration.
    """
    use_bond = _bond_bool_option(eli_param, "bond_channels", False)
    if not use_bond:
        stale = [name for name in _BOND_OPTION_KEYS if name in eli_param]
        if stale:
            logger.warning(
                "[eliashberg] %s set but bond_channels=false; the bond-channel "
                "options only apply to bond_channels=true and are ignored "
                "here.", ", ".join(stale))
        return False, None, None, None, {}, False

    bond_diagnostics = _bond_bool_option(eli_param, "bond_diagnostics", False)

    bond_green = eli_param.get("bond_green")
    if bond_green is not None:
        bond_green = str(bond_green)
        if not bond_green:
            raise ValueError(
                "[eliashberg] bond_green must be a path to a Green-function "
                "npz file; got an empty string. Omit the key to build the "
                "bare Green function from the transfer Hamiltonian.")

    max_shells = _bond_int_option(eli_param, "bond_max_shells")
    if max_shells is not None and max_shells < 0:
        raise ValueError(
            "[eliashberg] bond_max_shells must be >= 0 (shell 0 = the "
            "on-site Delta r = 0 point), got {}".format(max_shells))
    cap_gb = _bond_float_option(eli_param, "bond_memory_cap_gb",
                                _BOND_MEMORY_CAP_GB)
    if not np.isfinite(cap_gb) or cap_gb <= 0.0:
        raise ValueError(
            "[eliashberg] bond_memory_cap_gb must be a positive finite number, "
            "got {!r}".format(eli_param.get("bond_memory_cap_gb")))

    # Hermiticity-precondition knobs (review fix M7): the tolerance is a
    # RELATIVE asymmetry of the symmetrized kernel K~ (both the residual and
    # its scale are K~-referred, review fix I1), and dense_limit forces the
    # randomized probe estimator on a grid too large for the exact dense
    # residual.
    pre_opts = {}
    for key, opt in (("bond_precondition_atol", "atol"),
                     ("bond_precondition_rtol", "rtol")):
        val = _bond_float_option(eli_param, key)
        if val is not None:
            if not np.isfinite(val) or val < 0.0:
                raise ValueError(
                    "[eliashberg] {} must be a non-negative finite number, "
                    "got {!r}".format(key, eli_param[key]))
            pre_opts[opt] = val
    dense_limit = _bond_int_option(eli_param, "bond_precondition_dense_limit")
    if dense_limit is not None:
        if dense_limit < 0:
            raise ValueError(
                "[eliashberg] bond_precondition_dense_limit must be >= 0 "
                "(0 = always use the randomized probe estimator), got "
                "{}".format(dense_limit))
        pre_opts["dense_limit"] = dense_limit
    return True, bond_green, max_shells, cap_gb, pre_opts, bond_diagnostics


def _validate_bond_prereqs(chi0q_mode, norb, interactions,
                           imag_tol=_BOND_IMAG_TOL):
    """Top-level guards for ``bond_channels=true`` (spec S5).

    Every one of these is a "no silent wrong result" refusal: the internal
    physics functions (``bond_channels.bond_bubble`` / ``bare_bond_vertices``
    / ``make_bond_kernel``) stay general so the unit and reduction tests can
    exercise multi-orbital and complex-bond physics directly; only this
    top-level entry point restricts what a production run may request.
    """
    if chi0q_mode == "flex":
        raise ValueError(
            "[eliashberg] bond_channels=true is incompatible with "
            "chi0q_mode='flex'. That mode INGESTS A FLEX chi (chiq_s/chiq_c "
            "npz), which is scalar/orbital-only: it carries no bond-resolved "
            "chi-bar, and the bond path cannot dress one from it. Making "
            "flex.py output a bond-resolved chi-bar (the conserving-FLEX bond "
            "self-consistency) is a deferred follow-up spec. NOTE that this "
            "does NOT prevent you from using a FLEX-dressed GREEN function: "
            "the bond path builds its own chi-bar from the Green function, so "
            "set [eliashberg] bond_green = '<path to green.npz>' together "
            "with chi0q_mode='calc' or 'load' (which are then unused anyway) "
            "to feed an external/FLEX green.")

    _validate_bond_interactions(norb, interactions, imag_tol=imag_tol)


def _validate_bond_interactions(norb, interactions, imag_tol=_BOND_IMAG_TOL):
    """The chi0q-mode-INDEPENDENT half of :func:`_validate_bond_prereqs`.

    Split out so the DYNAMIC bond path can reuse the identical model guards
    (real V, CoulombInter declared, norb = 1) without the static path's
    ``chi0q_mode='flex'`` refusal: the dynamic bond branch never reads a chi
    file at all (spec S3.4 "chi-file bypass"), so ``chi0q_mode`` is
    warned-and-ignored there (:func:`_warn_dynamic_bond_ignores_chi`) instead
    of being an error.
    """
    # "Truly absent from the config" is the documented trigger (spec S5): a
    # DECLARED term whose values are all zero is explicitly allowed, so that a
    # V sweep keeps a constant channel topology (B does not jump 1 <-> 5) at
    # its V=0 point.
    if "CoulombInter" not in interactions:
        raise ValueError(
            "[eliashberg] bond_channels=true requires a CoulombInter term: "
            "bond channels carry the inter-site Fock/Cooper structure of "
            "V_ab(R), so asking for them with no inter-site interaction "
            "declared is a configuration mistake. Declare the CoulombInter "
            "file (or a combined `Coulomb` file containing R!=0 entries) -- "
            "a declared-but-zero V is allowed and keeps the channel "
            "topology fixed across a V sweep -- or use bond_channels=false.")

    nonfinite = [(irvec, orbvec, value)
                 for (irvec, orbvec), value in interactions["CoulombInter"].items()
                 if not np.isfinite(complex(value).real)
                 or not np.isfinite(complex(value).imag)]
    if nonfinite:
        irvec, orbvec, value = nonfinite[0]
        raise ValueError(
            "[eliashberg] bond_channels=true requires every CoulombInter "
            "value to be finite, but V_{}{}({}) = {} is not (real or "
            "imaginary part is NaN/inf; {} such entries). A NaN/inf "
            "interaction value is a data/configuration error, not a "
            "physical model, and must not silently reach the bond "
            "machinery.".format(
                orbvec[0], orbvec[1], tuple(irvec), value, len(nonfinite)))

    bad = [(irvec, orbvec, value)
           for (irvec, orbvec), value in interactions["CoulombInter"].items()
           if abs(complex(value).imag) > imag_tol]
    if bad:
        irvec, orbvec, value = bad[0]
        raise ValueError(
            "[eliashberg] bond_channels=true supports only REAL inversion-"
            "symmetric CoulombInter in v1, but V_{}{}({}) = {} has a nonzero "
            "imaginary part ({} such entries). Complex bonds need the "
            "Hermitian Q(D+D^dag)Q particle-particle form, whose end-to-end "
            "validation is a deferred follow-up.".format(
                orbvec[0], orbvec[1], tuple(irvec), value, len(bad)))

    if norb > 1:
        raise ValueError(
            "[eliashberg] bond_channels=true is validated for single-orbital "
            "models only (norb=1); this model has norb={}. The multi-orbital "
            "equations are implemented and unit-tested, but their end-to-end "
            "numerical validation is a deferred follow-up spec, so the "
            "top-level solver refuses to emit unvalidated numbers. Note also "
            "that the bond path CORRECTS the existing multi-orbital collapsed-"
            "bond Fock approximation (spec S4.3), so its results would differ "
            "from bond_channels=false by design.".format(norb))


# Number of (N_q, ND, ND)-shaped complex128 arrays that the bond-channel
# resource preflight (below) assumes are alive simultaneously: chi_bar,
# S_bond, C_bond, chi_s, chi_c, the fluctuation vertex F_q, its two
# matrix-product temporaries (chi_bar@S_bond, chi_bar@C_bond, or equivalent),
# and the real-space vertex Gamma_hat = ifftn(F_q). Named so the preflight's
# byte estimate and this list can never silently drift apart.
_BOND_N_Q_ARRAYS = 9

# Number of FREQUENCY-RESOLVED (N_q, nmat, ND, ND) arrays alive simultaneously
# at ``bond_channels.dress_bond_dynamic``'s high-water mark, read off its body:
# the input ``chi_bar_w``, the broadcast identity ``I_mat``, the two RPA
# denominators ``mat_s``/``mat_c``, and the two solves' outputs
# ``chi_s_w``/``chi_c_w``. (``S_full``/``C_full`` are broadcast VIEWS of the
# static S_bond/C_bond and cost nothing; they are budgeted at their own static
# size in ``q_bytes``.) Named here so the dynamic preflight and that function's
# buffer discipline cannot silently drift apart.
_BOND_DRESS_DYN_ARRAYS = 6


def _bond_memory_estimate(norb, bond_set, Nx, Ny, Nz, nmat, dynamic_nmat=None):
    """Byte budget for the bond path, broken down by buffer family.

    Split out of :func:`_bond_resource_preflight` so the accounting can be
    asserted directly against a MEASURED peak (``tests/test_sc_bond.py``):
    a preflight that undercounts is worse than no preflight, because a run it
    waves through then OOMs anyway.

    Everything is complex128 (16 B). ``unit = N_q * nmat * 16`` is one
    ``norb``-pair-free frequency-grid buffer.

    ``q_bytes``
        ``_BOND_N_Q_ARRAYS`` q-resolved ``(N_q, ND, ND)`` arrays alive
        simultaneously in the vertex assembly -- ``chi_bar``, ``S_bond``,
        ``C_bond``, ``chi_s``, ``chi_c``, the fluctuation vertex ``F_q``, its
        two matrix-product temporaries, and ``Gamma_hat = ifftn(F_q)``.
        ``bond_bubble``'s output ``chi_bar`` is the first of these, so it is
        NOT counted again in ``bubble_bytes``.
    ``bubble_bytes``
        ``bond_bubble``'s simultaneous working set at its high-water mark, the
        input ``green`` excluded (that is ``green_bytes``):
        ``BOND_BUBBLE_N4_BUFFERS`` ``norb**4``-sized plus
        ``BOND_BUBBLE_N2_BUFFERS`` ``norb**2``-sized frequency-grid buffers.
        The buffer-by-buffer list, and the ``del``/in-place discipline in
        ``bond_bubble`` that holds the real peak down to it, are documented
        beside those two constants in ``hwave.solver.bond_channels``; they are
        imported from there so the two sides cannot silently desync.
    ``green_bytes``
        the ``green_kw`` array itself, ``(norb, norb, Nx, Ny, Nz, nmat)`` --
        i.e. ``bond_bubble``'s input, alive throughout.
    ``chi_bar_bytes``
        one ``(N_q, ND, ND)`` array; reported separately (it is already inside
        ``q_bytes``) because it is the one q-resolved array ``bond_bubble``
        itself allocates, so ``bubble_bytes + chi_bar_bytes`` is the budget a
        measurement of ``bond_bubble`` alone must fit into.
    ``vertex_bytes``
        ``bare_bond_vertices``'s non-q-resolved ``ND x ND`` working set --
        ``BARE_VERTEX_ND2_BUFFERS`` buffers (``P``, ``D``, ``Dh``, ``Id``,
        ``Q_s``, ``Q_t``, ``B_s``, ``B_t``, ``Vpp_s``, ``Vpp_t``), NOT
        multiplied by ``N_q`` since none of them carry a q-axis (unlike
        ``S_bond``/``C_bond``, which are already inside ``q_bytes``). This is
        alive concurrently with ``q_bytes`` (``chi_bar`` from ``bond_bubble``
        is still held) and ``green_bytes`` (``green_kw`` is still held), so it
        is additive in ``peak``. The buffer-by-buffer list is documented
        beside the constant in ``hwave.solver.bond_channels``; imported from
        there so the two sides cannot silently desync.
    ``peak``
        ``q_bytes + bubble_bytes + vertex_bytes + green_bytes`` -- what the
        cap is applied to.

    ``dynamic_nmat``
        ``None`` (default) reproduces the static estimate above byte-for-
        byte -- IDENTICAL to omitting the keyword entirely, so every
        existing call site is unaffected.

        When set to the dynamic path's bosonic-frequency grid size, the
        budget is rebuilt for the DYNAMIC chain
        (``bond_bubble_dynamic -> dress_bond_dynamic ->
        make_bond_kernel_dynamic``), whose every large array carries the
        frequency axis. ``unit_dyn = N_q * dynamic_nmat * ND**2 * 16`` is one
        frequency-resolved q-array (``chi_bar_w`` / ``chi_s_w`` / ``F_q_w``
        / ``F_rt`` are all exactly this size); ``unit_mv = N_q *
        dynamic_nmat * ND * 16`` is one matvec-class buffer.

        The peak is a PHASE MAXIMUM, not a sum: the four phases below do not
        coexist (each releases its temporaries before the next allocates --
        that ``del`` discipline is the documented contract of the constants
        imported here), so summing them would over-refuse legitimate runs by
        roughly a factor two on the dominant term. What IS summed is what
        lives across every phase: the Green function, the two static
        ``(N_q, ND, ND)`` bare vertices ``S_bond``/``C_bond``, and
        ``bare_bond_vertices``'s ``ND x ND`` scratch.

        * ``bubble_phase`` -- ``bond_bubble_dynamic``'s working set
          (``BOND_BUBBLE_DYN_*`` buffers) plus its ``chi_bar_w`` output.
        * ``dress_phase`` -- ``dress_bond_dynamic``'s live set at its peak:
          ``_BOND_DRESS_DYN_ARRAYS`` frequency-resolved q-arrays.
        * ``build_phase`` -- the kernel builder's vertex hoist
          (``BOND_KERNEL_DYNAMIC_VERTEX_BUFFERS``) while the caller still
          holds ``chi_s_w``/``chi_c_w``, plus the persistent ``G2p_w``.
        * ``matvec_phase`` -- the persistent ``F_rt``
          (``BOND_KERNEL_DYNAMIC_VERTEX_PERSISTENT``) and ``G2p_w``, plus
          ``BOND_KERNEL_DYNAMIC_MATVEC_BUFFERS`` matvec-class buffers, with
          ``chi_s_w``/``chi_c_w`` still owned by the caller.

        ``g2_bytes`` (``G2p_w``, size ``nd**2 * N_q * dynamic_nmat``) is an
        EXPLICIT line rather than being folded into the matvec constant: its
        ratio to one matvec buffer is ``nd / B``, so at ``B = 1`` it is
        ``nd`` times a matvec buffer, not a rounding error (the note beside
        ``BOND_KERNEL_DYNAMIC_MATVEC_BUFFERS`` in ``bond_channels``).

        **IR branch.** The IR kernel (``make_bond_kernel_dynamic_ir``) has no
        memory contract of its own yet, so it is deliberately preflighted
        with this same UNIFORM formula: its sampling-node axes are never
        longer than the uniform window (``axB.n_freq``, ``axF.n_freq`` <=
        ``nmat``) and its buffer structure mirrors the uniform builder's, so
        the uniform figure is a conservative UPPER bound. A dedicated IR
        contract (with its own measured constants) is a deferred follow-up.
    """
    nd = norb * norb
    B = int(bond_set.n_channels)
    ND = nd * B
    Nq = int(Nx) * int(Ny) * int(Nz)
    itemsize = 16  # complex128
    unit = Nq * int(nmat) * itemsize

    chi_bar_bytes_static = Nq * ND * ND * itemsize
    vertex_bytes = bond_channels.BARE_VERTEX_ND2_BUFFERS * ND * ND * itemsize
    green_bytes = (norb ** 2) * unit

    if dynamic_nmat is None:
        chi_bar_bytes = chi_bar_bytes_static
        n4_buffers = bond_channels.BOND_BUBBLE_N4_BUFFERS
        n2_buffers = bond_channels.BOND_BUBBLE_N2_BUFFERS
        q_bytes = _BOND_N_Q_ARRAYS * chi_bar_bytes_static
        bubble_bytes = (n4_buffers * norb ** 4
                        + n2_buffers * norb ** 2) * unit
        return {"nd": nd, "B": B, "ND": ND, "Nq": Nq,
                "chi_bar_bytes": chi_bar_bytes,
                "q_bytes": q_bytes,
                "bubble_bytes": bubble_bytes,
                "vertex_bytes": vertex_bytes,
                "green_bytes": green_bytes,
                "peak": q_bytes + bubble_bytes + vertex_bytes + green_bytes}

    nmat_dyn = int(dynamic_nmat)
    unit_dyn = chi_bar_bytes_static * nmat_dyn
    unit_mv = Nq * nmat_dyn * ND * itemsize
    g2_bytes = nd * nd * Nq * nmat_dyn * itemsize
    bubble_bytes = (bond_channels.BOND_BUBBLE_DYN_N4_BUFFERS * norb ** 4
                    + bond_channels.BOND_BUBBLE_DYN_N2_BUFFERS
                    * norb ** 2) * unit
    # S_bond and C_bond -- the only q-resolved arrays that stay STATIC-sized
    # on the dynamic path (they carry no frequency axis, spec S3.2).
    q_bytes = 2 * chi_bar_bytes_static

    bubble_phase = bubble_bytes + unit_dyn
    dress_phase = _BOND_DRESS_DYN_ARRAYS * unit_dyn
    build_phase = (2 * unit_dyn
                   + bond_channels.BOND_KERNEL_DYNAMIC_VERTEX_BUFFERS
                   * unit_dyn + g2_bytes)
    matvec_phase = (2 * unit_dyn
                    + bond_channels.BOND_KERNEL_DYNAMIC_VERTEX_PERSISTENT
                    * unit_dyn
                    + bond_channels.BOND_KERNEL_DYNAMIC_MATVEC_BUFFERS
                    * unit_mv + g2_bytes)
    phase_peak = max(bubble_phase, dress_phase, build_phase, matvec_phase)
    return {"nd": nd, "B": B, "ND": ND, "Nq": Nq,
            "chi_bar_bytes": unit_dyn,
            "q_bytes": q_bytes,
            "bubble_bytes": bubble_bytes,
            "vertex_bytes": vertex_bytes,
            "green_bytes": green_bytes,
            "g2_bytes": g2_bytes,
            "dress_bytes": dress_phase,
            "kernel_build_bytes": build_phase,
            "kernel_matvec_bytes": matvec_phase,
            "peak": (green_bytes + q_bytes + vertex_bytes + phase_peak)}


def _bond_resource_preflight(norb, bond_set, Nx, Ny, Nz, nmat, cap_gb,
                             dynamic_nmat=None):
    """Estimate the peak memory of the bond path and refuse to exceed the cap.

    Spec S3.2 "Resource guard": only here are ``nd``, ``N_q``, ``B`` and the
    dtype all known, so this is where the ``ND = nd*B`` blow-up is caught --
    never a silent runaway allocation. Called BEFORE ``green_kw`` is
    allocated (review fix I-2: previously the bond-channel dispatch branch
    called ``_calc_green`` first, so a large-Nmat run could OOM before the
    cap was ever consulted); the estimate includes the ``green_kw``
    allocation itself so the preflight is not blind to it.

    See :func:`_bond_memory_estimate` for the buffer-by-buffer accounting.
    ``dynamic_nmat`` selects the DYNAMIC (frequency-resolved) budget -- pass
    the bosonic grid size the dynamic chain will really run on; ``None``
    keeps the static budget byte-for-byte.
    """
    est = _bond_memory_estimate(norb, bond_set, Nx, Ny, Nz, nmat,
                                dynamic_nmat=dynamic_nmat)
    B, ND, Nq = est["B"], est["ND"], est["Nq"]
    q_bytes = est["q_bytes"]
    bubble_bytes = est["bubble_bytes"]
    vertex_bytes = est["vertex_bytes"]
    green_bytes = est["green_bytes"]
    peak = est["peak"]

    logger.info(
        "Bond-channel preflight%s: B = %d channels, ND = nd*B = %d, "
        "N_q = %d, estimated peak memory = %.3f GB (cap %.3f GB)",
        "" if dynamic_nmat is None
        else " (dynamic, nmat = {})".format(int(dynamic_nmat)),
        B, ND, Nq, peak / 1.0e9, cap_gb)
    if dynamic_nmat is not None:
        logger.info(
            "Bond-channel preflight (dynamic) component budget: chi-bar/"
            "dressing %.3f GB, kernel build %.3f GB, kernel matvec %.3f GB, "
            "pair bubble G2 %.3f GB, green %.3f GB",
            est["dress_bytes"] / 1.0e9, est["kernel_build_bytes"] / 1.0e9,
            est["kernel_matvec_bytes"] / 1.0e9, est["g2_bytes"] / 1.0e9,
            green_bytes / 1.0e9)

    if peak > cap_gb * 1.0e9:
        detail = ("{:.3f} GB of q-resolved ND x ND arrays".format(
            q_bytes / 1.0e9) if dynamic_nmat is None else
            "{:.3f} GB of static q-resolved ND x ND vertices, a "
            "frequency-resolved phase peak of {:.3f} GB (dressing {:.3f}, "
            "kernel build {:.3f}, kernel matvec {:.3f}, pair bubble G2 "
            "{:.3f})".format(
                q_bytes / 1.0e9,
                max(est["dress_bytes"], est["kernel_build_bytes"],
                    est["kernel_matvec_bytes"]) / 1.0e9,
                est["dress_bytes"] / 1.0e9,
                est["kernel_build_bytes"] / 1.0e9,
                est["kernel_matvec_bytes"] / 1.0e9,
                est["g2_bytes"] / 1.0e9))
        raise ValueError(
            "[eliashberg] bond_channels: estimated peak memory {:.3f} GB "
            "exceeds bond_memory_cap_gb = {:.3f} GB. The cost is driven by "
            "B = {} bond channels (Delta r = {}) giving ND = nd*B = {} on "
            "N_q = {} q-points ({}) plus "
            "{:.3f} GB of bond-bubble work buffers plus {:.3f} GB of bare-"
            "vertex ND x ND temporaries plus {:.3f} GB for the Green "
            "function itself. Reduce bond_max_shells (fewer channels), "
            "reduce the k-grid/Nmat, or raise bond_memory_cap_gb.".format(
                peak / 1.0e9, cap_gb, B, list(bond_set.delta_r), ND, Nq,
                detail, bubble_bytes / 1.0e9, vertex_bytes / 1.0e9,
                green_bytes / 1.0e9))
    return peak


def _build_bond_m0_blocks(bond_set, interactions, inter_k, norb,
                          kx_array, ky_array, kz_array):
    """Build the Case-2-CORRECTED ``m=0`` spin/charge blocks (spec S4.3 star).

    The existing ``_build_sc_matrices_all_q`` places the FULL Fourier sum
    ``V_ab(q)`` into the Fock ``(ab,ab)`` element of the local block -- a
    collapsed-bond approximation that assigns the whole Fock sum to the local
    orbital-coherence operator. The bond path corrects this: only the on-site
    ``R=0`` component ``V_ab(R=0)`` stays in the ``m=0`` Fock element, while
    every ``R!=0`` component moves to its own bond channel
    (``bond_set.v_bond[m]``, consumed by ``bare_bond_vertices``). The Hartree
    ``(aa,bb) += 2 V_ab(q)`` stays in the ``m=0`` charge block (correct as is)
    and is built from the SAME filtered interaction as the bond blocks -- the
    ``ResolvedInteractionSet`` is the single source of truth (spec S3.0), so a
    ``bond_max_shells`` cutoff filters Hartree, Fock and Cooper consistently.

    The on-site ``R=0`` Fock component is read from ``bond_set.v_onsite``,
    NOT from the raw ``interactions["CoulombInter"]`` dict (review fix I-1):
    ``v_onsite`` is Hermiticity-closed by ``resolve_interactions`` the same
    way the ``m != 0`` bond channels are, so a user declaring only
    ``V_01(0)`` (not ``V_10(0)``) still gets a reversal-closed on-site Fock
    block instead of a silently asymmetric one.

    All other interaction types (CoulombIntra, Hund, Exchange, Ising, PairHop)
    ride in the ``m=0`` block exactly as today.

    Parameters
    ----------
    bond_set : bond_channels.ResolvedInteractionSet
        The resolved (reversal-closed, shell-filtered) bond topology; its
        ``v_onsite`` field supplies the ``R=0`` Fock component.
    interactions : dict
        Real-space interaction dicts (``_read_interaction_files``); no longer
        read for ``CoulombInter`` here (kept in the signature for backward
        compatibility with existing call sites/tests) -- the ``R=0`` value
        comes from ``bond_set.v_onsite`` instead.
    inter_k : dict
        k-space interactions (``_build_interaction_k``); the CoulombInter
        entry is NOT used -- ``V(q)`` is rebuilt from the filtered bond set.
    norb : int
    kx_array, ky_array, kz_array : ndarray
        k-point arrays (the q-grid of the S/C blocks).

    Returns
    -------
    S0_q, C0_q : ndarray, shape (Nx, Ny, Nz, norb**2, norb**2)
        The Case-2-corrected local blocks, ready for
        ``bond_channels.bare_bond_vertices``.
    """
    Nx, Ny, Nz = len(kx_array), len(ky_array), len(kz_array)

    # Everything EXCEPT CoulombInter through the existing builder; the
    # CoulombInter Hartree/Fock placement is redone below from the resolved
    # (filtered) interaction so the two can never disagree.
    inter_k_local = {k: v for k, v in inter_k.items() if k != "CoulombInter"}
    S0, C0 = _build_sc_matrices_all_q(inter_k_local, norb, Nx, Ny, Nz)

    # On-site (R=0) CoulombInter: the only Fock component that stays local.
    # Sourced from bond_set.v_onsite -- Hermiticity-closed by
    # resolve_interactions the same way the bond (R!=0) channels are (review
    # fix I-1) -- instead of reading interactions["CoulombInter"] raw, which
    # bypassed that closure and consistency check.
    V_r0 = np.array(bond_set.v_onsite, dtype=complex, copy=True)

    # Filtered V_ab(q) = V_ab(R=0) + sum_{m>=1} V_ab(Delta r_m) e^{-i q.Delta r_m}
    # -- same e^{-iqR} convention as _build_interaction_k._to_k, so the m=0
    # block is element-wise consistent with the rest of sc.py. (For the real
    # inversion-symmetric interactions v1 accepts, V(q) is real and the phase
    # sign is immaterial.)
    kx_mesh, ky_mesh, kz_mesh = np.meshgrid(
        kx_array, ky_array, kz_array, indexing='ij')
    V_q = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
    V_q += V_r0[:, :, np.newaxis, np.newaxis, np.newaxis]
    for m in range(1, bond_set.n_channels):
        Rx, Ry, Rz = bond_set.delta_r[m]
        phase = np.exp(-1j * (kx_mesh * Rx + ky_mesh * Ry + kz_mesh * Rz))
        V_q += bond_set.v_bond[m][:, :, np.newaxis, np.newaxis, np.newaxis] \
            * phase

    for a in range(norb):
        for b in range(norb):
            # Hartree (aa,bb): the FULL q-dependent 2 V_ab(q) (spec S4.3).
            C0[:, :, :, a * norb + a, b * norb + b] += 2.0 * V_q[a, b]
            # Fock (ab,ab): ONLY the R=0 component (the Case-2 correction).
            iab = a * norb + b
            S0[:, :, :, iab, iab] += V_r0[a, b]
            C0[:, :, :, iab, iab] -= V_r0[a, b]

    return S0, C0


def _build_bond_operator(bond_set, green_kw, interactions, inter_k, geom_info,
                         norb, kx_array, ky_array, kz_array, beta,
                         pairing_type, *, bond_max_shells, bond_memory_cap_gb,
                         precondition_opts=None, green_source=None):
    """Run the bond-resolved physics chain and return the Eliashberg operator.

    bond_bubble -> Case-2 corrected m=0 blocks -> bare_bond_vertices ->
    dress_bond -> make_bond_kernel (spec S4.2-S4.5).

    ``bond_set`` must already be resolved (``bond_channels.
    resolve_interactions``) and pre-flighted (``_bond_resource_preflight``) by
    the caller BEFORE ``green_kw`` was allocated (review fix I-2: the
    preflight must run before the potentially-large Green-function
    allocation, not after -- see ``calc_eliashberg``'s bond-channel dispatch
    branch). This function therefore no longer resolves the interaction or
    runs the preflight itself.

    ``precondition_opts`` (from ``_read_bond_config``) is forwarded verbatim
    to ``bond_channels.check_hermitian_preconditions`` -- the
    ``bond_precondition_atol/rtol/dense_limit`` knobs of ``[eliashberg]``
    (review fix M7).

    ``green_source`` is provenance only: the path ``green_kw`` was loaded from
    (``[eliashberg] bond_green``), or ``None`` when the caller built the bare
    Green function from the transfer Hamiltonian. It selects the wording of
    the recorded approximation level -- the two are physically different
    approximations and the record must say which one ran.

    ``bond_max_shells`` and ``bond_memory_cap_gb`` are likewise PROVENANCE
    ONLY here: the resource preflight and the shell cutoff itself already
    happened in the caller (``resolve_interactions``/
    ``_bond_resource_preflight``, before ``bond_set``/``green_kw`` were
    built), so by the time this function runs, neither value feeds into the
    computed operator -- they are only echoed into the returned
    ``provenance`` dict. They are keyword-only (review fix: honest
    signature) precisely so a call site cannot mistake them for
    physics-affecting positional arguments. (Past review fix: an earlier
    ``nmat`` parameter was dropped from this signature for the same reason --
    it was never read in the body, since ``green_kw.shape[-1]`` already fixes
    the Matsubara grid the bond bubble is built on, and ``bond_green`` can
    make that grid differ from the ``nmat`` the caller resolved the preflight
    with.)

    Returns
    -------
    (A, vec_size) : tuple
        The ``scipy.sparse.linalg.LinearOperator`` kernel and its vector size,
        in the same convention as ``_make_kernel_operator``.
    provenance : dict
        Additive output-provenance record (spec S5).
    attribution : dict
        Context for the ``lambda = lambda^pp + lambda^fl`` decomposition of
        spec S4.5: the two kernel parts (``pp``/``fl``), the full kernel and
        the ``sqrt(GG)`` metric (``weight``), consumed by
        ``bond_channels.attribute_lambda`` once a gap has converged.
    """
    Nx, Ny, Nz = len(kx_array), len(ky_array), len(kz_array)

    logger.info("Computing the bond-resolved bubble chi-bar...")
    chi_bar = bond_channels.bond_bubble(green_kw, bond_set, beta)

    S0_q, C0_q = _build_bond_m0_blocks(bond_set, interactions, inter_k, norb,
                                       kx_array, ky_array, kz_array)
    S_bond, C_bond, Vpp_s, Vpp_t = bond_channels.bare_bond_vertices(
        bond_set, S0_q, C0_q, norb)
    chi_s, chi_c = bond_channels.dress_bond(chi_bar, S_bond, C_bond)

    logger.info("Building the bond-resolved Eliashberg operator "
                "(pairing_type=%s)...", pairing_type)
    A, A_fl, A_pp, vec_size = bond_channels.make_bond_kernel_parts(
        chi_s, chi_c, S_bond, C_bond, Vpp_s, Vpp_t, green_kw, bond_set,
        pairing_type, beta)

    # v1 runtime preconditions of the Hermitian static path (spec S4.5): the
    # pair weight w = GG must be real (Hermitian) and >= 0 and the symmetrized
    # kernel K~ = -sqrt(w) Gamma sqrt(w) Hermitian, else the reported lambda is
    # not a real eigenvalue/Rayleigh quotient. Violations RAISE -- there is no
    # non-Hermitian biorthogonal fallback in v1.
    weight = bond_channels.pair_weight(green_kw, beta)
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
        "bond_n_channels": int(bond_set.n_channels),
        "bond_delta_r": [list(dr) for dr in bond_set.delta_r],
        "bond_max_shells": ("all" if bond_max_shells is None
                            else int(bond_max_shells)),
        "bond_memory_cap_gb": float(bond_memory_cap_gb),
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
    if bond_set.n_channels == 1:
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
        The nmat metadata from the npz file, if present.

    Returns
    -------
    int or None
        Position of the zero bosonic frequency along the stored axis, or
        None when the file carries no usable metadata -- the caller must
        then slice the CENTER of the frequency axis it actually uses (only
        the caller can identify that axis reliably, e.g. for the 6D
        reference format).
    """
    if freq_index is None:
        logger.warning(
            "chi0q file '{}' has no freq_index metadata; using the center "
            "of the stored frequency axis as the static slice.".format(
                file_name))
        return None

    freq_index = np.asarray(freq_index).ravel()
    if freq_index.size == 0:
        raise ValueError(
            "chi0q file '{}' records no frequencies (matsubara_frequency="
            "\"none\"?); it cannot provide the static susceptibility."
            .format(file_name))
    if freq_index.size != nfreq:
        # Legacy FLEX files store the FULL grid but a restricted freq_index
        # (FLEX never applied the matsubara_frequency filter): the DATA axis
        # is authoritative, so fall back to the pre-metadata behavior.
        logger.warning(
            "chi0q file '{}': freq_index length {} does not match the "
            "frequency axis length {}; ignoring the metadata and using the "
            "center of the stored frequency axis as the static slice."
            .format(file_name, freq_index.size, nfreq))
        return None

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


def _reject_ir_native(data, file_name, hint):
    """Fail fast when a uniform-grid reader is handed sparse-node data
    (design ir-matsubara-stage3.md Sec. 4). MUST run before any freq_index
    metadata interpretation: the legacy fallbacks (missing freq_index ->
    center row) would otherwise silently read a sparse node as the static
    slice."""
    if is_ir_native(data):
        raise ValueError(
            "file '{}' holds sparse-IR node data "
            "(frequency_grid=sparse_ir_nodes); {}".format(file_name, hint))


def _load_chi0q(input_dict):
    """Load chi0q from NPZ file produced by H-wave RPA solver.

    Parameters
    ----------
    input_dict : dict
        Parsed TOML input dictionary.

    Returns
    -------
    chi0q : ndarray
        Bare susceptibility array.
    static_index : int
        Position of the zero bosonic frequency along the frequency axis
        (determined from the freq_index metadata).
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
    enable_spin_orbital = input_dict.get("mode", {}).get(
        "enable_spin_orbital", False)
    validate_chi0q_index_convention(data, enable_spin_orbital, file_name)
    chi0q = data["chi0q"]
    logger.info("chi0q shape: {}".format(chi0q.shape))

    # coeff_tail provenance check (issue #80): the Matsubara tail correction
    # changes chi0q at O(1) (lambda differed by 50% in the reproduction), so
    # a config whose effective coeff_tail differs from the file's would
    # recompute DIFFERENT physics under chi0q_mode="calc". Warn so load-vs-
    # calc comparisons are not silently inconsistent. Files without the key
    # (older builds) load silently as before.
    if "coeff_tail" in data:
        file_tail = float(data["coeff_tail"])
        config_tail = float(input_dict.get("mode", {}).get("param", {})
                            .get("coeff_tail", 0.0))
        if file_tail != config_tail:
            logger.warning(
                "chi0q file '{}' was produced with coeff_tail = {} but this "
                "config's effective value is coeff_tail = {}; results are "
                "NOT comparable with a chi0q_mode=\"calc\" run under this "
                "config. Set [mode.param] coeff_tail = {} to match the "
                "file.".format(file_name, file_tail, config_tail, file_tail))

    freq_index, file_nmat = _read_freq_meta(data)
    # Identify the frequency axis from the array LAYOUT, never from the
    # freq_index length (a restricted freq_index can coincidentally match
    # an orbital axis).  4D raw: axis 0.  8D ref: last axis.  6D is either
    # raw (nmat, nvol, norb^4; the last four axes are equal) or ref
    # (norb, norb, Nx, Ny, Nz, nmat; the first two axes are equal) --
    # disambiguate structurally, with the freq_index length only as the
    # tiebreaker for degenerate shapes.  Without metadata
    # _static_freq_position returns None and the caller slices the center
    # of the axis it actually uses.
    if freq_index is not None:
        nfi = np.asarray(freq_index).size
        if chi0q.ndim == 4:
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
            elif nfi == chi0q.shape[0]:
                nfreq = chi0q.shape[0]
            else:
                nfreq = chi0q.shape[-1]
        else:
            nfreq = chi0q.shape[-1]
    else:
        nfreq = 0  # unused: no metadata -> _static_freq_position gives None
    config_nmat = input_dict.get("mode", {}).get("param", {}).get("Nmat",
                                                                _DEFAULT_NMAT)
    static_index = _static_freq_position(freq_index, nfreq, config_nmat,
                                         file_name, file_nmat=file_nmat)
    if static_index is not None and static_index != nfreq // 2:
        logger.info("static (zero bosonic frequency) slice at index {} "
                    "of {}".format(static_index, nfreq))
    return chi0q, static_index


def _calc_chi0q_internal(input_dict, chi0q_tensor="auto",
                         precomputed_mu=None):
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

    Returns
    -------
    chi0q : ndarray
        Bare susceptibility array. Shape depends on chi0q_tensor mode:
        - reduced: (nmat, nvol, norb, norb)
        - general: (nmat, nvol, norb, norb, norb, norb)
    """
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
    if chi0q_tensor == "auto":
        # Use "general" when inter-orbital interactions are present,
        # because their S/C matrices have off-diagonal elements that
        # couple to chi0q off-diagonal components.
        # With CoulombIntra only, S is block-diagonal and reduced is exact.
        files = info_inputfile.get("interaction", {})
        has_interorbital = any(k in files for k in
                              ["Hund", "Exchange", "CoulombInter",
                               "Ising", "PairHop"])
        if has_interorbital:
            calc_scheme = "general"
        else:
            calc_scheme = "reduced"
    elif chi0q_tensor == "general":
        calc_scheme = "general"
    else:
        calc_scheme = "reduced"

    info_mode_rpa = dict(info_mode)
    info_mode_rpa["calc_scheme"] = calc_scheme

    logger.info("Computing chi0q internally (calc_scheme={})...".format(calc_scheme))

    # Read input files via QLMSkInput
    read_io = read_input_k.QLMSkInput(info_inputfile)
    ham_info = read_io.get_param("ham")

    # Create RPA solver (sets up Lattice, Interaction, parameters)
    solver = sol_rpa.RPA(ham_info, info_log, info_mode_rpa)

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

    chi0q = solver._calc_chi0q(green0, green0_tail, beta)

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
    files = input_dict["file"]["input"]["interaction"]
    path_to_input = files.get("path_to_input", "")

    geom_file = os.path.join(path_to_input, files["Geometry"])
    geom_info = wan90.read_geom(geom_file)
    logger.info("norb = {}".format(geom_info["norb"]))

    transfer_file = os.path.join(path_to_input, files["Transfer"])
    hr = wan90.read_w90(transfer_file)

    interaction_types = ["CoulombIntra", "CoulombInter", "Hund", "Exchange",
                        "Ising", "PairLift", "PairHop"]
    interactions = {}
    for itype in interaction_types:
        if itype in files:
            f = os.path.join(path_to_input, files[itype])
            logger.info("Reading {} from {}".format(itype, f))
            interactions[itype] = wan90.read_w90(f)

    # combined 'Coulomb' input: same decomposition as UHFk/RPA
    # (wan90.split_coulomb; r=0 diagonal -> CoulombIntra, rest -> CoulombInter)
    if "Coulomb" in files:
        if "CoulombIntra" in interactions or "CoulombInter" in interactions:
            raise ValueError(
                "Coulomb cannot be specified together with "
                "CoulombIntra or CoulombInter")
        f = os.path.join(path_to_input, files["Coulomb"])
        logger.info("Reading Coulomb from {}".format(f))
        coulomb_intra, coulomb_inter = wan90.split_coulomb(wan90.read_w90(f))
        if coulomb_intra:
            interactions["CoulombIntra"] = coulomb_intra
        if coulomb_inter:
            interactions["CoulombInter"] = coulomb_inter

    return geom_info, hr, interactions


# ---------------------------------------------------------------------------
# k-space construction
# ---------------------------------------------------------------------------

def _build_hamiltonian_k(kx_array, ky_array, kz_array, hr, norb):
    """Build Hamiltonian epsilon(k) from real-space transfer integrals.

    Parameters
    ----------
    kx_array, ky_array, kz_array : ndarray
        k-point arrays.
    hr : dict
        Transfer parameters from wan90.read_w90().
    norb : int
        Number of orbitals.

    Returns
    -------
    epsilon_k : ndarray
        Hamiltonian in k-space, shape (norb, norb, Nx, Ny, Nz).
    """
    Nx, Ny, Nz = len(kx_array), len(ky_array), len(kz_array)
    epsilon_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
    kx_mesh, ky_mesh, kz_mesh = np.meshgrid(
        kx_array, ky_array, kz_array, indexing='ij'
    )
    # Solver-core convention (rpa.py _make_ham_trans: tab_r[R,orb1,orb2] +
    # fftn == e^{-ikR}): epsilon[a,b](k) = sum_R t_R[a,b] e^{-ikR}. This keeps
    # sc-built quantities element-wise consistent with arrays loaded from
    # FLEX/RPA output files. (The previous [orb2,orb1] + e^{+ikR} form is the
    # orbital transpose at -k; identical for real hoppings, different for
    # complex Hermitian ones.)
    for (irvec, orbvec), value in hr.items():
        orb1, orb2 = orbvec
        Rx, Ry, Rz = irvec
        epsilon_k[orb1, orb2, :, :, :] += value * np.exp(
            -1j * (kx_mesh * Rx + ky_mesh * Ry + kz_mesh * Rz)
        )
    return epsilon_k


def _build_interaction_k(kx_array, ky_array, kz_array, interactions, norb):
    """Build interaction matrices in k-space.

    Supports CoulombIntra (U), CoulombInter (U'), Hund (J), and Exchange (J').

    Parameters
    ----------
    kx_array, ky_array, kz_array : ndarray
        k-point arrays.
    interactions : dict
        Interaction parameters keyed by type.
    norb : int
        Number of orbitals.

    Returns
    -------
    inter_k : dict
        Dictionary of interactions in k-space.
        Keys: "CoulombIntra", "CoulombInter", "Hund", "Exchange".
        Values: ndarray of shape (norb, norb, Nx, Ny, Nz).
    """
    Nx, Ny, Nz = len(kx_array), len(ky_array), len(kz_array)
    kx_mesh, ky_mesh, kz_mesh = np.meshgrid(
        kx_array, ky_array, kz_array, indexing='ij'
    )

    def _to_k(value_r):
        # same solver-core convention as _build_hamiltonian_k:
        # V[a,b](q) = sum_R V_R[a,b] e^{-iqR}
        val_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        for (irvec, orbvec), value in value_r.items():
            orb1, orb2 = orbvec
            Rx, Ry, Rz = irvec
            val_k[orb1, orb2, :, :, :] += value * np.exp(
                -1j * (kx_mesh * Rx + ky_mesh * Ry + kz_mesh * Rz)
            )
        return val_k

    inter_k = {}
    for itype in ["CoulombIntra", "CoulombInter", "Hund", "Exchange",
                  "Ising", "PairLift", "PairHop"]:
        if itype in interactions:
            inter_k[itype] = _to_k(interactions[itype])

    return inter_k


# ---------------------------------------------------------------------------
# Green's function
# ---------------------------------------------------------------------------

def _calc_eigenvalues(epsilon_k):
    """Diagonalize epsilon(k) at every k-point.

    Parameters
    ----------
    epsilon_k : ndarray
        Hamiltonian, shape (norb, norb, Nx, Ny, Nz).

    Returns
    -------
    eigenvalues : ndarray
        Shape (Nx, Ny, Nz, norb).
    eigenvectors : ndarray
        Shape (Nx, Ny, Nz, norb, norb).
    """
    norb = epsilon_k.shape[0]
    Nx, Ny, Nz = epsilon_k.shape[2], epsilon_k.shape[3], epsilon_k.shape[4]
    eigenvalues = np.zeros((Nx, Ny, Nz, norb))
    eigenvectors = np.zeros((Nx, Ny, Nz, norb, norb), dtype=complex)

    # Batch diagonalization: transpose to (Nx, Ny, Nz, norb, norb) for vectorized eigh
    eps_batch = epsilon_k.transpose(2, 3, 4, 0, 1)  # (Nx, Ny, Nz, norb, norb)
    eigenvalues, eigenvectors = np.linalg.eigh(eps_batch)

    return eigenvalues, eigenvectors


def _determine_mu(eigenvalues, beta, n_target, norb):
    """Determine chemical potential using bisection.

    The filling convention follows the reference code: n_target is the
    number of electrons per orbital per spin channel. For example,
    n_target=0.75 means 3/4 filling per spin.

    Parameters
    ----------
    eigenvalues : ndarray
        Shape (Nx, Ny, Nz, norb).
    beta : float
        Inverse temperature.
    n_target : float
        Target filling per orbital per spin (e.g. 0.75 for 3/4 filling).
    norb : int
        Number of orbitals.

    Returns
    -------
    mu : float
        Chemical potential.
    """
    Nx, Ny, Nz = eigenvalues.shape[:3]
    nvol = Nx * Ny * Nz

    def _calc_n(mu):
        x = beta * (eigenvalues - mu)
        fermi = np.where(x > 100, 0.0, np.where(x < -100, 1.0, 1.0 / (1.0 + np.exp(x))))
        total_n = np.sum(fermi)
        return float(total_n / nvol - n_target * norb)

    emin = np.min(eigenvalues)
    emax = np.max(eigenvalues)
    mu = bisect(_calc_n, emin - 10.0, emax + 10.0)
    return float(mu)


def _calc_green(eigenvalues, eigenvectors, mu, beta, nmat):
    """Construct Green's function G(k, iwn).

    Parameters
    ----------
    eigenvalues : ndarray
        Shape (Nx, Ny, Nz, norb).
    eigenvectors : ndarray
        Shape (Nx, Ny, Nz, norb, norb).
    mu : float
        Chemical potential.
    beta : float
        Inverse temperature.
    nmat : int
        Number of Matsubara frequencies.

    Returns
    -------
    green_kw : ndarray
        Shape (norb, norb, Nx, Ny, Nz, nmat).
    """
    Nx, Ny, Nz, norb = eigenvalues.shape
    iomega = np.array([(2.0 * i + 1.0 - nmat) * np.pi for i in range(nmat)]) / beta

    # Vectorized Green's function construction:
    # G_{ij}(k, iwn) = sum_m U_{im}(k) U*_{jm}(k) / (iwn - (e_m(k) - mu))

    # factor[kx,ky,kz,i,j,m] = U[kx,ky,kz,i,m] * conj(U[kx,ky,kz,j,m])
    factor = np.einsum('...im,...jm->...ijm', eigenvectors, np.conj(eigenvectors))
    # factor shape: (Nx, Ny, Nz, norb, norb, norb)

    # denom[kx,ky,kz,m,w] = 1 / (iwn_w - (e_m(k) - mu))
    xi = eigenvalues - mu  # (Nx, Ny, Nz, norb)
    denom = 1.0 / (1j * iomega[None, None, None, None, :] - xi[:, :, :, :, None])
    # denom shape: (Nx, Ny, Nz, norb, nmat)

    # G[kx,ky,kz,i,j,w] = sum_m factor[...,i,j,m] * denom[...,m,w]
    # Batched GEMM over the flattened spatial axes (C-order, reshape-only):
    #   (nv, ij, m) @ (nv, m, w) -> (nv, ij, w). numpy.einsum does NOT lower
    #   the '...ijm,...mw->...ijw' form to BLAS GEMM, so reshape to matmul.
    nv = Nx * Ny * Nz
    G = factor.reshape(nv, norb * norb, norb) @ denom.reshape(nv, norb, nmat)
    green_kw_tmp = G.reshape(Nx, Ny, Nz, norb, norb, nmat)
    # shape: (Nx, Ny, Nz, norb, norb, nmat)

    # Transpose to output convention: (norb, norb, Nx, Ny, Nz, nmat)
    green_kw = green_kw_tmp.transpose(3, 4, 0, 1, 2, 5)

    return green_kw


# ---------------------------------------------------------------------------
# RPA vertices
# ---------------------------------------------------------------------------

def _build_sc_matrices_all_q(inter_k, norb, Nx, Ny, Nz):
    """Build spin (S) and charge (C) interaction matrices for all q-points at once.

    Follows Kuroki et al., Eq.(5) in arXiv:0902.3691:
        S_{l1l2,l3l4}, C_{l1l2,l3l4} for multi-orbital systems.

    Parameters
    ----------
    inter_k : dict
        Interactions in k-space from _build_interaction_k.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        Grid dimensions.

    Returns
    -------
    S_all : ndarray
        Spin interaction matrices, shape (Nx, Ny, Nz, norb^2, norb^2).
    C_all : ndarray
        Charge interaction matrices, shape (Nx, Ny, Nz, norb^2, norb^2).
    """
    nd = norb * norb
    S_all = np.zeros((Nx, Ny, Nz, nd, nd), dtype=complex)
    C_all = np.zeros((Nx, Ny, Nz, nd, nd), dtype=complex)

    def _get(itype):
        if itype in inter_k:
            return inter_k[itype]  # (norb, norb, Nx, Ny, Nz)
        return None

    U_mat = _get("CoulombIntra")
    Up_mat = _get("CoulombInter")
    J_mat = _get("Hund")
    Jp_mat = _get("Exchange")
    I_mat = _get("Ising")
    PH_mat = _get("PairHop")

    # Build using precomputed index arrays to avoid Python loops
    # for small norb (1-3), the loop overhead is negligible;
    # for larger norb the vectorized approach helps
    l1_arr, l2_arr, l3_arr, l4_arr = np.meshgrid(
        np.arange(norb), np.arange(norb), np.arange(norb), np.arange(norb),
        indexing='ij')
    l1f, l2f, l3f, l4f = l1_arr.ravel(), l2_arr.ravel(), l3_arr.ravel(), l4_arr.ravel()
    idx12 = l1f * norb + l2f
    idx34 = l3f * norb + l4f

    # Case 1: l1 == l2 == l3 == l4
    mask1 = (l1f == l2f) & (l2f == l3f) & (l3f == l4f)
    if U_mat is not None and np.any(mask1):
        for i in np.where(mask1)[0]:
            _l = l1f[i]
            S_all[:, :, :, idx12[i], idx34[i]] = U_mat[_l, _l]
            C_all[:, :, :, idx12[i], idx34[i]] = U_mat[_l, _l]

    # Case 2: l1==l3, l2==l4, l1!=l2
    mask2 = (l1f == l3f) & (l2f == l4f) & (l1f != l2f)
    for i in np.where(mask2)[0]:
        s_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        c_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        _l1, _l2 = l1f[i], l2f[i]
        if Up_mat is not None:
            s_q += Up_mat[_l1, _l2]
            c_q -= Up_mat[_l1, _l2]
        if I_mat is not None:
            s_q -= I_mat[_l1, _l2]
            c_q -= I_mat[_l1, _l2]
        if J_mat is not None:
            c_q += J_mat[_l1, _l2]
        S_all[:, :, :, idx12[i], idx34[i]] = s_q
        C_all[:, :, :, idx12[i], idx34[i]] = c_q

    # Case 3: l1==l2, l3==l4, l1!=l3
    mask3 = (l1f == l2f) & (l3f == l4f) & (l1f != l3f)
    for i in np.where(mask3)[0]:
        s_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        c_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        _l1, _l3 = l1f[i], l3f[i]
        if J_mat is not None:
            s_q += J_mat[_l1, _l3]
            c_q -= J_mat[_l1, _l3]
        if I_mat is not None:
            s_q -= 2.0 * I_mat[_l1, _l3]
        if Up_mat is not None:
            c_q += 2.0 * Up_mat[_l1, _l3]
        S_all[:, :, :, idx12[i], idx34[i]] = s_q
        C_all[:, :, :, idx12[i], idx34[i]] = c_q

    # Case 4: l1==l4, l2==l3, l1!=l2
    mask4 = (l1f == l4f) & (l2f == l3f) & (l1f != l2f)
    for i in np.where(mask4)[0]:
        s_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        _l1, _l2 = l1f[i], l2f[i]
        if Jp_mat is not None:
            s_q += Jp_mat[_l1, _l2]
        if PH_mat is not None:
            s_q += PH_mat[_l1, _l2]
        S_all[:, :, :, idx12[i], idx34[i]] = s_q
        C_all[:, :, :, idx12[i], idx34[i]] = s_q  # S = C for this channel

    return S_all, C_all


def _build_sc_matrices(inter_k, norb, ix, iy, iz):
    """Build spin (S) and charge (C) interaction matrices at a given q-point.

    Follows Kuroki et al., Eq.(5) in arXiv:0902.3691:
        S_{l1l2,l3l4}, C_{l1l2,l3l4} for multi-orbital systems.

    The composite index maps as (l1,l2) -> l1*norb + l2,
    giving norb^2 x norb^2 matrices.

    Parameters
    ----------
    inter_k : dict
        Interactions in k-space from _build_interaction_k.
    norb : int
        Number of orbitals.
    ix, iy, iz : int
        q-point indices.

    Returns
    -------
    S_mat : ndarray
        Spin interaction matrix, shape (norb^2, norb^2).
    C_mat : ndarray
        Charge interaction matrix, shape (norb^2, norb^2).
    """
    nd = norb * norb
    S_mat = np.zeros((nd, nd), dtype=complex)
    C_mat = np.zeros((nd, nd), dtype=complex)

    # Extract interaction values at this q-point
    def _get(itype):
        if itype in inter_k:
            return inter_k[itype][:, :, ix, iy, iz]
        return np.zeros((norb, norb), dtype=complex)

    U_mat = _get("CoulombIntra")    # U_mm (intra-orbital)
    Up_mat = _get("CoulombInter")   # U'_mm' (inter-orbital)
    J_mat = _get("Hund")            # J_mm' (Hund's coupling)
    Jp_mat = _get("Exchange")       # J'_mm' (pair-hopping)
    I_mat = _get("Ising")           # I_mm' (Ising S^z S^z)
    PH_mat = _get("PairHop")        # P_mm' (pair hopping)

    for l1 in range(norb):
        for l2 in range(norb):
            idx12 = l1 * norb + l2
            for l3 in range(norb):
                for l4 in range(norb):
                    idx34 = l3 * norb + l4

                    s_val = 0.0 + 0.0j
                    c_val = 0.0 + 0.0j

                    if l1 == l2 == l3 == l4:
                        # Same orbital: S = U, C = U
                        s_val = U_mat[l1, l1]
                        c_val = U_mat[l1, l1]
                    elif l1 == l3 and l2 == l4 and l1 != l2:
                        # l1=l3 != l2=l4 (cross): S = U' - I, C = -U' + J - I
                        s_val = Up_mat[l1, l2] - I_mat[l1, l2]
                        c_val = (-Up_mat[l1, l2] + J_mat[l1, l2]
                                 - I_mat[l1, l2])
                    elif l1 == l2 and l3 == l4 and l1 != l3:
                        # l1=l2 != l3=l4 (dens): S = J - 2I, C = 2U' - J
                        s_val = J_mat[l1, l3] - 2.0 * I_mat[l1, l3]
                        c_val = 2.0 * Up_mat[l1, l3] - J_mat[l1, l3]
                    elif l1 == l4 and l2 == l3 and l1 != l2:
                        # l1=l4 != l2=l3 (exch): S = J' + P, C = J' + P
                        s_val = Jp_mat[l1, l2] + PH_mat[l1, l2]
                        c_val = Jp_mat[l1, l2] + PH_mat[l1, l2]

                    S_mat[idx12, idx34] = s_val
                    C_mat[idx12, idx34] = c_val

    return S_mat, C_mat


def _compute_vertices(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                      pairing_type="singlet", static_index=None):
    """Compute effective pairing interaction V(q).

    Supports two modes:
    - Simple mode (only CoulombIntra/Inter, 2-index chi0q):
      Uses Wc=U+2V, Ws=-U formulation.
    - General mode (with Hund/Exchange, or 4-index chi0q):
      Uses S,C matrices from Kuroki et al.

    When 4-index chi0q is provided (8D array), general mode is always used
    because the full orbital tensor structure must be preserved.

    Pairing types:
    - singlet: V^s = (3/2) S chi_s S - (1/2) C chi_c C + (1/2)(S + C)
    - triplet: V^t = -(1/2) S chi_s S - (1/2) C chi_c C + (1/2)(C - S)

    Parameters
    ----------
    chi0q : ndarray
        Bare susceptibility.
        - 2-index: shape (norb, norb, Nx, Ny, Nz, nmat)
        - 4-index: shape (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
    inter_k : dict
        Interactions in k-space from _build_interaction_k.
    norb, Nx, Ny, Nz, nmat : int
        System parameters.
    pairing_type : str
        "singlet" or "triplet". Default "singlet".

    Returns
    -------
    Result depends on mode:
        Simple mode: tuple (Pc_q, Ps_q), each shape (norb, norb, Nx, Ny, Nz).
        General mode: Vs_q, shape (norb, norb, norb, norb, Nx, Ny, Nz).
    """
    has_interorbital_vertex = any(k in inter_k for k in
                                  ["Hund", "Exchange", "Ising", "PairHop"])
    chi0q_is_4index = (chi0q.ndim == 8)

    # PairLift does not enter the spin/charge (particle-hole) pairing vertex in
    # either mode (its S=C=0 contribution, verified against the full 4-index
    # RPA), so a configured PairLift term is silently inert. Warn the user.
    if "PairLift" in inter_k:
        logger.warning(
            "PairLift is configured but does not contribute to the S/C pairing "
            "vertex (S=C=0); it is ignored in the Eliashberg calculation.")

    if chi0q_is_4index or has_interorbital_vertex:
        # General mode: 4-index S,C matrices
        # Required when 4-index chi0q is available or Hund/Exchange present
        return _compute_vertices_general(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                                         pairing_type=pairing_type,
                                         static_index=static_index)
    else:
        # Simple mode: backward compatible with original implementation
        # Only used for 2-index chi0q without Hund/Exchange
        if "CoulombInter" in inter_k and norb > 1:
            logger.warning(
                "Multi-orbital CoulombInter in the simple (2-index) vertex mode "
                "uses the reduced density-density approximation: inter-orbital "
                "cross-channel contributions are dropped. Provide a 4-index "
                "chi0q (general mode) for the full Kuroki S/C treatment.")
        return _compute_vertices_simple(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                                        pairing_type=pairing_type,
                                        static_index=static_index)


def _compute_vertices_simple(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                             pairing_type="singlet", static_index=None):
    """Compute vertices using simple Wc=U+2V, Ws=-U formulation.

    For singlet:
        Pc = (Wc+Ws)/2 - (1/2) Wc chi_c Wc
        Ps = -Ws + (3/2) Ws chi_s Ws
    For triplet:
        Pc = (Wc-Ws)/2 - (1/2) Wc chi_c Wc   (charge part same sign change)
        Ps = Ws - (1/2) Ws chi_s Ws            (spin part sign flip)

    Returns Pc_q and Ps_q separately for backward compatibility.

    Returns
    -------
    Pc_q : ndarray
        Charge vertex, shape (norb, norb, Nx, Ny, Nz).
    Ps_q : ndarray
        Spin vertex, shape (norb, norb, Nx, Ny, Nz).
    """
    U_k = inter_k.get("CoulombIntra", np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex))
    V_k = inter_k.get("CoulombInter", np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex))

    # Transpose to batch dimension first: (Nx, Ny, Nz, norb, norb)
    Wc = (U_k + 2.0 * V_k).transpose(2, 3, 4, 0, 1).copy()
    Ws = (-U_k).transpose(2, 3, 4, 0, 1).copy()

    # chi0 at static limit: (Nx, Ny, Nz, norb, norb)
    # (zero bosonic frequency; nmat//2 only for a full frequency grid)
    si = nmat // 2 if static_index is None else static_index
    chi0_static = chi0q[:, :, :, :, :, si].transpose(2, 3, 4, 0, 1).copy()

    I_mat = np.broadcast_to(np.eye(norb, dtype=complex), (Nx, Ny, Nz, norb, norb)).copy()

    # Batched solve
    mat_s = I_mat + chi0_static @ Ws
    mat_c = I_mat + chi0_static @ Wc

    chis = np.linalg.solve(mat_s, chi0_static)
    chic = np.linalg.solve(mat_c, chi0_static)

    WsChisWs = Ws @ chis @ Ws
    WcChicWc = Wc @ chic @ Wc

    if pairing_type == "singlet":
        Pc_all = (Wc + Ws) / 2.0 - 0.5 * WcChicWc
        Ps_all = -Ws + 1.5 * WsChisWs
    elif pairing_type == "triplet":
        Pc_all = (Wc - Ws) / 2.0 - 0.5 * WcChicWc
        Ps_all = Ws - 0.5 * WsChisWs
    else:
        raise ValueError("Unknown pairing_type: '{}'. Use 'singlet' or 'triplet'.".format(
            pairing_type))

    # Transpose back: (Nx, Ny, Nz, norb, norb) -> (norb, norb, Nx, Ny, Nz)
    Pc_q = Pc_all.transpose(3, 4, 0, 1, 2)
    Ps_q = Ps_all.transpose(3, 4, 0, 1, 2)

    return Pc_q, Ps_q


def _compute_vertices_general(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                              pairing_type="singlet", static_index=None):
    """Compute effective pairing interaction using general S,C matrices.

    For singlet (Kuroki et al., arXiv:0902.3691, Eq.(6)):
        V^s = (3/2) S chi_s S - (1/2) C chi_c C + (1/2)(S + C)

    For triplet (Takimoto et al., PRB 69, 104504):
        V^t = -(1/2) S chi_s S - (1/2) C chi_c C + (1/2)(C - S)

    Supports both 2-index (reduced) and 4-index (general) chi0q:
    - 2-index: shape (norb, norb, Nx, Ny, Nz, nmat)
      Expanded to diagonal of (norb^2, norb^2) matrix.
    - 4-index: shape (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
      Used directly as (norb^2, norb^2) matrix. This is the correct
      treatment for multi-orbital systems.

    Returns
    -------
    V_q : ndarray
        Effective pairing interaction, shape (norb, norb, norb, norb, Nx, Ny, Nz).
    """
    nd = norb * norb
    chi0q_is_4index = (chi0q.ndim == 8)

    # Build S, C matrices for all q-points at once: (Nx, Ny, Nz, nd, nd)
    S_all, C_all = _build_sc_matrices_all_q(inter_k, norb, Nx, Ny, Nz)

    # Extract chi0 at static limit for all q-points
    # (zero bosonic frequency; nmat//2 only for a full frequency grid)
    si = nmat // 2 if static_index is None else static_index
    if chi0q_is_4index:
        # (norb, norb, norb, norb, Nx, Ny, Nz, nmat) -> (Nx, Ny, Nz, nd, nd)
        chi0_static = chi0q[:, :, :, :, :, :, :, si].reshape(
            nd, nd, Nx, Ny, Nz).transpose(2, 3, 4, 0, 1).copy()
    else:
        # (norb, norb, Nx, Ny, Nz, nmat) -> expand to (Nx, Ny, Nz, nd, nd)
        chi0_2d = chi0q[:, :, :, :, :, si].transpose(2, 3, 4, 0, 1).copy()
        # chi0_2d shape: (Nx, Ny, Nz, norb, norb)
        if norb == 1:
            chi0_static = chi0_2d.reshape(Nx, Ny, Nz, 1, 1)
        else:
            # Expand: chi0_{l1*norb+l2, l3*norb+l2} = chi0_2d[l1, l3]
            chi0_static = np.zeros((Nx, Ny, Nz, nd, nd), dtype=complex)
            for l2 in range(norb):
                chi0_static[:, :, :,
                            l2::norb,
                            l2::norb] = chi0_2d

    # Batched RPA solve for all q-points simultaneously
    # chi_s = [I - chi0 @ S]^{-1} @ chi0
    # chi_c = [I + chi0 @ C]^{-1} @ chi0
    # The solve itself lives in bond_channels.dress_bond, the SINGLE shared
    # dressing helper: the bond path runs the identical algebra at the enlarged
    # bond-major size ND = nd*B, so factoring it out lets the bond reduction
    # test exercise the real production code path instead of a mirror copy.
    # dress_bond performs exactly the same operations in the same order, so
    # this call is bit-identical to the historical inline body.
    #
    # cond_tol=None DISABLES dress_bond's near-singularity guard here: that
    # guard is a BOND-PATH-ONLY policy. This is the legacy default path, which
    # has always simply called np.linalg.solve; an invertible but poorly
    # conditioned RPA denominator must keep returning its (large but finite)
    # result rather than becoming a ValueError for existing users.
    chis, chic = bond_channels.dress_bond(chi0_static, S_all, C_all,
                                          cond_tol=None)

    SChisS = S_all @ chis @ S_all
    CChicC = C_all @ chic @ C_all

    if pairing_type == "singlet":
        Vs_all = 1.5 * SChisS - 0.5 * CChicC + 0.5 * (S_all + C_all)
    elif pairing_type == "triplet":
        Vs_all = -0.5 * SChisS - 0.5 * CChicC + 0.5 * (C_all - S_all)
    else:
        raise ValueError("Unknown pairing_type: '{}'. Use 'singlet' or 'triplet'.".format(
            pairing_type))

    # Reshape (Nx, Ny, Nz, nd, nd) -> (norb, norb, norb, norb, Nx, Ny, Nz)
    Vs_q = Vs_all.reshape(Nx, Ny, Nz, norb, norb, norb, norb).transpose(3, 4, 5, 6, 0, 1, 2)

    return Vs_q


def _compute_vertices_flex(chis, chic, inter_k, norb, Nx, Ny, Nz,
                           pairing_type="singlet", convention="kuroki"):
    """Compute pairing vertex from pre-computed FLEX susceptibilities.

    Uses the same formula as _compute_vertices_general but takes
    chi_s and chi_c directly instead of computing them from chi0q via RPA.
    This allows using dressed susceptibilities from the FLEX solver.

    ``convention`` selects the S/C interaction matrices applied to the
    susceptibilities: "kuroki" (default, arXiv:0902.3691, used by the reduced
    FLEX path and the rest of the Eliashberg solver) or "myo"
    (cond-mat/0407094). The two differ ONLY in the charge ``C(ab,ab) = -U'+2J``
    vs ``-U'+J`` element; they share the native [a,c,b,d] orbital-pair index
    layout. Susceptibilities produced by the general (full-vertex) FLEX path
    were computed with the MYO S/C and MUST be paired with ``convention='myo'``
    so the vertex stays self-consistent (mixing them with Kuroki S/C silently
    changes the physics in the C(ab,ab) channel). ``convention`` is a S/C-charge
    selector, not an index-layout flag: ``chis``/``chic`` are expected in the
    native [a,c,b,d] orbital-pair order regardless (issue #78).

    Parameters
    ----------
    chis : ndarray
        Spin susceptibility at static limit, shape (Nx, Ny, Nz, nd, nd)
        where nd = norb^2.
    chic : ndarray
        Charge susceptibility at static limit, same shape.
    inter_k : dict
        Interactions in k-space from _build_interaction_k.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        Grid dimensions.
    pairing_type : str
        "singlet" or "triplet".

    Returns
    -------
    Vs_q : ndarray
        Pairing vertex, shape (norb, norb, norb, norb, Nx, Ny, Nz).
    """
    nd = norb * norb

    if "PairLift" in inter_k:
        logger.warning(
            "PairLift is configured but does not contribute to the S/C pairing "
            "vertex (S=C=0); it is ignored in the Eliashberg calculation.")

    conv = convention.lower()
    if conv == "myo":
        from hwave.solver._sc_matrices_myo import build_sc_matrices_myo
        S_all, C_all = build_sc_matrices_myo(inter_k, norb, Nx, Ny, Nz)
    elif conv == "kuroki":
        S_all, C_all = _build_sc_matrices_all_q(inter_k, norb, Nx, Ny, Nz)
    else:
        raise ValueError(
            "Unknown convention: '{}'. Use 'kuroki' or 'myo'.".format(convention))

    SChisS = S_all @ chis @ S_all
    CChicC = C_all @ chic @ C_all

    if pairing_type == "singlet":
        Vs_all = 1.5 * SChisS - 0.5 * CChicC + 0.5 * (S_all + C_all)
    elif pairing_type == "triplet":
        Vs_all = -0.5 * SChisS - 0.5 * CChicC + 0.5 * (C_all - S_all)
    else:
        raise ValueError("Unknown pairing_type: '{}'".format(pairing_type))

    Vs_q = Vs_all.reshape(Nx, Ny, Nz, norb, norb, norb, norb).transpose(
        3, 4, 5, 6, 0, 1, 2)

    return Vs_q


def _resolve_flex_paths(input_dict):
    """Resolve the FLEX output directory and chi_s/chi_c/green file paths.

    Shared by `_load_flex_susceptibilities_full` (which reads the files) and
    `_load_flex_susceptibilities` (which re-derives the static frequency
    index from the same files' metadata), so the two never disagree about
    which files are in play.
    """
    file_input = input_dict.get("file", {}).get("input", {})
    flex_dir = file_input.get("path_to_flex_output",
                              input_dict.get("file", {}).get("output", {}).get(
                                  "path_to_output", "output"))

    eli_param = input_dict.get("eliashberg", {})
    chi_s_file = eli_param.get("flex_chi_s", "chiq_s.npz")
    chi_c_file = eli_param.get("flex_chi_c", "chiq_c.npz")
    green_file = eli_param.get("flex_green", "green.npz")

    chi_s_path = os.path.join(flex_dir, chi_s_file)
    chi_c_path = os.path.join(flex_dir, chi_c_file)
    green_path = os.path.join(flex_dir, green_file)
    return chi_s_path, chi_c_path, green_path


def _load_flex_susceptibilities_full(input_dict, norb, Nx, Ny, Nz,
                                     allow_ir=False):
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
        susceptibilities with (they differ only in the C(ab,ab) charge value)
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
        chi_s_raw, chi_c_raw, chi_convention = _read_flex_chi_raw(input_dict)
        ir_meta = None
        green_w = _load_flex_green(input_dict, norb, Nx, Ny, Nz)
    else:
        chi_s_raw, chi_c_raw, chi_convention, ir_meta = _read_flex_chi_raw(
            input_dict, allow_ir=True)
        green_w, green_meta = _load_flex_green(input_dict, norb, Nx, Ny, Nz,
                                               allow_ir=True)
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


def _read_flex_chi_raw(input_dict, allow_ir=False):
    """Read the raw FLEX chi_s / chi_c NPZ arrays and their orbital convention.

    Returns ``(chi_s_raw, chi_c_raw, chi_convention)`` in the H-wave layout
    ``(nfreq, nvol, nd, nd)`` -- no reshape/expansion, so callers that only
    need one static frequency can slice before expanding.

    With ``allow_ir=True`` (the dynamic Eliashberg caller) the return value
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
    if not allow_ir:
        _reject_ir_native(data_s, chi_s_path, _STATIC_IR_HINT)
    chi_s_raw = data_s["chiq_s"] if "chiq_s" in data_s else data_s["chiq"]
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
    if not allow_ir:
        _reject_ir_native(data_c, chi_c_path, _STATIC_IR_HINT)
    chi_c_raw = data_c["chiq_c"] if "chiq_c" in data_c else data_c["chiq"]
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
        return chi_s_raw, chi_c_raw, chi_convention

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
    return chi_s_raw, chi_c_raw, chi_convention, ir_meta


def _expand_flex_chi(chi_raw, norb, Nx, Ny, Nz, convention):
    """Reshape H-wave chi ``(nfreq, nvol, nd, nd)`` to
    ``(nfreq, Nx, Ny, Nz, nd, nd)`` in the ``nd = norb^2`` Eliashberg space,
    resolving the spin-orbital-vs-orbital-pair layout from ``convention``.

    The two FLEX conventions have DIFFERENT physical layouts that shape alone
    cannot always tell apart:

    - ``"myo"`` (general full-vertex FLEX) is already in orbital-pair space
      ``nd_chi = norb^2``; passed through unchanged.
    - ``"kuroki"`` (reduced / squashed FLEX) is in spin-orbital reduced space
      ``nd_chi = norb*ns`` (spin-block ordered ``s*norb + a``); the spin-up
      orbital block ``[:norb, :norb]`` is extracted and diagonally expanded to
      ``norb^2 x norb^2``.

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
    chi_full = chi_raw.reshape(nfreq, Nx, Ny, Nz, -1)
    nd_chi = int(np.sqrt(chi_full.shape[-1]))
    chi_full = chi_full.reshape(nfreq, Nx, Ny, Nz, nd_chi, nd_chi)

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
        raise ValueError(
            "FLEX chi dimension nd_chi={} matches neither the orbital-pair "
            "size norb^2={} nor the spin-orbital size norb*ns={}.".format(
                nd_chi, nd, nd_so))

    if not is_spin_orbital:
        return chi_full

    # Spin-orbital reduced (spin-block ordered s*norb+a): extract the spin-up
    # orbital block and scatter it diagonally into norb^2 x norb^2.
    chi_orb = chi_full[:, :, :, :, :norb, :norb]
    out = np.zeros((nfreq, Nx, Ny, Nz, nd, nd), dtype=complex)
    for l2 in range(norb):
        out[:, :, :, :, l2::norb, l2::norb] = chi_orb
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
    _, nfreq, nvol_g, norb1, norb2 = green_raw.shape
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
    unknown even to that internal dispatcher.
    """
    version = np.lib.format.read_magic(fh)
    reader = getattr(
        np.lib.format, "read_array_header_{}_{}".format(*version), None)
    if reader is not None:
        shape, _, _ = reader(fh)
        return shape
    generic = getattr(np.lib.format, "_read_array_header", None)
    if generic is None:
        raise _UnsupportedNpyHeaderVersion(
            "numpy.lib.format has no header reader for NPY format version "
            "{!r}".format(version))
    try:
        shape, _, _ = generic(fh, version)
    except ValueError as exc:
        raise _UnsupportedNpyHeaderVersion(
            "cannot parse NPY header version {!r}: {}".format(version, exc))
    return shape


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


def _validate_bond_green_nmat(nfreq, green_path, label="bond_green"):
    """Reject a bond Green file whose frequency count cannot define a grid.

    The file's ``nfreq`` OVERRIDES the configured ``Nmat`` on the bond path
    (the Green function defines the grid the bubble is built on), so it has to
    meet the same requirement the configured value does: POSITIVE and EVEN.
    ``bond_bubble`` builds a CENTERED fermionic grid
    ``iomega_n = (2n + 1 - nmat) pi/beta``; with an odd ``nmat`` the ``n =
    (nmat-1)/2`` point sits exactly at ``iomega = 0``, which is not a fermionic
    Matsubara frequency, and the bubble would process it silently -- an
    invalid susceptibility with no error at all. Mirrors the dynamic path's
    even-Nmat requirement (:func:`_validate_dynamic_prerequisites`).

    Called BEFORE the resource preflight and before the array is materialised,
    so a bad grid never sizes -- let alone allocates -- anything large.
    """
    nfreq = int(nfreq)
    if nfreq <= 0 or nfreq % 2 != 0:
        raise ValueError(
            "[eliashberg] {} {} carries Nmat={} Matsubara frequencies, but "
            "the bond path requires a POSITIVE EVEN count (the file's Nmat "
            "overrides [mode.param] Nmat, and the bond bubble is built on the "
            "centered fermionic grid iomega_n = (2n+1-Nmat)*pi/beta, which an "
            "odd count would place a spurious zero frequency on). Regenerate "
            "the Green function with an even Nmat.".format(
                label, green_path, nfreq))
    return nfreq


def _load_flex_green(input_dict, norb, Nx, Ny, Nz, allow_ir=False):
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
    if not allow_ir:
        _reject_ir_native(data_g, green_path, _STATIC_IR_HINT)
    green_raw = data_g["green"]
    # H-wave format: (nblock, nfreq, nvol, norb, norb)
    nblock, nmat_g, nvol, norb1, norb2 = green_raw.shape
    # Convert to sc.py format: (norb, norb, Nx, Ny, Nz, nfreq)
    green = green_raw[0].reshape(
        nmat_g, Nx, Ny, Nz, norb, norb
    ).transpose(4, 5, 1, 2, 3, 0).copy()
    if not allow_ir:
        return green
    meta = ir_native_meta(data_g) if is_ir_native(data_g) else None
    return green, meta


def _load_flex_susceptibilities(input_dict, norb, Nx, Ny, Nz):
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
    chi_s_raw, chi_c_raw, chi_convention = _read_flex_chi_raw(input_dict)

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
    chis = _expand_flex_chi(chi_s_raw[center_s:center_s + 1],
                            norb, Nx, Ny, Nz, chi_convention)[0]
    chic = _expand_flex_chi(chi_c_raw[center_c:center_c + 1],
                            norb, Nx, Ny, Nz, chi_convention)[0]

    green_dressed = _load_flex_green(input_dict, norb, Nx, Ny, Nz)
    return chis, chic, green_dressed, chi_convention


# ---------------------------------------------------------------------------
# G2 and Eliashberg kernel
# ---------------------------------------------------------------------------

def _calc_g2(green_kw, beta):
    """Calculate G2 = T * sum_n G(k, wn) G(-k+q, -wn).

    The temperature factor T = 1/beta is included so that the
    Eliashberg kernel correctly computes:
        lambda * sigma(k) = -(T/N_L) sum_{k',n'} P(k-k') G(k') G(-k') sigma(k')

    Parameters
    ----------
    green_kw : ndarray
        Green's function, shape (norb, norb, Nx, Ny, Nz, nmat).
    beta : float
        Inverse temperature.

    Returns
    -------
    G2 : ndarray
        Shape (norb, norb, norb, norb, Nx, Ny, Nz).

    Notes
    -----
    Thin alias for ``hwave.solver.bond_channels._g2_from_green`` (review fix:
    this used to be a verbatim ~20-line copy of that function; the two are
    bit-for-bit identical by construction now that only one implementation
    exists). ``sc.py`` already imports ``bond_channels``, and ``bond_channels``
    never imports ``sc``, so there is no circular-import obstacle to sharing
    it here -- kept as ``sc._calc_g2`` only for the existing public call
    sites/tests (``tests/test_sc_bond.py`` calls ``sc._calc_g2`` directly).
    """
    return bond_channels._g2_from_green(green_kw, beta)


def _eliashberg_kernel_fft(V_q, G2, sigma_old, norb):
    """Apply one iteration of the Eliashberg kernel using FFT convolution.

    Supports both simple (2-index) and general (4-index) pairing vertices.

    For the simple case (V_q shape: norb, norb, Nx, Ny, Nz):
        sigma_{il}(k) = -sum_q V_{ij}(q) * G2_{ijlm}(q) * sigma_{jm}(k)

    For the general case (V_q shape: norb, norb, norb, norb, Nx, Ny, Nz):
        sigma_{l1l4}(k) = -sum_q sum_{l2l3} V_{l1l2,l3l4}(q) * F_{l2l3}(q)
        where F_{l2l3}(q) = sum_{l5l6} G_{l2l5}(k-q) sigma_{l5l6}(k-q) G_{l3l6}(q-k)

    Parameters
    ----------
    V_q : ndarray
        Pairing vertex. Either shape (norb, norb, Nx, Ny, Nz) for simple mode,
        or shape (norb, norb, norb, norb, Nx, Ny, Nz) for general mode.
    G2 : ndarray
        Two-particle Green's function, shape (norb, norb, norb, norb, Nx, Ny, Nz).
    sigma_old : ndarray
        Previous gap function, shape (norb, norb, Nx, Ny, Nz).
    norb : int
        Number of orbitals.

    Returns
    -------
    sigma_new : ndarray
        Updated gap function, shape (norb, norb, Nx, Ny, Nz).
    """
    Nx, Ny, Nz = sigma_old.shape[-3], sigma_old.shape[-2], sigma_old.shape[-1]
    nvol = Nx * Ny * Nz

    # G2Sigma contraction: "ijlmpqs, jmpqs -> ilpqs"
    # = matrix-vector product treating (j,m) as contracted index
    # G2: (norb,norb,norb,norb,Nx,Ny,Nz) -> (norb, norb*norb, norb, nvol)
    # sigma: (norb,norb,Nx,Ny,Nz) -> (norb*norb, nvol)
    sigma_flat = sigma_old.reshape(norb * norb, nvol)  # (jm, site)
    G2_flat = G2.reshape(norb, norb * norb, norb, nvol)  # (i, jm, l, site) -- wrong
    # Need (i, l, jm, site) for matmul on jm
    # Actually einsum is: G2[i,j,l,m,s] * sigma[j,m,s] -> result[i,l,s]
    # Reshape G2 to (norb, norb, norb*norb, nvol): G2[i, l, jm, s]
    G2_r = G2.reshape(norb, norb, norb, norb, nvol)
    G2_r = G2_r.transpose(0, 2, 1, 3, 4).reshape(norb, norb, norb * norb, nvol)
    # G2_r[i, l, jm, s], sigma_flat[jm, s]
    # result[i,l,s] = sum_jm G2_r[i,l,jm,s] * sigma_flat[jm,s]
    G2Sigma = np.einsum('iljs,js->ils', G2_r, sigma_flat).reshape(norb, norb, Nx, Ny, Nz)

    if V_q.ndim == 5:
        # Simple mode: V_q is (norb, norb, Nx, Ny, Nz)
        P_r = ifftn(V_q, axes=(-3, -2, -1))
        G2Sigma_r = ifftn(G2Sigma, axes=(-3, -2, -1))

        Sigma_r = P_r * G2Sigma_r
        sigma_new = fftn(Sigma_r, axes=(-3, -2, -1))
        return -sigma_new

    else:
        # General mode: V_q is (norb, norb, norb, norb, Nx, Ny, Nz)
        F_q = G2Sigma  # (norb, norb, Nx, Ny, Nz) = (l2, l3, ...)

        V_r = ifftn(V_q, axes=(-3, -2, -1))
        F_r = ifftn(F_q, axes=(-3, -2, -1))

        # sigma[i,l,s] = sum_{j,k} V_r[i,j,k,l,s] * F_r[j,k,s]
        # Reshape for matmul: V_r -> (norb, norb^2, norb, nvol), F_r -> (norb^2, nvol)
        V_r_flat = V_r.reshape(norb, norb * norb, norb, nvol)
        F_r_flat = F_r.reshape(norb * norb, nvol)
        sigma_r = np.einsum('ijls,js->ils', V_r_flat, F_r_flat).reshape(
            norb, norb, Nx, Ny, Nz)

        sigma_new = fftn(sigma_r, axes=(-3, -2, -1))
        return -sigma_new


# ---------------------------------------------------------------------------
# Gap function initialization
# ---------------------------------------------------------------------------

def _initialize_gap(mode, norb, kx_array, ky_array, kz_array):
    """Initialize gap function sigma(k) with specified symmetry.

    Supports 2D and 3D lattice symmetries. The form factors are defined
    in terms of kx, ky, kz and work for any dimension (Nz=1 for 2D).

    Parameters
    ----------
    mode : str
        Initialization mode specifying the gap symmetry.
        Abbreviations cx=cos(kx), sx=sin(kx), etc.

        **Isotropic / s-wave:**
        - "cos"       : cos(kx+ky+kz) (original reference code)
        - "s"         : 1 (constant, k-independent)
        - "s_ext"     : cx*cy + cy*cz + cz*cx (extended s-wave, 3D)
        - "s_ext_2d"  : cx*cy (extended s-wave, 2D)

        **d-wave:**
        - "d_x2y2"    : cx - cy
        - "d_y2z2"    : cy - cz (in-plane d_{y^2-z^2}, even parity / singlet;
                        opposite-sign anti-nodes at (pi,0) and (0,pi))
        - "d_xy"      : sx*sy
        - "d_xz"      : sx*sz
        - "d_yz"      : sy*sz
        - "d_z2"      : 2*cz - cx - cy (d_{3z^2-r^2})

        **p-wave (odd parity):**
        - "p_x"       : sx
        - "p_y"       : sy
        - "p_z"       : sz

        **f-wave (odd parity, higher harmonic):**
        - "f_x"       : sx*(cx - cy)
        - "f_y"       : sy*(cy - cx)

        Note that "f_x"/"f_y" span an odd 2D representation together with
        "p_x"/"p_y" on a square lattice, so they are seeds (and the harmonic
        basis of the subspace tracking, see ``_odd_harmonic_basis``), not a
        point-group separation.

        **Other:**
        - "random"    : random (all symmetries mixed)
    norb : int
        Number of orbitals.
    kx_array, ky_array, kz_array : ndarray
        k-point arrays.

    Returns
    -------
    sigma : ndarray
        Normalized initial gap function, shape (norb, norb, Nx, Ny, Nz).
    """
    Nx, Ny, Nz = len(kx_array), len(ky_array), len(kz_array)
    kx_mesh, ky_mesh, kz_mesh = np.meshgrid(
        kx_array, ky_array, kz_array, indexing='ij'
    )
    I = np.identity(norb)
    cx, cy, cz = np.cos(kx_mesh), np.cos(ky_mesh), np.cos(kz_mesh)
    sx, sy, sz = np.sin(kx_mesh), np.sin(ky_mesh), np.sin(kz_mesh)

    form_factors = {
        # s-wave
        "cos":      lambda: np.cos(kx_mesh + ky_mesh + kz_mesh),
        "s":        lambda: np.ones((Nx, Ny, Nz)),
        "s_ext":    lambda: cx * cy + cy * cz + cz * cx,
        "s_ext_2d": lambda: cx * cy,
        # d-wave
        "d_x2y2":   lambda: cx - cy,
        "d_y2z2":   lambda: cy - cz,
        "d_xy":     lambda: sx * sy,
        "d_xz":     lambda: sx * sz,
        "d_yz":     lambda: sy * sz,
        "d_z2":     lambda: 2.0 * cz - cx - cy,
        # p-wave
        "p_x":      lambda: sx,
        "p_y":      lambda: sy,
        "p_z":      lambda: sz,
        # f-wave (odd; the higher odd harmonic of spec S7.7)
        "f_x":      lambda: sx * (cx - cy),
        "f_y":      lambda: sy * (cy - cx),
    }

    if mode in form_factors:
        f_k = form_factors[mode]()
        sigma = I[:, :, None, None, None] * f_k[None, None, :, :, :]
    elif mode == "random":
        rng = np.random.default_rng(12345)
        sigma = rng.standard_normal((norb, norb, Nx, Ny, Nz))
    else:
        raise ValueError(
            "Unknown init_gap mode: '{}'. Available: {}".format(
                mode, list(form_factors.keys()) + ["random"]))

    norm = np.linalg.norm(sigma)
    if norm > 0:
        sigma /= norm
    elif mode != "random":
        logger.warning(
            "init_gap='%s' evaluates to zero on the current k-grid: a form "
            "factor built from sin(k_a) vanishes when the a-axis is squashed "
            "(e.g. p_x/d_xy at Nx=1, or p_z/d_xz at Nz=1). The seed is all-"
            "zero -- choose a symmetry that is non-vanishing on this grid "
            "(e.g. p_y/p_z, d_yz or d_y2z2 for an in-plane [1, Ny, Nz] cell).",
            mode)
    return sigma


def _resolve_init_gap(init_gap, pairing_type):
    """Resolve the initial-gap symmetry, defaulting to the channel's parity.

    The Eliashberg kernel preserves parity (it commutes with k -> -k), so an
    even seed cannot reach the odd triplet solution and vice versa. When the
    user does not specify ``init_gap``, default to an even s-wave ('cos') seed
    for singlet and an odd p-wave ('p_x') seed for triplet.

    Parameters
    ----------
    init_gap : str or None
        User-specified gap symmetry, or None to use the default.
    pairing_type : str
        "singlet" or "triplet".

    Returns
    -------
    str
        The gap symmetry mode to pass to _initialize_gap.
    """
    if init_gap is not None:
        return init_gap
    return "p_x" if pairing_type == "triplet" else "cos"


# The odd (triplet-parity) harmonics the invariant-subspace tracker reports the
# tracked state on (spec S7.7): {sin kx, sin ky, sin kx (cos kx - cos ky),
# sin ky (cos ky - cos kx)}. They are exactly the odd `_initialize_gap` form
# factors, so the seed table stays the single source of truth.
_ODD_HARMONICS = ("p_x", "p_y", "f_x", "f_y")


def _odd_harmonic_basis(norb, kx_array, ky_array, kz_array):
    """Normalized odd harmonic basis for the subspace decomposition (S7.7).

    Returns ``{name: flat gap vector}`` for the odd form factors of
    ``_initialize_gap`` (``p_x``, ``p_y``, ``f_x``, ``f_y``), each raveled in
    the gap layout ``(norb, norb, Nx, Ny, Nz)``. Harmonics that vanish
    identically on the current k-grid (e.g. ``p_y`` when ``Ny == 1``) are
    omitted rather than returned as a zero vector, so a decomposition never
    divides by zero.
    """
    basis = {}
    for name in _ODD_HARMONICS:
        sigma = _initialize_gap(name, norb, kx_array, ky_array, kz_array)
        if np.linalg.norm(sigma) == 0.0:
            logger.debug(
                "odd harmonic '%s' vanishes on this k-grid; omitted from the "
                "subspace decomposition basis", name)
            continue
        basis[name] = sigma.ravel()
    return basis


def _bond_diagnostics_record(provenance, gap, eigenvalues, weight, norb,
                             kx_array, ky_array, kz_array,
                             deg_tol=_BOND_DIAGNOSTICS_DEG_TOL):
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

    if eigenvalues is not None and len(eigenvalues) > 0:
        clusters = bond_channels.cluster_eigenvalues(np.asarray(eigenvalues),
                                                     deg_tol=deg_tol)
        # The cluster of the RETURNED eigenpair, i.e. of row 0 -- NOT
        # clusters[0]. cluster_eigenvalues orders by descending Re lambda over
        # every computed eigenpair, while row 0 is the eigenpair the run
        # actually reports, which _reorder_eigenpairs_by_parity promoted out of
        # the channel's parity sector; an opposite-parity mode with a larger
        # Re lambda is demoted below it but still clustered here. Reporting
        # clusters[0]'s size as "the leading state's degeneracy" would then
        # describe a different state than the one the file returns.
        leading = next(c for c in clusters if 0 in c)
        provenance["bond_diagnostics_eigenvalue_clusters"] = str(clusters)
        provenance["bond_diagnostics_leading_cluster"] = str(leading)
        provenance["bond_diagnostics_leading_cluster_size"] = len(leading)
        provenance["bond_diagnostics_deg_tol"] = "{:.3e}".format(deg_tol)
        provenance["bond_diagnostics_clusters_note"] = (
            "near-degenerate groups of the COMPUTED eigenvalues (indices into "
            "the rows below, single-linkage on |lambda_i - lambda_j| <= "
            "deg_tol, ordered by descending Re lambda over ALL computed "
            "eigenpairs -- including any opposite-parity mode that the "
            "channel-parity promotion demoted below row 0). "
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

def _solve_iteration(green_kw, Vs_q, G2, sigma_init, norb,
                     max_iter=1000, alpha=0.5, tol=1.0e-5, pairing_type=None,
                     operator=None, spectral_shift=None, info_out=None):
    """Solve linearized Eliashberg equation by self-consistent iteration.

    Parameters
    ----------
    green_kw : ndarray
        Green's function, shape (norb, norb, Nx, Ny, Nz, nmat).
    Vs_q : ndarray
        Effective pairing interaction. Either shape (norb, norb, Nx, Ny, Nz)
        for simple mode or (norb, norb, norb, norb, Nx, Ny, Nz) for general.
    G2 : ndarray
        Two-particle Green's function.
    sigma_init : ndarray
        Initial gap function.
    norb : int
        Number of orbitals.
    max_iter : int
        Maximum iterations.
    alpha : float
        Mixing parameter (0 < alpha < 1).
    tol : float
        Convergence tolerance.
    pairing_type : str, optional
        "singlet" or "triplet". When given, every iterate is projected onto the
        channel's parity sector (even for singlet, odd for triplet). For a
        centrosymmetric system the kernel commutes with parity, so in exact
        arithmetic the iteration would stay in the seed's sector; the projection
        removes the wrong-parity component that numerical noise (FFT/GEMM
        round-off) reintroduces each step and that the power iteration would
        otherwise amplify toward the dominant eigenpair of the *other* channel.
        The projection is only legitimate when the kernel commutes with parity:
        a cross-sector leakage test (several random probes) is run up front and,
        for a non-centrosymmetric kernel (complex/SOC hopping, vertex odd in k,
        chiral gaps), projection is automatically disabled with a warning and the
        un-projected iteration is used. When None, no projection is applied (the
        historical behavior).
    operator : tuple, optional
        Pre-built ``(A, vec_size)`` kernel to iterate with, bypassing
        ``_make_kernel_operator(Vs_q, G2, ...)`` (``Vs_q``/``G2`` may then be
        None). Used by the bond-channel path, whose operator is built by
        ``bond_channels.make_bond_kernel``. When None (the default) the
        historical operator is built exactly as before.
    spectral_shift : float or "auto", optional
        Positive spectral shift sigma; the power iteration then runs on
        ``K + sigma*I`` and the returned eigenvalue is shifted back (see
        ``_solve_leading``). Needed whenever the kernel is repulsive-dominant
        (negative dominant eigenvalue), where the unshifted iterate flips sign
        every step and can never converge. When None (the default) the
        historical unshifted iteration runs unchanged.
    info_out : dict, optional
        When given, updated in place with ``_solve_leading``'s info dict --
        including, on a shifted run, the Rayleigh validation record
        (``eigenvalue_validated`` / ``rayleigh_residual`` / ...). The return
        arity is deliberately unchanged (many call sites unpack exactly four
        values), so this is the way a caller learns WHAT the returned number
        is.

    Returns
    -------
    sigma : ndarray
        Converged gap function.
    eigenvalue : float
        Leading eigenvalue: without ``spectral_shift`` the (unsigned) norm of
        sigma_new before normalization; with one, the signed eigenvalue when
        the Rayleigh check validates it, else the raw shifted iterate-norm
        estimate (see ``_solve_leading`` / ``info_out``).
    converged : bool
        Whether convergence was achieved.
    n_iter : int
        Number of iterations performed.
    """
    # Precompute the Eliashberg kernel operator once. The vertex IFFT
    # (V_r = ifftn(Vs_q)) and the G2 reshape/transpose are invariant across
    # iterations (only sigma_old changes); _make_kernel_operator hoists them out
    # of the per-matvec closure, saving one full vertex IFFT (and the G2
    # preprocessing) on every one of the up-to-max_iter iterations. The matvec
    # is numerically identical to _eliashberg_kernel_fft(Vs_q, G2, sigma, norb).
    Nx, Ny, Nz = sigma_init.shape[-3], sigma_init.shape[-2], sigma_init.shape[-1]
    A, vec_size = (_make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)
                   if operator is None else operator)
    shape = (norb, norb, Nx, Ny, Nz)

    # Parity projection is only legitimate when the kernel commutes with the
    # parity operator P (Delta_{ab}(k) -> Delta_{ba}(-k), i.e.
    # _reverse_k_and_orbital), which holds for centrosymmetric systems. For a
    # non-centrosymmetric kernel (complex/SOC hopping, a vertex odd in k, chiral
    # gaps) [A, P] != 0; projecting each iterate would then discard physical
    # components and steer the power iteration to the wrong eigenpair. We
    # measure this directly as the cross-sector leakage of the kernel: apply A
    # to even- and odd-projected random probes and ask how much of the result
    # lands in the *opposite* parity sector. Zero leakage <=> [A, P] = 0 <=>
    # projection is valid. Several independent probes are used (one probe could
    # accidentally lie near the nullspace of [A, P] and miss a non-commuting
    # kernel). If any probe leaks above threshold, warn and fall back to the
    # historical un-projected iteration.
    do_project = pairing_type is not None
    if do_project and pairing_type not in ("singlet", "triplet"):
        # Validate up front so an invalid channel always raises, even if the
        # centrosymmetry guard below later disables projection (in which case
        # _project_gap_parity -- which also validates -- is never reached).
        raise ValueError(
            "Unknown pairing_type: '{}'. Use 'singlet' or 'triplet'.".format(
                pairing_type))
    if do_project:
        rng = np.random.default_rng(0)
        leakage = 0.0
        for _ in range(_PARITY_GUARD_PROBES):
            x = (rng.standard_normal(vec_size)
                 + 1j * rng.standard_normal(vec_size)).reshape(shape)
            for sign in (1.0, -1.0):  # even (+) and odd (-) projected probes
                xp = 0.5 * (x + sign * _reverse_k_and_orbital(x))
                axp = A.matvec(xp.ravel()).reshape(shape)
                # component of A xp in the OPPOSITE parity sector
                opp = 0.5 * (axp - sign * _reverse_k_and_orbital(axp))
                denom = np.linalg.norm(axp) + 1.0e-300
                leakage = max(leakage, np.linalg.norm(opp) / denom)
        if leakage > 1.0e-8:
            logger.warning(
                "Eliashberg kernel does not commute with parity (cross-sector "
                "leakage = {:.2e}); the system appears non-centrosymmetric, so "
                "parity projection for the '{}' channel is disabled and the "
                "un-projected iteration is used.".format(leakage, pairing_type))
            do_project = False

    def _project(sigma):
        return (_project_gap_parity(sigma, pairing_type)
                if do_project else sigma)

    sigma_old = _project(sigma_init.copy())
    if do_project and np.linalg.norm(sigma_old) < 1.0e-12 * (
            np.linalg.norm(sigma_init) + 1.0e-300):
        raise ValueError(
            "Initial gap has no component in the '{}' parity sector; "
            "choose an init_gap of the matching parity (even for singlet, "
            "odd for triplet).".format(pairing_type))

    # The power-iteration loop itself (matvec, projection, normalize,
    # convergence check, mixing) is factored into _solve_leading so a future
    # dynamic-frequency kernel can reuse the identical driver. `A` was already
    # built above (for the leakage probe), so `make_operator` just hands back
    # that same operator rather than rebuilding it. The gap-parity projector
    # operates on the (norb, norb, Nx, Ny, Nz) tensor, so it is wrapped as a
    # flat-vector -> flat-vector `project_fn` for the generic driver; when no
    # projection applies (`do_project` False) `project_fn` is None, matching
    # the historical (unprojected) iteration exactly.
    make_operator = lambda: (A, vec_size)
    project_fn = ((lambda flat: _project(flat.reshape(shape)).ravel())
                 if do_project else None)

    eigenvalue, sigma_flat, info = _solve_leading(
        make_operator, vec_size, "iteration",
        max_iter=max_iter, convergence_tol=tol, init_vec=sigma_old.ravel(),
        alpha=alpha, project_fn=project_fn, spectral_shift=spectral_shift,
    )

    if info_out is not None:
        info_out.update(info)
    return sigma_flat.reshape(shape), eigenvalue, info["converged"], info["n_iter"]


def _make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz):
    """Create LinearOperator for the Eliashberg kernel.

    Parameters
    ----------
    Vs_q : ndarray
        Effective pairing interaction. Either shape (norb, norb, Nx, Ny, Nz)
        for simple mode or (norb, norb, norb, norb, Nx, Ny, Nz) for general.
    G2 : ndarray
        Two-particle Green's function.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        System dimensions.

    Returns
    -------
    A : LinearOperator
        Linear operator representing the Eliashberg kernel K.
    vec_size : int
        Size of the flattened vector.
    """
    vec_size = norb * norb * Nx * Ny * Nz
    nvol = Nx * Ny * Nz

    # Backend of the resident invariants (numpy on CPU, cupy when calc_eliashberg
    # parked them on the device). Everything the kernel builds/derives must live
    # on this same backend; assert Vs_q and G2 agree.
    xp = backend.array_module_of(Vs_q)
    if backend.array_module_of(G2) is not xp:
        raise ValueError(
            "_make_kernel_operator: Vs_q and G2 must be on the same backend "
            "(both host or both device).")

    # Precompute invariants that don't change per matvec call:
    # 1. V_r = IFFT(Vs_q) -- used every matvec
    V_r = backend.spatial_ifftn(Vs_q, axes=(-3, -2, -1))

    # 2. Precompute G2 reshaped for contraction
    G2_r = G2.reshape(norb, norb, norb, norb, nvol)
    G2_pre = G2_r.transpose(0, 2, 1, 3, 4).reshape(norb, norb, norb * norb, nvol)

    is_simple = (Vs_q.ndim == 5)

    # General (4-index) mode: the two contractions
    #   G2Sigma[i,l,s] = sum_j G2_pre[i,l,j,s] * sigma_flat[j,s]
    #   sigma_r[i,l,s] = sum_j V_r_flat[i,j,l,s] * F_r_flat[j,s]
    # are, per spatial point s, a small dense matvec over j (size norb*norb).
    # numpy.einsum does NOT lower these to batched BLAS GEMM, so we precompute
    # the operator tensors ONCE in GEMM-friendly (nvol, M=il, K=j) layout and
    # apply them with matmul (s, and the column axis z for matmat, as leading
    # batch axes). This is NOT bit-identical to the einsum (GEMM changes the
    # reduction order) but matches it to ~1e-13.
    if not is_simple:
        # G2_pre is (i, l, j, s); we need (s, (i,l), j).
        G2_gemm = G2_pre.transpose(3, 0, 1, 2).reshape(
            nvol, norb * norb, norb * norb)
        # V_r_flat is (i, j, l, s); the GEMM M-axis (i,l) needs i (axis0) and
        # l (axis2) adjacent, so transpose to (s, i, l, j) -> (s, (i,l), j).
        V_r_flat = V_r.reshape(norb, norb * norb, norb, nvol)
        V_r_gemm = V_r_flat.transpose(3, 0, 2, 1).reshape(
            nvol, norb * norb, norb * norb)

    def matvec(v):
        # scipy hands us a HOST numpy vector; move it onto the invariants'
        # backend before any device op (mixed numpy/cupy operands raise on real
        # CuPy). On the numpy backend this is a no-op.
        v = xp.asarray(v)
        sigma = v.reshape(norb, norb, Nx, Ny, Nz)
        sigma_flat = sigma.reshape(norb * norb, nvol)

        if is_simple:
            G2Sigma = xp.einsum('iljs,js->ils', G2_pre, sigma_flat).reshape(
                norb, norb, Nx, Ny, Nz)
            G2Sigma_r = backend.spatial_ifftn(G2Sigma, axes=(-3, -2, -1))
            Sigma_r = V_r * G2Sigma_r
            sigma_new = backend.spatial_fftn(Sigma_r, axes=(-3, -2, -1))
        else:
            sig = xp.moveaxis(sigma_flat, 1, 0)[..., np.newaxis]
            G2Sigma = (G2_gemm @ sig)[..., 0]            # (s, (i,l))
            G2Sigma = G2Sigma.reshape(nvol, norb, norb).transpose(
                1, 2, 0).reshape(norb, norb, Nx, Ny, Nz)

            F_r = backend.spatial_ifftn(G2Sigma, axes=(-3, -2, -1))
            F_r_flat = F_r.reshape(norb * norb, nvol)
            f = xp.moveaxis(F_r_flat, 1, 0)[..., np.newaxis]
            sigma_r = (V_r_gemm @ f)[..., 0]             # (s, (i,l))
            sigma_r = sigma_r.reshape(nvol, norb, norb).transpose(
                1, 2, 0).reshape(norb, norb, Nx, Ny, Nz)
            sigma_new = backend.spatial_fftn(sigma_r, axes=(-3, -2, -1))

        # Keep complex: for complex hopping (e.g. spin-orbit coupling),
        # non-centrosymmetric models, or chiral gaps the kernel is genuinely
        # complex; projecting to real would discard physical components. Return
        # a HOST array so scipy's ARPACK/power/subspace driver stays on numpy.
        return backend.to_host(-sigma_new).ravel()

    def matmat(B):
        # Batched form of matvec: apply the SAME kernel to every column of the
        # (vec_size, k) block B in a single FFT. The batch/column axis is the
        # LAST axis so the spatial FFT axes keep absolute positions (2, 3, 4)
        # and V_r / G2_pre (no column axis) broadcast over it.
        B = xp.asarray(B)
        if B.ndim == 1:
            return matvec(B)
        k = B.shape[1]
        sigma = B.reshape(norb, norb, Nx, Ny, Nz, k)
        sigma_flat = sigma.reshape(norb * norb, nvol, k)

        if is_simple:
            G2Sigma = xp.einsum('iljs,jsz->ilsz', G2_pre, sigma_flat).reshape(
                norb, norb, Nx, Ny, Nz, k)
            G2Sigma_r = backend.spatial_ifftn(G2Sigma, axes=(2, 3, 4))
            Sigma_r = V_r[..., np.newaxis] * G2Sigma_r
            sigma_new = backend.spatial_fftn(Sigma_r, axes=(2, 3, 4))
        else:
            sig = xp.moveaxis(sigma_flat, 1, 0)          # (nvol, norb*norb, k)
            G2Sigma = G2_gemm @ sig                      # (s, (i,l), k)
            G2Sigma = G2Sigma.reshape(nvol, norb, norb, k).transpose(
                1, 2, 0, 3).reshape(norb, norb, Nx, Ny, Nz, k)

            F_r = backend.spatial_ifftn(G2Sigma, axes=(2, 3, 4))
            F_r_flat = F_r.reshape(norb * norb, nvol, k)
            f = xp.moveaxis(F_r_flat, 1, 0)              # (nvol, norb*norb, k)
            sigma_r = V_r_gemm @ f                       # (s, (i,l), k)
            sigma_r = sigma_r.reshape(nvol, norb, norb, k).transpose(
                1, 2, 0, 3).reshape(norb, norb, Nx, Ny, Nz, k)
            sigma_new = backend.spatial_fftn(sigma_r, axes=(2, 3, 4))

        return backend.to_host(-sigma_new).reshape(vec_size, k)

    A = LinearOperator((vec_size, vec_size), matvec=matvec, matmat=matmat,
                       dtype=complex)
    return A, vec_size


def _order_by_seed_overlap(vals, vecs, seed_vec):
    """Order eigenpairs by descending overlap ``|<seed, vec>|`` with a seed
    eigenvector (columns already L2-normalized by ARPACK).

    Used for eigenvector continuation: when a converged eigenvector from a
    neighbouring parameter (e.g. the next temperature) is supplied as a seed,
    the physical branch is the eigenpair whose eigenvector maximally overlaps
    it -- NOT the algebraically largest one (which, near an exceptional point
    of the non-Hermitian kernel, can jump to a different branch). Ties and a
    zero seed fall back to real-part ordering.
    """
    s = np.asarray(seed_vec).ravel()
    ns = np.linalg.norm(s)
    if ns == 0:
        return _order_eigenpairs(vals, vecs)
    s = s / ns
    ov = np.abs(vecs.conj().T @ s) / (
        np.linalg.norm(vecs, axis=0) + 1.0e-300)
    # Primary key: descending overlap. Secondary key: descending real part, so
    # an exact overlap TIE deterministically falls back to the physical
    # real-part ordering (as the docstring promises) instead of the arbitrary
    # order np.linalg.eig / ARPACK happened to return. np.lexsort takes the
    # primary key last.
    idx = np.lexsort((-vals.real, -ov))
    return vals[idx], vecs[:, idx]


def _order_eigenpairs(vals, vecs):
    """Order eigenpairs by descending real part (largest first).

    The superconducting eigenvalue is the algebraically largest one
    (lambda -> 1 at Tc); only positive eigenvalues are physically relevant.
    The leading eigenpair must therefore be the largest by real part, NOT the
    largest by magnitude: a large-magnitude *negative* eigenvalue is an
    unphysical repulsive mode and must not be reported first (nor mask a
    smaller positive eigenvalue).

    Note: ARPACK still finds eigenvalues by magnitude (which='LM'), so a small
    positive eigenvalue masked by much larger negative ones may not be among
    the returned set; request more eigenvalues, or use a shift-invert method
    targeting a positive shift, to resolve such cases.

    Parameters
    ----------
    vals : ndarray
        Eigenvalues.
    vecs : ndarray
        Eigenvectors as columns, shape (vec_size, n).

    Returns
    -------
    vals, vecs ordered so that vals[0] has the largest real part.
    """
    idx = np.argsort(-vals.real)
    return vals[idx], vecs[:, idx]


def _reverse_k_and_orbital(gap):
    """Return the gap evaluated at (-k) with orbital indices transposed.

    For the pairing gap Delta_{ab}(k), the singlet/triplet parity condition is
    Delta_{ab}(k) = +/- Delta_{ba}(-k). This helper builds Delta_{ba}(-k):
    it reverses each spatial axis (k -> -k on the FFT grid, index i -> (N-i)%N)
    and swaps the two orbital indices.

    Parameters
    ----------
    gap : ndarray
        Gap function, shape (norb, norb, Nx, Ny, Nz).

    Returns
    -------
    ndarray
        Delta_{ba}(-k), same shape.
    """
    rev = gap[:, :, ::-1, ::-1, ::-1]
    rev = np.roll(rev, 1, axis=(2, 3, 4))   # index i -> (N - i) % N
    rev = np.swapaxes(rev, 0, 1)            # orbital transpose a <-> b
    return rev


def _project_gap_parity(gap, pairing_type):
    """Project a gap onto the parity sector of the pairing channel.

    Singlet gaps are even (Delta_{ab}(k) = +Delta_{ba}(-k)); triplet gaps are
    odd (Delta_{ab}(k) = -Delta_{ba}(-k)). With P the parity operator
    (P: Delta_{ab}(k) -> Delta_{ba}(-k), built by _reverse_k_and_orbital), the
    even/odd projectors are (1 +/- P)/2. The Eliashberg kernel commutes with P,
    so these projectors commute with it and the power iteration can be confined
    to the physical sector -- removing wrong-parity components that numerical
    noise would otherwise let the iteration amplify (see _solve_iteration).

    Parameters
    ----------
    gap : ndarray
        Gap function, shape (norb, norb, Nx, Ny, Nz).
    pairing_type : str
        "singlet" (even projector) or "triplet" (odd projector).

    Returns
    -------
    ndarray
        The component of ``gap`` in the channel's parity sector, same shape.
    """
    if pairing_type == "singlet":
        sign = 1.0
    elif pairing_type == "triplet":
        sign = -1.0
    else:
        raise ValueError(
            "Unknown pairing_type: '{}'. Use 'singlet' or 'triplet'.".format(
                pairing_type))
    return 0.5 * (gap + sign * _reverse_k_and_orbital(gap))


def _is_gap_parity(gap, pairing_type, tol=0.9):
    """Test whether a gap has the parity of the pairing channel.

    Singlet gaps are even (Delta_{ab}(k) = +Delta_{ba}(-k)); triplet gaps are
    odd (Delta_{ab}(k) = -Delta_{ba}(-k)). Returns True when the projection of
    the gap onto the requested parity sector retains at least ``tol`` of its
    norm.

    Parameters
    ----------
    gap : ndarray
        Gap function, shape (norb, norb, Nx, Ny, Nz).
    pairing_type : str
        "singlet" or "triplet".
    tol : float
        Minimum fraction of the norm that must survive the parity projection.

    Returns
    -------
    bool
    """
    proj = _project_gap_parity(gap, pairing_type)
    n = np.linalg.norm(gap)
    if n == 0:
        return False
    return np.linalg.norm(proj) / n >= tol


def _reorder_eigenpairs_by_parity(vals, gaps, pairing_type):
    """Stable-reorder eigenpairs so those matching the channel parity come first.

    The Eliashberg kernel preserves parity, so its eigenvectors split into even
    (singlet) and odd (triplet) sectors. When solving for a given channel, the
    physical solution is the leading eigenpair of the matching parity, which is
    not necessarily the globally leading one. This keeps every eigenpair (so the
    returned count is unchanged) but promotes the requested-parity ones,
    preserving their existing order.

    Parameters
    ----------
    vals : ndarray
        Eigenvalues, shape (n,).
    gaps : ndarray
        Eigenvectors as gaps, shape (n, norb, norb, Nx, Ny, Nz).
    pairing_type : str
        "singlet" or "triplet".

    Returns
    -------
    vals, gaps reordered with matching-parity eigenpairs first.
    """
    match = np.array([_is_gap_parity(g, pairing_type) for g in gaps])
    if not np.any(match):
        logger.warning(
            "No computed eigenpair matches the requested '%s' parity; the "
            "reported leading gap belongs to the other channel. Increase "
            "num_eigenvalues or check the pairing_type.", pairing_type)
    idx = np.concatenate([np.where(match)[0], np.where(~match)[0]])
    return vals[idx], gaps[idx]


def _shift_from_eigenvalues(vals, factor=0.9):
    """Estimate a shift-invert target near the physical (largest-real) eigenvalue.

    The superconducting eigenvalue is the algebraically largest one, so when a
    positive eigenvalue is present among the sampled values, aim the shift at
    the largest real part (NOT the largest magnitude, which could be a large
    negative repulsive mode that masks it). If no positive eigenvalue was
    sampled (no SC instability in range), fall back to the largest-magnitude
    eigenvalue so shift-invert still tracks the dominant mode.

    Parameters
    ----------
    vals : ndarray
        Sampled eigenvalues (e.g. a few largest-magnitude ARPACK values).
    factor : float
        Fraction of the target eigenvalue to use as the shift.

    Returns
    -------
    float
    """
    real = vals.real
    scale = float(np.max(np.abs(vals)))
    # "positive" means significantly positive relative to the spectrum scale,
    # so numerical-noise eigenvalues (~1e-16) do not pull the shift to zero.
    if scale > 0 and np.any(real > 1e-8 * scale):
        return float(np.max(real)) * factor
    # No significant positive eigenvalue: track the dominant magnitude. Note a
    # purely imaginary / all-zero sample yields a 0.0 shift, which sits at a
    # possible zero mode; physical FLEX inputs have a real dominant eigenvalue.
    return float(vals[np.argmax(np.abs(vals))].real) * factor


def _validate_spectral_shift(spectral_shift, solver_mode):
    """Validate the ``spectral_shift`` option and its mode compatibility.

    ``spectral_shift`` must be a positive finite number or the string
    ``"auto"``, and is only meaningful for the two drivers that can use it:
    the ARPACK largest-real path (``eigenvalue_method='arnoldi'``) and the
    power iteration (``solver_mode='iteration'``). Anything else raises.
    """
    if spectral_shift is None:
        return
    if isinstance(spectral_shift, str):
        if spectral_shift != "auto":
            raise ValueError(
                "[eliashberg] spectral_shift string must be \"auto\", got "
                "{!r}".format(spectral_shift))
    elif isinstance(spectral_shift, bool):
        # bool is a subclass of int, and float(True) == 1.0 succeeds -- so
        # without this explicit check a TOML `spectral_shift = true` would be
        # silently accepted as a numeric shift of 1.0 instead of being
        # rejected as the type error it is (review fix).
        raise ValueError(
            "[eliashberg] spectral_shift must be a positive finite number "
            "or \"auto\", got a bool: {!r}".format(spectral_shift))
    else:
        try:
            sv = float(spectral_shift)
        except (TypeError, ValueError, OverflowError):
            sv = float("nan")
        if not np.isfinite(sv) or sv <= 0.0:
            raise ValueError(
                "[eliashberg] spectral_shift must be a positive finite "
                "number or \"auto\", got {!r}".format(spectral_shift))
    if solver_mode not in ("arnoldi", "iteration"):
        raise ValueError(
            "[eliashberg] spectral_shift is only supported for "
            "eigenvalue_method='arnoldi' and solver_mode='iteration', not "
            "{!r}".format(solver_mode))


def _estimate_spectral_radius(A, probes, project_fn=None,
                              n_steps=_AUTO_SHIFT_PROBE_STEPS):
    """Estimate the spectral radius rho(A) with a few power-iteration steps.

    For a unit vector ``v``, ``||A v||`` is a lower bound on rho(A) that
    increases toward it as ``v`` is driven onto the dominant eigenvector; the
    largest value seen over all probes and steps is returned.

    ``project_fn`` (the gap-parity projector, when the caller iterates in a
    parity sector) is applied to every probe iterate, so the estimate is the
    spectral radius OF THAT SECTOR -- the only part of the spectrum the shifted
    iteration will actually see. Several probes are used because a physically
    motivated seed can be orthogonal to the dominant eigenvector (an f-wave
    seed and an s-wave dominant mode, say) and would then wildly underestimate
    rho on its own.
    """
    rho = 0.0
    for v0 in probes:
        v = np.asarray(v0, dtype=complex).ravel()
        if project_fn is not None:
            v = project_fn(v)
        nrm = np.linalg.norm(v)
        if nrm == 0.0 or not np.isfinite(nrm):
            continue
        v = v / nrm
        for _ in range(n_steps):
            w = A.matvec(v)
            if project_fn is not None:
                w = project_fn(w)
            nw = float(np.linalg.norm(w))
            if nw == 0.0 or not np.isfinite(nw):
                break
            rho = max(rho, nw)
            v = w / nw
    return rho


def _resolve_iteration_shift(A, spectral_shift, vec_size, init_vec,
                             project_fn=None):
    """Resolve the power-iteration spectral shift sigma > 0.

    A number is taken as given. ``"auto"`` estimates the spectral radius with
    ``_estimate_spectral_radius`` (probing from the caller's seed AND from a
    deterministic random vector) and returns
    ``_AUTO_SHIFT_MARGIN * rho_est + 1e-6`` -- i.e. slightly ABOVE the
    estimated radius, the smallest shift that pushes the whole spectrum into
    the right half-plane so that the algebraically largest eigenvalue of K
    becomes the dominant (and positive) eigenvalue of K + sigma*I.
    """
    if not isinstance(spectral_shift, str):
        return float(spectral_shift)
    # "auto" (validated by _validate_spectral_shift)
    rng = np.random.default_rng(0)
    probe = (rng.standard_normal(vec_size)
             + 1j * rng.standard_normal(vec_size))
    rho = _estimate_spectral_radius(A, (init_vec, probe), project_fn)
    sig = _AUTO_SHIFT_MARGIN * rho + 1.0e-6
    logger.info(
        "auto spectral_shift: spectral radius estimated as rho = %.6f from 2 "
        "probes x %d power-iteration steps; sigma = %.2f * rho + 1e-6 = %.6f.",
        rho, _AUTO_SHIFT_PROBE_STEPS, _AUTO_SHIFT_MARGIN, sig)
    if not np.isfinite(sig):
        raise ValueError(
            "auto spectral_shift overflowed to a non-finite value (spectral "
            "radius too large); pass an explicit shift.")
    return sig


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

# Absolute floor of the Rayleigh residual tolerance, and the factor applied to
# the caller's convergence_tol. The bound is scaled by (|lambda| + 1) so it is
# relative for a large eigenvalue and absolute for a tiny one.
_RAYLEIGH_RESIDUAL_ATOL = 1.0e-6
_RAYLEIGH_RESIDUAL_TOL_FACTOR = 100.0


def _signed_rayleigh(A, vec):
    """Signed Rayleigh quotient and relative residual of ``vec`` under ``A``.

    Returns ``(lambda, residual)`` with
    ``lambda = <v|A|v> / <v|v>`` (complex; it carries the SIGN of the
    eigenvalue, unlike a power-iterate norm) and
    ``residual = ||A v - lambda v|| / ||v||``. Costs exactly one matvec. A
    zero/non-finite vector or quotient yields ``residual = inf``, i.e. "not
    an eigenvector".
    """
    vec = np.asarray(vec).ravel()
    vnorm = float(np.linalg.norm(vec))
    if vnorm == 0.0 or not np.isfinite(vnorm):
        return complex(np.nan, np.nan), float("inf")
    Av = np.asarray(A.matvec(vec)).ravel()
    lam = complex(np.vdot(vec, Av) / (vnorm * vnorm))
    if not (np.isfinite(lam.real) and np.isfinite(lam.imag)):
        return lam, float("inf")
    return lam, float(np.linalg.norm(Av - lam * vec) / vnorm)


def _validate_shifted_eigenvalue(A_unshifted, vec, iterate_value, shift,
                                 convergence_tol):
    """Decide what a shifted power iteration is allowed to report.

    Parameters
    ----------
    A_unshifted : LinearOperator
        The ORIGINAL kernel K (not ``K + sigma*I``).
    vec : ndarray
        The vector the iteration returned.
    iterate_value : float
        ``||(K + sigma I) v|| - sigma``, the historical quantity.
    shift : float
        The sigma actually applied.
    convergence_tol : float
        The iteration's own tolerance; the residual bound is derived from it.

    Returns
    -------
    value : float
        ``lambda`` (the validated Rayleigh eigenvalue) when the residual check
        passes, otherwise ``iterate_value`` unchanged.
    record : dict
        ``eigenvalue_validated`` / ``eigenvalue_kind`` (``"eigenvalue"`` vs
        ``"shifted-iterate-norm-estimate"``) / ``rayleigh_eigenvalue`` /
        ``rayleigh_residual`` / ``rayleigh_residual_tol`` /
        ``shift_sufficient``. ``shift_sufficient`` is False when the mode the
        iteration locked onto has a NEGATIVE shifted eigenvalue -- the
        signature of a too-small sigma -- and None when nothing could be
        validated.
    """
    lam, resid = _signed_rayleigh(A_unshifted, vec)
    tol = max(_RAYLEIGH_RESIDUAL_ATOL,
              _RAYLEIGH_RESIDUAL_TOL_FACTOR * float(convergence_tol))
    scale = abs(lam)
    bound = tol * ((scale if np.isfinite(scale) else 0.0) + 1.0)
    validated = bool(np.isfinite(resid) and resid <= bound
                     and abs(lam.imag) <= bound)
    record = {
        "eigenvalue_validated": validated,
        "eigenvalue_kind": ("eigenvalue" if validated
                            else "shifted-iterate-norm-estimate"),
        "rayleigh_eigenvalue": (float(lam.real) if validated else None),
        "rayleigh_residual": resid,
        "rayleigh_residual_tol": bound,
        "shift_sufficient": None,
    }
    if not validated:
        return iterate_value, record
    value = float(lam.real)
    # A sufficient sigma puts the WHOLE spectrum in the right half-plane, so
    # the mode the iteration locks onto satisfies lambda + sigma > 0. A
    # negative implied shifted eigenvalue means the iterate was driven by
    # |lambda + sigma| with lambda + sigma < 0: sigma was too small, and the
    # mode found is the largest-MAGNITUDE one, not the algebraically largest.
    record["shift_sufficient"] = bool(value + shift > 0.0)
    return value, record


def _shifted_eigenvalue_note(solver_mode, spectral_shift, converged, n_iter,
                             iter_info, kernel_label="kernel"):
    """The ``eigenvalue.dat`` comment describing a SHIFTED iteration's number.

    Split out of ``calc_eliashberg`` so the static and dynamic entry points --
    and the tests -- share one wording, and so the "this is NOT an eigenvalue"
    case can never be silently dropped from the file.
    """
    info = iter_info or {}
    resid = info.get("rayleigh_residual")
    resid_txt = ("{:.3e}".format(resid) if isinstance(resid, float)
                 and np.isfinite(resid) else str(resid))
    tail = ("the iteration {} after {} steps."
            .format("converged" if converged else "did NOT converge", n_iter))
    if info.get("eigenvalue_validated"):
        note = (
            "solver_mode='{}' with spectral_shift={!r}: the power iteration "
            "ran on the shifted {} K + sigma*I and the value below is the "
            "SIGNED eigenvalue of K (the shift has been subtracted back), "
            "VALIDATED by the Rayleigh quotient <v|K|v>/<v|v> with residual "
            "||K v - lambda v||/||v|| = {}; {}"
            .format(solver_mode, spectral_shift, kernel_label, resid_txt,
                    tail))
        if info.get("shift_sufficient") is False:
            note += (" WARNING: spectral_shift appears INSUFFICIENT -- the "
                     "mode found has lambda + sigma < 0, so it is the "
                     "largest-MAGNITUDE eigenvalue of K + sigma*I, not the "
                     "algebraically largest (physical) one; re-run with a "
                     "larger spectral_shift or \"auto\".")
        return note
    return (
        "solver_mode='{}' with spectral_shift={!r}: the value below is a "
        "SHIFTED ITERATE-NORM ESTIMATE, ||(K + sigma*I) v|| - sigma. It is "
        "NOT an eigenvalue of K: the Rayleigh check <v|K|v>/<v|v> left a "
        "residual ||K v - lambda v||/||v|| = {} (tolerance {}), so no signed "
        "eigenvalue could be validated for this run -- do not quote it as "
        "lambda. {} Raise max_iter, raise spectral_shift, or use "
        "solver_mode=\"eigenvalue\"."
        .format(solver_mode, spectral_shift, resid_txt,
                info.get("rayleigh_residual_tol"), tail))


def _solve_leading(make_operator, vec_size, solver_mode, num_eigenvalues=10,
                   max_iter=1000, convergence_tol=1.0e-5, init_vec=None,
                   sigma_shift=None, alpha=0.5, project_fn=None, seed_vec=None,
                   spectral_shift=None):
    """Shared leading-eigenpair driver behind the static Eliashberg solvers.

    This holds the ARPACK/shift-invert eigen-selection-and-ordering body of
    ``_solve_eigenvalue`` and the power-iteration loop of ``_solve_iteration``,
    generalized to act on any ``make_operator``/``vec_size`` pair (not just the
    static ``_make_kernel_operator`` result) so a future dynamic-frequency
    kernel can reuse the exact same driver.

    Parameters
    ----------
    make_operator : callable
        Zero-argument callable returning ``(A, op_vec_size)``, i.e. the same
        convention as ``_make_kernel_operator`` (a
        ``scipy.sparse.linalg.LinearOperator`` acting on flat vectors of
        length ``vec_size``, plus that size). Called exactly once.
    vec_size : int
        Size of the flattened vector space that ``A`` acts on.
    solver_mode : str
        "arnoldi", "shift-invert-bicgstab", "shift-invert-gmres", or
        "shift-invert-lgmres" select the ARPACK/shift-invert eigenanalysis
        (the former body of ``_solve_eigenvalue``, excluding its "subspace"
        branch which stays in that wrapper). "iteration" selects the power
        loop (the former body of ``_solve_iteration``).
    num_eigenvalues : int
        Number of eigenvalues to request from ARPACK. Only used by the
        eigenvalue-family modes.
    max_iter : int
        Maximum number of power-iteration steps. Only used by "iteration".
    convergence_tol : float
        Power-iteration convergence tolerance on the normalized-iterate
        difference. Only used by "iteration".
    init_vec : ndarray, optional
        Flat initial vector, shape ``(vec_size,)``, for "iteration" mode
        (already projected onto the desired sector by the caller, if
        applicable). Required for "iteration" mode.
    sigma_shift : float, optional
        Shift-invert target for the "shift-invert-*" modes. If None, it is
        estimated from a preliminary Arnoldi pass, exactly as
        ``_solve_eigenvalue`` does.
    alpha : float
        Power-iteration mixing parameter. Only used by "iteration".
    project_fn : callable, optional
        Flat-vector -> flat-vector projector applied to every power-iteration
        iterate (e.g. the gap-parity projector in ``_solve_iteration``). Only
        used by "iteration"; when None, no projection is applied.
    spectral_shift : float or "auto", optional
        Positive spectral shift sigma. Supported by "arnoldi" (ARPACK asks for
        the largest-REAL eigenvalue of ``A + sigma*I``) and by "iteration"
        (the power loop runs on ``A + sigma*I``, whose dominant eigenvalue is
        positive even for a repulsive-dominant kernel, so it converges). In
        both cases sigma is subtracted back: the returned eigenvalue is the
        SIGNED eigenvalue of the original operator and the eigenvector is
        unchanged by the shift. ``"auto"`` estimates sigma from the spectral
        radius (see ``_resolve_iteration_shift`` for the iteration path).
        Rejected for every other mode.

    Returns
    -------
    leading_eigenvalue : complex or float
        The dominant eigenvalue (eigenvalue-family: largest real part, per
        ``_order_eigenpairs``). On the iteration path WITHOUT
        ``spectral_shift`` it is the converged/last iterate norm (unsigned),
        exactly as before. WITH a ``spectral_shift`` it is the signed Rayleigh
        eigenvalue ``<v|K|v>/<v|v>`` when that passes the residual check
        (``eig_analysis["eigenvalue_validated"]``), and otherwise the raw
        ``||(K + sigma I) v|| - sigma``, which is NOT an eigenvalue and must be
        reported as an estimate (see ``_validate_shifted_eigenvalue``).
    leading_eigenvector : ndarray
        Flat, shape ``(vec_size,)``.
    eig_analysis : dict
        Eigenvalue-family modes: ``{"eigenvalues": vals, "eigenvectors": vecs,
        "sigma_shift": sigma_shift}`` with ``vals`` ordered by descending real
        part (``_order_eigenpairs``) and ``vecs`` the matching eigenvectors as
        columns. "iteration" mode: ``{"converged": bool, "n_iter": int}`` --
        exactly the historical key set when no shift is requested -- plus,
        ONLY when a shift was actually applied (so the default path's dict is
        unchanged), ``spectral_shift`` (the sigma actually applied),
        ``eigenvalue_validated``, ``eigenvalue_kind``, ``rayleigh_eigenvalue``,
        ``rayleigh_residual``, ``rayleigh_residual_tol`` and
        ``shift_sufficient``.
    """
    # Validate spectral_shift up front -- before any branch, including the
    # iteration loop and the small-dense early return -- so an invalid value or
    # an incompatible mode fails fast everywhere.
    _validate_spectral_shift(spectral_shift, solver_mode)

    if solver_mode == "iteration":
        if init_vec is None:
            raise ValueError("init_vec is required for solver_mode='iteration'")
        A, _ = make_operator()

        # Spectral shift: iterate with K + sigma*I instead of K. A repulsive-
        # dominant kernel has a NEGATIVE dominant eigenvalue, which flips the
        # iterate's sign every step, so the convergence test
        # ||sigma_new/||sigma_new|| - sigma_old|| sits at ~2 and the loop can
        # never converge. Shifting the spectrum into the right half-plane makes
        # the ALGEBRAICALLY LARGEST eigenvalue dominant and positive (the same
        # physical selection as the arnoldi which='LR' path). The shift leaves
        # every eigenvector -- and hence every invariant subspace, including
        # the parity sectors that project_fn selects -- untouched; only the
        # eigenvalue moves, and it is shifted back below.
        shift = 0.0
        if spectral_shift is not None:
            shift = _resolve_iteration_shift(A, spectral_shift, vec_size,
                                             init_vec, project_fn)
            logger.info(
                "Power iteration with spectral shift sigma = %.6f (%s): "
                "iterating on K + sigma*I; the reported eigenvalue is "
                "un-shifted (lambda = lambda_shifted - sigma) and the "
                "eigenvector is unchanged by the shift.",
                shift,
                "auto" if isinstance(spectral_shift, str) else "explicit")
            A_unshifted = A
            A = LinearOperator(
                (vec_size, vec_size),
                matvec=lambda v: A_unshifted.matvec(v) + shift * v,
                dtype=complex)
        # the sigma actually applied (None when no shift was requested), so
        # callers can report it; `shift` itself is the arithmetic value.
        shift_used = shift if spectral_shift is not None else None

        def _finish(value, vec, converged, n_iter_done, rayleigh=True):
            """Assemble the return triple, validating the SHIFTED value.

            On the unshifted path this is exactly the historical dict --
            ``{"converged", "n_iter"}`` and NOTHING else, as at commit
            712b1a0; ``spectral_shift`` is added only when a shift was really
            applied, so a caller that compares or serializes the default
            dict sees the key set it always saw -- and no extra matvec is
            taken, so the default behaviour is bit-identical.
            """
            info = {"converged": converged, "n_iter": n_iter_done}
            if shift_used is None:
                return value, vec, info
            info["spectral_shift"] = shift_used
            if not rayleigh:
                # Nothing meaningful to validate (the zero-norm sentinel).
                info.update({"eigenvalue_validated": False,
                             "eigenvalue_kind": "not-computed",
                             "rayleigh_eigenvalue": None,
                             "rayleigh_residual": float("inf"),
                             "rayleigh_residual_tol": float("nan"),
                             "shift_sufficient": None})
                return value, vec, info
            value, record = _validate_shifted_eigenvalue(
                A_unshifted, vec, value, shift, convergence_tol)
            info.update(record)
            if record["eigenvalue_validated"]:
                logger.info(
                    "Shifted power iteration: reporting the SIGNED Rayleigh "
                    "eigenvalue lambda = %.8e of K (residual %.3e <= %.3e), "
                    "not the raw shifted iterate norm.",
                    value, record["rayleigh_residual"],
                    record["rayleigh_residual_tol"])
                if record["shift_sufficient"] is False:
                    logger.warning(
                        "spectral_shift = %.6f appears INSUFFICIENT: the "
                        "iteration locked onto a mode with lambda = %.8e, so "
                        "lambda + sigma = %.8e < 0 and the power iteration "
                        "selected the largest-MAGNITUDE eigenvalue of "
                        "K + sigma*I rather than the algebraically largest "
                        "(physical) one. Re-run with a larger spectral_shift "
                        "or spectral_shift=\"auto\".",
                        shift, value, value + shift)
            else:
                logger.warning(
                    "Shifted power iteration: the returned value %.8e is "
                    "||(K + sigma*I) v|| - sigma and is NOT an eigenvalue of "
                    "K -- the Rayleigh check <v|K|v>/<v|v> left a residual of "
                    "%.3e (tolerance %.3e), so no signed eigenvalue could be "
                    "validated. Do not quote it as lambda; raise max_iter, "
                    "raise spectral_shift, or use solver_mode=\"eigenvalue\".",
                    value, record["rayleigh_residual"],
                    record["rayleigh_residual_tol"])
            return value, vec, info

        sigma_old = init_vec
        eigenvalue = 0.0

        for iteration in range(max_iter):
            sigma_new = A.matvec(sigma_old)
            if project_fn is not None:
                sigma_new = project_fn(sigma_new)
            norm = np.linalg.norm(sigma_new)

            # A (possibly projected) iterate can collapse to the zero vector
            # when the kernel annihilates the requested sector; normalizing by
            # `norm` would then produce NaN. `np.linalg.norm` squares its
            # inputs, so it returns *exactly* 0.0 once the iterate underflows;
            # guard precisely on that, plus any non-finite norm, and report
            # non-convergence with the last finite iterate instead.
            if norm == 0.0 or not np.isfinite(norm):
                logger.warning(
                    "Eliashberg iterate collapsed to zero norm at iteration "
                    "{}; the kernel annihilates the requested sector, so the "
                    "eigenvalue cannot be normalized. Returning the previous "
                    "iterate as non-converged.".format(iteration))
                # 0.0 is the historical non-converged sentinel, NOT a shifted
                # eigenvalue, so it is returned as-is (un-shifting it would
                # dress a failure up as lambda = -sigma) and no Rayleigh
                # validation is attempted on it.
                return _finish(0.0, sigma_old, False, iteration + 1,
                               rayleigh=False)

            # `norm` is the eigenvalue of the (possibly shifted) operator;
            # subtract the shift back so every reported/logged value belongs to
            # the ORIGINAL operator K. `shift` is exactly 0.0 when no spectral
            # shift is in use, so the un-shifted path is bit-identical.
            eigenvalue = norm - shift
            diff = np.linalg.norm(sigma_new / norm - sigma_old)
            logger.info("Iteration {:4d}: eigenvalue = {:.6f}, diff = {:.6e}".format(
                iteration, eigenvalue, diff))

            if diff < convergence_tol:
                logger.info("Converged at iteration {}".format(iteration + 1))
                return _finish(eigenvalue, sigma_new / norm, True,
                               iteration + 1)

            sigma_old = (1.0 - alpha) * sigma_new / norm + alpha * sigma_old

        logger.warning("Failed to converge after {} iterations".format(max_iter))
        if shift_used is not None:
            logger.warning(
                "The power iteration ran on the shifted operator K + sigma*I "
                "with sigma = %.6f and still did not converge. Either sigma is "
                "too small (the dominant eigenvalue of K + sigma*I is still "
                "negative -- try a larger spectral_shift, or \"auto\") or it is "
                "so large that it compressed the relative eigenvalue gaps "
                "(try a smaller one, or solver_mode=\"eigenvalue\").", shift)
        return _finish(eigenvalue, sigma_old, False, max_iter)

    # Eigenvalue family: ARPACK Arnoldi or shift-invert.
    A, _ = make_operator()

    if not (solver_mode == "arnoldi" or solver_mode.startswith("shift-invert")):
        raise ValueError("Unknown eigenvalue method: {}".format(solver_mode))

    # (spectral_shift was validated at the top of the function, before every
    # branch, so an invalid value fails fast even on the small-dense early
    # return below.)

    # sigma_shift is the shift-invert target; it has no effect on the plain
    # arnoldi path (which uses spectral_shift instead). Warn rather than fail so
    # existing configs keep working.
    if sigma_shift is not None and solver_mode == "arnoldi":
        logger.warning(
            "[eliashberg] sigma_shift is ignored for eigenvalue_method="
            "'arnoldi' (it targets the shift-invert methods); use "
            "spectral_shift to bias the arnoldi selection.")

    if vec_size < 1:
        raise ValueError("Eliashberg operator has empty vector space")

    if vec_size <= 2:
        # scipy.sparse.linalg.eigs requires k < N - 1 for LinearOperator input,
        # so the smallest valid dynamic grid (e.g. norb=1, Nk=1, Nmat=2) cannot
        # go through ARPACK. Reconstruct the tiny dense operator directly.
        dense = np.empty((vec_size, vec_size), dtype=complex)
        basis = np.eye(vec_size, dtype=complex)
        for j in range(vec_size):
            dense[:, j] = A.matvec(basis[:, j])
        vals, vecs = np.linalg.eig(dense)
        # Match the ARPACK path below: with a seed eigenvector, track the branch
        # that overlaps it (eigenvector continuation); otherwise order by
        # largest real part (the physical SC eigenvalue), not magnitude.
        if seed_vec is not None:
            vals, vecs = _order_by_seed_overlap(vals, vecs, seed_vec)
        else:
            vals, vecs = _order_eigenpairs(vals, vecs)
        n_keep = min(max(1, num_eigenvalues), vec_size)
        vals = vals[:n_keep]
        vecs = vecs[:, :n_keep]
        return vals[0], vecs[:, 0], {"eigenvalues": vals,
                                     "eigenvectors": vecs,
                                     "sigma_shift": sigma_shift}

    max_ev = min(num_eigenvalues, vec_size - 2)
    if max_ev < 1:
        max_ev = 1

    logger.info("Computing {} eigenvalues with method='{}'...".format(
        max_ev, solver_mode))

    if solver_mode == "arnoldi":
        if spectral_shift is not None:
            # Select the LARGEST-REAL eigenvalue (the physical SC eigenvalue,
            # lambda -> 1 at Tc), which plain which='LM' misses when a small
            # positive lambda is masked by larger repulsive (negative)
            # eigenvalues. We ask ARPACK for the largest real part directly
            # (which='LR') -- unlike which='LM', this is the correct criterion
            # even for the non-Hermitian kernel's complex eigenvalues, where
            # |lambda+sigma| would otherwise be dominated by a large-|Im| mode.
            # The real spectral shift A -> A + sigma*I preserves the real-part
            # ordering (Re(lambda+sigma) = Re(lambda) + sigma) but moves the
            # spectrum into the right half-plane, which conditions ARPACK's LR
            # iteration; we subtract sigma afterwards.
            from scipy.sparse.linalg import LinearOperator as _LinOp
            # spectral_shift is validated (positive-finite / "auto", arnoldi)
            # near the top of _solve_leading.
            if isinstance(spectral_shift, str):  # "auto"
                k_pre = min(6, vec_size - 2)
                if k_pre >= 1:
                    vals_pre, _ = eigs(A, k=k_pre, which='LM')
                    sig = float(np.max(np.abs(vals_pre))) * 1.5 + 1.0e-6
                else:
                    sig = 1.0
                if not np.isfinite(sig):
                    raise ValueError(
                        "auto spectral_shift overflowed to a non-finite value "
                        "(spectral radius too large); pass an explicit shift.")
            else:
                sig = float(spectral_shift)
            logger.info("Spectral shift sigma={:.6f}: eigs(A+sigma*I, LR)".format(sig))
            A_sh = _LinOp(A.shape, matvec=lambda v: A.matvec(v) + sig * v,
                          dtype=A.dtype)
            vals, vecs = eigs(A_sh, k=max_ev, which='LR', v0=seed_vec)
            vals = vals - sig
        else:
            vals, vecs = eigs(A, k=max_ev, which='LM', v0=seed_vec)

    elif solver_mode.startswith("shift-invert"):
        if sigma_shift is None:
            # Estimate shift from a quick Arnoldi run. Sample a few
            # largest-magnitude eigenvalues and aim at the largest *real* part
            # (the physical SC eigenvalue), not the largest magnitude (which
            # can be a large negative repulsive mode).
            k_pre = min(6, vec_size - 2)
            if k_pre < 1:
                # Operator too small for a preliminary ARPACK pass (ARPACK
                # needs k < N-1); fall back to a neutral shift.
                sigma_shift = 0.0
            else:
                logger.info("Estimating shift with preliminary Arnoldi...")
                vals_pre, _ = eigs(A, k=k_pre, which='LM')
                sigma_shift = _shift_from_eigenvalues(vals_pre)
            logger.info("Using sigma_shift = {:.6f}".format(sigma_shift))
        vals, vecs = _eigs_shift_invert(
            A, vec_size, max_ev, solver_mode, sigma=sigma_shift,
            seed_vec=seed_vec
        )

    # With a seed eigenvector, track the branch that overlaps it (eigenvector
    # continuation); otherwise order by largest real part (the physical SC
    # eigenvalue), not magnitude.
    if seed_vec is not None:
        vals, vecs = _order_by_seed_overlap(vals, vecs, seed_vec)
    else:
        vals, vecs = _order_eigenpairs(vals, vecs)

    # A negative leading eigenvalue from plain which='LM' (no spectral_shift) is
    # often an artifact: a small positive lambda masked by larger repulsive
    # modes. Tip the user toward spectral_shift='auto'. Skip this when a
    # seed_vec is given -- there vals[0] is the seed-overlap continuation
    # branch, which may be intentionally negative, not a masked leading mode.
    # Require a meaningfully negative value (relative to the spectral scale)
    # so roundoff-scale negatives near a numerically-zero leading eigenvalue
    # do not trigger a misleading recommendation.
    if (solver_mode == "arnoldi" and spectral_shift is None
            and seed_vec is None and len(vals)):
        scale = float(np.max(np.abs(vals))) if len(vals) else 0.0
        neg_tol = 1.0e-8 * max(scale, 1.0)
        if vals[0].real < -neg_tol:
            logger.warning(
                "Leading eigenvalue Re(lambda)=%.4g is negative; if a positive "
                "(attractive) mode is expected, set [eliashberg] spectral_shift="
                "\"auto\" so the largest-REAL eigenvalue is selected instead of "
                "the largest-magnitude one.", vals[0].real)

    return vals[0], vecs[:, 0], {"eigenvalues": vals, "eigenvectors": vecs,
                                 "sigma_shift": sigma_shift}


def _solve_eigenvalue(Vs_q, G2, norb, Nx, Ny, Nz, num_eigenvalues=10,
                      method="arnoldi", sigma_shift=None, spectral_shift=None,
                      operator=None):
    """Solve linearized Eliashberg equation by eigenvalue analysis.

    Parameters
    ----------
    Vs_q : ndarray
        Effective pairing interaction. Either shape (norb, norb, Nx, Ny, Nz)
        for simple mode or (norb, norb, norb, norb, Nx, Ny, Nz) for general.
    G2 : ndarray
        Two-particle Green's function.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        System dimensions.
    num_eigenvalues : int
        Number of eigenvalues to compute.
    method : str
        Eigenvalue solver method:
        - "arnoldi" : Implicitly Restarted Arnoldi (ARPACK eigs).
            Standard Krylov subspace method for non-symmetric matrices.
            Directly finds the largest eigenvalues.
        - "shift-invert-bicgstab" : Shift-invert with BiCGSTAB.
            Solves (K - sigma*I)^{-1} eigenvalue problem using BiCGSTAB
            for the linear system. Efficient for finding eigenvalues near
            a target value sigma.
        - "shift-invert-gmres" : Shift-invert with GMRES.
            Same as above but uses GMRES for the linear system.
            More robust than BiCGSTAB for non-symmetric systems.
        - "shift-invert-lgmres" : Shift-invert with LGMRES.
            Uses LGMRES (loose GMRES) which can be faster for some problems.
        - "subspace" : Subspace iteration (block power method).
            Simultaneously propagates multiple vectors through the kernel
            with QR orthogonalization. Robust for degenerate eigenvalues
            and finds multiple eigenvalues with different symmetries.
    sigma_shift : float, optional
        Shift parameter for shift-invert methods. Eigenvalues near this
        value are found most efficiently. If None, a preliminary Arnoldi
        step estimates an appropriate shift value.
    operator : tuple, optional
        Pre-built ``(A, vec_size)`` kernel, bypassing
        ``_make_kernel_operator(Vs_q, G2, ...)`` (see ``_solve_iteration``).

    Returns
    -------
    eigenvalues : ndarray
        Leading eigenvalues ordered by descending real part (largest first).
    eigenvectors : ndarray
        Corresponding eigenvectors reshaped to (num_ev, norb, norb, Nx, Ny, Nz).
    """
    vec_size = norb * norb * Nx * Ny * Nz

    # Safety check (review fix): when a pre-built ``operator`` is supplied,
    # its own vec_size (``operator[1]``) is otherwise never read here -- every
    # downstream call (``_solve_leading``, ``_solve_subspace_iteration``) is
    # handed the LOCALLY recomputed ``norb*norb*Nx*Ny*Nz`` instead. The two
    # agree today for both callers of this path (the standard
    # ``_make_kernel_operator`` result and the bond-channel
    # ``bond_channels.make_bond_kernel*`` result both use the plain
    # orbital-pair x spatial vec_size), but nothing enforces that; a future
    # dynamic/enlarged bond operator with a DIFFERENT internal size would
    # silently reshape eigenvectors into the wrong shape instead of failing
    # loudly. Assert the invariant instead of trusting it silently.
    if operator is not None:
        assert operator[1] == vec_size, (
            "_solve_eigenvalue: the supplied operator's vec_size ({}) does "
            "not match norb**2*Nx*Ny*Nz ({}); eigenvectors would be reshaped "
            "with the wrong size.".format(operator[1], vec_size))

    if method == "subspace":
        # spectral_shift routes through the ARPACK largest-real path in
        # _solve_leading; the subspace driver never reaches it, so reject a
        # non-None spectral_shift here (same arnoldi-only contract) rather
        # than silently ignoring a misconfigured static input.
        if spectral_shift is not None:
            raise ValueError(
                "[eliashberg] spectral_shift is only supported for "
                "eigenvalue_method='arnoldi', not '{}'".format(method))
        # Subspace (block power) iteration has its own dedicated driver
        # (magnitude-based Ritz selection, not the ARPACK/shift-invert path),
        # so it is not routed through _solve_leading; call it directly, exactly
        # as before.
        max_ev = min(num_eigenvalues, vec_size - 2)
        if max_ev < 1:
            max_ev = 1
        logger.info("Computing {} eigenvalues with method='{}'...".format(
            max_ev, method))
        return _solve_subspace_iteration(
            Vs_q, G2, norb, Nx, Ny, Nz,
            num_eigenvalues=max_ev, operator=operator
        )

    # ARPACK Arnoldi / shift-invert: delegate the eigen-selection, shift
    # estimation, and descending-real-part ordering to the shared driver.
    make_operator = (lambda: _make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)) \
        if operator is None else (lambda: operator)
    _, _, eig_analysis = _solve_leading(
        make_operator, vec_size, method,
        num_eigenvalues=num_eigenvalues, sigma_shift=sigma_shift,
        spectral_shift=spectral_shift,
    )
    vals = eig_analysis["eigenvalues"]
    vecs = eig_analysis["eigenvectors"]

    eigenvectors = np.array([
        vecs[:, i].reshape(norb, norb, Nx, Ny, Nz)
        for i in range(len(vals))
    ])

    return vals, eigenvectors


def _eigs_shift_invert(A, vec_size, num_ev, method, sigma=0.0, rtol_linear=1e-8,
                       seed_vec=None):
    """Eigenvalue computation using shift-invert with iterative linear solver.

    Transforms the eigenvalue problem K*x = lambda*x into
    (K - sigma*I)^{-1}*x = nu*x where nu = 1/(lambda - sigma).
    Eigenvalues of the original problem near sigma become the largest
    eigenvalues of the transformed problem, which Arnoldi finds efficiently.

    The inverse (K - sigma*I)^{-1} is never formed explicitly; instead,
    each application solves (K - sigma*I)*y = x using an iterative solver.

    Parameters
    ----------
    A : LinearOperator
        The Eliashberg kernel operator K.
    vec_size : int
        Dimension of the vector space.
    num_ev : int
        Number of eigenvalues to compute.
    method : str
        One of "shift-invert-bicgstab", "shift-invert-gmres",
        "shift-invert-lgmres".
    sigma : float
        Shift value. Eigenvalues near sigma are found most efficiently.
        Default 0.0 finds eigenvalues nearest to zero.
    rtol_linear : float
        Relative tolerance for the iterative linear solver.

    Returns
    -------
    eigenvalues : ndarray
        Eigenvalues of the original problem K.
    eigenvectors : ndarray
        Corresponding eigenvectors.
    """
    solver_name = method.replace("shift-invert-", "")
    solver_map = {
        "bicgstab": bicgstab,
        "gmres": gmres,
        "lgmres": lgmres,
    }
    if solver_name not in solver_map:
        raise ValueError("Unknown linear solver: {}".format(solver_name))
    linear_solver = solver_map[solver_name]

    # (A - sigma*I)
    def shifted_matvec(v):
        return A.matvec(v) - sigma * v

    A_shifted = LinearOperator((vec_size, vec_size),
                               matvec=shifted_matvec, dtype=complex)

    solve_count = [0]
    fail_count = [0]

    # (A - sigma*I)^{-1} via iterative solver
    def inv_matvec(v):
        solve_count[0] += 1
        x, info = linear_solver(A_shifted, v, rtol=rtol_linear, maxiter=500)
        if info != 0:
            fail_count[0] += 1
            if fail_count[0] <= 3:
                logger.warning(
                    "Linear solver {} did not converge (info={}), "
                    "solve #{}, using approximate solution".format(
                        solver_name, info, solve_count[0]))
        return x

    A_inv = LinearOperator((vec_size, vec_size),
                           matvec=inv_matvec, dtype=complex)

    # eigs on (A - sigma*I)^{-1} finds eigenvalues nu = 1/(lambda - sigma)
    # largest |nu| correspond to lambda closest to sigma. A seed vector (v0)
    # biases the Arnoldi start toward the physical branch being tracked.
    nus, vecs = eigs(A_inv, k=num_ev, which='LM', v0=seed_vec)

    logger.info("Shift-invert: {} linear solves, {} failures".format(
        solve_count[0], fail_count[0]))

    # Convert back: lambda = 1/nu + sigma
    eigenvalues = 1.0 / nus + sigma

    return eigenvalues, vecs


def _solve_subspace_iteration(Vs_q, G2, norb, Nx, Ny, Nz,
                              num_eigenvalues=5, max_iter=300, tol=1e-6,
                              operator=None):
    """Find multiple eigenvalues by subspace iteration (block power method).

    Simultaneously propagates a block of vectors through the kernel,
    orthogonalizing via QR decomposition at each step. This naturally
    converges to the dominant invariant subspace.

    Parameters
    ----------
    Vs_q : ndarray
        Effective pairing interaction. Either shape (norb, norb, Nx, Ny, Nz)
        for simple mode or (norb, norb, norb, norb, Nx, Ny, Nz) for general.
    G2 : ndarray
        Two-particle Green's function.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        System dimensions.
    num_eigenvalues : int
        Number of eigenvalues to compute.
    max_iter : int
        Maximum iterations.
    tol : float
        Convergence tolerance for eigenvalues.
    operator : tuple, optional
        Pre-built ``(A, vec_size)`` kernel, bypassing
        ``_make_kernel_operator(Vs_q, G2, ...)`` (see ``_solve_iteration``).

    Returns
    -------
    eigenvalues : ndarray
        Converged eigenvalues ordered by descending magnitude.
    eigenvectors : ndarray
        Shape (num_eigenvalues, norb, norb, Nx, Ny, Nz).
    """
    A, vec_size = (_make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)
                   if operator is None else operator)
    num_ev = min(num_eigenvalues, vec_size)

    # Use extra vectors for better convergence (subspace padding)
    n_work = min(num_ev + max(num_ev, 5), vec_size)

    # Random initial subspace
    rng = np.random.default_rng(42)
    V = rng.standard_normal((vec_size, n_work))
    V, _ = np.linalg.qr(V)

    eigenvalues_old = np.zeros(num_ev, dtype=complex)

    for iteration in range(max_iter):
        # Apply kernel to all vectors at once: W = A @ V (kernel may be
        # complex). matmat batches all columns into a single FFT (numerically
        # identical to the per-column matvec).
        W = A.matmat(V)

        # Rayleigh quotient: H = V^dagger A V (conjugate transpose for complex)
        H = V.conj().T @ W

        # Eigendecomposition of small matrix
        evals_h, evecs_h = np.linalg.eig(H)
        idx = np.argsort(-np.abs(evals_h))
        evals_h = evals_h[idx]
        evecs_h = evecs_h[:, idx]

        # Ritz vectors: update subspace
        V = W @ evecs_h
        # Re-orthogonalize
        V, _ = np.linalg.qr(V)

        # Check convergence of the wanted eigenvalues (track the full complex
        # value so imaginary motion is not ignored for a complex kernel).
        eigenvalues_new = evals_h[:num_ev]
        diff = np.max(np.abs(eigenvalues_new - eigenvalues_old))

        if iteration % 10 == 0 or diff < tol:
            logger.info("Subspace iter {:4d}: eigenvalues = {}, diff = {:.2e}".format(
                iteration, np.array2string(eigenvalues_new, precision=4), diff))

        if diff < tol:
            logger.info("Subspace iteration converged at iteration {}".format(
                iteration + 1))
            break

        eigenvalues_old = eigenvalues_new.copy()

    # Extract final Ritz vectors for the wanted eigenvalues
    # Recompute from final V (batched matmat = per-column matvec)
    W = A.matmat(V)
    H = V.conj().T @ W
    evals_h, evecs_h = np.linalg.eig(H)
    # Subspace (block power) iteration is magnitude-based: it converges to the
    # dominant-magnitude invariant subspace, so report its modes by magnitude.
    # The physical SC selection (largest real) is applied in the default
    # arnoldi path; use that to read off the leading SC eigenvalue.
    idx = np.argsort(-np.abs(evals_h))

    eigenvalues = evals_h[idx[:num_ev]]
    ritz_vecs = V @ evecs_h[:, idx[:num_ev]]

    eigenvectors = np.array([
        ritz_vecs[:, i].reshape(norb, norb, Nx, Ny, Nz)
        for i in range(num_ev)
    ])

    return eigenvalues, eigenvectors


def _solve_shifted_bicg(Vs_q, G2, norb, Nx, Ny, Nz,
                        sigma_list, num_eigenvalues=3, tol_linear=1e-8):
    """Find eigenvalues near multiple target values using shifted BiCG.

    Note: this helper does not accept an ``operator=`` (bond-channel)
    argument -- it always builds its own kernel from ``Vs_q``/``G2`` via
    ``_make_kernel_operator``. It is unreachable from ``calc_eliashberg``
    when ``bond_channels=true`` (that path sets ``Vs_q = G2 = None`` and
    drives the pre-built bond operator through ``_solve_iteration``/
    ``_solve_eigenvalue`` instead; nothing in the bond dispatch calls this
    function).

    Parameters
    ----------
    Vs_q : ndarray
        Effective pairing interaction. Either shape (norb, norb, Nx, Ny, Nz)
        for simple mode or (norb, norb, norb, norb, Nx, Ny, Nz) for general.
    G2 : ndarray
        Two-particle Green's function.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        System dimensions.
    sigma_list : list of float
        Target shift values. Eigenvalues near each sigma are found.
    num_eigenvalues : int
        Number of eigenvalues to find near each sigma.
    tol_linear : float
        Tolerance for the linear solver.

    Returns
    -------
    all_eigenvalues : dict
        Dictionary mapping sigma -> array of eigenvalues.
    all_eigenvectors : dict
        Dictionary mapping sigma -> array of eigenvectors.
    """
    A, vec_size = _make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)

    max_ev = min(num_eigenvalues, vec_size - 2)
    if max_ev < 1:
        max_ev = 1

    all_eigenvalues = {}
    all_eigenvectors = {}

    # Use the first sigma as seed system
    sigma_seed = sigma_list[0]

    # Seed: build Krylov subspace for (A - sigma_seed I)
    # We use a manual BiCG implementation with shift tracking
    logger.info("Shifted BiCG: seed sigma = {:.6f}, {} additional shifts".format(
        sigma_seed, len(sigma_list) - 1))

    # For each shift, solve via shift-invert + eigs
    # The key optimization: reuse preconditioner from seed system
    for i, sigma in enumerate(sigma_list):
        logger.info("Processing shift sigma = {:.6f} ({}/{})".format(
            sigma, i + 1, len(sigma_list)))

        def shifted_matvec(v, s=sigma):
            return A.matvec(v) - s * v

        A_shifted = LinearOperator(
            (vec_size, vec_size), matvec=shifted_matvec, dtype=complex
        )

        def inv_matvec(v):
            x, info = gmres(A_shifted, v, rtol=tol_linear, maxiter=300)
            return x

        A_inv = LinearOperator(
            (vec_size, vec_size), matvec=inv_matvec, dtype=complex
        )

        try:
            nus, vecs = eigs(A_inv, k=max_ev, which='LM')
            eigenvalues = 1.0 / nus + sigma
            # Eigenvalues found near this shift; report by magnitude.
            idx = np.argsort(-np.abs(eigenvalues))
            eigenvalues = eigenvalues[idx]
            eigvecs = np.array([
                vecs[:, j].reshape(norb, norb, Nx, Ny, Nz)
                for j in idx
            ])
        except Exception as e:
            logger.warning("Shift sigma={:.6f} failed: {}".format(sigma, e))
            eigenvalues = np.array([])
            eigvecs = np.array([])

        all_eigenvalues[sigma] = eigenvalues
        all_eigenvectors[sigma] = eigvecs

        if len(eigenvalues) > 0:
            logger.info("  Found eigenvalues: {}".format(
                np.array2string(eigenvalues.real, precision=6)))

    return all_eigenvalues, all_eigenvectors


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _save_results(output_dir, sigma, eigenvalue, eigenvalues_eig, kx_array, ky_array, kz_array,
                  gap_file="gap.dat", eigenvalue_file="eigenvalue.dat",
                  eigenvalue_match=None, provenance=None,
                  eigenvalue_note=None):
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
    """
    os.makedirs(output_dir, exist_ok=True)

    if sigma is not None:
        norb = sigma.shape[0]
        Nx, Ny, Nz = len(kx_array), len(ky_array), len(kz_array)
        filepath = os.path.join(output_dir, gap_file)
        logger.info("Saving gap function to {}".format(filepath))

        with open(filepath, "w") as fw:
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
            if eigenvalue is not None:
                fw.write("# Iteration eigenvalue\n")
                if eigenvalue_note:
                    for line in str(eigenvalue_note).splitlines():
                        fw.write("# {}\n".format(line))
                fw.write("{:.8e}\n".format(eigenvalue))
            if eigenvalues_eig is not None:
                fw.write("# Eigenvalue analysis\n")
                # The optional trailing `match` column is 1 when the eigenvector
                # lies (dominantly) in the channel's parity sector (even for
                # singlet, odd for triplet) and 0 when it lies in the opposite
                # sector. The leading index/Re/Im/|ev| columns are unchanged, so
                # position-based parsers (e.g. reading column 1 as Re) are
                # unaffected; only the column count grows from 4 to 5.
                if eigenvalue_match is not None:
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
    """
    # --- Parse parameters ---
    mode_param = input_dict["mode"]["param"]
    T = mode_param["T"]
    beta = 1.0 / T
    cell_shape = mode_param["CellShape"]
    sub_shape = mode_param.get("SubShape", cell_shape)
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

    # Dispatch to dynamic Eliashberg if requested
    if _eliashberg_frequency(input_dict) == "dynamic":
        _validate_dynamic_prereqs(input_dict)
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

    solver_mode = eli_param.get("solver_mode", "iteration")
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
    gap_file = eli_param.get("output_gap", "gap.dat")
    eigenvalue_file = eli_param.get("output_eigenvalue", "eigenvalue.dat")

    output_dir = input_dict["file"]["output"]["path_to_output"]

    logger.info("=== Linearized Eliashberg Equation Solver ===")
    logger.info("T = {}, beta = {:.4f}".format(T, beta))
    logger.info("Cell: {}x{}x{}, Sub: {}x{}x{}, Grid: {}x{}x{}".format(
        Lx, Ly, Lz, Bx, By, Bz, Nx, Ny, Nz))
    logger.info("Nmat = {}, filling = {}".format(nmat, n_filling))
    logger.info("solver_mode = {}, chi0q_mode = {}".format(solver_mode, chi0q_mode))

    # --- Step 2: Read input files ---
    geom_info, hr, interactions = _read_interaction_files(input_dict)
    norb = geom_info["norb"]

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

    # --- Step 3: Setup k-mesh ---
    kx_array = np.linspace(0, 2.0 * np.pi, Nx, endpoint=False)
    ky_array = np.linspace(0, 2.0 * np.pi, Ny, endpoint=False)
    kz_array = np.linspace(0, 2.0 * np.pi, Nz, endpoint=False)

    # --- Step 4: Build Hamiltonian and diagonalize ---
    logger.info("Building Hamiltonian in k-space...")
    epsilon_k = _build_hamiltonian_k(kx_array, ky_array, kz_array, hr, norb)

    logger.info("Diagonalizing Hamiltonian...")
    eigenvalues, eigenvectors = _calc_eigenvalues(epsilon_k)

    # --- Step 5: Determine chemical potential ---
    logger.info("Determining chemical potential...")
    mu = _determine_mu(eigenvalues, beta, n_filling, norb)
    logger.info("mu = {:.6f}".format(mu))

    # --- Step 7: Build interaction in k-space ---
    logger.info("Building interactions in k-space...")
    inter_k = _build_interaction_k(kx_array, ky_array, kz_array, interactions, norb)

    # Pre-built Eliashberg operator + provenance record; only the bond path
    # sets them (every other path builds its operator from Vs_q/G2 as before).
    bond_operator = None
    bond_provenance = None
    bond_attribution = None

    if chi0q_mode == "flex":
        # --- FLEX mode: use pre-computed dressed susceptibilities ---
        logger.info("FLEX mode: loading dressed susceptibilities and Green's function")

        chis, chic, green_dressed, chi_convention = _load_flex_susceptibilities(
            input_dict, norb, Nx, Ny, Nz)
        logger.info("FLEX susceptibility convention: {}".format(chi_convention))

        # Use dressed Green's function if available, otherwise use bare
        if green_dressed is not None:
            green_kw = green_dressed
            logger.info("Using FLEX dressed Green's function")
        else:
            logger.info("Calculating bare Green's function G(k, iwn)...")
            green_kw = _calc_green(eigenvalues, eigenvectors, mu, beta, nmat)

        # Compute pairing vertex from FLEX susceptibilities
        logger.info("Computing FLEX vertices (pairing_type={})...".format(pairing_type))
        Vs_q = _compute_vertices_flex(chis, chic, inter_k, norb, Nx, Ny, Nz,
                                      pairing_type=pairing_type,
                                      convention=chi_convention)
        logger.info("FLEX vertex shape: {}".format(Vs_q.shape))

    elif use_bond_channels:
        # --- Bond-resolved mode: the enlarged (orbital-pair)x(bond) ladder ---
        # This path builds its own bubble (bond_bubble) directly from the
        # Green function, so no chi0q file/tensor is read: the m=m'=0 block of
        # chi-bar IS the ordinary chi0q, and the m!=0 blocks (which no chi0q
        # file carries) are what the bond channels add.
        #
        # Resolve the bond topology and run the resource preflight (S3.2)
        # BEFORE allocating green_kw (review fix I-2): green_kw's own bytes
        # are folded into the preflight's estimate, and a large-Nmat run must
        # be refused here, not after the (possibly large) Green-function
        # allocation has already happened.
        bond_set = bond_channels.resolve_interactions(
            interactions["CoulombInter"], geom_info["rvec"], norb,
            bond_max_shells=bond_max_shells)
        logger.info("Bond channels: B = %d, Delta r = %s",
                    bond_set.n_channels, list(bond_set.delta_r))
        if bond_set.dropped:
            logger.info("Bond channel provenance (dropped): %s",
                        list(bond_set.dropped))
        # The frequency grid the bond bubble will really be built on: for an
        # external bond_green that is the FILE's Nmat, which therefore has to
        # be known BEFORE the cap is applied -- otherwise an over-large
        # configured Nmat could refuse a run that fits comfortably (the file
        # is never reached). The header peek reads the npz member's shape
        # only, so the preflight still precedes every large allocation.
        nmat_bond = nmat
        if bond_green is not None:
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
        _bond_resource_preflight(norb, bond_set, Nx, Ny, Nz, nmat_bond,
                                 bond_memory_cap_gb)

        if bond_green is not None:
            # Externally supplied Green function (spec Goal): the milestone
            # green is FLEX-dressed and self-consistent, which the bare
            # transfer-Hamiltonian green is not. The bond path never reads a
            # chi file, so this is the only external input it takes.
            green_kw, nmat_green = _load_green_npz(
                bond_green, norb, Nx, Ny, Nz, label="bond_green")
            # Same requirement on the loaded count: when the header peek could
            # not resolve nfreq, this is the first place the file's real Nmat
            # is known, so the guard has to be repeated here rather than
            # trusted from above.
            _validate_bond_green_nmat(nmat_green, bond_green)
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
                                         nmat_green, bond_memory_cap_gb)
        else:
            logger.info("Calculating Green's function G(k, iwn)...")
            green_kw = _calc_green(eigenvalues, eigenvectors, mu, beta, nmat)
        bond_operator, bond_provenance, bond_attribution = _build_bond_operator(
            bond_set, green_kw, interactions, inter_k, geom_info, norb,
            kx_array, ky_array, kz_array, beta, pairing_type,
            bond_max_shells=bond_max_shells,
            bond_memory_cap_gb=bond_memory_cap_gb,
            precondition_opts=bond_precondition_opts,
            green_source=bond_green)
        Vs_q = None

    else:
        # --- Standard RPA mode ---

        # Step 6: Calculate bare Green's function
        logger.info("Calculating Green's function G(k, iwn)...")
        green_kw = _calc_green(eigenvalues, eigenvectors, mu, beta, nmat)

        # Step 1: Load or compute chi0q
        if chi0q_mode == "calc":
            chi0q_raw = _calc_chi0q_internal(input_dict, chi0q_tensor=chi0q_tensor,
                                                precomputed_mu=mu)
            # internally computed chi0q always carries the full frequency grid
            static_index = None
        else:
            chi0q_raw, static_index = _load_chi0q(input_dict)

        # Step 8: Convert chi0q format; the frequency axis is last in the
        # reference format, so read its length after the conversion (a 6D
        # input can be either raw H-wave (nmat first) or ref (nmat last))
        chi0q = _convert_chi0q_to_ref_format(chi0q_raw, norb, Nx, Ny, Nz)
        nmat_chi0q = chi0q.shape[-1]
        logger.info("chi0q converted to shape: {}".format(chi0q.shape))

        # Step 9: Compute RPA vertices
        logger.info("Computing RPA vertices (pairing_type={})...".format(pairing_type))
        vertex_result = _compute_vertices(chi0q, inter_k, norb, Nx, Ny, Nz, nmat_chi0q,
                                          pairing_type=pairing_type,
                                          static_index=static_index)
        if isinstance(vertex_result, tuple):
            Pc_q, Ps_q = vertex_result
            Vs_q = Pc_q + Ps_q
            logger.info("Simple mode: Pc + Ps vertex, shape {}".format(Vs_q.shape))
        else:
            Vs_q = vertex_result
            logger.info("General mode: 4-index V^s vertex, shape {}".format(Vs_q.shape))

    # --- Step 10: Compute G2 ---
    # The bond operator carries its own pair bubble (the same construction as
    # _calc_g2, applied inside bond_channels.make_bond_kernel), so G2 is not
    # needed -- and must not be built -- on that path.
    if use_bond_channels:
        G2 = None
    else:
        logger.info("Computing G2...")
        G2 = _calc_g2(green_kw, beta)

    # --- Step 11: Initialize gap function ---
    sigma_init = _initialize_gap(init_gap_mode, norb, kx_array, ky_array, kz_array)

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
    iteration_spectral_shift = (eli_param.get("spectral_shift")
                                if solver_mode in ("iteration", "both")
                                else None)

    # Rayleigh-validation record of a SHIFTED power iteration (empty on the
    # unshifted default path); decides whether the number written to
    # eigenvalue.dat may be called an eigenvalue at all.
    iteration_info = {}

    if solver_mode in ("iteration", "both"):
        logger.info("=== Self-consistent iteration ===")
        sigma_result, eigenvalue_iter, converged, n_iter = _solve_iteration(
            green_kw, Vs_q, G2, sigma_init, norb,
            max_iter=max_iter, alpha=alpha, tol=tol, pairing_type=pairing_type,
            operator=bond_operator,
            spectral_shift=iteration_spectral_shift,
            info_out=iteration_info,
        )
        logger.info("Iteration result: eigenvalue = {:.6f}, converged = {}, n_iter = {}".format(
            eigenvalue_iter, converged, n_iter))

    if solver_mode in ("eigenvalue", "both"):
        logger.info("=== Eigenvalue analysis ===")
        eigenvalues_eig, eigenvectors_eig = _solve_eigenvalue(
            Vs_q, G2, norb, Nx, Ny, Nz,
            num_eigenvalues=num_eigenvalues,
            method=eigenvalue_method,
            sigma_shift=eli_param.get("sigma_shift"),
            spectral_shift=eli_param.get("spectral_shift"),
            operator=bond_operator,
        )
        # The kernel preserves parity; promote the eigenpairs whose gap has the
        # requested channel parity (singlet even / triplet odd) so the reported
        # leading eigenpair is the physical solution for this channel.
        eigenvalues_eig, eigenvectors_eig = _reorder_eigenpairs_by_parity(
            eigenvalues_eig, eigenvectors_eig, pairing_type)
        # Flag whether each eigenvector lies (dominantly) in the channel's
        # parity sector (even for singlet, odd for triplet). _is_gap_parity is a
        # threshold test, so this is "dominantly the channel parity" rather than
        # an exact label; for a centrosymmetric model parity is an exact kernel
        # symmetry and the opposite-parity modes are unphysical for the channel
        # (e.g. an even mode in the triplet kernel is Pauli-forbidden).
        eigenvalue_match = np.array(
            [_is_gap_parity(g, pairing_type) for g in eigenvectors_eig])
        logger.info("Leading eigenvalues:")
        for i, ev in enumerate(eigenvalues_eig):
            tag = "" if eigenvalue_match[i] else "  [opposite-parity sector]"
            logger.info("  {:3d}: {:.6f} (|ev| = {:.6f}){}".format(
                i, ev.real, abs(ev), tag))

        # Use leading eigenvector as gap if no iteration result
        if sigma_result is None:
            sigma_result = eigenvectors_eig[0]
            eigenvalue_iter = eigenvalues_eig[0].real

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
        bond_provenance["lambda_rayleigh_solver_mode"] = solver_mode
        bond_provenance["lambda_rayleigh_converged"] = (
            str(bool(converged)) if solver_mode in ("iteration", "both")
            else "n/a")

        # The iteration path returns ||A sigma||, an UNSIGNED norm. On a
        # repulsive-dominant bond kernel the leading eigenvalue is NEGATIVE, so
        # the power iterate flips sign every step: diff ~ 2 forever and the
        # loop can never converge. That is structural, not a max_iter setting
        # -- warn loudly and point at the signed quantity (review fix C1b).
        # With a spectral_shift the iteration runs on K + sigma*I and the
        # returned value IS the signed eigenvalue, so this warning does not
        # apply (a mismatch there means the run simply did not converge, which
        # is already warned about in _solve_leading).
        if (iteration_spectral_shift is None
                and solver_mode in ("iteration", "both")
                and eigenvalue_iter is not None):
            lam_signed = attr["lambda"]
            if (np.sign(lam_signed) != np.sign(eigenvalue_iter)
                    or abs(lam_signed - eigenvalue_iter) > 1.0e-6):
                logger.warning(
                    "[eliashberg] solver_mode='%s' reports %.8e, which is the "
                    "UNSIGNED norm of the power iterate (||A sigma||), NOT a "
                    "signed eigenvalue; the signed physical quantity is "
                    "lambda_rayleigh = %.8e (converged = %s). Iteration mode "
                    "CANNOT converge on a repulsive-dominant bond kernel: its "
                    "dominant eigenvalue is negative, so the iterate flips "
                    "sign every step and the normalized difference stays "
                    "~2 regardless of max_iter. Set [eliashberg] "
                    "spectral_shift=\"auto\" to iterate on the shifted kernel "
                    "K + sigma*I (the reported lambda is then signed), or use "
                    "solver_mode=\"eigenvalue\" for the bond path.",
                    solver_mode, eigenvalue_iter, lam_signed, bool(converged))

    # --- Step 12c: opt-in bond diagnostics (S7.7 character analysis) ---
    # Additive comment lines only; off by default, so a run without the flag is
    # byte-identical to what it was before the diagnostics existed.
    if (bond_diagnostics and bond_provenance is not None
            and bond_attribution is not None):
        _bond_diagnostics_record(
            bond_provenance, sigma_result, eigenvalues_eig,
            bond_attribution["weight"], norb,
            kx_array, ky_array, kz_array)

    # --- Step 13: Save results ---
    # Label the iteration number in the file itself, so it can never be
    # silently confused with the signed lambda_rayleigh sitting next to it
    # (review fix C1c).
    #
    # The note is written ONLY when it carries information, i.e. when the run
    # opted into something that changes what the number means:
    #   * an explicit spectral_shift -- the value is then the SIGNED eigenvalue
    #     of K (or, if that could not be validated, a shifted iterate-norm
    #     ESTIMATE, see below), not the historical unsigned iterate norm; or
    #   * bond_channels -- eigenvalue.dat then also carries lambda_rayleigh,
    #     and the two numbers must be told apart.
    # With both inactive the meaning is exactly the historical one and there is
    # nothing to say, so eigenvalue.dat stays byte-for-byte the legacy file
    # (tests/test_sc_legacy_golden.py pins this against commit 712b1a0).
    eigenvalue_note = None
    if (solver_mode in ("iteration", "both") and eigenvalue_iter is not None
            and (iteration_spectral_shift is not None or use_bond_channels)):
        if iteration_spectral_shift is not None:
            # Shifted power iteration: the value is the SIGNED eigenvalue of
            # the original kernel only when the Rayleigh check validated it;
            # otherwise it is an iterate-norm ESTIMATE and must be labelled as
            # not-an-eigenvalue (see _validate_shifted_eigenvalue).
            eigenvalue_note = _shifted_eigenvalue_note(
                solver_mode, iteration_spectral_shift, converged, n_iter,
                iteration_info)
        else:
            note = ("solver_mode='{}': the value below is the UNSIGNED norm of "
                    "the power iterate (||A sigma||), not a signed eigenvalue; "
                    "the power iteration {} after {} steps."
                    .format(solver_mode,
                            "converged" if converged else "did NOT converge",
                            n_iter))
            if bond_provenance and "lambda_rayleigh" in bond_provenance:
                note += (" The signed physical eigenvalue of this run is "
                         "lambda_rayleigh = {} above. Iteration mode does not "
                         "converge for a repulsive-dominant bond kernel "
                         "(negative dominant eigenvalue) unless you set "
                         "[eliashberg] spectral_shift; otherwise use "
                         "solver_mode=\"eigenvalue\"."
                         .format(bond_provenance["lambda_rayleigh"]))
            eigenvalue_note = note

    _save_results(
        output_dir, sigma_result, eigenvalue_iter, eigenvalues_eig,
        kx_array, ky_array, kz_array,
        gap_file=gap_file, eigenvalue_file=eigenvalue_file,
        eigenvalue_match=eigenvalue_match,
        provenance=bond_provenance,
        eigenvalue_note=eigenvalue_note
    )

    logger.info("=== Done ===")


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

    calc_eliashberg(input_dict)


if __name__ == "__main__":
    main()
