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
from typing import NamedTuple

import numpy as np

from hwave.solver.vertex_table import sc_coefficients
from hwave.solver.kgrid import reverse_fft_axes
from hwave.solver.declarations import symmetrise_k
from hwave.solver import npy_header as _npy_header
from numpy.fft import fftn, ifftn
from scipy.optimize import bisect
from scipy.sparse.linalg import LinearOperator, eigs, bicgstab, gmres, lgmres

import hwave
import hwave.qlmsio.wan90 as wan90
from hwave.solver.rpa import validate_chi0q_index_convention, TAIL_ENDPOINT_CONVENTION
from hwave.solver.ir_axis import is_ir_native, ir_native_meta
from hwave.solver import backend
from hwave.solver import bond_channels
from hwave.solver import bubble
from hwave.solver import green as green_mod

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
    "bond_coeff_tail",
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


def _reject_bond_channels_dynamic(eli_param):
    """Reject ``bond_channels=true`` on the dynamic Eliashberg entry point.

    The bond-resolved vertex is implemented for the STATIC path only (spec S5
    guard / non-goal S2): the dynamic kernel would need the bond channels
    carried through the frequency-resolved vertex, which is a separate spec.
    Fail loudly rather than silently running the scalar dynamic vertex.

    The flag goes through the same validator as the static path, so a typo
    ("ture") is refused here too instead of quietly reading as false.
    """
    if _bond_bool_option(eli_param, "bond_channels", False):
        raise ValueError(
            "[eliashberg] bond_channels=true is not supported with "
            "frequency='dynamic': the bond-resolved pairing vertex is "
            "implemented for the STATIC linearized Eliashberg path only. "
            "Porting the bond kernel to eliashberg_dynamic.py is a deferred "
            "follow-up; use frequency='static' (or bond_channels=false).")


def _validate_dynamic_prereqs(input_dict):
    """Validate prerequisites for dynamic Eliashberg calculation.

    Parameters
    ----------
    input_dict : dict
        Parsed TOML configuration dictionary.

    Raises
    ------
    ValueError
        If bond_channels is requested, if chi0q_mode is not "flex", or if
        Nmat is odd.
    """
    eli = input_dict.get("eliashberg", {})
    _reject_bond_channels_dynamic(eli)
    if eli.get("chi0q_mode") != "flex":
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

def _resolve_bond_coeff_tail(eli_param):
    """Resolve ``[eliashberg] bond_coeff_tail`` (unified-bubble-kernel spec,
    "Tail default for the bond path").

    Resolution rule: key ABSENT -> requested value ``1.0`` (the physical
    asymptotic coefficient of a normal Green function; source
    ``"default"``); key PRESENT -> its value verbatim, validated exactly as
    RPA validates ``[mode.param] coeff_tail`` (rpa.py:1318-1332: a real
    scalar, type-strict -- booleans are rejected outright, since
    ``float(True) == 1.0`` would silently pass the check -- and finite; any
    finite real value is otherwise accepted, source ``"config"``).

    This intentionally differs from the RPA route's historical default
    ``coeff_tail = 0.0``: the bond path is pre-first-release code with no
    legacy users, so it defaults to the physically correct tail-on
    behaviour; the RPA route's own default is untouched by this function.

    Returns
    -------
    (requested, source) : (float, str)
        ``source`` is ``"default"`` or ``"config"``.
    """
    import numbers

    if "bond_coeff_tail" not in eli_param or eli_param["bond_coeff_tail"] is None:
        return 1.0, "default"

    raw = eli_param["bond_coeff_tail"]
    if isinstance(raw, (bool, np.bool_)) or not isinstance(raw, numbers.Real):
        raise ValueError(
            "[eliashberg] bond_coeff_tail must be a real number, got "
            "{!r}".format(raw))
    value = float(raw)
    if not np.isfinite(value):
        raise ValueError(
            "[eliashberg] bond_coeff_tail must be finite, got {}".format(value))
    return value, "config"


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
        simultaneously in the vertex assembly -- ``chi_bar``, ``S_bond``,
        ``C_bond``, ``chi_s``, ``chi_c``, the fluctuation vertex ``F_q``, its
        two matrix-product temporaries, and ``Gamma_hat = ifftn(F_q)``.
        ``bond_bubble``'s output ``chi_bar`` is the first of these, so it is
        NOT counted again in ``bubble_bytes``.
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
    ``vertex_bytes``
        ``bare_bond_vertices``'s non-q-resolved ``ND x ND`` working set --
        ``BARE_VERTEX_ND2_BUFFERS`` buffers (``P``, ``D``, ``Dh``, ``Id``,
        ``Q_s``, ``Q_t``, ``B_s``, ``B_t``, ``Vpp_s``, ``Vpp_t``), NOT
        multiplied by ``N_q`` since none of them carry a q-axis (unlike
        ``S_bond``/``C_bond``, which are already inside ``q_bytes``). This is
        alive concurrently with ``q_bytes`` (``chi_bar`` from ``bond_bubble``
        is still held) and ``carrier_bytes`` (``full_sc`` is still held), so
        it is additive in ``peak``. The buffer-by-buffer list is documented
        beside the constant in ``hwave.solver.bond_channels``; imported from
        there so the two sides cannot silently desync.
    ``peak``
        ``q_bytes + bubble_bytes + vertex_bytes + carrier_bytes`` -- the SUM
        combine rule preserved unchanged from before the tail correction;
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
    vertex_bytes = bond_channels.BARE_VERTEX_ND2_BUFFERS * ND * ND * itemsize
    carrier_bytes = (4 if tail_on else 2) * (norb ** 2) * unit
    return {"nd": nd, "B": B, "ND": ND, "Nq": Nq,
            "chi_bar_bytes": chi_bar_bytes,
            "q_bytes": q_bytes,
            "bubble_bytes": bubble_bytes,
            "vertex_bytes": vertex_bytes,
            "carrier_bytes": carrier_bytes,
            "peak": q_bytes + bubble_bytes + vertex_bytes + carrier_bytes}


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
    vertex_bytes = est["vertex_bytes"]
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
            "{:.3f} GB of bond-bubble work buffers plus {:.3f} GB of bare-"
            "vertex ND x ND temporaries plus {:.3f} GB for the Green "
            "carrier (tail_on={}). Reduce bond_max_shells (fewer channels), "
            "reduce the k-grid/Nmat, or raise bond_memory_cap_gb.".format(
                peak / 1.0e9, cap_gb, B, list(bond_set.delta_r), ND, Nq,
                q_bytes / 1.0e9, bubble_bytes / 1.0e9, vertex_bytes / 1.0e9,
                carrier_bytes / 1.0e9, tail_on))
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
        # Hermitize PER SHELL, BEFORE the q-phase multiply (Finding 1, #151
        # ED-adjudication campaign; tests/test_bond_vs_ed_oracle.py::
        # TestNullDirectionSolverSide). resolve_interactions guarantees,
        # EXACTLY, that bond_set.v_bond[reverse[m]][b, a] ==
        # conj(bond_set.v_bond[m][a, b]): the declared channel V_ab(R_m) and
        # its reversal partner V_ba(-R_m) multiply the SAME physical operator
        # sum_j n_{j,a} n_{j+R_m,b} (site relabeling j' = j+R_m turns the
        # partner's sum_j n_{j,b} n_{j-R_m,a} into the identical sum), so
        # Hermiticity forces their combined Hamiltonian coefficient
        # V_ab(R_m) + V_ba(-R_m) == v_bond[m][a,b] + conj(v_bond[m][a,b]) ==
        # 2*Re(v_bond[m][a,b]) -- real, and INDEPENDENT of q (this is the full
        # per-shell coefficient the Hartree term wants, unlike the ph Fock
        # slot in bare_bond_vertices which only wants HALF of it; see that
        # function's comment for why the two channels differ). This holds
        # identically for a genuinely complex, phase-carrying single-sided
        # declaration (not just a synthetic null-direction probe): a complex
        # CoulombInter value always collapses to twice its real part in this
        # density-density Hartree term.
        #
        # Before this fix, V_q summed the RAW (possibly complex) per-shell
        # value: the null-direction pair V_ab(+R)=V+i*eps / V_ba(-R)=V-i*eps
        # left a 2*eps-responsive residual in this Hartree block feeding
        # bare_bond_vertices' C_bond (that function's m!=0 Fock diagonal
        # carries the other, eps-magnitude half of the same finding).
        # Hermitizing per shell HERE, before multiplying by the complex
        # q-phase, keeps V_q genuinely q-dispersive/complex at q!=0 --
        # collapsing to .real only AFTER the phase sum (or on the full V_q)
        # would wrongly discard the physical q-dependent phase.
        # bond_channels.bare_bond_vertices' local-block hartree_rneq0
        # subtraction mirrors this exact Hermitization so the R!=0 Hartree it
        # removes from C0_loc still matches what is added here (keeps
        # Vpp_s/Vpp_t null-invariant).
        #
        # .real is bit-identical to the raw value whenever it is already real
        # (every direction adjudicated by this campaign, and every #82
        # golden-test coupling), so this is a provable no-op there.
        V_q += bond_set.v_bond[m][:, :, np.newaxis, np.newaxis, np.newaxis].real \
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


def _build_bond_operator(bond_set, green_carrier, interactions, inter_k,
                         geom_info, norb, kx_array, ky_array, kz_array, beta,
                         pairing_type, *, bond_max_shells, bond_memory_cap_gb,
                         precondition_opts=None, green_source=None,
                         g2_tail=False, bond_coeff_tail_source="default"):
    """Run the bond-resolved physics chain and return the Eliashberg operator.

    bond_bubble -> Case-2 corrected m=0 blocks -> bare_bond_vertices ->
    dress_bond -> make_bond_kernel (spec S4.2-S4.5).

    ``bond_set`` must already be resolved (``bond_channels.
    resolve_interactions``) and pre-flighted (``_bond_resource_preflight``) by
    the caller BEFORE ``green_carrier`` was built (review fix I-2: the
    preflight must run before the potentially-large Green-function
    allocation, not after -- see ``calc_eliashberg``'s bond-channel dispatch
    branch). This function therefore no longer resolves the interaction or
    runs the preflight itself.

    ``green_carrier`` : :class:`_BondGreen`
        The bubble reads ONLY ``deflated_kw``/``green0_tail`` (unified-
        bubble-kernel spec: "the bubble receives (deflated_kw, green0_tail);
        EVERY downstream consumer receives full_kw" -- here ``full_sc``, the
        sc-layout equivalent); every other consumer in this function
        (``make_bond_kernel_parts``, ``pair_weight``) reads
        ``green_carrier.full_sc`` -- the bare Green function on the same
        inputs, independent of the tail coefficient.

    ``precondition_opts`` (from ``_read_bond_config``) is forwarded verbatim
    to ``bond_channels.check_hermitian_preconditions`` -- the
    ``bond_precondition_atol/rtol/dense_limit`` knobs of ``[eliashberg]``
    (review fix M7).

    ``green_source`` is provenance only: the path the Green function was
    loaded from (``[eliashberg] bond_green``), or ``None`` when the caller
    built the bare Green function from the transfer Hamiltonian. It selects
    the wording of the recorded approximation level -- the two are
    physically different approximations and the record must say which one
    ran.

    ``bond_max_shells`` and ``bond_memory_cap_gb`` are likewise PROVENANCE
    ONLY here: the resource preflight and the shell cutoff itself already
    happened in the caller (``resolve_interactions``/
    ``_bond_resource_preflight``, before ``bond_set``/``green_carrier`` were
    built), so by the time this function runs, neither value feeds into the
    computed operator -- they are only echoed into the returned
    ``provenance`` dict. They are keyword-only (review fix: honest
    signature) precisely so a call site cannot mistake them for
    physics-affecting positional arguments. (Past review fix: an earlier
    ``nmat`` parameter was dropped from this signature for the same reason --
    it was never read in the body, since ``green_carrier.full_sc.shape[-1]``
    already fixes the Matsubara grid the bond bubble is built on, and
    ``bond_green`` can make that grid differ from the ``nmat`` the caller
    resolved the preflight with.)

    Returns
    -------
    (A, vec_size) : tuple
        The ``scipy.sparse.linalg.LinearOperator`` kernel and its vector size,
        in the same convention as ``_make_kernel_operator``.
    provenance : dict
        Additive output-provenance record (spec S5; unified-bubble-kernel
        spec's "Tail default for the bond path" adds the
        ``bond_coeff_tail_*``/``bond_endpoint_convention``/
        ``bond_green_source``/``external_full_green_assumed`` keys).
    attribution : dict
        Context for the ``lambda = lambda^pp + lambda^fl`` decomposition of
        spec S4.5: the two kernel parts (``pp``/``fl``), the full kernel and
        the ``sqrt(GG)`` metric (``weight``), consumed by
        ``bond_channels.attribute_lambda`` once a gap has converged.
    """
    Nx, Ny, Nz = len(kx_array), len(ky_array), len(kz_array)
    green_kw = green_carrier.full_sc

    logger.info("Computing the bond-resolved bubble chi-bar...")
    chi_bar = bubble.bond_bubble_static(
        green_carrier.deflated_kw, green_carrier.green0_tail, beta, bond_set,
        spatial_shape=(Nx, Ny, Nz))
    # bond_bubble_static returns the flattened-spatial (nvol, ND, ND)
    # canonical shape; every downstream consumer here (dress_bond,
    # bare_bond_vertices' S/C blocks, ...) still expects the unflattened
    # sc.py spatial layout (Nx, Ny, Nz, ND, ND) -- the same reshape
    # bond_channels.bond_bubble's compatibility wrapper applies.
    nD = chi_bar.shape[-1]
    chi_bar = chi_bar.reshape(Nx, Ny, Nz, nD, nD)

    S0_q, C0_q = _build_bond_m0_blocks(bond_set, interactions, inter_k, norb,
                                       kx_array, ky_array, kz_array)
    S_bond, C_bond, Vpp_s, Vpp_t = bond_channels.bare_bond_vertices(
        bond_set, S0_q, C0_q, norb)
    chi_s, chi_c = bond_channels.dress_bond(chi_bar, S_bond, C_bond)

    logger.info("Building the bond-resolved Eliashberg operator "
                "(pairing_type=%s)...", pairing_type)
    A, A_fl, A_pp, vec_size = bond_channels.make_bond_kernel_parts(
        chi_s, chi_c, S_bond, C_bond, Vpp_s, Vpp_t, green_kw, bond_set,
        pairing_type, beta, g2_tail=g2_tail)

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
        "bond_n_channels": int(bond_set.n_channels),
        "bond_delta_r": [list(dr) for dr in bond_set.delta_r],
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


def _load_chi0q(input_dict, norb=None):
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
    enable_spin_orbital = _resolve_spin_orbital_flag(input_dict)
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
        # CaseInsensitiveDict, like _read_interaction_files (round-5
        # review): an exact-case check classified a 'coulombinter' run as
        # reduced, silently omitting the very components the comment
        # above says inter-orbital interactions require. 'Coulomb' counts
        # too -- its split can produce a CoulombInter part.
        from requests.structures import CaseInsensitiveDict
        files = CaseInsensitiveDict(info_inputfile.get("interaction", {}))
        has_interorbital = any(k in files for k in
                              ["Hund", "Exchange", "CoulombInter",
                               "Ising", "PairHop", "Coulomb"])
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
    # Forward the RESOLVED boolean, not the raw value: RPA mixes truthiness
    # and == True checks on this flag, so a string "false" would diverge
    # internally (round-1 review of #83's guard).
    info_mode_rpa["enable_spin_orbital"] = _resolve_spin_orbital_flag(
        input_dict)

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
    from hwave.solver.declarations import validate_hermitian_closure

    for itype in interaction_types:
        if itype in files:
            f = os.path.join(path_to_input, files[itype])
            logger.info("Reading {} from {}".format(itype, f))
            tbl = wan90.read_w90(f)
            # issue #93: fail fast on non-Hermitian-closed declarations,
            # with the same rule and tolerance as the k-space reader
            validate_hermitian_closure(itype, tbl, source=f)
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
    # Solver-core convention (rpa.py _make_ham_trans, issue #133:
    # ifftn * nvol == e^{+ikR}, the documented Wannier90-style sign shared
    # with UHFk): epsilon[a,b](k) = sum_R t_R[a,b] e^{+ikR}. This keeps
    # sc-built quantities element-wise consistent with arrays loaded from
    # FLEX/RPA output files, which carry the same k labeling since the
    # #133 alignment. (The [orb1,orb2] placement is unchanged; only the
    # Fourier sign moved with the solver core.)
    for (irvec, orbvec), value in hr.items():
        orb1, orb2 = orbvec
        Rx, Ry, Rz = irvec
        epsilon_k[orb1, orb2, :, :, :] += value * np.exp(
            +1j * (kx_mesh * Rx + ky_mesh * Ry + kz_mesh * Rz)
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

    def _to_k(value_r, transpose=True):
        # Same Fourier phase as the solver core, e^{+iqR} since the #133
        # sign alignment, but the ORBITAL PAIR is stored transposed: an
        # entry (R, (a, b)) lands at [b, a].
        #
        # The interaction is a four-index object, and its MATRIX form carries a
        # pair-index transpose that the one-body Hamiltonian does not. H-wave's
        # own paper (arXiv:2308.00324) makes this explicit: Eq.(12) defines the
        # tensor W_ij^{aa'bb'} with the first pair at site i and the second at
        # site j, Eq.(16) puts the first pair at momentum +q, and Eq.(21) then
        # defines the matrix used in the RPA equation as
        #     [W(q)]^{ab} = W_q^{ba}
        # so for a density-density term W^{aabb} = V_ab(R) the matrix element is
        # [W(q)]^{(aa),(bb)} = V_ba(q), i.e. the matrix is V(q)^T.
        #
        # `rpa.py::_make_ham_inter` has always built that (it stores the entry's
        # first orbital in the SECOND pair slot); this builder did not, so the
        # susceptibility was solved with one orientation and the pairing vertex
        # assembled from the other (issue #96).
        #
        # Both consumers need the transpose, confirmed by exact
        # diagonalization on a bond set with V_ab(R) != V_ba(R):
        #   * the RPA ladder [I + chi0 W]^-1 chi0 -- solved here too, in
        #     _compute_vertices_simple -- selects V^T with a residual that
        #     scales linearly in V (pure O(V^2) truncation), against 98% for V;
        #   * the bare pair-scattering amplitude
        #     <k' a up, -k' b down| H_int |k a up, -k b down>, which is exact at
        #     first order, matches V^T to 6e-16, against 99% for V.
        #
        # NOTE: `_build_hamiltonian_k` is a one-body object and correctly keeps
        # the plain [a, b] placement. Do not "align" this builder to it.
        #
        # PairHop is EXCLUDED (transpose=False below), and unlike the rest that
        # is now a MEASURED result rather than a gap. It is the one type
        # `rpa.py` does not place through `_append_inter`: `_append_pairhop`
        # uses the slots (b, a, a, b) rather than the density-density
        # (b, b, a, a), and in the S/C matrices it lands in the pair-ANTIdiagonal
        # Case 4 rather than Case 1/3, so the density-density determination
        # above says nothing about it.
        #
        # Issue #100 settled it end to end (TestPairHopEndToEnd in
        # tests/test_vertex_orientation.py), with no index reasoning anywhere
        # between the input file and the answer: feed this builder an on-site
        # PairHop matrix, build S from it, form the solver's own linear-order
        # response -chi0 . S . chi0 with chi0 from `rpa.py`'s kernel and an
        # EXACT Green function, and compare against the cross-spin response of
        # the same Hamiltonian by exact diagonalization. The untransposed
        # placement -- i.e. S[(a,b),(b,a)] = P[a,b] -- reproduces it, and the
        # test locates the residual by checking chi0 against the exact bubble
        # in the same run, so what is left is the imaginary-time grid. Refining
        # that grid off-line, 192 -> 384 -> 768 points, takes the residual
        # 2.3e-3 -> 1.1e-3 -> 5.5e-4 (the test itself runs only the first);
        # transposing PairHop sits at 0.55 and does not move.
        #
        # Beware the plausible-looking derivation that says otherwise. The
        # vertex sits BETWEEN two chi0 factors, so its row index lives in
        # chi0's COLUMN convention and its column index in chi0's ROW
        # convention. chi0's row pair is reversed relative to its bilinear and
        # its column pair is not, so the reversal lands on the vertex's COLUMN
        # side. Applying it to the row side instead gives the transpose, and is
        # wrong.
        val_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        for (irvec, orbvec), value in value_r.items():
            orb1, orb2 = orbvec
            Rx, Ry, Rz = irvec
            i1, i2 = (orb2, orb1) if transpose else (orb1, orb2)
            val_k[i1, i2, :, :, :] += value * np.exp(
                +1j * (kx_mesh * Rx + ky_mesh * Ry + kz_mesh * Rz)
            )
        return val_k

    inter_k = {}
    for itype in ["CoulombIntra", "CoulombInter", "Hund", "Exchange",
                  "Ising", "PairLift", "PairHop"]:
        if itype in interactions:
            inter_k[itype] = _to_k(interactions[itype],
                                   transpose=(itype != "PairHop"))

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


def _green_sc_to_canonical(g):
    """Convert sc.py's Green-function layout to the shared kernel's
    canonical layout.

    ``(p, p, Nx, Ny, Nz, nmat) -> (1, nmat, Nx*Ny*Nz, p, p)``: transpose to
    ``(nmat, Nx, Ny, Nz, p, p)`` order, reshape the spatial axes to
    ``nvol`` in C order, and insert a leading length-1 block axis (the
    canonical layout's ``nblock`` -- sc.py's Green function never carries a
    block axis of its own). ``Nx, Ny, Nz`` are read directly off ``g``'s
    shape, so no separate ``spatial_shape`` argument is needed on this
    direction. ALWAYS returns a copy (``np.ascontiguousarray``) -- never a
    view aliasing ``g`` (unified-bubble-kernel spec, "Layout adapters").
    """
    g = np.asarray(g)
    if g.ndim != 6:
        raise ValueError(
            "_green_sc_to_canonical: expected ndim=6 (p, p, Nx, Ny, Nz, "
            "nmat), got ndim={} shape={}".format(g.ndim, g.shape))
    p, p2, Nx, Ny, Nz, nmat = g.shape
    if p2 != p:
        raise ValueError(
            "_green_sc_to_canonical: the two orbital axes must be equal; "
            "got shape {}".format(g.shape))
    transposed = g.transpose(5, 2, 3, 4, 0, 1)          # (nmat,Nx,Ny,Nz,p,p)
    canonical = transposed.reshape(nmat, Nx * Ny * Nz, p, p)[np.newaxis, ...]
    return np.ascontiguousarray(canonical)


def _green_canonical_to_sc(g, spatial_shape):
    """Exact inverse of :func:`_green_sc_to_canonical`.

    ``(1, nmat, nvol, p, p) -> (p, p, Nx, Ny, Nz, nmat)``. ``nvol`` alone
    does not determine its own ``(Nx, Ny, Nz)`` factorization, so this
    direction takes ``spatial_shape`` explicitly (the same way every other
    kernel entry point -- ``bond_bubble_static``, ``dense_bubble`` -- does).
    ALWAYS returns a copy (``np.ascontiguousarray``).
    """
    g = np.asarray(g)
    if g.ndim != 5 or g.shape[0] != 1:
        raise ValueError(
            "_green_canonical_to_sc: expected canonical shape (1, nmat, "
            "nvol, p, p); got ndim={} shape={}".format(g.ndim, g.shape))
    _, nmat, nvol, p, p2 = g.shape
    if p2 != p:
        raise ValueError(
            "_green_canonical_to_sc: the two orbital axes must be equal; "
            "got shape {}".format(g.shape))
    Nx, Ny, Nz = spatial_shape
    if Nx * Ny * Nz != nvol:
        raise ValueError(
            "_green_canonical_to_sc: spatial_shape {} does not match "
            "nvol={} (shape {})".format(spatial_shape, nvol, g.shape))
    reshaped = g[0].reshape(nmat, Nx, Ny, Nz, p, p)
    out = reshaped.transpose(4, 5, 1, 2, 3, 0)          # (p,p,Nx,Ny,Nz,nmat)
    return np.ascontiguousarray(out)


class _BondGreen(NamedTuple):
    """Carrier separating the FULL Green function -- fed, unchanged, to
    EVERY downstream consumer of the bond path (``pair_weight``,
    ``make_bond_kernel_parts``, the Eliashberg convolutions) -- from the
    DEFLATED Green/tail pair, fed ONLY to the bubble
    (``bubble.bond_bubble_static``). See the unified-bubble-kernel spec's
    "Green data flow in sc.py's bond branch".

    Attributes
    ----------
    full_sc : ndarray
        sc.py layout ``(p, p, Nx, Ny, Nz, nmat)`` -- the bare Green
        function ``G(k, iw_n)`` built from the eigen-decomposition,
        independent of ``coeff_tail`` (the tail term the bubble subtracts
        in frequency space is added straight back for ``full``; only
        ``deflated_kw`` depends on the coefficient).
    deflated_kw : ndarray
        Canonical ``(1, nmat, nvol, p, p)``, for the bubble only.
    green0_tail : ndarray or None
        Canonical, paired with ``deflated_kw``; ``None`` when the tail
        correction is off (``coeff_tail_applied == 0``).
    coeff_tail_requested : float
        The resolved ``[eliashberg] bond_coeff_tail`` config value.
    coeff_tail_applied : float
        What the build actually used: equals ``coeff_tail_requested`` for
        an internal build, ``0.0`` for an external Green function (the
        tail correction is not well defined for a file whose asymptote is
        unknown).
    source : str
        ``"internal"`` (built from sc.py's own eigen-decomposition) or
        ``"external_npz"`` (``[eliashberg] bond_green``).
    """
    full_sc: np.ndarray
    deflated_kw: np.ndarray
    green0_tail: np.ndarray | None
    coeff_tail_requested: float
    coeff_tail_applied: float
    source: str


def _build_bond_green(eigenvalues, eigenvectors, mu, beta, nmat,
                      coeff_tail_requested, bond_green_path):
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

    Returns
    -------
    _BondGreen
    """
    if bond_green_path is not None:
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

def _symmetrise_interactions_k(inter_k):
    """Reduce each interaction to its physical symmetric coefficient.

    Thin delegation: the reduction and its full derivation -- the
    momentum-space form of the reversed-bond partner, the PairHop
    Hermitian rule, and the UHFk-relation adjudication note -- are
    single-sourced in :func:`hwave.solver.declarations.symmetrise_k`
    (#108). Idempotent, so it is safe that both the all-q and per-q
    builders apply it.
    """
    return symmetrise_k(inter_k)


def _accumulate_coeff(dst, coeff, value):
    """dst += coeff * value with IEEE-preserving forms.

    NumPy evaluates ``1.0 * (Inf+0j)`` as ``Inf+NaNj`` (a full complex
    multiply), while a direct add preserves the zero imaginary part, so
    the +-1 coefficients of the vertex table must use direct add/subtract
    to keep the pre-table behavior for non-finite couplings. Zero
    coefficients suppress the contribution entirely (0.0 * Inf is NaN).
    ``dst`` must be an ndarray or writeable view; the update is in place.
    """
    if coeff == 0.0:
        return
    if coeff == 1.0:
        dst += value
    elif coeff == -1.0:
        dst -= value
    else:
        dst += coeff * value


def _build_sc_matrices_all_q(inter_k, norb, Nx, Ny, Nz,
                             _presymmetrised=False):
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

    # _presymmetrised is set by the per-q wrapper, which has already
    # symmetrised on the FULL grid: the -q partner of an off-site entry is
    # unreachable from a single-point slice, so re-symmetrising here would
    # average with the wrong (same-q) partner and corrupt off-site input.
    if not _presymmetrised:
        inter_k = _symmetrise_interactions_k(inter_k)

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
    # Accumulate rather than assign: this element is also reached by Case 3
    # below (which now includes l1 == l3), where an inter-site same-orbital
    # CoulombInter contributes 2 V_aa(q) to the charge channel.
    mask1 = (l1f == l2f) & (l2f == l3f) & (l3f == l4f)
    if U_mat is not None and np.any(mask1):
        sU, cU = sc_coefficients("CoulombIntra", "diag")
        for i in np.where(mask1)[0]:
            _l = l1f[i]
            _accumulate_coeff(S_all[:, :, :, idx12[i], idx34[i]],
                              sU, U_mat[_l, _l])
            _accumulate_coeff(C_all[:, :, :, idx12[i], idx34[i]],
                              cU, U_mat[_l, _l])

    # Case 2: l1==l3, l2==l4, l1!=l2
    # Coefficients come from the single adjudicated table
    # (hwave.solver.vertex_table). Signs, slots and per-type factors were
    # established against exact diagonalization in #113 (Ising sign, the
    # previously missing Hund S term, Exchange moved here from the
    # pair-antidiagonal Case 4); the SU(2) Kanamori consistency check --
    # Hund + Exchange at equal J giving S(ab,ab) without J and
    # C(ab,ab) = -U' + 2J -- is pinned in the tests.
    mask2 = (l1f == l3f) & (l2f == l4f) & (l1f != l2f)
    cross_terms = [(Up_mat, "CoulombInter"), (I_mat, "Ising"),
                   (J_mat, "Hund"), (Jp_mat, "Exchange")]
    for i in np.where(mask2)[0]:
        s_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        c_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        _l1, _l2 = l1f[i], l2f[i]
        for mat, itype in cross_terms:
            if mat is not None:
                sco, cco = sc_coefficients(itype, "cross")
                _accumulate_coeff(s_q, sco, mat[_l1, _l2])
                _accumulate_coeff(c_q, cco, mat[_l1, _l2])
        S_all[:, :, :, idx12[i], idx34[i]] = s_q
        C_all[:, :, :, idx12[i], idx34[i]] = c_q

    # Case 3: l1==l2, l3==l4 -- INCLUDING l1 == l3, but for CoulombInter ONLY.
    # The l1 != l3 exclusion dropped the inter-site same-orbital CoulombInter
    # from the charge diagonal C[(a,a),(a,a)], which must be U_a + 2 V_aa(q)
    # (issue #95); the simple two-index formulation used by chi0q_mode="load"
    # builds exactly that (`Wc = U_k + 2 V_k`, _compute_vertices_simple).
    # Case 1 above writes U_a into the same element, so both accumulate.
    # Hund and Ising stay restricted to l1 != l3: an orbital has no Hund or
    # Ising coupling with itself, and letting a stray diagonal entry through
    # here would silently move S as well.
    mask3 = (l1f == l2f) & (l3f == l4f)
    for i in np.where(mask3)[0]:
        s_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        c_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        _l1, _l3 = l1f[i], l3f[i]
        if _l1 != _l3:
            for mat, itype in ((J_mat, "Hund"), (I_mat, "Ising")):
                if mat is not None:
                    sco, cco = sc_coefficients(itype, "density")
                    _accumulate_coeff(s_q, sco, mat[_l1, _l3])
                    _accumulate_coeff(c_q, cco, mat[_l1, _l3])
        if Up_mat is not None:
            sco, cco = sc_coefficients("CoulombInter", "density")
            _accumulate_coeff(s_q, sco, Up_mat[_l1, _l3])
            _accumulate_coeff(c_q, cco, Up_mat[_l1, _l3])
        S_all[:, :, :, idx12[i], idx34[i]] += s_q
        C_all[:, :, :, idx12[i], idx34[i]] += c_q

    # Case 4: l1==l4, l2==l3, l1!=l2
    mask4 = (l1f == l4f) & (l2f == l3f) & (l1f != l2f)
    for i in np.where(mask4)[0]:
        s_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        c_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        _l1, _l2 = l1f[i], l2f[i]
        # Exchange used to sit here; exact diagonalization (issue #113) puts
        # its vertex on the pair-DIAGONAL slot family (Case 2), and end to end
        # the antidiagonal placement produced the right magnitude at the wrong
        # slots in both channels. Only PairHop belongs here (#100/#102).
        if PH_mat is not None:
            sco, cco = sc_coefficients("PairHop", "antidiag")
            _accumulate_coeff(s_q, sco, PH_mat[_l1, _l2])
            _accumulate_coeff(c_q, cco, PH_mat[_l1, _l2])
        S_all[:, :, :, idx12[i], idx34[i]] = s_q
        C_all[:, :, :, idx12[i], idx34[i]] = c_q

    return S_all, C_all


def _build_sc_matrices(inter_k, norb, ix, iy, iz):
    """Spin (S) and charge (C) matrices at a single q-point.

    Delegates to :func:`_build_sc_matrices_all_q` and slices, so there is ONE
    implementation of the S/C content. The previous hand-maintained copy had
    already drifted from the all-q builder once; after the issue #113 vertex
    corrections a second parallel copy would be a liability.
    """
    # Symmetrise on the FULL grid first -- the same-operator partner of an
    # off-site entry lives at -q, so slicing before symmetrising would discard
    # it -- then slice the (small, norb^2 per point) interaction arrays down to
    # the requested q and delegate on a 1x1x1 grid, so the (nd^2 per point) S/C
    # matrices are never built for the whole grid. Indexing goes through
    # np.arange so numpy negative indices keep their usual meaning and
    # out-of-range indices raise instead of returning empty slices.
    inter_sym = _symmetrise_interactions_k(inter_k)
    inter_1 = {}
    for t, M in inter_sym.items():
        jx = np.arange(M.shape[2])[ix]
        jy = np.arange(M.shape[3])[iy]
        jz = np.arange(M.shape[4])[iz]
        inter_1[t] = np.ascontiguousarray(
            M[:, :, jx:jx+1, jy:jy+1, jz:jz+1])
    S_all, C_all = _build_sc_matrices_all_q(inter_1, norb, 1, 1, 1,
                                            _presymmetrised=True)
    return S_all[0, 0, 0], C_all[0, 0, 0]

def _declarations_partner_closed(interactions, Nx, Ny, Nz, norb):
    """True iff every CoulombIntra/CoulombInter declaration table is
    EXACTLY closed under the reversed-bond partner (R,a,b) <-> (-R,b,a).

    Decided on the RAW tables, before any exponential transform, with
    bitwise comparison: a both-ends declaration carries identical
    literals, and 0.5*(x + x) == x exactly in IEEE, so closure detection
    is exact -- while the k-space arrays of a closed table are only
    ALGEBRAICALLY symmetric (the +R and -R exponentials are evaluated
    independently and differ by roundoff), which is why the decision
    cannot be made after the transform (PR #129 round 5: re-averaging a
    closed configuration changed saved artifacts at the 1e-18 level).
    Non-finite input reads as not closed (explicitly, before any
    comparison), which fails toward symmetrising. The closure test
    compares the table DIRECTLY with its reversed/transposed partner --
    no mean is formed, so a closed coefficient above float64.max / 2
    cannot overflow into a false 'open' (PR #129 round 6 reproduced
    exactly that with the earlier mean-based comparison)."""
    for itype in ("CoulombIntra", "CoulombInter"):
        tbl = interactions.get(itype) if interactions else None
        if not tbl:
            continue
        arr = np.zeros((Nx, Ny, Nz, norb, norb), dtype=complex)
        for (irvec, orbvec), v in tbl.items():
            arr[(irvec[0] % Nx, irvec[1] % Ny, irvec[2] % Nz,
                 *orbvec)] += v
        if not np.all(np.isfinite(arr)):
            return False
        rev = np.transpose(reverse_fft_axes(arr, (0, 1, 2)),
                           (0, 1, 2, 4, 3))
        if not (np.array_equal(arr.real, rev.real)
                and np.array_equal(arr.imag, rev.imag)):
            return False
    return True


def _compute_vertices(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                      pairing_type="singlet", static_index=None,
                      declarations_closed=False):
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
                                        static_index=static_index,
                                        declarations_closed=declarations_closed)


def _compute_vertices_simple(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                             pairing_type="singlet", static_index=None,
                             declarations_closed=False):
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
    # Symmetrise FIRST (PR #129 round 3): this path read the raw tables,
    # so a one-sided off-site declaration entered as v e^{+iqR} (the
    # documented sign, #133) while the
    # ring and the general S/C route read the same Hamiltonian as
    # v cos(qR) -- measured drift 0.7 at q = pi/2 for V(R=+x) = 0.7.
    # SKIPPED when the caller proved the raw declarations partner-closed
    # (round 5): the reduction is only ALGEBRAICALLY idempotent after the
    # exponential transform -- re-averaging a closed configuration mixed
    # independently rounded +R/-R exponentials and changed saved
    # artifacts at the 1e-18 level, violating byte parity for
    # previously-working symmetric runs.
    if not declarations_closed:
        inter_k = _symmetrise_interactions_k(
            {k: inter_k[k] for k in ("CoulombIntra", "CoulombInter")
             if k in inter_k})
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
        # The guard runs for EVERY 2-index input, including norb == 1: the
        # rejection of Exchange/PairHop applies there too (the one-orbital
        # builder gives them zero weight, so accepting them silently omits
        # the interaction -- the norb shortcut below bypassed the helper
        # entirely, round 8). Only the partial-dressing warning is
        # norb-gated, inside the helper.
        _warn_reduced_flex_missing_components(
            inter_k, norb, Nx, Ny, Nz, source="a reduced 2-index chi0q",
            sc_matrices=(S_all, C_all))
        if norb == 1:
            chi0_static = chi0_2d.reshape(Nx, Ny, Nz, 1, 1)
        else:
            # A 2-index (reduced) chi0q is the density-density diagonal
            # of the susceptibility: chi0_2d[a, b] IS chi0_{(a,a),(b,b)} (the
            # matching interaction reduction is einsum('kaabb->kab', ...) in
            # RPA._inflate_chi0q_and_ham).  With the orbital-pair flat index
            # (l1,l2) -> l1*norb + l2 used by _build_sc_matrices_all_q, it
            # therefore belongs at the density-pair positions:
            #
            #     chi0_{(a,a),(b,b)} = chi0_2d[a, b],  everything else zero.
            #
            # The historical placement chi0_{(l1,l2),(l3,l2)} = chi0_2d[l1,l3]
            # (a delta_{l2,l4} scatter, i.e. kron(chi0_2d, I_norb)) read the
            # density-pair index as an ordinary orbital index: it dropped the
            # inter-orbital density coupling chi0_{(0,0),(1,1)} and scattered
            # chi0_2d onto pair indices the reduced scheme never computed.  This
            # is the same defect that _expand_flex_chi carried on the
            # chi0q_mode="flex" route.
            #
            # Off-density rows/columns stay exactly zero.  For the interaction
            # terms that put S/C weight there (CoulombInter, Hund, Ising,
            # Exchange, PairHop) those channels are then undressed -- an honest
            # reflection of what a reduced chi0q contains, rather than a
            # fabricated dressing.
            chi0_static = np.zeros((Nx, Ny, Nz, nd, nd), dtype=complex)
            dens = np.arange(norb) * norb + np.arange(norb)
            chi0_static[..., dens[:, None], dens[None, :]] = chi0_2d

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


#: Interaction terms whose Kuroki S/C matrices have entries OUTSIDE the
#: density-pair block, keyed by the _build_sc_matrices_all_q case that puts
#: them there.  Case 2 is S/C[(a,b),(a,b)] and case 4 is S/C[(a,b),(b,a)], both
#: with a != b; a reduced run never computes chi on those pair indices,
#: so the corresponding fluctuation dressing is missing from the vertex.
#: Types whose vertex is PARTIALLY represented by a reduced chi: their
#: density-slot content is dressed, their cross-slot content is not.
_REDUCED_FLEX_PARTIAL = ("CoulombInter", "Hund", "Ising")
#: Types with NO density-diagonal vertex content at all
#: (hwave.solver.vertex_table): with a reduced chi nothing of them is
#: dressed, and since the #120 policy no reduced FLEX or RPA run
#: can even be produced with them -- such input is stale or mismatched,
#: and is REJECTED rather than warned about.
_REDUCED_FLEX_REJECTED = ("Exchange", "PairHop")



def _off_density_sc_weight(inter_k, norb, Nx, Ny, Nz, sc_matrices=None):
    """Largest |S| or |C| on the blocks a reduced chi cannot dress.

    The reduced susceptibility is zero on every pair index (a,b) with a != b, so
    the vertex is undressed exactly where S or C is nonzero there --
    S/C[(a,b),(a,b)] (case 2) and S/C[(a,b),(b,a)] (case 4) of
    :func:`_build_sc_matrices_all_q`.

    This inspects the COMBINED matrices rather than testing each configured term
    on its own. Those blocks are sums: under the adjudicated slot table
    (#113) case 2 mixes CoulombInter, Ising, Hund AND Exchange, and case 4
    carries PairHop alone. Terms can cancel there (equal CoulombInter and
    Hund cancel case 2 exactly in both channels), and a per-term test would
    then announce missing dressing that does not exist.
    """
    # Reuse the caller's matrices when it already has them: at the full grid
    # these are O(Nq * norb^4) each, and building a second pair while the first
    # is live doubles the allocation for what is only a diagnostic.
    if sc_matrices is not None:
        S_all, C_all = sc_matrices
    else:
        S_all, C_all = _build_sc_matrices_all_q(inter_k, norb, Nx, Ny, Nz)
    nd = norb * norb
    density = np.zeros(nd, dtype=bool)
    density[np.arange(norb) * norb + np.arange(norb)] = True
    off = ~density
    if not off.any():
        return 0.0
    weight = 0.0
    for M in (S_all, C_all):
        # anything outside the density x density sub-block, over all q
        weight = max(weight,
                     float(np.max(np.abs(M[..., off, :]))),
                     float(np.max(np.abs(M[..., :, off]))))
    return weight


def _build_vertex_sc_matrices(convention, inter_k, norb, Nx, Ny, Nz):
    """Build the (S, C) interaction matrices for the given orbital convention.

    Single dispatcher shared by :func:`_compute_vertices_flex` and its
    callers, so a dynamic run can build the frequency-INDEPENDENT matrices
    once and pass them to every per-frequency contraction instead of
    rebuilding ~nmat identical full-grid pairs.
    """
    conv = str(convention).lower()
    if conv == "myo":
        from hwave.solver._sc_matrices_myo import build_sc_matrices_myo
        return build_sc_matrices_myo(inter_k, norb, Nx, Ny, Nz)
    if conv == "kuroki":
        return _build_sc_matrices_all_q(inter_k, norb, Nx, Ny, Nz)
    raise ValueError(
        "Unknown convention: '{}'. Use 'kuroki' or 'myo'.".format(convention))


def _reject_reduced_flex_unsupported(inter_k, convention="kuroki",
                                     source=None):
    """Reject Exchange/PairHop wherever a Kuroki-convention vertex is built.

    This is the cheap half of :func:`_warn_reduced_flex_missing_components`:
    a dict scan with no S/C construction, so it can (and must) run at EVERY
    Kuroki vertex boundary -- :func:`_compute_vertices_flex` itself, hence
    also ``eliashberg_dynamic._instantaneous_vertex`` and every
    per-frequency dynamic contraction -- making the rejection an enforced
    invariant rather than an ordering artifact of which builder the
    orchestrator happens to call first (round-9 review: the IR
    instantaneous route silently accepted both terms and returned a zero
    vertex when called directly).

    Exchange and PairHop carry NO density-diagonal vertex content
    (hwave.solver.vertex_table), so a Kuroki (reduced-origin) chi dresses
    none of it and the result would silently omit the interaction; the
    general (myo) convention represents them fully and is not restricted.
    """
    if str(convention).lower() != "kuroki":
        return
    rejected = [k for k in _REDUCED_FLEX_REJECTED
                if k in inter_k and np.any(np.asarray(inter_k[k]) != 0)]
    if rejected:
        raise ValueError(
            "the Eliashberg vertex cannot be built from {} together with "
            "{}: those interactions have no density-diagonal vertex "
            "content at all (hwave.solver.vertex_table), so reduced "
            "(density-only) data dresses none of it and the result would "
            "silently omit the interaction. (H-wave's own reduced "
            "runs cannot even be produced with these terms since the "
            "unified scheme policy.) Provide a general (four-index) "
            "susceptibility instead -- for a FLEX source, re-run with "
            "calc_scheme='general'.".format(
                source or ("a REDUCED (calc_scheme='reduced') FLEX "
                           "susceptibility"),
                ", ".join(rejected)))


def _warn_reduced_flex_missing_components(inter_k, norb, Nx, Ny, Nz,
                                          convention="kuroki",
                                          source=None, sc_matrices=None):
    """Warn when a reduced (kuroki) FLEX chi cannot support the interaction.

    Call this ONCE per run, from the place that is about to build the pairing
    vertex -- NOT from inside ``_compute_vertices_flex``. That function is
    invoked once per bosonic Matsubara frequency by the dynamic kernel (so the
    warning would repeat ``Nmat`` times, ~1000 in production runs), and again by
    ``eliashberg_dynamic._instantaneous_vertex`` with chi = 0, where the
    message would be doubly misleading. The Exchange/PairHop REJECTION, by
    contrast, is enforced inside ``_compute_vertices_flex`` on every call
    (:func:`_reject_reduced_flex_unsupported`) -- it is a cheap scan and
    must hold on every route, not just the orchestrated one.

    A ``calc_scheme="reduced"`` FLEX run stores only the
    density-density diagonal chi_{(a,a),(b,b)} of the susceptibility.  The
    Kuroki S/C matrices built from CoulombIntra alone live entirely on that
    density-pair block, so the reduced treatment is exact.  CoulombInter,
    Hund and Ising also populate the off-density block S/C[(a,b),(a,b)]
    with a != b -- and there chi is identically zero simply because the
    reduced run never computed it.  Those channels then keep only the bare
    0.5*(S+C) term with no fluctuation dressing: a silent approximation
    rather than a solver error, hence the WARNING.  Exchange and PairHop
    carry NO density-diagonal vertex content at all (vertex_table), so a
    reduced chi dresses nothing of them; combined with the #120 scheme
    policy (no reduced run can be produced with them) that
    combination is REJECTED rather than warned about.

    This is a genuine limitation of the stored data, not of the loader: it
    cannot be repaired on the Eliashberg side.  The universal remedy is a
    general (four-index) susceptibility; for a FLEX source that means
    re-running with ``calc_scheme="general"`` (which stores the full
    orbital-pair chi).
    """
    if str(convention).lower() != "kuroki":
        # Only the reduced route stores a density-only chi; the
        # general (myo) path carries the full orbital-pair susceptibility.
        return
    # The rejection must run BEFORE the norb == 1 shortcut: a one-orbital
    # Exchange/PairHop has zero S/C weight in this builder too, so the
    # interaction would be silently omitted rather than represented --
    # exactly what the rejection exists to prevent (round-7 review; the
    # bypass returned zero vertices for both terms on both routes).
    _reject_reduced_flex_unsupported(inter_k, convention, source)
    if norb <= 1:
        # norb == 1 has no off-density pair index, so no PARTIAL dressing
        # can be missing (the rejection above still applies).
        return
    configured = [k for k in _REDUCED_FLEX_PARTIAL if k in inter_k]
    if not configured:
        return
    # Ask the assembled S/C matrices, not the individual terms: the off-density
    # blocks are sums and the configured terms can cancel there (e.g. equal
    # CoulombInter and Hund cancel case 2 in both channels under the
    # adjudicated table). Only a nonzero combined block means dressing is
    # actually missing.
    if _off_density_sc_weight(inter_k, norb, Nx, Ny, Nz,
                              sc_matrices) == 0.0:
        return
    # Decision above is on the assembled matrices; attribution here is per term,
    # so an inert term that happens to be configured is not named as a cause.
    # Attribute through the same assembled-and-symmetrised path as the
    # decision: a raw declaration that symmetrises to zero (#112) must not
    # be named as a cause.
    missing = [k for k in configured
               if _off_density_sc_weight({k: inter_k[k]}, norb,
                                         Nx, Ny, Nz) != 0.0] or configured
    logger.warning(
        "The Eliashberg vertex is being built from %s, which carries only the "
        "density-density components chi_{(a,a),(b,b)}, together with "
        "inter-orbital interaction(s) %s. Together those terms give the S/C "
        "matrices nonzero "
        "off-density components S/C[(a,b),(a,b)] and S/C[(a,b),(b,a)] (a != b) "
        "for which the reduced run computed no susceptibility, so those "
        "channels enter the pairing vertex undressed (bare 0.5*(S+C) only) and "
        "the resulting lambda is an approximation. calc_scheme='general' stores "
        "the full orbital-pair chi and removes this gap; note its off-site "
        "support is limited to same-orbital CoulombInter without sublattice "
        "folding -- for a model with off-site INTER-orbital interactions "
        "there is currently no FLEX-dressed vertex without this "
        "approximation.",
        source or ("a REDUCED (calc_scheme='reduced') FLEX "
                   "susceptibility"),
        ", ".join(missing))


def _compute_vertices_flex(chis, chic, inter_k, norb, Nx, Ny, Nz,
                           pairing_type="singlet", convention="kuroki",
                           sc_matrices=None):
    """Compute pairing vertex from pre-computed FLEX susceptibilities.

    Uses the same formula as _compute_vertices_general but takes
    chi_s and chi_c directly instead of computing them from chi0q via RPA.
    This allows using dressed susceptibilities from the FLEX solver.

    ``convention`` records which FLEX path produced the susceptibilities:
    "kuroki" (default; reduced path) or "myo" (general full-vertex path). The
    two S/C builders historically differed in the charge ``C(ab,ab)`` element
    (``-U'+2J`` vs ``-U'+J``); the exact-diagonalization adjudication of the
    per-type vertex content (issue #113: Hund ``+J`` and Exchange ``+J'``
    there) made them IDENTICAL, so the flag no longer selects different
    matrices -- it remains a provenance / layout discriminator, and legacy
    files are guarded by ``sc_vertex_version`` instead. ``chis``/``chic`` are
    expected in the native [a,c,b,d] orbital-pair order regardless
    (issue #78).

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
    sc_matrices : tuple of ndarray, optional
        Precomputed ``(S_all, C_all)`` for this convention
        (:func:`_build_vertex_sc_matrices`). The matrices are
        frequency-independent, so a dynamic run passes one pair to every
        per-frequency call instead of rebuilding ~nmat identical full-grid
        pairs (round-9 review).

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

    # Validate the convention BEFORE the cache branch: with precomputed
    # sc_matrices the dispatcher below is skipped, and an unknown tag would
    # otherwise be silently treated as "not kuroki" = unrestricted
    # (round-10 review reproduced convention="invalid" returning a vertex).
    conv = str(convention).lower()
    if conv not in ("kuroki", "myo"):
        raise ValueError(
            "Unknown convention: '{}'. Use 'kuroki' or 'myo'.".format(
                convention))

    # Enforced at THIS boundary (not only in the orchestrator): every Kuroki
    # vertex construction -- including chi = 0 via _instantaneous_vertex --
    # must reject Exchange/PairHop, or the interaction is silently omitted.
    _reject_reduced_flex_unsupported(inter_k, conv)

    if sc_matrices is not None:
        S_all, C_all = sc_matrices
    else:
        S_all, C_all = _build_vertex_sc_matrices(conv, inter_k,
                                                 norb, Nx, Ny, Nz)

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
        green_w = _load_flex_green(input_dict, norb, Nx, Ny, Nz)
    else:
        (chi_s_raw, chi_c_raw, chi_convention, ir_meta,
         legacy_tags) = _read_flex_chi_raw(
            input_dict, allow_ir=True, interactions=interactions)
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
        if affected or orientation or resym:
            versions = {}
            for data, path in ((data_s, chi_s_path), (data_c, chi_c_path)):
                if "sc_vertex_version" not in data:
                    reasons = []
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
                try:
                    arr = np.asarray(data["sc_vertex_version"])
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
                        "sc_vertex_version field ({!r}).".format(
                            path, data["sc_vertex_version"]))
            if len(set(versions.values())) != 1:
                raise ValueError(
                    "FLEX chi_s and chi_c carry different sc_vertex_version "
                    "values ({}); they must come from the same run.".format(
                        versions))
            ver = next(iter(versions.values()))
            if ver != 2:
                raise ValueError(
                    "FLEX susceptibility files carry sc_vertex_version = {} "
                    "but this code supports version 2 (the #113 exact-"
                    "diagonalization vertex content); regenerate the "
                    "susceptibilities with the current code.".format(ver))
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


#: Relative size, against the kept spin-up block, below which discarded spin
#: content is attributed to round-off rather than physics.
#:
#: A few hundred ulp of double precision. This is a NARROW MARGIN, not a
#: validated error bound: the error a linear solve can accumulate depends on
#: dimension, conditioning, backend and operation order, and near a
#: susceptibility pole a legitimate paramagnetic producer could exceed it. The
#: claim it rests on is only the measured one -- every producer exercised here
#: (CPU, uniform and IR axes, static and dynamic, plus production multi-orbital
#: output) is bit-exact, so the margin is never used in practice and exists so
#: that a backend whose solve is not bit-symmetric degrades to a warning rather
#: than aborting.
#:
#: It is NOT a claim about how weak a physical field can be. An earlier draft
#: used 1e-8, justified by a measured Zeeman case at ratio ~1.6 -- which sets no
#: lower bound at all, and would have relabelled a genuinely weak field as
#: round-off. If a supported backend is ever found to exceed this margin
#: legitimately, widen it deliberately with that measurement in hand rather than
#: treating the constant as already covering it.
_SPIN_DISCARD_ROUNDOFF_RATIO = 256 * np.finfo(float).eps


def _check_spin_block_discarded(chi_raw, norb, convention, label="chi",
                                legacy_up_block_only=False):
    """Reject, or at minimum report, spin content the embedding would drop.

    The reduced spin-orbital index is spin-block ordered ``s*norb + a``, and the
    embedding below keeps ONLY the spin-up block ``[:norb, :norb]``. Everything
    else -- the down block and both cross-spin blocks -- is dropped, because
    nothing downstream carries spin: the Kuroki S/C matrices and the
    singlet/triplet vertex formulas are norb^2-sized and paramagnetic.

    Dropping it is lossless exactly when it is redundant: the down block equals
    the up block and the cross blocks are zero. That holds bit-for-bit for a
    paramagnetic run -- the inflation writes the same array into both spin
    blocks and never touches the cross blocks
    (``FLEX._inflate_chi0q_and_ham``), and the channel vertices are spin-block
    diagonal, so the RPA solve preserves both properties.

    When it does not hold the discarded data is real and the eigenvalue is not a
    controlled approximation of anything, so this RAISES rather than warning:
    a spin-polarized model is outside what this formulation can express, and
    returning a number for it is worse than refusing.

    The test is on the DISCARDED DATA, not on the run's spin mode. Inferring
    "spin-polarized" from unequal diagonal blocks would be wrong in both
    directions: a spinful/spin-orbit run can have unequal blocks while still
    being time-reversal symmetric, and it can equally have equal diagonal blocks
    with nonzero cross-spin blocks -- discarded just the same. The stored npz
    records no ``spin_mode``, so the data is the only signal available.

    Two severities, because a false positive costs very different amounts:
    anything nonzero is reported, but only content above
    ``_SPIN_DISCARD_ROUNDOFF_RATIO`` -- a few hundred ulp of double precision --
    aborts. Exact bit-equality was confirmed on the uniform and IR axes, static
    and dynamic, and on production multi-orbital output, so the window is there
    purely so a backend whose solve is not bit-symmetric (a GPU batched solve,
    say) degrades to a warning instead of killing a legitimate paramagnetic run.
    It is deliberately anchored to machine precision and not to any assumption
    about how weak a physical field can be.
    """
    if str(convention).lower() != "kuroki":
        # Only the reduced layout has spin blocks to compare.
        return
    ns = 2
    chi = np.asarray(chi_raw)
    if norb <= 0:
        return
    if chi.ndim != 4 or chi.shape[-1] != chi.shape[-2]:
        # Rank AND squareness, before any dimension-based early return. Checking
        # only the last axis let a rectangular array such as (nfreq, nvol, 2, 8)
        # slip past -- its last dimension is not norb*ns, so the check returned,
        # and _expand_flex_chi then flattened 2*8 = 16, took sqrt, and
        # reinterpreted the rectangle as a 4x4 susceptibility.
        raise ValueError(
            "spin-block check expects a (nfreq, nvol, nd, nd) susceptibility "
            "with a square trailing matrix, got shape {}.".format(chi.shape))
    if chi.shape[-1] != norb * ns:
        # Not the reduced spin-orbital layout; _expand_flex_chi gives a better
        # diagnostic for this (it can recognise the spin-orbital-mode case).
        return
    up = chi[..., :norb, :norb]
    if up.size == 0:
        return
    down = chi[..., norb:, norb:]
    cross_ud = chi[..., :norb, norb:]
    cross_du = chi[..., norb:, :norb]

    # Accumulate per frequency rather than over the whole axis. Keeping the
    # whole-axis POLICY (see the loaders) does not require whole-axis
    # TEMPORARIES: abs(up), up-down and abs(up-down) at full size are each of
    # the order of the stored array, which on a production file would add
    # multi-GB allocations and could fail a valid file on memory alone.
    nfreq = chi.shape[0]

    # Per-frequency RATIOS, not a global difference over a global scale. The
    # ratio decides acceptance, and a single global pair lets a large redundant
    # frequency mask a small completely non-redundant one: with frequency 0
    # redundant at scale 1e14 and frequency 1 holding up=1, down=0, the true
    # mismatch is 100% while the global ratio is 1e-14 -- under the threshold.
    #
    # Non-finite values are rejected outright here rather than folded into the
    # maxima: max(0.0, nan) is 0.0 in Python, so a NaN in a discarded block
    # would leave the running maximum untouched and the array would be accepted
    # as exactly redundant.
    scale = d_down = d_cross = 0.0
    worst_ratio = 0.0
    worst_at = (0, 0.0, 0.0, 0.0)
    any_down = any_cross = False
    for w in range(nfreq):
        u, dn = up[w], down[w]
        cu, cd = cross_ud[w], cross_du[w]
        for name, blk in (("up", u), ("down", dn),
                          ("up-down cross", cu), ("down-up cross", cd)):
            if blk.size and not np.all(np.isfinite(blk)):
                raise ValueError(
                    "The reduced FLEX susceptibility {} has non-finite values "
                    "(NaN or Inf) in its {} block at frequency index {}. NaN "
                    "compares false against everything, so without this check "
                    "it would leave the redundancy maxima untouched and the "
                    "array would be accepted as exactly "
                    "redundant.".format(label, name, w))
        w_scale = float(np.max(np.abs(u)))
        w_down = float(np.max(np.abs(u - dn)))
        w_cross = max(float(np.max(np.abs(cu))), float(np.max(np.abs(cd))))
        scale = max(scale, w_scale)
        d_down = max(d_down, w_down)
        d_cross = max(d_cross, w_cross)
        w_ratio = max(w_down, w_cross) / max(w_scale, np.finfo(float).tiny)
        if w_ratio > worst_ratio:
            worst_ratio = w_ratio
            worst_at = (w, w_scale, w_down, w_cross)
        any_down = any_down or bool(np.any(dn))
        any_cross = any_cross or bool(np.any(cu)) or bool(np.any(cd))

    if d_down == 0.0 and d_cross == 0.0:
        return

    # An all-zero down/cross block is NOT a reliable marker of the legacy
    # single-populated representation: a saturated or projected spin sector, or
    # externally produced spin-resolved data, can look exactly the same. Absent
    # data and fully polarized data are indistinguishable from the values alone,
    # so the file has to SAY which it is. Accept it only when it does.
    if scale > 0.0 and not any_down and not any_cross:
        if legacy_up_block_only:
            # Name the actual source of the authorization: a file tag and a user
            # override carry different weight, and reporting one as the other
            # would misrepresent who vouched for the data.
            via = ("the [eliashberg] accept_up_block_only setting"
                   if legacy_up_block_only == "config_override"
                   else "the file's chi_spin_blocks='up_only' tag")
            logger.warning(
                "The reduced FLEX susceptibility %s has an all-zero down-spin "
                "block and all-zero cross-spin blocks. Per %s this is read as "
                "the layout in which only the block the Eliashberg step "
                "consumes was ever populated -- not as a spin-polarized run -- "
                "and the up-spin block is used.", label, via)
            return
        raise ValueError(
            "The reduced FLEX susceptibility {} has a nonzero up-spin block "
            "with an identically zero down-spin block and zero cross-spin "
            "blocks. That is either a file in which only the consumed block was "
            "ever populated, or a fully spin-polarized/projected susceptibility "
            "-- the values alone cannot tell them apart, and the second is not "
            "something this paramagnetic formulation can use. If you know it is "
            "the former, say so: set [eliashberg] accept_up_block_only = true, "
            "or tag the file itself with chi_spin_blocks='up_only'. Files "
            "written before this check existed carry no tag, so the "
            "configuration flag is the route for them.".format(label))

    # Describe the frequency that decided it: acceptance is a per-frequency
    # ratio, so quoting global maxima could pair a small tail's mismatch with a
    # huge unrelated scale from elsewhere on the axis.
    w_at, scale, d_down, d_cross = worst_at
    parts = []
    if d_down != 0.0:
        parts.append("the down-spin block differs from the up-spin block "
                     "by {:.3e}".format(d_down))
    if d_cross != 0.0:
        parts.append("the cross-spin blocks are nonzero "
                     "(max {:.3e})".format(d_cross))
    detail = " and ".join(parts)
    ratio = worst_ratio

    if ratio <= _SPIN_DISCARD_ROUNDOFF_RATIO:
        logger.warning(
            "The reduced FLEX susceptibility %s is not exactly spin-redundant: "
            "%s, against an up-block scale of %.3e at frequency index %d. "
            "That is %.1e of the kept block there -- within round-off of "
            "double precision, so it is treated as noise and the calculation "
            "continues using the up-spin block. Please report it: every tested "
            "path produces "
            "bit-identical spin blocks.", label, detail, scale, w_at, ratio)
        return

    raise ValueError(
        "The reduced FLEX susceptibility {} is SPIN-RESOLVED: {}, against an "
        "up-block scale of {:.3e} at frequency index {}. The Eliashberg "
        "pairing vertex here is "
        "paramagnetic -- the Kuroki S/C matrices carry no spin index and the "
        "singlet/triplet decomposition assumes spin-rotational symmetry -- so "
        "only the up-spin block could be used and the rest would be discarded. "
        "That is not a controlled approximation of the spin-resolved problem, "
        "so it is refused rather than silently returning a number. If you meant "
        "to run a paramagnetic calculation and a field such as coeff_extern is "
        "the only source of the splitting, drop it and re-run FLEX. If the spin "
        "structure is intrinsic to your model, this solver cannot describe it: "
        "a spin-resolved treatment needs an S_z-resolved vertex and the "
        "transverse susceptibility chi^(+-), which FLEX does not compute at all "
        "(calc_type='ring+ladder' is rejected on every FLEX scheme).".format(
            label, detail, scale, w_at))


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

    green_dressed = _load_flex_green(input_dict, norb, Nx, Ny, Nz)
    return chis, chic, green_dressed, chi_convention


# ---------------------------------------------------------------------------
# G2 and Eliashberg kernel
# ---------------------------------------------------------------------------

def _finite_positive_float64(value):
    """Convert a declared real-numeric scalar to float64 and validate the
    CONVERTED value; return it, or None if the value is structurally or
    numerically unacceptable.

    Conversion must happen BEFORE the finite/positive checks (round-9
    review): on platforms where np.longdouble is wider than binary64, a
    finite extended-precision value can pass an elementwise isfinite test
    and then narrow to inf (or a positive one to 0.0) in the float64 the
    computation actually uses -- recreating the fail-open comparisons these
    gates exist to prevent. Validating the narrowed float64 closes that.
    """
    if isinstance(value, (list, tuple)):
        # the contract is a SCALAR: [2.0] or [[2.0]] silently passing an
        # asarray-based gate would bless container-typed configs
        return None
    arr = np.asarray(value)
    # ndim == 0, not size == 1: np.array([2.0]) and np.array([[2.0]]) are
    # containers too, and NPZ metadata arrives as ndarray already
    if arr.ndim != 0 or arr.dtype.kind not in "iuf":
        return None
    v = float(arr.item())
    if not np.isfinite(v) or not v > 0:
        return None
    return v


# Numerically safe upper bound for beta in the pair bubble: the batched
# GEMM forms ~beta^2 intermediates before the division by beta, so a
# physically meaningless beta (T < 1e-75) would overflow them to NaN.
_G2_BETA_MAX = 1.0e75


def _coerce_run_beta(T):
    """Validate the run temperature and return beta = 1/T.

    The provenance gate compares the FILE's beta against the run's, and a
    NaN/Inf run temperature makes every tolerance comparison False -- i.e.
    an invalid run would slip through as 'matching' and calc_eliashberg
    would then build the whole Matsubara grid from beta = NaN (round-6
    review). Same real-numeric/finite/positive contract as the file-side
    beta gate.
    """
    t = _finite_positive_float64(T)
    if t is None:
        raise ValueError(
            "mode.param T must be a single finite positive real scalar, "
            "got {!r}.".format(T))
    beta = 1.0 / t
    # A subnormal T passes the checks above but its reciprocal overflows
    # to inf, recreating exactly the fail-open comparison this helper
    # exists to prevent (round-7 review) -- gate the RESULT too.
    if not np.isfinite(beta) or not beta > 0:
        raise ValueError(
            "mode.param T = {!r} yields a non-finite beta = 1/T; use a "
            "physically meaningful temperature.".format(T))
    return beta


def _coerce_config_bool(value, name):
    """Strictly parse a physics-relevant boolean config switch.

    backend.as_bool reads any unrecognized string ("ture", "garbage") as
    False and any nonzero integer as True -- for a switch that changes the
    physics (or gates a rejection), a spelling error must fail loudly
    instead of silently flipping the behavior. Accepted: real booleans,
    integers 0/1, and the strings true/false/yes/no/on/off/0/1 (case- and
    whitespace-insensitive).
    """
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, np.integer)) and value in (0, 1):
        return bool(value)
    if isinstance(value, str):
        v = value.strip().lower()
        if v in ("true", "yes", "on", "1"):
            return True
        if v in ("false", "no", "off", "0"):
            return False
    raise ValueError(
        "{} = {!r} is not a recognized boolean; use "
        "true/false (or yes/no, on/off, 0/1).".format(name, value))


def _coerce_g2_tail(value):
    """Strictly parse the [eliashberg] g2_tail switch (round-3 review)."""
    return _coerce_config_bool(value, "[eliashberg] g2_tail")


def _resolve_spin_orbital_flag(input_dict):
    """Resolve [mode] enable_spin_orbital to a real boolean, strictly.

    One resolver for every consumer in this module: the raw value used to
    be forwarded as-is, so a string "false" was truthy at the chi0q
    convention check and internally inconsistent inside RPA (its orbital
    counting uses truthiness while the transfer remap tests == True).
    Case-insensitive key lookup (config layers disagree on case handling;
    PR #128 sweep); unrecognized values raise.
    """
    from requests.structures import CaseInsensitiveDict
    return _coerce_config_bool(
        CaseInsensitiveDict(input_dict.get("mode", {})).get(
            "enable_spin_orbital", False),
        "[mode] enable_spin_orbital")


# Threshold for the asymptotic-regime diagnostic below: at the window edge
# |i wn G(k, i wn) - I| measures how far the Green function still is from
# its leading 1/(i wn) tail. 0.5 flags the regime where the neglected
# scales are at least comparable to the largest retained frequency.
_G2_TAIL_EDGE_DEV_THRESHOLD = 0.5


def _warn_if_g2_tail_outside_asymptotic_regime(green_kw, beta):
    """Warn when the tail correction's asymptotic premise looks violated.

    The correction subtracts the model I/(i wn) x -I/(i wn) inside the
    window and adds its exact full sum, which IMPROVES the result only when
    the largest retained frequency exceeds the relevant energy/self-energy
    scales -- otherwise the model does not describe the summand anywhere in
    the window and the added beta/4-complement can overshoot (measured:
    H = diag(1, -1), beta = 10, Nmat = 2 makes the corrected cross-orbital
    channel WORSE than the bare sum while staying PSD, so the positivity
    diagnostic alone cannot catch it). The same check flags a loaded Green
    function whose leading coefficient is not the identity (e.g. a scaled
    or truncation-damaged file): for those the fixed unit-coefficient
    correction is wrong by construction.

    Returns the measured edge deviation max_k |i wn G(k, i wn) - I| over
    both window endpoints.
    """
    norb = green_kw.shape[0]
    nmat = green_kw.shape[-1]
    if nmat < 1 or green_kw.size == 0:
        # Nothing to measure on an empty axis (or norb = 0); return
        # neutrally so the public ordering (this diagnostic, then
        # _calc_g2) surfaces the actionable even-grid ValueError instead
        # of an IndexError/max-of-empty here.
        return 0.0
    eye = np.eye(norb).reshape(norb, norb, 1, 1, 1)
    dev = 0.0
    for idx in (0, nmat - 1):
        wn = (2.0 * idx + 1.0 - nmat) * np.pi / beta
        dev = max(dev, float(np.abs(1j * wn * green_kw[..., idx]
                                    - eye).max()))
    if dev > _G2_TAIL_EDGE_DEV_THRESHOLD:
        logger.warning(
            "g2_tail: at the largest retained Matsubara frequency the Green "
            "function deviates from its asymptotic tail I/(i wn) by %.2f "
            "(threshold %.2f). The tail correction is asymptotic and can "
            "OVERSHOOT in this regime; the result may be less accurate than "
            "the bare sum, and passing the positivity check does not certify "
            "accuracy. Increase Nmat until the window edge exceeds the "
            "relevant energy scales, or check the loaded Green function's "
            "normalization/provenance.",
            dev, _G2_TAIL_EDGE_DEV_THRESHOLD)
    return dev


def _calc_g2(green_kw, beta, tail=True):
    """Calculate G2 = T * sum_n G(k, wn) G(-k+q, -wn).

    The temperature factor T = 1/beta is included so that the
    Eliashberg kernel correctly computes:
        lambda * sigma(k) = -(T/N_L) sum_{k',n'} P(k-k') G(k') G(-k') sigma(k')

    With ``tail=True`` (the default) the truncated Matsubara sum gets the
    analytic high-frequency tail correction (issue #86). The summand's exact
    leading tail is delta_ij delta_lm / wn^2 -- the 1/(i wn) coefficient of
    G is the identity by completeness of the eigenbasis, and a self-energy
    only enters at the next order, so this holds for dressed (FLEX) Green
    functions too, with no free coefficient (unlike chi0q's user-tunable
    ``coeff_tail``). Subtracting that model inside the window and adding its
    exact full sum (1/beta) sum_{n in Z} 1/wn^2 = beta/4 amounts to adding

        c = beta/4 - (1/beta) sum_{n in window} 1/wn^2   (> 0)

    to G2[i,i,l,l], i.e. c times the identity on the (il) gap space. This
    reduces the truncation error from O(1/Nmat) to the next order and in
    practice restores the positive semi-definiteness of G2 that makes the
    static Eliashberg kernel similar to a Hermitian matrix (real spectrum);
    the higher-order remainder can still leave a small exact eigenvalue
    negative, which _warn_if_g2_indefinite reports. The bare truncated sum
    can be slightly indefinite, which then injects spurious imaginary parts
    into the reported eigenvalues at small Nmat.

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
    """
    norb = green_kw.shape[0]
    Nx, Ny, Nz, nmat = (green_kw.shape[2], green_kw.shape[3],
                        green_kw.shape[4], green_kw.shape[5])
    nvol = Nx * Ny * Nz

    # beta gate (round-6 review): a non-finite or non-positive beta poisons
    # every frequency below, and a finite but astronomically large one
    # overflows the batched GEMM's ~beta^2 intermediates to NaN before the
    # division by beta -- measured: beta = 1e200 turns the whole pair
    # bubble into NaN while G itself is finite. beta values beyond
    # _G2_BETA_MAX have no physical meaning (T < 1e-75), so reject rather
    # than restructure the arithmetic (the tail=False stream must stay
    # bit-identical to the pre-correction implementation for valid input).
    beta_f = _finite_positive_float64(beta)
    if beta_f is None:
        raise ValueError(
            "beta must be a single finite positive real scalar, got "
            "{!r}.".format(beta))
    beta = beta_f
    if beta > _G2_BETA_MAX:
        raise ValueError(
            "beta = {!r} exceeds the numerically safe range (<= {!r}): "
            "the pair-bubble accumulation would overflow to NaN. Check "
            "the temperature.".format(beta, _G2_BETA_MAX))
    # Grid gate BEFORE any O(norb^2 nvol nmat) work (round-6 review: the
    # guard previously ran after the GEMM, which emitted overflow warnings
    # from data the call then rejected). The centered grid
    # wn = (2n + 1 - nmat) * pi / beta is only fermionic for EVEN nmat: an
    # odd nmat puts wn = 0 in the window, and nmat <= 0 would return the
    # bare analytic shift beta/4 with no Green-function samples at all.
    if tail and (nmat <= 0 or nmat % 2 != 0):
        raise ValueError(
            "the Matsubara tail correction (g2_tail) requires an even, "
            "positive number of frequencies on the centered fermionic "
            "grid; the Green function has nmat = {}. Fix the frequency "
            "axis of the input Green function (or Nmat), or set "
            "[eliashberg] g2_tail = false.".format(nmat))
    # Promote reduced-precision (e.g. complex64 file) input once so the
    # accumulation runs in complex128; a no-op for complex128 input, which
    # keeps the tail=False bit-parity with the pre-correction code.
    if green_kw.dtype != np.complex128:
        green_kw = np.ascontiguousarray(green_kw, dtype=np.complex128)

    # G(-k, -wn): centered-Matsubara flip on the frequency axis, then the
    # shared FFT-grid map k -> -k on the spatial axes.
    green_kw_inv = reverse_fft_axes(green_kw[..., ::-1], (2, 3, 4))
    # einsum("ijpqsk, lmpqsk -> ijlmpqs") sums over k (nmat dimension)
    # Reshape to use tensordot for BLAS: contract last axis (nmat) after
    # merging spatial dims
    # A: (norb, norb, nvol, nmat) -> (norb^2, nvol*nmat) -- but need per-site sum
    # Better: reshape to (norb*norb, nvol, nmat) then per-site outer product
    A = green_kw.reshape(norb * norb, nvol, nmat)       # (ij, site, k)
    B = green_kw_inv.reshape(norb * norb, nvol, nmat)   # (lm, site, k)
    # G2[ij, lm, site] = sum_k A[ij, site, k] * B[lm, site, k]
    # Per-site batched GEMM over the Matsubara axis n. numpy.einsum does NOT
    # lower 'isn,jsn->ijs' to BLAS GEMM; move site to the batch axis and matmul:
    #   As[s] @ Bs[s].T -> [i, j] per site; moveaxis(...,0,2) -> (i, j, s).
    As = np.moveaxis(A, 1, 0)             # (nvol, norb^2, nmat)
    Bs = np.moveaxis(B, 1, 0)             # (nvol, norb^2, nmat)
    G2 = np.moveaxis(As @ Bs.transpose(0, 2, 1), 0, 2)  # (norb^2, norb^2, nvol)
    G2 = G2.reshape(norb, norb, norb, norb, Nx, Ny, Nz)
    G2 = G2 / beta
    if tail:
        # Dimensionless form of the window complement: with
        # wn = r * pi / beta for the odd integers r = 2n + 1 - nmat,
        # beta/4 - (1/beta) sum 1/wn^2 = beta * (1/4 - (1/pi^2) sum 1/r^2).
        # Working in r keeps every intermediate O(1) regardless of beta
        # (the wn-based form underflows wn^2 to 0 for extreme beta).
        r = 2.0 * np.arange(nmat) + 1.0 - nmat
        coeff = beta * (0.25 - np.sum(1.0 / r**2) / np.pi**2)
        di = np.arange(norb)
        # G2[i, i, l, l] += coeff for every (i, l), all k
        G2[di[:, None], di[:, None], di[None, :], di[None, :]] += coeff
    # Fail-fast on a non-finite bubble (round-6 review): downstream
    # diagnostics may legitimately SKIP on size, and the solvers would then
    # propagate NaN into saved eigenvalues silently. Scan in bounded chunks
    # (round-9 review): a whole-tensor np.isfinite allocates a G2-sized
    # Boolean mask -- ~6% of the tensor again, enough to tip a large run
    # that otherwise fits into OOM.
    # ravel(order="K") flattens in MEMORY order and shares storage with
    # the (non-contiguous, ufunc-'K'-layout) G2; reshape(-1) would copy
    # the whole tensor here -- worse than the mask it replaced.
    flat = G2.ravel(order="K")
    step = 1 << 20
    for start in range(0, flat.size, step):
        if not np.all(np.isfinite(flat[start:start + step])):
            raise ValueError(
                "G2 (pair bubble) contains non-finite values; check the "
                "input Green function and beta.")
    return G2


# _warn_if_g2_indefinite skip thresholds (module-level so tests can patch
# them): eigvalsh work scales as (norb^2)^3 per k-point, i.e. norb^6 * nvol
# flops up to a constant. The peak memory is ~three full complex copies of
# the (nvol, norb^2, norb^2) view (the transpose/reshape materializes one,
# the residual |M - M^dag| evaluation another, the Hermitized matrix a
# third; tracemalloc-measured peak ~1.5x the two-copy figure), each
# 16 * norb^4 * nvol bytes.
_G2_CHECK_MAX_WORK = 20_000_000_000
_G2_CHECK_MAX_BYTES = 2_000_000_000


def _warn_if_g2_indefinite(G2, norb, tail_enabled):
    """Warn when the pair bubble G2 has a significantly negative eigenvalue.

    G2 viewed as a matrix on the (i,l) gap space, M[(i,l),(j,m)] = G2[i,j,l,m],
    is Hermitian positive semi-definite in the exact Matsubara sum; that is
    what makes the static kernel K = -Gamma W similar to a Hermitian matrix
    and its spectrum real. Truncation breaks positivity (issue #86), and
    users then see complex eigenvalues that look like a broken symmetry.
    Point them at Nmat / g2_tail instead of leaving that misread silent.

    A significantly non-Hermitian G2 gets its own warning first (a loaded
    Green function can be malformed, non-causal, or on the wrong grid;
    silently Hermitizing would hide exactly that defect), and the PSD test
    then runs on the Hermitian part. The check diagonalizes an
    (norb^2 x norb^2) block per k-point, so it is skipped (with a log line)
    when the eigensolver work estimate norb^6 * nvol or the temporaries'
    size would rival the solve itself.
    """
    nvol = int(np.prod(G2.shape[4:]))
    work = norb**6 * nvol
    nbytes = 3 * 16 * norb**4 * nvol
    if work > _G2_CHECK_MAX_WORK or nbytes > _G2_CHECK_MAX_BYTES:
        logger.info("G2 positivity check skipped (work estimate norb^6 * "
                    "nvol = %d, temporaries = %d bytes).", work, nbytes)
        return None
    M = G2.reshape(norb, norb, norb, norb, nvol)
    M = M.transpose(4, 0, 2, 1, 3).reshape(nvol, norb * norb, norb * norb)
    scale = float(np.abs(M).max()) if M.size else 0.0
    herm_residual = (float(np.abs(M - np.conj(M.transpose(0, 2, 1))).max())
                     if M.size else 0.0)
    if herm_residual > 1.0e-8 * max(scale, 1.0e-300):
        logger.warning(
            "G2 (pair bubble) is significantly non-Hermitian on the gap "
            "space: max |M - M^dag| = %.3e against max |M| = %.3e. The "
            "exact pair bubble is Hermitian, so check the input Green "
            "function's provenance (grid, temperature, causality); the "
            "positivity diagnostic below uses only the Hermitian part.",
            herm_residual, scale)
    M = 0.5 * (M + np.conj(M.transpose(0, 2, 1)))
    ev = np.linalg.eigvalsh(M)
    min_ev, max_ev = ev.min(), ev.max()
    if min_ev < -1.0e-6 * max(max_ev, 1.0e-300):
        hint = ("increase Nmat" if tail_enabled
                else "increase Nmat or re-enable the tail correction "
                     "([eliashberg] g2_tail = true)")
        logger.warning(
            "G2 (pair bubble) is not positive semi-definite: min eigenvalue "
            "= %.3e (max = %.3e). The most likely cause is Matsubara "
            "truncation error rather than a broken symmetry, and complex "
            "Eliashberg eigenvalues reported below are then spurious at "
            "this level -- %s. If the Green function was loaded from file, "
            "also check that its grid and temperature match this run.",
            min_ev, max_ev, hint)
    return min_ev


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
    rev = reverse_fft_axes(gap, (2, 3, 4))  # k -> -k on the FFT grid
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

def reject_spin_orbital_mode(input_dict):
    """Raise when [mode] enable_spin_orbital is set: the Eliashberg module
    does not support it, and until issue #83 it RAN anyway and printed
    eigenvalues built on four internal inconsistencies (norb never halved
    to the physical count; the internally computed chi0q in interleaved
    order against spin-block consumers -- measured 27% off, exactly a
    [0,2,1,3] permutation; interaction files read at the spin-orbital
    dimension, so U landed on (orb0, up)/(orb0, down) as two 'orbitals'
    and orbital 1 got none; paramagnetic orbital-space S/C rules applied
    to spin-orbital indices). A silent wrong result, not an approximation.

    Resolution is strict (shared _resolve_spin_orbital_flag): recognized
    false forms proceed as non-SO, unrecognized values raise -- a user who
    misspells "true" must not silently get a non-SO run that then
    misreads spin-orbital-formatted inputs.
    """
    if _resolve_spin_orbital_flag(input_dict):
        raise ValueError(
            "[mode] enable_spin_orbital = true is not supported by the "
            "Eliashberg module (hwave_sc): the pairing vertex is "
            "paramagnetic and the internal index/orbital-count "
            "conventions are not spin-orbital aware, so a run would "
            "produce silently wrong eigenvalues (issue #83). Use the "
            "UHFk/RPA/FLEX solvers for spin-orbital models.")


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
    reject_spin_orbital_mode(input_dict)
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
    from requests.structures import CaseInsensitiveDict
    g2_tail = _coerce_g2_tail(
        CaseInsensitiveDict(eli_param).get("g2_tail", True))
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
        # This path builds its own bubble (bond_bubble) directly from the
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
        bond_set = bond_channels.resolve_interactions(
            interactions["CoulombInter"], geom_info["rvec"], norb,
            bond_max_shells=bond_max_shells)
        logger.info("Bond channels: B = %d, Delta r = %s",
                    bond_set.n_channels, list(bond_set.delta_r))
        if bond_set.dropped:
            logger.info("Bond channel provenance (dropped): %s",
                        list(bond_set.dropped))

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
            bond_coeff_tail_requested, bond_green)
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
            bond_set, green_carrier, interactions, inter_k, geom_info, norb,
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
            chi0q_raw = _calc_chi0q_internal(input_dict, chi0q_tensor=chi0q_tensor,
                                                precomputed_mu=mu)
            # internally computed chi0q always carries the full frequency grid
            static_index = None
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
    # the finite-window part of _calc_g2, applied inside
    # bond_channels.make_bond_kernel), so G2 is not
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
