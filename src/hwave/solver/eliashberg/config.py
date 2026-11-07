"""Reading and validating the [eliashberg] configuration of hwave_sc.

Moved verbatim from hwave.sc (dev/sc-split-plan.md, step 1); the shell
re-exports every name here under its old spelling.
"""
from __future__ import annotations

import logging
import numpy as np

from hwave.solver.ir_axis import is_ir_native, ir_native_meta
from hwave.solver import bond_channels
from hwave.solver.eliashberg.vertex import (  # noqa: F401
    _NONBOND_FOCK_TYPES, _offsite_rows)

logger = logging.getLogger("hwave_sc")


# _warn_if_g2_indefinite skip thresholds (module-level so tests can patch
# them): eigvalsh work scales as (norb^2)^3 per k-point, i.e. norb^6 * nvol
# flops up to a constant. The peak memory is ~three full complex copies of
# the (nvol, norb^2, norb^2) view (the transpose/reshape materializes one,
# the residual |M - M^dag| evaluation another, the Hermitized matrix a
# third; tracemalloc-measured peak ~1.5x the two-copy figure), each
# 16 * norb^4 * nvol bytes.
_G2_CHECK_MAX_WORK = 20_000_000_000


_G2_CHECK_MAX_BYTES = 2_000_000_000


# Threshold for the asymptotic-regime diagnostic below: at the window edge
# |i wn G(k, i wn) - I| measures how far the Green function still is from
# its leading 1/(i wn) tail. 0.5 flags the regime where the neglected
# scales are at least comparable to the largest retained frequency.
_G2_TAIL_EDGE_DEV_THRESHOLD = 0.5


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


_NONBOND_REFUSED_OFFSITE = ("Exchange", "PairHop")


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


# Default peak-memory ceiling (GB) of the bond-channel preflight
# ([eliashberg] bond_memory_cap_gb; spec S3.2/S5).
_BOND_MEMORY_CAP_GB = 8.0


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


def _validate_dynamic_prereqs(input_dict):
    """Validate prerequisites for dynamic Eliashberg calculation.

    Parameters
    ----------
    input_dict : dict
        Parsed TOML configuration dictionary.

    Raises
    ------
    ValueError
        If a static-bond-only key is combined with ``bond_channels=true``, if
        chi0q_mode is not "flex", or if Nmat is odd.
    """
    eli = input_dict.get("eliashberg", {})
    # bond_channels = true on the DYNAMIC path selects the bond-resolved
    # pairing kernel, which takes its Green function and its bond topology
    # from the FLEX run, not from the static path's own keys. Refuse those
    # keys rather than accept and ignore them (a silently ignored bond_green
    # would leave the user believing a different propagator was used).
    if _bond_bool_option(eli, "bond_channels", False):
        for key in ("bond_green", "bond_max_shells"):
            if eli.get(key) is not None:
                raise ValueError(
                    "[eliashberg] {} is not used by the dynamic bond path "
                    "(frequency='dynamic' with bond_channels=true): the Green "
                    "function comes from path_to_flex_output and the topology "
                    "from the bond archive; remove the key".format(key))
        for key in ("zero_chi_s", "zero_chi_c"):
            # the same strict parser as bond_channels just above: read as
            # plain truthiness a typo like "ture" would slip past this
            # refusal and the unimplemented decomposition would be requested
            if _bond_bool_option(eli, key, False):
                raise ValueError(
                    "[eliashberg] {} is not implemented on the dynamic bond "
                    "path (frequency='dynamic' with bond_channels=true)"
                    .format(key))
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

    Every one of these is a "no silent wrong result" refusal: the shared
    bond modules (``bubble.bond_bubble_static``,
    ``bond_channels.build_sc_bond_channel``, ``eliashberg_bond``) stay
    general so the unit and reduction tests can exercise multi-orbital
    physics directly; only this top-level entry point restricts what a
    production run may request.
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

    # The bond channels are the off-site shells of every bond-resolved type
    # (CoulombInter, Hund, Ising, Exchange, PairHop). "Truly absent from the
    # config" is the documented trigger (spec S5): a DECLARED off-site row
    # whose value is zero is explicitly allowed, so that a V sweep keeps a
    # constant channel topology (B does not jump 1 <-> 5) at its V=0 point.
    def _rows(itype, offsite_only):
        return [(irvec, orbvec, value)
                for (irvec, orbvec), value in interactions.get(itype, {}).items()
                if not offsite_only or tuple(int(x) for x in irvec) != (0, 0, 0)]

    declared = [t for t in bond_channels._BOND_RESOLVABLE_TYPES
                if _rows(t, offsite_only=(t != "CoulombInter"))]
    if "CoulombInter" not in interactions and not declared:
        raise ValueError(
            "[eliashberg] bond_channels=true requires an off-site interaction "
            "(a CoulombInter term, or an off-site Hund, Ising, Exchange or "
            "PairHop row): bond channels carry the inter-site Fock/Cooper "
            "structure of the off-site interaction, so asking for them with "
            "no inter-site interaction declared is a configuration mistake. "
            "Declare the CoulombInter file (or a combined `Coulomb` file "
            "containing R!=0 entries) -- a declared-but-zero V is allowed "
            "and keeps the channel topology fixed across a V sweep -- or use "
            "bond_channels=false.")

    # CoulombInter is checked on every row (its R = 0 rows feed the on-site
    # block); the other bond-resolved types on their off-site rows, which are
    # the ones the bond channels carry.
    checked = [(t, _rows(t, offsite_only=(t != "CoulombInter")))
               for t in bond_channels._BOND_RESOLVABLE_TYPES]
    for itype, rows in checked:
        nonfinite = [(irvec, orbvec, value) for irvec, orbvec, value in rows
                     if not np.isfinite(complex(value).real)
                     or not np.isfinite(complex(value).imag)]
        if nonfinite:
            irvec, orbvec, value = nonfinite[0]
            raise ValueError(
                "[eliashberg] bond_channels=true requires every {} "
                "value to be finite, but V_{}{}({}) = {} is not (real or "
                "imaginary part is NaN/inf; {} such entries). A NaN/inf "
                "interaction value is a data/configuration error, not a "
                "physical model, and must not silently reach the bond "
                "machinery.".format(
                    itype, orbvec[0], orbvec[1], tuple(irvec), value,
                    len(nonfinite)))
    # An R = 0 CoulombInter row between an orbital and itself is a
    # same-orbital on-site density-density term. The bond-resolved vertex
    # (the locality-split channel-0 block shared with the FLEX bond gate)
    # carries such a row through its density (Hartree) slot only, while
    # exact diagonalization of the term also gives it the exchange crossing
    # of a Hubbard U (tests/test_bond_vs_ed_oracle.py, the "U" direction).
    # Refused rather than run with one of the two treatments: the on-site
    # same-orbital interaction belongs in CoulombIntra, where every path
    # treats it the same way.
    same_orbital = [(irvec, orbvec, value)
                    for irvec, orbvec, value in _rows("CoulombInter", False)
                    if tuple(int(x) for x in irvec) == (0, 0, 0)
                    and int(orbvec[0]) == int(orbvec[1]) and complex(value) != 0]
    if same_orbital:
        irvec, orbvec, value = same_orbital[0]
        raise ValueError(
            "[eliashberg] bond_channels=true does not accept an on-site "
            "same-orbital CoulombInter row, but V_{}{}({}) = {} is declared "
            "({} such entries). Declare an on-site interaction of an orbital "
            "with itself as CoulombIntra (a declared-zero row is "
            "accepted).".format(orbvec[0], orbvec[1], tuple(irvec), value,
                                len(same_orbital)))

    for itype, rows in checked:
        bad = [(irvec, orbvec, value) for irvec, orbvec, value in rows
               if abs(complex(value).imag) > imag_tol]
        if bad:
            irvec, orbvec, value = bad[0]
            raise ValueError(
                "[eliashberg] bond_channels=true supports only REAL "
                "inversion-symmetric {} in v1, but V_{}{}({}) = {} has a "
                "nonzero imaginary part ({} such entries). The bond-resolved "
                "vertex represents real off-site coefficients only; complex "
                "bonds are a deferred follow-up.".format(
                    itype, orbvec[0], orbvec[1], tuple(irvec), value,
                    len(bad)))


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


def _check_nonbond_offsite(interactions, where, site=None):
    """The off-site rows the non-bond (q-only) pairing vertex cannot represent.

    The non-bond paths keep only the Hartree (density) part of an ``R != 0``
    row (:func:`_build_interaction_k_for_sites`), as the RPA ring and the FLEX
    general path do; its exchange (Fock) crossing pairs a particle-hole pair on
    the bond and depends on two momenta, so it has no exact place in a vertex
    of ``q`` alone, and since issue #239 it is no longer collapsed into one.

    * An off-site ``Exchange`` or ``PairHop`` has no density part at all: it
      would contribute nothing, so it is REFUSED (``bond_channels = true``
      represents it exactly).
    * An off-site inter-orbital ``CoulombInter``, ``Hund`` or ``Ising`` keeps
      its Hartree part; ONE warning says that its Fock crossing is not
      included and points to the bond path.

    "Off-site" includes, with ``orbital_site`` labels ``site`` (#234), an
    ``R = 0`` row between orbitals of different sites."""
    refused = _offsite_rows(interactions, _NONBOND_REFUSED_OFFSITE, site)
    if refused:
        t = next(iter(refused))
        irvec, orbvec, v = refused[t][0]
        raise ValueError(
            "{}: an off-site {} row (R = {}, orbitals {}, value {}) has no density part, so "
            "the q-only pairing vertex of bond_channels = false cannot represent it (it would "
            "contribute nothing). Set [eliashberg] bond_channels = true, which treats every "
            "off-site row on its bond exactly, or remove the row.".format(
                where, t, tuple(int(x) for x in irvec), tuple(o + 1 for o in orbvec), v))
    fock = {t: [r for r in rows if r[1][0] != r[1][1]]
            for t, rows in _offsite_rows(interactions, _NONBOND_FOCK_TYPES, site).items()}
    fock = {t: rows for t, rows in fock.items() if rows}
    if fock:
        logger.warning(
            "%s: off-site inter-orbital %s declared. With bond_channels = false the pairing "
            "vertex keeps their Hartree (density) part only: their exchange (Fock) crossing "
            "is not included (earlier versions collapsed it into the q-only vertex, an "
            "approximation; issue #239). Set [eliashberg] bond_channels = true for the exact "
            "bond-resolved treatment.", where,
            ", ".join("{} ({} row(s))".format(t, len(r)) for t, r in fock.items()))
    return fock


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
    """Strictly parse a physics-relevant boolean config switch
    (:func:`hwave.solver.backend.strict_bool`)."""
    from hwave.solver.backend import strict_bool
    return strict_bool(value, name)


def _coerce_g2_tail(value):
    """Strictly parse the [eliashberg] g2_tail switch (round-3 review)."""
    return _coerce_config_bool(value, "[eliashberg] g2_tail")


def _resolve_spin_orbital_flag(input_dict):
    """Resolve [mode] enable_spin_orbital to a real boolean, strictly
    (:func:`hwave.solver.spin_orbital.resolve_flag`, the one resolver every
    solver uses; case-insensitive key lookup, unrecognized values raise)."""
    from hwave.solver import spin_orbital
    return spin_orbital.resolve_flag(input_dict.get("mode", {}))


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
