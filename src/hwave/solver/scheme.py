"""calc_scheme capability registry and the ``auto`` resolution rules (#167).

``calc_scheme="auto"`` promises a scheme that is EXACT for the declared
input, preferring the cheaper ``reduced`` scheme where exactness is
provable ("conservative exact"). Every consumer reads its OWN mode field
from :data:`CAPABILITIES` -- the modes are declared, never derived from
the vertex-family sets (a family-derived query misclassified PairLift).

Contents
--------
* :class:`Capability` / :data:`CAPABILITIES` -- one named record per
  reader interaction type, with explicit per-consumer modes.
* :func:`declared_types` -- fail-closed discovery over a reader table
  container (keys are canonicalized; an unknown type raises).
* :func:`flavour_conserved` -- the pre-fold structural predicate (RPA).
* :func:`resolve_rpa` / :func:`resolve_flex` -- the decision functions.
* :func:`legacy_1_0_resolution` -- the 1.0.x rule, kept verbatim as the
  WARNING baseline.
* :data:`RESOLUTION_TOKENS` -- the closed vocabulary stamped into npz.

This module imports no solver; it reads tables in the reader's single
value representation ``{((rx,ry,rz),(a,b)): complex}``.
"""
from __future__ import annotations

import math
from types import MappingProxyType
from typing import FrozenSet, Mapping, NamedTuple


class Capability(NamedTuple):
    """Per-type capability record.

    families : intrinsic vertex families (DOCUMENTATION only; the ring
        content of the type in hwave.solver.vertex_table terms).
    rpa_mode : "diag_ok" (never forces general on its own) |
        "conditional" (general iff flavour is not conserved) |
        "general_only" (unconditionally forces general).
    flex_forcing : FLEX auto resolves to general when True.
    sc_legacy_forcing : hwave_sc's auto rule counts the type as
        inter-orbital when True (frozen 1.0.x semantics).
    """
    families: FrozenSet[str]
    rpa_mode: str
    flex_forcing: bool
    sc_legacy_forcing: bool


RPA_MODES = frozenset({"diag_ok", "conditional", "general_only"})

_CAPABILITIES_RAW = {
    "CoulombIntra": Capability(frozenset({"diag"}), "diag_ok", False, False),
    "CoulombInter": Capability(frozenset({"density", "cross"}), "conditional", True, True),
    "Hund":         Capability(frozenset({"density", "cross"}), "conditional", True, True),
    "Ising":        Capability(frozenset({"density", "cross"}), "conditional", True, True),
    "Exchange":     Capability(frozenset({"cross"}), "general_only", True, True),
    "PairHop":      Capability(frozenset({"no_density_form"}), "general_only", True, True),
    # PairLift: rpa_mode conservative -- its ring content is real (#137,
    # ED-validated at first order; vertex_table.RING_SPIN_FLIP), so under
    # flavour conservation the reachability argument applies as to every
    # conditional type. flex_forcing=False: its FLEX particle-hole vertex is
    # exactly zero (inert). sc_legacy_forcing=False: frozen legacy behaviour.
    "PairLift":     Capability(frozenset({"cross"}), "conditional", False, False),
    # aggregate file (expands into Intra+Inter): conservatively cross-carrying.
    # Inspecting the expanded contents is a documented future relaxation.
    "Coulomb":      Capability(frozenset({"density", "cross"}), "conditional", True, True),
}
CAPABILITIES: Mapping[str, Capability] = MappingProxyType(dict(_CAPABILITIES_RAW))
del _CAPABILITIES_RAW

_LOWER_TO_CANONICAL = MappingProxyType({k.lower(): k for k in CAPABILITIES})


def _validate_registry():
    lowered = [k.lower() for k in CAPABILITIES]
    if len(set(lowered)) != len(lowered):
        raise ImportError("scheme.CAPABILITIES: canonical names collide case-insensitively")
    for name, cap in CAPABILITIES.items():
        if not isinstance(cap, Capability):
            raise ImportError("scheme.CAPABILITIES[{!r}] is not a Capability".format(name))
        if cap.rpa_mode not in RPA_MODES:
            raise ImportError("scheme.CAPABILITIES[{!r}]: rpa_mode {!r} not in {}".format(
                name, cap.rpa_mode, sorted(RPA_MODES)))
        if not isinstance(cap.flex_forcing, bool) or not isinstance(cap.sc_legacy_forcing, bool):
            raise ImportError("scheme.CAPABILITIES[{!r}]: consumer flags must be bool".format(name))


_validate_registry()

#: Sections a reader table container may hold that are NOT interactions.
#: Extern is a one-body term consumed by the predicate, never an interaction.
NON_INTERACTION_SECTIONS = frozenset(
    {"geometry", "transfer", "extern", "initial", "path_to_input"})


def canonical_name(name):
    """Canonical registry name for ``name`` (case-insensitive), or None."""
    return _LOWER_TO_CANONICAL.get(str(name).lower())


def is_structurally_nonzero(table, *, source="table"):
    """True iff any entry of ``table`` is nonzero. Every entry is inspected:
    a non-finite value (NaN / +-inf) raises even after a nonzero one."""
    nonzero = False
    for key, v in (table or {}).items():
        c = complex(v)
        if not (math.isfinite(c.real) and math.isfinite(c.imag)):
            raise ValueError(
                "{}: non-finite entry {!r} at {!r} -- calc_scheme resolution "
                "refuses to judge a table with NaN/inf".format(source, v, key))
        if c != 0:
            nonzero = True
    return nonzero


def declared_types(tables):
    """Canonical names of the structurally nonzero interaction tables in
    ``tables`` (a reader CaseInsensitiveDict or a plain dict). Fail-closed:
    a non-section key without a capability record raises."""
    out = set()
    for key in tables.keys():
        if str(key).lower() in NON_INTERACTION_SECTIONS:
            continue
        canon = canonical_name(key)
        if canon is None:
            raise ValueError(
                "no capability entry for interaction '{}': calc_scheme='auto' "
                "cannot judge it. Add a record to hwave.solver.scheme.CAPABILITIES "
                "or request calc_scheme explicitly.".format(key))
        if is_structurally_nonzero(tables[key], source=str(key)):
            out.add(canon)
    return frozenset(out)


def legacy_1_0_resolution(types, calc_type):
    """The 1.0.x auto rule (rpa.py ``_set_scheme`` before #167), kept
    VERBATIM as the baseline the promotion WARNING compares against.
    ``types`` are canonical names."""
    if calc_type == "ring+ladder":
        return "general"
    if "Exchange" in types or "PairHop" in types:
        return "general"
    return "reduced"


RESOLUTION_TOKENS = frozenset({
    "explicit",
    "auto:ring_ladder",
    "auto:general_only",
    "auto:no_discarded_content",
    "auto:exact:diagonal_transfer",
    "auto:exact:folded_diagonal",
    "auto:mixed:transfer",
    "auto:mixed:extern",
    "auto:mixed:trans_mod",
    "auto:mixed:green_init",
    "auto:flex_forcing",
})

CAUSE_PRECEDENCE = ("trans_mod", "green_init", "transfer", "extern")


def _find_key(tables, lower_name):
    for key in tables.keys():
        if str(key).lower() == lower_name:
            return key
    return None


def _has_nonzero_offdiagonal(tables, lower_name, *, index_limit=None):
    """True iff the section ``lower_name`` holds a nonzero entry with
    ``a != b``. Entries are judged individually (no cancellation credit).

    ``index_limit`` mirrors the CONSUMER's bounds: an entry carrying an
    index >= index_limit is fully invisible here, exactly as it is to
    ``rpa._make_ham_trans``, which skips it. The bounds filter therefore
    runs BEFORE the finiteness check -- a row that never enters H0(k)
    cannot promote the scheme, and cannot make the predicate refuse to
    judge the input either. Within the bounds, non-finite raises.
    """
    key = _find_key(tables, lower_name)
    if key is None:
        return False
    found = False
    for (irvec, orbvec), v in (tables[key] or {}).items():
        a, b = int(orbvec[0]), int(orbvec[1])
        if index_limit is not None and (a >= index_limit or b >= index_limit):
            continue
        c = complex(v)
        if not (math.isfinite(c.real) and math.isfinite(c.imag)):
            raise ValueError(
                "{}: non-finite entry {!r} at irvec={}, orbvec={}".format(
                    key, v, tuple(irvec), tuple(orbvec)))
        if c != 0 and a != b:
            found = True
    return found


def flavour_conserved(tables, *, norb_phys, coeff_extern,
                      trans_mod_present, green_init_present,
                      enable_spin_orbital=False):
    """Pre-fold structural flavour-conservation predicate (RPA).

    Returns ``(True, "diagonal_transfer")`` when the effective one-body
    Hamiltonian conserves the declared flavour, else ``(False, cause)`` with
    ``cause`` from :data:`CAUSE_PRECEDENCE` (first match wins; the
    precedence is fixed and tested). ``tables`` must be the PRE-FOLD
    container. In enable_spin_orbital mode the Transfer indices are the
    combined ``2*orb+spin`` index, so ``a == b`` is generalized-flavour
    conservation (a spin flip promotes).

    Both scans mirror their consumer's index bounds, so an entry H0(k)
    never sees is invisible to the decision:

    * Extern: ``_make_ham_trans`` reads physical orbital indices only, so
      the limit is ``norb_phys``.
    * Transfer: the normal branch of ``_make_ham_trans`` skips indices
      >= norb ("skip spin dependence"), while the spin-orbital branch
      consumes every index below the geometry norb = ``2*norb_phys``.
    """
    if not math.isfinite(float(coeff_extern)):
        raise ValueError("coeff_extern is non-finite: {!r}".format(coeff_extern))
    transfer_limit = 2 * int(norb_phys) if enable_spin_orbital else int(norb_phys)
    mixing = {
        "trans_mod": bool(trans_mod_present),
        "green_init": bool(green_init_present),
        "transfer": _has_nonzero_offdiagonal(tables, "transfer",
                                             index_limit=transfer_limit),
        "extern": (float(coeff_extern) != 0.0
                   and _has_nonzero_offdiagonal(tables, "extern",
                                                index_limit=int(norb_phys))),
    }
    for cause in CAUSE_PRECEDENCE:
        if mixing[cause]:
            return False, cause
    return True, "diagonal_transfer"


def resolve_rpa(types, calc_type, *, conserved, cause, has_sublattice):
    """RPA auto decision (deterministic precedence). ``types`` canonical."""
    if calc_type == "ring+ladder":
        return "general", "auto:ring_ladder"
    modes = {CAPABILITIES[t].rpa_mode for t in types}
    if "general_only" in modes:
        return "general", "auto:general_only"
    if "conditional" in modes:
        if conserved:
            return "reduced", ("auto:exact:folded_diagonal" if has_sublattice
                               else "auto:exact:diagonal_transfer")
        return "general", "auto:mixed:" + cause
    return "reduced", "auto:no_discarded_content"


def resolve_flex(types):
    """FLEX auto decision: general iff any declared type has flex_forcing."""
    if any(CAPABILITIES[t].flex_forcing for t in types):
        return "general", "auto:flex_forcing"
    return "reduced", "auto:no_discarded_content"


def estimate_chi_bytes(scheme, nmat, nvol, nd):
    """Size of the principal chi array for ``scheme``: complex128 over
    (nmat, nvol, nd^2) for reduced, (nmat, nvol, nd^4) for general.

    The name is matched case-insensitively: this is a pre-flight size
    estimate (the GPU path calls it before allocating), so an odd-cased but
    otherwise valid explicit request must size rather than crash. Anything
    that is not a RESOLVED scheme still raises."""
    scheme = str(scheme).lower()
    if scheme == "reduced":
        rank = int(nd) ** 2
    elif scheme == "general":
        rank = int(nd) ** 4
    else:
        raise ValueError("estimate_chi_bytes: unresolved scheme {!r}".format(scheme))
    return int(nmat) * int(nvol) * rank * 16
