"""The coverage report, the ranking and the summary.

* :func:`coverage_document` turns a :class:`CandidateSet` into the
  ``coverage.json`` document (``schema_version = 4``), with every key of
  the schema, and checks it with :func:`validate_coverage`.
* :func:`rank` ranks the collected rows on the free energy per primitive
  cell, ``FreeEnergy_Total / (Lx Ly Lz)`` (the UHFk free energy; not
  ``Energy_Total``, the internal energy), in one list per comparable group.
* :func:`order_parameter_records` reads the post-run order parameters of a
  run with :func:`hwave.results.read_uhfk_order_parameters`, per orbital
  and per site group.
* :func:`summary_text` renders the summary.

The only claim the report makes is the fixed sentence
(:data:`FIXED_SENTENCE`): the selected state is the lowest free energy
among the candidates tried.

**Row fields.** A row of the document has every key of :data:`ROW_FIELDS`:
the enumerator's fields (including ``pattern_label``, ``proposal_note``,
``duplicate_of`` and ``amplitude_requested``), the seed fields
(:data:`SEED_FIELDS`) and the post-run fields
(:data:`POST_RUN_FIELDS`); a key the row does not carry is null.

**Null rules** (:func:`validate_coverage`). The seed fields follow these
rules: para (tried) may set ``seed_base``, ``shell``,
``moments_h0`` and ``h_used = 0`` only; "no response (degenerate with
para)" ``seed_base``, ``shell``, ``moments_h0``, ``m0``, ``m_inf``; "no
response (constraint)" ``seed_base``, ``m0``; a writer refusal after
seeding whatever was computed; every other not-tried row (quota,
unrepresentable, duplicate, no response, a block-check or pattern refusal)
none; ``cap_reason`` only on a seeded row. ``amplitude_requested`` is the
candidate's requested amplitude, which the enumerator records on every row,
so it is exempt (it is not a result of seeding). The post-run fields are
null on every ``pre_run`` row and on every not-tried row; a tried
``post_run`` row has a ``status`` of :data:`STATUSES`.
"""
import dataclasses
import math
from fractions import Fraction
from typing import Literal, Mapping, Optional, Tuple

import numpy as np

from hwave.initial import cells, hints, labels
from hwave.initial.candidates import _thaw_json
from hwave.initial.facts import SITE_TOL, SiteGroups
from hwave.initial.seed import BLOCK_REASON, CONSTRAINT_REASON, DEGENERATE_REASON
from hwave.initial.symmetry import SZ_WARNING

__all__ = [
    "SCHEMA_VERSION", "FIXED_SENTENCE", "STRUCTURAL_EXCLUSIONS", "CHOSEN_BY",
    "DEGENERACY_TOL", "STATUSES", "SEED_FIELDS", "POST_RUN_FIELDS", "ROW_FIELDS",
    "PHASE_MARKERS", "FAILURE_KINDS", "coverage_document", "validate_coverage",
    "energy_fields",
    "RankEntry", "RankList", "Ranking", "rank", "failure_kind", "unranked_reason",
    "order_parameter_records",
    "orbital_centres", "summary_text",
]

#: ``coverage.json`` schema version. The first published version is 4 on
#: purpose: do not reset it to 1.
SCHEMA_VERSION = 4

#: the fixed sentence of every report, verbatim
FIXED_SENTENCE = ("The selected state is the lowest free energy among the candidates "
                  "tried. It is not a proof of the ground state.")

#: the structural exclusions: the kinds of state no candidate of this
#: version covers (they reach coverage.json and every summary)
STRUCTURAL_EXCLUSIONS = (
    "multi-q superpositions (double-Q, triple-Q, vortex crystals)",
    "sliding-phase variants beyond the two kept phases, including an Umklapp phase "
    "offset not pinned by symmetry",
    "combinations within a degenerate eigenspace other than the canonical basis",
    "rotated cells such as √3×√3, and cells outside subshapes",
    "orbital and bond patterns in the large tier (it has spin and charge patterns "
    "only), and bond patterns other than hints",
    "one amplitude per candidate",
    "hints describe the paramagnet only (static, at T_hint)",
    "spirals constructed only under a passed spin-rotation check",
    "Extern fields (unused by UHFk)",
    "q with a denominator above max_denominator",
)

#: who chose the tier; "auto (interim rule)" until the cost model exists
CHOSEN_BY = ("auto", "user", "auto (interim rule)")

#: the interim tier rule of this release
INTERIM_TIER_RULE = ("auto resolves to the large tier (no cost model in this release; "
                     "the small and medium tiers arrive later)")

#: the share of the host memory a tier may plan for
BUDGET_FRACTION = 0.25

#: near-degeneracy: |dF| per site (per cell when the grouping is ambiguous)
#: below this, in input energy units
DEGENERACY_TOL = 1e-4

#: the post-run statuses
STATUSES = ("converged", "not_converged", "convergence_unverified", "failed", "not_run",
            "seed_unverified")

PHASE_MARKERS = ("pre_run", "post_run")

#: the seed fields of a row
SEED_FIELDS = ("seed_base", "shell", "moments_h0", "m0", "m_inf", "m_target", "m_used",
               "h_hi", "h_used", "amplitude_requested", "amplitude_used", "cap_reason")

#: the post-run fields (filled by collect)
POST_RUN_FIELDS = ("run_dir", "status", "exit_code", "error_summary", "wall_time_s",
                   "iterations", "rest", "scf_residual", "scf_eps_tol", "scf_note",
                   "energy_total", "free_energy_total", "free_energy_per_cell",
                   "free_energy_per_site", "free_energy_per_site_reason",
                   "order_parameters", "final_label", "iteration_max", "eps", "mix",
                   "failure_kind")

#: why a ``failed`` row failed (``failure_kind``, set by collect; inferred
#: from ``iterations`` and ``exit_code`` when null): the solver crashed or
#: aborted, its outputs are missing or unreadable, or it performed no SCF
#: step (``IterationMax = 0``; the status file only)
FAILURE_KINDS = ("crash", "unreadable_output", "zero_steps")
_FAILURE_TEXT = {"crash": "solver crash or abort",
                 "unreadable_output": "outputs missing or unreadable",
                 "zero_steps": "zero SCF steps"}

#: every key of a row of the document
ROW_FIELDS = (
    "run_id", "candidate_id", "family", "sources", "q", "q_string", "phase", "sign",
    "channel", "pattern_label", "label", "rank", "lambda", "lambda_imag", "subshape",
    "cellshape", "representable", "reason", "proposed_subshape", "proposed_cellshape",
    "proposal_note", "equivalent_to", "duplicate_of", "citation", "seed",
) + SEED_FIELDS + ("tried", "not_tried_reason", "phase_marker") + POST_RUN_FIELDS

_MODEL_KEYS = ("input_sha256", "hwave_version", "hwave_initial_version", "spin_orbital",
               "norb", "norb_phys", "iflag_fock", "interacting_orbitals", "n_c",
               "site_groups", "Ncond", "Ncond_per_cell", "T", "2Sz", "cellshape", "su2",
               "warnings")
_HINT_KEYS = ("status", "skip_reason", "T_hint", "rule", "W0", "W_act", "W_all", "E_F",
              "D", "W_act_groups", "W_source", "q_status", "L_min", "t_hint_rule",
              "t_hint_source", "effect", "para_status", "para_free_energy",
              "para_residual", "para_polarisation", "para_eps", "gamma_self_adjoint",
              "method", "q_set")
#: null in a skipped hint record
_HINT_SKIPPED_NULL = ("para_free_energy", "para_residual", "para_polarisation",
                      "para_eps", "gamma_self_adjoint", "method")
_TIER_KEYS = ("name", "chosen_by", "estimate", "budget", "constants")
_LIMIT_KEYS = ("top_k", "max_candidates", "quota_instability", "precedence",
               "max_site_patterns", "n_random", "base_seed", "max_denominator",
               "amplitude", "T_hint", "structural_exclusions", "fixed_sentence")
_TOP_KEYS = ("schema_version", "phase_marker", "model", "tier", "hint", "cells", "rows",
             "limits", "ranking")
_SU2_STATUSES = ("passed", "failed", "not applicable (2Sz fixed)",
                 "not applicable (spin-orbital)")

_PARA_SEED = frozenset(("seed_base", "shell", "moments_h0", "h_used"))
_DEGENERATE_SEED = frozenset(("seed_base", "shell", "moments_h0", "m0", "m_inf"))
_CONSTRAINT_SEED = frozenset(("seed_base", "m0"))
_ALL_SEED = frozenset(SEED_FIELDS) - {"amplitude_requested"}


# --------------------------------------------------------------------------
# the document
# --------------------------------------------------------------------------

def _versions():
    import hwave
    # hwave.initial ships with H-wave (a preview until its documentation is
    # complete); it has no separate release number
    return hwave.__version__, "{} (preview)".format(hwave.__version__)


def _skipped_hint():
    out = {k: None for k in _HINT_KEYS}
    out.update(status="skipped", skip_reason="no hint record", q_status={}, q_set=[],
               E_F=[], D=[], W_act_groups=[])
    return out


def _tier_record(name, chosen_by):
    return {
        "name": name, "chosen_by": chosen_by,
        "estimate": {"bytes": None, "seconds": None, "class": None},
        "budget": {"budget_fraction": BUDGET_FRACTION, "bytes": None},
        "constants": {"budget_fraction": BUDGET_FRACTION},
        "rule": INTERIM_TIER_RULE if chosen_by == "auto (interim rule)" else None,
    }


def _t_hint_record(hint):
    return {
        "value": hint.get("T_hint"), "rule": hint.get("rule"),
        "source": hint.get("t_hint_source"), "W0": hint.get("W0"),
        "L_min": hint.get("L_min"),
        "constants": {"W_floor_act": hints.W_FLOOR_ACT,
                      "T_hint_rel_floor": hints.T_HINT_REL_FLOOR,
                      "structural_axis_max": hints.STRUCTURAL_AXIS_MAX,
                      "rule_mesh": hints.RULE_MESH, "rule_cap": hints.RULE_CAP},
    }


_GENERATED = ("baseline", "orbital", "instability")


def quota_boundary(rows) -> list:
    """Where a quota cut the q of one order unevenly (the generated rows of
    one ``SubShape``). The baseline order is round-robin over q by rank, so
    a quota that runs out part-way through a rank keeps some q of an order
    n and cuts others; lattice-equivalent q (such as the C4 partners
    (1/4, 1/4) and (1/4, 3/4)) can then be treated unequally, since no
    point-group map is used. One record per (SubShape, order) where a quota
    cut a row and the q of that order do not all have the same number of
    tried rows: ``{"subshape", "order", "tried_q", "untried_q", "counts"}``
    (duplicate rows are not counted)."""
    per = {}
    for r in rows:
        if r.get("family") not in _GENERATED or r.get("q") is None:
            continue
        why = r.get("not_tried_reason") or ""
        if why.startswith("duplicate of ") or r.get("duplicate_of"):
            continue
        if not r.get("tried") and not why.startswith("quota:"):
            continue
        B = tuple(r["subshape"])
        q = cells.to_qfrac(r["q"])
        rec = per.setdefault(B, {}).setdefault(cells.order(q), {})
        tried, cut = rec.get(cells.q_string(q), (0, 0))
        rec[cells.q_string(q)] = (tried + bool(r["tried"]),
                                  cut + (not r["tried"]))
    out = []
    for B, orders in per.items():
        for n in sorted(orders):
            counts = orders[n]
            if not any(c for _, c in counts.values()):
                continue
            if len({t for t, _ in counts.values()}) < 2:
                continue
            out.append({
                "subshape": list(B), "order": n,
                "tried_q": [q for q, (t, _) in counts.items() if t],
                "untried_q": [q for q, (t, _) in counts.items() if not t],
                "counts": {q: {"tried": t, "cut": c} for q, (t, c) in counts.items()},
            })
    return out


def _full_row(row, phase_marker):
    out = {k: None for k in ROW_FIELDS}
    out.update(_thaw_json(row))
    out["phase_marker"] = phase_marker
    return out


def coverage_document(cset, *, phase_marker: Literal["pre_run", "post_run"],
                      tier_chosen_by: str = "auto (interim rule)", ranking=None) -> dict:
    """The ``coverage.json`` document of ``cset``.

    ``phase_marker``: ``"pre_run"`` (prepare) or ``"post_run"`` (collect;
    the rows then carry the post-run fields). ``tier_chosen_by``: one of
    :data:`CHOSEN_BY` (``"user"`` when the tier was given explicitly).
    ``ranking``: a :class:`Ranking` stored under ``"ranking"`` (null when
    not given).

    The hint record moves from ``model["hint"]`` to the top-level
    ``hint``; ``limits`` gains ``T_hint`` (value, rule and constants), the
    label thresholds (:data:`hwave.initial.labels.THRESHOLDS`),
    ``degeneracy_tol``, ``site_tol``,
    ``quota_boundary`` (:func:`quota_boundary`), the structural exclusions
    and the fixed sentence. The document is checked by
    :func:`validate_coverage` before it is returned (``ValueError``)."""
    if phase_marker not in PHASE_MARKERS:
        raise ValueError("phase_marker: must be 'pre_run' or 'post_run', got {!r}"
                         .format(phase_marker))
    if tier_chosen_by not in CHOSEN_BY:
        raise ValueError("tier_chosen_by: must be one of {}, got {!r}".format(
            ", ".join(repr(c) for c in CHOSEN_BY), tier_chosen_by))
    model = _thaw_json(cset.model)
    hint = model.pop("hint", None) or _skipped_hint()
    hv, iv = _versions()
    model.setdefault("hwave_version", hv)
    model.setdefault("hwave_initial_version", iv)
    rows = [_full_row(r, phase_marker) for r in cset.rows]
    limits = _thaw_json(cset.limits)
    limits.update(
        T_hint=_t_hint_record(hint),
        label_thresholds=dict(labels.THRESHOLDS),
        degeneracy_tol=DEGENERACY_TOL, site_tol=SITE_TOL,
        quota_boundary=quota_boundary(rows),
        structural_exclusions=list(STRUCTURAL_EXCLUSIONS),
        fixed_sentence=FIXED_SENTENCE)
    doc = {
        "schema_version": SCHEMA_VERSION, "phase_marker": phase_marker,
        "model": model, "tier": _tier_record(cset.tier, tier_chosen_by), "hint": hint,
        "cells": _thaw_json(cset.cells), "rows": rows, "limits": limits,
        "ranking": None if ranking is None else ranking.to_json(),
    }
    validate_coverage(doc)
    return doc


# --------------------------------------------------------------------------
# validation
# --------------------------------------------------------------------------

def _missing(obj, keys, where, errs):
    if not isinstance(obj, Mapping):
        errs.append("{}: must be an object".format(where))
        return False
    ok = True
    for k in keys:
        if k not in obj:
            errs.append("{}.{}: required key missing".format(where, k))
            ok = False
    return ok


def _allowed_seed(row):
    """The seed fields ``row`` may set (the null rules of the module text)."""
    if row.get("tried"):
        if row.get("family") == "para" or row.get("channel") == "para":
            return _PARA_SEED
        return _ALL_SEED
    why = row.get("not_tried_reason") or ""
    if why == DEGENERATE_REASON:
        return _DEGENERATE_SEED
    if why == CONSTRAINT_REASON:
        return _CONSTRAINT_SEED
    if why.startswith("refused:") and not why.startswith(BLOCK_REASON) \
            and not why.startswith("refused: pattern"):
        return _ALL_SEED - {"cap_reason"}       # refused after seeding: as computed
    return frozenset()


def _check_row(i, row, phase_marker, errs):
    where = "rows[{}] ({})".format(i, row.get("run_id") if isinstance(row, Mapping) else "?")
    if not _missing(row, ROW_FIELDS, where, errs):
        return
    tried = row["tried"]
    if not isinstance(tried, bool):
        errs.append("{}.tried: must be a bool".format(where))
        return
    if tried != (row["not_tried_reason"] is None):
        errs.append("{}: not_tried_reason must be null exactly when tried (tried = {}, "
                    "reason {!r})".format(where, tried, row["not_tried_reason"]))
    if row["phase_marker"] != phase_marker:
        errs.append("{}.phase_marker: {!r} differs from the document's {!r}".format(
            where, row["phase_marker"], phase_marker))
    allowed = _allowed_seed(row)
    for k in _ALL_SEED:
        if k not in allowed and row.get(k) is not None:
            errs.append("{}.{}: must be null for this row ({})".format(
                where, k, row["not_tried_reason"] or row["family"]))
    if tried and allowed is _PARA_SEED and row.get("seed_base") is not None \
            and (isinstance(row.get("h_used"), bool) or row.get("h_used") != 0):
        # a para row whose seed is recorded (seed_base set) has h_used = 0
        # exactly; before the seed is merged every seed field is null
        errs.append("{}.h_used: a seeded para row has h_used = 0, got {!r}".format(
            where, row.get("h_used")))
    fk = row.get("failure_kind")
    if fk is not None and (row.get("status") != "failed" or fk not in FAILURE_KINDS):
        errs.append("{}.failure_kind: one of {} on a failed row only, got {!r} (status "
                    "{!r})".format(where, ", ".join(FAILURE_KINDS), fk, row.get("status")))
    if phase_marker == "pre_run" or not tried:
        for k in POST_RUN_FIELDS:
            if row.get(k) is not None:
                errs.append("{}.{}: post-run fields are null {}".format(
                    where, k, "before the runs" if phase_marker == "pre_run"
                    else "on a row that was not tried"))
    elif row["status"] not in STATUSES:
        errs.append("{}.status: a tried post-run row needs one of {}, got {!r}".format(
            where, ", ".join(STATUSES), row["status"]))


def _check_model(model, errs):
    if not _missing(model, _MODEL_KEYS, "model", errs):
        return
    _missing(model["site_groups"], ("source", "groups", "ambiguous"),
             "model.site_groups", errs)
    su2 = model["su2"]
    if _missing(su2, ("status", "max_residual", "method", "sz_residual"), "model.su2",
                errs):
        st = su2["status"]
        if st not in _SU2_STATUSES:
            errs.append("model.su2.status: must be one of {}, got {!r}".format(
                ", ".join(_SU2_STATUSES), st))
        if st.startswith("not applicable") if isinstance(st, str) else False:
            for k in ("max_residual", "method"):
                if su2[k] is not None:
                    errs.append("model.su2.{}: null when not applicable".format(k))
        if st != "not applicable (2Sz fixed)" and su2["sz_residual"] is not None:
            errs.append("model.su2.sz_residual: set only in the 2Sz case")


def _check_hint(hint, errs):
    if not _missing(hint, _HINT_KEYS, "hint", errs):
        return
    if hint["status"] not in ("computed", "skipped"):
        errs.append("hint.status: must be 'computed' or 'skipped', got {!r}".format(
            hint["status"]))
    if hint["status"] == "skipped":
        if not hint["skip_reason"]:
            errs.append("hint.skip_reason: a skipped hint names its reason")
        for k in _HINT_SKIPPED_NULL:
            if hint[k] is not None:
                errs.append("hint.{}: null when the hint is skipped".format(k))
        if hint["q_set"] != []:
            errs.append("hint.q_set: [] when the hint is skipped")
    elif hint["skip_reason"] is not None:
        errs.append("hint.skip_reason: null when the hint is computed")


def validate_coverage(doc) -> None:
    """Check a ``coverage.json`` document: every key of the schema
    present, ``schema_version = 4``, the null rules of the module text,
    ``tier.chosen_by`` in :data:`CHOSEN_BY`, the structural exclusions and
    the fixed sentence verbatim. ``ValueError`` listing every problem."""
    errs = []
    if not _missing(doc, _TOP_KEYS, "$", errs):
        raise ValueError("; ".join(errs))
    if doc.get("schema_version") != SCHEMA_VERSION:
        errs.append("schema_version: must be {}, got {!r}".format(
            SCHEMA_VERSION, doc.get("schema_version")))
    pm = doc.get("phase_marker")
    if pm not in PHASE_MARKERS:
        errs.append("phase_marker: must be 'pre_run' or 'post_run', got {!r}".format(pm))
    if "model" in doc:
        _check_model(doc["model"], errs)
    if "hint" in doc:
        _check_hint(doc["hint"], errs)
    tier = doc.get("tier")
    if _missing(tier, _TIER_KEYS, "tier", errs) and tier["chosen_by"] not in CHOSEN_BY:
        errs.append("tier.chosen_by: must be one of {}, got {!r}".format(
            ", ".join(repr(c) for c in CHOSEN_BY), tier["chosen_by"]))
    if "tier" in doc and isinstance(tier, Mapping) and "estimate" in tier:
        _missing(tier["estimate"], ("bytes", "seconds", "class"), "tier.estimate", errs)
    if not isinstance(doc.get("cells"), list):
        errs.append("cells: must be a list")
    rows = doc.get("rows")
    if not isinstance(rows, list):
        errs.append("rows: must be a list")
    else:
        for i, row in enumerate(rows):
            _check_row(i, row, pm, errs)
    limits = doc.get("limits")
    if _missing(limits, _LIMIT_KEYS, "limits", errs):
        if list(limits["structural_exclusions"] or ()) != list(STRUCTURAL_EXCLUSIONS):
            errs.append("limits.structural_exclusions: must be the ten v1 items verbatim")
        if limits["fixed_sentence"] != FIXED_SENTENCE:
            errs.append("limits.fixed_sentence: must be {!r}".format(FIXED_SENTENCE))
    if errs:
        raise ValueError("coverage document: " + "; ".join(errs))


# --------------------------------------------------------------------------
# energies and the ranking
# --------------------------------------------------------------------------

def _facts_record(facts) -> dict:
    """The fields the ranking reads, from a ModelFacts or a model record."""
    if isinstance(facts, Mapping):
        f = dict(facts)
        sg = f.get("site_groups") or {}
        if isinstance(sg, SiteGroups):
            sg = {"source": sg.source, "groups": sg.groups, "ambiguous": sg.ambiguous}
        return {"T": float(f["T"]), "Ncond_per_cell": f.get("Ncond_per_cell"),
                "2Sz": f.get("2Sz"), "cellshape": [int(x) for x in f["cellshape"]],
                "interacting_orbitals": [int(a) for a in f.get("interacting_orbitals", ())],
                "site_groups": {"source": sg.get("source"),
                                "groups": [list(g) for g in sg.get("groups", ())],
                                "ambiguous": list(sg.get("ambiguous", ()))}}
    sg = facts.site_groups
    npc = facts.ncond_per_cell
    return {"T": float(facts.T), "Ncond_per_cell": str(npc) if npc is not None else None,
            "2Sz": facts.two_sz, "cellshape": [int(x) for x in facts.cellshape],
            "interacting_orbitals": [int(a) for a in facts.interacting_orbitals],
            "site_groups": {"source": sg.source, "groups": [list(g) for g in sg.groups],
                            "ambiguous": list(sg.ambiguous)}}


def _sites_per_cell(f):
    """``(number of interacting site groups, reason when per site is null)``."""
    sg = f["site_groups"]
    if sg["ambiguous"] and sg["source"] == "clustered":
        return None, ("site grouping is ambiguous ({}); no user grouping is given, so "
                      "F per site is not defined and near-degeneracy is judged per cell"
                      .format("; ".join(sg["ambiguous"])))
    ic = set(f["interacting_orbitals"])
    n = sum(1 for g in sg["groups"] if any(int(a) in ic for a in g))
    if n == 0:
        return None, "no interacting site group (n_c = 0): F per site is not defined"
    return n, None


def _row_T(row, f):
    return float(row["T"]) if row.get("T") is not None else f["T"]


def energy_fields(row, *, facts) -> dict:
    """The ranking key of a collected row:
    ``free_energy_per_cell = FreeEnergy_Total / (Lx Ly Lz)`` with ``L`` the
    row's ``cellshape``; at T = 0 ``Energy_Total`` when the free energy is
    absent (an older H-wave); at T > 0 then no key (``refused`` names the
    reason). ``free_energy_per_site`` divides by the number of interacting
    site groups per primitive cell; null with ``free_energy_per_site_reason``
    when the grouping is ambiguous (and not given by the user) or there is
    no interacting site. The row's T is its ``"T"`` key when present, else
    the model's. ``facts``: a ModelFacts or the document's model record.

    Returns ``{"free_energy_per_cell", "free_energy_per_site",
    "free_energy_per_site_reason", "key_source", "refused"}``."""
    f = _facts_record(facts)
    vol = float(np.prod(row.get("cellshape") or f["cellshape"]))
    T = _row_T(row, f)
    out = {"free_energy_per_cell": None, "free_energy_per_site": None,
           "free_energy_per_site_reason": None, "key_source": None, "refused": None}
    F, E = row.get("free_energy_total"), row.get("energy_total")
    if F is not None:
        out.update(free_energy_per_cell=float(F) / vol, key_source="FreeEnergy_Total")
    elif T == 0.0 and E is not None:
        out.update(free_energy_per_cell=float(E) / vol, key_source="Energy_Total (T = 0)")
    elif E is not None:
        out["refused"] = ("FreeEnergy_Total is absent (an older H-wave) and T > 0: "
                          "Energy_Total is the internal energy, not the free energy")
    else:
        out["refused"] = "no FreeEnergy_Total or Energy_Total in energy.dat"
    key = out["free_energy_per_cell"]
    if key is not None and not math.isfinite(key):
        out.update(free_energy_per_cell=None, key_source=None,
                   refused="the energy in energy.dat is not finite ({})".format(key))
    n, reason = _sites_per_cell(f)
    if out["free_energy_per_cell"] is not None:
        if n is None:
            out["free_energy_per_site_reason"] = reason
        else:
            out["free_energy_per_site"] = out["free_energy_per_cell"] / n
    return out


@dataclasses.dataclass(frozen=True)
class RankEntry:
    """One ranked row. ``delta_per_cell``/``delta_per_site``: above the
    lowest of its list. ``flags``: ``"near-degenerate"``, ``"cap hit"``,
    ``"convergence unverified"``. ``near_degenerate_with``: the run ids
    within :data:`DEGENERACY_TOL` of it."""
    run_id: str
    candidate_id: Optional[str]
    pattern_label: Optional[str]
    final_label: Optional[str]
    status: str
    free_energy_per_cell: float
    free_energy_per_site: Optional[float]
    free_energy_per_site_reason: Optional[str]
    delta_per_cell: float
    delta_per_site: Optional[float]
    key_source: str
    flags: Tuple[str, ...]
    near_degenerate_with: Tuple[str, ...]
    group: Optional[Mapping] = None

    def to_json(self):
        out = dataclasses.asdict(self)
        out["flags"] = list(self.flags)
        out["near_degenerate_with"] = list(self.near_degenerate_with)
        out["group"] = None if self.group is None else dict(self.group)
        return out


@dataclasses.dataclass(frozen=True)
class RankList:
    """One ranking list. ``kind``: ``"converged"`` or
    ``"convergence_unverified"``. ``group``: ``{"T", "Ncond_per_cell",
    "cellshape", "2Sz"}`` shared by its rows (``None`` for a flagged list
    that mixes groups). ``degeneracy_scale``: ``"per site"`` or ``"per
    cell"``."""
    kind: str
    group: Optional[Mapping]
    entries: Tuple[RankEntry, ...]
    degeneracy_scale: str

    def to_json(self):
        return {"kind": self.kind, "group": None if self.group is None else dict(self.group),
                "degeneracy_scale": self.degeneracy_scale,
                "entries": [e.to_json() for e in self.entries]}


@dataclasses.dataclass(frozen=True)
class Ranking:
    """The ranking of the collected rows. ``lists``: the converged rows, one
    list per (T, Ncond per cell, CellShape, 2Sz), the model's group first.
    ``flagged_unverified``: the ``convergence_unverified`` rows on the same
    key. ``unranked``: every other tried row (``not_converged``,
    ``failed``, ``not_run``, ``seed_unverified``, not collected, or without
    a usable energy), each a copy of its row with ``unranked_reason``;
    listed, never dropped. ``warnings``: the group splits."""
    lists: Tuple[RankList, ...]
    flagged_unverified: RankList
    unranked: Tuple[Mapping, ...]
    warnings: Tuple[str, ...]
    degeneracy_tol: float = DEGENERACY_TOL

    def to_json(self):
        return {"lists": [l.to_json() for l in self.lists],
                "flagged_unverified": self.flagged_unverified.to_json(),
                "unranked": [dict(r) for r in self.unranked],
                "warnings": list(self.warnings), "degeneracy_tol": self.degeneracy_tol,
                "degeneracy_tol_unit": "input energy units"}


def _canon_ncond(v):
    """Ncond per cell as a canonical value: an exact ``Fraction`` for an
    int, a rational string or a float that is one (``1``, ``"1"`` and
    ``1.0`` agree), else the float."""
    if v is None:
        return None
    if isinstance(v, bool):
        raise ValueError("Ncond per cell: a number, got {!r}".format(v))
    if isinstance(v, Fraction):
        return v
    if isinstance(v, str):
        try:
            return Fraction(v.strip())
        except (ValueError, ZeroDivisionError):
            raise ValueError("Ncond per cell: a number, got {!r}".format(v)) from None
    x = float(v)
    fr = Fraction(x).limit_denominator(10 ** 6)
    return fr if float(fr) == x else Fraction(x)


def _canon_sz(v):
    """2Sz as ``None`` (free: unset or ``"free"`` in any case) or an int."""
    if v is None or (isinstance(v, str) and v.strip().lower() == "free"):
        return None
    if isinstance(v, bool):
        raise ValueError("2Sz: an integer or 'free', got {!r}".format(v))
    try:
        x = float(v)
    except (TypeError, ValueError):
        raise ValueError("2Sz: an integer or 'free', got {!r}".format(v)) from None
    if not x.is_integer():
        raise ValueError("2Sz: an integer or 'free', got {!r}".format(v))
    return int(x)


def _canon(T, ncond, cellshape, sz):
    return (float(T), _canon_ncond(ncond), tuple(int(x) for x in cellshape), _canon_sz(sz))


def _group(row, f):
    """The comparison group of a row, in canonical form (T a float, Ncond
    per cell a Fraction, 2Sz None or an int), so equal values of different
    types (``1`` and ``"1"``) are one group."""
    return _canon(_row_T(row, f),
                  row["Ncond_per_cell"] if row.get("Ncond_per_cell") is not None
                  else f["Ncond_per_cell"],
                  row.get("cellshape") or f["cellshape"],
                  row["2Sz"] if "2Sz" in row else f["2Sz"])


def _ncond_text(v):
    return None if v is None else str(v)


def _sz_text(v):
    return "free" if v is None else str(v)


def _group_json(g):
    return {"T": g[0], "Ncond_per_cell": _ncond_text(g[1]), "cellshape": list(g[2]),
            "2Sz": g[3]}


def _build_list(kind, items, group, tol):
    """``items``: (row, energy fields) of one list."""
    items = sorted(items, key=lambda it: (it[1]["free_energy_per_cell"], it[0]["run_id"]))
    per_site = all(e["free_energy_per_site"] is not None for _, e in items)
    scale = "per site" if per_site else "per cell"
    key = "free_energy_per_site" if per_site else "free_energy_per_cell"
    if not items:
        return RankList(kind, group, (), scale)
    lo_cell = items[0][1]["free_energy_per_cell"]
    lo_site = items[0][1]["free_energy_per_site"]
    near = {i: [] for i in range(len(items))}
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            if abs(items[j][1][key] - items[i][1][key]) < tol:
                near[i].append(items[j][0]["run_id"])
                near[j].append(items[i][0]["run_id"])
    entries = []
    for i, (row, e) in enumerate(items):
        flags = []
        if kind == "convergence_unverified":
            flags.append("convergence unverified")
        if near[i]:
            flags.append("near-degenerate")
        if row.get("cap_reason"):
            flags.append("cap hit")
        site = e["free_energy_per_site"]
        entries.append(RankEntry(
            run_id=row["run_id"], candidate_id=row.get("candidate_id"),
            pattern_label=row.get("pattern_label"), final_label=row.get("final_label"),
            status=row["status"], free_energy_per_cell=e["free_energy_per_cell"],
            free_energy_per_site=site,
            free_energy_per_site_reason=e["free_energy_per_site_reason"],
            delta_per_cell=e["free_energy_per_cell"] - lo_cell,
            delta_per_site=None if site is None or lo_site is None else site - lo_site,
            key_source=e["key_source"], flags=tuple(flags),
            near_degenerate_with=tuple(near[i]), group=group))
    return RankList(kind, group, tuple(entries), scale)


def _differences(g, main):
    out = []
    names = ("T", "Ncond per cell", "CellShape", "2Sz")
    for name, a, b in zip(names, g, main):
        if a != b:
            out.append((name, a, b))
    return out


def _split_warning(g, main, n):
    parts = []
    for name, a, b in _differences(g, main):
        if name == "CellShape":
            parts.append("CellShape {} differs from the model's {} (a k-mesh mismatch)"
                         .format(list(a), list(b)))
        elif name == "2Sz":
            parts.append("2Sz {} differs from the model's {}".format(_sz_text(a),
                                                                       _sz_text(b)))
        else:
            parts.append("{} = {} differs from the model's {}".format(
                name, _ncond_text(a) if name == "Ncond per cell" else a,
                _ncond_text(b) if name == "Ncond per cell" else b))
    return ("{} row(s) with {} are ranked in a separate list; they are not comparable "
            "with the model's list".format(n, "; ".join(parts)))


def failure_kind(row) -> Optional[str]:
    """The :data:`FAILURE_KINDS` entry of a ``failed`` row: its
    ``failure_kind`` when set, else ``"zero_steps"`` for ``iterations = 0``,
    ``"crash"`` for a non-zero ``exit_code``, else ``"unreadable_output"``;
    ``None`` for any other status."""
    if row.get("status") != "failed":
        return None
    if row.get("failure_kind") is not None:
        return row["failure_kind"]
    if row.get("iterations") == 0:
        return "zero_steps"
    if row.get("exit_code") not in (None, 0):
        return "crash"
    return "unreadable_output"


def unranked_reason(row) -> str:
    """Why a tried row is not ranked, specific to its status: the failure
    kind and ``error_summary`` of a failed row (a crash, unreadable outputs
    and zero SCF steps are told apart), the iterations and residual of a
    not-converged row, the check of a seed-unverified row."""
    status = row.get("status")
    summary = row.get("error_summary")
    if status is None:
        return "not collected (no post-run status)"
    if status == "failed":
        kind = failure_kind(row)
        text = "failed ({})".format(_FAILURE_TEXT[kind])
        if kind == "crash" and row.get("exit_code") not in (None, 0):
            text += ", exit code {}".format(row["exit_code"])
        if summary:
            text += ": {}".format(summary)
        elif kind == "zero_steps":
            text += ": no SCF step performed"
        return text
    if status == "not_converged":
        parts = []
        if row.get("iterations") is not None:
            parts.append("{} iterations{}".format(
                row["iterations"], "" if row.get("iteration_max") is None
                else " of IterationMax = {}".format(row["iteration_max"])))
        if row.get("scf_residual") is not None:
            parts.append("residual {}".format(_num(row["scf_residual"], "{:.3e}")))
        if row.get("scf_eps_tol") is not None:
            parts.append("tolerance {}".format(_num(row["scf_eps_tol"], "{:.3e}")))
        text = "not converged"
        if parts:
            text += " ({})".format(", ".join(parts))
        return text + (": {}".format(summary) if summary else "")
    if status == "not_run":
        return "not run (neither hwave.log nor an output exists)"
    if status == "seed_unverified":
        return "seed unverified: {}".format(
            summary or "an input or the initial file differs from the one prepared, or "
                       "the initial-state check refused it")
    return "status {}: not ranked".format(status)


def rank(rows, *, facts, degeneracy_tol: float = DEGENERACY_TOL) -> Ranking:
    """Rank the collected rows.

    Only ``converged`` rows enter a ranking list, ordered by
    :func:`energy_fields` ``free_energy_per_cell`` (ties by run id). Rows
    are ranked together only when T, Ncond per cell, CellShape and 2Sz
    agree (a row's own ``"T"``, ``"Ncond_per_cell"`` and ``"2Sz"`` keys
    override the model's; CellShape is the row's ``cellshape``); every other
    group is a separate list, with a warning. ``convergence_unverified``
    rows form the flagged list on the same key. Near-degeneracy: ``|dF|``
    per site below ``degeneracy_tol`` (per cell when a row of the list has
    no per-site value) flags both rows. A cap-hit row (``cap_reason`` set)
    is ranked and flagged. Not-tried rows are not runs and are left out
    (the coverage table lists them); every other tried row is
    ``unranked`` with its reason.

    ``facts``: a ModelFacts or the document's model record."""
    tol = float(degeneracy_tol)
    if not math.isfinite(tol) or tol < 0:
        raise ValueError("degeneracy_tol: a finite number >= 0, got {!r}".format(
            degeneracy_tol))
    f = _facts_record(facts)
    main = _canon(f["T"], f["Ncond_per_cell"], f["cellshape"], f["2Sz"])
    groups, unverified, unranked = {}, {}, []
    for row in rows:
        if not row.get("tried"):
            continue
        status = row.get("status")
        if status not in ("converged", "convergence_unverified"):
            unranked.append(dict(row, unranked_reason=unranked_reason(row)))
            continue
        e = energy_fields(row, facts=f)
        if e["refused"]:
            unranked.append(dict(row, unranked_reason=e["refused"]))
            continue
        if status == "converged":
            groups.setdefault(_group(row, f), []).append((row, e))
        else:
            unverified.setdefault(_group(row, f), []).append((row, e))
    order = sorted(groups, key=lambda g: (g != main, repr(g)))
    lists, warnings = [], []
    for g in order:
        lists.append(_build_list("converged", groups[g], _group_json(g), tol))
        if g != main:
            warnings.append(_split_warning(g, main, len(groups[g])))
    uorder = sorted(unverified, key=lambda g: (g != main, repr(g)))
    if len(uorder) <= 1:
        g = uorder[0] if uorder else main
        flagged = _build_list("convergence_unverified", unverified.get(g, []),
                              _group_json(g), tol)
    else:
        # deltas and near-degeneracy per comparison group; each entry names its group
        parts = [_build_list("convergence_unverified", unverified[g], _group_json(g), tol)
                 for g in uorder]
        scale = ("per site" if all(p.degeneracy_scale == "per site" for p in parts)
                 else "per cell")
        flagged = RankList("convergence_unverified", None,
                           tuple(e for p in parts for e in p.entries), scale)
        warnings.append("the convergence-unverified rows span {} comparison groups; their "
                        "deltas and near-degeneracy flags are computed within each group"
                        .format(len(uorder)))
    return Ranking(tuple(lists), flagged, tuple(unranked), tuple(warnings), tol)


# --------------------------------------------------------------------------
# post-run order parameters
# --------------------------------------------------------------------------

def orbital_centres(facts) -> np.ndarray:
    """The fractional centre of each physical orbital of the primitive cell
    (``(norb_phys, 3)``; in spin-orbital mode the shared centre of the two
    spin-orbitals), for :func:`order_parameter_records`."""
    from hwave.initial.facts import _physical_centers
    return np.array(_physical_centers(facts.problem)[0], dtype=float)


def _c(z):
    z = complex(z)
    return [float(z.real), float(z.imag)]


def _qs(q):
    return None if q is None else cells.q_string(q)


def _arr(a):
    return None if a is None else np.asarray(a, dtype=float).tolist()


def _op_json(op) -> dict:
    """An :class:`hwave.results.OrderParameters` as JSON (complex numbers as
    ``[re, im]``, q as ``[0, 1)`` strings)."""
    return {
        "site_of_orbital": np.asarray(op.site_of_orbital).astype(int).tolist(),
        "site_positions": _arr(op.site_positions),
        "charge": _arr(op.charge), "spin": _arr(op.spin), "moment": _arr(op.moment),
        "uniform": _arr(op.uniform), "uniform_norm": float(op.uniform_norm),
        "mean_local_moment": float(op.mean_local_moment),
        "q": None if op.q_frac is None else [cells.q_string(q) for q in op.q_frac],
        "spin_fourier": None if op.spin_fourier is None else
        [[_c(z) for z in v] for v in op.spin_fourier],
        "spin_amplitude": _arr(op.spin_amplitude),
        "spin_ordering_vector": _qs(op.spin_ordering_vector),
        "charge_fourier": None if op.charge_fourier is None else
        [_c(z) for z in op.charge_fourier],
        "charge_amplitude": _arr(op.charge_amplitude),
        "charge_ordering_vector": _qs(op.charge_ordering_vector),
        "charge_order_amplitude": None if op.charge_order_amplitude is None
        else float(op.charge_order_amplitude),
        "staggered": None if op.staggered is None else [_c(z) for z in op.staggered],
        "staggered_norm": None if op.staggered_norm is None else float(op.staggered_norm),
        "cell_shape": None if op.cell_shape is None else list(op.cell_shape),
        "lattice_shape": None if op.lattice_shape is None else list(op.lattice_shape),
        "n_cell": int(op.n_cell), "transverse_available": bool(op.transverse_available),
        "units": dict(op.units),
    }


def order_parameter_records(run_or_toml, *, site_groups, orbital_centres,
                            base_dir=None) -> dict:
    """The post-run order parameters of one UHFk run (the spin and charge
    amplitudes M(Q) and dn(Q), per orbital), read twice with
    :func:`hwave.results.read_uhfk_order_parameters`:

    * ``per_orbital``: ``site_of_orbital`` the identity, each physical
      orbital its own site;
    * ``per_site``: the site groups (``site_groups``, a :class:`SiteGroups`
      or the model record ``{"source", "groups", "ambiguous"}``); null, with
      ``per_site_reason``, when a clustered grouping is ambiguous.

    Both calls pass the orbital centres as ``site_positions``
    (``orbital_centres``, ``(norb_phys, 3)`` fractional, from
    :func:`orbital_centres`); a site group sits at the centre of its first
    orbital (its members lie within ``site_tol``, which the reader's 1e-8
    centre check would refuse). ``base_dir`` is passed to the reader (a
    collected WORK passes the run directory, so no run record is consulted).

    Returns ``{"per_orbital", "per_site", "per_site_reason"}`` (JSON)."""
    from hwave.results import read_uhfk_order_parameters
    centres = np.asarray(orbital_centres, dtype=float)
    if centres.ndim != 2 or centres.shape[1] != 3 or not np.all(np.isfinite(centres)):
        raise ValueError("orbital_centres: a finite (norb_phys, 3) array, got shape {}"
                         .format(centres.shape))
    n = centres.shape[0]
    per_orbital = read_uhfk_order_parameters(
        run_or_toml, site_of_orbital=np.arange(n), site_positions=centres,
        base_dir=base_dir)
    if isinstance(site_groups, SiteGroups):
        source, groups, ambiguous = (site_groups.source, site_groups.groups,
                                     site_groups.ambiguous)
    else:
        source = site_groups["source"]
        groups, ambiguous = site_groups["groups"], site_groups["ambiguous"]
    out = {"per_orbital": _op_json(per_orbital), "per_site": None, "per_site_reason": None}
    if ambiguous and source == "clustered":
        out["per_site_reason"] = ("the site grouping is ambiguous ({}); give the grouping "
                                  "with orbital_site for per-site order parameters"
                                  .format("; ".join(ambiguous)))
        return out
    sof = np.full(n, -1, dtype=int)
    for g, members in enumerate(groups):
        for a in members:
            sof[int(a)] = g
    if np.any(sof < 0):
        raise ValueError("site_groups: orbitals {} are in no group".format(
            np.flatnonzero(sof < 0).tolist()))
    positions = np.array([centres[int(min(m))] for m in groups])
    per_site = read_uhfk_order_parameters(
        run_or_toml, site_of_orbital=sof, site_positions=positions, base_dir=base_dir)
    out["per_site"] = _op_json(per_site)
    return out


# --------------------------------------------------------------------------
# the summary
# --------------------------------------------------------------------------

def _num(v, fmt="{:.8g}"):
    if v is None:
        return "-"
    if isinstance(v, float) and not math.isfinite(v):
        return str(v)
    return fmt.format(v) if isinstance(v, (int, float)) else str(v)


def _table(header, rows):
    rows = [[str(c) for c in r] for r in rows]
    widths = [max(len(h), *(len(r[i]) for r in rows)) if rows else len(h)
              for i, h in enumerate(header)]
    line = lambda cells_: "  " + " | ".join(c.ljust(w) for c, w in zip(cells_, widths))
    return [line(header), "  " + "-+-".join("-" * w for w in widths)] + \
        [line(r) for r in rows]


def _reason_key(reason):
    """The not-tried count key: the full reason, with the referent and map of
    a duplicate left out (``"duplicate"``, and ``"quota: <family>
    (duplicate of a row cut by its quota)"`` for a promoted duplicate)."""
    if reason.startswith("duplicate of "):
        return "duplicate"
    if reason.startswith("quota:") and "(duplicate of " in reason:
        return reason.split(" (duplicate of ")[0] + " (duplicate of a row cut by its quota)"
    return reason


def _near_clusters(lst):
    """The near-degenerate entries of a list as clusters (connected by the
    near-degeneracy relation), each once, in rank order."""
    ids = [e.run_id for e in lst.entries]
    adj = {e.run_id: set(e.near_degenerate_with) for e in lst.entries}
    seen, out = set(), []
    for rid in ids:
        if rid in seen or not adj[rid]:
            continue
        stack, comp = [rid], set()
        while stack:
            x = stack.pop()
            if x in comp:
                continue
            comp.add(x)
            stack.extend(adj.get(x, ()))
        seen |= comp
        out.append([i for i in ids if i in comp])
    return out


def _rank_rows(lst):
    out = []
    for e in lst.entries:
        d = e.delta_per_site if lst.degeneracy_scale == "per site" else e.delta_per_cell
        out.append([e.run_id, e.pattern_label or "-", _num(e.free_energy_per_cell),
                    _num(d, "{:.3e}"), e.status, e.final_label or "-",
                    ", ".join(e.flags) or "-"])
    return out


def summary_text(doc, ranking=None) -> str:
    """The summary text: the model and cost line; the table
    (id, label, F per cell, dF per site, status, final label); the flags;
    the structural exclusions; the fixed sentence (always the last line).
    ``ranking``: a :class:`Ranking` (post-run); without one the tried rows
    are listed as planned runs."""
    m, tier, lim = doc["model"], doc["tier"], doc["limits"]
    hint = doc.get("hint") or {}
    out = ["Candidate initial states: coverage report ({}, preview)".format(
        doc.get("phase_marker"))]
    su2 = (m.get("su2") or {}).get("status")
    out.append("Model: norb_phys = {}, n_c = {}, spin-orbital {}, Fock {}, Ncond = {} "
               "({} per cell), T = {}, 2Sz = {}, CellShape {}, SU(2) check {}".format(
                   m.get("norb_phys"), m.get("n_c"),
                   "on" if m.get("spin_orbital") else "off",
                   "on" if m.get("iflag_fock") else "off", m.get("Ncond"),
                   m.get("Ncond_per_cell"), m.get("T"), _sz_text(_canon_sz(m.get("2Sz"))),
                   m.get("cellshape"), su2))
    est = tier.get("estimate") or {}
    cost = ("no cost estimate in this release" if est.get("bytes") is None
            and est.get("seconds") is None else "estimate {} bytes, {} s ({})".format(
                est.get("bytes"), est.get("seconds"), est.get("class")))
    out.append("Tier: {} (chosen by {}); {}".format(tier.get("name"), tier.get("chosen_by"),
                                                     cost))
    th = lim.get("T_hint") or {}
    out.append("T_hint = {} ({}); hints {}{}".format(
        _num(th.get("value")), th.get("rule"), hint.get("status"),
        ": " + hint["skip_reason"] if hint.get("skip_reason") else ""))
    lt = lim.get("label_thresholds") or {}
    out.append("Labels name the components holding >= {:.0%} of the weight; TR even/odd "
               "when the other parity holds < {:.0%}, else mixed".format(
                   lt.get("dominant_weight", 0.1), lt.get("tr_minor_weight", 0.1)))
    rows = doc.get("rows") or []
    tried = [r for r in rows if r.get("tried")]
    reasons = {}
    for r in rows:
        if not r.get("tried"):
            key = _reason_key(r.get("not_tried_reason") or "")
            reasons[key] = reasons.get(key, 0) + 1
    out.append("Rows: {} ({} tried); not tried: {}".format(
        len(rows), len(tried), ", ".join("{} {}".format(n, k) for k, n in
                                         sorted(reasons.items())) or "none"))
    out.append("")
    header = ["id", "label", "F per cell", "dF per site", "status", "final label", "flags"]
    flags = []
    if ranking is None:
        out.append("Planned runs:")
        out.extend(_table(["id", "label", "SubShape"],
                          [[r["run_id"], r.get("pattern_label") or "-", r.get("subshape")]
                           for r in tried]))
    else:
        for i, lst in enumerate(ranking.lists):
            g = lst.group or {}
            out.append("Ranking {} (converged; T = {}, Ncond per cell {}, CellShape {}, "
                       "2Sz {}; near-degeneracy judged {}):".format(
                           i + 1, g.get("T"), g.get("Ncond_per_cell"), g.get("cellshape"),
                           _sz_text(g.get("2Sz")), lst.degeneracy_scale))
            hdr = list(header)
            if lst.degeneracy_scale == "per cell":
                hdr[3] = "dF per cell"
            out.extend(_table(hdr, _rank_rows(lst)))
            out.append("")
        if ranking.flagged_unverified.entries:
            out.append("Flagged: convergence unverified (ranked on the same key, kept "
                       "apart):")
            out.extend(_table(header, _rank_rows(ranking.flagged_unverified)))
            out.append("")
        if ranking.unranked:
            out.append("Not ranked (listed, never dropped):")
            out.extend(_table(["id", "label", "status", "iterations", "final label",
                               "reason"],
                              [[r["run_id"], r.get("pattern_label") or "-",
                                r.get("status") or "-", _num(r.get("iterations")),
                                r.get("final_label") or "-",
                                r.get("unranked_reason")] for r in ranking.unranked]))
            out.append("")
        for lst in tuple(ranking.lists) + (ranking.flagged_unverified,):
            for e in lst.entries:
                if "cap hit" in e.flags:
                    flags.append("cap hit: {} (the seed field reached the cap without the "
                                 "target amplitude)".format(e.run_id))
            for cluster in _near_clusters(lst):
                flags.append("near-degenerate: {} (|dF| {} < {:g} in input energy units)"
                             .format(", ".join(cluster), lst.degeneracy_scale,
                                     ranking.degeneracy_tol))
        flags.extend(ranking.warnings)
    for qb in lim.get("quota_boundary") or ():
        counts = "; ".join("{} {} of {}".format(q, c["tried"], c["tried"] + c["cut"])
                           for q, c in qb["counts"].items())
        flags.append("quota reached inside q-order n = {} on SubShape {}: rows tried per q: "
                     "{}; lattice-equivalent q may be treated unequally (no point-group "
                     "maps are used)".format(qb["order"], qb["subshape"], counts))
    warns = list(m.get("warnings") or ())
    out.append("Flags:" if flags or warns else "Flags: none")
    for w in warns:
        prefix = "WARNING: " if w.startswith(SZ_WARNING) else "warning: "
        out.append("  " + prefix + w)
    for fl in flags:
        out.append("  " + fl)
    out.append("")
    out.append("Limits: UHF overestimates order and ordering temperatures; the k mesh and "
               "the cells are finite; the hints are static, for the paramagnet at T_hint; "
               "the result is the lowest among the candidates tried.")
    out.append("Structural exclusions (no candidate covers these):")
    for i, item in enumerate(lim.get("structural_exclusions") or STRUCTURAL_EXCLUSIONS, 1):
        out.append("  ({}) {}".format(i, item))
    out.append("")
    out.append(FIXED_SENTENCE)
    return "\n".join(out) + "\n"
