"""The RPA/FLEX equivalence-table cell registry.

This module is the SIDE-EFFECT-FREE typed record registry described in
``docs/superpowers/specs/2026-08-18-equivalence-table-design.md`` (Step
4 of the RPA/FLEX consolidation sequence) and
``docs/superpowers/plans/2026-08-18-equivalence-table.md`` (Task 1 of
that plan). It defines the record types every later task builds on
(comparators, executors, the 38-cell inventory, and a rendered docs
page) plus the schema validators that keep the registry internally
consistent.

Hard constraints (binding on every future edit to this file):

  * NO ``unittest`` import here, and no import of any solver code. The
    module must be importable by both the generated-test module and the
    docs renderer without running anything. Validation TESTS (including
    ``SupplementaryLink`` existence/discoverability) live in
    ``tests/test_rpa_flex_equivalence_table.py``; this module only
    checks record SHAPE.
  * ``from __future__ import annotations`` -- Python 3.9 compatibility,
    no runtime ``X | None`` annotations.
  * Every dataclass here is ``frozen``. Nested mappings (``dict``-typed
    fields) are normalized to ``types.MappingProxyType`` in
    ``__post_init__`` so importers (tests, the renderer) cannot mutate a
    shared record in place.
  * Dataclass constructors never raise on BUSINESS-RULE violations (e.g.
    an unknown ``ExecuteReject.exc_type``, a mismatched
    ``required_observables`` bundle). Those invariants are reported by
    ``validate_registry`` so tests can construct intentionally-invalid
    records and assert on the returned error list. Construction only
    normalizes mapping fields; it never inspects cross-field business
    rules.

``CELLS`` and ``COVERAGE_OBLIGATIONS`` are intentionally EMPTY in this
task -- this module ships the machinery (types + validators), not the
38-cell inventory (Task 5) or its coverage predicates.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from enum import Enum
from types import MappingProxyType
from typing import Callable, Dict, List, Mapping, NamedTuple, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Enums
# ---------------------------------------------------------------------------


class Status(Enum):
    """Per-solver capability outcome for a cell."""

    SUPPORTED = "supported"
    REJECT = "reject"
    NOT_APPLICABLE = "not_applicable"


class Site(Enum):
    """Where an ``ExecuteReject`` proof's public call raises."""

    CONSTRUCTOR = "constructor"
    SOLVE = "solve"


# ---------------------------------------------------------------------------
# Proof steps -- the CLOSED PROOF_STEP_TYPES union
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ExecuteRun:
    """The cell's test runs the solver's ``solve()`` successfully.

    The single run-class proof for an ordinary SUPPORTED cell that also
    discharges a comparison (the run-once-per-cell rule).
    """


@dataclass(frozen=True)
class ExecuteConstruct:
    """The cell's test constructs the solver and asserts the resolved
    ``calc_scheme`` -- used for construction-only ``auto`` resolution
    rows (no solve; ``expected_resolved_scheme`` names the scheme
    ``_set_scheme`` is expected to resolve to).
    """

    expected_resolved_scheme: str


# The exc_type CLOSED set. Membership is validator-enforced (see
# ``validate_registry``), never checked at construction time.
_CLOSED_EXC_TYPES: Tuple[str, ...] = ("ValueError", "RuntimeError")


@dataclass(frozen=True)
class ExecuteReject:
    """The cell's test asserts the named public call (constructor |
    solve) raises ``exc_type`` with a stable, project-owned ``fragment``
    substring in its message.

    ``exc_type`` must be a member of the CLOSED set
    ``{"ValueError", "RuntimeError"}`` -- validator-enforced by
    ``validate_registry``; the executor (a later task) maps it to the
    concrete exception class via a fixed dict keyed on this same set.
    """

    site: Site
    exc_type: str
    fragment: str


@dataclass(frozen=True)
class PairedInvarianceRun:
    """Two solves of the SAME solver (documented exception to the
    run-once rule) asserting an invariance property -- e.g. the FLEX
    side of ``chi0q_init`` reuse: with vs. without the option, exhaustive
    key-set + per-key ``np.array_equal`` over ``green_info``.
    """


@dataclass(frozen=True)
class ExecuteChiqInitReuse:
    """Two RPA solves -- the ``chi0q_init`` reuse oracle (cell 33): a
    first run without ``chi0q_init`` capturing ``chi0q_A``/``chiq_A``,
    then a second run consuming ``chi0q_A`` via the public
    ``chi0q_init`` mechanism, asserting bitwise-identical
    ``chi0q_B``/``chiq_B``.
    """


# The CLOSED union. ``validate_registry`` rejects any step object whose
# type is not a member of this tuple -- unknown ProofStep types can
# never silently enter the registry.
PROOF_STEP_TYPES: Tuple[type, ...] = (
    ExecuteRun,
    ExecuteConstruct,
    ExecuteReject,
    PairedInvarianceRun,
    ExecuteChiqInitReuse,
)

# The run-class steps that discharge a SUPPORTED status.
_RUN_CLASS_STEP_TYPES: Tuple[type, ...] = (ExecuteRun, ExecuteConstruct, ExecuteChiqInitReuse)

# The run-class steps that additionally imply a comparison is expected
# (ExecuteConstruct is deliberately excluded -- construction-only rows
# never carry a comparison).
_SOLVE_RUNNING_STEP_TYPES: Tuple[type, ...] = (ExecuteRun, ExecuteChiqInitReuse)


# ---------------------------------------------------------------------------
# Supplementary links (documentation only -- zero proof authority)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SupplementaryLink:
    """A non-proof documentation link to an existing test.

    ``test_id`` must have the shape
    ``<dotted.module>::<TestCase>::<method>`` (schema-validated FORMAT
    only, here); existence/discoverability of the referenced module,
    class, and method is a TEST-module concern
    (``TestSupplementaryLinks``, via ``unittest.TestLoader``), not this
    registry's.
    """

    test_id: str
    claim: str


# ---------------------------------------------------------------------------
# Per-solver proof
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SolverProof:
    """One solver side's (RPA or FLEX) proof for a cell.

    Status<->step composition (EXHAUSTIVE, validator-enforced -- see
    ``validate_registry``):

      * SUPPORTED: exactly ONE run-class step
        (``ExecuteRun`` | ``ExecuteConstruct`` | ``ExecuteChiqInitReuse``)
        and NO other step types.
      * REJECT: exactly ONE ``ExecuteReject`` and NOTHING else.
      * NOT_APPLICABLE: exactly one ``PairedInvarianceRun`` OR zero
        steps (the narrowed claim -- then ``reason`` must contain the
        substring ``"no corresponding"``).

    ``reason`` is mandatory IFF ``status is Status.NOT_APPLICABLE``
    (biconditional: any other status must leave it empty).
    """

    status: Status
    steps: tuple
    links: tuple = ()
    reason: str = ""


# ---------------------------------------------------------------------------
# Observable specs and the ComparisonProof sum type (Equiv | Diverges)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ObservableSpec:
    """A single observable's comparator + checked-in tolerance."""

    comparator: str
    atol: float
    provenance: str


@dataclass(frozen=True)
class DivergingSpec:
    """A single observable that is classified Diverges: a two-sided
    asserted bracket ``(policy_ceiling, regression_bound]`` -- strict
    lower edge (``ceiling < residual``, so equality never classifies as
    Diverges), inclusive upper (``residual <= regression_bound``).
    """

    comparator: str
    ceiling: float
    regression_bound: float
    provenance: str


def _freeze_mapping(mapping) -> MappingProxyType:
    """Normalize a plain mapping into an immutable
    ``types.MappingProxyType`` view, so records cannot be mutated in
    place by importers (tests, the renderer). Idempotent on an
    already-frozen mapping.
    """

    if isinstance(mapping, MappingProxyType):
        return mapping
    return MappingProxyType(dict(mapping))


@dataclass(frozen=True)
class Equiv:
    """The ComparisonProof variant for a (SUPPORTED, SUPPORTED) cell
    whose every required observable meets its policy ceiling.

    Schema invariant (validator-enforced): ``observables.keys()`` must
    equal ``Cell.required_observables`` exactly, and every entry's
    ``atol`` must be <= its mapped ``POLICY_CEILINGS`` entry (equality
    allowed).
    """

    observables: Mapping[str, ObservableSpec]

    def __post_init__(self) -> None:
        object.__setattr__(self, "observables", _freeze_mapping(self.observables))


@dataclass(frozen=True)
class Diverges:
    """The ComparisonProof variant for a (SUPPORTED, SUPPORTED) cell
    with at least one observable above its policy ceiling.

    Schema invariants (validator-enforced): ``diverging`` is non-empty;
    ``diverging`` and ``others`` are disjoint; their union equals
    ``Cell.required_observables`` exactly; every ``DivergingSpec.ceiling``
    equals its mapped ``POLICY_CEILINGS`` entry, with
    ``ceiling < regression_bound``; every ``others`` entry obeys the
    ordinary ``Equiv``-style ceiling rule; ``step5_issue`` matches
    ``^#[1-9][0-9]*$`` (a real issue number -- placeholders never enter
    the registry).
    """

    diverging: Mapping[str, DivergingSpec]
    others: Mapping[str, ObservableSpec]
    step5_issue: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "diverging", _freeze_mapping(self.diverging))
        object.__setattr__(self, "others", _freeze_mapping(self.others))


# ---------------------------------------------------------------------------
# Fixtures and cells
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FixtureSpec:
    """An immutable input description for one cell.

    ``mu`` and ``filling`` are mutually exclusive (validator-enforced):
    fixed-mu cells set ``mu``; mu-coupled cells set ``filling``.
    ``extern`` names a committed Zeeman/Extern file in ``input_dir``
    (the public spin-diag route); ``geom``/``transfer`` override the
    default ``geom.dat``/``transfer.dat`` names for the SO family
    (``geom_so.dat``/``transfer_spinful.dat``).
    """

    input_dir: str
    interactions: Mapping[str, str]
    T: float
    mu: Optional[float]
    filling: Optional[float]
    CellShape: tuple
    SubShape: tuple
    Nmat: int
    extra_params: Mapping
    calc_type: str
    requested_scheme: str
    enable_spin_orbital: bool
    extern: Optional[str]
    geom: str = "geom.dat"
    transfer: str = "transfer.dat"

    def __post_init__(self) -> None:
        object.__setattr__(self, "interactions", _freeze_mapping(self.interactions))
        object.__setattr__(self, "extra_params", _freeze_mapping(self.extra_params))


@dataclass(frozen=True)
class Cell:
    """One row of the equivalence table."""

    cell_id: str
    fixture: FixtureSpec
    resolved_scheme: str
    expected_spin_mode: str
    rpa: SolverProof
    flex: SolverProof
    comparison: object  # Equiv | Diverges | None
    required_observables: tuple
    interaction_class: str
    notes: str


# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

# The Global-Constraints values, verbatim.
POLICY_CEILINGS: dict = {
    "mu_diag": 1e-10,
    "green_diag": 1e-10,
    "chi0q_mu": 1e-10,
    "chiq_mu": 1e-10,
    "chi0q_fixed": 1e-12,
    "chiq_fixed": 1e-12,
}

# Schema-validated at creation ({"source_sha": None, "run_ids": (),
# "status": "candidate"}); the freeze commit's ONLY registry change.
PROVENANCE: dict = {
    "source_sha": None,
    "run_ids": (),
    "status": "candidate",
}

# The 38-cell inventory -- populated by Task 5. Empty this task.
CELLS: tuple = ()

# Named predicates over ``CELLS`` (or any cells tuple passed to
# ``validate_registry``) -- populated by Task 5. Empty this task.
COVERAGE_OBLIGATIONS: Dict[str, Callable[[tuple], bool]] = {}


# ---------------------------------------------------------------------------
# method_name
# ---------------------------------------------------------------------------

_METHOD_NAME_TRANSLATION = {ord("."): "_", ord("+"): "_"}


def method_name(cell_id: str) -> str:
    """``"test_cell__" + cell_id`` with ``.`` and ``+`` mapped to ``_``.

    A registry test asserts the mapped names are unique AND valid
    Python identifiers across the whole registry (a duplicate would let
    ``setattr`` silently overwrite a generated proof method in Task 3).
    """

    return "test_cell__" + cell_id.translate(_METHOD_NAME_TRANSLATION)


# ---------------------------------------------------------------------------
# relationship() -- the spec's TOTAL status-pair lookup
# ---------------------------------------------------------------------------

# | rpa \ flex      | SUPPORTED                      | REJECT                          | NOT_APPLICABLE |
# | SUPPORTED       | see comparison outcome         | FLEX-REJECT.RPA-SUPPORTED       | RPA-ONLY        |
# | REJECT          | RPA-REJECT.FLEX-SUPPORTED      | BOTH-REJECT                     | RPA-REJECT      |
# | NOT_APPLICABLE  | FLEX-ONLY                      | FLEX-REJECT                     | N/A (no cell)   |
_STATUS_PAIR_RELATIONSHIP: Dict[Tuple[Status, Status], str] = {
    (Status.SUPPORTED, Status.REJECT): "FLEX-REJECT·RPA-SUPPORTED",
    (Status.SUPPORTED, Status.NOT_APPLICABLE): "RPA-ONLY",
    (Status.REJECT, Status.SUPPORTED): "RPA-REJECT·FLEX-SUPPORTED",
    (Status.REJECT, Status.REJECT): "BOTH-REJECT",
    (Status.REJECT, Status.NOT_APPLICABLE): "RPA-REJECT",
    (Status.NOT_APPLICABLE, Status.SUPPORTED): "FLEX-ONLY",
    (Status.NOT_APPLICABLE, Status.REJECT): "FLEX-REJECT",
}


def _construction_step(proof: SolverProof) -> Optional[ExecuteConstruct]:
    """Return the ``ExecuteConstruct`` step if ``proof`` is a
    single-step SUPPORTED construction-only proof, else ``None``.
    """

    if (
        proof.status is Status.SUPPORTED
        and len(proof.steps) == 1
        and isinstance(proof.steps[0], ExecuteConstruct)
    ):
        return proof.steps[0]
    return None


def relationship(cell: Cell) -> str:
    """The spec's TOTAL lookup over ``(rpa.status, flex.status)``, plus
    the construction-only variant.

    Construction-only rows (both sides carry a single ``ExecuteConstruct``
    step with matching ``expected_resolved_scheme`` values, equal to
    ``cell.resolved_scheme``) return ``"AUTO-RESOLVES(<scheme>)"`` --
    the status-pair lookup below applies to every other row.

    ``(SUPPORTED, SUPPORTED)`` (non-construction) rows resolve via the
    cell's comparison outcome: ``"EQUIV"`` for ``Equiv``,
    ``"DIVERGES(<step5_issue>)"`` for ``Diverges``.

    ``(NOT_APPLICABLE, NOT_APPLICABLE)`` has no relationship: the model
    declares this status pair nonexistent, so this raises ``ValueError``
    (mirrors ``validate_registry``'s outright rejection of the same
    pair).
    """

    rpa_construct = _construction_step(cell.rpa)
    flex_construct = _construction_step(cell.flex)
    if rpa_construct is not None and flex_construct is not None:
        if (
            rpa_construct.expected_resolved_scheme == flex_construct.expected_resolved_scheme
            and rpa_construct.expected_resolved_scheme == cell.resolved_scheme
        ):
            return "AUTO-RESOLVES({})".format(cell.resolved_scheme)

    rpa_status, flex_status = cell.rpa.status, cell.flex.status

    if rpa_status is Status.SUPPORTED and flex_status is Status.SUPPORTED:
        if isinstance(cell.comparison, Equiv):
            return "EQUIV"
        if isinstance(cell.comparison, Diverges):
            return "DIVERGES({})".format(cell.comparison.step5_issue)
        raise ValueError(
            "cell {!r}: (SUPPORTED, SUPPORTED) without a matching "
            "construction-only pair requires a comparison outcome "
            "(Equiv or Diverges), got {!r}".format(cell.cell_id, cell.comparison)
        )

    if rpa_status is Status.NOT_APPLICABLE and flex_status is Status.NOT_APPLICABLE:
        raise ValueError(
            "cell {!r}: (NOT_APPLICABLE, NOT_APPLICABLE) has no "
            "relationship -- the model declares this status pair "
            "nonexistent".format(cell.cell_id)
        )

    return _STATUS_PAIR_RELATIONSHIP[(rpa_status, flex_status)]


# ---------------------------------------------------------------------------
# validate_registry -- the schema validators
# ---------------------------------------------------------------------------

_LINK_RE = re.compile(
    r"^[A-Za-z_][A-Za-z0-9_]*(?:\.[A-Za-z_][A-Za-z0-9_]*)*"
    r"::[A-Za-z_][A-Za-z0-9_]*::[A-Za-z_][A-Za-z0-9_]*$"
)
_ISSUE_RE = re.compile(r"^#[1-9][0-9]*$")


def _ceiling_key(observable: str, fixture: FixtureSpec) -> Optional[str]:
    """The explicit observable -> ``POLICY_CEILINGS`` key mapping,
    keyed on observable name + the fixture's mu-mode: mu-coupled
    fixtures (``filling`` set) map to the ``*_mu`` ceilings; fixed-mu
    fixtures (``mu`` set) map to the ``*_fixed`` ceilings. Returns
    ``None`` when the fixture sets neither or both (the mu/filling
    exclusivity check reports that separately).
    """

    if fixture.filling is not None and fixture.mu is None:
        return "{}_mu".format(observable)
    if fixture.mu is not None and fixture.filling is None:
        return "{}_fixed".format(observable)
    return None


def _validate_side(cell_id: str, side_name: str, proof: SolverProof) -> List[str]:
    errors: List[str] = []

    for step in proof.steps:
        if not isinstance(step, PROOF_STEP_TYPES):
            errors.append(
                "cell {!r} {}: step {!r} is not a member of the closed "
                "PROOF_STEP_TYPES union".format(cell_id, side_name, step)
            )
    if errors:
        # A malformed step makes the composition checks below meaningless.
        return errors

    if proof.status is Status.NOT_APPLICABLE:
        if not proof.reason:
            errors.append(
                "cell {!r} {}: reason is mandatory when status is "
                "NOT_APPLICABLE".format(cell_id, side_name)
            )
    elif proof.reason:
        errors.append(
            "cell {!r} {}: reason must be empty unless status is "
            "NOT_APPLICABLE, got status={!r} reason={!r}".format(
                cell_id, side_name, proof.status, proof.reason
            )
        )

    if proof.status is Status.SUPPORTED:
        if len(proof.steps) != 1 or not isinstance(proof.steps[0], _RUN_CLASS_STEP_TYPES):
            errors.append(
                "cell {!r} {}: SUPPORTED requires exactly one run-class "
                "step (ExecuteRun | ExecuteConstruct | ExecuteChiqInitReuse) "
                "and no other step types, got {!r}".format(cell_id, side_name, proof.steps)
            )
    elif proof.status is Status.REJECT:
        if len(proof.steps) != 1 or not isinstance(proof.steps[0], ExecuteReject):
            errors.append(
                "cell {!r} {}: REJECT requires exactly one ExecuteReject "
                "step and nothing else, got {!r}".format(cell_id, side_name, proof.steps)
            )
        else:
            exc_type = proof.steps[0].exc_type
            if exc_type not in _CLOSED_EXC_TYPES:
                errors.append(
                    "cell {!r} {}: ExecuteReject.exc_type {!r} is not in "
                    "the closed set {!r}".format(cell_id, side_name, exc_type, _CLOSED_EXC_TYPES)
                )
    elif proof.status is Status.NOT_APPLICABLE:
        if len(proof.steps) == 1 and isinstance(proof.steps[0], PairedInvarianceRun):
            pass
        elif len(proof.steps) == 0:
            if "no corresponding" not in proof.reason:
                errors.append(
                    "cell {!r} {}: NOT_APPLICABLE with zero steps (the "
                    "narrowed claim) requires reason to contain "
                    "'no corresponding', got {!r}".format(cell_id, side_name, proof.reason)
                )
        else:
            errors.append(
                "cell {!r} {}: NOT_APPLICABLE requires exactly one "
                "PairedInvarianceRun step or zero steps, got {!r}".format(
                    cell_id, side_name, proof.steps
                )
            )
    else:
        errors.append(
            "cell {!r} {}: unknown status {!r}".format(cell_id, side_name, proof.status)
        )

    errors.extend(_validate_links(cell_id, side_name, proof.links))
    return errors


def _validate_links(cell_id: str, side_name: str, links: tuple) -> List[str]:
    errors: List[str] = []
    seen = set()
    for link in links:
        if not _LINK_RE.match(link.test_id):
            errors.append(
                "cell {!r} {}: SupplementaryLink.test_id {!r} does not "
                "match '<dotted.module>::<TestCase>::<method>'".format(
                    cell_id, side_name, link.test_id
                )
            )
        if link.test_id in seen:
            errors.append(
                "cell {!r} {}: duplicate SupplementaryLink.test_id {!r}".format(
                    cell_id, side_name, link.test_id
                )
            )
        seen.add(link.test_id)
    return errors


def _validate_equiv(cell: Cell) -> List[str]:
    errors: List[str] = []
    equiv = cell.comparison
    observed_keys = set(equiv.observables.keys())
    required_keys = set(cell.required_observables)
    if observed_keys != required_keys:
        errors.append(
            "cell {!r}: Equiv.observables keys {!r} must equal "
            "required_observables {!r} exactly".format(
                cell.cell_id, sorted(observed_keys), sorted(required_keys)
            )
        )
    for name, spec in equiv.observables.items():
        key = _ceiling_key(name, cell.fixture)
        if key is None or key not in POLICY_CEILINGS:
            errors.append(
                "cell {!r}: observable {!r} has no explicit "
                "POLICY_CEILINGS mapping for this fixture's mu-mode "
                "(mapped key {!r})".format(cell.cell_id, name, key)
            )
            continue
        ceiling = POLICY_CEILINGS[key]
        if spec.atol > ceiling:
            errors.append(
                "cell {!r}: observable {!r} atol {!r} exceeds its policy "
                "ceiling {!r} ({!r})".format(cell.cell_id, name, spec.atol, ceiling, key)
            )
    return errors


def _validate_diverges(cell: Cell) -> List[str]:
    errors: List[str] = []
    diverges = cell.comparison
    diverging_keys = set(diverges.diverging.keys())
    others_keys = set(diverges.others.keys())
    required_keys = set(cell.required_observables)

    if not diverging_keys:
        errors.append("cell {!r}: Diverges.diverging must be non-empty".format(cell.cell_id))
    overlap = diverging_keys & others_keys
    if overlap:
        errors.append(
            "cell {!r}: Diverges.diverging and .others must be disjoint, "
            "overlap {!r}".format(cell.cell_id, sorted(overlap))
        )
    if (diverging_keys | others_keys) != required_keys:
        errors.append(
            "cell {!r}: Diverges.diverging union others {!r} must equal "
            "required_observables {!r} exactly".format(
                cell.cell_id, sorted(diverging_keys | others_keys), sorted(required_keys)
            )
        )
    if not _ISSUE_RE.match(diverges.step5_issue):
        errors.append(
            "cell {!r}: Diverges.step5_issue {!r} does not match "
            "'^#[1-9][0-9]*$'".format(cell.cell_id, diverges.step5_issue)
        )

    for name, spec in diverges.diverging.items():
        key = _ceiling_key(name, cell.fixture)
        if key is None or key not in POLICY_CEILINGS:
            errors.append(
                "cell {!r}: diverging observable {!r} has no explicit "
                "POLICY_CEILINGS mapping for this fixture's mu-mode "
                "(mapped key {!r})".format(cell.cell_id, name, key)
            )
            continue
        ceiling = POLICY_CEILINGS[key]
        if spec.ceiling != ceiling:
            errors.append(
                "cell {!r}: DivergingSpec.ceiling {!r} for observable "
                "{!r} must equal the mapped policy ceiling {!r} "
                "({!r})".format(cell.cell_id, spec.ceiling, name, ceiling, key)
            )
        if not (spec.ceiling < spec.regression_bound):
            errors.append(
                "cell {!r}: DivergingSpec for observable {!r} must "
                "satisfy ceiling < regression_bound, got ceiling={!r} "
                "regression_bound={!r}".format(
                    cell.cell_id, name, spec.ceiling, spec.regression_bound
                )
            )

    for name, spec in diverges.others.items():
        key = _ceiling_key(name, cell.fixture)
        if key is None or key not in POLICY_CEILINGS:
            errors.append(
                "cell {!r}: observable {!r} (Diverges.others) has no "
                "explicit POLICY_CEILINGS mapping for this fixture's "
                "mu-mode (mapped key {!r})".format(cell.cell_id, name, key)
            )
            continue
        ceiling = POLICY_CEILINGS[key]
        if spec.atol > ceiling:
            errors.append(
                "cell {!r}: observable {!r} (Diverges.others) atol {!r} "
                "exceeds its policy ceiling {!r} ({!r})".format(
                    cell.cell_id, name, spec.atol, ceiling, key
                )
            )
    return errors


def _validate_cell(cell: Cell) -> List[str]:
    errors: List[str] = []
    errors.extend(_validate_side(cell.cell_id, "rpa", cell.rpa))
    errors.extend(_validate_side(cell.cell_id, "flex", cell.flex))

    # The (NOT_APPLICABLE, NOT_APPLICABLE) status pair is rejected
    # outright -- the relationship model declares it nonexistent.
    if cell.rpa.status is Status.NOT_APPLICABLE and cell.flex.status is Status.NOT_APPLICABLE:
        errors.append(
            "cell {!r}: (NOT_APPLICABLE, NOT_APPLICABLE) is not a valid "
            "status pair -- the relationship model declares it "
            "nonexistent".format(cell.cell_id)
        )

    # mu/filling mutually exclusive.
    if cell.fixture.mu is not None and cell.fixture.filling is not None:
        errors.append(
            "cell {!r}: fixture.mu and fixture.filling are mutually "
            "exclusive, both are set (mu={!r}, filling={!r})".format(
                cell.cell_id, cell.fixture.mu, cell.fixture.filling
            )
        )

    # required_observables empty IFF comparison is None.
    has_required = len(cell.required_observables) > 0
    has_comparison = cell.comparison is not None
    if has_required != has_comparison:
        errors.append(
            "cell {!r}: required_observables must be empty iff "
            "comparison is None (required_observables={!r}, "
            "comparison={!r})".format(cell.cell_id, cell.required_observables, cell.comparison)
        )

    # comparison present IFF both sides SUPPORTED via a solve-running
    # step (ExecuteRun | ExecuteChiqInitReuse); None for
    # construction-only rows.
    rpa_solve_running = (
        cell.rpa.status is Status.SUPPORTED
        and len(cell.rpa.steps) == 1
        and isinstance(cell.rpa.steps[0], _SOLVE_RUNNING_STEP_TYPES)
    )
    flex_solve_running = (
        cell.flex.status is Status.SUPPORTED
        and len(cell.flex.steps) == 1
        and isinstance(cell.flex.steps[0], _SOLVE_RUNNING_STEP_TYPES)
    )
    should_have_comparison = rpa_solve_running and flex_solve_running
    if has_comparison != should_have_comparison:
        errors.append(
            "cell {!r}: comparison must be present iff both solver "
            "proofs are SUPPORTED via a solve-running step (ExecuteRun "
            "/ ExecuteChiqInitReuse); got comparison={!r}, rpa "
            "steps={!r}, flex steps={!r}".format(
                cell.cell_id, cell.comparison, cell.rpa.steps, cell.flex.steps
            )
        )

    # Construction-only rows: both sides ExecuteConstruct with equal
    # expected_resolved_scheme == Cell.resolved_scheme.
    rpa_construct = _construction_step(cell.rpa)
    flex_construct = _construction_step(cell.flex)
    if (rpa_construct is None) != (flex_construct is None):
        errors.append(
            "cell {!r}: construction-only rows require BOTH sides to "
            "carry a single ExecuteConstruct step; only one side "
            "does".format(cell.cell_id)
        )
    elif rpa_construct is not None and flex_construct is not None:
        if not (
            rpa_construct.expected_resolved_scheme
            == flex_construct.expected_resolved_scheme
            == cell.resolved_scheme
        ):
            errors.append(
                "cell {!r}: construction-only ExecuteConstruct."
                "expected_resolved_scheme must agree on both sides and "
                "equal Cell.resolved_scheme (rpa={!r}, flex={!r}, "
                "resolved_scheme={!r})".format(
                    cell.cell_id,
                    rpa_construct.expected_resolved_scheme,
                    flex_construct.expected_resolved_scheme,
                    cell.resolved_scheme,
                )
            )
        if cell.comparison is not None or cell.required_observables:
            errors.append(
                "cell {!r}: construction-only rows must have "
                "comparison=None and required_observables=()".format(cell.cell_id)
            )

    # ComparisonProof completeness + the observable -> POLICY_CEILINGS
    # mapping.
    if isinstance(cell.comparison, Equiv):
        errors.extend(_validate_equiv(cell))
    elif isinstance(cell.comparison, Diverges):
        errors.extend(_validate_diverges(cell))
    elif cell.comparison is not None:
        errors.append(
            "cell {!r}: comparison must be an Equiv, a Diverges, or "
            "None, got {!r}".format(cell.cell_id, type(cell.comparison))
        )

    return errors


def validate_registry(cells) -> list:
    """Validate a collection of ``Cell`` records; return a list of
    human-readable violation strings (empty iff the registry is
    internally consistent). Never raises on data violations -- every
    check reports into the returned list so callers (tests) can assert
    on the full set of violations for a deliberately-invalid registry.
    """

    cells = tuple(cells)
    errors: List[str] = []

    # UNIQUENESS: cell_ids unique.
    ids_seen: Dict[str, int] = {}
    for cell in cells:
        ids_seen[cell.cell_id] = ids_seen.get(cell.cell_id, 0) + 1
    for cell_id, count in ids_seen.items():
        if count > 1:
            errors.append(
                "duplicate cell_id {!r} ({} occurrences)".format(cell_id, count)
            )

    # UNIQUENESS: method_name outputs unique AND valid identifiers.
    methods_seen: Dict[str, List[str]] = {}
    for cell in cells:
        name = method_name(cell.cell_id)
        if not name.isidentifier():
            errors.append(
                "cell {!r}: method_name {!r} is not a valid Python "
                "identifier".format(cell.cell_id, name)
            )
        methods_seen.setdefault(name, []).append(cell.cell_id)
    for name, cell_ids in methods_seen.items():
        if len(cell_ids) > 1:
            errors.append(
                "duplicate method_name {!r} produced by cell_ids "
                "{!r}".format(name, cell_ids)
            )

    for cell in cells:
        errors.extend(_validate_cell(cell))

    # Coverage obligations (empty this task; a no-op until Task 5
    # populates COVERAGE_OBLIGATIONS).
    for name, predicate in COVERAGE_OBLIGATIONS.items():
        try:
            satisfied = predicate(cells)
        except Exception as exc:  # pragma: no cover - defensive
            errors.append(
                "coverage obligation {!r} raised {!r} while "
                "evaluating".format(name, exc)
            )
            continue
        if not satisfied:
            errors.append("coverage obligation {!r} is not satisfied".format(name))

    return errors


# ---------------------------------------------------------------------------
# Output bundles and the named comparator registry (Task 2)
# ---------------------------------------------------------------------------
#
# This section is still SIDE-EFFECT-FREE: ``numpy`` is the only import it
# adds (explicitly allowed by the plan's Task 2 interface -- "the
# comparators live in the REGISTRY module (pure numpy ... still no
# unittest/solver imports)"). ``extract_bundle`` -- the function that
# actually touches a solver object's ``green_info`` -- lives in
# ``tests/test_rpa_flex_equivalence_table.py``, not here.
#
# Comparator policy (Global Constraints, binding): elementwise
# ``abs(a - b) <= atol`` on complex values, NO rtol; every comparator
# checks exact shape then all-finite FIRST; on failure the diagnostic
# reports both the max-|diff| VALUE and its INDEX.


class OutputBundle(NamedTuple):
    """One solver side's post-solve observable bundle -- the input to
    every ``Comparator.map``.

    NO ``mu`` field (spec amendment, 2026-08-18, recorded in
    ``docs/superpowers/specs/2026-08-18-equivalence-table-design.md``):
    RPA never retains its solved chemical potential as a public
    post-solve attribute (it is a local variable of the solve) and
    FLEX's stored ``mu`` is the post-mix DRESSED value -- a different
    algorithmic state. The mu/Green seam is owned exclusively by the
    Task-6 divergence diagnostic via ``scalar_residual`` /
    ``assert_scalar_within`` below, never a per-cell ``COMPARATORS``
    entry.

    ``chiq`` is populated by RPA (both ``general`` and ``reduced``
    schemes, ``src/hwave/solver/rpa.py:2110``) and left ``None`` by
    FLEX, which never writes a combined ``"chiq"`` key. ``chiq_s`` /
    ``chiq_c`` are populated by FLEX (``src/hwave/solver/flex.py:
    746-748``) and left ``None`` by RPA.
    """

    chi0q: np.ndarray
    chiq: Optional[np.ndarray]
    chiq_s: Optional[np.ndarray]
    chiq_c: Optional[np.ndarray]


def _require_field(bundle: OutputBundle, field: str, side: str) -> np.ndarray:
    value = getattr(bundle, field)
    if value is None:
        raise ValueError(
            "comparator: {} bundle field {!r} is unset (None)".format(side, field)
        )
    return np.asarray(value)


def _map_identity(observable: str, rpa: OutputBundle, flex: OutputBundle) -> Tuple[np.ndarray, np.ndarray]:
    """chi0q from both bundles: the RPA and FLEX bare bubble is the
    identical Lindhard object -- no dressing at zeroth order, so the two
    solvers' ``green_info["chi0q"]`` must match under the fixed/mu
    comparator policy verbatim.
    """

    a = _require_field(rpa, observable, "rpa")
    b = _require_field(flex, observable, "flex")
    return a, b


def _map_general_from_flex_channels(
    observable: str, rpa: OutputBundle, flex: OutputBundle
) -> Tuple[np.ndarray, np.ndarray]:
    """RPA's ``chiq`` vs the FLEX channel same/diff reconstruction,
    copied VERBATIM from
    ``tests/test_rpa_flex_oneshot_equivalence.py:88-99``
    (``TestGeneralSchemeOneShot.test_multiorbital_onsite_cells_match_flex``)::

        chiq = np.asarray(gr['chiq'])
        cs = np.asarray(gf['chiq_s'])
        cc = np.asarray(gf['chiq_c'])
        recon = np.zeros_like(chiq)
        same = 0.5 * (cc + cs)
        diff = 0.5 * (cc - cs)
        for s1 in (0, 1):
            for s2 in (0, 1):
                blk = same if s1 == s2 else diff
                recon[:, :,
                      s1*norb:(s1+1)*norb, s1*norb:(s1+1)*norb,
                      s2*norb:(s2+1)*norb, s2*norb:(s2+1)*norb] = blk

    ``norb`` is read from ``flex.chiq_c``'s trailing axis (the FLEX
    channel arrays are already exactly ``norb``-sized on every axis --
    unlike RPA's ``chiq``, which carries the full ``2*norb`` spin-orbital
    axes this reconstruction rebuilds).
    """

    if observable != "chiq":
        raise ValueError(
            "general_from_flex_channels comparator only maps 'chiq', got "
            "observable={!r}".format(observable)
        )
    chiq = _require_field(rpa, "chiq", "rpa")
    cs = _require_field(flex, "chiq_s", "flex")
    cc = _require_field(flex, "chiq_c", "flex")

    norb = cc.shape[-1]
    recon = np.zeros_like(chiq)
    same = 0.5 * (cc + cs)
    diff = 0.5 * (cc - cs)
    for s1 in (0, 1):
        for s2 in (0, 1):
            blk = same if s1 == s2 else diff
            recon[
                :, :,
                s1 * norb:(s1 + 1) * norb, s1 * norb:(s1 + 1) * norb,
                s2 * norb:(s2 + 1) * norb, s2 * norb:(s2 + 1) * norb,
            ] = blk
    return chiq, recon


def _map_reduced_blocks(
    observable: str, rpa: OutputBundle, flex: OutputBundle
) -> Tuple[np.ndarray, np.ndarray]:
    """RPA's ``chiq`` uu/ud blocks vs FLEX's spin/charge channels,
    copied VERBATIM from ``tests/test_rpa_flex_oneshot_equivalence.py``'s
    ``TestReducedOneShot`` (``_blocks`` + ``test_matrix_cells``,
    lines 105-134)::

        def _blocks(self, cq, norb):
            uu = cq[:, :, :norb, :norb]
            ud = cq[:, :, :norb, norb:]
            return uu, ud
        ...
        uu, ud = self._blocks(np.asarray(gr['chiq']), norb)
        cs = np.asarray(gf['chiq_s'])[:, :, :norb, :norb]
        cc = np.asarray(gf['chiq_c'])[:, :, :norb, :norb]
        np.testing.assert_allclose(uu - ud, cs, rtol=0.0, atol=1e-12)
        np.testing.assert_allclose(uu + ud, cc, rtol=0.0, atol=1e-12)

    Both equations are folded into ONE elementwise comparison by
    stacking the two block pairs along a new leading axis: index 0 is
    the spin channel (``uu - ud`` vs ``cs``), index 1 is the charge
    channel (``uu + ud`` vs ``cc``). ``norb`` is read from ``rpa.chiq``'s
    trailing axis (``2 * norb`` -- the uu/ud split point), matching how
    the oneshot suite derives it from the fixture rather than from the
    (already ``norb``-sized) FLEX channel arrays.
    """

    if observable != "chiq":
        raise ValueError(
            "reduced_blocks comparator only maps 'chiq', got "
            "observable={!r}".format(observable)
        )
    cq = _require_field(rpa, "chiq", "rpa")
    cs_full = _require_field(flex, "chiq_s", "flex")
    cc_full = _require_field(flex, "chiq_c", "flex")

    if cq.shape[-1] % 2 != 0:
        raise ValueError(
            "reduced_blocks: rpa.chiq's trailing axis {!r} is not even "
            "(cannot split into uu/ud norb halves)".format(cq.shape[-1])
        )
    norb = cq.shape[-1] // 2
    uu = cq[:, :, :norb, :norb]
    ud = cq[:, :, :norb, norb:]
    cs = cs_full[:, :, :norb, :norb]
    cc = cc_full[:, :, :norb, :norb]

    a = np.stack([uu - ud, uu + ud])
    b = np.stack([cs, cc])
    return a, b


def _check_comparable(a: np.ndarray, b: np.ndarray) -> None:
    """Exact-shape then all-finite FIRST -- binding comparator policy
    order. Raises ``ValueError`` (a structural precondition failure, not
    a numeric-magnitude failure -- that distinction is what lets
    ``assert_within`` reserve ``AssertionError`` for genuine tolerance
    violations).
    """

    if a.shape != b.shape:
        raise ValueError(
            "comparator: shape mismatch {!r} vs {!r}".format(a.shape, b.shape)
        )
    if not np.all(np.isfinite(a)):
        raise ValueError("comparator: the first array contains non-finite values (NaN/Inf)")
    if not np.all(np.isfinite(b)):
        raise ValueError("comparator: the second array contains non-finite values (NaN/Inf)")


class Comparator:
    """A named, pure-numpy elementwise comparator.

    ``map`` extracts the two arrays to compare from a pair of
    ``OutputBundle`` records for one named ``observable``. ``residual``
    and ``assert_within`` are the SHARED comparison primitives every
    named comparator uses -- exact shape then all-finite first (see
    ``_check_comparable``), then the elementwise complex
    ``abs(a - b)`` max (comparator policy: NO rtol). ``assert_within``
    reports both the max-|diff| VALUE and its INDEX on failure.
    """

    def __init__(self, name: str, mapper: Callable[[str, OutputBundle, OutputBundle], Tuple[np.ndarray, np.ndarray]]) -> None:
        self.name = name
        self._mapper = mapper

    def map(self, observable: str, rpa: OutputBundle, flex: OutputBundle) -> Tuple[np.ndarray, np.ndarray]:
        return self._mapper(observable, rpa, flex)

    def residual(self, a, b) -> float:
        a = np.asarray(a)
        b = np.asarray(b)
        _check_comparable(a, b)
        return float(np.max(np.abs(a - b)))

    def assert_within(self, a, b, atol: float) -> None:
        a = np.asarray(a)
        b = np.asarray(b)
        _check_comparable(a, b)
        diff = np.abs(a - b)
        raw_idx = np.unravel_index(int(np.argmax(diff)), diff.shape)
        idx = tuple(int(i) for i in raw_idx)
        max_diff = float(diff[idx])
        if max_diff > atol:
            raise AssertionError(
                "{}: max |diff| = {!r} at index {!r} exceeds atol "
                "{!r}".format(self.name, max_diff, idx, atol)
            )


# The CLOSED comparator registry -- exactly the three keys the plan's
# Appendix A cell matrix references (no "scalar" key: the diagnostic's
# scalar mu checkpoints use ``scalar_residual``/``assert_scalar_within``
# below, never a ``COMPARATORS`` entry).
COMPARATORS: Dict[str, Comparator] = {
    "identity": Comparator("identity", _map_identity),
    "general_from_flex_channels": Comparator(
        "general_from_flex_channels", _map_general_from_flex_channels
    ),
    "reduced_blocks": Comparator("reduced_blocks", _map_reduced_blocks),
}


# ---------------------------------------------------------------------------
# Registry-level scalar utilities -- the Task-6 diagnostic's mu
# checkpoints (never a COMPARATORS entry: no per-cell mu, spec amendment).
# ---------------------------------------------------------------------------


def scalar_residual(a: float, b: float) -> float:
    """The scalar analogue of ``Comparator.residual``: ``abs(a - b)``,
    all-finite first. Used by the Task-6 divergence diagnostic's mu/
    Green checkpoints, which are NOT per-cell ``COMPARATORS`` entries.
    """

    if not (math.isfinite(a) and math.isfinite(b)):
        raise ValueError(
            "scalar_residual requires finite inputs, got a={!r} b={!r}".format(a, b)
        )
    return abs(a - b)


def assert_scalar_within(a: float, b: float, atol: float) -> None:
    """The scalar analogue of ``Comparator.assert_within``: raises
    ``AssertionError`` iff ``scalar_residual(a, b) > atol`` (equality
    passes -- comparator policy: ``abs(a - b) <= atol``).
    """

    residual = scalar_residual(a, b)
    if residual > atol:
        raise AssertionError(
            "scalar residual {!r} (a={!r}, b={!r}) exceeds atol "
            "{!r}".format(residual, a, b, atol)
        )


def assert_diverges_bracket(residual: float, ceiling: float, regression_bound: float) -> None:
    """Assert ``residual`` lies in the ``Diverges`` two-sided bracket
    ``(ceiling, regression_bound]`` -- the assertion-helper-level
    counterpart of the schema invariant ``DivergingSpec`` enforces
    structurally (``ceiling < regression_bound``). STRICT lower edge:
    ``residual > ceiling``, so equality with the ceiling never
    classifies as Diverges -- a residual that HEALED back to (or below)
    the ceiling means the cell no longer diverges and must be
    recalibrated back to ``Equiv`` (the spec's truth invariant,
    ``docs/superpowers/specs/2026-08-18-equivalence-table-design.md``);
    the raised message says so explicitly. INCLUSIVE upper edge:
    ``residual <= regression_bound``, a further-drift regression guard.
    """

    if residual <= ceiling:
        raise AssertionError(
            "residual {!r} healed to <= ceiling {!r}; recalibrate this "
            "cell back to Equiv".format(residual, ceiling)
        )
    if residual > regression_bound:
        raise AssertionError(
            "residual {!r} exceeds the regression_bound {!r}; this cell "
            "has regressed further than the calibrated bound".format(
                residual, regression_bound
            )
        )
