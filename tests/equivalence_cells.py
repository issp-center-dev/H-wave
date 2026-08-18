"""The RPA/FLEX equivalence-table cell registry.

This module is the SIDE-EFFECT-FREE typed record registry behind the
RPA/FLEX support and equivalence table. It defines the record types
everything else builds on (comparators, executors, the cell inventory,
and the rendered docs page) plus the schema validators that keep the
registry internally consistent.

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

``CELLS`` carries the full cell inventory this registry defines: rows
1-36 plus the conditioning row, cell 38 (there is no row 37, so the
registry holds 37 cells total). The conditioning row depends on the
divergence-diagnostic apparatus (``TestMuGreenDivergenceDiagnostic`` in
the test module): both its Equiv comparison and its coverage obligation
are stated in terms of that apparatus, so the two must be kept in step.
``COVERAGE_OBLIGATIONS`` is populated with every named coverage
predicate, including "the conditioning row exists (38)".

Registry docstring maintenance checklist: ``tests/
equivalence_benchmark.md`` records ONE per-cell timing row per cell in
``CELLS``, grouped into calibration sections (the development-machine
pass first, then one section per continuous-integration calibration) --
``TestBenchmarkRegistryTie`` (in the test module) asserts the benchmark
file's MOST RECENT section's cell_id column is an EXACT multiset match
against ``CELLS``. Any future edit
to this module that ADDS a cell, REMOVES a cell, RENAMES a
``cell_id``, or ENLARGES an existing cell's fixture (a change that
could plausibly move its measured wall time) therefore also requires
appending a NEW section to ``tests/equivalence_benchmark.md`` (never
editing a past section in place -- each section is an immutable record
of one calibration event, exactly like ``tests/equivalence_calibration_
log.md``'s events) with a fresh measurement pass, or
``TestBenchmarkRegistryTie`` fails. The same edits also invalidate any
`Equiv`/`Diverges` atol/ceiling calibrated against the changed cell's
old fixture -- re-run ``python -m tests.equivalence_measure`` and
update the affected ``ObservableSpec``/``DivergingSpec`` provenance
too, not just the benchmark file.
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

# ---------------------------------------------------------------------------
# Candidate calibration: the two development-machine runner descriptors
# and the measured source commit shared by every ``_measured_equiv`` call
# in this module. One machine is macOS arm64, the other Linux x86_64
# (CPU only, no CUDA involved). Both are ADVISORY pre-freeze proxies:
# the gating continuous-integration runners (ubuntu-latest x Python
# 3.9-3.12) are what ultimately FREEZES these bounds, once a calibration
# run there closes the loop. See ``tests/equivalence_calibration_
# log.md`` (Event 2) for the full per-cell residual tables this
# module's literals were read off of, and
# ``tests/equivalence_benchmark.md`` for the timing side of the same
# sweep. One finding from the Linux x86_64 machine is recorded here
# because it also shapes the CI workflow
# (``.github/workflows/equivalence-calibration.yml`` pins
# ``OPENBLAS_NUM_THREADS``/``OMP_NUM_THREADS``/``MKL_NUM_THREADS``/
# ``NUMEXPR_NUM_THREADS`` to 1 for this reason): with the BLAS backend's
# default thread pool, that machine SEGFAULTS partway through a single
# ``equivalence_measure`` process (repeated small-matrix solver
# constructions across ~27 cells accumulate OpenBLAS worker threads until
# pthread_create fails) -- pinning every BLAS/OMP thread-count env var to
# 1 before invoking the interpreter eliminates it (verified: 3/3 clean
# runs afterward, 0/3 before).
# ---------------------------------------------------------------------------

_CANDIDATE_SOURCE_SHA = "b922a1c13b85b2d319bc65ee8c45183dc6ab2a47"
_CANDIDATE_LOCAL_RUNNER = (
    "macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, "
    "numpy 1.26.4, scipy 1.17.1)"
)
_CANDIDATE_CHITA_RUNNER = (
    "Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, "
    "Python 3.11.15, numpy 1.26.4, scipy 1.13.1)"
)


def _round_up_pow10(x: float) -> float:
    """Round ``x`` (``x`` > 0) UP to the nearest power of ten. Pure
    arithmetic, no I/O -- keeps this module side-effect-free. ``round``
    on the log guards against float noise landing an exact power of ten
    (e.g. ``1e-14``) just above its own exponent boundary and rounding
    up one extra decade.
    """

    if x <= 0.0:
        return 0.0
    exponent = math.ceil(round(math.log10(x), 12))
    return 10.0 ** exponent


def _candidate_atol(local_residual: float, chita_residual: float, ceiling: float) -> float:
    """The tolerance bound rule, in its candidate variant (the two
    development machines in place of the MIN/MAX over the gating
    continuous-integration runners a freeze uses): 10x the larger of the
    two measured residuals, floored at ``max(residual, 1e-15)`` BEFORE
    the 10x multiply, then rounded UP to a power of ten; capped at
    ``ceiling`` (the validator allows ``atol == ceiling``). This
    function only computes the number;
    the STOP-and-report (raw residual > ceiling) and margin-insufficient
    (rounded bound > ceiling) decisions were made by hand while reading
    ``tests/equivalence_calibration_log.md``'s Event 2 table off the
    measured JSON -- neither case arose in this calibration pass (every
    residual measured in the 1e-17..1e-12 range, comfortably under every
    ``POLICY_CEILINGS`` entry).
    """

    raw_max = max(local_residual, chita_residual)
    floored = max(raw_max, 1e-15)
    bound = _round_up_pow10(10.0 * floored)
    return min(bound, ceiling)


def _measured_equiv(
    chi0q_ceiling_key: str,
    chi0q_local: float,
    chi0q_chita: float,
    chiq_comparator: str,
    chiq_ceiling_key: str,
    chiq_local: float,
    chiq_chita: float,
) -> "Equiv":
    """Build the common ``Equiv(chi0q, chiq)`` shape every comparison
    cell in the registry uses: ``chi0q`` always via ``identity``;
    ``chiq`` via the caller's comparator (``general_from_flex_channels``
    or ``reduced_blocks``). This is the candidate-calibration pass:
    ``atol`` is the ``_candidate_atol`` bound (10x the larger of the two
    development-machine residuals, rounded up to a power of ten, capped
    at the mapped policy ceiling); ``provenance`` records BOTH measured
    residuals, both runner descriptors, and the measured source commit
    -- status stays "candidate", and only a closed calibration loop on
    the continuous-integration runners flips the module-level
    ``PROVENANCE`` record to "frozen".
    """

    def _prov(local_residual: float, chita_residual: float, atol: float, ceiling: float) -> str:
        return (
            "measured on two development machines at commit "
            "{}: {} max|diff| {:.3e}; {} max|diff| {:.3e}. Resulting "
            "atol {:.1e} (10x the larger residual, floored at 1e-15, "
            "rounded up to a power of ten; policy ceiling "
            "{:.1e}).".format(
                _CANDIDATE_SOURCE_SHA[:12],
                _CANDIDATE_LOCAL_RUNNER,
                local_residual,
                _CANDIDATE_CHITA_RUNNER,
                chita_residual,
                atol,
                ceiling,
            )
        )

    chi0q_ceiling = POLICY_CEILINGS[chi0q_ceiling_key]
    chiq_ceiling = POLICY_CEILINGS[chiq_ceiling_key]
    chi0q_atol = _candidate_atol(chi0q_local, chi0q_chita, chi0q_ceiling)
    chiq_atol = _candidate_atol(chiq_local, chiq_chita, chiq_ceiling)

    return Equiv(
        observables={
            "chi0q": ObservableSpec(
                comparator="identity",
                atol=chi0q_atol,
                provenance=_prov(chi0q_local, chi0q_chita, chi0q_atol, chi0q_ceiling),
            ),
            "chiq": ObservableSpec(
                comparator=chiq_comparator,
                atol=chiq_atol,
                provenance=_prov(chiq_local, chiq_chita, chiq_atol, chiq_ceiling),
            ),
        }
    )


# The cell inventory. It opens with a 2-cell BOOTSTRAP (rows 1 and 17),
# which exercises the generated-test machinery (build_solver / run_cell /
# the class-body generation loop) end to end on one Equiv comparison cell
# and one BOTH-REJECT cell; the remaining rows follow the same shapes.

_BOOTSTRAP_CELL_1 = Cell(
    cell_id="general.ring.onsite_coulombintra.fixedmu",
    fixture=FixtureSpec(
        input_dir="tests/equivalence_input/orb2",
        interactions={"CoulombIntra": "coulombintra.dat"},
        T=2.0,
        mu=0.0,
        filling=None,
        CellShape=(4, 4, 1),
        SubShape=(1, 1, 1),
        Nmat=32,
        extra_params={},
        calc_type="ring",
        requested_scheme="general",
        enable_spin_orbital=False,
        extern=None,
    ),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    comparison=_measured_equiv(
        "chi0q_fixed", 1.110233872252785e-16, 9.714456569222345e-17,
        "general_from_flex_channels", "chiq_fixed", 2.4980368208390374e-16, 1.942891400167647e-16,
    ),
    required_observables=("chi0q", "chiq"),
    interaction_class="onsite",
    notes=(
        "Bootstrap cell. CoulombIntra U=4.0, orbitals (1,1) and "
        "(2,2) -- tests/equivalence_input/orb2/coulombintra.dat, copied "
        "verbatim from tests/rpa/input_2orb/coulombintra.dat."
    ),
)

_BOOTSTRAP_CELL_17 = Cell(
    cell_id="reduced.ring.onsite_exchange.reject",
    fixture=FixtureSpec(
        input_dir="tests/equivalence_input/orb2",
        interactions={"Exchange": "exchange_onsite.dat"},
        T=2.0,
        mu=0.0,
        filling=None,
        CellShape=(4, 4, 1),
        SubShape=(1, 1, 1),
        Nmat=32,
        extra_params={},
        calc_type="ring",
        requested_scheme="reduced",
        enable_spin_orbital=False,
        extern=None,
    ),
    resolved_scheme="reduced",
    expected_spin_mode="spin-free",
    rpa=SolverProof(
        status=Status.REJECT,
        steps=(
            ExecuteReject(
                site=Site.CONSTRUCTOR,
                exc_type="ValueError",
                fragment="reduced",
            ),
        ),
    ),
    flex=SolverProof(
        status=Status.REJECT,
        steps=(
            ExecuteReject(
                site=Site.CONSTRUCTOR,
                exc_type="ValueError",
                fragment="reduced",
            ),
        ),
    ),
    comparison=None,
    required_observables=(),
    interaction_class="onsite",
    notes=(
        "Bootstrap cell. BOTH-REJECT at the shared _set_scheme "
        "check (src/hwave/solver/rpa.py:1302-1316): calc_scheme='reduced' "
        "solves with the density-diagonal part of the interaction, and "
        "Exchange's particle-hole vertex has no density-diagonal content. "
        "FLEX.__init__ calls super().__init__ (src/hwave/solver/flex.py:"
        "193), so the same constructor-time ValueError fires on both "
        "solvers. tests/equivalence_input/orb2/exchange_onsite.dat, "
        "J=0.20, on-site orbital pair (1,2)."
    ),
)

_SO_SPINFUL_FIXTURE_KWARGS = dict(
    input_dir="tests/equivalence_input/spin",
    interactions={"CoulombIntra": "coulombintra_so.dat"},
    T=2.0,
    mu=0.0,
    filling=None,
    CellShape=(4, 4, 1),
    SubShape=(1, 1, 1),
    Nmat=32,
    extra_params={},
    calc_type="ring",
    enable_spin_orbital=True,
    extern=None,
    geom="geom_so.dat",
    transfer="transfer_spinful.dat",
)

_CELL_35_SO_GENERAL_CONSTRUCTION_REJECT = Cell(
    cell_id="so.general.construction.reject",
    fixture=FixtureSpec(requested_scheme="general", **_SO_SPINFUL_FIXTURE_KWARGS),
    resolved_scheme="general",
    expected_spin_mode="spinful",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(
        status=Status.REJECT,
        steps=(
            ExecuteReject(
                site=Site.CONSTRUCTOR,
                exc_type="ValueError",
                fragment="enable_spin_orbital",
            ),
        ),
    ),
    comparison=None,
    required_observables=(),
    interaction_class="onsite",
    notes=(
        "E4 spinful fixture: tests/equivalence_input/spin/"
        "geom_so.dat (1 physical orbital) + transfer_spinful.dat "
        "(diagonal t=0.7 at R=+-1, on-site spin-mixing 0.35+0.15j / "
        "conjugate) + coulombintra_so.dat (U=4.0), enable_spin_orbital="
        "True. RPA supports SO end to end (verified: constructs and "
        "solves calc_scheme='general', spin_mode resolves to 'spinful'). "
        "FLEX general REJECTS at CONSTRUCTOR -- verified site: "
        "src/hwave/solver/flex.py:220-224 "
        "(FLEX._init_flex_param, scheme=='general' branch), "
        "'calc_scheme=\\'general\\' FLEX (v1) does not support "
        "enable_spin_orbital; deferred to the generalized FLEX solver.' "
        "-- matches the recorded expectation (CONSTRUCTOR, ValueError); "
        "no STOP-and-amend needed."
    ),
)

_CELL_36_SO_REDUCED_CONSTRUCTION_REJECT = Cell(
    cell_id="so.reduced.construction.reject",
    fixture=FixtureSpec(requested_scheme="reduced", **_SO_SPINFUL_FIXTURE_KWARGS),
    resolved_scheme="reduced",
    expected_spin_mode="spinful",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(
        status=Status.REJECT,
        steps=(
            ExecuteReject(
                site=Site.CONSTRUCTOR,
                exc_type="ValueError",
                fragment="enable_spin_orbital",
            ),
        ),
    ),
    comparison=None,
    required_observables=(),
    interaction_class="onsite",
    notes=(
        "Same E4 spinful fixture as cell 35, "
        "calc_scheme='reduced'. RPA supports SO end to end under "
        "reduced too (verified: constructs and solves, spin_mode "
        "resolves to 'spinful'; RPA reduced's spin_mode=='spinful' "
        "branch -- src/hwave/solver/rpa.py:1945-2006 -- runs regardless "
        "of how spin_mode was determined). FLEX reduced REJECTS at "
        "CONSTRUCTOR via the guard added alongside this cell: "
        "src/hwave/solver/flex.py:227-231 "
        "(FLEX._init_flex_param, scheme=='reduced' branch, mirroring "
        "the general-scheme guard at lines 220-224), "
        "'calc_scheme=\\'reduced\\' FLEX does not support "
        "enable_spin_orbital; deferred to the generalized FLEX solver.' "
        "-- before that guard, calc_scheme='reduced' had NO "
        "enable_spin_orbital guard at all (the audit-identified silent "
        "gap it closes)."
    ),
)

# ---------------------------------------------------------------------------
# The remaining cells: rows 2-16 and 18-34 (rows 1, 17, 35 and 36 are
# the bootstrap above; row 38 is the conditioning cell below). Every
# ``_measured_equiv`` call below (the helper is defined near
# ``POLICY_CEILINGS`` above, so the bootstrap cells can use it too) sets
# a CANDIDATE atol from the development-machine calibration sweep -- see
# ``tests/equivalence_calibration_log.md`` (Event 2) and
# ``tests/equivalence_benchmark.md`` for the full residual tables.
# ---------------------------------------------------------------------------


# -- G1 shared fixture shape: E2, calc_scheme='general', fixed mu=0.0 -------# -- G1 shared fixture shape: E2, calc_scheme='general', fixed mu=0.0 -------

_E2_GENERAL_FIXEDMU_KWARGS = dict(
    input_dir="tests/equivalence_input/orb2",
    T=2.0,
    mu=0.0,
    filling=None,
    CellShape=(4, 4, 1),
    SubShape=(1, 1, 1),
    Nmat=32,
    extra_params={},
    calc_type="ring",
    requested_scheme="general",
    enable_spin_orbital=False,
    extern=None,
)

# -- G1 shared fixture shape: E2, calc_scheme='general', mu-coupled --------

_E2_GENERAL_MU_KWARGS = dict(_E2_GENERAL_FIXEDMU_KWARGS)
_E2_GENERAL_MU_KWARGS.update(mu=None, filling=0.5)

_CELL_2_ONSITE_COULOMBINTER_FIXEDMU = Cell(
    cell_id="general.ring.onsite_coulombinter.fixedmu",
    fixture=FixtureSpec(
        interactions={"CoulombInter": "onsite_inter.dat"},
        **_E2_GENERAL_FIXEDMU_KWARGS,
    ),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    comparison=_measured_equiv(
        "chi0q_fixed", 1.110233872252785e-16, 9.714456569222345e-17,
        "general_from_flex_channels", "chiq_fixed", 1.2490109573171623e-16, 1.1102239035450868e-16,
    ),
    required_observables=("chi0q", "chiq"),
    interaction_class="onsite",
    notes=(
        "E2 onsite_inter.dat -- on-site inter-orbital CoulombInter "
        "V=0.7, orbital pair (1,2)."
    ),
)

_CELL_3_ONSITE_HUND_FIXEDMU = Cell(
    cell_id="general.ring.onsite_hund.fixedmu",
    fixture=FixtureSpec(
        interactions={"Hund": "hund_onsite.dat"},
        **_E2_GENERAL_FIXEDMU_KWARGS,
    ),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    comparison=_measured_equiv(
        "chi0q_fixed", 1.110233872252785e-16, 9.714456569222345e-17,
        "general_from_flex_channels", "chiq_fixed", 1.1102343486699909e-16, 1.1102304692257097e-16,
    ),
    required_observables=("chi0q", "chiq"),
    interaction_class="onsite",
    notes="E2 hund_onsite.dat -- on-site Hund J=1.0, orbital pair (1,2).",
)

_CELL_4_ONSITE_ISING_FIXEDMU = Cell(
    cell_id="general.ring.onsite_ising.fixedmu",
    fixture=FixtureSpec(
        interactions={"Ising": "ising_onsite.dat"},
        **_E2_GENERAL_FIXEDMU_KWARGS,
    ),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    comparison=_measured_equiv(
        "chi0q_fixed", 1.110233872252785e-16, 9.714456569222345e-17,
        "general_from_flex_channels", "chiq_fixed", 1.2490106872585652e-16, 1.1102283837830983e-16,
    ),
    required_observables=("chi0q", "chiq"),
    interaction_class="onsite",
    notes=(
        "E2 ising_onsite.dat -- on-site Ising J=0.15, orbital pair "
        "(1,2)."
    ),
)

_CELL_5_ONSITE_EXCHANGE_FIXEDMU = Cell(
    cell_id="general.ring.onsite_exchange.fixedmu",
    fixture=FixtureSpec(
        interactions={"Exchange": "exchange_onsite.dat"},
        **_E2_GENERAL_FIXEDMU_KWARGS,
    ),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    comparison=_measured_equiv(
        "chi0q_fixed", 1.110233872252785e-16, 9.714456569222345e-17,
        "general_from_flex_channels", "chiq_fixed", 1.110233872252785e-16, 9.714456569222345e-17,
    ),
    required_observables=("chi0q", "chiq"),
    interaction_class="onsite",
    notes=(
        "E2 exchange_onsite.dat -- on-site Exchange J=0.20, orbital "
        "pair (1,2)."
    ),
)

_CELL_6_ONSITE_PAIRHOP_FIXEDMU = Cell(
    cell_id="general.ring.onsite_pairhop.fixedmu",
    fixture=FixtureSpec(
        interactions={"PairHop": "pairhop_onsite.dat"},
        **_E2_GENERAL_FIXEDMU_KWARGS,
    ),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    comparison=_measured_equiv(
        "chi0q_fixed", 1.110233872252785e-16, 9.714456569222345e-17,
        "general_from_flex_channels", "chiq_fixed", 1.110233872252785e-16, 9.714538392232159e-17,
    ),
    required_observables=("chi0q", "chiq"),
    interaction_class="onsite",
    notes=(
        "E2 pairhop_onsite.dat -- on-site PairHop P=0.10, orbital "
        "pair (1,2)."
    ),
)

_CELL_7_ONSITE_PAIRLIFT_FIXEDMU = Cell(
    cell_id="general.ring.onsite_pairlift.fixedmu",
    fixture=FixtureSpec(
        interactions={"PairLift": "pairlift_onsite.dat"},
        **_E2_GENERAL_FIXEDMU_KWARGS,
    ),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    comparison=_measured_equiv(
        "chi0q_fixed", 1.110233872252785e-16, 9.714456569222345e-17,
        "general_from_flex_channels", "chiq_fixed", 1.110233872252785e-16, 9.714456569222345e-17,
    ),
    required_observables=("chi0q", "chiq"),
    interaction_class="onsite",
    notes=(
        "E2 pairlift_onsite.dat -- on-site PairLift L=0.15, orbital "
        "pair (1,2). FLEX warns (not rejects) that PairLift's "
        "particle-hole vertex is exactly zero (S=C=0, "
        "flex.py:1941-1945); it is physically inert but not an error."
    ),
)

_CELL_8_ONSITE_U_V_HUND_MU = Cell(
    cell_id="general.ring.onsite_u_v_hund.mu",
    fixture=FixtureSpec(
        interactions={
            "CoulombIntra": "coulombintra.dat",
            "CoulombInter": "onsite_inter.dat",
            "Hund": "hund_onsite.dat",
        },
        **_E2_GENERAL_MU_KWARGS,
    ),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    comparison=_measured_equiv(
        "chi0q_mu", 1.110233872252785e-16, 9.714456569222345e-17,
        "general_from_flex_channels", "chiq_mu", 3.0531310938257943e-16, 2.77555831224212e-16,
    ),
    required_observables=("chi0q", "chiq"),
    interaction_class="onsite",
    notes=(
        "E2 CoulombIntra U=4.0 + onsite_inter.dat V=0.7 + "
        "hund_onsite.dat J=1.0 -- the oneshot suite's own U+V+J "
        "composite (tests/test_rpa_flex_oneshot_equivalence.py's "
        "TestGeneralSchemeOneShot.CASES); filling=0.5."
    ),
)

_CELL_9_ONSITE_FULL_KANAMORI_MU = Cell(
    cell_id="general.ring.onsite_full_kanamori.mu",
    fixture=FixtureSpec(
        interactions={
            "CoulombIntra": "coulombintra.dat",
            "CoulombInter": "kanamori_interorb.dat",
            "Hund": "kanamori_hund.dat",
            "PairHop": "pairhop_onsite.dat",
            "Ising": "ising_onsite.dat",
            "Exchange": "exchange_onsite.dat",
        },
        **_E2_GENERAL_MU_KWARGS,
    ),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    comparison=_measured_equiv(
        "chi0q_mu", 1.110233872252785e-16, 9.714456569222345e-17,
        "general_from_flex_channels", "chiq_mu", 2.4980382961466737e-16, 1.942901993007593e-16,
    ),
    required_observables=("chi0q", "chiq"),
    interaction_class="onsite",
    notes=(
        "FROZEN full-Kanamori couplings (tests/test_flex_general.py): "
        "CoulombIntra U=4.0 (coulombintra.dat), CoulombInter V=1.0 "
        "(kanamori_interorb.dat), Hund J=0.30 (kanamori_hund.dat), "
        "PairHop P=0.10 (pairhop_onsite.dat), Ising J=0.15 "
        "(ising_onsite.dat), Exchange J=0.20 (exchange_onsite.dat); "
        "filling=0.5."
    ),
)

# -- G2 shared fixture shape: E2, calc_scheme='reduced', mu-coupled --------

_E2_REDUCED_MU_KWARGS = dict(
    input_dir="tests/equivalence_input/orb2",
    T=2.0,
    mu=None,
    filling=0.5,
    CellShape=(4, 4, 1),
    SubShape=(1, 1, 1),
    Nmat=32,
    extra_params={},
    calc_type="ring",
    requested_scheme="reduced",
    enable_spin_orbital=False,
    extern=None,
)

_CELL_10_REDUCED_COULOMBINTRA_SPINFREE_MU = Cell(
    cell_id="reduced.ring.onsite_coulombintra.spinfree.mu",
    fixture=FixtureSpec(
        interactions={"CoulombIntra": "coulombintra.dat"},
        **_E2_REDUCED_MU_KWARGS,
    ),
    resolved_scheme="reduced",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    comparison=_measured_equiv(
        "chi0q_mu", 1.110233872252785e-16, 9.714456569222345e-17,
        "reduced_blocks", "chiq_mu", 3.885849304444534e-16, 3.0531155156219455e-16,
    ),
    required_observables=("chi0q", "chiq"),
    interaction_class="onsite",
    notes=(
        "E2 coulombintra.dat, U=4.0, reduced/spin-free; filling=0.5. "
        "chi0q FINDING: FLEX writes the FULL spin-block-"
        "inflated bare bubble to green_info['chi0q'] under "
        "calc_scheme='reduced' (shape (nmat,nvol,2*norb,2*norb), "
        "matching chiq_s/chiq_c's convention), while RPA's reduced "
        "chi0q stays in its already-reduced density-diagonal shape "
        "(spin-free: (nmat,nvol,norb,norb)). Never caught before -- "
        "the existing oneshot suite (TestReducedOneShot) only "
        "compares chiq under calc_scheme='reduced'. extract_bundle "
        "(test module) now slices FLEX's reduced chi0q to its "
        "up-up norb x norb diagonal block before the 'identity' "
        "comparator runs, mirroring the SAME block extraction "
        "'reduced_blocks' already applies to chiq_s/chiq_c; measured "
        "equal to round-off on every reduced fixture checked here "
        "(E1 norb=1, E2 norb=2 spin-free, E4 spin-diag)."
    ),
)

_CELL_11_REDUCED_COULOMBINTER_SPINFREE_MU = Cell(
    cell_id="reduced.ring.onsite_coulombinter.spinfree.mu",
    fixture=FixtureSpec(
        interactions={"CoulombInter": "onsite_inter.dat"},
        **_E2_REDUCED_MU_KWARGS,
    ),
    resolved_scheme="reduced",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    comparison=_measured_equiv(
        "chi0q_mu", 1.110233872252785e-16, 9.714456569222345e-17,
        "reduced_blocks", "chiq_mu", 1.2490113201110995e-16, 1.249002003256928e-16,
    ),
    required_observables=("chi0q", "chiq"),
    interaction_class="onsite",
    notes=(
        "E2 onsite_inter.dat, V=0.7, reduced/spin-free; filling=0.5. "
        "See cell 10's notes for the chi0q block-extraction finding "
        "(applies uniformly to every reduced cell)."
    ),
)

_CELL_12_REDUCED_HUND_SPINFREE_MU = Cell(
    cell_id="reduced.ring.onsite_hund.spinfree.mu",
    fixture=FixtureSpec(
        interactions={"Hund": "hund_onsite.dat"},
        **_E2_REDUCED_MU_KWARGS,
    ),
    resolved_scheme="reduced",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    comparison=_measured_equiv(
        "chi0q_mu", 1.110233872252785e-16, 9.714456569222345e-17,
        "reduced_blocks", "chiq_mu", 1.1102343486699909e-16, 1.1102304692257097e-16,
    ),
    required_observables=("chi0q", "chiq"),
    interaction_class="onsite",
    notes=(
        "E2 hund_onsite.dat, J=1.0, reduced/spin-free; filling=0.5. "
        "See cell 10's notes for the chi0q block-extraction finding."
    ),
)

_CELL_13_REDUCED_ISING_SPINFREE_MU = Cell(
    cell_id="reduced.ring.onsite_ising.spinfree.mu",
    fixture=FixtureSpec(
        interactions={"Ising": "ising_onsite.dat"},
        **_E2_REDUCED_MU_KWARGS,
    ),
    resolved_scheme="reduced",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    comparison=_measured_equiv(
        "chi0q_mu", 1.110233872252785e-16, 9.714456569222345e-17,
        "reduced_blocks", "chiq_mu", 1.2490105326199038e-16, 1.2490057062540458e-16,
    ),
    required_observables=("chi0q", "chiq"),
    interaction_class="onsite",
    notes=(
        "E2 ising_onsite.dat, J=0.15, reduced/spin-free; filling=0.5. "
        "See cell 10's notes for the chi0q block-extraction finding."
    ),
)

_CELL_14_REDUCED_PAIRLIFT_SPINFREE_MU = Cell(
    cell_id="reduced.ring.onsite_pairlift.spinfree.mu",
    fixture=FixtureSpec(
        interactions={"PairLift": "pairlift_onsite.dat"},
        **_E2_REDUCED_MU_KWARGS,
    ),
    resolved_scheme="reduced",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    comparison=_measured_equiv(
        "chi0q_mu", 1.110233872252785e-16, 9.714456569222345e-17,
        "reduced_blocks", "chiq_mu", 1.110233872252785e-16, 9.714456569222345e-17,
    ),
    required_observables=("chi0q", "chiq"),
    interaction_class="onsite",
    notes=(
        "E2 pairlift_onsite.dat, L=0.15, reduced/spin-free; "
        "filling=0.5. PairLift's acceptance under reduced is a "
        "RECORDED FACT of the current resolver (its particle-hole "
        "vertex is exactly zero, so 'reduced' has nothing to drop -- "
        "unlike Exchange/PairHop, whose nonzero cross-slot vertex the "
        "density-diagonal projection WOULD drop and cell 17/18 "
        "reject for that reason). FLEX warns (not rejects) that "
        "PairLift is physically inert (rpa.py:1317-1322 / "
        "flex.py:1941-1945). See cell 10's notes for the chi0q "
        "block-extraction finding."
    ),
)

# -- G2 shared fixture shape: E4, calc_scheme='reduced', spin-diag ---------

_E4_REDUCED_SPINDIAG_KWARGS = dict(
    input_dir="tests/equivalence_input/spin",
    T=2.0,
    mu=None,
    filling=0.5,
    CellShape=(4, 4, 1),
    SubShape=(1, 1, 1),
    Nmat=32,
    extra_params={"coeff_extern": 1.0},
    calc_type="ring",
    requested_scheme="reduced",
    enable_spin_orbital=False,
    extern="extern_zeeman.dat",
)

_CELL_15_REDUCED_COULOMBINTRA_SPINDIAG_MU = Cell(
    cell_id="reduced.ring.onsite_coulombintra.spindiag.mu",
    fixture=FixtureSpec(
        interactions={"CoulombIntra": "coulombintra.dat"},
        **_E4_REDUCED_SPINDIAG_KWARGS,
    ),
    resolved_scheme="reduced",
    expected_spin_mode="spin-diag",
    rpa=SolverProof(
        status=Status.SUPPORTED,
        steps=(ExecuteRun(),),
        links=(
            SupplementaryLink(
                test_id=(
                    "tests.test_rpa_flex_equivalence_table::"
                    "TestReducedSpinDiagAcceptance::"
                    "test_reduced_rpa_also_solves_with_extern"
                ),
                claim=(
                    "RPA accepts the Extern/Zeeman spin-diagonal route "
                    "(the tests/equivalence_input/spin fixture) under "
                    "calc_scheme='reduced'."
                ),
            ),
        ),
    ),
    flex=SolverProof(
        status=Status.SUPPORTED,
        steps=(ExecuteRun(),),
        links=(
            SupplementaryLink(
                test_id=(
                    "tests.test_rpa_flex_equivalence_table::"
                    "TestReducedSpinDiagAcceptance::"
                    "test_reduced_flex_solves_with_extern_coulombintra"
                ),
                claim=(
                    "FLEX accepts the Extern/Zeeman spin-diagonal route "
                    "(the tests/equivalence_input/spin fixture) under "
                    "calc_scheme='reduced' -- the recorded expectation "
                    "for this row is that both solvers accept it."
                ),
            ),
        ),
    ),
    comparison=_measured_equiv(
        "chi0q_mu", 1.2490992897474145e-16, 1.249078513798569e-16,
        "reduced_blocks", "chiq_mu", 4.163877857486877e-16, 4.718522149908535e-16,
    ),
    required_observables=("chi0q", "chiq"),
    interaction_class="onsite",
    notes=(
        "E4 spin-diag route: coulombintra.dat (U=4.0, copied from "
        "E2) + extern_zeeman.dat (hz=0.1, coeff_extern=1.0); "
        "filling=0.5. chi0q comparison: RPA's spin-diag chi0q is "
        "5-dim ((2,nmat,nvol,norb,norb), one block per spin channel); "
        "extract_bundle stacks FLEX's inflated chi0q's up-up and "
        "down-down diagonal blocks into the same shape -- measured "
        "equal to round-off, and the off-diagonal (spin-mixing) "
        "blocks of FLEX's chi0q are exactly zero, as expected at the "
        "bare-bubble level."
    ),
)

_CELL_16_REDUCED_COULOMBINTER_SPINDIAG_MU = Cell(
    cell_id="reduced.ring.onsite_coulombinter.spindiag.mu",
    fixture=FixtureSpec(
        interactions={"CoulombInter": "onsite_inter.dat"},
        **_E4_REDUCED_SPINDIAG_KWARGS,
    ),
    resolved_scheme="reduced",
    expected_spin_mode="spin-diag",
    rpa=SolverProof(
        status=Status.SUPPORTED,
        steps=(ExecuteRun(),),
        links=(
            SupplementaryLink(
                test_id=(
                    "tests.test_rpa_flex_equivalence_table::"
                    "TestReducedSpinDiagAcceptance::"
                    "test_reduced_rpa_also_solves_with_extern"
                ),
                claim=(
                    "RPA accepts the Extern/Zeeman spin-diagonal route "
                    "(the tests/equivalence_input/spin fixture) under "
                    "calc_scheme='reduced'."
                ),
            ),
        ),
    ),
    flex=SolverProof(
        status=Status.SUPPORTED,
        steps=(ExecuteRun(),),
        links=(
            SupplementaryLink(
                test_id=(
                    "tests.test_rpa_flex_equivalence_table::"
                    "TestReducedSpinDiagAcceptance::"
                    "test_reduced_flex_solves_with_extern_coulombinter"
                ),
                claim=(
                    "FLEX accepts the Extern/Zeeman spin-diagonal route "
                    "(the tests/equivalence_input/spin fixture) under "
                    "calc_scheme='reduced' with CoulombInter."
                ),
            ),
        ),
    ),
    comparison=_measured_equiv(
        "chi0q_mu", 1.2490992897474145e-16, 1.249078513798569e-16,
        "reduced_blocks", "chiq_mu", 1.5266161502687854e-16, 1.3878173264691723e-16,
    ),
    required_observables=("chi0q", "chiq"),
    interaction_class="onsite",
    notes=(
        "E4 spin-diag route: onsite_inter.dat (V=0.7, copied from "
        "E2) + extern_zeeman.dat; filling=0.5. See cell 15's notes "
        "for the chi0q stacking construction."
    ),
)

_CELL_18_REDUCED_PAIRHOP_REJECT = Cell(
    cell_id="reduced.ring.onsite_pairhop.reject",
    fixture=FixtureSpec(
        interactions={"PairHop": "pairhop_onsite.dat"},
        **dict(_E2_GENERAL_FIXEDMU_KWARGS, requested_scheme="reduced"),
    ),
    resolved_scheme="reduced",
    expected_spin_mode="spin-free",
    rpa=SolverProof(
        status=Status.REJECT,
        steps=(ExecuteReject(site=Site.CONSTRUCTOR, exc_type="ValueError", fragment="reduced"),),
    ),
    flex=SolverProof(
        status=Status.REJECT,
        steps=(ExecuteReject(site=Site.CONSTRUCTOR, exc_type="ValueError", fragment="reduced"),),
    ),
    comparison=None,
    required_observables=(),
    interaction_class="onsite",
    notes=(
        "E2 pairhop_onsite.dat, P=0.10. BOTH-REJECT at the shared "
        "_set_scheme check (rpa.py:1302-1316), same shape as cell "
        "17: PairHop's particle-hole vertex has no density-diagonal "
        "content."
    ),
)

# ---------------------------------------------------------------------------
# Group G3 -- off-site (E1/E2), cells 19-28.
# ---------------------------------------------------------------------------

_CELL_19_OFFSITE_COULOMBINTER_SAMEORB_MU = Cell(
    cell_id="general.ring.offsite_coulombinter_sameorb.mu",
    fixture=FixtureSpec(
        input_dir="tests/equivalence_input/orb1",
        interactions={"CoulombInter": "coulombinter.dat"},
        T=2.0, mu=None, filling=0.75,
        CellShape=(4, 4, 1), SubShape=(1, 1, 1), Nmat=32,
        extra_params={}, calc_type="ring",
        requested_scheme="general", enable_spin_orbital=False, extern=None,
    ),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(
        status=Status.SUPPORTED,
        steps=(ExecuteRun(),),
        links=(
            SupplementaryLink(
                test_id=(
                    "tests.test_flex_offsite_general::"
                    "TestOffsiteGeneralFLEX::"
                    "test_one_orbital_offsite_v_matches_the_rpa_ring_exactly"
                ),
                claim=(
                    "The same fixture (tests/rpa/input's "
                    "coulombinter.dat, filling=0.75) element-complete "
                    "equal between RPA and FLEX for off-site "
                    "same-orbital CoulombInter."
                ),
            ),
        ),
    ),
    comparison=_measured_equiv(
        "chi0q_mu", 4.884981379996393e-15, 4.898859387959854e-15,
        "general_from_flex_channels", "chiq_mu", 2.0206059344954128e-14, 2.0261571408243727e-14,
    ),
    required_observables=("chi0q", "chiq"),
    interaction_class="offsite",
    notes=(
        "E1 coulombinter.dat -- off-site CoulombInter, same-orbital "
        "(a==b), R=+-x/+-y bonds, V=1.0; filling=0.75. The audit's "
        "exact-match case: FLEX's general path is measured equal to "
        "the RPA ring for this entry class (flex.py:1947-1970's own "
        "documented rationale)."
    ),
)

_CELL_20_REDUCED_OFFSITE_COULOMBINTER_MU = Cell(
    cell_id="reduced.ring.offsite_coulombinter.mu",
    fixture=FixtureSpec(
        input_dir="tests/equivalence_input/orb1",
        interactions={"CoulombInter": "coulombinter.dat"},
        T=2.0, mu=None, filling=0.75,
        CellShape=(4, 4, 1), SubShape=(1, 1, 1), Nmat=32,
        extra_params={}, calc_type="ring",
        requested_scheme="reduced", enable_spin_orbital=False, extern=None,
    ),
    resolved_scheme="reduced",
    expected_spin_mode="spin-free",
    rpa=SolverProof(
        status=Status.SUPPORTED,
        steps=(ExecuteRun(),),
        links=(
            SupplementaryLink(
                test_id=(
                    "tests.test_rpa_flex_oneshot_equivalence::"
                    "TestReducedOneShot::test_matrix_cells"
                ),
                claim=(
                    "An existing one-shot comparison over this same "
                    "fixture and filling pins the chiq uu-ud/uu+ud "
                    "block equivalence under calc_scheme='reduced'."
                ),
            ),
        ),
    ),
    flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    comparison=_measured_equiv(
        "chi0q_mu", 4.884981379996393e-15, 4.898859387959854e-15,
        "reduced_blocks", "chiq_mu", 3.552713730991191e-14, 3.5610405641548485e-14,
    ),
    required_observables=("chi0q", "chiq"),
    interaction_class="offsite",
    notes=(
        "E1 coulombinter.dat, reduced; filling=0.75 -- the oneshot "
        "suite's TestReducedOneShot case (a) "
        "('tests/rpa/input', CoulombInter only, filling=0.75, norb=1), "
        "which only ever compared chiq. See cell 10's notes for the "
        "chi0q block-extraction finding this cell adds on top."
    ),
)

_CELL_21_OFFSITE_COULOMBINTER_INTERORB_FLEXREJECT = Cell(
    cell_id="general.ring.offsite_coulombinter_interorb.flexreject",
    fixture=FixtureSpec(
        interactions={"CoulombInter": "offsite_coulombinter_interorb.dat"},
        **_E2_GENERAL_MU_KWARGS,
    ),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(
        status=Status.REJECT,
        steps=(
            ExecuteReject(
                site=Site.SOLVE, exc_type="ValueError",
                fragment="interaction 'CoulombInter'",
            ),
        ),
        links=(
            SupplementaryLink(
                test_id=(
                    "tests.test_flex_offsite_general::"
                    "TestOffsiteGeneralFLEX::test_rejected_offsite_classes"
                ),
                claim=(
                    "FLEX rejects off-site CoulombInter with a != b "
                    "(inter-orbital) under calc_scheme='general'."
                ),
            ),
        ),
    ),
    comparison=None,
    required_observables=(),
    interaction_class="offsite",
    notes=(
        "E2 offsite_coulombinter_interorb.dat -- off-site CoulombInter, "
        "inter-orbital (1,2), R=(1,0,0), coupling 0.2; filling=0.5. "
        "RPA=SUPPORTED (off-site inter-orbital CoulombInter has no "
        "known equivalence gap on the RPA side). FLEX SOLVE-time "
        "reject: flex.py:2020-2043 (_flex_compute_veff_general -> "
        "_inflate_chi0q_and_ham_general), the off-site guard's "
        "a==b-only CoulombInter exemption."
    ),
)

_CELL_22_OFFSITE_HUND_FLEXREJECT = Cell(
    cell_id="general.ring.offsite_hund.flexreject",
    fixture=FixtureSpec(
        interactions={"Hund": "offsite_hund.dat"},
        **_E2_GENERAL_MU_KWARGS,
    ),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(
        status=Status.REJECT,
        steps=(
            ExecuteReject(
                site=Site.SOLVE, exc_type="ValueError",
                fragment="interaction 'Hund'",
            ),
        ),
        links=(
            SupplementaryLink(
                test_id=(
                    "tests.test_flex_offsite_general::"
                    "TestOffsiteGeneralFLEX::test_rejected_offsite_classes"
                ),
                claim="FLEX rejects off-site Hund under calc_scheme='general'.",
            ),
        ),
    ),
    comparison=None,
    required_observables=(),
    interaction_class="offsite",
    notes=(
        "E2 offsite_hund.dat, R=(1,0,0), coupling 0.2; filling=0.5. "
        "FLEX SOLVE-time reject: flex.py:2020-2043 -- off-site Hund "
        "is not representable by a q-only vertex (non-local "
        "particle-hole pair)."
    ),
)

_CELL_23_OFFSITE_ISING_FLEXREJECT = Cell(
    cell_id="general.ring.offsite_ising.flexreject",
    fixture=FixtureSpec(
        interactions={"Ising": "offsite_ising.dat"},
        **_E2_GENERAL_MU_KWARGS,
    ),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(
        status=Status.REJECT,
        steps=(
            ExecuteReject(
                site=Site.SOLVE, exc_type="ValueError",
                fragment="interaction 'Ising'",
            ),
        ),
        links=(
            SupplementaryLink(
                test_id=(
                    "tests.test_flex_offsite_general::"
                    "TestOffsiteGeneralFLEX::test_rejected_offsite_classes"
                ),
                claim="FLEX rejects off-site Ising under calc_scheme='general'.",
            ),
        ),
    ),
    comparison=None,
    required_observables=(),
    interaction_class="offsite",
    notes=(
        "E2 offsite_ising.dat, R=(1,0,0), coupling 0.2; filling=0.5. "
        "FLEX SOLVE-time reject: flex.py:2020-2043."
    ),
)

_CELL_24_OFFSITE_EXCHANGE_FLEXREJECT = Cell(
    cell_id="general.ring.offsite_exchange.flexreject",
    fixture=FixtureSpec(
        interactions={"Exchange": "offsite_exchange.dat"},
        **_E2_GENERAL_MU_KWARGS,
    ),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(
        status=Status.REJECT,
        steps=(
            ExecuteReject(
                site=Site.SOLVE, exc_type="ValueError",
                fragment="interaction 'Exchange'",
            ),
        ),
        links=(
            SupplementaryLink(
                test_id=(
                    "tests.test_flex_offsite_general::"
                    "TestOffsiteGeneralFLEX::test_rejected_offsite_classes"
                ),
                claim="FLEX rejects off-site Exchange under calc_scheme='general'.",
            ),
        ),
    ),
    comparison=None,
    required_observables=(),
    interaction_class="offsite",
    notes=(
        "E2 offsite_exchange.dat, R=(1,0,0), coupling 0.2; "
        "filling=0.5. FLEX SOLVE-time reject: flex.py:2020-2043 -- "
        "non-local particle-hole pair off-site."
    ),
)

_CELL_25_OFFSITE_PAIRHOP_FLEXREJECT = Cell(
    cell_id="general.ring.offsite_pairhop.flexreject",
    fixture=FixtureSpec(
        interactions={"PairHop": "offsite_pairhop.dat"},
        **_E2_GENERAL_MU_KWARGS,
    ),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(
        status=Status.REJECT,
        steps=(
            ExecuteReject(
                site=Site.SOLVE, exc_type="ValueError",
                fragment="interaction 'PairHop'",
            ),
        ),
        links=(
            SupplementaryLink(
                test_id=(
                    "tests.test_flex_offsite_general::"
                    "TestOffsiteGeneralFLEX::test_rejected_offsite_classes"
                ),
                claim=(
                    "tests/test_flex_offsite_general.py's general "
                    "rejection-classes coverage (this exact PairHop "
                    "case is not one of its enumerated subTest rows -- "
                    "the module's docstring names off-site PairHop as "
                    "rejected for the same non-local-pair reason as "
                    "Exchange; linked for the surrounding module "
                    "context, not a per-type pin)."
                ),
            ),
        ),
    ),
    comparison=None,
    required_observables=(),
    interaction_class="offsite",
    notes=(
        "E2 offsite_pairhop.dat, R=(1,0,0), coupling 0.2; "
        "filling=0.5. RPA reads off-site PairHop but silently drops "
        "it (issue #157, on-site-only representation) -- SUPPORTED "
        "still holds since RPA does not raise. FLEX SOLVE-time "
        "reject: flex.py:2020-2043 -- non-local particle-hole pair "
        "off-site, same as Exchange."
    ),
)

_CELL_26_OFFSITE_COULOMBINTRA_LITERALKEY_REJECT = Cell(
    cell_id="general.ring.offsite_coulombintra_literalkey.reject",
    fixture=FixtureSpec(
        interactions={"CoulombIntra": "offsite_coulombintra.dat"},
        **_E2_GENERAL_MU_KWARGS,
    ),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(
        status=Status.REJECT,
        steps=(
            ExecuteReject(
                site=Site.CONSTRUCTOR, exc_type="ValueError",
                fragment="the documented operator is on-site and same-orbital",
            ),
        ),
    ),
    flex=SolverProof(
        status=Status.REJECT,
        steps=(
            ExecuteReject(
                site=Site.CONSTRUCTOR, exc_type="ValueError",
                fragment="the documented operator is on-site and same-orbital",
            ),
        ),
    ),
    comparison=None,
    required_observables=(),
    interaction_class="offsite",
    notes=(
        "AMENDED (2026-08-18, Task-3 finding, reviewer-confirmed): "
        "originally intended as a FLEX-solve-time rejection; the "
        "literal CoulombIntra key rejects ANY off-site/inter-orbital "
        "entry at READ time instead (declarations.py:195-211, "
        "validate_hermitian_closure, called from "
        "read_input_k.py:125-128 -- runs inside the fixture read, "
        "i.e. before either solver constructs). VERIFIED: "
        "reading tests/equivalence_input/orb2/offsite_coulombintra.dat "
        "under the literal 'CoulombIntra' key raises ValueError: "
        "\"CoulombIntra ... declares R=(1, 0, 0) orbitals (1, 1) = "
        "(0.2+0j): the documented operator is on-site and "
        "same-orbital (R = 0, a == b); off-site or inter-orbital "
        "entries belong in CoulombInter\" -- both solvers share this "
        "read, so BOTH-REJECT at the CONSTRUCTOR site (the read "
        "happens inside run_cell's construction phase, before "
        "RPA()/FLEX() are ever called). The aggregate 'Coulomb' route "
        "instead classifies this file's same-orbital off-site content "
        "as CoulombInter (cell 19's physics, FLEX-accepted) -- so the "
        "originally intended row has no reachable public "
        "configuration; this row now documents the input-layer "
        "rejection shared by both solvers."
    ),
)

_CELL_27_OFFSITE_COULOMBINTER_SAMEORB_SUBSHAPE = Cell(
    cell_id="general.ring.offsite_coulombinter_sameorb.subshape",
    fixture=FixtureSpec(
        input_dir="tests/equivalence_input/orb1",
        interactions={"CoulombInter": "coulombinter.dat"},
        T=2.0, mu=0.0, filling=None,
        CellShape=(4, 4, 1), SubShape=(2, 1, 1), Nmat=32,
        extra_params={}, calc_type="ring",
        requested_scheme="general", enable_spin_orbital=False, extern=None,
    ),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(
        status=Status.REJECT,
        steps=(
            ExecuteReject(
                site=Site.SOLVE, exc_type="ValueError",
                fragment="with sublattice folding",
            ),
        ),
    ),
    comparison=None,
    required_observables=(),
    interaction_class="offsite",
    notes=(
        "Cell 19's fixture (E1 coulombinter.dat, same-orbital off-site "
        "CoulombInter) with SubShape=(2,1,1); mu=0.0. STOP-and-amend "
        "audit VERIFIED: per the audit, FLEX REJECTS folded "
        "off-site entries even in the otherwise-accepted a==b "
        "CoulombInter class -- confirmed empirically: RPA solves "
        "(chi0q shape (32,8,2,2,2,2)); FLEX raises ValueError at "
        "SOLVE time (flex.py:2020-2043's has_fold branch), message "
        "'...interaction 'CoulombInter' has an off-site entry "
        "irvec=(1, 0, 0), orbvec=(0, 1), with sublattice folding...'. "
        "Matches the recorded expectation exactly -- NO STOP-and-amend "
        "triggered; FLEX-REJECT.RPA-SUPPORTED as predicted."
    ),
)

_CELL_28_ONSITE_COULOMBINTER_SUBSHAPE_MU = Cell(
    cell_id="general.ring.onsite_coulombinter.subshape.mu",
    fixture=FixtureSpec(
        interactions={"CoulombInter": "onsite_inter.dat"},
        **dict(_E2_GENERAL_MU_KWARGS, SubShape=(2, 1, 1)),
    ),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    comparison=_measured_equiv(
        "chi0q_mu", 1.6653615627141079e-16, 9.714691283011063e-17,
        "general_from_flex_channels", "chiq_mu", 1.9429144066570713e-16, 1.1102444068327482e-16,
    ),
    required_observables=("chi0q", "chiq"),
    interaction_class="onsite",
    notes=(
        "E2 onsite_inter.dat, SubShape=(2,1,1), filling=0.5 -- the "
        "audit's untested on-site folding path: sublattice folding "
        "combined with an on-site (not off-site) interaction. Equiv "
        "at round-off, unlike cell 27's off-site+folding combination."
    ),
)

# ---------------------------------------------------------------------------
# Group G4 -- options, cells 29-34.
# ---------------------------------------------------------------------------

_CELL_29_ONSITE_COULOMBINTER_COEFFTAIL_MU = Cell(
    cell_id="general.ring.onsite_coulombinter.coefftail.mu",
    fixture=FixtureSpec(
        interactions={"CoulombInter": "onsite_inter.dat"},
        **dict(_E2_GENERAL_MU_KWARGS, extra_params={"coeff_tail": 1.0}),
    ),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    comparison=_measured_equiv(
        "chi0q_mu", 1.110223918911846e-16, 1.1102377660297755e-16,
        "general_from_flex_channels", "chiq_mu", 1.3877796215672092e-16, 1.387790961379121e-16,
    ),
    required_observables=("chi0q", "chiq"),
    interaction_class="onsite",
    notes=(
        "E2 onsite_inter.dat, coeff_tail=1.0, filling=0.5 -- the "
        "audit's untested combination (the high-frequency tail "
        "correction, active) with a mu-coupled on-site cell. Equiv "
        "at round-off."
    ),
)

_CELL_30_AUTO_DENSITY_RESOLUTION = Cell(
    cell_id="auto.density.resolution",
    fixture=FixtureSpec(
        interactions={"CoulombInter": "onsite_inter.dat"},
        **dict(_E2_GENERAL_MU_KWARGS, requested_scheme="auto"),
    ),
    resolved_scheme="reduced",
    expected_spin_mode="spin-free",
    rpa=SolverProof(
        status=Status.SUPPORTED,
        steps=(ExecuteConstruct(expected_resolved_scheme="reduced"),),
    ),
    flex=SolverProof(
        status=Status.SUPPORTED,
        steps=(ExecuteConstruct(expected_resolved_scheme="reduced"),),
    ),
    comparison=None,
    required_observables=(),
    interaction_class="onsite",
    notes=(
        "E2 onsite_inter.dat (CoulombInter only), calc_scheme='auto'. "
        "VERIFIED: rpa.py:1285-1289 -- density-diagonal-only "
        "content (no Exchange/PairHop, calc_type='ring') resolves to "
        "'reduced'; both RPA and FLEX (which inherits _set_scheme) "
        "resolve identically. Construction-only: no solve."
    ),
)

_CELL_31_AUTO_EXCHANGE_RESOLUTION = Cell(
    cell_id="auto.exchange.resolution",
    fixture=FixtureSpec(
        interactions={"Exchange": "exchange_onsite.dat"},
        **dict(_E2_GENERAL_MU_KWARGS, requested_scheme="auto"),
    ),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(
        status=Status.SUPPORTED,
        steps=(ExecuteConstruct(expected_resolved_scheme="general"),),
    ),
    flex=SolverProof(
        status=Status.SUPPORTED,
        steps=(ExecuteConstruct(expected_resolved_scheme="general"),),
    ),
    comparison=None,
    required_observables=(),
    interaction_class="onsite",
    notes=(
        "E2 exchange_onsite.dat (Exchange only), calc_scheme='auto'. "
        "VERIFIED: rpa.py:1275-1284 -- Exchange has no "
        "density-diagonal vertex content, so auto resolves to "
        "'general' (both solvers). Construction-only: no solve."
    ),
)

_CELL_32_AUTO_PAIRHOP_RESOLUTION = Cell(
    cell_id="auto.pairhop.resolution",
    fixture=FixtureSpec(
        interactions={"PairHop": "pairhop_onsite.dat"},
        **dict(_E2_GENERAL_MU_KWARGS, requested_scheme="auto"),
    ),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(
        status=Status.SUPPORTED,
        steps=(ExecuteConstruct(expected_resolved_scheme="general"),),
    ),
    flex=SolverProof(
        status=Status.SUPPORTED,
        steps=(ExecuteConstruct(expected_resolved_scheme="general"),),
    ),
    comparison=None,
    required_observables=(),
    interaction_class="onsite",
    notes=(
        "E2 pairhop_onsite.dat (PairHop only), calc_scheme='auto'. "
        "VERIFIED: rpa.py:1275-1284 -- PairHop has no "
        "density-diagonal vertex content either, so auto resolves to "
        "'general' (both solvers). Construction-only: no solve."
    ),
)

_CELL_33_CHI0Q_INIT_REUSE = Cell(
    cell_id="chi0q_init.reuse",
    fixture=FixtureSpec(
        input_dir="tests/equivalence_input/orb1",
        interactions={"CoulombInter": "coulombinter.dat"},
        T=2.0, mu=0.0, filling=None,
        CellShape=(4, 4, 1), SubShape=(1, 1, 1), Nmat=32,
        extra_params={}, calc_type="ring",
        requested_scheme="general", enable_spin_orbital=False, extern=None,
    ),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteChiqInitReuse(),)),
    flex=SolverProof(
        status=Status.NOT_APPLICABLE,
        steps=(PairedInvarianceRun(),),
        reason=(
            "accepted; no corresponding option semantics (flex.py:408) "
            "-- FLEX starts every SCF loop from zero self-energy and "
            "recomputes chi0q from the dressed Green's function each "
            "iteration (flex.py:446-451's own docstring), so a "
            "chi0q_init entry loaded by the inherited RPA read_init is "
            "never consumed."
        ),
    ),
    comparison=None,
    required_observables=(),
    interaction_class="offsite",
    notes=(
        "E1 coulombinter.dat, mu=0.0. RPA side (ExecuteChiqInitReuse): "
        "run A (no chi0q_init) captures chi0q_A/chiq_A; run B loads "
        "chi0q_A via the PUBLIC chi0q_init mechanism (run A's own "
        "save_results() writes the file, run B's read_init() loads "
        "it) -- VERIFIED: np.array_equal(chi0q_B, chi0q_A) "
        "and np.array_equal(chiq_B, chiq_A) both True (chi0q passes "
        "through solve() untouched when supplied). FLEX side "
        "(PairedInvarianceRun): the SAME two-solve construction, "
        "exhaustive green_info key-set match + per-key "
        "np.array_equal -- VERIFIED: identical key sets "
        "({'chi0q','chiq_c','chiq_s','green','physics','sigma'}) and "
        "every array bitwise equal between the with/without-"
        "chi0q_init runs, proving (not merely asserting) FLEX's "
        "output is invariant to the option's presence."
    ),
)

_CELL_34_RINGLADDER_GENERAL_ONSITE_COULOMBINTRA = Cell(
    cell_id="ringladder.general.onsite_coulombintra",
    fixture=FixtureSpec(
        input_dir="tests/equivalence_input/mini",
        interactions={"CoulombIntra": "coulombintra.dat"},
        T=2.0, mu=0.0, filling=None,
        CellShape=(2, 1, 1), SubShape=(1, 1, 1), Nmat=8,
        extra_params={}, calc_type="ring+ladder",
        requested_scheme="general", enable_spin_orbital=False, extern=None,
    ),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(
        status=Status.REJECT,
        steps=(
            ExecuteReject(
                site=Site.CONSTRUCTOR, exc_type="ValueError",
                fragment="ring+ladder",
            ),
        ),
        links=(
            SupplementaryLink(
                test_id=(
                    "tests.test_spinful_transverse_ed::"
                    "TestSpinfulVertexExchangeOptOutRingLadderRejection::"
                    "test_ring_ladder_rejects_the_opt_out_combination"
                ),
                claim=(
                    "The transverse ring+ladder pipeline's own "
                    "rejection contract (a companion instance to "
                    "FLEX's blanket calc_type='ring+ladder' rejection)."
                ),
            ),
            SupplementaryLink(
                test_id=(
                    "tests.test_bond_transverse_w::"
                    "TestGateW1OnsiteReduction::"
                    "test_collapsed_output_matches_plain_chiq_pm_at_omega_zero"
                ),
                claim=(
                    "RPA's ring+ladder transverse channel resolves "
                    "correctly on the collapsed (non-bond) path -- the "
                    "capability FLEX does not implement."
                ),
            ),
        ),
    ),
    comparison=None,
    required_observables=(),
    interaction_class="onsite",
    notes=(
        "E3 (2x1x1, 1 orbital, Nmat=8), coulombintra.dat U=2.0; "
        "mu=0.0, calc_type='ring+ladder' (requires "
        "calc_scheme='general' -- rpa.py:1323-1325). RPA=SUPPORTED "
        "(VERIFIED: solves, green_info carries 'chiq_pm', the "
        "transverse channel). FLEX=REJECT at CONSTRUCTOR: "
        "flex.py:208-219 (_init_flex_param, scheme=='general' branch, "
        "checked BEFORE the enable_spin_orbital guard), 'FLEX does "
        "not support calc_type='ring+ladder' (the transverse ladder "
        "channel); FLEX general is ring-only....' -- VERIFIED "
        "via build_solver's construction path (calc_type "
        "must be threaded to FLEX's info dict too, not just RPA's -- "
        "see build_solver's docstring update)."
    ),
)

# ---------------------------------------------------------------------------
# Group G6 -- conditioning. FC = the one-orbital fixture
# files (tests/equivalence_input/orb1) with T=0.2, Nmat=256, filling=0.5
# overridden on top of that fixture's usual T=2.0/Nmat=32 (the param
# override lives in FixtureSpec, not in new files).
# ---------------------------------------------------------------------------

_CELL_38_CONDITIONING_FC_KWARGS = dict(
    input_dir="tests/equivalence_input/orb1",
    interactions={"CoulombInter": "coulombinter.dat"},
    T=0.2,
    mu=None,
    filling=0.5,
    CellShape=(4, 4, 1),
    SubShape=(1, 1, 1),
    Nmat=256,
    extra_params={},
    calc_type="ring",
    requested_scheme="general",
    enable_spin_orbital=False,
    extern=None,
)

_CELL_38_CONDITIONING_MU = Cell(
    cell_id="general.ring.offsite_coulombinter.conditioning.mu",
    fixture=FixtureSpec(**_CELL_38_CONDITIONING_FC_KWARGS),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    comparison=_measured_equiv(
        "chi0q_mu", 6.431799537609854e-13, 6.432077093507842e-13,
        "general_from_flex_channels", "chiq_mu", 3.3066883022456364e-12, 3.373345155219826e-12,
    ),
    required_observables=("chi0q", "chiq"),
    interaction_class="offsite",
    notes=(
        "FC (Group G6): the one-orbital fixture files "
        "(tests/equivalence_input/orb1) -- geom.dat/transfer.dat "
        "(1 orbital; NN hopping t=1.0 along +-x/+-y plus a "
        "single-diagonal-direction t'=0.5 term at R=+-(1,1)) + "
        "coulombinter.dat (OFF-SITE, same-orbital CoulombInter V=1.0, "
        "R=(+-1,0,0)/(0,+-1,0) -- the SAME file cell 19 uses), with "
        "T=0.2, Nmat=256, filling=0.5 overriding that fixture's usual "
        "T=2.0/Nmat=32. The interaction is off-site, and the cell_id "
        "and interaction_class both say so: a single orbital cannot "
        "carry an on-site inter-orbital coupling, so the off-site "
        "same-orbital file is the only CoulombInter this fixture has. "
        "van-Hove-shoulder verification of that filling: on the "
        "actual CellShape=(4,4,1) 16-k-point mesh, H0's eigenvalues are "
        "[-3, -2,-2,-2,-2, -1,-1,-1,-1, 1,1, 2,2,2,2, 5]; filling=0.5 "
        "(one-spin target 8 of 16 states) lands the Fermi level EXACTLY "
        "inside the 4-fold-degenerate eps=-1.0 cluster -- a coarse-mesh "
        "analogue of a van Hove shoulder (confirmed by direct "
        "H0_eigenvalue inspection; a finer 64x64 k-mesh on the same "
        "dispersion shows the true DOS peak near eps~-0.93, consistent "
        "with this cluster). filling=0.5 is kept as the stated starting "
        "value -- no re-selection was needed (see "
        "tests/equivalence_calibration_log.md). "
        "AMPLIFICATION (owned by TestConditioningAmplification, not this "
        "row): the mu/Green divergence diagnostic run in ISOLATION "
        "(Sigma=0, no solve()) on this fixture vs the benign cell 8 "
        "fixture shows >=10x amplification on the mu-seam residual "
        "(assertion 2) at the very first candidate T=0.2 -- no halve-T "
        "re-choice needed (single logged calibration attempt). "
        "THIS ROW's own (chi0q, chiq) comparison instead uses the FULL "
        "one-shot solve() pipeline, which is NOT the same computation "
        "the diagnostic isolates: FLEX's post-mix stored mu (after one "
        "SCF iteration re-solves mu against its own iteration-1 "
        "self-energy, flex.py:741-743) differs from the diagnostic's "
        "isolated Sigma=0 mu, so a close full-pipeline chi0q/chiq match "
        "here does not contradict the diagnostic's amplified seam-level "
        "residuals in isolation -- exactly the OutputBundle-documented "
        "reason mu is never a per-cell observable. Measured on the macOS "
        "arm64 development machine: chi0q max|diff| = 0.0 "
        "(bit-identical) and chiq "
        "max|diff| = 3.3066883022456364e-12, both far inside their "
        "mu-coupled ceilings. Wall time on the same machine: rpa solve() "
        "~0.15s + flex solve() ~0.01s (~0.16s combined) -- negligible "
        "despite Nmat=256, because this fixture is 1-orbital/16-k-point "
        "(nd=1); no CI-budget threat from this cell."
    ),
)

CELLS: tuple = (
    _BOOTSTRAP_CELL_1,
    _CELL_2_ONSITE_COULOMBINTER_FIXEDMU,
    _CELL_3_ONSITE_HUND_FIXEDMU,
    _CELL_4_ONSITE_ISING_FIXEDMU,
    _CELL_5_ONSITE_EXCHANGE_FIXEDMU,
    _CELL_6_ONSITE_PAIRHOP_FIXEDMU,
    _CELL_7_ONSITE_PAIRLIFT_FIXEDMU,
    _CELL_8_ONSITE_U_V_HUND_MU,
    _CELL_9_ONSITE_FULL_KANAMORI_MU,
    _CELL_10_REDUCED_COULOMBINTRA_SPINFREE_MU,
    _CELL_11_REDUCED_COULOMBINTER_SPINFREE_MU,
    _CELL_12_REDUCED_HUND_SPINFREE_MU,
    _CELL_13_REDUCED_ISING_SPINFREE_MU,
    _CELL_14_REDUCED_PAIRLIFT_SPINFREE_MU,
    _CELL_15_REDUCED_COULOMBINTRA_SPINDIAG_MU,
    _CELL_16_REDUCED_COULOMBINTER_SPINDIAG_MU,
    _BOOTSTRAP_CELL_17,
    _CELL_18_REDUCED_PAIRHOP_REJECT,
    _CELL_19_OFFSITE_COULOMBINTER_SAMEORB_MU,
    _CELL_20_REDUCED_OFFSITE_COULOMBINTER_MU,
    _CELL_21_OFFSITE_COULOMBINTER_INTERORB_FLEXREJECT,
    _CELL_22_OFFSITE_HUND_FLEXREJECT,
    _CELL_23_OFFSITE_ISING_FLEXREJECT,
    _CELL_24_OFFSITE_EXCHANGE_FLEXREJECT,
    _CELL_25_OFFSITE_PAIRHOP_FLEXREJECT,
    _CELL_26_OFFSITE_COULOMBINTRA_LITERALKEY_REJECT,
    _CELL_27_OFFSITE_COULOMBINTER_SAMEORB_SUBSHAPE,
    _CELL_28_ONSITE_COULOMBINTER_SUBSHAPE_MU,
    _CELL_29_ONSITE_COULOMBINTER_COEFFTAIL_MU,
    _CELL_30_AUTO_DENSITY_RESOLUTION,
    _CELL_31_AUTO_EXCHANGE_RESOLUTION,
    _CELL_32_AUTO_PAIRHOP_RESOLUTION,
    _CELL_33_CHI0Q_INIT_REUSE,
    _CELL_34_RINGLADDER_GENERAL_ONSITE_COULOMBINTRA,
    _CELL_35_SO_GENERAL_CONSTRUCTION_REJECT,
    _CELL_36_SO_REDUCED_CONSTRUCTION_REJECT,
    _CELL_38_CONDITIONING_MU,
)


# ---------------------------------------------------------------------------
# COVERAGE_OBLIGATIONS -- named predicates over CELLS: the full list of
# properties the inventory must have, including "the conditioning row
# exists (38)".
# ---------------------------------------------------------------------------


def _cell_ids(cells) -> set:
    return {c.cell_id for c in cells}


def _obligation_g1_every_interaction_type_appears(cells) -> bool:
    """Every interaction type the general scheme supports appears
    somewhere in G1 (cell_ids ``general.ring.onsite_*.fixedmu`` /
    ``.mu``)."""

    required = {
        "CoulombIntra", "CoulombInter", "Hund", "Ising", "Exchange",
        "PairHop", "PairLift",
    }
    seen: set = set()
    for cell in cells:
        if cell.cell_id.startswith("general.ring.onsite_") and (
            cell.cell_id.endswith(".fixedmu") or cell.cell_id.endswith(".mu")
        ):
            seen.update(cell.fixture.interactions.keys())
    return required <= seen


def _obligation_g2_both_reduced_spin_modes_appear(cells) -> bool:
    ids = _cell_ids(cells)
    has_spin_free = any(cid.endswith(".spinfree.mu") for cid in ids)
    has_spin_diag = any(cid.endswith(".spindiag.mu") for cid in ids)
    return has_spin_free and has_spin_diag


def _obligation_every_resolver_outcome_appears(cells) -> bool:
    ids = _cell_ids(cells)
    return {
        "auto.density.resolution",
        "auto.exchange.resolution",
        "auto.pairhop.resolution",
    } <= ids


def _obligation_both_so_guard_sites_appear(cells) -> bool:
    ids = _cell_ids(cells)
    return {
        "so.general.construction.reject",
        "so.reduced.construction.reject",
    } <= ids


def _obligation_full_kanamori_row_exists(cells) -> bool:
    return "general.ring.onsite_full_kanamori.mu" in _cell_ids(cells)


def _obligation_conditioning_row_exists(cells) -> bool:
    return "general.ring.offsite_coulombinter.conditioning.mu" in _cell_ids(cells)


def _obligation_at_least_one_both_reject(cells) -> bool:
    return any(
        cell.rpa.status is Status.REJECT and cell.flex.status is Status.REJECT
        for cell in cells
    )


def _obligation_at_least_one_flex_reject_rpa_supported(cells) -> bool:
    return any(
        cell.rpa.status is Status.SUPPORTED and cell.flex.status is Status.REJECT
        for cell in cells
    )


def _obligation_at_least_one_rpa_only(cells) -> bool:
    return any(
        cell.rpa.status is Status.SUPPORTED
        and cell.flex.status is Status.NOT_APPLICABLE
        for cell in cells
    )


COVERAGE_OBLIGATIONS: Dict[str, Callable[[tuple], bool]] = {
    "g1_every_interaction_type_appears": _obligation_g1_every_interaction_type_appears,
    "g2_both_reduced_spin_modes_appear": _obligation_g2_both_reduced_spin_modes_appear,
    "every_resolver_outcome_appears": _obligation_every_resolver_outcome_appears,
    "both_so_guard_sites_appear": _obligation_both_so_guard_sites_appear,
    "full_kanamori_row_exists": _obligation_full_kanamori_row_exists,
    "conditioning_row_exists": _obligation_conditioning_row_exists,
    "at_least_one_both_reject": _obligation_at_least_one_both_reject,
    "at_least_one_flex_reject_rpa_supported": _obligation_at_least_one_flex_reject_rpa_supported,
    "at_least_one_rpa_only": _obligation_at_least_one_rpa_only,
}


# ---------------------------------------------------------------------------
# method_name
# ---------------------------------------------------------------------------

_METHOD_NAME_TRANSLATION = {ord("."): "_", ord("+"): "_"}


def method_name(cell_id: str) -> str:
    """``"test_cell__" + cell_id`` with ``.`` and ``+`` mapped to ``_``.

    A registry test asserts the mapped names are unique AND valid
    Python identifiers across the whole registry (a duplicate would let
    ``setattr`` silently overwrite a generated proof method).
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

    # Coverage obligations (a no-op while COVERAGE_OBLIGATIONS is
    # empty).
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
# Output bundles and the named comparator registry
# ---------------------------------------------------------------------------
#
# This section is still SIDE-EFFECT-FREE: ``numpy`` is the only import
# it adds. The comparators deliberately live in this REGISTRY module
# (pure numpy, still no unittest and no solver imports), so the docs
# renderer can import them. ``extract_bundle`` -- the function that
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

    NO ``mu`` field: RPA never retains its solved chemical potential as
    a public post-solve attribute (it is a local variable of the solve)
    and FLEX's stored ``mu`` is the post-mix DRESSED value -- a
    different algorithmic state. The mu/Green seam is owned exclusively
    by the divergence diagnostic via ``scalar_residual`` /
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


# The CLOSED comparator registry -- exactly the three keys the cell
# inventory references (no "scalar" key: the diagnostic's scalar mu
# checkpoints use ``scalar_residual``/``assert_scalar_within`` below,
# never a ``COMPARATORS`` entry).
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
    recalibrated back to ``Equiv`` -- a row must always describe what
    the code actually does; the raised message says so explicitly.
    INCLUSIVE upper edge: ``residual <= regression_bound``, a
    further-drift regression guard.
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
