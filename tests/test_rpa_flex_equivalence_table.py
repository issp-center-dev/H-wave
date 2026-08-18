"""Validation tests for the RPA/FLEX equivalence-table cell registry
(``tests/equivalence_cells.py``, Task 1 of the equivalence-table plan).

This module tests the REGISTRY MACHINERY -- the typed records, the
schema validators (``validate_registry``), ``relationship()``'s total
status-pair lookup, and ``method_name`` -- NOT any cell content
(``CELLS``/``COVERAGE_OBLIGATIONS`` are intentionally empty until Task
5). Every test constructs its own minimal ``Cell``/``SolverProof``
fixtures via the small builders below.
"""

from __future__ import annotations

import unittest
from types import MappingProxyType

import numpy as np

from tests.equivalence_cells import (
    CELLS,
    COMPARATORS,
    COVERAGE_OBLIGATIONS,
    POLICY_CEILINGS,
    PROOF_STEP_TYPES,
    PROVENANCE,
    Cell,
    Comparator,
    DivergingSpec,
    Diverges,
    Equiv,
    ExecuteChiqInitReuse,
    ExecuteConstruct,
    ExecuteReject,
    FixtureSpec,
    ObservableSpec,
    OutputBundle,
    PairedInvarianceRun,
    ExecuteRun,
    Site,
    SolverProof,
    Status,
    SupplementaryLink,
    assert_diverges_bracket,
    assert_scalar_within,
    method_name,
    relationship,
    scalar_residual,
    validate_registry,
)


# ---------------------------------------------------------------------------
# extract_bundle -- TEST-module side (it touches solver objects): builds
# an ``OutputBundle`` from one solver's post-solve ``green_info``. The
# registry module (``tests/equivalence_cells.py``) never imports solver
# code; this function is the seam where a concrete ``solver_obj``/
# ``green_info`` pair (Task 3's ``build_solver``, or -- for these unit
# tests -- a small stand-in) becomes the comparator-facing bundle type.
# ---------------------------------------------------------------------------


def extract_bundle(solver_obj, green_info, solver_kind: str) -> OutputBundle:
    """Build the ``OutputBundle`` a comparator maps over.

    ``solver_kind`` selects the extraction shape:

    * ``"rpa"`` reads ``green_info["chi0q"]``/``green_info["chiq"]`` --
      RPA never populates the channel-decomposed keys
      (``src/hwave/solver/rpa.py:2110`` sets ``chiq``; there is no
      ``chiq_s``/``chiq_c`` write anywhere in ``rpa.py``), so those two
      bundle fields are left ``None``.
    * ``"flex"`` reads ``green_info["chi0q"]``/``green_info["chiq_s"]``/
      ``green_info["chiq_c"]`` -- FLEX never populates a combined
      ``"chiq"`` key (``src/hwave/solver/flex.py:746-748`` sets exactly
      ``chi0q``/``chiq_s``/``chiq_c``), so ``chiq`` is left ``None``.

    ``solver_obj`` is accepted for parity with the run-time callers
    (Task 3's ``build_solver`` produces the ``(solver_obj, green_info)``
    pair together) and is not otherwise inspected here -- every value
    this function returns comes from ``green_info``.
    """

    if solver_kind == "rpa":
        return OutputBundle(
            chi0q=np.asarray(green_info["chi0q"]),
            chiq=np.asarray(green_info["chiq"]),
            chiq_s=None,
            chiq_c=None,
        )
    if solver_kind == "flex":
        return OutputBundle(
            chi0q=np.asarray(green_info["chi0q"]),
            chiq=None,
            chiq_s=np.asarray(green_info["chiq_s"]),
            chiq_c=np.asarray(green_info["chiq_c"]),
        )
    raise ValueError(
        "extract_bundle: solver_kind must be 'rpa' or 'flex', got "
        "{!r}".format(solver_kind)
    )


# ---------------------------------------------------------------------------
# Small builders -- keep individual tests focused on the one thing they
# violate.
# ---------------------------------------------------------------------------


def _fixture(**overrides):
    defaults = dict(
        input_dir="tests/equivalence_input/orb1",
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
    )
    defaults.update(overrides)
    return FixtureSpec(**defaults)


def _run_proof(status=Status.SUPPORTED, step=None, links=(), reason=""):
    step = ExecuteRun() if step is None else step
    return SolverProof(status=status, steps=(step,), links=links, reason=reason)


def _reject_proof(site=Site.CONSTRUCTOR, exc_type="ValueError", fragment="reduced", links=()):
    return SolverProof(
        status=Status.REJECT,
        steps=(ExecuteReject(site=site, exc_type=exc_type, fragment=fragment),),
        links=links,
    )


def _equiv_observables(mu_mode="fixed"):
    suffix = "_mu" if mu_mode == "mu" else "_fixed"
    return {
        "chi0q": ObservableSpec(
            comparator="identity", atol=POLICY_CEILINGS["chi0q" + suffix], provenance="measured"
        ),
        "chiq": ObservableSpec(
            comparator="general_from_flex_channels",
            atol=POLICY_CEILINGS["chiq" + suffix],
            provenance="measured",
        ),
    }


def _equiv_cell(cell_id="equiv.cell", mu_mode="fixed", **overrides):
    fixture = overrides.pop("fixture", None)
    if fixture is None:
        fixture = _fixture(mu=0.0, filling=None) if mu_mode == "fixed" else _fixture(mu=None, filling=0.5)
    comparison = overrides.pop("comparison", Equiv(observables=_equiv_observables(mu_mode)))
    defaults = dict(
        cell_id=cell_id,
        fixture=fixture,
        resolved_scheme="general",
        expected_spin_mode="spin-free",
        rpa=_run_proof(),
        flex=_run_proof(),
        comparison=comparison,
        required_observables=("chi0q", "chiq"),
        interaction_class="onsite",
        notes="",
    )
    defaults.update(overrides)
    return Cell(**defaults)


def _diverges_cell(cell_id="diverges.cell", mu_mode="fixed", step5_issue="#123", **overrides):
    fixture = overrides.pop("fixture", None)
    if fixture is None:
        fixture = _fixture(mu=0.0, filling=None) if mu_mode == "fixed" else _fixture(mu=None, filling=0.5)
    suffix = "_mu" if mu_mode == "mu" else "_fixed"
    ceiling = POLICY_CEILINGS["chiq" + suffix]
    comparison = overrides.pop(
        "comparison",
        Diverges(
            diverging={
                "chiq": DivergingSpec(
                    comparator="general_from_flex_channels",
                    ceiling=ceiling,
                    regression_bound=ceiling * 10,
                    provenance="measured",
                )
            },
            others={
                "chi0q": ObservableSpec(
                    comparator="identity",
                    atol=POLICY_CEILINGS["chi0q" + suffix],
                    provenance="measured",
                )
            },
            step5_issue=step5_issue,
        ),
    )
    defaults = dict(
        cell_id=cell_id,
        fixture=fixture,
        resolved_scheme="general",
        expected_spin_mode="spin-free",
        rpa=_run_proof(),
        flex=_run_proof(),
        comparison=comparison,
        required_observables=("chi0q", "chiq"),
        interaction_class="onsite",
        notes="",
    )
    defaults.update(overrides)
    return Cell(**defaults)


def _reject_reject_cell(cell_id="both.reject.cell", **overrides):
    defaults = dict(
        cell_id=cell_id,
        fixture=_fixture(),
        resolved_scheme="reduced",
        expected_spin_mode="spin-free",
        rpa=_reject_proof(),
        flex=_reject_proof(),
        comparison=None,
        required_observables=(),
        interaction_class="onsite",
        notes="",
    )
    defaults.update(overrides)
    return Cell(**defaults)


def _construction_cell(cell_id="auto.resolution.cell", scheme="reduced", **overrides):
    defaults = dict(
        cell_id=cell_id,
        fixture=_fixture(requested_scheme="auto"),
        resolved_scheme=scheme,
        expected_spin_mode="spin-free",
        rpa=_run_proof(step=ExecuteConstruct(expected_resolved_scheme=scheme)),
        flex=_run_proof(step=ExecuteConstruct(expected_resolved_scheme=scheme)),
        comparison=None,
        required_observables=(),
        interaction_class="onsite",
        notes="",
    )
    defaults.update(overrides)
    return Cell(**defaults)


def _not_applicable_narrowed_cell(cell_id="na.narrowed.cell", **overrides):
    defaults = dict(
        cell_id=cell_id,
        fixture=_fixture(),
        resolved_scheme="general",
        expected_spin_mode="spin-free",
        rpa=_run_proof(),
        flex=SolverProof(
            status=Status.NOT_APPLICABLE,
            steps=(),
            reason="accepted; no corresponding option semantics (flex.py:408)",
        ),
        comparison=None,
        required_observables=(),
        interaction_class="onsite",
        notes="",
    )
    defaults.update(overrides)
    return Cell(**defaults)


def _not_applicable_paired_cell(cell_id="na.paired.cell", **overrides):
    defaults = dict(
        cell_id=cell_id,
        fixture=_fixture(),
        resolved_scheme="general",
        expected_spin_mode="spin-free",
        rpa=SolverProof(
            status=Status.NOT_APPLICABLE,
            steps=(PairedInvarianceRun(),),
            reason="invariance under chi0q_init reuse",
        ),
        flex=_run_proof(),
        comparison=None,
        required_observables=(),
        interaction_class="onsite",
        notes="",
    )
    defaults.update(overrides)
    return Cell(**defaults)


class TestRegistrySchema(unittest.TestCase):
    """One test per violation class enforced by ``validate_registry``,
    plus the multiplicity/uniqueness rules and a positive control.
    """

    # -- module-level constants -------------------------------------------

    def test_policy_ceilings_values(self):
        self.assertEqual(
            dict(POLICY_CEILINGS),
            {
                "mu_diag": 1e-10,
                "green_diag": 1e-10,
                "chi0q_mu": 1e-10,
                "chiq_mu": 1e-10,
                "chi0q_fixed": 1e-12,
                "chiq_fixed": 1e-12,
            },
        )

    def test_provenance_initial_state(self):
        self.assertEqual(
            dict(PROVENANCE),
            {"source_sha": None, "run_ids": (), "status": "candidate"},
        )

    def test_cells_and_coverage_obligations_empty_this_task(self):
        self.assertEqual(CELLS, ())
        self.assertEqual(dict(COVERAGE_OBLIGATIONS), {})

    def test_proof_step_types_is_the_closed_five_member_union(self):
        self.assertEqual(
            set(PROOF_STEP_TYPES),
            {ExecuteRun, ExecuteConstruct, ExecuteReject, PairedInvarianceRun, ExecuteChiqInitReuse},
        )

    # -- method_name --------------------------------------------------------

    def test_method_name_translates_dot_and_plus(self):
        self.assertEqual(
            method_name("general.ring.onsite_u_v_hund.mu"),
            "test_cell__general_ring_onsite_u_v_hund_mu",
        )
        self.assertEqual(method_name("auto+density.resolution"), "test_cell__auto_density_resolution")
        self.assertTrue(method_name("general.ring.onsite_ising.fixedmu").isidentifier())

    # -- immutability ---------------------------------------------------------

    def test_nested_mappings_are_frozen(self):
        fixture = _fixture(interactions={"CoulombIntra": "coulombintra.dat"}, extra_params={"coeff_tail": 1.0})
        self.assertIsInstance(fixture.interactions, MappingProxyType)
        self.assertIsInstance(fixture.extra_params, MappingProxyType)
        with self.assertRaises(TypeError):
            fixture.interactions["CoulombIntra"] = "other.dat"

        equiv = Equiv(observables=_equiv_observables("fixed"))
        self.assertIsInstance(equiv.observables, MappingProxyType)
        with self.assertRaises(TypeError):
            equiv.observables["chiq"] = None

        diverges = _diverges_cell().comparison
        self.assertIsInstance(diverges.diverging, MappingProxyType)
        self.assertIsInstance(diverges.others, MappingProxyType)
        with self.assertRaises(TypeError):
            diverges.diverging["chiq"] = None

    # -- closed step union ----------------------------------------------------

    def test_unknown_step_type_is_rejected(self):
        class RogueStep:
            pass

        cell = _equiv_cell(rpa=SolverProof(status=Status.SUPPORTED, steps=(RogueStep(),)))
        errors = validate_registry([cell])
        self.assertTrue(any("closed PROOF_STEP_TYPES union" in e for e in errors), errors)

    # -- SUPPORTED composition -------------------------------------------------

    def test_supported_rejects_more_than_one_step(self):
        cell = _equiv_cell(rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(), ExecuteRun())))
        errors = validate_registry([cell])
        self.assertTrue(any("SUPPORTED requires exactly one run-class step" in e for e in errors), errors)

    def test_supported_rejects_a_non_run_class_step(self):
        cell = _equiv_cell(
            rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteReject(Site.SOLVE, "ValueError", "x"),))
        )
        errors = validate_registry([cell])
        self.assertTrue(any("SUPPORTED requires exactly one run-class step" in e for e in errors), errors)

    # -- REJECT composition -----------------------------------------------------

    def test_reject_rejects_zero_steps(self):
        cell = _reject_reject_cell(rpa=SolverProof(status=Status.REJECT, steps=()))
        errors = validate_registry([cell])
        self.assertTrue(any("REJECT requires exactly one ExecuteReject" in e for e in errors), errors)

    def test_reject_rejects_a_non_reject_step(self):
        cell = _reject_reject_cell(rpa=SolverProof(status=Status.REJECT, steps=(ExecuteRun(),)))
        errors = validate_registry([cell])
        self.assertTrue(any("REJECT requires exactly one ExecuteReject" in e for e in errors), errors)

    def test_reject_exc_type_must_be_in_closed_set(self):
        cell = _reject_reject_cell(rpa=_reject_proof(exc_type="TypeError"))
        errors = validate_registry([cell])
        self.assertTrue(any("is not in the closed set" in e for e in errors), errors)

    # -- NOT_APPLICABLE composition ------------------------------------------

    def test_not_applicable_zero_steps_requires_narrowed_reason_phrase(self):
        cell = _not_applicable_narrowed_cell(
            flex=SolverProof(status=Status.NOT_APPLICABLE, steps=(), reason="accepted, ignored")
        )
        errors = validate_registry([cell])
        self.assertTrue(any("'no corresponding'" in e for e in errors), errors)

    def test_not_applicable_zero_steps_with_narrowed_phrase_is_valid(self):
        cell = _not_applicable_narrowed_cell()
        errors = validate_registry([cell])
        self.assertEqual(errors, [])

    def test_not_applicable_single_step_must_be_paired_invariance_run(self):
        cell = _not_applicable_paired_cell(
            rpa=SolverProof(status=Status.NOT_APPLICABLE, steps=(ExecuteRun(),), reason="x")
        )
        errors = validate_registry([cell])
        self.assertTrue(
            any("NOT_APPLICABLE requires exactly one PairedInvarianceRun step or zero steps" in e for e in errors),
            errors,
        )

    def test_not_applicable_paired_invariance_run_is_valid(self):
        cell = _not_applicable_paired_cell()
        errors = validate_registry([cell])
        self.assertEqual(errors, [])

    def test_reason_mandatory_for_not_applicable(self):
        cell = _not_applicable_narrowed_cell(
            flex=SolverProof(status=Status.NOT_APPLICABLE, steps=(), reason="")
        )
        errors = validate_registry([cell])
        self.assertTrue(any("reason is mandatory when status is NOT_APPLICABLE" in e for e in errors), errors)

    def test_reason_forbidden_outside_not_applicable(self):
        cell = _equiv_cell(rpa=_run_proof(reason="should not be here"))
        errors = validate_registry([cell])
        self.assertTrue(any("reason must be empty unless status is NOT_APPLICABLE" in e for e in errors), errors)

    def test_not_applicable_pair_rejected_outright(self):
        cell = _reject_reject_cell(
            rpa=SolverProof(status=Status.NOT_APPLICABLE, steps=(), reason="no corresponding thing"),
            flex=SolverProof(status=Status.NOT_APPLICABLE, steps=(), reason="no corresponding thing"),
        )
        errors = validate_registry([cell])
        self.assertTrue(any("(NOT_APPLICABLE, NOT_APPLICABLE) is not a valid status pair" in e for e in errors), errors)

    # -- fixture invariants -----------------------------------------------------

    def test_mu_and_filling_are_mutually_exclusive(self):
        cell = _equiv_cell(fixture=_fixture(mu=0.0, filling=0.5))
        errors = validate_registry([cell])
        self.assertTrue(any("mutually exclusive" in e for e in errors), errors)

    # -- required_observables / comparison presence ------------------------

    def test_required_observables_must_be_empty_when_comparison_is_none(self):
        cell = _reject_reject_cell(required_observables=("chi0q",))
        errors = validate_registry([cell])
        self.assertTrue(any("required_observables must be empty iff comparison is None" in e for e in errors), errors)

    def test_required_observables_must_be_nonempty_when_comparison_present(self):
        cell = _equiv_cell(required_observables=())
        errors = validate_registry([cell])
        self.assertTrue(any("required_observables must be empty iff comparison is None" in e for e in errors), errors)

    def test_comparison_forbidden_without_both_sides_solve_running(self):
        cell = _equiv_cell(flex=_reject_proof(), required_observables=("chi0q", "chiq"))
        errors = validate_registry([cell])
        self.assertTrue(any("comparison must be present iff both solver proofs" in e for e in errors), errors)

    def test_comparison_required_when_both_sides_solve_running(self):
        cell = _equiv_cell(comparison=None, required_observables=())
        errors = validate_registry([cell])
        self.assertTrue(any("comparison must be present iff both solver proofs" in e for e in errors), errors)

    def test_construction_only_rows_forbid_a_comparison(self):
        cell = _construction_cell(
            comparison=Equiv(observables=_equiv_observables("fixed")),
            required_observables=("chi0q", "chiq"),
        )
        errors = validate_registry([cell])
        self.assertTrue(any("comparison=None and required_observables=()" in e for e in errors), errors)

    # -- construction-only scheme agreement ---------------------------------

    def test_construction_only_requires_matching_expected_schemes(self):
        cell = _construction_cell(
            rpa=_run_proof(step=ExecuteConstruct(expected_resolved_scheme="reduced")),
            flex=_run_proof(step=ExecuteConstruct(expected_resolved_scheme="general")),
        )
        errors = validate_registry([cell])
        self.assertTrue(any("expected_resolved_scheme must agree on both sides" in e for e in errors), errors)

    def test_construction_only_valid_row_has_no_errors(self):
        cell = _construction_cell()
        errors = validate_registry([cell])
        self.assertEqual(errors, [])

    # -- Equiv completeness / ceilings ---------------------------------------

    def test_equiv_observable_keys_must_match_required_observables(self):
        cell = _equiv_cell(required_observables=("chi0q",))
        errors = validate_registry([cell])
        self.assertTrue(any("Equiv.observables keys" in e for e in errors), errors)

    def test_equiv_atol_may_equal_the_ceiling(self):
        cell = _equiv_cell(mu_mode="fixed")
        errors = validate_registry([cell])
        self.assertEqual(errors, [])

    def test_equiv_atol_above_ceiling_is_rejected(self):
        observables = dict(_equiv_observables("fixed"))
        observables["chiq"] = ObservableSpec(comparator="general_from_flex_channels", atol=1e-6, provenance="measured")
        cell = _equiv_cell(comparison=Equiv(observables=observables))
        errors = validate_registry([cell])
        self.assertTrue(any("exceeds its policy ceiling" in e for e in errors), errors)

    # -- Diverges completeness / ceilings -------------------------------------

    def test_diverges_valid_row_has_no_errors(self):
        cell = _diverges_cell()
        errors = validate_registry([cell])
        self.assertEqual(errors, [])

    def test_diverges_diverging_must_be_nonempty(self):
        cell = _diverges_cell(
            comparison=Diverges(
                diverging={},
                others={
                    "chi0q": ObservableSpec("identity", POLICY_CEILINGS["chi0q_fixed"], "measured"),
                    "chiq": ObservableSpec("general_from_flex_channels", POLICY_CEILINGS["chiq_fixed"], "measured"),
                },
                step5_issue="#1",
            )
        )
        errors = validate_registry([cell])
        self.assertTrue(any("Diverges.diverging must be non-empty" in e for e in errors), errors)

    def test_diverges_diverging_and_others_must_be_disjoint(self):
        ceiling = POLICY_CEILINGS["chiq_fixed"]
        cell = _diverges_cell(
            comparison=Diverges(
                diverging={"chiq": DivergingSpec("general_from_flex_channels", ceiling, ceiling * 10, "measured")},
                others={"chiq": ObservableSpec("general_from_flex_channels", ceiling, "measured")},
                step5_issue="#1",
            ),
            required_observables=("chiq",),
        )
        errors = validate_registry([cell])
        self.assertTrue(any("must be disjoint" in e for e in errors), errors)

    def test_diverges_union_must_equal_required_observables(self):
        cell = _diverges_cell(required_observables=("chiq",))
        errors = validate_registry([cell])
        self.assertTrue(any("must equal required_observables" in e for e in errors), errors)

    def test_diverges_ceiling_must_equal_mapped_policy_ceiling(self):
        cell = _diverges_cell(
            comparison=Diverges(
                diverging={"chiq": DivergingSpec("general_from_flex_channels", 1e-6, 1e-5, "measured")},
                others={"chi0q": ObservableSpec("identity", POLICY_CEILINGS["chi0q_fixed"], "measured")},
                step5_issue="#1",
            )
        )
        errors = validate_registry([cell])
        self.assertTrue(any("must equal the mapped policy ceiling" in e for e in errors), errors)

    def test_diverges_ceiling_must_be_strictly_below_regression_bound(self):
        ceiling = POLICY_CEILINGS["chiq_fixed"]
        cell = _diverges_cell(
            comparison=Diverges(
                diverging={"chiq": DivergingSpec("general_from_flex_channels", ceiling, ceiling, "measured")},
                others={"chi0q": ObservableSpec("identity", POLICY_CEILINGS["chi0q_fixed"], "measured")},
                step5_issue="#1",
            )
        )
        errors = validate_registry([cell])
        self.assertTrue(any("ceiling < regression_bound" in e for e in errors), errors)

    def test_diverges_step5_issue_must_match_issue_regex(self):
        for bad in ("1234", "#0", "#01", "#issue", ""):
            with self.subTest(step5_issue=bad):
                cell = _diverges_cell(cell_id="diverges.bad." + repr(bad), step5_issue=bad)
                errors = validate_registry([cell])
                self.assertTrue(any("step5_issue" in e for e in errors), errors)

    # -- SupplementaryLink format / duplicates -------------------------------

    def test_supplementary_link_format_is_validated(self):
        cell = _equiv_cell(
            rpa=_run_proof(links=(SupplementaryLink(test_id="not-a-valid-link", claim="x"),))
        )
        errors = validate_registry([cell])
        self.assertTrue(any("does not match" in e for e in errors), errors)

    def test_supplementary_link_valid_format_has_no_error(self):
        cell = _equiv_cell(
            rpa=_run_proof(
                links=(
                    SupplementaryLink(
                        test_id="tests.test_flex_offsite_general::TestFlexOffsiteGeneral::test_reject",
                        claim="rejects off-site general",
                    ),
                )
            )
        )
        errors = validate_registry([cell])
        self.assertEqual(errors, [])

    def test_duplicate_supplementary_link_within_a_cell_is_rejected(self):
        link = SupplementaryLink(test_id="tests.mod::TestCase::test_x", claim="x")
        cell = _equiv_cell(rpa=_run_proof(links=(link, link)))
        errors = validate_registry([cell])
        self.assertTrue(any("duplicate SupplementaryLink.test_id" in e for e in errors), errors)

    # -- uniqueness -------------------------------------------------------------

    def test_duplicate_cell_id_is_rejected(self):
        errors = validate_registry([_equiv_cell(cell_id="dup"), _equiv_cell(cell_id="dup")])
        self.assertTrue(any("duplicate cell_id" in e for e in errors), errors)

    def test_duplicate_method_name_from_distinct_cell_ids_is_rejected(self):
        errors = validate_registry(
            [_equiv_cell(cell_id="a.b"), _equiv_cell(cell_id="a+b")]
        )
        self.assertTrue(any("duplicate method_name" in e for e in errors), errors)

    def test_method_name_that_is_not_a_valid_identifier_is_rejected(self):
        # "-" survives method_name's translate() (only "." and "+" are
        # mapped to "_"), so a cell_id containing a hyphen yields a
        # non-identifier method name (e.g. "test_cell__bad-cell_id").
        cell = _equiv_cell(cell_id="bad-cell.id")
        errors = validate_registry([cell])
        self.assertTrue(
            any("is not a valid Python identifier" in e for e in errors), errors
        )

    # -- positive control ---------------------------------------------------

    def test_a_mixed_valid_registry_has_no_errors(self):
        cells = [
            _equiv_cell(cell_id="g1.equiv", mu_mode="fixed"),
            _equiv_cell(cell_id="g1.equiv.mu", mu_mode="mu"),
            _diverges_cell(cell_id="g1.diverges"),
            _reject_reject_cell(cell_id="g2.bothreject"),
            _construction_cell(cell_id="g4.auto"),
            _not_applicable_narrowed_cell(cell_id="g4.chi0q_init.narrowed"),
            _not_applicable_paired_cell(cell_id="g4.chi0q_init.paired"),
        ]
        errors = validate_registry(cells)
        self.assertEqual(errors, [])


class TestRelationship(unittest.TestCase):
    """``relationship()`` over the full ``(rpa.status, flex.status)``
    matrix, plus the construction-only ``AUTO-RESOLVES`` variant.
    """

    def test_supported_supported_equiv(self):
        self.assertEqual(relationship(_equiv_cell()), "EQUIV")

    def test_supported_supported_diverges(self):
        self.assertEqual(relationship(_diverges_cell(step5_issue="#42")), "DIVERGES(#42)")

    def test_supported_supported_without_comparison_raises(self):
        cell = _equiv_cell(comparison=None, required_observables=())
        with self.assertRaises(ValueError):
            relationship(cell)

    def test_supported_reject(self):
        cell = _equiv_cell(flex=_reject_proof(), comparison=None, required_observables=())
        self.assertEqual(relationship(cell), "FLEX-REJECT·RPA-SUPPORTED")

    def test_supported_not_applicable(self):
        cell = _not_applicable_narrowed_cell()
        self.assertEqual(relationship(cell), "RPA-ONLY")

    def test_reject_supported(self):
        cell = _equiv_cell(rpa=_reject_proof(), comparison=None, required_observables=())
        self.assertEqual(relationship(cell), "RPA-REJECT·FLEX-SUPPORTED")

    def test_reject_reject(self):
        self.assertEqual(relationship(_reject_reject_cell()), "BOTH-REJECT")

    def test_reject_not_applicable(self):
        cell = _reject_reject_cell(
            flex=SolverProof(status=Status.NOT_APPLICABLE, steps=(), reason="no corresponding option")
        )
        self.assertEqual(relationship(cell), "RPA-REJECT")

    def test_not_applicable_supported(self):
        cell = _not_applicable_paired_cell()
        self.assertEqual(relationship(cell), "FLEX-ONLY")

    def test_not_applicable_reject(self):
        cell = _reject_reject_cell(
            rpa=SolverProof(status=Status.NOT_APPLICABLE, steps=(), reason="no corresponding option")
        )
        self.assertEqual(relationship(cell), "FLEX-REJECT")

    def test_not_applicable_not_applicable_raises(self):
        cell = _reject_reject_cell(
            rpa=SolverProof(status=Status.NOT_APPLICABLE, steps=(), reason="no corresponding option"),
            flex=SolverProof(status=Status.NOT_APPLICABLE, steps=(), reason="no corresponding option"),
        )
        with self.assertRaises(ValueError):
            relationship(cell)

    def test_construction_only_auto_resolves(self):
        cell = _construction_cell(scheme="reduced")
        self.assertEqual(relationship(cell), "AUTO-RESOLVES(reduced)")

    def test_construction_only_mismatched_schemes_falls_back_and_raises(self):
        cell = _construction_cell(
            rpa=_run_proof(step=ExecuteConstruct(expected_resolved_scheme="reduced")),
            flex=_run_proof(step=ExecuteConstruct(expected_resolved_scheme="general")),
        )
        with self.assertRaises(ValueError):
            relationship(cell)


class TestExtractBundle(unittest.TestCase):
    """``extract_bundle`` -- fake ``green_info`` dicts, no real solver."""

    def test_rpa_bundle_reads_chi0q_and_chiq_leaves_channels_none(self):
        green_info = {"chi0q": np.array([[1.0 + 2.0j]]), "chiq": np.array([[3.0 - 4.0j]])}
        bundle = extract_bundle(object(), green_info, "rpa")
        np.testing.assert_array_equal(bundle.chi0q, green_info["chi0q"])
        np.testing.assert_array_equal(bundle.chiq, green_info["chiq"])
        self.assertIsNone(bundle.chiq_s)
        self.assertIsNone(bundle.chiq_c)

    def test_flex_bundle_reads_channels_leaves_chiq_none(self):
        green_info = {
            "chi0q": np.array([[1.0 + 2.0j]]),
            "chiq_s": np.array([[5.0 + 6.0j]]),
            "chiq_c": np.array([[7.0 - 8.0j]]),
        }
        bundle = extract_bundle(object(), green_info, "flex")
        np.testing.assert_array_equal(bundle.chi0q, green_info["chi0q"])
        np.testing.assert_array_equal(bundle.chiq_s, green_info["chiq_s"])
        np.testing.assert_array_equal(bundle.chiq_c, green_info["chiq_c"])
        self.assertIsNone(bundle.chiq)

    def test_unknown_solver_kind_raises(self):
        with self.assertRaises(ValueError):
            extract_bundle(object(), {"chi0q": np.zeros((1,))}, "flux")


class TestScalarUtilities(unittest.TestCase):
    """``scalar_residual`` / ``assert_scalar_within`` -- the Task-6
    diagnostic's mu checkpoints (never a COMPARATORS entry).
    """

    def test_scalar_residual_exact_match(self):
        self.assertEqual(scalar_residual(2.0, 2.0), 0.0)

    def test_scalar_residual_computes_abs_diff(self):
        self.assertEqual(scalar_residual(2.0, 5.0), 3.0)
        self.assertEqual(scalar_residual(5.0, 2.0), 3.0)

    def test_scalar_residual_rejects_nan(self):
        with self.assertRaises(ValueError):
            scalar_residual(float("nan"), 1.0)
        with self.assertRaises(ValueError):
            scalar_residual(1.0, float("inf"))

    def test_assert_scalar_within_passes_below_atol(self):
        assert_scalar_within(1.0, 1.2, atol=0.5)  # must not raise

    def test_assert_scalar_within_passes_at_exact_atol_boundary(self):
        assert_scalar_within(1.0, 1.5, atol=0.5)  # equality passes (<=)

    def test_assert_scalar_within_fails_above_atol(self):
        with self.assertRaises(AssertionError):
            assert_scalar_within(1.0, 1.51, atol=0.5)


class TestDivergesBracket(unittest.TestCase):
    """``assert_diverges_bracket`` -- the two-sided ``(ceiling,
    regression_bound]`` semantics at the assertion-helper level.
    """

    CEILING = 1e-10
    BOUND = 1e-9

    def test_equality_at_ceiling_fails_the_strict_lower_check(self):
        with self.assertRaises(AssertionError) as ctx:
            assert_diverges_bracket(self.CEILING, self.CEILING, self.BOUND)
        self.assertIn("recalibrate", str(ctx.exception))

    def test_below_ceiling_also_heals_and_asks_to_recalibrate(self):
        with self.assertRaises(AssertionError) as ctx:
            assert_diverges_bracket(self.CEILING * 0.5, self.CEILING, self.BOUND)
        self.assertIn("recalibrate", str(ctx.exception))

    def test_strictly_above_ceiling_and_below_bound_passes(self):
        assert_diverges_bracket(self.CEILING * 5, self.CEILING, self.BOUND)  # must not raise

    def test_equality_at_regression_bound_passes(self):
        assert_diverges_bracket(self.BOUND, self.CEILING, self.BOUND)  # must not raise

    def test_above_regression_bound_fails(self):
        with self.assertRaises(AssertionError) as ctx:
            assert_diverges_bracket(self.BOUND * 1.0001, self.CEILING, self.BOUND)
        self.assertIn("regression_bound", str(ctx.exception))


class TestComparatorRegistryShape(unittest.TestCase):
    """``COMPARATORS`` has exactly the three keys the plan's Appendix A
    cell matrix references -- no ``"scalar"`` key.
    """

    def test_exact_key_set(self):
        self.assertEqual(
            set(COMPARATORS.keys()),
            {"identity", "general_from_flex_channels", "reduced_blocks"},
        )

    def test_every_value_is_a_comparator(self):
        for name, comparator in COMPARATORS.items():
            with self.subTest(name=name):
                self.assertIsInstance(comparator, Comparator)


class TestIdentityComparator(unittest.TestCase):
    """``identity`` on chi0q -- hand-enumerated small asymmetric arrays."""

    def _bundle(self, chi0q):
        return OutputBundle(chi0q=chi0q, chiq=None, chiq_s=None, chiq_c=None)

    def test_map_returns_the_raw_chi0q_pair(self):
        rpa_chi0q = np.array([[1.0 + 2.0j, -3.0 + 0.5j], [0.25 - 1.5j, 7.0 + 7.0j]])
        flex_chi0q = np.array([[1.0 + 2.0j, -3.0 + 0.5j], [0.25 - 1.5j, 7.0 + 7.0j]])
        comparator = COMPARATORS["identity"]
        a, b = comparator.map("chi0q", self._bundle(rpa_chi0q), self._bundle(flex_chi0q))
        np.testing.assert_array_equal(a, rpa_chi0q)
        np.testing.assert_array_equal(b, flex_chi0q)
        self.assertEqual(comparator.residual(a, b), 0.0)
        comparator.assert_within(a, b, atol=0.0)  # must not raise

    def test_max_diff_value_and_index_are_reported_on_failure(self):
        rpa_chi0q = np.array([[1.0 + 2.0j, -3.0 + 0.5j], [0.25 - 1.5j, 7.0 + 7.0j]])
        flex_chi0q = rpa_chi0q.copy()
        flex_chi0q[1, 0] = rpa_chi0q[1, 0] + 0.75  # a single, known perturbation
        comparator = COMPARATORS["identity"]
        a, b = comparator.map("chi0q", self._bundle(rpa_chi0q), self._bundle(flex_chi0q))
        self.assertAlmostEqual(comparator.residual(a, b), 0.75)
        with self.assertRaises(AssertionError) as ctx:
            comparator.assert_within(a, b, atol=0.5)
        message = str(ctx.exception)
        self.assertIn("0.75", message)
        self.assertIn("(1, 0)", message)
        comparator.assert_within(a, b, atol=0.75)  # equality passes (<=)

    def test_shape_mismatch_is_rejected(self):
        comparator = COMPARATORS["identity"]
        a = np.zeros((2, 2), dtype=complex)
        b = np.zeros((2, 3), dtype=complex)
        with self.assertRaises(ValueError):
            comparator.residual(a, b)
        with self.assertRaises(ValueError):
            comparator.assert_within(a, b, atol=1.0)

    def test_nan_is_rejected(self):
        comparator = COMPARATORS["identity"]
        a = np.array([[1.0 + 0.0j, float("nan") + 0.0j]])
        b = np.array([[1.0 + 0.0j, 0.0 + 0.0j]])
        with self.assertRaises(ValueError):
            comparator.residual(a, b)
        with self.assertRaises(ValueError):
            comparator.assert_within(a, b, atol=1.0)

    def test_unset_bundle_field_is_rejected(self):
        comparator = COMPARATORS["identity"]
        rpa_bundle = OutputBundle(chi0q=np.zeros((1, 1)), chiq=None, chiq_s=None, chiq_c=None)
        flex_bundle = OutputBundle(chi0q=None, chiq=None, chiq_s=None, chiq_c=None)
        with self.assertRaises(ValueError):
            comparator.map("chi0q", rpa_bundle, flex_bundle)


class TestGeneralFromFlexChannelsComparator(unittest.TestCase):
    """``general_from_flex_channels`` -- the 0.5*(cc+-cs) same/diff
    reconstruction, verbatim from ``tests/test_rpa_flex_oneshot_
    equivalence.py:88-99``, on a hand-enumerated norb=2 case.

    ``cc``/``cs`` are built so that, writing ``n = p*8 + q*4 + r*2 + u``
    for the local (norb=2)^4 index tuple (p, q, r, u):

        cc[n] = 2*(n+1) + (200+2n)j        cs[n] = (100+2n) + (-2n-2)j

    so the reconstruction collapses to closed forms computed BY HAND
    (not by calling the comparator under test):

        same[n] = 0.5*(cc[n]+cs[n]) = (2n+51) + 99j
        diff[n] = 0.5*(cc[n]-cs[n]) = -49 + (2n+101)j
    """

    def _cc_cs(self):
        cc_flat = [complex(2 * n + 2, 200 + 2 * n) for n in range(16)]
        cs_flat = [complex(100 + 2 * n, -2 * n - 2) for n in range(16)]
        cc = np.array(cc_flat, dtype=complex).reshape(1, 1, 2, 2, 2, 2)
        cs = np.array(cs_flat, dtype=complex).reshape(1, 1, 2, 2, 2, 2)
        return cc, cs

    def _recon_expected(self):
        same_flat = [complex(2 * n + 51, 99) for n in range(16)]
        diff_flat = [complex(-49, 2 * n + 101) for n in range(16)]
        same = np.array(same_flat, dtype=complex).reshape(2, 2, 2, 2)
        diff = np.array(diff_flat, dtype=complex).reshape(2, 2, 2, 2)
        recon = np.zeros((1, 1, 4, 4, 4, 4), dtype=complex)
        recon[0, 0, 0:2, 0:2, 0:2, 0:2] = same
        recon[0, 0, 0:2, 0:2, 2:4, 2:4] = diff
        recon[0, 0, 2:4, 2:4, 0:2, 0:2] = diff
        recon[0, 0, 2:4, 2:4, 2:4, 2:4] = same
        return recon

    def test_exact_match_hand_enumerated(self):
        cc, cs = self._cc_cs()
        recon_expected = self._recon_expected()
        rpa_bundle = OutputBundle(chi0q=np.zeros((1, 1)), chiq=recon_expected.copy(), chiq_s=None, chiq_c=None)
        flex_bundle = OutputBundle(chi0q=np.zeros((1, 1)), chiq=None, chiq_s=cs, chiq_c=cc)

        comparator = COMPARATORS["general_from_flex_channels"]
        a, b = comparator.map("chiq", rpa_bundle, flex_bundle)
        np.testing.assert_array_equal(a, recon_expected)
        np.testing.assert_array_equal(b, recon_expected)
        self.assertEqual(comparator.residual(a, b), 0.0)
        comparator.assert_within(a, b, atol=0.0)  # must not raise

    def test_known_perturbation_is_caught_with_value_and_index(self):
        cc, cs = self._cc_cs()
        recon_expected = self._recon_expected()
        perturbed = recon_expected.copy()
        # index (0, 0, 1, 1, 3, 3): s1=0 (local p=q=1), s2=1 (local r=u=1)
        # -> n = 1*8+1*4+1*2+1 = 15 -> the "diff" block, value -49+131j.
        self.assertEqual(perturbed[0, 0, 1, 1, 3, 3], complex(-49, 131))
        perturbed[0, 0, 1, 1, 3, 3] = perturbed[0, 0, 1, 1, 3, 3] + 0.5

        rpa_bundle = OutputBundle(chi0q=np.zeros((1, 1)), chiq=perturbed, chiq_s=None, chiq_c=None)
        flex_bundle = OutputBundle(chi0q=np.zeros((1, 1)), chiq=None, chiq_s=cs, chiq_c=cc)
        comparator = COMPARATORS["general_from_flex_channels"]
        a, b = comparator.map("chiq", rpa_bundle, flex_bundle)
        self.assertAlmostEqual(comparator.residual(a, b), 0.5)
        with self.assertRaises(AssertionError) as ctx:
            comparator.assert_within(a, b, atol=0.49)
        message = str(ctx.exception)
        self.assertIn("0.5", message)
        self.assertIn("(0, 0, 1, 1, 3, 3)", message)
        comparator.assert_within(a, b, atol=0.5)  # equality passes (<=)

    def test_rejects_an_observable_other_than_chiq(self):
        cc, cs = self._cc_cs()
        rpa_bundle = OutputBundle(chi0q=np.zeros((1, 1)), chiq=self._recon_expected(), chiq_s=None, chiq_c=None)
        flex_bundle = OutputBundle(chi0q=np.zeros((1, 1)), chiq=None, chiq_s=cs, chiq_c=cc)
        with self.assertRaises(ValueError):
            COMPARATORS["general_from_flex_channels"].map("chi0q", rpa_bundle, flex_bundle)

    def test_rejects_unset_flex_channels(self):
        rpa_bundle = OutputBundle(chi0q=np.zeros((1, 1)), chiq=self._recon_expected(), chiq_s=None, chiq_c=None)
        flex_bundle = OutputBundle(chi0q=np.zeros((1, 1)), chiq=None, chiq_s=None, chiq_c=None)
        with self.assertRaises(ValueError):
            COMPARATORS["general_from_flex_channels"].map("chiq", rpa_bundle, flex_bundle)


class TestReducedBlocksComparator(unittest.TestCase):
    """``reduced_blocks`` -- the uu +- ud mapping from
    ``TestReducedOneShot`` (``tests/test_rpa_flex_oneshot_equivalence.
    py:109-134``), on a hand-enumerated norb=2 case.

    Writing ``n = i*2 + j`` for the local (norb=2)^2 index pair:

        uu[n] = (10+n) + (200+n)j          ud[n] = (1+n) + -(n+1)j

    so, computed BY HAND (not by calling the comparator under test):

        cs[n] = uu[n] - ud[n] = 9 + (201+2n)j
        cc[n] = uu[n] + ud[n] = (11+2n) + 199j
    """

    def _cq_cs_cc(self):
        uu_flat = [complex(10 + n, 200 + n) for n in range(4)]
        ud_flat = [complex(1 + n, -(n + 1)) for n in range(4)]
        cs_flat = [complex(9, 201 + 2 * n) for n in range(4)]
        cc_flat = [complex(11 + 2 * n, 199) for n in range(4)]

        uu = np.array(uu_flat, dtype=complex).reshape(2, 2)
        ud = np.array(ud_flat, dtype=complex).reshape(2, 2)
        cs = np.array(cs_flat, dtype=complex).reshape(1, 1, 2, 2)
        cc = np.array(cc_flat, dtype=complex).reshape(1, 1, 2, 2)

        cq = np.zeros((1, 1, 4, 4), dtype=complex)
        cq[0, 0, 0:2, 0:2] = uu
        cq[0, 0, 0:2, 2:4] = ud
        # rows 2:4 are never read by this mapper -- left zero deliberately.
        return cq, cs, cc, uu, ud

    def test_exact_match_hand_enumerated(self):
        cq, cs, cc, uu, ud = self._cq_cs_cc()
        rpa_bundle = OutputBundle(chi0q=np.zeros((1, 1)), chiq=cq, chiq_s=None, chiq_c=None)
        flex_bundle = OutputBundle(chi0q=np.zeros((1, 1)), chiq=None, chiq_s=cs, chiq_c=cc)

        comparator = COMPARATORS["reduced_blocks"]
        a, b = comparator.map("chiq", rpa_bundle, flex_bundle)

        expected_a = np.stack([(uu - ud).reshape(1, 1, 2, 2), (uu + ud).reshape(1, 1, 2, 2)])
        expected_b = np.stack([cs, cc])
        np.testing.assert_array_equal(a, expected_a)
        np.testing.assert_array_equal(b, expected_b)
        self.assertEqual(comparator.residual(a, b), 0.0)
        comparator.assert_within(a, b, atol=0.0)  # must not raise

    def test_known_perturbation_is_caught(self):
        cq, cs, cc, uu, ud = self._cq_cs_cc()
        perturbed = cq.copy()
        # 0.25 is exactly representable in binary floating point, so the
        # resulting residual is exact (no rounding noise at the atol
        # boundary below).
        perturbed[0, 0, 1, 3] = perturbed[0, 0, 1, 3] + 0.25  # ud[1,1] (n=3)

        rpa_bundle = OutputBundle(chi0q=np.zeros((1, 1)), chiq=perturbed, chiq_s=None, chiq_c=None)
        flex_bundle = OutputBundle(chi0q=np.zeros((1, 1)), chiq=None, chiq_s=cs, chiq_c=cc)
        comparator = COMPARATORS["reduced_blocks"]
        a, b = comparator.map("chiq", rpa_bundle, flex_bundle)
        self.assertEqual(comparator.residual(a, b), 0.25)
        with self.assertRaises(AssertionError):
            comparator.assert_within(a, b, atol=0.24)
        comparator.assert_within(a, b, atol=0.25)  # equality passes (<=)

    def test_rejects_an_observable_other_than_chiq(self):
        cq, cs, cc, _uu, _ud = self._cq_cs_cc()
        rpa_bundle = OutputBundle(chi0q=np.zeros((1, 1)), chiq=cq, chiq_s=None, chiq_c=None)
        flex_bundle = OutputBundle(chi0q=np.zeros((1, 1)), chiq=None, chiq_s=cs, chiq_c=cc)
        with self.assertRaises(ValueError):
            COMPARATORS["reduced_blocks"].map("chi0q", rpa_bundle, flex_bundle)

    def test_odd_trailing_axis_is_rejected(self):
        rpa_bundle = OutputBundle(
            chi0q=np.zeros((1, 1)), chiq=np.zeros((1, 1, 3, 3), dtype=complex), chiq_s=None, chiq_c=None
        )
        flex_bundle = OutputBundle(
            chi0q=np.zeros((1, 1)),
            chiq=None,
            chiq_s=np.zeros((1, 1, 1, 1), dtype=complex),
            chiq_c=np.zeros((1, 1, 1, 1), dtype=complex),
        )
        with self.assertRaises(ValueError):
            COMPARATORS["reduced_blocks"].map("chiq", rpa_bundle, flex_bundle)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
