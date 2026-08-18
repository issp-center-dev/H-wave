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

from tests.equivalence_cells import (
    CELLS,
    COVERAGE_OBLIGATIONS,
    POLICY_CEILINGS,
    PROOF_STEP_TYPES,
    PROVENANCE,
    Cell,
    DivergingSpec,
    Diverges,
    Equiv,
    ExecuteChiqInitReuse,
    ExecuteConstruct,
    ExecuteReject,
    FixtureSpec,
    ObservableSpec,
    PairedInvarianceRun,
    ExecuteRun,
    Site,
    SolverProof,
    Status,
    SupplementaryLink,
    method_name,
    relationship,
    validate_registry,
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


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
