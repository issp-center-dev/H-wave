"""Validation tests for the RPA/FLEX equivalence-table cell registry
(``tests/equivalence_cells.py``).

This module tests the REGISTRY MACHINERY -- the typed records, the
schema validators (``validate_registry``), ``relationship()``'s total
status-pair lookup, and ``method_name`` -- separately from the cell
content in ``CELLS``/``COVERAGE_OBLIGATIONS``. Every machinery test
constructs its own minimal ``Cell``/``SolverProof`` fixtures via the
small builders below.
"""

from __future__ import annotations

import logging
import os
import re
import tempfile
import time
import unittest
from contextlib import ExitStack
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

logger = logging.getLogger("qlms").getChild("test_equivalence_table")


# ---------------------------------------------------------------------------
# extract_bundle -- TEST-module side (it touches solver objects): builds
# an ``OutputBundle`` from one solver's post-solve ``green_info``. The
# registry module (``tests/equivalence_cells.py``) never imports solver
# code; this function is the seam where a concrete ``solver_obj``/
# ``green_info`` pair (from ``build_solver`` below, or -- for the
# machinery tests -- a small stand-in) becomes the comparator-facing
# bundle type.
# ---------------------------------------------------------------------------


def _reduce_flex_chi0q_for_reduced_scheme(solver_obj, chi0q_full: np.ndarray) -> np.ndarray:
    """FLEX FINDING: under ``calc_scheme='reduced'`` FLEX writes
    the FULL spin-block-INFLATED bare bubble to ``green_info["chi0q"]``
    -- shape ``(nmat, nvol, 2*norb, 2*norb)``, matching the
    ``chiq_s``/``chiq_c`` convention -- while RPA's reduced ``chi0q``
    stays in its already-reduced, density-diagonal shape: spin-free/
    spinful ``(nmat, nvol, norb, norb)`` (4-dim), or spin-diag
    ``(2, nmat, nvol, norb, norb)`` (5-dim, one block per spin channel).

    This had never been caught: the pre-existing oneshot suite
    (``tests/test_rpa_flex_oneshot_equivalence.py``'s
    ``TestReducedOneShot``) only ever compares ``chiq`` under
    ``calc_scheme='reduced'``, never ``chi0q`` -- and the required-
    observables freeze (Global Constraints) now mandates
    ``("chi0q", "chiq")`` for EVERY comparison cell.

    Measured on the macOS arm64 development machine (E1 norb=1 / E2
    norb=2 spin-free / E4 spin-diag): FLEX's up-up (and, for spin-diag,
    down-down) ``norb x norb`` diagonal block of its inflated ``chi0q``
    equals RPA's
    reduced ``chi0q`` to round-off (~1e-16 to ~5e-15), and the
    off-diagonal (spin-mixing) blocks are exactly zero at the
    bare-bubble level -- this is a REPRESENTATION-CONVENTION
    difference, not a numerical divergence. This mirrors EXACTLY the
    block extraction the ``reduced_blocks`` comparator already applies
    to ``chiq_s``/``chiq_c`` (``TestReducedOneShot._blocks``); ``chi0q``
    needed the identical treatment. Resolving it here (at the
    bundle-extraction seam, which is test-module-owned and NOT part of
    the closed ``COMPARATORS`` registry) keeps the ``identity``
    comparator's semantics -- and the registry's "chi0q -> identity"
    contract -- intact: by the time ``identity``
    runs, both bundle arrays already denote the SAME physical quantity.

    Only called when ``solver_obj.calc_scheme == "reduced"``.
    """

    norb = solver_obj.norb
    if solver_obj.spin_mode == "spin-diag":
        up = chi0q_full[:, :, :norb, :norb]
        down = chi0q_full[:, :, norb:, norb:]
        return np.stack([up, down])
    return chi0q_full[:, :, :norb, :norb]


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
      When ``solver_obj.calc_scheme == "reduced"``, the extracted
      ``chi0q`` is additionally block-reduced to RPA's shape -- see
      ``_reduce_flex_chi0q_for_reduced_scheme`` (a Task-5 finding).

    ``solver_obj`` is inspected for its ``calc_scheme``/``norb``/
    ``spin_mode`` attributes on the FLEX reduced-scheme path above;
    otherwise every value this function returns comes from
    ``green_info``.
    """

    if solver_kind == "rpa":
        return OutputBundle(
            chi0q=np.asarray(green_info["chi0q"]),
            chiq=np.asarray(green_info["chiq"]),
            chiq_s=None,
            chiq_c=None,
        )
    if solver_kind == "flex":
        chi0q = np.asarray(green_info["chi0q"])
        if getattr(solver_obj, "calc_scheme", None) == "reduced":
            chi0q = _reduce_flex_chi0q_for_reduced_scheme(solver_obj, chi0q)
        return OutputBundle(
            chi0q=chi0q,
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


def _schema_errors(errors):
    """Filter out "coverage obligation ..." messages from a
    ``validate_registry`` result.

    ``COVERAGE_OBLIGATIONS`` is evaluated against
    WHATEVER ``cells`` sequence is passed in -- including the small,
    deliberately-synthetic single/few-cell lists the ``TestRegistrySchema``
    tests below build to isolate ONE structural rule at a time. Those
    synthetic lists have no reason to satisfy registry-WIDE coverage
    predicates (e.g. "every interaction type appears in G1"), so tests
    that assert "this one schema rule produces zero errors" filter the
    unrelated coverage-obligation messages out first. Full-registry
    coverage is separately pinned by ``test_full_registry_satisfies_
    every_coverage_obligation`` against the real ``CELLS``.
    """

    return [e for e in errors if not e.startswith("coverage obligation ")]


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

    def test_cells_holds_the_full_37_cell_inventory(self):
        # CELLS carries rows 1-36 plus the G6 conditioning row, cell
        # 38 -- there is no row 37, so the registry holds 37 cells
        # total.
        self.assertEqual(len(CELLS), 37)
        ids = [cell.cell_id for cell in CELLS]
        self.assertEqual(len(ids), len(set(ids)), "duplicate cell_id in CELLS")
        for bootstrap_id in (
            "general.ring.onsite_coulombintra.fixedmu",
            "reduced.ring.onsite_exchange.reject",
            "so.general.construction.reject",
            "so.reduced.construction.reject",
            "general.ring.offsite_coulombinter.conditioning.mu",
        ):
            self.assertIn(bootstrap_id, ids)

    def test_coverage_obligations_populated_including_the_conditioning_predicate(self):
        self.assertEqual(
            sorted(COVERAGE_OBLIGATIONS.keys()),
            sorted(
                [
                    "g1_every_interaction_type_appears",
                    "g2_both_reduced_spin_modes_appear",
                    "every_resolver_outcome_appears",
                    "both_so_guard_sites_appear",
                    "full_kanamori_row_exists",
                    "conditioning_row_exists",
                    "at_least_one_both_reject",
                    "at_least_one_flex_reject_rpa_supported",
                    "at_least_one_rpa_only",
                ]
            ),
        )

    def test_bootstrap_cells_pass_validate_registry(self):
        errors = validate_registry(CELLS)
        self.assertEqual(errors, [])

    def test_full_registry_satisfies_every_coverage_obligation(self):
        for name, predicate in COVERAGE_OBLIGATIONS.items():
            with self.subTest(obligation=name):
                self.assertTrue(predicate(CELLS), "coverage obligation {!r} not satisfied by CELLS".format(name))

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
        self.assertEqual(_schema_errors(errors), [])

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
        self.assertEqual(_schema_errors(errors), [])

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
        self.assertEqual(_schema_errors(errors), [])

    # -- Equiv completeness / ceilings ---------------------------------------

    def test_equiv_observable_keys_must_match_required_observables(self):
        cell = _equiv_cell(required_observables=("chi0q",))
        errors = validate_registry([cell])
        self.assertTrue(any("Equiv.observables keys" in e for e in errors), errors)

    def test_equiv_atol_may_equal_the_ceiling(self):
        cell = _equiv_cell(mu_mode="fixed")
        errors = validate_registry([cell])
        self.assertEqual(_schema_errors(errors), [])

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
        self.assertEqual(_schema_errors(errors), [])

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
        self.assertEqual(_schema_errors(errors), [])

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
        self.assertEqual(_schema_errors(errors), [])


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


class _FakeReducedFlexSolverStub:
    """Minimal solver stand-in exposing exactly the three attributes
    ``_reduce_flex_chi0q_for_reduced_scheme`` reads (``calc_scheme``,
    ``norb``, ``spin_mode``) -- deliberately NOT the run_cell-oriented
    ``_FakeSolver`` below (that one models ``solve()``/``save_results()``/
    ``read_init()`` for the executor solve-count pins; this is a bare
    attribute stub for an isolated function call).
    """

    def __init__(self, norb, spin_mode):
        self.calc_scheme = "reduced"
        self.norb = norb
        self.spin_mode = spin_mode


class TestReduceFlexChi0qForReducedScheme(unittest.TestCase):
    """``_reduce_flex_chi0q_for_reduced_scheme`` -- ISOLATED unit
    coverage (a Task-5 review should-fix). ``TestExtractBundle``'s FLEX
    case passes a bare ``object()`` as ``solver_obj``, so
    ``getattr(solver_obj, "calc_scheme", None)`` is ``None`` there and
    this function is never reached from that test -- prior coverage of
    the reduction itself was end-to-end only, through the real E1/E2/E4
    fixtures (``TestBootstrapCellsRun``/the generated ``TestEquivalenceCells``
    methods). These tests call the function DIRECTLY against
    ``_FakeReducedFlexSolverStub`` with a SMALL, hand-enumerated,
    asymmetric complex ``norb=2`` array (every entry a distinct literal
    -- never generated by the mapping under test, and never sliced with
    ``[:, :, :norb, :norb]``-shaped code in this test file, so a
    transcription bug in the function under test cannot also be baked
    into the expected value).

    The shared input, ``nmat=1, nvol=1, nd=4`` (``norb=2``), written out
    row by row::

        row0 (a=0): 1+2j,  3-1j,  5+0j, -2+4j
        row1 (a=1): 0-3j,  7+7j, -4+1j,  2+2j
        row2 (a=2): 6-5j,  1+1j, -3+8j,  9-9j
        row3 (a=3): 4+4j, -1-1j,  8+3j, -6+6j

    so the two norb=2 diagonal blocks, read off directly (not sliced by
    code -- just the rows/columns 0-1 and 2-3 above, transcribed by
    hand):

        up-up   (rows/cols 0-1): [[1+2j, 3-1j], [0-3j, 7+7j]]
        down-down (rows/cols 2-3): [[-3+8j, 9-9j], [8+3j, -6+6j]]
    """

    def _hand_array(self):
        return np.array(
            [
                [
                    [1 + 2j, 3 - 1j, 5 + 0j, -2 + 4j],
                    [0 - 3j, 7 + 7j, -4 + 1j, 2 + 2j],
                    [6 - 5j, 1 + 1j, -3 + 8j, 9 - 9j],
                    [4 + 4j, -1 - 1j, 8 + 3j, -6 + 6j],
                ]
            ]
        ).reshape(1, 1, 4, 4)

    def _expected_up_up(self):
        return np.array([[1 + 2j, 3 - 1j], [0 - 3j, 7 + 7j]]).reshape(1, 1, 2, 2)

    def _expected_down_down(self):
        return np.array([[-3 + 8j, 9 - 9j], [8 + 3j, -6 + 6j]]).reshape(1, 1, 2, 2)

    def test_spin_free_returns_the_up_up_block(self):
        chi0q_full = self._hand_array()
        solver = _FakeReducedFlexSolverStub(norb=2, spin_mode="spin-free")
        result = _reduce_flex_chi0q_for_reduced_scheme(solver, chi0q_full)
        self.assertEqual(result.shape, (1, 1, 2, 2))
        np.testing.assert_array_equal(result, self._expected_up_up())

    def test_spinful_also_returns_the_up_up_block(self):
        # The function's branch is "spin-diag vs everything else" --
        # spin-free and spinful share the same up-up-only extraction
        # (reduced+spinful has no registry row, but the function itself
        # does not special-case it; pin the actual branch condition).
        chi0q_full = self._hand_array()
        solver = _FakeReducedFlexSolverStub(norb=2, spin_mode="spinful")
        result = _reduce_flex_chi0q_for_reduced_scheme(solver, chi0q_full)
        np.testing.assert_array_equal(result, self._expected_up_up())

    def test_spin_diag_stacks_up_up_and_down_down(self):
        chi0q_full = self._hand_array()
        solver = _FakeReducedFlexSolverStub(norb=2, spin_mode="spin-diag")
        result = _reduce_flex_chi0q_for_reduced_scheme(solver, chi0q_full)
        self.assertEqual(result.shape, (2, 1, 1, 2, 2))
        np.testing.assert_array_equal(result[0], self._expected_up_up())
        np.testing.assert_array_equal(result[1], self._expected_down_down())

    def test_spin_free_ignores_the_down_down_block_content(self):
        # Discrimination check: mutate the down-down block so it
        # DIFFERS from up-up (and from its original value), then assert
        # the spin-free result is still exactly the (unmutated) up-up
        # block -- not the mutated down-down block, and not some blend
        # of the two. A function that accidentally read the wrong
        # sub-block, or averaged the two, would fail this.
        chi0q_full = self._hand_array()
        chi0q_full[0, 0, 2:, 2:] = np.array([[123 - 45j, 67 + 8j], [-9 + 10j, 11 - 12j]])
        solver = _FakeReducedFlexSolverStub(norb=2, spin_mode="spin-free")
        result = _reduce_flex_chi0q_for_reduced_scheme(solver, chi0q_full)
        np.testing.assert_array_equal(result, self._expected_up_up())
        self.assertFalse(np.array_equal(result, chi0q_full[:, :, 2:, 2:]))


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
    """``COMPARATORS`` has exactly the three keys the cell inventory
    references -- no ``"scalar"`` key.
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


# ---------------------------------------------------------------------------
# Proof executors -- build_solver / run_cell + the generated per-cell
# TestCase methods.
#
# ``build_solver`` mirrors ``tests/test_rpa_flex_oneshot_equivalence.py``'s
# ``_read``/``_run`` construction verbatim (RPA's
# ``{'mode': 'RPA', 'param': par, 'calc_scheme': ..., 'calc_type': ...}``;
# FLEX's ``{'mode': 'FLEX', 'param': par + {'IterationMax': 1, 'Mix': 1.0,
# 'EPS': 1}, ...}``, no ``calc_type`` key) but with a FRESH
# ``tempfile.TemporaryDirectory()`` owned by the caller's ``ExitStack`` --
# ``tests/rpa/output`` is never touched by these cells.
# ---------------------------------------------------------------------------

# The ExecuteReject.exc_type CLOSED-set -> exception class mapping (fixed,
# per the Task-1 interface contract).
_EXC_TYPE_MAP = {"ValueError": ValueError, "RuntimeError": RuntimeError}


def _read_fixture(spec: FixtureSpec):
    """Read one fixture via ``QLMSkInput``, mirroring the oneshot
    suite's ``_read`` -- ``spec.geom``/``spec.transfer`` override the
    default ``geom.dat``/``transfer.dat`` names (the SO family);
    ``spec.extern`` names a committed Zeeman/Extern file when set (the
    public spin-diag route).
    """

    import hwave.qlmsio.read_input_k as read_input_k

    idict = {
        "path_to_input": spec.input_dir,
        "Geometry": spec.geom,
        "Transfer": spec.transfer,
    }
    idict.update(dict(spec.interactions))
    if spec.extern is not None:
        idict["Extern"] = spec.extern
    return read_input_k.QLMSkInput(
        {"path_to_input": spec.input_dir, "interaction": idict}
    )


def build_solver(spec: FixtureSpec, solver_kind: str, stack: ExitStack):
    """Construct one solver instance for a fixture. Returns
    ``(solver_obj, green_info, out_dir)``. ``out_dir`` is a FRESH
    ``tempfile.TemporaryDirectory()`` registered on ``stack`` (cleaned
    up when the caller's ``ExitStack`` closes) -- ``tests/rpa/output``
    is never used by cell fixtures.

    Mirrors ``tests/test_rpa_flex_oneshot_equivalence.py::_run``'s
    solver construction, with ONE Task-5 addition: FLEX ALSO gets
    ``calc_type`` (the oneshot suite never varies it off the "ring"
    default, so it never needed to). FLEX has no ``calc_type``
    parameter of its own, but it inherits ``_set_scheme`` from RPA
    (``FLEX(RPA)``), which reads ``info_mode.get("calc_type", "ring")``
    identically for both solvers -- omitting the key here left FLEX's
    ``self.calc_type`` silently defaulted to ``"ring"`` regardless of
    ``spec.calc_type``, which would have made cell 34's
    ``calc_type='ring+ladder'`` CONSTRUCTOR-time rejection
    unreachable through this executor. Passing ``"ring"`` explicitly
    for every OTHER cell is a no-op (same default value, verified no
    behavior change). FLEX additionally sets
    ``IterationMax=1, Mix=1.0, EPS=1`` (the one-shot / bare-bubble
    limit the whole equivalence table measures against).
    """

    import hwave.solver.flex as flex_mod
    import hwave.solver.rpa as rpa_mod

    reader = _read_fixture(spec)
    ham = reader.get_param("ham")

    par = dict(spec.extra_params)
    par["T"] = spec.T
    if spec.mu is not None:
        par["mu"] = spec.mu
    if spec.filling is not None:
        par["filling"] = spec.filling
    par["CellShape"] = list(spec.CellShape)
    par["SubShape"] = list(spec.SubShape)
    par["Nmat"] = spec.Nmat

    out_dir = stack.enter_context(tempfile.TemporaryDirectory())

    if solver_kind == "rpa":
        info = {
            "mode": "RPA",
            "param": par,
            "enable_spin_orbital": spec.enable_spin_orbital,
            "calc_scheme": spec.requested_scheme,
            "calc_type": spec.calc_type,
        }
        solver_obj = rpa_mod.RPA(ham, {}, info)
    elif solver_kind == "flex":
        pf = dict(par)
        pf.update({"IterationMax": 1, "Mix": 1.0, "EPS": 1})
        info = {
            "mode": "FLEX",
            "param": pf,
            "enable_spin_orbital": spec.enable_spin_orbital,
            "calc_scheme": spec.requested_scheme,
            "calc_type": spec.calc_type,
        }
        solver_obj = flex_mod.FLEX(ham, {}, info)
    else:
        raise ValueError(
            "build_solver: solver_kind must be 'rpa' or 'flex', got "
            "{!r}".format(solver_kind)
        )

    green_info = reader.get_param("green")
    return solver_obj, green_info, out_dir


def _construct_and_maybe_solve(fixture, solver_kind, step, stack, solver_factory=build_solver):
    """Execute ONE run-class proof step for one solver side.

    ``ExecuteConstruct``: build only -- ZERO solves (construction-only
    rows never call ``solve()``). ``ExecuteRun``: build then solve --
    EXACTLY ONE solve (the run-once-per-cell discharge rule). Returns
    ``(solver_obj, green_info, out_dir)`` in both cases.

    ``solver_factory`` defaults to ``build_solver`` and is overridable
    so the executor solve-count contract can be pinned against
    instrumented fake solvers (``TestExecutorSolveCounts`` below)
    without touching real fixtures or solver code.
    """

    if isinstance(step, ExecuteConstruct):
        return solver_factory(fixture, solver_kind, stack)
    if isinstance(step, ExecuteRun):
        solver_obj, green_info, out_dir = solver_factory(fixture, solver_kind, stack)
        solver_obj.solve(green_info, out_dir)
        return solver_obj, green_info, out_dir
    raise ValueError(
        "_construct_and_maybe_solve: unsupported run-class step {!r} -- "
        "PairedInvarianceRun/ExecuteChiqInitReuse are dedicated "
        "two-solve exceptions with their own executors (cell "
        "33)".format(step)
    )


def _execute_chiq_init_reuse(fixture, solver_kind, stack, solver_factory=build_solver):
    """The RPA side of cell 33 (``chi0q_init.reuse``): TWO RPA solves --
    the documented ``ExecuteChiqInitReuse`` exception to the run-once
    rule.

    Run A (no ``chi0q_init``) captures ``chi0q_A``/``chiq_A``. Run B
    consumes ``chi0q_A`` via the PUBLIC ``chi0q_init`` mechanism: run
    A's own ``save_results()`` writes it to a file (the same public
    entry point ``RPA.save_results`` documents), and run B's
    ``read_init()`` loads it back (``tests/test_rpa_flex_oneshot_
    equivalence.py::TestSpinModeConsistency`` establishes this same
    save-then-reload pattern, though it hand-builds the npz rather than
    routing through ``save_results``).

    The ORACLE (discharged internally, like ``_assert_reject``):
    ``np.array_equal(chi0q_B, chi0q_A)`` AND
    ``np.array_equal(chiq_B, chiq_A)`` -- reuse must reproduce the same
    dressed result bitwise. Exactly 2 solves total (one per run);
    ``TestExecutorSolveCounts`` pins this against instrumented fakes.
    """

    solver_a, green_a, out_a = solver_factory(fixture, solver_kind, stack)
    solver_a.solve(green_a, out_a)
    chi0q_a = np.array(green_a["chi0q"], copy=True)
    chiq_a = np.array(green_a["chiq"], copy=True)

    solver_a.save_results({"path_to_output": out_a, "chi0q": "chi0q_reuse.npz"}, green_a)

    solver_b, green_b, out_b = solver_factory(fixture, solver_kind, stack)
    green_b.update(solver_b.read_init({"path_to_input": out_a, "chi0q_init": "chi0q_reuse.npz"}))
    solver_b.solve(green_b, out_b)
    chi0q_b = np.asarray(green_b["chi0q"])
    chiq_b = np.asarray(green_b["chiq"])

    if not np.array_equal(chi0q_b, chi0q_a):
        raise AssertionError(
            "ExecuteChiqInitReuse: chi0q_B != chi0q_A bitwise -- "
            "chi0q_init reuse did not reproduce the same bare bubble"
        )
    if not np.array_equal(chiq_b, chiq_a):
        raise AssertionError(
            "ExecuteChiqInitReuse: chiq_B != chiq_A bitwise -- "
            "chi0q_init reuse did not reproduce the same dressed chiq"
        )


def _execute_paired_invariance_run(fixture, solver_kind, stack, solver_factory=build_solver):
    """The FLEX side of cell 33 (``chi0q_init.reuse``): TWO FLEX solves
    -- the documented ``PairedInvarianceRun`` exception to the
    run-once rule, proving the NOT_APPLICABLE claim rather than merely
    asserting it.

    FLEX starts every SCF loop from zero self-energy and recomputes
    ``chi0q`` from the dressed Green's function each iteration
    (``src/hwave/solver/flex.py:446-451``'s own docstring), so a
    ``chi0q_init`` entry loaded by the inherited RPA ``read_init`` is
    NEVER consumed. This executor drives the SAME two-solve
    construction ``_execute_chiq_init_reuse`` uses (run A plain, run B
    with a ``chi0q_init`` entry injected via the identical public
    ``save_results``/``read_init`` route) and asserts EXHAUSTIVE
    invariance: the ``green_info`` key sets match exactly, and every
    key's array is bitwise ``np.array_equal`` between the two runs --
    the option's presence changes NOTHING about FLEX's output.
    """

    solver_a, green_a, out_a = solver_factory(fixture, solver_kind, stack)
    solver_a.solve(green_a, out_a)

    solver_a.save_results({"path_to_output": out_a, "chi0q": "chi0q_reuse.npz"}, green_a)

    solver_b, green_b, out_b = solver_factory(fixture, solver_kind, stack)
    green_b.update(solver_b.read_init({"path_to_input": out_a, "chi0q_init": "chi0q_reuse.npz"}))
    solver_b.solve(green_b, out_b)

    keys_a = set(green_a.keys())
    keys_b = set(green_b.keys())
    if keys_a != keys_b:
        raise AssertionError(
            "PairedInvarianceRun: green_info key sets differ with vs "
            "without chi0q_init -- A={!r} B={!r}".format(sorted(keys_a), sorted(keys_b))
        )
    for key in keys_a:
        a_val = np.asarray(green_a[key])
        b_val = np.asarray(green_b[key])
        if not np.array_equal(a_val, b_val):
            raise AssertionError(
                "PairedInvarianceRun: green_info[{!r}] differs with vs "
                "without chi0q_init (FLEX is claimed invariant to this "
                "option)".format(key)
            )


def _assert_reject(fixture, solver_kind, step, stack, solver_factory=build_solver):
    """``ExecuteReject`` at the recorded ``site``: ``assertRaises`` the
    mapped exception class, then asserts ``step.fragment`` is a
    substring of the raised message.

    CONSTRUCTOR site: the solver never comes into existence -- ZERO
    solves. SOLVE site: construction succeeds, then ``solve()`` raises
    -- EXACTLY ONE ATTEMPTED solve.
    """

    exc_class = _EXC_TYPE_MAP[step.exc_type]
    tc = unittest.TestCase()
    if step.site is Site.CONSTRUCTOR:
        with tc.assertRaises(exc_class) as ctx:
            solver_factory(fixture, solver_kind, stack)
    elif step.site is Site.SOLVE:
        solver_obj, green_info, out_dir = solver_factory(fixture, solver_kind, stack)
        with tc.assertRaises(exc_class) as ctx:
            solver_obj.solve(green_info, out_dir)
    else:
        raise ValueError("_assert_reject: unknown Site {!r}".format(step.site))

    message = str(ctx.exception)
    if step.fragment not in message:
        raise AssertionError(
            "{} {}: ExecuteReject.fragment {!r} not found in the raised "
            "{} message {!r}".format(
                solver_kind, step.site, step.fragment, step.exc_type, message
            )
        )


def _run_side(fixture, solver_kind, proof, stack, solver_factory=build_solver):
    """Execute ONE solver side of a cell per its ``SolverProof``.

    Returns ``(solver_obj, green_info, out_dir)`` for an ``ExecuteRun``/
    ``ExecuteConstruct`` SUPPORTED proof; ``None`` for REJECT
    (discharged via ``_assert_reject``, itself), the narrowed
    zero-step NOT_APPLICABLE proof (nothing to execute -- the option
    has no corresponding semantics on this solver), and the two
    multirun exceptions -- ``ExecuteChiqInitReuse`` (SUPPORTED) and
    ``PairedInvarianceRun`` (NOT_APPLICABLE) -- which discharge their
    own bitwise oracle internally (like ``_assert_reject``) via
    ``_execute_chiq_init_reuse``/``_execute_paired_invariance_run``.
    """

    if proof.status is Status.REJECT:
        _assert_reject(fixture, solver_kind, proof.steps[0], stack, solver_factory)
        return None
    if proof.status is Status.NOT_APPLICABLE:
        if len(proof.steps) == 0:
            return None
        step = proof.steps[0]
        if isinstance(step, PairedInvarianceRun):
            _execute_paired_invariance_run(fixture, solver_kind, stack, solver_factory)
            return None
        raise NotImplementedError(  # pragma: no cover -- unreachable: validate_registry's closed step union
            "_run_side: unsupported NOT_APPLICABLE step {!r}".format(step)
        )
    # Status.SUPPORTED
    step = proof.steps[0]
    if isinstance(step, ExecuteChiqInitReuse):
        _execute_chiq_init_reuse(fixture, solver_kind, stack, solver_factory)
        return None
    return _construct_and_maybe_solve(fixture, solver_kind, step, stack, solver_factory)


def _is_construction_only(cell) -> bool:
    return (
        cell.rpa.status is Status.SUPPORTED
        and len(cell.rpa.steps) == 1
        and isinstance(cell.rpa.steps[0], ExecuteConstruct)
        and cell.flex.status is Status.SUPPORTED
        and len(cell.flex.steps) == 1
        and isinstance(cell.flex.steps[0], ExecuteConstruct)
    )


def _assert_resolved_scheme(cell, solver_kind, solver_obj) -> None:
    if solver_obj.calc_scheme != cell.resolved_scheme:
        raise AssertionError(
            "cell {!r} {}: constructed calc_scheme {!r} != "
            "Cell.resolved_scheme {!r}".format(
                cell.cell_id, solver_kind, solver_obj.calc_scheme, cell.resolved_scheme
            )
        )


def _assert_comparison(cell, rpa_bundle, flex_bundle) -> None:
    comparison = cell.comparison
    if isinstance(comparison, Equiv):
        for observable, spec in comparison.observables.items():
            comparator = COMPARATORS[spec.comparator]
            a, b = comparator.map(observable, rpa_bundle, flex_bundle)
            comparator.assert_within(a, b, spec.atol)
        return
    if isinstance(comparison, Diverges):
        for observable, spec in comparison.diverging.items():
            comparator = COMPARATORS[spec.comparator]
            a, b = comparator.map(observable, rpa_bundle, flex_bundle)
            residual = comparator.residual(a, b)
            assert_diverges_bracket(residual, spec.ceiling, spec.regression_bound)
        for observable, spec in comparison.others.items():
            comparator = COMPARATORS[spec.comparator]
            a, b = comparator.map(observable, rpa_bundle, flex_bundle)
            comparator.assert_within(a, b, spec.atol)
        return
    raise AssertionError(
        "cell {!r}: comparison must be Equiv or Diverges, got "
        "{!r}".format(cell.cell_id, comparison)
    )


def run_cell(cell, solver_factory=build_solver) -> None:
    """Execute every proof a cell records, owning its own
    ``ExitStack`` spanning construction -> solve -> comparison ->
    cleanup (Global Constraints: a fresh temp output dir and fresh
    ``green_info`` per cell run; ``tests/rpa/output`` untouched).

    Run-once discharge: for a (SUPPORTED, SUPPORTED) comparison cell,
    ``_run_side`` is called exactly once per solver side -- the single
    resulting run pair feeds every observable assertion.
    ``ExecuteReject`` is dispatched via ``assertRaises`` at the
    recorded site inside ``_run_side``/``_assert_reject``.
    ``ExecuteConstruct`` asserts the resolved ``calc_scheme`` on the
    publicly constructed solver, with no solve.

    ``solver_factory`` is a test-only override point (default
    ``build_solver``) letting ``TestExecutorSolveCounts`` pin the exact
    solve count per proof-step type against instrumented fake solvers.
    """

    with ExitStack() as stack:
        try:
            rpa_result = _run_side(cell.fixture, "rpa", cell.rpa, stack, solver_factory)
            flex_result = _run_side(cell.fixture, "flex", cell.flex, stack, solver_factory)
        except KeyboardInterrupt:
            raise
        except SystemExit as exc:
            # Several input-reader/solver construction paths (e.g.
            # src/hwave/solver/rpa.py's Lattice._init_lattice on an
            # incompatible CellShape/SubShape pairing) call
            # ``sys.exit(...)`` on malformed
            # input instead of raising. Left uncaught, that call
            # terminates the WHOLE unittest process silently (no
            # traceback, remaining cells/tests never run) -- a single
            # corrupted fixture must surface as one normal test
            # failure, not kill the runner. Converted to an
            # AssertionError with the exit code recorded; both solver-
            # execution paths above (rpa AND flex) are covered by this
            # one try/except.
            raise AssertionError(
                "cell {!r}: a solver-execution path called "
                "sys.exit({!r}) instead of raising -- this must "
                "surface as a normal test failure, not terminate the "
                "process".format(cell.cell_id, exc.code)
            ) from exc

        if cell.comparison is not None:
            rpa_obj, rpa_green, _rpa_out = rpa_result
            flex_obj, flex_green, _flex_out = flex_result
            rpa_bundle = extract_bundle(rpa_obj, rpa_green, "rpa")
            flex_bundle = extract_bundle(flex_obj, flex_green, "flex")
            _assert_comparison(cell, rpa_bundle, flex_bundle)
        elif _is_construction_only(cell):
            rpa_obj, _rpa_green, _rpa_out = rpa_result
            flex_obj, _flex_green, _flex_out = flex_result
            _assert_resolved_scheme(cell, "rpa", rpa_obj)
            _assert_resolved_scheme(cell, "flex", flex_obj)
        # Otherwise: every proof this cell records was already fully
        # discharged by _run_side (e.g. BOTH-REJECT; a REJECT side
        # paired with a SUPPORTED side that carries no comparison; or
        # cell 33's shape -- ExecuteChiqInitReuse/PairedInvarianceRun
        # each assert their own bitwise oracle internally, and their
        # (SUPPORTED, NOT_APPLICABLE) status pair never carries a
        # comparison either).


# ---------------------------------------------------------------------------
# The generated per-cell TestCase methods -- one ``test_cell__<id>``
# method per ``CELLS`` entry, attached via ``setattr`` (never a loop
# variable captured by reference -- ``_make`` binds ``cell`` in its own
# closure per call).
# ---------------------------------------------------------------------------


class TestEquivalenceCells(unittest.TestCase):
    """Every generated per-cell proof for the current ``CELLS``
    registry. Discoverable via ``unittest.TestLoader().loadTestsFromModule``
    like any other ``TestCase`` -- see
    ``TestGeneratedCellDiscovery`` below.
    """


def _make(cell):
    def test(self):
        run_cell(cell)

    test.__doc__ = "Generated proof for cell {!r} ({} <-> {}).".format(
        cell.cell_id, cell.rpa.status.value, cell.flex.status.value
    )
    return test


for _cell in CELLS:
    setattr(TestEquivalenceCells, method_name(_cell.cell_id), _make(_cell))
del _cell


class TestGeneratedCellDiscovery(unittest.TestCase):
    """The class-body generation loop actually produces discoverable
    unittest methods -- verified via ``unittest.TestLoader``, not by
    reaching into ``TestEquivalenceCells.__dict__`` directly.
    """

    def test_every_cell_has_a_discoverable_generated_method(self):
        import sys

        module = sys.modules[__name__]
        suite = unittest.TestLoader().loadTestsFromModule(module)
        discovered = set()
        stack = [suite]
        while stack:
            item = stack.pop()
            if isinstance(item, unittest.TestSuite):
                stack.extend(item)
            else:
                discovered.add("{}.{}".format(item.__class__.__name__, item._testMethodName))

        for cell in CELLS:
            expected = "TestEquivalenceCells.{}".format(method_name(cell.cell_id))
            self.assertIn(expected, discovered)

    def test_generated_method_count_matches_cells(self):
        generated = [
            name
            for name in vars(TestEquivalenceCells)
            if name.startswith("test_cell__")
        ]
        self.assertEqual(len(generated), len(CELLS))
        self.assertEqual(set(generated), {method_name(cell.cell_id) for cell in CELLS})


class TestBootstrapCellsRun(unittest.TestCase):
    """A direct (non-generated) run of the two bootstrap cells, so a
    failure here points straight at ``run_cell``/``build_solver``
    rather than requiring a trip through ``setattr`` first.
    """

    def test_bootstrap_equiv_cell_via_run_cell(self):
        cell = next(c for c in CELLS if c.cell_id == "general.ring.onsite_coulombintra.fixedmu")
        run_cell(cell)  # must not raise

    def test_bootstrap_reject_cell_via_run_cell(self):
        cell = next(c for c in CELLS if c.cell_id == "reduced.ring.onsite_exchange.reject")
        run_cell(cell)  # must not raise (both sides assertRaises internally)

    def test_so_general_construction_reject_cell_via_run_cell(self):
        cell = next(c for c in CELLS if c.cell_id == "so.general.construction.reject")
        run_cell(cell)  # must not raise -- RPA runs, FLEX assertRaises internally

    def test_so_reduced_construction_reject_cell_via_run_cell(self):
        cell = next(c for c in CELLS if c.cell_id == "so.reduced.construction.reject")
        run_cell(cell)  # must not raise -- RPA runs, FLEX assertRaises internally


# ---------------------------------------------------------------------------
# TestExecutorSolveCounts -- the instrumented FAKE-SOLVER unit test:
# stub solvers proving the EXACT solve count per proof-step type,
# independent of real fixtures/solver code (Global Constraints: "an
# executor unit test with instrumented fake solvers proves the EXACT
# solve count per step type").
# ---------------------------------------------------------------------------


class _FakeSolver:
    """A minimal stand-in for RPA/FLEX: records ``solve()`` calls and
    optionally raises on the Nth one (``fail_at_solve``), so tests can
    pin EXACT solve counts without running real physics.
    ``save_results``/``read_init`` are no-ops -- the two-solve
    executors (``_execute_chiq_init_reuse``/
    ``_execute_paired_invariance_run``) call them unconditionally, so a
    fake solver must accept the calls even though it does not model
    the underlying file round-trip (the fake tests inject fixed
    ``green_info`` content via the factory instead).
    """

    def __init__(self, fail_at_solve=None):
        self.solve_calls = 0
        self.calc_scheme = "general"
        self._fail_at_solve = fail_at_solve

    def solve(self, green_info, out_dir):
        self.solve_calls += 1
        if self._fail_at_solve is not None:
            raise self._fail_at_solve

    def save_results(self, info_outputfile, green_info):
        pass

    def read_init(self, info_inputfile):
        return {}


def _fake_factory(solver, fail_at_construct=None):
    """A ``solver_factory``-shaped callable returning a pre-built fake
    solver (or raising ``fail_at_construct`` before ever returning
    one -- ZERO solves possible in that branch).
    """

    def factory(fixture, solver_kind, stack):
        if fail_at_construct is not None:
            raise fail_at_construct
        return solver, {"chi0q": None}, "/dev/null"

    return factory


def _fake_fixture():
    # solver_factory bypasses real construction entirely -- this
    # FixtureSpec's field values are never read, only its shape.
    return FixtureSpec(
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


class TestExecutorSolveCounts(unittest.TestCase):
    """Pins how many solves each proof-step shape performs: ExecuteRun
    (1 solve), ExecuteConstruct (0 solves), constructor-site
    ExecuteReject (0 solves), solve-site ExecuteReject (exactly 1
    ATTEMPTED solve). The two multirun exceptions,
    ``PairedInvarianceRun``/``ExecuteChiqInitReuse``, are each pinned at
    EXACTLY 2 solves (see the tests under "the two multirun exception
    executors" below).
    """

    def test_execute_run_solves_exactly_once(self):
        solver = _FakeSolver()
        with ExitStack() as stack:
            obj, _green, _out = _construct_and_maybe_solve(
                _fake_fixture(), "rpa", ExecuteRun(), stack, solver_factory=_fake_factory(solver)
            )
        self.assertIs(obj, solver)
        self.assertEqual(solver.solve_calls, 1)

    def test_execute_construct_solves_zero_times(self):
        solver = _FakeSolver()
        with ExitStack() as stack:
            obj, _green, _out = _construct_and_maybe_solve(
                _fake_fixture(),
                "flex",
                ExecuteConstruct(expected_resolved_scheme="general"),
                stack,
                solver_factory=_fake_factory(solver),
            )
        self.assertIs(obj, solver)
        self.assertEqual(solver.solve_calls, 0)

    def test_constructor_site_reject_solves_zero_times(self):
        solver = _FakeSolver()
        step = ExecuteReject(site=Site.CONSTRUCTOR, exc_type="ValueError", fragment="boom")
        with ExitStack() as stack:
            _assert_reject(
                _fake_fixture(),
                "rpa",
                step,
                stack,
                solver_factory=_fake_factory(solver, fail_at_construct=ValueError("boom at construction")),
            )
        self.assertEqual(solver.solve_calls, 0)

    def test_solve_site_reject_attempts_exactly_one_solve(self):
        solver = _FakeSolver(fail_at_solve=ValueError("boom at solve"))
        step = ExecuteReject(site=Site.SOLVE, exc_type="ValueError", fragment="boom")
        with ExitStack() as stack:
            _assert_reject(
                _fake_fixture(), "flex", step, stack, solver_factory=_fake_factory(solver)
            )
        self.assertEqual(solver.solve_calls, 1)

    def test_reject_fragment_mismatch_is_reported(self):
        solver = _FakeSolver(fail_at_solve=ValueError("a totally different message"))
        step = ExecuteReject(site=Site.SOLVE, exc_type="ValueError", fragment="boom")
        with self.assertRaises(AssertionError) as ctx:
            with ExitStack() as stack:
                _assert_reject(
                    _fake_fixture(), "flex", step, stack, solver_factory=_fake_factory(solver)
                )
        self.assertIn("boom", str(ctx.exception))
        self.assertEqual(solver.solve_calls, 1)

    def test_run_cell_solves_each_side_exactly_once_on_equiv(self):
        rpa_solver = _FakeSolver()
        flex_solver = _FakeSolver()
        chi0q = np.array([[1.0 + 0.5j]])

        def factory(fixture, solver_kind, stack):
            solver = rpa_solver if solver_kind == "rpa" else flex_solver
            green = {"chi0q": chi0q.copy(), "chiq": chi0q.copy()} if solver_kind == "rpa" else {
                "chi0q": chi0q.copy(),
                "chiq_s": np.zeros((1, 1), dtype=complex),
                "chiq_c": np.zeros((1, 1), dtype=complex),
            }
            return solver, green, "/dev/null"

        cell = Cell(
            cell_id="fake.equiv.cell",
            fixture=_fake_fixture(),
            resolved_scheme="general",
            expected_spin_mode="spin-free",
            rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
            flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
            comparison=Equiv(
                observables={
                    "chi0q": ObservableSpec(comparator="identity", atol=0.0, provenance="fake"),
                }
            ),
            required_observables=("chi0q",),
            interaction_class="onsite",
            notes="",
        )
        run_cell(cell, solver_factory=factory)
        self.assertEqual(rpa_solver.solve_calls, 1)
        self.assertEqual(flex_solver.solve_calls, 1)

    def test_run_cell_solves_each_side_exactly_once_on_diverges(self):
        rpa_solver = _FakeSolver()
        flex_solver = _FakeSolver()
        ceiling = POLICY_CEILINGS["chi0q_fixed"]
        rpa_chi0q = np.array([[1.0 + 0.5j]])
        flex_chi0q = rpa_chi0q + 5.0 * ceiling  # inside (ceiling, ceiling*10]

        def factory(fixture, solver_kind, stack):
            if solver_kind == "rpa":
                return rpa_solver, {"chi0q": rpa_chi0q.copy(), "chiq": rpa_chi0q.copy()}, "/dev/null"
            return (
                flex_solver,
                {
                    "chi0q": flex_chi0q.copy(),
                    "chiq_s": np.zeros((1, 1), dtype=complex),
                    "chiq_c": np.zeros((1, 1), dtype=complex),
                },
                "/dev/null",
            )

        cell = Cell(
            cell_id="fake.diverges.cell",
            fixture=_fake_fixture(),
            resolved_scheme="general",
            expected_spin_mode="spin-free",
            rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
            flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
            comparison=Diverges(
                diverging={
                    "chi0q": DivergingSpec(
                        comparator="identity",
                        ceiling=ceiling,
                        regression_bound=ceiling * 10,
                        provenance="fake",
                    )
                },
                others={},
                step5_issue="#1",
            ),
            required_observables=("chi0q",),
            interaction_class="onsite",
            notes="",
        )
        run_cell(cell, solver_factory=factory)
        self.assertEqual(rpa_solver.solve_calls, 1)
        self.assertEqual(flex_solver.solve_calls, 1)

    def test_run_cell_solves_zero_times_on_construction_only(self):
        rpa_solver = _FakeSolver()
        flex_solver = _FakeSolver()

        def factory(fixture, solver_kind, stack):
            solver = rpa_solver if solver_kind == "rpa" else flex_solver
            solver.calc_scheme = "reduced"
            return solver, {}, "/dev/null"

        cell = Cell(
            cell_id="fake.construct.cell",
            fixture=_fake_fixture(),
            resolved_scheme="reduced",
            expected_spin_mode="spin-free",
            rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteConstruct(expected_resolved_scheme="reduced"),)),
            flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteConstruct(expected_resolved_scheme="reduced"),)),
            comparison=None,
            required_observables=(),
            interaction_class="onsite",
            notes="",
        )
        run_cell(cell, solver_factory=factory)
        self.assertEqual(rpa_solver.solve_calls, 0)
        self.assertEqual(flex_solver.solve_calls, 0)

    def test_run_cell_solves_zero_times_on_both_reject_at_constructor(self):
        def factory(fixture, solver_kind, stack):
            raise ValueError("reduced: rejected at construction")

        cell = Cell(
            cell_id="fake.bothreject.cell",
            fixture=_fake_fixture(),
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
            notes="",
        )
        run_cell(cell, solver_factory=factory)  # must not raise -- both sides reject as expected

    # -- the two multirun exception executors (cell 33) -------------------

    def test_execute_chiq_init_reuse_solves_exactly_twice(self):
        solvers = []

        def factory(fixture, solver_kind, stack):
            solver = _FakeSolver()
            solvers.append(solver)
            green = {
                "chi0q": np.array([[2.0 - 1.0j]]),
                "chiq": np.array([[3.0 + 4.0j]]),
            }
            return solver, green, "/dev/null"

        with ExitStack() as stack:
            _execute_chiq_init_reuse(_fake_fixture(), "rpa", stack, solver_factory=factory)
        self.assertEqual(len(solvers), 2)
        self.assertEqual([s.solve_calls for s in solvers], [1, 1])

    def test_execute_chiq_init_reuse_detects_a_mismatch(self):
        call_n = [0]

        def factory(fixture, solver_kind, stack):
            call_n[0] += 1
            solver = _FakeSolver()
            if call_n[0] == 1:
                green = {"chi0q": np.array([[1.0]]), "chiq": np.array([[2.0]])}
            else:
                green = {"chi0q": np.array([[1.0]]), "chiq": np.array([[999.0]])}
            return solver, green, "/dev/null"

        with self.assertRaises(AssertionError) as ctx:
            with ExitStack() as stack:
                _execute_chiq_init_reuse(_fake_fixture(), "rpa", stack, solver_factory=factory)
        self.assertIn("chiq_B != chiq_A", str(ctx.exception))
        self.assertEqual(call_n[0], 2)

    def test_paired_invariance_run_solves_exactly_twice(self):
        solvers = []

        def factory(fixture, solver_kind, stack):
            solver = _FakeSolver()
            solvers.append(solver)
            green = {
                "chi0q": np.array([[1.0 + 0.5j]]),
                "chiq_s": np.array([[0.1j]]),
                "chiq_c": np.array([[0.2j]]),
            }
            return solver, green, "/dev/null"

        with ExitStack() as stack:
            _execute_paired_invariance_run(_fake_fixture(), "flex", stack, solver_factory=factory)
        self.assertEqual(len(solvers), 2)
        self.assertEqual([s.solve_calls for s in solvers], [1, 1])

    def test_paired_invariance_run_detects_a_key_set_mismatch(self):
        call_n = [0]

        def factory(fixture, solver_kind, stack):
            call_n[0] += 1
            solver = _FakeSolver()
            green = {"chi0q": np.array([[1.0]])}
            if call_n[0] == 2:
                green["chiq_s"] = np.array([[0.0]])
            return solver, green, "/dev/null"

        with self.assertRaises(AssertionError) as ctx:
            with ExitStack() as stack:
                _execute_paired_invariance_run(_fake_fixture(), "flex", stack, solver_factory=factory)
        self.assertIn("key sets differ", str(ctx.exception))

    def test_paired_invariance_run_detects_a_value_mismatch(self):
        call_n = [0]

        def factory(fixture, solver_kind, stack):
            call_n[0] += 1
            solver = _FakeSolver()
            value = 1.0 if call_n[0] == 1 else 2.0
            green = {"chi0q": np.array([[value]])}
            return solver, green, "/dev/null"

        with self.assertRaises(AssertionError) as ctx:
            with ExitStack() as stack:
                _execute_paired_invariance_run(_fake_fixture(), "flex", stack, solver_factory=factory)
        self.assertIn("differs with vs without chi0q_init", str(ctx.exception))


class TestSystemExitEscape(unittest.TestCase):
    """Several input-reader/solver construction paths call
    ``sys.exit(...)`` on malformed input instead of raising (e.g.
    ``src/hwave/solver/rpa.py``'s
    ``Lattice._init_lattice``, which ``sys.exit(1)``s when ``SubShape``
    does not evenly divide ``CellShape`` -- rpa.py:571-573). Reading
    the input files themselves never reaches a ``sys.exit`` site in
    this codebase (``qlmsio`` raises ordinary exceptions on bad file
    content); the reachable ``sys.exit`` sites all sit in solver
    *construction*, validating already-parsed parameters against each
    other. This class pins ``run_cell``'s handling of that escape route
    using a deliberately-broken FixtureSpec ("corrupted input" at the
    parameter layer) built on a REAL, otherwise-valid committed input
    directory (``tests/equivalence_input/orb1``) driven through the
    REAL ``build_solver`` -- not a mocked factory -- so the test proves
    the actual integration, not just the try/except wrapper in
    isolation. Before this fix, ``sys.exit(1)`` from
    ``Lattice._init_lattice`` would kill the whole unittest process:
    every cell/test after the broken one would silently never run,
    with no traceback pointing at the actual cause.
    """

    def _broken_subshape_fixture(self):
        # CellShape=(4,4,1) is not evenly divisible by SubShape=(3,1,1)
        # -- hits rpa.py's Lattice._init_lattice sys.exit(1) site
        # ("SubShape is not compatible with CellShape") at the very
        # start of RPA.__init__/FLEX.__init__, before any interaction
        # file is even touched. tests/equivalence_input/orb1 is a real,
        # self-contained, otherwise-valid committed fixture directory.
        return FixtureSpec(
            input_dir="tests/equivalence_input/orb1",
            interactions={"CoulombIntra": "coulombintra.dat"},
            T=2.0,
            mu=0.0,
            filling=None,
            CellShape=(4, 4, 1),
            SubShape=(3, 1, 1),
            Nmat=32,
            extra_params={},
            calc_type="ring",
            requested_scheme="general",
            enable_spin_orbital=False,
            extern=None,
        )

    def test_run_cell_converts_system_exit_to_an_assertion_error(self):
        cell = Cell(
            cell_id="fake.systemexit.cell",
            fixture=self._broken_subshape_fixture(),
            resolved_scheme="general",
            expected_spin_mode="spin-free",
            rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
            flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
            comparison=None,
            required_observables=(),
            interaction_class="onsite",
            notes="",
        )
        with self.assertRaises(AssertionError) as ctx:
            run_cell(cell)
        message = str(ctx.exception)
        self.assertIn("fake.systemexit.cell", message)
        self.assertIn("sys.exit(1)", message)

    def test_process_continues_after_a_system_exit_cell(self):
        # If SystemExit escaped run_cell uncaught, this test method
        # itself would never reach this second assertion (the whole
        # process would already be gone) -- so simply completing this
        # test after catching the first cell's escape IS the evidence
        # that the process continues.
        cell = Cell(
            cell_id="fake.systemexit.cell",
            fixture=self._broken_subshape_fixture(),
            resolved_scheme="general",
            expected_spin_mode="spin-free",
            rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
            flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
            comparison=None,
            required_observables=(),
            interaction_class="onsite",
            notes="",
        )
        with self.assertRaises(AssertionError):
            run_cell(cell)

        # A second, ordinary cell runs normally afterward -- the
        # process was never killed.
        healthy_cell = Cell(
            cell_id="fake.systemexit.followup",
            fixture=_fake_fixture(),
            resolved_scheme="general",
            expected_spin_mode="spin-free",
            rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
            flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
            comparison=None,
            required_observables=(),
            interaction_class="onsite",
            notes="",
        )
        solver = _FakeSolver()
        run_cell(healthy_cell, solver_factory=_fake_factory(solver))
        self.assertEqual(solver.solve_calls, 2)  # both rpa and flex sides


# ---------------------------------------------------------------------------
# The oneshot canary + both bootstrap cells against the REAL fixture
# directories under tests/equivalence_input (orb1, orb2, mini) -- a
# direct sanity check that those .dat files parse and solve,
# independent of the registry/executor plumbing above.
# ---------------------------------------------------------------------------


class TestFixtureDirectoriesParse(unittest.TestCase):
    """Every fixture file under ``tests/equivalence_input/`` parses via
    ``QLMSkInput`` -- a structural smoke test, not a numerical claim
    (the registry cells are what actually exercise most of these
    files)."""

    ORB1_CASES = [
        {"CoulombIntra": "coulombintra.dat"},
        {"CoulombInter": "coulombinter.dat"},
    ]
    ORB2_CASES = [
        {"CoulombIntra": "coulombintra.dat"},
        {"CoulombInter": "onsite_inter.dat"},
        {"Hund": "hund_onsite.dat"},
        {"Ising": "ising_onsite.dat"},
        {"Exchange": "exchange_onsite.dat"},
        {"PairHop": "pairhop_onsite.dat"},
        {"PairLift": "pairlift_onsite.dat"},
        {"CoulombInter": "kanamori_interorb.dat"},
        {"Hund": "kanamori_hund.dat"},
        {"CoulombInter": "offsite_coulombinter_interorb.dat"},
        {"Hund": "offsite_hund.dat"},
        {"Ising": "offsite_ising.dat"},
        {"Exchange": "offsite_exchange.dat"},
        {"PairHop": "offsite_pairhop.dat"},
        {"Coulomb": "offsite_coulombintra.dat"},
    ]
    MINI_CASES = [
        {"CoulombIntra": "coulombintra.dat"},
    ]

    def _assert_parses(self, input_dir, interactions):
        spec = FixtureSpec(
            input_dir=input_dir,
            interactions=interactions,
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
        reader = _read_fixture(spec)
        ham = reader.get_param("ham")
        self.assertIsNotNone(ham)

    def test_orb1_fixtures_parse(self):
        for interactions in self.ORB1_CASES:
            with self.subTest(interactions=interactions):
                self._assert_parses("tests/equivalence_input/orb1", interactions)

    def test_orb2_fixtures_parse(self):
        for interactions in self.ORB2_CASES:
            with self.subTest(interactions=interactions):
                self._assert_parses("tests/equivalence_input/orb2", interactions)

    def test_mini_fixtures_parse(self):
        for interactions in self.MINI_CASES:
            with self.subTest(interactions=interactions):
                self._assert_parses("tests/equivalence_input/mini", interactions)


class TestSpinFixtureDirectoryParses(unittest.TestCase):
    """``tests/equivalence_input/spin/`` -- the E4 fixture family, which
    carries two routes: (a) the spin-diag route (E2's geom/transfer
    copies + ``extern_zeeman.dat``, hz=0.1, plus CoulombIntra/
    CoulombInter copied from E2) and (b) the spinful route
    (``geom_so.dat``/``transfer_spinful.dat`` + ``coulombintra_so.dat``,
    ``enable_spin_orbital=True``). Structural smoke test only -- the
    cells that exercise these numerically are cells 35-36 (below) and
    cells 15-16 (the spin-diag Equiv comparison rows)."""

    def _assert_parses(self, spec):
        reader = _read_fixture(spec)
        self.assertIsNotNone(reader.get_param("ham"))

    def test_spin_diag_route_coulombintra_parses(self):
        self._assert_parses(FixtureSpec(
            input_dir="tests/equivalence_input/spin",
            interactions={"CoulombIntra": "coulombintra.dat"},
            T=2.0, mu=None, filling=0.5,
            CellShape=(4, 4, 1), SubShape=(1, 1, 1), Nmat=32,
            extra_params={"coeff_extern": 1.0},
            calc_type="ring", requested_scheme="reduced",
            enable_spin_orbital=False, extern="extern_zeeman.dat",
        ))

    def test_spin_diag_route_coulombinter_parses(self):
        self._assert_parses(FixtureSpec(
            input_dir="tests/equivalence_input/spin",
            interactions={"CoulombInter": "onsite_inter.dat"},
            T=2.0, mu=None, filling=0.5,
            CellShape=(4, 4, 1), SubShape=(1, 1, 1), Nmat=32,
            extra_params={"coeff_extern": 1.0},
            calc_type="ring", requested_scheme="reduced",
            enable_spin_orbital=False, extern="extern_zeeman.dat",
        ))

    def test_spinful_route_parses(self):
        self._assert_parses(FixtureSpec(
            input_dir="tests/equivalence_input/spin",
            interactions={"CoulombIntra": "coulombintra_so.dat"},
            T=2.0, mu=0.0, filling=None,
            CellShape=(4, 4, 1), SubShape=(1, 1, 1), Nmat=32,
            extra_params={},
            calc_type="ring", requested_scheme="general",
            enable_spin_orbital=True, extern=None,
            geom="geom_so.dat", transfer="transfer_spinful.dat",
        ))


class TestReducedSpinDiagAcceptance(unittest.TestCase):
    """A spin-diag reduced FLEX fixture must SOLVE SUCCESSFULLY -- the
    E4 Extern/Zeeman route (``tests/equivalence_input/spin/
    extern_zeeman.dat``, hz=0.1) that cells 15-16 depend on via a
    ``SupplementaryLink`` to this test module. The tests below verify
    that FLEX consumes Extern exactly
    like RPA and completes the (one-shot, ``IterationMax=1``) SCF loop
    without raising; had FLEX rejected or ignored Extern, or had
    reduced not supported spin-diag, those two cells could not have
    been recorded as they are.
    """

    def _spin_diag_fixture(self, interactions):
        return FixtureSpec(
            input_dir="tests/equivalence_input/spin",
            interactions=interactions,
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

    def test_reduced_flex_solves_with_extern_coulombintra(self):
        fixture = self._spin_diag_fixture({"CoulombIntra": "coulombintra.dat"})
        with ExitStack() as stack:
            solver_obj, green_info, out_dir = build_solver(fixture, "flex", stack)
            solver_obj.solve(green_info, out_dir)  # must not raise
        self.assertEqual(solver_obj.spin_mode, "spin-diag")
        self.assertIn("chi0q", green_info)
        self.assertIn("chiq_s", green_info)
        self.assertIn("chiq_c", green_info)

    def test_reduced_flex_solves_with_extern_coulombinter(self):
        fixture = self._spin_diag_fixture({"CoulombInter": "onsite_inter.dat"})
        with ExitStack() as stack:
            solver_obj, green_info, out_dir = build_solver(fixture, "flex", stack)
            solver_obj.solve(green_info, out_dir)  # must not raise
        self.assertEqual(solver_obj.spin_mode, "spin-diag")

    def test_reduced_rpa_also_solves_with_extern(self):
        # Cells 15-16 compare RPA against FLEX on this same fixture;
        # confirm RPA also accepts it (the recorded expectation for the
        # spin-diag route is that both solvers accept it).
        fixture = self._spin_diag_fixture({"CoulombIntra": "coulombintra.dat"})
        with ExitStack() as stack:
            solver_obj, green_info, out_dir = build_solver(fixture, "rpa", stack)
            solver_obj.solve(green_info, out_dir)  # must not raise
        self.assertEqual(solver_obj.spin_mode, "spin-diag")
        self.assertIn("chiq", green_info)


class TestReducedSpinfulGuard(unittest.TestCase):
    """Pins the SOLVE-time ``FLEX._check_reduced_rejects_spinful``
    guard (``src/hwave/solver/flex.py``).

    ``test_direct_semantic_state_pin*``: construct the semantic state
    directly on a bare (not-``__init__``-ed) ``FLEX`` instance --
    setting ``_flex_general``/``spin_mode`` by hand and calling the
    extracted guard method directly, bypassing every public input
    route. This is the unit-level pin for that guard.

    ``test_public_trans_mod_route_reaches_spinful_and_is_rejected``:
    an INVESTIGATION FINDING, not defense-in-depth for a state believed
    unreachable. The G5 group was recorded on the expectation that
    reaching ``spin_mode == 'spinful'`` under reduced WITHOUT
    ``enable_spin_orbital`` is not possible through public inputs (a
    spin-mixing H0 requires the SO flag, which cell 36's construction
    guard already rejects), which is why that route has no cell of its
    own. That expectation does NOT hold: the public
    ``trans_mod`` input (``[file.input] trans_mod = "..."``, read by
    ``RPA.read_init``/``FLEX.read_init`` and documented in
    ``FLEX._solve_impl``'s own docstring -- "``green_init`` and
    ``trans_mod``... ARE consumed") sets ``do_spin_orbital = True`` in
    ``RPA._calc_epsilon_k`` UNCONDITIONALLY, independent of
    ``enable_spin_orbital`` (src/hwave/solver/rpa.py:2967-2979); the
    resulting ``spin_mode`` is then determined purely from the supplied
    H0(k)'s block structure (rpa.py:2997-3013). A user-supplied
    ``trans_mod.npz`` with a genuinely spin-mixing (off-diagonal,
    asymmetric) on-site block therefore reaches
    ``calc_scheme='reduced'`` + ``spin_mode == 'spinful'`` with
    ``enable_spin_orbital=False`` throughout -- never touching the
    CONSTRUCTOR-time guard cell 36 exercises. No registry cell is
    added for this route: recording one would also require
    characterizing RPA's own reduced+spinful-via-``trans_mod``
    behaviour, which is outside what this table measures, so it is left
    as an open question for a maintainer. The guard pinned here closes
    the gap regardless of how ``spin_mode`` came to be ``'spinful'``,
    as this test demonstrates.
    """

    def test_direct_semantic_state_pin_rejects_spinful(self):
        import hwave.solver.flex as flex_mod

        solver = flex_mod.FLEX.__new__(flex_mod.FLEX)
        solver._flex_general = False
        solver.spin_mode = "spinful"
        with self.assertRaises(ValueError) as ctx:
            solver._check_reduced_rejects_spinful()
        self.assertIn("spin_mode='spinful'", str(ctx.exception))

    def test_direct_semantic_state_pin_allows_spin_free_and_spin_diag(self):
        import hwave.solver.flex as flex_mod

        for mode in ("spin-free", "spin-diag"):
            with self.subTest(spin_mode=mode):
                solver = flex_mod.FLEX.__new__(flex_mod.FLEX)
                solver._flex_general = False
                solver.spin_mode = mode
                solver._check_reduced_rejects_spinful()  # must not raise

    def test_direct_semantic_state_pin_general_scheme_unaffected(self):
        # This method only guards the reduced scheme; the general
        # scheme's own (pre-existing) spin-free-only check lives
        # separately in _solve_impl.
        import hwave.solver.flex as flex_mod

        solver = flex_mod.FLEX.__new__(flex_mod.FLEX)
        solver._flex_general = True
        solver.spin_mode = "spinful"
        solver._check_reduced_rejects_spinful()  # must not raise

    def test_public_trans_mod_route_reaches_spinful_and_is_rejected(self):
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as flex_mod

        input_dir = "tests/equivalence_input/orb1"
        interaction = {
            "path_to_input": input_dir,
            "Geometry": "geom.dat",
            "Transfer": "transfer.dat",
        }
        reader = read_input_k.QLMSkInput(
            {"path_to_input": input_dir, "interaction": interaction}
        )
        ham = reader.get_param("ham")

        with ExitStack() as stack:
            tmp_dir = stack.enter_context(tempfile.TemporaryDirectory())
            tm_path = os.path.join(tmp_dir, "trans_mod.npz")
            # orb1's CellShape [4,4,1] -> cellvol=16; ns*norb_orig=2*1=2.
            # A genuinely spin-mixing (asymmetric off-diagonal) on-site
            # block -- NOT reachable by declaring enable_spin_orbital
            # alone, since this fixture never sets that flag.
            cellvol = 16
            tab = np.zeros((cellvol, 2, 2), dtype=complex)
            tab[:, 0, 0] = 0.3
            tab[:, 1, 1] = -0.2
            tab[:, 0, 1] = 0.1 + 0.05j
            tab[:, 1, 0] = 0.1 - 0.05j
            np.savez(tm_path, trans_mod=tab)

            param = {
                "T": 2.0, "mu": 0.0,
                "CellShape": [4, 4, 1], "SubShape": [1, 1, 1], "Nmat": 32,
            }
            info_mode = {
                "mode": "FLEX",
                "param": dict(param, IterationMax=1, Mix=1.0, EPS=1),
                "enable_spin_orbital": False,
                "calc_scheme": "reduced",
            }
            # Constructs cleanly: enable_spin_orbital=False, so the
            # CONSTRUCTOR-time guard (cell 36) never fires.
            solver = flex_mod.FLEX(ham, {}, info_mode)

            green_info = reader.get_param("green")
            green_info.update(
                solver.read_init({"path_to_input": input_dir, "trans_mod": tm_path})
            )
            out_dir = stack.enter_context(tempfile.TemporaryDirectory())
            with self.assertRaises(ValueError) as ctx:
                solver.solve(green_info, out_dir)

        self.assertEqual(solver.spin_mode, "spinful")
        self.assertIn("spin_mode='spinful'", str(ctx.exception))


# ---------------------------------------------------------------------------
# TestMuGreenDivergenceDiagnostic + TestConditioningAmplification:
# the mu/Green implementation-seam diagnostic -- the audit's top finding
# made permanently visible as a gating apparatus, plus cell 38's
# amplification obligation. Tracking issue: #160 (any future Diverges
# classification arising from this diagnostic would reference it; none
# does today -- both fixtures exercised here stay well inside their
# ceilings).
#
# ADAPTER CONTRACT (read verbatim before touching this section): a
# measurement apparatus using DIRECT internal-API invocation on
# CONSTRUCTED solver objects -- public solve() cannot expose RPA's
# solve-local mu (RPA never retains it as a public post-solve attribute;
# see OutputBundle's docstring in equivalence_cells.py) nor isolate the
# Sigma=0 seam FLEX's own solve() immediately moves past (FLEX re-solves
# mu against its own post-mix, generally NONZERO self-energy at the end
# of every SCF iteration -- flex.py:741-743). Both solvers are
# constructed via ``build_solver`` (NO ``solve()`` call) from the SAME
# ``FixtureSpec``, then ``_calc_epsilon_k(green_info)`` is invoked
# directly on each -- normally solve()'s own first step (RPA.solve
# rpa.py:1835, FLEX._solve_impl flex.py:466) -- to populate
# ``H0_eigenvalue``/``H0_eigenvector``/``spin_mode`` before any mu/Green
# call. Every fixture this diagnostic uses is mu-coupled
# (``fixture.filling`` set, ``calc_mu=True`` on both solver sides).
#
# Exact call signatures/pre-state (read flex.py ~1448-1584 and rpa.py's
# ``_find_mu``/``_calc_green``/``_calc_chi0q`` first; this docstring
# records the findings):
#
#   RPA._find_mu(Ncond, T) -> (dist, mu)                     (rpa.py:3065)
#     Ncond is the ONE-SPIN target: halved from ``self.Ncond`` when
#     ``spin_mode == "spin-free"`` (rpa.py:1837-1841), used AS-IS
#     otherwise. Pre-state: only needs ``self.H0_eigenvalue`` (set by
#     ``_calc_epsilon_k``) and ``self.ene_cutoff``.
#
#   FLEX._find_mu_dressed(sigma, beta, Ncond) -> mu           (flex.py:1489)
#     ``sigma`` shape (nblock, nmat, nvol, nd, nd), dtype complex128 --
#     the FULL self-energy array ``_calc_dressed_green`` consumes
#     (flex.py:1148-1149). ``Ncond`` uses the SAME one-spin-target
#     halving convention as ``RPA._find_mu`` (flex.py:483-486, its own
#     docstring cites ``RPA._find_mu`` explicitly). Pre-state: needs
#     ``self.H0_eigenvalue``/``H0_eigenvector`` (``_calc_epsilon_k``) --
#     internally calls ``_matsubara_number_operator(sigma, beta)`` once,
#     then root-finds via safeguarded Newton (NOT the scipy
#     bisect/Newton ``RPA._find_mu`` uses -- an independent numerical
#     method, the point of assertion 2 below). A ZERO ``sigma`` array
#     (Sigma=0) is FLEX's own pre-loop convention for "the first,
#     non-interacting-Sigma mu search" (flex.py:610-615's ``sigma =
#     xp.zeros(expected, dtype=np.complex128)`` before the SCF loop).
#
#   FLEX._calc_dressed_green(beta, mu, sigma) -> green         (flex.py:1139)
#     green_inv = (iwn+mu)*I - H0(k) - sigma; green = inv(green_inv) --
#     DIRECT matrix inversion, one nd x nd system per (block, freq, k).
#     Same shape/dtype ``sigma`` contract as above. At sigma=0 this is
#     algorithmically INDEPENDENT of RPA's eigenbasis builder below (no
#     shared code path) while representing the identical physical
#     quantity G0(k, iwn).
#
#   RPA._calc_green(beta, mu) -> (green, green_tail)           (rpa.py:3114)
#     Eigenbasis construction: green_mod.build_green(ew, ev, mu, beta,
#     nmat, coeff_tail, want_full=False) -- reconstructs G0 via the H0
#     spectral decomposition (V diag(1/(iwn+mu-ew)) V^H), never forming
#     G^{-1} or inverting a matrix. ``green_tail`` is an all-zero
#     shape-matched array when ``coeff_tail == 0`` (every fixture this
#     diagnostic uses sets ``extra_params={}``, so coeff_tail defaults to
#     0 -- ``green_tail``/the FLEX-side zero tail below are therefore
#     genuinely equal, not merely shape-compatible placeholders).
#
#   RPA._calc_chi0q(green_kw, green0_tail, beta) -> chi0q       (rpa.py:3158)
#     The SHARED bubble kernel (flex.py never overrides it -- FLEX
#     inherits this exact bound method). Called here on ``rpa_obj`` for
#     BOTH assertion-3 Green arrays, making the "shared kernel" claim
#     literal: one bound method, two Green inputs. ``green0_tail`` must
#     shape-match ``green_kw`` (bubble.py's own pairing guard); the
#     FLEX-side Green gets an explicit ``np.zeros_like`` tail since
#     ``_calc_dressed_green`` never returns one.
#
# NORMALIZATION: every residual here is a SCALAR (mu) or a raw elementwise
# max-|diff| (Green/chi0q) -- NOT a ``COMPARATORS`` entry (no bundle
# mapping; ``scalar_residual``/``assert_scalar_within`` for the two mu
# checkpoints, and ``_assert_max_diff_within`` below -- a local helper
# mirroring ``Comparator.assert_within``'s max-|diff|+index diagnostic --
# for the three array checkpoints).
# ---------------------------------------------------------------------------


def _assert_max_diff_within(a, b, atol, label):
    """The array-checkpoint counterpart of ``assert_scalar_within``:
    mirrors ``Comparator.assert_within``'s policy (exact shape then
    all-finite first, elementwise complex ``abs(a-b)``, the failure
    message reports both the max-|diff| VALUE and its INDEX) without
    requiring an ``OutputBundle`` -- the diagnostic compares raw
    Green/chi0q arrays straight off the solver internals, never through
    ``extract_bundle``.
    """

    a = np.asarray(a)
    b = np.asarray(b)
    if a.shape != b.shape:
        raise AssertionError(
            "{}: shape mismatch {!r} vs {!r}".format(label, a.shape, b.shape)
        )
    if not np.all(np.isfinite(a)):
        raise AssertionError("{}: the first array contains non-finite values".format(label))
    if not np.all(np.isfinite(b)):
        raise AssertionError("{}: the second array contains non-finite values".format(label))
    diff = np.abs(a - b)
    idx = tuple(int(i) for i in np.unravel_index(int(np.argmax(diff)), diff.shape))
    max_diff = float(diff[idx])
    if max_diff > atol:
        raise AssertionError(
            "{}: max |diff| = {!r} at index {!r} exceeds atol {!r} (CURRENT "
            "measured value, not a frozen number -- see "
            "tests/equivalence_calibration_log.md for the calibration "
            "history)".format(label, max_diff, idx, atol)
        )


def _diagnostic_residuals(fixture, solver_factory=build_solver):
    """Run the ORDERED five-assertion mu/Green divergence diagnostic on
    ONE fixture and return every checkpoint value + residual.

    Both solvers are constructed (never solved) from ``fixture``, then
    ``_calc_epsilon_k`` is invoked directly on each (see the module
    docstring above the call-signature contract this function
    implements). Returns a dict with the raw checkpoint values (so
    callers can build rich, dynamic failure messages via
    ``assert_scalar_within``/``_assert_max_diff_within``) AND the five
    named scalar residuals (``"assertion1"``..``"assertion5"``, all
    ``float``) for the amplification ratio (``TestConditioningAmplification``).
    """

    with ExitStack() as stack:
        rpa_obj, rpa_green, _rpa_out = solver_factory(fixture, "rpa", stack)
        flex_obj, flex_green, _flex_out = solver_factory(fixture, "flex", stack)

        # Normally solve()'s own first step (rpa.py:1835 / flex.py:466);
        # invoked directly here since neither solver is ever solved.
        rpa_obj._calc_epsilon_k(rpa_green)
        flex_obj._calc_epsilon_k(flex_green)

        beta = 1.0 / fixture.T

        # Same one-spin-target halving convention both solve() paths use
        # (rpa.py:1837-1841, flex.py:483-486).
        ncond_rpa = rpa_obj.Ncond / 2 if rpa_obj.spin_mode == "spin-free" else rpa_obj.Ncond
        ncond_flex = flex_obj.Ncond / 2 if flex_obj.spin_mode == "spin-free" else flex_obj.Ncond

        # ---- assertion 1: shared-mu identity control ----
        # RPA._find_mu invoked on the RPA instance vs the SAME method
        # (inherited unchanged) invoked on the FLEX instance -- expect
        # round-off agreement; a nonzero residual here means the two
        # constructed solver objects/pre-state differ, not the
        # algorithms (they run the identical code path).
        _dist_rpa, mu_rpa = rpa_obj._find_mu(ncond_rpa, rpa_obj.T)
        _dist_flex, mu_flex_instance = flex_obj._find_mu(ncond_flex, flex_obj.T)
        assertion1 = scalar_residual(mu_rpa, mu_flex_instance)

        # ---- assertion 2: the mu seam ----
        # RPA._find_mu's non-interacting root vs FLEX._find_mu_dressed's
        # Sigma=0 root -- two INDEPENDENT implementations of the same
        # physical mu.
        nblock, nvol, nd = flex_obj.H0_eigenvalue.shape
        sigma_zero = np.zeros((nblock, flex_obj.nmat, nvol, nd, nd), dtype=np.complex128)
        mu_flex_dressed_zero = flex_obj._find_mu_dressed(sigma_zero, beta, ncond_flex)
        assertion2 = scalar_residual(mu_rpa, mu_flex_dressed_zero)

        # ---- assertion 3: the Green-builder seam (bare), both at the
        # SAME shared mu (assertion 1's mu_rpa) ----
        green_rpa_at_shared_mu, tail_rpa_at_shared_mu = rpa_obj._calc_green(beta, mu_rpa)
        green_flex_at_shared_mu = flex_obj._calc_dressed_green(beta, mu_rpa, sigma_zero)
        assertion3 = float(np.max(np.abs(
            np.asarray(green_rpa_at_shared_mu) - np.asarray(green_flex_at_shared_mu)
        )))

        # ---- assertion 4: the composed seam -- the same two builders at
        # _find_mu_dressed's Sigma=0 mu (assertion 2's FLEX value) ----
        green_rpa_at_dressed_zero_mu, _tail_rpa_at_dressed_zero_mu = rpa_obj._calc_green(
            beta, mu_flex_dressed_zero
        )
        green_flex_at_dressed_zero_mu = flex_obj._calc_dressed_green(
            beta, mu_flex_dressed_zero, sigma_zero
        )
        assertion4 = float(np.max(np.abs(
            np.asarray(green_rpa_at_dressed_zero_mu) - np.asarray(green_flex_at_dressed_zero_mu)
        )))

        # ---- assertion 5: downstream chi0q from assertion 3's TWO
        # Greens through the shared bubble kernel (rpa_obj._calc_chi0q,
        # bound but identical on both instances -- see the module
        # docstring) ----
        chi0q_from_rpa_green = rpa_obj._calc_chi0q(
            green_rpa_at_shared_mu, tail_rpa_at_shared_mu, beta
        )
        zero_tail_for_flex_green = np.zeros_like(green_flex_at_shared_mu)
        chi0q_from_flex_green = rpa_obj._calc_chi0q(
            green_flex_at_shared_mu, zero_tail_for_flex_green, beta
        )
        assertion5 = float(np.max(np.abs(
            np.asarray(chi0q_from_rpa_green) - np.asarray(chi0q_from_flex_green)
        )))

    return {
        "mu_rpa": mu_rpa,
        "mu_flex_instance": mu_flex_instance,
        "mu_flex_dressed_zero": mu_flex_dressed_zero,
        "green_rpa_at_shared_mu": green_rpa_at_shared_mu,
        "green_flex_at_shared_mu": green_flex_at_shared_mu,
        "green_rpa_at_dressed_zero_mu": green_rpa_at_dressed_zero_mu,
        "green_flex_at_dressed_zero_mu": green_flex_at_dressed_zero_mu,
        "chi0q_from_rpa_green": chi0q_from_rpa_green,
        "chi0q_from_flex_green": chi0q_from_flex_green,
        "assertion1": assertion1,
        "assertion2": assertion2,
        "assertion3": assertion3,
        "assertion4": assertion4,
        "assertion5": assertion5,
    }


# The benign fixture is cell 8's (general.ring.onsite_u_v_hund.mu, E2,
# T=2.0); the FC conditioning fixture is cell 38's -- both are looked
# up from CELLS rather than duplicated, so they stay
# tied to the registry's actual fixture values.
_DIAGNOSTIC_BENIGN_FIXTURE = next(
    c for c in CELLS if c.cell_id == "general.ring.onsite_u_v_hund.mu"
).fixture
_DIAGNOSTIC_FC_FIXTURE = next(
    c for c in CELLS if c.cell_id == "general.ring.offsite_coulombinter.conditioning.mu"
).fixture


class TestMuGreenDivergenceDiagnostic(unittest.TestCase):
    """The five ORDERED assertions of the mu/Green divergence
    diagnostic, run on the BENIGN fixture (cell 8's) -- proving the
    diagnostic apparatus itself stays inside its ``*_diag``/``chi0q_mu``
    ceilings under ordinary conditions. The amplified (FC) case is
    ``TestConditioningAmplification`` below; see this module's own
    section docstring (above ``_assert_max_diff_within``) for the full
    adapter contract, exact call signatures, and Step-5 issue (#160).

    Five SEPARATE, independently-discoverable test methods (one per
    assertion) sharing ONE diagnostic run via ``setUpClass`` -- the
    apparatus is deterministic and side-effect-free per call, so
    computing it once and asserting five times is equivalent to five
    independent runs, without repeating the (cheap but non-trivial)
    solver construction five times.
    """

    @classmethod
    def setUpClass(cls):
        cls.values = _diagnostic_residuals(_DIAGNOSTIC_BENIGN_FIXTURE)

    def test_1_shared_mu_identity_control(self):
        v = self.values
        assert_scalar_within(v["mu_rpa"], v["mu_flex_instance"], POLICY_CEILINGS["mu_diag"])

    def test_2_the_mu_seam(self):
        v = self.values
        assert_scalar_within(
            v["mu_rpa"], v["mu_flex_dressed_zero"], POLICY_CEILINGS["mu_diag"]
        )

    def test_3_green_builder_seam_bare(self):
        v = self.values
        _assert_max_diff_within(
            v["green_rpa_at_shared_mu"], v["green_flex_at_shared_mu"],
            POLICY_CEILINGS["green_diag"],
            "assertion 3 (Green-builder seam, bare, at the shared mu)",
        )

    def test_4_composed_seam_dressed_zero_sigma(self):
        v = self.values
        _assert_max_diff_within(
            v["green_rpa_at_dressed_zero_mu"], v["green_flex_at_dressed_zero_mu"],
            POLICY_CEILINGS["green_diag"],
            "assertion 4 (composed seam, at _find_mu_dressed's Sigma=0 mu)",
        )

    def test_5_downstream_chi0q(self):
        v = self.values
        _assert_max_diff_within(
            v["chi0q_from_rpa_green"], v["chi0q_from_flex_green"],
            POLICY_CEILINGS["chi0q_mu"],
            "assertion 5 (downstream chi0q from assertion 3's two Greens)",
        )


class TestConditioningAmplification(unittest.TestCase):
    """Cell 38's amplification obligation: the SAME diagnostic
    apparatus (``_diagnostic_residuals``),
    run independently on the benign fixture (cell 8's) and the FC
    conditioning fixture (cell 38's), must show >= 10x amplification on
    AT LEAST ONE of the three seam residuals -- assertion 2 (the mu
    seam), assertion 4 (the composed seam), assertion 5 (downstream
    chi0q) -- with ``max(residual, 1e-15)`` denominators (an
    exactly-zero benign residual, as measured on assertion 2 here,
    would otherwise divide by zero).

    MONOTONICITY DERIVATION (why halving T is the calibration lever):
    reducing T sharpens the Fermi step (``RPA._find_mu``'s
    ``_fermi`` / ``FLEX._fermi_occupation``) against a near-degenerate
    H0 eigenvalue cluster. On the FC fixture's actual CellShape=(4,4,1)
    16-k-point mesh, H0's eigenvalues are
    ``[-3,-2,-2,-2,-2,-1,-1,-1,-1,1,1,2,2,2,2,5]``; filling=0.5 (one-spin
    target 8 of 16 states) pins the Fermi level EXACTLY inside the
    4-fold-degenerate eps=-1.0 cluster -- a coarse-mesh analogue of a van
    Hove shoulder (a finer 64x64 mesh on the same dispersion shows the
    true DOS peak near eps~-0.93, consistent with this cluster; see cell
    38's notes in equivalence_cells.py). dN/dmu grows as T shrinks near
    a near-degenerate manifold, so the two INDEPENDENT mu-root-finders
    (RPA's scipy bisect+Newton on the non-interacting Fermi count vs
    FLEX's safeguarded-Newton eigenvalue root on Sigma=0,
    ``_find_mu_dressed``) settle at their own numerical fixed points
    inside an increasingly narrow-but-nonzero window -- amplifying their
    disagreement as T shrinks. This predicts assertion 2's residual
    grows monotonically as T is halved; measured on the macOS arm64
    development machine (diagnostic-only -- not gated here) at
    T=0.2/0.1/0.05 on this fixture: 1.070e-12 -> 1.514e-12 ->
    1.687e-12, confirming the trend.

    T=0.2, the first temperature tried, already clears the >=10x bar on
    assertion 2 ALONE by roughly three orders of magnitude, so no
    halve-T re-choice was needed -- the single logged attempt is in
    ``tests/equivalence_calibration_log.md``.
    """

    @classmethod
    def setUpClass(cls):
        cls.benign = _diagnostic_residuals(_DIAGNOSTIC_BENIGN_FIXTURE)
        cls.fc = _diagnostic_residuals(_DIAGNOSTIC_FC_FIXTURE)

    def test_at_least_one_seam_amplifies_at_least_tenfold(self):
        benign = self.benign
        fc = self.fc
        keys = ("assertion2", "assertion4", "assertion5")
        ratios = {key: fc[key] / max(benign[key], 1e-15) for key in keys}
        self.assertTrue(
            any(ratio >= 10.0 for ratio in ratios.values()),
            "conditioning amplification criterion (cell 38, group G6) "
            "failed: none of the seam residuals amplified >= 10x "
            "from the benign fixture to FC -- CURRENT measured ratios "
            "{!r} (benign residuals {!r}, FC residuals {!r})".format(
                ratios,
                {key: benign[key] for key in keys},
                {key: fc[key] for key in keys},
            ),
        )


# ---------------------------------------------------------------------------
# TestSupplementaryLinks: import-based existence/discoverability
# of every SupplementaryLink recorded anywhere in CELLS. SupplementaryLink
# carries ZERO proof authority (Global Constraints) -- schema-only FORMAT
# validation (the "<dotted.module>::<TestCase>::<method>" shape, duplicate
# rejection) lives in equivalence_cells.validate_registry; THIS is the
# stale-pointer guard: the module actually imports, the class actually
# resolves via getattr, and the method actually exists and is callable.
# ---------------------------------------------------------------------------


class TestSupplementaryLinks(unittest.TestCase):
    def _all_links(self):
        links = []
        for cell in CELLS:
            links.extend(cell.rpa.links)
            links.extend(cell.flex.links)
        return links

    def test_every_supplementary_link_resolves(self):
        import importlib

        links = self._all_links()
        self.assertTrue(links, "expected at least one SupplementaryLink in CELLS")
        for link in links:
            with self.subTest(test_id=link.test_id):
                module_path, class_name, method = link.test_id.split("::")
                module = importlib.import_module(module_path)
                cls = getattr(module, class_name)
                self.assertTrue(
                    issubclass(cls, unittest.TestCase),
                    "{}::{} is not a unittest.TestCase subclass".format(module_path, class_name),
                )
                bound = getattr(cls, method)
                self.assertTrue(
                    callable(bound),
                    "{}::{}::{} is not callable".format(module_path, class_name, method),
                )

    def test_every_cell_carrying_links_appears_in_the_all_links_sweep(self):
        # A cheap sanity check that _all_links() actually walks CELLS
        # (not e.g. an accidentally-empty fixture) -- the real
        # existence/discoverability assertion is
        # test_every_supplementary_link_resolves above. Cell-to-cell
        # link REUSE (the same test cited as evidence for more than
        # one cell, e.g. cells 21-25 all citing
        # TestOffsiteGeneralFLEX::test_rejected_offsite_classes) is
        # legitimate -- SupplementaryLink carries zero proof authority,
        # so citing one existing test from several cells is not a
        # defect.
        cells_with_links = [
            cell for cell in CELLS if cell.rpa.links or cell.flex.links
        ]
        self.assertTrue(cells_with_links)
        self.assertGreaterEqual(len(self._all_links()), len(cells_with_links))


class TestBenchmarkRegistryTie(unittest.TestCase):
    """``tests/equivalence_benchmark.md``'s MOST RECENT section must
    cover EXACTLY the current ``CELLS`` inventory -- an EXACT multiset
    match on cell_id (duplicates detected, not masked by set
    semantics), per the registry docstring's maintenance checklist
    (``tests/equivalence_cells.py``'s module docstring) and the
    benchmark file's own header. A cell added to (or removed from) the
    registry without a matching new benchmark section fails this test.
    """

    _CELL_ID_ROW_RE = re.compile(r'^\|\s*`([^`]+)`\s*\|\s*[-\d.eE]+\s*\|\s*$')

    @classmethod
    def _most_recent_section_cell_ids(cls):
        """Parse ``tests/equivalence_benchmark.md``'s LAST
        ``cell_id | max_seconds_over_runs_and_runners``-shaped table
        (the most recently appended calibration section -- the file
        grows one section per stage, oldest first, exactly like
        ``tests/equivalence_calibration_log.md``) and return the list
        of cell_ids in that table, IN ORDER, WITH duplicates preserved.
        """

        path = os.path.join(os.path.dirname(__file__), "equivalence_benchmark.md")
        with open(path, "r", encoding="utf-8") as f:
            lines = f.readlines()

        # Every "cell_id | max_seconds_over_runs_and_runners" table in
        # the file starts with this exact header + separator pair;
        # collect the row-blocks that follow EACH occurrence, then keep
        # only the LAST block (the most recent section).
        header_indices = [
            i for i, line in enumerate(lines)
            if line.strip() == "| cell_id | max_seconds_over_runs_and_runners |"
        ]
        if not header_indices:
            raise AssertionError(
                "tests/equivalence_benchmark.md: no "
                "'cell_id | max_seconds_over_runs_and_runners' table found"
            )
        start = header_indices[-1] + 2  # skip the header line + its '|---|---|' separator
        cell_ids = []
        for line in lines[start:]:
            stripped = line.rstrip("\n")
            if not stripped.startswith("|"):
                break  # the table block ends at the first non-row line
            m = cls._CELL_ID_ROW_RE.match(stripped)
            if not m:
                break
            cell_ids.append(m.group(1))
        return cell_ids

    def test_most_recent_section_is_an_exact_multiset_match_with_cells(self):
        from collections import Counter

        benchmark_ids = self._most_recent_section_cell_ids()
        self.assertTrue(benchmark_ids, "parsed zero cell_id rows from the most recent benchmark section")

        registry_ids = [cell.cell_id for cell in CELLS]
        self.assertEqual(
            Counter(benchmark_ids), Counter(registry_ids),
            "tests/equivalence_benchmark.md's most recent section's cell_id "
            "column is not an exact multiset match with CELLS -- a cell was "
            "added/removed/renamed/duplicated without a matching benchmark "
            "update",
        )


# ---------------------------------------------------------------------------
# TestRenderedEquivalenceTableDrift: the committed user-facing
# page ``docs/en/source/algorithm/rpa_flex_equivalence.rst`` must be
# EXACTLY what ``docs/tools/render_equivalence_table.py`` renders from
# this registry. The Sphinx build never imports the registry (the RST is
# static); this in-memory drift test is what keeps the static page and
# the measured records from separating. It never writes the repository.
#
# ``docs/tools/`` is deliberately NOT an importable package (it holds a
# developer script, not shipped code), so the renderer is loaded BY FILE
# PATH via ``importlib.util``. A missing renderer or a missing RST is a
# FAILURE, never a skip -- an absent artifact is precisely the drift this
# test exists to catch.
# ---------------------------------------------------------------------------

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_RENDERER_PATH = os.path.join(
    _REPO_ROOT, "docs", "tools", "render_equivalence_table.py"
)
_RENDERED_RST_PATH = os.path.join(
    _REPO_ROOT, "docs", "en", "source", "algorithm", "rpa_flex_equivalence.rst"
)
_REGENERATE_HINT = "python docs/tools/render_equivalence_table.py --write"


def _load_renderer_module():
    """Load the generator script by file path (``docs/tools`` is not a
    package). Raises ``AssertionError`` -- a test FAILURE -- when the
    script is absent, so a deleted generator can never pass as a skip.
    """

    import importlib.util

    if not os.path.isfile(_RENDERER_PATH):
        raise AssertionError(
            "the equivalence-table generator is missing: expected "
            "docs/tools/render_equivalence_table.py"
        )
    spec = importlib.util.spec_from_file_location(
        "hwave_docs_render_equivalence_table", _RENDERER_PATH
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestRenderedEquivalenceTableDrift(unittest.TestCase):
    def _render(self) -> str:
        return _load_renderer_module().render()

    def _committed(self) -> str:
        if not os.path.isfile(_RENDERED_RST_PATH):
            raise AssertionError(
                "the committed equivalence-table page is missing: expected "
                "docs/en/source/algorithm/rpa_flex_equivalence.rst -- "
                "generate it with `{}`".format(_REGENERATE_HINT)
            )
        # newline="" disables universal-newline translation, so what is
        # read is exactly what is on disk; the CRLF fold that follows then
        # makes the comparison platform-INDEPENDENT, since a checkout with
        # core.autocrlf=true (there is no .gitattributes pinning *.rst to
        # LF) would otherwise fail this test spuriously. The write path
        # stays LF-exact, so the committed file itself is unaffected.
        with open(_RENDERED_RST_PATH, "r", encoding="utf-8", newline="") as f:
            return f.read().replace("\r\n", "\n")

    def test_committed_page_matches_the_registry_rendering(self):
        rendered = self._render()
        committed = self._committed()
        if rendered != committed:
            rendered_lines = rendered.splitlines()
            committed_lines = committed.splitlines()
            first = None
            for i, (a, b) in enumerate(zip(committed_lines, rendered_lines)):
                if a != b:
                    first = i
                    break
            if first is None:
                first = min(len(committed_lines), len(rendered_lines))
            detail = (
                "first difference at line {}:\n  committed: {!r}\n  "
                "rendered:  {!r}".format(
                    first + 1,
                    committed_lines[first] if first < len(committed_lines) else "<end of file>",
                    rendered_lines[first] if first < len(rendered_lines) else "<end of file>",
                )
            )
            self.maxDiff = 2000
            self.assertEqual(
                rendered,
                committed,
                "docs/en/source/algorithm/rpa_flex_equivalence.rst has drifted "
                "from tests/equivalence_cells.py -- regenerate it with `{}` "
                "and commit the result.\n{}".format(_REGENERATE_HINT, detail),
            )

    def test_rendering_is_deterministic(self):
        # The generator is a pure function of the registry: no
        # timestamps, no live measurements, no environment-dependent
        # content, so two renderings are byte-identical.
        module = _load_renderer_module()
        self.assertEqual(module.render(), module.render())

    _FROZEN_PROVENANCE = {
        "source_sha": "0123456789abcdef0123456789abcdef01234567",
        "run_ids": ("1234567890",),
        "status": "frozen",
    }

    def test_freezing_the_provenance_only_changes_the_status_paragraph(self):
        # PROVENANCE feeds no proof, no tolerance and no fixture, so
        # flipping it from "candidate" to "frozen" must regenerate the
        # page as a metadata-only change: the measured content is
        # untouched and only the status paragraph moves.
        module = _load_renderer_module()
        candidate = module.render()
        frozen = module.render(provenance=self._FROZEN_PROVENANCE)

        self.assertNotEqual(candidate, frozen)
        self.assertIn(self._FROZEN_PROVENANCE["source_sha"], frozen)
        self.assertIn(self._FROZEN_PROVENANCE["run_ids"][0], frozen)

        def _outside_the_status_section(text):
            before, rest = text.split("Status of these records", 1)
            _, after = rest.split("How to read the table", 1)
            return before, after

        self.assertEqual(
            _outside_the_status_section(candidate),
            _outside_the_status_section(frozen),
            "flipping PROVENANCE changed page content outside the status "
            "section -- the freeze regeneration is no longer metadata-only",
        )


# ---------------------------------------------------------------------------
# Non-gating module wall-time timer. The CI budget gate is the
# dedicated `python -m unittest tests.test_rpa_flex_equivalence_table`
# timing step in the calibration workflow, taken as the MAX over the
# four gating runners; this module-level timer is a local, non-gating
# early warning only.
# ---------------------------------------------------------------------------

_MODULE_START_TIME = None
_MODULE_TIME_BUDGET_SECONDS = 120.0


def setUpModule():
    global _MODULE_START_TIME
    _MODULE_START_TIME = time.monotonic()


def tearDownModule():
    elapsed = time.monotonic() - _MODULE_START_TIME
    logger.info("test_rpa_flex_equivalence_table wall time: %.3fs", elapsed)
    if elapsed > _MODULE_TIME_BUDGET_SECONDS:
        logger.warning(
            "*** test_rpa_flex_equivalence_table took %.3fs, over the "
            "%.0fs CI budget (NON-GATING here -- the freeze-time gate "
            "is the calibration workflow's dedicated `python -m unittest "
            "tests.test_rpa_flex_equivalence_table` timing step, MAX "
            "over the four gating runners; this local timer is an "
            "early warning only) ***",
            elapsed,
            _MODULE_TIME_BUDGET_SECONDS,
        )


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
