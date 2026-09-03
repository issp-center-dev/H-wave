"""Registry of the OPT-IN ("heavy") tests, and the decorator that gates them.

The suite runs in two shapes:

* the FAST gate -- what every push/PR runs. ``HWAVE_FULL_TESTS`` is unset,
  and every test named in :data:`HEAVY_TESTS` below is skipped.
* the FULL suite -- run with ``HWAVE_FULL_TESTS=1``. Nothing is skipped by
  this module.

    HWAVE_FULL_TESTS=1 python -m unittest discover tests/ -v
    HWAVE_FULL_TESTS=1 pytest tests/ -v

Why this exists: a measured full ``pytest tests/`` run was 882 s over 2252
tests, and 631 s of that (71%) sat in 86 tests spread over the four
exact-diagonalization (ED) oracle modules, plus a single 171 s case in
``test_bubble_kernel.py``. Those tests re-derive numerical CONTENT from a
dense Lehmann/ED sum; that content only moves when a vertex table or a
bubble kernel is edited deliberately, which is exactly the change a
pre-merge FULL run is there to catch. Everyday structural edits do not
move it, so paying 10+ minutes for them on every push buys nothing.

THE SELECTION RULE
------------------
**A test whose OWN runtime is roughly 5 s or more is opt-in.** Not the
module's runtime, not the class's -- the individual test method's. That is
why the entries below are per-method: several of these classes hold both a
30 s ED granule and a sub-second guard, and the guard stays in the fast
gate.

THE TWO DOCUMENTED EXCEPTIONS
-----------------------------
Four tests exceed the 5 s rule and STAY in the fast gate anyway. They are
NOT decorated, and they are deliberately absent from :data:`HEAVY_TESTS`:

* ``test_spinful_transverse_ed.TestGateS0ExtractionAdjudication``
  (``test_step3_dressed_extraction_matches_totalned``,
  ``test_step4_g2_g4_swap_mutation_must_fail``), ~31 s
* ``test_spinful_transverse_ed.TestGateS2ProductionRouting``
  (``test_old_slice_then_solve_fails_the_s0_granule``,
  ``test_production_passes_the_s0_granule``), ~31 s

Reason: everything in :data:`HEAVY_TESTS` pins per-type numerical CONTENT.
These four pin the production ROUTING instead -- that ``solve()`` takes the
dressed-extraction path at all, and that the old slice-then-solve path
still fails the S0 granule. Routing is what an everyday structural change
actually breaks, and these four are also what fixes issue #110 in place;
losing them from the per-push gate would be losing the guard, not deferring
a re-measurement. Their cost is accepted knowingly.

CONSISTENCY IS ENFORCED, NOT TRUSTED
------------------------------------
``tests/test_heavy_tests_registry.py`` fails the FAST gate when this file
and reality disagree: an entry naming a test that does not exist, a test
carrying :func:`heavy` but missing from :data:`HEAVY_TESTS`, or an entry
whose test is not actually skipped with ``HWAVE_FULL_TESTS`` unset.

The gate is a ``unittest``-level skip (``unittest.skipUnless``) on purpose:
the gating CI runner is ``python -m unittest discover``, where a
module-level ``pytest.skip`` degrades to a collection ERROR rather than a
skip.
"""

import os
import unittest

#: Environment variable that opts into the full suite. Any value other than
#: the falsey spellings below (and other than "unset") enables it.
FULL_SUITE_ENV = "HWAVE_FULL_TESTS"

#: Spellings of "off". Everything else -- "1", "true", "yes", "on", and any
#: other non-empty string -- is truthy. Kept deliberately small so a typo
#: like HWAVE_FULL_TESTS=ture errs towards RUNNING the tests rather than
#: silently skipping them in a run the operator believed was full.
_FALSEY = frozenset({"", "0", "false", "no", "off"})


def full_suite_enabled():
    """True when ``HWAVE_FULL_TESTS`` asks for the full suite.

    Read at import time by :func:`heavy` (that is when ``unittest`` needs
    the skip decision), and again by the consistency test, which reports
    which shape of the suite it is checking.
    """
    return os.environ.get(FULL_SUITE_ENV, "").strip().lower() not in _FALSEY


#: The opt-in set, as data: ``(module, class, method, reason)``. ``module``
#: is the bare module name under ``tests/`` (no ``tests.`` prefix and no
#: ``.py``), matching how the consistency test and the CI discovery command
#: both name modules. Keep it sorted by module then class then method.
#:
#: To add a test here: measure it, confirm it is >= ~5 s on its own, add
#: the row, and decorate the method with :func:`heavy`. Doing only one of
#: the two fails ``test_heavy_tests_registry.py``.
HEAVY_TESTS = (
    # --- tests/test_bond_transverse_ed.py -- W_pm bond ED granules -------
    ("test_bond_transverse_ed", "TestTask6ProductionPipelineGranules",
     "test_granule_multiorbital_offsite_coulombinter_regression",
     "multi-orbital off-site CoulombInter production granule vs ED"),
    ("test_bond_transverse_ed", "TestTask6ProductionPipelineGranules",
     "test_granule_multiorbital_offsite_exchange_both_orientations",
     "multi-orbital off-site Exchange production granule, both orientations"),
    ("test_bond_transverse_ed", "TestW0Granules",
     "test_granule_a_multiorbital_offsite_exchange",
     "W0 granule A: multi-orbital off-site Exchange vs ED"),
    # --- tests/test_bond_vs_ed_oracle.py -- bond-channel ED oracle -------
    ("test_bond_vs_ed_oracle", "TestCaseM", "test_joint_ray_ph_1_1",
     "case-M joint-ray ph adjudication at (1,1)"),
    ("test_bond_vs_ed_oracle", "TestCaseM", "test_joint_ray_pp_1_1",
     "case-M joint-ray pp adjudication at (1,1)"),
    ("test_bond_vs_ed_oracle", "TestPpInterOrbitalSanityFx3", "test_g_minus",
     "fx3 inter-orbital pp sanity, g- direction"),
    ("test_bond_vs_ed_oracle", "TestPpInterOrbitalSanityFx3", "test_g_plus",
     "fx3 inter-orbital pp sanity, g+ direction"),
    # --- tests/test_bubble_kernel.py -------------------------------------
    ("test_bubble_kernel", "TestIrBondOracle",
     "test_matches_direct_lehmann_sum",
     "IR bond bubble vs a direct Lehmann sum -- 171 s on its own, the "
     "single most expensive test in the suite"),
    # --- tests/test_offsite_exchange_ed_longitudinal.py -- #181 Tier 2 ---
    ("test_offsite_exchange_ed_longitudinal",
     "TestOffsiteExchangeLongitudinalControls",
     "test_one_orbital_controls_and_the_fock_residual_identity",
     "off-site Exchange vs density-bond controls, one-orbital ring (ED)"),
    ("test_offsite_exchange_ed_longitudinal",
     "TestOffsiteExchangeLongitudinalControls",
     "test_two_orbital_interorbital_exchange_and_controls",
     "off-site inter-orbital Exchange longitudinal content vs ED, L=3 norb=2"),
    # --- tests/test_rpa_vs_ed_oracle.py ----------------------------------
    ("test_rpa_vs_ed_oracle", "TestRPAGeneralVsED", "test_first_order",
     "general-scheme RPA first-order response vs ED"),
    # --- tests/test_spinful_transverse_ed.py -- t_so granule campaign ----
    ("test_spinful_transverse_ed", "TestTask4GranuleCampaign",
     "test_norb2_coulombinter",
     "t_so granule campaign, norb=2 CoulombInter"),
    ("test_spinful_transverse_ed", "TestTask4GranuleCampaign",
     "test_norb2_exchange", "t_so granule campaign, norb=2 Exchange"),
    ("test_spinful_transverse_ed", "TestTask4GranuleCampaign",
     "test_norb2_hund", "t_so granule campaign, norb=2 Hund"),
    ("test_spinful_transverse_ed", "TestTask4GranuleCampaign",
     "test_norb2_pairhop_complex",
     "t_so granule campaign, norb=2 PairHop (complex phase)"),
    ("test_spinful_transverse_ed", "TestTask4GranuleCampaign",
     "test_norb2_pairlift", "t_so granule campaign, norb=2 PairLift"),
    ("test_spinful_transverse_ed", "TestTask4SVDSensitivityGate",
     "test_norb2_active_types_direction_set",
     "SVD sensitivity gate over the norb=2 active-type direction set"),
    ("test_spinful_transverse_ed", "TestTask4SVDSensitivityGate",
     "test_norb2_tso_mediated_direction_set",
     "SVD sensitivity gate over the t_so-mediated direction set"),
    ("test_spinful_transverse_ed", "TestTask4TsoZeroControl",
     "test_hund_returns_to_pass_zero_at_tso_zero",
     "t_so=0 isolation control for Hund"),
    ("test_spinful_transverse_ed", "TestTask4TsoZeroControl",
     "test_pairlift_returns_to_pass_zero_at_tso_zero",
     "t_so=0 isolation control for PairLift"),
)

#: The four fast-gate tests that exceed the 5 s rule by documented
#: exception (see the module docstring). Data, not prose, so the
#: consistency test can assert they are NOT skipped in the fast gate --
#: silently decorating one of them would otherwise go unnoticed.
FAST_GATE_EXCEPTIONS = (
    ("test_spinful_transverse_ed", "TestGateS0ExtractionAdjudication",
     "test_step3_dressed_extraction_matches_totalned"),
    ("test_spinful_transverse_ed", "TestGateS0ExtractionAdjudication",
     "test_step4_g2_g4_swap_mutation_must_fail"),
    ("test_spinful_transverse_ed", "TestGateS2ProductionRouting",
     "test_old_slice_then_solve_fails_the_s0_granule"),
    ("test_spinful_transverse_ed", "TestGateS2ProductionRouting",
     "test_production_passes_the_s0_granule"),
)

#: Attribute stamped on every method :func:`heavy` decorates. The
#: consistency test scans the discovered suite for it, which is how a
#: decorated-but-unregistered test gets caught.
HEAVY_MARKER_ATTR = "__hwave_heavy__"

SKIP_REASON = (
    "heavy opt-in test (>= ~5 s on its own); set {}=1 to run the full "
    "suite".format(FULL_SUITE_ENV)
)


def heavy(test_method):
    """Mark ``test_method`` opt-in: skipped unless ``HWAVE_FULL_TESTS``.

    ``unittest.skipUnless`` (NOT ``pytest.mark.skipif``) because
    ``python -m unittest discover`` is the gating runner. The returned
    object carries :data:`HEAVY_MARKER_ATTR` in both suite shapes -- the
    marker records the DECLARATION, whereas ``__unittest_skip__`` records
    the current run's decision, and the consistency test needs both.
    """
    decorated = unittest.skipUnless(full_suite_enabled(), SKIP_REASON)(
        test_method)
    setattr(decorated, HEAVY_MARKER_ATTR, True)
    return decorated


def heavy_test_ids():
    """The registry as ``"module.Class.method"`` strings (a frozenset)."""
    return frozenset(
        "{}.{}.{}".format(mod, cls, meth)
        for mod, cls, meth, _reason in HEAVY_TESTS)
