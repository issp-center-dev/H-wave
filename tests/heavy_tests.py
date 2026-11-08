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
(The ``test_bubble_kernel.py`` case, ``TestIrBondOracle``, has since left
the registry: its oracle's 2e9-term Matsubara sum became the closed form,
and the test runs in well under a second.)

THE SELECTION RULE
------------------
**A test whose OWN runtime is roughly 5 s or more is opt-in.** Not the
module's runtime, not the class's -- the individual test method's. That is
why the entries below are per-method: several of these classes hold both a
30 s ED granule and a sub-second guard, and the guard stays in the fast
gate.

THE DOCUMENTED EXCEPTIONS (NONE AT PRESENT)
-------------------------------------------
:data:`FAST_GATE_EXCEPTIONS` lists tests that exceed the 5 s rule and STAY
in the fast gate anyway (not decorated, absent from :data:`HEAVY_TESTS`).
It is empty now. Until 2026-11-04 it held the four production-ROUTING
gates of ``test_spinful_transverse_ed`` (``TestGateS0ExtractionAdjudication``
and ``TestGateS2ProductionRouting``, two methods each): they pin that
``solve()`` takes the dressed-extraction path at all and that the old
slice-then-solve path still fails the S0 granule -- what an everyday
structural change actually breaks, and the in-place fix of issue #110 --
and at the ~31 s each they cost when carved out, keeping them on every
push was judged worth it. MEASURED on the CI runners on 2026-11-04 they
had grown to ~175 s EACH (12 minutes of the fast gate, paid in both of its
passes), and the FULL job now runs on every pull request too, so moving
them into :data:`HEAVY_TESTS` defers nothing past a merge: the guard still
runs before one, in `test-full`. The same measurement moved
``test_bond_vs_ed_oracle.TestCaseM``'s ``test_ph_g_plus`` /
``test_ph_g_minus`` (162 s each, never carved out -- simply over the rule).

And the inverse: two tests BELOW the 5 s rule are registered anyway --
``test_flex_second_order_ed_chain.TestG4``'s ``test_coulombinter_chain``
(~2.6 s) and ``test_hund_ising_chain_recorded`` (~1.9 s). They are the same
exact-diagonalization content, on the same two chain fixtures and at the
same Matsubara working point, as their ~18 s sibling
``test_interorbital_bond_chain``, which is over the rule; splitting the
three between the two shapes of the suite would mean a fast gate that
re-derives most of a chain ED without ever reaching the verdict the module
exists for. Their cost is deferred knowingly, with the sibling. The same
holds for ``test_so_pairing_ed.TestZeemanSzConserving``'s
``test_totaled_equals_sectored_L3`` (3.7 s measured on an idle machine,
10 s or more under load): the 12-mode engine cross-check of gate 0-S2
(ii), deferred with its ~64 s sibling; its L = 2 variant (0.03 s) keeps
the engine check in the fast gate.

Every OTHER entry below is over the rule on its own, measured. Two of them
sit close to it and carry the measurement in their reason string, so that a
reader who wants to re-check the selection does not have to re-time them
blind.

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
    # --- tests/test_bond_half_grid.py -- the half bosonic grid (#223) -----
    ("test_bond_half_grid", "TestPairingOracleIR", "test_mixing_and_batches",
     "twelve two-orbital bond-gate FLEX solves with the in-process IR pairing "
     "(full and half grid, both mixing schemes, three pinned batches), 6.4 s "
     "measured; the fixture sweep of the same class stays in the fast gate"),
    ("test_bond_half_grid", "TestPairingOracleUniform", "test_mixing_and_batches",
     "sixteen two-orbital bond-gate FLEX solves with the in-process uniform "
     "pairing of both channels (full and half grid, both mixing schemes, four "
     "pinned batches), 5.4 s measured; the fixture sweep of the same class "
     "stays in the fast gate"),
    # --- tests/test_bond_longitudinal_ed.py -- longitudinal bond ED granules
    ("test_bond_longitudinal_ed", "TestLongitudinalBondGranulesCaseM",
     "test_coulombinter_interorbital",
     "case-M off-site CoulombInter bond blocks vs ED (S and C)"),
    ("test_bond_longitudinal_ed", "TestLongitudinalBondGranulesCaseM",
     "test_exchange_interorbital",
     "case-M off-site Exchange bond blocks vs ED (S and C)"),
    ("test_bond_longitudinal_ed", "TestLongitudinalBondGranulesCaseM",
     "test_g3_bond_blocks_are_load_bearing_interorbital",
     "G3 on case M: bond blocks reproduce ED; Hartree-only misses by >= 10 tol"),
    ("test_bond_longitudinal_ed", "TestLongitudinalBondGranulesCaseM",
     "test_hund_interorbital",
     "case-M off-site Hund bond blocks vs ED (S and C)"),
    ("test_bond_longitudinal_ed", "TestLongitudinalBondGranulesCaseM",
     "test_ising_interorbital",
     "case-M off-site Ising bond blocks vs ED (S and C)"),
    ("test_bond_longitudinal_ed", "TestLongitudinalBondGranulesCaseM",
     "test_pairhop_interorbital",
     "case-M off-site PairHop bond blocks vs ED (S and C)"),
    ("test_bond_longitudinal_ed", "TestLongitudinalBondGranulesFx5",
     "test_coulombinter_g1_g2",
     "fx5 off-site CoulombInter bond blocks vs ED, shells 1 and 2"),
    ("test_bond_longitudinal_ed", "TestLongitudinalBondGranulesFx5",
     "test_exchange_g1",
     "fx5 off-site Exchange bond blocks vs ED"),
    ("test_bond_longitudinal_ed", "TestLongitudinalBondGranulesFx5",
     "test_g3_bond_blocks_are_load_bearing",
     "G3: bond blocks reproduce ED; Hartree-only misses by >= 10 tol"),
    ("test_bond_longitudinal_ed", "TestLongitudinalBondGranulesFx5",
     "test_hund_g1",
     "fx5 off-site Hund bond blocks vs ED"),
    ("test_bond_longitudinal_ed", "TestLongitudinalBondGranulesFx5",
     "test_ising_g1",
     "fx5 off-site Ising bond blocks vs ED"),
    ("test_bond_longitudinal_ed", "TestLongitudinalBondGranulesFx5",
     "test_pairhop_g1",
     "fx5 off-site PairHop bond blocks vs ED"),
    # --- tests/test_bond_pairing_ed.py -- bond pairing vertex vs ED (#204) -
    ("test_bond_pairing_ed", "TestPairingVertexAgainstED",
     "test_first_and_second_order",
     "pair susceptibility O(lambda), O(lambda^2) vs ED, both spin sectors"),
    ("test_bond_pairing_ed", "TestPairingVertexAgainstED",
     "test_free_pair_susceptibility",
     "free pair susceptibility vs ED (shares the ED fixture)"),
    ("test_bond_pairing_ed", "TestPairingVertexAgainstED",
     "test_onsite_interorbital_row",
     "O(lambda), O(lambda^2) vs ED with an on-site inter-orbital row (#262)"),
    ("test_bond_pairing_ed", "TestPairingVertexAgainstED",
     "test_opposite_placement_is_missed",
     "opposite-placement ED misses the production ladder"),
    ("test_bond_pairing_ed", "TestPairingVertexAgainstED",
     "test_second_order_needs_the_fluctuation_vertex",
     "ladder without the fluctuation vertex misses at O(lambda^2)"),
    ("test_bond_pairing_ed", "TestPairingVertexAgainstED",
     "test_spin_coefficients_are_not_interchangeable",
     "interchanged singlet/triplet coefficients miss at O(lambda^2)"),
    # --- tests/test_bond_supercell.py -- multi-orbital supercell equivalence (#239)
    ("test_bond_supercell", "TestDynamicEntries",
     "test_in_process_and_post_processing",
     "two FLEX bond-gate runs plus the post-processing of the supercell archive"),
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
    ("test_bond_vs_ed_oracle", "TestCaseM", "test_ph_g_minus",
     "case-M ph adjudication, g- direction -- 162 s on CI (2026-11-04)"),
    ("test_bond_vs_ed_oracle", "TestCaseM", "test_ph_g_plus",
     "case-M ph adjudication, g+ direction -- 162 s on CI (2026-11-04)"),
    ("test_bond_vs_ed_oracle", "TestPpInterOrbitalSanityFx3", "test_g_minus",
     "fx3 inter-orbital pp sanity, g- direction"),
    ("test_bond_vs_ed_oracle", "TestPpInterOrbitalSanityFx3", "test_g_plus",
     "fx3 inter-orbital pp sanity, g+ direction"),
    # --- tests/test_eliashberg_bond_entry.py -- the two pairing entries ---
    ("test_eliashberg_bond_entry", "TestPostProcessingRuns",
     "test_ir_lambda_arm_matches_uniform",
     "the IR arm of the entry-level agreement case: two real bond-gate FLEX "
     "runs (4x4, Nmat 64 and 128) plus eight hwave_sc solves on their "
     "archives, one per channel per Matsubara basis, comparing lambda_ir with "
     "lambda_uniform on the same archive -- it also needs sparse-ir, which "
     "the fast gate does not require. Its uniform sibling "
     "(test_uniform_lambda_and_outputs, 2.2 s) and the two-entry agreement "
     "case (TestInProcess.test_inprocess_equals_postprocessing, 1.2 s) are "
     "deliberately NOT registered: they are the only end-to-end runs of the "
     "hwave_sc bond entry, and a fast gate without them would exercise its "
     "refusals and nothing else"),
    # --- tests/test_eliashberg_dynamic_so.py -- spin-orbital dynamic (#284)
    ("test_eliashberg_dynamic_so", "TestKernelAgainstValidatedKernel",
     "test_exact_interacting_green_function",
     "the production spin-orbital dynamic operator equals the ED-validated "
     "kernel on the exact interacting Green function of the L = 3 spin-mixing "
     "chain (one ED solve of the fixture shared with test_so_pairing_ed, 13 s "
     "measured when it is the first method to build it, free when that "
     "module ran first in the process); the free-G "
     "sibling stays in the fast gate"),
    ("test_eliashberg_dynamic_so", "TestLMRetry",
     "test_retry_matches_the_dense_reduced_spectrum",
     "the which='LM' retry of the spin-orbital dynamic solver against the dense "
     "reduced spectrum at Nmat = 32 (a 4096 x 4096 dense eigenproblem, ~22 s); "
     "the ARPACK-referenced sibling stays in the fast gate"),
    # --- tests/test_flex_bond_device_store.py -- device block store (#275) -
    ("test_flex_bond_device_store", "TestKernelMatvecsAgree", "test_ir",
     "twelve bond-gate FLEX solves with the in-process IR pairing on the fake "
     "device (host vs device store; singlet, triplet, both; full and half "
     "grid), 6.2 s measured; the uniform-grid sibling stays in the fast gate"),
    # --- tests/test_flex_bond_gpu.py -- CPU vs GPU equivalence -----------
    ("test_flex_bond_gpu", "TestEndToEndEquivalence",
     "test_two_orbital_bond_chain",
     "two whole bond-gate FLEX solves of the two-orbital fixture (4x4, "
     "Nmat = 8, EPS = 1e-12, 45 iterations each), one per array backend, "
     "compared end to end -- 4.5 s measured on an A100, so it sits right "
     "at the rule; on a machine WITHOUT a CUDA device the whole class is "
     "skipped from setUpClass anyway, and the fast gate pays nothing"),
    # --- tests/test_flex_bond_memory.py -- #181 Phase B memory validation --
    ("test_flex_bond_memory", "TestSubprocessRSS",
     "test_subprocess_rss_below_estimate_plus_allowance",
     "bond-resolved FLEX peak RSS of a subprocess vs the preflight estimate "
     "plus the calibrated interpreter baseline (6x6, Nmat=64, depth 8, "
     "full output)"),
    # --- tests/test_flex_bond_onari_trend.py -- #181 Phase B trend milestone --
    ("test_flex_bond_onari_trend", "TestFlexBondOnariTrend",
     "test_regenerated_greens_reproduce_the_pinned_lambda",
     "regenerate the Hartree-Fock bond FLEX greens of the Onari-type trend "
     "milestone (hours; HWAVE_RUN_SLOW_FIXTURES=1) and re-derive lambda_t"),
    # --- tests/test_flex_second_order_ed_chain.py -- #181 second-order G4 --
    ("test_flex_second_order_ed_chain", "TestG4",
     "test_coulombinter_chain",
     "G4: exact diagonalisation of the L = 4 single-orbital chain pins the "
     "off-site second order of CoulombInter (O(V^2), O(U V), O(U^2)) for "
     "the local path plus the dropped class and for the bond gate"),
    ("test_flex_second_order_ed_chain", "TestG4",
     "test_hund_ising_chain_recorded",
     "G4 (recorded): off-site Hund and Ising on the L = 4 chain against "
     "exact diagonalisation, printed rather than asserted"),
    ("test_flex_second_order_ed_chain", "TestG4",
     "test_hund_times_offsite_v_chain",
     "G4: the mixed O(J V) coefficient of on-site Hund against the "
     "asymmetric inter-orbital off-site bond on the L = 3 two-orbital "
     "chain, and the vanishing of the mixed dropped class"),
    ("test_flex_second_order_ed_chain", "TestG4",
     "test_interorbital_bond_chain",
     "G4: the L = 3 two-orbital chain with an orbital-asymmetric "
     "inter-orbital off-site bond adjudicates the off-site vertex "
     "orientation and the bond gate, both asserted -- the gate against a "
     "floor the test measures itself and against the exact oracle"),
    # --- tests/test_flex_second_order_ed_onsite.py -- #181 second-order G3 --
    ("test_flex_second_order_ed_onsite", "TestG3",
     "test_every_type_and_pair",
     "G3: exact diagonalisation of the single-site three-orbital model pins "
     "the on-site second order of all seven interaction types and seven "
     "mixed pairs (~35 production maps at Nmat = 4096)"),
    # --- tests/test_flex_second_order_sopt.py -- #181 second-order gates --
    ("test_flex_second_order_sopt", "TestG2Heavy",
     "test_b_covering_set_equals_the_local_oracle",
     "G2 (b): end-to-end second-order coefficients of the general path over "
     "the covering set of 10 pure types and 12 mixed pairs, against the "
     "independent real-space oracle (Richardson ladder, ~600 maps)"),
    ("test_flex_second_order_sopt", "TestG2Heavy",
     "test_b_dropped_class_load_bearing_for_v",
     "G2 (b): the off/off uncrossed class the local weighting drops is not "
     "negligible for the off-site V, so the covering-set gate really pins "
     "the local weighting"),
    ("test_flex_second_order_sopt", "TestG2Heavy",
     "test_c_gate_on_reproduces_the_full_oracle",
     "G2 (c): with the bond gate on the second order is the FULL oracle for "
     "off-site CoulombInter; the off-site Hund/Ising outcome is recorded"),
    ("test_flex_sigma_hermitian", "TestColdStartReproduction",
     "test_cold_start_remnant_and_split_seed",
     "issue #281: two cold-started FLEX runs of the dimer past the Stoner "
     "instability (6x6, Nmat 32, up to 1000 linear-mixing iterations each), "
     "with and without the Hermitian projection, then the split conversion "
     "and a flex_hartree_fock run on the seed -- 4.3 s measured, at the rule; "
     "the fast gate keeps the deterministic seeded tests of the same module"),
    # --- tests/test_flex_spin_orbital_ed.py -- #83 spin-orbital FLEX --
    ("test_flex_spin_orbital_ed", "TestSecondOrderWithSOC",
     "test_every_type",
     "exact diagonalisation of the single-site three-orbital model with a "
     "spin-mixing one-body block pins the second order of the spin-orbital "
     "FLEX for all seven interaction types (28 production maps at "
     "Nmat = 4096, ~20 s)"),
    ("test_flex_spin_orbital_ed", "TestZeemanSecondOrder",
     "test_zeeman_without_and_with_soc",
     "gate 0-F2 (#88): the same second-order ED comparison with a Zeeman field, "
     "without and with the spin-mixing block, the field through Transfer and "
     "through Extern (all seven types, Nmat = 4096, 58 s measured); the "
     "three-type Nmat = 512 smoke variant stays in the fast gate"),
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
    # --- tests/test_sc_spin_orbital_bond_full.py -- #270 spin-orbital bond --
    ("test_sc_spin_orbital_bond_full", "TestAdmission", "test_phase_rows_iteration",
     "the traced per-phase memory rows of the static spin-orbital bond path in "
     "the iteration solver mode at Nmat 32 (S_z-packed Zeeman fixture, up to "
     "1000 power-iteration steps); below the rule, kept out of the fast set by "
     "its time budget; the fast module keeps the Nmat 8 smoke; macOS 0.43 s, "
     "Linux 1.0 s"),
    ("test_sc_spin_orbital_bond_full", "TestLayoutSelectionThroughEntry",
     "test_layout_selection_through_entry",
     "the public-entry layout pair of the static spin-orbital bond path (the "
     "S_z-conserving Zeeman fixture in the automatic S_z-packed and the forced "
     "dense layout, full spectrum); below the rule, kept out of the fast set by "
     "its time budget; macOS 0.10 s, Linux 0.14 s"),
    ("test_sc_spin_orbital_bond_full", "TestReportingThroughEntry", "test_state_rows",
     "the forced sector-solver states of the static spin-orbital bond path through "
     "the public entry (rows 2-7 S_z-packed, rows 3 and 7 dense, three solver modes "
     "each, 24 runs); below the rule, kept out of the fast set by its time budget; "
     "macOS 1.0 s, Linux 2.4 s"),
    ("test_sc_spin_orbital_bond_full", "TestZeroFieldThroughEntry",
     "test_zero_field_limit_through_entry",
     "the zero-field limit of the static spin-orbital bond path through the public "
     "entry (4x4, B = 1 and B = 5, one and two orbitals, both layouts, five fields, "
     "complete listings; 40 runs, plus the Extern/Transfer route pair and the "
     "operator at the first two fields for the isolation of the selected clusters); "
     "the fast module keeps the 4x2 operator-level variant; macOS 2.7 s, "
     "Linux 6.3 s"),
    # --- tests/test_so_bond_chi_ed.py -- spinful bond correlators vs ED ---
    ("test_so_bond_chi_ed", "TestChiFirstOrderED", "test_fixture_ii",
     "the first-order longitudinal and transverse bond correlators of the "
     "two-orbital three-cell chain vs exact diagonalization (12 modes, "
     "largest sector 924, two steps), 19 sketches, 15.6 s measured on macOS, "
     "19.0 s on Linux (numpy 1.26); the five-cell chain and its Zeeman variant "
     "stay in the fast gate"),
    # --- tests/test_so_bond_conventions.py -- spinful bond vertex vs ED ----
    ("test_so_bond_conventions", "TestFilingAgainstED", "test_fixture_i_per_type",
     "the spinful bond vertex of each interaction type alone vs the exact "
     "first-order vertex part on the five-cell chain, every q and all 14 "
     "sketches (seven exact diagonalizations at two steps), 5.2 s measured on "
     "macOS, 10.4 s on Linux (numpy 1.26); the combined fixture stays in the fast gate"),
    ("test_so_bond_conventions", "TestFilingAgainstED", "test_fixture_ii_combined",
     "the spinful bond vertex of the two-orbital three-cell chain vs the exact "
     "first-order vertex part (12 modes, largest sector 924, two steps), all "
     "14 sketches, 11.9 s measured on macOS, 13.9 s on Linux (numpy 1.26)"),
    # --- tests/test_so_bond_pairing_ed.py -- spinful bond pairing vs ED ---
    ("test_so_bond_pairing_ed", "TestSOBondPairingAgainstED",
     "test_first_and_second_order_two_sided",
     "the spinful bond pair susceptibility at O(lambda) and O(lambda^2) vs "
     "exact diagonalization of the two-orbital three-cell chain with on-site "
     "and nearest-neighbour interactions (ten ED solves, production at "
     "nmat = 128 and 256 for every reading), two-sided, 98.7 s measured on "
     "macOS, 233 s on Linux (numpy 1.26); the one-orbital smoke stays in the fast gate"),
    # --- tests/test_so_bond_pairing_oracle.py -- spinful bond pairing, 1st order
    ("test_so_bond_pairing_oracle", "TestFirstOrderSOBondPairing", "test_all_types_norb_1_2_3",
     "the static spinful bond pairing kernel at first order vs the independent "
     "spin-orbital construction, every off-site type at both orientations and "
     "CoulombIntra, norb_phys 1, 2 and 3 (39 cases, four production chains "
     "each), 3.7 s measured on macOS, 7.7 s on Linux (numpy 1.26); heavy because "
     "of the fast-gate budget, not its own runtime; the two-orbital subset "
     "stays in the fast gate"),
    # --- tests/test_so_pairing_ed.py -- spin-orbital pairing vertex vs ED (#83)
    ("test_so_pairing_ed", "TestSOPairingVertexAgainstED",
     "test_first_and_second_order",
     "spin-orbital pair susceptibility O(lambda), O(lambda^2) vs a spin-mixing ED"),
    ("test_so_pairing_ed", "TestSOPairingVertexAgainstED",
     "test_free_pair_susceptibility",
     "free spin-orbital pair susceptibility vs ED (shares the ED fixture)"),
    ("test_so_pairing_ed", "TestSOPairingVertexAgainstED",
     "test_opposite_reading_is_missed",
     "opposite reading of the fluctuation term misses at O(lambda^2)"),
    ("test_so_pairing_ed", "TestSOPairingVertexAgainstED",
     "test_second_order_needs_the_fluctuation_vertex",
     "spin-orbital ladder without the fluctuation term misses at O(lambda^2)"),
    ("test_so_pairing_ed", "TestZeemanSpinFlipModel",
     "test_first_and_second_order",
     "gate 0-S2 (i, #88): the spin-mixing chain with a Zeeman field, pair "
     "susceptibility O(lambda), O(lambda^2) vs TotalNED (five new ED solves, "
     "63 s measured)"),
    ("test_so_pairing_ed", "TestZeemanSzConserving",
     "test_first_and_second_order",
     "gate 0-S2 (ii, #88): the S_z-conserving chain with a Zeeman field, pair "
     "susceptibility O(lambda), O(lambda^2) vs ED (five new ED solves, 64 s "
     "measured)"),
    ("test_so_pairing_ed", "TestZeemanSzConserving",
     "test_totaled_equals_sectored_L3",
     "gate 0-S2 (ii, #88): TotalNED equals SectorED on the 12-mode L = 3 "
     "Zeeman chain (3.7 s measured; registered with its ED siblings, the "
     "L = 2 variant stays in the fast gate)"),
    # --- tests/test_spinful_transverse_ed.py -- production routing gates --
    # the former FAST_GATE_EXCEPTIONS (see the module docstring): ~175 s
    # each on CI (2026-11-04); test-full runs on every pull request, so
    # the routing guard still runs before a merge. Since the four share
    # one cached ED reference (_s0_ed_derivative), the first of them to
    # run pays ~19 s locally and the others under 1 s; each one ALONE
    # still costs ~19 s, over the rule, so all four stay registered.
    ("test_spinful_transverse_ed", "TestGateS0ExtractionAdjudication",
     "test_step3_dressed_extraction_matches_totalned",
     "gate S0: dressed extraction matches the total-N ED (routing guard)"),
    ("test_spinful_transverse_ed", "TestGateS0ExtractionAdjudication",
     "test_step4_g2_g4_swap_mutation_must_fail",
     "gate S0: the g2/g4 swap mutation must fail (routing guard)"),
    ("test_spinful_transverse_ed", "TestGateS2ProductionRouting",
     "test_old_slice_then_solve_fails_the_s0_granule",
     "gate S2: the old slice-then-solve path fails the S0 granule"),
    ("test_spinful_transverse_ed", "TestGateS2ProductionRouting",
     "test_production_passes_the_s0_granule",
     "gate S2: production solve() passes the S0 granule (issue #110)"),
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
    # --- tests/test_uhfk_orientation.py -- #193 end-to-end UHFk gate -----
    ("test_uhfk_orientation", "TestUHFkOrientation", "test_all_types_bond",
     "three UHFk runs per case x two inter-orbital Ising/Exchange/PairLift "
     "cases (normal and spin-orbital) against the reference revision on "
     "the reversed declaration (6.1-6.4 s measured)"),
    ("test_uhfk_orientation", "TestUHFkOrientation", "test_bond_variants",
     "three UHFk runs per case x three cases (Fock off, spin-orbital) "
     "against the reference revision on the reversed declaration "
     "(5.7 s measured -- the closest entry to the rule)"),
    ("test_uhfk_orientation", "TestUHFkOrientation", "test_pairhop_complex",
     "three UHFk runs per case x eight complex-PairHop cases (norb 2 and "
     "1, Fock on and off, normal and spin-orbital) against the reference "
     "revision on the reversed declaration"),
)

#: Fast-gate tests that exceed the 5 s rule by documented exception (see
#: the module docstring). Data, not prose, so the consistency test can
#: assert they are NOT skipped in the fast gate -- silently decorating one
#: of them would otherwise go unnoticed. Empty since 2026-11-04: the four
#: routing gates it held moved into HEAVY_TESTS (~175 s each on CI).
FAST_GATE_EXCEPTIONS = ()

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


#: The registry's own consistency gate, which the FULL job runs alongside
#: the registered tests (it is the one fast test whose meaning depends on
#: the full-suite switch: it drives both shapes in child processes).
REGISTRY_TEST_FILE = "tests/test_heavy_tests_registry.py"


def heavy_node_ids(shard=None):
    """The registry as pytest node ids, ``tests/<module>.py::Class::method``,
    in registry order. The FULL CI job selects exactly these (plus
    :data:`REGISTRY_TEST_FILE`), so the list is derived here at run time
    and never written out by hand.

    ``shard`` (1-based) restricts the list to the files of that shard of
    ``tests/heavy_shards.py``; a registered file that is in no shard is an
    error, so it can never be left out of every job."""
    ids = [("tests/{}.py::{}::{}".format(mod, cls, meth), mod)
           for mod, cls, meth, _reason in HEAVY_TESTS]
    if shard is None:
        return [node for node, _mod in ids]
    from tests import heavy_shards

    if not 1 <= shard <= heavy_shards.N_SHARDS:
        raise ValueError("shard {} out of 1..{}".format(shard, heavy_shards.N_SHARDS))
    where = heavy_shards.shard_of_module()
    unplaced = sorted({mod for _node, mod in ids if mod not in where})
    if unplaced:
        raise ValueError(
            "registered heavy test file(s) in no shard of tests/heavy_shards.py: "
            "{}".format(unplaced))
    return [node for node, mod in ids if where[mod] == shard]


def run_list(shard=None):
    """What one FULL job runs: :func:`heavy_node_ids` of the shard, plus
    :data:`REGISTRY_TEST_FILE` in shard 1 (and in the unsharded run)."""
    out = heavy_node_ids(shard)
    if shard in (None, 1):
        out.append(REGISTRY_TEST_FILE)
    return out


def missing_from_collection(expected, collect_output):
    """The ``expected`` node ids that do NOT appear in ``collect_output``,
    the text ``pytest --collect-only -q`` printed (one node id per line).

    A registered test that was renamed or removed without updating the
    registry must fail the FULL job rather than drop out of it silently.
    pytest does exit non-zero on a node id it cannot resolve, but that
    behaviour is the runner's, not this project's, and its message names
    only the first miss; the comparison here is explicit and names every
    missing id."""
    collected = set(line.strip() for line in collect_output.splitlines())
    return [node for node in expected if node not in collected]


def _main(argv=None):
    """Command line used by the FULL CI job.

    ``python -m tests.heavy_tests node-ids [--shard N]``
        print the node ids to run, one per line (:func:`run_list`).
    ``python -m tests.heavy_tests check-collected FILE [--shard N]``
        exit non-zero, naming them, if any registered node id (of shard N)
        is absent from FILE (the output of ``pytest --collect-only -q``).
    """
    import argparse
    import sys

    parser = argparse.ArgumentParser(prog="python -m tests.heavy_tests")
    sub = parser.add_subparsers(dest="cmd", required=True)
    ids_p = sub.add_parser("node-ids")
    ids_p.add_argument("--shard", type=int, default=None)
    chk = sub.add_parser("check-collected")
    chk.add_argument("collect_output")
    chk.add_argument("--shard", type=int, default=None)
    args = parser.parse_args(argv)

    if args.cmd == "node-ids":
        for node in run_list(args.shard):
            print(node)
        return 0
    ids = heavy_node_ids(args.shard)
    with open(args.collect_output) as fh:
        text = fh.read()
    missing = missing_from_collection(ids, text)
    if missing:
        sys.stderr.write(
            "registered heavy test(s) NOT collected by pytest -- renamed or "
            "removed without updating tests/heavy_tests.py:\n  "
            + "\n  ".join(missing) + "\n")
        return 1
    print("all {} registered heavy tests collected".format(len(ids)))
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
