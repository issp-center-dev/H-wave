"""Unit tests for ``tests/equivalence_freeze_check.py`` -- the
read-only aggregation/validation collector that turns the calibration
workflow's measurement artifacts into a freeze decision.

Every test builds its own small, HAND-BUILT sample set (2 fake cells --
one ``Equiv``, one ``Diverges`` -- instead of the real 37-cell/21-Equiv
registry) so each validation-failure class and each reducer can be
pinned in isolation, fast, without depending on
``tests.equivalence_measure`` actually having run. Real-registry usage
is exercised end to end whenever the calibration workflow itself runs.
"""

from __future__ import annotations

import json
import os
import tempfile
import unittest

from tests.equivalence_cells import (
    Cell,
    Diverges,
    DivergingSpec,
    Equiv,
    ExecuteRun,
    FixtureSpec,
    ObservableSpec,
    SolverProof,
    Status,
)
from tests.equivalence_freeze_check import (
    DIAGNOSTIC_METRICS,
    GATING_RUNNERS,
    Sample,
    assert_amplification_holds,
    build_report,
    load_samples_from_dir,
    max_cell_seconds,
    max_module_total_seconds,
    max_residual,
    min_residual,
    paired_amplification_ratios,
    parse_sample_filename,
    unittest_gate_seconds,
    validate_samples,
)

SOURCE_SHA = "a" * 40


def _fixture() -> FixtureSpec:
    return FixtureSpec(
        input_dir="tests/equivalence_input/orb1",
        interactions={},
        T=1.0,
        mu=0.0,
        filling=None,
        CellShape=(1, 1, 1),
        SubShape=(1, 1, 1),
        Nmat=8,
        extra_params={},
        calc_type="ring",
        requested_scheme="general",
        enable_spin_orbital=False,
        extern=None,
    )


_EQUIV_CELL = Cell(
    cell_id="fake.equiv",
    fixture=_fixture(),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    comparison=Equiv(
        observables={"chi0q": ObservableSpec(comparator="identity", atol=1e-10, provenance="x")}
    ),
    required_observables=("chi0q",),
    interaction_class="onsite",
    notes="",
)

_DIVERGES_CELL = Cell(
    cell_id="fake.diverges",
    fixture=_fixture(),
    resolved_scheme="general",
    expected_spin_mode="spin-free",
    rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
    comparison=Diverges(
        diverging={
            "chiq": DivergingSpec(
                comparator="identity", ceiling=1e-10, regression_bound=1e-9, provenance="x"
            )
        },
        others={"chi0q": ObservableSpec(comparator="identity", atol=1e-10, provenance="x")},
        step5_issue="#1",
    ),
    required_observables=("chi0q", "chiq"),
    interaction_class="onsite",
    notes="",
)

FAKE_CELLS = (_EQUIV_CELL, _DIVERGES_CELL)


def _diag_records(residual=1e-16):
    return [
        {"diagnostic": metric, "fixture": fixture, "residual": residual}
        for metric in DIAGNOSTIC_METRICS
        for fixture in ("benign", "fc", "geev")
    ]


def _measurement_records(
    runner_python="3.10.18",
    source_sha=SOURCE_SHA,
    chi0q_equiv=1e-16,
    chi0q_diverges=1e-16,
    chiq_diverges=1e-16,
    diag_residual=1e-16,
    platform="Linux-x86_64",
):
    records = [
        {"cell": "fake.equiv", "seconds": 0.01},
        {"cell": "fake.diverges", "seconds": 0.01},
        {"cell": "fake.equiv", "observable": "chi0q", "residual": chi0q_equiv},
        {"cell": "fake.diverges", "observable": "chi0q", "residual": chi0q_diverges},
        {"cell": "fake.diverges", "observable": "chiq", "residual": chiq_diverges},
    ]
    records.extend(_diag_records(diag_residual))
    records.append(
        {
            "module_total_seconds": 0.5,
            "source_sha": source_sha,
            "runner": {
                "platform": platform,
                "python": runner_python,
                "numpy": "1.26.4",
                "scipy": "1.13.1",
            },
        }
    )
    return tuple(records)


def _good_sample_set(source_sha=SOURCE_SHA, **measurement_kwargs):
    samples = []
    for runner in GATING_RUNNERS:
        for n in (1, 2, 3):
            samples.append(
                Sample(
                    runner=runner,
                    invocation=n,
                    kind="measurement",
                    records=_measurement_records(
                        runner_python=runner + ".0", source_sha=source_sha, **measurement_kwargs
                    ),
                )
            )
        samples.append(
            Sample(
                runner=runner,
                invocation=1,
                kind="unittest",
                records=({"unittest_module_process_seconds": 10.0},),
            )
        )
    return samples


class TestParseSampleFilename(unittest.TestCase):
    def test_measurement_filename(self):
        self.assertEqual(parse_sample_filename("calib-3.9-1.json"), ("3.9", 1, "measurement"))
        self.assertEqual(parse_sample_filename("/a/b/calib-3.12-3.json"), ("3.12", 3, "measurement"))

    def test_unittest_filename(self):
        self.assertEqual(parse_sample_filename("unittest-3.10.json"), ("3.10", 1, "unittest"))

    def test_unrecognized_filename_raises(self):
        with self.assertRaises(ValueError):
            parse_sample_filename("something-else.json")


class TestValidateSamplesHappyPath(unittest.TestCase):
    def test_a_complete_valid_set_has_no_errors(self):
        errors = validate_samples(_good_sample_set(), SOURCE_SHA, FAKE_CELLS)
        self.assertEqual(errors, [])


class TestValidateSamplesCompleteness(unittest.TestCase):
    def test_missing_measurement_sample_reports_the_multiplicity(self):
        samples = _good_sample_set()
        samples.pop(0)  # drop one measurement sample -> one short
        errors = validate_samples(samples, SOURCE_SHA, FAKE_CELLS)
        expected = "expected {} measurement samples".format(3 * len(GATING_RUNNERS))
        self.assertTrue(any(expected in e for e in errors))

    def test_missing_unittest_sample_reports_the_multiplicity(self):
        samples = [s for s in _good_sample_set()
                   if not (s.kind == "unittest" and s.runner == GATING_RUNNERS[0])]
        errors = validate_samples(samples, SOURCE_SHA, FAKE_CELLS)
        expected = "expected {} unittest-timing samples".format(len(GATING_RUNNERS))
        self.assertTrue(any(expected in e for e in errors))

    def test_missing_cell_timing_line_reported(self):
        samples = _good_sample_set()
        broken = samples[0]
        pruned_records = tuple(r for r in broken.records if r.get("cell") != "fake.diverges" or "seconds" not in r or "observable" in r)
        samples[0] = Sample(runner=broken.runner, invocation=broken.invocation, kind=broken.kind, records=pruned_records)
        errors = validate_samples(samples, SOURCE_SHA, FAKE_CELLS)
        self.assertTrue(any("missing timing line(s)" in e and "fake.diverges" in e for e in errors))

    def test_missing_observable_residual_reported(self):
        samples = _good_sample_set()
        broken = samples[0]
        pruned_records = tuple(
            r for r in broken.records
            if not (r.get("cell") == "fake.diverges" and r.get("observable") == "chiq")
        )
        samples[0] = Sample(runner=broken.runner, invocation=broken.invocation, kind=broken.kind, records=pruned_records)
        errors = validate_samples(samples, SOURCE_SHA, FAKE_CELLS)
        self.assertTrue(any("missing observable residual" in e for e in errors))

    def test_missing_diagnostic_checkpoint_reported(self):
        samples = _good_sample_set()
        broken = samples[0]
        pruned_records = tuple(
            r for r in broken.records
            if not (r.get("diagnostic") == "assertion3" and r.get("fixture") == "fc")
        )
        samples[0] = Sample(runner=broken.runner, invocation=broken.invocation, kind=broken.kind, records=pruned_records)
        errors = validate_samples(samples, SOURCE_SHA, FAKE_CELLS)
        self.assertTrue(any("missing diagnostic checkpoint" in e for e in errors))

    def test_duplicate_diagnostic_record_reported(self):
        samples = _good_sample_set()
        broken = samples[0]
        dupe_record = next(
            r for r in broken.records
            if r.get("diagnostic") == "assertion3" and r.get("fixture") == "fc"
        )
        bumped_records = broken.records + (dupe_record,)
        samples[0] = Sample(runner=broken.runner, invocation=broken.invocation, kind=broken.kind, records=bumped_records)
        errors = validate_samples(samples, SOURCE_SHA, FAKE_CELLS)
        self.assertTrue(any("duplicate diagnostic record" in e for e in errors))


class TestValidateSamplesUniqueness(unittest.TestCase):
    def test_duplicate_runner_invocation_is_detected(self):
        samples = _good_sample_set()
        # Duplicate invocation 1 of the first gating runner as one extra
        # measurement sample (also breaks the count -- both errors are
        # expected together, this test only asserts the duplicate one is
        # present).
        dupe = next(s for s in samples if s.kind == "measurement"
                    and s.runner == GATING_RUNNERS[0] and s.invocation == 1)
        samples.append(dupe)
        errors = validate_samples(samples, SOURCE_SHA, FAKE_CELLS)
        self.assertTrue(any("duplicate (runner, invocation)" in e for e in errors))

    def test_runner_outside_the_gating_set_is_reported(self):
        samples = _good_sample_set()
        rogue = Sample(
            runner="3.13",
            invocation=1,
            kind="measurement",
            records=_measurement_records(runner_python="3.13.0"),
        )
        samples.append(rogue)
        errors = validate_samples(samples, SOURCE_SHA, FAKE_CELLS)
        self.assertTrue(any("is not one of the gating runners" in e for e in errors))


class TestValidateSamplesSourceSha(unittest.TestCase):
    def test_source_sha_mismatch_is_rejected(self):
        samples = _good_sample_set()
        broken = samples[0]
        bad_records = tuple(
            {**r, "source_sha": "wrong" * 8} if "source_sha" in r else r for r in broken.records
        )
        samples[0] = Sample(runner=broken.runner, invocation=broken.invocation, kind=broken.kind, records=bad_records)
        errors = validate_samples(samples, SOURCE_SHA, FAKE_CELLS)
        self.assertTrue(any("source_sha" in e and "!=" in e for e in errors))


class TestValidateSamplesErrorRecords(unittest.TestCase):
    def test_an_error_record_is_rejected(self):
        samples = _good_sample_set()
        broken = samples[0]
        bad_records = broken.records + ({"error": "boom", "cell": "fake.equiv", "phase": "cell"},)
        samples[0] = Sample(runner=broken.runner, invocation=broken.invocation, kind=broken.kind, records=bad_records)
        errors = validate_samples(samples, SOURCE_SHA, FAKE_CELLS)
        self.assertTrue(any("'error' record(s) present" in e for e in errors))


class TestValidateSamplesRunnerMetadata(unittest.TestCase):
    def test_inconsistent_runner_metadata_across_invocations_is_detected(self):
        samples = _good_sample_set()
        broken = samples[1]  # first gating runner, invocation 2
        assert broken.runner == GATING_RUNNERS[0] and broken.invocation == 2
        bad_records = tuple(
            {**r, "runner": {**r["runner"], "platform": "a-totally-different-platform"}}
            if "runner" in r else r
            for r in broken.records
        )
        samples[1] = Sample(runner=broken.runner, invocation=broken.invocation, kind=broken.kind, records=bad_records)
        errors = validate_samples(samples, SOURCE_SHA, FAKE_CELLS)
        self.assertTrue(any("inconsistent runner metadata" in e for e in errors))

    def test_runner_python_not_matching_the_declared_label_is_detected(self):
        samples = _good_sample_set()
        broken = samples[0]
        bad_records = tuple(
            {**r, "runner": {**r["runner"], "python": "2.7.18"}} if "runner" in r else r
            for r in broken.records
        )
        samples[0] = Sample(runner=broken.runner, invocation=broken.invocation, kind=broken.kind, records=bad_records)
        errors = validate_samples(samples, SOURCE_SHA, FAKE_CELLS)
        self.assertTrue(any("does not start with the declared runner label" in e for e in errors))


class TestLoadSamplesFromDir(unittest.TestCase):
    """``gh run download`` (no ``-n``) writes one SUBDIRECTORY per
    artifact, named after the artifact, containing that artifact's
    single file -- e.g. ``<dir>/calib-3.10-1/calib-3.10-1.json``. This
    pins ``load_samples_from_dir`` against that exact real-world shape
    (not just a flat directory of files).
    """

    def test_loads_the_gh_run_download_nested_layout(self):
        with tempfile.TemporaryDirectory() as tmp:
            for runner in (GATING_RUNNERS[0],):
                for n in (1, 2, 3):
                    sub = os.path.join(tmp, "calib-{}-{}".format(runner, n))
                    os.makedirs(sub)
                    with open(os.path.join(sub, "calib-{}-{}.json".format(runner, n)), "w") as f:
                        f.write(json.dumps({"module_total_seconds": 1.0, "source_sha": "x", "runner": {}}) + "\n")
                sub = os.path.join(tmp, "unittest-{}".format(runner))
                os.makedirs(sub)
                with open(os.path.join(sub, "unittest-{}.json".format(runner)), "w") as f:
                    f.write(json.dumps({"unittest_module_process_seconds": 5.0}) + "\n")

            samples = load_samples_from_dir(tmp)

        self.assertEqual(len(samples), 4)
        measurement = sorted((s.runner, s.invocation) for s in samples if s.kind == "measurement")
        r0 = GATING_RUNNERS[0]
        self.assertEqual(measurement, [(r0, 1), (r0, 2), (r0, 3)])
        unittest_ones = [s for s in samples if s.kind == "unittest"]
        self.assertEqual(len(unittest_ones), 1)
        self.assertEqual(unittest_ones[0].runner, r0)


class TestReducers(unittest.TestCase):
    def test_max_residual_is_the_global_max_across_runners_and_invocations(self):
        samples = _good_sample_set()
        # Bump one sample's fake.equiv/chi0q residual above the rest.
        target = samples[5]
        bumped = tuple(
            {**r, "residual": 7.5e-13} if r.get("cell") == "fake.equiv" and r.get("observable") == "chi0q" else r
            for r in target.records
        )
        samples[5] = Sample(runner=target.runner, invocation=target.invocation, kind=target.kind, records=bumped)
        self.assertEqual(max_residual(samples, "fake.equiv", "chi0q"), 7.5e-13)

    def test_max_residual_raises_on_no_matching_records(self):
        with self.assertRaises(ValueError):
            max_residual(_good_sample_set(), "no.such.cell", "chi0q")

    def test_min_residual_is_the_global_min(self):
        samples = _good_sample_set()
        target = next(s for s in samples if s.kind == "measurement"
                      and s.runner == GATING_RUNNERS[0] and s.invocation == 3)
        idx = samples.index(target)
        lowered = tuple(
            {**r, "residual": 2.0e-18} if r.get("cell") == "fake.diverges" and r.get("observable") == "chiq" else r
            for r in target.records
        )
        samples[idx] = Sample(runner=target.runner, invocation=target.invocation, kind=target.kind, records=lowered)
        self.assertEqual(min_residual(samples, "fake.diverges", "chiq"), 2.0e-18)

    def test_min_residual_raises_on_no_matching_records(self):
        with self.assertRaises(ValueError):
            min_residual(_good_sample_set(), "no.such.cell", "chiq")

    def test_max_cell_seconds(self):
        samples = _good_sample_set()
        target = next(s for s in samples if s.kind == "measurement" and s.runner == "3.11" and s.invocation == 2)
        idx = samples.index(target)
        bumped = tuple(
            {**r, "seconds": 12.5} if r.get("cell") == "fake.equiv" and "observable" not in r else r
            for r in target.records
        )
        samples[idx] = Sample(runner=target.runner, invocation=target.invocation, kind=target.kind, records=bumped)
        result = max_cell_seconds(samples)
        self.assertEqual(result["fake.equiv"], 12.5)
        self.assertEqual(result["fake.diverges"], 0.01)

    def test_max_module_total_seconds(self):
        samples = _good_sample_set()
        target = samples[2]
        bumped = tuple(
            {**r, "module_total_seconds": 99.0} if "module_total_seconds" in r else r
            for r in target.records
        )
        samples[2] = Sample(runner=target.runner, invocation=target.invocation, kind=target.kind, records=bumped)
        self.assertEqual(max_module_total_seconds(samples), 99.0)

    def test_unittest_gate_seconds_requires_one_sample_per_gating_runner(self):
        samples = _good_sample_set()
        self.assertEqual(unittest_gate_seconds(samples), 10.0)

        fewer = [s for s in samples
                 if not (s.kind == "unittest" and s.runner == GATING_RUNNERS[0])]
        with self.assertRaises(ValueError):
            unittest_gate_seconds(fewer)

    def test_unittest_gate_seconds_is_the_max_over_gating_runners(self):
        samples = _good_sample_set()
        target = next(s for s in samples if s.kind == "unittest" and s.runner == "3.12")
        idx = samples.index(target)
        samples[idx] = Sample(
            runner=target.runner, invocation=target.invocation, kind=target.kind,
            records=({"unittest_module_process_seconds": 118.0},),
        )
        self.assertEqual(unittest_gate_seconds(samples), 118.0)


class TestPairedAmplificationRatio(unittest.TestCase):
    def _samples_with_diag(self, per_runner_fc, per_runner_benign):
        samples = []
        for runner in GATING_RUNNERS:
            for n, (fc_val, benign_val) in enumerate(zip(per_runner_fc[runner], per_runner_benign[runner]), start=1):
                records = [
                    {"diagnostic": "assertion2", "fixture": "fc", "residual": fc_val},
                    {"diagnostic": "assertion2", "fixture": "benign", "residual": benign_val},
                ]
                samples.append(Sample(runner=runner, invocation=n, kind="measurement", records=tuple(records)))
        return samples

    def test_ratio_is_min_fc_over_max_benign_per_runner(self):
        per_runner_fc = {r: [1.0e-11, 2.0e-11, 3.0e-11] for r in GATING_RUNNERS}
        per_runner_benign = {r: [1.0e-13, 2.0e-13, 5.0e-13] for r in GATING_RUNNERS}
        samples = self._samples_with_diag(per_runner_fc, per_runner_benign)
        ratios = paired_amplification_ratios(samples, metric="assertion2")
        for runner in GATING_RUNNERS:
            expected = 1.0e-11 / max(5.0e-13, 1e-15)  # MIN(fc) / MAX(benign)
            self.assertAlmostEqual(ratios[runner], expected)

    def test_ratio_floors_the_denominator_at_1e_minus_15(self):
        per_runner_fc = {r: [1.0e-11] for r in GATING_RUNNERS}
        per_runner_benign = {r: [0.0] for r in GATING_RUNNERS}
        samples = self._samples_with_diag(per_runner_fc, per_runner_benign)
        ratios = paired_amplification_ratios(samples, metric="assertion2")
        for runner in GATING_RUNNERS:
            self.assertAlmostEqual(ratios[runner], 1.0e-11 / 1e-15)

    def test_assert_amplification_holds_passes_when_every_runner_clears_the_bar(self):
        ratios = {r: 50.0 for r in GATING_RUNNERS}
        assert_amplification_holds(ratios)  # must not raise

    def test_assert_amplification_holds_flags_a_failing_runner(self):
        ratios = {r: 50.0 for r in GATING_RUNNERS}
        ratios["3.11"] = 3.0  # below the 10x threshold
        with self.assertRaises(ValueError) as ctx:
            assert_amplification_holds(ratios)
        self.assertIn("3.11", str(ctx.exception))

    def test_assert_amplification_holds_flags_a_missing_runner(self):
        ratios = {r: 50.0 for r in GATING_RUNNERS if r != "3.10"}
        with self.assertRaises(ValueError) as ctx:
            assert_amplification_holds(ratios)
        self.assertIn("3.10", str(ctx.exception))
        self.assertIn("no amplification ratio", str(ctx.exception))


class TestBuildReport(unittest.TestCase):
    def test_invalid_set_reports_validation_failures_and_skips_reducers(self):
        samples = _good_sample_set()
        samples.pop(0)
        report = build_report(samples, SOURCE_SHA, FAKE_CELLS)
        self.assertIn("VALIDATION FAILED", report)
        self.assertIn(
            "expected {} measurement samples".format(3 * len(GATING_RUNNERS)),
            report)
        self.assertNotIn("fake.equiv / chi0q", report)

    def test_valid_set_reports_every_cell_observable_and_the_gate(self):
        report = build_report(_good_sample_set(), SOURCE_SHA, FAKE_CELLS)
        self.assertIn("VALIDATION OK", report)
        self.assertIn("fake.equiv / chi0q", report)
        self.assertIn("fake.diverges / chi0q", report)
        self.assertIn("fake.diverges / chiq", report)
        self.assertIn("120s freeze-time budget", report)
        self.assertIn("Conditioning amplification", report)


if __name__ == "__main__":
    unittest.main()
