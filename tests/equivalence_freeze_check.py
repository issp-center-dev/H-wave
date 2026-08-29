"""The read-only aggregation/validation collector for the RPA/FLEX
equivalence-table calibration workflow
(``.github/workflows/equivalence-calibration.yml``).

The calibration workflow (``.github/workflows/equivalence-calibration.
yml``) uploads, per gating-runner job (``ubuntu-latest`` x Python
3.9/3.10/3.11/3.12): three ``python -m tests.equivalence_measure``
invocations (``calib-<py>-<n>.json``, ``n`` in 1..3) plus one
``python -m unittest tests.test_rpa_flex_equivalence_table`` timing
invocation (``unittest-<py>.json``). Whoever runs the calibration
downloads all sixteen artifacts (``gh run download``) for a single
successful run matching a source commit ``S``. This module parses
those files (or, for testing, hand-built in-memory equivalents),
VALIDATES completeness against the aggregation rules implemented
below, and computes every reducer those rules define. Its printed
report (``build_report``) is the SOLE input to any freeze decision --
a bound is never hand-computed from the raw artifacts.

This module is read-only: it never asserts numerical equivalence
itself (that already happened when ``equivalence_measure`` ran) and it
writes nothing. It DOES import ``tests.equivalence_cells`` (a pure,
side-effect-free data module -- not a violation of that module's own
"no unittest, no solver imports" constraint, which is about what
``equivalence_cells.py`` itself may import, not who may import it).
"""

from __future__ import annotations

import glob
import json
import os
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from tests.equivalence_cells import CELLS, Diverges, Equiv

# The four gating runners (Global Constraints: "the classification set
# is STRICTLY the four gating runners"). Used both as the expected
# ``Sample.runner`` label set and to size the multiplicity checks.
GATING_RUNNERS: Tuple[str, ...] = ("3.9", "3.10", "3.11", "3.12")

# The diagnostic apparatus records one named scalar per metric per
# fixture: 11 metrics x 3 fixtures = 33 records per measurement sample
# (spec 2026-08-28-mu-green-seam-160: every scalar is its OWN record).
DIAGNOSTIC_METRICS: Tuple[str, ...] = (
    "assertion1", "assertion2", "assertion3", "assertion4", "assertion5",
    "counter_cross_at_mu_rpa", "counter_cross_at_mu_flex",
    "number_residual_rpa", "number_residual_flex",
    "dyson_residual_eigenbasis", "dyson_residual_inv",
)
DIAGNOSTIC_FIXTURES: Tuple[str, ...] = ("benign", "fc", "geev")

# Global Constraints: the conditioning amplification threshold.
AMPLIFICATION_THRESHOLD = 10.0


# ---------------------------------------------------------------------------
# Sample: one (runner, invocation) artifact's parsed JSON-lines content.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Sample:
    """One downloaded artifact: ``kind="measurement"`` for a
    ``calib-<py>-<n>.json`` (one of the 3 per-runner ``equivalence_
    measure`` invocations) or ``kind="unittest"`` for a
    ``unittest-<py>.json`` (the single per-runner unittest-timing
    invocation). ``runner`` is the gating Python-version label
    (``"3.9"``..``"3.12"``); ``invocation`` is 1-based (always 1 for
    ``kind="unittest"``, since the workflow runs that step once per
    job). ``records`` is the tuple of parsed JSON objects, one per
    line, in file order.
    """

    runner: str
    invocation: int
    kind: str  # "measurement" | "unittest"
    records: Tuple[dict, ...]


_FILENAME_RE = re.compile(
    r'^calib-(?P<py>[0-9.]+)-(?P<n>[0-9]+)\.json$|'
    r'^unittest-(?P<upy>[0-9.]+)\.json$'
)


def parse_sample_filename(filename: str) -> Tuple[str, int, str]:
    """Parse ``calib-<py>-<n>.json`` -> ``(py, n, "measurement")`` or
    ``unittest-<py>.json`` -> ``(py, 1, "unittest")``, matching the
    exact naming the workflow's ``if: always()`` upload steps produce.
    Raises ``ValueError`` on any other name -- a naming drift between
    the workflow and this parser must fail loudly, not silently skip an
    artifact.
    """

    base = os.path.basename(filename)
    m = _FILENAME_RE.match(base)
    if not m:
        raise ValueError(
            "parse_sample_filename: {!r} does not match "
            "'calib-<py>-<n>.json' or 'unittest-<py>.json'".format(base)
        )
    if m.group("py") is not None:
        return m.group("py"), int(m.group("n")), "measurement"
    return m.group("upy"), 1, "unittest"


def _read_json_lines(path: str) -> Tuple[dict, ...]:
    with open(path, "r", encoding="utf-8") as f:
        return tuple(json.loads(line) for line in f if line.strip())


def load_sample_file(path: str, runner: Optional[str] = None,
                      invocation: Optional[int] = None,
                      kind: Optional[str] = None) -> Sample:
    """Load one artifact file into a ``Sample``. ``runner``/
    ``invocation``/``kind`` are inferred from the filename
    (``parse_sample_filename``) when not given explicitly -- tests
    that build synthetic files under an arbitrary name can still
    override them.
    """

    inferred_runner, inferred_invocation, inferred_kind = parse_sample_filename(path)
    return Sample(
        runner=runner if runner is not None else inferred_runner,
        invocation=invocation if invocation is not None else inferred_invocation,
        kind=kind if kind is not None else inferred_kind,
        records=_read_json_lines(path),
    )


def load_samples_from_dir(directory: str) -> List[Sample]:
    """Load every ``calib-*.json``/``unittest-*.json`` file found
    UNDER ``directory``, at any depth. Covers both the flat-directory
    case and the shape ``gh run download`` actually produces by
    default: one SUBDIRECTORY per uploaded artifact (named after the
    artifact, e.g. ``<directory>/calib-3.9-1/calib-3.9-1.json``), each
    containing that artifact's single file. Raises ``ValueError`` (via
    ``load_sample_file`` -> ``parse_sample_filename``) on any matched
    file whose name does not fit the expected pattern -- a stray file
    must fail loudly, not be silently skipped or silently included.
    """

    paths = sorted(
        set(
            glob.glob(os.path.join(directory, "**", "calib-*.json"), recursive=True)
            + glob.glob(os.path.join(directory, "**", "unittest-*.json"), recursive=True)
        )
    )
    return [load_sample_file(p) for p in paths]


# ---------------------------------------------------------------------------
# Expected shape, derived from the registry (pure data, no side effects).
# ---------------------------------------------------------------------------


def _observable_pairs(cells: Sequence) -> set:
    pairs = set()
    for cell in cells:
        comparison = cell.comparison
        if isinstance(comparison, Equiv):
            for observable in comparison.observables:
                pairs.add((cell.cell_id, observable))
        elif isinstance(comparison, Diverges):
            for observable in comparison.diverging:
                pairs.add((cell.cell_id, observable))
            for observable in comparison.others:
                pairs.add((cell.cell_id, observable))
    return pairs


def _diverging_pairs(cells: Sequence) -> set:
    pairs = set()
    for cell in cells:
        if isinstance(cell.comparison, Diverges):
            for observable in cell.comparison.diverging:
                pairs.add((cell.cell_id, observable))
    return pairs


def _diagnostic_pairs() -> set:
    return {(m, f) for m in DIAGNOSTIC_METRICS for f in DIAGNOSTIC_FIXTURES}


# ---------------------------------------------------------------------------
# Validation: every completeness/uniqueness/consistency rule the
# artifact set must satisfy. Returns a list of human-readable
# violation strings (empty list == valid); never raises on a bad
# artifact set -- callers decide what "invalid" means for their flow
# (``build_report`` refuses to compute reducers when this is non-empty).
# ---------------------------------------------------------------------------


def validate_samples(samples: Sequence[Sample], expected_source_sha: str,
                      cells: Sequence = CELLS) -> List[str]:
    errors: List[str] = []

    measurement = [s for s in samples if s.kind == "measurement"]
    unittest_samples = [s for s in samples if s.kind == "unittest"]

    # --- multiplicities: 12 measurement samples, 4 unittest samples ---
    if len(measurement) != 12:
        errors.append(
            "expected 12 measurement samples (3 invocations x 4 gating "
            "runners), got {}".format(len(measurement))
        )
    if len(unittest_samples) != 4:
        errors.append(
            "expected 4 unittest-timing samples (1 per gating runner), "
            "got {}".format(len(unittest_samples))
        )

    # --- (runner, invocation) uniqueness, within each kind ---
    for kind, group in (("measurement", measurement), ("unittest", unittest_samples)):
        seen: Dict[Tuple[str, int], int] = Counter()
        for s in group:
            seen[(s.runner, s.invocation)] += 1
        for key, count in seen.items():
            if count > 1:
                errors.append(
                    "duplicate (runner, invocation) identity {!r} "
                    "appears {} times among {} samples".format(key, count, kind)
                )

    # --- runner-set shape: exactly the 4 gating runners, 3 measurement
    #     invocations {1,2,3} and exactly 1 unittest invocation each ---
    measurement_counts = Counter(s.runner for s in measurement)
    for runner in GATING_RUNNERS:
        count = measurement_counts.get(runner, 0)
        if count != 3:
            errors.append(
                "runner {!r}: expected 3 measurement samples, got {}".format(runner, count)
            )
    for runner in measurement_counts:
        if runner not in GATING_RUNNERS:
            errors.append(
                "measurement sample runner {!r} is not one of the gating "
                "runners {!r}".format(runner, GATING_RUNNERS)
            )
    for runner_group, group in (("measurement", measurement),):
        invocations_by_runner: Dict[str, set] = defaultdict(set)
        for s in group:
            invocations_by_runner[s.runner].add(s.invocation)
        for runner, invocations in invocations_by_runner.items():
            if runner in GATING_RUNNERS and invocations != {1, 2, 3}:
                errors.append(
                    "runner {!r}: measurement invocations {!r} != "
                    "{{1, 2, 3}}".format(runner, sorted(invocations))
                )

    unittest_counts = Counter(s.runner for s in unittest_samples)
    for runner in GATING_RUNNERS:
        count = unittest_counts.get(runner, 0)
        if count != 1:
            errors.append(
                "runner {!r}: expected 1 unittest-timing sample, got {}".format(runner, count)
            )
    for runner in unittest_counts:
        if runner not in GATING_RUNNERS:
            errors.append(
                "unittest-timing sample runner {!r} is not one of the "
                "gating runners {!r}".format(runner, GATING_RUNNERS)
            )

    # --- per-sample content: zero error records ---
    for s in samples:
        error_records = [rec for rec in s.records if "error" in rec]
        if error_records:
            errors.append(
                "{} sample (runner={!r}, invocation={}): {} 'error' "
                "record(s) present -- {}".format(
                    s.kind, s.runner, s.invocation, len(error_records),
                    "; ".join(str(rec.get("error")) for rec in error_records[:3]),
                )
            )

    # --- measurement samples: exactly one summary line; source_sha ==
    #     expected; runner metadata present and self-consistent across
    #     a runner's own 3 invocations ---
    runner_metadata: Dict[str, dict] = {}
    for s in measurement:
        summaries = [rec for rec in s.records if "module_total_seconds" in rec]
        if len(summaries) != 1:
            errors.append(
                "measurement sample (runner={!r}, invocation={}): "
                "expected exactly 1 summary line ('module_total_seconds' "
                "present), got {}".format(s.runner, s.invocation, len(summaries))
            )
            continue
        summary = summaries[0]
        sha = summary.get("source_sha")
        if sha != expected_source_sha:
            errors.append(
                "measurement sample (runner={!r}, invocation={}): "
                "source_sha {!r} != expected {!r}".format(
                    s.runner, s.invocation, sha, expected_source_sha
                )
            )
        meta = summary.get("runner")
        if not isinstance(meta, dict) or not meta.get("python"):
            errors.append(
                "measurement sample (runner={!r}, invocation={}): missing "
                "or malformed 'runner' metadata block".format(s.runner, s.invocation)
            )
            continue
        if not str(meta["python"]).startswith(s.runner):
            errors.append(
                "measurement sample (runner={!r}, invocation={}): summary "
                "runner.python {!r} does not start with the declared "
                "runner label {!r}".format(s.runner, s.invocation, meta["python"], s.runner)
            )
        prior = runner_metadata.get(s.runner)
        if prior is not None and prior != meta:
            errors.append(
                "runner {!r}: inconsistent runner metadata across its own "
                "invocations -- {!r} vs {!r}".format(s.runner, prior, meta)
            )
        runner_metadata[s.runner] = meta

    # --- required cell/observable/diagnostic multiplicities, per
    #     measurement sample ---
    expected_cell_ids = {c.cell_id for c in cells}
    expected_obs_pairs = _observable_pairs(cells)
    expected_diag_pairs = _diagnostic_pairs()

    for s in measurement:
        timing_cells = {
            rec["cell"] for rec in s.records
            if "cell" in rec and "seconds" in rec and "observable" not in rec and "error" not in rec
        }
        missing_cells = expected_cell_ids - timing_cells
        extra_cells = timing_cells - expected_cell_ids
        if missing_cells:
            errors.append(
                "measurement sample (runner={!r}, invocation={}): missing "
                "timing line(s) for cell(s) {!r}".format(
                    s.runner, s.invocation, sorted(missing_cells)
                )
            )
        if extra_cells:
            errors.append(
                "measurement sample (runner={!r}, invocation={}): "
                "unexpected timing line(s) for cell(s) {!r}".format(
                    s.runner, s.invocation, sorted(extra_cells)
                )
            )

        obs_pairs = {
            (rec["cell"], rec["observable"]) for rec in s.records
            if "observable" in rec and "error" not in rec
        }
        missing_obs = expected_obs_pairs - obs_pairs
        extra_obs = obs_pairs - expected_obs_pairs
        if missing_obs:
            errors.append(
                "measurement sample (runner={!r}, invocation={}): missing "
                "observable residual(s) {!r}".format(
                    s.runner, s.invocation, sorted(missing_obs)
                )
            )
        if extra_obs:
            errors.append(
                "measurement sample (runner={!r}, invocation={}): "
                "unexpected observable residual(s) {!r}".format(
                    s.runner, s.invocation, sorted(extra_obs)
                )
            )

        diag_records = [
            (rec["diagnostic"], rec["fixture"]) for rec in s.records
            if "diagnostic" in rec and "error" not in rec
        ]
        diag_counts = Counter(diag_records)
        for pair, count in diag_counts.items():
            if count > 1:
                errors.append(
                    "measurement sample (runner={!r}, invocation={}): "
                    "duplicate diagnostic record for (metric, fixture) "
                    "{!r} appears {} times".format(
                        s.runner, s.invocation, pair, count
                    )
                )
        diag_pairs = set(diag_counts)
        missing_diag = expected_diag_pairs - diag_pairs
        extra_diag = diag_pairs - expected_diag_pairs
        if missing_diag:
            errors.append(
                "measurement sample (runner={!r}, invocation={}): missing "
                "diagnostic checkpoint(s) (metric, fixture) {!r}".format(
                    s.runner, s.invocation, sorted(missing_diag)
                )
            )
        if extra_diag:
            errors.append(
                "measurement sample (runner={!r}, invocation={}): "
                "unexpected diagnostic checkpoint(s) (metric, fixture) "
                "{!r}".format(
                    s.runner, s.invocation, sorted(extra_diag)
                )
            )

    # --- unittest samples: exactly one 'unittest_module_process_seconds' ---
    for s in unittest_samples:
        values = [rec for rec in s.records if "unittest_module_process_seconds" in rec]
        if len(values) != 1:
            errors.append(
                "unittest sample (runner={!r}): expected exactly 1 "
                "'unittest_module_process_seconds' record, got {}".format(
                    s.runner, len(values)
                )
            )

    return errors


# ---------------------------------------------------------------------------
# Reducers -- the aggregation rules a freeze decision is computed with.
# Each raises ValueError if the requested (cell, observable)/checkpoint
# has zero matching records -- silently returning e.g. ``-inf`` would
# let a missing quantity masquerade as "passed the bound", which is
# exactly the failure mode this module exists to prevent. Callers
# normally run
# ``validate_samples`` first and only call these when it returned no
# errors.
# ---------------------------------------------------------------------------


def max_residual(samples: Sequence[Sample], cell_id: str, observable: str) -> float:
    """Upper-bound reducer (atol / diagnostic-ceiling checks): per-
    runner MAX over its invocations, then cross-runner MAX. (The two-
    stage description is the definitional one; a flat max over every
    matching record is numerically identical to max-of-maxes and is
    what this implementation computes.)
    """

    values = [
        rec["residual"] for s in samples if s.kind == "measurement"
        for rec in s.records
        if rec.get("cell") == cell_id and rec.get("observable") == observable and "error" not in rec
    ]
    if not values:
        raise ValueError(
            "max_residual: no matching records for cell={!r} "
            "observable={!r}".format(cell_id, observable)
        )
    return max(values)


def min_residual(samples: Sequence[Sample], cell_id: str, observable: str) -> float:
    """Diverges-classification reducer: per-runner MIN over
    invocations, then cross-runner MIN -- the strict-lower assertion
    (``residual > ceiling``) must hold on EVERY ordinary execution, so
    a single low sample forbids classifying the cell as Diverges.
    """

    per_runner: Dict[str, List[float]] = defaultdict(list)
    for s in samples:
        if s.kind != "measurement":
            continue
        for rec in s.records:
            if rec.get("cell") == cell_id and rec.get("observable") == observable and "error" not in rec:
                per_runner[s.runner].append(rec["residual"])
    if not per_runner:
        raise ValueError(
            "min_residual: no matching records for cell={!r} "
            "observable={!r}".format(cell_id, observable)
        )
    runner_mins = [min(vals) for vals in per_runner.values()]
    return min(runner_mins)


def max_cell_seconds(samples: Sequence[Sample]) -> Dict[str, float]:
    """MAX ``seconds`` per cell, over every measurement sample (per-
    runner MAX then cross-runner MAX, as ``max_residual``)."""

    per_cell: Dict[str, List[float]] = defaultdict(list)
    for s in samples:
        if s.kind != "measurement":
            continue
        for rec in s.records:
            if "cell" in rec and "seconds" in rec and "observable" not in rec and "error" not in rec:
                per_cell[rec["cell"]].append(rec["seconds"])
    return {cell: max(vals) for cell, vals in per_cell.items()}


def max_module_total_seconds(samples: Sequence[Sample]) -> float:
    values = [
        rec["module_total_seconds"] for s in samples if s.kind == "measurement"
        for rec in s.records if "module_total_seconds" in rec
    ]
    if not values:
        raise ValueError("max_module_total_seconds: no summary lines found")
    return max(values)


def unittest_gate_seconds(samples: Sequence[Sample]) -> float:
    """The 120s freeze-time budget's own reducer: MAX over EXACTLY the
    four ``unittest_module_process_seconds`` samples (one per gating
    runner). Raises if the count is not exactly 4 -- the budget
    decision must never silently run on a partial sample set.
    """

    values = [
        rec["unittest_module_process_seconds"]
        for s in samples if s.kind == "unittest"
        for rec in s.records if "unittest_module_process_seconds" in rec
    ]
    if len(values) != 4:
        raise ValueError(
            "unittest_gate_seconds: expected exactly 4 "
            "unittest_module_process_seconds samples (1 per gating "
            "runner), got {}".format(len(values))
        )
    return max(values)


def paired_amplification_ratios(samples: Sequence[Sample], metric: str,
                                 benign_fixture: str = "benign",
                                 fc_fixture: str = "fc") -> Dict[str, float]:
    """The conditioning amplification reducer, per Global Constraints:
    per runner, ``MIN(FC invocations) / MAX(benign invocations)`` with
    ``max(denominator, 1e-15)``. Returns one ratio per runner that
    reported BOTH fixtures for this ``metric`` (a ``DIAGNOSTIC_METRICS``
    name, e.g. ``"assertion2"``); a runner reporting only one side is
    silently omitted here (``assert_amplification_holds`` below is what
    enforces "every runner must clear the bar" and will catch a runner
    with no ratio at all as a missing-runner error).
    """

    fc_vals: Dict[str, List[float]] = defaultdict(list)
    benign_vals: Dict[str, List[float]] = defaultdict(list)
    for s in samples:
        if s.kind != "measurement":
            continue
        for rec in s.records:
            if rec.get("diagnostic") != metric or "error" in rec:
                continue
            if rec.get("fixture") == fc_fixture:
                fc_vals[s.runner].append(rec["residual"])
            elif rec.get("fixture") == benign_fixture:
                benign_vals[s.runner].append(rec["residual"])

    ratios = {}
    for runner in fc_vals:
        if runner not in benign_vals:
            continue
        fc_min = min(fc_vals[runner])
        benign_max = max(benign_vals[runner])
        ratios[runner] = fc_min / max(benign_max, 1e-15)
    return ratios


def assert_amplification_holds(ratios: Dict[str, float],
                                threshold: float = AMPLIFICATION_THRESHOLD,
                                expected_runners: Sequence[str] = GATING_RUNNERS) -> None:
    """Global Constraints: ">=10x threshold must hold on EVERY
    runner". Raises ``ValueError`` naming every runner that either
    never reported a ratio or fell below the threshold.
    """

    problems = []
    for runner in expected_runners:
        if runner not in ratios:
            problems.append("{!r}: no amplification ratio computed".format(runner))
        elif ratios[runner] < threshold:
            problems.append("{!r}: {:.3g}x < {:.3g}x".format(runner, ratios[runner], threshold))
    if problems:
        raise ValueError(
            "assert_amplification_holds: threshold not met -- " + "; ".join(problems)
        )


# ---------------------------------------------------------------------------
# The report: the SOLE input to any freeze decision.
# ---------------------------------------------------------------------------


def build_report(samples: Sequence[Sample], expected_source_sha: str,
                  cells: Sequence = CELLS) -> str:
    """Validate ``samples`` and, if valid, compute every reducer this
    module defines for every comparison cell/observable in ``cells``.
    Returns a human-readable report string. If validation fails, the
    report states so and lists every violation -- it does NOT attempt
    to compute reducers over an incomplete/inconsistent sample set.
    """

    lines: List[str] = []
    errors = validate_samples(samples, expected_source_sha, cells)
    if errors:
        lines.append("VALIDATION FAILED ({} issue(s)):".format(len(errors)))
        for e in errors:
            lines.append("  - {}".format(e))
        return "\n".join(lines)

    lines.append("VALIDATION OK (12 measurement + 4 unittest samples, source_sha {})".format(
        expected_source_sha
    ))
    lines.append("")
    lines.append("Per-cell/observable MAX residual (candidate/frozen atol bound source):")
    for cell in cells:
        comparison = cell.comparison
        observables: List[str] = []
        if isinstance(comparison, Equiv):
            observables = list(comparison.observables)
        elif isinstance(comparison, Diverges):
            observables = list(comparison.diverging) + list(comparison.others)
        for observable in observables:
            value = max_residual(samples, cell.cell_id, observable)
            lines.append("  {} / {}: {:.6e}".format(cell.cell_id, observable, value))

    lines.append("")
    lines.append("120s freeze-time budget (MAX over the 4 unittest_module_process_seconds samples):")
    lines.append("  {:.3f}s".format(unittest_gate_seconds(samples)))

    lines.append("")
    lines.append("Conditioning amplification (assertion2), per runner:")
    for runner, ratio in sorted(
        paired_amplification_ratios(samples, metric="assertion2").items()
    ):
        lines.append("  {}: {:.3g}x".format(runner, ratio))

    return "\n".join(lines)


def main() -> int:  # pragma: no cover -- thin CLI wrapper, exercised manually
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", help="directory containing the downloaded calib-*/unittest-* artifacts")
    parser.add_argument("--source-sha", required=True, help="the expected (measured) source commit S")
    args = parser.parse_args()

    samples = load_samples_from_dir(args.directory)
    print(build_report(samples, args.source_sha))
    return 0


if __name__ == "__main__":  # pragma: no cover
    import sys

    sys.exit(main())
