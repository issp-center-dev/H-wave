"""The read-only measurement entry point for the RPA/FLEX equivalence
table. It produces the calibration artifacts
``tests/equivalence_freeze_check.py`` aggregates.

``python -m tests.equivalence_measure`` executes and TIMES every cell in
``tests.equivalence_cells.CELLS`` (every proof shape -- rejections,
construction-only rows, and the two multirun exceptions all included,
via the SAME executor (``tests.test_rpa_flex_equivalence_table._run_side``)
``run_cell`` itself dispatches through), then runs the mu/Green
divergence diagnostic (``_diagnostic_residuals``) on both its
fixtures. It prints one JSON object per line (never a JSON array -- so
a partial/interrupted run's already-printed lines stay individually
parseable) to stdout:

  * one ``{"cell": <cell_id>, "seconds": <float>}`` line per cell (every
    cell, regardless of proof shape);
  * one ``{"cell": <cell_id>, "observable": <name>, "residual": <float>}``
    line per (comparison-cell, observable) -- only for cells whose
    ``comparison`` is not ``None`` (``Equiv`` or ``Diverges``);
  * one ``{"diagnostic": 1..5, "fixture": "benign"|"fc", "residual":
    <float>}`` line per diagnostic checkpoint per fixture (10 lines
    total: 5 checkpoints x 2 fixtures);
  * a final ``{"module_total_seconds": <float>, "source_sha": <str>,
    "runner": {"platform": ..., "python": ..., "numpy": ..., "scipy":
    ...}}`` line.

``source_sha`` is ``git rev-parse HEAD`` (run against this repository's
root), falling back to the ``GITHUB_SHA`` environment variable when git
is unavailable (e.g. a tarball checkout) or the call fails.

This module NEVER asserts numerical equivalence -- residuals are
reported, not judged against any ceiling (that judgment is
``tests.equivalence_freeze_check``'s job, downstream). It FAILS CLOSED
on structural problems instead: an unexpected exception while executing
a cell's proofs (an ``ExecuteRun``/``ExecuteConstruct`` solve raising
when it should not, a REJECT proof's expected exception NOT being
raised, a multirun oracle's bitwise/exhaustive check failing), a
comparator's own structural guard firing (shape mismatch, non-finite
values -- ``Comparator.residual``'s ``_check_comparable`` already
raises ``ValueError`` for both), a missing ``green_info`` key or a
non-negligible discarded spin block in the FLEX reduced-scheme chi0q
reduction (``extract_bundle`` raising ``KeyError``/``ValueError``/
``AssertionError``), or a JSON-serialization failure of an emitted
record all produce an
``{"error": <message>, ...}`` line (never silently dropped) and set the
process's exit code to 1 -- but processing CONTINUES to the next cell/
checkpoint first (best-effort telemetry: one broken cell must not hide
every other cell's timing/residual data from the CI artifact).

Reuses ``tests.test_rpa_flex_equivalence_table``'s own builders
(``build_solver``, ``_run_side``, ``extract_bundle``, ``COMPARATORS``)
and the mu/Green divergence diagnostic (``_diagnostic_residuals``,
``_DIAGNOSTIC_BENIGN_FIXTURE``, ``_DIAGNOSTIC_FC_FIXTURE``) -- it does
not reimplement any solver-construction or comparison logic. Every
solver construction happens inside a fresh ``contextlib.ExitStack`` per
cell (mirroring ``run_cell``'s own lifetime contract: a fresh temp
output dir and fresh ``green_info`` per cell, ``tests/rpa/output``
untouched). Writes nothing to the repository -- every artifact is a
line printed to stdout; callers (a shell redirect, the calibration
workflow) own persistence. Normal ``unittest``/``pytest`` runs never
import or invoke this module.
"""

from __future__ import annotations

import json
import os
import platform
import subprocess
import sys
import time
from contextlib import ExitStack

import numpy as np
import scipy

from tests.equivalence_cells import CELLS, COMPARATORS, Diverges, Equiv
from tests.test_rpa_flex_equivalence_table import (
    _DIAGNOSTIC_BENIGN_FIXTURE,
    _DIAGNOSTIC_FC_FIXTURE,
    _cell_chi0q_atol,
    _diagnostic_residuals,
    _run_side,
    build_solver,
    extract_bundle,
)

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _suppress_perf_db_exit_report() -> None:
    """Every ``@do_profile``-wrapped solver method (``src/hwave/solver/
    perf.py``) records into a module-level ``PerfDB`` singleton whose
    ``__del__`` unconditionally prints a "Statistics" table to STDOUT
    when the process tears down -- unrelated project infrastructure,
    deliberately left unmodified here. Left alone, that table's lines
    would land AFTER this module's own final JSON line and break the
    "one JSON object per line, nothing else, ever" stdout contract for
    any downstream consumer that does not know to ignore trailing
    non-JSON text.
    Clearing the singleton's counters here (its own ``__del__`` returns
    immediately when they are empty) suppresses the report without
    touching ``perf.py`` itself.
    """

    import hwave.solver.perf as perf_mod

    perf_mod._perf_db_data._db_count.clear()
    perf_mod._perf_db_data._db_value.clear()


def _source_sha() -> str:
    """``git rev-parse HEAD`` against the repository root; falls back to
    ``$GITHUB_SHA`` when git is unavailable or the call fails (e.g. a
    tarball checkout with no ``.git`` directory)."""

    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=_REPO_ROOT,
            capture_output=True,
            text=True,
            check=True,
        )
        sha = result.stdout.strip()
        if sha:
            return sha
    except Exception:
        pass
    return os.environ.get("GITHUB_SHA", "")


def _emit(record: dict) -> bool:
    """Print one JSON line for ``record``. Returns ``True`` on success;
    on a serialization failure (a structural problem this module must
    fail closed on), prints a minimal fallback error line built ONLY
    from primitive fields (never re-attempts to serialize the
    offending value) and returns ``False``.
    """

    try:
        print(json.dumps(record), flush=True)
        return True
    except (TypeError, ValueError) as exc:
        fallback = {
            "error": "serialization failure: {}".format(exc),
            "phase": "emit",
            "record_keys": sorted(str(k) for k in record.keys()),
        }
        print(json.dumps(fallback), flush=True)
        return False


def _observable_specs(comparison) -> dict:
    if isinstance(comparison, Equiv):
        return dict(comparison.observables)
    if isinstance(comparison, Diverges):
        merged = dict(comparison.diverging)
        merged.update(comparison.others)
        return merged
    raise ValueError(
        "_observable_specs: comparison must be Equiv or Diverges, got "
        "{!r}".format(type(comparison))
    )


def _measure_cell(cell) -> bool:
    """Execute + time every proof ``cell`` records; emit its timing line
    and (if it carries a comparison) its per-observable residual lines.

    Returns ``True`` iff nothing structurally wrong was observed.
    """

    ok = True
    t0 = time.monotonic()
    try:
        with ExitStack() as stack:
            rpa_result = _run_side(cell.fixture, "rpa", cell.rpa, stack, build_solver)
            flex_result = _run_side(cell.fixture, "flex", cell.flex, stack, build_solver)
            elapsed = time.monotonic() - t0
            ok &= _emit({"cell": cell.cell_id, "seconds": float(elapsed)})

            if cell.comparison is not None:
                rpa_obj, rpa_green, _rpa_out = rpa_result
                flex_obj, flex_green, _flex_out = flex_result
                chi0q_atol = _cell_chi0q_atol(cell)
                rpa_bundle = extract_bundle(rpa_obj, rpa_green, "rpa", chi0q_atol)
                flex_bundle = extract_bundle(flex_obj, flex_green, "flex", chi0q_atol)
                for observable, spec in _observable_specs(cell.comparison).items():
                    comparator = COMPARATORS[spec.comparator]
                    a, b = comparator.map(observable, rpa_bundle, flex_bundle)
                    residual = comparator.residual(a, b)
                    ok &= _emit({
                        "cell": cell.cell_id,
                        "observable": observable,
                        "residual": float(residual),
                    })
    except KeyboardInterrupt:
        raise
    except (Exception, SystemExit) as exc:
        # A bare ``except Exception`` would let a ``sys.exit(...)``
        # from an input-reader/solver-construction path (SystemExit
        # derives from BaseException, not Exception) escape uncaught
        # and kill this process -- silently dropping every
        # cell/checkpoint measured after the broken one from the
        # calibration artifact. Explicitly catching SystemExit
        # alongside Exception keeps the "one broken cell must not hide
        # every other cell's data" contract; the exit code is recorded
        # on the emitted error line.
        elapsed = time.monotonic() - t0
        record = {
            "error": str(exc),
            "cell": cell.cell_id,
            "phase": "cell",
            "seconds": float(elapsed),
        }
        if isinstance(exc, SystemExit):
            code = exc.code
            record["exit_code"] = code if isinstance(code, (int, str)) or code is None else str(code)
        _emit(record)
        ok = False
    return ok


def _measure_diagnostics() -> bool:
    """Run the mu/Green divergence diagnostic on both its fixtures;
    emit 5 checkpoint lines per fixture (10 total)."""

    ok = True
    for fixture_name, fixture in (("benign", _DIAGNOSTIC_BENIGN_FIXTURE),
                                   ("fc", _DIAGNOSTIC_FC_FIXTURE)):
        try:
            values = _diagnostic_residuals(fixture)
        except Exception as exc:
            _emit({"error": str(exc), "phase": "diagnostic", "fixture": fixture_name})
            ok = False
            continue
        for checkpoint in range(1, 6):
            residual = values["assertion{}".format(checkpoint)]
            ok &= _emit({
                "diagnostic": checkpoint,
                "fixture": fixture_name,
                "residual": float(residual),
            })
    return ok


def main() -> int:
    module_start = time.monotonic()
    ok = True

    for cell in CELLS:
        ok &= _measure_cell(cell)

    ok &= _measure_diagnostics()

    module_total_seconds = time.monotonic() - module_start
    ok &= _emit({
        "module_total_seconds": float(module_total_seconds),
        "source_sha": _source_sha(),
        "runner": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
    })

    _suppress_perf_db_exit_report()

    if not ok:
        _emit({
            "error": (
                "one or more cells/diagnostic checkpoints raised a "
                "structural problem or failed to serialize; see the "
                "preceding 'error' records"
            ),
        })
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
