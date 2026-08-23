"""Unit tests for ``tests/equivalence_measure.py`` -- the read-only
measurement entry point behind the RPA/FLEX equivalence-table
calibration.

This module pins ONE behavior: ``_measure_cell`` must catch
``SystemExit`` alongside ``Exception`` so a ``sys.exit(...)`` reached
through a cell's solver-construction path (e.g.
``src/hwave/solver/rpa.py``'s ``Lattice._init_lattice`` on an
incompatible ``CellShape``/``SubShape`` pairing -- rpa.py:571-573)
cannot silently kill the whole measurement process. Every other
``equivalence_measure`` behavior (timing, residual emission, the
diagnostic pass, ``main()``'s exit-code contract) is exercised end to
end whenever ``python -m tests.equivalence_measure`` itself is run
during a calibration sweep; this module is deliberately narrow.
"""

from __future__ import annotations

import io
import json
import unittest
from contextlib import redirect_stdout

from tests.equivalence_cells import Cell, ExecuteRun, FixtureSpec, SolverProof, Status
from tests.equivalence_measure import _measure_cell


def _broken_subshape_fixture() -> FixtureSpec:
    # CellShape=(4,4,1) is not evenly divisible by SubShape=(3,1,1) --
    # hits rpa.py's Lattice._init_lattice sys.exit(1) site ("SubShape
    # is not compatible with CellShape") at the very start of
    # RPA.__init__/FLEX.__init__, before any interaction file is even
    # touched. tests/equivalence_input/orb1 is a real, self-contained,
    # otherwise-valid committed fixture directory -- only the
    # CellShape/SubShape pairing is deliberately broken (qlmsio's file
    # readers never reach a sys.exit site on bad file content in this
    # codebase; every reachable sys.exit site sits in solver
    # construction, validating already-parsed parameters).
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


def _system_exit_cell(cell_id: str) -> Cell:
    return Cell(
        cell_id=cell_id,
        fixture=_broken_subshape_fixture(),
        resolved_scheme="general",
        expected_spin_mode="spin-free",
        rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
        flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
        comparison=None,
        required_observables=(),
        interaction_class="onsite",
        notes="",
    )


def _healthy_run_cell(cell_id: str) -> Cell:
    # The broken fixture's sibling with a COMPATIBLE SubShape -- a
    # real ExecuteRun/ExecuteRun solve on both sides, same input_dir,
    # so this test stays independent of the rest of the registry's
    # fixture directories.
    fixture = FixtureSpec(
        input_dir="tests/equivalence_input/orb1",
        interactions={"CoulombIntra": "coulombintra.dat"},
        T=2.0,
        mu=0.0,
        filling=None,
        CellShape=(4, 4, 1),
        SubShape=(1, 1, 1),
        Nmat=8,
        extra_params={},
        calc_type="ring",
        requested_scheme="general",
        enable_spin_orbital=False,
        extern=None,
    )
    return Cell(
        cell_id=cell_id,
        fixture=fixture,
        resolved_scheme="general",
        expected_spin_mode="spin-free",
        rpa=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
        flex=SolverProof(status=Status.SUPPORTED, steps=(ExecuteRun(),)),
        comparison=None,
        required_observables=(),
        interaction_class="onsite",
        notes="",
    )


class TestMeasureCellSystemExitEscape(unittest.TestCase):
    """``_measure_cell`` must fail closed on a ``SystemExit`` raised
    while executing a cell's proofs -- an ``{"error": ...}`` JSON line
    with the exit code recorded, ``False`` returned (never an
    uncaught, process-killing exception), and processing must CONTINUE
    to the next cell.
    """

    def test_system_exit_emits_an_error_line_and_returns_false(self):
        buf = io.StringIO()
        with redirect_stdout(buf):
            ok = _measure_cell(_system_exit_cell("fake.measure.systemexit"))

        self.assertFalse(ok)
        lines = [json.loads(line) for line in buf.getvalue().splitlines() if line.strip()]
        self.assertEqual(len(lines), 1)
        record = lines[0]
        self.assertIn("error", record)
        self.assertEqual(record["cell"], "fake.measure.systemexit")
        self.assertEqual(record["phase"], "cell")
        self.assertEqual(record["exit_code"], 1)
        self.assertIn("seconds", record)

    def test_process_continues_to_the_next_cell_after_a_system_exit(self):
        buf = io.StringIO()
        with redirect_stdout(buf):
            ok_broken = _measure_cell(_system_exit_cell("fake.measure.systemexit"))
            ok_healthy = _measure_cell(_healthy_run_cell("fake.measure.followup"))

        self.assertFalse(ok_broken)
        self.assertTrue(ok_healthy)

        lines = [json.loads(line) for line in buf.getvalue().splitlines() if line.strip()]
        self.assertEqual(len(lines), 2)
        self.assertIn("error", lines[0])
        self.assertEqual(lines[0]["cell"], "fake.measure.systemexit")
        self.assertNotIn("error", lines[1])
        self.assertEqual(lines[1]["cell"], "fake.measure.followup")


if __name__ == "__main__":
    unittest.main()
