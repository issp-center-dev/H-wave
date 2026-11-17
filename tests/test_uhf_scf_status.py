"""``scf_status.dat``: the convergence status that UHFr and UHFk record.

Format version 1: one ``"{key} = {value}\\n"`` line per key, in the order of
:data:`KEYS`. ``[file.output] scf_status`` names the file (default
``"scf_status.dat"``); ``""`` disables it. With ``IterationMax = 0`` the
file is the only result written (a configured ``print_check`` file exists
and is empty), and one WARNING says so. The key is matched exactly; a
value that is not a string is refused before the SCF starts.
"""
import copy
import logging
import math
import os
import re
import shutil
import tempfile
import unittest

import tomli

from hwave import qlms

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

KEYS = ("format_version", "solver", "converged", "residual", "eps", "iterations",
        "iteration_max", "residual_definition", "note")

_DEFINITION = {"UHFk": "two_frobenius_over_size:uhfk_all_r",
               "UHFr": "two_frobenius_over_size:uhfr_full"}

_CASE = {"UHFk": "tests/uhfk/CoulombIntra", "UHFr": "tests/uhfr/coulombintra"}

_ZERO_WARNING = "no SCF step was performed (IterationMax = 0)"

def _input(solver, out, iteration_max=None, random_start=False, print_check=None,
           output=True):
    """The committed input of ``solver``'s case, with absolute input paths
    and the output in ``out``."""
    case = os.path.join(_REPO, _CASE[solver])
    with open(os.path.join(case, "input.toml"), "rb") as f:
        inp = tomli.load(f)
    inp["log"]["print_level"] = 0
    if print_check is not None:
        inp["log"]["print_check"] = print_check
    finput = inp["file"]["input"]
    finput["path_to_input"] = case
    if "interaction" in finput and solver == "UHFk":
        finput["interaction"]["path_to_input"] = case
    if random_start:
        finput.pop("initial", None)
        finput.pop("Initial", None)
        if solver == "UHFk":
            finput["initial_mode"] = "random"
    if iteration_max is not None:
        inp["mode"]["param"]["IterationMax"] = iteration_max
    if output:
        inp["file"]["output"]["path_to_output"] = out
    else:
        del inp["file"]["output"]
    return inp


def read_status(path):
    """``{key: value string}`` of a status file, and its raw text."""
    with open(path, "rb") as f:
        raw = f.read().decode()
    table = {}
    for line in raw.splitlines():
        key, _, value = line.partition(" = ")
        table[key] = value
    return table, raw


def check_invariants(test, table):
    eps = float(table["eps"])
    iterations = int(table["iterations"])
    iteration_max = int(table["iteration_max"])
    residual = float(table["residual"])
    test.assertTrue(math.isfinite(eps) and eps > 0)
    test.assertIn(table["converged"], ("true", "false"))
    test.assertTrue(0 <= iterations <= iteration_max)
    if iterations > 0:
        if table["residual"] == "nan":
            test.assertEqual(table["note"],
                             "residual not finite after step {}".format(iterations))
            test.assertEqual(table["converged"], "false")
        else:
            test.assertTrue(math.isfinite(residual) and residual >= 0)
    if table["converged"] == "true":
        test.assertGreater(iterations, 0)
        test.assertLess(residual, eps)
    if iterations == 0:
        test.assertEqual(table["converged"], "false")
        test.assertEqual(table["residual"], "nan")


class _Base(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self._level = logging.getLogger().level
        logging.getLogger().setLevel(logging.ERROR)

    def tearDown(self):
        logging.getLogger().setLevel(self._level)
        shutil.rmtree(self.d, ignore_errors=True)

    def _out(self, name):
        return os.path.join(self.d, name)

    def _run(self, inp):
        qlms.run(input_dict=copy.deepcopy(inp))

    def _status(self, out, name="scf_status.dat"):
        table, raw = read_status(os.path.join(out, name))
        check_invariants(self, table)
        return table, raw


class TestConvergence(_Base):

    def test_uhfk_converged(self):
        out = self._out("k")
        self._run(_input("UHFk", out))
        t, _ = self._status(out)
        self.assertEqual(t["format_version"], "1")
        self.assertEqual(t["solver"], "UHFk")
        self.assertEqual(t["converged"], "true")
        self.assertEqual(t["iterations"], "1")
        self.assertEqual(t["iteration_max"], "1000")
        self.assertLess(float(t["residual"]), 1e-12)
        self.assertEqual(t["eps"], repr(1e-12))
        self.assertEqual(t["residual_definition"], _DEFINITION["UHFk"])
        self.assertEqual(t["note"], "")

    def test_uhfk_not_converged(self):
        out = self._out("k")
        self._run(_input("UHFk", out, iteration_max=5, random_start=True))
        t, _ = self._status(out)
        self.assertEqual(t["converged"], "false")
        self.assertEqual(t["iterations"], "5")
        self.assertEqual(t["iteration_max"], "5")
        self.assertGreater(float(t["residual"]), float(t["eps"]))
        # the other results are still written
        self.assertTrue(os.path.exists(os.path.join(out, "energy.dat")))

    def test_uhfr_converged(self):
        out = self._out("r")
        self._run(_input("UHFr", out))
        t, _ = self._status(out)
        self.assertEqual(t["solver"], "UHFr")
        self.assertEqual(t["converged"], "true")
        self.assertGreater(int(t["iterations"]), 0)
        self.assertLess(float(t["residual"]), float(t["eps"]))
        self.assertEqual(t["eps"], repr(1e-8))
        self.assertEqual(t["residual_definition"], _DEFINITION["UHFr"])
        self.assertEqual(t["note"], "")

    def test_uhfr_not_converged(self):
        out = self._out("r")
        self._run(_input("UHFr", out, iteration_max=5, random_start=True))
        t, _ = self._status(out)
        self.assertEqual(t["converged"], "false")
        self.assertEqual(t["iterations"], "5")
        self.assertGreater(float(t["residual"]), float(t["eps"]))
        self.assertTrue(os.path.exists(os.path.join(out, "energy.dat")))

    def test_attributes_after_steps(self):
        """After a 5-step run the solver's attributes match the last line of
        the ``print_check`` trace (Rest to 8 significant digits)."""
        from hwave.solver import uhfk, uhfr
        for solver, cls in (("UHFk", uhfk.UHFk), ("UHFr", uhfr.UHFr)):
            with self.subTest(solver=solver):
                seen = []
                orig = cls.save_results

                def spy(this, *args, **kw):
                    seen.append((this.scf_converged, this.scf_iterations,
                                 this.scf_residual))
                    return orig(this, *args, **kw)

                out = self._out(solver + "_attr")
                cls.save_results = spy
                try:
                    self._run(_input(solver, out, iteration_max=5, random_start=True,
                                     print_check="check.dat"))
                finally:
                    cls.save_results = orig
                (converged, iterations, residual), = seen
                self.assertIs(converged, False)
                self.assertEqual(iterations, 5)
                with open(os.path.join(out, "check.dat")) as f:
                    last = f.read().splitlines()[-1].split(",")
                self.assertEqual(int(last[0]), 4)
                self.assertEqual("{:.8g}".format(residual), last[1].strip())
                t, _ = self._status(out)
                self.assertEqual(t["residual"], repr(residual))


class TestIterationMaxZero(_Base):

    def _zero_run(self, solver, **kw):
        out = self._out(solver)
        inp = _input(solver, out, iteration_max=0, **kw)
        with self.assertLogs("qlms", level="WARNING") as logs:
            self._run(inp)
        warnings = [r for r in logs.records if _ZERO_WARNING in r.getMessage()]
        self.assertEqual(len(warnings), 1, logs.output)
        return out, inp, warnings[0].getMessage()

    def _assert_zero_status(self, t, solver):
        self.assertEqual(t["solver"], solver)
        self.assertEqual(t["converged"], "false")
        self.assertEqual(t["iterations"], "0")
        self.assertEqual(t["iteration_max"], "0")
        self.assertEqual(t["residual"], "nan")
        self.assertEqual(t["note"], "no SCF step performed")

    def test_iteration_max_zero_via_run(self):
        for solver in ("UHFk", "UHFr"):
            with self.subTest(solver=solver):
                out, inp, msg = self._zero_run(solver)
                # energy, eigen and green are all configured
                for key in ("energy", "eigen", "green"):
                    self.assertIn(key, inp["file"]["output"])
                self.assertEqual(sorted(os.listdir(out)), ["scf_status.dat"])
                self.assertIn("no result file other than scf_status.dat is written", msg)
                t, _ = self._status(out)
                self._assert_zero_status(t, solver)

    def test_iteration_max_zero_with_print_check(self):
        for solver in ("UHFk", "UHFr"):
            with self.subTest(solver=solver):
                out, _, _ = self._zero_run(solver, print_check="check.dat")
                self.assertEqual(sorted(os.listdir(out)), ["check.dat", "scf_status.dat"])
                self.assertEqual(os.path.getsize(os.path.join(out, "check.dat")), 0)
                self._assert_zero_status(self._status(out)[0], solver)

    def test_iteration_max_zero_opted_out(self):
        for solver in ("UHFk", "UHFr"):
            with self.subTest(solver=solver):
                out = self._out(solver)
                inp = _input(solver, out, iteration_max=0)
                inp["file"]["output"]["scf_status"] = ""
                with self.assertLogs("qlms", level="WARNING") as logs:
                    self._run(inp)
                msgs = [r.getMessage() for r in logs.records
                        if _ZERO_WARNING in r.getMessage()]
                self.assertEqual(len(msgs), 1, logs.output)
                self.assertIn("no result file is written", msgs[0])
                self.assertEqual(os.listdir(out), [])

    def test_iteration_max_zero_direct(self):
        """The solver API without ``qlms.run``, built as qlms.run builds it."""
        from hwave.qlmsio.read_input import QLMSInput
        from requests.structures import CaseInsensitiveDict
        from hwave.qlmsio.read_input_k import QLMSkInput
        from hwave.solver.base import ParamDict
        from hwave.solver.uhfk import UHFk
        from hwave.solver.uhfr import UHFr
        for solver in ("UHFk", "UHFr"):
            with self.subTest(solver=solver):
                out = self._out(solver + "_direct")
                os.makedirs(out)
                inp = _input(solver, out, iteration_max=0)
                info_log = dict(inp["log"], print_step=1)
                info_mode = dict(inp["mode"])
                info_mode["param"] = ParamDict(info_mode["param"])
                finput = inp["file"]["input"]
                if solver == "UHFk":
                    read_io = QLMSkInput(finput)
                    cls = UHFk
                else:
                    files = CaseInsensitiveDict(
                        {k: os.path.join(finput["path_to_input"], v)
                         for k, v in finput["interaction"].items()})
                    for key in ("Initial", "OneBodyG"):
                        if key in finput:
                            files[key] = os.path.join(finput["path_to_input"], finput[key])
                    read_io = QLMSInput(files)
                    cls = UHFr
                green_info = read_io.get_param("green")
                sol = cls(read_io.get_param("ham"), info_log, info_mode)
                sol.solve(green_info, out)
                self.assertEqual((sol.scf_converged, sol.scf_iterations), (False, 0))
                self.assertTrue(math.isnan(sol.scf_residual))
                sol.save_results(inp["file"]["output"], green_info)
                self.assertEqual(sorted(os.listdir(out)), ["scf_status.dat"])
                self._assert_zero_status(self._status(out)[0], solver)


class TestFileFormat(_Base):

    def test_line_format_and_order(self):
        for solver in ("UHFk", "UHFr"):
            with self.subTest(solver=solver):
                out = self._out(solver)
                self._run(_input(solver, out))
                _, raw = self._status(out)
                lines = raw.splitlines(keepends=True)
                self.assertEqual([ln.split(" = ")[0] for ln in lines], list(KEYS))
                for line in lines:
                    self.assertRegex(line, r"^[a-z_]+ = .*\n$")
                self.assertEqual(lines[-1], "note = \n")
                self.assertNotIn("np.", raw)

    def test_default_without_file_output(self):
        """No ``[file.output]`` table: qlms.run writes into ``output``
        relative to the working directory."""
        here = os.getcwd()
        work = self._out("cwd")
        os.makedirs(work)
        try:
            os.chdir(work)
            self._run(_input("UHFk", None, output=False))
        finally:
            os.chdir(here)
        t, _ = self._status(os.path.join(work, "output"))
        self.assertEqual(t["converged"], "true")

    def test_custom_name(self):
        out = self._out("k")
        inp = _input("UHFk", out)
        inp["file"]["output"]["scf_status"] = "status.txt"
        self._run(inp)
        self.assertFalse(os.path.exists(os.path.join(out, "scf_status.dat")))
        self.assertEqual(self._status(out, "status.txt")[0]["converged"], "true")

    def test_disabled(self):
        for solver in ("UHFk", "UHFr"):
            with self.subTest(solver=solver):
                out = self._out(solver)
                inp = _input(solver, out)
                inp["file"]["output"]["scf_status"] = ""
                self._run(inp)
                self.assertFalse(os.path.exists(os.path.join(out, "scf_status.dat")))
                self.assertTrue(os.path.exists(os.path.join(out, "energy.dat")))


class TestEdgeCases(_Base):
    """A non-finite residual, a negative ``IterationMax``, a misspelled or
    non-string ``scf_status`` and the warnings of a zero-step run."""

    def _inject(self, cls, value, from_step):
        """Make ``cls._calc_phys`` report ``value`` as the residual from the
        ``from_step``-th call (1-based) on; returns the restore callable."""
        orig = cls._calc_phys
        calls = []

        def patched(this, *args, **kw):
            result = orig(this, *args, **kw)
            calls.append(1)
            if len(calls) >= from_step:
                this.physics["Rest"] = value
            return result

        cls._calc_phys = patched
        return lambda: setattr(cls, "_calc_phys", orig)

    def test_non_finite_residual(self):
        from hwave.solver import uhfk, uhfr
        for solver, cls, value in (("UHFk", uhfk.UHFk, float("nan")),
                                   ("UHFr", uhfr.UHFr, float("nan")),
                                   ("UHFk", uhfk.UHFk, float("inf"))):
            with self.subTest(solver=solver, value=value):
                out = self._out(solver + repr(value))
                restore = self._inject(cls, value, from_step=3)
                try:
                    self._run(_input(solver, out, iteration_max=5, random_start=True))
                finally:
                    restore()
                t, _ = self._status(out)
                self.assertEqual(t["iterations"], "5")
                self.assertEqual(t["residual"], "nan")
                self.assertEqual(t["converged"], "false")
                self.assertEqual(t["note"], "residual not finite after step 5")
                # the other results are written as after any capped run
                self.assertTrue(os.path.exists(os.path.join(out, "energy.dat")))

    def test_negative_iteration_max(self):
        """Reachable only with relaxed parameter checks."""
        for solver in ("UHFk", "UHFr"):
            with self.subTest(solver=solver):
                out = self._out(solver)
                inp = _input(solver, out, iteration_max=-3)
                inp["mode"]["param"]["trustme_interaction_range"] = True
                with self.assertLogs("qlms", level="WARNING") as logs:
                    self._run(inp)
                msgs = [r.getMessage() for r in logs.records
                        if "no SCF step was performed" in r.getMessage()]
                self.assertEqual(len(msgs), 1, logs.output)
                self.assertIn("(IterationMax = -3)", msgs[0])
                self.assertEqual(sorted(os.listdir(out)), ["scf_status.dat"])
                t, _ = self._status(out)
                self.assertEqual(t["iteration_max"], "0")
                self.assertEqual(t["iterations"], "0")
                self.assertEqual(t["residual"], "nan")
                self.assertEqual(t["note"], "IterationMax -3 treated as 0")

    def test_key_matched_exactly(self):
        out = self._out("k")
        inp = _input("UHFk", out)
        inp["file"]["output"]["SCF_Status"] = "other.dat"
        self._run(inp)
        self.assertFalse(os.path.exists(os.path.join(out, "other.dat")))
        self.assertEqual(self._status(out)[0]["converged"], "true")

    def test_non_string_value_refused_before_solve(self):
        from hwave.solver.base import scf_status_file_name
        for solver in ("UHFk", "UHFr"):
            for value in (False, 0, None):
                with self.subTest(solver=solver, value=value):
                    out = self._out(solver + repr(value))
                    inp = _input(solver, out, print_check="check.dat")
                    inp["file"]["output"]["scf_status"] = value
                    with self.assertRaisesRegex(ValueError, "scf_status"):
                        self._run(inp)
                    # solve creates the print_check file first: it never ran
                    self.assertEqual(os.listdir(out), [])
        with self.assertRaises(ValueError):
            scf_status_file_name({"scf_status": 1})
        self.assertEqual(scf_status_file_name({}), "scf_status.dat")
        self.assertEqual(scf_status_file_name({"scf_status": ""}), "")

    def test_zero_steps_single_warning(self):
        """After zero steps the one zero-step WARNING is the only warning
        (UHFr no longer adds "calculation is failed: rest=1.0")."""
        for solver in ("UHFk", "UHFr"):
            with self.subTest(solver=solver):
                out = self._out(solver)
                inp = _input(solver, out, iteration_max=0)
                # the committed UHFr case starts from a file and does not
                # read RndSeed, which qlms.run reports as an unused key
                del inp["mode"]["param"]["RndSeed"]
                with self.assertLogs("qlms", level="WARNING") as logs:
                    self._run(inp)
                self.assertEqual(len(logs.records), 1, logs.output)
                self.assertIn(_ZERO_WARNING, logs.records[0].getMessage())


class TestInvariants(_Base):
    """``test_invariants`` re-reads every status file of a set of runs (a
    converged and a capped run of each solver, and the zero-step runs) and
    checks the invariants of format version 1 on each. The runs are made
    here, so the test does not depend on the order in which the others
    run."""

    def test_invariants(self):
        files = []
        for solver in ("UHFk", "UHFr"):
            for name, kw in (("conv", {}),
                             ("cap", {"iteration_max": 5, "random_start": True}),
                             ("zero", {"iteration_max": 0})):
                out = self._out(solver + name)
                self._run(_input(solver, out, **kw))
                files.append(os.path.join(out, "scf_status.dat"))
        for path in files:
            with self.subTest(path=os.path.basename(os.path.dirname(path))):
                table, raw = read_status(path)
                self.assertEqual(tuple(table), KEYS)
                check_invariants(self, table)


if __name__ == "__main__":
    unittest.main()
