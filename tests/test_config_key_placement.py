"""Keys of ``[mode]`` / ``[mode.param]`` that a run would ignore are reported
(GitHub issue #236): a ``[mode]`` key placed in ``[mode.param]``, an unknown
or mis-cased ``[mode]`` key, and a ``[mode.param]`` key no solver reads."""
import copy
import logging
import pathlib
import shutil
import tempfile
import unittest

from hwave.qlms import MODE_TABLE_KEYS, check_mode_tables
from hwave.solver.base import ParamDict


class TestCheckModeTables(unittest.TestCase):

    def test_well_placed_tables_pass(self):
        info = {"mode": "RPA", "calc_scheme": "general", "calc_type": "ring",
                "enable_spin_orbital": False, "flag_fock": True,
                "param": {"T": 0.1, "CellShape": [4, 4, 1]}}
        self.assertEqual(check_mode_tables(info), [])

    def test_mode_key_in_param_is_refused(self):
        for key in ("calc_scheme", "Calc_Type", "enable_spin_orbital", "flag_fock",
                    "mode", "Mode", "MODE"):
            with self.subTest(key=key):
                msgs = check_mode_tables({"mode": "RPA", "param": {key: "x"}})
                self.assertEqual(len(msgs), 1)
                self.assertIn("[mode.param] {!r}".format(key), msgs[0])
                self.assertIn("move it to [mode]", msgs[0])

    def test_nested_param_table_is_refused(self):
        msgs = check_mode_tables({"mode": "RPA", "param": {"param": {"T": 1.0}}})
        self.assertEqual(len(msgs), 1)
        self.assertIn("no nested param table", msgs[0])

    def test_mapping_tables_are_checked_too(self):
        from requests.structures import CaseInsensitiveDict
        info = CaseInsensitiveDict({"mode": "RPA",
                                    "param": CaseInsensitiveDict({"calc_scheme": "general"})})
        msgs = check_mode_tables(info)
        self.assertEqual(len(msgs), 1)
        self.assertIn("move it to [mode]", msgs[0])

    def test_unknown_mode_key_is_refused(self):
        msgs = check_mode_tables({"mode": "RPA", "T": 0.1, "param": {}})
        self.assertEqual(len(msgs), 1)
        self.assertIn("[mode] 'T' is not a [mode] key", msgs[0])
        self.assertIn("[mode.param]", msgs[0])

    def test_mis_cased_mode_key_is_refused(self):
        msgs = check_mode_tables({"mode": "RPA", "Calc_Scheme": "general", "param": {}})
        self.assertEqual(len(msgs), 1)
        self.assertIn("case-sensitive", msgs[0])
        self.assertIn("'calc_scheme'", msgs[0])

    def test_every_listed_key_is_accepted_in_mode(self):
        info = {k: None for k in MODE_TABLE_KEYS}
        info["param"] = {}
        self.assertEqual(check_mode_tables(info), [])


class TestParamDict(unittest.TestCase):

    def test_reads_are_recorded_case_insensitively(self):
        p = ParamDict({"CellShape": [2, 2, 1], "T": 1.0, "Mix": 0.5})
        _ = p["cellshape"]
        _ = p.get("t")
        _ = "MIX" in p
        _ = "absent" in p
        self.assertEqual(p.accessed, {"cellshape", "t", "mix", "absent"})

    def test_iteration_is_not_a_read(self):
        p = ParamDict({"A": 1, "B": 2})
        self.assertEqual(sorted(k for k, _ in p.items()), ["A", "B"])
        self.assertEqual(sorted(p.values()), [1, 2])
        self.assertEqual(sorted(p), ["A", "B"])
        self.assertEqual(p.accessed, set())

    def test_views_follow_later_changes(self):
        p = ParamDict({"A": 1})
        items, values = p.items(), p.values()
        p["B"] = 2
        self.assertEqual(sorted(items), [("A", 1), ("B", 2)])
        self.assertEqual(sorted(values), [1, 2])
        self.assertIn(("a", 1), items)
        self.assertNotIn(("A", 2), items)
        self.assertEqual(len(items), 2)
        self.assertEqual(p.accessed, set())

    def test_copy_method_shares_the_record(self):
        p = ParamDict({"A": 1})
        q = p.copy()
        self.assertIsInstance(q, ParamDict)
        _ = q["a"]
        self.assertEqual(p.accessed, {"a"})

    def test_a_copy_shares_the_record_and_is_not_a_read(self):
        p = ParamDict({"A": 1, "B": 2})
        q = ParamDict(p)
        self.assertEqual(p.accessed, set())
        _ = q["b"]
        self.assertEqual(p.accessed, {"b"})
        self.assertEqual(q["A"], 1)


def _tiny_flex_input(d, extra_mode=None, extra_param=None):
    d = pathlib.Path(d)
    (d / "geom.dat").write_text("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n1\n0.0 0.0 0.0\n")
    (d / "transfer.dat").write_text(
        "transfer\n1\n4\n 1 1 1 1\n"
        "  1  0  0  1  1  -1.0  0.0\n -1  0  0  1  1  -1.0  0.0\n"
        "  0  1  0  1  1  -1.0  0.0\n  0 -1  0  1  1  -1.0  0.0\n")
    (d / "coulombintra.dat").write_text("coulombintra\n1\n1\n 1\n0 0 0 1 1 4.0 0.0\n")
    param = {"T": 4.0, "CellShape": [2, 2, 1], "SubShape": [1, 1, 1], "Nmat": 16,
             "filling": 0.5, "IterationMax": 2, "Mix": 0.5, "EPS": 6}
    param.update(extra_param or {})
    mode = {"mode": "FLEX", "calc_scheme": "reduced", "param": param}
    mode.update(extra_mode or {})
    return {
        "log": {"print_level": 0},
        "mode": mode,
        "file": {"input": {"path_to_input": "",
                           "interaction": {"path_to_input": str(d),
                                           "Geometry": "geom.dat",
                                           "Transfer": "transfer.dat",
                                           "CoulombIntra": "coulombintra.dat"}},
                 "output": {"path_to_output": str(d / "output")}},
    }


def _with_mode(inp, mode, **param):
    inp = copy.deepcopy(inp)
    inp["mode"]["mode"] = mode
    inp["mode"].pop("calc_scheme", None)
    inp["mode"]["param"].update(param)
    return inp


class TestRunReports(unittest.TestCase):
    """End to end through ``hwave.qlms.run``."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _run(self, inp):
        import hwave.qlms as qlms
        return qlms.run(input_dict=copy.deepcopy(inp))

    def test_misspelled_param_key_is_reported(self):
        inp = _tiny_flex_input(self.tmp, extra_param={"Tempreature": 0.1})
        with self.assertLogs("qlms", level=logging.WARNING) as logs:
            self._run(inp)
        hits = [m for m in logs.output if "not read by this FLEX run" in m]
        self.assertEqual(len(hits), 1, logs.output)
        self.assertIn("Tempreature", hits[0])
        # every key the solver does read stays out of the report
        for key in ("CellShape", "Nmat", "filling", "IterationMax"):
            self.assertNotIn(key, hits[0])

    def test_clean_input_reports_nothing(self):
        inp = _tiny_flex_input(self.tmp)
        logger = logging.getLogger("qlms")
        with self.assertLogs(logger, level=logging.WARNING) as logs:
            logger.warning("sentinel")         # assertLogs needs one record
            self._run(inp)
        self.assertFalse([m for m in logs.output if "not read by" in m], logs.output)

    def _reports(self, inp):
        logger = logging.getLogger("qlms")
        with self.assertLogs(logger, level=logging.WARNING) as logs:
            logger.warning("sentinel")         # assertLogs needs one record
            self._run(inp)
        return [m for m in logs.output if "not read by" in m]

    def test_other_modes_report_only_the_misspelled_key(self):
        base = _tiny_flex_input(self.tmp)
        cases = {
            "RPA": _with_mode(base, "RPA"),
            # UHFk checks the interaction range against the cell
            "UHFk": _with_mode(base, "UHFk", Ncond=16, CellShape=[4, 4, 1]),
        }
        for mode, inp in cases.items():
            # keep only the keys this solver reads (RPA is not iterative;
            # UHFk has no Matsubara axis and takes Ncond)
            drop = {"RPA": ("IterationMax", "Mix", "EPS"),
                    "UHFk": ("Nmat", "filling")}[mode]
            for key in drop:
                inp["mode"]["param"].pop(key)
            with self.subTest(mode=mode, clean=True):
                self.assertEqual(self._reports(inp), [])
            with self.subTest(mode=mode, clean=False):
                bad = copy.deepcopy(inp)
                bad["mode"]["param"]["CellShap"] = [2, 2, 1]
                hits = self._reports(bad)
                self.assertEqual(len(hits), 1, hits)
                self.assertIn("CellShap", hits[0])
                self.assertIn("this {} run".format(mode), hits[0])

    def test_uhfr_reports_only_the_misspelled_key(self):
        import os
        import tomli
        here = os.path.dirname(os.path.abspath(__file__))
        work = os.path.join(self.tmp, "uhfr")
        shutil.copytree(os.path.join(here, "uhfr", "coulombintra"), work)
        with open(os.path.join(work, "input.toml"), "rb") as f:
            inp = tomli.load(f)
        inp["log"]["print_level"] = 0
        inp["mode"]["param"]["Mixx"] = 0.5
        cwd = os.getcwd()
        os.chdir(work)
        try:
            hits = self._reports(inp)
        finally:
            os.chdir(cwd)
        self.assertEqual(len(hits), 1, hits)
        self.assertIn("Mixx", hits[0])
        # RndSeed is not read when the initial Green function is loaded
        self.assertIn("RndSeed", hits[0])

    def test_misplaced_calc_scheme_stops_the_run(self):
        inp = _tiny_flex_input(self.tmp, extra_param={"calc_scheme": "general"})
        with self.assertLogs("qlms", level=logging.ERROR) as logs, \
                self.assertRaises(SystemExit):
            self._run(inp)
        self.assertTrue(any("move it to [mode]" in m for m in logs.output), logs.output)

    def test_the_caller_input_is_not_modified(self):
        inp = _tiny_flex_input(self.tmp)
        import hwave.qlms as qlms
        qlms.run(input_dict=inp)
        self.assertIs(type(inp["mode"]["param"]), dict)


if __name__ == "__main__":
    unittest.main()
