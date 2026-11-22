"""Tests of ``hwave.results.run``: opening a run, output names, run records."""
import json
import logging
import os
import shutil
import tempfile
import types
import unittest
from pathlib import Path

import numpy as np

import hwave
import hwave.results as results
from hwave.results import (
    MissingConfigError, MissingOutputError, ResultsError, Run, RunRecordError,
    UnsupportedLayoutError, as_run, open_run, read_run_record, run_record_path,
    write_run_record)
from tests.hf_square_oracle import neel_gap_equation
from tests.results_fixtures import (
    _freeze, assert_deep_frozen, copy_fixture, fixture)


def _write_toml(directory, text, name="input.toml"):
    path = Path(directory) / name
    path.write_text(text)
    return path


UHFK_HEAD = '[mode]\nmode = "UHFk"\n[mode.param]\nNcond = 2\n'
UHFR_HEAD = '[mode]\nmode = "UHFr"\n[mode.param]\nNsite = 2\nNcond = 2\n'


class _TmpDir(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)


class TestExports(unittest.TestCase):
    def test_all_names(self):
        for name in ("READER_VERSION", "Run", "open_run", "as_run",
                     "read_run_record", "write_run_record", "run_record_path",
                     "ResultsError", "UnsupportedLayoutError",
                     "MissingOutputError", "MissingConfigError",
                     "AmbiguousRunError", "RunRecordError"):
            self.assertIn(name, results.__all__)
            self.assertTrue(hasattr(results, name), name)
        self.assertEqual(results.READER_VERSION, "1")

    def test_error_classes(self):
        for cls in (results.UnsupportedLayoutError, results.MissingOutputError,
                    results.MissingConfigError, results.AmbiguousRunError,
                    results.RunRecordError):
            self.assertTrue(issubclass(cls, ResultsError))
            err = cls("message", path="/x/y", hint="do this")
            self.assertEqual(err.path, Path("/x/y"))
            self.assertEqual(err.hint, "do this")
            self.assertIn("message", str(err))
            self.assertIn("do this", str(err))
        err = ResultsError("plain")
        self.assertIsNone(err.path)
        self.assertIsNone(err.hint)
        self.assertEqual(str(err), "plain")


class TestSolverRuns(unittest.TestCase):
    """Output names of real solver runs (the two tutorial models)."""

    @classmethod
    def setUpClass(cls):
        cls.uhfk = fixture("uhfk_cdw")
        cls.uhfr = fixture("uhfr_tutorial")

    def test_uhfk_outputs(self):
        run = open_run(self.uhfk)
        self.assertEqual(run.mode, "UHFk")
        self.assertFalse(run.spin_orbital)
        self.assertEqual(run.path("eigen").name, "eigen.npz")
        self.assertEqual(run.path("green").name, "green.npz")
        self.assertEqual(run.path("scf_status").name, "scf_status.dat")
        self.assertEqual(run.path("energy").name, "energy.dat")
        self.assertEqual(run.output_dir, self.uhfk.parent.resolve() / "output")
        for key in ("eigen", "green", "scf_status", "energy"):
            self.assertEqual(run.path(key).parent, run.output_dir)
        self.assertEqual(run.input_path, self.uhfk.resolve())
        self.assertEqual(run.base_dir_source, "toml_directory")
        self.assertIsNone(run.run_record)

    def test_uhfr_eigen_blocks(self):
        run = open_run(self.uhfr)
        self.assertEqual(run.mode, "UHFr")
        path = run.path("eigen:spin-up")
        self.assertEqual(path.name, "spin-up_eigen.npz")
        self.assertTrue(path.is_file())
        self.assertTrue(run.path("eigen:spin-down").is_file())
        with self.assertRaises(MissingOutputError) as cm:
            run.path("eigen:sz-free")
        self.assertIn("eigen:sz-free", str(cm.exception))
        self.assertIn("not configured", str(cm.exception))
        self.assertIn("eigen:spin-up, eigen:spin-down", cm.exception.hint)
        self.assertNotIn("eigen", run.outputs)
        self.assertEqual(run.path("green").name, "green.dat")
        self.assertEqual(run.path("scf_status").name, "scf_status.dat")

    def test_root_for_each_source(self):
        toml = copy_fixture("uhfk_cdw")
        tree = toml.parent
        other = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, other, ignore_errors=True)

        # 1. base_dir= (a relative one is resolved)
        cwd = os.getcwd()
        try:
            os.chdir(other)
            run = open_run(toml, base_dir=".")
        finally:
            os.chdir(cwd)
        self.assertEqual(run.root, other.resolve())
        self.assertEqual(run.base_dir_source, "argument")
        self.assertEqual(run.output_dir, other.resolve() / "output")

        # 2. a run record
        record_cwd = other / "recorded"
        write_run_record(toml, cwd=record_cwd, argv=["hwave", "input.toml"])
        run = open_run(toml)
        self.assertEqual(run.root, Path(read_run_record(toml)["cwd"]))
        self.assertEqual(run.root, record_cwd.resolve())
        self.assertEqual(run.base_dir_source, "run_record")
        self.assertEqual(run.run_record["argv"], ["hwave", "input.toml"])

        # 3. the TOML directory
        run_record_path(toml).unlink()
        run = open_run(toml)
        self.assertEqual(run.root, tree.resolve())
        self.assertEqual(run.base_dir_source, "toml_directory")
        self.assertIsNone(run.run_record)

    def test_root_choice_logged_at_debug(self):
        with self.assertLogs("hwave.results", level="DEBUG") as cm:
            open_run(self.uhfk)
        text = "\n".join(cm.output)
        self.assertIn("toml_directory", text)
        self.assertIn(str(self.uhfk.parent.resolve()), text)

    def test_missing_output_message(self):
        toml = copy_fixture("uhfk_cdw")
        green = toml.parent / "output" / "green.npz"
        green.unlink()
        run = open_run(toml)
        with self.assertRaises(MissingOutputError) as cm:
            run.path("green")
        err = cm.exception
        self.assertEqual(err.path, green.resolve())
        msg = str(err)
        self.assertIn(str(green.resolve()), msg)
        self.assertIn("toml_directory", msg)
        self.assertIn("base_dir=", msg)
        # a key that is not configured at all
        with self.assertRaises(MissingOutputError) as cm:
            run.path("onebodyg")
        msg = str(cm.exception)
        self.assertIn("onebodyg", msg)
        self.assertIn("not configured", msg)
        self.assertIn("energy", msg)          # the configured keys are listed

    def test_run_record_cwd(self):
        # the TOML sits in a sub-directory and the run was started from the
        # parent: only the run record knows where the outputs are
        toml = copy_fixture("uhfk_cdw")
        work = toml.parent
        conf = work / "conf"
        conf.mkdir()
        moved = conf / toml.name
        shutil.move(str(toml), str(moved))
        with self.assertRaises(MissingOutputError):
            open_run(moved).path("energy")
        record_file = write_run_record(moved, cwd=work, argv=["hwave", "conf/input.toml"])
        self.assertEqual(record_file, run_record_path(moved))
        run = open_run(moved)
        self.assertEqual(run.base_dir_source, "run_record")
        self.assertEqual(run.root, work.resolve())
        self.assertTrue(run.path("energy").is_file())
        self.assertEqual(run.meta["run_record_path"], str(record_file))


class TestOutputNames(_TmpDir):
    def test_npz_suffix(self):
        toml = _write_toml(self.tmp, UHFK_HEAD + (
            '[file.output]\neigen = "eigen.npz"\ngreen = "green.dat"\n'
            'rpa = "trans"\nonebodyg = "greenone.dat"\nenergy = "energy.dat"\n'))
        run = open_run(toml)
        out = self.tmp.resolve() / "output"
        self.assertEqual(run.outputs["eigen"], out / "eigen.npz")
        self.assertEqual(run.outputs["green"], out / "green.dat.npz")
        self.assertEqual(run.outputs["rpa"], out / "trans.npz")
        self.assertEqual(run.outputs["onebodyg"], out / "greenone.dat")   # text
        self.assertEqual(run.outputs["energy"], out / "energy.dat")       # text
        toml = _write_toml(self.tmp, UHFR_HEAD + (
            '[file.output]\neigen = "eig.npz"\ngreen = "green.dat"\n'))
        run = open_run(toml)
        self.assertEqual(run.outputs["eigen:sz-free"], out / "sz-free_eig.npz")
        self.assertEqual(run.outputs["green"], out / "green.dat")         # text

    def test_uhfk_defaults(self):
        toml = _write_toml(self.tmp, UHFK_HEAD)
        run = open_run(toml)
        out = self.tmp.resolve() / "output"
        self.assertEqual(run.outputs, {"green": out / "green.npz",
                                       "scf_status": out / "scf_status.dat"})
        toml = _write_toml(self.tmp, UHFK_HEAD + (
            '[file.output]\npath_to_output = "res"\ngreen = ""\nscf_status = ""\n'))
        self.assertEqual(open_run(toml).outputs, {})
        self.assertEqual(open_run(toml).output_dir, self.tmp.resolve() / "res")

    def test_uhfr_blocks_and_fij(self):
        cases = [("", ["sz-free"], "sz-free_f.dat"),
                 ("2Sz = \"free\"\n", ["sz-free"], "sz-free_f.dat"),
                 ("2Sz = 0\n", ["spin-up", "spin-down"], "Sz0_f.dat"),
                 ("2Sz = 2\n", ["spin-up", "spin-down"], None)]
        for extra, blocks, fij in cases:
            with self.subTest(extra=extra):
                toml = _write_toml(self.tmp, UHFR_HEAD + extra + (
                    '[file.output]\neigen = "eigen"\nfij = "f.dat"\n'
                    'initial = "init.dat"\n'))
                run = open_run(toml)
                out = self.tmp.resolve() / "output"
                eig = sorted(k for k in run.outputs if k.startswith("eigen:"))
                self.assertEqual(eig, sorted("eigen:" + b for b in blocks))
                for b in blocks:
                    self.assertEqual(run.outputs["eigen:" + b],
                                     out / (b + "_eigen.npz"))
                if fij is None:
                    self.assertNotIn("fij", run.outputs)
                else:
                    self.assertEqual(run.outputs["fij"], out / fij)
                self.assertEqual(run.outputs["initial"], out / "init.dat")
                self.assertNotIn("green", run.outputs)     # no UHFr default

    def test_print_check(self):
        for head in (UHFK_HEAD, UHFR_HEAD):
            toml = _write_toml(self.tmp, '[log]\nprint_check = "check.dat"\n' + head +
                               '[file.output]\npath_to_output = "o"\n')
            run = open_run(toml)
            self.assertEqual(run.outputs["print_check"], self.tmp.resolve() / "o" / "check.dat")
            # the [log] key is exact-case, as the solvers read it
            toml = _write_toml(self.tmp, '[log]\nPrint_Check = "check.dat"\n' + head)
            run = open_run(toml)
            self.assertNotIn("print_check", run.outputs)
            self.assertIn("[log] Print_Check", run.meta["ignored_output_keys"])

    def test_output_keys_exact_case(self):
        # the solvers read [file.output] keys exactly (uhfk.py:1796, uhfr.py:1271):
        # a key in other case is not written, so it is not configured here;
        # it is reported instead of being dropped without a word
        toml = _write_toml(self.tmp, UHFK_HEAD + (
            '[file.output]\nEnergy = "energy.dat"\nSCF_Status = "s.dat"\n'))
        run = open_run(toml)
        self.assertNotIn("energy", run.outputs)
        self.assertEqual(run.outputs["scf_status"].name, "scf_status.dat")
        self.assertEqual(sorted(run.meta["ignored_output_keys"]),
                         ["[file.output] Energy", "[file.output] SCF_Status"])

    def test_invalid_output_values(self):
        # a non-string or empty name is not configured: the reader does not
        # resolve it, although the solver may still write a file for an
        # empty name (UHFk eigen = "" gives output/.npz)
        toml = _write_toml(self.tmp, '[log]\nprint_check = 3\n' + UHFK_HEAD + (
            '[file.output]\nenergy = ""\neigen = 3\nscf_status = 1\n'))
        run = open_run(toml)
        self.assertNotIn("energy", run.outputs)
        self.assertNotIn("eigen", run.outputs)
        self.assertNotIn("scf_status", run.outputs)
        self.assertNotIn("print_check", run.outputs)
        self.assertEqual(sorted(run.meta["invalid_output_keys"]),
                         ["[file.output] eigen", "[file.output] scf_status",
                          "[log] print_check"])

    def test_names_match_solver(self):
        """Run.outputs and the solvers' own list of written files agree."""
        from hwave.solver.uhfk import UHFk
        from hwave.solver.uhfr import UHFr
        info = {"path_to_output": "out", "energy": "e.dat", "eigen": "eig",
                "green": "g.dat", "initial": "i.dat", "fij": "f.dat",
                "onebodyg": "o.dat", "rpa": "r.npz"}
        output = '[file.output]\n' + "".join(
            '{} = "{}"\n'.format(k, v) for k, v in info.items())
        for two_sz, blocks in ((None, ["sz-free"]), (0, ["spin-up", "spin-down"]),
                               (2, ["spin-up", "spin-down"])):
            fake = types.SimpleNamespace(
                green_list={b: {} for b in blocks}, param_mod={"2Sz": two_sz})
            expected = sorted(str(self.tmp.resolve() / p)
                              for _, p in UHFr._output_artifacts(fake, info))
            extra = "" if two_sz is None else "2Sz = {}\n".format(two_sz)
            run = open_run(_write_toml(self.tmp, UHFR_HEAD + extra + output))
            got = sorted(str(p) for k, p in run.outputs.items() if k != "scf_status")
            self.assertEqual(got, expected)
        fake = types.SimpleNamespace(param_ham={})
        expected = sorted(str(self.tmp.resolve() / p)
                          for _, p in UHFk._output_artifacts(fake, info))
        run = open_run(_write_toml(self.tmp, UHFK_HEAD + output))
        got = sorted(str(p) for k, p in run.outputs.items() if k != "scf_status")
        self.assertEqual(got, expected)

    def test_unsupported_mode_opens_empty(self):
        toml = _write_toml(self.tmp, '[mode]\nmode = "RPA"\n[file.output]\n'
                                     'chi0q = "chi0q"\n')
        run = open_run(toml)
        self.assertEqual(run.mode, "RPA")
        self.assertEqual(run.outputs, {})
        self.assertEqual(run.meta["unsupported_mode"], "RPA")
        with self.assertRaises(MissingOutputError) as cm:
            run.path("chi0q")
        self.assertIn("RPA", str(cm.exception))


class TestOpenRunRefusals(_TmpDir):
    def test_directory_refused(self):
        with self.assertRaises(MissingConfigError) as cm:
            open_run(self.tmp)
        err = cm.exception
        self.assertEqual(err.hint, "pass the input TOML file; directory discovery "
                                   "is not part of this API")
        self.assertIn("directory discovery", str(err))
        self.assertEqual(err.path, self.tmp)

    def test_missing_file(self):
        with self.assertRaises(MissingConfigError) as cm:
            open_run(self.tmp / "nope.toml")
        self.assertEqual(cm.exception.path, self.tmp / "nope.toml")

    def test_no_mode_table(self):
        toml = _write_toml(self.tmp, '[log]\nprint_level = 1\n')
        with self.assertRaises(UnsupportedLayoutError) as cm:
            open_run(toml)
        self.assertIn("[mode]", str(cm.exception))
        for text in ('mode = "UHFk"\n', '[mode]\nparam = {}\n',
                     '[mode]\nmode = 3\n', '[mode\n'):
            with self.subTest(text=text):
                with self.assertRaises(UnsupportedLayoutError):
                    open_run(_write_toml(self.tmp, text))

    def test_bad_tables(self):
        for text in ('file = 3\n' + UHFK_HEAD,
                     UHFK_HEAD + '[file]\noutput = 3\n',
                     UHFK_HEAD + '[file.output]\npath_to_output = 3\n'):
            with self.subTest(text=text):
                with self.assertRaises(UnsupportedLayoutError):
                    open_run(_write_toml(self.tmp, text))

    def test_as_run_rejects_base_dir_for_run(self):
        toml = _write_toml(self.tmp, UHFK_HEAD)
        run = open_run(toml)
        self.assertIs(as_run(run), run)
        with self.assertRaises(MissingConfigError) as cm:
            as_run(run, base_dir=self.tmp)
        self.assertIn("base_dir", str(cm.exception))
        self.assertIn("open_run(toml, base_dir=", cm.exception.hint)
        other = as_run(toml, base_dir=self.tmp)
        self.assertIsInstance(other, Run)
        self.assertEqual(other.base_dir_source, "argument")
        self.assertEqual(as_run(str(toml)).input_path, toml.resolve())
        with self.assertRaises(MissingConfigError):
            as_run(self.tmp)

    def test_run_is_frozen(self):
        run = open_run(_write_toml(self.tmp, UHFK_HEAD))
        with self.assertRaises(Exception):
            run.root = Path("/")
        # outputs and meta are read-only
        with self.assertRaises(TypeError):
            run.outputs["green"] = Path("/elsewhere")
        with self.assertRaises(TypeError):
            del run.outputs["green"]
        with self.assertRaises(TypeError):
            run.meta["note"] = "x"
        self.assertEqual(run.outputs["green"], self.tmp.resolve() / "output" / "green.npz")
        # every open builds fresh mappings
        again = open_run(run.input_path)
        self.assertIsNot(run.outputs, again.outputs)
        self.assertIsNot(run.meta, again.meta)
        self.assertIsNot(run.config, again.config)


class TestConfigReadOnly(unittest.TestCase):
    """Run.config is read-only at every depth: the readers consult it after
    open_run, so a caller cannot change what they see."""

    def test_config_cannot_be_changed(self):
        toml = copy_fixture("uhfk_u8")
        out = toml.parent / "output"
        text = toml.read_text()
        self.assertIn("[file.output]\n", text)
        toml.write_text(text.replace('[file.output]\n',
                                     '[file.output]\nscf_status = ""\n'))
        (out / "scf_status.dat").unlink()
        run = open_run(toml)
        expected = results.read_energy(run).values
        cfg = run.config
        for mutate in (lambda: cfg.__setitem__("mode", {}),
                       lambda: cfg["mode"]["param"].__setitem__("IterationMax", 0),
                       lambda: cfg["mode"]["param"].__setitem__("iterationmax", 0),
                       lambda: cfg["mode"]["param"].pop("IterationMax"),
                       lambda: cfg["file"]["output"].__setitem__("energy", "x.dat"),
                       lambda: cfg["file"]["input"]["interaction"].clear()):
            with self.assertRaises((TypeError, AttributeError)):
                mutate()
        self.assertIsInstance(cfg["mode"]["param"]["CellShape"], tuple)
        # the case rules survive behind the read-only view
        self.assertEqual(cfg["mode"]["param"]["iterationmax"], 1000)
        self.assertNotIn("Mode", cfg)
        self.assertEqual(run.param("IterationMax"), 1000)
        self.assertEqual(dict(results.read_energy(run).values), dict(expected))
        # a mutable copy is the caller's own
        copy = dict(cfg)
        copy["mode"] = {"param": {"IterationMax": 0}}
        self.assertEqual(run.param("IterationMax"), 1000)
        results.read_energy(run)


class TestSpinOrbitalFlag(_TmpDir):
    def test_flag(self):
        for value, expected in (("true", True), ("false", False), ("1.0", True)):
            toml = _write_toml(self.tmp, '[mode]\nmode = "UHFk"\nenable_spin_orbital = {}\n'
                               .format(value))
            self.assertIs(open_run(toml).spin_orbital, expected)
        self.assertIs(open_run(_write_toml(self.tmp, UHFK_HEAD)).spin_orbital, False)

    def test_bad_flag(self):
        toml = _write_toml(self.tmp, '[mode]\nmode = "UHFk"\nenable_spin_orbital = "maybe"\n')
        with self.assertRaises(UnsupportedLayoutError) as cm:
            open_run(toml)
        self.assertIn("enable_spin_orbital", str(cm.exception))
        # a float is accepted for UHFk only, as the solvers do
        toml = _write_toml(self.tmp, '[mode]\nmode = "UHFr"\nenable_spin_orbital = 1.0\n')
        with self.assertRaises(UnsupportedLayoutError):
            open_run(toml)


class TestCaseInsensitiveConfig(_TmpDir):
    """Run.config reads each layer as the solver does: the top-level tables,
    [mode] and [file.*] are exact-case (qlms.py:98-113, 26-53; uhfk.py and
    uhfr.py save_results; read_input_k.py:152-190), [mode.param] and
    [file.input.interaction] are case-insensitive (base.py ParamDict,
    qlms.py:137; qlms.py:146-149, read_input_k.py:93)."""
    MIXED = (
        '[mode]\nmode = "UHFr"\n'
        '[mode.param]\nNSITE = 2\nncond = 2\n2sz = 0\neps = 9\nMIX = 0.25\n'
        '[file]\n[file.input]\nPath_To_Input = "in"\nonebodyg = "g.def"\n'
        '[file.input.interaction]\nTRANS = "trans.def"\n'
        '[file.output]\nenergy = "energy.dat"\neigen = "eigen"\n')

    def test_mixed_case_toml(self):
        with self.assertLogs("hwave.results.run", level="WARNING") as logs:
            run = open_run(_write_toml(self.tmp, self.MIXED))
        cfg = run.config
        # case-insensitive layers: [mode.param] and the interaction table
        self.assertEqual(cfg["mode"]["param"]["EPS"], 9)
        self.assertEqual(cfg["mode"]["param"]["Nsite"], 2)
        self.assertEqual(cfg["file"]["input"]["interaction"]["Trans"], "trans.def")
        self.assertEqual(run.param("EPS"), 9)
        self.assertEqual(run.param("mix"), 0.25)
        self.assertEqual(run.param("MIX", default=1.0), 0.25)
        self.assertEqual(run.param("2Sz"), 0)
        self.assertEqual(run.param("NSite"), 2)
        # exact-case layers: the tables, [file.input] and [file.output]
        for lookup in (lambda: cfg["MODE"], lambda: cfg["File"],
                       lambda: cfg["mode"]["Param"],
                       lambda: cfg["file"]["INPUT"],
                       lambda: cfg["file"]["input"]["path_to_input"],
                       lambda: cfg["file"]["output"]["ENERGY"]):
            with self.assertRaises(KeyError):
                lookup()
        self.assertEqual(cfg["file"]["input"]["Path_To_Input"], "in")
        self.assertEqual(cfg["file"]["output"]["energy"], "energy.dat")
        # UHFr reads [file.input] path_to_input with its exact case
        # (qlms.py:108), so the variant is reported and logged
        self.assertEqual(run.meta["ignored_keys"], ("[file.input] Path_To_Input",))
        self.assertTrue(any("Path_To_Input" in line for line in logs.output))
        # the lower-case 2sz selects the two spin blocks, as in the solver
        self.assertIn("eigen:spin-up", run.outputs)
        self.assertIn("eigen:spin-down", run.outputs)

    def test_uhfk_input_key_in_other_case_not_honoured(self):
        """UHFk reads [file.input] keys exactly (read_input_k.py:155):
        ``Initial`` is not a start, so it is reported and logged."""
        toml = _write_toml(self.tmp, UHFK_HEAD + (
            '[file.input]\nInitial = "start.npz"\ninitial_mode = "random"\n'))
        with self.assertLogs("hwave.results.run", level="WARNING") as logs:
            run = open_run(toml)
        self.assertNotIn("initial", run.config["file"]["input"])
        self.assertEqual(run.config["file"]["input"]["initial_mode"], "random")
        self.assertEqual(run.meta["ignored_keys"], ("[file.input] Initial",))
        self.assertEqual(len(logs.output), 1)
        self.assertIn("[file.input] Initial", logs.output[0])
        self.assertIn("exact case", logs.output[0])

    def test_uhfr_initial_and_onebodyg_any_case(self):
        """UHFr matches [file.input] initial and onebodyg in any case
        (qlms.py:155-159): no report for those."""
        toml = _write_toml(self.tmp, UHFR_HEAD + (
            '[file.input]\nInitial = "init.dat"\nOneBodyG = "g.def"\n'))
        logger = logging.getLogger("hwave.results.run")
        with self.assertRaises(AssertionError):
            with self.assertLogs(logger, level="WARNING"):
                run = open_run(toml)
        run = open_run(toml)
        self.assertNotIn("ignored_keys", run.meta)
        self.assertEqual(run.config["file"]["input"]["OneBodyG"], "g.def")

    def test_tables_and_mode_keys_in_other_case(self):
        # a [Mode] table is not read by the solver (qlms.py:103, 122-124)
        toml = _write_toml(self.tmp, '[Mode]\nmode = "UHFk"\n')
        with self.assertRaisesRegex(UnsupportedLayoutError, r"\[Mode\].*exact case"):
            open_run(toml)
        # [File] and [Log] are ignored by the solver, a [mode] key in other
        # case is refused by it (qlms.py:26-53); the reader reports them
        toml = _write_toml(self.tmp, (
            '[mode]\nmode = "UHFk"\nEnable_Spin_Orbital = true\n'
            '[mode.param]\nNcond = 2\n'
            '[File.output]\nenergy = "e.dat"\n[Log]\nprint_check = "c.dat"\n'
            '[file.Input]\ninitial = "s.npz"\n'
            '[file.input]\nInteraction = {Geometry = "geom.dat"}\n'))
        with self.assertLogs("hwave.results.run", level="WARNING") as logs:
            run = open_run(toml)
        self.assertEqual(sorted(run.meta["ignored_keys"]), [
            "[File]", "[Log]", "[file.input] Interaction", "[file] Input",
            "[mode] Enable_Spin_Orbital"])
        self.assertEqual(len(logs.output), 5)
        self.assertFalse(run.spin_orbital)
        self.assertNotIn("energy", run.outputs)
        self.assertNotIn("print_check", run.outputs)

    def test_output_keys_also_in_ignored_keys(self):
        toml = _write_toml(self.tmp, UHFK_HEAD + '[file.output]\nEnergy = "e.dat"\n')
        run = open_run(toml)
        self.assertEqual(run.meta["ignored_output_keys"], ("[file.output] Energy",))
        self.assertEqual(run.meta["ignored_keys"], ("[file.output] Energy",))

    def test_defaults(self):
        run = open_run(_write_toml(self.tmp, UHFR_HEAD))
        self.assertEqual(run.param("Mix"), 0.5)
        self.assertEqual(run.param("EPS"), 6)
        self.assertEqual(run.param("IterationMax"), 20000)
        self.assertIsNone(run.param("2Sz"))
        self.assertIsNone(run.param("CellShape"))
        self.assertEqual(run.param("CellShape", default=[1, 1, 1]), [1, 1, 1])
        # a default given by the caller wins over the solver default ...
        self.assertEqual(run.param("Mix", default=0.3), 0.3)
        self.assertIsNone(run.param("EPS", default=None))
        # ... but not over a value in the file
        self.assertEqual(run.param("ncond", default=99), 2)
        run = open_run(_write_toml(self.tmp, UHFR_HEAD + '2Sz = "free"\n'))
        self.assertIsNone(run.param("2sz"))

    def test_defaults_match_solver(self):
        """The defaults Run.param uses are those of solver_base (base.py)."""
        from hwave.results.run import SOLVER_PARAM_DEFAULTS
        from hwave.solver.base import solver_base

        class _Bare(solver_base):
            name = "results_test"

        bare = _Bare({}, {}, {"mode": "UHFr", "param": {"Ncond": 1}})
        solver = {k: v for k, v in bare.param_mod.items()}
        self.assertEqual(sorted(k.lower() for k in solver),
                         sorted(k.lower() for k in SOLVER_PARAM_DEFAULTS))
        for key, value in SOLVER_PARAM_DEFAULTS.items():
            if key == "EPS":            # the solver stores 10**-EPS
                self.assertEqual(solver[key], 10 ** -value)
            elif key == "Ncond":
                continue                # given above (the solver requires it)
            else:
                self.assertEqual(solver[key], value, key)


class TestRunRecords(_TmpDir):
    def test_path_and_schema(self):
        toml = _write_toml(self.tmp, UHFK_HEAD, name="model.a.toml")
        self.assertEqual(run_record_path(toml), self.tmp / "hwave-run.model.a.json")
        work = self.tmp / "work"
        (self.tmp / "sub").mkdir()
        link_dir = self.tmp / "sub" / ".."
        path = write_run_record(link_dir / "model.a.toml", cwd=work,
                                argv=["hwave", Path("model.a.toml")])
        self.assertEqual(path, run_record_path(link_dir / "model.a.toml"))
        with open(path) as f:
            record = json.load(f)
        self.assertEqual(record, {
            "format": "hwave-run-record", "format_version": 1,
            "cwd": str(work.resolve()), "input": str(toml.resolve()),
            "argv": ["hwave", "model.a.toml"],
            "hwave_version": hwave.__version__})
        self.assertEqual(read_run_record(toml), record)
        self.assertIsNone(read_run_record(self.tmp / "other.toml"))
        with self.assertRaises(MissingConfigError):
            write_run_record(self.tmp / "other.toml", cwd=work, argv=[])

    def test_run_record_bad_version(self):
        toml = _write_toml(self.tmp, UHFK_HEAD)
        good = {"format": "hwave-run-record", "format_version": 1,
                "cwd": str(self.tmp.resolve()), "input": str(toml.resolve()),
                "argv": [], "hwave_version": "2.0.0"}
        bad = [dict(good, format_version=2), dict(good, format_version="1"),
               dict(good, format="other"), {"cwd": str(self.tmp)},
               dict(good, cwd="relative/dir"), dict(good, cwd=3),
               dict(good, input=None), [1, 2]]
        for record in bad:
            with self.subTest(record=record):
                run_record_path(toml).write_text(json.dumps(record))
                with self.assertRaises(RunRecordError) as cm:
                    read_run_record(toml)
                self.assertEqual(cm.exception.path, run_record_path(toml))
                with self.assertRaises(RunRecordError):
                    open_run(toml)
        run_record_path(toml).write_text("{not json")
        with self.assertRaises(RunRecordError):
            open_run(toml)
        # with base_dir= the record is not read at all
        run = open_run(toml, base_dir=self.tmp)
        self.assertIsNone(run.run_record)
        self.assertIn("not used", run.meta["run_record_note"])

    def test_run_record_mismatch_and_bypass(self):
        first = self.tmp / "first"
        first.mkdir()
        toml = _write_toml(first, UHFK_HEAD)
        write_run_record(toml, cwd=first, argv=["hwave", "input.toml"])
        self.assertEqual(open_run(toml).base_dir_source, "run_record")
        moved = self.tmp / "moved"
        shutil.move(str(first), str(moved))
        toml = moved / "input.toml"
        with self.assertRaises(RunRecordError) as cm:
            open_run(toml)
        msg = str(cm.exception)
        self.assertIn(str((first / "input.toml").resolve()), msg)
        self.assertIn(str(toml.resolve()), msg)
        self.assertIn("base_dir=", msg)
        run = open_run(toml, base_dir=moved)
        self.assertEqual(run.base_dir_source, "argument")
        self.assertEqual(run.root, moved.resolve())
        self.assertIsNone(run.run_record)
        note = run.meta["run_record_note"]
        self.assertIn(str(run_record_path(toml)), note)
        self.assertIn("base_dir=", note)
        # no note when there is no record
        run_record_path(toml).unlink()
        self.assertNotIn("run_record_note", open_run(toml, base_dir=moved).meta)

    def test_record_input_differs_only_in_case(self):
        toml = _write_toml(self.tmp, UHFK_HEAD)
        write_run_record(toml, cwd=self.tmp, argv=["hwave", "input.toml"])
        other_case = self.tmp / "INPUT.toml"
        if not other_case.exists():
            self.skipTest("the file system is case-sensitive")
        # the same file: not reported as moved
        self.assertEqual(read_run_record(other_case)["input"], str(toml.resolve()))
        self.assertEqual(open_run(other_case).base_dir_source, "run_record")


class TestFreezeHelpers(unittest.TestCase):
    def test_freeze(self):
        src = np.arange(4.0)
        nested = {"a": src, "b": [np.ones(2), (np.zeros(1), 3)], "c": "x"}
        frozen = _freeze(nested)
        assert_deep_frozen(frozen)
        self.assertIsNot(frozen, nested)
        self.assertIsNot(frozen["a"], src)
        self.assertTrue(src.flags.writeable)      # the input is not touched
        src[0] = 9.0
        self.assertEqual(frozen["a"][0], 0.0)     # a private copy
        self.assertEqual(frozen["c"], "x")
        with self.assertRaises(ValueError):
            frozen["b"][1][0][0] = 1.0

    def test_freeze_case_insensitive_dict(self):
        from requests.structures import CaseInsensitiveDict
        src = np.arange(3.0)
        mapping = CaseInsensitiveDict({"Green": src, "inner": {"a": [src]}})
        frozen = _freeze(mapping)
        self.assertIsInstance(frozen, CaseInsensitiveDict)
        self.assertIsNot(frozen["green"], src)
        assert_deep_frozen(frozen)
        with self.assertRaises(ValueError):
            frozen["GREEN"][0] = 1.0
        self.assertTrue(src.flags.writeable)
        with self.assertRaises(AssertionError):
            assert_deep_frozen(mapping)
        from types import MappingProxyType
        with self.assertRaises(AssertionError):
            assert_deep_frozen(MappingProxyType({"x": src}))
        self.assertEqual(type(_freeze(MappingProxyType({"x": src}))), dict)

    def test_assert_deep_frozen_detects(self):
        with self.assertRaises(AssertionError):
            assert_deep_frozen({"x": [np.zeros(2)]})
        base = np.zeros(4)
        view = base[:2]
        view.flags.writeable = False
        with self.assertRaises(AssertionError):
            assert_deep_frozen(view)               # the base is still writeable
        import dataclasses

        @dataclasses.dataclass
        class Loose:
            a: int
        with self.assertRaises(AssertionError):
            assert_deep_frozen(Loose(1))

        @dataclasses.dataclass(frozen=True)
        class Tight:
            a: object
        assert_deep_frozen(Tight(_freeze(np.zeros(2))))
        with self.assertRaises(AssertionError):
            assert_deep_frozen(Tight(np.zeros(2)))


class TestSquareOracle(unittest.TestCase):
    def test_golden_values(self):
        m, e = neel_gap_equation(4, 8.0)
        self.assertAlmostEqual(m, 0.4474190440247483, delta=1e-12)
        self.assertAlmostEqual(e, -7.38962582304832, delta=1e-10)
        m, e = neel_gap_equation(16, 2.0)
        self.assertAlmostEqual(m, 0.1917429137967485, delta=1e-12)
        self.assertAlmostEqual(e, -291.4365477135228, delta=1e-9)


if __name__ == "__main__":
    unittest.main()
