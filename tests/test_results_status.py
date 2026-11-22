"""Tests of ``hwave.results.status`` (the convergence status) and of the
``energy.dat`` reader of ``hwave.results.uhf``."""
import math
import re
import shutil
import unittest

import hwave.results as results
from hwave.results import (
    Energy, MissingOutputError, ScfStatus, UnsupportedLayoutError, open_run,
    read_energy, read_scf_status)
from hwave.results.status import NOT_RUN_NOTE, PREDATES_NOTE
from tests.results_fixtures import assert_deep_frozen, copy_fixture, fixture

#: a valid converged status file (format version 1), as UHFk writes it
GOOD = (
    ("format_version", "1"),
    ("solver", "UHFk"),
    ("converged", "true"),
    ("residual", "8.537948717029953e-16"),
    ("eps", "1e-12"),
    ("iterations", "1"),
    ("iteration_max", "1000"),
    ("residual_definition", "two_frobenius_over_size:uhfk_all_r"),
    ("note", ""),
)


def status_text(rows=GOOD, drop=(), extra=(), **changes):
    """The text of a status file: ``rows`` with ``changes`` applied, the keys
    in ``drop`` left out and the ``(key, value)`` pairs of ``extra``
    appended."""
    lines = []
    for key, value in rows:
        if key in drop:
            continue
        lines.append("{} = {}".format(key, changes.get(key, value)))
    lines.extend("{} = {}".format(k, v) for k, v in extra)
    return "\n".join(lines) + "\n"


def file_numbers(path):
    """``{key: float}`` of a ``key = value`` file, parsed independently of
    the reader (the solver's values differ in the last digits between
    platforms, so the tests compare the reader with the file and the file
    with the physics, to a tolerance)."""
    table = {}
    for line in path.read_text().splitlines():
        key, _, value = line.partition("=")
        table[key.strip()] = float(value)
    return table


def names_key(key):
    """The pattern of an error message that names ``key``."""
    return r": {!r}: ".format(key)


class TestExports(unittest.TestCase):
    def test_names(self):
        for name in ("Energy", "read_energy", "ScfStatus", "read_scf_status"):
            self.assertIn(name, results.__all__)
            self.assertTrue(hasattr(results, name))


class TestStatusFromFixtures(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.u8 = fixture("uhfk_u8").resolve()
        cls.not_converged = fixture("uhfk_u8_not_converged").resolve()
        cls.zero = fixture("uhfk_u8_zero_steps").resolve()

    def test_converged(self):
        st = read_scf_status(self.u8)
        self.assertIsInstance(st, ScfStatus)
        self.assertEqual(st.status, "converged")
        self.assertIs(st.converged, True)
        line, = [l for l in st.source.read_text().splitlines()
                 if l.startswith("residual = ")]
        self.assertEqual(st.residual, float(line.split("=")[1]))
        self.assertLess(st.residual, 1e-14)
        self.assertEqual(st.eps, 1e-12)
        self.assertEqual(st.iterations, 1)
        self.assertEqual(st.iteration_max, 1000)
        self.assertEqual(st.residual_definition, "two_frobenius_over_size:uhfk_all_r")
        self.assertEqual(st.source_kind, "scf_status")
        self.assertEqual(st.note, "")
        self.assertEqual(st.meta["note_kind"], "")
        self.assertEqual(st.meta["format_version"], 1)
        self.assertEqual(st.meta["solver"], "UHFk")
        self.assertEqual(st.source, self.u8.parent / "output" / "scf_status.dat")
        assert_deep_frozen(st)

    def test_accepts_a_run(self):
        run = open_run(self.u8)
        self.assertEqual(read_scf_status(run).status, "converged")

    def test_fresh_meta(self):
        a, b = read_scf_status(self.u8), read_scf_status(self.u8)
        self.assertIsNot(a.meta, b.meta)

    def test_not_converged(self):
        st = read_scf_status(self.not_converged)
        self.assertEqual(st.status, "not_converged")
        self.assertIs(st.converged, False)
        self.assertEqual(st.iterations, 5)
        self.assertEqual(st.iteration_max, 5)
        self.assertGreater(st.residual, st.eps)

    def test_not_converged_started_from_random(self):
        # the first step differs from that of a zero start (built next to the
        # fixture by its builder, which also refuses the fixture otherwise)
        # and from the committed start, which converges in one step at the
        # Neel energy
        from tests.results_fixtures import _first_check_line
        run_dir = self.not_converged.parent
        first = _first_check_line(run_dir)
        zero = _first_check_line(run_dir.parent / (run_dir.name + "-zero-start"))
        self.assertNotEqual(first[1:], zero[1:])
        self.assertGreater(abs(float(first[2]) + 7.38962582), 1e-3)
        self.assertNotIn("initial", open_run(self.not_converged).config["file"]["input"])

    def test_zero_steps(self):
        st = read_scf_status(self.zero)
        self.assertEqual(st.status, "not_converged")
        self.assertIs(st.converged, False)
        self.assertEqual(st.iterations, 0)
        self.assertEqual(st.iteration_max, 0)
        self.assertTrue(math.isnan(st.residual))
        self.assertEqual(st.note, "no SCF step performed")
        self.assertEqual(st.meta["note_kind"], "no_step")

    def test_uhfr(self):
        st = read_scf_status(fixture("uhfr_tutorial"))
        self.assertEqual(st.status, "converged")
        self.assertEqual(st.iterations, 34)
        self.assertEqual(st.eps, 1e-8)
        self.assertEqual(st.residual_definition, "two_frobenius_over_size:uhfr_full")

    def test_legacy_is_unknown(self):
        # energy.dat, but no status file and no print_check
        st = read_scf_status(fixture("uhfk_cdw_legacy"))
        self.assertEqual(st.status, "unknown")
        self.assertIsNone(st.converged)
        self.assertEqual(st.source_kind, "none")
        self.assertIsNone(st.source)
        self.assertEqual(st.note, PREDATES_NOTE)
        self.assertIn("predates the convergence status file", st.note)


class TestWriterForms(unittest.TestCase):
    """The NaN-residual and negative-IterationMax files as the solver writes
    them."""

    def test_nan_residual(self):
        toml = fixture("uhfk_u8_nan_residual").resolve()
        text = (toml.parent / "output" / "scf_status.dat").read_text()
        self.assertIn("residual = nan\n", text)
        self.assertIn("note = residual not finite after step 5\n", text)
        st = read_scf_status(toml)
        self.assertEqual(st.status, "not_converged")
        self.assertIs(st.converged, False)
        self.assertIsNone(st.residual)
        self.assertEqual(st.iterations, 5)
        self.assertEqual(st.meta["note_kind"], "residual_not_finite")
        self.assertEqual(st.meta["nonfinite_step"], 5)
        self.assertIs(st.meta["residual_not_finite"], True)

    def test_nan_residual_print_check(self):
        toml = copy_fixture("uhfk_u8_nan_residual").resolve()
        self.addCleanup(shutil.rmtree, toml.parent, ignore_errors=True)
        _without(toml, "scf_status.dat")
        st = read_scf_status(toml)
        self.assertEqual(st.status, "not_converged")
        self.assertEqual(st.iterations, 5)
        self.assertIsNone(st.residual)
        self.assertIs(st.meta["residual_not_finite"], True)

    def test_negative_iteration_max(self):
        toml = fixture("uhfk_u8_negative_max").resolve()
        out = toml.parent / "output"
        self.assertEqual(sorted(p.name for p in out.iterdir()),
                         ["check.dat", "scf_status.dat"])
        self.assertIn("note = IterationMax -3 treated as 0\n",
                      (out / "scf_status.dat").read_text())
        st = read_scf_status(toml)
        self.assertEqual(st.status, "not_converged")
        self.assertEqual(st.iterations, 0)
        self.assertEqual(st.iteration_max, 0)
        self.assertTrue(math.isnan(st.residual))
        self.assertEqual(st.meta["note_kind"], "iteration_max_negative")
        self.assertEqual(st.meta["iteration_max_requested"], -3)

    def test_negative_iteration_max_print_check(self):
        toml = copy_fixture("uhfk_u8_negative_max").resolve()
        self.addCleanup(shutil.rmtree, toml.parent, ignore_errors=True)
        _without(toml, "scf_status.dat")
        st = read_scf_status(toml)
        self.assertEqual(st.status, "not_converged")
        self.assertEqual(st.iterations, 0)
        self.assertEqual(st.iteration_max, 0)
        self.assertEqual(st.meta["iteration_max_requested"], -3)
        self.assertEqual(st.note, "no SCF step performed; IterationMax -3 treated as 0")


class _StatusFile(unittest.TestCase):
    """Each test writes its own status file into a private copy of
    ``uhfk_u8``."""

    @classmethod
    def setUpClass(cls):
        cls.toml = copy_fixture("uhfk_u8").resolve()
        cls.path = cls.toml.parent / "output" / "scf_status.dat"

    def read(self, text):
        self.path.write_text(text)
        return read_scf_status(self.toml)

    def refused(self, text, key):
        self.path.write_text(text)
        with self.assertRaises(UnsupportedLayoutError) as cm:
            read_scf_status(self.toml)
        self.assertRegex(str(cm.exception), names_key(key))
        self.assertEqual(cm.exception.path, self.path)
        return cm.exception


class TestStatusFileValues(_StatusFile):
    def test_good(self):
        self.assertEqual(self.read(status_text()).status, "converged")

    def test_false_is_not_converged(self):
        st = self.read(status_text(converged="false", residual="0.5"))
        self.assertEqual(st.status, "not_converged")
        self.assertIs(st.converged, False)
        self.assertEqual(st.residual, 0.5)

    def test_missing_keys(self):
        for key, _ in GOOD:
            with self.subTest(key=key):
                self.refused(status_text(drop=(key,)), key)

    def test_duplicate_key(self):
        err = self.refused(status_text(extra=(("iterations", "1"),)), "iterations")
        self.assertIn("line 10", str(err))

    def test_unknown_key(self):
        self.refused(status_text(extra=(("comment", "x"),)), "comment")

    def test_other_version(self):
        for version in ("2", "1.0", "", "one"):
            with self.subTest(version=version):
                self.refused(status_text(format_version=version), "format_version")

    def test_other_version_before_other_keys(self):
        # a newer format may have other keys: the version is named first
        text = "format_version = 2\nsomething = new\n"
        self.refused(text, "format_version")

    def test_bad_boolean(self):
        for value in ("True", "yes", "1", ""):
            with self.subTest(value=value):
                self.refused(status_text(converged=value), "converged")

    def test_bad_integers(self):
        for key in ("iterations", "iteration_max"):
            for value in ("1.0", "x", "", "-1"):
                with self.subTest(key=key, value=value):
                    self.refused(status_text(**{key: value}), key)

    def test_bad_solver(self):
        self.refused(status_text(solver="UHFr",
                                 residual_definition="two_frobenius_over_size:uhfr_full"),
                     "solver")
        self.refused(status_text(solver="RPA"), "solver")

    def test_bad_residual_definition(self):
        self.refused(status_text(residual_definition="two_frobenius_over_size:uhfr_full"),
                     "residual_definition")

    def test_bad_eps(self):
        for value in ("nan", "inf", "-inf", "0", "0.0", "-1e-12", "x"):
            with self.subTest(value=value):
                self.refused(status_text(eps=value), "eps")

    def test_bad_residual(self):
        for value in ("nan", "inf", "-1e-3", "x"):
            with self.subTest(value=value):
                self.refused(status_text(converged="false", iterations="7",
                                         residual=value), "residual")


class TestStatusFileInvariants(_StatusFile):
    def test_iterations_above_iteration_max(self):
        self.refused(status_text(iterations="1001"), "iterations")

    def test_converged_needs_residual_below_eps(self):
        self.refused(status_text(residual="1e-12"), "converged")
        self.refused(status_text(residual="0.5"), "converged")

    def test_converged_needs_a_step(self):
        self.refused(status_text(iterations="0", residual="nan",
                                 note="no SCF step performed"), "converged")

    def test_zero_steps_need_nan_residual(self):
        self.refused(status_text(converged="false", iterations="0", residual="0.0",
                                 note="no SCF step performed"), "residual")

    def test_zero_steps_need_a_note(self):
        self.refused(status_text(converged="false", iterations="0", residual="nan"),
                     "note")

    def test_no_step_note_with_steps(self):
        self.refused(status_text(note="no SCF step performed"), "note")


class TestStatusFileNotes(_StatusFile):
    NOT_FINITE = dict(converged="false", iterations="7", residual="nan",
                      note="residual not finite after step 7")

    def test_residual_not_finite_note(self):
        st = self.read(status_text(**self.NOT_FINITE))
        self.assertEqual(st.status, "not_converged")
        self.assertIs(st.converged, False)
        self.assertIsNone(st.residual)
        self.assertEqual(st.iterations, 7)
        self.assertEqual(st.note, "residual not finite after step 7")
        self.assertEqual(st.meta["note_kind"], "residual_not_finite")
        self.assertEqual(st.meta["nonfinite_step"], 7)
        self.assertIs(st.meta["residual_not_finite"], True)

    def test_residual_not_finite_note_with_converged(self):
        self.refused(status_text(**dict(self.NOT_FINITE, converged="true")), "converged")

    def test_residual_not_finite_note_other_step(self):
        self.refused(status_text(**dict(self.NOT_FINITE,
                                        note="residual not finite after step 6")),
                     "note")

    def test_residual_not_finite_note_finite_residual(self):
        self.refused(status_text(**dict(self.NOT_FINITE, residual="0.5")), "note")

    def test_negative_iteration_max_note(self):
        st = self.read(status_text(converged="false", iterations="0", residual="nan",
                                   iteration_max="0",
                                   note="IterationMax -3 treated as 0"))
        self.assertEqual(st.status, "not_converged")
        self.assertIs(st.converged, False)
        self.assertEqual(st.iterations, 0)
        self.assertEqual(st.iteration_max, 0)
        self.assertEqual(st.meta["note_kind"], "iteration_max_negative")
        self.assertEqual(st.meta["iteration_max_requested"], -3)

    def test_negative_iteration_max_note_needs_zero(self):
        self.refused(status_text(converged="false", iterations="0", residual="nan",
                                 iteration_max="5",
                                 note="IterationMax -3 treated as 0"), "iteration_max")
        self.refused(status_text(converged="false", iterations="0", residual="nan",
                                 iteration_max="0",
                                 note="IterationMax 3 treated as 0"), "note")

    def test_unknown_note(self):
        for note in ("converged nicely", "no SCF step performed.",
                     "residual not finite after step x"):
            with self.subTest(note=note):
                self.refused(status_text(note=note), "note")


def _without(toml, *names):
    for name in names:
        (toml.parent / "output" / name).unlink()


class TestPrintCheck(unittest.TestCase):
    """``scf_status.dat`` deleted: the ``print_check`` file decides."""

    def setUp(self):
        self.toml = copy_fixture("uhfk_u8").resolve()
        self.addCleanup(shutil.rmtree, self.toml.parent, ignore_errors=True)
        _without(self.toml, "scf_status.dat")
        self.check = self.toml.parent / "output" / "check.dat"

    def write_check(self, *rests):
        self.check.write_text("".join(
            "{}, {}, -7.3896258, 16, -3.555e-13\n".format(i, r)
            for i, r in enumerate(rests)))

    def test_fixture_trace(self):
        st = read_scf_status(self.toml)
        self.assertEqual(st.status, "converged")
        self.assertEqual(st.source_kind, "print_check")
        self.assertEqual(st.source, self.check)
        self.assertEqual(st.iterations, 1)
        self.assertEqual(st.residual,
                         float(self.check.read_text().split(",")[1]))
        self.assertLess(st.residual, 1e-14)
        self.assertEqual(st.eps, 1e-12)
        self.assertEqual(st.iteration_max, 1000)

    def test_band(self):
        for rest, status in (("9.999998e-13", "converged"),
                             ("1.0000002e-12", "not_converged"),
                             ("1e-12", "unknown"),
                             ("9.9999999e-13", "unknown"),
                             ("1.00000001e-12", "unknown")):
            with self.subTest(rest=rest):
                self.write_check("0.5", "0.25", rest)
                st = read_scf_status(self.toml)
                self.assertEqual(st.status, status)
                self.assertEqual(st.iterations, 3)
                self.assertEqual(st.residual, float(rest))
                self.assertEqual(st.converged,
                                 {"converged": True, "not_converged": False,
                                  "unknown": None}[status])

    def test_nan_rest(self):
        self.write_check("0.5", "nan")
        st = read_scf_status(self.toml)
        self.assertEqual(st.status, "not_converged")
        self.assertIsNone(st.residual)
        self.assertIs(st.meta["residual_not_finite"], True)

    def test_eps_key_case_insensitive(self):
        text = self.toml.read_text()
        self.assertIn("EPS = 12", text)
        self.toml.write_text(text.replace("EPS = 12", "eps = 10"))
        self.write_check("5e-11")
        self.assertEqual(read_scf_status(self.toml).status, "converged")

    def test_disabled_status_file_is_ignored(self):
        # with scf_status = "" a stale status file is not the run's
        (self.toml.parent / "output" / "scf_status.dat").write_text(status_text())
        text = self.toml.read_text().replace(
            "[file.output]\n", '[file.output]\n  scf_status = ""\n')
        self.toml.write_text(text)
        self.write_check("0.5")
        st = read_scf_status(self.toml)
        self.assertEqual(st.source_kind, "print_check")
        self.assertEqual(st.status, "not_converged")

    def test_malformed(self):
        for text, line in (("0, 0.5, -7, 16, 0\n1, 0.25, -7, 16\n", 2),
                           ("0, x, -7, 16, 0\n", 1),
                           ("0, 0.5, -7, 16, 0\n\n1, 0.5, -7, 16, 0\n", 2)):
            with self.subTest(text=text):
                self.check.write_text(text)
                with self.assertRaises(UnsupportedLayoutError) as cm:
                    read_scf_status(self.toml)
                self.assertIn("line {}".format(line), str(cm.exception))
                self.assertEqual(cm.exception.path, self.check)

    def test_uhfr_complex_energy_column(self):
        # UHFr's energy is complex: its column is complex text, written with
        # the format string of uhfr.py:609-614; only Rest is read
        line = "{}, {:.8g}, {:.8g}, {:.4g}, {:.4g}\n"
        energy = complex(-3.4081901344465244, 0.0)
        text = (line.format(0, 0.5, energy, 8.0, -2.2e-16)
                + line.format(1, 5e-13, energy, 8.0, -2.2e-16))
        self.assertIn("-3.4081901+0j", text)
        self.check.write_text(text)
        self.assertEqual(read_scf_status(self.toml).status, "converged")

    def test_incomplete_last_line(self):
        # a hard kill can cut the last line (the trace is buffered)
        for text, steps in (("0, 0.5, -7, 16, 0\n1, 0.25, -7, 16, 0\n2, 0.1", 2),
                            ("0, 0.5, -7, 16, 0\n1, 1e-13, -7, 16, 0", 1),
                            ("0, 0.5, -7", 0)):
            with self.subTest(text=text):
                self.check.write_text(text)
                st = read_scf_status(self.toml)
                self.assertEqual(st.status, "not_converged")
                self.assertIs(st.converged, False)
                self.assertEqual(st.iterations, steps)
                incomplete = steps + 1
                self.assertEqual(st.meta["incomplete_line"], incomplete)
                self.assertIn("incomplete line {}".format(incomplete), st.note)
                if steps:
                    self.assertEqual(st.residual, float(
                        text.splitlines()[steps - 1].split(",")[1]))
                else:
                    self.assertTrue(math.isnan(st.residual))

    def test_negative_iteration_max_config(self):
        text = self.toml.read_text()
        self.assertIn("IterationMax = 1000", text)
        self.toml.write_text(text.replace("IterationMax = 1000", "IterationMax = -3"))
        self.check.write_text("")
        st = read_scf_status(self.toml)
        self.assertEqual(st.status, "not_converged")
        self.assertEqual(st.iteration_max, 0)
        self.assertEqual(st.meta["iteration_max_requested"], -3)
        self.assertIn("IterationMax -3 treated as 0", st.note)

    def test_energy_missing_is_failed(self):
        _without(self.toml, "energy.dat")
        st = read_scf_status(self.toml)
        self.assertEqual(st.status, "failed")
        self.assertIsNone(st.converged)
        self.assertEqual(st.source_kind, "print_check")
        self.assertEqual(st.iterations, 1)
        self.assertEqual(st.meta["scf_verdict"], "converged")
        self.assertIn("energy.dat", st.note)


class TestEmptyPrintCheck(unittest.TestCase):
    def test_zero_steps(self):
        toml = copy_fixture("uhfk_u8_zero_steps").resolve()
        self.addCleanup(shutil.rmtree, toml.parent, ignore_errors=True)
        _without(toml, "scf_status.dat")
        st = read_scf_status(toml)
        self.assertEqual(st.status, "not_converged")
        self.assertIs(st.converged, False)
        self.assertEqual(st.iterations, 0)
        self.assertTrue(math.isnan(st.residual))
        self.assertEqual(st.note, "no SCF step performed")
        self.assertEqual(st.source_kind, "print_check")


class TestNoRecords(unittest.TestCase):
    def setUp(self):
        self.toml = copy_fixture("uhfk_u8").resolve()
        self.addCleanup(shutil.rmtree, self.toml.parent, ignore_errors=True)
        self.out = self.toml.parent / "output"

    def test_empty_output_directory(self):
        for path in self.out.iterdir():
            path.unlink()
        st = read_scf_status(self.toml)
        self.assertEqual(st.status, "not_run")
        self.assertIsNone(st.converged)
        self.assertEqual(st.source_kind, "none")
        self.assertEqual(st.note, NOT_RUN_NOTE)
        self.assertIn("qlms.run creates the output directory", st.note)

    def test_no_output_directory(self):
        shutil.rmtree(self.out)
        self.assertEqual(read_scf_status(self.toml).status, "not_run")

    def test_energy_missing_green_present(self):
        _without(self.toml, "scf_status.dat", "check.dat", "energy.dat")
        self.assertTrue((self.out / "green.dat.npz").is_file())
        st = read_scf_status(self.toml)
        self.assertEqual(st.status, "failed")
        self.assertIsNone(st.converged)
        self.assertEqual(st.source_kind, "none")
        self.assertIn("energy.dat", st.note)

    def test_unsupported_mode(self):
        self.toml.write_text(self.toml.read_text().replace('"UHFk"', '"RPA"'))
        with self.assertRaises(UnsupportedLayoutError):
            read_scf_status(self.toml)


class TestEnergy(unittest.TestCase):
    def test_uhfk(self):
        toml = fixture("uhfk_u8").resolve()
        e = read_energy(toml)
        self.assertIsInstance(e, Energy)
        self.assertEqual(e.solver, "UHFk")
        self.assertEqual(e.source, toml.parent / "output" / "energy.dat")
        self.assertEqual(e.values, file_numbers(e.source))
        values = e.values
        for field, key in (("total", "Energy_Total"), ("band", "Energy_Band"),
                           ("free_total", "FreeEnergy_Total"),
                           ("free_band", "FreeEnergy_Band"), ("ncond", "NCond"),
                           ("sz", "Sz"), ("mu", "ChemicalPotential"), ("sx", "Sx"),
                           ("sy", "Sy")):
            self.assertEqual(getattr(e, field), values[key], field)
        # the golden U = 8 Neel state (tests/hf_square_oracle.py, (4, 8))
        self.assertAlmostEqual(e.total, -7.38962582304832, delta=1e-10)
        self.assertAlmostEqual(e.free_total, e.total, delta=1e-12)
        self.assertAlmostEqual(e.ncond, 16.0, delta=1e-10)
        self.assertAlmostEqual(e.mu, 4.0, delta=1e-10)
        self.assertAlmostEqual(e.total, e.band + values["Energy_CoulombIntra"],
                               delta=1e-12)
        self.assertEqual(list(e.values)[:3],
                         ["Energy_Total", "Energy_Band", "Energy_CoulombIntra"])
        self.assertEqual(e.meta["missing_keys"], ())
        assert_deep_frozen(e)
        self.assertIsNot(e.values, read_energy(toml).values)

    def test_uhfr_block_chemical_potentials(self):
        e = read_energy(fixture("uhfr_tutorial"))
        self.assertEqual(e.solver, "UHFr")
        values = file_numbers(e.source)
        self.assertEqual(e.values, values)
        for field, key in (("total", "Energy_total"), ("band", "Energy_band"),
                           ("free_total", "FreeEnergy_total"),
                           ("free_band", "FreeEnergy_band"), ("ncond", "NCond"),
                           ("sz", "Sz")):
            self.assertEqual(getattr(e, field), values[key], field)
        self.assertEqual(e.mu, {"spin-up": values["ChemicalPotential_spin-up"],
                                "spin-down": values["ChemicalPotential_spin-down"]})
        self.assertAlmostEqual(e.total, -3.40819013444652, delta=1e-9)
        self.assertAlmostEqual(e.ncond, 8.0, delta=1e-10)
        for mu in e.mu.values():
            self.assertAlmostEqual(mu, 4.0, delta=1e-6)
        self.assertEqual((e.sx, e.sy), (0.0, 0.0))
        self.assertEqual(e.meta["missing_keys"], ())

    def test_legacy(self):
        e = read_energy(fixture("uhfk_cdw_legacy"))
        self.assertEqual(e.total, -2.2661791934436106)
        self.assertEqual(e.sz, -4.562322741819003e-16)
        self.assertIsNone(e.mu)
        self.assertIsNone(e.sx)
        self.assertIsNone(e.sy)
        self.assertIsNone(e.free_total)
        self.assertEqual(e.meta["missing_keys"],
                         ("free_total", "free_band", "mu", "sx", "sy"))


class TestEnergyFile(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.toml = copy_fixture("uhfk_u8").resolve()
        cls.path = cls.toml.parent / "output" / "energy.dat"
        cls.good = cls.path.read_text()

    def refused(self, text, key):
        self.path.write_text(text)
        with self.assertRaises(UnsupportedLayoutError) as cm:
            read_energy(self.toml)
        self.assertRegex(str(cm.exception), names_key(key))
        return cm.exception

    def test_bad_value(self):
        # the NCond digits differ between platforms: replace the whole line
        text, count = re.subn(r"(?m)^NCond\s*=.*$", "NCond = sixteen", self.good)
        self.assertEqual(count, 1)
        self.refused(text, "NCond")

    def test_duplicate_key(self):
        self.refused(self.good + "Sz = 0.0\n", "Sz")
        self.refused(self.good + "energy_total = 0.0\n", "energy_total")

    def test_plain_and_suffixed_chemical_potential(self):
        self.refused(self.good + "ChemicalPotential_0 = 1.0\n", "ChemicalPotential_0")

    def test_missing_file(self):
        self.path.unlink()
        try:
            with self.assertRaises(MissingOutputError):
                read_energy(self.toml)
        finally:
            self.path.write_text(self.good)

    def tearDown(self):
        self.path.write_text(self.good)


if __name__ == "__main__":
    unittest.main()
