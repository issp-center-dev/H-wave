"""Tests of ``hwave.results._io``: npz, text and geometry access, the legacy
table and the input-file paths."""
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

from hwave.results import (
    MissingOutputError, UnsupportedLayoutError, open_run)
from hwave.results._io import (
    MARKER_KINDS, Geometry, check_markers, input_file_path, member_shapes,
    open_npz, read_geometry, read_key_value_text, scalar_members)
from tests.results_fixtures import REPO, assert_deep_frozen, fixture


class _TmpDir(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def write(self, name, text):
        path = self.tmp / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path


def _fail(*args, **kwargs):
    raise AssertionError("an array member was loaded")


class TestKeyValueText(_TmpDir):
    def test_padded_and_unpadded(self):
        path = self.write("energy.dat",
                          "Energy_Total        = -2.25\n"
                          "NCond=16.0\n"
                          "\n"
                          "   \n"
                          "  note =  \n"
                          "residual_definition = a=b\n")
        self.assertEqual(read_key_value_text(path), [
            ("Energy_Total", "-2.25", 1),
            ("NCond", "16.0", 2),
            ("note", "", 5),
            ("residual_definition", "a=b", 6),
        ])

    def test_duplicates_are_kept_in_order(self):
        path = self.write("s.dat", "a = 1\nb = 2\na = 3\n")
        self.assertEqual(read_key_value_text(path),
                         [("a", "1", 1), ("b", "2", 2), ("a", "3", 3)])

    def test_line_without_equals(self):
        path = self.write("s.dat", "a = 1\njust text\n")
        with self.assertRaises(UnsupportedLayoutError) as cm:
            read_key_value_text(path)
        self.assertIn("line 2", str(cm.exception))
        self.assertEqual(cm.exception.path, path)

    def test_empty_key(self):
        path = self.write("s.dat", " = 1\n")
        with self.assertRaises(UnsupportedLayoutError):
            read_key_value_text(path)

    def test_missing_file(self):
        with self.assertRaises(MissingOutputError):
            read_key_value_text(self.tmp / "absent.dat")


class TestNpz(_TmpDir):
    def setUp(self):
        super().setUp()
        self.path = self.tmp / "a.npz"
        np.savez(self.path, big=np.zeros((3, 4, 5), dtype=complex),
                 marker=np.array("e_plus_ikR"), plain="text", count=np.int64(7),
                 x=np.float64(0.5), flag=np.bool_(True), raw=np.array(b"bytes"),
                 vec=np.arange(2))

    def test_member_shapes_header_only(self):
        with mock.patch.object(np.lib.npyio.NpzFile, "__getitem__", _fail), \
                mock.patch.object(np.lib.format, "read_array", _fail), \
                mock.patch("numpy.load", _fail):
            shapes = member_shapes(self.path)
        self.assertEqual(shapes, {"big": (3, 4, 5), "marker": (), "plain": (),
                                  "count": (), "x": (), "flag": (), "raw": (),
                                  "vec": (2,)})

    def test_scalar_members(self):
        with mock.patch.object(np.lib.npyio.NpzFile, "__getitem__", _fail), \
                mock.patch("numpy.load", _fail):
            scalars = scalar_members(self.path)
        self.assertEqual(scalars, {"marker": "e_plus_ikR", "plain": "text",
                                   "count": 7, "x": 0.5, "flag": True,
                                   "raw": "bytes"})
        self.assertIs(type(scalars["marker"]), str)
        self.assertIs(type(scalars["count"]), int)
        self.assertIs(type(scalars["x"]), float)
        self.assertIs(type(scalars["flag"]), bool)

    def test_open_npz_closes(self):
        with open_npz(self.path) as data:
            self.assertEqual(data["vec"].tolist(), [0, 1])
            handle = data.zip.fp
        self.assertIsNone(data.zip)
        self.assertTrue(handle is None or handle.closed)

    def test_open_npz_closes_on_error(self):
        with self.assertRaises(KeyError):
            with open_npz(self.path) as data:
                data["no_such_member"]
        self.assertIsNone(data.zip)

    def test_open_npz_no_pickle(self):
        path = self.tmp / "obj.npz"
        np.savez(path, obj=np.array([{"a": 1}], dtype=object))
        with open_npz(path) as data:
            with self.assertRaises(ValueError):
                data["obj"]
        with self.assertRaises(UnsupportedLayoutError):
            scalar_members(self._object_scalar())

    def _object_scalar(self):
        path = self.tmp / "obj0.npz"
        np.savez(path, obj=np.array({"a": 1}, dtype=object))
        return path

    def test_not_an_npz(self):
        bad = self.write("bad.npz", "not a zip file\n")
        for func in (member_shapes, scalar_members):
            with self.assertRaises(UnsupportedLayoutError):
                func(bad)
        with self.assertRaises(UnsupportedLayoutError):
            with open_npz(bad):
                pass

    def test_missing_npz(self):
        for func in (member_shapes, scalar_members):
            with self.assertRaises(MissingOutputError):
                func(self.tmp / "absent.npz")
        with self.assertRaises(MissingOutputError):
            with open_npz(self.tmp / "absent.npz"):
                pass


class TestGeometry(_TmpDir):
    def test_tests_uhfk_coulombintra(self):
        path = REPO / "tests" / "uhfk" / "CoulombIntra" / "geom.dat"
        geom = read_geometry(path)
        self.assertIsInstance(geom, Geometry)
        np.testing.assert_array_equal(geom.lattice_vectors, np.eye(3))
        self.assertEqual(geom.centres.shape, (1, 3))
        np.testing.assert_array_equal(geom.centres, np.zeros((1, 3)))
        self.assertEqual(geom.norb, 1)
        self.assertEqual(geom.source, path)
        assert_deep_frozen(geom)

    def test_two_orbitals(self):
        path = self.write("geom.dat",
                          "2 0 0\n0 1 0\n0 0 1\n2\n0.25 0.5 0\n0.75 0.5 0.0\n\n")
        geom = read_geometry(path)
        np.testing.assert_array_equal(geom.lattice_vectors, np.diag([2.0, 1, 1]))
        np.testing.assert_array_equal(geom.centres,
                                      [[0.25, 0.5, 0], [0.75, 0.5, 0]])

    def test_too_few_centres(self):
        path = self.write("geom.dat", "1 0 0\n0 1 0\n0 0 1\n2\n0 0 0\n")
        with self.assertRaises(UnsupportedLayoutError) as cm:
            read_geometry(path)
        self.assertIn("2 orbital", str(cm.exception))

    def test_extra_lines(self):
        path = self.write("geom.dat", "1 0 0\n0 1 0\n0 0 1\n1\n0 0 0\n0 0 0\n")
        with self.assertRaises(UnsupportedLayoutError):
            read_geometry(path)

    def test_bad_values(self):
        for text in ("1 0\n0 1 0\n0 0 1\n1\n0 0 0\n",
                     "1 0 0\n0 1 0\n0 0 1\nx\n0 0 0\n",
                     "1 0 0\n0 1 0\n0 0 1\n0\n",
                     "1 0 0\n0 1 0\n0 0 1\n1\n0 a 0\n"):
            path = self.write("geom.dat", text)
            with self.assertRaises(UnsupportedLayoutError):
                read_geometry(path)

    def test_missing(self):
        with self.assertRaises(MissingOutputError):
            read_geometry(self.tmp / "geom.dat")


class TestLegacyTable(unittest.TestCase):
    """The per-marker legacy table on the git-tracked tutorial outputs."""

    @classmethod
    def setUpClass(cls):
        cls.output = fixture("uhfk_cdw_legacy").parent / "output"

    def test_fixture_is_legacy(self):
        for name in ("green.npz", "eigen.npz"):
            shapes = member_shapes(self.output / name)
            self.assertNotIn("momentum_convention", shapes)
            self.assertNotIn("green_convention", shapes)
        self.assertIn("green_sublattice", member_shapes(self.output / "green.npz"))

    def test_kinds(self):
        self.assertEqual(set(MARKER_KINDS),
                         {"uhfk_eigen", "uhfk_green", "uhfk_green_sublattice"})
        with self.assertRaises(ValueError):
            check_markers(self.output / "eigen.npz", "energy", legacy=True)

    def test_refused_without_legacy(self):
        for name, kind in (("eigen.npz", "uhfk_eigen"),
                           ("green.npz", "uhfk_green_sublattice"),
                           ("green.npz", "uhfk_green")):
            with self.assertRaises(UnsupportedLayoutError) as cm:
                check_markers(self.output / name, kind, legacy=False)
            self.assertIn("momentum_convention", str(cm.exception))
            self.assertIn("legacy=True", cm.exception.hint)
            self.assertEqual(cm.exception.path, self.output / name)

    def test_legacy_eigen(self):
        notes = check_markers(self.output / "eigen.npz", "uhfk_eigen", legacy=True)
        self.assertEqual(len(notes), 1)
        self.assertEqual(notes[0]["member"], "momentum_convention")
        self.assertEqual(notes[0]["assumed"], "e_plus_ikR")
        self.assertIn("UHFk", notes[0]["reason"])

    def test_legacy_green_refuses_deflated_green(self):
        notes = check_markers(self.output / "green.npz", "uhfk_green_sublattice",
                              legacy=True)
        by_member = {note["member"]: note for note in notes}
        self.assertEqual(set(by_member), {"momentum_convention", "green"})
        self.assertEqual(by_member["momentum_convention"]["assumed"], "e_plus_ikR")
        self.assertIsNone(by_member["green"]["assumed"])
        self.assertIn("green_sublattice", by_member["green"]["reason"])
        self.assertIn("green_convention", by_member["green"]["reason"])


class TestMarkers(_TmpDir):
    """The marker rows on synthetic files."""

    def npz(self, **members):
        path = self.tmp / "f.npz"
        np.savez(path, **members)
        return path

    def test_current_files_need_no_assumption(self):
        g = np.zeros((1, 2, 1, 2, 1))
        eigen = self.npz(eigenvalue=np.zeros((1, 2)), momentum_convention="e_plus_ikR")
        self.assertEqual(check_markers(eigen, "uhfk_eigen", legacy=False), [])
        plain = self.npz(green=g, momentum_convention="e_plus_ikR")
        self.assertEqual(check_markers(plain, "uhfk_green", legacy=False), [])
        sub = self.npz(green=g, green_sublattice=g,
                       green_convention=np.array("green_slot_first"),
                       momentum_convention="e_plus_ikR")
        self.assertEqual(check_markers(sub, "uhfk_green_sublattice", legacy=False), [])

    def test_wrong_momentum_marker_is_refused(self):
        for value in ("e_minus_ikR", np.array(["e_plus_ikR", "e_plus_ikR"])):
            path = self.npz(green=np.zeros(1), momentum_convention=value)
            for legacy in (False, True):
                with self.assertRaises(UnsupportedLayoutError) as cm:
                    check_markers(path, "uhfk_green", legacy=legacy)
                self.assertNotIn("legacy=True", str(cm.exception))

    def test_wrong_green_convention_is_refused(self):
        g = np.zeros((1, 2, 1, 2, 1))
        path = self.npz(green=g, green_sublattice=g, green_convention="other",
                        momentum_convention="e_plus_ikR")
        with self.assertRaises(UnsupportedLayoutError) as cm:
            check_markers(path, "uhfk_green_sublattice", legacy=True)
        self.assertNotIn("legacy=True", str(cm.exception))

    def test_sublattice_without_green_sublattice_is_never_assumed(self):
        path = self.npz(green=np.zeros((1, 2, 1, 2, 1)))
        for legacy in (False, True):
            with self.assertRaises(UnsupportedLayoutError) as cm:
                check_markers(path, "uhfk_green_sublattice", legacy=legacy)
            self.assertIn("green_sublattice", str(cm.exception))
            self.assertNotIn("legacy=True", str(cm.exception))

    def test_marked_sublattice_file_without_green_convention(self):
        # momentum marker present, fold marker absent: the deflated green is
        # refused also without legacy=True (nothing is assumed)
        g = np.zeros((1, 2, 1, 2, 1))
        path = self.npz(green=g, green_sublattice=g, momentum_convention="e_plus_ikR")
        notes = check_markers(path, "uhfk_green_sublattice", legacy=False)
        self.assertEqual([n["member"] for n in notes], ["green"])
        self.assertIsNone(notes[0]["assumed"])

    def test_missing_member_payload(self):
        path = self.npz(momentum_convention="e_plus_ikR")
        with self.assertRaises(UnsupportedLayoutError):
            check_markers(path, "uhfk_eigen", legacy=False)
        with self.assertRaises(UnsupportedLayoutError):
            check_markers(path, "uhfk_green", legacy=False)

    def test_no_array_is_loaded(self):
        g = np.zeros((1, 2, 1, 2, 1))
        path = self.npz(green=g, green_sublattice=g, momentum_convention="e_plus_ikR",
                        green_convention="green_slot_first")
        with mock.patch.object(np.lib.npyio.NpzFile, "__getitem__", _fail), \
                mock.patch("numpy.load", _fail):
            self.assertEqual(check_markers(path, "uhfk_green_sublattice",
                                           legacy=False), [])


class TestInputFilePath(_TmpDir):
    def run_of(self, text):
        return open_run(self.write("input.toml", text))

    def test_uhfk_interaction_dir_case_insensitive(self):
        self.write("model/geom.dat", "x\n")
        run = self.run_of('[mode]\nmode = "UHFk"\n'
                          '[file.input]\npath_to_input = "ignored"\n'
                          '[file.input.interaction]\nPath_To_Input = "model"\n'
                          'GEOMETRY = "geom.dat"\nTransfer = "transfer.dat"\n')
        self.assertEqual(input_file_path(run, "Geometry"),
                         self.tmp.resolve() / "model" / "geom.dat")
        self.assertEqual(input_file_path(run, "geometry"),
                         input_file_path(run, "GeoMetry"))
        # existence is checked only when asked
        self.assertEqual(input_file_path(run, "transfer", must_exist=False),
                         self.tmp.resolve() / "model" / "transfer.dat")
        with self.assertRaises(MissingOutputError) as cm:
            input_file_path(run, "transfer")
        self.assertIn("transfer.dat", str(cm.exception))

    def test_uhfr_uses_file_input_dir(self):
        self.write("defs/trans.def", "x\n")
        run = self.run_of('[mode]\nmode = "UHFr"\n'
                          '[file.input]\npath_to_input = "defs"\n'
                          '[file.input.interaction]\nTRANS = "trans.def"\n')
        self.assertEqual(input_file_path(run, "Trans"),
                         self.tmp.resolve() / "defs" / "trans.def")

    def test_uhfr_dir_key_is_exact_case(self):
        # qlms.py:108 reads [file.input] path_to_input with its exact case
        self.write("trans.def", "x\n")
        run = self.run_of('[mode]\nmode = "UHFr"\n'
                          '[file.input]\nPath_To_Input = "defs"\n'
                          '[file.input.interaction]\nTrans = "trans.def"\n')
        self.assertEqual(input_file_path(run, "trans"),
                         self.tmp.resolve() / "trans.def")

    def test_base_dir(self):
        other = self.tmp / "elsewhere"
        (other / "model").mkdir(parents=True)
        (other / "model" / "geom.dat").write_text("x\n")
        toml = self.write("input.toml",
                          '[mode]\nmode = "UHFk"\n[file.input.interaction]\n'
                          'path_to_input = "model"\nGeometry = "geom.dat"\n')
        run = open_run(toml, base_dir=other)
        self.assertEqual(input_file_path(run, "Geometry"),
                         other.resolve() / "model" / "geom.dat")

    def test_absolute_name(self):
        target = self.write("abs/geom.dat", "x\n").resolve()
        run = self.run_of('[mode]\nmode = "UHFk"\n[file.input.interaction]\n'
                          'path_to_input = "model"\nGeometry = "{}"\n'.format(target))
        self.assertEqual(input_file_path(run, "Geometry"), target)

    def test_not_configured(self):
        run = self.run_of('[mode]\nmode = "UHFk"\n[file.input.interaction]\n'
                          'Geometry = "geom.dat"\n')
        with self.assertRaises(MissingOutputError) as cm:
            input_file_path(run, "CoulombIntra")
        self.assertIn("CoulombIntra", str(cm.exception))
        self.assertIn("[file.input.interaction]", cm.exception.hint)
        with self.assertRaises(MissingOutputError):
            input_file_path(run, "path_to_input")

    def test_non_string_name(self):
        run = self.run_of('[mode]\nmode = "UHFk"\n[file.input.interaction]\n'
                          'Geometry = 3\n')
        with self.assertRaises(UnsupportedLayoutError):
            input_file_path(run, "Geometry")

    def test_interaction_not_a_table(self):
        run = self.run_of('[mode]\nmode = "UHFk"\n[file.input]\ninteraction = "x"\n')
        with self.assertRaises(UnsupportedLayoutError) as cm:
            input_file_path(run, "Geometry")
        self.assertIn("[file.input.interaction]", str(cm.exception))

    def test_table_case_follows_solver(self):
        # qlms.run reads [file] and [file.input] with their exact case
        run = self.run_of('[mode]\nmode = "UHFk"\n[File.Input.Interaction]\n'
                          'Geometry = "geom.dat"\n')
        with self.assertRaises(MissingOutputError):
            input_file_path(run, "Geometry", must_exist=False)

    def test_unsupported_mode(self):
        run = self.run_of('[mode]\nmode = "Other"\n[file.input.interaction]\n'
                          'Geometry = "geom.dat"\n')
        with self.assertRaises(UnsupportedLayoutError):
            input_file_path(run, "Geometry")

    def test_tutorial_geometry(self):
        run = open_run(fixture("uhfk_cdw_legacy"))
        path = input_file_path(run, "geometry")
        self.assertEqual(path.name, "geom.dat")
        self.assertEqual(read_geometry(path).norb, 1)


if __name__ == "__main__":
    unittest.main()
