"""hwave_sc in non-spin-orbital mode refuses a spin-orbital archive (#271).

A spin-orbital susceptibility archive over N physical orbitals has nd = 2N,
i.e. the shape of a non-spin-orbital general archive over 2N orbitals. Read
as the latter it would be solved with the wrong index interpretation and no
error, so a non-spin-orbital run must refuse any archive that carries the
spin-orbital marker (``spin_orbital = True``, or a FLEX
``chi_convention = "spin_orbital"``) before any heavy work. Archives
without the marker keep loading as before.

The archives are REAL ones: a one-orbital spin-orbital RPA / FLEX run
through ``qlms.run`` (nd = 2, the shape of a two-orbital model).
"""

import os
import shutil
import tempfile
import unittest
from unittest import mock

import numpy as np

import hwave.sc as sc
from hwave import qlms
from hwave.solver import eliashberg_dynamic
from tests.test_sc_spin_orbital import (_input, _quiet, _run_flex,
                                        _write_one_orbital, _wfile, _resave)


def _write_two_orbital_non_so(d):
    """The same hopping as the one-orbital spin-orbital model, read as two
    non-spin-orbital orbitals (geometry, transfer copied verbatim)."""
    _wfile(os.path.join(d, "coulombintra.dat"), "CoulombIntra", 2,
           [(0, 0, 0, 1, 1, 2.0, 0.0), (0, 0, 0, 2, 2, 2.0, 0.0)])
    return {"CoulombIntra": "coulombintra.dat"}


def _so_rpa(d, out):
    files = _write_one_orbital(d)
    inp = _input(d, files, True, {}, nmat=32, out=out)
    inp["mode"]["calc_scheme"] = "general"
    inp["file"]["input"]["path_to_input"] = d
    del inp["eliashberg"]
    inp["log"] = {"print_level": 0}
    with _quiet():
        qlms.run(input_dict=inp)
    return os.path.join(d, out)


class _Base(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = tempfile.mkdtemp(prefix="sc_nonso_")
        cls.rpa = _so_rpa(cls.d, "so_rpa")
        cls.flex = _run_flex(cls.d, {"CoulombIntra": "coulombintra.dat"},
                             "so_flex", nmat=32)
        # the non-spin-orbital twin of the same shape lives in its own dir
        cls.n = os.path.join(cls.d, "nso")
        os.makedirs(cls.n)
        for f in ("geom.dat", "transfer.dat"):
            shutil.copy(os.path.join(cls.d, f), cls.n)
        cls.nfiles = _write_two_orbital_non_so(cls.n)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.d, True)

    def _inp(self, eli, out_dir=None, **kw):
        inp = _input(self.n, self.nfiles, False, eli, nmat=32,
                     out=out_dir or self.rpa)
        inp["file"]["input"]["path_to_input"] = self.n
        return inp


class TestMarkerIsReal(_Base):
    def test_archives_carry_the_markers(self):
        with np.load(os.path.join(self.rpa, "chi0q.npz")) as z:
            self.assertTrue(bool(z["spin_orbital"]))
            self.assertEqual(z["chi0q"].shape[-4:], (2, 2, 2, 2))
        with np.load(os.path.join(self.flex, "chiq.npz")) as z:
            self.assertTrue(bool(z["spin_orbital"]))
            self.assertEqual(str(z["chi_convention"]), "spin_orbital")


class TestLoadMode(_Base):
    def _boom(self):
        return (mock.patch.object(sc, "_convert_chi0q_to_ref_format",
                                  side_effect=AssertionError("heavy work")),
                mock.patch.object(sc, "_compute_vertices",
                                  side_effect=AssertionError("heavy work")))

    def test_loader_refuses_and_names_file(self):
        inp = self._inp({"chi0q_mode": "load"})
        with self.assertRaises(ValueError) as cm:
            sc._load_chi0q(inp, norb=2)
        msg = str(cm.exception)
        self.assertIn(os.path.join(self.rpa, "chi0q.npz"), msg)
        self.assertIn("enable_spin_orbital", msg)
        self.assertIn("non-spin-orbital run", msg)

    def test_calc_eliashberg_refuses_before_vertex_work(self):
        inp = self._inp({"chi0q_mode": "load", "solver_mode": "eigenvalue",
                         "num_eigenvalues": 4})
        p1, p2 = self._boom()
        with _quiet(), p1, p2:
            with self.assertRaisesRegex(ValueError, "spin_orbital"):
                sc.calc_eliashberg(inp)

    def test_refusal_precedes_diagonalization_and_green(self):
        inp = self._inp({"chi0q_mode": "load", "solver_mode": "eigenvalue",
                         "num_eigenvalues": 4})
        with _quiet(), \
                mock.patch.object(sc, "_calc_eigenvalues",
                                  side_effect=AssertionError("heavy work")), \
                mock.patch.object(sc, "_build_bond_green",
                                  side_effect=AssertionError("heavy work")):
            with self.assertRaisesRegex(ValueError, r"spin_orbital = True"):
                sc.calc_eliashberg(inp)

    def test_unmarked_archive_still_loads(self):
        d = os.path.join(self.d, "unmarked")
        os.makedirs(d)
        _resave(os.path.join(self.rpa, "chi0q.npz"),
                os.path.join(d, "chi0q.npz"), drop=("spin_orbital",))
        chi0q, _ = sc._load_chi0q(self._inp({"chi0q_mode": "load"}, d), norb=2)
        self.assertEqual(chi0q.shape[-4:], (2, 2, 2, 2))

    def test_false_marker_loads(self):
        d = os.path.join(self.d, "falsemark")
        os.makedirs(d)
        _resave(os.path.join(self.rpa, "chi0q.npz"),
                os.path.join(d, "chi0q.npz"), spin_orbital=np.array(False))
        sc._load_chi0q(self._inp({"chi0q_mode": "load"}, d), norb=2)

    def test_malformed_marker_refused(self):
        for label, val in (("int", np.array(1)), ("vec", np.array([True])),
                           ("str", np.array("True"))):
            with self.subTest(marker=label):
                d = os.path.join(self.d, "bad_" + label)
                os.makedirs(d)
                _resave(os.path.join(self.rpa, "chi0q.npz"),
                        os.path.join(d, "chi0q.npz"), spin_orbital=val)
                with self.assertRaisesRegex(ValueError, "malformed"):
                    sc._load_chi0q(self._inp({"chi0q_mode": "load"}, d),
                                   norb=2)

    def test_so_mode_still_loads_the_marked_archive(self):
        inp = _input(self.d, _write_one_orbital(self.d), True,
                     {"chi0q_mode": "load"}, nmat=32, out="so_rpa")
        chi0q, _ = sc._load_chi0q(inp, norb=2)
        self.assertEqual(chi0q.shape[-4:], (2, 2, 2, 2))


class TestFlexMode(_Base):
    def _flex_inp(self, dynamic=False, chi=None):
        chi = chi or os.path.join(self.flex, "chiq.npz")
        eli = {"chi0q_mode": "flex", "solver_mode": "eigenvalue",
               "num_eigenvalues": 4, "flex_chi_s": chi, "flex_chi_c": chi}
        if dynamic:
            eli["frequency"] = "dynamic"
        inp = self._inp(eli, self.flex)
        inp["file"]["input"]["path_to_flex_output"] = self.flex
        return inp

    def test_raw_reader_refuses(self):
        inp = self._flex_inp()
        with self.assertRaises(ValueError) as cm:
            sc._read_flex_chi_raw(inp)
        self.assertIn("chiq.npz", str(cm.exception))
        self.assertIn("enable_spin_orbital", str(cm.exception))

    def test_static_flex_refuses_before_heavy_work(self):
        inp = self._flex_inp()
        with _quiet(), \
                mock.patch.object(sc, "_compute_vertices_flex",
                                  side_effect=AssertionError("heavy work")), \
                mock.patch.object(sc, "_load_flex_green",
                                  side_effect=AssertionError("heavy work")):
            with self.assertRaisesRegex(ValueError, r"spin_orbital = True"):
                sc.calc_eliashberg(inp)

    def test_static_flex_refusal_precedes_diagonalization(self):
        inp = self._flex_inp()
        with _quiet(), \
                mock.patch.object(sc, "_calc_eigenvalues",
                                  side_effect=AssertionError("heavy work")):
            with self.assertRaisesRegex(ValueError, r"spin_orbital = True"):
                sc.calc_eliashberg(inp)

    def test_dynamic_refusal_precedes_interaction_build(self):
        inp = self._flex_inp(dynamic=True)
        with _quiet(), \
                mock.patch.object(sc, "_build_interaction_k_for_sites",
                                  side_effect=AssertionError("heavy work")):
            with self.assertRaisesRegex(ValueError, r"spin_orbital = True"):
                eliashberg_dynamic.solve_dynamic(inp)

    def test_dynamic_refuses_before_heavy_work(self):
        inp = self._flex_inp(dynamic=True)
        with _quiet(), \
                mock.patch.object(sc, "_load_flex_green",
                                  side_effect=AssertionError("heavy work")):
            with self.assertRaisesRegex(ValueError, r"spin_orbital = True"):
                eliashberg_dynamic.solve_dynamic(inp)

    def test_convention_tag_alone_is_refused(self):
        d = os.path.join(self.d, "tagonly")
        os.makedirs(d)
        p = os.path.join(d, "chiq.npz")
        _resave(os.path.join(self.flex, "chiq.npz"), p,
                drop=("spin_orbital",))
        with self.assertRaisesRegex(ValueError, "chi_convention"):
            sc._read_flex_chi_raw(self._flex_inp(chi=p))

    def test_malformed_marker_refused(self):
        d = os.path.join(self.d, "flexbad")
        os.makedirs(d)
        p = os.path.join(d, "chiq.npz")
        _resave(os.path.join(self.flex, "chiq.npz"), p,
                spin_orbital=np.array([True]))
        with self.assertRaisesRegex(ValueError, "malformed"):
            sc._read_flex_chi_raw(self._flex_inp(chi=p))

    def test_non_so_flex_archive_still_reads(self):
        d = os.path.join(self.d, "nso_flex")
        os.makedirs(d)
        shutil.copy(os.path.join(self.n, "geom.dat"), d)
        shutil.copy(os.path.join(self.n, "transfer.dat"), d)
        files = _write_two_orbital_non_so(d)
        out = _run_flex(d, files, "flex", so=False, nmat=32)
        inp = _input(d, files, False,
                     {"chi0q_mode": "flex", "flex_chi_s": "chiq_s.npz",
                      "flex_chi_c": "chiq_c.npz"}, nmat=32, out="flex")
        inp["file"]["input"]["path_to_flex_output"] = out
        s, c, conv, _ = sc._read_flex_chi_raw(inp)
        self.assertEqual(conv, "myo")


if __name__ == "__main__":
    unittest.main()
