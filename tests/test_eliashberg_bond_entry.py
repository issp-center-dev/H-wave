"""Entries of the bond-resolved dynamic pairing kernel (spec 10.2)."""
import os
import shutil
import tempfile
import unittest

import numpy as np

from tests.test_flex_bond_gate import _flex


def _run_gate(tmp, output_full=True, extra=None, inter=None):
    """A tiny real bond-gate FLEX solve (4x4, Nmat 8, two orbitals, CoulombInter) writing every
    artifact into tmp; returns (solver, green_info)."""
    par = {"longitudinal_bond_output_full": output_full}
    par.update(extra or {})
    s, r = _flex(par, inter=inter)
    gi = r.get_param("green")
    s.solve(gi, tmp)
    outputs = {"path_to_output": tmp, "chiq_s": "chiq_s", "chiq_c": "chiq_c", "green": "green",
               "sigma": "sigma", "longitudinal_bond": "longitudinal_bond"}
    s.validate_output_paths(outputs, tmp)
    s.save_results(outputs, gi)
    return s, gi


class TestArchiveSchema2(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_archive_carries_vertices(self):
        s, gi = _run_gate(self.tmp)
        with np.load(os.path.join(self.tmp, "longitudinal_bond.npz")) as d:
            self.assertEqual(int(d["bond_archive_schema"]), 2)
            self.assertEqual(int(d["norb"]), s.norb)
            np.testing.assert_array_equal(d["S_bond"], s._bond_S)
            np.testing.assert_array_equal(d["C_bond"], s._bond_C)
            self.assertNotIn("longitudinal_bond_S", d.files)
        np.testing.assert_array_equal(gi["longitudinal_bond_S"], s._bond_S)

    def test_no_vertices_published_without_output_full(self):
        s, gi = _run_gate(self.tmp, output_full=False)
        self.assertNotIn("longitudinal_bond_S", gi)
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "longitudinal_bond.npz")))


class TestLoadBondArchive(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.s, self.gi = _run_gate(self.tmp)
        self.path = os.path.join(self.tmp, "longitudinal_bond.npz")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _load(self, **kw):
        from hwave.solver.eliashberg_bond_io import load_bond_archive
        args = dict(norb=self.s.norb, nmat_expected=self.s.nmat,
                    cell_shape_expected=tuple(int(x) for x in self.s.lattice.shape),
                    beta_expected=1.0 / self.s.T)
        args.update(kw)
        return load_bond_archive(self.path, **args)

    def test_loads_and_members_on_demand(self):
        a = self._load()
        self.assertEqual(a.ND, a.B * self.s.norb ** 2)
        np.testing.assert_array_equal(a.S_bond, self.s._bond_S)
        self.assertEqual(a.chi_shape, (self.s.nmat, self.s.lattice.nvol, a.ND, a.ND))
        chi = a.member("chi_s_w")
        np.testing.assert_array_equal(chi, self.gi["longitudinal_bond_chi_s_w"])

    def test_header_without_loading(self):
        from hwave.solver.eliashberg_bond_io import _npz_member_header
        shape, dtype = _npz_member_header(self.path, "chi_c_w")
        self.assertEqual(shape, (self.s.nmat, self.s.lattice.nvol, self.s._bond_S.shape[-1], self.s._bond_S.shape[-1]))
        self.assertEqual(np.dtype(dtype), np.dtype(np.complex128))

    def _rewrite(self, **changes):
        with np.load(self.path) as d:
            members = {k: d[k] for k in d.files}
        members.update(changes)
        np.savez(self.path, **members)

    def test_schema1_refused(self):                                # 10.2.5
        self._rewrite(bond_archive_schema=np.int64(1))
        with self.assertRaisesRegex(ValueError, "longitudinal_bond_output_full = true"):
            self._load()
        self._rewrite(bond_archive_schema=np.int64(3))
        with self.assertRaisesRegex(ValueError, "unsupported archive schema"):
            self._load()

    def test_validations(self):
        with self.assertRaisesRegex(ValueError, "nmat"):
            self._load(nmat_expected=self.s.nmat + 2)
        with self.assertRaisesRegex(ValueError, "cell_shape"):
            self._load(cell_shape_expected=(2, 2, 1))
        with self.assertRaisesRegex(ValueError, "beta"):
            self._load(beta_expected=3.3)
        with self.assertRaisesRegex(ValueError, "norb"):
            self._load(norb=self.s.norb + 1)
        self._rewrite(index_order=np.str_("I = something else"))
        with self.assertRaisesRegex(ValueError, "index_order"):
            self._load()
        self._rewrite(index_order=np.str_("I = m*norb**2 + l1*norb + l2"),
                      reverse=np.array([0, 2, 1, 3][:len(np.load(self.path)["reverse"])]))
        with self.assertRaisesRegex(ValueError, "reverse"):
            self._load()

    def test_delta_r_aliasing_refused(self):
        with np.load(self.path) as d:
            dr = d["delta_r"].copy()
        dr[1] = dr[1] + np.array([4, 0, 0])      # aliases modulo the 4x4x1 cell
        self._rewrite(delta_r=dr)
        with self.assertRaisesRegex(ValueError, "modulo"):
            self._load()


class TestPairingControls(unittest.TestCase):
    def test_defaults_and_validation(self):
        from hwave.solver.eliashberg_bond_io import PairingControls
        c = PairingControls.from_param({}, pairing_types=("singlet",))
        self.assertEqual((c.solver_mode, c.eigenvalue_method, c.num_eigenvalues, c.max_iter, c.alpha,
                          c.convergence_tol, c.matsubara_basis, c.ir_tol, c.ir_fit_tol, c.fft_workers),
                         ("iteration", "arnoldi", 10, 1000, 0.5, 1e-5, "uniform", 1e-8, 0.1, 1))
        self.assertIsNone(c.init_gap); self.assertIsNone(c.bond_memory_cap_gb); self.assertFalse(c.ir_keep_static_chi)
        c2 = PairingControls.from_param({"matsubara_basis": "IR", "ir_fit_tol": 0, "num_eigenvalues": 3,
                                         "ir_keep_static_chi": "true"}, pairing_types=("singlet", "triplet"))
        self.assertEqual((c2.matsubara_basis, c2.ir_fit_tol, c2.num_eigenvalues, c2.ir_keep_static_chi),
                         ("ir", 0.0, 3, True))
        for bad in ({"matsubara_basis": "x"}, {"ir_fit_tol": -1}, {"num_eigenvalues": 0},
                    {"solver_mode": "z"}, {"bond_memory_cap_gb": -2}):
            with self.assertRaises(ValueError):
                PairingControls.from_param(bad, pairing_types=("singlet",))
        with self.assertRaises(ValueError):
            PairingControls.from_param({}, pairing_types=("x",))
