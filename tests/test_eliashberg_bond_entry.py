"""Entries of the bond-resolved dynamic pairing kernel (spec 10.2)."""
import os
import shutil
import tempfile
import unittest
from unittest import mock

import numpy as np

from tests.heavy_tests import heavy
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


def _npz_to_sidecar(index_path):
    """Rewrite a schema-2 npz archive at ``index_path`` into the sidecar
    layout from the SAME data (an independent oracle for the writer): the two
    channel members move to ``.npy`` files next to the index, which then names
    them and carries bond_archive_layout = "sidecar"."""
    with np.load(index_path) as d:
        members = {k: d[k] for k in d.files}
    chi_s = members.pop("chi_s_w")
    chi_c = members.pop("chi_c_w")
    stem = os.path.basename(index_path)[:-4] if index_path.endswith(".npz") \
        else os.path.basename(index_path)
    s_base = "{}_chi_s_w.npy".format(stem)
    c_base = "{}_chi_c_w.npy".format(stem)
    out_dir = os.path.dirname(index_path)
    np.save(os.path.join(out_dir, s_base), chi_s)
    np.save(os.path.join(out_dir, c_base), chi_c)
    members["bond_archive_layout"] = np.str_("sidecar")
    members["chi_s_w_file"] = np.str_(s_base)
    members["chi_c_w_file"] = np.str_(c_base)
    np.savez(index_path, **members)


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

    def test_default_npz_layout_stamps_the_layout_field(self):
        _run_gate(self.tmp)
        with np.load(os.path.join(self.tmp, "longitudinal_bond.npz")) as d:
            self.assertEqual(str(d["bond_archive_layout"]), "npz")
            self.assertIn("chi_s_w", d.files)
            self.assertIn("chi_c_w", d.files)
        # the default layout writes no sidecar .npy files
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "longitudinal_bond_chi_s_w.npy")))

    def test_default_npz_member_order_appends_layout_last(self):
        # every pre-existing member keeps its position; bond_archive_layout is
        # appended LAST, so removing it reproduces exactly the schema-2 order
        _run_gate(self.tmp)
        with np.load(os.path.join(self.tmp, "longitudinal_bond.npz")) as d:
            files = list(d.files)
        self.assertEqual(files.count("bond_archive_layout"), 1)
        self.assertEqual(files[-1], "bond_archive_layout")
        # the leading positional prefix is unchanged (schema, then the two
        # channels, then S_bond ...)
        self.assertEqual(files[:5],
                         ["bond_archive_schema", "chi_s_w", "chi_c_w", "S_bond", "C_bond"])


class TestSidecarWriter(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_sidecar_index_drops_the_channels_and_names_the_files(self):
        s, gi = _run_gate(self.tmp, extra={"longitudinal_bond_output_layout": "sidecar"})
        index = os.path.join(self.tmp, "longitudinal_bond.npz")
        s_npy = os.path.join(self.tmp, "longitudinal_bond_chi_s_w.npy")
        c_npy = os.path.join(self.tmp, "longitudinal_bond_chi_c_w.npy")
        self.assertTrue(os.path.exists(index) and os.path.exists(s_npy) and os.path.exists(c_npy))
        with np.load(index) as d:
            self.assertEqual(int(d["bond_archive_schema"]), 2)
            self.assertEqual(str(d["bond_archive_layout"]), "sidecar")
            self.assertNotIn("chi_s_w", d.files)
            self.assertNotIn("chi_c_w", d.files)
            self.assertEqual(str(d["chi_s_w_file"]), "longitudinal_bond_chi_s_w.npy")
            self.assertEqual(str(d["chi_c_w_file"]), "longitudinal_bond_chi_c_w.npy")
            # every other member the default archive carries is still present
            self.assertIn("S_bond", d.files)
            self.assertIn("delta_r", d.files)
        # the .npy sidecars hold exactly the channel data of the run
        np.testing.assert_array_equal(np.load(s_npy), gi["longitudinal_bond_chi_s_w"])
        np.testing.assert_array_equal(np.load(c_npy), gi["longitudinal_bond_chi_c_w"])

    def test_sidecar_index_carries_the_same_small_members_as_npz(self):
        tmp_npz = tempfile.mkdtemp()
        try:
            _run_gate(tmp_npz)
            _run_gate(self.tmp, extra={"longitudinal_bond_output_layout": "sidecar"})
            with np.load(os.path.join(tmp_npz, "longitudinal_bond.npz")) as dn, \
                    np.load(os.path.join(self.tmp, "longitudinal_bond.npz")) as ds:
                npz_only = set(dn.files) - set(ds.files)
                self.assertEqual(npz_only, {"chi_s_w", "chi_c_w"})
                sidecar_only = set(ds.files) - set(dn.files)
                self.assertEqual(sidecar_only, {"chi_s_w_file", "chi_c_w_file"})
                for k in set(dn.files) & set(ds.files):
                    if k == "bond_archive_layout":
                        continue
                    np.testing.assert_array_equal(dn[k], ds[k], err_msg=k)
        finally:
            shutil.rmtree(tmp_npz, ignore_errors=True)


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

    def test_non_finite_beta_refused(self):
        """``beta = nan`` slips through the relative comparison (every
        comparison with nan is false), and the run would then proceed on an
        archive whose temperature is unknown."""
        self._rewrite(beta=np.float64("nan"))
        with self.assertRaisesRegex(ValueError, "beta must be a finite positive number"):
            self._load()
        self._rewrite(beta=np.float64(0.0))
        with self.assertRaisesRegex(ValueError, "beta must be a finite positive number"):
            self._load()
        # the caller's own beta is held to the same rule
        self._rewrite(beta=np.float64(1.0 / self.s.T))
        with self.assertRaisesRegex(ValueError, "this run's beta"):
            self._load(beta_expected=float("inf"))

    def test_non_integral_topology_refused(self):
        """``delta_r`` / ``reverse`` used to be cast with ``astype(int64)``,
        which TRUNCATES: a float member carrying 0.5 would silently become a
        different bond topology."""
        with np.load(self.path) as d:
            dr = np.asarray(d["delta_r"], dtype=float)
            rev = np.asarray(d["reverse"], dtype=float)
            types = np.asarray(d["types"])
        bad = dr.copy()
        bad[1, 0] += 0.5
        self._rewrite(delta_r=bad)
        with self.assertRaisesRegex(ValueError, "delta_r"):
            self._load()
        self._rewrite(delta_r=dr)
        self._load()                      # the same values as floats are fine
        bad_rev = rev.copy()
        bad_rev[0] = 0.5
        self._rewrite(reverse=bad_rev)
        with self.assertRaisesRegex(ValueError, "reverse"):
            self._load()
        # ``types`` names the INTERACTION TYPES the vertices were built from,
        # not one entry per channel, so only emptiness is a contradiction
        self._rewrite(reverse=rev, types=types[:0])
        with self.assertRaisesRegex(ValueError, "types is empty"):
            self._load()
        self._rewrite(types=types)
        self._load()

    def test_delta_r_aliasing_refused(self):
        with np.load(self.path) as d:
            dr = d["delta_r"].copy()
            reverse = d["reverse"]
        (i,) = np.where((dr == np.array([1, 0, 0])).all(axis=1))
        i = int(i[0])
        j = int(reverse[i])

        # Two DIFFERENT channels that coincide modulo the 4x4x1 cell: (2,0,0)
        # and (-2,0,0) both reduce to (2,0,0) mod 4. The reversal involution
        # (delta_r[reverse] == -delta_r) still holds, so only the modulo
        # check can fire.
        dr_alias = dr.copy()
        dr_alias[i] = np.array([2, 0, 0])
        dr_alias[j] = np.array([-2, 0, 0])
        self._rewrite(delta_r=dr_alias)
        with self.assertRaisesRegex(ValueError, "modulo"):
            self._load()

        # A full-period shift of one row (its OWN residue unchanged, so it
        # does not collide with any other row) is NOT an alias and must be
        # accepted: the loader checks pairwise distinctness modulo
        # cell_shape, not that rows are already in canonical form.
        dr_shifted = dr.copy()
        dr_shifted[i] = dr[i] + np.array([4, 0, 0])
        dr_shifted[j] = -dr_shifted[i]
        self._rewrite(delta_r=dr_shifted)
        self._load()      # must not raise

    def test_malformed_archive_reports_the_contract(self):
        """A truncated / garbage file, or one missing a required member, must
        surface as a ValueError naming the archive contract, not a bare
        KeyError / BadZipFile / EOFError from deep in numpy."""
        from hwave.solver.eliashberg_bond_io import load_bond_archive

        def call(p):
            return load_bond_archive(
                p, norb=self.s.norb, nmat_expected=self.s.nmat,
                cell_shape_expected=tuple(int(x) for x in self.s.lattice.shape),
                beta_expected=1.0 / self.s.T)

        # a nonexistent path (the existing explicit message, naming the path)
        gone = os.path.join(self.tmp, "does_not_exist.npz")
        with self.assertRaisesRegex(ValueError, r"does_not_exist\.npz"):
            call(gone)
        # a few random bytes at a .npz path
        junk = os.path.join(self.tmp, "junk.npz")
        with open(junk, "wb") as f:
            f.write(os.urandom(41))
        with self.assertRaisesRegex(ValueError, r"junk\.npz.*not a valid bond archive"):
            call(junk)
        # a completely empty file
        empty = os.path.join(self.tmp, "empty.npz")
        open(empty, "wb").close()
        with self.assertRaisesRegex(ValueError, "not a valid bond archive"):
            call(empty)
        # a valid npz that is missing one required member
        with np.load(self.path) as d:
            members = {k: d[k] for k in d.files if k != "index_order"}
        sub = os.path.join(self.tmp, "sub.npz")
        np.savez(sub, **members)
        with self.assertRaisesRegex(ValueError, r"sub\.npz.*not a valid bond archive"):
            call(sub)
        # a missing chi member (read by header, not np.load) is caught too
        with np.load(self.path) as d:
            members = {k: d[k] for k in d.files if k != "chi_c_w"}
        no_chi = os.path.join(self.tmp, "no_chi.npz")
        np.savez(no_chi, **members)
        with self.assertRaisesRegex(ValueError, "not a valid bond archive"):
            call(no_chi)
        # the valid archive still loads
        self._load()


class TestSidecarLoader(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.s, self.gi = _run_gate(self.tmp)
        self.path = os.path.join(self.tmp, "longitudinal_bond.npz")
        self.s_npy = os.path.join(self.tmp, "longitudinal_bond_chi_s_w.npy")
        self.c_npy = os.path.join(self.tmp, "longitudinal_bond_chi_c_w.npy")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _load(self, **kw):
        from hwave.solver.eliashberg_bond_io import load_bond_archive
        args = dict(norb=self.s.norb, nmat_expected=self.s.nmat,
                    cell_shape_expected=tuple(int(x) for x in self.s.lattice.shape),
                    beta_expected=1.0 / self.s.T)
        args.update(kw)
        return load_bond_archive(self.path, **args)

    def test_sidecar_member_is_a_memmap_npz_member_is_plain(self):
        # the default npz layout returns a plain, fully-resident array
        a_npz = self._load()
        m_npz = a_npz.member("chi_s_w")
        self.assertNotIsInstance(m_npz, np.memmap)
        # the sidecar layout returns a lazily-read memmap
        _npz_to_sidecar(self.path)
        a = self._load()
        m = a.member("chi_s_w")
        self.assertIsInstance(m, np.memmap)
        self.assertEqual(m.shape, a.chi_shape)
        self.assertEqual(m.dtype, np.dtype(np.complex128))
        np.testing.assert_array_equal(np.asarray(m), self.gi["longitudinal_bond_chi_s_w"])

    def test_old_archive_without_layout_field_loads_as_npz(self):
        # an archive written before this field existed has no
        # bond_archive_layout member and must still load as the npz layout
        with np.load(self.path) as d:
            self.assertIn("bond_archive_layout", d.files)
            members = {k: d[k] for k in d.files if k != "bond_archive_layout"}
        np.savez(self.path, **members)
        a = self._load()
        m = a.member("chi_c_w")
        self.assertNotIsInstance(m, np.memmap)
        np.testing.assert_array_equal(m, self.gi["longitudinal_bond_chi_c_w"])

    def test_missing_sidecar_is_refused_with_the_contract(self):
        _npz_to_sidecar(self.path)
        os.remove(self.s_npy)
        with self.assertRaisesRegex(ValueError, "not a valid bond archive"):
            self._load()

    def test_garbage_sidecar_is_refused_with_the_contract(self):
        _npz_to_sidecar(self.path)
        with open(self.c_npy, "wb") as fh:
            fh.write(b"not a real npy file")
        with self.assertRaisesRegex(ValueError, "not a valid bond archive"):
            self._load()

    def test_wrong_shape_sidecar_is_refused(self):
        _npz_to_sidecar(self.path)
        np.save(self.s_npy, np.zeros((2, 2), dtype=np.complex128))
        with self.assertRaisesRegex(ValueError, "chi_s_w"):
            self._load()

    def test_wrong_dtype_sidecar_is_refused(self):
        _npz_to_sidecar(self.path)
        np.save(self.s_npy, np.zeros(self.s.nmat * self.s.lattice.nvol, dtype=np.float64)
                .reshape(self.s.nmat, self.s.lattice.nvol, 1, 1))
        with self.assertRaisesRegex(ValueError, "chi_s_w"):
            self._load()

    def test_truncated_sidecar_payload_is_refused(self):
        # a .npy whose header is intact but whose data body is cut short must
        # be caught during validation, not surface as a bare exception later
        _npz_to_sidecar(self.path)
        from numpy.lib import format as npfmt
        with open(self.s_npy, "rb") as fh:
            npfmt.read_magic(fh)
            npfmt.read_array_header_1_0(fh)
            data_offset = fh.tell()
        with open(self.s_npy, "r+b") as fh:
            fh.truncate(data_offset + 16)      # header + a few bytes only
        with self.assertRaisesRegex(ValueError, "not a valid bond archive"):
            self._load()


class TestPairingControls(unittest.TestCase):
    def test_defaults_and_validation(self):
        from hwave.solver.eliashberg_bond_io import PairingControls
        c = PairingControls.from_param({}, pairing_types=("singlet",))
        self.assertEqual((c.solver_mode, c.eigenvalue_method, c.num_eigenvalues, c.max_iter, c.alpha,
                          c.convergence_tol, c.matsubara_basis, c.ir_tol, c.ir_fit_tol, c.fft_workers),
                         ("iteration", "arnoldi", 10, 1000, 0.5, 1e-5, "uniform", 1e-8, 0.5, 1))
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

    def test_strict_parsing_of_the_remaining_keys(self):
        """A malformed value must be refused, not silently read as something
        else: ``"ture"`` is not false, ``fft_workers = 1.5`` is not 1, an
        unknown ``eigenvalue_method`` is not arnoldi, and ``alpha = true`` is
        not 1.0."""
        from hwave.solver.eliashberg_bond_io import PairingControls
        for bad in ({"ir_keep_static_chi": "ture"}, {"ir_keep_static_chi": 2.5},
                    {"ir_keep_static_chi": []},
                    {"fft_workers": 1.5}, {"fft_workers": 0}, {"fft_workers": True},
                    {"eigenvalue_method": "foo"},
                    {"alpha": True}, {"convergence_tol": True}):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError):
                    PairingControls.from_param(bad, pairing_types=("singlet",))
        for spelling, expected in (("TRUE", True), ("on", True), ("1", True), (True, True),
                                   ("No", False), ("off", False), ("0", False), (False, False)):
            c = PairingControls.from_param({"ir_keep_static_chi": spelling},
                                           pairing_types=("singlet",))
            self.assertIs(c.ir_keep_static_chi, expected, spelling)
        for method in ("arnoldi", "shift-invert-bicgstab", "shift-invert-gmres",
                       "shift-invert-lgmres"):
            c = PairingControls.from_param({"eigenvalue_method": method.upper()},
                                           pairing_types=("singlet",))
            self.assertEqual(c.eigenvalue_method, method)
        self.assertEqual(PairingControls.from_param({"fft_workers": 4},
                                                    pairing_types=("singlet",)).fft_workers, 4)
        # the documented "all cores" spelling of the dynamic solver stays valid,
        # down to the scipy.fft floor of -os.cpu_count(); one below is refused
        # here rather than by the first FFT after the FLEX solve
        import os
        floor = -(os.cpu_count() or 1)
        self.assertEqual(PairingControls.from_param({"fft_workers": -1},
                                                    pairing_types=("singlet",)).fft_workers, -1)
        self.assertEqual(PairingControls.from_param({"fft_workers": floor},
                                                    pairing_types=("singlet",)).fft_workers, floor)
        with self.assertRaisesRegex(ValueError, "fft_workers"):
            PairingControls.from_param({"fft_workers": floor - 1}, pairing_types=("singlet",))

    def test_parity_leakage_tol_default_depends_on_the_basis(self):
        """``parity_leakage_tol`` is optional; ``None`` resolves to 1e-8 on the
        uniform grid and 2e-2 with ``matsubara_basis = "ir"`` (spec 4.4: the IR
        representation of a uniform-FFT archive carries a parity asymmetry
        decaying as Nmat^-2, measured 7.6e-3 / 1.3e-3 at Nmat 64 / 128)."""
        from hwave.solver.eliashberg_bond_io import PairingControls
        c = PairingControls.from_param({}, pairing_types=("singlet",))
        self.assertIsNone(c.parity_leakage_tol)
        self.assertEqual(c.resolved_parity_leakage_tol, 1.0e-8)
        cir = PairingControls.from_param({"matsubara_basis": "ir"},
                                         pairing_types=("singlet",))
        self.assertIsNone(cir.parity_leakage_tol)
        self.assertEqual(cir.resolved_parity_leakage_tol, 2.0e-2)
        # an explicit value wins on either basis, zero included
        for basis in ("uniform", "ir"):
            c2 = PairingControls.from_param(
                {"matsubara_basis": basis, "parity_leakage_tol": 3.5e-3},
                pairing_types=("singlet",))
            self.assertEqual(c2.resolved_parity_leakage_tol, 3.5e-3)
        c0 = PairingControls.from_param({"parity_leakage_tol": 0},
                                        pairing_types=("singlet",))
        self.assertEqual(c0.resolved_parity_leakage_tol, 0.0)
        for bad in ({"parity_leakage_tol": -1e-3}, {"parity_leakage_tol": float("nan")}):
            with self.assertRaises(ValueError):
                PairingControls.from_param(bad, pairing_types=("singlet",))


class TestParityLeakageTolerance(unittest.TestCase):
    """The configurable refusal threshold of the parity probe (spec 4.4).

    The operator is the deliberately non-commuting one of the shared
    eigen-driver's unit tests (a dense random complex matrix): its
    cross-sector leakage is O(1), so a tolerance above it must let the solve
    through (with the warning of the ``[0.1 * tol, tol)`` band) and a
    tolerance below it must refuse, naming the key and the tolerance.
    """

    def setUp(self):
        from scipy.sparse.linalg import LinearOperator
        from hwave.solver import eliashberg_dynamic as ed
        self.ed = ed
        self.gap_shape = (1, 1, 2, 2, 1, 4)
        n = int(np.prod(self.gap_shape))
        rng = np.random.default_rng(11)
        M = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
        self.matvec = lambda x: M @ x
        kx = np.array([0.0, np.pi])
        ky = np.array([0.0, np.pi])
        kz = np.array([0.0])
        self.phi0, self.seed = ed.build_seed({}, "singlet", 1, kx, ky, kz,
                                             self.gap_shape, False, None, 4)
        self.leak = ed._parity_leakage(
            LinearOperator((n, n), matvec=self.matvec, dtype=complex),
            self.gap_shape, "singlet")

    def _run(self, **kw):
        return self.ed.run_leading_eigenproblem(
            self.matvec, self.gap_shape, {"solver_mode": "iteration", "max_iter": 5},
            "singlet", phi0=self.phi0, seed_vec=self.seed, use_ir=False, axF=None,
            nmat=4, parity_leakage_policy="refuse", **kw)

    def test_tolerance_above_the_leakage_completes_and_warns(self):
        self.assertGreater(self.leak, 1.0e-8)
        with self.assertLogs("qlms.eliashberg_dynamic", level="WARNING") as cm:
            out = self._run(parity_leakage_tol=1.5 * self.leak)
        self.assertEqual(len(out), 9)
        self.assertIsInstance(out[5], float)
        self.assertAlmostEqual(out[5], self.leak)
        # the 7th element is the returned gap's sector composition, the 8th
        # the sector the iteration projected onto, the 9th which eigenvalue
        # selection criterion produced the leading value
        self.assertAlmostEqual(sum(out[6].values()), 1.0, places=9)
        self.assertIn(out[7], ("channel", "combined_parity", "none"))
        self.assertEqual(out[8], "iteration")
        self.assertTrue(any("parity_leakage_tol" in m for m in cm.output), cm.output)

    def test_tolerance_below_the_leakage_refuses_naming_the_key(self):
        with self.assertRaisesRegex(ValueError, "parity_leakage_tol") as cm:
            self._run(parity_leakage_tol=0.5 * self.leak)
        self.assertIn("cross-sector leakage", str(cm.exception))

    def test_silent_below_a_tenth_of_the_tolerance(self):
        """The lower edge of the warning band: a leakage below ``0.1 * tol`` is
        accepted SILENTLY, so the warning of the sibling case above is a signal
        and not something every accepted run emits."""
        with self.assertNoLogs("qlms.eliashberg_dynamic", level="WARNING"):
            out = self._run(parity_leakage_tol=100.0 * self.leak)
        self.assertAlmostEqual(out[5], self.leak)

    def test_default_tolerance_is_1e_minus_8(self):
        import inspect
        # the default itself, not only its consequence on this leaky operator
        self.assertEqual(
            inspect.signature(self.ed.run_leading_eigenproblem)
            .parameters["parity_leakage_tol"].default, 1.0e-8)
        with self.assertRaisesRegex(ValueError, "cross-sector leakage"):
            self._run()

    def test_invalid_tolerance_refused(self):
        for bad in (-1.0, float("nan"), None):
            with self.assertRaisesRegex(ValueError, "parity_leakage_tol"):
                self._run(parity_leakage_tol=bad)


# ---------------------------------------------------------------------------
# Post-processing entry through hwave_sc (spec 6, 10.2.1 / 10.2.2)
# ---------------------------------------------------------------------------

_IN1 = "tests/rpa/input"
_IN2 = "tests/rpa/input_2orb"


def _sc_input(flex_dir, out_dir, T, nmat, cell, **eli):
    inp = {"mode": {"param": {"T": T, "CellShape": list(cell), "SubShape": [1, 1, 1],
                              "Nmat": nmat, "filling": 0.5}},
           "file": {"input": {"path_to_flex_output": flex_dir,
                              "interaction": {"path_to_input": _IN1, "Geometry": "geom.dat",
                                              "Transfer": "transfer.dat",
                                              "CoulombIntra": "coulombintra.dat",
                                              "CoulombInter": "coulombinter.dat"}},
                    "output": {"path_to_output": out_dir}},
           "eliashberg": {"chi0q_mode": "flex", "frequency": "dynamic", "bond_channels": True,
                          "pairing_type": "singlet", "solver_mode": "eigenvalue",
                          "num_eigenvalues": 3}}
    inp["eliashberg"].update(eli)
    return inp


def _run_gate_1orb(tmp, T=0.5, nmat=64, cell=(4, 4, 1), output_full=True, pairing=None, eli=None):
    """A real single-band (U = 4, V = 1) bond-gate FLEX run through ``qlms.run``,
    writing green.npz and the bond archive into ``tmp``."""
    import hwave.qlms as qlms
    par = {"T": T, "filling": 0.5, "CellShape": list(cell), "SubShape": [1, 1, 1], "Nmat": nmat,
           "IterationMax": 200, "Mix": 0.3, "EPS": 8, "flex_hartree_fock": True,
           "longitudinal_bond_channels": True, "longitudinal_bond_output_full": output_full,
           "mixing_scheme": "linear"}
    if pairing is not None:
        par["longitudinal_bond_pairing"] = pairing
    inp = {"mode": {"mode": "FLEX", "calc_scheme": "general", "param": par},
           "file": {"input": {"path_to_input": "",
                              "interaction": {"path_to_input": _IN1, "Geometry": "geom.dat",
                                              "Transfer": "transfer.dat",
                                              "CoulombIntra": "coulombintra.dat",
                                              "CoulombInter": "coulombinter.dat"}},
                    "output": {"path_to_output": tmp, "chiq_s": "chiq_s", "chiq_c": "chiq_c",
                               "green": "green", "longitudinal_bond": "longitudinal_bond"}}}
    if eli is not None:
        inp["eliashberg"] = eli
    qlms.run(input_dict=inp)
    return inp


class TestPostProcessingGuards(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _refused(self, needle, **eli):
        import hwave.sc as sc
        with self.assertRaisesRegex(ValueError, needle):
            sc.calc_eliashberg(_sc_input(self.tmp, self.tmp, 0.5, 8, (4, 4, 1), **eli))

    def test_refusals(self):                                # 10.2.1 (post-processing half)
        self._refused("bond_green", bond_green="x.npz")
        self._refused("bond_max_shells", bond_max_shells=1)
        self._refused("zero_chi", zero_chi_s=True)
        self._refused("zero_chi", zero_chi_c=True)
        import hwave.sc as sc
        inp = _sc_input(self.tmp, self.tmp, 0.5, 9, (4, 4, 1))
        with self.assertRaisesRegex(ValueError, "even Nmat"):
            sc.calc_eliashberg(inp)
        inp = _sc_input(self.tmp, self.tmp, 0.5, 8, (4, 4, 1))
        inp["eliashberg"]["chi0q_mode"] = "rpa"
        with self.assertRaisesRegex(ValueError, "chi0q_mode='flex'"):
            sc.calc_eliashberg(inp)
        # the archive is missing in tmp: the message names the FLEX prerequisite
        self._refused("longitudinal_bond_output_full = true")

    def test_direct_solve_dynamic_call_reaches_the_bond_path(self):
        """A DIRECT ``solve_dynamic`` call (bypassing ``calc_eliashberg``'s
        dispatch) must not silently run the scalar on-site vertex with
        ``bond_channels = true`` ignored: it delegates to the bond entry, which
        here gets as far as the missing archive."""
        from hwave.solver import eliashberg_dynamic
        inp = _sc_input(self.tmp, self.tmp, 0.5, 8, (4, 4, 1))
        with self.assertRaisesRegex(ValueError, "longitudinal_bond_output_full = true"):
            eliashberg_dynamic.solve_dynamic(inp)
        # and the strict reader still applies on that route
        inp["eliashberg"]["bond_channels"] = "ture"
        with self.assertRaisesRegex(ValueError, "bond_channels"):
            eliashberg_dynamic.solve_dynamic(inp)

    def test_static_bond_guards_unchanged(self):
        import hwave.sc as sc
        inp = _sc_input(self.tmp, self.tmp, 0.5, 8, (4, 4, 1))
        inp["eliashberg"]["frequency"] = "static"
        # the STATIC bond path still refuses chi0q_mode = "flex" (unchanged)
        with self.assertRaisesRegex(ValueError, "(?i)flex"):
            sc.calc_eliashberg(inp)


class TestFrequencyBatch(unittest.TestCase):
    """``_frequency_batch`` is a pure function of the two caps (spec 8): one
    dressing batch of ``7 * nb * nvol * ND**2 * 16`` bytes lives on the host AND
    on the device, so it must fit a quarter of BOTH caps."""

    def _nb(self, host_gib, device_gib=None, nmat=64, nvol=16, ND=45):
        from hwave.solver.eliashberg_bond_io import _frequency_batch
        gib = 1024 ** 3
        return _frequency_batch(nmat, nvol, ND, host_gib * gib,
                                device_cap=(None if device_gib is None else device_gib * gib))

    def test_host_only_and_device_limited(self):
        per_l = 7 * 16 * 45 * 45 * 16                      # 3.629 MB per frequency
        gib = 1024 ** 3
        # one address space (numpy): the host cap alone decides
        self.assertEqual(self._nb(1.0), min(64, int(0.25 * gib // per_l)))  # 73 -> nmat
        self.assertEqual(self._nb(1.0, None), self._nb(1.0))
        self.assertEqual(self._nb(0.1), int(0.25 * 0.1 * gib // per_l))     # 7
        # a device backend with a roomy host: the DEVICE cap is what binds, and
        # passing it must shrink nb (this is the whole point of the argument)
        big_host = self._nb(64.0)
        self.assertEqual(big_host, 64)                                  # nmat-capped
        self.assertLess(self._nb(64.0, 0.05), big_host)
        self.assertEqual(self._nb(64.0, 0.05), int(0.25 * 0.05 * gib // per_l))
        # the smaller of the two always wins, either way round
        self.assertEqual(self._nb(0.05, 64.0), self._nb(0.05))
        # never zero, never above nmat
        self.assertEqual(self._nb(1e-9, 1e-9), 1)
        self.assertEqual(self._nb(1e6, 1e6, nmat=8), 8)


class TestDeviceProbeFallback(unittest.TestCase):
    """Both entries budget device bytes before they admit a channel, and
    ``backend.device_available_bytes()`` is allowed to have no answer --
    it returns ``None`` when the query FAILS as well as when there is no
    device. Multiplying that would end a post-processing run with a
    ``TypeError`` about ``NoneType`` and turn an in-process channel's error
    into the same, so the entries fall back to the host cap and warn."""

    def test_entry_admission_falls_back_to_the_host_cap(self):
        from hwave.solver import backend as bk
        from hwave.solver.eliashberg_bond_io import _device_cap
        gib = 1024.0 ** 3
        with mock.patch.object(bk, "device_available_bytes", return_value=None):
            with self.assertLogs("qlms.solver.eliashberg_bond", "WARNING") as log:
                cap = _device_cap(3.0 * gib, "bond pairing admission")
        self.assertEqual(cap, 3.0 * gib)
        self.assertIn("device memory probe", "\n".join(log.output))
        self.assertIn("bond pairing admission", "\n".join(log.output))
        # the warning belongs to the entry module's logger, not the kernel's,
        # so a user filtering on qlms.solver sees it
        self.assertTrue(all(r.startswith("WARNING:qlms.solver.eliashberg_bond")
                            for r in log.output), log.output)
        with mock.patch.object(bk, "device_available_bytes", return_value=2000):
            self.assertEqual(_device_cap(3.0 * gib, "bond pairing admission"), 1800.0)


class TestPostProcessingRuns(unittest.TestCase):
    def test_sidecar_layout_matches_npz_layout(self):      # 10.2 (issue #205)
        """The bond archive written in the sidecar layout, memory-mapped by
        post-processing, gives the IDENTICAL leading eigenvalue and gap as the
        default single-file layout on the SAME data. One FLEX run; its archive
        is solved once as npz and once after rewriting it into the sidecar
        layout (the two big members moved to .npy files next to the index),
        on BOTH the uniform and the IR pairing paths -- the IR residual pass
        re-reads each sidecar member, so it exercises repeated memmap
        open/close.

        The leading eigenvalue matches to 1e-13 (measured |diff| 5.3e-15). The
        output gap is already gauge-fixed (_fix_gauge), so it is compared
        directly at rtol = 0, atol = 1e-11; the measured gap max|diff| is
        3.2e-13 (uniform) / 1.9e-14 (IR) at Nmat 64 -- the eigensolver's
        reduction-order noise, since a memmap can take a different BLAS path.
        """
        import hwave.sc as sc
        bases = [dict(tag="uniform", kw={})]
        try:
            import sparse_ir  # noqa: F401
            bases.append(dict(tag="ir", kw=dict(matsubara_basis="ir", ir_tol=1e-8,
                                                ir_keep_static_chi=True)))
        except ImportError:
            pass
        flex_dir = tempfile.mkdtemp()
        try:
            _run_gate_1orb(flex_dir, nmat=64)
            index = os.path.join(flex_dir, "longitudinal_bond.npz")
            # solve every basis on the intact npz archive FIRST (the sidecar
            # rewrite below is one-way)
            ref = {}
            for b in bases:
                out = tempfile.mkdtemp()
                try:
                    ref[b["tag"]] = self._solve_lambda_gap(sc, flex_dir, out, b["kw"])
                finally:
                    shutil.rmtree(out, ignore_errors=True)
            with np.load(index) as d:
                self.assertEqual(str(d["bond_archive_layout"]), "npz")
                self.assertIn("chi_s_w", d.files)
            # rewrite the archive in place into the sidecar layout
            _npz_to_sidecar(index)
            self.assertTrue(os.path.exists(
                os.path.join(flex_dir, "longitudinal_bond_chi_s_w.npy")))
            with np.load(index) as d:
                self.assertEqual(str(d["bond_archive_layout"]), "sidecar")
                self.assertNotIn("chi_s_w", d.files)
            for b in bases:
                out = tempfile.mkdtemp()
                try:
                    lam_side, gap_side = self._solve_lambda_gap(sc, flex_dir, out, b["kw"])
                finally:
                    shutil.rmtree(out, ignore_errors=True)
                lam_npz, gap_npz = ref[b["tag"]]
                self.assertTrue(np.isfinite(lam_npz))
                self.assertAlmostEqual(lam_side, lam_npz,
                                       delta=1e-13 * max(1.0, abs(lam_npz)), msg=b["tag"])
                maxdiff = float(np.max(np.abs(gap_side - gap_npz)))
                np.testing.assert_allclose(gap_side, gap_npz, rtol=0, atol=1e-11,
                                           err_msg="{}: gap max|diff|={:.3e}".format(b["tag"], maxdiff))
        finally:
            shutil.rmtree(flex_dir, ignore_errors=True)

    @staticmethod
    def _solve_lambda_gap(sc, flex_dir, out, kw):
        lam = sc.calc_eliashberg(_sc_input(flex_dir, out, 0.5, 64, (4, 4, 1), **kw))
        with np.load(os.path.join(out, "gap_dynamic.npz")) as d:
            gap = np.array(d["gap"])
        return lam, gap

    def test_uniform_lambda_and_outputs(self):              # 10.2.2 (uniform half)
        """A real bond-gate FLEX run, then ``hwave_sc`` on its archive for both
        pairing channels: the leading eigenvalue is a finite real number, the
        documented file set is written with the bond provenance, and the answer
        is stable in the FLEX run's own Nmat (64 vs 128, the discretization the
        post-processing inherits rather than one it chooses).

        The IR half of spec 10.2.2 -- the same lambda from
        ``matsubara_basis = "ir"`` -- is
        :meth:`test_ir_lambda_arm_matches_uniform` below, on this same working
        point.
        """
        import hwave.sc as sc
        lam = {}
        for nmat in (64, 128):
            flex_dir = tempfile.mkdtemp()
            try:
                _run_gate_1orb(flex_dir, nmat=nmat)
                for eta in ("singlet", "triplet"):
                    out = tempfile.mkdtemp()
                    try:
                        lam[(eta, nmat)] = sc.calc_eliashberg(_sc_input(
                            flex_dir, out, 0.5, nmat, (4, 4, 1), pairing_type=eta))
                        self.assertTrue(np.isfinite(lam[(eta, nmat)]))
                        for name in ("eigenvalue.dat", "gap.dat", "gap_dynamic.npz"):
                            self.assertTrue(os.path.exists(os.path.join(out, name)), name)
                        head = open(os.path.join(out, "eigenvalue.dat")).read()
                        self.assertIn("# bond_channels=true", head)
                        self.assertIn("# residency=", head)
                        with np.load(os.path.join(out, "gap_dynamic.npz")) as d:
                            self.assertTrue(bool(d["bond_channels"]))
                            self.assertEqual(str(d["pairing_type"]), eta)
                            self.assertEqual(d["gap"].shape, (1, 1, 4, 4, 1, nmat))
                            # the bond provenance: B = on-site + the four
                            # nearest neighbours of the V = 1 fixture
                            self.assertEqual(d["bond_delta_r"].shape, (5, 3))
                            self.assertEqual(tuple(d["bond_delta_r"][0]), (0, 0, 0))
                            self.assertEqual(d["bond_reverse"].shape, (5,))
                            self.assertTrue(str(d["bond_archive"]).endswith(
                                "longitudinal_bond.npz"))
                            self.assertIn(str(d["bond_residency"]),
                                          ("device", "host", "stream"))
                            self.assertEqual(d["gap_bond_projection"].shape,
                                             (5, 1, 1, nmat))
                    finally:
                        shutil.rmtree(out, ignore_errors=True)
            finally:
                shutil.rmtree(flex_dir, ignore_errors=True)
        for eta in ("singlet", "triplet"):
            # The reported eigenvalue is now the channel-projected leading value
            # (issue #202), which equals the dense even-frequency channel
            # eigenvalue exactly; on this 4x4 fixture its Nmat 64 -> 128
            # convergence is ~2.6% (triplet), so the coarse cross-Nmat check
            # uses a 3% band. Both Nmat take the same projected path.
            self.assertAlmostEqual(lam[(eta, 128)], lam[(eta, 64)],
                                   delta=3e-2 * abs(lam[(eta, 64)]), msg=eta)

    @heavy
    def test_ir_lambda_arm_matches_uniform(self):           # 10.2.2 (IR half)
        """The IR arm of 10.2.2 on the SAME real bond-gate archives: the IR
        basis must reproduce the uniform-grid eigenvalue, and the discrepancy
        must be a discretization one.

        Both bases are solved on the ONE archive per Nmat, so the only
        difference is the Matsubara representation. At Nmat 128 the two agree
        to the spec's 2e-2 relative band, and the disagreement at Nmat 128 is
        smaller than at Nmat 64 for both channels -- measured

            singlet  -0.1464292 (uniform) vs -0.1547410 (IR)  5.7 % at Nmat  64
                     -0.1455735           vs -0.1480647       1.7 % at Nmat 128
            triplet  -0.1987357           vs -0.2047572       3.0 % at Nmat  64
                     -0.1984555           vs -0.2002307       0.9 % at Nmat 128

        so the residual gap is the IR representability of a uniform-FFT
        archive, not a defect of the kernel algebra (which is exact to 1e-12
        on an IR-representable fixture:
        ``test_eliashberg_bond_kernel.TestKernelIR.test_ir_parity_commutation``).

        That same representability is what ``[eliashberg] parity_leakage_tol``
        exists for (spec 4.4): the measured cross-sector leakage is 7.56e-03 /
        9.70e-03 at Nmat 64 and 1.29e-03 / 1.68e-03 at Nmat 128, i.e. roughly
        ``Nmat^-2`` and below the IR default 2e-2, where the uniform grid sits
        at 1e-15. The Task 7b convergence guard is kept here as the third
        assertion. The componentwise IR fit residual with the constant
        retained is 0.161 / 0.157, below the ``ir_fit_tol`` default 0.5 (and
        inside its warning band); both numbers are recorded in the npz, which
        the test also asserts.
        """
        try:
            import sparse_ir  # noqa: F401
        except ImportError:
            raise unittest.SkipTest("sparse-ir not installed")
        import hwave.sc as sc

        lam, leak, resid = {}, {}, {}
        for nmat in (64, 128):
            flex_dir = tempfile.mkdtemp()
            try:
                _run_gate_1orb(flex_dir, nmat=nmat)
                for eta in ("singlet", "triplet"):
                    for basis in ("uniform", "ir"):
                        kw = dict(pairing_type=eta)
                        if basis == "ir":
                            # the defaults of spec rev 8 (ir_fit_tol 0.5,
                            # parity_leakage_tol 2e-2) must accept these real
                            # archives: nothing is loosened here
                            kw.update(matsubara_basis="ir", ir_tol=1e-8,
                                      ir_keep_static_chi=True)
                        out = tempfile.mkdtemp()
                        try:
                            lam[(basis, eta, nmat)] = sc.calc_eliashberg(_sc_input(
                                flex_dir, out, 0.5, nmat, (4, 4, 1), **kw))
                            self.assertTrue(np.isfinite(lam[(basis, eta, nmat)]))
                            with np.load(os.path.join(out, "gap_dynamic.npz")) as d:
                                # the measured parity leakage is recorded on
                                # BOTH bases; the IR fit residual only on IR
                                self.assertIn("bond_parity_leakage", d.files)
                                leak[(basis, eta, nmat)] = float(d["bond_parity_leakage"])
                                if basis == "ir":
                                    self.assertIn("bond_ir_fit_residual_rel", d.files)
                                    resid[(eta, nmat)] = float(
                                        np.asarray(d["bond_ir_fit_residual_rel"]).max())
                                else:
                                    self.assertNotIn("bond_ir_fit_residual_rel", d.files)
                            head = open(os.path.join(out, "eigenvalue.dat")).read()
                            self.assertIn("# parity_leakage=", head)
                        finally:
                            shutil.rmtree(out, ignore_errors=True)
            finally:
                shutil.rmtree(flex_dir, ignore_errors=True)

        for eta in ("singlet", "triplet"):
            rel = {n: abs(lam[("ir", eta, n)] - lam[("uniform", eta, n)])
                   / abs(lam[("uniform", eta, n)]) for n in (64, 128)}
            self.assertLessEqual(
                rel[128], 2e-2,
                "{}: lambda_ir {:.9f} vs lambda_uniform {:.9f} at Nmat 128 "
                "({:.3e} relative)".format(eta, lam[("ir", eta, 128)],
                                           lam[("uniform", eta, 128)], rel[128]))
            self.assertLess(rel[128], rel[64],
                            "{}: no Nmat convergence of lambda_ir, {:.3e} -> {:.3e}"
                            .format(eta, rel[64], rel[128]))
            # the Task 7b guard: the IR parity leakage of the archive converges
            # (and the uniform grid is at machine precision either way)
            self.assertLess(leak[("ir", eta, 64)], 5e-2, eta)
            self.assertLess(leak[("ir", eta, 128)], 0.6 * leak[("ir", eta, 64)],
                            "{}: no Nmat convergence of the parity leakage "
                            "{:.3e} -> {:.3e}".format(eta, leak[("ir", eta, 64)],
                                                      leak[("ir", eta, 128)]))
            self.assertLess(leak[("uniform", eta, 64)], 1e-10, eta)
            # and the IR fit residual stays under the shipped ir_fit_tol
            for n in (64, 128):
                self.assertLess(resid[(eta, n)], 0.5, (eta, n))


class TestQlmsForwarding(unittest.TestCase):
    """``qlms.run`` hands the ``[eliashberg]`` table of the SAME input file to
    the FLEX solver (spec 7 / 9).

    Fast end-to-end check of the forwarding seam only: one ``IterationMax = 1``
    two-orbital solve, asserting that the constructed solver carries the parsed
    :class:`PairingControls` (the new ``parity_leakage_tol`` included) and that
    the requested channel produced an eigenvalue artifact.
    """

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_run_forwards_the_eliashberg_table(self):
        import hwave.qlms as qlms
        par = {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1],
               "Nmat": 8, "IterationMax": 1, "Mix": 0.5, "EPS": 1e-12,
               "flex_hartree_fock": True, "mixing_scheme": "linear",
               "longitudinal_bond_channels": True,
               "longitudinal_bond_pairing": "singlet"}
        inp = {"mode": {"mode": "FLEX", "calc_scheme": "general", "param": par},
               "file": {"input": {"path_to_input": "",
                                  "interaction": {"path_to_input": _IN2,
                                                  "Geometry": "geom.dat",
                                                  "Transfer": "transfer.dat",
                                                  "CoulombInter": "coulombinter.dat"}},
                        "output": {"path_to_output": self.tmp, "green": "green"}},
               "eliashberg": {"solver_mode": "eigenvalue", "num_eigenvalues": 3,
                              "parity_leakage_tol": 5.0e-3}}
        captured = []
        real_init = qlms.sol_flex.FLEX.__init__

        def capture_init(solver, *a, **kw):
            real_init(solver, *a, **kw)
            captured.append(solver)

        with mock.patch.object(qlms.sol_flex.FLEX, "__init__", capture_init):
            qlms.run(input_dict=inp)
        self.assertEqual(len(captured), 1)
        ctl = captured[0]._pairing_controls
        self.assertIsNotNone(ctl)
        self.assertEqual(ctl.pairing_types, ("singlet",))
        self.assertEqual((ctl.solver_mode, ctl.num_eigenvalues), ("eigenvalue", 3))
        # the forwarded key reaches the resolved tolerance, and the default
        # ir_fit_tol of spec 4.5 is the one the table did not override
        self.assertEqual(ctl.parity_leakage_tol, 5.0e-3)
        self.assertEqual(ctl.resolved_parity_leakage_tol, 5.0e-3)
        self.assertEqual(ctl.ir_fit_tol, 0.5)
        # the in-process half of "the measured leakage is recorded in the
        # outputs" (spec 4.4): the eigenvalue header line and the npz member
        ev = os.path.join(self.tmp, "eigenvalue_bond_singlet.dat")
        self.assertTrue(os.path.exists(ev))
        with open(ev) as f:
            txt = f.read()
        self.assertIn("# parity_leakage=", txt)
        with np.load(os.path.join(self.tmp, "eliashberg_bond_singlet.npz")) as d:
            self.assertIn("bond_parity_leakage", d.files)
            leakage = float(d["bond_parity_leakage"])
        self.assertTrue(np.isfinite(leakage))
        self.assertGreaterEqual(leakage, 0.0)
        # the header line and the npz member are the SAME measurement
        self.assertAlmostEqual(
            float(txt.split("# parity_leakage=")[1].split()[0]), leakage,
            delta=1e-6 * max(leakage, 1e-12))


class TestInProcessHookBoundary(unittest.TestCase):
    """``run_inprocess_pairing`` promises never to raise: a pairing failure
    must cost the pairing results only, never the FLEX ones. The promise has
    to cover the handler's OWN setup too, not just the worker call."""

    def _solver(self):
        import types as _t
        return _t.SimpleNamespace(
            _pairing_controls=_t.SimpleNamespace(pairing_types=("singlet",)))

    def test_worker_failure_is_recorded_not_raised(self):
        from hwave.solver import eliashberg_bond_io as io
        gi = {}
        with mock.patch.object(io, "_run_inprocess_pairing",
                               side_effect=RuntimeError("worker boom")):
            with self.assertLogs("qlms.solver.eliashberg_bond", "ERROR"):
                io.run_inprocess_pairing(self._solver(), None, None, None, 1.0, gi)
        self.assertIn("worker boom", gi["pairing_singlet_error"])

    def test_a_failure_in_the_handlers_own_setup_is_recorded_too(self):
        """The backend import and the out-of-memory type tuple are built
        INSIDE the guard; built outside, a failure there would propagate and
        take the FLEX results with it."""
        from hwave.solver import backend as bk, eliashberg_bond_io as io
        gi = {}
        with mock.patch.object(bk, "_oom_error_types", side_effect=RuntimeError("probe boom")):
            with self.assertLogs("qlms.solver.eliashberg_bond", "ERROR"):
                io.run_inprocess_pairing(self._solver(), None, None, None, 1.0, gi)
        self.assertIn("probe boom", gi["pairing_singlet_error"])


class TestTemporaryNames(unittest.TestCase):
    """A fixed ``x.tmp.npz`` collides between two runs writing the same
    directory (a parameter sweep, a restarted job): one run's rename can
    publish the other's half-written file. The temporary carries the pid and
    a random token, and starts with a dot so it is not mistaken for output."""

    def test_tmp_name_is_unique_hidden_and_local(self):
        from hwave.solver.eliashberg_bond_io import _tmp_name
        path = os.path.join("/some", "where", "gap_bond_singlet.dat")
        a = _tmp_name(path)
        b = _tmp_name(path)
        self.assertEqual(os.path.dirname(a), os.path.dirname(path))
        self.assertTrue(os.path.basename(a).startswith("."), a)
        self.assertTrue(a.endswith(".dat"), a)
        self.assertIn("gap_bond_singlet", os.path.basename(a))
        self.assertNotEqual(a, b)
        npz = _tmp_name(os.path.join("/some", "where", "eliashberg_bond_singlet.npz"))
        self.assertTrue(npz.endswith(".npz"), npz)


class TestInProcess(unittest.TestCase):
    """The in-process entry at the end of the FLEX solve (spec 7, 10.2)."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _outputs(self, tmp):
        return {"path_to_output": tmp, "chiq_s": "chiq_s", "chiq_c": "chiq_c", "green": "green",
                "sigma": "sigma", "longitudinal_bond": "longitudinal_bond"}

    def test_outputs_written(self):                                # 10.2.6
        s, r = _flex({"longitudinal_bond_pairing": "both", "IterationMax": 3})
        gi = r.get_param("green")
        s.solve(gi, self.tmp)
        for eta in ("singlet", "triplet"):
            self.assertIn("pairing_{}_eigenvalue".format(eta), gi)
            self.assertNotIn("pairing_{}_error".format(eta), gi)
        s.validate_output_paths(self._outputs(self.tmp), self.tmp)
        s.save_results(self._outputs(self.tmp), gi)
        for eta in ("singlet", "triplet"):
            npz = os.path.join(self.tmp, "eliashberg_bond_{}.npz".format(eta))
            with np.load(npz) as d:
                for k in ("gap", "eigenvalue", "scf_converged", "state", "gap_bond_projection",
                          "bond_delta_r", "gap_sector_weights", "gap_sector_labels",
                          "sector_selection"):
                    self.assertIn(k, d.files, k)
                self.assertEqual(str(d["pairing_type"]), eta)
                # the in-process entry records the same two keys, under the
                # same names, as the file entry; the projected iteration
                # returns a gap in the channel's even-frequency sector
                labels = [str(x) for x in d["gap_sector_labels"]]
                self.assertEqual(labels, ["even_k_even_w", "odd_k_even_w",
                                          "even_k_odd_w", "odd_k_odd_w"])
                w = dict(zip(labels, np.asarray(d["gap_sector_weights"], dtype=float)))
                self.assertAlmostEqual(sum(w.values()), 1.0, places=9)
                own = "even_k_even_w" if eta == "singlet" else "odd_k_even_w"
                self.assertGreater(w[own], 0.999, (eta, w))
                # this two-orbital kernel conserves only the combined parity,
                # so the iteration is NOT restricted to the channel sector --
                # yet its leading eigenvector lands there anyway, which is
                # exactly what the recorded weights are for
                self.assertEqual(str(d["sector_selection"]),
                                 "combined_parity", eta)
            with open(os.path.join(self.tmp, "eigenvalue_bond_{}.dat".format(eta))) as f:
                txt = f.read()
            self.assertIn("scf_converged=", txt)
            self.assertIn("state=", txt)
            self.assertIn("# gap_sector_weights even_k_even_w=", txt)
            self.assertIn("# sector_selection=combined_parity", txt)
            self.assertTrue(os.path.exists(os.path.join(self.tmp, "gap_bond_{}.dat".format(eta))))
            for stem in ("eliashberg_bond_{}.tmp.npz", "gap_bond_{}.tmp.dat",
                         "eigenvalue_bond_{}.tmp.dat"):
                self.assertFalse(os.path.exists(os.path.join(self.tmp, stem.format(eta))))
        from hwave.solver import flex_bond
        with self.assertRaisesRegex(ValueError, "resolve to the same file"):
            flex_bond.resolve_output_paths(
                dict(self._outputs(self.tmp), eliashberg_bond_singlet="chiq_s"),
                self.tmp, ("chiq_s", "eliashberg_bond_singlet"))

    def test_cpu_run_never_touches_the_device_pool(self):
        """``free_device_pool`` swallows everything it raises, so on a CPU run
        the cupy import is merely a pointless attempt per channel -- but on a
        machine where cupy IS installed while the run is on numpy it would
        free another consumer's blocks. It must not be called at all."""
        from hwave.solver import backend as bk
        s, r = _flex({"longitudinal_bond_pairing": "both", "IterationMax": 3})
        gi = r.get_param("green")
        with mock.patch.object(bk, "free_device_pool",
                               side_effect=AssertionError("must not be called")) as spy, \
                mock.patch.object(bk, "_import_cupy",
                                  side_effect=AssertionError("must not be called")):
            s.solve(gi, self.tmp)
        for eta in ("singlet", "triplet"):
            self.assertIn("pairing_{}_eigenvalue".format(eta), gi)
            self.assertNotIn("pairing_{}_error".format(eta), gi)
        self.assertEqual(spy.call_count, 0)

    def test_configured_subdirectory_name_is_published(self):
        """A configured output name may carry a subdirectory. The temporary is
        written next to its target (so the rename stays on one filesystem), the
        directory is created, and nothing temporary survives."""
        s, r = _flex({"longitudinal_bond_pairing": "singlet", "IterationMax": 3})
        gi = r.get_param("green")
        s.solve(gi, self.tmp)
        self.assertNotIn("pairing_singlet_error", gi)
        outputs = dict(self._outputs(self.tmp),
                       gap_bond_singlet=os.path.join("sub", "gap_bond_singlet.dat"))
        s.validate_output_paths(outputs, self.tmp)
        s.save_results(outputs, gi)
        self.assertNotIn("pairing_singlet_error", gi)
        self.assertTrue(os.path.exists(os.path.join(self.tmp, "sub", "gap_bond_singlet.dat")))
        leftovers = [os.path.join(d, f) for d, _, fs in os.walk(self.tmp) for f in fs
                     if ".tmp" in f]
        self.assertEqual(leftovers, [])

    def _dir_snapshot(self, tmp):
        return {os.path.join(d, f) for d, _, fs in os.walk(tmp) for f in fs}

    def test_publication_rolls_back_over_previous_outputs(self):
        """A rename that fails partway must leave the PREVIOUS outputs of the
        channel intact -- neither a mix of new and old files nor a hole."""
        s, r = _flex({"longitudinal_bond_pairing": "singlet", "IterationMax": 3})
        gi = r.get_param("green")
        s.solve(gi, self.tmp)
        self.assertNotIn("pairing_singlet_error", gi)
        s.validate_output_paths(self._outputs(self.tmp), self.tmp)
        s.save_results(self._outputs(self.tmp), gi)
        self.assertNotIn("pairing_singlet_error", gi)
        targets = {k: os.path.join(self.tmp, fn) for k, fn in (
            ("eliashberg_bond", "eliashberg_bond_singlet.npz"),
            ("gap_bond", "gap_bond_singlet.dat"),
            ("eigenvalue_bond", "eigenvalue_bond_singlet.dat"))}
        # overwrite the three published files with distinctive sentinels so that
        # a file which is NOT rolled back would be detectably different
        sentinels = {}
        for i, (k, p) in enumerate(targets.items()):
            data = ("previous-{}-{}\n".format(k, i)).encode() * (i + 1)
            with open(p, "wb") as f:
                f.write(data)
            sentinels[k] = data
        before = self._dir_snapshot(self.tmp)
        real = os.replace
        state = {"failed": False}

        def fake_replace(src, dst):
            if dst == targets["gap_bond"] and not state["failed"]:
                state["failed"] = True
                raise OSError("disk full on the second publication")
            return real(src, dst)

        with mock.patch("hwave.solver.eliashberg_bond_io.os.replace", side_effect=fake_replace):
            s.save_results(self._outputs(self.tmp), gi)
        # every previous file is byte-identical to its sentinel (restored)
        for k, p in targets.items():
            with open(p, "rb") as f:
                self.assertEqual(f.read(), sentinels[k], k)
        # no temp or backup file survives; the directory is exactly as before
        self.assertEqual(self._dir_snapshot(self.tmp), before)
        self.assertIn("pairing_singlet_error", gi)
        self.assertIn("restored", gi["pairing_singlet_error"])

    def test_publication_failure_on_first_publish_leaves_no_partials(self):
        """No previous files: a failure on the second rename must leave NONE of
        the channel's three output files behind, and record the error."""
        s, r = _flex({"longitudinal_bond_pairing": "singlet", "IterationMax": 3})
        gi = r.get_param("green")
        s.solve(gi, self.tmp)
        self.assertNotIn("pairing_singlet_error", gi)
        s.validate_output_paths(self._outputs(self.tmp), self.tmp)
        names = ("eliashberg_bond_singlet.npz", "gap_bond_singlet.dat",
                 "eigenvalue_bond_singlet.dat")
        gap_target = os.path.join(self.tmp, "gap_bond_singlet.dat")
        real = os.replace
        state = {"failed": False}

        def fake_replace(src, dst):
            if dst == gap_target and not state["failed"]:
                state["failed"] = True
                raise OSError("disk full on the second publication")
            return real(src, dst)

        with mock.patch("hwave.solver.eliashberg_bond_io.os.replace", side_effect=fake_replace):
            s.save_results(self._outputs(self.tmp), gi)
        for fn in names:
            self.assertFalse(os.path.exists(os.path.join(self.tmp, fn)), fn)
        leftovers = [f for f in os.listdir(self.tmp)
                     if ".tmp" in f or ".bak" in f]
        self.assertEqual(leftovers, [])
        self.assertIn("pairing_singlet_error", gi)

    def test_pairing_failure_keeps_flex_outputs(self):             # 10.2.7
        from unittest import mock
        from hwave.solver import eliashberg_bond as eb, eliashberg_dynamic as ed, \
            flex_bond as fb, backend as bk
        oom = bk._oom_error_types()[0] if bk._oom_error_types() else MemoryError

        def oom_instance():
            """The backend's device out-of-memory type where a device backend
            is installed (CuPy's ``OutOfMemoryError(size, total)`` takes the
            byte counts, not a message), ``MemoryError`` otherwise."""
            try:
                return oom("boom")
            except TypeError:
                return oom(1, 1)

        real_dress = fb._dress

        def dress_fails_only_for_the_pairing(*a, **k):
            """The FLEX map and the pairing re-dressing share ``_dress``; the
            pairing passes ``iteration=None`` (argument 8), the map its own
            iteration number. Only the pairing call fails, so the FLEX solve
            itself still produces every artifact."""
            iteration = a[7] if len(a) > 7 else k.get("iteration")
            if iteration is None:
                raise ValueError("denominator singular")
            return real_dress(*a, **k)

        shared = [("accumulator", eb.PairVertexAccumulator, "add_dressed", oom_instance()),
                  ("dress", fb, "_dress", dress_fails_only_for_the_pairing)]
        per_channel = [("solver", ed, "run_leading_eigenproblem", RuntimeError("no convergence")),
                       ("outputs", ed, "write_dynamic_outputs", OSError("disk full"))]
        for label, target, attr, exc in shared + per_channel:
            tmp = tempfile.mkdtemp()
            try:
                s, r = _flex({"longitudinal_bond_pairing": "both", "IterationMax": 3})
                gi = r.get_param("green")
                with mock.patch.object(target, attr, side_effect=exc):
                    s.solve(gi, tmp)                              # never raises
                    for k in ("sigma", "green", "chiq_s", "longitudinal_bond_chi_s"):
                        self.assertIn(k, gi, label)
                    s.validate_output_paths(self._outputs(tmp), tmp)
                    s.save_results(self._outputs(tmp), gi)
                for fn in ("sigma.npz", "green.npz", "chiq_s.npz"):
                    self.assertTrue(os.path.exists(os.path.join(tmp, fn)), label)
                if label in ("accumulator", "dress"):
                    for eta in ("singlet", "triplet"):
                        self.assertIn("pairing_{}_error".format(eta), gi, label)
                        self.assertFalse(os.path.exists(
                            os.path.join(tmp, "eliashberg_bond_{}.npz".format(eta))), label)
                    self.assertEqual(gi["pairing_singlet_error"], gi["pairing_triplet_error"])
                else:
                    # every channel fails the same way for a per-channel stage patched
                    # globally; the failure is per channel and nothing of that channel
                    # is published
                    for eta in ("singlet", "triplet"):
                        self.assertIn("pairing_{}_error".format(eta), gi, label)
                        self.assertFalse(os.path.exists(
                            os.path.join(tmp, "eliashberg_bond_{}.npz".format(eta))), label)
                        self.assertFalse(os.path.exists(
                            os.path.join(tmp, "eliashberg_bond_{}.tmp.npz".format(eta))), label)
            finally:
                shutil.rmtree(tmp, ignore_errors=True)
        # per-channel isolation: patch the solver to fail on the FIRST call only
        tmp = tempfile.mkdtemp()
        try:
            s, r = _flex({"longitudinal_bond_pairing": "both", "IterationMax": 3})
            gi = r.get_param("green")
            real = ed.run_leading_eigenproblem
            calls = []

            def flaky(*a, **k):
                calls.append(1)
                if len(calls) == 1:
                    raise RuntimeError("first channel fails")
                return real(*a, **k)

            with mock.patch.object(ed, "run_leading_eigenproblem", side_effect=flaky):
                s.solve(gi, tmp)
            self.assertIn("pairing_singlet_error", gi)
            self.assertIn("pairing_triplet_eigenvalue", gi)
            s.validate_output_paths(self._outputs(tmp), tmp)
            s.save_results(self._outputs(tmp), gi)
            self.assertTrue(os.path.exists(os.path.join(tmp, "eliashberg_bond_triplet.npz")))
            self.assertFalse(os.path.exists(os.path.join(tmp, "eliashberg_bond_singlet.npz")))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_solve_entry_preflight_refusal(self):                  # 10.2.8
        s, r = _flex({"longitudinal_bond_pairing": "singlet", "IterationMax": 3})
        s._pairing_controls = s._pairing_controls.__class__(
            **dict(s._pairing_controls.__dict__, bond_memory_cap_gb=1e-9))
        gi = r.get_param("green")
        with self.assertRaisesRegex(ValueError, "matsubara_basis = 'ir'"):
            s.solve(gi, self.tmp)

    def test_flag_off_is_deterministic(self):                      # 10.2.4 (in-tree half)
        """Two runs of THIS tree with the pairing absent write byte-identical
        artifacts -- run-to-run determinism, which is what a comparison
        against a second run of the same code can prove and all it can prove.

        The backward-compatibility half of spec 10.2.4 -- identity with the
        revision this branch was cut from -- is
        :class:`TestFlagOffMatchesTheReferenceRevision`, which needs a second
        source tree and therefore skips wherever one is not provisioned. This
        test needs nothing, so the gate keeps a flag-off assertion even
        there."""
        ref = tempfile.mkdtemp()
        try:
            s1, r1 = _flex({"IterationMax": 3, "longitudinal_bond_output_full": True})
            gi1 = r1.get_param("green")
            s1.solve(gi1, ref)
            s1.validate_output_paths(self._outputs(ref), ref)
            s1.save_results(self._outputs(ref), gi1)
            s2, r2 = _flex({"IterationMax": 3, "longitudinal_bond_output_full": True})
            gi2 = r2.get_param("green")
            s2.solve(gi2, self.tmp)
            s2.validate_output_paths(self._outputs(self.tmp), self.tmp)
            s2.save_results(self._outputs(self.tmp), gi2)
            for fn in os.listdir(ref):
                a, b = os.path.join(ref, fn), os.path.join(self.tmp, fn)
                if fn.endswith(".npz"):
                    with np.load(a) as da, np.load(b) as db:
                        self.assertEqual(set(da.files), set(db.files), fn)
                        for k in da.files:
                            np.testing.assert_array_equal(da[k], db[k],
                                                          err_msg="{} {}".format(fn, k))
                else:
                    with open(a, "rb") as fa, open(b, "rb") as fb:
                        self.assertEqual(fa.read(), fb.read(), fn)
            self.assertFalse(any(fn.startswith("eliashberg_bond")
                                 for fn in os.listdir(self.tmp)))
        finally:
            shutil.rmtree(ref, ignore_errors=True)

    def test_inprocess_equals_postprocessing(self):                # 10.2.3
        import hwave.sc as sc
        flex_dir = tempfile.mkdtemp()
        try:
            _run_gate_1orb(flex_dir, pairing="both",
                           eli={"solver_mode": "eigenvalue", "num_eigenvalues": 3})
            for eta in ("singlet", "triplet"):
                with np.load(os.path.join(flex_dir,
                                          "eliashberg_bond_{}.npz".format(eta))) as d:
                    lam_ip, gap_ip, evs = float(d["eigenvalue"]), d["gap"], d["eigenvalues_all"]
                self.assertGreater(abs(evs[0] - evs[1]), 1e-3 * abs(evs[0]),
                                   "isolated leading eigenvalue")
                out = tempfile.mkdtemp()
                try:
                    lam_pp = sc.calc_eliashberg(_sc_input(flex_dir, out, 0.5, 64, (4, 4, 1),
                                                          pairing_type=eta))
                    with np.load(os.path.join(out, "gap_dynamic.npz")) as d:
                        gap_pp = d["gap"]
                finally:
                    shutil.rmtree(out, ignore_errors=True)
                self.assertAlmostEqual(lam_ip, lam_pp, delta=1e-10 * abs(lam_pp), msg=eta)
                ov = abs(np.vdot(gap_ip, gap_pp)) / (np.linalg.norm(gap_ip)
                                                     * np.linalg.norm(gap_pp))
                self.assertAlmostEqual(ov, 1.0, delta=1e-8, msg=eta)
        finally:
            shutil.rmtree(flex_dir, ignore_errors=True)


#: The revision the flag-off comparison of spec 10.2.4 is against: the merge
#: this branch was cut from (``develop`` at that point, the merge of the
#: bond-gate GPU work). It is deliberately NOT
#: ``test_flex_second_order_compat.DEVELOP_COMMIT``, which names an earlier
#: merge and belongs to the second-order comparisons: the claim here is that
#: a run WITHOUT the pairing table reproduces the tree this branch started
#: from, so the revision has to be that tree and no other.
BOND_REFERENCE_COMMIT = "65719b742c7c7d5cd3aada9d6e533c608563dcfe"

#: The environment variables that provision and require the reference
#: checkout of :data:`BOND_REFERENCE_COMMIT`. A pair of its own, separate
#: from the second-order harnesses': the two references are at different
#: revisions, so one checkout cannot serve both, and a "required" flag shared
#: between them would turn a missing bond reference into a failure of a
#: comparison whose own reference is sitting right there.
BOND_REFERENCE_PATH_ENV = "HWAVE_DEVELOP_CHECKOUT_BOND"
BOND_REFERENCE_REQUIRE_ENV = "HWAVE_REQUIRE_DEVELOP_COMPARISON_BOND"

#: The archive members this branch is allowed to add or change with the
#: pairing absent: the two bond vertices, the orbital count the loader needs
#: in order to check them, and the schema stamp that announces all three
#: (1 -> 2). Every other member of the archive -- and every other output
#: file -- must come out of both revisions bit for bit.
_ARCHIVE_DELTA = frozenset(("S_bond", "C_bond", "norb", "bond_archive_schema",
                            "longitudinal_bond_cond_tol", "bond_archive_layout"))

#: The provenance members the guard-policy work (issue #199) adds to EVERY
#: archive of a Hartree-Fock / bond-gate run, the bond archive included:
#: the requested policy and the number of tolerated violations. With the
#: default policy they are the constant pair ``("refuse", 0)`` -- asserted
#: below -- and they carry no numerical content, so their presence is not a
#: compatibility break; every member the reference wrote must still come out
#: bit for bit. ``longitudinal_bond_cond_tol`` (the conditioning floor the
#: bond archive now records, the historical ``1e-3``) is the same kind of
#: addition and is listed in :data:`_ARCHIVE_DELTA` with the bond members.
_PROVENANCE_DELTA = frozenset(("flex_guard_policy", "flex_guard_violations"))

#: The driver both revisions run: one ``qlms.run`` on an input dict read from
#: a file. A file rather than an argument because the reference tree is
#: driven through a subprocess and the dict carries paths.
_QLMS_RUN = r'''
import json, sys
import hwave.qlms as qlms
with open(sys.argv[1]) as fh:
    qlms.run(input_dict=json.load(fh))
'''


def _bond_compat_input(out_dir):
    """The run both revisions execute for the flag-off comparison: the
    single-band (U = 4, V = 1) bond gate with the full archive and NO
    ``[eliashberg]`` table -- i.e. exactly the configuration an existing user
    of the bond gate has today. ``path_to_input`` stays relative so that each
    tree reads its own copy of the fixture (the two copies are identical; a
    shared absolute path would instead compare one revision's run against a
    foreign input)."""
    return {"mode": {"mode": "FLEX", "calc_scheme": "general",
                     "param": {"T": 0.5, "filling": 0.5, "CellShape": [4, 4, 1],
                               "SubShape": [1, 1, 1], "Nmat": 16, "IterationMax": 200,
                               "Mix": 0.3, "EPS": 8, "flex_hartree_fock": True,
                               "mixing_scheme": "linear",
                               "longitudinal_bond_channels": True,
                               "longitudinal_bond_output_full": True}},
            "file": {"input": {"path_to_input": "",
                               "interaction": {"path_to_input": _IN1, "Geometry": "geom.dat",
                                               "Transfer": "transfer.dat",
                                               "CoulombIntra": "coulombintra.dat",
                                               "CoulombInter": "coulombinter.dat"}},
                     "output": {"path_to_output": out_dir, "chiq_s": "chiq_s",
                                "chiq_c": "chiq_c", "green": "green", "energy": "energy.dat",
                                "longitudinal_bond": "longitudinal_bond"}}}


class TestFlagOffMatchesTheReferenceRevision(unittest.TestCase):
    """spec 10.2.4: with the pairing absent, this branch reproduces the
    revision it was cut from (:data:`BOND_REFERENCE_COMMIT`) file by file --
    the archive's two new vertex members, ``norb``, the schema stamp and the
    recorded conditioning floor (:data:`_ARCHIVE_DELTA`), and the provenance
    pair of the guard policy (:data:`_PROVENANCE_DELTA`) excepted.

    The comparison needs a SECOND source tree, provisioned by CI as a
    detached worktree and pointed at by :data:`BOND_REFERENCE_PATH_ENV`;
    without one the test skips with the reason (and, where
    :data:`BOND_REFERENCE_REQUIRE_ENV` is set, fails instead of skipping, so
    a comparison that stops running in CI is reported rather than silently
    dropped)."""

    def setUp(self):
        self.ref_out = tempfile.mkdtemp()
        self.new_out = tempfile.mkdtemp()
        self.work = tempfile.mkdtemp()

    def tearDown(self):
        for d in (self.ref_out, self.new_out, self.work):
            shutil.rmtree(d, ignore_errors=True)

    def _run(self, checkout, out_dir, name):
        import json
        import subprocess
        import sys
        from tests.test_flex_second_order_compat import reference_subprocess_env
        path = os.path.join(self.work, "input_{}.json".format(name))
        with open(path, "w") as fh:
            json.dump(_bond_compat_input(out_dir), fh)
        if checkout == os.getcwd():
            # this tree's own run: keep it measured by the coverage run
            env = dict(os.environ,
                       PYTHONPATH=os.path.join(checkout, "src") + ":" + checkout)
        else:
            env = reference_subprocess_env(checkout)
        subprocess.run([sys.executable, "-B", "-c", _QLMS_RUN, path],
                       env=env, check=True, capture_output=True, cwd=checkout)

    @staticmethod
    def _identical(a, b):
        return a.shape == b.shape and a.dtype == b.dtype and np.array_equal(a, b)

    def test_outputs_match_the_reference_revision(self):           # 10.2.4
        from tests.test_flex_second_order_compat import develop_checkout
        dev, why = develop_checkout(expected_commit=BOND_REFERENCE_COMMIT,
                                    path_env=BOND_REFERENCE_PATH_ENV,
                                    require_env=BOND_REFERENCE_REQUIRE_ENV)
        if dev is None:
            raise unittest.SkipTest(why)
        self._run(dev, self.ref_out, "reference")
        self._run(os.getcwd(), self.new_out, "branch")
        ref_files = sorted(os.listdir(self.ref_out))
        self.assertEqual(ref_files, sorted(os.listdir(self.new_out)),
                         "the two revisions wrote different file sets")
        self.assertIn("longitudinal_bond.npz", ref_files)
        self.assertIn("energy.dat", ref_files)
        for fn in ref_files:
            a, b = os.path.join(self.ref_out, fn), os.path.join(self.new_out, fn)
            if not fn.endswith(".npz"):
                with open(a, "rb") as fa, open(b, "rb") as fb:
                    self.assertEqual(fa.read(), fb.read(), fn)
                continue
            with np.load(a) as da, np.load(b) as db:
                if fn != "longitudinal_bond.npz":
                    # every member the reference wrote, bit for bit; the only
                    # members this side may add are the provenance pair
                    self.assertEqual(set(da.files) - set(db.files), set(),
                                     "{}: members the reference wrote are missing".format(fn))
                    added = set(db.files) - set(da.files)
                    # ... plus, where the bond members are written (the
                    # chiq_s / chiq_c archives carry them too), the recorded
                    # conditioning floor
                    self.assertTrue(_PROVENANCE_DELTA <= added, (fn, added))
                    self.assertTrue(added <= _PROVENANCE_DELTA | {"longitudinal_bond_cond_tol"},
                                    (fn, added))
                    for k in da.files:
                        np.testing.assert_array_equal(da[k], db[k],
                                                      err_msg="{} {}".format(fn, k))
                    self._assert_default_provenance(db, fn)
                    if "longitudinal_bond_cond_tol" in added:
                        self.assertEqual(float(db["longitudinal_bond_cond_tol"]), 1.0e-3, fn)
                    continue
                # the one file this branch may change, and only in the named
                # members
                self.assertEqual(set(da.files) - set(db.files), set(),
                                 "the archive dropped members the reference wrote")
                added = set(db.files) - set(da.files)
                changed = {k for k in da.files if not self._identical(da[k], db[k])}
                self.assertEqual(added | changed, set(_ARCHIVE_DELTA | _PROVENANCE_DELTA))
                for k in set(da.files) - _ARCHIVE_DELTA:
                    np.testing.assert_array_equal(da[k], db[k],
                                                  err_msg="{} {}".format(fn, k))
                self.assertEqual(int(da["bond_archive_schema"]), 1)
                self.assertEqual(int(db["bond_archive_schema"]), 2)
                self.assertEqual(float(db["longitudinal_bond_cond_tol"]), 1.0e-3)
                self._assert_default_provenance(db, fn)

    def _assert_default_provenance(self, db, fn):
        """The provenance pair of issue #199 at its default: the reference
        revision has no guard policy, so this side must report the refusing
        policy and no tolerated violation -- anything else would mean the
        comparison ran under a different guard than the reference did."""
        self.assertEqual(str(db["flex_guard_policy"]), "refuse", fn)
        self.assertEqual(int(db["flex_guard_violations"]), 0, fn)
