"""Entries of the bond-resolved dynamic pairing kernel (spec 10.2)."""
import os
import shutil
import tempfile
import unittest

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


# ---------------------------------------------------------------------------
# Post-processing entry through hwave_sc (spec 6, 10.2.1 / 10.2.2)
# ---------------------------------------------------------------------------

_IN1 = "tests/rpa/input"


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


class TestPostProcessingRuns(unittest.TestCase):
    @heavy
    def test_uniform_lambda_and_outputs(self):              # 10.2.2 (uniform half)
        """A real bond-gate FLEX run, then ``hwave_sc`` on its archive for both
        pairing channels: the leading eigenvalue is a finite real number, the
        documented file set is written with the bond provenance, and the answer
        is stable in the FLEX run's own Nmat (64 vs 128, the discretization the
        post-processing inherits rather than one it chooses).

        The IR half of spec 10.2.2 (the same lambda from
        ``matsubara_basis = "ir"``) is NOT asserted here: the IR dynamic kernel
        does not commute with the combined parity once the instantaneous vertex
        term is present, so every IR run is refused by the
        ``parity_leakage_policy = "refuse"`` this entry is required to use.
        Measured on an exactly IR-representable fixture (IR fit residual
        5.8e-10, so not a data-quality effect): leakage 3.5e-12 without the
        instantaneous term and 3.52e-01 with it, unchanged from Nmat 128 to
        256, against 1e-15 for the uniform kernel in both cases. That is a
        property of the shared IR flat-term evaluation, not of the bond path
        (the on-site dynamic IR kernel leaks 3.3e-01 the same way), and it has
        to be fixed there before the comparison can be asserted.
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
            self.assertAlmostEqual(lam[(eta, 128)], lam[(eta, 64)],
                                   delta=2e-2 * abs(lam[(eta, 64)]), msg=eta)
