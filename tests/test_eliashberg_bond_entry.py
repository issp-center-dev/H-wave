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
        ``matsubara_basis = "ir"``) is still not a lambda comparison: see
        :meth:`test_ir_arm_parity_leakage_converges` below for what the IR arm
        does assert on this same working point, and why.
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

    @heavy
    def test_ir_arm_parity_leakage_converges(self):         # 10.2.2 (IR half)
        """The IR arm of 10.2.2 on the SAME real bond-gate archives.

        The flat-term fix (the instantaneous vertex multiplies the
        unregularized Matsubara sum, i.e. the midpoint of the tau = 0 jump,
        not ``F(0^+)``) makes the IR kernel algebra exact: on an exactly
        IR-representable fixture the leakage is 1e-12 with both flat terms
        live (``test_eliashberg_bond_kernel.TestKernelIR
        .test_ir_parity_commutation``). On REAL output it drops by two orders
        of magnitude -- singlet 2.80e-01 -> 7.56e-03 and triplet 1.69e-01 ->
        9.70e-03 at Nmat 64, 2.98e-01 -> 1.29e-03 and 1.29e-01 -> 1.68e-03 at
        Nmat 128 -- and the remainder CONVERGES as roughly ``Nmat^-2``
        (3.1e-04 at Nmat 256), which the O(1), Nmat-independent flat-term
        defect did not.

        What is still NOT asserted, and why there is no ``lambda_ir`` vs
        ``lambda_uniform`` comparison here: that remainder is above the
        entry's fixed ``parity_leakage_policy = "refuse"`` threshold of 1e-8,
        so the IR run is refused and produces no eigenvalue. It is NOT the
        flat term -- forcing the instantaneous vertex to zero AND dropping the
        retained constant leaves the same floor (7.69e-03 at Nmat 64, 2.40e-03
        at Nmat 128) -- but the IR representability of the uniform-FFT bond
        archive itself, and no reachable ``ir_wmax`` / ``ir_tol`` / ``Nmat``
        setting brings it to 1e-8 (measured over wmax in {auto, 40, 200},
        ir_tol in {1e-8, 1e-12}, Nmat in {64, 128, 256}: best 2.45e-04). That
        is a separate, pre-existing data-quality issue of the archive, not of
        the pairing kernel.
        """
        try:
            import sparse_ir  # noqa: F401
        except ImportError:
            raise unittest.SkipTest("sparse-ir not installed")
        import re
        import hwave.sc as sc

        def ir_run(flex_dir, out, nmat, eta, **extra):
            # ir_fit_tol = 0: on this real output the componentwise IR fit
            # residual WITH the constant retained is 1.609e-01 (Nmat 64) /
            # 1.568e-01 (Nmat 128), above the default ir_fit_tol = 0.1. That
            # default is deliberately NOT changed -- retaining the constant
            # already brings the residual down from 1.64 / 0.39, and the rest
            # is the same archive representability the docstring describes.
            kw = dict(pairing_type=eta, matsubara_basis="ir", ir_tol=1e-8,
                      ir_keep_static_chi=True, ir_fit_tol=0.0)
            kw.update(extra)
            return sc.calc_eliashberg(_sc_input(flex_dir, out, 0.5, nmat,
                                                (4, 4, 1), **kw))

        leak = {}
        for nmat in (64, 128):
            flex_dir = tempfile.mkdtemp()
            try:
                _run_gate_1orb(flex_dir, nmat=nmat)
                for eta in ("singlet", "triplet"):
                    out = tempfile.mkdtemp()
                    try:
                        with self.assertRaisesRegex(
                                ValueError, "cross-sector leakage") as cm:
                            ir_run(flex_dir, out, nmat, eta)
                        m = re.search(r"cross-sector leakage ([0-9.eE+-]+)",
                                      str(cm.exception))
                        self.assertIsNotNone(m)
                        leak[(eta, nmat)] = float(m.group(1))
                    finally:
                        shutil.rmtree(out, ignore_errors=True)
                if nmat == 64:
                    # and the default ir_fit_tol refuses EARLIER, with the
                    # measured residual named in the message
                    out = tempfile.mkdtemp()
                    try:
                        with self.assertRaisesRegex(
                                ValueError, "exceeds ir_fit_tol") as cm:
                            ir_run(flex_dir, out, 64, "singlet", ir_fit_tol=0.1)
                        m = re.search(r"relative residual ([0-9.eE+-]+)",
                                      str(cm.exception))
                        self.assertIsNotNone(m)
                        self.assertLess(float(m.group(1)), 0.5)   # was 1.64
                    finally:
                        shutil.rmtree(out, ignore_errors=True)
            finally:
                shutil.rmtree(flex_dir, ignore_errors=True)
        for eta in ("singlet", "triplet"):
            # two orders of magnitude below the pre-fix O(0.1) flat-term
            # leakage, and converging in the FLEX run's own Nmat
            self.assertLess(leak[(eta, 64)], 5e-2, eta)
            self.assertLess(leak[(eta, 128)], 0.6 * leak[(eta, 64)],
                            "{}: no Nmat convergence {:.3e} -> {:.3e}"
                            .format(eta, leak[(eta, 64)], leak[(eta, 128)]))


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
                          "bond_delta_r"):
                    self.assertIn(k, d.files, k)
                self.assertEqual(str(d["pairing_type"]), eta)
            with open(os.path.join(self.tmp, "eigenvalue_bond_{}.dat".format(eta))) as f:
                txt = f.read()
            self.assertIn("scf_converged=", txt)
            self.assertIn("state=", txt)
            self.assertTrue(os.path.exists(os.path.join(self.tmp, "gap_bond_{}.dat".format(eta))))
            for stem in ("eliashberg_bond_{}.tmp.npz", "gap_bond_{}.tmp.dat",
                         "eigenvalue_bond_{}.tmp.dat"):
                self.assertFalse(os.path.exists(os.path.join(self.tmp, stem.format(eta))))
        from hwave.solver import flex_bond
        with self.assertRaisesRegex(ValueError, "resolve to the same file"):
            flex_bond.resolve_output_paths(
                dict(self._outputs(self.tmp), eliashberg_bond_singlet="chiq_s"),
                self.tmp, ("chiq_s", "eliashberg_bond_singlet"))

    def test_pairing_failure_keeps_flex_outputs(self):             # 10.2.7
        from unittest import mock
        from hwave.solver import eliashberg_bond as eb, eliashberg_dynamic as ed, \
            flex_bond as fb, backend as bk
        oom = bk._oom_error_types()[0] if bk._oom_error_types() else MemoryError
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

        shared = [("accumulator", eb.PairVertexAccumulator, "add_dressed", oom("boom")),
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

    def test_flag_off_regression(self):                            # 10.2.4
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

    @heavy
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
