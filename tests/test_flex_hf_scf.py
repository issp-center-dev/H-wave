"""The standalone Hartree-Fock-renormalised FLEX loop (spec 2026-09-06
sections 1, 2, 4.2-4.3, 5.6; gates G0-off, G5, D7): gate-off inputs are
byte-identical to develop, the single-band Hubbard HF term is a chemical-
potential shift, the mean-field seed follows the initial-value
trajectory, the W-disabled fixed point equals paramagnetic UHFk, the
reset set, IterationMax = 0 and the provenance block."""
import os
import subprocess
import sys
import tempfile
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k

_IN2 = "tests/rpa/input_2orb"
_IN1 = "tests/rpa/input"


def _build(param_extra=None, path=_IN2, interactions=None, mixing="linear"):
    import hwave.solver.flex as flex_mod
    idict = {"path_to_input": path, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
    idict.update({"CoulombInter": "coulombinter.dat"} if interactions is None else interactions)
    r = read_input_k.QLMSkInput({"path_to_input": path, "interaction": idict})
    par = {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1],
           "Nmat": 32, "IterationMax": 4, "Mix": 0.5, "EPS": 8, "mixing_scheme": mixing}
    par.update(param_extra or {})
    info = {"mode": "FLEX", "param": par, "enable_spin_orbital": False, "calc_scheme": "general"}
    return flex_mod.FLEX(r.get_param("ham"), {}, info), r


_G0_SCRIPT = r'''
import os, sys, json, logging, numpy as np
sys.path.insert(0, sys.argv[1] + "/src"); sys.path.insert(0, sys.argv[1])
import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.flex as flex_mod
out = sys.argv[2]
logging.basicConfig(level=logging.DEBUG, filename=os.path.join(out, "log.txt"), filemode="w",
                    format="%(name)s %(levelname)s %(message)s")
p = "tests/rpa/input_2orb"
idict = {"path_to_input": p, "Geometry": "geom.dat", "Transfer": "transfer.dat", "CoulombInter": "coulombinter.dat"}
r = read_input_k.QLMSkInput({"path_to_input": p, "interaction": idict})
par = {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1], "Nmat": 32,
       "IterationMax": 3, "Mix": 0.5, "EPS": 8, "mixing_scheme": "anderson"}
s = flex_mod.FLEX(r.get_param("ham"), {}, {"mode": "FLEX", "param": par, "enable_spin_orbital": False, "calc_scheme": "general"})
gi = r.get_param("green")
s.solve(gi, out)
s.save_results({"path_to_output": out, "sigma": "sigma.npz", "green": "green.npz", "chi0q": "chi0q.npz", "chiq": "chiq.npz"}, gi)
'''


class TestG0Off(unittest.TestCase):
    """Gate-off inputs: every archive member and the complete log stream
    equal those of develop (the main checkout)."""

    def test_byte_identity_against_develop(self):
        here = os.getcwd()
        develop = os.environ.get("HWAVE_DEVELOP_CHECKOUT",
                                 os.path.abspath(os.path.join(here, "..", "..", "..")))
        if not os.path.isdir(os.path.join(develop, "src", "hwave")):
            self.skipTest("develop checkout not found at {}".format(develop))
        with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
            for root, out in ((here, a), (develop, b)):
                subprocess.run([sys.executable, "-B", "-c", _G0_SCRIPT, root, out], check=True, cwd=here)
            for name in ("sigma.npz", "green.npz", "chi0q.npz", "chiq.npz"):
                da, db = np.load(os.path.join(a, name), allow_pickle=True), np.load(os.path.join(b, name), allow_pickle=True)
                self.assertEqual(da.files, db.files, name)
                for k in da.files:
                    self.assertEqual(da[k].dtype, db[k].dtype, (name, k))
                    self.assertTrue(np.array_equal(da[k], db[k]), (name, k))
            la = open(os.path.join(a, "log.txt")).read().replace(here, "<root>").replace(a, "<out>")
            lb = open(os.path.join(b, "log.txt")).read().replace(develop, "<root>").replace(b, "<out>")
            self.assertEqual(la, lb)


class TestStandaloneHF(unittest.TestCase):

    def test_single_band_hubbard_hf_is_a_mu_shift(self):
        """norb=1, on-site U only: Sigma_HF is a constant absorbed by mu, so
        the dressed Green function equals the HF-off one."""
        s_on, r1 = _build({"flex_hartree_fock": True}, path=_IN1, interactions={"CoulombIntra": "coulombintra.dat"})
        s_off, r2 = _build({}, path=_IN1, interactions={"CoulombIntra": "coulombintra.dat"})
        g_on, g_off = r1.get_param("green"), r2.get_param("green")
        os.makedirs("tests/flex/output", exist_ok=True)
        s_on.solve(g_on, "tests/flex/output"); s_off.solve(g_off, "tests/flex/output")
        st = np.asarray(g_on["sigma_static"])
        self.assertLess(np.abs(st - st[0, 0, 0, 0, 0]).max(), 1e-12)     # k-independent constant
        np.testing.assert_allclose(np.asarray(g_on["green"]), np.asarray(g_off["green"]), rtol=0, atol=1e-10)
        self.assertAlmostEqual(g_on["physics"]["mu"] - g_off["physics"]["mu"], st[0, 0, 0, 0, 0].real, places=9)
        self.assertAlmostEqual(g_on["physics"]["NCond"], g_off["physics"]["NCond"], places=8)

    def test_d7_trajectory_delta_h_is_never_readded(self):
        """A trans_mod seed Delta H with the HF map forced to zero and
        linear mixing m: the static state is (1-m)^k Delta H after k maps."""
        import hwave.solver.flex_hf as flex_hf
        s, r = _build({"flex_hartree_fock": True, "IterationMax": 2, "Mix": 0.25})
        nvol, norb = s.lattice.nvol, s.norb
        shift = np.diag([0.3, -0.2])
        H = np.zeros((nvol, 2, norb, 2, norb), complex)
        base = np.asarray(s.ham_info.ham_trans_q).reshape(nvol, norb, norb)
        H[:, 0, :, 0, :] = base + shift; H[:, 1, :, 1, :] = base + shift
        gi = r.get_param("green"); gi["trans_mod"] = H.reshape(nvol, 2 * norb, 2 * norb)
        real = flex_hf.hf_map
        flex_hf.hf_map = lambda rho_r, tables, shape, norb_: np.zeros((rho_r.shape[0], norb_, norb_), complex)
        seen = []
        s._iteration_hook = lambda d: seen.append(d)
        try:
            s.solve(gi, "tests/flex/output")
        finally:
            flex_hf.hf_map = real
        self.assertEqual(len(seen), 2)
        self.assertEqual(np.abs(seen[0]["static_new"]).max(), 0.0)              # the map is zero
        np.testing.assert_allclose(np.asarray(gi["sigma_static"])[0, 0],
                                   0.75 ** 2 * np.broadcast_to(shift, (nvol, norb, norb)), rtol=0, atol=1e-12)
        self.assertIn("trans_mod", gi)                                          # input preserved

    def test_reset_set_iteration_max_zero_and_provenance(self):
        s, r = _build({"flex_hartree_fock": True})
        gi = r.get_param("green")
        os.makedirs("tests/flex/output", exist_ok=True)
        s.solve(gi, "tests/flex/output")
        self.assertEqual(s.scf_iterations, 4)
        for k in ("sigma", "sigma_static", "sigma_fluct", "green", "physics", "chi0q", "chiq_s", "chiq_c"):
            self.assertIn(k, gi)
        np.testing.assert_allclose(gi["sigma"], gi["sigma_static"] + gi["sigma_fluct"], atol=1e-14)
        # reuse: a second solver on the same container with IterationMax = 0
        s0, _ = _build({"flex_hartree_fock": True, "IterationMax": 0})
        gi["chiq_s"] = np.zeros(3)                                             # stale, must vanish
        s0.solve(gi, "tests/flex/output")
        for k in ("chi0q", "chiq_s", "chiq_c"):
            self.assertNotIn(k, gi)
        for k in ("sigma", "green", "physics"):
            self.assertIn(k, gi)
        self.assertEqual((s0.scf_iterations, s0.map_iteration, s0.state_iteration), (0, 0, 0))
        self.assertTrue(np.isnan(s0.scf_green_residual))
        with tempfile.TemporaryDirectory() as tmp:
            s0.save_results({"path_to_output": tmp, "sigma": "s.npz", "green": "g.npz", "chiq": "c.npz"}, gi)
            d = np.load(os.path.join(tmp, "s.npz"))
            for k in ("sigma_convention", "sigma_static", "sigma_fluct", "scf_converged", "scf_iterations",
                      "map_iteration", "state_iteration", "scf_sigma_residual", "scf_green_residual",
                      "scf_component_residual", "payload_kind", "hf_density_error", "hf_density_source",
                      "density_target_enforced"):
                self.assertIn(k, d.files, k)
            self.assertEqual(str(d["payload_kind"]), "final_state")
            self.assertEqual(str(d["sigma_convention"]), "split")
            self.assertFalse(os.path.exists(os.path.join(tmp, "c.npz")) and "chiq_s" in np.load(os.path.join(tmp, "c.npz")).files)
        # failure atomicity: a refused re-solve leaves nothing produced behind
        s_bad, _ = _build({"flex_hartree_fock": True})
        gi["trans_mod"] = np.zeros((s.lattice.nvol, 2 * s.norb, 2 * s.norb), complex)
        gi["sigma_init_envelope"] = s._read_sigma.__self__ and None
        gi.pop("sigma_init_envelope")
        import hwave.solver.flex_hf as flex_hf
        env = flex_hf.SigmaSeedEnvelope(sigma=np.zeros(1), sigma_static=None, sigma_fluct=None,
                                        marker="split", ir_meta=None, file_name="x")
        gi["sigma_init_envelope"] = env
        with self.assertRaises(ValueError):
            s_bad.solve(gi, "tests/flex/output")
        for k in ("sigma", "green", "physics", "chi0q"):
            self.assertNotIn(k, gi)
        self.assertFalse(hasattr(s_bad, "physics"))

    def test_hf_off_does_not_write_new_members(self):
        s, r = _build({})
        gi = r.get_param("green")
        os.makedirs("tests/flex/output", exist_ok=True)
        s.solve(gi, "tests/flex/output")
        with tempfile.TemporaryDirectory() as tmp:
            s.save_results({"path_to_output": tmp, "sigma": "s.npz", "chi0q": "c.npz"}, gi)
            self.assertNotIn("sigma_convention", np.load(os.path.join(tmp, "s.npz")).files)
            self.assertNotIn("payload_kind", np.load(os.path.join(tmp, "c.npz")).files)


class TestG5UHFkFixedPoint(unittest.TestCase):

    def test_w_disabled_fixed_point_equals_paramagnetic_uhfk(self):
        """Fluctuations disabled: the converged FLEX Hartree-Fock state equals
        a paramagnetic UHFk run (flag_fock=true, 2Sz=0, same T, per-cell
        normalisation) at the fixed point."""
        import hwave.solver.flex as flex_mod
        from hwave.solver.uhfk import UHFk
        # FLEX: seed sigma_static = HF map of the bare paramagnetic density (through trans_mod)
        s, r = _build({"flex_hartree_fock": True, "IterationMax": 400, "Mix": 0.5, "EPS": 10},
                      path=_IN2, interactions={"CoulombInter": "onsite_inter.dat", "CoulombIntra": "coulombintra.dat"})
        real = flex_mod.FLEX._calc_self_energy_general
        flex_mod.FLEX._calc_self_energy_general = lambda self_, g, v, b: np.zeros_like(g)
        gi = r.get_param("green")
        try:
            s.solve(gi, "tests/flex/output")
        finally:
            flex_mod.FLEX._calc_self_energy_general = real
        self.assertTrue(s.scf_converged)
        # UHFk on the same input, paramagnetic (2Sz = 0), T = 2.0, Ncond = filling * 2 * norb * nvol
        nvol, norb = s.lattice.nvol, s.norb
        idict = {"path_to_input": _IN2, "Geometry": "geom.dat", "Transfer": "transfer.dat",
                 "CoulombInter": "onsite_inter.dat", "CoulombIntra": "coulombintra.dat"}
        r2 = read_input_k.QLMSkInput({"path_to_input": _IN2, "interaction": idict})
        par = {"T": 2.0, "Ncond": int(round(0.5 * 2 * norb * nvol)), "2Sz": 0, "flag_fock": True,
               "CellShape": [4, 4, 1], "SubShape": [1, 1, 1], "IterationMax": 2000, "EPS": 12, "Mix": 0.5,
               "RndSeed": 1}
        u = UHFk(r2.get_param("ham"), {"print_level": 1, "print_step": 1}, {"mode": "UHFk", "param": par})
        with tempfile.TemporaryDirectory() as tmp:
            u.solve(r2.get_param("green"), tmp)
        G = u.Green                                                   # (nvol, 2, norb, 2, norb)
        # the FLEX final-state density
        from hwave.solver import hartree_fock as hf
        heff = hf.heff_eigenpairs(np.asarray(s.H0_k)[0], gi["sigma_static"][0, 0])
        dens = hf.equal_time_density(gi["green"], heff, s.mu, 1.0 / 2.0, (4, 4, 1))
        np.testing.assert_allclose(G[:, 0, :, 0, :], dens.rho_r, rtol=0, atol=1e-8)
        np.testing.assert_allclose(G[:, 1, :, 1, :], dens.rho_r, rtol=0, atol=1e-8)
        self.assertLess(np.abs(G[:, 0, :, 1, :]).max(), 1e-10)
        h_uhf = np.asarray(u.ham).reshape(nvol, 2, norb, 2, norb)[:, 0, :, 0, :]
        h_flex = np.asarray(s.H0_k)[0] + gi["sigma_static"][0, 0]
        np.testing.assert_allclose(h_flex, h_uhf, rtol=0, atol=1e-8)


if __name__ == "__main__":
    unittest.main()
