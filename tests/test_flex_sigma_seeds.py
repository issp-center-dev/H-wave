"""Self-energy seeds of the Phase B design (spec 2026-09-06 section 2.3,
D7, D10, D15): the seed envelope, the split-file validation, the copied
mean-field assembly, and the seed matrix under HF off / HF on."""
import os
import tempfile
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k

_IN = "tests/rpa/input_2orb"


def _build(param_extra=None, calc_scheme="general"):
    import hwave.solver.flex as flex_mod
    idict = {"path_to_input": _IN, "Geometry": "geom.dat", "Transfer": "transfer.dat",
             "CoulombInter": "coulombinter.dat"}
    r = read_input_k.QLMSkInput({"path_to_input": _IN, "interaction": idict})
    par = {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1],
           "Nmat": 32, "IterationMax": 1, "Mix": 1.0, "EPS": 1}
    par.update(param_extra or {})
    info = {"mode": "FLEX", "param": par, "enable_spin_orbital": False, "calc_scheme": calc_scheme}
    return flex_mod.FLEX(r.get_param("ham"), {}, info), r


class _SeedFiles:
    """A legacy sigma.npz written by a one-iteration FLEX run, plus split /
    total / unmarked variants derived from it."""

    def __init__(self, tmp):
        self.tmp = tmp
        s, r = _build({})
        gi = r.get_param("green")
        s.solve(gi, tmp)
        s.save_results({"path_to_output": tmp, "sigma": "legacy.npz"}, gi)
        self.legacy = os.path.join(tmp, "legacy.npz")
        d = dict(np.load(self.legacy, allow_pickle=True))
        self.sigma = d["sigma"]
        self.meta = {k: v for k, v in d.items() if k != "sigma"}
        rng = np.random.default_rng(0)
        st = rng.normal(size=self.sigma.shape[2:]) + 1j * rng.normal(size=self.sigma.shape[2:])
        st = 0.5 * (st + np.conj(np.swapaxes(st, -1, -2)))
        self.sigma_static = st[None, None]                          # (1, 1, nvol, norb, norb)
        self.sigma_fluct = self.sigma - self.sigma_static
        self.split = os.path.join(tmp, "split.npz")
        np.savez(self.split, sigma=self.sigma, sigma_static=self.sigma_static,
                 sigma_fluct=self.sigma_fluct, sigma_convention="split", **self.meta)
        self.total = os.path.join(tmp, "total.npz")
        np.savez(self.total, sigma=self.sigma, sigma_convention="total", **self.meta)

    def variant(self, name, **changes):
        d = dict(np.load(self.split, allow_pickle=True))
        d.update(changes)
        path = os.path.join(self.tmp, name)
        np.savez(path, **d)
        return path


class TestSeedEnvelope(unittest.TestCase):

    def test_read_sigma_returns_an_unvalidated_read_only_envelope(self):
        with tempfile.TemporaryDirectory() as tmp:
            f = _SeedFiles(tmp)
            s, _ = _build({})
            env = s._read_sigma(f.split)
            self.assertEqual(env.marker, "split")
            self.assertEqual(env.file_name, f.split)
            self.assertFalse(env.sigma.flags.writeable)
            self.assertFalse(env.sigma_static.flags.writeable)
            self.assertEqual(env.sigma_static.shape, (1, 1, 16, 2, 2))
            env = s._read_sigma(f.legacy)
            self.assertIsNone(env.marker)
            self.assertIsNone(env.sigma_static)
            env = s._read_sigma(f.total)
            self.assertEqual(env.marker, "total")

    def test_read_init_keeps_the_legacy_array_and_adds_the_envelope(self):
        with tempfile.TemporaryDirectory() as tmp:
            f = _SeedFiles(tmp)
            s, _ = _build({})
            info = s.read_init({"path_to_input": tmp, "sigma_init": "split.npz"})
            np.testing.assert_array_equal(info["sigma_init"], f.sigma)
            self.assertEqual(info["sigma_init_envelope"].marker, "split")

    def test_validate_split_seed_rules(self):
        from hwave.solver import flex_hf
        with tempfile.TemporaryDirectory() as tmp:
            f = _SeedFiles(tmp)
            s, _ = _build({})
            ok = s._read_sigma(f.split)
            st, fl = flex_hf.validate_split_seed(ok, f.sigma.shape)
            np.testing.assert_array_equal(st, f.sigma_static)
            np.testing.assert_array_equal(fl, f.sigma_fluct)
            bad_sum = f.variant("bad_sum.npz", sigma_fluct=f.sigma_fluct + 1e-6)
            bad_shape = f.variant("bad_shape.npz", sigma_static=f.sigma_static[:, :, :8])
            nonherm = f.sigma_static.copy(); nonherm[0, 0, :, 0, 1] += 0.1
            bad_herm = f.variant("bad_herm.npz", sigma_static=nonherm, sigma_fluct=f.sigma - nonherm)
            bad_freq = f.variant("bad_freq.npz", sigma_static=np.repeat(f.sigma_static, 2, axis=1))
            nan = f.sigma_fluct.copy(); nan[0, 0, 0, 0, 0] = np.nan
            bad_nan = f.variant("bad_nan.npz", sigma_fluct=nan)
            unknown = f.variant("unknown.npz", sigma_convention="weird")
            for path in (bad_sum, bad_shape, bad_herm, bad_freq, bad_nan, unknown):
                with self.subTest(path=os.path.basename(path)):
                    with self.assertRaises(ValueError):
                        flex_hf.validate_split_seed(s._read_sigma(path), f.sigma.shape)


class TestMeanFieldSeed(unittest.TestCase):

    def test_calc_trans_mod_copy_matches_and_does_not_mutate(self):
        s, r = _build({})
        nd = 2 * s.norb
        rng = np.random.default_rng(2)
        g0 = rng.normal(size=(1, nd, nd)) + 1j * rng.normal(size=(1, nd, nd))
        g0 = 0.5 * (g0 + np.conj(np.swapaxes(g0, -1, -2)))
        before_r = s.ham_info.ham_trans_r.copy()
        before_q = s.ham_info.ham_trans_q.copy()
        ref = s._calc_trans_mod(g0.copy())            # the legacy helper mutates ham_trans_r
        s.ham_info.ham_trans_r[...] = before_r
        got = s._calc_trans_mod_copy(g0)
        np.testing.assert_array_equal(got, ref)
        np.testing.assert_array_equal(s.ham_info.ham_trans_r, before_r)
        np.testing.assert_array_equal(s.ham_info.ham_trans_q, before_q)

    def _paramagnetic_trans_mod(self, s, shift):
        nvol, norb = s.lattice.nvol, s.norb
        H = np.zeros((nvol, 2, norb, 2, norb), complex)
        base = np.asarray(s.ham_info.ham_trans_q).reshape(nvol, norb, norb)
        H[:, 0, :, 0, :] = base + shift
        H[:, 1, :, 1, :] = base + shift
        return H.reshape(nvol, 2 * norb, 2 * norb)

    def test_seed_matrix_hf_on(self):
        from hwave.solver import flex_hf
        with tempfile.TemporaryDirectory() as tmp:
            f = _SeedFiles(tmp)
            s, r = _build({"flex_hartree_fock": True})
            nvol, norb, nmat = s.lattice.nvol, s.norb, s.nmat
            shift = np.diag([0.3, -0.2])
            expected = (1, nmat, nvol, norb, norb)
            # none: zero
            st, fl = s._assemble_static_seed({}, expected)
            self.assertEqual(st.shape, (1, 1, nvol, norb, norb)); self.assertEqual(np.abs(st).max(), 0.0)
            self.assertEqual(np.abs(fl).max(), 0.0)
            # trans_mod -> Delta H (bare band stays bare)
            gi = {"trans_mod": self._paramagnetic_trans_mod(s, shift)}
            st, fl = s._assemble_static_seed(gi, expected)
            np.testing.assert_allclose(st[0, 0], np.broadcast_to(shift, (nvol, norb, norb)), atol=1e-12)
            self.assertIn("trans_mod", gi)                       # inputs preserved
            # split seed: components as stored
            env = s._read_sigma(f.split)
            st, fl = s._assemble_static_seed({"sigma_init_envelope": env}, expected)
            np.testing.assert_array_equal(st, f.sigma_static); np.testing.assert_array_equal(fl, f.sigma_fluct)
            # split + mean field: refused (D15)
            with self.assertRaises(ValueError) as cm:
                s._assemble_static_seed({"sigma_init_envelope": env, "trans_mod": gi["trans_mod"]}, expected)
            self.assertIn("trans_mod", str(cm.exception))
            # total / unmarked: refused naming the converter
            for path in (f.total, f.legacy):
                with self.assertRaises(ValueError) as cm:
                    s._assemble_static_seed({"sigma_init_envelope": s._read_sigma(path)}, expected)
                self.assertIn("hwave_sigma_split", str(cm.exception))
            # spin-diag trans_mod: refused at step 4 (spin classification)
            H = gi["trans_mod"].reshape(nvol, 2, norb, 2, norb).copy(); H[:, 1, :, 1, :] += 0.5
            with self.assertRaises(ValueError) as cm:
                s._assemble_static_seed({"trans_mod": H.reshape(nvol, 2 * norb, 2 * norb)}, expected)
            self.assertIn("spin", str(cm.exception))
            # trans_mod and green_init: trans_mod wins (INFO)
            g0 = np.zeros((1, 2 * norb, 2 * norb), complex)
            with self.assertLogs("hwave.solver.flex", level="INFO") as lg:
                st2, _ = s._assemble_static_seed({"trans_mod": gi["trans_mod"], "green_init": g0}, expected)
            np.testing.assert_array_equal(st2, st if False else s._assemble_static_seed(gi, expected)[0])
            self.assertTrue(any("trans_mod" in m and "green_init" in m for m in lg.output))

    def test_hf_off_paths_unchanged(self):
        """With HF off the legacy array contract is what the loop consumes:
        a total seed is one array, a split seed's total is used."""
        with tempfile.TemporaryDirectory() as tmp:
            f = _SeedFiles(tmp)
            s, _ = _build({})
            info = s.read_init({"path_to_input": tmp, "sigma_init": "total.npz"})
            np.testing.assert_array_equal(info["sigma_init"], f.sigma)
            info = s.read_init({"path_to_input": tmp, "sigma_init": "split.npz"})
            np.testing.assert_array_equal(info["sigma_init"], f.sigma)


if __name__ == "__main__":
    unittest.main()
