"""``hwave_sigma_split`` (spec 2026-09-06 section 2.3): argument coupling,
the static.npz contract, --zero-static, the UHFk trans_mod recipe against
the solver's own mean-field seed path, --uhfk-spin-major, overwrite
refusal, marker admissibility and the output's validity as a Phase B seed."""
import os
import tempfile
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k

_IN2 = "tests/rpa/input_2orb"
_INTER = {"CoulombInter": "onsite_inter.dat", "CoulombIntra": "coulombintra.dat"}
_PAR = {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1], "Nmat": 8,
        "IterationMax": 1, "Mix": 1.0, "EPS": 1}


def _reader():
    idict = {"path_to_input": _IN2, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
    idict.update(_INTER)
    return read_input_k.QLMSkInput({"path_to_input": _IN2, "interaction": idict})


def _flex(hf, extra=None):
    import hwave.solver.flex as flex_mod
    r = _reader()
    par = dict(_PAR, flex_hartree_fock=hf)
    par.update(extra or {})
    info = {"mode": "FLEX", "param": par, "enable_spin_orbital": False, "calc_scheme": "general"}
    return flex_mod.FLEX(r.get_param("ham"), {}, info), r


class _Files:
    """total.npz (legacy, unmarked), split.npz (a Phase B archive) and a UHFk
    trans_mod.npz on the same lattice, written once per class."""

    def __init__(self, d):
        self.d = d
        s, r = _flex(False)
        gi = r.get_param("green")
        s.solve(gi, d)
        s.save_results({"path_to_output": d, "sigma": "total.npz"}, gi)
        s2, r2 = _flex(True)
        gi2 = r2.get_param("green")
        s2.solve(gi2, d)
        s2.save_results({"path_to_output": d, "sigma": "split.npz"}, gi2)
        self.total = os.path.join(d, "total.npz")
        self.split = os.path.join(d, "split.npz")
        # UHFk mean field on the same input (paramagnetic, T = 2)
        from hwave.solver.uhfk import UHFk
        nvol, norb = s.lattice.nvol, s.norb
        par = {"T": 2.0, "Ncond": int(round(0.5 * 2 * norb * nvol)), "2Sz": 0, "flag_fock": True,
               "CellShape": [4, 4, 1], "SubShape": [1, 1, 1], "IterationMax": 200, "EPS": 8, "Mix": 0.5,
               "RndSeed": 1}
        u = UHFk(r.get_param("ham"), {"print_level": 1, "print_step": 1}, {"mode": "UHFk", "param": par})
        u.solve(r.get_param("green"), d)
        self.trans_mod = os.path.join(d, "trans_mod.npz")
        u._save_trans_mod(self.trans_mod)
        self.norb, self.nvol = norb, nvol


def _run(argv):
    from hwave.sigma_split import main
    return main(argv)


class TestSigmaSplitCLI(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        cls.f = _Files(cls._tmp.name)

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def _out(self, name):
        return os.path.join(self.f.d, name)

    def _is_valid_seed(self, path):
        from hwave.solver.flex_hf import make_seed_envelope, validate_split_seed
        data = np.load(path)
        env = make_seed_envelope(data, path, data["sigma"], None)
        self.assertEqual(env.marker, "split")
        st, fl = validate_split_seed(env, data["sigma"].shape)
        np.testing.assert_allclose(st + fl, data["sigma"], rtol=0, atol=1e-13)
        return data

    def test_argument_coupling(self):
        f = self.f
        for argv in ([f.total, self._out("a.npz")],
                     [f.total, self._out("a.npz"), "--zero-static", "--static", f.split],
                     [f.total, self._out("a.npz"), "--uhfk-trans-mod", f.trans_mod],
                     [f.total, self._out("a.npz"), "--zero-static", "--bare-transfer", "x"],
                     [f.total, self._out("a.npz"), "--zero-static", "--uhfk-spin-major"]):
            with self.subTest(argv=argv):
                self.assertEqual(_run(argv), 1)
                self.assertFalse(os.path.exists(self._out("a.npz")))

    def test_zero_static(self):
        out = self._out("zero.npz")
        self.assertEqual(_run([self.f.total, out, "--zero-static"]), 0)
        data = self._is_valid_seed(out)
        self.assertEqual(np.abs(data["sigma_static"]).max(), 0.0)
        np.testing.assert_array_equal(data["sigma_fluct"], data["sigma"])
        total = np.load(self.f.total)
        for k in total.files:
            self.assertIn(k, data.files)
        self.assertNotIn("scf_converged", data.files)          # no synthesised provenance
        # overwrite refused without --force
        self.assertEqual(_run([self.f.total, out, "--zero-static"]), 1)
        self.assertEqual(_run([self.f.total, out, "--zero-static", "--force"]), 0)
        # the seed is accepted by a Phase B run
        s, r = _flex(True)
        gi = r.get_param("green")
        gi.update(s.read_init({"path_to_input": self.f.d, "sigma_init": "zero.npz"}))
        s.solve(gi, self.f.d)

    def test_static_contract(self):
        f = self.f
        split = np.load(f.split)
        # rank 5 (from a Phase B archive) and the sigma fallback, rank-4 promotion
        st5 = os.path.join(f.d, "st5.npz"); st4 = os.path.join(f.d, "st4.npz"); fb = os.path.join(f.d, "fb.npz")
        meta = {k: split[k] for k in ("cell_shape", "momentum_convention", "wavevector_unit", "wavevector_index")}
        np.savez(st5, sigma_static=split["sigma_static"], **meta)
        np.savez(st4, sigma_static=split["sigma_static"][:, 0], **meta)
        np.savez(fb, sigma=split["sigma_static"], **meta)
        for src in (st5, st4, fb):
            out = self._out("from_static.npz")
            self.assertEqual(_run([f.total, out, "--static", src, "--force"]), 0)
            data = self._is_valid_seed(out)
            np.testing.assert_array_equal(data["sigma_static"], split["sigma_static"])
        # metadata mismatch / missing
        bad = os.path.join(f.d, "bad.npz")
        m2 = dict(meta); m2["cell_shape"] = np.array([2, 8, 1])
        np.savez(bad, sigma_static=split["sigma_static"], **m2)
        self.assertEqual(_run([f.total, self._out("bad_out.npz"), "--static", bad]), 1)
        m3 = dict(meta); del m3["wavevector_index"]
        np.savez(bad, sigma_static=split["sigma_static"], **m3)
        self.assertEqual(_run([f.total, self._out("bad_out.npz"), "--static", bad]), 1)
        # wrong shape
        np.savez(bad, sigma_static=split["sigma_static"][:, :, :8], **meta)
        self.assertEqual(_run([f.total, self._out("bad_out.npz"), "--static", bad]), 1)
        # non-Hermitian
        nh = np.array(split["sigma_static"]); nh[0, 0, :, 0, 1] += 0.1
        np.savez(bad, sigma_static=nh, **meta)
        self.assertEqual(_run([f.total, self._out("bad_out.npz"), "--static", bad]), 1)
        self.assertFalse(os.path.exists(self._out("bad_out.npz")))

    def test_spin_major_static(self):
        f = self.f
        split = np.load(f.split)
        meta = {k: split[k] for k in ("cell_shape", "momentum_convention", "wavevector_unit", "wavevector_index")}
        st = split["sigma_static"][0, 0]
        sm = np.einsum('kab,st->ksatb', st, np.eye(2)).reshape(f.nvol, 2 * f.norb, 2 * f.norb)
        path = os.path.join(f.d, "sm.npz")
        np.savez(path, sigma_static=sm, **meta)
        out = self._out("sm_out.npz")
        self.assertEqual(_run([f.total, out, "--static", path, "--uhfk-spin-major"]), 0)
        np.testing.assert_array_equal(np.load(out)["sigma_static"], split["sigma_static"])
        self.assertEqual(_run([f.total, out, "--static", path, "--force"]), 1)     # rank/shape refused without the flag
        sm[:, 0, f.norb] = 0.3                                                       # spin mixing
        np.savez(path, sigma_static=sm, **meta)
        self.assertEqual(_run([f.total, out, "--static", path, "--uhfk-spin-major", "--force"]), 1)

    def test_uhfk_trans_mod_recipe_equals_the_solver_seed(self):
        f = self.f
        out = self._out("uhfk.npz")
        self.assertEqual(_run([f.total, out, "--uhfk-trans-mod", f.trans_mod, "--bare-transfer",
                               os.path.join(_IN2, "transfer.dat")]), 0)
        data = self._is_valid_seed(out)
        # the in-solver D7 path: trans_mod through read_init -> Delta H
        s, r = _flex(True)
        gi = r.get_param("green")
        gi.update(s.read_init({"path_to_input": f.d, "trans_mod": "trans_mod.npz"}))
        s._phase_b_reset(gi)
        s._phase_b_preflight(gi)
        static, fluct = s._phase_b_seed
        np.testing.assert_allclose(data["sigma_static"], static, rtol=0, atol=1e-12)
        self.assertGreater(np.abs(static).max(), 1e-3)
        # the .npz form of the bare transfer
        tr = np.load(f.trans_mod)
        npz_transfer = os.path.join(f.d, "transfer.npz")
        np.savez(npz_transfer, Transfer=np.asarray(s.ham_info.ham_trans_r).reshape(f.nvol, f.norb, f.norb))
        out2 = self._out("uhfk2.npz")
        self.assertEqual(_run([f.total, out2, "--uhfk-trans-mod", f.trans_mod, "--bare-transfer", npz_transfer]), 0)
        np.testing.assert_allclose(np.load(out2)["sigma_static"], static, rtol=0, atol=1e-12)

    def test_marker_admissibility(self):
        f = self.f
        self.assertEqual(_run([f.split, self._out("x.npz"), "--zero-static"]), 1)   # split input refused
        total = np.load(f.total)
        marked = os.path.join(f.d, "marked.npz")
        np.savez(marked, sigma_convention=np.str_("total"), **{k: total[k] for k in total.files})
        self.assertEqual(_run([marked, self._out("x.npz"), "--zero-static"]), 0)
        np.savez(marked, sigma_convention=np.str_("weird"), **{k: total[k] for k in total.files})
        self.assertEqual(_run([marked, self._out("y.npz"), "--zero-static"]), 1)
        self.assertFalse(os.path.exists(self._out("y.npz")))


if __name__ == "__main__":
    unittest.main()
