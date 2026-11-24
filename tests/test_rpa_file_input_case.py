#!/usr/bin/env python3
"""[file.input] keys are case-insensitive in hwave (#279).

Run through the public entry ``hwave.qlms.run``: a mixed-case spelling of
trans_mod / green_init / chi0q_init (RPA) or sigma_init (FLEX) must give the
same result as the lower-case one (and differ from the run without it), and
two spellings of the same key must be refused.
"""

import os
import shutil
import tempfile
import unittest

import numpy as np

import hwave.qlms as qlms
from tests.test_rpa_trans_mod import INPUT_DIR, _make_trans_mod_npz


def _config(outdir, extra, mode="RPA", param=None, output=None):
    inp = {
        "path_to_input": INPUT_DIR,
        "interaction": {
            "path_to_input": INPUT_DIR,
            "Geometry": "geom_2orb.dat",
            "Transfer": "transfer_nonso_2orb.dat",
            "CoulombIntra": "coulombintra_2orb.dat",
        },
    }
    inp.update(extra)
    p = {"T": 2.0, "filling": 0.5, "CellShape": [8, 1, 1],
         "SubShape": [1, 1, 1], "Nmat": 32}
    p.update(param or {})
    out = {"path_to_output": outdir}
    out.update(output or {"chiq": "chiq.npz"})
    return {
        "log": {"print_level": 0},
        "mode": {"mode": mode, "enable_spin_orbital": False,
                 "calc_scheme": "general", "param": p},
        "file": {"input": inp, "output": out},
    }


class TestFileInputCase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.mkdtemp()
        self.tm = os.path.join(self._tmp, "tm.npz")
        _make_trans_mod_npz(self.tm)
        # spin-polarized on-site green: shifts the Hartree band of the run
        self.gi = os.path.join(self._tmp, "gi.npz")
        green = np.zeros((8, 4, 4), dtype=np.complex128)
        green[0] = np.diag([0.9, 0.9, 0.1, 0.1])
        np.savez(self.gi, green=green)

    def tearDown(self):
        shutil.rmtree(self._tmp, ignore_errors=True)

    def _chiq(self, name, extra):
        out = os.path.join(self._tmp, name)
        qlms.run(input_dict=_config(out, extra))
        return np.load(os.path.join(out, "chiq.npz"))["chiq"]

    def _check_key(self, lower_key, mixed_key, value):
        lower = self._chiq("lower", {lower_key: value})
        mixed = self._chiq("mixed", {mixed_key: value})
        bare = self._chiq("bare", {})
        np.testing.assert_allclose(mixed, lower, rtol=0, atol=1e-12)
        self.assertGreater(np.max(np.abs(lower - bare)), 1e-6)

    def test_mixed_case_trans_mod_is_honoured(self):
        self._check_key("trans_mod", "Trans_Mod", self.tm)

    def test_mixed_case_green_init_is_honoured(self):
        self._check_key("green_init", "Green_Init", self.gi)

    def test_mixed_case_chi0q_init_is_honoured(self):
        # chi0q of a run with trans_mod, reloaded into a run without it
        out = os.path.join(self._tmp, "src")
        qlms.run(input_dict=_config(out, {"trans_mod": self.tm},
                                    output={"chi0q": "chi0q.npz"}))
        c0 = os.path.join(out, "chi0q.npz")
        self._check_key("chi0q_init", "Chi0q_Init", c0)

    def test_duplicate_spellings_refused(self):
        for lower, mixed in (("trans_mod", "Trans_Mod"),
                             ("green_init", "Green_Init"),
                             ("chi0q_init", "Chi0q_Init")):
            with self.subTest(key=lower):
                cfg = _config(os.path.join(self._tmp, "dup"),
                              {lower: self.tm, mixed: self.tm})
                with self.assertRaisesRegex(ValueError, "more than once"):
                    qlms.run(input_dict=cfg)


class TestFlexFileInputCase(unittest.TestCase):
    PARAM = {"IterationMax": 300, "Mix": 0.5, "EPS": 10}

    def setUp(self):
        self._tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._tmp, ignore_errors=True)

    def _flex(self, name, extra, param):
        out = os.path.join(self._tmp, name)
        p = dict(self.PARAM)
        p.update(param)
        cfg = _config(out, extra, mode="FLEX", param=p,
                      output={"sigma": "sigma.npz"})
        qlms.run(input_dict=cfg)
        return os.path.join(out, "sigma.npz")

    def test_mixed_case_sigma_init_is_honoured(self):
        seed = self._flex("seed", {}, {})
        one = {"IterationMax": 1, "Mix": 1.0, "EPS": 1}
        sig = lambda f: np.load(f)["sigma"]
        lower = sig(self._flex("lower", {"sigma_init": seed}, one))
        mixed = sig(self._flex("mixed", {"Sigma_Init": seed}, one))
        cold = sig(self._flex("cold", {}, one))
        np.testing.assert_allclose(mixed, lower, rtol=0, atol=1e-12)
        self.assertGreater(np.max(np.abs(lower - cold)), 1e-6)

    def test_duplicate_sigma_init_refused(self):
        cfg = _config(os.path.join(self._tmp, "dup"),
                      {"sigma_init": "a.npz", "Sigma_Init": "a.npz"},
                      mode="FLEX", param=self.PARAM,
                      output={"sigma": "sigma.npz"})
        with self.assertRaisesRegex(ValueError, "more than once"):
            qlms.run(input_dict=cfg)


if __name__ == "__main__":
    unittest.main()
