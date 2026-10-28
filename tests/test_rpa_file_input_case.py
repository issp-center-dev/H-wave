#!/usr/bin/env python3
"""[file.input] trans_mod / green_init are case-insensitive in hwave (#279).

Run through the public entry ``hwave.qlms.run``: a mixed-case spelling must
give the same chi0q as the lower-case one (and differ from the bare band),
and two spellings of the same key must be refused.
"""

import os
import shutil
import tempfile
import unittest

import numpy as np

import hwave.qlms as qlms
from tests.test_rpa_trans_mod import INPUT_DIR, _make_trans_mod_npz


def _config(outdir, tm_key, tm_file):
    inp = {
        "path_to_input": INPUT_DIR,
        "interaction": {
            "path_to_input": INPUT_DIR,
            "Geometry": "geom_2orb.dat",
            "Transfer": "transfer_nonso_2orb.dat",
            "CoulombIntra": "coulombintra_2orb.dat",
        },
    }
    if tm_key is not None:
        inp[tm_key] = tm_file
    return {
        "log": {"print_level": 0},
        "mode": {"mode": "RPA", "enable_spin_orbital": False,
                 "calc_scheme": "general",
                 "param": {"T": 2.0, "filling": 0.5, "CellShape": [8, 1, 1],
                           "SubShape": [1, 1, 1], "Nmat": 32}},
        "file": {"input": inp, "output": {"path_to_output": outdir,
                                  "chiq": "chiq.npz"}},
    }


class TestFileInputCase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.mkdtemp()
        self.tm = os.path.join(self._tmp, "tm.npz")
        _make_trans_mod_npz(self.tm)

    def tearDown(self):
        shutil.rmtree(self._tmp, ignore_errors=True)

    def _chiq(self, name, key):
        out = os.path.join(self._tmp, name)
        qlms.run(input_dict=_config(out, key, self.tm))
        return np.load(os.path.join(out, "chiq.npz"))["chiq"]

    def test_mixed_case_trans_mod_is_honoured(self):
        lower = self._chiq("lower", "trans_mod")
        mixed = self._chiq("mixed", "Trans_Mod")
        bare = self._chiq("bare", None)
        np.testing.assert_allclose(mixed, lower, atol=1e-12)
        self.assertGreater(np.max(np.abs(lower - bare)), 1e-6)

    def test_duplicate_spellings_refused(self):
        cfg = _config(os.path.join(self._tmp, "dup"), "trans_mod", self.tm)
        cfg["file"]["input"]["Trans_Mod"] = self.tm
        with self.assertRaises(ValueError):
            qlms.run(input_dict=cfg)


if __name__ == "__main__":
    unittest.main()
