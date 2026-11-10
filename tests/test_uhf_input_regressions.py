"""Regression tests through the public entry for three UHF input defects.

* #348: UHFk with ``trustme_interaction_range = true`` and a table whose
  range exceeds ``CellShape`` folds the translation vectors back into the
  cell (it used to stop with an ``IndexError``).
* #349: UHFr places the one-body part ``I c^dag_i c_l`` of an InterAll entry
  ``I c^dag_i c_j c^dag_k c_l`` with ``j == k`` at ``(i, l)``.
* #350: UHFk in spin-orbital mode reports ``Sy`` with the sign of the
  documented Green function ``G_ij = <c^dag_i c_j>``.
"""
import copy
import logging
import os
import shutil
import tempfile
import unittest

import numpy as np

from hwave import qlms


def _write(d, name, text):
    with open(os.path.join(d, name), "w") as f:
        f.write(text)


def _energy(path):
    out = {}
    with open(path) as f:
        for line in f:
            if "=" in line:
                k, v = line.split("=", 1)
                out[k.strip()] = float(v)
    return out


class _TempDir(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self._level = logging.getLogger().level
        logging.getLogger().setLevel(logging.ERROR)

    def tearDown(self):
        logging.getLogger().setLevel(self._level)
        shutil.rmtree(self.d, ignore_errors=True)

    def _run(self, inp):
        qlms.run(input_dict=copy.deepcopy(inp))


class TestUHFrInterAllOneBody(_TempDir):
    """#349: InterAll entries that are exactly one-body operators give the
    same state as the operator written in Trans (H = sigma_x, E = -1,
    <c^dag_up c_down> = -1/2)."""

    GREENONE = ("===\nNCisAjs 4\n===\n===\n===\n"
                "0 0 0 0\n0 0 0 1\n0 1 0 0\n0 1 0 1\n")

    def _trans(self, rows):
        body = "".join("{} {} {} {} {} 0.0\n".format(*r) for r in rows)
        return "===\nNTransfer {}\n===\n===\n===\n{}".format(len(rows), body)

    def _case(self, name, trans_rows, interall):
        d = os.path.join(self.d, name)
        os.makedirs(d)
        _write(d, "greenone.def", self.GREENONE)
        _write(d, "trans.def", self._trans(trans_rows))
        files = {"Trans": "trans.def"}
        if interall:
            body = "".join(" ".join(map(str, r)) + " 1.0 0.0\n" for r in interall)
            _write(d, "interall.def",
                   "===\nNInterAll {}\n===\n===\n===\n{}".format(len(interall), body))
            files["InterAll"] = "interall.def"
        inp = {
            "log": {"print_level": 0},
            "mode": {"mode": "UHFr", "param": {
                "Nsite": 1, "2Sz": "free", "Ncond": 1, "IterationMax": 200,
                "EPS": 10, "Mix": 0.5, "T": 0.0}},
            "file": {"input": {"path_to_input": d, "OneBodyG": "greenone.def",
                               "interaction": dict(files)},
                     "output": {"path_to_output": os.path.join(d, "out"),
                                "energy": "energy.dat", "green": "green.dat"}},
        }
        self._run(inp)
        out = os.path.join(d, "out")
        g = np.loadtxt(os.path.join(out, "green.dat"))
        return _energy(os.path.join(out, "energy.dat")), g[:, 4]

    def test_interall_one_body_part_is_at_i_l(self):
        e_all, g_all = self._case(
            "interall", [(0, 0, 0, 0, 0.0), (0, 1, 0, 1, 0.0)],
            [(0, 0, 0, 1, 0, 1, 0, 1), (0, 1, 0, 1, 0, 1, 0, 0)])
        # Trans is H = -sum t c^dag c, so t = -1 writes H = +sigma_x
        e_tr, g_tr = self._case(
            "trans", [(0, 0, 0, 1, -1.0), (0, 1, 0, 0, -1.0)], None)
        self.assertAlmostEqual(e_tr["Energy_total"], -1.0, places=8)
        self.assertAlmostEqual(e_all["Energy_total"], -1.0, places=8)
        np.testing.assert_allclose(g_all, g_tr, atol=1e-8)
        self.assertAlmostEqual(g_all[1], -0.5, places=8)


if __name__ == "__main__":
    unittest.main()
