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


_CHAIN_GEOM = "1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n1\n0.0 0.0 0.0\n"


def _wan(rows):
    lines = ["table", "1", str(len(rows)), " ".join(["1"] * len(rows))]
    lines += ["{} 0 0 1 1 {} 0.0".format(r, v) for r, v in rows]
    return "\n".join(lines) + "\n"


class TestUHFkTrustmeFolding(_TempDir):
    """#348: a hopping or interaction vector outside the cell is folded
    modulo CellShape at accumulation; for an off-site table the result
    equals the pre-folded table's exactly. The type of each row is still
    judged on the vector as written (an aggregate Coulomb row at R = L
    stays an inter-site row), so folding never reclassifies a row."""

    INTER_TYPES = ("CoulombInter", "Hund", "Ising", "Exchange", "PairLift", "PairHop")

    def _case(self, name, transfer, inter, subshape, inter_type="CoulombInter",
              trustme=True):
        d = os.path.join(self.d, name)
        os.makedirs(d)
        _write(d, "geom.dat", _CHAIN_GEOM)
        _write(d, "transfer.dat", _wan(transfer))
        _write(d, "inter.dat", _wan(inter))
        inp = {
            "log": {"print_level": 0},
            "mode": {"mode": "UHFk", "param": {
                "CellShape": [4, 1, 1], "SubShape": subshape,
                "trustme_interaction_range": trustme, "T": 0.0, "Ncond": 4,
                "IterationMax": 50}},
            "file": {"input": {"path_to_input": d, "interaction": {
                "path_to_input": d, "Geometry": "geom.dat",
                "Transfer": "transfer.dat", inter_type: "inter.dat"}},
                "output": {"path_to_output": os.path.join(d, "out"),
                           "energy": "energy", "green": "green"}},
        }
        self._run(inp)
        out = os.path.join(d, "out")
        return (_energy(os.path.join(out, "energy")),
                np.load(os.path.join(out, "green.npz"))["green"])

    def test_folded_equals_prefolded(self):
        far_t = [(1, -1.0), (-1, -1.0), (5, -0.1), (-5, -0.1)]
        near_t = [(1, -1.1), (-1, -1.1), (0, 0.0)]
        far_v = [(5, 0.3), (-5, 0.3)]
        near_v = [(1, 0.3), (-1, 0.3)]
        for sub in ([1, 1, 1], [2, 1, 1]):
            with self.subTest(subshape=sub):
                tag = "".join(map(str, sub))
                e_far, g_far = self._case("far" + tag, far_t, far_v, sub)
                e_near, g_near = self._case("near" + tag, near_t, near_v, sub)
                self.assertEqual(e_far, e_near)
                np.testing.assert_array_equal(g_far, g_near)

    def test_every_offsite_type(self):
        t = [(1, -1.0), (-1, -1.0)]
        for typ in self.INTER_TYPES:
            with self.subTest(type=typ):
                e_far, g_far = self._case("f" + typ, t, [(5, 0.3), (-5, 0.3)],
                                          [1, 1, 1], typ)
                e_near, g_near = self._case("n" + typ, t, [(1, 0.3), (-1, 0.3)],
                                            [1, 1, 1], typ)
                self.assertEqual(e_far, e_near)
                np.testing.assert_array_equal(g_far, g_near)

    def test_out_of_range_still_refused_without_trustme(self):
        with self.assertRaises(SystemExit):
            self._case("refused", [(1, -1.0), (-1, -1.0), (5, -0.1), (-5, -0.1)],
                       [(1, 0.3), (-1, 0.3)], [1, 1, 1], trustme=False)


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


class TestUHFkSOSpinSign(_TempDir):
    """#350: one spin-orbital site in a field H = -sigma_alpha has
    <S_alpha> = +1/2 for alpha = x, y, z (Transfer entry (a, b) is the
    coefficient of c^dag_a c_b, index 2*orb + spin + 1)."""

    FIELDS = {
        "Sx": [(1, 2, -1.0, 0.0), (2, 1, -1.0, 0.0)],
        "Sy": [(1, 2, 0.0, 1.0), (2, 1, 0.0, -1.0)],
        "Sz": [(1, 1, -1.0, 0.0), (2, 2, 1.0, 0.0)],
    }

    def _case(self, name, rows):
        d = os.path.join(self.d, name)
        os.makedirs(d)
        _write(d, "geom.dat", "1 0 0\n0 1 0\n0 0 1\n2\n0 0 0\n0 0 0\n")
        body = "".join("0 0 0 {} {} {} {}\n".format(*r) for r in rows)
        _write(d, "transfer.dat", "field\n2\n1\n 1\n" + body)
        inp = {
            "log": {"print_level": 0},
            "mode": {"mode": "UHFk", "enable_spin_orbital": True, "param": {
                "CellShape": [1, 1, 1], "SubShape": [1, 1, 1], "T": 0.0,
                "Ncond": 1, "IterationMax": 3}},
            "file": {"input": {"path_to_input": d, "interaction": {
                "path_to_input": d, "Geometry": "geom.dat",
                "Transfer": "transfer.dat"}},
                "output": {"path_to_output": os.path.join(d, "out"),
                           "energy": "energy", "green": "green"}},
        }
        self._run(inp)
        return _energy(os.path.join(d, "out", "energy"))

    def test_field_along_each_axis(self):
        for comp, rows in self.FIELDS.items():
            with self.subTest(component=comp):
                e = self._case(comp, rows)
                self.assertAlmostEqual(e["Energy_Total"], -1.0, places=10)
                for other in ("Sx", "Sy", "Sz"):
                    want = 0.5 if other == comp else 0.0
                    self.assertAlmostEqual(e[other], want, places=10, msg=other)


if __name__ == "__main__":
    unittest.main()
