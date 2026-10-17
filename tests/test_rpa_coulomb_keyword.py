"""Regression test: the combined 'Coulomb' interaction keyword in RPA/FLEX.

'Coulomb' is accepted by the shared k-space input reader and handled by UHFk
(split into the r=0 diagonal part = CoulombIntra and the rest = CoulombInter).
RPA's Interaction must apply the same decomposition instead of silently
ignoring the term (which yielded chiq without any Coulomb interaction).
"""
import unittest

import numpy as np

from hwave.solver.rpa import Interaction


class _FakeLattice:
    def __init__(self, shape):
        self.shape = shape
        self.nvol = shape[0] * shape[1] * shape[2]


def _make_bare_interaction(norb, shape, param_ham):
    obj = object.__new__(Interaction)
    obj.lattice = _FakeLattice(shape)
    obj.norb = norb
    obj.param_ham = param_ham
    obj._has_interaction = False
    obj._has_interaction_coulomb = False
    obj._has_interaction_exchange = False
    obj._has_interaction_pairhop = False
    return obj


class TestRPACoulombKeyword(unittest.TestCase):
    def test_coulomb_equals_intra_plus_inter(self):
        norb = 2
        shape = (2, 1, 1)
        coulomb = {
            ((0, 0, 0), (0, 0)): 2.0,   # r=0 diagonal -> CoulombIntra
            ((0, 0, 0), (1, 1)): 1.5,   # r=0 diagonal -> CoulombIntra
            ((0, 0, 0), (0, 1)): 0.3,   # r=0 off-diagonal -> CoulombInter
            ((0, 0, 0), (1, 0)): 0.3,
            ((1, 0, 0), (0, 0)): 0.5,   # r!=0 -> CoulombInter
            ((1, 0, 0), (1, 1)): 0.5,
        }
        intra = {k: v for k, v in coulomb.items()
                 if k[0] == (0, 0, 0) and k[1][0] == k[1][1]}
        inter = {k: v for k, v in coulomb.items() if k not in intra}

        obj_c = _make_bare_interaction(norb, shape, {"Coulomb": dict(coulomb)})
        obj_c._make_ham_inter()

        obj_ref = _make_bare_interaction(
            norb, shape, {"CoulombIntra": intra, "CoulombInter": inter})
        obj_ref._make_ham_inter()

        self.assertTrue(obj_c.has_interaction(),
                        "'Coulomb' must register as an interaction")
        self.assertTrue(obj_c.has_interaction_coulomb())
        self.assertGreater(np.max(np.abs(obj_c.ham_inter_q)), 1e-12)
        np.testing.assert_allclose(
            obj_c.ham_inter_q, obj_ref.ham_inter_q, rtol=0.0, atol=1e-12,
            err_msg="'Coulomb' must equal CoulombIntra + CoulombInter")
        # the Fierz correction (#104) must also agree between the two routes:
        # the fixture has an on-site inter-orbital entry, so it is nonempty
        self.assertIsNotNone(obj_c.ham_fierz_q)
        self.assertIsNotNone(obj_ref.ham_fierz_q)
        np.testing.assert_allclose(
            obj_c.ham_fierz_q, obj_ref.ham_fierz_q, rtol=0.0, atol=1e-12,
            err_msg="aggregate and explicit Fierz tensors must agree")

    def test_coulomb_conflicts_with_explicit_terms(self):
        norb = 1
        shape = (1, 1, 1)
        obj = _make_bare_interaction(
            norb, shape,
            {"Coulomb": {((0, 0, 0), (0, 0)): 2.0},
             "CoulombIntra": {((0, 0, 0), (0, 0)): 1.0}})
        with self.assertRaises(SystemExit):
            obj._make_ham_inter()


class TestHwaveScCoulombKeyword(unittest.TestCase):
    """hwave_sc must consume the combined 'Coulomb' keyword like the solvers
    (silently dropping it would build the Eliashberg vertex without the
    interaction the same TOML just used for chi0q/chiq)."""

    def _write_w90(self, filename, entries, norb):
        with open(filename, "w") as f:
            f.write("comment\n{}\n{}\n".format(norb, len(entries)))
            f.write(" ".join(["1"] * len(entries)) + "\n")
            for (r, a, b, v) in entries:
                f.write("  {:3d} {:3d} {:3d} {:3d} {:3d}  {:.6f}  0.0\n".format(
                    r[0], r[1], r[2], a + 1, b + 1, v))

    def test_coulomb_keyword_split_into_intra_and_inter(self):
        import os
        import shutil
        import tempfile
        from hwave.sc import _read_interaction_files

        d = tempfile.mkdtemp(prefix="sc_coulomb_")
        try:
            with open(os.path.join(d, "geom.dat"), "w") as f:
                f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n")
                f.write("2\n0.0 0.0 0.0\n0.0 0.0 0.0\n")
            self._write_w90(os.path.join(d, "transfer.dat"),
                            [((1, 0, 0), 0, 0, 1.0), ((-1, 0, 0), 0, 0, 1.0)],
                            norb=2)
            # Hermitian-closed declarations (issue #93: the readers now
            # reject one-sided/unclosed tables at read time)
            self._write_w90(os.path.join(d, "coulomb.dat"),
                            [((0, 0, 0), 0, 0, 2.0),   # r=0 diag -> intra
                             ((0, 0, 0), 0, 1, 0.3),   # r=0 offdiag -> inter
                             ((0, 0, 0), 1, 0, 0.3),
                             ((1, 0, 0), 0, 0, 0.5),   # r!=0 -> inter
                             ((-1, 0, 0), 0, 0, 0.5)],
                            norb=2)
            input_dict = {
                "file": {"input": {"interaction": {
                    "path_to_input": d,
                    "Geometry": "geom.dat",
                    "Transfer": "transfer.dat",
                    "Coulomb": "coulomb.dat",
                }}},
            }
            geom, hr, interactions = _read_interaction_files(input_dict)
            self.assertIn("CoulombIntra", interactions)
            self.assertIn("CoulombInter", interactions)
            self.assertEqual(
                interactions["CoulombIntra"][((0, 0, 0), (0, 0))], 2.0)
            self.assertEqual(
                interactions["CoulombInter"][((0, 0, 0), (0, 1))], 0.3)
            self.assertEqual(
                interactions["CoulombInter"][((1, 0, 0), (0, 0))], 0.5)
        finally:
            shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
