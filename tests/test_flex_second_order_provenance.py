# tests/test_flex_second_order_provenance.py
"""Provenance members flex_second_order / flex_second_order_schema in every
general-scheme archive (spec section 3), absent on reduced archives, the
seed reader's WARNING on a kernel mismatch and INFO when absent, and the
off-site warning wording per kernel."""
import os
import tempfile
import unittest

import numpy as np

from tests.test_flex_second_order_config import _build


def _solve_save(s, r, out, files):
    gi = r.get_param("green")
    s.solve(gi, out)
    s.save_results(dict({"path_to_output": out}, **files), gi)
    return gi


class TestProvenance(unittest.TestCase):

    def test_members_in_every_general_archive_and_absent_on_reduced(self):
        with tempfile.TemporaryDirectory() as out:
            for so in ("local", "takimoto"):
                s, r = _build({"flex_second_order": so, "Nmat": 8})
                _solve_save(s, r, out, {"chi0q": "chi0q", "chiq": "chiq", "sigma": "sigma", "green": "green"})
                for f in ("chi0q.npz", "chiq.npz", "chiq_s.npz", "chiq_c.npz", "sigma.npz", "green.npz"):
                    z = np.load(os.path.join(out, f))
                    self.assertEqual(str(z["flex_second_order"]), so, f)
                    self.assertEqual(int(z["flex_second_order_schema"]), 1, f)
                    self.assertEqual(z["flex_second_order"].dtype.kind, "U")
            s, r = _build({"Nmat": 8}, calc_scheme="reduced")
            _solve_save(s, r, out, {"sigma": "sigma", "green": "green"})
            for f in ("sigma.npz", "green.npz", "chiq_s.npz"):
                self.assertNotIn("flex_second_order", np.load(os.path.join(out, f)).files, f)

    def test_seed_warning_and_info(self):
        with tempfile.TemporaryDirectory() as out:
            s, r = _build({"flex_second_order": "takimoto", "Nmat": 8})
            _solve_save(s, r, out, {"sigma": "sigma"})
            s2, _ = _build({"flex_second_order": "local", "Nmat": 8})
            with self.assertLogs("hwave.solver.flex", level="WARNING") as cm:
                s2.read_init({"path_to_input": out, "sigma_init": "sigma.npz"})
            self.assertTrue(any("seed computed with flex_second_order = takimoto; this run uses local" in m
                                for m in cm.output))
            z = dict(np.load(os.path.join(out, "sigma.npz")))
            z.pop("flex_second_order"); z.pop("flex_second_order_schema")
            np.savez(os.path.join(out, "legacy.npz"), **z)
            with self.assertLogs("hwave.solver.flex", level="INFO") as cm:
                s2.read_init({"path_to_input": out, "sigma_init": "legacy.npz"})
            self.assertTrue(any("not recorded" in m and "reduced scheme or pre-2.1" in m for m in cm.output))
            self.assertFalse(any("WARNING" in m and "seed computed" in m for m in cm.output))

    def test_offsite_warning_wording(self):
        for so, needle in (("local", "exact at second order"), ("takimoto", "is omitted (the same approximation")):
            with self.assertLogs("hwave.solver.flex", level="WARNING") as cm:
                s, r = _build({"flex_second_order": so, "Nmat": 8})
                s._calc_epsilon_k({})
                G = s._calc_dressed_green(0.5, 0.1, np.zeros((1, 8, 16, 2, 2), complex))
                s._flex_compute_veff_general(s._calc_chi0q(G, np.zeros_like(G), 0.5)[0], s.ham_info.ham_inter_q)
            self.assertTrue(any(needle in m for m in cm.output), (so, cm.output))


if __name__ == "__main__":
    unittest.main()
