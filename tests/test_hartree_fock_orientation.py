"""The mean-field tables read an off-site row (r, a, b, v) in the documented
orientation (spec 2026-09-16 D-1): after the reader's closure every off-site
displacement table is conjugate-transposed, which equals the closed table at
-r. On-site tables (r = 0) are untouched; the INFO line of D-8 is emitted
exactly when a table changed."""
import logging
import unittest

import numpy as np

from hwave.solver import hartree_fock as hf


def _closed_only(param_ham, norb, shape):
    """The tables BEFORE the orientation step: rebuilt here from the same
    per-type recipe (closure only) so the identity new[r] == closed[-r] is
    asserted against an independent construction, not against the function
    under test."""
    nx, ny, nz = shape
    out = {}
    sign = {"Hund": -1.0, "Exchange": -1.0}
    for t in ("CoulombInter", "Hund", "Ising", "PairLift", "Exchange", "PairHop"):
        if t not in param_ham:
            continue
        tab = np.zeros((nx, ny, nz, norb, norb), dtype=np.complex128)
        for (irvec, orbvec), v in param_ham[t].items():
            tab[(*irvec, *orbvec)] += v
        out[t] = sign.get(t, 1.0) * hf._reverse_closed(tab)
    return out


def _reversed_r(tab):
    """closed[-r]: the grid reversal of every displacement axis."""
    return np.roll(np.flip(tab, axis=(0, 1, 2)), shift=1, axis=(0, 1, 2))


class TestOrientation(unittest.TestCase):
    shape = (3, 3, 1)          # every axis odd: no displacement is its own reverse
    norb = 2

    def _ham(self):
        return {
            "CoulombIntra": {((0, 0, 0), (0, 0)): 2.0, ((0, 0, 0), (1, 1)): 1.5},
            "CoulombInter": {((1, 0, 0), (0, 1)): 0.4, ((-1, 0, 0), (1, 0)): 0.4,
                             ((0, 1, 0), (0, 0)): 0.2, ((0, -1, 0), (0, 0)): 0.2,
                             ((0, 0, 0), (0, 1)): 0.7, ((0, 0, 0), (1, 0)): 0.7},
            "Hund": {((1, 0, 0), (0, 1)): 0.1, ((-1, 0, 0), (1, 0)): 0.1},
            "Ising": {((1, 0, 0), (0, 1)): 0.15, ((-1, 0, 0), (1, 0)): 0.15},
            "PairLift": {((1, 0, 0), (0, 1)): 0.05, ((-1, 0, 0), (1, 0)): 0.05},
            "Exchange": {((0, 1, 0), (0, 1)): 0.12, ((0, -1, 0), (1, 0)): 0.12},
            "PairHop": {((1, 0, 0), (0, 1)): 0.08 + 0.03j, ((-1, 0, 0), (1, 0)): 0.08 - 0.03j,
                        ((0, 0, 0), (0, 1)): 0.2 + 0.1j, ((0, 0, 0), (1, 0)): 0.2 - 0.1j},
        }

    def test_offsite_tables_equal_the_closed_table_at_minus_r(self):
        ham = self._ham()
        tabs = hf.build_interaction_tables(ham, self.norb, self.shape)
        closed = _closed_only(ham, self.norb, self.shape)
        for t, ref in closed.items():
            with self.subTest(type=t):
                got = tabs.inter_table[t]
                want = _reversed_r(ref)
                want[0, 0, 0] = ref[0, 0, 0]            # r = 0 untouched
                np.testing.assert_array_equal(got, want)
                # non-vacuity: the orientation step changed something off-site
                self.assertFalse(np.array_equal(got, ref), t)

    def test_onsite_tables_are_untouched_including_complex_pairhop(self):
        ham = self._ham()
        tabs = hf.build_interaction_tables(ham, self.norb, self.shape)
        closed = _closed_only(ham, self.norb, self.shape)
        for t, ref in closed.items():
            np.testing.assert_array_equal(tabs.inter_table[t][0, 0, 0], ref[0, 0, 0], t)
        np.testing.assert_array_equal(tabs.inter_table["CoulombIntra"][0, 0, 0],
                                      np.diag([2.0, 1.5]).astype(complex))

    def test_orient_documented_is_an_involution_on_closed_tables(self):
        ham = self._ham()
        closed = _closed_only(ham, self.norb, self.shape)
        for t, ref in closed.items():
            twice = hf.orient_documented(hf.orient_documented(ref))
            np.testing.assert_array_equal(twice, ref, t)

    def test_real_single_orbital_and_real_diagonal_rows_are_invariant(self):
        ham = {"CoulombInter": {((1, 0, 0), (0, 0)): 0.3, ((-1, 0, 0), (0, 0)): 0.3}}
        tabs = hf.build_interaction_tables(ham, 1, self.shape)
        closed = _closed_only(ham, 1, self.shape)
        np.testing.assert_array_equal(tabs.inter_table["CoulombInter"], closed["CoulombInter"])

    def test_complex_single_orbital_pairhop_is_not_invariant(self):
        ham = {"PairHop": {((1, 0, 0), (0, 0)): 0.1 + 0.05j, ((-1, 0, 0), (0, 0)): 0.1 - 0.05j}}
        tabs = hf.build_interaction_tables(ham, 1, self.shape)
        closed = _closed_only(ham, 1, self.shape)
        self.assertFalse(np.array_equal(tabs.inter_table["PairHop"], closed["PairHop"]))
        np.testing.assert_array_equal(tabs.inter_table["PairHop"], _reversed_r(closed["PairHop"]))

    def test_info_line_only_when_a_table_changed(self):
        onsite = {"CoulombIntra": {((0, 0, 0), (0, 0)): 2.0},
                  "CoulombInter": {((0, 0, 0), (0, 1)): 0.7, ((0, 0, 0), (1, 0)): 0.7}}
        logger = logging.getLogger("qlms.solver.hartree_fock")
        records = []
        handler = logging.Handler()
        handler.setLevel(logging.INFO)
        handler.emit = records.append
        logger.addHandler(handler)
        old_level = logger.level
        logger.setLevel(logging.INFO)
        try:
            hf.build_interaction_tables(onsite, 2, self.shape)
            self.assertFalse([r for r in records if "documented orientation" in r.getMessage()])
            records.clear()
            hf.build_interaction_tables(self._ham(), 2, self.shape)
            hits = [r for r in records if "documented orientation" in r.getMessage()]
            self.assertEqual(len(hits), 1)
            # (type, r) entries changed: CoulombInter +-x (2), Hund 2, Ising 2, PairLift 2,
            # Exchange +-y (2), PairHop +-x (2); the +-y CoulombInter rows are real and
            # orbital-diagonal -> unchanged
            self.assertIn("12 off-site displacement table(s)", hits[0].getMessage())
            records.clear()
            diag = {"PairHop": {((1, 0, 0), (0, 0)): 0.1 + 0.05j, ((-1, 0, 0), (0, 0)): 0.1 - 0.05j}}
            hf.build_interaction_tables(diag, 1, self.shape)
            hits = [r for r in records if "documented orientation" in r.getMessage()]
            self.assertEqual(len(hits), 1)
            self.assertIn("2 off-site displacement table(s)", hits[0].getMessage())
        finally:
            logger.removeHandler(handler)
            logger.setLevel(old_level)


if __name__ == "__main__":
    unittest.main()
