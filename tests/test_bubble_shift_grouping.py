"""One bond-bubble block per distinct channel shift (issue #224).

Every ``(m, m')`` block of the bond bubble depends on the pair only through
the shift ``Delta r_m - Delta r_m'``, so the bubble computes one block per
distinct shift and hands it to every pair of that shift. The reference here
is the per-pair loop itself: ``_pairs_by_shift`` patched to one group per
pair in row-major order (the previous iteration) must give BIT-identical
blocks on every entry point that groups -- the dynamic and static dense
longitudinal bubble, the IR static bubble and the transverse static bubble.
"""
import unittest
from unittest import mock

import numpy as np

from hwave.solver import bond_channels as bc
from hwave.solver import bubble
from tests.eliashberg_bond_fixtures import _green, _k_even_hopping

SHAPE, NORB, NMAT, BETA = (4, 3, 1), 2, 16, 2.0
#: first and second neighbours in x, first in y: 7 channels, 49 pairs
SHELLS = [(1, 0, 0), (-1, 0, 0), (2, 0, 0), (-2, 0, 0), (0, 1, 0), (0, -1, 0)]


def _ungrouped(delta_r_row, delta_r_col):
    out = []
    for m in range(len(delta_r_row)):
        for mp in range(len(delta_r_col)):
            shift = tuple(int(a) - int(b) for a, b in zip(delta_r_row[m], delta_r_col[mp]))
            out.append((shift, [(m, mp)]))
    return out


def _view():
    tbl = {"CoulombInter": {(R, (a, b)): 0.3 + 0j for R in SHELLS
                            for a in range(NORB) for b in range(NORB)}}
    topo = bc.resolve_bond_topology(tbl, np.eye(3), NORB,
                                    active_types=bc._LONGITUDINAL_ACTIVE_TYPES)
    return topo, bc.BondSetView(topo)


def _green_kw():
    return _green(_k_even_hopping(NORB, SHAPE, np.random.default_rng(3)), BETA, NMAT, mu=0.2)


def _n_shifts(delta_r):
    return len({tuple(int(x) for x in (a - b)) for a in np.asarray(delta_r)
                for b in np.asarray(delta_r)})


class TestPairsByShift(unittest.TestCase):

    def test_every_pair_once_grouped_by_shift(self):
        dr = np.array([(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0)])
        groups = bubble._pairs_by_shift(dr, dr)
        pairs = [p for _, ps in groups for p in ps]
        self.assertEqual(sorted(pairs), [(m, mp) for m in range(4) for mp in range(4)])
        for shift, ps in groups:
            self.assertEqual(ps, sorted(ps))
            for m, mp in ps:
                self.assertEqual(shift, tuple(int(x) for x in dr[m] - dr[mp]))
        self.assertEqual(len(groups), _n_shifts(dr))
        self.assertEqual(groups[0][0], (0, 0, 0))     # first appearance: (0, 0)


class TestLongitudinal(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.topo, cls.view = _view()
        cls.g = _green_kw()

    def _dynamic(self):
        return {k: np.array(b) for k, b in bubble._iter_bond_dynamic(
            self.g, None, BETA, self.view, spatial_shape=SHAPE)}

    def test_dynamic_bit_identical_and_one_block_per_shift(self):
        real = bubble._bond_pair_full_block
        with mock.patch.object(bubble, "_bond_pair_full_block", side_effect=real) as spy:
            got = self._dynamic()
        self.assertEqual(spy.call_count, _n_shifts(self.view.delta_r))
        self.assertLess(spy.call_count, self.view.n_channels ** 2)
        with mock.patch.object(bubble, "_pairs_by_shift", _ungrouped):
            want = self._dynamic()
        self.assertEqual(sorted(got), sorted(want))
        for k in want:
            np.testing.assert_array_equal(got[k], want[k], err_msg=str(k))

    def test_shared_block_is_read_only(self):
        it = bubble._iter_bond_dynamic(self.g, None, BETA, self.view, spatial_shape=SHAPE)
        _, block = next(it)
        self.assertFalse(block.flags.writeable)
        with self.assertRaises(ValueError):
            block[...] = 0.0

    def test_static_bit_identical(self):
        got = bubble.bond_bubble_static(self.g, None, BETA, self.view, spatial_shape=SHAPE)
        with mock.patch.object(bubble, "_pairs_by_shift", _ungrouped):
            want = bubble.bond_bubble_static(self.g, None, BETA, self.view,
                                             spatial_shape=SHAPE)
        np.testing.assert_array_equal(got, want)

    def test_ir_static_bit_identical(self):
        try:
            import sparse_ir  # noqa: F401
        except ImportError:
            self.skipTest("sparse_ir not installed")
        from hwave.solver.ir_axis import IRAxis
        axF = IRAxis(BETA, 10.0, 1e-8, "F")
        axB = IRAxis(BETA, 10.0, 1e-8, "B")
        h = _k_even_hopping(NORB, SHAPE, np.random.default_rng(4))
        iw = 1j * np.asarray(axF.freq_n) * np.pi / BETA
        g = np.stack([np.linalg.inv((w + 0.2) * np.eye(NORB)[None] - h) for w in iw])[None]
        got = bubble._ir_bond_static(g, axF, axB, self.view, spatial_shape=SHAPE)
        with mock.patch.object(bubble, "_pairs_by_shift", _ungrouped):
            want = bubble._ir_bond_static(g, axF, axB, self.view, spatial_shape=SHAPE)
        np.testing.assert_array_equal(got, want)
        self.assertGreater(float(np.abs(got).max()), 0.0)


class TestTransverse(unittest.TestCase):

    def test_transverse_static_bit_identical(self):
        from tests.test_bond_transverse_w import _transverse_topo_1d
        L = 5
        dr = [0, 1, -1, 2, -2]
        topo = _transverse_topo_1d(dr, [0, 2, 1, 4, 3], 1)
        rng = np.random.default_rng(5)
        g = np.stack([_green(_k_even_hopping(1, (L, 1, 1), rng), BETA, NMAT)[0]
                      for _ in range(2)])
        got = bubble.transverse_bond_bubble_static(g, None, BETA, topo, spatial_shape=(L, 1, 1))
        with mock.patch.object(bubble, "_pairs_by_shift", _ungrouped):
            want = bubble.transverse_bond_bubble_static(g, None, BETA, topo,
                                                        spatial_shape=(L, 1, 1))
        np.testing.assert_array_equal(got, want)


if __name__ == "__main__":
    unittest.main()
