"""Characterization tests for UHFk._find_dist_group_zero_t T=0 occupation scatter.

These tests pin down the exact occupation distribution produced by the T=0
occupation routine so that a vectorized rewrite of the per-element scatter can be
proven numerically identical.
"""
import unittest

import numpy as np

from hwave.solver.uhfk import UHFk


def _make_bare_uhfk(shape):
    """Build a UHFk instance with only the attributes _find_dist_group_zero_t needs."""
    solver = object.__new__(UHFk)
    nx, ny, nz = shape
    solver.shape = (nx, ny, nz)
    solver.nvol = nx * ny * nz
    solver.T = 0
    return solver


def _reference_zero_t(solver, ws_list, ncond):
    """Reference implementation: the original per-element scatter loop.

    Mirrors the pre-optimization _find_dist_group_zero_t body so we can assert
    the vectorized form is identical.
    """
    def _ksq_table(width):
        nx, ny, nz = solver.shape
        nvol = solver.nvol
        kx = np.roll((np.arange(nx) - nx // 2), -nx // 2) ** 2
        ky = np.roll((np.arange(ny) - ny // 2), -ny // 2) ** 2
        kz = np.roll((np.arange(nz) - nz // 2), -nz // 2) ** 2
        rr = np.zeros((nx, ny, nz))
        rr += np.broadcast_to(kx.reshape(nx, 1, 1), (nx, ny, nz))
        rr += np.broadcast_to(ky.reshape(1, ny, 1), (nx, ny, nz))
        rr += np.broadcast_to(kz.reshape(1, 1, nz), (nx, ny, nz))
        return np.broadcast_to(rr.reshape(nvol, 1), (nvol, width))

    all_ww = []
    all_ksq = []
    all_block_idx = []
    all_local_idx = []
    for b, w in enumerate(ws_list):
        k_sq = _ksq_table(w.shape[1]).flatten()
        ww = w.flatten()
        all_ww.append(ww)
        all_ksq.append(k_sq)
        all_block_idx.append(np.full(ww.size, b, dtype=int))
        all_local_idx.append(np.arange(ww.size))

    all_ww = np.concatenate(all_ww)
    all_ksq = np.concatenate(all_ksq)
    all_block_idx = np.concatenate(all_block_idx)
    all_local_idx = np.concatenate(all_local_idx)

    ev_idx = np.lexsort((all_ksq, all_ww))[0:ncond]

    dists = [np.zeros(w.size) for w in ws_list]
    for idx in ev_idx:
        b = all_block_idx[idx]
        li = all_local_idx[idx]
        dists[b][li] = 1.0
    dists = [d.reshape(w.shape) for d, w in zip(dists, ws_list)]
    return dists


class TestZeroTOccupation(unittest.TestCase):
    def _assert_matches_reference(self, shape, ws_list, ncond):
        solver = _make_bare_uhfk(shape)
        ref = _reference_zero_t(solver, ws_list, ncond)
        got, mu = solver._find_dist_group_zero_t(ws_list, None, ncond)

        self.assertEqual(len(got), len(ref))
        for g, r in zip(got, ref):
            self.assertEqual(g.shape, r.shape)
            self.assertEqual(g.dtype, r.dtype)
            np.testing.assert_array_equal(g, r)
        # occupation count must equal ncond
        total_occ = sum(g.sum() for g in got)
        self.assertEqual(total_occ, float(ncond))
        return got, mu

    def test_single_block_lowest_filled(self):
        # nvol=1 so each block is a single k-point; eigenvalues directly sortable
        rng = np.random.default_rng(0)
        w = rng.standard_normal((1, 6))
        ws_list = [w]
        got, mu = self._assert_matches_reference((1, 1, 1), ws_list, ncond=3)
        # explicit expected: the 3 lowest of the 6 eigenvalues are occupied
        order = np.argsort(w.flatten())
        expected = np.zeros(6)
        expected[order[:3]] = 1.0
        np.testing.assert_array_equal(got[0].flatten(), expected)

    def test_multi_block(self):
        rng = np.random.default_rng(1)
        ws_list = [
            rng.standard_normal((1, 4)),
            rng.standard_normal((1, 3)),
            rng.standard_normal((1, 5)),
        ]
        self._assert_matches_reference((1, 1, 1), ws_list, ncond=7)

    def test_multi_kpoint_multi_block(self):
        rng = np.random.default_rng(2)
        # nvol = 2*2*1 = 4 k-points
        ws_list = [
            rng.standard_normal((4, 3)),
            rng.standard_normal((4, 2)),
        ]
        self._assert_matches_reference((2, 2, 1), ws_list, ncond=11)

    def test_ncond_zero(self):
        rng = np.random.default_rng(3)
        ws_list = [rng.standard_normal((1, 4))]
        got, mu = self._assert_matches_reference((1, 1, 1), ws_list, ncond=0)
        self.assertEqual(got[0].sum(), 0.0)

    def test_ncond_full(self):
        rng = np.random.default_rng(4)
        ws_list = [rng.standard_normal((2, 3)), rng.standard_normal((2, 2))]
        # nvol = 2*1*1 = 2; total states = 2*3 + 2*2 = 10
        got, mu = self._assert_matches_reference((2, 1, 1), ws_list, ncond=10)
        for g in got:
            np.testing.assert_array_equal(g, np.ones_like(g))


if __name__ == "__main__":
    unittest.main()
