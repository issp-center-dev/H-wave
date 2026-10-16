"""Thread-parallel batched linear algebra of the bond dressing: every helper
must reproduce the serial numpy result bit for bit, and dress_batch /
dress_and_build_w must not move under ``workers``."""
import unittest

import numpy as np

from hwave.solver import bond_channels as bc
from hwave.solver import threaded as tl
from tests.test_flex_bond_dressing import _problem

needs_pool = unittest.skipUnless(tl.available(), "threadpoolctl is not installed")


def _stack(n, nd, seed=0):
    rng = np.random.default_rng(seed)
    a = rng.standard_normal((n, nd, nd)) + 1j * rng.standard_normal((n, nd, nd))
    return a + 3.0 * np.eye(nd)


class TestHelpersEqualSerial(unittest.TestCase):

    def test_solve_is_exact_for_every_worker_count(self):
        a, b = _stack(700, 6, 1), _stack(700, 6, 2)
        ref = np.linalg.solve(a, b)
        for w in (None, 1, 2, 3, 8):
            with self.subTest(workers=w):
                np.testing.assert_array_equal(tl.solve(a, b, workers=w, min_rows=16), ref)

    def test_svdvals_is_exact(self):
        a = _stack(333, 5, 3)
        ref = np.linalg.svd(a, compute_uv=False)
        for w in (1, 4, 7):
            with self.subTest(workers=w):
                np.testing.assert_array_equal(tl.svdvals(a, workers=w, min_rows=8), ref)

    def test_matmul_broadcast_on_either_batch_axis(self):
        rng = np.random.default_rng(4)
        x = rng.standard_normal((5, 9, 4, 4)) + 1j * rng.standard_normal((5, 9, 4, 4))
        w = rng.standard_normal((1, 9, 4, 4)) + 0j
        ref = x @ w
        for workers in (1, 2, 4, 16):
            with self.subTest(workers=workers):
                np.testing.assert_array_equal(tl.matmul(x, w, workers=workers, min_rows=1), ref)
        # a batch axis of 1 on the left operand and a full one on the right
        np.testing.assert_array_equal(tl.matmul(w, x, workers=3, min_rows=1), w @ x)

    def test_sandwich_equals_the_left_associated_triple_product(self):
        rng = np.random.default_rng(5)
        x = rng.standard_normal((6, 7, 3, 3)) + 1j * rng.standard_normal((6, 7, 3, 3))
        s = rng.standard_normal((7, 3, 3)) + 0j
        ref = s[None] @ x @ s[None]
        for workers in (1, 2, 5):
            with self.subTest(workers=workers):
                np.testing.assert_array_equal(
                    tl.sandwich(s, x, s, workers=workers, min_rows=1), ref)

    def test_a_singular_block_raises_linalgerror_from_a_worker(self):
        a = _stack(200, 4, 6)
        a[137] = 0.0
        with self.assertRaises(np.linalg.LinAlgError):
            tl.solve(a, a, workers=4, min_rows=8)

    def test_invalid_worker_count_is_refused(self):
        a = _stack(8, 3)
        with self.assertRaises(ValueError):
            tl.solve(a, a, workers=-2)

    def test_small_stacks_stay_in_one_chunk(self):
        self.assertEqual(len(tl._row_slices(10, 8, 256)), 1)
        self.assertEqual(len(tl._row_slices(1024, 8, 256)), 4)
        sl = tl._row_slices(1000, 3, 1)
        self.assertEqual(sl[0].start, 0)
        self.assertEqual(sl[-1].stop, 1000)
        self.assertEqual([s.stop for s in sl[:-1]], [s.start for s in sl[1:]])


@needs_pool
class TestDressBatchWorkers(unittest.TestCase):

    def setUp(self):
        self._min_rows = tl._MIN_ROWS
        tl._MIN_ROWS = 1               # the fixtures are tiny: force real chunking

    def tearDown(self):
        tl._MIN_ROWS = self._min_rows

    def _run(self, cb, W, channel, workers, **kw):
        v = []
        try:
            chi, cond = bc.dress_batch(cb, W, channel, workers=workers, violations=v, **kw)
            return ("ok", chi, cond, tuple(sorted((k, r[k]) for r in v for k in sorted(r))))
        except bc.BondConditioningError as e:
            return ("refused", str(e), e.channel, e.l, e.iq, e.q, e.worst, e.ratio, e.pole)

    def _same(self, cb, W, channel, **kw):
        a = self._run(cb, W, channel, 1, **kw)
        b = self._run(cb, W, channel, 4, **kw)
        self.assertEqual(a[0], b[0])
        if a[0] == "ok":
            np.testing.assert_array_equal(a[1], b[1])
            self.assertEqual(a[2:], b[2:])
        else:
            self.assertEqual(a[1:], b[1:])

    def test_dress_batch_is_independent_of_workers(self):
        chi_bar, S, C = _problem(nmat=6, nvol=8)
        nmat, nvol = chi_bar.shape[:2]
        for channel, W in (("spin", S), ("charge", C)):
            for guard_freqs in ("all", "static"):
                for method in ("svd", "interval"):
                    with self.subTest(channel=channel, guard_freqs=guard_freqs, method=method):
                        self._same(chi_bar, W, channel, l0=0, nmat=nmat,
                                   spatial_shape=(nvol, 1, 1), guard_freqs=guard_freqs,
                                   guard_method=method)

    def test_refusals_are_independent_of_workers(self):
        chi_bar, S, _C = _problem(nmat=6, nvol=8)
        nmat, nvol = chi_bar.shape[:2]
        for scale in (4.0, 16.0):
            for policy in ("refuse", "warn"):
                with self.subTest(scale=scale, policy=policy):
                    self._same(chi_bar, scale * S, "spin", l0=0, nmat=nmat,
                               spatial_shape=(nvol, 1, 1), guard_policy=policy)


@needs_pool
class TestDressAndBuildWWorkers(unittest.TestCase):

    def setUp(self):
        self._min_rows = tl._MIN_ROWS
        tl._MIN_ROWS = 1

    def tearDown(self):
        tl._MIN_ROWS = self._min_rows

    def test_effective_interaction_is_independent_of_workers(self):
        import types
        from hwave.solver.flex_bond import BondBlockStore, BondDeviceContext, dress_and_build_w
        from tests.test_flex_bond_dressing import _onsite_factors
        chi_bar, S, C = _problem(nmat=6, nvol=4, nd=4, B=3)
        nmat, nvol, ND = chi_bar.shape[:3]
        nd = 4
        out = {}
        for workers in (1, 3):
            with BondBlockStore(nmat, nvol, ND, nd, ("chibar", "W", "chi_s_w", "chi_c_w")) as store:
                store.put_freq_batch("chibar", 0, nmat, chi_bar)
                S_on = np.ascontiguousarray(S[:1, :nd, :nd]).repeat(nvol, axis=0) * 0.7
                C_on = np.ascontiguousarray(C[:1, :nd, :nd]).repeat(nvol, axis=0) * 0.3
                view = types.SimpleNamespace(n_channels=ND // nd)
                with BondDeviceContext.for_view(np, S, C, S_on, C_on, view, 2) as dev:
                    res = dress_and_build_w(store, dev, nb=4, output_full=True, nmat=nmat,
                                            nvol=nvol, nd=nd, spatial_shape=(4, 1, 1),
                                            factors=_onsite_factors(), workers=workers)
                out[workers] = (store.get_freq_batch("W", 0, nmat).copy(), res.collapse_s,
                                res.collapse_c, res.static_s, res.static_c, res.cond_min_s)
        for a, b in zip(out[1], out[3]):
            np.testing.assert_array_equal(a, b)


if __name__ == "__main__":
    unittest.main()
