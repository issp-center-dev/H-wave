"""Thread-parallel batched linear algebra of the bond dressing: every helper
must reproduce the serial numpy result to rounding, and dress_batch /
dress_and_build_w must not move under ``workers``. Not bit for bit: the
serial reference may run on multithreaded BLAS while the workers pin BLAS to
one thread, so the summation order (last digit) can differ."""
import unittest

import numpy as np

from hwave.solver import bond_channels as bc
from hwave.solver import threaded as tl
from tests.test_flex_bond_dressing import _problem

needs_pool = unittest.skipUnless(tl.available(), "threadpoolctl is not installed")


RTOL, ATOL = 1e-10, 1e-13


def _approx_equal(tc, a, b):
    if isinstance(a, (tuple, list)):
        tc.assertEqual(len(a), len(b))
        for x, y in zip(a, b):
            _approx_equal(tc, x, y)
    elif isinstance(a, (float, np.floating, np.ndarray, complex)):
        np.testing.assert_allclose(a, b, rtol=RTOL, atol=ATOL)
    else:
        tc.assertEqual(a, b)


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
                np.testing.assert_allclose(tl.solve(a, b, workers=w, min_rows=16), ref, rtol=RTOL, atol=ATOL)

    def test_svdvals_is_exact(self):
        a = _stack(333, 5, 3)
        ref = np.linalg.svd(a, compute_uv=False)
        for w in (1, 4, 7):
            with self.subTest(workers=w):
                np.testing.assert_allclose(tl.svdvals(a, workers=w, min_rows=8), ref, rtol=RTOL, atol=ATOL)

    def test_matmul_broadcast_on_either_batch_axis(self):
        rng = np.random.default_rng(4)
        x = rng.standard_normal((5, 9, 4, 4)) + 1j * rng.standard_normal((5, 9, 4, 4))
        w = rng.standard_normal((1, 9, 4, 4)) + 0j
        ref = x @ w
        for workers in (1, 2, 4, 16):
            with self.subTest(workers=workers):
                np.testing.assert_allclose(tl.matmul(x, w, workers=workers, min_rows=1), ref, rtol=RTOL, atol=ATOL)
        # a batch axis of 1 on the left operand and a full one on the right
        np.testing.assert_allclose(tl.matmul(w, x, workers=3, min_rows=1), w @ x, rtol=RTOL, atol=ATOL)

    def test_sandwich_equals_the_left_associated_triple_product(self):
        rng = np.random.default_rng(5)
        x = rng.standard_normal((6, 7, 3, 3)) + 1j * rng.standard_normal((6, 7, 3, 3))
        s = rng.standard_normal((7, 3, 3)) + 0j
        ref = s[None] @ x @ s[None]
        for workers in (1, 2, 5):
            with self.subTest(workers=workers):
                np.testing.assert_allclose(
                    tl.sandwich(s, x, s, workers=workers, min_rows=1), ref, rtol=RTOL, atol=ATOL)

    def test_sandwich_cuts_the_second_batch_axis(self):
        # one frequency and many sites: the chunks run along axis 1
        rng = np.random.default_rng(10)
        x = rng.standard_normal((1, 40, 3, 3)) + 1j * rng.standard_normal((1, 40, 3, 3))
        s = rng.standard_normal((40, 3, 3)) + 0j
        c = rng.standard_normal((40, 3, 3)) + 0j
        np.testing.assert_allclose(tl.sandwich(s, x, c, workers=4, min_rows=1),
                                   s[None] @ x @ c[None], rtol=RTOL, atol=ATOL)

    def test_stacks_shorter_than_the_worker_count(self):
        for n in (0, 1, 3):
            with self.subTest(n=n):
                a, b = _stack(n, 4, 11), _stack(n, 4, 12)
                np.testing.assert_allclose(tl.solve(a, b, workers=8, min_rows=1),
                                           np.linalg.solve(a, b), rtol=RTOL, atol=ATOL)
                np.testing.assert_allclose(tl.svdvals(a, workers=8, min_rows=1),
                                           np.linalg.svd(a, compute_uv=False),
                                           rtol=RTOL, atol=ATOL)

    def test_all_cores(self):
        a, b = _stack(300, 5, 13), _stack(300, 5, 14)
        np.testing.assert_allclose(tl.solve(a, b, workers=-1, min_rows=1),
                                   np.linalg.solve(a, b), rtol=RTOL, atol=ATOL)

    def test_a_singular_block_raises_linalgerror_from_a_worker(self):
        a = _stack(200, 4, 6)
        a[137] = 0.0
        with self.assertRaises(np.linalg.LinAlgError):
            tl.solve(a, a, workers=4, min_rows=8)

    def test_matmul_with_a_vector_operand_matches_numpy(self):
        rng = np.random.default_rng(8)
        m = rng.standard_normal((6, 4, 4)) + 0j
        v = rng.standard_normal(4) + 0j
        for workers in (1, 3):
            with self.subTest(workers=workers):
                np.testing.assert_allclose(tl.matmul(m, v, workers=workers, min_rows=1), m @ v,
                                           rtol=RTOL, atol=ATOL)
                np.testing.assert_allclose(tl.matmul(v, m, workers=workers, min_rows=1), v @ m,
                                           rtol=RTOL, atol=ATOL)

    @needs_pool
    def test_worker_failure_joins_every_task_and_restores_the_blas_limit(self):
        import threading
        import time
        from threadpoolctl import threadpool_info
        before = [(d["internal_api"], d["num_threads"]) for d in threadpool_info()]
        running = []

        def slow():
            running.append(1)
            time.sleep(0.3)
            running.pop()
            return 0

        def boom():
            raise RuntimeError("boom")

        with self.assertRaises(RuntimeError):
            tl._run([boom, slow, slow], 3)
        self.assertEqual(running, [])             # nothing still running after _run raised
        self.assertEqual(before, [(d["internal_api"], d["num_threads"]) for d in threadpool_info()])

    @needs_pool
    def test_concurrent_callers_are_serialized_and_restore_the_limit(self):
        import threading
        from threadpoolctl import threadpool_info
        before = [(d["internal_api"], d["num_threads"]) for d in threadpool_info()]
        a = _stack(64, 4, 9)
        ref = np.linalg.solve(a, a)
        out, errs = [], []

        def call():
            try:
                out.append(tl.solve(a, a, workers=2, min_rows=8))
            except Exception as e:                # pragma: no cover
                errs.append(e)

        ts = [threading.Thread(target=call) for _ in range(4)]
        [t.start() for t in ts]
        [t.join() for t in ts]
        self.assertEqual(errs, [])
        for r in out:
            np.testing.assert_allclose(r, ref, rtol=RTOL, atol=ATOL)
        self.assertEqual(before, [(d["internal_api"], d["num_threads"]) for d in threadpool_info()])

    def test_missing_threadpoolctl_warns_once_with_the_install_hint(self):
        saved = (tl.threadpool_limits, tl._warned)
        tl.threadpool_limits, tl._warned = None, False
        try:
            a = _stack(8, 3)
            with self.assertLogs("qlms.solver.threaded", level="WARNING") as logs:
                np.testing.assert_allclose(tl.solve(a, a, workers=4), np.linalg.solve(a, a),
                                           rtol=RTOL, atol=ATOL)
                tl.solve(a, a, workers=4)
            self.assertEqual(len(logs.records), 1)
            self.assertIn('pip install "hwave[threads]"', logs.output[0])
        finally:
            tl.threadpool_limits, tl._warned = saved

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
            _approx_equal(self, a[1:], b[1:])
        else:
            _approx_equal(self, a[2:], b[2:])      # the message text embeds rounded numbers

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
            np.testing.assert_allclose(a, b, rtol=RTOL, atol=ATOL)


if __name__ == "__main__":
    unittest.main()
