"""Pinned staging and threaded host copies of the bond block store (GitHub
issue #275, part 2).

A host store on the GPU path moves device blocks and batches through one
bounded pinned buffer, chunk by chunk along the frequency axis, and splits
its pair-strided host copies over threads; a device store leaves the device
through the same buffer. These tests force tiny chunks (several per array,
the last one partial) and threaded copies of tiny arrays, on the fake device
of tests/fake_device.py, and check that every byte arrives exactly as on the
previous path -- the plain ``asnumpy`` / ``asarray`` / assignment, which is
also what runs when the pinned allocation fails.
"""
import contextlib
import logging
import unittest
from unittest import mock

import numpy as np

from tests.fake_device import DevArray, fake_cupy, is_device


def _host(x):
    return np.array(x.view(np.ndarray) if isinstance(x, DevArray) else x, copy=True)


def _rand(shape, seed):
    rng = np.random.default_rng(seed)
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


@contextlib.contextmanager
def _threaded_tiny_copies():
    """Every host copy of the staging module runs threaded, however small."""
    from hwave.solver import staging
    with mock.patch.object(staging, "PARALLEL_COPY_MIN_BYTES", 0):
        yield


class TestParallelCopy(unittest.TestCase):

    def test_strided_bitwise_and_more_workers_than_rows(self):
        from hwave.solver import staging
        big = _rand((5, 3, 8, 8), 1)
        for workers in (1, 2, 3, 8):
            with self.subTest(workers=workers):
                dst = np.zeros_like(big)
                staging.parallel_copy(dst[:, :, 4:8, 0:4], big[:, :, 0:4, 4:8], workers,
                                      min_bytes=0)
                ref = np.zeros_like(big)
                ref[:, :, 4:8, 0:4] = big[:, :, 0:4, 4:8]
                np.testing.assert_array_equal(dst, ref)

    def test_overlapping_operands_are_one_assignment(self):
        from hwave.solver import staging
        x = np.arange(40.0).reshape(10, 4)
        ref = x.copy()
        ref[1:] = ref[:-1].copy()
        with mock.patch.object(staging, "_executor",
                               side_effect=AssertionError("threads on overlapping operands")):
            staging.parallel_copy(x[1:], x[:-1], 4, min_bytes=0)
        np.testing.assert_array_equal(x, ref)

    def test_the_real_threshold(self):
        """Just below PARALLEL_COPY_MIN_BYTES: one assignment; at and above it:
        threaded pieces. Bitwise either way."""
        from hwave.solver import staging
        limit = staging.PARALLEL_COPY_MIN_BYTES
        self.assertEqual(limit, 16 * 1024 ** 2)
        row = 64 * 1024                                  # bytes per leading-axis row
        real = staging._executor
        for nbytes, threaded in ((limit - row, False), (limit, True), (limit + row, True)):
            with self.subTest(nbytes=nbytes):
                n = nbytes // row
                src = np.arange(n * row // 8, dtype=np.float64).reshape(n, row // 8)
                dst = np.zeros_like(src)
                self.assertEqual(dst.nbytes, nbytes)
                with mock.patch.object(staging, "_executor", side_effect=real) as ex:
                    staging.parallel_copy(dst, src, 4)
                self.assertEqual(ex.called, threaded)
                np.testing.assert_array_equal(dst, src)

    def test_concurrent_callers_share_one_pool(self):
        """Callers asking for different worker counts at the same time never
        meet a pool shut down under them."""
        import threading
        from hwave.solver import staging
        errors = []

        def run(w):
            try:
                for _ in range(20):
                    src = _rand((16, 4, 4), w)
                    dst = np.empty_like(src)
                    staging.parallel_copy(dst, src, w, min_bytes=0)
                    np.testing.assert_array_equal(dst, src)
            except Exception as exc:              # noqa: BLE001
                errors.append(exc)
        ts = [threading.Thread(target=run, args=(w,)) for w in (2, 3, 5, 7, 9, 16)]
        for t in ts:
            t.start()
        for t in ts:
            t.join()
        self.assertEqual(errors, [])

    def test_default_workers_follow_the_affinity_mask(self):
        from hwave.solver import staging
        with mock.patch.object(staging.os, "sched_getaffinity", create=True,
                               return_value={0, 1, 2}), \
                mock.patch.object(staging.os, "cpu_count", return_value=64):
            self.assertEqual(staging.default_copy_workers(), 3)
        with mock.patch.object(staging.os, "sched_getaffinity", create=True,
                               return_value=set(range(32))):
            self.assertEqual(staging.default_copy_workers(), staging.COPY_WORKERS)
        with mock.patch.object(staging.os, "sched_getaffinity", create=True,
                               side_effect=AttributeError), \
                mock.patch.object(staging.os, "cpu_count", return_value=2):
            self.assertEqual(staging.default_copy_workers(), 2)
        with mock.patch.object(staging.os, "sched_getaffinity", create=True,
                               side_effect=OSError("no affinity")), \
                mock.patch.object(staging.os, "cpu_count", return_value=None):
            self.assertEqual(staging.default_copy_workers(), 1)
        saved = getattr(staging.os, "sched_getaffinity", None)
        try:
            if saved is not None:
                del staging.os.sched_getaffinity      # absent, as on macOS
            with mock.patch.object(staging.os, "cpu_count", return_value=5):
                self.assertEqual(staging.default_copy_workers(), 5)
            with mock.patch.object(staging.os, "cpu_count", return_value=None):
                self.assertEqual(staging.default_copy_workers(), 1)
        finally:
            if saved is not None:
                staging.os.sched_getaffinity = saved

    def test_shape_mismatch_refused(self):
        from hwave.solver import staging
        with self.assertRaises(ValueError):
            staging.parallel_copy(np.zeros((2, 3)), np.zeros((3, 2)), 2, min_bytes=0)


class TestPinnedStager(unittest.TestCase):

    def test_chunks_through_one_reused_buffer(self):
        from hwave.solver import staging
        shape = (7, 3, 4, 4)                       # one row: 3 * 16 * 16 = 768 bytes
        row = 3 * 4 * 4 * 16
        x = _rand(shape, 2)
        with fake_cupy() as cp, _threaded_tiny_copies():
            st = staging.PinnedStager(cp, max_bytes=2 * row + 5, workers=3)
            d = cp.asarray(x)
            cp.reset()
            # device -> strided host view
            big = np.zeros((7, 3, 8, 8), complex)
            st.copy_to_host(d, big[:, :, 4:8, 0:4])
            np.testing.assert_array_equal(big[:, :, 4:8, 0:4], x)
            self.assertEqual(np.count_nonzero(big[:, :, :4]), 0)
            self.assertEqual(st.chunks, 4)                         # 2 + 2 + 2 + 1 rows
            self.assertEqual([s[0] for _d, s, _n in cp.transfers], [2, 2, 2, 1])
            self.assertEqual(sum(n for _d, _s, n in cp.transfers), x.nbytes)
            # strided host -> new device array
            back = st.to_device(big[:, :, 4:8, 0:4])
            self.assertTrue(is_device(back) and back.flags.c_contiguous)
            np.testing.assert_array_equal(_host(back), x)
            # new host array
            h = st.to_host(d)
            self.assertFalse(is_device(h))
            np.testing.assert_array_equal(h, x)
            self.assertEqual(st.chunks, 12)
            # every chunk went through the pinned buffer, allocated once
            self.assertEqual(len(cp.pinned_io), 12)
            self.assertTrue(all(p for _d, p in cp.pinned_io))
            self.assertEqual(cp.pinned_allocs, [2 * row])
            self.assertEqual(st.nbytes, 2 * row)
            st.release()
            self.assertEqual(st.nbytes, 0)
            self.assertEqual(cp.pinned_frees, 1)

    def test_asynchronous_copies_are_synchronized_per_chunk(self):
        """With copies that only land at their stream's synchronize (an
        asynchronous DMA from / to pinned memory), several chunks through
        one buffer still give every byte: each chunk's stream is synchronized
        before the buffer is read or overwritten."""
        from hwave.solver import staging
        shape = (7, 3, 4, 4)
        row = 3 * 4 * 4 * 16
        x = _rand(shape, 6)
        with fake_cupy() as cp, _threaded_tiny_copies():
            d = cp.asarray(x)
            cp.deferred = True
            st = staging.PinnedStager(cp, max_bytes=2 * row, workers=2)
            big = np.zeros((7, 3, 8, 8), complex)
            st.copy_to_host(d, big[:, :, 0:4, 4:8])                   # D2H, 4 chunks
            for i in range(7):
                np.testing.assert_array_equal(big[i, :, 0:4, 4:8], x[i], err_msg=str(i))
            np.testing.assert_array_equal(st.to_host(d), x)
            back = st.to_device(big[:, :, 0:4, 4:8])                  # H2D, 4 chunks
            self.assertEqual(cp.stream.pending, [])
            for i in range(7):
                np.testing.assert_array_equal(_host(back[i]), x[i], err_msg=str(i))
            self.assertEqual(st.chunks, 12)
            self.assertGreaterEqual(cp.stream.syncs, 12)
            # a non-default current stream is the one used and synchronized
            other = cp.cuda.Stream(non_blocking=True)
            with mock.patch.object(cp.cuda, "get_current_stream", lambda: other):
                back = st.to_device(x[::-1])
            self.assertEqual(other.syncs, 4)
            np.testing.assert_array_equal(_host(back), x[::-1])

    def test_a_slice_larger_than_the_buffer_is_cut_within(self):
        """One frequency slice (512 bytes) larger than the buffer (200): the
        chunks run along the first axis whose sub-arrays fit (axis 2, rows
        of 64 bytes, 3 per chunk) and the pinned buffer stays within
        max_bytes; copies deferred to the stream, host side strided."""
        from hwave.solver import staging
        x = _rand((3, 2, 4, 4), 3)
        with fake_cupy() as cp, _threaded_tiny_copies():
            d = cp.asarray(x)
            cp.deferred = True
            st = staging.PinnedStager(cp, max_bytes=200, workers=2)
            cp.reset()
            big = np.zeros((3, 2, 8, 4), complex)
            st.copy_to_host(d, big[:, :, 2:6])
            np.testing.assert_array_equal(big[:, :, 2:6], x)
            self.assertEqual(np.count_nonzero(big[:, :, :2]) + np.count_nonzero(big[:, :, 6:]), 0)
            self.assertEqual([s for _d, s, _n in cp.transfers], [(3, 4), (1, 4)] * 6)
            back = st.to_device(big[:, :, 2:6])
            np.testing.assert_array_equal(_host(back), x)
            self.assertEqual(st.chunks, 24)
            self.assertEqual(cp.pinned_allocs, [3 * 4 * 16])
            self.assertLessEqual(st.nbytes, 200)
            self.assertTrue(all(n <= 200 for n in cp.pinned_allocs))
            # a buffer smaller than one element (8 < 16 bytes): the unstaged
            # path, no pinned memory at all
            cp.pinned_allocs.clear()
            st1 = staging.PinnedStager(cp, max_bytes=8)
            np.testing.assert_array_equal(st1.to_host(d), x)
            dst = np.zeros((3, 2, 4, 8), complex)
            st1.copy_to_host(d, dst[..., ::2])
            np.testing.assert_array_equal(dst[..., ::2], x)
            np.testing.assert_array_equal(_host(st1.to_device(dst[..., ::2])), x)
            self.assertEqual(st1.chunks, 0)
            self.assertEqual(st1.nbytes, 0)
            self.assertEqual(cp.pinned_allocs, [])

    def test_concurrent_callers_share_one_stager(self):
        """The stager's lock: threads moving different arrays through one
        buffer (deferred copies, several chunks each) all get their bytes."""
        import threading
        from hwave.solver import staging
        with fake_cupy() as cp:
            xs = [_rand((6, 2, 4, 4), 30 + i) for i in range(6)]
            ds = [cp.asarray(x) for x in xs]
            cp.deferred = True
            st = staging.PinnedStager(cp, max_bytes=2 * 2 * 4 * 4 * 16)
            out, errors = [None] * 6, []

            def run(i):
                try:
                    for _ in range(10):
                        out[i] = (st.to_host(ds[i]), _host(st.to_device(xs[i])))
                except Exception as exc:          # noqa: BLE001
                    errors.append(exc)
            ts = [threading.Thread(target=run, args=(i,)) for i in range(6)]
            for t in ts:
                t.start()
            for t in ts:
                t.join()
            self.assertEqual(errors, [])
            for x, (h, d) in zip(xs, out):
                np.testing.assert_array_equal(h, x)
                np.testing.assert_array_equal(d, x)

    def test_release_during_threaded_shared_copies(self):
        """Copies through the global worker pool (threshold 0) on one shared
        stager while another thread keeps releasing its buffer: no deadlock,
        every byte right (release waits for the copy in flight)."""
        import threading
        from hwave.solver import staging
        with fake_cupy() as cp, _threaded_tiny_copies():
            xs = [_rand((8, 2, 4, 4), 40 + i) for i in range(4)]
            ds = [cp.asarray(x) for x in xs]
            cp.deferred = True
            st = staging.PinnedStager(cp, max_bytes=3 * 2 * 4 * 4 * 16, workers=3)
            errors, stop = [], threading.Event()
            real = staging._executor
            used = []

            def spy():
                used.append(1)
                return real()

            def run(i):
                try:
                    for _ in range(8):
                        np.testing.assert_array_equal(st.to_host(ds[i]), xs[i])
                        big = np.zeros((8, 2, 4, 8), complex)
                        st.copy_to_host(ds[i], big[..., 1::2])
                        np.testing.assert_array_equal(big[..., 1::2], xs[i])
                        np.testing.assert_array_equal(_host(st.to_device(big[..., 1::2])), xs[i])
                except Exception as exc:          # noqa: BLE001
                    errors.append(exc)

            def releaser():
                while not stop.is_set():
                    st.release()
            with mock.patch.object(staging, "_executor", spy):
                ts = [threading.Thread(target=run, args=(i,)) for i in range(4)]
                r = threading.Thread(target=releaser)
                r.start()
                for t in ts:
                    t.start()
                for t in ts:
                    t.join(timeout=60)
                stop.set()
                r.join(timeout=60)
            self.assertFalse(any(t.is_alive() for t in ts) or r.is_alive(), "deadlock")
            self.assertEqual(errors, [])
            self.assertTrue(used)                 # the copies used the worker pool

    def test_release_waits_for_the_copy_in_flight(self):
        """Deterministically: a copy holding the stager is paused inside its
        host-side copy; release() from another thread demonstrably tries the
        stager's lock meanwhile and must not complete until the copy
        finishes; the copied bytes are right."""
        import threading
        from hwave.solver import staging

        class _Watched:
            """The stager's RLock, signalling when the release thread tries it."""

            def __init__(self, lock):
                self.lock, self.owner, self.tried = lock, None, threading.Event()

            def __enter__(self):
                if threading.current_thread() is self.owner:
                    self.tried.set()
                return self.lock.__enter__()

            def __exit__(self, *exc):
                return self.lock.__exit__(*exc)

        with fake_cupy() as cp:
            x = _rand((6, 2, 4, 4), 50)
            d = cp.asarray(x)
            st = staging.PinnedStager(cp, max_bytes=2 * 2 * 4 * 4 * 16, workers=1)
            st._lock = watched = _Watched(st._lock)
            inside, go = threading.Event(), threading.Event()
            real = staging.parallel_copy
            out = {}

            def paused(dst, src, *a, **k):
                if not inside.is_set():
                    inside.set()
                    go.wait(timeout=30)
                return real(dst, src, *a, **k)

            with mock.patch.object(staging, "parallel_copy", paused):
                c = threading.Thread(target=lambda: out.setdefault("h", st.to_host(d)))
                r = threading.Thread(target=st.release)
                watched.owner = r
                try:
                    c.start()
                    self.assertTrue(inside.wait(timeout=30))
                    r.start()
                    self.assertTrue(watched.tried.wait(timeout=30),
                                    "release did not try the stager's lock")
                    r.join(timeout=0.2)
                    waited = r.is_alive()
                finally:
                    go.set()
                    c.join(timeout=30)
                    r.join(timeout=30)
            self.assertTrue(waited, "release did not wait for the copy")
            self.assertFalse(c.is_alive() or r.is_alive(), "deadlock")
            np.testing.assert_array_equal(out["h"], x)

    def test_failed_pinned_allocation_takes_the_plain_path(self):
        from hwave.solver import staging
        x = _rand((5, 2, 4, 4), 4)
        with fake_cupy() as cp:
            cp.fail_pinned = True
            st = staging.PinnedStager(cp, max_bytes=64)
            d = cp.asarray(x)
            cp.reset()
            with self.assertLogs("qlms.solver.staging", logging.DEBUG) as log:
                h = st.to_host(d)
            self.assertIn("unstaged path", "\n".join(log.output))
            np.testing.assert_array_equal(h, x)
            dst = np.zeros_like(x)
            st.copy_to_host(d, dst)
            np.testing.assert_array_equal(dst, x)
            back = st.to_device(x[:, :, ::-1])
            np.testing.assert_array_equal(_host(back), x[:, :, ::-1])
            self.assertEqual(st.chunks, 0)
            self.assertEqual(cp.pinned_io, [])
            # whole-array copies, as before
            self.assertEqual([(t, s) for t, s, _n in cp.transfers],
                             [("d2h", x.shape), ("d2h", x.shape), ("h2d", x.shape)])
            self.assertFalse(st.available)
            # each call tries again (and says so)
            with self.assertLogs("qlms.solver.staging", logging.DEBUG) as log:
                st.copy_to_host(d, dst)
                st.to_device(x)
            self.assertEqual(sum("unstaged path" in r.getMessage() for r in log.records), 2)
            # once pinned memory is available again, the next call is staged
            cp.fail_pinned = False
            np.testing.assert_array_equal(st.to_host(d), x)
            self.assertTrue(st.available)
            self.assertGreater(st.chunks, 0)
            self.assertTrue(all(p for _d, p in cp.pinned_io))


class TestStoreStaging(unittest.TestCase):
    """Every store method that crosses host <-> device gives the previous
    path's bytes with several chunks and threaded copies."""

    def _stores(self, cp, half, staging_bytes):
        from hwave.solver.flex_bond import BondBlockStore
        nmat, nvol, nd, B = 8, 3, 4, 2
        kw = dict(half_grid=half)
        new = BondBlockStore(nmat, nvol, B * nd, nd, ("chibar", "W"),
                             staging_bytes=staging_bytes, copy_workers=3, **kw)
        old = BondBlockStore(nmat, nvol, B * nd, nd, ("chibar", "W"), **kw)
        return new, old

    def test_host_store_bitwise_against_plain_copies(self):
        from hwave.solver import half_grid as hg
        from tests.test_bond_half_grid import _hermitian_in_frequency
        nmat, nvol, nd, B = 8, 3, 4, 2
        row_pair = nvol * nd * nd * 16
        for half in (False, True):
            with self.subTest(half=half), fake_cupy() as cp, _threaded_tiny_copies():
                new, old = self._stores(cp, half, staging_bytes=2 * row_pair)
                with new, old:
                    ns = new.nslices
                    W = _hermitian_in_frequency(nmat, nvol, B * nd, seed=5)
                    Wh = hg.half_from_full(W, nmat) if half else W
                    # put_freq_batch of a device batch (two calls, partial chunks)
                    for l0, l1 in ((0, 3), (3, ns)):
                        new.put_freq_batch("W", l0, l1, cp.asarray(Wh[l0:l1]))
                        old._arr("W")[l0:l1] = Wh[l0:l1]
                    # put_pair of device blocks into chibar
                    for a in range(B):
                        for b in range(B):
                            blk = _rand((ns, nvol, nd, nd), 10 * a + b)
                            new.put_pair("chibar", a, b, cp.asarray(blk))
                            old._arr("chibar")[:, :, a * nd:(a + 1) * nd,
                                               b * nd:(b + 1) * nd] = blk
                    for name in ("chibar", "W"):
                        np.testing.assert_array_equal(new._arr(name), old._arr(name))
                    self.assertGreater(new.stager(cp).chunks, 2 * B * B)
                    # pair blocks to the device and to the host
                    for a in range(B):
                        for b in range(B):
                            ref = np.ascontiguousarray(old.get_pair("W", a, b)) if not half \
                                else hg.expand(old.get_pair("W", a, b), nmat,
                                               ba=old.get_pair("W", b, a))
                            out = new.get_pair_full("W", a, b, xp=cp)
                            self.assertTrue(is_device(out) and out.flags.c_contiguous)
                            np.testing.assert_array_equal(_host(out), ref)
                            out_h = new.get_pair_full("W", a, b)
                            self.assertFalse(is_device(out_h))
                            self.assertTrue(out_h.flags.c_contiguous)
                            np.testing.assert_array_equal(out_h, ref)
                    # frequency batches to the device; on the store's own
                    # module the view itself
                    b_dev = new.get_freq_batch_on("chibar", 1, ns, cp)
                    self.assertTrue(is_device(b_dev))
                    np.testing.assert_array_equal(_host(b_dev), old._arr("chibar")[1:ns])
                    v = new.get_freq_batch_on("chibar", 1, 3, np)
                    self.assertTrue(np.shares_memory(v, new._arr("chibar")))
                    self.assertTrue(all(p for _d, p in cp.pinned_io))
                    # bounded by staging_bytes, although a batch's frequency
                    # slice (4 * row_pair) is larger
                    self.assertLessEqual(new.stager(cp).nbytes, 2 * row_pair)
                    self.assertTrue(all(n <= 2 * row_pair for n in cp.pinned_allocs))
                    # the cast / broadcast path is the plain assignment
                    new.put_freq_batch("W", 0, 1, 0.5)
                    self.assertTrue(np.all(new._arr("W")[0] == 0.5))
                self.assertIsNone(new._stager)
                self.assertGreaterEqual(cp.pinned_frees, 1)

    def test_device_store_leaves_in_chunks(self):
        from hwave.solver import backend as bk
        from hwave.solver.flex_bond import BondBlockStore
        nmat, nvol, nd, B = 9, 3, 4, 2
        ND = B * nd
        slot = nvol * ND * ND * 16
        for fail in (False, True):
            with self.subTest(fail_pinned=fail), fake_cupy() as cp, _threaded_tiny_copies():
                cp.fail_pinned = fail
                with BondBlockStore(nmat, nvol, ND, nd, ("chibar", "W"), xp=cp,
                                    staging_bytes=4 * slot, copy_workers=2) as st:
                    for i, name in enumerate(("chibar", "W")):
                        st.put_freq_batch(name, 0, nmat, _rand((nmat, nvol, ND, ND), i))
                    ref = {n: bk.to_host(st._arr(n)) for n in ("chibar", "W")}
                    cp.reset()
                    st.to_host()
                    self.assertEqual(st.resident, "host")
                    for n in ref:
                        a = st._arr(n)
                        self.assertFalse(is_device(a))
                        self.assertTrue(a.flags.c_contiguous)
                        np.testing.assert_array_equal(a, ref[n])
                    self.assertEqual(sum(b for d, _s, b in cp.transfers if d == "d2h"),
                                     2 * nmat * slot)
                    rows = [s[0] for _d, s, _b in cp.transfers]
                    self.assertEqual(rows, [nmat, nmat] if fail else [4, 4, 1, 4, 4, 1])


class TestStoreLeavesTheDeviceAsynchronously(unittest.TestCase):

    def test_to_host_with_deferred_copies(self):
        from hwave.solver.flex_bond import BondBlockStore
        nmat, nvol, nd, B = 9, 3, 4, 2
        ND = B * nd
        with fake_cupy() as cp:
            with BondBlockStore(nmat, nvol, ND, nd, ("chibar", "W"), xp=cp,
                                staging_bytes=2 * nvol * ND * ND * 16) as st:
                ref = {}
                for i, name in enumerate(("chibar", "W")):
                    ref[name] = _rand((nmat, nvol, ND, ND), 20 + i)
                    st.put_freq_batch(name, 0, nmat, ref[name])
                cp.deferred = True
                st.to_host()
                for name in ref:
                    np.testing.assert_array_equal(st._arr(name), ref[name])


class TestCopyWorkersFollowTheBackend(unittest.TestCase):
    """A CPU solve copies in one thread (as before); a GPU solve with a host
    store gets the threaded copies."""

    def test_cpu_one_thread_gpu_host_store_threaded(self):
        from hwave.solver import flex_bond
        from tests.test_bond_half_grid import _run
        seen = []
        real = flex_bond.BondBlockStore

        class Spy(real):
            def __init__(self, *a, **k):
                super().__init__(*a, **k)
                seen.append((self.copy_workers, self.resident))
        with mock.patch.object(flex_bond, "BondBlockStore", Spy):
            _run("UV", {}, half=False)
            cpu = list(seen)
            seen.clear()
            with fake_cupy():
                _run("UV", {"gpu": True, "longitudinal_bond_store": "host"}, half=False)
        self.assertEqual(cpu, [(1, "host")])
        # the store's own default is the single plain assignment
        self.assertEqual(real(2, 1, 2, 1, ("W",)).copy_workers, 1)
        self.assertEqual(seen, [(None, "host")])


class TestRealCupyStreams(unittest.TestCase):
    """On a CUDA device: forced multi-chunk host <-> device staging on a
    NON-default stream, every chunk verified bit for bit (the pinned copies
    are asynchronous there)."""

    @classmethod
    def setUpClass(cls):
        try:
            import cupy
            if cupy.cuda.runtime.getDeviceCount() < 1:
                raise RuntimeError("no device")
        except Exception as exc:                  # noqa: BLE001
            raise unittest.SkipTest("CuPy with a CUDA device is needed ({})".format(exc))
        cls.cp = cupy

    def test_multi_chunk_on_a_non_default_stream(self):
        from hwave.solver import staging
        cp = self.cp
        n, rest = 33, (64, 16, 16)                       # 256 KiB per row
        row = int(np.prod(rest)) * 16
        x = _rand((n,) + rest, 7)
        st = staging.PinnedStager(cp, max_bytes=3 * row, workers=4)
        try:
            with mock.patch.object(staging, "PARALLEL_COPY_MIN_BYTES", 0), \
                    cp.cuda.Stream(non_blocking=True) as s:
                # H2D from a strided host view
                big = np.zeros((n, 64, 32, 16), complex)
                big[:, :, 8:24, :] = x
                d = st.to_device(big[:, :, 8:24, :])
                # queue device work on the stream so the D2H must wait for it
                d2 = d * 2.0
                d2 -= d
                host = st.to_host(d2)
                strided = np.zeros((n, 64, 16, 32), complex)
                st.copy_to_host(d2, strided[:, :, :, 4:20])
                s.synchronize()
            # three transfers of n rows, at most three rows per chunk
            self.assertGreaterEqual(st.chunks, 3 * -(-n // 3))
            for i in range(n):
                np.testing.assert_array_equal(cp.asnumpy(d[i]), x[i], err_msg="h2d %d" % i)
                np.testing.assert_array_equal(host[i], x[i], err_msg="d2h %d" % i)
                np.testing.assert_array_equal(strided[i, :, :, 4:20], x[i],
                                              err_msg="strided d2h %d" % i)
        finally:
            st.release()

    def test_slices_larger_than_the_buffer_on_a_non_default_stream(self):
        from hwave.solver import staging
        cp = self.cp
        n, rest = 5, (64, 16, 16)                        # 256 KiB per slice
        row = int(np.prod(rest)) * 16
        x = _rand((n,) + rest, 8)
        st = staging.PinnedStager(cp, max_bytes=row // 3, workers=4)
        try:
            with mock.patch.object(staging, "PARALLEL_COPY_MIN_BYTES", 0), \
                    cp.cuda.Stream(non_blocking=True) as s:
                d = st.to_device(x[:, :, :, ::-1])
                host = st.to_host(d)
                strided = np.zeros((n, 64, 16, 32), complex)
                st.copy_to_host(d, strided[:, :, :, 1::2])        # non-contiguous dst
                s.synchronize()
            self.assertLessEqual(st.nbytes, row // 3)
            self.assertGreaterEqual(st.chunks, 2 * 3 * n)
            for i in range(n):
                np.testing.assert_array_equal(host[i], x[i, :, :, ::-1], err_msg="slice %d" % i)
                np.testing.assert_array_equal(strided[i, :, :, 1::2], x[i, :, :, ::-1],
                                              err_msg="strided slice %d" % i)
                self.assertEqual(np.count_nonzero(strided[i, :, :, 0::2]), 0)
        finally:
            st.release()


class TestSolvesOldVersusNew(unittest.TestCase):
    """Whole bond-gate solves on the fake device: tiny pinned chunks with
    threaded copies give the outputs of the unstaged path (pinned allocation
    refused) bit for bit, for a host store and for a device store that
    leaves the device before the pairing kernel."""

    def _run(self, store, pairing, half, staged):
        from hwave.solver import staging
        from tests.test_flex_bond_device_store import _solve
        seen = []

        def setup(cp):
            cp.fail_pinned = not staged
            seen.append(cp)
        with mock.patch.object(staging, "PINNED_STAGING_BYTES", 4096), \
                mock.patch.object(staging, "PARALLEL_COPY_MIN_BYTES", 0):
            out = _solve(store, pairing, half, setup=setup)
        return out, seen[0]

    def _check(self, store, pairing, half):
        from tests.test_flex_bond_device_store import _assert_same_outputs
        new, cp_new = self._run(store, pairing, half, True)
        old, cp_old = self._run(store, pairing, half, False)
        self.assertEqual(new[0]._bond_store_resident, store)
        _assert_same_outputs(self, old, new)
        staged = [p for _d, p in cp_new.pinned_io]
        self.assertGreater(len(staged), 0)
        self.assertTrue(all(staged))
        self.assertEqual(cp_old.pinned_io, [])

    def test_host_store_full_grid(self):
        self._check("host", "none", False)

    def test_host_store_half_grid_singlet(self):
        self._check("host", "singlet", True)

    def test_device_store_singlet(self):
        self._check("device", "singlet", False)


if __name__ == "__main__":
    unittest.main()
