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

    def test_a_row_larger_than_the_buffer_is_one_chunk(self):
        from hwave.solver import staging
        x = _rand((3, 2, 4, 4), 3)
        with fake_cupy() as cp:
            st = staging.PinnedStager(cp, max_bytes=16)
            np.testing.assert_array_equal(st.to_host(cp.asarray(x)), x)
            self.assertEqual(st.chunks, 3)
            self.assertEqual(cp.pinned_allocs, [2 * 4 * 4 * 16])

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
                    self.assertEqual(len(cp.pinned_allocs), 1)
                    # bounded by max(staging_bytes, one frequency row of a batch)
                    self.assertEqual(new.stager(cp).nbytes, max(2 * row_pair, B * B * row_pair))
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
