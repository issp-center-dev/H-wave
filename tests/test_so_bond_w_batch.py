"""The FLEX effective-interaction batch of the spin-orbital bond path and
two small helpers of the SO bond FLEX loop.

* ``so_bond.so_bond_w_batch(X, Y, gamma_b, plan=)``: ``W = P_B[Gamma_B
  (Y - X/2) Gamma_B]`` per slice and momentum, from a bare bubble batch
  ``X`` and its dressing ``Y``. Dense: an explicit-index reference.
  S_z-packed: the sector products ``M@s = Gamma@s (Y@s - X@s/2) Gamma@s``
  followed by ``so_bond.p_b``, which reads block ``s`` of ``W`` from block
  ``-s`` of ``M``; it equals the dense result restricted to the test-side
  sector lists. Two deliberate breaks (a same-sector ``p_b`` and the 1/2
  of ``Y - X/2`` set to 1, both patched at call time) miss by >= 1e-3.
  The inputs are never written; an ``out`` sharing memory with an input,
  a batch above the plan's capacity and non-finite inputs are refused.
* ``so_bond.channel0_block(X, nd=)``: a new dense copy of the pair block
  of channel 0 (from the sectors for a packed ``X``).
* ``so_bond.inverse_per_sector(A, norb_phys)``: the inverse of an S_z
  block-diagonal ``(..., 2n, 2n)`` stack, block by block, with exact zeros
  between the blocks; a wrong shape or a nonzero off-sector entry is
  refused.
* ``so_bond.static_pole_indicator(X, gamma_b, shape)``: the spinful pole
  indicator of a static bubble slice, without logging; its entry is the
  one the ``hwave_sc`` static bond check reports, in either layout.
"""
import logging
import unittest
from unittest import mock

import numpy as np

from hwave.solver import so_bond
from hwave.solver.hartree_fock import NonFiniteError
from hwave.solver.so_bond import PackedArrays, SOBondCaps, plan_so_bond
from tests.so_bond_fixtures import (BUBBLE_NMAT, bubble_green, chain_gamma, chain_lists,
                                    chain_plan, copy_packed, dressed, static_bubble)
from tests.so_bond_oracles import FIXTURE_LZSZ, rel, restrict

NOCAP = SOBondCaps(None, None, None)

# the dense random case: B = 2 channels, nd = 2 modes (ND = 8), 3 k points,
# a batch of 2 slices
B, ND_MODES, NVOL, NB = 2, 2, 3, 2
ND = B * ND_MODES ** 2


def _rand(rng, shape):
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


def _dense_case(seed=11):
    rng = np.random.default_rng(seed)
    X = _rand(rng, (NB, NVOL, ND, ND))
    Y = _rand(rng, (NB, NVOL, ND, ND))
    G = _rand(rng, (NVOL, ND, ND))
    plan = plan_so_bond(nmat=4, nvol=NVOL, B=B, nd=ND_MODES, consumers=("flex_scf",),
                        caps=NOCAP, nb=NB)
    return X, Y, G, plan


def _reference(X, Y, G, nd, half=0.5):
    """The explicit-index ``W[l, q, (m; x, y), (m'; z, w)] = M[l, q, (m; y,
    x), (m'; w, z)]`` of ``M[l, q, I, J] = sum_{K, L} G[q, I, K] (Y - half
    X)[l, q, K, L] G[q, L, J]``."""
    nb, nvol, n, _ = X.shape
    npair = nd * nd
    T = Y - half * X
    M = np.zeros_like(T)
    for l in range(nb):
        for q in range(nvol):
            for I in range(n):
                for J in range(n):
                    M[l, q, I, J] = sum(G[q, I, K] * T[l, q, K, L] * G[q, L, J]
                                        for K in range(n) for L in range(n))
    W = np.empty_like(M)
    for I in range(n):
        m, r = divmod(I, npair)
        x, y = divmod(r, nd)
        for J in range(n):
            m2, r2 = divmod(J, npair)
            z, w = divmod(r2, nd)
            W[..., I, J] = M[..., m * npair + y * nd + x, m2 * npair + w * nd + z]
    return W


def _lzsz():
    """``(fx, Xd, Xp, Gd, Gp, Yd, Yp, plan_dense, plan_packed)``: the dressed
    batch of four slices of FIXTURE_LZSZ in both layouts, and the FLEX plans."""
    fx = FIXTURE_LZSZ
    xd, xp, gd, gp, yd, yp = dressed(fx, "batch")
    plans = [chain_plan(fx, layout, ("flex_scf",), BUBBLE_NMAT[fx])
             for layout in ("dense", "sz_packed")]
    return (fx, xd, xp, gd, gp, yd, yp) + tuple(plans)


def _same_sector_p_b(M, nd, *, out=None):
    """The rejected reading of ``P_B`` in the packed layout: block ``s`` of
    the result read from block ``s`` of ``M`` through the swap permutation
    (instead of block ``-s``)."""
    sm = M.sector_map
    blocks = tuple(np.take(np.take(b, cross.perm, axis=-2), cross.perm, axis=-1)
                   for b, cross in zip(M.blocks, sm.pb_cross))
    if out is None:
        return PackedArrays(sm, blocks)
    for o, b in zip(out.blocks, blocks):
        o[...] = b
    return out


class TestWBatch(unittest.TestCase):

    def test_dense_matches_reference(self):
        X, Y, G, plan = _dense_case()
        want = _reference(X, Y, G, ND_MODES)
        got = so_bond.so_bond_w_batch(X, Y, G, plan=plan)
        self.assertEqual((got.shape, got.dtype), (X.shape, np.complex128))
        self.assertTrue(np.isfinite(got).all())
        self.assertLessEqual(rel(got, want), 1e-13)
        # into a caller's out (C-contiguous and strided)
        out = np.empty_like(X)
        self.assertIs(so_bond.so_bond_w_batch(X, Y, G, plan=plan, out=out), out)
        self.assertLessEqual(rel(out, want), 1e-13)
        big = np.zeros((NB, NVOL, ND, 2 * ND), dtype=np.complex128)
        strided = big[..., ::2]
        self.assertFalse(strided.flags.c_contiguous)
        self.assertIs(so_bond.so_bond_w_batch(X, Y, G, plan=plan, out=strided), strided)
        self.assertLessEqual(rel(strided, want), 1e-13)

    def test_packed_equals_dense(self):
        fx, xd, xp, gd, gp, yd, yp, pd, pp = _lzsz()
        self.assertEqual(xd.shape[0], 4)
        M = so_bond._w_kernel(xp, yp, gp, 0.5)
        scale = max(float(np.abs(b).max()) for b in M.blocks)
        self.assertGreaterEqual(float(np.abs(M[+1] - M[-1]).max()), 1e-2 * scale)
        wd = so_bond.so_bond_w_batch(xd, yd, gd, plan=pd)
        self.assertLessEqual(rel(wd, _reference_matmul(xd, yd, gd, fx.nd)), 1e-13)
        wp = so_bond.so_bond_w_batch(xp, yp, gp, plan=pp)
        self.assertIsInstance(wp, PackedArrays)
        self.assertIs(wp.sector_map, xp.sector_map)
        self.assertEqual(wp.lead_shape, xp.lead_shape)
        lists = chain_lists(fx)
        for s, (got, want) in enumerate(zip(wp.blocks, restrict(wd, lists))):
            self.assertGreater(np.abs(want).max(), 0.0)
            self.assertLessEqual(rel(got, want), 1e-13, s)
        # the dense W of block-diagonal inputs is exactly zero off the sectors
        mask = np.ones(wd.shape[-2:], dtype=bool)
        for l in lists:
            mask[l[:, None], l] = False
        self.assertTrue(np.all(wd[..., mask] == 0))
        # into a caller's packed out
        out = PackedArrays.zeros(xp.sector_map, xp.lead_shape)
        self.assertIs(so_bond.so_bond_w_batch(xp, yp, gp, plan=pp, out=out), out)
        for a, b in zip(out.blocks, wp.blocks):
            self.assertTrue(np.array_equal(a, b))

    def test_cross_sector_mutation(self):
        fx, xd, xp, gd, gp, yd, yp, pd, pp = _lzsz()
        lists = chain_lists(fx)
        want = restrict(_reference_matmul(xd, yd, gd, fx.nd), lists)   # before the patch
        spy = mock.Mock(wraps=_same_sector_p_b)
        with mock.patch.object(so_bond, "p_b", spy):
            self.assertEqual(spy.call_count, 0)
            wrong = so_bond.so_bond_w_batch(xp, yp, gp, plan=pp)
            self.assertGreaterEqual(spy.call_count, 1)
        self.assertTrue(all(np.isfinite(b).all() for b in wrong.blocks))
        for s in (1, 2):                                   # sectors +1 and -1
            self.assertGreaterEqual(rel(wrong.blocks[s], want[s]), 1e-3, s)
        # the unpatched call meets the reference
        right = so_bond.so_bond_w_batch(xp, yp, gp, plan=pp)
        for s in range(3):
            self.assertLessEqual(rel(right.blocks[s], want[s]), 1e-13, s)

    def test_half_mutation(self):
        fx, xd, xp, gd, gp, yd, yp, pd, pp = _lzsz()
        lists = chain_lists(fx)
        want = _reference_matmul(xd, yd, gd, fx.nd)                   # before the patch
        kernel = so_bond._w_kernel

        def full(X, Y, gamma_b, half):
            return kernel(X, Y, gamma_b, 2.0 * half)                   # Y - X for half = 1/2

        for layout, args, plan in (("dense", (xd, yd, gd), pd), ("sz_packed", (xp, yp, gp), pp)):
            with self.subTest(layout=layout):
                spy = mock.Mock(wraps=full)
                with mock.patch.object(so_bond, "_w_kernel", spy):
                    self.assertEqual(spy.call_count, 0)
                    got = so_bond.so_bond_w_batch(*args, plan=plan)
                    self.assertEqual(spy.call_count, 1)
                self.assertEqual(spy.call_args[0][3], 0.5)
                dense = got if layout == "dense" else got.to_dense()
                self.assertTrue(np.isfinite(dense).all())
                self.assertGreaterEqual(rel(dense, want), 1e-3)
                if layout == "sz_packed":
                    for s, w in enumerate(restrict(want, lists)):
                        self.assertGreaterEqual(rel(got.blocks[s], w), 1e-3, s)

    def test_inputs_not_written(self):
        X, Y, G, plan = _dense_case()
        copies = [a.copy() for a in (X, Y, G)]
        so_bond.so_bond_w_batch(X, Y, G, plan=plan)
        so_bond.so_bond_w_batch(X, Y, G, plan=plan, out=np.empty_like(X))
        for a, c in zip((X, Y, G), copies):
            self.assertTrue(np.array_equal(a, c))
        _, _, xp, _, gp, _, yp, _, pp = _lzsz()
        xp, yp, gp = (copy_packed(p) for p in (xp, yp, gp))   # private copies of the cache
        copies = [copy_packed(p) for p in (xp, yp, gp)]
        so_bond.so_bond_w_batch(xp, yp, gp, plan=pp)
        so_bond.so_bond_w_batch(xp, yp, gp, plan=pp,
                                out=PackedArrays.zeros(xp.sector_map, xp.lead_shape))
        for a, c in zip((xp, yp, gp), copies):
            for x, y in zip(a.blocks, c.blocks):
                self.assertTrue(np.array_equal(x, y))

    def test_out_alias_refused(self):
        X, Y, G, plan = _dense_case()
        buf = np.zeros((NB, NVOL, ND, ND), dtype=np.complex128)
        cases = (("X", (X, Y, G), X), ("Y", (X, Y, G), Y), ("gamma_b", (X, Y, buf[0]), buf))
        for name, (x, y, g), out in cases:
            with self.subTest(layout="dense", name=name):
                before = out.copy()
                with self.assertRaises(ValueError) as cm:
                    so_bond.so_bond_w_batch(x, y, g, plan=plan, out=out)
                self.assertEqual(str(cm.exception),
                                 "so_bond_w_batch: out shares memory with {}".format(name))
                self.assertTrue(np.array_equal(out, before))
        _, _, xp, _, gp, _, yp, _, pp = _lzsz()
        xp, yp = copy_packed(xp), copy_packed(yp)
        pbuf = PackedArrays.zeros(xp.sector_map, xp.lead_shape)
        g0 = PackedArrays(xp.sector_map, tuple(b[0] for b in pbuf.blocks))
        for b, g in zip(g0.blocks, gp.blocks):
            b[...] = g
        cases = (("X", (xp, yp, gp), xp), ("Y", (xp, yp, gp), yp),
                 ("gamma_b", (xp, yp, g0), pbuf))
        for name, (x, y, g), out in cases:
            with self.subTest(layout="sz_packed", name=name):
                with self.assertRaises(ValueError) as cm:
                    so_bond.so_bond_w_batch(x, y, g, plan=pp, out=out)
                self.assertEqual(str(cm.exception),
                                 "so_bond_w_batch: out shares memory with {}".format(name))

    def test_plan_and_input_refusals(self):
        X, Y, G, plan = _dense_case()
        # the plan's boundary first: a batch above its capacity, a plan without flex_scf
        small = plan_so_bond(nmat=4, nvol=NVOL, B=B, nd=ND_MODES, consumers=("flex_scf",),
                             caps=NOCAP, nb=1)
        with self.assertRaisesRegex(ValueError, r"^so_bond_w_batch: a batch of 2 slices "
                                                r"exceeds the plan's dress_capacity = 1 "):
            so_bond.so_bond_w_batch(X, Y, G, plan=small)
        dyn = plan_so_bond(nmat=4, nvol=NVOL, B=B, nd=ND_MODES, consumers=("chi_dynamic",),
                           caps=NOCAP, nb=2)
        with self.assertRaisesRegex(ValueError, r"^so_bond_w_batch: the admission plan has no "
                                                r"row 'w_batch_workspace'"):
            so_bond.so_bond_w_batch(X, Y, G, plan=dyn)
        with self.assertRaisesRegex(ValueError, r"^so_bond_w_batch: plan must be a SOBondPlan"):
            so_bond.so_bond_w_batch(X, Y, G, plan=None)
        bad = X.copy()
        bad[1, 2, 3, 4] = np.nan
        for name, args in (("X", (bad, Y, G)), ("Y", (X, bad, G)),
                           ("gamma_b", (X, Y, bad[1].copy()))):
            with self.subTest(name=name):
                with self.assertRaises(NonFiniteError) as cm:
                    so_bond.so_bond_w_batch(*args, plan=plan)
                self.assertEqual(str(cm.exception),
                                 "so_bond_w_batch: {} has non-finite entries".format(name))
        with self.assertRaisesRegex(ValueError, r"^so_bond_w_batch: X and Y must both be "):
            so_bond.so_bond_w_batch(X, Y[:, :2], G, plan=plan)
        with self.assertRaisesRegex(ValueError, r"^so_bond_w_batch: out must be "):
            so_bond.so_bond_w_batch(X, Y, G, plan=plan, out=np.empty(X.shape))
        # layouts: packed inputs with a dense plan, mixed layouts
        _, _, xp, _, gp, _, yp, pd, pp = _lzsz()
        with self.assertRaisesRegex(ValueError, r"^so_bond_w_batch: the admission plan has "
                                                r"layout = dense but the call has layout = "
                                                r"sz_packed"):
            so_bond.so_bond_w_batch(xp, yp, gp, plan=pd)
        with self.assertRaisesRegex(ValueError, r"^so_bond_w_batch: X, Y and gamma_b must all "):
            so_bond.so_bond_w_batch(xp, yp, gp.to_dense(), plan=pp)
        with self.assertRaisesRegex(ValueError, r"^so_bond_w_batch: X, Y and gamma_b must all "
                                                r".*\(and out too\)$"):
            so_bond.so_bond_w_batch(xp, yp, gp, plan=pp, out=xp.to_dense())
        with self.assertRaisesRegex(ValueError, r"^so_bond_w_batch: out must be PackedArrays on "
                                                r"the sector map of X"):
            so_bond.so_bond_w_batch(xp, yp, gp, plan=pp,
                                    out=PackedArrays.zeros(xp.sector_map, (1,) + xp.lead_shape))


def _reference_matmul(X, Y, G, nd):
    """The dense W of the fixture batch by plain matrix products and the
    explicit slot swap (the explicit-index sums of :func:`_reference` are
    too slow at ND = 80)."""
    M = G[None] @ (Y - 0.5 * X) @ G[None]
    n = X.shape[-1]
    npair = nd * nd
    I = np.arange(n)
    m, r = np.divmod(I, npair)
    x, y = np.divmod(r, nd)
    f = m * npair + y * nd + x
    return M[..., f[:, None], f]


class TestReferences(unittest.TestCase):

    def test_matmul_reference_is_the_index_reference(self):
        X, Y, G, _ = _dense_case(seed=5)
        self.assertLessEqual(rel(_reference_matmul(X, Y, G, ND_MODES),
                                 _reference(X, Y, G, ND_MODES)), 1e-13)


class TestChannel0(unittest.TestCase):

    def test_dense_and_packed(self):
        fx, xd, xp, *_ = _lzsz()
        P = fx.nd ** 2
        dense = so_bond.channel0_block(xd, nd=fx.nd)
        self.assertEqual(dense.shape, xd.shape[:-2] + (P, P))
        self.assertTrue(np.array_equal(dense, xd[..., :P, :P]))
        self.assertFalse(np.shares_memory(dense, xd))
        self.assertTrue(dense.flags.c_contiguous)
        packed = so_bond.channel0_block(xp)
        self.assertEqual((packed.dtype, packed.shape), (np.complex128, dense.shape))
        self.assertTrue(np.array_equal(packed, dense))
        self.assertTrue(np.array_equal(so_bond.channel0_block(xp, nd=fx.nd), dense))
        self.assertFalse(any(np.shares_memory(packed, b) for b in xp.blocks))
        # off-sector entries of the channel-0 block are exact zeros, on-sector not all
        lists = [l[l < P] for l in chain_lists(fx)]
        mask = np.ones((P, P), dtype=bool)
        for l in lists:
            mask[l[:, None], l] = False
        self.assertTrue(np.all(packed[..., mask] == 0))
        self.assertGreater(np.abs(packed[..., ~mask]).min(), 0.0)
        # a static slice (one leading axis)
        one = PackedArrays(xp.sector_map, tuple(b[0] for b in xp.blocks))
        self.assertTrue(np.array_equal(so_bond.channel0_block(one), dense[0]))

    def test_refusals(self):
        fx, xd, xp, *_ = _lzsz()
        with self.assertRaisesRegex(ValueError, r"^channel0_block: a dense X needs nd "):
            so_bond.channel0_block(xd)
        with self.assertRaisesRegex(ValueError, r"^channel0_block: the sector map of X has "
                                                r"nd = 4, the call has nd = 2"):
            so_bond.channel0_block(xp, nd=2)
        with self.assertRaisesRegex(ValueError, r"^channel0_block: X must be a complex128 "):
            so_bond.channel0_block(xd[..., :5], nd=fx.nd)


class TestInversePerSector(unittest.TestCase):

    @staticmethod
    def _blockdiag(rng, lead, n):
        A = np.zeros(lead + (2 * n, 2 * n), dtype=np.complex128)
        eye = np.eye(n)
        A[..., :n, :n] = _rand(rng, lead + (n, n)) + 3 * eye
        A[..., n:, n:] = _rand(rng, lead + (n, n)) + 3 * eye
        return A

    def test_matches_inv(self):
        rng = np.random.default_rng(3)
        for n, lead in ((1, (5,)), (2, (3, 4)), (3, (2, 1, 3))):
            with self.subTest(n=n, lead=lead):
                A = self._blockdiag(rng, lead, n)
                before = A.copy()
                got = so_bond.inverse_per_sector(A, n)
                self.assertTrue(np.array_equal(A, before))
                self.assertEqual((got.shape, got.dtype), (A.shape, np.complex128))
                self.assertLessEqual(rel(got, np.linalg.inv(A)), 1e-13)
                self.assertTrue(np.all(got[..., :n, n:] == 0))
                self.assertTrue(np.all(got[..., n:, :n] == 0))
        # a single matrix (no leading axis)
        A = self._blockdiag(rng, (), 2)
        self.assertLessEqual(rel(so_bond.inverse_per_sector(A, 2), np.linalg.inv(A)), 1e-13)

    def test_shape_and_leak_refused(self):
        rng = np.random.default_rng(4)
        A = self._blockdiag(rng, (3,), 2)
        for bad in (A[..., :3, :3], A[..., :3], A[0, 0]):
            with self.subTest(shape=bad.shape):
                with self.assertRaises(ValueError) as cm:
                    so_bond.inverse_per_sector(bad, 2)
                self.assertEqual(str(cm.exception), (
                    "inverse_per_sector: A must end in shape (2*norb_phys, 2*norb_phys) = "
                    "(4, 4), got shape {}".format(bad.shape)))
        for n in (0, -1, 1.5, True):
            with self.subTest(norb_phys=n):
                with self.assertRaisesRegex(ValueError, r"^inverse_per_sector: norb_phys must "
                                                        r"be a positive integer"):
                    so_bond.inverse_per_sector(A, n)
        for i, j in ((0, 3), (2, 1)):
            with self.subTest(entry=(i, j)):
                leak = A.copy()
                leak[1, i, j] = 1e-14
                want = 1e-14 / float(np.abs(leak).max())
                with self.assertRaises(ValueError) as cm:
                    so_bond.inverse_per_sector(leak, 2)
                self.assertEqual(str(cm.exception), "inverse_per_sector: off-sector leak "
                                 "max|A_off| / max|A| = %.3e" % want)


class TestPoleIndicator(unittest.TestCase):

    def test_matches_sc(self):
        from hwave.solver.eliashberg import bond as sc_bond
        fx = FIXTURE_LZSZ
        g = bubble_green(fx)
        shape = tuple(fx.spatial_shape)
        entries = {}
        for layout in ("dense", "sz_packed"):
            with self.subTest(layout=layout):
                X = static_bubble(fx, layout, g)
                gamma = chain_gamma(fx, layout)
                with self.assertNoLogs(level=logging.DEBUG):
                    entry, per, worst = so_bond.static_pole_indicator(X, gamma, shape)
                report = sc_bond._so_bond_static_pole(X, gamma, shape)
                self.assertEqual(entry, report.indicators["spinful"])
                self.assertTrue(np.isfinite(entry[0]))
                if layout == "dense":
                    self.assertIsNone(per)
                    self.assertIsNone(worst)
                else:
                    self.assertEqual(set(per) - set(X.sector_map.sector_order), set())
                    self.assertIn(worst, per)
                    self.assertEqual(per[worst], entry)
                    self.assertEqual(entry[0], max(v[0] for v in per.values()))
                entries[layout] = entry
        # the sector spectra together are the dense block's
        self.assertAlmostEqual(entries["sz_packed"][0], entries["dense"][0], delta=1e-10)


if __name__ == "__main__":
    unittest.main()
