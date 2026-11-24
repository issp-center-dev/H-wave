"""S_z sectors of the spin-orbital bond objects: the sector map and packed arrays.

* ``so_bond.sz_sector_map`` lists the code slots ``(m; x, y)`` of each ph
  sector ``sigma = S_z(x) - S_z(y)`` in the order ``("0", "+1", "-1")``
  (ascending dense index inside a sector), the cross-sector permutations of
  ``P_B`` (``(m; x, y) -> (m; y, x)``) and of the exchange image
  (``(m; x, y) -> (reverse[m]; y, x)``), and the gather tables of the
  pp-sector pairing vertex.
* ``so_bond.PackedArrays`` holds the three diagonal sector blocks of a
  code-slot array; ``from_dense`` requires exact zeros off the sectors.

Every map is checked against an independent test-side enumeration
(``tests.so_bond_oracles``) or against the dense operation it stands for.
"""
import dataclasses
import itertools
import types
import unittest

import numpy as np

from hwave.solver import so_bond
from hwave.solver.so_bond import PackedArrays, SectorMap, sz_sector_map
from tests.so_bond_fixtures import exact
from tests.so_bond_oracles import (PH_SIGMAS, PP_ORDER, ph_sector_lists, pp_sector_lists,
                                   restrict, sz_of_mode)

#: on-site, +-x, +-y channels; ``reverse`` pairs each bond with its partner
_DR = ((0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0))
_REV = (0, 2, 1, 4, 3)


def _topo(B):
    return types.SimpleNamespace(delta_r=np.array(_DR[:B], dtype=np.int64),
                                 reverse=np.array(_REV[:B], dtype=np.int64))


def _random_block_diagonal(lists, lead, ND, seed=0):
    """A random complex code-slot array whose off-sector entries are exactly 0."""
    rng = np.random.default_rng(seed)
    M = np.zeros(tuple(lead) + (ND, ND), dtype=np.complex128)
    for l in lists:
        n = len(l)
        M[..., l[:, None], l] = (rng.normal(size=tuple(lead) + (n, n))
                                 + 1j * rng.normal(size=tuple(lead) + (n, n)))
    return M


def _slot(m, x, y, nd):
    return m * nd * nd + x * nd + y


NORB, B = 2, 3
SM = sz_sector_map(_topo(B), 2 * NORB)
LISTS = ph_sector_lists(B, 2 * NORB)


class TestSectorMap(unittest.TestCase):
    """Construction check of the sector map: sizes, offsets, index lists,
    inverse maps, dtypes, read-only arrays, the cache and the refusals."""

    def test_lists_sizes_offsets(self):
        for norb, B_ in itertools.product((1, 2, 3), (1, 3, 5)):
            with self.subTest(norb=norb, B=B_):
                nd = 2 * norb
                sm = sz_sector_map(_topo(B_), nd)
                lists = ph_sector_lists(B_, nd)
                self.assertEqual(sm.sizes, (B_ * 2 * norb ** 2, B_ * norb ** 2, B_ * norb ** 2))
                self.assertEqual((sm.B, sm.nd, sm.norb_phys, sm.ND), (B_, nd, norb, B_ * nd * nd))
                self.assertEqual(sm.offsets.tolist(), [0, sm.sizes[0], sm.sizes[0] + sm.sizes[1],
                                                       sm.ND])
                for s in range(3):
                    self.assertTrue(np.array_equal(
                        sm.index[sm.offsets[s]:sm.offsets[s + 1]], lists[s]))
                # the inverse maps: sector sigma and local position of every dense slot
                for s, sigma in enumerate(PH_SIGMAS):
                    l = lists[s]
                    self.assertTrue(np.all(sm.sector_of[l] == sigma))
                    self.assertTrue(np.array_equal(sm.local[l], np.arange(len(l))))

    def test_fields_and_dtypes(self):
        self.assertEqual(SM.sector_order, ("0", "+1", "-1"))
        self.assertEqual(SM.sigmas, (0, 1, -1))
        self.assertEqual(SM.pp_order, (1, 0, -1))
        for name, dt, shape in (("offsets", np.int64, (4,)), ("index", np.int64, (SM.ND,)),
                                ("sector_of", np.int8, (SM.ND,)),
                                ("local", np.int64, (SM.ND,))):
            arr = getattr(SM, name)
            self.assertEqual((arr.dtype, arr.shape), (np.dtype(dt), shape), name)
            self.assertFalse(arr.flags.writeable, name)
        arrays = [c.perm for c in SM.pb_cross + SM.gp_cross] + list(SM.pp_pairs)
        for g in SM.pp_gather:
            arrays += [g.sector, g.row_local, g.col_local]
        self.assertTrue(all(not a.flags.writeable for a in arrays))
        with self.assertRaises(dataclasses.FrozenInstanceError):
            SM.B = 4

    def test_cache_by_content(self):
        nd = 2 * NORB
        self.assertIs(sz_sector_map(_topo(B), nd), SM)
        other = types.SimpleNamespace(delta_r=np.array(_DR[:3]), reverse=[0, 2, 1])
        self.assertIs(sz_sector_map(other, nd), SM)        # list reverse, same content
        self.assertIsNot(sz_sector_map(_topo(B), 2), SM)
        self.assertIsNot(sz_sector_map(_topo(5), nd), sz_sector_map(_topo(3), nd))

    def test_refusals(self):
        for nd in (0, 3, -2, 2.0, "4", True):
            with self.subTest(nd=nd), self.assertRaisesRegex(ValueError, exact(
                    "sz_sector_map: nd must be a positive even integer (nd = 2*norb_phys), "
                    "got {!r}".format(nd))):
                sz_sector_map(_topo(3), nd)
        for rev in ([0, 1, 1], [0, 2, 3], [[0, 2, 1]], [0.0, 2.0, 1.0]):
            with self.subTest(reverse=rev), self.assertRaisesRegex(ValueError, exact(
                    "sz_sector_map: topo.reverse must be an integer (B,) involution with "
                    "entries in [0, B), got {}".format(np.asarray(rev).tolist()))):
                sz_sector_map(types.SimpleNamespace(delta_r=np.zeros((3, 3), int),
                                                    reverse=np.asarray(rev)), 4)

    def test_consumer_map_is_identity(self):
        self.assertEqual(so_bond._MAP_CONSUMERS, ("gamma_packer", "gamma_p", "store_scatter",
                                                  "dress_sectors", "w_cross", "kernel_gather"))
        for name in so_bond._MAP_CONSUMERS:
            self.assertIs(so_bond._consumer_map(name, SM), SM)
        with self.assertRaisesRegex(ValueError, exact(
                "_consumer_map: unknown consumer 'gamma', expected one of ('gamma_packer', "
                "'gamma_p', 'store_scatter', 'dress_sectors', 'w_cross', 'kernel_gather')")):
            so_bond._consumer_map("gamma", SM)


class TestSectorMapGates(unittest.TestCase):
    """Gate P-map: the packed layout and every index map reproduce the dense
    operation they stand for, exactly (packed == dense bit for bit)."""

    def test_round_trip_and_off_sector_zero(self):
        M = _random_block_diagonal(LISTS, (4, 5), SM.ND)
        P = PackedArrays.from_dense(M, SM, source="test")
        ref = restrict(M, LISTS)
        for s in range(3):
            self.assertTrue(np.array_equal(P.blocks[s], ref[s]))
        self.assertTrue(np.array_equal(P.to_dense(), M))
        self.assertEqual(P.lead_shape, (4, 5))

    def test_storage_factor_from_allocations(self):
        for norb, B_ in ((1, 1), (2, 3), (3, 5)):
            sm = sz_sector_map(_topo(B_), 2 * norb)
            self.assertEqual(16 * PackedArrays.zeros(sm, (4, 5)).nbytes,
                             6 * np.zeros((4, 5, sm.ND, sm.ND), complex).nbytes)

    def test_pb_cross_is_pair_swap(self):
        """``W@sigma = P_B(M@-sigma)``: ``P_B[M][(m;x,y),(m';z,w)] =
        M[(m;y,x),(m';w,z)]`` with no phase."""
        for norb, B_ in ((1, 3), (2, 3), (2, 5)):
            nd = 2 * norb
            sm = sz_sector_map(_topo(B_), nd)
            lists = ph_sector_lists(B_, nd)
            sw = np.array([_slot(m, y, x, nd) for m in range(B_) for x in range(nd)
                           for y in range(nd)])
            M = _random_block_diagonal(lists, (2,), sm.ND, seed=1)
            W = M[..., sw[:, None], sw]
            P = PackedArrays.from_dense(M, sm, source="test")
            ref = restrict(W, lists)
            for t, sigma in enumerate(PH_SIGMAS):
                src, perm = sm.pb_cross[t]
                self.assertEqual(src, -sigma)
                self.assertTrue(np.array_equal(P[src][..., perm[:, None], perm], ref[t]))

    def test_gp_cross_is_reversed_transpose(self):
        """``Gamma^P[(m;x,y),(m';z,w)] = Gamma[(rev m'; w,z),(rev m; y,x)]``
        (phase and momentum aside): the slot map ``(m;x,y) -> (rev m; y,x)``
        and a transpose."""
        for norb, B_ in ((1, 3), (2, 5)):
            nd = 2 * norb
            sm = sz_sector_map(_topo(B_), nd)
            lists = ph_sector_lists(B_, nd)
            g = np.array([_slot(_REV[m], y, x, nd) for m in range(B_) for x in range(nd)
                          for y in range(nd)])
            G = _random_block_diagonal(lists, (3,), sm.ND, seed=2)
            GP = G[..., g[None, :], g[:, None]]
            P = PackedArrays.from_dense(G, sm, source="test")
            ref = restrict(GP, lists)
            for t, sigma in enumerate(PH_SIGMAS):
                src, perm = sm.gp_cross[t]
                self.assertEqual(src, -sigma)
                got = np.swapaxes(P[src][..., perm[:, None], perm], -1, -2)
                self.assertTrue(np.array_equal(got, ref[t]))

    def test_pp_pairs_and_gather(self):
        """``V^p_{ab}[i, j] = V[(a; A, B), (b; C, D)]`` for the output pair
        ``i = (A, D)`` and the input pair ``j = (B, C)`` of pp sector ``p``,
        read from the packed sector blocks; every in-sector entry is read by
        exactly one pp sector."""
        for norb, B_ in ((1, 1), (1, 3), (2, 3)):
            nd = 2 * norb
            sm = sz_sector_map(_topo(B_), nd)
            lists = ph_sector_lists(B_, nd)
            pp = pp_sector_lists(nd)
            V = _random_block_diagonal(lists, (), sm.ND, seed=3)
            P = PackedArrays.from_dense(V, sm, source="test")
            hits = [np.zeros((n, n), dtype=np.int64) for n in sm.sizes]
            for k, p in enumerate(PP_ORDER):
                pairs = sm.pp_pairs[k]
                self.assertTrue(np.array_equal(pairs, pp[k]))
                self.assertEqual(pairs.dtype, np.int64)
                n = len(pairs)
                self.assertEqual(n, norb ** 2 * (2 if p == 0 else 1))
                g = sm.pp_gather[k]
                self.assertEqual((g.sector.shape, g.sector.dtype), ((n, n), np.dtype(np.int8)))
                self.assertEqual(g.row_local.shape, (B_, n, n))
                self.assertEqual(g.col_local.shape, (B_, n, n))
                expect = np.empty((B_, B_, n, n), dtype=np.complex128)
                for a, b, i, j in itertools.product(range(B_), range(B_), range(n), range(n)):
                    A, D = divmod(int(pairs[i]), nd)
                    Bm, C = divmod(int(pairs[j]), nd)
                    self.assertEqual(sz_of_mode(A, norb) + sz_of_mode(D, norb), p)
                    expect[a, b, i, j] = V[_slot(a, A, Bm, nd), _slot(b, C, D, nd)]
                got = np.empty_like(expect)
                for s, sigma in enumerate(PH_SIGMAS):
                    mask = g.sector == sigma
                    for a, b in itertools.product(range(B_), range(B_)):
                        r, c = g.row_local[a][mask], g.col_local[b][mask]
                        got[a, b][mask] = P[sigma][r, c]
                        np.add.at(hits[s], (r, c), 1)
                self.assertTrue(np.array_equal(got, expect))
            for s in range(3):
                self.assertTrue(np.all(hits[s] == 1), (norb, B_, s))


class TestPackedArrays(unittest.TestCase):
    """Construction check of ``PackedArrays``: the invariants, item access,
    ownership, ``to_dense(out=)`` and the ``from_dense`` refusals."""

    def _blocks(self, lead=(2,), sm=SM):
        return tuple(np.zeros(tuple(lead) + (n, n), dtype=np.complex128) for n in sm.sizes)

    def test_packed_arrays_invariants(self):
        b0, bp, bm = self._blocks()
        n0, n1 = SM.sizes[0], SM.sizes[1]
        cases = (
            ((bp, b0, bm), "PackedArrays: block 0 (sector '0') must end in shape ({0}, {0}), "
                           "got shape (2, {1}, {1})".format(n0, n1)),
            ((b0, bp.real.copy(), bm), "PackedArrays: block 1 (sector '+1') must be a "
                                       "complex128 ndarray, got float64"),
            ((b0, bp, bm.tolist()), "PackedArrays: block 2 (sector '-1') must be a complex128 "
                                    "ndarray, got list"),
            ((b0, bp, bm[..., :-1]), "PackedArrays: block 2 (sector '-1') must end in shape "
                                     "({0}, {0}), got shape (2, {0}, {1})".format(n1, n1 - 1)),
            ((b0, np.zeros((3, n1, n1), complex), bm),
             "PackedArrays: block 1 (sector '+1') has leading shape (3,), block 0 has (2,)"),
            ((b0, bp), "PackedArrays: blocks must be a tuple of exactly 3 arrays (sectors '0', "
                       "'+1', '-1'), got 2 blocks"),
            ([b0, bp, bm], "PackedArrays: blocks must be a tuple of exactly 3 arrays (sectors "
                           "'0', '+1', '-1'), got list"),
        )
        for blocks, msg in cases:
            with self.subTest(msg=msg), self.assertRaisesRegex(ValueError, exact(msg)):
                PackedArrays(SM, blocks)
        with self.assertRaisesRegex(ValueError, exact(
                "PackedArrays: sector_map must be a SectorMap, got SimpleNamespace")):
            PackedArrays(types.SimpleNamespace(), (b0, bp, bm))
        for field, value, msg in (
                ("sector_order", ("0", "+1"), "PackedArrays: sector_map.sector_order must be "
                 "('0', '+1', '-1'), got ('0', '+1')"),
                ("sector_order", ("0", "-1", "+1"), "PackedArrays: sector_map.sector_order "
                 "must be ('0', '+1', '-1'), got ('0', '-1', '+1')"),
                ("sigmas", (0, -1, 1), "PackedArrays: sector_map.sigmas must be (0, 1, -1), "
                 "got (0, -1, 1)"),
                ("sizes", SM.sizes[:2], "PackedArrays: sector_map.sizes must have 3 entries, "
                 "got {!r}".format(SM.sizes[:2]))):
            bad = dataclasses.replace(SM, **{field: value})
            with self.subTest(field=field, value=value), \
                    self.assertRaisesRegex(ValueError, exact(msg)):
                PackedArrays(bad, (b0, bp, bm))
        P = PackedArrays(SM, (b0, bp, bm))                 # valid: holds the given arrays
        self.assertIs(P.blocks[1], bp)
        self.assertEqual(P.lead_shape, (2,))
        self.assertEqual(P.nbytes, b0.nbytes + bp.nbytes + bm.nbytes)

    def test_item_access(self):
        P = PackedArrays(SM, self._blocks())
        for k, s in ((0, 0), ("0", 0), (1, 1), (+1, 1), ("+1", 1), (-1, 2), ("-1", 2)):
            self.assertIs(P[k], P.blocks[s])
        for k in ("1", 2, True, 1.0, None):
            with self.subTest(key=k), self.assertRaisesRegex(KeyError, "PackedArrays: sector"):
                P[k]

    def test_zeros_and_from_dense_own_their_blocks(self):
        Z = PackedArrays.zeros(SM, (2, 3))
        self.assertEqual(Z.lead_shape, (2, 3))
        for b, n in zip(Z.blocks, SM.sizes):
            self.assertEqual((b.shape, b.dtype), ((2, 3, n, n), np.dtype(np.complex128)))
            self.assertTrue(b.flags.owndata and b.flags.writeable and not b.any())
        M = _random_block_diagonal(LISTS, (2,), SM.ND)
        P = PackedArrays.from_dense(M, SM, source="test")
        for b in P.blocks:
            self.assertTrue(b.flags.owndata and b.flags.c_contiguous)
            self.assertFalse(np.shares_memory(b, M))
        P.blocks[0][...] = 0                               # never reaches the source
        self.assertTrue(np.array_equal(P.to_dense()[..., LISTS[1][:, None], LISTS[1]],
                                       M[..., LISTS[1][:, None], LISTS[1]]))
        R = PackedArrays.from_dense(M.real, SM, source="test")   # real input is promoted
        self.assertEqual(R.blocks[0].dtype, np.complex128)

    def test_to_dense_out(self):
        M = _random_block_diagonal(LISTS, (2,), SM.ND)
        P = PackedArrays.from_dense(M, SM, source="test")
        D1, D2 = P.to_dense(), P.to_dense()
        self.assertFalse(np.shares_memory(D1, D2))
        out = np.full((2, SM.ND, SM.ND), 7.0 + 1j)
        self.assertIs(P.to_dense(out=out), out)
        self.assertTrue(np.array_equal(out, M))
        for bad in (np.zeros((2, SM.ND, SM.ND)), np.zeros((3, SM.ND, SM.ND), complex)):
            with self.assertRaisesRegex(ValueError, exact(
                    "PackedArrays.to_dense: out must be a complex128 array of shape "
                    "(2, {0}, {0}), got {1} {2}".format(SM.ND, bad.dtype, bad.shape))):
                P.to_dense(out=bad)

    def test_from_dense_shape_refusal(self):
        with self.assertRaisesRegex(ValueError, exact(
                "PackedArrays.from_dense: 'chi' must end in shape (ND, ND) = ({0}, {0}), "
                "got shape (2, {0}, {1})".format(SM.ND, SM.ND - 1))):
            PackedArrays.from_dense(np.zeros((2, SM.ND, SM.ND - 1), complex), SM, source="chi")

    def test_from_dense_refusals(self):
        I, J = int(LISTS[0][3]), int(LISTS[1][2])          # an off-sector entry (0, +1)
        for eps, cls in ((1e-14, "numerical leakage from a non-canonical producer"),
                         (1e-12, "numerical leakage from a non-canonical producer"),
                         (1e-6, "incompatible Hamiltonian (produced from a different, "
                                "S_z-breaking input)")):
            M = 0.5 * _random_block_diagonal(LISTS, (2,), SM.ND, seed=4) / 10.0
            M[0, int(LISTS[0][0]), int(LISTS[0][1])] = 1.0          # max |M| = 1
            M[1, I, J] = eps
            leak = np.max(np.abs(M[1, I, J])) / np.max(np.abs(M))
            msg = ("'chi_w.npy' is dense with nonzero S_z-off-sector entries (measured leak "
                   "= %.3e; packing requires exactly 0) although this run conserves S_z: "
                   "%s. Regenerate it with this version." % (leak, cls))
            with self.subTest(eps=eps), self.assertRaisesRegex(ValueError, exact(msg)):
                PackedArrays.from_dense(M, SM, source="chi_w.npy")
        M = _random_block_diagonal(LISTS, (), SM.ND)
        nan_msg = ("'x' is dense with nonzero S_z-off-sector entries (measured leak = nan; "
                   "packing requires exactly 0) although this run conserves S_z: incompatible "
                   "Hamiltonian (produced from a different, S_z-breaking input). Regenerate it "
                   "with this version.")
        # a non-finite off-sector entry in the rows of each sector
        for r, c in ((I, J), (int(LISTS[1][0]), int(LISTS[2][1])),
                     (int(LISTS[2][3]), int(LISTS[0][2]))):
            for bad in (np.nan, np.inf):
                M[r, c] = bad
                with self.subTest(row=r, col=c, value=bad), \
                        self.assertRaisesRegex(ValueError, exact(nan_msg)):
                    PackedArrays.from_dense(M, SM, source="x")
                M[r, c] = 0.0
        # a finite off-sector leak with a non-finite in-sector entry: leak undefined
        M[int(LISTS[0][0]), int(LISTS[0][0])] = np.inf
        M[I, J] = 1e-14
        with self.assertRaisesRegex(ValueError, exact(nan_msg)):
            PackedArrays.from_dense(M, SM, source="x")
        M = _random_block_diagonal(LISTS, (), SM.ND)
        M[J, I] = -0.0                                     # an exact (negative) zero is a zero
        PackedArrays.from_dense(M, SM, source="test")


if __name__ == "__main__":
    unittest.main()
