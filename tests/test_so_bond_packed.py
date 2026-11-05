"""The S_z-sector packed layout of the spinful bond vertex.

* ``so_bond.so_bond_vertex`` with a plan of layout ``"sz_packed"`` returns
  ``PackedArrays`` ``(nvol, n_sigma, n_sigma)``: the diagonal ph-sector
  blocks (``sigma = S_z(x) - S_z(y)`` of the code slot ``(m; x, y)``) of
  the dense ``Gamma_B``. A record that couples two sectors is refused
  before anything is allocated.
* ``so_bond.gamma_p`` takes packed input and output: the exchange image
  ``(m; x, y) -> (reverse[m]; y, x)`` maps sector ``sigma`` onto ``-sigma``.
* ``so_bond.so_bond_bubble`` writes the static bubble into a packed sink and
  every bosonic slice into a ``so_bond.PackedBondBlockStore``
  (``so_bond.open_store``), from a Green function built per spin block
  (``tests.so_bond_fixtures.canonical_green``); a Green function with an
  entry between the spin blocks is refused.
* ``so_bond.so_bond_dress`` dresses packed bubbles sector by sector; its
  conditioning guard scores the extreme singular values over all sectors
  (the dense block's) and a refusal names the sectors where they lie.
* ``so_bond.so_bond_pairing_vertex_parts`` and
  ``so_bond.so_bond_pairing_vertex_batch`` build the pairing vertex
  ``1/4 (Gamma + Gamma^P) - Gamma chi Gamma`` sector by sector (static
  parts; one dressed frequency batch written into the caller's ``out``).
* ``so_bond.p_b`` reads block ``sigma`` of ``W`` from block ``-sigma`` of
  ``M``; the FLEX ``M = Gamma (Y - X/2) Gamma`` of a dressed batch has
  different +1 and -1 blocks, so that reading is visible.

Every packed result is compared with the dense production result restricted
to the test-side sector lists (``tests.so_bond_oracles.ph_sector_lists``).
"""
import unittest
from unittest import mock

import numpy as np

import re
import tracemalloc
import types

from hwave.solver import bond_channels, flex_bond, so_bond
from hwave.solver.hartree_fock import NonFiniteError
from hwave.solver.so_bond import (PackedArrays, PackedBondBlockStore, _file_records, gamma_p,
                                  open_store, plan_so_bond, slot_geometry, so_bond_bubble,
                                  so_bond_topology, so_bond_vertex, sz_sector_map)
from hwave.solver.so_records import compile_interaction
from tests.so_bond_fixtures import BUBBLE_NMAT, NOCAP, assert_recomputed
from tests.so_bond_fixtures import FIXTURE_NAMES as NAMES
from tests.so_bond_fixtures import batch_plan as _batch_plan
from tests.so_bond_fixtures import bubble_green as _green
from tests.so_bond_fixtures import chain_gamma as _gamma
from tests.so_bond_fixtures import chain_lists as _lists
from tests.so_bond_fixtures import chain_plan as _plan
from tests.so_bond_fixtures import chain_sector_map as _sm
from tests.so_bond_fixtures import chain_setup as _setup
from tests.so_bond_fixtures import copy_packed as _copy_packed
from tests.so_bond_fixtures import dress as _dress
from tests.so_bond_fixtures import dress_inputs as _dress_inputs
from tests.so_bond_fixtures import dress_plan as _dress_plan
from tests.so_bond_fixtures import dressed as _dressed
from tests.so_bond_fixtures import every_slice_bubble as _every_slice
from tests.so_bond_fixtures import exact
from tests.so_bond_fixtures import packed_lead as _packed_lead
from tests.so_bond_fixtures import pairing_parts as _parts
from tests.so_bond_fixtures import parts_plan as _parts_plan
from tests.so_bond_fixtures import production_bubbles as _bubbles
from tests.so_bond_fixtures import random_packed as _random_packed
from tests.so_bond_fixtures import static_bubble as _static
from tests.so_bond_fixtures import vertex_batch as _batch
from tests.so_bond_oracles import (ALL_CONSERVING_FIXTURES, CELL, FIXTURE_I, FIXTURE_LZSZ,
                                   FIXTURE_SZ_GRID, rel, restrict, sz_of_mode)


def _off_sector(M, lists):
    """The entries of ``M`` outside the diagonal sector blocks."""
    mask = np.ones(M.shape[-2:], dtype=bool)
    for l in lists:
        mask[l[:, None], l] = False
    return M[..., mask]


class TestPackedGamma(unittest.TestCase):
    """Gate P-map and P-alg for Gamma_B and its exchange image: the packed
    blocks equal the dense production result restricted to the test-side
    sector lists (<= 1e-13 relative, max-norm); the +1 and -1 blocks agree
    to 1e-15 relative (their records accumulate in different orders)."""

    def test_packed_equals_dense_restricted(self):
        for fx in ALL_CONSERVING_FIXTURES:
            with self.subTest(fixture=NAMES[fx]):
                d, p = _gamma(fx, "dense"), _gamma(fx, "sz_packed")
                self.assertIsInstance(d, np.ndarray)
                self.assertIsInstance(p, PackedArrays)
                geom = _setup(fx)[2]
                self.assertIs(p.sector_map, sz_sector_map(geom, fx.nd))
                self.assertEqual(p.lead_shape, (fx.nvol,))
                lists = _lists(fx)
                r = restrict(d, lists)
                for s in range(3):
                    self.assertEqual(p.blocks[s].shape, (fx.nvol, len(lists[s]), len(lists[s])))
                    self.assertGreater(np.abs(r[s]).max(), 0.0)
                    self.assertLessEqual(rel(p.blocks[s], r[s]), 1e-13)
                self.assertEqual(np.abs(_off_sector(d, lists)).max(), 0.0)

    def test_plus_minus_blocks_equal(self):
        for fx in ALL_CONSERVING_FIXTURES:
            with self.subTest(fixture=NAMES[fx]):
                g = _gamma(fx, "sz_packed")
                self.assertGreater(np.abs(g[+1]).max(), 0.0)
                self.assertLessEqual(rel(g[+1], g[-1]), 1e-15)

    def test_gamma_p_packed(self):
        for fx in ALL_CONSERVING_FIXTURES:
            with self.subTest(fixture=NAMES[fx]):
                geom = _setup(fx)[2]
                d, p = _gamma(fx, "dense"), _gamma(fx, "sz_packed")
                lists = _lists(fx)
                dp = restrict(gamma_p(d, geom, out=np.empty_like(d)), lists)
                out = PackedArrays.zeros(p.sector_map, p.lead_shape)
                self.assertIs(gamma_p(p, geom, out=out), out)
                for s in range(3):
                    self.assertLessEqual(rel(out.blocks[s], dp[s]), 1e-13)
                    self.assertLessEqual(rel(out.blocks[s], p.blocks[s]), 1e-13)

    def test_gamma_p_packed_random(self):
        """A random block-diagonal input, whose exchange image differs from
        it and whose +1 and -1 blocks differ: the cross-sector reading of
        the packed image is visible."""
        for fx in (FIXTURE_SZ_GRID, FIXTURE_LZSZ):
            with self.subTest(fixture=NAMES[fx]):
                geom = _setup(fx)[2]
                sm = sz_sector_map(geom, fx.nd)
                x = _random_packed(sm, fx.nvol, 4)
                self.assertGreaterEqual(rel(x[+1], x[-1]), 1e-2)
                dense = x.to_dense()
                lists = _lists(fx)
                ref = gamma_p(dense, geom, out=np.empty_like(dense))
                self.assertEqual(np.abs(_off_sector(ref, lists)).max(), 0.0)
                ref = restrict(ref, lists)
                out = gamma_p(x, geom, out=PackedArrays.zeros(sm, (fx.nvol,)))
                for s in range(3):
                    self.assertGreaterEqual(rel(ref[s], x.blocks[s]), 1e-2)
                    self.assertLessEqual(rel(out.blocks[s], ref[s]), 1e-13)

    def test_sector_coupling_record_refused(self):
        """The first record that couples two ph sectors is named (norb_phys 1:
        the seven-type fixture; norb_phys 2: an on-site inter-orbital
        PairLift); nothing is allocated."""
        onsite_pl = {"CoulombIntra": {((0, 0, 0), (0, 0)): 1.0, ((0, 0, 0), (1, 1)): 0.8},
                     "PairLift": {((0, 0, 0), (0, 1)): 0.3, ((0, 0, 0), (1, 0)): 0.3}}
        name = {0: "0", 1: "+1", -1: "-1"}
        for tables, norb in ((FIXTURE_I.tables_dict, 1), (onsite_pl, 2)):
            with self.subTest(norb_phys=norb):
                c = compile_interaction(tables, norb)
                topo = so_bond_topology(c, CELL)
                rec = next(r for r in c.records
                           if sum(sgn * sz_of_mode(x, norb)
                                  for sgn, x in zip((1, -1, 1, -1), r["modes"])))
                x1, x2, x3, x4 = (int(x) for x in rec["modes"])
                s_row = name[int(round(sz_of_mode(x1, norb) - sz_of_mode(x2, norb)))]
                s_col = name[int(round(sz_of_mode(x4, norb) - sz_of_mode(x3, norb)))]
                self.assertNotEqual(s_row, s_col)
                plan = plan_so_bond(nmat=2, nvol=5, B=len(topo.delta_r), nd=2 * norb,
                                    consumers=("vertex",), caps=NOCAP, layout="sz_packed",
                                    layout_reason="so bond layout: sz_packed (forced for a test)")
                msg = ("so_bond_vertex: a {} record couples the S_z sectors {} and {} in the "
                       "code slots; the packed layout needs an S_z-conserving interaction "
                       "(select_layout)".format(str(rec["type"]), s_row, s_col))
                with mock.patch.object(so_bond.PackedArrays, "zeros") as zeros:
                    with self.assertRaisesRegex(ValueError, exact(msg)):
                        so_bond_vertex(c, topo, (5, 1, 1), plan=plan)
                zeros.assert_not_called()
                self.assertEqual(str(rec["type"]), "PairLift")


class TestPackedGammaRefusals(unittest.TestCase):
    """Construction check of the packed filing and exchange image: the
    refusals of mixed or mismatched inputs (exact texts)."""

    @classmethod
    def setUpClass(cls):
        cls.fx = FIXTURE_SZ_GRID
        cls.geom = _setup(cls.fx)[2]
        cls.sm = sz_sector_map(cls.geom, cls.fx.nd)
        cls.p = _gamma(cls.fx, "sz_packed")

    def test_file_records_packed_coupling(self):
        c = compile_interaction(FIXTURE_I.tables_dict, 1)
        geom = slot_geometry(so_bond_topology(c, CELL), (5, 1, 1))
        out = PackedArrays.zeros(sz_sector_map(geom, 2), (5,))
        with self.assertRaisesRegex(ValueError, exact(
                "_file_records: a PairLift record couples the S_z sectors +1 and -1 in the code "
                "slots; the packed layout needs an S_z-conserving interaction "
                "(select_layout)")):
            _file_records(out, c.records, geom)

    def test_file_records_packed_shape(self):
        out = PackedArrays.zeros(self.sm, (4,))
        with self.assertRaisesRegex(ValueError, exact(
                "_file_records: the packed out has leading shape (4,) and B = 5; the geometry "
                "needs (15,) and B = 5")):
            _file_records(out, _setup(self.fx)[0].records, self.geom)

    def test_gamma_p_mixed(self):
        with self.assertRaisesRegex(ValueError, exact(
                "gamma_p: gamma_b and out must both be PackedArrays or both dense arrays")):
            gamma_p(self.p, self.geom, out=np.zeros((15, 20, 20), complex))
        with self.assertRaisesRegex(ValueError, exact(
                "gamma_p: gamma_b and out must both be PackedArrays or both dense arrays")):
            gamma_p(self.p.to_dense(), self.geom, out=PackedArrays.zeros(self.sm, (15,)))

    def test_gamma_p_other_map(self):
        other = sz_sector_map(_setup(FIXTURE_LZSZ)[2], 4)
        with self.assertRaisesRegex(ValueError, exact(
                "gamma_p: the sector map of gamma_b (B = 3, nd = 4) is not the map of the "
                "geometry (B = 5)")):
            gamma_p(PackedArrays.zeros(other, (15,)), self.geom,
                    out=PackedArrays.zeros(other, (15,)))
        with self.assertRaisesRegex(ValueError, exact(
                "gamma_p: the sector map of out is not the sector map of gamma_b")):
            gamma_p(self.p, self.geom, out=PackedArrays.zeros(other, (15,)))

    def test_gamma_p_packed_shape_and_alias(self):
        with self.assertRaisesRegex(ValueError, exact(
                "gamma_p: the packed out has leading shape (3,); the geometry needs (15,)")):
            gamma_p(self.p, self.geom, out=PackedArrays.zeros(self.sm, (3,)))
        with self.assertRaisesRegex(ValueError, exact(
                "gamma_p: the packed gamma_b has leading shape (3,); the geometry needs (15,)")):
            gamma_p(PackedArrays.zeros(self.sm, (3,)), self.geom,
                    out=PackedArrays.zeros(self.sm, (15,)))
        with self.assertRaisesRegex(ValueError, exact("gamma_p: out shares memory with gamma_b")):
            gamma_p(self.p, self.geom, out=self.p)


# ---------------------------------------------------------------------------
# The packed bubble and the packed block store
# ---------------------------------------------------------------------------


def _noncanonical_green(fx):
    """The canonical G with ONE entry between the spin blocks set to 1e-14."""
    g = _green(fx)
    g[0, 1, 2, 0, fx.norb_phys] = 1e-14
    return g


def _leak(dense, lists):
    return float(np.abs(_off_sector(dense, lists)).max() / np.abs(dense).max())


class TestPackedBubble(unittest.TestCase):
    """Gate P-alg of the packed bubble and the packed store's readers: from
    a Green function built per spin block the dense production bubble has
    exactly zero off-sector entries, and its packed form (static sink and
    every bosonic slice in the packed store) equals the dense result
    restricted to the test-side sector lists (<= 1e-13 relative); the
    store's dense pair blocks and the self-energy transport on a packed
    ``W`` store equal the dense store's bitwise."""

    def test_static_packed_equals_dense_restricted(self):
        for fx in BUBBLE_NMAT:
            with self.subTest(fixture=NAMES[fx]):
                d, p = _bubbles(fx, "static")
                lists = _lists(fx)
                self.assertEqual(np.abs(_off_sector(d, lists)).max(), 0.0)
                self.assertIsInstance(p, PackedArrays)
                r = restrict(d, lists)
                for s in range(3):
                    self.assertGreater(np.abs(r[s]).max(), 0.0)
                    self.assertLessEqual(rel(p.blocks[s], r[s]), 1e-13)

    def test_every_slice_packed_equals_dense_restricted(self):
        for fx in BUBBLE_NMAT:
            with self.subTest(fixture=NAMES[fx]):
                (d, dpairs), (p, ppairs) = _bubbles(fx, "all")
                lists = _lists(fx)
                self.assertEqual(d.shape[0], BUBBLE_NMAT[fx])
                self.assertEqual(np.abs(_off_sector(d, lists)).max(), 0.0)
                self.assertEqual(p.lead_shape, (BUBBLE_NMAT[fx], fx.nvol))
                r = restrict(d, lists)
                for s in range(3):
                    self.assertGreater(np.abs(r[s]).max(), 0.0)
                    self.assertLessEqual(rel(p.blocks[s], r[s]), 1e-13)
                self.assertEqual(set(ppairs), set(dpairs))
                for key, blk in dpairs.items():
                    assert_recomputed(self, ppairs[key], blk, key)

    def test_self_energy_on_packed_w(self):
        for fx in BUBBLE_NMAT:
            with self.subTest(fixture=NAMES[fx]):
                geom, sm, nmat = _setup(fx)[2], _sm(fx), BUBBLE_NMAT[fx]
                rng = np.random.default_rng(31)
                w = PackedArrays(sm, tuple(
                    rng.normal(size=(nmat, fx.nvol, n, n))
                    + 1j * rng.normal(size=(nmat, fx.nvol, n, n)) for n in sm.sizes))
                g = _green(fx)
                sig = []
                with flex_bond.BondBlockStore(nmat, fx.nvol, sm.ND, fx.nd ** 2, ("W",)) as ds, \
                        PackedBondBlockStore(nmat, fx.nvol, sm, ("W",)) as ps:
                    ds.put_freq_batch("W", 0, nmat, w.to_dense())
                    ps.put_freq_batch("W", 0, nmat, w)
                    for st in (ds, ps):
                        sig.append(flex_bond.calc_self_energy_bond(
                            st, g, fx.beta, geom, fx.spatial_shape, fx.nd, 1))
                self.assertGreater(np.abs(sig[0]).max(), 0.0)
                assert_recomputed(self, sig[1], sig[0])

    def test_static_packed_peak_below_dense(self):
        """The traced peak of the static path above the bubble itself (the
        production bubble replaced by a copy of its result), sink included:
        packed (sink 6/16 of dense, gathered per momentum) stays below
        dense and within the dense bubble plus the packed sink plus 1/4
        (the gather through a dense round trip would add about 1.3)."""
        fx = FIXTURE_SZ_GRID
        topo, sm, g = _setup(fx)[1], _sm(fx), _green(fx)
        pre = _bubbles(fx, "static")[0].copy()
        D = pre.nbytes
        peaks = {}
        for layout in ("dense", "sz_packed"):
            plan = _plan(fx, layout, ("chi_static",), BUBBLE_NMAT[fx])
            with mock.patch.object(so_bond.bubble, "bond_bubble_static",
                                   side_effect=lambda *a, **k: pre.copy()):
                tracemalloc.start()
                try:
                    sink = (np.empty(pre.shape, dtype=np.complex128) if layout == "dense"
                            else PackedArrays.zeros(sm, (fx.nvol,)))
                    so_bond_bubble(g, None, fx.beta, topo, fx.spatial_shape, plan=plan,
                                   freq="static", sink=sink)
                    peaks[layout] = tracemalloc.get_traced_memory()[1]
                finally:
                    tracemalloc.stop()
            if layout == "sz_packed":
                r = restrict(pre, _lists(fx))
                self.assertTrue(all(np.array_equal(sink.blocks[s], r[s]) for s in range(3)))
        self.assertGreaterEqual(peaks["dense"], 2 * D)
        self.assertLess(peaks["sz_packed"], peaks["dense"])
        self.assertLessEqual(peaks["sz_packed"], (1 + 6 / 16 + 1 / 4) * D)


#: the refusal of a pair block with off-sector entries (``%s`` the pair and slot)
PAIR_REFUSAL = ("PackedBondBlockStore: the pair block ({}, {}) of '{}' has S_z-off-sector "
                "entries (max |entry| = {}): the Green function it was built from is not "
                "block-diagonal in spin; build it per spin sector (canonical_band_eigen) or "
                "project a borrowed one (project_borrowed_green)")
STATIC_REFUSAL = ("so_bond_bubble: the static bond bubble has S_z-off-sector entries (leak = "
                  "{}): the Green function is not block-diagonal in spin; build it per spin "
                  "sector (canonical_band_eigen) or project a borrowed one "
                  "(project_borrowed_green)")


class TestPackedStore(unittest.TestCase):
    """Construction check of the packed block store (protocol, sizes, views
    and copies, refusals), of ``open_store`` and of the packed bubble's
    refusals (a Green function that is not block-diagonal in spin, a sink of
    the wrong layout or map)."""

    fx = FIXTURE_LZSZ

    def _store(self, names=("chibar",), **kw):
        return PackedBondBlockStore(4, self.fx.nvol, _sm(self.fx), names, **kw)

    def test_protocol_lists(self):
        present = ("nmat", "nslices", "nvol", "ND", "nd", "B", "half_grid", "resident",
                   "layout", "sector_map", "slot_names", "nbytes", "released", "put_pair",
                   "get_pair", "get_pair_full", "put_freq_batch", "get_freq_batch",
                   "get_freq_batch_on", "packed", "grid_index", "detach", "release_slot",
                   "release", "to_host", "__enter__", "__exit__")
        absent = ("put_pair_fullgrid", "expand_full", "stager", "copy_workers")
        with self._store() as st:
            for name in present:
                self.assertTrue(hasattr(st, name), name)
            for name in absent:
                self.assertFalse(hasattr(st, name), name)
            sm = _sm(self.fx)
            self.assertEqual((st.nmat, st.nslices, st.nvol, st.ND, st.nd, st.B),
                             (4, 4, self.fx.nvol, sm.ND, self.fx.nd ** 2, sm.B))
            self.assertEqual(st.sector_map.nd, self.fx.nd)       # modes, not the pair size
            self.assertEqual((st.half_grid, st.resident, st.layout),
                             (False, "host", "sz_packed"))
            self.assertIs(st.sector_map, sm)
            self.assertEqual(st.slot_names, ("chibar",))
            self.assertTrue(np.array_equal(st.grid_index(1, 3), [1, 2]))
            self.assertIsNone(st.to_host())
        self.assertTrue(st.released)

    def test_nbytes(self):
        sm = _sm(self.fx)
        plan = _plan(self.fx, "sz_packed", ("chi_dynamic",), 4)
        with open_store(plan, sm) as st, \
                flex_bond.BondBlockStore(4, self.fx.nvol, sm.ND, self.fx.nd ** 2,
                                         ("chibar",)) as ds:
            self.assertEqual(16 * st.nbytes, 6 * ds.nbytes)
            self.assertEqual(st.nbytes, dict((n, b) for _, n, b in plan.rows)["chibar"])

    def test_half_grid_and_device_refused(self):
        sm = _sm(self.fx)
        with self.assertRaisesRegex(ValueError, exact(
                "PackedBondBlockStore: the half bosonic grid is not available on the "
                "S_z-packed store")):
            PackedBondBlockStore(4, 5, sm, ("chibar",), half_grid=True)
        fake = types.ModuleType("cupy")
        with self.assertRaisesRegex(ValueError, exact(
                "PackedBondBlockStore: the S_z-packed store is host only (xp must be numpy), "
                "got cupy")):
            PackedBondBlockStore(4, 5, sm, ("chibar",), xp=fake)
        with self._store() as st:
            for meth, args in (("get_pair_full", ("chibar", 0, 0)),
                               ("get_freq_batch_on", ("chibar", 0, 1))):
                with self.subTest(method=meth):
                    with self.assertRaisesRegex(ValueError, exact(
                            "PackedBondBlockStore.{}: the S_z-packed store is host only (xp "
                            "must be None or numpy), got cupy".format(meth))):
                        getattr(st, meth)(*args, xp=fake)
            self.assertTrue(np.shares_memory(st.get_freq_batch_on("chibar", 0, 1, np).blocks[0],
                                             st.packed("chibar").blocks[0]))

    def test_views_and_copies(self):
        sm = _sm(self.fx)
        rng = np.random.default_rng(5)
        x = PackedArrays(sm, tuple(rng.normal(size=(2, 5, n, n)) + 0j for n in sm.sizes))
        with self._store(("chibar", "chi"), aliases={"W": "chibar"}) as st:
            self.assertEqual(st.slot_names, ("chibar", "chi", "W"))
            st.put_freq_batch("chibar", 1, 3, x)
            v = st.get_freq_batch("W", 1, 3)
            self.assertTrue(all(np.array_equal(a, b) for a, b in zip(v.blocks, x.blocks)))
            v.blocks[1][0, 0, 0, 0] = 7.0                     # a view: reaches the store
            self.assertEqual(st.packed("chibar")[+1][1, 0, 0, 0], 7.0)
            dense = st.packed("chibar").to_dense()
            self.assertTrue(np.array_equal(st.get_pair("chibar", 1, 2),
                                           dense[:, :, 16:32, 32:48]))
            pair = st.get_pair("chibar", 0, 0)                 # a new array
            pair[...] = 0
            dense[...] = 0                                     # to_dense copies
            self.assertEqual(st.packed("chibar")[+1][1, 0, 0, 0], 7.0)
            with self.assertRaisesRegex(ValueError, exact(
                    "PackedBondBlockStore: cannot detach 'W', an alias of 'chibar'")):
                st.detach("W")
            got = st.detach("chi")
            self.assertIsInstance(got, PackedArrays)
            self.assertEqual(got.lead_shape, (4, 5))
            self.assertEqual(st.slot_names, ("chibar", "W"))
            with self.assertRaisesRegex(KeyError, "PackedBondBlockStore: slot 'chi' was released"):
                st.packed("chi")
            with self.assertRaisesRegex(KeyError, "PackedBondBlockStore: no slot 'nope'"):
                st.packed("nope")
            st.release_slot("chibar")
            st.release_slot("chibar")                          # tolerated
            with self.assertRaisesRegex(KeyError, "slot 'chibar' was released"):
                st.packed("W")
            self.assertEqual(st.nbytes, 0)
        with self.assertRaisesRegex(RuntimeError, exact("PackedBondBlockStore: released")):
            st.packed("chibar")

    def test_batch_refusals(self):
        sm = _sm(self.fx)
        other = _sm(FIXTURE_SZ_GRID)
        with self._store() as st:
            with self.assertRaisesRegex(ValueError, exact(
                    "PackedBondBlockStore.put_freq_batch: the batch of 'chibar' must be "
                    "PackedArrays of the store's sector map with leading shape (2, 5)")):
                st.put_freq_batch("chibar", 0, 2, np.zeros((2, 5, sm.ND, sm.ND), complex))
            for bad in (PackedArrays.zeros(other, (2, 5)), PackedArrays.zeros(sm, (3, 5))):
                with self.assertRaisesRegex(ValueError, "^PackedBondBlockStore.put_freq_batch"):
                    st.put_freq_batch("chibar", 0, 2, bad)
            with self.assertRaisesRegex(ValueError, exact(
                    "PackedBondBlockStore: the slice range [3, 5) must lie in [0, 4) and be "
                    "non-empty")):
                st.get_freq_batch("chibar", 3, 5)
            with self.assertRaisesRegex(ValueError, exact(
                    "PackedBondBlockStore.put_pair: the block of 'chibar' must have shape "
                    "(4, 5, 16, 16), got (4, 5, 16, 4)")):
                st.put_pair("chibar", 0, 0, np.zeros((4, 5, 16, 4), complex))
            with self.assertRaisesRegex(ValueError, exact(
                    "PackedBondBlockStore: channel indices must be integers in [0, 3), got "
                    "(0, 3)")):
                st.get_pair("chibar", 0, 3)

    def test_put_pair_off_sector_refused(self):
        with self._store() as st:
            blk = np.zeros((4, 5, 16, 16), dtype=np.complex128)
            blk[2, 1, 1, 0] = 1e-14          # pair 1 = (0, 1): sigma 0; pair 0 = (0, 0): 0
            st.put_pair("chibar", 1, 2, blk)                   # same sector: stored
            self.assertEqual(st.get_pair("chibar", 1, 2)[2, 1, 1, 0], 1e-14)
            blk[2, 1, 2, 0] = 1e-14          # pair 2 = (0, 2): up-down, sigma +1
            with self.assertRaisesRegex(ValueError, exact(PAIR_REFUSAL.format(
                    1, 2, "chibar", "1.000e-14"))):
                st.put_pair("chibar", 1, 2, blk)
            for bad, text in ((np.nan, "nan"), (np.inf, "inf")):
                blk[...] = 0
                blk[0, 0, 2, 0] = bad
                with self.subTest(entry=text):
                    with self.assertRaisesRegex(ValueError, exact(PAIR_REFUSAL.format(
                            0, 1, "chibar", text))):
                        st.put_pair("chibar", 0, 1, blk)
            self.assertEqual(np.abs(st.get_pair("chibar", 0, 1)).max(), 0.0)

    def test_noncanonical_green_refused(self):
        for fx in BUBBLE_NMAT:
            with self.subTest(fixture=NAMES[fx]):
                g = _noncanonical_green(fx)
                dense = _static(fx, "dense", g)
                leak = _leak(dense, _lists(fx))
                self.assertGreater(leak, 0.0)
                sink = PackedArrays.zeros(_sm(fx), (fx.nvol,))
                with self.assertRaisesRegex(ValueError, exact(STATIC_REFUSAL.format(
                        "%.3e" % leak))):
                    _static(fx, "sz_packed", g, sink)
                self.assertEqual(max(np.abs(b).max() for b in sink.blocks), 0.0)
                # the first pair block the bubble hands over is (0, 0); its
                # off-sector maximum read from the dense store of the same G
                blk = _every_slice(fx, "dense", g)[1][(0, 0)]
                psig = _sm(fx).sector_of[:fx.nd ** 2]
                off = np.abs(blk[..., psig[:, None] != psig[None, :]]).max()
                self.assertGreater(off, 0.0)
                with self.assertRaisesRegex(ValueError, exact(PAIR_REFUSAL.format(
                        0, 0, "chibar", "%.3e" % off))):
                    _every_slice(fx, "sz_packed", g)

    def test_static_sink_refusals(self):
        fx = FIXTURE_SZ_GRID
        topo = _setup(fx)[1]
        g = _green(fx)
        plan = _plan(fx, "sz_packed", ("chi_static",), BUBBLE_NMAT[fx])
        with self.assertRaisesRegex(ValueError, "^" + re.escape(
                "so_bond_bubble: the admission plan has layout = sz_packed but the call has "
                "layout = dense (")):
            so_bond_bubble(g, None, fx.beta, topo, fx.spatial_shape, plan=plan, freq="static",
                           sink=np.empty((15, 20, 20), complex))
        map_text = ("so_bond_bubble: the sector map of the packed {} is not the map of the "
                    "bond channels (B = 5, nd = 2, reverse = [0, 4, 3, 2, 1])")
        self.assertEqual(_setup(fx)[2].reverse.tolist(), [0, 4, 3, 2, 1])
        # another fixture's map, and a map of the same B and nd for another reverse
        other_rev = sz_sector_map(types.SimpleNamespace(reverse=np.arange(5)), 2)
        for bad in (PackedArrays.zeros(_sm(self.fx), (15,)), PackedArrays.zeros(other_rev, (15,))):
            with self.assertRaisesRegex(ValueError, exact(map_text.format("static sink"))):
                so_bond_bubble(g, None, fx.beta, topo, fx.spatial_shape, plan=plan,
                               freq="static", sink=bad)
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_bubble: the packed static sink must have leading shape (15,), got "
                "(3,)")):
            so_bond_bubble(g, None, fx.beta, topo, fx.spatial_shape, plan=plan, freq="static",
                           sink=PackedArrays.zeros(_sm(fx), (3,)))
        pdyn = _plan(fx, "sz_packed", ("chi_dynamic",), BUBBLE_NMAT[fx])
        with PackedBondBlockStore(8, 15, other_rev, ("chibar",)) as st:
            with self.assertRaisesRegex(ValueError, exact(map_text.format("store"))):
                so_bond_bubble(g, None, fx.beta, topo, fx.spatial_shape, plan=pdyn, freq="all",
                               sink=st)
        dplan = _plan(fx, "dense", ("chi_dynamic",), BUBBLE_NMAT[fx])
        with PackedBondBlockStore(8, 15, _sm(fx), ("chibar",)) as st:
            with self.assertRaisesRegex(ValueError, "^" + re.escape(
                    "so_bond_bubble: the admission plan has layout = dense but the call has "
                    "layout = sz_packed (")):
                so_bond_bubble(g, None, fx.beta, topo, fx.spatial_shape, plan=dplan, freq="all",
                               sink=st)

    def test_open_store(self):
        fx = self.fx
        sm = _sm(fx)
        dplan = _plan(fx, "dense", ("chi_dynamic",), 4)
        with open_store(dplan) as st:
            self.assertIs(type(st), flex_bond.BondBlockStore)
            self.assertEqual((st.nslices, st.nvol, st.ND, st.nd, st.slot_names),
                             (4, fx.nvol, sm.ND, fx.nd ** 2, ("chibar",)))
        with self.assertRaisesRegex(ValueError, exact(
                "open_store: a dense plan takes no sector map")):
            open_store(dplan, sm)
        pplan = _plan(fx, "sz_packed", ("chi_dynamic",), 4)
        with open_store(pplan, sm) as st:
            self.assertIs(type(st), PackedBondBlockStore)
            self.assertEqual((st.nslices, st.nvol, st.slot_names), (4, fx.nvol, ("chibar",)))
        with self.assertRaisesRegex(ValueError, exact(
                "open_store: an sz_packed plan needs the sector map (sz_sector_map)")):
            open_store(pplan)
        with self.assertRaisesRegex(ValueError, exact(
                "open_store: the sector map (B = 5, nd = 2) is not the plan's (B = 3, nd = 4)")):
            open_store(pplan, _sm(FIXTURE_SZ_GRID))


# ---------------------------------------------------------------------------
# The packed dressing and its conditioning guard
# ---------------------------------------------------------------------------

#: the (fixture, vertex scale) of the near-singular refusal input (a local
#: minimum of the guard score over the scale; score 4.3e-4)
NEAR_SINGULAR = (FIXTURE_SZ_GRID, -8.35)
#: the smallest relative gap of every round-off-sensitive choice (see _isolation)
ISOLATED = 1e-6


def _isolation(xd, gd, lists):
    """Test side: the relative gaps that make the guard's choices
    independent of round-off: ``(second-smallest block score / smallest
    - 1, the same for the sector sigma_min at the worst block, largest
    sector sigma_max / second largest - 1 there)``."""
    smin, smax = _sector_svd(xd, gd, lists)
    glob = np.sort(_scores(smin.min(0), smax.max(0)))
    i = int(np.argmin(_scores(smin.min(0), smax.max(0))))
    lo, hi = np.sort(smin[:, i]), np.sort(smax[:, i])
    return glob[1] / glob[0] - 1.0, lo[1] / lo[0] - 1.0, hi[-1] / hi[-2] - 1.0


def _sector_svd(xd, gd, lists):
    """Test side: per ph sector the extreme singular values of ``1 + X_s
    Gamma_s`` (``X_s``, ``Gamma_s`` the restrictions of the dense inputs),
    flattened over ``(l, q)``: ``(smin, smax)``, each ``(3, nb * nvol)``."""
    smin, smax = [], []
    for x, g in zip(restrict(xd, lists), restrict(gd, lists)):
        n = x.shape[-1]
        sv = np.linalg.svd((np.eye(n) + x @ g[None]).reshape(-1, n, n), compute_uv=False)
        smin.append(sv[:, -1])
        smax.append(sv[:, 0])
    smin, smax = np.array(smin), np.array(smax)
    return smin, smax


def _scores(smin, smax):
    """The guard score of one block: min(sigma_min/sigma_max, sigma_min/max(1, sigma_max))."""
    return np.minimum(smin / smax, smin / np.maximum(1.0, smax))


def _worst(xd, gd, lists):
    """Test side: ``(global score, per-sector minimum score at the worst
    block, name of the sigma_min sector, name of the sigma_max sector,
    largest sector condition number)``. The global score of a block uses
    the extreme singular values over all three sectors (the dense
    block's)."""
    smin, smax = _sector_svd(xd, gd, lists)
    glob = _scores(smin.min(0), smax.max(0))
    i = int(np.argmin(glob))
    names = so_bond.SECTOR_ORDER
    return (float(glob[i]), float(_scores(smin[:, i], smax[:, i]).min()),
            names[int(np.argmin(smin[:, i]))], names[int(np.argmax(smax[:, i]))],
            float((smax / smin).max()))


class TestPackedDress(unittest.TestCase):
    """Gate P-alg of the packed dressing: ``chi_s = [1 + X_s Gamma_s]^{-1}
    X_s`` per ph sector equals the dense production dressing restricted to
    the test-side sector lists (<= 1e-13 relative, for the static bubble
    and a batch of four bosonic slices), and the dense result has exactly
    zero off-sector entries. The conditioning guard runs over all sectors
    at once: its score equals the dense guard's (<= 1e-13 relative) and
    the reported sectors of sigma_min / sigma_max are the test-side argmin
    / argmax over the sector singular values. The comparison inputs are
    pre-asserted well conditioned (every sector condition number <= 1e3,
    guard score in [1e-2, 1]); the refusal uses a separate near-singular
    input and quotes the dense refusal followed by the two sectors. The
    bubbles carry a seeded sector-preserving perturbation, and every input
    is pre-asserted free of near ties (worst block, and the sectors of
    sigma_min / sigma_max there, each ahead of the runner-up by >= 1e-6
    relative), so no choice depends on round-off."""

    def test_packed_equals_dense_restricted(self):
        for fx in BUBBLE_NMAT:
            for kind in ("static", "batch"):
                with self.subTest(fixture=NAMES[fx], kind=kind):
                    xd, xp, gd, gp, l0 = _dress_inputs(fx, kind)
                    lists = _lists(fx)
                    glob, _, _, _, cond = _worst(xd, gd, lists)
                    self.assertGreaterEqual(min(_isolation(xd, gd, lists)), ISOLATED)
                    self.assertLessEqual(cond, 1e3)
                    self.assertTrue(1e-2 <= glob <= 1.0, glob)
                    chi_d, diag_d = _dress(fx, "dense", xd, gd, l0, kind)
                    chi_p, diag_p = _dress(fx, "sz_packed", xp, gp, l0, kind)
                    self.assertEqual(xd.shape[0], 1 if kind == "static" else 4)
                    # the dense dressing of block-diagonal inputs stays block diagonal
                    self.assertTrue(np.all(_off_sector(chi_d, lists) == 0))
                    self.assertIsInstance(chi_p, PackedArrays)
                    self.assertIs(chi_p.sector_map, xp.sector_map)
                    self.assertEqual(chi_p.lead_shape, xd.shape[:2])
                    for s, (got, want) in enumerate(zip(chi_p.blocks, restrict(chi_d, lists))):
                        self.assertGreater(np.abs(want).max(), 0.0)
                        self.assertLessEqual(rel(got, want), 1e-13, s)
                    # the dressing is not the identity on any sector
                    for got, x in zip(chi_p.blocks, xp.blocks):
                        self.assertGreater(rel(got, x), 1e-3)
                    self.assertEqual(set(diag_d), {"cond_min"})
                    self.assertEqual(set(diag_p),
                                     {"cond_min", "sigma_min_sector", "sigma_max_sector"})
                    self.assertLessEqual(abs(diag_p["cond_min"] - diag_d["cond_min"]),
                                         1e-13 * diag_d["cond_min"])
                    self.assertLessEqual(abs(diag_p["cond_min"] - glob), 1e-13 * glob)

    def test_named_sectors(self):
        for fx in BUBBLE_NMAT:
            for kind in ("static", "batch"):
                with self.subTest(fixture=NAMES[fx], kind=kind):
                    xd, xp, gd, gp, l0 = _dress_inputs(fx, kind)
                    glob, per, s_min, s_max, _ = _worst(xd, gd, _lists(fx))
                    self.assertGreaterEqual(min(_isolation(xd, gd, _lists(fx))), ISOLATED)
                    if fx is FIXTURE_SZ_GRID:
                        # non-vacuity: the two sectors differ, and so the
                        # global score is not the minimum of the sector scores
                        self.assertNotEqual(s_min, s_max)
                        self.assertGreaterEqual(abs(per - glob), 1e-3 * glob)
                    _, diag = _dress(fx, "sz_packed", xp, gp, l0, kind)
                    self.assertEqual((diag["sigma_min_sector"], diag["sigma_max_sector"]),
                                     (s_min, s_max))

    def test_guard_disabled(self):
        fx = FIXTURE_SZ_GRID
        xd, xp, gd, gp, l0 = _dress_inputs(fx, "static")
        chi_d, diag_d = _dress(fx, "dense", xd, gd, l0, "static", cond_tol=None)
        chi_p, diag_p = _dress(fx, "sz_packed", xp, gp, l0, "static", cond_tol=None)
        self.assertEqual(diag_d, {"cond_min": None})
        self.assertEqual(diag_p, {"cond_min": None, "sigma_min_sector": None,
                                  "sigma_max_sector": None})
        for got, want in zip(chi_p.blocks, restrict(chi_d, _lists(fx))):
            self.assertLessEqual(rel(got, want), 1e-13)

    def _refusals(self, fx, scale, cond_tol):
        xd, xp, gd, gp, l0 = _dress_inputs(fx, "static", scale)
        kw = {} if cond_tol is None else {"cond_tol": cond_tol}
        with self.assertRaises(bond_channels.BondConditioningError) as dense:
            _dress(fx, "dense", xd, gd, l0, "static", **kw)
        with self.assertRaises(bond_channels.BondConditioningError) as packed:
            _dress(fx, "sz_packed", xp, gp, l0, "static", **kw)
        self.assertGreaterEqual(min(_isolation(xd, gd, _lists(fx))), ISOLATED)
        d, p = dense.exception, packed.exception
        self.assertIsInstance(d.__cause__, bond_channels.BondConditioningError)
        self.assertIsInstance(p.__cause__, bond_channels.BondConditioningError)
        self.assertEqual(str(p.__cause__), str(d.__cause__))
        return _worst(xd, gd, _lists(fx)), d, p

    def test_refusal_names_the_sectors(self):
        fx, scale = NEAR_SINGULAR
        (glob, per, s_min, s_max, cond), d, p = self._refusals(fx, scale, None)
        # the input is near singular, and the guard's score is not the
        # minimum of the sector scores there
        self.assertLess(glob, 1e-3)
        self.assertNotEqual(s_min, s_max)
        self.assertGreaterEqual(abs(per - glob), 1e-3 * glob)
        self.assertEqual(str(p), str(d) + " (S_z sector of sigma_min: {}; of sigma_max: {})"
                         .format(s_min, s_max))
        for name in ("channel", "iq", "q", "l", "cond_tol"):
            self.assertEqual(getattr(p, name), getattr(d, name), name)
        # sigma_min carries an absolute error ~ u sigma_max, i.e. ~ u kappa
        # relative: the tolerance scales with the sector condition number
        # (kappa ~ 1.4e3 here)
        self.assertLessEqual(cond, 1e4)
        for name in ("worst", "ratio", "pole", "smin", "smax"):
            self.assertLessEqual(abs(getattr(p, name) - getattr(d, name)),
                                 1e-14 * cond * abs(getattr(d, name)), name)

    def test_refusal_one_sector(self):
        # a raised cond_tol refuses a block whose sigma_min and sigma_max lie
        # in the same sector: the second clause is omitted
        fx = FIXTURE_LZSZ
        (glob, _, s_min, s_max, _), d, p = self._refusals(fx, 1.0, 0.5)
        self.assertLess(glob, 0.5)
        self.assertEqual(s_min, s_max)
        self.assertEqual(str(p), str(d) + " (S_z sector of sigma_min: {})".format(s_min))


class TestPackedDressRefusals(unittest.TestCase):
    """Construction check: the refusals of the packed dressing."""

    def test_batch_above_capacity(self):
        fx = FIXTURE_SZ_GRID
        xp = _dress_inputs(fx, "batch")[1]
        plan = _dress_plan(fx, "sz_packed", "static")
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_dress: a batch of 4 slices exceeds the plan's dress_capacity = 1 "
                "({})".format(plan.reason))):
            so_bond.so_bond_dress(xp, _gamma(fx, "sz_packed"), plan=plan, l0=2,
                                  spatial_shape=fx.spatial_shape)

    def test_non_finite(self):
        fx = FIXTURE_SZ_GRID
        _, xp, _, gp, l0 = _dress_inputs(fx, "static")
        plan = _dress_plan(fx, "sz_packed", "static")
        for arg in ("X", "gamma_b"):
            for s in range(3):
                for bad in (np.nan, np.inf):
                    with self.subTest(arg=arg, sector=s, bad=bad):
                        x = _packed_lead(xp, np.s_[...])
                        g = _packed_lead(gp, np.s_[...])
                        (x if arg == "X" else g).blocks[s][..., 0, 0] = bad
                        for cond_tol in (None, 1e-12):
                            with self.assertRaisesRegex(NonFiniteError, exact(
                                    "so_bond_dress: {} has non-finite entries".format(arg))):
                                so_bond.so_bond_dress(x, g, plan=plan, l0=l0,
                                                      spatial_shape=fx.spatial_shape,
                                                      cond_tol=cond_tol)

    def test_mixed_layouts_and_maps(self):
        fx = FIXTURE_SZ_GRID
        xd, xp, gd, gp, l0 = _dress_inputs(fx, "static")
        pplan, dplan = _dress_plan(fx, "sz_packed", "static"), _dress_plan(fx, "dense", "static")
        both = "so_bond_dress: X and gamma_b must both be PackedArrays or both dense arrays"
        with self.assertRaisesRegex(ValueError, exact(both)):
            so_bond.so_bond_dress(xp, gd, plan=pplan, l0=l0, spatial_shape=fx.spatial_shape)
        with self.assertRaisesRegex(ValueError, exact(both)):
            so_bond.so_bond_dress(xd, gp, plan=dplan, l0=l0, spatial_shape=fx.spatial_shape)
        other = _random_packed(_sm(FIXTURE_LZSZ), fx.nvol, 3)
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_dress: the sector map of gamma_b is not the sector map of X")):
            so_bond.so_bond_dress(xp, other, plan=pplan, l0=l0, spatial_shape=fx.spatial_shape)
        short = _packed_lead(gp, np.s_[:2])
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_dress: the packed gamma_b has leading shape (2,); X needs ({},)"
                .format(fx.nvol))):
            so_bond.so_bond_dress(xp, short, plan=pplan, l0=l0, spatial_shape=fx.spatial_shape)
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_dress: prod(spatial_shape) != nvol")):
            so_bond.so_bond_dress(xp, gp, plan=pplan, l0=l0, spatial_shape=(1, 1, 1))


# ---------------------------------------------------------------------------
# The pairing vertex: static parts, the frequency batch, the FLEX W map
# ---------------------------------------------------------------------------

#: the fixtures of the pairing comparisons
PAIRING_FIXTURES = (FIXTURE_SZ_GRID, FIXTURE_LZSZ)


def _assert_restricted(tc, packed, dense, lists, tol=1e-13):
    """``packed`` equals ``restrict(dense)`` (each sector nonzero) and the
    dense array is exactly zero off the sectors."""
    tc.assertIsInstance(packed, PackedArrays)
    tc.assertEqual(np.abs(_off_sector(dense, lists)).max(), 0.0)
    for s, (got, want) in enumerate(zip(packed.blocks, restrict(dense, lists))):
        tc.assertGreater(np.abs(want).max(), 0.0, s)
        tc.assertLessEqual(rel(got, want), tol, s)


def _packed_equal(a, b):
    return all(np.array_equal(x, y) for x, y in zip(a.blocks, b.blocks))


class TestPackedParts(unittest.TestCase):
    """Gate P-alg of the static pairing parts: ``V_bare = 1/4 (Gamma +
    Gamma^P)`` and ``V_fluct = -Gamma chi Gamma`` computed sector by sector
    from the packed vertex and the packed dressed static chi equal the
    dense production parts restricted to the test-side sector lists
    (<= 1e-13 relative, every sector nonzero), and the dense parts are
    exactly zero off the sectors."""

    def test_packed_equals_dense_restricted(self):
        for fx in PAIRING_FIXTURES:
            with self.subTest(fixture=NAMES[fx]):
                _, _, gd, gp, chi_d, chi_p = _dressed(fx, "static")
                g0, c0 = _copy_packed(gp), _copy_packed(chi_p)
                vb_d, vf_d = _parts(fx, "dense", gd, chi_d)
                vb_p, vf_p = _parts(fx, "sz_packed", gp, chi_p)
                lists = _lists(fx)
                for v in (vb_p, vf_p):
                    self.assertIs(v.sector_map, gp.sector_map)
                    self.assertEqual(v.lead_shape, (fx.nvol,))
                _assert_restricted(self, vb_p, vb_d, lists)
                _assert_restricted(self, vf_p, vf_d, lists)
                # the inputs are not written
                self.assertTrue(_packed_equal(gp, g0) and _packed_equal(chi_p, c0))


class TestPackedPartsRefusals(unittest.TestCase):
    """Construction check: the refusals of the packed pairing parts."""

    @classmethod
    def setUpClass(cls):
        cls.fx = fx = FIXTURE_SZ_GRID
        _, _, cls.gd, cls.gp, cls.cd, cls.cp = _dressed(fx, "static")
        cls.geom = _setup(fx)[2]
        cls.plan = _parts_plan(fx, "sz_packed")

    def _call(self, g, c, plan=None):
        return so_bond.so_bond_pairing_vertex_parts(g, c, self.geom,
                                                    plan=self.plan if plan is None else plan)

    def test_mixed_layouts(self):
        msg = ("so_bond_pairing_vertex_parts: gamma_b and chi_b must both be PackedArrays or "
               "both dense arrays")
        for g, c in ((self.gp, self.cd), (self.gd, self.cp)):
            with self.assertRaisesRegex(ValueError, exact(msg)):
                self._call(g, c)

    def test_packed_with_dense_plan(self):
        plan = _parts_plan(self.fx, "dense")
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_pairing_vertex_parts: the admission plan has layout = dense but the "
                "call has layout = sz_packed ({})".format(plan.reason))):
            self._call(self.gp, self.cp, plan)
        # and a dense call with a packed plan
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_pairing_vertex_parts: the admission plan has layout = sz_packed but "
                "the call has layout = dense ({})".format(self.plan.reason))):
            self._call(self.gd, self.cd)

    def test_maps_and_shapes(self):
        fx = self.fx
        other = _random_packed(_sm(FIXTURE_LZSZ), fx.nvol, 5)
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_pairing_vertex_parts: the sector map of chi_b is not the sector map "
                "of gamma_b")):
            self._call(self.gp, other)
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_pairing_vertex_parts: the packed chi_b has leading shape (2,); "
                "gamma_b has ({},)".format(fx.nvol))):
            self._call(self.gp, _packed_lead(self.cp, np.s_[:2]))
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_pairing_vertex_parts: the packed gamma_b must have leading shape "
                "(nvol,), got (1, {})".format(fx.nvol))):
            self._call(_packed_lead(self.gp, np.s_[None]), self.cp)
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_pairing_vertex_parts: chi_b shares memory with gamma_b")):
            self._call(self.gp, self.gp)

    def test_non_finite_and_exchange(self):
        for arg in ("gamma_b", "chi_b"):
            for s in range(3):
                with self.subTest(arg=arg, sector=s):
                    g, c = _copy_packed(self.gp), _copy_packed(self.cp)
                    (g if arg == "gamma_b" else c).blocks[s][0, 0, 0] = np.nan
                    with self.assertRaisesRegex(ValueError, exact(
                            "so_bond_pairing_vertex_parts: {} has non-finite entries"
                            .format(arg))):
                        self._call(g, c)
        # a vertex that is not its exchange image (only one sector broken)
        g = _copy_packed(self.gp)
        g.blocks[1][...] += 0.5 * np.abs(g.blocks[1]).max() * np.eye(g.blocks[1].shape[-1])[
            None] * np.arange(1, self.fx.nvol + 1)[:, None, None]
        with self.assertRaisesRegex(ValueError, r"^so_bond_pairing_vertex_parts: Gamma_B "
                                    r"violates the particle-exchange symmetry \(\|Gamma\^P - "
                                    r"Gamma\| = \S+ > 1e-13 \|Gamma\|\)$"):
            self._call(g, self.cp)


class TestVertexBatch(unittest.TestCase):
    """Gate P-alg of the per-frequency pairing vertex batch ``V[l] = 1/4
    (Gamma + Gamma^P) - Gamma Y[l] Gamma`` for a batch of four dressed
    bosonic slices (``so_bond_dress`` output): the dense batch equals the
    sum of the static parts of each slice (<= 1e-13 relative), and the
    packed batch equals the dense batch restricted to the test-side sector
    lists (<= 1e-13 relative), whose off-sector entries are exactly zero.
    Neither ``gamma_b`` nor ``Y`` is written."""

    def test_dense_batch_equals_parts(self):
        for fx in PAIRING_FIXTURES:
            with self.subTest(fixture=NAMES[fx]):
                _, _, gd, _, Y, _ = _dressed(fx, "batch")
                self.assertEqual(Y.shape[0], 4)
                g0, y0 = gd.copy(), Y.copy()
                V = _batch(fx, "dense", gd, Y)
                self.assertTrue(np.array_equal(gd, g0) and np.array_equal(Y, y0))
                for l in range(Y.shape[0]):
                    vb, vf = _parts(fx, "dense", gd, Y[l].copy())
                    self.assertGreater(rel(vf, vb), 1e-3)          # both terms count
                    self.assertLessEqual(rel(V[l], vb + vf), 1e-13, l)
                # the slices differ: the batch is not one slice repeated
                self.assertGreater(rel(V[0], V[3]), 1e-3)

    def test_packed_batch_equals_dense_restricted(self):
        for fx in PAIRING_FIXTURES:
            with self.subTest(fixture=NAMES[fx]):
                _, _, gd, gp, Yd, Yp = _dressed(fx, "batch")
                g0, y0 = _copy_packed(gp), _copy_packed(Yp)
                Vd = _batch(fx, "dense", gd, Yd)
                Vp = _batch(fx, "sz_packed", gp, Yp)
                self.assertTrue(_packed_equal(gp, g0) and _packed_equal(Yp, y0))
                self.assertEqual(Vp.lead_shape, Yp.lead_shape)
                _assert_restricted(self, Vp, Vd, _lists(fx))

    def test_out_overwritten(self):
        """``out`` holding garbage is overwritten completely."""
        fx = FIXTURE_SZ_GRID
        _, _, gd, gp, Yd, Yp = _dressed(fx, "batch")
        ref = _batch(fx, "sz_packed", gp, Yp)
        out = _random_packed(Yp.sector_map, fx.nvol, 9)
        out = PackedArrays(out.sector_map, tuple(np.repeat(b[None], 4, axis=0) * 1e6
                                                 for b in out.blocks))
        assert_recomputed(self, _batch(fx, "sz_packed", gp, Yp, out=out), ref)
        dense = np.full(Yd.shape, 1e6 + 1e6j)
        assert_recomputed(self, _batch(fx, "dense", gd, Yd, out=dense),
                          _batch(fx, "dense", gd, Yd))


class TestVertexBatchRefusals(unittest.TestCase):
    """Construction check: the refusals of the pairing vertex batch."""

    @classmethod
    def setUpClass(cls):
        cls.fx = fx = FIXTURE_SZ_GRID
        _, _, cls.gd, cls.gp, cls.Yd, cls.Yp = _dressed(fx, "batch")
        cls.geom = _setup(fx)[2]

    def _call(self, g, Y, out, plan):
        return so_bond.so_bond_pairing_vertex_batch(g, Y, self.geom, plan=plan, out=out)

    def test_out_aliases_Y(self):
        fx = self.fx
        msg = "so_bond_pairing_vertex_batch: out shares memory with Y"
        Yd, Yp = self.Yd.copy(), _copy_packed(self.Yp)
        with self.assertRaisesRegex(ValueError, exact(msg)):
            self._call(self.gd, Yd, Yd, _batch_plan(fx, "dense"))
        with self.assertRaisesRegex(ValueError, exact(msg)):
            self._call(self.gp, Yp, Yp, _batch_plan(fx, "sz_packed"))
        # one shared sector block suffices
        out = PackedArrays(Yp.sector_map, (np.zeros_like(Yp.blocks[0]), Yp.blocks[1],
                                           np.zeros_like(Yp.blocks[2])))
        with self.assertRaisesRegex(ValueError, exact(msg)):
            self._call(self.gp, Yp, out, _batch_plan(fx, "sz_packed"))
        self.assertTrue(np.array_equal(Yd, self.Yd) and _packed_equal(Yp, self.Yp))

    def test_plan_refusals(self):
        fx = self.fx
        out = np.empty_like(self.Yd)
        plan = _plan(fx, "dense", ("chi_dynamic",), BUBBLE_NMAT[fx])
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_pairing_vertex_batch: the admission plan has no row "
                "'vertex_batch_workspace' (closed consumers: vertex, chi_dynamic; {})"
                .format(plan.reason))):
            self._call(self.gd, self.Yd, out, plan)
        plan = _batch_plan(fx, "dense", nb=2)
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_pairing_vertex_batch: a batch of 4 slices exceeds the plan's batch "
                "capacity (dress_capacity = 2) ({})".format(plan.reason))):
            self._call(self.gd, self.Yd, out, plan)
        plan = _batch_plan(fx, "dense")
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_pairing_vertex_batch: the admission plan has layout = dense but the "
                "call has layout = sz_packed ({})".format(plan.reason))):
            self._call(self.gp, self.Yp, PackedArrays.zeros(self.Yp.sector_map,
                                                            self.Yp.lead_shape), plan)

    def test_argument_refusals(self):
        fx = self.fx
        dplan, pplan = _batch_plan(fx, "dense"), _batch_plan(fx, "sz_packed")
        name = "so_bond_pairing_vertex_batch: "
        with self.assertRaisesRegex(ValueError, exact(
                name + "gamma_b, Y and out must all be PackedArrays or all dense arrays")):
            self._call(self.gp, self.Yp, np.empty_like(self.Yd), pplan)
        with self.assertRaisesRegex(ValueError, exact(
                name + "out must be a complex128 array of shape {}, got complex128 {}".format(
                    self.Yd.shape, self.Yd[:2].shape))):
            self._call(self.gd, self.Yd, np.empty_like(self.Yd[:2]), dplan)
        with self.assertRaisesRegex(ValueError, exact(
                name + "Y must be a complex128 (nb, nvol, ND, ND) array matching gamma_b {}, "
                "got complex128 {}".format(self.gd.shape, self.Yd[0].shape))):
            self._call(self.gd, self.Yd[0], np.empty_like(self.Yd[0]), dplan)
        with self.assertRaisesRegex(ValueError, exact(
                name + "the packed out has leading shape (2, {}); Y has (4, {})".format(
                    fx.nvol, fx.nvol))):
            self._call(self.gp, self.Yp, PackedArrays.zeros(self.Yp.sector_map, (2, fx.nvol)),
                       pplan)
        other = _random_packed(_sm(FIXTURE_LZSZ), fx.nvol, 4)
        other = PackedArrays(other.sector_map, tuple(np.repeat(b[None], 4, axis=0)
                                                     for b in other.blocks))
        with self.assertRaisesRegex(ValueError, exact(
                name + "the sector map of Y is not the sector map of gamma_b")):
            self._call(self.gp, other, PackedArrays.zeros(self.Yp.sector_map,
                                                          self.Yp.lead_shape), pplan)
        Y = _copy_packed(self.Yp)
        Y.blocks[2][1, 0, 0, 0] = np.inf
        with self.assertRaisesRegex(ValueError, exact(name + "Y has non-finite entries")):
            self._call(self.gp, Y, PackedArrays.zeros(Y.sector_map, Y.lead_shape), pplan)
        Y = self.Yd.copy()
        Y[3, 0, 0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, exact(name + "Y has non-finite entries")):
            self._call(self.gd, Y, np.empty_like(Y), dplan)

    def test_more_argument_refusals(self):
        fx = self.fx
        dplan, pplan = _batch_plan(fx, "dense"), _batch_plan(fx, "sz_packed")
        name = "so_bond_pairing_vertex_batch: "
        sm = self.Yp.sector_map
        # out shares memory with gamma_b (dense and packed)
        out = np.zeros_like(self.Yd)
        out[0] = self.gd
        with self.assertRaisesRegex(ValueError, exact(name + "out shares memory with gamma_b")):
            self._call(out[0], self.Yd, out, dplan)
        pout = PackedArrays.zeros(sm, self.Yp.lead_shape)
        g = PackedArrays(sm, tuple(b[0] for b in pout.blocks))
        for b, src in zip(g.blocks, self.gp.blocks):
            b[...] = src
        with self.assertRaisesRegex(ValueError, exact(name + "out shares memory with gamma_b")):
            self._call(g, self.Yp, pout, pplan)
        # a non-finite packed gamma_b
        g = _copy_packed(self.gp)
        g.blocks[1][3, 0, 0] = np.inf
        with self.assertRaisesRegex(ValueError, exact(name + "gamma_b has non-finite entries")):
            self._call(g, self.Yp, PackedArrays.zeros(sm, self.Yp.lead_shape), pplan)
        # the sector map of out, and that of gamma_b against the geometry
        other = sz_sector_map(types.SimpleNamespace(reverse=np.arange(self.geom.B)), fx.nd)
        self.assertEqual(other.ND, sm.ND)
        with self.assertRaisesRegex(ValueError, exact(
                name + "the sector map of out is not the sector map of gamma_b")):
            self._call(self.gp, self.Yp, PackedArrays.zeros(other, self.Yp.lead_shape), pplan)
        with self.assertRaisesRegex(ValueError, exact(
                name + "the sector map of gamma_b (B = 5, nd = 2) is not the map of the "
                "geometry (B = 5)")):
            self._call(PackedArrays(other, tuple(b.copy() for b in self.gp.blocks)), self.Yp,
                       PackedArrays.zeros(sm, self.Yp.lead_shape), pplan)
        # the leading shape of a packed Y; the dtype of a dense gamma_b
        Y1 = _packed_lead(self.Yp, 0)
        with self.assertRaisesRegex(ValueError, exact(
                name + "the packed Y has leading shape ({},); gamma_b needs (nb, {})".format(
                    fx.nvol, fx.nvol))):
            self._call(self.gp, Y1, PackedArrays.zeros(sm, Y1.lead_shape), pplan)
        with self.assertRaisesRegex(ValueError, exact(
                name + "gamma_b must be a complex128 (nvol, ND, ND) array")):
            self._call(self.gd.real.copy(), self.Yd, np.empty_like(self.Yd), dplan)

    def test_exchange_refusal(self):
        """``Gamma_B`` that is not its exchange image is refused as by the
        static parts; ``out[0]`` then holds ``Gamma^P`` and the other slices
        are untouched."""
        fx = self.fx
        regex = (r"^so_bond_pairing_vertex_batch: Gamma_B violates the particle-exchange "
                 r"symmetry \(\|Gamma\^P - Gamma\| = \S+ > 1e-13 \|Gamma\|\)$")
        g = _copy_packed(self.gp)
        g.blocks[1][...] += 0.5 * np.abs(g.blocks[1]).max() * np.eye(g.blocks[1].shape[-1])[
            None] * np.arange(1, fx.nvol + 1)[:, None, None]
        gd = g.to_dense()
        for layout, gam, Y in (("sz_packed", g, self.Yp), ("dense", gd, self.Yd)):
            with self.subTest(layout=layout):
                if layout == "dense":
                    out = np.full(Y.shape, 7.0 + 7.0j)
                else:
                    out = PackedArrays(Y.sector_map, tuple(np.full(b.shape, 7.0 + 7.0j)
                                                           for b in Y.blocks))
                with self.assertRaisesRegex(ValueError, regex):
                    self._call(gam, Y, out, _batch_plan(fx, layout))
                if layout == "dense":
                    head, rest = out[0], out[1:]
                    img = gamma_p(gd, self.geom, out=np.empty_like(gd))
                else:
                    head = PackedArrays(out.sector_map, tuple(b[0] for b in out.blocks))
                    rest = np.concatenate([b[1:].ravel() for b in out.blocks])
                    img = gamma_p(g, self.geom,
                                  out=PackedArrays.zeros(out.sector_map, (fx.nvol,)))
                self.assertGreater(rel(img.to_dense() if layout == "sz_packed" else img, gd),
                                   1e-3)
                assert_recomputed(self, head, img)
                self.assertTrue(np.all(np.asarray(rest) == 7.0 + 7.0j))


class TestWCrossMap(unittest.TestCase):
    """Gate P-alg of the FLEX effective interaction in the packed layout:
    ``M = Gamma (Y - X/2) Gamma`` of a dressed batch of four slices
    (test-side products per sector and dense) has +1 and -1 blocks that
    differ by >= 1e-2 of its largest entry (pre-asserted, so reading a
    block from the wrong sector is visible); ``p_b`` of the packed ``M``
    equals ``p_b`` of the dense ``M`` restricted to the test-side sector
    lists (<= 1e-13 relative), the dense ``W`` is exactly zero off the
    sectors, and the same-sector reading misses by >= 1e-3."""

    def test_w_cross_sector(self):
        fx = FIXTURE_LZSZ
        xd, xp, gd, gp, Yd, Yp = _dressed(fx, "batch")
        Md = gd[None] @ (Yd - 0.5 * xd) @ gd[None]
        Mp = PackedArrays(xp.sector_map, tuple(
            g[None] @ (y - 0.5 * x) @ g[None] for g, y, x in zip(gp.blocks, Yp.blocks,
                                                               xp.blocks)))
        scale = max(np.abs(b).max() for b in Mp.blocks)
        self.assertGreaterEqual(np.abs(Mp[+1] - Mp[-1]).max(), 1e-2 * scale)
        Wd = so_bond.p_b(Md, fx.nd)
        Wp = so_bond.p_b(Mp, fx.nd)
        lists = _lists(fx)
        _assert_restricted(self, Wp, Wd, lists)
        # the same-sector reading (block sigma permuted from block sigma) misses
        sm = Mp.sector_map
        for s, cross in enumerate(sm.pb_cross):
            if cross.source == sm.sigmas[s]:
                continue
            perm = cross.perm
            wrong = np.take(np.take(Mp.blocks[s], perm, axis=-2), perm, axis=-1)
            self.assertGreaterEqual(rel(wrong, restrict(Wd, lists)[s]), 1e-3, s)


if __name__ == "__main__":
    unittest.main()
