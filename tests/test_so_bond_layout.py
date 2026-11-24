"""Layout selection of the spin-orbital bond objects and the geometry view.

* ``so_bond.parse_layout_key`` reads the ``so_bond_layout`` key: ``"auto"``
  or ``"dense"`` (case-insensitive); ``"sz_packed"`` cannot be forced.
* ``so_bond.band_conserves_sz`` decides from the Transfer table alone
  whether the band conserves S_z (every spin-off-diagonal entry exactly
  zero); Extern (``H1 (x) sigma_z``) never breaks it.
* ``so_bond.records_conserve_sz`` names the first compiled record with
  ``S_z(1) - S_z(2) + S_z(3) - S_z(4) != 0``.
* ``so_bond.select_layout`` combines the override, the band and the records
  into ``(layout, reason)`` with exact reason texts.
* ``BondSlotGeometry.n_channels`` lets the geometry feed the bond pairing
  kernel and the bond self-energy transport as a view; the public helpers
  equal their internal geometry-first adapters.
"""
import dataclasses
import inspect
import unittest

import numpy as np

from hwave.solver import bond_channels, flex_bond, so_bond
from hwave.solver.so_bond import (BandSzCheck, SOBondCaps, band_conserves_sz, parse_layout_key,
                                  plan_so_bond, records_conserve_sz, select_layout,
                                  slot_geometry, so_bond_bubble, so_bond_dress, so_bond_topology,
                                  so_bond_vertex)
from hwave.solver.so_records import compile_interaction
from tests.so_bond_fixtures import assert_recomputed, exact, pair_kernel, quiet
from tests.so_bond_oracles import (BASE, CELL, CONVENTION_DR, CONVENTION_REV, CONVENTION_SHAPE,
                                   closed_row, random_green)

K = "so_bond_layout"
PACKED_INFO = ("so bond layout: sz_packed (S_z conserved: no spin-flip transfer, "
               "no S_z-changing interaction)")
FORCED_INFO = 'so bond layout: dense (reason: forced by so_bond_layout = "dense")'
SZ_PACKED_MSG = ('so_bond_layout = "sz_packed" cannot be forced: the S_z-packed layout is '
                 'chosen automatically when the input conserves S_z (no spin-flip Transfer '
                 'entry, no S_z-changing interaction). Use "auto" (the default) or "dense".')


def other_msg(value):
    return 'so_bond_layout must be "auto" or "dense", got {!r}'.format(value)


def _spin_diagonal_transfer(h):
    """norb_phys = 1 chain: hopping -1 on both spins, Zeeman ``+-h`` on site
    (interleaved indices ``2 orb + spin``)."""
    tr = {}
    for s, sign in ((0, 1.0), (1, -1.0)):
        tr[((0, 0, 0), (s, s))] = sign * h
        tr[((1, 0, 0), (s, s))] = -1.0
        tr[((-1, 0, 0), (s, s))] = -1.0
    return tr


#: Zeeman field through the Transfer table
ZEEMAN_TRANSFER = _spin_diagonal_transfer(0.15)
#: no field in the Transfer table (Zeeman through Extern)
PLAIN_TRANSFER = _spin_diagonal_transfer(0.0)
#: an Extern table in physical orbital indices (``coeff_extern H1 (x) sigma_z``)
EXTERN_Z = {((0, 0, 0), (0, 0)): 0.2}

_INTRA = {"CoulombIntra": {((0, 0, 0), (0, 0)): BASE["CoulombIntra"]}}
#: every S_z-conserving type, one table each (norb_phys = 1)
_CONSERVING_TYPES = {t: {t: closed_row(1, 0, 0, BASE[t])}
                     for t in ("CoulombInter", "Hund", "Ising", "Exchange", "PairHop")}
_CONSERVING_TYPES["CoulombIntra"] = _INTRA
ALL_CONSERVING = {k: v for tbl in _CONSERVING_TYPES.values() for k, v in tbl.items()}
WITH_PAIRLIFT = {**_INTRA, "PairLift": closed_row(1, 0, 0, BASE["PairLift"])}
ZERO_PAIRLIFT = {**_INTRA, "PairLift": closed_row(1, 0, 0, 0.0)}
COMPILED_SZ = compile_interaction(ALL_CONSERVING, 1)

CONSERVING_CASES = {}
for _t, _tables in _CONSERVING_TYPES.items():
    CONSERVING_CASES["transfer zeeman, " + _t] = (ZEEMAN_TRANSFER, None, _tables)
    CONSERVING_CASES["extern zeeman, " + _t] = (PLAIN_TRANSFER, EXTERN_Z, _tables)
CONSERVING_CASES["transfer zeeman, every type"] = (ZEEMAN_TRANSFER, None, ALL_CONSERVING)
CONSERVING_CASES["extern zeeman, every type"] = (PLAIN_TRANSFER, EXTERN_Z, ALL_CONSERVING)


class TestSelection(unittest.TestCase):
    """Gate P-sel: the layout is selected exactly from the declared inputs."""

    def test_conserving_inputs_select_packed(self):
        for name, (transfer, extern, tables) in CONSERVING_CASES.items():
            with self.subTest(name):
                lay, why = select_layout(compile_interaction(tables, 1), transfer,
                                         override="auto")
                self.assertEqual((lay, why), ("sz_packed", PACKED_INFO))
                self.assertEqual(band_conserves_sz(transfer, extern, 1), BandSzCheck(True, None))
                self.assertIsNone(records_conserve_sz(compile_interaction(tables, 1)))

    def test_tiny_soc_selects_dense(self):
        tr = {**ZEEMAN_TRANSFER, ((1, 0, 0), (0, 1)): 1e-14}
        self.assertEqual(
            select_layout(COMPILED_SZ, tr, override="auto"),
            ("dense", "so bond layout: dense (reason: spin-flip transfer entry R=(1, 0, 0) "
                      "(1, 2) = (1e-14+0j))"))
        self.assertEqual(band_conserves_sz(tr, None, 1),
                         BandSzCheck(False, ((1, 0, 0), (0, 1), 1e-14 + 0j)))

    def test_exact_zero_spin_flip_entry_conserves(self):
        tr = {**ZEEMAN_TRANSFER, ((1, 0, 0), (0, 1)): 0j, ((0, 0, 0), (1, 0)): -0.0}
        self.assertEqual(select_layout(COMPILED_SZ, tr, override="auto")[0], "sz_packed")

    def test_pairlift(self):
        self.assertEqual(select_layout(compile_interaction(WITH_PAIRLIFT, 1), ZEEMAN_TRANSFER,
                                       override="auto"),
                         ("dense", "so bond layout: dense (reason: PairLift declared)"))
        self.assertEqual(select_layout(compile_interaction(ZERO_PAIRLIFT, 1), ZEEMAN_TRANSFER,
                                       override="auto"),
                         ("sz_packed", PACKED_INFO))
        rec = records_conserve_sz(compile_interaction(WITH_PAIRLIFT, 1))
        self.assertEqual(str(rec["type"]), "PairLift")

    def test_fallback_record_reason(self):
        """A non-PairLift S_z-changing record (not produced by the compiler)
        is named by its type."""
        recs = COMPILED_SZ.records.copy()
        i = int(np.flatnonzero(recs["type"] == "Exchange")[0])
        recs["modes"][i] = (0, 1, 0, 1)                     # S_z: +1/2 - (-1/2) twice
        recs.flags.writeable = False
        bad = dataclasses.replace(COMPILED_SZ, records=recs)
        self.assertEqual(
            select_layout(bad, ZEEMAN_TRANSFER, override="auto"),
            ("dense", "so bond layout: dense (reason: S_z-changing interaction record of "
                      "type Exchange)"))

    def test_precedence(self):
        """Forced dense, then the Transfer entry, then the record."""
        tr = {**ZEEMAN_TRANSFER, ((0, 0, 0), (0, 1)): 0.1, ((0, 0, 0), (1, 0)): 0.1}
        pl = compile_interaction(WITH_PAIRLIFT, 1)
        self.assertEqual(select_layout(pl, tr, override="dense"), ("dense", FORCED_INFO))
        self.assertEqual(select_layout(pl, tr, override="auto")[1],
                         "so bond layout: dense (reason: spin-flip transfer entry R=(0, 0, 0) "
                         "(1, 2) = (0.1+0j))")

    def test_canonical_predicate_is_independent(self):
        # PairLift on a spin-diagonal band: layout dense, the band still conserves S_z
        self.assertEqual(select_layout(compile_interaction(WITH_PAIRLIFT, 1), ZEEMAN_TRANSFER,
                                       override="auto")[0], "dense")
        self.assertEqual(band_conserves_sz(ZEEMAN_TRANSFER, None, 1), BandSzCheck(True, None))
        # forced dense on a conserving input: layout dense, the band conserves S_z
        self.assertEqual(select_layout(COMPILED_SZ, ZEEMAN_TRANSFER, override="dense"),
                         ("dense", FORCED_INFO))
        self.assertTrue(band_conserves_sz(ZEEMAN_TRANSFER, EXTERN_Z, 1).conserving)
        self.assertEqual(list(inspect.signature(band_conserves_sz).parameters),
                         ["transfer", "extern", "norb_phys"])

    def test_key_values(self):
        self.assertEqual([parse_layout_key(v, key=K) for v in ("auto", "AUTO", "Dense")],
                         ["auto", "auto", "dense"])
        for v, msg in (("sz_packed", SZ_PACKED_MSG), ("SZ_Packed", SZ_PACKED_MSG),
                       ("sparse", other_msg("sparse")), (1, other_msg(1)),
                       (None, other_msg(None))):
            with self.subTest(v=v):
                with self.assertRaisesRegex(ValueError, exact(msg)):
                    parse_layout_key(v, key=K)
        with self.assertRaisesRegex(ValueError, exact(SZ_PACKED_MSG)):
            select_layout(COMPILED_SZ, ZEEMAN_TRANSFER, override="sz_packed")

    def test_compiled_checked_in_every_mode(self):
        for override in ("auto", "dense"):
            with self.subTest(override=override), self.assertRaisesRegex(ValueError, exact(
                    "select_layout: compiled must be a CompiledInteraction "
                    "(so_records.compile_interaction), got dict")):
                select_layout({}, ZEEMAN_TRANSFER, override=override)

    def test_first_offending_entry_order(self):
        tr = {**PLAIN_TRANSFER, ((1, 0, 0), (0, 1)): 0.3, ((0, 0, 0), (1, 0)): 0.2 + 0.1j,
              ((-1, 0, 0), (1, 0)): 0.3}
        self.assertEqual(band_conserves_sz(tr, None, 1),
                         BandSzCheck(False, ((-1, 0, 0), (1, 0), 0.3 + 0j)))
        self.assertEqual(select_layout(COMPILED_SZ, tr, override="auto")[1],
                         "so bond layout: dense (reason: spin-flip transfer entry "
                         "R=(-1, 0, 0) (2, 1) = (0.3+0j))")

    def test_norb2_interleaved_reading(self):
        """norb_phys = 2: interleaved index ``2 orb + spin``; an inter-orbital
        same-spin entry conserves, an inter-orbital spin-flip one does not."""
        same = {((0, 0, 0), (0, 2)): 0.4, ((0, 0, 0), (2, 0)): 0.4,
                ((0, 0, 0), (1, 3)): 0.4, ((0, 0, 0), (3, 1)): 0.4}
        self.assertEqual(band_conserves_sz(same, None, 2), BandSzCheck(True, None))
        flip = {**same, ((1, 0, 0), (2, 1)): 0.5j}
        self.assertEqual(band_conserves_sz(flip, None, 2),
                         BandSzCheck(False, ((1, 0, 0), (2, 1), 0.5j)))


class TestGeometryView(unittest.TestCase):
    """Construction check: the geometry feeds the bond kernel and the bond
    self-energy as a view; the public helpers equal their adapters."""

    SHAPE = CONVENTION_SHAPE
    ND_MODES = 2

    @classmethod
    def setUpClass(cls):
        cls.geom = so_bond.BondSlotGeometry.from_arrays(CONVENTION_DR, CONVENTION_REV, cls.SHAPE)
        cls.bond_set = bond_channels.BondSetView(
            bond_channels.BondTopology(CONVENTION_DR, CONVENTION_REV, {}))

    def test_n_channels(self):
        self.assertEqual(self.geom.n_channels, self.geom.B)
        self.assertEqual(self.geom.n_channels, len(CONVENTION_DR))
        self.assertIsInstance(self.geom.n_channels, int)
        with self.assertRaises(AttributeError):
            self.geom.n_channels = 3

    def test_n_channels_feeds_kernel(self):
        nd, B, nvol = self.ND_MODES, self.geom.B, self.geom.nvol
        ND = B * nd * nd
        rng = np.random.default_rng(31)
        blocks = rng.normal(size=(1, nvol, ND, ND)) + 1j * rng.normal(size=(1, nvol, ND, ND))
        g2_shape = (nd,) * 4 + self.SHAPE + (1,)
        G2 = rng.normal(size=g2_shape) + 1j * rng.normal(size=g2_shape)
        k_geom = pair_kernel(blocks, G2, self.geom, nd, self.SHAPE, 2.0)
        k_view = pair_kernel(blocks, G2, self.bond_set, nd, self.SHAPE, 2.0)
        phi = rng.normal(size=int(np.prod(k_geom.gap_shape))) + 0j
        assert_recomputed(self, k_geom.matvec(phi), k_view.matvec(phi))

    def test_n_channels_feeds_self_energy(self):
        nd, B, nvol, nmat = self.ND_MODES, self.geom.B, self.geom.nvol, 4
        P2 = nd * nd
        rng = np.random.default_rng(32)
        W = rng.normal(size=(nmat, nvol, B * P2, B * P2)) * (1 + 0.5j)

        class Store:
            def get_pair(self, name, a, b):
                assert name == "W"
                return W[:, :, a * P2:(a + 1) * P2, b * P2:(b + 1) * P2]

        g = random_green(33, nmat, nvol, nd)[None]
        s_geom = flex_bond.calc_self_energy_bond(Store(), g, 2.0, self.geom, self.SHAPE, nd, 1)
        s_view = flex_bond.calc_self_energy_bond(Store(), g, 2.0, self.bond_set, self.SHAPE,
                                                 nd, 1)
        assert_recomputed(self, s_geom, s_view)

    def test_public_equals_adapter(self):
        """Construction check: the public helpers equal their geometry
        adapters on the same inputs."""
        shape, nmat, beta = (5, 1, 1), 4, 2.0
        tables = {**_INTRA, "CoulombInter": closed_row(1, 0, 0, BASE["CoulombInter"])}
        c = compile_interaction(tables, 1)
        topo = so_bond_topology(c, CELL)
        geom = slot_geometry(topo, shape)
        nvol, nd = geom.nvol, 2
        ND = geom.B * nd * nd
        with quiet():
            plan = plan_so_bond(nmat=nmat, nvol=nvol, B=geom.B, nd=nd, consumers=("chi_static",),
                                caps=SOBondCaps(None, None, None))
        gb = so_bond_vertex(c, topo, shape, plan=plan)
        assert_recomputed(self, gb, so_bond._vertex_on_geom(c, geom, plan=plan))
        g = 0.1 * random_green(34, nmat, nvol, nd)[None]
        x_pub = so_bond_bubble(g, None, beta, topo, shape, plan=plan, freq="static",
                               sink=np.empty((nvol, ND, ND), np.complex128))
        x_ad = so_bond._bubble_on_geom(g, None, beta, geom, plan=plan, freq="static",
                                       sink=np.empty((nvol, ND, ND), np.complex128))
        assert_recomputed(self, x_pub, x_ad)
        diag = {}
        chi_pub = so_bond_dress(x_pub[None], gb, plan=plan, l0=nmat // 2, spatial_shape=shape)
        chi_ad = so_bond._dress_impl(x_pub[None], gb, plan=plan, l0=nmat // 2,
                                     spatial_shape=shape,
                                     cond_tol=bond_channels._BOND_COND_FLOOR, diagnostics=diag)
        assert_recomputed(self, chi_pub, chi_ad)
        _, cond_min = bond_channels.dress_batch(x_pub[None], gb, "spinful", l0=nmat // 2,
                                                nmat=nmat, spatial_shape=shape)
        self.assertEqual(diag, {"cond_min": cond_min})


if __name__ == "__main__":
    unittest.main()
