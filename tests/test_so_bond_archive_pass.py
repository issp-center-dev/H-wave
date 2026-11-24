"""The susceptibility pass of the spin-orbital bond archive
(``flex_so_bond.so_bond_archive_pass``), at the level of the block store.

After the self-consistency loop the pass computes ``chi_B`` of the final
Green function: the bubble of ``green_scf`` (with its tail) into the store
slot ``chibar`` and, per frequency batch, the static pole indicator of the
``nu = 0`` slice (before that batch is dressed) and the dressing ``chi_B =
[1 + chibar_B Gamma_B]^{-1} chibar_B`` written over the batch of
``chibar``. These tests check the slot against a test-side composition
(``np.linalg.solve``), the bubble against the discrete label sum of
``tests/so_bond_oracles.py``, the pole indicator against the one of the FLEX
map on the same Green function, the trace scopes and the refusals.

Every reference is built before a patch is entered.
"""
import functools
import math
import unittest
from unittest import mock

import numpy as np

from hwave.solver import flex_bond, flex_so_bond, so_bond
from hwave.solver.bond_channels import BondConditioningError
from hwave.solver.so_bond import (PackedArrays, plan_so_bond, select_layout, slot_geometry,
                                  so_bond_topology, so_bond_vertex, sz_sector_map)
from hwave.solver.so_records import compile_interaction
from tests.so_bond_fixtures import C1A_NMAT_FAST, NOCAP, PB_TABLES_SZ, c1a_green
from tests.so_bond_oracles import (CELL, FIXTURE_I, FIXTURE_SZ_RING, bubble_label_to_code,
                                   free_green, lindhard_discrete, rel)

NMAT = C1A_NMAT_FAST
NB = 3
#: gate tolerance of the slot and of the bubble (fixed)
TOL = 1e-12
#: the floor of a deliberate break and of a non-vacuity check
MISS = 1e-3
#: two runs of the same production arithmetic (numpy 1.26.4 OpenBLAS is not
#: run-to-run reproducible, see tests/test_flex_so_bond_map.py)
REPRO_TOL = 1e-14
#: the guard of the unpatched passes: the dressing's default
COND_TOL = so_bond.bond_channels._BOND_COND_FLOOR
#: the label view of the slot convention (as tests/test_so_bond_chi_ring.py)
FLAGS = so_bond.SLOT_CONVENTION.label_view()

#: cases: (label, fixture, tables key, layout)
DENSE = ("I", FIXTURE_I, "fx", "dense")
PACKED = ("I S_z ring", FIXTURE_SZ_RING, ("pb", False), "sz_packed")
CASES = (DENSE, PACKED)


def _tables(fx, key):
    return fx.tables_dict if key == "fx" else PB_TABLES_SZ[key[1]]


@functools.lru_cache(maxsize=None)
def _setup(fx, key, layout, keep_chi, nmat=NMAT):
    """``(plan, gamma_b, geom, sector_map)`` of a case (consumer
    ``flex_scf``, ``nb = 3``)."""
    c = compile_interaction(_tables(fx, key), fx.norb_phys)
    topo = so_bond_topology(c, CELL)
    geom = slot_geometry(topo, fx.spatial_shape)
    reason, sm = None, None
    if layout == "sz_packed":
        selected, reason = select_layout(c, fx.transfer_dict, override="auto")
        assert selected == "sz_packed", reason
        sm = sz_sector_map(topo, fx.nd)
    plan = plan_so_bond(nmat=nmat, nvol=fx.nvol, B=len(topo.delta_r), nd=fx.nd,
                        consumers=("flex_scf",), caps=NOCAP, layout=layout,
                        layout_reason=reason, nb=NB, keep_chi=keep_chi)
    gamma = so_bond_vertex(c, topo, fx.spatial_shape, plan=plan)
    return plan, gamma, geom, sm


def _green(fx, layout, nmat=NMAT):
    """The "band" Green function ``(1, nmat, nvol, nd, nd)``; for the packed
    layout with the S_z off-sector entries exactly zero."""
    g = c1a_green(fx, "band", nmat)
    if layout == "sz_packed":
        n = fx.norb_phys
        g[..., :n, n:] = 0
        g[..., n:, :n] = 0
    return np.ascontiguousarray(g)[None].copy()


def _dense(a):
    return a.to_dense() if isinstance(a, PackedArrays) else np.asarray(a)


def _pass(case, G, *, tail=None, cond_tol=COND_TOL, keep_chi=True, plan=None, nmat=NMAT,
          trace=None):
    """The pass of ``case`` in a fresh store: ``(result, the detached
    chibar slot as a dense array)``."""
    _, fx, key, layout = case
    p, gamma, geom, sm = _setup(fx, key, layout, keep_chi, nmat)
    plan = p if plan is None else plan
    with so_bond.open_store(plan, sm) as st:
        res = flex_so_bond.so_bond_archive_pass(st, gamma, G, tail, fx.beta, geom, plan=plan,
                                                cond_tol=cond_tol, trace=trace)
        slot = _dense(st.detach("chibar"))
    return res, slot


def _composition(case, G, tail=None):
    """The test-side pass, dense: ``chibar`` of ``G`` (and ``tail``) in a
    separate store, ``chi = solve(1 + chibar Gamma, chibar)`` per ``(l,
    q)``. Returns ``(chi, chibar)``."""
    _, fx, key, _ = case
    c = compile_interaction(_tables(fx, key), fx.norb_phys)
    topo = so_bond_topology(c, CELL)
    plan = plan_so_bond(nmat=NMAT, nvol=fx.nvol, B=len(topo.delta_r), nd=fx.nd,
                        consumers=("chi_dynamic",), caps=NOCAP, layout="dense")
    gamma = so_bond_vertex(c, topo, fx.spatial_shape, plan=plan)
    ND = plan.ND
    with flex_bond.BondBlockStore(NMAT, fx.nvol, ND, fx.nd ** 2, ("chibar",)) as st:
        so_bond.so_bond_bubble(G.copy(), tail, fx.beta, topo, fx.spatial_shape, plan=plan,
                               freq="all", sink=st)
        chibar = st.detach("chibar")
    chi = np.linalg.solve(np.eye(ND) + np.matmul(chibar, gamma[None]), chibar)
    return chi, chibar


def _finite(tc, *values):
    for v in values:
        tc.assertTrue(np.isfinite(v), "non-finite metric")


def _same_entry(tc, a, b, what):
    """Pole-indicator entries ``(lambda, q, max_re)``."""
    tc.assertEqual(a[1], b[1], what)
    for x, y in ((a[0], b[0]), (a[2], b[2])):
        tc.assertEqual(np.isnan(x), np.isnan(y), what)
        if not np.isnan(y):
            tc.assertLessEqual(abs(x - y), REPRO_TOL * max(1.0, abs(y)), what)


def _same_pole(tc, a, b):
    """``(entry, per_sector, worst)`` of :func:`so_bond.static_pole_indicator`."""
    _same_entry(tc, a[0], b[0], "pole entry")
    tc.assertEqual(a[2], b[2], "worst sector")
    if b[1] is None:
        tc.assertIsNone(a[1])
        return
    tc.assertEqual(sorted(a[1]), sorted(b[1]))
    for k in b[1]:
        _same_entry(tc, a[1][k], b[1][k], "pole sector " + k)


class TestArchivePass(unittest.TestCase):

    def test_matches_composition(self):
        """The slot after the pass is ``chi_B`` of the Green function; with
        the dressing patched to the identity it is the bubble (miss)."""
        rows = []
        for case in CASES:
            label, fx, key, layout = case
            packed = layout == "sz_packed"
            with self.subTest(case=label):
                G = _green(fx, layout)
                chi, chibar = _composition(case, G)
                res, slot = _pass(case, G)
                dev = rel(slot, chi)
                nonvac = rel(chi, chibar)
                _finite(self, dev, nonvac, res.cond_min)
                self.assertLessEqual(dev, TOL)
                self.assertGreaterEqual(nonvac, MISS)       # the dressing is not the identity
                # break: the dressing returns its input
                spy = mock.Mock(side_effect=lambda X, gamma_b, **kw:
                                (PackedArrays(X.sector_map, tuple(b.copy() for b in X.blocks))
                                 if packed else X.copy()))
                with mock.patch.object(so_bond, "so_bond_dress", spy):
                    _, broken = _pass(case, G)
                self.assertEqual(spy.call_count, math.ceil(NMAT / NB))
                miss = rel(broken, chi)
                _finite(self, miss)
                self.assertGreaterEqual(miss, MISS)
                rows.append((label, dev, nonvac, miss, res.cond_min))
        print("\narchive pass vs composition (dev, |chi - chibar|, identity-dressing miss, "
              "cond_min):\n" + "\n".join("  %-12s %.3e  %.3e  %.3e  %.3e" % r for r in rows))

    def test_tail(self):
        """With a tail (``coeff_tail = 1``) the bubble runs on ``G - tail``
        with ``green0_tail``, as in the FLEX map; the tail matters."""
        _, fx, _, layout = DENSE
        coeff, beta, nd = 1.0, fx.beta, fx.nd
        # the tail-subtracted Green function G - tail_w of G = G_band +
        # coeff / (i w) is G_band; the tail argument is the one of the map
        # tests (tests/test_flex_so_bond_map.py::TestFullMap::test_tail)
        Gs = _green(fx, layout)
        tail = np.broadcast_to(coeff * beta / 2 * np.eye(nd), Gs.shape).astype(np.complex128)
        chi, _ = _composition(DENSE, Gs, tail)
        no_tail, _ = _composition(DENSE, Gs)
        res, slot = _pass(DENSE, Gs, tail=tail)
        dev, effect = rel(slot, chi), rel(no_tail, chi)
        print("\narchive pass tail (dev, tail effect): %.3e %.3e" % (dev, effect))
        _finite(self, dev, effect, res.cond_min)
        self.assertLessEqual(dev, TOL)
        self.assertGreaterEqual(effect, MISS)

    def test_requires_keep_chi(self):
        _, fx, key, layout = DENSE
        G = _green(fx, layout)
        with self.assertRaises(ValueError) as cm:
            _pass(DENSE, G, keep_chi=False)
        self.assertEqual(str(cm.exception),
                         "so_bond_archive_pass: the plan must be built with keep_chi = True "
                         "(plan_so_bond(..., keep_chi=True))")
        other = plan_so_bond(nmat=NMAT, nvol=fx.nvol, B=_setup(fx, key, layout, True)[0].B,
                             nd=fx.nd, consumers=("chi_dynamic",), caps=NOCAP, layout="dense",
                             keep_chi=True)
        with self.assertRaises(ValueError) as cm:
            _pass(DENSE, G, plan=other)
        self.assertEqual(str(cm.exception),
                         "so_bond_archive_pass: the plan needs the consumer flex_scf (closed "
                         "consumers: vertex, chi_dynamic)")
        with self.assertRaises(ValueError) as cm:
            _pass(DENSE, G[:, :, :, :1])
        self.assertTrue(str(cm.exception).startswith(
            "so_bond_archive_pass: green_scf must have shape"), str(cm.exception))

    def test_pole_matches_map(self):
        """The pole indicator of the pass (``keep_chi = True`` plan) equals
        the map's (``keep_chi = False`` plan) on the same Green function."""
        for case in CASES:
            label, fx, key, layout = case
            with self.subTest(case=label):
                G = _green(fx, layout)
                plan, gamma, geom, sm = _setup(fx, key, layout, False)
                self.assertFalse(plan.keep_chi)
                with so_bond.open_store(plan, sm) as st:
                    ref = flex_so_bond.so_bond_flex_map(st, gamma, G, G, None, fx.beta, geom,
                                                        plan=plan, cond_tol=COND_TOL)
                res, _ = _pass(case, G)
                self.assertTrue(_setup(fx, key, layout, True)[0].keep_chi)
                _same_pole(self, res.pole, ref.pole)
                self.assertTrue(np.isfinite(res.pole[0][0]))
                if layout == "sz_packed":
                    self.assertIsNotNone(res.pole[1])
                _finite(self, res.cond_min)
                self.assertLessEqual(abs(res.cond_min - ref.cond_min),
                                     REPRO_TOL * abs(ref.cond_min))

    def test_bubble_independent(self):
        """The production bubble of the pass (the dressing's input) equals
        the bubble-rule image of the discrete label sum (FIXTURE_I, Nmat 4,
        dense, no tail)."""
        nmat = 4
        _, fx, key, layout = DENSE
        geom = _setup(fx, key, layout, True, nmat)[2]
        self.assertTrue(np.array_equal(fx.geom.delta_r, geom.delta_r))
        self.assertTrue(np.array_equal(fx.geom.reverse, geom.reverse))
        g = free_green(fx.band[0], fx.band[1], fx.beta, nmat)
        want = bubble_label_to_code(lindhard_discrete(g, fx.beta, fx.geom), fx.geom, fx.nd,
                                    FLAGS)
        real = so_bond.so_bond_dress
        seen = []

        def dress(X, gamma_b, **kw):
            seen.append((kw["l0"], X.copy()))
            return real(X, gamma_b, **kw)
        spy = mock.Mock(side_effect=dress)
        with mock.patch.object(so_bond, "so_bond_dress", spy):
            _pass(DENSE, np.ascontiguousarray(g)[None].copy(), nmat=nmat)
        self.assertEqual(spy.call_count, math.ceil(nmat / NB))
        self.assertEqual([l0 for l0, _ in seen], list(range(0, nmat, NB)))
        got = np.concatenate([x for _, x in seen])
        self.assertEqual(got.shape, want.shape)
        dev = rel(got, want)
        print("\narchive pass bubble vs discrete label sum: %.3e" % dev)
        _finite(self, dev)
        self.assertLessEqual(dev, TOL)

    def test_trace_scopes(self):
        _, fx, _, layout = DENSE
        names = []

        def trace(name):
            names.append(name)
            return mock.MagicMock()
        _pass(DENSE, _green(fx, layout), trace=trace)
        self.assertEqual(names, ["archive_bubble"]
                         + ["archive_dressing"] * math.ceil(NMAT / NB))

    def test_cond_refusal_propagates(self):
        _, fx, _, layout = DENSE
        with self.assertRaises(BondConditioningError):
            _pass(DENSE, _green(fx, layout), cond_tol=0.99)


if __name__ == "__main__":
    unittest.main()
