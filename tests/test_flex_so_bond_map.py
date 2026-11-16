"""One map of the spin-orbital bond FLEX (``flex_so_bond.so_bond_flex_map``,
gate C1b-M parts (a) and (b)).

The map is the bubble of the (tail-subtracted) Green function into the store
slot ``chibar``, per frequency batch the dressing, ``W = P_B[Gamma (chi -
chibar/2) Gamma]`` written over the batch (or into a second slot), and the
bond transport of ``W`` with the full Green function. Its parts are the
C1a-validated production helpers; these tests check the wiring.

* (a) truncations: with ``so_bond.so_bond_dress`` patched to return its
  input (``chi := chibar``) the map is the second-order chain
  ``c1a_sigma_prod(..., reading="p_b")``; patched to ``X - X Gamma X`` it is
  ``c1a(ring 2) - c1a(ring_order=3)`` (the physical third order is the
  NEGATIVE of the ``ring_order = 3`` chain). FIXTURE_I and FIXTURE_II with
  the "band" and "rand" Green functions (dense), the S_z-packed layout on
  ``PB_TABLES_SZ`` / FIXTURE_SZ_RING. The rejected C1a readings, patched in
  for ``so_bond.p_b``, miss.
* (b) full map, "band" Green function only (the random one is refused by
  the conditioning guard on FIXTURE_II): against a test-side composition (a separate
  store, ``np.linalg.solve(1 + chibar Gamma, chibar)``, ``p_b``, the
  transport), also with a tail (the bubble on ``G - tail``, the transport on
  ``G``); ``W`` in place equals ``W`` in a second slot; the transport's
  arguments; the guard minimum and the pole indicator; the channel-0 blocks;
  the caller's output buffers.

Every reference is built before a patch is entered.
"""
import functools
import math
import tracemalloc
import unittest
from unittest import mock

import numpy as np

from hwave.solver import eliashberg_bond as eb
from hwave.solver import flex_bond, flex_so_bond, so_bond
from hwave.solver.so_bond import (PackedArrays, plan_so_bond, select_layout, slot_geometry,
                                  so_bond_topology, so_bond_vertex, sz_sector_map)
from hwave.solver.so_records import compile_interaction
from tests.so_bond_fixtures import (C1A_NMAT_FAST, NOCAP, PB_TABLES_SZ, READINGS, c1a_green,
                                    c1a_sigma_prod, copy_packed)
from tests.so_bond_oracles import (CELL, FIXTURE_I, FIXTURE_II, FIXTURE_SZ_RING, OracleGeometry,
                                   geom_of, rel)

NMAT = C1A_NMAT_FAST
NB = 3
#: gate tolerances (fixed) and the floor of a deliberate break
TOL_SIGMA = 1e-12
TOL_CHI = 1e-13
MISS = 1e-3
#: two runs of the same production arithmetic (another store, another call):
#: bitwise on numpy 2 (macOS); the numpy 1.26.4 OpenBLAS build (DYNAMIC_ARCH)
#: of the local test environment is not run-to-run reproducible, also with
#: one thread (measured: 1.4e-16 on the bubble alone, 4.4e-16 on sigma)
REPRO_TOL = 1e-14
#: the guard of the unpatched maps: the dressing's default
COND_TOL = so_bond.bond_channels._BOND_COND_FLOOR

#: cases: (label, fixture, tables key, layout)
DENSE_CASES = (("I", FIXTURE_I, "fx", "dense"), ("II", FIXTURE_II, "fx", "dense"))
PACKED_CASES = (("SZ ring", FIXTURE_SZ_RING, ("pb", False), "sz_packed"),
                ("SZ ring nnn", FIXTURE_SZ_RING, ("pb", True), "sz_packed"))
CASES = DENSE_CASES + PACKED_CASES
#: unpatched ("band") rows: one dense case per fixture and the packed ring
FULL_CASES = (DENSE_CASES[0], DENSE_CASES[1], PACKED_CASES[0])


def _tables(fx, key):
    return fx.tables_dict if key == "fx" else PB_TABLES_SZ[key[1]]


def _setup(fx, key, layout, w_in_place=None):
    """``(plan, gamma_b, geom, sector_map)`` of a case (``nb = 3``), built
    once per argument set."""
    return _setup_cached(fx, key, layout, w_in_place)


@functools.lru_cache(maxsize=None)
def _setup_cached(fx, key, layout, w_in_place):
    c = compile_interaction(_tables(fx, key), fx.norb_phys)
    topo = so_bond_topology(c, CELL)
    geom = slot_geometry(topo, fx.spatial_shape)
    reason, sm = None, None
    if layout == "sz_packed":
        selected, reason = select_layout(c, fx.transfer_dict, override="auto")
        assert selected == "sz_packed", reason
        sm = sz_sector_map(topo, fx.nd)
    plan = plan_so_bond(nmat=NMAT, nvol=fx.nvol, B=len(topo.delta_r), nd=fx.nd,
                        consumers=("flex_scf",), caps=NOCAP, layout=layout,
                        layout_reason=reason, nb=NB, w_in_place=w_in_place)
    gamma = so_bond_vertex(c, topo, fx.spatial_shape, plan=plan)
    return plan, gamma, geom, sm


def _green(fx, kind, layout):
    """``c1a_green`` ``(nmat, nvol, nd, nd)``; for the packed layout with the
    spin-off-diagonal (S_z off-sector) entries set to exactly zero."""
    g = c1a_green(fx, kind, NMAT)
    if layout == "sz_packed":
        n = fx.norb_phys
        g[..., :n, n:] = 0
        g[..., n:, :n] = 0
    return g


def _map(case, green_kw, green_scf=None, tail=None, *, cond_tol=None, w_in_place=None,
         workers=1, **kw):
    """The map of ``case`` in a fresh store (``green_scf`` defaults to
    ``green_kw``; both ``(nmat, nvol, nd, nd)``; ``tail`` rank 5 or None)."""
    _, fx, key, layout = case
    plan, gamma, geom, sm = _setup(fx, key, layout, w_in_place)
    gk = np.ascontiguousarray(green_kw)[None].copy()
    gs = gk if green_scf is None else np.ascontiguousarray(green_scf)[None].copy()
    with so_bond.open_store(plan, sm) as st:
        return flex_so_bond.so_bond_flex_map(st, gamma, gk, gs, tail, fx.beta, geom, plan=plan,
                                             cond_tol=cond_tol, workers=workers, **kw)


def _dress_identity(packed):
    def dress(X, gamma_b, **kw):
        return copy_packed(X) if packed else X.copy()
    return dress


def _dress_third(packed):
    """``X - X Gamma X`` (staged two-operand matmul), per sector when packed."""
    def one(x, g):
        t = np.matmul(x, g[None])
        return x - np.matmul(t, x)

    def dress(X, gamma_b, **kw):
        if packed:
            return PackedArrays(X.sector_map, tuple(one(x, g) for x, g in
                                                    zip(X.blocks, gamma_b.blocks)))
        return one(X, gamma_b)
    return dress


def _composition(case, g_bubble, tail, g_transport, *, swap_c0=False):
    """The test-side map, dense: ``chibar`` of ``g_bubble`` (and ``tail``) in
    a separate store, ``chi = solve(1 + chibar Gamma, chibar)`` per ``(l,
    q)``, ``W = p_b(Gamma (chi - chibar/2) Gamma)``, the transport of ``W``
    with ``g_transport``. Returns ``(sigma (nmat, nvol, nd, nd), chi0_c0,
    chi_c0, chibar)``."""
    _, fx, key, _ = case
    c = compile_interaction(_tables(fx, key), fx.norb_phys)
    topo = so_bond_topology(c, CELL)
    geom = slot_geometry(topo, fx.spatial_shape)
    plan = plan_so_bond(nmat=NMAT, nvol=fx.nvol, B=len(topo.delta_r), nd=fx.nd,
                        consumers=("chi_dynamic",), caps=NOCAP, layout="dense")
    gamma = so_bond_vertex(c, topo, fx.spatial_shape, plan=plan)
    ND, P = plan.ND, fx.nd ** 2
    with flex_bond.BondBlockStore(NMAT, fx.nvol, ND, P, ("chibar",)) as st:
        so_bond.so_bond_bubble(np.ascontiguousarray(g_bubble)[None].copy(), tail, fx.beta, topo,
                               fx.spatial_shape, plan=plan, freq="all", sink=st)
        chibar = st.detach("chibar")
    g = gamma[None]
    chi = np.linalg.solve(np.eye(ND) + np.matmul(chibar, g), chibar)
    M = np.matmul(np.matmul(g, chi - 0.5 * chibar), g)
    W = so_bond.p_b(M, fx.nd)
    sigma = flex_bond.calc_self_energy_bond(eb.ArrayBlockSource({"W": W}, P),
                                            np.ascontiguousarray(g_transport)[None].copy(),
                                            fx.beta, geom, fx.spatial_shape, fx.nd, 1)[0]
    return sigma, chibar[..., :P, :P].copy(), chi[..., :P, :P].copy(), chibar


def _same_run(tc, a, b, what):
    """``a`` and ``b`` are the same production arithmetic (:data:`REPRO_TOL`)."""
    a, b = np.asarray(a), np.asarray(b)
    tc.assertEqual(a.shape, b.shape, what)
    if not np.any(b):
        tc.assertFalse(np.any(a), what)
        return
    tc.assertLessEqual(rel(a, b), REPRO_TOL, what)


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


def _finite(tc, *values):
    for v in values:
        tc.assertTrue(np.isfinite(v), "non-finite metric")


class TestTruncations(unittest.TestCase):
    """C1b-M (a): the map with the dressing truncated at second and third
    order equals the C1a production chain."""

    def test_second_order(self):
        rows = []
        for case in CASES:
            label, fx, key, layout = case
            packed = layout == "sz_packed"
            for kind in ("band", "rand"):
                with self.subTest(case=label, kind=kind):
                    G = _green(fx, kind, layout)
                    spy = mock.Mock(side_effect=_dress_identity(packed))
                    ref = c1a_sigma_prod(fx, _tables(fx, key), G, reading="p_b")
                    self.assertEqual(spy.call_count, 0)
                    with mock.patch.object(so_bond, "so_bond_dress", spy):
                        res = _map(case, G)
                    self.assertEqual(spy.call_count, math.ceil(NMAT / NB))
                    dev = rel(res.sigma_fluct[0], ref)
                    rows.append((label, kind, dev))
                    _finite(self, dev)
                    self.assertLessEqual(dev, TOL_SIGMA)
        print("\nC1b-M (a) second order:\n" + "\n".join(
            "  %-12s %-4s %.3e" % r for r in rows))

    def test_third_order(self):
        rows = []
        for case in CASES:
            label, fx, key, layout = case
            packed = layout == "sz_packed"
            for kind in ("band", "rand"):
                with self.subTest(case=label, kind=kind):
                    G = _green(fx, kind, layout)
                    spy = mock.Mock(side_effect=_dress_third(packed))
                    tables = _tables(fx, key)
                    s2 = c1a_sigma_prod(fx, tables, G)
                    s3 = c1a_sigma_prod(fx, tables, G, ring_order=3)
                    self.assertEqual(spy.call_count, 0)
                    ref = s2 - s3
                    size3 = float(np.abs(s3).max() / np.abs(s2).max())
                    with mock.patch.object(so_bond, "so_bond_dress", spy):
                        res = _map(case, G)
                    self.assertEqual(spy.call_count, math.ceil(NMAT / NB))
                    dev = rel(res.sigma_fluct[0], ref)
                    rows.append((label, kind, dev, size3))
                    _finite(self, dev, size3)
                    self.assertGreaterEqual(size3, MISS)        # non-vacuity
                    self.assertLessEqual(dev, TOL_SIGMA)
        print("\nC1b-M (a) third order (dev, max|c1a(3)|/max|c1a(2)|):\n" + "\n".join(
            "  %-12s %-4s %.3e  %.3e" % r for r in rows))

    def test_wrong_readings_miss(self):
        case = DENSE_CASES[0]
        _, fx, key, layout = case
        geom = _setup(fx, key, layout)[2]
        G = _green(fx, "rand", layout)
        ref = c1a_sigma_prod(fx, _tables(fx, key), G, reading="p_b")
        rows = []
        for name, reading in sorted(READINGS.items()):
            if name == "p_b":
                continue
            with self.subTest(reading=name):
                def wrong(M, nd, *, out=None, reading=reading):
                    W = np.array(reading(M, geom, nd), dtype=np.complex128)
                    if out is None:
                        return W
                    out[...] = W
                    return out
                spy = mock.Mock(side_effect=wrong)
                with mock.patch.object(so_bond, "so_bond_dress", _dress_identity(False)), \
                        mock.patch.object(so_bond, "p_b", spy):
                    res = _map(case, G)
                self.assertEqual(spy.call_count, math.ceil(NMAT / NB))     # applied
                dev = rel(res.sigma_fluct[0], ref)
                rows.append((name, dev))
                _finite(self, dev)
                self.assertGreaterEqual(dev, MISS)
        print("\nC1b-M (a) rejected readings (FIXTURE_I rand):\n" + "\n".join(
            "  %-18s %.3e" % r for r in rows))


class TestFullMap(unittest.TestCase):
    """C1b-M (b): the unpatched map ("band" Green function)."""

    def test_matches_composition(self):
        rows = []
        for case in FULL_CASES:
            label, fx, key, layout = case
            with self.subTest(case=label):
                G = _green(fx, "band", layout)
                sigma, c0, c, _ = _composition(case, G, None, G)
                res = _map(case, G, cond_tol=COND_TOL)
                d = (rel(res.sigma_fluct[0], sigma), rel(res.chi0_c0, c0), rel(res.chi_c0, c))
                rows.append((label,) + d + (res.cond_min,))
                _finite(self, *d)
                self.assertLessEqual(d[0], TOL_SIGMA)
                self.assertLessEqual(d[1], TOL_CHI)
                self.assertLessEqual(d[2], TOL_CHI)
                # the dressing is not the identity here
                self.assertGreaterEqual(rel(c, c0), MISS)
        print("\nC1b-M (b) full map (sigma, chi0_c0, chi_c0, cond_min):\n" + "\n".join(
            "  %-12s %.3e  %.3e  %.3e  %.3e" % r for r in rows))

    def test_tail(self):
        case = DENSE_CASES[0]
        _, fx, key, layout = case
        coeff, beta, nd = 1.0, fx.beta, fx.nd
        G_band = _green(fx, "band", layout)
        w = (2 * np.arange(NMAT) + 1 - NMAT) * np.pi / beta
        eye = np.eye(nd)
        tail_w = np.broadcast_to((coeff / (1j * w))[:, None, None, None] * eye,
                                 G_band.shape).copy()
        tail = np.broadcast_to(coeff * beta / 2 * eye, (1,) + G_band.shape).astype(np.complex128)
        G = G_band + tail_w
        Gs = G - tail_w
        sigma, c0, c, _ = _composition(case, Gs, tail, G)
        res = _map(case, G, Gs, tail, cond_tol=COND_TOL)
        d = (rel(res.sigma_fluct[0], sigma), rel(res.chi0_c0, c0), rel(res.chi_c0, c))
        swapped = _map(case, Gs, G, tail, cond_tol=COND_TOL)
        miss = rel(swapped.sigma_fluct[0], sigma)
        # the tail matters: the bubble of G without the tail differs
        no_tail = _composition(case, Gs, None, G)[0]
        tail_effect = rel(no_tail, sigma)
        print("\nC1b-M (b) tail (sigma, chi0_c0, chi_c0; swapped miss; tail effect): "
              "%.3e %.3e %.3e; %.3e; %.3e" % (d + (miss, tail_effect)))
        _finite(self, *d, miss, tail_effect)
        self.assertLessEqual(d[0], TOL_SIGMA)
        self.assertLessEqual(d[1], TOL_CHI)
        self.assertLessEqual(d[2], TOL_CHI)
        self.assertGreaterEqual(miss, MISS)
        self.assertGreaterEqual(tail_effect, MISS)

    def test_in_place_equals_two_slots(self):
        for case in (DENSE_CASES[0], PACKED_CASES[0]):
            label, fx, key, layout = case
            with self.subTest(case=label):
                self.assertEqual(_setup(fx, key, layout)[0].aliases, (("W", "chibar"),))
                two = _setup(fx, key, layout, False)[0]
                self.assertEqual((two.slots, two.aliases), (("chibar", "W"), ()))
                G = _green(fx, "band", layout)
                a = _map(case, G, cond_tol=COND_TOL)
                b = _map(case, G, cond_tol=COND_TOL, w_in_place=False)
                for f in ("sigma_fluct", "chi0_c0", "chi_c0"):
                    _same_run(self, getattr(a, f), getattr(b, f), f)
                _same_run(self, a.cond_min, b.cond_min, "cond_min")
                self.assertEqual((a.sigma_min_sector, a.sigma_max_sector),
                                 (b.sigma_min_sector, b.sigma_max_sector))
                _same_pole(self, a.pole, b.pole)

    def test_transport_arguments(self):
        case = DENSE_CASES[0]
        _, fx, key, layout = case
        _, _, geom, _ = _setup(fx, key, layout)
        G = _green(fx, "band", layout)
        ref = _map(case, G, cond_tol=COND_TOL).sigma_fluct
        real = flex_bond.calc_self_energy_bond
        spy = mock.Mock(side_effect=real)
        with mock.patch.object(flex_bond, "calc_self_energy_bond", spy):
            res = _map(case, G, cond_tol=COND_TOL, workers=2)
        self.assertEqual(spy.call_count, 1)
        args = spy.call_args[0]
        self.assertEqual(len(args), 7)
        self.assertIs(args[3], geom)
        self.assertEqual(args[4], fx.spatial_shape)
        self.assertEqual(args[5], fx.nd)
        self.assertEqual(args[6], 2)
        self.assertEqual(args[2], fx.beta)
        self.assertTrue(np.array_equal(args[1], G[None]))
        self.assertLessEqual(rel(res.sigma_fluct, ref), TOL_SIGMA)
        # the label geometry of FIXTURE_I with channels 1 and 2 swapped (the
        # label geometry itself coincides with the slot one) misses
        label = geom_of(fx)
        self.assertTrue(np.array_equal(label.delta_r, geom.delta_r))
        perm = np.array([0, 2, 1, 3, 4])
        inv = np.argsort(perm)
        swapped = OracleGeometry(np.asarray(label.delta_r)[perm],
                                 inv[np.asarray(label.reverse)[perm]], fx.spatial_shape)

        def on_label(store, green, beta, view, shape, norb, workers):
            return real(store, green, beta, swapped, shape, norb, workers)
        spy2 = mock.Mock(side_effect=on_label)
        with mock.patch.object(flex_bond, "calc_self_energy_bond", spy2):
            bad = _map(case, G, cond_tol=COND_TOL)
        self.assertEqual(spy2.call_count, 1)
        miss = rel(bad.sigma_fluct, ref)
        print("\nC1b-M (b) transport: workers 2 vs 1 %.3e; label geometry miss %.3e"
              % (rel(res.sigma_fluct, ref), miss))
        _finite(self, miss)
        self.assertGreaterEqual(miss, MISS)

    def test_cond_and_pole(self):
        for case in (DENSE_CASES[0], DENSE_CASES[1], PACKED_CASES[0]):
            label, fx, key, layout = case
            with self.subTest(case=label):
                plan, gamma, geom, sm = _setup(fx, key, layout)
                G = _green(fx, "band", layout)
                chibar = _composition(case, G, None, G)[3]
                static = chibar[NMAT // 2]
                if sm is not None:
                    static = PackedArrays.from_dense(static, sm, source="reference")
                want = so_bond.static_pole_indicator(static, gamma, fx.spatial_shape)
                real = so_bond.so_bond_dress
                seen = []

                def dress(*a, **kw):
                    out = real(*a, **kw)
                    seen.append(dict(kw["diagnostics"]))
                    return out
                with mock.patch.object(so_bond, "so_bond_dress", side_effect=dress):
                    res = _map(case, G, cond_tol=COND_TOL)
                self.assertEqual(len(seen), math.ceil(NMAT / NB))
                worst = min(seen, key=lambda d: d["cond_min"])
                self.assertEqual(res.cond_min, worst["cond_min"])
                self.assertEqual(res.sigma_min_sector, worst.get("sigma_min_sector"))
                self.assertEqual(res.sigma_max_sector, worst.get("sigma_max_sector"))
                if sm is None:
                    self.assertIsNone(res.sigma_min_sector)
                else:
                    self.assertIn(res.sigma_min_sector, ("0", "+1", "-1"))
                _same_pole(self, res.pole, want)
                self.assertTrue(np.isfinite(res.pole[0][0]))
                # without the guard
                off = _map(case, G, cond_tol=None)
                self.assertTrue(math.isnan(off.cond_min))
                self.assertIsNone(off.sigma_min_sector)
                self.assertIsNone(off.sigma_max_sector)

    def test_mutations(self):
        """``channel0_block`` of ``X`` and of ``Y`` swapped: ``chi_c0`` (and
        ``chi0_c0``) miss."""
        case = DENSE_CASES[0]
        label, fx, key, layout = case
        _, gamma, _, _ = _setup(fx, key, layout)
        G = _green(fx, "band", layout)
        _, c0, c, _ = _composition(case, G, None, G)
        real = so_bond.channel0_block
        ND = gamma.shape[-1]
        calls = []

        def swapped(A, *, nd=None):
            calls.append(A)
            if len(calls) % 2:                         # the X call gets Y's block
                Y = np.linalg.solve(np.eye(ND) + np.matmul(A, gamma[None]), A)
                return real(Y, nd=nd)
            return real(calls[-2], nd=nd)               # the Y call gets X's block
        spy = mock.Mock(side_effect=swapped)
        with mock.patch.object(so_bond, "channel0_block", spy):
            res = _map(case, G, cond_tol=COND_TOL)
        self.assertEqual(spy.call_count, 2 * math.ceil(NMAT / NB))
        m = (rel(res.chi_c0, c), rel(res.chi0_c0, c0))
        print("\nC1b-M (b) channel-0 swap miss (chi_c0, chi0_c0): %.3e %.3e" % m)
        _finite(self, *m)
        self.assertGreaterEqual(m[0], MISS)
        self.assertGreaterEqual(m[1], MISS)

    def test_out_buffers(self):
        case = DENSE_CASES[1]                          # FIXTURE_II: nd = 4
        _, fx, key, layout = case
        plan, _, _, _ = _setup(fx, key, layout)
        G = _green(fx, "band", layout)
        P = fx.nd ** 2
        shape = (NMAT, fx.nvol, P, P)
        C = 16 * int(np.prod(shape))
        chi0_out = np.full(shape, np.nan + 0j)
        chi_out = np.full(shape, np.nan + 0j)
        ref = _map(case, G, cond_tol=COND_TOL)            # warm the caches
        res = _map(case, G, cond_tol=COND_TOL, chi0_out=chi0_out, chi_out=chi_out)
        self.assertIs(res.chi0_c0, chi0_out)
        self.assertIs(res.chi_c0, chi_out)
        for f in ("sigma_fluct", "chi0_c0", "chi_c0"):
            _same_run(self, getattr(res, f), getattr(ref, f), f)
        S = ref.sigma_fluct.nbytes

        def retained(**kw):
            tracemalloc.start()
            try:
                before = tracemalloc.get_traced_memory()[0]
                out = _map(case, G, cond_tol=COND_TOL, **kw)
                after = tracemalloc.get_traced_memory()[0]
            finally:
                tracemalloc.stop()
            return after - before, out
        with_out, _ = retained(chi0_out=chi0_out, chi_out=chi_out)
        without, _ = retained()
        print("\nC1b-M (b) retained bytes of a map: with buffers %d, without %d "
              "(sigma %d, one channel-0 block C = %d)" % (with_out, without, S, C))
        # with the buffers only sigma (and small objects) survives the map:
        # no channel-0 block of its own
        self.assertLess(with_out, S + C // 2)
        # without them both blocks are allocated (the measurement sees them)
        self.assertGreaterEqual(without, S + 2 * C)

    def test_buffer_refusals(self):
        case = DENSE_CASES[0]
        _, fx, key, layout = case
        P = fx.nd ** 2
        shape = (NMAT, fx.nvol, P, P)
        G = _green(fx, "band", layout)
        with self.assertRaisesRegex(ValueError, r"^so_bond_flex_map: chi0_out must be a "
                                    r"complex128 array of shape"):
            _map(case, G, chi0_out=np.zeros(shape))
        shared = np.zeros(shape, dtype=np.complex128)
        with self.assertRaisesRegex(ValueError, "^so_bond_flex_map: chi0_out and chi_out "
                                    "share memory$"):
            _map(case, G, chi0_out=shared, chi_out=shared)


if __name__ == "__main__":
    unittest.main()
