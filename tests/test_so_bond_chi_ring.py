"""The spinful bond bubble and its RPA dressing (``so_bond.so_bond_bubble``,
``so_bond.so_bond_dress``) against label-form oracles.

Code slots follow ``so_bond.SLOT_CONVENTION``. The bubble is
``chibar_B(q, i nu)[(m;x,y),(m';z,w)] = e^{-iq.dr_m}
chibar(-q, +i nu)[(+dr_m; y,x), (-dr_m'; z,w)]`` (the row operator
daggered), and the dressing ``chi_B = [1 + chibar_B Gamma_B]^{-1} chibar_B``
equals the same bubble-rule image of the label ring
``[1 + chibar(q) Gamma(-q)]^{-1} chibar(q)``, with no further momentum flip.
"""
import functools
import tracemalloc
import unittest
from types import SimpleNamespace
from typing import NamedTuple

import numpy as np

from hwave.solver import flex_bond
from hwave.solver.hartree_fock import NonFiniteError
from hwave.solver.bond_channels import BondSetView, BondTopology
from hwave.solver.so_bond import (SLOT_CONVENTION, SOBondCaps, plan_so_bond, so_bond_bubble,
                                  so_bond_dress, so_bond_topology, so_bond_vertex)
from hwave.solver.so_records import compile_interaction
from tests.so_bond_fixtures import production_gamma_b
from tests.so_bond_oracles import (CELL, FIXTURE_GRID_5X3, FIXTURE_I, _fermi, bubble_label_to_code,
                                   closed_row, expand_tables, file_gamma_label,
                                   free_green, lindhard_discrete, lindhard_matsubara,
                                   lindhard_static, matsubara_pair_kernel, mesh_negate,
                                   random_green, rel, ring_label)

NOCAP = SOBondCaps(None, None, None)
RING_TOL = 1e-12
FLAGS = SLOT_CONVENTION.label_view()


class Grid(NamedTuple):
    """A fixture (band and tables) on its mesh with a Matsubara count."""
    fixture: object
    nmat: int

    @property
    def geom(self):
        return self.fixture.geom

    @property
    def nd(self):
        return self.fixture.nd

    @property
    def beta(self):
        return self.fixture.beta

    @property
    def shape(self):
        return self.fixture.spatial_shape

    @property
    def band(self):
        return self.fixture.band


#: fixture (i) on the 5x1x1 chain, beta = 2, nmat = 16
FX_I_GRID = Grid(FIXTURE_I, 16)
#: the 5x3 grid with x, y and diagonal bonds, nmat = 8
GRID_5x3 = Grid(FIXTURE_GRID_5X3, 8)


def _view(geom):
    return BondSetView(BondTopology(geom.delta_r, geom.reverse, {}))


@functools.lru_cache(maxsize=None)
def _green(fx):
    """The free Green function of the band, ``(nmat, nvol, nd, nd)``."""
    return free_green(fx.band[0], fx.band[1], fx.beta, fx.nmat)


def _to_code(chi_label, fx, flags=FLAGS):
    """The bubble-rule image of a label array (static ``(nvol, NL, NL)`` or
    ``(nl, nvol, NL, NL)``)."""
    if chi_label.ndim == 3:
        return bubble_label_to_code(chi_label[None], fx.geom, fx.nd, flags)[0]
    return bubble_label_to_code(chi_label, fx.geom, fx.nd, flags)


def _flip_q(flags):
    """The reading with the label momentum flipped once more."""
    f = list(flags)
    f[2] = -f[2]
    return tuple(f)


def _prod_bubble(fx, freq):
    """The production bubble of the free Green function: ``(nvol, ND, ND)``
    (static) or ``(nmat, nvol, ND, ND)`` (every slice)."""
    g = _green(fx)[None].copy()
    nvol, nd, B = fx.fixture.nvol, fx.nd, fx.geom.B
    ND = B * nd * nd
    cons = ("chi_static",) if freq == "static" else ("chi_dynamic",)
    plan = plan_so_bond(nmat=fx.nmat, nvol=nvol, B=B, nd=nd, consumers=cons, caps=NOCAP)
    if freq == "static":
        sink = np.empty((nvol, ND, ND), dtype=np.complex128)
        out = so_bond_bubble(g, None, fx.beta, _view(fx.geom), fx.shape, plan=plan,
                             freq="static", sink=sink)
        assert out is sink
        return sink
    with flex_bond.BondBlockStore(fx.nmat, nvol, ND, nd * nd, ("chibar",)) as store:
        self_out = so_bond_bubble(g, None, fx.beta, _view(fx.geom), fx.shape, plan=plan,
                                  freq="all", sink=store)
        assert self_out is None
        return flex_bond.freq_batch_on(store, "chibar", 0, fx.nmat, np).copy()


@functools.lru_cache(maxsize=None)
def _gamma_b(fx):
    """The production ``Gamma_B`` of the fixture's tables on its own label
    channels (``fx.geom``)."""
    g = production_gamma_b(fx.fixture.tables_dict, fx.fixture.norb_phys, fx.geom, fx.shape)
    g.flags.writeable = False
    return g


@functools.lru_cache(maxsize=None)
def _gamma_label(fx):
    return file_gamma_label(expand_tables(fx.fixture.tables_dict, fx.fixture.norb_phys),
                            fx.geom, fx.nd)


def _prod_dress(fx, X):
    """The production dressing of a static ``(nvol, ND, ND)`` or a batched
    ``(nb, nvol, ND, ND)`` code bubble with the production ``Gamma_B``."""
    static = X.ndim == 3
    cons = ("chi_static",) if static else ("chi_dynamic",)
    plan = plan_so_bond(nmat=fx.nmat, nvol=fx.fixture.nvol, B=fx.geom.B, nd=fx.nd,
                        consumers=cons, caps=NOCAP)
    l0 = fx.nmat // 2 if static else fx.nmat // 2 + 1
    got = so_bond_dress(X[None] if static else X, _gamma_b(fx), plan=plan, l0=l0,
                        spatial_shape=fx.shape)
    return got[0] if static else got


@functools.lru_cache(maxsize=None)
def _chibar_label(fx, ls):
    out = (lindhard_static(*fx.band, fx.geom) if ls is None
           else lindhard_matsubara(*fx.band, fx.geom, ls))
    out.flags.writeable = False
    return out


def _chibar_code(fx):
    return _to_code(_chibar_label(fx, None), fx)


def _sector(arr, s):
    """The block of slots ``(m; x, y)`` with ``S_z`` transfer
    ``s_y - s_x = s`` (spin of the spin-block mode ``X = s * norb + a``,
    norb_phys = 1) in both the row and the column."""
    ND = arr.shape[-1]
    idx = np.arange(ND)
    x, y = (idx // 2) % 2, idx % 2
    sel = np.nonzero((y - x) == s)[0]
    return arr[..., sel[:, None], sel[None, :]]


# ---------------------------------------------------------------------------
# Oracles used below
# ---------------------------------------------------------------------------

class TestOracles(unittest.TestCase):
    def test_matsubara_kernel_against_the_sum(self):
        beta, xa, xb = 2.0, 0.37, -0.81
        n = np.arange(-200000, 200000)
        w = (2 * n + 1) * np.pi / beta
        for m in (1, -2, 3):
            nu = 2 * np.pi * m / beta
            brute = -(1.0 / beta) * np.sum(1.0 / ((1j * w - xa) * (1j * w - 1j * nu - xb)))
            want = complex(matsubara_pair_kernel(xa, xb, beta, nu))
            self.assertLessEqual(abs(brute - want) / abs(want), 1e-5, m)
        self.assertEqual(complex(matsubara_pair_kernel(xa, xb, beta, 0.0)),
                         complex((_fermi(xb, beta) - _fermi(xa, beta)) / (xa - xb)))

    def test_matsubara_zero_is_static(self):
        fx = FX_I_GRID
        self.assertLessEqual(rel(_chibar_label(fx, (0,))[0], _chibar_label(fx, None)), 1e-14)

    def test_ring_label_is_the_definition(self):
        fx = FX_I_GRID
        chi, gam = _chibar_label(fx, None), _gamma_label(fx)
        got = ring_label(chi, gam, fx.shape)
        neg = mesh_negate(fx.shape)
        for q in range(fx.fixture.nvol):
            want = np.linalg.inv(np.eye(chi.shape[1]) + chi[q] @ gam[neg[q]]) @ chi[q]
            self.assertLessEqual(rel(got[q], want), 1e-13)


# ---------------------------------------------------------------------------
# Entry checks
# ---------------------------------------------------------------------------

G_2x1 = random_green(5, 2, 2, 2)[None].copy()


def _plan_for(view, shape):
    return plan_so_bond(nmat=2, nvol=int(np.prod(shape)), B=len(view.delta_r), nd=2,
                        consumers=("chi_static",), caps=NOCAP)


class TestBubbleEntry(unittest.TestCase):
    def test_mesh_coincident_refused_before_allocation(self):
        view = SimpleNamespace(delta_r=((0, 0, 0), (1, 0, 0), (-1, 0, 0)), reverse=(0, 2, 1),
                               n_channels=3)
        plan = _plan_for(view, (2, 1, 1))
        tracemalloc.start()
        self.addCleanup(tracemalloc.stop)
        with self.assertRaises(ValueError) as cm:
            so_bond_bubble(G_2x1, None, 2.0, view, (2, 1, 1), plan=plan, freq="static",
                           sink=None)
        self.assertEqual(
            str(cm.exception),
            "validate_topology_against_mesh: channels 1 (delta_r=(1, 0, 0)) and 2 "
            "(delta_r=(-1, 0, 0)) collide on the mesh [2, 1, 1] -- both wrap to the same "
            "displacement (1, 0, 0) modulo the mesh (an alias such as R vs R+k*L on some "
            "axis, or the even-mesh self-reversal +L/2 == -L/2); a topology whose channels "
            "are not INJECTIVE modulo the mesh cannot be used at this mesh resolution.")
        self.assertLess(tracemalloc.get_traced_memory()[1], 100 * 1024)

    def test_plan_mismatch_refused_before_work(self):
        view = SimpleNamespace(delta_r=((0, 0, 0), (1, 0, 0), (-1, 0, 0)), reverse=(0, 2, 1),
                               n_channels=3)
        plan = plan_so_bond(nmat=4, nvol=5, B=3, nd=2, consumers=("chi_static",), caps=NOCAP)
        g = random_green(5, 2, 5, 2)[None].copy()
        sink = np.zeros((5, 12, 12), dtype=np.complex128)
        with self.assertRaises(ValueError) as cm:
            so_bond_bubble(g, None, 2.0, view, (5, 1, 1), plan=plan, freq="static", sink=sink)
        self.assertEqual(str(cm.exception),
                         "so_bond_bubble: the admission plan has nmat = 4 but the call has "
                         "nmat = 2 ({})".format(plan.reason))
        self.assertEqual(np.count_nonzero(sink), 0)

    def test_sink_refusals(self):
        view = SimpleNamespace(delta_r=((0, 0, 0), (1, 0, 0), (-1, 0, 0)), reverse=(0, 2, 1),
                               n_channels=3)
        g = random_green(5, 2, 5, 2)[None].copy()
        plan = plan_so_bond(nmat=2, nvol=5, B=3, nd=2, consumers=("chi_static",), caps=NOCAP)
        with self.assertRaises(ValueError) as cm:
            so_bond_bubble(g, None, 2.0, view, (5, 1, 1), plan=plan, freq="static",
                           sink=np.zeros((5, 12, 11), dtype=np.complex128))
        self.assertEqual(str(cm.exception),
                         "so_bond_bubble: the static sink must be a complex128 array of shape "
                         "(5, 12, 12), got complex128 (5, 12, 11)")
        dplan = plan_so_bond(nmat=2, nvol=5, B=3, nd=2, consumers=("chi_dynamic",), caps=NOCAP)
        with flex_bond.BondBlockStore(2, 5, 12, 4, ("chi",)) as store:
            with self.assertRaises(ValueError) as cm:
                so_bond_bubble(g, None, 2.0, view, (5, 1, 1), plan=dplan, freq="all",
                               sink=store)
        self.assertEqual(str(cm.exception),
                         "so_bond_bubble: the store slots ('chi',) do not include 'chibar' "
                         "({})".format(dplan.reason))

    def test_non_finite_static_bubble_not_published(self):
        """A non-finite Green function gives a non-finite static bubble; it
        is refused with NonFiniteError (as the every-slice path refuses a
        non-finite block) and the sink is left untouched."""
        view = SimpleNamespace(delta_r=((0, 0, 0), (1, 0, 0), (-1, 0, 0)), reverse=(0, 2, 1),
                               n_channels=3)
        plan = plan_so_bond(nmat=2, nvol=5, B=3, nd=2, consumers=("chi_static",), caps=NOCAP)
        for bad in (np.nan, np.inf):
            with self.subTest(bad=bad):
                g = random_green(5, 2, 5, 2)[None].copy()
                g[0, 1, 2, 0, 1] = bad
                sink = np.zeros((5, 12, 12), dtype=np.complex128)
                with self.assertRaises(NonFiniteError) as cm, np.errstate(invalid="ignore"):
                    so_bond_bubble(g, None, 2.0, view, (5, 1, 1), plan=plan, freq="static",
                                   sink=sink)
                self.assertEqual(str(cm.exception),
                                 "so_bond_bubble: non-finite static bond bubble")
                self.assertEqual(np.count_nonzero(sink), 0)


class TestDressFinite(unittest.TestCase):
    def test_non_finite_inputs_refused_before_dressing(self):
        """A NaN or inf in the bubble batch or the vertex is refused with
        NonFiniteError naming the argument, with or without the
        conditioning guard, before the dressing runs."""
        plan = plan_so_bond(nmat=8, nvol=5, B=3, nd=2, consumers=("chi_static",), caps=NOCAP)
        for arg in ("X", "gamma_b"):
            for bad in (np.nan, np.inf):
                for cond_tol in (None, 1e-12):
                    with self.subTest(arg=arg, bad=bad, cond_tol=cond_tol):
                        X = np.zeros((1, 5, 12, 12), complex)
                        G = np.zeros((5, 12, 12), complex)
                        (X[0] if arg == "X" else G)[2, 3, 4] = bad
                        with self.assertRaises(NonFiniteError) as cm:
                            so_bond_dress(X, G, plan=plan, l0=4, spatial_shape=(5, 1, 1),
                                          cond_tol=cond_tol)
                        self.assertEqual(str(cm.exception),
                                         "so_bond_dress: {} has non-finite entries".format(arg))


class TestDressCapacity(unittest.TestCase):
    def test_multi_slice_under_static_plan_refused(self):
        plan = plan_so_bond(nmat=8, nvol=5, B=3, nd=2, consumers=("chi_static",), caps=NOCAP)
        X = np.zeros((4, 5, 12, 12), complex)
        G = np.zeros((5, 12, 12), complex)
        tracemalloc.start()
        self.addCleanup(tracemalloc.stop)
        with self.assertRaises(ValueError) as cm:
            so_bond_dress(X, G, plan=plan, l0=0, spatial_shape=(5, 1, 1))
        self.assertEqual(str(cm.exception),
                         "so_bond_dress: a batch of 4 slices exceeds the plan's dress_capacity "
                         "= 1 ({})".format(plan.reason))
        self.assertLess(tracemalloc.get_traced_memory()[1], 100 * 1024)

    def test_batch_shape_refused(self):
        plan = plan_so_bond(nmat=8, nvol=5, B=3, nd=2, consumers=("chi_static",), caps=NOCAP)
        with self.assertRaises(ValueError) as cm:
            so_bond_dress(np.zeros((5, 12, 12), complex), np.zeros((5, 12, 12), complex),
                          plan=plan, l0=0, spatial_shape=(5, 1, 1))
        self.assertEqual(str(cm.exception),
                         "so_bond_dress: X must be a batch (nb, nvol, ND, ND), got shape "
                         "(5, 12, 12)")


# ---------------------------------------------------------------------------
# The production bubble and dressing against the label ring
# ---------------------------------------------------------------------------

class TestBubblePart(unittest.TestCase):
    """Gate A1-chi ring, bubble part: the production bubble of the free
    Green function equals the bubble-rule image of the discrete label sum,
    on every bosonic slice and, for the static bubble, on ``nu = 0``."""

    def test_production_bubble_vs_discrete_sum(self):
        for fx in (FX_I_GRID, GRID_5x3):
            want = _to_code(lindhard_discrete(_green(fx), fx.beta, fx.geom), fx)
            for freq in ("static", "all"):
                with self.subTest(L=fx.shape, freq=freq):
                    got = _prod_bubble(fx, freq)
                    ref = want[fx.nmat // 2] if freq == "static" else want
                    self.assertLessEqual(rel(got, ref), RING_TOL)


class TestDressingPart(unittest.TestCase):
    """Gate A1-chi ring, dressing part: the production dressing of the
    bubble-rule image of the analytic label bubble, with the production
    ``Gamma_B``, equals the bubble-rule image of the label ring built with
    the independently filed label ``Gamma``; reading the ring with one more
    momentum flip misses."""

    def test_static_and_matsubara(self):
        for fx in (FX_I_GRID, GRID_5x3):
            for ls in (None, (1, 3)):
                with self.subTest(L=fx.shape, ls=ls):
                    chibar_L = _chibar_label(fx, ls)
                    ring = ring_label(chibar_L, _gamma_label(fx), fx.shape)
                    X = _to_code(chibar_L, fx)
                    got = _prod_dress(fx, X)
                    want = _to_code(ring, fx)
                    self.assertLessEqual(rel(got, want), RING_TOL)
                    # the dressing is not the identity at this coupling
                    self.assertGreaterEqual(rel(got, X), 1e-2)
                    # one more momentum flip of the ring misses
                    self.assertGreaterEqual(rel(got, _to_code(ring, fx, _flip_q(FLAGS))), 1e-3)

    def test_longitudinal_and_transverse_slots_nonvacuous(self):
        got = _prod_dress(FX_I_GRID, _chibar_code(FX_I_GRID))
        for s in (0, +1, -1):
            with self.subTest(sector=s):
                self.assertGreaterEqual(np.abs(_sector(got, s)).max(), 1e-3 * np.abs(got).max())

    def test_one_band_hubbard_channels(self):
        c, zz, pm = _hubbard_channels(_dress_scalar(0.31, 1.1))
        self.assertAlmostEqual(c, 0.31 / (1 + 1.1 * 0.31), places=13)
        self.assertAlmostEqual(zz, 0.31 / (1 - 1.1 * 0.31), places=13)
        self.assertAlmostEqual(pm, zz, places=13)

    def test_all_declared_zero(self):
        tables = _all_zero_tables()
        with self.assertLogs("hwave.solver.so_records", "INFO") as log:
            c = compile_interaction(tables, 1)
        self.assertEqual(len(c.records), 0)
        self.assertEqual([r for r in log.output if "bond shell" in r],
                         ["INFO:hwave.solver.so_records:so_records: bond shell(s) "
                          + _shell_list(tables) + " carry only couplings with "
                          "|v| <= 1e-14 * max|v|"])
        topo = so_bond_topology(c, CELL)
        self.assertEqual(len(topo.delta_r), 5)                 # every declared shell is kept
        plan = plan_so_bond(nmat=4, nvol=5, B=5, nd=2, consumers=("chi_static",), caps=NOCAP)
        g = so_bond_vertex(c, topo, (5, 1, 1), plan=plan)
        self.assertEqual(np.count_nonzero(g), 0)
        X = np.empty((5, 20, 20), dtype=np.complex128)
        so_bond_bubble(random_green(9, 4, 5, 2)[None].copy(), None, 2.0, topo, (5, 1, 1),
                       plan=plan, freq="static", sink=X)
        self.assertGreater(np.count_nonzero(X), 0)
        got = so_bond_dress(X[None], g, plan=plan, l0=2, spatial_shape=(5, 1, 1))[0]
        self.assertTrue(np.array_equal(got, X))


def _dress_scalar(chi0, U):
    """The production dressing of the one-band Hubbard bubble ``chibar_B =
    chi0 * 1`` (a spin-diagonal Green function; on-site channel only, one
    momentum) with the vertex of ``CoulombIntra = U``. Returns ``(4, 4)``."""
    c = compile_interaction({"CoulombIntra": {((0, 0, 0), (0, 0)): U}}, 1)
    topo = so_bond_topology(c, CELL)
    plan = plan_so_bond(nmat=2, nvol=1, B=1, nd=2, consumers=("chi_static",), caps=NOCAP)
    g = so_bond_vertex(c, topo, (1, 1, 1), plan=plan)
    X = chi0 * np.eye(4, dtype=np.complex128)[None, None]
    return so_bond_dress(X, g, plan=plan, l0=1, spatial_shape=(1, 1, 1))[0, 0]


def _hubbard_channels(chi):
    """``(charge, zz, +-)`` of a dressed one-band ``(4, 4)`` code matrix:
    slots ``(x, y)`` at ``2 x + y``; charge and zz on ``(up up) +- (dn dn)``
    (normalized), ``+-`` on the slot ``(up, dn)``."""
    u = np.array([1.0, 0.0, 0.0, 1.0]) / np.sqrt(2.0)
    v = np.array([1.0, 0.0, 0.0, -1.0]) / np.sqrt(2.0)
    return complex(u @ chi @ u), complex(v @ chi @ v), complex(chi[1, 1])


def _all_zero_tables():
    return {"CoulombIntra": {((0, 0, 0), (0, 0)): 0.0},
            "CoulombInter": {**closed_row(1, 0, 0, 0.0), **closed_row(2, 0, 0, 0.0)},
            "Exchange": closed_row(1, 0, 0, 0.0)}


def _shell_list(tables):
    """Every declared off-site shell ``(-d, d)`` of the tables, sorted, in the
    format of the INFO line."""
    shells = set()
    for tbl in tables.values():
        for R, _ in tbl:
            if any(R):
                shells.add(tuple(sorted({tuple(R), tuple(-x for x in R)})))
    return str(sorted(shells))


if __name__ == "__main__":
    unittest.main()
