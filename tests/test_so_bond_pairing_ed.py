"""Exact diagonalization of the spinful bond pairing vertex at second order.

The comparison of ``tests/test_so_pairing_ed.py`` (its pair basis, ED
observables, ladder, stencils and error measure are imported from there)
with bond channels and off-site interactions. The observable is the static
pair susceptibility of the spin-orbital pair operators of a chain,
``chi_{(A, B, d'), (E, F, d)} = <Delta_(A, B, d'); Delta_(E, F, d)^dagger>``
(``TotalNED.pair_correlator``), and the production side is the
Bethe-Salpeter ladder through second order, ``chi_prod eta = 2 X[eta + P K
eta + P K P K eta]``, with ``G`` the exact interacting Green function of the
same ED (``TotalNED.green``) and ``K`` the bond pairing kernel
(``BondPairKernel``, ``B = 3`` channels ``0, +-1``, ``nfreq = nmat``, no
instantaneous vertex) of the production chain at every bosonic slice ``l``:

1. ``chibar_B(l)`` = ``so_bond_bubble(freq="all")`` of ``G``;
2. ``chi_B = so_bond_dress(chibar_B[l0:l1], lam Gamma_B)`` per frequency
   batch of ``plan.nb`` slices;
3. ``V(l) = V_bare + V_fluct(l)`` from
   ``so_bond_pairing_vertex_parts(lam Gamma_B, chi_B(l))``:
   ``V_bare = 0.25 (Gamma_B + Gamma_B^P)`` on every channel block and
   ``V_fluct = -Gamma_B chi_B Gamma_B``.

Every interaction is scaled by ``lam``; the dressing changes ``V_fluct`` at
O(lam^3) only, so the ``lam^1`` and ``lam^2`` coefficients (five-point
stencils, ``h = 0.05``) of production and ED must agree up to the Matsubara
truncation (Richardson over ``nmat = 128, 256``), with the tolerances of the
on-site comparison (2e-5 and 2e-4). Measured: 6.8e-6 (``c1``) and 5.8e-5
(``c2``).

Model: the two-orbital three-site chain with complex spin-flip hopping of
``tests/test_so_pairing_ed.py`` (``so_bond_oracles.FIXTURE_II``'s band),
every on-site interaction of that test (``TYPES``) plus the nearest-neighbour
inter-orbital rows of every off-site type (``X_01(+1) = v``, ``X_01(-1) =
0.6 v``: mirror-asymmetric, complex PairHop, PairLift), the ED terms from
``so_bond_oracles.build_terms``.

The comparison is two-sided: the fluctuation term read in the opposite
orientation (``V[(rev m'; d, c), (rev m; b, a)](-q)`` with the phase
``e^{+iq.(dr_m' - dr_m)}``, i.e. the exchange image ``gamma_p`` of
``V_fluct``) and the ladder without the fluctuation term each miss the
``lam^2`` coefficient by more than 1e-2 (measured through the same
Richardson extrapolation as the implemented reading, 5.8e-5: 4.4e-2 and
1.2e-1; a q-flipped ``chi_B`` would miss by 2.5e-2). The off-site rows
carry more than a tenth of the ED ``lam^2`` coefficient (measured 0.48,
against the on-site-only model). With the on-site rows alone (all three
channels kept) and the bubble undressed -- the one-channel production does
not dress, and the dressing differs from it at O(lam^3), which the
coefficients do not see but the susceptibilities do -- the bond chain
reproduces the one-channel production of ``tests/test_so_pairing_ed.py`` on
the same ``G`` to 1e-12 (measured 6.5e-17).

The fast smoke runs the same chain on a one-orbital three-site chain
(``FIXTURE_I``'s band with spin-flip hopping; CoulombInter, Exchange and
complex PairHop on the nearest-neighbour bond; ``nmat = 64, 128``) and is
two-sided too: the opposite reading, no fluctuation term and a q-flipped
``chi_B``, each Richardson-extrapolated as the implemented reading is, miss
the ``lam^2`` coefficient by more than three times its deviation.

Smoke tolerance. ``SMOKE_TOL = min(SMOKE_TOL_CEILING, 3 x the measured c1
deviation)``. The ceiling 2e-2 is an a-priori bound fixed before the
measurement (the same number serves as the floor of the spin-orbital FLEX
smoke tolerance of ``tests/test_flex_spin_orbital_ed.py``, a different
check); the measured deviations are 2.7e-5 (``c1``) and 9.8e-5 (``c2``), and
the three wrong readings miss ``c2`` by 4.8e-2, 1.5e-1 and 3.1e-2.
"""
import functools
import unittest

import numpy as np

import tests.test_so_pairing_ed as E
from hwave.solver import bond_channels as bc
from hwave.solver import flex_bond
from hwave.solver.eliashberg_dynamic import calc_g2_dynamic
from hwave.solver.so_bond import (SOBondCaps, gamma_p, plan_so_bond, slot_geometry,
                                  so_bond_bubble, so_bond_dress, so_bond_pairing_vertex_parts,
                                  so_bond_topology)
from hwave.solver.so_records import compile_interaction
from tests.ed_oracle_util import EDFixture, TotalNED
from tests.heavy_tests import heavy
from tests.so_bond_fixtures import green_kw, onsite_decl, pair_kernel, production_gamma_b
from tests.so_bond_oracles import (BASE, CELL, FIXTURE_I, FIXTURE_II, FixtureSpec, build_terms,
                                   closed_row, freeze_tables, merge_tables, mesh_negate, rel)

NOCAP = SOBondCaps(None, None, None)
#: the readings of the fluctuation term: as implemented, the opposite
#: orientation, none, and the dressed susceptibility read at -q
READINGS = ("direct", "opposite", "none", "chi_q_flipped")
SMOKE_TOL_CEILING = 2e-2
#: 3 x the measured smoke c1 deviation (2.7e-5), below the ceiling
SMOKE_TOL = min(SMOKE_TOL_CEILING, 8.1e-5)


#: the on-site interaction of ``tests/test_so_pairing_ed.py`` (``TYPES``)
ONSITE = merge_tables(*(onsite_decl(t, v) for t, v in E.TYPES))
#: the nearest-neighbour rows of every off-site type (fixture (ii))
NN = {t: tbl for t, tbl in FIXTURE_II.tables_dict.items() if t != "CoulombIntra"}
#: the two-orbital chain with the on-site and the nearest-neighbour rows
BOND_MODEL = FIXTURE_II._replace(tables=freeze_tables(merge_tables(ONSITE, NN)))
#: the same chain with the on-site rows alone
ONSITE_MODEL = FIXTURE_II._replace(tables=freeze_tables(ONSITE))
#: the one-orbital smoke chain: CoulombInter, Exchange and complex PairHop on
#: the nearest-neighbour bond, spin-flip hopping
SMOKE_MODEL = FixtureSpec.make(
    3, 1, FIXTURE_I.beta, FIXTURE_I.mu, FIXTURE_I.transfer_dict,
    {t: closed_row(1, 0, 0, BASE[t]) for t in ("CoulombInter", "Exchange", "PairHop")})
#: the Matsubara counts of each model (Richardson pair)
NMATS = {BOND_MODEL: E.NMATS, SMOKE_MODEL: (64, 128)}


def _pairs(model):
    return E._pair_basis(model.nd, model.L)


# ---------------------------------------------------------------------------
# exact diagonalization
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=None)
def _ed_series(model, nmat):
    """``((chi, G), ...)`` over ``E.LAMS``: the pair susceptibility and the
    Green function (on ``nmat`` points; none for ``nmat = 0``) of ``model``
    with every interaction scaled by ``lam`` (cached by content)."""
    fx = EDFixture(model.L, model.norb_phys, {}, [0.0] * model.norb_phys, 1.0 / model.beta,
                   model.mu)
    terms = build_terms(model.tables_dict, model.L, model.norb_phys)
    h1 = model.h1
    return tuple(E._ed_observables(TotalNED(fx, [t[:4] + (lam * t[4],) for t in terms], h1=h1),
                                   _pairs(model), nmat)
                 for lam in E.LAMS)


def _ed_chi(model):
    nmat = max(NMATS[model]) if model in NMATS else 0
    return tuple(chi for chi, _ in _ed_series(model, nmat))


def _ed_green(model, nmat):
    return tuple(E._green_sc(G, nmat, model.nd, model.L)
                 for _, G in _ed_series(model, max(NMATS[model])))


# ---------------------------------------------------------------------------
# production chain
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=None)
def _topology(model):
    return so_bond_topology(compile_interaction(model.tables_dict, model.norb_phys), CELL)


@functools.lru_cache(maxsize=None)
def _gamma_b(model, tables):
    """``Gamma_B`` of ``tables`` (frozen) on the channels of ``model``
    (read-only)."""
    out = production_gamma_b({t: {(R, ab): v for R, ab, v in tbl} for t, tbl in tables},
                             model.norb_phys, _topology(model), model.spatial_shape)
    out.flags.writeable = False
    return out


def _vertex_blocks(model, g, gamma, nmat, reading, dress):
    """``V(l) = V_bare + V_fluct(l)`` ``(nmat, nvol, ND, ND)`` of the
    production chain on the sc-layout Green function ``g``; ``reading``
    selects the fluctuation term (:data:`READINGS`); ``dress=False`` feeds
    the bare bubble in place of ``chi_B``."""
    shape, nvol, nd = model.spatial_shape, model.nvol, model.nd
    topo = _topology(model)
    geom = slot_geometry(topo, shape)
    ND = geom.B * nd * nd
    plan = plan_so_bond(nmat=nmat, nvol=nvol, B=geom.B, nd=nd,
                        consumers=("chi_dynamic", "pairing_static"), caps=NOCAP)
    neg = mesh_negate(shape)
    blocks = np.empty((nmat, nvol, ND, ND), dtype=np.complex128)
    with flex_bond.BondBlockStore(nmat, nvol, ND, nd * nd, ("chibar",)) as store:
        so_bond_bubble(green_kw(g), None, model.beta, topo, shape, plan=plan, freq="all",
                       sink=store)
        for l0 in range(0, nmat, plan.nb):
            l1 = min(l0 + plan.nb, nmat)
            chibar = flex_bond.freq_batch_on(store, "chibar", l0, l1, np)
            chi = (so_bond_dress(chibar, gamma, plan=plan, l0=l0, spatial_shape=shape)
                   if dress else np.array(chibar))
            for i in range(l1 - l0):
                c = chi[i][neg] if reading == "chi_q_flipped" else chi[i]
                vb, vf = so_bond_pairing_vertex_parts(gamma, np.ascontiguousarray(c), geom,
                                                      plan=plan)
                if reading == "opposite":
                    vf = gamma_p(vf, geom, out=np.empty_like(vf))
                elif reading == "none":
                    vf = 0.0
                blocks[l0 + i] = vb + vf
    return blocks


def _chain_chi(model, g, lam, nmat, reading="direct", tables=None, dress=True):
    """The production ladder through ``K K`` (pair basis of ``model``) at
    coupling ``lam`` with the vertex of ``tables`` (default: the model's)."""
    gamma = lam * _gamma_b(model, model.tables if tables is None else tables)
    blocks = _vertex_blocks(model, g, gamma, nmat, reading, dress)
    view = bc.BondSetView(_topology(model))
    K = pair_kernel(blocks, calc_g2_dynamic(g, model.beta), view, model.nd,
                     model.spatial_shape, model.beta)
    return E._ladder_chi(K, g, nmat, _pairs(model), model.beta)


@functools.lru_cache(maxsize=None)
def _production(model, nmat, reading="direct"):
    return tuple(_chain_chi(model, g, lam, nmat, reading)
                 for g, lam in zip(_ed_green(model, nmat), E.LAMS))


def _errors_vs_ed(fluct, model=BOND_MODEL):
    """``{order: deviation}`` of the Richardson production coefficients
    (``2 chi(nmat_hi) - chi(nmat_lo)``) of the reading ``fluct`` from the ED
    ones."""
    lo, hi = (_production(model, n, fluct) for n in NMATS[model])
    return E._errors(_ed_chi(model), tuple(2 * b - a for a, b in zip(lo, hi)), _pairs(model),
                     model.L)


def _coefficients(smoke=False):
    """``(dev_c1, dev_c2)`` of the implemented reading (:func:`_errors_vs_ed`,
    the path of every control reading)."""
    dev = _errors_vs_ed("direct", SMOKE_MODEL if smoke else BOND_MODEL)
    return dev[1], dev[2]


def _smoke_errors(reading):
    return _errors_vs_ed(reading, SMOKE_MODEL)


def _offsite_share_of_c2():
    """``|c2(model) - c2(on-site rows only)| / |c2(model)|`` of the ED, on
    the antisymmetric pair space."""
    P = E._antisymmetric_projector(_pairs(BOND_MODEL), BOND_MODEL.L)
    full = P @ E._stencil(_ed_chi(BOND_MODEL))[1] @ P
    onsite = P @ E._stencil(_ed_chi(ONSITE_MODEL))[1] @ P
    return np.linalg.norm(full - onsite) / np.linalg.norm(full)


def _bond_production_onsite_only(nmat):
    """The three-channel chain with the on-site rows alone and the bubble
    undressed, on the Green functions of :data:`BOND_MODEL`."""
    return tuple(_chain_chi(BOND_MODEL, g, lam, nmat, tables=ONSITE_MODEL.tables, dress=False)
                 for g, lam in zip(_ed_green(BOND_MODEL, nmat), E.LAMS))


def _harness_production(nmat):
    """The one-channel production of ``tests/test_so_pairing_ed.py`` (its
    on-site ``TYPES``, undressed bubble) on the same Green functions."""
    return tuple(E._production_chi(g, lam, nmat)
                 for g, lam in zip(_ed_green(BOND_MODEL, nmat), E.LAMS))


class TestSOBondPairingAgainstED(unittest.TestCase):
    """Gate A1-E: the spinful bond pairing vertex (bubble, dressing, static
    pairing parts, bond kernel) against exact diagonalization through
    second order, two-sided, with a fast two-sided smoke."""

    def _spin_mixing(self, model):
        """The spin-flip hopping is load bearing: the Green function has
        sizeable spin-off-diagonal blocks."""
        g = _ed_green(model, NMATS[model][0])[2]
        n = model.norb_phys
        self.assertGreater(np.linalg.norm(g[:n, n:]), 1e-1 * np.linalg.norm(g))

    @heavy
    def test_first_and_second_order_two_sided(self):
        self._spin_mixing(BOND_MODEL)
        c1, c2 = _coefficients()
        self.assertLess(c1, E.TOL_C1, "order 1: {:.3e}".format(c1))
        self.assertLess(c2, E.TOL_C2, "order 2: {:.3e}".format(c2))
        share = _offsite_share_of_c2()
        self.assertGreater(share, 1e-1, "off-site share of c2: {:.3e}".format(share))
        # the controls go through the same Richardson path, which the
        # implemented reading passes well below the miss threshold
        direct = _errors_vs_ed("direct")[2]
        self.assertLess(direct, E.MISS / 10, "direct: {:.3e}".format(direct))
        for reading in ("opposite", "none"):
            dev = _errors_vs_ed(reading)[2]
            self.assertGreater(dev, E.MISS, "{}: {:.3e}".format(reading, dev))
        dev = rel(np.array(_bond_production_onsite_only(NMATS[BOND_MODEL][0])),
                   np.array(_harness_production(NMATS[BOND_MODEL][0])))
        self.assertLessEqual(dev, 1e-12, "on-site rows vs one channel: {:.3e}".format(dev))

    def test_smoke_two_sided(self):
        self._spin_mixing(SMOKE_MODEL)
        dev = _coefficients(smoke=True)
        self.assertLess(dev[0], SMOKE_TOL, "order 1: {:.3e}".format(dev[0]))
        self.assertLess(dev[1], 10 * SMOKE_TOL, "order 2: {:.3e}".format(dev[1]))
        # the controls go through the same Richardson path as the implemented
        # reading (whose order-2 deviation is dev[1]) and miss by more than 3x it
        threshold = 3 * dev[1]
        for reading in ("opposite", "none", "chi_q_flipped"):
            miss = _smoke_errors(reading)[2]
            self.assertGreater(miss, threshold, "{}: {:.3e}".format(reading, miss))


if __name__ == "__main__":
    unittest.main()
