"""GitHub issue #181, Tier 3 Phase A -- gates G2 and G3: the bond-resolved
LONGITUDINAL (spin/charge) vertex against exact diagonalization, first
order in the coupling, with the #151 bond-channel oracle apparatus
(``tests/test_bond_vs_ed_oracle.py``: HF-subtracted bond correlators,
Richardson at CAMPAIGN_V1 and CAMPAIGN_V1/2, the derived ED-to-solver
index map, ``adjudicate_granule``'s tolerance rule).

Solver side: ``dS``/``dC`` = ``build_sc_bond_channel`` at unit declared
amplitude minus the same declared topology at zero (same ND on both
sides), on the channel-0 blocks of the production locality split
(``offsite.sc_matrices_from_split``); ``D_pred = (+chi_bar dS chi_bar,
-chi_bar dC chi_bar)`` (``pred_first_order``), exact in the coupling.
ED side: ``X_S = X^uu - X^ud``, ``X_C = X^uu + X^ud`` of the HF-subtracted
bond correlator, mapped into the solver's ``(q, ND, ND)`` frame.

Verdicts (spec "Gates"), recorded from the run that froze this module
(2026-09-05; NMAT = 1024, CAMPAIGN_V1 = 0.005; tol = 10 * max(delta_rich,
delta_nmat)):

    fx5 (L=5, norb=1, a == b bonds, shells +-1 and +-2 declared)
      CoulombInter R=1   S PASS (tol 7.3e-4, signal 0.135)  C PASS (2.4e-3, 0.437)
      CoulombInter R=2   S PASS (7.3e-4, 0.135)             C PASS (2.0e-3, 0.275)
      Hund         R=1   S PASS (1.2e-3, 0.219)             C PASS (1.2e-3, 0.219)
      Ising        R=1   S PASS (2.4e-3, 0.437)             C PASS (7.3e-4, 0.135)
      G3 (bond blocks zeroed): CoulombInter/Hund/Ising all FAIL with
        95-125 failing cells each -- the bond blocks are load-bearing.
    fx3 (L=3, norb=2, inter-orbital (0,1) bond, shell +-1)
      CoulombInter       S PASS (7.3e-4, 0.138)             C PASS (1.4e-3, 0.246)
      Hund               S PASS (7.3e-4, 0.137)             C PASS (7.3e-4, 0.137)
      Ising              S PASS (1.4e-3, 0.246)             C PASS (7.3e-4, 0.138)

The Exchange direction is ORACLE-ONLY -- production refuses an off-site
Exchange declaration -- and its verdict is recorded, not asserted: the
(+1, +1) bond-diagonal rule for a REAL off-site J is PASS on both fx5
(a == b) and fx3 (a != b), S and C (same tolerances as CoulombInter).
That is the recorded input for the follow-up "Exchange promotion" of
issue #181; it does not change Tier 2's finding (no q-ONLY content).
"""
import functools
import unittest

import numpy as np

from hwave.solver import bond_channels as bc
from hwave.solver import bubble
from hwave.solver.offsite import split_locality, sc_matrices_from_split
from tests import ed_oracle_util
from tests.heavy_tests import heavy
from tests.test_bond_transverse_ed import (
    _terms_exchange_offsite, _terms_hund_offsite, _terms_ising_offsite)
from tests.test_bond_vs_ed_oracle import (
    CAMPAIGN_V1, NMAT, ed_to_solver_bond_map, free_green, fx3, fx5,
    pred_first_order)

_ALL = ("CoulombInter", "Hund", "Ising", "Exchange")


def _decls(t, R, a, b, v, shells):
    """Two-sided declaration of type ``t`` at amplitude ``v`` on a FIXED
    topology: every shell in ``shells`` is declared (zero where it is
    not the direction), so ND is the same at every amplitude."""
    d = {}
    for S in shells:
        val = v if S == R else 0.0
        d[((S, 0, 0), (a, b))] = val
        d[((-S, 0, 0), (b, a))] = np.conj(val)
    return {t: d}


def _terms(t, fx, a, b, R, v):
    if t == "CoulombInter":
        return ed_oracle_util.canonical_density_terms(fx, [(a, b, R, v)])
    return {"Hund": _terms_hund_offsite, "Ising": _terms_ising_offsite,
            "Exchange": _terms_exchange_offsite}[t](fx, a, b, R, v)


class _Side:
    """Solver-side objects for one coupling direction on fixture ``fx``."""

    def __init__(self, fx, t, a, b, R, shells):
        self.fx, self.t, self.a, self.b, self.R = fx, t, a, b, R
        self.shells = tuple(shells)

    def topo(self, v):
        return bc.resolve_bond_topology(
            _decls(self.t, self.R, self.a, self.b, v, self.shells),
            np.eye(3), self.fx.norb, active_types=_ALL)

    def w0(self, v):
        decl = _decls(self.t, self.R, self.a, self.b, v, self.shells)

        class H:
            param_ham = decl
            param_ham_orig = None

            def _reshape_interaction(self, tbl, so):
                return tbl

        class Lat:
            subshape = (1, 1, 1)

        S0, C0 = sc_matrices_from_split(
            split_locality(H(), Lat()), ("CoulombInter", "Hund", "Ising"),
            self.fx.norb, self.fx.L, 1, 1)
        nd = self.fx.norb ** 2
        return S0.reshape(self.fx.L, nd, nd), C0.reshape(self.fx.L, nd, nd)

    def dW(self, bond_blocks=True):
        S1, C1 = self.w0(1.0)
        S0, C0 = self.w0(0.0)
        t1, t0 = self.topo(1.0), self.topo(0.0)
        types = (self.t,) if bond_blocks else ()
        dS = (bc.build_sc_bond_channel(t1, S1, "S", types=types)
              - bc.build_sc_bond_channel(t0, S0, "S", types=types))
        dC = (bc.build_sc_bond_channel(t1, C1, "C", types=types)
              - bc.build_sc_bond_channel(t0, C0, "C", types=types))
        return dS, dC

    def view(self):
        return bc.BondSetView(self.topo(0.0))

    def chibar(self, nmat):
        return _chibar(self.fx, self.shells, nmat)

    def channels_and_map(self):
        view = self.view()
        norb = self.fx.norb
        channels = []
        for m in range(view.n_channels):
            R = view.delta_r[m][0]
            for l1 in range(norb):
                for l2 in range(norb):
                    channels.append((R,) if norb == 1 else (R, l1, l2))
        return channels, ed_to_solver_bond_map(channels, _MapSet(view, norb))

    def ed_D(self, v1):
        channels, smap = self.channels_and_map()
        fx = self.fx

        @functools.lru_cache(maxsize=None)
        def chi_hf_sub(v):
            terms = _terms(self.t, fx, self.a, self.b, self.R, v)
            full = ed_oracle_util.SectorED(
                fx, terms=terms).bond_correlator(channels)
            hf = ed_oracle_util.SectorED(
                fx, terms=(), h1=ed_oracle_util.hf_h1_from_terms(fx, terms)
            ).bond_correlator(channels)
            return full - hf

        def mapped(v, s1, s2):
            return chi_hf_sub(v)[:, s1, s2][:, smap][:, :, smap]

        def X_S(v):
            return mapped(v, 0, 0) - mapped(v, 0, 1)

        def X_C(v):
            return mapped(v, 0, 0) + mapped(v, 0, 1)

        return (ed_oracle_util.richardson(X_S, v1),
                ed_oracle_util.richardson(X_C, v1))


class _MapSet:
    """What ``ed_to_solver_bond_map`` reads from a bond set: the channel
    list and ``norb`` (through ``v_bond[0].shape``)."""

    def __init__(self, view, norb):
        self.delta_r = view.delta_r
        self.n_channels = view.n_channels
        self.v_bond = [np.zeros((norb, norb), complex)] * view.n_channels


@functools.lru_cache(maxsize=None)
def _chibar(fx, shells, nmat):
    """The bare bond bubble on the fixed topology (independent of the
    direction and of its amplitude)."""
    topo = bc.resolve_bond_topology(
        _decls("CoulombInter", shells[0], 0, 0, 0.0, shells), np.eye(3),
        fx.norb, active_types=_ALL)
    view = bc.BondSetView(topo)
    green = free_green(fx, nmat)
    return bc.bond_bubble(green, view, beta=fx.beta)[:, 0, 0]


@functools.lru_cache(maxsize=None)
def _ed_side(fx_name, t, a, b, R, shells, v1):
    """The ED-side Richardson derivative, computed ONCE per direction and
    step (the G3 no-bond comparison reuses the same ED data)."""
    fx = {"fx5": fx5(), "fx3": fx3()}[fx_name]
    return _Side(fx, t, a, b, R, shells).ed_D(v1)


@functools.lru_cache(maxsize=None)
def adjudicate(fx_name, t, a, b, R, shells, bond_blocks=True):
    """The canonical S and C records for one direction, plus the raw
    finer-estimate arrays (for the G3 claim)."""
    fx = {"fx5": fx5(), "fx3": fx3()}[fx_name]
    side = _Side(fx, t, a, b, R, shells)
    label = "lb/{}/{}/R{}_{}{}{}".format(
        fx_name, t, R, a, b, "" if bond_blocks else "/nobond")
    dS, dC = side.dW(bond_blocks)
    Dn_S, Dn_C = pred_first_order(side.chibar(NMAT), dS, dC)
    D2_S, D2_C = pred_first_order(side.chibar(2 * NMAT), dS, dC)
    E1_S, E1_C = _ed_side(fx_name, t, a, b, R, shells, CAMPAIGN_V1)
    Eh_S, Eh_C = _ed_side(fx_name, t, a, b, R, shells, CAMPAIGN_V1 / 2)
    free = np.ones(dS.shape[1:], bool)
    zS = ed_oracle_util.projected_structural_zero_mask(
        free, np.any(np.abs(dS) > 0, axis=0))
    zC = ed_oracle_util.projected_structural_zero_mask(
        free, np.any(np.abs(dC) > 0, axis=0))
    rec_S = ed_oracle_util.adjudicate_granule(
        E1_S, Eh_S, Dn_S, D2_S, zS, label + "/S")
    rec_C = ed_oracle_util.adjudicate_granule(
        E1_C, Eh_C, Dn_C, D2_C, zC, label + "/C")
    return rec_S, rec_C, (Eh_S, Eh_C, D2_S, D2_C)


class TestFramePins(unittest.TestCase):
    """The bubble used here (``bond_channels.bond_bubble``, the #151
    oracle's) and the production bubble (``bubble.bond_bubble_static``)
    are the SAME frame -- so the ED map derived for the former applies
    to the vertex built for the latter."""

    def test_static_bubble_equals_the_oracle_bubble_static_slice(self):
        fx = fx5()
        shells = (1, 2)
        nmat = 64
        cb = _chibar(fx, shells, nmat)                      # (L, ND, ND)
        topo = bc.resolve_bond_topology(
            _decls("CoulombInter", 1, 0, 0, 0.0, shells), np.eye(3),
            fx.norb, active_types=_ALL)
        g = free_green(fx, nmat)[:, :, :, 0, 0, :]          # (norb, norb, L, nmat)
        g = np.moveaxis(g, (0, 1, 2, 3), (2, 3, 1, 0))[None]  # (1, nmat, L, norb, norb)
        cs = bubble.bond_bubble_static(
            np.ascontiguousarray(g), None, fx.beta, bc.BondSetView(topo),
            spatial_shape=(fx.L, 1, 1))
        np.testing.assert_allclose(np.asarray(cs), cb, rtol=0, atol=1e-12)

    def test_dw_support_is_the_documented_pattern(self):
        side = _Side(fx5(), "Hund", 0, 0, 1, (1, 2))
        dS, dC = side.dW()
        supp_S = np.any(np.abs(dS) > 0, axis=0)
        supp_C = np.any(np.abs(dC) > 0, axis=0)
        view = side.view()
        mp = view.delta_r.index((1, 0, 0))
        mm = view.delta_r.index((-1, 0, 0))
        exp_S = np.zeros((5, 5), bool)
        exp_S[mp, mp] = exp_S[mm, mm] = True   # the bond-diagonal Fock slots
        exp_C = exp_S.copy()
        exp_S[0, 0] = True                     # Hund's density content (+1, -1):
        exp_C[0, 0] = True                     # the Hartree slot in BOTH channels
        self.assertTrue(np.array_equal(supp_S, exp_S), (supp_S, exp_S))
        self.assertTrue(np.array_equal(supp_C, exp_C), (supp_C, exp_C))


class TestLongitudinalBondGranulesFx5(unittest.TestCase):
    """G2 on the one-orbital ring (a == b bonds): CoulombInter g1/g2,
    Hund g1, Ising g1 -- PASS asserted; Exchange g1 recorded."""

    def _assert_pass(self, rec_S, rec_C):
        self.assertEqual(rec_S["status"], "PASS", rec_S)
        self.assertEqual(rec_C["status"], "PASS", rec_C)

    @heavy
    def test_coulombinter_g1_g2(self):
        for R in (1, 2):
            with self.subTest(R=R):
                rec_S, rec_C, _ = adjudicate("fx5", "CoulombInter", 0, 0, R, (1, 2))
                self._assert_pass(rec_S, rec_C)

    @heavy
    def test_hund_g1(self):
        rec_S, rec_C, _ = adjudicate("fx5", "Hund", 0, 0, 1, (1, 2))
        self._assert_pass(rec_S, rec_C)

    @heavy
    def test_ising_g1(self):
        rec_S, rec_C, _ = adjudicate("fx5", "Ising", 0, 0, 1, (1, 2))
        self._assert_pass(rec_S, rec_C)

    @heavy
    def test_g3_bond_blocks_are_load_bearing(self):
        """G3, the falsifiable claim: with the bond blocks the first-order
        response is reproduced at the granule tolerance; with them zeroed
        (Tier 1's Hartree-only reading embedded in the bond frame) the
        SAME ED data is missed by at least 10 x tol."""
        for t in ("CoulombInter", "Hund", "Ising"):
            with self.subTest(t=t):
                rec_S, rec_C, (E_S, E_C, _, _) = adjudicate("fx5", t, 0, 0, 1, (1, 2))
                self._assert_pass(rec_S, rec_C)
                _, _, (_, _, P_S, P_C) = adjudicate("fx5", t, 0, 0, 1, (1, 2),
                                                    bond_blocks=False)
                miss = max(np.max(np.abs(E_S - P_S)), np.max(np.abs(E_C - P_C)))
                tol = max(rec_S["tol"], rec_C["tol"])
                self.assertGreaterEqual(miss, 10.0 * tol, (t, miss, tol))

    @heavy
    def test_exchange_g1_oracle_only_recorded(self):
        rec_S, rec_C, _ = adjudicate("fx5", "Exchange", 0, 0, 1, (1, 2))
        print("ORACLE-ONLY Exchange fx5 R=1 (0,0):",
              rec_S["status"], rec_C["status"],
              "tol", rec_S["tol"], rec_C["tol"],
              "max_signal", rec_S["max_signal"], rec_C["max_signal"])
        self.assertIn(rec_S["status"], ("PASS", "PASS-ZERO", "FAIL", "INCONCLUSIVE"))


class TestLongitudinalBondGranulesCaseM(unittest.TestCase):
    """G2 on the two-orbital case-M ring (a != b bond, R = +1, single
    shell): CoulombInter, Hund, Ising -- PASS asserted; Exchange recorded."""

    def _assert_pass(self, rec_S, rec_C):
        self.assertEqual(rec_S["status"], "PASS", rec_S)
        self.assertEqual(rec_C["status"], "PASS", rec_C)

    @heavy
    def test_coulombinter_interorbital(self):
        rec_S, rec_C, _ = adjudicate("fx3", "CoulombInter", 0, 1, 1, (1,))
        self._assert_pass(rec_S, rec_C)

    @heavy
    def test_hund_interorbital(self):
        rec_S, rec_C, _ = adjudicate("fx3", "Hund", 0, 1, 1, (1,))
        self._assert_pass(rec_S, rec_C)

    @heavy
    def test_ising_interorbital(self):
        rec_S, rec_C, _ = adjudicate("fx3", "Ising", 0, 1, 1, (1,))
        self._assert_pass(rec_S, rec_C)

    @heavy
    def test_g3_bond_blocks_are_load_bearing_interorbital(self):
        for t in ("CoulombInter", "Hund", "Ising"):
            with self.subTest(t=t):
                rec_S, rec_C, (E_S, E_C, _, _) = adjudicate("fx3", t, 0, 1, 1, (1,))
                self._assert_pass(rec_S, rec_C)
                _, _, (_, _, P_S, P_C) = adjudicate("fx3", t, 0, 1, 1, (1,),
                                                    bond_blocks=False)
                miss = max(np.max(np.abs(E_S - P_S)), np.max(np.abs(E_C - P_C)))
                tol = max(rec_S["tol"], rec_C["tol"])
                self.assertGreaterEqual(miss, 10.0 * tol, (t, miss, tol))

    @heavy
    def test_exchange_interorbital_oracle_only_recorded(self):
        rec_S, rec_C, _ = adjudicate("fx3", "Exchange", 0, 1, 1, (1,))
        print("ORACLE-ONLY Exchange fx3 R=1 (0,1):",
              rec_S["status"], rec_C["status"],
              "tol", rec_S["tol"], rec_C["tol"],
              "max_signal", rec_S["max_signal"], rec_C["max_signal"])
        self.assertIn(rec_S["status"], ("PASS", "PASS-ZERO", "FAIL", "INCONCLUSIVE"))


if __name__ == "__main__":
    unittest.main()
