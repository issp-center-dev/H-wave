#!/usr/bin/env python3
"""Off-site Exchange (and PairLift) carry NO q-representable longitudinal
spin/charge content: the #181 Tier-2 exact-diagonalization adjudication.

Question. The general (full-vertex) FLEX path solves the LONGITUDINAL
(spin-conserving) spin/charge channels with a q-only pair-space vertex.
Issue #181's investigation suggested that an off-site Exchange term,
``J c+_{ia,up} c_{jb,up} c+_{jb,dn} c_{ia,dn}`` (+ h.c.), might need a new
density-family ``(S, C)`` entry in :mod:`hwave.solver.vertex_table`, because
it regroups into ``-J S+_{ia} S-_{jb}``, a product of two LOCAL bilinears.
Those bilinears are spin-FLIP operators, though, so that regrouping lives
in the transverse channel; the longitudinal question has to be settled by
exact diagonalization.

Method (the #113 / bond-campaign recipe, tests/ed_oracle_util.py): on a
periodic ring, ``chi0 = to_solver_slots(chi_connected(fx))``; the
first-order vertex part of the response is ``vertex(v) =
chi_full(v) - chi_HF-only(v)`` (the O(v) self-energy removed exactly by the
Hartree-Fock insertion), differentiated at v -> 0 by Richardson. The
generalized pair space ``(a, c) x (b, d)`` splits into the spin-conserving
(longitudinal) and spin-flip sectors; in the longitudinal sector the
effective q-only vertex ``Gamma(q) = -chi0(q)^-1 D(q) chi0(q)^-1`` (the ring's
``chi = [1 + chi0 W]^-1 chi0`` convention) is read off per slot family and
per channel, ``S = W_cross - W_same``, ``C = W_cross + W_same``.

Findings (2026-09-03, recorded in the assertions below):

* CONTROLS reproduce the adjudicated table: an off-site density bond gives
  CoulombInter ``(S, C) = (0, +2) V(q)`` and Hund ``(+1, -1) V(q)`` on the
  density slots, with ``V(q) = 2 v cos q`` for the two-sided one-orbital
  bond and ``v e^{+iq}`` for the one-sided two-orbital bond. Measured
  discrepancy: 0.2% on the two-orbital density slots (1.9964 vs 2), 1-3%
  on the one-orbital ring (3.949 vs 4; 1.949 vs 2), where the bond's own
  exchange (Fock) crossing -- a NON-local pair that a q-only vertex cannot
  carry (the part Tier 1 omits with a warning; #181 Tier 3) -- lands on
  the same single slot. The ENFORCED tolerances below (0.05-0.06 absolute
  per unit coupling) are looser than the measurement: they bound the
  Fock residual, not the round-off.
* OFF-SITE EXCHANGE: the longitudinal first-order response is 40x smaller
  than the spin-flip one (1-orbital: 3.4e-3 vs 2.2e-1; 2-orbital: 4.8e-3
  vs 1.9e-1 per unit coupling) and carries NO density-family content
  (2-orbital density slots 1e-3 vs the controls' 2.0 / 1.0); what remains
  coincides, number for number, with the controls' Fock-crossing residual
  (1-orbital: S = C = 0.0513 / 0.0297 / 0.0083 over q against the
  CoulombInter control's S-channel residual 0.0509 / 0.0293 / 0.0076) --
  i.e. it is the same non-local-pair content, not a q-only vertex. For a
  purely imaginary J the longitudinal response is numerically zero
  (1.5e-7 on the one-orbital fixture, bound 1e-6; no analytic identity is
  claimed) while the spin-flip response stays at 1.6e-1.
* OFF-SITE PAIRLIFT (Sz-breaking; dense route): longitudinal response
  8e-9, spin-flip 2.2e-1 -- inert in the longitudinal channels off-site
  exactly as on-site.

Consequence: :data:`vertex_table.ADJUDICATED_SC` needs no off-site entry for
Exchange (its ``cross`` content is on-site only) or PairLift; the general
FLEX path's rejection of off-site Exchange stands, now on adjudicated
grounds (the physics of the term is transverse, outside the spin-free
longitudinal solve); off-site PairLift stays accepted-and-inert.
"""

import unittest

import numpy as np

from tests import ed_oracle_util as U
from tests.heavy_tests import heavy


V1 = 0.005


def _fixture_1orb():
    # the bond campaign's fx5: odd ring (non-self-inverse q), complex hopping
    return U.EDFixture(L=5, norb=1, t={(0, 0): 0.7 * np.exp(0.3j)},
                       eps=(0.0,), T=0.5, mu=0.2)


def _fixture_2orb():
    # the bond campaign's fx3 / case M
    return U.EDFixture(L=3, norb=2,
                       t={(0, 0): 0.5 + 0.2j, (1, 1): 0.35 - 0.15j,
                          (0, 1): 0.1 + 0.05j, (1, 0): 0.1 + 0.05j},
                       eps=(0.05, -0.03), T=0.45, mu=0.1)


def _fixture_pairlift():
    # dense route (PairLift breaks Sz): small ring
    return U.EDFixture(L=3, norb=1, t={(0, 0): 0.7 * np.exp(0.3j)},
                       eps=(0.0,), T=0.5, mu=0.2)


def _chi(fx, terms=None, h1=None, dense=False):
    if dense:
        hint = None if terms is None else U.h_int_from_terms(fx, terms)
        return U.to_solver_slots(U.chi_connected(fx, hint=hint, h1=h1))
    if terms is None:
        return U.to_solver_slots(U.SectorED(fx, terms=(), h1=h1).chi_connected())
    return U.to_solver_slots(U.SectorED(fx, terms).chi_connected())


def _first_order(fx, terms_of_v, v1=V1, dense=False):
    """Richardson derivative at v -> 0 of the HF-subtracted vertex part."""
    shape = (fx.L,) + (2 * fx.norb,) * 4

    def vertex(v):
        if v == 0.0:
            return np.zeros(shape, dtype=complex)
        t = terms_of_v(v)
        return (_chi(fx, terms=t, dense=dense)
                - _chi(fx, h1=U.hf_h1_from_terms(fx, t), dense=dense))

    return U.richardson(vertex, v1)


def _sectors(fx, D):
    """max|D| over the longitudinal (spin-conserving) and spin-flip pair
    sectors, per q."""
    nd = 2 * fx.norb
    npair = nd * nd
    Dm = D.reshape(fx.L, npair, npair)
    long_mask = np.array([(a // fx.norb) == (c // fx.norb)
                          for a in range(nd) for c in range(nd)])
    d_long = np.array([np.max(np.abs(Dm[q][np.ix_(long_mask, long_mask)]))
                       for q in range(fx.L)])
    d_flip = np.array([np.max(np.abs(Dm[q][np.ix_(~long_mask, ~long_mask)]))
                       for q in range(fx.L)])
    return d_long, d_flip


def _longitudinal_sc(fx, D):
    """Effective longitudinal q-only vertex in S/C form, per q: dicts
    ``S[q][(l1 l2, l3 l4)]``, ``C[q][...]`` over orbital-pair slots."""
    nd = 2 * fx.norb
    norb = fx.norb
    npair = nd * nd
    chi0 = _chi(fx).reshape(fx.L, npair, npair)
    Dm = D.reshape(fx.L, npair, npair)
    long_idx = [a * nd + c for a in range(nd) for c in range(nd)
                if (a // norb) == (c // norb)]
    pairs = [((a // norb), a % norb, c % norb) for a in range(nd) for c in range(nd)
             if (a // norb) == (c // norb)]
    S_all, C_all = [], []
    for q in range(fx.L):
        c0 = chi0[q][np.ix_(long_idx, long_idx)]
        dq = Dm[q][np.ix_(long_idx, long_idx)]
        G = -np.linalg.solve(c0, np.linalg.solve(c0.conj().T, dq.conj().T).conj().T)
        n = norb * norb
        W_same = np.zeros((n, n), dtype=complex)
        W_cross = np.zeros((n, n), dtype=complex)
        for i, (s1, o1, o2) in enumerate(pairs):
            for j, (s2, o3, o4) in enumerate(pairs):
                if s1 == 0 and s2 == 0:
                    W_same[o1 * norb + o2, o3 * norb + o4] = G[i, j]
                elif s1 == 0 and s2 == 1:
                    W_cross[o1 * norb + o2, o3 * norb + o4] = G[i, j]
        S_all.append(W_cross - W_same)
        C_all.append(W_cross + W_same)
    return np.array(S_all), np.array(C_all)


def _terms_exchange(fx, a, b, R, J):
    from tests.test_bond_transverse_ed import _terms_exchange_offsite
    return _terms_exchange_offsite(fx, a, b, R, J)


def _terms_hund(fx, a, b, R, v):
    from tests.test_bond_transverse_ed import _terms_hund_offsite
    return _terms_hund_offsite(fx, a, b, R, v)


def _terms_pairlift(fx, a, b, R, v):
    from tests.test_bond_transverse_ed import _terms_pairlift_offsite
    return _terms_pairlift_offsite(fx, a, b, R, v)


class TestOffsiteExchangeLongitudinalOneOrbital(unittest.TestCase):
    """Fast (< 3 s): the one-orbital ring, SectorED route."""

    @classmethod
    def setUpClass(cls):
        cls.fx = _fixture_1orb()

    def test_real_exchange_has_no_q_only_longitudinal_content(self):
        fx = self.fx
        D = _first_order(fx, lambda v: _terms_exchange(fx, 0, 0, 1, v))
        d_long, d_flip = _sectors(fx, D)
        # the term acts in the spin-flip sector (measured 3.4e-3 vs 2.2e-1)
        self.assertGreater(d_flip.max(), 0.1)
        self.assertLess(d_long.max(), 0.05 * d_flip.max())
        # what little longitudinal response exists is not a q-only density
        # entry: an entry x on the (00,00) slot would read x * 2 cos q, i.e.
        # 2x / 0.62x / -1.6x over this ring; the measured S = C = 0.051 /
        # 0.030 / 0.008 has the wrong shape (all positive) and is < 3% of
        # the controls' unit content (2.0 / 1.0). Bound: 0.06 per unit J.
        S, C = _longitudinal_sc(fx, D)
        self.assertLess(np.max(np.abs(S)), 0.06)
        self.assertLess(np.max(np.abs(C)), 0.06)
        cosq = 2.0 * np.cos(2.0 * np.pi * np.arange(fx.L) / fx.L)
        s_q = S[:, 0, 0].real
        x = np.dot(s_q, cosq) / np.dot(cosq, cosq)     # best q-only fit
        self.assertGreater(np.max(np.abs(s_q - x * cosq)), 0.5 * np.max(np.abs(s_q)))

    def test_imaginary_exchange_longitudinal_response_vanishes(self):
        fx = self.fx
        D = _first_order(fx, lambda v: _terms_exchange(fx, 0, 0, 1, 1j * v))
        d_long, d_flip = _sectors(fx, D)
        self.assertGreater(d_flip.max(), 0.1)          # measured 1.6e-1
        self.assertLess(d_long.max(), 1e-6)            # measured 1.5e-7

    def test_offsite_pairlift_is_longitudinally_inert(self):
        fx = _fixture_pairlift()
        D = _first_order(fx, lambda v: _terms_pairlift(fx, 0, 0, 1, v),
                         dense=True)
        d_long, d_flip = _sectors(fx, D)
        self.assertGreater(d_flip.max(), 0.05)
        self.assertLess(d_long.max(), 1e-6)


class TestOffsiteExchangeLongitudinalControls(unittest.TestCase):
    """Heavy: the controls that calibrate the bounds above (one-orbital)
    and the two-orbital inter-orbital case with its controls."""

    @heavy
    def test_one_orbital_controls_and_the_fock_residual_identity(self):
        fx = _fixture_1orb()
        cosq = 2.0 * np.cos(2.0 * np.pi * np.arange(fx.L) / fx.L)
        # CoulombInter: two-sided bond via canonical_density_terms (one
        # unordered class), density C = 2 V(q) = 2 * (2 v cos q) / v
        D_ci = _first_order(fx, lambda v: U.canonical_density_terms(fx, [(0, 0, 1, v)]))
        S_ci, C_ci = _longitudinal_sc(fx, D_ci)
        np.testing.assert_allclose(C_ci[:, 0, 0].real, 2.0 * cosq, atol=0.06)
        self.assertLess(np.max(np.abs(S_ci)), 0.06)
        # Hund: (S, C) = (+1, -1) V(q)
        D_h = _first_order(fx, lambda v: _terms_hund(fx, 0, 0, 1, v))
        S_h, C_h = _longitudinal_sc(fx, D_h)
        np.testing.assert_allclose(S_h[:, 0, 0].real, cosq, atol=0.06)
        np.testing.assert_allclose(C_h[:, 0, 0].real, -cosq, atol=0.06)
        # the Exchange term's whole longitudinal content equals the density
        # bond's Fock-crossing residual (the S channel of CoulombInter,
        # which has no q-only S content): the same non-local pair
        D_x = _first_order(fx, lambda v: _terms_exchange(fx, 0, 0, 1, v))
        S_x, C_x = _longitudinal_sc(fx, D_x)
        np.testing.assert_allclose(S_x[:, 0, 0], S_ci[:, 0, 0], atol=3e-3)
        np.testing.assert_allclose(C_x[:, 0, 0], S_ci[:, 0, 0], atol=3e-3)
        self.assertGreater(np.max(np.abs(S_ci[:, 0, 0])), 0.03)   # not vacuous
        # Richardson stability: the step 0.005 is not what makes the
        # Exchange content small (measured |D(v1) - D(v1/2)| = 2.7e-6)
        D_x_half = _first_order(fx, lambda v: _terms_exchange(fx, 0, 0, 1, v),
                                v1=V1 / 2)
        self.assertLess(np.max(np.abs(D_x_half - D_x)), 1e-4)

    @heavy
    def test_two_orbital_interorbital_exchange_and_controls(self):
        fx = _fixture_2orb()
        eq = np.exp(1j * 2.0 * np.pi * np.arange(fx.L) / fx.L)
        # slots at norb=2: 0=(00), 1=(01), 2=(10), 3=(11); density (00,11)=(0,3).
        # The one-sided bond (orbital 0 at j, orbital 1 at j+1) is V_01(q)
        # = v e^{+iq}; the (00,11) slot reads the pair-transposed entry
        # (hwave.sc's _build_interaction_k convention, issue #96), i.e.
        # e^{-iq}, and (11,00) reads e^{+iq}.
        D_ci = _first_order(fx, lambda v: U.canonical_density_terms(fx, [(0, 1, 1, v)]))
        S_ci, C_ci = _longitudinal_sc(fx, D_ci)
        np.testing.assert_allclose(C_ci[:, 0, 3], 2.0 * eq.conj(), atol=0.05)
        np.testing.assert_allclose(C_ci[:, 3, 0], 2.0 * eq, atol=0.05)
        self.assertLess(np.max(np.abs(S_ci[:, 0, 3])), 0.02)
        # the same placement as the solver's k-space builder: the (00,11)
        # slot reads mat[0, 1] of hwave.sc._build_interaction_k, i.e. the
        # pair-transposed entry (issue #96)
        import hwave.sc as sc
        kx = 2.0 * np.pi * np.arange(fx.L) / fx.L
        ik = sc._build_interaction_k(kx, np.array([0.0]), np.array([0.0]),
                                     {"CoulombInter": {((1, 0, 0), (0, 1)): 1.0,
                                                       ((-1, 0, 0), (1, 0)): 1.0}},
                                     2)["CoulombInter"]
        np.testing.assert_allclose(C_ci[:, 0, 3], 2.0 * ik[0, 1, :, 0, 0], atol=0.05)
        np.testing.assert_allclose(C_ci[:, 3, 0], 2.0 * ik[1, 0, :, 0, 0], atol=0.05)
        D_h = _first_order(fx, lambda v: _terms_hund(fx, 0, 1, 1, v))
        S_h, C_h = _longitudinal_sc(fx, D_h)
        np.testing.assert_allclose(S_h[:, 0, 3], eq.conj(), atol=0.05)
        np.testing.assert_allclose(C_h[:, 0, 3], -eq.conj(), atol=0.05)
        for label, phase in (("real", 1.0), ("imag", 1j)):
            with self.subTest(J=label):
                D_x = _first_order(fx, lambda v, p=phase: _terms_exchange(fx, 0, 1, 1, p * v))
                d_long, d_flip = _sectors(fx, D_x)
                self.assertGreater(d_flip.max(), 0.1)      # measured 1.9e-1
                self.assertLess(d_long.max(), 0.05 * d_flip.max())   # 4.8e-3
                S_x, C_x = _longitudinal_sc(fx, D_x)
                # density family: nothing (measured <= 3.6e-3 vs 2.0 / 1.0)
                for (r, c) in ((0, 3), (3, 0)):
                    self.assertLess(abs(S_x[:, r, c]).max(), 0.01)
                    self.assertLess(abs(C_x[:, r, c]).max(), 0.01)
                # everything else is at the controls' non-local residual
                # level (cross slots 0.011 / 0.031, same numbers as the
                # CoulombInter control's cross residual)
                self.assertLess(np.max(np.abs(S_x)), 0.05)
                self.assertLess(np.max(np.abs(C_x)), 0.05)
                if label == "real":
                    # ... and for real J it IS that residual, element for
                    # element on the cross slots (01,01), (10,10): the
                    # Exchange term's (S, C) both equal the CoulombInter
                    # control's S-channel Fock residual there (measured
                    # 0.0113 at q=0, 0.0307/0.0038 at q != 0) -- a nonzero
                    # identity, so an all-zero implementation cannot pass
                    for (r, c) in ((1, 1), (2, 2)):
                        np.testing.assert_allclose(S_x[:, r, c], S_ci[:, r, c], atol=3e-3)
                        np.testing.assert_allclose(C_x[:, r, c], S_ci[:, r, c], atol=3e-3)
                    self.assertGreater(np.max(np.abs(S_ci[:, 1, 1])), 0.02)


if __name__ == "__main__":
    unittest.main()
