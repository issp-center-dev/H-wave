#!/usr/bin/env python3

"""``TotalNED`` contract and calibration tests (Phase S, Task 1 of the
bond-transverse campaign, issue #110's validation groundwork).

See ``docs/superpowers/specs/2026-08-15-bond-transverse-design.md``
("Phase S -- the full spinful transverse") and
``docs/superpowers/plans/2026-08-17-bond-transverse-phase-s.md`` for the
binding spec/plan this module works under. ``TotalNED``
(``tests/ed_oracle_util.py``) is the total-N-sector ED oracle Task 1
built to represent Sz-BREAKING one-body Hamiltonians (e.g. a spin-orbit
hop ``t_so * c^dagger_{...,up} c_{...,dn}``), which ``SectorED``'s
(N_up, N_dn) blocking cannot host at all.

This module covers:

- construction: on a Sz-breaking ``h1``, ``SectorED`` must reject it
  (``TestTotalNEDConstructionContract``);
- the ``bond_correlator_transverse`` operator-construction contract
  (shape/dtype, R-wrapping duplicate rejection, the mixed-spin sector
  delta) -- mirroring Phase H's
  ``TestBondCorrelatorTransverseContract`` list, but for the total-N
  delta (always 0, never SectorED's ``(+1, -1)``) --
  (``TestBondCorrelatorTransverseContractTotalN``);
- the calibration gate: on a Sz-CONSERVING fixture (with interactions on
  and a R != 0 channel), ``TotalNED.bond_correlator_transverse`` must
  equal ``SectorED.bond_correlator_transverse`` at round-off -- both are
  exact ED on the identical Hilbert space, so any discrepancy is a real
  bug in one of the two engines (``TestTotalNEDCalibration``);
- the Sz-breaking discrimination check: on a genuinely Sz-breaking
  fixture, the mixed-spin (transverse) correlator blocks TotalNED
  computes are independently, measurably NONZERO -- confirming the
  fixture class Gate S0/S2 will build on actually exercises spin mixing,
  not an accidental near-zero response (``TestTotalNEDSzBreakingDiscrimination``).

Tests must be run from the repository root.
"""

import unittest

import numpy as np

from tests import ed_oracle_util
from tests.approx_util import assert_approx_array


def _h1_with_spin_orbit(fx, t_so):
    """``fx.build_h1()`` plus an on-site spin-flip hop
    ``t_so * c^dagger_{j,o,up} c_{j,o,dn} + h.c.``, per site/orbital -- a
    Sz-BREAKING one-body term (spin-orbit-like) that ``SectorED`` cannot
    host (its (N_up, N_dn) sector blocking assumes ``h1`` conserves each
    separately) but ``TotalNED`` can (its blocking is by TOTAL N only;
    see ``TotalNED``'s class docstring for the sector-algebra
    derivation: a bilinear of shape ``c^dagger_p c_q`` always conserves
    total N, whatever the spin parities of ``p`` and ``q`` are)."""
    h1 = fx.build_h1().copy()
    for j in range(fx.L):
        for o in range(fx.norb):
            up = fx.mode(j, o, 0)
            dn = fx.mode(j, o, 1)
            h1[up, dn] += t_so
            h1[dn, up] += np.conj(t_so)
    return h1


class TestTotalNEDConstructionContract(unittest.TestCase):
    """Construction on a small spinful fixture with a Sz-BREAKING ``h1``:
    ``SectorED`` must reject it (``_sector_hamiltonian``'s dropped-element
    floor check fires, since a spin-off-diagonal ``h1`` entry moves
    probability between (N_up, N_dn) sectors); ``TotalNED`` must accept it
    (total N is still conserved for the very same bilinear -- see
    ``TotalNED``'s class docstring)."""

    @staticmethod
    def _fx():
        return ed_oracle_util.EDFixture(
            L=2, norb=1, t={(0, 0): 0.4 - 0.1j}, eps=(0.05,), T=0.5, mu=0.1)

    def test_sectored_rejects_sz_breaking_h1(self):
        fx = self._fx()
        h1 = _h1_with_spin_orbit(fx, 0.3 + 0.05j)
        with self.assertRaises(AssertionError):
            ed_oracle_util.SectorED(fx, h1=h1)

    def test_totalned_accepts_sz_breaking_h1(self):
        fx = self._fx()
        h1 = _h1_with_spin_orbit(fx, 0.3 + 0.05j)
        te = ed_oracle_util.TotalNED(fx, h1=h1)   # must not raise
        # Sanity: the sector bookkeeping actually built something usable
        # (every total-N key from 0 to fx.nmode present, weights sum to 1).
        self.assertEqual(sorted(te._sector_states.keys()),
                          list(range(fx.nmode + 1)))
        total_w = sum(w.sum() for w in te._w.values())
        assert_approx_array(total_w, 1.0, rel=0, abs=1e-12)


class TestBondCorrelatorTransverseContractTotalN(unittest.TestCase):
    """Operator-construction contract for
    ``TotalNED.bond_correlator_transverse``, mirroring the Phase-H list
    (``test_bond_transverse_ed.py::TestBondCorrelatorTransverseContract``)
    against ``SectorED``'s IDENTICAL method: shape/dtype, R-wrapping
    duplicate rejection, and the mixed-spin sector displacement -- which
    is ``0`` here (not SectorED's ``(+1, -1)``; see ``TotalNED``'s class
    docstring's sector algebra)."""

    @staticmethod
    def _small_fx():
        # Same fixture shape as Phase H's contract-test fixture (2-site,
        # two-orbital, spinful) -- fast, disposable, norb=2 so channels
        # use the (R, a, b) 3-tuple form ``_channel_orbitals`` expects.
        t = {(0, 0): 0.3 + 0.1j, (1, 1): 0.2 - 0.05j,
             (0, 1): 0.1j, (1, 0): 0.1j}
        return ed_oracle_util.EDFixture(
            L=2, norb=2, t=t, eps=(0.05, -0.02), T=0.5, mu=0.1)

    def test_output_shape_and_dtype(self):
        fx = self._small_fx()
        te = ed_oracle_util.TotalNED(fx)
        channels = [(0, 0, 0), (1, 0, 0)]
        out = te.bond_correlator_transverse(channels)
        self.assertEqual(out.shape, (fx.L, len(channels), len(channels)))
        self.assertEqual(out.dtype, np.complex128)

    def test_duplicate_channels_after_wrapping_rejected(self):
        fx = self._small_fx()
        te = ed_oracle_util.TotalNED(fx)
        # (1, 0, 0) and (1 - L, 0, 0) both wrap (mod L=2) to R=1 -- same
        # channel after wrapping, must be rejected rather than silently
        # double-counted.
        channels = [(1, 0, 0), (1 - fx.L, 0, 0)]
        with self.assertRaises(ValueError):
            te.bond_correlator_transverse(channels)

    def test_mixed_spin_delta_is_zero(self):
        # Probe the sector displacement via the same operator-
        # construction API bond_correlator_transverse uses internally: a
        # bilinear whose creation leg is spin UP (mode parity 0) and
        # annihilation leg is spin DOWN (mode parity 1) must impart
        # delta = 0 under TOTAL-N blocking -- unlike SectorED, where the
        # SAME bilinear gives delta = (+1, -1) (see the class docstring's
        # sector-algebra derivation: total N = N_up + N_dn is invariant
        # under any single creation-plus-annihilation bilinear).
        fx = self._small_fx()
        te = ed_oracle_util.TotalNED(fx)
        mode_terms = [(fx.mode(j, 0, 0), fx.mode(j, 0, 1), 1.0)
                      for j in range(fx.L)]
        delta, _op = te._build_operator(mode_terms)
        self.assertEqual(delta, 0)


class TestTotalNEDCalibration(unittest.TestCase):
    """The calibration gate (Task 1, Step 4): on a Sz-CONSERVING fixture
    WITH interactions on (a real off-site CoulombInter term, not just the
    free chain) and including a R != 0 channel, ``TotalNED.
    bond_correlator_transverse`` must equal ``SectorED.
    bond_correlator_transverse`` at ``assert_approx_array(rel=0,
    abs=1e-12)`` -- both engines diagonalize the IDENTICAL Hamiltonian on
    the IDENTICAL full Hilbert space (TotalNED's total-N blocks are a
    coarser partition of the SAME Fock space SectorED's (N_up, N_dn)
    blocks refine), so any discrepancy beyond round-off is a real bug
    in one of the two sector-block engines, not a physics difference.

    Fixture: L=3, norb=1 ring (dim=64; TotalNED's largest total-N block
    is N=3, dim=C(6,3)=20 -- small and fast), the SAME (t, eps, T, mu) as
    ``test_bond_transverse_ed.py``'s
    ``test_granule_pairlift_offsite_passzero_control_dense`` fixture (a
    known-good parameter set already exercised in this campaign).
    Interaction: real off-site CoulombInter V=0.4 at R=1 (via
    ``canonical_density_terms`` -- V kept real, since a complex
    coefficient on a Hermitian density-density product would break
    Hermiticity of H). Channels: R=0, 1, 2 (norb=1's ``(R,)`` shorthand),
    so R != 0 is exercised on both legs.
    """

    def test_matches_sectored_with_interactions_and_offsite_channel(self):
        fx = ed_oracle_util.EDFixture(
            L=3, norb=1, t={(0, 0): 0.55 * np.exp(-0.15j)}, eps=(0.0,),
            T=0.55, mu=0.1)
        V = 0.4
        terms = ed_oracle_util.canonical_density_terms(fx, [(0, 0, 1, V)])
        channels = [(0,), (1,), (2,)]

        se = ed_oracle_util.SectorED(fx, terms=terms)
        te = ed_oracle_util.TotalNED(fx, terms=terms)

        want = se.bond_correlator_transverse(channels)
        got = te.bond_correlator_transverse(channels)

        self.assertEqual(got.shape, want.shape)
        assert_approx_array(got, want, rel=0, abs=1e-12)


class TestTotalNEDSzBreakingDiscrimination(unittest.TestCase):
    """The Sz-breaking discrimination check (Task 1, Step 5): on a
    genuinely Sz-breaking fixture (on-site spin-orbit hop ``t_so``, see
    ``_h1_with_spin_orbit``), the mixed-spin transverse correlator blocks
    ``TotalNED.bond_correlator_transverse`` computes must be
    independently, measurably NONZERO -- confirming this fixture class
    genuinely exercises spin mixing (rather than landing, by accident, on
    a near-zero response that would make a later Gate S0/S2 pin vacuous).

    Fixture: L=2, norb=1 (dim=16, tiny/fast), t=0.4-0.1j, eps=0.05,
    T=0.5, mu=0.1, on-site spin-orbit t_so=0.3+0.05j (same fixture and
    t_so as ``TestTotalNEDConstructionContract``, so the two tests probe
    the identical Hamiltonian from two angles: constructibility, then
    physical content).

    Measured lower bound (recorded per the brief's "assert a lower bound
    measured and recorded in the docstring"): the R=0 diagonal entry
    ``out[0, 0, 0]`` (q=0, self channel) has ``abs(...) ~ 0.30`` at these
    parameters -- an O(1) fraction of the O(1) coupling scale, not a
    round-off-level residual. The assertion below uses a conservative
    floor of ``1e-3`` (orders of magnitude above any Richardson/round-off
    noise floor used elsewhere in this campaign, e.g. adjudicate_
    granule's ``1e-13`` epsilon), so it is robust to the exact measured
    value while still being a genuine, non-vacuous discrimination check.
    """

    def test_mixed_spin_blocks_are_nonzero(self):
        fx = ed_oracle_util.EDFixture(
            L=2, norb=1, t={(0, 0): 0.4 - 0.1j}, eps=(0.05,), T=0.5, mu=0.1)
        h1 = _h1_with_spin_orbit(fx, 0.3 + 0.05j)
        te = ed_oracle_util.TotalNED(fx, h1=h1)
        channels = [(0,), (1,)]
        out = te.bond_correlator_transverse(channels)

        self.assertEqual(out.shape, (fx.L, len(channels), len(channels)))
        magnitudes = np.abs(out)
        max_mag = float(magnitudes.max())
        # Not vacuous: at least one cell clears a floor far above any
        # round-off/noise scale.
        self.assertGreaterEqual(max_mag, 1e-3, out)
        # Every q slice individually carries a nonzero mixed-spin signal
        # (not just one lucky cell) -- the "blocks are independently
        # nonzero" requirement.
        for qi in range(fx.L):
            self.assertGreaterEqual(float(magnitudes[qi].max()), 1e-3,
                                     (qi, out[qi]))


if __name__ == "__main__":
    unittest.main()
