#!/usr/bin/env python3

"""Contract tests and the H1 calibration gate for
``ed_oracle_util.SectorED.bond_correlator_transverse`` (Phase H of the
bond-transverse campaign, issue #110's validation groundwork).

See ``docs/superpowers/specs/2026-08-15-bond-transverse-design.md``
("Phase H -- the transverse ED harness") for the binding contract this
module pins, and ``docs/superpowers/plans/2026-08-15-bond-transverse-
phase-h.md`` for the task-by-task build order. This module (Task 1)
covers:

- the operator-construction contract (shape/dtype, R-wrapping duplicate
  rejection, the mixed-spin sector delta) via
  ``TestBondCorrelatorTransverseContract``;
- H1: the R == 0 channels of ``bond_correlator_transverse`` reproduce the
  bespoke ``_chi_pm_ed`` dense-ED reference from
  ``tests/test_rpa_vs_ed_oracle.py::TestTransverseComplexPairHop`` on the
  SAME 2-site spinful fixture, at round-off (both sides are exact ED on
  the same Hilbert space -- ``assert_approx_array(rel=0, abs=1e-12)``).

H2 (on-site 7-type re-adjudication) and H3 (the frame-map derivation and
the V=0 bond pin) are later tasks in this same campaign and are not
covered here.

Tests must be run from the repository root.
"""

import unittest

import numpy as np

import tests.test_rpa_vs_ed_oracle as oracle
from tests import ed_oracle_util
from tests.approx_util import assert_approx_array

# NOTE: imported as a module (``oracle``), never as bare names -- binding
# ``oracle.TestTransverseComplexPairHop`` at module scope here would make
# unittest's discovery pick it up as a SECOND, duplicate TestCase under
# this module too.


class TestBondCorrelatorTransverseContract(unittest.TestCase):
    """Operator-construction contract, independent of any physics
    comparison: shape/dtype, R-wrapping duplicate rejection, and the
    mixed-spin sector displacement the new method relies on."""

    @staticmethod
    def _small_fx():
        # 2-site, two-orbital, spinful -- fast, disposable fixture for the
        # contract tests (distinct from the H1 fixture below, which must
        # be the EXACT 2-site/2-orbital TestTransverseComplexPairHop
        # fixture). norb=2 so channels use the (R, a, b) 3-tuple form
        # ``_channel_orbitals`` expects (the norb == 1 bare-``(R,)``
        # shorthand is not what the brief's example channels use).
        t = {(0, 0): 0.3 + 0.1j, (1, 1): 0.2 - 0.05j,
             (0, 1): 0.1j, (1, 0): 0.1j}
        return ed_oracle_util.EDFixture(
            L=2, norb=2, t=t, eps=(0.05, -0.02), T=0.5, mu=0.1)

    def test_output_shape_and_dtype(self):
        fx = self._small_fx()
        se = ed_oracle_util.SectorED(fx)
        channels = [(0, 0, 0), (1, 0, 0)]
        out = se.bond_correlator_transverse(channels)
        self.assertEqual(out.shape, (fx.L, len(channels), len(channels)))
        self.assertEqual(out.dtype, np.complex128)

    def test_duplicate_channels_after_wrapping_rejected(self):
        fx = self._small_fx()
        se = ed_oracle_util.SectorED(fx)
        # (1, 0, 0) and (1 - L, 0, 0) both wrap (mod L=2) to R=1 -- same
        # channel after wrapping, must be rejected rather than silently
        # double-counted.
        channels = [(1, 0, 0), (1 - fx.L, 0, 0)]
        with self.assertRaises(ValueError):
            se.bond_correlator_transverse(channels)

    def test_mixed_spin_delta_is_plus_minus_one(self):
        # Probe the sector displacement via the same operator-construction
        # API bond_correlator_transverse uses internally: a bilinear whose
        # creation leg is spin UP (mode parity 0) and annihilation leg is
        # spin DOWN (mode parity 1) must impart delta = (+1, -1) -- the
        # up-sector gains one particle, the down-sector loses one.
        fx = self._small_fx()
        se = ed_oracle_util.SectorED(fx)
        mode_terms = [(fx.mode(j, 0, 0), fx.mode(j, 0, 1), 1.0)
                      for j in range(fx.L)]
        delta, _op = se._build_operator(mode_terms)
        self.assertEqual(delta, (1, -1))


class TestH1AgainstBespokeChiPmEd(unittest.TestCase):
    """H1: on R == 0 channels, ``bond_correlator_transverse`` must
    reproduce ``TestTransverseComplexPairHop._chi_pm_ed`` exactly (both
    are exact diagonalization of the SAME Hamiltonian on the SAME 2-site,
    2-orbital spinful Hilbert space -- round-off agreement, not a
    physics-tolerance one)."""

    def test_r0_matches_bespoke_chi_pm_ed(self):
        fx, _C, _CD, h1 = oracle._fx2_state()
        se = ed_oracle_util.SectorED(fx, h1=h1)

        # channels ordered (a outer, b inner) so the flat (I, J) axes
        # reshape directly onto the bespoke reference's (a, c, b, d) axes
        # below: channel (0, a, b) <-> creation orbital a (up),
        # annihilation orbital b (down).
        channels = [(0, a, b) for a in range(oracle.NORB)
                    for b in range(oracle.NORB)]
        got = se.bond_correlator_transverse(channels)
        got = got.reshape(oracle.LX, oracle.NORB, oracle.NORB,
                           oracle.NORB, oracle.NORB)

        bespoke = oracle.TestTransverseComplexPairHop()
        expected = bespoke._chi_pm_ed()  # base Hamiltonian, no interaction

        assert_approx_array(got, expected, rel=0, abs=1e-12)


if __name__ == "__main__":
    unittest.main()
