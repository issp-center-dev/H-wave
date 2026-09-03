#!/usr/bin/env python3
"""The locality-split S/C builder used by the general FLEX path (#181, Tier 1).

`hwave.sc._build_sc_matrices_all_q` applies the on-site Kanamori slot map to
whatever k-space interaction it is handed, with no locality bookkeeping: an
off-site entry then lands in the cross (ab,ab) / antidiag (ab,ba) slot
families, whose particle-hole pair is NON-local for R != 0 and not a function
of q alone; and the density-family gate `l1 != l3` (right for on-site Hund /
Ising) deletes the physical same-orbital off-site Hund / Ising.

`build_sc_matrices_locality_split` takes the on-site and off-site parts of the
interaction SEPARATELY (split on the PRE-fold table by the caller): the
on-site part goes through the full slot map unchanged; the off-site part is
written into the density family only (Hartree vertex V(q) on (aa,bb),
including a == b), which is exactly what the RPA ring carries for it. The
exchange (Fock) crossing of an off-site term is NOT representable here and is
deliberately absent (bond-resolved treatment: #181 Tier 3).
"""

import unittest

import numpy as np


def _grid(n):
    k = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    return k, k, np.array([0.0])


def _bond(a, b, v, R=(1, 0, 0)):
    """Hermitian-closed +-R bond between orbitals a and b."""
    mR = tuple(-x for x in R)
    return {(R, (a, b)): v, (mR, (b, a)): v}


class TestLocalitySplitBuilder(unittest.TestCase):

    def _build(self, onsite, offsite, norb, n=4):
        import hwave.sc as sc
        from hwave.solver._sc_matrices_myo import (
            build_sc_matrices_locality_split)

        kx, ky, kz = _grid(n)
        ik_on = sc._build_interaction_k(kx, ky, kz, onsite, norb)
        ik_off = sc._build_interaction_k(kx, ky, kz, offsite, norb)
        return build_sc_matrices_locality_split(ik_on, ik_off, norb,
                                                n, n, 1)

    def test_onsite_only_input_is_bitwise_the_full_builder(self):
        """With no off-site part the split builder IS the existing one:
        every on-site slot family (diag, cross, density, antidiag) keeps
        its adjudicated content bit for bit."""
        import hwave.sc as sc

        onsite = {
            "CoulombIntra": {((0, 0, 0), (0, 0)): 2.0, ((0, 0, 0), (1, 1)): 1.5},
            "CoulombInter": {((0, 0, 0), (0, 1)): 0.8, ((0, 0, 0), (1, 0)): 0.8},
            "Hund": {((0, 0, 0), (0, 1)): 0.3, ((0, 0, 0), (1, 0)): 0.3},
            "Exchange": {((0, 0, 0), (0, 1)): 0.3, ((0, 0, 0), (1, 0)): 0.3},
            "Ising": {((0, 0, 0), (0, 1)): 0.1, ((0, 0, 0), (1, 0)): 0.1},
            "PairHop": {((0, 0, 0), (0, 1)): 0.2, ((0, 0, 0), (1, 0)): 0.2},
        }
        kx, ky, kz = _grid(4)
        ik = sc._build_interaction_k(kx, ky, kz, onsite, 2)
        S_ref, C_ref = sc._build_sc_matrices_all_q(ik, 2, 4, 4, 1)
        S, C = self._build(onsite, {}, 2)
        self.assertTrue(np.array_equal(S, S_ref))
        self.assertTrue(np.array_equal(C, C_ref))
        # anti-vacuity: every family is populated in this fixture
        for slot in ((0, 0), (1, 1), (0, 3), (1, 2)):
            self.assertNotEqual(np.abs(S_ref[..., slot[0], slot[1]]).max(), 0.0)

    def test_offsite_content_lands_on_density_slots_only(self):
        """Off-site inter-orbital CoulombInter / Hund / Ising: the density
        slots (aa,bb) carry the adjudicated (S, C) per unit coupling times
        V_ab(q); the cross and antidiag slots stay exactly zero."""
        from hwave.solver.vertex_table import sc_coefficients

        n, v = 4, 0.2
        kx = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
        # The declaration (+x, (0, 1)) is the operator n_0(i) n_1(i+x); its
        # reversal partner (-x, (1, 0)) is the same operator. The k-space
        # matrix stores the orbital pair TRANSPOSED (hwave.sc's
        # _build_interaction_k, issue #96), so the (00,11) density slot
        # reads mat[0, 1] = v e^{-iqx} and (11,00) reads mat[1, 0] =
        # v e^{+iqx} -- the same orientation the shared builder's Case 3
        # uses for on-site input, and the ring's, measured element-complete
        # equal end to end (tests/test_flex_offsite_general.py).
        vq = {(0, 3): v * np.exp(-1j * kx), (3, 0): v * np.exp(+1j * kx)}
        for itype in ("CoulombInter", "Hund", "Ising"):
            with self.subTest(itype=itype):
                S, C = self._build({}, {itype: _bond(0, 1, v)}, 2, n=n)
                s, c = sc_coefficients(itype, "density")
                # density slots: (00,11) and (11,00); pair index a*norb+a
                for (i, j), vq_ij in vq.items():
                    np.testing.assert_allclose(S[:, 0, 0, i, j], s * vq_ij,
                                              atol=1e-15)
                    np.testing.assert_allclose(C[:, 0, 0, i, j], c * vq_ij,
                                              atol=1e-15)
                # cross (01,01), (10,10) and antidiag (01,10), (10,01)
                for (i, j) in ((1, 1), (2, 2), (1, 2), (2, 1)):
                    self.assertEqual(np.abs(S[..., i, j]).max(), 0.0)
                    self.assertEqual(np.abs(C[..., i, j]).max(), 0.0)
                # same-orbital density slots untouched by an a != b bond
                for (i, j) in ((0, 0), (3, 3)):
                    self.assertEqual(np.abs(S[..., i, j]).max(), 0.0)
                    self.assertEqual(np.abs(C[..., i, j]).max(), 0.0)

    def test_offsite_same_orbital_hund_and_ising_are_kept(self):
        """The on-site density gate (l1 != l3) must NOT apply off-site: a
        one-orbital off-site Hund / Ising is a physical density-density
        term and lands on the (00,00) slot with the density (S, C)."""
        from hwave.solver.vertex_table import sc_coefficients

        n, v = 4, 0.5
        kx = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
        vq = 2.0 * v * np.cos(kx)
        for itype in ("Hund", "Ising", "CoulombInter"):
            with self.subTest(itype=itype):
                S, C = self._build({}, {itype: _bond(0, 0, v)}, 1, n=n)
                s, c = sc_coefficients(itype, "density")
                np.testing.assert_allclose(S[:, 0, 0, 0, 0], s * vq, atol=1e-15)
                np.testing.assert_allclose(C[:, 0, 0, 0, 0], c * vq, atol=1e-15)
                self.assertNotEqual(abs(s) + abs(c), 0.0)

    def test_offsite_one_sided_declaration_reads_as_its_symmetric_part(self):
        """A one-sided off-site table (+R only, internal tables) is reduced
        to the same symmetric coefficient the on-site builder and the ring
        use: v e^{iqR} -> v cos(qR)."""
        n, v = 4, 0.4
        kx = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
        S, C = self._build({}, {"CoulombInter": {((1, 0, 0), (0, 0)): v}},
                           1, n=n)
        np.testing.assert_allclose(C[:, 0, 0, 0, 0], 2.0 * v * np.cos(kx),
                                  atol=1e-15)

    def test_offsite_exchange_or_pairhop_in_the_offsite_part_fails_closed(self):
        """Types with no density-family content (Exchange, PairHop,
        PairLift) and CoulombIntra have NO q-only off-site representation
        in this builder; handing them in the off-site part is a caller
        error and must raise rather than be silently dropped."""
        for itype in ("Exchange", "PairHop", "PairLift", "CoulombIntra"):
            with self.subTest(itype=itype):
                with self.assertRaises(ValueError) as cm:
                    self._build({}, {itype: _bond(0, 1, 0.2)}, 2)
                self.assertIn(itype, str(cm.exception))

    def test_onsite_and_offsite_parts_add(self):
        """Mixed input: on-site Kanamori plus an off-site bond of the same
        type; the result is the sum of the two single-part builds."""
        onsite = {"CoulombInter": {((0, 0, 0), (0, 1)): 0.8,
                                   ((0, 0, 0), (1, 0)): 0.8},
                  "Hund": {((0, 0, 0), (0, 1)): 0.3, ((0, 0, 0), (1, 0)): 0.3}}
        offsite = {"CoulombInter": _bond(0, 1, 0.2), "Hund": _bond(1, 1, 0.1)}
        S, C = self._build(onsite, offsite, 2)
        S1, C1 = self._build(onsite, {}, 2)
        S2, C2 = self._build({}, offsite, 2)
        np.testing.assert_allclose(S, S1 + S2, atol=0.0, rtol=0.0)
        np.testing.assert_allclose(C, C1 + C2, atol=0.0, rtol=0.0)
        self.assertNotEqual(np.abs(S2).max(), 0.0)


if __name__ == "__main__":
    unittest.main()
