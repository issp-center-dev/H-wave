#!/usr/bin/env python3
"""The adjudicated S/C vertex table, pinned per interaction type.

Every value below was established by exact diagonalization of H-wave's
documented file operators (issue #113; methodology in issue #107): the
O(lambda) self-energy is removed exactly via the Hartree-Fock response, the
extraction is validated by q-independence for on-site input, the label map is
anchored computationally on the end-to-end-adjudicated CoulombInter case, and
the per-channel scale is calibrated once on the CoulombIntra control.

S enters the spin channel as [1 - chi0 S]^-1 chi0 and C the charge channel as
[1 + chi0 C]^-1 chi0. Slots are (l1*norb+l2, l3*norb+l4).

Internal consistency checks carried by the table itself:
  * SU(2): for Hund + Exchange at equal J (the rotationally symmetric
    combination), S(ab,ab) has no J at all and C(ab,ab) = -U' + 2J -- the
    standard Kanamori literature values, reproduced from the corrected
    per-type split (the old bookkeeping assigned the whole 2J to Hund and
    nothing to Exchange, which is also why the "MYO vs Kuroki convention
    difference" dissolved: it was a per-type attribution error).
  * The transverse channel measured in the #105 series is consistent with
    this table through the SU(2) identity for the same combination.
"""

import unittest

import numpy as np


def _single(itype, val, norb=2):
    import hwave.sc as sc

    k = np.array([0.0])
    if itype == "CoulombIntra":
        entries = {((0, 0, 0), (a, a)): val for a in range(norb)}
    else:
        entries = {((0, 0, 0), (0, 1)): val, ((0, 0, 0), (1, 0)): val}
    ik = sc._build_interaction_k(k, k, k, {itype: entries}, norb)
    S, C = sc._build_sc_matrices_all_q(ik, norb, 1, 1, 1)
    return S[0, 0, 0], C[0, 0, 0]


class TestAdjudicatedVertexTable(unittest.TestCase):

    def test_per_type_sc_content(self):
        J = 0.7
        # slots: 0 = (0,0), 1 = (0,1), 2 = (1,0), 3 = (1,1)
        expected = {
            #            S nonzero slots            C nonzero slots
            "CoulombIntra": ({(0, 0): J, (3, 3): J}, {(0, 0): J, (3, 3): J}),
            "CoulombInter": ({(1, 1): J, (2, 2): J},
                             {(1, 1): -J, (2, 2): -J,
                              (0, 3): 2*J, (3, 0): 2*J}),
            "Hund":         ({(0, 3): J, (3, 0): J,
                              (1, 1): -J, (2, 2): -J},
                             {(0, 3): -J, (3, 0): -J,
                              (1, 1): J, (2, 2): J}),
            "Exchange":     ({(1, 1): J, (2, 2): J},
                             {(1, 1): J, (2, 2): J}),
            "Ising":        ({(0, 3): -2*J, (3, 0): -2*J,
                              (1, 1): J, (2, 2): J},
                             {(1, 1): -J, (2, 2): -J}),
            "PairLift":     ({}, {}),
            "PairHop":      ({(1, 2): J, (2, 1): J},
                             {(1, 2): J, (2, 1): J}),
        }
        for itype, (s_want, c_want) in expected.items():
            with self.subTest(interaction=itype):
                S, C = _single(itype, J)
                for name, M, want in (("S", S, s_want), ("C", C, c_want)):
                    built = np.zeros((4, 4), dtype=complex)
                    for (i, j), v in want.items():
                        built[i, j] = v
                    np.testing.assert_allclose(
                        M, built, atol=1e-12,
                        err_msg="%s matrix of %s does not match the "
                                "adjudicated table" % (name, itype))

    def test_su2_combination_recovers_the_standard_kanamori_values(self):
        """Hund + Exchange at equal J: S(ab,ab) carries no J and
        C(ab,ab) = +2J (with U' = 0 here) -- the standard literature values,
        obtained from the corrected per-type split."""
        import hwave.sc as sc

        J = 0.7
        k = np.array([0.0])
        both = {"Hund": {((0, 0, 0), (0, 1)): J, ((0, 0, 0), (1, 0)): J},
                "Exchange": {((0, 0, 0), (0, 1)): J, ((0, 0, 0), (1, 0)): J}}
        ik = sc._build_interaction_k(k, k, k, both, 2)
        S, C = sc._build_sc_matrices_all_q(ik, 2, 1, 1, 1)
        self.assertAlmostEqual(S[0, 0, 0, 1, 1].real, 0.0, places=12)
        self.assertAlmostEqual(C[0, 0, 0, 1, 1].real, 2 * J, places=12)


if __name__ == "__main__":
    unittest.main()
