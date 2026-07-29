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

import os
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




class TestDeclarationSymmetrisation(unittest.TestCase):
    """The S/C builders symmetrise the two redundant declarations of each
    orbital pair with the mean (plain transpose for every type whose two
    entries multiply the same operator; conjugated for PairHop, whose entries
    are Hermitian partners and whose complex phase is physical, #100/#105).

    Without this, an asymmetric Exchange declaration put unequal values at the
    two (ab,ab) slots -- diverging from the transverse channel, which
    symmetrises -- and an antisymmetric declaration, an identically zero
    Hamiltonian, produced a nonzero vertex.
    """

    def _sc(self, itype, entries):
        import hwave.sc as sc

        k = np.array([0.0])
        ik = sc._build_interaction_k(k, k, k, {itype: entries}, 2)
        S, C = sc._build_sc_matrices_all_q(ik, 2, 1, 1, 1)
        return S[0, 0, 0], C[0, 0, 0]

    def test_antisymmetric_declarations_give_zero(self):
        entries = {((0, 0, 0), (0, 1)): 1.0, ((0, 0, 0), (1, 0)): -1.0}
        for itype in ("Exchange", "CoulombInter", "Ising", "Hund"):
            with self.subTest(interaction=itype):
                S, C = self._sc(itype, entries)
                self.assertLess(float(np.max(np.abs(S))), 1e-12)
                self.assertLess(float(np.max(np.abs(C))), 1e-12)

    def test_asymmetric_exchange_uses_the_mean(self):
        S, C = self._sc("Exchange",
                        {((0, 0, 0), (0, 1)): 1.0, ((0, 0, 0), (1, 0)): 0.35})
        mean = (1.0 + 0.35) / 2
        for M in (S, C):
            self.assertAlmostEqual(M[1, 1].real, mean, places=12)
            self.assertAlmostEqual(M[2, 2].real, mean, places=12)

    def test_complex_hermitian_pairhop_keeps_its_phase(self):
        p = 0.7 + 0.4j
        S, C = self._sc("PairHop",
                        {((0, 0, 0), (0, 1)): p,
                         ((0, 0, 0), (1, 0)): np.conj(p)})
        # Hermitian-closed input: the conjugated mean is an identity, so the
        # full complex value survives (a plain mean would collapse it to Re p)
        self.assertGreater(float(np.max(np.abs(S.imag))), 0.39)


class TestLegacyFlexFileGuard(unittest.TestCase):
    """Susceptibility files that predate the #113 vertex corrections must not
    be silently paired with the corrected S/C matrices when the interaction
    set contains an affected type (Hund / Exchange / Ising). U/V-only inputs
    are provably unchanged, so legacy files stay usable there.
    """

    def _write(self, d, **extra):
        nd = 4
        arr = np.zeros((4, 4, nd, nd), dtype=complex)
        np.savez(os.path.join(d, "chiq_s.npz"), chiq_s=arr, **extra)
        np.savez(os.path.join(d, "chiq_c.npz"), chiq_c=arr, **extra)
        return {"file": {"output": {"path_to_output": d}}, "eliashberg": {}}

    def test_legacy_file_with_affected_type_is_rejected(self):
        import tempfile
        import hwave.sc as sc

        d = tempfile.mkdtemp()
        inp = self._write(d)               # no sc_vertex_version
        with self.assertRaises(ValueError) as cm:
            sc._load_flex_susceptibilities(
                inp, 2, 2, 2, 1, interactions={"Hund": {}})
        self.assertIn("sc_vertex_version", str(cm.exception))

    def test_legacy_file_with_uv_only_is_accepted(self):
        import tempfile
        import hwave.sc as sc

        d = tempfile.mkdtemp()
        inp = self._write(d)
        sc._load_flex_susceptibilities(
            inp, 2, 2, 2, 1,
            interactions={"CoulombIntra": {}, "CoulombInter": {}})

    def test_tagged_file_with_affected_type_is_accepted(self):
        import tempfile
        import hwave.sc as sc

        d = tempfile.mkdtemp()
        inp = self._write(d, sc_vertex_version=2)
        sc._load_flex_susceptibilities(
            inp, 2, 2, 2, 1, interactions={"Exchange": {}})


if __name__ == "__main__":
    unittest.main()
