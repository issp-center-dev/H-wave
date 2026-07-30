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


class TestOffsiteVqPreservation(unittest.TestCase):
    """The declaration symmetrisation must not touch off-site structure.

    The partner of an off-site entry (R, (a, b)) is (-R, (b, a)), which in
    momentum space is the orbital transpose AT -q. Averaging with the same-q
    transpose instead collapses V(q) = v e^{-iqR} to v cos(qR) -- at q = pi/2
    the S/C entry vanished outright. Checked at a q with sin(q) != 0 so any
    cos-collapse or conjugation error changes the answer.
    """

    def _sc(self):
        import hwave.sc as sc

        return sc

    def _interactions(self):
        # one-dimensional bond declared from both ends, as an interaction
        # file does; V(q) = 0.7 e^{-iq} must survive verbatim
        return {"CoulombInter": {((1, 0, 0), (0, 1)): 0.7,
                                 ((-1, 0, 0), (1, 0)): 0.7}}

    def test_offsite_vq_is_preserved_with_full_complex_phase(self):
        sc = self._sc()
        kx = np.linspace(0, 2 * np.pi, 4, endpoint=False)
        kz = np.array([0.0])
        ik = sc._build_interaction_k(kx, kx, kz, self._interactions(), 2)
        S, C = sc._build_sc_matrices_all_q(ik, 2, 4, 4, 1)
        q1 = kx[1]                          # pi/2: sin(q) = 1, maximal test
        # builder stores entry (R, (a, b)) at [b, a] with e^{-iqR}: the two
        # cross slots carry 0.7 e^{+iq} (from [0,1]) and 0.7 e^{-iq} ([1,0])
        want = {(1, 1): 0.7 * np.exp(+1j * q1),
                (2, 2): 0.7 * np.exp(-1j * q1)}
        for (i, j), w in want.items():
            for name, M in (("S", S), ("C", C)):
                got = M[1, 0, 0, i, j]
                mag = abs(got)
                self.assertAlmostEqual(mag, 0.7, places=12,
                                       msg="%s[%d,%d] magnitude" % (name, i, j))
                self.assertAlmostEqual(
                    abs(got - (w if name == "S" else -w)), 0.0, places=12,
                    msg="%s[%d,%d] phase" % (name, i, j))

    def test_per_q_builder_matches_all_q_including_negative_index(self):
        sc = self._sc()
        kx = np.linspace(0, 2 * np.pi, 4, endpoint=False)
        kz = np.array([0.0])
        ik = sc._build_interaction_k(kx, kx, kz, self._interactions(), 2)
        S, C = sc._build_sc_matrices_all_q(ik, 2, 4, 4, 1)
        for i in range(4):
            Sp, Cp = sc._build_sc_matrices(ik, 2, i, 0, 0)
            np.testing.assert_allclose(Sp, S[i, 0, 0], atol=1e-13)
            np.testing.assert_allclose(Cp, C[i, 0, 0], atol=1e-13)
        Sm, Cm = sc._build_sc_matrices(ik, 2, -1, 0, 0)   # numpy semantics
        np.testing.assert_allclose(Sm, S[3, 0, 0], atol=1e-13)
        with self.assertRaises(IndexError):
            sc._build_sc_matrices(ik, 2, 7, 0, 0)

    def test_offsite_hermitian_pairhop_phase_is_preserved(self):
        # PairHop's partner is the HERMITIAN entry (-R, (b,a), conj(p)); the
        # same-q conjugated transpose must be an identity for that closure
        # off-site too, keeping the full complex phase
        sc = self._sc()
        kx = np.linspace(0, 2 * np.pi, 4, endpoint=False)
        kz = np.array([0.0])
        p = 0.7 * np.exp(0.3j)
        ik = sc._build_interaction_k(
            kx, kx, kz,
            {"PairHop": {((1, 0, 0), (0, 1)): p,
                         ((-1, 0, 0), (1, 0)): np.conjugate(p)}}, 2)
        raw = ik["PairHop"].copy()
        sym = sc._symmetrise_interactions_k(ik)["PairHop"]
        np.testing.assert_allclose(sym, raw, atol=1e-13)
        # PairHop is stored UNtransposed (measured in #100): the entry
        # (R, (a, b)) lands at [a, b], so [0, 1] carries p e^{-iqR}
        q1 = kx[1]
        self.assertAlmostEqual(
            abs(sym[0, 1, 1, 0, 0] - p * np.exp(-1j * q1)), 0.0, places=12)

    def test_reversal_on_mixed_odd_even_grid(self):
        # (Nx, Ny, Nz) = (5, 4, 3): odd sizes pair i <-> n-i, even sizes make
        # 0 and n/2 self-partners; the both-ends identity must hold per axis
        sc = self._sc()
        kx = np.linspace(0, 2 * np.pi, 5, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, 4, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, 3, endpoint=False)
        ik = sc._build_interaction_k(
            kx, ky, kz,
            {"CoulombInter": {((1, 2, 0), (0, 1)): 0.7,
                              ((-1, -2, 0), (1, 0)): 0.7,
                              ((0, 0, 1), (0, 1)): 0.4,
                              ((0, 0, -1), (1, 0)): 0.4}}, 2)
        raw = ik["CoulombInter"].copy()
        sym = sc._symmetrise_interactions_k(ik)["CoulombInter"]
        np.testing.assert_allclose(sym, raw, atol=1e-13)
        S, C = sc._build_sc_matrices_all_q(ik, 2, 5, 4, 3)
        for idx in ((1, 0, 0), (2, 3, 1), (4, 2, 2), (0, 2, 0)):
            Sp, Cp = sc._build_sc_matrices(ik, 2, *idx)
            np.testing.assert_allclose(Sp, S[idx], atol=1e-13)
            np.testing.assert_allclose(Cp, C[idx], atol=1e-13)

    def test_one_sided_offsite_declaration_gets_the_mean(self):
        # a single-direction declaration means (jab + jba)/2 in this codebase:
        # half the weight at e^{-iq}, half at e^{+iq} on the transposed slot
        sc = self._sc()
        kx = np.linspace(0, 2 * np.pi, 4, endpoint=False)
        kz = np.array([0.0])
        ik = sc._build_interaction_k(
            kx, kx, kz, {"CoulombInter": {((1, 0, 0), (0, 1)): 0.7}}, 2)
        sym = sc._symmetrise_interactions_k(ik)["CoulombInter"]
        q1 = kx[1]
        self.assertAlmostEqual(
            abs(sym[1, 0, 1, 0, 0] - 0.35 * np.exp(-1j * q1)), 0.0, places=12)
        self.assertAlmostEqual(
            abs(sym[0, 1, 1, 0, 0] - 0.35 * np.exp(+1j * q1)), 0.0, places=12)


class TestLegacyFlexFileGuard(unittest.TestCase):
    """Susceptibility files that predate the #113 vertex corrections must not
    be silently paired with the corrected S/C matrices when the interaction
    set contains an affected type (Hund / Exchange / Ising). U/V-only inputs
    are provably unchanged, so legacy files stay usable there. The guard
    requires the version to DECODE to exactly 2 in BOTH files -- a version-1
    tag, a mismatched pair, or a malformed field is as wrong as no tag.
    Activation is by interaction CONTENT: an explicitly configured but empty
    Hund file contributes nothing and must not trip the guard.
    """

    HUND = {((0, 0, 0), (0, 1)): 0.5}

    def _write(self, d, extra_s=None, extra_c=None):
        nd = 4
        arr = np.zeros((4, 4, nd, nd), dtype=complex)
        np.savez(os.path.join(d, "chiq_s.npz"), chiq_s=arr, **(extra_s or {}))
        np.savez(os.path.join(d, "chiq_c.npz"), chiq_c=arr, **(extra_c or {}))
        return {"file": {"output": {"path_to_output": d}}, "eliashberg": {}}

    def _load(self, inp, interactions):
        import hwave.sc as sc

        return sc._load_flex_susceptibilities(
            inp, 2, 2, 2, 1, interactions=interactions)

    def test_legacy_file_with_affected_type_is_rejected(self):
        import tempfile

        d = tempfile.mkdtemp()
        inp = self._write(d)               # no sc_vertex_version
        with self.assertRaises(ValueError) as cm:
            self._load(inp, {"Hund": self.HUND})
        self.assertIn("sc_vertex_version", str(cm.exception))

    def test_legacy_file_with_uv_only_is_accepted(self):
        import tempfile

        d = tempfile.mkdtemp()
        inp = self._write(d)
        self._load(inp, {"CoulombIntra": {(0, 0): 4.0},
                         "CoulombInter": {((0, 0, 0), (0, 1)): 0.7}})

    def test_legacy_file_with_empty_affected_type_is_accepted(self):
        # a configured-but-empty Hund file contributes no interaction, so the
        # vertex content cannot differ; key PRESENCE must not trip the guard
        import tempfile

        d = tempfile.mkdtemp()
        inp = self._write(d)
        self._load(inp, {"Hund": {}, "Exchange": {}, "Ising": {}})

    def test_tagged_file_with_affected_type_is_accepted(self):
        import tempfile

        d = tempfile.mkdtemp()
        v = {"sc_vertex_version": 2}
        inp = self._write(d, extra_s=v, extra_c=v)
        self._load(inp, {"Exchange": self.HUND})

    def test_version_1_is_rejected(self):
        import tempfile

        d = tempfile.mkdtemp()
        v = {"sc_vertex_version": 1}
        inp = self._write(d, extra_s=v, extra_c=v)
        with self.assertRaises(ValueError) as cm:
            self._load(inp, {"Hund": self.HUND})
        self.assertIn("version 2", str(cm.exception))

    def test_mismatched_versions_are_rejected(self):
        import tempfile

        d = tempfile.mkdtemp()
        inp = self._write(d, extra_s={"sc_vertex_version": 2},
                          extra_c={"sc_vertex_version": 1})
        with self.assertRaises(ValueError) as cm:
            self._load(inp, {"Ising": self.HUND})
        self.assertIn("different", str(cm.exception))

    def test_fractional_version_is_rejected_not_truncated(self):
        # int() would truncate 2.9 to 2 and accept it
        import tempfile

        d = tempfile.mkdtemp()
        v = {"sc_vertex_version": 2.9}
        inp = self._write(d, extra_s=v, extra_c=v)
        with self.assertRaises(ValueError) as cm:
            self._load(inp, {"Hund": self.HUND})
        self.assertIn("malformed", str(cm.exception))

    def test_nonfinite_version_is_rejected_as_valueerror(self):
        # int(inf) raises OverflowError, which must not escape undocumented
        import tempfile

        d = tempfile.mkdtemp()
        v = {"sc_vertex_version": np.inf}
        inp = self._write(d, extra_s=v, extra_c=v)
        with self.assertRaises(ValueError) as cm:
            self._load(inp, {"Hund": self.HUND})
        self.assertIn("malformed", str(cm.exception))

    def test_nonscalar_version_is_rejected(self):
        import tempfile

        d = tempfile.mkdtemp()
        v = {"sc_vertex_version": np.array([2, 2])}
        inp = self._write(d, extra_s=v, extra_c=v)
        with self.assertRaises(ValueError) as cm:
            self._load(inp, {"Hund": self.HUND})
        self.assertIn("malformed", str(cm.exception))

    def test_integral_float_version_is_accepted(self):
        # 2.0 decodes to exactly 2; only truncation is forbidden
        import tempfile

        d = tempfile.mkdtemp()
        v = {"sc_vertex_version": 2.0}
        inp = self._write(d, extra_s=v, extra_c=v)
        self._load(inp, {"Hund": self.HUND})

    def test_full_frequency_loader_activates_the_same_guard(self):
        # the dynamic-Eliashberg path loads through
        # _load_flex_susceptibilities_full, which must forward the interaction
        # content into the same version guard as the static loader
        import tempfile
        import hwave.sc as sc

        d = tempfile.mkdtemp()
        inp = self._write(d)               # untagged
        with self.assertRaises(ValueError) as cm:
            sc._load_flex_susceptibilities_full(
                inp, 2, 2, 2, 1, interactions={"Hund": self.HUND})
        self.assertIn("sc_vertex_version", str(cm.exception))

    def test_malformed_version_is_rejected(self):
        import tempfile

        d = tempfile.mkdtemp()
        v = {"sc_vertex_version": np.array(["two"])}
        inp = self._write(d, extra_s=v, extra_c=v)
        with self.assertRaises(ValueError) as cm:
            self._load(inp, {"Hund": self.HUND})
        self.assertIn("malformed", str(cm.exception))

if __name__ == "__main__":
    unittest.main()
