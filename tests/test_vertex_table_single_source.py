#!/usr/bin/env python3
"""Both vertex builders derive from the one adjudicated table.

The pair-space S/C matrices (hwave.sc) and the ring's spin-resolved
cross-slot correction (hwave.solver.rpa) encode the same physics in two
representations, tied by the channel decomposition W_same = (C - S)/2,
W_cross = (C + S)/2. Historically the two drifted (#96, #104, #113);
they now read one table (hwave.solver.vertex_table). These tests measure
the consistency from the ASSEMBLED objects, not the table -- so a future
edit that bypasses the table on either side still fails here.
"""

import unittest

import numpy as np


class TestAssembledConsistency(unittest.TestCase):

    def _sc_cross(self, itype, J=0.7):
        """S and C at the (01, 01) cross slot from the sc builders."""
        import hwave.sc as sc

        k = np.array([0.0])
        entries = {((0, 0, 0), (0, 1)): J, ((0, 0, 0), (1, 0)): J}
        ik = sc._build_interaction_k(k, k, k, {itype: entries}, 2)
        S, C = sc._build_sc_matrices_all_q(ik, 2, 1, 1, 1)
        return S[0, 0, 0, 1, 1], C[0, 0, 0, 1, 1]

    def _rpa_fierz(self, itype, J=0.7):
        """Same-spin and cross-spin cross-slot entries of the ring's
        Fierz tensor, read from the assembled solver."""
        import hwave.solver.rpa as rpa_mod

        param_ham = {'Geometry': {'norb': 2, 'rvec': np.eye(3),
                                  'center': [[0, 0, 0]] * 2},
                     'Transfer': {((0, 0, 0), (0, 0)): 1.0},
                     itype: {((0, 0, 0), (0, 1)): J,
                             ((0, 0, 0), (1, 0)): J}}
        info = {'mode': 'RPA', 'param': {'T': 2.0, 'filling': 0.5,
                'CellShape': [1, 1, 1], 'SubShape': [1, 1, 1], 'Nmat': 4},
                'enable_spin_orbital': False, 'calc_scheme': 'general',
                'calc_type': 'ring'}
        sv = rpa_mod.RPA(param_ham, {}, info)
        f = sv.ham_info.ham_fierz_q
        if f is None:
            return 0.0, 0.0
        f = np.asarray(f).reshape(-1, 4, 4, 4, 4)[0]
        # spin-major nd = 4: up = {0: 0, 1: 1}, down = {0: 2, 1: 3};
        # entries sit at (a, b, a, b) per spin block (see the appender)
        same = f[0, 1, 0, 1]          # (up a, up b, up a, up b)
        cross = f[2, 3, 0, 1]         # (dn a, dn b, up a, up b)
        return same, cross

    def test_ring_fierz_is_the_channel_decomposition_of_sc(self):
        """For every type with cross-slot content: measure S/C from the
        sc builders and the same-spin/cross-spin entries from the ring,
        and require W_same == (C - S)/2, W_cross == (C + S)/2."""
        J = 0.7
        expected_sc = {"CoulombInter": (+J, -J), "Hund": (-J, +J),
                       "Ising": (+J, -J), "Exchange": (+J, +J)}
        for itype in ("CoulombInter", "Hund", "Ising", "Exchange"):
            with self.subTest(interaction=itype):
                S, C = self._sc_cross(itype, J=J)
                # guard against a vacuous pass: if both builders dropped
                # the cross content entirely, the identity below would
                # hold on zeros -- pin the expected NONZERO values first
                want_s, want_c = expected_sc[itype]
                self.assertAlmostEqual(complex(S), complex(want_s),
                                       places=12, msg="%s S" % itype)
                self.assertAlmostEqual(complex(C), complex(want_c),
                                       places=12, msg="%s C" % itype)
                same, cross = self._rpa_fierz(itype, J=J)
                self.assertAlmostEqual(
                    complex(same), complex((C - S) / 2.0), places=12,
                    msg="%s same-spin" % itype)
                self.assertAlmostEqual(
                    complex(cross), complex((C + S) / 2.0), places=12,
                    msg="%s cross-spin" % itype)

    def test_types_without_cross_content_have_no_fierz(self):
        for itype in ("PairHop", "PairLift", "CoulombIntra"):
            with self.subTest(interaction=itype):
                entries = ({((0, 0, 0), (0, 0)): 0.7,
                            ((0, 0, 0), (1, 1)): 0.7}
                           if itype == "CoulombIntra" else
                           {((0, 0, 0), (0, 1)): 0.7,
                            ((0, 0, 0), (1, 0)): 0.7})
                import hwave.solver.rpa as rpa_mod

                param_ham = {'Geometry': {'norb': 2, 'rvec': np.eye(3),
                                          'center': [[0, 0, 0]] * 2},
                             'Transfer': {((0, 0, 0), (0, 0)): 1.0},
                             itype: entries}
                info = {'mode': 'RPA',
                        'param': {'T': 2.0, 'filling': 0.5,
                                  'CellShape': [1, 1, 1],
                                  'SubShape': [1, 1, 1], 'Nmat': 4},
                        'enable_spin_orbital': False,
                        'calc_scheme': 'general', 'calc_type': 'ring'}
                sv = rpa_mod.RPA(param_ham, {}, info)
                self.assertIsNone(sv.ham_info.ham_fierz_q)


class TestTableContent(unittest.TestCase):

    def test_table_matches_the_adjudicated_values(self):
        """Pin the exported table itself, independent of the builders."""
        from hwave.solver.vertex_table import ADJUDICATED_SC

        self.assertEqual(ADJUDICATED_SC, {
            "CoulombIntra": {"diag": (+1.0, +1.0)},
            "CoulombInter": {"cross": (+1.0, -1.0), "density": (0.0, +2.0)},
            "Hund":         {"cross": (-1.0, +1.0), "density": (+1.0, -1.0)},
            "Exchange":     {"cross": (+1.0, +1.0)},
            "Ising":        {"cross": (+1.0, -1.0), "density": (-2.0, 0.0)},
            "PairLift":     {},
            "PairHop":      {"antidiag": (+1.0, +1.0)},
        })


class TestRingBaseContent(unittest.TestCase):
    """The ring's base spin tables are now DERIVED from the adjudicated
    table (density slots via the channel decomposition, spin-flip slots
    from the #105 transverse data). Pin the derived output against the
    hand-coded per-type tables the ring carried since its adjudication --
    these literals ARE the adjudicated base content and must never change
    as a side effect of table or derivation edits."""

    def test_derived_tables_equal_the_adjudicated_hand_coded_ones(self):
        from hwave.solver.vertex_table import ring_spin_table

        adjudicated = {
            'CoulombIntra': {(0, 0, 1, 1): 1, (1, 1, 0, 0): 1},
            'CoulombInter': {(0, 0, 0, 0): 1, (1, 1, 1, 1): 1,
                             (0, 0, 1, 1): 1, (1, 1, 0, 0): 1},
            'Hund':         {(0, 0, 0, 0): -1, (1, 1, 1, 1): -1},
            'Ising':        {(0, 0, 0, 0): 1, (1, 1, 1, 1): 1,
                             (0, 0, 1, 1): -1, (1, 1, 0, 0): -1},
            'PairLift':     {(0, 1, 0, 1): 1, (1, 0, 1, 0): 1},
            'Exchange':     {(0, 1, 1, 0): -1, (1, 0, 0, 1): -1},
            'PairHop':      {(0, 0, 1, 1): 1, (1, 1, 0, 0): 1},
        }
        for itype, want in adjudicated.items():
            with self.subTest(interaction=itype):
                got = ring_spin_table(itype)
                self.assertEqual(got, want)
                # same keys AND exactly representable weights: the derived
                # floats must multiply bit-identically to the old ints
                for k, w in got.items():
                    self.assertEqual(float(w), float(want[k]))

    def test_unknown_type_yields_empty_table(self):
        from hwave.solver.vertex_table import ring_spin_table

        self.assertEqual(ring_spin_table("NoSuchType"), {})

    def test_decomposition_is_algebraic_not_a_seven_case_accident(self):
        """A synthetic type with non-integral density content (S, C) =
        (x, y) must decompose to W_same = (y - x)/2, W_cross = (y + x)/2
        by construction -- the derivation is the algebraic channel
        identity, not something that merely reproduces the seven
        integer-valued adjudicated tables."""
        import hwave.solver.vertex_table as vt
        from types import MappingProxyType
        from unittest import mock

        x, y = -3.25, 1.75
        synthetic = dict(vt.ADJUDICATED_SC)
        synthetic["Synthetic"] = MappingProxyType({"density": (x, y)})
        with mock.patch.object(vt, "ADJUDICATED_SC",
                               MappingProxyType(synthetic)):
            got = vt.ring_spin_table("Synthetic")
        self.assertEqual(got, {
            (0, 0, 0, 0): (y - x) / 2.0, (1, 1, 1, 1): (y - x) / 2.0,
            (0, 0, 1, 1): (y + x) / 2.0, (1, 1, 0, 0): (y + x) / 2.0,
        })

    def test_spin_flip_data_is_frozen(self):
        import hwave.solver.vertex_table as vt

        with self.assertRaises(TypeError):
            vt.RING_SPIN_FLIP["Exchange"] = {}
        with self.assertRaises(TypeError):
            vt.RING_SPIN_FLIP["Exchange"][(0, 1, 1, 0)] = 0.0
        self.assertFalse(hasattr(vt, "_RING_SPIN_FLIP_RAW"),
                         "no mutable backing name may be retained")
        # ring_spin_table returns a fresh dict; mutating it must not leak
        t = vt.ring_spin_table("Exchange")
        t[(0, 1, 1, 0)] = 99.0
        self.assertEqual(vt.ring_spin_table("Exchange")[(0, 1, 1, 0)], -1.0)


class TestTableImmutability(unittest.TestCase):

    def test_both_mapping_levels_are_frozen(self):
        """The table-equality test would still pass if the proxies were
        removed; pin the freeze itself."""
        import hwave.solver.vertex_table as vt

        with self.assertRaises(TypeError):
            vt.ADJUDICATED_SC["Hund"] = {}
        with self.assertRaises(TypeError):
            vt.ADJUDICATED_SC["Hund"]["cross"] = (0.0, 0.0)
        self.assertFalse(hasattr(vt, "_ADJUDICATED_SC_RAW"),
                         "no mutable backing name may be retained")
        self.assertIsInstance(vt.sc_coefficients("Hund", "cross"), tuple)
        self.assertIsInstance(vt.sc_coefficients("Hund", "nope"), tuple)


class TestZeroCoefficientSuppression(unittest.TestCase):

    def test_nonfinite_input_does_not_leak_into_absent_channels(self):
        """0.0 * Inf is NaN: a channel with no content for a type must
        stay exactly zero even for non-finite couplings. Diagonal
        CoulombInter has no S content; density-slot Ising has no C
        content."""
        import warnings

        import hwave.sc as sc

        ctx = warnings.catch_warnings()
        ctx.__enter__()
        self.addCleanup(ctx.__exit__, None, None, None)
        warnings.simplefilter("ignore", RuntimeWarning)

        k = np.array([0.0])
        ik = sc._build_interaction_k(
            k, k, k, {"CoulombInter": {((0, 0, 0), (0, 0)): np.inf,
                                       ((0, 0, 0), (1, 1)): np.inf}}, 2)
        S, C = sc._build_sc_matrices_all_q(ik, 2, 1, 1, 1)
        self.assertTrue(np.all(S[0, 0, 0] == 0.0),
                        "diagonal CoulombInter must not touch S")

        ik = sc._build_interaction_k(
            k, k, k, {"Ising": {((0, 0, 0), (0, 1)): np.inf,
                                ((0, 0, 0), (1, 0)): np.inf}}, 2)
        S, C = sc._build_sc_matrices_all_q(ik, 2, 1, 1, 1)
        # the density slots (aa, bb) of Ising carry S only
        self.assertEqual(complex(C[0, 0, 0, 0, 3]), 0.0,
                         "density-slot Ising must not touch C")


class TestIEEEParity(unittest.TestCase):

    def test_plus_minus_two_coefficients_preserve_structure(self):
        """The +-2 coefficients: value doubled for finite input, and the
        legacy component structure preserved (numpy promotes the scalar
        multiply to a COMPLEX multiply, matching the pre-refactor
        2.0 * V form; a finite-only check could not distinguish the
        forms)."""
        import warnings

        import hwave.sc as sc

        ctx = warnings.catch_warnings()
        ctx.__enter__()
        self.addCleanup(ctx.__exit__, None, None, None)
        warnings.simplefilter("ignore", RuntimeWarning)

        m = np.zeros((2, 2, 1, 1, 1), dtype=complex)
        m[0, 1] = m[1, 0] = 0.7
        S, C = sc._build_sc_matrices_all_q(
            {"CoulombInter": m}, 2, 1, 1, 1, _presymmetrised=True)
        self.assertAlmostEqual(complex(C[0, 0, 0, 0, 3]), 1.4, places=13)

        # numpy promotes scalar * complex array to a complex multiply, so
        # 2.0 * (Inf+0j) is Inf+NaNj -- IDENTICAL to the pre-refactor
        # 2.0 * V form, and distinguishable from a replace-by-two-adds
        # implementation (which would give Inf+0j)
        m = np.zeros((2, 2, 1, 1, 1), dtype=complex)
        m[0, 1] = m[1, 0] = complex(np.inf, 0.0)
        S, C = sc._build_sc_matrices_all_q(
            {"CoulombInter": m}, 2, 1, 1, 1, _presymmetrised=True)
        got = C[0, 0, 0, 0, 3]          # +2 slot
        self.assertEqual(got.real, np.inf)
        self.assertTrue(np.isnan(got.imag),
                        "the legacy multiply form yields Inf+NaNj here")
        S, C = sc._build_sc_matrices_all_q(
            {"Ising": m}, 2, 1, 1, 1, _presymmetrised=True)
        got = S[0, 0, 0, 0, 3]          # -2 slot
        self.assertEqual(got.real, -np.inf)
        self.assertTrue(np.isnan(got.imag),
                        "the legacy multiply form yields -Inf+NaNj here")

    def test_plus_minus_one_coefficients_preserve_complex_infinities(self):
        """numpy's 1.0 * (Inf+0j) is Inf+NaNj (full complex multiply);
        the +-1 coefficients must use direct add/subtract so non-finite
        values keep the pre-refactor component structure. The interaction
        arrays are built by hand and passed pre-symmetrised: the public
        Fourier builder itself multiplies by complex phases and turns Inf
        into Inf+NaNj before the S/C builder is reached (identically in
        the old and new code), which would mask the accumulator
        semantics under test."""
        import warnings

        import hwave.sc as sc

        ctx = warnings.catch_warnings()
        ctx.__enter__()
        self.addCleanup(ctx.__exit__, None, None, None)
        warnings.simplefilter("ignore", RuntimeWarning)

        def mat(entries):
            m = np.zeros((2, 2, 1, 1, 1), dtype=complex)
            for (a, b), v in entries.items():
                m[a, b, 0, 0, 0] = v
            return m

        inf = complex(np.inf, 0.0)
        cases = [
            ("CoulombIntra", {(0, 0): inf, (1, 1): inf},
             [("S", (0, 0), np.inf), ("C", (0, 0), np.inf)]),
            ("CoulombInter", {(0, 1): inf, (1, 0): inf},
             [("S", (1, 1), np.inf), ("C", (1, 1), -np.inf)]),
            ("PairHop", {(0, 1): inf, (1, 0): inf},
             [("S", (1, 2), np.inf), ("C", (1, 2), np.inf)]),
            # Case 3's +-1 density coefficients (Hund): S at (aa, bb)
            ("Hund", {(0, 1): inf, (1, 0): inf},
             [("S", (0, 3), np.inf), ("C", (0, 3), -np.inf)]),
        ]
        for itype, entries, expected in cases:
            with self.subTest(interaction=itype):
                ik = {itype: mat(entries)}
                S, C = sc._build_sc_matrices_all_q(
                    ik, 2, 1, 1, 1, _presymmetrised=True)
                M = {"S": S[0, 0, 0], "C": C[0, 0, 0]}
                for ch, (i, j), want in expected:
                    got = M[ch][i, j]
                    self.assertEqual(got.real, want,
                                     "%s %s[%d,%d] real" % (itype, ch, i, j))
                    self.assertEqual(got.imag, 0.0,
                                     "%s %s[%d,%d] imag must stay exactly "
                                     "zero, not NaN" % (itype, ch, i, j))


class TestMixedTypeGolden(unittest.TestCase):

    def test_mixed_types_norb3_complex_slots(self):
        """Permanent legacy-parity protection: a mixed-type norb=3 complex
        input, with selected slots of every type and family pinned (the
        pre-refactor and refactored builders were verified bit-for-bit
        equal on this class of input during review; these values freeze
        that behavior)."""
        import hwave.sc as sc

        norb = 3

        def mat(entries):
            m = np.zeros((norb, norb, 1, 1, 1), dtype=complex)
            for (a, b), v in entries.items():
                m[a, b, 0, 0, 0] = v
            return m

        ik = {
            "CoulombIntra": mat({(a, a): 2.0 + 0.1 * a for a in range(3)}),
            "CoulombInter": mat({(0, 1): 0.7, (1, 0): 0.7,
                                 (0, 2): 0.4, (2, 0): 0.4}),
            "Hund": mat({(0, 1): 0.3, (1, 0): 0.3}),
            "Ising": mat({(1, 2): 0.2, (2, 1): 0.2}),
            "Exchange": mat({(0, 2): 0.25, (2, 0): 0.25}),
            "PairHop": mat({(0, 1): 0.15 + 0.05j, (1, 0): 0.15 - 0.05j}),
        }
        S, C = sc._build_sc_matrices_all_q(
            ik, norb, 1, 1, 1, _presymmetrised=True)
        S = S[0, 0, 0]
        C = C[0, 0, 0]
        p = lambda a, b: a * norb + b
        checks = [
            # diag: U_a at (aa, aa); CoulombInter aa density adds nothing here
            ("S", p(0, 0), p(0, 0), 2.0), ("C", p(0, 0), p(0, 0), 2.0),
            ("S", p(1, 1), p(1, 1), 2.1),
            # density (aa, bb): C = 2 V; Hund S = +J, C = -J; Ising S = -2 I
            ("C", p(0, 0), p(1, 1), 2 * 0.7 - 0.3),
            ("S", p(0, 0), p(1, 1), 0.3),
            ("S", p(1, 1), p(2, 2), -2 * 0.2),
            ("C", p(1, 1), p(2, 2), 0.0),
            ("C", p(0, 0), p(2, 2), 2 * 0.4),
            # cross (ab, ab): V (+, -), Hund (-, +), Ising (+, -), Exch (+, +)
            ("S", p(0, 1), p(0, 1), 0.7 - 0.3),
            ("C", p(0, 1), p(0, 1), -0.7 + 0.3),
            ("S", p(1, 2), p(1, 2), 0.2),
            ("C", p(1, 2), p(1, 2), -0.2),
            ("S", p(0, 2), p(0, 2), 0.4 + 0.25),
            ("C", p(0, 2), p(0, 2), -0.4 + 0.25),
            # antidiag (ab, ba): PairHop, complex phase preserved, S = C
            ("S", p(0, 1), p(1, 0), 0.15 + 0.05j),
            ("C", p(0, 1), p(1, 0), 0.15 + 0.05j),
            ("S", p(1, 0), p(0, 1), 0.15 - 0.05j),
        ]
        M = {"S": S, "C": C}
        for ch, i, j, want in checks:
            self.assertAlmostEqual(
                complex(M[ch][i, j]), complex(want), places=12,
                msg="%s[%d,%d]" % (ch, i, j))


if __name__ == "__main__":
    unittest.main()
