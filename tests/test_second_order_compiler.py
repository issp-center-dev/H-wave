"""The on-site interaction compiler (spec 2026-09-08 section 2.2): the
normative per-type records as literal fixtures, the first-derivative
contract against the Hartree-Fock kernel, the density-slot contract against
the adjudicated S/C table, Hermiticity, and the D7 refusals."""
import itertools
import unittest

import numpy as np

_TYPES = ("CoulombIntra", "CoulombInter", "Hund", "Ising", "Exchange", "PairHop", "PairLift")


def _tbl(itype, rows):
    """{type: {((0,0,0),(a,b)): v}} -- one type, on-site rows."""
    return {itype: {((0, 0, 0), (a, b)): complex(v) for (a, b, v) in rows}}


def _gi(a, s, norb):
    return s * norb + a


def _expected_gamma(itype, norb, rows):
    """Hand-written normal-ordered records of spec 2.2, antisymmetrised.
    H = 1/2 sum V_{pq,rs} c+_p c+_q c_s c_r ; a monomial x c+_p c+_q c_s c_r
    adds x to V[p,q,r,s] and V[q,p,s,r]."""
    M = 2 * norb
    V = np.zeros((M, M, M, M), complex)

    def add(p, q, s_, r, x):
        V[p, q, r, s_] += x
        V[q, p, s_, r] += x
    up, dn = 0, 1
    for (a, b, v) in rows:
        v = complex(v)
        # The reversal closure represents a bond by BOTH ordered entries
        # (a, b) and (b, a); for the types whose transposed row repeats the
        # SAME operator, each ordered entry therefore carries half the bond.
        # This is the Hartree-Fock kernel's convention, pinned by
        # test_first_derivative_equals_accumulate_hf. PairHop is excluded:
        # its (b, a) row is the reverse pair hop, a different operator.
        if itype in ("CoulombInter", "Hund", "Ising", "Exchange", "PairLift"):
            v = 0.5 * v
        g = lambda o, sp: _gi(o, sp, norb)
        if itype == "CoulombIntra":
            add(g(a, up), g(a, dn), g(a, dn), g(a, up), v)
        elif itype == "CoulombInter":
            for s1, s2 in itertools.product((up, dn), repeat=2):
                add(g(a, s1), g(b, s2), g(b, s2), g(a, s1), v)
        elif itype == "Hund":
            for s1 in (up, dn):
                add(g(a, s1), g(b, s1), g(b, s1), g(a, s1), -v)
        elif itype == "Ising":
            for s1, s2 in itertools.product((up, dn), repeat=2):
                add(g(a, s1), g(b, s2), g(b, s2), g(a, s1), v if s1 == s2 else -v)
        elif itype == "Exchange":
            for s1 in (up, dn):
                ms = 1 - s1
                add(g(a, s1), g(b, ms), g(b, s1), g(a, ms), -v)
        elif itype == "PairHop":
            add(g(a, up), g(a, dn), g(b, up), g(b, dn), -v)
        elif itype == "PairLift":
            add(g(a, up), g(b, up), g(a, dn), g(b, dn), -v)
            add(g(b, dn), g(a, dn), g(b, up), g(a, up), -v)
    return V - V.transpose(0, 1, 3, 2)


class TestLiteralFixtures(unittest.TestCase):

    def test_hubbard_cells(self):
        from hwave.solver.second_order import compile_onsite
        G = compile_onsite(_tbl("CoulombIntra", [(0, 0, 1.0)]), 1)
        self.assertAlmostEqual(G[0, 1, 0, 1], 1.0)      # (a up)(a dn),(a up)(a dn)
        self.assertAlmostEqual(G[0, 1, 1, 0], -1.0)     # swapped placement
        self.assertAlmostEqual(G[1, 0, 1, 0], 1.0)
        self.assertEqual(np.count_nonzero(G), 4)

    def test_every_type_matches_the_hand_written_records(self):
        from hwave.solver.second_order import compile_onsite
        norb = 2
        for itype in _TYPES:
            rows = [(0, 0, 1.0)] if itype == "CoulombIntra" else [(0, 1, 1.0), (1, 0, 1.0)]
            with self.subTest(itype=itype):
                G = compile_onsite(_tbl(itype, rows), norb)
                np.testing.assert_allclose(G, _expected_gamma(itype, norb, rows), atol=1e-14)
                self.assertGreater(np.abs(G).max(), 0.0)

    def test_complex_rows(self):
        """Same-operator types fold a complex coefficient to its real part
        over the closed orbit; PairHop keeps its phase (mirrored row v*)."""
        from hwave.solver.second_order import compile_onsite
        norb, v = 2, 1.0 + 0.5j
        for itype in ("Exchange", "PairLift"):
            with self.subTest(itype=itype):
                G = compile_onsite(_tbl(itype, [(0, 1, v), (1, 0, np.conj(v))]), norb)
                # closure: (v + conj(conj(v)))/2 = v per row; the two
                # half-weighted ordered rows sum to Re v on the orbit
                exp = _expected_gamma(itype, norb, [(0, 1, v), (1, 0, np.conj(v))])
                np.testing.assert_allclose(G, exp, atol=1e-14)
                np.testing.assert_allclose(G, G.conj().transpose(2, 3, 0, 1), atol=1e-14)
        G = compile_onsite(_tbl("PairHop", [(0, 1, v), (1, 0, np.conj(v))]), norb)
        exp = _expected_gamma("PairHop", norb, [(0, 1, v), (1, 0, np.conj(v))])
        np.testing.assert_allclose(G, exp, atol=1e-14)
        self.assertGreater(np.abs(G.imag).max(), 0.1)              # phase preserved
        np.testing.assert_allclose(G, G.conj().transpose(2, 3, 0, 1), atol=1e-14)


class TestContracts(unittest.TestCase):

    def _random_rho(self, norb, seed):
        """Generic Hermitian spin-block-diagonal density with orbital coherences.
        This is the sector the solver actually runs in."""
        rng = np.random.default_rng(seed)
        M = 2 * norb
        rho = np.zeros((M, M), complex)
        for s in (0, 1):
            X = rng.normal(size=(norb, norb)) + 1j * rng.normal(size=(norb, norb))
            blk = 0.5 * (X + X.conj().T) * 0.1 + 0.4 * np.eye(norb)
            rho[s * norb:(s + 1) * norb, s * norb:(s + 1) * norb] = blk
        return rho

    def _general_rho(self, norb, seed):
        """Generic Hermitian density with spin OFF-DIAGONAL blocks. Without
        these the kernel's spin-flip channels are identically zero and the
        Exchange (Hartree) and PairLift subtests compare 0 against 0."""
        rng = np.random.default_rng(seed)
        M = 2 * norb
        X = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M))
        return 0.1 * (X + X.conj().T) + 0.4 * np.eye(M)

    def test_first_derivative_equals_accumulate_hf(self):
        from hwave.solver import hartree_fock as hf
        from hwave.solver.second_order import compile_onsite, compile_onsite_v, hf_first_order
        norb = 2
        shape = (1, 1, 1)
        for itype in _TYPES:
            rows = [(0, 0, 0.7), (1, 1, 0.3)] if itype == "CoulombIntra" else [(0, 1, 0.7), (1, 0, 0.7)]
            param_ham = {itype: {((0, 0, 0), (a, b)): v for (a, b, v) in rows}}
            G = compile_onsite(_tbl(itype, rows), norb)
            V = compile_onsite_v(_tbl(itype, rows), norb)
            tabs = hf.build_interaction_tables(param_ham, norb, shape)
            for dens, rho in (("block", self._random_rho(norb, 3)),
                              ("general", self._general_rho(norb, 5))):
                # kernel density: rho_so[r, s, a, t, b] = <c+_{sa}(0) c_{tb}(r)>
                rho_so = rho.reshape(2, norb, 2, norb)[None]
                for fock in (True, False):
                    with self.subTest(itype=itype, dens=dens, fock=fock):
                        out = np.zeros((1, 2 * norb, 2 * norb), complex)
                        hf.accumulate_hf(out, rho_so, tabs.inter_table, tabs.spin_table, shape,
                                         include_fock=fock)
                        if dens == "general":
                            # anti-vacuity: a general density leaves no channel
                            # of any type in the kernel's zero sector
                            self.assertGreater(np.abs(out).max(), 0.0)
                        # Hartree + Fock from Gamma; the Hartree part alone from
                        # the non-antisymmetrised V (same contraction)
                        sig = hf_first_order(G if fock else V, rho)
                        np.testing.assert_allclose(sig, out[0], atol=1e-14,
                                                   err_msg="if this fails only by complex conjugation, "
                                                           "the kernel's density is <c+_{tb} c_{sa}>: "
                                                           "transpose rho in THIS test, not Gamma")

    def test_density_slots_reproduce_the_adjudicated_sc_table(self):
        from hwave.solver.second_order import compile_onsite, density_slots
        from hwave.solver.vertex_table import ADJUDICATED_SC
        norb = 2
        for itype, fam in (("CoulombIntra", "diag"), ("CoulombInter", "density"),
                           ("Hund", "density"), ("Ising", "density")):
            rows = [(0, 0, 1.0)] if itype == "CoulombIntra" else [(0, 1, 1.0), (1, 0, 1.0)]
            G = compile_onsite(_tbl(itype, rows), norb)
            d = density_slots(G, norb)
            a, b = (0, 0) if itype == "CoulombIntra" else (0, 1)
            S = d[0, 1, a, b] - d[0, 0, a, b]
            C = d[0, 0, a, b] + d[0, 1, a, b]
            s_ref, c_ref = ADJUDICATED_SC[itype][fam]
            with self.subTest(itype=itype):
                self.assertAlmostEqual(S, s_ref, places=12)
                self.assertAlmostEqual(C, c_ref, places=12)

    def test_hermiticity_guard_and_d7(self):
        from hwave.solver.second_order import compile_onsite, DegenerateRowError
        norb = 2
        with self.assertRaises(ValueError) as cm:
            compile_onsite(_tbl("Exchange", [(0, 1, 1.0 + 0.5j), (1, 0, 1.0 + 0.5j)]), norb,
                           closed=True)   # bypass the closure: not Hermitian
        self.assertIn("Hermitian", str(cm.exception))
        for itype, hint in (("CoulombInter", "CoulombIntra"), ("Hund", "one-body"),
                            ("Ising", "CoulombIntra"), ("Exchange", "CoulombIntra"),
                            ("PairHop", "2 Re"), ("PairLift", "identically zero")):
            with self.subTest(itype=itype):
                with self.assertRaises(DegenerateRowError) as cm:
                    compile_onsite(_tbl(itype, [(0, 0, 0.5)]), norb)
                self.assertIn(itype, str(cm.exception))
                self.assertIn(hint, str(cm.exception))



class TestRowGuards(unittest.TestCase):
    """Refusals that name the offending ROW.

    These are the diagnostics for a caller that reaches the compiler
    without the reader -- a programmatic pipeline, a fixture, a future
    scheme. Every one of them is a route to a silently wrong ``Gamma``
    rather than a crash, which is why they are refusals and not comments."""

    def test_non_finite_coefficients_are_refused_before_the_hermiticity_guard(self):
        """A NaN or Inf coupling would otherwise pass BOTH guards: the
        Hermiticity test is ``dev > herm_tol * scale``, and a comparison
        against NaN is FALSE, so a NaN table compiled without complaint and
        the non-finite values reached the kernel."""
        from hwave.solver.second_order import compile_onsite, compile_onsite_v
        for bad in (float("nan"), float("inf"), -float("inf"), complex(1.0, float("nan"))):
            for itype, rows in (("CoulombIntra", [(0, 0, bad)]),
                                ("CoulombInter", [(0, 1, bad), (1, 0, bad)])):
                with self.subTest(value=bad, itype=itype):
                    for fn in (compile_onsite, compile_onsite_v):
                        with self.assertRaises(ValueError) as cm:
                            fn(_tbl(itype, rows), 2)
                        msg = str(cm.exception)
                        self.assertIn("non-finite", msg)
                        self.assertIn(itype, msg)
                        self.assertIn("orbitals", msg)

    def test_a_nan_table_would_otherwise_survive_the_tensor_guard(self):
        """The reason the check above exists, measured: with the row check
        bypassed (``closed=True`` still runs it, so this reproduces the
        comparison directly) a NaN deviation does not exceed any tolerance."""
        dev = float("nan")
        self.assertFalse(dev > 1e-12 * 1.0)

    def test_orbital_indices_outside_the_model_are_refused(self):
        from hwave.solver.second_order import compile_onsite
        for rows in ([(0, 5, 1.0), (5, 0, 1.0)], [(-1, 0, 1.0), (0, -1, 1.0)]):
            with self.subTest(rows=rows):
                with self.assertRaises(ValueError) as cm:
                    compile_onsite(_tbl("CoulombInter", rows), 2)
                msg = str(cm.exception)
                self.assertIn("CoulombInter", msg)
                self.assertIn("outside", msg)

    def test_a_nonpositive_norb_is_refused(self):
        from hwave.solver.second_order import compile_onsite
        for norb in (0, -1):
            with self.subTest(norb=norb):
                with self.assertRaises(ValueError) as cm:
                    compile_onsite(_tbl("CoulombIntra", [(0, 0, 1.0)]), norb)
                self.assertIn("norb", str(cm.exception))

    def test_a_row_that_is_not_hermitian_closed_is_named(self):
        """The closure SYMMETRISES whatever it is given, so a table whose
        transposed row is not the conjugate is not rejected by it -- it is
        quietly replaced by its Hermitian part. The per-row guard names the
        type and the row; the tensor-wide guard, which stays, cannot."""
        from hwave.solver.second_order import compile_onsite
        v = 1.0 + 0.5j
        for itype in _TYPES:
            rows = [(0, 0, v)] if itype == "CoulombIntra" else [(0, 1, v), (1, 0, v)]
            with self.subTest(itype=itype):
                with self.assertRaises(ValueError) as cm:
                    compile_onsite(_tbl(itype, rows), 2)
                msg = str(cm.exception)
                self.assertIn("Hermitian-closed", msg)
                self.assertIn(itype, msg)
                self.assertIn("orbitals", msg)

    def test_a_lone_row_is_accepted_and_closed(self):
        """Only a DECLARED pair of ordered rows is checked: a lone row is
        legal input and the closure supplies its conjugate partner, which is
        what ``close_onsite_rows`` documents. Without this exemption the
        guard would refuse a table the solver has always accepted."""
        from hwave.solver.second_order import compile_onsite, close_onsite_rows
        v = 1.0 + 0.5j
        tbl = _tbl("CoulombInter", [(0, 1, v)])
        closed = close_onsite_rows(tbl)["CoulombInter"]
        self.assertAlmostEqual(closed[(0, 1)], 0.5 * v)
        self.assertAlmostEqual(closed[(1, 0)], 0.5 * np.conj(v))
        G = compile_onsite(tbl, 2)
        self.assertGreater(np.abs(G).max(), 0.0)                   # anti-vacuity
        np.testing.assert_allclose(G, G.conj().transpose(2, 3, 0, 1), atol=1e-14)


if __name__ == "__main__":
    unittest.main()
