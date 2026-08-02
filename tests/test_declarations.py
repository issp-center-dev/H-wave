"""The single-sourced declaration symmetrisation (#108 increment 2).

The module owns the reduction in both representations its consumers use
(real-space dense for the ring assembly, momentum-space for the S/C
builders); these tests pin the algebra once -- most importantly the
FFT-commutation identity that makes the two forms the SAME reduction --
plus bit-parity of each delegation against the pre-refactor inline
expressions.
"""

import os
import unittest

import numpy as np

from hwave.solver.declarations import symmetrise_dense, symmetrise_k
from hwave.solver.kgrid import reverse_fft_axes


def _rand(shape, seed):
    rng = np.random.default_rng(seed)
    return rng.standard_normal(shape) + 1j * rng.standard_normal(shape)


def _fft_ab_q(arr):
    """(nx, ny, nz, a, b) real-space table -> (a, b, nx, ny, nz) with the
    e^{-iqR} phase convention (numpy's forward FFT)."""
    return np.fft.fftn(arr, axes=(0, 1, 2)).transpose(3, 4, 0, 1, 2)


class TestFftCommutation(unittest.TestCase):
    """fft(symmetrise_dense(T)) == symmetrise_k(fft(T)): the identity
    that makes the two representations one reduction. Odd and even
    lattice dimensions, all-complex tables, both rules."""

    def test_same_operator_rule_commutes(self):
        arr = _rand((4, 3, 5, 3, 3), 1)
        lhs = _fft_ab_q(symmetrise_dense(arr))
        rhs = symmetrise_k({"Hund": _fft_ab_q(arr)})["Hund"]
        np.testing.assert_allclose(lhs, rhs, atol=1e-12)

    def test_pairhop_hermitian_rule_commutes(self):
        arr = _rand((4, 3, 5, 3, 3), 2)
        lhs = _fft_ab_q(symmetrise_dense(arr, hermitian=True))
        rhs = symmetrise_k({"PairHop": _fft_ab_q(arr)})["PairHop"]
        np.testing.assert_allclose(lhs, rhs, atol=1e-12)


class TestDenseForm(unittest.TestCase):

    def test_matches_the_pre_refactor_inline_expression(self):
        """Bit-parity: the ring assembly's historical inline expression
        (reverse + orbital transpose (+ conjugate for PairHop), mean)."""
        for herm, seed in ((False, 3), (True, 4)):
            with self.subTest(hermitian=herm):
                arr = _rand((4, 3, 5, 2, 2), seed)
                rev = np.transpose(reverse_fft_axes(arr, (0, 1, 2)),
                                   (0, 1, 2, 4, 3))
                if herm:
                    rev = np.conjugate(rev)
                old = 0.5 * (arr + rev)
                np.testing.assert_array_equal(
                    symmetrise_dense(arr, hermitian=herm), old)

    def test_one_sided_bond_splits_between_both_directions(self):
        """A single +x declaration means (T_ab + T_ba)/2: half the weight
        stays, half lands at the wrapped -x cell on the transposed
        orbital pair -- the reading the ring adjudicated (a one-sided
        off-site V has vertex v cos(qR), not v e^{-iqR})."""
        nx = 4
        arr = np.zeros((nx, 1, 1, 2, 2), dtype=complex)
        arr[1, 0, 0, 0, 1] = 1.0
        sym = symmetrise_dense(arr)
        self.assertAlmostEqual(sym[1, 0, 0, 0, 1].real, 0.5, places=15)
        self.assertAlmostEqual(sym[nx - 1, 0, 0, 1, 0].real, 0.5,
                               places=15)
        self.assertAlmostEqual(float(np.abs(sym).sum()), 1.0, places=15)

    def test_antisymmetric_declaration_vanishes(self):
        arr = np.zeros((1, 1, 1, 2, 2), dtype=complex)
        arr[0, 0, 0, 0, 1] = 1.0
        arr[0, 0, 0, 1, 0] = -1.0
        self.assertEqual(float(np.abs(symmetrise_dense(arr)).max()), 0.0)

    def test_hermitian_closed_pairhop_keeps_its_phase(self):
        p = 0.7 + 0.4j
        arr = np.zeros((1, 1, 1, 2, 2), dtype=complex)
        arr[0, 0, 0, 0, 1] = p
        arr[0, 0, 0, 1, 0] = np.conj(p)
        sym = symmetrise_dense(arr, hermitian=True)
        self.assertAlmostEqual(sym[0, 0, 0, 0, 1], p, places=15)


class TestKForm(unittest.TestCase):

    def test_matches_the_pre_refactor_inline_expression(self):
        """Bit-parity of the moved k-space body: fancy-gather reversal +
        orbital transpose mean; PairHop conjugate-transpose at fixed q."""
        M = _rand((2, 2, 4, 3, 5), 5)
        nx, ny, nz = 4, 3, 5
        Mrev = M[:, :, (-np.arange(nx)) % nx][:, :, :, (-np.arange(ny)) % ny][
            :, :, :, :, (-np.arange(nz)) % nz]
        old_plain = 0.5 * (M + Mrev.transpose(1, 0, 2, 3, 4))
        old_ph = 0.5 * (M + np.conj(M.transpose(1, 0, 2, 3, 4)))
        out = symmetrise_k({"Exchange": M, "PairHop": M})
        np.testing.assert_array_equal(out["Exchange"], old_plain)
        np.testing.assert_array_equal(out["PairHop"], old_ph)

    def test_idempotent(self):
        M = _rand((2, 2, 4, 3, 5), 6)
        once = symmetrise_k({"Hund": M, "PairHop": M})
        twice = symmetrise_k(once)
        for k in once:
            np.testing.assert_allclose(twice[k], once[k], atol=1e-15)


class TestProductionAnchoring(unittest.TestCase):
    """Anchor the convention to PRODUCTION transforms, not to a
    test-local FFT choice (round 1: a sign-flipped helper would pin the
    wrong identity symmetrically). One-direction bond, norb=2, nx=4 so
    q = pi/2 is not self-inverse."""

    def test_dense_route_equals_the_production_k_route(self):
        import hwave.sc as sc
        from hwave.solver.declarations import symmetrise_dense

        nx, norb = 4, 2
        kx = np.linspace(0, 2 * np.pi, nx, endpoint=False)
        kz = np.array([0.0])
        tbl = {"CoulombInter": {((1, 0, 0), (0, 1)): 0.7}}
        prod = sc._symmetrise_interactions_k(
            sc._build_interaction_k(kx, kz, kz, tbl, norb))["CoulombInter"]

        # dense route with _build_interaction_k's layout: an entry
        # (R, (a, b)) lands at [R, b, a] (the vertex pair transpose).
        # R -> q with the documented e^{+iqR} (#133): ifftn * nvol.
        arr = np.zeros((nx, 1, 1, norb, norb), dtype=complex)
        arr[1, 0, 0, 1, 0] = 0.7
        dense = (np.fft.ifftn(symmetrise_dense(arr), axes=(0, 1, 2))
                 * nx).transpose(3, 4, 0, 1, 2)
        np.testing.assert_allclose(dense, prod, atol=1e-13)
        # and the value itself is the adjudicated one-sided reading:
        # half weight on e^{+iq}, half on e^{-iq} at the transposed slot
        q1 = kx[1]
        self.assertAlmostEqual(
            prod[1, 0, 1, 0, 0], 0.35 * np.exp(+1j * q1), places=13)

    def test_rpa_ring_consumes_the_shared_dense_core(self):
        """Delegation spy: the ring assembly must reach
        declarations.symmetrise_dense (the single-sourcing claim, not
        just formula parity)."""
        import tempfile
        from unittest import mock
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.rpa as rpa_mod
        from hwave.solver import declarations

        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("  1.0 0.0 0.0\n  0.0 1.0 0.0\n  0.0 0.0 1.0\n"
                    "1\n 0.0 0.0 0.0\n")
        with open(os.path.join(d, "transfer.dat"), "w") as f:
            f.write("Transfer\n1\n1\n 1\n"
                    "   1    0    0    1    1   1.000000   0.0\n")
        with open(os.path.join(d, "coulombintra.dat"), "w") as f:
            f.write("CoulombIntra\n1\n1\n 1\n"
                    "   0    0    0    1    1   2.000000   0.0\n")
        read_io = read_input_k.QLMSkInput(
            {"interaction": {"path_to_input": d, "Geometry": "geom.dat",
                             "Transfer": "transfer.dat",
                             "CoulombIntra": "coulombintra.dat"}})
        info_mode = {'mode': 'RPA',
                     'param': {'T': 1.0, 'filling': 0.5,
                               'CellShape': [4, 1, 1],
                               'SubShape': [1, 1, 1], 'Nmat': 4}}
        with mock.patch.object(
                rpa_mod, "symmetrise_dense",
                side_effect=declarations.symmetrise_dense) as spy:
            rpa_mod.RPA(read_io.get_param("ham"), {}, info_mode)
        self.assertGreater(spy.call_count, 0)

    def test_sc_builders_consume_the_shared_k_form(self):
        """Delegation spy for both S/C builder routes."""
        from unittest import mock
        import hwave.sc as sc
        from hwave.solver import declarations

        M = _rand((2, 2, 2, 2, 1), 9)
        with mock.patch.object(
                sc, "symmetrise_k",
                side_effect=declarations.symmetrise_k) as spy:
            sc._build_sc_matrices_all_q({"Hund": M}, 2, 2, 2, 1)
            self.assertEqual(spy.call_count, 1)
            sc._build_sc_matrices({"Hund": M}, 2, 0, 0, 0)
        self.assertEqual(spy.call_count, 2)


class TestSimplePathSymmetrisation(unittest.TestCase):
    """The one-orbital simple/load vertex path read RAW tables (round 3:
    a one-sided off-site V entered as v e^{-iqR} while every other
    consumer read v cos(qR); measured drift 0.7 at q = pi/2). Pin the
    fix through the production _compute_vertices dispatch."""

    def test_one_sided_bond_reads_as_cos_q(self):
        import hwave.sc as sc

        nx, norb, nmat = 4, 1, 4
        kx = np.linspace(0, 2 * np.pi, nx, endpoint=False)
        kz = np.array([0.0])
        inter_k = sc._build_interaction_k(
            kx, kz, kz, {"CoulombInter": {((1, 0, 0), (0, 0)): 0.7}}, norb)
        chi0q = np.zeros((norb, norb, nx, 1, 1, nmat), dtype=complex)
        V = sc._compute_vertices(chi0q, inter_k, norb, nx, 1, 1, nmat)
        # with chi0 = 0 the vertex is the bare (Wc + Ws)/2 = V(q); the
        # symmetrised one-sided declaration is 0.7 cos(q): ~0 at pi/2,
        # NOT the raw -0.7i
        v_pi2 = complex(np.asarray(V)[..., 1, 0, 0].flat[0])
        self.assertLess(abs(v_pi2), 1e-12)
        v_0 = complex(np.asarray(V)[..., 0, 0, 0].flat[0])
        self.assertAlmostEqual(v_0.real, 0.7, places=12)

    def test_unknown_interaction_type_raises(self):
        with self.assertRaises(ValueError) as cm:
            symmetrise_k({"PairHpo": np.zeros((1, 1, 2, 1, 1),
                                              dtype=complex)})
        self.assertIn("PairHpo", str(cm.exception))

    def test_dense_form_is_idempotent(self):
        for herm, seed in ((False, 12), (True, 13)):
            with self.subTest(hermitian=herm):
                arr = _rand((4, 3, 5, 2, 2), seed)
                once = symmetrise_dense(arr, hermitian=herm)
                twice = symmetrise_dense(once, hermitian=herm)
                np.testing.assert_allclose(twice, once, atol=1e-15)

    def test_ring_and_sc_agree_on_an_asymmetric_declaration(self):
        """The drift class this module exists to prevent, committed as a
        measured comparison: an asymmetric on-site declaration must
        reduce to the SAME mean coefficient in the ring's dense route
        and the S/C k-route."""
        import hwave.sc as sc
        from hwave.solver.declarations import symmetrise_dense

        norb = 2
        k1 = np.array([0.0])
        tbl = {"Exchange": {((0, 0, 0), (0, 1)): 1.0,
                            ((0, 0, 0), (1, 0)): 0.35}}
        mean = (1.0 + 0.35) / 2
        # S/C k-route
        S, C = sc._build_sc_matrices_all_q(
            sc._build_interaction_k(k1, k1, k1, tbl, norb), norb, 1, 1, 1)
        self.assertAlmostEqual(S[0, 0, 0, 1, 1].real, mean, places=12)
        # ring dense route
        arr = np.zeros((1, 1, 1, norb, norb), dtype=complex)
        arr[0, 0, 0, 0, 1] = 1.0
        arr[0, 0, 0, 1, 0] = 0.35
        sym = symmetrise_dense(arr)
        self.assertAlmostEqual(sym[0, 0, 0, 0, 1].real, mean, places=12)
        self.assertAlmostEqual(sym[0, 0, 0, 1, 0].real, mean, places=12)


class TestClosureGate(unittest.TestCase):
    """Round 5: re-averaging an ALREADY-closed configuration mixed
    independently rounded +R/-R exponentials and changed saved artifacts
    at the 1e-18 level. The simple path now decides on the RAW tables
    (exact: 0.5*(x+x) == x in IEEE) and skips the k-space averaging for
    closed input, restoring byte parity."""

    def test_closure_detection(self):
        import hwave.sc as sc

        closed = {"CoulombInter": {((1, 0, 0), (0, 1)): 0.7,
                                   ((-1, 0, 0), (1, 0)): 0.7}}
        one_sided = {"CoulombInter": {((1, 0, 0), (0, 1)): 0.7}}
        herm = {"CoulombInter": {((0, 0, 0), (0, 1)): 0.5 + 0.2j,
                                 ((0, 0, 0), (1, 0)): 0.5 - 0.2j}}
        self.assertTrue(
            sc._declarations_partner_closed(closed, 4, 1, 1, 2))
        self.assertFalse(
            sc._declarations_partner_closed(one_sided, 4, 1, 1, 2))
        # Hermitian-closed complex is NOT same-operator closed: the mean
        # folds it to Re(p), so it must be (re)symmetrised
        self.assertFalse(
            sc._declarations_partner_closed(herm, 4, 1, 1, 2))
        self.assertTrue(sc._declarations_partner_closed({}, 4, 1, 1, 2))

    def test_huge_closed_coefficient_detected_without_overflow(self):
        """Round 6: the mean-based comparison overflowed for closed
        coefficients above float64.max / 2, reading them as open (with
        RuntimeWarnings). The direct partner comparison forms no mean."""
        import warnings
        import hwave.sc as sc

        big = np.finfo(np.float64).max * 0.75
        closed = {"CoulombInter": {((1, 0, 0), (0, 1)): big,
                                   ((-1, 0, 0), (1, 0)): big}}
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            self.assertTrue(
                sc._declarations_partner_closed(closed, 4, 1, 1, 2))
        nonfinite = {"CoulombInter": {((1, 0, 0), (0, 1)): np.inf,
                                      ((-1, 0, 0), (1, 0)): np.inf}}
        self.assertFalse(
            sc._declarations_partner_closed(nonfinite, 4, 1, 1, 2))

    def test_e2e_symmetric_run_takes_the_raw_pipeline(self):
        """End-to-end capture through calc_eliashberg (round 6): on a
        symmetric fixture the closure detector must return True in the
        production flow and the simple path must never call the k-space
        symmetriser -- i.e. the executed pipeline IS the pre-fix one, so
        saved artifacts are byte-identical by construction."""
        import tempfile
        from unittest import mock
        import hwave.sc as sc

        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("  1.0 0.0 0.0\n  0.0 1.0 0.0\n  0.0 0.0 1.0\n"
                    "1\n 0.0 0.0 0.0\n")
        with open(os.path.join(d, "transfer.dat"), "w") as f:
            f.write("Transfer\n1\n2\n 1 1\n"
                    "   1    0    0    1    1   1.000000   0.0\n"
                    "  -1    0    0    1    1   1.000000   0.0\n")
        # CoulombIntra only: the auto tensor selector (PR #128) chooses
        # the general scheme whenever CoulombInter is present, so the
        # SIMPLE path is reachable from the calc route only for U-only
        # input (V reaches it via chi0q_mode='load' with a 2-index file).
        # On-site U is partner-closed by construction.
        with open(os.path.join(d, "coulombintra.dat"), "w") as f:
            f.write("CoulombIntra\n1\n1\n 1\n"
                    "   0    0    0    1    1   2.000000   0.0\n")
        out = os.path.join(d, "out")
        os.makedirs(out, exist_ok=True)
        inp = {
            "mode": {"param": {"T": 2.0, "filling": 0.5,
                               "CellShape": [4, 4, 1],
                               "SubShape": [1, 1, 1], "Nmat": 8}},
            "file": {"input": {"path_to_input": d,
                               "interaction": {"path_to_input": d,
                                               "Geometry": "geom.dat",
                                               "Transfer": "transfer.dat",
                                               "CoulombIntra":
                                                   "coulombintra.dat"}},
                     "output": {"path_to_output": out}},
            "eliashberg": {"chi0q_mode": "calc", "frequency": "static",
                           "pairing_type": "singlet", "init_gap": "cos",
                           "solver_mode": "eigenvalue",
                           "eigenvalue_method": "arnoldi",
                           "num_eigenvalues": 1,
                           "output_eigenvalue": "eig.dat",
                           "output_gap": "gap.dat"},
        }
        real_closed = sc._declarations_partner_closed
        real_symm = sc._symmetrise_interactions_k
        real_simple = sc._compute_vertices_simple
        with mock.patch.object(sc, "_declarations_partner_closed",
                               side_effect=real_closed) as det_spy, \
                mock.patch.object(sc, "_symmetrise_interactions_k",
                                  side_effect=real_symm) as sym_spy, \
                mock.patch.object(sc, "_compute_vertices_simple",
                                  side_effect=real_simple) as simple_spy:
            sc.calc_eliashberg(inp)
        # the run really took the simple path, the detector fired and
        # returned True, and the k-space symmetriser was never invoked:
        # the executed pipeline is byte-identical to the pre-fix one
        self.assertEqual(simple_spy.call_count, 1)
        self.assertEqual(det_spy.call_count, 1)
        self.assertEqual(sym_spy.call_count, 0)
        self.assertTrue(os.path.exists(os.path.join(out, "gap.dat")))

    def test_closed_input_keeps_the_raw_bits(self):
        """With declarations_closed=True the vertex is computed from the
        UNTOUCHED k-tables: bitwise equal to the pre-fix raw formula."""
        import hwave.sc as sc

        nx, norb, nmat = 4, 1, 4
        kx = np.linspace(0, 2 * np.pi, nx, endpoint=False)
        kz = np.array([0.0])
        tbl = {"CoulombInter": {((1, 0, 0), (0, 0)): 0.7,
                                ((-1, 0, 0), (0, 0)): 0.7}}
        inter_k = sc._build_interaction_k(kx, kz, kz, tbl, norb)
        chi0q = np.zeros((norb, norb, nx, 1, 1, nmat), dtype=complex)
        V = sc._compute_vertices(chi0q, inter_k, norb, nx, 1, 1, nmat,
                                 declarations_closed=True)
        # pre-fix raw formula at chi0 = 0: (Wc + Ws)/2 = V_k, untouched
        raw = np.asarray(inter_k["CoulombInter"]).transpose(2, 3, 4, 0, 1)
        got = np.asarray(V[0] + V[1]) if isinstance(V, tuple) else np.asarray(V)
        flat_raw = raw.reshape(got.shape) if got.shape != raw.shape else raw
        np.testing.assert_array_equal(got.real, flat_raw.real)
        np.testing.assert_array_equal(got.imag, flat_raw.imag)


class TestIEEESpecialParity(unittest.TestCase):
    """Committed parity for non-finite inputs (verified adversarially in
    review; pinned here so a refactor cannot change the propagation):
    both forms reproduce the pre-refactor expressions bitwise, NaN and
    Inf included."""

    def test_dense_with_nonfinite_entries(self):
        arr = _rand((2, 2, 2, 2, 2), 10)
        arr[0, 0, 0, 0, 1] = np.nan
        arr[1, 0, 1, 1, 0] = np.inf + 1j
        for herm in (False, True):
            with self.subTest(hermitian=herm):
                rev = np.transpose(reverse_fft_axes(arr, (0, 1, 2)),
                                   (0, 1, 2, 4, 3))
                if herm:
                    rev = np.conjugate(rev)
                old = 0.5 * (arr + rev)
                new = symmetrise_dense(arr, hermitian=herm)
                np.testing.assert_array_equal(new.real, old.real)
                np.testing.assert_array_equal(new.imag, old.imag)

    def test_k_with_nonfinite_entries(self):
        M = _rand((2, 2, 2, 2, 1), 11)
        M[0, 1, 0, 0, 0] = np.nan
        M[1, 0, 1, 1, 0] = -np.inf
        nx, ny, nz = 2, 2, 1
        Mrev = M[:, :, (-np.arange(nx)) % nx][:, :, :, (-np.arange(ny)) % ny][
            :, :, :, :, (-np.arange(nz)) % nz]
        old_plain = 0.5 * (M + Mrev.transpose(1, 0, 2, 3, 4))
        old_ph = 0.5 * (M + np.conj(M.transpose(1, 0, 2, 3, 4)))
        out = symmetrise_k({"Ising": M, "PairHop": M})
        for got, want in ((out["Ising"], old_plain),
                          (out["PairHop"], old_ph)):
            np.testing.assert_array_equal(got.real, want.real)
            np.testing.assert_array_equal(got.imag, want.imag)




class TestHermitianClosureValidation(unittest.TestCase):
    """Read-time rejection of unclosed declaration tables (issue #93)."""

    def test_disagreeing_pair_raises_naming_everything(self):
        from hwave.solver.declarations import validate_hermitian_closure

        tbl = {((0, 0, 0), (0, 1)): 0.3, ((0, 0, 0), (1, 0)): 0.9}
        with self.assertRaises(ValueError) as cm:
            validate_hermitian_closure("Hund", tbl, source="hund.dat")
        msg = str(cm.exception)
        for frag in ("hund.dat", "0.3", "0.9", "Hermitian-closed"):
            self.assertIn(frag, msg)

    def test_one_sided_declaration_raises(self):
        from hwave.solver.declarations import validate_hermitian_closure

        with self.assertRaises(ValueError):
            validate_hermitian_closure(
                "CoulombInter", {((1, 0, 0), (0, 1)): 0.7})

    def test_closed_tables_pass(self):
        from hwave.solver.declarations import validate_hermitian_closure

        validate_hermitian_closure(
            "CoulombInter", {((1, 0, 0), (0, 1)): 0.7,
                             ((-1, 0, 0), (1, 0)): 0.7})
        p = 0.7 + 0.4j
        validate_hermitian_closure(
            "PairHop", {((0, 0, 0), (0, 1)): p,
                        ((0, 0, 0), (1, 0)): np.conj(p)})
        validate_hermitian_closure("CoulombIntra",
                                   {((0, 0, 0), (0, 0)): 4.0})

    def test_coulombintra_rules(self):
        from hwave.solver.declarations import validate_hermitian_closure

        with self.assertRaises(ValueError):
            validate_hermitian_closure(
                "CoulombIntra", {((1, 0, 0), (0, 0)): 1.0})
        with self.assertRaises(ValueError):
            validate_hermitian_closure(
                "CoulombIntra", {((0, 0, 0), (0, 1)): 1.0})


if __name__ == "__main__":
    unittest.main()
