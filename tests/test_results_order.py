"""Tests of ``hwave.results.order``: the order parameters of a one-cell
density, their exact Fourier definitions and the sum-rule tolerance.

The input ``rho[o, s, t]`` is the stored Green slot ``<c^dag_{o s} c_{o t}>``
(G_ab = <c^dag_a c_b>, as both UHF solvers write it). For a spin along +y
that slot is ``[[1/2, i/2], [-i/2, 1/2]]``; the transposed (operator density
matrix) reading would flip the sign of S^y, which ``test_sy_sign_convention``
and the comparisons with the solver's own totals guard against.
"""
import cmath
import math
import unittest
from fractions import Fraction

import numpy as np

from hwave.results import (
    OrderParameters, cosine_amplitude, open_run, order_parameters_from_density,
    read_energy, read_scf_status, sum_rule_tolerance)
from tests.results_fixtures import assert_deep_frozen, fixture

SIGMA = (np.array([[0, 1], [1, 0]], dtype=complex),
         np.array([[0, -1j], [1j, 0]]),
         np.array([[1, 0], [0, -1]], dtype=complex))


def slot(n, s):
    """The Green slot ``G_st = <c^dag_s c_t>`` of one orbital with charge
    ``n`` and spin ``s``: the transpose of ``n/2 + s.sigma``."""
    rho_op = 0.5 * n * np.eye(2) + sum(c * m for c, m in zip(s, SIGMA))
    return rho_op.T.copy()


def density(charges, spins):
    return np.array([slot(n, s) for n, s in zip(charges, spins)])


def q_index(result, q):
    return result.q_frac.index(tuple(Fraction(x) for x in q))


def square_cell():
    """The four sites (x, y) of a 2x2 cell, site index x + 2y."""
    return np.array([[x, y, 0] for y in (0, 1) for x in (0, 1)], dtype=float)


def random_rotation(rng):
    a = rng.normal(size=(3, 3))
    q, r = np.linalg.qr(a)
    q = q * np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


class TestDefinitions(unittest.TestCase):

    def test_sy_sign_convention(self):
        """``rho`` is the Green slot <c^dag_s c_t>: <c^dag_up c_down> = i/2
        is a spin along +y (the transposed reading gives -y)."""
        rho = np.array([[[0.5, 0.5j], [-0.5j, 0.5]]])
        r = order_parameters_from_density(
            rho, site_positions=[[0, 0, 0]], cell_shape=(1, 1, 1),
            lattice_shape=(1, 1, 1))
        np.testing.assert_allclose(r.spin[0], [0, 0.5, 0], rtol=0, atol=1e-15)
        np.testing.assert_allclose(r.charge, [1.0], rtol=0, atol=1e-15)
        # and the helper of this module builds the same slot
        np.testing.assert_array_equal(slot(1.0, (0, 0.5, 0)), rho[0])

    def test_site_values(self):
        spins = [(0.1, -0.2, 0.3), (-0.05, 0.25, 0.0)]
        r = order_parameters_from_density(
            density([0.7, 1.4], spins), site_positions=[[0, 0, 0], [1, 0, 0]],
            cell_shape=(2, 1, 1), lattice_shape=(2, 1, 1))
        self.assertIsInstance(r, OrderParameters)
        np.testing.assert_allclose(r.spin, spins, rtol=0, atol=1e-15)
        np.testing.assert_allclose(r.charge, [0.7, 1.4], rtol=0, atol=1e-15)
        np.testing.assert_allclose(
            r.moment, np.linalg.norm(spins, axis=1), rtol=0, atol=1e-15)
        np.testing.assert_allclose(r.uniform, np.mean(spins, axis=0),
                                   rtol=0, atol=1e-15)
        np.testing.assert_allclose(r.spin_fourier[q_index(r, (0, 0, 0))],
                                   r.uniform, rtol=0, atol=1e-15)
        self.assertEqual(r.charge_amplitude[q_index(r, (0, 0, 0))], 0.0)
        self.assertEqual(r.units, {
            "spin": "hbar (same convention as Sz in energy.dat)",
            "charge": "electrons"})
        self.assertIsNone(r.source)
        self.assertTrue(r.transverse_available)
        self.assertEqual(r.n_cell, 1)

    def test_formula_fractional_positions(self):
        pos = [(0.0, 0.0, 0.0), (1 / 3, 0.0, 0.0), (0.25, 0.5, 0.0)]
        spins = [(0.1, 0.2, -0.3), (-0.4, 0.05, 0.15), (0.2, -0.1, 0.35)]
        charges = [1.2, 0.9, 0.6]
        r = order_parameters_from_density(
            density(charges, spins), site_positions=pos, cell_shape=(3, 2, 1),
            lattice_shape=(6, 2, 1))
        nbar = sum(charges) / 3
        for q in ((Fraction(1, 3), 0, 0), (Fraction(-1, 3), Fraction(1, 2), 0)):
            want_m = [0j, 0j, 0j]
            want_n = 0j
            for rs, ss, ns in zip(pos, spins, charges):
                phase = cmath.exp(2j * math.pi * sum(float(qa) * ra
                                                     for qa, ra in zip(q, rs)))
                for a in range(3):
                    want_m[a] += ss[a] * phase / 3
                want_n += (ns - nbar) * phase / 3
            i = q_index(r, q)
            np.testing.assert_allclose(r.spin_fourier[i], want_m, rtol=0, atol=1e-14)
            np.testing.assert_allclose(r.charge_fourier[i], want_n, rtol=0, atol=1e-14)
            self.assertAlmostEqual(r.spin_amplitude[i],
                                   math.sqrt(sum(abs(m) ** 2 for m in want_m)),
                                   delta=1e-14)
            self.assertAlmostEqual(r.charge_amplitude[i], abs(want_n), delta=1e-14)
            # the tau = 0.25 site makes the amplitudes non-real
            self.assertGreater(np.max(np.abs(np.imag(r.spin_fourier[i]))), 1e-3)
            self.assertGreater(abs(np.imag(r.charge_fourier[i])), 1e-3)

    def test_q_set(self):
        r = order_parameters_from_density(
            density([1] * 12, [(0, 0, 0)] * 12),
            site_positions=[[x, y, 0] for y in range(3) for x in range(4)],
            cell_shape=(4, 3, 1), lattice_shape=(4, 3, 1))
        F = Fraction
        qx = [F(-1, 4), F(0), F(1, 4), F(1, 2)]
        qy = [F(-1, 3), F(0), F(1, 3)]
        self.assertEqual(r.q_frac, tuple((a, b, F(0)) for a in qx for b in qy))
        for q in r.q_frac:
            for c in q:
                self.assertIsInstance(c, Fraction)
        self.assertEqual(r.spin_fourier.shape, (12, 3))
        self.assertEqual(r.charge_fourier.shape, (12,))

    def test_neel_any_axis(self):
        rng = np.random.default_rng(359)
        u = rng.normal(size=3)
        u /= np.linalg.norm(u)
        eps = [(-1) ** (x + y) for y in (0, 1) for x in (0, 1)]
        r = order_parameters_from_density(
            density([1] * 4, [0.5 * e * u for e in eps]),
            site_positions=square_cell(), cell_shape=(2, 2, 1),
            lattice_shape=(4, 4, 1))
        self.assertAlmostEqual(r.staggered_norm, 0.5, delta=1e-15)
        self.assertAlmostEqual(r.uniform_norm, 0.0, delta=1e-15)
        np.testing.assert_allclose(r.staggered.imag, 0, rtol=0, atol=1e-15)
        np.testing.assert_allclose(r.staggered.real, 0.5 * u, rtol=0, atol=1e-15)
        self.assertEqual(r.spin_ordering_vector,
                         (Fraction(1, 2), Fraction(1, 2), Fraction(0)))
        self.assertAlmostEqual(r.mean_local_moment, 0.5, delta=1e-15)

    def test_su2_rotation_invariance(self):
        """A global SU(2) rotation of the Green slots rotates every spin by
        the same SO(3) matrix and leaves every norm unchanged."""
        rng = np.random.default_rng(20261122)
        spins = rng.normal(size=(4, 3)) * 0.2
        charges = rng.uniform(0.5, 1.5, size=4)
        rho = density(charges, spins)
        pos = square_cell() + [0.1, 0.3, 0]
        kw = dict(site_positions=pos, cell_shape=(2, 2, 1), lattice_shape=(2, 2, 1))
        base = order_parameters_from_density(rho, **kw)
        # U = exp(-i theta n.sigma / 2); rho_op -> U rho_op U^dag, so the
        # slot G = rho_op^T -> conj(U) G U^T
        n = rng.normal(size=3)
        n /= np.linalg.norm(n)
        theta = 1.234
        U = (math.cos(theta / 2) * np.eye(2)
             - 1j * math.sin(theta / 2) * sum(c * m for c, m in zip(n, SIGMA)))
        rot = order_parameters_from_density(
            np.array([U.conj() @ g @ U.T for g in rho]), **kw)
        R = np.array([[0.5 * np.trace(SIGMA[a] @ U @ SIGMA[b] @ U.conj().T).real
                       for b in range(3)] for a in range(3)])
        np.testing.assert_allclose(rot.spin, base.spin @ R.T, rtol=0, atol=1e-14)
        self.assertGreater(np.max(np.abs(rot.spin - base.spin)), 1e-2)
        for name in ("moment", "spin_amplitude", "charge", "charge_amplitude"):
            np.testing.assert_allclose(getattr(rot, name), getattr(base, name),
                                       rtol=0, atol=1e-14, err_msg=name)
        for name in ("uniform_norm", "staggered_norm", "mean_local_moment",
                     "charge_order_amplitude"):
            self.assertAlmostEqual(getattr(rot, name), getattr(base, name),
                                   delta=1e-14, msg=name)
        self.assertEqual(rot.spin_ordering_vector, base.spin_ordering_vector)
        self.assertEqual(rot.charge_ordering_vector, base.charge_ordering_vector)

    def test_ferromagnet(self):
        u = np.array([0.3, -0.4, 0.5])
        u /= np.linalg.norm(u)
        r = order_parameters_from_density(
            density([1] * 4, [0.5 * u] * 4), site_positions=square_cell(),
            cell_shape=(2, 2, 1), lattice_shape=(2, 2, 1))
        self.assertAlmostEqual(r.uniform_norm, 0.5, delta=1e-15)
        self.assertAlmostEqual(r.staggered_norm, 0.0, delta=1e-15)
        self.assertEqual(r.spin_ordering_vector, (0, 0, 0))
        self.assertIsNone(r.charge_ordering_vector)
        self.assertAlmostEqual(r.charge_order_amplitude, 0.0, delta=1e-15)

    def test_charge_checkerboard(self):
        charges = [2 if (x + y) % 2 == 0 else 0 for y in (0, 1) for x in (0, 1)]
        r = order_parameters_from_density(
            density(charges, [(0, 0, 0)] * 4), site_positions=square_cell(),
            cell_shape=(2, 2, 1), lattice_shape=(2, 2, 1))
        self.assertAlmostEqual(r.charge_order_amplitude, 1.0, delta=1e-15)
        self.assertEqual(r.charge_ordering_vector,
                         (Fraction(1, 2), Fraction(1, 2), Fraction(0)))
        self.assertIsNone(r.spin_ordering_vector)
        np.testing.assert_allclose(r.charge, charges, rtol=0, atol=0)

    def test_spiral_120(self):
        spins = [0.5 * np.array([math.cos(2 * math.pi * j / 3),
                                 math.sin(2 * math.pi * j / 3), 0]) for j in range(3)]
        r = order_parameters_from_density(
            density([1] * 3, spins), site_positions=[[j, 0, 0] for j in range(3)],
            cell_shape=(3, 1, 1), lattice_shape=(3, 1, 1))
        for q in ((Fraction(1, 3), 0, 0), (Fraction(-1, 3), 0, 0)):
            self.assertAlmostEqual(r.spin_amplitude[q_index(r, q)],
                                   math.sqrt(2) / 4, delta=1e-15)
        self.assertAlmostEqual(r.spin_amplitude[q_index(r, (0, 0, 0))], 0,
                               delta=1e-15)
        self.assertAlmostEqual(float(np.sum(r.spin_amplitude ** 2)), 0.25, delta=1e-15)
        self.assertEqual(r.spin_ordering_vector, (Fraction(-1, 3), 0, 0))
        # local moment = sqrt(2) |M(Q)| for a single-Q coplanar spiral
        np.testing.assert_allclose(
            r.moment, math.sqrt(2) * r.spin_amplitude[q_index(r, (Fraction(1, 3), 0, 0))],
            rtol=0, atol=1e-15)
        # L = (3, 1, 1): B_x = 3 is odd, so there is no staggered moment
        self.assertIsNone(r.staggered)
        self.assertIsNone(r.staggered_norm)

    def test_ordering_tie_tolerance(self):
        """Amplitudes within 1e-12 tie and go to the smallest |q|^2, then to
        the lexicographically smallest q; a larger gap decides by size."""
        def result(c):
            d = 0.4
            spins = [(d, 0, c), (0, 0, -c), (-d, 0, c), (0, 0, -c)]
            return order_parameters_from_density(
                density([1] * 4, spins), site_positions=[[j, 0, 0] for j in range(4)],
                cell_shape=(4, 1, 1), lattice_shape=(4, 1, 1))
        r = result(0.2 + 5e-13)
        i_half, i_quarter = q_index(r, (0.5, 0, 0)), q_index(r, (0.25, 0, 0))
        self.assertGreater(r.spin_amplitude[i_half], r.spin_amplitude[i_quarter])
        self.assertEqual(r.spin_ordering_vector, (Fraction(-1, 4), 0, 0))
        self.assertEqual(result(0.2 + 1e-9).spin_ordering_vector,
                         (Fraction(1, 2), 0, 0))

    def test_ordering_floor(self):
        """No ordering vector when the largest amplitude is <= 1e-10."""
        def result(m):
            return order_parameters_from_density(
                density([1 + m, 1 - m], [(0, 0, m), (0, 0, -m)]),
                site_positions=[[0, 0, 0], [1, 0, 0]], cell_shape=(2, 1, 1),
                lattice_shape=(2, 1, 1))
        small, large = result(5e-11), result(1e-9)
        self.assertIsNone(small.spin_ordering_vector)
        self.assertIsNone(small.charge_ordering_vector)
        self.assertAlmostEqual(small.charge_order_amplitude, 5e-11, delta=1e-16)
        self.assertEqual(large.spin_ordering_vector, (Fraction(1, 2), 0, 0))
        self.assertEqual(large.charge_ordering_vector, (Fraction(1, 2), 0, 0))

    def test_ordering_tie_prefers_small_q_over_lexicographic(self):
        """A tie between (+-1/4, 1/2, 0) and (0, 1/2, 0): the smaller |q|^2
        wins although (-1/4, 1/2, 0) is lexicographically smaller."""
        c, d = 0.2 + 5e-13, 0.4
        pos, spins = [], []
        for y in (0, 1):
            for x in range(4):
                pos.append([x, y, 0])
                spins.append(((-1) ** y * d * [1, 0, -1, 0][x], 0, (-1) ** y * c))
        r = order_parameters_from_density(
            density([1] * 8, spins), site_positions=pos, cell_shape=(4, 2, 1),
            lattice_shape=(4, 2, 1))
        amp = r.spin_amplitude
        self.assertAlmostEqual(amp[q_index(r, (-0.25, 0.5, 0))], 0.2, delta=1e-15)
        self.assertAlmostEqual(amp[q_index(r, (0, 0.5, 0))], 0.2, delta=1e-12)
        self.assertEqual(r.spin_ordering_vector, (0, Fraction(1, 2), 0))

    def test_stripe_not_staggered(self):
        r = order_parameters_from_density(
            density([1, 1], [(0, 0, 0.5), (0, 0, -0.5)]),
            site_positions=[[0, 0, 0], [1, 0, 0]], cell_shape=(2, 1, 1),
            lattice_shape=(4, 4, 1))
        self.assertIsNone(r.staggered)
        self.assertIsNone(r.staggered_norm)
        self.assertEqual(r.spin_ordering_vector, (Fraction(1, 2), 0, 0))

    def test_staggered_only_on_extended_axes(self):
        """q_AF has 1/2 only where L > 1: a chain along x in a (2, 1, 1)
        lattice is staggered at (1/2, 0, 0)."""
        r = order_parameters_from_density(
            density([1, 1], [(0, 0, 0.5), (0, 0, -0.5)]),
            site_positions=[[0, 0, 0], [1, 0, 0]], cell_shape=(2, 1, 1),
            lattice_shape=(4, 1, 1))
        np.testing.assert_allclose(r.staggered, [0, 0, 0.5], rtol=0, atol=1e-15)
        self.assertAlmostEqual(r.staggered_norm, 0.5, delta=1e-15)

    def test_cosine_precondition(self):
        spins = np.array([(0.1, 0, 0.4), (0.2, -0.3, -0.1)])
        q = (Fraction(1, 2), 0, 0)
        for tau, equal in ((0.0, True), (0.25, False)):
            pos = np.array([[0, 0, 0], [1, 0, 0]], dtype=float) + [tau, 0, 0]
            r = order_parameters_from_density(
                density([1, 1], spins), site_positions=pos, cell_shape=(2, 1, 1),
                lattice_shape=(2, 1, 1))
            cos = cosine_amplitude(spins, pos, q)
            want = np.mean([math.cos(2 * math.pi * 0.5 * p[0]) * s
                            for p, s in zip(pos, spins)], axis=0)
            np.testing.assert_allclose(cos, want, rtol=0, atol=1e-15)
            same = np.allclose(cos, r.spin_fourier[q_index(r, q)], rtol=0, atol=1e-12)
            self.assertEqual(same, equal, msg="tau = {}".format(tau))
        # a scalar per site works as well
        self.assertAlmostEqual(
            float(cosine_amplitude([2.0, 0.0], [[0, 0, 0], [1, 0, 0]], q)), 1.0,
            delta=1e-15)

    def test_cosine_rejects_complex_s(self):
        pos = [[0, 0, 0], [1, 0, 0]]
        q = (Fraction(1, 2), 0, 0)
        with self.assertRaisesRegex(ValueError, "S must be real"):
            cosine_amplitude([0.1 + 1e-3j, 0.2], pos, q)
        # an imaginary part within the stated tolerance is dropped
        got = cosine_amplitude(np.array([0.1 + 1e-14j, 0.2]), pos, q)
        self.assertAlmostEqual(float(got), (0.1 - 0.2) / 2, delta=1e-15)
        self.assertFalse(np.iscomplexobj(got))
        with self.assertRaisesRegex(ValueError, "S must be"):
            cosine_amplitude(["a", "b"], pos, q)

    def test_cosine_rejects_bad_q(self):
        pos = [[0, 0, 0], [1, 0, 0]]
        for q in (0.5, (0.5, 0), (0.5, 0, 0, 0), ("x", 0, 0), (float("nan"), 0, 0),
                  (1j, 0, 0), None):
            with self.subTest(q=q):
                with self.assertRaisesRegex(ValueError, "q must"):
                    cosine_amplitude([1.0, 0.0], pos, q)

    def test_site_of_orbital_required(self):
        rho = density([1, 1], [(0, 0, 0.5), (0, 0, -0.5)])
        for kw in (dict(site_positions=[[0, 0, 0]], cell_shape=(1, 1, 1),
                        lattice_shape=(1, 1, 1)),
                   dict(cell_shape=(1, 1, 1), lattice_shape=(1, 1, 1))):
            with self.assertRaisesRegex(
                    ValueError, "site_of_orbital is required when a cell holds "
                                "several orbitals"):
                order_parameters_from_density(rho, **kw)
        r = order_parameters_from_density(
            rho, site_of_orbital=[0, 0], site_positions=[[0, 0, 0]],
            cell_shape=(1, 1, 1), lattice_shape=(1, 1, 1))
        self.assertEqual(r.spin.shape, (1, 3))

    def test_site_of_orbital_default(self):
        rho = density([1, 1], [(0, 0, 0.5), (0, 0, -0.5)])
        r = order_parameters_from_density(
            rho, site_positions=[[0, 0, 0], [1, 0, 0]], cell_shape=(2, 1, 1),
            lattice_shape=(2, 1, 1))
        np.testing.assert_array_equal(r.site_of_orbital, [0, 1])
        # from the cell shape when there are no positions
        r = order_parameters_from_density(rho, cell_shape=(2, 1, 1),
                                          lattice_shape=(2, 1, 1))
        np.testing.assert_array_equal(r.site_of_orbital, [0, 1])

    def test_site_sum(self):
        r = order_parameters_from_density(
            density([1, 1], [(0.1, 0.2, 0.5), (-0.1, -0.2, -0.5)]),
            site_of_orbital=[0, 0], site_positions=[[0, 0, 0]],
            cell_shape=(1, 1, 1), lattice_shape=(2, 2, 1))
        np.testing.assert_allclose(r.spin, [[0, 0, 0]], rtol=0, atol=1e-15)
        np.testing.assert_allclose(r.charge, [2.0], rtol=0, atol=1e-15)
        self.assertAlmostEqual(r.mean_local_moment, 0.0, delta=1e-15)
        np.testing.assert_allclose(r.moment, [0.0], rtol=0, atol=1e-15)

    def test_mean_local_moment_inside_cell(self):
        """An antiferromagnet inside one cell (two sites at the same cell
        position, no lattice Q captures it): only mean_local_moment sees it."""
        r = order_parameters_from_density(
            density([1, 1], [(0, 0, 0.5), (0, 0, -0.5)]),
            site_of_orbital=[0, 1], site_positions=[[0.25, 0, 0], [0.75, 0, 0]],
            cell_shape=(1, 1, 1), lattice_shape=(4, 1, 1))
        self.assertAlmostEqual(r.uniform_norm, 0.0, delta=1e-15)
        self.assertIsNone(r.spin_ordering_vector)
        self.assertAlmostEqual(r.mean_local_moment, 0.5, delta=1e-15)

    def test_positions_none(self):
        r = order_parameters_from_density(
            density([1, 1], [(0, 0, 0.5), (0, 0, 0.1)]), cell_shape=(2, 1, 1),
            lattice_shape=(2, 1, 1))
        self._assert_no_fourier(r)
        np.testing.assert_allclose(r.uniform, [0, 0, 0.3], rtol=0, atol=1e-15)
        self.assertEqual(r.cell_shape, (2, 1, 1))

    def test_without_lattice_shape(self):
        """No lattice: the per-site and uniform (Q = 0) values only."""
        r = order_parameters_from_density(
            density([1, 0.5], [(0.2, 0, 0.3), (0, -0.1, 0.1)]),
            site_of_orbital=[0, 1], site_positions=[[0, 0, 0], [1, 0, 0]])
        self._assert_no_fourier(r)
        self.assertIsNone(r.cell_shape)
        self.assertIsNone(r.lattice_shape)
        np.testing.assert_allclose(r.uniform, [0.1, -0.05, 0.2], rtol=0, atol=1e-15)
        self.assertAlmostEqual(r.uniform_norm, math.sqrt(0.01 + 0.0025 + 0.04),
                               delta=1e-15)
        np.testing.assert_allclose(r.charge, [1, 0.5], rtol=0, atol=1e-15)
        self.assertGreater(r.mean_local_moment, 0.2)

    def _assert_no_fourier(self, r):
        for name in ("q_frac", "spin_fourier", "spin_amplitude",
                     "spin_ordering_vector", "charge_fourier", "charge_amplitude",
                     "charge_ordering_vector", "charge_order_amplitude",
                     "staggered", "staggered_norm"):
            self.assertIsNone(getattr(r, name), msg=name)

    def test_transverse_flag_and_n_cell(self):
        r = order_parameters_from_density(
            density([1], [(0, 0, 0.5)]), cell_shape=(1, 1, 1),
            lattice_shape=(4, 4, 1), n_cell=16, transverse_available=False)
        self.assertFalse(r.transverse_available)
        self.assertEqual(r.n_cell, 16)
        self.assertEqual(r.lattice_shape, (4, 4, 1))
        # M is per site; N_site * n_cell * M is the total
        np.testing.assert_allclose(r.uniform, [0, 0, 0.5], rtol=0, atol=1e-15)

    def test_n_cell_against_shapes(self):
        """With both shapes, n_cell is prod(L) / prod(B): derived when not
        given, and a different value is refused (e.g. the number of lattice
        sites instead of the number of magnetic cells)."""
        rho = density([1] * 4, [(0, 0, 0.5)] * 4)
        kw = dict(site_positions=square_cell(), cell_shape=(2, 2, 1),
                  lattice_shape=(4, 4, 1))
        self.assertEqual(order_parameters_from_density(rho, **kw).n_cell, 4)
        self.assertEqual(order_parameters_from_density(rho, n_cell=4, **kw).n_cell, 4)
        with self.assertRaisesRegex(ValueError, "n_cell 16 .*4 magnetic cells"):
            order_parameters_from_density(rho, n_cell=16, **kw)
        # without both shapes there is nothing to check; the default is 1
        r = order_parameters_from_density(rho, site_positions=square_cell(),
                                          cell_shape=(2, 2, 1))
        self.assertEqual(r.n_cell, 1)
        r = order_parameters_from_density(rho, site_of_orbital=[0, 1, 2, 3],
                                          n_cell=7)
        self.assertEqual(r.n_cell, 7)

    def test_input_validation(self):
        good = density([1, 1], [(0, 0, 0.5), (0, 0, -0.5)])
        pos = [[0, 0, 0], [1, 0, 0]]
        shapes = dict(cell_shape=(2, 1, 1), lattice_shape=(2, 1, 1))
        bad_herm = good.copy()
        bad_herm[0, 0, 1] += 1e-9
        bad_nan = good.copy()
        bad_nan[1, 0, 0] = np.nan
        cases = [
            ("rho", dict(rho=good[:, 0])),
            ("rho", dict(rho=np.zeros((2, 3, 3)))),
            ("rho", dict(rho=np.zeros((0, 2, 2)))),
            ("rho", dict(rho=bad_herm)),
            ("rho", dict(rho=bad_nan)),
            ("rho", dict(rho=np.array([["a", "b"], ["c", "d"]]))),
            ("site_of_orbital", dict(site_of_orbital=[0])),
            ("site_of_orbital", dict(site_of_orbital=[0, 2],
                                     site_positions=[[0, 0, 0]] * 3)),
            ("site_of_orbital", dict(site_of_orbital=[-1, 0])),
            ("site_of_orbital", dict(site_of_orbital=[0, 0.5])),
            ("site_of_orbital", dict(site_of_orbital=[[0, 1]])),
            ("site_positions", dict(site_positions=[[0, 0], [1, 0]])),
            ("site_positions", dict(site_positions=[[0, 0, 0]] * 3)),
            ("site_positions", dict(site_positions=[[0, 0, 0], [np.inf, 0, 0]])),
            ("cell_shape", dict(cell_shape=(0, 1, 1))),
            ("cell_shape", dict(cell_shape=(2, 1))),
            ("cell_shape", dict(cell_shape=(2.5, 1, 1))),
            ("cell_shape", dict(cell_shape=(True, 1, 1))),
            ("lattice_shape", dict(lattice_shape=(3, 1, 1))),
            ("lattice_shape", dict(lattice_shape=(-2, 1, 1))),
            ("n_cell", dict(n_cell=0)),
            ("n_cell", dict(n_cell=1.5)),
            ("n_cell", dict(n_cell=True)),
        ]
        for name, change in cases:
            kw = dict(rho=good, site_positions=pos, **shapes)
            kw.update(change)
            rho = kw.pop("rho")
            with self.subTest(name=name, change=sorted(change)):
                with self.assertRaisesRegex(ValueError, name):
                    order_parameters_from_density(rho, **kw)
        # integer-valued arrays are accepted
        r = order_parameters_from_density(
            good, site_of_orbital=np.array([0.0, 1.0]), site_positions=pos,
            cell_shape=np.array([2, 1, 1]), lattice_shape=[6, 1, 1],
            n_cell=np.int64(3))
        self.assertEqual(r.cell_shape, (2, 1, 1))
        self.assertEqual(r.n_cell, 3)
        self.assertIsInstance(r.n_cell, int)

    def test_arrays_copied_and_frozen(self):
        rho = density([1, 1, 0.5], [(0, 0, 0.5), (0.1, 0, -0.5), (0, 0.2, 0)])
        pos = np.array([[0, 0, 0], [1, 0, 0], [0.5, 0.5, 0]])
        sites = np.array([0, 1, 2])
        r = order_parameters_from_density(
            rho, site_of_orbital=sites, site_positions=pos, cell_shape=(2, 2, 1),
            lattice_shape=(2, 2, 1))
        assert_deep_frozen(r)
        before = (r.spin.copy(), r.site_positions.copy())

        def arrays(obj):
            if isinstance(obj, np.ndarray):
                yield obj
            elif isinstance(obj, dict):
                for v in obj.values():
                    yield from arrays(v)
            elif isinstance(obj, (list, tuple)):
                for v in obj:
                    yield from arrays(v)

        import dataclasses
        found = 0
        for f in dataclasses.fields(r):
            for arr in arrays(getattr(r, f.name)):
                found += 1
                for inp in (rho, pos, sites):
                    self.assertFalse(np.shares_memory(arr, inp), msg=f.name)
        self.assertGreater(found, 10)
        rho[:] = 0
        pos[:] = 7
        sites[:] = 0
        np.testing.assert_array_equal(r.spin, before[0])
        np.testing.assert_array_equal(r.site_positions, before[1])
        np.testing.assert_array_equal(r.site_of_orbital, [0, 1, 2])
        with self.assertRaises(dataclasses.FrozenInstanceError):
            r.uniform_norm = 1.0


class TestSumRuleTolerance(unittest.TestCase):

    def test_tolerance_formula(self):
        got = sum_rule_tolerance(1e-7, mix=0.5, g_size=1024, nvol=16, nd=8)
        want = 0.5 * 16 * math.sqrt(8) * 1e-7 * 1024 / 2 + 16 * 8 * 1e-12 + 1e-12
        self.assertAlmostEqual(got, want, delta=1e-20)
        # the mixing term vanishes at mix = 1 and doubles from 0.5 to 0
        floor = sum_rule_tolerance(1e-7, mix=1.0, g_size=1024, nvol=16, nd=8)
        self.assertAlmostEqual(floor, 16 * 8 * 1e-12 + 1e-12, delta=1e-24)
        zero = sum_rule_tolerance(1e-7, mix=0.0, g_size=1024, nvol=16, nd=8)
        self.assertAlmostEqual(zero - floor, 2 * (got - floor), delta=1e-18)
        # UHFr: g_size = 4 N^2, nvol = 1, nd = 2 N
        n = 16
        got = sum_rule_tolerance(3e-9, mix=0.3, g_size=4 * n * n, nvol=1, nd=2 * n)
        want = 0.7 * math.sqrt(2 * n) * 3e-9 * 4 * n * n / 2 + 2 * n * 1e-12 + 1e-12
        self.assertAlmostEqual(got, want, delta=1e-22)

    def test_tolerance_rejects_bad_residual(self):
        for residual in (float("nan"), float("inf"), -1e-3, None):
            with self.subTest(residual=residual):
                with self.assertRaisesRegex(ValueError, "residual"):
                    sum_rule_tolerance(residual, mix=0.5, g_size=4, nvol=1, nd=2)


def _stub_calc_phys(green, enable_spin_orbital, nd, norb, ns):
    """The solver's own totals (``UHFk._calc_phys``) on a one-cell stub."""
    import hwave.solver.uhfk as uhfk_module
    stub = object.__new__(uhfk_module.UHFk)
    stub.shape = (1, 1, 1)
    stub.nvol, stub.nd, stub.norb, stub.ns = 1, nd, norb, ns
    stub.enable_spin_orbital = enable_spin_orbital
    stub.Green = np.asarray(green, dtype=complex)
    stub.Green_prev = stub.Green.copy()
    stub.param_mod = {"Mix": 1.0}
    stub.threshold = 1e-12
    stub.physics = {}
    stub._calc_phys()
    return stub.physics


class TestAgainstSolverTotals(unittest.TestCase):
    """N_site * n_cell * M and n_cell * sum(n_s) against the totals that the
    solver writes to energy.dat, on a random transverse state in both index
    conventions, and on a real run along a generic axis."""

    def _check(self, rho, phys):
        r = order_parameters_from_density(rho, cell_shape=(len(rho), 1, 1),
                                          lattice_shape=(len(rho), 1, 1))
        n_site = len(r.charge)
        np.testing.assert_allclose(n_site * r.uniform,
                                   [phys["Sx"], phys["Sy"], phys["Sz"]],
                                   rtol=0, atol=1e-13)
        self.assertAlmostEqual(float(np.sum(r.charge)), phys["NCond"], delta=1e-13)
        self.assertGreater(abs(phys["Sy"]), 1e-2)
        self.assertGreater(abs(phys["Sx"]), 1e-2)

    def _random_spin_block(self, norb, seed):
        rng = np.random.default_rng(seed)
        a = rng.normal(size=(2 * norb, 2 * norb)) + 1j * rng.normal(size=(2 * norb, 2 * norb))
        return (a + a.conj().T) / 2          # index s * norb + orb

    def test_normal_mode(self):
        norb = 3
        mat = self._random_spin_block(norb, 359)
        green = mat.reshape(1, 2, norb, 2, norb)
        phys = _stub_calc_phys(green, False, nd=2 * norb, norb=norb, ns=2)
        rho = np.array([green[0][:, a, :, a] for a in range(norb)])
        self._check(rho, phys)

    def test_spin_orbital_mode(self):
        from hwave.solver import spin_orbital
        norb = 3
        mat = self._random_spin_block(norb, 1122)
        g0 = spin_orbital.spinblock_to_interleaved(mat, norb)
        phys = _stub_calc_phys(g0.reshape(1, 2 * norb, 2 * norb), True,
                               nd=2 * norb, norb=2 * norb, ns=1)
        rho = np.array([g0[2 * j:2 * j + 2, 2 * j:2 * j + 2] for j in range(norb)])
        self._check(rho, phys)

    def test_ferromagnet_generic_axis_run(self):
        """The saved Green function of a UHFk ferromagnet along a generic
        axis u (fixture ``uhfk_fm_axis``) against Sx, Sy, Sz and NCond of
        its energy.dat, within ``sum_rule_tolerance``."""
        toml = fixture("uhfk_fm_axis")
        run = open_run(toml)
        with np.load(run.path("green")) as npz:
            green = npz["green"]
        nvol = green.shape[0]
        self.assertEqual(green.shape, (16, 2, 1, 2, 1))
        rho = green[0][:, 0, :, 0][None]
        r = order_parameters_from_density(
            rho, site_positions=[[0, 0, 0]], cell_shape=(1, 1, 1),
            lattice_shape=(4, 4, 1), n_cell=nvol)
        energy = read_energy(run)
        status = read_scf_status(run)
        self.assertEqual(status.status, "converged")
        tol = sum_rule_tolerance(status.residual, mix=0.5, g_size=green.size,
                                 nvol=nvol, nd=2)
        self.assertLess(tol, 1e-6)
        total = len(r.charge) * r.n_cell * r.uniform
        np.testing.assert_allclose(total, [energy.sx, energy.sy, energy.sz],
                                   rtol=0, atol=tol)
        self.assertAlmostEqual(r.n_cell * float(np.sum(r.charge)), 5.0, delta=tol)
        u = np.array([0.6, 0.2537, 0.7734])
        u /= np.linalg.norm(u)
        np.testing.assert_allclose(r.uniform / r.uniform_norm, u, rtol=0, atol=1e-8)
        self.assertGreater(energy.sy, 0.5)
        self.assertEqual(r.spin_ordering_vector, (0, 0, 0))


if __name__ == "__main__":
    unittest.main()
