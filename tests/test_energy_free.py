#!/usr/bin/env python3
"""Finite-temperature energy output: internal energy E and free energy F.

At finite temperature the mean-field code minimizes the Helmholtz free energy
F = E - T*S.  The solver must expose BOTH the internal energy E (= sum_n eps_n
f(eps_n) - <V>/2) and the free energy F.  At T=0 the two coincide (S=0), so the
historical "Energy_Total" semantics are preserved.

These tests drive _calc_energy with a stubbed solver state (no interactions, so
the total energy is just the band energy) and check the band-energy
decomposition analytically.
"""

import unittest

import numpy as np

from hwave.solver.uhfk import UHFk


def _make_uhfk_stub(T, eigenvalues, mu):
    """Build a minimal UHFk object exercising only the band-energy path."""
    solver = UHFk.__new__(UHFk)

    nl = len(eigenvalues)
    solver.shape = (1, 1, 1)
    solver.nvol = 1
    solver.norb = nl
    solver.ns = 1
    solver.nd = nl
    solver.T = T
    solver.ene_cutoff = 100.0
    solver.enable_spin_orbital = False
    solver.iflag_fock = True

    ws = np.array(eigenvalues, dtype=float).reshape(1, nl)
    vs = np.eye(nl, dtype=complex).reshape(1, nl, nl)
    solver._green_list = {
        "eigenvalue": [ws],
        "eigenvector": [vs],
        "mu": [mu],
    }
    solver.block_to_group = [0]
    solver.group_nconds = [nl // 2 if nl > 1 else 1]

    solver.Green = np.zeros((1, 1, nl, 1, nl), dtype=complex)
    solver.inter_table = {}
    solver.spin_table = {}
    solver.physics = {}
    return solver


def _fermi(T, mu, eps):
    return 1.0 / (1.0 + np.exp((np.asarray(eps) - mu) / T))


class TestFiniteTempEnergyAndFreeEnergy(unittest.TestCase):
    def test_both_energy_and_free_energy_are_output(self):
        T = 0.5
        eigenvalues = [-1.0, 1.0]
        mu = 0.0
        solver = _make_uhfk_stub(T, eigenvalues, mu)

        solver._calc_energy()

        # Internal energy and free energy must both be present.
        self.assertIn("Ene", solver.physics)
        self.assertIn("FreeEne", solver.physics)
        self.assertIn("Total", solver.physics["Ene"])
        self.assertIn("Total", solver.physics["FreeEne"])

    def test_internal_band_energy_is_sum_eps_f(self):
        T = 0.5
        eigenvalues = np.array([-1.0, 1.0])
        mu = 0.0
        solver = _make_uhfk_stub(T, eigenvalues, mu)

        solver._calc_energy()

        f = _fermi(T, mu, eigenvalues)
        e_band_internal = np.sum(eigenvalues * f)  # -0.761594...
        self.assertAlmostEqual(
            solver.physics["Ene"]["Band"].real, e_band_internal, places=10
        )

    def test_free_band_energy_is_grand_potential(self):
        T = 0.5
        eigenvalues = np.array([-1.0, 1.0])
        mu = 0.0
        solver = _make_uhfk_stub(T, eigenvalues, mu)

        solver._calc_energy()

        f = _fermi(T, mu, eigenvalues)
        nn = np.sum(f)
        omega0 = -T * np.sum(np.log1p(np.exp(-(eigenvalues - mu) / T)))
        e_band_free = mu * nn + omega0  # -1.126928...
        self.assertAlmostEqual(
            solver.physics["FreeEne"]["Band"].real, e_band_free, places=10
        )

    def test_free_energy_not_above_internal_energy(self):
        # F = E - T*S with S >= 0, so F <= E at finite temperature.
        T = 0.5
        eigenvalues = np.array([-1.0, 1.0])
        mu = 0.0
        solver = _make_uhfk_stub(T, eigenvalues, mu)

        solver._calc_energy()

        E = solver.physics["Ene"]["Total"].real
        F = solver.physics["FreeEne"]["Total"].real
        self.assertLessEqual(F, E + 1.0e-12)

    def test_zero_temperature_energy_equals_free_energy(self):
        # At T=0, S=0 so E == F and Energy_Total semantics are preserved.
        T = 0.0
        eigenvalues = np.array([-1.0, 1.0])
        mu = 0.0
        solver = _make_uhfk_stub(T, eigenvalues, mu)

        solver._calc_energy()

        E = solver.physics["Ene"]["Total"].real
        F = solver.physics["FreeEne"]["Total"].real
        self.assertAlmostEqual(E, F, places=12)


class TestUHFrFiniteTempEnergyAndFreeEnergy(unittest.TestCase):
    """UHFr (real-space) must expose the same internal/free energy split."""

    def _make_uhfr_stub(self, T, eigenvalues, mu, occupied):
        from hwave.solver.uhfr import UHFr

        solver = UHFr.__new__(UHFr)
        nl = len(eigenvalues)
        solver.Nsize = nl
        solver.T = T
        solver.ene_cutoff = 100.0

        ev = np.array(eigenvalues, dtype=float)
        vec = np.eye(nl, dtype=complex)
        solver.green_list = {
            "all": {
                "eigenvalue": ev,
                "eigenvector": vec,
                "mu": mu,
                "occupied": occupied,
            }
        }
        # No interaction: Ham_local = 0 so InterAll vanishes.
        solver.Green = np.zeros((2 * nl, 2 * nl), dtype=complex)
        solver.Ham_local = np.zeros(((2 * nl) ** 2, (2 * nl) ** 2), dtype=complex)
        solver.physics = {}
        return solver

    def test_internal_and_free_band_energy(self):
        T = 0.5
        eigenvalues = np.array([-1.0, 1.0])
        mu = 0.0
        solver = self._make_uhfr_stub(T, eigenvalues, mu, occupied=1)

        solver._calc_energy()

        f = _fermi(T, mu, eigenvalues)
        e_band_internal = np.sum(eigenvalues * f)
        nn = np.sum(f)
        omega0 = -T * np.sum(np.log1p(np.exp(-(eigenvalues - mu) / T)))
        e_band_free = mu * nn + omega0

        self.assertAlmostEqual(
            solver.physics["Ene"]["band"].real, e_band_internal, places=10
        )
        self.assertIn("FreeEne", solver.physics)
        self.assertAlmostEqual(
            solver.physics["FreeEne"]["band"].real, e_band_free, places=10
        )
        self.assertLessEqual(
            solver.physics["FreeEne"]["Total"].real,
            solver.physics["Ene"]["Total"].real + 1.0e-12,
        )

    def test_zero_temperature_equivalence(self):
        T = 0.0
        eigenvalues = np.array([-1.0, 1.0])
        mu = 0.0
        solver = self._make_uhfr_stub(T, eigenvalues, mu, occupied=1)

        solver._calc_energy()

        self.assertAlmostEqual(
            solver.physics["Ene"]["Total"].real,
            solver.physics["FreeEne"]["Total"].real,
            places=12,
        )


class TestChemicalPotentialZeroT(unittest.TestCase):
    """At T=0 a real chemical potential (Fermi level) is computed and output.

    The Fermi level is taken as the midpoint of the HOMO-LUMO gap, which is the
    T -> 0+ limit of the finite-temperature chemical-potential equation.
    """

    def test_uhfk_zero_t_fermi_level_is_gap_midpoint(self):
        solver = UHFk.__new__(UHFk)
        solver.shape = (1, 1, 1)
        solver.nvol = 1
        ws = [np.array([[-3.0, -1.0, 2.0, 5.0]])]
        vs = [np.eye(4, dtype=complex).reshape(1, 4, 4)]

        _dists, mu = solver._find_dist_group_zero_t(ws, vs, 2)

        self.assertAlmostEqual(mu, 0.5, places=12)  # (-1 + 2) / 2

    def test_uhfk_zero_t_fermi_level_edge_cases(self):
        solver = UHFk.__new__(UHFk)
        solver.shape = (1, 1, 1)
        solver.nvol = 1
        ws = [np.array([[-3.0, -1.0, 2.0, 5.0]])]
        vs = [np.eye(4, dtype=complex).reshape(1, 4, 4)]

        _d, mu_empty = solver._find_dist_group_zero_t(ws, vs, 0)
        _d, mu_full = solver._find_dist_group_zero_t(ws, vs, 4)
        self.assertAlmostEqual(mu_empty, -3.0, places=12)
        self.assertAlmostEqual(mu_full, 5.0, places=12)

    def test_uhfr_zero_t_fermi_level(self):
        from hwave.solver.uhfr import UHFr

        solver = UHFr.__new__(UHFr)
        w = np.array([-3.0, -1.0, 2.0, 5.0])
        self.assertAlmostEqual(solver._fermi_level_zero_t(w, 2), 0.5, places=12)
        self.assertAlmostEqual(solver._fermi_level_zero_t(w, 0), -3.0, places=12)
        self.assertAlmostEqual(solver._fermi_level_zero_t(w, 4), 5.0, places=12)


if __name__ == "__main__":
    unittest.main()
