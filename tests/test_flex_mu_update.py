#!/usr/bin/env python3

"""Tests for the self-consistent chemical-potential update in FLEX.

Standard self-consistent FLEX re-solves the chemical potential mu every SCF
iteration so that the DRESSED Green's function reproduces the target particle
number N = Ncond.  If mu is frozen at its non-interacting value, the growing
self-energy makes the actual particle number drift away from Ncond, so the
converged solution describes a different filling than the user requested.

These tests pin:
1. The dressed particle-number formula reduces to the analytic Fermi count at
   Sigma = 0 (correctness of the tail-corrected Matsubara sum).
2. After a full FLEX run with calc_mu=True and U > 0 on a DOPED (non
   particle-hole-symmetric) filling, the dressed Green's function carries
   exactly Ncond electrons.
3. When mu is fixed (calc_mu=False), it is NOT re-solved (behavior unchanged).
"""

import os
import unittest
import numpy as np


def _make_solver(fixing, Lx=8, Ly=8, Nmat=128, T=1.0, U=3.0,
                 iteration_max=1, mix=1.0):
    """Create a FLEX solver.  ``fixing`` is either {'mu': value} or
    {'Ncond': value} to select fixed-mu vs particle-number modes."""
    info_log = {}
    param = {
        'T': T,
        'CellShape': [Lx, Ly, 1],
        'SubShape': [1, 1, 1],
        'Nmat': Nmat,
        'IterationMax': iteration_max,
        'Mix': mix,
    }
    param.update(fixing)
    info_mode = {'mode': 'FLEX', 'param': param, 'calc_scheme': 'reduced'}
    info_file_input = {
        'path_to_input': 'tests/rpa/input',
        'interaction': {
            'path_to_input': 'tests/rpa/input',
            'Geometry': 'geom.dat',
            'Transfer': 'transfer.dat',
            'CoulombIntra': 'coulombintra.dat',
        },
    }

    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as solver_flex
    read_io = read_input_k.QLMSkInput(info_file_input)
    ham_info = read_io.get_param("ham")
    green_info = read_io.get_param("green")
    solver = solver_flex.FLEX(ham_info, info_log, info_mode)
    return solver, green_info


def _fermi(T, mu, ev, cutoff=1.0e2):
    w = (ev - mu) / T
    mask = w < cutoff
    w1 = np.where(mask, w, 0.0)
    return np.where(mask, 1.0 / (1.0 + np.exp(w1)), 0.0)


def _number_from_green(solver, green_kw, mu, beta):
    """Reference particle count from a stored dressed Green's function:
    N = sum_a f(eps_a - mu) + (1/beta) sum_{k,n} Tr[G - G0(mu)]."""
    nmat = green_kw.shape[1]
    ew = solver.H0_eigenvalue
    n_ref = _fermi(1.0 / beta, mu, ew).sum()
    iomega = (np.arange(nmat) * 2 + 1 - nmat) * np.pi / beta
    trG = np.einsum('bnkaa->bnk', green_kw)
    trG0 = (1.0 / ((1j * iomega)[np.newaxis, :, np.newaxis, np.newaxis]
                   + (mu - ew)[:, np.newaxis, :, :])).sum(axis=-1)
    return n_ref + ((trG - trG0).sum() / beta).real


class TestFlexMuUpdate(unittest.TestCase):

    def test_dressed_number_reduces_to_fermi_at_zero_sigma(self):
        """_calc_number_dressed with Sigma=0 must equal the analytic Fermi
        count that _find_mu uses (the tail-corrected sum is exact there)."""
        solver, green_info = _make_solver({'Ncond': 48.0})
        beta = 1.0 / solver.T
        solver._calc_epsilon_k(green_info)
        nblock, nvol, nd_block = solver.H0_eigenvalue.shape
        sigma0 = np.zeros((nblock, solver.nmat, nvol, nd_block, nd_block),
                          dtype=np.complex128)
        mu = 0.37

        n_dressed = solver._calc_number_dressed(sigma0, mu, beta)
        n_fermi = _fermi(solver.T, mu, solver.H0_eigenvalue).sum()

        self.assertAlmostEqual(n_dressed, n_fermi, places=10)

    def test_eig_number_matches_dressed_number(self):
        """The eigenvalue-accelerated particle count used inside the mu search
        must equal the reference inversion-based _calc_number_dressed to
        machine precision, for a NON-trivial (frequency-dependent, complex)
        self-energy.  This pins the optimization: G^{-1}(mu) = M + mu*I with M
        fixed, so Tr[G(mu)] = sum_j 1/(lam_j(M) + mu) needs the eigenvalues of
        M only once, not one inversion per trial mu."""
        solver, green_info = _make_solver({'Ncond': 40.0}, Nmat=64)
        beta = 1.0 / solver.T
        solver._calc_epsilon_k(green_info)
        nblock, nvol, nd_block = solver.H0_eigenvalue.shape

        rng = np.random.default_rng(1)
        sigma = 0.3 * (rng.standard_normal((nblock, solver.nmat, nvol,
                                            nd_block, nd_block))
                       + 1j * rng.standard_normal((nblock, solver.nmat, nvol,
                                                   nd_block, nd_block)))
        mu = 0.42

        n_ref = solver._calc_number_dressed(sigma, mu, beta)
        lam, ew = solver._matsubara_number_operator(sigma, beta)
        n_eig = solver._number_from_eigs(lam, ew, mu, beta)

        self.assertAlmostEqual(n_ref, n_eig, places=10)

    def test_number_derivative_matches_finite_difference(self):
        """dN/dmu returned by _number_from_eigs(with_deriv=True) must match a
        central finite difference -- Sigma is fixed during the mu search, so
        the derivative is analytic and drives a Newton solve."""
        # T != 1 so a wrong 1/T coefficient in dN_ref/dmu cannot hide
        solver, green_info = _make_solver({'Ncond': 40.0}, Nmat=64, T=2.0)
        beta = 1.0 / solver.T
        solver._calc_epsilon_k(green_info)
        nblock, nvol, nd_block = solver.H0_eigenvalue.shape

        rng = np.random.default_rng(3)
        sigma = 0.3 * (rng.standard_normal((nblock, solver.nmat, nvol,
                                            nd_block, nd_block))
                       + 1j * rng.standard_normal((nblock, solver.nmat, nvol,
                                                   nd_block, nd_block)))
        sigma = 0.5 * (sigma + sigma.conj().swapaxes(-1, -2))
        lam, ew = solver._matsubara_number_operator(sigma, beta)

        h = 1.0e-6
        for mu in (-2.0, -0.5, 0.3, 1.5):
            _, dN = solver._number_from_eigs(lam, ew, mu, beta,
                                             with_deriv=True)
            fd = (solver._number_from_eigs(lam, ew, mu + h, beta)
                  - solver._number_from_eigs(lam, ew, mu - h, beta)) / (2 * h)
            self.assertAlmostEqual(dN, fd, places=5,
                                   msg="dN/dmu mismatch at mu={}".format(mu))

    def test_particle_number_conserved_after_scf_doped(self):
        """After a converging FLEX run at a doped filling with U>0, the stored
        dressed Green's function must carry exactly Ncond electrons -- the
        whole point of re-solving mu each iteration."""
        Ncond = 40.0   # Nstate = 8*8*2 = 128 -> filling 0.3125 (doped)
        solver, green_info = _make_solver({'Ncond': Ncond}, U=3.0,
                                          iteration_max=60, mix=0.4)
        os.makedirs('tests/flex/output', exist_ok=True)
        solver.solve(green_info, 'tests/flex/output')

        beta = 1.0 / solver.T
        green_kw = green_info["green"]
        # spin-free counts one spin -> target is Ncond/2
        n_actual = _number_from_green(solver, green_kw, solver.mu, beta)

        self.assertAlmostEqual(n_actual, Ncond / 2.0, places=4)

    def test_fixed_mu_is_not_resolved(self):
        """calc_mu=False: mu stays at the user value through the whole run."""
        mu_in = 0.5
        solver, green_info = _make_solver({'mu': mu_in}, U=3.0,
                                          iteration_max=5, mix=0.5)
        os.makedirs('tests/flex/output', exist_ok=True)
        solver.solve(green_info, 'tests/flex/output')

        self.assertEqual(solver.mu, mu_in)

    def test_mu_search_in_charge_gap_is_flat_safe(self):
        """When a charge gap makes N(mu) flat (dN/dmu -> 0) at the target
        filling, the Newton-based mu search must not blow up (f/df with a
        vanishing derivative): it must still return a finite mu inside the gap
        with N(mu) == Ncond.  Built on a synthetic gapped spectrum (half the
        k-points at -1, half at +1 -> spectral gap (-1, +1)) filled to the gap
        edge, at low T so the plateau is sharp."""
        solver, green_info = _make_solver({'Ncond': 40.0}, Nmat=64, T=0.02)
        solver._calc_epsilon_k(green_info)
        nb, nv, nd = solver.H0_eigenvalue.shape

        ew = np.zeros((nb, nv, nd))
        ew[0, :nv // 2, 0] = -1.0
        ew[0, nv // 2:, 0] = +1.0
        solver.H0_eigenvalue = ew
        solver.H0_eigenvector = np.ones((nb, nv, nd, nd), dtype=np.complex128)

        beta = 1.0 / solver.T
        sigma = np.zeros((nb, solver.nmat, nv, nd, nd), dtype=np.complex128)
        Ncond_target = float(nv // 2)   # fill the lower group -> mu in the gap

        lam, ewr = solver._matsubara_number_operator(sigma, beta)
        _, slope = solver._number_from_eigs(lam, ewr, 0.0, beta,
                                            with_deriv=True)
        self.assertLess(abs(slope), 1.0e-6)   # confirm the gap is actually flat

        mu = solver._find_mu_dressed(sigma, beta, Ncond_target)
        n_at_mu = solver._number_from_eigs(lam, ewr, mu, beta)

        self.assertTrue(np.isfinite(mu))
        self.assertTrue(-1.0 < mu < 1.0, "mu={} not inside the gap".format(mu))
        self.assertAlmostEqual(n_at_mu, Ncond_target, places=6)


if __name__ == '__main__':
    unittest.main()
