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


if __name__ == '__main__':
    unittest.main()
