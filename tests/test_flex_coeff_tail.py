#!/usr/bin/env python3

"""Tests for the Matsubara high-frequency tail (coeff_tail) handling in FLEX.

RPA's tail acceleration is a two-step contract:
  1. _calc_green subtracts aa/(i w_n) in frequency space, and
  2. _calc_chi0q subtracts the analytic tau-space constant green0_tail.
FLEX builds its Green's function via _calc_dressed_green instead of
_calc_green, so step 1 must be applied by FLEX itself before handing the
Green's function to _calc_chi0q / the self-energy convolution.  These tests
pin that contract: with coeff_tail != 0 the FLEX results must remain a small
tail-acceleration *correction*, not a change of physics.
"""

import os
import unittest
import numpy as np


def _make_solver(mode_cls, Lx=8, Ly=8, Nmat=64, T=2.0, mu=0.0,
                 coeff_tail=0.0, iteration_max=1, mix=1.0):
    """Create an RPA or FLEX solver on the shared 1-orbital Hubbard fixture."""
    info_log = {}
    info_mode = {
        'mode': mode_cls,
        'param': {
            'T': T,
            'mu': mu,
            'CellShape': [Lx, Ly, 1],
            'SubShape': [1, 1, 1],
            'Nmat': Nmat,
            'coeff_tail': coeff_tail,
            'IterationMax': iteration_max,
            'Mix': mix,
        },
        'calc_scheme': 'reduced',
    }
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
    read_io = read_input_k.QLMSkInput(info_file_input)
    ham_info = read_io.get_param("ham")
    green_info = read_io.get_param("green")

    if mode_cls == "FLEX":
        import hwave.solver.flex as solver_flex
        solver = solver_flex.FLEX(ham_info, info_log, info_mode)
    else:
        import hwave.solver.rpa as solver_rpa
        solver = solver_rpa.RPA(ham_info, info_log, info_mode)

    return solver, green_info


def _run_flex(Nmat, coeff_tail):
    """Run one FLEX iteration (Sigma_in = 0) and return green_info."""
    solver, green_info = _make_solver("FLEX", Nmat=Nmat,
                                      coeff_tail=coeff_tail,
                                      iteration_max=1, mix=1.0)
    os.makedirs('tests/flex/output', exist_ok=True)
    solver.solve(green_info, 'tests/flex/output')
    return green_info


class TestFlexCoeffTail(unittest.TestCase):
    """coeff_tail must be a small acceleration correction in FLEX, as in RPA."""

    def test_chi0q_first_iteration_matches_rpa(self):
        """With Sigma=0 the FLEX chi0q must equal the RPA chi0q for the SAME
        coeff_tail: both describe the identical bare bubble, so any deviation
        beyond matrix-inversion roundoff means the tail contract is broken."""
        coeff_tail = 1.0

        rpa, rpa_info = _make_solver("RPA", coeff_tail=coeff_tail)
        beta = 1.0 / rpa.T
        rpa._calc_epsilon_k(rpa_info)
        green0, green0_tail = rpa._calc_green(beta, rpa.mu_value)
        chi0q_rpa = rpa._calc_chi0q(green0, green0_tail, beta)[0]

        flex_info = _run_flex(Nmat=64, coeff_tail=coeff_tail)
        # FLEX outputs the spin-inflated chi0q; the orbital block sits on the
        # spin-diagonal (1 orbital here: element [0, 0] of the 2x2 block).
        chi0q_flex = flex_info["chi0q"][..., 0:1, 0:1]

        np.testing.assert_allclose(chi0q_flex, chi0q_rpa,
                                   rtol=0.0, atol=1.0e-10)

    def test_chi0q_coeff_tail_is_small_correction(self):
        """chi0q(coeff_tail=1) vs chi0q(coeff_tail=0) must differ only by
        the finite-Nmat truncation error (~1e-2 here), never by O(1).

        Since the #134 endpoint fix the corrected result sits essentially
        on the exact value, so this difference IS coeff_tail=0's full
        truncation error rather than a partial cancellation with the
        old (broken) tail -- the bound reflects that scale."""
        chi0q_0 = _run_flex(Nmat=64, coeff_tail=0.0)["chi0q"]
        chi0q_1 = _run_flex(Nmat=64, coeff_tail=1.0)["chi0q"]

        scale = np.abs(chi0q_0).max()
        diff = np.abs(chi0q_1 - chi0q_0).max()
        self.assertLess(diff, 0.15 * scale,
                        "coeff_tail changed chi0q by {:.3e} (scale {:.3e}); "
                        "tail handling is inconsistent".format(diff, scale))

    def test_sigma_converges_to_same_limit_regardless_of_coeff_tail(self):
        """coeff_tail must not change the physical self-energy.  The self-energy
        convolution keeps the full physical G (no tail reconstruction), so at a
        large Nmat the two conventions must agree to the finite-Nmat truncation
        error -- i.e. coeff_tail leaves Sigma essentially unchanged, never O(1)
        wrong (the pre-fix bug shifted G(tau) by -aa/2 for chi0q; this guards
        that the sigma path is not similarly corrupted)."""
        n_big = 512
        sigma_0 = _run_flex(Nmat=n_big, coeff_tail=0.0)["sigma"]
        sigma_1 = _run_flex(Nmat=n_big, coeff_tail=1.0)["sigma"]

        scale = np.abs(sigma_0).max()
        diff = np.abs(sigma_1 - sigma_0).max()
        self.assertLess(diff, 0.02 * scale,
                        "coeff_tail changed the self-energy by {:.3e} "
                        "(scale {:.3e}); the sigma path is tail-corrupted"
                        .format(diff, scale))


if __name__ == '__main__':
    unittest.main()
