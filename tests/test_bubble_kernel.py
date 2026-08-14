"""Gates for the unified particle-hole bubble kernel (Step 2 of the H-wave
bubble-unification series).

Binding spec: ``docs/superpowers/specs/2026-08-14-unified-bubble-kernel-design.md``.
This module grows through the series' tasks; it starts with the Green
builder (``hwave.solver.green.build_green``, extracted from
``RPA._calc_green`` in Task 1).
"""

import unittest

import numpy as np

from tests.approx_util import assert_approx_array
from hwave.solver import green as green_mod


def _tiny_eig(nvol=4, p=2, seed=7):
    rng = np.random.default_rng(seed)
    ev = rng.uniform(-2.0, 2.0, size=(nvol, p))
    # random unitary per k (QR of a random complex matrix)
    m = rng.normal(size=(nvol, p, p)) + 1j * rng.normal(size=(nvol, p, p))
    q, _ = np.linalg.qr(m)
    return ev, q


def _make_rpa_solver(coeff_tail=0.0, Lx=4, Ly=4, Nmat=8):
    """Build a minimal 1-orbital RPA solver from the ``tests/rpa/input``
    fixture (same construction pattern as
    ``tests/test_chi0q_coeff_tail_provenance.py``'s ``_make_rpa`` --
    ``tests/test_rpa_chi0q_guards.py`` only builds a bare
    ``object.__new__`` stub for its shape-guard tests, so this module
    copies the fuller construction pattern used elsewhere in the suite
    for tests that need a real eigen-decomposition)."""
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.rpa as sol_rpa

    info_input = {
        'path_to_input': 'tests/rpa/input',
        'interaction': {
            'path_to_input': 'tests/rpa/input',
            'Geometry': 'geom.dat',
            'Transfer': 'transfer.dat',
            'CoulombIntra': 'coulombintra.dat',
        },
    }
    param = {'T': 2.0, 'mu': 0.0,
             'CellShape': [Lx, Ly, 1], 'SubShape': [1, 1, 1], 'Nmat': Nmat,
             'coeff_tail': coeff_tail}
    info_mode = {'mode': 'RPA', 'param': param, 'calc_scheme': 'general'}
    ham = read_input_k.QLMSkInput(info_input).get_param("ham")
    solver = sol_rpa.RPA(ham, {}, info_mode)
    solver._init_wavevec()
    return solver


class TestBuildGreen(unittest.TestCase):

    def test_formula_pin_against_direct_sum(self):
        ev, V = _tiny_eig()
        beta, nmat, mu, ct = 2.0, 8, 0.1, 1.0
        full, defl, tail = green_mod.build_green(ev, V, mu, beta, nmat, ct)
        self.assertEqual(full.shape, (1, nmat, ev.shape[0], 2, 2))
        n, k = 3, 1
        iw = 1j * (2 * n + 1 - nmat) * np.pi / beta
        g_direct = (V[k] * (1.0 / (iw - (ev[k] - mu)))) @ V[k].conj().T
        assert_approx_array(full[0, n, k], g_direct, rel=0, abs=1e-13)
        assert_approx_array(full[0, n, k] - defl[0, n, k],
                             ct * (V[k] @ V[k].conj().T) / iw, rel=0, abs=1e-13)
        assert_approx_array(tail[0, n, k],
                             V[k] @ V[k].conj().T * ct * beta / 2, rel=0, abs=1e-13)

    def test_coeff_tail_zero_aliases(self):
        ev, V = _tiny_eig()
        full, defl, tail = green_mod.build_green(ev, V, 0.0, 2.0, 8, 0.0)
        self.assertIs(full, defl)
        self.assertIsNone(tail)

    def test_matches_rpa_solver_calc_green(self):
        """Oracle: RPA._calc_green on the same eigen-decomposition must
        equal build_green's (deflated, tail) pair exactly (same math,
        transcribed -- not just close)."""
        solver = _make_rpa_solver(coeff_tail=0.5)
        solver._calc_epsilon_k({})
        beta = 1.0 / solver.T
        mu = solver.mu_value
        green_old, tail_old = solver._calc_green(beta, mu)

        full_new, defl_new, tail_new = green_mod.build_green(
            solver.H0_eigenvalue, solver.H0_eigenvector, mu, beta,
            solver.nmat, solver.coeff_tail)

        assert_approx_array(defl_new, green_old, rel=0, abs=1e-13)
        assert_approx_array(tail_new, tail_old, rel=0, abs=1e-13)

    def test_rejects_unflattened_shapes(self):
        """sc.py's genuinely UNFLATTENED (Nx, Ny, Nz, p) /
        (Nx, Ny, Nz, p, p) layout (ndim 4 / 5) must be rejected -- the
        caller is required to flatten to (nvol, p) / (nvol, p, p) first.

        (A reshape that stays within the accepted ndim pairs -- e.g. to
        (nblock, nvol, p) / (nblock, nvol, p, p) -- is not a case this
        test can use: for nvol=4, p=2 such a reshape is numerically
        self-consistent as a genuine 2-block decomposition of the same
        data (verified: each reshaped (block, k) slice reproduces one
        original (k) slice exactly), so it is legitimately ACCEPTED, not
        a shape to reject. The real invalid case is the unflattened
        spatial layout build_green documents as out of contract.)"""
        ev, V = _tiny_eig()
        with self.assertRaises(ValueError):
            green_mod.build_green(ev.reshape(2, 2, 1, 2), V.reshape(2, 2, 1, 2, 2),
                                   0.0, 2.0, 8, 0.0)

    def test_rejects_mismatched_p(self):
        ev, V = _tiny_eig(nvol=4, p=2)
        bad_V = np.zeros((4, 3, 3), dtype=complex)
        with self.assertRaises(ValueError):
            green_mod.build_green(ev, bad_V, 0.0, 2.0, 8, 0.0)

    def test_rejects_nonpositive_beta(self):
        ev, V = _tiny_eig()
        with self.assertRaises(ValueError):
            green_mod.build_green(ev, V, 0.0, 0.0, 8, 0.0)
        with self.assertRaises(ValueError):
            green_mod.build_green(ev, V, 0.0, -1.0, 8, 0.0)

    def test_rejects_odd_nmat(self):
        ev, V = _tiny_eig()
        with self.assertRaises(ValueError):
            green_mod.build_green(ev, V, 0.0, 2.0, 7, 0.0)

    def test_rejects_nonfinite_coeff_tail(self):
        ev, V = _tiny_eig()
        with self.assertRaises(ValueError):
            green_mod.build_green(ev, V, 0.0, 2.0, 8, float('nan'))
        with self.assertRaises(ValueError):
            green_mod.build_green(ev, V, 0.0, 2.0, 8, float('inf'))


if __name__ == "__main__":
    unittest.main()
