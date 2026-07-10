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
    # The fixture's CoulombIntra file pins U=4; override every on-site entry so
    # the requested U actually drives the interacting behavior under test
    # (otherwise the U argument would be silently ignored).
    if "CoulombIntra" in ham_info:
        ham_info["CoulombIntra"] = {k: complex(U)
                                    for k in ham_info["CoulombIntra"]}
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

    def test_output_triple_is_consistent(self):
        """The stored (mu, green, sigma) must satisfy the Dyson equation:
        rebuilding G from green_info['sigma'] at solver.mu must reproduce the
        stored green AND carry Ncond electrons -- even for a short (non-fully-
        converged) run.  Guards the final-consistency rebuild."""
        Ncond = 40.0
        solver, green_info = _make_solver({'Ncond': Ncond}, U=3.0,
                                          iteration_max=3, mix=0.4)
        os.makedirs('tests/flex/output', exist_ok=True)
        solver.solve(green_info, 'tests/flex/output')

        beta = 1.0 / solver.T
        green_from_sigma = solver._calc_dressed_green(beta, solver.mu,
                                                      green_info["sigma"])
        # stored green equals G rebuilt from the stored sigma at the stored mu
        np.testing.assert_allclose(green_from_sigma, green_info["green"],
                                   rtol=0.0, atol=1.0e-12)
        # and that consistent triple carries the target particle number
        n_from_sigma = _number_from_green(solver, green_from_sigma, solver.mu,
                                          beta)
        self.assertAlmostEqual(n_from_sigma, Ncond / 2.0, places=4)

    def test_find_mu_brackets_root_beyond_initial_guess(self):
        """A near-empty target pushes mu far below the initial bracket
        [band_min - pad, ...]; _find_mu_dressed must EXPAND the bracket and
        still return a finite mu with N(mu) == target, not sys.exit.  This
        exercises the adaptive-bracket path that also covers multi-orbital /
        off-diagonal self-energies whose root moves beyond max|Re Sigma|."""
        solver, green_info = _make_solver({'Ncond': 40.0}, Nmat=64, T=1.0)
        solver._calc_epsilon_k(green_info)
        beta = 1.0 / solver.T
        nb, nv, nd = solver.H0_eigenvalue.shape
        sigma = np.zeros((nb, solver.nmat, nv, nd, nd), dtype=np.complex128)

        target = 1.0e-3   # nearly empty -> mu far below band_min - 1
        mu = solver._find_mu_dressed(sigma, beta, target)
        n_at_mu = solver._number_from_eigs(
            *solver._matsubara_number_operator(sigma, beta), mu, beta)

        self.assertTrue(np.isfinite(mu))
        self.assertLess(mu, float(solver.H0_eigenvalue.min()) - 1.0)
        self.assertAlmostEqual(n_at_mu, target, places=6)

    def test_find_mu_unbracketable_target_raises(self):
        """A target above the total number of states can never be bracketed
        (N(mu) -> Nstate < target as mu -> +inf).  _find_mu_dressed must raise
        a catchable RuntimeError (not sys.exit, which would tear down the whole
        interpreter and break sweeps/notebooks)."""
        solver, green_info = _make_solver({'Ncond': 40.0}, Nmat=64, T=1.0)
        solver._calc_epsilon_k(green_info)
        beta = 1.0 / solver.T
        nb, nv, nd = solver.H0_eigenvalue.shape
        sigma = np.zeros((nb, solver.nmat, nv, nd, nd), dtype=np.complex128)

        nstate = int(np.prod(solver.H0_eigenvalue.shape))
        with self.assertRaises(RuntimeError):
            solver._find_mu_dressed(sigma, beta, float(nstate) * 10.0)

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


class TestFlexPhysicsOutput(unittest.TestCase):
    """N (particle number) and Sz reporting from the dressed Green function."""

    def test_occupation_sum_matches_number_dressed(self):
        """Per-orbital occupation summed over orbitals must equal the traced
        _calc_number_dressed to machine precision (nonzero sigma): the
        occupation is just the trace resolved per orbital."""
        solver, green_info = _make_solver({'Ncond': 40.0}, Nmat=64)
        beta = 1.0 / solver.T
        solver._calc_epsilon_k(green_info)
        nblock, nvol, nd_block = solver.H0_eigenvalue.shape

        rng = np.random.default_rng(7)
        sigma = 0.3 * (rng.standard_normal((nblock, solver.nmat, nvol,
                                            nd_block, nd_block))
                       + 1j * rng.standard_normal((nblock, solver.nmat, nvol,
                                                   nd_block, nd_block)))
        mu = 0.35
        green_kw = solver._calc_dressed_green(beta, mu, sigma)

        nocc = solver._calc_occupation_dressed(green_kw, mu, beta)
        n_ref = solver._calc_number_dressed(sigma, mu, beta)

        self.assertAlmostEqual(nocc.sum(), n_ref, places=10)

    def test_physics_spin_free_N_and_zero_Sz(self):
        """Spin-free run: NCond == 2 * one-spin count, Sz == 0 by construction."""
        solver, green_info = _make_solver({'Ncond': 40.0}, U=3.0,
                                          iteration_max=5, mix=0.4)
        os.makedirs('tests/flex/output', exist_ok=True)
        solver.solve(green_info, 'tests/flex/output')

        beta = 1.0 / solver.T
        physics = green_info["physics"]
        n_one_spin = solver._calc_number_dressed(green_info["sigma"],
                                                 solver.mu, beta)

        self.assertAlmostEqual(physics["NCond"], 2.0 * n_one_spin, places=8)
        self.assertEqual(physics["Sz"], 0.0)
        self.assertEqual(physics["mu"], solver.mu)

    def test_physics_N_equals_Ncond_after_calc_mu(self):
        """Converged calc_mu doped run: reported NCond equals the target Ncond."""
        Ncond = 40.0
        solver, green_info = _make_solver({'Ncond': Ncond}, U=3.0,
                                          iteration_max=60, mix=0.4)
        os.makedirs('tests/flex/output', exist_ok=True)
        solver.solve(green_info, 'tests/flex/output')

        self.assertAlmostEqual(green_info["physics"]["NCond"], Ncond, places=4)

    def test_fixed_mu_outputs_N_to_energy_file(self):
        """calc_mu=False (fixed mu) writes energy.dat with a finite NCond and
        ChemicalPotential == the input mu -- one point of the mu-N curve."""
        mu_in = 0.5
        solver, green_info = _make_solver({'mu': mu_in}, U=3.0,
                                          iteration_max=5, mix=0.5)
        out_dir = 'tests/flex/output'
        os.makedirs(out_dir, exist_ok=True)
        solver.solve(green_info, out_dir)
        solver.save_results({'path_to_output': out_dir, 'energy': 'energy.dat'},
                            green_info)

        with open(os.path.join(out_dir, 'energy.dat')) as fr:
            lines = dict(line.split('=', 1) for line in fr
                         if '=' in line)
        ncond = float(lines['NCond '])
        mu_out = float(lines['ChemicalPotential '])

        self.assertTrue(np.isfinite(ncond))
        self.assertGreater(ncond, 0.0)
        self.assertAlmostEqual(mu_out, mu_in, places=12)

    def test_physics_spin_diag_Sz(self):
        """Spin-diag: Sz == (N_up - N_down)/2 with N_up/N_down the per-block
        counts.  Built on a synthetic magnetized 2-block spectrum (up/down bands
        split) so N_up != N_down, checked against the analytic Fermi counts."""
        solver, green_info = _make_solver({'Ncond': 40.0}, Nmat=64, T=0.5)
        solver._calc_epsilon_k(green_info)
        _, nvol, _ = solver.H0_eigenvalue.shape

        # spin-diag: 2 blocks (up, down), one orbital, split bands -> magnetized
        solver.spin_mode = "spin-diag"
        solver.norb = 1
        eps_up = np.linspace(-2.0, 2.0, nvol) - 0.5     # up band shifted down
        eps_down = np.linspace(-2.0, 2.0, nvol) + 0.5   # down band shifted up
        ew = np.stack([eps_up[:, None], eps_down[:, None]])   # (2, nvol, 1)
        solver.H0_eigenvalue = ew
        solver.H0_eigenvector = np.ones((2, nvol, 1, 1), dtype=np.complex128)

        beta = 1.0 / solver.T
        mu = 0.0
        sigma = np.zeros((2, solver.nmat, nvol, 1, 1), dtype=np.complex128)
        green_kw = solver._calc_dressed_green(beta, mu, sigma)

        physics = solver._calc_physics_dressed(green_kw, mu, beta)

        # analytic reference at Sigma=0: N_spin = sum_k f(eps_spin - mu)
        n_up = _fermi(solver.T, mu, eps_up).sum()
        n_down = _fermi(solver.T, mu, eps_down).sum()
        self.assertAlmostEqual(physics["NCond"], n_up + n_down, places=8)
        self.assertAlmostEqual(physics["Sz"], 0.5 * (n_up - n_down), places=8)
        self.assertNotAlmostEqual(physics["Sz"], 0.0, places=3)  # truly magnetized


class TestEigvalsSmall(unittest.TestCase):
    """Closed-form batched eigenvalues for the mu-search operator M.

    The mu search only consumes the eigenvalues as an (unordered) set --
    N(mu) = sum_j 1/(lam_j + mu) -- so the closed-form nd<=2 path must
    reproduce LAPACK geev's spectra as multisets, on the input's array
    backend (avoiding both the host round-trip and the batched-geev cost)."""

    def _assert_same_spectra(self, lam_a, lam_b, nd):
        a = np.sort_complex(np.asarray(lam_a).reshape(-1, nd))
        b = np.sort_complex(np.asarray(lam_b).reshape(-1, nd))
        np.testing.assert_allclose(a, b, atol=1e-12)

    def test_closed_form_1x1_matches_geev(self):
        from hwave.solver.flex import _eigvals_small
        rng = np.random.default_rng(3)
        M = (rng.standard_normal((2, 5, 7, 1, 1))
             + 1j * rng.standard_normal((2, 5, 7, 1, 1)))
        self._assert_same_spectra(_eigvals_small(M), np.linalg.eigvals(M), 1)

    def test_closed_form_2x2_matches_geev(self):
        from hwave.solver.flex import _eigvals_small
        rng = np.random.default_rng(4)
        M = (rng.standard_normal((3, 4, 6, 2, 2))
             + 1j * rng.standard_normal((3, 4, 6, 2, 2)))
        self._assert_same_spectra(_eigvals_small(M), np.linalg.eigvals(M), 2)

    def test_geev_fallback_3x3(self):
        from hwave.solver.flex import _eigvals_small
        rng = np.random.default_rng(5)
        M = (rng.standard_normal((2, 3, 3, 3))
             + 1j * rng.standard_normal((2, 3, 3, 3)))
        self._assert_same_spectra(_eigvals_small(M), np.linalg.eigvals(M), 3)

    def test_number_from_eigs_returns_plain_floats(self):
        """The mu-search closure consumes N (and dN/dmu) in host-side float
        comparisons, so _number_from_eigs must return plain Python floats."""
        solver, _ = _make_solver({'Ncond': 0.8})
        nvol = solver.lattice.nvol
        rng = np.random.default_rng(6)
        solver.H0_eigenvalue = rng.standard_normal((1, nvol, 1))
        beta = 1.0 / solver.T
        sigma = 0.1 * (rng.standard_normal((1, solver.nmat, nvol, 1, 1))
                       + 1j * rng.standard_normal((1, solver.nmat, nvol, 1, 1)))
        # keep H0_eigenvector consistent for _matsubara_number_operator
        solver.H0_eigenvector = np.ones((1, nvol, 1, 1), dtype=np.complex128)
        lam, ew = solver._matsubara_number_operator(sigma, beta)
        n = solver._number_from_eigs(lam, ew, 0.1, beta)
        n2, dn = solver._number_from_eigs(lam, ew, 0.1, beta, with_deriv=True)
        self.assertIsInstance(n, float)
        self.assertIsInstance(n2, float)
        self.assertIsInstance(dn, float)
        self.assertEqual(n, n2)


def test_mu_search_gpu_matches_cpu():
    """On a CUDA device the whole dressed mu search (closed-form eigenvalues
    + device-resident N(mu) evaluations) must reproduce the CPU mu."""
    import pytest
    cupy = pytest.importorskip("cupy")
    try:
        cupy.zeros(1)
    except Exception:
        pytest.skip("cupy installed but no usable CUDA device")

    solver, _ = _make_solver({'Ncond': 0.8}, Lx=4, Ly=4, Nmat=32)
    nvol = solver.lattice.nvol
    rng = np.random.default_rng(7)
    solver.H0_eigenvalue = rng.standard_normal((1, nvol, 1))
    solver.H0_eigenvector = np.ones((1, nvol, 1, 1), dtype=np.complex128)
    beta = 1.0 / solver.T
    sigma = 0.2 * (rng.standard_normal((1, solver.nmat, nvol, 1, 1))
                   + 1j * rng.standard_normal((1, solver.nmat, nvol, 1, 1)))

    mu_cpu = solver._find_mu_dressed(sigma, beta, 0.8)

    solver.H0_eigenvalue = cupy.asarray(solver.H0_eigenvalue)
    solver.H0_eigenvector = cupy.asarray(solver.H0_eigenvector)
    mu_gpu = solver._find_mu_dressed(cupy.asarray(sigma), beta, 0.8)
    assert np.isclose(mu_gpu, mu_cpu, atol=1e-10)


if __name__ == '__main__':
    unittest.main()
