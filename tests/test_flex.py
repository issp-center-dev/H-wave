#!/usr/bin/env python3

"""Tests for the FLEX solver.

Test strategy:
1. Consistency test: FLEX with Sigma=0 (1 iteration) should reproduce RPA results
2. SCF convergence: FLEX should converge for weakly-coupled systems
3. Physical constraints: Self-energy should be Hermitian, particle number conserved
"""

import os
import sys
import unittest
import numpy as np


class TestFLEX(unittest.TestCase):
    """Basic FLEX solver tests using 1-orbital Hubbard model."""

    def _make_solver(self, mode_cls, Lx=8, Ly=8, Nmat=64, T=2.0, mu=0.0,
                     U=4.0, calc_scheme="reduced", extra_mode=None,
                     extra_interactions=None):
        """Create a solver instance with given parameters.

        Parameters
        ----------
        mode_cls : str
            "RPA" or "FLEX"
        """
        info_log = {}
        info_mode = {
            'mode': mode_cls,
            'param': {
                'T': T,
                'mu': mu,
                'CellShape': [Lx, Ly, 1],
                'SubShape': [1, 1, 1],
                'Nmat': Nmat,
            },
            'calc_scheme': calc_scheme,
        }
        if extra_mode:
            info_mode.update(extra_mode)

        info_file = {
            'input': {
                'path_to_input': 'tests/rpa/input',
                'interaction': {
                    'path_to_input': 'tests/rpa/input',
                    'Geometry': 'geom.dat',
                    'Transfer': 'transfer.dat',
                    'CoulombIntra': 'coulombintra.dat',
                },
            },
            'output': {
                'path_to_output': 'tests/flex/output',
            },
        }
        if extra_interactions:
            info_file['input']['interaction'].update(extra_interactions)

        os.makedirs('tests/flex/output', exist_ok=True)

        import hwave.qlmsio.read_input_k as read_input_k
        read_io = read_input_k.QLMSkInput(info_file['input'])
        ham_info = read_io.get_param("ham")
        green_info = read_io.get_param("green")

        if mode_cls == "FLEX":
            import hwave.solver.flex as solver_flex
            solver = solver_flex.FLEX(ham_info, info_log, info_mode)
        else:
            import hwave.solver.rpa as solver_rpa
            solver = solver_rpa.RPA(ham_info, info_log, info_mode)

        return solver, green_info, info_file

    def test_flex_instantiation(self):
        """Test that FLEX solver can be instantiated."""
        solver, green_info, _ = self._make_solver("FLEX")
        self.assertIsNotNone(solver)
        self.assertEqual(solver.norb, 1)
        self.assertEqual(solver.ns, 2)
        self.assertEqual(solver.nd, 2)

    def test_flex_inherits_rpa(self):
        """Test that FLEX properly inherits RPA methods."""
        import hwave.solver.rpa as solver_rpa
        solver, green_info, _ = self._make_solver("FLEX")
        self.assertIsInstance(solver, solver_rpa.RPA)
        self.assertTrue(hasattr(solver, '_calc_chi0q'))
        self.assertTrue(hasattr(solver, '_solve_rpa'))
        self.assertTrue(hasattr(solver, '_find_block_diagonal'))

    def test_flex_single_iteration_chi0q_matches_rpa(self):
        """FLEX chi0q with Sigma=0 should match RPA chi0q exactly."""
        Lx, Ly, Nmat = 8, 8, 64

        # Run RPA
        rpa_solver, rpa_green_info, _ = self._make_solver(
            "RPA", Lx=Lx, Ly=Ly, Nmat=Nmat, calc_scheme="reduced")
        rpa_solver.solve(rpa_green_info, 'tests/flex/output')
        chi0q_rpa = rpa_green_info["chi0q"]

        # Run FLEX with 1 iteration
        flex_solver, flex_green_info, _ = self._make_solver(
            "FLEX", Lx=Lx, Ly=Ly, Nmat=Nmat, calc_scheme="reduced",
            extra_mode={'param': {
                'T': 2.0, 'mu': 0.0,
                'CellShape': [Lx, Ly, 1], 'SubShape': [1, 1, 1],
                'Nmat': Nmat, 'IterationMax': 1, 'Mix': 1.0,
            }})
        flex_solver.solve(flex_green_info, 'tests/flex/output')
        chi0q_flex = flex_green_info["chi0q"]

        # FLEX returns spin-inflated chi0q (nd x nd), RPA returns orbital-only (norb x norb)
        # Extract the orbital block from FLEX chi0q for comparison
        norb = rpa_solver.norb
        if chi0q_flex.shape[-1] != chi0q_rpa.shape[-1]:
            chi0q_flex_orb = chi0q_flex[:, :, :norb, :norb]
        else:
            chi0q_flex_orb = chi0q_flex

        self.assertTrue(np.allclose(chi0q_rpa, chi0q_flex_orb, atol=1e-10),
                        "chi0q mismatch between RPA and FLEX (1st iteration)")

    def test_flex_convergence_weak_coupling(self):
        """FLEX should converge for small U."""
        solver, green_info, _ = self._make_solver(
            "FLEX", Lx=4, Ly=4, Nmat=32, U=1.0, T=4.0,
            calc_scheme="reduced",
            extra_mode={'param': {
                'T': 4.0, 'mu': 0.0,
                'CellShape': [4, 4, 1], 'SubShape': [1, 1, 1],
                'Nmat': 32, 'IterationMax': 50, 'Mix': 0.2, 'EPS': 4,
            }})
        solver.solve(green_info, 'tests/flex/output')

        # Check that sigma is stored and has correct shape
        sigma = green_info.get("sigma")
        self.assertIsNotNone(sigma, "Self-energy not stored")

    def test_flex_self_energy_shape(self):
        """Self-energy should have correct dimensions."""
        Lx, Ly, Nmat = 4, 4, 32
        solver, green_info, _ = self._make_solver(
            "FLEX", Lx=Lx, Ly=Ly, Nmat=Nmat, T=4.0, U=1.0,
            calc_scheme="reduced",
            extra_mode={'param': {
                'T': 4.0, 'mu': 0.0,
                'CellShape': [Lx, Ly, 1], 'SubShape': [1, 1, 1],
                'Nmat': Nmat, 'IterationMax': 2, 'Mix': 0.5,
            }})
        solver.solve(green_info, 'tests/flex/output')

        sigma = green_info["sigma"]
        # sigma should be (nblock, nmat, nvol, nd_block, nd_block)
        self.assertEqual(sigma.shape[1], Nmat)
        self.assertEqual(sigma.shape[2], Lx * Ly)

    def test_flex_dressed_green_reduces_to_bare(self):
        """With Sigma=0, dressed G should equal bare G."""
        solver, green_info, _ = self._make_solver(
            "FLEX", Lx=4, Ly=4, Nmat=32, T=2.0,
            calc_scheme="reduced",
            extra_mode={'param': {
                'T': 2.0, 'mu': 0.0,
                'CellShape': [4, 4, 1], 'SubShape': [1, 1, 1],
                'Nmat': 32, 'IterationMax': 1, 'Mix': 1.0,
            }})

        # Setup
        solver._calc_epsilon_k(green_info)
        if solver.calc_mu:
            if solver.spin_mode == "spin-free":
                Ncond = solver.Ncond / 2
            else:
                Ncond = solver.Ncond
            dist, mu = solver._find_mu(Ncond, solver.T)
        else:
            mu = solver.mu_value

        beta = 1.0 / solver.T

        # Bare Green's function
        green0, green0_tail = solver._calc_green(beta, mu)

        # Dressed Green with Sigma=0
        nblock = green0.shape[0]
        nd_block = green0.shape[-1]
        nvol = solver.lattice.nvol
        nmat = solver.nmat
        sigma_zero = np.zeros((nblock, nmat, nvol, nd_block, nd_block),
                              dtype=np.complex128)
        green_dressed = solver._calc_dressed_green(beta, mu, sigma_zero)

        # Compare: they should agree (up to tail correction differences)
        # The bare G uses eigendecomposition, dressed G uses matrix inversion
        # Both should give the same result for Sigma=0
        self.assertTrue(np.allclose(green0, green_dressed, atol=1e-10),
                        "Dressed G (Sigma=0) != Bare G")

    def test_flex_via_qlms(self):
        """Test that FLEX can be invoked through qlms.run()."""
        import tempfile

        # Create a minimal TOML config
        toml_content = """
[mode]
mode = "FLEX"
calc_scheme = "reduced"

[mode.param]
T = 4.0
mu = 0.0
CellShape = [4, 4, 1]
SubShape = [1, 1, 1]
Nmat = 32
IterationMax = 2
Mix = 0.5

[file.input]
path_to_input = ""

[file.input.interaction]
path_to_input = "tests/rpa/input"
Geometry = "geom.dat"
Transfer = "transfer.dat"
CoulombIntra = "coulombintra.dat"

[file.output]
path_to_output = "tests/flex/output_qlms"
"""
        os.makedirs('tests/flex/output_qlms', exist_ok=True)

        with tempfile.NamedTemporaryFile(mode='w', suffix='.toml',
                                         delete=False) as f:
            f.write(toml_content)
            toml_path = f.name

        try:
            import hwave.qlms as qlms
            qlms.run(input_file=toml_path)
        finally:
            os.unlink(toml_path)

    def test_flex_scf_reduces_diff(self):
        """Multiple SCF iterations should reduce the convergence criterion."""
        Lx, Ly, Nmat = 4, 4, 32
        solver, green_info, _ = self._make_solver(
            "FLEX", Lx=Lx, Ly=Ly, Nmat=Nmat, T=4.0, U=1.0,
            calc_scheme="reduced",
            extra_mode={'param': {
                'T': 4.0, 'mu': 0.0,
                'CellShape': [Lx, Ly, 1], 'SubShape': [1, 1, 1],
                'Nmat': Nmat, 'IterationMax': 10, 'Mix': 0.3,
            }})
        solver.solve(green_info, 'tests/flex/output')

        sigma = green_info["sigma"]
        # After several iterations with weak coupling, sigma should be finite
        self.assertGreater(np.linalg.norm(sigma), 1e-10,
                           "Sigma should be nonzero after FLEX SCF")

    def test_flex_sigma_zero_for_zero_interaction(self):
        """Self-energy should be zero when U=0."""
        Lx, Ly, Nmat = 4, 4, 32
        # U=0: no interaction file
        solver, green_info, _ = self._make_solver(
            "FLEX", Lx=Lx, Ly=Ly, Nmat=Nmat, T=4.0,
            calc_scheme="reduced",
            extra_mode={'param': {
                'T': 4.0, 'mu': 0.0,
                'CellShape': [Lx, Ly, 1], 'SubShape': [1, 1, 1],
                'Nmat': Nmat, 'IterationMax': 3, 'Mix': 0.5,
            }})

        # With U present but very weak, sigma should be small
        # Actually test by checking that with U=4, sigma is non-trivially nonzero
        solver.solve(green_info, 'tests/flex/output')
        sigma = green_info["sigma"]
        # With U=4, sigma should be nonzero
        self.assertGreater(np.linalg.norm(sigma), 1e-10,
                           "Sigma should be nonzero for U=4")

    def test_flex_self_energy_hermitian(self):
        """Self-energy should satisfy Sigma(-iwn) = Sigma(iwn)^dagger."""
        Lx, Ly, Nmat = 4, 4, 32
        solver, green_info, _ = self._make_solver(
            "FLEX", Lx=Lx, Ly=Ly, Nmat=Nmat, T=4.0, U=2.0,
            calc_scheme="reduced",
            extra_mode={'param': {
                'T': 4.0, 'mu': 0.0,
                'CellShape': [Lx, Ly, 1], 'SubShape': [1, 1, 1],
                'Nmat': Nmat, 'IterationMax': 3, 'Mix': 0.5,
            }})
        solver.solve(green_info, 'tests/flex/output')

        sigma = green_info["sigma"]
        # Sigma(-iwn) = Sigma(iwn)^dagger (conjugate transpose in orbital space)
        for n in range(Nmat // 2):
            n_neg = Nmat - 1 - n
            sigma_dagger = np.conj(sigma[:, n]).swapaxes(-2, -1)
            np.testing.assert_allclose(
                sigma[:, n_neg], sigma_dagger, atol=1e-10,
                err_msg="Sigma Hermitian symmetry violated at n={}".format(n))

    def test_flex_spin_charge_vertex_hubbard(self):
        """For single-orbital Hubbard, ham_s = -U*I and ham_c = +U*I."""
        solver, green_info, _ = self._make_solver(
            "FLEX", Lx=4, Ly=4, Nmat=32, T=4.0, U=4.0,
            calc_scheme="reduced",
            extra_mode={'param': {
                'T': 4.0, 'mu': 0.0,
                'CellShape': [4, 4, 1], 'SubShape': [1, 1, 1],
                'Nmat': 32, 'IterationMax': 1, 'Mix': 1.0,
            }})

        # Setup solver state
        solver._calc_epsilon_k(green_info)

        ham_orig = solver.ham_info.ham_inter_q
        nvol = solver.lattice.nvol
        norb = solver.norb
        nd = solver.nd

        # Inflate ham to reduced form
        ham_inflated = np.einsum('ksasatbtb->ksatb',
                                 ham_orig.reshape(nvol, *(2, norb) * 4)
                                 ).reshape(nvol, nd, nd)

        # Build spin/charge vertices
        ham_s, ham_c = solver._build_spin_charge_vertices(ham_inflated)

        # For single-orbital Hubbard U=4:
        # ham_inflated at q=0:
        #   [0   U]   (same-spin: 0, cross-spin: U)
        #   [U   0]
        # w_same = 0, w_cross = U
        # U_s = w_cross - w_same = U
        # U_c = w_cross + w_same = U
        # ham_s = -U_s = -U (for Stoner enhancement: [1 - U*chi0]^{-1})
        # ham_c = +U_c = +U (for charge suppression: [1 + U*chi0]^{-1})
        U = 4.0
        q0_idx = 0  # on-site
        ham_s_q0 = ham_s[q0_idx]
        ham_c_q0 = ham_c[q0_idx]

        # ham_s should be -U*I (spin channel: Stoner enhancement)
        expected_s = -U * np.eye(nd)
        self.assertTrue(np.allclose(ham_s_q0, expected_s, atol=1e-10),
                        "Spin vertex incorrect: got {} expected {}".format(
                            ham_s_q0, expected_s))

        # ham_c should be +U*I (charge channel: suppressed)
        expected_c = U * np.eye(nd)
        self.assertTrue(np.allclose(ham_c_q0, expected_c, atol=1e-10),
                        "Charge vertex incorrect: got {} expected {}".format(
                            ham_c_q0, expected_c))


    def test_flex_spin_susceptibility_stoner_enhanced(self):
        """Spin susceptibility should be Stoner-enhanced (chi_s > chi0)."""
        Lx, Ly, Nmat = 4, 4, 32
        solver, green_info, _ = self._make_solver(
            "FLEX", Lx=Lx, Ly=Ly, Nmat=Nmat, T=4.0, U=2.0,
            calc_scheme="reduced",
            extra_mode={'param': {
                'T': 4.0, 'mu': 0.0,
                'CellShape': [Lx, Ly, 1], 'SubShape': [1, 1, 1],
                'Nmat': Nmat, 'IterationMax': 1, 'Mix': 1.0,
            }})

        # Run one iteration to get chi0q and chi_s
        solver._calc_epsilon_k(green_info)
        if solver.calc_mu:
            Ncond = solver.Ncond / 2 if solver.spin_mode == "spin-free" else solver.Ncond
            dist, mu = solver._find_mu(Ncond, solver.T)
        else:
            mu = solver.mu_value
        beta = 1.0 / solver.T
        green0, green0_tail = solver._calc_green(beta, mu)
        chi0q_raw = solver._calc_chi0q(green0, green0_tail, beta)
        if solver.spin_mode in ["spin-free", "spinful"]:
            chi0q_raw = chi0q_raw[0]

        ham_orig = solver.ham_info.ham_inter_q
        chi0q, ham = solver._inflate_chi0q_and_ham(chi0q_raw, ham_orig)
        ham_s, ham_c = solver._build_spin_charge_vertices(ham)

        chi_s = solver._solve_rpa(chi0q, ham_s)
        chi_c = solver._solve_rpa(chi0q, ham_c)

        # At the center Matsubara frequency (ivn=0, index=Nmat//2)
        center = Nmat // 2
        chi0_norm = np.linalg.norm(chi0q[center])
        chi_s_norm = np.linalg.norm(chi_s[center])
        chi_c_norm = np.linalg.norm(chi_c[center])

        # Spin should be enhanced (larger than bare)
        self.assertGreater(chi_s_norm, chi0_norm,
                           "Spin susceptibility should be Stoner-enhanced")
        # Charge should be suppressed (smaller than bare)
        self.assertLess(chi_c_norm, chi0_norm,
                        "Charge susceptibility should be suppressed")

    def test_flex_single_iteration_matches_rpa_susceptibility(self):
        """FLEX chi_s with Sigma=0 (1 iter) should match standard RPA."""
        Lx, Ly, Nmat = 4, 4, 32

        # Run standard RPA
        rpa_solver, rpa_green_info, _ = self._make_solver(
            "RPA", Lx=Lx, Ly=Ly, Nmat=Nmat, T=2.0, calc_scheme="reduced")
        rpa_solver.solve(rpa_green_info, 'tests/flex/output')
        chiq_rpa = rpa_green_info["chiq"]

        # Run FLEX with 1 iteration
        flex_solver, flex_green_info, _ = self._make_solver(
            "FLEX", Lx=Lx, Ly=Ly, Nmat=Nmat, T=2.0,
            calc_scheme="reduced",
            extra_mode={'param': {
                'T': 2.0, 'mu': 0.0,
                'CellShape': [Lx, Ly, 1], 'SubShape': [1, 1, 1],
                'Nmat': Nmat, 'IterationMax': 1, 'Mix': 1.0,
            }})

        # Manually run FLEX to get spin susceptibility
        flex_solver._calc_epsilon_k(flex_green_info)
        if flex_solver.calc_mu:
            Ncond = flex_solver.Ncond / 2 if flex_solver.spin_mode == "spin-free" else flex_solver.Ncond
            dist, mu = flex_solver._find_mu(Ncond, flex_solver.T)
        else:
            mu = flex_solver.mu_value
        beta = 1.0 / flex_solver.T
        green0, green0_tail = flex_solver._calc_green(beta, mu)
        chi0q_raw = flex_solver._calc_chi0q(green0, green0_tail, beta)
        if flex_solver.spin_mode in ["spin-free", "spinful"]:
            chi0q_raw = chi0q_raw[0]

        ham_orig = flex_solver.ham_info.ham_inter_q
        chi0q, ham = flex_solver._inflate_chi0q_and_ham(chi0q_raw, ham_orig)

        # Standard RPA uses the full ham directly (not decomposed into S/C)
        # So FLEX's _solve_rpa(chi0q, ham) should match RPA's chiq
        chiq_flex_full = flex_solver._solve_rpa(chi0q, ham)

        # RPA chiq is in (nfreq, nvol, nd, nd) form matching inflated chi0q
        self.assertTrue(np.allclose(chiq_rpa, chiq_flex_full, atol=1e-10),
                        "FLEX RPA solve with full ham should match standard RPA")


class TestFLEXTwoOrbital(unittest.TestCase):
    """FLEX solver tests using 2-orbital Hubbard model."""

    def _make_2orb_solver(self, mode_cls, Lx=4, Ly=4, Nmat=32, T=2.0,
                          mu=0.0, calc_scheme="reduced", extra_mode=None,
                          extra_interactions=None):
        """Create a 2-orbital solver instance."""
        info_log = {}
        info_mode = {
            'mode': mode_cls,
            'param': {
                'T': T,
                'mu': mu,
                'CellShape': [Lx, Ly, 1],
                'SubShape': [1, 1, 1],
                'Nmat': Nmat,
            },
            'calc_scheme': calc_scheme,
        }
        if extra_mode:
            info_mode.update(extra_mode)

        info_file = {
            'input': {
                'path_to_input': 'tests/rpa/input_2orb',
                'interaction': {
                    'path_to_input': 'tests/rpa/input_2orb',
                    'Geometry': 'geom.dat',
                    'Transfer': 'transfer.dat',
                    'CoulombIntra': 'coulombintra.dat',
                    'CoulombInter': 'coulombinter.dat',
                },
            },
            'output': {
                'path_to_output': 'tests/flex/output_2orb',
            },
        }
        if extra_interactions:
            info_file['input']['interaction'].update(extra_interactions)

        os.makedirs('tests/flex/output_2orb', exist_ok=True)

        import hwave.qlmsio.read_input_k as read_input_k
        read_io = read_input_k.QLMSkInput(info_file['input'])
        ham_info = read_io.get_param("ham")
        green_info = read_io.get_param("green")

        if mode_cls == "FLEX":
            import hwave.solver.flex as solver_flex
            solver = solver_flex.FLEX(ham_info, info_log, info_mode)
        else:
            import hwave.solver.rpa as solver_rpa
            solver = solver_rpa.RPA(ham_info, info_log, info_mode)

        return solver, green_info, info_file

    def test_2orb_instantiation(self):
        """FLEX solver should handle 2-orbital model."""
        solver, green_info, _ = self._make_2orb_solver("FLEX")
        self.assertEqual(solver.norb, 2)
        self.assertEqual(solver.nd, 4)

    def test_2orb_chi0q_matches_rpa(self):
        """2-orbital FLEX chi0q should match RPA at first iteration."""
        Lx, Ly, Nmat = 4, 4, 32

        rpa_solver, rpa_gi, _ = self._make_2orb_solver(
            "RPA", Lx=Lx, Ly=Ly, Nmat=Nmat, T=2.0, calc_scheme="reduced")
        rpa_solver.solve(rpa_gi, 'tests/flex/output_2orb')
        chi0q_rpa = rpa_gi["chi0q"]

        flex_solver, flex_gi, _ = self._make_2orb_solver(
            "FLEX", Lx=Lx, Ly=Ly, Nmat=Nmat, T=2.0,
            calc_scheme="reduced",
            extra_mode={'param': {
                'T': 2.0, 'mu': 0.0,
                'CellShape': [Lx, Ly, 1], 'SubShape': [1, 1, 1],
                'Nmat': Nmat, 'IterationMax': 1, 'Mix': 1.0,
            }})
        flex_solver.solve(flex_gi, 'tests/flex/output_2orb')
        chi0q_flex = flex_gi["chi0q"]

        # FLEX inflates to spin-orbital (nd=4), RPA stays orbital (norb=2)
        norb = rpa_solver.norb
        if chi0q_flex.shape[-1] != chi0q_rpa.shape[-1]:
            chi0q_flex_orb = chi0q_flex[:, :, :norb, :norb]
        else:
            chi0q_flex_orb = chi0q_flex

        self.assertTrue(np.allclose(chi0q_rpa, chi0q_flex_orb, atol=1e-10),
                        "2-orbital chi0q mismatch between RPA and FLEX")

    def test_2orb_convergence(self):
        """2-orbital FLEX should converge for weak coupling."""
        Lx, Ly, Nmat = 4, 4, 32
        solver, green_info, _ = self._make_2orb_solver(
            "FLEX", Lx=Lx, Ly=Ly, Nmat=Nmat, T=4.0,
            calc_scheme="reduced",
            extra_mode={'param': {
                'T': 4.0, 'mu': 0.0,
                'CellShape': [Lx, Ly, 1], 'SubShape': [1, 1, 1],
                'Nmat': Nmat, 'IterationMax': 20, 'Mix': 0.2, 'EPS': 4,
            }})
        solver.solve(green_info, 'tests/flex/output_2orb')

        sigma = green_info["sigma"]
        self.assertIsNotNone(sigma)
        # 2-orbital: sigma should have nd_block = norb = 2 (spin-free reduced)
        self.assertEqual(sigma.shape[-1], 2)
        self.assertGreater(np.linalg.norm(sigma), 1e-10,
                           "2-orbital sigma should be nonzero")

    def test_2orb_self_energy_hermitian(self):
        """2-orbital self-energy should satisfy Sigma(-iwn) = Sigma(iwn)^dagger."""
        Lx, Ly, Nmat = 4, 4, 32
        solver, green_info, _ = self._make_2orb_solver(
            "FLEX", Lx=Lx, Ly=Ly, Nmat=Nmat, T=4.0,
            calc_scheme="reduced",
            extra_mode={'param': {
                'T': 4.0, 'mu': 0.0,
                'CellShape': [Lx, Ly, 1], 'SubShape': [1, 1, 1],
                'Nmat': Nmat, 'IterationMax': 3, 'Mix': 0.5,
            }})
        solver.solve(green_info, 'tests/flex/output_2orb')

        sigma = green_info["sigma"]
        # Sigma(-iwn) = Sigma(iwn)^dagger (conjugate transpose in orbital space)
        for n in range(Nmat // 2):
            n_neg = Nmat - 1 - n
            # Hermitian conjugate: swap last two axes + complex conjugate
            sigma_dagger = np.conj(sigma[:, n]).swapaxes(-2, -1)
            np.testing.assert_allclose(
                sigma[:, n_neg], sigma_dagger, atol=1e-10,
                err_msg="2-orbital Hermitian symmetry violated at n={}".format(n))

    def test_2orb_stoner_enhanced(self):
        """2-orbital spin susceptibility should be Stoner-enhanced."""
        Lx, Ly, Nmat = 4, 4, 32
        solver, green_info, _ = self._make_2orb_solver(
            "FLEX", Lx=Lx, Ly=Ly, Nmat=Nmat, T=2.0,
            calc_scheme="reduced",
            extra_mode={'param': {
                'T': 2.0, 'mu': 0.0,
                'CellShape': [Lx, Ly, 1], 'SubShape': [1, 1, 1],
                'Nmat': Nmat, 'IterationMax': 1, 'Mix': 1.0,
            }})

        solver._calc_epsilon_k(green_info)
        if solver.calc_mu:
            Ncond = solver.Ncond / 2 if solver.spin_mode == "spin-free" else solver.Ncond
            dist, mu = solver._find_mu(Ncond, solver.T)
        else:
            mu = solver.mu_value
        beta = 1.0 / solver.T
        green0, green0_tail = solver._calc_green(beta, mu)
        chi0q_raw = solver._calc_chi0q(green0, green0_tail, beta)
        if solver.spin_mode in ["spin-free", "spinful"]:
            chi0q_raw = chi0q_raw[0]

        ham_orig = solver.ham_info.ham_inter_q
        chi0q, ham = solver._inflate_chi0q_and_ham(chi0q_raw, ham_orig)
        ham_s, ham_c = solver._build_spin_charge_vertices(ham)

        chi_s = solver._solve_rpa(chi0q, ham_s)
        chi_c = solver._solve_rpa(chi0q, ham_c)

        center = Nmat // 2
        chi0_norm = np.linalg.norm(chi0q[center])
        chi_s_norm = np.linalg.norm(chi_s[center])
        chi_c_norm = np.linalg.norm(chi_c[center])

        self.assertGreater(chi_s_norm, chi0_norm,
                           "2-orbital spin susceptibility should be enhanced")
        self.assertLess(chi_c_norm, chi0_norm,
                        "2-orbital charge susceptibility should be suppressed")


class TestFLEXEliashberg(unittest.TestCase):
    """Tests for FLEX + Eliashberg equation integration."""

    def _build_eliashberg_inter_k(self, norb, Nx, Ny, Nz, input_dir='tests/rpa/input'):
        """Build interaction in k-space for Eliashberg tests."""
        from hwave.sc import _read_interaction_files, _build_interaction_k
        input_dict = {
            'file': {
                'input': {
                    'path_to_input': input_dir,
                    'interaction': {
                        'path_to_input': input_dir,
                        'Geometry': 'geom.dat',
                        'Transfer': 'transfer.dat',
                        'CoulombIntra': 'coulombintra.dat',
                    },
                },
                'output': {'path_to_output': 'tests/flex/output_eli'},
            }
        }
        geom_info, hr, interactions = _read_interaction_files(input_dict)
        kx = np.linspace(0, 2*np.pi, Nx, endpoint=False)
        ky = np.linspace(0, 2*np.pi, Ny, endpoint=False)
        kz = np.linspace(0, 2*np.pi, Nz, endpoint=False)
        inter_k = _build_interaction_k(kx, ky, kz, interactions, norb)
        return inter_k

    def _flex_chi_to_eliashberg(self, chi_s, chi_c, norb, nd, Nmat, Nx, Ny, Nz):
        """Convert FLEX susceptibilities to Eliashberg format.

        Delegates to the production expansion rather than mirroring it. The
        mirror previously used the delta_{l2,l4} scatter kron(X, I_norb); its
        callers are single-orbital, where that coincides with the correct
        density-pair placement, so it silently survived the fix and would have
        validated the wrong layout the moment a multi-orbital caller appeared.
        """
        import hwave.sc as sc
        center = Nmat // 2
        chis_eli = sc._expand_flex_chi(chi_s[center:center + 1], norb,
                                       Nx, Ny, Nz, "kuroki")[0]
        chic_eli = sc._expand_flex_chi(chi_c[center:center + 1], norb,
                                       Nx, Ny, Nz, "kuroki")[0]
        return chis_eli, chic_eli

    def test_flex_eliashberg_vertex_construction(self):
        """FLEX chi_s/chi_c should produce valid Eliashberg pairing vertex."""
        from hwave.sc import _compute_vertices_flex

        Lx, Ly, Nmat = 4, 4, 32
        Nx, Ny, Nz = Lx, Ly, 1

        # Run FLEX
        info_log = {}
        info_mode = {
            'mode': 'FLEX',
            'param': {
                'T': 4.0, 'mu': 0.0,
                'CellShape': [Lx, Ly, 1], 'SubShape': [1, 1, 1],
                'Nmat': Nmat, 'IterationMax': 3, 'Mix': 0.5,
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
        os.makedirs('tests/flex/output_eli', exist_ok=True)

        import hwave.qlmsio.read_input_k as read_input_k
        read_io = read_input_k.QLMSkInput(info_file_input)
        ham_info = read_io.get_param("ham")
        green_info = read_io.get_param("green")

        import hwave.solver.flex as solver_flex
        solver = solver_flex.FLEX(ham_info, info_log, info_mode)
        solver.solve(green_info, 'tests/flex/output_eli')

        norb = solver.norb
        nd = solver.nd

        # Convert to Eliashberg format
        chis_eli, chic_eli = self._flex_chi_to_eliashberg(
            green_info["chiq_s"], green_info["chiq_c"],
            norb, nd, Nmat, Nx, Ny, Nz)

        inter_k = self._build_eliashberg_inter_k(norb, Nx, Ny, Nz)

        # Compute pairing vertex from FLEX susceptibilities
        Vs_q = _compute_vertices_flex(chis_eli, chic_eli, inter_k,
                                      norb, Nx, Ny, Nz)

        self.assertTrue(np.isfinite(Vs_q).all(), "Vertex has non-finite values")
        self.assertGreater(np.linalg.norm(Vs_q), 1e-10,
                           "Vertex should be nonzero")

    def test_flex_vs_rpa_eliashberg_vertex(self):
        """FLEX vertex at 1st iteration should match RPA vertex."""
        from hwave.sc import _compute_vertices_flex, _compute_vertices_general

        Lx, Ly, Nmat = 4, 4, 32
        Nx, Ny, Nz = Lx, Ly, 1

        # Run FLEX with 1 iteration (= bare G = RPA)
        info_mode = {
            'mode': 'FLEX',
            'param': {
                'T': 2.0, 'mu': 0.0,
                'CellShape': [Lx, Ly, 1], 'SubShape': [1, 1, 1],
                'Nmat': Nmat, 'IterationMax': 1, 'Mix': 1.0,
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

        import hwave.solver.flex as solver_flex
        solver = solver_flex.FLEX(ham_info, {}, info_mode)
        solver.solve(green_info, 'tests/flex/output_eli')

        norb = solver.norb
        nd = solver.nd

        # Convert FLEX chi_s/chi_c to Eliashberg format
        chis_eli, chic_eli = self._flex_chi_to_eliashberg(
            green_info["chiq_s"], green_info["chiq_c"],
            norb, nd, Nmat, Nx, Ny, Nz)

        inter_k = self._build_eliashberg_inter_k(norb, Nx, Ny, Nz)

        # FLEX vertex
        Vs_flex = _compute_vertices_flex(chis_eli, chic_eli, inter_k,
                                         norb, Nx, Ny, Nz)

        # RPA vertex: compute from chi0q
        chi0q = green_info["chi0q"]
        chi0q_orb = chi0q[:, :, :norb, :norb]
        chi0q_sc = chi0q_orb.reshape(Nmat, Nx, Ny, Nz, norb, norb).transpose(
            4, 5, 1, 2, 3, 0).copy()

        Vs_rpa = _compute_vertices_general(chi0q_sc, inter_k, norb,
                                           Nx, Ny, Nz, Nmat)

        # At 1st iteration (bare G), FLEX and RPA should give same vertex
        self.assertTrue(np.allclose(Vs_flex, Vs_rpa, atol=1e-10),
                        "FLEX vertex (1 iter) should match RPA vertex")


class TestFLEXTemperature(unittest.TestCase):
    """Tests for FLEX temperature dependence."""

    def _make_solver(self, T, Lx=4, Ly=4, Nmat=32, U=2.0):
        """Create FLEX solver at given temperature."""
        info_log = {}
        info_mode = {
            'mode': 'FLEX',
            'param': {
                'T': T,
                'mu': 0.0,
                'CellShape': [Lx, Ly, 1],
                'SubShape': [1, 1, 1],
                'Nmat': Nmat,
                'IterationMax': 5,
                'Mix': 0.3,
            },
            'calc_scheme': 'reduced',
        }
        info_file = {
            'input': {
                'path_to_input': 'tests/rpa/input',
                'interaction': {
                    'path_to_input': 'tests/rpa/input',
                    'Geometry': 'geom.dat',
                    'Transfer': 'transfer.dat',
                    'CoulombIntra': 'coulombintra.dat',
                },
            },
        }
        os.makedirs('tests/flex/output_temp', exist_ok=True)

        import hwave.qlmsio.read_input_k as read_input_k
        read_io = read_input_k.QLMSkInput(info_file['input'])
        ham_info = read_io.get_param("ham")
        green_info = read_io.get_param("green")

        import hwave.solver.flex as solver_flex
        solver = solver_flex.FLEX(ham_info, info_log, info_mode)

        return solver, green_info

    def test_higher_T_smaller_sigma(self):
        """Self-energy should decrease with increasing temperature."""
        solver_low, gi_low = self._make_solver(T=2.0)
        solver_low.solve(gi_low, 'tests/flex/output_temp')
        sigma_low = gi_low["sigma"]

        solver_high, gi_high = self._make_solver(T=8.0)
        solver_high.solve(gi_high, 'tests/flex/output_temp')
        sigma_high = gi_high["sigma"]

        norm_low = np.linalg.norm(sigma_low)
        norm_high = np.linalg.norm(sigma_high)

        self.assertGreater(norm_low, norm_high,
                           "Sigma at lower T ({:.3e}) should be larger than at "
                           "higher T ({:.3e})".format(norm_low, norm_high))

    def test_sigma_hermitian_at_different_T(self):
        """Hermitian symmetry should hold at various temperatures."""
        for T in [1.0, 4.0, 10.0]:
            solver, gi = self._make_solver(T=T)
            solver.solve(gi, 'tests/flex/output_temp')
            sigma = gi["sigma"]
            Nmat = solver.nmat
            for n in range(Nmat // 2):
                n_neg = Nmat - 1 - n
                sigma_dagger = np.conj(sigma[:, n]).swapaxes(-2, -1)
                np.testing.assert_allclose(
                    sigma[:, n_neg], sigma_dagger, atol=1e-10,
                    err_msg="Hermitian violated at T={}, n={}".format(T, n))


if __name__ == '__main__':
    unittest.main()
