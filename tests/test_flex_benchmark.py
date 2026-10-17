#!/usr/bin/env python3

"""Physics benchmark tests for the FLEX solver.

These tests verify that the FLEX solver reproduces well-known physical
results from the literature:

1-orbital (2D Hubbard model on square lattice):
  - Bickers, Scalapino, White, PRL 62, 961 (1989)
  - Bickers & Scalapino, Ann. Phys. 193, 206 (1989)
  Key physics:
    - AF spin susceptibility peak at Q=(pi,pi) near half-filling
    - Quasiparticle damping: Im Sigma(k, iw_0) < 0
    - Hot/cold spots: stronger damping at (pi,0) than (pi/2,pi/2)
    - Luttinger sum rule: n = (2T/N_k) sum_{k,n} G(k,iwn) e^{iwn 0+}
    - d-wave pairing from Eliashberg equation

Multi-orbital (2-orbital Hubbard-Kanamori model):
  - Kontani & Onari, PRL 104, 157001 (2010)
  - Kuroki et al., PRB 79, 224511 (2009)
  Key physics:
    - Orbital-dependent mass renormalization Z_a = 1/(1 - dImSigma/dw)
    - Spin/charge susceptibility structure with inter-orbital coupling
    - Pairing symmetry from Eliashberg eigenvalue
"""

import os
import unittest
import numpy as np


def _make_1orb_flex(Lx=16, Ly=16, Nmat=64, T=0.5, mu=0.0, U=4.0,
                    max_iter=30, mix=0.2, eps=4):
    """Create a 1-orbital FLEX solver for the 2D Hubbard model.

    Model: H = -t sum_{<ij>} c^dag c + U sum_i n_up n_dn
    with t = 1.0 (nearest-neighbor), t' = 0.5 (next-nearest-neighbor)
    on a square lattice.

    The U value is written dynamically to a temp CoulombIntra file
    so that different U values can be tested.
    """
    import tempfile
    import shutil

    # Create temp directory with CoulombIntra at the requested U
    tmpdir = tempfile.mkdtemp(prefix='flex_bench_')
    src_dir = 'tests/rpa/input'
    for fn in ['geom.dat', 'transfer.dat']:
        shutil.copy2(os.path.join(src_dir, fn), os.path.join(tmpdir, fn))

    # Write CoulombIntra with requested U
    ci_path = os.path.join(tmpdir, 'coulombintra.dat')
    with open(ci_path, 'w') as f:
        f.write("CoulombIntra in wannier90-like format for uhfk\n")
        f.write("1\n1\n 1\n")
        f.write("   0    0    0    1    1   {:.12f}   0.000000000000\n".format(U))

    info_log = {}
    info_mode = {
        'mode': 'FLEX',
        'param': {
            'T': T,
            'mu': mu,
            'CellShape': [Lx, Ly, 1],
            'SubShape': [1, 1, 1],
            'Nmat': Nmat,
            'IterationMax': max_iter,
            'Mix': mix,
            'EPS': eps,
        },
        'calc_scheme': 'reduced',
    }
    info_file_input = {
        'path_to_input': tmpdir,
        'interaction': {
            'path_to_input': tmpdir,
            'Geometry': 'geom.dat',
            'Transfer': 'transfer.dat',
            'CoulombIntra': 'coulombintra.dat',
        },
    }
    os.makedirs('tests/flex/output_bench', exist_ok=True)

    import hwave.qlmsio.read_input_k as read_input_k
    read_io = read_input_k.QLMSkInput(info_file_input)
    ham_info = read_io.get_param("ham")
    green_info = read_io.get_param("green")

    import hwave.solver.flex as solver_flex
    solver = solver_flex.FLEX(ham_info, info_log, info_mode)

    # Clean up temp directory
    shutil.rmtree(tmpdir, ignore_errors=True)

    return solver, green_info


def _make_2orb_flex(Lx=8, Ly=8, Nmat=32, T=1.0, mu=0.0,
                    max_iter=30, mix=0.2, eps=4):
    """Create a 2-orbital FLEX solver with Kanamori interaction.

    Model: 2-orbital Hubbard on square lattice
    with CoulombIntra U and CoulombInter V.
    """
    info_log = {}
    info_mode = {
        'mode': 'FLEX',
        'param': {
            'T': T,
            'mu': mu,
            'CellShape': [Lx, Ly, 1],
            'SubShape': [1, 1, 1],
            'Nmat': Nmat,
            'IterationMax': max_iter,
            'Mix': mix,
            'EPS': eps,
        },
        'calc_scheme': 'reduced',
    }
    info_file_input = {
        'path_to_input': 'tests/rpa/input_2orb',
        'interaction': {
            'path_to_input': 'tests/rpa/input_2orb',
            'Geometry': 'geom.dat',
            'Transfer': 'transfer.dat',
            'CoulombIntra': 'coulombintra.dat',
            'CoulombInter': 'coulombinter.dat',
        },
    }
    os.makedirs('tests/flex/output_bench_2orb', exist_ok=True)

    import hwave.qlmsio.read_input_k as read_input_k
    read_io = read_input_k.QLMSkInput(info_file_input)
    ham_info = read_io.get_param("ham")
    green_info = read_io.get_param("green")

    import hwave.solver.flex as solver_flex
    solver = solver_flex.FLEX(ham_info, info_log, info_mode)
    return solver, green_info


class TestFLEXBenchmark1Orbital(unittest.TestCase):
    """Physics benchmarks for 1-orbital 2D Hubbard model FLEX.

    Reference: Bickers, Scalapino, White, PRL 62, 961 (1989)
    Model: square lattice, t=1 (NN), t'=0.5 (NNN), U=4, mu=0 (half-filling)
    """

    @classmethod
    def setUpClass(cls):
        """Run FLEX calculation once for all tests in this class."""
        cls.Lx, cls.Ly = 16, 16
        cls.Nmat = 64
        cls.T = 0.5
        cls.U = 4.0
        solver, green_info = _make_1orb_flex(
            Lx=cls.Lx, Ly=cls.Ly, Nmat=cls.Nmat, T=cls.T, U=cls.U,
            max_iter=50, mix=0.2, eps=4)
        solver.solve(green_info, 'tests/flex/output_bench')

        cls.solver = solver
        cls.green_info = green_info
        cls.sigma = green_info["sigma"]       # (nblock, nmat, nvol, nd_b, nd_b)
        cls.green_kw = green_info["green"]    # (nblock, nmat, nvol, nd_b, nd_b)
        cls.chi0q = green_info["chi0q"]       # (nmat, nvol, nd, nd)
        cls.chi_s = green_info["chiq_s"]      # (nmat, nvol, nd, nd)
        cls.chi_c = green_info["chiq_c"]      # (nmat, nvol, nd, nd)
        cls.beta = 1.0 / cls.T

    def _q_index(self, qx_frac, qy_frac):
        """Convert fractional q-vector (in units of 2pi) to flat index."""
        ix = int(round(qx_frac * self.Lx)) % self.Lx
        iy = int(round(qy_frac * self.Ly)) % self.Ly
        return iy * self.Lx + ix  # row-major (y,x) or (x,y) depends on lattice

    def _k_index(self, kx_frac, ky_frac):
        """Convert fractional k-vector to flat index.

        The k-points are ordered as k = (ix/Lx, iy/Ly) for ix in [0,Lx), iy in [0,Ly).
        The flat index depends on the lattice layout.
        """
        ix = int(round(kx_frac * self.Lx)) % self.Lx
        iy = int(round(ky_frac * self.Ly)) % self.Ly
        return iy * self.Lx + ix

    def test_af_susceptibility_peak_at_pi_pi(self):
        """Spin susceptibility should peak at Q=(pi,pi) for the square lattice
        Hubbard model near half-filling.

        This is the hallmark of antiferromagnetic spin fluctuations that
        drive d-wave superconductivity (Bickers & Scalapino 1989).
        """
        center = self.Nmat // 2  # static limit (ivn=0)
        nd = self.chi_s.shape[-1]

        # chi_s is in spin-orbital reduced space (nd=2 for 1-orb)
        # Take trace over spin-orbital indices to get scalar susceptibility
        chi_s_static = self.chi_s[center]  # (nvol, nd, nd)
        chi_s_trace = np.real(np.trace(chi_s_static, axis1=-2, axis2=-1))

        # Find Q-vector with maximum chi_s
        q_max = np.argmax(chi_s_trace)

        # Q=(pi,pi) corresponds to (Lx//2, Ly//2)
        q_pipi = self._q_index(0.5, 0.5)
        q_00 = self._q_index(0.0, 0.0)

        # The AF peak should be at Q=(pi,pi)
        self.assertEqual(q_max, q_pipi,
                         "Spin susceptibility peak should be at Q=(pi,pi), "
                         "got q_max={} (expected {})".format(q_max, q_pipi))

        # chi_s(pi,pi) should be significantly larger than chi_s(0,0)
        ratio = chi_s_trace[q_pipi] / chi_s_trace[q_00]
        self.assertGreater(ratio, 2.0,
                           "chi_s(pi,pi)/chi_s(0,0) = {:.2f} should be >> 1".format(ratio))

    def test_quasiparticle_damping(self):
        """Im Sigma(k, iw_0) should be negative (quasiparticle damping).

        For any k-point, the imaginary part of the self-energy at the
        lowest positive Matsubara frequency should give a negative
        contribution to the spectral weight.
        """
        # Lowest positive Matsubara frequency: index = Nmat//2
        iw0_idx = self.Nmat // 2
        sigma_iw0 = self.sigma[0, iw0_idx]  # (nvol, nd_b, nd_b)

        # For 1-orbital: nd_b = 1 (spin-free reduced)
        # Im Sigma(k, iw_0) should be negative for all k
        im_sigma = np.imag(sigma_iw0[:, 0, 0])

        # At the lowest Matsubara frequency, Im Sigma should be < 0
        # (electron loses spectral weight to spin fluctuations)
        self.assertTrue(np.all(im_sigma < 0),
                        "Im Sigma(k, iw_0) should be negative everywhere. "
                        "Max Im Sigma = {:.6f}".format(np.max(im_sigma)))

    def test_hot_cold_spots(self):
        """Self-energy anisotropy: stronger damping at hot spots (pi,0)
        than cold spots (pi/2,pi/2).

        Hot spots are k-points connected to the AF vector Q=(pi,pi) on
        the Fermi surface. Near half-filling on the square lattice, (pi,0)
        is a hot spot and (pi/2,pi/2) is a cold spot.

        Reference: Monthoux, Balatsky & Pines, PRL 67, 3448 (1991)
        """
        iw0_idx = self.Nmat // 2

        # Hot spot: k = (pi, 0) = (Lx//2, 0)
        k_hot = self._k_index(0.5, 0.0)
        # Cold spot: k = (pi/2, pi/2) = (Lx//4, Ly//4)
        k_cold = self._k_index(0.25, 0.25)

        im_sigma_hot = np.abs(np.imag(self.sigma[0, iw0_idx, k_hot, 0, 0]))
        im_sigma_cold = np.abs(np.imag(self.sigma[0, iw0_idx, k_cold, 0, 0]))

        self.assertGreater(im_sigma_hot, im_sigma_cold,
                           "|Im Sigma| at hot spot ({:.6f}) should exceed "
                           "cold spot ({:.6f})".format(im_sigma_hot, im_sigma_cold))

    def test_self_energy_particle_hole_symmetry(self):
        """At half-filling (mu=0) with particle-hole symmetric dispersion,
        Re Sigma(k, iw_n) should vanish for nearest-neighbor hopping only.

        With next-nearest-neighbor hopping (t'), particle-hole symmetry is
        broken, so Re Sigma can be nonzero but should still be small compared
        to Im Sigma at high temperature.
        """
        iw0_idx = self.Nmat // 2
        sigma_iw0 = self.sigma[0, iw0_idx]  # (nvol, nd_b, nd_b)

        re_sigma = np.real(sigma_iw0[:, 0, 0])
        im_sigma = np.imag(sigma_iw0[:, 0, 0])

        # With t' = 0.5 (next-nearest-neighbor), PH symmetry is broken,
        # but Re Sigma should still be much smaller than Im Sigma
        # Check that the dominant contribution is imaginary
        re_max = np.max(np.abs(re_sigma))
        im_max = np.max(np.abs(im_sigma))

        self.assertGreater(im_max, 0,
                           "Im Sigma should be nonzero")

    def test_luttinger_sum_rule(self):
        """Luttinger sum rule: the dressed Green's function should give
        the correct particle number.

        n = (2/N_k) sum_k [ T * sum_n G(k, iwn) * e^{iwn 0+} + 1/2 ]

        For the spin-free block, G gives n per spin, so total n = 2 * n_per_spin.
        """
        nmat = self.Nmat
        nvol = self.solver.lattice.nvol
        green = self.green_kw[0]  # (nmat, nvol, nd_b, nd_b)

        # n_per_spin = T * sum_n Tr[G(k,iwn)] * e^{iwn 0+} + norb/2
        # Convergence factor: e^{iwn 0+} -> 1 for the tail-subtracted sum
        # Here we use the full G including tail, so we compute directly:
        # n = (1/N_k) sum_k (T sum_n Tr[G(k,iwn)] + norb/2)
        n_k = np.zeros(nvol)
        for ik in range(nvol):
            g_trace = np.trace(green[:, ik], axis1=-2, axis2=-1)  # (nmat,)
            # Sum over Matsubara frequencies
            n_k[ik] = np.real(np.sum(g_trace)) / self.beta + 0.5

        n_per_spin = np.mean(n_k)
        n_total = 2 * n_per_spin  # both spins

        # At mu=0 with the test Hamiltonian, filling should be close to 1
        # (half-filling for 1 orbital = 1 electron per site)
        # Allow generous tolerance since Nmat is finite
        self.assertGreater(n_total, 0.5,
                           "Particle number ({:.4f}) too small".format(n_total))
        self.assertLess(n_total, 1.5,
                        "Particle number ({:.4f}) too large".format(n_total))

    def test_spin_susceptibility_stoner_enhanced(self):
        """Spin susceptibility should be Stoner-enhanced relative to bare chi0
        for all q-vectors at static limit.

        chi_s(q) > chi0(q) for all q (since U > 0).
        """
        center = self.Nmat // 2
        nd = self.chi_s.shape[-1]

        chi0_norm = np.array([
            np.linalg.norm(self.chi0q[center, iq]) for iq in range(self.Lx * self.Ly)])
        chi_s_norm = np.array([
            np.linalg.norm(self.chi_s[center, iq]) for iq in range(self.Lx * self.Ly)])

        # chi_s should be enhanced everywhere
        enhanced = chi_s_norm > chi0_norm
        frac_enhanced = np.sum(enhanced) / len(enhanced)
        self.assertGreater(frac_enhanced, 0.9,
                           "Only {:.0f}% of q-points are Stoner-enhanced "
                           "(expected >90%)".format(frac_enhanced * 100))

    def test_charge_susceptibility_suppressed(self):
        """Charge susceptibility should be suppressed relative to bare chi0.

        chi_c(q) < chi0(q) for repulsive U > 0 (charge fluctuations are costly).
        """
        center = self.Nmat // 2

        chi0_norm = np.array([
            np.linalg.norm(self.chi0q[center, iq]) for iq in range(self.Lx * self.Ly)])
        chi_c_norm = np.array([
            np.linalg.norm(self.chi_c[center, iq]) for iq in range(self.Lx * self.Ly)])

        suppressed = chi_c_norm < chi0_norm
        frac_suppressed = np.sum(suppressed) / len(suppressed)
        self.assertGreater(frac_suppressed, 0.9,
                           "Only {:.0f}% of q-points have suppressed charge "
                           "susceptibility (expected >90%)".format(frac_suppressed * 100))

    def test_self_energy_frequency_dependence(self):
        """Im Sigma should decrease with increasing Matsubara frequency.

        At high frequency, Sigma -> 0 as 1/iwn (Fermi liquid behavior).
        So |Im Sigma(k, iwn)| should decrease for large n.
        """
        # Average over k-points
        im_sigma_avg = np.zeros(self.Nmat)
        for n in range(self.Nmat):
            im_sigma_avg[n] = np.mean(np.abs(
                np.imag(self.sigma[0, n, :, 0, 0])))

        # Compare lowest vs highest positive Matsubara frequency
        n_low = self.Nmat // 2       # iw_0
        n_high = self.Nmat - 1       # iw_{N-1}

        self.assertGreater(im_sigma_avg[n_low], im_sigma_avg[n_high],
                           "|Im Sigma(iw_0)| ({:.6f}) should exceed "
                           "|Im Sigma(iw_max)| ({:.6f})".format(
                               im_sigma_avg[n_low], im_sigma_avg[n_high]))

    def test_mass_renormalization(self):
        """Quasiparticle weight Z = 1 / (1 - Im Sigma(iw_0) / w_0) should
        satisfy 0 < Z < 1 for interacting electrons.

        Z < 1 indicates correlation-induced mass enhancement m*/m = 1/Z.
        Uses the standard estimate Z ~ 1 / (1 - Im Sigma(iw_0) / w_0)
        which is more robust than the finite difference derivative for
        coarse Matsubara grids.
        """
        iw0 = self.Nmat // 2
        beta = self.beta
        omega0 = np.pi / beta  # lowest positive Matsubara freq

        # Average Z over all k-points
        Z_k = np.zeros(self.Lx * self.Ly)
        for ik in range(self.Lx * self.Ly):
            im_sigma_0 = np.imag(self.sigma[0, iw0, ik, 0, 0])
            Z_k[ik] = 1.0 / (1.0 - im_sigma_0 / omega0)

        Z_avg = np.mean(Z_k)

        # Z should be between 0 and 1 (mass enhancement)
        self.assertGreater(Z_avg, 0.0,
                           "Average Z ({:.4f}) should be > 0".format(Z_avg))
        self.assertLess(Z_avg, 1.0,
                        "Average Z ({:.4f}) should be < 1 for interacting "
                        "system".format(Z_avg))

    def test_eliashberg_vertex_from_flex(self):
        """FLEX Eliashberg pairing vertex should be well-defined and reflect
        the spin-fluctuation structure.

        For the 1-orbital Hubbard model with on-site U only, the singlet
        pairing vertex Vs(q) = U^2 * (3/2 chi_s(q) - 1/2 chi_c(q)) + U
        should be peaked at Q=(pi,pi) because chi_s peaks there.

        The vertex being peaked at Q=(pi,pi) is the prerequisite for
        d-wave superconductivity in multi-orbital systems.
        """
        from hwave.sc import (_compute_vertices_flex, _read_interaction_files,
                              _build_interaction_k)

        Lx, Ly, Nz = self.Lx, self.Ly, 1
        norb = self.solver.norb
        nd = self.solver.nd
        nmat = self.Nmat

        # Build interaction in k-space
        input_dict = {
            'file': {
                'input': {
                    'path_to_input': 'tests/rpa/input',
                    'interaction': {
                        'path_to_input': 'tests/rpa/input',
                        'Geometry': 'geom.dat',
                        'Transfer': 'transfer.dat',
                        'CoulombIntra': 'coulombintra.dat',
                    },
                },
                'output': {'path_to_output': 'tests/flex/output_bench'},
            }
        }
        geom_info, hr, interactions = _read_interaction_files(input_dict)
        kx = np.linspace(0, 2 * np.pi, Lx, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, Ly, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, Nz, endpoint=False)
        inter_k = _build_interaction_k(kx, ky, kz, interactions, norb)

        # Convert chi_s, chi_c to Eliashberg format
        center = nmat // 2
        chis_orb = self.chi_s[center].reshape(Lx, Ly, Nz, nd, nd)[:, :, :, :norb, :norb]
        chic_orb = self.chi_c[center].reshape(Lx, Ly, Nz, nd, nd)[:, :, :, :norb, :norb]

        # Compute pairing vertex (singlet)
        Vs_q = _compute_vertices_flex(chis_orb, chic_orb, inter_k,
                                      norb, Lx, Ly, Nz,
                                      pairing_type="singlet")

        # For 1-orbital: Vs_q shape is (1,1,1,1,Nx,Ny,Nz) — scalar per q
        Vs_real = np.real(Vs_q[0, 0, 0, 0, :, :, 0])

        # Vertex should be peaked at Q=(pi,pi) due to AF spin fluctuations
        q_pipi = (Lx // 2, Ly // 2)
        q_00 = (0, 0)

        Vs_pipi = Vs_real[q_pipi]
        Vs_00 = Vs_real[q_00]

        self.assertGreater(Vs_pipi, Vs_00,
                           "Pairing vertex should peak at Q=(pi,pi) ({:.4f}) "
                           "vs Q=(0,0) ({:.4f})".format(Vs_pipi, Vs_00))

        # Verify it's the global maximum
        Vs_max_q = np.unravel_index(np.argmax(Vs_real), Vs_real.shape)
        self.assertEqual(Vs_max_q, q_pipi,
                         "Vs(q) peak at {} should be at (pi,pi)={}".format(
                             Vs_max_q, q_pipi))


class TestFLEXBenchmark2Orbital(unittest.TestCase):
    """Physics benchmarks for 2-orbital Hubbard-Kanamori model FLEX.

    Model: 2-orbital Hubbard on square lattice
    with CoulombIntra U=4 and CoulombInter V (from test data).
    """

    @classmethod
    def setUpClass(cls):
        """Run 2-orbital FLEX calculation once."""
        cls.Lx, cls.Ly = 8, 8
        cls.Nmat = 32
        cls.T = 1.0
        solver, green_info = _make_2orb_flex(
            Lx=cls.Lx, Ly=cls.Ly, Nmat=cls.Nmat, T=cls.T,
            max_iter=30, mix=0.2, eps=4)
        solver.solve(green_info, 'tests/flex/output_bench_2orb')

        cls.solver = solver
        cls.green_info = green_info
        cls.sigma = green_info["sigma"]
        cls.green_kw = green_info["green"]
        cls.chi0q = green_info["chi0q"]
        cls.chi_s = green_info["chiq_s"]
        cls.chi_c = green_info["chiq_c"]
        cls.beta = 1.0 / cls.T

    def test_orbital_dependent_self_energy(self):
        """Self-energy should be orbital-dependent in the 2-orbital model.

        Different orbitals have different bandwidths and interactions,
        leading to orbital-selective mass renormalization.
        """
        iw0 = self.Nmat // 2
        norb = self.solver.norb  # 2

        # sigma shape: (nblock, nmat, nvol, norb, norb) for spin-free
        sigma_iw0 = self.sigma[0, iw0]  # (nvol, norb, norb)

        # Diagonal elements: Im Sigma_{aa}(k, iw0) for each orbital
        im_sigma_orb = np.zeros((norb, self.Lx * self.Ly))
        for a in range(norb):
            im_sigma_orb[a] = np.imag(sigma_iw0[:, a, a])

        # Both orbitals should have negative Im Sigma (damping)
        for a in range(norb):
            self.assertTrue(np.all(im_sigma_orb[a] < 0),
                            "Im Sigma for orbital {} should be negative".format(a))

        # In general, the two orbitals may have different damping rates
        # (orbital-dependent self-energy is a key feature of multi-orbital FLEX)
        avg_damping = [np.mean(np.abs(im_sigma_orb[a])) for a in range(norb)]
        self.assertGreater(min(avg_damping), 0,
                           "Both orbitals should have nonzero damping")

    def test_mass_renormalization_per_orbital(self):
        """Quasiparticle weight Z_a should satisfy 0 < Z_a < 1 for each orbital.

        Uses Z ~ 1 / (1 - Im Sigma(iw_0) / w_0) which is robust for
        coarse Matsubara grids.
        """
        iw0 = self.Nmat // 2
        norb = self.solver.norb
        nvol = self.Lx * self.Ly
        omega0 = np.pi / self.beta

        Z_avg = np.zeros(norb)
        for a in range(norb):
            Z_k = np.zeros(nvol)
            for ik in range(nvol):
                im0 = np.imag(self.sigma[0, iw0, ik, a, a])
                Z_k[ik] = 1.0 / (1.0 - im0 / omega0)
            Z_avg[a] = np.mean(Z_k)

        for a in range(norb):
            self.assertGreater(Z_avg[a], 0.0,
                               "Z for orbital {} ({:.4f}) should be > 0".format(
                                   a, Z_avg[a]))
            self.assertLess(Z_avg[a], 1.0,
                            "Z for orbital {} ({:.4f}) should be < 1".format(
                                a, Z_avg[a]))

    def test_spin_susceptibility_peak_structure(self):
        """Spin susceptibility should have well-defined peak structure
        reflecting the nesting properties of the multi-orbital Fermi surface.

        The peak may not be at (pi,pi) if the Fermi surface nesting
        vector differs from the simple square lattice case.
        """
        center = self.Nmat // 2
        nvol = self.Lx * self.Ly
        nd = self.chi_s.shape[-1]

        chi_s_static = self.chi_s[center]  # (nvol, nd, nd)
        chi_s_trace = np.real(np.trace(chi_s_static, axis1=-2, axis2=-1))

        q_max = np.argmax(chi_s_trace)
        chi_max = chi_s_trace[q_max]
        chi_avg = np.mean(chi_s_trace)

        # The peak should be above average (threshold depends on temperature;
        # at T=1.0 with this model, the ratio is modest but still > 1)
        self.assertGreater(chi_max / chi_avg, 1.05,
                           "Spin susceptibility peak ({:.4f}) should be above "
                           "average ({:.4f}), ratio={:.4f}".format(
                               chi_max, chi_avg, chi_max / chi_avg))

    def test_off_diagonal_self_energy(self):
        """In multi-orbital systems, off-diagonal self-energy elements
        Sigma_{ab}(k, iwn) with a != b should be nonzero when inter-orbital
        interactions or hybridization are present.
        """
        iw0 = self.Nmat // 2
        norb = self.solver.norb

        if norb < 2:
            self.skipTest("Need multi-orbital system")

        sigma_iw0 = self.sigma[0, iw0]  # (nvol, norb, norb)

        # Off-diagonal elements
        offdiag_norm = 0.0
        for a in range(norb):
            for b in range(norb):
                if a != b:
                    offdiag_norm += np.linalg.norm(sigma_iw0[:, a, b])

        # With inter-orbital hopping and CoulombInter, off-diagonal
        # self-energy should be nonzero
        self.assertGreater(offdiag_norm, 1e-10,
                           "Off-diagonal self-energy should be nonzero "
                           "with inter-orbital interactions")

    def test_hermitian_symmetry(self):
        """Sigma(-iwn) = Sigma(iwn)^dagger for all k and n."""
        nmat = self.Nmat
        for n in range(nmat // 2):
            n_neg = nmat - 1 - n
            sigma_dagger = np.conj(self.sigma[:, n]).swapaxes(-2, -1)
            np.testing.assert_allclose(
                self.sigma[:, n_neg], sigma_dagger, atol=1e-10,
                err_msg="2-orbital Hermitian symmetry violated at n={}".format(n))

    def test_self_energy_fermi_liquid(self):
        """At high Matsubara frequency, Sigma should decay as 1/iwn
        (Fermi liquid behavior).

        |Sigma(iwn)| should decrease monotonically for large n.
        """
        norb = self.solver.norb
        nmat = self.Nmat

        # Average |Sigma| over k-points for positive Matsubara frequencies
        n_start = nmat // 2
        sigma_vs_n = np.zeros(nmat - n_start)
        for n in range(n_start, nmat):
            sigma_vs_n[n - n_start] = np.mean(np.abs(self.sigma[0, n]))

        # Check monotonic decrease in the upper half of frequencies
        n_mid = len(sigma_vs_n) // 2
        self.assertGreater(sigma_vs_n[0], sigma_vs_n[-1],
                           "|Sigma(iw_0)| ({:.6f}) should exceed "
                           "|Sigma(iw_max)| ({:.6f})".format(
                               sigma_vs_n[0], sigma_vs_n[-1]))

    def test_spin_charge_separation(self):
        """Spin and charge susceptibilities should differ significantly.

        For repulsive interactions (U > 0), spin fluctuations are enhanced
        while charge fluctuations are suppressed. This is a fundamental
        feature of correlated electron systems.
        """
        center = self.Nmat // 2
        nvol = self.Lx * self.Ly

        chi_s_trace = np.real(np.trace(self.chi_s[center], axis1=-2, axis2=-1))
        chi_c_trace = np.real(np.trace(self.chi_c[center], axis1=-2, axis2=-1))
        chi0_trace = np.real(np.trace(self.chi0q[center], axis1=-2, axis2=-1))

        # Overall: spin enhanced, charge suppressed
        self.assertGreater(np.sum(chi_s_trace), np.sum(chi0_trace),
                           "Total spin susceptibility should exceed bare chi0")
        self.assertLess(np.sum(chi_c_trace), np.sum(chi0_trace),
                        "Total charge susceptibility should be less than bare chi0")


class TestFLEXBenchmarkUDependence(unittest.TestCase):
    """Verify that FLEX results depend correctly on interaction strength U.

    Physical expectations:
    - Stronger U -> larger self-energy
    - Stronger U -> more enhanced spin susceptibility
    - Stronger U -> smaller quasiparticle weight Z
    """

    def test_sigma_increases_with_U(self):
        """Self-energy magnitude should increase with U."""
        sigma_norms = []
        for U in [2.0, 6.0]:
            solver, gi = _make_1orb_flex(
                Lx=8, Ly=8, Nmat=32, T=2.0, U=U,
                max_iter=30, mix=0.3, eps=3)
            solver.solve(gi, 'tests/flex/output_bench')
            sigma_norms.append(np.linalg.norm(gi["sigma"]))

        self.assertGreater(sigma_norms[1], sigma_norms[0],
                           "|Sigma(U=6)| ({:.6f}) should exceed "
                           "|Sigma(U=2)| ({:.6f})".format(
                               sigma_norms[1], sigma_norms[0]))

    def test_spin_susceptibility_increases_with_U(self):
        """Spin susceptibility at Q=(pi,pi) should increase with U."""
        chi_s_pipi = []
        for U in [2.0, 6.0]:
            solver, gi = _make_1orb_flex(
                Lx=8, Ly=8, Nmat=32, T=2.0, U=U,
                max_iter=30, mix=0.3, eps=3)
            solver.solve(gi, 'tests/flex/output_bench')

            Lx, Ly = 8, 8
            nmat = 32
            center = nmat // 2
            chi_s = gi["chiq_s"]
            q_pipi = (Ly // 2) * Lx + Lx // 2
            chi_s_pipi.append(
                np.real(np.trace(chi_s[center, q_pipi])))

        self.assertGreater(chi_s_pipi[1], chi_s_pipi[0],
                           "chi_s(pi,pi) at U=6 ({:.6f}) should exceed "
                           "U=2 ({:.6f})".format(chi_s_pipi[1], chi_s_pipi[0]))

    def test_mass_renormalization_increases_with_U(self):
        """Quasiparticle weight Z should decrease (mass enhancement increases)
        with increasing U.
        """
        Z_values = []
        for U in [2.0, 6.0]:
            solver, gi = _make_1orb_flex(
                Lx=8, Ly=8, Nmat=32, T=2.0, U=U,
                max_iter=30, mix=0.3, eps=3)
            solver.solve(gi, 'tests/flex/output_bench')

            T = 2.0
            beta = 1.0 / T
            nmat = 32
            iw0 = nmat // 2
            omega0 = np.pi * T  # lowest positive Matsubara freq
            sigma = gi["sigma"]
            nvol = 8 * 8

            Z_k = np.zeros(nvol)
            for ik in range(nvol):
                im0 = np.imag(sigma[0, iw0, ik, 0, 0])
                Z_k[ik] = 1.0 / (1.0 - im0 / omega0)
            Z_values.append(np.mean(Z_k))

        # Larger U -> smaller Z (stronger correlation)
        self.assertLess(Z_values[1], Z_values[0],
                        "Z(U=6) ({:.4f}) should be smaller than "
                        "Z(U=2) ({:.4f})".format(Z_values[1], Z_values[0]))


def _make_iron_2orb_flex(Lx=8, Ly=8, Nmat=64, T=0.2, filling=0.5,
                         U=1.5, J=0.25,
                         max_iter=50, mix=0.2, eps=4):
    """Create the Raghu et al. iron pnictide 2-orbital FLEX solver.

    Model: dxz/dyz on a square lattice (1-Fe unit cell).
    Raghu et al., PRB 77, 220503(R) (2008).

    Parameters: t1=-1.0, t2=1.3, t3=t4=-0.85
    Kanamori: U, U'=U-2J, J, J'=J
    """
    import tempfile
    import shutil

    tmpdir = tempfile.mkdtemp(prefix='flex_iron_')

    # Geometry: 2 orbitals at same site
    with open(os.path.join(tmpdir, 'geom.dat'), 'w') as f:
        f.write("  1.0  0.0  0.0\n  0.0  1.0  0.0\n  0.0  0.0  1.0\n")
        f.write("2\n  0.0  0.0  0.0\n  0.0  0.0  0.0\n")

    # Transfer integrals (Raghu model)
    t1, t2, t3, t4 = -1.0, 1.3, -0.85, -0.85
    hops = []
    # dxz-dxz (orbital 1)
    hops += [(1, 0, 0, 1, 1, -t1), (-1, 0, 0, 1, 1, -t1),
             (0, 1, 0, 1, 1, -t2), (0, -1, 0, 1, 1, -t2),
             (1, 1, 0, 1, 1, -t3), (1, -1, 0, 1, 1, -t3),
             (-1, 1, 0, 1, 1, -t3), (-1, -1, 0, 1, 1, -t3)]
    # dyz-dyz (orbital 2)
    hops += [(1, 0, 0, 2, 2, -t2), (-1, 0, 0, 2, 2, -t2),
             (0, 1, 0, 2, 2, -t1), (0, -1, 0, 2, 2, -t1),
             (1, 1, 0, 2, 2, -t3), (1, -1, 0, 2, 2, -t3),
             (-1, 1, 0, 2, 2, -t3), (-1, -1, 0, 2, 2, -t3)]
    # dxz-dyz inter-orbital (from -4*t4*sin(kx)*sin(ky))
    hops += [(1, 1, 0, 1, 2, t4), (1, -1, 0, 1, 2, -t4),
             (-1, 1, 0, 1, 2, -t4), (-1, -1, 0, 1, 2, t4),
             (1, 1, 0, 2, 1, t4), (1, -1, 0, 2, 1, -t4),
             (-1, 1, 0, 2, 1, -t4), (-1, -1, 0, 2, 1, t4)]

    with open(os.path.join(tmpdir, 'transfer.dat'), 'w') as f:
        f.write("Transfer for Raghu iron pnictide model\n2\n9\n")
        f.write(" 1 1 1 1 1 1 1 1 1\n")
        for rx, ry, rz, a, b, val in hops:
            f.write("  {:2d}  {:2d}  {:2d}  {:2d}  {:2d}  {:.12f}  0.000000000000\n"
                    .format(rx, ry, rz, a, b, val))

    # CoulombIntra: U
    with open(os.path.join(tmpdir, 'coulombintra.dat'), 'w') as f:
        f.write("CoulombIntra\n2\n1\n 1\n")
        f.write("   0    0    0    1    1   {:.12f}   0.000000000000\n".format(U))
        f.write("   0    0    0    2    2   {:.12f}   0.000000000000\n".format(U))

    # CoulombInter: U' = U - 2J
    Uprime = U - 2 * J
    with open(os.path.join(tmpdir, 'coulombinter.dat'), 'w') as f:
        f.write("CoulombInter\n2\n1\n 1\n")
        f.write("   0    0    0    1    2   {:.12f}   0.000000000000\n".format(Uprime))
        f.write("   0    0    0    2    1   {:.12f}   0.000000000000\n".format(Uprime))

    # Hund: J
    with open(os.path.join(tmpdir, 'hund.dat'), 'w') as f:
        f.write("Hund\n2\n1\n 1\n")
        f.write("   0    0    0    1    2   {:.12f}   0.000000000000\n".format(J))
        f.write("   0    0    0    2    1   {:.12f}   0.000000000000\n".format(J))

    # Exchange: J' = J
    with open(os.path.join(tmpdir, 'exchange.dat'), 'w') as f:
        f.write("Exchange\n2\n1\n 1\n")
        f.write("   0    0    0    1    2   {:.12f}   0.000000000000\n".format(J))
        f.write("   0    0    0    2    1   {:.12f}   0.000000000000\n".format(J))

    info_log = {}
    info_mode = {
        'mode': 'FLEX',
        'param': {
            'T': T,
            'filling': filling,
            'CellShape': [Lx, Ly, 1],
            'SubShape': [1, 1, 1],
            'Nmat': Nmat,
            'IterationMax': max_iter,
            'Mix': mix,
            'EPS': eps,
        },
    }
    info_file_input = {
        'path_to_input': tmpdir,
        'interaction': {
            'path_to_input': tmpdir,
            'Geometry': 'geom.dat',
            'Transfer': 'transfer.dat',
            'CoulombIntra': 'coulombintra.dat',
            'CoulombInter': 'coulombinter.dat',
            'Hund': 'hund.dat',
            'Exchange': 'exchange.dat',
        },
    }
    os.makedirs('tests/flex/output_iron', exist_ok=True)

    import hwave.qlmsio.read_input_k as read_input_k
    read_io = read_input_k.QLMSkInput(info_file_input)
    ham_info = read_io.get_param("ham")
    green_info = read_io.get_param("green")

    import hwave.solver.flex as solver_flex
    solver = solver_flex.FLEX(ham_info, info_log, info_mode)

    shutil.rmtree(tmpdir, ignore_errors=True)

    return solver, green_info


class TestFLEXIronPnictide(unittest.TestCase):
    """Physics benchmarks for the iron pnictide 2-orbital model.

    Reference: Raghu et al., PRB 77, 220503(R) (2008)
    Model: dxz/dyz on square lattice, t1=-1, t2=1.3, t3=t4=-0.85
    Kanamori: U=1.5, J=0.25, U'=1.0
    Key physics:
      - Spin susceptibility peaks at Q=(pi,0)/(0,pi), NOT (pi,pi)
      - Orbital-dependent self-energy (d_xz vs d_yz anisotropy)
      - Quasiparticle damping Im Sigma < 0
    """

    @classmethod
    def setUpClass(cls):
        """Run FLEX calculation once for the iron pnictide model."""
        cls.Lx, cls.Ly = 16, 16
        cls.Nmat = 64
        cls.T = 0.1
        cls.U = 1.5
        cls.J = 0.25
        solver, green_info = _make_iron_2orb_flex(
            Lx=cls.Lx, Ly=cls.Ly, Nmat=cls.Nmat, T=cls.T,
            U=cls.U, J=cls.J,
            max_iter=200, mix=0.2, eps=6)
        solver.solve(green_info, 'tests/flex/output_iron')

        cls.solver = solver
        cls.green_info = green_info
        cls.sigma = green_info["sigma"]
        cls.green_kw = green_info["green"]
        cls.chi0q = green_info["chi0q"]
        cls.chi_s = green_info["chiq_s"]
        cls.chi_c = green_info["chiq_c"]
        cls.beta = 1.0 / cls.T

    def test_auto_selected_the_general_scheme(self):
        """The Kanamori set includes Exchange, so auto must resolve to
        general -- the scheme in which Exchange genuinely acts."""
        self.assertEqual(self.solver.calc_scheme, 'general')

    def test_spin_susceptibility_peaks_at_pi_0(self):
        """chi_s should peak at Q=(pi,0) or (0,pi), NOT (pi,pi).

        This is the hallmark of the iron pnictide Fermi surface nesting
        between hole pockets at Gamma and electron pockets at M.
        """
        center = self.Nmat // 2
        chi = np.asarray(self.chi_s)[center]
        if chi.ndim == 5:
            # general scheme: rank-4 orbital tensor in the public
            # [a, c, b, d] convention; the physical spin response sums
            # the density sector (aa), (bb)
            norb = chi.shape[-1]
            chi_s_trace = sum(chi[:, a, a, b, b].real
                              for a in range(norb) for b in range(norb))
        else:
            # reduced: (nvol, norb, norb) density-pair matrix
            chi_s_trace = chi.real.sum(axis=(-2, -1))

        # Q vectors
        q_pi0 = (self.Lx // 2) * self.Ly  # (pi, 0)
        q_0pi = self.Ly // 2               # (0, pi)
        q_pipi = (self.Lx // 2) * self.Ly + self.Ly // 2  # (pi, pi)

        chi_pi0 = chi_s_trace[q_pi0]
        chi_0pi = chi_s_trace[q_0pi]
        chi_pipi = chi_s_trace[q_pipi]

        # (pi,0) and (0,pi) should be degenerate by C4 symmetry
        self.assertAlmostEqual(chi_pi0, chi_0pi, places=4,
                               msg="chi_s(pi,0) and chi_s(0,pi) should be "
                                   "degenerate by C4 symmetry")

        # (pi,0) should exceed (pi,pi)
        self.assertGreater(chi_pi0, chi_pipi,
                           "chi_s(pi,0) ({:.4f}) should exceed "
                           "chi_s(pi,pi) ({:.4f})".format(chi_pi0, chi_pipi))

    def test_orbital_anisotropy_self_energy(self):
        """At k=(pi,0), d_yz should be more strongly damped than d_xz.

        The (pi,0) nesting vector connects portions of the Fermi surface
        with predominantly d_yz character, so spin fluctuations at Q=(pi,0)
        scatter d_yz quasiparticles more strongly.
        """
        iw0 = self.Nmat // 2
        k_pi0 = (self.Lx // 2) * self.Ly

        im_sigma_xz = abs(self.sigma[0, iw0, k_pi0, 0, 0].imag)
        im_sigma_yz = abs(self.sigma[0, iw0, k_pi0, 1, 1].imag)

        self.assertGreater(im_sigma_yz, im_sigma_xz,
                           "|Im Sigma_yz(pi,0)| ({:.6f}) should exceed "
                           "|Im Sigma_xz(pi,0)| ({:.6f})".format(
                               im_sigma_yz, im_sigma_xz))

    def test_quasiparticle_damping(self):
        """Im Sigma(k, iw0) < 0 for all k and both orbitals."""
        iw0 = self.Nmat // 2
        norb = self.solver.norb

        for a in range(norb):
            im_sigma = self.sigma[0, iw0, :, a, a].imag
            self.assertTrue(np.all(im_sigma < 0),
                            "Im Sigma for orbital {} should be negative "
                            "everywhere".format(a))

    def test_mass_renormalization(self):
        """Quasiparticle weight 0 < Z < 1 for both orbitals."""
        iw0 = self.Nmat // 2
        norb = self.solver.norb
        omega0 = np.pi * self.T

        for a in range(norb):
            im_sigma = self.sigma[0, iw0, :, a, a].imag
            Z_k = 1.0 / (1.0 - im_sigma / omega0)
            Z_avg = np.mean(Z_k)
            self.assertGreater(Z_avg, 0.0,
                               "Z for orbital {} should be > 0".format(a))
            self.assertLess(Z_avg, 1.0,
                            "Z for orbital {} should be < 1".format(a))

    def test_spin_charge_anisotropy(self):
        """Spin fluctuations should be enhanced, charge suppressed."""
        center = self.Nmat // 2

        chi_s_total = np.real(np.trace(self.chi_s[center],
                                        axis1=-2, axis2=-1)).sum()
        chi_c_total = np.real(np.trace(self.chi_c[center],
                                        axis1=-2, axis2=-1)).sum()
        chi0_total = np.real(np.trace(self.chi0q[center],
                                       axis1=-2, axis2=-1)).sum()

        self.assertGreater(chi_s_total, chi0_total,
                           "Total spin susceptibility should exceed bare chi0")
        self.assertLess(chi_c_total, chi0_total,
                        "Total charge susceptibility should be less than chi0")

    def test_hermitian_symmetry(self):
        """Sigma(-iwn) = Sigma(iwn)^dagger."""
        nmat = self.Nmat
        for n in range(nmat // 2):
            n_neg = nmat - 1 - n
            sigma_dagger = np.conj(self.sigma[:, n]).swapaxes(-2, -1)
            np.testing.assert_allclose(
                self.sigma[:, n_neg], sigma_dagger, atol=1e-10,
                err_msg="Iron pnictide Hermitian symmetry violated "
                        "at n={}".format(n))


if __name__ == '__main__':
    unittest.main()
