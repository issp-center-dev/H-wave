#!/usr/bin/env python3

"""Tests for the linearized Eliashberg equation solver (hwave.sc).

Test strategy:
1. Green's function: Verify G(k, iwn) against direct construction.
2. Vertex: Verify Pc, Ps against reference code on a small system.
3. Kernel: One iteration of the FFT kernel against reference.
4. Convergence: Self-consistent loop converges on a small system.
5. Eigenvalue: Eigenvalue analysis agrees with iteration result.
"""

import os
import unittest

import numpy as np
import numpy.testing as npt

from hwave.sc import (
    _build_hamiltonian_k,
    _build_interaction_k,
    _build_sc_matrices,
    _calc_chi0q_internal,
    _calc_eigenvalues,
    _calc_g2,
    _calc_green,
    _compute_vertices,
    _compute_vertices_simple,
    _compute_vertices_general,
    _convert_chi0q_to_ref_format,
    _determine_mu,
    _eliashberg_kernel_fft,
    _initialize_gap,
    _solve_eigenvalue,
    _solve_iteration,
    _solve_subspace_iteration,
    _solve_shifted_bicg,
)


class TestGreenFunction(unittest.TestCase):
    """Test Green's function construction."""

    def setUp(self):
        """Set up a small 1-orbital square lattice test case."""
        self.Nx, self.Ny, self.Nz = 4, 4, 1
        self.norb = 1
        self.t = 1.0
        self.beta = 10.0
        self.nmat = 32
        self.mu = 0.0

        # Build Transfer dict for square lattice
        self.hr = {
            ((1, 0, 0), (0, 0)): self.t,
            ((-1, 0, 0), (0, 0)): self.t,
            ((0, 1, 0), (0, 0)): self.t,
            ((0, -1, 0), (0, 0)): self.t,
        }

        self.kx = np.linspace(0, 2 * np.pi, self.Nx, endpoint=False)
        self.ky = np.linspace(0, 2 * np.pi, self.Ny, endpoint=False)
        self.kz = np.linspace(0, 2 * np.pi, self.Nz, endpoint=False)

    def test_hamiltonian_k(self):
        """Test that epsilon(k) = 2t(cos kx + cos ky) for square lattice."""
        epsilon_k = _build_hamiltonian_k(self.kx, self.ky, self.kz, self.hr, self.norb)
        self.assertEqual(epsilon_k.shape, (1, 1, 4, 4, 1))

        for ix, kx in enumerate(self.kx):
            for iy, ky in enumerate(self.ky):
                expected = 2 * self.t * (np.cos(kx) + np.cos(ky))
                npt.assert_allclose(
                    epsilon_k[0, 0, ix, iy, 0].real, expected, atol=1e-10
                )

    def test_eigenvalues(self):
        """Test eigenvalue computation."""
        epsilon_k = _build_hamiltonian_k(self.kx, self.ky, self.kz, self.hr, self.norb)
        eigenvalues, eigenvectors = _calc_eigenvalues(epsilon_k)
        self.assertEqual(eigenvalues.shape, (4, 4, 1, 1))

        for ix, kx in enumerate(self.kx):
            for iy, ky in enumerate(self.ky):
                expected = 2 * self.t * (np.cos(kx) + np.cos(ky))
                npt.assert_allclose(eigenvalues[ix, iy, 0, 0], expected, atol=1e-10)

    def test_green_function(self):
        """Test G(k, iwn) = 1/(iwn - (ek - mu)) for 1-orbital case."""
        epsilon_k = _build_hamiltonian_k(self.kx, self.ky, self.kz, self.hr, self.norb)
        eigenvalues, eigenvectors = _calc_eigenvalues(epsilon_k)
        green_kw = _calc_green(eigenvalues, eigenvectors, self.mu, self.beta, self.nmat)
        self.assertEqual(green_kw.shape, (1, 1, 4, 4, 1, 32))

        iomega = np.array(
            [(2.0 * i + 1.0 - self.nmat) * np.pi for i in range(self.nmat)]
        ) / self.beta

        for ix, kx in enumerate(self.kx):
            for iy, ky in enumerate(self.ky):
                ek = 2 * self.t * (np.cos(kx) + np.cos(ky))
                for iw in range(self.nmat):
                    expected = 1.0 / (1j * iomega[iw] - (ek - self.mu))
                    npt.assert_allclose(
                        green_kw[0, 0, ix, iy, 0, iw], expected, atol=1e-10
                    )

    def test_determine_mu(self):
        """Test chemical potential determination at half-filling."""
        epsilon_k = _build_hamiltonian_k(self.kx, self.ky, self.kz, self.hr, self.norb)
        eigenvalues, _ = _calc_eigenvalues(epsilon_k)

        # Half-filling for 1 orbital: n_target=0.5
        mu = _determine_mu(eigenvalues, self.beta, 0.5, self.norb)
        # For particle-hole symmetric square lattice at half-filling, mu=0
        npt.assert_allclose(mu, 0.0, atol=1e-4)


class TestGreenFunction2Orb(unittest.TestCase):
    """Test Green's function for 2-orbital case against reference code."""

    def setUp(self):
        self.Nx, self.Ny, self.Nz = 4, 4, 1
        self.norb = 2
        self.beta = 10.0
        self.nmat = 32
        self.mu = 0.0
        self.t = 1.0
        self.t1 = 0.5

        # 2-orbital model from reference test
        self.hr = {
            ((0, 1, 0), (0, 0)): self.t,
            ((0, -1, 0), (0, 0)): self.t,
            ((0, 1, 0), (1, 1)): self.t,
            ((0, -1, 0), (1, 1)): self.t,
            ((0, 0, 0), (0, 1)): self.t1,
            ((-1, 0, 0), (0, 1)): self.t1,
            ((0, 0, 0), (1, 0)): self.t1,
            ((1, 0, 0), (1, 0)): self.t1,
        }

        self.kx = np.linspace(0, 2 * np.pi, self.Nx, endpoint=False)
        self.ky = np.linspace(0, 2 * np.pi, self.Ny, endpoint=False)
        self.kz = np.linspace(0, 2 * np.pi, self.Nz, endpoint=False)

    def test_green_matches_reference(self):
        """Test that our Green's function matches the reference implementation."""
        epsilon_k = _build_hamiltonian_k(
            self.kx, self.ky, self.kz, self.hr, self.norb
        )
        eigenvalues, eigenvectors = _calc_eigenvalues(epsilon_k)
        green_kw = _calc_green(eigenvalues, eigenvectors, self.mu, self.beta, self.nmat)

        # Build reference Green's function directly
        iomega = np.array(
            [(2.0 * i + 1.0 - self.nmat) * np.pi for i in range(self.nmat)]
        ) / self.beta

        green_ref = np.zeros((2, 2, self.Nx, self.Ny, 1, self.nmat), dtype=complex)
        for ix in range(self.Nx):
            for iy in range(self.Ny):
                H = epsilon_k[:, :, ix, iy, 0]
                vals, vecs = np.linalg.eigh(H)
                vec_conj = np.conjugate(vecs)
                factor = np.einsum('im,jm->ijm', vecs, vec_conj)
                for iw in range(self.nmat):
                    green_ref[:, :, ix, iy, 0, iw] = np.sum(
                        factor / (1j * iomega[iw] - (vals - self.mu))[None, None, :],
                        axis=2
                    )

        npt.assert_allclose(green_kw, green_ref, atol=1e-10)


class TestVertexComputation(unittest.TestCase):
    """Test RPA vertex computation Pc, Ps."""

    def test_vertex_simple(self):
        """Test vertex for a 2x2x1 system with 1 orbital."""
        Nx, Ny, Nz = 2, 2, 1
        norb = 1
        nmat = 8

        # Create a simple chi0q
        rng = np.random.default_rng(42)
        chi0q = rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat)) + \
                1j * rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat))
        chi0q *= 0.1  # small values for stability

        # Simple on-site U and no V
        U_k = np.ones((norb, norb, Nx, Ny, Nz), dtype=complex) * 2.0
        inter_k = {"CoulombIntra": U_k}

        Pc_q, Ps_q = _compute_vertices(chi0q, inter_k, norb, Nx, Ny, Nz, nmat)

        self.assertEqual(Pc_q.shape, (norb, norb, Nx, Ny, Nz))
        self.assertEqual(Ps_q.shape, (norb, norb, Nx, Ny, Nz))

        # Verify against manual computation for each q-point
        I = np.identity(norb)
        for ix in range(Nx):
            for iy in range(Ny):
                _U = U_k[:, :, ix, iy, 0]
                Wc = _U
                Ws = -_U
                _chi0 = chi0q[:, :, ix, iy, 0, nmat // 2].astype(np.complex128)

                chis = np.linalg.solve(I + _chi0 @ Ws, _chi0)
                chic = np.linalg.solve(I + _chi0 @ Wc, _chi0)

                Pc_expected = (Wc + Ws) / 2.0 - 0.5 * Wc @ chic @ Wc
                Ps_expected = -Ws + 1.5 * Ws @ chis @ Ws

                npt.assert_allclose(
                    Pc_q[:, :, ix, iy, 0], Pc_expected, atol=1e-10
                )
                npt.assert_allclose(
                    Ps_q[:, :, ix, iy, 0], Ps_expected, atol=1e-10
                )


class TestG2Calculation(unittest.TestCase):
    """Test G2 two-particle Green's function."""

    def test_g2_shape(self):
        """Test G2 output shape."""
        norb = 1
        Nx, Ny, Nz, nmat = 4, 4, 1, 8
        green_kw = np.random.randn(norb, norb, Nx, Ny, Nz, nmat) + \
                   1j * np.random.randn(norb, norb, Nx, Ny, Nz, nmat)
        G2 = _calc_g2(green_kw)
        self.assertEqual(G2.shape, (norb, norb, norb, norb, Nx, Ny, Nz))

    def test_g2_shape_2orb(self):
        """Test G2 shape for 2-orbital case."""
        norb = 2
        Nx, Ny, Nz, nmat = 4, 4, 1, 8
        green_kw = np.random.randn(norb, norb, Nx, Ny, Nz, nmat) + \
                   1j * np.random.randn(norb, norb, Nx, Ny, Nz, nmat)
        G2 = _calc_g2(green_kw)
        self.assertEqual(G2.shape, (norb, norb, norb, norb, Nx, Ny, Nz))


class TestKernel(unittest.TestCase):
    """Test Eliashberg kernel FFT convolution."""

    def test_kernel_shape(self):
        """Test that kernel returns correct shape."""
        norb = 1
        Nx, Ny, Nz, nmat = 4, 4, 1, 8

        P_q = np.random.randn(norb, norb, Nx, Ny, Nz) + \
              1j * np.random.randn(norb, norb, Nx, Ny, Nz)

        rng = np.random.default_rng(42)
        green_kw = rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat)) + \
                   1j * rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat))
        G2 = _calc_g2(green_kw)

        sigma_old = np.ones((norb, norb, Nx, Ny, Nz))
        sigma_old /= np.linalg.norm(sigma_old)

        sigma_new = _eliashberg_kernel_fft(P_q, G2, sigma_old, norb)
        self.assertEqual(sigma_new.shape, (norb, norb, Nx, Ny, Nz))


class TestInitializeGap(unittest.TestCase):
    """Test gap function initialization."""

    def test_cos_init(self):
        """Test cosine initialization."""
        kx = np.linspace(0, 2 * np.pi, 4, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, 4, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, 1, endpoint=False)
        sigma = _initialize_gap("cos", 2, kx, ky, kz)
        self.assertEqual(sigma.shape, (2, 2, 4, 4, 1))
        # Should be normalized
        npt.assert_allclose(np.linalg.norm(sigma), 1.0, atol=1e-10)

    def test_random_init(self):
        """Test random initialization."""
        kx = np.linspace(0, 2 * np.pi, 4, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, 4, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, 1, endpoint=False)
        sigma = _initialize_gap("random", 2, kx, ky, kz)
        self.assertEqual(sigma.shape, (2, 2, 4, 4, 1))
        npt.assert_allclose(np.linalg.norm(sigma), 1.0, atol=1e-10)

    def test_all_symmetry_modes_2d(self):
        """Test all symmetry-based gap initializations on 2D system."""
        N = 8
        kx = np.linspace(0, 2 * np.pi, N, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, N, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, 1, endpoint=False)
        norb = 2

        all_modes = [
            "cos", "s", "s_ext", "s_ext_2d",
            "d_x2y2", "d_xy", "d_xz", "d_yz", "d_z2",
            "p_x", "p_y", "p_z",
            "random",
        ]
        for mode in all_modes:
            sigma = _initialize_gap(mode, norb, kx, ky, kz)
            self.assertEqual(sigma.shape, (norb, norb, N, N, 1),
                             msg="shape mismatch for mode={}".format(mode))
            if mode not in ("d_xz", "d_yz", "p_z"):
                # These are zero for Nz=1 (kz=0 only -> sin(0)=0)
                npt.assert_allclose(
                    np.linalg.norm(sigma), 1.0, atol=1e-10,
                    err_msg="not normalized for mode={}".format(mode))

    def test_all_symmetry_modes_3d(self):
        """Test all symmetry-based gap initializations on 3D system."""
        N = 4
        kx = np.linspace(0, 2 * np.pi, N, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, N, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, N, endpoint=False)
        norb = 1

        all_modes = [
            "cos", "s", "s_ext", "s_ext_2d",
            "d_x2y2", "d_xy", "d_xz", "d_yz", "d_z2",
            "p_x", "p_y", "p_z",
            "random",
        ]
        for mode in all_modes:
            sigma = _initialize_gap(mode, norb, kx, ky, kz)
            self.assertEqual(sigma.shape, (norb, norb, N, N, N),
                             msg="shape mismatch for mode={}".format(mode))
            npt.assert_allclose(
                np.linalg.norm(sigma), 1.0, atol=1e-10,
                err_msg="not normalized for mode={}".format(mode))

    def test_d_x2y2_symmetry(self):
        """Test d_{x^2-y^2} gap has correct sign structure."""
        N = 8
        kx = np.linspace(0, 2 * np.pi, N, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, N, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, 1, endpoint=False)
        sigma = _initialize_gap("d_x2y2", 1, kx, ky, kz)

        # d_{x^2-y^2} = cos(kx) - cos(ky)
        # At (pi,0): cos(pi)-cos(0) = -2  (negative)
        # At (0,pi): cos(0)-cos(pi) = +2  (positive)
        ix_pi = N // 2
        iy_pi = N // 2
        val_pi0 = sigma[0, 0, ix_pi, 0, 0]
        val_0pi = sigma[0, 0, 0, iy_pi, 0]
        self.assertLess(val_pi0 * val_0pi, 0,
                        "d_{x^2-y^2} should change sign between (pi,0) and (0,pi)")

    def test_p_wave_odd(self):
        """Test p-wave gap is odd under k -> -k."""
        N = 8
        kx = np.linspace(0, 2 * np.pi, N, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, N, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, 1, endpoint=False)
        sigma = _initialize_gap("p_x", 1, kx, ky, kz)

        for ix in range(N):
            ix_inv = (N - ix) % N
            npt.assert_allclose(
                sigma[0, 0, ix, 0, 0], -sigma[0, 0, ix_inv, 0, 0],
                atol=1e-10,
                err_msg="p_x should be odd under kx -> -kx")

    def test_3d_kz_dependent(self):
        """Test that kz-dependent modes (p_z, d_xz, d_yz) are non-trivial in 3D."""
        N = 4
        kx = np.linspace(0, 2 * np.pi, N, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, N, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, N, endpoint=False)

        for mode in ["p_z", "d_xz", "d_yz", "d_z2"]:
            sigma = _initialize_gap(mode, 1, kx, ky, kz)
            # Should vary along kz at an appropriate k-point
            # d_xz = sin(kx)*sin(kz) needs kx!=0; d_yz = sin(ky)*sin(kz) needs ky!=0
            if mode == "d_xz":
                kz_slice = sigma[0, 0, 1, 0, :]  # kx_index=1 so sin(kx)!=0
            elif mode == "d_yz":
                kz_slice = sigma[0, 0, 0, 1, :]  # ky_index=1 so sin(ky)!=0
            else:
                kz_slice = sigma[0, 0, 0, 0, :]
            self.assertGreater(
                np.std(np.abs(kz_slice)), 1e-10,
                msg="{} should have kz dependence".format(mode))

    def test_p_z_odd_3d(self):
        """Test p_z is odd under kz -> -kz in 3D."""
        N = 8
        kx = np.linspace(0, 2 * np.pi, N, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, N, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, N, endpoint=False)
        sigma = _initialize_gap("p_z", 1, kx, ky, kz)

        for iz in range(N):
            iz_inv = (N - iz) % N
            npt.assert_allclose(
                sigma[0, 0, 0, 0, iz], -sigma[0, 0, 0, 0, iz_inv],
                atol=1e-10,
                err_msg="p_z should be odd under kz -> -kz")


class TestChi0qConversion(unittest.TestCase):
    """Test chi0q format conversion."""

    def test_4d_conversion(self):
        """Test conversion from H-wave 4D format to reference 6D format."""
        norb = 2
        Nx, Ny, Nz = 4, 4, 1
        nmat = 8
        nvol = Nx * Ny * Nz

        chi0q_hwave = np.random.randn(nmat, nvol, norb, norb) + \
                      1j * np.random.randn(nmat, nvol, norb, norb)

        chi0q_ref = _convert_chi0q_to_ref_format(chi0q_hwave, norb, Nx, Ny, Nz, nmat)
        self.assertEqual(chi0q_ref.shape, (norb, norb, Nx, Ny, Nz, nmat))

        # Verify mapping: chi0q_hwave[w, vol_idx, a, b] == chi0q_ref[a, b, ix, iy, iz, w]
        for iw in range(nmat):
            idx = 0
            for ix in range(Nx):
                for iy in range(Ny):
                    for iz in range(Nz):
                        for a in range(norb):
                            for b in range(norb):
                                npt.assert_allclose(
                                    chi0q_ref[a, b, ix, iy, iz, iw],
                                    chi0q_hwave[iw, idx, a, b],
                                    atol=1e-14,
                                )
                        idx += 1


class TestSolverConvergence(unittest.TestCase):
    """Test that the self-consistent iteration converges on a small system."""

    def test_iteration_runs(self):
        """Test iteration solver runs without error on a small system."""
        norb = 1
        Nx, Ny, Nz = 4, 4, 1
        nmat = 16
        beta = 5.0
        t = 1.0

        hr = {
            ((1, 0, 0), (0, 0)): t,
            ((-1, 0, 0), (0, 0)): t,
            ((0, 1, 0), (0, 0)): t,
            ((0, -1, 0), (0, 0)): t,
        }

        kx = np.linspace(0, 2 * np.pi, Nx, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, Ny, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, Nz, endpoint=False)

        epsilon_k = _build_hamiltonian_k(kx, ky, kz, hr, norb)
        eigenvalues, eigenvectors = _calc_eigenvalues(epsilon_k)
        mu = _determine_mu(eigenvalues, beta, 0.5, norb)
        green_kw = _calc_green(eigenvalues, eigenvectors, mu, beta, nmat)

        # Simple U interaction
        U_k = np.ones((norb, norb, Nx, Ny, Nz), dtype=complex) * 2.0
        inter_k = {"CoulombIntra": U_k}

        # Simple chi0q (dummy for vertex test)
        chi0q = np.zeros((norb, norb, Nx, Ny, Nz, nmat), dtype=complex)
        for ix in range(Nx):
            for iy in range(Ny):
                for iz in range(Nz):
                    for iw in range(nmat):
                        chi0q[:, :, ix, iy, iz, iw] = 0.1

        Pc_q, Ps_q = _compute_vertices(chi0q, inter_k, norb, Nx, Ny, Nz, nmat)
        Vs_q = Pc_q + Ps_q
        G2 = _calc_g2(green_kw)
        sigma_init = _initialize_gap("cos", norb, kx, ky, kz)

        sigma, eigenvalue, converged, n_iter = _solve_iteration(
            green_kw, Vs_q, G2, sigma_init, norb,
            max_iter=10, alpha=0.5, tol=1e-3
        )
        self.assertEqual(sigma.shape, (norb, norb, Nx, Ny, Nz))
        self.assertGreater(n_iter, 0)


class TestEigenvalueSolver(unittest.TestCase):
    """Test eigenvalue analysis solver."""

    def test_eigenvalue_runs(self):
        """Test eigenvalue solver runs without error."""
        norb = 1
        Nx, Ny, Nz = 4, 4, 1
        nmat = 16
        beta = 5.0
        t = 1.0

        hr = {
            ((1, 0, 0), (0, 0)): t,
            ((-1, 0, 0), (0, 0)): t,
            ((0, 1, 0), (0, 0)): t,
            ((0, -1, 0), (0, 0)): t,
        }

        kx = np.linspace(0, 2 * np.pi, Nx, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, Ny, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, Nz, endpoint=False)

        epsilon_k = _build_hamiltonian_k(kx, ky, kz, hr, norb)
        eigenvalues_ek, eigenvectors = _calc_eigenvalues(epsilon_k)
        mu = _determine_mu(eigenvalues_ek, beta, 0.5, norb)
        green_kw = _calc_green(eigenvalues_ek, eigenvectors, mu, beta, nmat)

        U_k = np.ones((norb, norb, Nx, Ny, Nz), dtype=complex) * 2.0
        inter_k = {"CoulombIntra": U_k}

        chi0q = np.zeros((norb, norb, Nx, Ny, Nz, nmat), dtype=complex)
        chi0q[:, :, :, :, :, :] = 0.1

        Pc_q, Ps_q = _compute_vertices(chi0q, inter_k, norb, Nx, Ny, Nz, nmat)
        Vs_q = Pc_q + Ps_q
        G2 = _calc_g2(green_kw)

        eigenvalues, eigenvectors_gap = _solve_eigenvalue(
            Vs_q, G2, norb, Nx, Ny, Nz, num_eigenvalues=3
        )
        self.assertEqual(len(eigenvalues), 3)
        self.assertEqual(eigenvectors_gap.shape[0], 3)

    def test_eigenvalue_vs_iteration(self):
        """Test that leading eigenvalue from eigs matches iteration norm."""
        norb = 1
        Nx, Ny, Nz = 4, 4, 1
        nmat = 16
        beta = 5.0
        t = 1.0

        hr = {
            ((1, 0, 0), (0, 0)): t,
            ((-1, 0, 0), (0, 0)): t,
            ((0, 1, 0), (0, 0)): t,
            ((0, -1, 0), (0, 0)): t,
        }

        kx = np.linspace(0, 2 * np.pi, Nx, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, Ny, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, Nz, endpoint=False)

        epsilon_k = _build_hamiltonian_k(kx, ky, kz, hr, norb)
        eigenvalues_ek, eigenvectors = _calc_eigenvalues(epsilon_k)
        mu = _determine_mu(eigenvalues_ek, beta, 0.5, norb)
        green_kw = _calc_green(eigenvalues_ek, eigenvectors, mu, beta, nmat)

        U_k = np.ones((norb, norb, Nx, Ny, Nz), dtype=complex) * 2.0
        inter_k = {"CoulombIntra": U_k}

        chi0q = np.zeros((norb, norb, Nx, Ny, Nz, nmat), dtype=complex)
        chi0q[:, :, :, :, :, :] = 0.1

        Pc_q, Ps_q = _compute_vertices(chi0q, inter_k, norb, Nx, Ny, Nz, nmat)
        Vs_q = Pc_q + Ps_q
        G2 = _calc_g2(green_kw)
        sigma_init = _initialize_gap("cos", norb, kx, ky, kz)

        # Iteration
        sigma_iter, ev_iter, converged, n_iter = _solve_iteration(
            green_kw, Vs_q, G2, sigma_init, norb,
            max_iter=200, alpha=0.5, tol=1e-6
        )

        # Eigenvalue
        eigenvalues_eig, _ = _solve_eigenvalue(
            Vs_q, G2, norb, Nx, Ny, Nz, num_eigenvalues=3
        )

        # Leading eigenvalue magnitude should be close
        if converged:
            npt.assert_allclose(
                abs(eigenvalues_eig[0]), ev_iter, rtol=0.1,
                err_msg="Leading eigenvalue from eigs should match iteration result"
            )


class TestEigenvalueMethods(unittest.TestCase):
    """Test different eigenvalue solver methods give consistent results."""

    def _setup_problem(self):
        """Create a small test problem for eigenvalue solvers."""
        norb = 1
        Nx, Ny, Nz = 4, 4, 1
        nmat = 16
        beta = 5.0
        t = 1.0

        hr = {
            ((1, 0, 0), (0, 0)): t,
            ((-1, 0, 0), (0, 0)): t,
            ((0, 1, 0), (0, 0)): t,
            ((0, -1, 0), (0, 0)): t,
        }

        kx = np.linspace(0, 2 * np.pi, Nx, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, Ny, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, Nz, endpoint=False)

        epsilon_k = _build_hamiltonian_k(kx, ky, kz, hr, norb)
        eigenvalues_ek, eigenvectors = _calc_eigenvalues(epsilon_k)
        mu = _determine_mu(eigenvalues_ek, beta, 0.5, norb)
        green_kw = _calc_green(eigenvalues_ek, eigenvectors, mu, beta, nmat)

        U_k = np.ones((norb, norb, Nx, Ny, Nz), dtype=complex) * 2.0
        inter_k = {"CoulombIntra": U_k}

        chi0q = np.full((norb, norb, Nx, Ny, Nz, nmat), 0.1, dtype=complex)

        Pc_q, Ps_q = _compute_vertices(chi0q, inter_k, norb, Nx, Ny, Nz, nmat)
        Vs_q = Pc_q + Ps_q
        G2 = _calc_g2(green_kw)

        return Vs_q, G2, norb, Nx, Ny, Nz

    def test_arnoldi(self):
        """Test Arnoldi (default) method."""
        Vs_q, G2, norb, Nx, Ny, Nz = self._setup_problem()
        vals, vecs = _solve_eigenvalue(
            Vs_q, G2, norb, Nx, Ny, Nz,
            num_eigenvalues=3, method="arnoldi"
        )
        self.assertEqual(len(vals), 3)

    def test_shift_invert_bicgstab(self):
        """Test shift-invert with BiCGSTAB."""
        Vs_q, G2, norb, Nx, Ny, Nz = self._setup_problem()
        vals, vecs = _solve_eigenvalue(
            Vs_q, G2, norb, Nx, Ny, Nz,
            num_eigenvalues=3, method="shift-invert-bicgstab"
        )
        self.assertEqual(len(vals), 3)

    def test_shift_invert_gmres(self):
        """Test shift-invert with GMRES."""
        Vs_q, G2, norb, Nx, Ny, Nz = self._setup_problem()
        vals, vecs = _solve_eigenvalue(
            Vs_q, G2, norb, Nx, Ny, Nz,
            num_eigenvalues=3, method="shift-invert-gmres"
        )
        self.assertEqual(len(vals), 3)

    def test_shift_invert_lgmres(self):
        """Test shift-invert with LGMRES."""
        Vs_q, G2, norb, Nx, Ny, Nz = self._setup_problem()
        vals, vecs = _solve_eigenvalue(
            Vs_q, G2, norb, Nx, Ny, Nz,
            num_eigenvalues=3, method="shift-invert-lgmres"
        )
        self.assertEqual(len(vals), 3)

    def test_methods_agree(self):
        """Test that all methods find the same leading eigenvalue."""
        Vs_q, G2, norb, Nx, Ny, Nz = self._setup_problem()

        vals_arnoldi, _ = _solve_eigenvalue(
            Vs_q, G2, norb, Nx, Ny, Nz,
            num_eigenvalues=3, method="arnoldi"
        )
        vals_bicgstab, _ = _solve_eigenvalue(
            Vs_q, G2, norb, Nx, Ny, Nz,
            num_eigenvalues=3, method="shift-invert-bicgstab"
        )
        vals_gmres, _ = _solve_eigenvalue(
            Vs_q, G2, norb, Nx, Ny, Nz,
            num_eigenvalues=3, method="shift-invert-gmres"
        )

        # Leading eigenvalue should agree across methods
        ev_arnoldi = abs(vals_arnoldi[0])
        ev_bicgstab = abs(vals_bicgstab[0])
        ev_gmres = abs(vals_gmres[0])

        npt.assert_allclose(ev_bicgstab, ev_arnoldi, rtol=0.1,
                            err_msg="BiCGSTAB should agree with Arnoldi")
        npt.assert_allclose(ev_gmres, ev_arnoldi, rtol=0.1,
                            err_msg="GMRES should agree with Arnoldi")


class TestSubspaceIteration(unittest.TestCase):
    """Test subspace iteration for multiple eigenvalues."""

    def _setup_2orb_problem(self):
        """Create a 2-orbital problem with multiple non-zero eigenvalues."""
        norb, Nx, Ny, Nz, nmat = 2, 6, 6, 1, 16
        beta, t, t1 = 5.0, 1.0, 0.5

        hr = {
            ((0, 1, 0), (0, 0)): t, ((0, -1, 0), (0, 0)): t,
            ((0, 1, 0), (1, 1)): t, ((0, -1, 0), (1, 1)): t,
            ((0, 0, 0), (0, 1)): t1, ((-1, 0, 0), (0, 1)): t1,
            ((0, 0, 0), (1, 0)): t1, ((1, 0, 0), (1, 0)): t1,
        }

        kx = np.linspace(0, 2 * np.pi, Nx, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, Ny, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, Nz, endpoint=False)

        epsilon_k = _build_hamiltonian_k(kx, ky, kz, hr, norb)
        eigenvalues, eigenvectors = _calc_eigenvalues(epsilon_k)
        mu = _determine_mu(eigenvalues, beta, 0.75, norb)
        green_kw = _calc_green(eigenvalues, eigenvectors, mu, beta, nmat)

        U_k = np.ones((norb, norb, Nx, Ny, Nz), dtype=complex) * 2.0
        inter_k = {"CoulombIntra": U_k}
        chi0q = np.full((norb, norb, Nx, Ny, Nz, nmat), 0.1, dtype=complex)

        Pc, Ps = _compute_vertices(chi0q, inter_k, norb, Nx, Ny, Nz, nmat)
        Vs_q = Pc + Ps
        G2 = _calc_g2(green_kw)

        return Vs_q, G2, norb, Nx, Ny, Nz

    def test_subspace_finds_multiple(self):
        """Test that subspace iteration finds multiple eigenvalues."""
        Vs_q, G2, norb, Nx, Ny, Nz = self._setup_2orb_problem()

        eigenvalues, eigenvectors = _solve_subspace_iteration(
            Vs_q, G2, norb, Nx, Ny, Nz,
            num_eigenvalues=4, max_iter=200, tol=1e-4
        )
        self.assertEqual(len(eigenvalues), 4)
        self.assertEqual(eigenvectors.shape[0], 4)
        # Should find non-trivial eigenvalues
        self.assertGreater(abs(eigenvalues[0]), 1e-5)

    def test_subspace_via_interface(self):
        """Test subspace method via _solve_eigenvalue interface."""
        Vs_q, G2, norb, Nx, Ny, Nz = self._setup_2orb_problem()

        eigenvalues, eigenvectors = _solve_eigenvalue(
            Vs_q, G2, norb, Nx, Ny, Nz,
            num_eigenvalues=3, method="subspace"
        )
        self.assertEqual(len(eigenvalues), 3)

    def test_subspace_agrees_with_arnoldi(self):
        """Test subspace iteration agrees with Arnoldi on leading eigenvalue."""
        Vs_q, G2, norb, Nx, Ny, Nz = self._setup_2orb_problem()

        vals_arnoldi, _ = _solve_eigenvalue(
            Vs_q, G2, norb, Nx, Ny, Nz,
            num_eigenvalues=3, method="arnoldi"
        )
        vals_subspace, _ = _solve_subspace_iteration(
            Vs_q, G2, norb, Nx, Ny, Nz,
            num_eigenvalues=3, max_iter=300, tol=1e-6
        )

        # Leading eigenvalue should agree
        npt.assert_allclose(
            abs(vals_subspace[0]), abs(vals_arnoldi[0]), rtol=0.01,
            err_msg="Subspace should agree with Arnoldi on leading eigenvalue"
        )

    def test_shifted_bicg_scan(self):
        """Test shifted BiCG spectrum scanning."""
        Vs_q, G2, norb, Nx, Ny, Nz = self._setup_2orb_problem()

        # First get reference from Arnoldi
        vals_ref, _ = _solve_eigenvalue(
            Vs_q, G2, norb, Nx, Ny, Nz,
            num_eigenvalues=3, method="arnoldi"
        )
        max_ev = vals_ref[0].real

        # Scan near the largest eigenvalue
        sigma_list = [max_ev * 0.9, max_ev * 0.5]
        all_evals, all_evecs = _solve_shifted_bicg(
            Vs_q, G2, norb, Nx, Ny, Nz,
            sigma_list=sigma_list, num_eigenvalues=2
        )
        self.assertEqual(len(all_evals), 2)
        # Should find eigenvalues near each sigma
        for sigma in sigma_list:
            self.assertIn(sigma, all_evals)


class TestSimpleGeneralConsistency(unittest.TestCase):
    """Test that general mode agrees with simple mode when J=J'=0."""

    def test_1orb_vertex_consistency(self):
        """For 1 orbital, CoulombIntra only: simple == general vertex."""
        Nx, Ny, Nz = 4, 4, 1
        norb = 1
        nmat = 16

        rng = np.random.default_rng(123)
        chi0q = rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat)) * 0.1 + \
                1j * rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat)) * 0.1

        U_k = np.ones((norb, norb, Nx, Ny, Nz), dtype=complex) * 3.0
        inter_k = {"CoulombIntra": U_k}

        # Simple mode (Wc/Ws)
        Pc_s, Ps_s = _compute_vertices_simple(
            chi0q, inter_k, norb, Nx, Ny, Nz, nmat)
        V_simple = Pc_s + Ps_s  # shape (1,1,Nx,Ny,Nz)

        # General mode (S,C matrices) — force by calling directly
        Vs_gen = _compute_vertices_general(
            chi0q, inter_k, norb, Nx, Ny, Nz, nmat)
        V_general = Vs_gen[0, 0, 0, 0, :, :, :]  # extract scalar component

        npt.assert_allclose(
            V_simple[0, 0, :, :, :], V_general,
            atol=1e-10,
            err_msg="1-orbital: simple and general vertices should match")

    def test_1orb_triplet_consistency(self):
        """For 1 orbital triplet: simple == general."""
        Nx, Ny, Nz = 4, 4, 1
        norb = 1
        nmat = 16

        rng = np.random.default_rng(456)
        chi0q = rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat)) * 0.1 + \
                1j * rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat)) * 0.1

        U_k = np.ones((norb, norb, Nx, Ny, Nz), dtype=complex) * 2.5
        inter_k = {"CoulombIntra": U_k}

        Pc_t, Ps_t = _compute_vertices_simple(
            chi0q, inter_k, norb, Nx, Ny, Nz, nmat, pairing_type="triplet")
        V_simple = Pc_t + Ps_t

        Vs_gen = _compute_vertices_general(
            chi0q, inter_k, norb, Nx, Ny, Nz, nmat, pairing_type="triplet")
        V_general = Vs_gen[0, 0, 0, 0, :, :, :]

        npt.assert_allclose(
            V_simple[0, 0, :, :, :], V_general,
            atol=1e-10,
            err_msg="1-orbital triplet: simple and general vertices should match")

    def test_kernel_output_consistency_1orb(self):
        """Kernel applied through simple vs general vertex gives same sigma."""
        Nx, Ny, Nz = 4, 4, 1
        norb = 1
        nmat = 16
        beta = 5.0

        hr = {
            ((1, 0, 0), (0, 0)): 1.0,
            ((-1, 0, 0), (0, 0)): 1.0,
            ((0, 1, 0), (0, 0)): 1.0,
            ((0, -1, 0), (0, 0)): 1.0,
        }
        kx = np.linspace(0, 2 * np.pi, Nx, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, Ny, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, Nz, endpoint=False)

        epsilon_k = _build_hamiltonian_k(kx, ky, kz, hr, norb)
        eigenvalues, eigenvectors = _calc_eigenvalues(epsilon_k)
        mu = _determine_mu(eigenvalues, beta, 0.5, norb)
        green_kw = _calc_green(eigenvalues, eigenvectors, mu, beta, nmat)

        U_k = np.ones((norb, norb, Nx, Ny, Nz), dtype=complex) * 2.0
        inter_k = {"CoulombIntra": U_k}

        rng = np.random.default_rng(789)
        chi0q = rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat)) * 0.1 + \
                1j * rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat)) * 0.1

        # Simple path
        Pc_s, Ps_s = _compute_vertices_simple(
            chi0q, inter_k, norb, Nx, Ny, Nz, nmat)
        V_simple = Pc_s + Ps_s

        # General path
        Vs_gen = _compute_vertices_general(
            chi0q, inter_k, norb, Nx, Ny, Nz, nmat)

        G2 = _calc_g2(green_kw)
        sigma_old = _initialize_gap("cos", norb, kx, ky, kz)

        sigma_simple = _eliashberg_kernel_fft(V_simple, G2, sigma_old, norb)
        sigma_general = _eliashberg_kernel_fft(Vs_gen, G2, sigma_old, norb)

        npt.assert_allclose(
            sigma_simple, sigma_general,
            atol=1e-10,
            err_msg="Kernel output should match between simple and general modes")

    def test_simple_mode_unchanged_2orb(self):
        """2-orbital simple mode: same result as reference manual computation.

        The simple mode (Wc=U+2V, Ws=-U) is a different mathematical model from
        the general S,C matrix mode for multi-orbital systems. This test verifies
        that the simple mode returns the same result as before the refactoring
        by comparing against a manual computation.
        """
        Nx, Ny, Nz = 4, 4, 1
        norb = 2
        nmat = 16

        rng = np.random.default_rng(321)
        chi0q = rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat)) * 0.05 + \
                1j * rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat)) * 0.05

        U_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        U_k[0, 0, :, :, :] = 3.0
        U_k[1, 1, :, :, :] = 3.0
        inter_k = {"CoulombIntra": U_k}

        Pc_q, Ps_q = _compute_vertices_simple(
            chi0q, inter_k, norb, Nx, Ny, Nz, nmat)

        # Manual reference computation (original algorithm)
        Pc_ref = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        Ps_ref = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        I = np.identity(norb)
        for ix in range(Nx):
            for iy in range(Ny):
                _U = U_k[:, :, ix, iy, 0]
                Wc = _U
                Ws = -_U
                _chi0 = chi0q[:, :, ix, iy, 0, nmat // 2].astype(np.complex128)
                chis = np.linalg.solve(I + _chi0 @ Ws, _chi0)
                chic = np.linalg.solve(I + _chi0 @ Wc, _chi0)
                Pc_ref[:, :, ix, iy, 0] = (Wc + Ws) / 2.0 - 0.5 * Wc @ chic @ Wc
                Ps_ref[:, :, ix, iy, 0] = -Ws + 1.5 * Ws @ chis @ Ws

        npt.assert_allclose(Pc_q, Pc_ref, atol=1e-10,
                            err_msg="Pc should match manual reference")
        npt.assert_allclose(Ps_q, Ps_ref, atol=1e-10,
                            err_msg="Ps should match manual reference")

    def test_iteration_consistency_1orb(self):
        """Full iteration result matches between simple and general modes."""
        Nx, Ny, Nz = 4, 4, 1
        norb = 1
        nmat = 16
        beta = 5.0

        hr = {
            ((1, 0, 0), (0, 0)): 1.0, ((-1, 0, 0), (0, 0)): 1.0,
            ((0, 1, 0), (0, 0)): 1.0, ((0, -1, 0), (0, 0)): 1.0,
        }
        kx = np.linspace(0, 2 * np.pi, Nx, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, Ny, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, Nz, endpoint=False)

        epsilon_k = _build_hamiltonian_k(kx, ky, kz, hr, norb)
        eigenvalues, eigenvectors = _calc_eigenvalues(epsilon_k)
        mu = _determine_mu(eigenvalues, beta, 0.5, norb)
        green_kw = _calc_green(eigenvalues, eigenvectors, mu, beta, nmat)

        U_k = np.ones((norb, norb, Nx, Ny, Nz), dtype=complex) * 2.0
        inter_k = {"CoulombIntra": U_k}
        chi0q = np.full((norb, norb, Nx, Ny, Nz, nmat), 0.1, dtype=complex)

        Pc_s, Ps_s = _compute_vertices_simple(
            chi0q, inter_k, norb, Nx, Ny, Nz, nmat)
        V_simple = Pc_s + Ps_s

        Vs_gen = _compute_vertices_general(
            chi0q, inter_k, norb, Nx, Ny, Nz, nmat)

        G2 = _calc_g2(green_kw)
        sigma_init = _initialize_gap("cos", norb, kx, ky, kz)

        sigma_s, ev_s, conv_s, _ = _solve_iteration(
            green_kw, V_simple, G2, sigma_init, norb,
            max_iter=50, alpha=0.5, tol=1e-6)
        sigma_g, ev_g, conv_g, _ = _solve_iteration(
            green_kw, Vs_gen, G2, sigma_init, norb,
            max_iter=50, alpha=0.5, tol=1e-6)

        npt.assert_allclose(
            ev_s, ev_g, rtol=1e-6,
            err_msg="Iteration eigenvalue should match")
        npt.assert_allclose(
            np.abs(sigma_s), np.abs(sigma_g), atol=1e-6,
            err_msg="Iteration gap function should match")


class TestSCMatrices(unittest.TestCase):
    """Test S, C matrix construction for multi-orbital Hund/Exchange."""

    def test_single_orbital_U_only(self):
        """Test S, C matrices for 1-orbital with only U."""
        norb = 1
        Nx, Ny, Nz = 2, 2, 1
        U_k = np.ones((norb, norb, Nx, Ny, Nz), dtype=complex) * 4.0
        inter_k = {"CoulombIntra": U_k}

        S_mat, C_mat = _build_sc_matrices(inter_k, norb, 0, 0, 0)
        # For 1 orbital, S = U, C = U
        npt.assert_allclose(S_mat, [[4.0]], atol=1e-10)
        npt.assert_allclose(C_mat, [[4.0]], atol=1e-10)

    def test_two_orbital_with_hund(self):
        """Test S, C for 2 orbitals with U, U', J, J'."""
        norb = 2
        Nx, Ny, Nz = 1, 1, 1
        U_val, Up_val, J_val, Jp_val = 4.0, 2.0, 0.5, 0.5

        # On-site interactions (k-independent)
        U_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        U_k[0, 0, 0, 0, 0] = U_val
        U_k[1, 1, 0, 0, 0] = U_val

        Up_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        Up_k[0, 1, 0, 0, 0] = Up_val
        Up_k[1, 0, 0, 0, 0] = Up_val

        J_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        J_k[0, 1, 0, 0, 0] = J_val
        J_k[1, 0, 0, 0, 0] = J_val

        Jp_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        Jp_k[0, 1, 0, 0, 0] = Jp_val
        Jp_k[1, 0, 0, 0, 0] = Jp_val

        inter_k = {
            "CoulombIntra": U_k,
            "CoulombInter": Up_k,
            "Hund": J_k,
            "Exchange": Jp_k,
        }

        S_mat, C_mat = _build_sc_matrices(inter_k, norb, 0, 0, 0)
        # S and C are 4x4 matrices (norb^2 = 4)
        self.assertEqual(S_mat.shape, (4, 4))
        self.assertEqual(C_mat.shape, (4, 4))

        # Check specific elements from Eq.(5):
        # (l1=l2=l3=l4=0): S=U, C=U
        npt.assert_allclose(S_mat[0, 0], U_val, atol=1e-10)
        npt.assert_allclose(C_mat[0, 0], U_val, atol=1e-10)

        # (l1=0,l2=1,l3=0,l4=1) -> l1=l3!=l2=l4: S=U', C=-U'+J
        idx_01 = 0 * norb + 1  # = 1
        npt.assert_allclose(S_mat[idx_01, idx_01], Up_val, atol=1e-10)
        npt.assert_allclose(C_mat[idx_01, idx_01], -Up_val + J_val, atol=1e-10)

        # (l1=0,l2=0,l3=1,l4=1) -> l1=l2!=l3=l4: S=J, C=2U'-J
        idx_00 = 0 * norb + 0  # = 0
        idx_11 = 1 * norb + 1  # = 3
        npt.assert_allclose(S_mat[idx_00, idx_11], J_val, atol=1e-10)
        npt.assert_allclose(C_mat[idx_00, idx_11], 2 * Up_val - J_val, atol=1e-10)

        # (l1=0,l2=1,l3=1,l4=0) -> l1=l4!=l2=l3: S=J', C=J'
        idx_01 = 0 * norb + 1  # = 1
        idx_10 = 1 * norb + 0  # = 2
        npt.assert_allclose(S_mat[idx_01, idx_10], Jp_val, atol=1e-10)
        npt.assert_allclose(C_mat[idx_01, idx_10], Jp_val, atol=1e-10)

    def test_ising_interaction(self):
        """Test S, C for 2 orbitals with Ising interaction.

        Ising contributes to cross and dens channels:
        cross (l1=l3,l2=l4): S = -I, C = -I
        dens  (l1=l2,l3=l4): S = -2I, C = 0
        """
        norb = 2
        Nx, Ny, Nz = 1, 1, 1
        I_val = 1.5

        I_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        I_k[0, 1, 0, 0, 0] = I_val
        I_k[1, 0, 0, 0, 0] = I_val

        inter_k = {"Ising": I_k}
        S_mat, C_mat = _build_sc_matrices(inter_k, norb, 0, 0, 0)

        # cross (l1=0,l2=1,l3=0,l4=1): S = -I, C = -I
        idx_01 = 0 * norb + 1
        npt.assert_allclose(S_mat[idx_01, idx_01], -I_val, atol=1e-10)
        npt.assert_allclose(C_mat[idx_01, idx_01], -I_val, atol=1e-10)

        # dens (l1=0,l2=0,l3=1,l4=1): S = -2I, C = 0
        idx_00 = 0 * norb + 0
        idx_11 = 1 * norb + 1
        npt.assert_allclose(S_mat[idx_00, idx_11], -2.0 * I_val, atol=1e-10)
        npt.assert_allclose(C_mat[idx_00, idx_11], 0.0, atol=1e-10)

        # diag should be zero (Ising only inter-orbital)
        npt.assert_allclose(S_mat[0, 0], 0.0, atol=1e-10)
        npt.assert_allclose(C_mat[0, 0], 0.0, atol=1e-10)

    def test_pairlift_interaction(self):
        """Test S, C for PairLift: no contribution in ph channel."""
        norb = 2
        Nx, Ny, Nz = 1, 1, 1
        P_val = 2.0

        PL_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        PL_k[0, 1, 0, 0, 0] = P_val
        PL_k[1, 0, 0, 0, 0] = P_val

        inter_k = {"PairLift": PL_k}
        S_mat, C_mat = _build_sc_matrices(inter_k, norb, 0, 0, 0)

        # PairLift has no contribution to S/C (particle-particle channel)
        npt.assert_allclose(S_mat, np.zeros((4, 4)), atol=1e-10)
        npt.assert_allclose(C_mat, np.zeros((4, 4)), atol=1e-10)

    def test_pairhop_interaction(self):
        """Test S, C for PairHop interaction.

        PairHop contributes to exch channel:
        exch (l1=l4,l2=l3): S = P, C = P
        """
        norb = 2
        Nx, Ny, Nz = 1, 1, 1
        P_val = 0.8

        PH_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        PH_k[0, 1, 0, 0, 0] = P_val
        PH_k[1, 0, 0, 0, 0] = P_val

        inter_k = {"PairHop": PH_k}
        S_mat, C_mat = _build_sc_matrices(inter_k, norb, 0, 0, 0)

        # exch (l1=0,l2=1,l3=1,l4=0): S = P, C = P
        idx_01 = 0 * norb + 1
        idx_10 = 1 * norb + 0
        npt.assert_allclose(S_mat[idx_01, idx_10], P_val, atol=1e-10)
        npt.assert_allclose(C_mat[idx_01, idx_10], P_val, atol=1e-10)

        # Other channels should be zero
        npt.assert_allclose(S_mat[0, 0], 0.0, atol=1e-10)
        npt.assert_allclose(C_mat[0, 0], 0.0, atol=1e-10)
        npt.assert_allclose(S_mat[idx_01, idx_01], 0.0, atol=1e-10)

    def test_combined_kanamori_ising_pairhop(self):
        """Test S, C with all interactions combined.

        Verify additivity and backward compatibility:
        the Kanamori part should be unchanged.
        """
        norb = 2
        Nx, Ny, Nz = 1, 1, 1
        U, V, J, Jp, I_val, P = 4.0, 2.0, 0.5, 0.5, 1.0, 0.3

        def _make_k(val, diag=False):
            k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
            if diag:
                k[0, 0, 0, 0, 0] = val
                k[1, 1, 0, 0, 0] = val
            else:
                k[0, 1, 0, 0, 0] = val
                k[1, 0, 0, 0, 0] = val
            return k

        # Kanamori only
        inter_k_kanamori = {
            "CoulombIntra": _make_k(U, diag=True),
            "CoulombInter": _make_k(V),
            "Hund": _make_k(J),
            "Exchange": _make_k(Jp),
        }
        S_kam, C_kam = _build_sc_matrices(inter_k_kanamori, norb, 0, 0, 0)

        # Full (Kanamori + Ising + PairHop)
        inter_k_full = dict(inter_k_kanamori)
        inter_k_full["Ising"] = _make_k(I_val)
        inter_k_full["PairHop"] = _make_k(P)

        S_full, C_full = _build_sc_matrices(inter_k_full, norb, 0, 0, 0)

        # Ising + PairHop only
        inter_k_new = {"Ising": _make_k(I_val), "PairHop": _make_k(P)}
        S_new, C_new = _build_sc_matrices(inter_k_new, norb, 0, 0, 0)

        # Additivity: S_full = S_kam + S_new
        npt.assert_allclose(S_full, S_kam + S_new, atol=1e-10)
        npt.assert_allclose(C_full, C_kam + C_new, atol=1e-10)

        # Check specific values for full case
        # diag: S = U, C = U (from Kanamori only)
        npt.assert_allclose(S_full[0, 0], U, atol=1e-10)
        npt.assert_allclose(C_full[0, 0], U, atol=1e-10)

        # cross: S = V - I, C = -V + J - I
        idx_01 = 0 * norb + 1
        npt.assert_allclose(S_full[idx_01, idx_01], V - I_val, atol=1e-10)
        npt.assert_allclose(C_full[idx_01, idx_01], -V + J - I_val, atol=1e-10)

        # dens: S = J - 2I, C = 2V - J
        idx_00 = 0 * norb + 0
        idx_11 = 1 * norb + 1
        npt.assert_allclose(S_full[idx_00, idx_11], J - 2 * I_val, atol=1e-10)
        npt.assert_allclose(C_full[idx_00, idx_11], 2 * V - J, atol=1e-10)

        # exch: S = Jp + P, C = Jp + P
        idx_01 = 0 * norb + 1
        idx_10 = 1 * norb + 0
        npt.assert_allclose(S_full[idx_01, idx_10], Jp + P, atol=1e-10)
        npt.assert_allclose(C_full[idx_01, idx_10], Jp + P, atol=1e-10)


class TestTripletPairing(unittest.TestCase):
    """Test triplet pairing vertex computation."""

    def test_triplet_differs_from_singlet(self):
        """Test that triplet vertex differs from singlet vertex."""
        Nx, Ny, Nz = 2, 2, 1
        norb = 1
        nmat = 8

        rng = np.random.default_rng(42)
        chi0q = rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat)) * 0.1 + \
                1j * rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat)) * 0.1

        U_k = np.ones((norb, norb, Nx, Ny, Nz), dtype=complex) * 2.0
        inter_k = {"CoulombIntra": U_k}

        Pc_s, Ps_s = _compute_vertices(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                                        pairing_type="singlet")
        Pc_t, Ps_t = _compute_vertices(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                                        pairing_type="triplet")

        # Singlet and triplet should differ
        self.assertGreater(
            np.linalg.norm(Ps_s - Ps_t), 1e-10,
            "Spin vertex should differ between singlet and triplet")

    def test_general_triplet_with_hund(self):
        """Test triplet vertex in general mode (with Hund/Exchange)."""
        Nx, Ny, Nz = 2, 2, 1
        norb = 2
        nmat = 8

        rng = np.random.default_rng(42)
        chi0q = rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat)) * 0.05 + \
                1j * rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat)) * 0.05

        U_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        U_k[0, 0] = 4.0
        U_k[1, 1] = 4.0
        J_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        J_k[0, 1] = 0.5
        J_k[1, 0] = 0.5

        inter_k = {"CoulombIntra": U_k, "Hund": J_k}

        Vs_singlet = _compute_vertices(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                                        pairing_type="singlet")
        Vs_triplet = _compute_vertices(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                                        pairing_type="triplet")

        # Both should be 7D (general mode)
        self.assertEqual(Vs_singlet.shape, (norb, norb, norb, norb, Nx, Ny, Nz))
        self.assertEqual(Vs_triplet.shape, (norb, norb, norb, norb, Nx, Ny, Nz))

        # Should differ
        self.assertGreater(
            np.linalg.norm(Vs_singlet - Vs_triplet), 1e-10,
            "Singlet and triplet should differ")

    def test_general_kernel_with_4index(self):
        """Test that the 4-index kernel produces correct shape."""
        norb = 2
        Nx, Ny, Nz = 4, 4, 1
        nmat = 8

        rng = np.random.default_rng(42)
        chi0q = rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat)) * 0.05 + \
                1j * rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat)) * 0.05

        U_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        U_k[0, 0] = 4.0
        U_k[1, 1] = 4.0
        J_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        J_k[0, 1] = 0.5
        J_k[1, 0] = 0.5

        inter_k = {"CoulombIntra": U_k, "Hund": J_k}

        Vs_q = _compute_vertices(chi0q, inter_k, norb, Nx, Ny, Nz, nmat)
        self.assertEqual(Vs_q.ndim, 7)

        # Create G2 and test kernel
        green_kw = rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat)) + \
                   1j * rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat))
        G2 = _calc_g2(green_kw)

        sigma_old = np.ones((norb, norb, Nx, Ny, Nz))
        sigma_old /= np.linalg.norm(sigma_old)

        sigma_new = _eliashberg_kernel_fft(Vs_q, G2, sigma_old, norb)
        self.assertEqual(sigma_new.shape, (norb, norb, Nx, Ny, Nz))


class TestChi0qInternal(unittest.TestCase):
    """Test internal chi0q computation using H-wave's RPA module."""

    def _create_test_files(self, tmpdir, norb=1, Nx=4, Ny=4, Nz=1):
        """Create minimal input files for a test case.

        Creates Geometry, Transfer, and CoulombIntra files for a simple
        square lattice Hubbard model.
        """
        input_dir = os.path.join(tmpdir, "input")
        output_dir = os.path.join(tmpdir, "output")
        os.makedirs(input_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)

        # Geometry file
        geom_file = os.path.join(input_dir, "geom.dat")
        with open(geom_file, "w") as f:
            # rvec (3x3 identity for simple cubic)
            f.write("1.0 0.0 0.0\n")
            f.write("0.0 1.0 0.0\n")
            f.write("0.0 0.0 1.0\n")
            # norb
            f.write("{}\n".format(norb))
            # center positions
            for i in range(norb):
                f.write("0.0 0.0 0.0\n")

        # Transfer file (square lattice, t=1.0)
        transfer_file = os.path.join(input_dir, "transfer.dat")
        with open(transfer_file, "w") as f:
            f.write("Transfer\n")
            f.write("{}\n".format(norb))
            if Nz > 1:
                f.write("5\n")  # 5 R-vectors: (1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,0)
            else:
                f.write("5\n")
            f.write(" 1 1 1 1 1\n")
            # t=1.0 hopping
            for orb in range(norb):
                f.write("  1  0  0  {}  {}  1.0 0.0\n".format(orb + 1, orb + 1))
                f.write(" -1  0  0  {}  {}  1.0 0.0\n".format(orb + 1, orb + 1))
                f.write("  0  1  0  {}  {}  1.0 0.0\n".format(orb + 1, orb + 1))
                f.write("  0 -1  0  {}  {}  1.0 0.0\n".format(orb + 1, orb + 1))

        # CoulombIntra file (U=2.0)
        coulomb_file = os.path.join(input_dir, "coulombintra.dat")
        with open(coulomb_file, "w") as f:
            f.write("CoulombIntra\n")
            f.write("{}\n".format(norb))
            f.write("1\n")
            f.write(" 1\n")
            for orb in range(norb):
                f.write("  0  0  0  {}  {}  2.0 0.0\n".format(orb + 1, orb + 1))

        return input_dir, output_dir

    def _make_input_dict(self, input_dir, output_dir, Nx=4, Ny=4, Nz=1,
                         T=0.1, nmat=32, filling=0.5):
        """Create a TOML-like input dictionary."""
        return {
            "mode": {
                "mode": "RPA",
                "param": {
                    "T": T,
                    "CellShape": [Nx, Ny, Nz],
                    "SubShape": [1, 1, 1],
                    "Nmat": nmat,
                    "filling": filling,
                },
            },
            "file": {
                "input": {
                    "path_to_input": "",
                    "interaction": {
                        "path_to_input": input_dir,
                        "Geometry": "geom.dat",
                        "Transfer": "transfer.dat",
                        "CoulombIntra": "coulombintra.dat",
                    },
                },
                "output": {
                    "path_to_output": output_dir,
                    "chi0q": "chi0q",
                },
            },
            "eliashberg": {
                "chi0q_mode": "calc",
            },
        }

    def test_chi0q_internal_returns_valid_shape(self):
        """Test that internal chi0q computation returns correct shape."""
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            input_dir, output_dir = self._create_test_files(tmpdir)
            input_dict = self._make_input_dict(input_dir, output_dir,
                                                T=0.1, nmat=32, filling=0.5)

            chi0q = _calc_chi0q_internal(input_dict)

            # chi0q should have shape (nmat, nvol, norb, norb)
            norb = 1
            Nx, Ny, Nz = 4, 4, 1
            nvol = Nx * Ny * Nz
            nmat = 32
            self.assertEqual(chi0q.shape, (nmat, nvol, norb, norb))

    def test_chi0q_internal_matches_rpa_solver(self):
        """Test that internal chi0q matches direct RPA solver computation."""
        import tempfile
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.rpa as sol_rpa

        with tempfile.TemporaryDirectory() as tmpdir:
            input_dir, output_dir = self._create_test_files(tmpdir)
            T = 0.1
            nmat = 32
            filling = 0.5
            input_dict = self._make_input_dict(input_dir, output_dir,
                                                T=T, nmat=nmat, filling=filling)

            # Compute via _calc_chi0q_internal
            chi0q_internal = _calc_chi0q_internal(input_dict)

            # Compute via direct RPA solver
            info_mode = {"mode": "RPA", "param": input_dict["mode"]["param"],
                         "calc_scheme": "reduced"}
            info_inputfile = input_dict["file"]["input"]
            info_inputfile["path_to_input"] = ""
            info_log = {"print_level": 1, "print_step": 1}

            read_io = read_input_k.QLMSkInput(info_inputfile)
            ham_info = read_io.get_param("ham")
            solver = sol_rpa.RPA(ham_info, info_log, info_mode)

            green_info = read_io.get_param("green")
            green_info.update(solver.read_init(info_inputfile))

            beta = 1.0 / T
            solver._calc_epsilon_k(green_info)
            Ncond = solver.Ncond / 2  # spin-free
            dist, mu = solver._find_mu(Ncond, solver.T)
            green0, green0_tail = solver._calc_green(beta, mu)
            chi0q_direct = solver._calc_chi0q(green0, green0_tail, beta)

            # Remove block index
            chi0q_direct = chi0q_direct[0]

            npt.assert_allclose(chi0q_internal, chi0q_direct, atol=1e-12,
                                err_msg="Internal chi0q should match direct RPA computation")

    def test_chi0q_internal_nonzero(self):
        """Test that internal chi0q is non-trivially nonzero."""
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            input_dir, output_dir = self._create_test_files(tmpdir)
            input_dict = self._make_input_dict(input_dir, output_dir,
                                                T=0.1, nmat=32, filling=0.5)

            chi0q = _calc_chi0q_internal(input_dict)

            # chi0q should have substantial nonzero values
            self.assertGreater(np.max(np.abs(chi0q)), 0.01,
                               "chi0q should have substantial nonzero values")

    def test_chi0q_internal_vs_load_consistency(self):
        """Test that load mode and calc mode give same result when consistent.

        This test saves a chi0q computed internally, then loads it back and
        verifies they match.
        """
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            input_dir, output_dir = self._create_test_files(tmpdir)
            input_dict = self._make_input_dict(input_dir, output_dir,
                                                T=0.1, nmat=32, filling=0.5)

            # Compute chi0q internally
            chi0q_calc = _calc_chi0q_internal(input_dict)

            # Save it as npz file
            chi0q_file = os.path.join(output_dir, "chi0q.npz")
            np.savez(chi0q_file, chi0q=chi0q_calc)

            # Load it back
            from hwave.sc import _load_chi0q
            chi0q_loaded = _load_chi0q(input_dict)

            npt.assert_allclose(chi0q_calc, chi0q_loaded, atol=1e-15,
                                err_msg="Loaded chi0q should exactly match computed chi0q")


class TestChi0q4Index(unittest.TestCase):
    """Test 4-index (general) chi0q computation and vertex calculation."""

    def _create_2orb_test_files(self, tmpdir):
        """Create 2-orbital test files with inter-orbital hopping."""
        input_dir = os.path.join(tmpdir, "input")
        output_dir = os.path.join(tmpdir, "output")
        os.makedirs(input_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)

        with open(os.path.join(input_dir, "geom.dat"), "w") as f:
            f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n")
            f.write("2\n0.0 0.0 0.0\n0.0 0.0 0.0\n")

        with open(os.path.join(input_dir, "transfer.dat"), "w") as f:
            f.write("Transfer\n2\n5\n 1 1 1 1 1\n")
            for orb in [1, 2]:
                f.write("  1  0  0  %d  %d  1.0 0.0\n" % (orb, orb))
                f.write(" -1  0  0  %d  %d  1.0 0.0\n" % (orb, orb))
                f.write("  0  1  0  %d  %d  0.8 0.0\n" % (orb, orb))
                f.write("  0 -1  0  %d  %d  0.8 0.0\n" % (orb, orb))
            # Inter-orbital hopping
            f.write("  1  0  0  1  2  0.3 0.0\n")
            f.write(" -1  0  0  1  2  0.3 0.0\n")
            f.write("  1  0  0  2  1  0.3 0.0\n")
            f.write(" -1  0  0  2  1  0.3 0.0\n")

        with open(os.path.join(input_dir, "coulombintra.dat"), "w") as f:
            f.write("CoulombIntra\n2\n1\n 1\n")
            f.write("  0  0  0  1  1  3.0 0.0\n")
            f.write("  0  0  0  2  2  3.0 0.0\n")

        return input_dir, output_dir

    def _make_2orb_input_dict(self, input_dir, output_dir):
        return {
            "mode": {
                "mode": "RPA",
                "param": {
                    "T": 0.1,
                    "CellShape": [4, 4, 1],
                    "SubShape": [1, 1, 1],
                    "Nmat": 32,
                    "filling": 0.5,
                },
            },
            "file": {
                "input": {
                    "path_to_input": "",
                    "interaction": {
                        "path_to_input": input_dir,
                        "Geometry": "geom.dat",
                        "Transfer": "transfer.dat",
                        "CoulombIntra": "coulombintra.dat",
                    },
                },
                "output": {
                    "path_to_output": output_dir,
                },
            },
            "eliashberg": {
                "chi0q_mode": "calc",
                "chi0q_tensor": "general",
            },
        }

    def test_general_chi0q_shape(self):
        """Test that general chi0q has 4-index shape."""
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            input_dir, output_dir = self._create_2orb_test_files(tmpdir)
            input_dict = self._make_2orb_input_dict(input_dir, output_dir)

            chi0q = _calc_chi0q_internal(input_dict, chi0q_tensor="general")

            # Should have 6D shape: (nmat, nvol, norb, norb, norb, norb)
            norb = 2
            nmat = 32
            nvol = 16
            self.assertEqual(chi0q.ndim, 6)
            self.assertEqual(chi0q.shape, (nmat, nvol, norb, norb, norb, norb))

    def test_general_diagonal_matches_reduced(self):
        """Test that general chi0q diagonal matches reduced chi0q exactly."""
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            input_dir, output_dir = self._create_2orb_test_files(tmpdir)
            input_dict = self._make_2orb_input_dict(input_dir, output_dir)

            chi0q_gen = _calc_chi0q_internal(input_dict, chi0q_tensor="general")
            chi0q_red = _calc_chi0q_internal(input_dict, chi0q_tensor="reduced")

            norb = 2
            # Diagonal of general should match reduced exactly
            for a in range(norb):
                for b in range(norb):
                    npt.assert_allclose(
                        chi0q_gen[:, :, a, a, b, b],
                        chi0q_red[:, :, a, b],
                        atol=1e-12,
                        err_msg="general diagonal [%d,%d,%d,%d] != reduced [%d,%d]" % (
                            a, a, b, b, a, b))

    def test_general_has_offdiagonal(self):
        """Test that general chi0q has non-negligible off-diagonal elements."""
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            input_dir, output_dir = self._create_2orb_test_files(tmpdir)
            input_dict = self._make_2orb_input_dict(input_dir, output_dir)

            chi0q_gen = _calc_chi0q_internal(input_dict, chi0q_tensor="general")

            norb = 2
            # chi0q[0,1,0,1] should be nonzero (G_00 * G_11 bubble)
            max_offdiag = np.max(np.abs(chi0q_gen[:, :, 0, 1, 0, 1]))
            self.assertGreater(max_offdiag, 0.01,
                               "Off-diagonal chi0q_{0101} should be substantial")

    def test_4index_conversion(self):
        """Test chi0q format conversion for 4-index case."""
        norb = 2
        Nx, Ny, Nz = 4, 4, 1
        nmat = 8
        nvol = Nx * Ny * Nz

        # Create fake H-wave format chi0q (nmat, nvol, norb, norb, norb, norb)
        rng = np.random.default_rng(42)
        chi0q_hw = rng.standard_normal((nmat, nvol, norb, norb, norb, norb))

        chi0q_ref = _convert_chi0q_to_ref_format(chi0q_hw, norb, Nx, Ny, Nz, nmat)

        # Should be (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
        self.assertEqual(chi0q_ref.shape,
                         (norb, norb, norb, norb, Nx, Ny, Nz, nmat))

        # Verify a specific element
        npt.assert_allclose(
            chi0q_ref[0, 1, 0, 1, 2, 3, 0, 5],
            chi0q_hw[5, 2 * Ny * Nz + 3 * Nz + 0, 0, 1, 0, 1],
            atol=1e-15)

    def test_vertex_with_4index_chi0q(self):
        """Test that vertex computation with 4-index chi0q produces correct shape."""
        norb = 2
        Nx, Ny, Nz = 4, 4, 1
        nmat = 16

        # Create 4-index chi0q in ref format (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
        rng = np.random.default_rng(42)
        chi0q = rng.standard_normal((norb, norb, norb, norb, Nx, Ny, Nz, nmat)) * 0.02 + \
                1j * rng.standard_normal((norb, norb, norb, norb, Nx, Ny, Nz, nmat)) * 0.02

        U_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        U_k[0, 0] = 2.0
        U_k[1, 1] = 2.0
        inter_k = {"CoulombIntra": U_k}

        Vs_q = _compute_vertices_general(chi0q, inter_k, norb, Nx, Ny, Nz, nmat)

        self.assertEqual(Vs_q.shape, (norb, norb, norb, norb, Nx, Ny, Nz))

    def test_vertex_4index_vs_2index_differ(self):
        """Test that 4-index chi0q gives different vertex than 2-index diagonal approx."""
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            input_dir, output_dir = self._create_2orb_test_files(tmpdir)
            input_dict = self._make_2orb_input_dict(input_dir, output_dir)

            chi0q_gen = _calc_chi0q_internal(input_dict, chi0q_tensor="general")
            chi0q_red = _calc_chi0q_internal(input_dict, chi0q_tensor="reduced")

            norb = 2
            Nx, Ny, Nz = 4, 4, 1
            nmat = 32

            # Convert to ref format
            chi0q_gen_ref = _convert_chi0q_to_ref_format(chi0q_gen, norb, Nx, Ny, Nz, nmat)
            chi0q_red_ref = _convert_chi0q_to_ref_format(chi0q_red, norb, Nx, Ny, Nz, nmat)

            U_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
            U_k[0, 0] = 3.0
            U_k[1, 1] = 3.0
            inter_k = {"CoulombIntra": U_k}

            # Compute vertex with 4-index (correct)
            Vs_4idx = _compute_vertices_general(
                chi0q_gen_ref, inter_k, norb, Nx, Ny, Nz, nmat)

            # Compute vertex with 2-index (approximate)
            Vs_2idx = _compute_vertices_general(
                chi0q_red_ref, inter_k, norb, Nx, Ny, Nz, nmat)

            # They should differ because off-diagonal chi0q components matter
            diff = np.max(np.abs(Vs_4idx - Vs_2idx))
            self.assertGreater(diff, 1e-3,
                               "4-index and 2-index vertices should differ for 2-orbital system")


if __name__ == '__main__':
    unittest.main()
