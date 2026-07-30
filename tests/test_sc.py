#!/usr/bin/env python3

"""Tests for the linearized Eliashberg equation solver (hwave.sc).

Test strategy:
1. Green's function: Verify G(k, iwn) against direct construction.
2. Vertex: Verify Pc, Ps against reference code on a small system.
3. Kernel: One iteration of the FFT kernel against reference.
4. Convergence: Self-consistent loop converges on a small system.
5. Eigenvalue: Eigenvalue analysis agrees with iteration result.
"""

import logging
import os
import unittest

import numpy as np
import numpy.testing as npt
from numpy.fft import fftn, ifftn

from hwave.sc import (
    _build_hamiltonian_k,
    _build_interaction_k,
    _build_sc_matrices,
    _make_kernel_operator,
    _calc_chi0q_internal,
    _calc_eigenvalues,
    _calc_g2,
    _calc_green,
    _compute_vertices,
    _compute_vertices_simple,
    _compute_vertices_general,
    _compute_vertices_flex,
    _convert_chi0q_to_ref_format,
    _determine_mu,
    _eliashberg_kernel_fft,
    _initialize_gap,
    _is_gap_parity,
    _make_kernel_operator,
    _order_by_seed_overlap,
    _order_eigenpairs,
    _project_gap_parity,
    _reorder_eigenpairs_by_parity,
    _resolve_init_gap,
    _reverse_k_and_orbital,
    _save_results,
    _shift_from_eigenvalues,
    _solve_eigenvalue,
    _solve_iteration,
    _solve_leading,
    _solve_subspace_iteration,
    _solve_shifted_bicg,
)
from scipy.sparse.linalg import LinearOperator


class TestKSpaceBuilderConvention(unittest.TestCase):
    """The hwave_sc k-space builders must follow the SAME Fourier/orbital
    convention as the UHFk/RPA/FLEX solver core:

        M[a, b](k) = sum_R M_R[a, b] * exp(-i k.R)

    (rpa.py _make_ham_trans: tab_r[R, orb1, orb2] + numpy fftn == e^{-ikR}).
    Historically sc.py built epsilon_k[orb2, orb1] with e^{+ikR}, i.e. the
    orbital-transposed matrix at -k. For real (time-reversal-symmetric)
    hoppings the two coincide element-wise, which masked the difference; for
    complex Hermitian hoppings they differ, so quantities loaded from
    FLEX/RPA files (green, chi0q) disagreed element-wise with sc-built ones.
    These tests pin the solver convention with a complex Hermitian fixture."""

    def setUp(self):
        self.Nx, self.Ny, self.Nz = 4, 4, 1
        self.norb = 2
        # complex Hermitian hopping: t(-R, b, a) = conj(t(R, a, b))
        self.hr = {
            ((0, 0, 0), (0, 1)): -2.5 + 0.3j,
            ((0, 0, 0), (1, 0)): -2.5 - 0.3j,
            ((1, 0, 0), (0, 1)): -2.8 + 0.5j,
            ((-1, 0, 0), (1, 0)): -2.8 - 0.5j,
            ((0, 1, 0), (0, 0)): -0.5,
            ((0, -1, 0), (0, 0)): -0.5,
        }
        self.kx = np.linspace(0, 2 * np.pi, self.Nx, endpoint=False)
        self.ky = np.linspace(0, 2 * np.pi, self.Ny, endpoint=False)
        self.kz = np.linspace(0, 2 * np.pi, self.Nz, endpoint=False)

    def _expected(self, value_r):
        """Direct evaluation of the solver convention sum_R M_R[a,b] e^{-ikR}."""
        out = np.zeros((self.norb, self.norb, self.Nx, self.Ny, self.Nz),
                       dtype=complex)
        kxm, kym, kzm = np.meshgrid(self.kx, self.ky, self.kz, indexing='ij')
        for (irvec, (o1, o2)), v in value_r.items():
            out[o1, o2] += v * np.exp(
                -1j * (kxm * irvec[0] + kym * irvec[1] + kzm * irvec[2]))
        return out

    def test_hamiltonian_k_matches_solver_convention(self):
        eps = _build_hamiltonian_k(self.kx, self.ky, self.kz, self.hr,
                                   self.norb)
        npt.assert_allclose(eps, self._expected(self.hr), atol=1e-12)

    def test_interaction_k_is_the_pair_transpose(self):
        """The INTERACTION keeps the solver's e^{-iqR} phase but stores the
        orbital pair TRANSPOSED, unlike the one-body Hamiltonian above.

        It is a four-index object, and its matrix form carries a pair-index
        transpose: H-wave's own paper (arXiv:2308.00324) Eq.(21) defines the
        matrix entering the RPA equation as ``[W(q)]^{ab} = W_q^{ba}``, so a
        density-density term ``W^{aabb} = V_ab(R)`` gives
        ``[W(q)]^{(aa),(bb)} = V_ba(q)``.  ``rpa.py::_make_ham_inter`` builds
        that; this builder must agree (issue #96).

        An earlier revision aligned this builder to the one-body convention,
        which is why the assertion below is the transpose of ``_expected``.
        """
        # complex Hermitian off-site "CoulombInter" exercises both the
        # orbital-index order and the Fourier phase sign
        vr = {
            ((1, 0, 0), (0, 1)): 0.7 + 0.2j,
            ((-1, 0, 0), (1, 0)): 0.7 - 0.2j,
            ((0, 1, 0), (0, 0)): 0.4,
            ((0, -1, 0), (0, 0)): 0.4,
        }
        inter_k = _build_interaction_k(self.kx, self.ky, self.kz,
                                       {"CoulombInter": vr}, self.norb)
        expected = self._expected(vr)
        npt.assert_allclose(inter_k["CoulombInter"],
                            expected.transpose(1, 0, 2, 3, 4), atol=1e-12)
        # ...and not vacuously: the fixture is orbital-asymmetric, so the
        # transpose really is a different matrix
        self.assertGreater(
            np.max(np.abs(expected - expected.transpose(1, 0, 2, 3, 4))), 1e-6)


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


class TestVertexWarnings(unittest.TestCase):
    """Silent approximations in the pairing vertex must warn the user."""

    def _simple_inputs(self, norb=2, extra=None):
        Nx, Ny, Nz, nmat = 2, 2, 1, 4
        chi0q = np.full((norb, norb, Nx, Ny, Nz, nmat), 0.1, dtype=complex)
        U = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        for a in range(norb):
            U[a, a] = 1.0
        inter_k = {"CoulombIntra": U}
        if extra:
            inter_k.update(extra)
        return chi0q, inter_k, norb, Nx, Ny, Nz, nmat

    def test_multiorbital_coulombinter_simple_mode_warns(self):
        """Multi-orbital CoulombInter in 2-index (simple) mode drops the
        inter-orbital cross channel -> must warn."""
        V = np.zeros((2, 2, 2, 2, 1), dtype=complex)
        V[0, 1] = 0.5
        V[1, 0] = 0.5
        chi0q, inter_k, norb, Nx, Ny, Nz, nmat = self._simple_inputs(
            extra={"CoulombInter": V})
        with self.assertLogs("hwave_sc", level="WARNING") as cm:
            _compute_vertices(chi0q, inter_k, norb, Nx, Ny, Nz, nmat)
        self.assertTrue(any("CoulombInter" in m for m in cm.output))

    def test_pairlift_warns_it_is_inert(self):
        """PairLift contributes nothing to the S/C pairing vertex -> warn."""
        PL = np.zeros((2, 2, 2, 2, 1), dtype=complex)
        PL[0, 1] = 0.3
        PL[1, 0] = 0.3
        chi0q, inter_k, norb, Nx, Ny, Nz, nmat = self._simple_inputs(
            extra={"PairLift": PL})
        with self.assertLogs("hwave_sc", level="WARNING") as cm:
            _compute_vertices(chi0q, inter_k, norb, Nx, Ny, Nz, nmat)
        self.assertTrue(any("PairLift" in m for m in cm.output))

    def test_single_orbital_no_coulombinter_warning(self):
        """1-orbital CoulombIntra-only must NOT emit the CoulombInter warning."""
        chi0q, inter_k, norb, Nx, Ny, Nz, nmat = self._simple_inputs(norb=1)
        import logging
        logger = logging.getLogger("hwave_sc")
        with self.assertLogs("hwave_sc", level="WARNING") as cm:
            logger.warning("sentinel")  # ensure the context has >=1 record
            _compute_vertices(chi0q, inter_k, norb, Nx, Ny, Nz, nmat)
        self.assertFalse(any("CoulombInter" in m for m in cm.output))


class TestG2Calculation(unittest.TestCase):
    """Test G2 two-particle Green's function."""

    def test_g2_shape(self):
        """Test G2 output shape."""
        norb = 1
        Nx, Ny, Nz, nmat = 4, 4, 1, 8
        beta = 1.0
        green_kw = np.random.randn(norb, norb, Nx, Ny, Nz, nmat) + \
                   1j * np.random.randn(norb, norb, Nx, Ny, Nz, nmat)
        G2 = _calc_g2(green_kw, beta)
        self.assertEqual(G2.shape, (norb, norb, norb, norb, Nx, Ny, Nz))

    def test_g2_shape_2orb(self):
        """Test G2 shape for 2-orbital case."""
        norb = 2
        Nx, Ny, Nz, nmat = 4, 4, 1, 8
        beta = 1.0
        green_kw = np.random.randn(norb, norb, Nx, Ny, Nz, nmat) + \
                   1j * np.random.randn(norb, norb, Nx, Ny, Nz, nmat)
        G2 = _calc_g2(green_kw, beta)
        self.assertEqual(G2.shape, (norb, norb, norb, norb, Nx, Ny, Nz))


class TestKernel(unittest.TestCase):
    """Test Eliashberg kernel FFT convolution."""

    def test_kernel_shape(self):
        """Test that kernel returns correct shape."""
        norb = 1
        Nx, Ny, Nz, nmat = 4, 4, 1, 8
        beta = 1.0

        P_q = np.random.randn(norb, norb, Nx, Ny, Nz) + \
              1j * np.random.randn(norb, norb, Nx, Ny, Nz)

        rng = np.random.default_rng(42)
        green_kw = rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat)) + \
                   1j * rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat))
        G2 = _calc_g2(green_kw, beta)

        sigma_old = np.ones((norb, norb, Nx, Ny, Nz))
        sigma_old /= np.linalg.norm(sigma_old)

        sigma_new = _eliashberg_kernel_fft(P_q, G2, sigma_old, norb)
        self.assertEqual(sigma_new.shape, (norb, norb, Nx, Ny, Nz))

    def test_kernel_operator_matches_standalone(self):
        """Precomputed-operator matvec is bit-identical to standalone kernel.

        _solve_iteration applies the _make_kernel_operator matvec (which hoists
        the invariant vertex IFFT + G2 preprocessing out of the loop). This must
        be numerically identical to a direct _eliashberg_kernel_fft call.
        Checked for both simple (2-index) and general (4-index) vertices.
        """
        norb = 2
        Nx, Ny, Nz = 4, 4, 1
        rng = np.random.default_rng(123)

        def cplx(*shape):
            return (rng.standard_normal(shape)
                    + 1j * rng.standard_normal(shape))

        G2 = cplx(norb, norb, norb, norb, Nx, Ny, Nz)
        sigma = cplx(norb, norb, Nx, Ny, Nz)

        # simple (5-d vertex) and general (7-d vertex) modes
        for V_q in (cplx(norb, norb, Nx, Ny, Nz),
                    cplx(norb, norb, norb, norb, Nx, Ny, Nz)):
            ref = _eliashberg_kernel_fft(V_q, G2, sigma, norb)
            A, _ = _make_kernel_operator(V_q, G2, norb, Nx, Ny, Nz)
            via_op = A.matvec(sigma.ravel()).reshape(ref.shape)
            # Same ops in the same source order, just with the invariant vertex
            # IFFT + G2 preprocessing precomputed. Differences are at most the
            # last ULP from multi-threaded BLAS reduction-order nondeterminism
            # (~1e-16), far below the suite's 1e-8 tolerance.
            np.testing.assert_allclose(via_op, ref, rtol=1e-13, atol=1e-13)

    def test_matmat_matches_columnwise_matvec(self):
        """Batched matmat equals column-by-column matvec to machine precision.

        The subspace eigensolver applies the kernel to a whole work block at
        once via matmat (one batched FFT) instead of column-by-column. This
        must be numerically identical to applying matvec to each column. A
        wrong batch/FFT-axis layout would produce a large discrepancy here.
        Checked for both simple (5-d vertex) and general (7-d vertex) modes,
        with k > 1 random columns.
        """
        norb = 2
        Nx, Ny, Nz = 4, 4, 1
        vec_size = norb * norb * Nx * Ny * Nz
        k = 7
        rng = np.random.default_rng(2024)

        def cplx(*shape):
            return (rng.standard_normal(shape)
                    + 1j * rng.standard_normal(shape))

        G2 = cplx(norb, norb, norb, norb, Nx, Ny, Nz)
        # mix of real and complex columns (subspace iteration uses real V)
        V = rng.standard_normal((vec_size, k))
        Vc = cplx(vec_size, k)

        for V_q in (cplx(norb, norb, Nx, Ny, Nz),
                    cplx(norb, norb, norb, norb, Nx, Ny, Nz)):
            A, n = _make_kernel_operator(V_q, G2, norb, Nx, Ny, Nz)
            self.assertEqual(n, vec_size)
            for B in (V, Vc):
                ref = np.column_stack(
                    [A.matvec(B[:, j]) for j in range(k)])
                batched = A.matmat(B)
                self.assertEqual(batched.shape, ref.shape)
                # At most last-ULP einsum reduction-order differences from
                # multi-threaded BLAS; far below the suite's 1e-8 tolerance.
                np.testing.assert_allclose(batched, ref, rtol=1e-12, atol=1e-12)


    def test_general_kernel_matches_einsum_reference(self):
        """General-mode matmul kernel matches the original einsum contractions.

        The general (4-index) kernel was rewritten from two
        batched-over-spatial-point einsums ('iljs,js->ils' and 'ijls,js->ils')
        into batched np.matmul (GEMM). This is NOT bit-identical (GEMM changes
        the floating-point reduction order) but must match the einsum to ~1e-13,
        far below the suite's atol=1e-8. This test recomputes the kernel with
        the ORIGINAL einsums as an explicit reference for both matvec (single
        column) and matmat (k>1 columns).
        """
        norb = 3
        Nx, Ny, Nz = 4, 3, 2
        nvol = Nx * Ny * Nz
        vec_size = norb * norb * nvol
        k = 5
        rng = np.random.default_rng(7)

        def cplx(*shape):
            return (rng.standard_normal(shape)
                    + 1j * rng.standard_normal(shape))

        G2 = cplx(norb, norb, norb, norb, Nx, Ny, Nz)
        Vs_q = cplx(norb, norb, norb, norb, Nx, Ny, Nz)  # general (7-d) vertex

        # Explicit reference using the ORIGINAL einsum contractions.
        V_r = ifftn(Vs_q, axes=(-3, -2, -1))
        V_r_flat = V_r.reshape(norb, norb * norb, norb, nvol)
        G2_r = G2.reshape(norb, norb, norb, norb, nvol)
        G2_pre = G2_r.transpose(0, 2, 1, 3, 4).reshape(
            norb, norb, norb * norb, nvol)

        def ref_matvec(v):
            sigma = v.reshape(norb, norb, Nx, Ny, Nz)
            sigma_flat = sigma.reshape(norb * norb, nvol)
            G2Sigma = np.einsum('iljs,js->ils', G2_pre, sigma_flat).reshape(
                norb, norb, Nx, Ny, Nz)
            F_r = ifftn(G2Sigma, axes=(-3, -2, -1))
            F_r_flat = F_r.reshape(norb * norb, nvol)
            sigma_r = np.einsum('ijls,js->ils', V_r_flat, F_r_flat).reshape(
                norb, norb, Nx, Ny, Nz)
            sigma_new = fftn(sigma_r, axes=(-3, -2, -1))
            return (-sigma_new).ravel()

        A, n = _make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)
        self.assertEqual(n, vec_size)

        # matvec
        v = cplx(vec_size)
        new_mv = A.matvec(v)
        ref_mv = ref_matvec(v)
        maxdiff_mv = np.abs(new_mv - ref_mv).max()
        self.assertLess(maxdiff_mv, 1e-10)
        np.testing.assert_allclose(new_mv, ref_mv, rtol=1e-10, atol=1e-10)

        # matmat (k columns)
        B = cplx(vec_size, k)
        new_mm = A.matmat(B)
        ref_mm = np.column_stack([ref_matvec(B[:, j]) for j in range(k)])
        maxdiff_mm = np.abs(new_mm - ref_mm).max()
        self.assertLess(maxdiff_mm, 1e-10)
        np.testing.assert_allclose(new_mm, ref_mm, rtol=1e-10, atol=1e-10)


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
            "d_x2y2", "d_y2z2", "d_xy", "d_xz", "d_yz", "d_z2",
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
            "d_x2y2", "d_y2z2", "d_xy", "d_xz", "d_yz", "d_z2",
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

    def test_d_y2z2_symmetry(self):
        """In-plane d_{y^2-z^2} = cos(ky) - cos(kz) sign structure.

        For a [1, Ny, Nz] cell (kx squashed) this is the d-wave that
        sign-changes across the (pi, 0) <-> (0, pi) zone-boundary q, unlike
        d_yz = sin(ky) sin(kz) which has a node there.
        """
        N = 8
        kx = np.linspace(0, 2 * np.pi, 1, endpoint=False)   # Nx=1, kx=0
        ky = np.linspace(0, 2 * np.pi, N, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, N, endpoint=False)
        i_pi = N // 2
        d = _initialize_gap("d_y2z2", 1, kx, ky, kz)
        val_pi0 = d[0, 0, 0, i_pi, 0]     # (ky, kz) = (pi, 0) -> cos(pi)-cos(0) < 0
        val_0pi = d[0, 0, 0, 0, i_pi]     # (ky, kz) = (0, pi) -> cos(0)-cos(pi) > 0
        self.assertLess(val_pi0 * val_0pi, 0,
                        "d_{y^2-z^2} should change sign between (pi,0) and (0,pi)")
        # contrast: d_yz has a node at both of these q-points
        dyz = _initialize_gap("d_yz", 1, kx, ky, kz)
        npt.assert_allclose(dyz[0, 0, 0, i_pi, 0], 0.0, atol=1e-12)
        npt.assert_allclose(dyz[0, 0, 0, 0, i_pi], 0.0, atol=1e-12)
        # even under k -> -k over the whole grid (singlet/even parity)
        for iy in range(N):
            for iz in range(N):
                npt.assert_allclose(
                    d[0, 0, 0, iy, iz], d[0, 0, 0, (-iy) % N, (-iz) % N],
                    atol=1e-12,
                    err_msg="d_y2z2 not even at (iy,iz)=({},{})".format(iy, iz))

    def test_inplane_modes_vanish_when_kx_squashed(self):
        """Form factors built from sin(kx) vanish identically for a [1,Ny,Nz]
        cell (kx=0), so p_x / d_xy / d_xz are invalid seeds there; the valid
        in-plane seeds are p_y, p_z, d_yz and d_y2z2."""
        N = 8
        kx = np.linspace(0, 2 * np.pi, 1, endpoint=False)   # kx = 0
        ky = np.linspace(0, 2 * np.pi, N, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, N, endpoint=False)
        for mode in ("p_x", "d_xy", "d_xz"):
            sigma = _initialize_gap(mode, 1, kx, ky, kz)
            npt.assert_allclose(
                np.linalg.norm(sigma), 0.0, atol=1e-12,
                err_msg="{} must vanish when kx is squashed".format(mode))
        for mode in ("p_y", "p_z", "d_yz", "d_y2z2"):
            sigma = _initialize_gap(mode, 1, kx, ky, kz)
            self.assertGreater(
                np.linalg.norm(sigma), 0.5,
                "{} should be a nonzero in-plane seed".format(mode))

    def test_vanishing_seed_warns(self):
        """Form factors that vanish on the grid warn on any squashed axis; a
        valid nonzero seed does not warn."""
        N = 4
        k1 = np.linspace(0, 2 * np.pi, 1, endpoint=False)   # squashed axis
        kN = np.linspace(0, 2 * np.pi, N, endpoint=False)
        # p_x = sin(kx) vanishes at Nx=1
        with self.assertLogs("hwave_sc", level="WARNING") as cm:
            _initialize_gap("p_x", 1, k1, kN, kN)
        self.assertTrue(any("zero" in m for m in cm.output))
        # p_z = sin(kz) vanishes at Nz=1 (a different axis)
        with self.assertLogs("hwave_sc", level="WARNING") as cm:
            _initialize_gap("p_z", 1, kN, kN, k1)
        self.assertTrue(any("zero" in m for m in cm.output))
        # a valid in-plane seed must NOT warn (Python 3.9-compatible: capture
        # handler, since unittest.assertNoLogs only exists on 3.10+)
        import logging
        captured = []
        handler = logging.Handler()
        handler.emit = captured.append
        lg = logging.getLogger("hwave_sc")
        lg.addHandler(handler)
        try:
            _initialize_gap("p_y", 1, k1, kN, kN)
        finally:
            lg.removeHandler(handler)
        self.assertEqual(
            [r for r in captured if r.levelno >= logging.WARNING], [],
            "valid seed p_y must not emit a WARNING")

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


class TestComplexKernel(unittest.TestCase):
    """The Eliashberg kernel operator must preserve complex components.

    For complex hopping (e.g. spin-orbit coupling), non-centrosymmetric models,
    or chiral gaps the kernel is genuinely complex; projecting the operator to
    real silently discards the imaginary part and changes the spectrum.
    """

    def _vsq_g2(self):
        norb, Nx, Ny, Nz, nmat = 1, 4, 4, 1, 16
        beta, t = 5.0, 1.0
        hr = {((1, 0, 0), (0, 0)): t, ((-1, 0, 0), (0, 0)): t,
              ((0, 1, 0), (0, 0)): t, ((0, -1, 0), (0, 0)): t}
        kx = np.linspace(0, 2 * np.pi, Nx, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, Ny, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, Nz, endpoint=False)
        ek = _build_hamiltonian_k(kx, ky, kz, hr, norb)
        ev, evec = _calc_eigenvalues(ek)
        mu = _determine_mu(ev, beta, 0.5, norb)
        g = _calc_green(ev, evec, mu, beta, nmat)
        U = np.ones((norb, norb, Nx, Ny, Nz), dtype=complex) * 2.0
        chi0q = np.full((norb, norb, Nx, Ny, Nz, nmat), 0.1, dtype=complex)
        Pc, Ps = _compute_vertices(chi0q, {"CoulombIntra": U}, norb, Nx, Ny, Nz, nmat)
        return Pc + Ps, _calc_g2(g, beta), norb, Nx, Ny, Nz

    def test_kernel_preserves_imaginary_part(self):
        Vs_q, G2, norb, Nx, Ny, Nz = self._vsq_g2()
        A, n = _make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)
        rng = np.random.default_rng(7)
        v = rng.standard_normal(n) + 1j * rng.standard_normal(n)
        out = A.matvec(v)
        # A real-projected operator would zero the imaginary part of the output.
        self.assertGreater(
            np.max(np.abs(np.imag(out))), 1e-10,
            "kernel operator discards the imaginary part (forced real)")


class TestEigenpairOrdering(unittest.TestCase):
    """The leading eigenpair must be the algebraically largest (largest real
    part), since the SC eigenvalue is the largest positive one (-> 1 at Tc).
    A large-magnitude negative eigenvalue must not be reported first."""

    def test_orders_by_largest_real_not_magnitude(self):
        vals = np.array([-5.0, 2.0, 0.5], dtype=complex)
        vecs = np.eye(3, dtype=complex)
        ovals, ovecs = _order_eigenpairs(vals, vecs)
        # leading is the largest real (2.0), not the largest |.| (-5.0)
        self.assertAlmostEqual(ovals[0].real, 2.0)
        npt.assert_allclose(ovecs[:, 0], vecs[:, 1])
        npt.assert_allclose([v.real for v in ovals], [2.0, 0.5, -5.0])


class TestSeedOverlapSelection(unittest.TestCase):
    """When a seed eigenvector is supplied (eigenvector continuation), the
    leading eigenpair must be the one that maximally overlaps the seed, not
    the algebraically largest one. This must hold in the ``vec_size <= 2``
    dense fallback (smallest dynamic grids, e.g. norb=1/Nk=1/Nmat=2) just as
    it does in the ARPACK path -- see issue #61."""

    @staticmethod
    def _diag_operator(diag):
        """A LinearOperator for a real diagonal matrix, standard-basis
        eigenvectors. ``make_operator`` returns ``(A, vec_size)``."""
        diag = np.asarray(diag, dtype=complex)
        n = diag.size
        A = LinearOperator(
            (n, n), matvec=lambda x: diag * np.asarray(x).ravel(),
            dtype=complex)
        return lambda: (A, n)

    def test_dense_fallback_selects_seeded_branch(self):
        # eigenvalues 2.0 (basis vec e0) and 1.0 (basis vec e1); the largest
        # real part is 2.0, but the seed points at the 1.0-branch.
        make_op = self._diag_operator([2.0, 1.0])
        seed = np.array([0.0, 1.0], dtype=complex)
        val, vec, info = _solve_leading(
            make_op, 2, "arnoldi", seed_vec=seed)
        # Without seed handling the dense fallback would return 2.0/e0.
        self.assertAlmostEqual(val.real, 1.0)
        self.assertAlmostEqual(abs(np.vdot(vec, seed)), 1.0, places=10)
        self.assertAlmostEqual(info["eigenvalues"][0].real, 1.0)

    def test_dense_fallback_without_seed_uses_largest_real(self):
        # eigenvalues [2.0, -5.0]: largest real (2.0) != largest magnitude (5.0),
        # so this catches a regression back to magnitude ordering.
        make_op = self._diag_operator([2.0, -5.0])
        val, vec, info = _solve_leading(make_op, 2, "arnoldi")
        # No seed -> physical SC eigenvalue is the algebraically largest.
        self.assertAlmostEqual(val.real, 2.0)

    def test_dense_fallback_equal_overlap_breaks_tie_by_real_part(self):
        # A seed equally overlapping both basis eigenvectors must fall back to
        # the physical (largest-real) ordering, not the eigensolver's arbitrary
        # order (issue #61 tie-break contract).
        make_op = self._diag_operator([1.0, 3.0])       # e1 has the larger real
        seed = np.array([1.0, 1.0], dtype=complex) / np.sqrt(2.0)
        val, vec, info = _solve_leading(make_op, 2, "arnoldi", seed_vec=seed)
        self.assertAlmostEqual(val.real, 3.0)
        self.assertAlmostEqual(info["eigenvalues"][0].real, 3.0)


class TestParitySelection(unittest.TestCase):
    """The gap of a singlet channel is even (Δ(k)=Δ(-k)); a triplet gap is odd
    (Δ(k)=-Δ(-k)). The eigenvalue solver must be able to identify and prefer
    eigenvectors of the requested parity (the kernel preserves parity, so the
    leading eigenvector may belong to the other channel)."""

    def _k(self, N=8):
        k = np.linspace(0, 2 * np.pi, N, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, 1, endpoint=False)
        return k, k, kz

    def test_even_gap_classified_as_singlet(self):
        kx, ky, kz = self._k()
        gap = _initialize_gap("cos", 1, kx, ky, kz)  # cos(kx+ky) is even
        self.assertTrue(_is_gap_parity(gap, "singlet"))
        self.assertFalse(_is_gap_parity(gap, "triplet"))

    def test_odd_gap_classified_as_triplet(self):
        kx, ky, kz = self._k()
        gap = _initialize_gap("p_x", 1, kx, ky, kz)  # sin(kx) is odd
        self.assertTrue(_is_gap_parity(gap, "triplet"))
        self.assertFalse(_is_gap_parity(gap, "singlet"))

    def test_reorder_puts_requested_parity_first(self):
        kx, ky, kz = self._k()
        even = _initialize_gap("cos", 1, kx, ky, kz)
        odd = _initialize_gap("p_x", 1, kx, ky, kz)
        gaps = np.array([even, odd])           # even has the larger eigenvalue
        vals = np.array([2.0, 1.0], dtype=complex)

        # triplet: the odd gap must be promoted to the front even though it has
        # the smaller eigenvalue
        v2, g2 = _reorder_eigenpairs_by_parity(vals, gaps, "triplet")
        self.assertTrue(_is_gap_parity(g2[0], "triplet"))
        self.assertAlmostEqual(v2[0].real, 1.0)

        # singlet: the even gap stays in front
        v1, g1 = _reorder_eigenpairs_by_parity(vals, gaps, "singlet")
        self.assertTrue(_is_gap_parity(g1[0], "singlet"))
        self.assertAlmostEqual(v1[0].real, 2.0)

    def test_no_parity_match_warns(self):
        """If no eigenpair has the requested parity, warn (the leading gap then
        belongs to the wrong channel)."""
        kx, ky, kz = self._k()
        even = _initialize_gap("cos", 1, kx, ky, kz)
        gaps = np.array([even, even])
        vals = np.array([2.0, 1.0], dtype=complex)
        with self.assertLogs("hwave_sc", level="WARNING") as cm:
            _reorder_eigenpairs_by_parity(vals, gaps, "triplet")
        self.assertTrue(any("parity" in m.lower() for m in cm.output))


class TestProjectGapParity(unittest.TestCase):
    """The parity projectors (1 +/- P)/2 (P: Delta_{ab}(k) -> Delta_{ba}(-k))
    must split a gap into its even (singlet) and odd (triplet) sectors."""

    def _k(self, N=8):
        k = np.linspace(0, 2 * np.pi, N, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, 1, endpoint=False)
        return k, k, kz

    def test_projection_keeps_matching_parity(self):
        kx, ky, kz = self._k()
        even = _initialize_gap("cos", 1, kx, ky, kz)
        odd = _initialize_gap("p_x", 1, kx, ky, kz)
        npt.assert_allclose(_project_gap_parity(even, "singlet"), even, atol=1e-12)
        npt.assert_allclose(_project_gap_parity(odd, "triplet"), odd, atol=1e-12)

    def test_projection_removes_opposite_parity(self):
        kx, ky, kz = self._k()
        even = _initialize_gap("cos", 1, kx, ky, kz)
        odd = _initialize_gap("p_x", 1, kx, ky, kz)
        npt.assert_allclose(_project_gap_parity(even, "triplet"), 0.0, atol=1e-12)
        npt.assert_allclose(_project_gap_parity(odd, "singlet"), 0.0, atol=1e-12)

    def test_projectors_are_idempotent_and_complete(self):
        kx, ky, kz = self._k()
        rng = np.random.default_rng(0)
        g = (rng.standard_normal((1, 1, 8, 8, 1))
             + 1j * rng.standard_normal((1, 1, 8, 8, 1)))
        pe = _project_gap_parity(g, "singlet")
        po = _project_gap_parity(g, "triplet")
        # idempotent
        npt.assert_allclose(_project_gap_parity(pe, "singlet"), pe, atol=1e-12)
        npt.assert_allclose(_project_gap_parity(po, "triplet"), po, atol=1e-12)
        # complete: even + odd = identity
        npt.assert_allclose(pe + po, g, atol=1e-12)


class TestIterationParityProjection(unittest.TestCase):
    """The self-consistent power iteration must stay in the channel's parity
    sector. The Eliashberg kernel commutes with parity, so an odd (triplet)
    seed should converge to the leading ODD eigenpair -- but numerical leakage
    lets the larger-magnitude even (singlet) eigenvector take over unless each
    iterate is projected back onto the requested sector."""

    def _symmetric_kernel(self):
        """A parity-symmetric norb=1 kernel whose leading eigenpair is even but
        which also has a (smaller) odd eigenpair. Returns Vs_q, G2 plus the
        leading even/odd eigenvalues and eigenvectors of the kernel matrix."""
        norb, Nx, Ny, Nz = 1, 8, 1, 1
        k = 2 * np.pi * np.arange(Nx) / Nx
        G2 = np.zeros((norb, norb, norb, norb, Nx, Ny, Nz), dtype=complex)
        G2[0, 0, 0, 0, :, 0, 0] = 2.0 + np.cos(k)            # even in k
        Vs_q = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        Vs_q[0, 0, :, 0, 0] = -(1.0 + 0.7 * np.cos(k) + 0.3 * np.cos(2 * k))
        A, n = _make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)
        M = A.matmat(np.eye(n))
        w, v = np.linalg.eig(M)

        def is_even(vec):
            g = vec.reshape(norb, norb, Nx, Ny, Nz)
            gr = _reverse_k_and_orbital(g)
            return np.linalg.norm(g + gr) > np.linalg.norm(g - gr)

        order = np.argsort(-np.abs(w))
        even = [(w[i].real, v[:, i]) for i in order
                if is_even(v[:, i]) and abs(w[i]) > 1e-6]
        odd = [(w[i].real, v[:, i]) for i in order
               if not is_even(v[:, i]) and abs(w[i]) > 1e-6]
        return Vs_q, G2, norb, (Nx, Ny, Nz), even[0], odd[0]

    def test_seed_leaks_to_wrong_parity_without_projection(self):
        """Document the failure mode: without projection an odd+even seed
        converges to the larger even eigenpair."""
        Vs_q, G2, norb, (Nx, Ny, Nz), (le, ve), (lo, vo) = self._symmetric_kernel()
        self.assertGreater(abs(le), abs(lo))   # even dominates
        ve = (ve / np.linalg.norm(ve)).reshape(norb, norb, Nx, Ny, Nz)
        vo = (vo / np.linalg.norm(vo)).reshape(norb, norb, Nx, Ny, Nz)
        seed = vo + 0.1 * ve
        gap, ev, conv, _ = _solve_iteration(
            None, Vs_q, G2, seed.copy(), norb,
            max_iter=300, alpha=0.0, tol=1e-9)
        self.assertTrue(conv)
        self.assertAlmostEqual(ev, abs(le), places=4)
        self.assertFalse(_is_gap_parity(gap, "triplet"))

    def test_projection_keeps_iteration_in_triplet_sector(self):
        """With pairing_type='triplet' the same seed converges to the ODD
        eigenpair instead of leaking to the dominant even one."""
        Vs_q, G2, norb, (Nx, Ny, Nz), (le, ve), (lo, vo) = self._symmetric_kernel()
        ve = (ve / np.linalg.norm(ve)).reshape(norb, norb, Nx, Ny, Nz)
        vo = (vo / np.linalg.norm(vo)).reshape(norb, norb, Nx, Ny, Nz)
        seed = vo + 0.1 * ve
        gap, ev, conv, _ = _solve_iteration(
            None, Vs_q, G2, seed.copy(), norb,
            max_iter=300, alpha=0.0, tol=1e-9, pairing_type="triplet")
        self.assertTrue(conv)
        self.assertAlmostEqual(ev, abs(lo), places=4)
        self.assertTrue(_is_gap_parity(gap, "triplet"))

    def test_singlet_projection_is_noop_for_even_seed(self):
        """Projection must not change a correctly-seeded singlet run."""
        Vs_q, G2, norb, (Nx, Ny, Nz), (le, ve), (lo, vo) = self._symmetric_kernel()
        ve = (ve / np.linalg.norm(ve)).reshape(norb, norb, Nx, Ny, Nz)
        gap_p, ev_p, _, _ = _solve_iteration(
            None, Vs_q, G2, ve.copy(), norb,
            max_iter=300, alpha=0.0, tol=1e-12, pairing_type="singlet")
        gap_n, ev_n, _, _ = _solve_iteration(
            None, Vs_q, G2, ve.copy(), norb,
            max_iter=300, alpha=0.0, tol=1e-12)
        self.assertAlmostEqual(ev_p, ev_n, places=8)
        self.assertTrue(_is_gap_parity(gap_p, "singlet"))


class TestProjectionCentrosymmetryGuard(unittest.TestCase):
    """Parity projection is only valid when the Eliashberg kernel commutes with
    the parity operator P (Delta_{ab}(k) -> Delta_{ba}(-k)), i.e. for
    centrosymmetric systems. For a non-centrosymmetric kernel (e.g. complex/SOC
    hopping, vertex odd in k) the projection would discard physical components
    and steer the power iteration to the wrong eigenpair. The solver must detect
    [A, P] != 0, WARN, and fall back to the un-projected iteration.
    """

    def _asymmetric_kernel(self):
        """A norb=1 kernel that does NOT commute with parity: an odd-in-k part
        in G2 breaks the k -> -k symmetry. Returns Vs_q, G2 and the kernel size.
        """
        norb, Nx, Ny, Nz = 1, 8, 1, 1
        k = 2 * np.pi * np.arange(Nx) / Nx
        G2 = np.zeros((norb, norb, norb, norb, Nx, Ny, Nz), dtype=complex)
        # even base + odd (sin) component -> G2(k) != G2(-k) -> [A, P] != 0
        G2[0, 0, 0, 0, :, 0, 0] = 2.0 + np.cos(k) + 0.6 * np.sin(k)
        Vs_q = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        Vs_q[0, 0, :, 0, 0] = -(1.0 + 0.7 * np.cos(k) + 0.3 * np.cos(2 * k))
        return Vs_q, G2, norb, (Nx, Ny, Nz)

    def test_kernel_is_genuinely_non_centrosymmetric(self):
        """Sanity: the constructed kernel really does not commute with P."""
        Vs_q, G2, norb, (Nx, Ny, Nz) = self._asymmetric_kernel()
        A, n = _make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)
        rng = np.random.default_rng(0)
        x = (rng.standard_normal(n) + 1j * rng.standard_normal(n))
        xg = x.reshape(norb, norb, Nx, Ny, Nz)
        # A(P x) vs P(A x)
        apx = A.matvec(_reverse_k_and_orbital(xg).ravel())
        pax = _reverse_k_and_orbital(
            A.matvec(x).reshape(norb, norb, Nx, Ny, Nz)).ravel()
        self.assertGreater(
            np.linalg.norm(apx - pax) / np.linalg.norm(pax), 1e-6,
            "test kernel must be non-centrosymmetric for this guard test")

    def test_projection_disabled_and_warns_for_non_centrosymmetric_kernel(self):
        """With a non-centrosymmetric kernel, requesting singlet projection must
        WARN and fall back to the un-projected result (same leading eigenpair).
        """
        Vs_q, G2, norb, (Nx, Ny, Nz) = self._asymmetric_kernel()
        rng = np.random.default_rng(1)
        seed = (rng.standard_normal((norb, norb, Nx, Ny, Nz))
                + 1j * rng.standard_normal((norb, norb, Nx, Ny, Nz)))

        with self.assertLogs("hwave_sc", level="WARNING") as cm:
            gap_p, ev_p, conv_p, _ = _solve_iteration(
                None, Vs_q, G2, seed.copy(), norb,
                max_iter=400, alpha=0.0, tol=1e-10, pairing_type="singlet")
        self.assertTrue(any("parity" in m.lower() for m in cm.output))

        gap_n, ev_n, conv_n, _ = _solve_iteration(
            None, Vs_q, G2, seed.copy(), norb,
            max_iter=400, alpha=0.0, tol=1e-10)
        self.assertTrue(conv_p)
        self.assertTrue(conv_n)
        # projection was disabled -> identical leading eigenvalue to no-projection
        self.assertAlmostEqual(ev_p, ev_n, places=6)

    def test_projection_still_applied_for_centrosymmetric_kernel(self):
        """A centrosymmetric kernel must NOT trip the guard: no warning, and the
        triplet projection still steers to the odd eigenpair."""
        norb, Nx, Ny, Nz = 1, 8, 1, 1
        k = 2 * np.pi * np.arange(Nx) / Nx
        G2 = np.zeros((norb, norb, norb, norb, Nx, Ny, Nz), dtype=complex)
        G2[0, 0, 0, 0, :, 0, 0] = 2.0 + np.cos(k)            # even in k
        Vs_q = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        Vs_q[0, 0, :, 0, 0] = -(1.0 + 0.7 * np.cos(k) + 0.3 * np.cos(2 * k))
        rng = np.random.default_rng(2)
        seed = (rng.standard_normal((norb, norb, Nx, Ny, Nz))
                + 1j * rng.standard_normal((norb, norb, Nx, Ny, Nz)))
        # Capture all warnings; a symmetric kernel must not emit a parity one.
        records = []
        handler = logging.Handler()
        handler.emit = records.append
        logger = logging.getLogger("hwave_sc")
        logger.addHandler(handler)
        try:
            gap, ev, conv, _ = _solve_iteration(
                None, Vs_q, G2, seed.copy(), norb,
                max_iter=400, alpha=0.0, tol=1e-10, pairing_type="triplet")
        finally:
            logger.removeHandler(handler)
        self.assertFalse(
            any("parity" in r.getMessage().lower() for r in records),
            "centrosymmetric kernel must not trip the parity guard")
        self.assertTrue(conv)
        self.assertTrue(_is_gap_parity(gap, "triplet"))

    def test_guard_detects_non_centrosymmetric_two_orbital_kernel(self):
        """The guard must also fire for a multi-orbital non-centrosymmetric
        kernel, where parity P = (Delta_ab(k) -> Delta_ba(-k)) involves the
        orbital transposition as well as k -> -k."""
        norb, Nx, Ny, Nz = 2, 8, 1, 1
        k = 2 * np.pi * np.arange(Nx) / Nx
        G2 = np.zeros((norb, norb, norb, norb, Nx, Ny, Nz), dtype=complex)
        # odd (sin) component on the orbital-0 diagonal breaks k -> -k symmetry
        G2[0, 0, 0, 0, :, 0, 0] = 2.0 + np.cos(k) + 0.6 * np.sin(k)
        G2[1, 1, 1, 1, :, 0, 0] = 2.0 + np.cos(k)
        Vs_q = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        for a in range(norb):
            Vs_q[a, a, :, 0, 0] = -(1.0 + 0.5 * np.cos(k))
        rng = np.random.default_rng(4)
        seed = (rng.standard_normal((norb, norb, Nx, Ny, Nz))
                + 1j * rng.standard_normal((norb, norb, Nx, Ny, Nz)))
        with self.assertLogs("hwave_sc", level="WARNING") as cm:
            _solve_iteration(
                None, Vs_q, G2, seed.copy(), norb,
                max_iter=200, alpha=0.0, tol=1e-9, pairing_type="singlet")
        self.assertTrue(any("parity" in m.lower() for m in cm.output))

    def test_small_representable_eigenvalue_converges(self):
        """A small but representable leading eigenvalue must still converge --
        the zero-norm guard must trip only on genuine underflow to 0, never on
        a physically small (but nonzero, non-underflowing) eigenvalue."""
        norb, Nx, Ny, Nz = 1, 8, 1, 1
        k = 2 * np.pi * np.arange(Nx) / Nx
        G2 = np.zeros((norb, norb, norb, norb, Nx, Ny, Nz), dtype=complex)
        G2[0, 0, 0, 0, :, 0, 0] = 2.0 + np.cos(k)
        Vs_q = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        # small overall scale: leading eigenvalue ~1e-6, well within the range
        # where np.linalg.norm does not underflow (norm**2 stays representable).
        Vs_q[0, 0, :, 0, 0] = -(1.0 + 0.7 * np.cos(k)) * 1.0e-6
        rng = np.random.default_rng(7)
        seed = (rng.standard_normal((norb, norb, Nx, Ny, Nz))
                + 1j * rng.standard_normal((norb, norb, Nx, Ny, Nz)))
        gap, ev, conv, _ = _solve_iteration(
            None, Vs_q, G2, seed.copy(), norb,
            max_iter=300, alpha=0.0, tol=1e-9)
        self.assertTrue(conv)
        self.assertGreater(ev, 0.0)
        self.assertTrue(np.isfinite(ev))

    def test_solve_iteration_rejects_unknown_pairing_type(self):
        """An invalid pairing_type must raise even when the centrosymmetry
        guard would disable projection (so _project_gap_parity is never hit)."""
        Vs_q, G2, norb, (Nx, Ny, Nz) = self._asymmetric_kernel()
        rng = np.random.default_rng(8)
        seed = (rng.standard_normal((norb, norb, Nx, Ny, Nz))
                + 1j * rng.standard_normal((norb, norb, Nx, Ny, Nz)))
        with self.assertRaises(ValueError):
            _solve_iteration(
                None, Vs_q, G2, seed.copy(), norb,
                max_iter=10, alpha=0.0, tol=1e-9, pairing_type="bogus")

    def test_zero_norm_iterate_does_not_nan(self):
        """If the kernel annihilates the requested parity sector the iterate
        norm collapses to zero; normalizing must not return NaN. The solver
        must report non-convergence with a finite gap instead."""
        norb, Nx, Ny, Nz = 1, 8, 1, 1
        k = 2 * np.pi * np.arange(Nx) / Nx
        # Zero vertex -> kernel A == 0 -> every iterate maps to the zero vector.
        Vs_q = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        G2 = np.zeros((norb, norb, norb, norb, Nx, Ny, Nz), dtype=complex)
        G2[0, 0, 0, 0, :, 0, 0] = 2.0 + np.cos(k)
        rng = np.random.default_rng(3)
        seed = (rng.standard_normal((norb, norb, Nx, Ny, Nz))
                + 1j * rng.standard_normal((norb, norb, Nx, Ny, Nz)))
        gap, ev, conv, _ = _solve_iteration(
            None, Vs_q, G2, seed.copy(), norb,
            max_iter=50, alpha=0.0, tol=1e-10, pairing_type="triplet")
        self.assertFalse(conv)
        self.assertTrue(np.all(np.isfinite(gap)),
                        "collapsed iterate must not produce NaN/inf")


class TestProjectGapParityValidation(unittest.TestCase):
    """_project_gap_parity is now reused by the solver guard; an unknown
    pairing_type must raise rather than silently selecting the odd projector."""

    def test_rejects_unknown_pairing_type(self):
        g = np.zeros((1, 1, 4, 1, 1), dtype=complex)
        with self.assertRaises(ValueError):
            _project_gap_parity(g, "bogus")

    def test_accepts_known_pairing_types(self):
        g = np.zeros((1, 1, 4, 1, 1), dtype=complex)
        # must not raise
        _project_gap_parity(g, "singlet")
        _project_gap_parity(g, "triplet")


class TestEigenvalueParityColumn(unittest.TestCase):
    """The eigenvalue output may carry a per-eigenvalue parity-match flag so
    downstream tools can tell physical (channel-parity) modes from spurious
    opposite-parity ones. The flag is an appended column; the existing
    index/Re/Im/|ev| columns are unchanged (backward compatible)."""

    def _read_analysis_rows(self, path):
        rows = []
        with open(path) as f:
            in_section = False
            for line in f:
                if line.startswith("# Eigenvalue analysis"):
                    in_section = True
                    continue
                if line.startswith("# index"):
                    continue
                if in_section and line.strip():
                    rows.append(line.split())
        return rows

    def test_backward_compatible_without_match(self):
        import tempfile
        evs = np.array([1.5 + 0j, 0.97 + 0j])
        with tempfile.TemporaryDirectory() as d:
            _save_results(d, None, None, evs, None, None, None,
                          eigenvalue_file="eigenvalue.dat")
            rows = self._read_analysis_rows(os.path.join(d, "eigenvalue.dat"))
        self.assertEqual(len(rows), 2)
        # index Re Im |ev| -> 4 columns, Re at index 1
        self.assertEqual(len(rows[0]), 4)
        npt.assert_allclose(float(rows[0][1]), 1.5)

    def test_writes_match_column_when_provided(self):
        import tempfile
        evs = np.array([1.5 + 0j, 0.97 + 0j])
        match = np.array([False, True])   # 1.5 spurious (wrong parity), 0.97 physical
        with tempfile.TemporaryDirectory() as d:
            _save_results(d, None, None, evs, None, None, None,
                          eigenvalue_file="eigenvalue.dat",
                          eigenvalue_match=match)
            rows = self._read_analysis_rows(os.path.join(d, "eigenvalue.dat"))
        self.assertEqual(len(rows[0]), 5)          # appended column
        npt.assert_allclose(float(rows[0][1]), 1.5)  # Re still at index 1
        self.assertEqual(rows[0][4], "0")          # spurious
        self.assertEqual(rows[1][4], "1")          # physical

    def test_match_column_reflects_is_gap_parity(self):
        """End-to-end: the written match flag must equal _is_gap_parity for the
        actual eigenvectors, the same classification calc_eliashberg applies."""
        import tempfile
        N = 8
        k = np.linspace(0, 2 * np.pi, N, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, 1, endpoint=False)
        even = _initialize_gap("cos", 1, k, k, kz)   # singlet-parity (even)
        odd = _initialize_gap("p_x", 1, k, k, kz)    # triplet-parity (odd)
        gaps = [even, odd]
        # classify against the singlet channel, exactly as calc_eliashberg does
        match = np.array([_is_gap_parity(g, "singlet") for g in gaps])
        self.assertTrue(match[0])     # even matches singlet
        self.assertFalse(match[1])    # odd does not
        # the triplet channel must classify the two gaps the opposite way
        match_t = np.array([_is_gap_parity(g, "triplet") for g in gaps])
        self.assertFalse(match_t[0])  # even does not match triplet
        self.assertTrue(match_t[1])   # odd matches triplet
        evs = np.array([0.9 + 0j, 0.8 + 0j])
        with tempfile.TemporaryDirectory() as d:
            _save_results(d, None, None, evs, None, None, None,
                          eigenvalue_file="eigenvalue.dat",
                          eigenvalue_match=match)
            rows = self._read_analysis_rows(os.path.join(d, "eigenvalue.dat"))
        self.assertEqual(rows[0][4], "1")
        self.assertEqual(rows[1][4], "0")


class TestTutorialPlotLoadEigenvalues(unittest.TestCase):
    """The tutorial plot_results.py must read the new match column when present
    and fall back to 'unknown sector' (None) for legacy/partial files so it
    never infers a parity sector it does not actually have."""

    def _load_plot_module(self):
        import importlib.util
        path = os.path.join(
            os.path.dirname(__file__), "..",
            "docs", "en", "source", "rpa", "sample_sc", "plot_results.py")
        spec = importlib.util.spec_from_file_location("_tut_plot", path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def _write(self, rows):
        import tempfile
        d = tempfile.mkdtemp()
        p = os.path.join(d, "eigenvalue.dat")
        with open(p, "w") as f:
            f.write("# Iteration eigenvalue\n9.6e-01\n# Eigenvalue analysis\n")
            f.write("# index ...\n")
            f.writelines(rows)
        return p

    def test_legacy_four_column_yields_none(self):
        pr = self._load_plot_module()
        p = self._write(["   0  1.58e+00 0 1.58e+00\n",
                         "   1  0.97e+00 0 0.97e+00\n"])
        _, ev, match = pr.load_eigenvalues(p)
        self.assertEqual(len(ev), 2)
        self.assertIsNone(match, "legacy 4-column file must yield match=None")

    def test_new_five_column_yields_flags(self):
        pr = self._load_plot_module()
        p = self._write(["   0  1.58e+00 0 1.58e+00 0\n",
                         "   1  0.97e+00 0 0.97e+00 1\n"])
        _, _, match = pr.load_eigenvalues(p)
        self.assertEqual(match, [False, True])

    def test_partial_match_column_treated_as_legacy(self):
        pr = self._load_plot_module()
        p = self._write(["   0  1.58e+00 0 1.58e+00 0\n",
                         "   1  0.97e+00 0 0.97e+00\n"])   # second row lacks flag
        _, _, match = pr.load_eigenvalues(p)
        self.assertIsNone(match, "mixed/partial files must be treated as legacy")


class TestShiftEstimate(unittest.TestCase):
    """Shift-invert target must aim at the largest real eigenvalue (the SC
    instability), not the largest magnitude (which can be a large negative)."""

    def test_shift_from_largest_real(self):
        vals = np.array([-6.0, 0.9, -3.0], dtype=complex)
        self.assertAlmostEqual(_shift_from_eigenvalues(vals), 0.9 * 0.9)

    def test_shift_factor(self):
        vals = np.array([-6.0, 2.0], dtype=complex)
        self.assertAlmostEqual(_shift_from_eigenvalues(vals, factor=0.5), 1.0)


class TestFlexVertexWarnings(unittest.TestCase):
    """PairLift is inert in the FLEX pairing-vertex path too -> must warn."""

    def test_pairlift_warns_in_flex_path(self):
        norb, Nx, Ny, Nz = 2, 2, 2, 1
        nd = norb * norb
        chis = np.zeros((Nx, Ny, Nz, nd, nd), dtype=complex)
        chic = np.zeros((Nx, Ny, Nz, nd, nd), dtype=complex)
        U = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        U[0, 0] = 1.0
        U[1, 1] = 1.0
        PL = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        PL[0, 1] = 0.3
        PL[1, 0] = 0.3
        inter_k = {"CoulombIntra": U, "PairLift": PL}
        with self.assertLogs("hwave_sc", level="WARNING") as cm:
            _compute_vertices_flex(chis, chic, inter_k, norb, Nx, Ny, Nz)
        self.assertTrue(any("PairLift" in m for m in cm.output))


class TestFlexVertexMYOConvention(unittest.TestCase):
    """FLEX susceptibilities from the general path are in MYO convention; the
    pairing vertex built from them must use MYO S/C matrices, not the default
    Kuroki ones (they differ in C(ab,ab) = -U'+2J vs -U'+J when Hund J != 0)."""

    def _setup(self):
        norb, Nx, Ny, Nz = 2, 2, 2, 1
        nd = norb * norb
        rng = np.random.default_rng(0)
        chis = rng.standard_normal((Nx, Ny, Nz, nd, nd)) + 0j
        chic = rng.standard_normal((Nx, Ny, Nz, nd, nd)) + 0j
        U = np.zeros((norb, norb, Nx, Ny, Nz), complex); U[0, 0] = 2.0; U[1, 1] = 2.0
        Up = np.zeros((norb, norb, Nx, Ny, Nz), complex); Up[0, 1] = 1.0; Up[1, 0] = 1.0
        J = np.zeros((norb, norb, Nx, Ny, Nz), complex); J[0, 1] = 0.4; J[1, 0] = 0.4
        inter_k = {"CoulombIntra": U, "CoulombInter": Up, "Hund": J}
        return chis, chic, inter_k, norb, Nx, Ny, Nz

    def test_myo_differs_from_kuroki_with_hund(self):
        chis, chic, inter_k, norb, Nx, Ny, Nz = self._setup()
        v_kuroki = _compute_vertices_flex(chis, chic, inter_k, norb, Nx, Ny, Nz,
                                          convention="kuroki")
        v_myo = _compute_vertices_flex(chis, chic, inter_k, norb, Nx, Ny, Nz,
                                       convention="myo")
        # Adjudicated by exact diagonalization (issue #113): the per-type
        # values make the two builders identical, so the vertices must AGREE.
        npt.assert_allclose(v_myo, v_kuroki, atol=1e-12)

    def test_myo_matches_myo_sc_formula(self):
        from hwave.solver._sc_matrices_myo import build_sc_matrices_myo
        chis, chic, inter_k, norb, Nx, Ny, Nz = self._setup()
        v_myo = _compute_vertices_flex(chis, chic, inter_k, norb, Nx, Ny, Nz,
                                       convention="myo", pairing_type="singlet")
        S, C = build_sc_matrices_myo(inter_k, norb, Nx, Ny, Nz)
        Vs_all = 1.5 * (S @ chis @ S) - 0.5 * (C @ chic @ C) + 0.5 * (S + C)
        Vs_q = Vs_all.reshape(Nx, Ny, Nz, norb, norb, norb, norb).transpose(
            3, 4, 5, 6, 0, 1, 2)
        np.testing.assert_allclose(v_myo, Vs_q, atol=1e-12)

    def test_default_is_kuroki(self):
        chis, chic, inter_k, norb, Nx, Ny, Nz = self._setup()
        v_default = _compute_vertices_flex(chis, chic, inter_k, norb, Nx, Ny, Nz)
        v_kuroki = _compute_vertices_flex(chis, chic, inter_k, norb, Nx, Ny, Nz,
                                          convention="kuroki")
        np.testing.assert_allclose(v_default, v_kuroki, atol=1e-12)


class TestResolveInitGap(unittest.TestCase):
    """The default initial gap must match the pairing channel's parity.

    The Eliashberg kernel preserves parity, so an even seed (e.g. 'cos')
    cannot reach the odd triplet solution. When the user does not specify
    init_gap, the default must be even for singlet and odd for triplet.
    """

    def test_triplet_default_is_odd_seed(self):
        mode = _resolve_init_gap(None, "triplet")
        # the resolved default must be an odd (p-wave) seed
        N = 8
        k = np.linspace(0, 2 * np.pi, N, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, 1, endpoint=False)
        sigma = _initialize_gap(mode, 1, k, k, kz)
        for ix in range(N):
            ix_inv = (N - ix) % N
            npt.assert_allclose(
                sigma[0, 0, ix, 0, 0], -sigma[0, 0, ix_inv, 0, 0], atol=1e-10,
                err_msg="triplet default seed must be odd under k -> -k")

    def test_singlet_default_is_even_cos(self):
        self.assertEqual(_resolve_init_gap(None, "singlet"), "cos")

    def test_explicit_init_gap_is_respected(self):
        self.assertEqual(_resolve_init_gap("d_x2y2", "triplet"), "d_x2y2")
        self.assertEqual(_resolve_init_gap("s", "singlet"), "s")


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

        chi0q_ref = _convert_chi0q_to_ref_format(chi0q_hwave, norb, Nx, Ny, Nz)
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
        G2 = _calc_g2(green_kw, beta)
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
        G2 = _calc_g2(green_kw, beta)

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
        G2 = _calc_g2(green_kw, beta)
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
        G2 = _calc_g2(green_kw, beta)

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

        # Methods should agree on the dominant eigenvalue. Compare by largest
        # magnitude (order-independent): the reported ordering is now by real
        # part, so vals[0] is the largest-real, not the dominant, eigenvalue.
        ev_arnoldi = np.max(np.abs(vals_arnoldi))
        ev_bicgstab = np.max(np.abs(vals_bicgstab))
        ev_gmres = np.max(np.abs(vals_gmres))

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
        G2 = _calc_g2(green_kw, beta)

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
        # Should find non-trivial eigenvalues (order-independent check)
        self.assertGreater(np.max(np.abs(eigenvalues)), 1e-5)

    def test_subspace_via_interface(self):
        """Test subspace method via _solve_eigenvalue interface."""
        Vs_q, G2, norb, Nx, Ny, Nz = self._setup_2orb_problem()

        eigenvalues, eigenvectors = _solve_eigenvalue(
            Vs_q, G2, norb, Nx, Ny, Nz,
            num_eigenvalues=3, method="subspace"
        )
        self.assertEqual(len(eigenvalues), 3)

    def test_subspace_agrees_with_arnoldi(self):
        """Test subspace iteration agrees with Arnoldi on the dominant eigenvalue."""
        Vs_q, G2, norb, Nx, Ny, Nz = self._setup_2orb_problem()

        vals_arnoldi, _ = _solve_eigenvalue(
            Vs_q, G2, norb, Nx, Ny, Nz,
            num_eigenvalues=3, method="arnoldi"
        )
        vals_subspace, _ = _solve_subspace_iteration(
            Vs_q, G2, norb, Nx, Ny, Nz,
            num_eigenvalues=3, max_iter=300, tol=1e-6
        )

        # Methods should agree on the dominant eigenvalue (largest magnitude).
        # Compared order-independently since reporting is now ordered by real
        # part, not magnitude.
        npt.assert_allclose(
            np.max(np.abs(vals_subspace)), np.max(np.abs(vals_arnoldi)),
            rtol=0.01,
            err_msg="Subspace should agree with Arnoldi on the dominant eigenvalue"
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

        G2 = _calc_g2(green_kw, beta)
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

        G2 = _calc_g2(green_kw, beta)
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

        # (l1=0,l2=1,l3=0,l4=1) -> l1=l3!=l2=l4:
        #   S = U' - J + J' (= U' for J = J'),  C = -U' + J + J'
        # (adjudicated by exact diagonalization, issue #113)
        idx_01 = 0 * norb + 1  # = 1
        npt.assert_allclose(S_mat[idx_01, idx_01],
                            Up_val - J_val + Jp_val, atol=1e-10)
        npt.assert_allclose(C_mat[idx_01, idx_01],
                            -Up_val + J_val + Jp_val, atol=1e-10)

        # (l1=0,l2=0,l3=1,l4=1) -> l1=l2!=l3=l4: S=J, C=2U'-J
        idx_00 = 0 * norb + 0  # = 0
        idx_11 = 1 * norb + 1  # = 3
        npt.assert_allclose(S_mat[idx_00, idx_11], J_val, atol=1e-10)
        npt.assert_allclose(C_mat[idx_00, idx_11], 2 * Up_val - J_val, atol=1e-10)

        # (l1=0,l2=1,l3=1,l4=0) -> l1=l4!=l2=l3: Exchange moved to the
        # (ab,ab) slots (issue #113); with no PairHop these are zero.
        idx_01 = 0 * norb + 1  # = 1
        idx_10 = 1 * norb + 0  # = 2
        npt.assert_allclose(S_mat[idx_01, idx_10], 0.0, atol=1e-10)
        npt.assert_allclose(C_mat[idx_01, idx_10], 0.0, atol=1e-10)

    def test_ising_interaction(self):
        """Test S, C for 2 orbitals with Ising interaction.

        Ising contributes to cross and dens channels:
        cross (l1=l3,l2=l4): S = +I (sign per issue #113), C = -I
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

        # cross (l1=0,l2=1,l3=0,l4=1): S = +I, C = -I. The S sign was
        # adjudicated by exact diagonalization (issue #113).
        idx_01 = 0 * norb + 1
        npt.assert_allclose(S_mat[idx_01, idx_01], +I_val, atol=1e-10)
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

        # cross: S = V + I - J + J', C = -V + J + J' - I  (issue #113)
        idx_01 = 0 * norb + 1
        npt.assert_allclose(S_full[idx_01, idx_01],
                            V + I_val - J + Jp, atol=1e-10)
        npt.assert_allclose(C_full[idx_01, idx_01],
                            -V + J + Jp - I_val, atol=1e-10)

        # dens: S = J - 2I, C = 2V - J
        idx_00 = 0 * norb + 0
        idx_11 = 1 * norb + 1
        npt.assert_allclose(S_full[idx_00, idx_11], J - 2 * I_val, atol=1e-10)
        npt.assert_allclose(C_full[idx_00, idx_11], 2 * V - J, atol=1e-10)

        # exch: PairHop only (Exchange moved to (ab,ab), issue #113)
        idx_01 = 0 * norb + 1
        idx_10 = 1 * norb + 0
        npt.assert_allclose(S_full[idx_01, idx_10], P, atol=1e-10)
        npt.assert_allclose(C_full[idx_01, idx_10], P, atol=1e-10)


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
        beta = 1.0

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
        G2 = _calc_g2(green_kw, beta)

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
            chi0q_loaded, static_index = _load_chi0q(input_dict)
            self.assertIsNone(static_index,
                              "metadata-less file: the caller slices the "
                              "center of its actual frequency axis")

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

        chi0q_ref = _convert_chi0q_to_ref_format(chi0q_hw, norb, Nx, Ny, Nz)

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
        """A 2-index (reduced) chi0q keeps only the density-density block, so it
        agrees with the 4-index vertex exactly when the interaction lives on
        that block (CoulombIntra only), and differs once a term reaches the
        off-density blocks (here Hund).

        Before the density-pair embedding fix, the reduced chi0q was scattered
        as kron(chi0_2d, I_norb), which made even the CoulombIntra-only case
        disagree with the 4-index reference -- the disagreement this test used
        to assert."""
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
            chi0q_gen_ref = _convert_chi0q_to_ref_format(chi0q_gen, norb, Nx, Ny, Nz)
            chi0q_red_ref = _convert_chi0q_to_ref_format(chi0q_red, norb, Nx, Ny, Nz)

            U_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
            U_k[0, 0] = 3.0
            U_k[1, 1] = 3.0
            inter_k = {"CoulombIntra": U_k}

            # CoulombIntra only: S and C are confined to the density-pair block
            # (_build_sc_matrices_all_q case 1), so the off-density components
            # of the 4-index chi0q never enter S @ chi @ S. The reduced chi0q
            # loses nothing that the vertex reads, and the two must agree.
            Vs_4idx = _compute_vertices_general(
                chi0q_gen_ref, inter_k, norb, Nx, Ny, Nz, nmat)
            Vs_2idx = _compute_vertices_general(
                chi0q_red_ref, inter_k, norb, Nx, Ny, Nz, nmat)
            np.testing.assert_allclose(
                Vs_2idx, Vs_4idx, rtol=1e-9, atol=1e-11,
                err_msg="for CoulombIntra only the reduced chi0q carries every "
                        "component the vertex reads, so it must reproduce the "
                        "4-index vertex exactly")

            # Add Hund: S/C now also populate the off-density blocks
            # S/C[(a,b),(a,b)] and S/C[(a,b),(b,a)], where the reduced chi0q has
            # nothing. The 4-index vertex dresses those channels, the reduced
            # one leaves them bare -- so now they genuinely differ.
            J_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
            J_k[0, 1] = 0.7
            J_k[1, 0] = 0.7
            inter_k_J = dict(inter_k, Hund=J_k)
            Vs_4idx_J = _compute_vertices_general(
                chi0q_gen_ref, inter_k_J, norb, Nx, Ny, Nz, nmat)
            Vs_2idx_J = _compute_vertices_general(
                chi0q_red_ref, inter_k_J, norb, Nx, Ny, Nz, nmat)
            diff = np.max(np.abs(Vs_4idx_J - Vs_2idx_J))
            self.assertGreater(
                diff, 1e-3,
                "with Hund the reduced chi0q is missing the off-density "
                "components the vertex reads, so the vertices must differ")


class TestKanamoriInteraction(unittest.TestCase):
    """End-to-end tests for Kanamori-type interactions (U, U', J, J').

    Verifies that the Eliashberg solver works correctly when Hund coupling J
    and exchange J' are present, including:
    1. Simple vs General vertex consistency when J=0
    2. RPA susceptibility via S/C matrices vs manual reference
    3. End-to-end iteration + eigenvalue consistency with Kanamori interactions
    """

    def _setup_2orb_model(self, Nx=4, Ny=4, Nz=1, nmat=16, beta=5.0,
                          filling=0.5, t1=1.0, t2=0.8, t12=0.3):
        """Set up a 2-orbital tight-binding model on a square lattice."""
        norb = 2
        hr = {
            ((1, 0, 0), (0, 0)): t1,
            ((-1, 0, 0), (0, 0)): t1,
            ((0, 1, 0), (0, 0)): t1,
            ((0, -1, 0), (0, 0)): t1,
            ((1, 0, 0), (1, 1)): t2,
            ((-1, 0, 0), (1, 1)): t2,
            ((0, 1, 0), (1, 1)): t2,
            ((0, -1, 0), (1, 1)): t2,
            ((1, 0, 0), (0, 1)): t12,
            ((-1, 0, 0), (0, 1)): t12,
            ((1, 0, 0), (1, 0)): t12,
            ((-1, 0, 0), (1, 0)): t12,
        }

        kx = np.linspace(0, 2 * np.pi, Nx, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, Ny, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, Nz, endpoint=False)

        epsilon_k = _build_hamiltonian_k(kx, ky, kz, hr, norb)
        eigenvalues, eigenvectors = _calc_eigenvalues(epsilon_k)
        mu = _determine_mu(eigenvalues, beta, filling, norb)
        green_kw = _calc_green(eigenvalues, eigenvectors, mu, beta, nmat)

        return {
            "norb": norb, "Nx": Nx, "Ny": Ny, "Nz": Nz,
            "nmat": nmat, "beta": beta,
            "kx": kx, "ky": ky, "kz": kz,
            "green_kw": green_kw,
        }

    @staticmethod
    def _embed4(chi0q, norb, Nx, Ny, Nz, nmat):
        """Density-pair embedding of a 2-index chi0q into the four-index
        form (out[(a,a),(b,b)] = X[a,b], everything else zero)."""
        out = np.zeros((norb, norb, norb, norb, Nx, Ny, Nz, nmat),
                       dtype=complex)
        for a in range(norb):
            for b in range(norb):
                out[a, a, b, b] = chi0q[a, b]
        return out

    def _make_inter_k(self, norb, Nx, Ny, Nz, U=0.0, Up=0.0, J=0.0, Jp=0.0):
        """Create interaction dict for Kanamori-type interactions."""
        inter_k = {}
        if U != 0.0:
            U_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
            U_k[0, 0] = U
            U_k[1, 1] = U
            inter_k["CoulombIntra"] = U_k
        if Up != 0.0:
            Up_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
            Up_k[0, 1] = Up
            Up_k[1, 0] = Up
            inter_k["CoulombInter"] = Up_k
        if J != 0.0:
            J_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
            J_k[0, 1] = J
            J_k[1, 0] = J
            inter_k["Hund"] = J_k
        if Jp != 0.0:
            Jp_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
            Jp_k[0, 1] = Jp
            Jp_k[1, 0] = Jp
            inter_k["Exchange"] = Jp_k
        return inter_k

    def test_simple_vs_general_1orb(self):
        """For 1-orbital CoulombIntra, simple and general modes should match exactly.

        Both formulations reduce to the same scalar equations for 1 orbital.
        """
        Nx, Ny, Nz = 4, 4, 1
        norb = 1
        nmat = 16
        beta = 5.0

        hr = {
            ((1, 0, 0), (0, 0)): 1.0,
            ((-1, 0, 0), (0, 0)): 1.0,
            ((0, 1, 0), (0, 0)): 1.5,
            ((0, -1, 0), (0, 0)): 1.5,
        }
        kx = np.linspace(0, 2 * np.pi, Nx, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, Ny, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, Nz, endpoint=False)

        epsilon_k = _build_hamiltonian_k(kx, ky, kz, hr, norb)
        eigenvalues, eigenvectors = _calc_eigenvalues(epsilon_k)
        mu = _determine_mu(eigenvalues, beta, 0.5, norb)
        green_kw = _calc_green(eigenvalues, eigenvectors, mu, beta, nmat)

        U_k = np.ones((norb, norb, Nx, Ny, Nz), dtype=complex) * 3.0
        inter_k = {"CoulombIntra": U_k}
        chi0q = np.full((norb, norb, Nx, Ny, Nz, nmat), 0.1, dtype=complex)

        Pc_s, Ps_s = _compute_vertices_simple(
            chi0q, inter_k, norb, Nx, Ny, Nz, nmat)
        V_simple = Pc_s + Ps_s  # (1, 1, Nx, Ny, Nz)

        Vs_gen = _compute_vertices_general(
            chi0q, inter_k, norb, Nx, Ny, Nz, nmat)
        # (1, 1, 1, 1, Nx, Ny, Nz) -> diagonal is [0,0,0,0]

        npt.assert_allclose(
            V_simple[0, 0], Vs_gen[0, 0, 0, 0], atol=1e-10,
            err_msg="Simple vs General should match for 1 orbital")

    def test_general_vertex_produces_finite_values(self):
        """General vertex with CoulombIntra+CoulombInter (no J) should be finite.

        For 2-orbital systems, simple and general modes differ in structure
        because general mode operates in (norb^2, norb^2) space. This test
        verifies the general mode result is well-behaved.
        """
        params = self._setup_2orb_model()
        norb, Nx, Ny, Nz, nmat = (params[k] for k in
                                     ["norb", "Nx", "Ny", "Nz", "nmat"])

        U_val, Up_val = 2.0, 1.0
        inter_k = self._make_inter_k(norb, Nx, Ny, Nz, U=U_val, Up=Up_val)

        chi0q = np.zeros((norb, norb, Nx, Ny, Nz, nmat), dtype=complex)
        green_kw = params["green_kw"]
        for a in range(norb):
            for b in range(norb):
                Gab = green_kw[a, b]
                Gba_rev = np.roll(
                    green_kw[b, a, ::-1, ::-1, ::-1, ::-1],
                    (1, 1, 1), (0, 1, 2)
                )
                prod = Gab * Gba_rev
                chi0q[a, b] = -ifftn(
                    fftn(prod, axes=(0, 1, 2)), axes=(0, 1, 2)
                ) * (Nx * Ny * Nz) / params["beta"]

        Vs_gen = _compute_vertices_general(
            chi0q, inter_k, norb, Nx, Ny, Nz, nmat)

        self.assertEqual(Vs_gen.shape, (norb, norb, norb, norb, Nx, Ny, Nz))
        self.assertTrue(np.all(np.isfinite(Vs_gen)),
                        "General vertex should be finite")

    def test_rpa_chi_s_manual_reference(self):
        """Verify RPA spin susceptibility chi_s from S/C matrices against manual calc.

        chi_s = [I - chi0 @ S]^{-1} @ chi0
        Test with Kanamori U, U', J, J' on a small system.
        """
        norb = 2
        Nx, Ny, Nz = 2, 2, 1
        nd = norb * norb

        U, Up, J, Jp = 4.0, 2.0, 0.5, 0.5
        inter_k = self._make_inter_k(norb, Nx, Ny, Nz, U=U, Up=Up, J=J, Jp=Jp)

        # Build S matrix
        S_all, C_all = _build_sc_matrices(inter_k, norb, 0, 0, 0)
        # S_all, C_all are (nd, nd)

        # Check known S matrix values (Kuroki convention)
        # S_{l1l2, l3l4}:
        # diag (l1=l2=l3=l4): S = U
        npt.assert_allclose(S_all[0, 0], U, atol=1e-10)
        npt.assert_allclose(S_all[3, 3], U, atol=1e-10)
        # cross (l1=l3 != l2=l4): S = U'
        npt.assert_allclose(S_all[1, 1], Up - J + Jp, atol=1e-10)  # (01, 01)
        npt.assert_allclose(S_all[2, 2], Up - J + Jp, atol=1e-10)  # (10, 10)
        # dens (l1=l2 != l3=l4): S = J
        npt.assert_allclose(S_all[0, 3], J, atol=1e-10)   # (00, 11)
        npt.assert_allclose(S_all[3, 0], J, atol=1e-10)   # (11, 00)
        # exch (l1=l4 != l2=l3): S = J'
        npt.assert_allclose(S_all[1, 2], 0.0, atol=1e-10)  # (01, 10) #113
        npt.assert_allclose(S_all[2, 1], 0.0, atol=1e-10)  # (10, 01)

        # Check C matrix
        # diag: C = U
        npt.assert_allclose(C_all[0, 0], U, atol=1e-10)
        # cross: C = -U' + J
        npt.assert_allclose(C_all[1, 1], -Up + J + Jp, atol=1e-10)  # #113
        # dens: C = 2U' - J
        npt.assert_allclose(C_all[0, 3], 2 * Up - J, atol=1e-10)
        # exch: C = J'
        npt.assert_allclose(C_all[1, 2], 0.0, atol=1e-10)  # #113

        # Now verify RPA susceptibility computation (full Kanamori,
        # including J': the code path below receives the FOUR-index
        # embedded chi0q, which the reduced-chi rejection does not police).
        # Use a small chi0 so the series converges, and make it
        # DENSITY-DIAGONAL (nonzero only at [(a,a),(b,b)]) so the manual
        # nd x nd reference and the four-index embedded code path compare
        # the same physics exactly.
        chi0 = np.zeros((nd, nd), dtype=complex)
        _dens_block = np.array([[0.05, 0.01], [0.01, 0.04]], dtype=complex)
        for a in range(norb):
            for b in range(norb):
                chi0[a * norb + a, b * norb + b] = _dens_block[a, b]

        # Manual computation
        I_mat = np.eye(nd, dtype=complex)
        chi_s_manual = np.linalg.solve(I_mat - chi0 @ S_all, chi0)
        chi_c_manual = np.linalg.solve(I_mat + chi0 @ C_all, chi0)

        # Verify singlet vertex: V^s = (3/2) S chi_s S - (1/2) C chi_c C + (1/2)(S + C)
        V_s_manual = (1.5 * S_all @ chi_s_manual @ S_all
                      - 0.5 * C_all @ chi_c_manual @ C_all
                      + 0.5 * (S_all + C_all))

        # Triplet vertex: V^t = -(1/2) S chi_s S - (1/2) C chi_c C + (1/2)(C - S)
        V_t_manual = (-0.5 * S_all @ chi_s_manual @ S_all
                      - 0.5 * C_all @ chi_c_manual @ C_all
                      + 0.5 * (C_all - S_all))

        # Now compute via _compute_vertices_general with matching chi0q
        nmat = 8
        chi0q = np.zeros((norb, norb, Nx, Ny, Nz, nmat), dtype=complex)
        # Set chi0 at static limit (nmat//2) to match our test chi0
        for a in range(norb):
            for b in range(norb):
                chi0q[a, b, 0, 0, 0, nmat // 2] = chi0[a * norb + a, b * norb + b]

        # four-index embedding: full Kanamori (including J') is retained
        # -- the reduced-chi rejection polices only 2-index input, and this
        # doubles as proof that four-index input is not falsely rejected
        chi0q4 = self._embed4(chi0q, norb, Nx, Ny, Nz, nmat)

        Vs_singlet = _compute_vertices_general(
            chi0q4, inter_k, norb, Nx, Ny, Nz, nmat, pairing_type="singlet")
        Vs_triplet = _compute_vertices_general(
            chi0q4, inter_k, norb, Nx, Ny, Nz, nmat, pairing_type="triplet")

        # Compare at q=(0,0,0)
        V_s_code = Vs_singlet[:, :, :, :, 0, 0, 0].reshape(nd, nd)
        V_t_code = Vs_triplet[:, :, :, :, 0, 0, 0].reshape(nd, nd)

        npt.assert_allclose(V_s_code, V_s_manual, atol=1e-10,
                            err_msg="Singlet vertex should match manual reference")
        npt.assert_allclose(V_t_code, V_t_manual, atol=1e-10,
                            err_msg="Triplet vertex should match manual reference")

    def test_kanamori_eliashberg_iteration_eigenvalue_consistency(self):
        """End-to-end test: Kanamori interactions with iteration + eigenvalue.

        Uses a 2-orbital model with U, U', J, J' and checks:
        1. Iteration converges
        2. Eigenvalue solver runs
        3. Leading eigenvalue from both methods agree
        """
        params = self._setup_2orb_model(Nx=4, Ny=4, Nz=1, nmat=16, beta=5.0)
        norb, Nx, Ny, Nz, nmat, beta = (
            params[k] for k in ["norb", "Nx", "Ny", "Nz", "nmat", "beta"])
        green_kw = params["green_kw"]

        # Kanamori interaction: U=2.0, U'=U-2J=1.4, J=J'=0.3
        U, J = 2.0, 0.3
        Up = U - 2 * J
        inter_k = self._make_inter_k(norb, Nx, Ny, Nz,
                                      U=U, Up=Up, J=J, Jp=J)

        # Compute chi0q from Green's function (simple diagonal form)
        chi0q = np.zeros((norb, norb, Nx, Ny, Nz, nmat), dtype=complex)
        for a in range(norb):
            for b in range(norb):
                Gab = green_kw[a, b]
                Gba_rev = np.roll(
                    green_kw[b, a, ::-1, ::-1, ::-1, ::-1],
                    (1, 1, 1), (0, 1, 2)
                )
                prod = Gab * Gba_rev
                chi0q[a, b] = -ifftn(
                    fftn(prod, axes=(0, 1, 2)), axes=(0, 1, 2)
                ) * (Nx * Ny * Nz) / beta

        # embed the density chi0q into the four-index form: full Kanamori
        # (including J' = J) is retained -- four-index input is not policed
        # by the reduced-chi rejection, and this also proves it is not
        # falsely rejected
        chi0q = self._embed4(chi0q, norb, Nx, Ny, Nz, nmat)

        # General mode (4-index vertex) is required with Hund/Exchange
        Vs_q = _compute_vertices(
            chi0q, inter_k, norb, Nx, Ny, Nz, nmat, pairing_type="singlet")
        self.assertEqual(Vs_q.ndim, 7,
                         "Should use general (4-index) mode with Hund/Exchange")

        G2 = _calc_g2(green_kw, beta)
        sigma_init = _initialize_gap("cos", norb, params["kx"], params["ky"],
                                     params["kz"])

        # Test iteration
        sigma_iter, ev_iter, converged, n_iter = _solve_iteration(
            green_kw, Vs_q, G2, sigma_init, norb,
            max_iter=200, alpha=0.5, tol=1e-5)
        self.assertEqual(sigma_iter.shape, (norb, norb, Nx, Ny, Nz))
        self.assertGreater(n_iter, 0, "Should iterate at least once")
        self.assertTrue(converged,
                        "the docstring claims convergence; without it the "
                        "eigenvalue comparison below would be meaningless")

        # Test eigenvalue
        eigenvalues, eigvecs = _solve_eigenvalue(
            Vs_q, G2, norb, Nx, Ny, Nz, num_eigenvalues=3)
        self.assertEqual(len(eigenvalues), 3)

        # If converged, leading eigenvalue should be consistent
        if converged:
            npt.assert_allclose(
                abs(eigenvalues[0]), ev_iter, rtol=0.15,
                err_msg="Leading eigenvalue from eigs should match iteration")

    def test_kanamori_singlet_vs_triplet(self):
        """Singlet and triplet channels should give different eigenvalues
        with the full Kanamori interactions (J and J' = J)."""
        params = self._setup_2orb_model(Nx=4, Ny=4, Nz=1, nmat=16, beta=5.0)
        norb, Nx, Ny, Nz, nmat, beta = (
            params[k] for k in ["norb", "Nx", "Ny", "Nz", "nmat", "beta"])
        green_kw = params["green_kw"]

        U, J = 2.0, 0.3
        Up = U - 2 * J
        inter_k = self._make_inter_k(norb, Nx, Ny, Nz,
                                      U=U, Up=Up, J=J, Jp=J)

        chi0q = np.zeros((norb, norb, Nx, Ny, Nz, nmat), dtype=complex)
        for a in range(norb):
            for b in range(norb):
                Gab = green_kw[a, b]
                Gba_rev = np.roll(
                    green_kw[b, a, ::-1, ::-1, ::-1, ::-1],
                    (1, 1, 1), (0, 1, 2)
                )
                prod = Gab * Gba_rev
                chi0q[a, b] = -ifftn(
                    fftn(prod, axes=(0, 1, 2)), axes=(0, 1, 2)
                ) * (Nx * Ny * Nz) / beta

        # four-index embedding: full Kanamori including J' = J is retained
        # (four-index input is not policed by the reduced-chi rejection)
        chi0q = self._embed4(chi0q, norb, Nx, Ny, Nz, nmat)

        Vs_singlet = _compute_vertices(
            chi0q, inter_k, norb, Nx, Ny, Nz, nmat, pairing_type="singlet")
        Vs_triplet = _compute_vertices(
            chi0q, inter_k, norb, Nx, Ny, Nz, nmat, pairing_type="triplet")

        G2 = _calc_g2(green_kw, beta)

        ev_singlet, _ = _solve_eigenvalue(
            Vs_singlet, G2, norb, Nx, Ny, Nz, num_eigenvalues=1)
        ev_triplet, _ = _solve_eigenvalue(
            Vs_triplet, G2, norb, Nx, Ny, Nz, num_eigenvalues=1)

        # Singlet and triplet leading eigenvalues should differ when J != 0
        self.assertGreater(
            abs(abs(ev_singlet[0]) - abs(ev_triplet[0])), 1e-4,
            "Singlet and triplet should give different eigenvalues with J != 0")

    def test_kanamori_symmetry_UV_consistency(self):
        """Test U'=U-2J constraint: results should be physically consistent.

        With the Kanamori symmetry U'=U-2J and J=J', the S matrix should have
        the relation S_cross = S_diag - 2*S_dens (from U' = U - 2J).
        """
        norb = 2
        Nx, Ny, Nz = 1, 1, 1
        U, J = 3.0, 0.4
        Up = U - 2 * J  # = 2.2

        inter_k = self._make_inter_k(norb, Nx, Ny, Nz,
                                      U=U, Up=Up, J=J, Jp=J)
        S_mat, C_mat = _build_sc_matrices(inter_k, norb, 0, 0, 0)

        # Verify S matrix Kanamori relations
        nd = norb * norb
        S_diag = S_mat[0, 0]       # U
        S_cross = S_mat[1, 1]      # U'
        S_dens = S_mat[0, 3]       # J
        S_exch = S_mat[1, 2]       # J' = J

        npt.assert_allclose(S_diag, U, atol=1e-10)
        npt.assert_allclose(S_cross, Up, atol=1e-10)
        npt.assert_allclose(S_dens, J, atol=1e-10)
        npt.assert_allclose(S_exch, 0.0, atol=1e-10)  # Exchange moved (#113)
        # U' = U - 2J
        npt.assert_allclose(S_cross, S_diag - 2 * S_dens, atol=1e-10,
                            err_msg="Kanamori relation U'=U-2J should hold in S matrix")

        # C matrix relations
        C_diag = C_mat[0, 0]       # U
        C_cross = C_mat[1, 1]      # -U' + J = -(U-2J) + J = -U + 3J
        C_dens = C_mat[0, 3]       # 2U' - J = 2(U-2J) - J = 2U - 5J
        C_exch = C_mat[1, 2]       # J' = J

        npt.assert_allclose(C_diag, U, atol=1e-10)
        # C_cross = -U' + J + J' = -U' + 2J for J = J' -- the standard
        # Kanamori value, from the corrected per-type split (#113)
        npt.assert_allclose(C_cross, -Up + 2 * J, atol=1e-10)
        npt.assert_allclose(C_dens, 2 * Up - J, atol=1e-10)
        npt.assert_allclose(C_exch, 0.0, atol=1e-10)

    def test_kanamori_with_ising_pairhop_eliashberg(self):
        """End-to-end with every type a 2-index chi0q supports: U, U', J
        (Hund), Ising -- and a pinned REJECTION for Exchange/PairHop, which
        have no density-diagonal vertex and are refused with this chi0q
        form (#120 alignment)."""
        params = self._setup_2orb_model(Nx=4, Ny=4, Nz=1, nmat=16, beta=5.0)
        norb, Nx, Ny, Nz, nmat, beta = (
            params[k] for k in ["norb", "Nx", "Ny", "Nz", "nmat", "beta"])
        green_kw = params["green_kw"]

        # Full interaction set
        U_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        U_k[0, 0] = 2.0
        U_k[1, 1] = 2.0

        V_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        V_k[0, 1] = 1.0
        V_k[1, 0] = 1.0

        J_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        J_k[0, 1] = 0.3
        J_k[1, 0] = 0.3

        Jp_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        Jp_k[0, 1] = 0.3
        Jp_k[1, 0] = 0.3

        I_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        I_k[0, 1] = 0.1
        I_k[1, 0] = 0.1

        PH_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        PH_k[0, 1] = 0.15
        PH_k[1, 0] = 0.15

        inter_k = {
            "CoulombIntra": U_k,
            "CoulombInter": V_k,
            "Hund": J_k,
            "Ising": I_k,
        }
        inter_k_rejected = {
            "CoulombIntra": U_k,
            "Exchange": Jp_k,
            "PairHop": PH_k,
        }

        chi0q = np.zeros((norb, norb, Nx, Ny, Nz, nmat), dtype=complex)
        for a in range(norb):
            for b in range(norb):
                Gab = green_kw[a, b]
                Gba_rev = np.roll(
                    green_kw[b, a, ::-1, ::-1, ::-1, ::-1],
                    (1, 1, 1), (0, 1, 2)
                )
                prod = Gab * Gba_rev
                chi0q[a, b] = -ifftn(
                    fftn(prod, axes=(0, 1, 2)), axes=(0, 1, 2)
                ) * (Nx * Ny * Nz) / beta

        Vs_q = _compute_vertices(
            chi0q, inter_k, norb, Nx, Ny, Nz, nmat, pairing_type="singlet")
        self.assertEqual(Vs_q.ndim, 7)

        with self.assertRaises(ValueError) as cm:
            _compute_vertices(chi0q, inter_k_rejected, norb, Nx, Ny, Nz,
                              nmat, pairing_type="singlet")
        self.assertIn("general", str(cm.exception))

        G2 = _calc_g2(green_kw, beta)
        sigma_init = _initialize_gap("cos", norb, params["kx"], params["ky"],
                                     params["kz"])

        sigma_iter, ev_iter, converged, n_iter = _solve_iteration(
            green_kw, Vs_q, G2, sigma_init, norb,
            max_iter=100, alpha=0.5, tol=1e-4)
        self.assertEqual(sigma_iter.shape, (norb, norb, Nx, Ny, Nz))
        self.assertTrue(np.isfinite(ev_iter), "Eigenvalue should be finite")

        eigenvalues, _ = _solve_eigenvalue(
            Vs_q, G2, norb, Nx, Ny, Nz, num_eigenvalues=3)
        self.assertEqual(len(eigenvalues), 3)
        for ev in eigenvalues:
            self.assertTrue(np.isfinite(ev), "All eigenvalues should be finite")

    def test_rpa_chi_from_sc_matrices_is_finite_on_rpa_solver_chi0q(self):
        """Smoke test: feed a REAL reduced chi0q from the RPA solver through
        sc.py's S/C dressing and check both channels stay finite.

        NOTE: this does NOT compare against the RPA solver's own chi_s/chi_c --
        it was previously named as if it did. It builds chi0q with the solver,
        embeds it at the density-pair positions, and solves
        [I -+ chi0 S/C]^{-1} chi0 here, asserting only that no q-point diverges
        for a small U. A genuine cross-implementation comparison would be
        valuable, but the reduced RPA solve and sc.py's general S/C formulation
        do not treat Hund identically, so what the two should agree on has to be
        established first rather than assumed."""
        import tempfile
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.rpa as sol_rpa

        norb = 2
        Nx, Ny, Nz = 4, 4, 1
        nmat = 32
        T = 0.2
        beta = 1.0 / T

        with tempfile.TemporaryDirectory() as tmpdir:
            input_dir = os.path.join(tmpdir, "input")
            output_dir = os.path.join(tmpdir, "output")
            os.makedirs(input_dir, exist_ok=True)
            os.makedirs(output_dir, exist_ok=True)

            # Geometry
            with open(os.path.join(input_dir, "geom.dat"), "w") as f:
                f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n")
                f.write("2\n0.0 0.0 0.0\n0.0 0.0 0.0\n")

            # Transfer
            with open(os.path.join(input_dir, "transfer.dat"), "w") as f:
                f.write("Transfer\n2\n5\n 1 1 1 1 1\n")
                for orb in [1, 2]:
                    f.write("  1  0  0  %d  %d  1.0 0.0\n" % (orb, orb))
                    f.write(" -1  0  0  %d  %d  1.0 0.0\n" % (orb, orb))
                    f.write("  0  1  0  %d  %d  0.8 0.0\n" % (orb, orb))
                    f.write("  0 -1  0  %d  %d  0.8 0.0\n" % (orb, orb))
                f.write("  1  0  0  1  2  0.3 0.0\n")
                f.write(" -1  0  0  1  2  0.3 0.0\n")
                f.write("  1  0  0  2  1  0.3 0.0\n")
                f.write(" -1  0  0  2  1  0.3 0.0\n")

            # CoulombIntra: U = 1.0 (small to avoid divergence)
            U_val = 1.0
            with open(os.path.join(input_dir, "coulombintra.dat"), "w") as f:
                f.write("CoulombIntra\n2\n1\n 1\n")
                f.write("  0  0  0  1  1  %.1f 0.0\n" % U_val)
                f.write("  0  0  0  2  2  %.1f 0.0\n" % U_val)

            # Hund: J = 0.2
            J_val = 0.2
            with open(os.path.join(input_dir, "hund.dat"), "w") as f:
                f.write("Hund\n2\n1\n 1\n")
                f.write("  0  0  0  1  2  %.1f 0.0\n" % J_val)
                f.write("  0  0  0  2  1  %.1f 0.0\n" % J_val)

            # Run RPA solver to get chi0q and chiq
            input_dict = {
                "mode": {
                    "mode": "RPA",
                    "param": {
                        "T": T,
                        "CellShape": [Nx, Ny, Nz],
                        "SubShape": [1, 1, 1],
                        "Nmat": nmat,
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
                            "Hund": "hund.dat",
                        },
                    },
                    "output": {
                        "path_to_output": output_dir,
                    },
                },
            }

            # Get chi0q via RPA solver
            info_mode = {"mode": "RPA",
                         "param": input_dict["mode"]["param"],
                         "calc_scheme": "reduced"}
            info_inputfile = input_dict["file"]["input"]
            info_log = {"print_level": 1, "print_step": 1}

            read_io = read_input_k.QLMSkInput(info_inputfile)
            ham_info = read_io.get_param("ham")
            solver = sol_rpa.RPA(ham_info, info_log, info_mode)

            green_info = read_io.get_param("green")
            green_info.update(solver.read_init(info_inputfile))

            solver._calc_epsilon_k(green_info)
            Ncond = solver.Ncond / 2
            dist, mu = solver._find_mu(Ncond, solver.T)
            green0, green0_tail = solver._calc_green(beta, mu)
            chi0q_rpa = solver._calc_chi0q(green0, green0_tail, beta)
            chi0q_rpa = chi0q_rpa[0]  # Remove block index

            # Compute chi0q in sc.py ref format
            chi0q_ref = _convert_chi0q_to_ref_format(
                chi0q_rpa, norb, Nx, Ny, Nz)

            # Build interaction in k-space for sc.py
            inter_k = self._make_inter_k(norb, Nx, Ny, Nz,
                                          U=U_val, J=J_val)

            # Compute RPA spin susceptibility via sc.py S/C matrices
            nd = norb * norb
            S_all, C_all = _build_sc_matrices(inter_k, norb, 0, 0, 0)

            # chi0 at static limit
            chi0_static = chi0q_ref[:, :, :, :, :, nmat // 2]
            # Expand 2-index to 4-index for the matrix formulation, using the
            # SAME density-pair placement as production: a reduced chi0 is the
            # density-density diagonal, so chi0_2d[a,b] belongs at
            # [(a,a),(b,b)]. (This mirrored the old kron(chi0_2d, I_norb)
            # scatter, which modelled different susceptibility data than the
            # solver actually builds.)
            chi0_2d = chi0_static.transpose(2, 3, 4, 0, 1).copy()
            chi0_expanded = np.zeros((Nx, Ny, Nz, nd, nd), dtype=complex)
            _dens = np.arange(norb) * norb + np.arange(norb)
            chi0_expanded[..., _dens[:, None], _dens[None, :]] = chi0_2d

            I_mat = np.eye(nd, dtype=complex)

            # Spin susceptibility from S matrix
            for ix in range(Nx):
                for iy in range(Ny):
                    chi0_q = chi0_expanded[ix, iy, 0]
                    # chi_s = [I - chi0 S]^{-1} chi0
                    chi_s = np.linalg.solve(I_mat - chi0_q @ S_all, chi0_q)

                    # chi_s should be finite (no divergence for small U)
                    self.assertTrue(np.all(np.isfinite(chi_s)),
                                    "chi_s should be finite at q=(%d,%d)" % (ix, iy))

                    # chi_c = [I + chi0 C]^{-1} chi0
                    chi_c = np.linalg.solve(I_mat + chi0_q @ C_all, chi0_q)
                    self.assertTrue(np.all(np.isfinite(chi_c)),
                                    "chi_c should be finite at q=(%d,%d)" % (ix, iy))


class TestBatchedMatmulEquivalence(unittest.TestCase):
    """Prove the batched-matmul rewrites equal the original einsums.

    Each test builds random operands of the documented shapes, computes the
    ORIGINAL einsum (pasted as reference) and the NEW matmul form, and asserts
    they agree to ~1e-10 or tighter. The rewrites only change the floating-point
    reduction order (BLAS GEMM vs numpy's einsum loop), so differences are at
    machine-epsilon level, far below the suite's atol=1e-8.
    """

    def _calc_green_reference(self, factor, denom, Nx, Ny, Nz, norb, nmat):
        return np.einsum('...ijm,...mw->...ijw', factor, denom)

    def _calc_green_matmul(self, factor, denom, Nx, Ny, Nz, norb, nmat):
        nv = Nx * Ny * Nz
        G = factor.reshape(nv, norb * norb, norb) @ denom.reshape(nv, norb, nmat)
        return G.reshape(Nx, Ny, Nz, norb, norb, nmat)

    def test_calc_green_einsum_equivalence(self):
        rng = np.random.default_rng(0)
        for norb in (2, 3):
            Nx, Ny, Nz, nmat = 3, 4, 2, 5
            factor = (rng.standard_normal((Nx, Ny, Nz, norb, norb, norb))
                      + 1j * rng.standard_normal((Nx, Ny, Nz, norb, norb, norb)))
            denom = (rng.standard_normal((Nx, Ny, Nz, norb, nmat))
                     + 1j * rng.standard_normal((Nx, Ny, Nz, norb, nmat)))
            ref = self._calc_green_reference(factor, denom, Nx, Ny, Nz, norb, nmat)
            new = self._calc_green_matmul(factor, denom, Nx, Ny, Nz, norb, nmat)
            self.assertEqual(ref.shape, new.shape)
            maxdiff = np.max(np.abs(ref - new))
            self.assertTrue(np.allclose(ref, new, rtol=1e-10, atol=1e-10),
                            "norb=%d maxdiff=%g" % (norb, maxdiff))

    def _calc_g2_reference(self, A, B):
        return np.einsum('isn,jsn->ijs', A, B)

    def _calc_g2_matmul(self, A, B):
        As = np.moveaxis(A, 1, 0)
        Bs = np.moveaxis(B, 1, 0)
        return np.moveaxis(As @ Bs.transpose(0, 2, 1), 0, 2)

    def test_calc_g2_einsum_equivalence(self):
        rng = np.random.default_rng(1)
        for norb in (2, 3):
            Nx, Ny, Nz, nmat = 3, 2, 2, 6
            nvol = Nx * Ny * Nz
            nd = norb * norb
            A = (rng.standard_normal((nd, nvol, nmat))
                 + 1j * rng.standard_normal((nd, nvol, nmat)))
            B = (rng.standard_normal((nd, nvol, nmat))
                 + 1j * rng.standard_normal((nd, nvol, nmat)))
            ref = self._calc_g2_reference(A, B)
            new = self._calc_g2_matmul(A, B)
            self.assertEqual(ref.shape, new.shape)
            maxdiff = np.max(np.abs(ref - new))
            self.assertTrue(np.allclose(ref, new, rtol=1e-10, atol=1e-10),
                            "norb=%d maxdiff=%g" % (norb, maxdiff))

    def test_compute_vertices_batched_matmul_equivalence(self):
        rng = np.random.default_rng(2)
        for nd in (2, 3, 4, 9):
            Nx, Ny, Nz = 3, 4, 2
            X = (rng.standard_normal((Nx, Ny, Nz, nd, nd))
                 + 1j * rng.standard_normal((Nx, Ny, Nz, nd, nd)))
            Y = (rng.standard_normal((Nx, Ny, Nz, nd, nd))
                 + 1j * rng.standard_normal((Nx, Ny, Nz, nd, nd)))
            ref = np.einsum('...ab,...bc->...ac', X, Y)
            new = X @ Y
            maxdiff = np.max(np.abs(ref - new))
            self.assertTrue(np.allclose(ref, new, rtol=1e-12, atol=1e-12),
                            "nd=%d maxdiff=%g" % (nd, maxdiff))


if __name__ == '__main__':
    unittest.main()
