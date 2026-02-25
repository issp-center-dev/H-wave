#!/usr/bin/env python3
"""Tests for block matrix optimization.

Verifies that _solve_rpa with block-diagonal detection produces
identical results to the original full-matrix solve.
"""

import unittest
import numpy as np


class TestRPABlockDiagonal(unittest.TestCase):
    """Test block-diagonal optimization in RPA._solve_rpa."""

    def _make_solver_stub(self, nvol):
        """Create a minimal RPA-like object with _solve_rpa and _find_block_diagonal."""
        import hwave.solver.rpa as rpa_module

        class LatticeStub:
            pass

        stub = object.__new__(rpa_module.RPA)
        stub.lattice = LatticeStub()
        stub.lattice.nvol = nvol
        return stub

    def _solve_rpa_full(self, chi0q, ham, nvol):
        """Reference full-matrix solve (original algorithm)."""
        nmat = chi0q.shape[0]
        chi_shape = chi0q.shape
        ndx = np.prod(chi_shape[2:2 + (len(chi_shape) - 2) // 2])

        chi0q_2d = chi0q.reshape(nmat, nvol, ndx, ndx)
        ham_2d = ham.reshape(nvol, ndx, ndx)

        mat = np.tile(np.eye(ndx, dtype=np.complex128), (nmat, nvol, 1, 1))
        mat += np.einsum('lkab,kbc->lkac', chi0q_2d, ham_2d)
        sol = np.linalg.solve(mat, chi0q_2d)
        return sol.reshape(chi_shape)

    def test_block_diagonal_2x2_reduced(self):
        """Test block-diagonal case: 2 spin blocks, reduced scheme (nd x nd)."""
        nmat, nvol, norb, ns = 4, 8, 2, 2
        nd = norb * ns

        # Build block-diagonal chi0q: spin-up and spin-down blocks independent
        chi0q = np.zeros((nmat, nvol, nd, nd), dtype=np.complex128)
        rng = np.random.RandomState(42)
        for s in range(ns):
            blk = rng.randn(nmat, nvol, norb, norb) + 1j * rng.randn(nmat, nvol, norb, norb)
            # Make Hermitian per (k) for physical realism
            blk = (blk + blk.conj().transpose(0, 1, 3, 2)) / 2
            chi0q[:, :, s * norb:(s + 1) * norb, s * norb:(s + 1) * norb] = blk

        # Build block-diagonal ham
        ham = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for s in range(ns):
            blk = rng.randn(nvol, norb, norb) + 1j * rng.randn(nvol, norb, norb)
            blk = (blk + blk.conj().transpose(0, 2, 1)) / 2
            ham[:, s * norb:(s + 1) * norb, s * norb:(s + 1) * norb] = blk

        solver = self._make_solver_stub(nvol)
        sol_block = solver._solve_rpa(chi0q, ham)
        sol_full = self._solve_rpa_full(chi0q, ham, nvol)

        self.assertTrue(np.allclose(sol_block, sol_full, atol=1e-12),
                        "Block-diagonal solve must match full solve")

    def test_block_diagonal_general_scheme(self):
        """Test block-diagonal case: general scheme (nd^2 x nd^2)."""
        nmat, nvol, norb, ns = 2, 4, 1, 2
        nd = norb * ns

        # General scheme: chi0q shape (nmat, nvol, nd, nd, nd, nd)
        chi0q = np.zeros((nmat, nvol, nd, nd, nd, nd), dtype=np.complex128)
        rng = np.random.RandomState(123)
        # Fill spin-diagonal blocks only (block structure in nd^2 space)
        for s1 in range(ns):
            for s2 in range(ns):
                if s1 == s2:
                    a = s1 * norb
                    b = (s1 + 1) * norb
                    blk = rng.randn(nmat, nvol, norb, norb, norb, norb) \
                        + 1j * rng.randn(nmat, nvol, norb, norb, norb, norb)
                    chi0q[:, :, a:b, a:b, a:b, a:b] = blk

        ham = np.zeros((nvol, nd, nd, nd, nd), dtype=np.complex128)
        for s1 in range(ns):
            a = s1 * norb
            b = (s1 + 1) * norb
            blk = rng.randn(nvol, norb, norb, norb, norb) \
                + 1j * rng.randn(nvol, norb, norb, norb, norb)
            ham[:, a:b, a:b, a:b, a:b] = blk

        solver = self._make_solver_stub(nvol)
        sol_block = solver._solve_rpa(chi0q, ham)
        sol_full = self._solve_rpa_full(chi0q, ham, nvol)

        self.assertTrue(np.allclose(sol_block, sol_full, atol=1e-12),
                        "Block-diagonal solve (general) must match full solve")

    def test_no_block_structure(self):
        """Test that full matrices without block structure still work correctly."""
        nmat, nvol, nd = 3, 4, 4
        rng = np.random.RandomState(99)

        chi0q = rng.randn(nmat, nvol, nd, nd) + 1j * rng.randn(nmat, nvol, nd, nd)
        chi0q = (chi0q + chi0q.conj().transpose(0, 1, 3, 2)) / 2

        ham = rng.randn(nvol, nd, nd) + 1j * rng.randn(nvol, nd, nd)
        ham = (ham + ham.conj().transpose(0, 2, 1)) / 2

        solver = self._make_solver_stub(nvol)
        sol_block = solver._solve_rpa(chi0q, ham)
        sol_full = self._solve_rpa_full(chi0q, ham, nvol)

        self.assertTrue(np.allclose(sol_block, sol_full, atol=1e-12),
                        "Full matrix (no blocks) must produce identical results")

    def test_find_block_diagonal_detection(self):
        """Test _find_block_diagonal correctly identifies block structure."""
        import hwave.solver.rpa as rpa_module

        class LatticeStub:
            nvol = 2

        solver = object.__new__(rpa_module.RPA)
        solver.lattice = LatticeStub()

        # 4x4 matrix with 2x2 block-diagonal structure
        ham = np.zeros((2, 4, 4), dtype=np.complex128)
        ham[:, 0:2, 0:2] = np.random.randn(2, 2, 2)
        ham[:, 2:4, 2:4] = np.random.randn(2, 2, 2)

        blocks = solver._find_block_diagonal(ham)
        self.assertIsNotNone(blocks)
        self.assertEqual(len(blocks), 2)
        block_sizes = sorted([len(b) for b in blocks])
        self.assertEqual(block_sizes, [2, 2])

    def test_find_block_diagonal_no_blocks(self):
        """Test _find_block_diagonal returns None for full matrix."""
        import hwave.solver.rpa as rpa_module

        class LatticeStub:
            nvol = 2

        solver = object.__new__(rpa_module.RPA)
        solver.lattice = LatticeStub()

        ham = np.random.randn(2, 4, 4) + 1j * np.random.randn(2, 4, 4)
        blocks = solver._find_block_diagonal(ham)
        self.assertIsNone(blocks)

    def test_three_blocks(self):
        """Test with 3 independent blocks of different sizes."""
        nmat, nvol = 2, 4
        # Block sizes: 2, 3, 1 -> ndx = 6
        ndx = 6
        rng = np.random.RandomState(77)

        chi0q = np.zeros((nmat, nvol, ndx, ndx), dtype=np.complex128)
        ham = np.zeros((nvol, ndx, ndx), dtype=np.complex128)

        block_ranges = [(0, 2), (2, 5), (5, 6)]
        for a, b in block_ranges:
            sz = b - a
            c = rng.randn(nmat, nvol, sz, sz) + 1j * rng.randn(nmat, nvol, sz, sz)
            chi0q[:, :, a:b, a:b] = (c + c.conj().transpose(0, 1, 3, 2)) / 2
            h = rng.randn(nvol, sz, sz) + 1j * rng.randn(nvol, sz, sz)
            ham[:, a:b, a:b] = (h + h.conj().transpose(0, 2, 1)) / 2

        solver = self._make_solver_stub(nvol)
        sol_block = solver._solve_rpa(chi0q, ham)
        sol_full = self._solve_rpa_full(chi0q, ham, nvol)

        self.assertTrue(np.allclose(sol_block, sol_full, atol=1e-12),
                        "3-block solve must match full solve")


class TestRPAEighBlockDecomposition(unittest.TestCase):
    """Test that RPA _calc_epsilon_k block decomposition gives correct eigenvalues."""

    def _make_rpa_stub(self, nvol):
        """Create a minimal RPA-like object for eigenvalue tests."""
        import hwave.solver.rpa as rpa_module

        class LatticeStub:
            pass

        stub = object.__new__(rpa_module.RPA)
        stub.lattice = LatticeStub()
        stub.lattice.nvol = nvol
        return stub

    def test_eigh_block_diagonal_matches_full(self):
        """Block-diagonal eigh should produce same eigenvalues as full eigh.

        H0 shape (nblock_spin, nvol, nd, nd) with block-diagonal structure.
        """
        nblock_spin, nvol, nd = 1, 8, 6
        rng = np.random.RandomState(200)

        # Build block-diagonal H0: blocks of size 3 and 3
        H0 = np.zeros((nblock_spin, nvol, nd, nd), dtype=np.complex128)
        for a, b in [(0, 3), (3, 6)]:
            sz = b - a
            blk = rng.randn(nblock_spin, nvol, sz, sz) + \
                1j * rng.randn(nblock_spin, nvol, sz, sz)
            blk = (blk + blk.conj().transpose(0, 1, 3, 2)) / 2
            H0[:, :, a:b, a:b] = blk

        # Full eigh
        w_full, v_full = np.linalg.eigh(H0)

        # Block eigh via RPA stub
        stub = self._make_rpa_stub(nvol)
        blocks = stub._find_block_diagonal(
            H0.reshape(nblock_spin * nvol, nd, nd))
        self.assertIsNotNone(blocks)
        self.assertEqual(len(blocks), 2)

        # Reconstruct block eigenvalues
        w_block = np.zeros((nblock_spin, nvol, nd), dtype=np.float64)
        col = 0
        for blk_idx in blocks:
            idx = np.array(blk_idx)
            ix = np.ix_(idx, idx)
            wb, _ = np.linalg.eigh(H0[:, :, ix[0], ix[1]])
            nb = len(idx)
            w_block[:, :, col:col + nb] = wb
            col += nb

        # Sort eigenvalues per k-point for comparison
        w_full_sorted = np.sort(w_full.reshape(nblock_spin, nvol, -1), axis=-1)
        w_block_sorted = np.sort(w_block.reshape(nblock_spin, nvol, -1), axis=-1)

        self.assertTrue(np.allclose(w_full_sorted, w_block_sorted, atol=1e-12),
                        "Block-decomposed eigh eigenvalues must match full eigh")

    def test_eigh_spin_diag_with_orbital_blocks(self):
        """Spin-diagonal H0 (nblock=2) with orbital blocks within each spin.

        nblock=2, nvol=4, norb=4.
        Each spin block has 2+2 orbital structure.
        """
        nblock_spin, nvol, norb = 2, 4, 4
        rng = np.random.RandomState(210)

        H0 = np.zeros((nblock_spin, nvol, norb, norb), dtype=np.complex128)
        for g in range(nblock_spin):
            for a, b in [(0, 2), (2, 4)]:
                sz = b - a
                blk = rng.randn(nvol, sz, sz) + 1j * rng.randn(nvol, sz, sz)
                blk = (blk + blk.conj().transpose(0, 2, 1)) / 2
                H0[g, :, a:b, a:b] = blk

        w_full, _ = np.linalg.eigh(H0)
        w_full_sorted = np.sort(w_full.reshape(nblock_spin, nvol, -1), axis=-1)

        stub = self._make_rpa_stub(nvol)
        blocks = stub._find_block_diagonal(
            H0.reshape(nblock_spin * nvol, norb, norb))

        self.assertIsNotNone(blocks)
        self.assertEqual(len(blocks), 2)

        w_block = np.zeros_like(w_full)
        col = 0
        for blk_idx in blocks:
            idx = np.array(blk_idx)
            ix = np.ix_(idx, idx)
            wb, _ = np.linalg.eigh(H0[:, :, ix[0], ix[1]])
            nb = len(idx)
            w_block[:, :, col:col + nb] = wb
            col += nb

        w_block_sorted = np.sort(w_block.reshape(nblock_spin, nvol, -1), axis=-1)
        self.assertTrue(np.allclose(w_full_sorted, w_block_sorted, atol=1e-12),
                        "Spin-diag + orbital block eigh must match full eigh")

    def test_eigh_no_blocks_unchanged(self):
        """Full matrix (no blocks) should produce identical results."""
        nblock_spin, nvol, nd = 1, 4, 4
        rng = np.random.RandomState(220)

        H0 = rng.randn(nblock_spin, nvol, nd, nd) + \
            1j * rng.randn(nblock_spin, nvol, nd, nd)
        H0 = (H0 + H0.conj().transpose(0, 1, 3, 2)) / 2

        stub = self._make_rpa_stub(nvol)
        blocks = stub._find_block_diagonal(
            H0.reshape(nblock_spin * nvol, nd, nd))
        self.assertIsNone(blocks,
                          "Full matrix should not detect block structure")


class TestRPASpinOrbitalEquivalence(unittest.TestCase):
    """Test that RPA _solve_rpa and _find_block_diagonal produce equivalent
    results under spin-orbital index reordering.

    Normal ordering: flat index = s * norb + a  (spin-block first)
    SO ordering:     flat index = 2 * a + s      (interleaved spin-orbital)

    Since _solve_rpa computes chiq = [1 + chi0q * W]^(-1) * chi0q,
    the result must be invariant under simultaneous permutation of all indices.
    """

    def _make_solver_stub(self, nvol):
        """Create a minimal RPA stub with _solve_rpa and _find_block_diagonal."""
        import hwave.solver.rpa as rpa_module

        class LatticeStub:
            pass

        stub = object.__new__(rpa_module.RPA)
        stub.lattice = LatticeStub()
        stub.lattice.nvol = nvol
        return stub

    def _build_perm(self, norb):
        """Build permutation: perm[i_so] = i_normal.

        Normal: s*norb+a, SO: 2*a+s.
        """
        nd = 2 * norb
        perm = np.zeros(nd, dtype=int)
        for a in range(norb):
            for s in range(2):
                perm[2 * a + s] = s * norb + a
        return perm

    def _normal_to_so_2d(self, mat, perm):
        """Reorder last 2 indices from normal to SO ordering."""
        return mat[..., perm[:, None], perm[None, :]]

    def _so_to_normal_2d(self, mat, perm):
        """Reorder last 2 indices from SO to normal ordering."""
        inv_perm = np.argsort(perm)
        return mat[..., inv_perm[:, None], inv_perm[None, :]]

    def _normal_to_so_4d(self, mat, perm):
        """Reorder last 4 indices from normal to SO ordering."""
        ix = np.ix_(perm, perm, perm, perm)
        return mat[..., ix[0], ix[1], ix[2], ix[3]]

    def _so_to_normal_4d(self, mat, perm):
        """Reorder last 4 indices from SO to normal ordering."""
        inv_perm = np.argsort(perm)
        ix = np.ix_(inv_perm, inv_perm, inv_perm, inv_perm)
        return mat[..., ix[0], ix[1], ix[2], ix[3]]

    # ----------------------------------------------------------------
    # Block detection tests
    # ----------------------------------------------------------------
    def test_find_block_diagonal_so_ordering(self):
        """Block detection finds equivalent blocks in SO ordering.

        norb=3, nd=6.
        Normal: spin blocks [0,1,2] and [3,4,5].
        SO: equivalent blocks [0,2,4] and [1,3,5] (interleaved).
        """
        norb, nvol = 3, 4
        nd = 2 * norb
        rng = np.random.RandomState(310)

        ham = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for s in range(2):
            blk = rng.randn(nvol, norb, norb) + 1j * rng.randn(nvol, norb, norb)
            ham[:, s*norb:(s+1)*norb, s*norb:(s+1)*norb] = blk

        solver = self._make_solver_stub(nvol)
        blocks_normal = solver._find_block_diagonal(ham)

        perm = self._build_perm(norb)
        ham_so = self._normal_to_so_2d(ham, perm)
        blocks_so = solver._find_block_diagonal(ham_so)

        self.assertIsNotNone(blocks_normal)
        self.assertIsNotNone(blocks_so)
        self.assertEqual(len(blocks_normal), len(blocks_so))

        # Verify SO blocks map back to same normal blocks
        # perm[i_so] = i_normal, so perm maps SO indices to normal indices
        blocks_so_as_normal = [sorted(perm[b].tolist()) for b in blocks_so]
        blocks_normal_sorted = [sorted(b) for b in blocks_normal]
        self.assertEqual(sorted(map(tuple, blocks_so_as_normal)),
                         sorted(map(tuple, blocks_normal_sorted)))

    def test_find_block_diagonal_so_mixed_blocks(self):
        """SO block detection with orbital sub-blocks.

        norb=3, nd=6. Orbitals 0,1 coupled; orbital 2 isolated.
        Normal: blocks {0,1}, {2}, {3,4}, {5} (4 blocks).
        SO: blocks {0,2}, {4}, {1,3}, {5} (same structure, non-contiguous).
        """
        norb, nvol = 3, 4
        nd = 2 * norb
        rng = np.random.RandomState(320)

        ham = np.zeros((nvol, nd, nd), dtype=np.complex128)
        # Spin-up: orbitals 0-1 coupled
        ham[:, 0, 0] = rng.randn(nvol)
        ham[:, 1, 1] = rng.randn(nvol)
        ham[:, 0, 1] = rng.randn(nvol)
        ham[:, 1, 0] = rng.randn(nvol)
        ham[:, 2, 2] = rng.randn(nvol)  # orbital 2 isolated
        # Spin-down: same structure
        ham[:, 3, 3] = rng.randn(nvol)
        ham[:, 4, 4] = rng.randn(nvol)
        ham[:, 3, 4] = rng.randn(nvol)
        ham[:, 4, 3] = rng.randn(nvol)
        ham[:, 5, 5] = rng.randn(nvol)

        solver = self._make_solver_stub(nvol)
        blocks_normal = solver._find_block_diagonal(ham)

        perm = self._build_perm(norb)
        ham_so = self._normal_to_so_2d(ham, perm)
        blocks_so = solver._find_block_diagonal(ham_so)

        self.assertIsNotNone(blocks_normal)
        self.assertIsNotNone(blocks_so)
        self.assertEqual(len(blocks_normal), 4)
        self.assertEqual(len(blocks_so), 4)

        blocks_so_as_normal = [sorted(perm[b].tolist()) for b in blocks_so]
        blocks_normal_sorted = [sorted(b) for b in blocks_normal]
        self.assertEqual(sorted(map(tuple, blocks_so_as_normal)),
                         sorted(map(tuple, blocks_normal_sorted)))

    # ----------------------------------------------------------------
    # _solve_rpa equivalence tests (reduced scheme: nd x nd)
    # ----------------------------------------------------------------
    def test_solve_rpa_so_reduced_spin_diag(self):
        """Spin-diagonal chi0q/ham: normal vs SO in reduced scheme.

        Block structure: normal has contiguous blocks [0,1],[2,3];
        SO has interleaved blocks [0,2],[1,3].
        """
        norb, nmat, nvol = 2, 4, 8
        nd = 2 * norb
        rng = np.random.RandomState(300)

        chi0q = np.zeros((nmat, nvol, nd, nd), dtype=np.complex128)
        for s in range(2):
            blk = rng.randn(nmat, nvol, norb, norb) + \
                  1j * rng.randn(nmat, nvol, norb, norb)
            blk = (blk + blk.conj().transpose(0, 1, 3, 2)) / 2
            chi0q[:, :, s*norb:(s+1)*norb, s*norb:(s+1)*norb] = blk

        ham = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for s in range(2):
            blk = rng.randn(nvol, norb, norb) + 1j * rng.randn(nvol, norb, norb)
            blk = (blk + blk.conj().transpose(0, 2, 1)) / 2
            ham[:, s*norb:(s+1)*norb, s*norb:(s+1)*norb] = blk

        solver = self._make_solver_stub(nvol)
        sol_normal = solver._solve_rpa(chi0q, ham)

        perm = self._build_perm(norb)
        chi0q_so = self._normal_to_so_2d(chi0q, perm)
        ham_so = self._normal_to_so_2d(ham, perm)
        sol_so = solver._solve_rpa(chi0q_so, ham_so)
        sol_so_as_normal = self._so_to_normal_2d(sol_so, perm)

        np.testing.assert_allclose(
            sol_so_as_normal, sol_normal, atol=1e-12,
            err_msg="RPA solve (reduced, spin-diag) must be invariant under SO reordering"
        )

    def test_solve_rpa_so_reduced_full_matrix(self):
        """Full (non-block-diagonal) chi0q/ham: normal vs SO in reduced scheme."""
        norb, nmat, nvol = 2, 3, 4
        nd = 2 * norb
        rng = np.random.RandomState(301)

        chi0q = rng.randn(nmat, nvol, nd, nd) + 1j * rng.randn(nmat, nvol, nd, nd)
        chi0q = (chi0q + chi0q.conj().transpose(0, 1, 3, 2)) / 2

        ham = rng.randn(nvol, nd, nd) + 1j * rng.randn(nvol, nd, nd)
        ham = (ham + ham.conj().transpose(0, 2, 1)) / 2

        solver = self._make_solver_stub(nvol)
        sol_normal = solver._solve_rpa(chi0q, ham)

        perm = self._build_perm(norb)
        chi0q_so = self._normal_to_so_2d(chi0q, perm)
        ham_so = self._normal_to_so_2d(ham, perm)
        sol_so = solver._solve_rpa(chi0q_so, ham_so)
        sol_so_as_normal = self._so_to_normal_2d(sol_so, perm)

        np.testing.assert_allclose(
            sol_so_as_normal, sol_normal, atol=1e-12,
            err_msg="RPA solve (reduced, full matrix) must be invariant under SO reordering"
        )

    def test_solve_rpa_so_reduced_larger_system(self):
        """Larger system: norb=4, nd=8, spin-diagonal, reduced scheme."""
        norb, nmat, nvol = 4, 4, 4
        nd = 2 * norb
        rng = np.random.RandomState(302)

        chi0q = np.zeros((nmat, nvol, nd, nd), dtype=np.complex128)
        for s in range(2):
            blk = rng.randn(nmat, nvol, norb, norb) + \
                  1j * rng.randn(nmat, nvol, norb, norb)
            blk = (blk + blk.conj().transpose(0, 1, 3, 2)) / 2
            chi0q[:, :, s*norb:(s+1)*norb, s*norb:(s+1)*norb] = blk

        ham = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for s in range(2):
            blk = rng.randn(nvol, norb, norb) + 1j * rng.randn(nvol, norb, norb)
            blk = (blk + blk.conj().transpose(0, 2, 1)) / 2
            ham[:, s*norb:(s+1)*norb, s*norb:(s+1)*norb] = blk

        solver = self._make_solver_stub(nvol)
        sol_normal = solver._solve_rpa(chi0q, ham)

        perm = self._build_perm(norb)
        sol_so = solver._solve_rpa(
            self._normal_to_so_2d(chi0q, perm),
            self._normal_to_so_2d(ham, perm))
        sol_so_as_normal = self._so_to_normal_2d(sol_so, perm)

        np.testing.assert_allclose(
            sol_so_as_normal, sol_normal, atol=1e-12,
            err_msg="RPA solve (reduced, norb=4) must be invariant under SO reordering"
        )

    # ----------------------------------------------------------------
    # _solve_rpa equivalence tests (general scheme: nd^2 x nd^2)
    # ----------------------------------------------------------------
    def test_solve_rpa_so_general_spin_diag(self):
        """Spin-diagonal chi0q/ham: normal vs SO in general scheme.

        chi0q shape: (nmat, nvol, nd, nd, nd, nd)
        ham shape: (nvol, nd, nd, nd, nd)
        """
        norb, nmat, nvol = 2, 2, 4
        nd = 2 * norb
        rng = np.random.RandomState(400)

        chi0q = np.zeros((nmat, nvol, nd, nd, nd, nd), dtype=np.complex128)
        for s in range(2):
            a = s * norb
            b = (s + 1) * norb
            blk = rng.randn(nmat, nvol, norb, norb, norb, norb) + \
                  1j * rng.randn(nmat, nvol, norb, norb, norb, norb)
            chi0q[:, :, a:b, a:b, a:b, a:b] = blk

        ham = np.zeros((nvol, nd, nd, nd, nd), dtype=np.complex128)
        for s in range(2):
            a = s * norb
            b = (s + 1) * norb
            blk = rng.randn(nvol, norb, norb, norb, norb) + \
                  1j * rng.randn(nvol, norb, norb, norb, norb)
            ham[:, a:b, a:b, a:b, a:b] = blk

        solver = self._make_solver_stub(nvol)
        sol_normal = solver._solve_rpa(chi0q, ham)

        perm = self._build_perm(norb)
        chi0q_so = self._normal_to_so_4d(chi0q, perm)
        ham_so = self._normal_to_so_4d(ham, perm)
        sol_so = solver._solve_rpa(chi0q_so, ham_so)
        sol_so_as_normal = self._so_to_normal_4d(sol_so, perm)

        np.testing.assert_allclose(
            sol_so_as_normal, sol_normal, atol=1e-12,
            err_msg="RPA solve (general, spin-diag) must be invariant under SO reordering"
        )

    def test_solve_rpa_so_general_full_matrix(self):
        """Full chi0q/ham: normal vs SO in general scheme."""
        norb, nmat, nvol = 2, 2, 4
        nd = 2 * norb
        rng = np.random.RandomState(401)

        chi0q = rng.randn(nmat, nvol, nd, nd, nd, nd) + \
                1j * rng.randn(nmat, nvol, nd, nd, nd, nd)
        ham = rng.randn(nvol, nd, nd, nd, nd) + \
              1j * rng.randn(nvol, nd, nd, nd, nd)

        solver = self._make_solver_stub(nvol)
        sol_normal = solver._solve_rpa(chi0q, ham)

        perm = self._build_perm(norb)
        chi0q_so = self._normal_to_so_4d(chi0q, perm)
        ham_so = self._normal_to_so_4d(ham, perm)
        sol_so = solver._solve_rpa(chi0q_so, ham_so)
        sol_so_as_normal = self._so_to_normal_4d(sol_so, perm)

        np.testing.assert_allclose(
            sol_so_as_normal, sol_normal, atol=1e-12,
            err_msg="RPA solve (general, full matrix) must be invariant under SO reordering"
        )


class TestUHFkDetectBlocks(unittest.TestCase):
    """Test _detect_blocks in UHFk for various transfer/interaction patterns.

    Index convention: nd = norb * ns.
    For ns=2:  indices 0..norb-1 are spin-up, norb..2*norb-1 are spin-down.
    The (s, a) pair maps to index s*norb + a.
    """

    def _make_uhfk_stub(self, norb, ns, nvol, ham_trans,
                        inter_table=None, spin_table=None,
                        iflag_fock=True, sz_free=True, Nconds=None,
                        enable_spin_orbital=False):
        """Create a minimal UHFk-like object for _detect_blocks testing.

        Parameters
        ----------
        norb : int
            Number of orbitals.
        ns : int
            Number of spin species (2 for normal, 1 for spin-orbital).
        nvol : int
            Number of k-points (only shape matters for detection).
        ham_trans : ndarray (nvol, nd, nd)
            Transfer Hamiltonian in k-space.
        inter_table : dict or None
            Interaction tables {name: ndarray(nvol, norb, norb) or None}.
        spin_table : dict or None
            Spin combination tables {name: ndarray(ns, ns, ns, ns)}.
        iflag_fock : bool
            Whether Fock term is enabled.
        sz_free : bool
            Whether 2Sz is free (True) or fixed (False).
        Nconds : list or None
            Electron counts per block. Defaults to [norb*ns] if sz_free,
            or [norb//2, norb//2] if not sz_free.
        enable_spin_orbital : bool
            Whether spin-orbital mode is enabled.
        """
        import hwave.solver.uhfk as uhfk_module

        nd = norb * ns
        stub = object.__new__(uhfk_module.UHFk)
        stub.norb = norb
        stub.ns = ns
        stub.nd = nd
        stub.nvol = nvol
        stub.ham_trans = ham_trans
        stub.inter_table = inter_table if inter_table is not None else {}
        stub.spin_table = spin_table if spin_table is not None else {}
        stub.iflag_fock = iflag_fock
        stub.sz_free = sz_free
        stub.enable_spin_orbital = enable_spin_orbital
        if enable_spin_orbital:
            stub.norb_phys = norb // 2
        else:
            stub.norb_phys = norb
        if Nconds is not None:
            stub.Nconds = Nconds
        elif sz_free:
            stub.Nconds = [nd]
        else:
            stub.Nconds = [norb // 2, norb - norb // 2]
        return stub

    # ----------------------------------------------------------------
    # Pattern 1: Spin does NOT mix, orbitals DO mix
    #   Transfer: T[(up,a),(up,b)] != 0 and T[(dn,a),(dn,b)] != 0
    #             T[(up,a),(dn,b)] = 0
    #   Expected: 2 blocks of size norb (spin-up, spin-down)
    # ----------------------------------------------------------------
    def test_spin_separate_orbitals_mixed(self):
        """Transfer couples orbitals within each spin sector but not across spins.

        norb=3, ns=2 -> nd=6.  Indices: up=[0,1,2], dn=[3,4,5].
        Expected: 2 blocks [0,1,2] and [3,4,5].
        """
        norb, ns, nvol = 3, 2, 4
        nd = norb * ns
        rng = np.random.RandomState(10)

        ham_trans = np.zeros((nvol, nd, nd), dtype=np.complex128)
        # Fill spin-up block (0:3, 0:3)
        blk = rng.randn(nvol, norb, norb) + 1j * rng.randn(nvol, norb, norb)
        ham_trans[:, :norb, :norb] = blk
        # Fill spin-down block (3:6, 3:6)
        blk = rng.randn(nvol, norb, norb) + 1j * rng.randn(nvol, norb, norb)
        ham_trans[:, norb:, norb:] = blk

        stub = self._make_uhfk_stub(norb, ns, nvol, ham_trans, sz_free=False,
                                    Nconds=[2, 2])
        stub._detect_blocks()

        blocks = stub.block_info
        self.assertEqual(len(blocks), 2, "Should detect 2 spin blocks")
        sizes = sorted([len(b) for b in blocks])
        self.assertEqual(sizes, [3, 3])
        # Verify blocks correspond to pure spin sectors
        for blk in blocks:
            all_up = all(i < norb for i in blk)
            all_dn = all(i >= norb for i in blk)
            self.assertTrue(all_up or all_dn,
                            "Each block should be a pure spin sector")

    # ----------------------------------------------------------------
    # Pattern 2: Spin DOES mix (spin-orbit coupling via transfer)
    #   Transfer: T[(up,a),(dn,b)] != 0
    #   Expected: 1 single block of size nd
    # ----------------------------------------------------------------
    def test_spin_mixed_via_transfer(self):
        """Transfer has off-diagonal spin elements (spin-orbit coupling).

        norb=2, ns=2 -> nd=4.
        All indices connected via transfer -> single block [0,1,2,3].
        """
        norb, ns, nvol = 2, 2, 4
        nd = norb * ns
        rng = np.random.RandomState(20)

        # Full nd x nd transfer (all spin-orbital pairs connected)
        ham_trans = rng.randn(nvol, nd, nd) + 1j * rng.randn(nvol, nd, nd)

        stub = self._make_uhfk_stub(norb, ns, nvol, ham_trans, sz_free=True,
                                    Nconds=[nd])
        stub._detect_blocks()

        blocks = stub.block_info
        self.assertEqual(len(blocks), 1, "Spin mixing -> single block")
        self.assertEqual(len(blocks[0]), nd)

    # ----------------------------------------------------------------
    # Pattern 3: Orbitals do NOT mix, spins do NOT mix
    #   Transfer: T[(s,a),(s,a)] only (diagonal)
    #   No interactions.
    #   Expected: nd separate blocks of size 1 each.
    # ----------------------------------------------------------------
    def test_orbitals_separate_spins_separate(self):
        """Diagonal transfer only (no orbital or spin mixing).

        norb=3, ns=2 -> nd=6.
        Expected: 6 blocks of size 1.
        """
        norb, ns, nvol = 3, 2, 4
        nd = norb * ns
        rng = np.random.RandomState(30)

        ham_trans = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for i in range(nd):
            ham_trans[:, i, i] = rng.randn(nvol) + 1j * rng.randn(nvol)

        stub = self._make_uhfk_stub(norb, ns, nvol, ham_trans, sz_free=True,
                                    Nconds=[nd])
        stub._detect_blocks()

        blocks = stub.block_info
        self.assertEqual(len(blocks), nd,
                         "Fully diagonal transfer -> {} separate blocks".format(nd))
        for blk in blocks:
            self.assertEqual(len(blk), 1)

    # ----------------------------------------------------------------
    # Pattern 4: Orbitals DO mix within spin, plus interaction that
    #   does NOT mix spins (e.g. Hund coupling)
    #   Hund spin_table: spin[0,0,0,0]=1, spin[1,1,1,1]=1
    #   -> Fock term couples (s,a)-(s,b) for same spin only
    #   Expected: 2 blocks (spin-up, spin-down)
    # ----------------------------------------------------------------
    def test_hund_interaction_no_spin_mixing(self):
        """Hund interaction couples orbitals within same spin sector.

        norb=2, ns=2, nd=4.
        Transfer: diagonal only per orbital.
        Hund: J[a,b] != 0 for a != b, spin_table only same-spin.
        Expected: 2 blocks [0,1] and [2,3].
        """
        norb, ns, nvol = 2, 2, 4
        nd = norb * ns

        # Diagonal transfer (no orbital mixing from transfer alone)
        ham_trans = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for i in range(nd):
            ham_trans[:, i, i] = 1.0

        # Hund interaction: J[a,b] for a != b
        jab = np.zeros((nvol, norb, norb), dtype=np.complex128)
        jab[:, 0, 1] = 1.0
        jab[:, 1, 0] = 1.0

        spin = np.zeros((2, 2, 2, 2), dtype=int)
        spin[0, 0, 0, 0] = 1
        spin[1, 1, 1, 1] = 1

        inter_table = {"Hund": jab}
        spin_table = {"Hund": spin}

        stub = self._make_uhfk_stub(norb, ns, nvol, ham_trans,
                                    inter_table=inter_table,
                                    spin_table=spin_table,
                                    iflag_fock=True, sz_free=False,
                                    Nconds=[1, 1])
        stub._detect_blocks()

        blocks = stub.block_info
        self.assertEqual(len(blocks), 2,
                         "Hund (same-spin only) should give 2 spin blocks")
        sizes = sorted([len(b) for b in blocks])
        self.assertEqual(sizes, [2, 2])

    # ----------------------------------------------------------------
    # Pattern 5: Exchange interaction mixes spins
    #   Exchange spin_table: spin[0,1,0,1]=1, spin[1,0,1,0]=1
    #   Fock term: couples (s,a)-(t,b) with spin[s,u,t,v]
    #   -> spin_fock = max over (u,v) of |spin[s,u,t,v]|
    #   For Exchange: spin_fock[0,0]=|spin[0,1,0,1]|=1 via (u=1,v=1)
    #     Actually spin[s,u,t,v]: spin_fock[s,t] = max_{u,v} |spin[s,u,t,v]|
    #     spin[0,1,0,1]=1 -> spin_fock: max over (1,3) axes of |spin|
    #   -> Exchange Fock mixes spins -> single block
    # ----------------------------------------------------------------
    def test_exchange_interaction_mixes_spins(self):
        """Exchange interaction couples different spins via Fock term.

        norb=2, ns=2, nd=4.
        Transfer: spin-diagonal (no mixing).
        Exchange: J[a,b] with spin[0,1,0,1]=1, spin[1,0,1,0]=1.
        With Fock enabled, this couples (up,a)-(down,b) -> single block.
        """
        norb, ns, nvol = 2, 2, 4
        nd = norb * ns

        # Spin-diagonal transfer
        ham_trans = np.zeros((nvol, nd, nd), dtype=np.complex128)
        ham_trans[:, :norb, :norb] = np.eye(norb)
        ham_trans[:, norb:, norb:] = np.eye(norb)

        # Exchange: J[0,1] and J[1,0] nonzero
        jab = np.zeros((nvol, norb, norb), dtype=np.complex128)
        jab[:, 0, 1] = 0.5
        jab[:, 1, 0] = 0.5

        spin = np.zeros((2, 2, 2, 2), dtype=int)
        spin[0, 1, 0, 1] = 1
        spin[1, 0, 1, 0] = 1

        inter_table = {"Exchange": jab}
        spin_table = {"Exchange": spin}

        stub = self._make_uhfk_stub(norb, ns, nvol, ham_trans,
                                    inter_table=inter_table,
                                    spin_table=spin_table,
                                    iflag_fock=True, sz_free=True,
                                    Nconds=[nd])
        stub._detect_blocks()

        blocks = stub.block_info
        self.assertEqual(len(blocks), 1,
                         "Exchange Fock term mixes spins -> single block")
        self.assertEqual(len(blocks[0]), nd)

    # ----------------------------------------------------------------
    # Pattern 6: Exchange without Fock -> only Hartree term
    #   Hartree couples (s,a)-(t,a) (same orbital, cross spin)
    #   But Exchange spin_table only has spin[0,1,0,1], spin[1,0,1,0]
    #   Hartree: spin_mat_h[s,t] = any(spin[s,:,:,t])
    #     spin_mat_h[0,1] = any(spin[0,:,:,1]) = spin[0,1,0,1]=1  -> True
    #   So Hartree still couples (up,a)-(down,a) -> mixes spins
    # ----------------------------------------------------------------
    def test_exchange_without_fock_hartree_pairs(self):
        """Exchange Hartree term couples spins per orbital (no orbital mixing).

        norb=2, ns=2, nd=4. Fock disabled. Transfer diagonal.
        Exchange Hartree couples (up,a)-(down,a) for each a.
        Since transfer is diagonal, orbitals stay separate.
        Expected: 2 blocks {0,2} and {1,3} (orbital-spin pairs).
        """
        norb, ns, nvol = 2, 2, 4
        nd = norb * ns

        ham_trans = np.zeros((nvol, nd, nd), dtype=np.complex128)
        ham_trans[:, :norb, :norb] = np.eye(norb)
        ham_trans[:, norb:, norb:] = np.eye(norb)

        jab = np.zeros((nvol, norb, norb), dtype=np.complex128)
        jab[:, 0, 1] = 0.5
        jab[:, 1, 0] = 0.5

        spin = np.zeros((2, 2, 2, 2), dtype=int)
        spin[0, 1, 0, 1] = 1
        spin[1, 0, 1, 0] = 1

        inter_table = {"Exchange": jab}
        spin_table = {"Exchange": spin}

        stub = self._make_uhfk_stub(norb, ns, nvol, ham_trans,
                                    inter_table=inter_table,
                                    spin_table=spin_table,
                                    iflag_fock=False, sz_free=True,
                                    Nconds=[nd])
        stub._detect_blocks()

        blocks = stub.block_info
        # Hartree couples (up,a)-(down,a) but orbitals stay separate (no Fock)
        self.assertEqual(len(blocks), 2,
                         "Exchange Hartree (no Fock) -> 2 orbital-spin pair blocks")
        sizes = sorted([len(b) for b in blocks])
        self.assertEqual(sizes, [2, 2])
        # Each block should pair up/down of same orbital
        for blk in blocks:
            indices = sorted(blk.tolist())
            self.assertEqual(indices[1] - indices[0], norb,
                             "Block should pair (up_a, dn_a)")

    # ----------------------------------------------------------------
    # Pattern 7: CoulombIntra (same-orbital, cross-spin Hartree)
    #   spin_table: spin[0,1,1,0]=1, spin[1,0,0,1]=1
    #   This is purely Hartree (no Fock for CoulombIntra typically)
    #   Couples (up,a)-(down,a) via Hartree
    #   With diagonal transfer -> each orbital pair (up_a, dn_a) forms block
    # ----------------------------------------------------------------
    def test_coulomb_intra_pairs_spins_per_orbital(self):
        """CoulombIntra couples up/down on same orbital.

        norb=3, ns=2, nd=6. Diagonal transfer.
        CoulombIntra: U[a,a] nonzero, spin[0,1,1,0]=1, spin[1,0,0,1]=1.
        Hartree couples (up,a)-(down,a) for each a.
        Expected: 3 blocks of size 2: {0,3}, {1,4}, {2,5}.
        """
        norb, ns, nvol = 3, 2, 4
        nd = norb * ns

        # Diagonal transfer
        ham_trans = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for i in range(nd):
            ham_trans[:, i, i] = 1.0

        # CoulombIntra: U on diagonal only
        uab = np.zeros((nvol, norb, norb), dtype=np.complex128)
        for a in range(norb):
            uab[:, a, a] = 2.0

        spin = np.zeros((2, 2, 2, 2), dtype=int)
        spin[0, 1, 1, 0] = 1
        spin[1, 0, 0, 1] = 1

        inter_table = {"CoulombIntra": uab}
        spin_table = {"CoulombIntra": spin}

        stub = self._make_uhfk_stub(norb, ns, nvol, ham_trans,
                                    inter_table=inter_table,
                                    spin_table=spin_table,
                                    iflag_fock=True, sz_free=True,
                                    Nconds=[nd])
        stub._detect_blocks()

        blocks = stub.block_info
        self.assertEqual(len(blocks), norb,
                         "CoulombIntra + diagonal transfer -> {} orbital-pair blocks".format(norb))
        for blk in blocks:
            self.assertEqual(len(blk), 2,
                             "Each block should be (up_a, dn_a) pair")
            # Verify it's a spin pair for same orbital
            indices = sorted(blk)
            self.assertEqual(indices[1] - indices[0], norb,
                             "Block should pair up/down of same orbital")

    # ----------------------------------------------------------------
    # Pattern 8: CoulombInter with Fock -> mixes orbitals within spin
    #   CoulombInter spin_table: [0,0,0,0]=1, [1,1,1,1]=1,
    #                            [0,1,1,0]=1, [1,0,0,1]=1
    #   Fock: spin_fock[s,t] = max_{u,v} |spin[s,u,t,v]|
    #     spin_fock[0,0] = max(|spin[0,0,0,0]|,|spin[0,1,0,1]|,...) = 1
    #     spin_fock[1,1] = 1, spin_fock[0,1] = |spin[0,0,1,0]|... check all
    #   Actually: spin[0,1,1,0]=1 -> in transpose(0,2,1,3) = spin_t[0,1,1,0]
    #     spin_fock[s,t] is computed from max over axes (1,3) of |spin|
    #   Need to check: spin_fock = max over (u,v) of |spin[s,u,t,v]|
    #   For CoulombInter: spin_fock[0,1] = max_u,v |spin[0,u,1,v]|
    #     = max(|spin[0,0,1,0]|, |spin[0,0,1,1]|, |spin[0,1,1,0]|, |spin[0,1,1,1]|)
    #     = |spin[0,1,1,0]| = 1
    #   So Fock term of CoulombInter DOES couple spins! -> single block
    # ----------------------------------------------------------------
    def test_coulomb_inter_fock_mixes_spins(self):
        """CoulombInter Fock term couples different spins.

        norb=2, ns=2, nd=4.
        spin_fock[0,1] != 0 because spin[0,1,1,0]=1 -> single block.
        """
        norb, ns, nvol = 2, 2, 4
        nd = norb * ns

        ham_trans = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for i in range(nd):
            ham_trans[:, i, i] = 1.0

        # CoulombInter with inter-orbital terms
        vab = np.zeros((nvol, norb, norb), dtype=np.complex128)
        vab[:, 0, 1] = 1.0
        vab[:, 1, 0] = 1.0

        spin = np.zeros((2, 2, 2, 2), dtype=int)
        spin[0, 0, 0, 0] = 1
        spin[1, 1, 1, 1] = 1
        spin[0, 1, 1, 0] = 1
        spin[1, 0, 0, 1] = 1

        inter_table = {"CoulombInter": vab}
        spin_table = {"CoulombInter": spin}

        stub = self._make_uhfk_stub(norb, ns, nvol, ham_trans,
                                    inter_table=inter_table,
                                    spin_table=spin_table,
                                    iflag_fock=True, sz_free=True,
                                    Nconds=[nd])
        stub._detect_blocks()

        blocks = stub.block_info
        self.assertEqual(len(blocks), 1,
                         "CoulombInter Fock mixes spins and orbitals -> single block")

    # ----------------------------------------------------------------
    # Pattern 9: Transfer couples some orbitals, not others
    #   3 orbitals, ns=2: orb 0 and 1 coupled, orb 2 separate
    #   No interactions.
    #   Expected: 4 blocks: {up_0, up_1}, {dn_0, dn_1}, {up_2}, {dn_2}
    # ----------------------------------------------------------------
    def test_partial_orbital_mixing_transfer_only(self):
        """Transfer couples orbitals 0-1 but not orbital 2.

        norb=3, ns=2, nd=6.
        Expected: 4 blocks: {0,1}, {3,4}, {2}, {5}.
        """
        norb, ns, nvol = 3, 2, 4
        nd = norb * ns
        rng = np.random.RandomState(40)

        ham_trans = np.zeros((nvol, nd, nd), dtype=np.complex128)
        # Spin-up: orb 0-1 coupled
        ham_trans[:, 0, 1] = rng.randn(nvol) + 1j * rng.randn(nvol)
        ham_trans[:, 1, 0] = ham_trans[:, 0, 1].conj()
        ham_trans[:, 0, 0] = 1.0
        ham_trans[:, 1, 1] = 1.0
        ham_trans[:, 2, 2] = 1.0  # orb 2 isolated
        # Spin-down: same pattern
        ham_trans[:, 3, 4] = rng.randn(nvol) + 1j * rng.randn(nvol)
        ham_trans[:, 4, 3] = ham_trans[:, 3, 4].conj()
        ham_trans[:, 3, 3] = 1.0
        ham_trans[:, 4, 4] = 1.0
        ham_trans[:, 5, 5] = 1.0  # orb 2 isolated

        stub = self._make_uhfk_stub(norb, ns, nvol, ham_trans, sz_free=True,
                                    Nconds=[nd])
        stub._detect_blocks()

        blocks = stub.block_info
        self.assertEqual(len(blocks), 4,
                         "Partial orbital mixing -> 4 blocks")
        sizes = sorted([len(b) for b in blocks])
        self.assertEqual(sizes, [1, 1, 2, 2])

    # ----------------------------------------------------------------
    # Pattern 10: Combined - partial orbital mixing + CoulombIntra
    #   Transfer couples orb 0-1 within each spin.
    #   CoulombIntra couples (up,a)-(dn,a) for all a.
    #   Expected: 2 blocks: {up_0, up_1, dn_0, dn_1} and {up_2, dn_2}
    # ----------------------------------------------------------------
    def test_partial_orbital_plus_coulomb_intra(self):
        """Transfer couples orbs 0-1, CoulombIntra couples up/down.

        norb=3, ns=2, nd=6.
        Transfer: orb 0-1 coupled within each spin, orb 2 separate.
        CoulombIntra: (up,a)-(dn,a) for all a.
        Expected: 2 blocks: {0,1,3,4} and {2,5}.
        """
        norb, ns, nvol = 3, 2, 4
        nd = norb * ns

        ham_trans = np.zeros((nvol, nd, nd), dtype=np.complex128)
        # Spin-up: orb 0-1 coupled
        ham_trans[:, 0, 1] = 0.5
        ham_trans[:, 1, 0] = 0.5
        for i in range(nd):
            ham_trans[:, i, i] = 1.0
        # Spin-down: orb 0-1 coupled
        ham_trans[:, 3, 4] = 0.5
        ham_trans[:, 4, 3] = 0.5

        # CoulombIntra: U on all diagonal orbitals
        uab = np.zeros((nvol, norb, norb), dtype=np.complex128)
        for a in range(norb):
            uab[:, a, a] = 2.0

        spin = np.zeros((2, 2, 2, 2), dtype=int)
        spin[0, 1, 1, 0] = 1
        spin[1, 0, 0, 1] = 1

        inter_table = {"CoulombIntra": uab}
        spin_table = {"CoulombIntra": spin}

        stub = self._make_uhfk_stub(norb, ns, nvol, ham_trans,
                                    inter_table=inter_table,
                                    spin_table=spin_table,
                                    iflag_fock=True, sz_free=True,
                                    Nconds=[nd])
        stub._detect_blocks()

        blocks = stub.block_info
        self.assertEqual(len(blocks), 2,
                         "Partial orb mixing + CoulombIntra -> 2 blocks")
        sizes = sorted([len(b) for b in blocks])
        self.assertEqual(sizes, [2, 4])
        # The larger block should contain indices {0,1,3,4}
        big_block = [b for b in blocks if len(b) == 4][0]
        self.assertEqual(sorted(big_block.tolist()), [0, 1, 3, 4])
        small_block = [b for b in blocks if len(b) == 2][0]
        self.assertEqual(sorted(small_block.tolist()), [2, 5])

    # ----------------------------------------------------------------
    # Pattern 11: Spin-orbital mode (ns=1)
    #   When enable_spin_orbital=True, ns=1 and norb already includes spin.
    #   Block structure is purely from transfer connectivity.
    # ----------------------------------------------------------------
    def test_spin_orbital_mode_block_from_transfer(self):
        """Spin-orbital mode (ns=1): blocks determined by transfer only.

        norb=4 (includes spin), ns=1, nd=4.
        Transfer: orb 0-1 coupled, orb 2-3 coupled.
        Expected: 2 blocks {0,1} and {2,3}.
        """
        norb, ns, nvol = 4, 1, 4
        nd = norb * ns  # 4

        ham_trans = np.zeros((nvol, nd, nd), dtype=np.complex128)
        # Group 1: orb 0-1
        ham_trans[:, 0, 0] = 1.0
        ham_trans[:, 1, 1] = 1.0
        ham_trans[:, 0, 1] = 0.3
        ham_trans[:, 1, 0] = 0.3
        # Group 2: orb 2-3
        ham_trans[:, 2, 2] = 1.0
        ham_trans[:, 3, 3] = 1.0
        ham_trans[:, 2, 3] = 0.3
        ham_trans[:, 3, 2] = 0.3

        stub = self._make_uhfk_stub(norb, ns, nvol, ham_trans,
                                    sz_free=True, Nconds=[2],
                                    enable_spin_orbital=True)
        stub._detect_blocks()

        blocks = stub.block_info
        self.assertEqual(len(blocks), 2, "Spin-orbital: 2 groups from transfer")
        sizes = sorted([len(b) for b in blocks])
        self.assertEqual(sizes, [2, 2])

    # ----------------------------------------------------------------
    # Pattern 12: PairLift interaction (cross-spin, same-orbital Fock)
    #   PairLift: spin[0,0,1,1]=1, spin[1,1,0,0]=1
    #   Fock: spin_fock[s,t] = max_{u,v} |spin[s,u,t,v]|
    #     spin_fock[0,1] = max(|spin[0,0,1,0]|, |spin[0,0,1,1]|,
    #                          |spin[0,1,1,0]|, |spin[0,1,1,1]|) = 1
    #   So Fock connects (up,a)-(down,b) if J[a,b] != 0
    # ----------------------------------------------------------------
    def test_pairlift_fock_mixes_spins(self):
        """PairLift Fock term mixes spins.

        norb=2, ns=2, nd=4.
        PairLift with J[0,1] != 0 -> Fock couples up-down.
        Expected: single block.
        """
        norb, ns, nvol = 2, 2, 4
        nd = norb * ns

        ham_trans = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for i in range(nd):
            ham_trans[:, i, i] = 1.0

        jab = np.zeros((nvol, norb, norb), dtype=np.complex128)
        jab[:, 0, 1] = 1.0
        jab[:, 1, 0] = 1.0

        spin = np.zeros((2, 2, 2, 2), dtype=int)
        spin[0, 0, 1, 1] = 1
        spin[1, 1, 0, 0] = 1

        inter_table = {"PairLift": jab}
        spin_table = {"PairLift": spin}

        stub = self._make_uhfk_stub(norb, ns, nvol, ham_trans,
                                    inter_table=inter_table,
                                    spin_table=spin_table,
                                    iflag_fock=True, sz_free=True,
                                    Nconds=[nd])
        stub._detect_blocks()

        blocks = stub.block_info
        self.assertEqual(len(blocks), 1,
                         "PairLift Fock mixes spins -> single block")

    # ----------------------------------------------------------------
    # Pattern 13: Multiple interactions combined
    #   Hund (same-spin) + CoulombIntra (cross-spin same-orbital)
    #   norb=2, ns=2, nd=4
    #   Hund Fock: couples (s,a)-(s,b) -> merges orbs within spin
    #   CoulombIntra Hartree: couples (up,a)-(dn,a) -> merges spins per orbital
    #   Together: all connected -> single block
    # ----------------------------------------------------------------
    def test_hund_plus_coulomb_intra_single_block(self):
        """Hund + CoulombIntra together connect all indices.

        Hund merges orbitals within each spin.
        CoulombIntra bridges spins on same orbital.
        Combined: single block.
        """
        norb, ns, nvol = 2, 2, 4
        nd = norb * ns

        ham_trans = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for i in range(nd):
            ham_trans[:, i, i] = 1.0

        # Hund
        hund_jab = np.zeros((nvol, norb, norb), dtype=np.complex128)
        hund_jab[:, 0, 1] = 1.0
        hund_jab[:, 1, 0] = 1.0
        hund_spin = np.zeros((2, 2, 2, 2), dtype=int)
        hund_spin[0, 0, 0, 0] = 1
        hund_spin[1, 1, 1, 1] = 1

        # CoulombIntra
        uab = np.zeros((nvol, norb, norb), dtype=np.complex128)
        uab[:, 0, 0] = 2.0
        uab[:, 1, 1] = 2.0
        ci_spin = np.zeros((2, 2, 2, 2), dtype=int)
        ci_spin[0, 1, 1, 0] = 1
        ci_spin[1, 0, 0, 1] = 1

        inter_table = {"Hund": hund_jab, "CoulombIntra": uab}
        spin_table = {"Hund": hund_spin, "CoulombIntra": ci_spin}

        stub = self._make_uhfk_stub(norb, ns, nvol, ham_trans,
                                    inter_table=inter_table,
                                    spin_table=spin_table,
                                    iflag_fock=True, sz_free=True,
                                    Nconds=[nd])
        stub._detect_blocks()

        blocks = stub.block_info
        self.assertEqual(len(blocks), 1,
                         "Hund + CoulombIntra -> fully connected -> single block")

    # ----------------------------------------------------------------
    # Pattern 14: Ncond distribution for sz_free blocks
    # ----------------------------------------------------------------
    def test_ncond_distribution_sz_free(self):
        """Check Ncond is correctly distributed for sz_free mode.

        norb=3, ns=2, nd=6. Transfer fully couples all orbs within each spin.
        Creates 2 spin blocks of size 3 each.
        Total Ncond=4, should split 2+2.
        """
        norb, ns, nvol = 3, 2, 4
        nd = norb * ns
        rng = np.random.RandomState(50)

        ham_trans = np.zeros((nvol, nd, nd), dtype=np.complex128)
        # Fully couple all orbitals within spin-up
        blk = rng.randn(nvol, norb, norb) + 1j * rng.randn(nvol, norb, norb)
        ham_trans[:, :norb, :norb] = blk
        # Fully couple all orbitals within spin-down
        blk = rng.randn(nvol, norb, norb) + 1j * rng.randn(nvol, norb, norb)
        ham_trans[:, norb:, norb:] = blk

        stub = self._make_uhfk_stub(norb, ns, nvol, ham_trans,
                                    sz_free=True, Nconds=[4])
        stub._detect_blocks()

        blocks = stub.block_info
        self.assertEqual(len(blocks), 2,
                         "Spin-diagonal transfer -> 2 spin blocks")
        self.assertEqual(sum(stub.Nconds), 4,
                         "Total Ncond must be preserved")
        # Equal-size blocks should get equal Ncond
        self.assertEqual(stub.Nconds[0], 2)
        self.assertEqual(stub.Nconds[1], 2)

    # ----------------------------------------------------------------
    # Pattern 15: Ncond distribution for 2Sz fixed (sz_free=False)
    # ----------------------------------------------------------------
    def test_ncond_distribution_sz_fixed(self):
        """Check Ncond is correctly set for 2Sz fixed mode.

        norb=2, ns=2, nd=4. Spin-diagonal transfer -> 2 spin blocks.
        Nconds=[2, 1] (2 up electrons, 1 down electron).
        """
        norb, ns, nvol = 2, 2, 4
        nd = norb * ns

        ham_trans = np.zeros((nvol, nd, nd), dtype=np.complex128)
        ham_trans[:, :norb, :norb] = np.eye(norb)
        ham_trans[:, norb:, norb:] = np.eye(norb)
        ham_trans[:, 0, 1] = 0.3
        ham_trans[:, 1, 0] = 0.3
        ham_trans[:, 2, 3] = 0.3
        ham_trans[:, 3, 2] = 0.3

        stub = self._make_uhfk_stub(norb, ns, nvol, ham_trans,
                                    sz_free=False, Nconds=[2, 1])
        stub._detect_blocks()

        blocks = stub.block_info
        self.assertEqual(len(blocks), 2,
                         "Spin-diagonal transfer -> 2 spin blocks")

        # Check Ncond assignment matches spin sectors
        for i, blk in enumerate(blocks):
            all_up = all(idx < norb for idx in blk)
            all_dn = all(idx >= norb for idx in blk)
            if all_up:
                self.assertEqual(stub.Nconds[i], 2,
                                 "Spin-up block Ncond should be 2")
            elif all_dn:
                self.assertEqual(stub.Nconds[i], 1,
                                 "Spin-down block Ncond should be 1")


    # ----------------------------------------------------------------
    # Pattern 16: Partial orbital mixing via interaction only
    #   Transfer: diagonal (no orbital coupling).
    #   Hund (same-spin Fock): J[0,1] != 0 couples orb 0-1 within each spin.
    #   Orbital 2 has no interaction -> stays separate.
    #   Expected: 4 blocks: {up_0, up_1}, {up_2}, {dn_0, dn_1}, {dn_2}
    # ----------------------------------------------------------------
    def test_partial_orbital_mixing_via_interaction(self):
        """Hund interaction couples orbitals 0-1 but not orbital 2.

        norb=3, ns=2, nd=6. Diagonal transfer.
        Hund J[0,1] != 0: Fock couples (s,0)-(s,1) for same spin.
        Expected: 4 blocks: {0,1}, {2}, {3,4}, {5}.
        """
        norb, ns, nvol = 3, 2, 4
        nd = norb * ns

        ham_trans = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for i in range(nd):
            ham_trans[:, i, i] = 1.0

        # Hund: only J[0,1] and J[1,0] nonzero
        hund_jab = np.zeros((nvol, norb, norb), dtype=np.complex128)
        hund_jab[:, 0, 1] = 1.0
        hund_jab[:, 1, 0] = 1.0

        hund_spin = np.zeros((2, 2, 2, 2), dtype=int)
        hund_spin[0, 0, 0, 0] = 1
        hund_spin[1, 1, 1, 1] = 1

        inter_table = {"Hund": hund_jab}
        spin_table = {"Hund": hund_spin}

        stub = self._make_uhfk_stub(norb, ns, nvol, ham_trans,
                                    inter_table=inter_table,
                                    spin_table=spin_table,
                                    iflag_fock=True, sz_free=True,
                                    Nconds=[nd])
        stub._detect_blocks()

        blocks = stub.block_info
        self.assertEqual(len(blocks), 4,
                         "Hund on orbs 0-1 only + diagonal transfer -> 4 blocks")
        sizes = sorted([len(b) for b in blocks])
        self.assertEqual(sizes, [1, 1, 2, 2])

    # ----------------------------------------------------------------
    # Pattern 17: Chain-like orbital connectivity
    #   Transfer: orb 0-1 coupled and orb 1-2 coupled -> all 3 connected
    #   Orb 3 isolated.
    #   Expected: 2 blocks per spin (size 3 and 1), total 4 blocks.
    # ----------------------------------------------------------------
    def test_chain_orbital_connectivity(self):
        """Transfer creates chain: orb0-orb1-orb2, orb3 isolated.

        norb=4, ns=2, nd=8.
        Expected: 4 blocks: {0,1,2}, {3}, {4,5,6}, {7}.
        """
        norb, ns, nvol = 4, 2, 4
        nd = norb * ns

        ham_trans = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for i in range(nd):
            ham_trans[:, i, i] = 1.0
        # Spin-up chain: 0-1-2
        ham_trans[:, 0, 1] = 0.3
        ham_trans[:, 1, 0] = 0.3
        ham_trans[:, 1, 2] = 0.3
        ham_trans[:, 2, 1] = 0.3
        # Spin-down chain: 4-5-6
        ham_trans[:, 4, 5] = 0.3
        ham_trans[:, 5, 4] = 0.3
        ham_trans[:, 5, 6] = 0.3
        ham_trans[:, 6, 5] = 0.3

        stub = self._make_uhfk_stub(norb, ns, nvol, ham_trans,
                                    sz_free=True, Nconds=[nd])
        stub._detect_blocks()

        blocks = stub.block_info
        self.assertEqual(len(blocks), 4,
                         "Chain connectivity + isolated orb -> 4 blocks")
        sizes = sorted([len(b) for b in blocks])
        self.assertEqual(sizes, [1, 1, 3, 3])

    # ----------------------------------------------------------------
    # Pattern 18: Partial orbital mixing via interaction bridging transfer gaps
    #   Transfer: orb 0 and orb 1 separate (diagonal).
    #   CoulombInter (Fock): J[0,1] != 0 -> couples (s,0)-(s,1) within same spin.
    #   CoulombIntra: U[a,a] -> couples (up,a)-(dn,a).
    #   Together: Fock bridges orb 0-1, CoulombIntra bridges spins
    #   -> single block.
    # ----------------------------------------------------------------
    def test_interaction_bridges_orbital_gap(self):
        """CoulombInter Fock bridges orbital gap that transfer doesn't.

        norb=2, ns=2, nd=4. Diagonal transfer.
        CoulombInter Fock: J[0,1] couples (s,0)-(s,1) per spin.
        CoulombIntra: couples (up,a)-(dn,a).
        Combined: all indices connected -> single block.
        """
        norb, ns, nvol = 2, 2, 4
        nd = norb * ns

        ham_trans = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for i in range(nd):
            ham_trans[:, i, i] = 1.0

        # CoulombInter: J[0,1] != 0
        vab = np.zeros((nvol, norb, norb), dtype=np.complex128)
        vab[:, 0, 1] = 1.0
        vab[:, 1, 0] = 1.0
        ci_spin = np.zeros((2, 2, 2, 2), dtype=int)
        ci_spin[0, 0, 0, 0] = 1
        ci_spin[1, 1, 1, 1] = 1
        ci_spin[0, 1, 1, 0] = 1
        ci_spin[1, 0, 0, 1] = 1

        # CoulombIntra: U on diagonal
        uab = np.zeros((nvol, norb, norb), dtype=np.complex128)
        uab[:, 0, 0] = 2.0
        uab[:, 1, 1] = 2.0
        cu_spin = np.zeros((2, 2, 2, 2), dtype=int)
        cu_spin[0, 1, 1, 0] = 1
        cu_spin[1, 0, 0, 1] = 1

        inter_table = {"CoulombInter": vab, "CoulombIntra": uab}
        spin_table = {"CoulombInter": ci_spin, "CoulombIntra": cu_spin}

        stub = self._make_uhfk_stub(norb, ns, nvol, ham_trans,
                                    inter_table=inter_table,
                                    spin_table=spin_table,
                                    iflag_fock=True, sz_free=True,
                                    Nconds=[nd])
        stub._detect_blocks()

        blocks = stub.block_info
        self.assertEqual(len(blocks), 1,
                         "CoulombInter Fock + CoulombIntra bridge all gaps -> single block")

    # ----------------------------------------------------------------
    # Pattern 19: Large system with multiple independent orbital groups
    #   norb=6, ns=2, nd=12.
    #   Transfer: orb {0,1,2} coupled, orb {3,4} coupled, orb {5} isolated.
    #   No interactions -> 6 blocks (3 per spin).
    #   With CoulombIntra -> 3 blocks (each spin pair merged).
    # ----------------------------------------------------------------
    def test_multiple_orbital_groups_with_coulomb_intra(self):
        """Multiple orbital groups merged by CoulombIntra across spins.

        norb=6, ns=2, nd=12.
        Transfer groups: {0,1,2}, {3,4}, {5} per spin.
        CoulombIntra pairs up/down -> 3 blocks: sizes 6, 4, 2.
        """
        norb, ns, nvol = 6, 2, 4
        nd = norb * ns

        ham_trans = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for i in range(nd):
            ham_trans[:, i, i] = 1.0
        # Spin-up: group {0,1,2}
        ham_trans[:, 0, 1] = 0.5
        ham_trans[:, 1, 0] = 0.5
        ham_trans[:, 1, 2] = 0.5
        ham_trans[:, 2, 1] = 0.5
        # Spin-up: group {3,4}
        ham_trans[:, 3, 4] = 0.5
        ham_trans[:, 4, 3] = 0.5
        # Spin-down: group {6,7,8}
        ham_trans[:, 6, 7] = 0.5
        ham_trans[:, 7, 6] = 0.5
        ham_trans[:, 7, 8] = 0.5
        ham_trans[:, 8, 7] = 0.5
        # Spin-down: group {9,10}
        ham_trans[:, 9, 10] = 0.5
        ham_trans[:, 10, 9] = 0.5

        # CoulombIntra on all orbitals
        uab = np.zeros((nvol, norb, norb), dtype=np.complex128)
        for a in range(norb):
            uab[:, a, a] = 2.0
        cu_spin = np.zeros((2, 2, 2, 2), dtype=int)
        cu_spin[0, 1, 1, 0] = 1
        cu_spin[1, 0, 0, 1] = 1

        inter_table = {"CoulombIntra": uab}
        spin_table = {"CoulombIntra": cu_spin}

        stub = self._make_uhfk_stub(norb, ns, nvol, ham_trans,
                                    inter_table=inter_table,
                                    spin_table=spin_table,
                                    iflag_fock=True, sz_free=True,
                                    Nconds=[nd])
        stub._detect_blocks()

        blocks = stub.block_info
        self.assertEqual(len(blocks), 3,
                         "3 orbital groups, each merged across spins -> 3 blocks")
        sizes = sorted([len(b) for b in blocks])
        self.assertEqual(sizes, [2, 4, 6])

    # ----------------------------------------------------------------
    # Pattern 20: Partial spin-orbit coupling in transfer
    #   Only some orbitals have spin mixing in transfer; others don't.
    #   norb=3, ns=2, nd=6.
    #   Transfer: orb 0 has spin mixing (up_0 - dn_0 coupled).
    #             orb 1 and 2 are spin-diagonal.
    #   Expected: block {0, 3} (mixed spin), {1}, {2}, {4}, {5}
    #             But orb 1 and 2 are connected within spin if transfer
    #             couples them. Let's make orb 1-2 coupled within spin.
    #   Expected: {0, 3}, {1, 2}, {4, 5}
    # ----------------------------------------------------------------
    def test_partial_spin_orbit_in_transfer(self):
        """Only orbital 0 has spin-orbit mixing in transfer.

        norb=3, ns=2, nd=6.
        Transfer: orb 0 mixes spins (up_0 <-> dn_0).
                  orb 1-2 coupled within spin only.
        Expected: 3 blocks: {0,3} (spin-mixed), {1,2} (up), {4,5} (dn).
        """
        norb, ns, nvol = 3, 2, 4
        nd = norb * ns

        ham_trans = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for i in range(nd):
            ham_trans[:, i, i] = 1.0
        # Orb 0: spin mixing (up_0 <-> dn_0)
        ham_trans[:, 0, 3] = 0.5  # up_0 to dn_0
        ham_trans[:, 3, 0] = 0.5  # dn_0 to up_0
        # Orb 1-2: coupled within spin-up
        ham_trans[:, 1, 2] = 0.3
        ham_trans[:, 2, 1] = 0.3
        # Orb 1-2: coupled within spin-down
        ham_trans[:, 4, 5] = 0.3
        ham_trans[:, 5, 4] = 0.3

        stub = self._make_uhfk_stub(norb, ns, nvol, ham_trans,
                                    sz_free=True, Nconds=[nd])
        stub._detect_blocks()

        blocks = stub.block_info
        self.assertEqual(len(blocks), 3,
                         "Partial spin-orbit -> 3 blocks")
        sizes = sorted([len(b) for b in blocks])
        self.assertEqual(sizes, [2, 2, 2])

        # Find the spin-mixed block (contains both up and dn indices)
        for blk in blocks:
            indices = sorted(blk.tolist())
            if 0 in indices:
                self.assertIn(3, indices,
                              "Orb 0 up/dn should be in same block")
            if 1 in indices:
                self.assertIn(2, indices,
                              "Orbs 1,2 in spin-up should be together")
            if 4 in indices:
                self.assertIn(5, indices,
                              "Orbs 1,2 in spin-down should be together")


class TestUHFrDetectBlocks(unittest.TestCase):
    """Test _detect_blocks in UHFr for various Hamiltonian patterns.

    UHFr uses flat indices 0..2*Nsize-1.
    For spin-up/down split: up=[0..Nsize-1], down=[Nsize..2*Nsize-1].
    """

    def _make_uhfr_stub(self, Nsize, Ham_trans, Ham_local,
                        TwoSz=None, Ncond=None):
        """Create a minimal UHFr-like object for _detect_blocks testing.

        Parameters
        ----------
        Nsize : int
            Number of sites.
        Ham_trans : ndarray (2*Nsize, 2*Nsize)
            Transfer Hamiltonian.
        Ham_local : ndarray (2*Nsize, 2*Nsize, 2*Nsize, 2*Nsize)
            Interaction matrix (before reshape to (4N^2, 4N^2)).
        TwoSz : int or None
            If None, sz-free mode. Otherwise spin-up/down split.
        Ncond : int or None
            Electron count. Defaults to Nsize.
        """
        import hwave.solver.uhfr as uhfr_module

        nd = 2 * Nsize
        stub = object.__new__(uhfr_module.UHFr)
        stub.Nsize = Nsize
        stub.Ham_trans = Ham_trans
        stub.Ham_local = Ham_local.reshape(nd ** 2, nd ** 2)

        if Ncond is None:
            Ncond = Nsize

        if TwoSz is None:
            stub.green_list = {
                "sz-free": {
                    "label": list(range(nd)),
                    "occupied": Ncond
                }
            }
        else:
            stub.green_list = {
                "spin-up": {
                    "label": list(range(Nsize)),
                    "value": 0.5,
                    "occupied": int((Ncond + TwoSz) / 2)
                },
                "spin-down": {
                    "label": list(range(Nsize, nd)),
                    "value": -0.5,
                    "occupied": int((Ncond - TwoSz) / 2)
                }
            }
        return stub

    def test_spin_diagonal_transfer(self):
        """Transfer only couples within spin sectors.

        Nsize=3, nd=6. Transfer: up-up and dn-dn coupled, no cross-spin.
        TwoSz specified -> already 2 spin blocks from green_list.
        No further splitting expected since orbitals all connected within spin.
        """
        Nsize = 3
        nd = 2 * Nsize
        rng = np.random.RandomState(300)

        Ham_trans = np.zeros((nd, nd), dtype=complex)
        # Spin-up: all sites connected
        blk = rng.randn(Nsize, Nsize) + 1j * rng.randn(Nsize, Nsize)
        Ham_trans[:Nsize, :Nsize] = (blk + blk.conj().T) / 2
        # Spin-down: all sites connected
        blk = rng.randn(Nsize, Nsize) + 1j * rng.randn(Nsize, Nsize)
        Ham_trans[Nsize:, Nsize:] = (blk + blk.conj().T) / 2

        Ham_local = np.zeros((nd, nd, nd, nd), dtype=complex)

        stub = self._make_uhfr_stub(Nsize, Ham_trans, Ham_local,
                                    TwoSz=1, Ncond=3)
        stub._detect_blocks()

        # Should remain as 2 spin blocks
        self.assertEqual(len(stub.green_list), 2)

    def test_diagonal_transfer_orbital_split(self):
        """Diagonal transfer splits into individual site blocks.

        Nsize=3, nd=6, sz-free mode. Diagonal transfer only.
        No interactions -> 6 independent blocks.
        """
        Nsize = 3
        nd = 2 * Nsize

        Ham_trans = np.diag(np.arange(1, nd + 1, dtype=complex))
        Ham_local = np.zeros((nd, nd, nd, nd), dtype=complex)

        stub = self._make_uhfr_stub(Nsize, Ham_trans, Ham_local,
                                    TwoSz=None, Ncond=3)
        stub._detect_blocks()

        self.assertEqual(len(stub.green_list), nd,
                         "Diagonal transfer -> {} blocks".format(nd))

    def test_interaction_connects_spins(self):
        """Interaction term connects spin-up and spin-down sites.

        Nsize=2, nd=4. Diagonal transfer.
        Ham_local couples (0,2) and (1,3) -> merges spin sectors per site.
        TwoSz=None -> starts as single block.
        Expected: 2 blocks {0,2} and {1,3}.
        """
        Nsize = 2
        nd = 2 * Nsize

        Ham_trans = np.diag(np.ones(nd, dtype=complex))
        Ham_local = np.zeros((nd, nd, nd, nd), dtype=complex)
        # Couple site 0 (up) with site 2 (down) via interaction
        Ham_local[0, 0, 2, 2] = 1.0
        Ham_local[2, 2, 0, 0] = 1.0
        # Couple site 1 (up) with site 3 (down) via interaction
        Ham_local[1, 1, 3, 3] = 1.0
        Ham_local[3, 3, 1, 1] = 1.0

        stub = self._make_uhfr_stub(Nsize, Ham_trans, Ham_local,
                                    TwoSz=None, Ncond=2)
        stub._detect_blocks()

        self.assertEqual(len(stub.green_list), 2,
                         "Interaction connects spin pairs -> 2 blocks")
        for k, info in stub.green_list.items():
            self.assertEqual(len(info["label"]), 2)

    def test_partial_transfer_coupling(self):
        """Transfer couples sites 0-1 within spin-up, 2-3 within spin-down.

        Nsize=3, nd=6. Site 2 (up) and 5 (dn) are isolated.
        TwoSz specified -> start with 2 spin blocks.
        Expected: 4 blocks: {0,1}, {2}, {3,4}, {5}.
        """
        Nsize = 3
        nd = 2 * Nsize

        Ham_trans = np.zeros((nd, nd), dtype=complex)
        # Spin-up: site 0-1 coupled
        Ham_trans[0, 0] = 1.0
        Ham_trans[1, 1] = 1.0
        Ham_trans[2, 2] = 1.0
        Ham_trans[0, 1] = 0.5
        Ham_trans[1, 0] = 0.5
        # Spin-down: site 3-4 coupled
        Ham_trans[3, 3] = 1.0
        Ham_trans[4, 4] = 1.0
        Ham_trans[5, 5] = 1.0
        Ham_trans[3, 4] = 0.5
        Ham_trans[4, 3] = 0.5

        Ham_local = np.zeros((nd, nd, nd, nd), dtype=complex)

        stub = self._make_uhfr_stub(Nsize, Ham_trans, Ham_local,
                                    TwoSz=1, Ncond=3)
        stub._detect_blocks()

        self.assertEqual(len(stub.green_list), 4,
                         "Partial coupling -> 4 sub-blocks")
        sizes = sorted([len(v["label"]) for v in stub.green_list.values()])
        self.assertEqual(sizes, [1, 1, 2, 2])

    def test_full_coupling_single_block(self):
        """Fully coupled transfer -> no splitting.

        Nsize=3, nd=6. All sites connected via transfer.
        sz-free mode -> single block.
        """
        Nsize = 3
        nd = 2 * Nsize
        rng = np.random.RandomState(310)

        Ham_trans = rng.randn(nd, nd) + 1j * rng.randn(nd, nd)
        Ham_trans = (Ham_trans + Ham_trans.conj().T) / 2
        Ham_local = np.zeros((nd, nd, nd, nd), dtype=complex)

        stub = self._make_uhfr_stub(Nsize, Ham_trans, Ham_local,
                                    TwoSz=None, Ncond=3)
        stub._detect_blocks()

        self.assertEqual(len(stub.green_list), 1,
                         "Fully coupled -> single block")

    def test_ncond_preserved_after_split(self):
        """Total occupied count is preserved after block splitting.

        Nsize=4, nd=8. sz-free mode with Ncond=6.
        Diagonal transfer -> 8 blocks with Ncond summing to 6.
        """
        Nsize = 4
        nd = 2 * Nsize

        Ham_trans = np.diag(np.ones(nd, dtype=complex))
        Ham_local = np.zeros((nd, nd, nd, nd), dtype=complex)

        stub = self._make_uhfr_stub(Nsize, Ham_trans, Ham_local,
                                    TwoSz=None, Ncond=6)
        stub._detect_blocks()

        total_occ = sum(v["occupied"] for v in stub.green_list.values())
        self.assertEqual(total_occ, 6,
                         "Total Ncond must be preserved after splitting")


class TestUHFkSpinOrbitalInteraction(unittest.TestCase):
    """Test spin-orbital mode with interaction terms.

    Verifies that running UHFk in spin-orbital mode (ns=1, nd=2*norb_phys)
    with interactions produces results equivalent to normal mode (ns=2, norb=norb_phys).

    In spin-orbital mode, the orbital index uses the convention:
        index = 2 * physical_orbital + spin  (spin=0 for up, 1 for down)
    """

    def _make_uhfk_stub_full(self, norb_phys, nvol, ham_trans_normal,
                             inter_table=None, spin_table=None,
                             iflag_fock=True, sz_free=True, Nconds=None,
                             green_init=None):
        """Create a UHFk-like stub with both normal and spin-orbital modes.

        Parameters
        ----------
        norb_phys : int
            Number of physical orbitals.
        nvol : int
            Number of k-points.
        ham_trans_normal : ndarray (nvol, 2*norb_phys, 2*norb_phys)
            Transfer Hamiltonian in normal mode (s*norb + a ordering).
        inter_table : dict or None
            Interaction tables {name: ndarray(nvol, norb_phys, norb_phys) or None}.
        spin_table : dict or None
            Spin combination tables {name: ndarray(2,2,2,2)}.
        iflag_fock : bool
            Whether Fock term is enabled.
        green_init : ndarray or None
            Initial Green function in normal mode (nvol, 2, norb_phys, 2, norb_phys).

        Returns
        -------
        stub_normal, stub_so : pair of UHFk-like stubs
        """
        import hwave.solver.uhfk as uhfk_module

        nd = 2 * norb_phys

        # --- Normal mode stub ---
        stub_n = object.__new__(uhfk_module.UHFk)
        stub_n.norb = norb_phys
        stub_n.ns = 2
        stub_n.nd = nd
        stub_n.nvol = nvol
        stub_n.shape = (nvol, 1, 1)  # (nx, ny, nz) for FFT
        stub_n.ham_trans = ham_trans_normal
        # Ensure all interaction types exist in table (defaulting to None)
        all_types = ['CoulombIntra', 'CoulombInter', 'Hund', 'Ising',
                     'PairLift', 'Exchange', 'PairHop']
        _inter = {t: None for t in all_types}
        _spin = {}
        if inter_table is not None:
            _inter.update(inter_table)
        if spin_table is not None:
            _spin.update(spin_table)
        stub_n.inter_table = _inter
        stub_n.spin_table = _spin
        stub_n.iflag_fock = iflag_fock
        stub_n.sz_free = sz_free
        stub_n.enable_spin_orbital = False
        stub_n.norb_phys = norb_phys
        stub_n.block_info = [np.arange(nd)]
        stub_n.Nconds = Nconds if Nconds is not None else [nd]
        stub_n.threshold = 1e-12

        # --- Spin-orbital mode stub ---
        # Convert ham_trans from normal ordering (s*norb+a) to SO ordering (2*a+s)
        ham_trans_so = np.zeros_like(ham_trans_normal)
        for s in range(2):
            for t in range(2):
                # Normal: s*norb+a -> SO: 2*a+s
                ham_trans_so[:, s::2, t::2] = ham_trans_normal[:, s*norb_phys:(s+1)*norb_phys,
                                                                t*norb_phys:(t+1)*norb_phys]

        stub_so = object.__new__(uhfk_module.UHFk)
        stub_so.norb = 2 * norb_phys  # nd in SO mode
        stub_so.ns = 1
        stub_so.nd = nd
        stub_so.nvol = nvol
        stub_so.shape = (nvol, 1, 1)
        stub_so.ham_trans = ham_trans_so
        _inter_so = {t: None for t in all_types}
        _spin_so = {}
        if inter_table is not None:
            _inter_so.update(inter_table)
        if spin_table is not None:
            _spin_so.update(spin_table)
        stub_so.inter_table = _inter_so
        stub_so.spin_table = _spin_so
        stub_so.iflag_fock = iflag_fock
        stub_so.sz_free = sz_free
        stub_so.enable_spin_orbital = True
        stub_so.norb_phys = norb_phys
        stub_so.block_info = [np.arange(nd)]
        stub_so.Nconds = Nconds if Nconds is not None else [nd]
        stub_so.threshold = 1e-12

        return stub_n, stub_so

    def _normal_green_to_so(self, green_normal, norb_phys):
        """Convert Green function from normal to spin-orbital ordering.

        normal: (nvol, 2, norb_phys, 2, norb_phys)
        so: (nvol, 1, 2*norb_phys, 1, 2*norb_phys)
        """
        nvol = green_normal.shape[0]
        nd = 2 * norb_phys
        G_so = np.zeros((nvol, 1, nd, 1, nd), dtype=np.complex128)
        for s in range(2):
            for t in range(2):
                G_so[:, 0, s::2, 0, t::2] = green_normal[:, s, :, t, :]
        return G_so

    def _so_ham_to_normal(self, ham_so, norb_phys):
        """Convert Hamiltonian from spin-orbital to normal ordering.

        so: (nvol, nd, nd) with index = 2*a+s
        normal: (nvol, nd, nd) with index = s*norb+a
        """
        nvol = ham_so.shape[0]
        nd = 2 * norb_phys
        H_n = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for s in range(2):
            for t in range(2):
                H_n[:, s*norb_phys:(s+1)*norb_phys,
                     t*norb_phys:(t+1)*norb_phys] = ham_so[:, s::2, t::2]
        return H_n

    # ----------------------------------------------------------------
    # Test 1: CoulombIntra produces identical Hamiltonian in both modes
    # ----------------------------------------------------------------
    def test_coulomb_intra_equivalence(self):
        """CoulombIntra interaction: normal vs spin-orbital mode.

        norb_phys=2, nvol=1.
        CoulombIntra U[0,0]=4.0, U[1,1]=3.0
        """
        norb_phys, nvol = 2, 1
        nd = 2 * norb_phys
        rng = np.random.RandomState(42)

        # Transfer: simple diagonal + small off-diagonal
        ham_t = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for i in range(nd):
            ham_t[:, i, i] = float(i) + 1.0
        ham_t[:, 0, 1] = 0.3
        ham_t[:, 1, 0] = 0.3
        ham_t[:, 2, 3] = 0.3
        ham_t[:, 3, 2] = 0.3

        # CoulombIntra
        uab = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        uab[:, 0, 0] = 4.0
        uab[:, 1, 1] = 3.0
        cu_spin = np.zeros((2, 2, 2, 2), dtype=int)
        cu_spin[0, 1, 1, 0] = 1
        cu_spin[1, 0, 0, 1] = 1

        inter_table = {"CoulombIntra": uab}
        spin_table = {"CoulombIntra": cu_spin}

        stub_n, stub_so = self._make_uhfk_stub_full(
            norb_phys, nvol, ham_t,
            inter_table=inter_table, spin_table=spin_table,
            iflag_fock=True
        )

        # Random Green function in normal mode
        green_n = 0.01 * (rng.randn(nvol, 2, norb_phys, 2, norb_phys)
                          + 1j * rng.randn(nvol, 2, norb_phys, 2, norb_phys))
        # Make it roughly Hermitian-like for diagonal
        for s in range(2):
            for a in range(norb_phys):
                green_n[:, s, a, s, a] = 0.25 + 0.01 * rng.randn(nvol)

        stub_n.Green = green_n
        stub_so.Green = self._normal_green_to_so(green_n, norb_phys)

        # Build Hamiltonians
        stub_n._make_ham()
        stub_so._make_ham()

        # Convert SO Hamiltonian to normal ordering for comparison
        ham_so_in_normal = self._so_ham_to_normal(stub_so.ham, norb_phys)

        np.testing.assert_allclose(
            ham_so_in_normal, stub_n.ham, atol=1e-12,
            err_msg="CoulombIntra: SO and normal mode Hamiltonians should match"
        )

    # ----------------------------------------------------------------
    # Test 2: CoulombInter with Fock term
    # ----------------------------------------------------------------
    def test_coulomb_inter_fock_equivalence(self):
        """CoulombInter interaction with Fock: normal vs spin-orbital mode.

        norb_phys=2, nvol=2.
        CoulombInter V[0,1]=2.0, V[1,0]=2.0
        """
        norb_phys, nvol = 2, 2
        nd = 2 * norb_phys
        rng = np.random.RandomState(123)

        ham_t = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for i in range(nd):
            ham_t[:, i, i] = float(i) + 1.0

        vab = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        vab[:, 0, 1] = 2.0
        vab[:, 1, 0] = 2.0
        ci_spin = np.zeros((2, 2, 2, 2), dtype=int)
        ci_spin[0, 0, 0, 0] = 1
        ci_spin[1, 1, 1, 1] = 1
        ci_spin[0, 1, 1, 0] = 1
        ci_spin[1, 0, 0, 1] = 1

        inter_table = {"CoulombInter": vab}
        spin_table = {"CoulombInter": ci_spin}

        stub_n, stub_so = self._make_uhfk_stub_full(
            norb_phys, nvol, ham_t,
            inter_table=inter_table, spin_table=spin_table,
            iflag_fock=True
        )

        green_n = 0.01 * (rng.randn(nvol, 2, norb_phys, 2, norb_phys)
                          + 1j * rng.randn(nvol, 2, norb_phys, 2, norb_phys))
        for s in range(2):
            for a in range(norb_phys):
                green_n[:, s, a, s, a] = 0.25 + 0.01 * rng.randn(nvol)

        stub_n.Green = green_n
        stub_so.Green = self._normal_green_to_so(green_n, norb_phys)

        stub_n._make_ham()
        stub_so._make_ham()

        ham_so_in_normal = self._so_ham_to_normal(stub_so.ham, norb_phys)
        np.testing.assert_allclose(
            ham_so_in_normal, stub_n.ham, atol=1e-12,
            err_msg="CoulombInter Fock: SO and normal mode Hamiltonians should match"
        )

    # ----------------------------------------------------------------
    # Test 3: Hund interaction (same-spin coupling)
    # ----------------------------------------------------------------
    def test_hund_equivalence(self):
        """Hund interaction: normal vs spin-orbital mode.

        norb_phys=2, nvol=1.
        """
        norb_phys, nvol = 2, 1
        nd = 2 * norb_phys
        rng = np.random.RandomState(77)

        ham_t = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for i in range(nd):
            ham_t[:, i, i] = float(i) + 1.0

        jab = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        jab[:, 0, 1] = -1.5
        jab[:, 1, 0] = -1.5

        hund_spin = np.zeros((2, 2, 2, 2), dtype=int)
        hund_spin[0, 0, 0, 0] = 1
        hund_spin[1, 1, 1, 1] = 1

        inter_table = {"Hund": jab}
        spin_table = {"Hund": hund_spin}

        stub_n, stub_so = self._make_uhfk_stub_full(
            norb_phys, nvol, ham_t,
            inter_table=inter_table, spin_table=spin_table,
            iflag_fock=True
        )

        green_n = 0.01 * (rng.randn(nvol, 2, norb_phys, 2, norb_phys)
                          + 1j * rng.randn(nvol, 2, norb_phys, 2, norb_phys))
        for s in range(2):
            for a in range(norb_phys):
                green_n[:, s, a, s, a] = 0.25 + 0.01 * rng.randn(nvol)

        stub_n.Green = green_n
        stub_so.Green = self._normal_green_to_so(green_n, norb_phys)

        stub_n._make_ham()
        stub_so._make_ham()

        ham_so_in_normal = self._so_ham_to_normal(stub_so.ham, norb_phys)
        np.testing.assert_allclose(
            ham_so_in_normal, stub_n.ham, atol=1e-12,
            err_msg="Hund: SO and normal mode Hamiltonians should match"
        )

    # ----------------------------------------------------------------
    # Test 4: Exchange interaction
    # ----------------------------------------------------------------
    def test_exchange_equivalence(self):
        """Exchange interaction: normal vs spin-orbital mode."""
        norb_phys, nvol = 2, 1
        nd = 2 * norb_phys
        rng = np.random.RandomState(88)

        ham_t = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for i in range(nd):
            ham_t[:, i, i] = float(i) + 1.0

        jab = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        jab[:, 0, 1] = -2.0
        jab[:, 1, 0] = -2.0

        ex_spin = np.zeros((2, 2, 2, 2), dtype=int)
        ex_spin[0, 1, 0, 1] = 1
        ex_spin[1, 0, 1, 0] = 1

        inter_table = {"Exchange": jab}
        spin_table = {"Exchange": ex_spin}

        stub_n, stub_so = self._make_uhfk_stub_full(
            norb_phys, nvol, ham_t,
            inter_table=inter_table, spin_table=spin_table,
            iflag_fock=True
        )

        green_n = 0.01 * (rng.randn(nvol, 2, norb_phys, 2, norb_phys)
                          + 1j * rng.randn(nvol, 2, norb_phys, 2, norb_phys))
        for s in range(2):
            for a in range(norb_phys):
                green_n[:, s, a, s, a] = 0.25 + 0.01 * rng.randn(nvol)

        stub_n.Green = green_n
        stub_so.Green = self._normal_green_to_so(green_n, norb_phys)

        stub_n._make_ham()
        stub_so._make_ham()

        ham_so_in_normal = self._so_ham_to_normal(stub_so.ham, norb_phys)
        np.testing.assert_allclose(
            ham_so_in_normal, stub_n.ham, atol=1e-12,
            err_msg="Exchange: SO and normal mode Hamiltonians should match"
        )

    # ----------------------------------------------------------------
    # Test 5: PairHop interaction
    # ----------------------------------------------------------------
    def test_pairhop_equivalence(self):
        """PairHop interaction: normal vs spin-orbital mode."""
        norb_phys, nvol = 2, 2
        nd = 2 * norb_phys
        rng = np.random.RandomState(99)

        ham_t = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for i in range(nd):
            ham_t[:, i, i] = float(i) + 1.0

        jab = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        jab[:, 0, 1] = 1.0
        jab[:, 1, 0] = 1.0

        ph_spin = np.zeros((2, 2, 2, 2), dtype=int)
        ph_spin[0, 1, 1, 0] = 1
        ph_spin[1, 0, 0, 1] = 1

        inter_table = {"PairHop": jab}
        spin_table = {"PairHop": ph_spin}

        stub_n, stub_so = self._make_uhfk_stub_full(
            norb_phys, nvol, ham_t,
            inter_table=inter_table, spin_table=spin_table,
            iflag_fock=True
        )

        green_n = 0.01 * (rng.randn(nvol, 2, norb_phys, 2, norb_phys)
                          + 1j * rng.randn(nvol, 2, norb_phys, 2, norb_phys))
        for s in range(2):
            for a in range(norb_phys):
                green_n[:, s, a, s, a] = 0.25 + 0.01 * rng.randn(nvol)

        stub_n.Green = green_n
        stub_so.Green = self._normal_green_to_so(green_n, norb_phys)

        stub_n._make_ham()
        stub_so._make_ham()

        ham_so_in_normal = self._so_ham_to_normal(stub_so.ham, norb_phys)
        np.testing.assert_allclose(
            ham_so_in_normal, stub_n.ham, atol=1e-12,
            err_msg="PairHop: SO and normal mode Hamiltonians should match"
        )

    # ----------------------------------------------------------------
    # Test 5b: Ising interaction
    # ----------------------------------------------------------------
    def test_ising_equivalence(self):
        """Ising interaction: normal vs spin-orbital mode.

        Ising spin table: [0,0,0,0]=1, [1,1,1,1]=1, [0,1,1,0]=-1, [1,0,0,1]=-1
        This represents J Sz_i Sz_j (Hartree + Fock).
        """
        norb_phys, nvol = 2, 2
        nd = 2 * norb_phys
        rng = np.random.RandomState(110)

        ham_t = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for i in range(nd):
            ham_t[:, i, i] = float(i) + 1.0

        jab = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        jab[:, 0, 1] = 2.0
        jab[:, 1, 0] = 2.0

        ising_spin = np.zeros((2, 2, 2, 2), dtype=int)
        ising_spin[0, 0, 0, 0] = 1
        ising_spin[1, 1, 1, 1] = 1
        ising_spin[0, 1, 1, 0] = -1
        ising_spin[1, 0, 0, 1] = -1

        inter_table = {"Ising": jab}
        spin_table = {"Ising": ising_spin}

        stub_n, stub_so = self._make_uhfk_stub_full(
            norb_phys, nvol, ham_t,
            inter_table=inter_table, spin_table=spin_table,
            iflag_fock=True
        )

        green_n = 0.01 * (rng.randn(nvol, 2, norb_phys, 2, norb_phys)
                          + 1j * rng.randn(nvol, 2, norb_phys, 2, norb_phys))
        for s in range(2):
            for a in range(norb_phys):
                green_n[:, s, a, s, a] = 0.25 + 0.01 * rng.randn(nvol)

        stub_n.Green = green_n
        stub_so.Green = self._normal_green_to_so(green_n, norb_phys)

        stub_n._make_ham()
        stub_so._make_ham()

        ham_so_in_normal = self._so_ham_to_normal(stub_so.ham, norb_phys)
        np.testing.assert_allclose(
            ham_so_in_normal, stub_n.ham, atol=1e-12,
            err_msg="Ising: SO and normal mode Hamiltonians should match"
        )

    # ----------------------------------------------------------------
    # Test 5c: PairLift interaction
    # ----------------------------------------------------------------
    def test_pairlift_equivalence(self):
        """PairLift interaction: normal vs spin-orbital mode.

        PairLift spin table: [0,0,1,1]=1, [1,1,0,0]=1
        This is a Coulomb-type interaction (Hartree + Fock).
        """
        norb_phys, nvol = 2, 2
        nd = 2 * norb_phys
        rng = np.random.RandomState(120)

        ham_t = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for i in range(nd):
            ham_t[:, i, i] = float(i) + 1.0

        jab = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        jab[:, 0, 1] = 1.5
        jab[:, 1, 0] = 1.5

        pl_spin = np.zeros((2, 2, 2, 2), dtype=int)
        pl_spin[0, 0, 1, 1] = 1
        pl_spin[1, 1, 0, 0] = 1

        inter_table = {"PairLift": jab}
        spin_table = {"PairLift": pl_spin}

        stub_n, stub_so = self._make_uhfk_stub_full(
            norb_phys, nvol, ham_t,
            inter_table=inter_table, spin_table=spin_table,
            iflag_fock=True
        )

        green_n = 0.01 * (rng.randn(nvol, 2, norb_phys, 2, norb_phys)
                          + 1j * rng.randn(nvol, 2, norb_phys, 2, norb_phys))
        for s in range(2):
            for a in range(norb_phys):
                green_n[:, s, a, s, a] = 0.25 + 0.01 * rng.randn(nvol)

        stub_n.Green = green_n
        stub_so.Green = self._normal_green_to_so(green_n, norb_phys)

        stub_n._make_ham()
        stub_so._make_ham()

        ham_so_in_normal = self._so_ham_to_normal(stub_so.ham, norb_phys)
        np.testing.assert_allclose(
            ham_so_in_normal, stub_n.ham, atol=1e-12,
            err_msg="PairLift: SO and normal mode Hamiltonians should match"
        )

    # ----------------------------------------------------------------
    # Test 6: Multiple interactions combined
    # ----------------------------------------------------------------
    def test_combined_interactions_equivalence(self):
        """Multiple interactions: CoulombIntra + CoulombInter + Exchange.

        norb_phys=2, nvol=2. Tests that the combined Hamiltonian matches.
        """
        norb_phys, nvol = 2, 2
        nd = 2 * norb_phys
        rng = np.random.RandomState(55)

        ham_t = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for i in range(nd):
            ham_t[:, i, i] = float(i) + 1.0
        ham_t[:, 0, 1] = 0.5
        ham_t[:, 1, 0] = 0.5
        ham_t[:, 2, 3] = 0.5
        ham_t[:, 3, 2] = 0.5

        # CoulombIntra
        uab = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        uab[:, 0, 0] = 4.0
        uab[:, 1, 1] = 3.0
        cu_spin = np.zeros((2, 2, 2, 2), dtype=int)
        cu_spin[0, 1, 1, 0] = 1
        cu_spin[1, 0, 0, 1] = 1

        # CoulombInter
        vab = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        vab[:, 0, 1] = 2.0
        vab[:, 1, 0] = 2.0
        ci_spin = np.zeros((2, 2, 2, 2), dtype=int)
        ci_spin[0, 0, 0, 0] = 1
        ci_spin[1, 1, 1, 1] = 1
        ci_spin[0, 1, 1, 0] = 1
        ci_spin[1, 0, 0, 1] = 1

        # Exchange
        jex = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        jex[:, 0, 1] = -1.5
        jex[:, 1, 0] = -1.5
        ex_spin = np.zeros((2, 2, 2, 2), dtype=int)
        ex_spin[0, 1, 0, 1] = 1
        ex_spin[1, 0, 1, 0] = 1

        inter_table = {
            "CoulombIntra": uab,
            "CoulombInter": vab,
            "Exchange": jex,
        }
        spin_table = {
            "CoulombIntra": cu_spin,
            "CoulombInter": ci_spin,
            "Exchange": ex_spin,
        }

        stub_n, stub_so = self._make_uhfk_stub_full(
            norb_phys, nvol, ham_t,
            inter_table=inter_table, spin_table=spin_table,
            iflag_fock=True
        )

        green_n = 0.01 * (rng.randn(nvol, 2, norb_phys, 2, norb_phys)
                          + 1j * rng.randn(nvol, 2, norb_phys, 2, norb_phys))
        for s in range(2):
            for a in range(norb_phys):
                green_n[:, s, a, s, a] = 0.25 + 0.01 * rng.randn(nvol)

        stub_n.Green = green_n
        stub_so.Green = self._normal_green_to_so(green_n, norb_phys)

        stub_n._make_ham()
        stub_so._make_ham()

        ham_so_in_normal = self._so_ham_to_normal(stub_so.ham, norb_phys)
        np.testing.assert_allclose(
            ham_so_in_normal, stub_n.ham, atol=1e-12,
            err_msg="Combined interactions: SO and normal mode Hamiltonians should match"
        )

    # ----------------------------------------------------------------
    # Test 7: Energy calculation equivalence
    # ----------------------------------------------------------------
    def test_energy_equivalence(self):
        """Interaction energy: normal vs spin-orbital mode.

        Uses CoulombIntra + CoulombInter with Fock.
        """
        norb_phys, nvol = 2, 1
        nd = 2 * norb_phys
        rng = np.random.RandomState(33)

        ham_t = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for i in range(nd):
            ham_t[:, i, i] = float(i) + 1.0

        # CoulombIntra
        uab = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        uab[:, 0, 0] = 4.0
        uab[:, 1, 1] = 3.0
        cu_spin = np.zeros((2, 2, 2, 2), dtype=int)
        cu_spin[0, 1, 1, 0] = 1
        cu_spin[1, 0, 0, 1] = 1

        # CoulombInter
        vab = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        vab[:, 0, 1] = 2.0
        vab[:, 1, 0] = 2.0
        ci_spin = np.zeros((2, 2, 2, 2), dtype=int)
        ci_spin[0, 0, 0, 0] = 1
        ci_spin[1, 1, 1, 1] = 1
        ci_spin[0, 1, 1, 0] = 1
        ci_spin[1, 0, 0, 1] = 1

        inter_table = {"CoulombIntra": uab, "CoulombInter": vab}
        spin_table = {"CoulombIntra": cu_spin, "CoulombInter": ci_spin}

        stub_n, stub_so = self._make_uhfk_stub_full(
            norb_phys, nvol, ham_t,
            inter_table=inter_table, spin_table=spin_table,
            iflag_fock=True
        )

        green_n = 0.01 * (rng.randn(nvol, 2, norb_phys, 2, norb_phys)
                          + 1j * rng.randn(nvol, 2, norb_phys, 2, norb_phys))
        for s in range(2):
            for a in range(norb_phys):
                green_n[:, s, a, s, a] = 0.25 + 0.01 * rng.randn(nvol)

        stub_n.Green = green_n
        stub_so.Green = self._normal_green_to_so(green_n, norb_phys)

        # Build Hamiltonians (needed for band energy setup)
        stub_n._make_ham()
        stub_so._make_ham()

        # Setup minimal eigenvalue data for band energy
        w_n, v_n = np.linalg.eigh(stub_n.ham)
        w_so, v_so = np.linalg.eigh(stub_so.ham)
        stub_n._green_list = {
            "eigenvalue": [w_n], "eigenvector": [v_n], "mu": np.array([0.0])
        }
        stub_so._green_list = {
            "eigenvalue": [w_so], "eigenvector": [v_so], "mu": np.array([0.0])
        }
        stub_n.T = 0
        stub_so.T = 0
        stub_n.ene_cutoff = 1e2
        stub_so.ene_cutoff = 1e2
        stub_n.physics = {"Ene": {}}
        stub_so.physics = {"Ene": {}}
        stub_n.norb = norb_phys
        stub_so.norb = 2 * norb_phys

        stub_n._calc_energy()
        stub_so._calc_energy()

        # Interaction energies should match
        for itype in inter_table:
            if itype in stub_n.physics["Ene"] and itype in stub_so.physics["Ene"]:
                np.testing.assert_allclose(
                    stub_so.physics["Ene"][itype].real,
                    stub_n.physics["Ene"][itype].real,
                    atol=1e-12,
                    err_msg=f"{itype} energy: SO and normal mode should match"
                )

    # ----------------------------------------------------------------
    # Test 8: Block detection with interaction in spin-orbital mode
    # ----------------------------------------------------------------
    def test_block_detection_so_with_coulomb_intra(self):
        """Block detection in spin-orbital mode with CoulombIntra.

        norb_phys=2, nd=4. Diagonal transfer + CoulombIntra.
        CoulombIntra couples (up,a)-(dn,a): in SO, couples 2a with 2a+1.
        Expected: 2 blocks: {0,1} and {2,3} (each orbital pairs up/dn).
        """
        import hwave.solver.uhfk as uhfk_module

        norb_phys = 2
        nd = 2 * norb_phys
        nvol = 1

        # Diagonal transfer in SO ordering
        ham_trans = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for i in range(nd):
            ham_trans[:, i, i] = float(i) + 1.0

        # CoulombIntra
        uab = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        uab[:, 0, 0] = 4.0
        uab[:, 1, 1] = 3.0
        cu_spin = np.zeros((2, 2, 2, 2), dtype=int)
        cu_spin[0, 1, 1, 0] = 1
        cu_spin[1, 0, 0, 1] = 1

        stub = object.__new__(uhfk_module.UHFk)
        stub.norb = nd
        stub.ns = 1
        stub.nd = nd
        stub.nvol = nvol
        stub.ham_trans = ham_trans
        stub.inter_table = {"CoulombIntra": uab}
        stub.spin_table = {"CoulombIntra": cu_spin}
        stub.iflag_fock = True
        stub.sz_free = True
        stub.enable_spin_orbital = True
        stub.norb_phys = norb_phys
        stub.Nconds = [nd]
        stub._detect_blocks()

        blocks = stub.block_info
        self.assertEqual(len(blocks), 2,
                         "CoulombIntra in SO mode: 2 blocks (one per physical orbital)")
        sizes = sorted([len(b) for b in blocks])
        self.assertEqual(sizes, [2, 2])

        # Each block should contain one up (even) and one down (odd) index
        for blk in blocks:
            indices = sorted(blk.tolist())
            self.assertEqual(len(indices), 2)
            self.assertEqual(indices[1] - indices[0], 1,
                             "Block should contain consecutive up/dn pair")

    # ----------------------------------------------------------------
    # Test 9: Block detection SO with Fock connects orbitals
    # ----------------------------------------------------------------
    def test_block_detection_so_fock_bridges_orbitals(self):
        """In SO mode, CoulombInter Fock bridges physical orbitals.

        norb_phys=2, nd=4. CoulombInter J[0,1]!=0 with Fock.
        Fock connects (s,0)-(s,1) -> in SO: 0-2 and 1-3.
        CoulombIntra or Hartree connects (up,a)-(dn,a).
        Combined: all indices connected -> single block.
        """
        import hwave.solver.uhfk as uhfk_module

        norb_phys = 2
        nd = 2 * norb_phys
        nvol = 1

        ham_trans = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for i in range(nd):
            ham_trans[:, i, i] = float(i) + 1.0

        # CoulombInter (Fock bridges orbitals + Hartree bridges spins)
        vab = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        vab[:, 0, 1] = 2.0
        vab[:, 1, 0] = 2.0
        ci_spin = np.zeros((2, 2, 2, 2), dtype=int)
        ci_spin[0, 0, 0, 0] = 1
        ci_spin[1, 1, 1, 1] = 1
        ci_spin[0, 1, 1, 0] = 1
        ci_spin[1, 0, 0, 1] = 1

        stub = object.__new__(uhfk_module.UHFk)
        stub.norb = nd
        stub.ns = 1
        stub.nd = nd
        stub.nvol = nvol
        stub.ham_trans = ham_trans
        stub.inter_table = {"CoulombInter": vab}
        stub.spin_table = {"CoulombInter": ci_spin}
        stub.iflag_fock = True
        stub.sz_free = True
        stub.enable_spin_orbital = True
        stub.norb_phys = norb_phys
        stub.Nconds = [nd]
        stub._detect_blocks()

        blocks = stub.block_info
        self.assertEqual(len(blocks), 1,
                         "CoulombInter Fock + Hartree in SO mode: single block")

    # ----------------------------------------------------------------
    # Test 10: Sz calculation in spin-orbital mode
    # ----------------------------------------------------------------
    def test_sz_calculation_so(self):
        """Verify Sz calculation works correctly in spin-orbital mode."""
        import hwave.solver.uhfk as uhfk_module

        norb_phys = 2
        nd = 2 * norb_phys
        nvol = 1

        stub = object.__new__(uhfk_module.UHFk)
        stub.norb = nd
        stub.ns = 1
        stub.nd = nd
        stub.nvol = nvol
        stub.shape = (1, 1, 1)
        stub.enable_spin_orbital = True
        stub.norb_phys = norb_phys
        stub.threshold = 1e-12
        stub.param_mod = {"Mix": 0.5}
        stub.physics = {"Ene": {"Total": 0.0, "Band": 0.0}, "NCond": 0.0,
                         "Sz": 0.0, "Rest": 1.0}

        # Green function: G[2a, 2a] = n_{a,up}, G[2a+1, 2a+1] = n_{a,dn}
        green_so = np.zeros((nvol, 1, nd, 1, nd), dtype=np.complex128)
        # Orbital 0: n_up=0.7, n_dn=0.3
        green_so[0, 0, 0, 0, 0] = 0.7  # orb0, up
        green_so[0, 0, 1, 0, 1] = 0.3  # orb0, dn
        # Orbital 1: n_up=0.4, n_dn=0.6
        green_so[0, 0, 2, 0, 2] = 0.4  # orb1, up
        green_so[0, 0, 3, 0, 3] = 0.6  # orb1, dn

        stub.Green = green_so
        stub.Green_prev = np.zeros_like(green_so)

        stub._calc_phys()

        # Expected Sz = 0.5 * ((0.7 - 0.3) + (0.4 - 0.6)) = 0.5 * 0.2 = 0.1
        expected_sz = 0.5 * ((0.7 - 0.3) + (0.4 - 0.6))
        np.testing.assert_allclose(
            stub.physics["Sz"], expected_sz, atol=1e-12,
            err_msg="Sz calculation in SO mode should be correct"
        )


    # ================================================================
    # Integrated SCF cycle tests:
    # Run make_ham -> diag -> green -> calc_energy in both normal and
    # spin-orbital modes with spin-diagonal transfer + interactions,
    # verify all results match after multiple iterations.
    # ================================================================

    def _run_scf_cycle(self, stub, n_steps=3, mix=0.5):
        """Run n_steps of the SCF cycle: make_ham -> diag -> green -> calc_energy -> calc_phys.

        Initializes the Green function from the non-interacting ground state
        (transfer Hamiltonian only) to ensure a physical starting point.
        A physical Green function satisfies G_{ab}(r) = G_{ba}(-r)^*, which
        guarantees the mean-field Hamiltonian is Hermitian at every step.

        Parameters
        ----------
        stub : UHFk stub
            Must have ham_trans, inter_table, etc. set up.
        n_steps : int
            Number of SCF iterations.
        mix : float
            Mixing parameter for Green function update.
        """
        stub.param_mod = {"Mix": mix}
        stub.physics = {"Ene": {"Total": 0.0, "Band": 0.0},
                         "NCond": 0.0, "Sz": 0.0, "Rest": 1.0}
        stub._green_list = {}

        # Initialize Green function from non-interacting ground state
        # Set zero Green as placeholder (needed by _green which saves Green_prev)
        stub.Green = np.zeros((stub.nvol, stub.ns, stub.norb, stub.ns, stub.norb),
                              dtype=np.complex128)
        stub.ham = stub.ham_trans.copy()
        stub._diag()
        stub._green()

        for _ in range(n_steps):
            stub._make_ham()
            stub._diag()
            stub._green()
            stub._calc_energy()
            stub._calc_phys()

    def _so_green_to_normal(self, green_so, norb_phys):
        """Convert Green from SO (nvol, 1, nd, 1, nd) to normal (nvol, 2, norb_phys, 2, norb_phys)."""
        nvol = green_so.shape[0]
        nd = 2 * norb_phys
        G_flat = green_so.reshape(nvol, nd, nd)
        G_n = np.zeros((nvol, 2, norb_phys, 2, norb_phys), dtype=np.complex128)
        for s in range(2):
            for t in range(2):
                G_n[:, s, :, t, :] = G_flat[:, s::2, t::2]
        return G_n

    def test_scf_cycle_coulomb_intra(self):
        """Full SCF cycle: CoulombIntra, normal vs spin-orbital mode.

        Sets up spin-diagonal transfer (no spin-orbit coupling in transfer)
        with CoulombIntra interaction. Runs 5 SCF steps in both modes
        and verifies that Green functions, energies, and physical quantities match.
        """
        norb_phys = 2
        nd = 2 * norb_phys
        nvol = 2

        # Spin-diagonal transfer in normal ordering
        ham_t = np.zeros((nvol, nd, nd), dtype=np.complex128)
        # Spin-up block
        ham_t[:, 0, 0] = 1.0
        ham_t[:, 1, 1] = 2.0
        ham_t[:, 0, 1] = 0.3
        ham_t[:, 1, 0] = 0.3
        # Spin-down block (same hopping)
        ham_t[:, 2, 2] = 1.0
        ham_t[:, 3, 3] = 2.0
        ham_t[:, 2, 3] = 0.3
        ham_t[:, 3, 2] = 0.3

        # CoulombIntra: U[0,0]=4.0, U[1,1]=3.0
        uab = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        uab[:, 0, 0] = 4.0
        uab[:, 1, 1] = 3.0
        cu_spin = np.zeros((2, 2, 2, 2), dtype=int)
        cu_spin[0, 1, 1, 0] = 1
        cu_spin[1, 0, 0, 1] = 1

        inter_table = {"CoulombIntra": uab}
        spin_table = {"CoulombIntra": cu_spin}

        stub_n, stub_so = self._make_uhfk_stub_full(
            norb_phys, nvol, ham_t,
            inter_table=inter_table, spin_table=spin_table,
            iflag_fock=True, Nconds=[3]
        )

        # Set additional required attributes for full cycle
        for stub in [stub_n, stub_so]:
            stub.T = 0
            stub.ene_cutoff = 1e2

        # Run 5 SCF steps (Green initialized from non-interacting ground state)
        self._run_scf_cycle(stub_n, n_steps=5)
        self._run_scf_cycle(stub_so, n_steps=5)

        # Compare Green functions (convert SO back to normal ordering)
        green_so_as_normal = self._so_green_to_normal(stub_so.Green, norb_phys)
        np.testing.assert_allclose(
            green_so_as_normal, stub_n.Green, atol=1e-10,
            err_msg="Green functions should match after SCF cycle"
        )

        # Compare energies
        for key in stub_n.physics["Ene"]:
            en = stub_n.physics["Ene"][key]
            eso = stub_so.physics["Ene"][key]
            e_n = en.real if hasattr(en, 'real') else en
            e_so = eso.real if hasattr(eso, 'real') else eso
            np.testing.assert_allclose(
                e_so, e_n, atol=1e-10,
                err_msg=f"Energy '{key}' should match after SCF cycle"
            )

        # Compare physical quantities
        np.testing.assert_allclose(
            stub_so.physics["NCond"], stub_n.physics["NCond"], atol=1e-10,
            err_msg="NCond should match"
        )
        np.testing.assert_allclose(
            stub_so.physics["Sz"], stub_n.physics["Sz"], atol=1e-10,
            err_msg="Sz should match"
        )

    def test_scf_cycle_coulomb_inter_exchange(self):
        """Full SCF cycle: CoulombInter + Exchange, normal vs spin-orbital.

        More complex test with two interaction types and Fock term active.
        Runs 5 SCF steps and compares all outputs.
        """
        norb_phys = 2
        nd = 2 * norb_phys
        nvol = 2

        # Spin-diagonal transfer with different k-point values
        ham_t = np.zeros((nvol, nd, nd), dtype=np.complex128)
        # k=0
        ham_t[0, 0, 0] = 1.0;  ham_t[0, 1, 1] = 2.5
        ham_t[0, 0, 1] = 0.4;  ham_t[0, 1, 0] = 0.4
        ham_t[0, 2, 2] = 1.0;  ham_t[0, 3, 3] = 2.5
        ham_t[0, 2, 3] = 0.4;  ham_t[0, 3, 2] = 0.4
        # k=1 (different hopping)
        ham_t[1, 0, 0] = 1.5;  ham_t[1, 1, 1] = 2.0
        ham_t[1, 0, 1] = -0.2; ham_t[1, 1, 0] = -0.2
        ham_t[1, 2, 2] = 1.5;  ham_t[1, 3, 3] = 2.0
        ham_t[1, 2, 3] = -0.2; ham_t[1, 3, 2] = -0.2

        # CoulombInter: V[0,1]=2.0
        vab = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        vab[:, 0, 1] = 2.0
        vab[:, 1, 0] = 2.0
        ci_spin = np.zeros((2, 2, 2, 2), dtype=int)
        ci_spin[0, 0, 0, 0] = 1
        ci_spin[1, 1, 1, 1] = 1
        ci_spin[0, 1, 1, 0] = 1
        ci_spin[1, 0, 0, 1] = 1

        # Exchange: J[0,1]=-1.5
        jex = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        jex[:, 0, 1] = -1.5
        jex[:, 1, 0] = -1.5
        ex_spin = np.zeros((2, 2, 2, 2), dtype=int)
        ex_spin[0, 1, 0, 1] = 1
        ex_spin[1, 0, 1, 0] = 1

        inter_table = {"CoulombInter": vab, "Exchange": jex}
        spin_table = {"CoulombInter": ci_spin, "Exchange": ex_spin}

        stub_n, stub_so = self._make_uhfk_stub_full(
            norb_phys, nvol, ham_t,
            inter_table=inter_table, spin_table=spin_table,
            iflag_fock=True, Nconds=[3]
        )

        for stub in [stub_n, stub_so]:
            stub.T = 0
            stub.ene_cutoff = 1e2

        self._run_scf_cycle(stub_n, n_steps=5)
        self._run_scf_cycle(stub_so, n_steps=5)

        green_so_as_normal = self._so_green_to_normal(stub_so.Green, norb_phys)
        np.testing.assert_allclose(
            green_so_as_normal, stub_n.Green, atol=1e-10,
            err_msg="Green functions should match (CoulombInter + Exchange)"
        )

        for key in stub_n.physics["Ene"]:
            en = stub_n.physics["Ene"][key]
            eso = stub_so.physics["Ene"][key]
            e_n = en.real if hasattr(en, 'real') else en
            e_so = eso.real if hasattr(eso, 'real') else eso
            np.testing.assert_allclose(
                e_so, e_n, atol=1e-10,
                err_msg=f"Energy '{key}' should match (CoulombInter + Exchange)"
            )

    def test_scf_cycle_all_interactions(self):
        """Full SCF cycle: CoulombIntra + CoulombInter + Hund + Exchange.

        Comprehensive test with all Coulomb-type interactions.
        3 physical orbitals, nvol=2, 5 SCF steps.
        """
        norb_phys = 3
        nd = 2 * norb_phys
        nvol = 2

        # Spin-diagonal transfer
        ham_t = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for k in range(nvol):
            for a in range(norb_phys):
                ham_t[k, a, a] = float(a) + 1.0 + 0.5 * k
                ham_t[k, norb_phys + a, norb_phys + a] = float(a) + 1.0 + 0.5 * k
            ham_t[k, 0, 1] = 0.3; ham_t[k, 1, 0] = 0.3
            ham_t[k, 1, 2] = 0.2; ham_t[k, 2, 1] = 0.2
            ham_t[k, norb_phys+0, norb_phys+1] = 0.3
            ham_t[k, norb_phys+1, norb_phys+0] = 0.3
            ham_t[k, norb_phys+1, norb_phys+2] = 0.2
            ham_t[k, norb_phys+2, norb_phys+1] = 0.2

        # CoulombIntra
        uab = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        uab[:, 0, 0] = 4.0; uab[:, 1, 1] = 3.5; uab[:, 2, 2] = 3.0
        cu_spin = np.zeros((2, 2, 2, 2), dtype=int)
        cu_spin[0, 1, 1, 0] = 1; cu_spin[1, 0, 0, 1] = 1

        # CoulombInter
        vab = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        vab[:, 0, 1] = 1.5; vab[:, 1, 0] = 1.5
        vab[:, 1, 2] = 1.0; vab[:, 2, 1] = 1.0
        ci_spin = np.zeros((2, 2, 2, 2), dtype=int)
        ci_spin[0, 0, 0, 0] = 1; ci_spin[1, 1, 1, 1] = 1
        ci_spin[0, 1, 1, 0] = 1; ci_spin[1, 0, 0, 1] = 1

        # Hund
        jhund = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        jhund[:, 0, 1] = -0.8; jhund[:, 1, 0] = -0.8
        hund_spin = np.zeros((2, 2, 2, 2), dtype=int)
        hund_spin[0, 0, 0, 0] = 1; hund_spin[1, 1, 1, 1] = 1

        # Exchange
        jex = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        jex[:, 0, 1] = -0.5; jex[:, 1, 0] = -0.5
        ex_spin = np.zeros((2, 2, 2, 2), dtype=int)
        ex_spin[0, 1, 0, 1] = 1; ex_spin[1, 0, 1, 0] = 1

        inter_table = {
            "CoulombIntra": uab, "CoulombInter": vab,
            "Hund": jhund, "Exchange": jex,
        }
        spin_table = {
            "CoulombIntra": cu_spin, "CoulombInter": ci_spin,
            "Hund": hund_spin, "Exchange": ex_spin,
        }

        stub_n, stub_so = self._make_uhfk_stub_full(
            norb_phys, nvol, ham_t,
            inter_table=inter_table, spin_table=spin_table,
            iflag_fock=True, Nconds=[4]
        )

        for stub in [stub_n, stub_so]:
            stub.T = 0
            stub.ene_cutoff = 1e2

        self._run_scf_cycle(stub_n, n_steps=5)
        self._run_scf_cycle(stub_so, n_steps=5)

        green_so_as_normal = self._so_green_to_normal(stub_so.Green, norb_phys)
        np.testing.assert_allclose(
            green_so_as_normal, stub_n.Green, atol=1e-10,
            err_msg="Green functions should match (all interactions)"
        )

        for key in stub_n.physics["Ene"]:
            en = stub_n.physics["Ene"][key]
            eso = stub_so.physics["Ene"][key]
            e_n = en.real if hasattr(en, 'real') else en
            e_so = eso.real if hasattr(eso, 'real') else eso
            np.testing.assert_allclose(
                e_so, e_n, atol=1e-10,
                err_msg=f"Energy '{key}' should match (all interactions)"
            )

        np.testing.assert_allclose(
            stub_so.physics["NCond"], stub_n.physics["NCond"], atol=1e-10,
            err_msg="NCond should match (all interactions)"
        )
        np.testing.assert_allclose(
            stub_so.physics["Sz"], stub_n.physics["Sz"], atol=1e-10,
            err_msg="Sz should match (all interactions)"
        )

    def test_scf_cycle_pairhop(self):
        """Full SCF cycle: PairHop interaction, normal vs spin-orbital.

        PairHop has a different contraction structure in _make_ham,
        so it needs a separate test.
        """
        norb_phys = 2
        nd = 2 * norb_phys
        nvol = 2

        ham_t = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for k in range(nvol):
            for a in range(norb_phys):
                ham_t[k, a, a] = float(a) + 1.0 + 0.3 * k
                ham_t[k, norb_phys + a, norb_phys + a] = float(a) + 1.0 + 0.3 * k

        jab = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        jab[:, 0, 1] = 1.0; jab[:, 1, 0] = 1.0
        ph_spin = np.zeros((2, 2, 2, 2), dtype=int)
        ph_spin[0, 1, 1, 0] = 1; ph_spin[1, 0, 0, 1] = 1

        inter_table = {"PairHop": jab}
        spin_table = {"PairHop": ph_spin}

        stub_n, stub_so = self._make_uhfk_stub_full(
            norb_phys, nvol, ham_t,
            inter_table=inter_table, spin_table=spin_table,
            iflag_fock=True, Nconds=[3]
        )

        for stub in [stub_n, stub_so]:
            stub.T = 0
            stub.ene_cutoff = 1e2

        self._run_scf_cycle(stub_n, n_steps=5)
        self._run_scf_cycle(stub_so, n_steps=5)

        green_so_as_normal = self._so_green_to_normal(stub_so.Green, norb_phys)
        np.testing.assert_allclose(
            green_so_as_normal, stub_n.Green, atol=1e-10,
            err_msg="Green functions should match (PairHop)"
        )

        for key in stub_n.physics["Ene"]:
            en = stub_n.physics["Ene"][key]
            eso = stub_so.physics["Ene"][key]
            e_n = en.real if hasattr(en, 'real') else en
            e_so = eso.real if hasattr(eso, 'real') else eso
            np.testing.assert_allclose(
                e_so, e_n, atol=1e-10,
                err_msg=f"Energy '{key}' should match (PairHop)"
            )

    def test_scf_cycle_no_fock(self):
        """Full SCF cycle without Fock term: normal vs spin-orbital.

        Tests that Hartree-only mode also works correctly.
        """
        norb_phys = 2
        nd = 2 * norb_phys
        nvol = 1

        ham_t = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for a in range(norb_phys):
            ham_t[:, a, a] = float(a) + 1.0
            ham_t[:, norb_phys + a, norb_phys + a] = float(a) + 1.0
        ham_t[:, 0, 1] = 0.5; ham_t[:, 1, 0] = 0.5
        ham_t[:, 2, 3] = 0.5; ham_t[:, 3, 2] = 0.5

        # CoulombIntra + CoulombInter (Hartree only)
        uab = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        uab[:, 0, 0] = 4.0; uab[:, 1, 1] = 3.0
        cu_spin = np.zeros((2, 2, 2, 2), dtype=int)
        cu_spin[0, 1, 1, 0] = 1; cu_spin[1, 0, 0, 1] = 1

        vab = np.zeros((nvol, norb_phys, norb_phys), dtype=np.complex128)
        vab[:, 0, 1] = 2.0; vab[:, 1, 0] = 2.0
        ci_spin = np.zeros((2, 2, 2, 2), dtype=int)
        ci_spin[0, 0, 0, 0] = 1; ci_spin[1, 1, 1, 1] = 1
        ci_spin[0, 1, 1, 0] = 1; ci_spin[1, 0, 0, 1] = 1

        inter_table = {"CoulombIntra": uab, "CoulombInter": vab}
        spin_table = {"CoulombIntra": cu_spin, "CoulombInter": ci_spin}

        stub_n, stub_so = self._make_uhfk_stub_full(
            norb_phys, nvol, ham_t,
            inter_table=inter_table, spin_table=spin_table,
            iflag_fock=False, Nconds=[3]  # Fock disabled
        )

        for stub in [stub_n, stub_so]:
            stub.T = 0
            stub.ene_cutoff = 1e2

        self._run_scf_cycle(stub_n, n_steps=5)
        self._run_scf_cycle(stub_so, n_steps=5)

        green_so_as_normal = self._so_green_to_normal(stub_so.Green, norb_phys)
        np.testing.assert_allclose(
            green_so_as_normal, stub_n.Green, atol=1e-10,
            err_msg="Green functions should match (Hartree only, no Fock)"
        )

        for key in stub_n.physics["Ene"]:
            en = stub_n.physics["Ene"][key]
            eso = stub_so.physics["Ene"][key]
            e_n = en.real if hasattr(en, 'real') else en
            e_so = eso.real if hasattr(eso, 'real') else eso
            np.testing.assert_allclose(
                e_so, e_n, atol=1e-10,
                err_msg=f"Energy '{key}' should match (Hartree only)"
            )


if __name__ == '__main__':
    unittest.main()
