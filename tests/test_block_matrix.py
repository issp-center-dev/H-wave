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


if __name__ == '__main__':
    unittest.main()
