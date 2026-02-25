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


class TestUHFkDetectBlocks(unittest.TestCase):
    """Test _detect_blocks in UHFk for various transfer/interaction patterns.

    Index convention: nd = norb * ns.
    For ns=2:  indices 0..norb-1 are spin-up, norb..2*norb-1 are spin-down.
    The (s, a) pair maps to index s*norb + a.
    """

    def _make_uhfk_stub(self, norb, ns, nvol, ham_trans,
                        inter_table=None, spin_table=None,
                        iflag_fock=True, sz_free=True, Nconds=None):
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
                                    sz_free=True, Nconds=[2])
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


if __name__ == '__main__':
    unittest.main()
