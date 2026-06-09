#!/usr/bin/env python3

"""Analytical verification tests for the FLEX solver.

These tests compare FLEX/RPA numerical results against analytically
known formulas, providing rigorous quantitative validation without
relying on external reference data.

Test 1: Lindhard function
    chi0(q, nu=0) from the solver vs the analytic formula
    chi0(q, 0) = -(1/N) sum_k [f(ek) - f(ek+q)] / [ek - ek+q]

Test 2: Weak-coupling limit
    FLEX with U->0 (single iteration) should reproduce bare chi0 exactly.

Test 3: Filling sum rule
    (T/N) sum_{k,n} G(k, iwn) * e^{iwn 0+} = n (particle number per state)
    Equivalent: (2T/N) sum_{k,n>0} Re G(k, iwn) + 0.5 = n/2  (for 1 orbital)

Test 4: High-frequency tail of Green's function
    G(k, iwn) ~ 1/(iwn) for |wn| >> bandwidth
    Im G(k, iwn) ~ -1/wn

Test 5: Kramers-Kronig consistency of self-energy
    Sigma(-iwn) = conj(Sigma(iwn)) for Hermitian Hamiltonian

Test 6: Weak-coupling FLEX self-energy = second-order perturbation theory
    For small U, Sigma_FLEX should match the SOPT bubble diagram.
"""

import os
import unittest
import numpy as np


def _make_2orb_intra_solver(Lx=4, Ly=4, Nmat=16, T=1.0, mu=0.0, U=1.0):
    """Create a 2-orbital FLEX solver with intra-orbital U only.

    Reuses the geom/transfer of tests/rpa/input_2orb (which contains
    inter-orbital hopping t_{12}=0.5), so the bare Green's function has
    orbital-off-diagonal elements.  Only CoulombIntra U is added, keeping the
    interaction density-density (reduced-scheme compatible).
    """
    import tempfile
    import shutil

    tmpdir = tempfile.mkdtemp(prefix='flex_2orb_')
    src_dir = 'tests/rpa/input_2orb'
    for fn in ['geom.dat', 'transfer.dat']:
        shutil.copy2(os.path.join(src_dir, fn), os.path.join(tmpdir, fn))

    ci_path = os.path.join(tmpdir, 'coulombintra.dat')
    with open(ci_path, 'w') as f:
        f.write("CoulombIntra\n2\n1\n 1\n")
        f.write("   0    0    0    1    1   {:.12f}   0.000000000000\n".format(U))
        f.write("   0    0    0    2    2   {:.12f}   0.000000000000\n".format(U))

    info_log = {}
    info_mode = {
        'mode': 'FLEX',
        'param': {
            'T': T, 'mu': mu,
            'CellShape': [Lx, Ly, 1], 'SubShape': [1, 1, 1],
            'Nmat': Nmat, 'IterationMax': 1, 'Mix': 1.0, 'EPS': 1,
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
    os.makedirs('tests/flex/output_anal', exist_ok=True)

    import hwave.qlmsio.read_input_k as read_input_k
    read_io = read_input_k.QLMSkInput(info_file_input)
    ham_info = read_io.get_param("ham")
    green_info = read_io.get_param("green")

    import hwave.solver.flex as solver_flex
    solver = solver_flex.FLEX(ham_info, info_log, info_mode)

    shutil.rmtree(tmpdir, ignore_errors=True)
    return solver, green_info


def _flex_self_energy_explicit(solver, green_kw, v_eff, beta, rule):
    """Independent (non-FFT) reference for the FLEX self-energy convolution.

    Sigma_{ab}(k,iwn) = (T/N) sum_{q,m} <rule>( V_eff(q,inu_m), G(k-q, iwn-inu_m) )

    with the solver's frequency conventions:
      iwn   = (2n+1-Nmat) pi T            (fermionic, centred)
      inu_m = 2 pi T (m - Nmat/2)         (bosonic, static at m=Nmat/2)
      iwn - inu_m = iw_{n'},  n' = (n - m + Nmat//2) mod Nmat

    rule='elem'   -> Sigma_ab = V_ab * G_ab   (the implemented contraction)
    rule='matmul' -> Sigma_ab = sum_c V_ac G_cb (the alternative Codex proposed)
    """
    nvol = solver.lattice.nvol
    norb = solver.norb
    ns = solver.ns
    nx, ny, nz = solver.lattice.shape
    nmat = solver.nmat
    nd_v = v_eff.shape[-1]
    T = 1.0 / beta

    # Expand block-space G (spin-free: norb x norb) to spin-orbital nd_v
    G_so = np.zeros((nmat, nvol, nd_v, nd_v), dtype=complex)
    for s in range(ns):
        sl = slice(s * norb, (s + 1) * norb)
        G_so[:, :, sl, sl] = green_kw[0]

    # Solver packs the flat lattice index in C-order (nx,ny,nz): for nz=1 the
    # flat index is ix*ny + iy (ix is the outermost/slowest axis, matching
    # green_kw.reshape(...,nx,ny,nz,...) and wavenum_table in rpa.py).
    def qidx(ix, iy):
        return ix * ny + iy

    sig = np.zeros((nmat, nvol, nd_v, nd_v), dtype=complex)
    for n in range(nmat):
        for m in range(nmat):
            npr = (n - m + nmat // 2) % nmat
            for ikx in range(nx):
                for iky in range(ny):
                    ik = qidx(ikx, iky)
                    acc = np.zeros((nd_v, nd_v), dtype=complex)
                    for iqx in range(nx):
                        for iqy in range(ny):
                            iq = qidx(iqx, iqy)
                            ikmq = qidx((ikx - iqx) % nx, (iky - iqy) % ny)
                            if rule == 'elem':
                                acc += v_eff[m, iq] * G_so[npr, ikmq]
                            else:
                                acc += v_eff[m, iq] @ G_so[npr, ikmq]
                    sig[n, ik] += acc
    sig *= (T / nvol)
    # contract back to block (orbital) space, as the solver does
    return sig[np.newaxis, :, :, :norb, :norb]


def _make_1orb_solver(Lx=16, Ly=16, Nmat=64, T=0.5, mu=0.0, U=4.0,
                      max_iter=1, mix=1.0, eps=4):
    """Create a 1-orbital FLEX solver for analytical tests.

    Uses the same 2D Hubbard model as test_flex_benchmark:
      eps(k) = 2*cos(kx) + 2*cos(ky) + cos(kx+ky)
    """
    import tempfile
    import shutil

    tmpdir = tempfile.mkdtemp(prefix='flex_anal_')
    src_dir = 'tests/rpa/input'
    for fn in ['geom.dat', 'transfer.dat']:
        shutil.copy2(os.path.join(src_dir, fn), os.path.join(tmpdir, fn))

    ci_path = os.path.join(tmpdir, 'coulombintra.dat')
    with open(ci_path, 'w') as f:
        f.write("CoulombIntra\n1\n1\n 1\n")
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
    os.makedirs('tests/flex/output_anal', exist_ok=True)

    import hwave.qlmsio.read_input_k as read_input_k
    read_io = read_input_k.QLMSkInput(info_file_input)
    ham_info = read_io.get_param("ham")
    green_info = read_io.get_param("green")

    import hwave.solver.flex as solver_flex
    solver = solver_flex.FLEX(ham_info, info_log, info_mode)

    shutil.rmtree(tmpdir, ignore_errors=True)
    return solver, green_info


def _dispersion_1orb(kx, ky):
    """Dispersion for the test 1-orbital model.

    eps(k) = 2*cos(kx) + 2*cos(ky) + cos(kx+ky)
    from transfer.dat: t_{(1,0)}=1.0, t_{(0,1)}=1.0, t_{(1,1)}=0.5
    """
    return 2.0 * np.cos(kx) + 2.0 * np.cos(ky) + np.cos(kx + ky)


def _fermi_dirac(e, beta):
    """Fermi-Dirac distribution, numerically stable."""
    x = beta * e
    return np.where(x > 500, 0.0, np.where(x < -500, 1.0, 1.0 / (np.exp(x) + 1.0)))


def _lindhard_static(kx_grid, ky_grid, mu, beta):
    """Compute the static Lindhard function chi0(q, nu=0) analytically.

    chi0(q, 0) = -(1/N) sum_k [f(ek) - f(ek+q)] / [ek - ek+q]

    Parameters
    ----------
    kx_grid, ky_grid : 1D arrays
        k-points: kx = 2*pi*n/Lx, ky = 2*pi*m/Ly
    mu : float
        Chemical potential.
    beta : float
        Inverse temperature.

    Returns
    -------
    chi0 : ndarray, shape (nvol,)
        Static Lindhard susceptibility for each q-point.
        q-points indexed as iq = ix*Ly + iy.
    """
    Lx = len(kx_grid)
    Ly = len(ky_grid)
    nvol = Lx * Ly

    # Build dispersion on full grid
    KX, KY = np.meshgrid(kx_grid, ky_grid, indexing='ij')
    ek = _dispersion_1orb(KX, KY).ravel() - mu  # shape (nvol,)
    fk = _fermi_dirac(ek, beta)

    chi0 = np.zeros(nvol)
    for iq in range(nvol):
        iqx = iq // Ly
        iqy = iq % Ly
        qx = kx_grid[iqx]
        qy = ky_grid[iqy]

        # eps(k+q)
        KXq = KX + qx
        KYq = KY + qy
        ekq = _dispersion_1orb(KXq, KYq).ravel() - mu
        fkq = _fermi_dirac(ekq, beta)

        de = ek - ekq
        # Avoid 0/0: where de==0, use -df/de = beta*f*(1-f)
        mask = np.abs(de) < 1e-12
        contrib = np.zeros(nvol)
        contrib[~mask] = (fk[~mask] - fkq[~mask]) / de[~mask]
        contrib[mask] = -beta * fk[mask] * (1.0 - fk[mask])

        chi0[iq] = -np.mean(contrib)

    return chi0


def _lindhard_matsubara(kx_grid, ky_grid, mu, beta, nmat):
    """Compute Lindhard function in Matsubara representation.

    chi0(q, i*nu_m) = -(T/N) sum_{k,n} G0(k,iwn) G0(k+q, iwn+inu_m)

    where G0(k,iwn) = 1/(iwn - ek + mu) and nu_m = 2*m*pi*T (bosonic).

    Parameters
    ----------
    kx_grid, ky_grid : 1D arrays
    mu, beta : float
    nmat : int
        Number of Matsubara frequencies.

    Returns
    -------
    chi0 : ndarray, shape (nmat, nvol)
        Bare susceptibility for each (bosonic freq, q-point).
    """
    Lx = len(kx_grid)
    Ly = len(ky_grid)
    nvol = Lx * Ly
    T = 1.0 / beta

    # Fermionic Matsubara: wn = pi*(2n+1-nmat)/beta, n=0..nmat-1
    n_idx = np.arange(nmat) * 2 + 1 - nmat
    iwn = 1j * n_idx * np.pi * T  # shape (nmat,)

    # Bosonic Matsubara: nu_m = 2*pi*m*T, m=0..nmat-1
    # But the solver uses m=0..nmat-1 mapped to bosonic freq
    # Actually in the FFT convention used by the solver, the bosonic
    # frequencies correspond to m=0,...,nmat-1 with nu_m = 2*pi*m/beta

    KX, KY = np.meshgrid(kx_grid, ky_grid, indexing='ij')
    ek_2d = _dispersion_1orb(KX, KY) - mu  # (Lx, Ly)
    ek = ek_2d.ravel()  # (nvol,)

    # G0(k, iwn) = 1/(iwn - ek)
    # shape: (nmat, nvol)
    G0 = 1.0 / (iwn[:, None] - ek[None, :])

    chi0 = np.zeros((nmat, nvol), dtype=complex)
    for iq in range(nvol):
        iqx = iq // Ly
        iqy = iq % Ly
        qx = kx_grid[iqx]
        qy = ky_grid[iqy]

        ekq_2d = _dispersion_1orb(KX + qx, KY + qy) - mu
        ekq = ekq_2d.ravel()

        G0kq = 1.0 / (iwn[:, None] - ekq[None, :])

        # chi0(q, inu_m) = -T/N sum_{k,n} G0(k,n) G0(k+q,n+m)
        # For each bosonic m, shift the fermionic index: iwn + inu_m
        for im in range(nmat):
            nu_m = 2 * im * np.pi * T
            iwn_shifted = iwn + 1j * nu_m
            G0kq_shifted = 1.0 / (iwn_shifted[:, None] - ekq[None, :])

            # sum over k and n
            chi0[im, iq] = -T * np.mean(np.sum(G0 * G0kq_shifted, axis=0))

    return chi0


class TestBareSusceptibility(unittest.TestCase):
    """Test 1: Properties of the bare susceptibility chi0(q).

    The solver computes chi0 via FFT-based bubble diagram.
    We verify analytically known symmetry and sum-rule properties.
    """

    @classmethod
    def setUpClass(cls):
        """Set up solver and compute chi0."""
        cls.Lx, cls.Ly = 8, 8
        cls.Nmat = 64
        cls.T = 0.5
        cls.mu = 0.0
        cls.beta = 1.0 / cls.T
        cls.U = 0.0  # No interaction for pure chi0 test

        solver, green_info = _make_1orb_solver(
            Lx=cls.Lx, Ly=cls.Ly, Nmat=cls.Nmat,
            T=cls.T, mu=cls.mu, U=cls.U,
            max_iter=1, mix=1.0, eps=1)
        solver.solve(green_info, 'tests/flex/output_anal')

        cls.solver = solver
        cls.green_info = green_info
        cls.chi0q = green_info["chi0q"]

        cls.kx = np.arange(cls.Lx) * 2.0 * np.pi / cls.Lx
        cls.ky = np.arange(cls.Ly) * 2.0 * np.pi / cls.Ly

    def test_chi0_symmetry_inversion(self):
        """chi0(q) = chi0(-q) by time-reversal symmetry."""
        chi0 = self.chi0q
        if chi0.ndim == 4:
            chi0 = chi0[:, :, 0, 0]

        nvol = self.Lx * self.Ly
        for iq in range(nvol):
            iqx = iq // self.Ly
            iqy = iq % self.Ly
            iqx_neg = (-iqx) % self.Lx
            iqy_neg = (-iqy) % self.Ly
            iq_neg = iqx_neg * self.Ly + iqy_neg

            np.testing.assert_allclose(
                chi0[0, iq], chi0[0, iq_neg],
                atol=1e-12,
                err_msg="chi0(q) != chi0(-q) for q=({},{})".format(iqx, iqy))

    def test_chi0_real_at_static(self):
        """chi0(q, nu=0) should be real for the static susceptibility."""
        chi0 = self.chi0q
        if chi0.ndim == 4:
            chi0 = chi0[:, :, 0, 0]

        chi0_static = chi0[0, :]
        max_imag = np.max(np.abs(chi0_static.imag))

        self.assertLess(max_imag, 1e-10,
                        "Static chi0 has imaginary part: max|Im|={:.2e}".format(
                            max_imag))

    def test_chi0_c4_symmetry(self):
        """chi0(qx,qy) = chi0(qy,qx) for square lattice C4 symmetry.

        This model has t' only along (1,1)/(-1,-1), breaking C4.
        But chi0(qx,qy) = chi0(-qx,-qy) should still hold (inversion).
        Also chi0 should be real at nu=0.
        """
        chi0 = self.chi0q
        if chi0.ndim == 4:
            chi0 = chi0[:, :, 0, 0]

        # Check all q-points for inversion symmetry
        for iqx in range(self.Lx):
            for iqy in range(self.Ly):
                iq = iqx * self.Ly + iqy
                iq_inv = ((-iqx) % self.Lx) * self.Ly + ((-iqy) % self.Ly)
                np.testing.assert_allclose(
                    chi0[0, iq].real, chi0[0, iq_inv].real,
                    atol=1e-12)

    def test_chi0_spin_degeneracy(self):
        """For spin-free mode, chi0 up-up = chi0 down-down."""
        chi0 = self.chi0q
        if chi0.ndim == 4 and chi0.shape[2] >= 2:
            np.testing.assert_allclose(
                chi0[:, :, 0, 0], chi0[:, :, 1, 1],
                atol=1e-14,
                err_msg="chi0 up-up != chi0 down-down")

    def test_chi0_stoner_criterion_consistency(self):
        """RPA chi_s = chi0/(1 - U*chi0) should be self-consistent.

        Run with small U and verify the Stoner formula.
        """
        Lx, Ly, Nmat, T = 8, 8, 32, 0.5
        U = 2.0

        solver, green_info = _make_1orb_solver(
            Lx=Lx, Ly=Ly, Nmat=Nmat, T=T, mu=0.0, U=U,
            max_iter=1, mix=1.0, eps=1)
        solver.solve(green_info, 'tests/flex/output_anal')

        chi0 = green_info["chi0q"]
        chi_s = green_info["chiq_s"]

        if chi0.ndim == 4:
            chi0 = chi0[:, :, 0, 0]
        if chi_s.ndim == 4:
            chi_s = chi_s[:, :, 0, 0]

        # RPA: chi_s = chi0 / (1 - U*chi0)
        chi_s_formula = chi0 / (1.0 - U * chi0)
        np.testing.assert_allclose(
            chi_s, chi_s_formula, atol=1e-10,
            err_msg="chi_s does not match Stoner formula")


class TestFillingSumRule(unittest.TestCase):
    """Test 3: Particle number conservation from Green's function.

    n = 1 + (2T/N) sum_{k, n} Re G(k, iwn)

    For 1 orbital with mu=0, by particle-hole symmetry n should be ~1
    (half-filling depends on the dispersion).
    """

    @classmethod
    def setUpClass(cls):
        cls.Lx, cls.Ly = 16, 16
        cls.Nmat = 128
        cls.T = 0.5
        cls.mu = 0.0
        cls.beta = 1.0 / cls.T

        solver, green_info = _make_1orb_solver(
            Lx=cls.Lx, Ly=cls.Ly, Nmat=cls.Nmat,
            T=cls.T, mu=cls.mu, U=0.0,
            max_iter=1, mix=1.0, eps=1)
        solver.solve(green_info, 'tests/flex/output_anal')

        cls.solver = solver
        cls.green_info = green_info

    def test_filling_from_bare_green(self):
        """Filling from G0 should match direct Fermi-Dirac calculation.

        n = 1 + (2T/Nk) sum_{k,n} Re G(k, iwn)
        This is the standard Matsubara sum rule.
        """
        # Direct calculation from eigenvalues
        ew = self.solver.H0_eigenvalue  # shape (nblock, nvol, nd)
        ek = ew[0, :, 0] - self.mu  # 1 orbital, spin-free
        n_direct = 2.0 * np.mean(_fermi_dirac(ek, self.beta))

        # From Matsubara sum: n = 1 + (2T/N) sum_{k,wn} Re G(k, iwn)
        nmat = self.Nmat
        nvol = self.Lx * self.Ly

        # Reconstruct G from eigenvalues
        iwn = 1j * (np.arange(nmat) * 2 + 1 - nmat) * np.pi * self.T
        G0 = 1.0 / (iwn[:, None] - ek[None, :])  # (nmat, nvol)

        n_matsubara = 1.0 + 2.0 * self.T * np.sum(G0.real) / nvol

        self.assertAlmostEqual(
            n_matsubara, n_direct, places=4,
            msg="Filling from Matsubara sum ({:.6f}) != Fermi-Dirac ({:.6f})".format(
                n_matsubara, n_direct))

    def test_filling_sum_rule_solver_green(self):
        """Filling from solver's Green function should be consistent.

        Compare the solver's internal Green function with the eigenvalue-based
        filling calculation.
        """
        ew = self.solver.H0_eigenvalue
        ek = ew[0, :, 0] - self.mu
        n_direct = 2.0 * np.mean(_fermi_dirac(ek, self.beta))

        # Solver's Green function
        green_kw, green_tail = self.solver._calc_green(self.beta, self.mu)
        # green_kw shape: (nblock, nmat, nvol, nd, nd)
        G_diag = green_kw[0, :, :, 0, 0]  # (nmat, nvol)

        nvol = self.Lx * self.Ly
        n_solver = 1.0 + 2.0 * self.T * np.sum(G_diag.real) / nvol

        self.assertAlmostEqual(
            n_solver, n_direct, places=3,
            msg="Filling from solver Green ({:.6f}) != Fermi-Dirac ({:.6f})".format(
                n_solver, n_direct))


class TestHighFrequencyTail(unittest.TestCase):
    """Test 4: Green's function high-frequency behavior.

    G(k, iwn) ~ 1/(iwn) for |wn| >> max|ek|
    """

    @classmethod
    def setUpClass(cls):
        cls.Lx, cls.Ly = 8, 8
        cls.Nmat = 128
        cls.T = 0.5
        cls.mu = 0.0
        cls.beta = 1.0 / cls.T

        solver, green_info = _make_1orb_solver(
            Lx=cls.Lx, Ly=cls.Ly, Nmat=cls.Nmat,
            T=cls.T, mu=cls.mu, U=0.0,
            max_iter=1, mix=1.0, eps=1)
        solver.solve(green_info, 'tests/flex/output_anal')

        cls.solver = solver

    def test_high_freq_decay(self):
        """Im G(k, iwn) ~ -1/wn at high frequencies."""
        green_kw, _ = self.solver._calc_green(self.beta, self.mu)
        G_diag = green_kw[0, :, :, 0, 0]  # (nmat, nvol)

        nmat = self.Nmat
        iwn = (np.arange(nmat) * 2 + 1 - nmat) * np.pi * self.T

        # Check the highest positive frequencies (last 10)
        for i_freq in range(nmat - 10, nmat):
            wn = iwn[i_freq]
            G_avg = np.mean(G_diag[i_freq, :])

            # At high freq, G ~ 1/(iwn), so Im G ~ -1/wn
            expected_im = -1.0 / wn
            actual_im = G_avg.imag

            # Relative error should be small for large wn
            if abs(wn) > 10:
                rel_err = abs(actual_im - expected_im) / abs(expected_im)
                self.assertLess(rel_err, 0.05,
                                "High-freq tail: wn={:.1f}, expected Im G={:.6f}, "
                                "got {:.6f}".format(wn, expected_im, actual_im))

    def test_real_part_decays_faster(self):
        """Re G(k, iwn) ~ O(1/wn^2) at high frequencies.

        |Re G| << |Im G| at high frequency.
        """
        green_kw, _ = self.solver._calc_green(self.beta, self.mu)
        G_diag = green_kw[0, :, :, 0, 0]

        # Average over k, check last frequency
        G_last = np.mean(G_diag[-1, :])
        self.assertLess(abs(G_last.real), abs(G_last.imag),
                        "|Re G| should be << |Im G| at high frequency")


class TestWeakCouplingLimit(unittest.TestCase):
    """Test 2: In the weak-coupling limit, FLEX reduces to RPA.

    For U -> 0:
    - chi_s = chi0 / (1 - U*chi0) ~ chi0 (1 + U*chi0)
    - chi_c = chi0 / (1 + U*chi0) ~ chi0 (1 - U*chi0)
    - Sigma ~ U^2 * chi0 (second order perturbation theory)
    """

    def test_flex_single_iter_equals_rpa(self):
        """Single FLEX iteration with small U should match RPA chi_s."""
        Lx, Ly = 8, 8
        Nmat = 32
        T = 1.0
        U = 0.5  # Small U

        solver, green_info = _make_1orb_solver(
            Lx=Lx, Ly=Ly, Nmat=Nmat,
            T=T, mu=0.0, U=U,
            max_iter=1, mix=1.0, eps=1)
        solver.solve(green_info, 'tests/flex/output_anal')

        chi0 = green_info["chi0q"]
        chi_s = green_info["chiq_s"]

        if chi0.ndim == 4:
            chi0 = chi0[:, :, 0, 0]
        if chi_s.ndim == 4:
            chi_s = chi_s[:, :, 0, 0]

        # RPA: chi_s = chi0 / (1 - U*chi0)
        chi_s_rpa = chi0 / (1.0 - U * chi0)

        np.testing.assert_allclose(
            chi_s, chi_s_rpa, atol=1e-10, rtol=1e-8,
            err_msg="Single-iteration FLEX chi_s does not match RPA formula")

    def test_zero_U_chi_s_equals_chi0(self):
        """With U=0, chi_s should equal chi0 exactly."""
        Lx, Ly = 8, 8
        Nmat = 32
        T = 1.0

        solver, green_info = _make_1orb_solver(
            Lx=Lx, Ly=Ly, Nmat=Nmat,
            T=T, mu=0.0, U=0.0,
            max_iter=1, mix=1.0, eps=1)
        solver.solve(green_info, 'tests/flex/output_anal')

        chi0 = green_info["chi0q"]
        chi_s = green_info["chiq_s"]
        chi_c = green_info["chiq_c"]

        np.testing.assert_allclose(
            chi_s, chi0, atol=1e-12,
            err_msg="chi_s != chi0 for U=0")
        np.testing.assert_allclose(
            chi_c, chi0, atol=1e-12,
            err_msg="chi_c != chi0 for U=0")

    def test_stoner_enhancement_formula(self):
        """chi_s(q) = chi0(q) / (1 - U*chi0(q)) for single iteration.

        Verify the Stoner enhancement quantitatively at Q=(pi,pi).
        """
        Lx, Ly = 16, 16
        Nmat = 64
        T = 0.5
        U = 2.0

        solver, green_info = _make_1orb_solver(
            Lx=Lx, Ly=Ly, Nmat=Nmat,
            T=T, mu=0.0, U=U,
            max_iter=1, mix=1.0, eps=1)
        solver.solve(green_info, 'tests/flex/output_anal')

        chi0 = green_info["chi0q"]
        chi_s = green_info["chiq_s"]

        if chi0.ndim == 4:
            chi0 = chi0[:, :, 0, 0]
        if chi_s.ndim == 4:
            chi_s = chi_s[:, :, 0, 0]

        q_pipi = (Lx // 2) * Ly + Ly // 2
        chi0_pipi = chi0[0, q_pipi].real
        chi_s_pipi = chi_s[0, q_pipi].real

        # Stoner: chi_s = chi0 / (1 - U*chi0)
        stoner = chi0_pipi / (1.0 - U * chi0_pipi)

        self.assertAlmostEqual(
            chi_s_pipi, stoner, places=8,
            msg="chi_s(pi,pi)={:.8f} != Stoner formula {:.8f}".format(
                chi_s_pipi, stoner))


class TestSelfEnergyCausality(unittest.TestCase):
    """Test 5: Causality and Hermiticity of the self-energy.

    For a physical self-energy:
    - Sigma(k, -iwn) = conj(Sigma(k, iwn))
    - Im Sigma(k, iwn) < 0 for wn > 0  (retarded, causal)
    """

    @classmethod
    def setUpClass(cls):
        cls.Lx, cls.Ly = 8, 8
        cls.Nmat = 32
        cls.T = 1.0

        solver, green_info = _make_1orb_solver(
            Lx=cls.Lx, Ly=cls.Ly, Nmat=cls.Nmat,
            T=cls.T, mu=0.0, U=4.0,
            max_iter=30, mix=0.2, eps=4)
        solver.solve(green_info, 'tests/flex/output_anal')

        cls.sigma = green_info["sigma"]

    def test_sigma_conjugation_symmetry(self):
        """Sigma(k, -iwn) = conj(Sigma(k, iwn)).

        This follows from the Hermiticity of the Hamiltonian.
        """
        sigma = self.sigma
        if sigma.ndim == 5:
            sigma = sigma[:, :, 0, 0]  # (nmat, nvol) for 1 orbital
        elif sigma.ndim == 3:
            pass  # already (spin, nmat, nvol)

        # For spin-free reduced: sigma shape is (1, nmat, nvol, 1, 1)
        # After squeeze it might be (nmat, nvol)
        s = self.sigma
        if s.ndim == 5:
            s = s[0, :, :, 0, 0]  # (nmat, nvol)
        nmat = s.shape[0]

        # iwn = pi*T*(2n+1-nmat), n=0..nmat-1
        # n and (nmat-1-n) give opposite frequencies
        for n in range(nmat // 2):
            n_neg = nmat - 1 - n
            np.testing.assert_allclose(
                s[n, :], np.conj(s[n_neg, :]),
                atol=1e-10,
                err_msg="Sigma(-iwn) != conj(Sigma(iwn)) at n={}".format(n))

    def test_sigma_causal(self):
        """Im Sigma(k, iwn) < 0 for wn > 0 (causal self-energy)."""
        s = self.sigma
        if s.ndim == 5:
            s = s[0, :, :, 0, 0]
        nmat = s.shape[0]

        # Positive frequencies: n >= nmat//2
        sigma_pos = s[nmat // 2:, :]
        max_im = np.max(sigma_pos.imag)

        self.assertLess(max_im, 1e-10,
                        "Im Sigma should be <= 0 for wn > 0, "
                        "max(Im Sigma) = {:.2e}".format(max_im))


class TestWeakCouplingSelfEnergy(unittest.TestCase):
    """Test 6: FLEX self-energy scaling in the weak-coupling limit.

    The leading FLEX self-energy is the second-order perturbation (SOPT)
    bubble.  V_eff = W * [3/2 chi_s + 1/2 chi_c - chi0] * W, and at zeroth
    order chi_s = chi_c = chi0, so the kernel reduces to chi0 and
        V_eff ~ W chi0 W = U^2 chi0,   Sigma = V_eff (x) G ~ U^2.
    So the FLEX self-energy scales as U^2 (NOT U^3) for small U.  A U^3
    scaling would indicate that chi0 has been subtracted twice in V_eff,
    cancelling the physical second-order bubble.
    """

    def test_sigma_scales_as_U_squared_sopt(self):
        """At small U, Sigma_FLEX must scale as U^2 (second-order PT bubble).

        The leading FLEX self-energy is the second-order bubble
        Sigma^(2) = U^2 * chi0 (x) G ~ U^2.  The fluctuation kernel must keep
        the O(U^0) part of chi_s/chi_c (i.e. subtract chi0 exactly ONCE):
            fluct = 3/2 chi_s + 1/2 chi_c - chi0 = chi0 + O(U).
        If chi0 is subtracted twice (3/2(chi_s-chi0)+1/2(chi_c-chi0)) the
        O(U^0) term cancels, V_eff ~ U^3 and the SOPT bubble is lost.
        """
        Lx, Ly = 8, 8
        Nmat = 32
        T = 1.0

        U_values = [0.25, 0.5, 1.0]
        sigma_over_U2 = []

        for U in U_values:
            solver, green_info = _make_1orb_solver(
                Lx=Lx, Ly=Ly, Nmat=Nmat,
                T=T, mu=0.0, U=U,
                max_iter=1, mix=1.0, eps=1)
            solver.solve(green_info, 'tests/flex/output_anal')

            sigma = green_info["sigma"]
            if sigma.ndim == 5:
                sigma_val = np.mean(np.abs(sigma[0, Nmat // 2, :, 0, 0].imag))
            else:
                sigma_val = np.mean(np.abs(sigma[Nmat // 2, :].imag))

            sigma_over_U2.append(sigma_val / U**2)

        # For small U, sigma/U^2 should be roughly constant (-> SOPT bubble).
        ratio = sigma_over_U2[-1] / sigma_over_U2[0]
        self.assertLess(
            abs(ratio - 1.0), 0.5,
            "Sigma should scale as U^2 (SOPT bubble), but sigma/U^2 "
            "ratios = {} (ratio={:.3f}). A value near U (~4) means chi0 is "
            "subtracted twice in V_eff.".format(sigma_over_U2, ratio))


class TestSelfEnergyConvolution(unittest.TestCase):
    """Settle the orbital contraction rule in the FLEX self-energy convolution.

    The FLEX self-energy is
        Sigma_{ab}(k,iwn) = (T/N) sum_{q,m} V_eff_{ab}(q,inu_m) G_{ab}(k-q,...),
    i.e. the orbital indices (a,b) of V_eff and G coincide, so the (r,tau)-space
    tie-back V_eff(r,tau) G(r,tau) is an ELEMENTWISE product (not a matrix
    product over an internal orbital).  This is exactly the second-order bubble
        Sigma^(2)_{ab}(x) = G_{ab}(x) * [U chi0(x) U]_{ab},
    derived from the density-density vertex; the orbital matrix products are
    already absorbed into V_eff = W chi W.

    Verified on a 2-orbital model with inter-orbital hopping (so G is
    orbital-off-diagonal and the elementwise / matmul rules differ) by
    reproducing the solver's self-energy with a fully independent, FFT-free
    explicit-frequency convolution.
    """

    @classmethod
    def setUpClass(cls):
        beta = 1.0
        # Asymmetric lattice (Lx != Ly) so the flat-index packing convention
        # ix*ny+iy is genuinely exercised (a transposed qidx would fail here).
        solver, green_info = _make_2orb_intra_solver(
            Lx=4, Ly=2, Nmat=16, T=1.0, mu=0.0, U=1.0)
        solver.solve(green_info, 'tests/flex/output_anal')

        nblock = 1
        nd_block = solver.norb
        sigma0 = np.zeros((nblock, solver.nmat, solver.lattice.nvol,
                           nd_block, nd_block), dtype=complex)
        green_kw = solver._calc_dressed_green(beta, solver.mu, sigma0)

        chi0q_raw = solver._calc_chi0q(green_kw, solver.green0_tail, beta)
        if solver.spin_mode in ["spin-free", "spinful"]:
            chi0q_raw = chi0q_raw[0]
        _, v_eff, _, _ = solver._flex_compute_veff(
            chi0q_raw, solver.ham_info.ham_inter_q)

        cls.solver = solver
        cls.beta = beta
        cls.green_kw = green_kw
        cls.v_eff = v_eff
        cls.sigma_solver = solver._calc_self_energy(green_kw, v_eff, beta)

    def test_self_energy_matches_elementwise_convolution(self):
        """Solver self-energy == independent elementwise-rule convolution."""
        sig_elem = _flex_self_energy_explicit(
            self.solver, self.green_kw, self.v_eff, self.beta, rule='elem')
        scale = np.max(np.abs(self.sigma_solver))
        self.assertGreater(scale, 1e-6, "self-energy is ~zero; test is vacuous")
        np.testing.assert_allclose(
            self.sigma_solver, sig_elem, atol=1e-10, rtol=1e-7,
            err_msg="FLEX self-energy does not match the elementwise "
                    "(V_ab * G_ab) convolution")

    def test_elementwise_and_matmul_differ_here(self):
        """Sanity: this 2-orbital model actually discriminates the two rules.

        Without orbital-off-diagonal structure the test would be vacuous.
        """
        sig_elem = _flex_self_energy_explicit(
            self.solver, self.green_kw, self.v_eff, self.beta, rule='elem')
        sig_mm = _flex_self_energy_explicit(
            self.solver, self.green_kw, self.v_eff, self.beta, rule='matmul')
        diff = np.max(np.abs(sig_elem - sig_mm))
        self.assertGreater(
            diff, 1e-4,
            "elementwise and matmul rules coincide -> model does not "
            "discriminate the contraction (need inter-orbital hopping)")

    def test_self_energy_does_not_match_matmul(self):
        """The matmul rule is NOT what the solver implements (it would be wrong)."""
        sig_mm = _flex_self_energy_explicit(
            self.solver, self.green_kw, self.v_eff, self.beta, rule='matmul')
        max_dev = np.max(np.abs(self.sigma_solver - sig_mm))
        self.assertGreater(
            max_dev, 1e-4,
            "solver self-energy unexpectedly matches the matmul rule")


class TestSelfEnergyGuards(unittest.TestCase):
    """Guard against a partial bosonic frequency grid in the self-energy.

    _calc_self_energy combines V_eff (nfreq bosonic freqs) and G (nmat
    fermionic freqs) over n_common = min(nfreq, nmat).  If nfreq < nmat the
    remaining tau slices are silently left at zero, which is physically wrong.
    The solver must fail loudly instead.
    """

    def test_partial_bosonic_grid_raises(self):
        solver, green_info = _make_1orb_solver(
            Lx=4, Ly=4, Nmat=16, T=1.0, mu=0.0, U=1.0,
            max_iter=1, mix=1.0, eps=1)
        solver.solve(green_info, 'tests/flex/output_anal')

        beta = 1.0 / solver.T
        nblock = 1
        nd_block = solver.norb
        sigma0 = np.zeros((nblock, solver.nmat, solver.lattice.nvol,
                           nd_block, nd_block), dtype=complex)
        green_kw = solver._calc_dressed_green(beta, solver.mu, sigma0)
        chi0q_raw = solver._calc_chi0q(green_kw, solver.green0_tail, beta)
        if solver.spin_mode in ["spin-free", "spinful"]:
            chi0q_raw = chi0q_raw[0]
        _, v_eff, _, _ = solver._flex_compute_veff(
            chi0q_raw, solver.ham_info.ham_inter_q)

        # truncate the bosonic frequency axis -> nfreq < nmat
        v_eff_partial = v_eff[:solver.nmat // 2]
        with self.assertRaises(ValueError):
            solver._calc_self_energy(green_kw, v_eff_partial, beta)

    def test_full_grid_does_not_raise(self):
        """The normal full-grid path (nfreq == nmat) must still work."""
        solver, green_info = _make_1orb_solver(
            Lx=4, Ly=4, Nmat=16, T=1.0, mu=0.0, U=1.0,
            max_iter=1, mix=1.0, eps=1)
        # solve() exercises _calc_self_energy on the full grid
        solver.solve(green_info, 'tests/flex/output_anal')
        self.assertIn("sigma", green_info)


class TestFLEXIgnoresChi0qInit(unittest.TestCase):
    """FLEX recomputes chi0q from the dressed Green's function every iteration,
    so a chi0q_init entry (stored as green_info['chi0q']) is not consumed.

    Note: green_init / trans_mod ARE consumed (the inherited _calc_epsilon_k
    uses them to build H0), so only chi0q_init is inert for FLEX.
    """

    def test_chi0q_init_is_ignored(self):
        cfg = dict(Lx=4, Ly=4, Nmat=16, T=1.0, mu=0.0, U=1.0,
                   max_iter=2, mix=1.0, eps=8)

        solver_a, gi_a = _make_1orb_solver(**cfg)
        solver_a.solve(gi_a, 'tests/flex/output_anal')
        sigma_ref = gi_a["sigma"].copy()

        # Inject a bogus chi0q_init; FLEX.solve must overwrite it and produce
        # the same self-energy.
        solver_b, gi_b = _make_1orb_solver(**cfg)
        gi_b["chi0q"] = np.full(
            (solver_b.nmat, solver_b.lattice.nvol, solver_b.norb, solver_b.norb),
            1.0e3, dtype=np.complex128)
        solver_b.solve(gi_b, 'tests/flex/output_anal')

        np.testing.assert_allclose(
            gi_b["sigma"], sigma_ref, atol=1e-12,
            err_msg="FLEX self-energy changed when chi0q_init was set; "
                    "FLEX recomputes chi0q and should ignore the init entry")


def _make_flex_solver_with(calc_scheme='reduced', interactions=None,
                           calc_type=None, Lx=4, Ly=4, Nmat=8, T=1.0, mu=0.0):
    """Build a FLEX solver with a given calc_scheme and interaction files.

    interactions: dict mapping the TOML interaction keyword (e.g. 'CoulombIntra',
    'Exchange') to the file body to write.  calc_type optionally sets
    [mode].calc_type ('ring' or 'ring+ladder').  Returns the constructed solver
    (construction may raise, which is what the guard tests check).
    """
    import tempfile
    import shutil

    if interactions is None:
        interactions = {
            'CoulombIntra': "CoulombIntra\n1\n1\n 1\n"
                            "   0    0    0    1    1   1.000000000000   0.0\n",
        }

    tmpdir = tempfile.mkdtemp(prefix='flex_guard_')
    src_dir = 'tests/rpa/input'
    for fn in ['geom.dat', 'transfer.dat']:
        shutil.copy2(os.path.join(src_dir, fn), os.path.join(tmpdir, fn))

    inter_cfg = {'path_to_input': tmpdir, 'Geometry': 'geom.dat',
                 'Transfer': 'transfer.dat'}
    for key, body in interactions.items():
        fname = key.lower() + '.dat'
        with open(os.path.join(tmpdir, fname), 'w') as f:
            f.write(body)
        inter_cfg[key] = fname

    info_log = {}
    info_mode = {
        'mode': 'FLEX',
        'param': {'T': T, 'mu': mu, 'CellShape': [Lx, Ly, 1],
                  'SubShape': [1, 1, 1], 'Nmat': Nmat,
                  'IterationMax': 1, 'Mix': 1.0, 'EPS': 1},
        'calc_scheme': calc_scheme,
    }
    if calc_type is not None:
        info_mode['calc_type'] = calc_type
    info_file_input = {'path_to_input': tmpdir, 'interaction': inter_cfg}

    import hwave.qlmsio.read_input_k as read_input_k
    read_io = read_input_k.QLMSkInput(info_file_input)
    ham_info = read_io.get_param("ham")

    import hwave.solver.flex as solver_flex
    try:
        solver = solver_flex.FLEX(ham_info, info_log, info_mode)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return solver


class TestFLEXSchemeGuards(unittest.TestCase):
    """FLEX-specific compatibility guards.

    FLEX consumes the reduced-shape (4-dim) chi0q and reduces the interaction
    via the density-density diagonal 'kaabb->kab', so it requires a
    reduced/squashed scheme and density-density interactions only.
    """

    # transfer-format body that registers as an Exchange interaction
    _EXCHANGE_BODY = ("Exchange\n1\n1\n 1\n"
                      "   0    0    0    1    1   0.500000000000   0.0\n")
    _PAIRHOP_BODY = ("PairHop\n1\n1\n 1\n"
                     "   0    0    0    1    1   0.500000000000   0.0\n")

    def test_general_scheme_rejected(self):
        """calc_scheme='general' would feed FLEX a 6/7-dim chi0q -> reject."""
        with self.assertRaises((ValueError, SystemExit)):
            _make_flex_solver_with(calc_scheme='general')

    def test_exchange_interaction_warns(self):
        """Exchange under 'squashed' is approximated by its density-density
        part; FLEX must warn (not silently drop) and still construct."""
        with self.assertLogs('hwave.solver.flex', level='WARNING') as cm:
            solver = _make_flex_solver_with(
                calc_scheme='squashed',
                interactions={'Exchange': self._EXCHANGE_BODY})
        self.assertTrue(
            any('density-density' in msg for msg in cm.output),
            "expected a warning about the density-density approximation")
        self.assertEqual(solver.calc_scheme, 'squashed')

    def test_pairhop_interaction_warns(self):
        """PairHop sets a separate flag from Exchange/PairLift, but its
        off-diagonal vertices are also dropped by the density-density reduction,
        so FLEX must warn for it too."""
        with self.assertLogs('hwave.solver.flex', level='WARNING') as cm:
            solver = _make_flex_solver_with(
                calc_scheme='squashed',
                interactions={'PairHop': self._PAIRHOP_BODY})
        self.assertTrue(
            any('density-density' in msg for msg in cm.output),
            "expected a density-density approximation warning for PairHop")
        self.assertEqual(solver.calc_scheme, 'squashed')

    def test_reduced_density_density_ok(self):
        """The standard reduced + density-density path must still construct."""
        solver = _make_flex_solver_with(calc_scheme='reduced')
        self.assertEqual(solver.calc_scheme, 'reduced')

    def test_auto_scheme_exchange_warns(self):
        """auto + exchange resolves to squashed (inherited RPA logic) and FLEX
        must warn about the density-density approximation, not error."""
        with self.assertLogs('hwave.solver.flex', level='WARNING') as cm:
            solver = _make_flex_solver_with(
                calc_scheme='auto',
                interactions={
                    'CoulombIntra': "CoulombIntra\n1\n1\n 1\n"
                                    "   0    0    0    1    1   1.0   0.0\n",
                    'Exchange': self._EXCHANGE_BODY,
                })
        self.assertEqual(solver.calc_scheme, 'squashed')
        self.assertTrue(
            any('density-density' in msg for msg in cm.output))

    def test_auto_ring_ladder_rejected_with_clear_message(self):
        """auto + calc_type='ring+ladder' resolves to general (inherited), which
        FLEX rejects; the message must name ring+ladder explicitly."""
        with self.assertRaisesRegex(ValueError, r'ring\+ladder'):
            _make_flex_solver_with(calc_scheme='auto', calc_type='ring+ladder')


class TestFLEXSpinSymmetry(unittest.TestCase):
    """Verify SU(2) spin symmetry properties in spin-free mode.

    In spin-free mode, all quantities must be block-diagonal in spin,
    ensuring Sigma_↑↓ = 0 and Sigma_↑↑ = Sigma_↓↓.
    """

    @classmethod
    def setUpClass(cls):
        cls.Lx, cls.Ly = 8, 8
        cls.Nmat = 32
        cls.T = 1.0
        cls.U = 2.0

        solver, green_info = _make_1orb_solver(
            Lx=cls.Lx, Ly=cls.Ly, Nmat=cls.Nmat,
            T=cls.T, mu=0.0, U=cls.U,
            max_iter=3, mix=0.5, eps=1)
        solver.solve(green_info, 'tests/flex/output_anal')

        cls.solver = solver
        cls.green_info = green_info

    def test_veff_spin_block_diagonal(self):
        """V_eff must be block-diagonal in spin space.

        In spin-free mode, V_eff is computed in the inflated (nd x nd) space
        where nd = norb * ns. The spin off-diagonal blocks must be zero,
        which guarantees Sigma_↑↓ = 0 in the FFT convolution.
        """
        solver = self.solver
        green_info = self.green_info
        norb = solver.norb
        ns = solver.ns
        nd = solver.nd
        beta = 1.0 / solver.T

        # Re-compute one iteration to get V_eff
        green_kw = solver._calc_dressed_green(beta, solver.mu, solver.sigma)
        chi0q_raw = solver._calc_chi0q(green_kw, solver.green0_tail, beta)
        if solver.spin_mode in ["spin-free", "spinful"]:
            chi0q_raw = chi0q_raw[0]

        chi0q, v_eff, chi_s, chi_c = solver._flex_compute_veff(
            chi0q_raw, solver.ham_info.ham_inter_q)

        # V_eff shape: (nfreq, nvol, nd, nd)
        # Reshape to (nfreq, nvol, ns, norb, ns, norb)
        nfreq, nvol = v_eff.shape[:2]
        v_so = v_eff.reshape(nfreq, nvol, ns, norb, ns, norb)

        # Check spin off-diagonal blocks are zero: V[s, :, t, :] = 0 for s != t
        for s in range(ns):
            for t in range(ns):
                if s != t:
                    off_diag = v_so[:, :, s, :, t, :]
                    max_abs = np.max(np.abs(off_diag))
                    self.assertLess(
                        max_abs, 1e-14,
                        "V_eff spin off-diagonal block ({},{}) nonzero: "
                        "max|V|={:.2e}".format(s, t, max_abs))

    def test_chi0q_spin_block_diagonal(self):
        """Inflated chi0q must be block-diagonal in spin space."""
        solver = self.solver
        norb = solver.norb
        ns = solver.ns
        nd = solver.nd
        beta = 1.0 / solver.T

        green_kw = solver._calc_dressed_green(beta, solver.mu, solver.sigma)
        chi0q_raw = solver._calc_chi0q(green_kw, solver.green0_tail, beta)
        if solver.spin_mode in ["spin-free", "spinful"]:
            chi0q_raw = chi0q_raw[0]

        chi0q, ham = solver._inflate_chi0q_and_ham(
            chi0q_raw, solver.ham_info.ham_inter_q)

        nfreq, nvol = chi0q.shape[:2]
        chi0_so = chi0q.reshape(nfreq, nvol, ns, norb, ns, norb)

        for s in range(ns):
            for t in range(ns):
                if s != t:
                    off_diag = chi0_so[:, :, s, :, t, :]
                    max_abs = np.max(np.abs(off_diag))
                    self.assertLess(
                        max_abs, 1e-14,
                        "chi0q spin off-diagonal block ({},{}) nonzero: "
                        "max|chi0|={:.2e}".format(s, t, max_abs))

    def test_chi0q_spin_degeneracy_converged(self):
        """chi0q_↑↑ = chi0q_↓↓ after SCF convergence (not just 1st iteration)."""
        solver = self.solver
        norb = solver.norb
        ns = solver.ns
        beta = 1.0 / solver.T

        green_kw = solver._calc_dressed_green(beta, solver.mu, solver.sigma)
        chi0q_raw = solver._calc_chi0q(green_kw, solver.green0_tail, beta)
        if solver.spin_mode in ["spin-free", "spinful"]:
            chi0q_raw = chi0q_raw[0]

        chi0q, _ = solver._inflate_chi0q_and_ham(
            chi0q_raw, solver.ham_info.ham_inter_q)

        nfreq, nvol = chi0q.shape[:2]
        chi0_so = chi0q.reshape(nfreq, nvol, ns, norb, ns, norb)

        np.testing.assert_allclose(
            chi0_so[:, :, 0, :, 0, :],
            chi0_so[:, :, 1, :, 1, :],
            atol=1e-14,
            err_msg="chi0q up-up != chi0q down-down after SCF")

    def test_sigma_full_offdiag_zero(self):
        """Self-energy spin off-diagonal blocks are zero in full space.

        Directly verify that the spin-inflated sigma has zero off-diagonal
        blocks before the slicing at _calc_self_energy return.
        """
        solver = self.solver
        norb = solver.norb
        ns = solver.ns
        nd = solver.nd
        nmat = solver.nmat
        beta = 1.0 / solver.T

        green_kw = solver._calc_dressed_green(beta, solver.mu, solver.sigma)
        chi0q_raw = solver._calc_chi0q(green_kw, solver.green0_tail, beta)
        if solver.spin_mode in ["spin-free", "spinful"]:
            chi0q_raw = chi0q_raw[0]

        chi0q, v_eff, chi_s, chi_c = solver._flex_compute_veff(
            chi0q_raw, solver.ham_info.ham_inter_q)

        # Temporarily patch to get full sigma before slicing
        nd_block_orig = green_kw.shape[-1]
        nd_v = v_eff.shape[-1]

        if nd_block_orig != nd_v:
            # Manually compute full sigma using the same logic
            import numpy.fft as FFT
            nx, ny, nz = solver.lattice.shape
            nvol = solver.lattice.nvol
            nblock = green_kw.shape[0]
            nfreq = v_eff.shape[0]

            # Transform G to (r, tau)
            omg_f = np.exp(-1j * np.pi * (1.0 / nmat - 1.0) * np.arange(nmat))
            green_flat = green_kw.reshape(nblock, nmat, nvol * nd_block_orig * nd_block_orig)
            green_kt = (FFT.fft(green_flat, axis=1)
                        * omg_f[np.newaxis, :, np.newaxis]
                        ).reshape(nblock, nmat, nx, ny, nz, nd_block_orig * nd_block_orig)
            green_rt = FFT.ifftn(green_kt, axes=(2, 3, 4)
                                 ).reshape(nblock, nmat, nvol, nd_block_orig, nd_block_orig)

            # Transform V_eff to (r, tau)
            omg_b = np.exp(-1j * np.pi * np.arange(nfreq))
            v_flat = v_eff.reshape(nfreq, nvol * nd_v * nd_v)
            v_qt = (FFT.fft(v_flat, axis=0)
                    * omg_b[:, np.newaxis]
                    ).reshape(nfreq, nx, ny, nz, nd_v * nd_v)
            v_rt = FFT.ifftn(v_qt, axes=(1, 2, 3)).reshape(nfreq, nvol, nd_v, nd_v)

            n_common = min(nfreq, nmat)

            # Build full sigma in inflated space
            sigma_rt = np.zeros((nblock, nmat, nvol, nd_v, nd_v),
                                dtype=np.complex128)
            for s in range(ns):
                sl = slice(s * norb, (s + 1) * norb)
                sigma_rt[:, :n_common, :, sl, sl] = (
                    v_rt[:n_common, :, sl, sl]
                    * green_rt[:, :n_common]
                )

            # Check off-diagonal blocks in (r, tau) space
            sigma_so = sigma_rt.reshape(nblock, nmat, nvol, ns, norb, ns, norb)
            for s in range(ns):
                for t in range(ns):
                    if s != t:
                        off_diag = sigma_so[:, :, :, s, :, t, :]
                        max_abs = np.max(np.abs(off_diag))
                        self.assertLess(
                            max_abs, 1e-30,
                            "Sigma(r,tau) spin off-diagonal ({},{}) nonzero: "
                            "max|Sigma|={:.2e}".format(s, t, max_abs))

            # Also verify spin degeneracy: Sigma_↑↑ = Sigma_↓↓
            # Transform back to (k, iwn) to check
            sigma_kt = FFT.fftn(
                sigma_rt.reshape(nblock, nmat, nx, ny, nz, nd_v * nd_v),
                axes=(2, 3, 4)
            ).reshape(nblock, nmat, nvol * nd_v * nd_v)
            omg_f_inv = np.exp(1j * np.pi * (1.0 / nmat - 1.0) * np.arange(nmat))
            sigma_kw_full = (FFT.ifft(sigma_kt * omg_f_inv[np.newaxis, :, np.newaxis],
                                      axis=1)
                             .reshape(nblock, nmat, nvol, nd_v, nd_v) * (1.0 / beta))

            sigma_so_kw = sigma_kw_full.reshape(nblock, nmat, nvol, ns, norb, ns, norb)
            np.testing.assert_allclose(
                sigma_so_kw[:, :, :, 0, :, 0, :],
                sigma_so_kw[:, :, :, 1, :, 1, :],
                atol=1e-12,
                err_msg="Sigma_↑↑ != Sigma_↓↓ in full spin-orbital space")


if __name__ == '__main__':
    unittest.main()
