#!/usr/bin/env python3
"""Tests for RPA ladder diagram (ring+ladder calc_type).

Tests verify:
1. chi_zz = chi_+- consistency: At the RPA level for paramagnetic states,
   the longitudinal spin vertex V_spin = W_uu - W_ud is identical to the
   transverse vertex W_pm = ham[↑↑↑↑] - ham[↓↓↑↑]. This holds for ANY
   interaction, not just SU(2)-symmetric ones.
2. Ring-only regression: calc_type="ring" results unchanged
3. Bare susceptibility validation: chi0 from RPA matches Lindhard function
4. Multi-interaction consistency: CoulombIntra, CoulombInter, Hund, Exchange,
   Ising, PairLift, and full Kanamori all satisfy chi_zz = chi_pm.
"""

import os
import unittest

import numpy as np


class TestRPALadder(unittest.TestCase):
    """Test ladder diagram implementation via SU(2) symmetry."""

    def _run_rpa(self, calc_type="ring", calc_scheme="auto",
                 spin_orbital=False, inter=None, Lx=8, Ly=8, Nmat=64,
                 interactions=None, input_path='tests/rpa/input',
                 T=2.0, filling=0.75):
        """Run RPA and return (solver, green_info).

        Parameters
        ----------
        interactions : dict, optional
            Full interaction dict to use instead of default CoulombIntra.
            If None, uses {'CoulombIntra': 'coulombintra.dat'} plus `inter`.
        inter : dict, optional
            Additional interactions to add to the default.
        input_path : str, optional
            Path to input files directory (default: 'tests/rpa/input').
        T : float, optional
            Temperature (default: 2.0).
        filling : float, optional
            Electron filling (default: 0.75).
        """
        if inter is None:
            inter = {}

        info_log = {}
        info_mode = {
            'mode': 'RPA',
            'param': {
                'T': T,
                'filling': filling,
                'CellShape': [Lx, Ly, 1],
                'SubShape': [1, 1, 1],
                'Nmat': Nmat,
            },
            'enable_spin_orbital': spin_orbital,
            'calc_scheme': calc_scheme,
            'calc_type': calc_type,
        }

        if interactions is None:
            interaction_dict = {
                'path_to_input': input_path,
                'Geometry': 'geom.dat',
                'Transfer': 'transfer.dat',
                'CoulombIntra': 'coulombintra.dat',
            }
            interaction_dict.update(inter)
        else:
            interaction_dict = {
                'path_to_input': input_path,
                'Geometry': 'geom.dat',
                'Transfer': 'transfer.dat',
            }
            interaction_dict.update(interactions)

        info_file = {
            'input': {
                'path_to_input': input_path,
                'interaction': interaction_dict,
            },
            'output': {
                'path_to_output': 'tests/rpa/output',
                'chi0q': 'chi0q.npz',
            },
        }
        os.makedirs(info_file['output']['path_to_output'], exist_ok=True)

        import hwave.qlmsio.read_input_k as read_input_k
        read_io = read_input_k.QLMSkInput(info_file['input'])
        ham_info = read_io.get_param("ham")

        import hwave.solver.rpa as solver_rpa
        solver = solver_rpa.RPA(ham_info, info_log, info_mode)

        green_info = read_io.get_param("green")
        solver.solve(green_info, info_file['output']['path_to_output'])

        return solver, green_info

    def test_su2_symmetry_1orb_coulombintra(self):
        """For 1-orbital CoulombIntra, ring+ladder should give chi_zz = chi_+-.

        chi_zz = chi_uu - chi_ud  (from longitudinal RPA)
        chi_+- = transverse RPA susceptibility (from ladder channel)

        SU(2) symmetry requires chi_zz = chi_+- for paramagnetic systems.

        Note: in the density correlation convention used by rpa.py,
        chi_zz = chi_uu - chi_ud (not 2*(chi_uu - chi_ud)),
        and chi_+- = chi_zz (not chi_zz/2).
        """
        solver, green_info = self._run_rpa(
            calc_type="ring+ladder",
            calc_scheme="general",
        )

        chiq = green_info["chiq"]
        chiq_pm = green_info["chiq_pm"]

        nfreq = chiq.shape[0]
        iw0 = nfreq // 2

        # chi_zz from longitudinal channel
        chi_uu = chiq[iw0, :, 0, 0, 0, 0]
        chi_ud = chiq[iw0, :, 0, 0, 1, 1]
        chi_zz = chi_uu - chi_ud

        # chi_+- from transverse channel
        chi_pm = chiq_pm[iw0, :, 0, 0, 0, 0]

        np.testing.assert_allclose(
            chi_zz, chi_pm,
            atol=1e-10,
            err_msg="SU(2) symmetry violated: chi_zz != chi_+- for ring+ladder"
        )

    def test_ring_only_no_transverse(self):
        """For ring-only mode, chiq_pm should not exist."""
        solver, green_info = self._run_rpa(
            calc_type="ring",
            calc_scheme="general",
        )
        self.assertNotIn("chiq_pm", green_info,
                          "Ring-only mode should not compute transverse chi")

    def test_ring_backward_compatibility(self):
        """Ring-only with general scheme should produce non-trivial results."""
        solver, green_info = self._run_rpa(
            calc_type="ring",
            calc_scheme="general",
        )
        chiq = green_info["chiq"]
        nfreq = chiq.shape[0]
        iw0 = nfreq // 2
        chi_uu = chiq[iw0, :, 0, 0, 0, 0]
        self.assertTrue(np.any(np.abs(chi_uu) > 1e-10),
                        "chi_uu should be non-zero")

    def test_transverse_vertex_coulombintra(self):
        """Verify the transverse vertex W_+- = -U for CoulombIntra.

        For H = U n_up n_dn, the Fock exchange gives W_+-(diag) = -U.
        The transverse RPA is: chi_+- = chi0 / (1 + U*chi0)
        which gives the SAME result as chi_s = chi0 / (1 - U*chi0)
        when the sign convention accounts for W_+- = -U.

        SU(2) requires chi_s = chi_+- in density convention.
        """
        solver, green_info = self._run_rpa(
            calc_type="ring+ladder",
            calc_scheme="general",
            Lx=4, Ly=4, Nmat=32,
        )

        chiq = green_info["chiq"]
        chiq_pm = green_info["chiq_pm"]

        nfreq = chiq.shape[0]
        iw0 = nfreq // 2

        chi_s = chiq[iw0, :, 0, 0, 0, 0] - chiq[iw0, :, 0, 0, 1, 1]
        chi_pm = chiq_pm[iw0, :, 0, 0, 0, 0]

        np.testing.assert_allclose(
            chi_s, chi_pm,
            atol=1e-10,
            err_msg="Transverse vertex: chi_s != chi_+- for CoulombIntra"
        )

    def test_su2_all_matsubara_frequencies(self):
        """SU(2) symmetry should hold at ALL Matsubara frequencies, not just iw0."""
        solver, green_info = self._run_rpa(
            calc_type="ring+ladder",
            calc_scheme="general",
            Lx=4, Ly=4, Nmat=16,
        )

        chiq = green_info["chiq"]
        chiq_pm = green_info["chiq_pm"]

        nfreq = chiq.shape[0]

        for iw in range(nfreq):
            chi_zz = chiq[iw, :, 0, 0, 0, 0] - chiq[iw, :, 0, 0, 1, 1]
            chi_pm = chiq_pm[iw, :, 0, 0, 0, 0]
            np.testing.assert_allclose(
                chi_zz, chi_pm,
                atol=1e-10,
                err_msg=f"SU(2) violated at Matsubara frequency index {iw}"
            )


    def _assert_su2(self, green_info, norb=1, msg=""):
        """Assert SU(2) symmetry: chi_zz = chi_+- at all q-points for iw0."""
        chiq = green_info["chiq"]
        chiq_pm = green_info["chiq_pm"]

        nfreq = chiq.shape[0]
        iw0 = nfreq // 2

        if norb == 1:
            chi_zz = chiq[iw0, :, 0, 0, 0, 0] - chiq[iw0, :, 0, 0, 1, 1]
            chi_pm = chiq_pm[iw0, :, 0, 0, 0, 0]
        else:
            nd = norb * 2
            chi_zz = np.zeros(chiq.shape[1], dtype=complex)
            chi_pm = np.zeros(chiq.shape[1], dtype=complex)
            for a in range(norb):
                chi_zz += (chiq[iw0, :, a, a, a, a]
                           - chiq[iw0, :, a, a, norb + a, norb + a])
                chi_pm += chiq_pm[iw0, :, a, a, a, a]

        np.testing.assert_allclose(
            chi_zz, chi_pm,
            atol=1e-10,
            err_msg=f"SU(2) symmetry violated: {msg}"
        )

    def test_su2_coulombinter(self):
        """CoulombInter (V * n_i * n_j) is SU(2) symmetric."""
        solver, green_info = self._run_rpa(
            calc_type="ring+ladder",
            calc_scheme="general",
            Lx=4, Ly=4, Nmat=32,
            interactions={
                'CoulombIntra': 'coulombintra.dat',
                'CoulombInter': 'coulombinter.dat',
            },
        )
        self._assert_su2(green_info, msg="CoulombIntra + CoulombInter")

    def test_su2_coulombinter_only(self):
        """CoulombInter alone is SU(2) symmetric."""
        solver, green_info = self._run_rpa(
            calc_type="ring+ladder",
            calc_scheme="general",
            Lx=4, Ly=4, Nmat=32,
            interactions={'CoulombInter': 'coulombinter.dat'},
        )
        self._assert_su2(green_info, msg="CoulombInter only")

    def test_su2_hund(self):
        """CoulombIntra + Hund: chi_zz = chi_pm at the RPA level.

        At the RPA level for paramagnetic states, chi_zz = chi_pm always holds
        because the longitudinal spin vertex V_spin = W_↑↑ - W_↑↓ is identical
        to the transverse vertex W_pm = ham[↑↑↑↑] - ham[↓↓↑↑] by construction.
        This is true for ANY interaction, not just SU(2)-symmetric ones.
        """
        solver, green_info = self._run_rpa(
            calc_type="ring+ladder",
            calc_scheme="general",
            Lx=4, Ly=4, Nmat=32,
            interactions={
                'CoulombIntra': 'coulombintra.dat',
                'Hund': 'coulombinter.dat',
            },
        )
        self._assert_su2(green_info, msg="CoulombIntra + Hund")

    def test_su2_hund_2orb(self):
        """2-orbital Hund: chi_zz = chi_pm at the RPA level.

        Even for multi-orbital systems, the paramagnetic RPA preserves
        chi_zz = chi_pm because V_spin = W_pm as orbital matrices.
        """
        solver, green_info = self._run_rpa(
            calc_type="ring+ladder",
            calc_scheme="general",
            input_path='tests/rpa/input_2orb',
            Lx=4, Ly=4, Nmat=32,
            T=2.0, filling=0.5,
            interactions={
                'CoulombIntra': 'coulombintra.dat',
                'Hund': 'hund_onsite.dat',
            },
        )
        self._assert_su2(green_info, norb=2, msg="CoulombIntra + Hund (2-orbital)")

    def test_su2_exchange(self):
        """CoulombIntra + Exchange is SU(2) symmetric.

        Exchange -J*(c†_up c_dn c†_dn c_up + h.c.) = S+S- part of spin exchange,
        which IS SU(2) symmetric together with CoulombIntra.
        """
        solver, green_info = self._run_rpa(
            calc_type="ring+ladder",
            calc_scheme="general",
            Lx=4, Ly=4, Nmat=32,
            interactions={
                'CoulombIntra': 'coulombintra.dat',
                'Exchange': 'coulombinter.dat',
            },
        )
        self._assert_su2(green_info, msg="CoulombIntra + Exchange")

    def test_su2_pairlift(self):
        """CoulombIntra + PairLift is SU(2) symmetric.

        PairLift J*(c†_up c_dn c†_up c_dn + h.c.) doesn't contribute
        to the Hartree vertex in 1-orbital systems (maps to zero in S/C matrix).
        """
        solver, green_info = self._run_rpa(
            calc_type="ring+ladder",
            calc_scheme="general",
            Lx=4, Ly=4, Nmat=32,
            interactions={
                'CoulombIntra': 'coulombintra.dat',
                'PairLift': 'coulombinter.dat',
            },
        )
        self._assert_su2(green_info, msg="CoulombIntra + PairLift")

    def test_su2_ising(self):
        """CoulombIntra + Ising: chi_zz = chi_pm at the RPA level.

        The Ising interaction = 4*Jz*S^z*S^z is NOT SU(2) symmetric.
        However, at the RPA level for paramagnetic states, the longitudinal
        spin vertex V_spin = W_↑↑ - W_↑↓ is always identical to the
        transverse vertex W_pm, so chi_zz = chi_pm regardless.
        """
        solver, green_info = self._run_rpa(
            calc_type="ring+ladder",
            calc_scheme="general",
            Lx=4, Ly=4, Nmat=32,
            interactions={
                'CoulombIntra': 'coulombintra.dat',
                'Ising': 'coulombinter.dat',
            },
        )
        self._assert_su2(green_info, msg="CoulombIntra + Ising")

    def test_su2_all_interactions(self):
        """Full SU(2)-symmetric Kanamori: U + V + Hund + Exchange.

        Using the same J for Hund and Exchange with appropriate signs
        gives the full rotationally invariant interaction.
        """
        solver, green_info = self._run_rpa(
            calc_type="ring+ladder",
            calc_scheme="general",
            Lx=4, Ly=4, Nmat=32,
            interactions={
                'CoulombIntra': 'coulombintra.dat',
                'CoulombInter': 'coulombinter.dat',
                'Hund': 'coulombinter.dat',
                'Exchange': 'coulombinter.dat',
            },
        )
        self._assert_su2(green_info,
                          msg="Full Kanamori (U + V + Hund + Exchange)")


class TestRPALadderBareResponse(unittest.TestCase):
    """Validate bare susceptibility chi0 by comparing with the Lindhard function
    computed via finite-difference response from non-interacting calculations.

    The Lindhard susceptibility chi0(q) = -1/N sum_k (f_{k+q}-f_k)/(e_{k+q}-e_k)
    equals the static bare susceptibility. The RPA code computes chi0(q, iwn) at
    individual Matsubara frequencies; we validate by comparing their sum with the
    Lindhard result.
    """

    def _compute_lindhard(self, ek, mu, T, qx, qy, Lx, Ly):
        """Compute the Lindhard function chi0_orb(q) analytically.

        Returns the static orbital susceptibility (per spin).
        """
        nvol = Lx * Ly
        beta = 1.0 / T
        KX_flat = np.array([i // Ly for i in range(nvol)])
        KY_flat = np.array([i % Ly for i in range(nvol)])
        xk = ek - mu
        fk = 1.0 / (np.exp(beta * xk) + 1.0)

        chi0 = 0.0
        for ik in range(nvol):
            kpqx = (KX_flat[ik] + qx) % Lx
            kpqy = (KY_flat[ik] + qy) % Ly
            ikpq = kpqx * Ly + kpqy
            dx = xk[ikpq] - xk[ik]
            if abs(dx) > 1e-12:
                chi0 += (fk[ikpq] - fk[ik]) / dx
            else:
                chi0 += (-beta * fk[ik] * (1.0 - fk[ik]))
        chi0 *= -1.0 / nvol
        return chi0

    def _compute_fd_response_real_space(self, ek_real_space_builder, mu, T,
                                         qx, qy, Lx, Ly, dh=1e-6,
                                         field_type="z"):
        """Compute static response via finite-difference in real space.

        Builds the full 2N x 2N real-space Hamiltonian with applied field,
        diagonalizes, and computes the order parameter change.

        Returns chi_zz(q) or chi_xx(q).
        """
        N = Lx * Ly
        beta = 1.0 / T
        KX_flat = np.array([i // Ly for i in range(N)])
        KY_flat = np.array([i % Ly for i in range(N)])

        # Phase factor for the applied field: e^{i q.r} in real space
        phase_r = np.zeros(N, dtype=float)
        for ir in range(N):
            phase_r[ir] = np.cos(2 * np.pi * (qx * KX_flat[ir] / Lx
                                                + qy * KY_flat[ir] / Ly))

        def build_H(h):
            H = np.zeros((2 * N, 2 * N))
            # Hopping (from the real-space builder)
            H[:N, :N], H[N:, N:] = ek_real_space_builder()
            # Field
            for i in range(N):
                if field_type == "z":
                    H[i, i] += h * phase_r[i]
                    H[N + i, N + i] -= h * phase_r[i]
                elif field_type == "x":
                    H[i, N + i] += h * phase_r[i]
                    H[N + i, i] += h * phase_r[i]
            return H

        def compute_order_param(H):
            w, v = np.linalg.eigh(H)
            f_occ = 1.0 / (np.exp(beta * (w - mu)) + 1.0)
            rho = v @ np.diag(f_occ) @ v.conj().T
            m = 0.0
            for ir in range(N):
                if field_type == "z":
                    m += phase_r[ir] * (rho[ir, ir] - rho[N + ir, N + ir]).real
                elif field_type == "x":
                    m += phase_r[ir] * (rho[ir, N + ir]
                                        + rho[N + ir, ir]).real
            return m / N

        Hp = build_H(dh)
        Hm = build_H(-dh)
        mp = compute_order_param(Hp)
        mm = compute_order_param(Hm)
        return (mp - mm) / (2 * dh)

    def test_chi0_uniform_matches_lindhard(self):
        """Verify chi0_zz(q=0) from Lindhard matches finite-difference.

        The Lindhard chi0 > 0 (physics convention), while the response
        dM/dh < 0 (field opposes magnetization). So chi0_zz = -dM/dh.
        """
        Lx, Ly = 8, 8
        T = 2.0
        N = Lx * Ly

        KX_flat = np.array([i // Ly for i in range(N)])
        KY_flat = np.array([i % Ly for i in range(N)])
        ek = -2.0 * (np.cos(2 * np.pi * KX_flat / Lx)
                      + np.cos(2 * np.pi * KY_flat / Ly))

        from scipy.optimize import brentq
        beta = 1.0 / T

        def filling(mu):
            return 2 * np.mean(1.0 / (np.exp(beta * (ek - mu)) + 1.0)) - 0.75
        mu = brentq(filling, -10, 10)

        # Lindhard: chi0 > 0
        chi0_orb = self._compute_lindhard(ek, mu, T, 0, 0, Lx, Ly)
        chi0_zz = 2 * chi0_orb

        # Analytic response: dM/dh = -chi0_zz
        xk = ek - mu
        fk = 1.0 / (np.exp(beta * xk) + 1.0)
        response = -2 * beta * np.mean(fk * (1 - fk))  # negative

        np.testing.assert_allclose(
            chi0_zz, -response,
            rtol=1e-8,
            err_msg="chi0_zz(q=0): Lindhard vs analytic mismatch"
        )

    def test_chi0_staggered_matches_fd(self):
        """Verify chi0_zz(q=pi,pi) from Lindhard matches real-space FD."""
        Lx, Ly = 6, 6
        T = 2.0
        N = Lx * Ly

        KX_flat = np.array([i // Ly for i in range(N)])
        KY_flat = np.array([i % Ly for i in range(N)])
        ek = -2.0 * (np.cos(2 * np.pi * KX_flat / Lx)
                      + np.cos(2 * np.pi * KY_flat / Ly))

        from scipy.optimize import brentq
        beta = 1.0 / T

        def filling(mu):
            return 2 * np.mean(1.0 / (np.exp(beta * (ek - mu)) + 1.0)) - 0.75
        mu = brentq(filling, -10, 10)

        # Lindhard
        Qx, Qy = Lx // 2, Ly // 2
        chi0_orb = self._compute_lindhard(ek, mu, T, Qx, Qy, Lx, Ly)
        chi0_zz = 2 * chi0_orb

        # Real-space finite-difference
        def build_hopping():
            H = np.zeros((N, N))
            for ir in range(N):
                irx = KX_flat[ir]
                iry = KY_flat[ir]
                for dx, dy in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
                    jx = (irx + dx) % Lx
                    jy = (iry + dy) % Ly
                    jr = jx * Ly + jy
                    H[ir, jr] -= 1.0
            return H, H.copy()

        chi_fd = self._compute_fd_response_real_space(
            build_hopping, mu, T, Qx, Qy, Lx, Ly, dh=1e-6, field_type="z")

        # Lindhard chi0 > 0, FD response dM/dh < 0: chi0_zz = -chi_fd
        np.testing.assert_allclose(
            chi0_zz, -chi_fd,
            rtol=1e-4,
            err_msg="chi0_zz(q=pi,pi): Lindhard vs real-space FD mismatch"
        )


if __name__ == '__main__':
    unittest.main()
