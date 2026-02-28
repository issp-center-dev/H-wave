"""FLEX (Fluctuation Exchange Approximation) solver.

This module implements the FLEX approximation for itinerant electron systems.
FLEX extends RPA by using dressed (self-consistent) Green's functions instead
of bare ones, iterating G -> chi0 -> chi -> Sigma -> G until convergence.

The solver inherits from the RPA class to reuse infrastructure for:
- Lattice and interaction setup
- Diagonalization and chemical potential search
- Bare susceptibility calculation
- RPA equation solving
- Block-diagonal structure detection
"""

from __future__ import annotations
from typing import Optional

import sys
import os
import numpy as np
import numpy.fft as FFT
from requests.structures import CaseInsensitiveDict

try:
    from .perf import do_profile
except ImportError:
    from functools import wraps
    def do_profile(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)
        return wrapper

import logging
logger = logging.getLogger(__name__)

from .rpa import RPA, Lattice, Interaction


class FLEX(RPA):
    """FLEX solver that extends RPA with self-consistent Green's functions.

    The FLEX approximation computes the self-energy from spin and charge
    fluctuations, then updates the Green's function self-consistently.

    The SCF loop is:
        1. Compute G(k, iwn) from H0 and Sigma
        2. Compute chi0(q, ivn) = -T/Nk sum_k G(k+q) G(k)
        3. Inflate chi0 and ham to common spin-orbital space
        4. Decompose ham into spin/charge channels
        5. Compute chi_s, chi_c via RPA equations
        6. Construct V_eff from chi_s, chi_c, chi0
        7. Compute Sigma(k, iwn) = T/Nk sum_q V_eff(q) G(k-q)
        8. Check convergence; if not, go to 1

    Parameters
    ----------
    param_ham : dict
        Hamiltonian parameters.
    info_log : dict
        Logging configuration.
    info_mode : dict
        Calculation mode parameters. FLEX-specific parameters include:
        - IterationMax: Maximum SCF iterations (default: 100)
        - Mix: Mixing parameter for self-energy update (default: 0.2)
        - EPS: Convergence criterion exponent (default: 6, i.e., 1e-6)
    """

    @do_profile
    def __init__(self, param_ham, info_log, info_mode):
        logger.debug(">>> FLEX.__init__")

        # Initialize RPA infrastructure (lattice, interaction, params)
        super().__init__(param_ham, info_log, info_mode)

        # FLEX-specific parameters
        self._init_flex_param()

    def _init_flex_param(self):
        """Initialize FLEX-specific parameters."""
        logger.debug(">>> FLEX._init_flex_param")

        self.max_iter = int(self.param_mod.get("IterationMax", 100))
        self.mix = float(self.param_mod.get("Mix", 0.2))

        eps_exp = self.param_mod.get("EPS", 6)
        if isinstance(eps_exp, float) and eps_exp < 1.0:
            self.eps = eps_exp
        else:
            self.eps = 10.0 ** (-int(eps_exp))

        logger.info("FLEX parameters:")
        logger.info("    max_iter        = {}".format(self.max_iter))
        logger.info("    mix             = {}".format(self.mix))
        logger.info("    eps             = {:e}".format(self.eps))

    @do_profile
    def solve(self, green_info, path_to_output):
        """Solve the FLEX equations self-consistently.

        Parameters
        ----------
        green_info : dict
            Dictionary containing Green's function information.
        path_to_output : str
            Path to output directory.
        """
        logger.info("Start FLEX calculations")

        beta = 1.0 / self.T
        nvol = self.lattice.nvol
        nmat = self.nmat
        norb = self.norb
        ns = self.ns
        nd = self.nd

        # Step 1: Compute band structure and chemical potential
        self._calc_epsilon_k(green_info)

        if self.calc_mu:
            if self.spin_mode == "spin-free":
                Ncond = self.Ncond / 2
            else:
                Ncond = self.Ncond
            dist, mu = self._find_mu(Ncond, self.T)
        else:
            mu = self.mu_value

        self.mu = mu

        # Step 2: Compute bare Green's function G0(k, iwn)
        green0, green0_tail = self._calc_green(beta, mu)

        # Store for reference
        self.green0 = green0
        self.green0_tail = green0_tail

        # Initialize self-energy to zero
        # Shape: (nblock, nmat, nvol, nd_block, nd_block)
        nblock = green0.shape[0]
        nd_block = green0.shape[-1]
        sigma = np.zeros((nblock, nmat, nvol, nd_block, nd_block),
                         dtype=np.complex128)

        # Prepare interaction Hamiltonian (full spin-orbital space)
        ham_orig = self.ham_info.ham_inter_q

        # Main SCF loop
        converged = False
        for iteration in range(self.max_iter):
            logger.info("FLEX iteration {}/{}".format(iteration + 1, self.max_iter))

            # Step 3: Compute dressed Green's function G(k, iwn)
            green_kw = self._calc_dressed_green(beta, mu, sigma)

            # Step 4: Compute chi0(q, ivn) from dressed G
            chi0q_raw = self._calc_chi0q(green_kw, green0_tail, beta)

            # Remove spin block dimension
            if self.spin_mode in ["spin-free", "spinful"]:
                assert chi0q_raw.shape[0] == 1
                chi0q_raw = chi0q_raw[0]
            else:
                assert chi0q_raw.shape[0] == 2

            # Step 5: Inflate chi0q and ham to common reduced space,
            # then compute spin/charge susceptibilities and V_eff
            chi0q_out, v_eff, chi_s, chi_c = self._flex_compute_veff(
                chi0q_raw, ham_orig)

            # Step 6: Compute self-energy Sigma(k, iwn)
            sigma_new = self._calc_self_energy(green_kw, v_eff, beta)

            # Step 7: Mix and check convergence
            diff = self._calc_convergence(sigma, sigma_new)
            logger.info("  convergence: |dSigma|/|Sigma| = {:.3e}".format(diff))

            sigma = (1.0 - self.mix) * sigma + self.mix * sigma_new

            if diff < self.eps:
                logger.info("FLEX converged after {} iterations".format(
                    iteration + 1))
                converged = True
                break

        if not converged:
            logger.warning("FLEX did not converge after {} iterations "
                           "(diff={:.3e}, eps={:.3e})".format(
                               self.max_iter, diff, self.eps))

        # Store results
        self.sigma = sigma
        self.green_kw = green_kw
        self.chi_s = chi_s
        self.chi_c = chi_c

        # Store in green_info for output
        green_info["chi0q"] = chi0q_out
        green_info["chiq_s"] = chi_s
        green_info["chiq_c"] = chi_c
        green_info["sigma"] = sigma
        green_info["green"] = green_kw

        logger.info("End FLEX calculations")

    @do_profile
    def _calc_dressed_green(self, beta, mu, sigma):
        """Compute dressed Green's function G(k, iwn) = [G0^{-1} - Sigma]^{-1}.

        Parameters
        ----------
        beta : float
            Inverse temperature.
        mu : float
            Chemical potential.
        sigma : ndarray
            Self-energy, shape (nblock, nmat, nvol, nd, nd).

        Returns
        -------
        ndarray
            Dressed Green's function, shape (nblock, nmat, nvol, nd, nd).
        """
        logger.debug(">>> FLEX._calc_dressed_green")

        ew = self.H0_eigenvalue
        ev = self.H0_eigenvector

        nblock, nvol, nd = ew.shape
        nmat = self.nmat

        # Matsubara frequencies
        iomega = (np.arange(nmat) * 2 + 1 - nmat) * np.pi / beta

        # Reconstruct H0 in orbital basis from eigendecomposition
        # H0 = ev @ diag(ew) @ ev†, using matmul for BLAS efficiency
        H0_k = np.matmul(ev * ew[:, :, np.newaxis, :], np.conj(ev).swapaxes(-2, -1))

        # G^{-1}(k, iwn) = (iwn + mu) * I - H0(k) - Sigma(k, iwn)
        eye = np.eye(nd, dtype=np.complex128)

        # Vectorized construction of G^{-1} for all frequencies
        # iomega shape: (nmat,) -> broadcast to (1, nmat, 1, 1, 1)
        iw = 1j * iomega[np.newaxis, :, np.newaxis, np.newaxis, np.newaxis]
        # H0_k shape: (nblock, nvol, nd, nd) -> (nblock, 1, nvol, nd, nd)
        H0_exp = H0_k[:, np.newaxis, :, :, :]

        green_inv = (iw + mu) * eye - H0_exp - sigma

        # G(k, iwn) = [G^{-1}]^{-1}
        green = np.linalg.inv(green_inv)

        return green

    @do_profile
    def _flex_compute_veff(self, chi0q_raw, ham_orig):
        """Inflate chi0q, decompose into spin/charge channels, and compute V_eff.

        This method handles the spin inflation (same as RPA.solve()) and then
        decomposes the interaction into spin and charge channels for FLEX.

        Parameters
        ----------
        chi0q_raw : ndarray
            Raw chi0q from _calc_chi0q (before spin inflation).
        ham_orig : ndarray
            Original interaction Hamiltonian in full spin-orbital space.

        Returns
        -------
        chi0q_inflated : ndarray
            Chi0q after spin inflation (for output).
        v_eff : ndarray
            Effective FLEX interaction V_eff(q, ivn).
        chi_s : ndarray
            Spin susceptibility.
        chi_c : ndarray
            Charge susceptibility.
        """
        nvol = self.lattice.nvol
        norb = self.norb
        ns = self.ns
        nd = self.nd

        # Inflate chi0q to reduced spin-orbital space (same logic as RPA.solve())
        chi0q, ham = self._inflate_chi0q_and_ham(chi0q_raw, ham_orig)

        # Decompose ham into spin and charge channels
        ham_s, ham_c = self._build_spin_charge_vertices(ham)

        # Solve RPA for spin channel: chi_s = [1 - chi0*V_s]^{-1} chi0
        chi_s = self._solve_rpa(chi0q, ham_s)

        # Solve RPA for charge channel: chi_c = [1 + chi0*V_c]^{-1} chi0
        chi_c = self._solve_rpa(chi0q, ham_c)

        # Compute V_eff
        v_eff = self._calc_veff(chi0q, chi_s, chi_c, ham)

        return chi0q, v_eff, chi_s, chi_c

    @do_profile
    def _inflate_chi0q_and_ham(self, chi0q_raw, ham_orig):
        """Apply spin inflation to chi0q and ham, matching RPA.solve() logic.

        Takes raw chi0q (orbital-only for spin-free, spin-block for spin-diag,
        or spin-orbital for spinful) and inflates to the reduced space used
        for FLEX calculations.

        Parameters
        ----------
        chi0q_raw : ndarray
            Raw chi0q before spin inflation.
        ham_orig : ndarray
            Full interaction Hamiltonian.

        Returns
        -------
        chi0q : ndarray
            Inflated chi0q in reduced spin-orbital space.
        ham : ndarray
            Interaction Hamiltonian in matching reduced space.
        """
        nvol = self.lattice.nvol
        norb = self.norb
        ns = self.ns
        nd = self.nd

        if self.spin_mode == "spin-free":
            # chi0q_raw shape: (nmat, nvol, norb, norb) for reduced
            nfreq = chi0q_raw.shape[0]

            # Inflate to spin-orbital reduced space
            spin_tensor = np.identity(ns)
            chi0q = np.einsum('lkab,st->lksatb',
                              chi0q_raw.reshape(nfreq, nvol, norb, norb),
                              spin_tensor).reshape(nfreq, nvol, nd, nd)

            ham = np.einsum('ksasatbtb->ksatb',
                            ham_orig.reshape(nvol, *(ns, norb) * 4)
                            ).reshape(nvol, nd, nd)

        elif self.spin_mode == "spin-diag":
            # chi0q_raw shape: (nblock=2, nmat, nvol, norb, norb) for reduced
            nblock_s, nfreq, nvol_c, norb1, norb2 = chi0q_raw.shape

            spin_tensor = np.identity(ns)
            chi0q = np.einsum('glkab,gh->lkgahb',
                              chi0q_raw,
                              spin_tensor).reshape(nfreq, nvol, nd, nd)

            ham = np.einsum('ksasatbtb->ksatb',
                            ham_orig.reshape(nvol, *(ns, norb) * 4)
                            ).reshape(nvol, nd, nd)

        elif self.spin_mode == "spinful":
            # chi0q_raw shape: (nmat, nvol, nd, nd) for reduced
            chi0q = chi0q_raw

            ham = np.einsum('kaabb->kab',
                            ham_orig.reshape(nvol, *(nd,) * 4)
                            ).reshape(nvol, nd, nd)

        return chi0q, ham

    @do_profile
    def _build_spin_charge_vertices(self, ham_inflated):
        """Build spin and charge interaction vertices from inflated ham.

        The inflated ham is in reduced (nd x nd) form. We decompose it into
        spin and charge channels.

        For the Hubbard model with U on-site:
            ham_inflated has structure:
                [W_↑↑  W_↑↓]     where W_↑↑ = same-spin, W_↑↓ = cross-spin
                [W_↓↑  W_↓↓]

            Physical susceptibilities:
                chi_s = chi0 * [1 - U_s chi0]^{-1}  (spin, Stoner-enhanced)
                chi_c = chi0 * [1 + U_c chi0]^{-1}  (charge, suppressed)
            where:
                U_s = W_cross - W_same   (for Hubbard: U - 0 = U)
                U_c = W_cross + W_same   (for Hubbard: U + 0 = U)

            _solve_rpa convention: [1 + chi0 * ham]^{-1} chi0
                => ham_s = -U_s  (to get [1 - U_s chi0]^{-1})
                => ham_c = +U_c  (to get [1 + U_c chi0]^{-1})

        Parameters
        ----------
        ham_inflated : ndarray, shape (nvol, nd, nd)
            Interaction in reduced spin-orbital space.

        Returns
        -------
        ham_s : ndarray, shape (nvol, nd, nd)
            Spin channel vertex (with sign for _solve_rpa convention).
        ham_c : ndarray, shape (nvol, nd, nd)
            Charge channel vertex.
        """
        logger.debug(">>> FLEX._build_spin_charge_vertices")

        nvol = self.lattice.nvol
        norb = self.norb
        ns = self.ns
        nd = self.nd

        # ham_inflated is (nvol, nd, nd) where nd = norb * ns
        # Reshape to (nvol, ns, norb, ns, norb)
        ham_so = ham_inflated.reshape(nvol, ns, norb, ns, norb)

        # Same-spin block: (s, a) -> (s, a) same s
        # Cross-spin block: (s, a) -> (t, a) with s != t
        # Average over both spins for symmetry:
        w_same = 0.5 * (ham_so[:, 0, :, 0, :] + ham_so[:, 1, :, 1, :])
        w_cross = 0.5 * (ham_so[:, 1, :, 0, :] + ham_so[:, 0, :, 1, :])

        # U_s = w_cross - w_same (for Hubbard: U - 0 = U)
        # U_c = w_cross + w_same (for Hubbard: U + 0 = U)
        u_s = w_cross - w_same  # (nvol, norb, norb)
        u_c = w_cross + w_same  # (nvol, norb, norb)

        # For _solve_rpa [1 + chi0 * ham]^{-1} chi0:
        #   ham_s = -U_s  -> [1 - U_s chi0]^{-1} chi0  (Stoner enhancement)
        #   ham_c = +U_c  -> [1 + U_c chi0]^{-1} chi0  (charge suppression)
        spin_id = np.eye(ns)
        ham_s = np.einsum('kab,st->ksatb', -u_s, spin_id).reshape(nvol, nd, nd)
        ham_c = np.einsum('kab,st->ksatb', u_c, spin_id).reshape(nvol, nd, nd)

        return ham_s, ham_c

    @do_profile
    def _calc_veff(self, chi0q, chi_s, chi_c, ham_inflated):
        """Compute effective FLEX interaction V_eff(q, ivn).

        V_eff = W * [3/2*(chi_s - chi0) + 1/2*(chi_c - chi0)] * W

        The subtraction of chi0 removes the double-counted second-order diagram.

        Parameters
        ----------
        chi0q : ndarray
            Bare susceptibility (inflated), shape (nfreq, nvol, nd, nd).
        chi_s : ndarray
            Spin susceptibility, same shape.
        chi_c : ndarray
            Charge susceptibility, same shape.
        ham_inflated : ndarray
            Interaction in reduced space, shape (nvol, nd, nd).

        Returns
        -------
        ndarray
            Effective interaction V_eff, shape (nfreq, nvol, nd, nd).
        """
        logger.debug(">>> FLEX._calc_veff")

        nfreq = chi0q.shape[0]
        nvol = chi0q.shape[1]
        nd = chi0q.shape[-1]

        # Reshape to matrices
        chi0q_2d = chi0q.reshape(nfreq, nvol, nd, nd)
        chi_s_2d = chi_s.reshape(nfreq, nvol, nd, nd)
        chi_c_2d = chi_c.reshape(nfreq, nvol, nd, nd)
        ham_2d = ham_inflated.reshape(nvol, nd, nd)

        # Fluctuation susceptibility
        fluct_chi = 1.5 * (chi_s_2d - chi0q_2d) + 0.5 * (chi_c_2d - chi0q_2d)

        # V_eff = W * fluct_chi * W
        # Use batched matmul instead of einsum for better BLAS utilization
        # W @ fluct_chi: broadcast (nvol, nd, nd) @ (nfreq, nvol, nd, nd)
        tmp = np.matmul(ham_2d, fluct_chi)
        # tmp @ W: (nfreq, nvol, nd, nd) @ (nvol, nd, nd)
        v_eff = np.matmul(tmp, ham_2d)

        return v_eff

    @do_profile
    def _calc_self_energy(self, green_kw, v_eff, beta):
        """Compute self-energy Sigma(k, iwn) via FFT convolution.

        Sigma(k, iwn) = T/Nk * sum_{q,m} V_eff(q, ivm) * G(k-q, iwn-ivm)

        The convolution is computed efficiently using FFT:
        1. Transform G and V_eff to (r, tau) space
        2. Multiply in (r, tau) space: Sigma(r,tau) = V(r,tau) * G(r,tau)
        3. Transform back to (k, iwn) space

        Parameters
        ----------
        green_kw : ndarray
            Dressed Green's function, shape (nblock, nmat, nvol, nd_block, nd_block).
        v_eff : ndarray
            Effective interaction, shape (nfreq, nvol, nd, nd).
        beta : float
            Inverse temperature.

        Returns
        -------
        ndarray
            Self-energy, shape (nblock, nmat, nvol, nd_block, nd_block).
        """
        logger.debug(">>> FLEX._calc_self_energy")

        nx, ny, nz = self.lattice.shape
        nvol = self.lattice.nvol
        nmat = self.nmat
        nblock = green_kw.shape[0]
        nd_block = green_kw.shape[-1]
        nd_v = v_eff.shape[-1]
        nfreq = v_eff.shape[0]

        # --- Transform Green's function to (r, tau) space ---

        # Matsubara freq -> imaginary time for G (fermionic)
        # Use broadcasting instead of einsum for phase multiplication
        omg_f = np.exp(-1j * np.pi * (1.0 / nmat - 1.0) * np.arange(nmat))
        green_flat = green_kw.reshape(nblock, nmat, nvol * nd_block * nd_block)
        green_kt = (FFT.fft(green_flat, axis=1)
                     * omg_f[np.newaxis, :, np.newaxis]
                     ).reshape(nblock, nmat, nx, ny, nz, nd_block * nd_block)

        # k-space -> real-space for G
        green_rt = FFT.ifftn(green_kt, axes=(2, 3, 4)
                             ).reshape(nblock, nmat, nvol, nd_block, nd_block)

        # --- Transform V_eff to (r, tau) space ---

        # Bosonic Matsubara freq -> imaginary time
        # Bosonic phase: (-1)^j = exp(-i*pi*j)
        omg_b = np.exp(-1j * np.pi * np.arange(nfreq))
        v_flat = v_eff.reshape(nfreq, nvol * nd_v * nd_v)
        v_qt = (FFT.fft(v_flat, axis=0)
                * omg_b[:, np.newaxis]
                ).reshape(nfreq, nx, ny, nz, nd_v * nd_v)

        # q-space -> real-space for V_eff
        v_rt = FFT.ifftn(v_qt, axes=(1, 2, 3)).reshape(nfreq, nvol, nd_v, nd_v)

        # --- Compute Sigma(r, tau) ---
        n_common = min(nfreq, nmat)

        # If nd_block != nd_v, we need to expand G to match V_eff space
        # spin-free: nd_block = norb, nd_v = norb*ns = nd
        if nd_block != nd_v:
            # Expand Green's function to spin-orbital space using np.kron
            # For spin-free: G_{sa,tb} = G_{ab} * delta_{st}
            norb = self.norb
            ns = self.ns
            # Vectorized spin inflation: use Kronecker product with identity
            # green_rt[:, :n_common] has shape (nblock, n_common, nvol, norb, norb)
            # We need kron(eye(ns), G) for each (nblock, time, vol)
            spin_eye = np.eye(ns, dtype=np.complex128)
            # sigma_rt[g, t, r] = v_rt[t, r] * kron(I_s, G[g, t, r])
            # Compute directly without full expansion:
            # kron(I_s, G) is block-diagonal with G repeated ns times on diagonal
            # Hadamard with v_rt: only diagonal blocks of v_rt contribute
            sigma_rt = np.zeros((nblock, nmat, nvol, nd_v, nd_v),
                                dtype=np.complex128)
            for s in range(ns):
                sl = slice(s * norb, (s + 1) * norb)
                sigma_rt[:, :n_common, :, sl, sl] = (
                    v_rt[:n_common, :, sl, sl]
                    * green_rt[:, :n_common]
                )
        else:
            # Sigma(r,tau) = V_eff(r,tau) * G(r,tau)  (element-wise product)
            sigma_rt = np.zeros((nblock, nmat, nvol, nd_v, nd_v),
                                dtype=np.complex128)
            sigma_rt[:, :n_common] = v_rt[:n_common] * green_rt[:, :n_common]

        # --- Transform Sigma back to (k, iwn) space ---

        # Real-space -> k-space
        sigma_kt = FFT.fftn(
            sigma_rt.reshape(nblock, nmat, nx, ny, nz, nd_v * nd_v),
            axes=(2, 3, 4)
        ).reshape(nblock, nmat, nvol * nd_v * nd_v)

        # Imaginary time -> Matsubara freq (fermionic)
        # Use broadcasting instead of einsum for phase multiplication
        omg_f_inv = np.exp(1j * np.pi * (1.0 / nmat - 1.0) * np.arange(nmat))
        sigma_kw = (FFT.ifft(sigma_kt * omg_f_inv[np.newaxis, :, np.newaxis],
                             axis=1)
                    .reshape(nblock, nmat, nvol, nd_v, nd_v) * (1.0 / beta))

        # If we expanded to spin-orbital space, contract back to block space
        if nd_block != nd_v:
            return sigma_kw[:, :, :, :nd_block, :nd_block]

        return sigma_kw

    def _calc_convergence(self, sigma_old, sigma_new):
        """Calculate convergence criterion for self-energy.

        Parameters
        ----------
        sigma_old : ndarray
            Previous self-energy.
        sigma_new : ndarray
            New self-energy.

        Returns
        -------
        float
            Relative difference |sigma_new - sigma_old| / |sigma_new|.
        """
        diff = np.linalg.norm(sigma_new - sigma_old)
        norm = np.linalg.norm(sigma_new)
        if norm < 1.0e-30:
            return diff
        return diff / norm

    @do_profile
    def save_results(self, info_outputfile, green_info):
        """Save FLEX calculation results.

        Parameters
        ----------
        info_outputfile : dict
            Output file configuration.
        green_info : dict
            Calculation results.
        """
        logger.info("Save FLEX results")
        path_to_output = info_outputfile["path_to_output"]

        self._init_wavevec()

        # Save chi0q
        if "chi0q" in info_outputfile:
            file_name = os.path.join(path_to_output, info_outputfile["chi0q"])
            np.savez(file_name,
                     chi0q=green_info["chi0q"],
                     freq_index=self.freq_index,
                     wavevector_unit=self.kvec,
                     wavevector_index=self.wavenum_table)
            logger.info("save_results: save chi0q in file {}".format(file_name))

        # Save susceptibilities (spin and charge channels separately for Eliashberg)
        common_meta = dict(freq_index=self.freq_index,
                           wavevector_unit=self.kvec,
                           wavevector_index=self.wavenum_table)

        if "chiq_s" in green_info:
            file_name = os.path.join(path_to_output,
                                     info_outputfile.get("chiq_s", "chiq_s"))
            np.savez(file_name, chiq_s=green_info["chiq_s"], **common_meta)
            logger.info("save_results: save chiq_s in file {}".format(file_name))

        if "chiq_c" in green_info:
            file_name = os.path.join(path_to_output,
                                     info_outputfile.get("chiq_c", "chiq_c"))
            np.savez(file_name, chiq_c=green_info["chiq_c"], **common_meta)
            logger.info("save_results: save chiq_c in file {}".format(file_name))

        if "chiq" in info_outputfile:
            file_name = os.path.join(path_to_output, info_outputfile["chiq"])
            save_dict = dict(**common_meta)
            if "chiq_s" in green_info:
                save_dict["chiq_s"] = green_info["chiq_s"]
            if "chiq_c" in green_info:
                save_dict["chiq_c"] = green_info["chiq_c"]
            np.savez(file_name, **save_dict)
            logger.info("save_results: save chiq in file {}".format(file_name))

        # Save self-energy
        if "sigma" in info_outputfile:
            file_name = os.path.join(path_to_output, info_outputfile["sigma"])
            np.savez(file_name,
                     sigma=green_info.get("sigma"),
                     freq_index=self.freq_index,
                     wavevector_unit=self.kvec,
                     wavevector_index=self.wavenum_table)
            logger.info("save_results: save sigma in file {}".format(file_name))

        # Save Green's function
        if "green" in info_outputfile:
            file_name = os.path.join(path_to_output, info_outputfile["green"])
            np.savez(file_name,
                     green=green_info.get("green"),
                     freq_index=self.freq_index,
                     wavevector_unit=self.kvec,
                     wavevector_index=self.wavenum_table)
            logger.info("save_results: save green in file {}".format(file_name))
