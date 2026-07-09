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

        # FLEX consumes the reduced-shape (4-dim) chi0q and reduces the
        # interaction via the density-density diagonal ('kaabb->kab') on the
        # reduced/squashed path.  The 'general' scheme selects the paramagnetic
        # full-vertex multi-orbital path (v1: spin-free only; the spin_mode
        # guard is enforced in solve(), where spin_mode is determined).
        scheme = self.calc_scheme.lower()
        if scheme == "general":
            # FLEX general is the paramagnetic full-vertex path: it sums the
            # RPA *ring* (bubble) series with the full rank-4 vertex. The
            # transverse *ladder* channel (calc_type='ring+ladder') is not part
            # of this path, so reject it even though it forces scheme='general'.
            if getattr(self, "calc_type", "ring") == "ring+ladder":
                msg = ("FLEX does not support calc_type='ring+ladder' (the "
                       "transverse ladder channel); FLEX general is ring-only. "
                       "Use the default calc_type='ring' with "
                       "calc_scheme='general'.")
                logger.error(msg)
                raise ValueError(msg)
            if getattr(self.ham_info, "enable_spin_orbital", False):
                raise ValueError(
                    "calc_scheme='general' FLEX (v1) does not support "
                    "enable_spin_orbital; deferred to the generalized FLEX "
                    "solver.")
            self._flex_general = True
        elif scheme in ("reduced", "squashed"):
            self._flex_general = False
        else:
            if getattr(self, "calc_type", "ring") == "ring+ladder":
                msg = ("FLEX does not support calc_type='ring+ladder' (the "
                       "transverse ladder channel); use the default 'ring' "
                       "with calc_scheme='reduced'/'squashed'/'general'.")
            else:
                msg = ("FLEX requires calc_scheme='reduced', 'squashed', or "
                       "'general', got '{}'.".format(self.calc_scheme))
            logger.error(msg)
            raise ValueError(msg)
        if not self._flex_general and (
                self.ham_info.has_interaction_exchange()
                or self.ham_info.has_interaction_pairhop()):
            # FLEX reduces the vertex via the density-density diagonal
            # ('kaabb->kab'), so exchange/spin-flip/pair off-diagonal vertices
            # are dropped.  This is a deliberate (common) approximation, but
            # warn so it is not silent (cf. the inherited reduced+exchange
            # guard, which does not cover the squashed scheme).  Exchange and
            # PairLift set the exchange flag; PairHop sets a separate flag, so
            # both are checked here to cover all off-diagonal interaction types.
            logger.warning(
                "FLEX uses the density-density reduction; exchange- and "
                "pair-hopping-type interactions (Exchange, PairLift, PairHop) "
                "are approximated by their density-density part "
                "(off-diagonal vertices are dropped). "
                "Use calc_scheme='general' to keep them.")

        self.max_iter = int(self.param_mod.get("IterationMax", 100))
        self.mix = float(self.param_mod.get("Mix", 0.2))

        eps_exp = self.param_mod.get("EPS", 6)
        if isinstance(eps_exp, float) and eps_exp < 1.0:
            self.eps = eps_exp
        else:
            self.eps = 10.0 ** (-int(eps_exp))

        # Lazy cache for the q-independent MYO S/C matrices (general path only);
        # built once on first use in _inflate_chi0q_and_ham_general and reused
        # across SCF iterations. None until then.
        self._myo_sc_cache = None

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

        Notes
        -----
        FLEX starts the SCF loop from zero self-energy and recomputes chi0q from
        the dressed Green's function every iteration, so a `chi0q_init` entry
        loaded by the inherited RPA `read_init` is NOT consumed.  `green_init`
        and `trans_mod`, however, ARE consumed: the inherited `_calc_epsilon_k`
        uses them to build the (mean-field-shifted) transfer H0(k), exactly as
        in RPA.
        """
        logger.info("Start FLEX calculations")

        beta = 1.0 / self.T
        nvol = self.lattice.nvol
        nmat = self.nmat
        norb = self.norb
        ns = self.ns
        nd = self.nd

        # Step 1: Compute band structure and chemical potential.
        # _calc_epsilon_k determines self.spin_mode from H0(k); the general
        # full-vertex path (v1) is paramagnetic only, so guard here once
        # spin_mode is known.
        self._calc_epsilon_k(green_info)

        if self._flex_general and self.spin_mode != "spin-free":
            raise ValueError(
                "calc_scheme='general' FLEX (v1) supports spin_mode='spin-free' "
                "only, got '{}'. spin-diag/spinful are deferred to the "
                "generalized FLEX solver.".format(self.spin_mode))

        if self.calc_mu:
            # spin-free counts one spin, so the target is halved (as in
            # RPA._find_mu / RPA.solve).  Ncond_target is reused every SCF
            # iteration to re-solve mu from the DRESSED Green's function.
            if self.spin_mode == "spin-free":
                Ncond_target = self.Ncond / 2
            else:
                Ncond_target = self.Ncond
            # initial mu from the non-interacting bands (Sigma = 0)
            dist, mu = self._find_mu(Ncond_target, self.T)
        else:
            Ncond_target = None
            mu = self.mu_value

        self.mu = mu

        # Step 2: Compute bare Green's function G0(k, iwn)
        green0, green0_tail = self._calc_green(beta, mu)

        # Store for reference
        self.green0 = green0
        self.green0_tail = green0_tail

        # High-frequency tail contract (coeff_tail): RPA's tail acceleration
        # is a two-step pair -- _calc_green subtracts aa/(i w_n) in FREQUENCY
        # space, and _calc_chi0q subtracts the analytic tau-space constant
        # green0_tail (= VV† aa beta/2) after the Matsubara FFT.  The dressed
        # Green's function below comes from _calc_dressed_green, which returns
        # the FULL physical G, so the frequency-space term must be subtracted
        # here before handing G to _calc_chi0q; otherwise G(tau) is uniformly
        # shifted by -aa/2 and chi0q is O(1) wrong.
        #
        # The tail subtraction is applied ONLY to the chi0q transform (the one
        # paired with green0_tail).  The self-energy convolution keeps the
        # full physical G: its FFT pipeline is an exact cyclic frequency
        # convolution and needs no tau-space tail reconstruction -- measured
        # on the 8x8 Hubbard fixture, reconstructing the "true" G(tau) there
        # does not improve the Nmat convergence of Sigma.
        aa = self.coeff_tail
        if aa != 0.0:
            iomega = (np.arange(nmat) * 2 + 1 - nmat) * np.pi / beta
            ev = self.H0_eigenvector
            VVt = ev @ np.conj(ev).swapaxes(-2, -1)  # (nblock, nvol, nd, nd)
            green_tail_w = ((aa / (1j * iomega))[np.newaxis, :, np.newaxis,
                                                 np.newaxis, np.newaxis]
                            * VVt[:, np.newaxis])
        else:
            green_tail_w = None

        # Initialize self-energy to zero
        # Shape: (nblock, nmat, nvol, nd_block, nd_block)
        nblock = green0.shape[0]
        nd_block = green0.shape[-1]
        sigma = np.zeros((nblock, nmat, nvol, nd_block, nd_block),
                         dtype=np.complex128)

        # Prepare interaction Hamiltonian (full spin-orbital space)
        ham_orig = self.ham_info.ham_inter_q

        # Main SCF loop
        diff = float("inf")
        converged = False
        for iteration in range(self.max_iter):
            logger.info("FLEX iteration {}/{}".format(iteration + 1, self.max_iter))

            # Re-solve mu so the DRESSED G reproduces the target particle
            # number: as Sigma grows the frozen non-interacting mu no longer
            # yields Ncond electrons, so the run would otherwise converge to a
            # different filling than requested.  mu is solved for the current
            # sigma BEFORE building green_kw so the stored (mu, green_kw) pair
            # is self-consistent (N(green_kw) == Ncond).  A fixed mu
            # (calc_mu=False) is left untouched.
            if self.calc_mu:
                mu = self._find_mu_dressed(sigma, beta, Ncond_target)
                self.mu = mu

            # Step 3: Compute dressed Green's function G(k, iwn)
            green_kw = self._calc_dressed_green(beta, mu, sigma)

            # Tail-subtracted G for the Matsubara-FFT transforms (see the
            # coeff_tail contract above); identical to green_kw when aa == 0.
            if green_tail_w is not None:
                green_scf = green_kw - green_tail_w
            else:
                green_scf = green_kw

            # Step 4: Compute chi0(q, ivn) from dressed G
            chi0q_raw = self._calc_chi0q(green_scf, green0_tail, beta)

            # Remove spin block dimension
            if self.spin_mode in ["spin-free", "spinful"]:
                assert chi0q_raw.shape[0] == 1
                chi0q_raw = chi0q_raw[0]
            else:
                assert chi0q_raw.shape[0] == 2

            # Step 5: Inflate chi0q and ham, then compute spin/charge
            # susceptibilities and V_eff.  The general (paramagnetic full-vertex)
            # path keeps the full rank-4 orbital vertex; the reduced/squashed
            # path uses the density-density reduction.
            # Step 6: Compute self-energy Sigma(k, iwn).
            if self._flex_general:
                chi0q_out, v_eff, chi_s, chi_c = \
                    self._flex_compute_veff_general(chi0q_raw, ham_orig)
                sigma_new = self._calc_self_energy_general(green_kw, v_eff, beta)
            else:
                chi0q_out, v_eff, chi_s, chi_c = self._flex_compute_veff(
                    chi0q_raw, ham_orig)
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

        if self.max_iter == 0:
            # No SCF iteration ran: green_kw / chi_s / chi_c / chi0q_out were
            # never computed.  Warn-and-return instead of dereferencing them.
            logger.warning("FLEX IterationMax=0: no SCF step performed; "
                           "no results stored.")
            return

        # Final-output consistency: during the loop green_kw was built from the
        # PRE-mix sigma, while `sigma` below is the POST-mix estimate, so the
        # stored (mu, green, sigma) triple would not satisfy the Dyson equation
        # green = [G0^{-1} - sigma]^{-1} (noticeable for non-converged or
        # IterationMax=1 runs; negligible once converged).  Rebuild the dressed
        # G from the final stored sigma -- and, for calc_mu, re-solve mu for it
        # -- so the stored triple is mutually consistent and N(green) == Ncond.
        if self.calc_mu:
            mu = self._find_mu_dressed(sigma, beta, Ncond_target)
            self.mu = mu
        green_kw = self._calc_dressed_green(beta, mu, sigma)

        # Physical observables from the final (consistent) dressed G: particle
        # number N and spin Sz.  In fixed-mu mode this N is the mu-N single
        # point; in calc_mu mode it equals the target Ncond.
        physics = self._calc_physics_dressed(green_kw, mu, beta)
        self.physics = physics
        logger.info("FLEX: NCond = {}, Sz = {}, ChemicalPotential = {}".format(
            physics["NCond"], physics["Sz"], physics["mu"]))

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
        green_info["physics"] = physics

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
    def _calc_number_dressed(self, sigma, mu, beta):
        """Particle number carried by the dressed Green's function at ``mu``.

        The self-energy shifts and renormalizes the dressed G, so the actual
        electron count differs from the non-interacting Fermi count.  The count
        is evaluated with the standard tail-subtracted Matsubara sum: the
        non-interacting reference is summed analytically (Fermi function) and
        only the small, absolutely-convergent difference ``G - G0`` is summed
        over the finite Matsubara grid::

            N(mu) = sum_{block,k,a} f(eps_a(k) - mu)
                    + (1/beta) sum_{k,n} Tr[ G(k,iwn) - G0(k,iwn) ]

        Here ``G0(k,iwn) = [(iwn+mu)I - H0(k)]^{-1}`` (same mu), whose analytic
        occupation is exactly the Fermi term, and ``G - G0 = G0 Sigma G`` decays
        as ``1/(iwn)^2`` so the residual sum needs no e^{iwn 0+} convergence
        factor and is real.  At ``Sigma = 0`` the residual vanishes and ``N``
        reduces exactly to the ``_find_mu`` Fermi count.

        The total is compared against the SAME target convention as
        :meth:`RPA._find_mu`: the sum runs over all spin blocks, so for
        spin-free (one block) it is the one-spin count (target ``Ncond/2``).

        Parameters
        ----------
        sigma : ndarray
            Self-energy, shape (nblock, nmat, nvol, nd_block, nd_block).
        mu : float
            Trial chemical potential.
        beta : float
            Inverse temperature.

        Returns
        -------
        float
            Particle number N(mu).
        """
        nmat = self.nmat
        ew = self.H0_eigenvalue                       # (nblock, nvol, nd_block)

        # analytic non-interacting reference (Fermi function)
        n_ref = self._fermi_occupation(1.0 / beta, mu, ew).sum()

        # dressed correction Tr[G - G0], summed over the finite Matsubara grid
        green = self._calc_dressed_green(beta, mu, sigma)
        iomega = (np.arange(nmat) * 2 + 1 - nmat) * np.pi / beta
        trG = np.einsum('bnkaa->bnk', green)          # (nblock, nmat, nvol)
        trG0 = (1.0 / ((1j * iomega)[np.newaxis, :, np.newaxis, np.newaxis]
                       + (mu - ew)[:, np.newaxis, :, :])).sum(axis=-1)
        corr = (trG - trG0).sum() / beta

        return n_ref + corr.real

    @staticmethod
    def _fermi_occupation(t, mu, ev, ene_cutoff=1.0e2):
        """Fermi function with the same overflow guard as RPA._find_mu."""
        w = (ev - mu) / t
        mask = w < ene_cutoff
        w1 = np.where(mask, w, 0.0)
        v1 = 1.0 / (1.0 + np.exp(w1))
        return np.where(mask, v1, 0.0)

    @do_profile
    def _calc_occupation_dressed(self, green_kw, mu, beta):
        r"""k-summed occupation per (spin block, orbital) from the dressed G.

        The per-orbital analogue of :meth:`_calc_number_dressed`: the same
        tail-subtracted Matsubara sum, resolved on each orbital ``a`` instead of
        traced.  The non-interacting reference is projected onto the orbital
        basis with the H0 eigenvectors ``V``::

            n_a(k) = sum_j |V_aj|^2 f(eps_j(k) - mu)
                     + (1/beta) sum_n [ G_aa(k,iwn) - G0_aa(k,iwn) ]
            G0_aa(k,iwn) = sum_j |V_aj|^2 / (iwn + mu - eps_j(k))

        Summing the result over orbitals recovers ``_calc_number_dressed`` (by
        unitarity ``sum_a |V_aj|^2 = 1``); it is split per orbital here so N and
        Sz can be assembled from spin-resolved partial sums.

        Parameters
        ----------
        green_kw : ndarray
            Dressed Green's function, shape (nblock, nmat, nvol, nd_block, nd_block).
        mu : float
            Chemical potential.
        beta : float
            Inverse temperature.

        Returns
        -------
        ndarray
            Occupation summed over k, shape (nblock, nd_block).
        """
        nmat = self.nmat
        ew = self.H0_eigenvalue                       # (nblock, nvol, nd_block)
        ev = self.H0_eigenvector                      # (nblock, nvol, a, j)
        vsq = np.abs(ev) ** 2                          # |V_aj|^2

        # bare per-orbital reference: n0_a = sum_j |V_aj|^2 f(eps_j - mu)
        f = self._fermi_occupation(1.0 / beta, mu, ew)     # (nblock, nvol, j)
        n0 = np.einsum('bkaj,bkj->bka', vsq, f)            # (nblock, nvol, a)

        # dressed correction (1/beta) sum_n [G_aa - G0_aa]
        iomega = (np.arange(nmat) * 2 + 1 - nmat) * np.pi / beta
        g_aa = np.einsum('bnkaa->bnka', green_kw)          # (nblock, nmat, nvol, a)
        denom = ((1j * iomega)[np.newaxis, :, np.newaxis, np.newaxis]
                 + (mu - ew)[:, np.newaxis, :, :])         # (nblock, nmat, nvol, j)
        g0_aa = np.einsum('bkaj,bnkj->bnka', vsq, 1.0 / denom)
        corr = (g_aa - g0_aa).sum(axis=1) / beta           # sum over n -> (nblock, nvol, a)

        nocc = n0 + corr.real
        return nocc.sum(axis=1)                            # sum over k -> (nblock, nd_block)

    @do_profile
    def _calc_physics_dressed(self, green_kw, mu, beta):
        """Assemble the particle number N and spin Sz from the dressed G.

        Uses :meth:`_calc_occupation_dressed` and the FLEX spin-block orbital
        ordering (``s*norb + a``).  Conventions per spin mode:

        - spin-free (one block computed): ``N = 2 * sum(nocc)``, ``Sz = 0``.
        - spin-diag (block 0 = up, block 1 = down):
          ``N = N_up + N_down``, ``Sz = (N_up - N_down)/2``.
        - spinful (single block, ``nd = 2*norb``, up = orbitals ``[0, norb)``,
          down = ``[norb, 2*norb)``): ``N = sum(nocc)``,
          ``Sz = (N_up - N_down)/2``.

        Parameters
        ----------
        green_kw : ndarray
            Dressed Green's function (nblock, nmat, nvol, nd_block, nd_block).
        mu : float
            Chemical potential.
        beta : float
            Inverse temperature.

        Returns
        -------
        dict
            ``{"NCond": N_total, "Sz": Sz, "mu": mu}`` (plain floats).
        """
        nocc = self._calc_occupation_dressed(green_kw, mu, beta)
        norb = self.norb

        if self.spin_mode == "spin-free":
            n_total = 2.0 * nocc.sum()
            sz = 0.0
        elif self.spin_mode == "spin-diag":
            n_up = nocc[0].sum()
            n_down = nocc[1].sum()
            n_total = n_up + n_down
            sz = 0.5 * (n_up - n_down)
        else:  # spinful: spin folded into the orbital index (spin-block order)
            n_up = nocc[0, :norb].sum()
            n_down = nocc[0, norb:].sum()
            n_total = nocc.sum()
            sz = 0.5 * (n_up - n_down)

        return {"NCond": float(n_total.real if np.iscomplexobj(n_total)
                               else n_total),
                "Sz": float(sz.real if np.iscomplexobj(sz) else sz),
                "mu": float(mu)}

    @do_profile
    def _matsubara_number_operator(self, sigma, beta):
        r"""Eigenvalues of the mu-independent part of ``G^{-1}``, for the mu
        search.

        During the chemical-potential search ``sigma`` is held FIXED, so

            G^{-1}(k, iwn; mu) = (iwn + mu) I - H0(k) - Sigma(k, iwn)
                               = M(k, iwn) + mu I,
            M(k, iwn) = iwn I - H0(k) - Sigma(k, iwn)   (mu-independent).

        Diagonalizing ``M`` ONCE then gives ``Tr[G(mu)] = sum_j 1/(lam_j + mu)``
        for ANY trial ``mu`` -- one eigenvalue decomposition per iteration
        instead of one full matrix inversion per bisection step (the mu search
        was ~50% of solve()).  ``M`` is non-Hermitian, but its eigenvalues sit
        at ``lam_j ~ iwn - (real band+Sigma)``, i.e. ``|Im(lam_j)| ~ |wn| >=
        pi/beta > 0``, so ``lam_j + mu`` (mu real) never touches zero: the
        1/(lam+mu) sum is well conditioned.

        Parameters
        ----------
        sigma : ndarray
            Self-energy, shape (nblock, nmat, nvol, nd_block, nd_block).
        beta : float
            Inverse temperature.

        Returns
        -------
        lam : ndarray
            Eigenvalues of ``M``, shape (nblock, nmat, nvol, nd_block).
        ew : ndarray
            H0 band energies ``self.H0_eigenvalue`` (returned for convenience,
            used by :meth:`_number_from_eigs` for the analytic reference).
        """
        ew = self.H0_eigenvalue
        ev = self.H0_eigenvector
        nblock, nvol, nd = ew.shape
        nmat = self.nmat

        iomega = (np.arange(nmat) * 2 + 1 - nmat) * np.pi / beta
        H0_k = np.matmul(ev * ew[:, :, np.newaxis, :],
                         np.conj(ev).swapaxes(-2, -1))       # (nb, nvol, nd, nd)
        eye = np.eye(nd, dtype=np.complex128)
        iw = 1j * iomega[np.newaxis, :, np.newaxis, np.newaxis, np.newaxis]
        M = iw * eye - H0_k[:, np.newaxis, :, :, :] - sigma  # (nb,nmat,nvol,nd,nd)

        lam = np.linalg.eigvals(M)                            # (nb,nmat,nvol,nd)
        return lam, ew

    def _number_from_eigs(self, lam, ew, mu, beta, with_deriv=False):
        """Particle number N(mu) (and optionally dN/dmu) from precomputed
        ``M`` eigenvalues.

        Cheap closure evaluated for every trial ``mu`` in the mu search:
        identical formula and result as :meth:`_calc_number_dressed` (verified
        to machine precision in the tests), but reusing the eigenvalues from
        :meth:`_matsubara_number_operator` instead of re-inverting G::

            N(mu) = sum_{block,k,a} f(eps_a - mu)
                    + (1/beta) sum_{k,n} [ sum_j 1/(lam_j + mu)
                                           - sum_a 1/(iwn + mu - eps_a) ]

        Because ``sigma`` (hence ``lam``) is FIXED during the search, the
        derivative is analytic and cheap to evaluate from the same
        eigenvalues, which lets the search use Newton's method::

            dN/dmu = sum_a (1/T) f_a (1 - f_a)
                     + (1/beta) sum_{k,n} [ -sum_j 1/(lam_j + mu)^2
                                            + sum_a 1/(iwn + mu - eps_a)^2 ]

        Parameters
        ----------
        lam : ndarray
            ``M`` eigenvalues from :meth:`_matsubara_number_operator`,
            shape (nblock, nmat, nvol, nd_block).
        ew : ndarray
            H0 band energies, shape (nblock, nvol, nd_block).
        mu : float
            Trial chemical potential.
        beta : float
            Inverse temperature.
        with_deriv : bool, optional
            If True, return the tuple ``(N, dN/dmu)`` instead of ``N``.

        Returns
        -------
        float or tuple of float
            ``N(mu)``, or ``(N(mu), dN/dmu)`` when ``with_deriv`` is True.
        """
        nmat = self.nmat
        T = 1.0 / beta
        iomega = (np.arange(nmat) * 2 + 1 - nmat) * np.pi / beta

        f = self._fermi_occupation(T, mu, ew)
        n_ref = f.sum()

        denomG = lam + mu                                    # (nb, nmat, nvol, nd)
        denomG0 = ((1j * iomega)[np.newaxis, :, np.newaxis, np.newaxis]
                   + (mu - ew)[:, np.newaxis, :, :])
        trG = (1.0 / denomG).sum(axis=-1)                    # (nb, nmat, nvol)
        trG0 = (1.0 / denomG0).sum(axis=-1)
        n = n_ref + ((trG - trG0).sum() / beta).real

        if not with_deriv:
            return n

        # dN/dmu: the Fermi part is +f(1-f)/T; each 1/(x+mu) term contributes
        # -1/(x+mu)^2 (same trG - trG0 combination, both differentiated).
        dref = ((1.0 / T) * f * (1.0 - f)).sum()
        dtrG = (-1.0 / denomG ** 2).sum(axis=-1)
        dtrG0 = (-1.0 / denomG0 ** 2).sum(axis=-1)
        dn = dref + ((dtrG - dtrG0).sum() / beta).real
        return n, dn

    @do_profile
    def _find_mu_dressed(self, sigma, beta, Ncond):
        """Re-solve mu so the dressed Green's function carries ``Ncond``.

        Holds ``sigma`` fixed and solves ``N(mu) = Ncond``.  The ``M``
        eigenvalues are computed ONCE via :meth:`_matsubara_number_operator``,
        so both ``N(mu)`` and its analytic derivative ``dN/dmu`` are cheap sums
        over those eigenvalues (:meth:`_number_from_eigs`) -- no per-step matrix
        inversion.  ``N(mu)`` is smooth and monotonically increasing, so the
        root is found with a **safeguarded Newton** iteration (Numerical
        Recipes ``rtsafe``): a Newton step when it stays inside the current
        bracket and makes progress, otherwise a bisection step.  This gets
        Newton's quadratic convergence (a handful of evaluations from the
        previous iteration's mu) while the maintained bracket guarantees
        convergence even if a Newton step misbehaves.

        The initial bracket is the bare eigenvalue range widened by the
        self-energy scale, because Sigma can push spectral weight outside the
        bare band edges.

        Parameters
        ----------
        sigma : ndarray
            Current self-energy (fixed during the mu search).
        beta : float
            Inverse temperature.
        Ncond : float
            Target particle number (already halved for spin-free, matching
            :meth:`RPA._find_mu`).

        Returns
        -------
        float
            Chemical potential mu such that N(mu) = Ncond.
        """
        w = self.H0_eigenvalue

        # Diagonalize the mu-independent part of G^{-1} once; each trial mu is
        # then a cheap sum over eigenvalues (see _matsubara_number_operator).
        lam, ew = self._matsubara_number_operator(sigma, beta)

        def _delta_n(mu, with_deriv=False):
            r = self._number_from_eigs(lam, ew, mu, beta, with_deriv=with_deriv)
            if with_deriv:
                return r[0] - Ncond, r[1]
            return r - Ncond

        # Initial bracket: the bare band range widened by the self-energy scale.
        # This is only a starting guess -- a general (multi-orbital, off-diagonal
        # or non-Hermitian) Sigma can shift the dressed spectral weight, and
        # hence the root, by more than max|Re Sigma|.  So EXPAND the bracket
        # (doubling the pad) until N(mu) changes sign, instead of failing on the
        # first guess.  N(mu) -> 0 as mu -> -inf and -> Nstate as mu -> +inf, so
        # a finite root is always bracketed after enough expansion.
        pad = float(np.abs(sigma.real).max()) + 1.0
        lo = float(w.min()) - pad
        hi = float(w.max()) + pad
        f_lo = _delta_n(lo)
        f_hi = _delta_n(hi)
        for _ in range(60):
            if f_lo == 0.0:
                return lo
            if f_hi == 0.0:
                return hi
            if f_lo * f_hi < 0.0:
                break
            span = hi - lo
            lo -= span
            hi += span
            f_lo = _delta_n(lo)
            f_hi = _delta_n(hi)
        else:
            # 60 doublings span ~1e18 * initial width: a real root cannot be
            # this far out. Fail loudly rather than return garbage. Raise a
            # catchable exception (not sys.exit) so parameter sweeps, pipelines,
            # and notebooks can handle/aggregate the failure instead of having
            # the whole interpreter torn down.
            raise RuntimeError(
                "FLEX._find_mu_dressed: chemical-potential root not bracketed "
                "after expansion to [{}, {}] (N-Ncond = {}, {}). Check the "
                "target filling/Ncond and temperature.".format(
                    lo, hi, f_lo, f_hi))

        # orient the bracket so f(lo) < 0 < f(hi) (N increases with mu)
        if f_lo > 0.0:
            lo, hi = hi, lo

        mu = 0.5 * (lo + hi)
        dmu_old = abs(hi - lo)
        dmu = dmu_old
        f, df = _delta_n(mu, with_deriv=True)

        xtol = 1.0e-12
        for _ in range(100):
            # take bisection when the Newton step leaves the bracket, is not
            # shrinking the interval fast enough, or the derivative is flat;
            # otherwise take Newton.  df -> 0 is the CHARGE-GAP regime: N(mu) is
            # flat across the gap, so a Newton step is unreliable (and f/df ill
            # defined).  The rtsafe product test below already routes tiny df to
            # bisection, but guard df == 0 explicitly so f/df is never formed.
            newton_out = ((mu - hi) * df - f) * ((mu - lo) * df - f) > 0.0
            slow = abs(2.0 * f) > abs(dmu_old * df)
            if newton_out or slow or df == 0.0:
                dmu_old = dmu
                dmu = 0.5 * (hi - lo)
                mu = lo + dmu
            else:
                dmu_old = dmu
                dmu = f / df
                mu = mu - dmu

            if abs(dmu) < xtol:
                break

            f, df = _delta_n(mu, with_deriv=True)
            # keep the sign-change bracket [lo, hi] around the root
            if f < 0.0:
                lo = mu
            else:
                hi = mu

        logger.info("FLEX._find_mu_dressed: mu = {}".format(mu))
        return mu

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
    def _flex_compute_veff_general(self, chi0q_raw, ham_orig):
        """General (paramagnetic full-vertex) counterpart of _flex_compute_veff.

        Mirrors the return signature ``(chi0q_out, v_eff, chi_s, chi_c)`` of
        :meth:`_flex_compute_veff`, but uses the general-path methods that keep
        the full rank-4 orbital vertex (MYO-convention S/C matrices) instead of
        the density-density reduction.

        Parameters
        ----------
        chi0q_raw : ndarray
            Raw rank-6 chi0q ``(nmat, nvol, norb, norb, norb, norb)`` from
            :meth:`_calc_chi0q` after the spin block dimension has been stripped.
        ham_orig : ndarray
            Original interaction Hamiltonian in full spin-orbital space.

        Returns
        -------
        chi0q_out : ndarray
            chi0q for public output, in the RPA ``[a,c,b,d]`` orbital convention
            (rank-6), consistent with the reduced path -- NOT the internal MYO
            layout used for the S/C math.
        v_eff : ndarray
            Effective FLEX interaction ``(nmat, nvol, norb^2, norb^2)``.
        chi_s : ndarray
            Spin susceptibility in the general-path **MYO** orbital convention
            (rank-6 ``(nmat,nvol,m,n,mu,nu)``).
        chi_c : ndarray
            Charge susceptibility in the general-path **MYO** orbital convention
            (rank-6 ``(nmat,nvol,m,n,mu,nu)``).
        """
        chi0q, Us, Uc = self._inflate_chi0q_and_ham_general(chi0q_raw, ham_orig)
        chi_s, chi_c = self._solve_channels_general(chi0q, Us, Uc)
        v_eff = self._calc_veff_general(chi0q, chi_s, chi_c, Us, Uc)
        # Expose chi0q in the RPA [a,c,b,d] convention (consistent with the
        # reduced path) for the public output. The MYO transpose is an internal
        # detail of the S/C math; transposing back (the transpose is an
        # involution) recovers the input convention. chi_s/chi_c remain in the
        # general-path MYO convention (see _solve_channels_general).
        chi0q_out = chi0q.transpose(0, 1, 4, 5, 2, 3)
        return chi0q_out, v_eff, chi_s, chi_c

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

            # Inflate to spin-orbital reduced space.
            # Equivalent to a Kronecker product with I_ns: scatter chi0q onto
            # the spin-block diagonal (same spin block for row and col),
            # leaving off-diagonal spin blocks exactly zero. Bit-identical to
            # np.einsum('lkab,st->lksatb', chi0q, I_ns), but faster.
            chi0q_src = chi0q_raw.reshape(nfreq, nvol, norb, norb)
            chi0q = np.zeros((nfreq, nvol, nd, nd), dtype=chi0q_src.dtype)
            for s in range(ns):
                sl = slice(s * norb, (s + 1) * norb)
                chi0q[..., sl, sl] = chi0q_src

            ham = np.einsum('ksasatbtb->ksatb',
                            ham_orig.reshape(nvol, *(ns, norb) * 4)
                            ).reshape(nvol, nd, nd)

        elif self.spin_mode == "spin-diag":
            # chi0q_raw shape: (nblock=2, nmat, nvol, norb, norb) for reduced
            nblock_s, nfreq, nvol_c, norb1, norb2 = chi0q_raw.shape

            # Inflate per-spin-block chi0q to spin-orbital reduced space.
            # Source block g goes to spin block (g, g) on the diagonal, with
            # off-diagonal spin blocks exactly zero. Bit-identical to
            # np.einsum('glkab,gh->lkgahb', chi0q, I_ns), but faster.
            chi0q = np.zeros((nfreq, nvol, nd, nd), dtype=chi0q_raw.dtype)
            for s in range(ns):
                sl = slice(s * norb, (s + 1) * norb)
                chi0q[..., sl, sl] = chi0q_raw[s]

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
    def _inflate_chi0q_and_ham_general(self, chi0q_raw, ham_orig):
        """Convert chi0q to MYO convention; build full MYO S/C interaction matrices.

        This is the paramagnetic (spin-free) full-vertex general-path analogue
        of :meth:`_inflate_chi0q_and_ham`.  Unlike the reduced/squashed path,
        which collapses the vertex onto the density-density diagonal and works
        in the spin-orbital ``nd = norb * ns`` space, the MYO paramagnetic
        formalism (cond-mat/0407094 Eqs. 5/6) is purely in ORBITAL space and
        keeps the full rank-4 vertex.  Hence the working dimension here is
        ``no = self.norb`` (NOT ``self.nd``): chi0q stays orbital-only and the
        S/C matrices are ``norb^2 x norb^2`` MYO Kanamori blocks.

        This method ALSO converts chi0q from the RPA orbital-pair convention to
        the MYO convention via an orbital-pair transpose.  ``RPA._calc_chi0q``
        labels the rank-6 bare bubble as ``chi0[..., a, c, b, d]`` (rpa.py:
        ``chi0[a,c,b,d] = G[a,b] . G_rev[d,c]``); viewed as a
        ``(norb^2, norb^2)`` matrix with row ``(a,c)`` and column ``(b,d)`` this
        is the MYO susceptibility ``chi0_{mn,mu nu}`` TRANSPOSED (RPA's row pair
        ``(a,c)`` is MYO's COLUMN ``(mu,nu)``; RPA's column pair ``(b,d)`` is
        MYO's ROW ``(m,n)``).  The whole downstream general-path machinery -- the
        MYO S/C matrices, :meth:`_solve_channels_general`,
        :meth:`_calc_veff_general`, and :meth:`_sigma_orbital_contract` -- is
        written in the MYO convention, so we transpose the orbital-pair axes here
        to convert RPA -> MYO.  Verified against the physical brute-force
        pipeline (``tests/flex_bruteforce_ref``) to ~1e-13.

        Parameters
        ----------
        chi0q_raw : ndarray
            Bare susceptibility, shape (nmat,nvol,norb,norb,norb,norb)
            (6-dim: freq, vol, 4 orbital legs) from ``RPA._calc_chi0q`` for the
            spin-free general scheme, in RPA orbital-pair convention
            ``[a,c,b,d]``, after the (size-1) spin-block dimension has been
            stripped in ``solve()``.
        ham_orig : ignored
            Present for signature parity with ``_inflate_chi0q_and_ham``; the
            general path rebuilds the interaction from the real-space
            ``self.ham_info.param_ham`` rather than the pre-inflated
            ``ham_inter_q``.

        Returns
        -------
        chi0q : ndarray
            ``chi0q_raw`` with its orbital-pair axes transposed into MYO
            convention (shape unchanged: 6-dim ``(nmat,nvol,m,n,mu,nu)``).
        Us : ndarray
            MYO spin (S) interaction matrices, shape ``(nvol, norb^2, norb^2)``.
        Uc : ndarray
            MYO charge (C) interaction matrices, shape ``(nvol, norb^2, norb^2)``.

        Notes
        -----
        The S/C matrices are built via ``build_sc_matrices_myo`` from an
        ``inter_k`` dict assembled with ``hwave.sc._build_interaction_k``.  For
        on-site Kanamori interactions the S/C matrices are CONSTANT over q, so
        the k-array ordering used to build ``inter_k`` (a plain linspace grid)
        need not match the FFT q-grid for this v1 path; this is reshaped to
        ``(nvol, norb^2, norb^2)`` purely as ``Nx*Ny*Nz`` independent copies.
        The reshape yields the matrix-per-q form that the downstream
        channel-solver (``_solve_rpa``) consumes as ``ham``.  Because the S/C
        matrices are q-independent constants for on-site Kanamori, they are
        cached across SCF iterations; only the (per-iteration) chi0 transpose is
        recomputed each call.
        """
        logger.debug(">>> FLEX._inflate_chi0q_and_ham_general")

        # RPA's _calc_chi0q labels the rank-6 bare bubble as chi0[..., a, c, b, d]
        # (rpa.py: chi0[a,c,b,d] = G[a,b]·G_rev[d,c]); as a (norb^2, norb^2) matrix
        # (row=(a,c), col=(b,d)) this is the MYO susceptibility χ⁰_{mn,μν} TRANSPOSED
        # (RPA's (a,c) pair is MYO's column (μ,ν); (b,d) is MYO's row (m,n)). The whole
        # downstream general-path machinery (MYO S/C matrices, _solve_channels_general,
        # _calc_veff_general, _sigma_orbital_contract) is written in the MYO convention,
        # so transpose the orbital-pair axes here to convert RPA→MYO. Verified against
        # the physical brute-force pipeline (tests/flex_bruteforce_ref) to ~1e-13.
        chi0q = chi0q_raw.transpose(0, 1, 4, 5, 2, 3)   # (nmat,nvol,a,c,b,d) -> (nmat,nvol,m,n,μ,ν)
        assert chi0q.ndim == 6

        # MYO S/C matrices are q-independent constants for on-site Kanamori, so
        # build them once and cache across SCF iterations. The chi0 transpose
        # above stays OUTSIDE the cache (chi0 changes every iteration).
        cache = getattr(self, "_myo_sc_cache", None)
        if cache is None:
            from hwave.sc import _build_interaction_k
            from hwave.solver._sc_matrices_myo import build_sc_matrices_myo

            # PairLift does not contribute to the particle-hole spin/charge
            # (S/C) vertex: its contribution is S=C=0 (verified against the full
            # 4-index RPA; cf. the same treatment in hwave.sc), so the MYO S/C
            # builder correctly omits it. It is therefore physically inert here,
            # but warn so a configured PairLift term is not silently ignored --
            # matching the Eliashberg-path wording.
            if "PairLift" in self.ham_info.param_ham:
                logger.warning(
                    "PairLift is configured but does not contribute to the S/C "
                    "pairing vertex (S=C=0); it is ignored in the general FLEX "
                    "calculation.")

            # Fail-fast: the general (full-vertex) path builds the MYO S/C
            # matrices on a uniform k-grid that is NOT the FFT q-grid used by
            # _calc_chi0q.  For ON-SITE Kanamori interactions the S/C matrices
            # are q-independent constants so this is exact; but an OFF-SITE
            # interaction entry (irvec != (0,0,0)) makes them genuinely
            # q-dependent on the wrong grid -> silently wrong physics.  v1 of
            # the general path is on-site-only, so reject off-site entries.
            for itype in ("CoulombIntra", "CoulombInter", "Hund",
                          "Exchange", "PairHop", "Ising"):
                if itype in self.ham_info.param_ham:
                    for (irvec, orbvec) in self.ham_info.param_ham[itype]:
                        if tuple(irvec) != (0, 0, 0):
                            raise ValueError(
                                "FLEX calc_scheme='general' (v1) supports "
                                "on-site interactions only; interaction '{}' "
                                "has an off-site entry irvec={}. Off-site "
                                "two-body interactions are not yet supported "
                                "by the general full-vertex path.".format(
                                    itype, tuple(irvec)))

            no = self.norb
            nx, ny, nz = self.lattice.shape

            # Build k-space interactions from the raw real-space param_ham. The
            # k-array ordering is irrelevant for on-site Kanamori terms (constant
            # over q), so a simple uniform grid suffices for v1.
            kx = np.linspace(0, 2.0 * np.pi, nx, endpoint=False)
            ky = np.linspace(0, 2.0 * np.pi, ny, endpoint=False)
            kz = np.linspace(0, 2.0 * np.pi, nz, endpoint=False)
            inter_k = _build_interaction_k(kx, ky, kz,
                                           self.ham_info.param_ham, no)

            # MYO S/C matrices: (nx, ny, nz, norb^2, norb^2).
            Us, Uc = build_sc_matrices_myo(inter_k, no, nx, ny, nz)

            # Reshape to (nvol, norb^2, norb^2) for the downstream channel solver.
            nvol = self.lattice.nvol
            Us = Us.reshape(nvol, no * no, no * no)
            Uc = Uc.reshape(nvol, no * no, no * no)

            self._myo_sc_cache = (Us, Uc)
        else:
            Us, Uc = cache

        return chi0q, Us, Uc

    @do_profile
    def _solve_channels_general(self, chi0q, Us, Uc):
        """Solve the spin and charge RPA channels for the general MYO path.

        In the MYO / Takimoto-Hotta-Ueda paramagnetic full-vertex convention
        (cond-mat/0407094 Eq. 4) the channel susceptibilities are

            chi_s = [I - chi0 . Us]^{-1} . chi0      (spin   channel)
            chi_c = [I + chi0 . Uc]^{-1} . chi0      (charge channel)

        i.e. the spin channel enters with a MINUS sign (it is the Stoner /
        magnetic-instability channel that diverges as chi0.Us -> I) and the
        charge channel with a PLUS sign.  This mirrors the sign discussion in
        :meth:`_build_spin_charge_vertices` for the reduced path.

        The inherited matrix solver :meth:`RPA._solve_rpa` computes
        ``[I + chi0 . ham]^{-1} . chi0``.  To realize the MYO signs we therefore
        pass ``ham_s = -Us`` (so ``+chi0.(-Us) = -chi0.Us``) and
        ``ham_c = +Uc``.  ``_solve_rpa`` accepts the orbital-space ``chi0q`` --
        a 6-dimensional array of shape ``(nmat, nvol, norb, norb, norb, norb)``
        (four orbital legs ``m,n,mu,nu``) -- and the ``(nvol, norb^2, norb^2)``
        interaction matrices directly (it reshapes chi0q to
        ``(nmat, nvol, norb^2, norb^2)`` internally), so ``Us``/``Uc`` from
        :meth:`_inflate_chi0q_and_ham_general` are used as-is.

        Parameters
        ----------
        chi0q : ndarray
            Bare susceptibility, shape ``(nmat, nvol, norb, norb, norb, norb)``
            (six dimensions: frequency, volume, and four orbital legs).
        Us : ndarray
            MYO spin (S) interaction matrices ``(nvol, norb^2, norb^2)``.
        Uc : ndarray
            MYO charge (C) interaction matrices ``(nvol, norb^2, norb^2)``.

        Returns
        -------
        chi_s : ndarray
            Spin-channel RPA susceptibility, same shape as ``chi0q``.
        chi_c : ndarray
            Charge-channel RPA susceptibility, same shape as ``chi0q``.
        """
        logger.debug(">>> FLEX._solve_channels_general")

        # _solve_rpa computes [I + chi0 ham]^-1 chi0; pass -Us / +Uc to obtain
        # chi_s = [I - chi0 Us]^-1 chi0 and chi_c = [I + chi0 Uc]^-1 chi0.
        chi_s = self._solve_rpa(chi0q, -Us)
        chi_c = self._solve_rpa(chi0q, +Uc)
        return chi_s, chi_c

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
        # Inflate channel vertices to spin-orbital reduced space by scattering
        # onto the spin-block diagonal (Kronecker product with I_ns).
        # Bit-identical to np.einsum('kab,st->ksatb', u, I_ns), but faster.
        ham_s = np.zeros((nvol, nd, nd), dtype=u_s.dtype)
        ham_c = np.zeros((nvol, nd, nd), dtype=u_c.dtype)
        for s in range(ns):
            sl = slice(s * norb, (s + 1) * norb)
            ham_s[..., sl, sl] = -u_s
            ham_c[..., sl, sl] = u_c

        return ham_s, ham_c

    @do_profile
    def _calc_veff(self, chi0q, chi_s, chi_c, ham_inflated):
        """Compute effective FLEX interaction V_eff(q, ivn).

        V_eff = W * [3/2*chi_s + 1/2*chi_c - chi0] * W

        chi0 is subtracted exactly once: at zeroth order chi_s = chi_c = chi0,
        so the kernel reduces to chi0 and V_eff = W chi0 W reproduces the
        second-order (SOPT) bubble.  This single subtraction removes the
        double-counted second-order diagram contained in chi_s and chi_c.

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

        # Fluctuation susceptibility.
        # FLEX kernel: 3/2 chi_s + 1/2 chi_c - chi0.  chi0 must be subtracted
        # exactly ONCE: at zeroth order chi_s = chi_c = chi0, so the kernel is
        # 3/2 chi0 + 1/2 chi0 - chi0 = chi0, giving V_eff = W chi0 W ~ U^2 (the
        # second-order bubble).  Subtracting chi0 in *both* terms would cancel
        # the O(U^0) part and lose the leading SOPT self-energy.
        fluct_chi = 1.5 * chi_s_2d + 0.5 * chi_c_2d - chi0q_2d

        # V_eff = W * fluct_chi * W
        # Use batched matmul instead of einsum for better BLAS utilization
        # W @ fluct_chi: broadcast (nvol, nd, nd) @ (nfreq, nvol, nd, nd)
        tmp = np.matmul(ham_2d, fluct_chi)
        # tmp @ W: (nfreq, nvol, nd, nd) @ (nvol, nd, nd)
        v_eff = np.matmul(tmp, ham_2d)

        return v_eff

    @do_profile
    def _calc_veff_general(self, chi0q, chi_s, chi_c, Us, Uc):
        r"""Compute the MYO full-vertex effective interaction V_eff(q, ivn).

        Implements the fluctuation part of the paramagnetic full-vertex
        interaction (Takimoto-Hotta-Ueda / MYO convention, cond-mat/0407094
        Eq. 4)::

            V(q) = 3/2 Us.chi_s.Us + 1/2 Uc.chi_c.Uc
                   - 1/4 (Us+Uc).chi0.(Us+Uc)

        The first-order constant terms (``+3/2 Us - 1/2 Uc``) are intentionally
        EXCLUDED so that, like the reduced path :meth:`_calc_veff`, the
        second-order (SOPT) limit is reproduced: at zeroth order
        ``chi_s = chi_c = chi0`` and (with ``Us = Uc = U``) the kernel reduces
        to ``U.chi0.U``, the leading second-order bubble.

        All three susceptibilities are six-dimensional orbital-space arrays of
        shape ``(nmat, nvol, norb, norb, norb, norb)`` (four orbital legs).
        They are flattened to ``(nmat, nvol, ndx, ndx)`` matrices with
        ``ndx = norb^2`` BEFORE the matrix products, using the same flatten
        convention as the MYO S/C builder: the row index is
        ``(l1 * norb + l2)`` (``idx12``) and the column index is
        ``(l3 * norb + l4)`` (``idx34``), matching ``Us`` / ``Uc`` / ``chi0q``.
        ``Us`` / ``Uc`` are ``(nvol, ndx, ndx)`` and are broadcast over the
        frequency axis.

        Parameters
        ----------
        chi0q : ndarray
            Bare susceptibility, shape ``(nmat, nvol, norb, norb, norb, norb)``.
        chi_s : ndarray
            Spin-channel RPA susceptibility, same shape as ``chi0q``.
        chi_c : ndarray
            Charge-channel RPA susceptibility, same shape as ``chi0q``.
        Us : ndarray
            MYO spin (S) interaction matrices ``(nvol, norb^2, norb^2)``.
        Uc : ndarray
            MYO charge (C) interaction matrices ``(nvol, norb^2, norb^2)``.

        Returns
        -------
        ndarray
            Effective interaction V_eff, shape
            ``(nmat, nvol, norb^2, norb^2)``.
        """
        logger.debug(">>> FLEX._calc_veff_general")

        nmat, nvol = chi0q.shape[0], chi0q.shape[1]
        no = chi0q.shape[2]
        ndx = no * no

        # Flatten the four orbital legs to (row, col) matrices BEFORE matmul.
        chi0_2d = chi0q.reshape(nmat, nvol, ndx, ndx)
        chis_2d = chi_s.reshape(nmat, nvol, ndx, ndx)
        chic_2d = chi_c.reshape(nmat, nvol, ndx, ndx)

        # Broadcast the (nvol, ndx, ndx) interaction matrices over frequency.
        UsB = Us[np.newaxis]            # (1, nvol, ndx, ndx)
        UcB = Uc[np.newaxis]
        Uspc = (Us + Uc)[np.newaxis]

        term_s = 1.5 * (UsB @ chis_2d @ UsB)
        term_c = 0.5 * (UcB @ chic_2d @ UcB)
        term_0 = 0.25 * (Uspc @ chi0_2d @ Uspc)

        v_eff = term_s + term_c - term_0    # (nmat, nvol, ndx, ndx)
        return v_eff

    @do_profile
    def _calc_self_energy(self, green_kw, v_eff, beta):
        """Compute self-energy Sigma(k, iwn) via FFT convolution.

        Sigma(k, iwn) = T/Nk * sum_{q,m} V_eff(q, ivm) * G(k-q, iwn-ivm)

        The convolution is computed efficiently using FFT:
        1. Transform G and V_eff to (r, tau) space
        2. Multiply in (r, tau) space: Sigma(r,tau) = V(r,tau) * G(r,tau)
        3. Transform back to (k, iwn) space

        ``green_kw`` is the FULL physical Green's function (NOT the
        tail-subtracted one used for chi0q): the self-energy convolution is an
        exact cyclic frequency convolution and needs no tau-space tail
        reconstruction.

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

        # The self-energy convolution Sigma(r,tau) = V_eff(r,tau) * G(r,tau)
        # requires V_eff and G to share the same imaginary-time grid, i.e. the
        # bosonic frequency count must equal the fermionic one.  A smaller
        # nfreq would silently leave tau slices at zero (n_common below), so
        # fail loudly instead of returning a wrong self-energy.
        if nfreq != nmat:
            raise ValueError(
                "FLEX self-energy requires a full bosonic frequency grid: "
                "V_eff has nfreq={} but Nmat={}".format(nfreq, nmat))

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

        # If nd_block != nd_v, each Green block only sees its own spin slot
        # of V_eff (spin-free: nd_block = norb, nd_v = norb*ns = nd)
        if nd_block != nd_v:
            # The inflated V_eff carries the (raw) chi0q block g on its spin
            # diagonal slot (g, g) (see _inflate_chi0q_and_ham), so:
            #   spin-diag (nblock == ns): Green block g pairs with slot g
            #   spin-free (nblock == 1):  both slots are identical, use slot 0
            # Sigma stays in block space (nd_block x nd_block); the spin
            # off-diagonal slots of the inflated form are exactly zero.
            norb = self.norb
            sigma_rt = np.zeros((nblock, nmat, nvol, nd_block, nd_block),
                                dtype=np.complex128)
            for g in range(nblock):
                s = g if nblock == self.ns else 0
                sl = slice(s * norb, (s + 1) * norb)
                sigma_rt[g, :n_common] = (
                    v_rt[:n_common, :, sl, sl]
                    * green_rt[g, :n_common]
                )
        else:
            # Sigma(r,tau) = V_eff(r,tau) * G(r,tau)  (element-wise product)
            sigma_rt = np.zeros((nblock, nmat, nvol, nd_v, nd_v),
                                dtype=np.complex128)
            sigma_rt[:, :n_common] = v_rt[:n_common] * green_rt[:, :n_common]

        # --- Transform Sigma back to (k, iwn) space ---
        nd_sig = sigma_rt.shape[-1]

        # Real-space -> k-space
        sigma_kt = FFT.fftn(
            sigma_rt.reshape(nblock, nmat, nx, ny, nz, nd_sig * nd_sig),
            axes=(2, 3, 4)
        ).reshape(nblock, nmat, nvol * nd_sig * nd_sig)

        # Imaginary time -> Matsubara freq (fermionic)
        # Use broadcasting instead of einsum for phase multiplication
        omg_f_inv = np.exp(1j * np.pi * (1.0 / nmat - 1.0) * np.arange(nmat))
        sigma_kw = (FFT.ifft(sigma_kt * omg_f_inv[np.newaxis, :, np.newaxis],
                             axis=1)
                    .reshape(nblock, nmat, nvol, nd_sig, nd_sig) * (1.0 / beta))

        return sigma_kw

    @do_profile
    def _sigma_orbital_contract(self, v_rt, green_rt):
        r"""Per-(r,tau) rank-4 orbital contraction.

        Implements the orbital part of MYO (cond-mat/0407094) Eq. 3::

            Sigma_{mn} = sum_{mu,nu} V_{(mu m),(nu n)} G_{mu nu}

        Parameters
        ----------
        v_rt : ndarray, shape (nfreq, nvol, norb^2, norb^2)
            MYO vertex in (r, tau).  The last two axes flatten the orbital
            pairs (mu, m) [row] and (nu, n) [col], i.e.
            ``v_rt.reshape(..., norb, norb, norb, norb)`` has axes
            (freq, vol, mu, m, nu, n).
        green_rt : ndarray, shape (nblock, nfreq, nvol, norb, norb)
            Green's function in (r, tau), last two axes (mu, nu).

        Returns
        -------
        ndarray, shape (nblock, nfreq, nvol, norb, norb)
            Self-energy in (r, tau), last two axes (m, n).
        """
        nfreq, nvol = v_rt.shape[0], v_rt.shape[1]
        no = green_rt.shape[-1]
        v6 = v_rt.reshape(nfreq, nvol, no, no, no, no)   # (f, r, mu, m, nu, n)
        # Sigma[g,f,r,m,n] = sum_{mu,nu} v6[f,r,mu,m,nu,n] * green_rt[g,f,r,mu,nu]
        return np.einsum('frumvn,gfruv->gfrmn', v6, green_rt)

    @do_profile
    def _calc_self_energy_general(self, green_kw, v_eff, beta):
        r"""Compute Sigma(k, iwn) for the general (full-vertex) FLEX path.

        Implements MYO (cond-mat/0407094) Eq. 3::

            Sigma_{mn}(k, iwn)
              = T/Nk * sum_{q,iv} sum_{mu,nu}
                       V_{(mu m),(nu n)}(q, iv) G_{mu nu}(k-q, iwn-iv)

        The frequency/momentum FFT transport (G -> (r,tau), V -> (r,tau),
        per-(r,tau) combine, then back to (k, iwn)) is IDENTICAL to the reduced
        :meth:`_calc_self_energy`; see that method for the derivation.  The ONLY
        difference is the per-(r,tau) combine step: the reduced path uses an
        element-wise Hadamard ``V * G``, whereas here it is the rank-4 orbital
        contraction :meth:`_sigma_orbital_contract`.

        Because the general path is pure-orbital (``nd_block == nd_v == norb``),
        the spin-inflation branch of the reduced method is not needed and is
        dropped.

        Parameters
        ----------
        green_kw : ndarray, shape (nblock, nmat, nvol, norb, norb)
            Dressed Green's function (orbital axes (mu, nu)).
        v_eff : ndarray, shape (nfreq, nvol, norb^2, norb^2)
            MYO effective interaction; last two axes flatten the orbital pairs
            (mu, m) [row] and (nu, n) [col].
        beta : float
            Inverse temperature.

        Returns
        -------
        ndarray, shape (nblock, nmat, nvol, norb, norb)
            Self-energy (orbital axes (m, n)).
        """
        logger.debug(">>> FLEX._calc_self_energy_general")

        nx, ny, nz = self.lattice.shape
        nvol = self.lattice.nvol
        nmat = self.nmat
        nblock = green_kw.shape[0]
        norb = green_kw.shape[-1]
        ndx = v_eff.shape[-1]   # = norb^2
        nfreq = v_eff.shape[0]

        # Shared with _calc_self_energy: V_eff and G must share the imaginary-
        # time grid (bosonic freq count == fermionic), else tau slices would be
        # silently left at zero.  Fail loudly.
        if nfreq != nmat:
            raise ValueError(
                "FLEX self-energy requires a full bosonic frequency grid: "
                "V_eff has nfreq={} but Nmat={}".format(nfreq, nmat))

        # --- Transform Green's function to (r, tau) space ---
        # (transport identical to _calc_self_energy)
        omg_f = np.exp(-1j * np.pi * (1.0 / nmat - 1.0) * np.arange(nmat))
        green_flat = green_kw.reshape(nblock, nmat, nvol * norb * norb)
        green_kt = (FFT.fft(green_flat, axis=1)
                    * omg_f[np.newaxis, :, np.newaxis]
                    ).reshape(nblock, nmat, nx, ny, nz, norb * norb)
        green_rt = FFT.ifftn(green_kt, axes=(2, 3, 4)
                             ).reshape(nblock, nmat, nvol, norb, norb)

        # --- Transform V_eff to (r, tau) space ---
        # (transport identical to _calc_self_energy)
        omg_b = np.exp(-1j * np.pi * np.arange(nfreq))
        v_flat = v_eff.reshape(nfreq, nvol * ndx * ndx)
        v_qt = (FFT.fft(v_flat, axis=0)
                * omg_b[:, np.newaxis]
                ).reshape(nfreq, nx, ny, nz, ndx * ndx)
        v_rt = FFT.ifftn(v_qt, axes=(1, 2, 3)).reshape(nfreq, nvol, ndx, ndx)

        # --- Compute Sigma(r, tau): rank-4 orbital contraction (the only new
        # bug-prone step; the rest is reused transport). ---
        sigma_rt = self._sigma_orbital_contract(v_rt, green_rt)

        # --- Transform Sigma back to (k, iwn) space ---
        # (transport identical to _calc_self_energy)
        sigma_kt = FFT.fftn(
            sigma_rt.reshape(nblock, nmat, nx, ny, nz, norb * norb),
            axes=(2, 3, 4)
        ).reshape(nblock, nmat, nvol * norb * norb)

        omg_f_inv = np.exp(1j * np.pi * (1.0 / nmat - 1.0) * np.arange(nmat))
        sigma_kw = (FFT.ifft(sigma_kt * omg_f_inv[np.newaxis, :, np.newaxis],
                             axis=1)
                    .reshape(nblock, nmat, nvol, norb, norb) * (1.0 / beta))

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

        # FLEX never applies the matsubara_frequency filter to its outputs
        # (unlike RPA.solve): every stored frequency axis is the full nmat
        # grid, so the metadata must describe that grid -- writing the
        # restricted user option (self.freq_index) would make the strict
        # hwave_sc loader reject a valid file.
        full_freq_index = np.arange(self.nmat)
        if len(self.freq_index) < self.nmat:
            logger.warning(
                "matsubara_frequency is restricted but FLEX outputs always "
                "hold the full frequency grid; the option is ignored.")

        # Save physical observables (particle number, spin, chemical potential)
        # as a plain-text energy file, mirroring UHFr/UHFk.  In fixed-mu mode
        # the NCond line is the mu-N single point for the given mu.
        if "energy" in info_outputfile:
            physics = green_info.get("physics", getattr(self, "physics", None))
            if physics is None:
                logger.warning("save_results: no physics data to write to "
                               "'{}'".format(info_outputfile["energy"]))
            else:
                file_name = os.path.join(path_to_output,
                                         info_outputfile["energy"])
                with open(file_name, "w") as fw:
                    fw.write("NCond = {}\n".format(physics["NCond"]))
                    fw.write("Sz = {}\n".format(physics["Sz"]))
                    fw.write("ChemicalPotential = {}\n".format(physics["mu"]))
                logger.info("save_results: save energy in file {}".format(
                    file_name))

        # Save chi0q
        if "chi0q" in info_outputfile:
            file_name = os.path.join(path_to_output, info_outputfile["chi0q"])
            np.savez(file_name,
                     chi0q=green_info["chi0q"],
                     freq_index=full_freq_index,
                     # full grid size: lets consumers locate the zero bosonic
                     # frequency (index nmat//2) unambiguously
                     nmat=self.nmat,
                     wavevector_unit=self.kvec,
                     wavevector_index=self.wavenum_table,
                     # FLEX chi0q comes from the same spin-block-ordered RPA
                     # internals; tag it so the chi0q consumers (RPA read_chi0q,
                     # hwave_sc _load_chi0q) accept the file in SO mode.
                     index_convention="spin_block")
            logger.info("save_results: save chi0q in file {}".format(file_name))

        # Save susceptibilities (spin and charge channels separately for
        # Eliashberg). Tag the orbital convention so the downstream Eliashberg
        # consumer (_load_flex_susceptibilities / _compute_vertices_flex) pairs
        # them with the matching S/C matrices: the general (full-vertex) path
        # produces MYO-convention susceptibilities, the reduced path Kuroki.
        common_meta = dict(freq_index=full_freq_index,
                           nmat=self.nmat,
                           wavevector_unit=self.kvec,
                           wavevector_index=self.wavenum_table,
                           chi_convention=("myo" if self._flex_general
                                           else "kuroki"))

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
                     freq_index=full_freq_index,
                     wavevector_unit=self.kvec,
                     wavevector_index=self.wavenum_table)
            logger.info("save_results: save sigma in file {}".format(file_name))

        # Save Green's function
        if "green" in info_outputfile:
            file_name = os.path.join(path_to_output, info_outputfile["green"])
            np.savez(file_name,
                     green=green_info.get("green"),
                     freq_index=full_freq_index,
                     wavevector_unit=self.kvec,
                     wavevector_index=self.wavenum_table)
            logger.info("save_results: save green in file {}".format(file_name))
