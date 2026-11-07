"""The one-body band of hwave_sc: H(k), eigenpairs, the chemical potential and the mean-field band.

Moved verbatim from hwave.sc (dev/sc-split-plan.md, step 1); the shell
re-exports every name here under its old spelling.
"""
from __future__ import annotations

import logging
from typing import NamedTuple
import numpy as np
from scipy.optimize import brentq

logger = logging.getLogger("hwave_sc")


# Absolute tolerance (energy units) of the spin-independence test of the
# ASSEMBLED one-body band H0(k) -- trans_mod, or T + Sigma_HF[green_init], plus
# coeff_extern * Extern (issue #273). Its up-up and down-down blocks may
# differ, and its spin-mixing blocks may be nonzero, by at most this much.
# It only absorbs convergence round-off of a paramagnetic UHF solution; it
# is absolute so that a large common energy offset cannot hide a finite
# splitting (a relative test would).
SPIN_FREE_ATOL = 1.0e-10


# Interaction/Hamiltonian arrays of the RPA solver that _calc_epsilon_k can
# modify in place (the trans_mod branch adds Extern into the green_info
# array; the spin-orbital green_init branch of _calc_trans_mod used to add its
# mean field into a view of ham_trans_r, before issue #278). The check pass
# restores all of them.
_EPS_K_HAM_ARRAYS = ("ham_trans_r", "ham_trans_q", "ham_extern_r",
                     "ham_extern_q")


HERMITIAN_RTOL = 1.0e-6


# Hermiticity of the assembled band: max |H(k) - H(k)^dagger| above
# max(HERMITIAN_ATOL, HERMITIAN_RTOL * max |H|) is refused; between
# HERMITIAN_ATOL and that bound the Hermitian part is used, with a warning.
HERMITIAN_ATOL = 1.0e-10


_FLEX_BAND_REASON = ("the susceptibilities and the Green function are read "
                     "from the FLEX output, which already includes the "
                     "one-body band of that run")


# [file.input] keys that replace the bare Transfer band (RPA read_init).
_ONE_BODY_FILE_KEYS = ("trans_mod", "green_init")


class _MeanFieldBand(NamedTuple):
    """The one-body band of an hwave_sc run whose input carries
    ``trans_mod`` / ``green_init`` / an active ``Extern`` (issue #273).

    Attributes
    ----------
    epsilon_k : ndarray, shape (norb, norb, Nx, Ny, Nz)
        The RPA solver's assembled H0(k), in sc.py's layout.
    solver : RPA or None
        The RPA solver that assembled it (``_calc_epsilon_k`` already run),
        reused for chi0q in ``chi0q_mode = "calc"``; ``None`` once released.
    sources : tuple of str
        Which inputs define the band (for messages).
    same_as_bare : bool
        True when the assembled H0(k) is bit-for-bit the RPA reader's bare
        Transfer band (e.g. a trans_mod written for U = 0). The caller then
        keeps its own Transfer construction, so such a run is bit-identical
        to the run without the extra input.
    """
    epsilon_k: np.ndarray
    solver: object
    sources: tuple
    same_as_bare: bool = False


def _one_body_band_sources(input_dict):
    """Inputs that make the one-body band differ from the bare Transfer.

    ``trans_mod`` / ``green_init`` in ``[file.input]`` and an ``Extern``
    table in ``[file.input.interaction]``, whatever ``coeff_extern`` is: the
    RPA reader turns a declared Extern into a two-block H0 even at a zero
    coefficient, so that band has to go through the same assembly and
    check (a zero field then simply yields the Transfer band). Keys are
    matched case-insensitively.
    """
    from requests.structures import CaseInsensitiveDict
    fin = (input_dict.get("file", {}) or {}).get("input", {}) or {}
    lowered = {str(k).lower() for k in fin}
    sources = [key for key in _ONE_BODY_FILE_KEYS if key in lowered]
    inter = CaseInsensitiveDict(fin.get("interaction", {}) or {})
    if "Extern" in inter:
        sources.append("Extern")
    return tuple(sources)


def _warn_one_body_band_ignored(input_dict, where, reason=_FLEX_BAND_REASON):
    """WARNING for the paths whose chi and G do not come from hwave_sc's own
    band: the FLEX paths (the band of the FLEX run is already inside its
    outputs) and the bond path with an external ``bond_green``. A
    ``trans_mod`` / ``green_init`` / ``Extern`` in the hwave_sc input is not
    applied there (issue #273)."""
    sources = _one_body_band_sources(input_dict)
    if sources:
        logger.warning("%s is ignored by hwave_sc for %s: %s.",
                       " / ".join(sources), where, reason)


def _spin_dependence(arr, norb):
    """``(splitting, mixing)`` of a ``(..., 2*norb, 2*norb)`` spin-block
    array: max |up-up - down-down| and max |spin-mixing block|."""
    a = np.asarray(arr).reshape(-1, 2, norb, 2, norb)
    split = float(np.max(np.abs(a[:, 0, :, 0, :] - a[:, 1, :, 1, :])))
    mix = float(max(np.max(np.abs(a[:, 0, :, 1, :])),
                    np.max(np.abs(a[:, 1, :, 0, :]))))
    return split, mix


def _refuse_spin_dependent(what, detail, tol=SPIN_FREE_ATOL):
    raise ValueError(
        "hwave_sc: the one-body band defined by {} is spin-dependent ({}; "
        "[eliashberg] spin_free_tol = {:g}). The singlet/triplet pairing "
        "kernel of hwave_sc carries no spin index. If the field is meant to "
        "be paramagnetic, tighten the UHF convergence ([mode.param] EPS), "
        "start UHF from a spin-symmetric initial state, or raise "
        "[eliashberg] spin_free_tol (the band is then spin-averaged). A "
        "genuinely magnetic (spin-diagonal with unequal spins, or "
        "spin-mixing) mean field, or a finite Extern field: run it with "
        "[mode] enable_spin_orbital = true (see the spin-polarized section "
        "of the Eliashberg tutorial).".format(what, detail, tol))


class _AssembledH0Ready(Exception):
    """Internal sentinel: stops the check pass of :func:`_checked_band`
    once H0(k) is assembled, before RPA diagonalizes it."""


def _checked_band(solver, green_info, sources, spin_free_tol=SPIN_FREE_ATOL,
                  allow_average=True, allow_hermitize=True):
    """Assemble the one-body H0(k) ONCE, check it, and return its spin-free
    block ``(nvol, norb, norb)`` -- the band that is then used everywhere.

    The assembly is RPA's own ``_calc_epsilon_k`` (trans_mod, or T +
    Sigma_HF[green_init], plus coeff_extern * Extern), run with the solver's
    ``_keep_spinful`` switch so the two spin blocks stay uncollapsed (RPA's
    own collapse is a RELATIVE ``np.allclose`` test, which a large common
    offset defeats), and stopped before the eigendecomposition. The pass is
    transactional: ``green_info`` is copied and every Hamiltonian array the
    assembly can modify is restored afterwards.

    The checked H0 must be finite and Hermitian (:data:`HERMITIAN_ATOL`).
    Its up-up and down-down blocks are compared and its spin-mixing blocks
    measured (absolute values):

    * both within :data:`SPIN_FREE_ATOL`: the up-up block is returned
      (Extern included) -- the strict default;
    * either above ``spin_free_tol`` (``[eliashberg] spin_free_tol``,
      default :data:`SPIN_FREE_ATOL`): ``ValueError``;
    * otherwise (between the two): the spin-averaged block
      ``(H_up + H_down)/2`` is returned, the mixing blocks are discarded,
      and a WARNING reports both discarded magnitudes.
    """
    ham = solver.ham_info
    snapshot = {}
    for name in _EPS_K_HAM_ARRAYS:
        val = getattr(ham, name, None)
        snapshot[name] = None if val is None else np.array(val, copy=True)
    gi = green_info.copy()
    for key in ("trans_mod", "green_init"):
        if key in gi:
            gi[key] = np.array(gi[key], copy=True)

    def _stop(*args, **kwargs):
        raise _AssembledH0Ready()

    solver._keep_spinful = True
    solver._find_block_diagonal = _stop      # instance override, removed below
    try:
        # non-finite inputs are reported explicitly below
        with np.errstate(invalid="ignore", over="ignore"):
            solver._calc_epsilon_k(gi)
    except _AssembledH0Ready:
        pass
    finally:
        del solver._find_block_diagonal
        solver._keep_spinful = False
        for name, val in snapshot.items():
            setattr(ham, name, val)
    H0 = np.array(solver.H0_k, dtype=np.complex128, copy=True)
    norb = solver.norb
    what = " / ".join(sources)

    if not np.all(np.isfinite(H0)):
        raise ValueError(
            "hwave_sc: the one-body band assembled from {} has non-finite "
            "(NaN or Inf) entries. Check [mode.param] coeff_extern, the "
            "Extern file and the trans_mod / green_init file.".format(what))
    herm = float(np.max(np.abs(H0 - np.conj(np.swapaxes(H0, -1, -2)))))
    herm_bound = max(HERMITIAN_ATOL,
                     HERMITIAN_RTOL * float(np.max(np.abs(H0))))
    if herm > herm_bound:
        raise ValueError(
            "hwave_sc: the one-body band assembled from {} is not Hermitian "
            "(max |H(k) - H(k)^dagger| = {:.3e} > {:.3e} = max({:g}, {:g} * "
            "max |H|)). Check the trans_mod / green_init file and the "
            "Transfer / Extern tables (each R entry needs its conjugate "
            "partner at -R).".format(what, herm, herm_bound, HERMITIAN_ATOL,
                                     HERMITIAN_RTOL))
    if herm > HERMITIAN_ATOL and not allow_hermitize:
        raise ValueError(
            "hwave_sc: the one-body band assembled from {} is not exactly "
            "Hermitian (max |H(k) - H(k)^dagger| = {:.3e} > {:g}). With "
            "chi0q_mode = 'load' the band is not replaced by its Hermitian "
            "part: the loaded chi0q was computed by the RPA run on the raw "
            "band, so mu and G would not match it. Use chi0q_mode = 'calc', "
            "or fix the band (each R entry of trans_mod / Transfer / Extern "
            "needs its conjugate partner at -R).".format(
                what, herm, HERMITIAN_ATOL))
    if herm > HERMITIAN_ATOL:
        # round-off-level asymmetry (e.g. +-R rows rounded independently)
        logger.warning(
            "hwave_sc: the one-body band assembled from %s is not exactly "
            "Hermitian (max |H(k) - H(k)^dagger| = %.3e, within %.3e); using "
            "its Hermitian part (H + H^dagger)/2.", what, herm, herm_bound)
        H0 = 0.5 * (H0 + np.conj(np.swapaxes(H0, -1, -2)))

    if H0.shape[0] == 2:                       # Extern on the Transfer band
        up, dn = H0[0], H0[1]
        split, mix = float(np.max(np.abs(up - dn))), 0.0
    elif H0.shape[-1] == 2 * norb:             # (1, nvol, 2norb, 2norb)
        split, mix = _spin_dependence(H0[0], norb)
        blk = H0[0].reshape(-1, 2, norb, 2, norb)
        up, dn = blk[:, 0, :, 0, :], blk[:, 1, :, 1, :]
    else:                                      # already one spin-free block
        up = dn = H0[0]
        split, mix = 0.0, 0.0
    if not (np.isfinite(split) and np.isfinite(mix)):
        raise ValueError(
            "hwave_sc: the spin comparison of the band assembled from {} is "
            "not finite. Check [mode.param] coeff_extern, the Extern file "
            "and the trans_mod / green_init file.".format(what))
    detail = ("assembled H0(k): max |H_up - H_down| = {:.3e}, max "
              "|spin-mixing| = {:.3e}".format(split, mix))
    if split > spin_free_tol or mix > spin_free_tol:
        _refuse_spin_dependent(what, detail, spin_free_tol)
    if not allow_average and (split > SPIN_FREE_ATOL or mix > SPIN_FREE_ATOL):
        raise ValueError(
            "hwave_sc: the one-body band from {} is spin-dependent beyond "
            "SPIN_FREE_ATOL = {:g} ({}). chi0q_mode = 'load' reads a chi0q "
            "that an external RPA run computed on its own band, which cannot "
            "be the spin-averaged band; the spin-averaged band "
            "([eliashberg] spin_free_tol) is available only when hwave_sc "
            "computes chi0q itself (chi0q_mode = 'calc', or the bond_channels "
            "path). Use chi0q_mode = 'calc', or regenerate the chi0q file "
            "from a spin-symmetric band (tighter UHF convergence or a "
            "spin-symmetric UHF start).".format(what, SPIN_FREE_ATOL, detail))
    if split > SPIN_FREE_ATOL or mix > SPIN_FREE_ATOL:
        # within the user's spin_free_tol: spin-average, drop the mixing
        logger.warning(
            "hwave_sc: the one-body band from %s is spin-dependent beyond "
            "SPIN_FREE_ATOL = %g but within [eliashberg] spin_free_tol = %g "
            "(%s). Using the spin-averaged band (H_up + H_down)/2 and "
            "DISCARDING the spin-mixing blocks: discarded max |H_up - H_down|"
            "/2 = %.3e, discarded max |spin-mixing| = %.3e.",
            what, SPIN_FREE_ATOL, spin_free_tol, detail, split / 2.0, mix)
        return np.ascontiguousarray(0.5 * (up + dn))
    return np.ascontiguousarray(up)


def _check_green_init_spin_free(g, norb, spin_free_tol):
    """Refuse a spin-dependent green_init Green function.

    ``g``: the green_init array as RPA read it, ``(nvol, 2*norb, 2*norb)``
    in spin-block order. A magnetic G describes a magnetic state whatever
    band its mean field happens to give for the declared interactions (with
    no interaction acting on the moment, for instance, T + Sigma_HF[G] is
    spin-free although G is not), so G itself must be spin-independent in
    every cell: up-up equal to down-down and no spin-off-diagonal block,
    within max(SPIN_FREE_ATOL, spin_free_tol) on the dimensionless G
    entries. The assembled band is checked separately (:func:`_checked_band`).
    """
    g = np.asarray(g)
    if not np.all(np.isfinite(g)):
        raise ValueError("hwave_sc: the green_init Green function has "
                         "non-finite (NaN or Inf) entries.")
    split, mix = _spin_dependence(g, norb)
    tol = max(SPIN_FREE_ATOL, spin_free_tol)
    if split > tol or mix > tol:
        raise ValueError(
            "hwave_sc: green_init describes a spin-dependent state (max "
            "|G_up - G_down| = {:.3e}, max |spin-mixing G| = {:.3e}; "
            "tolerance max(1e-10, [eliashberg] spin_free_tol) = {:g} on the "
            "dimensionless Green-function entries). The singlet/triplet "
            "pairing kernel of hwave_sc carries no spin index; for a magnetic "
            "state run it with [mode] enable_spin_orbital = true (see the "
            "spin-polarized section of the Eliashberg tutorial). For a "
            "paramagnetic state, the trans_mod written by UHFk is the "
            "consistent input.".format(split, mix, tol))


def _build_hamiltonian_k(kx_array, ky_array, kz_array, hr, norb):
    """Build Hamiltonian epsilon(k) from real-space transfer integrals.

    Parameters
    ----------
    kx_array, ky_array, kz_array : ndarray
        k-point arrays.
    hr : dict
        Transfer parameters from wan90.read_w90().
    norb : int
        Number of orbitals.

    Returns
    -------
    epsilon_k : ndarray
        Hamiltonian in k-space, shape (norb, norb, Nx, Ny, Nz).
    """
    Nx, Ny, Nz = len(kx_array), len(ky_array), len(kz_array)
    epsilon_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
    kx_mesh, ky_mesh, kz_mesh = np.meshgrid(
        kx_array, ky_array, kz_array, indexing='ij'
    )
    # Solver-core convention (rpa.py _make_ham_trans, issue #133:
    # ifftn * nvol == e^{+ikR}, the documented Wannier90-style sign shared
    # with UHFk): epsilon[a,b](k) = sum_R t_R[a,b] e^{+ikR}. This keeps
    # sc-built quantities element-wise consistent with arrays loaded from
    # FLEX/RPA output files, which carry the same k labeling since the
    # #133 alignment. (The [orb1,orb2] placement is unchanged; only the
    # Fourier sign moved with the solver core.)
    for (irvec, orbvec), value in hr.items():
        orb1, orb2 = orbvec
        Rx, Ry, Rz = irvec
        epsilon_k[orb1, orb2, :, :, :] += value * np.exp(
            +1j * (kx_mesh * Rx + ky_mesh * Ry + kz_mesh * Rz)
        )
    return epsilon_k


def _calc_eigenvalues(epsilon_k):
    """Diagonalize epsilon(k) at every k-point.

    Parameters
    ----------
    epsilon_k : ndarray
        Hamiltonian, shape (norb, norb, Nx, Ny, Nz).

    Returns
    -------
    eigenvalues : ndarray
        Shape (Nx, Ny, Nz, norb).
    eigenvectors : ndarray
        Shape (Nx, Ny, Nz, norb, norb).
    """
    norb = epsilon_k.shape[0]
    Nx, Ny, Nz = epsilon_k.shape[2], epsilon_k.shape[3], epsilon_k.shape[4]
    eigenvalues = np.zeros((Nx, Ny, Nz, norb))
    eigenvectors = np.zeros((Nx, Ny, Nz, norb, norb), dtype=complex)

    # Batch diagonalization: transpose to (Nx, Ny, Nz, norb, norb) for vectorized eigh
    eps_batch = epsilon_k.transpose(2, 3, 4, 0, 1)  # (Nx, Ny, Nz, norb, norb)
    eigenvalues, eigenvectors = np.linalg.eigh(eps_batch)

    return eigenvalues, eigenvectors


def _determine_mu(eigenvalues, beta, n_target, norb):
    """Determine chemical potential, matching RPA's search to round-off.

    The filling convention follows the reference code: n_target is the
    number of electrons per orbital per spin channel. For example,
    n_target=0.75 means 3/4 filling per spin.

    Historically this used a plain ``scipy.optimize.bisect`` with default
    tolerances, whose particle-number residual sat well above round-off and
    disagreed with ``RPA._find_mu`` (#185).  It now uses the SAME machinery
    as RPA: ``brentq`` with ``xtol=1e-14`` bracketed on the eigenvalue span
    (falling back to the wide bracket only when that span is not
    sign-bracketed), followed by the transactional Newton polish
    ``rpa._polish_mu_root``, and it reuses ``rpa._masked_fermi_delta_n`` so
    the overflow-guarded Fermi arithmetic (and the particle-number residual
    it drives to round-off) is bit-for-bit identical to RPA's.  The former
    overflow cutoff of 100 becomes the ``ene_cutoff`` of the masked helper.

    Parameters
    ----------
    eigenvalues : ndarray
        Shape (Nx, Ny, Nz, norb).
    beta : float
        Inverse temperature.
    n_target : float
        Target filling per orbital per spin (e.g. 0.75 for 3/4 filling).
    norb : int
        Number of orbitals.

    Returns
    -------
    mu : float
        Chemical potential.
    """
    from hwave.solver import rpa as _rpa

    Nx, Ny, Nz = eigenvalues.shape[:3]
    nvol = Nx * Ny * Nz
    T = 1.0 / beta
    ene_cutoff = 100.0
    w = np.asarray(eigenvalues)

    # RPA counts an UNNORMALIZED sum over all (k, band) entries of one spin
    # block; sc.py's per-nvol target is n_target * norb, so the matching
    # unnormalized target is n_target * norb * nvol.  Same eigenvalues, same
    # physical per-spin filling as RPA._find_mu.
    target = n_target * norb * nvol

    def _delta_n_and_deriv(mu):
        return _rpa._masked_fermi_delta_n(w, T, mu, target, ene_cutoff)

    def _calc_delta_n(mu):
        return _delta_n_and_deriv(mu)[0]

    ev = np.sort(w.flatten())
    lo, hi = float(ev[0]), float(ev[-1])

    bracketed = (_calc_delta_n(lo) * _calc_delta_n(hi)) < 0.0
    if bracketed:
        mu = brentq(_calc_delta_n, lo, hi, xtol=1e-14)
    else:
        # fall back to the historical wide bracket
        mu = brentq(_calc_delta_n, lo - 10.0, hi + 10.0, xtol=1e-14)

    # Transactional Newton polish (#160): drives the particle-number residual
    # toward round-off; never returns a worse root than brentq's.
    mu = _rpa._polish_mu_root(_delta_n_and_deriv, mu,
                              bracket=(lo, hi) if bracketed else None)
    return float(mu)
