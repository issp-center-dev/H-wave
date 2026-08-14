"""INTERNAL module -- the shared non-interacting Green's function builder.

This module is NOT public API. It exists so ``RPA._calc_green`` and
sc.py's bond (Eliashberg) Green carrier can construct ``G0(k, iw_n)`` from
the same eigen-decomposition of ``H0(k)`` through one code path instead of
two hand-duplicated ones. Callers are solver classes and the bubble kernel
(``bubble.py``); the public entry points stay the documented solver
classes and ``bond_channels`` functions.

The binding contract is
``docs/superpowers/specs/2026-08-14-unified-bubble-kernel-design.md``
(section "Green data flow"). See that document for the rationale; this
module transcribes the numerics of ``RPA._calc_green`` (rpa.py, around
line 2877 at series start) without altering them.
"""

import math

import numpy as np

from . import backend as _bk


def build_green(eigenvalues, eigenvectors, mu, beta, nmat, coeff_tail):
    """Build the non-interacting Green's function from an eigen-decomposition.

    Parameters
    ----------
    eigenvalues : ndarray, shape (nvol, p) or (nblock, nvol, p)
        Real eigenvalues of ``H0(k)``, FLATTENED over k already (the
        caller reshapes any ``(Nx, Ny, Nz, p)`` spatial layout to
        ``(nvol, p)`` first, in C order).
    eigenvectors : ndarray, shape (nvol, p, p) or (nblock, nvol, p, p)
        The matching eigenvectors (``eigenvectors[..., k, :, i]`` is the
        i-th eigenvector at k-point k), same flattening as ``eigenvalues``.
    mu : float
        Chemical potential.
    beta : float
        Inverse temperature; must be > 0.
    nmat : int
        Number of (fermionic) Matsubara frequencies; must be a positive
        even integer.
    coeff_tail : float
        High-frequency tail coefficient (see below); must be a finite
        real number. ``0`` disables the tail correction.

    Returns
    -------
    full_kw : ndarray, complex128, shape (nblock, nmat, nvol, p, p)
        The bare Green's function G(k, iw_n).
    deflated_kw : ndarray, complex128, shape (nblock, nmat, nvol, p, p)
        ``full_kw`` with the ``coeff_tail / (iw_n)`` term subtracted in
        frequency space (faster tau-space decay for FFT use). When
        ``coeff_tail == 0`` this is the SAME object as ``full_kw`` (no
        copy) -- deliberate, per the spec's aliasing contract.
    green0_tail : ndarray, complex128, shape (nblock, nmat, nvol, p, p), or None
        The analytic imaginary-time tail contribution that recovers
        ``full_kw`` from ``deflated_kw`` (constant over the frequency
        axis). ``None`` when ``coeff_tail == 0``.

    Notes
    -----
    Centered fermionic Matsubara grid: ``iw_n = 1j*(2n+1-nmat)*pi/beta``.

    ``deflated[g,n,k] = V [1/(iw_n - (e-mu)) - coeff_tail/(iw_n)] V dagger``

    ``full[g,n,k] = deflated[g,n,k] + coeff_tail * (V V dagger)[g,k] / (iw_n)``

    ``green0_tail = V V dagger * coeff_tail * beta / 2`` (constant over n).

    The leading block axis is inserted (length 1) when the inputs carry
    none, so the output is always canonical ``(nblock, nmat, nvol, p, p)``.
    Accepted input shapes are FLATTENED ONLY: ``ndim`` must pair up as
    ``(2, 3)`` (no block axis) or ``(3, 4)`` (with a block axis); anything
    else -- including sc.py's unflattened ``(Nx, Ny, Nz, p)`` /
    ``(Nx, Ny, Nz, p, p)`` layout -- raises ``ValueError``.
    """
    beta, nmat, coeff_tail = _validate_scalars(beta, nmat, coeff_tail)
    ev, V, nblock, nvol, p = _validate_and_canonicalize_shapes(
        eigenvalues, eigenvectors)

    xp = _bk.array_module_of(ev)

    # Centered fermionic Matsubara grid.
    iomega = (xp.arange(nmat) * 2 + 1 - nmat) * np.pi / beta
    wn = 1j * iomega[np.newaxis, :, np.newaxis, np.newaxis]        # (1,nmat,1,1)
    ek = (ev - mu)[:, np.newaxis, :, :]                            # (nblock,1,nvol,p)

    aa = coeff_tail
    g_deflated = 1.0 / (wn - ek) - aa / wn                         # (nblock,nmat,nvol,p)

    # G_ab(k,iw_n) = sum_j V_{a,j} V*_{b,j} * g_j = (V * g) @ V dagger
    V_conj_t = xp.conj(V).swapaxes(-2, -1)                         # (nblock,nvol,p,p)
    Vg = V[:, np.newaxis, :, :, :] * g_deflated[:, :, :, np.newaxis, :]
    deflated_kw = (Vg @ V_conj_t[:, np.newaxis, :, :, :]).astype(
        np.complex128, copy=False)

    if aa == 0.0:
        return deflated_kw, deflated_kw, None

    # VV dagger built ONCE at (nblock,nvol,p,p) -- never materialized at
    # (nblock,nmat,nvol,p,p) before the division by iw_n; the division and
    # the addition back onto `deflated_kw` broadcast the frequency axis in
    # a single assembly step.
    VVt = V @ V_conj_t                                             # (nblock,nvol,p,p)
    wn5 = wn[:, :, :, :, np.newaxis]                               # (1,nmat,1,1,1)
    full_kw = (deflated_kw + aa * VVt[:, np.newaxis, :, :, :] / wn5).astype(
        np.complex128, copy=False)

    green0_tail = (VVt[:, np.newaxis, :, :, :]
                   * xp.ones((1, nmat, 1, 1, 1)) * aa * 0.5 * beta).astype(
        np.complex128, copy=False)

    return full_kw, deflated_kw, green0_tail


def _validate_scalars(beta, nmat, coeff_tail):
    """Validate ``beta``, ``nmat``, ``coeff_tail``; return normalized
    ``(float, int, float)``. All failures are ``ValueError`` (survives
    ``python -O``, unlike a bare ``assert``)."""
    if isinstance(beta, bool) or not isinstance(
            beta, (int, float, np.integer, np.floating)):
        raise ValueError(
            "build_green: beta must be a real number, got {!r}".format(beta))
    beta = float(beta)
    if not (beta > 0.0):
        raise ValueError(
            "build_green: beta must be > 0, got {}".format(beta))

    if isinstance(nmat, bool) or not isinstance(nmat, (int, np.integer)):
        raise ValueError(
            "build_green: nmat must be an integer, got {!r}".format(nmat))
    nmat = int(nmat)
    if nmat <= 0 or nmat % 2 != 0:
        raise ValueError(
            "build_green: nmat must be a positive even integer, got {}".format(nmat))

    if isinstance(coeff_tail, bool) or not isinstance(
            coeff_tail, (int, float, np.integer, np.floating)):
        raise ValueError(
            "build_green: coeff_tail must be a real number, got {!r}".format(
                coeff_tail))
    coeff_tail = float(coeff_tail)
    if not math.isfinite(coeff_tail):
        raise ValueError(
            "build_green: coeff_tail must be finite, got {}".format(coeff_tail))

    return beta, nmat, coeff_tail


def _validate_and_canonicalize_shapes(eigenvalues, eigenvectors):
    """Validate the FLATTENED-only shape contract and insert a leading
    length-1 block axis when the inputs carry none. Returns
    ``(eigenvalues, eigenvectors, nblock, nvol, p)`` with both arrays at
    ndim 3 / 4 respectively."""
    ev = eigenvalues
    V = eigenvectors

    if ev.ndim == 2:
        expected_v_ndim = 3
    elif ev.ndim == 3:
        expected_v_ndim = 4
    else:
        raise ValueError(
            "build_green: eigenvalues must be flattened over k, ndim in "
            "{{2, 3}} ((nvol, p) or (nblock, nvol, p)); got ndim={} "
            "shape={} -- unflattened (e.g. sc.py's (Nx, Ny, Nz, p)) "
            "eigenvalues must be reshaped by the caller first".format(
                ev.ndim, ev.shape))

    if V.ndim != expected_v_ndim:
        raise ValueError(
            "build_green: eigenvectors ndim ({}) does not pair with "
            "eigenvalues ndim ({}); expected eigenvectors ndim {} for "
            "eigenvalues shape {}".format(
                V.ndim, ev.ndim, expected_v_ndim, ev.shape))

    if ev.ndim == 2:
        ev = ev[np.newaxis, :, :]
        V = V[np.newaxis, :, :, :]

    nblock, nvol, p = ev.shape

    if p <= 0 or nvol <= 0 or nblock <= 0:
        raise ValueError(
            "build_green: nblock, nvol, p must all be positive, got "
            "nblock={}, nvol={}, p={}".format(nblock, nvol, p))

    if V.shape != (nblock, nvol, p, p):
        raise ValueError(
            "build_green: eigenvectors shape {} does not match eigenvalues "
            "-- expected {} (nblock={}, nvol={}, p={} from eigenvalues "
            "shape {})".format(
                V.shape, (nblock, nvol, p, p), nblock, nvol, p, ev.shape))

    return ev, V, nblock, nvol, p
