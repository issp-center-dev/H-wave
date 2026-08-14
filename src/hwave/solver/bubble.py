"""INTERNAL module -- shared particle-hole bubble orbital contractions and
the dense uniform-Matsubara-grid transport.

This module is NOT public API; the supported user-facing surface stays the
solver classes (``RPA``, ``FLEX``) and the documented ``bond_channels``
functions. ``bubble.py``'s own names carry no external compatibility
promise and may be renamed or restructured within the series without
notice to callers outside ``hwave.solver``.

The binding contract is
``docs/superpowers/specs/2026-08-14-unified-bubble-kernel-design.md``
(sections "Module layout", "Green/tail contracts", "Scheme x bond", and
"Output contracts"). This module transcribes the numerics of
``RPA._calc_chi0q`` (rpa.py, around line 2936 at series start) without
altering them -- including the equal-time endpoint-mean correction
(issue #134) and the data-driven tail gate. ``bubble.py`` deliberately
imports neither ``sparse_ir`` nor ``flex.py``, so the plain dense path
stays usable without the optional IR dependency.
"""

import numpy as np

from . import backend as _bk
from . import matsubara as _ms
from .kgrid import reverse_fft_axes

_SCHEMES = ("reduced", "general")


def contract_reduced(g_fwd, g_rev):
    """Reduced (diagonal-orbital-pair) particle-hole contraction.

    ``chi0[..., a, b] = g_fwd[..., a, b] * g_rev[..., b, a]``. Generic over
    any shared leading shape -- used both on the full
    ``(nblock, nmat, nvol, p, p)`` tensors and on the single-frequency
    equal-time endpoint slices (``(nblock, nvol, p, p)``).
    """
    return g_fwd * g_rev.swapaxes(-2, -1)


def contract_general(g_fwd, g_rev):
    """General (full orbital-pair) particle-hole contraction.

    ``chi0[..., a, c, b, d] = g_fwd[..., a, b] * g_rev[..., d, c]``, built
    by an outer product (broadcasting) followed by an axis swap -- NOT
    einsum, for parity with the dense general branch's existing
    performance choice (``RPA._calc_chi0q``'s general path,
    rpa.py:3091-3114). Generic over any shared leading shape, same as
    :func:`contract_reduced`.
    """
    # g_fwd[..., a, None, b, None] * g_rev[..., None, d, None, c]
    #   -> (..., a, d, b, c)
    x = (g_fwd[..., :, np.newaxis, :, np.newaxis]
         * g_rev[..., np.newaxis, :, np.newaxis, :])
    # (..., a, d, b, c) -> (..., a, c, b, d): swap the trailing d/c axes.
    return x.swapaxes(-3, -1)


def _validate_scheme(scheme):
    if scheme not in _SCHEMES:
        raise ValueError(
            "bubble kernel: scheme must be one of {!r}, got {!r}".format(
                _SCHEMES, scheme))


def _validate_beta(beta):
    if isinstance(beta, bool) or not isinstance(
            beta, (int, float, np.integer, np.floating)):
        raise ValueError(
            "bubble kernel: beta must be a real number, got {!r}".format(beta))
    beta = float(beta)
    if not (beta > 0.0):
        raise ValueError(
            "bubble kernel: beta must be > 0, got {}".format(beta))
    return beta


def _validate_spatial_shape(spatial_shape, nvol):
    try:
        shape = tuple(spatial_shape)
    except TypeError:
        raise ValueError(
            "bubble kernel: spatial_shape must be a 3-tuple of positive "
            "integers, got {!r}".format(spatial_shape))
    if len(shape) != 3 or any(
            isinstance(s, bool) or not isinstance(s, (int, np.integer)) or s <= 0
            for s in shape):
        raise ValueError(
            "bubble kernel: spatial_shape must be a 3-tuple of positive "
            "integers, got {!r}".format(spatial_shape))
    shape = tuple(int(s) for s in shape)
    prod = shape[0] * shape[1] * shape[2]
    if prod != nvol:
        raise ValueError(
            "bubble kernel: spatial_shape {} (volume {}) does not match "
            "the Green function's volume axis ({})".format(
                shape, prod, nvol))
    return shape


def _promote_complex(arr, xp, label):
    if arr.dtype.kind != 'c':
        raise ValueError(
            "bubble kernel: {} must be complex, got dtype {}".format(
                label, arr.dtype))
    # complex64 is promoted to complex128 (a copy); complex128 input is
    # returned unchanged (astype(..., copy=False) is then a no-op).
    return arr.astype(xp.complex128, copy=False)


def _validate_dense_inputs(green_kw, green0_tail, beta, spatial_shape, scheme):
    """Validate a ``dense_bubble`` call. Returns the normalized
    ``(green_kw, green0_tail, beta, nblock, nmat, nvol, nd, spatial_shape)``
    -- ``green_kw``/``green0_tail`` promoted to complex128 (copy only when
    needed). All failures are ``ValueError`` (survives ``python -O``)."""
    _validate_scheme(scheme)
    beta = _validate_beta(beta)

    xp = _bk.array_module_of(green_kw)

    if green_kw.ndim != 5:
        raise ValueError(
            "bubble kernel: green_kw must have ndim 5 "
            "(nblock, nmat, nvol, p, p), got ndim={} shape={}".format(
                green_kw.ndim, green_kw.shape))
    nblock, nmat, nvol, nd, nd2 = green_kw.shape
    if nblock not in (1, 2):
        # 1 = spin-free/spinful, 2 = spin-diag. See RPA._calc_chi0q's
        # identical guard (rpa.py issue #125) for why this is a ValueError,
        # not a bare assert, and why 0/3+ are rejected rather than
        # silently truncated by a caller's block handling.
        raise ValueError(
            "bubble kernel: Green's function block axis ({}) must be 1 "
            "(spin-free/spinful) or 2 (spin-diag)".format(nblock))
    if nmat <= 0 or nmat % 2 != 0:
        raise ValueError(
            "bubble kernel: Green's function frequency axis (nmat={}) "
            "must be a positive even integer".format(nmat))
    if nd != nd2 or nd < 1:
        raise ValueError(
            "bubble kernel: orbital axes must be square and nonempty, "
            "got ({}, {})".format(nd, nd2))

    spatial_shape = _validate_spatial_shape(spatial_shape, nvol)

    green_kw = _promote_complex(green_kw, xp, "green_kw")

    if green0_tail is not None:
        # Structural pairing check only (mirrors RPA._calc_chi0q's
        # green0_tail guard): shape/dtype/backend equality cannot prove
        # the tail came from the same Green build -- that provenance is
        # the CALLER's responsibility (see the module docstring / spec's
        # Green/tail contract).
        if green0_tail.shape != green_kw.shape:
            raise ValueError(
                "bubble kernel: green0_tail shape {} does not match the "
                "Green function {} -- the tail must be the paired array "
                "from the same Green build".format(
                    green0_tail.shape, green_kw.shape))
        if _bk.array_module_of(green0_tail) is not xp:
            raise ValueError(
                "bubble kernel: green0_tail and green_kw must live on the "
                "same array backend")
        green0_tail = _promote_complex(green0_tail, xp, "green0_tail")

    return (green_kw, green0_tail, beta, nblock, nmat, nvol, nd, spatial_shape)


class _DensePrepared:
    """Internal container for :func:`_prepare_dense`'s output. Not part of
    this module's public surface."""

    __slots__ = ("green_rt", "green_rev", "sgn", "tail_on",
                 "fwd0_p", "rev0_p", "jump_f", "jump_r_rev")

    def __init__(self, green_rt, green_rev, sgn, tail_on,
                 fwd0_p, rev0_p, jump_f, jump_r_rev):
        self.green_rt = green_rt
        self.green_rev = green_rev
        self.sgn = sgn
        self.tail_on = tail_on
        self.fwd0_p = fwd0_p
        self.rev0_p = rev0_p
        self.jump_f = jump_f
        self.jump_r_rev = jump_r_rev


def _prepare_dense(green_kw, green0_tail, beta, spatial_shape, workers):
    """Fermion-to-tau, conditional tau-space tail subtraction, spatial
    ifftn, the FFT-grid reversal (``kgrid.reverse_fft_axes``), the
    fermionic tau-domain sign, and the equal-time endpoint-branch metadata
    (issue #134). Transcribes ``RPA._calc_chi0q`` (rpa.py:3020-3071)
    verbatim.

    ``green_kw``/``green0_tail`` are assumed already validated and
    complex128-promoted (by :func:`_validate_dense_inputs`); this function
    performs no validation of its own. ``beta`` is accepted for interface
    parity with the IR transport's ``_prepare_ir`` (a later task in the
    series) -- the dense prefix computed here does not itself need it.
    """
    xp = _bk.array_module_of(green_kw)
    nblock, nmat, nvol, nd, _ = green_kw.shape
    nx, ny, nz = spatial_shape

    green_flat = green_kw.reshape(nblock, nmat, nvol * nd * nd)
    green_kt = _ms.fermion_to_tau(green_flat, axis=1)
    green_kt = green_kt.reshape(nblock, nmat, nx, ny, nz, nd, nd)
    if green0_tail is not None:
        green_kt = green_kt - green0_tail.reshape(
            nblock, nmat, nx, ny, nz, nd, nd)

    green_rt = _bk.spatial_ifftn(
        green_kt.reshape(nblock, nmat, nx, ny, nz, nd * nd),
        axes=(2, 3, 4), workers=workers)

    green_rev = reverse_fft_axes(green_rt, (1, 2, 3, 4)).reshape(
        nblock, nmat, nvol, nd, nd)

    # Equal-time endpoint correction (issue #134): see rpa.py:3034-3056 for
    # the full derivation. No-op -- including bitwise -- when the tail is
    # zero. The predicate inspects only frequency slice 0: that is the
    # ONLY slice the correction reads (the tail is frequency-constant by
    # construction), so reducing the full repeated tensor would force a
    # needless device synchronization on the GPU path.
    tail_on = False
    fwd0_p = rev0_p = jump_f = jump_r_rev = None
    if green0_tail is not None:
        tail_on = bool(xp.any(
            green0_tail.reshape(nblock, nmat, -1)[:, 0] != 0))
        if tail_on:
            tail_kt0 = green0_tail.reshape(
                nblock, nmat, nx, ny, nz, nd * nd)[:, 0:1]
            jump_f = _bk.spatial_ifftn(2.0 * tail_kt0, axes=(2, 3, 4),
                                        workers=workers)
            # the reversed factor lives at -r: its jump is the spatial
            # reversal of the forward one
            jump_r_rev = reverse_fft_axes(jump_f, (2, 3, 4))
            jump_f = jump_f.reshape(nblock, nvol, nd, nd)
            jump_r_rev = jump_r_rev.reshape(nblock, nvol, nd, nd)
            fwd0_p = green_rt.reshape(nblock, nmat, nvol, nd, nd)[:, 0].copy()
            rev0_p = green_rev[:, 0].copy()

    sgn = xp.full(nmat, -1)
    sgn[0] = 1

    green_rt5 = green_rt.reshape(nblock, nmat, nvol, nd, nd)

    return _DensePrepared(green_rt5, green_rev, sgn, tail_on,
                           fwd0_p, rev0_p, jump_f, jump_r_rev)


def _assemble_plain(prepped, scheme, beta, spatial_shape, workers):
    """Contract (reduced or general), apply the equal-time endpoint-mean
    correction, and transport back to (bosonic Matsubara frequency,
    k-space) -- the dense arm. Transcribes ``RPA._calc_chi0q``
    (rpa.py:3076-3126) verbatim."""
    nx, ny, nz = spatial_shape
    nblock, nmat, nvol, nd, _ = prepped.green_rt.shape

    contract = contract_reduced if scheme == "reduced" else contract_general
    nd_shape = (nd, nd) if scheme == "reduced" else (nd, nd, nd, nd)
    nds = nd * nd if scheme == "reduced" else nd ** 4

    chi0_rt = contract(prepped.green_rt, prepped.green_rev)
    sgn_bc = prepped.sgn.reshape((1, nmat) + (1,) * (chi0_rt.ndim - 2))
    chi0_rt = chi0_rt * sgn_bc

    if prepped.tail_on:
        # mean of the two equal-time branches (sgn[0] = +1)
        chi0_rt[:, 0] = 0.5 * (
            contract(prepped.fwd0_p, prepped.rev0_p + prepped.jump_r_rev)
            + contract(prepped.fwd0_p + prepped.jump_f, prepped.rev0_p))

    chi0_qt = _bk.spatial_fftn(
        chi0_rt.reshape(nblock, nmat, nx, ny, nz, nds),
        axes=(2, 3, 4), workers=workers)

    chi0_qt_flat = chi0_qt.reshape(nblock, nmat, nvol * nds)
    chi0_qw = _ms.tau_to_boson(chi0_qt_flat, axis=1).reshape(
        nblock, nmat, nvol, *nd_shape) * (-1.0 / beta)

    return chi0_qw


def dense_bubble(green_kw, green0_tail, beta, *, spatial_shape, scheme,
                  workers=None):
    """Dense uniform-Matsubara-grid particle-hole bubble chi0(q, iW_m).

    Parameters
    ----------
    green_kw : ndarray, complex, shape (nblock, nmat, nvol, p, p)
        Green's function in k-space and fermionic Matsubara frequency.
        When ``green0_tail`` is not None this is the DEFLATED Green
        function (``coeff_tail / (iw_n)`` already subtracted); when
        ``green0_tail`` is None it is the FULL Green function and no tail
        machinery runs. ``nblock`` is 1 (spin-free/spinful) or 2
        (spin-diag).
    green0_tail : ndarray or None, same shape/dtype family as ``green_kw``
        The paired tau-space tail add-back from the same Green build, or
        None to disable the tail correction. Shape/dtype/backend equality
        with ``green_kw`` is validated, but this cannot prove the pair
        came from the same Green build -- that provenance is the caller's
        responsibility (see the module docstring).
    beta : float
        Inverse temperature; must be > 0.
    spatial_shape : tuple of 3 positive ints (Nx, Ny, Nz)
        The lattice shape; ``prod(spatial_shape)`` must equal
        ``green_kw.shape[2]`` (``nvol``).
    scheme : {"reduced", "general"}
        ``"reduced"`` returns the diagonal-orbital-pair bubble
        ``(nblock, nmat, nvol, p, p)``; ``"general"`` returns the full
        orbital-pair bubble ``(nblock, nmat, nvol, p, p, p, p)``.
    workers : int or None, optional
        Forwarded to ``backend.spatial_ifftn``/``spatial_fftn`` exactly as
        ``RPA._calc_chi0q`` forwards its own ``fft_workers`` setting.

    Returns
    -------
    ndarray, complex128
        The bare susceptibility chi0(q, iW_m), bosonic Matsubara frequency
        on axis 1 (uniform grid, same ``nmat`` as the input), k-space on
        axis 2. Inputs are promoted to complex128 (complex64 is copied
        up); real-dtype inputs raise ``ValueError``.

    Raises
    ------
    ValueError
        On any shape/dtype/scheme mismatch (see the validation performed
        by :func:`_validate_dense_inputs`); these checks use ``ValueError``
        rather than a bare ``assert`` so they survive ``python -O``.
    """
    (green_kw, green0_tail, beta, nblock, nmat, nvol, nd, spatial_shape
     ) = _validate_dense_inputs(green_kw, green0_tail, beta, spatial_shape,
                                 scheme)

    prepped = _prepare_dense(green_kw, green0_tail, beta, spatial_shape,
                              workers)
    return _assemble_plain(prepped, scheme, beta, spatial_shape, workers)
