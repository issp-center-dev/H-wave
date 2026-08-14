"""INTERNAL module -- shared particle-hole bubble orbital contractions and
the dense uniform-Matsubara-grid / sparse-IR transports.

This module is NOT public API; the supported user-facing surface stays the
solver classes (``RPA``, ``FLEX``) and the documented ``bond_channels``
functions. ``bubble.py``'s own names carry no external compatibility
promise and may be renamed or restructured within the series without
notice to callers outside ``hwave.solver``.

The binding contract is
``docs/superpowers/specs/2026-08-14-unified-bubble-kernel-design.md``
(sections "Module layout", "Green/tail contracts", "Scheme x bond", "The
bond assembly", and "Output contracts"). This module transcribes the
numerics of ``RPA._calc_chi0q`` (rpa.py, around line 2936 at series start),
``FLEX._calc_chi0q_ir`` / ``FLEX._calc_chi0q_general_ir`` (flex.py, around
lines 818-917 at series start), and ``bond_channels.bond_bubble``
(bond_channels.py, around lines 1171-1315 at series start) without
altering them -- including the equal-time endpoint-mean correction (issue
#134) and the data-driven tail gate on the dense path. ``bubble.py``
deliberately imports neither ``sparse_ir`` nor ``flex.py`` -- the IR entry
points receive already-constructed ``IRAxis`` objects and duck-type their
attributes/methods (``statistics``, ``beta``, ``wmax``, ``eps``, ``tau``,
``n_tau``, ``n_freq``, ``freq_to_tau_points``, ``tau_to_freq``) -- so the
plain dense path stays usable without the optional IR dependency; the bond
entry point takes a ``bond_set`` (duck-typed ``delta_r``/``n_channels``)
and never imports ``bond_channels`` either.
"""

import numpy as np

from . import backend as _bk
from . import matsubara as _ms
from .kgrid import reverse_fft_axes

_SCHEMES = ("reduced", "general")

# Per-pair working-set buffer counts for the streaming bond assembly (spec:
# "Memory contract"). This module owns the canonical definition;
# ``bond_channels.py`` re-exports these names (unified-bubble-kernel series,
# Task 8) rather than keeping its own copy, so the two cannot silently
# desync.
#
# N4_BUFFERS is 4, not the legacy ``bond_bubble``'s 3: ``contract_general``
# returns a non-contiguous SWAPAXES VIEW (the ``(..., a, d, b, c) ->
# (..., a, c, b, d)`` permutation), and ``_bond_pair_full_block``'s
# immediately following ``.reshape(nmat, nx, ny, nz, npair, npair)`` cannot
# express that as a further view -- it silently COPIES, so the pre-reshape
# buffer and its reshaped copy briefly coexist as a 4th N4-sized array
# (measured: ``tests/test_sc_bond.py::TestPreflightMemoryBudget``). The
# legacy body avoided this because ``np.einsum`` wrote directly into the
# final contiguous layout; the shared kernel trades that one extra
# transient for sharing ``contract_general`` with the plain dense path.
BOND_BUBBLE_N4_BUFFERS = 4
BOND_BUBBLE_N2_BUFFERS = 3


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
    if nmat <= 0:
        # Evenness is deliberately NOT required here (spec: "Kernel input
        # validation", AMENDED 2026-08-14, Task 6): the centered grid
        # iw_n = 1j*(2n+1-nmat)*pi/beta is well-defined for odd nmat too,
        # the legacy _calc_chi0q body accepted odd nmat, and the issue-#91
        # orbital-order regression locks in test_flex_general.py run a
        # degenerate nmat=1 oracle through the public dispatch path. The
        # bond entry points (bond_bubble_static, _iter_bond_dynamic) DO
        # still require even nmat -- see _validate_even_nmat_for_bond --
        # because their Omega=0 identification (static_index = nmat // 2)
        # assumes the even centered grid.
        raise ValueError(
            "bubble kernel: Green's function frequency axis (nmat={}) "
            "must be a positive integer".format(nmat))
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


def _validate_even_nmat_for_bond(nmat):
    """The bond entry points (``bond_bubble_static``, ``_iter_bond_dynamic``)
    require an even ``nmat``, unlike the plain dense path (see the
    evenness note in :func:`_validate_dense_inputs`): their static
    ``Omega=0`` slice is identified by ``static_index = nmat // 2``, which
    only lands on the bosonic zero-frequency point for the even centered
    Matsubara grid. ``ValueError`` (survives ``python -O``)."""
    if nmat % 2 != 0:
        raise ValueError(
            "bubble kernel: the bond bubble entry points require an even "
            "nmat (got nmat={}) -- the static Omega=0 slice is identified "
            "by static_index = nmat // 2, which assumes the even centered "
            "Matsubara grid".format(nmat))


def _validate_bond_set(bond_set):
    """Structural guard on ``bond_set`` (a
    ``bond_channels.ResolvedInteractionSet``, duck-typed here -- this
    module never imports ``bond_channels``): ``n_channels``/``delta_r``
    present, ``len(delta_r) == n_channels``, and every displacement an
    integer 3-tuple (spec: "Kernel input validation"). Returns
    ``(n_channels, delta_r)``. ``ValueError`` on any failure (survives
    ``python -O``)."""
    n_channels = getattr(bond_set, "n_channels", None)
    delta_r = getattr(bond_set, "delta_r", None)
    if n_channels is None or delta_r is None:
        raise ValueError(
            "bubble kernel: bond_set must provide 'n_channels' and "
            "'delta_r' attributes (see "
            "bond_channels.resolve_interactions)")
    if len(delta_r) != n_channels:
        raise ValueError(
            "bubble kernel: bond_set.delta_r length ({}) does not match "
            "bond_set.n_channels ({})".format(len(delta_r), n_channels))
    for i, dr in enumerate(delta_r):
        try:
            ok = len(dr) == 3 and all(
                not isinstance(c, bool) and isinstance(c, (int, np.integer))
                for c in dr)
        except TypeError:
            ok = False
        if not ok:
            raise ValueError(
                "bubble kernel: bond_set.delta_r[{}] must be a 3-tuple of "
                "integers, got {!r}".format(i, dr))
    return n_channels, delta_r


def _validate_ir_axes(ax_fermi, ax_bose):
    """Axis-compatibility guard for the IR entry points (spec: "IR axis
    compatibility"). Duck-typed -- reads ``statistics``/``beta``/``wmax``/
    ``eps`` off the given objects without importing ``ir_axis`` or
    ``sparse_ir``. Validates PARAMETER equality only: separately
    constructed axes with equal parameters pass. ``ValueError`` names the
    offending field and both values."""
    fermi_stat = getattr(ax_fermi, "statistics", None)
    if fermi_stat != "F":
        raise ValueError(
            "bubble kernel: ax_fermi.statistics must be 'F', got {!r}"
            .format(fermi_stat))
    bose_stat = getattr(ax_bose, "statistics", None)
    if bose_stat != "B":
        raise ValueError(
            "bubble kernel: ax_bose.statistics must be 'B', got {!r}"
            .format(bose_stat))
    for field in ("beta", "wmax", "eps"):
        v_fermi = getattr(ax_fermi, field)
        v_bose = getattr(ax_bose, field)
        if v_fermi != v_bose:
            raise ValueError(
                "bubble kernel: ax_fermi.{f} ({a!r}) does not match "
                "ax_bose.{f} ({b!r}) -- the IR entry points require "
                "parameter-identical fermionic/bosonic axes"
                .format(f=field, a=v_fermi, b=v_bose))


def _validate_ir_inputs(green_kw, ax_fermi, ax_bose, spatial_shape, scheme):
    """Validate an ``ir_bubble`` call. Returns the normalized
    ``(green_kw, nblock, nw, nvol, nd, spatial_shape)`` -- ``green_kw``
    promoted to complex128 (copy only when needed). All failures are
    ``ValueError`` (survives ``python -O``). No tail argument -- the IR
    entry points take the FULL Green function at fermionic IR nodes (see
    the module docstring's Green/tail contract)."""
    _validate_scheme(scheme)
    _validate_ir_axes(ax_fermi, ax_bose)

    xp = _bk.array_module_of(green_kw)

    if green_kw.ndim != 5:
        raise ValueError(
            "bubble kernel: green_kw must have ndim 5 "
            "(nblock, n_freq_F, nvol, p, p), got ndim={} shape={}".format(
                green_kw.ndim, green_kw.shape))
    nblock, nw, nvol, nd, nd2 = green_kw.shape
    if nblock not in (1, 2):
        raise ValueError(
            "bubble kernel: Green's function block axis ({}) must be 1 "
            "(spin-free/spinful) or 2 (spin-diag)".format(nblock))
    if nd != nd2 or nd < 1:
        raise ValueError(
            "bubble kernel: orbital axes must be square and nonempty, "
            "got ({}, {})".format(nd, nd2))
    if nw != ax_fermi.n_freq:
        raise ValueError(
            "bubble kernel: green_kw frequency axis (nw={}) does not "
            "match ax_fermi.n_freq ({}) -- green_kw must be sampled at "
            "the fermionic IR nodes".format(nw, ax_fermi.n_freq))

    spatial_shape = _validate_spatial_shape(spatial_shape, nvol)

    green_kw = _promote_complex(green_kw, xp, "green_kw")

    return (green_kw, nblock, nw, nvol, nd, spatial_shape)


class _IRPrepared:
    """Internal container for :func:`_prepare_ir`'s output. Not part of
    this module's public surface."""

    __slots__ = ("green_rt", "green_rev")

    def __init__(self, green_rt, green_rev):
        self.green_rt = green_rt
        self.green_rev = green_rev


def _prepare_ir(green_kw, ax_fermi, ax_bose, spatial_shape, workers):
    """Fermionic-IR-node-to-bosonic-tau evaluation, spatial ifftn, and the
    IR tau/spatial reversal. Transcribes ``FLEX._calc_chi0q_ir`` /
    ``FLEX._calc_chi0q_general_ir`` (flex.py:842-856 / 889-902) verbatim.

    Unlike the dense transport's ``_prepare_dense`` (a roll+flip over the
    discrete FFT tau grid via ``kgrid.reverse_fft_axes`` on BOTH the tau
    and spatial axes at once), the IR tau reversal is a plain reflection
    ``flip`` of the (already reversal-symmetric) bosonic tau node array,
    kept as its own step distinct from the shared FFT-grid spatial
    reversal -- the two transports never share a reversal implementation
    (spec: "Module layout").

    ``green_kw`` is assumed already validated and complex128-promoted (by
    :func:`_validate_ir_inputs`); this function performs no validation of
    its own. ``ax_fermi``/``ax_bose`` are duck-typed ``IRAxis``-like
    objects (see the module docstring)."""
    xp = _bk.array_module_of(green_kw)
    nblock, nw, nvol, nd, _ = green_kw.shape
    nx, ny, nz = spatial_shape

    # G(k, tau_B): fermionic-node coefficients evaluated at the bosonic
    # axis's tau nodes.
    g = xp.moveaxis(green_kw.reshape(nblock, nw, nvol * nd * nd), 1, -1)
    g_tau = ax_fermi.freq_to_tau_points(g, ax_bose.tau)
    ntB = ax_bose.n_tau
    g_tau = xp.moveaxis(g_tau, -1, 1).reshape(
        nblock, ntB, nx, ny, nz, nd * nd)

    g_rt = _bk.spatial_ifftn(g_tau, axes=(2, 3, 4), workers=workers)

    # G(-r, -tau) = -G(-r, beta - tau): interior symmetric tau nodes ->
    # plain tau reversal (with the global fermionic -1); r -> -r is the
    # shared FFT-grid reversal (identical to the dense path).
    g_rev = -xp.flip(reverse_fft_axes(g_rt, (2, 3, 4)), axis=1)
    g_rt = g_rt.reshape(nblock, ntB, nvol, nd, nd)
    g_rev = g_rev.reshape(nblock, ntB, nvol, nd, nd)

    return _IRPrepared(g_rt, g_rev)


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
    # Never needed again (fully consumed by the ifftn above); released here
    # rather than left bound for the rest of the function -- the bond
    # path's memory budget (BOND_BUBBLE_N2_BUFFERS, sc._bond_memory_estimate)
    # assumes this buffer does not coexist with green_rev below.
    del green_kt

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


def _assemble_plain(prepped, scheme, beta, spatial_shape, workers,
                    ax_bose=None):
    """Contract (reduced or general) and transport back to (bosonic
    Matsubara frequency, k-space).

    Two arms, selected by ``ax_bose`` (``None`` -> dense, given -> IR):

    - Dense arm (``ax_bose is None``): applies the equal-time endpoint-mean
      correction and the per-tau ``sgn``, then ``spatial fftn ->
      tau_to_boson -> x(-1/beta)``. Transcribes ``RPA._calc_chi0q``
      (rpa.py:3076-3126) verbatim.
    - IR arm (``ax_bose`` given): the bubble sign is the single explicit
      ``-1`` on the contraction (no ``sgn`` array, no endpoint correction --
      the IR basis already represents the tail), then ``spatial fftn ->
      ax_bose.tau_to_freq`` with NO ``1/beta`` factor (the IR transforms
      are physical). Transcribes ``FLEX._calc_chi0q_ir`` /
      ``FLEX._calc_chi0q_general_ir`` (flex.py:858-868 / 904-917) verbatim.
    """
    nx, ny, nz = spatial_shape
    nblock, nt, nvol, nd, _ = prepped.green_rt.shape

    contract = contract_reduced if scheme == "reduced" else contract_general
    nd_shape = (nd, nd) if scheme == "reduced" else (nd, nd, nd, nd)
    nds = nd * nd if scheme == "reduced" else nd ** 4

    if ax_bose is not None:
        xp = _bk.array_module_of(prepped.green_rt)
        chi0_rt = -contract(prepped.green_rt, prepped.green_rev)

        chi0_qt = _bk.spatial_fftn(
            chi0_rt.reshape(nblock, nt, nx, ny, nz, nds),
            axes=(2, 3, 4), workers=workers)

        chi0_qt_flat = xp.moveaxis(
            chi0_qt.reshape(nblock, nt, nvol * nds), 1, -1)
        chi0_q = ax_bose.tau_to_freq(chi0_qt_flat)
        chi0_q = xp.moveaxis(chi0_q, -1, 1).reshape(
            nblock, ax_bose.n_freq, nvol, *nd_shape)
        return chi0_q

    chi0_rt = contract(prepped.green_rt, prepped.green_rev)
    sgn_bc = prepped.sgn.reshape((1, nt) + (1,) * (chi0_rt.ndim - 2))
    chi0_rt = chi0_rt * sgn_bc

    if prepped.tail_on:
        # mean of the two equal-time branches (sgn[0] = +1)
        chi0_rt[:, 0] = 0.5 * (
            contract(prepped.fwd0_p, prepped.rev0_p + prepped.jump_r_rev)
            + contract(prepped.fwd0_p + prepped.jump_f, prepped.rev0_p))

    chi0_qt = _bk.spatial_fftn(
        chi0_rt.reshape(nblock, nt, nx, ny, nz, nds),
        axes=(2, 3, 4), workers=workers)

    chi0_qt_flat = chi0_qt.reshape(nblock, nt, nvol * nds)
    chi0_qw = _ms.tau_to_boson(chi0_qt_flat, axis=1).reshape(
        nblock, nt, nvol, *nd_shape) * (-1.0 / beta)

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


def ir_bubble(green_kw, ax_fermi, ax_bose, *, spatial_shape, scheme,
              workers=None):
    """Sparse-IR particle-hole bubble chi0(q, iOmega_l) on the bosonic IR
    nodes, computed natively from the full Green function on the fermionic
    IR nodes.

    Parameters
    ----------
    green_kw : ndarray, complex, shape (nblock, n_freq_F, nvol, p, p)
        The FULL Green's function (no tail subtracted) sampled at
        ``ax_fermi``'s fermionic Matsubara nodes -- there is no
        ``green0_tail`` argument on the IR path (the IR basis represents
        the ``coeff_tail / (i w)`` asymptotics within ``ax_fermi.eps``).
        ``nblock`` is 1 (spin-free/spinful) or 2 (spin-diag).
    ax_fermi, ax_bose : IRAxis-like
        Already-constructed fermionic/bosonic IR axis objects (duck-typed:
        ``statistics``, ``beta``, ``wmax``, ``eps``, ``tau``, ``n_tau``,
        ``n_freq``, ``freq_to_tau_points``, ``tau_to_freq`` -- this module
        never imports ``sparse_ir`` or ``ir_axis`` itself). Must be
        parameter-compatible: ``ax_fermi.statistics == "F"``,
        ``ax_bose.statistics == "B"``, and exact float equality of
        ``beta``, ``wmax``, and ``eps`` between the two -- see
        :func:`_validate_ir_axes`.
    spatial_shape : tuple of 3 positive ints (Nx, Ny, Nz)
        The lattice shape; ``prod(spatial_shape)`` must equal
        ``green_kw.shape[2]`` (``nvol``).
    scheme : {"reduced", "general"}
        ``"reduced"`` returns the diagonal-orbital-pair bubble
        ``(nblock, n_freq_B, nvol, p, p)``; ``"general"`` returns the full
        orbital-pair bubble ``(nblock, n_freq_B, nvol, p, p, p, p)``.
    workers : int or None, optional
        Forwarded to ``backend.spatial_ifftn``/``spatial_fftn`` exactly as
        the FLEX IR methods forward their own ``fft_workers`` setting.

    Returns
    -------
    ndarray, complex128
        The bare susceptibility chi0(q, iOmega_l) on ``ax_bose``'s bosonic
        Matsubara nodes (axis 1, length ``ax_bose.n_freq``), k-space on
        axis 2. Inputs are promoted to complex128 (complex64 is copied up);
        real-dtype inputs raise ``ValueError``. No ``1/beta`` factor is
        applied -- the IR transforms are physical (the module docstring's
        "Green/tail contracts" section).

    Raises
    ------
    ValueError
        On any shape/dtype/scheme/axis-compatibility mismatch (see
        :func:`_validate_ir_inputs` / :func:`_validate_ir_axes`); these
        checks use ``ValueError`` rather than a bare ``assert`` so they
        survive ``python -O``.
    """
    (green_kw, nblock, nw, nvol, nd, spatial_shape
     ) = _validate_ir_inputs(green_kw, ax_fermi, ax_bose, spatial_shape,
                             scheme)

    prepped = _prepare_ir(green_kw, ax_fermi, ax_bose, spatial_shape,
                          workers)
    return _assemble_plain(prepped, scheme, None, spatial_shape, workers,
                           ax_bose=ax_bose)


# ===========================================================================
# bond-enlarged bubble  --  streaming (m, m') channel-pair assembly
# ===========================================================================
#
# Transcribes ``bond_channels.bond_bubble`` (bond_channels.py:1171-1315)
# onto the shared ``_prepare_dense``/``contract_general`` primitives: the
# roll-by-``Delta r_m - Delta r_m'`` of the reversed Green function, the
# general orbital contraction, and the reshape into the ``(npair, npair)``
# pair block, one channel pair at a time so the per-pair working set stays
# O(nmat * nvol * p**4) rather than O(nmat * nvol * ND**2) (spec: "The bond
# assembly"). The STATIC entry point (``bond_bubble_static``, Omega=0 only)
# and the DYNAMIC entry point (``_iter_bond_dynamic``, every bosonic
# Matsubara frequency) share this ONE per-pair pipeline
# (``_bond_pair_full_block`` via ``_iter_bond_channel_pairs``) -- static
# reads only the Omega=0 slice of each pair's full block, keeping its
# zero-materialization property at the ND**2 scale.


def _roll_spatial(arr, axis, spatial_shape, shift):
    """Circularly shift ``arr`` by ``shift`` (a 3-tuple of integers) along
    the flattened spatial axis at position ``axis`` -- temporarily splits
    that axis into ``(Nx, Ny, Nz)``, rolls, and reshapes back. Used for the
    per-channel-pair ``Delta r_m - Delta r_m'`` shift (the Fourier-shift
    realization of the bond phase ``e^{i k . (Delta r_m - Delta r_m')}``,
    see ``bond_channels.bond_bubble``'s docstring)."""
    xp = _bk.array_module_of(arr)
    nx, ny, nz = spatial_shape
    shape = arr.shape
    expanded = arr.reshape(shape[:axis] + (nx, ny, nz) + shape[axis + 1:])
    rolled = xp.roll(expanded, shift=shift, axis=(axis, axis + 1, axis + 2))
    return rolled.reshape(shape)


def _bond_pair_full_block(green_fwd_sgn, green_rev, tail_pack, beta,
                          spatial_shape, workers, shift):
    """Compute the FULL-frequency bond bubble block for ONE ``(m, m')``
    channel pair: ``(nmat, nvol, npair, npair)`` complex128, FRESHLY
    ALLOCATED (never a view into a reused buffer -- every array computed
    here, including the returned one, is a fresh allocation from an
    elementwise/FFT op, so no caller-side copy is needed).

    This is the ONE private per-pair pipeline shared by
    :func:`_assemble_bond_static` (which reads only the ``Omega=0`` slice,
    keeping its zero-materialization property at the ``ND**2`` scale) and
    :func:`_iter_bond_dynamic` (which yields the full block) -- spec:
    "static and dynamic share ONE private per-pair helper".

    Roll the reversed Green (and, when the tail correction is active, the
    reversed branch's equal-time jump term -- SAME shift, spec: "endpoint
    correction under a shift") by ``shift`` (``Delta r_m - Delta r_m'``),
    ``contract_general``, overwrite the equal-time (tau=0) slice with the
    endpoint-mean correction, reshape the trailing ``(a, c, b, d)`` axes
    into the ``(npair, npair)`` pair block (spec's "Scheme x bond and the
    pair-block mapping" -- ``contract_general``'s own axis order already
    IS the row/column split, no moveaxis needed), spatial ``fftn``,
    ``tau_to_boson``, scale by ``-1/beta``, and free every per-pair
    temporary before returning (mirrors ``bond_channels.bond_bubble``'s
    ``del`` discipline -- the ``BOND_BUBBLE_N4_BUFFERS`` /
    ``BOND_BUBBLE_N2_BUFFERS`` memory contract depends on it).

    ``green_fwd_sgn``/``green_rev``: ``(nmat, nvol, nd, nd)`` (block axis
    already dropped by the caller). ``tail_pack``: ``None`` (tail off) or
    the 4-tuple ``(fwd0_p, rev0_p, jump_f, jump_r_rev)``, each
    ``(nvol, nd, nd)`` (block axis already dropped)."""
    nmat, nvol, nd, _ = green_rev.shape
    nx, ny, nz = spatial_shape
    npair = nd * nd

    green_rev_shifted = _roll_spatial(green_rev, 1, spatial_shape, shift)

    # (nmat, nvol, a, c, b, d)
    chi0_rt = contract_general(green_fwd_sgn, green_rev_shifted)
    del green_rev_shifted

    if tail_pack is not None:
        fwd0_p, rev0_p, jump_f, jump_r_rev = tail_pack
        rev0_p_shifted = _roll_spatial(rev0_p, 0, spatial_shape, shift)
        jump_r_rev_shifted = _roll_spatial(
            jump_r_rev, 0, spatial_shape, shift)
        chi0_rt[0] = 0.5 * (
            contract_general(fwd0_p, rev0_p_shifted + jump_r_rev_shifted)
            + contract_general(fwd0_p + jump_f, rev0_p_shifted))
        del rev0_p_shifted, jump_r_rev_shifted

    # reshape to the pair block: contract_general's trailing axis order is
    # already (a, c, b, d) -- row = (a, c), col = (b, d) (spec's
    # pair-block mapping) -- and split nvol back to (Nx, Ny, Nz) for the
    # spatial FFT, in one reshape.
    chi0_rt = chi0_rt.reshape(nmat, nx, ny, nz, npair, npair)

    chi0_qt = _bk.spatial_fftn(chi0_rt, axes=(1, 2, 3), workers=workers)
    del chi0_rt
    chi0_qt = chi0_qt.reshape(nmat, nvol, npair, npair)

    chi0_qw = _ms.tau_to_boson(chi0_qt, axis=0)
    del chi0_qt

    block = chi0_qw * (-1.0 / beta)
    del chi0_qw

    return block


def _iter_bond_channel_pairs(prepped, bond_set_validated, beta,
                             spatial_shape, workers):
    """Shared per-pair generator: yields ``((m, m'), block)`` in row-major
    channel order (``m`` outer, ``m'`` inner), ``block`` the FULL-frequency
    ``(nmat, nvol, npair, npair)`` complex128 bond bubble for that channel
    pair (see :func:`_bond_pair_full_block`), on ``_prepare_dense``'s
    already-prepared tau/real-space tensors (``nblock`` assumed already
    validated as 1 by the caller). ``bond_set_validated`` is the
    ``(n_channels, delta_r)`` pair already returned by
    :func:`_validate_bond_set` -- callers validate ONCE at their own entry
    point and pass the result down (no repeated validation per pair)."""
    _, nmat, nvol, nd, _ = prepped.green_rt.shape
    n_channels, delta_r = bond_set_validated

    # nblock == 1 (validated by the caller): drop the block axis.
    green_rt = prepped.green_rt[0]        # (nmat, nvol, nd, nd)
    green_rev = prepped.green_rev[0]      # (nmat, nvol, nd, nd)
    sgn = prepped.sgn.reshape((nmat, 1, 1, 1))
    # Baked in ONCE (sgn is +-1, exact, so this commutes with the
    # contraction regardless of which factor carries it) -- mirrors
    # bond_bubble's G_fwd_sgn, computed outside the (m, m') loop.
    green_fwd_sgn = green_rt * sgn
    # green_rt is consumed: never read again below (every pair rolls a
    # fresh copy of green_rev instead). Release BOTH the local view and
    # prepped's own reference (a bare `del green_rt` would leave
    # prepped.green_rt holding the underlying buffer alive) -- this
    # generator is the sole consumer of a freshly built _DensePrepared
    # instance per bond_bubble_static/_iter_bond_dynamic call, so clearing
    # the attribute here cannot affect any other caller.
    del green_rt
    prepped.green_rt = None

    tail_pack = None
    if prepped.tail_on:
        tail_pack = (prepped.fwd0_p[0], prepped.rev0_p[0],
                     prepped.jump_f[0], prepped.jump_r_rev[0])

    for m in range(n_channels):
        dm = delta_r[m]
        for mp in range(n_channels):
            dmp = delta_r[mp]
            shift = (int(dm[0]) - int(dmp[0]),
                     int(dm[1]) - int(dmp[1]),
                     int(dm[2]) - int(dmp[2]))

            block = _bond_pair_full_block(
                green_fwd_sgn, green_rev, tail_pack, beta, spatial_shape,
                workers, shift)
            yield (m, mp), block


def _assemble_bond_static(prepped, beta, bond_set_validated, spatial_shape,
                          workers):
    """Stream the bond-enlarged static bubble
    ``chi_bar(q; Delta r_m, Delta r_m')`` one ``(m, m')`` channel pair at a
    time, reading only the ``Omega=0`` slice (bosonic Matsubara index
    ``nmat // 2``) of each pair's full block from
    :func:`_iter_bond_channel_pairs` -- this keeps the working set at
    O(per-pair general bubble) rather than O(nmat * ND**2) (spec: "The
    bond assembly"). ``bond_set_validated`` is the already-validated
    ``(n_channels, delta_r)`` pair (see :func:`_iter_bond_channel_pairs`)."""
    xp = _bk.array_module_of(prepped.green_rt)
    _, nmat, nvol, nd, _ = prepped.green_rt.shape
    npair = nd * nd
    n_channels, _ = bond_set_validated
    nD = npair * n_channels
    static_index = nmat // 2

    chi_bar = xp.zeros((nvol, nD, nD), dtype=xp.complex128)

    for (m, mp), block in _iter_bond_channel_pairs(
            prepped, bond_set_validated, beta, spatial_shape, workers):
        chi_bar[:, m * npair:(m + 1) * npair,
                mp * npair:(mp + 1) * npair] = block[static_index]
        del block

    return chi_bar


def _iter_bond_dynamic(green_kw, green0_tail, beta, bond_set, *,
                       spatial_shape, workers=None):
    """Bond-enlarged DYNAMIC (all bosonic Matsubara frequencies)
    particle-hole bubble, streamed one ``(m, m')`` channel pair at a time.

    Parameters
    ----------
    green_kw, green0_tail, beta, bond_set, spatial_shape, workers
        Identical semantics to :func:`bond_bubble_static` (same
        validation, same Green/tail contracts, same endpoint-mean
        correction under a spatial shift). Validation runs EAGERLY on
        entry -- this function is a plain function (not a generator
        function) that validates and prepares, then returns a generator
        for the per-pair iteration, so a caller that never advances the
        returned generator still gets immediate ``ValueError`` feedback on
        a bad input.

    Yields
    ------
    ((m, mp), block)
        ``(m, mp)`` in row-major channel order (``m`` outer, ``mp`` inner,
        both ``0 .. bond_set.n_channels - 1``). ``block`` is
        ``(nmat, nvol, npair, npair)`` complex128, FRESHLY ALLOCATED per
        pair (never a reused buffer). Frequency axis 0 is bosonic
        Matsubara index ``l``, related to the bosonic Matsubara index by
        ``2*l - nmat`` (zero at ``l = nmat // 2`` -- pinned to equal
        :func:`bond_bubble_static`'s output at that same channel pair, up
        to floating-point round-off, since both read the SAME per-pair
        pipeline).

    Raises
    ------
    ValueError
        Same conditions as :func:`bond_bubble_static` (shape/dtype/
        nblock/bond_set mismatches, and an odd ``nmat`` -- unlike the
        plain dense path, the bond entry points still require an even
        ``nmat``, see :func:`_validate_even_nmat_for_bond`); see
        :func:`_validate_dense_inputs` / :func:`_validate_bond_set`.
    """
    (green_kw, green0_tail, beta, nblock, nmat, nvol, nd, spatial_shape
     ) = _validate_dense_inputs(green_kw, green0_tail, beta, spatial_shape,
                                 "general")
    if nblock != 1:
        raise ValueError(
            "bubble kernel: _iter_bond_dynamic requires green_kw's block "
            "axis to be 1 (bond entry points take no scheme argument and "
            "no nblock=2 cross-block semantics are defined), got "
            "nblock={}".format(nblock))
    _validate_even_nmat_for_bond(nmat)
    bond_set_validated = _validate_bond_set(bond_set)

    prepped = _prepare_dense(green_kw, green0_tail, beta, spatial_shape,
                             workers)

    return _iter_bond_channel_pairs(prepped, bond_set_validated, beta,
                                    spatial_shape, workers)


def bond_bubble_static(green_kw, green0_tail, beta, bond_set, *,
                       spatial_shape, workers=None):
    """Bond-enlarged static (``Omega=0``) particle-hole bubble
    ``chi_bar(q; Delta r_m, Delta r_m')``.

    Parameters
    ----------
    green_kw : ndarray, complex, shape (1, nmat, nvol, p, p)
        Green's function in k-space and fermionic Matsubara frequency,
        canonical layout with a single propagator block (``nblock == 1``
        strictly -- see ``ValueError`` below). Same deflated/full
        Green/tail semantics as :func:`dense_bubble` (module docstring's
        "Green/tail contracts").
    green0_tail : ndarray or None, same shape/dtype family as ``green_kw``
        The paired tau-space tail add-back, or None to disable the tail
        correction. Same pairing-provenance caveat as :func:`dense_bubble`.
    beta : float
        Inverse temperature; must be > 0.
    bond_set : ResolvedInteractionSet-like
        Bond-channel topology (``delta_r``, ``n_channels``), e.g. from
        ``bond_channels.resolve_interactions``.
    spatial_shape : tuple of 3 positive ints (Nx, Ny, Nz)
        The lattice shape; ``prod(spatial_shape)`` must equal
        ``green_kw.shape[2]`` (``nvol``).
    workers : int or None, optional
        Forwarded to ``backend.spatial_ifftn``/``spatial_fftn``.

    Returns
    -------
    ndarray, complex128, shape (nvol, ND, ND)
        ``ND = npair * bond_set.n_channels`` with ``npair = p**2``; the
        static (``Omega=0``) bond-enlarged bubble, bond-major/orbital-minor
        index ``I = m * npair + (a * p + c)`` (spec's pair-block mapping,
        numerically identical to ``bond_channels.bond_bubble``'s
        ``(l1*norb+l2)`` convention).

    Raises
    ------
    ValueError
        On any shape/dtype/nblock/bond_set mismatch (see
        :func:`_validate_dense_inputs` / :func:`_validate_bond_set`); these
        checks use ``ValueError`` rather than a bare ``assert`` so they
        survive ``python -O``. In particular ``green_kw``'s block axis
        must be exactly 1 -- the bond entry points take no ``scheme``
        argument (the enlarged object requires the general orbital
        product) and no ``nblock == 2`` cross-block semantics are defined
        (spec: "Scheme x bond and the pair-block mapping") -- and
        ``nmat`` must be even (unlike the plain dense path): the static
        ``Omega=0`` slice is identified by ``nmat // 2``, which assumes
        the even centered Matsubara grid (see
        :func:`_validate_even_nmat_for_bond`).
    """
    (green_kw, green0_tail, beta, nblock, nmat, nvol, nd, spatial_shape
     ) = _validate_dense_inputs(green_kw, green0_tail, beta, spatial_shape,
                                 "general")
    if nblock != 1:
        raise ValueError(
            "bubble kernel: bond_bubble_static requires green_kw's block "
            "axis to be 1 (bond entry points take no scheme argument and "
            "no nblock=2 cross-block semantics are defined), got "
            "nblock={}".format(nblock))
    _validate_even_nmat_for_bond(nmat)
    bond_set_validated = _validate_bond_set(bond_set)

    prepped = _prepare_dense(green_kw, green0_tail, beta, spatial_shape,
                             workers)
    return _assemble_bond_static(prepped, beta, bond_set_validated,
                                 spatial_shape, workers)
