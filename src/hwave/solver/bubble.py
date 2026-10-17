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
``RPA._calc_chi0q_transverse`` (rpa.py, around line 3016 at the start of
the 2026-08-15-transverse-bubble-on-kernel series),
``FLEX._calc_chi0q_ir`` / ``FLEX._calc_chi0q_general_ir`` (flex.py, around
lines 818-917 at series start), and ``bond_channels.bond_bubble``
(bond_channels.py, around lines 1171-1315 at series start) without
altering them -- including the equal-time endpoint-mean correction (issue
#134) and the data-driven tail gate on the dense path. ``bubble.py``
deliberately imports neither ``sparse_ir`` nor ``flex.py`` -- the IR entry
points receive already-constructed ``IRAxis`` objects and duck-type their
attributes/methods (``statistics``, ``beta``, ``wmax``, ``eps``, ``tau``,
``n_tau``, ``n_freq``, ``freq_to_tau_points``, ``tau_to_freq``) -- so the
plain dense path stays usable without the optional IR dependency; the
LONGITUDINAL bond entry points (``bond_bubble_static``, ``_iter_bond_dynamic``)
take a ``bond_set`` (duck-typed ``delta_r``/``n_channels``) and never import
``bond_channels`` either.

The one exception is ``transverse_bond_bubble_static`` (Phase W Task 4,
``docs/superpowers/specs/2026-08-15-bond-transverse-design.md``): its
``topo`` argument is a ``bond_channels.TransverseTopology``, and the
spec requires ``validate_topology_against_mesh`` -- defined in
``bond_channels.py`` -- to run at entry. Since ``bond_channels.py``
imports THIS module at its own top level (``from . import bubble as
_bubble``), a module-level ``from . import bond_channels`` here would
cycle; ``transverse_bond_bubble_static`` therefore imports
``bond_channels`` LOCALLY, inside the function, instead.
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
# N4_BUFFERS = 3, matching the legacy ``bond_bubble``'s count: the
# ``contract_general``-then-reshape copy in :func:`_bond_pair_full_block`
# only ever holds 2 N4-sized buffers at once (below ``tau_to_boson``'s
# internal peak of 3, so it is not the binding constraint). The REAL risk
# to this budget is generator retention across ``yield``:
# :func:`_iter_bond_channel_pairs` is a generator, and a generator's local
# variables -- including ``block``, the value it just yielded -- stay bound
# in its suspended frame until it is explicitly cleared, REGARDLESS of
# whether the consumer drops its own reference (``_assemble_bond_static``'s
# ``del block`` only releases the CONSUMER's reference). Without an
# explicit ``del`` on the generator's side, the previous pair's block stays
# alive while the next pair's is being built, adding a 4th simultaneous
# N4-sized array. :func:`_iter_bond_channel_pairs` clears its own
# reference immediately after each ``yield`` to close this.
BOND_BUBBLE_N4_BUFFERS = 3
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
    if isinstance(n_channels, bool) or not isinstance(
            n_channels, (int, np.integer)) or n_channels <= 0:
        raise ValueError(
            "bubble kernel: bond_set.n_channels must be a positive integer, "
            "got {!r} -- a malformed (e.g. zero-channel) bond set would "
            "otherwise pass through and yield a plausible-looking but "
            "empty/zero-sized susceptibility".format(n_channels))
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


def _assemble_cross_block(prepped, fwd_block, rev_block, scheme, beta,
                          spatial_shape, workers):
    """Cross-block contraction for the transverse bubble: contract
    ``prepped.green_rt[fwd_block]`` (the forward factor) against
    ``prepped.green_rev[rev_block]`` (the reversed factor) via the SAME
    ``contract_reduced``/``contract_general`` primitives and the SAME
    per-tau ``sgn`` / equal-time endpoint-mean correction / spatial
    ``fftn`` / ``tau_to_boson`` / final ``x(-1/beta)`` transport as
    :func:`_assemble_plain`'s dense arm -- the only differences are that
    the two factors are read from DIFFERENT block indices
    (``fwd_block``/``rev_block``) rather than the same one, and the
    block axis is consumed by the crossing rather than carried through
    to the output. Transcribed from the legacy ``RPA._calc_chi0q_transverse``
    body (deleted at the end of the migration series; see git history)
    verbatim.

    The endpoint-mean correction threads ``fwd_block``/``rev_block``
    through EACH of its four independently-indexable terms
    (``fwd0_p[fwd_block]``, ``rev0_p[rev_block]``, ``jump_f[fwd_block]``,
    ``jump_r_rev[rev_block]``) -- the forward branch and its jump always
    read ``fwd_block``, the reversed branch and its jump always read
    ``rev_block``, never the same index on both sides and never a mixed
    pair beyond this pattern (spec: "The kernel entry")::

        chi0_endpoint = 0.5 * ( contract(fwd0_p[fwd], rev0_p[rev] + jump_r_rev[rev])
                              + contract(fwd0_p[fwd] + jump_f[fwd], rev0_p[rev]) )

    ``prepped`` is a :class:`_DensePrepared` instance (from
    :func:`_prepare_dense` on a two-block Green function); this function
    performs no validation of its own -- callers (:func:`transverse_bubble`)
    validate ``nblock == 2`` before building ``prepped``."""
    nx, ny, nz = spatial_shape
    _, nmat, nvol, nd, _ = prepped.green_rt.shape

    contract = contract_reduced if scheme == "reduced" else contract_general
    nd_shape = (nd, nd) if scheme == "reduced" else (nd, nd, nd, nd)
    nds = nd * nd if scheme == "reduced" else nd ** 4

    g_fwd = prepped.green_rt[fwd_block]
    g_rev = prepped.green_rev[rev_block]

    chi0_rt = contract(g_fwd, g_rev)
    sgn_bc = prepped.sgn.reshape((nmat,) + (1,) * (chi0_rt.ndim - 1))
    chi0_rt = chi0_rt * sgn_bc

    if prepped.tail_on:
        # mean of the two equal-time branches (sgn[0] = +1)
        chi0_rt[0] = 0.5 * (
            contract(prepped.fwd0_p[fwd_block],
                     prepped.rev0_p[rev_block] + prepped.jump_r_rev[rev_block])
            + contract(prepped.fwd0_p[fwd_block] + prepped.jump_f[fwd_block],
                       prepped.rev0_p[rev_block]))

    chi0_qt = _bk.spatial_fftn(
        chi0_rt.reshape(nmat, nx, ny, nz, nds),
        axes=(1, 2, 3), workers=workers)

    chi0_qt_flat = chi0_qt.reshape(nmat, nvol * nds)
    chi0_qw = _ms.tau_to_boson(chi0_qt_flat, axis=0).reshape(
        nmat, nvol, *nd_shape) * (-1.0 / beta)

    return chi0_qw


def transverse_bubble(green_kw, green0_tail, beta, *, spatial_shape,
                      scheme, workers=None):
    """Cross-block particle-hole bubble chi0_+-: contract block 0
    (forward) with block 1 (reversed).

    ``green_kw``/``green0_tail``: canonical ``(2, nmat, nvol, p, p)`` --
    ``nblock == 2`` REQUIRED (``ValueError`` otherwise, naming both the
    actual ``nblock`` and this contract). Green/tail semantics identical
    to :func:`dense_bubble` (deflated pair; ``None`` => full G, tail
    machinery off; conditional validation). ``scheme`` in
    ``{"reduced", "general"}``.

    This docstring states the SPIN-DIAG two-block contract only -- block
    0 is G_up, block 1 is G_down, and this entry point computes their
    cross term chi0_+-. The full spinful (spin-mixing) transverse bubble
    generalizes this -- that is a SEPARATE later spec/campaign (Step 3b),
    not a task within this series; it is not what this function computes.

    Returns
    -------
    ndarray, complex128
        ``(nmat, nvol, p, p)`` reduced / ``(nmat, nvol, p, p, p, p)``
        general -- NO block axis (the block dimension is consumed by the
        up/down crossing), bosonic uniform frequency axis 0 with the
        same ``l <-> 2l - nmat`` mapping as :func:`dense_bubble`. Sign
        and normalization identical to the legacy transverse body:
        per-tau ``sgn``, endpoint-mean correction, final ``x(-1/beta)``.

    Raises
    ------
    ValueError
        If ``green_kw``'s block axis is not exactly 2, or on any other
        shape/dtype/scheme mismatch (see :func:`_validate_dense_inputs`);
        these checks use ``ValueError`` rather than a bare ``assert`` so
        they survive ``python -O``.
    """
    if green_kw.ndim != 5:
        raise ValueError(
            "bubble kernel: green_kw must have ndim 5 "
            "(nblock, nmat, nvol, p, p), got ndim={} shape={}".format(
                green_kw.ndim, green_kw.shape))
    nblock = green_kw.shape[0]
    if nblock != 2:
        # Checked FIRST, before _validate_dense_inputs's own (looser,
        # nblock in {1, 2}) guard -- transverse_bubble's up/down crossing
        # requires exactly 2 blocks, a stricter contract than the generic
        # dense entry point's, so both nblock=1 and nblock=3 must surface
        # THIS message (naming the actual nblock and the required value)
        # rather than the generic one.
        raise ValueError(
            "bubble kernel: transverse_bubble requires green_kw's block "
            "axis to be exactly 2 (a spin-diag Green's function with "
            "G_up as block 0 and G_down as block 1 -- the up/down "
            "cross-block contraction this entry point computes has no "
            "meaning for any other block count), got nblock={} "
            "(required: 2)".format(nblock))

    (green_kw, green0_tail, beta, nblock, nmat, nvol, nd, spatial_shape
     ) = _validate_dense_inputs(green_kw, green0_tail, beta, spatial_shape,
                                scheme)

    prepped = _prepare_dense(green_kw, green0_tail, beta, spatial_shape,
                             workers)
    return _assemble_cross_block(prepped, 0, 1, scheme, beta,
                                 spatial_shape, workers)


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
            # A generator's locals -- including `block`, just yielded --
            # stay bound in its suspended frame across the yield regardless
            # of what the consumer does with ITS reference, so without this
            # the previous pair's block stays alive while the next pair's
            # is being built (a 4th N4-sized buffer beyond the documented
            # budget). Runs when the generator is resumed for the next
            # pair, i.e. AFTER the consumer's own `del block` (if any) has
            # already dropped its reference -- so this is the point the
            # buffer's refcount actually reaches zero.
            del block


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


# ===========================================================================
# bond-enlarged bubble  --  private IR x bond cell
# ===========================================================================
#
# _ir_bond_static computes the SAME static (Omega=0) bond-enlarged bubble as
# bond_bubble_static, natively on the sparse-IR bosonic node l=0 rather than
# via dense-Matsubara-grid truncation, reusing _prepare_ir (the same
# fermionic-node-to-bosonic-tau transport, spatial ifftn, and reversed
# propagator ir_bubble already uses) and the same per-pair
# roll -> contract_general -> pair-reshape -> spatial fftn pipeline as the
# dense bond assembly, with tau_to_freq_points(block, [0]) in place of the
# dense path's spatial-fftn -> tau_to_boson -> (-1/beta) tail (spec:
# "_ir_bond_static assembly sequence").
#
# Sign convention: _prepare_ir's g_rev already carries ONE explicit minus
# (the fermionic-antiperiodicity sign, spec: "IR Omega=0 semantics"). The
# per-pair contraction below applies a SECOND explicit minus at
# contract_general -- the closed-fermion-loop bubble sign -- exactly
# mirroring ir_bubble's own IR arm (_assemble_plain's
# ``chi0_rt = -contract(...)``, the direct counterpart of dense_bubble's
# ``sgn`` array plus its final ``-1/beta``, both of which also carry two
# independent sign contributions). Both signs are REQUIRED: at zero bond
# displacement (single channel, delta_r=(0,0,0)) this cell's diagonal block
# must equal ir_bubble's general-scheme Omega=0 slice EXACTLY (mirroring the
# dense bond assembly's own zero-displacement identity with dense_bubble's
# general scheme, verified bitwise during development) -- omitting the
# second minus produces the exact negative of that value, not a
# discretization-scale discrepancy, so TestIrBondOracle/TestIrBondVsDense
# would fail outright (not just miss a tolerance) rather than merely
# degrade, making this identity self-checking.


def _validate_ir_bond_inputs(green_kw, ax_fermi, ax_bose, spatial_shape):
    """Validate an ``_ir_bond_static`` call: IR axis compatibility
    (:func:`_validate_ir_axes`) plus everything :func:`_validate_ir_inputs`
    validates for the general-orbital contraction, MINUS the ``scheme``
    argument itself (the bond entry point takes none, mirroring
    :func:`bond_bubble_static`'s dense counterpart) -- and ``nblock`` must
    be exactly 1 (no ``scheme`` argument and no ``nblock == 2`` cross-block
    semantics are defined for the bond entry points). Returns the
    normalized ``(green_kw, nvol, nd, spatial_shape)``. ``ValueError`` on
    any failure (survives ``python -O``)."""
    _validate_ir_axes(ax_fermi, ax_bose)

    xp = _bk.array_module_of(green_kw)

    if green_kw.ndim != 5:
        raise ValueError(
            "bubble kernel: green_kw must have ndim 5 "
            "(nblock, n_freq_F, nvol, p, p), got ndim={} shape={}".format(
                green_kw.ndim, green_kw.shape))
    nblock, nw, nvol, nd, nd2 = green_kw.shape
    if nblock != 1:
        raise ValueError(
            "bubble kernel: _ir_bond_static requires green_kw's block "
            "axis to be 1 (bond entry points take no scheme argument and "
            "no nblock=2 cross-block semantics are defined), got "
            "nblock={}".format(nblock))
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

    return (green_kw, nvol, nd, spatial_shape)


def _ir_bond_pair_block(g_fwd, g_rev, ax_bose, spatial_shape, workers,
                        shift, freq_indices):
    """Compute the IR bond bubble block for ONE ``(m, mp)`` channel pair at
    the requested bosonic IR frequency indices: ``(nvol, npair, npair,
    len(freq_indices))`` complex128, FRESHLY ALLOCATED (mirrors
    :func:`_bond_pair_full_block`'s per-pair pipeline through the IR
    transport -- spec: "the IR path mirrors its structure but through the
    IR transport").

    Roll the reversed Green by ``shift`` (``Delta r_m - Delta r_m'``),
    ``contract_general`` (with the closed-loop sign, see the module
    comment above), reshape the trailing ``(a, c, b, d)`` axes into the
    ``(npair, npair)`` pair block and split ``nvol`` back to
    ``(Nx, Ny, Nz)`` in one reshape (identical to the dense pipeline's
    reshape), spatial ``fftn``, then ``ax_bose.tau_to_freq_points`` at the
    requested indices -- every per-pair temporary is freed before return,
    mirroring the dense path's ``del`` discipline.

    ``g_fwd``/``g_rev``: ``(n_tau_B, nvol, nd, nd)`` (block axis already
    dropped by the caller)."""
    nt, nvol, nd, _ = g_rev.shape
    nx, ny, nz = spatial_shape
    npair = nd * nd

    g_rev_shifted = _roll_spatial(g_rev, 1, spatial_shape, shift)

    # (n_tau_B, nvol, a, c, b, d) -- the closed-fermion-loop bubble sign
    # (see the module comment above: this is the SECOND of the two signs
    # the IR bond cell needs, on top of the one already embedded in g_rev
    # by _prepare_ir).
    chi0_rt = -contract_general(g_fwd, g_rev_shifted)
    del g_rev_shifted

    # pair-block reshape + split nvol back to (Nx, Ny, Nz), one reshape
    # (contract_general's trailing axis order (a, c, b, d) already IS the
    # row/column split -- identical to the dense pipeline).
    chi0_rt = chi0_rt.reshape(nt, nx, ny, nz, npair, npair)

    chi0_qt = _bk.spatial_fftn(chi0_rt, axes=(1, 2, 3), workers=workers)
    del chi0_rt
    chi0_qt = chi0_qt.reshape(nt, nvol * npair * npair)

    xp = _bk.array_module_of(chi0_qt)
    chi0_qt_flat = xp.moveaxis(chi0_qt, 0, -1)   # (nvol*npair*npair, n_tau_B)
    del chi0_qt

    val = ax_bose.tau_to_freq_points(chi0_qt_flat, freq_indices)
    del chi0_qt_flat

    block = val.reshape(nvol, npair, npair, freq_indices.size)
    del val

    return block


def _ir_bond_static(green_kw, ax_fermi, ax_bose, bond_set, *,
                    spatial_shape, workers=None):
    """Bond-enlarged static (bosonic Matsubara index ``n = 0``) IR
    particle-hole bubble ``chi_bar(q; Delta r_m, Delta r_m')``, computed
    natively on the sparse-IR bosonic axis rather than by dense-Matsubara
    truncation.

    Parameters
    ----------
    green_kw : ndarray, complex, shape (1, n_freq_F, nvol, p, p)
        The FULL Green's function (no tail subtracted) sampled at
        ``ax_fermi``'s fermionic Matsubara nodes -- no ``green0_tail``
        argument (same Green/tail contract as :func:`ir_bubble`; the IR
        basis already represents the ``coeff_tail / (i w)`` asymptotics
        within ``ax_fermi.eps``). ``nblock`` must be exactly 1.
    ax_fermi, ax_bose : IRAxis-like
        Already-constructed fermionic/bosonic IR axis objects; must be
        parameter-compatible (see :func:`_validate_ir_axes`).
    bond_set : ResolvedInteractionSet-like
        Bond-channel topology (``delta_r``, ``n_channels``).
    spatial_shape : tuple of 3 positive ints (Nx, Ny, Nz)
        The lattice shape; ``prod(spatial_shape)`` must equal
        ``green_kw.shape[2]`` (``nvol``).
    workers : int or None, optional
        Forwarded to ``backend.spatial_ifftn``/``spatial_fftn``.

    Returns
    -------
    ndarray, complex128, shape (nvol, ND, ND)
        ``ND = npair * bond_set.n_channels`` with ``npair = p**2``; the
        static (``iOmega_0 = 0``) bond-enlarged bubble, same
        bond-major/orbital-minor index convention as
        :func:`bond_bubble_static`.

    Raises
    ------
    ValueError
        On any shape/dtype/nblock/axis-compatibility/bond_set mismatch
        (see :func:`_validate_ir_bond_inputs` / :func:`_validate_bond_set`);
        these checks use ``ValueError`` rather than a bare ``assert`` so
        they survive ``python -O``.
    """
    (green_kw, nvol, nd, spatial_shape
     ) = _validate_ir_bond_inputs(green_kw, ax_fermi, ax_bose, spatial_shape)
    n_channels, delta_r = _validate_bond_set(bond_set)

    prepped = _prepare_ir(green_kw, ax_fermi, ax_bose, spatial_shape,
                          workers)
    # nblock == 1 (validated above): drop the block axis. No forward-side
    # sgn precompute here (unlike the dense path's green_fwd_sgn) -- the IR
    # transport has no per-tau sgn array (spec: "no sgn array"), so the
    # SAME g_rt is reused unmodified for every channel pair below.
    g_rt = prepped.green_rt[0]
    g_rev = prepped.green_rev[0]

    npair = nd * nd
    nD = npair * n_channels
    xp = _bk.array_module_of(g_rt)
    chi_bar = xp.zeros((nvol, nD, nD), dtype=xp.complex128)
    freq0 = np.array([0], dtype=np.int64)

    for m in range(n_channels):
        dm = delta_r[m]
        for mp in range(n_channels):
            dmp = delta_r[mp]
            shift = (int(dm[0]) - int(dmp[0]),
                     int(dm[1]) - int(dmp[1]),
                     int(dm[2]) - int(dmp[2]))

            block = _ir_bond_pair_block(g_rt, g_rev, ax_bose, spatial_shape,
                                        workers, shift, freq0)
            chi_bar[:, m * npair:(m + 1) * npair,
                    mp * npair:(mp + 1) * npair] = block[..., 0]
            del block

    return chi_bar


# ===========================================================================
# bond-enlarged bubble  --  static cross-block (transverse) assembly
# ===========================================================================
#
# transverse_bond_bubble_static computes the STATIC (Omega=0) bond-enlarged
# transverse particle-hole bubble chi0_+-(q; Delta r_m, Delta r_m'), the
# cross-block analogue of bond_bubble_static: instead of contracting one
# propagator block against itself channel-pair by channel-pair (the
# longitudinal bond path above), it contracts block 1 (DOWN, forward) against
# block 0 (UP, reversed) -- the SAME per-pair roll-by-(Delta r_m - Delta r_m')
# / contract_general / endpoint-mean / spatial-fftn / tau_to_boson pipeline
# (_bond_pair_full_block, unmodified) driven by ONE new block-aware generator
# (_iter_transverse_bond_channel_pairs) that plays the role
# _iter_bond_channel_pairs plays for the nblock=1 longitudinal path.
#
# fwd = block 1 (DOWN), rev = block 0 (UP): the AMENDED 2026-08-15 Phase-H
# H3 adjudication (spec, "The bond bubble (static)") -- the Wick contraction
# of the campaign's observable <S+;S+dagger> factorizes as
# -G_down(r,tau)*G_up(-r,-tau), so the forward propagator is the DOWN Green
# function. This is NOT the same role assignment as 3a's
# ``transverse_bubble(fwd=0, rev=1)`` (block 0 = up = forward there): that
# earlier entry point computes the CONJUGATE object <S-;S-dagger>, which
# coincides with <S+;S+dagger> only when G_up == G_down (spec's NOTE under
# "The bond bubble (static)"). The role swap is not absorbable in the frame
# map (a pure index permutation cannot express a spin-role swap): a naive
# fwd=up/rev=down assignment here was measured (Phase-H H3, an Nmat-independent
# effect) to leave an O(1) residual, not a finite-Nmat discretization error --
# see ``transverse_bond_bubble_static``'s docstring for the recorded
# mutation-check magnitude on THIS entry point.


def _iter_transverse_bond_channel_pairs(prepped, delta_r, beta,
                                        spatial_shape, workers):
    """Shared per-pair generator for the transverse (cross-block) bond
    bubble: yields ``((m, mp), block)`` in row-major channel order (``m``
    outer, ``mp`` inner), ``block`` the FULL-frequency ``(nmat, nvol,
    npair, npair)`` complex128 bond bubble for that channel pair, built by
    :func:`_bond_pair_full_block` -- the SAME per-pair pipeline the
    longitudinal path's :func:`_iter_bond_channel_pairs` uses, differing
    only in WHICH block of ``prepped`` supplies the forward/reversed
    factors: forward = block 1 (DOWN), reversed = block 0 (UP) (see the
    module comment above and ``transverse_bond_bubble_static``'s
    docstring). ``prepped`` is assumed built from a two-block
    (``nblock == 2``) Green function; this function performs no
    validation of its own -- callers (``transverse_bond_bubble_static``)
    validate ``nblock == 2`` before building ``prepped``.

    ``delta_r``: ``(B, 3)`` integer array (``topo.delta_r``, already
    validated by the caller).

    Mirrors :func:`_iter_bond_channel_pairs`'s del discipline: the
    unused halves of ``prepped.green_rt``/``prepped.green_rev`` (block 0
    of ``green_rt`` -- UP, never read as the forward factor here -- and
    block 1 of ``green_rev`` -- DOWN, never read as the reversed factor
    here) are released as soon as the needed block is extracted, rather
    than staying bound in ``prepped`` for the whole ``(m, mp)`` loop; and
    the generator clears its own ``block`` reference immediately after
    each ``yield`` (see :func:`_iter_bond_channel_pairs`'s comment on why
    a generator's suspended-frame locals need an explicit ``del`` beyond
    whatever the consumer does with its own reference)."""
    _, nmat, nvol, nd, _ = prepped.green_rt.shape
    B = delta_r.shape[0]

    sgn = prepped.sgn.reshape((nmat, 1, 1, 1))
    # Forward = block 1 (DOWN). `* sgn` is an elementwise product -> a
    # fresh array, decoupled from the two-block buffer underneath, so
    # dropping prepped.green_rt right after this line is safe.
    green_rt_full = prepped.green_rt
    green_fwd_sgn = green_rt_full[1] * sgn
    del green_rt_full
    prepped.green_rt = None

    # Reversed = block 0 (UP). Plain indexing is a VIEW into the two-block
    # buffer, so this is copied explicitly before the full buffer (with
    # its unused block 1) is released.
    green_rev_full = prepped.green_rev
    green_rev = green_rev_full[0].copy()
    del green_rev_full
    prepped.green_rev = None

    tail_pack = None
    if prepped.tail_on:
        # Endpoint threading (spec: "The bond bubble (static)" / mirrors
        # _assemble_cross_block's own threading): forward branch and its
        # jump always read block 1 (down), reversed branch and its jump
        # always read block 0 (up).
        tail_pack = (prepped.fwd0_p[1], prepped.rev0_p[0],
                     prepped.jump_f[1], prepped.jump_r_rev[0])

    for m in range(B):
        dm = delta_r[m]
        for mp in range(B):
            dmp = delta_r[mp]
            shift = (int(dm[0]) - int(dmp[0]),
                     int(dm[1]) - int(dmp[1]),
                     int(dm[2]) - int(dmp[2]))

            block = _bond_pair_full_block(
                green_fwd_sgn, green_rev, tail_pack, beta, spatial_shape,
                workers, shift)
            yield (m, mp), block
            # See _iter_bond_channel_pairs's identical comment: releases
            # the previous pair's block from this generator's own
            # suspended frame, keeping the per-pair working set bounded.
            del block


def transverse_bond_bubble_static(green_kw, green0_tail, beta, topo, *,
                                  spatial_shape, workers=None):
    """Bond-enlarged static (``Omega=0``) TRANSVERSE particle-hole bubble
    ``chi0_+-(q; Delta r_m, Delta r_m')``, the cross-block analogue of
    :func:`bond_bubble_static`.

    **Block roles (AMENDED 2026-08-15, Phase-H H3 adjudication -- spec
    "The bond bubble (static)"): fwd = block 1 (DOWN), rev = block 0
    (UP).** The Wick contraction of the observable ``<S+;S+dagger>``
    factorizes as ``-G_down(r,tau)*G_up(-r,-tau)``, i.e. the forward
    propagator is the DOWN Green function -- independently derived by the
    H3 reviewer and pinned by the V=0 gate (an un-swapped ``fwd=up,
    rev=down`` choice leaves an Nmat-INDEPENDENT residual, not a
    finite-Nmat discretization error; see the mutation-check note below).
    This is NOT the same convention as :func:`transverse_bubble`'s own
    ``fwd=block 0 (up), rev=block 1 (down)``: that entry point computes
    the CONJUGATE object ``<S-;S-dagger>``, which coincides with
    ``<S+;S+dagger>`` only when ``G_up == G_down`` (an open question this
    spec explicitly leaves to a later granule, not settled by this
    function). The role swap is NOT absorbable into the roll/frame-map
    machinery (a pure index permutation cannot express a spin-role swap).

    Parameters
    ----------
    green_kw : ndarray, complex, shape (2, nmat, nvol, p, p)
        Canonical two-block Green's function, ``nblock == 2`` REQUIRED
        (block 0 = G_up, block 1 = G_down -- ``ValueError`` naming both
        the actual ``nblock`` and the required value otherwise). Same
        deflated/full Green/tail semantics as :func:`dense_bubble`
        (module docstring's "Green/tail contracts"; deflated pair when
        ``green0_tail`` is given, ``None`` => full G and the tail
        machinery off).
    green0_tail : ndarray or None, same shape/dtype family as ``green_kw``
        The paired tau-space tail add-back, or ``None`` to disable the
        tail correction.
    beta : float
        Inverse temperature; must be > 0.
    topo : bond_channels.TransverseTopology
        The master transverse bond topology (``delta_r``/``reverse``/
        ``coeffs`` -- only ``delta_r`` is read here; ``coeffs`` belongs to
        the VERTEX, ``W_pm_bond``, not the bubble). Validated at entry via
        ``bond_channels.validate_topology_against_mesh`` (the ONE
        mesh-dependent lifecycle point shared with ``W_pm_bond`` -- spec,
        "The master transverse topology").
    spatial_shape : tuple of 3 positive ints (Nx, Ny, Nz)
        The lattice shape; ``prod(spatial_shape)`` must equal
        ``green_kw.shape[2]`` (``nvol``), and ``topo.delta_r`` must be
        INJECTIVE modulo this mesh (``validate_topology_against_mesh``).
    workers : int or None, optional
        Forwarded to ``backend.spatial_ifftn``/``spatial_fftn``.

    Returns
    -------
    ndarray, complex128, shape (nvol, ND, ND)
        ``ND = npair * B`` with ``npair = p**2``, ``B = len(topo.delta_r)``;
        the static (``Omega=0``) bond-enlarged transverse bubble,
        bond-major/orbital-minor index ``I = m * npair + (a * p + c)`` --
        the SAME index convention :func:`bond_bubble_static` uses. The
        ``Omega=0`` element is bosonic Matsubara index ``nmat // 2`` (even
        ``nmat`` REQUIRED, same as :func:`bond_bubble_static` -- see
        :func:`_validate_even_nmat_for_bond`).

    Raises
    ------
    ValueError
        If ``green_kw``'s block axis is not exactly 2 (naming the actual
        ``nblock`` and the required value 2); on any other shape/dtype
        mismatch (see :func:`_validate_dense_inputs`); if ``nmat`` is odd
        (see :func:`_validate_even_nmat_for_bond`); or if ``topo``/
        ``spatial_shape`` fail mesh-injectivity or ``topo``'s own
        constructor invariants (propagated from
        ``bond_channels.validate_topology_against_mesh``). All failures
        are ``ValueError`` (survives ``python -O``).

    Notes
    -----
    Mutation check (fwd/rev swapped to block 0 (up) forward / block 1
    (down) reversed, run manually during development against
    ``tests/test_bond_transverse_w.py``'s V=0 near-machine pin, then
    reverted -- ``git diff`` on this file carries no trace, per the
    established mutation-check discipline in
    ``tests/test_bubble_kernel.py``/``tests/test_bond_transverse_ed.py``'s
    H3): with the roles swapped, the pin's max-abs distance jumps from
    round-off (~1e-13) to a MEASURED 2.465e-2 (norb=1, L=5 fixture) /
    2.959e-2 (norb=2, L=3 fixture) -- an O(1) margin (the same ~2.5e-2
    class H3's own mutation check measures), confirming the role
    assignment above is load-bearing and not absorbable elsewhere.
    """
    if green_kw.ndim != 5:
        raise ValueError(
            "bubble kernel: green_kw must have ndim 5 "
            "(nblock, nmat, nvol, p, p), got ndim={} shape={}".format(
                green_kw.ndim, green_kw.shape))
    nblock = green_kw.shape[0]
    if nblock != 2:
        # Checked FIRST, before _validate_dense_inputs's own (looser,
        # nblock in {1, 2}) guard -- mirrors transverse_bubble's identical
        # early check: the cross-block down/up contraction requires
        # exactly 2 blocks, so both nblock=1 and nblock=3+ must surface
        # THIS message (naming the actual nblock and the required value)
        # rather than the generic one.
        raise ValueError(
            "bubble kernel: transverse_bond_bubble_static requires "
            "green_kw's block axis to be exactly 2 (a spin-diag Green's "
            "function with G_up as block 0 and G_down as block 1 -- the "
            "down/up cross-block contraction this entry point computes "
            "has no meaning for any other block count), got nblock={} "
            "(required: 2)".format(nblock))

    (green_kw, green0_tail, beta, nblock, nmat, nvol, nd, spatial_shape
     ) = _validate_dense_inputs(green_kw, green0_tail, beta, spatial_shape,
                                "general")
    _validate_even_nmat_for_bond(nmat)

    # The ONE mesh-dependent validation lifecycle point (spec, "The master
    # transverse topology"), shared with W_pm_bond -- LOCAL import (see the
    # module docstring's import-graph note: bond_channels imports this
    # module at its own top level, so a module-level import here would
    # cycle).
    from . import bond_channels as _bc
    _bc.validate_topology_against_mesh(topo, spatial_shape)
    # Re-fetch a freshly-validated, alias-safe topology for local use (the
    # same defensive re-construction validate_topology_against_mesh applies
    # internally -- topo.coeffs is a plain dict a caller could have
    # re-keyed between that call returning and this line; delta_r/reverse
    # are read-only arrays, but re-running the constructor costs little and
    # keeps this entry point's own invariants independent of that detail).
    topo = _bc.TransverseTopology(topo.delta_r, topo.reverse, topo.coeffs)

    delta_r = np.asarray(topo.delta_r)
    B = delta_r.shape[0]

    prepped = _prepare_dense(green_kw, green0_tail, beta, spatial_shape,
                             workers)

    xp = _bk.array_module_of(green_kw)
    npair = nd * nd
    ND = npair * B
    static_index = nmat // 2

    chi_bar = xp.zeros((nvol, ND, ND), dtype=xp.complex128)

    for (m, mp), block in _iter_transverse_bond_channel_pairs(
            prepped, delta_r, beta, spatial_shape, workers):
        chi_bar[:, m * npair:(m + 1) * npair,
                mp * npair:(mp + 1) * npair] = block[static_index]
        del block

    return chi_bar
