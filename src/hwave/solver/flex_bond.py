"""The dynamic bond pipeline of the Phase B FLEX solver (GitHub issue
#181; spec docs/superpowers/specs/2026-09-06-flex-bond-sigma-phase-b-181-design.md,
section 3): the dense bond-block store, the bubble assembly, the batched
dressing, the bond-aware self-energy transport, the memory preflight and
the output helpers.
"""
import numpy as np

from . import backend as _bk
from . import matsubara as _ms
from . import bubble as _bubble
from . import second_order as _so


class BondBlockStore:
    """Owner of the dense ``(nmat, nvol, ND, ND)`` bond arrays (spec 3.5).

    A context manager: construction allocates every named array (rolling
    back on a partial failure), ``release()`` runs on every exit of the
    ``with`` block. Pair access (``put_pair`` / ``get_pair``: a view
    ``(nmat, nvol, nd, nd)`` onto channel pair ``(alpha, beta)``) and
    frequency-batch access (``put_freq_batch`` / ``get_freq_batch``: a
    contiguous ``(nb, nvol, ND, ND)`` view) are the only ways in and out;
    ``detach(name)`` transfers ownership of one array to the caller.
    A future streaming store implements the same five methods.
    """

    def __init__(self, nmat, nvol, ND, nd, names):
        names = tuple(names)
        if len(set(names)) != len(names):
            raise ValueError("BondBlockStore: duplicate array names {}".format(names))
        if ND % nd != 0:
            raise ValueError("BondBlockStore: ND={} is not a multiple of nd={}".format(ND, nd))
        self.nmat, self.nvol, self.ND, self.nd = int(nmat), int(nvol), int(ND), int(nd)
        self.B = self.ND // self.nd
        self._arrays = {}
        self._released_slots = set()
        try:
            for name in names:
                self._arrays[name] = np.zeros((self.nmat, self.nvol, self.ND, self.ND),
                                              dtype=np.complex128)
        except MemoryError:
            self._arrays.clear()
            raise
        self._released = False

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.release()
        return False

    def _arr(self, name):
        if self._released:
            raise RuntimeError("BondBlockStore: released")
        if name in self._released_slots:
            raise KeyError("BondBlockStore: slot {!r} was released".format(name))
        return self._arrays[name]

    def put_pair(self, name, alpha, beta, block):
        a = self._arr(name)
        nd = self.nd
        a[:, :, alpha * nd:(alpha + 1) * nd, beta * nd:(beta + 1) * nd] = block

    def get_pair(self, name, alpha, beta):
        a = self._arr(name)
        nd = self.nd
        return a[:, :, alpha * nd:(alpha + 1) * nd, beta * nd:(beta + 1) * nd]

    def put_freq_batch(self, name, l0, l1, batch):
        self._arr(name)[l0:l1] = batch

    def get_freq_batch(self, name, l0, l1):
        return self._arr(name)[l0:l1]

    def detach(self, name):
        a = self._arr(name)
        del self._arrays[name]
        return a

    def release_slot(self, name):
        """Free ONE named array (the pairing step drops ``W`` before it
        allocates the IR coefficients, spec 7 step 1). Later access to the
        slot raises KeyError; ``release`` / ``__exit__`` tolerate it."""
        # after release every array is gone, so a valid slot name would look
        # "unknown": say what really happened instead
        if self._released:
            raise RuntimeError("BondBlockStore: released")
        if name not in self._arrays:
            if name in self._released_slots:
                return
            raise KeyError("BondBlockStore: no slot {!r}".format(name))
        del self._arrays[name]
        self._released_slots.add(name)

    def release(self):
        self._arrays.clear()
        self._released = True

    @property
    def released(self):
        return self._released


class BondDeviceContext:
    """Owner of the IMMUTABLE bond inputs for one ``_solve_phase_b`` call
    (spec 2026-09-17 section 4.2): the vertices ``S``, ``C``, the on-site
    ``S_on``, ``C_on``, their sum ``SpC_on`` (formed on the host here, once),
    the mixed-block pair permutation ``perm``, the block-weight ``mask``
    and the bubble's tau-space tail ``green0_tail`` (issue #196: the
    bubble computes on the module of the Green function, which the SCF
    loop already keeps on the device, so its tail belongs to the
    solve-scoped device set rather than to a per-iteration transfer).
    Construction transfers them to the array module ``xp`` exactly once
    (``_bk.to_device``; the identity on numpy, so the CPU path gets the
    host arrays themselves); every SCF iteration reuses them. A context
    manager like :class:`BondBlockStore`: ``release()`` drops the references
    (the cupy pool may then reuse the blocks). ``_solve_phase_b`` creates
    it AFTER the memory preflight (its allocation is part of the predicted
    device need) and nothing else creates device copies of the vertices."""

    _NAMES = ("S", "C", "S_on", "C_on", "SpC_on", "perm", "mask", "green0_tail")

    def __init__(self, xp, S, C, S_on=None, C_on=None, perm=None, mask=None,
                 green0_tail=None):
        self.xp = xp
        # SpC_on exists only when BOTH are present, so one of them alone would
        # be silently dropped into a context whose SpC_on is None
        if (S_on is None) != (C_on is None):
            raise ValueError("BondDeviceContext: S_on and C_on must be given together "
                             "(both None or both arrays); got S_on {}, C_on {}"
                             .format("None" if S_on is None else "an array",
                                     "None" if C_on is None else "an array"))
        SpC_on = None if (S_on is None or C_on is None) else np.asarray(S_on) + np.asarray(C_on)
        host = dict(S=S, C=C, S_on=S_on, C_on=C_on, SpC_on=SpC_on, perm=perm, mask=mask,
                    green0_tail=green0_tail)
        self._arrays = {k: (None if host[k] is None else _bk.to_device(host[k], xp))
                        for k in self._NAMES}
        self._released = False

    @classmethod
    def for_view(cls, xp, S, C, S_on, C_on, view, norb, green0_tail=None, *, perm=None):
        """The context for the bond view ``view`` (issue #198): derives the
        pair permutation from :func:`_mixed_pair_permutation` (looked up on
        the module, so a test may still replace it) and the block-weight
        mask from :func:`mixed_block_mask`, both from ``view.n_channels``
        and ``nd = norb * norb``. ``perm`` overrides the permutation (the
        identity-permutation control of the second-order tests). Refuses a
        vertex whose trailing two dimensions are not ``(n_channels * nd,
        n_channels * nd)``; the shapes are read off the arrays as given
        (host or device -- nothing is converted here)."""
        nd = norb * norb
        B = int(view.n_channels)
        ND = B * nd
        for label, V in (("S", S), ("C", C)):
            shape = tuple(np.shape(V))
            if shape[-2:] != (ND, ND):
                raise ValueError("BondDeviceContext.for_view: the vertex {} has trailing "
                                 "shape {}, not the (ND, ND) = ({}, {}) of view.n_channels "
                                 "* norb**2 = {} * {}".format(label, shape[-2:], ND, ND, B, nd))
        if perm is None:
            perm = _mixed_pair_permutation(B, nd, norb)
        return cls(xp, S, C, S_on, C_on, perm, mixed_block_mask(B, nd),
                   green0_tail=green0_tail)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.release()
        return False

    def __getattr__(self, name):
        # attribute access for the owned arrays; everything else is normal
        if name in BondDeviceContext._NAMES:
            arrays = self.__dict__.get("_arrays")
            if self.__dict__.get("_released", True) or arrays is None:
                raise RuntimeError("BondDeviceContext: released")
            return arrays[name]
        raise AttributeError(name)

    def release(self):
        self._arrays = {}
        self._released = True

    @property
    def released(self):
        return self._released


def assemble_bubble(store, green_scf, green0_tail, beta, view, spatial_shape, workers):
    """Fill ``store['chibar']`` pair by pair from ``bubble._iter_bond_dynamic``
    (spec 3.1). ``green_scf`` is the TAIL-SUBTRACTED single-block Green
    function the general bubble consumes; ``green0_tail`` its tail.

    The bubble computes on the ARRAY MODULE of ``green_scf`` (issue
    #196): handed a device-resident Green function it runs the whole
    per-pair pipeline on the device and only the finished block crosses
    to the host store below. ``green0_tail`` must therefore live on the
    same module as ``green_scf`` -- a mismatch is refused here, by name,
    rather than as the backend-equality message of
    ``bubble._validate_dense_inputs``. Nothing is transferred here: the
    caller owns the placement of both arrays."""
    from .hartree_fock import NonFiniteError
    if green0_tail is not None:
        xp_g = _bk.array_module_of(green_scf)
        xp_t = _bk.array_module_of(green0_tail)
        if xp_t is not xp_g:
            raise ValueError(
                "assemble_bubble: green0_tail lives on {} but green_scf on {}; the bond "
                "bubble computes on the array module of green_scf, so its paired tail has "
                "to be on that same module (transfer it once per solve, e.g. through "
                "BondDeviceContext(green0_tail=...))".format(
                    getattr(xp_t, "__name__", xp_t), getattr(xp_g, "__name__", xp_g)))
    for (alpha, beta_), block in _bubble._iter_bond_dynamic(
            green_scf, green0_tail, beta, view, spatial_shape=tuple(spatial_shape),
            workers=workers):
        block = _bk.to_host(block)
        if not np.all(np.isfinite(block)):
            raise NonFiniteError("non-finite bond bubble block ({}, {})".format(alpha, beta_))
        store.put_pair("chibar", alpha, beta_, block)
        del block


# =============================================================================
# Batched dressing, effective interaction and collapses (spec 3.2-3.3)
# =============================================================================

from dataclasses import dataclass as _dataclass

from . import bond_channels as _bc


from .hartree_fock import NonFiniteError as _NonFiniteError


def _at(iteration):
    return "" if iteration is None else " (SCF iteration {})".format(iteration)


def _dress(cb, V, channel, l0, nmat, spatial_shape, cond_tol, iteration, guard_freqs="all",
           guard_policy="refuse", violations=None, guard_method="auto", stats=None):
    try:
        chi_b, cond = _bc.dress_batch(cb, V, channel, l0=l0, nmat=nmat, spatial_shape=spatial_shape,
                                      cond_tol=cond_tol, guard_freqs=guard_freqs,
                                      guard_policy=guard_policy, violations=violations,
                                      # the warnings of a tolerated violation do not
                                      # pass through the refusal wrapper below, so the
                                      # iteration reaches them at the source
                                      iteration=iteration, guard_method=guard_method, stats=stats)
    except ValueError as exc:
        if iteration is None:
            raise
        message = "{}{}".format(exc, _at(iteration))
        if isinstance(exc, _bc.BondConditioningError):
            # keep the structured refusal (issue #198) on the exception the
            # production path exposes, with the iteration appended
            raise _bc.BondConditioningError(
                message, channel=exc.channel, iq=exc.iq, q=exc.q, l=exc.l, worst=exc.worst,
                ratio=exc.ratio, pole=exc.pole, smin=exc.smin, smax=exc.smax,
                cond_tol=exc.cond_tol) from exc
        raise ValueError(message) from exc
    xp = _bk.array_module_of(chi_b)
    if not bool(xp.all(xp.isfinite(chi_b))):
        raise _NonFiniteError("non-finite dressed {} channel in the frequency batch starting "
                              "at l={}{}".format(channel, l0, _at(iteration)))
    return chi_b, cond


def _mixed_pair_permutation(B, nd, norb):
    """Index vector on the ``(B nd)`` pair axis: the IDENTITY on channel 0,
    the orbital-pair transpose ``(l1, l2) -> (l2, l1)`` inside every
    ``m != 0`` block.

    Applied to BOTH pair axes of the masked second-order product, it reads
    the bond-side leg of the bubble and the bond vertex at the transposed
    orbital pair -- the exact mixed on-site x bond second order (spec
    2026-09-16 R3, issue #192). It is the identity at ``norb = 1`` and on
    orbital-diagonal bonds, so no single-orbital or orbital-diagonal result
    moves under it.

    ``B = ND // nd`` is the channel count and ``nd = norb * norb`` the pair
    dimension of one channel block; the pair index inside block ``m`` is
    ``m * nd + a * norb + b`` (the flattening
    :func:`~hwave.solver.bond_channels.build_sc_bond_channel` places the
    bond coefficients in)."""
    perm = np.arange(B * nd)
    for m in range(1, B):
        for l1 in range(norb):
            for l2 in range(norb):
                perm[m * nd + l1 * norb + l2] = m * nd + l2 * norb + l1
    return perm


def mixed_block_mask(B, nd):
    """Block-weight mask of the second-order MIXED term on the ``(B nd)``
    pair axis (issue #198: one home for what the caller and the tests used
    to spell out inline): ``0.5`` on the channel-0 row and column blocks
    (the mixed on-site x bond blocks, taken once as the exchange
    skeleton), ``0`` on the channel-0 / channel-0 block (the channel-0
    second order is added separately) and ``0`` on every bond / bond
    block (their ring second order is the direct skeleton again). Applied
    as ``A *= 0.5 * mask`` to ``S chibar S + C chibar C`` in
    :func:`dress_and_build_w`. ``float64`` ``(B nd, B nd)``."""
    ND = B * nd
    mask = np.zeros((ND, ND))
    mask[:nd, :] = 0.5
    mask[:, :nd] = 0.5
    mask[:nd, :nd] = 0.0
    return mask


@_dataclass(frozen=True)
class DressResult:
    collapse0: np.ndarray      # (nmat, nvol, nd, nd)  channel-0 block of chibar
    collapse_s: np.ndarray     # (nmat, nvol, nd, nd)  channel-0 block of chi_s
    collapse_c: np.ndarray
    static_s: np.ndarray       # (nvol, ND, ND) at Omega = 0
    static_c: np.ndarray
    cond_min_s: float
    cond_min_c: float
    #: violations tolerated under guard_policy = "warn" in THIS map, both
    #: channels and both guard kinds (GitHub issue #199); always 0 under
    #: the default "refuse" policy, which raises instead
    guard_violations: int = 0
    #: how many of those were findings of the reduced "static" mode's solve
    #: RESIDUAL rather than of the conditioning guard (the conditioning count
    #: is the difference), so a diagnostic can name the guard that spoke
    guard_residual_violations: int = 0
    #: blocks the conditioning guard decomposed exactly in this map, and the
    #: number it guarded (issue #197; equal under guard_method = "svd")
    guard_exact_blocks: int = 0
    guard_blocks: int = 0


def dress_and_build_w(store, dev, *, nb, output_full, nmat, nvol, nd, spatial_shape,
                      cond_tol=_bc._BOND_COND_FLOOR, iteration=None, factors=None,
                      second_order="takimoto", guard_freqs="all", guard_policy="refuse",
                      guard_method="auto"):
    """The spec 3.2-3.3 loop (rev 19): per frequency batch, dress spin then
    charge (one channel batch alive at a time) and consume each into the
    effective interaction

        W = 3/2 S (chi_s - chibar) S + 1/2 C (chi_c - chibar) C + W2

    whose second-order part ``W2`` is exact for the off-site content and
    reduces to the general path's for the on-site content:
    channel-0 block ``[3/2 S chibar S + 1/2 C chibar C]_00 - 1/4 (S_on +
    C_on) chibar_00 (S_on + C_on)`` (the general path's subtraction
    restricted to the ON-SITE vertices ``S_on``, ``C_on``: a spin-
    independent density interaction is counted once by the ring),
    mixed channel-0/bond blocks ``1/4 (S chibar S + C chibar C)`` read at
    the PERMUTED pair index on both axes (the exchange skeleton, once:
    :func:`_mixed_pair_permutation` -- identity on channel 0, the
    orbital-pair transpose inside every bond block, so that the bond-side
    leg of the bubble and the bond vertex are taken at the transposed
    orbital pair; identity at ``norb = 1`` and on orbital-diagonal bonds)
    and bond-bond blocks ``0`` (their ring second order is the direct
    skeleton again). Also collects the channel-0 collapses, the static
    slices (slice assignment into preallocated buffers) and, with
    ``output_full``, the store's ``chi_s_w``/``chi_c_w``.

    ``second_order`` selects the CHANNEL-0 second order only (spec
    2026-09-08 D5), exactly as ``flex_second_order`` does for the
    standalone general path:

    ``"takimoto"``
        the rev-19 expression quoted above -- unchanged, byte for byte.
    ``"local"``
        the exact local kernel of :mod:`hwave.solver.second_order`
        instead: the channel-0 block is the ring restricted to channel 0
        (already in ``W_b``, the bare bubble subtracted from both
        channels) plus ``W2(chibar_00)`` from
        :func:`~hwave.solver.second_order.accumulate_batch`, which needs
        the compiled ``factors`` pack (host-backed: this whole module is
        host-side). Since the general path's channel-0 flattening is the
        bond store's, the two agree wherever the off-site content is
        declared zero -- that is gate G0 under BOTH kernels.

    The mixed and bond-bond blocks are the same under both values: the
    gate resums the off/off exchange topology the local kernel does not
    carry, and that is the D4/D5 design.

    ``dev`` is the :class:`BondDeviceContext` owning the vertices (spec
    2026-09-17 section 4.2): each frequency batch is moved from the store
    (host) onto ``dev.xp`` once (``_bk.to_device``), dressed and reduced to
    ``W`` entirely on that module, then the collapses, the static slices
    and ``W`` itself are copied back to the host (``_bk.to_host``) before
    the ONE store write of the batch. On numpy (``dev.xp is np``) every
    transfer is the identity, so this is byte-for-byte the host loop.

    ``cond_tol`` is the conditioning floor handed to every
    :func:`~hwave.solver.bond_channels.dress_batch` call and
    ``guard_policy`` its policy (GitHub issue #199): under ``"warn"`` a
    violation of the conditioning guard (or, with
    ``guard_freqs = "static"``, of the solve residual) is logged and the
    map continues, and the number of such violations over BOTH channels
    and the whole frequency grid is returned as
    ``DressResult.guard_violations``, of which
    ``DressResult.guard_residual_violations`` were residual findings.

    ``guard_method`` (GitHub issue #197) is threaded to every
    :func:`~hwave.solver.bond_channels.dress_batch` call unchanged
    (``"auto"``, ``"svd"`` or ``"interval"``); the number of blocks the
    guard decomposed exactly and the number it guarded, summed over both
    channels and every frequency batch, are returned as
    ``DressResult.guard_exact_blocks`` and ``DressResult.guard_blocks``."""
    if second_order not in ("local", "takimoto"):
        raise ValueError("dress_and_build_w: second_order must be \"local\" or \"takimoto\", "
                         "got {!r}".format(second_order))
    if second_order == "local" and factors is None:
        raise ValueError("dress_and_build_w: flex_second_order = \"local\" needs the factors")
    xp = dev.xp
    S, C, SpC_on = dev.S, dev.C, dev.SpC_on
    ND = S.shape[-1]
    norb = int(round(nd ** 0.5))
    if norb * norb != nd:
        raise ValueError("dress_and_build_w: nd = {} is not a square of an orbital "
                         "count".format(nd))
    # the pair axis is a whole number of channel blocks, or the block layout
    # the permutation below assumes (pair index m * nd + a * norb + b) does
    # not describe this matrix and B = ND // nd would silently truncate.
    # BondBlockStore makes the same check on ITS ND; this one is on the
    # VERTEX's, which is a separate argument and need not be the store's
    if ND % nd != 0:
        raise ValueError("dress_and_build_w: the vertex pair dimension ND = {} is not a "
                         "multiple of the channel block size nd = {}, so it does not carry "
                         "whole bond-channel blocks".format(ND, nd))
    # the pair permutation and the block-weight mask, owned by the device
    # context (built once, ahead of every SCF iteration)
    perm, mask = dev.perm, dev.mask
    collapse0 = np.empty((nmat, nvol, nd, nd), dtype=np.complex128)      # HOST accumulators
    collapse_s = np.empty_like(collapse0)
    collapse_c = np.empty_like(collapse0)
    static_s = np.zeros((nvol, ND, ND), dtype=np.complex128)
    static_c = np.zeros((nvol, ND, ND), dtype=np.complex128)
    l_static = nmat // 2
    cond_s = cond_c = np.inf
    # one list for the whole map: every warned guard violation of either
    # channel and either guard kind (empty under guard_policy = "refuse")
    violations = []
    # accumulated over both channels and every frequency batch (issue #197)
    stats = {"guard_exact_blocks": 0, "guard_blocks": 0}
    for l0 in range(0, nmat, nb):
        l1 = min(nmat, l0 + nb)
        cb = _bk.to_device(store.get_freq_batch("chibar", l0, l1), xp)     # 1 H2D
        collapse0[l0:l1] = _bk.to_host(cb[:, :, :nd, :nd])
        chi_s_b, cs = _dress(cb, S, "spin", l0, nmat, spatial_shape, cond_tol, iteration,
                             guard_freqs, guard_policy, violations, guard_method, stats)
        cond_s = min(cond_s, cs if cs is not None else np.inf)
        collapse_s[l0:l1] = _bk.to_host(chi_s_b[:, :, :nd, :nd])
        if l0 <= l_static < l1:
            static_s[...] = _bk.to_host(chi_s_b[l_static - l0])
        if output_full:
            store.put_freq_batch("chi_s_w", l0, l1, _bk.to_host(chi_s_b))
        chi_s_b -= cb
        W_b = 1.5 * (S[None] @ chi_s_b @ S[None])
        del chi_s_b
        chi_c_b, cc = _dress(cb, C, "charge", l0, nmat, spatial_shape, cond_tol, iteration,
                             guard_freqs, guard_policy, violations, guard_method, stats)
        cond_c = min(cond_c, cc if cc is not None else np.inf)
        collapse_c[l0:l1] = _bk.to_host(chi_c_b[:, :, :nd, :nd])
        if l0 <= l_static < l1:
            static_c[...] = _bk.to_host(chi_c_b[l_static - l0])
        if output_full:
            store.put_freq_batch("chi_c_w", l0, l1, _bk.to_host(chi_c_b))
        chi_c_b -= cb
        W_b += 0.5 * (C[None] @ chi_c_b @ C[None])
        del chi_c_b
        # second-order term. A / Bc carry BOTH the channel-0 block of the
        # legacy kernel and the mixed blocks of either kernel, so they are
        # built once, before the channel-0 branch.
        A = S[None] @ cb @ S[None]
        Bc = C[None] @ cb @ C[None]
        if second_order == "local":
            # accumulate_batch is array-module generic; `factors` must carry arrays
            # of xp's module (the caller passes the device pack on the GPU)
            _so.accumulate_batch(W_b[:, :, :nd, :nd], cb[:, :, :nd, :nd], l0, factors)
        else:
            W_b[:, :, :nd, :nd] += 1.5 * A[:, :, :nd, :nd] + 0.5 * Bc[:, :, :nd, :nd]
            W_b[:, :, :nd, :nd] -= 0.25 * (SpC_on[None] @ cb[:, :, :nd, :nd] @ SpC_on[None])
        A += Bc
        del Bc
        A *= 0.5 * mask
        if norb > 1:
            # at norb = 1 the pair permutation is the identity, and the
            # fancy-index copy would allocate a whole batch-shaped temporary
            # to reproduce A
            A = A[:, :, perm[:, None], perm[None, :]]
        W_b += A
        del A
        if not bool(xp.all(xp.isfinite(W_b))):
            raise _NonFiniteError("non-finite effective interaction W in the frequency batch "
                                  "[{}, {}){}".format(l0, l1, _at(iteration)))
        store.put_freq_batch("W", l0, l1, _bk.to_host(W_b))          # written ONCE, last
        del W_b, cb
    return DressResult(collapse0=collapse0, collapse_s=collapse_s, collapse_c=collapse_c,
                       static_s=static_s, static_c=static_c,
                       cond_min_s=float(cond_s), cond_min_c=float(cond_c),
                       guard_violations=len(violations),
                       guard_residual_violations=sum(1 for v in violations
                                                     if v["kind"] == "residual"),
                       guard_exact_blocks=int(stats["guard_exact_blocks"]),
                       guard_blocks=int(stats["guard_blocks"]))


# =============================================================================
# Bond-aware self-energy transport (spec 3.4)
# =============================================================================



def calc_self_energy_bond(store, green_kw, beta, view, shape, norb, workers, xp=np):
    """Sigma_fluct(k, iw) from the bond-resolved effective interaction ``W``
    in ``store`` (spec 3.4 rev 19, normative equation):

        Sigma_ab(k) = T/N sum_{q,nu} sum_{alpha,beta,c,d}
            e^{+i(k-q).(R_beta - R_alpha)} W_{(alpha,c,a),(beta,d,b)}(q) G_cd(k-q)

    The bond form factor sits on the internal leg k - q at both ends of W,
    with the block (alpha, beta) of W paired with the bubble block
    (beta, alpha)'s phase: the bubble is chibar_{alpha beta}(q) = -(T/N)
    sum_k e^{i(k-q).(R_alpha - R_beta)} G(k) G(k-q) and both functional
    derivatives of the ring functional with respect to G carry the phase
    of the differentiated (transposed) block (pinned by the second-order
    skeletons: with this transport, AT norb = 1, the channel-0, mixed and
    bond-bond blocks of 1/2 (S chibar S + C chibar C) are the direct
    skeleton, twice the exchange skeleton and half the direct skeleton to
    1e-15.  At norb > 1 the "twice the exchange skeleton" reading of the
    mixed blocks holds only once their pair index is permuted -- see
    :func:`_mixed_pair_permutation`, which dress_and_build_w applies and
    which is the identity at norb = 1; with it the mixed blocks are the
    exact mixed on-site x bond second order on every bond).  With
    chibar(q, i nu)^dagger = chibar(q, -i nu) the symmetry
    Sigma(k, i w)^dagger = Sigma(k, -i w) holds to round-off.  In real
    space the phase is G(r + R_beta - R_alpha), a roll of G by R_alpha - R_beta.

    ``green_kw`` is rank five `(1, nmat, nvol, norb, norb)`; the result has
    the same shape.  With a single on-site channel this reproduces
    `FLEX._calc_self_energy_general` byte for byte.

    ``xp`` selects the array module the accumulation runs in (numpy or
    cupy; default numpy). ``green_kw`` may be host or device -- it is moved
    to ``xp`` here via :func:`_bk.to_device`, as is each host ``W`` block
    read from ``store``. The returned ``sigma`` is an ``xp`` array; the
    caller moves it to the host. With ``xp is np`` this is the previous
    code, line for line."""
    if not (xp is np or (hasattr(xp, "zeros") and _bk.array_module_of(xp.zeros(1)) is xp)):
        raise TypeError("calc_self_energy_bond: xp must be numpy or cupy")
    nx, ny, nz = (int(x) for x in shape)
    nvol = nx * ny * nz
    nmat = green_kw.shape[1]
    P = norb
    nd = norb * norb
    B = view.n_channels
    G_kw = _bk.to_device(green_kw, xp)[0]
    G_rt = _bk.spatial_ifftn(
        _ms.fermion_to_tau(G_kw.reshape(nmat, nvol * P * P), axis=0).reshape(nmat, nx, ny, nz, P * P),
        axes=(1, 2, 3), workers=workers).reshape(nmat, nx, ny, nz, P, P)
    Sigma_rt = xp.zeros((nmat, nvol, P, P), dtype=xp.complex128)
    axes = (1, 2, 3)
    for alpha in range(B):
        Ra = np.asarray(view.delta_r[alpha], dtype=int)
        for bt in range(B):
            Rb = np.asarray(view.delta_r[bt], dtype=int)
            shift = tuple(int(x) for x in (Ra - Rb))   # e^{+ik'.(R_b - R_a)} G(k') = G(r + R_b - R_a)
            if shift == (0, 0, 0):
                G_sh = G_rt.reshape(nmat, nvol, P, P)
            else:
                G_sh = xp.roll(G_rt, shift, axis=axes).reshape(nmat, nvol, P, P)
            blk = _bk.to_device(np.ascontiguousarray(store.get_pair("W", alpha, bt)), xp)   # 1 H2D per pair
            blk_qt = _ms.boson_to_tau(blk.reshape(nmat, nvol * nd * nd), axis=0)
            del blk
            Wab_rt = _bk.spatial_ifftn(blk_qt.reshape(nmat, nx, ny, nz, nd * nd),
                                       axes=axes, workers=workers).reshape(nmat, nvol, P, P, P, P)
            del blk_qt
            A = xp.einsum('frcadb,frcd->frab', Wab_rt, G_sh)
            del Wab_rt, G_sh
            Sigma_rt += A
            del A
    del G_rt
    tmp = _bk.spatial_fftn(Sigma_rt.reshape(nmat, nx, ny, nz, P * P), axes=axes, workers=workers)
    del Sigma_rt
    sigma = _ms.tau_to_fermion(tmp.reshape(nmat, nvol * P * P), axis=0)
    del tmp
    sigma = sigma.reshape(1, nmat, nvol, P, P)
    sigma *= 1.0 / beta
    return sigma


# =============================================================================
# Memory preflight (spec 3.6) and output-path resolution (spec 4.2)
# =============================================================================

import os as _os

_DRESSING_OPS_WARN = 1.0e12
_TRANSPORT_OPS_WARN = 1.0e11
_GIB = float(1024 ** 3)


#: the device phase rows, in the order a refusal names them when two of
#: them are exactly equal
_DEV_PHASE_ORDER = ("bubble", "dressing", "transport")


def _largest_dev_phase(rows):
    """Name of the largest device phase row in ``rows``; an exact tie goes
    to the earlier name of :data:`_DEV_PHASE_ORDER`."""
    return max(_DEV_PHASE_ORDER,
               key=lambda k: (rows[k], -_DEV_PHASE_ORDER.index(k)))


def dressing_ops(nmat, nvol, ND):
    """Operation count of the two dense ND x ND solves over all (l, q)."""
    return 2.0 * nmat * nvol * float(ND) ** 3


def transport_ops(B, nmat, nvol, norb):
    """Operation count of the bond self-energy transport (transforms plus
    the orbital contraction, P^3 = norb^6)."""
    P = norb ** 2
    return float(B * B) * nmat * nvol * (P ** 2 * np.log2(max(nmat * nvol, 2)) + P ** 3)


def estimate_bond_memory(*, nmat, nvol, norb, B, depth, output_full, split_seed, n_types,
                         freq_batch, cap_gb, mixing, factor_bytes=0, device_available=None):
    """The named-buffer lifetime table of spec 3.6 (every row raw), the
    batch selection and the admission decision against ``cap_gb`` (binary
    GiB). Returns a dict with ``persistent_rows``, ``phase_rows`` (at the
    selected ``nb``), ``persistent``, ``nb``, ``peak`` and the symbols;
    raises ``ValueError`` naming every row when even ``nb = 1`` exceeds
    the cap, or when ``freq_batch`` does.

    ``factor_bytes`` is the compiled second-order factor pack's
    ``SecondOrderFactors.nbytes`` (0 with the legacy
    ``flex_second_order = takimoto``, which compiles none): the pack is
    built once per solver and lives for the whole solve, so it is a PERSISTENT row
    (``second_order_factors``). Its batch temporaries are not a row of
    their own -- they are two ``(nb, nvol, nd, nd)`` arrays, i.e. 2/B^2
    of a single ``(nb, nvol, ND, ND)`` buffer, well inside the
    ``dressing`` row's six.

    ``device_available`` (bytes free on the GPU, measured by the caller
    at the entry of the Phase B solve, before the block store and the
    device vertex context are allocated) is optional; when given, a
    second, incremental DEVICE table is built alongside the host one and
    the returned ``nb`` is capped to what both agree on. The device
    table carries what the gated run keeps resident on the device for
    the whole Phase B solve --

    * ``vertices_static = 5 * S``: the static vertex stack of the device
      context, allocated right after the measurement;
    * ``flex_arrays = 5 * G``: the SCF loop's device Green functions and
      self-energies, allocated per iteration INSIDE the loop and so not
      yet part of the measured reading;
    * ``green0_tail = G``: the bubble's tau-space tail, transferred to the
      device by the vertex context right after the measurement (issue
      #196). It is counted unconditionally: the tail exists whenever
      ``coeff_tail != 0``, and counting it when there is none is a
      deliberate margin, as with ``second_order_factors``;
    * ``second_order_factors = factor_bytes``: the device mirror of the
      compiled factor pack, which exists exactly when the pack does
      (``flex_second_order = "local"``). The mirror is already live at
      the measurement point, so counting it is a deliberate margin
      rather than a missing allocation --

    plus three phase rows that are never simultaneously live:
    ``bubble = max(prep, pair)`` -- the CPU expression, taken over as a
    DEVICE row since issue #196 put the bond bubble on the solver's
    array module. It is an UPPER BOUND there: the streamed device kernel
    never holds the static ``S``-sized buffers the host row was derived
    with, so the real device allocation is smaller --
    ``dressing(nb) = 7 * nb * nvol * ND^2 * 16`` during the per-batch
    dressing solve and ``transport = 6 * C`` during the bond
    self-energy transport. ``guard(nb) = 3 * nb * nvol * ND^2 * 16`` is the
    conditioning guard's own device temporaries (issue #197):
    INCREMENTAL to the dressing row rather than a phase of its own (the
    guard runs inside the dressing phase, never simultaneously with the
    bubble or the transport), so it is added to ``dressing`` in the need
    below but never competes for the name a refusal gives the largest
    phase. The device need at a batch size is therefore
    ``1.25 * (vertices_static + flex_arrays + green0_tail +
    second_order_factors + max(bubble, dressing(nb) + guard(nb), transport))``
    against ``device_cap = 0.9 *
    device_available``. Refusal at ``nb = 1``
    names whichever of the three phase rows (``bubble``, ``dressing``,
    ``transport`` -- the RAW rows, not ``dressing + guard``) is the
    largest (an exact tie goes to the earlier of ``bubble``, ``dressing``,
    ``transport``); an explicit
    ``freq_batch`` is checked against both the host and the device
    table; otherwise the selected ``nb`` is ``min`` of the largest
    batch each table admits, and the dict gains ``device_rows``,
    ``device_need``, ``device_cap``, ``device_nb`` and ``device_table``.

    Giving ``device_available`` also changes ONE host row: the bubble no
    longer allocates its temporaries on the host, which then holds only
    the finished channel-pair block on its way into the store, so the
    host ``bubble`` row becomes ``C``. The first, host-only preflight
    call runs before the backend is known and keeps the CPU expression,
    which is the conservative reading of the two."""
    nmat, nvol, norb, B = int(nmat), int(nvol), int(norb), int(B)
    depth = max(1, int(depth)) if mixing == "anderson" else 0      # the mixer's effective depth
    it = 16
    P = norb * norb
    ND = B * P
    U = nmat * nvol * ND * ND * it
    G = nmat * nvol * P * it
    C = nmat * nvol * P * P * it
    S = nvol * ND * ND * it
    H = nvol * P * it
    persistent_rows = {
        "chibar_W": 2 * U,
        "chi_sc_w": 2 * U if output_full else 0,
        "vertices_static": 5 * S,
        "state": G + H,
        "seed_envelope": (2 * G + H) if split_seed else 0,
        "anderson_history": 2 * depth * 2 * G,
        "anderson_work": 2 * max(depth - 1, 0) * 2 * G,       # the mixer's retained dR/dX stacks
        "collapses": 3 * C,
        "eigenpairs_hf": 5 * H,
        "flex_arrays": 5 * G,
        "hf_tables": int(n_types) * H,
        "second_order_factors": int(factor_bytes),
    }
    persistent = sum(persistent_rows.values())
    prep = it * nvol * (ND ** 2 + 4 * P * (nmat + 2))
    pair = it * nvol * (2 * ND ** 2 + 3 * nmat * P ** 2 + 2 * nmat * P + 8 * P)
    per_batch = 6 * nvol * ND * ND * it
    # On the GPU path the bubble's temporaries are the DEVICE's (issue
    # #196); the host only holds the one finished channel-pair block that
    # is on its way into the store, i.e. C bytes.
    host_bubble = C if device_available is not None else max(prep, pair)

    def _phase_rows(nb):
        return {
            "green_mu": 6 * G + H,
            # the evaluator's G_ref and G - G_ref (2 G), the Heff eigenpairs and
            # their workspace (4 H), rho_k/rho_r (2 H), rho_so/out (8 H), the
            # kernel's spin-major temporaries (20 H), the owning Sigma_HF copy (H)
            "density_hf": 2 * G + 4 * H + 2 * H + 8 * H + 20 * H + H,
            "bubble": host_bubble,
            "dressing": per_batch * nb,
            "transport": 4 * G + 4 * C,
            # materialised new total, green_inv + inverse + workspace, the G
            # difference, the two frequency-broadcast static stacks of the
            # component residual (7 G), the new static pair (H)
            "convergence": 7 * G + H,
            # anderson: residual + output (the history and work stacks are
            # persistent); linear: the difference, the scaled step, the new pair
            "mixing": 4 * G if mixing == "anderson" else 3 * G + H,
            "post_scf": G + (U if output_full else C),
        }

    def _peak(nb):
        return 1.25 * (persistent + max(_phase_rows(nb).values()))

    cap_bytes = float(cap_gb) * _GIB

    def _table(nb):
        rows = ["  persistent {:>18s}: {:10.4f} GiB".format(k, v / _GIB)
                for k, v in persistent_rows.items()]
        rows += ["  phase      {:>18s}: {:10.4f} GiB".format(k, v / _GIB)
                 for k, v in _phase_rows(nb).items()]
        return "\n".join(rows)

    if freq_batch is not None:
        nb = int(freq_batch)
        if not 1 <= nb <= nmat:
            raise ValueError("longitudinal_bond_freq_batch must be in [1, Nmat={}], got {}"
                             .format(nmat, nb))
        if _peak(nb) > cap_bytes:
            raise ValueError(
                "[mode.param] longitudinal_bond_freq_batch = {} gives an estimated peak of "
                "{:.4f} GiB = 1.25 * (persistent + max phase row) above "
                "longitudinal_bond_memory_cap_gb = {:.4f} GiB (B = {}, ND = {}, nvol = {}, "
                "Nmat = {}); the rows are\n{}\nReduce the batch, the k mesh or Nmat, drop "
                "declared-zero outer shells with longitudinal_bond_max_shells, or raise the cap."
                .format(nb, _peak(nb) / _GIB, cap_bytes / _GIB, B, ND, nvol, nmat, _table(nb)))
    else:
        if _peak(1) > cap_bytes:
            raise ValueError(
                "[mode.param] longitudinal_bond_channels=true: the estimated peak host memory "
                "{:.4f} GiB (frequency batch 1) = 1.25 * (persistent + max phase row) exceeds "
                "longitudinal_bond_memory_cap_gb = {:.4f} GiB (B = {}, ND = {}, nvol = {}, "
                "Nmat = {}); the rows are\n{}\nReduce the k mesh or Nmat, drop declared-zero "
                "outer shells with longitudinal_bond_max_shells, disable "
                "longitudinal_bond_output_full, lower the Anderson depth, or raise the cap."
                .format(_peak(1) / _GIB, cap_bytes / _GIB, B, ND, nvol, nmat, _table(1)))
        nb = 1
        for cand in range(nmat, 0, -1):
            if _peak(cand) <= cap_bytes:
                nb = cand
                break

    device = {}
    if device_available is not None:
        vertices = 5 * S
        transport = 6 * C
        dev_persistent = {"vertices_static": vertices, "flex_arrays": 5 * G,
                          "green0_tail": G,
                          "second_order_factors": int(factor_bytes)}
        dev_persistent_sum = sum(dev_persistent.values())
        dev_bubble = max(prep, pair)
        def _dev_rows(n):
            rows = dict(dev_persistent)
            rows["bubble"] = dev_bubble
            rows["dressing"] = 7 * n * nvol * ND * ND * it
            # the conditioning guard's device temporaries (issue #197):
            # incremental to the dressing row, not a phase of its own --
            # see _largest_dev_phase, which never names it
            rows["guard"] = 3 * n * nvol * ND * ND * it
            rows["transport"] = transport
            return rows
        def _dev_need(n):
            r = _dev_rows(n)
            return 1.25 * (dev_persistent_sum
                           + max(r["bubble"], r["dressing"] + r["guard"], r["transport"]))
        dev_cap = 0.9 * float(device_available)
        def _dev_table(n):
            return "\n".join("  device     {:>18s}: {:10.4f} GiB".format(k, v / _GIB)
                             for k, v in _dev_rows(n).items())
        if _dev_need(1) > dev_cap:
            phase = _largest_dev_phase(_dev_rows(1))
            raise ValueError(
                "[mode.param] gpu=true with longitudinal_bond_channels: the estimated device need "
                "{:.4f} GiB (frequency batch 1, phase '{}') = 1.25 * (persistent rows + max phase "
                "row) exceeds 0.9 * the available device memory {:.4f} GiB; the rows are\n{}\nReduce "
                "the k mesh or Nmat, drop declared-zero outer shells with "
                "longitudinal_bond_max_shells, or run with gpu=false.".format(
                    _dev_need(1) / _GIB, phase, dev_cap / _GIB, _dev_table(1)))
        if freq_batch is not None and _dev_need(int(freq_batch)) > dev_cap:
            raise ValueError(
                "[mode.param] longitudinal_bond_freq_batch = {} gives an estimated device need of "
                "{:.4f} GiB above 0.9 * the available device memory {:.4f} GiB (host peak {:.4f} GiB "
                "against the host cap {:.4f} GiB); the device rows are\n{}".format(
                    int(freq_batch), _dev_need(int(freq_batch)) / _GIB, dev_cap / _GIB,
                    _peak(int(freq_batch)) / _GIB, cap_bytes / _GIB, _dev_table(int(freq_batch))))
        device_nb = 1
        for cand in range(nmat, 0, -1):
            if _dev_need(cand) <= dev_cap:
                device_nb = cand
                break
        nb = min(nb, device_nb)
        device = dict(device_rows=_dev_rows(nb), device_need=_dev_need(nb), device_cap=dev_cap,
                      device_nb=device_nb, device_table=_dev_table(nb))

    phase_rows = _phase_rows(nb)
    return dict(persistent_rows=persistent_rows, phase_rows=phase_rows, persistent=persistent,
                nb=nb, peak=_peak(nb), cap_bytes=cap_bytes, U=U, G_bytes=G, C_bytes=C,
                S_bytes=S, H_bytes=H, B=B, ND=ND, nvol=nvol, nmat=nmat, table=_table(nb),
                dressing_ops=dressing_ops(nmat, nvol, ND),
                transport_ops=transport_ops(B, nmat, nvol, norb),
                **device)


_NPZ_ARTIFACTS = ("chi0q", "chiq_s", "chiq_c", "chiq", "sigma", "green", "longitudinal_bond",
                  "eliashberg_bond_singlet", "eliashberg_bond_triplet")
_DEFAULT_FILES = {"chiq_s": "chiq_s", "chiq_c": "chiq_c",
                  "longitudinal_bond": "longitudinal_bond.npz",
                  "eliashberg_bond_singlet": "eliashberg_bond_singlet",
                  "eliashberg_bond_triplet": "eliashberg_bond_triplet",
                  "eigenvalue_bond_singlet": "eigenvalue_bond_singlet.dat",
                  "eigenvalue_bond_triplet": "eigenvalue_bond_triplet.dat",
                  "gap_bond_singlet": "gap_bond_singlet.dat",
                  "gap_bond_triplet": "gap_bond_triplet.dat"}


def resolve_output_paths(info_outputfile, path_to_output, active_outputs):
    """``{name: absolute path}`` for the artifacts in ``active_outputs``
    that will actually be written (spec 4.2): joins with ``path_to_output``,
    applies the ``.npz`` suffix rule of ``numpy.savez`` to the archive
    artifacts only, normalises to absolute paths (symlink equivalence is
    not checked) and refuses two artifacts resolving to one file."""
    out = {}
    for name in active_outputs:
        fn = info_outputfile.get(name, _DEFAULT_FILES.get(name))
        if fn is None:
            continue
        p = _os.path.join(str(path_to_output), str(fn))
        if name in _NPZ_ARTIFACTS and not p.endswith(".npz"):
            p += ".npz"
        out[name] = _os.path.abspath(_os.path.normpath(p))
    seen = {}
    for name, p in out.items():
        if p in seen:
            raise ValueError(
                "[file.output] the artifacts '{}' and '{}' resolve to the same file {}; "
                "give them distinct file names (numpy appends '.npz' to archive names "
                "without that suffix).".format(seen[p], name, p))
        seen[p] = name
    return out
