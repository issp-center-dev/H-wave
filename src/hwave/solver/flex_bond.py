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
from . import threaded as _threaded
from . import half_grid as _hg
from . import staging as _staging


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

    ``aliases`` maps extra names onto allocated ones, so that two roles
    share one array: without in-process pairing the FLEX map writes ``W``
    over the ``chibar`` batch it has just consumed (issue #222), with
    ``aliases={"W": "chibar"}``. An alias cannot be detached or released
    on its own.

    ``half_grid`` (issue #223): store only the non-negative half of the
    uniform bosonic grid, ``nslices = nmat // 2 + 1`` slices per array, slice
    ``h`` holding the grid index ``l(h) = (h + nmat // 2) mod nmat``
    (:mod:`hwave.solver.half_grid`). ``nmat`` keeps its meaning (the
    solver's full grid); frequency-batch access is in SLICE index (``h``
    under ``half_grid``, ``l`` otherwise; :meth:`grid_index` names the
    ``l`` of each slice), and ``get_pair`` is a view on the stored slices.
    :meth:`get_pair_full` and :meth:`expand_full` rebuild the full grid
    from the symmetry ``X(q, l) = X(q, nmat - l)^dagger``. On the full grid
    every one of them is the identity it replaces.

    ``xp`` (issue #275): the array module every slot is allocated on --
    numpy (the default, the host store) or cupy (the DEVICE store, chosen by
    the device admission when the store fits). The store owns placement:
    ``put_pair`` / ``put_freq_batch`` convert their argument to ``xp`` once
    (:func:`backend.to_device`, which is the identity for an argument already
    on ``xp`` and the host copy ``to_host`` for a device argument of a host
    store), while ``get_pair`` / ``get_freq_batch`` return views of ``xp``
    and ``get_pair_full`` a contiguous array of the module asked for. On a
    device store the dressing's frequency batch is therefore a VIEW of the
    slot, exactly as on the CPU path, and the per-iteration host <-> device
    copies of the store disappear. ``detach`` returns an array of ``xp``.

    A HOST store on the GPU path (issue #275, part 2) moves its bytes
    through :mod:`hwave.solver.staging`: device blocks and batches arrive,
    and leave for the device (``get_pair_full`` / :meth:`get_freq_batch_on`),
    in chunks through one pinned staging buffer of at most
    ``staging_bytes`` (default :data:`staging.PINNED_STAGING_BYTES`), and the
    pair-strided host copies are split over ``copy_workers`` threads (default
    1, the plain single assignment; None picks
    :func:`staging.default_copy_workers`, as the GPU solve does). A
    DEVICE store leaves the device (:meth:`to_host`) through the same
    buffer. Only the route of the bytes changes, never a byte.
    """

    def __init__(self, nmat, nvol, ND, nd, names, aliases=None, *, half_grid=False, xp=np,
                 staging_bytes=None, copy_workers=1):
        names = tuple(names)
        if len(set(names)) != len(names):
            raise ValueError("BondBlockStore: duplicate array names {}".format(names))
        aliases = dict(aliases or {})
        for alias, target in aliases.items():
            if alias in names or target not in names:
                raise ValueError("BondBlockStore: alias {!r} -> {!r} must name a new name "
                                 "for an allocated array {}".format(alias, target, names))
        self._aliases = aliases
        if ND % nd != 0:
            raise ValueError("BondBlockStore: ND={} is not a multiple of nd={}".format(ND, nd))
        self.nmat, self.nvol, self.ND, self.nd = int(nmat), int(nvol), int(ND), int(nd)
        self.B = self.ND // self.nd
        self.half_grid = bool(half_grid)
        self.nslices = _hg.nslices(self.nmat) if self.half_grid else self.nmat
        self.xp = np if xp is None else xp
        self._arrays = {}
        self._released_slots = set()
        self._staging_bytes = staging_bytes
        self.copy_workers = copy_workers
        self._stager = None
        try:
            for name in names:
                self._arrays[name] = self.xp.zeros((self.nslices, self.nvol, self.ND, self.ND),
                                                   dtype=np.complex128)
        except (MemoryError,) + tuple(_bk.oom_error_types()):
            # a partial allocation is rolled back (a device store hands its
            # blocks back to the pool) before the error propagates
            self._arrays.clear()
            raise
        self._released = False

    @property
    def resident(self):
        """``"host"`` (numpy slots) or ``"device"`` (any other module)."""
        return "host" if self.xp is np else "device"

    @property
    def slot_names(self):
        """The names the store answers to now, as a tuple: the allocated
        slots (released or detached ones excluded) in allocation order, then
        the aliases whose target is still allocated."""
        return tuple(self._arrays) + tuple(a for a, t in self._aliases.items()
                                           if t in self._arrays)

    @property
    def nbytes(self):
        """Bytes of the slots currently allocated (aliases counted once)."""
        return int(sum(int(a.nbytes) for a in self._arrays.values()))

    def _to_store(self, x):
        # the one conversion point of the store (issue #275): the identity
        # for an argument already on the store's module
        return _bk.to_device(x, self.xp)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.release()
        return False

    def _arr(self, name):
        if self._released:
            raise RuntimeError("BondBlockStore: released")
        name = self._aliases.get(name, name)
        if name in self._released_slots:
            raise KeyError("BondBlockStore: slot {!r} was released".format(name))
        return self._arrays[name]

    def stager(self, xp_dev):
        """The store's pinned stager for the device module ``xp_dev``
        (created on first use, dropped by :meth:`release`)."""
        if self._stager is None or self._stager.xp is not xp_dev:
            if self._stager is not None:
                self._stager.release()
            self._stager = _staging.PinnedStager(xp_dev, self._staging_bytes,
                                                 self.copy_workers)
        return self._stager

    def _put(self, dst, x):
        # dst: a view of a slot; x: a block or batch of either module
        if self.xp is not np:
            dst[...] = self._to_store(x)
            return
        xp_x = _bk.array_module_of(x)
        if tuple(getattr(x, "shape", ())) != tuple(dst.shape) or \
                getattr(x, "dtype", None) != dst.dtype:
            dst[...] = _bk.to_host(x)      # broadcasting / casting: the plain path
        elif xp_x is np:
            _staging.parallel_copy(dst, x, self.copy_workers)
        else:
            # a device block into a host store: chunked through the pinned
            # buffer straight into the (strided) slot view
            self.stager(xp_x).copy_to_host(x, dst)

    def put_pair(self, name, alpha, beta, block):
        a = self._arr(name)
        nd = self.nd
        self._put(a[:, :, alpha * nd:(alpha + 1) * nd, beta * nd:(beta + 1) * nd], block)

    def get_pair(self, name, alpha, beta):
        a = self._arr(name)
        nd = self.nd
        return a[:, :, alpha * nd:(alpha + 1) * nd, beta * nd:(beta + 1) * nd]

    def put_pair_fullgrid(self, name, alpha, beta, block_full):
        """Store a FULL-grid ``(nmat, nvol, nd, nd)`` pair block: its slices
        ``l(h)`` under ``half_grid``, the whole block otherwise."""
        if self.half_grid:
            block_full = _hg.half_from_full(block_full, self.nmat)
        self.put_pair(name, alpha, beta, block_full)

    def get_pair_full(self, name, alpha, beta, xp=None):
        """The pair block ``(alpha, beta)`` on the FULL grid as a contiguous
        ``(nmat, nvol, nd, nd)`` array of module ``xp`` (default: the
        store's own). Under ``half_grid`` the slices ``l in (0, nmat/2)``
        are the conjugate transposes of the stored slices of the transposed
        pair ``(beta, alpha)``; on the full grid this is exactly
        ``to_device(ascontiguousarray(get_pair(...)), xp)``. The gather and
        the mirror run on the store's module (a device kernel on a device
        store); only a module change crosses host <-> device."""
        ab = self.get_pair(name, alpha, beta)
        cross = not (xp is None or xp is self.xp)
        if not self.half_grid:
            if self.xp is not np:
                full = self.xp.ascontiguousarray(ab)
            elif cross:
                # host store -> device: gathered chunk by chunk into the
                # pinned buffer, no host temporary of the whole block
                return self.stager(xp).to_device(ab)
            else:
                full = _staging.parallel_copy(np.empty(ab.shape, dtype=ab.dtype), ab,
                                              self.copy_workers)
        else:
            full = _hg.expand(ab, self.nmat, ba=self.get_pair(name, beta, alpha))
        if not cross:
            return full
        return self._cross(full, xp)

    def _cross(self, x, xp):
        """``x`` (of the store's module) on module ``xp`` through the stager."""
        if self.xp is np:
            return self.stager(xp).to_device(x)
        if xp is np:
            return self.stager(self.xp).to_host(x)
        return _bk.to_device(x, xp)

    def get_freq_batch_on(self, name, l0, l1, xp):
        """:meth:`get_freq_batch` on module ``xp``: the view itself when
        ``xp`` is the store's module (on a device store the dressing then
        works on the slot, as on the CPU path), otherwise a staged copy."""
        batch = self.get_freq_batch(name, l0, l1)
        if xp is None or xp is self.xp:
            return batch
        return self._cross(batch, xp)

    def grid_index(self, s0, s1):
        """The grid index ``l`` of the stored slices ``s0 ... s1 - 1``."""
        if self.half_grid:
            return _hg.grid_index(self.nmat, s0, s1)
        return np.arange(int(s0), int(s1), dtype=np.int64)

    def expand_full(self, name, out=None, chunk=None):
        """The whole array ``name`` on the full grid ``(nmat, nvol, ND, ND)``
        (``out`` may be a memory map); the stored array itself on the full
        grid when ``out`` is None. An ``out`` of another array module than
        the store's receives a converted copy."""
        a = self._arr(name)
        if out is not None and _bk.array_module_of(out) is not self.xp:
            a = _bk.to_device(a, _bk.array_module_of(out))
        if not self.half_grid:
            if out is None:
                return a
            out[...] = a
            return out
        return _hg.expand(a, self.nmat, out=out, chunk=chunk)

    def put_freq_batch(self, name, l0, l1, batch):
        self._put(self._arr(name)[l0:l1], batch)

    def get_freq_batch(self, name, l0, l1):
        return self._arr(name)[l0:l1]

    def to_host(self):
        """Move every remaining slot to host memory (issue #275): after the
        SCF the in-process pairing reads the store only to build its vertex
        and the kernel blocks, and the device memory is better spent on the
        pairing kernel. Slot by slot, so at most one slot is held twice; the
        store's module becomes numpy. The identity for a host store.

        Each slot comes over in chunks through the store's pinned staging
        buffer straight into its new host array (issue #275, part 2: one
        ``asnumpy`` of a 10 GiB slot ran far below the bus rate); without a
        pinned buffer it is the plain ``asnumpy``."""
        if self.xp is np:
            return
        st = self.stager(self.xp)
        for name in list(self._arrays):
            self._arrays[name] = st.to_host(self._arrays[name])
        self.xp = np

    def _refuse_alias(self, name, what):
        if name in self._aliases:
            raise ValueError("BondBlockStore: cannot {} {!r}, an alias of {!r}".format(
                what, name, self._aliases[name]))

    def detach(self, name):
        self._refuse_alias(name, "detach")
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
        self._refuse_alias(name, "release")
        if name not in self._arrays:
            if name in self._released_slots:
                return
            raise KeyError("BondBlockStore: no slot {!r}".format(name))
        del self._arrays[name]
        self._released_slots.add(name)

    def release(self):
        self._arrays.clear()
        self._released = True
        if self._stager is not None:
            self._stager.release()
            self._stager = None

    @property
    def released(self):
        return self._released


def freq_batch_on(source, name, l0, l1, xp):
    """``source.get_freq_batch(name, l0, l1)`` on module ``xp``: through
    :meth:`BondBlockStore.get_freq_batch_on` (staged, or the view itself on
    the store's module) when the source has it, else ``to_device``."""
    get_on = getattr(source, "get_freq_batch_on", None)
    if get_on is not None:
        return get_on(name, l0, l1, xp)
    return _bk.to_device(source.get_freq_batch(name, l0, l1), xp)


class BondDeviceContext:
    """Owner of the IMMUTABLE bond inputs for one ``_solve_phase_b`` call
    (spec 2026-09-17 section 4.2): the vertices ``S``, ``C``,
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

    _NAMES = ("S", "C", "perm", "mask", "green0_tail")

    def __init__(self, xp, S, C, perm=None, mask=None, green0_tail=None):
        self.xp = xp
        host = dict(S=S, C=C, perm=perm, mask=mask, green0_tail=green0_tail)
        self._arrays = {k: (None if host[k] is None else _bk.to_device(host[k], xp))
                        for k in self._NAMES}
        self._released = False

    @classmethod
    def for_view(cls, xp, S, C, view, norb, green0_tail=None, *, perm=None):
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
        return cls(xp, S, C, perm, mixed_block_mask(B, nd),
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
    half = bool(getattr(store, "half_grid", False))
    for (alpha, beta_), block in _bubble._iter_bond_dynamic(
            green_scf, green0_tail, beta, view, spatial_shape=tuple(spatial_shape),
            workers=workers):
        # the whole block is checked for finiteness on its own module, before
        # the half grid (issue #223) drops the mirrored slices, so a
        # non-finite slice the half grid would drop is still refused
        xp_b = _bk.array_module_of(block)
        if not bool(xp_b.all(xp_b.isfinite(block))):
            raise NonFiniteError("non-finite bond bubble block ({}, {})".format(alpha, beta_))
        if half:
            # only the stored slices l(h) are kept; on the device the
            # selection runs before any transfer
            block = _hg.half_from_full(block, store.nmat)
        # the store converts (issue #275): a host copy for a host store, the
        # block itself for a device store on the block's module
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


def _sandwich(left, x, right, workers):
    """``left[None] @ x @ right[None]``, thread-chunked on the host."""
    if _bk.array_module_of(x) is np:
        return _threaded.sandwich(left, x, right, workers=workers)
    return left[None] @ x @ right[None]


def _dress(cb, V, channel, l0, nmat, spatial_shape, cond_tol, iteration, guard_freqs="all",
           guard_policy="refuse", violations=None, guard_method="auto", stats=None, workers=1,
           guard_exact="always", ls=None, guard_only=False):
    try:
        chi_b, cond = _bc.dress_batch(cb, V, channel, l0=l0, nmat=nmat, spatial_shape=spatial_shape,
                                      ls=ls, guard_only=guard_only,
                                      cond_tol=cond_tol, guard_freqs=guard_freqs,
                                      guard_policy=guard_policy, violations=violations,
                                      # the warnings of a tolerated violation do not
                                      # pass through the refusal wrapper below, so the
                                      # iteration reaches them at the source
                                      iteration=iteration, guard_method=guard_method, stats=stats,
                                      workers=workers, guard_exact=guard_exact)
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
    if chi_b is None:              # guard_only: the diagnostics alone
        return chi_b, cond
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
    #: ``None`` when ``cond_min_s`` is the exact minimum of the map (always
    #: the case under ``guard_exact = "always"``). Otherwise (issue #276:
    #: the interval alone proved the pass of some frequency batch, so they
    #: were not decomposed) ``cond_min_s`` is NaN and this is the bracket
    #: ``(low, up, index)`` of the exact minimum, ``low <= min <= up``, with
    #: ``index = l * nvol + q`` the flattened grid index of the smallest
    #: upper end (``None`` when that end is an exactly decomposed batch's)
    cond_bracket_s: tuple = None
    cond_bracket_c: tuple = None
    #: half grid only (issue #223): whether the mirrored slices were guarded
    #: on their own in this map (``False`` when the time-reversal relation of
    #: :mod:`hwave.solver.half_grid` held, so that every mirrored score is a
    #: stored one), and the relation's measured relative deviation; ``None``
    #: on the full grid
    mirror_guarded: bool = None
    time_reversal_deviation: float = None


class HalfGridDressing:
    """The conditioning of the mirrored slices of a half-grid dressing pass
    (issue #223), shared by the FLEX map (:func:`dress_and_build_w`) and the
    in-process pairing's re-dressing
    (:meth:`~hwave.solver.eliashberg_bond.PairVertexAccumulator.add_dressed`),
    so that both apply ONE rule.

    The mirrored slices ``l' = nmat - l`` are skipped only where that is
    provably exact (:mod:`hwave.solver.half_grid`, relations (2)-(3)): the
    guard refuses on every violation (``guard_policy = "refuse"``) over
    every frequency (``guard_freqs = "all"``) and time reversal makes every
    mirrored score a stored one -- checked here on the stored ``chibar``,
    which the caller passes BEFORE any ``W`` batch overwrites it. Even then
    a batch whose minimum score lies within :func:`half_grid.decision_margin`
    of ``cond_tol`` has its mirror guarded on its own. Under ``"warn"`` and
    under the ``"static"`` residual the mirror is always guarded.

    The caller supplies ``dress(cb, V, channel, ls, guard_only=False)``,
    its own :func:`_dress` call for the rows ``ls`` of a batch (the stored
    rows, or with ``guard_only`` the diagnostics of their mirrors)."""

    def __init__(self, chibar_half, S, C, *, nmat, spatial_shape, cond_tol, guard_freqs,
                 guard_policy):
        self.nmat = int(nmat)
        self.cond_tol = cond_tol
        trs = _hg.time_reversal_check(chibar_half, S, C, spatial_shape)
        self.deviation = trs["deviation"]
        self.mirror_guard = not (trs["holds"] and guard_freqs == "all"
                                 and guard_policy == "refuse" and cond_tol is not None)
        self.margin = _hg.decision_margin(trs, int(S.shape[-1]), cond_tol)
        self.mq = _hg.minus_q_index(spatial_shape)
        self.spatial_shape = tuple(spatial_shape)

    def mirror(self, dress, cb, ls, V, channel, cmin):
        """The dressing diagnostics of the mirrored slices ``l' = nmat - l``
        of a batch (``l = 0`` and ``nmat / 2`` are their own mirrors and are
        skipped); their scores join ``cmin`` when it is given."""
        nmat = self.nmat
        sel = np.nonzero((ls != 0) & (ls != nmat // 2))[0]
        if not sel.size:
            return
        xp = _bk.array_module_of(cb)
        cb_m = xp.ascontiguousarray(xp.conj(cb[xp.asarray(sel)]).swapaxes(-1, -2))
        _chi, cm = dress(cb_m, V, channel, nmat - ls[sel], guard_only=True)
        del _chi, cb_m
        if cmin is not None:
            cmin.add(cm)

    def _near(self, score):
        # within the Weyl margin of cond_tol: the mirrored partner's score
        # could fall on the other side
        return score is not None and abs(score - self.cond_tol) <= self.margin * self.cond_tol

    def stored(self, dress, cb, ls, V, channel, cmin):
        """Dress the stored rows ``ls`` of a batch and return the dressed
        channel. Under the skip, the mirror is guarded on its own only for a
        batch whose minimum is within the margin of ``cond_tol`` (if the
        minimum clears it, every block does), and a refusal at a stored
        block is reported at its time-reversal mirror ``(-q, nmat - l)``,
        which has the same score and which the full grid meets first."""
        nmat = self.nmat
        if self.mirror_guard:
            chi, c = dress(cb, V, channel, ls)
            if cmin is not None:
                cmin.add(c)
            return chi
        try:
            chi, c = dress(cb, V, channel, ls)
        except _bc.BondConditioningError as exc:
            if exc.l is None or exc.l in (0, nmat // 2):
                raise
            if self._near(exc.worst):
                # the mirror's own guard decides: it refuses at the mirrored
                # location if that block fails too (the full grid meets nu < 0
                # first); otherwise the stored refusal stands
                self.mirror(dress, cb, ls, V, channel, None)
                raise
            # the decision is the full grid's; with several violating blocks
            # the location the full grid reports first can differ
            new = _bc.relocate_refusal(exc, nmat - exc.l, int(self.mq[exc.iq]), nmat,
                                       self.spatial_shape)
            if new is None:
                raise
            raise new from exc
        if cmin is not None:
            cmin.add(c)
        low = c.low if isinstance(c, _bc.GuardBracket) else c
        if self._near(low):
            self.mirror(dress, cb, ls, V, channel, cmin)
        return chi


def dress_and_build_w(store, dev, *, nb, output_full, nmat, nvol, nd, spatial_shape,
                      cond_tol=_bc._BOND_COND_FLOOR, iteration=None, factors=None,
                      guard_freqs="all", guard_policy="refuse",
                      guard_method="auto", workers=1, bond_w2=None, guard_exact="always"):
    """The spec 3.2-3.3 loop (rev 19): per frequency batch, dress spin then
    charge (one channel batch alive at a time) and consume each into the
    effective interaction

        W = 3/2 S (chi_s - chibar) S + 1/2 C (chi_c - chibar) C + W2

    whose second-order part ``W2`` is exact for the off-site content and
    reduces to the general path's for the on-site content (spec 2026-09-08
    D5): the channel-0 block is the ring restricted to channel 0 (already
    in ``W_b``, the bare bubble subtracted from both channels) plus the
    exact local kernel ``W2(chibar_00)`` of
    :func:`hwave.solver.second_order.accumulate_batch`, which needs the
    compiled ``factors`` pack (of the array module ``dev.xp``: the caller
    passes the device mirror on the GPU). Since the general path's
    channel-0 flattening is the bond store's, the two agree wherever the
    off-site content is declared zero -- that is gate G0. The
    mixed channel-0/bond blocks are ``1/4 (S chibar S + C chibar C)`` read at
    the PERMUTED pair index on both axes (the exchange skeleton, once:
    :func:`_mixed_pair_permutation` -- identity on channel 0, the
    orbital-pair transpose inside every bond block, so that the bond-side
    leg of the bubble and the bond vertex are taken at the transposed
    orbital pair; identity at ``norb = 1`` and on orbital-diagonal bonds)
    and bond-bond blocks ``0`` (their ring second order is the direct
    skeleton again). Also collects the channel-0 collapses, the static
    slices (slice assignment into preallocated buffers) and, with
    ``output_full``, the store's ``chi_s_w``/``chi_c_w``.

    The mixed and bond-bond blocks resum the off/off exchange topology the
    local kernel does not carry (the D4/D5 design).

    ``dev`` is the :class:`BondDeviceContext` owning the vertices (spec
    2026-09-17 section 4.2): each frequency batch is moved from the store
    onto ``dev.xp`` once (``_bk.to_device``; the identity for a device store,
    issue #275), dressed and reduced to ``W`` entirely on that module, then
    the collapses and the static slices are copied back to the host
    (``_bk.to_host``) and ``W`` is handed to the ONE store write of the
    batch, which converts it to the store's module. On numpy (``dev.xp is
    np``) every transfer is the identity, so this is byte-for-byte the host
    loop.

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
    ``DressResult.guard_exact_blocks`` and ``DressResult.guard_blocks``.

    ``guard_exact`` (issue #276) is threaded to every
    :func:`~hwave.solver.bond_channels.dress_batch` call unchanged
    (``"always"``, the default, or ``"decide"``). Under ``"decide"`` a
    frequency batch whose interval alone proves the pass is not decomposed,
    and when the map's minimum then is not known exactly,
    ``DressResult.cond_min_s`` / ``cond_min_c`` are NaN and
    ``cond_bracket_s`` / ``cond_bracket_c`` carry its bracket.

    ``bond_w2`` (issue #181): the factors ``[(triple, A, B), ...]`` of
    :func:`hwave.solver.second_order.build_bond_w2_factors`, on ``dev.xp``.
    When given, the second order of every block is ``1/2 sum A chibar B``
    (exact in the bond-enlarged basis, off-site Exchange included) and the
    channel-0 kernel, the mixed-block rule and the zero bond-bond blocks
    above are not used; ``factors`` is then ignored.

    The ring's pair orientation (issue #251): with or without ``bond_w2``,
    the ring, third order and up, is taken with the two orbitals of every
    pair swapped inside its channel block -- the conserving orientation
    under :func:`calc_self_energy_bond`'s contraction."""
    if bond_w2 is not None:
        factors = None
    elif factors is None:
        raise ValueError("dress_and_build_w: the second-order kernel needs the factors "
                         "(or bond_w2)")
    xp = dev.xp
    S, C = dev.S, dev.C
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
    # the ring's pair orientation (issue #251): the orbitals of every pair are
    # swapped inside its channel block, (m, c, a) -> (m, a, c) -- the
    # orientation in which the ring is the functional derivative of
    # Tr (chibar U)^n under calc_self_energy_bond's contraction (pinned by
    # tests/test_flex_ring_orientation.py). The identity at norb = 1.
    pin = None
    if norb > 1:
        pin = xp.asarray(np.arange(ND).reshape(ND // nd, norb, norb).transpose(0, 2, 1).ravel())
    collapse0 = np.empty((nmat, nvol, nd, nd), dtype=np.complex128)      # HOST accumulators
    collapse_s = np.empty_like(collapse0)
    collapse_c = np.empty_like(collapse0)
    static_s = np.zeros((nvol, ND, ND), dtype=np.complex128)
    static_c = np.zeros((nvol, ND, ND), dtype=np.complex128)
    # issue #223: under the store's half grid the batches run over the
    # stored slices h = 0 ... nmat/2 (the static slice is h = 0), the
    # dressing reports the true grid index l(h), and the collapses are
    # mirrored onto the full grid before return
    half = bool(getattr(store, "half_grid", False))
    n_slices = store.nslices if half else nmat
    l_static = 0 if half else nmat // 2
    mirror_guard, trs_dev, hgd = None, None, None
    if half:
        # the mirrored slices' own conditioning (HalfGridDressing): checked on
        # the stored chibar BEFORE any W batch overwrites it
        hgd = HalfGridDressing(store.get_freq_batch("chibar", 0, n_slices), S, C, nmat=nmat,
                               spatial_shape=spatial_shape, cond_tol=cond_tol,
                               guard_freqs=guard_freqs, guard_policy=guard_policy)
        mirror_guard, trs_dev = hgd.mirror_guard, hgd.deviation
    cond_s, cond_c = _CondMin(), _CondMin()
    # one list for the whole map: every warned guard violation of either
    # channel and either guard kind (empty under guard_policy = "refuse")
    violations = []
    # accumulated over both channels and every frequency batch (issue #197)
    stats = {"guard_exact_blocks": 0, "guard_blocks": 0}

    def _dress_rows(cb, V_, ch, ls, guard_only=False):
        # one half-grid dressing call: the stored rows (statistics counted)
        # or, guard_only, the diagnostics of their mirrored rows
        return _dress(cb, V_, ch, int(ls[0]), nmat, spatial_shape, cond_tol, iteration,
                      guard_freqs, guard_policy, violations, guard_method,
                      None if guard_only else stats, workers, guard_exact, ls=ls,
                      guard_only=guard_only)

    for l0 in range(0, n_slices, nb):
        l1 = min(n_slices, l0 + nb)
        if half:
            ls = store.grid_index(l0, l1)
            rows, l_first = ls, int(ls[0])
            where = _hg.describe_slices(nmat, l0, l1)
        else:
            ls, rows, l_first = None, slice(l0, l1), l0
            where = "[{}, {})".format(l0, l1)
        # 1 H2D from a host store; on a device store (issue #275) the identity:
        # cb is then a VIEW of the slot that the W batch below overwrites
        # (with "W" aliased onto "chibar"), as on the CPU path
        cb = freq_batch_on(store, "chibar", l0, l1, xp)
        if mirror_guard:
            # time reversal is broken: the mirrored slices l' = nmat - l of
            # this batch (l = 0 and nmat/2 are their own mirrors and are not
            # repeated) go through the same dressing diagnostics as on the
            # full grid -- refusals, warn records, the exact-singular and the
            # "static" residual checks name the TRUE l'. Their scores join
            # cond_min (the minimum over the whole grid), not the guard
            # statistics, which count the stored slices. The temporaries
            # (one mirrored batch and the solve's) are below the stored
            # dressing's own, which runs after them.
            for V_, ch, cmin in ((S, "spin", cond_s), (C, "charge", cond_c)):
                hgd.mirror(_dress_rows, cb, ls, V_, ch, cmin)
        collapse0[rows] = _bk.to_host(cb[:, :, :nd, :nd])
        if half:
            chi_s_b = hgd.stored(_dress_rows, cb, ls, S, "spin", cond_s)
        else:
            chi_s_b, cs = _dress(cb, S, "spin", l_first, nmat, spatial_shape, cond_tol,
                                 iteration, guard_freqs, guard_policy, violations, guard_method,
                                 stats, workers, guard_exact, ls=ls)
            cond_s.add(cs)
        collapse_s[rows] = _bk.to_host(chi_s_b[:, :, :nd, :nd])
        if l0 <= l_static < l1:
            static_s[...] = _bk.to_host(chi_s_b[l_static - l0])
        if output_full:
            store.put_freq_batch("chi_s_w", l0, l1, chi_s_b)
        chi_s_b -= cb
        W_b = 1.5 * _sandwich(S, chi_s_b, S, workers)
        del chi_s_b
        if half:
            chi_c_b = hgd.stored(_dress_rows, cb, ls, C, "charge", cond_c)
        else:
            chi_c_b, cc = _dress(cb, C, "charge", l_first, nmat, spatial_shape, cond_tol,
                                 iteration, guard_freqs, guard_policy, violations, guard_method,
                                 stats, workers, guard_exact, ls=ls)
            cond_c.add(cc)
        collapse_c[rows] = _bk.to_host(chi_c_b[:, :, :nd, :nd])
        if l0 <= l_static < l1:
            static_c[...] = _bk.to_host(chi_c_b[l_static - l0])
        if output_full:
            store.put_freq_batch("chi_c_w", l0, l1, chi_c_b)
        chi_c_b -= cb
        W_b += 0.5 * _sandwich(C, chi_c_b, C, workers)
        del chi_c_b
        if pin is not None:
            W_b = W_b[:, :, pin[:, None], pin[None, :]]
        if bond_w2 is not None:
            for _t, Af, Bf in bond_w2:
                T = Af[None] @ cb
                T = T @ Bf[None]
                T *= 0.5
                W_b += T
                del T
            if not bool(xp.all(xp.isfinite(W_b))):
                raise _NonFiniteError("non-finite effective interaction W in the frequency "
                                      "batch {}{}".format(where, _at(iteration)))
            store.put_freq_batch("W", l0, l1, W_b)
            del W_b, cb
            continue
        # second-order term: A / Bc carry the mixed blocks
        A = _sandwich(S, cb, S, workers)
        Bc = _sandwich(C, cb, C, workers)
        # accumulate_batch is array-module generic; `factors` must carry arrays
        # of xp's module (the caller passes the device pack on the GPU)
        _so.accumulate_batch(W_b[:, :, :nd, :nd], cb[:, :, :nd, :nd], l0, factors)
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
                                  "{}{}".format(where, _at(iteration)))
        # written ONCE, last: when "W" aliases "chibar" (issue #222) this
        # overwrites the batch just consumed, and nothing above reads it again
        store.put_freq_batch("W", l0, l1, W_b)
        del W_b, cb
    if half:
        for coll in (collapse0, collapse_s, collapse_c):
            _hg.fill_mirror(coll, nmat)
    cond_min_s, bracket_s = cond_s.result()
    cond_min_c, bracket_c = cond_c.result()
    return DressResult(collapse0=collapse0, collapse_s=collapse_s, collapse_c=collapse_c,
                       static_s=static_s, static_c=static_c,
                       cond_min_s=cond_min_s, cond_min_c=cond_min_c,
                       cond_bracket_s=bracket_s, cond_bracket_c=bracket_c,
                       guard_violations=len(violations),
                       guard_residual_violations=sum(1 for v in violations
                                                     if v["kind"] == "residual"),
                       guard_exact_blocks=int(stats["guard_exact_blocks"]),
                       guard_blocks=int(stats["guard_blocks"]),
                       mirror_guarded=mirror_guard, time_reversal_deviation=trs_dev)


class _CondMin:
    """The per-map minimum of one channel's guard score over the frequency
    batches (issue #276): the exact minimum over the batches the guard
    decomposed and the bracket over the batches the interval alone passed
    (:class:`~hwave.solver.bond_channels.GuardBracket`)."""

    def __init__(self):
        self.exact = np.inf
        self.bracketed = False
        self.low = np.inf
        self.up = np.inf
        self.index = None

    def add(self, c):
        if c is None:
            return
        if isinstance(c, _bc.GuardBracket):
            self.low = min(self.low, c.low)
            if not self.bracketed or c.up < self.up:
                self.up, self.index = c.up, int(c.index)
            self.bracketed = True
        else:
            self.exact = min(self.exact, c)

    def result(self):
        """``(cond_min, bracket)``: the exact minimum and ``None`` when no
        batch was bracketed, or when no bracketed lower end is below the
        exact minimum -- the minimum is then attained in a decomposed batch
        (the interval encloses every bracketed block's exact score); else
        ``(nan, (low, up, index))``."""
        if not self.bracketed or self.exact <= self.low:
            return float(self.exact), None
        if self.exact < self.up:
            return float("nan"), (float(self.low), float(self.exact), None)
        return float("nan"), (float(self.low), float(self.up), self.index)


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
    to ``xp`` here via :func:`_bk.to_device`, as is each ``W`` block read
    from a host ``store`` (a device store hands its blocks over in place). The returned ``sigma`` is an ``xp`` array; the
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
    get_full = getattr(store, "get_pair_full", None)
    for alpha in range(B):
        Ra = np.asarray(view.delta_r[alpha], dtype=int)
        for bt in range(B):
            Rb = np.asarray(view.delta_r[bt], dtype=int)
            shift = tuple(int(x) for x in (Ra - Rb))   # e^{+ik'.(R_b - R_a)} G(k') = G(r + R_b - R_a)
            if shift == (0, 0, 0):
                G_sh = G_rt.reshape(nmat, nvol, P, P)
            else:
                G_sh = xp.roll(G_rt, shift, axis=axes).reshape(nmat, nvol, P, P)
            if get_full is None:
                src = store.get_pair("W", alpha, bt)
                blk = _bk.to_device(_bk.array_module_of(src).ascontiguousarray(src), xp)
            else:
                # 1 H2D per pair from a host store, none from a device store
                # (issue #275); under the store's half grid (issue #223) the
                # full-grid block is rebuilt from the stored slices of (alpha,
                # beta) and the mirror of those of (beta, alpha)
                blk = get_full("W", alpha, bt, xp=xp)
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


#: the batch floor of the device-store admission under
#: ``longitudinal_bond_store = "auto"`` (issue #275): the store goes on the
#: device only if a frequency batch of at least ``ceil(this / nvol)`` slices
#: (capped by the host table's batch) still fits next to it -- below that a
#: batched solve underfills the device and the host store with a larger
#: batch may be the faster choice. A provisional value, to be revisited with
#: measurements.
DEVICE_STORE_MIN_BLOCKS = 4096

#: the device phase rows, in the order a refusal names them when two of
#: them are exactly equal
_DEV_PHASE_ORDER = ("bubble", "dressing", "transport")


def _largest_dev_phase(rows):
    """Name of the largest device phase row in ``rows``; an exact tie goes
    to the earlier name of :data:`_DEV_PHASE_ORDER`."""
    return max(_DEV_PHASE_ORDER,
               key=lambda k: (rows[k], -_DEV_PHASE_ORDER.index(k)))


def dressing_ops(nmat, nvol, ND, w2_triples=0):
    """Operation count of the two dense ND x ND solves over all (l, q), plus
    the two ND x ND products per spin triple of the exact bond-enlarged
    second order (issue #181) when it is used."""
    return 2.0 * (1 + int(w2_triples)) * nmat * nvol * float(ND) ** 3


def transport_ops(B, nmat, nvol, norb):
    """Operation count of the bond self-energy transport (transforms plus
    the orbital contraction, P^3 = norb^6)."""
    P = norb ** 2
    return float(B * B) * nmat * nvol * (P ** 2 * np.log2(max(nmat * nvol, 2)) + P ** 3)


def estimate_bond_memory(*, nmat, nvol, norb, B, depth, output_full, split_seed, n_types,
                         freq_batch, cap_gb, mixing, factor_bytes=0, device_available=None,
                         guard_freqs="all", guard_enabled=True, w_in_place=False,
                         w2_triples=0, half_grid=False, store="host", pairing=False):
    """The named-buffer lifetime table of spec 3.6 (every row raw), the
    batch selection and the admission decision against ``cap_gb`` (binary
    GiB). Returns a dict with ``persistent_rows``, ``phase_rows`` (at the
    selected ``nb``), ``persistent``, ``nb``, ``peak`` and the symbols;
    raises ``ValueError`` naming every row when even ``nb = 1`` exceeds
    the cap, or when ``freq_batch`` does.

    ``w_in_place`` (issue #222): ``W`` is written over ``chibar`` in ONE
    store slot, so the ``chibar_W`` row is ``U`` instead of ``2 U``. The
    solver sets it exactly when no in-process pairing follows the solve
    (the pairing reads the last map's ``chibar`` back).

    ``factor_bytes`` is the compiled second-order factor pack's
    ``SecondOrderFactors.nbytes`` (0 when none is compiled): the pack is
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
      context, allocated right after the measurement. This is a
      CONSERVATIVE charge, not exact accounting: the resident vertex set is
      now ``S``, ``C`` (``S`` bytes each), the float64 block-weight
      ``mask`` (``ND x ND``) and the int64 pair permutation (``ND``) --
      the on-site vertices it was originally sized for were dropped with
      the retired ``"takimoto"`` kernel (issue #253) -- and on the host
      the same charge also covers the dressing's static slices
      (``static_s`` / ``static_c``, ``2 S``) and the Green-function and
      self-energy copies that the per-phase tracemalloc check of
      ``tests/test_flex_bond_memory.py`` measures. Re-deriving the table
      row by row is a follow-up;
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
    self-energy transport. ``guard(nb) = 4 * nb * nvol * ND^2 * 16`` under
    ``guard_freqs = "all"`` (``4 * nvol * ND^2 * 16``, independent of
    ``nb``, under ``"static"``, where only one frequency slice of ``nvol``
    blocks is ever guarded; ``0`` when the guard is disabled,
    ``guard_enabled=False``) is the conditioning guard's own device
    temporaries (issue #197), incremental to the dressing row rather than
    a phase of its own (the guard runs inside the dressing phase, never
    simultaneously with the bubble or the transport), so it is added to
    ``dressing`` in the need below but never competes for the name a
    refusal gives the largest phase. The factor is 4, not the 3 device
    temporaries active at any one instant, to keep headroom over the
    measured ~3x peak. The device need at a batch size is therefore
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
    which is the conservative reading of the two.

    ``half_grid`` (issue #223): the block store keeps ``nslices = nmat // 2
    + 1`` frequency slices instead of ``nmat``, so the store rows
    (``chibar_W``, ``chi_sc_w``) and the dressing operation count scale by
    ``nslices / nmat`` and the batch is chosen in ``[1, nslices]``. The
    transport still rebuilds one full-grid pair block at a time (its row is
    unchanged), and with ``output_full`` the publish step expands the two
    dynamic channels to the full grid one member at a time, so
    ``post_scf`` holds ``2 U - U_store`` (both full members, plus the half
    of the second while it is expanded) over the persistent rows.

    ``store`` (issue #275, with ``device_available`` only): where the bond
    block store lives -- ``"host"`` (the table above, unchanged), ``"device"``
    or ``"auto"``. A DEVICE store adds its slots as a persistent device row,
    ``block_store = n_slots * U_store`` (``n_slots`` = 1 with ``W`` in place,
    2 with the pairing's separate ``W``, plus 2 with ``output_full``), so the
    device need becomes ``1.25 * (persistent rows + block_store + max(bubble,
    dressing + guard, transport))`` against the same ``0.9 *
    device_available``. ``"auto"`` takes the device store when that need fits
    at the batch floor ``min(host-table nb, ceil(DEVICE_STORE_MIN_BLOCKS /
    nvol))`` (at the pinned ``freq_batch`` when one is given) and the host
    store otherwise; ``"device"`` takes it whenever the need fits at ``nb =
    1`` (or at the pinned batch) and raises ``ValueError`` with the device
    rows otherwise. The dict gains ``store_resident`` (``"host"`` /
    ``"device"``), ``store_bytes``, ``store_reason`` (one line for the log),
    ``device_nb_store`` (the largest batch the device-store table admits, 0
    when none) and, for a device store, ``device_rows["block_store"]``; the
    ``device_*`` keys then describe the device-store table. The host table
    is not changed by the residency (the host preflight runs before the
    backend is known, so it keeps charging the store on the host: a
    deliberately conservative reading).

    ``pairing`` (in-process pairing on): a device store also charges the
    pairing's re-dressing, ``pairing_redress = norb^4 nvol Nmat 16 + 2 nb
    nvol ND^2 16`` (the pair bubble plus two batch temporaries), a device
    row of the device-store table only."""
    if store not in ("host", "device", "auto"):
        raise ValueError("estimate_bond_memory: store must be 'host', 'device' or 'auto', "
                         "got {!r}".format(store))
    nmat, nvol, norb, B = int(nmat), int(nvol), int(norb), int(B)
    depth = max(1, int(depth)) if mixing == "anderson" else 0      # the mixer's effective depth
    it = 16
    P = norb * norb
    ND = B * P
    U = nmat * nvol * ND * ND * it
    # the frequency slices the block store keeps (issue #223)
    nslices = (nmat // 2 + 1) if half_grid else nmat
    U_store = nslices * nvol * ND * ND * it
    G = nmat * nvol * P * it
    C = nmat * nvol * P * P * it
    S = nvol * ND * ND * it
    H = nvol * P * it
    # the static vertex stack: a conservative 5 S, see the docstring (the
    # resident set is S, C, mask and perm; the margin also covers the static
    # slices and the Green-function copies; re-deriving it is a follow-up)
    V = 5 * S
    persistent_rows = {
        "chibar_W": (1 if w_in_place else 2) * U_store,
        "chi_sc_w": 2 * U_store if output_full else 0,
        "vertices_static": V,
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
            "post_scf": G + ((2 * U - U_store if half_grid else U) if output_full else C),
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
        # a batch is never longer than the stored frequency axis
        nb = min(nb, nslices)
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
        for cand in range(nslices, 0, -1):
            if _peak(cand) <= cap_bytes:
                nb = cand
                break

    device = {}
    if device_available is not None:
        vertices = V
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
            # see _largest_dev_phase, which never names it. Under
            # guard_freqs = "static" only one frequency slice of nvol
            # blocks is ever guarded per batch, so the row does not scale
            # with n there; a disabled guard (cond_tol = None) allocates
            # nothing.
            rows["guard"] = (4 * (n if guard_freqs == "all" else 1) * nvol * ND * ND * it
                             if guard_enabled else 0)
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
                "{:.4f} GiB (frequency batch 1, phase '{}') = 1.25 * (persistent rows + max(bubble, "
                "dressing + guard, transport)) exceeds 0.9 * the available device memory {:.4f} "
                "GiB; the rows are\n{}\nReduce "
                "the k mesh or Nmat, drop declared-zero outer shells with "
                "longitudinal_bond_max_shells, or run with gpu=false.".format(
                    _dev_need(1) / _GIB, phase, dev_cap / _GIB, _dev_table(1)))
        if freq_batch is not None and _dev_need(nb) > dev_cap:
            raise ValueError(
                "[mode.param] longitudinal_bond_freq_batch = {} gives an estimated device need of "
                "{:.4f} GiB above 0.9 * the available device memory {:.4f} GiB (host peak {:.4f} GiB "
                "against the host cap {:.4f} GiB); the device rows are\n{}".format(
                    int(freq_batch), _dev_need(nb) / _GIB, dev_cap / _GIB,
                    _peak(nb) / _GIB, cap_bytes / _GIB, _dev_table(nb)))
        device_nb = 1
        for cand in range(nslices, 0, -1):
            if _dev_need(cand) <= dev_cap:
                device_nb = cand
                break
        # the device-resident block store (issue #275): every slot of the
        # store is one persistent device row on top of the host-store table
        n_slots = (1 if w_in_place else 2) + (2 if output_full else 0)
        store_bytes = n_slots * U_store

        def _pair_redress(n):
            # the in-process pairing's re-dressing (issue #275): with a host
            # store it ran next to the device memory the store left free; with
            # the store resident it is charged here -- the pair bubble G2
            # (norb^4 nvol nmat, the uniform grid's frequency count, an upper
            # bound for the IR nodes) plus two batch temporaries (the S chi S
            # contribution and its scaled copy) on top of the dressing row
            if not pairing:
                return 0
            return norb ** 4 * nvol * nmat * it + 2 * n * nvol * ND * ND * it

        def _dev_rows_store(n):
            rows = _dev_rows(n)
            rows["block_store"] = store_bytes
            if pairing:
                rows["pairing_redress"] = _pair_redress(n)
            return rows

        def _dev_need_store(n):
            return _dev_need(n) + 1.25 * (store_bytes + _pair_redress(n))

        def _dev_table_store(n):
            return "\n".join("  device     {:>18s}: {:10.4f} GiB".format(k, v / _GIB)
                              for k, v in _dev_rows_store(n).items())
        device_nb_store = 0
        for cand in range(nslices, 0, -1):
            if _dev_need_store(cand) <= dev_cap:
                device_nb_store = cand
                break
        nb_host_store = min(nb, device_nb)
        resident, reason = "host", None
        if store == "host":
            reason = "longitudinal_bond_store = \"host\""
        else:
            if freq_batch is not None:
                at = nb
            elif store == "device":
                at = 1
            else:
                at = min(nb, max(1, -(-DEVICE_STORE_MIN_BLOCKS // nvol)))
            fits = _dev_need_store(at) <= dev_cap
            what = ("{:.4f} GiB = 1.25 * (persistent rows + block_store {:.4f} GiB + max phase) "
                    "at frequency batch {} against 0.9 * available {:.4f} GiB".format(
                        _dev_need_store(at) / _GIB, store_bytes / _GIB, at, dev_cap / _GIB))
            if fits:
                resident = "device"
                reason = "the device need with the store resident is " + what
            elif store == "device":
                raise ValueError(
                    "[mode.param] longitudinal_bond_store = \"device\": the bond block store "
                    "({:.4f} GiB, {} slot(s) of {:.4f} GiB) does not fit on the device: the "
                    "estimated device need is {}; the device rows are\n{}\nUse "
                    "longitudinal_bond_store = \"auto\" or \"host\", enable "
                    "longitudinal_bond_half_grid, or reduce the k mesh or Nmat.".format(
                        store_bytes / _GIB, n_slots, U_store / _GIB, what,
                        _dev_table_store(at)))
            else:
                reason = "the store does not fit on the device: the need with it resident " \
                         "would be " + what
        if resident == "device":
            nb = min(nb, device_nb_store)
            device = dict(device_rows=_dev_rows_store(nb), device_need=_dev_need_store(nb),
                          device_cap=dev_cap, device_nb=device_nb_store,
                          device_table=_dev_table_store(nb))
        else:
            nb = nb_host_store
            device = dict(device_rows=_dev_rows(nb), device_need=_dev_need(nb),
                          device_cap=dev_cap, device_nb=device_nb, device_table=_dev_table(nb))
        device.update(store_resident=resident, store_bytes=store_bytes, store_reason=reason,
                      device_nb_store=device_nb_store)

    phase_rows = _phase_rows(nb)
    return dict(persistent_rows=persistent_rows, phase_rows=phase_rows, persistent=persistent,
                nb=nb, peak=_peak(nb), cap_bytes=cap_bytes, U=U, G_bytes=G, C_bytes=C,
                S_bytes=S, H_bytes=H, B=B, ND=ND, nvol=nvol, nmat=nmat, table=_table(nb),
                nslices=nslices, U_store=U_store, half_grid=bool(half_grid),
                dressing_ops=dressing_ops(nslices, nvol, ND, w2_triples),
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
