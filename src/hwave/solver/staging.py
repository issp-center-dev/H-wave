"""Host-side copies and bounded pinned staging for the bond block store
(GitHub issue #275, part 2).

Two helpers move bytes without changing any of them:

* :func:`parallel_copy` -- ``dst[...] = src`` for host arrays, split over the
  leading axis among a few threads. numpy releases the GIL for a copy of a
  numeric dtype, so the pair-strided writes and reads of a host store
  (``put_pair`` / ``get_pair_full``) no longer run on one core. Small copies
  stay a single assignment.
* :class:`PinnedStager` -- host <-> device copies through ONE reusable
  page-locked (pinned) staging buffer of bounded size, in chunks: ranges
  along the leading axis, or -- when one leading-axis slice is larger than
  the buffer -- along the first later axis whose sub-arrays fit, with the
  axes before it fixed. A copy from pageable memory is staged by the driver through
  its own small pinned buffers (and ``cupy.asnumpy`` of a whole store slot
  also allocates and first-touches a full-size pageable result in one
  thread); a chunk at a time through a pinned buffer runs at the bus rate,
  and the host side of each chunk is a :func:`parallel_copy`. Every chunk
  writes straight into its destination, so no full-size temporary is made.

When the pinned buffer cannot be allocated (or the module has no pinned
allocator), or one array element is larger than the buffer, THAT copy takes
the previous path (``asnumpy`` / ``asarray`` of the whole array), with a
debug message for a failed allocation; the next copy tries again. The bytes
that arrive are the same either way.
"""
import logging
import os
import threading
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from . import backend as _bk

logger = logging.getLogger("qlms.solver.staging")

#: the default size of the pinned staging buffer (bytes): large enough that
#: a chunk's fixed cost is negligible, small enough that page-locking it is
#: harmless next to a store of tens of GiB
PINNED_STAGING_BYTES = 256 * 1024 ** 2

#: host copies below this many bytes are one plain assignment
PARALLEL_COPY_MIN_BYTES = 16 * 1024 ** 2

#: the default number of copy threads (memory bandwidth saturates early)
COPY_WORKERS = 8

_pool = None
_pool_lock = threading.Lock()


def _usable_cpus():
    """CPUs this process may run on: its affinity mask where the platform
    reports one (a batch job's cpuset), else ``os.cpu_count()``."""
    try:
        return len(os.sched_getaffinity(0)) or 1
    except (AttributeError, OSError):
        return os.cpu_count() or 1


def default_copy_workers():
    """``min(COPY_WORKERS, usable CPUs)``: a job given few CPUs is not
    oversubscribed."""
    return max(1, min(COPY_WORKERS, _usable_cpus()))


def _executor():
    """The one copy pool of the process (``COPY_WORKERS`` threads), created
    on first use and never shut down, so a concurrent caller cannot find it
    closed between obtaining it and submitting to it."""
    global _pool
    with _pool_lock:
        if _pool is None:
            _pool = ThreadPoolExecutor(max_workers=max(1, COPY_WORKERS),
                                       thread_name_prefix="hwave-copy")
        return _pool


def parallel_copy(dst, src, workers=None, min_bytes=None):
    """``dst[...] = src`` for HOST arrays of the same shape, the leading axis
    cut into up to ``workers`` contiguous pieces copied by threads. Bitwise
    the plain assignment (which it is for small or overlapping operands).
    Returns ``dst``."""
    if tuple(dst.shape) != tuple(src.shape):
        raise ValueError("parallel_copy: shapes differ: {} vs {}".format(dst.shape, src.shape))
    w = default_copy_workers() if workers is None else max(1, int(workers))
    limit = PARALLEL_COPY_MIN_BYTES if min_bytes is None else int(min_bytes)
    n = int(dst.shape[0]) if dst.ndim else 0
    w = min(w, n)
    # overlapping operands keep numpy's overlap-safe single assignment
    if w <= 1 or int(dst.nbytes) < limit or np.may_share_memory(dst, src):
        dst[...] = src
        return dst
    edges = [n * i // w for i in range(w + 1)]

    def piece(i):
        dst[edges[i]:edges[i + 1]] = src[edges[i]:edges[i + 1]]
    # list() re-raises the first worker exception here
    list(_executor().map(piece, range(w)))
    return dst


class PinnedStager:
    """Host <-> device copies through one reusable pinned buffer of at most
    ``max_bytes`` bytes (allocated on first use,
    grown only up to ``max_bytes``, freed by :meth:`release`).

    ``xp`` is the DEVICE module (cupy, or a test double): the buffer comes
    from ``xp.cuda.alloc_pinned_memory``, a device chunk moves to the host by
    ``chunk.get(out=buffer)`` and back by ``chunk.set(buffer)``. A chunk is a
    range along one axis ``k`` with every axis before it fixed, ``k`` being
    the first axis whose sub-arrays fit in ``max_bytes`` (axis 0 unless one
    frequency slice is larger than the buffer). For a contiguous device
    array the device side of each chunk is therefore contiguous; the host
    side may be any strided view (a channel-pair block of a host store).

    Every chunk copy is issued on the CURRENT stream (captured per call and
    passed explicitly) and that stream is synchronized after it, before the
    staging buffer is read (device -> host) or overwritten or released
    (host -> device): with a pinned source, ``set`` is truly asynchronous,
    and ``get`` is non-blocking on older CuPy or other streams.

    A failed pinned allocation sends THAT call down the unstaged path (with a
    debug message); the next call tries again. The copies hold a lock, so
    concurrent callers sharing one stager take turns with the buffer."""

    def __init__(self, xp, max_bytes=None, workers=None):
        self.xp = xp
        self.max_bytes = int(PINNED_STAGING_BYTES if max_bytes is None else max_bytes)
        if self.max_bytes < 1:
            raise ValueError("PinnedStager: max_bytes must be positive")
        self.workers = workers
        self._mem = None
        self._owner = None
        self._nbytes = 0
        self._failed = False     # the LAST allocation attempt failed
        self._lock = threading.RLock()
        self.chunks = 0          # chunk copies made (for tests and logs)

    @property
    def nbytes(self):
        """Bytes of pinned memory currently held."""
        return self._nbytes

    @property
    def available(self):
        """False when the last pinned allocation failed (it is retried)."""
        return not self._failed

    def release(self):
        """Drop the pinned buffer and hand cached pinned blocks back."""
        with self._lock:
            held = self._mem is not None
            self._mem = self._owner = None
            self._nbytes = 0
            if held:
                try:
                    self.xp.get_default_pinned_memory_pool().free_all_blocks()
                except Exception:            # noqa: BLE001 - only a cache hint
                    pass

    def _buffer(self, nbytes):
        """A pinned uint8 buffer of at least ``nbytes``, or None (fallback)."""
        if self._mem is not None and self._nbytes >= nbytes:
            return self._mem
        self.release()
        try:
            mem = self.xp.cuda.alloc_pinned_memory(int(nbytes))
            buf = np.ndarray((int(nbytes),), dtype=np.uint8, buffer=mem)
        except Exception as exc:         # noqa: BLE001 - any failure: old path
            self._failed = True
            logger.debug("pinned staging buffer of %d bytes unavailable (%s: %s); this host "
                         "<-> device copy of the bond block store takes the unstaged path",
                         int(nbytes), type(exc).__name__, exc)
            return None
        self._failed = False
        # keep the allocation alive with the view
        self._mem, self._nbytes, self._owner = buf, int(nbytes), mem
        return buf

    def _plan(self, shape, dtype):
        """``(chunks, buf)``: the index tuples of the chunks (None for an
        empty array) and the buffer (None: the unstaged path)."""
        shape = tuple(int(x) for x in shape)
        item = np.dtype(dtype).itemsize
        if not shape or 0 in shape:
            return None, None
        nd = len(shape)
        k = 0
        while k < nd - 1 and int(np.prod(shape[k + 1:], dtype=np.int64)) * item > self.max_bytes:
            k += 1
        inner = int(np.prod(shape[k + 1:], dtype=np.int64)) * item
        if inner > self.max_bytes:
            return None, None            # not even one element fits: unstaged
        per = max(1, min(shape[k], self.max_bytes // inner))
        buf = self._buffer(per * inner)
        if buf is None:
            return None, None
        chunks = [prefix + (slice(i0, min(shape[k], i0 + per)),)
                  for prefix in np.ndindex(*shape[:k])
                  for i0 in range(0, shape[k], per)]
        return chunks, buf

    @staticmethod
    def _view(buf, shape, dtype):
        cnt = int(np.prod(shape, dtype=np.int64))
        return buf[:cnt * np.dtype(dtype).itemsize].view(dtype).reshape(shape)

    def copy_to_host(self, src, dst):
        """``dst[...] = src`` for a DEVICE ``src`` and a HOST ``dst`` (any
        strides) of the same shape and dtype."""
        if tuple(src.shape) != tuple(dst.shape) or src.dtype != dst.dtype:
            raise ValueError("PinnedStager.copy_to_host: {} {} vs {} {}".format(
                src.shape, src.dtype, dst.shape, dst.dtype))
        if src.ndim == 0:
            dst[...] = _bk.to_host(src)
            return dst
        with self._lock:
            chunks, buf = self._plan(src.shape, src.dtype)
            if buf is None:
                parallel_copy(dst, _bk.to_host(src), self.workers)
                return dst
            stream = self.xp.cuda.get_current_stream()
            for idx in chunks:
                part = src[idx]
                stage = self._view(buf, part.shape, src.dtype)
                part.get(stream=stream, out=stage)
                stream.synchronize()             # the chunk has landed in stage
                parallel_copy(dst[idx], stage, self.workers)
                self.chunks += 1
            return dst

    def to_host(self, src):
        """A new host array equal to the device array ``src``."""
        if src.ndim == 0:
            return _bk.to_host(src)
        with self._lock:
            _chunks, buf = self._plan(src.shape, src.dtype)
            if buf is None:
                return _bk.to_host(src)          # no full-size double on the old path
            dst = np.empty(src.shape, dtype=src.dtype)
            return self.copy_to_host(src, dst)

    def to_device(self, src):
        """A new contiguous device array equal to the host array ``src`` (any
        strides): each chunk is gathered into the pinned buffer and set into
        its part of the result."""
        xp = self.xp
        if src.ndim == 0:
            return xp.asarray(src)
        with self._lock:
            chunks, buf = self._plan(src.shape, src.dtype)
            if buf is None:
                return xp.asarray(np.ascontiguousarray(src))
            dst = xp.empty(src.shape, dtype=src.dtype)
            stream = xp.cuda.get_current_stream()
            for idx in chunks:
                part = src[idx]
                stage = self._view(buf, part.shape, src.dtype)
                parallel_copy(stage, part, self.workers)
                dst[idx].set(stage, stream=stream)
                stream.synchronize()             # stage may be overwritten / freed now
                self.chunks += 1
            return dst
