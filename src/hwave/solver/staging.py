"""Host-side copies and bounded pinned staging for the bond block store
(GitHub issue #275, part 2).

Two helpers move bytes without changing any of them:

* :func:`parallel_copy` -- ``dst[...] = src`` for host arrays, split over the
  leading axis among a few threads. numpy releases the GIL for a copy of a
  numeric dtype, so the pair-strided writes and reads of a host store
  (``put_pair`` / ``get_pair_full``) no longer run on one core. Small copies
  stay a single assignment.
* :class:`PinnedStager` -- host <-> device copies through ONE reusable
  page-locked (pinned) staging buffer of bounded size, in chunks along the
  leading axis. A copy from pageable memory is staged by the driver through
  its own small pinned buffers (and ``cupy.asnumpy`` of a whole store slot
  also allocates and first-touches a full-size pageable result in one
  thread); a chunk at a time through a pinned buffer runs at the bus rate,
  and the host side of each chunk is a :func:`parallel_copy`. Every chunk
  writes straight into its destination, so no full-size temporary is made.

When the pinned buffer cannot be allocated (or the module has no pinned
allocator) the stager says so once at debug level and every copy takes the
previous path (``asnumpy`` / ``asarray`` of the whole array). The bytes that
arrive are the same either way.
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


def default_copy_workers():
    return max(1, min(COPY_WORKERS, os.cpu_count() or 1))


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
    ``max(max_bytes, one leading-axis row)`` bytes (allocated on first use,
    grown only when a row is larger, freed by :meth:`release`).

    ``xp`` is the DEVICE module (cupy, or a test double): the buffer comes
    from ``xp.cuda.alloc_pinned_memory``, a device chunk moves to the host by
    ``chunk.get(out=buffer)`` and back by ``chunk.set(buffer)``. Chunks run
    along axis 0 of the array, so the device side of each chunk is a
    contiguous slice of a contiguous array and the host side may be any
    strided view (a channel-pair block of a host store).

    Every chunk copy is issued on the CURRENT stream (captured per call and
    passed explicitly) and that stream is synchronized after it, before the
    staging buffer is read (device -> host) or overwritten or released
    (host -> device): with a pinned source, ``set`` is truly asynchronous,
    and ``get`` is non-blocking on older CuPy or other streams."""

    def __init__(self, xp, max_bytes=None, workers=None):
        self.xp = xp
        self.max_bytes = int(PINNED_STAGING_BYTES if max_bytes is None else max_bytes)
        if self.max_bytes < 1:
            raise ValueError("PinnedStager: max_bytes must be positive")
        self.workers = workers
        self._mem = None
        self._owner = None
        self._nbytes = 0
        self._failed = False
        self.chunks = 0          # chunk copies made (for tests and logs)

    @property
    def nbytes(self):
        """Bytes of pinned memory currently held."""
        return self._nbytes

    @property
    def available(self):
        return not self._failed

    def release(self):
        """Drop the pinned buffer and hand cached pinned blocks back."""
        held = self._mem is not None
        self._mem = self._owner = None
        self._nbytes = 0
        if held:
            try:
                self.xp.get_default_pinned_memory_pool().free_all_blocks()
            except Exception:            # noqa: BLE001 - only a cache hint
                pass

    def _rows(self, shape, dtype):
        row = int(np.prod(shape[1:], dtype=np.int64)) * np.dtype(dtype).itemsize
        return max(1, row), max(1, self.max_bytes // max(1, row))

    def _buffer(self, nbytes):
        """A pinned uint8 buffer of at least ``nbytes``, or None (fallback)."""
        if self._failed:
            return None
        if self._mem is not None and self._nbytes >= nbytes:
            return self._mem
        self.release()
        try:
            mem = self.xp.cuda.alloc_pinned_memory(int(nbytes))
            buf = np.ndarray((int(nbytes),), dtype=np.uint8, buffer=mem)
        except Exception as exc:         # noqa: BLE001 - any failure: old path
            self._failed = True
            logger.debug("pinned staging buffer of %d bytes unavailable (%s: %s); host <-> "
                         "device copies of the bond block store take the unstaged path",
                         int(nbytes), type(exc).__name__, exc)
            return None
        # keep the allocation alive with the view
        self._mem, self._nbytes, self._owner = buf, int(nbytes), mem
        return buf

    def _plan(self, shape, dtype):
        n = int(shape[0]) if len(shape) else 0
        row, per = self._rows(shape, dtype)
        per = min(per, max(1, n))
        buf = self._buffer(per * row) if n else None
        return n, per, buf

    @staticmethod
    def _view(buf, n, shape, dtype):
        cnt = n * int(np.prod(shape[1:], dtype=np.int64))
        return buf[:cnt * np.dtype(dtype).itemsize].view(dtype).reshape((n,) + tuple(shape[1:]))

    def copy_to_host(self, src, dst):
        """``dst[...] = src`` for a DEVICE ``src`` and a HOST ``dst`` (any
        strides) of the same shape and dtype."""
        if tuple(src.shape) != tuple(dst.shape) or src.dtype != dst.dtype:
            raise ValueError("PinnedStager.copy_to_host: {} {} vs {} {}".format(
                src.shape, src.dtype, dst.shape, dst.dtype))
        if src.ndim == 0:
            dst[...] = _bk.to_host(src)
            return dst
        n, per, buf = self._plan(src.shape, src.dtype)
        if buf is None:
            parallel_copy(dst, _bk.to_host(src), self.workers)
            return dst
        stream = self.xp.cuda.get_current_stream()
        for i0 in range(0, n, per):
            i1 = min(n, i0 + per)
            stage = self._view(buf, i1 - i0, src.shape, src.dtype)
            src[i0:i1].get(stream=stream, out=stage)
            stream.synchronize()             # the chunk has landed in stage
            parallel_copy(dst[i0:i1], stage, self.workers)
            self.chunks += 1
        return dst

    def to_host(self, src):
        """A new host array equal to the device array ``src``."""
        if src.ndim == 0:
            return _bk.to_host(src)
        _n, _per, buf = self._plan(src.shape, src.dtype)
        if buf is None:
            return _bk.to_host(src)          # no full-size double on the old path
        dst = np.empty(src.shape, dtype=src.dtype)
        return self.copy_to_host(src, dst)

    def to_device(self, src):
        """A new contiguous device array equal to the host array ``src`` (any
        strides): each chunk is gathered into the pinned buffer and set into
        its slice of the result."""
        xp = self.xp
        if src.ndim == 0:
            return xp.asarray(src)
        n, per, buf = self._plan(src.shape, src.dtype)
        if buf is None:
            return xp.asarray(np.ascontiguousarray(src))
        dst = xp.empty(src.shape, dtype=src.dtype)
        stream = xp.cuda.get_current_stream()
        for i0 in range(0, n, per):
            i1 = min(n, i0 + per)
            stage = self._view(buf, i1 - i0, src.shape, src.dtype)
            parallel_copy(stage, src[i0:i1], self.workers)
            dst[i0:i1].set(stage, stream=stream)
            stream.synchronize()             # stage may be overwritten / freed now
            self.chunks += 1
        return dst
