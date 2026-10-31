"""Thread-parallel stacked small-matrix linear algebra (bond dressing).

numpy's batched ``solve`` / ``svd`` / ``matmul`` loop over the stack in one
thread, and a multithreaded BLAS makes tiny matrices slower, not faster. Here
the stack is cut into contiguous chunks that worker threads process with the
BLAS pinned to one thread. Every matrix is independent, so the result equals
the serial numpy result up to rounding (the last digit can differ when the
serial reference runs on a multithreaded BLAS). ``workers`` of ``None`` or 1 is the plain
numpy call; without ``threadpoolctl`` a request for more falls back to it
(with one warning), since threads over a multithreaded BLAS are unsafe.
"""
import logging
import os
import threading
from concurrent.futures import ThreadPoolExecutor, wait

import numpy as np

try:
    from threadpoolctl import threadpool_limits
except ImportError:  # pragma: no cover - exercised by the environment
    threadpool_limits = None

logger = logging.getLogger("qlms.solver.threaded")

#: stacks shorter than this per chunk are not worth a thread
_MIN_ROWS = 256

_pools = {}
_pools_lock = threading.Lock()
# threadpool_limits is process-global: one caller at a time owns the pinned region
_run_lock = threading.Lock()
_warned = False
_warned_lock = threading.Lock()


def available():
    return threadpool_limits is not None


def _resolve_workers(workers):
    if workers is None:
        return 1
    workers = int(workers)
    if workers == -1:
        return os.cpu_count() or 1
    if workers < 1:
        raise ValueError("workers must be None, -1 or a positive integer; got {}".format(workers))
    return workers


def _effective_workers(workers):
    global _warned
    w = _resolve_workers(workers)
    if w > 1 and not available():
        with _warned_lock:
            first, _warned = not _warned, True
        if first:
            logger.warning("threaded linear algebra requested (workers=%d) but threadpoolctl "
                           "is not installed; running serially (install it with "
                           "pip install \"hwave[threads]\")", w)
        return 1
    return w


def _pool(w):
    with _pools_lock:
        if w not in _pools:
            _pools[w] = ThreadPoolExecutor(max_workers=w, thread_name_prefix="hwave-la")
        return _pools[w]


def _split(n, nchunks):
    nchunks = max(1, min(nchunks, n))
    edges = [n * i // nchunks for i in range(nchunks + 1)]
    return [slice(edges[i], edges[i + 1]) for i in range(nchunks)]


def _nchunks(rows, workers, min_rows):
    if min_rows is None:
        min_rows = _MIN_ROWS
    return max(1, min(workers, rows // max(1, int(min_rows))))


def _row_slices(n, workers, min_rows):
    return _split(n, _nchunks(n, workers, min_rows))


def _run(tasks, w):
    with _run_lock, threadpool_limits(limits=1, user_api="blas"):
        futures = [_pool(w).submit(t) for t in tasks]
        # join every future before the BLAS limit is restored, then re-raise the first error
        wait(futures)
        return [f.result() for f in futures]


def _take(arr, axis, sl):
    if arr.shape[axis] == 1:
        return arr
    idx = [slice(None)] * arr.ndim
    idx[axis] = sl
    return arr[tuple(idx)]


def _align(a, b):
    nd = max(a.ndim, b.ndim)
    return (a.reshape((1,) * (nd - a.ndim) + a.shape),
            b.reshape((1,) * (nd - b.ndim) + b.shape))


def _batch_axis(shape, workers):
    """Axis (among the leading batch axes) to cut: the first one long enough
    for every worker, else the longest."""
    lead = shape[:-2]
    for ax, n in enumerate(lead):
        if n >= workers:
            return ax
    return max(range(len(lead)), key=lambda ax: lead[ax])


def solve(a, b, *, workers=None, min_rows=None):
    w = _effective_workers(workers)
    a = np.asarray(a)
    b = np.asarray(b)
    if w == 1 or a.ndim != 3 or b.ndim != 3 or a.shape[0] != b.shape[0]:
        return np.linalg.solve(a, b)
    sl = _row_slices(a.shape[0], w, min_rows)
    if len(sl) == 1:
        return np.linalg.solve(a, b)
    parts = _run([(lambda s=s: np.linalg.solve(a[s], b[s])) for s in sl], w)
    return np.concatenate(parts, axis=0)


def svdvals(a, *, workers=None, min_rows=None):
    w = _effective_workers(workers)
    a = np.asarray(a)
    if w == 1 or a.ndim != 3:
        return np.linalg.svd(a, compute_uv=False)
    sl = _row_slices(a.shape[0], w, min_rows)
    if len(sl) == 1:
        return np.linalg.svd(a, compute_uv=False)
    parts = _run([(lambda s=s: np.linalg.svd(a[s], compute_uv=False)) for s in sl], w)
    return np.concatenate(parts, axis=0)


def matmul(a, b, *, workers=None, min_rows=None):
    w = _effective_workers(workers)
    a = np.asarray(a)
    b = np.asarray(b)
    if w == 1 or a.ndim < 2 or b.ndim < 2 or (a.ndim < 3 and b.ndim < 3):
        return a @ b
    a, b = _align(a, b)
    lead = np.broadcast_shapes(a.shape[:-2], b.shape[:-2])
    rows = int(np.prod(lead))
    nch = _nchunks(rows, w, min_rows)
    if nch == 1:
        return a @ b
    ax = _batch_axis(lead + (0, 0), nch)
    sl = _split(lead[ax], nch)
    out = np.empty(lead + (a.shape[-2], b.shape[-1]), dtype=np.result_type(a, b))

    def work(s):
        o = [slice(None)] * out.ndim
        o[ax] = s
        np.matmul(_take(a, ax, s), _take(b, ax, s), out=out[tuple(o)])

    _run([(lambda s=s: work(s)) for s in sl], w)
    return out


def sandwich(left, x, right, *, workers=None, min_rows=None):
    """``left[None] @ x @ right[None]`` (left-associated), ``left`` and
    ``right`` `(nvol, n, n)`, ``x`` `(nb, nvol, n, n)`."""
    w = _effective_workers(workers)
    left = np.asarray(left)
    x = np.asarray(x)
    right = np.asarray(right)
    L, R = left[None], right[None]
    nch = _nchunks(x.shape[0] * x.shape[1], w, min_rows) if w > 1 else 1
    if nch == 1:
        return L @ x @ R
    ax = 0 if x.shape[0] >= nch else 1
    sl = _split(x.shape[ax], nch)
    out = np.empty(x.shape, dtype=np.result_type(L, x, R))

    def work(s):
        o = [slice(None)] * out.ndim
        o[ax] = s
        out[tuple(o)] = _take(L, ax, s) @ x[tuple(o[:2])] @ _take(R, ax, s)

    _run([(lambda s=s: work(s)) for s in sl], w)
    return out
