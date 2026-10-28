"""A fake "cupy" for CPU-only tests of the device-resident code paths
(issue #275).

:func:`fake_cupy` patches :func:`hwave.solver.backend._import_cupy` so that
the backend sees a GPU: arrays of :class:`DevArray` (an ``ndarray`` subclass
whose class lives in a module named ``cupy...``) are "device" arrays for
:func:`hwave.solver.backend.array_module_of`, and the returned module
:class:`FakeCupy` proxies numpy, returning device arrays from every function.
The arithmetic is numpy's, so a device run and a host run of the same
kernels are bitwise comparable; what the fake adds is the BOUNDARY:

* every host <-> device copy is recorded in :attr:`FakeCupy.transfers`
  (``asarray`` / ``array`` / ``ascontiguousarray`` of a host array, and a
  host array assigned into a device array, are host-to-device;
  ``asnumpy`` is device-to-host);
* a host ``ndarray`` (of at least one dimension) mixed into a ufunc or a
  module function with device arrays raises ``TypeError``, as cupy does.

What it cannot catch: ``np.asarray`` / ``np.array`` of a device array and a
device array assigned into a host array succeed silently here (the subclass
IS an ndarray), where cupy refuses the implicit conversion. Such call sites
are caught by the transfer accounting (no ``asnumpy`` is recorded for them)
and by the real-GPU tests of ``tests/test_flex_bond_gpu.py``.
"""
import contextlib
import types
from unittest import mock

import numpy as np


class DevArray(np.ndarray):
    """A "device" array: a numpy array whose type is in a ``cupy`` module."""

    def __array_ufunc__(self, ufunc, method, *inputs, out=None, **kw):
        _refuse_host(inputs, "ufunc {}".format(ufunc.__name__))
        args = tuple(_base(x) for x in inputs)
        if out is not None:
            _refuse_host(out, "ufunc {} (out)".format(ufunc.__name__))
            kw["out"] = tuple(_base(o) for o in out)
        res = getattr(ufunc, method)(*args, **kw)
        if out is not None:
            return out[0] if len(out) == 1 else tuple(out)
        return _wrap(res)

    def __setitem__(self, key, value):
        if _is_host(value):
            _MODULE.record("h2d", value)
        np.ndarray.__setitem__(self, key, value)


DevArray.__module__ = "cupy._fake_device"


def _is_host(x):
    return isinstance(x, np.ndarray) and not isinstance(x, DevArray) and x.ndim > 0


def _base(x):
    return x.view(np.ndarray) if isinstance(x, DevArray) else x


def _refuse_host(args, what):
    for x in args:
        items = x if isinstance(x, (list, tuple)) else (x,)
        for y in items:
            if _is_host(y):
                raise TypeError("fake cupy: {} got a host numpy.ndarray operand of shape {} "
                                "(cupy refuses implicit host arrays)".format(what, y.shape))


def _wrap(res):
    if isinstance(res, DevArray):
        return res
    if isinstance(res, np.ndarray):
        return res.view(DevArray)
    if isinstance(res, tuple):
        items = [_wrap(r) for r in res]
        return type(res)(*items) if hasattr(res, "_fields") else tuple(items)
    if isinstance(res, list):
        return [_wrap(r) for r in res]
    return res


class _Namespace:
    """numpy (or one of its submodules) with every function's array results
    returned as device arrays and host array arguments refused."""

    def __init__(self, target, name):
        self._target = target
        self.__name__ = name

    def __getattr__(self, attr):
        obj = getattr(self._target, attr)
        if isinstance(obj, types.ModuleType):
            ns = _Namespace(obj, "{}.{}".format(self.__name__, attr))
            setattr(self, attr, ns)
            return ns
        if isinstance(obj, type) or not callable(obj):
            return obj

        def fn(*a, **k):
            _refuse_host(a, "{}.{}".format(self.__name__, attr))
            _refuse_host(tuple(v for v in k.values() if isinstance(v, np.ndarray)),
                         "{}.{}".format(self.__name__, attr))
            return _wrap(obj(*a, **k))
        fn.__name__ = attr
        return fn


class _OutOfMemoryError(MemoryError):
    """The fake ``cupy.cuda.memory.OutOfMemoryError``."""


class _Pool:
    def __init__(self, module):
        self._m = module

    def used_bytes(self):
        return 0

    def total_bytes(self):
        return 0

    def free_all_blocks(self):
        self._m.pool_frees += 1


class FakeCupy(_Namespace):
    """The fake module. ``free`` is what ``cuda.runtime.memGetInfo`` reports;
    ``fail_zeros`` makes the next ``zeros`` calls whose request is at least
    ``fail_zeros_min`` bytes raise the out-of-memory error (a count)."""

    def __init__(self, free=2 ** 40):
        super().__init__(np, "cupy")
        self.free = int(free)
        self.transfers = []
        self.pool_frees = 0
        self.fail_zeros = 0
        self.fail_zeros_min = 0
        self.zeros_calls = []
        runtime = types.SimpleNamespace(getDeviceCount=lambda: 1,
                                        memGetInfo=lambda: (self.free, 2 * self.free))
        memory = types.SimpleNamespace(OutOfMemoryError=_OutOfMemoryError)
        self.cuda = types.SimpleNamespace(runtime=runtime, memory=memory)
        self._pool = _Pool(self)

    def get_default_memory_pool(self):
        return self._pool

    # -- transfers ----------------------------------------------------------
    def record(self, direction, arr):
        self.transfers.append((direction, tuple(arr.shape), int(arr.nbytes)))

    def reset(self):
        self.transfers.clear()

    def asarray(self, a, dtype=None, order=None):
        if isinstance(a, DevArray):
            if (dtype is None or np.dtype(dtype) == a.dtype) and order in (None, "K", "A"):
                return a
            return np.asarray(a.view(np.ndarray), dtype=dtype, order=order).view(DevArray)
        if _is_host(a):
            self.record("h2d", a)
        return np.array(a, dtype=dtype, order=order, copy=True).view(DevArray)

    def array(self, a, dtype=None, copy=True, order="K", **kw):
        if isinstance(a, DevArray):
            return np.array(a.view(np.ndarray), dtype=dtype, copy=copy,
                            order=order).view(DevArray)
        if _is_host(a):
            self.record("h2d", a)
        return np.array(a, dtype=dtype, copy=True, order=order).view(DevArray)

    def ascontiguousarray(self, a, dtype=None):
        if isinstance(a, DevArray):
            return np.ascontiguousarray(a.view(np.ndarray), dtype=dtype).view(DevArray)
        if _is_host(a):
            self.record("h2d", a)
        return np.array(a, dtype=dtype, order="C", copy=True).view(DevArray)

    def asnumpy(self, a):
        if isinstance(a, DevArray):
            self.record("d2h", a)
            return np.array(a.view(np.ndarray), copy=True)
        return np.asarray(a)

    def zeros(self, shape, dtype=float, order="C"):
        nbytes = int(np.prod(shape)) * np.dtype(dtype).itemsize
        self.zeros_calls.append(nbytes)
        if self.fail_zeros and nbytes >= self.fail_zeros_min:
            self.fail_zeros -= 1
            raise _OutOfMemoryError("fake cupy: out of memory allocating {} bytes".format(nbytes))
        return np.zeros(shape, dtype=dtype, order=order).view(DevArray)


#: the module the current patch installed (DevArray.__setitem__ records into it)
_MODULE = None


@contextlib.contextmanager
def fake_cupy(free=2 ** 40):
    """Patch the backend so that ``gpu = true`` selects :class:`FakeCupy`;
    yields the module."""
    global _MODULE
    from hwave.solver import backend
    mod = FakeCupy(free=free)
    prev, _MODULE = _MODULE, mod
    try:
        with mock.patch.object(backend, "_import_cupy", lambda: mod):
            yield mod
    finally:
        _MODULE = prev


def is_device(x):
    return isinstance(x, DevArray)
