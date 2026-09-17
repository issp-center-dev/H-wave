"""backend helpers added for the bond-gate GPU path: to_device, device_available_bytes,
gpu_available, _oom_error_types, device_pool_used_bytes. All tests run without a GPU
(cupy is mocked)."""
import types
import unittest
from unittest import mock

import numpy as np

from hwave.solver import backend


class _Pool:
    def __init__(self, total, used):
        self._t, self._u = total, used
    def total_bytes(self):
        return self._t
    def used_bytes(self):
        return self._u


def _fake_cupy(free=10 * 2**30, total=40 * 2**30, pool_total=6 * 2**30, pool_used=2 * 2**30):
    cupy = types.SimpleNamespace()
    cupy.cuda = types.SimpleNamespace()
    cupy.cuda.runtime = types.SimpleNamespace(getDeviceCount=lambda: 1,
                                              memGetInfo=lambda: (free, total))
    cupy.get_default_memory_pool = lambda: _Pool(pool_total, pool_used)
    cupy.asarray = np.asarray
    return cupy


class TestToDevice(unittest.TestCase):
    def test_numpy_backend_is_identity(self):
        a = np.arange(3.0)
        self.assertIs(backend.to_device(a, np), a)

    def test_cupy_backend_calls_asarray(self):
        cupy = _fake_cupy()
        calls = []
        cupy.asarray = lambda a: calls.append(a) or a
        a = np.arange(3.0)
        backend.to_device(a, cupy)
        self.assertEqual(len(calls), 1)


class TestDeviceAvailableBytes(unittest.TestCase):
    def test_none_without_cupy(self):
        def _no():
            raise ImportError("no cupy")
        with mock.patch.object(backend, "_import_cupy", _no):
            self.assertIsNone(backend.device_available_bytes())

    def test_driver_free_plus_pool_cached(self):
        cupy = _fake_cupy(free=10 * 2**30, pool_total=6 * 2**30, pool_used=2 * 2**30)
        with mock.patch.object(backend, "_import_cupy", lambda: cupy):
            self.assertEqual(backend.device_available_bytes(), 14 * 2**30)

    def test_none_when_query_fails(self):
        cupy = _fake_cupy()
        def _boom():
            raise RuntimeError("cudaErrorNoDevice")
        cupy.cuda.runtime.memGetInfo = _boom
        with mock.patch.object(backend, "_import_cupy", lambda: cupy):
            self.assertIsNone(backend.device_available_bytes())


class TestGpuAvailable(unittest.TestCase):
    def test_false_without_cupy(self):
        def _no():
            raise ImportError("no cupy")
        with mock.patch.object(backend, "_import_cupy", _no):
            self.assertFalse(backend.gpu_available())

    def test_true_with_device(self):
        with mock.patch.object(backend, "_import_cupy", lambda: _fake_cupy()):
            self.assertTrue(backend.gpu_available())


class TestOomHelpers(unittest.TestCase):
    """The two helpers the solver's out-of-memory handler needs: a tuple that
    is safe to use in an ``except`` clause with or without cupy, and a pool
    reading that never raises."""

    def test_error_types_empty_without_cupy(self):
        def _no():
            raise ImportError("no cupy")
        with mock.patch.object(backend, "_import_cupy", _no):
            self.assertEqual(backend._oom_error_types(), ())

    def test_error_types_carry_the_cupy_exception(self):
        class _Oom(Exception):
            pass
        cupy = _fake_cupy()
        cupy.cuda.memory = types.SimpleNamespace(OutOfMemoryError=_Oom)
        with mock.patch.object(backend, "_import_cupy", lambda: cupy):
            self.assertEqual(backend._oom_error_types(), (_Oom,))
            # usable as an except clause
            try:
                raise _Oom("out of memory")
            except backend._oom_error_types():
                caught = True
            self.assertTrue(caught)

    def test_pool_used_bytes(self):
        cupy = _fake_cupy(pool_total=6 * 2**30, pool_used=2 * 2**30)
        with mock.patch.object(backend, "_import_cupy", lambda: cupy):
            self.assertEqual(backend.device_pool_used_bytes(), 2 * 2**30)

    def test_pool_used_bytes_zero_without_cupy(self):
        def _no():
            raise ImportError("no cupy")
        with mock.patch.object(backend, "_import_cupy", _no):
            self.assertEqual(backend.device_pool_used_bytes(), 0)


if __name__ == "__main__":
    unittest.main()
