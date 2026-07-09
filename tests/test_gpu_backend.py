"""Tests for hwave.solver.backend (numpy/cupy array-backend selection)."""
import logging

import numpy as np
import pytest


def test_get_backend_cpu_returns_numpy():
    from hwave.solver import backend
    xp, gpu_active = backend.get_backend(False)
    assert xp is np
    assert gpu_active is False


def test_get_backend_gpu_falls_back_without_cupy(monkeypatch, caplog):
    """gpu=True on a machine without (usable) CuPy must warn and fall back to
    numpy, not crash: the run is numerically identical, only slower."""
    from hwave.solver import backend

    def _no_cupy():
        raise ImportError("No module named 'cupy'")

    monkeypatch.setattr(backend, "_import_cupy", _no_cupy)
    logger = logging.getLogger("test_backend_fallback")
    with caplog.at_level(logging.WARNING, logger=logger.name):
        xp, gpu_active = backend.get_backend(True, logger=logger)
    assert xp is np
    assert gpu_active is False
    assert any("cupy" in rec.message.lower() for rec in caplog.records)


def test_get_backend_gpu_active_with_cupy():
    cupy = pytest.importorskip("cupy")
    try:
        ndev = cupy.cuda.runtime.getDeviceCount()
    except Exception:
        ndev = 0
    if ndev == 0:
        pytest.skip("cupy installed but no CUDA device")
    from hwave.solver import backend
    xp, gpu_active = backend.get_backend(True)
    assert xp is cupy
    assert gpu_active is True


def test_array_module_of_numpy():
    from hwave.solver import backend
    assert backend.array_module_of(np.zeros(3)) is np


def test_to_host_numpy_passthrough():
    from hwave.solver import backend
    a = np.arange(4.0)
    assert backend.to_host(a) is a


def test_array_module_and_to_host_cupy():
    cupy = pytest.importorskip("cupy")
    try:
        a = cupy.arange(4.0)
    except Exception:
        pytest.skip("cupy installed but no usable CUDA device")
    from hwave.solver import backend
    assert backend.array_module_of(a) is cupy
    back = backend.to_host(a)
    assert isinstance(back, np.ndarray)
    np.testing.assert_array_equal(back, np.arange(4.0))


def test_matsubara_transforms_cupy_match_numpy():
    """The matsubara transforms must give identical results (up to fp64
    round-off) on cupy arrays, since the dynamic Eliashberg GPU kernel calls
    them on device arrays."""
    cupy = pytest.importorskip("cupy")
    try:
        cupy.zeros(1)
    except Exception:
        pytest.skip("cupy installed but no usable CUDA device")
    from hwave.solver import matsubara as ms
    rng = np.random.default_rng(2)
    arr = rng.standard_normal((3, 8, 5)) + 1j * rng.standard_normal((3, 8, 5))
    for fn in (ms.fermion_to_tau, ms.boson_to_tau,
               ms.tau_to_fermion, ms.tau_to_boson):
        ref = fn(arr, axis=1)
        got = cupy.asnumpy(fn(cupy.asarray(arr), axis=1))
        assert isinstance(fn(cupy.asarray(arr), axis=1), cupy.ndarray)
        np.testing.assert_allclose(got, ref, atol=1e-12)
