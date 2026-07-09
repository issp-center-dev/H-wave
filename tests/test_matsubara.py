import numpy as np
import pytest
from hwave.solver import matsubara as ms


def test_fermion_tau_roundtrip():
    rng = np.random.default_rng(0)
    x = rng.standard_normal((3, 8, 5)) + 1j*rng.standard_normal((3, 8, 5))
    back = ms.tau_to_fermion(ms.fermion_to_tau(x, axis=1), axis=1)
    assert np.allclose(back, x, atol=1e-12)


def test_boson_tau_roundtrip():
    rng = np.random.default_rng(1)
    x = rng.standard_normal((8, 4)) + 1j*rng.standard_normal((8, 4))
    back = ms.tau_to_boson(ms.boson_to_tau(x, axis=0), axis=0)
    assert np.allclose(back, x, atol=1e-12)


def test_fermion_to_tau_matches_inline_formula():
    rng = np.random.default_rng(2)
    x = rng.standard_normal((2, 6)) + 1j*rng.standard_normal((2, 6))
    from hwave.solver.perf import FFT
    n = 6
    omg = np.exp(-1j*np.pi*(1.0/n - 1.0)*np.arange(n))
    expected = FFT.fft(x, axis=1) * omg[np.newaxis, :]
    assert np.allclose(ms.fermion_to_tau(x, axis=1), expected, atol=1e-14)
