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


def _reference(name, x, axis):
    """The transforms as written before #309 (layout aside)."""
    n = x.shape[axis]
    shape = [1] * x.ndim
    shape[axis] = n
    m = np.arange(n)
    phase = {"fermion_to_tau": np.exp(-1j * np.pi * (1.0 / n - 1.0) * m),
             "boson_to_tau": np.exp(-1j * np.pi * m),
             "tau_to_fermion": np.exp(1j * np.pi * (1.0 / n - 1.0) * m),
             "tau_to_boson": np.exp(1j * np.pi * (-1) * m)}[name].reshape(shape)
    if name.endswith("_to_tau"):
        return np.fft.fft(x, axis=axis) * phase
    return np.fft.ifft(x * phase, axis=axis)


@pytest.mark.parametrize("name", ["fermion_to_tau", "boson_to_tau",
                                  "tau_to_fermion", "tau_to_boson"])
@pytest.mark.parametrize("axis", [0, 1, 2, -1])
def test_transforms_return_c_ordered_values_unchanged(name, axis):
    # numpy < 2 leaves the transformed axis with unit stride; the transforms
    # return a C-ordered array (#309) with the same values, bit for bit
    rng = np.random.default_rng(3)
    x = rng.standard_normal((6, 8, 5)) + 1j * rng.standard_normal((6, 8, 5))
    out = getattr(ms, name)(x, axis=axis)
    assert out.flags.c_contiguous
    assert out.dtype == np.complex128
    assert np.array_equal(out, _reference(name, x, axis))


def test_transforms_accept_real_complex64_and_non_contiguous_input():
    rng = np.random.default_rng(4)
    z = rng.standard_normal((6, 8, 4)) + 1j * rng.standard_normal((6, 8, 4))
    inputs = {"real": rng.standard_normal((8, 6, 4)),
              "complex64": z.astype(np.complex64),
              "F-ordered": np.asfortranarray(z),
              "transposed": z.transpose(2, 0, 1),
              "sliced": z[:, ::2, :]}
    for name in ("fermion_to_tau", "boson_to_tau", "tau_to_fermion", "tau_to_boson"):
        for label, arr in inputs.items():
            for axis in (0, 1, 2):
                out = getattr(ms, name)(arr, axis=axis)
                ref = _reference(name, arr, axis)
                assert out.flags.c_contiguous, (name, label, axis)
                assert out.dtype == ref.dtype, (name, label, axis)
                assert np.array_equal(out, ref), (name, label, axis)
