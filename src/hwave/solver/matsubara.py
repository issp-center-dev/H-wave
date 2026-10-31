"""Fermionic/bosonic Matsubara <-> imaginary-time transforms.

Single source of truth for the phase/FFT conventions used by FLEX (chi0 = G*G,
Sigma = V*G) and the dynamic Eliashberg kernel. Verbatim code-motion of the
inline transforms in rpa._calc_chi0q / flex._calc_self_energy: the phase factors
are copied, not re-derived. The temperature factor 1/beta is NOT applied here
(callers apply it), so these are pure transforms.

The transforms are array-backend agnostic: the FFT and the phase factors are
taken from the module (numpy or cupy) that owns the input array, so device
arrays stay on the device.
"""
import numpy as np

from hwave.solver.backend import array_module_of


def _bcast(vec, ndim, axis):
    shape = [1] * ndim
    shape[axis] = vec.size
    return vec.reshape(shape)


# Every transform returns a C-ordered array. numpy < 2 transforms an inner
# axis by swapping it to the end and back, so its result has the transformed
# (frequency / time) axis with unit stride; the stacked matmul / solve calls
# downstream are several times slower with OpenBLAS on that layout (#309).
# The C-ordered result is written into a buffer the transform allocates
# anyway, so the peak memory does not grow.


def _forward(xp, arr, axis, phase):
    f = xp.fft.fft(arr, axis=axis)
    out = xp.empty(f.shape, dtype=xp.result_type(f, phase))
    xp.multiply(f, _bcast(phase, arr.ndim, axis), out=out)
    return out


def _backward(xp, arr, axis, phase):
    tmp = xp.empty(arr.shape, dtype=xp.result_type(arr, phase))
    xp.multiply(arr, _bcast(phase, arr.ndim, axis), out=tmp)
    res = xp.fft.ifft(tmp, axis=axis)
    if res.flags.c_contiguous:
        return res
    if res.dtype == tmp.dtype:
        tmp[...] = res
        return tmp
    del tmp
    return xp.ascontiguousarray(res)


def fermion_to_tau(arr, axis):
    xp = array_module_of(arr)
    n = arr.shape[axis]
    omg = xp.exp(-1j * np.pi * (1.0 / n - 1.0) * xp.arange(n))
    return _forward(xp, arr, axis, omg)


def boson_to_tau(arr, axis):
    xp = array_module_of(arr)
    n = arr.shape[axis]
    omg = xp.exp(-1j * np.pi * xp.arange(n))
    return _forward(xp, arr, axis, omg)


def tau_to_fermion(arr, axis):
    xp = array_module_of(arr)
    n = arr.shape[axis]
    omg_inv = xp.exp(1j * np.pi * (1.0 / n - 1.0) * xp.arange(n))
    return _backward(xp, arr, axis, omg_inv)


def tau_to_boson(arr, axis):
    xp = array_module_of(arr)
    n = arr.shape[axis]
    omg_inv = xp.exp(1j * np.pi * (-1) * xp.arange(n))
    return _backward(xp, arr, axis, omg_inv)


def tau_to_boson_row_weights(n, index, xp=np):
    """Weights ``w`` (length ``n``) of ONE output row of :func:`tau_to_boson`:
    ``tau_to_boson(arr, axis)[index] == sum_m w[m] * arr[m]`` along ``axis``
    (the ``1/n`` of the inverse FFT included). Lets a caller that needs a
    single bosonic frequency (the static ``nmat // 2`` row) accumulate the
    imaginary-time sum without materializing the whole frequency axis.

    The phase is :func:`tau_to_boson`'s own ``exp(-i pi m)``, evaluated by
    the same expression so the two agree to round-off; the inverse-DFT kernel
    ``exp(2 pi i index m / n)`` is taken EXACTLY as ``(-1)**m`` at the
    centre row ``index = n // 2`` of an even ``n`` (the zero bosonic
    frequency), and evaluated otherwise."""
    m = xp.arange(n)
    phase = xp.exp(1j * np.pi * (-1) * m)
    if n % 2 == 0 and index == n // 2:
        kernel = 1 - 2 * (m % 2)
    else:
        kernel = xp.exp(2j * np.pi * index * m / n)
    return phase * kernel / n
