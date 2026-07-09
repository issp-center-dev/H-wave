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


def fermion_to_tau(arr, axis):
    xp = array_module_of(arr)
    n = arr.shape[axis]
    omg = xp.exp(-1j * np.pi * (1.0 / n - 1.0) * xp.arange(n))
    return xp.fft.fft(arr, axis=axis) * _bcast(omg, arr.ndim, axis)


def boson_to_tau(arr, axis):
    xp = array_module_of(arr)
    n = arr.shape[axis]
    omg = xp.exp(-1j * np.pi * xp.arange(n))
    return xp.fft.fft(arr, axis=axis) * _bcast(omg, arr.ndim, axis)


def tau_to_fermion(arr, axis):
    xp = array_module_of(arr)
    n = arr.shape[axis]
    omg_inv = xp.exp(1j * np.pi * (1.0 / n - 1.0) * xp.arange(n))
    return xp.fft.ifft(arr * _bcast(omg_inv, arr.ndim, axis), axis=axis)


def tau_to_boson(arr, axis):
    xp = array_module_of(arr)
    n = arr.shape[axis]
    omg_inv = xp.exp(1j * np.pi * (-1) * xp.arange(n))
    return xp.fft.ifft(arr * _bcast(omg_inv, arr.ndim, axis), axis=axis)
