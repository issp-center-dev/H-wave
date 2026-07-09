"""Fermionic/bosonic Matsubara <-> imaginary-time transforms.

Single source of truth for the phase/FFT conventions used by FLEX (chi0 = G*G,
Sigma = V*G) and the dynamic Eliashberg kernel. Verbatim code-motion of the
inline transforms in rpa._calc_chi0q / flex._calc_self_energy: the phase factors
are copied, not re-derived. The temperature factor 1/beta is NOT applied here
(callers apply it), so these are pure transforms.
"""
import numpy as np

from hwave.solver.perf import FFT


def _bcast(vec, ndim, axis):
    shape = [1] * ndim
    shape[axis] = vec.size
    return vec.reshape(shape)


def fermion_to_tau(arr, axis):
    n = arr.shape[axis]
    omg = np.exp(-1j * np.pi * (1.0 / n - 1.0) * np.arange(n))
    return FFT.fft(arr, axis=axis) * _bcast(omg, arr.ndim, axis)


def boson_to_tau(arr, axis):
    n = arr.shape[axis]
    omg = np.exp(-1j * np.pi * np.arange(n))
    return FFT.fft(arr, axis=axis) * _bcast(omg, arr.ndim, axis)


def tau_to_fermion(arr, axis):
    n = arr.shape[axis]
    omg_inv = np.exp(1j * np.pi * (1.0 / n - 1.0) * np.arange(n))
    return FFT.ifft(arr * _bcast(omg_inv, arr.ndim, axis), axis=axis)


def tau_to_boson(arr, axis):
    n = arr.shape[axis]
    omg_inv = np.exp(1j * np.pi * (-1) * np.arange(n))
    return FFT.ifft(arr * _bcast(omg_inv, arr.ndim, axis), axis=axis)
