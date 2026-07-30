"""The FFT-grid reversal map, single-sourced (#108).

``reverse_fft_axes`` applies the index map ``i -> (N - i) % N`` along the
given axes: the momentum/position negation ``k -> -k`` (``r -> -r``) on an
FFT-ordered grid, where index 0 is the Gamma point / origin and stays put.

Before this module the map was spelled three different ways at eighteen
call sites -- ``roll(arr[::-1], +1)``, ``flip(roll(arr, -1))`` and a
fancy-index gather ``arr[(-arange(n)) % n]`` -- kept in sync only by
comments ("same reverse+roll map as the uniform path"). The spellings are
the same gather,
so consolidating them is bit-identical; having ONE implementation matters
because the map has a known trap: for an axis of length N <= 2 it is the
identity, so a misapplied reversal (wrong axes, or reversing an orbital
axis along with the lattice ones) is invisible on one- and two-site
fixtures and was exactly the defect history of the transverse channel
(see rpa.py's general-path note). The tests pin the map's algebra --
involution, index-0 fixed point, equivalence of all three historical
spellings -- once, here.

The helper is a GENERIC periodic-grid index reversal: it applies to any
axis whose index 0 is the origin of a periodic FFT-ordered coordinate.
That includes the uniform imaginary-time grid (tau_j = j*beta/N, so
-tau_j wraps to tau_{(N-j) % N}): rpa.py's chi0 kernels deliberately pass
their tau axis through this map, with the fermionic sign handled
separately. What stays OUT of this module is every frequency/time axis
whose reversal is NOT this map -- a centered fermionic Matsubara axis
reverses as ``n -> nmat - 1 - n`` with no roll, a symmetric IR node set
as a plain flip, and flex.py's tau-node sets likewise flip without
rolling. Callers compose those locally where the physics lives.
"""

import numpy as np

from hwave.solver import backend as _bk


def reverse_fft_axes(arr, axes):
    """Apply ``i -> (N - i) % N`` along ``axes`` of an FFT-ordered array.

    Works on numpy and cupy arrays alike (the module is taken from the
    array itself). The map is an involution and leaves index 0 in place.

    Parameters
    ----------
    arr : ndarray
        Array carrying FFT-ordered axes (index 0 = Gamma / origin).
    axes : tuple of int
        The axes to negate. Axes not listed are untouched.

    Returns
    -------
    ndarray
        New array with the listed axes reversed on the FFT grid.
    """
    xp = _bk.array_module_of(arr)
    axes = tuple(axes)
    return xp.roll(xp.flip(arr, axis=axes), 1, axis=axes)
