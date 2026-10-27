"""The non-negative half of the uniform bosonic Matsubara grid (issue #223).

Every bosonic object of the bond-resolved FLEX map -- the bubble ``chibar``,
the dressed spin and charge channels, the effective interaction ``W`` and
the channel-0 collapses -- satisfies, at the SAME wave vector ``q``,

    X(q, l) = X(q, nmat - l)^dagger                                  (1)

where ``^dagger`` is the conjugate transpose of the pair matrix (the
``ND x ND`` block, or ``nd x nd`` for a collapse) and ``l`` is the index of
the uniform grid ``nu_l = pi T (2 l - nmat)`` (``l = nmat // 2`` is the
static slice). It follows from ``G(k, -i w) = G(k, i w)^dagger`` and from
vertices that are Hermitian at every ``q``; time reversal is not used.

The half grid keeps ``nmat // 2 + 1`` slices: slice ``h`` holds the grid
index

    l(h) = (h + nmat // 2) mod nmat,    h = 0 ... nmat // 2,

so ``h = 0`` is the static slice and ``h = nmat // 2`` is the Nyquist slot
``l = 0``. Both endpoints are their own partners under (1) and must be
stored. Every other ``l`` in ``(0, nmat / 2)`` is rebuilt from its partner
``nmat - l``, which sits at slice ``h = nmat / 2 - l``.

For a channel-pair block ``(alpha, beta)`` of an ``ND x ND`` object, (1)
reads ``X_{alpha beta}(l) = X_{beta alpha}(nmat - l)^dagger``: the mirror of
a pair block reads the TRANSPOSED pair ``(beta, alpha)``.

This module is the single home of these index maps and of the mirror; the
helpers are array-module generic (numpy or cupy, by the module of the
input).
"""
import numpy as np

from . import backend as _bk


def nslices(nmat):
    """Number of stored slices of the half grid for an even ``nmat``."""
    nmat = int(nmat)
    if nmat < 2 or nmat % 2:
        raise ValueError("half_grid: nmat must be an even integer >= 2, got {}".format(nmat))
    return nmat // 2 + 1


def grid_index(nmat, s0=0, s1=None):
    """``l(h)`` for the slices ``h = s0 ... s1 - 1`` (int64 numpy array)."""
    ns = nslices(nmat)
    s1 = ns if s1 is None else int(s1)
    if not 0 <= int(s0) <= s1 <= ns:
        raise ValueError("half_grid: slice range [{}, {}) outside [0, {}]".format(s0, s1, ns))
    return (np.arange(int(s0), s1, dtype=np.int64) + int(nmat) // 2) % int(nmat)


def half_from_full(full, nmat):
    """The stored slices ``full[l(h)]``, ``h = 0 ... nmat // 2``, of a
    full-grid array (frequency on axis 0): ``l = nmat/2 ... nmat - 1`` then
    the Nyquist slot ``l = 0``. A new array of the module of ``full``."""
    xp = _bk.array_module_of(full)
    h0 = int(nmat) // 2
    if int(full.shape[0]) != int(nmat):
        raise ValueError("half_grid.half_from_full: axis 0 has {} slices, not nmat = {}"
                         .format(full.shape[0], nmat))
    return xp.concatenate((full[h0:], full[:1]), axis=0)


def _dagger(x):
    xp = _bk.array_module_of(x)
    return xp.conj(x).swapaxes(-1, -2)


def expand(ab, nmat, ba=None, out=None, chunk=None):
    """The full grid ``(nmat, ...)`` of a half-grid array ``ab`` (frequency
    on axis 0, ``nmat // 2 + 1`` slices; the last two axes are the matrix
    the dagger of (1) acts on).

    ``ba`` is the half-grid array the mirrored slices are read from: the
    TRANSPOSED channel pair ``(beta, alpha)`` when ``ab`` is the pair block
    ``(alpha, beta)``; ``None`` (the default) means ``ab`` itself, which is
    the case of a whole ``ND x ND`` object or of a channel-0 collapse.

    ``out`` (optional) receives the result, e.g. a memory map or a
    preallocated buffer; otherwise a new array of the module of ``ab`` is
    returned. ``chunk`` bounds the number of mirrored slices conjugated at
    once (``None``: all of them), so the temporary of the mirror is at most
    ``chunk`` slices."""
    nmat = int(nmat)
    h0 = nmat // 2
    ns = nslices(nmat)
    ba = ab if ba is None else ba
    if int(ab.shape[0]) != ns or int(ba.shape[0]) != ns:
        raise ValueError("half_grid.expand: the half-grid arrays must have nmat // 2 + 1 = {} "
                         "slices on axis 0, got {} and {}".format(ns, ab.shape[0], ba.shape[0]))
    xp = _bk.array_module_of(ab)
    if out is None:
        out = xp.empty((nmat,) + tuple(ab.shape[1:]), dtype=ab.dtype)
    out[h0:] = ab[:h0]            # l = nmat/2 ... nmat - 1   (h = 0 ... nmat/2 - 1)
    out[0] = ab[h0]               # Nyquist l = 0             (h = nmat/2)
    # l = 1 ... nmat/2 - 1: partner nmat - l at h = nmat/2 - l
    step = (h0 - 1) if chunk is None else max(1, int(chunk))
    for l_lo in range(1, h0, max(step, 1)):
        l_hi = min(h0, l_lo + step)
        # h = h0 - l for l = l_lo ... l_hi - 1, i.e. h descending from h0 - l_lo
        hs = slice(h0 - l_lo, h0 - l_hi, -1)
        out[l_lo:l_hi] = _dagger(ba[hs])
    return out


def fill_mirror(full, nmat):
    """In place: given a FULL-grid array whose stored slices ``l(h)`` are
    already set, write the mirrored slices ``l = 1 ... nmat/2 - 1`` as the
    conjugate transposes of ``full[nmat - l]``. Returns ``full``."""
    nmat = int(nmat)
    h0 = nmat // 2
    if int(full.shape[0]) != nmat:
        raise ValueError("half_grid.fill_mirror: axis 0 has {} slices, not nmat = {}"
                         .format(full.shape[0], nmat))
    if h0 > 1:
        # l = 1 ... h0 - 1 <- nmat - l = nmat - 1 ... h0 + 1
        full[1:h0] = _dagger(full[nmat - 1:h0:-1])
    return full


def describe_slices(nmat, s0, s1):
    """Message text naming the half-grid slices ``[s0, s1)`` by their grid
    indices (``l`` is contiguous modulo ``nmat``)."""
    ls = grid_index(nmat, s0, s1)
    return "of half-grid slices [{}, {}) (grid indices l = {} ... {})".format(
        int(s0), int(s1), int(ls[0]), int(ls[-1]))
