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


# ---------------------------------------------------------------------------
# When the mirrored slices need no conditioning guard of their own
# ---------------------------------------------------------------------------
#
# The guard scores the RPA denominator D(q, l) = 1 - s chibar(q, l) V(q)
# (s = -1 spin, +1 charge; V = S or C) by sigma_min/sigma_max and
# sigma_min/max(1, sigma_max). The half grid solves D only at the stored
# slices; at a mirrored slice l' = nmat - l it would be
#
#     D(q, l') = 1 - s chibar(q, l)^dagger V(q),                        (2)
#
# by (1). Its determinant is that of D(q, l)^dagger (push-through, V
# Hermitian), but its singular values are those of 1 - s V(q) chibar(q, l),
# NOT of D(q, l): the two products differ in order, so their scores can differ
# by up to cond(V) (measured: 0.3 % with complex hopping near an instability).
#
# If in addition, at every stored slice,
#
#     chibar(q, l)^T = chibar(-q, l)    and    V(q)^T = V(-q),          (3)
#
# then chibar(q, l)^dagger = conj(chibar(-q, l)) and V(q) = conj(V(-q)) (V
# Hermitian), so (2) is
#
#     D(q, l') = conj(1 - s chibar(-q, l) V(-q)) = conj(D(-q, l)):
#
# entrywise conjugation leaves the singular values unchanged, so the
# mirrored block at (q, l') has EXACTLY the score of the stored block at
# (-q, l). The stored slices then already contain every score of the
# mirrored ones, and the minimum over the whole grid is the minimum over the
# stored slices.
#
# (3) is time-reversal symmetry of the spin-free problem: for a Hamiltonian
# with real hopping, real interaction coefficients and a real seed,
# H(-k) = H(k)^T and Sigma(-k, i w) = Sigma(k, i w)^T, hence
# G(-k, i w) = G(k, i w)^T; inserted into the bubble
# chibar_{AB}(q) = -(T/N) sum_k e^{i(k-q)(R_A - R_B)} G(k) G(k - q) and
# substituting k -> q - k gives chibar(-q) = chibar(q)^T on the pair matrix,
# and the bond vertices built from real coefficients satisfy
# V(-q) = V(q)^T. A complex (time-reversal breaking) hopping breaks (3).
# Rather than infer (3) from the inputs, :func:`mirror_scores_are_stored`
# checks it directly on the stored chibar of the current map and on the
# vertices: a deviation delta of (3) moves a singular value of D by at most
# ||delta V|| (Weyl), so a tolerance at round-off level keeps the scores equal
# to round-off.

#: relative tolerance of (3): max |X(q)^T - X(-q)| / max |X|
TIME_REVERSAL_RTOL = 1.0e-12


def minus_q_index(spatial_shape):
    """The C-order flattened index of -q for every flattened q."""
    nx, ny, nz = (int(x) for x in spatial_shape)
    qi = np.arange(nx * ny * nz).reshape(nx, ny, nz)
    return qi[(-np.arange(nx)) % nx][:, (-np.arange(ny)) % ny][:, :, (-np.arange(nz)) % nz].ravel()


def _transpose_dev(x, mq):
    """max |x(q)^T - x(-q)| over the leading axes (q on axis -3)."""
    xp = _bk.array_module_of(x)
    return float(xp.abs(x.swapaxes(-1, -2) - x[..., mq, :, :]).max()) if x.size else 0.0


def mirror_scores_are_stored(chibar_half, S, C, spatial_shape, rtol=TIME_REVERSAL_RTOL):
    """``(holds, deviation)``: whether (3) holds for the stored slices
    ``chibar_half`` ``(nslices, nvol, ND, ND)`` and the vertices ``S``,
    ``C`` ``(nvol, ND, ND)``, to ``rtol`` relative to the largest entry of
    each. ``deviation`` is the largest relative deviation found (the scan
    stops at the first slice that fails)."""
    mq = minus_q_index(spatial_shape)
    worst = 0.0
    for V in (S, C):
        scale = max(1.0, float(np.abs(_bk.to_host(V)).max()) if np.size(V) else 1.0)
        worst = max(worst, _transpose_dev(V, mq) / scale)
        if worst > rtol:
            return False, worst
    xp = _bk.array_module_of(chibar_half)
    scale = float(xp.abs(chibar_half).max()) if chibar_half.size else 0.0
    if scale == 0.0:
        return True, worst
    for h in range(int(chibar_half.shape[0])):
        worst = max(worst, _transpose_dev(chibar_half[h], mq) / scale)
        if worst > rtol:
            return False, worst
    return True, worst


def seed_mirror_deviation(sigma_fluct):
    """Relative deviation of ``Sigma(k, -i w) = Sigma(k, i w)^dagger`` of a
    fermionic ``(nb, nmat, nvol, n, n)`` self-energy: the fermionic slot
    ``n`` (``w = pi T (2 n + 1 - nmat)``) pairs with ``nmat - 1 - n``."""
    x = np.asarray(sigma_fluct)
    if x.size == 0:
        return 0.0
    y = np.conj(x[:, ::-1]).swapaxes(-1, -2)
    return float(np.abs(x - y).max()) / max(1.0, float(np.abs(x).max()))
