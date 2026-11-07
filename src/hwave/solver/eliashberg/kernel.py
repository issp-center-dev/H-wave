"""Green's function layouts, the pair bubble G2 and the static Eliashberg kernel of hwave_sc.

Moved verbatim from hwave.sc (dev/sc-split-plan.md, step 1); the shell
re-exports every name here under its old spelling.
"""
from __future__ import annotations

import logging
from typing import NamedTuple
import numpy as np
from numpy.fft import fftn, ifftn
from scipy.sparse.linalg import LinearOperator, eigs, bicgstab, gmres, lgmres

from hwave.solver.kgrid import reverse_fft_axes
from hwave.solver import backend
from hwave.solver.eliashberg.config import (  # noqa: F401
    _finite_positive_float64)

logger = logging.getLogger("hwave_sc")


# Numerically safe upper bound for beta in the pair bubble: the batched
# GEMM forms ~beta^2 intermediates before the division by beta, so a
# physically meaningless beta (T < 1e-75) would overflow them to NaN.
_G2_BETA_MAX = 1.0e75


def _green_sc_to_canonical(g):
    """Convert sc.py's Green-function layout to the shared kernel's
    canonical layout.

    ``(p, p, Nx, Ny, Nz, nmat) -> (1, nmat, Nx*Ny*Nz, p, p)``: transpose to
    ``(nmat, Nx, Ny, Nz, p, p)`` order, reshape the spatial axes to
    ``nvol`` in C order, and insert a leading length-1 block axis (the
    canonical layout's ``nblock`` -- sc.py's Green function never carries a
    block axis of its own). ``Nx, Ny, Nz`` are read directly off ``g``'s
    shape, so no separate ``spatial_shape`` argument is needed on this
    direction. ALWAYS returns a copy (``np.ascontiguousarray``) -- never a
    view aliasing ``g`` (unified-bubble-kernel spec, "Layout adapters").
    """
    g = np.asarray(g)
    if g.ndim != 6:
        raise ValueError(
            "_green_sc_to_canonical: expected ndim=6 (p, p, Nx, Ny, Nz, "
            "nmat), got ndim={} shape={}".format(g.ndim, g.shape))
    p, p2, Nx, Ny, Nz, nmat = g.shape
    if p2 != p:
        raise ValueError(
            "_green_sc_to_canonical: the two orbital axes must be equal; "
            "got shape {}".format(g.shape))
    transposed = g.transpose(5, 2, 3, 4, 0, 1)          # (nmat,Nx,Ny,Nz,p,p)
    canonical = transposed.reshape(nmat, Nx * Ny * Nz, p, p)[np.newaxis, ...]
    return np.ascontiguousarray(canonical)


def _green_canonical_to_sc(g, spatial_shape):
    """Exact inverse of :func:`_green_sc_to_canonical`.

    ``(1, nmat, nvol, p, p) -> (p, p, Nx, Ny, Nz, nmat)``. ``nvol`` alone
    does not determine its own ``(Nx, Ny, Nz)`` factorization, so this
    direction takes ``spatial_shape`` explicitly (the same way every other
    kernel entry point -- ``bond_bubble_static``, ``dense_bubble`` -- does).
    ALWAYS returns a copy (``np.ascontiguousarray``).
    """
    g = np.asarray(g)
    if g.ndim != 5 or g.shape[0] != 1:
        raise ValueError(
            "_green_canonical_to_sc: expected canonical shape (1, nmat, "
            "nvol, p, p); got ndim={} shape={}".format(g.ndim, g.shape))
    _, nmat, nvol, p, p2 = g.shape
    if p2 != p:
        raise ValueError(
            "_green_canonical_to_sc: the two orbital axes must be equal; "
            "got shape {}".format(g.shape))
    Nx, Ny, Nz = spatial_shape
    if Nx * Ny * Nz != nvol:
        raise ValueError(
            "_green_canonical_to_sc: spatial_shape {} does not match "
            "nvol={} (shape {})".format(spatial_shape, nvol, g.shape))
    reshaped = g[0].reshape(nmat, Nx, Ny, Nz, p, p)
    out = reshaped.transpose(4, 5, 1, 2, 3, 0)          # (p,p,Nx,Ny,Nz,nmat)
    return np.ascontiguousarray(out)


class _BondGreen(NamedTuple):
    """Carrier separating the FULL Green function -- fed, unchanged, to
    EVERY downstream consumer of the bond path (``pair_weight``, the pair
    bubble ``G2`` of the bond kernel) -- from the
    DEFLATED Green/tail pair, fed ONLY to the bubble
    (``bubble.bond_bubble_static``). See the unified-bubble-kernel spec's
    "Green data flow in sc.py's bond branch".

    Attributes
    ----------
    full_sc : ndarray
        sc.py layout ``(p, p, Nx, Ny, Nz, nmat)`` -- the bare Green
        function ``G(k, iw_n)`` built from the eigen-decomposition,
        independent of ``coeff_tail`` (the tail term the bubble subtracts
        in frequency space is added straight back for ``full``; only
        ``deflated_kw`` depends on the coefficient).
    deflated_kw : ndarray
        Canonical ``(1, nmat, nvol, p, p)``, for the bubble only.
    green0_tail : ndarray or None
        Canonical, paired with ``deflated_kw``; ``None`` when the tail
        correction is off (``coeff_tail_applied == 0``).
    coeff_tail_requested : float
        The resolved ``[eliashberg] bond_coeff_tail`` config value.
    coeff_tail_applied : float
        What the build actually used: equals ``coeff_tail_requested`` for
        an internal build, ``0.0`` for an external Green function (the
        tail correction is not well defined for a file whose asymptote is
        unknown).
    source : str
        ``"internal"`` (built from sc.py's own eigen-decomposition) or
        ``"external_npz"`` (``[eliashberg] bond_green``).
    """
    full_sc: np.ndarray
    deflated_kw: np.ndarray
    green0_tail: np.ndarray | None
    coeff_tail_requested: float
    coeff_tail_applied: float
    source: str


def _calc_g2(green_kw, beta, tail=True):
    """Calculate G2 = T * sum_n G(k, wn) G(-k+q, -wn).

    The temperature factor T = 1/beta is included so that the
    Eliashberg kernel correctly computes:
        lambda * sigma(k) = -(T/N_L) sum_{k',n'} P(k-k') G(k') G(-k') sigma(k')

    With ``tail=True`` (the default) the truncated Matsubara sum gets the
    analytic high-frequency tail correction (issue #86). The summand's exact
    leading tail is delta_ij delta_lm / wn^2 -- the 1/(i wn) coefficient of
    G is the identity by completeness of the eigenbasis, and a self-energy
    only enters at the next order, so this holds for dressed (FLEX) Green
    functions too, with no free coefficient (unlike chi0q's user-tunable
    ``coeff_tail``). Subtracting that model inside the window and adding its
    exact full sum (1/beta) sum_{n in Z} 1/wn^2 = beta/4 amounts to adding

        c = beta/4 - (1/beta) sum_{n in window} 1/wn^2   (> 0)

    to G2[i,i,l,l], i.e. c times the identity on the (il) gap space. This
    reduces the truncation error from O(1/Nmat) to the next order and in
    practice restores the positive semi-definiteness of G2 that makes the
    static Eliashberg kernel similar to a Hermitian matrix (real spectrum);
    the higher-order remainder can still leave a small exact eigenvalue
    negative, which _warn_if_g2_indefinite reports. The bare truncated sum
    can be slightly indefinite, which then injects spurious imaginary parts
    into the reported eigenvalues at small Nmat.

    Parameters
    ----------
    green_kw : ndarray
        Green's function, shape (norb, norb, Nx, Ny, Nz, nmat).
    beta : float
        Inverse temperature.

    Returns
    -------
    G2 : ndarray
        Shape (norb, norb, norb, norb, Nx, Ny, Nz).
    """
    norb = green_kw.shape[0]
    Nx, Ny, Nz, nmat = (green_kw.shape[2], green_kw.shape[3],
                        green_kw.shape[4], green_kw.shape[5])
    nvol = Nx * Ny * Nz

    # beta gate (round-6 review): a non-finite or non-positive beta poisons
    # every frequency below, and a finite but astronomically large one
    # overflows the batched GEMM's ~beta^2 intermediates to NaN before the
    # division by beta -- measured: beta = 1e200 turns the whole pair
    # bubble into NaN while G itself is finite. beta values beyond
    # _G2_BETA_MAX have no physical meaning (T < 1e-75), so reject rather
    # than restructure the arithmetic (the tail=False stream must stay
    # bit-identical to the pre-correction implementation for valid input).
    beta_f = _finite_positive_float64(beta)
    if beta_f is None:
        raise ValueError(
            "beta must be a single finite positive real scalar, got "
            "{!r}.".format(beta))
    beta = beta_f
    if beta > _G2_BETA_MAX:
        raise ValueError(
            "beta = {!r} exceeds the numerically safe range (<= {!r}): "
            "the pair-bubble accumulation would overflow to NaN. Check "
            "the temperature.".format(beta, _G2_BETA_MAX))
    # Grid gate BEFORE any O(norb^2 nvol nmat) work (round-6 review: the
    # guard previously ran after the GEMM, which emitted overflow warnings
    # from data the call then rejected). The centered grid
    # wn = (2n + 1 - nmat) * pi / beta is only fermionic for EVEN nmat: an
    # odd nmat puts wn = 0 in the window, and nmat <= 0 would return the
    # bare analytic shift beta/4 with no Green-function samples at all.
    if tail and (nmat <= 0 or nmat % 2 != 0):
        raise ValueError(
            "the Matsubara tail correction (g2_tail) requires an even, "
            "positive number of frequencies on the centered fermionic "
            "grid; the Green function has nmat = {}. Fix the frequency "
            "axis of the input Green function (or Nmat), or set "
            "[eliashberg] g2_tail = false.".format(nmat))
    # Promote reduced-precision (e.g. complex64 file) input once so the
    # accumulation runs in complex128; a no-op for complex128 input, which
    # keeps the tail=False bit-parity with the pre-correction code.
    if green_kw.dtype != np.complex128:
        green_kw = np.ascontiguousarray(green_kw, dtype=np.complex128)

    # G(-k, -wn): centered-Matsubara flip on the frequency axis, then the
    # shared FFT-grid map k -> -k on the spatial axes.
    green_kw_inv = reverse_fft_axes(green_kw[..., ::-1], (2, 3, 4))
    # einsum("ijpqsk, lmpqsk -> ijlmpqs") sums over k (nmat dimension)
    # Reshape to use tensordot for BLAS: contract last axis (nmat) after
    # merging spatial dims
    # A: (norb, norb, nvol, nmat) -> (norb^2, nvol*nmat) -- but need per-site sum
    # Better: reshape to (norb*norb, nvol, nmat) then per-site outer product
    A = green_kw.reshape(norb * norb, nvol, nmat)       # (ij, site, k)
    B = green_kw_inv.reshape(norb * norb, nvol, nmat)   # (lm, site, k)
    # G2[ij, lm, site] = sum_k A[ij, site, k] * B[lm, site, k]
    # Per-site batched GEMM over the Matsubara axis n. numpy.einsum does NOT
    # lower 'isn,jsn->ijs' to BLAS GEMM; move site to the batch axis and matmul:
    #   As[s] @ Bs[s].T -> [i, j] per site; moveaxis(...,0,2) -> (i, j, s).
    As = np.moveaxis(A, 1, 0)             # (nvol, norb^2, nmat)
    Bs = np.moveaxis(B, 1, 0)             # (nvol, norb^2, nmat)
    G2 = np.moveaxis(As @ Bs.transpose(0, 2, 1), 0, 2)  # (norb^2, norb^2, nvol)
    G2 = G2.reshape(norb, norb, norb, norb, Nx, Ny, Nz)
    G2 = G2 / beta
    if tail:
        # Dimensionless form of the window complement: with
        # wn = r * pi / beta for the odd integers r = 2n + 1 - nmat,
        # beta/4 - (1/beta) sum 1/wn^2 = beta * (1/4 - (1/pi^2) sum 1/r^2).
        # Working in r keeps every intermediate O(1) regardless of beta
        # (the wn-based form underflows wn^2 to 0 for extreme beta).
        r = 2.0 * np.arange(nmat) + 1.0 - nmat
        coeff = beta * (0.25 - np.sum(1.0 / r**2) / np.pi**2)
        di = np.arange(norb)
        # G2[i, i, l, l] += coeff for every (i, l), all k
        G2[di[:, None], di[:, None], di[None, :], di[None, :]] += coeff
    # Fail-fast on a non-finite bubble (round-6 review): downstream
    # diagnostics may legitimately SKIP on size, and the solvers would then
    # propagate NaN into saved eigenvalues silently. Scan in bounded chunks
    # (round-9 review): a whole-tensor np.isfinite allocates a G2-sized
    # Boolean mask -- ~6% of the tensor again, enough to tip a large run
    # that otherwise fits into OOM.
    # ravel(order="K") flattens in MEMORY order and shares storage with
    # the (non-contiguous, ufunc-'K'-layout) G2; reshape(-1) would copy
    # the whole tensor here -- worse than the mask it replaced.
    flat = G2.ravel(order="K")
    step = 1 << 20
    for start in range(0, flat.size, step):
        if not np.all(np.isfinite(flat[start:start + step])):
            raise ValueError(
                "G2 (pair bubble) contains non-finite values; check the "
                "input Green function and beta.")
    return G2


def _eliashberg_kernel_fft(V_q, G2, sigma_old, norb):
    """Apply one iteration of the Eliashberg kernel using FFT convolution.

    Supports both simple (2-index) and general (4-index) pairing vertices.

    For the simple case (V_q shape: norb, norb, Nx, Ny, Nz):
        sigma_{il}(k) = -sum_q V_{ij}(q) * G2_{ijlm}(q) * sigma_{jm}(k)

    For the general case (V_q shape: norb, norb, norb, norb, Nx, Ny, Nz):
        sigma_{l1l4}(k) = -sum_q sum_{l2l3} V_{l1l2,l3l4}(q) * F_{l2l3}(q)
        where F_{l2l3}(q) = sum_{l5l6} G_{l2l5}(k-q) sigma_{l5l6}(k-q) G_{l3l6}(q-k)

    Parameters
    ----------
    V_q : ndarray
        Pairing vertex. Either shape (norb, norb, Nx, Ny, Nz) for simple mode,
        or shape (norb, norb, norb, norb, Nx, Ny, Nz) for general mode.
    G2 : ndarray
        Two-particle Green's function, shape (norb, norb, norb, norb, Nx, Ny, Nz).
    sigma_old : ndarray
        Previous gap function, shape (norb, norb, Nx, Ny, Nz).
    norb : int
        Number of orbitals.

    Returns
    -------
    sigma_new : ndarray
        Updated gap function, shape (norb, norb, Nx, Ny, Nz).
    """
    Nx, Ny, Nz = sigma_old.shape[-3], sigma_old.shape[-2], sigma_old.shape[-1]
    nvol = Nx * Ny * Nz

    # G2Sigma contraction: "ijlmpqs, jmpqs -> ilpqs"
    # = matrix-vector product treating (j,m) as contracted index
    # G2: (norb,norb,norb,norb,Nx,Ny,Nz) -> (norb, norb*norb, norb, nvol)
    # sigma: (norb,norb,Nx,Ny,Nz) -> (norb*norb, nvol)
    sigma_flat = sigma_old.reshape(norb * norb, nvol)  # (jm, site)
    G2_flat = G2.reshape(norb, norb * norb, norb, nvol)  # (i, jm, l, site) -- wrong
    # Need (i, l, jm, site) for matmul on jm
    # Actually einsum is: G2[i,j,l,m,s] * sigma[j,m,s] -> result[i,l,s]
    # Reshape G2 to (norb, norb, norb*norb, nvol): G2[i, l, jm, s]
    G2_r = G2.reshape(norb, norb, norb, norb, nvol)
    G2_r = G2_r.transpose(0, 2, 1, 3, 4).reshape(norb, norb, norb * norb, nvol)
    # G2_r[i, l, jm, s], sigma_flat[jm, s]
    # result[i,l,s] = sum_jm G2_r[i,l,jm,s] * sigma_flat[jm,s]
    G2Sigma = np.einsum('iljs,js->ils', G2_r, sigma_flat).reshape(norb, norb, Nx, Ny, Nz)

    if V_q.ndim == 5:
        # Simple mode: V_q is (norb, norb, Nx, Ny, Nz)
        P_r = ifftn(V_q, axes=(-3, -2, -1))
        G2Sigma_r = ifftn(G2Sigma, axes=(-3, -2, -1))

        Sigma_r = P_r * G2Sigma_r
        sigma_new = fftn(Sigma_r, axes=(-3, -2, -1))
        return -sigma_new

    else:
        # General mode: V_q is (norb, norb, norb, norb, Nx, Ny, Nz)
        F_q = G2Sigma  # (norb, norb, Nx, Ny, Nz) = (l2, l3, ...)

        V_r = ifftn(V_q, axes=(-3, -2, -1))
        F_r = ifftn(F_q, axes=(-3, -2, -1))

        # sigma[i,l,s] = sum_{j,k} V_r[i,j,k,l,s] * F_r[j,k,s]
        # Reshape for matmul: V_r -> (norb, norb^2, norb, nvol), F_r -> (norb^2, nvol)
        V_r_flat = V_r.reshape(norb, norb * norb, norb, nvol)
        F_r_flat = F_r.reshape(norb * norb, nvol)
        sigma_r = np.einsum('ijls,js->ils', V_r_flat, F_r_flat).reshape(
            norb, norb, Nx, Ny, Nz)

        sigma_new = fftn(sigma_r, axes=(-3, -2, -1))
        return -sigma_new


def _make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz):
    """Create LinearOperator for the Eliashberg kernel.

    Parameters
    ----------
    Vs_q : ndarray
        Effective pairing interaction. Either shape (norb, norb, Nx, Ny, Nz)
        for simple mode or (norb, norb, norb, norb, Nx, Ny, Nz) for general.
    G2 : ndarray
        Two-particle Green's function.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        System dimensions.

    Returns
    -------
    A : LinearOperator
        Linear operator representing the Eliashberg kernel K.
    vec_size : int
        Size of the flattened vector.
    """
    vec_size = norb * norb * Nx * Ny * Nz
    nvol = Nx * Ny * Nz

    # Backend of the resident invariants (numpy on CPU, cupy when calc_eliashberg
    # parked them on the device). Everything the kernel builds/derives must live
    # on this same backend; assert Vs_q and G2 agree.
    xp = backend.array_module_of(Vs_q)
    if backend.array_module_of(G2) is not xp:
        raise ValueError(
            "_make_kernel_operator: Vs_q and G2 must be on the same backend "
            "(both host or both device).")

    # Precompute invariants that don't change per matvec call:
    # 1. V_r = IFFT(Vs_q) -- used every matvec
    V_r = backend.spatial_ifftn(Vs_q, axes=(-3, -2, -1))

    # 2. Precompute G2 reshaped for contraction
    G2_r = G2.reshape(norb, norb, norb, norb, nvol)
    G2_pre = G2_r.transpose(0, 2, 1, 3, 4).reshape(norb, norb, norb * norb, nvol)

    is_simple = (Vs_q.ndim == 5)

    # General (4-index) mode: the two contractions
    #   G2Sigma[i,l,s] = sum_j G2_pre[i,l,j,s] * sigma_flat[j,s]
    #   sigma_r[i,l,s] = sum_j V_r_flat[i,j,l,s] * F_r_flat[j,s]
    # are, per spatial point s, a small dense matvec over j (size norb*norb).
    # numpy.einsum does NOT lower these to batched BLAS GEMM, so we precompute
    # the operator tensors ONCE in GEMM-friendly (nvol, M=il, K=j) layout and
    # apply them with matmul (s, and the column axis z for matmat, as leading
    # batch axes). This is NOT bit-identical to the einsum (GEMM changes the
    # reduction order) but matches it to ~1e-13.
    if not is_simple:
        # G2_pre is (i, l, j, s); we need (s, (i,l), j).
        G2_gemm = G2_pre.transpose(3, 0, 1, 2).reshape(
            nvol, norb * norb, norb * norb)
        # V_r_flat is (i, j, l, s); the GEMM M-axis (i,l) needs i (axis0) and
        # l (axis2) adjacent, so transpose to (s, i, l, j) -> (s, (i,l), j).
        V_r_flat = V_r.reshape(norb, norb * norb, norb, nvol)
        V_r_gemm = V_r_flat.transpose(3, 0, 2, 1).reshape(
            nvol, norb * norb, norb * norb)

    def matvec(v):
        # scipy hands us a HOST numpy vector; move it onto the invariants'
        # backend before any device op (mixed numpy/cupy operands raise on real
        # CuPy). On the numpy backend this is a no-op.
        v = xp.asarray(v)
        sigma = v.reshape(norb, norb, Nx, Ny, Nz)
        sigma_flat = sigma.reshape(norb * norb, nvol)

        if is_simple:
            G2Sigma = xp.einsum('iljs,js->ils', G2_pre, sigma_flat).reshape(
                norb, norb, Nx, Ny, Nz)
            G2Sigma_r = backend.spatial_ifftn(G2Sigma, axes=(-3, -2, -1))
            Sigma_r = V_r * G2Sigma_r
            sigma_new = backend.spatial_fftn(Sigma_r, axes=(-3, -2, -1))
        else:
            sig = xp.moveaxis(sigma_flat, 1, 0)[..., np.newaxis]
            G2Sigma = (G2_gemm @ sig)[..., 0]            # (s, (i,l))
            G2Sigma = G2Sigma.reshape(nvol, norb, norb).transpose(
                1, 2, 0).reshape(norb, norb, Nx, Ny, Nz)

            F_r = backend.spatial_ifftn(G2Sigma, axes=(-3, -2, -1))
            F_r_flat = F_r.reshape(norb * norb, nvol)
            f = xp.moveaxis(F_r_flat, 1, 0)[..., np.newaxis]
            sigma_r = (V_r_gemm @ f)[..., 0]             # (s, (i,l))
            sigma_r = sigma_r.reshape(nvol, norb, norb).transpose(
                1, 2, 0).reshape(norb, norb, Nx, Ny, Nz)
            sigma_new = backend.spatial_fftn(sigma_r, axes=(-3, -2, -1))

        # Keep complex: for complex hopping (e.g. spin-orbit coupling),
        # non-centrosymmetric models, or chiral gaps the kernel is genuinely
        # complex; projecting to real would discard physical components. Return
        # a HOST array so scipy's ARPACK/power/subspace driver stays on numpy.
        return backend.to_host(-sigma_new).ravel()

    def matmat(B):
        # Batched form of matvec: apply the SAME kernel to every column of the
        # (vec_size, k) block B in a single FFT. The batch/column axis is the
        # LAST axis so the spatial FFT axes keep absolute positions (2, 3, 4)
        # and V_r / G2_pre (no column axis) broadcast over it.
        B = xp.asarray(B)
        if B.ndim == 1:
            return matvec(B)
        k = B.shape[1]
        sigma = B.reshape(norb, norb, Nx, Ny, Nz, k)
        sigma_flat = sigma.reshape(norb * norb, nvol, k)

        if is_simple:
            G2Sigma = xp.einsum('iljs,jsz->ilsz', G2_pre, sigma_flat).reshape(
                norb, norb, Nx, Ny, Nz, k)
            G2Sigma_r = backend.spatial_ifftn(G2Sigma, axes=(2, 3, 4))
            Sigma_r = V_r[..., np.newaxis] * G2Sigma_r
            sigma_new = backend.spatial_fftn(Sigma_r, axes=(2, 3, 4))
        else:
            sig = xp.moveaxis(sigma_flat, 1, 0)          # (nvol, norb*norb, k)
            G2Sigma = G2_gemm @ sig                      # (s, (i,l), k)
            G2Sigma = G2Sigma.reshape(nvol, norb, norb, k).transpose(
                1, 2, 0, 3).reshape(norb, norb, Nx, Ny, Nz, k)

            F_r = backend.spatial_ifftn(G2Sigma, axes=(2, 3, 4))
            F_r_flat = F_r.reshape(norb * norb, nvol, k)
            f = xp.moveaxis(F_r_flat, 1, 0)              # (nvol, norb*norb, k)
            sigma_r = V_r_gemm @ f                       # (s, (i,l), k)
            sigma_r = sigma_r.reshape(nvol, norb, norb, k).transpose(
                1, 2, 0, 3).reshape(norb, norb, Nx, Ny, Nz, k)
            sigma_new = backend.spatial_fftn(sigma_r, axes=(2, 3, 4))

        return backend.to_host(-sigma_new).reshape(vec_size, k)

    A = LinearOperator((vec_size, vec_size), matvec=matvec, matmat=matmat,
                       dtype=complex)
    return A, vec_size
