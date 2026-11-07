"""The eigenvalue layer of hwave_sc: the leading-eigenpair driver, ARPACK and shift helpers, parity sectors and the static solver drivers.

Moved verbatim from hwave.sc (dev/sc-split-plan.md, step 2); the shell
re-exports every name here under its old spelling.
"""
from __future__ import annotations

import logging
import numpy as np
from scipy.sparse.linalg import LinearOperator, eigs, bicgstab, gmres, lgmres

from hwave.solver.kgrid import reverse_fft_axes
from hwave.solver.eliashberg.config import _validate_spectral_shift
from hwave.solver.eliashberg.kernel import _make_kernel_operator

logger = logging.getLogger("hwave_sc")


# Seed test, on K itself: the seed is an eigenvector of K when its residual
# ||K v - lambda v|| / (spectral scale * ||v||) is below _SEED_EIGENVECTOR_TOL,
# and is
# annihilated by K when ||K v|| <= _SEED_ANNIHILATED_REL * (spectral scale) *
# ||v||. 1e-8 is far below any seed that carries a genuine admixture of other
# modes (those reach the residual of their admixture), yet far above the
# round-off residual (~1e-14) of an analytic form factor that happens to be
# an exact eigenvector by symmetry; between the two the rate-aware stopping
# rule takes over.
_SEED_EIGENVECTOR_TOL = 1.0e-8


# Largest operator dimension for which a failed ARPACK largest-real-part solve
# falls back to a dense eigendecomposition (vec_size matvecs plus an O(n^3)
# eig with vectors: several seconds and a few hundred MB at n = 2000). Above
# it the leading eigenvalue is NaN.
_DENSE_FALLBACK_MAX = 2000


# A sign-alternating two-cycle is RECORDED (for the diagnosis of a run that
# reaches max_iter) when, over the last _TWO_CYCLE_WINDOW steps, every flip
# residual ||u_n + u_{n-1}|| stays below convergence_tol and the even and odd
# iterates have individually converged (||u_n - u_{n-2}|| < convergence_tol).
# It never stops the run: no finite window rules out a tiny component of a
# positive dominant mode that is still growing.
_TWO_CYCLE_WINDOW = 8


# A converged shifted iteration whose eigenvalue on K is below this fraction
# of the shift (|lambda| <= _COLLAPSE_REL * sigma) found a vector K annihilates
# (to round-off), not a pairing eigenvalue: reported as NaN, not converged.
_COLLAPSE_REL = 1.0e-10


_NULL_PROBE_COUNT = 3


_DIFF_ROUNDOFF = 1.0e-13


_SEED_ANNIHILATED_REL = 1.0e-10


_RATE_CLIP = 0.9999


_AUTO_SHIFT_MARGIN = 1.05


# fixed seed of the deterministic probe vector of the spectral scale
_SEED_ADMIX_SEED = 20261101


# Rate-aware stopping of the power iteration: the remaining distance to the
# fixed point is estimated as diff / (1 - rate), with rate the largest ratio
# of consecutive iterate differences over the last _RATE_WINDOW steps,
# clipped at _RATE_CLIP; a rate >= 1 (differences not shrinking: a drift)
# never stops the run, and a difference below _DIFF_ROUNDOFF is accepted
# outright (round-off level, where the ratios are noise).
_RATE_WINDOW = 5


_RAYLEIGH_RESIDUAL_TOL_FACTOR = 100.0


# Absolute floor of the Rayleigh residual tolerance, and the factor applied to
# the caller's convergence_tol. The bound is scaled by (|lambda| + 1) so it is
# relative for a large eigenvalue and absolute for a tiny one.
_RAYLEIGH_RESIDUAL_ATOL = 1.0e-6


# spectral_shift="auto" on the POWER-ITERATION path (see
# _resolve_iteration_shift): number of probe matvecs per probe vector used to
# estimate the spectral radius rho(K), and the safety factor applied to it.
# sigma = _AUTO_SHIFT_MARGIN * rho_est + 1e-6 sits just above the estimated
# radius, which is the smallest shift that makes every eigenvalue of K + sigma*I
# have a positive real part (so the power iteration converges to the
# algebraically largest eigenvalue of K). The margin absorbs the fact that a
# power-iteration estimate approaches rho FROM BELOW; it is deliberately small
# because a larger shift compresses the relative eigenvalue gaps and therefore
# slows the iteration down.
_AUTO_SHIFT_PROBE_STEPS = 15


# Number of independent random probes used by _solve_iteration to test whether
# the Eliashberg kernel commutes with parity before enabling gap-parity
# projection. More than one guards against a single probe accidentally lying
# near the nullspace of the commutator and missing a non-centrosymmetric kernel.
_PARITY_GUARD_PROBES = 3


#: Seed of the deterministic ARPACK start vector (GitHub issue #255).
_ARPACK_START_SEED = 255


# The odd (triplet-parity) harmonics the invariant-subspace tracker reports the
# tracked state on (spec S7.7): {sin kx, sin ky, sin kx (cos kx - cos ky),
# sin ky (cos ky - cos kx)}. They are exactly the odd `_initialize_gap` form
# factors, so the seed table stays the single source of truth.
_ODD_HARMONICS = ("p_x", "p_y", "f_x", "f_y")


# An operator is taken to be exactly zero only when ARPACK stops with error
# -9 ("starting vector is zero") AND it maps _NULL_PROBE_COUNT independent
# random vectors (a generator of its own, so an operator that merely
# annihilates the deterministic ARPACK start vector does not qualify) to
# exactly zero, elementwise, with no tolerance: a tiny but nonzero operator
# is a failed solve like any other.
_NULL_PROBE_SEED = 917253


class _NonFiniteEigenvalues(ArithmeticError):
    """An ARPACK call returned NaN or infinite eigenvalues; treated like a
    solve that did not converge."""


def _is_arpack_zero_start_error(exc):
    """True for ARPACK error -9 ("starting vector is zero"). The auto
    shift's radius probe wraps the ARPACK error in _ShiftProbeFailed, whose
    message is the wrapped error's, so the message test covers it too."""
    import re
    m = re.match(r"ARPACK error (-?\d+)", str(exc))
    return bool(m) and int(m.group(1)) == -9


def _null_flag(record):
    """The ``null_operator`` entry of a failed-solve record, if any."""
    return ({"null_operator": True} if record.get("null_operator")
            else {})


def _operator_is_exactly_zero(A, n):
    """True when ``A`` maps _NULL_PROBE_COUNT independent random vectors to
    exactly zero, elementwise."""
    rng = np.random.default_rng(_NULL_PROBE_SEED)
    for _ in range(_NULL_PROBE_COUNT):
        x = rng.standard_normal(n) + 1j * rng.standard_normal(n)
        if np.any(np.asarray(A.matvec(x)).ravel() != 0):
            return False
    return True


def _initialize_gap(mode, norb, kx_array, ky_array, kz_array):
    """Initialize gap function sigma(k) with specified symmetry.

    Supports 2D and 3D lattice symmetries. The form factors are defined
    in terms of kx, ky, kz and work for any dimension (Nz=1 for 2D).

    Parameters
    ----------
    mode : str
        Initialization mode specifying the gap symmetry.
        Abbreviations cx=cos(kx), sx=sin(kx), etc.

        **Isotropic / s-wave:**
        - "cos"       : cos(kx+ky+kz) (original reference code)
        - "s"         : 1 (constant, k-independent)
        - "s_ext"     : cx*cy + cy*cz + cz*cx (extended s-wave, 3D)
        - "s_ext_2d"  : cx*cy (extended s-wave, 2D)

        **d-wave:**
        - "d_x2y2"    : cx - cy
        - "d_y2z2"    : cy - cz (in-plane d_{y^2-z^2}, even parity / singlet;
                        opposite-sign anti-nodes at (pi,0) and (0,pi))
        - "d_xy"      : sx*sy
        - "d_xz"      : sx*sz
        - "d_yz"      : sy*sz
        - "d_z2"      : 2*cz - cx - cy (d_{3z^2-r^2})

        **p-wave (odd parity):**
        - "p_x"       : sx
        - "p_y"       : sy
        - "p_z"       : sz

        **f-wave (odd parity, higher harmonic):**
        - "f_x"       : sx*(cx - cy)
        - "f_y"       : sy*(cy - cx)

        Note that "f_x"/"f_y" span an odd 2D representation together with
        "p_x"/"p_y" on a square lattice, so they are seeds (and the harmonic
        basis of the subspace tracking, see ``_odd_harmonic_basis``), not a
        point-group separation.

        **Other:**
        - "random"    : random (all symmetries mixed)
    norb : int
        Number of orbitals.
    kx_array, ky_array, kz_array : ndarray
        k-point arrays.

    Returns
    -------
    sigma : ndarray
        Normalized initial gap function, shape (norb, norb, Nx, Ny, Nz).
    """
    Nx, Ny, Nz = len(kx_array), len(ky_array), len(kz_array)
    kx_mesh, ky_mesh, kz_mesh = np.meshgrid(
        kx_array, ky_array, kz_array, indexing='ij'
    )
    I = np.identity(norb)
    cx, cy, cz = np.cos(kx_mesh), np.cos(ky_mesh), np.cos(kz_mesh)
    sx, sy, sz = np.sin(kx_mesh), np.sin(ky_mesh), np.sin(kz_mesh)

    form_factors = {
        # s-wave
        "cos":      lambda: np.cos(kx_mesh + ky_mesh + kz_mesh),
        "s":        lambda: np.ones((Nx, Ny, Nz)),
        "s_ext":    lambda: cx * cy + cy * cz + cz * cx,
        "s_ext_2d": lambda: cx * cy,
        # d-wave
        "d_x2y2":   lambda: cx - cy,
        "d_y2z2":   lambda: cy - cz,
        "d_xy":     lambda: sx * sy,
        "d_xz":     lambda: sx * sz,
        "d_yz":     lambda: sy * sz,
        "d_z2":     lambda: 2.0 * cz - cx - cy,
        # p-wave
        "p_x":      lambda: sx,
        "p_y":      lambda: sy,
        "p_z":      lambda: sz,
        # f-wave (odd; the higher odd harmonic of spec S7.7)
        "f_x":      lambda: sx * (cx - cy),
        "f_y":      lambda: sy * (cy - cx),
    }

    if mode in form_factors:
        f_k = form_factors[mode]()
        sigma = I[:, :, None, None, None] * f_k[None, None, :, :, :]
    elif mode == "random":
        rng = np.random.default_rng(12345)
        sigma = rng.standard_normal((norb, norb, Nx, Ny, Nz))
    else:
        raise ValueError(
            "Unknown init_gap mode: '{}'. Available: {}".format(
                mode, list(form_factors.keys()) + ["random"]))

    norm = np.linalg.norm(sigma)
    if norm > 0:
        sigma /= norm
    elif mode != "random":
        logger.warning(
            "init_gap='%s' evaluates to zero on the current k-grid: a form "
            "factor built from sin(k_a) vanishes when the a-axis is squashed "
            "(e.g. p_x/d_xy at Nx=1, or p_z/d_xz at Nz=1). The seed is all-"
            "zero -- choose a symmetry that is non-vanishing on this grid "
            "(e.g. p_y/p_z, d_yz or d_y2z2 for an in-plane [1, Ny, Nz] cell).",
            mode)
    return sigma


def _odd_harmonic_basis(norb, kx_array, ky_array, kz_array):
    """Normalized odd harmonic basis for the subspace decomposition (S7.7).

    Returns ``{name: flat gap vector}`` for the odd form factors of
    ``_initialize_gap`` (``p_x``, ``p_y``, ``f_x``, ``f_y``), each raveled in
    the gap layout ``(norb, norb, Nx, Ny, Nz)``. Harmonics that vanish
    identically on the current k-grid (e.g. ``p_y`` when ``Ny == 1``) are
    omitted rather than returned as a zero vector, so a decomposition never
    divides by zero.
    """
    basis = {}
    for name in _ODD_HARMONICS:
        sigma = _initialize_gap(name, norb, kx_array, ky_array, kz_array)
        if np.linalg.norm(sigma) == 0.0:
            logger.debug(
                "odd harmonic '%s' vanishes on this k-grid; omitted from the "
                "subspace decomposition basis", name)
            continue
        basis[name] = sigma.ravel()
    return basis


def _solve_iteration(green_kw, Vs_q, G2, sigma_init, norb,
                     max_iter=1000, alpha=0.5, tol=1.0e-5, pairing_type=None,
                     operator=None, spectral_shift=None, info_out=None):
    """Solve linearized Eliashberg equation by self-consistent iteration.

    Parameters
    ----------
    green_kw : ndarray
        Green's function, shape (norb, norb, Nx, Ny, Nz, nmat).
    Vs_q : ndarray
        Effective pairing interaction. Either shape (norb, norb, Nx, Ny, Nz)
        for simple mode or (norb, norb, norb, norb, Nx, Ny, Nz) for general.
    G2 : ndarray
        Two-particle Green's function.
    sigma_init : ndarray
        Initial gap function.
    norb : int
        Number of orbitals.
    max_iter : int
        Maximum iterations.
    alpha : float
        Mixing parameter (0 < alpha < 1).
    tol : float
        Convergence tolerance.
    pairing_type : str, optional
        "singlet" or "triplet". When given, every iterate is projected onto the
        channel's parity sector (even for singlet, odd for triplet). For a
        centrosymmetric system the kernel commutes with parity, so in exact
        arithmetic the iteration would stay in the seed's sector; the projection
        removes the wrong-parity component that numerical noise (FFT/GEMM
        round-off) reintroduces each step and that the power iteration would
        otherwise amplify toward the dominant eigenpair of the *other* channel.
        The projection is only legitimate when the kernel commutes with parity:
        a cross-sector leakage test (several random probes) is run up front and,
        for a non-centrosymmetric kernel (complex/SOC hopping, vertex odd in k,
        chiral gaps), projection is automatically disabled with a warning and the
        un-projected iteration is used. When None, no projection is applied (the
        historical behavior).
    operator : tuple, optional
        Pre-built ``(A, vec_size)`` kernel to iterate with, bypassing
        ``_make_kernel_operator(Vs_q, G2, ...)`` (``Vs_q``/``G2`` may then be
        None). Used by the bond-channel path, whose operator is built by
        ``_build_bond_operator``. When None (the default) the
        historical operator is built exactly as before.
    spectral_shift : float or "auto", optional
        Positive spectral shift sigma; the power iteration then runs on
        ``K + sigma*I`` and the returned eigenvalue is shifted back (see
        ``_solve_leading``). Needed whenever the kernel is repulsive-dominant
        (negative dominant eigenvalue), where the unshifted iterate flips sign
        every step and can never converge. When None (the default of this
        internal function) the unshifted iteration runs; the ``[eliashberg]``
        configuration defaults to ``"auto"`` (``_resolve_spectral_shift``).
    info_out : dict, optional
        When given, updated in place with ``_solve_leading``'s info dict --
        including, on a shifted run, the Rayleigh validation record
        (``eigenvalue_validated`` / ``rayleigh_residual`` / ...). The return
        arity is deliberately unchanged (many call sites unpack exactly four
        values), so this is the way a caller learns WHAT the returned number
        is.

    Returns
    -------
    sigma : ndarray
        Converged gap function.
    eigenvalue : float
        Leading eigenvalue: the SIGNED Rayleigh eigenvalue of the converged
        iterate, or NaN when the iteration did not converge (see
        ``_solve_leading`` / ``info_out``).
    converged : bool
        Whether convergence was achieved (the ``eigen_converged`` verdict).
    n_iter : int
        Number of iterations performed.
    """
    # Precompute the Eliashberg kernel operator once. The vertex IFFT
    # (V_r = ifftn(Vs_q)) and the G2 reshape/transpose are invariant across
    # iterations (only sigma_old changes); _make_kernel_operator hoists them out
    # of the per-matvec closure, saving one full vertex IFFT (and the G2
    # preprocessing) on every one of the up-to-max_iter iterations. The matvec
    # is numerically identical to _eliashberg_kernel_fft(Vs_q, G2, sigma, norb).
    Nx, Ny, Nz = sigma_init.shape[-3], sigma_init.shape[-2], sigma_init.shape[-1]
    A, vec_size = (_make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)
                   if operator is None else operator)
    shape = (norb, norb, Nx, Ny, Nz)

    # Parity projection is only legitimate when the kernel commutes with the
    # parity operator P (Delta_{ab}(k) -> Delta_{ba}(-k), i.e.
    # _reverse_k_and_orbital), which holds for centrosymmetric systems. For a
    # non-centrosymmetric kernel (complex/SOC hopping, a vertex odd in k, chiral
    # gaps) [A, P] != 0; projecting each iterate would then discard physical
    # components and steer the power iteration to the wrong eigenpair. We
    # measure this directly as the cross-sector leakage of the kernel: apply A
    # to even- and odd-projected random probes and ask how much of the result
    # lands in the *opposite* parity sector. Zero leakage <=> [A, P] = 0 <=>
    # projection is valid. Several independent probes are used (one probe could
    # accidentally lie near the nullspace of [A, P] and miss a non-commuting
    # kernel). If any probe leaks above threshold, warn and fall back to the
    # historical un-projected iteration.
    do_project = pairing_type is not None
    if do_project and pairing_type not in ("singlet", "triplet"):
        # Validate up front so an invalid channel always raises, even if the
        # centrosymmetry guard below later disables projection (in which case
        # _project_gap_parity -- which also validates -- is never reached).
        raise ValueError(
            "Unknown pairing_type: '{}'. Use 'singlet' or 'triplet'.".format(
                pairing_type))
    if do_project:
        rng = np.random.default_rng(0)
        leakage = 0.0
        for _ in range(_PARITY_GUARD_PROBES):
            x = (rng.standard_normal(vec_size)
                 + 1j * rng.standard_normal(vec_size)).reshape(shape)
            for sign in (1.0, -1.0):  # even (+) and odd (-) projected probes
                xp = 0.5 * (x + sign * _reverse_k_and_orbital(x))
                axp = A.matvec(xp.ravel()).reshape(shape)
                # component of A xp in the OPPOSITE parity sector
                opp = 0.5 * (axp - sign * _reverse_k_and_orbital(axp))
                denom = np.linalg.norm(axp) + 1.0e-300
                leakage = max(leakage, np.linalg.norm(opp) / denom)
        if leakage > 1.0e-8:
            logger.warning(
                "Eliashberg kernel does not commute with parity (cross-sector "
                "leakage = {:.2e}); the system appears non-centrosymmetric, so "
                "parity projection for the '{}' channel is disabled and the "
                "un-projected iteration is used.".format(leakage, pairing_type))
            do_project = False

    def _project(sigma):
        return (_project_gap_parity(sigma, pairing_type)
                if do_project else sigma)

    sigma_old = _project(sigma_init.copy())
    if do_project and np.linalg.norm(sigma_old) < 1.0e-12 * (
            np.linalg.norm(sigma_init) + 1.0e-300):
        raise ValueError(
            "Initial gap has no component in the '{}' parity sector; "
            "choose an init_gap of the matching parity (even for singlet, "
            "odd for triplet).".format(pairing_type))

    # The power-iteration loop itself (matvec, projection, normalize,
    # convergence check, mixing) is factored into _solve_leading so a future
    # dynamic-frequency kernel can reuse the identical driver. `A` was already
    # built above (for the leakage probe), so `make_operator` just hands back
    # that same operator rather than rebuilding it. The gap-parity projector
    # operates on the (norb, norb, Nx, Ny, Nz) tensor, so it is wrapped as a
    # flat-vector -> flat-vector `project_fn` for the generic driver; when no
    # projection applies (`do_project` False) `project_fn` is None, matching
    # the historical (unprojected) iteration exactly.
    make_operator = lambda: (A, vec_size)
    project_fn = ((lambda flat: _project(flat.reshape(shape)).ravel())
                 if do_project else None)

    eigenvalue, sigma_flat, info = _solve_leading(
        make_operator, vec_size, "iteration",
        max_iter=max_iter, convergence_tol=tol, init_vec=sigma_old.ravel(),
        alpha=alpha, project_fn=project_fn, spectral_shift=spectral_shift,
    )

    if info_out is not None:
        info_out.update(info)
    return sigma_flat.reshape(shape), eigenvalue, info["converged"], info["n_iter"]


def _arpack_start(A, v0=None):
    """The start vector of an ARPACK (``scipy.sparse.linalg.eigs``) call on
    the operator ``A``: ``v0`` when the caller has one (an
    eigenvector-continuation seed), else a FIXED-seed pseudo-random vector
    of ``A``'s dtype (GitHub issue #255).

    Without a start vector ARPACK draws a random one on every call, so the
    converged eigenvalue moved by ~1e-9 relative from run to run on the
    same input. A constant vector would be reproducible too, but it carries
    the full lattice symmetry, and the Krylov space it starts stays in that
    sector -- it could miss, e.g., a d-wave leading eigenvalue -- so the
    vector is generic, only its seed is fixed."""
    if v0 is not None:
        return v0
    n = int(A.shape[0])
    dtype = np.dtype(A.dtype)
    rng = np.random.default_rng(_ARPACK_START_SEED)
    if np.issubdtype(dtype, np.complexfloating):
        v = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    else:
        v = rng.standard_normal(n)
    if np.issubdtype(dtype, np.inexact):
        v = v.astype(dtype)
    # (an integer operator keeps the float64 vector: eigs promotes it to double)
    return v / np.linalg.norm(v)


def _order_by_seed_overlap(vals, vecs, seed_vec):
    """Order eigenpairs by descending overlap ``|<seed, vec>|`` with a seed
    eigenvector (columns already L2-normalized by ARPACK).

    Used for eigenvector continuation: when a converged eigenvector from a
    neighbouring parameter (e.g. the next temperature) is supplied as a seed,
    the physical branch is the eigenpair whose eigenvector maximally overlaps
    it -- NOT the algebraically largest one (which, near an exceptional point
    of the non-Hermitian kernel, can jump to a different branch). Ties and a
    zero seed fall back to real-part ordering.
    """
    s = np.asarray(seed_vec).ravel()
    ns = np.linalg.norm(s)
    if ns == 0:
        return _order_eigenpairs(vals, vecs)
    s = s / ns
    ov = np.abs(vecs.conj().T @ s) / (
        np.linalg.norm(vecs, axis=0) + 1.0e-300)
    # Primary key: descending overlap. Secondary key: descending real part, so
    # an exact overlap TIE deterministically falls back to the physical
    # real-part ordering (as the docstring promises) instead of the arbitrary
    # order np.linalg.eig / ARPACK happened to return. np.lexsort takes the
    # primary key last.
    idx = np.lexsort((-vals.real, -ov))
    return vals[idx], vecs[:, idx]


def _order_eigenpairs(vals, vecs):
    """Order eigenpairs by descending real part (largest first).

    The superconducting eigenvalue is the algebraically largest one
    (lambda -> 1 at Tc); only positive eigenvalues are physically relevant.
    The leading eigenpair must therefore be the largest by real part, NOT the
    largest by magnitude: a large-magnitude *negative* eigenvalue is an
    unphysical repulsive mode and must not be reported first (nor mask a
    smaller positive eigenvalue).

    Note: ARPACK still finds eigenvalues by magnitude (which='LM'), so a small
    positive eigenvalue masked by much larger negative ones may not be among
    the returned set; request more eigenvalues, or use a shift-invert method
    targeting a positive shift, to resolve such cases.

    Parameters
    ----------
    vals : ndarray
        Eigenvalues.
    vecs : ndarray
        Eigenvectors as columns, shape (vec_size, n).

    Returns
    -------
    vals, vecs ordered so that vals[0] has the largest real part.
    """
    idx = np.argsort(-vals.real)
    return vals[idx], vecs[:, idx]


def _reverse_k_and_orbital(gap):
    """Return the gap evaluated at (-k) with orbital indices transposed.

    For the pairing gap Delta_{ab}(k), the singlet/triplet parity condition is
    Delta_{ab}(k) = +/- Delta_{ba}(-k). This helper builds Delta_{ba}(-k):
    it reverses each spatial axis (k -> -k on the FFT grid, index i -> (N-i)%N)
    and swaps the two orbital indices.

    Parameters
    ----------
    gap : ndarray
        Gap function, shape (norb, norb, Nx, Ny, Nz).

    Returns
    -------
    ndarray
        Delta_{ba}(-k), same shape.
    """
    rev = reverse_fft_axes(gap, (2, 3, 4))  # k -> -k on the FFT grid
    rev = np.swapaxes(rev, 0, 1)            # orbital transpose a <-> b
    return rev


def _project_gap_parity(gap, pairing_type):
    """Project a gap onto the parity sector of the pairing channel.

    Singlet gaps are even (Delta_{ab}(k) = +Delta_{ba}(-k)); triplet gaps are
    odd (Delta_{ab}(k) = -Delta_{ba}(-k)). With P the parity operator
    (P: Delta_{ab}(k) -> Delta_{ba}(-k), built by _reverse_k_and_orbital), the
    even/odd projectors are (1 +/- P)/2. The Eliashberg kernel commutes with P,
    so these projectors commute with it and the power iteration can be confined
    to the physical sector -- removing wrong-parity components that numerical
    noise would otherwise let the iteration amplify (see _solve_iteration).

    Parameters
    ----------
    gap : ndarray
        Gap function, shape (norb, norb, Nx, Ny, Nz).
    pairing_type : str
        "singlet" (even projector) or "triplet" (odd projector).

    Returns
    -------
    ndarray
        The component of ``gap`` in the channel's parity sector, same shape.
    """
    if pairing_type == "singlet":
        sign = 1.0
    elif pairing_type == "triplet":
        sign = -1.0
    else:
        raise ValueError(
            "Unknown pairing_type: '{}'. Use 'singlet' or 'triplet'.".format(
                pairing_type))
    return 0.5 * (gap + sign * _reverse_k_and_orbital(gap))


def _is_gap_parity(gap, pairing_type, tol=0.9):
    """Test whether a gap has the parity of the pairing channel.

    Singlet gaps are even (Delta_{ab}(k) = +Delta_{ba}(-k)); triplet gaps are
    odd (Delta_{ab}(k) = -Delta_{ba}(-k)). Returns True when the projection of
    the gap onto the requested parity sector retains at least ``tol`` of its
    norm.

    Parameters
    ----------
    gap : ndarray
        Gap function, shape (norb, norb, Nx, Ny, Nz).
    pairing_type : str
        "singlet" or "triplet".
    tol : float
        Minimum fraction of the norm that must survive the parity projection.

    Returns
    -------
    bool
    """
    proj = _project_gap_parity(gap, pairing_type)
    n = np.linalg.norm(gap)
    if n == 0:
        return False
    return np.linalg.norm(proj) / n >= tol


def _reorder_eigenpairs_by_parity(vals, gaps, pairing_type, warn=True):
    """Stable-reorder eigenpairs so those matching the channel parity come first.

    The Eliashberg kernel preserves parity, so its eigenvectors split into even
    (singlet) and odd (triplet) sectors. When solving for a given channel, the
    physical solution is the leading eigenpair of the matching parity, which is
    not necessarily the globally leading one. This keeps every eigenpair (so the
    returned count is unchanged) but promotes the requested-parity ones,
    preserving their existing order.

    Parameters
    ----------
    vals : ndarray
        Eigenvalues, shape (n,).
    gaps : ndarray
        Eigenvectors as gaps, shape (n, norb, norb, Nx, Ny, Nz).
    pairing_type : str
        "singlet" or "triplet".

    Returns
    -------
    vals, gaps reordered with matching-parity eigenpairs first.
    """
    match = np.array([_is_gap_parity(g, pairing_type) for g in gaps])
    if warn and not np.any(match):
        logger.warning(
            "No computed eigenpair matches the requested '%s' parity; the "
            "reported leading gap belongs to the other channel. Increase "
            "num_eigenvalues or check the pairing_type.", pairing_type)
    idx = np.concatenate([np.where(match)[0], np.where(~match)[0]])
    return vals[idx], gaps[idx]


def _parity_sector_basis(shape, pairing_type):
    """Orthonormal coordinates of the channel's parity sector.

    The parity map ``R: (a, b, k) -> (b, a, -k)`` (``_reverse_k_and_orbital``)
    is an involution of the flat gap index. Even (singlet) sector: one
    coordinate ``(e_i + e_Ri)/sqrt(2)`` per orbit ``i < R i`` plus ``e_i`` for
    each self-paired ``i = R i``; odd (triplet) sector: ``(e_i - e_Ri)/
    sqrt(2)`` per orbit (self-paired elements vanish). Returns ``(to_full,
    to_reduced, dim)``."""
    n = int(np.prod(shape))
    partner = np.ascontiguousarray(
        _reverse_k_and_orbital(np.arange(n).reshape(shape))).ravel()
    flat = np.arange(n)
    rep = flat[flat < partner]
    rep_p = partner[rep]
    self_idx = flat[flat == partner] if pairing_type == "singlet" else \
        np.zeros(0, dtype=int)
    sign = 1.0 if pairing_type == "singlet" else -1.0
    npair, nself = len(rep), len(self_idx)
    r2 = np.sqrt(0.5)

    def to_full(c):
        c = np.asarray(c).ravel()
        g = np.zeros(n, dtype=complex)
        g[rep] = r2 * c[:npair]
        g[rep_p] = sign * r2 * c[:npair]
        g[self_idx] = c[npair:]
        return g

    def to_reduced(g):
        g = np.asarray(g).ravel()
        return np.concatenate([r2 * (g[rep] + sign * g[rep_p]), g[self_idx]])

    return to_full, to_reduced, npair + nself


def _listing_with_leader(vals, gaps, lead, lead_gap):
    """The eigen-analysis table with the reported channel leader as row 0.

    The sector leader becomes row 0 (its vector has the requested parity, so
    its match is 1); rows 1.. are the global eigen-analysis set exactly as
    the solver returned it (``num_eigenvalues`` rows). Nothing is merged:
    when the leader's state is also among the global rows it is listed
    twice. Every emitted vector therefore stays paired with the eigenvalue
    the solver computed for it."""
    vals = np.asarray(vals)
    gaps = np.asarray(gaps)
    shape = gaps.shape[1:]
    out_vals = np.concatenate([np.asarray([lead], dtype=np.result_type(
        vals.dtype, np.asarray(lead).dtype)), vals])
    out_gaps = np.concatenate([np.asarray(lead_gap).reshape((1,) + shape),
                               gaps])
    return out_vals, out_gaps


def _channel_eigenpairs(vals, gaps, operator, shape, pairing_type,
                        num_eigenvalues, spectral_shift=None):
    """The reported leader of the requested pairing channel (static,
    non-spin-orbital).

    The global eigen-analysis set (``vals``/``gaps``) is kept for the
    listing, reordered so the eigenpairs of the requested parity come first.
    It does not by itself identify the channel's leader: the leading
    ``num_eigenvalues`` eigenpairs of the whole kernel may contain a channel
    eigenpair of large magnitude but not the channel's leader. When the
    kernel commutes with parity (probed as in ``_solve_iteration``), the
    leader is therefore ALWAYS obtained from a solve on the channel's own
    parity coordinates (``_parity_sector_basis``: no complement, no penalty),
    with the shared largest-real-part contract of ``_solve_leading``
    (LR re-solve, dense fallback, NaN sentinel). When it does not commute,
    the sector cannot be solved on its own; the leader is then the retained
    eigenpair that is dominantly of the requested parity, or NaN if there is
    none.

    The sector solve uses the user's ``spectral_shift`` (largest real part
    of the shifted sector operator) and otherwise the plain Arnoldi rules; it
    is always an Arnoldi solve, whatever ``eigenvalue_method`` the global
    listing used.

    Returns ``(vals_r, gaps_r, lead, lead_gap, found, note, sector_vals)``:
    the reordered listing, the leading eigenvalue and gap, whether a channel
    leader was obtained, a one-line description of where it came from (or
    None), and the sector solve's eigenvalues (None without a sector
    solve)."""
    from scipy.sparse.linalg import LinearOperator
    A, n = operator
    rng = np.random.default_rng(0)
    leakage = 0.0
    for _ in range(_PARITY_GUARD_PROBES):
        x = (rng.standard_normal(n) + 1j * rng.standard_normal(n)).reshape(shape)
        opp_max, img_max = 0.0, 0.0
        for sign in (1.0, -1.0):
            xp = 0.5 * (x + sign * _reverse_k_and_orbital(x))
            axp = np.asarray(A.matvec(xp.ravel())).reshape(shape)
            opp = 0.5 * (axp - sign * _reverse_k_and_orbital(axp))
            opp_max = max(opp_max, float(np.linalg.norm(opp)))
            img_max = max(img_max, float(np.linalg.norm(axp)))
        # normalised by the LARGER sector image, so an analytically null (or
        # tiny) opposite sector does not read as leakage 1
        leakage = max(leakage, opp_max / (img_max + 1.0e-300))
    nan = complex(np.nan, np.nan)
    commuting = leakage <= 1.0e-8
    vals_r, gaps_r = _reorder_eigenpairs_by_parity(vals, gaps, pairing_type,
                                                   warn=not commuting)
    if not commuting:
        matched = bool(len(gaps_r) and _is_gap_parity(gaps_r[0], pairing_type))
        logger.warning(
            "The Eliashberg kernel does not commute with parity (leakage "
            "%.2e), so the '%s' channel cannot be solved on its own parity "
            "sector; %s", leakage, pairing_type,
            "the leading value is the computed eigenpair that is dominantly of "
            "the requested parity." if matched else
            "no computed eigenpair has the requested parity: no eigenvalue is "
            "reported (NaN, eigen_converged = false). Increase "
            "num_eigenvalues.")
        if matched:
            return vals_r, gaps_r, vals_r[0], gaps_r[0], True, None, None
        return vals_r, gaps_r, nan, None, False, None, None
    to_full, to_reduced, dim = _parity_sector_basis(shape, pairing_type)
    if dim == 0:
        logger.warning("The '%s' parity sector is empty on this grid.",
                       pairing_type)
        return vals_r, gaps_r, nan, None, False, None, None
    red = LinearOperator(
        (dim, dim), matvec=lambda c: to_reduced(A.matvec(to_full(c))),
        dtype=complex)
    lead, cvec, info = _solve_leading(lambda: (red, dim), dim, "arnoldi",
                                      num_eigenvalues=num_eigenvalues,
                                      spectral_shift=spectral_shift)
    if not np.isfinite(lead):
        return vals_r, gaps_r, nan, None, False, None, None
    lead_gap = to_full(cvec).reshape(shape)
    note = ("leading value: the largest real part of the '{}' parity sector, "
            "solved on that sector's own coordinates (selection {}); it is "
            "row 0 of the table, followed by the global eigen-analysis "
            "set".format(
                pairing_type, info.get("selection")))
    return (vals_r, gaps_r, complex(lead), lead_gap, True, note,
            np.asarray(info.get("eigenvalues")))


def _shift_from_eigenvalues(vals, factor=0.9):
    """Estimate a shift-invert target near the physical (largest-real) eigenvalue.

    The superconducting eigenvalue is the algebraically largest one, so when a
    positive eigenvalue is present among the sampled values, aim the shift at
    the largest real part (NOT the largest magnitude, which could be a large
    negative repulsive mode that masks it). If no positive eigenvalue was
    sampled (no SC instability in range), fall back to the largest-magnitude
    eigenvalue so shift-invert still tracks the dominant mode.

    Parameters
    ----------
    vals : ndarray
        Sampled eigenvalues (e.g. a few largest-magnitude ARPACK values).
    factor : float
        Fraction of the target eigenvalue to use as the shift.

    Returns
    -------
    float
    """
    real = vals.real
    scale = float(np.max(np.abs(vals)))
    # "positive" means significantly positive relative to the spectrum scale,
    # so numerical-noise eigenvalues (~1e-16) do not pull the shift to zero.
    if scale > 0 and np.any(real > 1e-8 * scale):
        return float(np.max(real)) * factor
    # No significant positive eigenvalue: track the dominant magnitude. Note a
    # purely imaginary / all-zero sample yields a 0.0 shift, which sits at a
    # possible zero mode; physical FLEX inputs have a real dominant eigenvalue.
    return float(vals[np.argmax(np.abs(vals))].real) * factor


def _normalize_spectral_shift(spectral_shift):
    """``None`` for "no shift": ``None`` itself or the string ``"none"``
    (any case; TOML has no null, so ``spectral_shift = "none"`` is how an
    input asks for the unshifted solver). Anything else is returned as is
    (and validated by ``_validate_spectral_shift``)."""
    if isinstance(spectral_shift, str) and spectral_shift.strip().lower() == "none":
        return None
    return spectral_shift


def _resolve_spectral_shift(spectral_shift, solver_mode):
    """The ``[eliashberg] spectral_shift`` a solver leg actually uses.

    For the power iteration (``solver_mode = "iteration"``) an absent value
    (``None``) means ``"auto"``: an unshifted power iteration cannot converge
    on a repulsive-dominant kernel (its dominant eigenvalue is negative) and,
    even where it converges, follows the largest-MAGNITUDE eigenvalue, not the
    largest real part. ``"none"`` keeps the unshifted iteration. For every
    other solver an absent value stays ``None`` (no shift)."""
    if spectral_shift is None and str(solver_mode).lower() == "iteration":
        return "auto"
    return _normalize_spectral_shift(spectral_shift)


class _ShiftProbeFailed(Exception):
    """The which='LM' spectral-radius probe of the auto spectral_shift did not
    converge (and returned no usable eigenvalue)."""


def _estimate_spectral_radius(A, probes, project_fn=None,
                              n_steps=_AUTO_SHIFT_PROBE_STEPS):
    """Estimate the spectral radius rho(A) with a few power-iteration steps.

    For a unit vector ``v``, ``||A v||`` is a lower bound on rho(A) that
    increases toward it as ``v`` is driven onto the dominant eigenvector; the
    largest value seen over all probes and steps is returned.

    ``project_fn`` (the gap-parity projector, when the caller iterates in a
    parity sector) is applied to every probe iterate, so the estimate is the
    spectral radius OF THAT SECTOR -- the only part of the spectrum the shifted
    iteration will actually see. Several probes are used because a physically
    motivated seed can be orthogonal to the dominant eigenvector (an f-wave
    seed and an s-wave dominant mode, say) and would then wildly underestimate
    rho on its own.
    """
    rho = 0.0
    for v0 in probes:
        v = np.asarray(v0, dtype=complex).ravel()
        if project_fn is not None:
            v = project_fn(v)
        nrm = np.linalg.norm(v)
        if nrm == 0.0 or not np.isfinite(nrm):
            continue
        v = v / nrm
        for _ in range(n_steps):
            w = A.matvec(v)
            if project_fn is not None:
                w = project_fn(w)
            nw = float(np.linalg.norm(w))
            if nw == 0.0 or not np.isfinite(nw):
                break
            rho = max(rho, nw)
            v = w / nw
    return rho


def _resolve_iteration_shift(A, spectral_shift, vec_size, init_vec,
                             project_fn=None):
    """Resolve the power-iteration spectral shift sigma > 0.

    A number is taken as given. ``"auto"`` estimates the spectral radius with
    ``_estimate_spectral_radius`` (probing from the caller's seed AND from a
    deterministic random vector) and returns
    ``_AUTO_SHIFT_MARGIN * rho_est + 1e-6`` -- i.e. slightly ABOVE the
    estimated radius, the smallest shift that pushes the whole spectrum into
    the right half-plane so that the algebraically largest eigenvalue of K
    becomes the dominant (and positive) eigenvalue of K + sigma*I.
    """
    if not isinstance(spectral_shift, str):
        return float(spectral_shift)
    # "auto" (validated by _validate_spectral_shift)
    rng = np.random.default_rng(0)
    probe = (rng.standard_normal(vec_size)
             + 1j * rng.standard_normal(vec_size))
    rho = _estimate_spectral_radius(A, (init_vec, probe), project_fn)
    sig = _AUTO_SHIFT_MARGIN * rho + 1.0e-6
    logger.info(
        "auto spectral_shift: spectral radius estimated as rho = %.6f from 2 "
        "probes x %d power-iteration steps; sigma = %.2f * rho + 1e-6 = %.6f.",
        rho, _AUTO_SHIFT_PROBE_STEPS, _AUTO_SHIFT_MARGIN, sig)
    if not np.isfinite(sig):
        raise ValueError(
            "auto spectral_shift overflowed to a non-finite value (spectral "
            "radius too large); pass an explicit shift.")
    return sig


def _signed_rayleigh(A, vec):
    """Signed Rayleigh quotient and relative residual of ``vec`` under ``A``.

    Returns ``(lambda, residual)`` with
    ``lambda = <v|A|v> / <v|v>`` (complex; it carries the SIGN of the
    eigenvalue, unlike a power-iterate norm) and
    ``residual = ||A v - lambda v|| / ||v||``. Costs exactly one matvec. A
    zero/non-finite vector or quotient yields ``residual = inf``, i.e. "not
    an eigenvector".
    """
    vec = np.asarray(vec).ravel()
    vnorm = float(np.linalg.norm(vec))
    if vnorm == 0.0 or not np.isfinite(vnorm):
        return complex(np.nan, np.nan), float("inf")
    Av = np.asarray(A.matvec(vec)).ravel()
    lam = complex(np.vdot(vec, Av) / (vnorm * vnorm))
    if not (np.isfinite(lam.real) and np.isfinite(lam.imag)):
        return lam, float("inf")
    return lam, float(np.linalg.norm(Av - lam * vec) / vnorm)


def _validate_shifted_eigenvalue(A_unshifted, vec, iterate_value, shift,
                                 convergence_tol, spectral_scale=None,
                                 rho_max=None):
    """Decide what a shifted power iteration is allowed to report.

    Parameters
    ----------
    A_unshifted : LinearOperator
        The ORIGINAL kernel K (not ``K + sigma*I``).
    vec : ndarray
        The vector the iteration returned.
    iterate_value : float
        ``||(K + sigma I) v|| - sigma``, the historical quantity.
    shift : float
        The sigma actually applied.
    convergence_tol : float
        The iteration's own tolerance; the residual bound is derived from it.
    spectral_scale : float, optional
        A spectral scale of K (the auto shift's radius estimate, or a probe
        norm). When given, the residual bound is ``tol * spectral_scale``
        instead of ``tol * (|lambda| + 1)``: a bound that grows with the
        lambda being checked lets a non-normal kernel's off-spectrum
        Rayleigh quotients pass.
    rho_max : float, optional
        An upper bound of the spectral radius; a Rayleigh value beyond it is
        not an eigenvalue and is rejected.

    Returns
    -------
    value : float
        ``lambda`` (the validated Rayleigh eigenvalue) when the residual check
        passes, otherwise ``iterate_value`` unchanged.
    record : dict
        ``eigenvalue_validated`` / ``eigenvalue_kind`` (``"eigenvalue"`` vs
        ``"shifted-iterate-norm-estimate"``) / ``rayleigh_eigenvalue`` /
        ``rayleigh_residual`` / ``rayleigh_residual_tol`` /
        ``shift_sufficient``. ``shift_sufficient`` is False when the mode the
        iteration locked onto has a NEGATIVE shifted eigenvalue -- the
        signature of a too-small sigma -- and None when nothing could be
        validated.
    """
    lam, resid = _signed_rayleigh(A_unshifted, vec)
    tol = max(_RAYLEIGH_RESIDUAL_ATOL,
              _RAYLEIGH_RESIDUAL_TOL_FACTOR * float(convergence_tol))
    if spectral_scale is not None and spectral_scale > 0.0:
        bound = tol * float(spectral_scale)
    else:
        scale = abs(lam)
        bound = tol * ((scale if np.isfinite(scale) else 0.0) + 1.0)
    validated = bool(np.isfinite(resid) and resid <= bound
                     and abs(lam.imag) <= bound
                     and (rho_max is None or abs(lam) <= rho_max))
    record = {
        "eigenvalue_validated": validated,
        "eigenvalue_kind": ("eigenvalue" if validated
                            else "shifted-iterate-norm-estimate"),
        "rayleigh_eigenvalue": (float(lam.real) if validated else None),
        "rayleigh_residual": resid,
        "rayleigh_residual_tol": bound,
        "shift_sufficient": None,
    }
    if not validated:
        return iterate_value, record
    value = float(lam.real)
    # A sufficient sigma puts the WHOLE spectrum in the right half-plane, so
    # the mode the iteration locks onto satisfies lambda + sigma > 0. A
    # negative implied shifted eigenvalue means the iterate was driven by
    # |lambda + sigma| with lambda + sigma < 0: sigma was too small, and the
    # mode found is the largest-MAGNITUDE one, not the algebraically largest.
    record["shift_sufficient"] = bool(value + shift > 0.0)
    return value, record


def _stall_advice(diffs, tol, max_iter, window=20):
    """Advice for a power iteration that ran out of ``max_iter``.

    When the iterate difference was still decreasing geometrically over the
    last ``window`` steps, the observed rate gives the number of further steps
    needed to reach ``tol``; otherwise the iteration was not converging at all
    (a complex leading eigenvalue -- a complex-conjugate pair -- or two
    eigenvalues of equal modulus, which a power iteration never resolves)."""
    d = [x for x in diffs[-(window + 1):] if np.isfinite(x) and x > 0.0]
    if len(d) >= 3:
        k = len(d) - 1
        rate = (d[-1] / d[0]) ** (1.0 / k)
        # the stopping rule needs diff / (1 - rate) < tol (rate-aware)
        target = tol * (1.0 - rate) if 0.0 < rate < 1.0 else tol
        if 0.0 < rate < 1.0 and d[-1] > target:
            more = int(np.ceil(np.log(target / d[-1]) / np.log(rate)))
            return ("The iterate difference was still decreasing by a factor "
                    "{:.4f} per step (now {:.2e}; the stopping rule needs "
                    "diff / (1 - rate) < convergence_tol {:.1e}): "
                    "about {} more steps are needed, i.e. max_iter of about "
                    "{} -- or use solver_mode=\"eigenvalue\", which does not "
                    "depend on the eigenvalue gap this way."
                    .format(rate, d[-1], tol, more, max_iter + more))
        return ("The iterate difference was not decreasing (factor {:.4f} per "
                "step over the last {} steps): the iteration is not converging "
                "at all. A complex leading eigenvalue (a complex-conjugate "
                "pair) or two eigenvalues of equal modulus can never be "
                "resolved by a power iteration; use solver_mode=\"eigenvalue\"."
                .format(rate, k))
    return ("Raise max_iter, or use solver_mode=\"eigenvalue\" (a complex "
            "leading eigenvalue can never be reached by a power iteration).")


def _solve_leading(make_operator, vec_size, solver_mode, num_eigenvalues=10,
                   max_iter=1000, convergence_tol=1.0e-5, init_vec=None,
                   sigma_shift=None, alpha=0.5, project_fn=None, seed_vec=None,
                   spectral_shift=None, warn_negative_leading=True):
    """Shared leading-eigenpair driver behind the static Eliashberg solvers.

    This holds the ARPACK/shift-invert eigen-selection-and-ordering body of
    ``_solve_eigenvalue`` and the power-iteration loop of ``_solve_iteration``,
    generalized to act on any ``make_operator``/``vec_size`` pair (not just the
    static ``_make_kernel_operator`` result) so a future dynamic-frequency
    kernel can reuse the exact same driver.

    Parameters
    ----------
    make_operator : callable
        Zero-argument callable returning ``(A, op_vec_size)``, i.e. the same
        convention as ``_make_kernel_operator`` (a
        ``scipy.sparse.linalg.LinearOperator`` acting on flat vectors of
        length ``vec_size``, plus that size). Called exactly once.
    vec_size : int
        Size of the flattened vector space that ``A`` acts on.
    solver_mode : str
        "arnoldi", "shift-invert-bicgstab", "shift-invert-gmres", or
        "shift-invert-lgmres" select the ARPACK/shift-invert eigenanalysis
        (the former body of ``_solve_eigenvalue``, excluding its "subspace"
        branch which stays in that wrapper). "iteration" selects the power
        loop (the former body of ``_solve_iteration``).
    num_eigenvalues : int
        Number of eigenvalues to request from ARPACK. Only used by the
        eigenvalue-family modes.
    max_iter : int
        Maximum number of power-iteration steps. Only used by "iteration".
    convergence_tol : float
        Power-iteration convergence tolerance on the normalized-iterate
        difference. Only used by "iteration".
    init_vec : ndarray, optional
        Flat initial vector, shape ``(vec_size,)``, for "iteration" mode
        (already projected onto the desired sector by the caller, if
        applicable). Required for "iteration" mode.
    sigma_shift : float, optional
        Shift-invert target for the "shift-invert-*" modes. If None, it is
        estimated from a preliminary Arnoldi pass, exactly as
        ``_solve_eigenvalue`` does.
    alpha : float
        Power-iteration mixing parameter. Only used by "iteration", and only
        without a ``spectral_shift`` (with a shift it is not applied; see
        ``info["alpha_applied"]``).
    project_fn : callable, optional
        Flat-vector -> flat-vector projector applied to every power-iteration
        iterate (e.g. the gap-parity projector in ``_solve_iteration``). Only
        used by "iteration"; when None, no projection is applied.
    spectral_shift : float or "auto", optional
        Positive spectral shift sigma. Supported by "arnoldi" (ARPACK asks for
        the largest-REAL eigenvalue of ``A + sigma*I``) and by "iteration"
        (the power loop runs on ``A + sigma*I``, whose dominant eigenvalue is
        positive even for a repulsive-dominant kernel, so it converges). In
        both cases sigma is subtracted back: the returned eigenvalue is the
        SIGNED eigenvalue of the original operator and the eigenvector is
        unchanged by the shift. ``"auto"`` estimates sigma from the spectral
        radius (see ``_resolve_iteration_shift`` for the iteration path).
        Rejected for every other mode.
    warn_negative_leading : bool, optional
        When True (default), the plain ``arnoldi`` (no ``spectral_shift``, no
        ``seed_vec``) path logs an advisory if the leading eigenvalue comes out
        negative. Pass False when the caller detects and handles that case
        itself (e.g. the dynamic solver's automatic largest-real re-solve), so
        the advisory is not emitted twice.

    Returns
    -------
    leading_eigenvalue : complex or float
        Eigenvalue family: the eigenvalue of largest real part, per
        ``_order_eigenpairs``; an unshifted ``which='LM'`` set that cannot
        contain the largest real part (its largest real part is below its
        smallest magnitude) is re-solved for it (``selection = "LR_retry"``,
        unseeded runs only). Iteration path: the SIGNED Rayleigh eigenvalue
        ``<v|K|v>/<v|v>`` of the converged iterate on the original K (with or
        without ``spectral_shift``), and ``NaN`` whenever the iteration did not
        converge -- by ``max_iter``, by a detected sign-alternating two-cycle
        (a negative dominant eigenvalue of the iterated operator) or by a
        failed residual check. The iterate norm is never returned.
    leading_eigenvector : ndarray
        Flat, shape ``(vec_size,)``.
    eig_analysis : dict
        Eigenvalue-family modes: ``{"eigenvalues": vals, "eigenvectors": vecs,
        "sigma_shift": sigma_shift, "selection": <str>}`` with ``vals`` ordered
        by descending real part (``_order_eigenpairs``), ``vecs`` the matching
        eigenvectors as columns, and ``selection`` the criterion actually used
        -- ``"LM"`` (largest magnitude), ``"LR"`` (largest real part, a
        ``spectral_shift`` was applied), ``"dense-LR"`` (the tiny-operator dense
        eigendecomposition, ordered by real part) or ``"shift-invert"``.
        "iteration" mode: ``converged`` (the ``eigen_converged`` verdict:
        iterate difference below ``convergence_tol`` AND Rayleigh residual
        check passed), ``diff_converged`` (the bare difference test),
        ``n_iter``, ``two_cycle``, ``eigenvalue_validated``,
        ``eigenvalue_kind`` (``"eigenvalue"``, ``"not-converged"``,
        ``"not-computed"`` or ``"shifted-iterate-norm-estimate"``),
        ``rayleigh_eigenvalue``, ``rayleigh_residual``,
        ``rayleigh_residual_tol``, ``shift_sufficient``, on a two-cycle
        ``cycle_eigenvalue`` (the Rayleigh quotient of the cycling iterate on
        K, i.e. the negative dominant eigenvalue), and ``spectral_shift``
        (the sigma actually applied) only when a shift was requested. The
        eigenvalue family's dict may also carry ``lm_max_real`` /
        ``lm_min_abs`` (the bound that triggered the re-solve),
        ``lr_retry_failed``, and ``lm_solve_failed`` (the unshifted or
        preliminary which='LM' solve did not converge; the selection is then
        ``"dense-LR"`` or the NaN sentinel ``"not-converged"``). An exactly
        zero operator (ARPACK error -9 confirmed by independent probes)
        reports eigenvalues 0 with the selection ``"null-operator"`` and
        ``null_operator = True``.
    """
    # Validate spectral_shift up front -- before any branch, including the
    # iteration loop and the small-dense early return -- so an invalid value or
    # an incompatible mode fails fast everywhere.
    spectral_shift = _normalize_spectral_shift(spectral_shift)
    _validate_spectral_shift(spectral_shift, solver_mode)

    if solver_mode == "iteration":
        if init_vec is None:
            raise ValueError("init_vec is required for solver_mode='iteration'")
        A, _ = make_operator()

        # Spectral shift: iterate with K + sigma*I instead of K. A repulsive-
        # dominant kernel has a NEGATIVE dominant eigenvalue, which flips the
        # iterate's sign every step, so the convergence test
        # ||sigma_new/||sigma_new|| - sigma_old|| sits at ~2 and the loop can
        # never converge. Shifting the spectrum into the right half-plane makes
        # the ALGEBRAICALLY LARGEST eigenvalue dominant and positive (the same
        # physical selection as the arnoldi which='LR' path). The shift leaves
        # every eigenvector -- and hence every invariant subspace, including
        # the parity sectors that project_fn selects -- untouched; only the
        # eigenvalue moves, and it is shifted back below.
        shift = 0.0
        if spectral_shift is not None:
            shift = _resolve_iteration_shift(A, spectral_shift, vec_size,
                                             init_vec, project_fn)
            logger.info(
                "Power iteration with spectral shift sigma = %.6f (%s): "
                "iterating on K + sigma*I; the reported eigenvalue is "
                "un-shifted (lambda = lambda_shifted - sigma) and the "
                "eigenvector is unchanged by the shift.",
                shift,
                "auto" if isinstance(spectral_shift, str) else "explicit")
            A_unshifted = A
            A = LinearOperator(
                (vec_size, vec_size),
                matvec=lambda v: A_unshifted.matvec(v) + shift * v,
                dtype=complex)
        # the sigma actually applied (None when no shift was requested), so
        # callers can report it; `shift` itself is the arithmetic value.
        shift_used = shift if spectral_shift is not None else None
        # the operator whose eigenvalue is reported: the ORIGINAL kernel K
        A_report = A_unshifted if shift_used is not None else A
        # The mixing v <- (1 - alpha) u + alpha v damps the iterate like a
        # second spectral shift: it halves the contraction (1 - rate) of the
        # dominant mode at alpha = 0.5. With a spectral shift the eigenvalues
        # of K + sigma*I are already positive and well separated in sign, so
        # the mixing only slows the run; it is applied only without a shift
        # (where it damps the sign flips of a negative dominant mode).
        alpha_applied = float(alpha) if shift_used is None else 0.0
        if shift_used is not None and alpha:
            logger.info("Power iteration with a spectral shift: the mixing "
                        "parameter alpha = %g is not applied.", alpha)

        def _finish(vec, diff_converged, n_iter_done, two_cycle=False,
                    collapsed=False):
            """Assemble the return triple.

            The reported value is the SIGNED Rayleigh quotient
            ``<v|K|v>/<v|v>`` of the returned vector on the original kernel K
            (one extra matvec), and only when the iteration converged AND that
            quotient passes the residual check; in every other case it is NaN
            -- a non-converged power iteration has no eigenvalue to report,
            and its last iterate norm (``|lambda_-| (1 - alpha) / (1 + alpha)``
            on a sign-alternating two-cycle) is not one. ``info["converged"]``
            is that verdict (the ``eigen_converged`` flag of the outputs);
            ``info["diff_converged"]`` is the bare iterate-difference test.
            """
            info = {"converged": False, "diff_converged": bool(diff_converged),
                    "n_iter": n_iter_done, "two_cycle": bool(two_cycle),
                    "alpha_applied": alpha_applied}
            if shift_used is not None:
                info["spectral_shift"] = shift_used
            if collapsed:
                # nothing meaningful to validate (the zero-norm iterate)
                info.update({"eigenvalue_validated": False,
                             "eigenvalue_kind": "not-computed",
                             "rayleigh_eigenvalue": None,
                             "rayleigh_residual": float("inf"),
                             "rayleigh_residual_tol": float("nan"),
                             "shift_sufficient": None})
                return float("nan"), vec, info
            value, record = _validate_shifted_eigenvalue(
                A_report, vec, float("nan"), shift, convergence_tol,
                spectral_scale=spec_scale, rho_max=rho_max)
            if shift_used is None:
                record["shift_sufficient"] = None
            info.update(record)
            if not diff_converged:
                # the vector may still be an eigenvector (a two-cycle iterate
                # is one of the dominant NEGATIVE mode), but it is not the
                # eigenpair the iteration was asked for: report nothing
                info["cycle_eigenvalue"] = (
                    float(value) if two_cycle and record["eigenvalue_validated"]
                    else None)
                info.update({"eigenvalue_validated": False,
                             "eigenvalue_kind": "not-converged",
                             "rayleigh_eigenvalue": None})
                return float("nan"), vec, info
            if not record["eigenvalue_validated"]:
                logger.warning(
                    "Power iteration: the iterate difference fell below "
                    "convergence_tol, but the Rayleigh check <v|K|v>/<v|v> "
                    "left a residual of %.3e (tolerance %.3e), so the vector "
                    "is not an eigenvector of K; no eigenvalue is reported "
                    "(lambda = NaN, eigen_converged = false). Tighten "
                    "convergence_tol or use solver_mode=\"eigenvalue\".",
                    record["rayleigh_residual"], record["rayleigh_residual_tol"])
                return float("nan"), vec, info
            if (shift_used is not None
                    and abs(value) <= _COLLAPSE_REL * max(abs(shift), 1.0e-300)):
                # K annihilates the vector (to round-off): the shifted
                # iteration "converged" on the eigenvalue sigma of K + sigma*I
                # that a null vector of K carries, not on a pairing mode
                logger.warning(
                    "Power iteration: the converged vector is annihilated by "
                    "the kernel (lambda = %.3e, below %.0e of the shift "
                    "%.6f); this is not a pairing eigenvalue, so none is "
                    "reported (lambda = NaN, eigen_converged = false). The "
                    "iteration's sector may carry no other mode; use "
                    "solver_mode=\"eigenvalue\" or another init_gap.",
                    value, _COLLAPSE_REL, shift)
                info.update({"eigenvalue_validated": False,
                             "eigenvalue_kind": "annihilated",
                             "rayleigh_eigenvalue": None})
                return float("nan"), vec, info
            info["converged"] = True
            logger.info(
                "Power iteration: reporting the SIGNED Rayleigh eigenvalue "
                "lambda = %.8e of K (residual %.3e <= %.3e).",
                value, record["rayleigh_residual"],
                record["rayleigh_residual_tol"])
            if record["shift_sufficient"] is False:
                logger.warning(
                    "spectral_shift = %.6f appears INSUFFICIENT: the "
                    "iteration locked onto a mode with lambda = %.8e, so "
                    "lambda + sigma = %.8e < 0 and the power iteration "
                    "selected the largest-MAGNITUDE eigenvalue of "
                    "K + sigma*I rather than the algebraically largest "
                    "(physical) one. Re-run with a larger spectral_shift "
                    "or spectral_shift=\"auto\".",
                    shift, value, value + shift)
            return value, vec, info

        def _project(v):
            return project_fn(v) if project_fn is not None else v

        sigma_old = np.asarray(init_vec)

        # Spectral scale of K, for the seed test and the residual check: the
        # auto shift's own radius estimate, else one deterministic probe
        # ||K g|| / ||g|| (a typical magnitude of K in the iteration's sector).
        rng = np.random.default_rng(_SEED_ADMIX_SEED)
        probe = _project(rng.standard_normal(vec_size)
                         + 1j * rng.standard_normal(vec_size))
        pn = float(np.linalg.norm(probe))
        probe_scale = (float(np.linalg.norm(_project(A_report.matvec(probe))))
                       / pn if pn > 0.0 else 0.0)
        rho_est = None
        if shift_used is not None and isinstance(spectral_shift, str):
            rho_est = max((shift - 1.0e-6) / _AUTO_SHIFT_MARGIN, 0.0)
        # (the probe can exceed the spectral radius by orders of magnitude on
        # a strongly non-normal kernel, so it is used only without the
        # auto-shift estimate)
        spec_scale = rho_est if rho_est else probe_scale
        # no eigenvalue of K can lie beyond the spectral radius; the auto
        # shift sits 5 % above its estimate, so allow a further 5 %
        rho_max = (_AUTO_SHIFT_MARGIN * shift if rho_est is not None
                   else None)

        # The seed test, on K itself: a seed that is an eigenvector of K (or
        # that K annihilates) can never leave its mode -- the iteration then
        # cannot search the seed's symmetry sector for the largest real part.
        # It is NOT perturbed (a perturbation projected only on the channel
        # parity would leave the lattice-symmetry sector of init_gap); the run
        # reports no eigenvalue instead.
        seed_norm = float(np.linalg.norm(sigma_old))
        if seed_norm == 0.0 or not np.isfinite(seed_norm):
            raise ValueError("the initial gap of the power iteration is zero "
                             "or not finite")
        k0 = _project(A_report.matvec(sigma_old))
        k0_norm = float(np.linalg.norm(k0))
        lam0 = complex(np.vdot(sigma_old, k0) / (seed_norm * seed_norm))
        # relative to the spectral scale, not to ||K v||: a seed on a mode of
        # small |lambda| is as much an eigenvector as one on a large mode
        res_scale = (spec_scale * seed_norm if spec_scale > 0.0
                     else max(k0_norm, 1.0e-300))
        seed_resid = float(np.linalg.norm(k0 - lam0 * sigma_old)) / res_scale
        if k0_norm <= _SEED_ANNIHILATED_REL * spec_scale * seed_norm:
            seed_cause = "annihilated"
        elif seed_resid < _SEED_EIGENVECTOR_TOL:
            seed_cause = "eigenvector"
        else:
            seed_cause = None
        if seed_cause is not None:
            what = ("is annihilated by the kernel (K v = 0)"
                    if seed_cause == "annihilated" else
                    "is an eigenvector of the kernel (lambda = {:.8g})"
                    .format(lam0.real))
            logger.warning(
                "Power iteration: the initial gap %s, so the power iteration "
                "cannot search its symmetry sector for the largest real part; "
                "no eigenvalue is reported (lambda = NaN, eigen_converged = "
                "false). Use solver_mode=\"eigenvalue\" or a different "
                "init_gap.", what)
            info = {"converged": False, "diff_converged": False, "n_iter": 0,
                    "two_cycle": False, "seed_is_eigenvector": True,
                    "alpha_applied": alpha_applied,
                    "eigenvalue_validated": False,
                    "eigenvalue_kind": "seed-" + seed_cause,
                    "seed_eigenvalue": (0.0 if seed_cause == "annihilated"
                                        else float(lam0.real)),
                    "rayleigh_eigenvalue": None,
                    "rayleigh_residual": float("nan"),
                    "rayleigh_residual_tol": float("nan"),
                    "shift_sufficient": None}
            if shift_used is not None:
                info["spectral_shift"] = shift_used
            return float("nan"), sigma_old / seed_norm, info

        u_prev = None
        u_prev2 = None
        flips = []          # flip residuals ||u_n + u_{n-1}||
        diffs = []          # iterate differences
        cycling = False     # a sustained sign-alternating two-cycle right now

        for iteration in range(max_iter):
            sigma_new = _project(A.matvec(sigma_old))
            norm = np.linalg.norm(sigma_new)

            # A (possibly projected) iterate can collapse to the zero vector
            # when the kernel annihilates the requested sector; normalizing by
            # `norm` would then produce NaN. Report non-convergence with the
            # last finite iterate instead.
            if norm == 0.0 or not np.isfinite(norm):
                logger.warning(
                    "Eliashberg iterate collapsed to zero norm at iteration "
                    "{}; the kernel annihilates the requested sector, so the "
                    "eigenvalue cannot be normalized. Returning the previous "
                    "iterate as non-converged (lambda = NaN)."
                    .format(iteration))
                return _finish(sigma_old, False, iteration + 1, collapsed=True)

            u = sigma_new / norm
            diff = float(np.linalg.norm(u - sigma_old))
            diffs.append(diff)
            # `norm` is the iterate norm of the (possibly shifted) operator:
            # a progress indicator only, never the reported eigenvalue
            logger.info("Iteration {:4d}: iterate norm = {:.6f}, diff = {:.6e}"
                        .format(iteration, norm - shift, diff))

            # Rate-aware stopping: a small step is not a small distance to the
            # fixed point when the iterate contracts slowly (a seed close to a
            # subdominant mode drifts away from it by tiny steps). The
            # remaining distance is estimated as diff / (1 - rate), with the
            # contraction rate the LARGEST ratio of consecutive differences
            # over the last _RATE_WINDOW steps (clipped below 1), and the run
            # stops only when that estimate is below convergence_tol. A
            # difference at round-off level is accepted outright. The window
            # also replaces the earlier fixed minimum number of steps.
            if diff < convergence_tol and len(diffs) > _RATE_WINDOW:
                recent = diffs[-(_RATE_WINDOW + 1):]
                ratios = [b / a for a, b in zip(recent[:-1], recent[1:])
                          if a > 0.0]
                rate = max(ratios) if ratios else 1.0
                # a difference that is not shrinking is a drift, not a
                # converged iterate: no finite distance estimate
                remaining = (diff / (1.0 - min(rate, _RATE_CLIP))
                             if rate < 1.0 else float("inf"))
                if diff <= _DIFF_ROUNDOFF or remaining < convergence_tol:
                    logger.info("Converged at iteration {} (contraction rate "
                                "{:.4f}, estimated remaining distance {:.2e})"
                                .format(iteration + 1, rate, remaining))
                    return _finish(u, True, iteration + 1)

            # Sign-alternating two-cycle (a NEGATIVE dominant eigenvalue of
            # the iterated operator: only without a sufficient shift). It is
            # recorded for the diagnosis but does NOT stop the run: no finite
            # window can rule out a tiny component of a positive dominant mode
            # that is still growing.
            if u_prev is not None:
                flips.append(float(np.linalg.norm(u + u_prev)))
            cycling = bool(
                u_prev2 is not None
                and len(flips) >= _TWO_CYCLE_WINDOW
                and max(flips[-_TWO_CYCLE_WINDOW:]) < convergence_tol
                and np.linalg.norm(u - u_prev2) < convergence_tol)
            u_prev2 = u_prev
            u_prev = u

            sigma_old = ((1.0 - alpha_applied) * u + alpha_applied * sigma_old
                         if alpha_applied else u)

        if cycling:
            value, vec, info = _finish(u_prev, False, max_iter, two_cycle=True)
            cyc = info.get("cycle_eigenvalue")
            logger.warning(
                "Failed to converge after %d iterations: the iterate is in a "
                "sign-alternating two-cycle (A v = -|A v| v'), i.e. the "
                "dominant (largest-magnitude) eigenvalue of the iterated "
                "operator%s is NEGATIVE%s. The power iteration cannot "
                "converge to the leading pairing eigenvalue here, so no "
                "eigenvalue is reported (lambda = NaN, eigen_converged = "
                "false). Use solver_mode=\"eigenvalue\", or spectral_shift="
                "\"auto\" (the default for solver_mode=\"iteration\")%s.",
                max_iter,
                "" if shift_used is None
                else " K + sigma*I (sigma = {:.6f}: the spectral_shift is "
                     "INSUFFICIENT)".format(shift),
                "" if cyc is None
                else " (Rayleigh quotient on K: {:.6f})".format(cyc),
                "" if shift_used is None else " or a larger spectral_shift")
            return value, vec, info

        logger.warning(
            "Failed to converge after %d iterations: no eigenvalue is "
            "reported (lambda = NaN, eigen_converged = false). %s",
            max_iter, _stall_advice(diffs, convergence_tol, max_iter))
        if shift_used is not None:
            logger.warning(
                "The power iteration ran on the shifted operator K + sigma*I "
                "with sigma = %.6f and still did not converge. Either sigma is "
                "too small (the dominant eigenvalue of K + sigma*I is still "
                "negative -- try a larger spectral_shift, or \"auto\") or it is "
                "so large that it compressed the relative eigenvalue gaps "
                "(try a smaller one, or solver_mode=\"eigenvalue\").", shift)
        return _finish(sigma_old, False, max_iter)

    # Eigenvalue family: ARPACK Arnoldi or shift-invert.
    A, _ = make_operator()

    if not (solver_mode == "arnoldi" or solver_mode.startswith("shift-invert")):
        raise ValueError("Unknown eigenvalue method: {}".format(solver_mode))

    # (spectral_shift was validated at the top of the function, before every
    # branch, so an invalid value fails fast even on the small-dense early
    # return below.)

    # sigma_shift is the shift-invert target; it has no effect on the plain
    # arnoldi path (which uses spectral_shift instead). Warn rather than fail so
    # existing configs keep working.
    if sigma_shift is not None and solver_mode == "arnoldi":
        logger.warning(
            "[eliashberg] sigma_shift is ignored for eigenvalue_method="
            "'arnoldi' (it targets the shift-invert methods); use "
            "spectral_shift to bias the arnoldi selection.")

    if vec_size < 1:
        raise ValueError("Eliashberg operator has empty vector space")

    if vec_size <= 2:
        # scipy.sparse.linalg.eigs requires k < N - 1 for LinearOperator input,
        # so the smallest valid dynamic grid (e.g. norb=1, Nk=1, Nmat=2) cannot
        # go through ARPACK. Reconstruct the tiny dense operator directly.
        dense = np.empty((vec_size, vec_size), dtype=complex)
        basis = np.eye(vec_size, dtype=complex)
        for j in range(vec_size):
            dense[:, j] = A.matvec(basis[:, j])
        try:
            if not np.all(np.isfinite(dense)):
                # checked before LAPACK sees the matrix (it would print an
                # "illegal value" message); the handler below is the backstop
                raise np.linalg.LinAlgError(
                    "the operator has non-finite entries")
            vals, vecs = np.linalg.eig(dense)
        except np.linalg.LinAlgError as exc:
            # e.g. an operator that returns NaN: the NaN sentinel of a
            # failed solve, with a finite placeholder vector
            logger.warning(
                "the %d-dimensional operator cannot be diagonalized (%s): no "
                "leading eigenvalue is reported (NaN, eigen_converged = "
                "false).", vec_size, exc)
            n_keep = min(max(1, num_eigenvalues), vec_size)
            v = np.full(n_keep, np.nan + 0j)
            w = np.eye(vec_size, n_keep, dtype=complex)
            return v[0], w[:, 0], {"eigenvalues": v, "eigenvectors": w,
                                   "sigma_shift": sigma_shift,
                                   "selection": "not-converged"}
        # Match the ARPACK path below: with a seed eigenvector, track the branch
        # that overlaps it (eigenvector continuation); otherwise order by
        # largest real part (the physical SC eigenvalue), not magnitude.
        if seed_vec is not None:
            vals, vecs = _order_by_seed_overlap(vals, vecs, seed_vec)
        else:
            vals, vecs = _order_eigenpairs(vals, vecs)
        n_keep = min(max(1, num_eigenvalues), vec_size)
        vals = vals[:n_keep]
        vecs = vecs[:, :n_keep]
        return vals[0], vecs[:, 0], {"eigenvalues": vals,
                                     "eigenvectors": vecs,
                                     "sigma_shift": sigma_shift,
                                     "selection": "dense-LR"}

    max_ev = min(num_eigenvalues, vec_size - 2)
    if max_ev < 1:
        max_ev = 1

    logger.info("Computing {} eigenvalues with method='{}'...".format(
        max_ev, solver_mode))

    retry_record = {}
    selection_label = None

    def _largest_real_solve(shift_option):
        """``eigs(A + sigma*I, which='LR')`` with sigma subtracted back."""
        from scipy.sparse.linalg import LinearOperator as _LinOp
        if isinstance(shift_option, str):  # "auto"
            from scipy.sparse.linalg import ArpackError, ArpackNoConvergence
            k_pre = min(6, vec_size - 2)
            if k_pre >= 1:
                try:
                    vals_pre, _ = eigs(A, k=k_pre, which='LM',
                                       v0=_arpack_start(A))
                except ArpackNoConvergence as exc:
                    vals_pre = np.asarray(exc.eigenvalues)
                    vals_pre = vals_pre[np.isfinite(vals_pre)]
                    if not vals_pre.size:
                        raise _ShiftProbeFailed(exc) from exc
                    logger.warning(
                        "auto spectral_shift: the spectral-radius probe "
                        "(which='LM') did not fully converge (%s); the %d "
                        "converged eigenvalues set the shift.", exc,
                        vals_pre.size)
                except ArpackError as exc:
                    # any other ARPACK error (e.g. -9 on a null operator):
                    # no radius, handled like a failed solve
                    raise _ShiftProbeFailed(exc) from exc
                if not np.all(np.isfinite(vals_pre)):
                    # a NaN/Inf radius cannot size the shift: a failed
                    # probe, handled like a failed solve (not the overflow
                    # error below, which is for a finite but huge radius)
                    raise _ShiftProbeFailed(_NonFiniteEigenvalues(
                        "the probe returned non-finite eigenvalues"))
                sig = float(np.max(np.abs(vals_pre))) * 1.5 + 1.0e-6
            else:
                sig = 1.0
            if not np.isfinite(sig):
                raise ValueError(
                    "auto spectral_shift overflowed to a non-finite value "
                    "(spectral radius too large); pass an explicit shift.")
        else:
            sig = float(shift_option)
        logger.info("Spectral shift sigma={:.6f}: eigs(A+sigma*I, LR)".format(sig))
        A_sh = _LinOp(A.shape, matvec=lambda v: A.matvec(v) + sig * v,
                      dtype=A.dtype)
        vals_lr, vecs_lr = eigs(A_sh, k=max_ev, which='LR',
                                v0=_arpack_start(A_sh, seed_vec))
        if not np.all(np.isfinite(vals_lr)):
            raise _NonFiniteEigenvalues("non-finite eigenvalues returned")
        return vals_lr - sig, vecs_lr

    def _dense_largest_real():
        """The complete spectrum by dense eigendecomposition (vec_size
        matvecs), the largest real parts first, max_ev of them."""
        dense = np.empty((vec_size, vec_size), dtype=complex)
        basis_col = np.zeros(vec_size, dtype=complex)
        for j in range(vec_size):
            basis_col[j] = 1.0
            dense[:, j] = A.matvec(basis_col)
            basis_col[j] = 0.0
        if not np.all(np.isfinite(dense)):
            # never hand LAPACK a non-finite matrix (it would print an
            # "illegal value" message); the caller's LinAlgError handler
            # returns the NaN sentinel
            raise np.linalg.LinAlgError("the operator has non-finite entries")
        dvals, dvecs = np.linalg.eig(dense)
        dvals, dvecs = _order_eigenpairs(dvals, dvecs)
        return dvals[:max_ev], dvecs[:, :max_ev]

    def _largest_real_or_fallback(shift_option, lm=None):
        """``_largest_real_solve``; when ARPACK does not converge, a dense
        eigendecomposition for vec_size <= _DENSE_FALLBACK_MAX, else NaN
        eigenvalues (the vectors of ``lm`` when given) -- never the
        largest-magnitude set as if it were the largest real part.
        Returns ``(vals, vecs, label_suffix_or_None, record)``."""
        from scipy.sparse.linalg import ArpackError
        try:
            v, w = _largest_real_solve(shift_option)
            return v, w, None, {}
        except (ArpackError, _ShiftProbeFailed, _NonFiniteEigenvalues) as exc:
            what = ("the spectral-radius probe (which='LM') of the auto "
                    "spectral_shift" if isinstance(exc, _ShiftProbeFailed)
                    else "the largest-real-part solve (which='LR')")
            return _failed_solve(exc, what, lm=lm)

    def _failed_solve(exc, what, lm=None):
        """The contract after an ARPACK solve did not converge: a dense
        eigendecomposition (largest real parts first) for vec_size <=
        _DENSE_FALLBACK_MAX, else NaN eigenvalues with finite placeholder
        vectors and the label "not-converged". The one exception is an
        operator that is exactly zero (ARPACK error -9 confirmed by
        _operator_is_exactly_zero): every eigenvalue is 0, reported with
        orthonormal random vectors and the label "null-operator" at any
        dimension.
        Returns ``(vals, vecs, label, record)``."""
        from scipy.sparse.linalg import ArpackError
        if (_is_arpack_zero_start_error(exc)
                and _operator_is_exactly_zero(A, vec_size)):
            logger.info(
                "arnoldi: %s stopped with ARPACK error -9 (starting vector "
                "is zero) and the operator maps %d independent random vectors "
                "to exactly zero: it is the zero operator, every eigenvalue "
                "is 0.", what, _NULL_PROBE_COUNT)
            rng_null = np.random.default_rng(_SEED_ADMIX_SEED)
            w, _ = np.linalg.qr(
                rng_null.standard_normal((vec_size, max_ev))
                + 1j * rng_null.standard_normal((vec_size, max_ev)))
            return (np.zeros(max_ev, dtype=complex), w, "null-operator",
                    {"null_operator": True})
        if vec_size <= _DENSE_FALLBACK_MAX:
            logger.warning(
                "arnoldi: %s did not converge (%s); the operator "
                "(dimension %d <= %d) is diagonalized densely instead.",
                what, exc, vec_size, _DENSE_FALLBACK_MAX)
            try:
                v, w = _dense_largest_real()
                return v, w, "dense-LR", {"lr_retry_failed": True}
            except np.linalg.LinAlgError as dense_exc:
                # e.g. an operator that returns NaN: no eigendecomposition
                # exists, so the NaN sentinel below is the answer
                logger.warning(
                    "arnoldi: the dense diagonalization failed as well (%s).",
                    dense_exc)
                logger.warning(
                    "arnoldi: no leading eigenvalue is reported (NaN, "
                    "eigen_converged = false).")
        else:
            logger.warning(
                "arnoldi: %s did not converge (%s) "
                "and the operator (dimension %d > %d) is too large to "
                "diagonalize densely: no leading eigenvalue is reported (NaN, "
                "eigen_converged = false).%s", what, exc, vec_size,
                _DENSE_FALLBACK_MAX,
                "" if solver_mode.startswith("shift-invert") else
                " Set an explicit spectral_shift or raise num_eigenvalues.")
        # finite placeholder vectors (NOT eigenvectors; the NaN
        # eigenvalues say so): the which='LM' vectors, which every
        # downstream symmetry check accepts (they lie in the operator's
        # space -- for the spin-orbital solvers, the antisymmetric one)
        partial = getattr(exc, "eigenvectors", None)
        if lm is not None:
            w = lm[1]
        elif (partial is not None and np.ndim(partial) == 2
              and partial.shape[0] == vec_size and partial.shape[1] > 0
              and np.all(np.isfinite(partial))):
            # the partially converged vectors of the failed solve
            w = np.asarray(partial)
        else:
            # after a non-finite result another ARPACK pass on the same
            # operator would only repeat it: go straight to random vectors
            cause = exc.args[0] if (isinstance(exc, _ShiftProbeFailed)
                                    and exc.args) else exc
            try:
                if isinstance(cause, (_NonFiniteEigenvalues,
                                      np.linalg.LinAlgError)):
                    raise _NonFiniteEigenvalues("non-finite operator output")
                _, w = eigs(A, k=max_ev, which='LM', v0=_arpack_start(A))
                if not np.all(np.isfinite(w)):
                    raise _NonFiniteEigenvalues("non-finite vectors")
            except (ArpackError, _NonFiniteEigenvalues):
                rng_ph = np.random.default_rng(_SEED_ADMIX_SEED)
                w = (rng_ph.standard_normal((vec_size, max_ev))
                     + 1j * rng_ph.standard_normal((vec_size, max_ev)))
                w /= np.linalg.norm(w, axis=0)
        v = np.full(w.shape[1], np.nan + 0j)
        return v, w, "not-converged", {"lr_retry_failed": True}

    if solver_mode == "arnoldi":
        if spectral_shift is not None:
            # Select the LARGEST-REAL eigenvalue (the physical SC eigenvalue,
            # lambda -> 1 at Tc), which plain which='LM' misses when a small
            # positive lambda is masked by larger repulsive (negative)
            # eigenvalues. We ask ARPACK for the largest real part directly
            # (which='LR') -- unlike which='LM', this is the correct criterion
            # even for the non-Hermitian kernel's complex eigenvalues, where
            # |lambda+sigma| would otherwise be dominated by a large-|Im| mode.
            # The real spectral shift A -> A + sigma*I preserves the real-part
            # ordering (Re(lambda+sigma) = Re(lambda) + sigma) but moves the
            # spectrum into the right half-plane, which conditions ARPACK's LR
            # iteration; we subtract sigma afterwards.
            # (spectral_shift is validated near the top of _solve_leading.)
            vals, vecs, fallback, retry_record = _largest_real_or_fallback(
                spectral_shift)
            selection_label = fallback or "LR"
        else:
            from scipy.sparse.linalg import ArpackError
            try:
                vals, vecs = eigs(A, k=max_ev, which='LM',
                                  v0=_arpack_start(A, seed_vec))
                if not np.all(np.isfinite(vals)):
                    raise _NonFiniteEigenvalues(
                        "non-finite eigenvalues returned")
                selection_label = "LM"
            except (ArpackError, _NonFiniteEigenvalues) as exc:
                # the unshifted solve itself did not converge: the same
                # contract as a failed largest-real-part solve (dense
                # fallback, else NaN "not-converged"), never an exception
                vals, vecs, selection_label, _rec = _failed_solve(
                    exc, "the unshifted solve (which='LM')")
                retry_record = dict({"lm_solve_failed": True},
                                    **_null_flag(_rec))
            # The LM set holds every eigenvalue of magnitude above the
            # smallest one it contains, m = min |lambda_i|; any eigenvalue
            # outside it has real part <= m. So the largest real part FOUND
            # is the largest real part of the whole spectrum whenever it is
            # >= m. Otherwise (e.g. a repulsive-dominant kernel whose few
            # largest-magnitude modes are all negative) the physical,
            # largest-real eigenvalue may be missing: re-solve for the
            # largest real part. A seeded run tracks its continuation branch
            # on purpose and is left alone.
            if selection_label == "LM" and seed_vec is None and len(vals):
                re_max = float(np.max(vals.real))
                m_min = float(np.min(np.abs(vals)))
                tol_c = 1.0e-10 * max(1.0, float(np.max(np.abs(vals))))
                if re_max < m_min - tol_c:
                    logger.warning(
                        "arnoldi (which='LM', %d eigenvalues): the largest "
                        "real part found, %.6g, is below the smallest "
                        "magnitude found, %.6g, so an eigenvalue with a "
                        "larger real part may lie outside the computed set; "
                        "re-solving for the largest real part "
                        "(spectral_shift=\"auto\"). Set [eliashberg] "
                        "spectral_shift=\"auto\" to skip the first solve.",
                        max_ev, re_max, m_min)
                    vals, vecs, fallback, rec = _largest_real_or_fallback(
                        "auto", lm=(vals, vecs))
                    selection_label = fallback or "LR_retry"
                    retry_record = dict({"lm_max_real": re_max,
                                         "lm_min_abs": m_min}, **rec)

    elif solver_mode.startswith("shift-invert"):
        if sigma_shift is None:
            # Estimate shift from a quick Arnoldi run. Sample a few
            # largest-magnitude eigenvalues and aim at the largest *real* part
            # (the physical SC eigenvalue), not the largest magnitude (which
            # can be a large negative repulsive mode).
            k_pre = min(6, vec_size - 2)
            if k_pre < 1:
                # Operator too small for a preliminary ARPACK pass (ARPACK
                # needs k < N-1); fall back to a neutral shift.
                sigma_shift = 0.0
            else:
                logger.info("Estimating shift with preliminary Arnoldi...")
                from scipy.sparse.linalg import ArpackError
                try:
                    vals_pre, _ = eigs(A, k=k_pre, which='LM',
                                       v0=_arpack_start(A))
                    if not np.all(np.isfinite(vals_pre)):
                        raise _NonFiniteEigenvalues(
                            "non-finite eigenvalues returned")
                    sigma_shift = _shift_from_eigenvalues(vals_pre)
                except (ArpackError, _NonFiniteEigenvalues) as exc:
                    # no target for the shift-invert solve: the same
                    # contract as a failed solve (dense fallback, else NaN)
                    vals, vecs, selection_label, _rec = _failed_solve(
                        exc, "the preliminary which='LM' solve that sets "
                             "the shift-invert target")
                    retry_record = dict({"lm_solve_failed": True},
                                        **_null_flag(_rec))
                    sigma_shift = None
        if selection_label is None:
            logger.info("Using sigma_shift = {:.6f}".format(sigma_shift))
            from scipy.sparse.linalg import ArpackError
            try:
                vals, vecs = _eigs_shift_invert(
                    A, vec_size, max_ev, solver_mode, sigma=sigma_shift,
                    seed_vec=seed_vec
                )
                if not np.all(np.isfinite(vals)):
                    # e.g. a singular K - sigma*I (a null kernel with a zero
                    # target): no eigenvalue, the failed-solve contract
                    raise _NonFiniteEigenvalues(
                        "non-finite eigenvalues returned")
                selection_label = "shift-invert"
            except (ArpackError, _NonFiniteEigenvalues) as exc:
                vals, vecs, selection_label, _rec = _failed_solve(
                    exc, "the shift-invert solve")
                retry_record = dict({"shift_invert_failed": True},
                                    **_null_flag(_rec))

    # With a seed eigenvector, track the branch that overlaps it (eigenvector
    # continuation); otherwise order by largest real part (the physical SC
    # eigenvalue), not magnitude.
    if seed_vec is not None:
        vals, vecs = _order_by_seed_overlap(vals, vecs, seed_vec)
    else:
        vals, vecs = _order_eigenpairs(vals, vecs)

    # A negative leading eigenvalue from plain which='LM' (no spectral_shift) is
    # often an artifact: a small positive lambda masked by larger repulsive
    # modes. Tip the user toward spectral_shift='auto'. Skip this when a
    # seed_vec is given -- there vals[0] is the seed-overlap continuation
    # branch, which may be intentionally negative, not a masked leading mode.
    # Require a meaningfully negative value (relative to the spectral scale)
    # so roundoff-scale negatives near a numerically-zero leading eigenvalue
    # do not trigger a misleading recommendation.
    if (warn_negative_leading and solver_mode == "arnoldi"
            and spectral_shift is None and seed_vec is None and len(vals)):
        scale = float(np.max(np.abs(vals))) if len(vals) else 0.0
        neg_tol = 1.0e-8 * max(scale, 1.0)
        if vals[0].real < -neg_tol:
            logger.warning(
                "Leading eigenvalue Re(lambda)=%.4g is negative; if a positive "
                "(attractive) mode is expected, set [eliashberg] spectral_shift="
                "\"auto\" so the largest-REAL eigenvalue is selected instead of "
                "the largest-magnitude one.", vals[0].real)

    return vals[0], vecs[:, 0], dict({"eigenvalues": vals, "eigenvectors": vecs,
                                      "sigma_shift": sigma_shift,
                                      "selection": selection_label},
                                     **retry_record)


def _solve_eigenvalue(Vs_q, G2, norb, Nx, Ny, Nz, num_eigenvalues=10,
                      method="arnoldi", sigma_shift=None, spectral_shift=None,
                      operator=None):
    """Solve linearized Eliashberg equation by eigenvalue analysis.

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
    num_eigenvalues : int
        Number of eigenvalues to compute.
    method : str
        Eigenvalue solver method:
        - "arnoldi" : Implicitly Restarted Arnoldi (ARPACK eigs).
            Standard Krylov subspace method for non-symmetric matrices.
            Directly finds the largest eigenvalues.
        - "shift-invert-bicgstab" : Shift-invert with BiCGSTAB.
            Solves (K - sigma*I)^{-1} eigenvalue problem using BiCGSTAB
            for the linear system. Efficient for finding eigenvalues near
            a target value sigma.
        - "shift-invert-gmres" : Shift-invert with GMRES.
            Same as above but uses GMRES for the linear system.
            More robust than BiCGSTAB for non-symmetric systems.
        - "shift-invert-lgmres" : Shift-invert with LGMRES.
            Uses LGMRES (loose GMRES) which can be faster for some problems.
        - "subspace" : Subspace iteration (block power method).
            Simultaneously propagates multiple vectors through the kernel
            with QR orthogonalization. Robust for degenerate eigenvalues
            and finds multiple eigenvalues with different symmetries.
    sigma_shift : float, optional
        Shift parameter for shift-invert methods. Eigenvalues near this
        value are found most efficiently. If None, a preliminary Arnoldi
        step estimates an appropriate shift value.
    operator : tuple, optional
        Pre-built ``(A, vec_size)`` kernel, bypassing
        ``_make_kernel_operator(Vs_q, G2, ...)`` (see ``_solve_iteration``).

    Returns
    -------
    eigenvalues : ndarray
        Leading eigenvalues ordered by descending real part (largest first).
    eigenvectors : ndarray
        Corresponding eigenvectors reshaped to (num_ev, norb, norb, Nx, Ny, Nz).
    """
    vec_size = norb * norb * Nx * Ny * Nz

    # Safety check (review fix): when a pre-built ``operator`` is supplied,
    # its own vec_size (``operator[1]``) is otherwise never read here -- every
    # downstream call (``_solve_leading``, ``_solve_subspace_iteration``) is
    # handed the LOCALLY recomputed ``norb*norb*Nx*Ny*Nz`` instead. The two
    # agree today for both callers of this path (the standard
    # ``_make_kernel_operator`` result and the bond-channel
    # ``_build_bond_operator`` result both use the plain
    # orbital-pair x spatial vec_size), but nothing enforces that; a future
    # dynamic/enlarged bond operator with a DIFFERENT internal size would
    # silently reshape eigenvectors into the wrong shape instead of failing
    # loudly. Assert the invariant instead of trusting it silently.
    if operator is not None:
        assert operator[1] == vec_size, (
            "_solve_eigenvalue: the supplied operator's vec_size ({}) does "
            "not match norb**2*Nx*Ny*Nz ({}); eigenvectors would be reshaped "
            "with the wrong size.".format(operator[1], vec_size))

    spectral_shift = _normalize_spectral_shift(spectral_shift)
    if method == "subspace":
        # spectral_shift routes through the ARPACK largest-real path in
        # _solve_leading; the subspace driver never reaches it, so reject a
        # non-None spectral_shift here (same arnoldi-only contract) rather
        # than silently ignoring a misconfigured static input.
        if spectral_shift is not None:
            raise ValueError(
                "[eliashberg] spectral_shift is only supported for "
                "eigenvalue_method='arnoldi', not '{}'".format(method))
        # Subspace (block power) iteration has its own dedicated driver
        # (magnitude-based Ritz selection, not the ARPACK/shift-invert path),
        # so it is not routed through _solve_leading; call it directly, exactly
        # as before.
        max_ev = min(num_eigenvalues, vec_size - 2)
        if max_ev < 1:
            max_ev = 1
        logger.info("Computing {} eigenvalues with method='{}'...".format(
            max_ev, method))
        return _solve_subspace_iteration(
            Vs_q, G2, norb, Nx, Ny, Nz,
            num_eigenvalues=max_ev, operator=operator
        )

    # ARPACK Arnoldi / shift-invert: delegate the eigen-selection, shift
    # estimation, and descending-real-part ordering to the shared driver.
    make_operator = (lambda: _make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)) \
        if operator is None else (lambda: operator)
    _, _, eig_analysis = _solve_leading(
        make_operator, vec_size, method,
        num_eigenvalues=num_eigenvalues, sigma_shift=sigma_shift,
        spectral_shift=spectral_shift,
    )
    vals = eig_analysis["eigenvalues"]
    vecs = eig_analysis["eigenvectors"]

    eigenvectors = np.array([
        vecs[:, i].reshape(norb, norb, Nx, Ny, Nz)
        for i in range(len(vals))
    ])

    return vals, eigenvectors


def _eigs_shift_invert(A, vec_size, num_ev, method, sigma=0.0, rtol_linear=1e-8,
                       seed_vec=None):
    """Eigenvalue computation using shift-invert with iterative linear solver.

    Transforms the eigenvalue problem K*x = lambda*x into
    (K - sigma*I)^{-1}*x = nu*x where nu = 1/(lambda - sigma).
    Eigenvalues of the original problem near sigma become the largest
    eigenvalues of the transformed problem, which Arnoldi finds efficiently.

    The inverse (K - sigma*I)^{-1} is never formed explicitly; instead,
    each application solves (K - sigma*I)*y = x using an iterative solver.

    Parameters
    ----------
    A : LinearOperator
        The Eliashberg kernel operator K.
    vec_size : int
        Dimension of the vector space.
    num_ev : int
        Number of eigenvalues to compute.
    method : str
        One of "shift-invert-bicgstab", "shift-invert-gmres",
        "shift-invert-lgmres".
    sigma : float
        Shift value. Eigenvalues near sigma are found most efficiently.
        Default 0.0 finds eigenvalues nearest to zero.
    rtol_linear : float
        Relative tolerance for the iterative linear solver.

    Returns
    -------
    eigenvalues : ndarray
        Eigenvalues of the original problem K.
    eigenvectors : ndarray
        Corresponding eigenvectors.
    """
    solver_name = method.replace("shift-invert-", "")
    solver_map = {
        "bicgstab": bicgstab,
        "gmres": gmres,
        "lgmres": lgmres,
    }
    if solver_name not in solver_map:
        raise ValueError("Unknown linear solver: {}".format(solver_name))
    linear_solver = solver_map[solver_name]

    # (A - sigma*I)
    def shifted_matvec(v):
        return A.matvec(v) - sigma * v

    A_shifted = LinearOperator((vec_size, vec_size),
                               matvec=shifted_matvec, dtype=complex)

    solve_count = [0]
    fail_count = [0]

    # (A - sigma*I)^{-1} via iterative solver
    def inv_matvec(v):
        solve_count[0] += 1
        x, info = linear_solver(A_shifted, v, rtol=rtol_linear, maxiter=500)
        if info != 0:
            fail_count[0] += 1
            if fail_count[0] <= 3:
                logger.warning(
                    "Linear solver {} did not converge (info={}), "
                    "solve #{}, using approximate solution".format(
                        solver_name, info, solve_count[0]))
        return x

    A_inv = LinearOperator((vec_size, vec_size),
                           matvec=inv_matvec, dtype=complex)

    # eigs on (A - sigma*I)^{-1} finds eigenvalues nu = 1/(lambda - sigma)
    # largest |nu| correspond to lambda closest to sigma. A seed vector (v0)
    # biases the Arnoldi start toward the physical branch being tracked.
    nus, vecs = eigs(A_inv, k=num_ev, which='LM', v0=_arpack_start(A_inv, seed_vec))

    logger.info("Shift-invert: {} linear solves, {} failures".format(
        solve_count[0], fail_count[0]))

    # Convert back: lambda = 1/nu + sigma
    eigenvalues = 1.0 / nus + sigma

    return eigenvalues, vecs


def _solve_subspace_iteration(Vs_q, G2, norb, Nx, Ny, Nz,
                              num_eigenvalues=5, max_iter=300, tol=1e-6,
                              operator=None):
    """Find multiple eigenvalues by subspace iteration (block power method).

    Simultaneously propagates a block of vectors through the kernel,
    orthogonalizing via QR decomposition at each step. This naturally
    converges to the dominant invariant subspace.

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
    num_eigenvalues : int
        Number of eigenvalues to compute.
    max_iter : int
        Maximum iterations.
    tol : float
        Convergence tolerance for eigenvalues.
    operator : tuple, optional
        Pre-built ``(A, vec_size)`` kernel, bypassing
        ``_make_kernel_operator(Vs_q, G2, ...)`` (see ``_solve_iteration``).

    Returns
    -------
    eigenvalues : ndarray
        Converged eigenvalues ordered by descending magnitude.
    eigenvectors : ndarray
        Shape (num_eigenvalues, norb, norb, Nx, Ny, Nz).
    """
    A, vec_size = (_make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)
                   if operator is None else operator)
    num_ev = min(num_eigenvalues, vec_size)

    # Use extra vectors for better convergence (subspace padding)
    n_work = min(num_ev + max(num_ev, 5), vec_size)

    # Random initial subspace
    rng = np.random.default_rng(42)
    V = rng.standard_normal((vec_size, n_work))
    V, _ = np.linalg.qr(V)

    eigenvalues_old = np.zeros(num_ev, dtype=complex)

    for iteration in range(max_iter):
        # Apply kernel to all vectors at once: W = A @ V (kernel may be
        # complex). matmat batches all columns into a single FFT (numerically
        # identical to the per-column matvec).
        W = A.matmat(V)

        # Rayleigh quotient: H = V^dagger A V (conjugate transpose for complex)
        H = V.conj().T @ W

        # Eigendecomposition of small matrix
        evals_h, evecs_h = np.linalg.eig(H)
        idx = np.argsort(-np.abs(evals_h))
        evals_h = evals_h[idx]
        evecs_h = evecs_h[:, idx]

        # Ritz vectors: update subspace
        V = W @ evecs_h
        # Re-orthogonalize
        V, _ = np.linalg.qr(V)

        # Check convergence of the wanted eigenvalues (track the full complex
        # value so imaginary motion is not ignored for a complex kernel).
        eigenvalues_new = evals_h[:num_ev]
        diff = np.max(np.abs(eigenvalues_new - eigenvalues_old))

        if iteration % 10 == 0 or diff < tol:
            logger.info("Subspace iter {:4d}: eigenvalues = {}, diff = {:.2e}".format(
                iteration, np.array2string(eigenvalues_new, precision=4), diff))

        if diff < tol:
            logger.info("Subspace iteration converged at iteration {}".format(
                iteration + 1))
            break

        eigenvalues_old = eigenvalues_new.copy()

    # Extract final Ritz vectors for the wanted eigenvalues
    # Recompute from final V (batched matmat = per-column matvec)
    W = A.matmat(V)
    H = V.conj().T @ W
    evals_h, evecs_h = np.linalg.eig(H)
    # Subspace (block power) iteration is magnitude-based: it converges to the
    # dominant-magnitude invariant subspace, so report its modes by magnitude.
    # The physical SC selection (largest real) is applied in the default
    # arnoldi path; use that to read off the leading SC eigenvalue.
    idx = np.argsort(-np.abs(evals_h))

    eigenvalues = evals_h[idx[:num_ev]]
    ritz_vecs = V @ evecs_h[:, idx[:num_ev]]

    eigenvectors = np.array([
        ritz_vecs[:, i].reshape(norb, norb, Nx, Ny, Nz)
        for i in range(num_ev)
    ])

    return eigenvalues, eigenvectors


def _solve_shifted_bicg(Vs_q, G2, norb, Nx, Ny, Nz,
                        sigma_list, num_eigenvalues=3, tol_linear=1e-8):
    """Find eigenvalues near multiple target values using shifted BiCG.

    Note: this helper does not accept an ``operator=`` (bond-channel)
    argument -- it always builds its own kernel from ``Vs_q``/``G2`` via
    ``_make_kernel_operator``. It is unreachable from ``calc_eliashberg``
    when ``bond_channels=true`` (that path sets ``Vs_q = G2 = None`` and
    drives the pre-built bond operator through ``_solve_iteration``/
    ``_solve_eigenvalue`` instead; nothing in the bond dispatch calls this
    function).

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
    sigma_list : list of float
        Target shift values. Eigenvalues near each sigma are found.
    num_eigenvalues : int
        Number of eigenvalues to find near each sigma.
    tol_linear : float
        Tolerance for the linear solver.

    Returns
    -------
    all_eigenvalues : dict
        Dictionary mapping sigma -> array of eigenvalues.
    all_eigenvectors : dict
        Dictionary mapping sigma -> array of eigenvectors.
    """
    A, vec_size = _make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)

    max_ev = min(num_eigenvalues, vec_size - 2)
    if max_ev < 1:
        max_ev = 1

    all_eigenvalues = {}
    all_eigenvectors = {}

    # Use the first sigma as seed system
    sigma_seed = sigma_list[0]

    # Seed: build Krylov subspace for (A - sigma_seed I)
    # We use a manual BiCG implementation with shift tracking
    logger.info("Shifted BiCG: seed sigma = {:.6f}, {} additional shifts".format(
        sigma_seed, len(sigma_list) - 1))

    # For each shift, solve via shift-invert + eigs
    # The key optimization: reuse preconditioner from seed system
    for i, sigma in enumerate(sigma_list):
        logger.info("Processing shift sigma = {:.6f} ({}/{})".format(
            sigma, i + 1, len(sigma_list)))

        def shifted_matvec(v, s=sigma):
            return A.matvec(v) - s * v

        A_shifted = LinearOperator(
            (vec_size, vec_size), matvec=shifted_matvec, dtype=complex
        )

        def inv_matvec(v):
            x, info = gmres(A_shifted, v, rtol=tol_linear, maxiter=300)
            return x

        A_inv = LinearOperator(
            (vec_size, vec_size), matvec=inv_matvec, dtype=complex
        )

        try:
            nus, vecs = eigs(A_inv, k=max_ev, which='LM', v0=_arpack_start(A_inv))
            eigenvalues = 1.0 / nus + sigma
            # Eigenvalues found near this shift; report by magnitude.
            idx = np.argsort(-np.abs(eigenvalues))
            eigenvalues = eigenvalues[idx]
            eigvecs = np.array([
                vecs[:, j].reshape(norb, norb, Nx, Ny, Nz)
                for j in idx
            ])
        except Exception as e:
            logger.warning("Shift sigma={:.6f} failed: {}".format(sigma, e))
            eigenvalues = np.array([])
            eigvecs = np.array([])

        all_eigenvalues[sigma] = eigenvalues
        all_eigenvectors[sigma] = eigvecs

        if len(eigenvalues) > 0:
            logger.info("  Found eigenvalues: {}".format(
                np.array2string(eigenvalues.real, precision=6)))

    return all_eigenvalues, all_eigenvectors


class _SOAntisymmetricBasis:
    """Orthonormal basis of the fermion-antisymmetric spin-orbital gaps
    ``Delta_AB(k) = -Delta_BA(-k)`` on ``shape = (nd, nd, Nx, Ny, Nz)``.

    The map ``R: (A, B, k) -> (B, A, -k)`` is an involution of the flat gap
    index. Each orbit ``{i, R i}`` with ``i != R i`` carries one coordinate,
    the unit vector ``(e_i - e_{R i}) / sqrt(2)`` (``i`` the smaller index of
    the orbit); a self-paired element ``i = R i`` (``A = B`` and ``k = -k``)
    is forced to zero by the antisymmetry and carries none. ``E``
    (:meth:`to_full`) is therefore an isometry onto range(P), ``P = (1 -
    R) / 2``, with ``E^dagger E = 1`` and ``E E^dagger = P``; ``dim`` is
    ``(n - n_self) / 2``. The spin-orbital eigenproblem is solved on
    ``E^dagger K E`` (:meth:`reduce`), where no Pauli-forbidden mode exists,
    so no shift or selection can return one.

    ``reverse`` is the involution ``R`` as a map on index arrays of
    ``shape``; the default is the static ``(A, B, k) -> (B, A, -k)``
    (:func:`_reverse_k_and_orbital`). The dynamic solver (issue #284) passes
    ``eliashberg_dynamic._reverse_kw_and_orbital`` on ``(nd, nd, Nx, Ny, Nz,
    nmat)``, ``(A, B, k, n) -> (B, A, -k, nmat - 1 - n)``: a fermionic
    frequency is never zero, so that map has no self-paired element."""

    def __init__(self, shape, reverse=None):
        self.shape = tuple(int(x) for x in shape)
        n = int(np.prod(self.shape))
        idx = np.arange(n).reshape(self.shape)
        rev = _reverse_k_and_orbital if reverse is None else reverse
        self.partner = np.ascontiguousarray(rev(idx)).ravel()
        if self.partner.shape != (n,) or not np.array_equal(
                self.partner[self.partner], np.arange(n)):
            raise ValueError("the reverse map of the antisymmetric basis must "
                             "be an involution of the flat gap index")
        flat = np.arange(n)
        self.rep = flat[flat < self.partner]
        self.rep_partner = self.partner[self.rep]
        self.n_self = int(np.sum(flat == self.partner))
        self.n = n
        self.dim = len(self.rep)

    #: the full flat gap indices a subset basis acts on (:meth:`from_subset`);
    #: ``None`` for the basis of the whole gap
    indices = None

    @classmethod
    def from_subset(cls, shape, indices, reverse=None):
        """The basis restricted to the flat gap indices ``indices`` of
        ``shape`` (strictly ascending, closed under the involution ``R``).

        The sector vector is ``v = g[indices]`` of a flat gap ``g``;
        :meth:`to_full`, :meth:`to_reduced` and :meth:`reduce` act on sector
        vectors (length ``n = len(indices)``) and :meth:`scatter` places one
        in the full flat gap. Because the subset is closed, every orbit
        ``{i, R i}`` lies inside or outside it, and the coordinates are those
        of the full basis on the orbits inside (same representative, same
        sign). A subset that is not closed is refused (``ValueError`` naming
        the first index whose partner is missing).

        ``shape`` stays the FULL gap shape, so a sector vector cannot be
        reshaped to it: with a subset basis call
        :func:`_so_solve_eigenvalue` with ``map_back = False`` and place each
        vector with ``scatter(to_full(c))`` (never :func:`_so_solve_iteration`,
        which maps back by reshaping)."""
        shape = tuple(int(x) for x in shape)
        n_full = int(np.prod(shape))
        idx = np.asarray(indices)
        if idx.ndim != 1 or not (idx.size == 0 or np.issubdtype(idx.dtype, np.integer)):
            raise ValueError("the subset of the antisymmetric basis must be a one-dimensional "
                             "integer index array")
        idx = idx.astype(np.int64)
        if idx.size and (idx[0] < 0 or idx[-1] >= n_full or np.any(np.diff(idx) <= 0)):
            raise ValueError("the subset of the antisymmetric basis must be strictly ascending "
                             "flat gap indices in [0, {})".format(n_full))
        rev = _reverse_k_and_orbital if reverse is None else reverse
        partner = np.ascontiguousarray(
            rev(np.arange(n_full).reshape(shape))).ravel()
        if partner.shape != (n_full,) or not np.array_equal(
                partner[partner], np.arange(n_full)):
            raise ValueError("the reverse map of the antisymmetric basis must "
                             "be an involution of the flat gap index")
        pos = np.full(n_full, -1, dtype=np.int64)
        pos[idx] = np.arange(idx.size)
        local = pos[partner[idx]]
        missing = np.flatnonzero(local < 0)
        if missing.size:
            i = int(idx[missing[0]])
            raise ValueError("the subset of the antisymmetric basis is not closed under the "
                             "involution (index {} maps to {})".format(i, int(partner[i])))
        del pos, partner
        self = cls.__new__(cls)
        self.shape = shape
        idx.flags.writeable = False
        self.indices = idx
        self.partner = local
        flat = np.arange(idx.size)
        self.rep = flat[flat < local]
        self.rep_partner = local[self.rep]
        self.n_self = int(np.sum(flat == local))
        self.n = int(idx.size)
        self.dim = len(self.rep)
        return self

    def scatter(self, v):
        """A sector vector ``(n,)`` / ``(n, k)`` placed in the full flat gap
        (zeros elsewhere; a copy for the basis of the whole gap)."""
        v = np.asarray(v)
        if self.indices is None:
            return np.array(v, dtype=complex)
        out = np.zeros((int(np.prod(self.shape)),) + v.shape[1:], dtype=complex)
        out[self.indices] = v
        return out

    def to_full(self, c):
        """``E c``: reduced coordinates ``(dim,)`` or ``(dim, k)`` -> flat
        antisymmetric gaps ``(n,)`` / ``(n, k)``."""
        c = np.asarray(c)
        out = np.zeros((self.n,) + c.shape[1:], dtype=complex)
        h = c / np.sqrt(2.0)
        out[self.rep] = h
        out[self.rep_partner] = -h
        return out

    def to_reduced(self, v):
        """``E^dagger v``: flat gaps ``(n,)`` / ``(n, k)`` -> reduced
        coordinates."""
        v = np.asarray(v)
        return (v[self.rep] - v[self.rep_partner]) / np.sqrt(2.0)

    def reduce(self, A, *, adjoint=False):
        """``E^dagger K E`` of the flat kernel operator ``A`` (``E^dagger P
        = E^dagger``, so this is also ``E^dagger (P K) E``).

        ``adjoint=True`` also attaches ``rmatvec = E^dagger K^dagger E``
        (``A`` must provide ``rmatvec``); by default the operator has the
        forward product only."""
        def matvec(c):
            v = self.to_full(np.asarray(c).ravel())
            return self.to_reduced(np.asarray(A.matvec(v)).ravel())

        def matmat(C):
            C = np.asarray(C)
            if C.ndim == 1:
                return matvec(C)
            return self.to_reduced(np.asarray(A.matmat(self.to_full(C))))
        if not adjoint:
            return LinearOperator((self.dim, self.dim), matvec=matvec,
                                  matmat=matmat, dtype=complex)

        def rmatvec(c):
            v = self.to_full(np.asarray(c).ravel())
            return self.to_reduced(np.asarray(A.rmatvec(v)).ravel())
        return LinearOperator((self.dim, self.dim), matvec=matvec,
                              matmat=matmat, rmatvec=rmatvec, dtype=complex)
