"""Dynamic Eliashberg equation solver (frequency-resolved).

This module implements the full-frequency Eliashberg equation solver for
analyzing superconducting instabilities with frequency-dependent kernels.
"""

import logging
import os

import numpy as np

from hwave.solver import backend
from hwave.solver import matsubara as ms
from hwave.solver.kgrid import reverse_fft_axes
# Shared spatial-FFT helpers (scipy-parallel on CPU, cuFFT on GPU) live in
# backend.py so RPA/FLEX use the same implementations; keep the module-local
# names used throughout this file and by the tests.
from hwave.solver.backend import _SFFT  # noqa: F401
from hwave.solver.backend import spatial_fftn as _spatial_fftn
from hwave.solver.backend import spatial_ifftn as _spatial_ifftn

logger = logging.getLogger("qlms").getChild("eliashberg_dynamic")


def _gpu_requested(eli_param):
    """Whether ``[eliashberg] gpu`` requests GPU execution.

    Accepts a real TOML boolean, and (for programmatic dict configs) a truthy
    string such as ``"true"``/``"1"``/``"yes"`` so a stray ``"false"`` string is
    not treated as True by ``bool("false")``. Delegates to the shared
    ``backend.as_bool`` (also used by the RPA/FLEX/UHF gpu-flag reads).
    """
    return backend.as_bool(eli_param.get("gpu", False))


def _gpu_required_requested(eli_param):
    """Whether ``[eliashberg] gpu_required`` requests strict GPU mode: fail
    fast (get_backend raises) instead of falling back to CPU when CuPy/CUDA is
    unusable. Like ``_gpu_requested``, coerced via ``backend.as_bool`` so a
    programmatic string "false"/"0" is not read as True. Default false.
    """
    return backend.as_bool(eli_param.get("gpu_required", False))


def _ir_keep_static_requested(eli_param):
    """Whether ``[eliashberg] ir_keep_static_chi`` requests retaining the
    IR-compression static (frequency-independent) susceptibility component.

    Like ``_gpu_requested``, delegates to ``backend.as_bool`` so a programmatic
    dict config with a string ``"false"``/``"off"``/``"0"`` is not read as True
    by plain ``bool("false")`` -- which would silently retain the static
    component and change the numerical path. TOML booleans are unaffected.
    """
    return backend.as_bool(eli_param.get("ir_keep_static_chi", False))


# ---------------------------------------------------------------------------
# Channel sectors of the frequency-resolved gap
#
# The static path filters eigenpairs by the (k, orbital) parity operator
# Delta_{ab}(k) -> Delta_{ba}(-k) (sc._reverse_k_and_orbital), so a singlet
# request never reports an odd-parity (triplet-sector) mode. The dynamic gap
# carries an extra fermionic Matsubara axis, and fermion antisymmetry makes the
# COMBINED operator
#     P: Delta_{ab}(k, iw_n) -> Delta_{ba}(-k, -iw_n)
# the antisymmetry CONSTRAINT: even (+) for singlet, odd (-) for triplet.
# That constraint alone is not a channel, because each of its eigenspaces
# still contains two frequency parities:
#     singlet (P = +1): (even k, even w) OR (odd k, odd w)
#     triplet (P = -1): (odd k, even w)  OR (even k, odd w)
# The conventional channels reported in the literature are the EVEN-FREQUENCY
# sectors, so the projectors below select
#     singlet -> (even k, even w),   triplet -> (odd k, even w).
#
# P factorises into two commuting involutions,
#     Pk: Delta_{ab}(k, iw_n) -> Delta_{ba}(-k, iw_n)   (momentum parity,
#                                                        carrying the orbital
#                                                        transpose, as on the
#                                                        static path)
#     Pw: Delta_{ab}(k, iw_n) -> Delta_{ab}(k, -iw_n)   (frequency parity)
# with P = Pk Pw. EVEN FREQUENCY is invariance under Pw alone: an even-frequency
# singlet obeys Delta_{ab}(k, iw) = Delta_{ab}(k, -iw), and together with the
# antisymmetry constraint that gives Delta_{ba}(-k, iw) = Delta_{ab}(k, iw),
# i.e. even momentum parity as well.
#
# Odd-frequency pairing (the other two sectors) is NOT solved for by this
# module; a returned gap's weight in each of the four sectors is recorded by
# ``gap_sector_weights`` so an impure result is visible in the outputs
# (issue #209).
#
# Neither half is by itself a symmetry of a general multi-orbital pairing
# kernel -- only the product P is (measured: a two-orbital FLEX solve leaks
# 0.46 under Pk and 0.54 under Pw while P is conserved to 3e-16). The channel
# projection is therefore applied to the power iteration only where
# ``_channel_leakage`` says the kernel preserves the channel sector; otherwise
# the iteration falls back to the combined-parity sector, so the reported
# eigenvalue is always an eigenvalue of the kernel and never of a restriction
# of it. Which sector a run selected in is recorded as ``sector_selection``.
# ---------------------------------------------------------------------------

def _reverse_kw_and_orbital(gap_w):
    """Return Delta_{ba}(-k, -iw_n) for a dynamic gap.

    Parameters
    ----------
    gap_w : ndarray
        Gap, shape ``(norb, norb, Nx, Ny, Nz, nmat)``.

    Returns
    -------
    ndarray
        ``Delta_{ba}(-k, -iw_n)``, same shape. Spatial ``k -> -k`` is the
        shared FFT-grid reversal ``i -> (N - i) % N``
        (``kgrid.reverse_fft_axes``, matching
        ``sc._reverse_k_and_orbital``); the centered fermionic Matsubara
        partner of index ``n`` is ``nmat - 1 - n`` (a plain reversal, no roll);
        the two orbital indices are swapped.
    """
    rev = reverse_fft_axes(gap_w, (2, 3, 4))  # k -> -k on the FFT grid
    rev = rev[..., ::-1]                     # iw_n -> -iw_n (centered fermionic)
    rev = np.swapaxes(rev, 0, 1)             # orbital transpose a <-> b
    return rev


def _reverse_k_and_orbital_dynamic(gap_w):
    """Return ``Delta_{ba}(-k, iw_n)``: the momentum half ``Pk`` of the
    combined parity operator, frequency untouched.

    The same ``k -> -k`` FFT-grid reversal ``i -> (N - i) % N`` and orbital
    transpose as ``sc._reverse_k_and_orbital`` on the static path: the momentum
    parity of a pair amplitude carries the exchange of the two orbitals.
    """
    return np.swapaxes(reverse_fft_axes(gap_w, (2, 3, 4)), 0, 1)


def _reverse_w_dynamic(gap_w):
    """Return ``Delta_{ab}(k, -iw_n)``: the frequency half ``Pw`` of the
    combined parity operator, momentum and orbitals untouched.

    A plain reversal of the last axis -- correct both on the centered fermionic
    uniform grid (the partner of index ``n`` is ``nmat - 1 - n``) and on the
    symmetric IR fermionic nodes. Invariance under it is what "even frequency"
    means.
    """
    return gap_w[..., ::-1]


def _project_parity_dynamic(gap_w, pairing_type):
    """Project a dynamic gap onto the channel's EVEN-FREQUENCY sector.

    The combined parity ``P = Pk Pw`` is the fermion-antisymmetry constraint
    (``+`` singlet, ``-`` triplet), but each of its eigenspaces still holds an
    even-frequency and an odd-frequency sector. The conventional channels are
    the even-frequency ones, so this is the product of the momentum-parity
    projector and the even-frequency projector,

        ``0.25 * (g + s Pk g + Pw g + s Pk Pw g)``,  ``s = +1`` singlet,
        ``s = -1`` triplet,

    i.e. (even k, even w) for singlet and (odd k, even w) for triplet. It is
    idempotent, and the singlet and triplet projectors are orthogonal. The two
    odd-frequency sectors are annihilated by both: odd-frequency pairing is not
    solved for here (issue #209).
    """
    if pairing_type == "singlet":
        sign = 1.0
    elif pairing_type == "triplet":
        sign = -1.0
    else:
        raise ValueError(
            "Unknown pairing_type: '{}'. Use 'singlet' or 'triplet'.".format(
                pairing_type))
    # Written as the composition (1 + Pw)/2 . (1 + s Pk)/2 rather than the
    # expanded four-term sum: Pk and Pw commute, so the two are the same
    # projector, but on a grid where Pk is the identity (every k its own
    # inverse and a single orbital, e.g. a 2-point axis) the SINGLET
    # composition reduces BIT-identically to the historical (1 + P)/2, so such
    # runs are unchanged to the last digit. For the triplet on such a grid the
    # two do NOT agree and cannot: the odd-k sector is empty there, so this
    # projector is identically zero while the historical (1 - P)/2 was not.
    half = 0.5 * (gap_w + sign * _reverse_k_and_orbital_dynamic(gap_w))
    return 0.5 * (half + _reverse_w_dynamic(half))


# labels and (Pk, Pw) signs of the four sectors, in the order they are written
# to the outputs
_SECTOR_LABELS = ("even_k_even_w", "odd_k_even_w", "even_k_odd_w",
                  "odd_k_odd_w")
_SECTOR_SIGNS = ((1.0, 1.0), (-1.0, 1.0), (1.0, -1.0), (-1.0, -1.0))


def gap_sector_weights(gap_w):
    """Squared-norm fraction of ``gap_w`` in each of the four (momentum,
    frequency) parity sectors of ``Pk`` and ``Pw``.

    Returns a dict keyed by ``even_k_even_w`` (the conventional singlet),
    ``odd_k_even_w`` (the conventional triplet), ``even_k_odd_w`` and
    ``odd_k_odd_w`` (the two odd-frequency sectors this module does not solve
    for). ``even_w`` / ``odd_w`` is the parity under the plain frequency
    reversal ``Pw``; ``even_k`` / ``odd_k`` that under ``Pk``, which carries
    the orbital transpose. The four projectors resolve the identity, so the
    weights sum to 1; a zero gap reports all zeros instead of dividing by zero.
    """
    gap_w = np.asarray(gap_w)
    denom = float(np.linalg.norm(gap_w)) ** 2
    if denom <= 0.0 or not np.isfinite(denom):
        return {label: 0.0 for label in _SECTOR_LABELS}
    gw = _reverse_w_dynamic(gap_w)
    pk = _reverse_k_and_orbital_dynamic(gap_w)
    pkw = _reverse_k_and_orbital_dynamic(gw)
    out = {}
    for label, (sk, sw) in zip(_SECTOR_LABELS, _SECTOR_SIGNS):
        comp = 0.25 * (gap_w + sk * pk + sw * gw + sk * sw * pkw)
        out[label] = float(np.linalg.norm(comp) ** 2 / denom)
    return out


def _is_parity_dynamic(gap_w, pairing_type, tol=0.9):
    """True when ``gap_w`` retains at least ``tol`` of its norm under the
    channel's even-frequency sector projection (mirrors ``sc._is_gap_parity``)."""
    proj = _project_parity_dynamic(gap_w, pairing_type)
    n = np.linalg.norm(gap_w)
    if n == 0:
        return False
    return np.linalg.norm(proj) / n >= tol


#: the combined-parity sign of each channel (the antisymmetry constraint)
_CHANNEL_PARITY_SIGN = {"singlet": 1.0, "triplet": -1.0}


def _project_combined_parity_dynamic(gap_w, pairing_type_or_sign):
    """``(1 + s P)/2`` with the COMBINED parity
    ``P: Delta_ab(k, iw) -> Delta_ba(-k, -iw)``, ``s = +1`` singlet /
    ``-1`` triplet (a channel name or the sign itself is accepted).

    The fermion-antisymmetry constraint itself, which is what the pairing
    kernel actually commutes with (exactly, for every model tested: it is a
    property of the pair bubble, not of the lattice). ``_parity_leakage``
    measures against THIS operator, and the power iteration falls back to it
    when the narrower channel sector is not preserved.
    """
    if isinstance(pairing_type_or_sign, str):
        try:
            sign = _CHANNEL_PARITY_SIGN[pairing_type_or_sign]
        except KeyError:
            raise ValueError(
                "Unknown pairing_type: '{}'. Use 'singlet' or 'triplet'."
                .format(pairing_type_or_sign))
    else:
        sign = float(pairing_type_or_sign)
    return 0.5 * (gap_w + sign * _reverse_kw_and_orbital(gap_w))


def _is_combined_parity_dynamic(gap_w, pairing_type, tol=0.9):
    """True when ``gap_w`` retains at least ``tol`` of its norm under the
    channel's COMBINED-parity projection -- the fallback classification used
    when no eigenpair lies in the channel's even-frequency sector."""
    proj = _project_combined_parity_dynamic(gap_w, pairing_type)
    n = np.linalg.norm(gap_w)
    if n == 0:
        return False
    return np.linalg.norm(proj) / n >= tol


def _parity_leakage(A, gap_shape, pairing_type, n_probe=None, seed=0):
    """Max fraction of ``A x`` that lands in the OPPOSITE combined-parity
    sector.

    Zero (to numerical precision) iff the kernel commutes with the combined
    parity ``P``, i.e. iff the direct-term pairing kernel is the physical
    kernel. Mirrors ``sc._solve_iteration``'s centrosymmetry guard verbatim:
    for each of ``n_probe`` random vectors and each parity sign, project the
    probe into that sector, apply ``A``, and measure the norm fraction of the
    image that lands in the opposite sector (denominator
    ``||A xp|| + 1e-300``). Uses the same probe count as the static path
    (``sc._PARITY_GUARD_PROBES``).

    It is deliberately measured against the COMBINED parity and not against
    the narrower even-frequency channel sector: this number answers "is the
    direct-term kernel the physical kernel?", which is a statement about
    ``P``, and ``P`` is the one thing every kernel conserves. Whether the
    channel sector itself is preserved -- the condition for projecting the
    power iteration onto it -- is a separate question, measured by
    ``_channel_leakage``.
    """
    # pairing_type is unused by design: the guard probes BOTH parity sectors
    # (even and odd), exactly as sc._solve_iteration does, so the leakage
    # measure is independent of the requested channel.
    import hwave.sc as sc
    if n_probe is None:
        n_probe = sc._PARITY_GUARD_PROBES
    rng = np.random.default_rng(seed)
    leakage = 0.0
    for _ in range(n_probe):
        x = (rng.standard_normal(gap_shape)
             + 1j * rng.standard_normal(gap_shape))
        for sign in (1.0, -1.0):
            xp = _project_combined_parity_dynamic(x, sign)   # even (+)/odd (-)
            axp = A.matvec(xp.ravel()).reshape(gap_shape)
            denom = np.linalg.norm(axp) + 1.0e-300
            leakage = max(
                leakage,
                np.linalg.norm(_project_combined_parity_dynamic(axp, -sign))
                / denom)
    return leakage


def _channel_leakage(A, gap_shape, pairing_type, n_probe=None, seed=0):
    """Max fraction of ``A x`` that leaves the CHANNEL sector ``x`` lies in.

    Same probe construction as ``_parity_leakage``, but the probe is projected
    onto the requested channel's even-frequency sector and the measured
    quantity is the norm fraction of the image that is NOT in that sector --
    the other channel and both odd-frequency sectors alike, since all of them
    spoil the projection in the same way.

    Zero (to numerical precision) iff the channel sector is an invariant
    subspace of ``A``, which is exactly the condition under which projecting
    every power iterate onto it still yields an eigenpair of ``A``. It is a
    strictly stronger requirement than ``_parity_leakage``: a kernel can
    commute with the combined parity and still mix the two frequency parities
    inside each of its eigenspaces (a general multi-orbital kernel does).

    Returns ``None`` -- not ``0.0`` -- when the channel sector is EMPTY on this
    grid, so that "nothing to probe" cannot be read as "perfectly preserved".
    That happens for the triplet whenever every k is its own inverse (two
    points or fewer per periodic axis) and there is a single orbital: no
    odd-k function exists at all.
    """
    import hwave.sc as sc
    if n_probe is None:
        n_probe = sc._PARITY_GUARD_PROBES
    rng = np.random.default_rng(seed)
    leakage = None
    for _ in range(n_probe):
        x = (rng.standard_normal(gap_shape)
             + 1j * rng.standard_normal(gap_shape))
        xp = _project_parity_dynamic(x, pairing_type)
        if np.linalg.norm(xp) == 0.0:
            continue                     # empty sector: nothing to probe
        axp = A.matvec(xp.ravel()).reshape(gap_shape)
        denom = np.linalg.norm(axp) + 1.0e-300
        out_of_sector = axp - _project_parity_dynamic(axp, pairing_type)
        this = np.linalg.norm(out_of_sector) / denom
        leakage = this if leakage is None else max(leakage, this)
    return leakage


def _project_seed_dynamic(phi0, pairing_type):
    """Project the seed gap onto the channel sector (no renormalization).

    Mirrors the static guard in ``sc._solve_iteration``: if the seed has no
    component in the requested parity sector, raise instead of silently
    iterating from a (near-)zero vector (which would return a bogus
    ``lambda~0`` non-converged result). Like the static path, the projected
    seed is returned as-is (the power iteration normalizes each iterate, so the
    initial scale does not affect the converged eigenpair).
    """
    proj = _project_parity_dynamic(phi0, pairing_type)
    if np.linalg.norm(proj) < 1.0e-12 * (np.linalg.norm(phi0) + 1.0e-300):
        raise ValueError(
            "Initial gap has no component in the '{}' channel sector (even k "
            "and even frequency for singlet, odd k and even frequency for "
            "triplet); choose an init_gap of the matching k parity, and note "
            "that a purely odd-frequency seed belongs to neither channel."
            .format(pairing_type))
    return proj


def _project_seed_combined_dynamic(phi0, pairing_type):
    """``_project_seed_dynamic`` for the combined-parity fallback sector."""
    proj = _project_combined_parity_dynamic(phi0, pairing_type)
    if np.linalg.norm(proj) < 1.0e-12 * (np.linalg.norm(phi0) + 1.0e-300):
        raise ValueError(
            "Initial gap has no component in the '{}' combined-parity sector "
            "(Delta_ab(k, iw) = {} Delta_ba(-k, -iw)); choose an init_gap of "
            "the matching parity.".format(
                pairing_type,
                "+" if _CHANNEL_PARITY_SIGN[pairing_type] > 0 else "-"))
    return proj


def _reorder_eigenpairs_by_parity_dynamic(vals, vecs, gap_shape, pairing_type):
    """Promote eigenpairs whose gap lies in the channel's even-frequency sector.

    Mirrors ``sc._reorder_eigenpairs_by_parity`` for the frequency-resolved
    gap so the reported leading dynamic eigenpair is the physical solution for
    the requested channel, not ARPACK's raw largest-|lambda| (which can be an
    opposite-parity mode).

    Parameters
    ----------
    vals : ndarray
        Eigenvalues, shape ``(k,)``, already ordered by descending real part.
    vecs : ndarray
        Eigenvectors as columns, shape ``(vec_size, k)``.
    gap_shape : tuple
        ``(norb, norb, Nx, Ny, Nz, nmat)`` to reshape each column.
    pairing_type : str
        "singlet" or "triplet".

    Matching runs in two stages. First the channel's even-frequency sector.
    Only if NO eigenpair lies in it -- which is what happens when the kernel
    does not conserve the frequency parity, so its eigenvectors are
    even/odd-frequency mixtures of definite combined parity -- the channel's
    COMBINED-parity sector is tried instead, with a warning. The returned
    ``match`` array belongs to whichever stage was used.

    Returns
    -------
    vals, vecs, match, selection
        Reordered with the matching eigenpairs first (order preserved within
        each group), the reordered boolean match array, and the stage that
        produced it: ``"channel"``, ``"combined_parity"``, or ``"none"`` when
        neither matched anything (``match`` is then all False). The caller must
        carry ``selection`` into the outputs: the same column of ones means a
        different thing in each case.
    """
    def _flag(test):
        return np.array([test(vecs[:, i].reshape(gap_shape), pairing_type)
                         for i in range(vecs.shape[1])])

    selection = "channel"
    match = _flag(_is_parity_dynamic)
    if not np.any(match):
        fallback = _flag(_is_combined_parity_dynamic)
        if np.any(fallback):
            selection = "combined_parity"
            logger.warning(
                "No dynamic eigenpair lies in the '%s' even-frequency sector; "
                "the frequency parity is not a symmetry of this kernel, so the "
                "eigenpairs of the channel's combined-parity sector are "
                "promoted instead. See gap_sector_weights for the composition "
                "of the reported gap.", pairing_type)
            match = fallback
        else:
            selection = "none"
            logger.warning(
                "No dynamic eigenpair matches the requested '%s' parity; the "
                "reported leading gap belongs to the other channel. Increase "
                "num_eigenvalues or check pairing_type.", pairing_type)
    idx = np.concatenate([np.where(match)[0], np.where(~match)[0]])
    return vals[idx], vecs[:, idx], match[idx], selection


def _channel_projection_is_valid(channel_leakage, parity_leakage_tol):
    """Whether the channel's even-frequency sector may be projected onto.

    The same gate the projected power iteration uses: the sector must exist on
    this grid (``channel_leakage`` is not ``None``) and the kernel must commute
    with the channel projector within tolerance, so that restricting the
    eigenproblem to the sector still yields eigenpairs of the kernel.
    """
    return channel_leakage is not None and channel_leakage <= parity_leakage_tol


def _solve_channel_projected(make_operator, matvec, vec_size, gap_shape,
                             pairing_type, num_eigenvalues):
    """Largest-real eigenpair of the kernel restricted to the channel sector.

    Issue #202: a plain ``which='LM'`` ARPACK set asks for the eigenvalues of
    largest MAGNITUDE, so on a kernel that CONSERVES the channel's
    even-frequency sector the channel's (small, positive) Tc eigenvalue can be
    masked by larger-real opposite-sector modes and fall beyond the requested
    ``num_eigenvalues`` -- making the reported value depend on that count. When
    the sector is a valid invariant subspace (see
    ``_channel_projection_is_valid``) we instead solve the projected operator
    ``P K P`` with the channel projector ``P`` (``_project_parity_dynamic``, the
    even-frequency projector of PR #211), whose spectrum on the sector is
    exactly the channel's, and select the largest real part -- which is
    ``num_eigenvalues``-independent.

    The orthogonal complement (everything ``P`` annihilates) is mapped to a
    large NEGATIVE eigenvalue rather than left at zero, so the null space cannot
    win the largest-real selection and mask a repulsive-dominant sector whose
    own leading eigenvalue is negative; the reported value is then the least
    repulsive channel eigenvalue, not a spurious zero.

    Returns the same ``(eigenvalue, eigenvector, info)`` triple as
    ``sc._solve_leading`` (with the complement eigenvalues subtracted back out
    of range), so the caller reorders and gauge-fixes it exactly as the
    unprojected result.
    """
    import hwave.sc as sc
    from scipy.sparse.linalg import LinearOperator, eigs

    # scale of the full kernel, to size the complement penalty relative to the
    # sector spectrum (a cheap preliminary, like the auto-shift estimate)
    A0, _ = make_operator()
    k_pre = min(6, vec_size - 2)
    rho = 1.0
    if k_pre >= 1:
        try:
            pre_vals, _ = eigs(A0, k=k_pre, which='LM')
            rho = float(np.max(np.abs(pre_vals)))
        except Exception:
            rho = 1.0
    if not np.isfinite(rho) or rho <= 0.0:
        rho = 1.0
    penalty = 4.0 * rho + 1.0

    def _P(flat):
        return _project_parity_dynamic(
            np.asarray(flat).reshape(gap_shape), pairing_type).ravel()

    def projected_matvec(v):
        v = np.asarray(v).ravel()
        pv = _P(v)
        kpv = _P(np.asarray(matvec(pv)).ravel())
        return kpv - penalty * (v - pv)

    def make_projected():
        return (LinearOperator((vec_size, vec_size), matvec=projected_matvec,
                               dtype=complex),
                vec_size)

    # An explicit spectral shift large enough that A + sigma I has an
    # all-positive-real spectrum (sector eigenvalues >= -rho, complement
    # -penalty), so ARPACK's which='LR' iteration is well conditioned; it is
    # subtracted back inside _solve_leading. The tiny-operator dense path
    # ignores the shift and orders by real part directly, which is equivalent.
    sigma = 2.0 * penalty + rho + 1.0
    return sc._solve_leading(
        make_projected, vec_size, "arnoldi",
        num_eigenvalues=num_eigenvalues, spectral_shift=sigma, seed_vec=None)


def estimate_memory_bytes(norb, Nk, nmat):
    v = 16 * norb**4 * Nk * nmat          # vertex
    g2 = v                                # G2, same shape
    vec = 16 * norb**2 * Nk * nmat        # one gap/Arnoldi vector
    return int((v + g2 + 4 * vec) * 1.5)  # + ~4 work vectors, 1.5x FFT headroom


def _available_ram_bytes():
    try:
        import psutil
        return int(psutil.virtual_memory().available)
    except Exception:
        try:
            with open("/proc/meminfo") as fh:
                for line in fh:
                    if line.startswith("MemAvailable:"):
                        return int(line.split()[1]) * 1024
        except Exception:
            pass
    return 8 * 2**30  # fallback: 8 GiB


def check_memory(norb, Nk, nmat, mem_limit_gb=None):
    need = estimate_memory_bytes(norb, Nk, nmat)
    if mem_limit_gb == 0:
        logger.info("dynamic Eliashberg est. peak %.2f GiB (guard disabled)",
                    need / 2**30)
        return
    limit = ((0.8 * _available_ram_bytes()) if mem_limit_gb is None
             else mem_limit_gb * 2**30)
    logger.info("dynamic Eliashberg est. peak %.2f GiB (limit %.2f GiB)",
                need / 2**30, limit / 2**30)
    if need > limit:
        raise MemoryError(
            "dynamic Eliashberg estimated peak {:.1f} GiB exceeds the limit "
            "{:.1f} GiB; reduce Nmat/Nk or set [eliashberg] mem_limit_gb "
            "(0 disables).".format(need / 2**30, limit / 2**30))


# NOTE: do NOT write a private chi loader here — path/name/convention/
# spin-orbital-expansion logic lives ONLY in sc._load_flex_susceptibilities_full
# (Task 1c). load_flex_chi_dynamic delegates to it to avoid contract drift.
# The module also needs the FFT alias for the Task-6 kernel:
#   from hwave.solver.perf import FFT

def _npz_freq_size(path, keys, axis):
    """Frequency-axis length of an array in an ``.npz`` read from its header
    only (no data load).

    Parameters
    ----------
    path : str
        Path to the ``.npz`` file.
    keys : tuple of str
        Candidate array names to try, in order (e.g. ``("chiq_s", "chiq")``).
    axis : int
        Axis holding the Matsubara frequency in the H-wave layout (chi: 0;
        green: 1).

    Returns
    -------
    int or None
        The size of ``axis`` for the first present key, or ``None`` if the file
        is absent/unreadable or none of the keys is in the archive (the caller
        then falls back to the config grid and lets the loader raise the proper
        missing-file error).
    """
    import zipfile

    from hwave.solver import npy_header as _npy_header

    # Best-effort header probe. This is a pure optimization/safety pre-check --
    # the loader below is the authoritative path -- so ANY failure returns None
    # and lets the loader raise the existing, clearer error. The broad catch is
    # deliberate: besides unreadable/malformed/truncated files and a missing
    # axis, it also covers a future NPY header version the shared reader in
    # hwave.solver.npy_header does not know, which must degrade gracefully
    # rather than crash.
    try:
        with zipfile.ZipFile(path) as z:
            names = set(z.namelist())
            for key in keys:
                member = key + ".npy"
                if member in names:
                    with z.open(member) as f:
                        shape = _npy_header.read_npy_header_shape(f)
                    return int(shape[axis])
    except Exception:
        return None
    return None


def _ir_validate_native_nodes(arr, meta, ax, label, beta):
    """Validate native-node metadata without fitting an unused channel."""
    statistics = str(meta["statistics"])
    if statistics != ax.statistics:
        raise ValueError(
            "IR-native {}: ir_statistics={!r} does not match the run's "
            "{} axis ({!r}).".format(
                label,
                statistics,
                "fermionic" if ax.statistics == "F" else "bosonic",
                ax.statistics,
            )
        )
    file_beta = float(meta["beta"])
    if not np.isclose(file_beta, beta, rtol=1e-9, atol=1e-9 * beta):
        raise ValueError(
            "IR-native {}: ir_beta {:.12g} differs from this run's beta "
            "{:.12g} (rel {:.3e}); the susceptibilities/Green function are "
            "physics input and must match. Re-run FLEX at this temperature."
            .format(label, file_beta, beta,
                    abs(file_beta - beta) / abs(beta)))
    freq_n = np.asarray(meta["freq_n"], dtype=np.int64)
    if arr.shape[-1] != freq_n.size:
        raise ValueError(
            "IR-native {}: stored frequency-axis length {} differs from "
            "len(ir_freq_n)={}.".format(label, arr.shape[-1], freq_n.size)
        )
    return freq_n


def _ir_refit_nodes(arr, meta, ax, label, beta):
    """Bring IR-native node values (..., n_file_nodes) onto the RUN's axis
    nodes (design ir-matsubara-stage3.md Sec. 4.1). Returns NODE VALUES.

    Exact node-set equality is a pure pass-through (the common case: same
    beta, auto wmax, same sparse-ir version) -- a fit/eval round trip would
    not be zero-cost and could perturb the values.
    """
    freq_n = _ir_validate_native_nodes(arr, meta, ax, label, beta)
    if np.array_equal(freq_n, ax.freq_n):
        logger.info("IR-native %s: file node set equals the run basis "
                    "(%d nodes); stored values used directly.", label,
                    freq_n.size)
        return arr
    if ax.wmax < meta["wmax"] or ax.eps > meta["tol"]:
        logger.warning(
            "IR-native %s: the run basis is weaker than the writer's "
            "(wmax %.4g vs %.4g, eps %.1e vs %.1e); refit quality is "
            "empirical -- consider matching [eliashberg] ir_wmax/ir_tol to "
            "the FLEX run.", label, ax.wmax, meta["wmax"], ax.eps,
            meta["tol"])
    coeffs = ax.fit_from_freq_points(arr, freq_n)
    resid = float(np.abs(
        ax.eval_to_freq_points(coeffs, freq_n) - arr).max())
    scale = float(np.abs(arr).max()) or 1.0
    logger.info("IR-native %s: refit %d file nodes -> %d run nodes (L=%d), "
                "max residual at file nodes %.3e (rel %.3e)", label,
                freq_n.size, ax.n_freq, ax.L, resid, resid / scale)
    if resid > 0.05 * scale:
        logger.warning(
            "IR-native %s: refit residual is large (>5%% of the data "
            "scale); the run basis cannot represent the file's content -- "
            "raise [eliashberg] ir_wmax or tighten [eliashberg] ir_tol.",
            label)
    return ax.eval_to_freq(coeffs)


def _npz_is_ir_native(path):
    """Header-level D-1 discriminator (design ir-matsubara-stage3.md): a
    lazily-opened npz reads only the two small metadata keys, so this stays
    a preflight (no large-array load). Missing/unreadable files return
    False -- the authoritative loader raises the clearer error."""
    from hwave.solver.ir_axis import is_ir_native
    try:
        with np.load(path) as data:
            return is_ir_native(data)
    except Exception:
        return False


def load_flex_chi_dynamic(input_dict, norb, Nx, Ny, Nz, allow_ir=False,
                          interactions=None):
    import hwave.sc as sc
    cfg_nmat = int(input_dict["mode"]["param"].get("Nmat", 1024))
    if cfg_nmat % 2 != 0:
        raise ValueError("dynamic Eliashberg requires even Nmat; got {}".format(cfg_nmat))

    # Fail BEFORE allocating: read the STORED frequency-axis sizes from the NPZ
    # headers (no data load), validate them against the config Nmat, and size
    # the memory guard on the actual files. Using only the config Nmat would let
    # a larger on-disk grid slip past the guard and OOM inside the loader before
    # the post-load grid check could fire.
    chi_s_path, chi_c_path, green_path = sc._resolve_flex_paths(input_dict)
    # Stage 3, D-1 FIRST (design Sec. 4.1): for IR-native files the stored
    # frequency-axis length is the sparse node count, not Nmat, so the
    # uniform-grid equality check below must not see them.
    native = {"chis": _npz_is_ir_native(chi_s_path),
              "chic": _npz_is_ir_native(chi_c_path)}
    if os.path.exists(green_path):
        native["green"] = _npz_is_ir_native(green_path)
    any_native = any(native.values())
    if any_native and not allow_ir:
        first = next(k for k, v in native.items() if v)
        paths = {"chis": chi_s_path, "chic": chi_c_path, "green": green_path}
        raise ValueError(
            "file '{}' holds sparse-IR node data "
            "(frequency_grid=sparse_ir_nodes); set [eliashberg] "
            "matsubara_basis = \"ir\" to consume it, or re-run FLEX with "
            "[mode.param] write_densified = true.".format(paths[first]))
    stored = {"chis": _npz_freq_size(chi_s_path, ("chiq_s", "chiq"), axis=0),
              "chic": _npz_freq_size(chi_c_path, ("chiq_c", "chiq"), axis=0)}
    if os.path.exists(green_path):
        stored["green"] = _npz_freq_size(green_path, ("green",), axis=1)
    if not any_native:
        for name, n in stored.items():
            if n is not None and n != cfg_nmat:
                raise ValueError(
                    "dynamic Eliashberg grid mismatch: {} nmat={} differs "
                    "from config Nmat={}".format(name, n, cfg_nmat))
        # sizes agree with the config here, but guard on the stored value so
        # the estimate always reflects what is on disk.
        file_nmat = max([n for n in stored.values() if n is not None]
                        + [cfg_nmat])
    else:
        # IR-native: every large tensor downstream lives on the node axis
        # (design Sec. 4.1) -- size the guard from the node counts. Read
        # them from the tiny ir_freq_n member itself rather than the
        # best-effort header probe above: on numpy >= 2.4 the probe's
        # private-API path degrades to None (by design), which would
        # silently undersize the guard here. Unreadable files fall back to
        # the conservative config Nmat.
        counts = []
        for path, is_native in ((chi_s_path, native.get("chis")),
                                (chi_c_path, native.get("chic")),
                                (green_path, native.get("green"))):
            if is_native:
                try:
                    with np.load(path) as d:
                        counts.append(int(np.asarray(d["ir_freq_n"]).size))
                except Exception:
                    pass
        file_nmat = max(counts) if counts else cfg_nmat
    check_memory(norb, Nx * Ny * Nz, file_nmat,
                 input_dict["eliashberg"].get("mem_limit_gb"))

    # reuse the exact static path/name/convention/spin-orbital-expansion logic
    if allow_ir:
        chis_w, chic_w, green_w, chi_convention, ir_meta = \
            sc._load_flex_susceptibilities_full(input_dict, norb, Nx, Ny, Nz,
                                                allow_ir=True,
                                                interactions=interactions)
    else:
        chis_w, chic_w, green_w, chi_convention = \
            sc._load_flex_susceptibilities_full(input_dict, norb, Nx, Ny, Nz,
                                                interactions=interactions)
        ir_meta = None
    # Belt-and-suspenders: the header check above already rejected a mismatch,
    # but keep a post-load assertion in case the loader reshapes unexpectedly.
    if ir_meta is None:
        green_nmat = green_w.shape[-1] if green_w is not None else cfg_nmat
        if not (chis_w.shape[-1] == chic_w.shape[-1] == green_nmat
                == cfg_nmat):
            raise ValueError(
                "dynamic Eliashberg grid mismatch: nmat differs — chis={}, "
                "chic={}, green={}, config Nmat={}".format(
                    chis_w.shape[-1], chic_w.shape[-1], green_nmat, cfg_nmat))
    else:
        for label, arr, key in (("chiq_s", chis_w, "chis"),
                                ("chiq_c", chic_w, "chic"),
                                ("green", green_w, "green")):
            if arr is not None and arr.shape[-1] != ir_meta[key][
                    "freq_n"].size:
                raise ValueError(
                    "IR-native {}: stored axis length {} does not match its "
                    "own ir_freq_n count {} -- the file is inconsistent."
                    .format(label, arr.shape[-1],
                            ir_meta[key]["freq_n"].size))
    if not allow_ir:
        return chis_w, chic_w, green_w, chi_convention
    return chis_w, chic_w, green_w, chi_convention, ir_meta


def compute_vertices_flex_dynamic(chis_w, chic_w, inter_k, norb,
                                  Nx, Ny, Nz, pairing_type, convention,
                                  sc_matrices=None):
    """Full-frequency pairing vertex: apply sc._compute_vertices_flex per
    bosonic Matsubara frequency and stack along the trailing axis.

    Parameters
    ----------
    chis_w, chic_w : ndarray
        Spin/charge susceptibilities, shape (Nx, Ny, Nz, nd, nd, nmat)
        with nd = norb**2.
    inter_k : dict
        Interactions in k-space from sc._build_interaction_k.
    norb, Nx, Ny, Nz : int
        Orbital count and grid dimensions.
    pairing_type : str
        "singlet" or "triplet".
    convention : str
        "kuroki" or "myo" — MUST match the orbital convention chis_w/chic_w
        were produced in (chi_convention from the FLEX loader); forwarded
        unchanged to sc._compute_vertices_flex so the matching S/C matrices
        are used at every frequency.

    Returns
    -------
    Vs_q_w : ndarray
        Pairing vertex, shape (norb, norb, norb, norb, Nx, Ny, Nz, nmat).
    """
    import hwave.sc as sc
    nmat = chis_w.shape[-1]

    # Reject BEFORE the S/C build: the pair is O(Nq * norb^4) and an
    # unsupported calculation must not allocate heavily on its way to the
    # validation error (round-10 review).
    sc._reject_reduced_flex_unsupported(inter_k, convention)

    # The S/C matrices are frequency-independent: build ONE pair for the
    # whole run (or take the caller's) and hand it to the diagnostic and to
    # every per-frequency contraction. Previously each of the ~nmat calls
    # rebuilt an identical full-grid pair (round-9 review).
    if sc_matrices is None:
        sc_matrices = sc._build_vertex_sc_matrices(convention, inter_k,
                                                   norb, Nx, Ny, Nz)

    # Once for the whole run: _compute_vertices_flex below is called per
    # frequency, so warning inside it would repeat nmat (~1000) times.
    sc._warn_reduced_flex_missing_components(
        inter_k, norb, Nx, Ny, Nz, convention,
        sc_matrices=(sc_matrices if str(convention).lower() == "kuroki"
                     else None))

    def _one(l):
        return sc._compute_vertices_flex(
            chis_w[..., l], chic_w[..., l], inter_k, norb, Nx, Ny, Nz,
            pairing_type=pairing_type, convention=convention,
            sc_matrices=sc_matrices)

    v0 = _one(0)
    out = np.empty(v0.shape + (nmat,), dtype=v0.dtype)
    out[..., 0] = v0
    for l in range(1, nmat):
        out[..., l] = _one(l)
    return out


def calc_g2_dynamic(green_kw, beta):
    """Frequency-resolved pair bubble: identical to sc._calc_g2 except the
    fermionic Matsubara sum is NOT taken (the frequency axis is kept).

    sc._calc_g2 computes G2[i,j,l,m,x,y,z] = (1/beta) * sum_n
    green_kw[i,j,x,y,z,n] * green_kw_inv[l,m,x,y,z,n], where green_kw_inv is
    G(-k,-wn) built via the shared FFT-grid reversal (kgrid.reverse_fft_axes,
    i -> (N - i) % N). This function drops the sum over n and
    returns the per-frequency summand, so calc_g2_dynamic(...).sum(axis=-1)
    reproduces sc._calc_g2(..., tail=False) to machine precision (see
    tests/test_eliashberg_dynamic.py::test_g2_dynamic_sums_to_static).
    The static path's Matsubara tail correction (issue #86) is a property of
    taking the frequency sum, so it has no per-frequency counterpart here;
    the dynamic kernel keeps the bare summand.

    Parameters
    ----------
    green_kw : ndarray
        Green's function, shape (norb, norb, Nx, Ny, Nz, nmat).
    beta : float
        Inverse temperature.

    Returns
    -------
    G2_w : ndarray
        Shape (norb, norb, norb, norb, Nx, Ny, Nz, nmat).
    """
    norb = green_kw.shape[0]
    Nx, Ny, Nz, nmat = green_kw.shape[2], green_kw.shape[3], green_kw.shape[4], green_kw.shape[5]
    nvol = Nx * Ny * Nz

    # G(-k, -wn) via the shared FFT-grid reversal -- SAME construction
    # as sc._calc_g2.
    green_kw_inv = reverse_fft_axes(green_kw[..., ::-1], (2, 3, 4))
    # Same reshape/index layout as sc._calc_g2's A/B (ij, site, n) and
    # (lm, site, n), but keep the per-frequency product instead of summing
    # (matmul-)contracting over n.
    A = green_kw.reshape(norb * norb, nvol, nmat)      # (ij, site, n)
    B = green_kw_inv.reshape(norb * norb, nvol, nmat)  # (lm, site, n)
    # G2[ij, lm, site, n] = A[ij, site, n] * B[lm, site, n]  (no sum over n)
    G2 = A[:, np.newaxis, :, :] * B[np.newaxis, :, :, :]
    G2 = G2.reshape(norb, norb, norb, norb, Nx, Ny, Nz, nmat)
    return G2 / beta


def vertex_qw_to_rt(Vs_q_w, workers=1):
    r"""Transform the pairing vertex from (q, i nu_l) to (r, tau).

    This is the vertex leg of ``eliashberg_kernel_dynamic``: spatial q->r via
    ifftn (which carries the single spatial fold's 1/N), frequency boson->tau.
    It does not depend on the trial gap, so callers that apply the kernel
    repeatedly (power iteration / Arnoldi) should compute it once and pass the
    result via the kernel's ``Vs_rt`` argument.

    Parameters
    ----------
    Vs_q_w : ndarray
        Pairing vertex, shape (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
        on the bosonic Matsubara axis. A cupy array is transformed on the
        device.
    workers : int, optional
        FFT worker threads for the spatial transform on the numpy backend
        (see :func:`_spatial_ifftn`); ignored on the cupy backend.

    Returns
    -------
    V_rt : ndarray
        Same shape, in (r, tau), on the same backend as the input.
    """
    return _spatial_ifftn(ms.boson_to_tau(Vs_q_w, axis=-1),
                          axes=(4, 5, 6), workers=workers)


def eliashberg_kernel_dynamic(Vs_q_w, G2_w, phi_w, norb, beta, Vs_rt=None,
                              workers=1):
    r"""Apply the frequency-resolved (tau-product) Eliashberg kernel.

    Implements one action of the linearized Eliashberg operator on a trial
    gap ``phi``, keeping the full fermionic Matsubara axis. The structure
    mirrors the static kernel (``sc._make_kernel_operator``) and the FLEX
    self-energy (``flex._calc_self_energy``): an orbital contraction with the
    pair bubble G2, then an imaginary-time PRODUCT with the pairing vertex
    (NOT a circular frequency convolution) done via matsubara transforms.

    Orbital contraction (matches the static kernel):
        F_{l2,l3}(k,m)     = sum_{l5,l6} G2_{l2,l5,l3,l6}(k,m) phi_{l5,l6}(k,m)
        phi_out_{l1,l4}    = - sum_{l2,l3} V_{l1,l2,l3,l4}(r,tau) F_{l2,l3}(r,tau)

    The temperature factor 1/beta lives inside G2 (``calc_g2_dynamic``); the
    kernel carries only the spatial -(1/N) via the ifftn (numpy) convention.
    ``beta`` is accepted for signature symmetry with the static path but is
    not applied here.

    Parameters
    ----------
    Vs_q_w : ndarray or None
        Pairing vertex, shape (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
        on the bosonic Matsubara axis. May be None when ``Vs_rt`` is given.
    G2_w : ndarray
        Pair bubble, shape (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
        on the fermionic Matsubara axis (already divided by beta).
    phi_w : ndarray
        Trial gap, shape (norb, norb, Nx, Ny, Nz, nmat).
    norb : int
        Number of orbitals (for signature symmetry; inferred from shapes).
    beta : float
        Inverse temperature (unused here; T is inside G2).
    Vs_rt : ndarray, optional
        Precomputed ``vertex_qw_to_rt(Vs_q_w)``. The vertex transform does not
        depend on ``phi_w``, so iterative solvers pass it once instead of
        paying the (norb^4 x Nvol x nmat)-sized transform on every matvec.
    workers : int, optional
        FFT worker threads for the spatial transforms on the numpy backend
        (see :func:`_spatial_ifftn`); ignored on the cupy backend.

    Returns
    -------
    phi_out_w : ndarray
        Same shape as ``phi_w``, on the same backend as ``G2_w``.

    Notes
    -----
    The array backend follows ``G2_w``: when ``G2_w`` (and ``Vs_rt``) are cupy
    arrays the whole kernel runs on the GPU, and a numpy ``phi_w`` (as handed
    over by the host-side iterative solvers) is transferred to the device
    here.
    """
    xp = backend.array_module_of(G2_w)
    phi_w = xp.asarray(phi_w)
    # F_{l2,l3}(k, m) = sum_{l5,l6} G2_{l2,l5,l3,l6}(k,m) phi_{l5,l6}(k,m)
    F = xp.einsum('iljmxyzn,lmxyzn->ijxyzn', G2_w, phi_w)
    # spatial k->r on F (per orbital pair, per fermionic freq); freq fermion->tau
    F_rt = _spatial_ifftn(ms.fermion_to_tau(F, axis=-1),
                          axes=(2, 3, 4), workers=workers)
    # V(q, iv_l) -> (r, tau): spatial q->r, freq boson->tau. Normalize the
    # vertex to G2's backend: a caller may pass a host (numpy) Vs_rt/Vs_q_w with
    # a device (cupy) G2_w, and the einsum below must not mix host/device arrays.
    V_rt = vertex_qw_to_rt(Vs_q_w, workers=workers) if Vs_rt is None else Vs_rt
    V_rt = xp.asarray(V_rt)
    # phi_out_{l1,l4}(r,tau) = - sum_{l2,l3} V_{l1,l2,l3,l4}(r,tau) F_{l2,l3}(r,tau)
    prod = -xp.einsum('abcdxyzt,bcxyzt->adxyzt', V_rt, F_rt)
    # back: spatial r->k (fftn), freq tau->fermion. The single spatial fold's
    # -(1/N) is already carried by the ifftn above (numpy divides by N on the
    # inverse transform; the r->k fftn does not re-multiply the pair bubble).
    phi_out = ms.tau_to_fermion(_spatial_fftn(prod, axes=(2, 3, 4),
                                              workers=workers), axis=-1)
    return phi_out


def frequency_inner(a, b):
    """Inner product over ALL components (orbital, k, and frequency axes)."""
    return np.vdot(a, b)


def _fix_gauge(phi_w):
    r"""Fix the eigenvector gauge deterministically.

    An eigenvector of the (complex, non-Hermitian) Eliashberg kernel is only
    defined up to an overall complex scale. To make the reported gap
    reproducible across runs / linear-algebra backends we pin two freedoms:

    1. **Magnitude** — L2-normalize over *all* components (orbital, k, and the
       full Matsubara axis) so ``||phi|| = 1``.
    2. **Global phase** — multiply by the single phase that makes the
       largest-``|magnitude|`` component real and positive. Ties in the
       magnitude are broken by the first component in lexicographic
       ``(orb1, orb2, kx, ky, kz, iomega)`` order, which is exactly C-order of
       the ``(norb, norb, Nx, Ny, Nz, nmat)`` array (so ``np.argmax`` on the
       raveled magnitudes -- which returns the first maximizing index -- gives
       the tie-break for free).

    The magnitudes are invariant under a global phase, so applying an arbitrary
    phase before ``_fix_gauge`` reproduces the identical array (see
    ``tests/test_eliashberg_dynamic.py::test_gauge_deterministic``).

    Parameters
    ----------
    phi_w : ndarray
        Gap / eigenvector, shape (norb, norb, Nx, Ny, Nz, nmat).

    Returns
    -------
    ndarray
        Same shape as ``phi_w``, L2-normalized and phase-fixed.
    """
    phi = np.asarray(phi_w).astype(complex, copy=True)
    nrm = np.linalg.norm(phi)
    if nrm > 0:
        phi /= nrm
    flat = phi.ravel()  # C-order == lexicographic (orb1,orb2,kx,ky,kz,iomega)
    pivot = flat[int(np.argmax(np.abs(flat)))]
    if pivot != 0:
        phi /= (pivot / abs(pivot))  # rotate so the pivot becomes real-positive
    return phi


def write_dynamic_outputs(output_dir, gap_w, eigenvalue, T, pairing_type,
                          kx_array, ky_array, kz_array, beta,
                          gap_file="gap.dat", npz_file="gap_dynamic.npz",
                          extra_meta=None, sector_weights=None,
                          selection=None, eigenvalue_selection=None):
    r"""Write the dynamic-Eliashberg gap outputs.

    Produces two files under ``output_dir``:

    * ``gap_dynamic.npz`` -- the full frequency-resolved gap plus metadata.
      Keys: ``gap`` (norb, norb, Nx, Ny, Nz, nmat), ``iomega`` (the centered
      fermionic Matsubara frequencies :math:`\omega_n=(2n+1-N_{mat})\pi T`),
      ``T``, ``pairing_type``, ``frequency`` (== "dynamic"), ``eigenvalue``,
      ``axis_order``, and ``normalization`` (documenting the ``_fix_gauge``
      convention).
    * ``gap.dat`` -- a plain-text slice at a single Matsubara frequency (the
      smallest positive :math:`\omega_n`, index ``nmat//2``), mirroring the
      static ``sc._save_results`` column layout (kx ky kz then Re/Im per orbital
      pair). Its FIRST line is a ``#``-prefixed header that carries
      ``frequency=dynamic`` together with the slice index and its
      :math:`\omega_n` value.

    Parameters
    ----------
    output_dir : str
        Destination directory (created if absent).
    gap_w : ndarray
        Gauge-fixed gap, shape (norb, norb, Nx, Ny, Nz, nmat).
    eigenvalue : float
        Leading eigenvalue lambda.
    T : float
        Temperature.
    pairing_type : str
        "singlet" or "triplet".
    kx_array, ky_array, kz_array : ndarray
        k-point arrays.
    beta : float
        Inverse temperature (accepted for signature symmetry; iomega uses T).
    gap_file, npz_file : str
        Output filenames.
    sector_weights : dict, optional
        ``gap_sector_weights(gap_w)``. When given, the npz also carries
        ``gap_sector_weights`` (the four fractions, in the order of
        ``gap_sector_labels``) and ``gap_sector_labels``. A caller that already
        put those two keys into ``extra_meta`` (the in-process bond writer)
        must not pass this as well; ``extra_meta`` wins.
    selection : str, optional
        ``run_leading_eigenproblem``'s ``sector_selection``, written under that
        name. As with ``sector_weights``, ``extra_meta`` wins when it already
        carries the key.
    eigenvalue_selection : str, optional
        ``run_leading_eigenproblem``'s ``eigenvalue_selection`` (issue #202),
        written under that name when not ``None``. As with ``selection``,
        ``extra_meta`` wins when it already carries the key.
    """
    os.makedirs(output_dir, exist_ok=True)
    norb = gap_w.shape[0]
    Nx, Ny, Nz, nmat = gap_w.shape[2:6]

    # centered fermionic Matsubara frequencies: iw_n = (2n + 1 - Nmat) pi T
    n_idx = np.arange(nmat)
    iomega = (2.0 * n_idx + 1.0 - nmat) * np.pi * T

    axis_order = "(orb1, orb2, kx, ky, kz, iomega)"
    normalization = ("L2-normalized over all (orb,k,iomega) components; global "
                     "phase fixes the largest-|magnitude| component real-"
                     "positive (lexicographic (orb1,orb2,kx,ky,kz,iomega) "
                     "tie-break)")

    meta = dict(extra_meta or {})
    if sector_weights is not None:
        meta.setdefault("gap_sector_weights",
                        np.array([float(sector_weights[label])
                                  for label in _SECTOR_LABELS]))
        meta.setdefault("gap_sector_labels", np.array(_SECTOR_LABELS))
    if selection is not None:
        meta.setdefault("sector_selection", str(selection))
    if eigenvalue_selection is not None:
        meta.setdefault("eigenvalue_selection", str(eigenvalue_selection))

    np.savez(
        os.path.join(output_dir, npz_file),
        gap=gap_w,
        iomega=iomega,
        T=T,
        pairing_type=pairing_type,
        frequency="dynamic",
        eigenvalue=eigenvalue,
        axis_order=axis_order,
        normalization=normalization,
        # Fourier-sign provenance (issue #133): the gap is k-resolved
        momentum_convention="e_plus_ikR",  # = rpa.MOMENTUM_CONVENTION
        **meta,
    )

    # gap.dat: the fermionic slice nearest omega = 0^+ (smallest positive w_n).
    n0 = nmat // 2
    logger.info("Saving dynamic gap slice (index %d, iw_n=%.6e) to %s",
                n0, iomega[n0], os.path.join(output_dir, gap_file))
    with open(os.path.join(output_dir, gap_file), "w") as fw:
        header = ["# frequency=dynamic",
                  "index={}".format(n0),
                  "iomega_n={:.8e}".format(iomega[n0]),
                  "pairing_type={}".format(pairing_type),
                  "eigenvalue={:.8e}".format(eigenvalue),
                  "T={:.8e}".format(T)]
        if extra_meta and (
            extra_meta.get("zero_chi_s") or extra_meta.get("zero_chi_c")
        ):
            header.extend(
                [
                    "zero_chi_s={}".format(
                        str(bool(extra_meta.get("zero_chi_s", False))).lower()
                    ),
                    "zero_chi_c={}".format(
                        str(bool(extra_meta.get("zero_chi_c", False))).lower()
                    ),
                ]
            )
        fw.write("  ".join(header) + "\n")
        cols = ["# kx", "ky", "kz"]
        for i in range(norb):
            for j in range(norb):
                cols.append("Re(sigma_{}{})".format(i, j))
                cols.append("Im(sigma_{}{})".format(i, j))
        fw.write(" ".join(cols) + "\n")
        for ix in range(Nx):
            kx = kx_array[ix]
            if kx > np.pi:
                kx -= 2.0 * np.pi
            for iy in range(Ny):
                ky = ky_array[iy]
                if ky > np.pi:
                    ky -= 2.0 * np.pi
                for iz in range(Nz):
                    kz = kz_array[iz]
                    if kz > np.pi:
                        kz -= 2.0 * np.pi
                    parts = ["{:.8f}".format(kx), "{:.8f}".format(ky),
                             "{:.8f}".format(kz)]
                    for i in range(norb):
                        for j in range(norb):
                            val = gap_w[i, j, ix, iy, iz, n0]
                            parts.append("{:.8e}".format(val.real))
                            parts.append("{:.8e}".format(val.imag))
                    fw.write(" ".join(parts) + "\n")


def _ir_auto_wmax(hr, inter_k, norb, beta, mu=None, filling=None):
    """Heuristic default for ir_wmax (design Sec. 4): 3x the sum of the
    single-particle spectral half-range and the largest interaction scale.
    Aborts (ValueError) when the estimate cannot be formed, rather than
    silently defaulting.

    The spectral half-range is ``max|eps_k - mu|`` from the actual dispersion
    eps(k), diagonalized on a coarse k-mesh built from the real-space transfer
    integrals; this is the real-frequency extent the IR basis must cover. The
    chemical potential ``mu`` is used directly when given, else solved from
    ``filling`` via ``sc._determine_mu`` (mu is what sets where the spectral
    weight sits relative to zero -- ignoring it and using max|eps_k| would
    re-introduce any on-site energy offset). ``hr`` is the flat wan90 layout
    ``{((Rx,Ry,Rz),(orb1,orb2)): scalar}`` (see ``sc._build_hamiltonian_k``);
    it must NOT be summed as if each value were a per-R matrix -- doing so
    returns the grand total of |t| over every (R, orbital-pair), a large
    overestimate on realistic multi-hopping models (issue #57)."""
    try:
        import hwave.sc as sc

        # Even nk includes the zone boundary (k = pi); a coarse but tight
        # bound on the spectral range for the heuristic. The interaction adds
        # an extra spectral scale on top of the band.
        nk = 16
        kaxis = np.linspace(0.0, 2.0 * np.pi, nk, endpoint=False)
        eps_k = sc._build_hamiltonian_k(kaxis, kaxis, kaxis, hr, norb)
        # (norb, norb, Nx, Ny, Nz) -> (Nx, Ny, Nz, norb, norb); Hermitize
        # (guards against tiny asymmetry in the input) and diagonalize.
        hk = np.moveaxis(eps_k, (0, 1), (-2, -1))
        hk = 0.5 * (hk + np.conjugate(np.swapaxes(hk, -1, -2)))
        evals = np.linalg.eigvalsh(hk)
        if mu is None:
            mu = (sc._determine_mu(evals, beta, float(filling), norb)
                  if filling is not None else 0.0)
        u = 0.0
        for arr in inter_k.values():
            u = max(u, float(np.abs(np.asarray(arr)).max()))
    except Exception as exc:
        raise ValueError(
            "ir_wmax auto-estimate failed ({}); set [eliashberg] ir_wmax "
            "explicitly (a real-frequency bandwidth in the same energy "
            "units as the Hamiltonian).".format(exc))
    # The band-plus-interaction formula and its positivity check are the one
    # shared estimator (issue #184): ir_axis.auto_wmax, mu-aware, identical to
    # the FLEX side. Its ValueError already names [eliashberg] ir_wmax.
    from hwave.solver.ir_axis import auto_wmax
    return auto_wmax(evals, mu, u, param_hint="[eliashberg] ir_wmax")


def _ir_axes_for_run(eli_param, beta, hr, inter_k, norb, mu=None, filling=None):
    """Build the fermionic/bosonic IR axes for a dynamic-Eliashberg run."""
    from hwave.solver.ir_axis import IRAxis
    eps = float(eli_param.get("ir_tol", 1.0e-8))
    wmax = eli_param.get("ir_wmax")
    if wmax is None:
        wmax = _ir_auto_wmax(hr, inter_k, norb, beta, mu=mu, filling=filling)
        logger.info("IR: auto ir_wmax = %.6g (override with [eliashberg] "
                    "ir_wmax)", wmax)
    wmax = float(wmax)
    axF = IRAxis(beta=beta, wmax=wmax, eps=eps, statistics="F")
    axB = IRAxis(beta=beta, wmax=wmax, eps=eps, statistics="B")
    logger.info("IR: Lambda=%.3g eps=%.1e -> L_F=%d (nodes %d), L_B=%d "
                "(nodes %d)", beta * wmax, eps, axF.L, axF.n_freq,
                axB.L, axB.n_freq)
    return axF, axB


def _ir_compress(arr, ax, nmat, label, drop_constant=False,
                 keep_constant=False, error_on_large_constant=True,
                 max_chunk_bytes=1 << 28):
    """Fit a centered-uniform-grid array (..., nmat) to IR and return its
    values on the sparse frequency nodes (..., n_freq).

    The least-squares fit runs over ALL uniform frequencies (the uniform grid
    doubles as an oversampled residual check, design Sec. 3.2/5) and is
    applied in chunks so transient buffers stay bounded. The max residual
    over the full grid is logged; a residual large relative to the data
    scale warns with the remedy.

    ``drop_constant=True`` (used for the uniform-FFT susceptibilities)
    augments the fit with a frequency-independent constant. When it is the
    O(beta/Nmat) discretization artifact of the discrete tau -> i nu transform
    (a delta(tau) component the smooth IR basis cannot represent) it is small
    and DISCARDING it makes the IR representation closer to the continuum
    object than the raw uniform data.

    But when the susceptibility is *static-dominated* (large and nearly flat in
    nu within the sampled window -- the near-critical regime that matters for
    superconductivity), the constant column absorbs physical static weight, and
    dropping it silently corrupts the result (issue #57). A fitted constant
    that EXCEEDS the data scale cannot be the (small) O(beta/Nmat) artifact --
    it is the signature of an ill-conditioned fit on static-dominated data --
    so this raises ValueError unless ``keep_constant=True`` (retain the constant
    by adding it back onto every frequency node) or ``error_on_large_constant=
    False`` (drop it anyway -- for the kernel-algebra gate, which feeds both
    kernels the same densified data and asserts operator equivalence, not data
    fidelity). A constant above 5% of the data scale (but below it) still warns.
    The largest constant is logged.
    """
    lead = arr.reshape(-1, nmat)
    rows = max(1, int(max_chunk_bytes // max(1, arr.itemsize * nmat)))
    out = np.empty((lead.shape[0], ax.n_freq), dtype=np.complex128)
    fit_m, _ = ax.uniform_matrices(nmat, with_constant=drop_constant)
    resid = 0.0
    const_max = 0.0
    for s in range(0, lead.shape[0], rows):
        chunk = lead[s:s + rows]
        sol = chunk @ fit_m
        if drop_constant:
            coeffs = sol[..., :ax.L]
            const = sol[..., ax.L:ax.L + 1]
            const_max = max(const_max, float(np.abs(const).max()))
            resid = max(resid, float(np.abs(
                ax.eval_to_uniform(coeffs, nmat) + const - chunk).max()))
            node = ax.eval_to_freq(coeffs)
            if keep_constant:
                node = node + const
            out[s:s + rows] = node
        else:
            coeffs = sol
            resid = max(resid, float(
                np.abs(ax.eval_to_uniform(coeffs, nmat) - chunk).max()))
            out[s:s + rows] = ax.eval_to_freq(coeffs)
    scale = float(np.abs(lead).max()) or 1.0
    logger.info("IR compress %-8s: nmat=%d -> nodes=%d (L=%d), max uniform "
                "residual %.3e (rel %.3e)%s", label, nmat, ax.n_freq, ax.L,
                resid, resid / scale,
                ((", retained" if keep_constant else ", discarded")
                 + " frequency-independent constant %.3e" % const_max)
                if drop_constant else "")
    if drop_constant:
        if (not keep_constant and error_on_large_constant
                and const_max > scale):
            raise ValueError(
                "IR compress {}: the discarded frequency-independent component "
                "({:.3e}) exceeds the data scale ({:.3e}). It cannot be the "
                "small O(beta/Nmat) delta(tau) discretization artifact; this is "
                "an ill-conditioned fit on a static-dominated susceptibility, "
                "and dropping the constant gives an unusable result (issue "
                "#57). Lower [eliashberg] ir_wmax (Lambda = beta*wmax may be "
                "far too large), increase [mode.param] Nmat in the FLEX run, or "
                "set [eliashberg] ir_keep_static_chi = true to retain the "
                "static component.".format(label, const_max, scale))
        # Warn regardless of keep_constant/escape-hatch: a large constant is a
        # diagnostic signal (e.g. over-large ir_wmax) that must not be silenced
        # just because the caller opted to retain or tolerate it.
        if const_max > 0.05 * scale:
            logger.warning(
                "IR compress %s: the frequency-independent component "
                "(%.3e) is unusually large (>5%% of the data scale %.3e). The "
                "O(beta/Nmat) discretization constant should be small; a large "
                "value may indicate an unexpected constant in the input data -- "
                "check the FLEX output / increase [mode.param] Nmat in the FLEX "
                "run.%s", label, const_max, scale,
                " (retained via ir_keep_static_chi)" if keep_constant else "")
    if resid > 1.0e3 * ax.eps * scale:
        logger.warning(
            "IR fit residual for %s is large (rel %.3e > 1e3*ir_tol); the "
            "object may exceed the basis bandwidth -- raise ir_wmax or "
            "tighten ir_tol.", label, resid / scale)
    return out.reshape(arr.shape[:-1] + (ax.n_freq,))


def _ir_vertex_to_rtau(V_nodes, axB, axF, workers=1):
    """Pairing vertex on bosonic frequency nodes -> (r, fermionic tau nodes).

    The kernel's tau product lives on the FERMIONIC tau nodes (the result
    V*F is anti-periodic, so the fermionic fit applies); the bosonic
    coefficients are evaluated there exactly (design Sec. 3.2). This is the
    IR analogue of the hoisted ``Vs_rt`` invariant.
    """
    coeffs = axB.fit_from_freq(V_nodes)
    V_tau = axB.eval_to_tau_points(coeffs, axF.tau)
    return _spatial_ifftn(V_tau, axes=(4, 5, 6), workers=workers)


def _instantaneous_vertex(inter_k, norb, Nx, Ny, Nz, pairing_type,
                          convention, sc_matrices=None):
    """The frequency-INDEPENDENT part of the pairing vertex: the bare
    ``0.5*(S+C)``-type term of ``sc._compute_vertices_flex``, obtained by
    evaluating the vertex formula at chi_s = chi_c = 0. For a pure-Hubbard
    (CoulombIntra-only) model the Kuroki matrices have S = C = U on the
    intra-orbital element, so the SINGLET bare term ``0.5 (S + C)`` is the
    familiar +U while the TRIPLET one ``0.5 (C - S)`` vanishes; the triplet
    term is nonzero wherever the assembled matrices have C != S, which is
    the generic outcome of an inter-orbital U', a Hund / Ising term or an
    off-site density interaction (particular combinations can cancel, e.g.
    U' = J on the (ab,ab) element). On a frequency-EVEN pair amplitude the flat term is
    insensitive to the tau = 0 jump of F (its even-l IR coefficients vanish),
    which is why the issue-#57 defect stayed invisible on the shipped
    fixtures, whose gaps are flat in frequency.

    The Kuroki Exchange/PairHop rejection is enforced by the delegate
    itself (``sc._reject_reduced_flex_unsupported``), so this route is
    guarded even when called directly. ``sc_matrices`` forwards a
    precomputed (S, C) pair so the run's single build is reused."""
    import hwave.sc as sc
    zero = np.zeros((Nx, Ny, Nz, norb ** 2, norb ** 2), dtype=complex)
    return sc._compute_vertices_flex(zero, zero, inter_k, norb, Nx, Ny, Nz,
                                     pairing_type=pairing_type,
                                     convention=convention,
                                     sc_matrices=sc_matrices)


def eliashberg_kernel_ir(V_rt_tau, G2_nodes, phi_nodes, axF, beta,
                         V_inst_rt=None, workers=1):
    """Apply the dynamic Eliashberg kernel on sparse IR nodes.

    Mirrors ``eliashberg_kernel_dynamic`` with the phase-twisted FFT
    frequency transforms replaced by the IR node transforms (fused
    fit+evaluate matmuls). The IR transforms are PHYSICAL (G(tau) carries
    its 1/beta, the tau->freq step is the integral over tau), while the
    uniform-grid FFT chain carries one net factor beta; the explicit
    ``beta`` factor here restores the identical operator normalization
    (pinned by test_ir_matvec_matches_uniform_kernel).

    ``V_inst_rt`` (issue #57): the frequency-INDEPENDENT part of the
    pairing vertex (:func:`_instantaneous_vertex`), spatially transformed
    to r, shape (norb, norb, norb, norb, Nx, Ny, Nz). In imaginary time it
    is ``V_inst * delta(tau)`` -- OUT OF the bosonic IR basis, so it must
    never be fitted (``_ir_vertex_to_rtau`` would alias it into an
    uncontrolled smooth function); its tau integral is analytic instead:
    ``integral dtau e^{i w tau} V_inst delta(tau) F(tau) = V_inst F(0)``
    with ``F(0) = (1/beta) sum_nu F(i nu)`` -- the UNREGULARIZED Matsubara
    sum, i.e. the MIDPOINT ``0.5 * (F(0^+) + F(0^-))`` of the tau = 0 jump,
    evaluated exactly through the fermionic basis
    (``axF.u_matsubara_sum``).

    It must NOT be ``u_zero_plus``: that is the one-sided ``F(0^+)``, the
    midpoint PLUS half of the tau-jump. The jump half ANTI-commutes with
    the frequency reversal ``i w -> -i w`` (on IR coefficients the node
    reversal acts as ``c_l -> (-1)^(l+1) c_l``, and ``u_l(0^+)`` mixes both
    parities), so with ``u_zero_plus`` the kernel stops commuting with the
    combined parity operator as soon as the instantaneous vertex is nonzero
    -- i.e. for any model with a CoulombInter term -- and its flat term does
    not converge to the uniform-grid one at any Nmat. The earlier claim that
    "F ~ 1/nu^2, so the equal-time value is continuous" holds only for the
    particular smooth probes the shipped gates drive the kernel with; the
    kernel is a linear operator and must be right on every input.

    The uniform-grid kernel needs no such split: its dense tau grid
    represents the delta as a single bin, which IS the Matsubara sum
    truncated at Nmat.
    """
    xp = backend.array_module_of(G2_nodes)
    phi_nodes = xp.asarray(phi_nodes)
    F = xp.einsum('iljmxyzn,lmxyzn->ijxyzn', G2_nodes, phi_nodes)
    F_coeff = axF.fit_from_freq(F)
    F_rt = _spatial_ifftn(axF.eval_to_tau(F_coeff), axes=(2, 3, 4),
                          workers=workers)
    prod = -xp.einsum('abcdxyzt,bcxyzt->adxyzt', V_rt_tau, F_rt)
    out = axF.tau_to_freq(_spatial_fftn(prod, axes=(2, 3, 4),
                                        workers=workers))
    if V_inst_rt is not None:
        u0 = axF.u_matsubara_sum
        if xp is not np:
            u0 = xp.asarray(u0)
        F0_r = _spatial_ifftn(F_coeff @ u0, axes=(2, 3, 4), workers=workers)
        inst_r = -xp.einsum('abcdxyz,bcxyz->adxyz',
                            xp.asarray(V_inst_rt), F0_r)
        inst_k = _spatial_fftn(inst_r, axes=(2, 3, 4), workers=workers)
        out = out + inst_k[..., np.newaxis]
    return beta * out


def _load_seed_gap(eli_param, gap_shape, use_ir, axF, nmat):
    """Load an eigenvector-continuation seed from ``[eliashberg]
    seed_eigenvector`` (a ``gap_dynamic.npz`` written by a neighbouring run).

    Returns a flat complex vector matching ``gap_shape`` (C-order, the space
    the kernel operator acts on), or ``None`` when no seed is configured. The
    seed must be on the same uniform Matsubara grid and CellShape as this run
    (fail-fast otherwise); on the IR path the uniform-grid seed gap is refit
    onto the IR fermionic nodes so it lands in the same eigenvector space.
    """
    path = eli_param.get("seed_eigenvector")
    if not path:
        return None
    data = np.load(path)
    gap = np.asarray(data["gap"])   # (norb,norb,Nx,Ny,Nz,Nmat)
    # Fourier-sign provenance gate (issue #133): the seed's k labels must
    # match this run's convention; spatial axes are (2, 3, 4) BEFORE any
    # IR compression. Legacy unmarked seeds pass only if k-even.
    from hwave.solver.rpa import (validate_momentum_convention,
                                  check_momentum_marker)
    check_momentum_marker(data, path)   # unconditional (round-10)
    if gap.ndim != 6:
        raise ValueError(
            "seed_eigenvector '{}' has ndim = {} but a gap file is always "
            "(norb, norb, Nx, Ny, Nz, Nmat); the file is malformed or not "
            "a gap file.".format(path, gap.ndim))
    validate_momentum_convention(data, path, gap, (2, 3, 4),
                                 tuple(gap.shape[2:5]))
    if gap.shape[-1] != nmat:
        raise ValueError(
            "seed_eigenvector has Nmat={} but this run uses Nmat={}; "
            "eigenvector continuation requires the same uniform Matsubara grid "
            "and CellShape.".format(gap.shape[-1], nmat))
    if use_ir:
        gap = _ir_compress(gap, axF, nmat, "seed_gap", drop_constant=False)
    if gap.shape != gap_shape:
        raise ValueError(
            "seed_eigenvector shape {} does not match this run's {}."
            .format(gap.shape, gap_shape))
    return np.ascontiguousarray(gap).astype(complex).ravel()


def build_seed(eli_param, pairing_type, norb, kx, ky, kz, gap_shape, use_ir, axF, nmat):
    """The static init_gap form factor broadcast flat over the frequency
    axis (normalized), and the optional eigenvector-continuation seed.

    ``seed_vec`` is an optional eigenvector-continuation seed: a
    ``gap_dynamic.npz`` from a neighbouring run (e.g. the next
    temperature). Used as the ARPACK start vector AND to pick the
    eigenpair that overlaps it -- tracking one physical branch across an
    exceptional point of the non-Hermitian kernel, where the
    algebraically-largest eigenvalue can jump between a real and a
    complex branch. On the IR path the (uniform-grid) seed gap is refit
    onto the IR fermionic nodes so it lives in the same eigenvector
    space.
    """
    import hwave.sc as sc
    init_gap_mode = sc._resolve_init_gap(eli_param.get("init_gap"), pairing_type)
    sigma_static = sc._initialize_gap(init_gap_mode, norb, kx, ky, kz)
    phi0 = np.broadcast_to(sigma_static[..., np.newaxis], gap_shape).copy().astype(complex)
    n0 = np.linalg.norm(phi0)
    if n0 > 0:
        phi0 /= n0
    seed_vec = _load_seed_gap(eli_param, gap_shape, use_ir, axF, nmat)
    return phi0, seed_vec


def run_leading_eigenproblem(matvec, gap_shape, eli_param, pairing_type, *, phi0, seed_vec,
                             use_ir, axF, nmat, logger_label="dynamic kernel",
                             parity_leakage_policy="warn", parity_leakage_tol=1.0e-8):
    """Solve for the leading Eliashberg eigenpair of ``matvec`` and gauge-fix
    (and, on the IR path, densify) the resulting gap onto the uniform grid.

    Shared eigen-driver behind the on-site and bond-resolved dynamic
    solvers: the same ``scipy.sparse.linalg.LinearOperator`` +
    ``sc._solve_leading`` machinery, parity leakage probe / iterate
    projection / ARPACK eigenpair reordering, and gauge fix.

    ``parity_leakage_policy``: ``"warn"`` (default) reproduces the historical
    behaviour -- the cross-sector leakage probe runs only on the
    ``solver_mode="iteration"`` path, and a leaky kernel there falls back to
    an un-projected iteration with a logged warning. ``"refuse"`` runs the
    probe up front, before any solve, on every solver mode, and raises
    ``ValueError`` if the kernel does not commute with the channel's parity
    within tolerance.

    ``parity_leakage_tol`` (float >= 0) is that tolerance, the
    ``[eliashberg] parity_leakage_tol`` key of the bond entries. Its default
    ``1e-8`` is the historical hard-coded threshold, which is what the uniform
    grid reaches; the bond entries raise it to ``2e-2`` under
    ``matsubara_basis = "ir"``, where the IR representation of a uniform-FFT
    archive carries a parity asymmetry of its own decaying as ``Nmat^-2``.
    Under ``"refuse"`` a leakage in the band ``[0.1 * tol, tol)`` is accepted
    with a WARNING.

    Returns ``(lam, gap_w, eigenvalues_all, eigenvalue_match, eigenvalue_note,
    leakage, sector_weights, sector_selection, eigenvalue_selection)``, where
    ``eigenvalue_selection`` names the criterion that produced the reported
    leading eigenvalue (issue #202) -- ``"LR_projected"`` (a symmetry-valid
    projected solve on a conserved channel sector, the num_eigenvalues-
    independent Tc criterion), ``"LM"`` (plain arnoldi, largest magnitude),
    ``"LR"`` (user ``spectral_shift``, largest real part), ``"LR_retry"`` (the
    automatic largest-real re-solve when an unprojected ``"LM"`` set held no
    positive channel eigenvalue), ``"dense-LR"`` (the tiny-operator dense
    largest-real path), ``"shift-invert"``, or ``"iteration"`` -- and
    ``leakage`` is the measured
    combined-parity cross-sector leakage as a float, or ``None`` when no probe
    ran (the ``"warn"`` policy on the non-iteration solver modes),
    ``sector_weights`` is ``gap_sector_weights`` of the returned uniform-grid
    ``gap_w``, and ``sector_selection`` is the sector this run selected in --
    ``"channel"`` (the channel's even-frequency sector),
    ``"combined_parity"`` (the wider antisymmetry sector, used when the kernel
    does not conserve the frequency parity) or ``"none"``. On the iteration
    path it is the sector the iterates were projected onto; on the eigenvalue
    solver family, which never projects, it is the stage
    ``_reorder_eigenpairs_by_parity_dynamic`` matched with -- and therefore
    what the ``match`` column of ``eigenvalue.dat`` means.
    """
    from scipy.sparse.linalg import LinearOperator
    import hwave.sc as sc

    if parity_leakage_policy not in ("warn", "refuse"):
        raise ValueError("parity_leakage_policy must be 'warn' or 'refuse', got {!r}"
                         .format(parity_leakage_policy))
    try:
        tol_ok = parity_leakage_tol is not None \
            and np.isfinite(float(parity_leakage_tol)) and float(parity_leakage_tol) >= 0.0
    except (TypeError, ValueError):
        tol_ok = False
    if not tol_ok:
        raise ValueError("parity_leakage_tol must be a finite number >= 0, got {!r}"
                         .format(parity_leakage_tol))
    parity_leakage_tol = float(parity_leakage_tol)
    # the measured probe value, reported to the caller so the outputs can
    # record it; stays None when no probe runs on this path
    leakage = None

    vec_size = int(np.prod(gap_shape))
    solver_mode = eli_param.get("solver_mode", "iteration")
    eigenvalue_method = eli_param.get("eigenvalue_method", "arnoldi")
    num_eigenvalues = eli_param.get("num_eigenvalues", 10)
    max_iter = eli_param.get("max_iter", 1000)
    alpha = eli_param.get("alpha", 0.5)
    tol = eli_param.get("convergence_tol", 1.0e-5)

    def make_operator():
        return LinearOperator((vec_size, vec_size), matvec=matvec, dtype=complex), vec_size

    if parity_leakage_policy == "refuse":
        A_probe, _ = make_operator()
        leakage = _parity_leakage(A_probe, gap_shape, pairing_type)
        if leakage > parity_leakage_tol:
            raise ValueError(
                "the {} does not commute with the combined parity (cross-sector leakage "
                "{:.2e} > parity_leakage_tol = {:.2e}); the direct-term pairing kernel is the "
                "physical kernel only on a definite-parity subspace, so this run is refused. "
                "Likely causes: an asymmetric susceptibility archive, an IR fit error, or an "
                "inconsistent S/C and chi pair. Remedies: use matsubara_basis = 'uniform', "
                "tighten ir_tol / raise ir_wmax, regenerate the archive at a larger Nmat, or "
                "raise [eliashberg] parity_leakage_tol if this leakage is acceptable."
                .format(logger_label, leakage, parity_leakage_tol))
        if parity_leakage_tol > 0.0 and leakage >= 0.1 * parity_leakage_tol:
            logger.warning(
                "The %s commutes with the combined parity only to %.2e, within a decade of "
                "parity_leakage_tol = %.2e; the '%s' eigenvalue carries that much "
                "cross-sector contamination.",
                logger_label, leakage, parity_leakage_tol, pairing_type)

    # Map [eliashberg] controls to the _solve_leading solver_mode string,
    # exactly as calc_eliashberg does for the static path.
    eigenvalue_match = None
    # Mirrors sc.calc_eliashberg's eigenvalue_note (review fix I-2): when
    # spectral_shift is active on this power-iteration path, the value
    # written below is the SIGNED eigenvalue of the shifted-and-subtracted-
    # back kernel rather than the unshifted UNSIGNED iterate norm, and that
    # meaning change must be labelled in the output, not just logged.
    eigenvalue_note = None
    # The sector this run selected in: the one the power iteration projected
    # its iterates onto, or -- on the eigenvalue solver family, which never
    # projects -- the stage the eigenpair reordering matched with.
    sector_selection = "none"
    # Which selection criterion produced the reported leading eigenvalue, for
    # the outputs (issue #202): "LR_projected" (projected solve on a conserved
    # channel sector), "LM" (plain arnoldi, largest magnitude), "LR" (user
    # spectral_shift, largest real part), "LR_retry" (the automatic largest-real
    # re-solve when an unprojected LM set held no positive channel eigenvalue),
    # "dense-LR" (the tiny-operator dense path), "shift-invert", or "iteration".
    eigenvalue_selection = None
    if solver_mode == "iteration":
        # Mirror the static _solve_iteration: project every iterate onto the
        # requested sector so numerical noise cannot let the power iteration
        # drift into another one. A projection only yields an eigenpair of the
        # kernel when the sector it projects onto is an invariant subspace, so
        # probe before projecting, and use the widest sector that qualifies:
        #   channel (even-frequency) -> combined parity -> none.
        if leakage is None:              # not already measured by "refuse"
            A_probe, _ = make_operator()
            leakage = _parity_leakage(A_probe, gap_shape, pairing_type)
        A_probe, _ = make_operator()
        channel_leakage = _channel_leakage(A_probe, gap_shape, pairing_type)
        if channel_leakage is None:
            # the sector does not exist on this grid, so there is no gap to
            # iterate for; say why instead of failing later on an empty seed
            raise ValueError(
                "the '{}' even-frequency sector is empty on this k grid "
                "(every k is its own inverse); use a grid with at least one "
                "momentum axis of 3 or more points, or the singlet channel"
                .format(pairing_type))
        if _channel_projection_is_valid(channel_leakage, parity_leakage_tol):
            sector_selection = "channel"

            def project_fn(flat):
                return _project_parity_dynamic(
                    flat.reshape(gap_shape), pairing_type).ravel()
            # Project + normalize the seed, raising if it has no in-sector
            # component (matches sc._solve_iteration's guard).
            phi0 = _project_seed_dynamic(phi0, pairing_type)
        elif leakage <= parity_leakage_tol:
            sector_selection = "combined_parity"
            logger.warning(
                "frequency parity is not a symmetry of the kernel (channel "
                "leakage %.2e); the '%s' eigenvalue is the leading one of the "
                "combined-parity sector; see gap_sector_weights",
                channel_leakage, pairing_type)

            def project_fn(flat):
                return _project_combined_parity_dynamic(
                    flat.reshape(gap_shape), pairing_type).ravel()
            phi0 = _project_seed_combined_dynamic(phi0, pairing_type)
        else:
            sector_selection = "none"
            logger.warning(
                "Dynamic Eliashberg kernel does not commute with parity "
                "(cross-sector leakage %.2e); parity projection for the '%s' "
                "channel is disabled and the un-projected iteration is used.",
                leakage, pairing_type)
            project_fn = None
        # spectral_shift is honoured on the power-iteration path too: with a
        # repulsive-dominant kernel (negative dominant eigenvalue) the iterate
        # flips sign every step and the loop can never converge; iterating on
        # K + sigma*I fixes that and sc._solve_leading subtracts sigma back, so
        # the eigenvalue below is still the signed eigenvalue of K.
        iteration_spectral_shift = eli_param.get("spectral_shift")
        eigenvalue, sigma_flat, info = sc._solve_leading(
            make_operator, vec_size, "iteration",
            max_iter=max_iter, convergence_tol=tol, alpha=alpha,
            init_vec=phi0.ravel(), project_fn=project_fn,
            spectral_shift=iteration_spectral_shift)
        eigenvalues_all = None
        eigenvalue_selection = "iteration"
        if iteration_spectral_shift is not None:
            # Shifted power iteration: the value is the SIGNED eigenvalue of
            # the original dynamic kernel only when sc._solve_leading's
            # Rayleigh check validated it (<v|K|v>/<v|v> with a small
            # residual). Otherwise ||(K + sigma*I) v|| - sigma is NOT an
            # eigenvalue at all -- e.g. an insufficient sigma leaves the
            # dominant shifted eigenvalue negative -- and the shared note
            # labels it as an estimate so eigenvalue.dat cannot be misread.
            eigenvalue_note = sc._shifted_eigenvalue_note(
                "iteration", iteration_spectral_shift,
                info.get("converged"), info.get("n_iter"), info,
                kernel_label="dynamic kernel")
    else:
        # "eigenvalue" / "both": use the ARPACK/shift-invert eigen family.
        # Note: "both" degrades to eigenvalue-only here (the static path also
        # runs a power-iteration leg); the ARPACK leading pair is returned.
        if solver_mode == "both":
            logger.warning(
                "Dynamic solver_mode='both' runs the eigenvalue leg only; "
                "the power-iteration cross-check is skipped.")
        user_spectral_shift = eli_param.get("spectral_shift")

        # Issue #202: when the kernel CONSERVES the channel's even-frequency
        # sector, ARPACK's plain which='LM' (largest MAGNITUDE) set can miss the
        # channel's small positive Tc eigenvalue behind larger-real opposite-
        # sector modes, and the reported value would then depend on
        # num_eigenvalues. The fix is a symmetry-valid PROJECTED solve on that
        # sector (its largest-real eigenpair is num_eigenvalues-independent),
        # mirroring what the power iteration does with project_fn. Only the
        # plain arnoldi method benefits, and a seeded run tracks a continuation
        # branch on purpose, so neither shift-invert nor a seeded run projects.
        channel_valid = False
        if eigenvalue_method == "arnoldi" and seed_vec is None:
            A_probe, _ = make_operator()
            channel_leak = _channel_leakage(A_probe, gap_shape, pairing_type)
            channel_valid = _channel_projection_is_valid(
                channel_leak, parity_leakage_tol)

        if channel_valid:
            # Projected solve: P K P restricted to the channel sector, largest
            # real part. sector_selection is the channel by construction.
            eigenvalue, sigma_flat, info = _solve_channel_projected(
                make_operator, matvec, vec_size, gap_shape, pairing_type,
                num_eigenvalues)
            eigenvalues_all = info.get("eigenvalues")
            vecs_all = info.get("eigenvectors")
            if eigenvalues_all is not None and vecs_all is not None:
                eigenvalues_all, vecs_all, eigenvalue_match, sector_selection = \
                    _reorder_eigenpairs_by_parity_dynamic(
                        eigenvalues_all, vecs_all, gap_shape, pairing_type)
                eigenvalue = eigenvalues_all[0]
                sigma_flat = vecs_all[:, 0]
            eigenvalue_selection = "LR_projected"
        else:
            # Unprojected solve: the general multi-orbital case (the kernel does
            # not conserve the channel) or a seeded / shift-invert run. Here the
            # raw largest-real leading is already num_eigenvalues-independent,
            # and picking the first parity match is the documented best effort.
            # A plain which='LM' set that holds no positive channel eigenvalue
            # is re-solved ONCE with spectral_shift="auto" (which='LR'); a user
            # spectral_shift already used which='LR', and a seeded run is left
            # alone by design, so neither is re-solved.
            retry_eligible = (user_spectral_shift is None and seed_vec is None
                              and eigenvalue_method == "arnoldi")
            eigenvalue, sigma_flat, info = sc._solve_leading(
                make_operator, vec_size, eigenvalue_method,
                num_eigenvalues=num_eigenvalues,
                sigma_shift=eli_param.get("sigma_shift"),
                spectral_shift=user_spectral_shift,
                seed_vec=seed_vec,
                warn_negative_leading=not retry_eligible)
            eigenvalues_all = info.get("eigenvalues")
            vecs_all = info.get("eigenvectors")
            # Promote the eigenpair in the channel's sector so the reported
            # leading lambda is the physical singlet/triplet solution, not
            # ARPACK's raw largest-|lambda|. The stage that matched labels the
            # match column of eigenvalue.dat, so it travels with the outputs.
            if eigenvalues_all is not None and vecs_all is not None:
                eigenvalues_all, vecs_all, eigenvalue_match, sector_selection = \
                    _reorder_eigenpairs_by_parity_dynamic(
                        eigenvalues_all, vecs_all, gap_shape, pairing_type)
                eigenvalue = eigenvalues_all[0]
                sigma_flat = vecs_all[:, 0]
            # The criterion _solve_leading actually used (issue #202 review:
            # read it, do not infer -- the tiny-operator dense path is a
            # largest-real solve, not "LM").
            eigenvalue_selection = info.get("selection")

            if (eigenvalue_selection == "LM" and seed_vec is None
                    and eigenvalues_all is not None):
                vals_arr = np.asarray(eigenvalues_all)
                scale = float(np.max(np.abs(vals_arr))) if vals_arr.size else 0.0
                neg_tol = 1.0e-8 * max(1.0, scale)
                no_match = (sector_selection == "none"
                            or eigenvalue_match is None
                            or not bool(np.any(eigenvalue_match)))
                leading_re = float(np.real(eigenvalue))
                if no_match or leading_re < -neg_tol:
                    logger.warning(
                        "%s: arnoldi (which='LM', num_eigenvalues=%d) returned "
                        "no positive '%s' eigenvalue (leading Re(lambda) = "
                        "%.4g); re-solving with spectral_shift = \"auto\" "
                        "(largest real part) -- set [eliashberg] spectral_shift "
                        "= \"auto\" or use solver_mode = \"iteration\" to avoid "
                        "the extra solve",
                        logger_label, num_eigenvalues, pairing_type, leading_re)
                    eigenvalue, sigma_flat, info = sc._solve_leading(
                        make_operator, vec_size, eigenvalue_method,
                        num_eigenvalues=num_eigenvalues,
                        sigma_shift=eli_param.get("sigma_shift"),
                        spectral_shift="auto",
                        seed_vec=seed_vec)
                    eigenvalues_all = info.get("eigenvalues")
                    vecs_all = info.get("eigenvectors")
                    if eigenvalues_all is not None and vecs_all is not None:
                        eigenvalues_all, vecs_all, eigenvalue_match, \
                            sector_selection = \
                            _reorder_eigenpairs_by_parity_dynamic(
                                eigenvalues_all, vecs_all, gap_shape,
                                pairing_type)
                        eigenvalue = eigenvalues_all[0]
                        sigma_flat = vecs_all[:, 0]
                    eigenvalue_selection = "LR_retry"
                    retry_note = (
                        "eigenvalue selection: which='LM' with "
                        "num_eigenvalues={} returned no positive '{}' "
                        "eigenvalue (Re(lambda) = {:.4g}); the values below "
                        "come from a second solve with spectral_shift='auto' "
                        "(largest real part)".format(
                            num_eigenvalues, pairing_type, leading_re))
                    eigenvalue_note = (retry_note if not eigenvalue_note
                                       else eigenvalue_note + "\n" + retry_note)

    lam = float(np.real(eigenvalue))
    logger.info("%s leading eigenvalue lambda = %.6f", logger_label, lam)

    # --- Outputs ---
    # Gauge-fix the eigenvector (deterministic phase/normalization) so the
    # written gap is reproducible across runs and linear-algebra backends.
    gap_w = _fix_gauge(np.asarray(sigma_flat).reshape(gap_shape))
    if use_ir:
        # Densify the node-resolved gap back to the run's uniform grid so
        # the output format/metadata is IDENTICAL to the uniform path (the
        # gauge was fixed on nodes; re-fix after densification so the pivot
        # convention refers to the written array). The npz records the IR
        # provenance (design Sec. 3.2).
        gap_w = _fix_gauge(axF.eval_to_uniform(
            axF.fit_from_freq(gap_w), nmat))
    # Sector composition of what is actually returned: the solver only ever
    # SELECTS the channel's even-frequency sector (projected iteration) or
    # PREFERS it (eigenvalue reordering), so an impure result must be visible.
    sector_weights = gap_sector_weights(gap_w)
    channel_label = ("even_k_even_w" if pairing_type == "singlet"
                     else "odd_k_even_w")
    channel_weight = sector_weights[channel_label]
    logger.info("leading '%s' gap: %s weight %.6f",
                pairing_type, channel_label, channel_weight)
    if channel_weight < 0.999:
        logger.warning(
            "the returned gap is not a pure %s even-frequency state: %s weight "
            "%.6f (sector weights %s)", pairing_type, channel_label,
            channel_weight,
            "  ".join("{}={:.6f}".format(k, sector_weights[k])
                      for k in _SECTOR_LABELS))
    return (lam, gap_w, eigenvalues_all, eigenvalue_match, eigenvalue_note,
            None if leakage is None else float(leakage), sector_weights,
            sector_selection, eigenvalue_selection)


#: name of the per-eigenvalue match column, by the sector the run selected in.
#: ``None`` is the historical, unqualified label kept for callers that do not
#: say which sector they matched -- a ``1`` under it is ambiguous, which is
#: exactly what issue #209 removed from the solver.
_MATCH_COLUMN = {
    "channel": "match(1=channel even-frequency sector)",
    "combined_parity": "match(1=combined-parity sector; "
                       "no even-frequency eigenpair)",
    None: "match(1=channel-parity)",
}


def write_eigenvalue_file(path, lam, eigenvalues_all, eigenvalue_match, note,
                          header_lines=(), sector_weights=None,
                          selection=None, eigenvalue_selection=None):
    """Write the ``eigenvalue.dat`` leading-eigenvalue-and-spectrum file.

    ``header_lines`` are written as additional ``# ...`` lines right after
    the fixed first header line and before ``note`` -- e.g. the dynamic
    solver's ``zero_chi_s``/``zero_chi_c`` diagnostic-flag line.

    ``sector_weights`` (a ``gap_sector_weights`` dict) adds one more header
    line, ``# gap_sector_weights even_k_even_w=... ...``, right after them, and
    ``selection`` (``run_leading_eigenproblem``'s ``sector_selection``) adds
    ``# sector_selection=<channel|combined_parity|none>`` after that.
    ``eigenvalue_selection`` (issue #202), when not ``None``, adds
    ``# eigenvalue_selection: <LM|LR|LR_retry|shift-invert|iteration|subspace>``
    after that.

    ``selection`` also NAMES the per-eigenvalue ``match`` column, which means a
    different thing in each case: an even-frequency channel match, or only a
    combined-parity one because no even-frequency eigenpair existed. Callers
    that pass nothing keep the historical, unqualified label.
    """
    with open(path, "w") as fw:
        fw.write("# Dynamic Eliashberg leading eigenvalue\n")
        for line in header_lines:
            fw.write("# {}\n".format(line))
        if sector_weights is not None:
            fw.write("# gap_sector_weights {}\n".format(
                " ".join("{}={:.6f}".format(label,
                                            float(sector_weights[label]))
                         for label in _SECTOR_LABELS)))
        if selection is not None:
            fw.write("# sector_selection={}\n".format(selection))
        if eigenvalue_selection is not None:
            fw.write("# eigenvalue_selection: {}\n".format(eigenvalue_selection))
        if note:
            for line in str(note).splitlines():
                fw.write("# {}\n".format(line))
        fw.write("{:.8e}\n".format(lam))
        if eigenvalues_all is not None:
            if eigenvalue_match is not None:
                fw.write("# index  Re(eigenvalue)  Im(eigenvalue)  "
                         "|eigenvalue|  {}\n".format(_MATCH_COLUMN[selection]
                                                     if selection in _MATCH_COLUMN
                                                     else _MATCH_COLUMN[None]))
                for i, ev in enumerate(eigenvalues_all):
                    fw.write("{:4d} {:15.8e} {:15.8e} {:15.8e} {:d}\n".format(
                        i, ev.real, ev.imag, abs(ev),
                        int(bool(eigenvalue_match[i]))))
            else:
                fw.write("# index  Re(eigenvalue)  Im(eigenvalue)  "
                         "|eigenvalue|\n")
                for i, ev in enumerate(eigenvalues_all):
                    fw.write("{:4d} {:15.8e} {:15.8e} {:15.8e}\n".format(
                        i, ev.real, ev.imag, abs(ev)))


def solve_dynamic(input_dict):
    """Solve the dynamic (frequency-resolved) Eliashberg equation.

    Reads the FLEX outputs (full-frequency chi_s/chi_c and the dressed
    Green's function), builds the frequency-resolved pairing vertex and pair
    bubble, and finds the leading eigenpair of the tau-product Eliashberg
    kernel via the shared driver ``sc._solve_leading`` (the same eigenvalue
    ordering / shift-invert / iteration machinery as the static path).

    Parameters
    ----------
    input_dict : dict
        Parsed TOML configuration dictionary.

    Returns
    -------
    float
        The leading (largest real part) Eliashberg eigenvalue lambda.
    """
    import hwave.sc as sc

    # spin-orbital mode is unsupported here exactly as on the static path;
    # solve_dynamic is publicly callable, so guard both entries (issue #83)
    sc.reject_spin_orbital_mode(input_dict)
    mode_param = input_dict["mode"]["param"]
    T = mode_param["T"]
    # shared validated conversion (round-7 review): an unchecked 1/T here
    # let a subnormal T reach the dynamic solver as beta = inf
    beta = sc._coerce_run_beta(T)
    cell_shape = mode_param["CellShape"]
    # Resolve SubShape by the PACKAGE convention (documented default:
    # CellShape, i.e. the whole cell as one supercell) and guard the
    # RESOLVED value: a guard keyed on the explicit key alone let an
    # omitted SubShape default to a fully folded configuration and reach
    # the file reader with the very mismatch the guard exists to stop
    # (round-7 review).
    if isinstance(cell_shape, (list, tuple)):
        _cs = list(cell_shape)
    else:
        _cs = [cell_shape]
    while len(_cs) < 3:
        _cs.append(1)
    _ss = list(mode_param.get("SubShape", _cs))
    while len(_ss) < 3:
        _ss.append(1)
    if _ss != [1, 1, 1]:
        # supported nowhere in this module: the geometry and
        # interactions are consumed UNFOLDED here, so a folded
        # susceptibility mismatches the expected orbital count and an
        # off-site bond would fold onto an on-site supercell entry
        # (round-4 review); failing late produced an unhelpful shape
        # error instead of this actionable one
        raise ValueError(
            "SubShape (sublattice folding) is not supported by the "
            "Eliashberg module: fold the model into the unit cell "
            "yourself, or set SubShape = [1, 1, 1] explicitly. (Note: "
            "omitting SubShape defaults it to CellShape -- the whole "
            "cell as one supercell -- per the package convention, so it "
            "must be set explicitly here.)")
    sub_shape = _ss
    if isinstance(cell_shape, list):
        cell_shape = list(cell_shape)
        while len(cell_shape) < 3:
            cell_shape.append(1)
    Lx, Ly, Lz = cell_shape
    if isinstance(sub_shape, list):
        sub_shape = list(sub_shape)
        while len(sub_shape) < 3:
            sub_shape.append(1)
    Bx, By, Bz = sub_shape
    Nx, Ny, Nz = Lx // Bx, Ly // By, Lz // Bz
    Nk = Nx * Ny * Nz

    eli_param = input_dict.get("eliashberg", {})
    # bond_channels = true selects the bond-resolved pairing kernel, which is a
    # different vertex entirely (spec 6). sc.calc_eliashberg dispatches to it
    # FIRST, so that _validate_dynamic_prereqs runs before anything is read;
    # this second, equivalent route only covers a DIRECT solve_dynamic() call,
    # which must not be able to fall through to the scalar on-site vertex with
    # the flag silently ignored. The flag goes through the same strict reader as
    # every other bond option, so a typo is refused rather than read as false.
    if sc._bond_bool_option(eli_param, "bond_channels", False):
        from hwave.solver import eliashberg_bond_io
        return eliashberg_bond_io.solve_dynamic_bond(input_dict)
    pairing_type = eli_param.get("pairing_type", "singlet")
    use_gpu = _gpu_requested(eli_param)
    xp, gpu_active = backend.get_backend(
        use_gpu, logger=logger, required=_gpu_required_requested(eli_param))
    matsubara_basis = str(
        eli_param.get("matsubara_basis", "uniform")).lower()
    if matsubara_basis not in ("uniform", "ir"):
        raise ValueError(
            "matsubara_basis must be 'uniform' or 'ir', got '{}'."
            .format(matsubara_basis))
    use_ir = (matsubara_basis == "ir")
    zero_chi_c = backend.as_bool(eli_param.get("zero_chi_c", False))
    zero_chi_s = backend.as_bool(eli_param.get("zero_chi_s", False))

    # --- Geometry / interactions (norb from the geometry file) ---
    geom_info, hr, interactions = sc._read_interaction_files(input_dict)
    norb = geom_info["norb"]

    kx_array = np.linspace(0, 2.0 * np.pi, Nx, endpoint=False)
    ky_array = np.linspace(0, 2.0 * np.pi, Ny, endpoint=False)
    kz_array = np.linspace(0, 2.0 * np.pi, Nz, endpoint=False)
    inter_k = sc._build_interaction_k(kx_array, ky_array, kz_array,
                                      interactions, norb)

    # --- FLEX inputs (full frequency) ---
    if use_ir:
        chis_w, chic_w, green_w, chi_convention, ir_file_meta = \
            load_flex_chi_dynamic(input_dict, norb, Nx, Ny, Nz,
                                  allow_ir=True, interactions=interactions)
    else:
        chis_w, chic_w, green_w, chi_convention = load_flex_chi_dynamic(
            input_dict, norb, Nx, Ny, Nz, interactions=interactions)
        ir_file_meta = None
    if green_w is None:
        raise ValueError(
            "dynamic Eliashberg requires the dressed green.npz from the FLEX "
            "run (the pair bubble G2 is built from it); none was found. Check "
            "[file.input] path_to_flex_output / [eliashberg] flex_green.")
    nmat = chis_w.shape[-1]

    # --- IR frequency axis (design Sec. 3.2 / Stage-3 Sec. 4.1):
    # everything downstream (vertex assembly, pair bubble, kernel, parity
    # machinery) is per-frequency or reversal-based, so it operates on the
    # sparse symmetric node axis unchanged. The full uniform tensors of the
    # VERTEX and G2 are never built on the IR path.
    axF = axB = None
    if use_ir:
        axF, axB = _ir_axes_for_run(eli_param, beta, hr, inter_k, norb,
                                    mu=mode_param.get("mu"),
                                    filling=mode_param.get("filling"))
        if ir_file_meta is not None:
            # IR-native inputs (Stage 3): the files already hold node
            # values; refit each onto the run axes (pass-through when the
            # node sets coincide). No drop_constant -- node values carry no
            # uniform-FFT delta(tau) artifact.
            if zero_chi_s:
                _ir_validate_native_nodes(
                    chis_w, ir_file_meta["chis"], axB, "chiq_s", beta
                )
                chis_w = np.zeros(chis_w.shape[:-1] + (axB.n_freq,), dtype=chis_w.dtype)
            else:
                chis_w = _ir_refit_nodes(
                    chis_w, ir_file_meta["chis"], axB, "chiq_s", beta
                )
            if zero_chi_c:
                _ir_validate_native_nodes(
                    chic_w, ir_file_meta["chic"], axB, "chiq_c", beta
                )
                chic_w = np.zeros(chic_w.shape[:-1] + (axB.n_freq,), dtype=chic_w.dtype)
            else:
                chic_w = _ir_refit_nodes(
                    chic_w, ir_file_meta["chic"], axB, "chiq_c", beta
                )
            green_w = _ir_refit_nodes(green_w, ir_file_meta["green"], axF,
                                      "green", beta)
            # the uniform grid exists only as the OUTPUT grid here
            nmat = int(input_dict["mode"]["param"].get("Nmat", 1024))
        else:
            keep_static = _ir_keep_static_requested(eli_param)
            if zero_chi_s:
                chis_w = np.zeros(chis_w.shape[:-1] + (axB.n_freq,), dtype=chis_w.dtype)
            else:
                chis_w = _ir_compress(
                    chis_w,
                    axB,
                    nmat,
                    "chiq_s",
                    drop_constant=True,
                    keep_constant=keep_static,
                )
            if zero_chi_c:
                chic_w = np.zeros(chic_w.shape[:-1] + (axB.n_freq,), dtype=chic_w.dtype)
            else:
                chic_w = _ir_compress(
                    chic_w,
                    axB,
                    nmat,
                    "chiq_c",
                    drop_constant=True,
                    keep_constant=keep_static,
                )
            green_w = _ir_compress(green_w, axF, nmat, "green")
    nfreq_axis = axF.n_freq if use_ir else nmat

    # --- Diagnostic: optionally zero one fluctuation channel to decompose the
    #     pairing vertex into its spin (chi_s) and charge (chi_c) contributions.
    #     Both the singlet V = 1.5 S.chi_s.S - 0.5 C.chi_c.C + 0.5(S+C) and the
    #     triplet V = -0.5 S.chi_s.S - 0.5 C.chi_c.C + 0.5(C-S) vertices are
    #     linear in chi_s, chi_c, so this works for either pairing_type. Both
    #     flags default off, so the production vertex is unchanged. NOTE: the
    #     instantaneous bare term is retained in every case, and the linearized-gap
    #     eigenvalue problem is nonlinear in the vertex, so eigenvalues from
    #     separately zeroed runs are NOT additive
    #     (lambda_spin + lambda_charge != lambda_full in general).
    #     Booleans coerced via backend.as_bool (as for the gpu/ir flags) so a
    #     programmatic string "false" does not silently enable the diagnostic.
    if zero_chi_c and zero_chi_s:
        logger.warning(
            "zero_chi_c=zero_chi_s=True: both susceptibilities "
            "zeroed; bare (instantaneous) vertex only (diagnostic)."
        )
        chic_w[...] = 0
        chis_w[...] = 0
    else:
        if zero_chi_c:
            logger.warning(
                "zero_chi_c=True: charge susceptibility zeroed in "
                "the pairing vertex (spin+bare channel; diagnostic)."
            )
            chic_w[...] = 0
        if zero_chi_s:
            logger.warning(
                "zero_chi_s=True: spin susceptibility zeroed in the "
                "pairing vertex (charge+bare channel; diagnostic)."
            )
            chis_w[...] = 0

    # --- Vertex and pair bubble on the frequency axis ---
    logger.info("Computing dynamic FLEX pairing vertex (pairing_type=%s, "
                "convention=%s)...", pairing_type, chi_convention)
    # Reject before allocating the O(Nq * norb^4) pair (round-10 review).
    sc._reject_reduced_flex_unsupported(inter_k, chi_convention)
    # One S/C build for the whole solve: reused by every per-frequency
    # contraction and by the IR instantaneous vertex below (round-9 review).
    sc_mats = sc._build_vertex_sc_matrices(chi_convention, inter_k,
                                           norb, Nx, Ny, Nz)
    Vs_q_w = compute_vertices_flex_dynamic(
        chis_w, chic_w, inter_k, norb, Nx, Ny, Nz,
        pairing_type=pairing_type, convention=chi_convention,
        sc_matrices=sc_mats)
    logger.info("Computing frequency-resolved pair bubble G2...")
    G2_w = calc_g2_dynamic(green_w, beta)

    # --- Seed: static init_gap form factor, broadcast flat across omega ---
    gap_shape = (norb, norb, Nx, Ny, Nz, nfreq_axis)
    phi0, seed_vec = build_seed(eli_param, pairing_type, norb, kx_array, ky_array,
                                kz_array, gap_shape, use_ir, axF, nmat)

    vec_size = norb * norb * Nk * nfreq_axis
    assert phi0.size == vec_size

    # GPU path: park the two large invariants (pair bubble and vertex) on the
    # device once; every matvec then only moves the gap vector across PCIe.
    if gpu_active:
        logger.info("GPU backend active (CuPy): moving G2 and the pairing "
                    "vertex to the device (%.2f GB each).", G2_w.nbytes / 1e9)
        # Two resident tensors plus roughly one same-sized transform
        # workspace per matvec (the gap-sized arrays are norb^2 smaller).
        backend.warn_if_device_memory_short(
            3 * G2_w.nbytes, logger, label="the dynamic Eliashberg kernel")
        G2_w = xp.asarray(G2_w)
        Vs_q_w = xp.asarray(Vs_q_w)

    # Spatial-FFT parallelism for the CPU kernel (scipy.fft workers): the
    # default 1 keeps the serial numpy path (bit-compatible with previous
    # releases); -1 uses all cores. Opt-in so existing runs are unchanged and
    # concurrent solves do not oversubscribe against OMP/MKL threads. Ignored
    # on the GPU backend (cuFFT already runs on the device).
    fft_workers = eli_param.get("fft_workers", 1)

    # The vertex's (q, i nu) -> (r, tau) transform is phi-independent and
    # dominates the matvec cost, so do it once here; drop the (q, i nu) form
    # to keep the resident vertex memory unchanged. On the IR path the tau
    # grid is the fermionic node set (the product V*F is anti-periodic).
    V_inst_rt = None
    if use_ir:
        # Issue #57: split off the frequency-INDEPENDENT (bare 0.5*(S+C))
        # part of the vertex BEFORE the bosonic-basis fit -- in tau it is a
        # delta(tau), out of any IR basis, and fitting it aliases it into
        # an uncontrolled smooth function. The kernel handles it
        # analytically (see eliashberg_kernel_ir).
        V_inst = _instantaneous_vertex(inter_k, norb, Nx, Ny, Nz,
                                       pairing_type=pairing_type,
                                       convention=chi_convention,
                                       sc_matrices=sc_mats)
        inst_scale = float(np.abs(V_inst).max())
        if inst_scale > 0.0:
            logger.info("IR: instantaneous vertex part split off "
                        "analytically (max |V_inst| = %.6g).", inst_scale)
            # xp.asarray: on the GPU path Vs_q_w is already a device array
            # (moved above), while V_inst is host-built -- the subtraction
            # must not mix backends. Plain no-op cast on numpy.
            Vs_q_w = Vs_q_w - xp.asarray(V_inst)[..., np.newaxis]
            V_inst_rt = _spatial_ifftn(V_inst.astype(complex),
                                       axes=(4, 5, 6), workers=fft_workers)
            if gpu_active:
                V_inst_rt = xp.asarray(V_inst_rt)
        Vs_rt = _ir_vertex_to_rtau(Vs_q_w, axB, axF, workers=fft_workers)
    else:
        Vs_rt = vertex_qw_to_rt(Vs_q_w, workers=fft_workers)
    del Vs_q_w

    def _matvec(x):
        if use_ir:
            out = eliashberg_kernel_ir(
                Vs_rt, G2_w, x.reshape(gap_shape), axF, beta,
                V_inst_rt=V_inst_rt, workers=fft_workers)
        else:
            out = eliashberg_kernel_dynamic(
                None, G2_w, x.reshape(gap_shape), norb, beta, Vs_rt=Vs_rt,
                workers=fft_workers)
        return backend.to_host(out).ravel()

    # the on-site path keeps the historical "warn" policy and its default
    # tolerance; the measured leakage is not part of its output format
    lam, gap_w, eigenvalues_all, eigenvalue_match, dynamic_eigenvalue_note, \
        _leakage, sector_weights, sector_selection, eigenvalue_selection = \
        run_leading_eigenproblem(
            _matvec, gap_shape, eli_param, pairing_type, phi0=phi0,
            seed_vec=seed_vec, use_ir=use_ir, axF=axF, nmat=nmat,
            logger_label="Dynamic Eliashberg")

    # --- Outputs ---
    output_dir = input_dict["file"]["output"]["path_to_output"]
    os.makedirs(output_dir, exist_ok=True)
    eigenvalue_file = eli_param.get("output_eigenvalue", "eigenvalue.dat")
    write_eigenvalue_file(
        os.path.join(output_dir, eigenvalue_file), lam, eigenvalues_all,
        eigenvalue_match, dynamic_eigenvalue_note,
        header_lines=(
            ["zero_chi_s={}  zero_chi_c={}".format(
                str(zero_chi_s).lower(), str(zero_chi_c).lower())]
            if (zero_chi_s or zero_chi_c) else []),
        sector_weights=sector_weights, selection=sector_selection,
        eigenvalue_selection=eigenvalue_selection)

    gap_file = eli_param.get("output_gap", "gap.dat")
    # Provenance metadata is added ONLY on the opt-in IR path: the default
    # uniform output keeps its exact historical key set.
    extra_meta = {}
    if use_ir:
        extra_meta.update(
            {
                "matsubara_basis": "ir",
                "ir_tol": axF.eps,
                "ir_wmax": axF.wmax,
                "ir_L": axF.L,
            }
        )
    if zero_chi_s or zero_chi_c:
        extra_meta.update({"zero_chi_s": zero_chi_s, "zero_chi_c": zero_chi_c})
    write_dynamic_outputs(
        output_dir,
        gap_w,
        lam,
        T,
        pairing_type,
        kx_array,
        ky_array,
        kz_array,
        beta,
        gap_file=gap_file,
        extra_meta=extra_meta or None,
        sector_weights=sector_weights,
        selection=sector_selection,
        eigenvalue_selection=eigenvalue_selection,
    )

    return lam
