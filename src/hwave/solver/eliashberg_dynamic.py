"""Dynamic Eliashberg equation solver (frequency-resolved).

This module implements the full-frequency Eliashberg equation solver for
analyzing superconducting instabilities with frequency-dependent kernels.
"""

import os
import logging
import numpy as np

from hwave.solver import backend
from hwave.solver import matsubara as ms

# Shared spatial-FFT helpers (scipy-parallel on CPU, cuFFT on GPU) live in
# backend.py so RPA/FLEX use the same implementations; keep the module-local
# names used throughout this file and by the tests.
from hwave.solver.backend import (_SFFT,                       # noqa: F401
                                  spatial_fftn as _spatial_fftn,
                                  spatial_ifftn as _spatial_ifftn)

logger = logging.getLogger("qlms").getChild("eliashberg_dynamic")


def _gpu_requested(eli_param):
    """Whether ``[eliashberg] gpu`` requests GPU execution.

    Accepts a real TOML boolean, and (for programmatic dict configs) a truthy
    string such as ``"true"``/``"1"``/``"yes"`` so a stray ``"false"`` string is
    not treated as True by ``bool("false")``. Delegates to the shared
    ``backend.as_bool`` (also used by the RPA/FLEX/UHF gpu-flag reads).
    """
    return backend.as_bool(eli_param.get("gpu", False))


# ---------------------------------------------------------------------------
# Channel-parity of the frequency-resolved gap
#
# The static path filters eigenpairs by the (k, orbital) parity operator P:
# Delta_{ab}(k) -> Delta_{ba}(-k) (sc._reverse_k_and_orbital), so a singlet
# request never reports an odd-parity (triplet-sector) mode. The dynamic gap
# carries an extra fermionic Matsubara axis, and fermion antisymmetry makes the
# physical channel-parity the COMBINED operator
#     P: Delta_{ab}(k, iw_n) -> Delta_{ba}(-k, -iw_n),
# even (+) for singlet, odd (-) for triplet. This admits both a conventional
# (even-k, even-w) and an odd-frequency (odd-k, odd-w) singlet, and rejects the
# mixed even-w/odd-k mode that ARPACK's raw largest-|lambda| can surface.
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
        FFT-grid index map ``i -> (N - i) % N`` (reverse + roll-by-1, matching
        ``sc._reverse_k_and_orbital``); the centered fermionic Matsubara
        partner of index ``n`` is ``nmat - 1 - n`` (a plain reversal, no roll);
        the two orbital indices are swapped.
    """
    rev = gap_w[:, :, ::-1, ::-1, ::-1, :]
    rev = np.roll(rev, 1, axis=(2, 3, 4))   # k -> -k on the FFT grid
    rev = rev[..., ::-1]                     # iw_n -> -iw_n (centered fermionic)
    rev = np.swapaxes(rev, 0, 1)             # orbital transpose a <-> b
    return rev


def _project_parity_dynamic(gap_w, pairing_type):
    """Project a dynamic gap onto the channel's combined-parity sector.

    ``(1 +/- P)/2`` with ``+`` for singlet (even) and ``-`` for triplet (odd).
    """
    if pairing_type == "singlet":
        sign = 1.0
    elif pairing_type == "triplet":
        sign = -1.0
    else:
        raise ValueError(
            "Unknown pairing_type: '{}'. Use 'singlet' or 'triplet'.".format(
                pairing_type))
    return 0.5 * (gap_w + sign * _reverse_kw_and_orbital(gap_w))


def _is_parity_dynamic(gap_w, pairing_type, tol=0.9):
    """True when ``gap_w`` retains at least ``tol`` of its norm under the
    channel's combined-parity projection (mirrors ``sc._is_gap_parity``)."""
    proj = _project_parity_dynamic(gap_w, pairing_type)
    n = np.linalg.norm(gap_w)
    if n == 0:
        return False
    return np.linalg.norm(proj) / n >= tol


def _parity_leakage(A, gap_shape, pairing_type, n_probe=None, seed=0):
    """Max fraction of ``A x`` that lands in the OPPOSITE parity sector.

    Zero (to numerical precision) iff the kernel commutes with the parity
    operator, so parity projection of the power iteration is legitimate.
    Mirrors ``sc._solve_iteration``'s centrosymmetry guard verbatim: for each
    of ``n_probe`` random vectors and each parity sign, project the probe into
    that sector, apply ``A``, and measure the norm fraction of the image that
    lands in the opposite sector (denominator ``||A xp|| + 1e-300``). Uses the
    same probe count as the static path (``sc._PARITY_GUARD_PROBES``).
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
        for sign, pt in ((1.0, "singlet"), (-1.0, "triplet")):
            xp = _project_parity_dynamic(x, pt)          # even (+) / odd (-)
            axp = A.matvec(xp.ravel()).reshape(gap_shape)
            opp = "triplet" if pt == "singlet" else "singlet"
            denom = np.linalg.norm(axp) + 1.0e-300
            leakage = max(leakage,
                          np.linalg.norm(_project_parity_dynamic(axp, opp))
                          / denom)
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
            "Initial gap has no component in the '{}' parity sector; choose "
            "an init_gap of the matching parity (even for singlet, odd for "
            "triplet).".format(pairing_type))
    return proj


def _reorder_eigenpairs_by_parity_dynamic(vals, vecs, gap_shape, pairing_type):
    """Promote eigenpairs whose gap has the channel's combined parity.

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

    Returns
    -------
    vals, vecs, match : ndarray
        Reordered with channel-parity eigenpairs first (order preserved within
        each group), plus the reordered boolean match array.
    """
    match = np.array([_is_parity_dynamic(vecs[:, i].reshape(gap_shape),
                                         pairing_type)
                      for i in range(vecs.shape[1])])
    if not np.any(match):
        logger.warning(
            "No dynamic eigenpair matches the requested '%s' parity; the "
            "reported leading gap belongs to the other channel. Increase "
            "num_eigenvalues or check pairing_type.", pairing_type)
    idx = np.concatenate([np.where(match)[0], np.where(~match)[0]])
    return vals[idx], vecs[:, idx], match[idx]


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
    from numpy.lib import format as _npformat
    # Best-effort header probe. This is a pure optimization/safety pre-check --
    # the loader below is the authoritative path -- so ANY failure returns None
    # and lets the loader raise the existing, clearer error. The broad catch is
    # deliberate: besides unreadable/malformed/truncated files and a missing
    # axis, it also covers the numpy.lib.format internals used here being
    # renamed/removed in a future numpy (AttributeError), which must degrade
    # gracefully rather than crash.
    try:
        with zipfile.ZipFile(path) as z:
            names = set(z.namelist())
            for key in keys:
                member = key + ".npy"
                if member in names:
                    with z.open(member) as f:
                        version = _npformat.read_magic(f)
                        _npformat._check_version(version)
                        shape, _fortran, _dtype = _npformat._read_array_header(
                            f, version)
                    return int(shape[axis])
    except Exception:
        return None
    return None


def load_flex_chi_dynamic(input_dict, norb, Nx, Ny, Nz):
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
    stored = {"chis": _npz_freq_size(chi_s_path, ("chiq_s", "chiq"), axis=0),
              "chic": _npz_freq_size(chi_c_path, ("chiq_c", "chiq"), axis=0)}
    if os.path.exists(green_path):
        stored["green"] = _npz_freq_size(green_path, ("green",), axis=1)
    for name, n in stored.items():
        if n is not None and n != cfg_nmat:
            raise ValueError(
                "dynamic Eliashberg grid mismatch: {} nmat={} differs from "
                "config Nmat={}".format(name, n, cfg_nmat))
    # sizes agree with the config here, but guard on the stored value so the
    # estimate always reflects what is on disk.
    file_nmat = max([n for n in stored.values() if n is not None]
                    + [cfg_nmat])
    check_memory(norb, Nx * Ny * Nz, file_nmat,
                 input_dict["eliashberg"].get("mem_limit_gb"))

    # reuse the exact static path/name/convention/spin-orbital-expansion logic
    chis_w, chic_w, green_w, chi_convention = \
        sc._load_flex_susceptibilities_full(input_dict, norb, Nx, Ny, Nz)
    # Belt-and-suspenders: the header check above already rejected a mismatch,
    # but keep a post-load assertion in case the loader reshapes unexpectedly.
    green_nmat = green_w.shape[-1] if green_w is not None else cfg_nmat
    if not (chis_w.shape[-1] == chic_w.shape[-1] == green_nmat == cfg_nmat):
        raise ValueError(
            "dynamic Eliashberg grid mismatch: nmat differs — chis={}, chic={}, "
            "green={}, config Nmat={}".format(
                chis_w.shape[-1], chic_w.shape[-1], green_nmat, cfg_nmat))
    return chis_w, chic_w, green_w, chi_convention


def compute_vertices_flex_dynamic(chis_w, chic_w, inter_k, norb,
                                  Nx, Ny, Nz, pairing_type, convention):
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

    def _one(l):
        return sc._compute_vertices_flex(
            chis_w[..., l], chic_w[..., l], inter_k, norb, Nx, Ny, Nz,
            pairing_type=pairing_type, convention=convention)

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
    G(-k,-wn) built via roll+flip. This function drops the sum over n and
    returns the per-frequency summand, so calc_g2_dynamic(...).sum(axis=-1)
    reproduces sc._calc_g2(...) to machine precision (see
    tests/test_eliashberg_dynamic.py::test_g2_dynamic_sums_to_static).

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

    # G(-k, -wn) via roll+flip -- SAME construction as sc._calc_g2.
    green_kw_inv = np.roll(
        green_kw[:, :, ::-1, ::-1, ::-1, ::-1],
        (1, 1, 1), (2, 3, 4)
    )
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
                          extra_meta=None):
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
        **(extra_meta or {}),
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
        band = float(np.abs(evals - float(mu)).max())
        u = 0.0
        for arr in inter_k.values():
            u = max(u, float(np.abs(np.asarray(arr)).max()))
        wmax = 3.0 * (band + u)
    except Exception as exc:
        raise ValueError(
            "ir_wmax auto-estimate failed ({}); set [eliashberg] ir_wmax "
            "explicitly (a real-frequency bandwidth in the same energy "
            "units as the Hamiltonian).".format(exc))
    if not np.isfinite(wmax) or wmax <= 0.0:
        raise ValueError(
            "ir_wmax auto-estimate is not a positive finite number "
            "(spectral range + interaction scale gave {}); set [eliashberg] "
            "ir_wmax explicitly.".format(wmax))
    return wmax


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
                          convention):
    """The frequency-INDEPENDENT part of the pairing vertex: the bare
    ``0.5*(S+C)``-type term of ``sc._compute_vertices_flex``, obtained by
    evaluating the vertex formula at chi_s = chi_c = 0. For a pure-Hubbard
    (CoulombIntra-only) model in the Kuroki convention this cancels
    exactly (S+C = U-U = 0), which is why the issue-#57 defect stayed
    invisible on the CoulombIntra-only test fixtures."""
    import hwave.sc as sc
    zero = np.zeros((Nx, Ny, Nz, norb ** 2, norb ** 2), dtype=complex)
    return sc._compute_vertices_flex(zero, zero, inter_k, norb, Nx, Ny, Nz,
                                     pairing_type=pairing_type,
                                     convention=convention)


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
    with ``F(0) = (1/beta) sum_nu F(i nu)`` evaluated exactly through the
    fermionic basis (``u_zero_plus``; F ~ 1/nu^2, so the equal-time value
    is continuous and needs no 0^+ regularization). The uniform-grid
    kernel needs no such split: its dense tau grid represents the delta as
    a single bin.
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
        u0 = axF.u_zero_plus
        if xp is not np:
            u0 = xp.asarray(u0)
        F0_r = _spatial_ifftn(F_coeff @ u0, axes=(2, 3, 4), workers=workers)
        inst_r = -xp.einsum('abcdxyz,bcxyz->adxyz',
                            xp.asarray(V_inst_rt), F0_r)
        inst_k = _spatial_fftn(inst_r, axes=(2, 3, 4), workers=workers)
        out = out + inst_k[..., np.newaxis]
    return beta * out


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
    from scipy.sparse.linalg import LinearOperator

    mode_param = input_dict["mode"]["param"]
    T = mode_param["T"]
    beta = 1.0 / T
    cell_shape = mode_param["CellShape"]
    sub_shape = mode_param.get("SubShape", cell_shape)
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
    pairing_type = eli_param.get("pairing_type", "singlet")
    # Mirror calc_eliashberg's config -> _solve_leading string mapping.
    solver_mode = eli_param.get("solver_mode", "iteration")
    eigenvalue_method = eli_param.get("eigenvalue_method", "arnoldi")
    num_eigenvalues = eli_param.get("num_eigenvalues", 10)
    max_iter = eli_param.get("max_iter", 1000)
    alpha = eli_param.get("alpha", 0.5)
    tol = eli_param.get("convergence_tol", 1.0e-5)
    init_gap_mode = sc._resolve_init_gap(eli_param.get("init_gap"), pairing_type)
    use_gpu = _gpu_requested(eli_param)
    xp, gpu_active = backend.get_backend(use_gpu, logger=logger)
    matsubara_basis = str(
        eli_param.get("matsubara_basis", "uniform")).lower()
    if matsubara_basis not in ("uniform", "ir"):
        raise ValueError(
            "matsubara_basis must be 'uniform' or 'ir', got '{}'."
            .format(matsubara_basis))
    use_ir = (matsubara_basis == "ir")

    # --- Geometry / interactions (norb from the geometry file) ---
    geom_info, hr, interactions = sc._read_interaction_files(input_dict)
    norb = geom_info["norb"]

    kx_array = np.linspace(0, 2.0 * np.pi, Nx, endpoint=False)
    ky_array = np.linspace(0, 2.0 * np.pi, Ny, endpoint=False)
    kz_array = np.linspace(0, 2.0 * np.pi, Nz, endpoint=False)
    inter_k = sc._build_interaction_k(kx_array, ky_array, kz_array,
                                      interactions, norb)

    # --- FLEX inputs (full frequency) ---
    chis_w, chic_w, green_w, chi_convention = load_flex_chi_dynamic(
        input_dict, norb, Nx, Ny, Nz)
    if green_w is None:
        raise ValueError(
            "dynamic Eliashberg requires the dressed green.npz from the FLEX "
            "run (the pair bubble G2 is built from it); none was found. Check "
            "[file.input] path_to_flex_output / [eliashberg] flex_green.")
    nmat = chis_w.shape[-1]

    # --- Optional IR compression of the frequency axis (design Sec. 3.2):
    # everything downstream (vertex assembly, pair bubble, kernel, parity
    # machinery) is per-frequency or reversal-based, so it operates on the
    # sparse symmetric node axis unchanged. The full uniform tensors of the
    # VERTEX and G2 are never built on the IR path.
    axF = axB = None
    if use_ir:
        axF, axB = _ir_axes_for_run(eli_param, beta, hr, inter_k, norb,
                                    mu=mode_param.get("mu"),
                                    filling=mode_param.get("filling"))
        keep_static = bool(eli_param.get("ir_keep_static_chi", False))
        chis_w = _ir_compress(chis_w, axB, nmat, "chiq_s",
                              drop_constant=True, keep_constant=keep_static)
        chic_w = _ir_compress(chic_w, axB, nmat, "chiq_c",
                              drop_constant=True, keep_constant=keep_static)
        green_w = _ir_compress(green_w, axF, nmat, "green")
    nfreq_axis = axF.n_freq if use_ir else nmat

    # --- Vertex and pair bubble on the frequency axis ---
    logger.info("Computing dynamic FLEX pairing vertex (pairing_type=%s, "
                "convention=%s)...", pairing_type, chi_convention)
    Vs_q_w = compute_vertices_flex_dynamic(
        chis_w, chic_w, inter_k, norb, Nx, Ny, Nz,
        pairing_type=pairing_type, convention=chi_convention)
    logger.info("Computing frequency-resolved pair bubble G2...")
    G2_w = calc_g2_dynamic(green_w, beta)

    # --- Seed: static init_gap form factor, broadcast flat across omega ---
    gap_shape = (norb, norb, Nx, Ny, Nz, nfreq_axis)
    sigma_static = sc._initialize_gap(init_gap_mode, norb,
                                      kx_array, ky_array, kz_array)
    phi0 = np.broadcast_to(sigma_static[..., np.newaxis], gap_shape).copy()
    phi0 = phi0.astype(complex)
    n0 = np.linalg.norm(phi0)
    if n0 > 0:
        phi0 /= n0

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
                                       convention=chi_convention)
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

    def make_operator():
        op = LinearOperator((vec_size, vec_size), matvec=_matvec,
                            dtype=complex)
        return op, vec_size

    # Map [eliashberg] controls to the _solve_leading solver_mode string,
    # exactly as calc_eliashberg does for the static path.
    eigenvalue_match = None
    if solver_mode == "iteration":
        # Mirror the static _solve_iteration: project every iterate onto the
        # channel's combined-parity sector so numerical noise cannot let the
        # power iteration drift into the opposite-parity (e.g. triplet) mode.
        # The projection is legitimate only when the kernel commutes with P;
        # probe the cross-sector leakage first and fall back to the
        # un-projected iteration (with a warning) if it does not.
        A_probe, _ = make_operator()
        leak = _parity_leakage(A_probe, gap_shape, pairing_type)
        if leak <= 1.0e-8:
            def project_fn(flat):
                return _project_parity_dynamic(
                    flat.reshape(gap_shape), pairing_type).ravel()
            # Project + normalize the seed, raising if it has no in-sector
            # component (matches sc._solve_iteration's guard).
            phi0 = _project_seed_dynamic(phi0, pairing_type)
        else:
            logger.warning(
                "Dynamic Eliashberg kernel does not commute with parity "
                "(cross-sector leakage %.2e); parity projection for the '%s' "
                "channel is disabled and the un-projected iteration is used.",
                leak, pairing_type)
            project_fn = None
        eigenvalue, sigma_flat, info = sc._solve_leading(
            make_operator, vec_size, "iteration",
            max_iter=max_iter, convergence_tol=tol, alpha=alpha,
            init_vec=phi0.ravel(), project_fn=project_fn)
        eigenvalues_all = None
    else:
        # "eigenvalue" / "both": use the ARPACK/shift-invert eigen family.
        # Note: "both" degrades to eigenvalue-only here (the static path also
        # runs a power-iteration leg); the ARPACK leading pair is returned.
        if solver_mode == "both":
            logger.warning(
                "Dynamic solver_mode='both' runs the eigenvalue leg only; "
                "the power-iteration cross-check is skipped.")
        eigenvalue, sigma_flat, info = sc._solve_leading(
            make_operator, vec_size, eigenvalue_method,
            num_eigenvalues=num_eigenvalues)
        eigenvalues_all = info.get("eigenvalues")
        vecs_all = info.get("eigenvectors")
        # Promote the eigenpair with the channel's combined (k, orbital,
        # frequency) parity so the reported leading lambda is the physical
        # singlet/triplet solution -- not ARPACK's raw largest-|lambda|, which
        # on real FLEX data can be an opposite-parity (wrong-channel) mode.
        if eigenvalues_all is not None and vecs_all is not None:
            eigenvalues_all, vecs_all, eigenvalue_match = \
                _reorder_eigenpairs_by_parity_dynamic(
                    eigenvalues_all, vecs_all, gap_shape, pairing_type)
            eigenvalue = eigenvalues_all[0]
            sigma_flat = vecs_all[:, 0]

    lam = float(np.real(eigenvalue))
    logger.info("Dynamic Eliashberg leading eigenvalue lambda = %.6f", lam)

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
    output_dir = input_dict["file"]["output"]["path_to_output"]
    os.makedirs(output_dir, exist_ok=True)
    eigenvalue_file = eli_param.get("output_eigenvalue", "eigenvalue.dat")
    with open(os.path.join(output_dir, eigenvalue_file), "w") as fw:
        fw.write("# Dynamic Eliashberg leading eigenvalue\n")
        fw.write("{:.8e}\n".format(lam))
        if eigenvalues_all is not None:
            if eigenvalue_match is not None:
                fw.write("# index  Re(eigenvalue)  Im(eigenvalue)  "
                         "|eigenvalue|  match(1=channel-parity)\n")
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

    gap_file = eli_param.get("output_gap", "gap.dat")
    # Provenance metadata is added ONLY on the opt-in IR path: the default
    # uniform output keeps its exact historical key set.
    if use_ir:
        extra_meta = {"matsubara_basis": "ir", "ir_tol": axF.eps,
                      "ir_wmax": axF.wmax, "ir_L": axF.L}
    else:
        extra_meta = None
    write_dynamic_outputs(output_dir, gap_w, lam, T, pairing_type,
                          kx_array, ky_array, kz_array, beta,
                          gap_file=gap_file, extra_meta=extra_meta)

    return lam
