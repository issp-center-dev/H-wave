"""FLEX-side Hartree-Fock adapter (GitHub issue #181, Tier 3 Phase B; spec
docs/superpowers/specs/2026-09-06-flex-bond-sigma-phase-b-181-design.md,
section 2.2).

The spin-free FLEX density ``rho_ab(r)`` (one spin block, per spin) is
embedded paramagnetically into UHFk's spin-major layout and passed through
the shared kernel (:mod:`hwave.solver.hartree_fock`); the up block of the
result is the static Hartree-Fock self-energy ``Sigma_HF(k)``. Every
accepted interaction entry reaches the kernel (the table builder's
``discarded`` report is turned into a refusal here); a contraction may be
symmetry-zero (PairLift, D12).
"""
import numpy as np

from . import hartree_fock as _hf


def build_flex_hf_tables(param_ham, norb, shape):
    """The shared interaction tables for the FLEX HF map. Refuses
    (``ValueError``) when the per-type semantics would drop a declared
    entry, naming the first one, so that D4's "every accepted term" holds."""
    tabs = _hf.build_interaction_tables(param_ham, norb, shape)
    if tabs.discarded:
        t, irvec, orbvec, v = tabs.discarded[0]
        raise ValueError(
            "flex_hartree_fock: the '{}' entry at irvec={}, orbvec={} (value "
            "{!r}) has no Hartree-Fock representation in this solver ({} "
            "entries must be on-site and orbital-diagonal) and would be "
            "silently dropped; remove or correct the declaration ({} "
            "such entries).".format(t, irvec, orbvec, v, t, len(tabs.discarded)))
    return tabs


def hf_map(rho_r, tables, shape, norb, *, block_tol=1e-12, herm_tol=1e-10):
    """``Sigma_HF(k)`` `(nvol, norb, norb)` from the per-spin real-space
    density ``rho_r`` `(nvol, norb, norb)` (spec 2.2): paramagnetic
    embedding ``rho_so[:, s, a, t, b] = delta_st rho_ab``, the kernel with
    the Fock term, then the spin-block checks (up == down, spin-off-
    diagonal zero, both to ``block_tol`` relative -- a violation is a bug
    and is refused) and the k-space Hermiticity check (``herm_tol``).
    Returns an OWNING copy."""
    rho_r = np.asarray(rho_r)
    nvol = rho_r.shape[0]
    rho_so = np.zeros((nvol, 2, norb, 2, norb), dtype=np.complex128)
    rho_so[:, 0, :, 0, :] = rho_r
    rho_so[:, 1, :, 1, :] = rho_r
    out = np.zeros((nvol, 2 * norb, 2 * norb), dtype=np.complex128)
    _hf.accumulate_hf(out, rho_so, tables.inter_table, tables.spin_table,
                      tuple(shape), include_fock=True)
    if not np.all(np.isfinite(out)):
        raise _hf.NonFiniteError("hf_map: non-finite Hartree-Fock term")
    o = out.reshape(nvol, 2, norb, 2, norb)
    scale = max(1.0, float(np.max(np.abs(o))))
    d_blocks = float(np.max(np.abs(o[:, 0, :, 0, :] - o[:, 1, :, 1, :])))
    d_off = float(max(np.max(np.abs(o[:, 0, :, 1, :])), np.max(np.abs(o[:, 1, :, 0, :]))))
    if d_blocks > block_tol * scale or d_off > block_tol * scale:
        raise ValueError(
            "hf_map: the paramagnetic embedding produced unequal spin blocks "
            "({:.3e}) or a spin-off-diagonal term ({:.3e}); this is a "
            "kernel/embedding defect, not a physical signal".format(d_blocks, d_off))
    sigma = np.array(o[:, 0, :, 0, :], copy=True)
    ok, err = _hf.is_hermitian_batch(sigma, herm_tol)
    if not ok:
        raise ValueError(
            "hf_map: Sigma_HF(k) is not Hermitian (relative deviation {:.3e}); "
            "the density matrix violates rho_ab(r) = conj(rho_ba(-r))".format(err))
    return sigma


# =============================================================================
# Self-energy seeds (spec section 2.3)
# =============================================================================

from dataclasses import dataclass as _dataclass


@_dataclass(frozen=True)
class SigmaSeedEnvelope:
    """An UNVALIDATED self-energy seed as loaded from a ``sigma.npz``:
    every required member materialised and read-only, the convention
    marker (``"split"``, ``"total"`` or ``None`` for a legacy file), the
    IR metadata (``None`` for uniform files) and the file name for
    diagnostics. Semantic validation happens in ``solve`` at refusal
    precedence step 5 (:func:`validate_split_seed`)."""
    sigma: np.ndarray
    sigma_static: object
    sigma_fluct: object
    marker: object
    ir_meta: object
    file_name: str


def _readonly(a):
    a = np.array(a, dtype=np.complex128, copy=True)
    a.flags.writeable = False
    return a


def make_seed_envelope(data, file_name, sigma, ir_meta):
    """Build the envelope from an opened npz ``data`` (its ``sigma`` and
    ``ir_meta`` already extracted by the caller)."""
    marker = None
    if "sigma_convention" in data.files:
        marker = str(np.asarray(data["sigma_convention"]).ravel()[0])
    st = _readonly(data["sigma_static"]) if "sigma_static" in data.files else None
    fl = _readonly(data["sigma_fluct"]) if "sigma_fluct" in data.files else None
    return SigmaSeedEnvelope(sigma=_readonly(sigma), sigma_static=st, sigma_fluct=fl,
                             marker=marker, ir_meta=ir_meta, file_name=str(file_name))


def split_tail_diagnostic(fluct, nmat):
    """D11 tail check of a split seed: on the outer 25 % of the positive
    Matsubara window, the pair-even part ``e(n) = [s(n) + s(-n)]/2`` (zero for a
    pure ``1/iw`` tail, constant for a frequency-independent remainder) is
    compared with the pair-odd part ``o(n)``. Returns ``(r, e_max)`` with
    ``r = max_n |e(n)| / max(max_n |o(n)|, 1e-300)``, or ``None`` when
    ``nmat < 32`` (no window)."""
    nmat = int(nmat)
    if nmat < 32:
        return None
    f = np.asarray(fluct)
    n_outer = max(nmat // 8, 1)
    hi = np.arange(nmat - n_outer, nmat)
    lo = nmat - 1 - hi
    e = 0.5 * (f[:, hi] + f[:, lo])
    o = 0.5 * (f[:, hi] - f[:, lo])
    e_max = float(np.max(np.abs(e))) if e.size else 0.0
    o_max = float(np.max(np.abs(o))) if o.size else 0.0
    return e_max / max(o_max, 1e-300), e_max


def validate_split_seed(env, expected_shape, *, sum_tol=1e-12, herm_tol=1e-10):
    """Semantic validation of a ``"split"`` seed (spec 2.3): marker, exact
    shapes, finiteness, singleton frequency axis of ``sigma_static``, the
    sum identity, Hermiticity of ``sigma_static``. Returns the components
    EXACTLY as stored (copies)."""
    f = env.file_name
    if env.marker != "split":
        raise ValueError(
            "sigma_init '{}': sigma_convention must be \"split\" for a run with "
            "flex_hartree_fock=true (got {!r}); convert a total archive with "
            "hwave_sigma_split".format(f, env.marker))
    if env.sigma_static is None or env.sigma_fluct is None:
        raise ValueError("sigma_init '{}': a split archive needs sigma_static and sigma_fluct".format(f))
    expected = tuple(int(x) for x in expected_shape)
    nb, nmat, nvol, n1, n2 = expected
    if tuple(env.sigma.shape) != expected:
        raise ValueError("sigma_init '{}': sigma shape {} != expected {}".format(f, env.sigma.shape, expected))
    if tuple(env.sigma_fluct.shape) != expected:
        raise ValueError("sigma_init '{}': sigma_fluct shape {} != expected {}".format(f, env.sigma_fluct.shape, expected))
    if tuple(env.sigma_static.shape) != (nb, 1, nvol, n1, n2):
        raise ValueError(
            "sigma_init '{}': sigma_static shape {} != {} (it must be frequency independent, "
            "singleton frequency axis)".format(f, env.sigma_static.shape, (nb, 1, nvol, n1, n2)))
    for name, a in (("sigma", env.sigma), ("sigma_static", env.sigma_static), ("sigma_fluct", env.sigma_fluct)):
        if not np.all(np.isfinite(a)):
            raise ValueError("sigma_init '{}': {} has non-finite entries".format(f, name))
    scale = max(1.0, float(np.max(np.abs(env.sigma))))
    dev = float(np.max(np.abs(env.sigma - (env.sigma_static + env.sigma_fluct))))
    if dev > sum_tol * scale:
        raise ValueError(
            "sigma_init '{}': sigma != sigma_static + sigma_fluct (max deviation {:.3e})".format(f, dev))
    ok, err = _hf.is_hermitian_batch(env.sigma_static[:, 0], herm_tol)
    if not ok:
        raise ValueError("sigma_init '{}': sigma_static is not Hermitian (relative deviation {:.3e})".format(f, err))
    return (np.array(env.sigma_static, copy=True), np.array(env.sigma_fluct, copy=True))


# =============================================================================
# hwave_sigma_split converter (spec 2.3)
# =============================================================================

_SPLIT_FIELDS = ("sigma_convention", "sigma_static", "sigma_fluct")
_STATIC_META = ("cell_shape", "momentum_convention", "wavevector_unit", "wavevector_index")


def _spin_reduce(h, norb, label):
    """(nvol, 2 norb, 2 norb) spin-major -> (nvol, norb, norb); refuses a
    spin-mixing or unequal-block matrix (Phase B is spin-free)."""
    nvol = h.shape[0]
    h5 = h.reshape(nvol, 2, norb, 2, norb)
    scale = max(1.0, float(np.max(np.abs(h5))))
    off = max(float(np.max(np.abs(h5[:, 0, :, 1, :]))), float(np.max(np.abs(h5[:, 1, :, 0, :]))))
    diff = float(np.max(np.abs(h5[:, 0, :, 0, :] - h5[:, 1, :, 1, :])))
    if off > 1e-12 * scale or diff > 1e-12 * scale:
        raise ValueError("{}: the spin-major matrix is spin-mixing ({:.3e}) or has unequal spin "
                         "blocks ({:.3e}); a Phase B seed must be spin-free".format(label, off, diff))
    return np.ascontiguousarray(h5[:, 0, :, 0, :])


def _bare_transfer_k(transfer_file, cell_shape, norb):
    """H0_bare(k) from a Transfer input (Wannier90-style text or the .npz
    form) exactly as the k-space reader assembles it (e^{+ikR}, C-order
    momentum layout), (nvol, norb, norb)."""
    nx, ny, nz = (int(x) for x in cell_shape)
    nvol = nx * ny * nz
    if str(transfer_file).endswith(".npz"):
        data = np.load(transfer_file)
        if "Transfer" not in data:
            raise ValueError("--bare-transfer {}: the archive carries no 'Transfer' member"
                             .format(transfer_file))
        tab_r = np.asarray(data["Transfer"], dtype=np.complex128)
        if tab_r.shape not in ((nvol, norb, norb), (nx, ny, nz, norb, norb)):
            raise ValueError("--bare-transfer {}: Transfer shape {} does not match cell {} and "
                             "norb {}".format(transfer_file, tab_r.shape, list(cell_shape), norb))
        tab_r = tab_r.reshape(nx, ny, nz, norb, norb)
    else:
        from hwave.qlmsio.wan90 import read_w90
        table = read_w90(str(transfer_file))
        tab_r = np.zeros((nx, ny, nz, norb, norb), dtype=np.complex128)
        for (irvec, orbvec), v in table.items():
            if orbvec[0] < norb and orbvec[1] < norb:
                tab_r[(*irvec, *orbvec)] += v
            else:
                raise ValueError("--bare-transfer {}: orbital index {} exceeds norb {} of the "
                                 "total archive (a spin-dependent Transfer file is not a bare "
                                 "one-body transfer)".format(transfer_file, orbvec, norb))
    return (np.fft.ifftn(tab_r, axes=(0, 1, 2)) * nvol).reshape(nvol, norb, norb)


def _static_from_uhfk_trans_mod(trans_mod_file, transfer_file, cell_shape, norb):
    """Delta H = H_mod(k) - H0_bare(k) from a native UHFk ``trans_mod``
    archive (real-space, spin-major ``(nvol, 2 norb, 2 norb)``), with the
    solver's own Fourier convention (``RPA._read_trans_mod``)."""
    data = np.load(trans_mod_file)
    if "trans_mod" not in data:
        raise ValueError("--uhfk-trans-mod {}: no 'trans_mod' member".format(trans_mod_file))
    tab_r = np.asarray(data["trans_mod"], dtype=np.complex128)
    nx, ny, nz = (int(x) for x in cell_shape)
    nvol = nx * ny * nz
    nd = 2 * norb
    if tab_r.shape != (nvol, nd, nd):
        raise ValueError("--uhfk-trans-mod {}: trans_mod shape {} does not match (nvol, 2 norb, "
                         "2 norb) = {} of the total archive (a sublattice run must be converted "
                         "on its deflated cell)".format(trans_mod_file, tab_r.shape, (nvol, nd, nd)))
    h_mod = (np.fft.ifftn(tab_r.reshape(nx, ny, nz, nd, nd), axes=(0, 1, 2)) * nvol
             ).reshape(nvol, nd, nd)
    h_mod = _spin_reduce(h_mod, norb, "--uhfk-trans-mod")
    h0 = _bare_transfer_k(transfer_file, cell_shape, norb)
    return (h_mod - h0)[None, None]


def _static_from_file(static_file, total, norb, nvol, spin_major, log):
    data = np.load(static_file)
    if "sigma_static" in data:
        st = np.asarray(data["sigma_static"])
    elif "sigma" in data:
        log.info("--static {}: no 'sigma_static' member; using 'sigma' as the static correction"
                 .format(static_file))
        st = np.asarray(data["sigma"])
    else:
        raise ValueError("--static {}: neither 'sigma_static' nor 'sigma' is present".format(static_file))
    for key in _STATIC_META:
        if key not in data or key not in total:
            raise ValueError("--static {}: the mandatory metadata member '{}' must be present in "
                             "both the static and the total archive".format(static_file, key))
        a, b = data[key], total[key]
        same = (str(a) == str(b)) if key == "momentum_convention" else (
            np.asarray(a).shape == np.asarray(b).shape and np.array_equal(np.asarray(a), np.asarray(b)))
        if not same:
            raise ValueError("--static {}: metadata '{}' differs from the total archive ({!r} vs "
                             "{!r})".format(static_file, key, a, b))
    if spin_major:
        if st.shape != (nvol, 2 * norb, 2 * norb):
            raise ValueError("--static --uhfk-spin-major: expected shape {} , got {}".format(
                (nvol, 2 * norb, 2 * norb), st.shape))
        st = _spin_reduce(np.asarray(st, dtype=np.complex128), norb, "--static")[None, None]
    elif st.ndim == 4:
        if st.shape != (1, nvol, norb, norb):
            raise ValueError("--static: rank-4 sigma_static must be (1, nvol, norb, norb) = {}, got {}"
                             .format((1, nvol, norb, norb), st.shape))
        st = st[:, None]
    elif st.shape != (1, 1, nvol, norb, norb):
        raise ValueError("--static: sigma_static must be (1, 1, nvol, norb, norb) = {} (or rank 4 "
                         "without the frequency axis), got {}".format((1, 1, nvol, norb, norb), st.shape))
    return np.array(st, dtype=np.complex128, copy=True)


def sigma_split_convert(total_path, out_path, *, static_path=None, zero_static=False,
                        uhfk_trans_mod=None, bare_transfer=None, uhfk_spin_major=False,
                        force=False, logger=None):
    """Write a ``"split"`` self-energy archive from a total-form one (spec
    2.3). Exactly one static source: ``static_path`` (the static.npz
    contract), ``zero_static`` or ``uhfk_trans_mod`` + ``bare_transfer``.
    Returns a summary dict; raises ``ValueError`` / ``FileExistsError`` on
    any refusal (nothing is written then)."""
    import logging
    import os
    log = logger or logging.getLogger(__name__)
    n_src = int(static_path is not None) + int(zero_static) + int(uhfk_trans_mod is not None)
    if n_src != 1:
        raise ValueError("exactly one static source is required (--static, --zero-static or "
                         "--uhfk-trans-mod)")
    if (uhfk_trans_mod is not None) != (bare_transfer is not None):
        raise ValueError("--bare-transfer is required with --uhfk-trans-mod and accepted only with it")
    if uhfk_spin_major and static_path is None:
        raise ValueError("--uhfk-spin-major applies to --static only")
    out_path = str(out_path)
    if not out_path.endswith(".npz"):
        out_path += ".npz"          # numpy.savez appends it; check and report the real name
    if os.path.exists(out_path) and not force:
        raise FileExistsError("output '{}' exists; pass --force to overwrite".format(out_path))
    total = np.load(total_path)
    marker = str(total["sigma_convention"]) if "sigma_convention" in total else None
    if marker not in (None, "total"):
        raise ValueError("'{}' carries sigma_convention={!r}; only a total-form archive (no "
                         "marker or \"total\") can be converted".format(total_path, marker))
    if "sigma" not in total:
        raise ValueError("'{}' has no 'sigma' member".format(total_path))
    sigma = np.asarray(total["sigma"], dtype=np.complex128)
    if sigma.ndim != 5 or sigma.shape[0] != 1 or sigma.shape[-1] != sigma.shape[-2]:
        raise ValueError("'{}': sigma must be rank 5 (1, nmat, nvol, norb, norb), got {}".format(
            total_path, sigma.shape))
    if not np.all(np.isfinite(sigma)):
        raise ValueError("'{}': sigma is not finite".format(total_path))
    _, nmat, nvol, norb, _ = sigma.shape
    if "cell_shape" in total:
        cs = tuple(int(x) for x in total["cell_shape"])
        if int(np.prod(cs)) != nvol:
            raise ValueError("'{}': cell_shape {} does not match nvol {}".format(total_path, cs, nvol))
    if zero_static:
        static = np.zeros((1, 1, nvol, norb, norb), dtype=np.complex128)
    elif static_path is not None:
        static = _static_from_file(static_path, total, norb, nvol, uhfk_spin_major, log)
    else:
        if "cell_shape" not in total:
            raise ValueError("--uhfk-trans-mod needs the cell_shape member of the total archive")
        static = _static_from_uhfk_trans_mod(uhfk_trans_mod, bare_transfer, total["cell_shape"], norb)
    if not np.all(np.isfinite(static)):
        raise ValueError("the static correction is not finite")
    from hwave.solver.hartree_fock import is_hermitian_batch
    ok, err = is_hermitian_batch(static[:, 0], 1e-10)
    if not ok:
        raise ValueError("the static correction is not Hermitian (relative deviation {:.3e})".format(err))
    fluct = sigma - static
    if not np.all(np.isfinite(fluct)):
        raise ValueError("the fluctuation part sigma - sigma_static is not finite")
    # the solver's own seed validation, IN MEMORY, before anything is written
    env = SigmaSeedEnvelope(sigma=_readonly(sigma), sigma_static=_readonly(static),
                            sigma_fluct=_readonly(fluct), marker="split", ir_meta=None,
                            file_name=out_path)
    validate_split_seed(env, sigma.shape)
    members = {k: total[k] for k in total.files if k not in _SPLIT_FIELDS and k != "sigma"}
    members.update(sigma=sigma, sigma_convention=np.str_("split"), sigma_static=static,
                   sigma_fluct=fluct)
    np.savez(out_path, **members)
    # read-back check of the written archive
    data = np.load(out_path)
    validate_split_seed(make_seed_envelope(data, out_path, data["sigma"], None), sigma.shape)
    return dict(out=out_path, nmat=nmat, nvol=nvol, norb=norb,
                static_max=float(np.max(np.abs(static))), fluct_max=float(np.max(np.abs(fluct))))
