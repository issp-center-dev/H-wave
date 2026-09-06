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
