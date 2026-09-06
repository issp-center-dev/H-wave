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
