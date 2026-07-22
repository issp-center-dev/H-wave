"""APBC vs PBC internal consistency tests, driven via __new__ stubs.

We bypass the full UHFk constructor (which needs a complete param_ham and
loader-side state) and exercise _make_ham_trans() directly with just the
attributes it reads. This mirrors the stub pattern in
tests/test_uhfk_zero_t_occupation.py.
"""
import numpy as np
import pytest

from hwave.solver.uhfk import UHFk


def make_trans_stub(cellshape, transfer_dict,
                    boundary_theta=(0.0, 0.0, 0.0), norb=1):
    """Build a minimal UHFk for calling _make_ham_trans().

    Assumes SubShape = [1, 1, 1] (so shape == cellshape) and the
    non-spin-orbital branch.
    """
    s = object.__new__(UHFk)
    s.shape = tuple(cellshape)
    s.cellshape = tuple(cellshape)
    s.nvol = cellshape[0] * cellshape[1] * cellshape[2]
    s.norb = norb
    s.ns = 2
    s.nd = 2 * norb
    s.enable_spin_orbital = False
    s.param_ham = {"Transfer": dict(transfer_dict)}
    s.boundary_theta = tuple(boundary_theta)
    s.boundary_periodic = all(t == 0.0 for t in boundary_theta)
    return s


def nn_1d(t=1.0):
    """1D nearest-neighbor hopping (size-independent: dict keys are R vectors)."""
    return {((1, 0, 0), (0, 0)): -t, ((-1, 0, 0), (0, 0)): -t}


def _apply_pre_fold_phase(s):
    """Mimic _init_interaction's APBC pre-fold phase injection.

    Production wires APBC phase into Transfer in _init_interaction, before any
    sublattice fold, using the original signed irvec. Stub-based tests that
    bypass __init__ and call _make_ham_trans() directly must apply the phase
    themselves; this helper does that for the no-sublattice case.
    """
    if s.boundary_periodic:
        return
    from hwave.solver._apbc_phase import transfer_phase
    theta_arr = np.array(s.boundary_theta, dtype=np.float64)
    L_arr = np.array(s.cellshape, dtype=np.float64)
    s.param_ham["Transfer"] = {
        k: v * transfer_phase(np.asarray(k[0], dtype=np.float64), theta_arr, L_arr)
        for k, v in s.param_ham["Transfer"].items()
    }


def hamk(cellshape, transfer, theta):
    s = make_trans_stub(cellshape, transfer, boundary_theta=theta)
    _apply_pre_fold_phase(s)
    s._make_ham_trans()
    return s.ham_trans.copy()


def test_pbc_path_is_hermitian():
    """Sanity: with theta = 0 the phase code path is skipped and H(k) remains Hermitian."""
    h = hamk([4, 1, 1], nn_1d(), theta=(0.0, 0.0, 0.0))
    for hk in h:
        np.testing.assert_allclose(hk, hk.conj().T, atol=1e-14)


def test_apbc_x_differs_from_pbc():
    """APBC must produce a different H(k) than PBC for non-trivial hopping."""
    h_pbc = hamk([4, 1, 1], nn_1d(), theta=(0.0, 0.0, 0.0))
    h_apbc = hamk([4, 1, 1], nn_1d(), theta=(np.pi, 0.0, 0.0))
    assert not np.allclose(h_pbc, h_apbc)


def test_L_apbc_equals_2L_pbc_odd_k_subset():
    """1-body APBC spectrum of L sites = odd-k subset of 2L PBC spectrum."""
    L = 4
    # APBC on L
    h_apbc = hamk([L, 1, 1], nn_1d(), (np.pi, 0.0, 0.0))
    eigs_apbc = []
    for hk in h_apbc:
        eigs_apbc.extend(np.linalg.eigvalsh(hk).tolist())
    spec_apbc = np.unique(np.round(np.sort(eigs_apbc), 10))

    # PBC on 2L: collect eigenvalues at FFT indices n=1,3,...,2L-1 (odd k).
    h_pbc = hamk([2 * L, 1, 1], nn_1d(), (0.0, 0.0, 0.0))
    odd_eigs = []
    for n_super, hk in enumerate(h_pbc):
        if n_super % 2 == 1:
            odd_eigs.extend(np.linalg.eigvalsh(hk).tolist())
    spec_odd = np.unique(np.round(np.sort(odd_eigs), 10))

    np.testing.assert_allclose(spec_apbc, spec_odd, atol=1e-10)


def test_hk_is_hermitian_under_apbc():
    L = 6
    h = hamk([L, 1, 1], nn_1d(), (np.pi, 0.0, 0.0))
    for n, hk in enumerate(h):
        np.testing.assert_allclose(
            hk, hk.conj().T, atol=1e-12,
            err_msg=f"H(k_{n}) not Hermitian under APBC",
        )
