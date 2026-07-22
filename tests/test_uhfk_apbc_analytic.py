"""Analytic-spectrum tests for free-fermion APBC eigenvalues.

Same __new__ stub pattern as tests/test_uhfk_apbc_consistency.py — we drive
_make_ham_trans() directly and diagonalize each k-block.
"""
import numpy as np
import pytest

from hwave.solver.uhfk import UHFk


def make_trans_stub(cellshape, transfer_dict,
                    boundary_theta=(0.0, 0.0, 0.0), norb=1):
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
    return {((1, 0, 0), (0, 0)): -t, ((-1, 0, 0), (0, 0)): -t}


def nn_2d(t=1.0):
    return {
        ((1, 0, 0), (0, 0)): -t, ((-1, 0, 0), (0, 0)): -t,
        ((0, 1, 0), (0, 0)): -t, ((0, -1, 0), (0, 0)): -t,
    }


def nn_3d(t=1.0):
    return {
        ((1, 0, 0), (0, 0)): -t, ((-1, 0, 0), (0, 0)): -t,
        ((0, 1, 0), (0, 0)): -t, ((0, -1, 0), (0, 0)): -t,
        ((0, 0, 1), (0, 0)): -t, ((0, 0, -1), (0, 0)): -t,
    }


def spectrum_from_ham_trans(ham_trans):
    eigs = []
    for hk in ham_trans:
        eigs.extend(np.linalg.eigvalsh(hk).tolist())
    return np.sort(np.asarray(eigs, dtype=np.float64))


def expected_1d_apbc(L, t=1.0):
    spec = []
    for n in range(L):
        eps = -2.0 * t * np.cos((2 * n + 1) * np.pi / L)
        spec.extend([eps, eps])  # spin-degenerate
    return np.sort(np.asarray(spec))


def expected_2d(Lx, Ly, theta_x, theta_y, t=1.0):
    spec = []
    for nx in range(Lx):
        for ny in range(Ly):
            kx = (2 * np.pi * nx + theta_x) / Lx
            ky = (2 * np.pi * ny + theta_y) / Ly
            eps = -2.0 * t * (np.cos(kx) + np.cos(ky))
            spec.extend([eps, eps])
    return np.sort(np.asarray(spec))


def expected_3d(Lx, Ly, Lz, tx, ty, tz, t=1.0):
    spec = []
    for nx in range(Lx):
        for ny in range(Ly):
            for nz in range(Lz):
                kx = (2 * np.pi * nx + tx) / Lx
                ky = (2 * np.pi * ny + ty) / Ly
                kz = (2 * np.pi * nz + tz) / Lz
                eps = -2.0 * t * (np.cos(kx) + np.cos(ky) + np.cos(kz))
                spec.extend([eps, eps])
    return np.sort(np.asarray(spec))


def _apply_pre_fold_phase(s):
    """Mimic _init_interaction's APBC pre-fold phase injection (no sublattice)."""
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
    return s.ham_trans


@pytest.mark.parametrize("L", [4, 6, 8])
def test_1d_apbc_spectrum(L):
    got = spectrum_from_ham_trans(hamk([L, 1, 1], nn_1d(), (np.pi, 0.0, 0.0)))
    want = expected_1d_apbc(L)
    np.testing.assert_allclose(got, want, atol=1e-12)


def test_2d_mixed_apbc_spectrum():
    Lx, Ly = 4, 6
    got = spectrum_from_ham_trans(
        hamk([Lx, Ly, 1], nn_2d(), (np.pi, 0.0, 0.0))
    )
    want = expected_2d(Lx, Ly, np.pi, 0.0)
    np.testing.assert_allclose(got, want, atol=1e-12)


def test_3d_full_apbc_spectrum():
    Lx, Ly, Lz = 4, 4, 4
    got = spectrum_from_ham_trans(
        hamk([Lx, Ly, Lz], nn_3d(), (np.pi, np.pi, np.pi))
    )
    want = expected_3d(Lx, Ly, Lz, np.pi, np.pi, np.pi)
    np.testing.assert_allclose(got, want, atol=1e-12)


# ---- 2-orbital unit cell (norb_orig > 1) ------------------------------------
# Single-site/single-orbital chain tests above only exercise APBC against a
# trivially translation-invariant Hamiltonian. For models like honeycomb or
# kagome the unit cell already contains multiple inequivalent orbitals, and
# APBC stacks on top of that. Pin the simplest non-trivial case (a 2-orbital
# 1D chain) so multi-orbital APBC stays correct.

def two_orbital_chain_transfer(E0, E1, ta, tb, tau):
    """1D NN chain with 2 inequivalent orbitals per site.

    H = sum_x [
        E0 n_{x,0} + E1 n_{x,1}
        + tau (c^dag_{x,0} c_{x,1} + h.c.)
        - ta (c^dag_{x+1,0} c_{x,0} + h.c.)
        - tb (c^dag_{x+1,1} c_{x,1} + h.c.)
    ]
    """
    return {
        # on-site
        ((0, 0, 0), (0, 0)): complex(E0),
        ((0, 0, 0), (1, 1)): complex(E1),
        ((0, 0, 0), (0, 1)): complex(tau),
        ((0, 0, 0), (1, 0)): complex(tau),       # tau real -> hermitian
        # intra-orbital NN
        ((1, 0, 0), (0, 0)): complex(-ta), ((-1, 0, 0), (0, 0)): complex(-ta),
        ((1, 0, 0), (1, 1)): complex(-tb), ((-1, 0, 0), (1, 1)): complex(-tb),
    }


def expected_2orbital_chain_apbc(L, E0, E1, ta, tb, tau, t=1.0):
    """Analytic 2-band spectrum of the above H under APBC theta_x = pi.

    For each k in {(2n+1) pi / L : n = 0..L-1}, diagonalize the 2x2
        H(k) = [[E0 - 2 ta cos(k), tau          ],
                [tau,               E1 - 2 tb cos(k)]]
    """
    spec = []
    for n in range(L):
        k = (2 * n + 1) * np.pi / L
        diag0 = E0 - 2.0 * ta * np.cos(k)
        diag1 = E1 - 2.0 * tb * np.cos(k)
        Hk = np.array([[diag0, tau], [tau, diag1]], dtype=np.float64)
        eigs = np.linalg.eigvalsh(Hk)
        # spin-degenerate (no Zeeman): each k contributes 2 (orbital) * 2 (spin)
        spec.extend([eigs[0], eigs[0], eigs[1], eigs[1]])
    return np.sort(np.asarray(spec, dtype=np.float64))


def hamk_norb(cellshape, transfer, theta, norb):
    s = make_trans_stub(cellshape, transfer, boundary_theta=theta, norb=norb)
    s.nd = 2 * norb  # ns=2 spin sectors x norb orbitals
    _apply_pre_fold_phase(s)
    s._make_ham_trans()
    return s.ham_trans


def test_2orbital_chain_apbc_spectrum():
    """Multi-orbital unit cell (norb_orig=2) under APBC: H-wave spectrum
    must match the analytic 2-band model on the L=4 1D chain.

    This pins that the APBC gauge phase composes correctly with non-trivial
    orbital structure inside the unit cell (the situation in honeycomb,
    kagome, multi-orbital Hubbard, etc.).
    """
    L = 4
    E0, E1 = 0.5, -0.3
    ta, tb = 1.0, 0.5
    tau = 0.4

    got = spectrum_from_ham_trans(
        hamk_norb([L, 1, 1],
                  two_orbital_chain_transfer(E0, E1, ta, tb, tau),
                  (np.pi, 0.0, 0.0),
                  norb=2)
    )
    want = expected_2orbital_chain_apbc(L, E0, E1, ta, tb, tau)
    np.testing.assert_allclose(got, want, atol=1e-12)


def two_orbital_chain_transfer_peierls(E0, E1, ta, tb, tint, phi):
    """1D NN chain with 2 inequivalent orbitals AND a complex Peierls phase
    on the inter-orbital NN hop.

    Adds to the real-only model:
        + tint exp(+i phi) c^dag_{x+1, 0} c_{x, 1}  + h.c.
    Breaks time-reversal symmetry for nonzero phi, so the test exercises
    multi-orbital + complex-hopping APBC simultaneously.
    """
    trans = two_orbital_chain_transfer(E0, E1, ta, tb, tau=0.0)
    # T_{alpha=0, beta=1}(R=+1) = + tint * exp(+i phi)
    trans[((1, 0, 0), (0, 1))] = complex(tint * np.exp(+1j * phi))
    # Hermitian partner: T_{1, 0}(R=-1) = conj
    trans[((-1, 0, 0), (1, 0))] = complex(tint * np.exp(-1j * phi))
    return trans


def expected_2orbital_chain_apbc_peierls(L, E0, E1, ta, tb, tint, phi):
    """Analytic 2-band spectrum with the complex Peierls inter-orbital hop.

    H(k) = [[E0 - 2 ta cos(k),         tint * exp(+i (phi - k))],
            [conj(tint * exp(...)),   E1 - 2 tb cos(k)         ]]
    """
    spec = []
    for n in range(L):
        k = (2 * n + 1) * np.pi / L
        diag0 = E0 - 2.0 * ta * np.cos(k)
        diag1 = E1 - 2.0 * tb * np.cos(k)
        # Peierls hop in k-space: sum over R of T(R)_{01} exp(+i k R)
        # T(R=+1)_{01} = tint exp(+i phi), T(R=-1)_{01} = 0  (only the partner
        # at R=-1 with orbital flipped is set). So the off-diagonal is
        # tint exp(+i phi) * exp(+i k) = tint exp(+i(phi+k)).
        off = tint * np.exp(+1j * (phi + k))
        Hk = np.array([[diag0, off], [np.conj(off), diag1]], dtype=np.complex128)
        eigs = np.linalg.eigvalsh(Hk)
        spec.extend([eigs[0], eigs[0], eigs[1], eigs[1]])
    return np.sort(np.asarray(spec, dtype=np.float64))


def test_2orbital_chain_apbc_spectrum_complex_peierls():
    """Multi-orbital + complex (Peierls) hopping under APBC must also match
    the analytic 2-band spectrum.

    Combines the multi-orbital structure with a complex inter-orbital hop
    that breaks time-reversal symmetry, so any sign / phase mistake in the
    APBC + complex-hopping composition shows up.
    """
    L = 4
    E0, E1 = 0.5, -0.3
    ta, tb = 1.0, 0.5
    tint = 0.6
    phi = 0.4  # nonzero -> TRS-broken complex hop

    got = spectrum_from_ham_trans(
        hamk_norb([L, 1, 1],
                  two_orbital_chain_transfer_peierls(E0, E1, ta, tb, tint, phi),
                  (np.pi, 0.0, 0.0),
                  norb=2)
    )
    want = expected_2orbital_chain_apbc_peierls(L, E0, E1, ta, tb, tint, phi)
    np.testing.assert_allclose(got, want, atol=1e-12)
