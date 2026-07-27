import numpy as np
import pytest

from hwave.solver import eliashberg_dynamic as ed
from tests.oracle_bond_dynamic import (
    oracle_bond_matvec, oracle_bond_bubble, oracle_pair_bubble_X)

RNG = np.random.default_rng(20260728)


def _symmetric_green(norb, Nx, Ny, Nz, nmat, beta):
    """Random green in the v1 symmetry class: G(-k,iw)=G(k,iw),
    G(k,-iw)=conj(G(k,iw)).  Built from a real inversion-symmetric
    dispersion: G = 1/(i w_n - e(k))."""
    e = RNG.normal(size=(Nx, Ny, Nz))
    e = 0.5 * (e + np.roll(e[::-1, ::-1, ::-1], (1, 1, 1), (0, 1, 2)))
    n_t = np.arange(nmat) - nmat // 2
    iw = 1j * (2 * n_t + 1) * np.pi / beta
    G = 1.0 / (iw[None, None, None, :] - e[..., None])
    return G[None, None, ...].astype(complex)  # (1,1,Nx,Ny,Nz,nmat)


def test_oracle_matvec_equals_scalar_kernel_at_B1():
    Nx = Ny = 4; Nz = 1; nmat = 8; beta = 5.0
    green = _symmetric_green(1, Nx, Ny, Nz, nmat, beta)
    phi = (RNG.normal(size=(1, 1, Nx, Ny, Nz, nmat))
           + 1j * RNG.normal(size=(1, 1, Nx, Ny, Nz, nmat)))
    # scalar vertex on the bosonic axis -> bond F_w with B=1, ND=1
    Vs = (RNG.normal(size=(Nx, Ny, Nz, nmat))
          + 1j * RNG.normal(size=(Nx, Ny, Nz, nmat)))
    # enforce the bosonic reality class F(q,j) = conj(F(-q, -j))
    Vs = 0.5 * (Vs + np.conj(np.roll(
        Vs[::-1, ::-1, ::-1, ::-1], (1, 1, 1, 1), (0, 1, 2, 3))))
    F_w = Vs[..., None, None]                       # (Nx,Ny,Nz,nmat,1,1)
    Y_oracle = oracle_bond_matvec(phi, green, F_w, np.zeros((1, 1)),
                                  [(0, 0, 0)], beta)
    # production scalar kernel: Vs_q_w shape (1,1,1,1,Nx,Ny,Nz,nmat),
    # G2_w = calc_g2_dynamic (carries 1/beta)
    Vs_q_w = Vs[None, None, None, None]
    G2_w = ed.calc_g2_dynamic(green, beta)
    Y_prod = ed.eliashberg_kernel_dynamic(Vs_q_w, G2_w, phi, 1, beta)
    np.testing.assert_allclose(Y_oracle, Y_prod, rtol=1e-10, atol=1e-12)


def test_oracle_bubble_B1_zero_transfer_matches_static_chi0_shape():
    Nx = Ny = 4; Nz = 1; nmat = 8; beta = 5.0
    green = _symmetric_green(1, Nx, Ny, Nz, nmat, beta)
    chi = oracle_bond_bubble(green, [(0, 0, 0)], beta, [0])[0]
    assert chi.shape == (Nx, Ny, Nz, 1, 1)
    # particle-hole flip symmetry at zero transfer: chi(q) = chi(-q)*
    chi_neg = np.conj(np.roll(chi[::-1, ::-1, ::-1], (1, 1, 1), (0, 1, 2)))
    np.testing.assert_allclose(chi, chi_neg, rtol=1e-10, atol=1e-12)
