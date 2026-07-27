"""Independent direct-summation oracle for the dynamic bond kernel.

Normative reference for spec 2026-07-27-dynamic-bond-channels-design.md
sections 3.1/3.3. MUST NOT import hwave.solver.matsubara, hwave.solver
backend spatial helpers, or hwave.solver.ir_axis (oracle-independence
requirement, spec section 3.5).
"""
import numpy as np


def _neg_k(i, n):
    return (-i) % n


def _neg_n(n_idx, nmat):
    # omega_{n~} -> -omega_{n~} on the symmetric grid: array index nmat-1-n
    return nmat - 1 - n_idx


def _phase(kx, ky, kz, d, Nx, Ny, Nz):
    return np.exp(2j * np.pi * (kx * d[0] / Nx + ky * d[1] / Ny
                                + kz * d[2] / Nz))


def oracle_pair_bubble_X(green, phi, beta):
    """X_{(l3 l4)}(k', n') = sum_{l5 l6} G_{l3 l5}(k',n') phi_{l5 l6}(k',n')
    G_{l4 l6}(-k',-n').  Normative: NO 1/beta factor here (spec 3.3)."""
    norb = green.shape[0]
    Nx, Ny, Nz, nmat = green.shape[2:6]
    nd = norb * norb
    X = np.zeros((nd, Nx, Ny, Nz, nmat), dtype=complex)
    for kx in range(Nx):
        for ky in range(Ny):
            for kz in range(Nz):
                mkx, mky, mkz = _neg_k(kx, Nx), _neg_k(ky, Ny), _neg_k(kz, Nz)
                for n in range(nmat):
                    mn = _neg_n(n, nmat)
                    G = green[:, :, kx, ky, kz, n]
                    Gm = green[:, :, mkx, mky, mkz, mn]
                    P = phi[:, :, kx, ky, kz, n]
                    # X_{l3 l4} = sum_{l5 l6} G_{l3 l5} P_{l5 l6} Gm_{l4 l6}
                    Xm = G @ P @ Gm.T
                    X[:, kx, ky, kz, n] = Xm.reshape(nd)
    return X


def oracle_bond_bubble(green, delta_r, beta, j_list):
    """chi_bar at bosonic transfers j~ in j_list (spec 3.1), circular wrap.

    chi_bar_{(m,idx),(m',idx')}(q, i nu_j~)
      = -(T/N) sum_{k,n~} e^{i k.(dr_m - dr_m')}
          G_{l1 l3}(k+q, w_{n~+j~}) G_{l4 l2}(k, w_{n~})
    with the fermionic index n+j wrapped mod nmat (stored-window torus).
    idx = l1*norb + l2 (row), idx' = l3*norb + l4 (column).
    """
    norb = green.shape[0]
    Nx, Ny, Nz, nmat = green.shape[2:6]
    nd = norb * norb
    B = len(delta_r)
    ND = nd * B
    N = Nx * Ny * Nz
    T = 1.0 / beta
    out = {}
    for jt in j_list:
        chi = np.zeros((Nx, Ny, Nz, ND, ND), dtype=complex)
        for qx in range(Nx):
            for qy in range(Ny):
                for qz in range(Nz):
                    acc = np.zeros((ND, ND), dtype=complex)
                    for kx in range(Nx):
                        for ky in range(Ny):
                            for kz in range(Nz):
                                kqx, kqy, kqz = ((kx + qx) % Nx,
                                                 (ky + qy) % Ny,
                                                 (kz + qz) % Nz)
                                for n in range(nmat):
                                    nj = (n + jt) % nmat   # circular wrap
                                    Gkq = green[:, :, kqx, kqy, kqz, nj]
                                    Gk = green[:, :, kx, ky, kz, n]
                                    for m in range(B):
                                        for mp in range(B):
                                            ph = _phase(
                                                kx, ky, kz,
                                                np.subtract(delta_r[m],
                                                            delta_r[mp]),
                                                Nx, Ny, Nz)
                                            for l1 in range(norb):
                                                for l2 in range(norb):
                                                    for l3 in range(norb):
                                                        for l4 in range(norb):
                                                            I = m * nd + l1 * norb + l2
                                                            J = mp * nd + l3 * norb + l4
                                                            acc[I, J] += (
                                                                ph
                                                                * Gkq[l1, l3]
                                                                * Gk[l4, l2])
                    chi[qx, qy, qz] = -(T / N) * acc
        out[jt] = chi
    return out


def oracle_bond_matvec(phi, green, F_w, Vpp, delta_r, beta):
    """Normative matvec, spec 3.3 boxed equation.

    phi   : (norb, norb, Nx, Ny, Nz, nmat)
    F_w   : (Nx, Ny, Nz, nmat, ND, ND)  bosonic axis, zero transfer at
            array index nmat//2 (spec 3.1 map)
    Vpp   : (ND, ND) frequency-independent Cooper vertex
    returns Y, same shape as phi.
    """
    norb = green.shape[0]
    Nx, Ny, Nz, nmat = green.shape[2:6]
    nd = norb * norb
    B = len(delta_r)
    N = Nx * Ny * Nz
    T = 1.0 / beta
    X = oracle_pair_bubble_X(green, phi, beta)   # (nd, Nx,Ny,Nz, nmat)

    # --- instantaneous coefficients (external-frequency independent) -----
    # A_{m; cd} = (T/N) sum_{k', n'} e^{-i k'.dr_m} X_{cd}(k', n')
    A = np.zeros((B, nd), dtype=complex)
    for m in range(B):
        for kx in range(Nx):
            for ky in range(Ny):
                for kz in range(Nz):
                    ph = np.conj(_phase(kx, ky, kz, delta_r[m], Nx, Ny, Nz))
                    A[m] += ph * X[:, kx, ky, kz, :].sum(axis=-1)
    A *= (T / N)

    Y = np.zeros_like(phi)
    for kx in range(Nx):
        for ky in range(Ny):
            for kz in range(Nz):
                for n in range(nmat):
                    y = np.zeros(nd, dtype=complex)
                    # fluctuation: -(T/N) sum_{k',n'} F(q, j) with
                    # q = k - k' (mod), j = (n - n') + nmat//2 (mod)
                    for kpx in range(Nx):
                        for kpy in range(Ny):
                            for kpz in range(Nz):
                                qx, qy, qz = ((kx - kpx) % Nx,
                                              (ky - kpy) % Ny,
                                              (kz - kpz) % Nz)
                                for npr in range(nmat):
                                    j = (n - npr + nmat // 2) % nmat
                                    Fq = F_w[qx, qy, qz, j]      # (ND, ND)
                                    for m in range(B):
                                        phm = _phase(kx, ky, kz,
                                                     delta_r[m], Nx, Ny, Nz)
                                        for mp in range(B):
                                            phmp = _phase(kpx, kpy, kpz,
                                                          delta_r[mp],
                                                          Nx, Ny, Nz)
                                            blk = Fq[m * nd:(m + 1) * nd,
                                                     mp * nd:(mp + 1) * nd]
                                            y += (phm * phmp
                                                  * (blk @ X[:, kpx, kpy,
                                                             kpz, npr]))
                    y *= -(T / N)
                    # instantaneous: -(1/2) sum_{m,m'} Vpp e^{+i k.dr_m'} A_m
                    for m in range(B):
                        for mp in range(B):
                            phmp = _phase(kx, ky, kz, delta_r[mp],
                                          Nx, Ny, Nz)
                            blk = Vpp[m * nd:(m + 1) * nd,
                                      mp * nd:(mp + 1) * nd]
                            y += -0.5 * phmp * (blk @ A[m])
                    Y[:, :, kx, ky, kz, n] = y.reshape(norb, norb)
    return Y
