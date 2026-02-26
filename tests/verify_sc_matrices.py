#!/usr/bin/env python3
"""Verify S/C matrices for all interaction types including Ising and PairHop.

Verifies the S/C values in sc.py by comparing the Kuroki formula
with the full 4-index RPA for a paramagnetic system.

Convention (Kuroki et al., arXiv:0902.3691):
  chi_s = [1 - chi0 S]^{-1} chi0
  chi_c = [1 + chi0 C]^{-1} chi0

S/C contributions by interaction type (in Kuroki convention):
  CoulombIntra U: diag S=U, C=U
  CoulombInter V: cross S=V, C=-V+J; dens S=J, C=2V-J (Kanamori)
  Hund J: contributes to dens/cross (part of Kanamori)
  Exchange J': exch S=J', C=J'
  Ising I: cross S=-I, C=-I; dens S=-2I, C=0
  PairLift P: S=0, C=0 (no contribution in ph channel)
  PairHop P: exch S=P, C=P
"""

import numpy as np
import sys
sys.path.insert(0, 'src')


def build_gamma_so(interactions, norb):
    """Build 4-index RPA vertex in spin-orbital space (rpa.py convention)."""
    nd = 2 * norb
    spin_table = {
        'CoulombIntra': {(0,0,1,1): 1, (1,1,0,0): 1},
        'CoulombInter': {(0,0,0,0): 1, (1,1,1,1): 1, (0,0,1,1): 1, (1,1,0,0): 1},
        'Hund':         {(0,0,0,0): -1, (1,1,1,1): -1},
        'Ising':        {(0,0,0,0): 1, (1,1,1,1): 1, (0,0,1,1): -1, (1,1,0,0): -1},
        'PairLift':     {(0,1,0,1): 1, (1,0,1,0): 1},
        'Exchange':     {(0,1,1,0): -1, (1,0,0,1): -1},
        'PairHop':      {(0,0,1,1): 1, (1,1,0,0): 1},
    }
    Gamma = np.zeros((nd*nd, nd*nd), dtype=complex)
    for itype, coupling_matrix in interactions.items():
        spins = spin_table[itype]
        for a in range(norb):
            for b in range(norb):
                v = coupling_matrix[a, b]
                if abs(v) < 1e-15:
                    continue
                for (s1, s2, s3, s4), w in spins.items():
                    if itype == 'PairHop':
                        row = (s4*norb+b)*nd + (s3*norb+a)
                        col = (s1*norb+a)*nd + (s2*norb+b)
                    else:
                        row = (s4*norb+b)*nd + (s3*norb+b)
                        col = (s1*norb+a)*nd + (s2*norb+a)
                    Gamma[row, col] += v * w
    return Gamma


def build_sc_code(interactions, norb):
    """Build S/C as implemented in sc.py."""
    nd = norb * norb
    S = np.zeros((nd, nd), dtype=complex)
    C = np.zeros((nd, nd), dtype=complex)

    def _get(itype):
        return interactions.get(itype, np.zeros((norb, norb)))

    U_mat = _get("CoulombIntra")
    Up_mat = _get("CoulombInter")
    J_mat = _get("Hund")
    Jp_mat = _get("Exchange")
    I_mat = _get("Ising")
    PH_mat = _get("PairHop")

    for l1 in range(norb):
        for l2 in range(norb):
            idx12 = l1 * norb + l2
            for l3 in range(norb):
                for l4 in range(norb):
                    idx34 = l3 * norb + l4
                    s_val = 0.0
                    c_val = 0.0

                    if l1 == l2 == l3 == l4:
                        s_val = U_mat[l1, l1]
                        c_val = U_mat[l1, l1]
                    elif l1 == l3 and l2 == l4 and l1 != l2:
                        s_val = Up_mat[l1, l2] - I_mat[l1, l2]
                        c_val = (-Up_mat[l1, l2] + J_mat[l1, l2]
                                 - I_mat[l1, l2])
                    elif l1 == l2 and l3 == l4 and l1 != l3:
                        s_val = J_mat[l1, l3] - 2.0 * I_mat[l1, l3]
                        c_val = 2.0 * Up_mat[l1, l3] - J_mat[l1, l3]
                    elif l1 == l4 and l2 == l3 and l1 != l2:
                        s_val = Jp_mat[l1, l2] + PH_mat[l1, l2]
                        c_val = Jp_mat[l1, l2] + PH_mat[l1, l2]

                    S[idx12, idx34] = s_val
                    C[idx12, idx34] = c_val
    return S, C


def test_sc_via_rpa(name, interactions, norb, G_orb):
    """Test S/C by comparing Kuroki formula with full 4-index RPA."""
    nd = 2 * norb
    nd2 = norb * norb

    # Build paramagnetic spin-orbital chi0
    chi0_so = np.zeros((nd*nd, nd*nd), dtype=complex)
    for s in range(2):
        for la in range(norb):
            for lc in range(norb):
                for lb in range(norb):
                    for ld in range(norb):
                        row = (s*norb+la)*nd + (s*norb+lc)
                        col = (s*norb+lb)*nd + (s*norb+ld)
                        chi0_so[row, col] = G_orb[la, lb] * G_orb[ld, lc]

    # Full RPA
    Gamma = build_gamma_so(interactions, norb)
    I_so = np.eye(nd*nd)
    chi_full = np.linalg.solve(I_so + chi0_so @ Gamma, chi0_so)

    # Extract spin/charge chi
    chi_uu = np.zeros((nd2, nd2), dtype=complex)
    chi_ud = np.zeros((nd2, nd2), dtype=complex)
    for l1 in range(norb):
        for l2 in range(norb):
            for l3 in range(norb):
                for l4 in range(norb):
                    row = l1*nd + l2
                    col_uu = l3*nd + l4
                    col_dd = (norb+l3)*nd + (norb+l4)
                    chi_uu[l1*norb+l2, l3*norb+l4] = chi_full[row, col_uu]
                    chi_ud[l1*norb+l2, l3*norb+l4] = chi_full[row, col_dd]

    chi_s_full = chi_uu - chi_ud
    chi_c_full = chi_uu + chi_ud

    # Kuroki with sc.py S/C
    chi0_pair = np.zeros((nd2, nd2), dtype=complex)
    for l1 in range(norb):
        for l2 in range(norb):
            for l3 in range(norb):
                for l4 in range(norb):
                    chi0_pair[l1*norb+l2, l3*norb+l4] = G_orb[l1, l3] * G_orb[l4, l2]

    S, C = build_sc_code(interactions, norb)
    I_pair = np.eye(nd2)
    chi_s_kuroki = np.linalg.solve(I_pair - chi0_pair @ S, chi0_pair)
    chi_c_kuroki = np.linalg.solve(I_pair + chi0_pair @ C, chi0_pair)

    ds = np.max(np.abs(chi_s_full - chi_s_kuroki))
    dc = np.max(np.abs(chi_c_full - chi_c_kuroki))

    # Also compute exact S/C from inversion
    chi0_inv = np.linalg.inv(chi0_pair)
    S_exact = chi0_inv - np.linalg.inv(chi_s_full)
    C_exact = np.linalg.inv(chi_c_full) - chi0_inv

    ds_exact = np.max(np.abs(S - S_exact))
    dc_exact = np.max(np.abs(C - C_exact))

    print(f"  {name:30s}: |Δchi_s|={ds:.2e}, |Δchi_c|={dc:.2e}, "
          f"|ΔS_vs_exact|={ds_exact:.2e}, |ΔC_vs_exact|={dc_exact:.2e}")

    return ds, dc


def main():
    norb = 2
    np.random.seed(42)
    G_raw = np.random.randn(norb, norb) + 1j*np.random.randn(norb, norb)
    G_orb = (G_raw + G_raw.T.conj()) / (4 * norb)

    U_mat = np.diag([4.0, 4.0])
    V_mat = np.array([[0, 2.0], [2.0, 0]])
    J_mat = np.array([[0, 0.8], [0.8, 0]])
    Jp_mat = np.array([[0, 0.6], [0.6, 0]])
    I_mat = np.array([[0, 1.0], [1.0, 0]])

    print("S/C matrix verification: sc.py values vs full 4-index RPA")
    print("=" * 90)

    # Individual interactions
    test_sc_via_rpa("CoulombIntra U=4", {'CoulombIntra': U_mat}, norb, G_orb)
    test_sc_via_rpa("CoulombInter V=2", {'CoulombInter': V_mat}, norb, G_orb)
    test_sc_via_rpa("Hund J=0.8", {'Hund': J_mat}, norb, G_orb)
    test_sc_via_rpa("Exchange J'=0.6", {'Exchange': Jp_mat}, norb, G_orb)
    test_sc_via_rpa("Ising I=1.0", {'Ising': I_mat}, norb, G_orb)
    test_sc_via_rpa("PairLift P=1.0", {'PairLift': I_mat}, norb, G_orb)
    test_sc_via_rpa("PairHop P=1.0", {'PairHop': I_mat}, norb, G_orb)

    print()
    # Combined
    test_sc_via_rpa("Full Kanamori",
                    {'CoulombIntra': U_mat, 'CoulombInter': V_mat,
                     'Hund': J_mat, 'Exchange': Jp_mat}, norb, G_orb)
    test_sc_via_rpa("Kanamori + Ising",
                    {'CoulombIntra': U_mat, 'CoulombInter': V_mat,
                     'Hund': J_mat, 'Exchange': Jp_mat, 'Ising': I_mat}, norb, G_orb)
    test_sc_via_rpa("All interactions",
                    {'CoulombIntra': U_mat, 'CoulombInter': V_mat,
                     'Hund': J_mat, 'Exchange': Jp_mat,
                     'Ising': I_mat, 'PairLift': I_mat, 'PairHop': I_mat}, norb, G_orb)


if __name__ == '__main__':
    main()
