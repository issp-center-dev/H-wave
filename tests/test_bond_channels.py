import unittest

import numpy as np

from hwave.solver import bond_channels as bc
from tests.approx_util import ApproxTestCase


def _bond_topo(coulomb_inter, cell, norb, max_shells=None):
    """The shared bond topology of a CoulombInter declaration."""
    return bc.resolve_bond_topology({"CoulombInter": coulomb_inter}, cell, norb,
                                    max_shells=max_shells,
                                    active_types=("CoulombInter",))


def _bond_view(coulomb_inter, cell, norb, max_shells=None):
    """:class:`BondSetView` of :func:`_bond_topo` (the bubble's bond set)."""
    return bc.BondSetView(_bond_topo(coulomb_inter, cell, norb, max_shells))

def _nn_square(V=0.25):
    # single orbital, NN bonds ±x,±y on a square lattice
    return {((1,0,0),(0,0)): V, ((-1,0,0),(0,0)): V,
            ((0,1,0),(0,0)): V, ((0,-1,0),(0,0)): V}


def _toy_green_4x4(Nx=4, Ny=4, Nz=1, Nmat=8, t=1.0, mu=0.3, beta=2.0):
    """Single-orbital tight-binding k-space Green function on a tiny grid,
    in the sc.py layout (norb, norb, Nx, Ny, Nz, nmat) with the same
    Matsubara convention as hwave.sc._build_bond_green / rpa.py's _calc_green:
    iomega_n = (2n+1-nmat)*pi/beta.

    Returns (green, T, N) where T = 1/beta and N = Nx*Ny*Nz, matching the
    (T/N) prefactor convention of spec Eq.14.
    """
    kx = 2.0 * np.pi * np.arange(Nx) / Nx
    ky = 2.0 * np.pi * np.arange(Ny) / Ny
    KX, KY = np.meshgrid(kx, ky, indexing='ij')
    eps = -2.0 * t * (np.cos(KX) + np.cos(KY)) - mu  # (Nx,Ny)
    eps = np.broadcast_to(eps[:, :, np.newaxis], (Nx, Ny, Nz))

    iom = (2.0 * np.arange(Nmat) + 1.0 - Nmat) * np.pi / beta  # (Nmat,)

    green = np.empty((1, 1, Nx, Ny, Nz, Nmat), dtype=complex)
    green[0, 0] = 1.0 / (1j * iom[np.newaxis, np.newaxis, np.newaxis, :]
                         - eps[:, :, :, np.newaxis])

    T = 1.0 / beta
    N = Nx * Ny * Nz
    return green, T, N


def _chi0q_ref(green, T, N):
    """Ordinary (single-orbital) chi0q via a literal k-sum:
    chi0(q) = -(T/N) sum_{k,iwn} G(k+q,iwn) G(k,iwn).
    """
    _, _, Nx, Ny, Nz, Nmat = green.shape
    g = green[0, 0]  # (Nx,Ny,Nz,Nmat)
    chi0 = np.zeros((Nx, Ny, Nz), dtype=complex)
    for qx in range(Nx):
        for qy in range(Ny):
            for qz in range(Nz):
                s = 0.0 + 0.0j
                for kx in range(Nx):
                    for ky in range(Ny):
                        for kz in range(Nz):
                            kqx = (kx + qx) % Nx
                            kqy = (ky + qy) % Ny
                            kqz = (kz + qz) % Nz
                            s += np.sum(g[kqx, kqy, kqz, :] * g[kx, ky, kz, :])
                chi0[qx, qy, qz] = -(T / N) * s
    return chi0


def _toy_green_4x4_2orb(Nx=4, Ny=4, Nz=1, Nmat=4, t=1.0, mu1=0.3, mu2=-0.2,
                         t12=0.05 + 0.02j, beta=2.0):
    """Two-orbital k-space Green function on a tiny grid, same layout/
    convention as ``_toy_green_4x4`` but with a NON-symmetric inter-orbital
    2x2 Hamiltonian::

        H(k) = [[eps1(k), t12], [conj(t12), eps2(k)]]

    with distinct on-site dispersions (eps1 disperses along x, eps2 along y)
    and a small complex inter-orbital hopping ``t12``. This makes
    ``G_{01}(k) != G_{10}(k)`` (and both differ from ``G_{00}``/``G_{11}``),
    which is exactly what is needed to catch an index-transposition bug in
    ``bond_bubble``'s einsum/reshape that a diagonal (norb=1 or
    orbital-symmetric) Green function could not expose.

    G(k,iwn) = (iwn*I - H(k))^{-1}, computed via the explicit 2x2 inverse
    ``[[a,b],[c,d]]^{-1} = 1/(ad-bc) [[d,-b],[-c,a]]``.
    """
    kx = 2.0 * np.pi * np.arange(Nx) / Nx
    ky = 2.0 * np.pi * np.arange(Ny) / Ny
    KX, KY = np.meshgrid(kx, ky, indexing='ij')
    eps1 = -2.0 * t * np.cos(KX) - mu1
    eps2 = -2.0 * t * np.cos(KY) - mu2
    eps1 = np.broadcast_to(eps1[:, :, np.newaxis], (Nx, Ny, Nz))
    eps2 = np.broadcast_to(eps2[:, :, np.newaxis], (Nx, Ny, Nz))

    iom = (2.0 * np.arange(Nmat) + 1.0 - Nmat) * np.pi / beta  # (Nmat,)
    iomega = iom[np.newaxis, np.newaxis, np.newaxis, :]

    a = 1j * iomega - eps1[:, :, :, np.newaxis]
    d = 1j * iomega - eps2[:, :, :, np.newaxis]
    b = -t12
    c = -np.conj(t12)
    det = a * d - b * c

    green = np.empty((2, 2, Nx, Ny, Nz, Nmat), dtype=complex)
    green[0, 0] = d / det
    green[0, 1] = t12 / det
    green[1, 0] = np.conj(t12) / det
    green[1, 1] = a / det

    T = 1.0 / beta
    N = Nx * Ny * Nz
    return green, T, N


def _nn_square_2orb(V=0.25, Voff=0.05 + 0.02j):
    # 2-orbital NN bonds ±x,±y on a square lattice with a NON-symmetric
    # orbital-off-diagonal V_ab(R) on the ±x bonds (idx=l1*norb+l2 with
    # l1 != l2), so an index transposition changes the result.
    return {
        ((1, 0, 0), (0, 0)): V,        ((1, 0, 0), (1, 1)): 0.6 * V,
        ((1, 0, 0), (0, 1)): Voff,
        ((-1, 0, 0), (0, 0)): V,       ((-1, 0, 0), (1, 1)): 0.6 * V,
        ((-1, 0, 0), (1, 0)): np.conj(Voff),
        ((0, 1, 0), (0, 0)): V,        ((0, 1, 0), (1, 1)): 0.6 * V,
        ((0, -1, 0), (0, 0)): V,       ((0, -1, 0), (1, 1)): 0.6 * V,
    }


def _direct_bond_bubble(green, bond_set, T, N):
    """Literal nested-loop evaluation of spec Eq.14 (S4.2):

    chi_bar[(m,idx),(m',idx')](q) =
        -(T/N) sum_{k,iwn} e^{i k.(Delta r_m - Delta r_m')}
                            G_{l1,l3}(k+q,iwn) G_{l4,l2}(k,iwn)

    with idx = l1*norb+l2, idx' = l3*norb+l4 (sc.py bond-major convention,
    S4.1). Deliberately independent of the FFT machinery in bond_bubble.
    """
    norb, norb2, Nx, Ny, Nz, Nmat = green.shape
    assert norb2 == norb
    npair = norb * norb
    B = bond_set.n_channels
    nd_enlarged = npair * B

    kx_list = 2.0 * np.pi * np.arange(Nx) / Nx
    ky_list = 2.0 * np.pi * np.arange(Ny) / Ny
    kz_list = 2.0 * np.pi * np.arange(Nz) / Nz

    chi_bar = np.zeros((Nx, Ny, Nz, nd_enlarged, nd_enlarged), dtype=complex)

    for m in range(B):
        drm = bond_set.delta_r[m]
        for mp in range(B):
            drmp = bond_set.delta_r[mp]
            ddx = drm[0] - drmp[0]
            ddy = drm[1] - drmp[1]
            ddz = drm[2] - drmp[2]

            for qx in range(Nx):
                for qy in range(Ny):
                    for qz in range(Nz):
                        for l1 in range(norb):
                            for l2 in range(norb):
                                idx = l1 * norb + l2
                                I = m * npair + idx
                                for l3 in range(norb):
                                    for l4 in range(norb):
                                        idxp = l3 * norb + l4
                                        Ip = mp * npair + idxp

                                        s = 0.0 + 0.0j
                                        for kx in range(Nx):
                                            for ky in range(Ny):
                                                for kz in range(Nz):
                                                    phase = np.exp(1j * (
                                                        kx_list[kx] * ddx
                                                        + ky_list[ky] * ddy
                                                        + kz_list[kz] * ddz))
                                                    kqx = (kx + qx) % Nx
                                                    kqy = (ky + qy) % Ny
                                                    kqz = (kz + qz) % Nz
                                                    Gk_q = green[l1, l3, kqx, kqy, kqz, :]
                                                    Gk = green[l4, l2, kx, ky, kz, :]
                                                    s += phase * np.sum(Gk_q * Gk)
                                        chi_bar[qx, qy, qz, I, Ip] = -(T / N) * s
    return chi_bar


class TestBondTopologyAndBondBubble(ApproxTestCase):

    def test_bond_bubble_matches_direct_sum_tiny_grid(self):
        # 4x4x1 single-orbital tight-binding G; compare FFT bond bubble to
        # explicit k-sum of Eq.14 for a chosen (q, Δr, Δr').  (helper builds G.)
        from hwave.solver.bond_channels import bond_bubble
        green, T, N = _toy_green_4x4()                      # helper in test file
        beta = 1.0 / T
        bset = _bond_view(_nn_square(0.25), np.eye(3), norb=1)
        chib = bond_bubble(green, bset, beta=beta)          # (4,4,1,5,5)
        ref  = _direct_bond_bubble(green, bset, T, N)       # explicit Eq.14 loop
        assert chib.shape == (4, 4, 1, 5, 5)
        assert np.allclose(chib, ref, atol=1e-10)
        # m=0 block equals ordinary chi0q
        assert np.allclose(chib[..., 0, 0], _chi0q_ref(green, T, N), atol=1e-10)

    def test_bond_bubble_matches_direct_sum_tiny_grid_2orb(self):
        # Same oracle as test_bond_bubble_matches_direct_sum_tiny_grid but at
        # norb=2 with a non-symmetric inter-orbital Green function and an
        # orbital-off-diagonal bond interaction, so that a transposition bug in
        # the einsum string or the moveaxis/reshape block extraction (which the
        # norb=1 case cannot see, since l1=l2=l3=l4=0 there) would be caught.
        from hwave.solver.bond_channels import bond_bubble
        green, T, N = _toy_green_4x4_2orb(Nx=4, Ny=4, Nz=1, Nmat=4)
        beta = 1.0 / T
        bset = _bond_view(_nn_square_2orb(), np.eye(3), norb=2)
        assert bset.n_channels == 5
        chib = bond_bubble(green, bset, beta=beta)
        ref = _direct_bond_bubble(green, bset, T, N)
        npair = 2 * 2
        nd = npair * bset.n_channels
        assert chib.shape == (4, 4, 1, nd, nd)
        assert np.allclose(chib, ref, atol=1e-10)
        # m=0 block (idx,idxp both l1=l2=0 orbital-diagonal slice) equals the
        # ordinary multi-orbital chi0q for the (0,0) orbital pair.
        assert np.allclose(chib[..., 0, 0], _chi0q_ref(green[0:1, 0:1], T, N),
                            atol=1e-10)

    def test_square_nn_gives_five_channels(self):
        r = _bond_view(_nn_square(), np.eye(3), norb=1)
        assert r.n_channels == 5           # {0,+x,-x,+y,-y}
        assert r.delta_r[0] == (0,0,0)     # Δr=0 first
        # reverse map pairs +R with -R
        for m, dr in enumerate(r.delta_r):
            rm = r.delta_r[r.reverse[m]]
            assert tuple(-np.array(dr)) == rm

    def test_reverse_synthesis_from_plus_R_only(self):
        ci = {((1,0,0),(0,0)): 0.3}                  # +x only
        topo = _bond_topo(ci, np.eye(3), norb=1)
        dr = [tuple(int(x) for x in r) for r in topo.delta_r]
        assert len(dr) == 3                           # {0,+x,-x} synthesized
        v = topo.coeffs["CoulombInter"]
        # channel m files the row declared at -delta_r[m] (issue #232): the
        # declared +x row sits on the -x channel, its synthesized partner on +x
        self.assertApprox(v[dr.index((-1,0,0))][0,0], 0.3)
        self.assertApprox(v[dr.index((1,0,0))][0,0], 0.3)  # V_ab(R)=V_ba(-R)* real

    def test_v0_declared_topology_kept(self):
        ci = {((1,0,0),(0,0)): 0.0, ((0,1,0),(0,0)): 0.0}  # declared, zero-valued
        r = _bond_view(ci, np.eye(3), norb=1)
        assert r.n_channels == 5           # topology declared -> B stays 5 at V=0

    def test_cutoff_shell_zero_with_nonzero_V_errors(self):
        with self.assertRaises(ValueError):
            _bond_view(_nn_square(0.25), np.eye(3), norb=1, max_shells=0)

    def test_bond_reverse_mismatch_raises_for_nonzero_bond(self):
        """When BOTH V_01(+R) and V_10(-R) are declared but disagree with the
        Hermiticity constraint V_ab(R) = conj(V_ba(-R)) beyond tolerance, the
        resolver must raise -- not silently pick one side or average them."""
        norb = 2
        # V_01(+R) = 0.5+0.1j : declared. Consistency requires V_10(-R) =
        # conj(V_01(+R)) = 0.5-0.1j, but 0.1+0.5j is declared instead --
        # disagreement is not a rounding-level effect (|diff| = 0.72).
        ci = {((1, 0, 0), (0, 1)): 0.5 + 0.1j, ((-1, 0, 0), (1, 0)): 0.1 + 0.5j}
        with self.assertRaisesRegex(ValueError, r"resolve_bond_topology.*disagrees"):
            _bond_view(ci, np.eye(3), norb=norb)


# ---------------------------------------------------------------------------
# dress_bond -- the RPA-ladder dressing at the enlarged (ND x ND) index (spec
# S4.4); the non-bond static path calls it at ND = nd.
# ---------------------------------------------------------------------------

class TestDressBond(unittest.TestCase):

    def test_dress_bond_reduces_to_sc_dressing_at_B1(self):
        """Spec S4.4, B=1 reduction: at B=1 (ND=nd), dress_bond's batched solves
        must equal the existing sc.py inline dressing
        (sc.py:1098-1104, _compute_vertices_general) bit-for-bit (same algebra,
        same orientation: chi_s = solve(I - chi_bar@S, chi_bar),
        chi_c = solve(I + chi_bar@C, chi_bar))."""
        from hwave.solver.bond_channels import dress_bond

        rng = np.random.default_rng(1234)
        Nx, Ny, Nz = 2, 2, 1
        norb = 2
        nd = norb * norb  # ND = nd at B=1

        # random Hermitian chi_bar per q-point (a physical susceptibility block)
        A = (rng.normal(size=(Nx, Ny, Nz, nd, nd))
             + 1j * rng.normal(size=(Nx, Ny, Nz, nd, nd)))
        chi_bar = 0.5 * (A + np.conj(np.swapaxes(A, -1, -2)))

        # random real S, C (the Kuroki S/C matrices are real)
        S_bond = rng.normal(size=(Nx, Ny, Nz, nd, nd)).astype(complex)
        C_bond = rng.normal(size=(Nx, Ny, Nz, nd, nd)).astype(complex)

        chi_s, chi_c = dress_bond(chi_bar, S_bond, C_bond)
        assert chi_s.shape == (Nx, Ny, Nz, nd, nd)
        assert chi_c.shape == (Nx, Ny, Nz, nd, nd)

        # existing sc.py inline dressing, copied verbatim (sc.py:1098-1104)
        I_mat = np.broadcast_to(np.eye(nd, dtype=complex), (Nx, Ny, Nz, nd, nd)).copy()
        mat_s = I_mat - chi_bar @ S_bond
        mat_c = I_mat + chi_bar @ C_bond
        chi_s_ref = np.linalg.solve(mat_s, chi_bar)
        chi_c_ref = np.linalg.solve(mat_c, chi_bar)

        np.testing.assert_allclose(chi_s, chi_s_ref, atol=1e-12, rtol=1e-10)
        np.testing.assert_allclose(chi_c, chi_c_ref, atol=1e-12, rtol=1e-10)

    def test_dress_bond_enlarged_shape_and_mixes_bond_sectors(self):
        """dress_bond at B>1 returns full (Nx,Ny,Nz,ND,ND) arrays -- chi_bar is not
        block-diagonal in the bond index, so chi_s/chi_c generally have nonzero
        off-bond-block entries even though S_bond/C_bond are block-diagonal."""
        from hwave.solver.bond_channels import dress_bond

        rng = np.random.default_rng(7)
        Nx, Ny, Nz = 2, 1, 1
        norb = 1
        nd = norb * norb
        B = 3
        ND = nd * B

        A = (rng.normal(size=(Nx, Ny, Nz, ND, ND))
             + 1j * rng.normal(size=(Nx, Ny, Nz, ND, ND)))
        chi_bar = 0.5 * (A + np.conj(np.swapaxes(A, -1, -2)))
        # block-diagonal-in-m S_bond/C_bond, but off-diagonal chi_bar
        S_bond = np.zeros((Nx, Ny, Nz, ND, ND), dtype=complex)
        C_bond = np.zeros((Nx, Ny, Nz, ND, ND), dtype=complex)
        for m in range(B):
            S_bond[:, :, :, m, m] = 0.3 + 0.1 * m
            C_bond[:, :, :, m, m] = 0.2 - 0.05 * m

        chi_s, chi_c = dress_bond(chi_bar, S_bond, C_bond)
        assert chi_s.shape == (Nx, Ny, Nz, ND, ND)
        assert chi_c.shape == (Nx, Ny, Nz, ND, ND)
        # off-bond-block entries are generally nonzero (chi_bar mixes sectors)
        assert not np.allclose(chi_s[..., 0, 1], 0.0)


# ---------------------------------------------------------------------------
# The static bond Eliashberg operator (``sc._bond_kernel_parts``: fluctuation
# + instantaneous / bare part, built from the bond pairing pieces of
# ``eliashberg_bond``), spec S4.5. Reference-kernel oracle (S7.5) and the
# analytic single-band triplet fixture lambda_t = -1 (S7.9).
# ---------------------------------------------------------------------------

def _g2_direct(green, beta):
    """Independent literal G2 = (1/beta) sum_n G(k,wn) G(-k,-wn), matching the
    sc.py Matsubara/-k convention (iomega_n = (2n+1-nmat)*pi/beta, so the -wn
    partner of index n is index nmat-1-n, and -k is (N-k)%N).  Deliberately
    NOT the vectorized sc._calc_g2 -- a literal nested loop so the oracle is
    fully independent of the implementation under test.
    """
    norb, _, Nx, Ny, Nz, nmat = green.shape
    G2 = np.zeros((norb, norb, norb, norb, Nx, Ny, Nz), dtype=complex)
    for kx in range(Nx):
        for ky in range(Ny):
            for kz in range(Nz):
                krx, kry, krz = (-kx) % Nx, (-ky) % Ny, (-kz) % Nz
                for n in range(nmat):
                    nr = nmat - 1 - n
                    for i in range(norb):
                        for j in range(norb):
                            for l in range(norb):
                                for m in range(norb):
                                    G2[i, j, l, m, kx, ky, kz] += (
                                        green[i, j, kx, ky, kz, n]
                                        * green[l, m, krx, kry, krz, nr])
    return G2 / beta


def _parity_projected(phi, pairing_type):
    """Literal ``0.5 (phi_ab(k) +/- phi_ba(-k))`` (+ singlet, - triplet)."""
    norb, _, Nx, Ny, Nz = phi.shape
    sign = 1.0 if pairing_type == "singlet" else -1.0
    out = np.zeros_like(phi)
    for a in range(norb):
        for b in range(norb):
            for kx in range(Nx):
                for ky in range(Ny):
                    for kz in range(Nz):
                        out[a, b, kx, ky, kz] = 0.5 * (
                            phi[a, b, kx, ky, kz]
                            + sign * phi[b, a, (-kx) % Nx, (-ky) % Ny, (-kz) % Nz])
    return out


def _direct_bond_kernel(F_q, W0_q, green, bond_set, beta, phi, pairing_type):
    """Literal O(N_k^2) evaluation of the two-part bond Eliashberg operator
    (spec S4.5), independent of the production kernel.

    Returns ``(K^fl phi, K^pp phi)`` with

      K^fl phi(k)_{ad} = -(1/N) sum_{k'} sum_{bc} Veff_fl(k,k')_{abcd} X_bc(k')
      K^pp phi(k)_{ad} = -(1/N) sum_{k'} sum_{bc} W0(k-k')_{abcd} X^P_bc(k')

    where
      X_cd(k')      = sum_{ef} G2[c,e,d,f](k') phi_{ef}(k')      (cd = c*norb+d)
      Veff_fl_{abcd}(k,k') = sum_{m,m'} e^{-i k'.dr_m} F_q(k-k')_{(m,ab),(m',cd)}
                                          e^{-i k.dr_m'}
    and ``X^P`` the same pair bubble of the gap projected onto the parity
    sector of the pairing channel. The fluctuation part crosses the orbitals
    of the ph-layout F (``Y_ad <- F_(ab),(cd) X_bc``) and takes its input leg
    from the ROW channel and its output leg from the COLUMN channel; every
    pair leg carries ``e^{-ik.dr}``. This is the pair geometry issue #232
    fixed (a declared two-body row on the bond of its Transfer row), proven
    against an independent spin-orbital reference in
    tests/test_bond_pairing_placement.py. The bare part is the channel-0
    pair matrix ``W0 = 0.5 (S_0 + C_0)`` (singlet) / ``0.5 (C_0 - S_0)``
    (triplet), q = k - k' dependent through the off-site Hartree.
    """
    norb = green.shape[0]
    nd = norb * norb
    B = bond_set.n_channels
    Nx, Ny, Nz, nmat = green.shape[2:6]
    N = Nx * Ny * Nz
    dr = bond_set.delta_r

    G2 = _g2_direct(green, beta)

    def _bubble(p):
        X = np.zeros((nd, Nx, Ny, Nz), dtype=complex)
        for kx in range(Nx):
            for ky in range(Ny):
                for kz in range(Nz):
                    for c in range(norb):
                        for d in range(norb):
                            s = 0.0 + 0.0j
                            for e in range(norb):
                                for f in range(norb):
                                    s += (G2[c, e, d, f, kx, ky, kz]
                                          * p[e, f, kx, ky, kz])
                            X[c * norb + d, kx, ky, kz] = s
        return X

    X = _bubble(phi)
    XP = _bubble(_parity_projected(phi, pairing_type))

    kx_list = 2.0 * np.pi * np.arange(Nx) / Nx
    ky_list = 2.0 * np.pi * np.arange(Ny) / Ny
    kz_list = 2.0 * np.pi * np.arange(Nz) / Nz

    def kdot(ix, iy, iz, d):
        return kx_list[ix] * d[0] + ky_list[iy] * d[1] + kz_list[iz] * d[2]

    out_fl = np.zeros((nd, Nx, Ny, Nz), dtype=complex)
    out_pp = np.zeros((nd, Nx, Ny, Nz), dtype=complex)
    for kx in range(Nx):
        for ky in range(Ny):
            for kz in range(Nz):
                acc_fl = np.zeros(nd, dtype=complex)
                acc_pp = np.zeros(nd, dtype=complex)
                for px in range(Nx):
                    for py in range(Ny):
                        for pz in range(Nz):
                            qx, qy, qz = (kx - px) % Nx, (ky - py) % Ny, (kz - pz) % Nz
                            Fblk = F_q[qx, qy, qz]
                            Veff = np.zeros((norb,) * 4, dtype=complex)
                            for m in range(B):
                                e_in = np.exp(-1j * kdot(px, py, pz, dr[m]))
                                for mp in range(B):
                                    e_out = np.exp(-1j * kdot(kx, ky, kz, dr[mp]))
                                    Veff += e_in * e_out * Fblk[
                                        m * nd:(m + 1) * nd,
                                        mp * nd:(mp + 1) * nd].reshape((norb,) * 4)
                            acc_fl += np.einsum(
                                "abcd,bc->ad", Veff,
                                X[:, px, py, pz].reshape(norb, norb)).reshape(nd)
                            acc_pp += np.einsum(
                                "abcd,bc->ad",
                                W0_q[qx, qy, qz].reshape((norb,) * 4),
                                XP[:, px, py, pz].reshape(norb, norb)).reshape(nd)
                out_fl[:, kx, ky, kz] = -(1.0 / N) * acc_fl
                out_pp[:, kx, ky, kz] = -(1.0 / N) * acc_pp
    shape = (norb, norb, Nx, Ny, Nz)
    return out_fl.reshape(shape), out_pp.reshape(shape)


def _fluctuation_vertex(pairing_type, S_bond, chi_s, C_bond, chi_c):
    """F_eta(q) per spec S4.5 (batched matmul over q)."""
    SxS = S_bond @ chi_s @ S_bond
    CxC = C_bond @ chi_c @ C_bond
    if pairing_type == "singlet":
        return 1.5 * SxS - 0.5 * CxC
    return -0.5 * SxS - 0.5 * CxC


def _bare_channel0(pairing_type, S_bond, C_bond, nd):
    """The channel-0 bare pair matrix ``0.5 (S_0 +/- C_0)`` per q."""
    S0 = S_bond[..., :nd, :nd]
    C0 = C_bond[..., :nd, :nd]
    if pairing_type == "singlet":
        return 0.5 * (S0 + C0)
    return 0.5 * (C0 - S0)


class TestStaticBondKernel(unittest.TestCase):

    def _check_oracle(self, rng, green, bset, norb, shape, beta, n_vec):
        import hwave.sc as sc
        Nx, Ny, Nz = shape
        nd = norb * norb
        ND = nd * bset.n_channels

        def _herm(shp):
            A = rng.normal(size=shp) + 1j * rng.normal(size=shp)
            return 0.5 * (A + np.conj(np.swapaxes(A, -1, -2)))

        chi_s = _herm((Nx, Ny, Nz, ND, ND))
        chi_c = _herm((Nx, Ny, Nz, ND, ND))
        S_bond = rng.normal(size=(Nx, Ny, Nz, ND, ND)).astype(complex)
        C_bond = rng.normal(size=(Nx, Ny, Nz, ND, ND)).astype(complex)

        for pairing_type in ("singlet", "triplet"):
            with self.subTest(pairing_type=pairing_type):
                A, A_fl, A_pp, vec_size = sc._bond_kernel_parts(
                    S_bond, C_bond, [("spin", chi_s), ("charge", chi_c)],
                    green, bset, norb, shape, beta, pairing_type)
                assert vec_size == norb * norb * Nx * Ny * Nz
                F_q = _fluctuation_vertex(pairing_type, S_bond, chi_s, C_bond, chi_c)
                W0_q = _bare_channel0(pairing_type, S_bond, C_bond, nd)
                for _ in range(n_vec):
                    phi = (rng.normal(size=(norb, norb, Nx, Ny, Nz))
                           + 1j * rng.normal(size=(norb, norb, Nx, Ny, Nz)))
                    ref_fl, ref_pp = _direct_bond_kernel(
                        F_q, W0_q, green, bset, beta, phi, pairing_type)
                    got_fl = A_fl.matvec(phi.ravel()).reshape(phi.shape)
                    got_pp = A_pp.matvec(phi.ravel()).reshape(phi.shape)
                    got = A.matvec(phi.ravel()).reshape(phi.shape)
                    np.testing.assert_allclose(got_fl, ref_fl, atol=1e-10)
                    np.testing.assert_allclose(got_pp, ref_pp, atol=1e-10)
                    np.testing.assert_allclose(got, ref_fl + ref_pp, atol=1e-10)

    def test_bond_kernel_reference_oracle_singleband(self):
        """Spec S7.5: the operator equals a direct O(N_k^2) evaluation of the
        two-part V_eff(k,k') on a 4x4 grid, pinning the phase convention of
        the fluctuation part (e^{-i(k'.dr_m + k.dr_m')} with the row channel m
        on the input leg; issue #232), the q=k-k' argument of both parts, the
        parity projection of the bare part and the -(1/N) normalization."""
        rng = np.random.default_rng(20260726)
        norb = 1
        beta = 2.0
        green, T, N = _toy_green_4x4(Nx=4, Ny=4, Nz=1, Nmat=6, beta=beta)
        bset = _bond_view(_nn_square(0.3), np.eye(3), norb=norb)
        self._check_oracle(rng, green, bset, norb, (4, 4, 1), beta, 3)

    def test_bond_kernel_reference_oracle_two_orbital(self):
        """Spec S7.5, norb=2: the same oracle at two orbitals with a non-symmetric
        Green function and orbital-off-diagonal bonds, so an index transposition
        (l1<->l2 / row<->col in F or the bare block, or in the X = GGphi
        contraction) that the norb=1 case cannot expose is caught.

        Uses a 4x2x1 grid (bonds run along +-x/+-y): with >=3 points along the
        x axis, k_x takes values other than {0, pi}, so the bond phases
        e^{+-ik.dr} are genuinely complex (not just +-1). A 2x2x1 grid (k in
        {0, pi} on every axis) makes every bond phase real and the kernel
        invariant under flipping either phase sign, so it cannot catch a
        phase-sign bug coupled to the orbital-index handling."""
        rng = np.random.default_rng(11)
        norb = 2
        beta = 2.0
        green, T, N = _toy_green_4x4_2orb(Nx=4, Ny=2, Nz=1, Nmat=4, beta=beta)
        bset = _bond_view(_nn_square_2orb(), np.eye(3), norb=norb)
        self._check_oracle(rng, green, bset, norb, (4, 2, 1), beta, 2)

    def test_bond_kernel_triplet_lambda_minus_one(self):
        """Spec S7.9: single-band 1D 4-point grid, U=4, V=1, chi_s=chi_c=0 (only the
        instantaneous part), (T/N) G(k) G(-k) = 1/4.  The odd gap Delta=(0,1,0,-1)
        (= sin k) is an eigenvector with lambda_t = -1; the fluctuation piece is
        zero, so lambda_t^pp = -1 and lambda_t^fl = 0. The vertices come from
        the production chain (topology, locality split, bond vertices)."""
        import hwave.sc as sc

        norb = 1
        Nx, Ny, Nz = 4, 1, 1
        nmat = 1
        beta = 1.0
        # green == 1 everywhere, nmat=1  =>  G2[k] = (1/beta) * 1 * 1 = 1 for all k,
        # so X = G2 phi = phi and the 1/N (=1/4) fold gives (T/N)G G = 1/4.
        green = np.ones((norb, norb, Nx, Ny, Nz, nmat), dtype=complex)

        it = {"CoulombIntra": {((0, 0, 0), (0, 0)): 4.0},
              "CoulombInter": {((1, 0, 0), (0, 0)): 1.0, ((-1, 0, 0), (0, 0)): 1.0}}
        topo, view, split = sc._resolve_bond_topology(it, {"rvec": np.eye(3)}, norb, None)
        assert view.n_channels == 3
        S_bond, C_bond = sc._build_bond_vertices(topo, split, norb, (Nx, Ny, Nz))
        chi_zero = np.zeros_like(S_bond)

        A, A_fl, A_pp, vec_size = sc._bond_kernel_parts(
            S_bond, C_bond, [("spin", chi_zero), ("charge", chi_zero)],
            green, view, norb, (Nx, Ny, Nz), beta, "triplet")

        # odd gap Delta = sin k on k in {0, pi/2, pi, -pi/2}
        kx = 2.0 * np.pi * np.arange(Nx) / Nx
        phi = np.sin(kx).reshape(norb, norb, Nx, Ny, Nz).astype(complex)

        out = A.matvec(phi.ravel()).reshape(norb, norb, Nx, Ny, Nz)
        # K phi = -phi  ->  lambda_t = -1
        np.testing.assert_allclose(out, -phi, atol=1e-12)

        # attribution: Rayleigh quotients of the two parts on the (real)
        # eigenvector -- lambda_t^fl = 0, lambda_t^pp = -1.
        v = phi.ravel()
        denom = np.vdot(v, v)
        lam_fl = (np.vdot(v, A_fl.matvec(v)) / denom).real
        lam_pp = (np.vdot(v, A_pp.matvec(v)) / denom).real
        np.testing.assert_allclose(lam_fl, 0.0, rtol=0, atol=1e-12)
        np.testing.assert_allclose(lam_pp, -1.0, rtol=0, atol=1e-12)

    def test_bond_kernel_rejects_unknown_pairing_type(self):
        import hwave.sc as sc
        norb = 1
        beta = 2.0
        green, T, N = _toy_green_4x4(Nx=4, Ny=4, Nz=1, Nmat=4, beta=beta)
        bset = _bond_view(_nn_square(0.3), np.eye(3), norb=norb)
        ND = bset.n_channels
        z = np.zeros((4, 4, 1, ND, ND), dtype=complex)
        with self.assertRaises(ValueError):
            sc._bond_kernel_parts(z, z, [("spin", z), ("charge", z)], green, bset,
                                  norb, (4, 4, 1), beta, "nonsense")


# ---------------------------------------------------------------------------
# Boundary coverage of the resolver and of the enlarged RPA denominators
# (Codex phase-1 review: singular / near-singular I - chi_bar S and
# I + chi_bar C, an empty declared CoulombInter, the multi-shell cutoff, and a
# non-real on-site diagonal).
# ---------------------------------------------------------------------------

def _diag_denominator_fixture(bad_value, channel, Nx=2):
    """``(chi_bar, S_bond, C_bond)`` at ND = 2 whose ``channel`` denominator is
    ``diag(1, bad_value)`` at q = (1, 0, 0) and the identity everywhere else.

    ND must be >= 2 for the relative sigma_min/sigma_max criterion to mean
    anything: a 1x1 block always has ratio 1.
    """
    ND = 2
    shape = (Nx, 1, 1, ND, ND)
    chi_bar = np.zeros(shape, dtype=complex)
    S_bond = np.zeros(shape, dtype=complex)
    C_bond = np.zeros(shape, dtype=complex)
    for ix in range(Nx):
        chi_bar[ix, 0, 0] = np.eye(ND, dtype=complex)
    if channel == "spin":
        # I - chi_bar S: entry (1,1) becomes 1 - (1 - bad) = bad
        S_bond[1, 0, 0, 1, 1] = 1.0 - bad_value
    else:
        # I + chi_bar C: entry (1,1) becomes 1 - (1 - bad) = bad
        C_bond[1, 0, 0, 1, 1] = -(1.0 - bad_value)
    return chi_bar, S_bond, C_bond


class TestDressBondConditioningGuard(unittest.TestCase):

    def test_dress_bond_rejects_an_exactly_singular_denominator(self):
        """An exactly singular ``I -/+ chi_bar S/C`` must raise an ACTIONABLE
        error naming the channel and the q-point, not a bare LinAlgError."""
        from hwave.solver.bond_channels import dress_bond
        for channel in ["spin", "charge"]:
            with self.subTest(channel=channel):
                chi_bar, S_bond, C_bond = _diag_denominator_fixture(0.0, channel)
                with self.assertRaises(ValueError) as cm:
                    dress_bond(chi_bar, S_bond, C_bond)
                msg = str(cm.exception)
                assert channel in msg
                assert "(1, 0, 0)" in msg
                assert "Singular matrix" not in msg

    def test_dress_bond_rejects_a_nearly_singular_denominator(self):
        """A NEARLY singular denominator silently produces enormous, unreliable
        vertices; it must be refused with the conditioning ratio in the message."""
        from hwave.solver.bond_channels import dress_bond
        for channel in ["spin", "charge"]:
            with self.subTest(channel=channel):
                chi_bar, S_bond, C_bond = _diag_denominator_fixture(1.0e-4, channel)
                with self.assertRaises(ValueError) as cm:
                    dress_bond(chi_bar, S_bond, C_bond)
                msg = str(cm.exception)
                assert channel in msg
                assert "1.000e-04" in msg or "1.0e-04" in msg

    def test_dress_bond_accepts_a_well_conditioned_denominator(self):
        from hwave.solver.bond_channels import dress_bond
        for channel in ["spin", "charge"]:
            with self.subTest(channel=channel):
                chi_bar, S_bond, C_bond = _diag_denominator_fixture(0.5, channel)
                chi_s, chi_c = dress_bond(chi_bar, S_bond, C_bond)
                assert chi_s.shape == chi_bar.shape and chi_c.shape == chi_bar.shape

    def test_bond_conditioning_score_is_the_value_the_guard_publishes(self):
        """The scoring function and the raising wrapper must not drift: the
        warn-and-continue policy of GitHub issue #199 reports the score the
        refusal would have quoted."""
        import re
        from hwave.solver.bond_channels import (bond_conditioning_score,
                                                _check_bond_conditioning)
        chi_bar, S_bond, _C = _diag_denominator_fixture(0.5, "spin")
        mat = np.eye(chi_bar.shape[-1]) - chi_bar @ S_bond
        worst, iq, ratio_iq, pole_iq, smin_iq, smax_iq = bond_conditioning_score(mat)
        self.assertAlmostEqual(worst, _check_bond_conditioning("spin", mat, 1.0e-3), places=15)
        self.assertAlmostEqual(worst, min(ratio_iq, pole_iq), places=15)
        self.assertAlmostEqual(smin_iq / smax_iq, ratio_iq, places=12)
        chi_bar, S_bond, _C = _diag_denominator_fixture(1.0e-4, "spin")
        mat = np.eye(chi_bar.shape[-1]) - chi_bar @ S_bond
        worst, iq, ratio_iq, pole_iq, smin_iq, smax_iq = bond_conditioning_score(mat)
        with self.assertRaises(ValueError) as cm:
            _check_bond_conditioning("spin", mat, 1.0e-3)
        quoted = re.findall(r"= (-?[0-9.]+e[-+][0-9]+)", str(cm.exception))
        self.assertGreaterEqual(len(quoted), 2)
        self.assertAlmostEqual(min(float(quoted[0]), float(quoted[1])) / worst, 1.0, places=3)
        self.assertEqual(iq, 1)

    def test_dress_bond_conditioning_floor_is_tunable(self):
        """The floor is a knob, not a hard-coded refusal: an explicit lower
        ``cond_tol`` lets a deliberately stiff study through."""
        from hwave.solver.bond_channels import dress_bond
        chi_bar, S_bond, C_bond = _diag_denominator_fixture(1.0e-4, "spin")
        chi_s, _ = dress_bond(chi_bar, S_bond, C_bond, cond_tol=1.0e-8)
        assert np.all(np.isfinite(chi_s))


def _batch_denominator_fixture(bad_value, channel):
    """``(chi_bar_b, W)`` for :func:`dress_batch`: one frequency, two q
    points, ND = 2, whose ``channel`` denominator is ``diag(1, bad_value)``
    at q index 1 and the identity at q index 0."""
    chi_bar, S_bond, C_bond = _diag_denominator_fixture(bad_value, channel)
    cb = chi_bar.reshape(1, 2, 2, 2)                     # (nb=1, nvol=2, ND, ND)
    W = (S_bond if channel == "spin" else C_bond).reshape(2, 2, 2)
    return cb, W


class TestDressBatchGuardPolicy(unittest.TestCase):
    """``flex_guard_policy = "warn"`` (GitHub issue #199): a transient
    conditioning violation is logged and the batch continues."""

    def test_warn_logs_and_records_a_conditioning_violation(self):
        from hwave.solver import bond_channels as bc
        for channel in ("spin", "charge"):
            with self.subTest(channel=channel):
                cb, W = _batch_denominator_fixture(1.0e-4, channel)
                with self.assertRaises(ValueError):
                    bc.dress_batch(cb, W, channel, l0=0, nmat=2, spatial_shape=(2, 1, 1))
                violations = []
                with self.assertLogs("qlms.solver.bond_channels", level="WARNING") as cm:
                    chi, cond = bc.dress_batch(cb, W, channel, l0=0, nmat=2,
                                               spatial_shape=(2, 1, 1),
                                               guard_policy="warn", violations=violations)
                msg = "\n".join(cm.output)
                self.assertIn('flex_guard_policy = "warn"', msg)
                self.assertIn(channel, msg)
                self.assertIn("(1, 0, 0)", msg)
                self.assertTrue(np.all(np.isfinite(chi)))
                self.assertLess(cond, 1.0e-3)
                self.assertEqual(len(violations), 1)
                self.assertEqual(violations[0]["channel"], channel)
                self.assertEqual(violations[0]["kind"], "conditioning")
                self.assertEqual(violations[0]["q"], 1)
                self.assertEqual(violations[0]["l"], 0)
                self.assertLess(violations[0]["value"], 1.0e-3)

    def test_warn_does_not_change_a_well_conditioned_batch(self):
        from hwave.solver import bond_channels as bc
        cb, W = _batch_denominator_fixture(0.5, "spin")
        ref, c_ref = bc.dress_batch(cb, W, "spin", l0=0, nmat=2, spatial_shape=(2, 1, 1))
        violations = []
        out, c_out = bc.dress_batch(cb, W, "spin", l0=0, nmat=2, spatial_shape=(2, 1, 1),
                                    guard_policy="warn", violations=violations)
        np.testing.assert_array_equal(out, ref)
        self.assertEqual(c_out, c_ref)
        self.assertEqual(violations, [])

    def test_an_exactly_singular_block_still_raises_under_warn(self):
        """There is nothing finite to continue with: the policy tolerates a
        NEARLY singular denominator, not an exactly singular solve.

        Nothing may be reported as tolerated either -- the warning claims the
        dressing continues, and a run that ends here did not continue."""
        import logging
        from hwave.solver import bond_channels as bc
        cb, W = _batch_denominator_fixture(0.0, "spin")
        violations = []
        records = []
        logger = logging.getLogger("qlms.solver.bond_channels")
        h = logging.Handler(); h.emit = lambda rec: records.append(rec.getMessage())
        logger.addHandler(h)
        try:
            with self.assertRaises(ValueError) as cm:
                bc.dress_batch(cb, W, "spin", l0=0, nmat=2, spatial_shape=(2, 1, 1),
                               guard_policy="warn", violations=violations)
        finally:
            logger.removeHandler(h)
        self.assertIn("spin", str(cm.exception))
        self.assertEqual(violations, [])
        self.assertFalse([m for m in records if "continues" in m])

    def test_a_disabled_guard_scores_nothing_under_either_policy(self):
        """``cond_tol = None`` disables the conditioning guard itself: no
        score, no warning -- and an exactly singular solve still raises."""
        import logging
        from hwave.solver import bond_channels as bc
        for policy in ("refuse", "warn"):
            with self.subTest(policy=policy):
                cb, W = _batch_denominator_fixture(1.0e-4, "spin")
                violations = []
                records = []
                logger = logging.getLogger("qlms.solver.bond_channels")
                h = logging.Handler(); h.emit = lambda rec: records.append(rec.getMessage())
                logger.addHandler(h)
                try:
                    chi, cond = bc.dress_batch(cb, W, "spin", l0=0, nmat=2,
                                               spatial_shape=(2, 1, 1), cond_tol=None,
                                               guard_policy=policy, violations=violations)
                finally:
                    logger.removeHandler(h)
                self.assertIsNone(cond)
                self.assertTrue(np.all(np.isfinite(chi)))
                self.assertEqual(violations, [])
                self.assertEqual(records, [])
                cb, W = _batch_denominator_fixture(0.0, "spin")
                # with no guard the "refuse" path reaches the solver itself,
                # whose LinAlgError is not a ValueError in every numpy
                with self.assertRaises((ValueError, np.linalg.LinAlgError)):
                    bc.dress_batch(cb, W, "spin", l0=0, nmat=2, spatial_shape=(2, 1, 1),
                                   cond_tol=None, guard_policy=policy)

    def test_the_warnings_name_the_scf_iteration(self):
        import logging
        from hwave.solver import bond_channels as bc
        cb, W = _batch_denominator_fixture(1.0e-4, "spin")
        violations = []
        with self.assertLogs("qlms.solver.bond_channels", level="WARNING") as cm:
            bc.dress_batch(cb, W, "spin", l0=0, nmat=2, spatial_shape=(2, 1, 1),
                           guard_policy="warn", violations=violations, iteration=5)
        self.assertIn("(SCF iteration 5)", "\n".join(cm.output))
        self.assertEqual(violations[0]["iteration"], 5)
        # without one the message ends as before
        violations = []
        with self.assertLogs("qlms.solver.bond_channels", level="WARNING") as cm:
            bc.dress_batch(cb, W, "spin", l0=0, nmat=2, spatial_shape=(2, 1, 1),
                           guard_policy="warn", violations=violations)
        self.assertNotIn("SCF iteration", "\n".join(cm.output))
        self.assertIsNone(violations[0]["iteration"])
        # the residual guard of the reduced "static" mode says so too. Its
        # trigger -- the solve residual of a near-singular unchecked slice --
        # is BLAS-dependent (it can be exactly 0 on some LAPACK builds), so the
        # branch is exercised by fixing the residual monitor above the
        # tolerance rather than by a fixture that may or may not produce one.
        cb = np.zeros((1, 2, 2, 2), complex)
        Wr = np.zeros((2, 2, 2), complex)
        for q in range(2):
            cb[0, q] = np.eye(2)          # mat = identity: the solve is clean
        violations = []
        records = []
        logger = logging.getLogger("qlms.solver.bond_channels")
        h = logging.Handler(); h.emit = lambda rec: records.append(rec.getMessage())
        old_level = logger.level
        logger.addHandler(h); logger.setLevel(logging.WARNING)
        orig_residual = bc.solve_residual
        bc.solve_residual = lambda *a, **k: np.array([[0.0, 1.0e-3]])   # q = 1 over tol
        try:
            bc.dress_batch(cb, Wr, "spin", l0=0, nmat=2, spatial_shape=(2, 1, 1),
                           guard_freqs="static", residual_tol=1.0e-6,
                           guard_policy="warn", violations=violations, iteration=5)
        finally:
            bc.solve_residual = orig_residual
            logger.removeHandler(h); logger.setLevel(old_level)
        self.assertIn("(SCF iteration 5)", "\n".join(records))
        self.assertEqual([(v["kind"], v["iteration"]) for v in violations], [("residual", 5)])

    def test_invalid_guard_policy_is_refused(self):
        from hwave.solver import bond_channels as bc
        cb, W = _batch_denominator_fixture(0.5, "spin")
        for bad in ("ignore", "Warn", 1, None):
            with self.subTest(policy=bad):
                with self.assertRaises(ValueError) as cm:
                    bc.dress_batch(cb, W, "spin", l0=0, nmat=2, spatial_shape=(2, 1, 1),
                                   guard_policy=bad)
                self.assertIn("guard_policy", str(cm.exception))

    def test_warn_also_tolerates_the_static_residual_guard(self):
        """With guard_freqs = "static" an unchecked slice is judged by its
        solve residual; under "warn" that finding is logged, not raised, and
        under "refuse" it is raised.

        nmat = 2 puts the static slice at l = 1, outside this one-frequency
        batch, so only the residual guard speaks. The residual of a
        near-singular structured 2x2 solve is BLAS-dependent (it can be exactly
        0 on some LAPACK builds), so the guard's branch is exercised by fixing
        the residual monitor above the tolerance rather than by a fixture that
        may or may not produce one; the denominator itself is the identity, so
        the solve is clean and chi is finite."""
        import logging
        from hwave.solver import bond_channels as bc
        cb = np.zeros((1, 2, 2, 2), complex)
        W = np.zeros((2, 2, 2), complex)
        for q in range(2):
            cb[0, q] = np.eye(2)          # mat = identity: the solve is clean
        violations = []
        records = []
        logger = logging.getLogger("qlms.solver.bond_channels")
        h = logging.Handler(); h.emit = lambda rec: records.append(rec.getMessage())
        old_level = logger.level
        logger.addHandler(h); logger.setLevel(logging.WARNING)
        orig_residual = bc.solve_residual
        bc.solve_residual = lambda *a, **k: np.array([[0.0, 1.0e-3]])   # q = 1 over tol
        try:
            chi, _ = bc.dress_batch(cb, W, "spin", l0=0, nmat=2, spatial_shape=(2, 1, 1),
                                    guard_freqs="static", residual_tol=1.0e-6,
                                    guard_policy="warn", violations=violations)
            # and the same finding is refused under the default policy
            with self.assertRaises(ValueError):
                bc.dress_batch(cb, W, "spin", l0=0, nmat=2, spatial_shape=(2, 1, 1),
                               guard_freqs="static", residual_tol=1.0e-6)
        finally:
            bc.solve_residual = orig_residual
            logger.removeHandler(h); logger.setLevel(old_level)
        self.assertIn('flex_guard_policy = "warn"', "\n".join(records))
        self.assertTrue(np.all(np.isfinite(chi)))
        self.assertEqual([v["kind"] for v in violations], ["residual"])


# --- an empty declared CoulombInter ----------------------------------------

class TestEmptyDeclaredCoulombInter(unittest.TestCase):

    def test_resolver_on_an_empty_mapping(self):
        """A DECLARED but empty CoulombInter is legal (spec S5: a
        declared-but-zero V keeps the channel topology fixed across a V sweep);
        it must resolve to the on-site channel alone."""
        topo = _bond_topo({}, np.eye(3), norb=2)
        r = bc.BondSetView(topo)
        assert r.n_channels == 1
        assert r.delta_r == ((0, 0, 0),)
        assert r.reverse == (0,)
        assert np.array_equal(np.asarray(topo.coeffs["CoulombInter"]),
                              np.zeros((1, 2, 2)))

    def test_validate_bond_prereqs_accepts_a_declared_empty_coulomb_inter(self):
        import hwave.sc as sc
        sc._validate_bond_prereqs("calc", 1, {"CoulombInter": {}})


# --- FIX 3 (review round 3): non-finite CoulombInter values ----------------
#
# The real-interaction guard checked ``abs(complex(value).imag) > imag_tol``
# to reject complex interactions -- but NaN comparisons are always False and
# a real infinity has zero imaginary part, so NaN and +-inf entries used to
# pass the guard and reach the bond machinery instead of being rejected.

class TestNonFiniteCoulombInterRejected(unittest.TestCase):

    def test_validate_bond_prereqs_rejects_nan_coulomb_inter(self):
        import hwave.sc as sc
        ci = {((1, 0, 0), (0, 0)): float("nan")}
        with self.assertRaisesRegex(ValueError, "(?i)finite"):
            sc._validate_bond_prereqs("calc", 1, {"CoulombInter": ci})

    def test_validate_bond_prereqs_rejects_positive_infinity_coulomb_inter(self):
        import hwave.sc as sc
        ci = {((1, 0, 0), (0, 0)): float("inf")}
        with self.assertRaisesRegex(ValueError, "(?i)finite"):
            sc._validate_bond_prereqs("calc", 1, {"CoulombInter": ci})

    def test_validate_bond_prereqs_rejects_negative_infinity_coulomb_inter(self):
        import hwave.sc as sc
        ci = {((1, 0, 0), (0, 0)): float("-inf")}
        with self.assertRaisesRegex(ValueError, "(?i)finite"):
            sc._validate_bond_prereqs("calc", 1, {"CoulombInter": ci})

    def test_validate_bond_prereqs_rejects_nonfinite_imaginary_part(self):
        """A NaN/inf imaginary part alone (real part finite) must also be caught,
        not just a NaN/inf real part."""
        import hwave.sc as sc
        ci = {((1, 0, 0), (0, 0)): complex(0.3, float("nan"))}
        with self.assertRaisesRegex(ValueError, "(?i)finite"):
            sc._validate_bond_prereqs("calc", 1, {"CoulombInter": ci})

        ci2 = {((1, 0, 0), (0, 0)): complex(0.3, float("inf"))}
        with self.assertRaisesRegex(ValueError, "(?i)finite"):
            sc._validate_bond_prereqs("calc", 1, {"CoulombInter": ci2})

    def test_validate_bond_prereqs_error_names_offending_bond_and_value(self):
        """The error must name the offending (irvec, orbvec) and value, matching
        the style of the existing complex-value guard message."""
        import hwave.sc as sc
        ci = {((2, 0, 0), (0, 0)): float("nan")}
        with self.assertRaisesRegex(ValueError, r"\(2, 0, 0\)") as cm:
            sc._validate_bond_prereqs("calc", 1, {"CoulombInter": ci})
        assert "nan" in str(cm.exception).lower()


# --- the shell cutoff across MULTIPLE shells -------------------------------

def _two_shell_square(v1=0.3, v2=0.1):
    """NN (length 1) and NNN (length sqrt(2)) bonds on a square lattice."""
    ci = {}
    for R in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)):
        ci[(R, (0, 0))] = v1
    for R in ((1, 1, 0), (1, -1, 0), (-1, 1, 0), (-1, -1, 0)):
        ci[(R, (0, 0))] = v2
    return ci


class TestShellCutoffMultipleShells(unittest.TestCase):

    def test_max_shells_cuts_whole_shells(self):
        """``max_shells`` keeps whole |R|-shells, and refuses a cutoff that
        would drop a shell carrying declared nonzero content (it never
        discards declared coefficients silently)."""
        ci = _two_shell_square()
        nn_only = {(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)}
        nnn_only = {(1, 1, 0), (1, -1, 0), (-1, 1, 0), (-1, -1, 0)}

        full = _bond_view(ci, np.eye(3), norb=1)
        assert full.n_channels == 9                      # 1 + 4 + 4

        with self.assertRaisesRegex(ValueError, r"max_shells=1"):
            _bond_view(ci, np.eye(3), norb=1, max_shells=1)

        # the same cutoff is legal when the outer shell is declared at zero
        nn = _bond_view(_two_shell_square(v2=0.0), np.eye(3), norb=1,
                        max_shells=1)
        assert nn.n_channels == 5
        assert set(nn.delta_r[1:]) == nn_only
        # reversal closure survives the cutoff
        for m, dr in enumerate(nn.delta_r):
            assert nn.delta_r[nn.reverse[m]] == tuple(-x for x in dr)

        both = _bond_view(ci, np.eye(3), norb=1, max_shells=2)
        assert both.n_channels == 9
        assert set(both.delta_r[1:]) == nn_only | nnn_only

        # asking for more shells than exist is a no-op, not an error
        assert _bond_view(ci, np.eye(3), norb=1, max_shells=5).n_channels == 9


# ---------------------------------------------------------------------------
# Phase-1 round-2 review fixes:
#   R2-1  the conditioning guard must also measure the ABSOLUTE distance to the
#         pole (a 1x1 or a uniformly-small denominator has ratio == 1);
#   R2-2  a tolerance-ACCEPTED reverse pair must be projected onto an exactly
#         conjugate pair (the contract promises Hermiticity-closed output);
#   R2-3  max_shells validation (non-negative integral; the zero-shell
#         ambiguity guard keys on DECLARED nonzero V).
# ---------------------------------------------------------------------------

def _scalar_denominator_fixture(denom, channel, Nx=2):
    """``(chi_bar, S_bond, C_bond)`` at ND = 1 whose ``channel`` denominator is
    ``[[denom]]`` at q = (1, 0, 0) and ``[[1]]`` everywhere else.

    A 1x1 block has sigma_min/sigma_max == 1 for ANY nonzero value, so only an
    absolute pole-distance criterion can catch it.
    """
    ND = 1
    shape = (Nx, 1, 1, ND, ND)
    chi_bar = np.zeros(shape, dtype=complex)
    S_bond = np.zeros(shape, dtype=complex)
    C_bond = np.zeros(shape, dtype=complex)
    for ix in range(Nx):
        chi_bar[ix, 0, 0] = np.eye(ND, dtype=complex)
    if channel == "spin":
        S_bond[1, 0, 0, 0, 0] = 1.0 - denom
    else:
        C_bond[1, 0, 0, 0, 0] = -(1.0 - denom)
    return chi_bar, S_bond, C_bond


class TestPhase1RoundTwoReviewFixes(unittest.TestCase):

    def test_dress_bond_rejects_a_near_singular_scalar_denominator(self):
        """B = 1, norb = 1 (empty/local interaction set) gives ND = 1, where the
        relative sigma_min/sigma_max ratio is identically 1. A denominator of
        1e-12 still divides the vertices by 1e-12 and must be refused."""
        from hwave.solver.bond_channels import dress_bond
        for channel in ["spin", "charge"]:
            with self.subTest(channel=channel):
                chi_bar, S_bond, C_bond = _scalar_denominator_fixture(1.0e-12, channel)
                with self.assertRaises(ValueError) as cm:
                    dress_bond(chi_bar, S_bond, C_bond)
                msg = str(cm.exception)
                assert channel in msg
                assert "(1, 0, 0)" in msg

    def test_dress_bond_accepts_a_well_conditioned_scalar_denominator(self):
        from hwave.solver.bond_channels import dress_bond
        for channel in ["spin", "charge"]:
            with self.subTest(channel=channel):
                chi_bar, S_bond, C_bond = _scalar_denominator_fixture(0.5, channel)
                chi_s, chi_c = dress_bond(chi_bar, S_bond, C_bond)
                assert np.all(np.isfinite(chi_s)) and np.all(np.isfinite(chi_c))

    def test_dress_bond_rejects_a_uniformly_small_denominator(self):
        """``eps * I`` at ND > 1 has sigma_min/sigma_max == 1 as well, yet the
        dressed vertices blow up by 1/eps; the guard must catch it."""
        from hwave.solver.bond_channels import dress_bond
        for channel in ["spin", "charge"]:
            with self.subTest(channel=channel):
                eps = 1.0e-9
                ND, Nx = 3, 2
                shape = (Nx, 1, 1, ND, ND)
                chi_bar = np.zeros(shape, dtype=complex)
                S_bond = np.zeros(shape, dtype=complex)
                C_bond = np.zeros(shape, dtype=complex)
                for ix in range(Nx):
                    chi_bar[ix, 0, 0] = np.eye(ND, dtype=complex)
                if channel == "spin":
                    S_bond[1, 0, 0] = (1.0 - eps) * np.eye(ND, dtype=complex)
                else:
                    C_bond[1, 0, 0] = -(1.0 - eps) * np.eye(ND, dtype=complex)
                with self.assertRaises(ValueError) as cm:
                    dress_bond(chi_bar, S_bond, C_bond)
                assert channel in str(cm.exception)

    def test_bond_within_tolerance_pair_is_projected_to_exact_hermiticity(self):
        """A declared forward/reverse bond pair that agrees only WITHIN the
        reversal tolerance must still come back EXACTLY Hermiticity-closed:
        coeffs[m] must equal conj(coeffs[reverse[m]].T) exactly."""
        v = 0.5 + 0.1j
        eps = 5.0e-9
        ci = {((1, 0, 0), (0, 1)): v,
              ((-1, 0, 0), (1, 0)): np.conj(v) + eps}
        topo = _bond_topo(ci, np.eye(3), norb=2)
        arr = np.asarray(topo.coeffs["CoulombInter"])
        rev = np.asarray(topo.reverse)
        for m in range(arr.shape[0]):
            a = arr[m]
            b = arr[int(rev[m])]
            assert np.array_equal(a, b.conj().T), "channel {} not closed".format(m)

    def test_max_shells_rejects_a_negative_value(self):
        """A negative max_shells must not be clamped to zero (that would
        silently drop every inter-site bond and bypass the zero-shell
        ambiguity guard)."""
        with self.assertRaisesRegex(ValueError, r"max_shells"):
            _bond_view(_nn_square(0.25), np.eye(3), norb=1, max_shells=-1)

    def test_max_shells_rejects_a_non_integral_value(self):
        with self.assertRaisesRegex(ValueError, r"max_shells"):
            _bond_view(_nn_square(0.25), np.eye(3), norb=1, max_shells=1.5)

    def test_max_shells_zero_rejects_a_tiny_declared_nonzero_v(self):
        """The documented rule is "any DECLARED nonzero inter-site V"; a tiny
        value is still declared and nonzero, so max_shells=0 stays ambiguous
        rather than silently discarding it."""
        with self.assertRaisesRegex(ValueError, r"max_shells=0"):
            _bond_view(_nn_square(1.0e-15), np.eye(3), norb=1, max_shells=0)


# --- R2-4 (optional): validate the documented green layout / beta ----------

class TestGreenLayoutAndBetaValidation(unittest.TestCase):

    def test_pair_weight_validates_green_and_beta(self):
        from hwave.solver.bond_channels import pair_weight
        for bad in ["ndim", "orbital", "beta"]:
            with self.subTest(bad=bad):
                green = np.ones((1, 1, 2, 1, 1, 2), dtype=complex)
                beta = 1.0
                if bad == "ndim":
                    green = np.ones((1, 1, 2, 1, 1), dtype=complex)
                elif bad == "orbital":
                    green = np.ones((1, 2, 2, 1, 1, 2), dtype=complex)
                else:
                    beta = -1.0
                with self.assertRaises(ValueError):
                    pair_weight(green, beta)


if __name__ == "__main__":
    unittest.main()
