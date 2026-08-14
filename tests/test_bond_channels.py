import unittest

import numpy as np

from hwave.solver.bond_channels import resolve_interactions, bare_bond_vertices
from tests.approx_util import ApproxTestCase

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


class TestResolveInteractionsAndBondBubble(ApproxTestCase):

    def test_bond_bubble_matches_direct_sum_tiny_grid(self):
        # 4x4x1 single-orbital tight-binding G; compare FFT bond bubble to
        # explicit k-sum of Eq.14 for a chosen (q, Δr, Δr').  (helper builds G.)
        from hwave.solver.bond_channels import bond_bubble, resolve_interactions
        green, T, N = _toy_green_4x4()                      # helper in test file
        beta = 1.0 / T
        bset = resolve_interactions(_nn_square(0.25), np.eye(3), norb=1)
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
        from hwave.solver.bond_channels import bond_bubble, resolve_interactions
        green, T, N = _toy_green_4x4_2orb(Nx=4, Ny=4, Nz=1, Nmat=4)
        beta = 1.0 / T
        bset = resolve_interactions(_nn_square_2orb(), np.eye(3), norb=2)
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
        r = resolve_interactions(_nn_square(), np.eye(3), norb=1)
        assert r.n_channels == 5           # {0,+x,-x,+y,-y}
        assert r.delta_r[0] == (0,0,0)     # Δr=0 first
        # reverse map pairs +R with -R
        for m, dr in enumerate(r.delta_r):
            rm = r.delta_r[r.reverse[m]]
            assert tuple(-np.array(dr)) == rm

    def test_reverse_synthesis_from_plus_R_only(self):
        ci = {((1,0,0),(0,0)): 0.3}                  # +x only
        r = resolve_interactions(ci, np.eye(3), norb=1)
        assert r.n_channels == 3                      # {0,+x,-x} synthesized
        ix = r.delta_r.index((1,0,0)); imx = r.delta_r.index((-1,0,0))
        self.assertApprox(r.v_bond[ix][0,0], 0.3)
        self.assertApprox(r.v_bond[imx][0,0], 0.3)  # V_ab(R)=V_ba(-R)* real

    def test_v0_declared_topology_kept(self):
        ci = {((1,0,0),(0,0)): 0.0, ((0,1,0),(0,0)): 0.0}  # declared, zero-valued
        r = resolve_interactions(ci, np.eye(3), norb=1)
        assert r.n_channels == 5           # topology declared -> B stays 5 at V=0

    def test_cutoff_shell_zero_with_nonzero_V_errors(self):
        with self.assertRaises(ValueError):
            resolve_interactions(_nn_square(0.25), np.eye(3), norb=1, bond_max_shells=0)

    def test_bond_reverse_mismatch_raises_for_nonzero_bond(self):
        """The m != 0 (genuine bond, Delta r != 0) counterpart of
        test_v_onsite_inconsistent_orientations_raises: when BOTH V_01(+R) and
        V_10(-R) are declared but disagree with the Hermiticity constraint
        V_ab(R) = conj(V_ba(-R)) beyond tolerance, resolve_interactions must raise
        -- not silently pick one side or average them. This is the m != 0 variant
        of the consolidated consistency check (bond_channels.py resolve_interactions
        Step 2); only the on-site (R=0) mismatch was previously covered."""
        norb = 2
        # V_01(+R) = 0.5+0.1j : declared. Consistency requires V_10(-R) =
        # conj(V_01(+R)) = 0.5-0.1j, but 0.1+0.5j is declared instead --
        # disagreement is not a rounding-level effect (|diff| = 0.72).
        ci = {((1, 0, 0), (0, 1)): 0.5 + 0.1j, ((-1, 0, 0), (1, 0)): 0.1 + 0.5j}
        with self.assertRaisesRegex(ValueError, r"resolve_interactions.*disagrees"):
            resolve_interactions(ci, np.eye(3), norb=norb)


# ---------------------------------------------------------------------------
# I-1 (review fix): the on-site (Delta r=0) block must be Hermiticity-closed
# by resolve_interactions itself -- the SAME reverse/consistency synthesis
# applied to the bond (R!=0) channels -- via a new ``v_onsite`` field, instead
# of sc.py reading the raw dict directly (which bypassed the closure and could
# silently produce a non-Hermitian on-site Fock block for a multi-orbital
# model; spec S3.0).
# ---------------------------------------------------------------------------

class TestOnsiteHermiticityClosure(ApproxTestCase):

    def test_v_onsite_hermitian_closed_when_only_one_orientation_declared(self):
        norb = 2
        v = 0.7 + 0.2j
        ci = {((0, 0, 0), (0, 1)): v}          # only V_01(0) declared, not V_10(0)
        r = resolve_interactions(ci, np.eye(3), norb=norb)
        assert r.v_onsite.shape == (norb, norb)
        self.assertApprox(r.v_onsite[0, 1], v)
        self.assertApprox(r.v_onsite[1, 0], np.conj(v))   # synthesized

    def test_v_onsite_inconsistent_orientations_raises(self):
        norb = 2
        # V_01(0) and V_10(0) declared but NOT conjugate-consistent
        ci = {((0, 0, 0), (0, 1)): 0.5 + 0.1j, ((0, 0, 0), (1, 0)): 0.1 + 0.5j}
        with self.assertRaises(ValueError):
            resolve_interactions(ci, np.eye(3), norb=norb)

    def test_v_onsite_consistent_both_orientations_declared(self):
        norb = 2
        v = 0.4 - 0.05j
        ci = {((0, 0, 0), (0, 1)): v, ((0, 0, 0), (1, 0)): np.conj(v)}
        r = resolve_interactions(ci, np.eye(3), norb=norb)
        self.assertApprox(r.v_onsite[0, 1], v)
        self.assertApprox(r.v_onsite[1, 0], np.conj(v))

    def test_v_onsite_all_zero_when_undeclared(self):
        norb = 2
        r = resolve_interactions(_nn_square_2orb(), np.eye(3), norb=norb)
        assert np.allclose(r.v_onsite, 0.0)


# ---------------------------------------------------------------------------
# Task 3: bare_bond_vertices -- S_bond, C_bond, and the particle-particle
# (Cooper) vertices Vpp_s, Vpp_t. Analytic element fixtures (spec 7.4/7.9/7.11).
# ---------------------------------------------------------------------------

def _idx(norb, m, l1, l2):
    """Bond-major enlarged index I = m*nd + (l1*norb + l2)."""
    return m * (norb * norb) + l1 * norb + l2


def _build_two_orbital_case():
    """Spec 7.4 fixture. 2 orbitals, channels {0,+R,-R}, on-site U'=V01[R=0]=1,
    inter-site V01[R=+-R]=1/4. Intra U=0, no Hund/Exchange/PairHop.

    Returns (bond_set, S0_q, C0_q, norb) with S0_q/C0_q the Case-2-corrected
    m=0 ph blocks on a 1x1x1 grid (only R=0 Fock in the (ab,ab) element, the
    full q=0 Hartree 2*V01(q=0) in the (aa,bb) charge element).
    """
    norb = 2
    nd = norb * norb
    ci = {((1, 0, 0), (0, 1)): 0.25, ((1, 0, 0), (1, 0)): 0.25,
          ((-1, 0, 0), (0, 1)): 0.25, ((-1, 0, 0), (1, 0)): 0.25}
    bond_set = resolve_interactions(ci, np.eye(3), norb=norb)
    assert bond_set.n_channels == 3

    Up = 1.0                       # on-site U' = V01(R=0)
    Vq0 = Up + 0.25 + 0.25         # V01(q=0) = R=0 + (+R) + (-R) = 1.5
    S0 = np.zeros((1, 1, 1, nd, nd), dtype=complex)
    C0 = np.zeros((1, 1, 1, nd, nd), dtype=complex)

    def i(a, b):
        return a * norb + b

    # Fock (ab,ab): only on-site U' survives (Case-2-corrected)
    S0[0, 0, 0, i(0, 1), i(0, 1)] = +Up
    S0[0, 0, 0, i(1, 0), i(1, 0)] = +Up
    C0[0, 0, 0, i(0, 1), i(0, 1)] = -Up
    C0[0, 0, 0, i(1, 0), i(1, 0)] = -Up
    # Hartree (aa,bb): full q=0 Fourier sum 2*V01(q=0)
    C0[0, 0, 0, i(0, 0), i(1, 1)] = 2.0 * Vq0
    C0[0, 0, 0, i(1, 1), i(0, 0)] = 2.0 * Vq0
    return bond_set, S0, C0, norb


def _build_asymmetric_two_orbital_bond(VAB=0.4 + 0.3j, VBA=0.1 - 0.2j):
    """A genuinely asymmetric two-orbital bond, V_01(+R) = VAB != VBA =
    V_10(+R) (both orientations declared directly, so no synthesis is
    involved), with no on-site (R=0) CoulombInter declared.

    The two RAW "totals" (the R != 0 sum entering the Hartree correction)
    ``total_01 = V_01(+R) + V_01(-R)`` and ``total_10 = V_10(+R) + V_10(-R)``
    satisfy ``conj(total_10) == total_01`` (a structural consequence of
    reversal-Hermiticity, true for ANY reversal-closed set) but are
    themselves DIFFERENT complex numbers whenever either value has a nonzero
    imaginary part.

    UPDATE (Finding 1 fix, #151 ED-adjudication campaign): the Hartree sum
    entering ``sc._build_bond_m0_blocks``/``bare_bond_vertices`` no longer
    uses these RAW totals -- it Hermitizes each shell PER SLOT
    (``Re(V_ab(R_m))``) before summing, per the operator-algebra derivation
    at the Hermitization site (a declared channel and its
    ``resolve_interactions`` reversal partner multiply the SAME physical
    density operator, so Hermiticity forces their combined coefficient to be
    real). The HERMITIZED totals, ``Re(V_01(+R))+Re(V_01(-R))`` and
    ``Re(V_10(+R))+Re(V_10(-R))``, turn out to be EQUAL to each other for
    this (or any 2-orbital, single-shell) fixture -- both reduce to
    ``Re(VAB)+Re(VBA)`` -- which is why they no longer distinguish a
    Hartree-subtraction orbital-order transposition bug (the whole class of
    bug the raw totals' inequality used to catch is now structurally
    impossible: post-fix, this Hartree block is always Hermitian/real by
    construction). ``VAB``/``VBA`` are still asymmetric on their REAL PARTS
    (``Re(VAB)=0.4 != Re(VBA)=0.1``), so the fixture still exercises and
    pins the Fock-placement orbital order (S_bond/C_bond bond-diagonal
    element).
    """
    norb = 2
    ci = {
        ((1, 0, 0), (0, 1)): VAB,
        ((1, 0, 0), (1, 0)): VBA,
        ((-1, 0, 0), (0, 1)): np.conj(VBA),
        ((-1, 0, 0), (1, 0)): np.conj(VAB),
    }
    bond_set = resolve_interactions(ci, np.eye(3), norb=norb)
    assert bond_set.n_channels == 3
    return bond_set, norb, VAB, VBA


def _build_local_full_two_orbital(U=3.0, Up=1.0, J=0.5, Jp=0.2, PH=0.1, Vbond=0.25):
    """2-orbital local fixture with U,U',J,J',PH and one inter-site bond
    V01[+-R]=Vbond, for spec 7.11 (complete local Cooper vertex + B=1 crossing).
    """
    norb = 2
    nd = norb * norb
    ci = {((1, 0, 0), (0, 1)): Vbond, ((1, 0, 0), (1, 0)): Vbond,
          ((-1, 0, 0), (0, 1)): Vbond, ((-1, 0, 0), (1, 0)): Vbond}
    bond_set = resolve_interactions(ci, np.eye(3), norb=norb)

    Vq0 = Up + 2.0 * Vbond   # V01(q=0) = R=0 + (+R)+(-R)
    S0 = np.zeros((1, 1, 1, nd, nd), dtype=complex)
    C0 = np.zeros((1, 1, 1, nd, nd), dtype=complex)

    def i(a, b):
        return a * norb + b

    for a in range(norb):
        # Case 1 (aa,aa): S=C=U
        S0[0, 0, 0, i(a, a), i(a, a)] = U
        C0[0, 0, 0, i(a, a), i(a, a)] = U
    for a in range(norb):
        for b in range(norb):
            if a == b:
                continue
            # Case 2 (ab,ab): S=U', C=-U'+J
            S0[0, 0, 0, i(a, b), i(a, b)] = Up
            C0[0, 0, 0, i(a, b), i(a, b)] = -Up + J
            # Case 3 (aa,bb): S=J, C=2*V(q)-J
            S0[0, 0, 0, i(a, a), i(b, b)] = J
            C0[0, 0, 0, i(a, a), i(b, b)] = 2.0 * Vq0 - J
            # Case 4 (ab,ba): S=C=J'+PH
            S0[0, 0, 0, i(a, b), i(b, a)] = Jp + PH
            C0[0, 0, 0, i(a, b), i(b, a)] = Jp + PH
    return bond_set, S0, C0, norb, dict(U=U, Up=Up, J=J, Jp=Jp, PH=PH)


class TestBareBondVertices(ApproxTestCase):

    def test_bare_matrix_elements_two_orbital(self):
        """Spec 7.4: Fock signs & placement in S_bond/C_bond."""
        bond_set, S0, C0, norb = _build_two_orbital_case()
        S_bond, C_bond, _, _ = bare_bond_vertices(bond_set, S0, C0, norb)

        q0 = (0, 0, 0)
        # m=0 Fock (ab,ab): +U' in S, -U' in C
        self.assertApprox(
            S_bond[q0][_idx(norb, 0, 0, 1), _idx(norb, 0, 0, 1)], 1.0)
        self.assertApprox(
            C_bond[q0][_idx(norb, 0, 0, 1), _idx(norb, 0, 0, 1)], -1.0)
        # m=0 Hartree (aa,bb): 2*V01(q=0) = 3, only in the charge C_bond
        self.assertApprox(
            C_bond[q0][_idx(norb, 0, 0, 0), _idx(norb, 0, 1, 1)], 3.0)
        self.assertApprox(
            S_bond[q0][_idx(norb, 0, 0, 0), _idx(norb, 0, 1, 1)], 0.0)

        # m != 0 bond-diagonal Fock: +V01(R_m) in S, -V01(R_m) in C, for both
        # +R and -R channels (m=1 is -R, m=2 is +R for this bond_set).
        for m in (1, 2):
            self.assertApprox(
                S_bond[q0][_idx(norb, m, 0, 1), _idx(norb, m, 0, 1)], 0.25)
            self.assertApprox(
                C_bond[q0][_idx(norb, m, 0, 1), _idx(norb, m, 0, 1)], -0.25)
        # no (R,01) <-> (R,10) bare-vertex element (Fock is diagonal in the pair idx)
        self.assertApprox(
            S_bond[q0][_idx(norb, 1, 0, 1), _idx(norb, 1, 1, 0)], 0.0)
        self.assertApprox(
            C_bond[q0][_idx(norb, 1, 0, 1), _idx(norb, 1, 1, 0)], 0.0)

        # S_bond/C_bond are block-diagonal in the bond index m (m=0 and m!=0 blocks
        # never couple; distinct m!=0 channels never couple).
        nd = norb * norb
        for m in range(1, bond_set.n_channels):
            # off-diagonal-in-m coupling to the m=0 block is zero
            blk = S_bond[q0][m * nd:(m + 1) * nd, 0:nd]
            assert np.allclose(blk, 0.0)

    def test_bare_bond_vertices_asymmetric_V_orbital_order_pinned(self):
        """Task 3/6 (review); updated by the Finding 1 fix (#151 ED-adjudication
        campaign). Pins that bare_bond_vertices' Hartree subtraction (the Case-2
        correction inside the local Cooper block) and its m != 0 Fock placement
        (S_bond/C_bond bond-diagonal element) use the SAME orbital-pair index
        order -- a transposition in either would change the numeric result here
        (unlike the file's other bare_bond_vertices fixtures, which all declare
        a symmetric V_01 == V_10 and so cannot detect this).

        Both slots now Hermitize per shell before use (Re(V_ab(R_m)) -- see
        bare_bond_vertices' and sc._build_bond_m0_blocks' inline derivations):
        a declared channel and its resolve_interactions reversal partner
        multiply the SAME physical density operator after site relabeling, so
        Hermiticity forces their combined coefficient to be real. VAB/VBA's
        REAL PARTS (0.4 vs 0.1) are still different, so the Fock-placement
        check below still catches an orbital-order transposition directly (no
        cancellation needed). The Hartree totals, however, become EQUAL once
        Hermitized (see _build_asymmetric_two_orbital_bond's docstring) --
        that specific transposition class is now structurally impossible, so
        this test's Hartree section (2) below only confirms this fixture's
        production-matching cancellation still holds, not an orbital-order
        distinction.
        """
        bond_set, norb, VAB, VBA = _build_asymmetric_two_orbital_bond()
        nd = norb * norb

        def i(a, b):
            return a * norb + b

        # Hermitized totals (Finding 1 fix): Re(V_ab(+R)) + Re(V_ab(-R)),
        # computed directly from the DECLARED values (independent of
        # bond_set.v_bond, so this is not circular with the implementation).
        total_01 = VAB.real + np.conj(VBA).real   # Re(V_01(+R)) + Re(V_01(-R))
        total_10 = VBA.real + np.conj(VAB).real   # Re(V_10(+R)) + Re(V_10(-R))
        # Re(conj(z)) == Re(z), so both totals reduce to Re(VAB)+Re(VBA) and
        # are therefore equal -- a structural consequence of Hermitization,
        # not a fixture coincidence (see the docstring above / the builder's).
        self.assertApprox(total_01, total_10)
        self.assertApprox(total_01, VAB.real + VBA.real)
        assert abs(VAB.real - VBA.real) > 1e-6   # Fock check below is not a no-op

        # Case-2-corrected q=0 Hartree the caller (sc._build_bond_m0_blocks)
        # would hand to bare_bond_vertices: C0[aa,bb] = 2 * V_ab(q=0), with
        # V_ab(q=0) = V_ab(R=0) [undeclared -> 0] + total_ab (Hermitized).
        S0 = np.zeros((1, 1, 1, nd, nd), dtype=complex)
        C0 = np.zeros((1, 1, 1, nd, nd), dtype=complex)
        C0[0, 0, 0, i(0, 0), i(1, 1)] = 2.0 * total_01
        C0[0, 0, 0, i(1, 1), i(0, 0)] = 2.0 * total_10

        S_bond, C_bond, Vpp_s, Vpp_t = bare_bond_vertices(bond_set, S0, C0, norb)
        q0 = (0, 0, 0)

        # (1) Fock placement: the m != 0 bond-diagonal element (ab,ab) must carry
        # +/-Re(V_ab(R_m)), NOT +/-Re(V_ba(R_m)) -- Re(VAB) != Re(VBA) makes a
        # transposition visible directly (no cancellation needed).
        m_plus = bond_set.delta_r.index((1, 0, 0))
        self.assertApprox(
            S_bond[q0][_idx(norb, m_plus, 0, 1), _idx(norb, m_plus, 0, 1)],
            VAB.real)
        self.assertApprox(
            S_bond[q0][_idx(norb, m_plus, 1, 0), _idx(norb, m_plus, 1, 0)],
            VBA.real)
        self.assertApprox(
            C_bond[q0][_idx(norb, m_plus, 0, 1), _idx(norb, m_plus, 0, 1)],
            -VAB.real)
        self.assertApprox(
            C_bond[q0][_idx(norb, m_plus, 1, 0), _idx(norb, m_plus, 1, 0)],
            -VBA.real)

        # (2) Hartree subtraction (Case-2 correction): with no other local
        # interaction declared, the fully-corrected local block (read off
        # Vpp_s/Vpp_t's m=0 sub-block, which is built from S0_loc/C0_loc AFTER
        # the correction) must cancel EXACTLY to zero on BOTH the (0,0)<-(1,1)
        # and (1,1)<-(0,0) charge elements -- production's own hartree_rneq0
        # subtraction Hermitizes the identical way, so it matches this
        # fixture's C0 exactly regardless of which of the (now-equal) totals
        # lands in which slot.
        np.testing.assert_allclose(Vpp_s[0:nd, 0:nd], 0.0, atol=1e-12)
        np.testing.assert_allclose(Vpp_t[0:nd, 0:nd], 0.0, atol=1e-12)

    def test_pp_vertex_singleband_fixture(self):
        """Spec 7.9: 1D single band, U=4, V=1, chi=0. Vpp on basis {0,+R,-R}."""
        norb = 1
        ci = {((1, 0, 0), (0, 0)): 1.0, ((-1, 0, 0), (0, 0)): 1.0}
        bond_set = resolve_interactions(ci, np.eye(3), norb=norb)
        assert bond_set.n_channels == 3

        U, V = 4.0, 1.0
        # m=0 ph block: S=U, C=U+2*V(q).  V(q=0)=V(+R)+V(-R)=2 -> Hartree 2*V(q=0)=4
        S0 = np.full((1, 1, 1, 1, 1), U, dtype=complex)
        C0 = np.full((1, 1, 1, 1, 1), U + 2.0 * (2.0 * V) / 2.0 * 2.0, dtype=complex)
        # (be explicit: C0 = U + 2*V(q=0) = 4 + 2*2 = 8)
        C0 = np.full((1, 1, 1, 1, 1), 8.0, dtype=complex)

        _, _, Vpp_s, Vpp_t = bare_bond_vertices(bond_set, S0, C0, norb)

        # Hartree excluded: (0,0) = 2U = 8 (NOT 2U+2V); bond block D(1+-P).
        np.testing.assert_allclose(
            Vpp_s, [[8, 0, 0], [0, 1, 1], [0, 1, 1]], atol=1e-12)
        np.testing.assert_allclose(
            Vpp_t, [[0, 0, 0], [0, 1, -1], [0, -1, 1]], atol=1e-12)

    def test_pp_local_L_table_and_parity_subspace(self):
        """Spec 7.11: complete local L^{s/t} table entries, the m=0/m!=0 split
        (no double counting), and parity projection on the Q_eta subspace."""
        bond_set, S0, C0, norb, p = _build_local_full_two_orbital()
        nd = norb * norb
        _, _, Vpp_s, Vpp_t = bare_bond_vertices(bond_set, S0, C0, norb)

        L_s = Vpp_s[0:nd, 0:nd]
        L_t = Vpp_t[0:nd, 0:nd]

        def i(a, b):
            return a * norb + b

        U, Up, J, Jp, PH = p["U"], p["Up"], p["J"], p["Jp"], p["PH"]
        # L^s table
        self.assertApprox(L_s[i(0, 0), i(0, 0)], 2 * U)                  # aa<-aa
        self.assertApprox(L_s[i(0, 1), i(0, 1)], 2 * Up)                 # ab<-ab
        self.assertApprox(L_s[i(0, 1), i(1, 0)], J)                      # ab<-ba
        self.assertApprox(L_s[i(0, 0), i(1, 1)], 2 * (Jp + PH))          # aa<-bb
        # L^t table
        self.assertApprox(L_t[i(0, 0), i(0, 0)], 0.0)                    # aa<-aa
        self.assertApprox(L_t[i(0, 1), i(0, 1)], 2 * (Up - J))           # ab<-ab
        self.assertApprox(L_t[i(0, 1), i(1, 0)], -2 * Up + J)            # ab<-ba
        self.assertApprox(L_t[i(0, 0), i(1, 1)], 0.0)                    # aa<-bb

        # No double counting: the inter-site bond V01(R) does NOT leak into the
        # local (m=0) block -- L^s(01,01) stays exactly 2U' (on-site only).
        self.assertApprox(L_s[i(0, 1), i(0, 1)], 2 * Up)

        # Parity: build the local pair-reversal P_loc:(cd)->(dc), Q_eta=(I+-P)/2.
        P_loc = np.zeros((nd, nd), dtype=complex)
        for a in range(norb):
            for b in range(norb):
                P_loc[i(b, a), i(a, b)] = 1.0
        Q_s = 0.5 * (np.eye(nd) + P_loc)
        Q_t = 0.5 * (np.eye(nd) - P_loc)

        Vs_proj = Q_s @ L_s @ Q_s
        Vt_proj = Q_t @ L_t @ Q_t
        # projects onto its parity sector: P V = eta V
        np.testing.assert_allclose(P_loc @ Vs_proj, +1.0 * Vs_proj, atol=1e-12)
        np.testing.assert_allclose(P_loc @ Vt_proj, -1.0 * Vt_proj, atol=1e-12)
        # opposite-parity sector annihilated
        np.testing.assert_allclose(Vs_proj @ Q_t, 0.0, atol=1e-12)
        np.testing.assert_allclose(Vt_proj @ Q_s, 0.0, atol=1e-12)

        # Reduction to the existing sc.py instantaneous operator: the crossed L is
        # exactly the existing Vexist[a,b,c,d]=(S+C)_{(ab),(cd)} reinterpreted as a
        # pp-pair operator E[(a,d),(b,c)] -- i.e. L == crossing(S+C) with Hartree
        # excluded. Verify on the parity subspace they coincide.
        W_s = (S0[0, 0, 0] + C0[0, 0, 0]).copy()
        # remove the inter-site (R!=0) Hartree from (aa,bb) to get the local W_s
        for a in range(norb):
            for b in range(norb):
                W_s[i(a, a), i(b, b)] -= 2.0 * sum(
                    bond_set.v_bond[m][a, b] for m in range(1, bond_set.n_channels))
        Vexist4 = W_s.reshape(norb, norb, norb, norb)
        E_s = Vexist4.transpose(0, 3, 1, 2).reshape(nd, nd)
        np.testing.assert_allclose(L_s, E_s, atol=1e-12)


# ---------------------------------------------------------------------------
# Task 4: dress_bond -- enlarged (ND x ND) RPA-ladder dressing, plus the
# complex-interaction Hermiticity/parity oracle for bare_bond_vertices' Vpp
# (spec S4.4 / S4.5 / S7.12).
# ---------------------------------------------------------------------------

def _build_complex_two_orbital_bond():
    """2-orbital bond fixture with a genuinely complex, non-inversion-
    -symmetric-real value: V_01(+R) = 0.3 + 0.1j (only this entry declared;
    resolve_interactions synthesizes V_10(-R) = conj(V_01(+R)) automatically
    per the reversal-closure Hermiticity rule V_ab(R)=V_ba(-R)*).

    ``bare_bond_vertices`` assumes ``C0_q`` already carries the FULL q=0
    Hartree ``2*V_ab(q=0)`` (its Case-2-correction step subtracts the R!=0
    part of it to recover the local-only piece -- see its docstring/S4.3
    star); passing a bare zero C0 here would violate that precondition and
    make the subtraction inject spurious values into the local block. So we
    set the (aa,bb) Hartree entries to exactly ``2*V_ab(q=0)`` (with no
    on-site V declared, ``V_ab(q=0)`` is just the R!=0 bond sum), which the
    Case-2 correction then cancels back to zero -- leaving the local (m=0)
    block identically zero. Vpp_s/Vpp_t are then supported purely on the
    bond (m!=0) sector, isolating the Q_eta(D+D^dag)Q_eta term for the
    oracle.

    ``hartree_rneq0`` Hermitizes each shell (``.real``) before summing,
    matching ``sc._build_bond_m0_blocks``' Finding 1 fix (#151 ED-
    adjudication campaign) -- see that function's inline comment for the
    derivation. Without this, this fixture's caller-supplied C0 would no
    longer match what ``bare_bond_vertices``' own (also-fixed)
    ``hartree_rneq0`` subtracts, reintroducing a spurious residual into the
    local block this fixture is built to keep exactly zero.
    """
    norb = 2
    nd = norb * norb
    v = 0.3 + 0.1j
    ci = {((1, 0, 0), (0, 1)): v}
    bond_set = resolve_interactions(ci, np.eye(3), norb=norb)
    assert bond_set.n_channels == 3
    S0 = np.zeros((1, 1, 1, nd, nd), dtype=complex)
    C0 = np.zeros((1, 1, 1, nd, nd), dtype=complex)

    def i(a, b):
        return a * norb + b

    for a in range(norb):
        for b in range(norb):
            hartree_rneq0 = sum(
                bond_set.v_bond[m][a, b].real
                for m in range(1, bond_set.n_channels)
            )
            C0[0, 0, 0, i(a, a), i(b, b)] = 2.0 * hartree_rneq0

    return bond_set, S0, C0, norb, v


def _explicit_D_and_P(bond_set, norb):
    """Independent (test-side) construction of D_off and the parity
    permutation P, from the bond_set's public delta_r/v_bond/reverse -- built
    directly from the spec formula (S4.5), not by importing the
    implementation's internals."""
    nd = norb * norb
    B = bond_set.n_channels
    ND = nd * B
    D = np.zeros((ND, ND), dtype=complex)
    for m in range(1, B):
        vb = bond_set.v_bond[m]
        for l1 in range(norb):
            for l2 in range(norb):
                I = m * nd + l1 * norb + l2
                D[I, I] = vb[l1, l2]
    P = np.zeros((ND, ND), dtype=complex)
    for m in range(B):
        mr = bond_set.reverse[m]
        for l1 in range(norb):
            for l2 in range(norb):
                J = m * nd + l1 * norb + l2
                I = mr * nd + l2 * norb + l1
                P[I, J] = 1.0
    return D, P


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

    def test_complex_bond_pp_vertex_hermitian_and_parity_projected(self):
        """Spec S4.5 (complex paragraph) / S7.12: for a genuinely complex bond
        V_ab(R)=V_ba(-R)*, the derived Cooper vertex V^pp_eta = Q_eta(D+D^dag)Q_eta
        must be (a) Hermitian, (b) parity-projected: P V^pp_eta = eta V^pp_eta,
        and (c) exactly equal to the explicit Q_eta(D+D^dag)Q_eta built
        independently from the bond_set's public data (pins the formula, not just
        the invariants)."""
        bond_set, S0, C0, norb, v = _build_complex_two_orbital_bond()
        _, _, Vpp_s, Vpp_t = bare_bond_vertices(bond_set, S0, C0, norb)
        ND = Vpp_s.shape[0]

        # (a) Hermitian
        np.testing.assert_allclose(Vpp_s, Vpp_s.conj().T, atol=1e-12)
        np.testing.assert_allclose(Vpp_t, Vpp_t.conj().T, atol=1e-12)

        # (b) parity-projected: P @ Vpp_eta = eta * Vpp_eta
        D, P = _explicit_D_and_P(bond_set, norb)
        np.testing.assert_allclose(P @ Vpp_s, +1.0 * Vpp_s, atol=1e-12)
        np.testing.assert_allclose(P @ Vpp_t, -1.0 * Vpp_t, atol=1e-12)

        # (c) exact formula match: Q_eta (D + D^dag) Q_eta, built independently
        # from bond_set's public delta_r/v_bond/reverse (local block is zero here,
        # so this pins the bond sector exactly).
        Id = np.eye(ND, dtype=complex)
        Dh = D + D.conj().T
        Q_s = 0.5 * (Id + P)
        Q_t = 0.5 * (Id - P)
        np.testing.assert_allclose(Vpp_s, Q_s @ Dh @ Q_s, atol=1e-12)
        np.testing.assert_allclose(Vpp_t, Q_t @ Dh @ Q_t, atol=1e-12)

        # sanity: D itself is NOT Hermitian (this is a genuinely complex bond;
        # otherwise the Hermiticity assertion above would be a vacuous check that
        # even a buggy "V^pp = D(1+-P)" naive form could pass by accident here).
        assert not np.allclose(D, D.conj().T)

    def test_real_inversion_symmetric_bond_reduces_to_D_times_one_plus_minus_P(self):
        """Spec S4.5: for a REAL inversion-symmetric bond ([D,P]=0), the derived
        Q_eta(D+D^dag)Q_eta form collapses to the simpler D_off(1 +- P). Cross-
        checks the complex-oracle machinery against the existing real single-band
        fixture (test_pp_vertex_singleband_fixture, spec S7.9)."""
        norb = 1
        ci = {((1, 0, 0), (0, 0)): 1.0, ((-1, 0, 0), (0, 0)): 1.0}
        bond_set = resolve_interactions(ci, np.eye(3), norb=norb)
        U, V = 4.0, 1.0
        S0 = np.full((1, 1, 1, 1, 1), U, dtype=complex)
        C0 = np.full((1, 1, 1, 1, 1), 8.0, dtype=complex)  # U + 2*V(q=0) = 4 + 4

        _, _, Vpp_s, Vpp_t = bare_bond_vertices(bond_set, S0, C0, norb)
        ND = Vpp_s.shape[0]

        D, P = _explicit_D_and_P(bond_set, norb)
        # real inversion-symmetric bond: D is real and commutes with P
        assert np.allclose(D.imag, 0.0)
        assert np.allclose(D @ P, P @ D, atol=1e-12)

        # The reduction D_off(1+-P) applies to the bond (m!=0) sector only -- the
        # m=0 local block (here L_s(0,0)=2U=8, from S0=C0=U/8 at norb=1) is a
        # separate term (spec S4.5 "no double counting"), not part of D_off.
        nd = norb * norb
        Id = np.eye(ND, dtype=complex)
        D1P = D @ (Id + P)
        D1mP = D @ (Id - P)
        np.testing.assert_allclose(Vpp_s[nd:, nd:], D1P[nd:, nd:], atol=1e-12)
        np.testing.assert_allclose(Vpp_t[nd:, nd:], D1mP[nd:, nd:], atol=1e-12)
        # local (m=0) block untouched by the bond reduction
        np.testing.assert_allclose(Vpp_s[0:nd, 0:nd], [[8.0]], atol=1e-12)
        np.testing.assert_allclose(Vpp_t[0:nd, 0:nd], [[0.0]], atol=1e-12)

    def test_bare_bond_vertices_rejects_non_closed_set(self):
        """Internal-interface guard: a non-reversal-closed bond_set raises."""
        from hwave.solver.bond_channels import ResolvedInteractionSet
        norb = 1
        nd = 1
        # a hand-built set whose reverse map is NOT an involution / not closed
        bad = ResolvedInteractionSet(
            delta_r=((0, 0, 0), (1, 0, 0)),
            v_bond=(np.zeros((1, 1), dtype=complex), np.ones((1, 1), dtype=complex)),
            reverse=(0, 1),                       # (1,0,0) claims itself as its reverse
            n_channels=2,
        )
        S0 = np.zeros((1, 1, 1, nd, nd), dtype=complex)
        C0 = np.zeros((1, 1, 1, nd, nd), dtype=complex)
        with self.assertRaises(ValueError):
            bare_bond_vertices(bad, S0, C0, norb)


# ---------------------------------------------------------------------------
# Task 5: make_bond_kernel -- two-part Eliashberg operator (fluctuation +
# instantaneous / bare-pp), spec S4.5.  Reference-kernel oracle (S7.5) and the
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


def _direct_bond_kernel(F_q, Vpp, green, bond_set, beta, phi):
    """Literal O(N_k^2) evaluation of the two-part bond Eliashberg operator
    (spec S4.5), independent of make_bond_kernel.

    K phi(k)_{ab} = -(1/N) sum_{k'} [ Veff_fl(k,k') + Veff_pp(k,k') ]_{ab,cd}
                                     X_cd(k')
    with
      X_cd(k')      = sum_{ef} G2[c,e,d,f](k') phi_{ef}(k')      (cd = c*norb+d)
      Veff_fl_{ab,cd}(k,k') = sum_{m,m'} e^{i k.dr_m} F_q(k-k')_{(m,ab),(m',cd)}
                                          e^{i k'.dr_m'}
      Veff_pp_{ab,cd}(k,k') = 1/2 sum_{m,m'} Vpp_{(m,ab),(m',cd)}
                                          e^{i (k.dr_m' - k'.dr_m)}
    """
    norb = green.shape[0]
    nd = norb * norb
    B = bond_set.n_channels
    Nx, Ny, Nz, nmat = green.shape[2:6]
    N = Nx * Ny * Nz
    dr = bond_set.delta_r

    G2 = _g2_direct(green, beta)
    # X_cd(k) with cd = c*norb+d
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
                                      * phi[e, f, kx, ky, kz])
                        X[c * norb + d, kx, ky, kz] = s

    kx_list = 2.0 * np.pi * np.arange(Nx) / Nx
    ky_list = 2.0 * np.pi * np.arange(Ny) / Ny
    kz_list = 2.0 * np.pi * np.arange(Nz) / Nz

    def kdot(ix, iy, iz, d):
        return kx_list[ix] * d[0] + ky_list[iy] * d[1] + kz_list[iz] * d[2]

    out = np.zeros((nd, Nx, Ny, Nz), dtype=complex)
    for kx in range(Nx):
        for ky in range(Ny):
            for kz in range(Nz):
                acc = np.zeros(nd, dtype=complex)
                for px in range(Nx):
                    for py in range(Ny):
                        for pz in range(Nz):
                            qx, qy, qz = (kx - px) % Nx, (ky - py) % Ny, (kz - pz) % Nz
                            Fblk = F_q[qx, qy, qz]
                            Veff = np.zeros((nd, nd), dtype=complex)
                            for m in range(B):
                                em = np.exp(1j * kdot(kx, ky, kz, dr[m]))
                                for mp in range(B):
                                    emp = np.exp(1j * kdot(px, py, pz, dr[mp]))
                                    Veff += em * emp * Fblk[
                                        m * nd:(m + 1) * nd,
                                        mp * nd:(mp + 1) * nd]
                            Vpp_eff = np.zeros((nd, nd), dtype=complex)
                            for m in range(B):
                                for mp in range(B):
                                    ph = np.exp(1j * (
                                        kdot(kx, ky, kz, dr[mp])
                                        - kdot(px, py, pz, dr[m])))
                                    Vpp_eff += 0.5 * ph * Vpp[
                                        m * nd:(m + 1) * nd,
                                        mp * nd:(mp + 1) * nd]
                            acc += (Veff + Vpp_eff) @ X[:, px, py, pz]
                out[:, kx, ky, kz] = -(1.0 / N) * acc
    return out.reshape(norb, norb, Nx, Ny, Nz)


def _fluctuation_vertex(pairing_type, S_bond, chi_s, C_bond, chi_c):
    """F_eta(q) per spec S4.5 (batched matmul over q)."""
    SxS = S_bond @ chi_s @ S_bond
    CxC = C_bond @ chi_c @ C_bond
    if pairing_type == "singlet":
        return 1.5 * SxS - 0.5 * CxC
    return -0.5 * SxS - 0.5 * CxC


class TestMakeBondKernel(unittest.TestCase):

    def test_bond_kernel_reference_oracle_singleband(self):
        """Spec S7.5: the LinearOperator equals a direct O(N_k^2) evaluation of the
        two-part V_eff(k,k') on a 4x4 grid, pinning the two phase conventions
        (fluctuation e^{i(k.dr_m + k'.dr_m')}, instantaneous cross phase
        e^{i(k.dr_m' - k'.dr_m)}), the q=k-k' argument, and the -(1/N) / -1/2
        normalizations."""
        from hwave.solver.bond_channels import make_bond_kernel

        for pairing_type in ["singlet", "triplet"]:
            with self.subTest(pairing_type=pairing_type):
                rng = np.random.default_rng(20260726)
                norb = 1
                nd = norb * norb
                beta = 2.0
                green, T, N = _toy_green_4x4(Nx=4, Ny=4, Nz=1, Nmat=6, beta=beta)
                Nx, Ny, Nz = 4, 4, 1
                bset = resolve_interactions(_nn_square(0.3), np.eye(3), norb=norb)
                B = bset.n_channels
                ND = nd * B

                # random Hermitian-per-q chi_s/chi_c, real random S/C, random Vpp.
                def _herm(shape):
                    A = rng.normal(size=shape) + 1j * rng.normal(size=shape)
                    return 0.5 * (A + np.conj(np.swapaxes(A, -1, -2)))

                chi_s = _herm((Nx, Ny, Nz, ND, ND))
                chi_c = _herm((Nx, Ny, Nz, ND, ND))
                S_bond = rng.normal(size=(Nx, Ny, Nz, ND, ND)).astype(complex)
                C_bond = rng.normal(size=(Nx, Ny, Nz, ND, ND)).astype(complex)
                Vpp_s = rng.normal(size=(ND, ND)).astype(complex)
                Vpp_t = rng.normal(size=(ND, ND)).astype(complex)

                A, vec_size = make_bond_kernel(
                    chi_s, chi_c, S_bond, C_bond, Vpp_s, Vpp_t,
                    green, bset, pairing_type, beta)
                assert vec_size == norb * norb * Nx * Ny * Nz

                F_q = _fluctuation_vertex(pairing_type, S_bond, chi_s, C_bond, chi_c)
                Vpp = Vpp_s if pairing_type == "singlet" else Vpp_t

                for _ in range(3):
                    phi = (rng.normal(size=(norb, norb, Nx, Ny, Nz))
                           + 1j * rng.normal(size=(norb, norb, Nx, Ny, Nz)))
                    got = A.matvec(phi.ravel()).reshape(norb, norb, Nx, Ny, Nz)
                    ref = _direct_bond_kernel(F_q, Vpp, green, bset, beta, phi)
                    np.testing.assert_allclose(got, ref, atol=1e-10)

    def test_bond_kernel_reference_oracle_two_orbital(self):
        """Spec S7.5, norb=2: the same oracle at two orbitals with a non-symmetric
        Green function and orbital-off-diagonal bonds, so an index transposition
        (l1<->l2 / row<->col in F or Vpp, or in the X = GGphi contraction) that the
        norb=1 case cannot expose is caught.

        Uses a 4x2x1 grid (bonds run along +-x/+-y): with >=3 points along the
        x axis, k_x takes values other than {0, pi}, so the bond phases
        e^{+-ik.dr} are genuinely complex (not just +-1). A 2x2x1 grid (k in
        {0, pi} on every axis) makes every bond phase real and the kernel
        invariant under flipping either phase sign, so it cannot catch a
        phase-sign bug coupled to the orbital-index handling; see review of
        Task 5 (bond-channels-eliashberg)."""
        from hwave.solver.bond_channels import make_bond_kernel

        rng = np.random.default_rng(11)
        norb = 2
        nd = norb * norb
        beta = 2.0
        green, T, N = _toy_green_4x4_2orb(Nx=4, Ny=2, Nz=1, Nmat=4, beta=beta)
        Nx, Ny, Nz = 4, 2, 1
        bset = resolve_interactions(_nn_square_2orb(), np.eye(3), norb=norb)
        B = bset.n_channels
        ND = nd * B

        def _herm(shape):
            A = rng.normal(size=shape) + 1j * rng.normal(size=shape)
            return 0.5 * (A + np.conj(np.swapaxes(A, -1, -2)))

        chi_s = _herm((Nx, Ny, Nz, ND, ND))
        chi_c = _herm((Nx, Ny, Nz, ND, ND))
        S_bond = rng.normal(size=(Nx, Ny, Nz, ND, ND)).astype(complex)
        C_bond = rng.normal(size=(Nx, Ny, Nz, ND, ND)).astype(complex)
        Vpp_s = rng.normal(size=(ND, ND)).astype(complex)
        Vpp_t = rng.normal(size=(ND, ND)).astype(complex)

        for pairing_type in ("singlet", "triplet"):
            A, vec_size = make_bond_kernel(
                chi_s, chi_c, S_bond, C_bond, Vpp_s, Vpp_t,
                green, bset, pairing_type, beta)
            F_q = _fluctuation_vertex(pairing_type, S_bond, chi_s, C_bond, chi_c)
            Vpp = Vpp_s if pairing_type == "singlet" else Vpp_t
            for _ in range(2):
                phi = (rng.normal(size=(norb, norb, Nx, Ny, Nz))
                       + 1j * rng.normal(size=(norb, norb, Nx, Ny, Nz)))
                got = A.matvec(phi.ravel()).reshape(norb, norb, Nx, Ny, Nz)
                ref = _direct_bond_kernel(F_q, Vpp, green, bset, beta, phi)
                np.testing.assert_allclose(got, ref, atol=1e-10)

    def test_bond_kernel_triplet_lambda_minus_one(self):
        """Spec S7.9: single-band 1D 4-point grid, U=4, V=1, chi_s=chi_c=0 (only the
        instantaneous pp part), (T/N) G(k) G(-k) = 1/4.  The odd gap Delta=(0,1,0,-1)
        (= sin k) is an eigenvector with lambda_t = -1; the fluctuation piece is
        zero, so lambda_t^pp = -1 and lambda_t^fl = 0."""
        from hwave.solver.bond_channels import make_bond_kernel

        norb = 1
        Nx, Ny, Nz = 4, 1, 1
        nmat = 1
        beta = 1.0
        # green == 1 everywhere, nmat=1  =>  G2[k] = (1/beta) * 1 * 1 = 1 for all k,
        # so X = G2 phi = phi and the 1/N (=1/4) fold gives (T/N)G G = 1/4.
        green = np.ones((norb, norb, Nx, Ny, Nz, nmat), dtype=complex)

        ci = {((1, 0, 0), (0, 0)): 1.0, ((-1, 0, 0), (0, 0)): 1.0}
        bset = resolve_interactions(ci, np.eye(3), norb=norb)
        assert bset.n_channels == 3

        U, V = 4.0, 1.0
        S0 = np.full((1, 1, 1, 1, 1), U, dtype=complex)
        C0 = np.full((1, 1, 1, 1, 1), 8.0, dtype=complex)  # U + 2 V(q=0) = 4 + 4
        S_bond, C_bond, Vpp_s, Vpp_t = bare_bond_vertices(bset, S0, C0, norb)

        ND = S_bond.shape[-1]
        chi_zero = np.zeros((Nx, Ny, Nz, ND, ND), dtype=complex)

        A, vec_size = make_bond_kernel(
            chi_zero, chi_zero, S_bond, C_bond, Vpp_s, Vpp_t,
            green, bset, "triplet", beta)

        # odd gap Delta = sin k on k in {0, pi/2, pi, -pi/2}
        kx = 2.0 * np.pi * np.arange(Nx) / Nx
        phi = np.sin(kx).reshape(norb, norb, Nx, Ny, Nz).astype(complex)

        out = A.matvec(phi.ravel()).reshape(norb, norb, Nx, Ny, Nz)
        # K phi = -phi  ->  lambda_t = -1
        np.testing.assert_allclose(out, -phi, atol=1e-12)

        # attribution: split the kernel.  fluctuation-only kernel (Vpp=0) annihilates
        # the gap (lambda_t^fl = 0); pp-only kernel (chi=0 already) reproduces the
        # full result (lambda_t^pp = -1).  Rayleigh quotient on the (real) eigenvector.
        zeroVpp = np.zeros_like(Vpp_t)
        A_fl, _ = make_bond_kernel(
            chi_zero, chi_zero, S_bond, C_bond, zeroVpp, zeroVpp,
            green, bset, "triplet", beta)
        out_fl = A_fl.matvec(phi.ravel()).reshape(norb, norb, Nx, Ny, Nz)
        v = phi.ravel()
        denom = np.vdot(v, v)
        lam_fl = (np.vdot(v, out_fl.ravel()) / denom).real
        lam_pp = (np.vdot(v, out.ravel()) / denom).real - lam_fl
        np.testing.assert_allclose(lam_fl, 0.0, rtol=0, atol=1e-12)
        np.testing.assert_allclose(lam_pp, -1.0, rtol=0, atol=1e-12)

    def test_bond_kernel_rejects_unknown_pairing_type(self):
        from hwave.solver.bond_channels import make_bond_kernel
        norb = 1
        beta = 2.0
        green, T, N = _toy_green_4x4(Nx=4, Ny=4, Nz=1, Nmat=4, beta=beta)
        bset = resolve_interactions(_nn_square(0.3), np.eye(3), norb=norb)
        ND = bset.n_channels
        z = np.zeros((4, 4, 1, ND, ND), dtype=complex)
        Vpp = np.zeros((ND, ND), dtype=complex)
        with self.assertRaises(ValueError):
            make_bond_kernel(z, z, z, z, Vpp, Vpp, green, bset, "nonsense", beta)


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

    def test_dress_bond_conditioning_floor_is_tunable(self):
        """The floor is a knob, not a hard-coded refusal: an explicit lower
        ``cond_tol`` lets a deliberately stiff study through."""
        from hwave.solver.bond_channels import dress_bond
        chi_bar, S_bond, C_bond = _diag_denominator_fixture(1.0e-4, "spin")
        chi_s, _ = dress_bond(chi_bar, S_bond, C_bond, cond_tol=1.0e-8)
        assert np.all(np.isfinite(chi_s))


# --- the on-site (Delta r = 0) diagonal must be REAL -----------------------

class TestOnsiteDiagonalMustBeReal(ApproxTestCase):

    def test_resolve_interactions_rejects_a_complex_onsite_diagonal(self):
        """Hermiticity requires V_aa(0) = conj(V_aa(0)), i.e. real. Accepting an
        imaginary diagonal would break the function's Hermiticity-closed
        contract silently."""
        ci = {((0, 0, 0), (0, 0)): 0.5 + 0.2j}
        with self.assertRaisesRegex(ValueError, r"on-site"):
            resolve_interactions(ci, np.eye(3), norb=1)

    def test_resolve_interactions_closes_a_roundoff_imaginary_onsite_diagonal(self):
        """Within the reversal tolerance the imaginary part is round-off; the
        closed value must come out EXACTLY real."""
        ci = {((0, 0, 0), (0, 0)): 0.5 + 1.0e-14j}
        r = resolve_interactions(ci, np.eye(3), norb=1)
        assert r.v_onsite[0, 0].imag == 0.0
        self.assertApprox(r.v_onsite[0, 0].real, 0.5)

    def test_resolve_interactions_keeps_a_complex_offdiagonal_onsite(self):
        """Only the DIAGONAL has to be real; V_01(0) = conj(V_10(0)) may be
        genuinely complex and must still be accepted."""
        ci = {((0, 0, 0), (0, 1)): 0.3 + 0.1j}
        r = resolve_interactions(ci, np.eye(3), norb=2)
        self.assertApprox(r.v_onsite[0, 1], 0.3 + 0.1j)
        self.assertApprox(r.v_onsite[1, 0], 0.3 - 0.1j)


# --- an empty declared CoulombInter ----------------------------------------

class TestEmptyDeclaredCoulombInter(unittest.TestCase):

    def test_resolve_interactions_on_an_empty_mapping(self):
        """A DECLARED but empty CoulombInter is legal (spec S5: a
        declared-but-zero V keeps the channel topology fixed across a V sweep);
        it must resolve to the on-site channel alone."""
        r = resolve_interactions({}, np.eye(3), norb=2)
        assert r.n_channels == 1
        assert r.delta_r == ((0, 0, 0),)
        assert r.reverse == (0,)
        assert np.array_equal(np.asarray(r.v_bond[0]), np.zeros((2, 2)))
        assert np.array_equal(np.asarray(r.v_onsite), np.zeros((2, 2)))
        assert r.dropped == ()

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

    def test_bond_max_shells_cuts_whole_shells_and_records_them(self):
        ci = _two_shell_square()
        nn_only = {(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)}
        nnn_only = {(1, 1, 0), (1, -1, 0), (-1, 1, 0), (-1, -1, 0)}

        full = resolve_interactions(ci, np.eye(3), norb=1)
        assert full.n_channels == 9                      # 1 + 4 + 4

        nn = resolve_interactions(ci, np.eye(3), norb=1, bond_max_shells=1)
        assert nn.n_channels == 5
        assert set(nn.delta_r[1:]) == nn_only
        # cut shells are RECORDED, never silently discarded
        assert {d[0] for d in nn.dropped} == {"shell_cutoff"}
        assert {d[1] for d in nn.dropped} == nnn_only
        # reversal closure survives the cutoff
        for m, dr in enumerate(nn.delta_r):
            assert nn.delta_r[nn.reverse[m]] == tuple(-x for x in dr)

        both = resolve_interactions(ci, np.eye(3), norb=1, bond_max_shells=2)
        assert both.n_channels == 9
        assert set(both.delta_r[1:]) == nn_only | nnn_only
        assert both.dropped == ()

        # asking for more shells than exist is a no-op, not an error
        assert resolve_interactions(ci, np.eye(3), norb=1,
                                    bond_max_shells=5).n_channels == 9


# ---------------------------------------------------------------------------
# Phase-1 round-2 review fixes:
#   R2-1  the conditioning guard must also measure the ABSOLUTE distance to the
#         pole (a 1x1 or a uniformly-small denominator has ratio == 1);
#   R2-2  a tolerance-ACCEPTED reverse pair must be projected onto an exactly
#         conjugate pair (the contract promises Hermiticity-closed output);
#   R2-3  bond_max_shells validation (non-negative integral; the zero-shell
#         ambiguity guard keys on DECLARED nonzero V, not on tol_incl).
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

    def test_onsite_within_tolerance_pair_is_projected_to_exact_hermiticity(self):
        """A declared forward/reverse pair that agrees only WITHIN the reversal
        tolerance must still come back EXACTLY Hermiticity-closed -- keeping both
        declared values verbatim leaves a residual asymmetry in an object whose
        contract says it is closed."""
        v = 0.3 + 0.1j
        eps = 2.0e-9                     # inside atol + rtol*|v| for the defaults
        ci = {((0, 0, 0), (0, 1)): v,
              ((0, 0, 0), (1, 0)): np.conj(v) + eps}
        r = resolve_interactions(ci, np.eye(3), norb=2)
        vo = np.asarray(r.v_onsite)
        assert np.array_equal(vo, vo.conj().T)

    def test_bond_within_tolerance_pair_is_projected_to_exact_hermiticity(self):
        """Same for a genuine (Delta r != 0) bond: v_bond[m] must equal
        conj(v_bond[reverse[m]].T) exactly."""
        v = 0.5 + 0.1j
        eps = 5.0e-9
        ci = {((1, 0, 0), (0, 1)): v,
              ((-1, 0, 0), (1, 0)): np.conj(v) + eps}
        r = resolve_interactions(ci, np.eye(3), norb=2)
        for m in range(r.n_channels):
            a = np.asarray(r.v_bond[m])
            b = np.asarray(r.v_bond[r.reverse[m]])
            assert np.array_equal(a, b.conj().T), "channel {} not closed".format(m)

    def test_bond_max_shells_rejects_a_negative_value(self):
        """A negative bond_max_shells used to be clamped to zero, silently dropping
        every inter-site bond and bypassing the zero-shell ambiguity guard."""
        with self.assertRaisesRegex(ValueError, r"bond_max_shells"):
            resolve_interactions(_nn_square(0.25), np.eye(3), norb=1,
                                 bond_max_shells=-1)

    def test_bond_max_shells_rejects_a_non_integral_value(self):
        with self.assertRaisesRegex(ValueError, r"bond_max_shells"):
            resolve_interactions(_nn_square(0.25), np.eye(3), norb=1,
                                 bond_max_shells=1.5)

    def test_bond_max_shells_zero_rejects_a_tiny_declared_nonzero_v(self):
        """The documented rule is "any DECLARED nonzero inter-site V"; a value
        below tol_incl is still declared and nonzero, so bond_max_shells=0 stays
        ambiguous rather than silently discarding it."""
        with self.assertRaisesRegex(ValueError, r"bond_max_shells=0"):
            resolve_interactions(_nn_square(1.0e-15), np.eye(3), norb=1,
                                 bond_max_shells=0)


# --- R2-4 (optional): validate the documented green layout / beta ----------

def _tiny_kernel_args(Nx=2, nmat=2):
    """Minimal well-formed arguments for make_bond_kernel at norb = 1, B = 1."""
    from hwave.solver.bond_channels import bare_bond_vertices
    bset = resolve_interactions({}, np.eye(3), norb=1)
    ND = 1 * bset.n_channels
    shape = (Nx, 1, 1, ND, ND)
    chi = np.zeros(shape, dtype=complex)
    S = np.zeros(shape, dtype=complex)
    C = np.zeros(shape, dtype=complex)
    Vpp = np.zeros((ND, ND), dtype=complex)
    green = np.ones((1, 1, Nx, 1, 1, nmat), dtype=complex)
    return dict(chi_s=chi, chi_c=chi, S_bond=S, C_bond=C, Vpp_s=Vpp,
                Vpp_t=Vpp, green=green, bond_set=bset,
                pairing_type="singlet", beta=1.0)


def _call_kernel(fn, **over):
    kw = _tiny_kernel_args()
    kw.update(over)
    return fn(kw["chi_s"], kw["chi_c"], kw["S_bond"], kw["C_bond"],
              kw["Vpp_s"], kw["Vpp_t"], kw["green"], kw["bond_set"],
              kw["pairing_type"], kw["beta"])


class TestGreenLayoutAndBetaValidation(unittest.TestCase):

    def test_bond_kernel_validates_green_and_beta(self):
        import hwave.solver.bond_channels as bc
        for fn_name in ["make_bond_kernel", "make_bond_kernel_parts"]:
            for bad in ["ndim", "orbital", "beta"]:
                with self.subTest(fn_name=fn_name, bad=bad):
                    fn = getattr(bc, fn_name)
                    if bad == "ndim":
                        over = {"green": np.ones((1, 1, 2, 1, 1), dtype=complex)}
                    elif bad == "orbital":
                        over = {"green": np.ones((1, 2, 2, 1, 1, 2), dtype=complex)}
                    else:
                        over = {"beta": 0.0}
                    with self.assertRaises(ValueError):
                        _call_kernel(fn, **over)

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
