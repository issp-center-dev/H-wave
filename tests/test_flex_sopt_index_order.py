#!/usr/bin/env python3

"""Pin the orbital index order of the FLEX self-energy against an EXACT ground
truth that is independent of every convention chosen inside H-wave.

Why this file exists
--------------------
The general (full-vertex) FLEX path and its brute-force reference
(``tests/flex_bruteforce_ref.py``) share one convention: the orbital-pair index
order of the bare bubble.  Validating one against the other is therefore
circular, and it stayed circular long enough to ship a transposed self-energy
(issue #91): for a density-only interaction the reduced and general schemes
produced self-energies that agreed on the orbital diagonal but differed by
~20% of the diagonal scale off it, while both susceptibilities still matched to
machine precision.

The ground truth here is exact diagonalization.  A small interacting model is
solved exactly, its self-energy is expanded to O(U^2) by finite differences in
U, and the result is compared against the two candidate index orders of the
second-order skeleton.  Nothing in the comparison depends on how H-wave (or a
paper) labels an index.

The model
---------
One site, three orbitals, intra-orbital U on each.  Three orbitals with a
COMPLEX hopping loop (net flux) are required: with two orbitals the
off-diagonal ``G_{mn}(tau)`` carries a single global phase which factors out,
and the two candidate index orders become identically equal -- a two-orbital
model cannot detect this bug.

The discriminator
-----------------
In imaginary time the second-order skeleton for an intra-orbital U is a
product of three propagators.  The two candidate orbital index orders are

    A_{mn}(tau) = G_{mn}(tau) G_{mn}(tau) G_{nm}(beta - tau)   <- chi0_{mn} G_{mn}
    B_{mn}(tau) = G_{mn}(tau) G_{nm}(tau) G_{mn}(beta - tau)   <- chi0_{nm} G_{mn}

They coincide for m == n, so the diagonal only fixes the overall constant; the
OFF-diagonal decides.  For an intra-orbital U the first-order (Hartree)
self-energy is orbital-diagonal, so its O(U^2) piece (from dn/dU) contaminates
the diagonal with a frequency-independent constant -- the constant is therefore
fitted from diagonal DIFFERENCES between Matsubara frequencies, and the
off-diagonal (Hartree-free) is compared directly.
"""

import contextlib
import unittest
import numpy as np


@contextlib.contextmanager
def _quiet_matmul():
    """numpy's Accelerate-backed complex matmul emits spurious FP warnings on
    these exactly-representable +-1 operator matrices; results are unaffected."""
    with np.errstate(all="ignore"):
        yield

NORB = 3
NSPIN = 2
NSO = NORB * NSPIN
DIM = 1 << NSO

BETA = 4.0
MU = 0.35


def _hopping():
    """3-orbital on-site block with a complex hopping loop (net flux 0.7 rad)."""
    t01 = 0.6 * np.exp(0.2j)
    t12 = 0.5 * np.exp(0.3j)
    t20 = 0.4 * np.exp(0.2j)
    h = np.zeros((NORB, NORB), dtype=complex)
    h[0, 0], h[1, 1], h[2, 2] = 0.2, -0.3, 0.1
    h[0, 1], h[1, 0] = t01, np.conj(t01)
    h[1, 2], h[2, 1] = t12, np.conj(t12)
    h[2, 0], h[0, 2] = t20, np.conj(t20)
    return h


H0 = _hopping()


def _so(m, s):
    return m * NSPIN + s


def _build_annihilators():
    ops = []
    for p in range(NSO):
        c = np.zeros((DIM, DIM))
        for st in range(DIM):
            if not (st >> p) & 1:
                continue
            sign = (-1) ** bin(st & ((1 << p) - 1)).count("1")
            c[st ^ (1 << p), st] = sign
        ops.append(c.astype(complex))
    return ops


C = _build_annihilators()
CD = [c.conj().T for c in C]


def _hamiltonian(U):
    with _quiet_matmul():
        H = np.zeros((DIM, DIM), dtype=complex)
        for s in range(NSPIN):
            for m in range(NORB):
                for n in range(NORB):
                    if H0[m, n] != 0:
                        H = H + H0[m, n] * (CD[_so(m, s)] @ C[_so(n, s)])
        for m in range(NORB):
            nup = CD[_so(m, 0)] @ C[_so(m, 0)]
            ndn = CD[_so(m, 1)] @ C[_so(m, 1)]
            H = H + U * (nup @ ndn)
        N = sum(CD[p] @ C[p] for p in range(NSO))
    return H - MU * N


def _green_ed(U, iws, spin=0):
    """Lehmann representation of G_{mn}(i omega) for the exact model."""
    E, V = np.linalg.eigh(_hamiltonian(U))
    E = E - E.min()
    w = np.exp(-BETA * E)
    Z = w.sum()
    with _quiet_matmul():
        cm = [V.conj().T @ C[_so(m, spin)] @ V for m in range(NORB)]
    boltz = w[:, None] + w[None, :]
    dE = E[None, :] - E[:, None]
    G = np.zeros((len(iws), NORB, NORB), dtype=complex)
    for m in range(NORB):
        for n in range(NORB):
            num = cm[m] * np.conj(cm[n]) * boltz
            for i, iw in enumerate(iws):
                G[i, m, n] = (num / (iw - dE)).sum()
    return G / Z


def _green0_w(iws):
    eye = np.eye(NORB)
    return np.array([np.linalg.inv(iw * eye - (H0 - MU * eye)) for iw in iws])


def _green0_tau(taus):
    """Exact G0_{mn}(tau) on 0 < tau < beta for this finite model."""
    eps, W = np.linalg.eigh(H0 - MU * np.eye(NORB))
    ga = -np.exp(-np.outer(taus, eps)) / (1.0 + np.exp(-BETA * eps))[None, :]
    return np.einsum('ta,ma,na->tmn', ga, W, W.conj())


def _ft_tau_to_iw(f_tau, taus, iws):
    """Simpson quadrature of int_0^beta dtau exp(i omega tau) f(tau)."""
    dt = taus[1] - taus[0]
    wt = np.ones(len(taus))
    wt[1:-1:2] = 4.0
    wt[2:-1:2] = 2.0
    wt *= dt / 3.0
    out = np.zeros((len(iws),) + f_tau.shape[1:], dtype=complex)
    for i, iw in enumerate(iws):
        out[i] = np.tensordot(np.exp(iw * taus) * wt, f_tau, axes=(0, 0))
    return out


def _sigma2_exact(U, iws):
    """O(U^2) coefficient of the exact self-energy, by finite differences.

    Sigma(U) = U Sigma1 + U^2 Sigma2 + O(U^3), so
    (Sigma(2U) - 2 Sigma(U)) / (2 U^2) = Sigma2 + O(U).
    """
    G0w = _green0_w(iws)
    S = {}
    for u in (U, 2 * U):
        Gw = _green_ed(u, iws)
        S[u] = np.array([np.linalg.inv(G0w[i]) - np.linalg.inv(Gw[i])
                         for i in range(len(iws))])
    return (S[2 * U] - 2 * S[U]) / (2 * U * U)


def _candidates(iws, ntau=8001):
    taus = np.linspace(0.0, BETA, ntau)
    Gt = _green0_tau(taus)
    Gr = _green0_tau(BETA - taus)
    A = np.einsum('tmn,tmn,tnm->tmn', Gt, Gt, Gr)
    B = np.einsum('tmn,tnm,tmn->tmn', Gt, Gt, Gr)
    return (_ft_tau_to_iw(A, taus, iws), _ft_tau_to_iw(B, taus, iws))


def _fit_constant(sig2, Aw):
    """Overall constant from diagonal DIFFERENCES (cancels the Hartree shift)."""
    num = den = 0.0
    for m in range(NORB):
        d_ed = sig2[1:, m, m] - sig2[0, m, m]
        d_a = Aw[1:, m, m] - Aw[0, m, m]
        num += np.sum(d_ed * np.conj(d_a))
        den += np.sum(np.abs(d_a) ** 2)
    return num / den


class TestSecondOrderIndexOrder(unittest.TestCase):
    """Exact-diagonalization lock on the orbital index order (issue #91)."""

    @classmethod
    def setUpClass(cls):
        cls.iws = 1j * (2 * np.arange(8) + 1) * np.pi / BETA
        cls.Aw, cls.Bw = _candidates(cls.iws)

    def test_exact_diagonalization_reproduces_g0_at_zero_coupling(self):
        """Guard the ED machinery itself before trusting it as ground truth."""
        np.testing.assert_allclose(_green_ed(0.0, self.iws),
                                   _green0_w(self.iws), atol=1e-12)

    def test_three_orbital_flux_makes_the_candidates_distinguishable(self):
        """Without a net flux the two index orders coincide and the test is
        vacuous, so assert the model actually separates them."""
        scale = np.max(np.abs(self.Aw))
        self.assertGreater(np.max(np.abs(self.Aw - self.Bw)), 1e-3 * scale)

    def test_offdiagonal_selects_the_reduced_index_order(self):
        """A (chi0_{mn} G_{mn}) must match the exact O(U^2) self-energy off the
        orbital diagonal; B (chi0_{nm} G_{mn}) must not."""
        sig2 = _sigma2_exact(2.0e-3, self.iws)
        c = _fit_constant(sig2, self.Aw)
        # the constant is +1 up to the O(U) truncation of the finite difference
        self.assertAlmostEqual(c.real, 1.0, places=2)
        self.assertAlmostEqual(c.imag, 0.0, places=2)

        worst_a = worst_b = 0.0
        for m in range(NORB):
            for n in range(NORB):
                if m == n:
                    continue
                scale = np.max(np.abs(sig2[:, m, n]))
                self.assertGreater(scale, 1e-6)
                worst_a = max(worst_a, np.max(
                    np.abs(sig2[:, m, n] - c * self.Aw[:, m, n])) / scale)
                worst_b = max(worst_b, np.max(
                    np.abs(sig2[:, m, n] - c * self.Bw[:, m, n])) / scale)
        self.assertLess(worst_a, 2.0e-2)
        self.assertGreater(worst_b, 1.0e-1)

    def test_the_selected_order_converges_as_the_coupling_vanishes(self):
        """The residual of A must be pure O(U^3) truncation -- it has to shrink
        proportionally to U.  B's residual is structural and must not."""
        res = {}
        for u in (2.0e-3, 5.0e-4):
            sig2 = _sigma2_exact(u, self.iws)
            c = _fit_constant(sig2, self.Aw)
            ra = rb = 0.0
            for m in range(NORB):
                for n in range(NORB):
                    if m == n:
                        continue
                    scale = np.max(np.abs(sig2[:, m, n]))
                    ra = max(ra, np.max(
                        np.abs(sig2[:, m, n] - c * self.Aw[:, m, n])) / scale)
                    rb = max(rb, np.max(
                        np.abs(sig2[:, m, n] - c * self.Bw[:, m, n])) / scale)
            res[u] = (ra, rb)
        (ra_hi, rb_hi), (ra_lo, rb_lo) = res[2.0e-3], res[5.0e-4]
        # 4x smaller U -> ~4x smaller residual for the correct order
        self.assertLess(ra_lo, 0.4 * ra_hi)
        # ... and no improvement at all for the wrong one
        self.assertGreater(rb_lo, 0.9 * rb_hi)


class TestBruteForceRefOrbitalSlots(unittest.TestCase):
    """Anchor the brute-force reference's ORBITAL slots to the ED result.

    ``tests/flex_bruteforce_ref.py`` cannot be compared to exact
    diagonalization frequency by frequency: it wraps the Matsubara index
    cyclically (``(iw + iv) % nmat``), a deliberate toy that does not carry the
    real fermionic/bosonic phases.  What CAN be anchored -- and what issue #91
    was about -- is the orbital wiring, which is independent of the frequency
    convention.

    Exact diagonalization (above) established that the bubble multiplying
    ``G_{mn}`` is the one whose FORWARD propagator is also ``G_{mn}``.  Composed
    with ``sigma_bruteforce``'s ``V[mu, m, nu, n]`` indexing, a density-only
    vertex picks out ``chi0[m, m, n, n]``, so that element must be the bubble
    built from ``G_{mn}(k+q) G_{nm}(k)`` -- not its transpose.
    """

    def test_density_block_of_chi0_has_the_forward_propagator_first(self):
        from tests.flex_bruteforce_ref import chi0_bruteforce

        rng = np.random.default_rng(91)
        norb, Nk, nmat, T = 3, 3, 4, 1.7
        G = (rng.standard_normal((norb, norb, Nk, nmat))
             + 1j * rng.standard_normal((norb, norb, Nk, nmat)))

        chi0 = chi0_bruteforce(G, T=T, Nk=Nk)

        # The ED-proven density-block bubble, with the SAME cyclic wrap the
        # reference uses (so only the orbital slots are under test).
        for a in range(norb):
            for b in range(norb):
                for q in range(Nk):
                    for iv in range(nmat):
                        s = 0j
                        for k in range(Nk):
                            for iw in range(nmat):
                                s += (G[a, b, (k + q) % Nk, (iw + iv) % nmat]
                                      * G[b, a, k, iw])
                        np.testing.assert_allclose(
                            chi0[a, a, b, b, q, iv], -(T / Nk) * s, atol=1e-12)

    def test_density_only_vertex_contracts_onto_the_matching_green_element(self):
        """Sigma_{mn} must pick up chi0[m,m,n,n] times G_{mn} -- the pairing the
        ED test selected."""
        from tests.flex_bruteforce_ref import sigma_bruteforce

        rng = np.random.default_rng(92)
        norb, Nk, nmat, T = 3, 2, 2, 1.0
        G = (rng.standard_normal((norb, norb, Nk, nmat))
             + 1j * rng.standard_normal((norb, norb, Nk, nmat)))
        chi = (rng.standard_normal((norb, norb, Nk, nmat))
               + 1j * rng.standard_normal((norb, norb, Nk, nmat)))

        # density-only vertex: V[mu, m, nu, n] nonzero only for mu == m, nu == n
        V = np.zeros((norb, norb, norb, norb, Nk, nmat), dtype=complex)
        for m in range(norb):
            for n in range(norb):
                V[m, m, n, n] = chi[m, n]

        sig = sigma_bruteforce(G, V, T=T, Nk=Nk)

        ref = np.zeros((norb, norb, Nk, nmat), dtype=complex)
        for m in range(norb):
            for n in range(norb):
                for k in range(Nk):
                    for iw in range(nmat):
                        s = 0j
                        for q in range(Nk):
                            for iv in range(nmat):
                                s += (chi[m, n, q, iv]
                                      * G[m, n, (k - q) % Nk, (iw - iv) % nmat])
                        ref[m, n, k, iw] = (T / Nk) * s
        np.testing.assert_allclose(sig, ref, atol=1e-12)


if __name__ == '__main__':
    unittest.main()
