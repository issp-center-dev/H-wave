"""The Heff-referenced equal-time density evaluator of the Phase B design
(spec 2026-09-06 section 2.1): UHFk's density convention, exactness for a
frequency-independent self-energy at any Nmat, the O(nmat^-3) truncation
for a dynamic one, and the refusals."""
import unittest

import numpy as np


def _ring(L, norb, seed=0):
    """A complex-hopping ring: H0(k) = t e^{ik} + t^dag e^{-ik} + diag(eps)."""
    rng = np.random.default_rng(seed)
    t = rng.normal(size=(norb, norb)) + 1j * rng.normal(size=(norb, norb))
    eps = np.diag(rng.normal(size=norb))
    k = 2.0 * np.pi * np.arange(L) / L
    H0 = (t[None] * np.exp(1j * k)[:, None, None]
          + t.conj().T[None] * np.exp(-1j * k)[:, None, None] + eps[None])
    return H0                                       # (L, norb, norb)


def _green(H0, sigma_w, mu, beta, nmat):
    """G(k, iw_n) = [(iw_n + mu) - H0 - sigma(iw_n)]^-1, shape (nmat, L, norb, norb)."""
    L, norb = H0.shape[0], H0.shape[1]
    iw = 1j * (2 * np.arange(nmat) + 1 - nmat) * np.pi / beta
    eye = np.eye(norb)
    ginv = (iw[:, None, None, None] + mu) * eye - H0[None] - sigma_w
    return np.linalg.inv(ginv)


def _fermi(e, mu, beta):
    return 1.0 / (np.exp(beta * (e - mu)) + 1.0)


class TestEqualTimeDensity(unittest.TestCase):

    def setUp(self):
        self.L, self.norb, self.beta, self.mu, self.nmat = 6, 2, 4.0, 0.2, 64
        self.H0 = _ring(self.L, self.norb)
        self.shape = (self.L, 1, 1)

    def test_bare_matches_uhfk_green_convention(self):
        from hwave.solver import hartree_fock as hf
        e, V = np.linalg.eigh(self.H0)
        f = _fermi(e, self.mu, self.beta)
        # UHFk _green: rho_ab(k) = sum_j conj(V_aj) f_j V_bj ; rho(r) = fftn(rho_k, norm="forward")
        rho_k_ref = np.einsum('kaj,kj,kbj->kab', V.conj(), f, V)
        rho_r_ref = np.fft.fftn(rho_k_ref.reshape(self.L, 1, 1, self.norb, self.norb),
                                axes=(0, 1, 2), norm="forward").reshape(self.L, self.norb, self.norb)
        G = _green(self.H0, 0.0, self.mu, self.beta, self.nmat)[None]      # (1, nmat, L, norb, norb)
        heff = hf.heff_eigenpairs(self.H0, np.zeros_like(self.H0))
        res = hf.equal_time_density(G, heff, self.mu, self.beta, self.shape)
        np.testing.assert_allclose(res.rho_r, rho_r_ref, rtol=0, atol=1e-12)
        np.testing.assert_allclose(res.rho_k, rho_k_ref, rtol=0, atol=1e-12)
        self.assertFalse(res.rho_r.flags.writeable)
        self.assertFalse(res.rho_k.flags.writeable)
        self.assertIsInstance(res.n_per_spin, float)
        self.assertAlmostEqual(res.n_per_spin, float(f.sum()), places=12)

    def test_static_seed_exact_at_finite_nmat(self):
        from hwave.solver import hartree_fock as hf
        rng = np.random.default_rng(1)
        S = rng.normal(size=(self.L, self.norb, self.norb)) + 1j * rng.normal(size=(self.L, self.norb, self.norb))
        S = 0.5 * (S + S.conj().swapaxes(-1, -2))
        heff = hf.heff_eigenpairs(self.H0, S)
        for nmat in (16, 64):
            G = _green(self.H0, S[None], self.mu, self.beta, nmat)[None]
            res = hf.equal_time_density(G, heff, self.mu, self.beta, self.shape)
            e, _ = np.linalg.eigh(self.H0 + S)
            self.assertAlmostEqual(res.n_per_spin, float(_fermi(e, self.mu, self.beta).sum()), places=12)

    def test_dynamic_tail_exponent(self):
        """Single orbital, two-pole self-energy A/(iw - E): the exact density
        is the Fermi-weighted residue sum of the two poles of G."""
        from hwave.solver import hartree_fock as hf
        L, beta, mu = 8, 3.0, 0.1
        H0 = _ring(L, 1, seed=2)
        h = H0[:, 0, 0].real
        A, E = 0.3, -0.4
        disc = np.sqrt((h - mu - E) ** 2 + 4 * A)
        wp = 0.5 * ((h - mu + E) + disc)
        wm = 0.5 * ((h - mu + E) - disc)
        rp = (wp - E) / (wp - wm)
        rm = (wm - E) / (wm - wp)
        # poles of G measured from mu: G = (iw - E) / ((iw + mu - h)(iw - E) - A); poles at iw = w_pm
        n_exact = float(np.sum(rp * _fermi(wp, 0.0, beta) + rm * _fermi(wm, 0.0, beta)))
        heff = hf.heff_eigenpairs(H0, np.zeros_like(H0))
        errs = []
        for nmat in (64, 128, 256):
            iw = 1j * (2 * np.arange(nmat) + 1 - nmat) * np.pi / beta
            sig = (A / (iw - E))[:, None, None, None] * np.ones((1, L, 1, 1))
            G = _green(H0, sig, mu, beta, nmat)[None]
            res = hf.equal_time_density(G, heff, mu, beta, (L, 1, 1))
            errs.append(abs(res.n_per_spin - n_exact))
        slope = -np.polyfit(np.log([64, 128, 256]), np.log(errs), 1)[0]
        print("tail exponent", slope, errs)
        self.assertGreaterEqual(slope, 2.5)
        self.assertLessEqual(slope, 3.5)

    def test_refusals(self):
        from hwave.solver import hartree_fock as hf
        G = _green(self.H0, 0.0, self.mu, self.beta, self.nmat)[None]
        heff = hf.heff_eigenpairs(self.H0, np.zeros_like(self.H0))
        with self.assertRaises(ValueError):                    # non-Hermitian Heff
            bad = np.zeros_like(self.H0); bad[:, 0, 1] = 1.0
            hf.heff_eigenpairs(self.H0, bad)
        with self.assertRaises(hf.NonFiniteError):
            Gn = G.copy(); Gn[0, 0, 0, 0, 0] = np.nan
            hf.equal_time_density(Gn, heff, self.mu, self.beta, self.shape)
        with self.assertRaises(ValueError):                    # symmetry violation
            Gb = G.copy(); Gb[:, :, 0, 0, 1] += 0.5
            hf.equal_time_density(Gb, heff, self.mu, self.beta, self.shape)


if __name__ == "__main__":
    unittest.main()
