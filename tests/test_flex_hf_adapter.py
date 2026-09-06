"""The FLEX-side Hartree-Fock adapter (spec 2026-09-06 section 2.2): the
paramagnetic embedding of the spin-free density into the shared kernel,
the spin-block checks, the owning copy, the admissibility filter on the
table builder, and PairLift's symmetry-zero map (D12)."""
import unittest

import numpy as np


def _rho(nvol, norb, seed=0):
    rng = np.random.default_rng(seed)
    g = rng.normal(size=(nvol, norb, norb)) + 1j * rng.normal(size=(nvol, norb, norb))
    # Hermitian closure rho(r) = conj(rho(-r))^T on a (nvol,1,1) ring
    rev = np.conj(np.swapaxes(g, -1, -2))
    rev = np.roll(rev[::-1], 1, axis=0)
    return 0.5 * (g + rev)


_KANAMORI = {
    "CoulombIntra": {((0, 0, 0), (0, 0)): 2.0, ((0, 0, 0), (1, 1)): 1.5},
    "CoulombInter": {((0, 0, 0), (0, 1)): 0.7, ((0, 0, 0), (1, 0)): 0.7,
                     ((1, 0, 0), (0, 1)): 0.4, ((-1, 0, 0), (1, 0)): 0.4},
    "Hund": {((0, 0, 0), (0, 1)): 0.3, ((0, 0, 0), (1, 0)): 0.3,
             ((1, 0, 0), (0, 0)): 0.1, ((-1, 0, 0), (0, 0)): 0.1},
    "Ising": {((1, 0, 0), (0, 1)): 0.15, ((-1, 0, 0), (1, 0)): 0.15},
    "Exchange": {((0, 0, 0), (0, 1)): 0.3, ((0, 0, 0), (1, 0)): 0.3},
    "PairHop": {((0, 0, 0), (0, 1)): 0.3, ((0, 0, 0), (1, 0)): 0.3},
    "PairLift": {((1, 0, 0), (0, 1)): 0.05, ((-1, 0, 0), (1, 0)): 0.05},
}


class TestHfAdapter(unittest.TestCase):

    def setUp(self):
        self.shape, self.norb = (4, 1, 1), 2
        self.nvol = 4
        self.rho = _rho(self.nvol, self.norb)

    def test_map_is_the_up_block_of_the_spin_major_kernel(self):
        from hwave.solver import flex_hf, hartree_fock as hf
        tabs = flex_hf.build_flex_hf_tables(_KANAMORI, self.norb, self.shape)
        sig = flex_hf.hf_map(self.rho, tabs, self.shape, self.norb)
        rho_so = np.zeros((self.nvol, 2, self.norb, 2, self.norb), complex)
        rho_so[:, 0, :, 0, :] = self.rho
        rho_so[:, 1, :, 1, :] = self.rho
        out = np.zeros((self.nvol, 2 * self.norb, 2 * self.norb), complex)
        hf.accumulate_hf(out, rho_so, tabs.inter_table, tabs.spin_table, self.shape, include_fock=True)
        o = out.reshape(self.nvol, 2, self.norb, 2, self.norb)
        np.testing.assert_array_equal(sig, o[:, 0, :, 0, :])
        self.assertTrue(sig.flags.owndata)                     # owning copy, not a view
        self.assertEqual(np.abs(o[:, 0, :, 1, :]).max(), 0.0)  # paramagnetic: no spin off-diagonal
        np.testing.assert_array_equal(o[:, 0, :, 0, :], o[:, 1, :, 1, :])
        self.assertGreater(np.abs(sig).max(), 0.0)

    def test_hermiticity_in_k_space_is_checked(self):
        from hwave.solver import flex_hf
        tabs = flex_hf.build_flex_hf_tables(_KANAMORI, self.norb, self.shape)
        sig = flex_hf.hf_map(self.rho, tabs, self.shape, self.norb)
        np.testing.assert_allclose(sig, np.conj(np.swapaxes(sig, -1, -2)), atol=1e-12)
        bad = self.rho.copy(); bad[:, 0, 1] += 0.3                # breaks rho(r) = conj(rho(-r))^T
        with self.assertRaises(ValueError):
            flex_hf.hf_map(bad, tabs, self.shape, self.norb)

    def test_pairlift_map_is_zero_on_the_paramagnetic_embedding(self):
        from hwave.solver import flex_hf, hartree_fock as hf
        tabs = flex_hf.build_flex_hf_tables({"PairLift": _KANAMORI["PairLift"]}, self.norb, self.shape)
        sig = flex_hf.hf_map(self.rho, tabs, self.shape, self.norb)    # coherences and R != 0 present
        self.assertLessEqual(np.abs(sig).max(), 1e-14)
        # kernel-level: a spin-off-diagonal density DOES produce a PairLift map
        rho_so = np.zeros((self.nvol, 2, self.norb, 2, self.norb), complex)
        rho_so[:, 0, :, 1, :] = self.rho
        rho_so[:, 1, :, 0, :] = np.conj(np.swapaxes(self.rho, -1, -2))
        out = np.zeros((self.nvol, 4, 4), complex)
        hf.accumulate_hf(out, rho_so, tabs.inter_table, tabs.spin_table, self.shape, include_fock=True)
        self.assertGreater(np.abs(out).max(), 1e-6)

    def test_discarded_coulombintra_is_refused_naming_the_declaration(self):
        from hwave.solver import flex_hf
        ph = {"CoulombIntra": {((0, 0, 0), (0, 0)): 2.0, ((1, 0, 0), (0, 0)): 0.5}}
        with self.assertRaises(ValueError) as cm:
            flex_hf.build_flex_hf_tables(ph, self.norb, self.shape)
        self.assertIn("CoulombIntra", str(cm.exception))
        self.assertIn("(1, 0, 0)", str(cm.exception))
        ph = {"CoulombIntra": {((0, 0, 0), (0, 1)): 0.5}}
        with self.assertRaises(ValueError):
            flex_hf.build_flex_hf_tables(ph, self.norb, self.shape)

    def test_linear_in_each_coefficient(self):
        from hwave.solver import flex_hf
        for t in ("CoulombIntra", "CoulombInter", "Hund", "Ising", "Exchange", "PairHop"):
            with self.subTest(t=t):
                def sig_at(scale):
                    tbl = {k: (v * scale) for k, v in _KANAMORI[t].items()}
                    tabs = flex_hf.build_flex_hf_tables({t: tbl}, self.norb, self.shape)
                    return flex_hf.hf_map(self.rho, tabs, self.shape, self.norb)
                s1, s2 = sig_at(1.0), sig_at(2.0)
                self.assertGreater(np.abs(s1).max(), 1e-8, t)
                np.testing.assert_allclose(s2, 2.0 * s1, rtol=1e-12, atol=1e-14)


if __name__ == "__main__":
    unittest.main()
