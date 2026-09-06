"""Frequency-batched bond dressing, the effective interaction and the
collapses (spec 2026-09-06 sections 3.2-3.3)."""
import unittest

import numpy as np


def _problem(nmat=6, nvol=4, nd=4, B=3, seed=0):
    rng = np.random.default_rng(seed)
    ND = B * nd
    chi_bar = 0.05 * (rng.normal(size=(nmat, nvol, ND, ND)) + 1j * rng.normal(size=(nmat, nvol, ND, ND)))
    S = 0.4 * rng.normal(size=(nvol, ND, ND)); S = 0.5 * (S + np.swapaxes(S, 1, 2))
    C = 0.3 * rng.normal(size=(nvol, ND, ND)); C = 0.5 * (C + np.swapaxes(C, 1, 2))
    return chi_bar, S + 0j, C + 0j


class TestDressBatch(unittest.TestCase):

    def test_equals_per_frequency_dress_channel(self):
        from hwave.solver import bond_channels as bc
        chi_bar, S, C = _problem()
        nmat, nvol = chi_bar.shape[:2]
        for channel, W in (("spin", S), ("charge", C)):
            chi_b, cond = bc.dress_batch(chi_bar[1:4], W, channel, l0=1, nmat=nmat, spatial_shape=(4, 1, 1))
            for i, l in enumerate(range(1, 4)):
                ref, c_ref = bc.dress_channel(chi_bar[l], W, channel, spatial_shape=(4, 1, 1))
                np.testing.assert_array_equal(chi_b[i], ref)
            self.assertIsInstance(cond, float)

    def test_refusal_names_frequency_and_q(self):
        from hwave.solver import bond_channels as bc
        chi_bar, S, C = _problem(nmat=4, nvol=2, nd=1, B=1)
        chi_bar[:] = 0.0
        chi_bar[2, 1, 0, 0] = 1.0                            # l=2, q index 1: 1 - 1*1 = 0
        W = np.ones((2, 1, 1), complex)
        with self.assertRaises(ValueError) as cm:
            bc.dress_batch(chi_bar[2:4], W, "spin", l0=2, nmat=4, spatial_shape=(2, 1, 1))
        msg = str(cm.exception)
        self.assertIn("bosonic", msg)
        self.assertIn(str(2 * 2 - 4), msg)                     # 2l - nmat = 0
        self.assertIn("(1, 0, 0)", msg)
        self.assertIn("spin", msg)


class TestDressAndBuildW(unittest.TestCase):

    def test_w_and_collapses_and_static(self):
        from hwave.solver.flex_bond import BondBlockStore, dress_and_build_w
        chi_bar, S, C = _problem(nmat=6, nvol=4, nd=4, B=3)
        nmat, nvol, ND = chi_bar.shape[:3]
        nd = 4
        for output_full in (False, True):
            names = ("chibar", "W") + (("chi_s_w", "chi_c_w") if output_full else ())
            with BondBlockStore(nmat, nvol, ND, nd, names) as store:
                store.put_freq_batch("chibar", 0, nmat, chi_bar)
                S_on = np.ascontiguousarray(S[:1, :nd, :nd]).repeat(nvol, axis=0) * 0.7
                C_on = np.ascontiguousarray(C[:1, :nd, :nd]).repeat(nvol, axis=0) * 0.3
                res = dress_and_build_w(store, S, C, S_on=S_on, C_on=C_on, nb=4,
                                        output_full=output_full, nmat=nmat,
                                        nvol=nvol, nd=nd, spatial_shape=(4, 1, 1))
                I = np.eye(ND)
                chi_s = np.linalg.solve(I - chi_bar @ S, chi_bar)
                chi_c = np.linalg.solve(I + chi_bar @ C, chi_bar)
                # spec 3.3 rev 19: ring beyond second order + the exact second order
                W_ref = 1.5 * S @ (chi_s - chi_bar) @ S + 0.5 * C @ (chi_c - chi_bar) @ C
                A = S @ chi_bar @ S
                Bc = C @ chi_bar @ C
                W2 = np.zeros_like(W_ref)
                W2[:, :, :nd, :nd] = (1.5 * A + 0.5 * Bc)[:, :, :nd, :nd] \
                    - 0.25 * (S_on + C_on) @ chi_bar[:, :, :nd, :nd] @ (S_on + C_on)
                W2[:, :, :nd, nd:] = 0.25 * (A + Bc)[:, :, :nd, nd:]
                W2[:, :, nd:, :nd] = 0.25 * (A + Bc)[:, :, nd:, :nd]
                W_ref = W_ref + W2
                np.testing.assert_allclose(store.get_freq_batch("W", 0, nmat), W_ref, rtol=1e-12, atol=1e-13)
                np.testing.assert_allclose(res.collapse0, chi_bar[:, :, :nd, :nd], rtol=0, atol=1e-14)
                np.testing.assert_allclose(res.collapse_s, chi_s[:, :, :nd, :nd], atol=1e-12)
                np.testing.assert_allclose(res.collapse_c, chi_c[:, :, :nd, :nd], atol=1e-12)
                np.testing.assert_allclose(res.static_s, chi_s[nmat // 2], atol=1e-12)
                np.testing.assert_allclose(res.static_c, chi_c[nmat // 2], atol=1e-12)
                self.assertTrue(res.static_s.flags.owndata)
                self.assertGreater(res.cond_min_s, 0.0)
                if output_full:
                    np.testing.assert_allclose(store.get_freq_batch("chi_s_w", 0, nmat), chi_s, atol=1e-12)
                    np.testing.assert_allclose(store.get_freq_batch("chi_s_w", nmat // 2, nmat // 2 + 1)[0],
                                               res.static_s, atol=0)
                else:
                    with self.assertRaises(KeyError):
                        store.get_freq_batch("chi_s_w", 0, 1)


if __name__ == "__main__":
    unittest.main()
