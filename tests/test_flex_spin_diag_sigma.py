"""Regression test: FLEX spin-diag self-energy must pair each Green block
with its own spin slot of V_eff.

In spin-diag mode the up/down chi0q blocks are inflated onto distinct spin
diagonals of V_eff (block g -> slot (g, g)), so the self-energy of Green
block g must be built from V_eff slot g.  Truncating the spin-inflated
sigma with [..., :nd_block, :nd_block] instead selects the up-spin vertex
for BOTH spin blocks.

The check uses linearity: with identical up/down Green blocks and
V_down = 2 * V_up (elementwise), the down self-energy must be exactly
twice the up self-energy.
"""
import unittest

import numpy as np

from hwave.solver.flex import FLEX


def _make_flex_stub(shape, nmat, norb, ns):
    flex = FLEX.__new__(FLEX)

    class LatticeStub:
        pass

    flex.lattice = LatticeStub()
    flex.lattice.shape = shape
    flex.lattice.nvol = shape[0] * shape[1] * shape[2]
    flex.nmat = nmat
    flex.norb = norb
    flex.ns = ns
    flex.nd = norb * ns
    return flex


class TestFlexSpinDiagSelfEnergy(unittest.TestCase):
    def test_down_block_uses_down_spin_vertex(self):
        shape = (2, 1, 1)
        nmat, norb, ns = 4, 1, 2
        nvol = 2
        nd = norb * ns
        flex = _make_flex_stub(shape, nmat, norb, ns)

        rng = np.random.default_rng(42)
        g_block = (rng.standard_normal((nmat, nvol, norb, norb))
                   + 1j * rng.standard_normal((nmat, nvol, norb, norb)))
        green_kw = np.stack([g_block, g_block])  # identical up/down blocks

        v_up = (rng.standard_normal((nmat, nvol, norb, norb))
                + 1j * rng.standard_normal((nmat, nvol, norb, norb)))
        v_eff = np.zeros((nmat, nvol, nd, nd), dtype=np.complex128)
        v_eff[..., :norb, :norb] = v_up
        v_eff[..., norb:, norb:] = 2.0 * v_up

        sigma = flex._calc_self_energy(green_kw, v_eff, beta=10.0)

        self.assertEqual(sigma.shape, (2, nmat, nvol, norb, norb))
        self.assertGreater(np.max(np.abs(sigma[0])), 1e-12)
        np.testing.assert_allclose(
            sigma[1], 2.0 * sigma[0], rtol=0.0, atol=1e-12,
            err_msg="down-spin self-energy must use the down-spin V_eff slot")

    def test_spin_free_single_block_unchanged(self):
        # spin-free: one Green block, both V_eff spin slots identical;
        # the result must equal the elementwise product path.
        shape = (2, 1, 1)
        nmat, norb, ns = 4, 2, 2
        nvol = 2
        nd = norb * ns
        flex = _make_flex_stub(shape, nmat, norb, ns)

        rng = np.random.default_rng(7)
        green_kw = (rng.standard_normal((1, nmat, nvol, norb, norb))
                    + 1j * rng.standard_normal((1, nmat, nvol, norb, norb)))
        v_orb = (rng.standard_normal((nmat, nvol, norb, norb))
                 + 1j * rng.standard_normal((nmat, nvol, norb, norb)))
        v_eff = np.zeros((nmat, nvol, nd, nd), dtype=np.complex128)
        for s in range(ns):
            sl = slice(s * norb, (s + 1) * norb)
            v_eff[..., sl, sl] = v_orb

        sigma = flex._calc_self_energy(green_kw, v_eff, beta=10.0)

        # reference: same-shape elementwise path (nd_block == nd_v)
        flex_ref = _make_flex_stub(shape, nmat, norb, 1)
        sigma_ref = flex_ref._calc_self_energy(green_kw, v_orb, beta=10.0)

        np.testing.assert_allclose(sigma, sigma_ref, rtol=0.0, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
