"""The mixed-block weight mask and the pair permutation of the bond gate's
second-order term have ONE home (issue #198): ``flex_bond.mixed_block_mask``
and ``BondDeviceContext.for_view``. Before, the caller and about ten test
modules repeated the same three-line mask recipe inline."""

import types
import unittest
from unittest import mock

import numpy as np

from hwave.solver import flex_bond
from hwave.solver.flex_bond import BondDeviceContext, _mixed_pair_permutation


def _inline_mask(B, nd):
    """The recipe every call site used to spell out."""
    mask = np.zeros((B * nd,) * 2)
    mask[:nd, :] = 0.5
    mask[:, :nd] = 0.5
    mask[:nd, :nd] = 0.0
    return mask


def _vertices(B, norb, seed=0):
    nd = norb * norb
    ND = B * nd
    rng = np.random.default_rng(seed)
    S = rng.normal(size=(ND, ND))
    C = rng.normal(size=(ND, ND))
    S_on = rng.normal(size=(nd, nd))
    C_on = rng.normal(size=(nd, nd))
    return S, C, S_on, C_on


class TestMixedBlockMask(unittest.TestCase):

    def test_matches_the_inline_recipe(self):
        for B, nd in [(1, 1), (1, 4), (3, 1), (3, 4), (5, 9)]:
            with self.subTest(B=B, nd=nd):
                mask = flex_bond.mixed_block_mask(B, nd)
                np.testing.assert_array_equal(mask, _inline_mask(B, nd))
                self.assertEqual(mask.dtype, np.float64)

    def test_weights_are_the_documented_block_pattern(self):
        # channel-0 / channel-0 block 0, channel-0 row and column blocks
        # 0.5, bond / bond blocks 0
        B, nd = 3, 4
        mask = flex_bond.mixed_block_mask(B, nd)
        self.assertTrue(np.all(mask[:nd, :nd] == 0.0))
        self.assertTrue(np.all(mask[:nd, nd:] == 0.5))
        self.assertTrue(np.all(mask[nd:, :nd] == 0.5))
        self.assertTrue(np.all(mask[nd:, nd:] == 0.0))


class TestForView(unittest.TestCase):

    def test_equals_the_manual_construction(self):
        B, norb = 3, 2
        nd = norb * norb
        S, C, S_on, C_on = _vertices(B, norb)
        tail = np.arange(12.0).reshape(3, 4)
        view = types.SimpleNamespace(n_channels=B)
        with BondDeviceContext.for_view(np, S, C, S_on, C_on, view, norb,
                                        green0_tail=tail) as dev, \
                BondDeviceContext(np, S, C, S_on, C_on,
                                  _mixed_pair_permutation(B, nd, norb),
                                  _inline_mask(B, nd), green0_tail=tail) as ref:
            for name in BondDeviceContext._NAMES:
                with self.subTest(name=name):
                    np.testing.assert_array_equal(getattr(dev, name), getattr(ref, name))
            self.assertIs(dev.xp, np)

    def test_without_on_site_vertices(self):
        B, norb = 2, 1
        S, C, _S_on, _C_on = _vertices(B, norb)
        view = types.SimpleNamespace(n_channels=B)
        with BondDeviceContext.for_view(np, S, C, None, None, view, norb) as dev:
            self.assertIsNone(dev.S_on)
            self.assertIsNone(dev.SpC_on)
            np.testing.assert_array_equal(dev.mask, _inline_mask(B, 1))

    def test_uses_the_module_level_permutation(self):
        # the identity-permutation negative control of the second-order
        # tests patches flex_bond._mixed_pair_permutation; the factory must
        # keep looking it up on the module so that control still bites
        B, norb = 2, 2
        nd = norb * norb
        S, C, S_on, C_on = _vertices(B, norb)
        view = types.SimpleNamespace(n_channels=B)
        with mock.patch.object(flex_bond, "_mixed_pair_permutation",
                               lambda B_, nd_, norb_: np.arange(B_ * nd_)):
            with BondDeviceContext.for_view(np, S, C, S_on, C_on, view, norb) as dev:
                np.testing.assert_array_equal(dev.perm, np.arange(B * nd))
        with BondDeviceContext.for_view(np, S, C, S_on, C_on, view, norb) as dev:
            self.assertFalse(np.array_equal(dev.perm, np.arange(B * nd)))

    def test_accepts_a_permutation_override(self):
        B, norb = 2, 2
        nd = norb * norb
        S, C, S_on, C_on = _vertices(B, norb)
        view = types.SimpleNamespace(n_channels=B)
        ident = np.arange(B * nd)
        with BondDeviceContext.for_view(np, S, C, S_on, C_on, view, norb,
                                        perm=ident) as dev:
            np.testing.assert_array_equal(dev.perm, ident)

    def test_refuses_a_vertex_of_the_wrong_pair_dimension(self):
        B, norb = 2, 2
        S, C, S_on, C_on = _vertices(B + 1, norb)       # ND = 3 nd, view says 2
        view = types.SimpleNamespace(n_channels=B)
        with self.assertRaises(ValueError) as cm:
            BondDeviceContext.for_view(np, S, C, S_on, C_on, view, norb)
        self.assertIn("n_channels", str(cm.exception))

    def test_refuses_a_non_square_or_mismatched_vertex(self):
        B, norb = 2, 2
        nd = norb * norb
        S, C, S_on, C_on = _vertices(B, norb)
        view = types.SimpleNamespace(n_channels=B)
        with self.assertRaises(ValueError) as cm:      # S not square
            BondDeviceContext.for_view(np, S[:, :-1], C, S_on, C_on, view, norb)
        self.assertIn("vertex S", str(cm.exception))
        with self.assertRaises(ValueError) as cm:      # C of another pair dimension
            BondDeviceContext.for_view(np, S, np.eye(B * nd + nd), S_on, C_on, view, norb)
        self.assertIn("vertex C", str(cm.exception))

    def test_reads_shapes_without_converting_the_arrays(self):
        # a device array cannot be np.asarray'd; the factory must only look
        # at .shape (BondDeviceContext itself moves the arrays with
        # _bk.to_device). Stand-in: an object exposing shape and refusing
        # conversion, with to_device patched to accept it.
        B, norb = 1, 1
        view = types.SimpleNamespace(n_channels=B)

        class _Opaque:
            shape = (2, 1, 1)
            def __array__(self, *a, **k):
                raise TypeError("implicit conversion is refused")
        with mock.patch.object(flex_bond._bk, "to_device", lambda a, xp: a):
            with BondDeviceContext.for_view(np, _Opaque(), _Opaque(), None, None,
                                            view, norb) as dev:
                self.assertIsInstance(dev.S, _Opaque)


if __name__ == "__main__":
    unittest.main()
