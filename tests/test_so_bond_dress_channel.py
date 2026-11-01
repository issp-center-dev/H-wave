"""The ``"spinful"`` channel of the bond dressing.

The spinful (spin-orbital) bond vertex is dressed by the plain RPA
resummation ``chi = [1 + chibar W]^{-1} chibar`` with the full
antisymmetrized vertex ``W``: the channel sign is ``+1``, the same as the
charge channel. The spin and charge channels are unchanged.
"""
import unittest

import numpy as np

import hwave.solver.bond_channels as bc


def _dress_channel_args(channel):
    """Positional arguments of a minimal ``dress_channel`` call."""
    return np.zeros((1, 4, 4), complex), np.zeros((1, 4, 4), complex), channel


class TestDressChannels(unittest.TestCase):
    def test_spinful_sign(self):
        self.assertEqual(bc._DRESS_CHANNELS, {"spin": -1.0, "charge": 1.0, "spinful": 1.0})

    def test_messages_exact(self):
        with self.assertRaises(ValueError) as cm:
            bc.dress_batch(np.zeros((1, 1, 4, 4), complex), np.zeros((1, 4, 4), complex),
                           "spinless", l0=0, nmat=2, spatial_shape=(1, 1, 1))
        self.assertEqual(str(cm.exception),
                         "dress_batch: channel must be 'spin', 'charge' or 'spinful', "
                         "got 'spinless'")
        with self.assertRaises(ValueError) as cm:
            bc.dress_channel(*_dress_channel_args("spinless"), spatial_shape=(1, 1, 1))
        self.assertEqual(str(cm.exception),
                         "dress_channel: channel must be 'spin', 'charge' or 'spinful', "
                         "got 'spinless'")

    def test_spinful_equals_charge_arithmetic(self):
        """The spinful channel is ``[1 + X W]^{-1} X``: the same arithmetic as
        the charge channel, and an explicit inverse agrees."""
        rng = np.random.default_rng(7)
        X = 0.1 * (rng.standard_normal((2, 6, 4, 4)) + 1j * rng.standard_normal((2, 6, 4, 4)))
        W = 0.1 * (rng.standard_normal((6, 4, 4)) + 1j * rng.standard_normal((6, 4, 4)))
        got, _ = bc.dress_batch(X, W, "spinful", l0=0, nmat=2, spatial_shape=(3, 2, 1))
        chg, _ = bc.dress_batch(X, W, "charge", l0=0, nmat=2, spatial_shape=(3, 2, 1))
        self.assertTrue(np.array_equal(got, chg))
        want = np.linalg.inv(np.eye(4) + X @ W[None]) @ X
        self.assertLessEqual(np.abs(got - want).max() / np.abs(want).max(), 1e-13)
        one, _ = bc.dress_channel(X[0], W, "spinful", spatial_shape=(3, 2, 1))
        self.assertLessEqual(np.abs(one - want[0]).max() / np.abs(want[0]).max(), 1e-13)


if __name__ == "__main__":
    unittest.main()
