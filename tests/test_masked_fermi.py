"""``hartree_fock.masked_fermi``, the distribution function FLEX evaluates.

``_legacy_fermi_occupation`` is the body of ``FLEX._fermi_occupation`` as it
stood before FLEX started calling ``hartree_fock.masked_fermi`` (NumPy branch,
copied verbatim; do not edit it). The shared function must reproduce its
arithmetic bit for bit, including at the cutoff and for arguments that would
overflow ``exp``.
"""
import unittest
from unittest import mock

import numpy as np

from hwave.solver import backend
from hwave.solver import flex
from hwave.solver.hartree_fock import masked_fermi


def _legacy_fermi_occupation(t, mu, ev, ene_cutoff=1.0e2):
    xp = np
    w = (ev - mu) / t
    mask = w < ene_cutoff
    w1 = xp.where(mask, w, 0.0)
    v1 = 1.0 / (1.0 + xp.exp(w1))
    return xp.where(mask, v1, 0.0)


def _same_bits(a, b):
    return a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()


class TestMaskedFermiArithmetic(unittest.TestCase):

    def test_equals_the_former_flex_expression_bit_for_bit(self):
        rng = np.random.default_rng(20261125)
        ev = rng.uniform(-6.0, 6.0, size=(2, 30, 3))      # (nblock, nvol, nd)
        for t in (0.5, 0.1, 0.02, 1.0e-3):
            for mu in (-1.3, 0.0, 0.7):
                with self.subTest(t=t, mu=mu):
                    self.assertTrue(_same_bits(
                        masked_fermi(t, mu, ev),
                        _legacy_fermi_occupation(t, mu, ev)))

    def test_the_cutoff_argument_is_honoured(self):
        ev = np.linspace(-3.0, 3.0, 61)
        for cutoff in (5.0, 40.0, 1.0e2):
            with self.subTest(cutoff=cutoff):
                self.assertTrue(_same_bits(
                    masked_fermi(0.05, 0.2, ev, cutoff),
                    _legacy_fermi_occupation(0.05, 0.2, ev, cutoff)))

    def test_values(self):
        f = masked_fermi(0.5, 0.25, np.array([-0.75, 0.25, 1.25]))
        np.testing.assert_array_equal(
            f, 1.0 / (1.0 + np.exp(np.array([-2.0, 0.0, 2.0]))))
        self.assertEqual(f[1], 0.5)

    def test_zero_at_and_above_the_cutoff(self):
        # (ev - mu) / t = 99.5, 100, 100.5 with t = 1, mu = 0
        f = masked_fermi(1.0, 0.0, np.array([99.5, 100.0, 100.5]))
        self.assertGreater(f[0], 0.0)
        self.assertEqual(f[1], 0.0)
        self.assertEqual(f[2], 0.0)

    def test_no_overflow_beyond_the_cutoff(self):
        ev = np.array([-1.0e6, 0.0, 1.0e3, 1.0e6])
        with np.errstate(over="raise", invalid="raise"):
            f = masked_fermi(1.0e-3, 0.0, ev)
        np.testing.assert_array_equal(f, np.array([1.0, 0.5, 0.0, 0.0]))


class TestMaskedFermiArrayModule(unittest.TestCase):

    def test_a_numpy_array_gives_a_numpy_array(self):
        f = masked_fermi(0.5, 0.0, np.zeros((2, 3)))
        self.assertIs(type(f), np.ndarray)
        self.assertEqual(f.shape, (2, 3))

    def test_a_device_array_is_evaluated_by_its_own_module(self):
        """FLEX keeps ``H0_eigenvalue`` on the device under ``gpu = true``;
        the function must then call the device module, not NumPy."""
        calls = []

        class _DeviceArray(np.ndarray):
            pass
        _DeviceArray.__module__ = "cupy._core.core"

        class _FakeCupy:
            @staticmethod
            def where(*args):
                calls.append("where")
                return np.where(*args)

            @staticmethod
            def exp(x):
                calls.append("exp")
                return np.exp(x)

        ev = np.linspace(-2.0, 2.0, 9).view(_DeviceArray)
        with mock.patch.object(backend, "_import_cupy", lambda: _FakeCupy):
            f = masked_fermi(0.5, 0.1, ev)
        self.assertEqual(calls, ["where", "exp", "where"])
        np.testing.assert_array_equal(
            np.asarray(f),
            _legacy_fermi_occupation(0.5, 0.1, np.linspace(-2.0, 2.0, 9)))


class TestFlexUsesTheSharedFunction(unittest.TestCase):

    def test_flex_has_no_copy_of_its_own(self):
        self.assertIs(flex._masked_fermi, masked_fermi)
        self.assertFalse(hasattr(flex.FLEX, "_fermi_occupation"))


if __name__ == "__main__":
    unittest.main()
