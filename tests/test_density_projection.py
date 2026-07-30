#!/usr/bin/env python3
"""Unit pins for the shared density projections (#107 increment 2).

These are pure element selections, so exact equality is asserted
against the original einsum notations they replaced -- including the
spin-major layout identity between the factorized and combined-index
forms that lets every plain-reduced site share one function.
"""

import unittest

import numpy as np

from hwave.solver.density_projection import (
    project_density_pairs, project_density_squashed)


class TestDensityProjections(unittest.TestCase):

    def setUp(self):
        rng = np.random.default_rng(11)
        self.ns, self.norb, self.nvol = 2, 3, 5
        self.nd = self.ns * self.norb
        shape = (self.nvol,) + (self.nd,) * 4
        self.h = (rng.standard_normal(shape)
                  + 1j * rng.standard_normal(shape))

    def test_pairs_matches_the_combined_index_form(self):
        want = np.einsum(
            'kaabb->kab',
            self.h.reshape(self.nvol, *(self.nd,) * 4)
        ).reshape(self.nvol, self.nd, self.nd)
        got = project_density_pairs(self.h, self.nvol, self.nd, np)
        np.testing.assert_array_equal(got, want)

    def test_pairs_matches_the_factorized_form(self):
        """The spin-major layout identity: 'ksasatbtb->ksatb' through the
        (ns, norb) view selects the same elements in the same order."""
        want = np.einsum(
            'ksasatbtb->ksatb',
            self.h.reshape(self.nvol, *(self.ns, self.norb) * 4)
        ).reshape(self.nvol, self.nd, self.nd)
        got = project_density_pairs(self.h, self.nvol, self.nd, np)
        np.testing.assert_array_equal(got, want)

    def test_squashed_matches_its_original_form(self):
        want = np.einsum(
            'ksauatbvb->ksuatvb',
            self.h.reshape(self.nvol, *(self.ns, self.norb) * 4)
        ).reshape(self.nvol, *(self.ns, self.ns, self.norb) * 2)
        got = project_density_squashed(self.h, self.nvol, self.ns,
                                       self.norb, np)
        np.testing.assert_array_equal(got, want)


if __name__ == "__main__":
    unittest.main()
