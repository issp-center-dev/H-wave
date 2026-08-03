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


class TestLayoutsAndDegenerateShapes(unittest.TestCase):
    """The helpers are pure element selections and must be layout- and
    shape-agnostic: Fortran order, non-contiguous views, and degenerate
    dimensions all reduce to the same values as a fresh C-order copy."""

    def _check(self, ns, norb, nvol):
        rng = np.random.default_rng(ns * 100 + norb * 10 + nvol)
        nd = ns * norb
        shape = (nvol,) + (nd,) * 4
        h = (rng.standard_normal(shape)
             + 1j * rng.standard_normal(shape))
        want_p = project_density_pairs(h.copy(), nvol, nd, np)
        want_s = project_density_squashed(h.copy(), nvol, ns, norb, np)
        for label, variant in (
                ("fortran", np.asfortranarray(h)),
                ("noncontig", np.ascontiguousarray(
                    np.concatenate([h, h], axis=0))[:nvol]),
        ):
            with self.subTest(layout=label, ns=ns, norb=norb, nvol=nvol):
                np.testing.assert_array_equal(
                    project_density_pairs(variant, nvol, nd, np), want_p)
                np.testing.assert_array_equal(
                    project_density_squashed(variant, nvol, ns, norb, np),
                    want_s)

    def test_layout_and_shape_matrix(self):
        for ns, norb, nvol in ((2, 1, 1), (2, 1, 4), (2, 3, 2), (1, 2, 3)):
            self._check(ns, norb, nvol)


if __name__ == "__main__":
    unittest.main()
