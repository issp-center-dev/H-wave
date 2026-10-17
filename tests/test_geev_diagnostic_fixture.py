"""Qualification gates for the geev-path diagnostic fixture (spec
2026-08-28-mu-green-seam-160, Change 2d): (i) nd_block >= 3 so FLEX's
_matsubara_number_operator takes LAPACK geev, not _eigvals_small;
(ii) H0_k has no finer block structure under the project's own
_find_block_diagonal; (iii) plateau-free: dN/dmu >= 1 at both roots.
"""

import unittest
from contextlib import ExitStack

import numpy as np

from tests.equivalence_cells import GEEV_DIAGNOSTIC_FIXTURE
from tests.test_rpa_flex_equivalence_table import build_solver


class TestGeevDiagnosticFixtureQualification(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._stack = ExitStack()
        cls.rpa, cls.rpa_green, _ = build_solver(
            GEEV_DIAGNOSTIC_FIXTURE, "rpa", cls._stack)
        cls.flex, cls.flex_green, _ = build_solver(
            GEEV_DIAGNOSTIC_FIXTURE, "flex", cls._stack)
        cls.rpa._calc_epsilon_k(cls.rpa_green)
        cls.flex._calc_epsilon_k(cls.flex_green)

    @classmethod
    def tearDownClass(cls):
        cls._stack.close()

    def test_nd_block_takes_the_geev_path(self):
        # flex.py's _matsubara_number_operator dispatches on nd <= 2
        # (_eigvals_small closed form); nd >= 3 is the LAPACK geev path.
        nblock, nvol, nd = self.flex.H0_eigenvalue.shape
        self.assertGreaterEqual(nd, 3)

    def test_spin_mode_is_spin_free(self):
        self.assertEqual(self.rpa.spin_mode, "spin-free")

    def test_h0_has_no_finer_block_structure(self):
        # Criterion: the project's own block detector on the assembled
        # H0 (Task 2's stored self.H0_k) finds a single block of size nd.
        H0 = self.rpa.H0_k
        nblock, nvol, nd, _ = H0.shape
        blocks = self.rpa._find_block_diagonal(
            H0.reshape(nblock * nvol, nd, nd))
        # _find_block_diagonal returns None when no finer block
        # structure exists (single block) -- that is the passing case.
        if blocks is not None:
            self.assertEqual(len(blocks), 1)
            self.assertEqual(sorted(blocks[0]), list(range(nd)))
        self.assertTrue(
            blocks is None
            or (len(blocks) == 1 and sorted(blocks[0]) == list(range(nd))),
            "H0 has finer block structure: {!r}".format(blocks))

    def test_plateau_free_at_both_roots(self):
        fx = GEEV_DIAGNOSTIC_FIXTURE
        beta = 1.0 / fx.T
        ncond = (self.rpa.Ncond / 2
                 if self.rpa.spin_mode == "spin-free" else self.rpa.Ncond)
        _dist, mu_rpa = self.rpa._find_mu(ncond, self.rpa.T)
        nblock, nvol, nd = self.flex.H0_eigenvalue.shape
        sigma_zero = np.zeros((nblock, self.flex.nmat, nvol, nd, nd),
                              dtype=np.complex128)
        mu_flex = self.flex._find_mu_dressed(sigma_zero, beta, ncond)
        lam, ew = self.flex._matsubara_number_operator(sigma_zero, beta)
        for mu in (mu_rpa, mu_flex):
            _n, dn = self.flex._number_from_eigs(lam, ew, mu, beta,
                                                 with_deriv=True)
            self.assertGreaterEqual(dn, 1.0)


if __name__ == "__main__":
    unittest.main()
