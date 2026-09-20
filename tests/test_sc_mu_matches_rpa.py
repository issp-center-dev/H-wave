"""Regression test (#185): hwave_sc's chemical-potential search must agree
with RPA's to round-off.

``calc_eliashberg`` computes mu with ``_determine_mu`` (sc.py, called at the
``mu = _determine_mu(eigenvalues, beta, n_filling, norb)`` line) and passes it
downstream as ``precomputed_mu``, bypassing ``RPA._find_mu``.  RPA's search
(brentq xtol=1e-14 + Newton fallback + the transactional Newton polish
``_polish_mu_root``) drives the particle-number residual to round-off; a plain
bisection with default tolerances does not.  This test runs the SAME band
through both mu searches and checks that they agree.

Filling-convention mapping (the two entry points encode filling differently):
``_determine_mu`` targets ``n_target`` electrons per orbital PER SPIN, i.e. the
root of ``sum_f / nvol == n_target * norb`` over one spin block.  RPA's
``_find_mu`` root-finds an UNNORMALIZED sum over all (block, k, band) entries of
one spin block equal to ``Ncond``.  Feeding the identical eigenvalues, the same
physical per-spin filling is therefore ``Ncond = nvol * n_target * norb``.
"""
import unittest

import numpy as np

from hwave.sc import _determine_mu
from hwave.solver.rpa import RPA, _masked_fermi_delta_n


class _Host:
    """Minimal duck-typed RPA host: ``_find_mu`` reads only these two
    attributes (mirrors tests/test_rpa_find_mu.py)."""
    ene_cutoff = 1.0e2

    def __init__(self, levels):
        self.H0_eigenvalue = np.asarray(levels, dtype=np.float64)


def _cosine_band(Nx, Ny, norb):
    kx = np.linspace(0.0, 2.0 * np.pi, Nx, endpoint=False)
    ky = np.linspace(0.0, 2.0 * np.pi, Ny, endpoint=False)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    base = -2.0 * (np.cos(KX) + np.cos(KY))
    ev = np.zeros((Nx, Ny, 1, norb))
    for o in range(norb):
        ev[:, :, 0, o] = base + 0.3 * o  # split orbitals to avoid degeneracy
    return ev


class TestScMuMatchesRpa(unittest.TestCase):
    # (Nx, Ny, norb, n_target, beta, residual_ceiling).  The last row is a
    # large-energy-scale / low-T case: #160 pins the RPA residual ceiling at
    # 1e-10 for that regime, 1e-12 for the normal cases.
    CASES = [
        (8, 8, 1, 0.5, 10.0, 1.0e-12),
        (12, 12, 1, 0.375, 40.0, 1.0e-12),
        (6, 6, 2, 0.25, 20.0, 1.0e-12),
        (16, 16, 1, 0.5, 100.0, 1.0e-10),
    ]

    def test_mu_agrees_with_rpa_to_roundoff(self):
        for Nx, Ny, norb, n_target, beta, res_ceiling in self.CASES:
            with self.subTest(Nx=Nx, norb=norb, n_target=n_target, beta=beta):
                ev = _cosine_band(Nx, Ny, norb)
                nvol = Nx * Ny * 1
                T = 1.0 / beta

                mu_sc = _determine_mu(ev, beta, n_target, norb)

                # same eigenvalues, same physical per-spin filling
                Ncond = nvol * n_target * norb
                w = ev.reshape(1, nvol, norb)
                _dist, mu_rpa = RPA._find_mu(_Host(w), Ncond, T)

                # the two searches must land on the same root
                self.assertLessEqual(abs(mu_sc - mu_rpa), 1.0e-10)

                # and mu_sc must satisfy particle number to round-off, judged
                # with the RPA masked-Fermi counter at the matching target
                res_sc, _ = _masked_fermi_delta_n(w, T, mu_sc, Ncond, 100.0)
                self.assertLessEqual(abs(res_sc), res_ceiling)


if __name__ == "__main__":
    unittest.main()
