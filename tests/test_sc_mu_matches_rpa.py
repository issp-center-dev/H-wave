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


def _bisect_reference(w, T, target, lo, hi):
    """Plain bisection of the particle-number residual down to one ulp.

    An independent reference for the cases RPA's own search cannot solve
    (see the near-full case below); it shares only the masked Fermi counter,
    not the bracketing or polishing logic under test.
    """
    def f(mu):
        return _masked_fermi_delta_n(w, T, mu, target, 100.0)[0]

    assert f(lo) <= 0.0 <= f(hi)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid == lo or mid == hi:
            break
        if f(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _target_with_root_at(ev, beta, endpoint, norb):
    """Filling whose root sits EXACTLY on a band edge.

    ``delta_n`` at that edge is then exactly 0.0, so the strict
    ``product < 0`` bracketing test in ``_determine_mu`` rejects the
    eigenvalue span and routes the search through the wide-bracket
    fallback.  ``nvol * norb`` is a power of two on the fixtures used here,
    so ``n_target * norb * nvol`` reproduces the sum bit-for-bit.
    """
    nvol = ev.shape[0] * ev.shape[1] * ev.shape[2]
    w = ev.reshape(1, nvol, norb)
    edge = float(np.min(ev)) if endpoint == "lo" else float(np.max(ev))
    total, _ = _masked_fermi_delta_n(w, 1.0 / beta, edge, 0.0, 100.0)
    n_target = total / (norb * nvol)
    assert n_target * norb * nvol == total
    return n_target, edge


class TestScMuMatchesRpa(unittest.TestCase):
    # (Nx, Ny, norb, n_target, beta, mu_ceiling, residual_ceiling).  The last
    # row is a large-energy-scale / low-T case: #160 pins the RPA residual
    # ceiling at 1e-10 for that regime, 1e-12 for the normal cases, and the
    # two searches must agree at least as tightly.
    CASES = [
        (8, 8, 1, 0.5, 10.0, 1.0e-12, 1.0e-12),
        (12, 12, 1, 0.375, 40.0, 1.0e-12, 1.0e-12),
        (6, 6, 2, 0.25, 20.0, 1.0e-12, 1.0e-12),
        (16, 16, 1, 0.5, 100.0, 1.0e-10, 1.0e-10),
    ]

    def test_mu_agrees_with_rpa_to_roundoff(self):
        for (Nx, Ny, norb, n_target, beta, mu_ceiling,
             res_ceiling) in self.CASES:
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
                self.assertLessEqual(abs(mu_sc - mu_rpa), mu_ceiling)

                # and mu_sc must satisfy particle number to round-off, judged
                # with the RPA masked-Fermi counter at the matching target
                res_sc, _ = _masked_fermi_delta_n(w, T, mu_sc, Ncond, 100.0)
                self.assertLessEqual(abs(res_sc), res_ceiling)

    def test_wide_bracket_fallback_near_empty_matches_rpa(self):
        # A near-empty band: the root sits BELOW every eigenvalue, so the
        # eigenvalue span is not sign-bracketed and _determine_mu takes the
        # wide (lo-10, hi+10) bracket and the unbracketed polish.
        Nx, Ny, norb, n_target, beta = 8, 8, 1, 1.0e-6, 10.0
        ev = _cosine_band(Nx, Ny, norb)
        nvol = Nx * Ny
        T = 1.0 / beta
        w = ev.reshape(1, nvol, norb)
        Ncond = nvol * n_target * norb

        d_lo, _ = _masked_fermi_delta_n(w, T, float(np.min(ev)), Ncond, 100.0)
        d_hi, _ = _masked_fermi_delta_n(w, T, float(np.max(ev)), Ncond, 100.0)
        self.assertGreater(d_lo * d_hi, 0.0,
                           "fixture must NOT be bracketed by the span")

        mu_sc = _determine_mu(ev, beta, n_target, norb)
        self.assertLess(mu_sc, float(np.min(ev)))

        _dist, mu_rpa = RPA._find_mu(_Host(w), Ncond, T)
        self.assertLessEqual(abs(mu_sc - mu_rpa), 1.0e-12)

        res_sc, _ = _masked_fermi_delta_n(w, T, mu_sc, Ncond, 100.0)
        self.assertLessEqual(abs(res_sc), 1.0e-12)

    def test_wide_bracket_fallback_near_full_is_accurate(self):
        # The mirror case: a near-full band, root ABOVE every eigenvalue.
        # Before #218 RPA's fallback (a secant search started at the band
        # bottom) could not solve this one -- delta_n flattens onto its
        # horizontal asymptote above the band -- so the reference is an
        # independent bisection of the same particle-number residual;
        # since #218 RPA retries on a widened bracket and is compared too.
        Nx, Ny, norb, n_target, beta = 8, 8, 1, 1.0 - 1.0e-6, 10.0
        ev = _cosine_band(Nx, Ny, norb)
        nvol = Nx * Ny
        T = 1.0 / beta
        w = ev.reshape(1, nvol, norb)
        Ncond = nvol * n_target * norb
        lo, hi = float(np.min(ev)), float(np.max(ev))

        d_lo, _ = _masked_fermi_delta_n(w, T, lo, Ncond, 100.0)
        d_hi, _ = _masked_fermi_delta_n(w, T, hi, Ncond, 100.0)
        self.assertGreater(d_lo * d_hi, 0.0,
                           "fixture must NOT be bracketed by the span")

        mu_sc = _determine_mu(ev, beta, n_target, norb)
        self.assertGreater(mu_sc, hi)

        mu_ref = _bisect_reference(w, T, Ncond, hi, hi + 10.0)
        res_sc, dn_sc = _masked_fermi_delta_n(w, T, mu_sc, Ncond, 100.0)
        self.assertLessEqual(abs(res_sc), 1.0e-12)

        # Above the band the residual reaches exactly 0.0 over a PLATEAU of
        # mu values one counting quantum wide: ulp(Ncond) / dn, about 1.1e-11
        # on this fixture.  Any two searches that both drive the residual to
        # round-off can land anywhere in it, so the mu ceiling here is that
        # width (never tighter than the 1e-12 used elsewhere), not a float
        # comparison of two equally exact roots.
        plateau = float(np.spacing(Ncond) / dn_sc)
        self.assertLessEqual(abs(mu_sc - mu_ref),
                             max(1.0e-12, 2.0 * plateau))

        _dist, mu_rpa = RPA._find_mu(_Host(w), Ncond, T)
        self.assertGreater(mu_rpa, hi)
        res_rpa, _ = _masked_fermi_delta_n(w, T, mu_rpa, Ncond, 100.0)
        self.assertLessEqual(abs(res_rpa), 1.0e-12)
        self.assertLessEqual(abs(mu_sc - mu_rpa),
                             max(1.0e-12, 2.0 * plateau))

    def test_endpoint_root_at_band_bottom_matches_rpa(self):
        # delta_n(lo) == 0.0 EXACTLY: the strict product < 0 bracketing test
        # rejects the span, so the wide-bracket fallback must still find the
        # root -- which is the band bottom itself.
        Nx, Ny, norb, beta = 8, 8, 1, 10.0
        ev = _cosine_band(Nx, Ny, norb)
        n_target, edge = _target_with_root_at(ev, beta, "lo", norb)
        nvol = Nx * Ny
        T = 1.0 / beta
        w = ev.reshape(1, nvol, norb)
        Ncond = nvol * n_target * norb

        d_lo, _ = _masked_fermi_delta_n(w, T, edge, Ncond, 100.0)
        self.assertEqual(d_lo, 0.0)

        mu_sc = _determine_mu(ev, beta, n_target, norb)
        self.assertLessEqual(abs(mu_sc - edge), 1.0e-12)

        _dist, mu_rpa = RPA._find_mu(_Host(w), Ncond, T)
        self.assertLessEqual(abs(mu_sc - mu_rpa), 1.0e-12)

        res_sc, _ = _masked_fermi_delta_n(w, T, mu_sc, Ncond, 100.0)
        self.assertLessEqual(abs(res_sc), 1.0e-12)

    def test_endpoint_root_at_band_top(self):
        # The mirror endpoint: delta_n(hi) == 0.0 exactly.  The strict
        # product < 0 test rejects the span here too; before #218 RPA's
        # fallback could not solve an at-or-above-band root, so the primary
        # assertion is the root itself, which the construction of the
        # target pins exactly; since #218 RPA's widened bracket finds it
        # and is compared as well.
        Nx, Ny, norb, beta = 8, 8, 1, 10.0
        ev = _cosine_band(Nx, Ny, norb)
        n_target, edge = _target_with_root_at(ev, beta, "hi", norb)
        nvol = Nx * Ny
        T = 1.0 / beta
        w = ev.reshape(1, nvol, norb)
        Ncond = nvol * n_target * norb

        d_hi, _ = _masked_fermi_delta_n(w, T, edge, Ncond, 100.0)
        self.assertEqual(d_hi, 0.0)

        mu_sc = _determine_mu(ev, beta, n_target, norb)
        self.assertLessEqual(abs(mu_sc - edge), 1.0e-12)

        res_sc, dn_sc = _masked_fermi_delta_n(w, T, mu_sc, Ncond, 100.0)
        self.assertLessEqual(abs(res_sc), 1.0e-12)

        _dist, mu_rpa = RPA._find_mu(_Host(w), Ncond, T)
        res_rpa, _ = _masked_fermi_delta_n(w, T, mu_rpa, Ncond, 100.0)
        self.assertLessEqual(abs(res_rpa), 1.0e-12)
        plateau = float(np.spacing(Ncond) / dn_sc)
        self.assertLessEqual(abs(mu_rpa - edge), max(1.0e-12, 2.0 * plateau))


if __name__ == "__main__":
    unittest.main()
