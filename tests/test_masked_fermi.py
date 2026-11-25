"""``hartree_fock.masked_fermi`` and the particle-number counter built on it.

RPA, FLEX, hwave_sc and UHFk evaluate the Fermi function through
``hartree_fock.masked_fermi``. Before, the same expression was written out
in ``FLEX._fermi_occupation``, in a closure of ``RPA._find_mu``, in
``rpa._masked_fermi_delta_n`` and in two closures of UHFk.

``_legacy_masked_fermi`` is that expression and
``_legacy_masked_fermi_delta_n`` is the former body of the counter (both
copied verbatim from ``rpa.py``), and ``_legacy_find_dist_group_nonzero_t``
is the former body of the UHFk mu search with its closure (copied verbatim
from ``uhfk.py``). Do not edit them. The shared functions must reproduce
their arithmetic bit for bit, at every call site and with the cutoff that
site passes.
"""
import logging
import unittest
from unittest import mock

import numpy as np

from hwave.solver import hartree_fock
from hwave.solver.eliashberg.band import _determine_mu
from hwave.solver.eliashberg.so_solve import _so_own_state
from hwave.solver.hartree_fock import masked_fermi, masked_fermi_delta_n
from hwave.solver.rpa import RPA
from hwave.solver.uhfk import UHFk
from tests.fake_device import fake_cupy, is_device


def _legacy_masked_fermi(t, mu, ev_, ene_cutoff=1.0e2):
    w_ = (ev_ - mu) / t
    mask_ = w_ < ene_cutoff
    w1_ = np.where(mask_, w_, 0.0)
    v1_ = 1.0 / (1.0 + np.exp(w1_))
    return np.where(mask_, v1_, 0.0)


def _legacy_masked_fermi_delta_n(w, T, mu, target, ene_cutoff):
    x = (w - mu) / T
    mask = x < ene_cutoff
    x1 = np.where(mask, x, 0.0)
    f = np.where(mask, 1.0 / (1.0 + np.exp(x1)), 0.0)
    n = float(np.sum(f).real) - target
    dn = float(np.sum(f * (1.0 - f)).real / T)
    return n, dn


logger = logging.getLogger(__name__)


def _legacy_find_dist_group_nonzero_t(self, ws_list, vs_list, ncond):
    from scipy import optimize

    # Collect all eigenvalues for bracket search
    all_ev = np.sort(np.concatenate([w.flatten() for w in ws_list]))
    occupied_number = ncond

    def _fermi(t, mu, ev):
        w_ = (ev - mu) / t
        mask_ = w_ < self.ene_cutoff
        w1_ = np.where( mask_, w_, 0.0 )
        v1_ = 1.0 / (1.0 + np.exp(w1_))
        v_ = np.where( mask_, v1_, 0.0 )
        return v_

    def _calc_delta_n(mu):
        nn = 0.0
        for w, v in zip(ws_list, vs_list):
            ff = _fermi(self.T, mu, w)
            nn += np.einsum('kal,kl,kal->', np.conjugate(v), ff, v).real
        return nn - occupied_number

    # find mu s.t. <n>(mu) = N0
    is_converged = False
    if (_calc_delta_n(all_ev[0]) * _calc_delta_n(all_ev[-1])) < 0.0:
        logger.debug("+++ find mu: try bisection")
        mu, r = optimize.bisect(_calc_delta_n, all_ev[0], all_ev[-1], full_output=True, disp=False)
        is_converged = r.converged
    if not is_converged:
        logger.debug("+++ find mu: try newton")
        mu, r = optimize.newton(_calc_delta_n, all_ev[0], full_output=True)
        is_converged = r.converged
    if not is_converged:
        logger.error("find mu: not converged. abort")
        exit(1)

    dists = [_fermi(self.T, mu, w) for w in ws_list]
    return dists, mu


class _BitsMixin:

    def assertSameBits(self, got, expected):
        got, expected = np.asarray(got), np.asarray(expected)
        self.assertEqual(got.dtype, expected.dtype)
        self.assertEqual(got.shape, expected.shape)
        if got.tobytes() != expected.tobytes():
            flat_g, flat_e = got.ravel(), expected.ravel()
            where = [i for i in range(flat_g.size)
                     if flat_g[i:i + 1].tobytes() != flat_e[i:i + 1].tobytes()]
            i = where[0]
            self.fail("{} of {} entries differ in their bits; the first is entry {}: "
                      "{!r} != {!r}".format(len(where), flat_g.size, i,
                                            flat_g[i], flat_e[i]))


def _levels(shape=(2, 30, 3), seed=20261125):
    return np.random.default_rng(seed).uniform(-6.0, 6.0, size=shape)


class TestMaskedFermiArithmetic(_BitsMixin, unittest.TestCase):

    def test_equals_the_former_expression_bit_for_bit(self):
        ev = _levels()                                    # (nblock, nvol, nd)
        for t in (0.5, 0.1, 0.02, 1.0e-3):
            for mu in (-1.3, 0.0, 0.7):
                for cutoff in (5.0, 20.0, 1.0e2):
                    with self.subTest(t=t, mu=mu, cutoff=cutoff):
                        self.assertSameBits(
                            masked_fermi(t, mu, ev, cutoff),
                            _legacy_masked_fermi(t, mu, ev, cutoff))

    def test_the_default_cutoff_is_100(self):
        ev = _levels()
        self.assertSameBits(masked_fermi(0.02, 0.3, ev),
                            _legacy_masked_fermi(0.02, 0.3, ev, 1.0e2))

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


class TestMaskedFermiArrayModule(_BitsMixin, unittest.TestCase):

    def test_a_host_array_gives_a_host_array(self):
        f = masked_fermi(0.5, 0.0, np.zeros((2, 3)))
        self.assertIs(type(f), np.ndarray)
        self.assertEqual(f.shape, (2, 3))

    def test_a_device_array_stays_on_the_device(self):
        """FLEX keeps ``H0_eigenvalue`` on the device under ``gpu = true``.
        The fake device module refuses a host array mixed into a device
        operation and returns device arrays, so the result is a device
        array only if the function evaluated it with the device module."""
        host = np.linspace(-2.0, 2.0, 9)
        with fake_cupy() as xp:
            f = masked_fermi(0.5, 0.1, xp.asarray(host))
            self.assertTrue(is_device(f))
            got = xp.asnumpy(f)
        self.assertSameBits(got, _legacy_masked_fermi(0.5, 0.1, host))


class TestParticleNumberCounter(unittest.TestCase):

    def test_equals_the_former_counter_exactly(self):
        w = _levels()
        for T in (0.5, 0.1, 0.02, 1.0e-3):
            for mu in (-1.3, 0.0, 0.7):
                for cutoff in (20.0, 1.0e2):
                    with self.subTest(T=T, mu=mu, cutoff=cutoff):
                        # tuples of Python floats: equality is bitwise
                        self.assertEqual(
                            masked_fermi_delta_n(w, T, mu, 37.5, cutoff),
                            _legacy_masked_fermi_delta_n(w, T, mu, 37.5, cutoff))


class _RPAHost:
    """The two attributes ``RPA._find_mu`` reads."""

    def __init__(self, levels, ene_cutoff):
        self.H0_eigenvalue = levels
        self.ene_cutoff = ene_cutoff


class TestFindMuDistribution(_BitsMixin, unittest.TestCase):
    """``RPA._find_mu`` returns the distribution at the root. No caller
    under ``src/`` uses it, so no output file can show a change here."""

    def test_the_returned_distribution_has_the_former_bits(self):
        rng = np.random.default_rng(7)
        for shape, T, ncond, cutoff in (((1, 30, 2), 0.5, 30.0, 1.0e2),
                                        ((2, 30, 3), 0.1, 71.0, 1.0e2),
                                        ((1, 64, 1), 0.02, 20.0, 1.0e2),
                                        ((1, 30, 2), 0.05, 17.0, 20.0)):
            w = rng.uniform(-4.0, 4.0, size=shape)
            with self.subTest(shape=shape, T=T, cutoff=cutoff):
                dist, mu = RPA._find_mu(_RPAHost(w, cutoff), ncond, T)
                self.assertSameBits(dist, _legacy_masked_fermi(T, mu, w, cutoff))
                self.assertLess(abs(float(dist.sum()) - ncond), 1.0e-9)

    def test_the_cutoff_of_the_solver_is_used(self):
        # at T = 0.05 the levels reach (w - mu) / T = 80: masked at 20, kept at 100
        w = np.random.default_rng(7).uniform(-4.0, 4.0, size=(1, 30, 2))
        dist, mu = RPA._find_mu(_RPAHost(w, 20.0), 17.0, 0.05)
        self.assertNotEqual(dist.tobytes(),
                            _legacy_masked_fermi(0.05, mu, w, 1.0e2).tobytes())


class TestEliashbergUsesTheCounter(unittest.TestCase):
    """``hwave_sc`` roots and reports the count of
    ``hartree_fock.masked_fermi_delta_n`` with the cutoff 100."""

    def setUp(self):
        self.eig = np.random.default_rng(3).uniform(-2.0, 2.0, size=(6, 5, 1, 2))
        patcher = mock.patch.object(hartree_fock, "masked_fermi_delta_n",
                                    wraps=hartree_fock.masked_fermi_delta_n)
        self.spy = patcher.start()
        self.addCleanup(patcher.stop)

    def test_determine_mu(self):
        beta, n_target, norb = 10.0, 0.75, 2
        mu = _determine_mu(self.eig, beta, n_target, norb)
        self.assertGreater(self.spy.call_count, 2)
        for call in self.spy.call_args_list:
            self.assertEqual(call.args[4], 100.0)
        residual, _ = _legacy_masked_fermi_delta_n(
            self.eig, 1.0 / beta, mu, n_target * norb * 30, 100.0)
        self.assertLess(abs(residual), 1.0e-10)

    def test_so_own_state(self):
        beta, mu = 10.0, 0.1
        state = _so_own_state(self.eig, beta, mu, 8, 1.5)
        self.assertEqual(self.spy.call_count, 1)
        self.assertEqual(self.spy.call_args.args[4], 100.0)
        count, _ = _legacy_masked_fermi_delta_n(self.eig, 1.0 / beta, mu, 0.0, 100.0)
        self.assertEqual(state.density, count / 30)


class TestUHFkSites(_BitsMixin, unittest.TestCase):
    """The two UHFk sites, on a converged two-orbital chain with
    ``ene_cutoff = 20`` at ``T = 0.05``: levels with ``(w - mu) / T``
    between 20 and 100 are masked with the solver's cutoff and would not be
    with the default, so a site that dropped ``ene_cutoff`` fails here."""

    T = 0.05
    CUTOFF = 20.0

    @classmethod
    def setUpClass(cls):
        ham = {
            "Geometry": {"norb": 2, "rvec": np.eye(3),
                         "center": np.zeros((2, 3)), "degree": np.ones(2)},
            "Transfer": {
                ((0, 0, 0), (0, 0)): -0.4, ((0, 0, 0), (1, 1)): 0.6,
                ((0, 0, 0), (0, 1)): -0.3, ((0, 0, 0), (1, 0)): -0.3,
                ((1, 0, 0), (0, 1)): -1.0, ((-1, 0, 0), (1, 0)): -1.0,
            },
            "Exchange": {((0, 0, 0), (0, 1)): 0.4, ((0, 0, 0), (1, 0)): 0.4},
        }
        params = {"CellShape": [4, 1, 1], "SubShape": [1, 1, 1],
                  "Ncond": 8, "2Sz": None, "IterationMax": 1000,
                  "eps": 1e-10, "Mix": 0.5, "RndSeed": 7,
                  "T": cls.T, "ene_cutoff": cls.CUTOFF}
        cls.solver = UHFk(ham, {"print_level": 0},
                          {"mode": "UHFk", "flag_fock": True}, params)
        cls.solver.solve({"initial_mode": "random"}, ".")
        cls.ws = cls.solver._green_list["eigenvalue"]
        cls.vs = cls.solver._green_list["eigenvector"]

    def test_the_fixture_reaches_between_the_two_cutoffs(self):
        self.assertEqual(self.solver.ene_cutoff, self.CUTOFF)
        self.assertLess(self.solver.physics["Rest"], 1e-10)
        mu = self.solver._green_list["mu"][0]
        x = (np.concatenate([w.ravel() for w in self.ws]) - mu) / self.T
        self.assertTrue(np.any((x > self.CUTOFF) & (x < 1.0e2)))

    def test_the_mu_search(self):
        """The root and the distributions equal those of the former
        method; with the default cutoff the former method finds another
        root and other distributions."""
        dists, mu = self.solver._find_dist_group_nonzero_t(self.ws, self.vs, 8)
        ref_dists, ref_mu = _legacy_find_dist_group_nonzero_t(
            self.solver, self.ws, self.vs, 8)
        self.assertEqual(mu, ref_mu)
        self.assertEqual(len(dists), len(ref_dists))
        for dist, ref in zip(dists, ref_dists):
            self.assertSameBits(dist, ref)

        class _DefaultCutoff:
            T = self.T
            ene_cutoff = 1.0e2
        other_dists, other_mu = _legacy_find_dist_group_nonzero_t(
            _DefaultCutoff, self.ws, self.vs, 8)
        self.assertNotEqual(mu, other_mu)
        self.assertTrue(any(d.tobytes() != o.tobytes()
                            for d, o in zip(dists, other_dists)))

    def test_the_band_energy(self):
        solver = self.solver
        solver._calc_energy()
        mus = solver._green_list["mu"]

        def band(cutoff):
            e = 0.0
            for k, w in enumerate(self.ws):
                mu = mus[solver.block_to_group[k]]
                e += np.sum(w * _legacy_masked_fermi(self.T, mu, w, cutoff))
            return np.real(e)

        self.assertEqual(solver.physics["Ene"]["Band"], band(self.CUTOFF))
        self.assertNotEqual(band(self.CUTOFF), band(1.0e2))


if __name__ == "__main__":
    unittest.main()
