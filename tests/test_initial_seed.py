"""Field seeding (``hwave.initial.seed``; design section 7, validation R1).

The seed of a pattern O is the density of ``H_base - h O`` with the
solver's blocks and mu rule (``green_function``, uniform T = 0 shell
fill). ``m(h) = Re tr[O rho(h)] = Re sum_ab O_ab G_ab(h)`` per primitive
cell, with ``rho = G^T`` and G the stored Green function (``G_ab =
<c^dag_a c_b>``).
"""
import math
import os
import shutil
import tempfile
import time
import unittest
from fractions import Fraction
from unittest import mock

import numpy as np

from tests import initial_fixtures

_CHAIN = {((1, 0, 0), 0, 0): -1.0, ((-1, 0, 0), 0, 0): -1.0}
_SY = np.array([[0, -1j], [1j, 0]], dtype=complex)


def _spec(entries, q, subshape, *, phase=0.0, sign=1):
    from hwave.initial.patterns import PatternSpec
    return PatternSpec(q=tuple(str(Fraction(x)) for x in q), subshape=subshape, phase=phase,
                       sign=sign, modulation="cos", entries=tuple(entries))


def _e(a=0, b=None, r=(0, 0, 0), **c):
    from hwave.initial.patterns import Entry
    return Entry(a, a if b is None else b, r, **c)


def _green_k(G, shape):
    nvol = G.shape[0]
    nd = G.shape[1] * G.shape[2]
    return np.fft.ifftn(G.reshape(tuple(shape) + (nd, nd)), axes=(0, 1, 2),
                        norm="forward").reshape(nvol, nd, nd)


class _Tmp(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="initial_seed_")
        self._n = 0
        self._t0 = time.perf_counter()

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)
        # each fast test method stays under 5 s (tests/heavy_tests.py)
        self.assertLess(time.perf_counter() - self._t0, 5.0)

    def dest(self):
        self._n += 1
        return os.path.join(self.root, "m{}".format(self._n))

    def facts(self, toml):
        from hwave.initial import facts
        return facts.model_facts(toml)

    def chain(self, L, subshape, *, ncond=None, T=0.0, U=1.0, transfer=None, fock=True,
              two_sz=None, norb=1, centers=None, extra=None, so=False, intra=None):
        p = {"Ncond": L * norb if ncond is None else ncond, "T": T,
             "SubShape": list(subshape)}
        if two_sz is not None:
            p["2Sz"] = two_sz
        mode = {}
        if not fock:
            mode["flag_fock"] = False
        if so:
            mode["enable_spin_orbital"] = True
        nphys = norb // 2 if so else norb
        inter = {"CoulombIntra": intra if intra is not None else
                 {((0, 0, 0), a, a): U for a in range(nphys)}}
        toml = initial_fixtures.write_input(
            self.dest(), cellshape=(L, 1, 1), norb=norb,
            transfer=_CHAIN if transfer is None else transfer, centers=centers,
            interactions=inter, param=p, mode=mode or None, extra_files=extra)
        return self.facts(toml)

    def seed(self, f, spec, **kw):
        from hwave.initial import hints, patterns, seed
        p = f.problem
        O = None if spec is None else patterns.realise(spec, f, p.subshape)
        kw.setdefault("T_run", p.T)
        kw.setdefault("W0", hints.t_hint_bare(p).W0)
        kw.setdefault("seed_base", "bare")
        return seed.seed_pattern(f, p.ham_trans_k, O, **kw), O


# --------------------------------------------------------------------------
# m(h): G versus rho
# --------------------------------------------------------------------------

class TestSigmaY(_Tmp):
    def sigma_y_model(self):
        """One orbital, 2 x 1, Sz free, Fock on (CoulombIntra joins the
        spins in one block). An on-site Zeeman term ``-B sigma^z`` (B = 1/4,
        a 9-column Transfer; no bond fits a length-2 axis) makes
        ``sigma^y`` non-commuting with H_base: at half filling and T = 0
        the lower level of each k is filled, and ``<sigma^y> = h/sqrt(h^2 +
        B^2) = 2/sqrt(5) = 0.894`` at h = 0.5 for the field ``-h
        sigma^y``."""
        entries = {((0, 0, 0), 0, 0, 0, 0): -0.25}
        entries[((0, 0, 0), 0, 0, 1, 1)] = 0.25
        text = initial_fixtures.spin_transfer_text(1, entries)
        return self.chain(2, (1, 1, 1), ncond=2, extra={"transfer.dat": text})

    def test_sigma_y_green_vs_rho(self):
        from hwave.initial import seed
        f = self.sigma_y_model()
        p = f.problem
        self.assertTrue(p.spin_resolved_transfer)
        # the unnormalised sigma^y at every k (spin-block index s*norb + a)
        O = np.broadcast_to(_SY, (p.nvol, 2, 2)).copy()
        G, _ = seed.density_at(f, p.ham_trans_k, O, 0.5, T_run=0.0, W0=1.0)
        m = seed.m_value(O, G, problem=p)
        self.assertAlmostEqual(m, 2 / math.sqrt(5), delta=1e-12)
        self.assertAlmostEqual(m, 0.894, delta=1e-3)
        # the wrong form Re tr[O G] = <O^T> reverses the imaginary sigma^y
        Gk = _green_k(G, p.shape)
        wrong = float(np.real(np.einsum("kab,kba->", O, Gk))) / p.nvol
        self.assertAlmostEqual(wrong, -m, delta=1e-12)
        self.assertAlmostEqual(wrong, -0.894, delta=1e-3)
        # the normalised sigma^y seed: m(h) > 0 and increasing, target met
        res, On = self.seed(f, _spec([_e(cy=1.0)], ("0", "0", "0"), (1, 1, 1)), T_run=0.0)
        self.assertTrue(res.tried)
        self.assertIsNone(res.cap_reason)
        dm_sat = res.m_inf - res.m0
        self.assertLessEqual(abs((res.m_used - res.m0) - 0.5 * dm_sat), 1e-3 * dm_sat + 1e-15)
        ms = [seed.m_value(On, seed.density_at(f, p.ham_trans_k, On, h, T_run=0.0,
                                               W0=1.0)[0], problem=p)
              for h in (0.1, 0.3, 0.5, 1.0, 2.0)]
        self.assertGreater(ms[0], 0)
        self.assertTrue(all(b > a for a, b in zip(ms, ms[1:])), ms)
        # the seeded G read back as the PR 1a Green slot: S^y > 0
        self.assertGreater(seed.moments_of(res.G, f)["spin"][0][1], 0)

    def test_m_value_per_primitive_cell(self):
        from hwave.initial import patterns, seed
        f = self.chain(4, (2, 1, 1), ncond=4, T=0.1)
        p = f.problem
        O = patterns.realise(_spec([_e(c0=1.0)], ("0", "0", "0"), (1, 1, 1)), f, (2, 1, 1))
        G, _ = seed.density_at(f, p.ham_trans_k, O, 0.0, T_run=0.1, W0=4.0)
        # O = n / sqrt(2) per site: m = density per cell / sqrt(2) = 1/sqrt(2)
        # (to the electron count of the solver's T > 0 mu bisection, xtol 2e-12)
        self.assertAlmostEqual(seed.m_value(O, G, problem=p), 1 / math.sqrt(2), delta=1e-10)


# --------------------------------------------------------------------------
# m_inf
# --------------------------------------------------------------------------

class TestMInf(_Tmp):
    def test_m_inf_plateau(self):
        """Staggered charge at half filling: the cut lies in the O gap
        ``delta_O = sqrt(2) >= 0.1``; m(h) approaches m_inf like 1/h^2."""
        from hwave.initial import patterns, seed
        f = self.chain(4, (2, 1, 1), ncond=4, T=0.0)
        p = f.problem
        O = patterns.realise(_spec([_e(c0=1.0)], ("1/2", "0", "0"), (2, 1, 1)), f, (2, 1, 1))
        m_inf = seed.m_inf(O, p.blocks, bvol=2)
        self.assertAlmostEqual(m_inf, 1 / math.sqrt(2), delta=1e-14)
        delta_o = math.sqrt(2)
        w_floor = seed.w_floor(p.ham_trans_k)
        self.assertAlmostEqual(w_floor, 4.0, delta=1e-12)
        h_op = float(max(np.max(np.abs(np.linalg.eigvalsh(p.ham_trans_k))), 0))
        prev = -np.inf
        for fac in (1e4, 1e5, 1e6):
            h = fac * w_floor
            G, _ = seed.density_at(f, p.ham_trans_k, O, h, T_run=0.0, W0=4.0)
            m = seed.m_value(O, G, problem=p)
            self.assertGreaterEqual(m, prev)
            self.assertLessEqual(m_inf - m, 10 * h_op ** 2 / (h ** 2 * delta_o) + 1e-12)
            self.assertLessEqual(m, m_inf + 1e-12)
            prev = m

    def test_cut_degenerate_ferro(self):
        """sigma^z ferro at quarter filling: the cut lies inside the +1
        eigenspace. m(h) = m_inf (to 1e-12 at T = 0; at T > 0 to the electron
        count of the solver's mu bisection, xtol 2e-12), while the
        occupations inside the eigenspace stay selected by H_base: filled
        by energy at T = 0, thermal at T > 0, never uniform."""
        from hwave.initial import patterns, seed
        for T, atol in ((0.0, 1e-12), (0.1, 1e-10)):
            f = self.chain(8, (1, 1, 1), ncond=4, T=T)
            p = f.problem
            O = patterns.realise(_spec([_e(cz=1.0)], ("0", "0", "0"), (1, 1, 1)), f, (1, 1, 1))
            m_inf = seed.m_inf(O, p.blocks, bvol=1)
            self.assertAlmostEqual(m_inf, 4 / math.sqrt(2) / 8, delta=1e-14)
            w_floor = seed.w_floor(p.ham_trans_k)
            for fac in (1e2, 1e3, 1e4):
                G, _ = seed.density_at(f, p.ham_trans_k, O, fac * w_floor, T_run=T, W0=4.0)
                self.assertAlmostEqual(seed.m_value(O, G, problem=p), m_inf, delta=atol)
                occ = np.real(_green_k(G, p.shape)[:, 0, 0])      # up spin per k
                self.assertGreater(float(occ.max() - occ.min()), 0.5)
                np.testing.assert_allclose(np.sum(occ), 4, atol=1e-9)


# --------------------------------------------------------------------------
# the h search
# --------------------------------------------------------------------------

class TestSearch(_Tmp):
    def test_band_insulator_small_O(self):
        """Eight orbitals per cell with alternating on-site energies +-3 (a
        band insulator at half filling) and a ferro sigma^z over all of them:
        ``||O||_op = 1/4``; m(h) stays flat until h exceeds the gap scale, and
        the doubling still reaches the target."""
        from hwave.initial import patterns, seed
        n = 8
        t = {}
        for a in range(n):
            t[((0, 0, 0), a, a)] = 3.0 if a % 2 else -3.0
            b, r = (a + 1) % n, (1, 0, 0) if a == n - 1 else (0, 0, 0)
            t[(r, a, b)] = -1.0
            t[(tuple(-x for x in r), b, a)] = -1.0
        f = self.chain(4, (1, 1, 1), ncond=4 * n, T=0.0, norb=n, transfer=t,
                       centers=[(a / n, 0, 0) for a in range(n)])
        O_spec = _spec([_e(a, cz=1.0) for a in range(n)], ("0", "0", "0"), (1, 1, 1))
        O = patterns.realise(O_spec, f, (1, 1, 1))
        self.assertAlmostEqual(seed.op_norm(O), 0.25, delta=1e-14)
        res, _ = self.seed(f, O_spec)
        self.assertTrue(res.tried)
        self.assertIsNone(res.cap_reason)
        self.assertAlmostEqual(res.m0, 0.0, delta=1e-12)
        self.assertGreater(res.h_used, seed.w_floor(f.problem.ham_trans_k))
        dm_sat = res.m_inf - res.m0
        self.assertTrue(res.m_used - res.m0 >= 0.5 * dm_sat
                        or abs(res.m_used - res.m0 - 0.5 * dm_sat) <= 1e-3 * dm_sat)

    def test_full_filling_no_response(self):
        from hwave.initial import seed
        f = self.chain(4, (2, 1, 1), ncond=8, T=0.0)
        calls = []
        real = seed.density_at

        def counting(*a, **k):
            calls.append(a[3])
            return real(*a, **k)

        with mock.patch.object(seed, "density_at", counting):
            res, _ = self.seed(f, _spec([_e(c0=1.0)], ("1/2", "0", "0"), (2, 1, 1)))
        self.assertEqual([h for h in calls if h != 0], [])
        self.assertFalse(res.tried)
        self.assertEqual(res.not_tried_reason, "no response (degenerate with para)")
        self.assertAlmostEqual(res.m_inf, res.m0, delta=1e-12)

    def test_crystal_field_m0_above_target(self):
        """Orbital polarisation under a crystal field: m(0) exceeds
        ``amplitude * m_inf``, and the target is the change Delta m."""
        from hwave.initial import seed
        t = {((1, 0, 0), a, a): -1.0 for a in range(2)}
        t.update({((-1, 0, 0), a, a): -1.0 for a in range(2)})
        t[((0, 0, 0), 0, 0)] = -1.5
        t[((0, 0, 0), 1, 1)] = 1.5
        f = self.chain(4, (1, 1, 1), ncond=4, T=0.1, norb=2, transfer=t)
        res, _ = self.seed(f, _spec([_e(0, c0=1.0), _e(1, c0=-1.0)], ("0", "0", "0"),
                                    (1, 1, 1)))
        self.assertTrue(res.tried)
        self.assertGreater(res.m0, 0.5 * res.m_inf)
        dm_sat = res.m_inf - res.m0
        self.assertAlmostEqual(res.m_target, res.m0 + 0.5 * dm_sat, delta=1e-14)
        self.assertLessEqual(abs(res.m_used - res.m_target), 1e-3 * dm_sat + 1e-15)
        self.assertGreater(res.h_used, 0)

    def test_flat_band_w_floor(self):
        from hwave.initial import seed
        f = self.chain(4, (2, 1, 1), ncond=4, T=0.1, transfer={((0, 0, 0), 0, 0): 0.0})
        self.assertEqual(seed.w_floor(f.problem.ham_trans_k), 1e-3)
        res, _ = self.seed(f, _spec([_e(c0=1.0)], ("1/2", "0", "0"), (2, 1, 1)))
        self.assertTrue(res.tried)
        self.assertIsNone(res.cap_reason)
        self.assertGreater(res.h_used, 1e-3)
        self.assertGreater(res.h_hi, 1e-3)
        dm_sat = res.m_inf - res.m0
        self.assertLessEqual(abs(res.m_used - res.m_target), 1e-3 * dm_sat)

    def test_origin_invariance_seed(self):
        from hwave.initial import seed
        out = []
        for C in (0.0, 7.3):
            t = dict(_CHAIN)
            t[((0, 0, 0), 0, 0)] = C
            f = self.chain(8, (2, 1, 1), ncond=8, T=0.0, transfer=t)
            res, _ = self.seed(f, _spec([_e(cz=1.0)], ("1/2", "0", "0"), (2, 1, 1)))
            out.append((seed.w_floor(f.problem.ham_trans_k), res))
        (w1, r1), (w2, r2) = out
        self.assertLessEqual(abs(w1 - w2), 1e-12 * w1)
        self.assertEqual(r1.shell, r2.shell)
        self.assertTrue(r1.shell["h0"][0]["open_shell"])
        self.assertLessEqual(abs(r1.h_used - r2.h_used), 1e-12 * r1.h_used)
        self.assertLessEqual(float(np.max(np.abs(r1.G - r2.G))),
                             1e-12 * float(np.max(np.abs(r1.G))))
        self.assertLessEqual(abs(r1.m_used - r2.m_used), 1e-12)

    def test_t0_seed_nonpara_default_amplitude(self):
        from hwave.initial import facts, seed
        toml = initial_fixtures.square(self.dest(), 4, U=2.0, T=0.0, subshape=(2, 2, 1))
        f = facts.model_facts(toml)
        res, _ = self.seed(f, _spec([_e(cz=1.0)], ("1/2", "1/2", "0"), (2, 2, 1)))
        self.assertTrue(res.tried)
        self.assertEqual(res.amplitude_requested, 0.5)
        self.assertGreaterEqual(res.amplitude_used, 0.5 - 1e-3)
        self.assertTrue(res.shell["h0"][0]["open_shell"])
        para, _ = self.seed(f, None)
        self.assertGreater(float(np.max(np.abs(res.G - para.G))), 0.05)
        self.assertGreater(res.m_used - res.m0, 0.1)

    def test_shell_tolerance(self):
        """``shell_tol = 1e-10 * max(1, W0, h ||O||_op)``: the h term keeps
        the tolerance above the rounding of ``H_base - h O`` at large h."""
        from hwave.initial import patterns, seed
        from hwave.solver import initial_facade
        f = self.chain(4, (2, 1, 1), ncond=4, T=0.0)
        p = f.problem
        O = patterns.realise(_spec([_e(c0=1.0)], ("1/2", "0", "0"), (2, 1, 1)), f, (2, 1, 1))
        seen = []
        real = initial_facade.green_function

        def spy(*a, **k):
            seen.append(k["shell_tol"])
            return real(*a, **k)

        with mock.patch.object(initial_facade, "green_function", spy):
            for h, W0 in ((0.0, 0.5), (0.0, 3.0), (2.0, 3.0), (1e4, 3.0)):
                seed.density_at(f, p.ham_trans_k, O, h, T_run=0.0, W0=W0)
        op = seed.op_norm(O)
        self.assertAlmostEqual(op, 1 / math.sqrt(2), delta=1e-14)
        want = [1e-10, 3e-10, 3e-10, 1e-10 * 1e4 * op]
        np.testing.assert_allclose(seen, want, rtol=1e-14)

    def test_cap_doublings_takes_effect(self):
        """Patching the public CAP_DOUBLINGS moves the guard: with 1 the
        seed stops at 2 W_floor."""
        from hwave.initial import seed
        f = self.chain(4, (2, 1, 1), ncond=4, T=0.1)
        spec = _spec([_e(c0=1.0)], ("1/2", "0", "0"), (2, 1, 1))
        free, _ = self.seed(f, spec, amplitude=0.999)
        self.assertIsNone(free.cap_reason)
        self.assertGreater(free.h_used, 2 * seed.w_floor(f.problem.ham_trans_k))
        with mock.patch.object(seed, "CAP_DOUBLINGS", 1):
            res, _ = self.seed(f, spec, amplitude=0.999)
        self.assertEqual(res.cap_reason, "target not reached at the cap")
        self.assertEqual(res.h_used, 2 * seed.w_floor(f.problem.ham_trans_k))
        self.assertEqual(res.h_hi, res.h_used)

    def test_electron_count_at_large_h(self):
        """At T > 0 the solver's mu bisection has a tolerance relative to
        |mu| ~ h, so the raw count error grows with h; density_at refines mu
        and keeps every group's count within N_RTOL * max(1, Ncond), also at
        h far above the plateau (2^24 W_floor; beyond about 2^26 W_floor the
        solver's own mu search fails on this fixture, see density_at)."""
        from hwave.initial import patterns, seed
        f = self.chain(4, (2, 1, 1), ncond=4, T=0.1, two_sz=2)
        p = f.problem
        O = patterns.realise(_spec([_e(cz=1.0)], ("1/2", "0", "0"), (2, 1, 1)), f, (2, 1, 1))
        wf = seed.w_floor(p.ham_trans_k)
        refined = []
        for fac in (1e5, 1e6, 2.0 ** 24):
            G, info, meta = seed.density_at(f, p.ham_trans_k, O, fac * wf, T_run=0.1, W0=4.0,
                                            return_meta=True)
            counts = seed._group_counts(G, p)
            for n, n0 in zip(counts, p.blocks.group_nconds):
                self.assertLessEqual(abs(n - n0), 1e-9 * max(1, p.ncond))
            self.assertLessEqual(meta["n_residual"], 1e-9 * max(1, p.ncond))
            refined.append(meta["mu_refined"])
        self.assertTrue(refined[-1])          # the raw fill was off at large h
        # the seeded row records it
        res, _ = self.seed(f, _spec([_e(cz=1.0)], ("1/2", "0", "0"), (2, 1, 1)))
        self.assertEqual(sorted(res.meta), ["mu_refined", "n_residual"])

    def test_mu_search_failure_is_the_cap(self):
        """On the 2Sz = 2 chain the solver's own T > 0 mu search fails beyond
        about 2^26 W_floor. A target that is never met (m scaled down by a
        wrapper, so the doubling keeps going) stops at the last h that
        succeeded, with the cap reason, the meta entry and a WARNING."""
        from hwave.initial import seed
        f = self.chain(4, (2, 1, 1), ncond=4, T=0.1, two_sz=2)
        spec = _spec([_e(cz=1.0)], ("1/2", "0", "0"), (2, 1, 1))
        real = seed.m_value
        with mock.patch.object(seed, "m_value", lambda *a, **k: 1e-3 * real(*a, **k)), \
                self.assertLogs("qlms.initial.seed", level="WARNING") as cm:
            res, _ = self.seed(f, spec)
        wf = seed.w_floor(f.problem.ham_trans_k)
        self.assertTrue(res.tried)
        self.assertTrue(res.cap_reason.startswith("mu search failed at h="))
        failed = res.meta["mu_search_failed_at"]
        self.assertEqual(res.h_used, res.h_hi)
        self.assertEqual(failed, 2 * res.h_used)
        self.assertGreater(failed, 2 ** 20 * wf)
        self.assertLess(res.amplitude_used, 0.5)
        self.assertTrue(any("mu search failed" in line for line in cm.output))
        np.testing.assert_allclose(np.trace(res.G[0].reshape(4, 4)).real * 2, 4, atol=1e-9)
        # the solver's error itself is still raised by density_at
        from hwave.initial import patterns
        Ok = patterns.realise(spec, f, (2, 1, 1))
        with self.assertRaises(seed.MuSearchFailed):
            seed.density_at(f, f.problem.ham_trans_k, Ok, failed, T_run=0.1, W0=4.0)

    def test_other_runtime_errors_propagate(self):
        """Only the mu search failure is converted: a RuntimeError subclass
        from the fill (here NotImplementedError, RecursionError) propagates."""
        from hwave.initial import patterns, seed
        from hwave.solver import initial_facade
        f = self.chain(4, (2, 1, 1), ncond=4, T=0.1)
        O = patterns.realise(_spec([_e(cz=1.0)], ("1/2", "0", "0"), (2, 1, 1)), f, (2, 1, 1))
        for exc in (NotImplementedError("x"), RecursionError("y")):
            with mock.patch.object(initial_facade, "green_function", side_effect=exc):
                with self.assertRaises(type(exc)):
                    seed.density_at(f, f.problem.ham_trans_k, O, 1.0, T_run=0.1, W0=4.0)
        for exc in (RuntimeError("newton"), initial_facade.MuNotConvergedError("mu")):
            with mock.patch.object(initial_facade, "green_function", side_effect=exc):
                with self.assertRaises(seed.MuSearchFailed):
                    seed.density_at(f, f.problem.ham_trans_k, O, 1.0, T_run=0.1, W0=4.0)

    def test_small_h_unrefined_and_raw(self):
        """At moderate h the fill is the solver's, bit for bit, with no
        refinement."""
        from hwave.initial import patterns, seed
        from hwave.solver import initial_facade
        f = self.chain(4, (2, 1, 1), ncond=4, T=0.1, two_sz=2)
        p = f.problem
        O = patterns.realise(_spec([_e(cz=1.0)], ("1/2", "0", "0"), (2, 1, 1)), f, (2, 1, 1))
        for h in (0.0, 0.3, 3.0, 50.0):
            G, info, meta = seed.density_at(f, p.ham_trans_k, O, h, T_run=0.1, W0=4.0,
                                            return_meta=True)
            self.assertIs(meta["mu_refined"], False)
            tol = 1e-10 * max(1.0, 4.0, h * seed.op_norm(O))
            G_raw, info_raw = initial_facade.green_function(
                p.ham_trans_k - h * O, 0.1, p.blocks, shape=p.shape, ns=p.ns, norb=p.norb,
                fill="uniform", shell_tol=tol, ene_cutoff=p.ene_cutoff, return_info=True)
            np.testing.assert_array_equal(G, G_raw)
            self.assertEqual(info.mu, info_raw.mu)

    def test_mu_failure_in_bisection(self):
        """A mu search failure at a bisection midpoint: the kept point
        already meets the target, so there is no cap_reason; the failure is
        in meta only."""
        from hwave.initial import seed
        f = self.chain(4, (2, 1, 1), ncond=4, T=0.1)
        spec = _spec([_e(c0=1.0)], ("1/2", "0", "0"), (2, 1, 1))
        ref, _ = self.seed(f, spec)
        self.assertIsNone(ref.cap_reason)
        self.assertNotEqual(ref.h_used, ref.h_hi)          # the bisection ran
        # the first doubling point meets the target, so every 0 < h < h_hi
        # call below is a bisection midpoint
        self.assertEqual(ref.h_hi, seed.w_floor(f.problem.ham_trans_k))
        real = seed.density_at
        calls = []

        def failing(*a, **k):
            h = a[3]
            if 0 < h < ref.h_hi:
                calls.append(h)
                if len(calls) == 1:
                    raise seed.MuSearchFailed(h, RuntimeError("test"))
            return real(*a, **k)

        with mock.patch.object(seed, "density_at", failing):
            res, _ = self.seed(f, spec)
        self.assertTrue(res.tried)
        self.assertIsNone(res.cap_reason)
        self.assertEqual(res.meta["mu_search_failed_at"], calls[0])
        self.assertEqual(res.h_used, res.h_hi)
        self.assertGreaterEqual(res.m_used, res.m_target)

    def test_electron_count_warning(self):
        """A count still off after the refinement is logged with its size."""
        from hwave.initial import patterns, seed
        f = self.chain(4, (2, 1, 1), ncond=4, T=0.1)
        p = f.problem
        O = patterns.realise(_spec([_e(cz=1.0)], ("1/2", "0", "0"), (2, 1, 1)), f, (2, 1, 1))
        with mock.patch.object(seed, "N_RTOL", -1.0), \
                self.assertLogs("qlms.initial.seed", level="WARNING") as cm:
            _, _, meta = seed.density_at(f, p.ham_trans_k, O, 1.0, T_run=0.1, W0=4.0,
                                         return_meta=True)
        self.assertTrue(meta["mu_refined"])
        self.assertIn("electron count", cm.output[0])

    def test_amplitude_range(self):
        f = self.chain(4, (2, 1, 1), ncond=4, T=0.1)
        spec = _spec([_e(c0=1.0)], ("1/2", "0", "0"), (2, 1, 1))
        for amp in (1.0, 1.5, 0.0, -0.1):
            with self.assertRaisesRegex(ValueError, "amplitude"):
                self.seed(f, spec, amplitude=amp)
        with self.assertRaisesRegex(ValueError, "seed_base"):
            self.seed(f, spec, seed_base="para")


# --------------------------------------------------------------------------
# the seeded density
# --------------------------------------------------------------------------

class TestDensity(_Tmp):
    def check_density(self, f, res):
        p = f.problem
        G = res.G
        n = p.nvol * float(np.real(np.trace(G[0].reshape(p.nd, p.nd))))
        self.assertAlmostEqual(n, p.ncond, delta=1e-9)
        if p.two_sz is not None:
            nup = p.nvol * float(np.real(np.trace(G[0, 0, :, 0, :])))
            ndn = p.nvol * float(np.real(np.trace(G[0, 1, :, 1, :])))
            self.assertAlmostEqual(nup - ndn, p.two_sz, delta=1e-9)
        Gk = _green_k(G, p.shape)
        np.testing.assert_allclose(Gk, np.conj(np.swapaxes(Gk, -1, -2)), atol=1e-12)
        ev = np.linalg.eigvalsh(Gk)
        self.assertGreaterEqual(float(ev.min()), -1e-10)
        self.assertLessEqual(float(ev.max()), 1 + 1e-10)

    def test_ncond_2sz_exact(self):
        for T in (0.0, 0.1):
            f = self.chain(4, (2, 1, 1), ncond=4, T=T, two_sz=2)
            res, _ = self.seed(f, _spec([_e(cz=1.0)], ("1/2", "0", "0"), (2, 1, 1)))
            self.assertTrue(res.tried)
            self.check_density(f, res)

    def test_eigen_bounds(self):
        for T in (0.0, 0.1):
            f = self.chain(8, (2, 1, 1), ncond=6, T=T)
            for spec in (_spec([_e(cz=1.0)], ("1/2", "0", "0"), (2, 1, 1)),
                         _spec([_e(cx=1.0), _e(r=(1, 0, 0), c0=0.3j)], ("0", "0", "0"),
                               (1, 1, 1))):
                res, _ = self.seed(f, spec)
                self.assertTrue(res.tried)
                self.check_density(f, res)


# --------------------------------------------------------------------------
# constraint modes and block refusals
# --------------------------------------------------------------------------

class TestRefusals(_Tmp):
    def test_constraint_modes_marked(self):
        from hwave.initial import patterns, seed
        q0 = ("0", "0", "0")
        f = self.chain(4, (2, 1, 1), ncond=4, T=0.1)
        charge = patterns.realise(_spec([_e(c0=1.0)], q0, (1, 1, 1)), f, (2, 1, 1))
        spin = patterns.realise(_spec([_e(cz=1.0)], q0, (1, 1, 1)), f, (2, 1, 1))
        self.assertTrue(seed.constraint_mode(charge, f))
        self.assertFalse(seed.constraint_mode(spin, f))
        res, _ = self.seed(f, _spec([_e(c0=1.0)], q0, (1, 1, 1)))
        self.assertFalse(res.tried)
        self.assertEqual(res.not_tried_reason, "no response (constraint)")
        self.assertAlmostEqual(res.m0, 1 / math.sqrt(2), delta=1e-10)   # T > 0 mu xtol
        # 2Sz: uniform S^z is a constraint mode too (one mu per spin)
        f2 = self.chain(4, (2, 1, 1), ncond=4, T=0.1, two_sz=0)
        spin2 = patterns.realise(_spec([_e(cz=1.0)], q0, (1, 1, 1)), f2, (2, 1, 1))
        self.assertTrue(seed.constraint_mode(spin2, f2))
        neel = patterns.realise(_spec([_e(cz=1.0)], ("1/2", "0", "0"), (2, 1, 1)), f2,
                                (2, 1, 1))
        self.assertFalse(seed.constraint_mode(neel, f2))
        # a uniform charge on the interacting orbital only, with a ligand,
        # is not the electron number: it can respond
        t = {((1, 0, 0), 0, 0): -1.0, ((-1, 0, 0), 0, 0): -1.0,
             ((0, 0, 0), 0, 1): -0.5, ((0, 0, 0), 1, 0): -0.5}
        f3 = self.chain(4, (1, 1, 1), ncond=4, T=0.1, norb=2, transfer=t,
                        centers=[(0, 0, 0), (0.5, 0, 0)],
                        intra={((0, 0, 0), 0, 0): 2.0})
        c3 = patterns.realise(_spec([_e(c0=1.0)], q0, (1, 1, 1)), f3, (1, 1, 1))
        self.assertFalse(seed.constraint_mode(c3, f3))

    def test_block_refusals(self):
        from hwave.initial import patterns, seed
        q0 = ("0", "0", "0")
        # transverse spin in 2Sz mode
        f = self.chain(4, (1, 1, 1), ncond=4, T=0.1, two_sz=0)
        O = patterns.realise(_spec([_e(cx=1.0)], q0, (1, 1, 1)), f, (1, 1, 1))
        msgs = seed.block_check(O, f.problem.blocks, ns=2, norb=f.problem.norb)
        self.assertEqual(len(msgs), 1)
        self.assertIn("transverse", msgs[0])
        res, _ = self.seed(f, _spec([_e(cx=1.0)], q0, (1, 1, 1)))
        self.assertFalse(res.tried)
        self.assertTrue(res.not_tried_reason.startswith("refused: block check"))
        self.assertIn("transverse", res.not_tried_reason)
        # an out-of-block orbital entry: two orbitals Transfer never connects
        t = {((1, 0, 0), a, a): -1.0 for a in range(2)}
        t.update({((-1, 0, 0), a, a): -1.0 for a in range(2)})
        f2 = self.chain(4, (1, 1, 1), ncond=4, T=0.1, norb=2, transfer=t)
        O2 = patterns.realise(_spec([_e(0, 1, c0=1.0)], q0, (1, 1, 1)), f2, (1, 1, 1))
        msgs = seed.block_check(O2, f2.problem.blocks, ns=2, norb=f2.problem.norb)
        self.assertEqual(len(msgs), 1)
        self.assertIn("outside the solver blocks", msgs[0])
        self.assertNotIn("transverse", msgs[0])
        # noncollinear CoulombIntra seed: accepted with Fock on, refused off
        spec = _spec([_e(cx=1.0)], q0, (1, 1, 1))
        on = self.chain(4, (1, 1, 1), ncond=4, T=0.1)
        self.assertEqual(seed.block_check(patterns.realise(spec, on, (1, 1, 1)),
                                          on.problem.blocks), [])
        res, _ = self.seed(on, spec)
        self.assertTrue(res.tried)
        off = self.chain(4, (1, 1, 1), ncond=4, T=0.1, fock=False)
        self.assertNotEqual(seed.block_check(patterns.realise(spec, off, (1, 1, 1)),
                                             off.problem.blocks), [])
        res, _ = self.seed(off, spec)
        self.assertFalse(res.tried)
        self.assertTrue(res.not_tried_reason.startswith("refused: block check"))


# --------------------------------------------------------------------------
# the null table and the moments
# --------------------------------------------------------------------------

_FIELDS = ("G", "seed_base", "shell", "moments_h0", "m0", "m_inf", "m_target", "m_used",
           "h_hi", "h_used", "amplitude_requested", "amplitude_used", "cap_reason")


class TestNullTable(_Tmp):
    def assert_set(self, res, names):
        for name in _FIELDS:
            v = getattr(res, name)
            if name in names:
                self.assertIsNotNone(v, name)
            else:
                self.assertIsNone(v, name)

    def test_null_table(self):
        from hwave.initial import seed
        f = self.chain(4, (2, 1, 1), ncond=4, T=0.1)
        stag = _spec([_e(c0=1.0)], ("1/2", "0", "0"), (2, 1, 1))
        all_but_cap = set(_FIELDS) - {"cap_reason"}
        # seeded
        res, _ = self.seed(f, stag)
        self.assertTrue(res.tried)
        self.assertIsNone(res.not_tried_reason)
        self.assert_set(res, all_but_cap)
        self.assertEqual(sorted(res.shell), ["h0", "h_used"])
        # seeded at the cap
        with mock.patch.object(seed, "CAP_DOUBLINGS", 0):
            res, _ = self.seed(f, stag, amplitude=0.99)
        self.assertTrue(res.tried)
        self.assert_set(res, set(_FIELDS))
        self.assertEqual(res.cap_reason, "target not reached at the cap")
        self.assertEqual(res.h_used, res.h_hi)
        self.assertEqual(res.h_used, seed.w_floor(f.problem.ham_trans_k))
        self.assertLess(res.amplitude_used, 0.99)
        # para
        res, _ = self.seed(f, None)
        self.assertTrue(res.tried)
        self.assert_set(res, {"G", "seed_base", "shell", "moments_h0", "h_used"})
        self.assertEqual(res.h_used, 0.0)
        # no response (degenerate with para)
        full = self.chain(4, (2, 1, 1), ncond=8, T=0.1)
        res, _ = self.seed(full, stag)
        self.assertFalse(res.tried)
        self.assertEqual(res.not_tried_reason, "no response (degenerate with para)")
        self.assert_set(res, {"seed_base", "shell", "moments_h0", "m0", "m_inf"})
        # no response (constraint)
        res, _ = self.seed(f, _spec([_e(c0=1.0)], ("0", "0", "0"), (1, 1, 1)))
        self.assertEqual(res.not_tried_reason, "no response (constraint)")
        self.assert_set(res, {"seed_base", "m0"})
        # refused (block check)
        t2 = self.chain(4, (1, 1, 1), ncond=4, T=0.1, two_sz=0)
        res, _ = self.seed(t2, _spec([_e(cx=1.0)], ("0", "0", "0"), (1, 1, 1)))
        self.assertFalse(res.tried)
        self.assert_set(res, set())
        # seeding not performed (quota, unrepresentable, duplicate, null)
        res = seed.SeedResult.not_performed("quota: baseline")
        self.assertFalse(res.tried)
        self.assertEqual(res.not_tried_reason, "quota: baseline")
        self.assert_set(res, set())


class TestMoments(_Tmp):
    def test_moments_h0_use_green_slot_for_pr1a(self):
        """Review Focus 1: PR 1a's order_parameters_from_density reads the
        stored Green slot ``rho[o, s, t] = G[0, s, o, t, o]``. An H_base with
        an on-site ``-b sigma^y`` term gives S^y > 0; feeding rho = G^T
        instead would flip it."""
        from hwave.results import order_parameters_from_density
        entries = {((1, 0, 0), 0, 0, s, s): -1.0 for s in range(2)}
        entries[((0, 0, 0), 0, 0, 0, 1)] = 0.3j          # -0.3 sigma^y
        text = initial_fixtures.spin_transfer_text(1, entries)
        f = self.chain(4, (1, 1, 1), ncond=4, T=0.1, extra={"transfer.dat": text})
        res, _ = self.seed(f, None)
        G0 = res.G
        spin = np.array(res.moments_h0["spin"])
        self.assertEqual(spin.shape, (1, 3))
        self.assertGreater(spin[0, 1], 0.01)
        want = order_parameters_from_density(np.einsum("sata->ast", G0[0]),
                                             site_of_orbital=[0], n_cell=None)
        np.testing.assert_allclose(spin, want.spin, atol=1e-15)
        np.testing.assert_allclose(res.moments_h0["charge"], want.charge, atol=1e-15)
        wrong = order_parameters_from_density(np.einsum("sata->ats", G0[0]),
                                              site_of_orbital=[0], n_cell=None)
        self.assertLess(wrong.spin[0, 1], -0.01)

    def test_moments_per_site_group(self):
        from hwave.initial import seed
        f = self.chain(4, (2, 1, 1), ncond=4, T=0.1)
        res, _ = self.seed(f, _spec([_e(cz=1.0)], ("1/2", "0", "0"), (2, 1, 1)))
        m = seed.moments_of(res.G, f)
        self.assertEqual(len(m["spin"]), 2)            # one per site of the folded cell
        self.assertGreater(m["spin"][0][2], 0.05)
        self.assertAlmostEqual(m["spin"][0][2], -m["spin"][1][2], delta=1e-12)
        np.testing.assert_allclose(res.moments_h0["spin"], 0, atol=1e-12)
        np.testing.assert_allclose(res.moments_h0["charge"], [1, 1], atol=1e-12)

    def test_moments_h0_spin_orbital(self):
        """SO mode: the moments come from PR 1a's so_density_from_green (the
        normal-mode slot does not apply to ns = 1). An on-site Zeeman term
        on physical orbital 1 gives S^z > 0 there at h = 0, and a cz seed on
        physical orbital 1 polarises that orbital."""
        from hwave.initial import seed
        from hwave.results import order_parameters_from_density, so_density_from_green
        t = {}
        for a in range(2):
            for s in range(2):
                t[((1, 0, 0), 2 * a + s, 2 * a + s)] = -1.0
                t[((-1, 0, 0), 2 * a + s, 2 * a + s)] = -1.0
        t[((0, 0, 0), 2, 2)] = -0.2                   # orbital 1, spin up
        t[((0, 0, 0), 3, 3)] = 0.2                    # orbital 1, spin down
        f = self.chain(4, (1, 1, 1), ncond=8, T=0.1, norb=4, transfer=t, so=True,
                       centers=[(0, 0, 0), (0, 0, 0), (0.5, 0, 0), (0.5, 0, 0)])
        self.assertEqual(f.site_groups.groups, ((0,), (1,)))
        self.assertTrue(f.spin_orbital)
        para, _ = self.seed(f, None)
        rho = so_density_from_green(para.G[0, 0, :, 0, :], 2, 1)
        want = order_parameters_from_density(rho, site_of_orbital=[0, 1], n_cell=None)
        np.testing.assert_allclose(para.moments_h0["spin"], want.spin, atol=1e-15)
        self.assertAlmostEqual(para.moments_h0["spin"][0][2], 0.0, delta=1e-12)
        self.assertGreater(para.moments_h0["spin"][1][2], 0.01)
        res, _ = self.seed(f, _spec([_e(1, cz=1.0)], ("0", "0", "0"), (1, 1, 1)))
        self.assertTrue(res.tried)
        m = seed.moments_of(res.G, f)
        self.assertGreater(m["spin"][1][2], para.moments_h0["spin"][1][2] + 0.05)
        self.assertAlmostEqual(m["spin"][0][2], 0.0, delta=1e-10)


if __name__ == "__main__":
    unittest.main()
