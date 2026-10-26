"""The conditioning guard of the bond dressing skips the exact host SVD when
the interval alone already decides the pass (issue #276).

The interval guard (issue #197) runs on numpy as well as on cupy, so the
skip is exercised here on the host: ``guard_method = "interval"`` is passed
explicitly to :func:`dress_batch`, and the end-to-end FLEX runs map the
numpy entry of :data:`bond_channels._GUARD_AUTO` to ``"interval"`` (the
resolution ``"auto"`` makes on a device), so no GPU is needed. The exact
host decompositions are counted with a spy on
:func:`bond_channels.bond_conditioning_score`, which every exact path calls.
"""
import re
import tempfile
import unittest
from unittest import mock

import numpy as np

from hwave.solver import bond_channels as bc
from tests.test_flex_bond_dressing import _onsite_factors, _problem


def _spy():
    return mock.patch.object(bc, "bond_conditioning_score", wraps=bc.bond_conditioning_score)


def _solve_matrix(cb, W, channel):
    mat = cb @ W[None]
    if channel == "spin":
        mat = -mat
    return (mat + np.eye(mat.shape[-1])).reshape(-1, mat.shape[-1], mat.shape[-1])


def _exact_min(flat):
    """The exact score minimum of a flattened stack, from the full SVD."""
    return bc.bond_conditioning_score(flat.reshape(-1, 1, 1, flat.shape[-1], flat.shape[-1]))


def _run(cb, W, channel, exact, **kw):
    """Outcome of one dress_batch call under the interval guard, with the
    number of exact host decompositions it ran."""
    v = []
    with _spy() as spy:
        try:
            with mock.patch.object(bc.logger, "warning") as warn:
                chi, cond = bc.dress_batch(cb, W, channel, guard_method="interval",
                                           guard_exact=exact, violations=v, **kw)
            out = ("ok", chi, cond, v, [c.args for c in warn.call_args_list])
        except bc.BondConditioningError as e:
            out = ("refused", str(e), e.channel, e.l, e.iq, e.q, e.worst, e.ratio, e.pole,
                   e.smin, e.smax)
    return out, spy.call_count


class TestDecisionHelpers(unittest.TestCase):

    def _iv(self, low, up, valid=None):
        low = np.asarray(low, float)
        up = np.asarray(up, float)
        valid = np.ones(len(low), bool) if valid is None else np.asarray(valid, bool)
        return bc.GuardInterval(sigma_max_low=np.ones(len(low)), sigma_max_up=np.ones(len(low)),
                                sigma_min_low=low, sigma_min_up=up, rho=np.zeros(len(low)),
                                valid=valid)

    def test_threshold(self):
        self.assertEqual(bc._GUARD_SKIP_MARGIN, 2.0)
        self.assertEqual(bc.interval_skip_threshold(1e-3, 6),
                         2.0e-3 + bc._GUARD_ROUNDING * 36 * 2.0 ** -53)

    def test_decides_only_above_the_threshold(self):
        thr = bc.interval_skip_threshold(1e-3, 4)
        br = bc.interval_decides_pass(self._iv([0.5, 0.3, 0.4], [0.9, 0.35, 0.33]), 1e-3, 4)
        self.assertEqual(br, bc.GuardBracket(low=0.3, up=0.33, index=2, n_blocks=3))
        # straddling: the smallest lower end at or below the threshold
        self.assertIsNone(bc.interval_decides_pass(self._iv([0.5, thr], [0.9, 0.1]), 1e-3, 4))
        self.assertIsNone(bc.interval_decides_pass(self._iv([0.5, 1.5e-3], [0.9, 0.1]), 1e-3, 4))
        self.assertIsNotNone(bc.interval_decides_pass(
            self._iv([0.5, np.nextafter(thr, 1.0)], [0.9, 0.1]), 1e-3, 4))

    def test_any_invalid_block_or_empty_stack_defers_to_the_exact_path(self):
        self.assertIsNone(bc.interval_decides_pass(
            self._iv([0.5, 0.5], [0.9, 0.9], valid=[True, False]), 1e-3, 4))
        self.assertIsNone(bc.interval_decides_pass(self._iv([], []), 1e-3, 4))

    def test_unknown_modes_are_refused(self):
        cb, S, _ = _problem(nmat=2, nvol=2)
        with self.assertRaises(ValueError):
            bc.resolve_conditioning_guard(_solve_matrix(cb, S, "spin"), np, 1e-3, exact="sometimes")
        with self.assertRaises(ValueError) as cm:
            bc.dress_batch(cb, S, "spin", l0=0, nmat=2, spatial_shape=(2, 1, 1),
                           guard_exact="sometimes")
        self.assertIn("guard_exact", str(cm.exception))


class TestResolverSkip(unittest.TestCase):

    def test_well_conditioned_stack_is_not_decomposed_and_the_bracket_holds(self):
        cb, S, _ = _problem()
        flat = _solve_matrix(cb, S, "spin")
        ref = _exact_min(flat)
        with _spy() as spy:
            out = bc.resolve_conditioning_guard(flat, np, 1e-3, exact="decide")
        self.assertEqual(spy.call_count, 0)
        self.assertIsInstance(out, bc.GuardBracket)
        self.assertLessEqual(out.low, ref[0])
        self.assertLessEqual(ref[0], out.up)
        self.assertGreater(out.low, 1e-3)
        iv = bc.bond_conditioning_interval(flat, np)
        self.assertEqual(out.index, int(np.argmin(iv.score_up)))
        self.assertEqual(out.n_blocks, flat.shape[0])

    def test_straddling_stack_runs_the_exact_path_unchanged(self):
        cb, S, _ = _problem()
        flat = _solve_matrix(cb, S, "spin")
        iv = bc.bond_conditioning_interval(flat, np)
        tol = 0.75 * float(np.min(iv.score_low))      # 2 * tol above the smallest lower end
        with _spy() as spy:
            out = bc.resolve_conditioning_guard(flat, np, tol, exact="decide")
        self.assertEqual(spy.call_count, 1)
        self.assertEqual(out, bc.resolve_conditioning_guard(flat, np, tol))

    def test_always_never_skips(self):
        cb, S, _ = _problem()
        flat = _solve_matrix(cb, S, "spin")
        with _spy() as spy:
            out = bc.resolve_conditioning_guard(flat, np, 1e-3, exact="always")
        self.assertEqual(spy.call_count, 1)
        self.assertEqual(out[:6], _exact_min(flat))


class TestDressBatchDecisions(unittest.TestCase):
    """``guard_exact = "decide"`` against ``"always"`` on the dressing
    fixture: the same pass / refuse / warn outcome and outputs everywhere,
    the exact value whenever the exact path runs, and a bracket that holds
    the exact minimum whenever it is skipped."""

    def _compare(self, cb, W, channel, cond_tol, **kw):
        a, n_a = _run(cb, W, channel, "always", cond_tol=cond_tol, **kw)
        d, n_d = _run(cb, W, channel, "decide", cond_tol=cond_tol, **kw)
        self.assertGreaterEqual(n_a, 1)
        self.assertEqual(a[0], d[0])
        if a[0] == "refused":
            self.assertEqual(a[1:], d[1:])
            self.assertEqual(n_d, 1)
            return "refused"
        np.testing.assert_array_equal(a[1], d[1])
        self.assertEqual(a[3], d[3])                  # warn records
        self.assertEqual(a[4], d[4])                  # warn log lines
        self.assertIsInstance(a[2], float)
        if isinstance(d[2], bc.GuardBracket):
            self.assertEqual(n_d, 0)
            self.assertEqual(a[3], [])
            self.assertLessEqual(d[2].low, a[2])
            self.assertLessEqual(a[2], d[2].up)
            self.assertGreater(d[2].low, cond_tol)
            return "skipped"
        self.assertEqual(n_d, 1)
        self.assertEqual(a[2], d[2])                  # exact float equality
        return "warned" if a[3] else "exact"

    def test_tolerance_sweep_both_channels_and_policies(self):
        chi_bar, S, C = _problem()
        nmat, nvol = chi_bar.shape[:2]
        seen = set()
        for scale in (1.0, 4.0, 8.0, 16.0):
            for channel, W0 in (("spin", S), ("charge", C)):
                W = scale * W0
                exact = _exact_min(_solve_matrix(chi_bar, W, channel))[0]
                for factor in (0.01, 0.1, 0.4, 0.9, 1.0, 1.5):
                    tol = min(0.5, factor * exact)
                    for policy in ("refuse", "warn"):
                        with self.subTest(scale=scale, channel=channel, factor=factor,
                                          policy=policy):
                            seen.add(self._compare(chi_bar, W, channel, tol, l0=0, nmat=nmat,
                                                   spatial_shape=(nvol, 1, 1),
                                                   guard_policy=policy))
        # every branch of the decision was exercised
        self.assertEqual(seen, {"skipped", "exact", "refused", "warned"})

    def test_frequency_sub_batches_report_grid_indices(self):
        chi_bar, S, _ = _problem()
        nmat, nvol = chi_bar.shape[:2]
        for l0, nb in ((0, nmat), (1, 3), (5, 1)):
            with self.subTest(l0=l0, nb=nb):
                cb = chi_bar[l0:l0 + nb]
                _chi, br = bc.dress_batch(cb, S, "spin", l0=l0, nmat=nmat,
                                          spatial_shape=(nvol, 1, 1), guard_method="interval",
                                          guard_exact="decide")
                self.assertIsInstance(br, bc.GuardBracket)
                iv = bc.bond_conditioning_interval(_solve_matrix(cb, S, "spin"), np)
                self.assertEqual(br.index, l0 * nvol + int(np.argmin(iv.score_up)))

    def test_skip_only_under_guard_freqs_all_and_the_interval_method(self):
        chi_bar, S, _ = _problem()
        nmat, nvol = chi_bar.shape[:2]
        kw = dict(l0=0, nmat=nmat, spatial_shape=(nvol, 1, 1), guard_exact="decide")
        for method, freqs in (("svd", "all"), ("auto", "all"), ("interval", "static"),
                              ("svd", "static")):
            with self.subTest(method=method, freqs=freqs):
                ref_chi, ref = bc.dress_batch(chi_bar, S, "spin", guard_method=method,
                                              guard_freqs=freqs, l0=0, nmat=nmat,
                                              spatial_shape=(nvol, 1, 1))
                with _spy() as spy:
                    chi, cond = bc.dress_batch(chi_bar, S, "spin", guard_method=method,
                                               guard_freqs=freqs, **kw)
                self.assertEqual(spy.call_count, 1)
                self.assertIsInstance(cond, float)
                self.assertEqual(cond, ref)
                np.testing.assert_array_equal(chi, ref_chi)

    def test_counts(self):
        chi_bar, S, _ = _problem()
        nmat, nvol = chi_bar.shape[:2]
        stats = {"guard_exact_blocks": 0, "guard_blocks": 0}
        bc.dress_batch(chi_bar, S, "spin", l0=0, nmat=nmat, spatial_shape=(nvol, 1, 1),
                       guard_method="interval", guard_exact="decide", stats=stats)
        self.assertEqual(stats, {"guard_exact_blocks": 0, "guard_blocks": nmat * nvol})

    def test_default_is_always(self):
        chi_bar, S, _ = _problem()
        nmat, nvol = chi_bar.shape[:2]
        with _spy() as spy:
            _chi, cond = bc.dress_batch(chi_bar, S, "spin", l0=0, nmat=nmat,
                                        spatial_shape=(nvol, 1, 1), guard_method="interval")
        self.assertEqual(spy.call_count, 1)
        self.assertIsInstance(cond, float)


class TestDressAndBuildW(unittest.TestCase):

    def _map(self, exact, chi_bar, S, C, nb=3, cond_tol=bc._BOND_COND_FLOOR):
        import types
        from hwave.solver import flex_bond as fb
        factors = _onsite_factors()
        nmat, nvol, ND, nd = chi_bar.shape[0], chi_bar.shape[1], S.shape[-1], 4
        view = types.SimpleNamespace(n_channels=ND // nd)
        with fb.BondBlockStore(nmat, nvol, ND, nd, ("chibar", "W")) as store, \
                fb.BondDeviceContext.for_view(np, S, C, view, 2) as dev:
            store.put_freq_batch("chibar", 0, nmat, chi_bar)
            with _spy() as spy:
                res = fb.dress_and_build_w(store, dev, nb=nb, output_full=False, nmat=nmat,
                                           nvol=nvol, nd=nd, spatial_shape=(nvol, 1, 1),
                                           factors=factors, guard_method="interval",
                                           guard_exact=exact, cond_tol=cond_tol)
            return store.get_freq_batch("W", 0, nmat).copy(), res, spy.call_count

    def test_map_outputs_identical_and_bracket_holds_the_exact_minimum(self):
        chi_bar, S, C = _problem()
        W_a, r_a, n_a = self._map("always", chi_bar, S, C)
        W_d, r_d, n_d = self._map("decide", chi_bar, S, C)
        np.testing.assert_array_equal(W_d, W_a)
        for f in ("collapse0", "collapse_s", "collapse_c", "static_s", "static_c"):
            np.testing.assert_array_equal(getattr(r_d, f), getattr(r_a, f))
        self.assertEqual(n_a, 4)                       # two batches x two channels
        self.assertEqual(n_d, 0)
        self.assertIsNone(r_a.cond_bracket_s)
        self.assertIsNone(r_a.cond_bracket_c)
        self.assertEqual(r_d.guard_exact_blocks, 0)
        self.assertEqual(r_d.guard_blocks, r_a.guard_blocks)
        nvol = chi_bar.shape[1]
        for ch, W in (("s", S), ("c", C)):
            exact = getattr(r_a, "cond_min_" + ch)
            low, up, index = getattr(r_d, "cond_bracket_" + ch)
            self.assertTrue(np.isnan(getattr(r_d, "cond_min_" + ch)))
            self.assertLessEqual(low, exact)
            self.assertLessEqual(exact, up)
            # the smallest upper end over the whole grid, at its grid index
            iv = bc.bond_conditioning_interval(
                _solve_matrix(chi_bar, W, "spin" if ch == "s" else "charge"), np)
            self.assertEqual(up, float(np.min(iv.score_up)))
            self.assertEqual(index, int(np.argmin(iv.score_up)))
            self.assertLess(index, chi_bar.shape[0] * nvol)

    def test_mixed_map_reports_the_exact_minimum(self):
        # one frequency batch near the tolerance (exact path), the others
        # well conditioned (skipped): the minimum is the decomposed batch's
        chi_bar, S, C = _problem()
        chi_bar = chi_bar.copy()
        chi_bar[:3] *= 0.05                            # batch [0, 3): far from the floor
        exact_hot = _exact_min(_solve_matrix(chi_bar[3:], S, "spin"))[0]
        tol = 0.4 * exact_hot
        _W, r_a, _ = self._map("always", chi_bar, S, C, cond_tol=tol)
        _W, r_d, n_d = self._map("decide", chi_bar, S, C, cond_tol=tol)
        self.assertGreater(n_d, 0)
        self.assertLess(n_d, 4)
        self.assertIsNone(r_d.cond_bracket_s)
        self.assertEqual(r_d.cond_min_s, r_a.cond_min_s)

    def test_cond_min_aggregation(self):
        from hwave.solver.flex_bond import _CondMin
        m = _CondMin()
        m.add(None)
        self.assertEqual(m.result(), (np.inf, None))
        m.add(0.5)
        m.add(0.2)
        self.assertEqual(m.result(), (0.2, None))
        m.add(bc.GuardBracket(low=0.5, up=0.9, index=7, n_blocks=4))
        self.assertEqual(m.result(), (0.2, None))      # 0.2 <= every lower end
        m.add(bc.GuardBracket(low=0.2, up=0.6, index=3, n_blocks=4))
        self.assertEqual(m.result(), (0.2, None))      # a tie is still exact
        m.add(bc.GuardBracket(low=0.1, up=0.6, index=4, n_blocks=4))
        cond, br = m.result()
        self.assertTrue(np.isnan(cond))
        self.assertEqual(br, (0.1, 0.2, None))         # the exact batch holds the upper end
        b = _CondMin()
        b.add(bc.GuardBracket(low=0.4, up=np.inf, index=0, n_blocks=1))
        b.add(bc.GuardBracket(low=0.3, up=0.7, index=5, n_blocks=1))
        b.add(bc.GuardBracket(low=0.6, up=0.8, index=9, n_blocks=1))
        cond, br = b.result()
        self.assertTrue(np.isnan(cond))
        self.assertEqual(br, (0.3, 0.7, 5))


_COND_RE = re.compile(r"bond cond_min spin (.*) charge (.*) \(tol ")
_BRACKET_RE = re.compile(r"in bracket \[(\S+), (\S+)\] \(interval bound, no exact SVD; "
                         r"smallest upper end at (flat index \d+|n/a)\)")


def _parse(text):
    m = _BRACKET_RE.fullmatch(text)
    if m:
        return ("bracket", float(m.group(1)), float(m.group(2)))
    return ("exact", float(text))


class TestFlexEndToEnd(unittest.TestCase):
    """The bond-gate FLEX on the committed two-orbital input, with the
    interval guard the device path resolves to."""

    def _solve(self, extra):
        from tests.test_flex_bond_gate import _flex
        s, r = _flex(extra)
        gi = r.get_param("green")
        with mock.patch.dict(bc._GUARD_AUTO, {"numpy": "interval"}), _spy() as spy, \
                tempfile.TemporaryDirectory() as out, \
                self.assertLogs("hwave.solver.flex", level="INFO") as cm:
            s.solve(gi, out)
        lines = [m.groups() for m in (_COND_RE.search(x) for x in cm.output) if m]
        return s, gi, lines, spy.call_count

    def test_decide_matches_always_and_archives_the_exact_minimum(self):
        s_a, gi_a, lines_a, n_a = self._solve({"IterationMax": 3,
                                               "longitudinal_bond_guard_exact": "always"})
        s_d, gi_d, lines_d, n_d = self._solve({"IterationMax": 3})
        self.assertEqual(s_d.longitudinal_bond_guard_exact, "decide")
        self.assertEqual(len(lines_a), 3)
        self.assertEqual(len(lines_d), 3)
        for k in ("sigma", "green", "chiq_s", "chiq_c", "longitudinal_bond_chi_s",
                  "longitudinal_bond_chi_c", "longitudinal_bond_cond_min_s",
                  "longitudinal_bond_cond_min_c"):
            np.testing.assert_array_equal(np.asarray(gi_d[k]), np.asarray(gi_a[k]), err_msg=k)
        # maps 1 and 2 were passed on the interval; map 3 (IterationMax) is exact
        self.assertLess(n_d, n_a)
        for it, (la, ld) in enumerate(zip(lines_a, lines_d), start=1):
            for ta, td in zip(la, ld):
                ka, va = _parse(ta)
                self.assertEqual(ka, "exact")
                pd = _parse(td)
                if it < 3:
                    self.assertEqual(pd[0], "bracket", (it, td))
                    # the printed bracket holds the printed exact minimum
                    self.assertLessEqual(pd[1], va)
                    self.assertLessEqual(va, pd[2])
                else:
                    self.assertEqual(pd, ("exact", va))
        self.assertEqual(float(gi_d["longitudinal_bond_cond_min_s"]), s_a._bond_last.cond_min_s)
        self.assertIsNone(s_d._bond_last.cond_bracket_s)

    def test_converging_map_is_exact(self):
        # EPS = -3 (eps = 1e3): every map passes the residual test, so the
        # counter completes at map 3 of IterationMax = 50; that map must be exact
        s_a, gi_a, lines_a, _ = self._solve({"IterationMax": 50, "EPS": -3,
                                             "longitudinal_bond_guard_exact": "always"})
        s_d, gi_d, lines_d, _ = self._solve({"IterationMax": 50, "EPS": -3})
        self.assertTrue(s_d.scf_converged)
        self.assertEqual(s_d.scf_iterations, 3)
        self.assertEqual(len(lines_d), 3)
        self.assertEqual([_parse(t)[0] for t in lines_d[0] + lines_d[1]], ["bracket"] * 4)
        self.assertEqual([_parse(t)[0] for t in lines_d[2]], ["exact"] * 2)
        for ch in ("s", "c"):
            k = "longitudinal_bond_cond_min_" + ch
            self.assertEqual(float(gi_d[k]), float(gi_a[k]))

    def test_refusal_is_identical(self):
        # a floor above the fixture's minimum: both modes refuse in map 1
        # with the same structured error
        errs = []
        for exact in ("always", "decide"):
            from tests.test_flex_bond_gate import _flex
            s, r = _flex({"IterationMax": 3, "longitudinal_bond_cond_tol": 0.9,
                          "longitudinal_bond_guard_exact": exact})
            gi = r.get_param("green")
            with mock.patch.dict(bc._GUARD_AUTO, {"numpy": "interval"}), \
                    tempfile.TemporaryDirectory() as out:
                with self.assertRaises(bc.BondConditioningError) as cm:
                    s.solve(gi, out)
            e = cm.exception
            errs.append((str(e), e.channel, e.l, e.iq, e.worst, e.smin, e.smax))
        self.assertEqual(errs[0], errs[1])

    def test_numpy_auto_path_is_unchanged(self):
        # the CPU path resolves "auto" to "svd": "decide" changes nothing
        from tests.test_flex_bond_gate import _flex
        outs = []
        for exact in ("always", "decide"):
            s, r = _flex({"IterationMax": 2, "longitudinal_bond_guard_exact": exact})
            gi = r.get_param("green")
            with _spy() as spy, tempfile.TemporaryDirectory() as out:
                s.solve(gi, out)
            outs.append((gi, spy.call_count, s._bond_guard_exact_total))
        self.assertEqual(outs[0][1], outs[1][1])
        self.assertEqual(outs[0][2], outs[1][2])
        for k in ("sigma", "green", "longitudinal_bond_cond_min_s", "longitudinal_bond_cond_min_c"):
            np.testing.assert_array_equal(np.asarray(outs[1][0][k]), np.asarray(outs[0][0][k]))


class TestMapMode(unittest.TestCase):

    def _solver(self, exact="decide", max_iter=10):
        from hwave.solver.flex import FLEX
        s = FLEX.__new__(FLEX)
        s.longitudinal_bond_guard_exact = exact
        s.max_iter = max_iter
        return s

    def test_forced_on_every_map_that_can_be_the_last(self):
        s = self._solver()
        self.assertEqual(s._bond_guard_exact_for_map(1, 0, 3), "decide")
        self.assertEqual(s._bond_guard_exact_for_map(5, 1, 3), "decide")
        self.assertEqual(s._bond_guard_exact_for_map(5, 2, 3), "always")
        self.assertEqual(s._bond_guard_exact_for_map(10, 0, 3), "always")
        self.assertEqual(self._solver(max_iter=1)._bond_guard_exact_for_map(1, 0, 3), "always")
        a = self._solver("always")
        self.assertEqual(a._bond_guard_exact_for_map(1, 0, 3), "always")


class TestConfigKey(unittest.TestCase):

    def _raw(self, param, gate=True):
        from hwave.solver.flex import FLEX
        base = {"flex_hartree_fock": True, "longitudinal_bond_channels": gate, "Nmat": 8}
        base.update(param)
        return FLEX._parse_phase_b_keys({"param": base, "calc_scheme": "general"})

    def test_default_and_case_insensitive_values_and_key(self):
        self.assertEqual(self._raw({})["longitudinal_bond_guard_exact"], "decide")
        self.assertEqual(self._raw({"longitudinal_bond_guard_exact": " ALWAYS "})
                         ["longitudinal_bond_guard_exact"], "always")
        self.assertEqual(self._raw({"Longitudinal_Bond_Guard_Exact": "Decide"})
                         ["longitudinal_bond_guard_exact"], "decide")
        self.assertEqual(self._raw({"LONGITUDINAL_BOND_GUARD_EXACT": "always"})
                         ["longitudinal_bond_guard_exact"], "always")

    def test_invalid_values_are_refused(self):
        for bad in ("never", "", 1, True, None, ["always"]):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError) as cm:
                    self._raw({"longitudinal_bond_guard_exact": bad})
                self.assertIn("longitudinal_bond_guard_exact", str(cm.exception))

    def test_stale_without_the_gate(self):
        with self.assertLogs("hwave.solver.flex", level="WARNING") as cm:
            raw = self._raw({"longitudinal_bond_guard_exact": "always"}, gate=False)
        self.assertIn("longitudinal_bond_guard_exact", "\n".join(cm.output))
        self.assertEqual(raw["longitudinal_bond_guard_exact"], "decide")


if __name__ == "__main__":
    unittest.main()
