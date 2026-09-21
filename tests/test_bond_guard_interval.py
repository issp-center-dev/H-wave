"""The interval guard of the bond dressing (issue #197): certified
two-sided bounds of every block's conditioning score from batched device
operations, the exact set of blocks that could be the batch minimum, and
the resolver that reproduces the SVD guard's outputs from them."""
import unittest
import warnings
from unittest import mock

import numpy as np

from hwave.solver import bond_channels as bc

ND = 6


def _score_ref(blocks):
    """Reference per-block (sigma_max, sigma_min, score) from the production
    scorer, one block at a time so every block is scored (the scorer itself
    returns only the minimum)."""
    out = []
    for A in blocks:
        worst, _iq, _ratio, pole, smin, smax = bc.bond_conditioning_score(A[None, None, None])
        out.append((smax, smin, worst))
    return np.array(out)


def _near_identity(n, seed, scale=0.3):
    rng = np.random.default_rng(seed)
    G = rng.normal(size=(n, ND, ND)) + 1j * rng.normal(size=(n, ND, ND))
    return np.eye(ND) + scale * G / np.sqrt(ND)


def _with_singular_values(sv, seed=0):
    """A block with the prescribed singular values (random unitary factors)."""
    rng = np.random.default_rng(seed)
    U, _ = np.linalg.qr(rng.normal(size=(ND, ND)) + 1j * rng.normal(size=(ND, ND)))
    V, _ = np.linalg.qr(rng.normal(size=(ND, ND)) + 1j * rng.normal(size=(ND, ND)))
    return U @ np.diag(np.asarray(sv, dtype=float)) @ V.conj().T


def _assert_enclosed(test, iv, blocks, mask=None):
    ref = _score_ref(blocks)
    for j in range(len(blocks)):
        if mask is not None and not mask[j]:
            continue
        with test.subTest(block=j):
            test.assertTrue(iv.valid[j])
            smax, smin, score = ref[j]
            test.assertLessEqual(iv.sigma_max_low[j], smax)
            test.assertLessEqual(smax, iv.sigma_max_up[j])
            test.assertLessEqual(iv.sigma_min_low[j], smin)
            test.assertLessEqual(smin, iv.sigma_min_up[j])
            lo = iv.sigma_min_low[j] / max(1.0, iv.sigma_max_up[j])
            up = iv.sigma_min_up[j] / max(1.0, iv.sigma_max_low[j])
            test.assertLessEqual(lo, score)
            test.assertLessEqual(score, up)


class TestIntervalContainment(unittest.TestCase):

    def test_near_identity_random(self):
        blocks = _near_identity(64, 1)
        iv = bc.bond_conditioning_interval(blocks, np)
        _assert_enclosed(self, iv, blocks)
        self.assertEqual(iv.valid.shape, (64,))
        self.assertTrue(np.all(iv.rho < 1e-10))

    def test_ill_conditioned_scores(self):
        blocks = np.stack([_with_singular_values([1.0] * (ND - 1) + [s], seed=i)
                           for i, s in enumerate([1e-6, 1e-4, 1e-2, 1e-1])])
        _assert_enclosed(self, bc.bond_conditioning_interval(blocks, np), blocks)

    def test_repeated_singular_values(self):
        blocks = np.stack([0.5 * np.eye(ND), 3.0 * np.eye(ND), np.eye(ND)]).astype(complex)
        iv = bc.bond_conditioning_interval(blocks, np)
        _assert_enclosed(self, iv, blocks)
        # Rayleigh ends are exact here (every direction is dominant)
        np.testing.assert_allclose(iv.sigma_max_low, [0.5, 3.0, 1.0], rtol=1e-12)
        np.testing.assert_allclose(iv.sigma_min_up, [0.5, 3.0, 1.0], rtol=1e-12)

    def test_ends_are_rounded_outward(self):
        # both blocks are pure scalar multiples of the identity, so the
        # Rayleigh and norm ends coincide exactly with the reference singular
        # value in exact arithmetic; the outward-rounding term (_GUARD_ROUNDING)
        # must still open the interval strictly on both sides at the ulp level
        blocks = np.stack([3.0 * np.eye(ND), 0.5 * np.eye(ND)]).astype(complex)
        iv = bc.bond_conditioning_interval(blocks, np)
        for j, ref in enumerate([3.0, 0.5]):
            with self.subTest(block=j):
                self.assertLess(iv.sigma_max_low[j], ref)
                self.assertLess(ref, iv.sigma_max_up[j])
                self.assertLess(iv.sigma_min_low[j], ref)
                self.assertLess(ref, iv.sigma_min_up[j])
                widen = iv.sigma_max_up[j] / ref - 1.0
                self.assertGreater(widen, 0.0)
                self.assertLess(widen, 1e-11)

    def test_top_direction_orthogonal_to_the_start_vector(self):
        # diag(1, 1, ..., 1, 5) with the all-ones start: the power iteration
        # sees the top direction (its component is nonzero), and a block whose
        # top eigenvector has a ZERO component along the start (antisymmetric
        # under a permutation) must still be enclosed, only loosely
        d = np.ones(ND); d[-1] = 5.0
        A1 = np.diag(d).astype(complex)
        Q = np.eye(ND); Q[0, 0] = Q[1, 1] = 0; Q[0, 1] = 1; Q[1, 0] = -1   # rotates e0-e1
        A2 = (Q @ np.diag([1.0] + [0.2] * (ND - 1)) @ Q.T).astype(complex)  # top direction (e0 - e1)/sqrt2 ... orthogonal to ones
        blocks = np.stack([A1, A2])
        _assert_enclosed(self, bc.bond_conditioning_interval(blocks, np), blocks)

    def test_scales_inside_and_outside_the_scale_guard(self):
        base = _near_identity(4, 3)
        small, large = base * 1e-8, base * 1e8
        small90, large90 = base * 1e-90, base * 1e90
        inside = np.concatenate([small, large, small90, large90])
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            iv = bc.bond_conditioning_interval(inside, np)
        self.assertEqual(caught, [])
        _assert_enclosed(self, iv, inside)
        # 1e+-90 stays deep inside the scale guard (+-1e100): the Rayleigh
        # ends must stay finite and strictly positive there too, not just
        # near identity scale
        near_ends = slice(len(small) + len(large), len(inside))
        self.assertTrue(np.all(iv.sigma_max_low[near_ends] > 0))
        self.assertTrue(np.all(np.isfinite(iv.sigma_min_up[near_ends])))

        outside = np.concatenate([base * 1e-150, base * 1e150])
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            iv = bc.bond_conditioning_interval(outside, np)
        self.assertEqual(caught, [])
        self.assertFalse(np.any(iv.valid))
        self.assertTrue(np.all(iv.sigma_min_low == 0.0))
        self.assertTrue(np.all(np.isinf(iv.sigma_min_up)))

    def test_non_normal_jordan_block(self):
        J = np.eye(ND) + np.diag(np.full(ND - 1, 10.0), 1)     # strongly non-normal
        blocks = np.stack([J.astype(complex)])
        _assert_enclosed(self, bc.bond_conditioning_interval(blocks, np), blocks)

    def test_inaccurate_inverse_is_enclosed_or_invalid(self):
        # sigma_min ~ 1e-14: the inverse is inaccurate, rho is large; either
        # rho < 1 and the (loose) interval encloses, or the block is invalid
        A = _with_singular_values([1.0] * (ND - 1) + [1e-14], seed=7)
        iv = bc.bond_conditioning_interval(np.stack([A]), np)
        if iv.valid[0]:
            _assert_enclosed(self, iv, np.stack([A]))
        else:
            self.assertEqual(iv.sigma_min_low[0], 0.0)

    def test_exactly_singular_member_inside_a_valid_stack(self):
        blocks = _near_identity(5, 4)
        blocks[2] = 0.0
        iv = bc.bond_conditioning_interval(blocks, np)
        self.assertFalse(iv.valid[2])
        self.assertTrue(np.isnan(iv.rho[2]) or iv.rho[2] >= 1.0)
        mask = np.ones(5, bool); mask[2] = False
        _assert_enclosed(self, iv, blocks, mask=mask)

    def test_all_singular_stack(self):
        blocks = np.zeros((3, ND, ND), complex)
        iv = bc.bond_conditioning_interval(blocks, np)
        self.assertFalse(np.any(iv.valid))

    def test_corrupted_inverse_is_caught_by_the_residual(self):
        blocks = _near_identity(3, 5)
        real_inv = np.linalg.inv
        def bad_inv(a):
            out = real_inv(a)
            out[1] *= 1.1                      # wrong by 10 %: rho = 0.1 * sqrt(ND) ~ 0.245,
                                                # still < 1 -- a valid block whose inverse is off
                                                # by 10 %, still enclosed because the Neumann
                                                # bounds hold for any B with rho < 1
            out[2] += 10.0                     # garbage: rho > 1
            return out
        with mock.patch.object(np.linalg, "inv", bad_inv):
            iv = bc.bond_conditioning_interval(blocks, np)
        self.assertTrue(iv.valid[0]); self.assertTrue(iv.valid[1]); self.assertFalse(iv.valid[2])
        _assert_enclosed(self, iv, blocks, mask=np.array([True, True, False]))

    def test_property_based_against_the_scorer(self):
        for seed in range(20):
            blocks = _near_identity(256, 100 + seed, scale=0.9)
            _assert_enclosed(self, bc.bond_conditioning_interval(blocks, np), blocks)

    def test_refuses_other_modules_and_dtypes(self):
        with self.assertRaises(TypeError):
            bc.bond_conditioning_interval(np.eye(ND, dtype=np.complex64)[None], np)
        with self.assertRaises(TypeError):
            bc.bond_conditioning_interval(np.eye(ND)[None].astype(complex), object())

    def test_oom_types_propagate(self):
        class _Oom(Exception):
            pass
        with mock.patch.object(bc._bk, "oom_error_types", lambda: (_Oom,)), \
                mock.patch.object(np.linalg, "inv", mock.Mock(side_effect=_Oom("device"))):
            with self.assertRaises(_Oom):
                bc.bond_conditioning_interval(_near_identity(2, 6), np)


class TestExactSet(unittest.TestCase):

    def _interval(self, low, up, valid=None):
        low = np.asarray(low, float); up = np.asarray(up, float)
        valid = np.ones(len(low), bool) if valid is None else np.asarray(valid, bool)
        return bc.GuardInterval(sigma_max_low=np.ones(len(low)), sigma_max_up=np.ones(len(low)),
                                sigma_min_low=low, sigma_min_up=up,
                                rho=np.zeros(len(low)), valid=valid)

    def test_candidates_within_the_margin_of_the_smallest_upper_end(self):
        # score = sigma_min here (sigma_max ends are 1)
        iv = self._interval(low=[0.05, 0.30, 0.09, 0.50], up=[0.2, 0.6, 0.1, 0.9])
        # U = 0.1; margin 2 -> low <= 0.2: blocks 0 and 2
        np.testing.assert_array_equal(bc.select_exact_blocks(iv), [0, 2])
        np.testing.assert_array_equal(bc.select_exact_blocks(iv, margin=1.0), [0, 2])
        np.testing.assert_array_equal(bc.select_exact_blocks(iv, margin=4.0), [0, 1, 2])

    def test_invalid_blocks_are_always_included_and_ignored_for_u(self):
        iv = self._interval(low=[0.0, 0.30, 0.09], up=[np.inf, 0.6, 0.1], valid=[False, True, True])
        np.testing.assert_array_equal(bc.select_exact_blocks(iv), [0, 2])

    def test_all_invalid_selects_everything(self):
        iv = self._interval(low=[0, 0], up=[np.inf, np.inf], valid=[False, False])
        np.testing.assert_array_equal(bc.select_exact_blocks(iv), [0, 1])


class TestResolver(unittest.TestCase):

    def _stack(self, scores, seed=0):
        return np.stack([_with_singular_values([1.0] * (ND - 1) + [s], seed=seed + i)
                         for i, s in enumerate(scores)])

    def test_matches_the_full_stack_scorer(self):
        blocks = self._stack([0.5, 0.02, 0.3, 0.8, 0.05])
        ref = bc.bond_conditioning_score(blocks.reshape(5, 1, 1, ND, ND))
        out = bc.resolve_conditioning_guard(blocks, np, 1.0e-3)
        self.assertEqual(out[:6], ref)                 # exact float equality
        worst, i, _r, _p, _smin, _smax, n_exact, n_blocks = out
        self.assertEqual(i, 1)
        self.assertEqual(n_blocks, 5)
        self.assertLess(n_exact, 5)                    # 0.02 is far below 0.05 * ... the rest
        self.assertGreaterEqual(n_exact, 1)

    def test_tied_minima_report_the_first_index(self):
        A = _with_singular_values([1.0] * (ND - 1) + [0.01], seed=3)
        blocks = np.stack([A * 1.0, _with_singular_values([1.0] * ND, seed=4), A.copy()])
        ref = bc.bond_conditioning_score(blocks.reshape(3, 1, 1, ND, ND))
        out = bc.resolve_conditioning_guard(blocks, np, 1.0e-3)
        self.assertEqual(out[:6], ref)
        self.assertEqual(out[1], 0)

    def test_singular_member_is_decomposed_and_located(self):
        blocks = self._stack([0.5, 0.4, 0.3])
        blocks[1] = 0.0
        ref = bc.bond_conditioning_score(blocks.reshape(3, 1, 1, ND, ND))
        out = bc.resolve_conditioning_guard(blocks, np, 1.0e-3)
        self.assertEqual(out[:6], ref)
        self.assertEqual(out[1], 1)
        self.assertGreaterEqual(out[6], 1)

    def test_uniform_batch_decomposes_everything(self):
        blocks = _near_identity(32, 9, scale=0.05)
        out = bc.resolve_conditioning_guard(blocks, np, 1.0e-3)
        ref = bc.bond_conditioning_score(blocks.reshape(32, 1, 1, ND, ND))
        self.assertEqual(out[:6], ref)
        self.assertEqual(out[6], 32)                   # the documented worst case

    def test_property_based_against_the_full_scorer(self):
        for seed in range(20):
            rng = np.random.default_rng(seed)
            scores = 10.0 ** rng.uniform(-5, 0, size=64)
            blocks = self._stack(scores, seed=1000 * seed)
            ref = bc.bond_conditioning_score(blocks.reshape(64, 1, 1, ND, ND))
            out = bc.resolve_conditioning_guard(blocks, np, 1.0e-3)
            self.assertEqual(out[:6], ref, seed)


class TestDressBatchGuardMethod(unittest.TestCase):

    def _batch(self, sigma_min=1.0e-9, l_bad=1, q_bad=1, nb=3, nvol=2):
        cb = np.zeros((nb, nvol, 2, 2), dtype=complex)
        for l in range(nb):
            for q in range(nvol):
                c = 0.3 if (l, q) != (l_bad, q_bad) else 1.0 - sigma_min
                cb[l, q] = c * np.eye(2)
        W = np.broadcast_to(np.eye(2, dtype=complex), (nvol, 2, 2)).copy()
        return cb, W

    def test_rejects_an_unknown_method(self):
        cb, W = self._batch(sigma_min=0.5)
        with self.assertRaises(ValueError) as cm:
            bc.dress_batch(cb, W, "spin", l0=0, nmat=16, spatial_shape=(2, 1, 1), guard_method="fast")
        self.assertIn("guard_method", str(cm.exception))

    def test_auto_is_svd_on_numpy(self):
        self.assertEqual(bc._GUARD_AUTO["numpy"], "svd")
        cb, W = self._batch(sigma_min=0.5)
        with mock.patch.object(bc, "resolve_conditioning_guard",
                               side_effect=AssertionError("interval path must not run")):
            bc.dress_batch(cb, W, "spin", l0=0, nmat=16, spatial_shape=(2, 1, 1))

    def test_interval_reproduces_the_refusal(self):
        cb, W = self._batch()
        outs = []
        for method in ("svd", "interval"):
            with self.assertRaises(bc.BondConditioningError) as cm:
                bc.dress_batch(cb, W, "spin", l0=4, nmat=16, spatial_shape=(2, 1, 1),
                               guard_method=method)
            e = cm.exception
            outs.append((str(e), e.channel, e.l, e.iq, e.q, e.worst, e.ratio, e.pole, e.smin, e.smax))
        self.assertEqual(outs[0], outs[1])

    def test_interval_reproduces_a_passing_batch_and_counts(self):
        cb, W = self._batch(sigma_min=0.5)
        ref, c_ref = bc.dress_batch(cb, W, "charge", l0=0, nmat=16, spatial_shape=(2, 1, 1))
        stats = {"guard_exact_blocks": 0, "guard_blocks": 0}
        out, c_out = bc.dress_batch(cb, W, "charge", l0=0, nmat=16, spatial_shape=(2, 1, 1),
                                    guard_method="interval", stats=stats)
        np.testing.assert_array_equal(out, ref)
        self.assertEqual(c_out, c_ref)
        self.assertEqual(stats["guard_blocks"], 6)
        self.assertLessEqual(stats["guard_exact_blocks"], 6)
        self.assertGreaterEqual(stats["guard_exact_blocks"], 1)

    def test_static_mode_counts_only_the_slice(self):
        cb, W = self._batch(sigma_min=0.5)
        stats = {"guard_exact_blocks": 0, "guard_blocks": 0}
        bc.dress_batch(cb, W, "spin", l0=6, nmat=16, spatial_shape=(2, 1, 1),
                       guard_freqs="static", guard_method="interval", stats=stats)
        self.assertEqual(stats["guard_blocks"], 2)

    def test_warn_records_are_identical(self):
        cb, W = self._batch()
        recs = []
        for method in ("svd", "interval"):
            v = []
            with self.assertLogs("qlms.solver.bond_channels", level="WARNING") as cm:
                bc.dress_batch(cb, W, "spin", l0=4, nmat=16, spatial_shape=(2, 1, 1),
                               guard_policy="warn", violations=v, guard_method=method)
            recs.append((v, cm.output))
        self.assertEqual(recs[0], recs[1])

    def test_static_refusal_with_offset_is_identical(self):
        cb, W = self._batch(l_bad=2, q_bad=1)
        outs = []
        for method in ("svd", "interval"):
            with self.assertRaises(bc.BondConditioningError) as cm:
                bc.dress_batch(cb, W, "spin", l0=6, nmat=16, spatial_shape=(2, 1, 1),
                               guard_freqs="static", guard_method=method)
            e = cm.exception
            cause = e.__cause__
            outs.append((str(e), e.channel, e.l, e.iq, e.q, e.worst, e.ratio, e.pole,
                        e.smin, e.smax, cause.iq, cause.q))
        self.assertEqual(outs[0], outs[1])


class TestDressAndBuildWGuardMethod(unittest.TestCase):

    def test_map_is_identical_and_counts_are_reported(self):
        import types
        from hwave.solver import flex_bond as fb
        from tests.test_flex_bond_dressing import _problem
        chi_bar, S, C = _problem()
        nmat, nvol, ND, nd = chi_bar.shape[0], chi_bar.shape[1], S.shape[-1], 4
        S_on = np.ascontiguousarray(S[:, :nd, :nd]); C_on = np.ascontiguousarray(C[:, :nd, :nd])
        view = types.SimpleNamespace(n_channels=ND // nd)

        def run(method):
            with fb.BondBlockStore(nmat, nvol, ND, nd, ("chibar", "W")) as store, \
                    fb.BondDeviceContext.for_view(np, S, C, S_on, C_on, view, 2) as dev:
                store.put_freq_batch("chibar", 0, nmat, chi_bar)
                res = fb.dress_and_build_w(store, dev, nb=3, output_full=False, nmat=nmat,
                                           nvol=nvol, nd=nd, spatial_shape=(nvol, 1, 1),
                                           second_order="takimoto", guard_method=method)
                return store.get_freq_batch("W", 0, nmat).copy(), res

        W_ref, r_ref = run("svd")
        W, r = run("interval")
        np.testing.assert_array_equal(W, W_ref)
        for f in ("collapse0", "collapse_s", "collapse_c", "static_s", "static_c"):
            np.testing.assert_array_equal(getattr(r, f), getattr(r_ref, f))
        self.assertEqual((r.cond_min_s, r.cond_min_c), (r_ref.cond_min_s, r_ref.cond_min_c))
        self.assertEqual(r_ref.guard_blocks, 2 * nmat * nvol)          # both channels, all blocks
        self.assertEqual(r_ref.guard_exact_blocks, r_ref.guard_blocks)  # the svd path decomposes all
        self.assertEqual(r.guard_blocks, 2 * nmat * nvol)
        self.assertLessEqual(r.guard_exact_blocks, r.guard_blocks)

    def test_dress_result_keeps_positional_construction(self):
        from hwave.solver.flex_bond import DressResult
        z = np.zeros((1, 1, 1, 1))
        r = DressResult(z, z, z, z, z, 1.0, 1.0)
        self.assertEqual((r.guard_violations, r.guard_exact_blocks, r.guard_blocks), (0, 0, 0))

    def test_memory_table_has_an_incremental_guard_row(self):
        from hwave.solver import flex_bond as fb
        est = fb.estimate_bond_memory(nmat=8, nvol=4, norb=2, B=3, depth=2, output_full=False,
                                      split_seed=False, n_types=1, freq_batch=None, cap_gb=200.0,
                                      mixing="anderson", device_available=8 * 1024 ** 3)
        rows = est["device_rows"]
        self.assertIn("guard", rows)
        U_b = int(est["nb"]) * 4 * (3 * 4) ** 2 * 16
        self.assertEqual(rows["guard"], 3 * U_b)


class TestErrstateHelper(unittest.TestCase):
    """``bond_channels._errstate`` (issue #197 fix round 1): numpy's warning
    suppression for the norm / product / quotient arithmetic of the interval
    guard, a no-op off numpy -- cupy has no ``errstate`` and its device
    arithmetic raises no such warnings."""

    def test_errstate_helper_is_a_no_op_off_numpy(self):
        import contextlib
        import types
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            with bc._errstate(np):
                1.0 / np.zeros(1)
        self.assertEqual(caught, [])
        obj = types.SimpleNamespace()
        self.assertIsInstance(bc._errstate(obj), contextlib.nullcontext)
