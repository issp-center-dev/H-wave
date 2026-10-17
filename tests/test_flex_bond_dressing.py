"""Frequency-batched bond dressing, the effective interaction and the
collapses (spec 2026-09-06 sections 3.2-3.3)."""
import types
import unittest

import numpy as np


def _problem(nmat=6, nvol=4, nd=4, B=3, seed=0):
    rng = np.random.default_rng(seed)
    ND = B * nd
    chi_bar = 0.05 * (rng.normal(size=(nmat, nvol, ND, ND)) + 1j * rng.normal(size=(nmat, nvol, ND, ND)))
    S = 0.4 * rng.normal(size=(nvol, ND, ND)); S = 0.5 * (S + np.swapaxes(S, 1, 2))
    C = 0.3 * rng.normal(size=(nvol, ND, ND)); C = 0.5 * (C + np.swapaxes(C, 1, 2))
    return chi_bar, S + 0j, C + 0j


class TestDressBatch(unittest.TestCase):

    def test_equals_per_frequency_dress_channel(self):
        from hwave.solver import bond_channels as bc
        chi_bar, S, C = _problem()
        nmat, nvol = chi_bar.shape[:2]
        for channel, W in (("spin", S), ("charge", C)):
            chi_b, cond = bc.dress_batch(chi_bar[1:4], W, channel, l0=1, nmat=nmat, spatial_shape=(4, 1, 1))
            for i, l in enumerate(range(1, 4)):
                ref, c_ref = bc.dress_channel(chi_bar[l], W, channel, spatial_shape=(4, 1, 1))
                np.testing.assert_array_equal(chi_b[i], ref)
            self.assertIsInstance(cond, float)

    def test_accepts_an_array_like_input(self):
        """The batch argument is documented as an array, not specifically a
        numpy one: an array-like (nested list) must be coerced, not crash on
        a missing ``ndim``."""
        from hwave.solver import bond_channels as bc
        chi_bar, S, _C = _problem(nmat=2, nvol=2, nd=1, B=1, seed=3)
        ref, _ = bc.dress_batch(chi_bar, S, "spin", l0=0, nmat=2, spatial_shape=(2, 1, 1))
        out, _ = bc.dress_batch(chi_bar.tolist(), S, "spin", l0=0, nmat=2,
                                spatial_shape=(2, 1, 1))
        np.testing.assert_array_equal(out, ref)

    def test_refusal_names_frequency_and_q(self):
        from hwave.solver import bond_channels as bc
        chi_bar, S, C = _problem(nmat=4, nvol=2, nd=1, B=1)
        chi_bar[:] = 0.0
        chi_bar[2, 1, 0, 0] = 1.0                            # l=2, q index 1: 1 - 1*1 = 0
        W = np.ones((2, 1, 1), complex)
        with self.assertRaises(ValueError) as cm:
            bc.dress_batch(chi_bar[2:4], W, "spin", l0=2, nmat=4, spatial_shape=(2, 1, 1))
        msg = str(cm.exception)
        self.assertIn("bosonic", msg)
        self.assertIn(str(2 * 2 - 4), msg)                     # 2l - nmat = 0
        self.assertIn("(1, 0, 0)", msg)
        self.assertIn("spin", msg)


def _onsite_factors():
    """A small on-site-only factor pack (vpair = None, so it carries no
    nvol-dependent array and fits the synthetic nd = 4 problem above)."""
    from hwave.solver.second_order import build_factors
    from tests.test_second_order_factors import _split_for
    s, split = _split_for({"CoulombIntra": [(0, 0, 0, 1, 1, 1.3, 0.0),
                                            (0, 0, 0, 2, 2, 0.9, 0.0)],
                           "Exchange": [(0, 0, 0, 1, 2, 0.4, 0.0),
                                        (0, 0, 0, 2, 1, 0.4, 0.0)]})
    f = build_factors(split, s.lattice, s.norb)
    assert f.vpair is None and f.nd == 4
    return f


class TestDressAndBuildW(unittest.TestCase):

    def test_w_and_collapses_and_static(self):
        """The rev-19 channel-0 kernel ("takimoto", explicit) and the exact
        local one ("local"); the mixed and bond-bond blocks are the same
        under both.

        The mixed (channel-0 x bond) blocks are read at the PERMUTED pair
        index on both axes -- the orbital-pair transpose inside every bond
        block, spec 2026-09-16 R3 (issue #192) -- and the permutation is
        built here independently of the production helper, from the block
        layout alone. The synthetic problem has ``nd = 4`` (norb = 2) and a
        random, orbital-asymmetric ``chi_bar``/``S``/``C``, so it sees the
        permutation."""
        from hwave.solver.flex_bond import BondBlockStore, BondDeviceContext, dress_and_build_w
        from hwave.solver.second_order import dense_w2
        chi_bar, S, C = _problem(nmat=6, nvol=4, nd=4, B=3)
        nmat, nvol, ND = chi_bar.shape[:3]
        nd = 4
        factors = _onsite_factors()
        for output_full, second_order in ((False, "takimoto"), (True, "takimoto"),
                                          (False, "local"), (True, "local")):
            with self.subTest(output_full=output_full, second_order=second_order):
                names = ("chibar", "W") + (("chi_s_w", "chi_c_w") if output_full else ())
                with BondBlockStore(nmat, nvol, ND, nd, names) as store:
                    store.put_freq_batch("chibar", 0, nmat, chi_bar)
                    S_on = np.ascontiguousarray(S[:1, :nd, :nd]).repeat(nvol, axis=0) * 0.7
                    C_on = np.ascontiguousarray(C[:1, :nd, :nd]).repeat(nvol, axis=0) * 0.3
                    view = types.SimpleNamespace(n_channels=ND // nd)
                    norb = int(round(nd ** 0.5))
                    with BondDeviceContext.for_view(np, S, C, S_on, C_on, view, norb) as dev:
                        res = dress_and_build_w(store, dev, nb=4,
                                                output_full=output_full, nmat=nmat,
                                                nvol=nvol, nd=nd, spatial_shape=(4, 1, 1),
                                                factors=factors, second_order=second_order)
                    I = np.eye(ND)
                    chi_s = np.linalg.solve(I - chi_bar @ S, chi_bar)
                    chi_c = np.linalg.solve(I + chi_bar @ C, chi_bar)
                    # spec 3.3 rev 19: ring beyond second order + the exact second order
                    W_ref = 1.5 * S @ (chi_s - chi_bar) @ S + 0.5 * C @ (chi_c - chi_bar) @ C
                    A = S @ chi_bar @ S
                    Bc = C @ chi_bar @ C
                    W2 = np.zeros_like(W_ref)
                    if second_order == "takimoto":
                        W2[:, :, :nd, :nd] = (1.5 * A + 0.5 * Bc)[:, :, :nd, :nd] \
                            - 0.25 * (S_on + C_on) @ chi_bar[:, :, :nd, :nd] @ (S_on + C_on)
                    else:
                        # spec 2026-09-08 D5: the exact local kernel on channel 0
                        W2[:, :, :nd, :nd] = dense_w2(chi_bar[:, :, :nd, :nd], factors)
                    # the pair permutation of the mixed blocks: identity on
                    # channel 0 (never reached below, which only reads the
                    # bond half), the orbital-pair transpose (l1, l2) ->
                    # (l2, l1) inside every bond block
                    norb = int(round(nd ** 0.5))
                    pb = np.arange(nd, ND)
                    for m in range(1, ND // nd):
                        for l1 in range(norb):
                            for l2 in range(norb):
                                pb[m * nd + l1 * norb + l2 - nd] = m * nd + l2 * norb + l1
                    AB = 0.25 * (A + Bc)
                    W2[:, :, :nd, nd:] = AB[:, :, :nd, pb]
                    W2[:, :, nd:, :nd] = AB[:, :, pb, :nd]
                    W_ref = W_ref + W2
                    np.testing.assert_allclose(store.get_freq_batch("W", 0, nmat), W_ref, rtol=1e-12, atol=1e-13)
                    np.testing.assert_allclose(res.collapse0, chi_bar[:, :, :nd, :nd], rtol=0, atol=1e-14)
                    np.testing.assert_allclose(res.collapse_s, chi_s[:, :, :nd, :nd], rtol=0, atol=1e-12)
                    np.testing.assert_allclose(res.collapse_c, chi_c[:, :, :nd, :nd], rtol=0, atol=1e-12)
                    np.testing.assert_allclose(res.static_s, chi_s[nmat // 2], rtol=0, atol=1e-12)
                    np.testing.assert_allclose(res.static_c, chi_c[nmat // 2], rtol=0, atol=1e-12)
                    self.assertTrue(res.static_s.flags.owndata)
                    self.assertGreater(res.cond_min_s, 0.0)
                    if output_full:
                        np.testing.assert_allclose(store.get_freq_batch("chi_s_w", 0, nmat), chi_s, rtol=0, atol=1e-12)
                        np.testing.assert_allclose(store.get_freq_batch("chi_s_w", nmat // 2, nmat // 2 + 1)[0],
                                                   res.static_s, atol=0)
                    else:
                        with self.assertRaises(KeyError):
                            store.get_freq_batch("chi_s_w", 0, 1)


class TestDressAndBuildWDevice(unittest.TestCase):

    def test_numpy_context_matches_reference_batch_sizes(self):
        """The refactored loop with the numpy context equals itself for every batch
        size (1, non-divisor, boundary containing l_static, full) to round-off, and
        the store's W is written once per batch (spy)."""
        from hwave.solver import flex_bond as fb
        chi_bar, S, C = _problem(nmat=8, nvol=4, nd=4, B=3, seed=5)
        nmat, nvol, ND = chi_bar.shape[0], chi_bar.shape[1], S.shape[-1]
        nd = 4
        S_on = S[:, :nd, :nd].copy(); C_on = C[:, :nd, :nd].copy()
        view = types.SimpleNamespace(n_channels=ND // nd)
        outs = {}
        for nb in (1, 3, nmat // 2, nmat):
            with fb.BondBlockStore(nmat, nvol, ND, nd, ("chibar", "W")) as store, \
                    fb.BondDeviceContext.for_view(np, S, C, S_on, C_on, view, 2) as dev:
                store.put_freq_batch("chibar", 0, nmat, chi_bar)
                puts = []
                orig = store.put_freq_batch
                def _spy(name, l0, l1, batch, _o=orig):
                    puts.append((name, l0, l1)); return _o(name, l0, l1, batch)
                store.put_freq_batch = _spy
                res = fb.dress_and_build_w(store, dev, nb=nb, output_full=False, nmat=nmat,
                                           nvol=nvol, nd=nd, spatial_shape=(4, 1, 1),
                                           second_order="takimoto")
                outs[nb] = (store._arrays["W"].copy(), res.collapse0.copy(), res.static_s.copy(),
                            res.cond_min_s)
                self.assertEqual([p for p in puts if p[0] == "W"],
                                 [("W", l0, min(nmat, l0 + nb)) for l0 in range(0, nmat, nb)])
        for nb in (1, 3, nmat // 2):
            for a, b in zip(outs[nb], outs[nmat]):
                if isinstance(a, float):
                    self.assertAlmostEqual(a, b, places=12)
                else:
                    np.testing.assert_allclose(a, b, rtol=1e-12, atol=1e-14)


class TestDressGuardPolicyPassThrough(unittest.TestCase):
    """``flex_bond._dress`` owns the SCF iteration, so both the refusals and
    the warnings of a tolerated violation must name it (GitHub issue #199)."""

    @staticmethod
    def _near_singular(nb=1, nvol=2):
        """(chi_bar_b, W) whose spin denominator is nearly singular at q = 1
        in every frequency of the batch."""
        cb = np.zeros((nb, nvol, 2, 2), complex)
        W = np.zeros((nvol, 2, 2), complex)
        cb[:, :] = np.eye(2)
        W[1] = np.eye(2) - np.array([[1.0, 1.0], [1.0, 1.0 + 1.0e-13]], complex)
        return cb, W

    def test_the_conditioning_warning_names_the_iteration(self):
        from hwave.solver import flex_bond
        cb, W = self._near_singular()
        violations = []
        with self.assertLogs("qlms.solver.bond_channels", level="WARNING") as cm:
            chi, cond = flex_bond._dress(cb, W, "spin", 0, 2, (2, 1, 1), 1.0e-3, 3,
                                         "all", "warn", violations)
        self.assertIn("(SCF iteration 3)", "\n".join(cm.output))
        self.assertTrue(np.all(np.isfinite(chi)))
        self.assertEqual([v["iteration"] for v in violations], [3])

    def test_the_static_mode_warning_names_the_iteration(self):
        """``guard_freqs = "static"`` still SVD-checks the zero-frequency
        slice, and that warning carries the iteration too. (The residual
        guard of the unchecked slices is reached only through an explicit
        ``residual_tol``, which ``_dress`` does not expose; it is covered at
        the ``dress_batch`` level.)"""
        from hwave.solver import flex_bond
        cb, W = self._near_singular(nb=2)            # the batch holds l = 0 and l = 1
        violations = []
        with self.assertLogs("qlms.solver.bond_channels", level="WARNING") as cm:
            chi, _ = flex_bond._dress(cb, W, "spin", 0, 2, (2, 1, 1), 1.0e-3, 3,
                                      "static", "warn", violations)
        self.assertIn("(SCF iteration 3)", "\n".join(cm.output))
        self.assertTrue(np.all(np.isfinite(chi)))
        self.assertEqual([v["iteration"] for v in violations], [3])

    def test_a_refusal_still_names_the_iteration_once(self):
        from hwave.solver import flex_bond
        cb, W = self._near_singular()
        with self.assertRaises(ValueError) as cm:
            flex_bond._dress(cb, W, "spin", 0, 2, (2, 1, 1), 1.0e-3, 3, "all", "refuse", None)
        msg = str(cm.exception)
        self.assertEqual(msg.count("SCF iteration 3"), 1)


class TestDressBatchGuardModes(unittest.TestCase):

    def test_all_is_the_default_and_unchanged(self):
        from hwave.solver import bond_channels as bc
        chi_bar, S, C = _problem()
        nmat = chi_bar.shape[0]
        ref, c_ref = bc.dress_batch(chi_bar, S, "spin", l0=0, nmat=nmat, spatial_shape=(4, 1, 1))
        out, c_out = bc.dress_batch(chi_bar, S, "spin", l0=0, nmat=nmat, spatial_shape=(4, 1, 1),
                                    guard_freqs="all")
        np.testing.assert_array_equal(out, ref)
        self.assertEqual(c_out, c_ref)

    def test_static_checks_only_the_zero_frequency_slice(self):
        """A batch singular ONLY at a nonzero bosonic frequency: refused by "all",
        accepted by "static" when the solve is accurate (here the block is exactly
        singular, so the solve raises / the residual check refuses -- both name the
        frequency); a batch singular ONLY at l = nmat//2 is refused by both."""
        from hwave.solver import bond_channels as bc
        nmat, nvol = 4, 2
        W = np.ones((nvol, 1, 1), complex)
        # near-singular (not exactly) at l=3 (bosonic index 2), q=1: 1 - 0.999999 = 1e-6
        chi_bar = np.zeros((nmat, nvol, 1, 1), complex)
        chi_bar[3, 1, 0, 0] = 1.0 - 1e-6
        with self.assertRaises(ValueError) as cm:
            bc.dress_batch(chi_bar, W, "spin", l0=0, nmat=nmat, spatial_shape=(2, 1, 1),
                           guard_freqs="all")
        self.assertIn("bosonic Matsubara index 2", str(cm.exception))
        chi, cond = bc.dress_batch(chi_bar, W, "spin", l0=0, nmat=nmat, spatial_shape=(2, 1, 1),
                                   guard_freqs="static")
        self.assertTrue(np.all(np.isfinite(chi)))
        self.assertAlmostEqual(cond, 1.0)                    # the static slice is the identity
        # singular at the static slice: both modes refuse
        chi_bar2 = np.zeros((nmat, nvol, 1, 1), complex)
        chi_bar2[nmat // 2, 0, 0, 0] = 1.0
        for mode in ("all", "static"):
            with self.assertRaises(ValueError):
                bc.dress_batch(chi_bar2, W, "spin", l0=0, nmat=nmat, spatial_shape=(2, 1, 1),
                               guard_freqs=mode)

    def test_static_residual_check_refuses_an_inaccurate_unchecked_solve(self):
        from hwave.solver import bond_channels as bc
        nmat, nvol = 4, 1
        W = np.ones((nvol, 1, 1), complex)
        chi_bar = np.zeros((nmat, nvol, 1, 1), complex)
        chi_bar[3, 0, 0, 0] = 1.0                            # exactly singular at l=3 (unchecked)
        with self.assertRaises(ValueError) as cm:
            bc.dress_batch(chi_bar, W, "spin", l0=0, nmat=nmat, spatial_shape=(1, 1, 1),
                           guard_freqs="static")
        msg = str(cm.exception)
        self.assertIn("static", msg)
        self.assertIn("spin", msg)

    def test_invalid_guard_mode_refused(self):
        from hwave.solver import bond_channels as bc
        chi_bar, S, C = _problem()
        with self.assertRaises(ValueError):
            bc.dress_batch(chi_bar, S, "spin", l0=0, nmat=6, spatial_shape=(4, 1, 1),
                           guard_freqs="none")

    def test_solve_residual_is_zero_safe(self):
        from hwave.solver import bond_channels as bc
        mat = np.eye(3, dtype=complex)[None, None]
        cb = np.zeros((1, 1, 3, 3), complex)
        chi = np.zeros_like(cb)
        r = bc.solve_residual(mat, chi, cb)
        self.assertEqual(r.shape, (1, 1))
        self.assertEqual(float(r[0, 0]), 0.0)

    def test_static_mode_skips_the_guard_outside_the_static_slice(self):
        """A sub-batch that does not contain l = nmat // 2 gets no SVD guard
        under "static" (cond_min is None); the sub-batch that does contain it
        gets a float. The concatenated chi still matches the "all" result."""
        from hwave.solver import bond_channels as bc
        chi_bar, S, C = _problem()
        nmat = chi_bar.shape[0]
        l_static = nmat // 2
        full, _ = bc.dress_batch(chi_bar, S, "spin", l0=0, nmat=nmat, spatial_shape=(4, 1, 1))
        nb = 2
        parts, conds = [], []
        for l0 in range(0, nmat, nb):
            l1 = min(nmat, l0 + nb)
            p, c = bc.dress_batch(chi_bar[l0:l1], S, "spin", l0=l0, nmat=nmat,
                                  spatial_shape=(4, 1, 1), guard_freqs="static")
            parts.append(p); conds.append((l0, l1, c))
        saw_none = False
        for l0, l1, c in conds:
            if l0 <= l_static < l1:
                self.assertIsInstance(c, float)
            else:
                self.assertIsNone(c)
                saw_none = True
        self.assertTrue(saw_none)                             # at least one excluded sub-batch
        np.testing.assert_allclose(np.concatenate(parts), full, rtol=1e-12, atol=1e-14)

    def test_batch_size_invariance_of_chi_and_cond(self):
        from hwave.solver import bond_channels as bc
        chi_bar, S, C = _problem(nmat=8, nvol=4, nd=2, B=2, seed=3)
        nmat = chi_bar.shape[0]
        full, c_full = bc.dress_batch(chi_bar, S, "spin", l0=0, nmat=nmat, spatial_shape=(4, 1, 1))
        for nb in (1, 3, nmat // 2):
            parts, conds = [], []
            for l0 in range(0, nmat, nb):
                l1 = min(nmat, l0 + nb)
                p, c = bc.dress_batch(chi_bar[l0:l1], S, "spin", l0=l0, nmat=nmat,
                                      spatial_shape=(4, 1, 1))
                parts.append(p); conds.append(c)
            np.testing.assert_allclose(np.concatenate(parts), full, rtol=1e-12, atol=1e-14)
            self.assertAlmostEqual(min(conds), c_full, places=12)


if __name__ == "__main__":
    unittest.main()
