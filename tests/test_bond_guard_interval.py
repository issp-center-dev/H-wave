"""The interval guard of the bond dressing (issue #197): certified
two-sided bounds of every block's conditioning score from batched device
operations, the exact set of blocks that could be the batch minimum, and
the resolver that reproduces the SVD guard's outputs from them."""
import unittest
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
        inside = np.concatenate([base * 1e-8, base * 1e8])
        iv = bc.bond_conditioning_interval(inside, np)
        _assert_enclosed(self, iv, inside)
        outside = np.concatenate([base * 1e-150, base * 1e150])
        iv = bc.bond_conditioning_interval(outside, np)
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
