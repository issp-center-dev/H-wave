"""The provenance-split self-energy state and its mixing (spec 2026-09-06
sections 2.4 and 5.6): linear pair mixing equals the legacy total update,
the stacked Anderson mixer with REAL coefficients equals an independent
reference and applies one coefficient vector to both components, the
static block stays Hermitian, the cancellation regression, and the
three-residual pass counter."""
import unittest

import numpy as np


def _herm(rng, shape):
    a = rng.normal(size=shape) + 1j * rng.normal(size=shape)
    return 0.5 * (a + np.conj(np.swapaxes(a, -1, -2)))


def _state(rng, nmat=6, nvol=4, norb=2):
    from hwave.solver.flex_mixing import SplitState
    st = _herm(rng, (1, 1, nvol, norb, norb))
    fl = _herm(rng, (1, nmat, nvol, norb, norb))
    return SplitState(static=st, fluct=fl)


def _reference_anderson(xs, rs, mix):
    """Independent stacked-vector Anderson step with a REAL Gram system
    (spec 2.4): x_{k+1} = x + mix r - (dX + mix dR)^T gamma,
    gamma = solve(Re(dR dR^H) + lam I, Re(dR conj r))."""
    m = len(xs)
    x, r = xs[-1], rs[-1]
    if m == 1:
        return x + mix * r
    dR = np.stack([rs[i + 1] - rs[i] for i in range(m - 1)])
    dX = np.stack([xs[i + 1] - xs[i] for i in range(m - 1)])
    G = (dR.conj() @ dR.T).real
    b = (dR.conj() @ r).real
    lam = 1.0e-10 * max(float(np.trace(G)) / (m - 1), 1.0e-300)
    gamma = np.linalg.solve(G + lam * np.eye(m - 1), b)
    return x + mix * r - (dX + mix * dR).T @ gamma


class TestSplitMixing(unittest.TestCase):

    def test_linear_pair_equals_legacy_total(self):
        from hwave.solver.flex_mixing import linear_mix_pair
        rng = np.random.default_rng(0)
        s, n = _state(rng), _state(rng)
        out = linear_mix_pair(s, n, 0.3)
        legacy = (1.0 - 0.3) * s.total() + 0.3 * n.total()
        np.testing.assert_allclose(out.total(), legacy, rtol=0, atol=1e-14)
        self.assertEqual(out.static.shape, s.static.shape)

    def test_anderson_matches_reference_and_shares_coefficients(self):
        from hwave.solver.flex_mixing import StackedAndersonMixer, SplitState
        rng = np.random.default_rng(1)
        mix, depth = 0.3, 4
        mixer = StackedAndersonMixer(mix, depth)
        s = _state(rng)
        xs, rs = [], []
        for it in range(6):
            n = SplitState(static=s.static + 0.1 * _herm(rng, s.static.shape),
                           fluct=s.fluct + 0.1 * _herm(rng, s.fluct.shape))
            x = mixer.stack(s); r = mixer.stack(n) - x
            xs.append(x); rs.append(r)
            if len(xs) > depth:
                xs.pop(0); rs.pop(0)
            ref = _reference_anderson(xs, rs, mix)
            out = mixer.step(s, n)
            np.testing.assert_allclose(mixer.stack(out), ref, rtol=0, atol=1e-12)
            # one real coefficient vector for both components: the sum of the
            # component outputs equals the same affine combination of the summed histories
            gamma = mixer.last_gamma
            if gamma is not None:
                nb = s.static.size * s.fluct.shape[1]      # the static block is stacked BROADCAST
                tot_x = [xx[:nb].reshape(s.fluct.shape) + xx[nb:].reshape(s.fluct.shape) for xx in xs]
                tot_r = [rr[:nb].reshape(s.fluct.shape) + rr[nb:].reshape(s.fluct.shape) for rr in rs]
                dR = np.stack([tot_r[i + 1] - tot_r[i] for i in range(len(tot_r) - 1)])
                dX = np.stack([tot_x[i + 1] - tot_x[i] for i in range(len(tot_x) - 1)])
                exp_total = tot_x[-1] + mix * tot_r[-1] - np.tensordot(gamma, dX + mix * dR, axes=(0, 0))
                np.testing.assert_allclose(out.total(), exp_total, rtol=0, atol=1e-12)
            self.assertTrue(np.allclose(out.static, np.conj(np.swapaxes(out.static, -1, -2)), atol=1e-13))
            s = out

    def test_cancellation_is_not_convergence(self):
        from hwave.solver.flex_mixing import SplitState, residuals
        rng = np.random.default_rng(2)
        s = _state(rng)
        X = _herm(rng, s.static.shape)
        n = SplitState(static=s.static + X, fluct=s.fluct - np.repeat(X, s.fluct.shape[1], axis=1))
        G = rng.normal(size=s.fluct.shape) + 0j
        res_sigma, res_g, res_comp = residuals(s, n, G, G)
        self.assertLess(res_sigma, 1e-14)
        self.assertLess(res_g, 1e-300 + 1e-14)
        self.assertGreater(res_comp, 1e-2)

    def test_pass_counter_three_consecutive(self):
        from hwave.solver.flex_mixing import PassCounter
        pc = PassCounter(eps=1e-6, needed=3)
        self.assertEqual(pc.update(1e-7, 1e-7, 1e-7), (False, 1))
        self.assertEqual(pc.update(1e-7, 1e-3, 1e-7), (False, 0))     # any residual above eps resets
        for expect in ((False, 1), (False, 2), (True, 3)):
            self.assertEqual(pc.update(1e-8, 1e-8, 1e-8), expect)


if __name__ == "__main__":
    unittest.main()
