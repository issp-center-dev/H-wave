"""Numerical-equivalence pin for the _green per-k Green-function build rewrite.

The UHFk solver builds the per-block k-space Green function as

    G_ab(k) = sum_l conj(V_kal) * dist_kl * V_kbl

which was implemented with np.einsum('kal,kl,kbl->kab', conj(v), dist, v) and is
rewritten to an explicit batched matmul for performance. This test pins that the
matmul form is numerically equivalent (index order preserved) to the einsum.
"""

import numpy as np
import unittest


def green_einsum(v, dist):
    return np.einsum('kal, kl, kbl -> kab', np.conjugate(v), dist, v)


def green_matmul(v, dist):
    return np.matmul(np.conjugate(v) * dist[:, None, :], v.transpose(0, 2, 1))


class TestGreenMatmulEquivalence(unittest.TestCase):
    def _check(self, nvol, nd, seed):
        rng = np.random.default_rng(seed)
        v = (rng.standard_normal((nvol, nd, nd))
             + 1j * rng.standard_normal((nvol, nd, nd)))
        dist = rng.random((nvol, nd))

        ref = green_einsum(v, dist)
        got = green_matmul(v, dist)

        maxdiff = np.max(np.abs(ref - got))
        self.assertTrue(
            np.allclose(ref, got, rtol=1e-11, atol=1e-11),
            msg="maxdiff={:.3e} for nvol={}, nd={}".format(maxdiff, nvol, nd),
        )
        # tighter pin on the reduction-order difference
        self.assertLess(maxdiff, 1e-11)

    def test_small(self):
        self._check(nvol=4, nd=3, seed=0)

    def test_medium(self):
        self._check(nvol=27, nd=8, seed=1)

    def test_larger(self):
        self._check(nvol=64, nd=12, seed=2)

    def test_single_orbital(self):
        self._check(nvol=8, nd=1, seed=3)


if __name__ == "__main__":
    unittest.main()
