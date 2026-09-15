"""The batched local kernel (spec 2.3-2.5): the Hubbard identity
W2 = U^2 chibar, the V-V and U-V checks, batch independence, Hermiticity
W2(q, i nu)^dagger = W2(q, -i nu), in-place accumulation, and the non-finite
checkpoint."""
import unittest

import numpy as np

from tests.test_second_order_factors import _split_for


def _factors(rows_by_type):
    from hwave.solver.second_order import build_factors
    s, split = _split_for(rows_by_type)
    return s, build_factors(split, s.lattice, 2)


def _random_chibar(nmat, nvol, nd, seed, hermitian_pair=True):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(nmat, nvol, nd, nd)) + 1j * rng.normal(size=(nmat, nvol, nd, nd))
    if hermitian_pair:      # chibar(q, l)^dagger = chibar(q, nmat - l)
        Y = np.empty_like(X)
        for l in range(nmat):
            Y[l] = 0.5 * (X[l] + X[(nmat - l) % nmat].conj().swapaxes(-1, -2))
        return Y
    return X


class TestKernel(unittest.TestCase):

    def test_hubbard_identity(self):
        from hwave.solver.second_order import dense_w2
        s, f = _factors({"CoulombIntra": [(0, 0, 0, 1, 1, 0.7, 0.0), (0, 0, 0, 2, 2, 0.7, 0.0)]})
        cb = _random_chibar(8, 16, 4, 0)
        W2 = dense_w2(cb, f)
        U2 = np.zeros((4, 4)); U2[0, 0] = U2[3, 3] = 0.49       # density slots (aa),(aa)
        # explicit: W2[(aa),(bb)] = U^2 chibar[(aa),(bb)] for a == b, zero elsewhere
        for i in range(4):
            for j in range(4):
                if i in (0, 3) and j in (0, 3) and i == j:
                    np.testing.assert_allclose(W2[:, :, i, j], 0.49 * cb[:, :, i, j], atol=1e-13)
                elif i in (0, 3) and j in (0, 3):
                    np.testing.assert_allclose(W2[:, :, i, j], 0.49 * cb[:, :, i, j], atol=1e-13)
                else:
                    self.assertLess(np.abs(W2[:, :, i, j]).max(), 1e-13)

    def test_offsite_v_and_uv_checks(self):
        """One orbital active (orbital 2 undeclared): W2_00 = 2 V(q)^2 chibar_00 and
        the U-V cross term 2 U V(q) chibar_00 (spec 2.4 checks)."""
        from hwave.solver.second_order import dense_w2
        rows_v = [(1, 0, 0, 1, 1, 0.3, 0.0), (-1, 0, 0, 1, 1, 0.3, 0.0)]
        s, fv = _factors({"CoulombInter": rows_v})
        s, fuv = _factors({"CoulombInter": rows_v, "CoulombIntra": [(0, 0, 0, 1, 1, 0.5, 0.0)]})
        s, fu = _factors({"CoulombIntra": [(0, 0, 0, 1, 1, 0.5, 0.0)]})
        cb = _random_chibar(8, 16, 4, 1)
        vq = fv.vpair[0, 0, :, 0, 0]                                  # V(q) on (00),(00)
        Wv = dense_w2(cb, fv)
        np.testing.assert_allclose(Wv[:, :, 0, 0], 2.0 * vq[None] ** 2 * cb[:, :, 0, 0], atol=1e-13)
        Wuv = dense_w2(cb, fuv)
        Wu = dense_w2(cb, fu)
        cross = Wuv - Wu - Wv
        np.testing.assert_allclose(cross[:, :, 0, 0], 2.0 * 0.5 * vq[None] * cb[:, :, 0, 0], atol=1e-13)

    def test_batch_independence_and_in_place(self):
        from hwave.solver.second_order import accumulate_batch, dense_w2
        s, f = _factors({"CoulombIntra": [(0, 0, 0, 1, 1, 0.7, 0.0), (0, 0, 0, 2, 2, 0.4, 0.0)],
                         "Hund": [(0, 0, 0, 1, 2, 0.2, 0.0), (0, 0, 0, 2, 1, 0.2, 0.0)],
                         "CoulombInter": [(1, 0, 0, 1, 2, 0.3, 0.0), (-1, 0, 0, 2, 1, 0.3, 0.0)]})
        cb = _random_chibar(8, 16, 4, 2)
        ref = dense_w2(cb, f)
        for nb in (1, 3, 8):
            out = np.ones((8, 16, 4, 4), complex)          # pre-filled: accumulation adds
            for l0 in range(0, 8, nb):
                accumulate_batch(out[l0:l0 + nb], cb[l0:l0 + nb], l0, f)
            np.testing.assert_allclose(out - 1.0, ref, rtol=0, atol=1e-12)

    def test_hermiticity_without_inversion_symmetry(self):
        from hwave.solver.second_order import dense_w2
        s, f = _factors({"CoulombIntra": [(0, 0, 0, 1, 1, 0.7, 0.0)],
                         "Exchange": [(0, 0, 0, 1, 2, 0.2, 0.1), (0, 0, 0, 2, 1, 0.2, -0.1)],
                         "PairHop": [(0, 0, 0, 1, 2, 0.15, 0.05), (0, 0, 0, 2, 1, 0.15, -0.05)],
                         "CoulombInter": [(1, 0, 0, 1, 2, 0.3, 0.0), (-1, 0, 0, 2, 1, 0.3, 0.0)]})
        cb = _random_chibar(8, 16, 4, 3)
        W2 = dense_w2(cb, f)
        for l in range(8):
            np.testing.assert_allclose(W2[l].conj().swapaxes(-1, -2), W2[(8 - l) % 8], atol=1e-12)

    def test_non_finite_checkpoint(self):
        from hwave.solver.hartree_fock import NonFiniteError
        from hwave.solver.second_order import accumulate_batch
        s, f = _factors({"CoulombIntra": [(0, 0, 0, 1, 1, 0.7, 0.0)]})
        cb = _random_chibar(4, 16, 4, 4)
        cb[1, 2, 0, 0] = np.nan
        out = np.zeros_like(cb)
        with self.assertRaises(NonFiniteError) as cm:
            accumulate_batch(out, cb, 0, f)
        self.assertIn("second-order kernel W2", str(cm.exception))
        self.assertIn("[0, 4)", str(cm.exception))


if __name__ == "__main__":
    unittest.main()
