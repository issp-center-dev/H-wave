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


def _factors_norb1(rows_by_type, shape=(4, 4, 1)):
    """Factors of a ONE-orbital fixture (nd = 1: the density slot is the whole
    pair axis, the worst case for the kernel's slot handling) on a lattice of
    the given shape."""
    import os
    import shutil
    import tempfile

    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as flex_mod
    from hwave.solver.offsite import split_locality
    from hwave.solver.second_order import build_factors
    from tests.test_second_order_factors import _write_wan

    d = tempfile.mkdtemp()
    for f in ("geom.dat", "transfer.dat"):
        shutil.copy(os.path.join("tests/rpa/input", f), d)
    idict = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
    for t, rows in rows_by_type.items():
        _write_wan(os.path.join(d, t.lower() + ".dat"), t, 1, rows)
        idict[t] = t.lower() + ".dat"
    r = read_input_k.QLMSkInput({"path_to_input": d, "interaction": idict})
    par = {"T": 2.0, "filling": 0.5, "CellShape": list(shape), "SubShape": [1, 1, 1],
           "Nmat": 8, "IterationMax": 1, "Mix": 1.0, "EPS": 1, "flex_second_order": "takimoto"}
    s = flex_mod.FLEX(r.get_param("ham"), {}, {"mode": "FLEX", "param": par,
                                               "enable_spin_orbital": False,
                                               "calc_scheme": "general"})
    return build_factors(split_locality(s.ham_info, s.lattice), s.lattice, 1)


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
        # the (aa),(bb) density block is scaled by U^2, zero elsewhere
        for i in range(4):
            for j in range(4):
                if i in (0, 3) and j in (0, 3):
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

    def test_offsite_spin_selection_pinned_by_ising(self):
        """Ising's opposite-spin density weight (-1, vs CoulombInter's +1) pins
        the mixed term's spin selection sig = sq (spec 2.3-2.4): for Hubbard U
        the only surviving on-site triple under the accumulate_batch gate is
        (UP, DN, DN), so the mixed term uses vpair[UP, DN] = vpair[DN, UP],
        which for Ising is -I(q) (vpair[UP, UP] is +I(q)); a mutation
        sig = UP would instead read vpair[UP, UP] = +I(q) and flip the sign
        below. Hund has no opposite-spin density coupling, so its U-Hund
        cross term is zero -- no sign pressure from that type."""
        from hwave.solver.second_order import dense_w2
        rows_i = [(1, 0, 0, 1, 1, 0.3, 0.0), (-1, 0, 0, 1, 1, 0.3, 0.0)]
        s, fi = _factors({"Ising": rows_i})
        s, fui = _factors({"Ising": rows_i, "CoulombIntra": [(0, 0, 0, 1, 1, 0.5, 0.0)]})
        s, fu = _factors({"CoulombIntra": [(0, 0, 0, 1, 1, 0.5, 0.0)]})
        cb = _random_chibar(8, 16, 4, 5)
        iq = fi.vpair[0, 0, :, 0, 0]                                  # +I(q) on (00),(00)
        np.testing.assert_allclose(fi.vpair[0, 1, :, 0, 0], -iq, atol=1e-13)
        Wui = dense_w2(cb, fui)
        Wu = dense_w2(cb, fu)
        Wi = dense_w2(cb, fi)
        cross = Wui - Wu - Wi
        # note the NEGATIVE of the CoulombInter cross term of test_offsite_v_and_uv_checks
        np.testing.assert_allclose(cross[:, :, 0, 0], -2.0 * 0.5 * iq[None] * cb[:, :, 0, 0], atol=1e-13)

        s, fh = _factors({"Hund": rows_i})
        s, fuh = _factors({"Hund": rows_i, "CoulombIntra": [(0, 0, 0, 1, 1, 0.5, 0.0)]})
        Wuh = dense_w2(cb, fuh)
        Wh = dense_w2(cb, fh)
        cross_h = Wuh - Wu - Wh
        self.assertLess(np.abs(cross_h[:, :, 0, 0]).max(), 1e-13)

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

    def test_caller_supplied_work_buffers(self):
        """``work`` lends the kernel its two temporaries (spec 2.5: the whole
        general-path assembly peaks at two, not four). The result must be
        bit-identical to the self-allocating call, the buffers' contents on
        entry must not matter, and a mis-shaped buffer must be refused."""
        from hwave.solver.second_order import accumulate_batch
        s, f = _factors({"CoulombIntra": [(0, 0, 0, 1, 1, 0.7, 0.0), (0, 0, 0, 2, 2, 0.4, 0.0)],
                         "Hund": [(0, 0, 0, 1, 2, 0.2, 0.0), (0, 0, 0, 2, 1, 0.2, 0.0)],
                         "CoulombInter": [(1, 0, 0, 1, 2, 0.3, 0.0), (-1, 0, 0, 2, 1, 0.3, 0.0)]})
        cb = _random_chibar(8, 16, 4, 5)
        ref = np.zeros_like(cb)
        accumulate_batch(ref, cb, 0, f)
        self.assertGreater(np.abs(ref).max(), 1e-6)                    # anti-vacuity
        # dirty buffers on entry: the kernel must overwrite, never read them
        T1 = np.full_like(cb, 7.0 - 3.0j)
        T2 = np.full_like(cb, -11.0 + 2.0j)
        got = np.zeros_like(cb)
        accumulate_batch(got, cb, 0, f, work=(T1, T2))
        np.testing.assert_array_equal(got, ref)
        # a sliced view of a longer buffer is what the general path lends
        big1, big2 = np.empty((12, 16, 4, 4), complex), np.empty((12, 16, 4, 4), complex)
        got2 = np.zeros_like(cb)
        accumulate_batch(got2, cb, 0, f, work=(big1[:8], big2[:8]))
        np.testing.assert_array_equal(got2, ref)
        # wrong shape is refused, naming which buffer and both shapes
        for bad in ((np.empty((7, 16, 4, 4), complex), np.empty_like(cb)),
                    (np.empty_like(cb), np.empty((8, 16, 4, 2), complex))):
            with self.assertRaises(ValueError) as cm:
                accumulate_batch(np.zeros_like(cb), cb, 0, f, work=bad)
            self.assertIn("accumulate_batch", str(cm.exception))
            self.assertIn("expected chibar_b's (8, 16, 4, 4)", str(cm.exception))

    def test_allocation_peak_with_lent_buffers(self):
        """Spec 2.5's budget, measured: with ``work`` supplied the kernel
        allocates NOTHING that scales with the batch. Every product goes
        through ``matmul(..., out=)`` into the lent buffers and every
        density-slot gather and scatter is a strided VIEW (the (a, a) pair
        slots are the arithmetic sequence a * (norb + 1), so basic slicing
        selects them), so what is left is the boolean mask of the finiteness
        checkpoint -- one byte per complex element, a sixteenth of the batch
        -- plus numpy's own buffering for the strided in-place updates, which
        is capped by its buffer size (8192 elements per operand) and does NOT
        grow with the batch.

        Both facts are asserted: an absolute bound at one size, and -- the
        assertion that actually catches a batch-shaped allocation -- that
        quadrupling the batch grows the peak by no more than the mask does (a
        batch-shaped array would add sixteen times that).

        Measured around the ``accumulate_batch`` call ALONE: the output, the
        chibar and both work buffers are built before the window, and a
        warm-up call precedes it. The one-orbital cases are the worst case
        (nd = 1: a density slot IS the whole pair axis, so a slot-shaped copy
        would be full-size)."""
        import tracemalloc
        from hwave.solver.second_order import accumulate_batch, dense_w2
        _ALLOW = 1 << 20            # numpy's strided-operand buffers, constant
        rows_v = [(1, 0, 0, 1, 1, 0.3, 0.0), (-1, 0, 0, 1, 1, 0.3, 0.0)]
        intra1 = [(0, 0, 0, 1, 1, 0.5, 0.0)]
        _, f2 = _factors({"CoulombIntra": [(0, 0, 0, 1, 1, 0.7, 0.0), (0, 0, 0, 2, 2, 0.4, 0.0)],
                          "Hund": [(0, 0, 0, 1, 2, 0.2, 0.0), (0, 0, 0, 2, 1, 0.2, 0.0)],
                          "CoulombInter": [(1, 0, 0, 1, 2, 0.3, 0.0), (-1, 0, 0, 2, 1, 0.3, 0.0)]})
        cases = [
            ("norb=2, on-site + off-site", f2, 16, 4),
            ("norb=1, off-site only",
             _factors_norb1({"CoulombInter": rows_v}, shape=(16, 16, 1)), 256, 1),
            ("norb=1, off-site + U",
             _factors_norb1({"CoulombInter": rows_v, "CoulombIntra": intra1},
                            shape=(16, 16, 1)), 256, 1),
        ]

        def _peak(f, nmat, nvol, nd):
            cb = _random_chibar(nmat, nvol, nd, 7)
            ref = dense_w2(cb, f)                              # self-allocating call
            out = np.zeros_like(cb)
            T1, T2 = np.empty_like(cb), np.empty_like(cb)
            accumulate_batch(out, cb, 0, f, work=(T1, T2))      # warm-up
            out[...] = 0.0
            tracing = tracemalloc.is_tracing()
            if not tracing:
                tracemalloc.start()
            tracemalloc.reset_peak()
            base = tracemalloc.get_traced_memory()[0]
            accumulate_batch(out, cb, 0, f, work=(T1, T2))
            peak = tracemalloc.get_traced_memory()[1]
            if not tracing:
                tracemalloc.stop()
            self.assertGreater(np.abs(ref).max(), 1e-6)         # anti-vacuity
            np.testing.assert_array_equal(out, ref)             # work must not change it
            return peak - base, cb.nbytes, cb.nbytes // 16

        for name, f, nvol, nd in cases:
            with self.subTest(case=name):
                small, batch_s, mask_s = _peak(f, 128, nvol, nd)
                big, batch_b, mask_b = _peak(f, 512, nvol, nd)
                print("\naccumulate_batch [{}]: batch {} B peak extra {} B (mask {} B); "
                      "batch {} B peak extra {} B (mask {} B); growth {} B"
                      .format(name, batch_s, small, mask_s, batch_b, big, mask_b, big - small))
                self.assertLess(small, mask_s + _ALLOW)
                self.assertLess(big - small, (mask_b - mask_s) + (1 << 16))

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
