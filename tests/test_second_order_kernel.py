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

    with tempfile.TemporaryDirectory() as d:
        for f in ("geom.dat", "transfer.dat"):
            shutil.copy(os.path.join("tests/rpa/input", f), d)
        idict = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
        for t, rows in rows_by_type.items():
            _write_wan(os.path.join(d, t.lower() + ".dat"), t, 1, rows)
            idict[t] = t.lower() + ".dat"
        r = read_input_k.QLMSkInput({"path_to_input": d, "interaction": idict})
        par = {"T": 2.0, "filling": 0.5, "CellShape": list(shape), "SubShape": [1, 1, 1],
               "Nmat": 8, "IterationMax": 1, "Mix": 1.0, "EPS": 1,
               "flex_second_order": "takimoto"}
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

    def test_every_matmul_output_is_contiguous(self):
        """Portability of the ``out=`` targets: the natural sub-block of a
        scratch buffer (``T2[:, :, :, :norb]`` and friends) is STRIDED, which
        numpy accepts as ``matmul(out=)`` but another array module need not,
        and the GPU path does reach this kernel. Every ``out=`` the kernel
        passes must therefore be C-contiguous -- asserted here by wrapping the
        array module's ``matmul`` for the duration of one call -- while the
        density-slot gathers and scatters stay strided views (no ``out=``).

        A caller that lends NON-contiguous buffers of the right shape is not
        an error: the kernel ignores the loan and allocates its own, and the
        result must be the same."""
        from unittest import mock
        from hwave.solver import backend as _bk
        from hwave.solver.second_order import accumulate_batch, dense_w2
        s, f = _factors({"CoulombIntra": [(0, 0, 0, 1, 1, 0.7, 0.0), (0, 0, 0, 2, 2, 0.4, 0.0)],
                         "Hund": [(0, 0, 0, 1, 2, 0.2, 0.0), (0, 0, 0, 2, 1, 0.2, 0.0)],
                         "CoulombInter": [(1, 0, 0, 1, 2, 0.3, 0.0), (-1, 0, 0, 2, 1, 0.3, 0.0)]})
        self.assertGreater(f.norb, 1)          # nd > norb: the sub-blocks are strided
        self.assertIsNotNone(f.vpair)          # the off-site branch runs
        cb = _random_chibar(8, 16, 4, 5)
        ref = dense_w2(cb, f)
        self.assertGreater(np.abs(ref).max(), 1e-6)            # anti-vacuity

        seen = []
        real = np.matmul

        def checking_matmul(a, b, out=None, **kw):
            if out is not None:
                seen.append(bool(out.flags.c_contiguous))
                self.assertTrue(out.flags.c_contiguous,
                                "matmul out= is not contiguous: shape {} strides {}"
                                .format(out.shape, out.strides))
            return real(a, b, out=out, **kw) if out is not None else real(a, b, **kw)

        class _Proxy:
            matmul = staticmethod(checking_matmul)

            def __getattr__(self, name):
                return getattr(np, name)

        got = np.zeros_like(cb)
        T1, T2 = np.empty_like(cb), np.empty_like(cb)
        with mock.patch.object(_bk, "array_module_of", return_value=_Proxy()):
            accumulate_batch(got, cb, 0, f, work=(T1, T2))
        self.assertGreater(len(seen), 4)                       # the wrapper really ran
        self.assertTrue(all(seen))
        np.testing.assert_array_equal(got, ref)

        # a non-contiguous loan of the right shape: accepted, own buffers used
        big = np.empty((8, 16, 8, 4), complex)
        nc1, nc2 = big[:, :, ::2, :], big[:, :, 1::2, :]
        self.assertFalse(nc1.flags.c_contiguous)
        got2 = np.zeros_like(cb)
        with mock.patch.object(_bk, "array_module_of", return_value=_Proxy()):
            accumulate_batch(got2, cb, 0, f, work=(nc1, nc2))
        np.testing.assert_array_equal(got2, ref)
        # a mis-shaped loan is still refused
        with self.assertRaises(ValueError):
            accumulate_batch(np.zeros_like(cb), cb, 0, f,
                             work=(np.empty((7, 16, 4, 4), complex), np.empty_like(cb)))

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



class TestGuards(unittest.TestCase):
    """Preconditions of :func:`accumulate_batch`.

    Every one of these is a SILENT-corruption route rather than a crash:
    an array whose pair axes do not match the factor pack contracts the
    wrong slots, a mismatched output broadcasts instead of accumulating,
    and a lent buffer that overlaps an operand is overwritten while that
    operand is still being read."""

    def _pack(self):
        return _factors({"CoulombIntra": [(0, 0, 0, 1, 1, 0.7, 0.0), (0, 0, 0, 2, 2, 0.4, 0.0)],
                         "CoulombInter": [(1, 0, 0, 1, 2, 0.3, 0.0), (-1, 0, 0, 2, 1, 0.3, 0.0)]})[1]

    def test_shape_mismatches_are_named(self):
        from hwave.solver.second_order import accumulate_batch
        f = self._pack()
        cb = _random_chibar(8, 16, 4, 11)
        # out_b of a different batch length
        with self.assertRaises(ValueError) as cm:
            accumulate_batch(np.zeros((4, 16, 4, 4), complex), cb, 0, f)
        self.assertIn("out_b", str(cm.exception))
        # pair axes that do not match the factor pack
        for bad in ((8, 16, 2, 2), (8, 16, 4, 2)):
            with self.assertRaises(ValueError) as cm:
                accumulate_batch(np.zeros(bad, complex), np.zeros(bad, complex), 0, f)
            self.assertIn("chibar_b", str(cm.exception))
            self.assertIn("(nb, 16, 4, 4)", str(cm.exception))
        # a volume axis that does not match, with an off-site pack
        self.assertIsNotNone(f.vpair)
        with self.assertRaises(ValueError) as cm:
            accumulate_batch(np.zeros((8, 9, 4, 4), complex), np.zeros((8, 9, 4, 4), complex), 0, f)
        self.assertIn("(nb, 16, 4, 4)", str(cm.exception))
        # rank, not just extents
        with self.assertRaises(ValueError):
            accumulate_batch(np.zeros((16, 4, 4), complex), np.zeros((16, 4, 4), complex), 0, f)

    def test_shape_validation_runs_with_no_nonzero_triple(self):
        """A pack with nothing to contract still validates: otherwise the
        guard would be silently disabled exactly where the caller gets no
        other signal that the arrays are wrong."""
        from hwave.solver.second_order import accumulate_batch, build_factors
        s, split = _split_for({})
        f = build_factors(split, s.lattice, 2)
        self.assertEqual(f.triples, ())
        self.assertIsNone(f.vpair)
        with self.assertRaises(ValueError) as cm:
            accumulate_batch(np.zeros((8, 16, 3, 3), complex),
                             np.zeros((8, 16, 3, 3), complex), 0, f)
        self.assertIn("chibar_b", str(cm.exception))

    def test_empty_and_zero_packs_leave_the_output_untouched(self):
        """An accumulator adds; with nothing to add it must change nothing --
        for a pack with no factors at all AND for one whose off-site vertex
        is declared but zero."""
        from hwave.solver.second_order import accumulate_batch, build_factors
        cb = _random_chibar(8, 16, 4, 12)
        empty_s, empty_split = _split_for({})
        zero_s, zero_split = _split_for(
            {"CoulombInter": [(1, 0, 0, 1, 1, 0.0, 0.0), (-1, 0, 0, 1, 1, 0.0, 0.0)]})
        packs = [("no factors", build_factors(empty_split, empty_s.lattice, 2), False),
                 ("zero off-site rows", build_factors(zero_split, zero_s.lattice, 2), True)]
        for name, f, has_vpair in packs:
            with self.subTest(pack=name):
                self.assertEqual(f.triples, ())
                self.assertEqual(f.vpair is not None, has_vpair)
                prefilled = np.full((8, 16, 4, 4), 3.0 - 2.0j)
                out = prefilled.copy()
                accumulate_batch(out, cb, 0, f)
                np.testing.assert_array_equal(out, prefilled)

    def test_low_rank_inputs_are_named_not_indexed(self):
        """A rank-0/1/2 argument must reach the NAMED refusal, not an
        ``IndexError`` from the kernel's own ``chibar_b.shape[2]``: the
        validation runs before anything indexes an axis the array need not
        have."""
        from hwave.solver.second_order import accumulate_batch
        f = self._pack()
        good = _random_chibar(8, 16, 4, 16)
        for shape in ((), (4,), (16, 4)):
            bad = np.zeros(shape, complex)
            with self.subTest(shape=shape, which="chibar_b"):
                with self.assertRaises(ValueError) as cm:
                    accumulate_batch(np.zeros_like(bad), bad, 0, f)
                self.assertIn("chibar_b", str(cm.exception))
            with self.subTest(shape=shape, which="out_b"):
                with self.assertRaises(ValueError) as cm:
                    accumulate_batch(bad, good, 0, f)
                self.assertIn("out_b", str(cm.exception))

    def test_byte_extent_spans_the_region_a_view_reaches_across(self):
        """``_byte_extent`` is the device alias test's notion of "occupied",
        and it is NOT ``nbytes``.

        The arithmetic is array-module agnostic (shape and strides only), so
        it is exercised here on numpy stand-ins for the device views it is
        actually written for. The case that matters is the strided one: a
        view with half its base's ``nbytes`` still reaches across nearly the
        whole base, and an overlap test built on ``nbytes`` would call it
        disjoint from a buffer sitting in the gap it strides over."""
        from hwave.solver.second_order import _byte_extent
        base = np.zeros((4, 6), dtype=np.complex128)        # itemsize 16
        self.assertEqual(_byte_extent(base), (0, base.nbytes))
        strided = base[:, ::2]
        self.assertEqual(strided.nbytes, base.nbytes // 2)
        # spans from the first element to the last one it reaches
        self.assertEqual(_byte_extent(strided), (0, (3 * 6 + 4) * 16 + 16))
        self.assertGreater(_byte_extent(strided)[1], strided.nbytes)
        # a reversed view has a NEGATIVE low offset relative to its own base
        # pointer (which numpy places at the view's first element)
        rev = base[::-1]
        self.assertEqual(_byte_extent(rev), (-3 * 6 * 16, 16 * 6))
        # one row of the middle: offset is carried by the POINTER, not here
        row = base[2]
        self.assertEqual(_byte_extent(row), (0, 6 * 16))
        # an empty array spans nothing
        self.assertEqual(_byte_extent(base[:0]), (0, 0))
        self.assertEqual(_byte_extent(np.zeros((0, 3), complex)), (0, 0))

    def test_device_alias_rule_uses_the_extent_not_nbytes(self):
        """The device branch of ``_shares_storage`` on a numpy stand-in that
        exposes a CuPy-style ``.data.ptr``: two views that interleave inside
        one buffer must be reported as possibly sharing, which an
        ``nbytes``-based span would miss.

        ``_shares_storage`` short-circuits on real ``numpy.ndarray`` pairs
        (it can answer exactly there), so the stand-in below is what routes
        the call through the pointer-extent branch the device path takes."""
        from hwave.solver.second_order import _shares_storage, _byte_extent

        class _Ptr(object):
            def __init__(self, ptr):
                self.ptr = ptr

        class _DeviceLike(object):
            """Everything the device branch reads: a base pointer, a shape,
            strides and an itemsize."""

            def __init__(self, arr, base):
                self.shape, self.strides, self.itemsize = arr.shape, arr.strides, arr.itemsize
                self.nbytes = arr.nbytes
                off = arr.__array_interface__["data"][0] - base.__array_interface__["data"][0]
                self.data = _Ptr(1 << 20)                  # an arbitrary device address
                self.data.ptr += off

        base = np.zeros((4, 6), dtype=np.complex128)
        # a widely strided view (columns 0 and 5) and a small buffer sitting
        # INSIDE the region it strides over
        wide, inner = base[:, ::5], base[1, 2:3]
        a, b = _DeviceLike(wide, base), _DeviceLike(inner, base)
        # the nbytes rule calls them disjoint ...
        self.assertTrue(a.data.ptr + a.nbytes <= b.data.ptr
                        or b.data.ptr + b.nbytes <= a.data.ptr)
        # ... and the extent rule refuses, which is the safe direction for a
        # guard that cannot do element-level analysis on a device array
        self.assertTrue(_shares_storage(a, b))
        # a pair that genuinely shares elements is caught either way
        self.assertTrue(np.shares_memory(base, base[1]))
        self.assertTrue(_shares_storage(_DeviceLike(base, base),
                                        _DeviceLike(base[1], base)))
        # genuinely separate buffers are still reported disjoint
        other = np.zeros((4, 6), dtype=np.complex128)
        far = _DeviceLike(other, other)
        far.data.ptr = 1 << 24
        self.assertFalse(_shares_storage(_DeviceLike(base, base), far))
        # and an empty operand shares nothing
        empty = _DeviceLike(base[:0], base)
        self.assertEqual(_byte_extent(empty), (0, 0))
        self.assertFalse(_shares_storage(_DeviceLike(base, base), empty))
        # a module that exposes no base pointer cannot say
        class _Opaque(object):
            shape, strides, itemsize, nbytes = (1,), (16,), 16, 16
        self.assertIsNone(_shares_storage(_Opaque(), _Opaque()))

    def test_output_aliasing_the_bubble_is_refused(self):
        """``out_b`` is ACCUMULATED into while ``chibar_b`` is still being
        read (every spin triple reads the bubble again, and the off-site
        terms read its density slots afterwards), so the two sharing storage
        corrupts the result exactly as a lent buffer overlapping an operand
        does."""
        from hwave.solver.second_order import accumulate_batch
        f = self._pack()
        cb = _random_chibar(8, 16, 4, 17)
        for name, out in (("the same array", cb), ("a view of it", cb[...])):
            with self.subTest(case=name):
                with self.assertRaises(ValueError) as cm:
                    accumulate_batch(out, cb, 0, f)
                self.assertIn("out_b shares memory with chibar_b", str(cm.exception))
        # a separate output of the same shape is of course fine
        accumulate_batch(np.zeros_like(cb), cb, 0, f)

    def test_aliasing_of_the_lent_buffers_is_refused(self):
        from hwave.solver.second_order import accumulate_batch
        f = self._pack()
        cb = _random_chibar(8, 16, 4, 13)
        out = np.zeros_like(cb)
        T = np.empty_like(cb)
        # the same object twice
        with self.assertRaises(ValueError) as cm:
            accumulate_batch(out, cb, 0, f, work=(T, T))
        self.assertIn("same array", str(cm.exception))
        # two views of one buffer
        big = np.empty((8, 16, 4, 4), complex)
        with self.assertRaises(ValueError) as cm:
            accumulate_batch(out, cb, 0, f, work=(big[...], big[...]))
        self.assertIn("share memory", str(cm.exception))
        # a buffer that IS an operand, or a view of one
        for bad, which in (((cb, np.empty_like(cb)), "chibar_b"),
                           ((np.empty_like(cb), out), "out_b"),
                           ((out[...], np.empty_like(cb)), "out_b")):
            with self.assertRaises(ValueError) as cm:
                accumulate_batch(out, cb, 0, f, work=bad)
            self.assertIn("shares memory with {}".format(which), str(cm.exception))

    def test_non_contiguous_operands(self):
        """``out_b`` and ``chibar_b`` may be strided views -- the bond path
        passes exactly that (``W_b[:, :, :nd, :nd]``). Only the LENT buffers
        need contiguity, and a non-contiguous loan is discarded rather than
        refused."""
        from hwave.solver.second_order import accumulate_batch, dense_w2
        rows_v = [(1, 0, 0, 1, 2, 0.3, 0.0), (-1, 0, 0, 2, 1, 0.3, 0.0)]
        packs = [("on-site only",
                  _factors({"CoulombIntra": [(0, 0, 0, 1, 1, 0.7, 0.0),
                                             (0, 0, 0, 2, 2, 0.4, 0.0)]})[1]),
                 ("with off-site",
                  _factors({"CoulombIntra": [(0, 0, 0, 1, 1, 0.7, 0.0)],
                            "CoulombInter": rows_v})[1])]
        for name, f in packs:
            with self.subTest(pack=name):
                cb = _random_chibar(8, 16, 4, 14)
                ref = dense_w2(cb, f)
                self.assertGreater(np.abs(ref).max(), 1e-6)            # anti-vacuity
                # strided views of larger arrays, on BOTH operands
                cb_big = np.zeros((8, 16, 8, 8), complex)
                cb_big[:, :, :4, :4] = cb
                cb_view = cb_big[:, :, :4, :4]
                out_big = np.zeros((8, 16, 8, 8), complex)
                out_view = out_big[:, :, :4, :4]
                self.assertFalse(cb_view.flags.c_contiguous)
                self.assertFalse(out_view.flags.c_contiguous)
                accumulate_batch(out_view, cb_view, 0, f)
                np.testing.assert_allclose(out_view, ref, rtol=0, atol=1e-13)
                # and nothing outside the block was written
                out_big[:, :, :4, :4] = 0.0
                self.assertEqual(np.abs(out_big).max(), 0.0)

    def test_norb1_offsite_and_mixed_identities(self):
        """``norb = 1`` (nd = 1: the density slot IS the whole pair axis) --
        the two closed-form checks of spec 2.4 written directly on the
        kernel's output, as the ``norb = 2`` test does on its (00) block."""
        from hwave.solver.second_order import dense_w2
        rows_v = [(1, 0, 0, 1, 1, 0.3, 0.0), (-1, 0, 0, 1, 1, 0.3, 0.0)]
        intra = [(0, 0, 0, 1, 1, 0.5, 0.0)]
        fv = _factors_norb1({"CoulombInter": rows_v})
        fu = _factors_norb1({"CoulombIntra": intra})
        fuv = _factors_norb1({"CoulombInter": rows_v, "CoulombIntra": intra})
        nvol = fv.nvol
        cb = _random_chibar(8, nvol, 1, 15)
        vq = fv.vpair[0, 0, :, 0, 0]
        self.assertGreater(np.abs(vq).max(), 1e-3)                     # anti-vacuity
        Wv, Wu, Wuv = dense_w2(cb, fv), dense_w2(cb, fu), dense_w2(cb, fuv)
        # off-site only: 2 V(q)^2 chibar (the direct skeleton with its spin sum)
        np.testing.assert_allclose(Wv[:, :, 0, 0], 2.0 * vq[None] ** 2 * cb[:, :, 0, 0],
                                   rtol=0, atol=1e-13)
        # on-site only: U^2 chibar
        np.testing.assert_allclose(Wu[:, :, 0, 0], 0.25 * cb[:, :, 0, 0], rtol=0, atol=1e-13)
        # mixed: 2 U V(q) chibar
        np.testing.assert_allclose((Wuv - Wu - Wv)[:, :, 0, 0],
                                   2.0 * 0.5 * vq[None] * cb[:, :, 0, 0], rtol=0, atol=1e-13)


if __name__ == "__main__":
    unittest.main()
