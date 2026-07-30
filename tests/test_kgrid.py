"""The single-sourced FFT-grid reversal map (#108).

Nineteen call sites across sc.py, rpa.py, flex.py, eliashberg_dynamic.py
and uhfk.py used to spell the map ``i -> (N - i) % N`` three different
ways; these tests pin the map's algebra once and prove all historical
spellings are the same gather, so the consolidation is bit-identical by
construction.
"""

import unittest

import numpy as np
import pytest

from hwave.solver.kgrid import reverse_fft_axes


class TestReverseFftAxes(unittest.TestCase):

    def test_matches_the_index_formula(self):
        """out[i] == arr[(N - i) % N] on every listed axis, odd and even N."""
        rng = np.random.default_rng(3)
        for shape, axes in [((5,), (0,)), ((4,), (0,)),
                            ((3, 4, 5), (0, 1, 2)),
                            ((2, 3, 4, 5, 6), (2, 3, 4))]:
            with self.subTest(shape=shape, axes=axes):
                a = rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
                got = reverse_fft_axes(a, axes)
                want = a
                for ax in axes:
                    n = a.shape[ax]
                    want = np.take(want, (-np.arange(n)) % n, axis=ax)
                np.testing.assert_array_equal(got, want)

    def test_involution_and_gamma_fixed_point(self):
        rng = np.random.default_rng(4)
        a = rng.standard_normal((3, 5, 4, 7)) * 1j + rng.standard_normal((3, 5, 4, 7))
        r = reverse_fft_axes(a, (1, 2, 3))
        np.testing.assert_array_equal(reverse_fft_axes(r, (1, 2, 3)), a)
        # index 0 on every reversed axis (Gamma / origin) stays put
        np.testing.assert_array_equal(r[:, 0, 0, 0], a[:, 0, 0, 0])

    def test_all_three_historical_spellings_are_this_gather(self):
        """roll(flip), flip(roll(-1)) and the fancy-index gather were the
        three spellings in the tree; all are bit-identical to the shared
        implementation."""
        rng = np.random.default_rng(5)
        a = rng.standard_normal((4, 3, 5, 2)) + 1j * rng.standard_normal((4, 3, 5, 2))
        axes = (0, 1, 2)
        got = reverse_fft_axes(a, axes)
        # spelling 1: reverse then roll +1 (sc.py, eliashberg_dynamic.py)
        s1 = np.roll(a[::-1, ::-1, ::-1, :], (1, 1, 1), axes)
        # spelling 2: roll -1 then flip (rpa.py, flex.py, uhfk.py)
        s2 = np.flip(np.roll(a, -1, axis=axes), axis=axes)
        # spelling 3: fancy-index gather (sc._symmetrise_interactions_k)
        s3 = a[(-np.arange(4)) % 4][:, (-np.arange(3)) % 3][:, :, (-np.arange(5)) % 5]
        for i, s in enumerate((s1, s2, s3), 1):
            with self.subTest(spelling=i):
                np.testing.assert_array_equal(got, s)

    def test_short_axis_identity_trap(self):
        """For N <= 2 the map is the identity: (N - i) % N == i for
        i in {0, 1}. This is WHY a misapplied reversal is invisible on
        one- and two-site fixtures (the transverse-channel defect
        history) -- pinned here so the trap is executable knowledge."""
        rng = np.random.default_rng(6)
        for n in (1, 2):
            a = rng.standard_normal((n, 5)) + 0j
            np.testing.assert_array_equal(reverse_fft_axes(a, (0,)), a)
        a3 = rng.standard_normal((3, 5)) + 0j
        self.assertFalse(np.array_equal(reverse_fft_axes(a3, (0,)), a3))

    def test_composed_site_shapes_are_bit_identical_to_the_old_code(self):
        """The two composite forms that fold in a frequency/tau axis:
        sc._calc_g2's centered-Matsubara flip + spatial map, and flex's
        signed tau-flip + spatial map. Both must reproduce the exact old
        inline expressions."""
        rng = np.random.default_rng(7)
        # sc._calc_g2 / eliashberg_dynamic: (norb, norb, Nx, Ny, Nz, nmat)
        g = rng.standard_normal((2, 2, 4, 3, 2, 6)) \
            + 1j * rng.standard_normal((2, 2, 4, 3, 2, 6))
        old = np.roll(g[:, :, ::-1, ::-1, ::-1, ::-1], (1, 1, 1), (2, 3, 4))
        new = reverse_fft_axes(g[..., ::-1], (2, 3, 4))
        np.testing.assert_array_equal(new, old)
        # flex chi0 kernel: (nblock, ntau, nx, ny, nz, nd*nd)
        gt = rng.standard_normal((1, 6, 4, 3, 2, 4)) \
            + 1j * rng.standard_normal((1, 6, 4, 3, 2, 4))
        old = -np.flip(np.roll(gt, -1, axis=(2, 3, 4)), axis=(1, 2, 3, 4))
        new = -np.flip(reverse_fft_axes(gt, (2, 3, 4)), axis=1)
        np.testing.assert_array_equal(new, old)


class TestProductionCallSites(unittest.TestCase):
    """Call the cheap converted PRODUCTION functions directly against the
    exact pre-consolidation inline expressions, so a wrong axis tuple at a
    call site cannot hide behind the helper-level tests. (The heavy kernel
    sites -- the RPA/FLEX chi0 kernels and UHFk's interaction setup -- are
    covered end-to-end by their reference-value suites.)"""

    def test_sc_reverse_k_and_orbital(self):
        import hwave.sc as sc

        rng = np.random.default_rng(11)
        gap = rng.standard_normal((3, 3, 4, 3, 5)) \
            + 1j * rng.standard_normal((3, 3, 4, 3, 5))
        old = np.swapaxes(
            np.roll(gap[:, :, ::-1, ::-1, ::-1], 1, axis=(2, 3, 4)), 0, 1)
        np.testing.assert_array_equal(sc._reverse_k_and_orbital(gap), old)

    def test_eliashberg_dynamic_reverse_kw_and_orbital(self):
        from hwave.solver import eliashberg_dynamic as ed

        rng = np.random.default_rng(12)
        gap_w = rng.standard_normal((3, 3, 4, 3, 5, 6)) \
            + 1j * rng.standard_normal((3, 3, 4, 3, 5, 6))
        old = np.roll(gap_w[:, :, ::-1, ::-1, ::-1, :], 1, axis=(2, 3, 4))
        old = np.swapaxes(old[..., ::-1], 0, 1)
        np.testing.assert_array_equal(ed._reverse_kw_and_orbital(gap_w), old)

    def test_rpa_assemble_transverse_vertex(self):
        """The transverse vertex assembly's q -> -q used to be a
        coordinate-wise flat-index formula (the nineteenth spelling);
        pin the production method against a replica of that exact old
        gather, on a lattice with every dimension > 2 so the map is not
        the identity anywhere."""
        import types
        import hwave.solver.rpa as rpa_module

        rng = np.random.default_rng(14)
        norb, ns = 2, 2
        nx, ny, nz = 4, 3, 5
        nvol, nd = nx * ny * nz, norb * ns

        stub = object.__new__(rpa_module.RPA)
        stub.norb, stub.ns = norb, ns
        stub.lattice = types.SimpleNamespace(nvol=nvol, shape=(nx, ny, nz))

        ham_orig = rng.standard_normal((nvol, nd, nd, nd, nd)) \
            + 1j * rng.standard_normal((nvol, nd, nd, nd, nd))
        got = stub._assemble_transverse_vertex(ham_orig)

        # replica of the pre-consolidation assembly (coordinate-wise qrev)
        ham_4d = ham_orig.reshape(nvol, nd, nd, nd, nd)
        cross_block = ham_4d[:, norb:, norb:, :norb, :norb]
        spin_flip_block = ham_4d[:, :norb, norb:, :norb, norb:]
        pair_swap = (0, 3, 4, 1, 2)
        oi = np.arange(norb)
        palin = ((oi[:, None, None, None] == oi[None, :, None, None])
                 & (oi[None, None, :, None] == oi[None, None, None, :]))[None]
        _ix = (-np.arange(nx)) % nx
        _iy = (-np.arange(ny)) % ny
        _iz = (-np.arange(nz)) % nz
        qrev = ((_ix[:, None, None] * ny + _iy[None, :, None]) * nz
                + _iz[None, None, :]).reshape(-1)

        def _mean_old(block):
            swapped = block.transpose(*pair_swap)
            return 0.5 * (block
                          + np.where(palin, swapped[qrev], np.conj(swapped)))

        want = (-_mean_old(cross_block).transpose(0, 2, 4, 1, 3)
                + _mean_old(spin_flip_block))
        np.testing.assert_array_equal(got, want)

    def test_sc_symmetrise_interactions_k(self):
        import hwave.sc as sc

        rng = np.random.default_rng(13)
        nx, ny, nz = 4, 3, 5
        M = rng.standard_normal((3, 3, nx, ny, nz)) \
            + 1j * rng.standard_normal((3, 3, nx, ny, nz))
        Mrev_old = M[:, :, (-np.arange(nx)) % nx][:, :, :, (-np.arange(ny)) % ny][
            :, :, :, :, (-np.arange(nz)) % nz]
        want = 0.5 * (M + Mrev_old.transpose(1, 0, 2, 3, 4))
        got = sc._symmetrise_interactions_k({"Hund": M})["Hund"]
        np.testing.assert_array_equal(got, want)


class TestCupyPath(unittest.TestCase):
    """Device-array coverage; skips without a usable CUDA device.

    The numpy path cannot validate CuPy's own non-contiguous reshape
    inside the transverse assembly's unflatten -> map -> flatten round
    trip, so exercise the helper and that round trip on device."""

    def test_reverse_fft_axes_on_device(self):
        cupy = pytest.importorskip("cupy")
        try:
            cupy.zeros(1)
        except Exception:
            self.skipTest("cupy installed but no usable CUDA device")

        rng = np.random.default_rng(15)
        a = rng.standard_normal((60, 2, 2, 2, 2)) \
            + 1j * rng.standard_normal((60, 2, 2, 2, 2))
        d = cupy.asarray(a).transpose(0, 3, 4, 1, 2)   # non-contiguous
        h = a.transpose(0, 3, 4, 1, 2)
        got = reverse_fft_axes(
            d.reshape(4, 3, 5, *d.shape[1:]), (0, 1, 2)).reshape(d.shape)
        self.assertIsInstance(got, cupy.ndarray)
        want = reverse_fft_axes(
            h.reshape(4, 3, 5, *h.shape[1:]), (0, 1, 2)).reshape(h.shape)
        np.testing.assert_array_equal(cupy.asnumpy(got), want)


if __name__ == "__main__":
    unittest.main()
