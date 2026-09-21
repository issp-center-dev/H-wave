"""CPU vs GPU equivalence of the bond-resolved FLEX kernels, and one
end-to-end run of the bond gate on both backends (spec 2026-09-17 section 5).

Nothing here imports cupy at module import time: the whole module is
collected on a CPU-only machine and every case is SKIPPED from
``setUpClass`` (``unittest.SkipTest``) when
:func:`hwave.solver.backend.gpu_available` says there is no usable CUDA
device. That is a unittest-level skip on purpose -- ``python -m unittest
discover`` is the gating runner, where a module-level ``pytest.skip``
would degrade to a collection ERROR.

Tolerances. The two backends run the SAME algorithm on different BLAS /
LAPACK / FFT implementations, so the results agree to floating-point
round-off rather than bit for bit. The dressing tolerances are therefore
scaled by the condition number of the fixture's own solve matrix
``1 -/+ chibar V`` (:func:`_cond_max`) -- the factor by which either
backend's round-off may be amplified by the solve -- and the solve
residual of BOTH backends is checked against the same scaled bound, which
is an absolute statement about each result rather than a comparison of
the two. Set ``HWAVE_GPU_EQUIV_REPORT=1`` to print the measured maxima of
every comparison (off in a normal run, which stays quiet).
"""
import contextlib
import os
import tempfile
import types
import unittest
from unittest import mock

import numpy as np

from tests.heavy_tests import heavy
from tests.test_flex_bond_dressing import _problem

#: Print the measured deviations. Opt-in so the normal run is quiet.
_REPORT = os.environ.get("HWAVE_GPU_EQUIV_REPORT", "").strip().lower() \
    not in ("", "0", "false", "no", "off")


def _cond_max(mat):
    """Worst 2-norm condition number over the (l, q) blocks of ``mat``."""
    sv = np.linalg.svd(mat.reshape(-1, mat.shape[-1], mat.shape[-1]), compute_uv=False)
    return float(np.max(sv[:, 0] / sv[:, -1]))


def _deviation(label, out, ref):
    """Max absolute and max relative deviation of ``out`` from ``ref``,
    reported when ``HWAVE_GPU_EQUIV_REPORT`` asks for it."""
    out = np.asarray(out)
    ref = np.asarray(ref)
    amax = float(np.max(np.abs(out - ref))) if ref.size else 0.0
    scale = float(np.max(np.abs(ref))) if ref.size else 0.0
    rmax = amax / scale if scale > 0.0 else 0.0
    if _REPORT:
        print("[gpu-equivalence] {}: max|gpu - cpu| = {:.3e}, max|cpu| = {:.3e}, "
              "ratio = {:.3e}".format(label, amax, scale, rmax))
    return amax, rmax


class _GpuCase(unittest.TestCase):
    """Base: skip the whole class without a usable device, and import cupy
    only once that check has passed."""

    @classmethod
    def setUpClass(cls):
        from hwave.solver import backend
        if not backend.gpu_available():
            raise unittest.SkipTest("no usable CUDA device / cupy")
        import cupy
        cls.cupy = cupy


class TestDressBatchEquivalence(_GpuCase):
    """``bond_channels.dress_batch`` (spec 4.1) on numpy and on cupy."""

    def _check(self, chi_bar, V, channel, guard, label, cond_tol=None):
        from hwave.solver import bond_channels as bc
        nmat, nvol, ND = chi_bar.shape[0], chi_bar.shape[1], chi_bar.shape[-1]
        tol = {} if cond_tol is None else {"cond_tol": cond_tol}
        ref, c_ref = bc.dress_batch(chi_bar, V, channel, l0=0, nmat=nmat,
                                    spatial_shape=(nvol, 1, 1), guard_freqs=guard, **tol)
        out_d, c_out = bc.dress_batch(self.cupy.asarray(chi_bar), V, channel, l0=0, nmat=nmat,
                                      spatial_shape=(nvol, 1, 1), guard_freqs=guard, **tol)
        # the kernel returns the module of its input, never a host copy
        self.assertEqual(type(out_d).__module__.split(".")[0], "cupy")
        out = self.cupy.asnumpy(out_d)
        sign = 1.0 if channel == "spin" else -1.0
        mat = np.eye(ND)[None, None] - sign * (chi_bar @ V[None])
        cm = _cond_max(mat)
        _deviation("{} ({}, guard={})".format(label, channel, guard), out, ref)
        np.testing.assert_allclose(out, ref, rtol=1e-10 * cm,
                                   atol=1e-12 * np.max(np.abs(chi_bar)))
        # absolute statement about each backend's own solve, not a comparison
        for name, chi in (("cpu", ref), ("gpu", out)):
            r = float(np.max(bc.solve_residual(mat, chi, chi_bar)))
            if _REPORT:
                print("[gpu-equivalence] {} ({}, guard={}) {} solve residual = "
                      "{:.3e} (bound {:.3e})".format(label, channel, guard, name,
                                                     r, 1e-11 * cm))
            self.assertLessEqual(r, 1e-11 * cm)
        if c_ref is not None:
            self.assertIsNotNone(c_out)
            self.assertAlmostEqual(c_out, c_ref, delta=1e-10 * abs(c_ref))

    def test_random_general_complex(self):
        chi_bar, S, C = _problem(nmat=6, nvol=4, nd=4, B=3, seed=0)
        for ch, V in (("spin", S), ("charge", C)):
            for guard in ("all", "static"):
                with self.subTest(channel=ch, guard=guard):
                    self._check(chi_bar, V, ch, guard, "random general")

    def test_ill_conditioned_but_accepted(self):
        """A GENUINELY ill-conditioned denominator (condition number 1e5)
        that both backends still solve.

        The denominator is prescribed directly: with a q-independent vertex
        ``V`` and ``chibar = (1 - M) V^-1`` the kernel's solve matrix is
        exactly ``M``, so choosing ``M = U diag(1, ..., 1e-5) U^H`` with a
        random unitary ``U`` fixes its condition number at 1e5 -- the regime
        where the two backends' LAPACK and cuSOLVER round-off is amplified,
        which is what the ``cond_max``-scaled tolerances exist for. (An
        earlier fixture scaled the denominator uniformly, which leaves the
        condition number at 1 and tested nothing of the kind.)

        The default conditioning floor refuses anything with
        ``sigma_min / sigma_max <= 1e-3``, i.e. every condition number at or
        above 1e3, so the guard floor is lowered here deliberately -- the
        escape hatch the guard's own refusal message names. The guard still
        runs: it reports the 1e-5 minimum, and both backends must agree on
        it.
        """
        chi_bar, S, _C = _problem(nmat=4, nvol=2, nd=2, B=2, seed=7)
        nmat, nvol, ND = chi_bar.shape[0], chi_bar.shape[1], chi_bar.shape[-1]
        V = np.repeat(S[:1], nvol, axis=0)                    # q-independent vertex
        rng = np.random.default_rng(11)
        U, _r = np.linalg.qr(rng.normal(size=(ND, ND)) + 1j * rng.normal(size=(ND, ND)))
        d = np.ones(ND)
        d[-1] = 1.0e-5
        M = (U * d) @ U.conj().T                              # cond(M) = 1e5
        chi_bar = np.broadcast_to((np.eye(ND) - M) @ np.linalg.inv(V[0]),
                                  (nmat, nvol, ND, ND)).copy()
        # spin channel: mat = 1 - chibar V = M, by construction
        cm = _cond_max(np.eye(ND)[None, None] - chi_bar @ V[None])
        self.assertGreaterEqual(cm, 1.0e4)
        for guard in ("all", "static"):
            with self.subTest(guard=guard):
                self._check(chi_bar, V, "spin", guard, "ill-conditioned denominator",
                            cond_tol=1e-8)

    def test_near_zero_values(self):
        chi_bar, S, _C = _problem(nmat=4, nvol=2, nd=2, B=2, seed=8)
        chi_bar = chi_bar * 1e-14
        self._check(chi_bar, S, "spin", "all", "near-zero bubble")

    def test_interval_guard_matches_the_host_svd_guard(self):
        """Guard OUTPUTS of the device interval path equal the production host
        SVD path on the same device stack (issue #197): the minimum, its
        location and the refusal text; the dressed batch keeps the solve
        tolerance of the other cases."""
        from hwave.solver import bond_channels as bc
        chi_bar, S, C = _problem(nmat=6, nvol=8, nd=4, B=3, seed=11)
        nmat, nvol = chi_bar.shape[:2]
        for scale in (1.0, 6.0, 16.0):              # passing, near, refusing
            for channel, W in (("spin", S), ("charge", C)):
                with self.subTest(scale=scale, channel=channel):
                    cb_d = self.cupy.asarray(chi_bar); W_d = self.cupy.asarray(scale * W)
                    outs = []
                    for method in ("svd", "interval"):
                        try:
                            chi, cond = bc.dress_batch(cb_d, W_d, channel, l0=0, nmat=nmat,
                                                       spatial_shape=(nvol, 1, 1), guard_method=method)
                            outs.append(("ok", cond, self.cupy.asnumpy(chi)))
                        except bc.BondConditioningError as e:
                            outs.append(("refused", str(e), e.l, e.iq, e.q, e.worst, e.smin, e.smax))
                    self.assertEqual(outs[0][0], outs[1][0])
                    if outs[0][0] == "ok":
                        self.assertEqual(outs[0][1], outs[1][1])
                        np.testing.assert_array_equal(outs[0][2], outs[1][2])
                    else:
                        self.assertEqual(outs[0][1:], outs[1][1:])

    def test_interval_guard_device_memory_increment(self):
        # bond_conditioning_interval calls _norm_upper twice per call (once
        # on A, once on B after the inverse exists), so the spy below only
        # samples the pool at those two points -- the assertion is an upper
        # bound on the increment AT THOSE POINTS, not a measured peak; the
        # controller records the actual peak from the GPU run.
        from hwave.solver import bond_channels as bc
        chi_bar, S, _C = _problem(nmat=16, nvol=32, nd=4, B=3, seed=12)
        n, ND = 16 * 32, S.shape[-1]
        mat = self.cupy.asarray(np.eye(ND) - chi_bar.reshape(n, ND, ND) @ S[0])
        pool = self.cupy.get_default_memory_pool()
        self.cupy.cuda.Device().synchronize()
        before = pool.used_bytes()
        peak = [before]
        real = bc._norm_upper
        def spy(X, xp):
            peak.append(pool.used_bytes()); return real(X, xp)
        with mock.patch.object(bc, "_norm_upper", spy):
            bc.bond_conditioning_interval(mat, self.cupy)
        self.cupy.cuda.Device().synchronize()
        U_b = n * ND * ND * 16
        self.assertLessEqual(max(peak) - before, 3 * U_b)
        _deviation("guard memory (U_b units)", np.array([(max(peak) - before) / U_b]), np.array([0.0]))


def _bubble_fixture(norb=2, shape=(4, 4, 1), nmat=8, beta=2.0, seed=11):
    """A physical bond fixture (Hermitian k-even hopping, bare Green
    function, a three-channel bond set) plus a NONZERO frequency-constant
    tau-space tail, so that the equal-time endpoint branch of
    ``bubble._prepare_dense`` is exercised on both backends."""
    from tests.eliashberg_bond_fixtures import physical_fixture
    fx = physical_fixture(norb=norb, shape=shape, nmat=nmat, beta=beta, seed=seed)
    rng = np.random.default_rng(seed + 1)
    nvol = fx["nvol"]
    t = rng.normal(size=(nvol, norb, norb)) + 1j * rng.normal(size=(nvol, norb, norb))
    t = 0.25 * (t + np.conj(t).swapaxes(-2, -1))
    tail = np.broadcast_to(t[None, None], fx["green"].shape).copy()
    return fx, np.ascontiguousarray(fx["green"]), tail


class TestBubbleEquivalence(_GpuCase):
    """``flex_bond.assemble_bubble`` on numpy and on cupy (issue #196).

    The bubble kernel is array-module generic, so a device-resident Green
    function makes the whole per-pair pipeline -- the Matsubara
    transform, the spatial FFTs, the reversal and the pair contraction --
    run on the device; only the finished block crosses to the host store.
    This case pins that the two backends produce the same ``chibar``, and
    that the cupy run really took the device rather than falling back to
    a host copy somewhere on the way in."""

    def _assemble(self, xp, fx, green, tail, seen=None):
        from hwave.solver import bubble, flex_bond as fb
        from hwave.solver import backend as bk
        nmat, nvol, nd, ND = fx["nmat"], fx["nvol"], fx["nd"], fx["ND"]
        real_prepare = bubble._prepare_dense

        def _spy(green_kw, *a, **k):
            seen.append(bk.array_module_of(green_kw).__name__)
            return real_prepare(green_kw, *a, **k)
        store = fb.BondBlockStore(nmat, nvol, ND, nd, ("chibar",))
        ctx = fb.BondDeviceContext(xp, fx["S"], fx["C"], green0_tail=tail)
        with store, ctx:
            g = bk.to_device(green, xp)
            self.assertEqual(type(g).__module__.split(".")[0],
                             "cupy" if xp is not np else "numpy")
            patch = (mock.patch.object(bubble, "_prepare_dense", _spy)
                     if seen is not None else contextlib.nullcontext())
            with patch:
                fb.assemble_bubble(store, g, ctx.green0_tail, fx["beta"], fx["view"],
                                   fx["spatial_shape"], 1)
            return {(m, mp): store.get_pair("chibar", m, mp).copy()
                    for m in range(fx["B"]) for mp in range(fx["B"])}

    def test_chibar_matches_block_by_block(self):
        fx, green, tail = _bubble_fixture()
        seen_cpu, seen_gpu = [], []
        ref = self._assemble(np, fx, green, tail, seen_cpu)
        out = self._assemble(self.cupy, fx, green, tail, seen_gpu)
        # the bubble ran where the Green function was, with no host detour
        self.assertEqual(set(seen_cpu), {"numpy"})
        self.assertEqual(set(seen_gpu), {"cupy"})
        self.assertEqual(sorted(out), sorted(ref))
        for key in sorted(ref):
            _deviation("chibar {}".format(key), out[key], ref[key])
            np.testing.assert_allclose(out[key], ref[key], rtol=1e-10, atol=1e-12)

    def test_context_holds_the_tail_on_the_device_and_releases_it(self):
        from hwave.solver import flex_bond as fb
        fx, _green, tail = _bubble_fixture()
        with fb.BondDeviceContext(self.cupy, fx["S"], fx["C"], green0_tail=tail) as dev:
            self.assertEqual(type(dev.green0_tail).__module__.split(".")[0], "cupy")
            np.testing.assert_array_equal(self.cupy.asnumpy(dev.green0_tail), tail)
        self.assertTrue(dev.released)
        with self.assertRaises(RuntimeError):
            dev.green0_tail

    def test_a_host_tail_with_a_device_green_function_is_refused(self):
        from hwave.solver import backend as bk, flex_bond as fb
        fx, green, tail = _bubble_fixture()
        nmat, nvol, nd, ND = fx["nmat"], fx["nvol"], fx["nd"], fx["ND"]
        with fb.BondBlockStore(nmat, nvol, ND, nd, ("chibar",)) as store:
            with self.assertRaises(ValueError) as cm:
                fb.assemble_bubble(store, bk.to_device(green, self.cupy), tail,
                                   fx["beta"], fx["view"], fx["spatial_shape"], 1)
        msg = str(cm.exception)
        self.assertIn("green0_tail", msg)
        self.assertIn("cupy", msg)
        self.assertIn("numpy", msg)


class TestDressAndBuildWEquivalence(_GpuCase):
    """``flex_bond.dress_and_build_w`` (spec 4.2): the effective interaction,
    the collapses and the static slices, over several frequency batch sizes."""

    def test_w_collapses_static_and_batch_sizes(self):
        from hwave.solver import flex_bond as fb
        chi_bar, S, C = _problem(nmat=8, nvol=4, nd=4, B=3, seed=5)
        nmat, nvol, ND, nd = chi_bar.shape[0], chi_bar.shape[1], S.shape[-1], 4
        S_on = np.ascontiguousarray(S[:, :nd, :nd])
        C_on = np.ascontiguousarray(C[:, :nd, :nd])
        view = types.SimpleNamespace(n_channels=ND // nd)

        def run(xp, nb):
            with fb.BondBlockStore(nmat, nvol, ND, nd, ("chibar", "W")) as store, \
                    fb.BondDeviceContext.for_view(xp, S, C, S_on, C_on, view, 2) as dev:
                store.put_freq_batch("chibar", 0, nmat, chi_bar)
                res = fb.dress_and_build_w(store, dev, nb=nb, output_full=False, nmat=nmat,
                                           nvol=nvol, nd=nd, spatial_shape=(nvol, 1, 1),
                                           second_order="takimoto")
                return store.get_freq_batch("W", 0, nmat).copy(), res

        W_ref, r_ref = run(np, nmat)
        mat = np.eye(ND)[None, None] - chi_bar @ S[None]
        cm = _cond_max(mat)
        for nb in sorted({1, 3, nmat // 2, nmat}):
            with self.subTest(nb=nb):
                W, r = run(self.cupy, nb)
                _deviation("W (nb={})".format(nb), W, W_ref)
                _deviation("collapse_s (nb={})".format(nb), r.collapse_s, r_ref.collapse_s)
                _deviation("static_c (nb={})".format(nb), r.static_c, r_ref.static_c)
                np.testing.assert_allclose(W, W_ref, rtol=1e-10 * cm,
                                           atol=1e-12 * np.max(np.abs(W_ref)))
                np.testing.assert_allclose(r.collapse_s, r_ref.collapse_s,
                                           rtol=1e-10 * cm, atol=1e-13)
                np.testing.assert_allclose(r.collapse_c, r_ref.collapse_c,
                                           rtol=1e-10 * cm, atol=1e-13)
                np.testing.assert_allclose(r.static_s, r_ref.static_s,
                                           rtol=1e-10 * cm, atol=1e-13)
                np.testing.assert_allclose(r.static_c, r_ref.static_c,
                                           rtol=1e-10 * cm, atol=1e-13)
                self.assertAlmostEqual(r.cond_min_s, r_ref.cond_min_s, delta=1e-10)
                self.assertAlmostEqual(r.cond_min_c, r_ref.cond_min_c, delta=1e-10)

    def test_output_full_susceptibilities_match(self):
        """``output_full=True`` (issue #198 item 3): the batched dressed
        susceptibilities ``chi_s_w`` / ``chi_c_w`` the store writes alongside
        ``W`` -- otherwise covered on the GPU only through ``output_full``'s
        DressResult fields, never through the store's own full-frequency
        arrays."""
        from hwave.solver import flex_bond as fb
        chi_bar, S, C = _problem(nmat=8, nvol=4, nd=4, B=3, seed=5)
        nmat, nvol, ND, nd = chi_bar.shape[0], chi_bar.shape[1], S.shape[-1], 4
        S_on = np.ascontiguousarray(S[:, :nd, :nd])
        C_on = np.ascontiguousarray(C[:, :nd, :nd])
        view = types.SimpleNamespace(n_channels=ND // nd)

        def run(xp, nb):
            with fb.BondBlockStore(nmat, nvol, ND, nd,
                                   ("chibar", "W", "chi_s_w", "chi_c_w")) as store, \
                    fb.BondDeviceContext.for_view(xp, S, C, S_on, C_on, view, 2) as dev:
                store.put_freq_batch("chibar", 0, nmat, chi_bar)
                res = fb.dress_and_build_w(store, dev, nb=nb, output_full=True, nmat=nmat,
                                           nvol=nvol, nd=nd, spatial_shape=(nvol, 1, 1),
                                           second_order="takimoto")
                return (store.get_freq_batch("W", 0, nmat).copy(),
                        store.get_freq_batch("chi_s_w", 0, nmat).copy(),
                        store.get_freq_batch("chi_c_w", 0, nmat).copy(), res)

        W_ref, chi_s_ref, chi_c_ref, r_ref = run(np, nmat)
        mat = np.eye(ND)[None, None] - chi_bar @ S[None]
        cm = _cond_max(mat)
        for nb in sorted({1, 3, nmat // 2, nmat}):
            with self.subTest(nb=nb):
                W, chi_s_w, chi_c_w, r = run(self.cupy, nb)
                _deviation("W (nb={})".format(nb), W, W_ref)
                _deviation("chi_s_w (nb={})".format(nb), chi_s_w, chi_s_ref)
                _deviation("chi_c_w (nb={})".format(nb), chi_c_w, chi_c_ref)
                np.testing.assert_allclose(W, W_ref, rtol=1e-10 * cm,
                                           atol=1e-12 * np.max(np.abs(W_ref)))
                np.testing.assert_allclose(chi_s_w, chi_s_ref, rtol=1e-10 * cm,
                                           atol=1e-12 * np.max(np.abs(chi_s_ref)))
                np.testing.assert_allclose(chi_c_w, chi_c_ref, rtol=1e-10 * cm,
                                           atol=1e-12 * np.max(np.abs(chi_c_ref)))
                np.testing.assert_allclose(r.collapse_s, r_ref.collapse_s,
                                           rtol=1e-10 * cm, atol=1e-13)
                np.testing.assert_allclose(r.collapse_c, r_ref.collapse_c,
                                           rtol=1e-10 * cm, atol=1e-13)
                np.testing.assert_allclose(r.static_s, r_ref.static_s,
                                           rtol=1e-10 * cm, atol=1e-13)
                np.testing.assert_allclose(r.static_c, r_ref.static_c,
                                           rtol=1e-10 * cm, atol=1e-13)
                self.assertAlmostEqual(r.cond_min_s, r_ref.cond_min_s, delta=1e-10)
                self.assertAlmostEqual(r.cond_min_c, r_ref.cond_min_c, delta=1e-10)

    def test_interval_guard_whole_map(self):
        from hwave.solver import flex_bond as fb
        chi_bar, S, C = _problem(nmat=8, nvol=4, nd=4, B=3, seed=13)
        # Push ONE block toward the instability so its score is far below
        # every other block's (factor tuned by a host dry run, recorded in
        # the task report): at factor 100 the spin-channel score of block
        # (l=0, q=0) is ~0.0153 against a minimum of ~0.4532 over the other
        # 31 blocks (ratio ~0.034, well under the 1/20 the exact set is
        # pruned against), and the whole map still passes the guard on both
        # channels (min score over the map stays above cond_tol = 1e-3) --
        # the uniform seed=13 fixture alone gives every block nearly the
        # same score, which makes the pruning count assertions below
        # vacuous (interval decomposes the same 64 of 64 blocks as svd).
        chi_bar[0, 0] *= 100.0
        nmat, nvol, ND, nd = chi_bar.shape[0], chi_bar.shape[1], S.shape[-1], 4
        S_on = np.ascontiguousarray(S[:, :nd, :nd]); C_on = np.ascontiguousarray(C[:, :nd, :nd])
        view = types.SimpleNamespace(n_channels=ND // nd)
        def run(xp, method):
            with fb.BondBlockStore(nmat, nvol, ND, nd, ("chibar", "W")) as store, \
                    fb.BondDeviceContext.for_view(xp, S, C, S_on, C_on, view, 2) as dev:
                store.put_freq_batch("chibar", 0, nmat, chi_bar)
                res = fb.dress_and_build_w(store, dev, nb=3, output_full=False, nmat=nmat,
                                           nvol=nvol, nd=nd, spatial_shape=(nvol, 1, 1),
                                           second_order="takimoto", guard_method=method)
                return store.get_freq_batch("W", 0, nmat).copy(), res
        W_ref, r_ref = run(self.cupy, "svd")
        W, r = run(self.cupy, "interval")
        _W_host, r_host = run(np, "interval")
        np.testing.assert_array_equal(W, W_ref)               # same solves, same guard
        self.assertEqual((r.cond_min_s, r.cond_min_c), (r_ref.cond_min_s, r_ref.cond_min_c))
        # the svd path decomposes every guarded block, on both channels
        self.assertEqual(r_ref.guard_exact_blocks, r_ref.guard_blocks)
        # device and host intervals select the same exact set
        self.assertEqual(r.guard_exact_blocks, r_host.guard_exact_blocks)
        # the spread above gives the interval guard something to prune
        self.assertLess(r.guard_exact_blocks, r.guard_blocks)


class TestTransportEquivalence(_GpuCase):
    """``flex_bond.calc_self_energy_bond`` (spec 4.3) in both array modules."""

    def test_sigma_matches(self):
        from tests.test_flex_bond_sigma_transport import _transport_problem
        from hwave.solver import flex_bond as fb
        store, green_kw, beta, view, shape, norb = _transport_problem()
        with store:
            ref = fb.calc_self_energy_bond(store, green_kw, beta, view, shape, norb, 1)
            out_d = fb.calc_self_energy_bond(store, green_kw, beta, view, shape, norb, 1,
                                             xp=self.cupy)
        self.assertEqual(type(out_d).__module__.split(".")[0], "cupy")
        out = self.cupy.asnumpy(out_d)
        _deviation("sigma (bond transport)", out, ref)
        np.testing.assert_allclose(out, ref, rtol=1e-11, atol=1e-13 * np.max(np.abs(ref)))

    def test_sigma_matches_multichannel(self):
        """The B >= 3 counterpart (issue #198 item 2): the B = 1 fixture
        above never takes the rolled-block branch (``xp.roll`` for
        ``alpha != beta``) of ``calc_self_energy_bond``, since a single
        channel has no off-diagonal (alpha, beta) pair. This fixture does."""
        from tests.test_flex_bond_sigma_transport import _transport_problem_multichannel
        from hwave.solver import flex_bond as fb
        store, green_kw, beta, view, shape, norb = _transport_problem_multichannel()
        with store:
            ref = fb.calc_self_energy_bond(store, green_kw, beta, view, shape, norb, 1)
            out_d = fb.calc_self_energy_bond(store, green_kw, beta, view, shape, norb, 1,
                                             xp=self.cupy)
        self.assertEqual(type(out_d).__module__.split(".")[0], "cupy")
        out = self.cupy.asnumpy(out_d)
        _deviation("sigma (bond transport, multichannel)", out, ref)
        np.testing.assert_allclose(out, ref, rtol=1e-11, atol=1e-13 * np.max(np.abs(ref)))


class TestMuSearchDevicePath(_GpuCase):
    """The Phase B chemical-potential search under cupy.

    ``_find_mu_dressed`` replaces the band energies with the host-built
    Heff reference spectrum (``ew_ref``) whenever the bond path supplies
    one, while the eigenvalues it is combined with live on the solver's
    backend. Before the fix that went with this module the host array was
    used as-is and the very first trial mu of the very first GPU iteration
    raised ``TypeError: Unsupported type <class 'numpy.ndarray'>``, so no
    bond-gate solve could run on a device at all. One iteration is enough
    to reach it."""

    def test_phase_b_mu_search_matches_the_host(self):
        from tests.test_flex_bond_gate import _flex
        mus = {}
        for gpu in (False, True):
            s, r = _flex({"gpu": gpu, "IterationMax": 1})
            gi = r.get_param("green")
            with tempfile.TemporaryDirectory() as out:
                s.solve(gi, out)
            mus[gpu] = float(s.mu)
            if gpu:
                self.assertEqual(str(s._bond_xp_name), "cupy")
        self.assertTrue(np.isfinite(mus[True]))
        if _REPORT:
            print("[gpu-equivalence] mu after one Phase B iteration: cpu "
                  "{:.15e}, gpu {:.15e}".format(mus[False], mus[True]))
        self.assertAlmostEqual(mus[True], mus[False], delta=1e-10)


class TestEndToEndEquivalence(_GpuCase):
    """One whole bond-gate FLEX solve on each backend."""

    @heavy
    def test_two_orbital_bond_chain(self):
        """tests/rpa/input_2orb with the bond gate, 4x4, Nmat 8, T = 2 (far
        from any instability): the converged flag agrees, the iteration count
        to within 2, the chemical potential to 1e-9, and the self-energy and
        the bond spin susceptibility to rtol 1e-8 / atol 1e-10.

        The loop is driven to EPS = 1e-12 (45 linear iterations at Mix = 0.5,
        ~0.6 s per solve on CPU), NOT to the 1e-8 the comparison uses. At
        EPS = 1e-8 the run stops while the fixed point is still moving, and
        two equally valid ways of reaching it (measured here by varying Mix
        alone, on one backend) already separate by 2e-8 in Sigma -- above the
        comparison's own rtol. Comparing there would measure the stopping
        rule, not the two backends; converging two orders further makes the
        residual (2e-13) negligible against the tolerance, so what is left is
        the backends' round-off."""
        from tests.test_flex_bond_gate import _flex
        outs = {}
        for gpu in (False, True):
            s, r = _flex({"gpu": gpu, "IterationMax": 200, "EPS": 12, "Mix": 0.5})
            gi = r.get_param("green")
            with tempfile.TemporaryDirectory() as out:
                s.solve(gi, out)
                outs[gpu] = (bool(s.scf_converged), int(s.scf_iterations), float(s.mu),
                             np.array(gi["sigma"]), np.array(gi["longitudinal_bond_chi_s"]))
            if gpu:
                # the run really took the device path, so a silent numpy
                # fallback cannot pass this test by comparing CPU with CPU
                self.assertEqual(str(s._bond_xp_name), "cupy")
        (ca, ia, ma, sa, xa), (cb, ib, mb, sb, xb) = outs[False], outs[True]
        _deviation("end-to-end sigma", sb, sa)
        _deviation("end-to-end bond chi_s", xb, xa)
        if _REPORT:
            print("[gpu-equivalence] end-to-end: converged {}/{}, iterations {}/{}, "
                  "mu {:.15e}/{:.15e}".format(ca, cb, ia, ib, ma, mb))
        # a run that stopped at the iteration cap would compare two moving
        # states, where the tolerances below mean nothing
        self.assertTrue(ca, "the CPU reference run did not converge")
        self.assertEqual(ca, cb)
        self.assertLessEqual(abs(ia - ib), 2)
        self.assertAlmostEqual(ma, mb, delta=1e-9)
        np.testing.assert_allclose(sb, sa, rtol=1e-8, atol=1e-10)
        np.testing.assert_allclose(xb, xa, rtol=1e-8, atol=1e-10)


if __name__ == "__main__":
    unittest.main()
