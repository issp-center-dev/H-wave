"""The device-resident bond block store (issue #275), on the CPU.

A fake ``cupy`` (:mod:`tests.fake_device`) makes the backend see a GPU
whose arrays are numpy arrays of a subclass, so the device-store and the
host-store runs of the SAME kernels execute numpy's arithmetic and can be
compared BIT FOR BIT at a pinned frequency batch: what may change with the
residency is only where the store lives and which host <-> device copies
happen, and both are pinned here -- the values by ``assert_array_equal``,
the copies by the fake's transfer record.

Covered: the store's own contract on either module; the bubble, the
dressing and the transport against a host store; whole bond-gate FLEX
solves (pairing none and singlet, full and half grid); the residency
admission (``longitudinal_bond_store``) and its log; the fallback to the
host store when allocating the device store runs out of memory; the key's
parsing. The real-device counterparts are in ``tests/test_flex_bond_gpu.py``
(``TestDeviceStoreEquivalence``).
"""
import contextlib
import math
import tempfile
import types
import unittest
from unittest import mock

import numpy as np

from tests.fake_device import DevArray, fake_cupy, is_device
from tests.heavy_tests import heavy


def _host(x):
    return np.array(x.view(np.ndarray) if isinstance(x, DevArray) else x, copy=True)


class TestStoreOnTheDevice(unittest.TestCase):
    """``BondBlockStore(xp=...)``: placement, conversion, views."""

    def test_device_slots_views_and_conversion(self):
        from hwave.solver.flex_bond import BondBlockStore
        nmat, nvol, nd, B = 4, 3, 4, 2
        rng = np.random.default_rng(0)
        blk = rng.normal(size=(nmat, nvol, nd, nd)) + 1j * rng.normal(size=(nmat, nvol, nd, nd))
        with fake_cupy() as cp:
            with BondBlockStore(nmat, nvol, B * nd, nd, ("chibar", "W"), xp=cp) as st:
                self.assertEqual(st.resident, "device")
                self.assertIs(st.xp, cp)
                self.assertEqual(st.nbytes, 2 * nmat * nvol * (B * nd) ** 2 * 16)
                cp.reset()
                st.put_pair("chibar", 1, 0, blk)                     # host block: one H2D
                self.assertEqual(cp.transfers, [("h2d", blk.shape, blk.nbytes)])
                view = st.get_pair("chibar", 1, 0)
                self.assertTrue(is_device(view))
                np.testing.assert_array_equal(_host(view), blk)
                cp.reset()
                st.put_pair("W", 0, 1, cp.asarray(blk))              # device block
                cp.transfers.pop()                                    # the asarray above
                self.assertEqual(cp.transfers, [])
                batch = st.get_freq_batch("chibar", 1, 3)
                self.assertTrue(is_device(batch))
                self.assertTrue(batch.flags.c_contiguous)
                # a view: writing through it reaches the slot
                batch[0, 0, 0, 0] = 7.0
                self.assertEqual(complex(st.get_freq_batch("chibar", 1, 2)[0, 0, 0, 0]), 7.0)
                cp.reset()
                st.put_freq_batch("W", 0, 2, np.ones((2, nvol, B * nd, B * nd), complex))
                self.assertEqual([d for d, _s, _n in cp.transfers], ["h2d"])
                w = st.detach("W")
                self.assertTrue(is_device(w))

    def test_host_store_converts_a_device_block_once(self):
        """A host store on the GPU path receives device blocks and copies
        each to the host exactly once, in the store (issue #275: the call
        sites no longer transfer)."""
        from hwave.solver.flex_bond import BondBlockStore
        nmat, nvol, nd, B = 4, 3, 4, 2
        blk = np.arange(nmat * nvol * nd * nd, dtype=float).reshape(nmat, nvol, nd, nd) + 0j
        with fake_cupy() as cp:
            d = cp.asarray(blk)
            with BondBlockStore(nmat, nvol, B * nd, nd, ("chibar",)) as st:
                self.assertEqual(st.resident, "host")
                cp.reset()
                st.put_pair("chibar", 0, 1, d)
                self.assertEqual(cp.transfers, [("d2h", blk.shape, blk.nbytes)])
                cp.reset()
                st.put_freq_batch("chibar", 0, 2, cp.asarray(np.zeros((2, nvol, B * nd, B * nd),
                                                                       complex)))
                self.assertEqual([x[0] for x in cp.transfers], ["h2d", "d2h"])
                self.assertFalse(is_device(st.get_freq_batch("chibar", 0, 1)))
                np.testing.assert_array_equal(st.get_pair("chibar", 0, 1)[2:], blk[2:])

    def test_get_pair_full_on_either_module(self):
        from hwave.solver import half_grid as hg
        from hwave.solver.flex_bond import BondBlockStore
        from tests.test_bond_half_grid import _hermitian_in_frequency
        nmat, nvol, nd, B = 8, 3, 4, 2
        x = _hermitian_in_frequency(nmat, nvol, B * nd, seed=3)
        for half in (False, True):
            with self.subTest(half=half), fake_cupy() as cp:
                src = hg.half_from_full(x, nmat) if half else x
                with BondBlockStore(nmat, nvol, B * nd, nd, ("W",), half_grid=half) as h, \
                        BondBlockStore(nmat, nvol, B * nd, nd, ("W",), half_grid=half,
                                       xp=cp) as d:
                    h.put_freq_batch("W", 0, h.nslices, src)
                    d.put_freq_batch("W", 0, d.nslices, src)
                    for a in range(B):
                        for b in range(B):
                            ref = h.get_pair_full("W", a, b)
                            cp.reset()
                            out = d.get_pair_full("W", a, b, xp=cp)        # no copy across
                            self.assertEqual(cp.transfers, [])
                            self.assertTrue(is_device(out) and out.flags.c_contiguous)
                            np.testing.assert_array_equal(_host(out), ref)
                            out_h = d.get_pair_full("W", a, b, xp=np)       # one D2H
                            self.assertEqual([t[0] for t in cp.transfers], ["d2h"])
                            self.assertFalse(is_device(out_h))
                            np.testing.assert_array_equal(out_h, ref)
                            # a host store asked for the device: today's single H2D
                            cp.reset()
                            out_d = h.get_pair_full("W", a, b, xp=cp)
                            self.assertEqual([t[0] for t in cp.transfers], ["h2d"])
                            np.testing.assert_array_equal(_host(out_d), ref)
                    # expand_full into a host buffer from a device store
                    buf = np.empty((nmat, nvol, B * nd, B * nd), complex)
                    d.expand_full("W", out=buf)
                    np.testing.assert_array_equal(buf, h.expand_full("W", out=np.empty_like(buf)))

    def test_aliases_release_and_rollback_on_device(self):
        from hwave.solver.flex_bond import BondBlockStore
        with fake_cupy() as cp:
            with BondBlockStore(4, 2, 2, 1, ("chibar",), {"W": "chibar"}, xp=cp) as st:
                st.put_freq_batch("W", 0, 4, np.ones((4, 2, 2, 2), complex))
                self.assertEqual(complex(st.get_freq_batch("chibar", 0, 1)[0, 0, 0, 0]), 1.0)
            with BondBlockStore(4, 2, 2, 1, ("chibar", "W"), xp=cp) as st:
                st.release_slot("W")
                with self.assertRaises(KeyError):
                    st.get_pair("W", 0, 0)
            # an out-of-memory error on the SECOND slot: the first is rolled
            # back and the device error propagates
            real = cp.zeros
            first = []

            def zeros(*a, **k):
                if first:
                    raise cp.cuda.memory.OutOfMemoryError("boom")
                first.append(real(*a, **k))
                return first[0]
            with mock.patch.object(cp, "zeros", zeros):
                with self.assertRaises(cp.cuda.memory.OutOfMemoryError):
                    BondBlockStore(4, 2, 2, 1, ("chibar", "W"), xp=cp)
            self.assertEqual(len(first), 1)


# ---------------------------------------------------------------------------
# the kernels: a device store against a host store, on the same device
# ---------------------------------------------------------------------------

class TestBubbleIntoTheStore(unittest.TestCase):

    def _run(self, cp, fx, green, tail, store_xp, half):
        from hwave.solver import backend as bk, flex_bond as fb
        nmat, nvol, nd, ND = fx["nmat"], fx["nvol"], fx["nd"], fx["ND"]
        with fb.BondBlockStore(nmat, nvol, ND, nd, ("chibar",), half_grid=half,
                               xp=store_xp) as st, \
                fb.BondDeviceContext(cp, fx["S"], fx["C"], green0_tail=tail) as ctx:
            g = bk.to_device(green, cp)
            cp.reset()
            fb.assemble_bubble(st, g, ctx.green0_tail, fx["beta"], fx["view"],
                               fx["spatial_shape"], 1)
            moved = list(cp.transfers)
            blocks = {(a, b): _host(st.get_pair("chibar", a, b))
                      for a in range(fx["B"]) for b in range(fx["B"])}
        return blocks, moved

    def test_device_store_equals_host_store_bitwise(self):
        from tests.test_flex_bond_gpu import _bubble_fixture
        fx, green, tail = _bubble_fixture()
        for half in (False, True):
            with self.subTest(half=half), fake_cupy() as cp:
                ref, moved_h = self._run(cp, fx, green, tail, np, half)
                out, moved_d = self._run(cp, fx, green, tail, cp, half)
                self.assertEqual(sorted(out), sorted(ref))
                for k in ref:
                    np.testing.assert_array_equal(out[k], ref[k])
                # the host store brings every block back; the device store none
                self.assertEqual(sum(1 for d, s, _ in moved_h if d == "d2h" and len(s) == 4),
                                 fx["B"] ** 2)
                self.assertEqual([t for t in moved_d if len(t[1]) == 4], [])


def _dress(cp, store_xp, half, *, aliased, nb=3, output_full=False, bond_w2=False):
    """One dress_and_build_w on a store of module ``store_xp`` with the
    vertices on the fake device; returns (W, DressResult, transfers, extra)."""
    from hwave.solver import flex_bond as fb, half_grid as hg
    from tests.test_bond_half_grid import _hermitian_in_frequency
    from tests.test_flex_bond_dressing import _problem
    from tests.test_flex_bond_gpu import _onsite_factors_on
    nmat, nvol, nd, B = 8, 4, 4, 3
    ND = B * nd
    _x, S, C = _problem(nmat=nmat, nvol=nvol, nd=nd, B=B, seed=5)
    x = 0.05 * _hermitian_in_frequency(nmat, nvol, ND, seed=5)
    view = types.SimpleNamespace(n_channels=B)
    names = ("chibar",) if aliased else ("chibar", "W")
    if output_full:
        names += ("chi_s_w", "chi_c_w")
    aliases = {"W": "chibar"} if aliased else None
    kw = {}
    if bond_w2:
        rng = np.random.default_rng(251)
        Af = rng.normal(size=(nvol, ND, ND)) + 1j * rng.normal(size=(nvol, ND, ND))
        Bf = rng.normal(size=(nvol, ND, ND)) + 1j * rng.normal(size=(nvol, ND, ND))
        kw["bond_w2"] = [((0, 0, 0), cp.asarray(0.1 * Af), cp.asarray(0.1 * Bf))]
    else:
        kw["factors"] = _onsite_factors_on(cp)
    with fb.BondBlockStore(nmat, nvol, ND, nd, names, aliases, half_grid=half,
                           xp=store_xp) as st, \
            fb.BondDeviceContext.for_view(cp, S, C, view, 2) as dev:
        st.put_freq_batch("chibar", 0, st.nslices, hg.half_from_full(x, nmat) if half else x)
        cp.reset()
        res = fb.dress_and_build_w(st, dev, nb=nb, output_full=output_full, nmat=nmat,
                                   nvol=nvol, nd=nd, spatial_shape=(nvol, 1, 1), **kw)
        moved = list(cp.transfers)
        W = _host(st.get_freq_batch("W", 0, st.nslices))
        extra = {}
        if output_full:
            extra = {n: _host(st.get_freq_batch(n, 0, st.nslices)) for n in ("chi_s_w", "chi_c_w")}
    return W, res, moved, extra, ND


class TestDressingWithADeviceStore(unittest.TestCase):

    def _compare(self, **kw):
        for half in (False, True):
            for aliased in (True, False):
                with self.subTest(half=half, aliased=aliased, **kw), fake_cupy() as cp:
                    W_h, r_h, moved_h, ex_h, ND = _dress(cp, np, half, aliased=aliased, **kw)
                    W_d, r_d, moved_d, ex_d, _ = _dress(cp, cp, half, aliased=aliased, **kw)
                    np.testing.assert_array_equal(W_d, W_h)
                    for name in ("collapse0", "collapse_s", "collapse_c", "static_s",
                                 "static_c"):
                        np.testing.assert_array_equal(getattr(r_d, name), getattr(r_h, name))
                    self.assertEqual((r_d.cond_min_s, r_d.cond_min_c, r_d.cond_bracket_s,
                                      r_d.guard_exact_blocks, r_d.mirror_guarded),
                                     (r_h.cond_min_s, r_h.cond_min_c, r_h.cond_bracket_s,
                                      r_h.guard_exact_blocks, r_h.mirror_guarded))
                    for n in ex_h:
                        np.testing.assert_array_equal(ex_d[n], ex_h[n])
                    # the host store moves every chibar batch in and every W
                    # batch out; the device store moves none of them
                    big_h = [t for t in moved_h if len(t[1]) == 4 and t[1][-2:] == (ND, ND)]
                    big_d = [t for t in moved_d if len(t[1]) == 4 and t[1][-2:] == (ND, ND)]
                    self.assertGreater(len(big_h), 0)
                    self.assertEqual(big_d, [])

    def test_local_second_order(self):
        self._compare()

    def test_bond_w2_and_output_full(self):
        self._compare(bond_w2=True, output_full=True)


class TestTransportFromADeviceStore(unittest.TestCase):

    def test_sigma_bitwise_and_no_block_transfers(self):
        from hwave.solver import backend as bk, flex_bond as fb, half_grid as hg
        from tests.test_bond_half_grid import _hermitian_in_frequency
        from tests.test_flex_bond_sigma_transport import _transport_problem_multichannel
        st0, G, beta, view, shape, norb = _transport_problem_multichannel()
        with st0:
            nmat, nvol, ND, nd = st0.nmat, st0.nvol, st0.ND, st0.nd
        W = 0.1 * _hermitian_in_frequency(nmat, nvol, ND, seed=9)
        for half in (False, True):
            with self.subTest(half=half), fake_cupy() as cp:
                g = bk.to_device(G, cp)
                outs, moved = {}, {}
                for name, sx in (("host", np), ("device", cp)):
                    with fb.BondBlockStore(nmat, nvol, ND, nd, ("W",), half_grid=half,
                                           xp=sx) as st:
                        st.put_freq_batch("W", 0, st.nslices,
                                          hg.half_from_full(W, nmat) if half else W)
                        cp.reset()
                        sig = fb.calc_self_energy_bond(st, g, beta, view, shape, norb, 1, xp=cp)
                        moved[name] = list(cp.transfers)
                        self.assertTrue(is_device(sig))
                        outs[name] = _host(sig)
                np.testing.assert_array_equal(outs["device"], outs["host"])
                B = view.n_channels
                self.assertEqual(sum(1 for d, s, _ in moved["host"] if d == "h2d" and len(s) == 4),
                                 B * B)
                self.assertEqual([t for t in moved["device"] if len(t[1]) == 4], [])


# ---------------------------------------------------------------------------
# whole solves
# ---------------------------------------------------------------------------

def _solve(store, pairing="none", half=False, extra=None, cp_kw=None, setup=None):
    """One bond-gate FLEX solve of tests/test_bond_half_grid's "UV" fixture on
    the fake device (frequency batch pinned to 3 there); returns the solver,
    green_info, the archives and the fake's transfer record."""
    from tests.test_bond_half_grid import _ELI_UNIFORM, _run
    par = {"gpu": True, "longitudinal_bond_store": store,
           "longitudinal_bond_pairing": pairing}
    par.update(extra or {})
    with fake_cupy(**(cp_kw or {})) as cp:
        if setup is not None:
            setup(cp)
        s, gi, files = _run("UV", par, half=half,
                            eli=None if pairing == "none" else _ELI_UNIFORM)
        return s, gi, files, list(cp.transfers)


def _assert_same_outputs(test, a, b):
    """Every archive member and every array of green_info, bit for bit."""
    (_sa, gia, fa, _ta), (_sb, gib, fb, _tb) = a, b
    test.assertEqual(sorted(fa), sorted(fb))
    for fn in fa:
        test.assertEqual(sorted(fa[fn]), sorted(fb[fn]), fn)
        for m in fa[fn]:
            x, y = fa[fn][m], fb[fn][m]
            test.assertEqual(x.shape, y.shape, (fn, m))
            if x.dtype.kind in "fc":
                np.testing.assert_array_equal(x, y, err_msg="{}:{}".format(fn, m))
            else:
                test.assertTrue(np.array_equal(x, y), (fn, m))
    for k in set(gia) | set(gib):
        x, y = gia.get(k), gib.get(k)
        if isinstance(x, np.ndarray) or isinstance(y, np.ndarray):
            test.assertFalse(is_device(x) or is_device(y), k)
            np.testing.assert_array_equal(np.asarray(x), np.asarray(y), err_msg=str(k))


class TestEndToEnd(unittest.TestCase):
    """Whole solves: a device store gives the host store's outputs bit for
    bit at the pinned batch, without the store's copies."""

    def _check(self, pairing, half):
        host = _solve("host", pairing, half)
        dev = _solve("device", pairing, half)
        self.assertEqual(host[0]._bond_store_resident, "host")
        self.assertEqual(dev[0]._bond_store_resident, "device")
        self.assertEqual(host[0]._bond_nb, dev[0]._bond_nb)
        _assert_same_outputs(self, host, dev)
        if pairing != "none":
            self.assertTrue(any(f.startswith("eliashberg_bond_") for f in dev[2]))
            self.assertFalse(any("error" in str(k) for k in dev[1]), dev[1].keys())
        ND, nd = int(dev[0]._bond_est["ND"]), dev[0].norb ** 2
        # the store's copies: chibar / W batches and the pair blocks to the
        # device happen with the host store and never with the device store
        self.assertGreater(len(_big(host[3], ND, nd)), 0)
        self.assertEqual(_big(dev[3], ND, nd), [])

    def test_pairing_none_full_grid(self):
        self._check("none", False)

    def test_pairing_none_half_grid(self):
        self._check("none", True)

    def test_pairing_singlet_full_grid(self):
        self._check("singlet", False)

    def test_pairing_singlet_half_grid(self):
        self._check("singlet", True)

    def test_auto_takes_the_device_when_it_fits(self):
        s = _solve("auto")[0]
        self.assertEqual(s._bond_store_resident, "device")
        self.assertEqual(s._bond_est["store_resident"], "device")
        self.assertIn("block_store", s._bond_est["device_rows"])


def _big(transfers, ND, nd):
    """The recorded copies of the store's traffic: 4-D frequency batches
    ``(nb, nvol, ND, ND)`` either way, and 4-D channel-pair blocks ``(n,
    nvol, nd, nd)`` going to the device (the dressing's collapses come back
    in that shape, so the device-to-host direction of it is not the
    store's)."""
    return [(d, s) for d, s, _n in transfers
            if len(s) == 4 and (s[-2:] == (ND, ND) or (d == "h2d" and s[-2:] == (nd, nd)))]


class TestOutOfMemoryFallback(unittest.TestCase):
    """Allocating the device store fails: "auto" continues on the host store
    before any arithmetic; "device" refuses."""

    @staticmethod
    def _fail_store(cp):
        # every zeros request of at least one store slot fails once (the store
        # is the first allocation of that size)
        nslices, nvol, ND = 8, 16, 20
        cp.fail_zeros, cp.fail_zeros_min = 1, nslices * nvol * ND * ND * 16

    def test_auto_falls_back_to_the_host_store(self):
        ref = _solve("host")
        with self.assertLogs("hwave.solver.flex", level="WARNING") as cm:
            out = _solve("auto", setup=self._fail_store)
        s = out[0]
        self.assertEqual(s._bond_store_resident, "host")
        self.assertIs(s._bond_est, s._bond_est_host_store)
        self.assertEqual(s._bond_nb, int(s._bond_est_host_store["nb"]))
        self.assertTrue(any("allocating the bond block store on the device failed" in m
                            for m in cm.output), cm.output)
        _assert_same_outputs(self, ref, out)

    def test_device_refuses(self):
        from tests.fake_device import _OutOfMemoryError
        with self.assertLogs("hwave.solver.flex", level="ERROR") as cm:
            with self.assertRaises(_OutOfMemoryError):
                _solve("device", setup=self._fail_store)
        self.assertTrue(any("longitudinal_bond_store = \"device\"" in m for m in cm.output))

    def test_auto_fallback_restores_the_host_store_batch(self):
        """Without a pinned batch the device table can choose a smaller batch
        than the host-store table; the fallback goes back to the latter."""
        def setup(cp):
            self._fail_store(cp)
        s = _solve("auto", extra={"longitudinal_bond_freq_batch": None}, setup=setup)[0]
        self.assertEqual(s._bond_store_resident, "host")
        self.assertEqual(s._bond_nb, int(s._bond_est_host_store["nb"]))


# ---------------------------------------------------------------------------
# admission
# ---------------------------------------------------------------------------

class TestAdmission(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        from tests.test_flex_bond_gate import _flex
        cls.solvers = {}
        for key, extra in (("none", {}), ("singlet", {"longitudinal_bond_pairing": "singlet",
                                                      "IterationMax": 1}),
                           ("full", {"longitudinal_bond_output_full": True}),
                           ("half", {"longitudinal_bond_half_grid": True})):
            s, r = _flex(dict({"IterationMax": 0}, **extra), mixing="anderson")
            gi = r.get_param("green")
            with tempfile.TemporaryDirectory() as out:
                s.solve(gi, out)
            cls.solvers[key] = s

    def _est(self, key, available, store, **kw):
        from hwave.solver import flex_bond
        args = dict(self.solvers[key]._bond_est_kwargs, **kw)
        return flex_bond.estimate_bond_memory(device_available=available, store=store, **args)

    def test_store_row_and_slot_count(self):
        for key, slots in (("none", 1), ("singlet", 2), ("full", 3), ("half", 1)):
            with self.subTest(key=key):
                est = self._est(key, 2 ** 50, "device")
                self.assertEqual(est["store_resident"], "device")
                self.assertEqual(est["store_bytes"], slots * est["U_store"])
                self.assertEqual(est["device_rows"]["block_store"], est["store_bytes"])
                host = self._est(key, 2 ** 50, "host")
                self.assertNotIn("block_store", host["device_rows"])
                self.assertEqual(host["store_resident"], "host")
                # the host table does not depend on the residency
                self.assertEqual(est["persistent_rows"], host["persistent_rows"])
        # the half grid halves the store row
        self.assertLess(self._est("half", 2 ** 50, "device")["store_bytes"],
                        self._est("none", 2 ** 50, "device")["store_bytes"])

    def _need(self, key, nb):
        """1.25 (persistent device rows + store + max phase) at batch nb."""
        est = self._est(key, 2 ** 62, "device", freq_batch=nb)
        r = est["device_rows"]
        phases = ("bubble", "dressing", "transport", "guard")
        need = 1.25 * (sum(v for k, v in r.items() if k not in phases)
                       + max(r["bubble"], r["dressing"] + r["guard"], r["transport"]))
        self.assertAlmostEqual(need, est["device_need"], delta=1e-6 * need)
        return need

    def test_auto_follows_the_formula_at_the_batch_floor(self):
        from hwave.solver import flex_bond
        s = self.solvers["none"]
        nb_host = self._est("none", 2 ** 62, "host")["nb"]
        floor = min(nb_host, math.ceil(flex_bond.DEVICE_STORE_MIN_BLOCKS / s.lattice.nvol))
        need = self._need("none", floor)
        fit = self._est("none", int(need * 1.001 / 0.9), "auto")
        self.assertEqual(fit["store_resident"], "device")
        self.assertGreaterEqual(fit["nb"], floor)
        miss = self._est("none", int(need * 0.999 / 0.9), "auto")
        self.assertEqual(miss["store_resident"], "host")
        self.assertIn("does not fit", miss["store_reason"])
        self.assertNotIn("block_store", miss["device_rows"])

    def test_device_refuses_below_nb_1_and_admits_above(self):
        need1 = self._need("singlet", 1)
        est = self._est("singlet", int(need1 * 1.001 / 0.9), "device")
        self.assertEqual((est["store_resident"], est["nb"]), ("device", 1))
        with self.assertRaisesRegex(ValueError, "longitudinal_bond_store = \"device\".*"
                                                "block_store"):
            self._est("singlet", int(need1 * 0.999 / 0.9), "device")

    def test_a_pinned_batch_decides_at_that_batch(self):
        need3 = self._need("none", 3)
        est = self._est("none", int(need3 * 1.001 / 0.9), "auto", freq_batch=3)
        self.assertEqual((est["store_resident"], est["nb"]), ("device", 3))
        est = self._est("none", int(need3 * 0.999 / 0.9), "auto", freq_batch=3)
        self.assertEqual((est["store_resident"], est["nb"]), ("host", 3))

    def test_solver_admission_frees_the_pool_first_and_logs(self):
        import hwave.solver.flex as flex_mod
        s = self.solvers["none"]
        order = []
        fake_xp = types.SimpleNamespace(__name__="cupy")
        with mock.patch.object(flex_mod._bk, "free_device_pool",
                               lambda: order.append("free")), \
                mock.patch.object(flex_mod._bk, "device_available_bytes",
                                  lambda: (order.append("measure"), 2 ** 50)[1]):
            with self.assertLogs("hwave.solver.flex", level="INFO") as cm:
                s._bond_device_admission(fake_xp)
        self.assertEqual(order, ["free", "measure"])
        self.assertEqual(s._bond_store_resident, "device")
        joined = "\n".join(cm.output)
        self.assertIn("kept in DEVICE memory", joined)
        self.assertIn("block store on the device", joined)
        s.longitudinal_bond_store = "host"
        try:
            with mock.patch.object(flex_mod._bk, "device_available_bytes", lambda: 2 ** 50):
                with self.assertLogs("hwave.solver.flex", level="INFO") as cm:
                    s._bond_device_admission(fake_xp)
        finally:
            s.longitudinal_bond_store = "auto"
        self.assertEqual(s._bond_store_resident, "host")
        self.assertIn("kept in HOST memory", "\n".join(cm.output))

    def test_device_needs_a_device(self):
        s = self.solvers["none"]
        s.longitudinal_bond_store = "device"
        try:
            with self.assertRaisesRegex(ValueError, "needs the GPU backend"):
                s._bond_device_admission(np)
            import hwave.solver.flex as flex_mod
            with mock.patch.object(flex_mod._bk, "device_available_bytes", lambda: None):
                with self.assertRaisesRegex(ValueError, "could not be queried"):
                    s._bond_device_admission(types.SimpleNamespace(__name__="cupy"))
        finally:
            s.longitudinal_bond_store = "auto"

    def test_transfer_volume_with_a_device_store(self):
        s = self.solvers["none"]
        est = self._est("none", 2 ** 50, "device")
        with self.assertLogs("hwave.solver.flex", level="INFO") as cm:
            s._log_bond_transfer_volume(est)
        line = [m for m in cm.output if "transfer volume" in m][0]
        u = est["nslices"] * est["nvol"] * est["ND"] ** 2 * 16
        d2h = 2 * u + 2 * est["S_bytes"] + 3 * est["C_bytes"] + est["G_bytes"]
        self.assertIn("host to device 0.0000 GiB", line)
        self.assertIn("device to host {:.4f} GiB".format(d2h / 1024.0 ** 3), line)

    def test_cpu_run_is_unaffected(self):
        """On the numpy backend the store is the host store under every
        value but "device", which is refused."""
        from tests.test_flex_bond_gate import _flex
        for v in ("auto", "host"):
            s, r = _flex({"IterationMax": 1, "longitudinal_bond_store": v})
            with tempfile.TemporaryDirectory() as out:
                s.solve(r.get_param("green"), out)
            self.assertEqual(s._bond_store_resident, "host")
        s, r = _flex({"IterationMax": 1, "longitudinal_bond_store": "device"})
        with tempfile.TemporaryDirectory() as out:
            with self.assertRaisesRegex(ValueError, "needs the GPU backend"):
                s.solve(r.get_param("green"), out)


class TestKey(unittest.TestCase):

    def _raw(self, param, gate=True):
        from hwave.solver.flex import FLEX
        p = {"flex_hartree_fock": True, "longitudinal_bond_channels": gate, "Nmat": 8}
        p.update(param)
        return FLEX._parse_phase_b_keys({"param": p, "calc_scheme": "general"})

    def test_default_values_and_case(self):
        self.assertEqual(self._raw({})["longitudinal_bond_store"], "auto")
        for key in ("longitudinal_bond_store", "Longitudinal_Bond_Store"):
            for v, want in (("auto", "auto"), ("HOST", "host"), (" Device ", "device")):
                with self.subTest(key=key, v=v):
                    self.assertEqual(self._raw({key: v})["longitudinal_bond_store"], want)

    def test_strict_values(self):
        for v in ("gpu", "", True, 1, None, "pinned"):
            with self.subTest(v=v):
                with self.assertRaisesRegex(ValueError, "longitudinal_bond_store must be"):
                    self._raw({"longitudinal_bond_store": v})

    def test_stale_without_the_gate(self):
        from hwave.solver.flex import FLEX
        self.assertIn("longitudinal_bond_store", FLEX._PHASE_B_BOND_KEYS)
        with self.assertLogs("hwave.solver.flex", level="WARNING") as cm:
            out = self._raw({"longitudinal_bond_store": "garbage"}, gate=False)
        self.assertEqual(out["longitudinal_bond_store"], "auto")
        self.assertIn("longitudinal_bond_store", "\n".join(cm.output))
        with self.assertLogs("hwave.solver.flex", level="WARNING") as cm:
            out = FLEX._parse_phase_b_keys({"param": {"longitudinal_bond_store": "device"}})
        self.assertEqual(out["longitudinal_bond_store"], "auto")
        self.assertIn("longitudinal_bond_store", "\n".join(cm.output))

    def test_rpa_warns_it_is_flex_only(self):
        from tests.test_flex_bond_gate import _rpa
        with self.assertLogs("hwave.solver.rpa", level="WARNING") as cm:
            _rpa({"longitudinal_bond_store": "device"})
        self.assertTrue(any("longitudinal_bond_store" in m and "FLEX-only" in m
                            for m in cm.output), cm.output)


if __name__ == "__main__":
    unittest.main()
