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

#: the index of the fake's transfer record at which the in-process pairing
#: of the last :func:`_solve` started (None: no pairing step ran)
_PAIRING_MARK = []


def _solve(store, pairing="none", half=False, extra=None, cp_kw=None, setup=None):
    """One bond-gate FLEX solve of tests/test_bond_half_grid's "UV" fixture on
    the fake device (frequency batch pinned to 3 there); returns the solver,
    green_info, the archives and the fake's transfer record. The record
    position where the in-process pairing started is left in
    ``_PAIRING_MARK[0]``."""
    from hwave.solver import eliashberg_bond_io as io
    from tests.test_bond_half_grid import _ELI_UNIFORM, _run
    par = {"gpu": True, "longitudinal_bond_store": store,
           "longitudinal_bond_pairing": pairing}
    par.update(extra or {})
    _PAIRING_MARK[:] = [None]
    real = io.run_inprocess_pairing
    with fake_cupy(**(cp_kw or {})) as cp:
        def marked(*a, **k):
            _PAIRING_MARK[0] = len(cp.transfers)
            return real(*a, **k)
        if setup is not None:
            setup(cp)
        with mock.patch.object(io, "run_inprocess_pairing", marked):
            s, gi, files = _run("UV", par, half=half,
                                eli=None if pairing == "none" else _ELI_UNIFORM)
        return s, gi, files, list(cp.transfers)


#: the eigensolver's outputs: ARPACK follows its iterates, and on a real
#: BLAS a differently aligned operand can move them at round-off, so they are
#: compared to 1e-10 (the gap modulo its global phase), as in
#: tests/test_flex_bond_gpu.py; everything else bit for bit
_EIG_TOL = 1.0e-10


def _is_eig(name):
    return name.startswith("eliashberg_bond_") or name.startswith("pairing_") or \
        name.startswith("eigenvalue_bond_") or name.startswith("gap_bond_")


def _close_modulo_phase(test, x, y, label, phase=1.0):
    x, y = np.asarray(x), np.asarray(y)
    test.assertEqual(x.shape, y.shape, label)
    if x.dtype.kind not in "fc":
        test.assertTrue(np.array_equal(x, y), label)
        return
    if y.dtype.kind == "c":
        y = y * phase
    scale = float(np.max(np.abs(x))) if x.size else 0.0
    test.assertLessEqual(float(np.max(np.abs(x - y))) if x.size else 0.0,
                         _EIG_TOL * max(scale, 1e-300), label)


def _phase(x, y):
    ov = np.vdot(np.asarray(y).ravel(), np.asarray(x).ravel())
    return ov / abs(ov) if abs(ov) > 0 else 1.0


def _assert_same_outputs(test, a, b):
    """Every archive member and every array of green_info, bit for bit."""
    (_sa, gia, fa, _ta), (_sb, gib, fb, _tb) = a, b
    test.assertEqual(sorted(fa), sorted(fb))
    for fn in fa:
        test.assertEqual(sorted(fa[fn]), sorted(fb[fn]), fn)
        if _is_eig(fn):
            ph = _phase(fa[fn]["gap"], fb[fn]["gap"]) if "gap" in fa[fn] else 1.0
            for m in fa[fn]:
                _close_modulo_phase(test, fa[fn][m], fb[fn][m], (fn, m),
                                    ph if m.startswith("gap") and m != "gap_sector_weights"
                                    else 1.0)
            continue
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
            if _is_eig(str(k)):
                ph = (_phase(gia[str(k)], gib[str(k)]) if str(k).endswith("_gap") else 1.0)
                _close_modulo_phase(test, x, y, k, ph)
                continue
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
        mark = _PAIRING_MARK[0]
        scf, pair = (dev[3], []) if mark is None else (dev[3][:mark], dev[3][mark:])
        # the store's copies during the SCF: chibar / W batches and the pair
        # blocks to the device happen with the host store and never with the
        # device store
        self.assertGreater(len(_big(host[3], ND, nd)), 0)
        self.assertEqual(_big(scf, ND, nd), [])
        if pairing != "none":
            # before the kernel the store leaves the device once: the vertex
            # slot W comes to the host (chibar is released), as one slot
            st = dev[0]
            slot = (st._bond_est["nslices"], st.lattice.nvol, ND, ND)
            self.assertEqual([t for t in _big(pair, ND, nd) if t[1][-2:] == (ND, ND)],
                             [("d2h", slot)])

    def test_pairing_none_full_grid(self):
        self._check("none", False)

    def test_pairing_none_half_grid(self):
        self._check("none", True)

    def test_pairing_singlet_full_grid(self):
        self._check("singlet", False)

    def test_pairing_singlet_half_grid(self):
        self._check("singlet", True)

    def test_recorded_transfers_are_bounded_by_the_estimate(self):
        """The logged per-iteration transfer volume is an upper bound of what
        one bond map really copies (recorded by the fake), for either
        residency and direction."""
        import hwave.solver.flex as flex_mod
        for store in ("host", "device"):
            for half in (False, True):
                with self.subTest(store=store, half=half):
                    per_map = []
                    real = flex_mod.FLEX._phase_b_bond_map

                    def wrapped(self_, *a, **k):
                        from hwave.solver import backend as bk
                        cp = bk._import_cupy()
                        n0 = len(cp.transfers)
                        out = real(self_, *a, **k)
                        moved = cp.transfers[n0:]
                        # the per-map index vectors (the ring's pair
                        # permutation) are not part of the estimate
                        per_map.append((sum(b for d, sh, b in moved
                                            if d == "h2d" and len(sh) > 1),
                                        sum(b for d, _s, b in moved if d == "d2h")))
                        return out
                    with mock.patch.object(flex_mod.FLEX, "_phase_b_bond_map", wrapped):
                        s = _solve(store, half=half)[0]
                    h2d_est, d2h_est = s._bond_transfer_estimate
                    self.assertEqual(len(per_map), 3)
                    for h2d, d2h in per_map:
                        self.assertLessEqual(d2h, d2h_est)
                        self.assertGreater(d2h, 0)
                        if store == "host":
                            self.assertLessEqual(h2d, h2d_est)
                        else:
                            self.assertEqual(h2d_est, 0)
                            self.assertEqual(h2d, 0)

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


class TestStoreLeavesTheDeviceBeforeThePairingKernel(unittest.TestCase):
    """After the SCF a device store would hold device memory the pairing
    kernel needs (kugui, Nmat 1024: the kernel fell back to its host
    residency and the whole run got slower). The store leaves the device
    before the kernel's admission, so the kernel gets the device residency
    and the results are those of the host store."""

    def _device(self, release=True):
        """A solve whose device has room for the kernel only without the
        store: the device reading is large until the store exists, then a
        fixed capacity minus the store's live device slots."""
        from hwave.solver import backend as bk, eliashberg_bond as eb
        from hwave.solver import eliashberg_bond_io as io, flex_bond as fb
        table = eb.estimate_pair_memory(
            nmat=8, ntau=0, nfreq=8, nvol=16, norb=2, B=5, num_eigenvalues=3,
            residency="auto", ir=False, L_B=0, n_channels=1, in_process=True)
        d_need = table.need("device")[1]
        slots = 2 * 8 * 16 * 20 * 20 * 16          # chibar + W, full grid
        capacity = d_need / 0.9 + 0.5 * slots
        stores = []
        real_store = fb.BondBlockStore

        class Spy(real_store):
            def __init__(self, *a, **k):
                super().__init__(*a, **k)
                stores.append(self)

        def available():
            if not stores:
                return 2 ** 40
            live = sum(st.nbytes for st in stores
                       if not st.released and st.resident == "device")
            return int(capacity - live)
        patches = [mock.patch.object(fb, "BondBlockStore", Spy),
                   mock.patch.object(bk, "device_available_bytes", available)]
        if not release:
            patches.append(mock.patch.object(io, "_release_device_store",
                                             lambda *a, **k: None))
        return patches

    def _run(self, store, release=True):
        with contextlib.ExitStack() as st:
            for p in self._device(release):
                st.enter_context(p)
            cm = st.enter_context(self.assertLogs("qlms.solver.eliashberg_bond", level="INFO"))
            out = _solve(store, "singlet")
        res = [m.split("kernel residency ")[1] for m in cm.output if "kernel residency" in m]
        return out, res

    def test_kernel_on_the_device_and_results_unchanged(self):
        host, res_h = self._run("host")
        dev, res_d = self._run("device")
        self.assertEqual(dev[0]._bond_store_resident, "device")
        self.assertEqual(res_h, ["device"])
        self.assertEqual(res_d, ["device"])
        _assert_same_outputs(self, host, dev)
        # the control: with the store kept on the device the same reading
        # leaves the kernel no room, which is the situation the release fixes
        _kept, res_k = self._run("device", release=False)
        self.assertEqual(res_k, ["host"])


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
        """Without a pinned batch and with little device memory, the
        device-store table admits a SMALLER batch than the host-store table
        (the store takes memory the batch would use); after the fallback the
        host-store batch is the one in force. The batch floor is lowered so
        that "auto" admits the store at the small batch, and output_full
        makes the store (3 slots) larger than one batch row."""
        from hwave.solver import flex_bond
        extra = {"longitudinal_bond_freq_batch": None, "longitudinal_bond_output_full": True}
        probe = _solve("host", extra=extra)[0]
        kw = dict(probe._bond_est_kwargs)
        need = {n: flex_bond.estimate_bond_memory(device_available=2 ** 62, store="device",
                                                  **dict(kw, freq_batch=n))["device_need"]
                for n in range(1, 9)}
        chosen = None
        with mock.patch.object(flex_bond, "DEVICE_STORE_MIN_BLOCKS", 1):
            for n in range(1, 8):
                free = int(need[n] * 1.0001 / 0.9)
                dev = flex_bond.estimate_bond_memory(device_available=free, store="auto", **kw)
                host = flex_bond.estimate_bond_memory(device_available=free, store="host", **kw)
                if dev["store_resident"] == "device" and dev["nb"] < host["nb"]:
                    chosen = (free, dev["nb"], host["nb"])
                    break
        self.assertIsNotNone(chosen, "no device memory gives a smaller device-store batch")
        free, nb_dev, nb_host = chosen

        def setup(cp):
            self._fail_store(cp)
        with mock.patch.object(flex_bond, "DEVICE_STORE_MIN_BLOCKS", 1):
            s = _solve("auto", extra=extra, cp_kw={"free": free}, setup=setup)[0]
        self.assertEqual(s._bond_store_resident, "host")
        self.assertEqual(int(s._bond_est_host_store["nb"]), nb_host)
        self.assertEqual(s._bond_nb, nb_host)
        self.assertNotEqual(s._bond_nb, nb_dev)
        # and without the failure the device store runs at the smaller batch
        with mock.patch.object(flex_bond, "DEVICE_STORE_MIN_BLOCKS", 1):
            s = _solve("auto", extra=extra, cp_kw={"free": free})[0]
        self.assertEqual((s._bond_store_resident, s._bond_nb), ("device", nb_dev))


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

    #: an estimate-only problem at the issue's nvol (1024), so that
    #: nvol * (host batch) is far above DEVICE_STORE_MIN_BLOCKS
    _BIG = dict(nmat=64, nvol=1024, norb=2, B=9, depth=0, output_full=False,
                split_seed=False, n_types=1, freq_batch=None, cap_gb=1.0e6,
                mixing="linear", w_in_place=True)

    def _big(self, available, store, **kw):
        from hwave.solver import flex_bond
        return flex_bond.estimate_bond_memory(device_available=available, store=store,
                                              **dict(self._BIG, **kw))

    def test_batch_floor_is_decisive(self):
        """Available memory between the need at the floor ceil(4096 / nvol)
        and the need at the host batch: "auto" keeps the store on the device
        with a batch between the two; deciding at the host batch would not."""
        from hwave.solver import flex_bond
        nb_host = self._big(2 ** 62, "host")["nb"]
        floor = math.ceil(flex_bond.DEVICE_STORE_MIN_BLOCKS / self._BIG["nvol"])
        self.assertGreater(self._BIG["nvol"] * nb_host, flex_bond.DEVICE_STORE_MIN_BLOCKS)
        self.assertLess(floor, nb_host)
        n_mid = (floor + nb_host) // 2
        need_mid = self._big(2 ** 62, "device", freq_batch=n_mid)["device_need"]
        need_host = self._big(2 ** 62, "device", freq_batch=nb_host)["device_need"]
        self.assertLess(need_mid, need_host)
        est = self._big(int(need_mid * 1.0001 / 0.9), "auto")
        self.assertEqual(est["store_resident"], "device")
        self.assertTrue(floor <= est["nb"] < nb_host, (floor, est["nb"], nb_host))
        self.assertIn("at frequency batch {} ".format(floor), est["store_reason"])

    def test_device_only_batch_is_not_capped_by_the_host_table(self):
        """The host cap limits the batch below what the device admits:
        device_nb_store reports the device's own limit, the selected batch is
        the host's."""
        from hwave.solver import flex_bond
        peak5 = self._big(None, "host", freq_batch=5)["peak"]
        cap_gb = peak5 * 1.0001 / flex_bond._GIB
        est = self._big(2 ** 50, "auto", cap_gb=cap_gb)
        self.assertEqual(est["store_resident"], "device")
        self.assertEqual(est["nb"], 5)
        self.assertEqual(est["device_nb_store"], est["nslices"])
        self.assertEqual(est["device_nb"], est["nslices"])

    def test_pairing_redress_row(self):
        """With the in-process pairing a device store also charges the
        pairing's re-dressing (G2 plus two batch temporaries); the host-store
        table and the pairing-free table do not have the row."""
        on = self._big(2 ** 50, "device", w_in_place=False, pairing=True, freq_batch=4)
        off = self._big(2 ** 50, "device", w_in_place=False, pairing=False, freq_batch=4)
        host = self._big(2 ** 50, "host", w_in_place=False, pairing=True, freq_batch=4)
        it, nvol, ND = 16, self._BIG["nvol"], 9 * 4
        row = 2 ** 4 * nvol * 64 * it + 2 * 4 * nvol * ND * ND * it
        self.assertEqual(on["device_rows"]["pairing_redress"], row)
        self.assertNotIn("pairing_redress", off["device_rows"])
        self.assertNotIn("pairing_redress", host["device_rows"])
        self.assertAlmostEqual(on["device_need"] - off["device_need"], 1.25 * row, delta=1.0)
        # the solver passes the flag
        self.assertTrue(self.solvers["singlet"]._bond_est_kwargs["pairing"])
        self.assertFalse(self.solvers["none"]._bond_est_kwargs["pairing"])

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
        with self.assertRaisesRegex(ValueError, "needs gpu = true"):
            _flex({"IterationMax": 1, "longitudinal_bond_store": "device"})


class TestKey(unittest.TestCase):

    def _raw(self, param, gate=True):
        from hwave.solver.flex import FLEX
        p = {"flex_hartree_fock": True, "longitudinal_bond_channels": gate, "Nmat": 8,
             "gpu": True}
        p.update(param)
        return FLEX._parse_phase_b_keys({"param": p, "calc_scheme": "general"})

    def test_default_values_and_case(self):
        self.assertEqual(self._raw({})["longitudinal_bond_store"], "auto")
        for key in ("longitudinal_bond_store", "Longitudinal_Bond_Store"):
            for v, want in (("auto", "auto"), ("HOST", "host"), (" Device ", "device")):
                with self.subTest(key=key, v=v):
                    self.assertEqual(self._raw({key: v})["longitudinal_bond_store"], want)

    def test_device_needs_gpu_at_parse_time(self):
        for gpu in (False, "false", None):
            with self.subTest(gpu=gpu):
                p = {"longitudinal_bond_store": "device"}
                p["gpu"] = gpu
                if gpu is None:
                    del p["gpu"]
                    from hwave.solver.flex import FLEX
                    with self.assertRaisesRegex(ValueError, "needs gpu = true"):
                        FLEX._parse_phase_b_keys({"param": {
                            "flex_hartree_fock": True, "longitudinal_bond_channels": True,
                            "Nmat": 8, "longitudinal_bond_store": "device"},
                            "calc_scheme": "general"})
                else:
                    with self.assertRaisesRegex(ValueError, "needs gpu = true"):
                        self._raw(p)
        for v in ("auto", "host"):
            self._raw({"longitudinal_bond_store": v, "gpu": False})

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
