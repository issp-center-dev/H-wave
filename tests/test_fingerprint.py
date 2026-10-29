"""Golden tests for the canonical digests.  The expected bytes come from a
test-side encoder (struct.pack, hashlib), never from the module under test."""
import hashlib
import struct
import unittest
from unittest import mock

import numpy as np

from hwave.solver import fingerprint as fp

try:
    import cupy
    cupy.zeros(1)
    _HAS_CUPY = True
except Exception:  # pragma: no cover - no GPU
    _HAS_CUPY = False


def _s(x):
    b = x.encode("utf-8")
    return struct.pack("<I", len(b)) + b


def _b(v):
    return b"\x01" if v else b"\x00"


def _i(v):
    return struct.pack("<q", v)


def _c(z):
    return struct.pack("<dd", z.real + 0.0, z.imag + 0.0)


def _ref_digest(component, parts):
    h = hashlib.sha256()
    h.update(_s("hwave-fp") + _i(1) + _s(component))
    for p in parts:
        h.update(p)
    return h.hexdigest()


def _ref_table_digest(component, table):
    rows = sorted((tuple(k[0]) + tuple(k[1]), v) for k, v in table.items())
    parts = [_b(True), _i(len(rows))]
    for key, v in rows:
        parts.append(b"".join(_i(x) for x in key) + _c(complex(v)))
    return _ref_digest(component, parts)


def _ref_array_digest(component, a):
    a = np.asarray(a, dtype=complex)
    parts = [_s("c16"), _i(a.ndim)] + [_i(n) for n in a.shape]
    flat = a.reshape(-1)
    parts.append(b"".join(_c(z) for z in flat))
    return _ref_digest(component, parts)


def _rand(shape, seed=0):
    r = np.random.default_rng(seed)
    return r.standard_normal(shape) + 1j * r.standard_normal(shape)


def _to_sc(g):
    """(1, Nmat, nvol, nd, nd) with nvol = 6 -> (nd, nd, 3, 2, 1, Nmat)."""
    nmat, nd = g.shape[1], g.shape[3]
    return g[0].reshape(nmat, 3, 2, 1, nd, nd).transpose(4, 5, 1, 2, 3, 0)


T = {((0, 0, 0), (0, 0)): -1.0 + 0j, ((1, 0, 0), (0, 1)): 0.5 - 0.25j}
E = {((0, 0, 0), (1, 1)): 0.1 + 0j}


class TestPrimitives(unittest.TestCase):
    def test_int_real_bool_string(self):
        self.assertEqual(fp.encode_int(-3), struct.pack("<q", -3))
        self.assertEqual(fp.encode_real(-0.0), struct.pack("<d", 0.0))
        self.assertEqual(fp.encode_bool(True), b"\x01")
        self.assertEqual(fp.encode_bool(False), b"\x00")
        b = "µ-α".encode("utf-8")
        self.assertEqual(fp.encode_string("µ-α"), struct.pack("<I", len(b)) + b)

    def test_complex_signed_zero_both_parts(self):
        self.assertEqual(fp.encode_complex(complex(1.0, -0.0)), struct.pack("<dd", 1.0, 0.0))
        self.assertEqual(fp.encode_complex(complex(-0.0, 2.0)), struct.pack("<dd", 0.0, 2.0))

    def test_non_finite_refused(self):
        for v in (float("nan"), float("inf")):
            with self.assertRaisesRegex(ValueError, "non-finite"):
                fp.encode_real(v)

    def test_constants(self):
        self.assertEqual(fp.FINGERPRINT_VERSION, 1)
        self.assertEqual(fp.CHUNK_BYTES, 64 * 2**20)
        self.assertEqual(fp.PRODUCER_STAGES, (
            "rpa_so", "flex_so_onsite", "flex_so_bond", "sc_so_static",
            "sc_so_dynamic", "sc_so_bond_static", "sc_so_bond_dynamic"))


class TestDigests(unittest.TestCase):
    def test_table_layout_order_presence(self):
        t = {((1, 0, 0), (0, 1)): 0.5 - 0.25j, ((0, 0, 0), (0, 0)): -1.0 + 0j}
        self.assertEqual(fp.table_digest("transfer", t), _ref_table_digest("transfer", t))
        self.assertEqual(fp.table_digest("transfer", t),
                         fp.table_digest("transfer", dict(reversed(list(t.items())))))
        self.assertEqual(fp.table_digest("extern", None), _ref_digest("extern", [_b(False)]))
        self.assertNotEqual(fp.table_digest("extern", None), fp.table_digest("extern", {}))

    def test_table_signed_keys_sort_and_chunking(self):
        t = {((0, -1, 2), (1, 0)): 1j, ((-1, 3, 0), (0, 0)): 2 + 0j,
             ((0, -1, 2), (0, 1)): complex(-0.0, -0.0), ((0, 0, 0), (0, 0)): 1.5}
        ref = _ref_table_digest("transfer", t)
        self.assertEqual(fp.table_digest("transfer", t), ref)
        with mock.patch.object(fp, "CHUNK_BYTES", 56):
            self.assertEqual(fp.table_digest("transfer", t), ref)

    def test_table_signed_zero_canonicalized(self):
        neg = {((0, 0, 0), (0, 0)): complex(-0.0, 1.0), ((1, 0, 0), (0, 0)): complex(2.0, -0.0)}
        pos = {((0, 0, 0), (0, 0)): complex(0.0, 1.0), ((1, 0, 0), (0, 0)): complex(2.0, 0.0)}
        self.assertEqual(fp.table_digest("transfer", neg), fp.table_digest("transfer", pos))
        self.assertEqual(fp.table_digest("transfer", neg), _ref_table_digest("transfer", pos))

    def test_table_non_finite_refused(self):
        t = {((0, 0, 0), (0, 0)): complex(float("nan"), 0.0)}
        with self.assertRaisesRegex(ValueError, "non-finite.*transfer"):
            fp.table_digest("transfer", t)

    def test_domain_separation(self):
        t = {((0, 0, 0), (0, 0)): 1.0 + 0j}
        self.assertNotEqual(fp.table_digest("transfer", t), fp.table_digest("extern", t))

    def test_band_digests_reference(self):
        m = fp.BandMeta(0.25, True, "spin_block", 4, 2, (4, 4, 1), (2, 2, 1))
        d = fp.band_digests(T, E, m)
        self.assertEqual(set(d), {"fp_transfer", "fp_extern", "fp_meta", "fp_band"})
        meta = (struct.pack("<d", 0.25) + _b(True) + _s("spin_block") + _i(4) + _i(2)
                + b"".join(_i(x) for x in (4, 4, 1, 2, 2, 1)))
        self.assertEqual(d["fp_meta"], _ref_digest("meta", [meta]))
        self.assertEqual(d["fp_transfer"], _ref_table_digest("transfer", T))
        self.assertEqual(d["fp_extern"], _ref_table_digest("extern", E))
        self.assertEqual(d["fp_band"], _ref_digest(
            "band", [_s(d["fp_transfer"]), _s(d["fp_extern"]), _s(d["fp_meta"])]))
        bi = fp.band_inputs(T, E, m)
        self.assertEqual(bi.digests, d)
        self.assertIs(bi.transfer, T)
        self.assertIs(bi.meta, m)

    def test_band_digest_detects_coeff_extern(self):
        m0 = fp.BandMeta(0.0, True, "spin_block", 4, 2, (4, 4, 1), (4, 4, 1))
        self.assertNotEqual(fp.band_digests(T, E, m0)["fp_meta"],
                            fp.band_digests(T, E, m0._replace(coeff_extern=0.1))["fp_meta"])
        self.assertNotEqual(fp.band_digests(T, E, m0)["fp_band"],
                            fp.band_digests(T, None, m0)["fp_band"])

    def test_green_layouts_endianness_and_reference(self):
        g = _rand((1, 8, 6, 4, 4))
        g_sc = _to_sc(g)
        a = fp.green_digest(g, "b_f_k_o_o_v1")
        self.assertEqual(a, fp.green_digest(g_sc, "sc_o_o_xyz_f"))
        self.assertEqual(a, fp.green_digest(g.astype(">c16"), "b_f_k_o_o_v1"))
        self.assertEqual(a, _ref_array_digest("green", g[0]))

    def test_signed_zero_and_non_finite_green(self):
        g = _rand((1, 2, 2, 2, 2))
        h = g.copy()
        h[0, 0, 0, 0, 0] = complex(-0.0, -0.0)
        g[0, 0, 0, 0, 0] = 0j
        self.assertEqual(fp.green_digest(g, "b_f_k_o_o_v1"), fp.green_digest(h, "b_f_k_o_o_v1"))
        h[0, 1, 1, 1, 1] = complex(0.0, float("inf"))
        with self.assertRaisesRegex(ValueError, "non-finite"):
            fp.green_digest(h, "b_f_k_o_o_v1")

    def test_shape_header_through_module(self):
        z = np.zeros((1, 2, 3, 1, 1), complex)
        self.assertNotEqual(fp.green_digest(z, "b_f_k_o_o_v1"),
                            fp.green_digest(z.reshape(1, 3, 2, 1, 1), "b_f_k_o_o_v1"))

    def test_chunk_boundary_crossing(self):
        g = _rand((1, 8, 6, 4, 4))
        ref = fp.green_digest(g, "b_f_k_o_o_v1")
        for nbytes in (16, 6 * 4 * 4 * 16, 3 * 6 * 4 * 4 * 16 + 16):
            with mock.patch.object(fp, "CHUNK_BYTES", nbytes):
                self.assertEqual(fp.green_digest(g, "b_f_k_o_o_v1"), ref)
                self.assertEqual(fp.green_digest(_to_sc(g), "sc_o_o_xyz_f"), ref)

    def test_sc_noncontiguous_no_whole_copy(self):
        g = _rand((1, 8, 6, 4, 4))
        ref = fp.green_digest(g, "b_f_k_o_o_v1")
        # non-C-contiguous 6-D array: a strided view of a larger buffer
        big = np.zeros(_to_sc(g).shape[:-1] + (16,), complex)
        big[..., ::2] = _to_sc(g)
        sc = big[..., ::2]
        self.assertFalse(sc.flags.c_contiguous)
        whole = sc.size * 16
        seen = []
        real_host = fp.backend.to_host
        real_acg = np.ascontiguousarray

        def rec_host(a):
            seen.append(a.size * 16)
            return real_host(a)

        def rec_acg(a, *args, **kw):
            seen.append(np.asarray(a).size * 16)
            return real_acg(a, *args, **kw)

        for chunk in (fp.CHUNK_BYTES, 6 * 4 * 4 * 16, 16):
            seen.clear()
            with mock.patch.object(fp, "CHUNK_BYTES", chunk), \
                    mock.patch.object(fp.backend, "to_host", rec_host), \
                    mock.patch.object(fp.np, "ascontiguousarray", rec_acg):
                self.assertEqual(fp.green_digest(sc, "sc_o_o_xyz_f"), ref)
            self.assertTrue(seen)
            if chunk < whole:
                self.assertLess(max(seen), whole)
                self.assertLessEqual(max(seen), max(chunk, 4 * 4 * 16) * 16)

    def test_caller_array_not_mutated(self):
        g = _rand((1, 2, 2, 2, 2))
        g[0, 0, 0, 0, 0] = complex(-0.0, -0.0)
        before = g.copy()
        fp.green_digest(g, "b_f_k_o_o_v1")
        np.testing.assert_array_equal(np.signbit(g.view("f8")), np.signbit(before.view("f8")))

    def test_unknown_layout_refused(self):
        with self.assertRaises(ValueError):
            fp.green_digest(_rand((1, 2, 2, 2, 2)), "nope")

    @unittest.skipUnless(_HAS_CUPY, "cupy not available")
    def test_device_array_same_digest(self):
        g = _rand((1, 8, 6, 4, 4))
        self.assertEqual(fp.green_digest(cupy.asarray(g), "b_f_k_o_o_v1"),
                         fp.green_digest(g, "b_f_k_o_o_v1"))


if __name__ == "__main__":
    unittest.main()
