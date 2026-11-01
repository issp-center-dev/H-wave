"""Golden tests for the canonical digests.  The expected bytes come from a
test-side encoder (struct.pack, hashlib), never from the module under test."""
import hashlib
import struct
import unittest
from unittest import mock

import numpy as np

from hwave.solver import fingerprint as fp
from hwave.solver.so_records import compile_interaction
from tests.so_bond_oracles import COMPILER_TABLES, closed_row

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


def _ref_header(component, version):
    return _s("hwave-fp") + _i(version) + _s(component)


def _ref_records_digest(records, version):
    rows = sorted(
        (str(r["type"]).encode("utf-8"), tuple(int(x) for x in r["modes"]),
         tuple(int(x) for x in np.asarray(r["cells"]).reshape(-1)), complex(r["v"]))
        for r in records)
    h = hashlib.sha256()
    h.update(_ref_header("interaction", version))
    h.update(_i(len(rows)))
    for ty, m, c, v in rows:
        h.update(_s(ty.decode("utf-8")) + b"".join(_i(x) for x in m + c) + _c(v))
    return h.hexdigest()


def _ref_topology_digest(rows, version):
    rows = sorted(tuple(int(x) for x in r) for r in np.asarray(rows).reshape(-1, 3))
    h = hashlib.sha256()
    h.update(_ref_header("topology", version))
    h.update(_i(len(rows)))
    for r in rows:
        h.update(b"".join(_i(x) for x in r))
    return h.hexdigest()


def _plus_zero(r):
    r = r.copy()
    r["v"][0] = complex(r["v"][0].real + 0.0, r["v"][0].imag + 0.0)
    return r


def _plus_zero_shell(tables):
    t = {k: dict(v) for k, v in tables.items()}
    t["CoulombInter"].update(closed_row(3, 0, 1, 0.0))
    return t


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
        # Fortran-ordered 6-D array: merging the xyz axes needs a copy
        sc = np.asfortranarray(_to_sc(g))
        self.assertFalse(sc.flags.c_contiguous)
        self.assertFalse(np.shares_memory(sc.reshape(4, 4, 6, 8), sc))  # non-vacuity
        whole = sc.nbytes

        class Rec(np.ndarray):
            """Records the nbytes of every array that reshape/copy/astype
            materializes (i.e. that does not share memory with its source)."""
            log = []

            def _note(self, r):
                if not np.may_share_memory(r, self):
                    Rec.log.append(r.nbytes)
                return r

            def reshape(self, *a, **k):
                return self._note(super().reshape(*a, **k))

            def copy(self, *a, **k):
                return self._note(super().copy(*a, **k))

            def astype(self, *a, **k):
                return self._note(super().astype(*a, **k))

        rec = sc.view(Rec)
        for chunk in (6 * 4 * 4 * 16, 16):  # one frequency slice; below one k-row
            Rec.log.clear()
            with mock.patch.object(fp, "CHUNK_BYTES", chunk):
                self.assertEqual(fp.green_digest(rec, "sc_o_o_xyz_f"), ref)
            # no materialized temporary may approach the whole array, and
            # every one is at most max(chunk, one k-row of nd*nd complex)
            self.assertLess(whole, 8 * 6 * 4 * 4 * 16 + 1)
            for n in Rec.log:
                self.assertLess(n, whole)
                self.assertLessEqual(n, max(chunk, 4 * 4 * 16))

    def test_caller_array_not_mutated(self):
        g = _rand((1, 2, 2, 2, 2))
        g[0, 0, 0, 0, 0] = complex(-0.0, -0.0)
        before = g.copy()
        fp.green_digest(g, "b_f_k_o_o_v1")
        np.testing.assert_array_equal(np.signbit(g.view("f8")), np.signbit(before.view("f8")))

    def test_unknown_layout_refused(self):
        with self.assertRaises(ValueError):
            fp.green_digest(_rand((1, 2, 2, 2, 2)), "nope")

    def test_block_source_same_digest(self):
        """``green_digest_blocks`` equals ``green_digest`` of the
        materialized array at every chunking, and is only ever asked for
        in-range blocks of at most max(chunk, one k-row)."""
        g = _rand((1, 8, 6, 4, 4))
        ref = fp.green_digest(g, "b_f_k_o_o_v1")
        asked = []

        def block(f0, f1, k0, k1):
            asked.append((f1 - f0) * (k1 - k0) * 4 * 4 * 16)
            return g[0, f0:f1, k0:k1].copy()
        for chunk in (fp.CHUNK_BYTES, 16, 6 * 4 * 4 * 16, 3 * 6 * 4 * 4 * 16 + 16):
            asked.clear()
            with mock.patch.object(fp, "CHUNK_BYTES", chunk):
                self.assertEqual(fp.green_digest_blocks(8, 6, 4, block), ref)
            self.assertTrue(asked)
            for n in asked:
                self.assertLessEqual(n, max(chunk, 4 * 4 * 16))

    def test_block_source_shape_checked(self):
        with self.assertRaisesRegex(ValueError, "expected"):
            fp.green_digest_blocks(2, 2, 2, lambda f0, f1, k0, k1: np.zeros((1, 1, 2, 2), complex))

    @unittest.skipUnless(_HAS_CUPY, "cupy not available")
    def test_device_array_same_digest(self):
        g = _rand((1, 8, 6, 4, 4))
        self.assertEqual(fp.green_digest(cupy.asarray(g), "b_f_k_o_o_v1"),
                         fp.green_digest(g, "b_f_k_o_o_v1"))


class TestInteractionDigests(unittest.TestCase):
    def test_goldens(self):
        c = compile_interaction(COMPILER_TABLES, 2)
        self.assertEqual(fp.interaction_digest(c.records), _ref_records_digest(c.records, 2))
        self.assertEqual(fp.topology_digest(c.declared_topology),
                         _ref_topology_digest(c.declared_topology, 2))

    def test_pure_array_golden(self):
        from hwave.solver.so_records import RECORD_DTYPE
        r = np.zeros(2, dtype=RECORD_DTYPE)
        r["type"] = ["Hund", "CoulombIntra"]
        r["modes"] = [[0, 1, 1, 0], [1, 1, 1, 1]]
        r["cells"][0, 1] = [1, 0, 0]
        r["v"] = [0.5 - 0.25j, 2.0]
        self.assertEqual(fp.interaction_digest(r), _ref_records_digest(r, 2))
        self.assertEqual(fp.topology_digest(np.array([[1, 0, 0], [-1, 0, 0], [0, 0, 0]])),
                         _ref_topology_digest([[-1, 0, 0], [0, 0, 0], [1, 0, 0]], 2))

    def test_order_and_unsorted_input(self):
        c = compile_interaction(COMPILER_TABLES, 2)
        self.assertEqual(fp.interaction_digest(c.records[::-1]), fp.interaction_digest(c.records))
        rev = c.declared_topology[::-1].copy()
        self.assertEqual(fp.topology_digest(rev), fp.topology_digest(c.declared_topology))

    def test_shuffled_permutation(self):
        c = compile_interaction(COMPILER_TABLES, 2)
        rng = np.random.default_rng(1107)
        pr = rng.permutation(len(c.records))
        pt = rng.permutation(len(c.declared_topology))
        self.assertFalse((pr == np.arange(len(pr))).all())
        self.assertEqual(fp.interaction_digest(c.records[pr]), fp.interaction_digest(c.records))
        self.assertEqual(fp.topology_digest(c.declared_topology[pt]),
                         fp.topology_digest(c.declared_topology))

    def test_signed_zero_against_reference(self):
        c = compile_interaction(COMPILER_TABLES, 2)
        for part in ("real", "imag"):
            r = c.records.copy()
            z = r["v"][0]
            r["v"][0] = complex(-0.0, z.imag) if part == "real" else complex(z.real, -0.0)
            self.assertTrue(np.signbit(getattr(r["v"][0], part)))
            self.assertEqual(fp.interaction_digest(r), _ref_records_digest(r, 2))
            self.assertEqual(fp.interaction_digest(r), fp.interaction_digest(_plus_zero(r)))

    def test_non_finite_message_and_version(self):
        c = compile_interaction(COMPILER_TABLES, 2)
        self.assertNotEqual(fp.interaction_digest(c.records), _ref_records_digest(c.records, 1))
        for bad_value in (complex(np.nan, 0.0), complex(np.inf, 0.0), complex(-np.inf, 0.0),
                          complex(0.0, np.nan), complex(0.0, np.inf), complex(0.0, -np.inf)):
            bad = c.records.copy()
            bad["v"][3] = bad_value
            self.assertFalse(np.isfinite(bad["v"][3]))
            rec = bad[3]
            key = (str(rec["type"]), tuple(int(x) for x in rec["modes"]),
                   tuple(tuple(int(x) for x in row) for row in rec["cells"]))
            with self.assertRaises(ValueError) as cm:
                fp.interaction_digest(bad)
            self.assertEqual(str(cm.exception), f"non-finite value in interaction entry {key}")

    def test_zero_rows(self):
        a = compile_interaction(COMPILER_TABLES, 2)
        b = compile_interaction(_plus_zero_shell(COMPILER_TABLES), 2)
        self.assertEqual(fp.interaction_digest(a.records), fp.interaction_digest(b.records))
        self.assertNotEqual(fp.topology_digest(a.declared_topology),
                            fp.topology_digest(b.declared_topology))

    def test_version_one_unchanged(self):
        self.assertEqual(fp.FINGERPRINT_VERSION, 1)
        self.assertEqual(fp.INTERACTION_FINGERPRINT_VERSION, 2)
        self.assertEqual(fp.digest("x", [b"ab"]), _ref_digest("x", [b"ab"]))


if __name__ == "__main__":
    unittest.main()
