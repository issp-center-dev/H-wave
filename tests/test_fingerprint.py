"""Golden tests for the canonical digests.  The expected bytes come from a
test-side encoder (struct.pack, hashlib), never from the module under test."""
import hashlib
import inspect
import struct
import unittest
from unittest import mock

import numpy as np

from hwave.solver import fingerprint as fp
from hwave.solver.so_records import compile_interaction
from tests.so_bond_oracles import (BOND_MEMBERS, COMPILER_TABLES, EXPECTED_RPA_PRODUCER,
                                   PROVENANCE_MEMBERS, UNVERSIONED_REFUSAL_TEXT,
                                   V1_REFUSAL_TEXT, V2_ONLY_MEMBERS, closed_row)

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
    """The test-side digest at fingerprint version 2 (the version every
    spin-orbital archive is written at)."""
    h = hashlib.sha256()
    h.update(_s("hwave-fp") + _i(2) + _s(component))
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
        self.assertEqual(fp.FINGERPRINT_VERSION, 2)
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

    def test_one_version(self):
        # every digest is written at the one fingerprint version
        self.assertEqual(fp.FINGERPRINT_VERSION, 2)
        self.assertFalse(hasattr(fp, "INTERACTION_FINGERPRINT_VERSION"))
        self.assertEqual(fp.digest("x", [b"ab"]), _ref_digest("x", [b"ab"]))
        self.assertEqual(fp.interaction_digest(_COMPILED.records),
                         _ref_records_digest(_COMPILED.records, 2))



# ------------------------------------------------------ version-2 fixtures
def _ref_digest_v(component, parts, version):
    h = hashlib.sha256()
    h.update(_ref_header(component, version))
    for p in parts:
        h.update(p)
    return h.hexdigest()


def _ref_i8_array(a):
    a = np.asarray(a, dtype="<i8")
    return _s("i8") + _i(a.ndim) + b"".join(_i(n) for n in a.shape) + a.tobytes()


def _ref_resolved(delta_r, reverse, slot, version=2):
    return _ref_digest_v("resolved_topology",
                         [_ref_i8_array(delta_r), _ref_i8_array(reverse), _s(slot)], version)


META2 = fp.BandMeta(0.25, True, "spin_block", 4, 2, (4, 4, 1), (2, 2, 1))
STATE2 = fp.so_state(beta=2.0, nmat=8, mu=-0.3, density=0.8, calc_mu=True, target=0.8)
DR = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], dtype=np.int64)
REV = np.array([0, 2, 1], dtype=np.int64)
REV_SWAPPED = np.array([0, 1, 2], dtype=np.int64)
SLOT = "sr-:sc-:q-:nu+:rowrev1:colrev0:phr-:phc0"
BOND = fp.BondInputs(DR, REV, SLOT, 1)
GREEN_FP = "a" * 64
_COMPILED = compile_interaction(COMPILER_TABLES, 2)
INTER = fp.interaction_inputs(_COMPILED)

#: payload of each kind (small arrays of the recorded ranks)
_PAYLOAD = {
    "chi0q": {"chi0q": np.zeros((2, 1, 4, 4, 4, 4), complex)},
    "chiq": {"chiq": np.zeros((2, 1, 4, 4, 4, 4), complex)},
    "green": {"green": np.zeros((1, 2, 1, 4, 4), complex), "green_layout": np.str_("b_f_k_o_o_v1")},
    "sigma": {"sigma_static": np.zeros((1, 1, 1, 4, 4), complex),
              "sigma": np.zeros((1, 2, 1, 4, 4), complex)},
    "gap": {"gap": np.zeros((4, 4, 1, 1, 1), complex)},
}
_GAP_DYNAMIC = {"gap": np.zeros((4, 4, 1, 1, 1, 2), complex)}


def _v2_members(stage="rpa_so", **kw):
    args = dict(producer_stage=stage, producer_label="test producer", transfer=T, extern=E,
                meta=META2, fp_green=GREEN_FP, state=STATE2, interaction=INTER)
    if fp.STAGE_CAPABILITY[stage] == "so_bond":
        args["bond"] = BOND
    args.update(kw)
    return fp.so_provenance_members(**args)


def _archive(kind="chiq", stage="rpa_so", payload=None, **kw):
    d = dict(_v2_members(stage, **kw))
    d.update(_PAYLOAD[kind] if payload is None else payload)
    d["spin_orbital"] = np.bool_(True)
    d["index_convention"] = np.str_("spin_block")
    if stage == "sc_so_bond_static":
        d["block_layout"] = np.str_("dense")
    return d


def _malformed_text(path, check, stored, expected, label="test producer"):
    return ("'{}' is malformed or truncated: {} ({} vs {}). Regenerate it with {}."
            .format(path, check, stored, expected, label))


class TestVersionTwoDigests(unittest.TestCase):
    def test_header_version_changes_every_digest(self):
        g = _rand((1, 2, 3, 2, 2))
        digests = {
            "digest": (lambda v: fp.digest("x", [b"ab"], version=v),
                       lambda v: _ref_digest_v("x", [b"ab"], v)),
            "table": (lambda v: fp.table_digest("transfer", T, version=v), None),
            "band": (lambda v: fp.band_digests(T, E, META2, version=v)["fp_band"], None),
            "green": (lambda v: fp.green_digest(g, "b_f_k_o_o_v1", version=v), None),
            "green_blocks": (lambda v: fp.green_digest_blocks(
                2, 3, 2, lambda f0, f1, k0, k1: g[0, f0:f1, k0:k1], version=v), None),
            "interaction": (lambda v: fp.interaction_digest(_COMPILED.records, version=v),
                            lambda v: _ref_records_digest(_COMPILED.records, v)),
            "topology": (lambda v: fp.topology_digest(_COMPILED.declared_topology, version=v),
                         lambda v: _ref_topology_digest(_COMPILED.declared_topology, v)),
        }
        for name, (f, ref) in digests.items():
            with self.subTest(name=name):
                self.assertNotEqual(f(2), f(1))
                if ref is not None:
                    self.assertEqual(f(2), ref(2))
                    self.assertEqual(f(1), ref(1))
        # the band digests at version 2, from the test-side encoder
        d = fp.band_digests(T, E, META2, version=2)
        meta = (struct.pack("<d", 0.25) + _b(True) + _s("spin_block") + _i(4) + _i(2)
                + b"".join(_i(x) for x in (4, 4, 1, 2, 2, 1)))
        self.assertEqual(d["fp_meta"], _ref_digest_v("meta", [meta], 2))
        rows = sorted((tuple(k[0]) + tuple(k[1]), v) for k, v in T.items())
        parts = [_b(True), _i(len(rows))] + [b"".join(_i(x) for x in k) + _c(complex(v))
                                             for k, v in rows]
        self.assertEqual(d["fp_transfer"], _ref_digest_v("transfer", parts, 2))
        self.assertEqual(d["fp_band"], _ref_digest_v(
            "band", [_s(d["fp_transfer"]), _s(d["fp_extern"]), _s(d["fp_meta"])], 2))
        self.assertEqual(fp.green_digest(g, "b_f_k_o_o_v1", version=2),
                         _ref_digest_v("green", [_s("c16"), _i(4)] + [_i(n) for n in (2, 3, 2, 2)]
                                       + [b"".join(_c(z) for z in g[0].reshape(-1))], 2))
        # defaults: every digest at FINGERPRINT_VERSION (2)
        self.assertEqual(fp.band_digests(T, E, META2), d)
        self.assertEqual(fp.green_digest(g, "b_f_k_o_o_v1"),
                         fp.green_digest(g, "b_f_k_o_o_v1", version=2))
        self.assertEqual(fp.interaction_digest(_COMPILED.records),
                         fp.interaction_digest(_COMPILED.records, version=2))
        self.assertEqual(fp.band_inputs(T, E, META2, version=2).digests, d)

    def test_resolved_topology(self):
        self.assertEqual(fp.resolved_topology_digest(DR, REV, SLOT), _ref_resolved(DR, REV, SLOT))
        self.assertNotEqual(fp.resolved_topology_digest(DR, REV_SWAPPED, SLOT),
                            fp.resolved_topology_digest(DR, REV, SLOT))
        self.assertNotEqual(fp.resolved_topology_digest(DR, REV, SLOT, version=1),
                            fp.resolved_topology_digest(DR, REV, SLOT))
        from hwave.solver import so_bond
        self.assertEqual(so_bond.SLOT_CONVENTION.to_string(), SLOT)
        self.assertEqual(fp.resolved_topology_digest(DR, REV, so_bond.SLOT_CONVENTION),
                         _ref_resolved(DR, REV, SLOT))
        with self.assertRaisesRegex(ValueError, "delta_r"):
            fp.resolved_topology_digest(DR[:, :2], REV, SLOT)
        with self.assertRaisesRegex(ValueError, "reverse"):
            fp.resolved_topology_digest(DR, REV[:2], SLOT)

    def test_interaction_inputs(self):
        ii = fp.interaction_inputs(_COMPILED)
        self.assertEqual(ii.fp_interaction, _ref_records_digest(_COMPILED.records, 2))
        self.assertEqual(ii.fp_topology, _ref_topology_digest(_COMPILED.declared_topology, 2))
        self.assertIs(ii.records, _COMPILED.records)
        self.assertIs(ii.declared_topology, _COMPILED.declared_topology)


class TestProvenanceSchema(unittest.TestCase):
    """Construction check: table (a), literal member sets per kind and stage."""
    COMMON = {"fingerprint_version", "producer_stage", "producer_label", "producer_version",
              "fp_transfer", "fp_extern", "fp_meta", "fp_band", "fp_green", "fp_interaction",
              "fp_topology", "fp_resolved_topology", "fp_capability", "prov_transfer_keys",
              "prov_transfer_vals", "prov_extern_present", "prov_extern_keys", "prov_extern_vals",
              "prov_coeff_extern", "prov_enable_spin_orbital", "prov_index_convention",
              "prov_norb", "prov_norb_phys", "prov_cell_shape", "prov_sub_shape",
              "prov_rec_types", "prov_rec_modes", "prov_rec_cells", "prov_rec_v",
              "prov_topology", "state_beta", "state_nmat", "state_mu", "state_density",
              "state_calc_mu", "state_target", "state_mu_tol"}
    BOND = {"delta_r", "reverse", "slot_convention", "max_shells"}
    KINDS = {"chi0q": ("rpa_so", "flex_so_onsite"), "chiq": ("rpa_so", "flex_so_onsite"),
             "green": ("flex_so_onsite",), "sigma": ("flex_so_onsite",),
             "gap": ("sc_so_static", "sc_so_dynamic", "sc_so_bond_static")}

    def test_sets(self):
        self.assertEqual(len(self.COMMON), 37)
        self.assertEqual({k: tuple(v) for k, v in fp.KIND_STAGES.items()}, self.KINDS)
        for kind, stages in self.KINDS.items():
            for st in stages:
                with self.subTest(kind=kind, stage=st):
                    want = self.COMMON | (self.BOND if st == "sc_so_bond_static" else set())
                    self.assertEqual(fp.required_members(kind, st), frozenset(want))
                    self.assertEqual(set(_v2_members(st)), want)
        self.assertEqual(set(_v2_members(fp_green_source="b" * 64)),
                         self.COMMON | {"fp_green_source"})
        with self.assertRaises(ValueError) as cm:
            fp.required_members("green", "rpa_so")
        self.assertEqual(str(cm.exception),
                         "artifact kind 'green' has no producer_stage 'rpa_so' "
                         "(its stages: flex_so_onsite)")

    def test_writer_values(self):
        m = _v2_members()
        self.assertEqual(m["fingerprint_version"], 2)
        self.assertEqual(m["fingerprint_version"].dtype, np.int64)
        self.assertEqual(str(m["fp_capability"]), "so_onsite")
        self.assertEqual(str(m["fp_resolved_topology"]), "none:onsite")
        self.assertEqual(str(m["fp_interaction"]), INTER.fp_interaction)
        self.assertEqual(str(m["fp_topology"]), INTER.fp_topology)
        self.assertEqual(str(m["fp_band"]), fp.band_digests(T, E, META2, version=2)["fp_band"])
        n = len(_COMPILED.records)
        self.assertEqual(m["prov_rec_types"].shape, (n,))
        self.assertEqual(m["prov_rec_types"].dtype.kind, "U")
        self.assertEqual(m["prov_rec_modes"].shape, (n, 4))
        self.assertEqual(m["prov_rec_cells"].shape, (n, 4, 3))
        self.assertEqual(m["prov_rec_v"].dtype, np.complex128)
        np.testing.assert_array_equal(m["prov_topology"], _COMPILED.declared_topology)
        b = _v2_members("sc_so_bond_static")
        self.assertEqual(str(b["fp_capability"]), "so_bond")
        self.assertEqual(str(b["fp_resolved_topology"]), _ref_resolved(DR, REV, SLOT))
        self.assertEqual(str(b["slot_convention"]), SLOT)
        self.assertEqual(int(b["max_shells"]), 1)
        np.testing.assert_array_equal(b["delta_r"], DR)
        self.assertEqual(tuple(m), PROVENANCE_MEMBERS)
        self.assertEqual(tuple(b), PROVENANCE_MEMBERS + BOND_MEMBERS)
        # the writer takes no version: it writes the one version
        self.assertNotIn("version",
                         inspect.signature(fp.so_provenance_members).parameters)

    def test_writer_refusals(self):
        cases = (
            (dict(stage="sc_so_bond_static", bond=None),
             "producer_stage 'sc_so_bond_static' needs the resolved bond topology"),
            (dict(stage="rpa_so", bond=BOND),
             "producer_stage 'rpa_so' is not bond-resolved; it takes no resolved bond topology"),
            (dict(stage="rpa_so", interaction=None),
             "the spin-orbital provenance needs the interaction fingerprints (interaction=)"),
            (dict(stage="rpa_so", fp_green_source=GREEN_FP),
             "fp_green_source must differ from fp_green"),
        )
        for kw, text in cases:
            with self.subTest(text=text):
                stage = kw.pop("stage")
                with self.assertRaises(ValueError) as cm:
                    _v2_members(stage, **kw)
                self.assertEqual(str(cm.exception), text)

    def _malformed(self, name, bad, check, stored, expected, stage="rpa_so", kind="chiq"):
        d = _archive(kind, stage)
        d[name] = bad
        with self.assertRaises(fp.SOArchiveError) as cm:
            fp.validate_structure(d, "a.npz", kind=kind)
        self.assertEqual(str(cm.exception), _malformed_text("a.npz", check, stored, expected))

    def _rows(self, rows):
        for name, bad, check, stored, expected in rows:
            with self.subTest(name=name, check=check):
                self._malformed(name, bad, check, stored, expected)

    def test_shape_rules(self):
        """dtype, rank and shape rows of the provenance schema."""
        n = len(_COMPILED.records)
        self._rows((
            ("fp_capability", np.array(["so_onsite"]),
             "member fp_capability has the wrong dtype or rank", "<U9/1d", "str/0d"),
            ("prov_rec_modes", np.zeros((3, 3), np.int64),
             "member prov_rec_modes has the wrong shape", "(3, 3)", str((n, 4))),
            ("prov_rec_cells", np.zeros((3, 4), np.int64),
             "member prov_rec_cells has the wrong dtype or rank", "int64/2d", "int64/3d"),
            ("prov_rec_v", np.zeros(2, complex),
             "member prov_rec_v has the wrong shape", "(2,)", str((n,))),
            ("prov_rec_types", np.array([1, 2, 3]),
             "member prov_rec_types has the wrong dtype or rank", "int64/1d", "str/1d"),
            ("prov_topology", np.zeros((2, 2), np.int64),
             "member prov_topology has the wrong shape", "(2, 2)", "(2, 3)"),
            ("fingerprint_version", np.float64(2.0),
             "fingerprint_version is not an integer scalar", "float64/0d", "int64/0d"),
        ))

    def test_value_rules(self):
        """Value rows: capability, sentinel, digests, finiteness, duplicated fields."""
        n = len(_COMPILED.records)
        self._rows((
            ("fp_capability", np.str_("so_bond"),
             "fp_capability does not match producer_stage 'rpa_so'", "'so_bond'", "'so_onsite'"),
            ("fp_capability", np.str_("bond"),
             "fp_capability is not a known capability", "'bond'", "one of so_onsite, so_bond"),
            ("fp_resolved_topology", np.str_("0" * 64),
             "fp_resolved_topology of a so_onsite stage", repr("0" * 64), "'none:onsite'"),
            ("fp_interaction", np.str_("A" * 64),
             "fp_interaction is not a 64-character lowercase hex digest", repr("A" * 64),
             "64 hex characters"),
        ))
        v = np.full(n, 1.0 + 0j)
        v[1] = np.nan
        self._malformed("prov_rec_v", v, "member prov_rec_v is not finite",
                        "entry 1: (nan+0j)", "finite")
        self._malformed("fp_green_source", np.str_(GREEN_FP),
                        "fp_green_source equals fp_green", repr(GREEN_FP), "a different digest")
        self._malformed("delta_r", DR, "member delta_r is not written by a so_onsite stage",
                        "present", "absent")
        self._malformed("index_convention", np.str_("interleaved"),
                        "index_convention differs from prov_index_convention",
                        "'interleaved'", "'spin_block'")

    def test_bond_rules(self):
        kw = dict(stage="sc_so_bond_static", kind="gap")
        self._malformed("fp_resolved_topology", np.str_("none:onsite"),
                        "fp_resolved_topology is not a 64-character lowercase hex digest",
                        "'none:onsite'", "64 hex characters", **kw)
        self._malformed("slot_convention", np.str_("sr+"), "slot_convention does not parse",
                        "'sr+'", "the slot-convention grammar", **kw)
        other = "sr+:sc-:q-:nu+:rowrev1:colrev0:phr-:phc0"
        self._malformed("slot_convention", np.str_(other),
                        "slot_convention is not this version's", repr(other), repr(SLOT), **kw)
        self._malformed("reverse", np.zeros(2, np.int64), "member reverse has the wrong shape",
                        "(2,)", "(3,)", **kw)
        d = _archive("gap", "sc_so_bond_static")
        del d["max_shells"]
        with self.assertRaises(fp.SOArchiveError) as cm:
            fp.validate_structure(d, "a.npz", kind="gap")
        self.assertEqual(str(cm.exception), _malformed_text(
            "a.npz", "member max_shells is missing", "absent", "present"))
        self.assertEqual(fp.validate_structure(_archive("gap", "sc_so_bond_static"), "a.npz",
                                               kind="gap"), ("gap", "sc_so_bond_static", 2))

    def test_valid_archives(self):
        for kind, stages in self.KINDS.items():
            for st in stages:
                if st == "sc_so_bond_static":
                    continue
                with self.subTest(kind=kind, stage=st):
                    payload = _GAP_DYNAMIC if st == "sc_so_dynamic" else None
                    d = _archive(kind, st, payload=payload)
                    self.assertEqual(fp.validate_structure(d, "a.npz", kind=kind), (kind, st, 2))


class _NoPayload(dict):
    """A member mapping whose payload members must never be read."""
    def __getitem__(self, k):
        if k in ("chi0q", "chiq", "green", "sigma_static", "gap"):
            raise AssertionError("payload read")
        return super().__getitem__(k)


class TestPayloadTables(unittest.TestCase):
    def setUp(self):
        import tempfile, shutil
        self.d = tempfile.mkdtemp(prefix="fp_payload_")
        self.addCleanup(shutil.rmtree, self.d, True)

    def _save(self, d, name="a.npz"):
        import os
        path = os.path.join(self.d, name)
        np.savez(path, **d)
        return path

    def test_payload_names(self):
        d = _archive("sigma", "flex_so_onsite")
        del d["sigma_static"]
        with self.assertRaises(fp.SOArchiveError) as cm:
            fp.validate_structure(d, "sigma.npz", kind="sigma")
        self.assertEqual(str(cm.exception), _malformed_text(
            "sigma.npz", "payload member sigma_static of a sigma archive is missing",
            "absent", "present"))
        d = _archive("green", "flex_so_onsite")
        d["green_layout"] = np.str_("sc_o_o_xyz_f")
        with self.assertRaises(fp.SOArchiveError) as cm:
            fp.validate_structure(d, "green.npz", kind="green")
        self.assertEqual(str(cm.exception), _malformed_text(
            "green.npz", "green_layout is not the canonical layout", "'sc_o_o_xyz_f'",
            "'b_f_k_o_o_v1'"))
        d = _archive("chiq")
        d["spin_orbital"] = np.bool_(False)
        with self.assertRaises(fp.SOArchiveError) as cm:
            fp.validate_structure(d, "chiq.npz", kind="chiq")
        self.assertEqual(str(cm.exception), _malformed_text(
            "chiq.npz", "marker spin_orbital is not true", "False", "True"))
        d = _archive("chiq")
        del d["index_convention"]
        with self.assertRaises(fp.SOArchiveError) as cm:
            fp.validate_structure(d, "chiq.npz", kind="chiq")
        self.assertEqual(str(cm.exception), _malformed_text(
            "chiq.npz", "marker index_convention is missing", "absent", "present"))
        # the spin_orbital marker and a green archive's green_layout are required
        for kind, stage, member, check in (
                ("chiq", "rpa_so", "spin_orbital", "marker spin_orbital is missing"),
                ("green", "flex_so_onsite", "spin_orbital", "marker spin_orbital is missing"),
                ("green", "flex_so_onsite", "green_layout",
                 "payload member green_layout of a green archive is missing")):
            with self.subTest(kind=kind, member=member):
                d = _archive(kind, stage)
                del d[member]
                with self.assertRaises(fp.SOArchiveError) as cm:
                    fp.validate_structure(d, "x.npz", kind=kind)
                self.assertEqual(str(cm.exception), _malformed_text(
                    "x.npz", check, "absent", "present"))
        self.assertEqual(fp.PAYLOAD_MEMBERS["green"], ("green", "green_layout"))

    def test_header_only(self):
        import zipfile
        cases = (
            ("green", {"green": np.zeros((1, 2, 1, 4, 4)),
                       "green_layout": np.str_("b_f_k_o_o_v1")},
             "payload member green has the wrong dtype", "float64", "complex128"),
            ("green", {"green": np.zeros((2, 1, 4, 4), complex),
                       "green_layout": np.str_("b_f_k_o_o_v1")},
             "payload member green has the wrong rank", "4", "5"),
            ("green", {"green": np.asfortranarray(np.zeros((1, 2, 1, 4, 4), complex)),
                       "green_layout": np.str_("b_f_k_o_o_v1")},
             "payload member green is in Fortran order", "Fortran", "C"),
            ("chi0q", {"chi0q": np.zeros((2, 1, 4, 4, 4), complex)},
             "payload member chi0q has the wrong rank", "5", "4 or 6"),
        )
        real_open = zipfile.ZipFile.open
        for kind, payload, check, stored, expected in cases:
            with self.subTest(check=check, stored=stored):
                path = self._save(_archive(kind, "flex_so_onsite", payload=payload))
                pname = kind + ".npy"
                with zipfile.ZipFile(path) as z:
                    raw = z.read(pname)
                hlen = 10 + int.from_bytes(raw[8:10], "little")
                self.assertEqual(raw[6:8], b"\x01\x00")
                self.assertGreater(len(raw), hlen)        # the data section is there

                class Spy:
                    def __init__(self, fh):
                        self.fh, self.pos = fh, 0

                    def read(self, n=-1):
                        if n is None or n < 0 or self.pos + n > hlen:
                            raise AssertionError("payload data read")
                        b = self.fh.read(n)
                        self.pos += len(b)
                        return b

                    def tell(self):
                        return self.pos

                    def __enter__(self):
                        return self

                    def __exit__(self, *a):
                        self.fh.close()

                    def close(self):
                        self.fh.close()

                def spy_open(zself, name, *a, **k):
                    fh = real_open(zself, name, *a, **k)
                    return Spy(fh) if getattr(name, "filename", name) == pname else fh

                mapping = _NoPayload(_archive(kind, "flex_so_onsite", payload=payload))
                with mock.patch.object(zipfile.ZipFile, "open", spy_open):
                    with np.load(path) as npz:
                        for data, kw in ((mapping, {"archive": path}), (npz, {})):
                            with self.assertRaises(fp.SOArchiveError) as cm:
                                fp.validate_structure(data, "a.npz", kind=kind, **kw)
                            self.assertEqual(str(cm.exception),
                                             _malformed_text("a.npz", check, stored, expected))
        # the reduced (rank-4) spin-orbital chi0q passes: its layout refusal is the loader's
        path = self._save(_archive("chi0q", payload={"chi0q": np.zeros((2, 1, 16, 16), complex)}))
        with np.load(path) as npz:
            self.assertEqual(fp.validate_structure(npz, "a.npz", kind="chi0q"),
                             ("chi0q", "rpa_so", 2))
        # every layout a spin-orbital RPA writes passes (spin-diagonal bands
        # carry a leading spin axis: ranks 5 and 7); its layout refusal is the loader's
        for shape in ((2, 1, 16, 16), (2, 2, 1, 4, 4), (2, 1, 4, 4, 4, 4),
                      (2, 2, 1, 2, 2, 2, 2)):
            with self.subTest(rpa_chi0q_rank=len(shape)):
                path = self._save(_archive("chi0q", payload={"chi0q": np.zeros(shape, complex)}))
                with np.load(path) as npz:
                    self.assertEqual(fp.validate_structure(npz, "a.npz", kind="chi0q"),
                                     ("chi0q", "rpa_so", 2))
        self.assertEqual(fp.PAYLOAD_RANKS[("chi0q", "rpa_so")], (4, 5, 6, 7))
        self.assertEqual(fp.PAYLOAD_RANKS[("chiq", "rpa_so")], (4, 6))
        # a big-endian payload is complex128
        path = self._save(_archive("chiq", payload={"chiq": np.zeros((2, 1, 4, 4, 4, 4), ">c16")}))
        with np.load(path) as npz:
            self.assertEqual(fp.validate_structure(npz, "a.npz", kind="chiq"),
                             ("chiq", "rpa_so", 2))

    def test_member_size_matches_the_header(self):
        # a payload member whose stored size is not header + prod(shape) *
        # itemsize (truncated, or with trailing bytes) is refused from the
        # zip directory and the header alone
        import zipfile
        cases = (("chi0q", "rpa_so", "chi0q"), ("chiq", "rpa_so", "chiq"),
                 ("green", "flex_so_onsite", "green"), ("sigma", "flex_so_onsite", "sigma_static"),
                 ("gap", "sc_so_static", "gap"))
        for kind, stage, pname in cases:
            for cut, extra in ((16, b""), (0, b"\0" * 16)):
                with self.subTest(kind=kind, truncated=bool(cut)):
                    src = self._save(_archive(kind, stage))
                    path = self._save({}, name="bad.npz")
                    with zipfile.ZipFile(src) as zin, zipfile.ZipFile(path, "w") as zout:
                        for info in zin.infolist():
                            raw = zin.read(info.filename)
                            if info.filename == pname + ".npy":
                                hlen = 10 + int.from_bytes(raw[8:10], "little")
                                good = len(raw)
                                raw = raw[:len(raw) - cut] + extra
                            zout.writestr(info.filename, raw)
                    self.assertEqual(good - hlen, _PAYLOAD[kind][pname].nbytes)
                    want = _malformed_text(
                        "a.npz", "payload member {} has the wrong size".format(pname),
                        str(good - cut + len(extra)), str(good))
                    mapping = _NoPayload(_archive(kind, stage))
                    with np.load(path) as npz:
                        for data, kw in ((mapping, {"archive": path}), (npz, {})):
                            with self.assertRaises(fp.SOArchiveError) as cm:
                                fp.validate_structure(data, "a.npz", kind=kind, **kw)
                            self.assertEqual(str(cm.exception), want)
                    # the intact archive passes
                    with np.load(src) as npz:
                        self.assertEqual(fp.validate_structure(npz, "a.npz", kind=kind),
                                         (kind, stage, 2))

    def test_member_size_in_a_compressed_archive(self):
        # the stored size is the uncompressed one: a compressed archive
        # passes, and a truncated member of one is refused
        import os
        import zipfile
        src = os.path.join(self.d, "c.npz")
        np.savez_compressed(src, **_archive("chiq"))
        with zipfile.ZipFile(src) as z:
            self.assertEqual(z.getinfo("chiq.npy").compress_type, zipfile.ZIP_DEFLATED)
        with np.load(src) as npz:
            self.assertEqual(fp.validate_structure(npz, "c.npz", kind="chiq"),
                             ("chiq", "rpa_so", 2))
        bad = os.path.join(self.d, "bad.npz")
        with zipfile.ZipFile(src) as zin, \
                zipfile.ZipFile(bad, "w", compression=zipfile.ZIP_DEFLATED) as zout:
            for info in zin.infolist():
                raw = zin.read(info.filename)
                if info.filename == "chiq.npy":
                    good = len(raw)
                    raw = raw[:-16]
                zout.writestr(info.filename, raw)
        with np.load(bad) as npz:
            with self.assertRaises(fp.SOArchiveError) as cm:
                fp.validate_structure(npz, "c.npz", kind="chiq")
        self.assertEqual(str(cm.exception), _malformed_text(
            "c.npz", "payload member chiq has the wrong size", str(good - 16), str(good)))

    def test_mapping_skips_headers(self):
        d = _NoPayload(_archive("green", "flex_so_onsite",
                               payload={"green": np.zeros(3),
                                        "green_layout": np.str_("b_f_k_o_o_v1")}))
        self.assertEqual(fp.validate_structure(d, "a.npz", kind="green"),
                         ("green", "flex_so_onsite", 2))

    def test_kind_required(self):
        # the caller names the archive kind; it is never inferred
        for f in (fp.validate_structure, fp.check_so_archive):
            with self.subTest(f=f.__name__):
                p = inspect.signature(f).parameters["kind"]
                self.assertIs(p.default, inspect.Parameter.empty)
                self.assertEqual(p.kind, inspect.Parameter.KEYWORD_ONLY)
        with self.assertRaises(TypeError):
            fp.validate_structure(_archive("sigma", "flex_so_onsite"), "s.npz")
        with self.assertRaises(fp.SOArchiveError) as cm:
            fp.validate_structure(_archive("chiq", "rpa_so"), "x.npz", kind="green")
        self.assertEqual(str(cm.exception), _malformed_text(
            "x.npz", "producer_stage does not write a green archive", "'rpa_so'",
            "one of flex_so_onsite"))
        with self.assertRaises(ValueError) as cm:
            fp.validate_structure(_archive(), "x.npz", kind="Chiq")
        self.assertEqual(str(cm.exception), "unknown artifact kind 'Chiq'")

    def test_deferred_stages(self):
        texts = {
            "sc_so_bond_dynamic": "'g.npz' was written by producer_stage 'sc_so_bond_dynamic', "
                                  "which this version does not read yet (the dynamic "
                                  "spin-orbital bond path arrives in a later version).",
            "flex_so_bond": "'g.npz' was written by producer_stage 'flex_so_bond', which this "
                            "version does not read yet (the spin-orbital bond FLEX arrives in "
                            "a later version).",
        }
        for stage, text in texts.items():
            with self.subTest(stage=stage):
                d = _archive("gap", "sc_so_static")
                d["producer_stage"] = np.str_(stage)
                with self.assertRaises(fp.SOArchiveError) as cm:
                    fp.validate_structure(d, "g.npz", kind="gap")
                self.assertEqual(str(cm.exception), text)
        d = _archive("gap", "sc_so_static")
        d["producer_stage"] = np.str_("sc_so_later")
        with self.assertRaises(fp.SOArchiveError) as cm:
            fp.validate_structure(d, "g.npz", kind="gap")
        self.assertEqual(str(cm.exception), _malformed_text(
            "g.npz", "producer_stage is not a known stage", "'sc_so_later'",
            "one of " + ", ".join(fp.PRODUCER_STAGES)))


def _records(*rows):
    from hwave.solver.so_records import RECORD_DTYPE
    r = np.zeros(len(rows), dtype=RECORD_DTYPE)
    for i, (ty, modes, cell1, v) in enumerate(rows):
        r[i]["type"] = ty
        r[i]["modes"] = modes
        r[i]["cells"][1] = cell1
        r[i]["v"] = v
    return r


class _Compiled:
    def __init__(self, records):
        self.records = records
        self.declared_topology = np.array([[-1, 0, 0], [0, 0, 0], [1, 0, 0]], np.int64)


class TestVersionTwoChecks(unittest.TestCase):
    RUN_BAND = fp.band_inputs(T, E, META2)        # this run's tables
    RUN_STATE = fp.SOState(2.0, 8, float("nan"), float("nan"), True, 0.8, float("nan"))

    def check(self, d, path="chiq.npz", **kw):
        args = dict(kind="chiq", run_band=self.RUN_BAND, run_state=self.RUN_STATE,
                    nmat_rule="exact", run_interaction=INTER, need="so_onsite")
        args.update(kw)
        return fp.check_so_archive(d, path, **args)

    def test_accepted(self):
        info = self.check(_archive())
        self.assertEqual(info.version, 2)
        self.assertEqual(info.producer_fp_green, GREEN_FP)
        self.assertEqual(info.producer_label, "test producer")
        self.assertEqual(info.producer_state, fp._file_state(_archive()))

    def test_prov_table_disagreeing_is_not_refused(self):
        for name in ("prov_transfer_vals", "prov_rec_v"):
            with self.subTest(name=name):
                d = _archive()
                d[name] = d[name] + 1.0
                self.assertEqual(self.check(d).version, 2)

    def test_band_digest_alone_differs(self):
        """A file whose fp_band disagrees while its three component digests
        agree with the run's is refused (the check fails closed)."""
        d = _archive()
        d["fp_band"] = np.str_("0" * 64)
        run_band = fp.band_digests(T, E, META2, version=2)["fp_band"]
        with self.assertRaises(fp.SOArchiveError) as cm:
            self.check(d)
        self.assertEqual(str(cm.exception),
                         "'chiq.npz' (written by test producer) does not match this run: band "
                         "differs: fp_band: file 000000000000, run {} (fp_transfer, fp_extern "
                         "and fp_meta agree; the band digest is inconsistent). Archives are "
                         "reused only at the identical state point (not a warm start across a "
                         "sweep); re-run test producer with this input (the same "
                         "Transfer/Extern/coeff_extern).".format(run_band[:12]))

    def test_band_and_state_rules(self):
        T2 = dict(T)
        T2[((0, 0, 0), (0, 0))] = -1.5 + 0j
        with self.assertRaises(fp.SOArchiveError) as cm:
            self.check(_archive(), run_band=fp.band_inputs(T2, E, META2))
        self.assertEqual(str(cm.exception),
                         "'chiq.npz' (written by test producer) does not match this run: band "
                         "differs: Transfer: file 2 rows, run 2 rows; first differing row "
                         "(R=(0, 0, 0), a=1, b=1): file (-1+0j), run (-1.5+0j). Archives are "
                         "reused only at the identical state point (not a warm start across a "
                         "sweep); re-run test producer with this input (the same "
                         "Transfer/Extern/coeff_extern).")
        with self.assertRaises(fp.SOArchiveError) as cm:
            self.check(_archive(), run_state=self.RUN_STATE._replace(beta=4.0))
        self.assertEqual(str(cm.exception),
                         "'chiq.npz' (written by test producer) does not match this run: state "
                         "differs: beta 2.0 (file) vs 4.0 (run). Archives are reused only at "
                         "the identical state point; re-run test producer at this filling and "
                         "temperature.")
        g = _rand((1, 8, 1, 4, 4))
        d = _archive(fp_green=fp.green_digest(g, "b_f_k_o_o_v1", version=2))
        self.assertEqual(self.check(d, green=g, green_layout="b_f_k_o_o_v1").version, 2)
        v1_fp = fp.green_digest(g, "b_f_k_o_o_v1", version=1)
        with self.assertRaises(fp.SOArchiveError) as cm:
            self.check(_archive(fp_green=v1_fp), green=g, green_layout="b_f_k_o_o_v1")
        self.assertEqual(str(cm.exception),
                         "'chiq.npz' (written by test producer): fp_green mismatch: the Green "
                         "function this run loaded has digest {}..., the archive records {}... "
                         "(both at fingerprint_version 2); the Green function and the "
                         "susceptibility do not come from the same run. Regenerate both with "
                         "test producer.".format(
                             fp.green_digest(g, "b_f_k_o_o_v1", version=2)[:12], v1_fp[:12]))

    def test_interaction_mismatch_diff(self):
        file_rec = _records(("Exchange", (0, 1, 1, 0), (1, 0, 0), 0.29),
                            ("Hund", (0, 1, 1, 0), (0, 0, 0), 0.1))
        run_rec = file_rec.copy()
        run_rec[0]["v"] = 0.31
        d = _archive(interaction=fp.interaction_inputs(_Compiled(file_rec)))
        with self.assertRaises(fp.SOArchiveError) as cm:
            self.check(d, run_interaction=fp.interaction_inputs(_Compiled(run_rec)))
        self.assertEqual(str(cm.exception),
                         "'chiq.npz' (written by test producer) does not match this run: "
                         "interaction differs: file 2 records, run 2 records; first differing "
                         "record Exchange modes (0, 1, 1, 0) cells ((0, 0, 0), (1, 0, 0), "
                         "(0, 0, 0), (0, 0, 0)): file (0.29+0j), run (0.31+0j). Archives are "
                         "reused only at the identical state point (not a warm start across a "
                         "sweep); re-run test producer with this input (the same interaction "
                         "definitions).")
        with self.assertRaises(fp.SOArchiveError) as cm:
            self.check(d, run_interaction=fp.interaction_inputs(_Compiled(file_rec[1:])))
        self.assertIn("interaction differs: file 2 records, run 1 records; first differing "
                      "record Exchange modes (0, 1, 1, 0) cells ((0, 0, 0), (1, 0, 0), "
                      "(0, 0, 0), (0, 0, 0)): file (0.29+0j), run absent.", str(cm.exception))

    def test_capability_message(self):
        with self.assertRaises(fp.SOArchiveError) as cm:
            self.check(_archive(), need="so_bond")
        self.assertEqual(str(cm.exception), "chiq.npz was produced by test producer "
                                            "(so_onsite), but this run needs a so_bond artifact")

    def test_bond_consumer_topology(self):
        d = _archive("gap", "sc_so_bond_static")
        kw = dict(path="gap.npz", kind="gap", need="so_bond")
        self.assertEqual(self.check(d, run_bond=BOND, **kw).version, 2)
        with self.assertRaises(fp.SOArchiveError) as cm:
            self.check(d, run_bond=BOND._replace(reverse=REV_SWAPPED), **kw)
        self.assertEqual(str(cm.exception),
                         "'gap.npz' (written by test producer) does not match this run: resolved "
                         "bond topology differs: file delta_r [[0, 0, 0], [1, 0, 0], [-1, 0, 0]] "
                         "reverse [0, 2, 1] (max_shells 1), run delta_r [[0, 0, 0], [1, 0, 0], "
                         "[-1, 0, 0]] reverse [0, 1, 2] (max_shells 1). Archives are reused only "
                         "at the identical state point (not a warm start across a sweep); re-run "
                         "test producer with this input (the same bond settings).")
        # another shell cutoff keeping the same channels is compatible
        self.assertEqual(self.check(d, run_bond=BOND._replace(max_shells=2), **kw).version, 2)
        other = fp.interaction_inputs(_Compiled(_COMPILED.records))
        self.assertEqual(other.fp_interaction, INTER.fp_interaction)
        with self.assertRaises(fp.SOArchiveError) as cm:
            self.check(d, run_bond=BOND, run_interaction=other, **kw)
        self.assertEqual(str(cm.exception),
                         "'gap.npz' (written by test producer) does not match this run: declared "
                         "topology differs: file {} displacements, run 3 displacements. Archives "
                         "are reused only at the identical state point (not a warm start across a "
                         "sweep); re-run test producer with this input (the same interaction "
                         "definitions).".format(len(_COMPILED.declared_topology)))
        with self.assertRaises(ValueError) as cm:
            self.check(d, **kw)
        self.assertEqual(str(cm.exception), "a so_bond consumer needs run_bond")
        # an on-site consumer of a bond-resolved archive
        with self.assertRaises(fp.SOArchiveError) as cm:
            self.check(d, path="gap.npz", kind="gap")
        self.assertEqual(str(cm.exception), "gap.npz was produced by test producer "
                                            "(so_bond), but this run needs a so_onsite artifact")

    def test_caller_arguments(self):
        with self.assertRaises(ValueError) as cm:
            self.check(_archive(), run_interaction=None)
        self.assertEqual(str(cm.exception), "check_so_archive needs run_interaction")
        with self.assertRaises(ValueError) as cm:
            self.check(_archive(), need=None)
        self.assertEqual(str(cm.exception),
                         "need must be one of ('so_onsite', 'so_bond'), got None")

    def test_expected_producer_names_an_unreadable_label(self):
        d = _archive()
        d["producer_label"] = np.str_("")
        d["fp_band"] = np.str_("z" * 64)
        with self.assertRaises(fp.SOArchiveError) as cm:
            self.check(d, expected_producer="the spin-orbital RPA run")
        self.assertTrue(str(cm.exception).endswith("Regenerate it with the spin-orbital RPA run."))

    def test_old_versions_refused(self):
        """A version-1 archive (without the interaction fingerprints) and an
        unversioned one are refused before any other check, naming the
        producer to re-run; another version is malformed."""
        v1 = {k: v for k, v in _archive().items() if k not in V2_ONLY_MEMBERS}
        v1["fingerprint_version"] = np.int64(1)
        v1["fp_band"] = np.str_("z")                  # never reached
        unversioned = {k: v for k, v in _archive().items()
                       if k not in PROVENANCE_MEMBERS}
        v1_text = V1_REFUSAL_TEXT.replace("{path}", "chiq.npz").replace("{producer}", "{}")
        old_text = UNVERSIONED_REFUSAL_TEXT.replace("{path}", "chiq.npz").replace(
            "{producer}", "{}")
        expected = EXPECTED_RPA_PRODUCER
        unreadable = dict(v1, producer_label=np.str_(" "))
        cases = (
            (v1, {}, v1_text.format("test producer")),
            (v1, {"expected_producer": expected}, v1_text.format("test producer")),
            (unreadable, {"expected_producer": expected}, v1_text.format(expected)),
            (unreadable, {}, v1_text.format(fp.UNKNOWN_PRODUCER)),
            (unversioned, {"expected_producer": expected}, old_text.format(expected)),
            (unversioned, {}, old_text.format(fp.UNKNOWN_PRODUCER)),
        )
        for d, kw, text in cases:
            for name, call in (("check_so_archive", lambda: self.check(_NoPayload(d), **kw)),
                               ("validate_structure", lambda: fp.validate_structure(
                                   _NoPayload(d), "chiq.npz", kind="chiq", **kw))):
                with self.subTest(text=text[:60], kw=kw, via=name):
                    with self.assertRaises(fp.SOArchiveError) as cm:
                        call()
                    self.assertEqual(str(cm.exception), text)
        d = _archive()
        d["fingerprint_version"] = np.int64(3)
        with self.assertRaises(fp.SOArchiveError) as cm:
            self.check(d)
        self.assertEqual(str(cm.exception), _malformed_text(
            "chiq.npz", "unsupported fingerprint_version", 3, 2))

if __name__ == "__main__":
    unittest.main()
