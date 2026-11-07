"""The layout members of bond-resolved spin-orbital archives.

* A bond-resolved archive records its slot layout in ``block_layout``
  (``"dense"`` or ``"sz_packed"``); a ``sc_so_bond_static`` gap carries it
  (the solve layout; the gap itself is always stored full) and no member
  with the reserved prefix ``sz_sector_``. An on-site stage carries none of
  these members. Other unknown members are accepted (forward compatible).
* The slot-archive machinery (``check_layout_members`` with a slot kind and
  ``check_slot_payloads``) is admitted for no archive kind in this version;
  it is exercised here through a test-local rules table naming a synthetic
  slot kind ``bondtest`` on the ``flex_so_bond`` stage, on synthetic
  archives written by the builder below, both as zip files and as plain
  mappings (whose arrays are checked, never skipped).
"""
import os
import tempfile
import types
import unittest
import zipfile

import numpy as np

from hwave.solver import fingerprint as fp
from hwave.solver import so_bond
from hwave.solver.fingerprint import SlotRules
from hwave.solver.so_records import RECORD_DTYPE

PATH = "slots.npz"
LABEL = "test producer"
GREEN_FP = "a" * 64
SLOT = "sr-:sc-:q-:nu+:rowrev1:colrev0:phr-:phc0"
MARKERS = {"spin_orbital", "index_convention"}

TEST_RULES = SlotRules(
    frozenset({"bondtest"}),
    {"bondtest": ("Gamma_bond", "chi_w")},
    {"bondtest": ("chibar_w",)},
    {("bondtest", "flex_so_bond"): {"Gamma_bond": (3,), "chi_w": (4,), "chibar_w": (4,)}})

#: the provenance members of a fingerprint_version 2 bond-resolved archive
PROVENANCE = {
    "fingerprint_version", "producer_stage", "producer_label", "producer_version",
    "fp_transfer", "fp_extern", "fp_meta", "fp_band", "fp_green",
    "prov_transfer_keys", "prov_transfer_vals", "prov_extern_present",
    "prov_extern_keys", "prov_extern_vals",
    "prov_coeff_extern", "prov_enable_spin_orbital", "prov_index_convention",
    "prov_norb", "prov_norb_phys", "prov_cell_shape", "prov_sub_shape",
    "state_beta", "state_nmat", "state_mu", "state_density", "state_calc_mu",
    "state_target", "state_mu_tol",
    "fp_interaction", "fp_topology", "fp_resolved_topology", "fp_capability",
    "prov_rec_types", "prov_rec_modes", "prov_rec_cells", "prov_rec_v", "prov_topology",
    "delta_r", "reverse", "slot_convention", "max_shells",
}


def _interaction():
    r = np.zeros(1, dtype=RECORD_DTYPE)
    r[0]["type"] = "CoulombIntra"
    r[0]["modes"] = (0, 1, 1, 0)
    r[0]["v"] = 1.0
    return fp.interaction_inputs(types.SimpleNamespace(
        records=r, declared_topology=np.zeros((1, 3), np.int64)))


INTER = _interaction()


def _bond(B):
    """``B`` (odd) channels along x: 0, +1, -1, +2, -2, ... and their reverse."""
    rows = [0] + [s * k for k in range(1, B // 2 + 1) for s in (1, -1)]
    dr = np.zeros((B, 3), np.int64)
    dr[:, 0] = rows
    rev = np.array([0] + [i + 1 if i % 2 else i - 1 for i in range(1, B)], np.int64)
    return fp.BondInputs(dr, rev, SLOT, 1)


def _sector_map(B=3, norb_phys=1):
    return so_bond.sz_sector_map(types.SimpleNamespace(reverse=_bond(B).reverse),
                                 2 * norb_phys)


def _block_diagonal(sm, lead, seed):
    """A random code-slot array ``lead + (ND, ND)`` without off-sector entries."""
    rng = np.random.default_rng(seed)
    M = np.zeros(tuple(lead) + (sm.ND, sm.ND), np.complex128)
    for s in range(3):
        idx = sm.index[sm.offsets[s]:sm.offsets[s + 1]]
        n = len(idx)
        M[..., idx[:, None], idx] = (rng.normal(size=tuple(lead) + (n, n))
                                     + 1j * rng.normal(size=tuple(lead) + (n, n)))
    return M


_LEADS = {"Gamma_bond": (2,), "chi_w": (2, 2), "chibar_w": (2, 2)}


def _archive(path_or_none, *, stage, kind, layout, B=3, norb_phys=1, nvol=4, byteorder="=",
             drop=(), extra={}):
    """A synthetic spin-orbital archive of ``kind`` written by ``stage``.

    ``kind`` ``"gap"`` carries ``gap``; ``"chiq"`` carries ``chiq``;
    ``"bondtest"`` carries ``Gamma_bond`` and ``chi_w`` in ``layout`` (dense:
    ``<p>``; packed: ``<p>@0``, ``<p>@+1``, ``<p>@-1`` and the sector-map
    members) from the same block-diagonal data. ``layout`` ``None`` writes
    no layout member. The payload blocks are stored in ``byteorder``. The
    members ``drop`` are removed and ``extra`` added (or replaced). Writes
    ``np.savez`` to ``path_or_none`` when given; returns the member dict."""
    nd = 2 * norb_phys
    bond = _bond(B) if fp.STAGE_CAPABILITY[stage] == "so_bond" else None
    meta = fp.BandMeta(0.0, True, "spin_block", nd, norb_phys, (nvol, 1, 1), (nvol, 1, 1))
    state = fp.so_state(beta=2.0, nmat=8, mu=-0.3, density=0.8, calc_mu=True, target=0.8)
    d = dict(fp.so_provenance_members(
        producer_stage=stage, producer_label=LABEL,
        transfer={((0, 0, 0), (0, 0)): -1.0 + 0j}, extern=None, meta=meta, fp_green=GREEN_FP,
        state=state, interaction=INTER, bond=bond))
    d["spin_orbital"] = np.bool_(True)
    d["index_convention"] = np.str_("spin_block")
    payload = {}
    if kind == "gap":
        payload["gap"] = np.zeros((nd, nd, nvol, 1, 1), np.complex128)
    elif kind == "chiq":
        payload["chiq"] = np.zeros((1, nvol, nd, nd, nd, nd), np.complex128)
    elif kind == "bondtest":
        sm = so_bond.sz_sector_map(types.SimpleNamespace(reverse=bond.reverse), nd)
        for seed, p in enumerate(("Gamma_bond", "chi_w")):
            M = _block_diagonal(sm, _LEADS[p], seed)
            if layout == "sz_packed":
                payload.update(so_bond.packed_members(
                    p, so_bond.PackedArrays.from_dense(M, sm, source=p)))
            else:
                payload[p] = M
        if layout == "sz_packed":
            d.update(so_bond.sector_map_members(sm))
    else:
        raise ValueError(kind)
    dt = np.dtype(np.complex128).newbyteorder(byteorder)
    d.update({k: v.astype(dt) for k, v in payload.items()})
    if layout is not None:
        d["block_layout"] = np.array(layout, dtype="<U16")
    for k in drop:
        del d[k]
    d.update(extra)
    if path_or_none is not None:
        np.savez(path_or_none, **d)
    return d


def _malformed_text(check, stored, expected, path=PATH):
    return ("'{}' is malformed or truncated: {} ({} vs {}). Regenerate it with {}."
            .format(path, check, stored, expected, LABEL))


class _TmpDir(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self._n = 0

    def path(self):
        self._n += 1
        return os.path.join(self._tmp.name, "a{}.npz".format(self._n))


class TestGapLayoutMembers(unittest.TestCase):
    """The layout members of a bond gap and of an on-site archive, checked
    through ``validate_structure``."""

    def _gap(self, layout="dense", **kw):
        return _archive(None, stage="sc_so_bond_static", kind="gap", layout=layout, **kw)

    def _refused(self, d, check, stored, expected, kind="gap"):
        with self.assertRaises(fp.SOArchiveError) as cm:
            fp.validate_structure(d, "gap.npz", kind=kind)
        self.assertEqual(str(cm.exception),
                         _malformed_text(check, stored, expected, path="gap.npz"))

    def test_accepted(self):
        for layout in ("dense", "sz_packed"):
            with self.subTest(layout=layout):
                self.assertEqual(fp.validate_structure(self._gap(layout), "gap.npz", kind="gap"),
                                 ("gap", "sc_so_bond_static", 2))
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "gap.npz")
            _archive(p, stage="sc_so_bond_static", kind="gap", layout="sz_packed")
            with np.load(p) as z:
                self.assertEqual(fp.validate_structure(z, "gap.npz", kind="gap"),
                                 ("gap", "sc_so_bond_static", 2))

    def test_block_layout_refused(self):
        self._refused(self._gap(drop=("block_layout",)), "member block_layout is missing",
                      "absent", "present")
        for value in ("Dense", "sz-packed"):
            with self.subTest(value=value):
                self._refused(self._gap(value), "block_layout is not a known layout",
                              repr(value), "'dense' or 'sz_packed'")
        self._refused(self._gap(extra={"block_layout": np.array(["dense"])}),
                      "block_layout is not a known layout", "['dense']",
                      "'dense' or 'sz_packed'")

    def test_reserved_prefix_refused(self):
        for name in ("sz_sector_order", "sz_sector_indx", "sz_sector_"):
            for layout in ("dense", "sz_packed"):
                with self.subTest(name=name, layout=layout):
                    self._refused(self._gap(layout, extra={name: np.arange(3)}),
                                  "member {} is not written by a gap archive".format(name),
                                  "present", "absent")

    def test_unknown_member_accepted(self):
        d = self._gap(extra={"notes_extra": np.str_("a note")})
        self.assertEqual(fp.validate_structure(d, "gap.npz", kind="gap"),
                         ("gap", "sc_so_bond_static", 2))

    def test_onsite_stage_refuses_layout_members(self):
        for name, value in (("block_layout", np.array("dense", dtype="<U16")),
                            ("sz_sector_indx", np.arange(3))):
            with self.subTest(name=name):
                d = _archive(None, stage="rpa_so", kind="chiq", layout=None,
                             extra={name: value})
                self._refused(d, "member {} is not written by a so_onsite stage".format(name),
                              "present", "absent", kind="chiq")
        d = _archive(None, stage="rpa_so", kind="chiq", layout=None)
        self.assertEqual(fp.validate_structure(d, "c.npz", kind="chiq"), ("chiq", "rpa_so", 2))

    def test_tables_unchanged(self):
        self.assertEqual(fp.DEFERRED_STAGES, {
            "sc_so_bond_dynamic": "the dynamic spin-orbital bond path",
            "flex_so_bond": "the spin-orbital bond FLEX"})
        self.assertEqual(fp.KIND_STAGES, {
            "chi0q": ("rpa_so", "flex_so_onsite"),
            "chiq": ("rpa_so", "flex_so_onsite"),
            "green": ("flex_so_onsite",),
            "sigma": ("flex_so_onsite",),
            "gap": ("sc_so_static", "sc_so_dynamic", "sc_so_bond_static")})
        self.assertEqual(fp.PRODUCTION_SLOT_RULES, SlotRules(frozenset(), {}, {}, {}))
        d = _archive(None, stage="flex_so_bond", kind="bondtest", layout="sz_packed")
        with self.assertRaises(fp.SOArchiveError) as cm:
            fp.validate_structure(d, "b.npz", kind="chiq")
        self.assertEqual(str(cm.exception),
                         "'b.npz' was written by producer_stage 'flex_so_bond', which this "
                         "version does not read yet (the spin-orbital bond FLEX arrives in a "
                         "later version).")


def _check(data, *, zf, rules=TEST_RULES, B=3, nd=2, stage="flex_so_bond", kind="bondtest"):
    """Both slot checks of ``data``, as a consumer of the slot kind runs them."""
    names = fp._members(data)
    malformed = fp._malformed_raiser(PATH, LABEL)
    layout = fp.check_layout_members(data, names, kind=kind, stage=stage, B=B, nd=nd,
                                     malformed=malformed, path=PATH, label=LABEL, rules=rules)
    fp.check_slot_payloads(data, names, zf, kind=kind, stage=stage, layout=layout, B=B, nd=nd,
                           sector_dims=_sector_map(B, nd // 2).sizes, malformed=malformed,
                           rules=rules)
    return layout


class TestSlotMachinery(_TmpDir):
    """The slot-archive checks on synthetic archives, each as a zip and as a
    plain mapping."""

    SM = _sector_map()

    def _each(self, layout, **kw):
        """``(form, data, zf)`` for the zip file and the plain mapping."""
        p = self.path()
        d = _archive(p, stage="flex_so_bond", kind="bondtest", layout=layout, **kw)
        z = np.load(p)
        self.addCleanup(z.close)
        return (("zip", z, z.zip), ("mapping", d, None))

    def _refused(self, layout, check, stored, expected, forms=("zip", "mapping"), **kw):
        for form, data, zf in self._each(layout, **kw):
            if form not in forms:
                continue
            with self.subTest(form=form, check=check), \
                    self.assertRaises(fp.SOArchiveError) as cm:
                _check(data, zf=zf)
            self.assertEqual(str(cm.exception), _malformed_text(check, stored, expected))

    def test_accepted(self):
        for layout in ("dense", "sz_packed"):
            for form, data, zf in self._each(layout):
                with self.subTest(layout=layout, form=form):
                    self.assertEqual(_check(data, zf=zf), layout)

    def test_payload_names(self):
        g = np.zeros((2, 3, 3), np.complex128)
        c = np.zeros((2, 2, 3, 3), np.complex128)
        self._refused("dense", "payload member Gamma_bond does not follow block_layout dense",
                      "['Gamma_bond', 'Gamma_bond@0']", "['Gamma_bond']",
                      extra={"Gamma_bond@0": g})
        self._refused("sz_packed", "payload member chi_w does not follow block_layout sz_packed",
                      "['chi_w@+1', 'chi_w@0']", "['chi_w@0', 'chi_w@+1', 'chi_w@-1']",
                      drop=("chi_w@-1",))
        self._refused("sz_packed", "payload member chi_w does not follow block_layout sz_packed",
                      "['chi_w', 'chi_w@+1', 'chi_w@-1', 'chi_w@0']",
                      "['chi_w@0', 'chi_w@+1', 'chi_w@-1']",
                      extra={"chi_w": np.zeros((2, 2, 12, 12), np.complex128)})
        self._refused("sz_packed",
                      "payload member chibar_w does not follow block_layout sz_packed",
                      "['chibar_w@0']", "['chibar_w@0', 'chibar_w@+1', 'chibar_w@-1']",
                      extra={"chibar_w@0": c})

    def test_optional_payload_accepted(self):
        n = self.SM.sizes
        extra = {"chibar_w@" + s: np.zeros((2, 2, k, k), np.complex128)
                 for s, k in zip(so_bond.SECTOR_ORDER, n)}
        for form, data, zf in self._each("sz_packed", extra=extra):
            with self.subTest(form=form):
                self.assertEqual(_check(data, zf=zf), "sz_packed")

    def test_reserved_names(self):
        for name in ("sz_sector_index", "sz_sector_indx"):
            self._refused("dense",
                          "member {} is not written by a dense slot archive".format(name),
                          "present", "absent", extra={name: np.arange(12)})
        self._refused("sz_packed", "member sz_sector_indx is not a sector-map member (the "
                      "prefix sz_sector_ is reserved)", "present", "absent",
                      extra={"sz_sector_indx": np.arange(12)})
        for form, data, zf in self._each("sz_packed", extra={"extra_note": np.str_("x")}):
            with self.subTest(form=form):
                self.assertEqual(_check(data, zf=zf), "sz_packed")

    def test_permuted_sector_map(self):
        idx = np.array(self.SM.index, np.int64)
        idx[[1, 4]] = idx[[4, 1]]
        text = ("'{}' is malformed or truncated: sz_sector_index differs from the map of "
                "(B = 3, nd = 2) at position 1 ({} vs {}). Regenerate it with {}.".format(
                    PATH, idx[1], self.SM.index[1], LABEL))
        for form, data, zf in self._each("sz_packed", extra={"sz_sector_index": idx}):
            with self.subTest(form=form), self.assertRaises(fp.SOArchiveError) as cm:
                _check(data, zf=zf)
            self.assertEqual(str(cm.exception), text)

    def test_block_headers(self):
        n = self.SM.sizes[1]
        good = np.zeros((2, n, n), np.complex128)
        self._refused("sz_packed", "payload member Gamma_bond@+1 has the wrong slot shape",
                      str((n, n + 1)), str((n, n)),
                      extra={"Gamma_bond@+1": np.zeros((2, n, n + 1), np.complex128)})
        self._refused("dense", "payload member Gamma_bond has the wrong slot shape",
                      "(12, 13)", "(12, 12)",
                      extra={"Gamma_bond": np.zeros((2, 12, 13), np.complex128)})
        self._refused("sz_packed", "payload member Gamma_bond@+1 has the wrong rank", "2", "3",
                      extra={"Gamma_bond@+1": good[0]})
        self._refused("sz_packed", "payload member Gamma_bond@+1 has the wrong dtype",
                      "complex64", "complex128",
                      extra={"Gamma_bond@+1": good.astype(np.complex64)})
        self._refused("sz_packed", "payload member Gamma_bond@+1 is in Fortran order",
                      "Fortran", "C",
                      extra={"Gamma_bond@+1": np.asfortranarray(np.zeros((2, n, n), complex))})

    def test_truncated(self):
        p = self.path()
        _archive(p, stage="flex_so_bond", kind="bondtest", layout="sz_packed")
        q = self.path()
        with zipfile.ZipFile(p) as src, zipfile.ZipFile(q, "w") as dst:
            for info in src.infolist():
                raw = src.read(info)
                if info.filename == "chi_w@-1.npy":
                    raw, size = raw[:-16], len(raw)
                dst.writestr(info.filename, raw)
        with np.load(q) as z, self.assertRaises(fp.SOArchiveError) as cm:
            _check(z, zf=z.zip)
        self.assertEqual(str(cm.exception), _malformed_text(
            "payload member chi_w@-1 has the wrong size", str(size - 16), str(size)))

    def test_mapping_never_skips_shapes(self):
        """A plain mapping whose block has the wrong slot shape passes every
        other check, and is refused by its shape."""
        n = self.SM.sizes[0]
        d = _archive(None, stage="flex_so_bond", kind="bondtest", layout="sz_packed",
                     extra={"chi_w@0": np.zeros((2, 2, n + 1, n + 1), np.complex128)})
        names = fp._members(d)
        malformed = fp._malformed_raiser(PATH, LABEL)
        self.assertEqual(fp.check_layout_members(
            d, names, kind="bondtest", stage="flex_so_bond", B=3, nd=2, malformed=malformed,
            path=PATH, label=LABEL, rules=TEST_RULES), "sz_packed")
        with self.assertRaises(fp.SOArchiveError) as cm:
            fp.check_slot_payloads(d, names, None, kind="bondtest", stage="flex_so_bond",
                                   layout="sz_packed", B=3, nd=2, sector_dims=self.SM.sizes,
                                   malformed=malformed, rules=TEST_RULES)
        self.assertEqual(str(cm.exception), _malformed_text(
            "payload member chi_w@0 has the wrong slot shape", str((n + 1, n + 1)),
            str((n, n))))

    def test_byte_order_check(self):
        for layout in ("dense", "sz_packed"):
            for form, data, zf in self._each(layout, byteorder=">"):
                with self.subTest(layout=layout, form=form):
                    self.assertEqual(data["chi_w" if layout == "dense" else "chi_w@0"]
                                     .dtype.byteorder, ">")
                    self.assertEqual(_check(data, zf=zf), layout)


class _Spy(dict):
    """A member mapping recording every ``__getitem__``."""
    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.reads = []

    def __getitem__(self, k):
        self.reads.append(k)
        return super().__getitem__(k)


class TestSlotRules(unittest.TestCase):
    def test_incomplete_rules_fail_closed(self):
        r = TEST_RULES
        ranks = dict(r.ranks[("bondtest", "flex_so_bond")])
        del ranks["chibar_w"]
        cases = (
            (r._replace(payloads={}), "SlotRules: slot kind bondtest has no payloads entry"),
            (r._replace(ranks={}),
             "SlotRules: slot kind bondtest has no ranks for stage flex_so_bond entry"),
            (r._replace(ranks={("bondtest", "flex_so_bond"): ranks}),
             "SlotRules: slot kind bondtest has no rank for payload chibar_w on stage "
             "flex_so_bond entry"),
        )
        d = _archive(None, stage="flex_so_bond", kind="bondtest", layout="sz_packed")
        for rules, text in cases:
            spy = _Spy(d)
            names = set(d)
            malformed = fp._malformed_raiser(PATH, LABEL)
            with self.subTest(text=text):
                with self.assertRaises(ValueError) as cm:
                    fp.check_layout_members(spy, names, kind="bondtest", stage="flex_so_bond",
                                            B=3, nd=2, malformed=malformed, path=PATH,
                                            label=LABEL, rules=rules)
                self.assertIs(type(cm.exception), ValueError)
                self.assertEqual(str(cm.exception), text)
                with self.assertRaises(ValueError) as cm:
                    fp.check_slot_payloads(spy, names, None, kind="bondtest",
                                           stage="flex_so_bond", layout="sz_packed", B=3, nd=2,
                                           sector_dims=_sector_map().sizes,
                                           malformed=malformed, rules=rules)
                self.assertIs(type(cm.exception), ValueError)
                self.assertEqual(str(cm.exception), text)
                self.assertEqual(spy.reads, [])

    def test_production_rules_admit_no_slot_kind(self):
        self.assertEqual(fp.PRODUCTION_SLOT_RULES.slot_kinds, frozenset())
        fp._check_rules(fp.PRODUCTION_SLOT_RULES, "gap", "sc_so_bond_static")


class TestLayoutTables(unittest.TestCase):
    def test_no_new_digest(self):
        dense = _archive(None, stage="flex_so_bond", kind="bondtest", layout="dense")
        packed = _archive(None, stage="flex_so_bond", kind="bondtest", layout="sz_packed")
        fps = {k for k in dense if k.startswith("fp_")}
        self.assertEqual(fps, {k for k in packed if k.startswith("fp_")})
        for k in fps:
            self.assertEqual(str(dense[k]), str(packed[k]))
        for d in (dense, packed):
            payload = {k for k in d if k.split("@")[0] in ("Gamma_bond", "chi_w")}
            self.assertEqual(set(d) - payload - set(fp.LAYOUT_MEMBERS) - MARKERS, PROVENANCE)

    def test_byte_order_read(self):
        """Big-endian blocks are read as native complex128, equal to the
        native archive's, in both layouts and both directions."""
        sm = _sector_map()
        for layout in ("dense", "sz_packed"):
            native = _archive(None, stage="flex_so_bond", kind="bondtest", layout=layout)
            big = _archive(None, stage="flex_so_bond", kind="bondtest", layout=layout,
                           byteorder=">")
            for want in ("sz_packed", "dense"):
                with self.subTest(layout=layout, want=want):
                    ref = so_bond.read_slot_member(native, "chi_w", want=want, sector_map=sm,
                                                   file=PATH)
                    got = so_bond.read_slot_member(big, "chi_w", want=want, sector_map=sm,
                                                   file=PATH)
                    pairs = (zip(got.blocks, ref.blocks) if want == "sz_packed"
                             else ((got, ref),))
                    for g, r in pairs:
                        self.assertEqual(g.dtype, np.dtype(np.complex128))
                        self.assertTrue(g.dtype.isnative)
                        np.testing.assert_array_equal(g, r)

    def test_layout_members_literal(self):
        self.assertEqual(fp.LAYOUT_MEMBERS, ("block_layout",) + so_bond._MAP_MEMBERS)
        self.assertEqual(fp.RESERVED_PREFIX, "sz_sector_")
        for name in so_bond._MAP_MEMBERS:
            self.assertTrue(name.startswith(fp.RESERVED_PREFIX))


if __name__ == "__main__":
    unittest.main()
