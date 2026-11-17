"""The layout members of bond-resolved spin-orbital archives.

* A bond-resolved archive records its slot layout in ``block_layout``
  (``"dense"`` or ``"sz_packed"``); a ``sc_so_bond_static`` gap carries it
  (the solve layout; the gap itself is always stored full) and no member
  with the reserved prefix ``sz_sector_``. An on-site stage carries none of
  these members. Other unknown members are accepted (forward compatible).
* The slot-archive machinery (``check_layout_members`` with a slot kind and
  ``check_slot_payloads``) reads the kind ``"bond"`` of the ``flex_so_bond``
  stage (the production table); every row runs both on that table and on a
  test-local table naming a synthetic slot kind ``bondtest``, on synthetic
  archives written by the builder below, both as zip files and as plain
  mappings (whose arrays are checked, never skipped).
* A slot block is stored inline (``<block>``) or, with a sidecar directory,
  as a sidecar (``<block>_file`` naming a ``.npy`` file in that directory
  and ``<block>_nbytes`` its size); the npy header of a sidecar is read,
  never its data.
* The other kinds of a bond-resolved stage (``chi0q``, ``chiq``, ``green``,
  ``sigma``) carry no member with the reserved prefix.
"""
import io
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

    ``kind`` ``"gap"`` carries ``gap``; ``"chiq"``, ``"chi0q"``, ``"green"``
    and ``"sigma"`` their payloads; ``"bondtest"`` and ``"bond"`` carry
    ``Gamma_bond`` and ``chi_w`` in ``layout`` (dense:
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
    elif kind in ("chiq", "chi0q"):
        payload[kind] = np.zeros((1, nvol, nd, nd, nd, nd), np.complex128)
    elif kind == "green":
        payload["green"] = np.zeros((1, 2, nvol, nd, nd), np.complex128)
        d["green_layout"] = np.str_(fp.GREEN_LAYOUT)
    elif kind == "sigma":
        payload["sigma_static"] = np.zeros((1, 1, nvol, nd, nd), np.complex128)
    elif kind in ("bondtest", "bond"):
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


def _sidecar(d, directory, blocks, nonce):
    """Move the ``blocks`` of the member dict ``d`` to sidecar files
    ``<block>.<nonce>.npy`` in ``directory`` (``<block>_file`` the base name,
    ``<block>_nbytes`` the file size); returns ``{block: sidecar path}``."""
    out = {}
    for b in blocks:
        name = "{}.{}.npy".format(b, nonce)
        full = os.path.join(directory, name)
        np.save(full, d.pop(b))
        d[b + "_file"] = np.str_(name)
        d[b + "_nbytes"] = np.int64(os.path.getsize(full))
        out[b] = full
    return out


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
    plain mapping, on the test-local rules (``bondtest``)."""

    SM = _sector_map()
    RULES = TEST_RULES
    KIND = "bondtest"

    def _check(self, data, zf):
        return _check(data, zf=zf, rules=self.RULES, kind=self.KIND)

    def _each(self, layout, **kw):
        """``(form, data, zf)`` for the zip file and the plain mapping."""
        p = self.path()
        d = _archive(p, stage="flex_so_bond", kind=self.KIND, layout=layout, **kw)
        z = np.load(p)
        self.addCleanup(z.close)
        return (("zip", z, z.zip), ("mapping", d, None))

    def _refused(self, layout, check, stored, expected, forms=("zip", "mapping"), **kw):
        for form, data, zf in self._each(layout, **kw):
            if form not in forms:
                continue
            with self.subTest(form=form, check=check), \
                    self.assertRaises(fp.SOArchiveError) as cm:
                self._check(data, zf)
            self.assertEqual(str(cm.exception), _malformed_text(check, stored, expected))

    def test_accepted(self):
        for layout in ("dense", "sz_packed"):
            for form, data, zf in self._each(layout):
                with self.subTest(layout=layout, form=form):
                    self.assertEqual(self._check(data, zf), layout)

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
                self.assertEqual(self._check(data, zf), "sz_packed")

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
                self.assertEqual(self._check(data, zf), "sz_packed")

    def test_permuted_sector_map(self):
        idx = np.array(self.SM.index, np.int64)
        idx[[1, 4]] = idx[[4, 1]]
        text = ("'{}' is malformed or truncated: sz_sector_index differs from the map of "
                "(B = 3, nd = 2) at position 1 ({} vs {}). Regenerate it with {}.".format(
                    PATH, idx[1], self.SM.index[1], LABEL))
        for form, data, zf in self._each("sz_packed", extra={"sz_sector_index": idx}):
            with self.subTest(form=form), self.assertRaises(fp.SOArchiveError) as cm:
                self._check(data, zf)
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
        _archive(p, stage="flex_so_bond", kind=self.KIND, layout="sz_packed")
        q = self.path()
        with zipfile.ZipFile(p) as src, zipfile.ZipFile(q, "w") as dst:
            for info in src.infolist():
                raw = src.read(info)
                if info.filename == "chi_w@-1.npy":
                    raw, size = raw[:-16], len(raw)
                dst.writestr(info.filename, raw)
        with np.load(q) as z, self.assertRaises(fp.SOArchiveError) as cm:
            self._check(z, z.zip)
        self.assertEqual(str(cm.exception), _malformed_text(
            "payload member chi_w@-1 has the wrong size", str(size - 16), str(size)))

    def test_mapping_never_skips_shapes(self):
        """A plain mapping whose block has the wrong slot shape passes every
        other check, and is refused by its shape."""
        n = self.SM.sizes[0]
        d = _archive(None, stage="flex_so_bond", kind=self.KIND, layout="sz_packed",
                     extra={"chi_w@0": np.zeros((2, 2, n + 1, n + 1), np.complex128)})
        names = fp._members(d)
        malformed = fp._malformed_raiser(PATH, LABEL)
        self.assertEqual(fp.check_layout_members(
            d, names, kind=self.KIND, stage="flex_so_bond", B=3, nd=2, malformed=malformed,
            path=PATH, label=LABEL, rules=self.RULES), "sz_packed")
        with self.assertRaises(fp.SOArchiveError) as cm:
            fp.check_slot_payloads(d, names, None, kind=self.KIND, stage="flex_so_bond",
                                   layout="sz_packed", B=3, nd=2, sector_dims=self.SM.sizes,
                                   malformed=malformed, rules=self.RULES)
        self.assertEqual(str(cm.exception), _malformed_text(
            "payload member chi_w@0 has the wrong slot shape", str((n + 1, n + 1)),
            str((n, n))))

    def test_byte_order_check(self):
        for layout in ("dense", "sz_packed"):
            for form, data, zf in self._each(layout, byteorder=">"):
                with self.subTest(layout=layout, form=form):
                    self.assertEqual(data["chi_w" if layout == "dense" else "chi_w@0"]
                                     .dtype.byteorder, ">")
                    self.assertEqual(self._check(data, zf), layout)


class TestSlotMachineryProduction(TestSlotMachinery):
    """The same rows on the production rules (kind ``bond``)."""
    RULES = fp.PRODUCTION_SLOT_RULES
    KIND = "bond"


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

    def test_production_rules_complete(self):
        self.assertEqual(fp.PRODUCTION_SLOT_RULES.slot_kinds, frozenset({"bond"}))
        fp._check_rules(fp.PRODUCTION_SLOT_RULES, "bond", "flex_so_bond")
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



class TestProductionRules(unittest.TestCase):
    """The tables of this version: the bond FLEX stage is read for its five
    archive kinds; only the dynamic bond stage is deferred."""

    def test_tables(self):
        self.assertEqual(fp.ARTIFACT_KINDS, ("chi0q", "chiq", "green", "sigma", "gap", "bond"))
        self.assertEqual(fp.KIND_STAGES, {
            "chi0q": ("rpa_so", "flex_so_onsite", "flex_so_bond"),
            "chiq": ("rpa_so", "flex_so_onsite", "flex_so_bond"),
            "green": ("flex_so_onsite", "flex_so_bond"),
            "sigma": ("flex_so_onsite", "flex_so_bond"),
            "gap": ("sc_so_static", "sc_so_dynamic", "sc_so_bond_static"),
            "bond": ("flex_so_bond",)})
        self.assertEqual(fp.DEFERRED_STAGES,
                         {"sc_so_bond_dynamic": "the dynamic spin-orbital bond path"})
        self.assertEqual(fp.PAYLOAD_MEMBERS, {
            "chi0q": ("chi0q",), "chiq": ("chiq",), "green": ("green", "green_layout"),
            "sigma": ("sigma_static",), "gap": ("gap",), "bond": ()})
        self.assertEqual(fp.PAYLOAD_RANKS, {
            ("chi0q", "rpa_so"): (4, 5, 6, 7), ("chi0q", "flex_so_onsite"): (4, 6),
            ("chi0q", "flex_so_bond"): (6,),
            ("chiq", "rpa_so"): (4, 6), ("chiq", "flex_so_onsite"): (6,),
            ("chiq", "flex_so_bond"): (6,),
            ("green", "flex_so_onsite"): (5,), ("green", "flex_so_bond"): (5,),
            ("sigma", "flex_so_onsite"): (5,), ("sigma", "flex_so_bond"): (5,),
            ("gap", "sc_so_static"): (5,), ("gap", "sc_so_dynamic"): (6,),
            ("gap", "sc_so_bond_static"): (5,)})
        self.assertEqual(fp.PRODUCTION_SLOT_RULES, SlotRules(
            frozenset({"bond"}), {"bond": ("Gamma_bond", "chi_w")}, {"bond": ("chibar_w",)},
            {("bond", "flex_so_bond"): {"Gamma_bond": (3,), "chi_w": (4,),
                                        "chibar_w": (4,)}}))

    def test_required_members(self):
        for kind in ("chi0q", "chiq", "green", "sigma", "bond"):
            with self.subTest(kind=kind):
                self.assertEqual(fp.required_members(kind, "flex_so_bond"),
                                 frozenset(PROVENANCE))
        with self.assertRaises(ValueError) as cm:
            fp.required_members("bond", "rpa_so")
        self.assertEqual(str(cm.exception),
                         "artifact kind 'bond' has no producer_stage 'rpa_so' (its stages: "
                         "flex_so_bond)")

    def test_bond_archive_validates(self):
        with tempfile.TemporaryDirectory() as tmp:
            for layout in ("dense", "sz_packed"):
                p = os.path.join(tmp, layout + ".npz")
                d = _archive(p, stage="flex_so_bond", kind="bond", layout=layout)
                with np.load(p) as z:
                    for form, data in (("zip", z), ("mapping", d)):
                        with self.subTest(layout=layout, form=form):
                            self.assertEqual(fp.validate_structure(data, PATH, kind="bond"),
                                             ("bond", "flex_so_bond", 2))

    def test_signatures(self):
        import inspect
        for f, names in ((fp.validate_structure, ("kind", "archive", "expected_producer",
                                                  "sidecar_dir", "sources_out")),
                         (fp.check_slot_payloads, ("kind", "stage", "layout", "B", "nd",
                                                   "sector_dims", "malformed", "rules",
                                                   "sidecar_dir", "sources_out"))):
            params = inspect.signature(f).parameters
            with self.subTest(f=f.__name__):
                for n in names:
                    self.assertEqual(params[n].kind, inspect.Parameter.KEYWORD_ONLY)
                for n in ("sidecar_dir", "sources_out"):
                    self.assertIsNone(params[n].default)


class TestReservedNonSlot(unittest.TestCase):
    """The non-slot kinds of the bond FLEX stage carry block_layout and no
    member with the reserved prefix; a gap and an on-site stage keep their
    own texts."""

    def _refused(self, d, kind, check, stored, expected):
        with self.assertRaises(fp.SOArchiveError) as cm:
            fp.validate_structure(d, PATH, kind=kind)
        self.assertEqual(str(cm.exception), _malformed_text(check, stored, expected))

    def test_each_kind(self):
        with tempfile.TemporaryDirectory() as tmp:
            for kind in ("chi0q", "chiq", "green", "sigma"):
                for layout in ("dense", "sz_packed"):
                    p = os.path.join(tmp, "{}_{}.npz".format(kind, layout))
                    d = _archive(p, stage="flex_so_bond", kind=kind, layout=layout)
                    with self.subTest(kind=kind, layout=layout):
                        with np.load(p) as z:
                            self.assertEqual(fp.validate_structure(z, PATH, kind=kind),
                                             (kind, "flex_so_bond", 2))
                        self.assertEqual(fp.validate_structure(d, PATH, kind=kind),
                                         (kind, "flex_so_bond", 2))
                        for name in ("sz_sector_order", "sz_sector_indx", "sz_sector_"):
                            bad = dict(d, **{name: np.arange(3)})
                            self._refused(bad, kind, "member {} is not written by a {} archive "
                                          "of a bond-resolved stage".format(name, kind),
                                          "present", "absent")
                        self._refused({k: v for k, v in d.items() if k != "block_layout"},
                                      kind, "member block_layout is missing", "absent",
                                      "present")
                        # unknown members stay accepted
                        self.assertEqual(fp.validate_structure(
                            dict(d, notes_extra=np.str_("x")), PATH, kind=kind),
                            (kind, "flex_so_bond", 2))

    def test_sector_map_on_a_chiq(self):
        """The sector-map members of a packed bond archive are refused on a
        chiq (the first reserved name in sorted order)."""
        d = _archive(None, stage="flex_so_bond", kind="bond", layout="sz_packed")
        d["chiq"] = np.zeros((1, 4, 2, 2, 2, 2), np.complex128)
        self._refused(d, "chiq", "member sz_sector_index is not written by a chiq archive of a "
                      "bond-resolved stage", "present", "absent")

    def test_gap_and_onsite_texts_unchanged(self):
        d = _archive(None, stage="sc_so_bond_static", kind="gap", layout="dense",
                     extra={"sz_sector_order": np.arange(3)})
        with self.assertRaises(fp.SOArchiveError) as cm:
            fp.validate_structure(d, PATH, kind="gap")
        self.assertEqual(str(cm.exception), _malformed_text(
            "member sz_sector_order is not written by a gap archive", "present", "absent"))
        d = _archive(None, stage="flex_so_onsite", kind="chiq", layout=None,
                     extra={"sz_sector_order": np.arange(3)})
        with self.assertRaises(fp.SOArchiveError) as cm:
            fp.validate_structure(d, PATH, kind="chiq")
        self.assertEqual(str(cm.exception), _malformed_text(
            "member sz_sector_order is not written by a so_onsite stage", "present", "absent"))


NONCE = "0123456789abcdef" * 2


class TestSidecar(unittest.TestCase):
    """Sidecar blocks of a bond archive: accepted in both layouts, as a zip
    and as a plain mapping, with the sources map filled; each malformed
    variant is built from an archive this class first accepts."""

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        cls.base = {}
        for layout in ("dense", "sz_packed"):
            sub = os.path.join(cls._tmp.name, layout)
            os.makedirs(sub)
            d = _archive(None, stage="flex_so_bond", kind="bond", layout=layout)
            blocks = sorted(n for n in d if n.split("@")[0] == "chi_w")
            shapes = {n: tuple(d[n].shape) for n in d if n.split("@")[0] in
                      ("Gamma_bond", "chi_w")}
            files = _sidecar(d, sub, blocks, NONCE)
            p = os.path.join(sub, "bond.npz")
            np.savez(p, **d)
            cls.base[layout] = types.SimpleNamespace(dir=sub, path=p, data=d, files=files,
                                                     shapes=shapes, blocks=blocks)
        cls._n = 0

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def _forms(self, d):
        """``(form, data)``: the plain mapping and a zip written from it."""
        type(self)._n += 1
        p = os.path.join(self._tmp.name, "v{}.npz".format(self._n))
        np.savez(p, **d)
        z = np.load(p)
        self.addCleanup(z.close)
        return (("mapping", d), ("zip", z))

    def _file(self, layout, tag, array=None, raw=None):
        """A sidecar file ``<tag>.npy`` in the layout's directory."""
        b = self.base[layout]
        name = tag + ".npy"
        full = os.path.join(b.dir, name)
        if raw is None:
            np.save(full, array)
        else:
            with open(full, "wb") as fh:
                fh.write(raw)
        return name

    def _validate(self, data, sidecar_dir, sources=None):
        return fp.validate_structure(data, PATH, kind="bond", sidecar_dir=sidecar_dir,
                                     sources_out=sources)

    def _refused(self, layout, d, check, stored, expected, sidecar_dir="default"):
        b = self.base[layout]
        sidecar_dir = b.dir if sidecar_dir == "default" else sidecar_dir
        for form, data in self._forms(d):
            with self.subTest(layout=layout, check=check, form=form):
                with self.assertRaises(fp.SOArchiveError) as cm:
                    self._validate(data, sidecar_dir)
                self.assertEqual(str(cm.exception), _malformed_text(check, stored, expected))

    def _variant(self, layout, drop=(), **extra):
        d = {k: v for k, v in self.base[layout].data.items() if k not in drop}
        d.update(extra)
        return d

    def test_accepted(self):
        for layout, b in self.base.items():
            want = {n: (b.files[n] if n in b.files else "inline", shape)
                    for n, shape in b.shapes.items()}
            with np.load(b.path) as z:
                for form, data in (("zip", z), ("mapping", b.data)):
                    with self.subTest(layout=layout, form=form):
                        sources = {}
                        self.assertEqual(self._validate(data, b.dir, sources),
                                         ("bond", "flex_so_bond", 2))
                        self.assertEqual(sources, want)
                        self.assertEqual(self._validate(data, b.dir),
                                         ("bond", "flex_so_bond", 2))

    def test_every_block_sidecar_and_inline_only(self):
        """The generic check admits any block as a sidecar (the schema split
        is the archive reader's), and an all-inline archive with a sidecar
        directory reports every block inline."""
        with tempfile.TemporaryDirectory() as tmp:
            for layout in ("dense", "sz_packed"):
                d = _archive(None, stage="flex_so_bond", kind="bond", layout=layout)
                inline = dict(d)
                blocks = sorted(n for n in d if n.split("@")[0] in ("Gamma_bond", "chi_w"))
                shapes = {n: tuple(d[n].shape) for n in blocks}
                files = _sidecar(d, tmp, blocks, NONCE)
                for data, want in ((d, {n: (files[n], shapes[n]) for n in blocks}),
                                   (inline, {n: ("inline", shapes[n]) for n in blocks})):
                    with self.subTest(layout=layout, sidecars=data is d):
                        sources = {}
                        self.assertEqual(self._validate(data, tmp, sources),
                                         ("bond", "flex_so_bond", 2))
                        self.assertEqual(sources, want)

    def test_member_rows(self):
        b = self.base["dense"]
        name = str(b.data["chi_w_file"])
        n = int(b.data["chi_w_nbytes"])
        inline = np.load(b.files["chi_w"])
        self._refused("dense", self._variant("dense", chi_w=inline),
                      "payload member chi_w is stored both inline and as a sidecar",
                      "both", "one of them")
        self._refused("dense", self._variant("dense", drop=("chi_w_file",), chi_w=inline),
                      "payload member chi_w is stored both inline and as a sidecar",
                      "both", "one of them")
        self._refused("dense", self._variant("dense", drop=("chi_w_nbytes",)),
                      "sidecar member chi_w_nbytes is missing", "absent", "present")
        self._refused("dense", self._variant("dense", drop=("chi_w_file",)),
                      "sidecar member chi_w_file is missing", "absent", "present")
        for bad in (np.array([name]), np.int64(3), np.bytes_(name.encode())):
            a = np.asarray(bad)
            self._refused("dense", self._variant("dense", chi_w_file=bad),
                          "sidecar member chi_w_file is not a 0-d string",
                          "{}/{}d".format(a.dtype, a.ndim), "str/0d")
        for bad in (np.float64(n), np.array([n]), np.bool_(True), np.str_(str(n))):
            a = np.asarray(bad)
            self._refused("dense", self._variant("dense", chi_w_nbytes=bad),
                          "sidecar member chi_w_nbytes is not a 0-d integer",
                          "{}/{}d".format(a.dtype, a.ndim), "int64/0d")
        self._refused("dense", self._variant("dense"),
                      "payload member chi_w is a sidecar but no sidecar directory was given",
                      "sidecar", "inline", sidecar_dir=None)
        for bad in ("sub/" + name, "../" + name, os.path.join(b.dir, name), "a\\b.npy",
                    "", ".", ".."):
            self._refused("dense", self._variant("dense", chi_w_file=np.str_(bad)),
                          "sidecar file name of chi_w is not a plain file name", repr(bad),
                          "a file name without a directory part")
        self._refused("dense", self._variant("dense", chi_w_file=np.str_("gone.npy")),
                      "sidecar member chi_w file 'gone.npy' is missing", "absent", "present")
        os.makedirs(os.path.join(b.dir, "adir.npy"), exist_ok=True)
        self._refused("dense", self._variant("dense", chi_w_file=np.str_("adir.npy")),
                      "sidecar member chi_w file 'adir.npy' is missing", "absent", "present")

    def test_size_rows(self):
        b = self.base["dense"]
        n = int(b.data["chi_w_nbytes"])
        with open(b.files["chi_w"], "rb") as fh:
            raw = fh.read()
        self.assertEqual(len(raw), n)
        self._refused("dense", self._variant("dense", chi_w_nbytes=np.int64(n + 1)),
                      "sidecar member chi_w has the wrong size",
                      "{} bytes, {} recorded".format(n, n + 1), "{} bytes".format(n))
        for tag, cut in (("short", raw[:-16]), ("long", raw + b"\0" * 8)):
            name = self._file("dense", "chi_w." + tag, raw=cut)
            self._refused("dense", self._variant("dense", chi_w_file=np.str_(name),
                                                 chi_w_nbytes=np.int64(len(cut))),
                          "sidecar member chi_w has the wrong size",
                          "{} bytes, {} recorded".format(len(cut), len(cut)),
                          "{} bytes".format(n))

    def test_header_rows(self):
        from hwave.solver import npy_header
        b = self.base["dense"]
        good = np.load(b.files["chi_w"])
        garbage = b"this is not an npy file"
        try:
            npy_header.read_npy_header(io.BytesIO(garbage))
        except ValueError as exc:
            reason = repr(str(exc))
        name = self._file("dense", "chi_w.garbage", raw=garbage)
        # the header is read before the size is compared
        for nbytes in (len(garbage), len(garbage) + 1):
            self._refused("dense", self._variant("dense", chi_w_file=np.str_(name),
                                                 chi_w_nbytes=np.int64(nbytes)),
                          "payload member chi_w has an unreadable npy header", reason,
                          "a valid npy header")
        rows = (
            ("c64", good.astype(np.complex64), "payload member chi_w has the wrong dtype",
             "complex64", "complex128"),
            ("fortran", np.asfortranarray(good), "payload member chi_w is in Fortran order",
             "Fortran", "C"),
            ("rank", good[0], "payload member chi_w has the wrong rank", "3", "4"),
            ("slot", np.zeros((2, 2, 12, 13), np.complex128),
             "payload member chi_w has the wrong slot shape", "(12, 13)", "(12, 12)"),
        )
        for tag, a, check, stored, expected in rows:
            name = self._file("dense", "chi_w." + tag, array=a)
            size = os.path.getsize(os.path.join(b.dir, name))
            self._refused("dense", self._variant("dense", chi_w_file=np.str_(name),
                                                 chi_w_nbytes=np.int64(size)),
                          check, stored, expected)

    def test_packed_rows(self):
        b = self.base["sz_packed"]
        n_s = b.shapes["chi_w@+1"][-1]
        self._refused("sz_packed", self._variant("sz_packed", drop=("chi_w@+1_nbytes",)),
                      "sidecar member chi_w@+1_nbytes is missing", "absent", "present")
        name = self._file("sz_packed", "chi_w@+1.slot",
                          array=np.zeros((2, 2, n_s, n_s + 1), np.complex128))
        size = os.path.getsize(os.path.join(b.dir, name))
        self._refused("sz_packed", self._variant("sz_packed", **{
            "chi_w@+1_file": np.str_(name), "chi_w@+1_nbytes": np.int64(size)}),
            "payload member chi_w@+1 has the wrong slot shape", str((n_s, n_s + 1)),
            str((n_s, n_s)))
        self._refused("sz_packed",
                      self._variant("sz_packed", drop=("chi_w@-1_file", "chi_w@-1_nbytes")),
                      "payload member chi_w does not follow block_layout sz_packed",
                      "['chi_w@+1', 'chi_w@0']", "['chi_w@0', 'chi_w@+1', 'chi_w@-1']")
        # an orphan nbytes alone still names the block
        self._refused("sz_packed", self._variant("sz_packed", drop=("chi_w@-1_file",)),
                      "sidecar member chi_w@-1_file is missing", "absent", "present")

    def test_sidecar_name_on_a_dense_archive(self):
        """A packed block name stored as a sidecar is a stray block of a
        dense archive (the suffixes are stripped before the name check)."""
        b = self.base["dense"]
        self._refused("dense", self._variant("dense", **{
            "chi_w@0_file": b.data["chi_w_file"], "chi_w@0_nbytes": b.data["chi_w_nbytes"]}),
            "payload member chi_w does not follow block_layout dense",
            "['chi_w', 'chi_w@0']", "['chi_w']")

    def test_no_payload_data_read(self):
        """Neither the inline blocks of a mapping-free zip nor the sidecars
        are loaded: the sidecar check opens the files for their header
        only."""
        b = self.base["sz_packed"]
        opened = []
        real_open = open

        def spy_open(f, mode="r", *a, **kw):
            opened.append((f, mode))
            return real_open(f, mode, *a, **kw)

        from unittest import mock
        with np.load(b.path) as z, mock.patch("numpy.load") as load, \
                mock.patch("builtins.open", spy_open):
            self.assertEqual(self._validate(z, b.dir), ("bond", "flex_so_bond", 2))
        load.assert_not_called()
        self.assertEqual(sorted(f for f, _ in opened), sorted(b.files.values()))
        self.assertTrue(all(m == "rb" for _, m in opened))


if __name__ == "__main__":
    unittest.main()
