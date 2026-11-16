"""The schema-4 spin-orbital bond archive: its writer.

* ``write_so_bond_archive`` writes the members of the design's archive
  table (section 6.2) plus the markers of this version, dense or S_z-packed,
  as one npz or as an index with ``.npy`` sidecars named
  ``<block>.<nonce>.npy``; the result passes the fingerprint structure check
  of kind ``"bond"``. Every check here reads the raw members (``np.load``)
  and that structure check.
* No store-sized copy: a sidecar is written straight from the array, an
  inline npz member costs one zip chunk.
* Atomic index: on a failure before the index commit the directory is left
  as it was; a failure after the companion commit removes the previous
  index; a committed archive is never damaged by its own clean-up.
* The sidecars of the archive an index replaces are removed unless another
  index in the directory (or an unreadable or malformed candidate) may
  still name them.
* Preconditions, the ``green_file`` rule and the pole group are refused
  before anything is written.
"""
import contextlib
import logging
import os
import shutil
import tempfile
import tracemalloc
import unittest
from unittest import mock

import numpy as np

from hwave.solver import fingerprint as fp
from hwave.solver import so_bond
from hwave.solver import so_bond_archive as sa
from tests import so_bond_archive_fixtures as fx

LOG = logging.getLogger("hwave.test.so_bond_archive")

_REAL_REPLACE = os.replace
_REAL_REMOVE = os.remove
_REAL_SAVE = np.save

BASE_MEMBERS = {
    "bond_archive_schema", "spin_orbital", "index_convention", "index_order",
    "block_layout", "freq_axis", "bond_archive_layout", "norb_phys", "nd", "nmat", "B",
    "cell_shape", "beta", "types", "green_file", "archive_nonce", "green_layout",
    "momentum_convention", "wavevector_unit", "wavevector_index", "scf_converged",
} | set(fp.SO_PROVENANCE_MEMBERS) | set(fp.SO_BOND_MEMBERS)
MAP_MEMBERS = {"sz_sector_order", "sz_sector_offsets", "sz_sector_index"}
POLE = {"rpa_past_pole": np.bool_(True),
        "rpa_pole_indicator_spinful": np.float64(1.2),
        "rpa_pole_q_spinful": np.array([0.5, 0.0, 0.0])}


def _blocks(p, block_layout):
    return [p] if block_layout == "dense" else [p + "@0", p + "@+1", p + "@-1"]


def expected_members(block_layout, output_layout, chibar=False):
    names = set(BASE_MEMBERS)
    if block_layout == "sz_packed":
        names |= MAP_MEMBERS
    names |= set(_blocks("Gamma_bond", block_layout))
    for p in ["chi_w"] + (["chibar_w"] if chibar else []):
        for b in _blocks(p, block_layout):
            if output_layout == "sidecar":
                names |= {b + "_file", b + "_nbytes"}
            else:
                names.add(b)
    return names


def _listing(d):
    out = {}
    for n in sorted(os.listdir(d)):
        full = os.path.join(d, n)
        if os.path.isfile(full):
            with open(full, "rb") as fh:
                out[n] = fh.read()
        else:
            out[n] = None
    return out


def _sidecar_names(path):
    with np.load(path) as z:
        return sorted(str(z[n][()]) for n in z.files
                      if n.endswith("_file") and n != "green_file")


def _old_info(path, removed=(), kept=(), kept_by=(), nonso=(), replaced=True):
    text = sa.T_OLD_SIDECARS if replaced else sa.T_OLD_SIDECARS_NO_INDEX
    j = lambda v: ", ".join(sorted(v)) if v else "none"  # noqa: E731
    tail = sa.T_NONSO_TAIL.format(names=", ".join(sorted(nonso))) if nonso else ""
    return "INFO:{}:{}".format(LOG.name, text.format(
        path=path, n=len(removed), removed=j(removed), k=len(kept), kept=j(kept),
        kept_by=j(kept_by), nonso=tail))


def _selective(real, target, exc):
    """``real`` except for a call naming the path ``target`` (any positional
    argument: the file of ``os.remove``, the destination of ``os.replace``),
    which raises ``exc``."""
    def fake(*a, **kw):
        if any(isinstance(x, str) and os.path.abspath(x) == os.path.abspath(target)
               for x in a):
            raise exc
        return real(*a, **kw)
    return fake


class _Dir(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = tmp.name

    def sub(self, name):
        d = os.path.join(self.dir, name)
        os.makedirs(d)
        return d

    def write(self, name="lb.npz", *, d=None, inputs=None, output_layout="sidecar",
              block_layout="dense", **kw):
        if inputs is None:
            inputs = fx.archive_inputs(block_layout=block_layout)
        args = dict(inputs.kwargs)
        args.update(kw)
        path = os.path.join(d or self.dir, name)
        out = sa.write_so_bond_archive(path, output_layout=output_layout, log=LOG, **args)
        return path, out

    def assert_valid(self, path, block_layout=None, output_layout=None):
        src = {}
        with np.load(path) as z:
            got = fp.validate_structure(z, path, kind="bond",
                                        sidecar_dir=os.path.dirname(path), sources_out=src)
            self.assertEqual(got, ("bond", "flex_so_bond", 2))
            if block_layout is not None:
                self.assertEqual(str(z["block_layout"][()]), block_layout)
            if output_layout is not None:
                self.assertEqual(str(z["bond_archive_layout"][()]), output_layout)
        for block, (where, _shape) in src.items():
            if block.split("@")[0] == "Gamma_bond":
                self.assertEqual(where, "inline", block)
        return src


# ----------------------------------------------------------------- members
class TestWriter(_Dir):
    def test_members(self):
        for block_layout in ("dense", "sz_packed"):
            for output_layout in ("npz", "sidecar"):
                for chibar in (False, True):
                    with self.subTest(block_layout=block_layout, output_layout=output_layout,
                                      chibar=chibar):
                        self._members(block_layout, output_layout, chibar)

    def _members(self, block_layout, output_layout, chibar):
        d = self.sub("{}-{}-{}".format(block_layout, output_layout, chibar))
        inputs = fx.archive_inputs(block_layout=block_layout, chibar=chibar)
        path, out = self.write(d=d, inputs=inputs, output_layout=output_layout)
        kw = inputs.kwargs
        with np.load(path) as z:
            self.assertEqual(set(z.files), expected_members(block_layout, output_layout,
                                                            chibar))
            nonce = str(z["archive_nonce"][()])
            fixed = {"bond_archive_schema": (np.int64, 4), "spin_orbital": (np.bool_, True),
                     "index_convention": ("<U16", "spin_block"),
                     "index_order": ("<U32", "bond_major_pair_code_slot"),
                     "freq_axis": ("<U32", "uniform_boson_centered"),
                     "bond_archive_layout": ("<U8", output_layout),
                     "block_layout": ("<U16", block_layout),
                     "green_layout": ("<U32", "b_f_k_o_o_v1"),
                     "archive_nonce": ("<U32", nonce),
                     "green_file": (None, "green.npz"),
                     "momentum_convention": (None, "e_plus_ikR"),
                     "norb_phys": (np.int64, 1), "nd": (np.int64, 2), "nmat": (np.int64, 4),
                     "B": (np.int64, 3), "beta": (np.float64, 2.0),
                     "scf_converged": (np.bool_, True)}
            for k, (dt, v) in fixed.items():
                a = z[k]
                self.assertEqual(a.ndim, 0, k)
                if dt is not None:
                    self.assertEqual(a.dtype, np.dtype(dt), k)
                self.assertEqual(a[()], v, k)
            self.assertEqual(len(nonce), 32)
            self.assertTrue(all(c in "0123456789abcdef" for c in nonce))
            self.assertEqual(z["cell_shape"].dtype, np.int64)
            np.testing.assert_array_equal(z["cell_shape"], [4, 1, 1])
            self.assertEqual(z["types"].dtype, np.dtype("<U12"))
            np.testing.assert_array_equal(z["types"], ["CoulombIntra"])
            self.assertEqual(z["wavevector_unit"].dtype, np.float64)
            np.testing.assert_array_equal(z["wavevector_unit"],
                                          kw["momentum"]["wavevector_unit"])
            self.assertEqual(z["wavevector_index"].dtype, np.int64)
            np.testing.assert_array_equal(z["wavevector_index"],
                                          kw["momentum"]["wavevector_index"])
            for k, v in kw["provenance"].items():
                np.testing.assert_array_equal(z[k], v, err_msg=k)
                self.assertEqual(z[k].dtype, np.asarray(v).dtype, k)
            if block_layout == "sz_packed":
                for k, v in so_bond.sector_map_members(inputs.sector_map).items():
                    np.testing.assert_array_equal(z[k], v, err_msg=k)
            arrays = {"Gamma_bond": kw["gamma_b"], "chi_w": kw["chi_w"]}
            if chibar:
                arrays["chibar_w"] = kw["chibar_w"]
            sidecars = []
            for p, arr in arrays.items():
                parts = ({p: arr} if block_layout == "dense"
                         else so_bond.packed_members(p, arr))
                for b, want in parts.items():
                    if output_layout == "sidecar" and p != "Gamma_bond":
                        name = str(z[b + "_file"][()])
                        self.assertEqual(name, "{}.{}.npy".format(b, nonce))
                        full = os.path.join(d, name)
                        self.assertEqual(int(z[b + "_nbytes"][()]), os.path.getsize(full))
                        self.assertEqual(z[b + "_nbytes"].dtype, np.int64)
                        got = np.load(full)
                        sidecars.append(name)
                    else:
                        got = z[b]
                    self.assertEqual(got.dtype, np.complex128, b)
                    np.testing.assert_array_equal(got, want, err_msg=b)
        src = self.assert_valid(path, block_layout, output_layout)
        for block, (where, _shape) in src.items():
            want_side = output_layout == "sidecar" and block.split("@")[0] != "Gamma_bond"
            self.assertEqual(where != "inline", want_side, block)
        self.assertEqual(sorted(os.listdir(d)), sorted(["lb.npz"] + sidecars))
        self.assertEqual(out, {"path": path, "archive_nonce": nonce,
                               "bond_archive_layout": output_layout,
                               "block_layout": block_layout, "sidecars": sorted(sidecars),
                               "removed": [], "kept": [], "kept_by": [], "nonso": []})

    def test_extra_members_written(self):
        path, _ = self.write(output_layout="npz", extra_members=dict(POLE))
        with np.load(path) as z:
            for k, v in POLE.items():
                np.testing.assert_array_equal(z[k], v)
                self.assertEqual(z[k].dtype, np.asarray(v).dtype)
        self.assert_valid(path)

    def test_scf_converged_false(self):
        path, _ = self.write(output_layout="npz", scf_converged=False)
        with np.load(path) as z:
            self.assertIs(bool(z["scf_converged"][()]), False)


class TestNoCopy(_Dir):
    """A 64 MiB ``chi_w``: the sidecar write peaks at <= 1 MiB of traced
    allocations, the npz write at <= min(M, 16 MiB) + 1 MiB with M the
    largest inline member."""

    def _peak(self, output_layout, inputs, name):
        d = self.sub(name)
        np.save(os.path.join(d, "warm.npy"), np.zeros(8))
        np.savez(os.path.join(d, "warm.npz"), a=np.zeros(8))
        os.remove(os.path.join(d, "warm.npy"))
        os.remove(os.path.join(d, "warm.npz"))
        tracemalloc.start()
        try:
            tracemalloc.reset_peak()
            base = tracemalloc.get_traced_memory()[0]
            path, _ = self.write(d=d, inputs=inputs, output_layout=output_layout)
            peak = tracemalloc.get_traced_memory()[1] - base
        finally:
            tracemalloc.stop()
        return path, peak

    def test_no_copy(self):
        cell = (8, 8, 8)
        nmat = 512
        inputs = fx.archive_inputs(B=1, cell=cell, nmat=nmat, chi_fill=1.0 + 0.5j)
        self.assertEqual(inputs.chi_w.nbytes, 64 * 2 ** 20)
        path, peak = self._peak("sidecar", inputs, "sidecar")
        self.assertLessEqual(peak, 2 ** 20, peak)
        self.assert_valid(path)
        M = sa.inline_max_bytes(nmat=nmat, nvol=512, B=1, nd=2, block_layout="dense",
                                output_layout="npz", n_transfer=1, n_rec=1, n_topology=1)
        self.assertEqual(M, 64 * 2 ** 20)
        path, peak = self._peak("npz", inputs, "npz")
        self.assertLessEqual(peak, min(M, 16 * 2 ** 20) + 2 ** 20, peak)
        self.assert_valid(path)


# ------------------------------------------------------------ atomic index
class TestAtomicIndex(_Dir):
    def _failures(self, path):
        exc = OSError("injected replace")

        def replace():
            return mock.patch.object(os, "replace", _selective(_REAL_REPLACE, path, exc)), exc

        calls = {"n": 0}
        exc2 = OSError("injected save")

        def partial_save(file, arr, *a, **kw):
            calls["n"] += 1
            if calls["n"] == 2:
                with open(file, "wb") as fh:
                    fh.write(b"\x93NUMPY partial")
                raise exc2
            return _REAL_SAVE(file, arr, *a, **kw)

        exc3 = RuntimeError("injected companion")

        def companion():
            raise exc3

        return [("replace", replace()[0], {}, exc),
                ("second_sector_save", mock.patch.object(np, "save", partial_save), {}, exc2),
                ("companion", contextlib.nullcontext(),
                 {"commit_companion": companion}, exc3)], calls

    def test_atomic_index(self):
        for previous in (False, True):
            d = self.sub("prev-{}".format(previous))
            path = os.path.join(d, "lb.npz")
            inputs = fx.archive_inputs(block_layout="sz_packed")
            if previous:
                self.write(d=d, inputs=inputs)
            rows, calls = self._failures(path)
            for name, patch, extra, exc in rows:
                with self.subTest(previous=previous, failure=name):
                    calls["n"] = 0
                    before = _listing(d)
                    with patch:
                        with self.assertRaises(type(exc)) as cm:
                            self.write(d=d, inputs=inputs, **extra)
                    self.assertIs(cm.exception, exc)
                    self.assertEqual(_listing(d), before)


class TestPostCommit(_Dir):
    def test_post_commit_failure_keeps_archive(self):
        path, first = self.write()
        with mock.patch.object(sa, "_scan_candidates",
                               side_effect=RuntimeError("injected scan")):
            with self.assertLogs(LOG, "WARNING") as cm:
                path2, out = self.write()
        self.assertEqual(path2, path)
        self.assertEqual(cm.output, ["WARNING:{}:{}".format(
            LOG.name, sa.T_OLD_SIDECARS_SCAN_FAILED.format(path=path,
                                                            error="injected scan"))])
        self.assert_valid(path)
        for n in out["sidecars"]:
            self.assertTrue(os.path.isfile(os.path.join(self.dir, n)))


class TestCompanion(_Dir):
    def test_commit_order(self):
        """``commit_companion`` runs after the temporary index is written and
        before the index replaces the previous one."""
        seen = []
        path = os.path.join(self.dir, "lb.npz")

        def companion():
            tmps = [n for n in os.listdir(self.dir) if n.endswith(".tmp.npz")]
            seen.append((len(tmps), os.path.exists(path)))

        self.write(commit_companion=companion)
        self.assertEqual(seen, [(1, False)])
        self.assertFalse([n for n in os.listdir(self.dir) if n.endswith(".tmp.npz")])

    def test_companion_failure_removes_old_index(self):
        for row in ("old_index", "no_old_index", "removal_fails"):
            with self.subTest(row=row):
                d = self.sub(row)
                path = os.path.join(d, "lb.npz")
                green = os.path.normpath(os.path.join(d, "green.npz"))
                if row != "no_old_index":
                    self.write(d=d)
                old = _sidecar_names(path) if row != "no_old_index" else []
                exc = OSError("injected index commit")
                committed = []
                patches = [mock.patch.object(os, "replace",
                                             _selective(_REAL_REPLACE, path, exc))]
                if row == "removal_fails":
                    exc_rm = PermissionError("injected removal")
                    patches.append(mock.patch.object(os, "remove",
                                                     _selective(_REAL_REMOVE, path, exc_rm)))
                with patches[0], (patches[1] if len(patches) > 1
                                  else contextlib.nullcontext()):
                    if row == "no_old_index":
                        with self.assertNoLogs(LOG):
                            with self.assertRaises(OSError) as cm:
                                self.write(d=d, commit_companion=lambda: committed.append(1))
                    else:
                        with self.assertLogs(LOG, "WARNING") as lg:
                            with self.assertRaises(OSError) as cm:
                                self.write(d=d, commit_companion=lambda: committed.append(1))
                self.assertIs(cm.exception, exc)
                self.assertEqual(committed, [1])
                if row == "old_index":
                    self.assertEqual(lg.output, ["WARNING:{}:{}".format(
                        LOG.name, sa.T_OLD_INDEX_REMOVED.format(path=path, green=green,
                                                                error=exc))])
                    self.assertFalse(os.path.exists(path))
                elif row == "removal_fails":
                    self.assertEqual(lg.output, ["WARNING:{}:{}".format(
                        LOG.name, sa.T_OLD_INDEX_KEPT.format(path=path, green=green,
                                                             error=exc_rm))])
                    self.assertTrue(os.path.exists(path))
                # the new sidecars and the temporary index are gone; the old
                # sidecars stay
                self.assertEqual(sorted(os.listdir(d)),
                                 sorted(old + (["lb.npz"] if row == "removal_fails" else [])))


class TestCleanupErrors(_Dir):
    def test_cleanup_error_does_not_mask(self):
        exc = OSError("injected save")
        target = {}

        def failing_save(file, arr, *a, **kw):
            target["path"] = file
            with open(file, "wb") as fh:
                fh.write(b"partial")
            raise exc

        def failing_remove(p, *a, **kw):
            if "path" in target and os.path.abspath(p) == os.path.abspath(target["path"]):
                raise PermissionError("injected removal")
            return _REAL_REMOVE(p, *a, **kw)

        with mock.patch.object(np, "save", failing_save), \
                mock.patch.object(os, "remove", failing_remove):
            with self.assertLogs(LOG, "WARNING") as lg:
                with self.assertRaises(OSError) as cm:
                    self.write()
        self.assertIs(cm.exception, exc)
        self.assertEqual(lg.output, ["WARNING:{}:{}".format(
            LOG.name, sa.T_CLEANUP_FAILED.format(name=target["path"],
                                                 error="injected removal"))])
        self.assertEqual(sorted(os.listdir(self.dir)), [os.path.basename(target["path"])])


# ---------------------------------------------------------- old sidecars
class TestOldSidecars(_Dir):
    def test_old_sidecars_removed(self):
        for block_layout in ("dense", "sz_packed"):
            with self.subTest(block_layout=block_layout):
                d = self.sub(block_layout)
                inputs = fx.archive_inputs(block_layout=block_layout)
                path, first = self.write(d=d, inputs=inputs)
                for extra in ("x.npy", "chi_w.{}.npy".format("0" * 31)):
                    np.save(os.path.join(d, extra), np.zeros(2))
                with self.assertLogs(LOG, "INFO") as lg:
                    _, second = self.write(d=d, inputs=inputs)
                self.assertEqual(lg.output, [_old_info(path, removed=first["sidecars"])])
                self.assertEqual(second["removed"], first["sidecars"])
                self.assertEqual(sorted(os.listdir(d)), sorted(
                    ["lb.npz", "x.npy", "chi_w.{}.npy".format("0" * 31)]
                    + second["sidecars"]))
                self.assert_valid(path)

    def test_old_names_follow_the_nonce(self):
        """Only the replaced index's names ``<block>.<its nonce>.npy`` are old
        sidecars: a name of another nonce, or not a sidecar name at all, is
        never removed."""
        path, first = self.write()
        other = "chi_w.{}.npy".format("b" * 32)
        for n in (other, "x.npy"):
            np.save(os.path.join(self.dir, n), np.zeros(2))
        with np.load(path) as z:
            d = {k: z[k] for k in z.files}
        d["chi_w_file"] = np.str_(other)
        d["chibar_w_file"] = np.str_("x.npy")
        np.savez(path, **d)
        with self.assertNoLogs(LOG):
            _, second = self.write()
        self.assertEqual(second["removed"], [])
        self.assertEqual(sorted(os.listdir(self.dir)), sorted(
            ["lb.npz", other, "x.npy"] + first["sidecars"] + second["sidecars"]))

    def test_name_the_new_index_carries_is_kept(self):
        for block_layout in ("dense", "sz_packed"):
            with self.subTest(block_layout=block_layout):
                d = self.sub(block_layout)
                inputs = fx.archive_inputs(block_layout=block_layout)
                path, first = self.write(d=d, inputs=inputs)
                with mock.patch.object(sa, "_new_nonce", return_value=first["archive_nonce"]):
                    with self.assertLogs(LOG, "INFO") as lg:
                        _, second = self.write(d=d, inputs=inputs)
                self.assertEqual(second["sidecars"], first["sidecars"])
                self.assertEqual(lg.output, [_old_info(path, kept=first["sidecars"])])
                self.assertEqual(sorted(os.listdir(d)), sorted(["lb.npz"] + first["sidecars"]))
                self.assert_valid(path)

    def test_candidate_matrix(self):
        """The other ``*.npz`` files in the directory: a hidden partial
        temporary is skipped; a candidate is an index only with a 0-d
        integer ``bond_archive_schema == 4``; an unreadable or malformed
        candidate keeps every old name and is named in ``kept by``; a valid
        index keeps the names it carries."""
        inputs = fx.archive_inputs(block_layout="sz_packed")

        def schema(value, file_value=None):
            def make(d, old):
                m = {"bond_archive_schema": value,
                     "chi_w@0_file": np.str_(old[0]) if file_value is None else file_value}
                np.savez(os.path.join(d, "other.npz"), **m)
                return "other.npz"
            return make

        def raw(name, data):
            def make(d, old):
                with open(os.path.join(d, name), "wb") as fh:
                    fh.write(data)
                return name
            return make

        rows = [
            ("hidden_tmp", raw(".x.1-ab.tmp.npz", b"PK\x03\x04partial"), "removed"),
            ("junk", raw("junk.npz", b"not a zip file at all"), "all"),
            ("schema_float", schema(np.float64(4.0)), "removed"),
            ("schema_str", schema(np.str_("4")), "removed"),
            ("schema_list", schema(np.array([4])), "removed"),
            ("schema_bool", schema(np.bool_(True)), "removed"),
            ("nonstring_file", schema(np.int64(4), np.int64(3)), "all"),
            ("index", schema(np.int64(4)), "one"),
        ]
        for name, make, expect in rows:
            with self.subTest(candidate=name):
                d = self.sub(name)
                path, first = self.write(d=d, inputs=inputs)
                old = sorted(first["sidecars"])
                cand = make(d, old)
                with self.assertLogs(LOG, "INFO") as lg:
                    _, second = self.write(d=d, inputs=inputs)
                if expect == "removed":
                    removed, kept, by = old, [], []
                elif expect == "all":
                    removed, kept, by = [], old, [cand]
                else:
                    removed, kept, by = old[1:], old[:1], [cand]
                self.assertEqual(lg.output, [_old_info(path, removed, kept, by)])
                self.assertEqual((second["removed"], second["kept"], second["kept_by"]),
                                 (removed, kept, by))
                for n in kept + second["sidecars"]:
                    self.assertTrue(os.path.isfile(os.path.join(d, n)), n)
                for n in removed:
                    self.assertFalse(os.path.exists(os.path.join(d, n)), n)
                self.assertTrue(os.path.exists(os.path.join(d, cand)))
                self.assert_valid(path)

    def test_copied_index_keeps_sidecars(self):
        a, first = self.write("A.npz")
        c = os.path.join(self.dir, "C.npz")
        shutil.copy(a, c)
        with self.assertLogs(LOG, "INFO") as lg:
            _, second = self.write("C.npz")
        self.assertEqual(lg.output, [_old_info(c, kept=first["sidecars"], kept_by=["A.npz"])])
        self.assert_valid(a)
        self.assert_valid(c)

    def test_sidecar_then_npz(self):
        path, first = self.write()
        with self.assertLogs(LOG, "INFO") as lg:
            self.write(output_layout="npz")
        self.assertEqual(lg.output, [_old_info(path, removed=first["sidecars"])])
        self.assertEqual(os.listdir(self.dir), ["lb.npz"])
        self.assert_valid(path, output_layout="npz")

    def test_nonso_sidecars_listed(self):
        for previous in ("none", "nonso_index", "so_index"):
            with self.subTest(previous=previous):
                d = self.sub(previous)
                path = os.path.join(d, "lb.npz")
                nonso = ["lb_chi_c_w.npy", "lb_chi_s_w.npy"]
                for n in nonso:
                    np.save(os.path.join(d, n), np.zeros(2))
                np.save(os.path.join(d, "other_chi_s_w.npy"), np.zeros(2))
                removed = []
                if previous == "nonso_index":
                    np.savez(path, bond_archive_schema=np.int64(3),
                             bond_archive_layout=np.str_("sidecar"),
                             chi_s_w_file=np.str_(nonso[1]), chi_c_w_file=np.str_(nonso[0]))
                elif previous == "so_index":
                    removed = self.write(d=d)[1]["sidecars"]
                with self.assertLogs(LOG, "INFO") as lg:
                    self.write(d=d, output_layout="npz")
                self.assertEqual(lg.output, [_old_info(path, removed=removed, nonso=nonso,
                                                       replaced=previous != "none")])
                for n in nonso:
                    self.assertTrue(os.path.isfile(os.path.join(d, n)))

    def test_unlink_failure_warns(self):
        path, first = self.write()
        target = os.path.join(self.dir, first["sidecars"][0])
        exc = PermissionError("injected unlink")
        with mock.patch.object(os, "remove", _selective(_REAL_REMOVE, target, exc)):
            with self.assertLogs(LOG, "INFO") as lg:
                _, second = self.write()
        self.assertEqual(lg.output, ["WARNING:{}:{}".format(
            LOG.name, sa.T_OLD_SIDECARS_FAILED.format(path=path, name=target, error=exc))])
        self.assertEqual(second["removed"], [])
        self.assertTrue(os.path.isfile(target))
        self.assert_valid(path)

    def test_no_log_when_nothing_to_do(self):
        """A fresh write, and an npz archive replacing an npz archive, log nothing."""
        for output_layout in ("sidecar", "npz"):
            with self.subTest(output_layout=output_layout):
                d = self.sub(output_layout)
                with self.assertNoLogs(LOG):
                    self.write(d=d, output_layout=output_layout)
        with self.assertNoLogs(LOG):
            self.write(d=os.path.join(self.dir, "npz"), output_layout="npz")


# ------------------------------------------------------------ preconditions
class TestRefusals(_Dir):
    @classmethod
    def setUpClass(cls):
        cls.dense = fx.archive_inputs()
        cls.packed = fx.archive_inputs(block_layout="sz_packed")

    def refused(self, text, *, inputs=None, name="lb.npz", output_layout="sidecar",
                **changes):
        inputs = inputs or self.dense
        args = dict(inputs.kwargs)
        args.update(changes)
        with self.assertRaises(ValueError) as cm:
            sa.write_so_bond_archive(os.path.join(self.dir, name),
                                     output_layout=output_layout, log=LOG, **args)
        self.assertEqual(str(cm.exception), text)
        self.assertEqual([n for n in os.listdir(self.dir) if n not in ("sub", "g.npz")], [])

    def test_green_file_refused(self):
        os.makedirs(os.path.join(self.dir, "sub"))
        absolute = os.path.join(self.dir, "g.npz")
        for g in ("", ".", "sub", "sub/", absolute, "a/../g.npz", "../g.npz", "./"):
            with self.subTest(green_file=g):
                self.refused("write_so_bond_archive: green_file must be a relative file path "
                             "without '..', got {!r}".format(g), green_file=g)

    def test_green_file_rule(self):
        os.makedirs(os.path.join(self.dir, "sub"))
        ok = ["green.npz", "sub/g.npz", "./g.npz", "a b.npz"]
        bad = ["", ".", "sub", "sub/", "/abs/g.npz", "a/../g.npz", "..", "./", 3, None]
        for g in ok:
            self.assertTrue(sa.green_file_ok(g, self.dir), g)
        for g in bad:
            self.assertFalse(sa.green_file_ok(g, self.dir), g)

    def test_green_file_below_accepted(self):
        path, _ = self.write(green_file="sub/g.npz", output_layout="npz")
        with np.load(path) as z:
            self.assertEqual(str(z["green_file"][()]), "sub/g.npz")

    def test_path_must_end_npz(self):
        for name in ("lb.npy", "lb", "lb.NPZ"):
            with self.subTest(name=name):
                self.refused("write_so_bond_archive: path must end in .npz, got {!r}".format(
                    os.path.join(self.dir, name)), name=name)

    def test_duplicate_scalar_refused(self):
        for key in ("delta_r", "reverse", "slot_convention", "max_shells", "state_beta"):
            with self.subTest(key=key):
                s = dict(self.dense.kwargs["scalars"])
                s[key] = self.dense.kwargs["provenance"][key]
                self.refused("write_so_bond_archive: scalars must not carry {} (it comes "
                             "from provenance only)".format(key), scalars=s)

    def test_preconditions(self):
        kw = self.dense.kwargs
        prov = kw["provenance"]

        def scal(**ch):
            s = dict(kw["scalars"])
            for k, v in ch.items():
                if v is _DROP:
                    del s[k]
                else:
                    s[k] = v
            return s

        def provc(**ch):
            p = dict(prov)
            for k, v in ch.items():
                if v is _DROP:
                    del p[k]
                else:
                    p[k] = v
            return p

        def mom(**ch):
            m = dict(kw["momentum"])
            m.update(ch)
            return m

        W = "write_so_bond_archive: "
        chi = kw["chi_w"]
        rows = [
            ("block_layout", dict(block_layout="packed"),
             W + "block_layout must be 'dense' or 'sz_packed', got 'packed'"),
            ("output_layout", dict(output_layout="zip"),
             W + "output_layout must be 'npz' or 'sidecar', got 'zip'"),
            ("scalar keys", dict(scalars=scal(types=_DROP)),
             W + "scalars must have exactly the keys norb_phys, nd, nmat, B, cell_shape, "
             "beta, types, got B, beta, cell_shape, nd, nmat, norb_phys"),
            ("norb_phys", dict(scalars=scal(norb_phys=0)),
             W + "scalars['norb_phys'] must be an integer >= 1, got 0"),
            ("norb_phys bool", dict(scalars=scal(norb_phys=True)),
             W + "scalars['norb_phys'] must be an integer >= 1, got True"),
            ("nd", dict(scalars=scal(nd=4)),
             W + "scalars['nd'] must be an integer equal to 2 * norb_phys, got 4"),
            ("nmat", dict(scalars=scal(nmat=5)),
             W + "scalars['nmat'] must be an integer > 0 and even, got 5"),
            ("nmat float", dict(scalars=scal(nmat=4.0)),
             W + "scalars['nmat'] must be an integer > 0 and even, got 4.0"),
            ("B", dict(scalars=scal(B=0)), W + "scalars['B'] must be an integer >= 1, got 0"),
            ("beta", dict(scalars=scal(beta=float("inf"))),
             W + "scalars['beta'] must be a finite float > 0, got inf"),
            ("beta int", dict(scalars=scal(beta=2)),
             W + "scalars['beta'] must be a finite float > 0, got 2"),
            ("cell_shape", dict(scalars=scal(cell_shape=np.array([4, 1]))),
             W + "scalars['cell_shape'] must be an integer array of shape (3,) with entries "
             ">= 1, got [4, 1]"),
            ("types dtype", dict(scalars=scal(types=np.array(["CoulombIntra"], "<U20"))),
             W + "scalars['types'] must be a sorted 1-d <U12 array of distinct names, got "
             "['CoulombIntra']"),
            ("types order", dict(scalars=scal(types=np.array(["Hund", "CoulombIntra"],
                                                             "<U12"))),
             W + "scalars['types'] must be a sorted 1-d <U12 array of distinct names, got "
             "['Hund', 'CoulombIntra']"),
            ("provenance keys", dict(provenance=provc(max_shells=_DROP, extra=np.int64(1))),
             W + "provenance must carry the spin-orbital and bond provenance members; "
             "missing max_shells; unexpected extra"),
            ("stage", dict(provenance=provc(producer_stage=np.str_("flex_so_onsite"))),
             W + "provenance producer_stage must be 'flex_so_bond', got 'flex_so_onsite'"),
            ("index convention", dict(provenance=provc(
                prov_index_convention=np.str_("interleaved"))),
             W + "provenance prov_index_convention must be 'spin_block', got 'interleaved'"),
            ("relation B", dict(scalars=scal(B=1)),
             W + "scalars disagree with provenance: B = 1 but len(delta_r) = 3"),
            ("relation norb_phys", dict(scalars=scal(norb_phys=2, nd=4)),
             W + "scalars disagree with provenance: norb_phys = 2 but prov_norb_phys = 1"),
            ("relation nmat", dict(scalars=scal(nmat=6)),
             W + "scalars disagree with provenance: nmat = 6 but state_nmat = 4"),
            ("relation beta", dict(scalars=scal(beta=3.0)),
             W + "scalars disagree with provenance: beta = 3.0 but state_beta = 2.0"),
            ("relation cell", dict(scalars=scal(cell_shape=np.array([2, 2, 1]))),
             W + "scalars disagree with provenance: cell_shape = [2, 2, 1] but "
             "prov_cell_shape = [4, 1, 1]"),
            ("dense map", dict(sector_map=self.packed.sector_map),
             W + "a dense archive takes sector_map = None, got SectorMap"),
            ("packed map", dict(sector_map=None, inputs=self.packed),
             W + "a sz_packed archive needs the sector map of the provenance's (reverse, "
             "nd), got NoneType"),
            ("packed other map", dict(inputs=self.packed, sector_map=so_bond.sz_sector_map(
                fx.bond_inputs(1), 2)),
             W + "a sz_packed archive needs the sector map of the provenance's (reverse, "
             "nd), got another SectorMap"),
            ("dense type", dict(chi_w=list(chi.shape)),
             W + "chi_w of a dense archive must be a numpy array, got list"),
            ("packed type", dict(inputs=self.packed, chi_w=self.packed.chi_w),
             W + "chi_w of a sz_packed archive must be a PackedArrays, got ndarray"),
            ("packed on other map", dict(inputs=self.packed,
                                         gamma_b=so_bond.PackedArrays.zeros(
                                             so_bond.sz_sector_map(fx.bond_inputs(1), 2), (4,))),
             W + "gamma_b is packed on another sector map"),
            ("dtype", dict(chi_w=chi.astype(np.complex64)),
             W + "chi_w must be complex128, got complex64"),
            ("byte order", dict(chi_w=chi.astype(">c16")),
             W + "chi_w must be complex128, got >c16"),
            ("order", dict(chi_w=np.asfortranarray(chi)), W + "chi_w must be C-contiguous"),
            ("shape", dict(chi_w=chi[:2]),
             W + "chi_w must have shape (4, 4, 12, 12), got (2, 4, 12, 12)"),
            ("gamma shape", dict(gamma_b=np.ascontiguousarray(kw["gamma_b"][:, :6, :6])),
             W + "gamma_b must have shape (4, 12, 12), got (4, 6, 6)"),
            ("chibar shape", dict(chibar_w=np.ascontiguousarray(chi[:, :2])),
             W + "chibar_w must have shape (4, 4, 12, 12), got (4, 2, 12, 12)"),
            ("packed block order", dict(inputs=self.packed, chi_w=so_bond.PackedArrays(
                self.packed.sector_map, tuple(np.asfortranarray(b)
                                              for b in self.packed.kwargs["chi_w"].blocks))),
             W + "chi_w@0 must be C-contiguous"),
            ("packed lead", dict(inputs=self.packed, gamma_b=so_bond.PackedArrays.zeros(
                self.packed.sector_map, (3,))),
             W + "gamma_b@0 must have shape (4, 6, 6), got (3, 6, 6)"),
            ("momentum keys", dict(momentum={"wavevector_unit": np.eye(3)}),
             W + "momentum must have exactly the keys momentum_convention, wavevector_unit, "
             "wavevector_index, got wavevector_unit"),
            ("momentum convention", dict(momentum=mom(momentum_convention="e_minus_ikR")),
             W + "momentum_convention must be 'e_plus_ikR', got 'e_minus_ikR'"),
            ("wavevector_unit", dict(momentum=mom(wavevector_unit=np.eye(2))),
             W + "wavevector_unit must be a finite real array of shape (3, 3), got float64 "
             "(2, 2)"),
            ("wavevector_index", dict(momentum=mom(wavevector_index=np.zeros((3, 3), int))),
             W + "wavevector_index must be an integer array of shape (4, 3), got int64 (3, 3)"),
            ("scf_converged", dict(scf_converged=1),
             W + "scf_converged must be a bool, got 1"),
            ("commit_companion", dict(commit_companion="yes"),
             W + "commit_companion must be callable or None, got str"),
            ("extra collision", dict(extra_members={"nd": np.int64(2)}),
             W + "extra_members key 'nd' collides with an archive member"),
        ]
        for name, changes, text in rows:
            with self.subTest(row=name):
                changes = dict(changes)
                inputs = changes.pop("inputs", None)
                output_layout = changes.pop("output_layout", "sidecar")
                self.refused(text, inputs=inputs, output_layout=output_layout, **changes)

    def test_extra_members_group(self):
        W = "write_so_bond_archive: extra_members is not the pole group: "
        rows = [
            ("partial", {"rpa_past_pole": np.bool_(True)},
             "pole member rpa_pole_indicator_spinful is missing (absent vs present)"),
            ("past false", dict(POLE, rpa_past_pole=np.bool_(False)),
             "rpa_past_pole is not true (False vs True)"),
            ("past dtype", dict(POLE, rpa_past_pole=np.int64(1)),
             "rpa_past_pole is not true (1 vs True)"),
            ("indicator dtype", dict(POLE, rpa_pole_indicator_spinful=np.float32(1.2)),
             "rpa_pole_indicator_spinful is not a finite float64 scalar (float32/0d vs "
             "float64/0d)"),
            ("indicator finite", dict(POLE, rpa_pole_indicator_spinful=np.float64("nan")),
             "rpa_pole_indicator_spinful is not a finite float64 scalar (nan vs finite)"),
            ("q shape", dict(POLE, rpa_pole_q_spinful=np.zeros(2)),
             "rpa_pole_q_spinful is not a finite float64 (3,) array (float64/(2,) vs "
             "float64/(3,))"),
            ("q finite", dict(POLE, rpa_pole_q_spinful=np.array([0.0, np.inf, 0.0])),
             "rpa_pole_q_spinful is not a finite float64 (3,) array (inf vs finite)"),
            ("extra key", dict(POLE, rpa_pole_indicator_charge=np.float64(0.5)),
             "member rpa_pole_indicator_charge is not a pole member (present vs absent)"),
        ]
        for name, extra, tail in rows:
            with self.subTest(row=name):
                self.refused(W + tail, extra_members=extra)

    def test_pole_group_rule(self):
        self.assertIsNone(sa.pole_group_problem({}))
        self.assertIsNone(sa.pole_group_problem(dict(POLE)))
        self.assertEqual(sa.pole_group_problem({"rpa_pole_q_spinful": POLE["rpa_pole_q_spinful"]}),
                         ("pole member rpa_past_pole is missing", "absent", "present"))


_DROP = object()


# ------------------------------------------------------------------ layout
class TestArchiveLayout(unittest.TestCase):
    MB = {"chi_w@0": 600, "chi_w@+1": 300, "chi_w@-1": 300}

    def test_rows(self):
        L = sa.archive_layout
        self.assertEqual(L(member_bytes=self.MB, small_array_bytes=1000, request=None),
                         sa.ArchiveLayout("npz", False, "chi_w@0", 600, 1000))
        self.assertEqual(L(member_bytes=self.MB, small_array_bytes=500, request=None),
                         sa.ArchiveLayout("sidecar", True, "chi_w@0", 600, 500))
        self.assertEqual(L(member_bytes=self.MB, small_array_bytes=600, request=None).layout,
                         "npz")
        for req in ("sidecar", " SideCar "):
            self.assertEqual(L(member_bytes=self.MB, small_array_bytes=500, request=req),
                             sa.ArchiveLayout("sidecar", False, "chi_w@0", 600, 500))
            self.assertEqual(L(member_bytes=self.MB, small_array_bytes=1000, request=req),
                             sa.ArchiveLayout("sidecar", False, "chi_w@0", 600, 1000))
        for req in ("npz", "NPZ", " Npz"):
            self.assertEqual(L(member_bytes=self.MB, small_array_bytes=1000, request=req),
                             sa.ArchiveLayout("npz", False, "chi_w@0", 600, 1000))
        # a tie: the first name in sorted order
        self.assertEqual(L(member_bytes={"chi_w@-1": 5, "chi_w@+1": 5}, small_array_bytes=1,
                           request=None).member, "chi_w@+1")

    def test_npz_refused(self):
        big = {"chi_w": 3 * 2 ** 30}
        with self.assertRaises(ValueError) as cm:
            sa.archive_layout(member_bytes=big, small_array_bytes=float(2 ** 30),
                              request=" NPZ ")
        self.assertEqual(str(cm.exception), sa.T_NPZ_TOO_LARGE.format(
            member="chi_w", size="3.00 GiB", bound="1.00 GiB"))
        self.assertEqual(str(cm.exception), (
            "[mode] enable_spin_orbital = true with [mode.param] longitudinal_bond_channels "
            "= true (FLEX): longitudinal_bond_output_layout = \"npz\" cannot hold the bond "
            "archive member chi_w (3.00 GiB, above the small-array bound 1.00 GiB = "
            "min(1 GiB, cap / 8)); set longitudinal_bond_output_layout = \"sidecar\" or "
            "remove it (the sidecar layout is then chosen)."))

    def test_bad_request(self):
        for req in ("zip", 1, ""):
            with self.subTest(request=req):
                with self.assertRaises(ValueError) as cm:
                    sa.archive_layout(member_bytes=self.MB, small_array_bytes=1, request=req)
                self.assertEqual(str(cm.exception), "archive_layout: request must be None, "
                                 "'npz' or 'sidecar', got {!r}".format(req))
        with self.assertRaises(ValueError):
            sa.archive_layout(member_bytes={}, small_array_bytes=1, request=None)

    def test_format_bytes(self):
        self.assertEqual(sa.format_bytes(0), "0.00 MiB")
        self.assertEqual(sa.format_bytes(3 * 2 ** 19), "1.50 MiB")
        self.assertEqual(sa.format_bytes(2 ** 30 - 1), "1024.00 MiB")
        self.assertEqual(sa.format_bytes(2 ** 30), "1.00 GiB")
        self.assertEqual(sa.format_bytes(float(2 ** 30) * 2.5), "2.50 GiB")

    def test_slot_member_bytes(self):
        self.assertEqual(sa.slot_member_bytes(nmat=8, nvol=9, B=5, nd=4, layout="dense"),
                         {"chi_w": 16 * 8 * 9 * 80 ** 2})
        self.assertEqual(sa.slot_member_bytes(nmat=8, nvol=9, B=5, nd=4, layout="sz_packed"),
                         {"chi_w@0": 16 * 8 * 9 * 40 ** 2, "chi_w@+1": 16 * 8 * 9 * 20 ** 2,
                          "chi_w@-1": 16 * 8 * 9 * 20 ** 2})
        self.assertEqual(list(sa.slot_member_bytes(nmat=2, nvol=1, B=1, nd=2,
                                                   layout="sz_packed")),
                         ["chi_w@0", "chi_w@+1", "chi_w@-1"])
        with self.assertRaises(ValueError):
            sa.slot_member_bytes(nmat=2, nvol=1, B=1, nd=2, layout="npz")

    def test_inline_max_bytes(self):
        f = sa.inline_max_bytes
        # npz: the slot member; sidecar: Gamma_bond
        self.assertEqual(f(nmat=8, nvol=9, B=5, nd=4, block_layout="dense",
                           output_layout="npz"), 16 * 8 * 9 * 80 ** 2)
        self.assertEqual(f(nmat=8, nvol=9, B=5, nd=4, block_layout="dense",
                           output_layout="sidecar"), 16 * 9 * 80 ** 2)
        self.assertEqual(f(nmat=8, nvol=9, B=5, nd=4, block_layout="sz_packed",
                           output_layout="sidecar"), 16 * 9 * 40 ** 2)
        self.assertEqual(f(nmat=8, nvol=9, B=5, nd=4, block_layout="sz_packed",
                           output_layout="npz"), 16 * 8 * 9 * 40 ** 2)
        # tiny blocks: the record and transfer tables
        self.assertEqual(f(nmat=2, nvol=100, B=1, nd=2, block_layout="dense",
                           output_layout="sidecar"), 16 * 100 * 16)
        self.assertEqual(f(nmat=2, nvol=1000, B=1, nd=2, block_layout="sz_packed",
                           output_layout="sidecar"), 16 * 1000 * 2 ** 2)
        self.assertEqual(f(nmat=2, nvol=1, B=1, nd=2, block_layout="dense",
                           output_layout="sidecar", n_rec=1000), 96 * 1000)
        self.assertEqual(f(nmat=2, nvol=1, B=1, nd=2, block_layout="dense",
                           output_layout="sidecar", n_transfer=1000), 40 * 1000)
        self.assertIsInstance(f(nmat=2, nvol=1, B=1, nd=2, block_layout="dense",
                                output_layout="npz"), int)

    def test_written_members_bounded(self):
        """``inline_max_bytes`` bounds every inline member of a written archive."""
        with tempfile.TemporaryDirectory() as d:
            for block_layout in ("dense", "sz_packed"):
                for output_layout in ("npz", "sidecar"):
                    inputs = fx.archive_inputs(block_layout=block_layout)
                    p = os.path.join(d, "{}-{}.npz".format(block_layout, output_layout))
                    sa.write_so_bond_archive(p, output_layout=output_layout, log=LOG,
                                             **inputs.kwargs)
                    prov = inputs.kwargs["provenance"]
                    M = sa.inline_max_bytes(
                        nmat=4, nvol=4, B=3, nd=2, block_layout=block_layout,
                        output_layout=output_layout,
                        n_transfer=len(prov["prov_transfer_vals"]),
                        n_extern=len(prov["prov_extern_vals"]),
                        n_rec=len(prov["prov_rec_v"]),
                        n_topology=len(prov["prov_topology"]))
                    with np.load(p) as z:
                        biggest = max(z[n].nbytes for n in z.files)
                    self.assertEqual(biggest, M, (block_layout, output_layout))


class TestMessages(unittest.TestCase):
    def test_constants(self):
        H = ("[mode] enable_spin_orbital = true with [mode.param] longitudinal_bond_channels "
             "= true (FLEX)")
        self.assertEqual(sa.SO_BOND_FLEX_HEAD, H)
        self.assertEqual((sa.SCHEMA, sa.INDEX_ORDER, sa.FREQ_AXIS, sa.PRODUCER_LABEL),
                         (4, "bond_major_pair_code_slot", "uniform_boson_centered",
                          "hwave FLEX (spin-orbital, bond channels)"))
        self.assertEqual(sa.PRODUCER_LABEL, fx.LABEL)
        want = {
            "T_NOT_SO_ARCHIVE": "'{path}' is not a spin-orbital bond archive "
            "(bond_archive_schema {schema}, spin_orbital {so}); write it with FLEX, [mode] "
            "enable_spin_orbital = true, [mode.param] longitudinal_bond_channels = true and "
            "longitudinal_bond_output_full = true.",
            "T_COMPANION_MISSING": "companion green file '{green}' named by '{archive}' is "
            "missing (keep the archive, its sidecars and the green file together, or re-run "
            "FLEX)",
            "T_NEED_GREEN": H + ": longitudinal_bond_output_full = true needs [file.output] "
            "green (the bond archive names the Green function its susceptibility was built "
            "from as its companion green_file; [file.output] keys are case-sensitive); add "
            "green to [file.output].",
            "T_GREEN_DIR": H + ": the bond archive '{archive}' and its companion green file "
            "'{green}' must be in the same directory or the green file below it (green_file "
            "is stored as a relative path without '..').",
            "T_FORCED": "spin-orbital bond archive: sidecar layout (the member {member} of "
            "{size} is above the small-array bound {bound}); the susceptibility is written as "
            ".npy files next to the index.",
            "T_PASS": "spin-orbital bond archive: chi_w is chi_B of the final Green function "
            "(one extra susceptibility pass after the self-consistency).",
            "T_PASS_REFUSAL_NOTE": " This happened in the final susceptibility pass of the "
            "bond archive; set longitudinal_bond_output_full = false to skip that pass.",
            "T_OLD_SIDECARS": "spin-orbital bond archive '{path}': removed {n} sidecar file(s) "
            "of the archive it replaced: {removed}; kept {k} (named by another index, a "
            "malformed candidate or this one): {kept}; kept by: {kept_by}{nonso}",
            "T_OLD_SIDECARS_FAILED": "spin-orbital bond archive '{path}': could not remove the "
            "old sidecar file '{name}' ({error}); remove it by hand if no other archive uses "
            "it.",
            "T_ARCHIVE_KEY_IGNORED": "[file.output] longitudinal_bond is written only with "
            "[mode.param] longitudinal_bond_output_full = true; no bond archive is written by "
            "this run.",
            "T_NO_ARCHIVE_ITER0": "spin-orbital bond archive: IterationMax = 0 runs no map, so "
            "no bond archive is written.",
            "T_REMEDY_ARCHIVE": "longitudinal_bond_output_full = true adds the archive phases "
            "(archive_bubble, archive_dressing, archive_write); setting it to false removes "
            "them.",
            "T_OLD_INDEX_REMOVED": "spin-orbital bond archive '{path}': writing the new index "
            "failed after the companion green file '{green}' was replaced; the previous index "
            "was removed so that no archive names a different Green function ({error}).",
            "T_OLD_INDEX_KEPT": "spin-orbital bond archive '{path}': writing the new index "
            "failed after the companion green file '{green}' was replaced, and the previous "
            "index could not be removed ({error}); it names a different Green function, so "
            "load_so_bond_archive refuses it; remove it or re-run FLEX.",
            "T_CLEANUP_FAILED": "spin-orbital bond archive: could not remove '{name}' while "
            "cleaning up after a failed write ({error}); remove it by hand.",
            "T_NONSO_TAIL": "; the non-spin-orbital sidecars {names} at this stem are left in "
            "place",
        }
        for name, text in want.items():
            self.assertEqual(getattr(sa, name), text, name)
        self.assertEqual(sa.T_OLD_SIDECARS_NO_INDEX,
                         sa.T_OLD_SIDECARS.replace(" of the archive it replaced", ""))

    def test_no_flex_import(self):
        import subprocess
        import sys
        code = ("import sys; import hwave.solver.so_bond_archive; "
                "print('hwave.solver.flex' in sys.modules)")
        env = dict(os.environ)
        r = subprocess.run([sys.executable, "-B", "-c", code], capture_output=True,
                           text=True, env=env)
        self.assertEqual(r.stdout.strip(), "False", r.stderr)

    def test_sidecar_name_rule(self):
        n = "a" * 32
        self.assertEqual(sa.sidecar_name("chi_w@+1", n), "chi_w@+1.{}.npy".format(n))
        self.assertEqual(sa.parse_sidecar_name("chi_w@+1.{}.npy".format(n)), ("chi_w@+1", n))
        self.assertEqual(sa.parse_sidecar_name("Gamma_bond.{}.npy".format(n)),
                         ("Gamma_bond", n))
        for bad in ("chi_w.{}.npy".format("A" * 32), "chi_w@2.{}.npy".format(n),
                    "chi.{}.npy".format(n), "sub/chi_w.{}.npy".format(n),
                    "chi_w.{}.npz".format(n), "chi_w.{}.npy".format(n[:31])):
            self.assertIsNone(sa.parse_sidecar_name(bad), bad)


if __name__ == "__main__":
    unittest.main()
