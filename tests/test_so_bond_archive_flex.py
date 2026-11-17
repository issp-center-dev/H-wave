"""The schema-4 bond archive written by the spin-orbital bond FLEX (issue
#270): ``[mode] enable_spin_orbital = true`` with ``[mode.param]
longitudinal_bond_channels = true`` and ``longitudinal_bond_output_full =
true``.

The archive's ``chi_w`` is ``chi_B`` of the FINAL Green function (the one in
``green.npz``, named by ``green_file``), computed by one extra
susceptibility pass after the self-consistency. These tests run the solver
and ``save_results`` at the public entry and read the archive back with
``hwave.solver.so_bond_archive.load_so_bond_archive``:

* :class:`TestRoundTrip` -- every block layout and output layout is
  accepted, and the slot payloads read back equal the arrays handed to the
  writer (captured by a spy on the writer, never compared with the reader
  itself);
* :class:`TestContent` -- ``chi_w`` against a test-local oracle built from
  the companion ``green.npz`` (the production bubble and vertex, and an
  independent ``np.linalg.solve`` per ``(l, q)``), with and without a tail, its
  non-vacuity against ``chiq.npz`` and a break row;
* :class:`TestCompanion`, :class:`TestSidecar`, :class:`TestStructure`,
  :class:`TestCompatibility` -- the companion pair, the sidecar layout, the
  structural and compatibility refusals on copies of real archives;
* :class:`TestCrossRefusals` -- the non-spin-orbital loader, the on-site
  hwave_sc and the capability backstop on the real outputs;
* :class:`TestWiring` -- the pass, the pole check, the payload's lifetime,
  the output paths, the notices and the failure paths of the green commit.

The real archives are built once per module (``setUpModule``) and copied
by every row that changes one. Mutations patch module attributes looked up
at call time; every reference is computed before its patch is entered.
"""
import contextlib
import glob
import logging
import os
import shutil
import tempfile
import unittest
from typing import NamedTuple
from unittest import mock

import numpy as np

from tests import so_bond_entry_fixtures as F
from tests import so_bond_flex_fixtures as X
from tests.test_so_bond_archive import _DROP, _M, _Spies, _bounded, _members, _rewrite

LABEL = "hwave FLEX (spin-orbital, bond channels)"
H = ("[mode] enable_spin_orbital = true with [mode.param] longitudinal_bond_channels = true "
     "(FLEX)")
#: the gate tolerance of the archive's content and of the layout comparison
TOL = 1e-12
#: the floor of a deliberate break and of the non-vacuity of the tail
MISS = 1e-3
#: chi_w (final G) differs from chiq.npz (last map's G) at least this much
NONVAC = 1e-6

T_PASS = ("spin-orbital bond archive: chi_w is chi_B of the final Green function (one extra "
          "susceptibility pass after the self-consistency).")
T_PASS_REFUSAL_NOTE = (" This happened in the final susceptibility pass of the bond archive; "
                       "set longitudinal_bond_output_full = false to skip that pass.")
REMEDY = (" (remedy: lower [mode.param] longitudinal_bond_cond_tol, a finite number in (0, 1), "
          "default 1.0e-3; the dressing guard is not relaxed by flex_guard_policy)")
INFO_NO_ARCHIVE = (
    "spin-orbital bond FLEX: chi0q and chiq are the channel-0 blocks of chibar_B and chi_B of "
    "the last SCF iteration; the bond-resolved arrays are written only with "
    "longitudinal_bond_output_full = true.")
INFO_ARCHIVE = (
    "spin-orbital bond FLEX: chi0q and chiq are the channel-0 blocks of chibar_B and chi_B of "
    "the last SCF iteration; the bond archive holds chi_B of the final Green function.")
POLE_WHERE = ("FLEX longitudinal_bond_channels (spin-orbital, bond archive of the final Green "
              "function)")
POLE_CONSEQUENCE = ("The bond archive's chi_w is built on that susceptibility; it describes no "
                    "physical state and is still written.")
MAP_POLE_WHERE = "FLEX longitudinal_bond_channels (spin-orbital, final map)"
T_NO_ARCHIVE_ITER0 = ("spin-orbital bond archive: IterationMax = 0 runs no map, so no bond "
                      "archive is written.")
OLD_ITER0 = ("save_results: no map was executed (IterationMax=0); the last-map archives chi0q, "
             "chiq_s, chiq_c, chiq and the bond archive are omitted.")
T_FORCED = ("spin-orbital bond archive: sidecar layout (the member {member} of {size} is above "
            "the small-array bound {bound}); the susceptibility is written as .npy files next "
            "to the index.")
T_COMPANION_MISSING = ("companion green file '{green}' named by '{archive}' is missing (keep the "
                       "archive, its sidecars and the green file together, or re-run FLEX)")
T_CLEANUP_FAILED = ("spin-orbital bond archive: could not remove '{name}' while cleaning up after "
                    "a failed write ({error}); remove it by hand.")
T_OLD_INDEX_REMOVED = (
    "spin-orbital bond archive '{path}': writing the new index failed after the companion "
    "green file '{green}' was replaced; the previous index was removed so that no archive "
    "names a different Green function ({error}).")
T_OLD_INDEX_KEPT = (
    "spin-orbital bond archive '{path}': writing the new index failed after the companion "
    "green file '{green}' was replaced, and the previous index could not be removed ({error}); "
    "it names a different Green function, so load_so_bond_archive refuses it; remove it or "
    "re-run FLEX.")
FP_GREEN_MISMATCH = (
    "'{path}' (written by {label}): fp_green mismatch: the Green function this run loaded has "
    "digest {run}..., the archive records {file}... (both at fingerprint_version 2); the Green "
    "function and the susceptibility do not come from the same run. Regenerate both with "
    "{label}.")
T_NONSO_LOADER = (
    "bond archive {path}: this is a spin-orbital bond archive (written by the spin-orbital "
    "FLEX with bond channels); the non-spin-orbital bond consumers cannot read it, and no "
    "consumer of it exists in this version.")
T_ONSITE_BOND = (
    "[mode] enable_spin_orbital = true (hwave_sc, issue #83): '{path}' was written by the "
    "spin-orbital FLEX with bond channels ([mode.param] longitudinal_bond_channels = true): "
    "the chi0q, chiq and green of that run describe the bond-resolved run, which [eliashberg] "
    "chi0q_mode = \"flex\" cannot use. Re-run FLEX with longitudinal_bond_channels = false "
    "for this mode.")
EQUAL_PATHS = ("[file.output] the artifacts '{}' and '{}' resolve to the same file {}; give "
               "them distinct file names (numpy appends '.npz' to archive names without that "
               "suffix).")

#: the [file.output] keys of every run here
OUTPUTS = {"green": "green", "sigma": "sigma", "chi0q": "chi0q", "chiq": "chiq"}
#: the SCF settings of the shared runs: two maps, Mix 0.5 (the final G is
#: not the last map's)
BASE_PARAM = {"IterationMax": 2, "Mix": 0.5, "Longitudinal_Bond_Output_Full": True}
DENSE = {"SO_Bond_Layout": "Dense"}
SIDECAR = {"Longitudinal_Bond_Output_Layout": "SideCar"}


# --------------------------------------------------------------------------
# running and capturing
# --------------------------------------------------------------------------

class Run(NamedTuple):
    """One solve + save_results: the solver, its green_info, the output
    directory, the log records and the writer captures (one per call)."""
    solver: object
    green_info: dict
    out: str
    records: tuple
    written: tuple


def _copy_arrays(a):
    from hwave.solver import so_bond
    if isinstance(a, so_bond.PackedArrays):
        return so_bond.PackedArrays(a.sector_map, tuple(np.array(b, copy=True) for b in a.blocks))
    return np.array(a, copy=True)


@contextlib.contextmanager
def writer_spy():
    """Patch ``so_bond_archive.write_so_bond_archive``: every call records
    its keyword arguments, with copies of ``chi_w`` and ``gamma_b``
    (``chi_w_obj`` / ``gamma_b_obj`` keep the objects themselves)."""
    from hwave.solver import so_bond_archive as sa
    calls = []
    real = sa.write_so_bond_archive

    def spy(path, **kw):
        rec = dict(kw, path=path, chi_w_obj=kw["chi_w"], gamma_b_obj=kw["gamma_b"])
        rec["chi_w"] = _copy_arrays(kw["chi_w"])
        rec["gamma_b"] = _copy_arrays(kw["gamma_b"])
        calls.append(rec)
        return real(path, **kw)
    with mock.patch.object(sa, "write_so_bond_archive", spy):
        yield calls


class _Records(logging.Handler):
    def __init__(self):
        super().__init__(level=logging.INFO)
        self.records = []

    def emit(self, record):
        self.records.append(record)


def run_case(d, fixture, param=None, *, output=None, solver_hook=None):
    """Write the FLEX case of ``fixture`` into ``d`` (``[file.output]`` =
    :data:`OUTPUTS` updated by ``output``), validate the paths, solve and
    save as ``hwave.qlms.run`` does; returns a :class:`Run`."""
    out_map = dict(OUTPUTS)
    out_map.update(output or {})
    os.makedirs(d, exist_ok=True)
    inp = X.write_flex_case(d, fixture, param=param, output=out_map)
    info = inp["file"]["output"]
    out = info["path_to_output"]
    os.makedirs(out, exist_ok=True)
    solver, green_info = X.build_solver(inp)
    if solver_hook is not None:
        solver_hook(solver)
    root = logging.getLogger()
    handler = _Records()
    level = root.level
    root.addHandler(handler)
    root.setLevel(logging.INFO)
    try:
        with writer_spy() as calls:
            solver.validate_output_paths(info, out)
            solver.solve(green_info, out)
            solver.save_results(info, green_info)
    finally:
        root.removeHandler(handler)
        root.setLevel(level)
    return Run(solver, green_info, out, tuple(handler.records), tuple(calls))


def messages(records, level=None):
    return [r.getMessage() for r in records if level is None or r.levelname == level]


def load(path, solver, **kw):
    from hwave.solver import so_bond_archive as sa
    return sa.load_so_bond_archive(path, **X.run_inputs(solver), **kw)


def dense(a):
    from hwave.solver import so_bond
    return a.to_dense() if isinstance(a, so_bond.PackedArrays) else np.asarray(a)


def rel(a, b):
    return float(np.abs(a - b).max()) / float(np.abs(b).max())


# --------------------------------------------------------------------------
# the shared real archives
# --------------------------------------------------------------------------

#: name -> (fixture, extra [mode.param]); the archives of setUpModule
BASES = {
    "packed_npz": (X.ZEEMAN_SZ_ALL, {}),
    "packed_sidecar": (X.ZEEMAN_SZ_ALL, SIDECAR),
    "dense_npz": (X.ZEEMAN_SZ_ALL, DENSE),
    "dense_sidecar": (X.ZEEMAN_SZ_ALL, dict(DENSE, **SIDECAR)),
    "soc_npz": (F.SOC_ALL_TYPES, {}),
    "soc_sidecar": (F.SOC_ALL_TYPES, SIDECAR),
}
#: name -> (block layout, output layout)
LAYOUTS = {"packed_npz": ("sz_packed", "npz"), "packed_sidecar": ("sz_packed", "sidecar"),
           "dense_npz": ("dense", "npz"), "dense_sidecar": ("dense", "sidecar"),
           "soc_npz": ("dense", "npz"), "soc_sidecar": ("dense", "sidecar")}
_TMP = []
RUNS = {}


def setUpModule():
    _TMP.append(tempfile.TemporaryDirectory(prefix="so_bond_archive_flex_"))
    for name, (fx, extra) in BASES.items():
        d = os.path.join(_TMP[0].name, name)
        os.makedirs(d)
        RUNS[name] = run_case(d, fx, dict(BASE_PARAM, **extra))


def tearDownModule():
    RUNS.clear()
    while _TMP:
        _TMP.pop().cleanup()


def archive_path(run):
    return os.path.join(run.out, "longitudinal_bond.npz")


class _Copies(unittest.TestCase):
    """A per-test directory; :meth:`copy` copies a shared base directory
    into it and returns the archive path of the copy."""

    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="so_bond_archive_flex_row_")
        self.addCleanup(shutil.rmtree, self.dir, True)

    def copy(self, base="dense_npz"):
        """``<dir>/c<n>/out/longitudinal_bond.npz``: a copy of the output
        directory of the shared run ``base`` (a run of :func:`run_case`
        into ``<dir>/c<n>`` writes into the same directory)."""
        d = os.path.join(self.dir, "c{}".format(len(os.listdir(self.dir))), "out")
        shutil.copytree(RUNS[base].out, d)
        return os.path.join(d, "longitudinal_bond.npz")

    def refused(self, path, text, base="dense_npz", **kw):
        """The reader refuses ``path`` with ``text`` and reads no slot
        payload (no slot block of the npz, no memory map)."""
        with _Spies() as spy:
            with self.assertRaises(ValueError) as cm:
                load(path, RUNS[base].solver, **kw)
        self.assertEqual(str(cm.exception), text)
        self.assertEqual(spy.slot_reads(), [])
        self.assertEqual(spy.mmap, [])


# --------------------------------------------------------------------------
# round trip
# --------------------------------------------------------------------------

class TestRoundTrip(unittest.TestCase):
    def test_accepts_and_reads(self):
        from hwave.solver import so_bond
        for name, (bl, ol) in LAYOUTS.items():
            with self.subTest(base=name):
                run = RUNS[name]
                self.assertEqual(len(run.written), 1)
                w = run.written[0]
                a = load(archive_path(run), run.solver)
                self.assertEqual((a.layout, a.bond_archive_layout), (bl, ol))
                self.assertEqual((w["block_layout"], w["output_layout"]), (bl, ol))
                self.assertEqual((a.B, a.nd, a.nmat, a.nvol), (5, 2, X.FLEX_NMAT, 16))
                self.assertEqual(a.green_path, os.path.join(run.out, "green.npz"))
                for slot, want in (("chi_w", w["chi_w"]), ("Gamma_bond", w["gamma_b"])):
                    got = a.read_slot(slot)
                    if bl == "sz_packed":
                        self.assertIsInstance(got, so_bond.PackedArrays)
                        for g, e in zip(got.blocks, want.blocks):
                            np.testing.assert_array_equal(g, e)
                    else:
                        np.testing.assert_array_equal(got, want)
                self.assertFalse(a.has_slot("chibar_w"))
                # two maps do not converge: the archive records the run's flag
                self.assertIs(run.solver.scf_converged, False)
                self.assertIs(bool(_members(archive_path(run))["scf_converged"]), False)
                sidecars = sorted(glob.glob(os.path.join(run.out, "*.npy")))
                self.assertEqual(len(sidecars), {"npz": 0, "sidecar": 1 if bl == "dense"
                                                 else 3}[ol], sidecars)

    def test_packed_equals_dense(self):
        """The packed archive read dense equals the dense run's chi_w; the
        dense run's off-sector entries are exactly zero."""
        for ol in ("npz", "sidecar"):
            with self.subTest(output_layout=ol):
                p = load(archive_path(RUNS["packed_" + ol]), RUNS["packed_" + ol].solver)
                dn = load(archive_path(RUNS["dense_" + ol]), RUNS["dense_" + ol].solver)
                pd, dd = p.read_slot("chi_w", want="dense"), dn.read_slot("chi_w")
                dev = rel(pd, dd)
                self.assertTrue(np.isfinite(dev))
                self.assertLessEqual(dev, TOL)
                sm = p.sector_map
                sector = np.full(dd.shape[-1], -1)
                for s in range(3):
                    sector[np.asarray(sm.index)[sm.offsets[s]:sm.offsets[s + 1]]] = s
                self.assertTrue(np.all(sector >= 0))
                off = sector[:, None] != sector[None, :]
                self.assertTrue(off.any())
                self.assertEqual(float(np.abs(dd[..., off]).max()), 0.0)
                self.assertGreater(float(np.abs(dd[..., ~off]).max()), 0.1)


# --------------------------------------------------------------------------
# content
# --------------------------------------------------------------------------

def oracle_chi(run, tail=None):
    """chi_B of the companion green.npz of ``run``: the production bubble
    of ``G - tail_w`` (``tail_w = coeff / (i w_n)``, ``tail`` the paired
    ``(coeff, beta / 2)`` constant) in a separate dense store, the vertex
    rebuilt test-side, ``np.linalg.solve(I + chibar Gamma, chibar)`` per
    ``(l, q)``. ``tail`` is ``None`` or ``coeff``."""
    from hwave.solver import flex_bond, so_bond
    s = run.solver
    with np.load(os.path.join(run.out, "green.npz"), allow_pickle=False) as z:
        G = np.array(z["green"])
    compiled = s.ham_info.so_context.compiled
    topo = so_bond.so_bond_topology(compiled, np.eye(3))
    shape = tuple(int(x) for x in s.lattice.shape)
    nmat, nvol, nd = G.shape[1], G.shape[2], G.shape[-1]
    beta = 1.0 / s.T
    plan = so_bond.plan_so_bond(nmat=nmat, nvol=nvol, B=len(topo.delta_r), nd=nd,
                                consumers=("chi_dynamic",),
                                caps=so_bond.SOBondCaps(None, None, None), layout="dense")
    gamma = so_bond.so_bond_vertex(compiled, topo, shape, plan=plan)
    g0_tail = None
    if tail is not None:
        iw = (np.arange(nmat) * 2 + 1 - nmat) * np.pi / beta
        G = G - (tail / (1j * iw))[None, :, None, None, None] * np.eye(nd)
        g0_tail = np.broadcast_to(tail * beta / 2 * np.eye(nd), G.shape).astype(np.complex128)
    ND = plan.ND
    with flex_bond.BondBlockStore(nmat, nvol, ND, nd * nd, ("chibar",)) as st:
        so_bond.so_bond_bubble(np.ascontiguousarray(G), g0_tail, beta, topo, shape, plan=plan,
                               freq="all", sink=st)
        chibar = st.detach("chibar")
    return np.linalg.solve(np.eye(ND) + np.matmul(chibar, gamma[None]), chibar)


class TestContent(_Copies):
    def test_oracle(self):
        for name in ("dense_npz", "packed_sidecar", "soc_npz"):
            with self.subTest(base=name):
                run = RUNS[name]
                chi = oracle_chi(run)
                a = load(archive_path(run), run.solver)
                dev = rel(a.read_slot("chi_w", want="dense"), chi)
                print("archive chi_w vs oracle %-15s %.3e" % (name, dev))
                self.assertTrue(np.isfinite(dev))
                self.assertLessEqual(dev, TOL)

    def test_tail(self):
        """coeff_tail = 1: the archive is chi_B of G - tail with its paired
        tail (the tail of tests/test_so_bond_archive_pass.py::test_tail)."""
        run = run_case(self.dir, X.ZEEMAN_SZ_ALL, dict(BASE_PARAM, coeff_tail=1.0))
        a = load(archive_path(run), run.solver)
        got = a.read_slot("chi_w", want="dense")
        chi, no_tail = oracle_chi(run, tail=1.0), oracle_chi(run)
        dev, effect = rel(got, chi), rel(no_tail, chi)
        print("archive chi_w tail (dev, tail effect): %.3e %.3e" % (dev, effect))
        self.assertTrue(np.isfinite(dev) and np.isfinite(effect))
        self.assertLessEqual(dev, TOL)
        self.assertGreaterEqual(effect, MISS)

    def test_differs_from_chiq(self):
        """chi_w is built from the final G, chiq.npz from the last map's
        (Mix 0.5): their channel-0 blocks differ."""
        for name in ("dense_npz", "soc_npz"):
            with self.subTest(base=name):
                run = RUNS[name]
                chi = load(archive_path(run), run.solver).read_slot("chi_w")
                P = 4
                with np.load(os.path.join(run.out, "chiq.npz"), allow_pickle=False) as z:
                    chiq = np.array(z["chiq"]).reshape(chi.shape[0], chi.shape[1], P, P)
                d = rel(chi[..., :P, :P], chiq)
                print("archive chi_w channel 0 vs chiq.npz %-10s %.3e" % (name, d))
                self.assertTrue(np.isfinite(d))
                self.assertGreaterEqual(d, NONVAC)

    def test_break_row(self):
        """The dressing of the pass replaced by the identity (inside a
        wrapper of the pass only): the archive misses the oracle, the maps
        dress as often as without the break."""
        from hwave.solver import flex_so_bond, so_bond
        real_pass, real_dress = flex_so_bond.so_bond_archive_pass, so_bond.so_bond_dress
        inside = []

        def counting(counts, identity):
            def dress(X_, *a, **kw):
                counts["pass" if inside else "maps"] += 1
                return X_ if (inside and identity) else real_dress(X_, *a, **kw)
            return dress

        def wrapped(*a, **kw):
            inside.append(1)
            try:
                return real_pass(*a, **kw)
            finally:
                inside.pop()

        runs, counts = {}, {}
        for name, identity in (("ref", False), ("broken", True)):
            counts[name] = {"pass": 0, "maps": 0}
            with mock.patch.object(flex_so_bond, "so_bond_archive_pass", wrapped), \
                    mock.patch.object(so_bond, "so_bond_dress",
                                      counting(counts[name], identity)):
                runs[name] = run_case(os.path.join(self.dir, name), X.ZEEMAN_SZ_ALL,
                                      dict(BASE_PARAM, **DENSE))
            if name == "ref":
                chi = oracle_chi(runs["ref"])
        broken = runs["broken"]
        self.assertGreater(counts["broken"]["pass"], 0)
        self.assertEqual(counts["broken"], counts["ref"])
        self.assertLessEqual(rel(load(archive_path(runs["ref"]), runs["ref"].solver)
                                 .read_slot("chi_w"), chi), TOL)
        miss = rel(load(archive_path(broken), broken.solver).read_slot("chi_w"), chi)
        print("archive break row miss %.3e" % miss)
        self.assertTrue(np.isfinite(miss))
        self.assertGreaterEqual(miss, MISS)


# --------------------------------------------------------------------------
# companion
# --------------------------------------------------------------------------

class TestCompanion(_Copies):
    def test_archive_fp_green_equals_green_npz(self):
        from hwave.solver import fingerprint
        for name in LAYOUTS:
            with self.subTest(base=name):
                run = RUNS[name]
                g = _members(os.path.join(run.out, "green.npz"))
                arch = _members(archive_path(run))
                digest = fingerprint.green_digest(g["green"], fingerprint.GREEN_LAYOUT)
                self.assertEqual(str(arch["fp_green"]), str(g["fp_green"]))
                self.assertEqual(str(g["fp_green"]), digest)
                self.assertEqual(str(arch["green_file"]), "green.npz")
                self.assertEqual(str(arch["producer_stage"]), "flex_so_bond")

    def test_refusals(self):
        # missing
        path = self.copy()
        os.remove(os.path.join(os.path.dirname(path), "green.npz"))
        self.refused(path, T_COMPANION_MISSING.format(
            green=os.path.join(os.path.dirname(path), "green.npz"), archive=path))
        # a companion of another G: the consistent green.npz of a run of the
        # same case stopped after one map (a genuinely different G)
        other_run = run_case(os.path.join(self.dir, "other"), X.ZEEMAN_SZ_ALL,
                             dict(BASE_PARAM, IterationMax=1, **DENSE))
        other = os.path.join(other_run.out, "green.npz")
        path = self.copy()
        self.assertGreaterEqual(rel(_members(other)["green"],
                                    _members(os.path.join(os.path.dirname(path),
                                                          "green.npz"))["green"]), NONVAC)
        shutil.copy(other, os.path.join(os.path.dirname(path), "green.npz"))
        self.assertNotEqual(str(_members(other)["fp_green"]),
                            str(_members(path)["fp_green"]))
        with self.assertRaises(ValueError) as cm:
            load(path, RUNS["dense_npz"].solver)
        self.assertEqual(str(cm.exception), FP_GREEN_MISMATCH.format(
            path=path, label=LABEL, run=str(_members(other)["fp_green"])[:12],
            file=str(_members(path)["fp_green"])[:12]))
        # green_file absolute / with ..
        for g in (os.path.join(self.dir, "green.npz"), "../green.npz", "sub/../../green.npz"):
            with self.subTest(green_file=g):
                path = self.copy()
                _rewrite(path, dict(green_file=np.str_(g)))
                self.refused(path, _M(path, "green_file is not a relative file path",
                                      _bounded(repr(g)), "a file name or a relative path "
                                      "below the archive directory", LABEL))

    def test_suffixless_names(self):
        run = run_case(self.dir, X.ZEEMAN_SZ_ALL, BASE_PARAM,
                       output={"green": "g", "longitudinal_bond": "lb"})
        self.assertTrue(os.path.isfile(os.path.join(run.out, "g.npz")))
        self.assertFalse(os.path.exists(os.path.join(run.out, "longitudinal_bond.npz")))
        a = load(os.path.join(run.out, "lb.npz"), run.solver)
        self.assertEqual(a.green_path, os.path.join(run.out, "g.npz"))
        self.assertEqual(str(_members(os.path.join(run.out, "lb.npz"))["green_file"]), "g.npz")

    def test_green_below_the_archive(self):
        run = run_case(self.dir, X.ZEEMAN_SZ_ALL, BASE_PARAM,
                       output={"green": "sub/g", "longitudinal_bond": "lb"},
                       solver_hook=lambda s: os.makedirs(os.path.join(self.dir, "out", "sub")))
        self.assertEqual(str(_members(os.path.join(run.out, "lb.npz"))["green_file"]),
                         os.path.join("sub", "g.npz"))
        load(os.path.join(run.out, "lb.npz"), run.solver)


# --------------------------------------------------------------------------
# sidecar / npz
# --------------------------------------------------------------------------

class TestSidecar(_Copies):
    def test_forced_sidecar(self):
        """No layout requested and a member above the small-array bound:
        the sidecar layout, one INFO."""
        from hwave.solver import so_bond
        with mock.patch.object(so_bond, "small_array_bound", return_value=float(2 ** 17)):
            run = run_case(self.dir, X.ZEEMAN_SZ_ALL, BASE_PARAM)
        info = T_FORCED.format(member="chi_w@0", size="0.39 MiB", bound="0.12 MiB")
        self.assertEqual(messages(run.records).count(info), 1)
        a = load(archive_path(run), run.solver)
        self.assertEqual(a.bond_archive_layout, "sidecar")
        self.assertEqual(len(glob.glob(os.path.join(run.out, "chi_w@*.npy"))), 3)
        self.assertFalse([m for m in messages(RUNS["packed_npz"].records)
                          if m.startswith("spin-orbital bond archive: sidecar layout")])

    def _sidecar(self, path, block="chi_w"):
        nonce = str(_members(path)["archive_nonce"])
        return os.path.join(os.path.dirname(path), "{}.{}.npy".format(block, nonce))

    def test_refusals(self):
        base = "dense_sidecar"
        # one byte longer / shorter
        for delta in (1, -1):
            with self.subTest(size=delta):
                path = self.copy(base)
                full = self._sidecar(path)
                size = os.path.getsize(full)
                with open(full, "r+b") as fh:
                    if delta > 0:
                        fh.seek(0, 2)
                        fh.write(b"\0")
                    else:
                        fh.truncate(size - 1)
                self.refused(path, _M(path, "sidecar member chi_w has the wrong size",
                                      "{} bytes, {} recorded".format(size + delta, size),
                                      "{} bytes".format(size), LABEL), base=base)
        chi = np.load(self._sidecar(self.copy(base)))
        # wrong dtype, wrong slot shape (the recorded size updated)
        for row, arr, check, stored, expected in (
                ("dtype", chi.astype(np.complex64), "payload member chi_w has the wrong dtype",
                 "complex64", "complex128"),
                ("shape", np.ascontiguousarray(chi.reshape(chi.shape[:2] + (10, 40))),
                 "payload member chi_w has the wrong slot shape", "(10, 40)", "(20, 20)")):
            with self.subTest(row=row):
                path = self.copy(base)
                full = self._sidecar(path)
                np.save(full, arr)
                _rewrite(path, dict(chi_w_nbytes=np.int64(os.path.getsize(full))))
                self.refused(path, _M(path, check, stored, expected, LABEL), base=base)
        # a block both inline and as a sidecar
        path = self.copy(base)
        _rewrite(path, dict(chi_w=chi))
        self.refused(path, _M(path, "payload member chi_w is stored both inline and as a "
                              "sidecar", "both", "one of them", LABEL), base=base)

    def test_other_archive_sidecar(self):
        """Two archives in one directory; A's chi_w_file pointed at B's
        sidecar (right size and header) is refused by the nonce rule."""
        out = os.path.join(self.dir, "two")
        a = run_case(out, X.ZEEMAN_SZ_ALL, dict(BASE_PARAM, **DENSE, **SIDECAR),
                     output={"longitudinal_bond": "a", "green": "ga"})
        b = run_case(out, X.ZEEMAN_SZ_ALL, dict(BASE_PARAM, **DENSE, **SIDECAR,
                                                IterationMax=1),
                     output={"longitudinal_bond": "b", "green": "gb"})
        pa, pb = os.path.join(a.out, "a.npz"), os.path.join(b.out, "b.npz")
        load(pa, a.solver)
        nb, na = str(_members(pb)["chi_w_file"]), str(_members(pa)["chi_w_file"])
        self.assertTrue(os.path.isfile(os.path.join(out, "out", nb)))
        _rewrite(pa, dict(chi_w_file=np.str_(nb)))
        with _Spies() as spy:
            with self.assertRaises(ValueError) as cm:
                load(pa, a.solver)
        self.assertEqual(str(cm.exception), _M(
            pa, "sidecar file name of chi_w does not carry the archive_nonce", repr(nb),
            repr(na), LABEL))
        self.assertEqual((spy.slot_reads(), spy.mmap), ([], []))

    def test_rerun_removes_old_sidecars(self):
        path = self.copy("packed_sidecar")
        d = os.path.dirname(path)
        old = sorted(os.path.basename(p) for p in glob.glob(os.path.join(d, "*.npy")))
        self.assertEqual(len(old), 3)
        run = run_case(os.path.dirname(d), X.ZEEMAN_SZ_ALL,
                       dict(BASE_PARAM, **SIDECAR, IterationMax=1))
        self.assertEqual(run.out, d)
        new = sorted(os.path.basename(p) for p in glob.glob(os.path.join(d, "*.npy")))
        self.assertEqual(len(new), 3)
        self.assertFalse(set(old) & set(new))
        info = [m for m in messages(run.records)
                if m.startswith("spin-orbital bond archive '{}': removed 3".format(path))]
        self.assertEqual(len(info), 1, messages(run.records))
        load(path, run.solver)


# --------------------------------------------------------------------------
# structure and compatibility
# --------------------------------------------------------------------------

class TestStructure(_Copies):
    """Step 1 of the reader on copies of the real dense npz archive (B = 5,
    nd = 2, Nmat = 16, beta = 4, a 4 x 4 x 1 mesh)."""

    def test_bases_accepted(self):
        for name in LAYOUTS:
            with self.subTest(base=name):
                load(self.copy(name), RUNS[name].solver)

    def test_rows(self):
        from hwave.solver import so_bond_archive as sa
        path = self.copy()
        np.savez(path, bond_archive_schema=np.int64(3), S_bond=np.zeros((2, 2)))
        self.refused(path, sa.T_NOT_SO_ARCHIVE.format(path=path, schema="3", so="absent"))
        rows = [
            ("schema 3", dict(bond_archive_schema=np.int64(3)), None),
            ("nmat odd", dict(nmat=np.int64(15)),
             ("member nmat is out of range", "15", "> 0 and even")),
            ("nd", dict(nd=np.int64(4)), ("nd differs from 2 * norb_phys", "4", "2")),
            ("B", dict(B=np.int64(4)), ("B differs from len(delta_r)", "4", "5")),
            ("nmat", dict(nmat=np.int64(8)), ("nmat differs from state_nmat", "8", "16")),
            ("beta", dict(beta=np.float64(2.0)),
             ("beta differs from state_beta", "2.0", "4.0")),
            ("cell_shape", dict(cell_shape=np.array([2, 8, 1])),
             ("cell_shape differs from prov_cell_shape", "[2, 8, 1]", "[4, 4, 1]")),
            ("index_order", dict(index_order=np.str_("pair_major")),
             ("index_order is not this version's", "'pair_major'",
              "'bond_major_pair_code_slot'")),
            ("freq_axis", dict(freq_axis=np.str_("ir")),
             ("freq_axis is not this version's", "'ir'", "'uniform_boson_centered'")),
            ("green_layout", dict(green_layout=np.str_("b_k_f")),
             ("green_layout is not the canonical layout", "'b_k_f'", "'b_f_k_o_o_v1'")),
            ("scf_converged", dict(scf_converged=_DROP),
             ("member scf_converged is missing", "absent", "present")),
            ("sidecar layout of an npz archive", dict(bond_archive_layout=np.str_("sidecar")),
             ("payload member chi_w does not follow bond_archive_layout sidecar", "inline",
              "sidecar")),
        ]
        for name, ch, row in rows:
            with self.subTest(row=name):
                path = self.copy()
                _rewrite(path, ch)
                text = (sa.T_NOT_SO_ARCHIVE.format(path=path, schema="3", so="True")
                        if row is None else _M(path, *row, LABEL))
                self.refused(path, text)


class TestCompatibility(_Copies):
    """Step 2: the run differs from the archive's (and its companion's)
    provenance. The companion is checked first, so the refusal names it."""

    def _refused(self, base, run, text_of):
        path = self.copy(base)
        green = os.path.join(os.path.dirname(path), "green.npz")
        from hwave.solver import so_bond_archive as sa
        with self.assertRaises(ValueError) as cm:
            sa.load_so_bond_archive(path, **run)
        self.assertEqual(str(cm.exception), text_of(green))

    def test_rows(self):
        base = "packed_npz"
        s = RUNS[base].solver
        good = X.run_inputs(s)
        point = "the identical state point"
        st = good["run_state"]
        meta = good["run_band"].meta
        from hwave.solver import fingerprint
        band = fingerprint.band_inputs(good["run_band"].transfer, good["run_band"].extern,
                                       meta._replace(coeff_extern=meta.coeff_extern + 0.25))
        rows = (
            ("Extern", dict(good, run_band=band), lambda g: (
                "'{}' (written by {}) does not match this run: band differs: coeff_extern: "
                "file {!r}, run {!r}. Archives are reused only at {} (not a warm start across "
                "a sweep); re-run {} with this input (the same Transfer/Extern/coeff_extern)."
                .format(g, LABEL, meta.coeff_extern, meta.coeff_extern + 0.25, point, LABEL))),
            ("beta", dict(good, run_state=st._replace(beta=2.0)), lambda g: (
                "'{}' (written by {}) does not match this run: state differs: beta {!r} (file) "
                "vs 2.0 (run). Archives are reused only at {}; re-run {} at this filling and "
                "temperature.".format(g, LABEL, st.beta, point, LABEL))),
            ("filling", dict(good, run_state=st._replace(mu=st.mu + 0.5)), lambda g: (
                "'{}' (written by {}) does not match this run: state differs: fixed mu {!r} "
                "(file) vs {!r} (run) energy units. Archives are reused only at {}; re-run {} "
                "at this filling and temperature.".format(g, LABEL, st.mu, st.mu + 0.5, point,
                                                          LABEL))),
        )
        for name, run, text_of in rows:
            with self.subTest(row=name):
                self._refused(base, run, text_of)

    def test_topology(self):
        """Equal channels with another max_shells are accepted; other
        channels are refused."""
        base = "packed_npz"
        good = X.run_inputs(RUNS[base].solver)
        bond = good["run_bond"]
        from hwave.solver import so_bond_archive as sa
        sa.load_so_bond_archive(self.copy(base), **dict(
            good, run_bond=bond._replace(max_shells=int(bond.max_shells) + 3)))
        dr, rv = np.asarray(bond.delta_r), np.asarray(bond.reverse)
        perm = np.array([0, 2, 1, 4, 3]) if len(dr) == 5 else None
        self.assertIsNotNone(perm)
        other = bond._replace(delta_r=dr[perm], reverse=np.argsort(perm)[rv[perm]])
        path = self.copy(base)
        green = os.path.join(os.path.dirname(path), "green.npz")
        with self.assertRaises(ValueError) as cm:
            sa.load_so_bond_archive(path, **dict(good, run_bond=other))
        self.assertEqual(str(cm.exception), (
            "'{}' (written by {}) does not match this run: resolved bond topology differs: "
            "file delta_r {} reverse {} (max_shells {}), run delta_r {} reverse {} (max_shells "
            "{}). Archives are reused only at the identical state point (not a warm start "
            "across a sweep); re-run {} with this input (the same bond settings).").format(
                green, LABEL, dr.tolist(), rv.tolist(), int(bond.max_shells),
                np.asarray(other.delta_r).tolist(), np.asarray(other.reverse).tolist(),
                int(bond.max_shells), LABEL))


# --------------------------------------------------------------------------
# cross-refusals on the real outputs (the on-site hwave_sc rows
# of B = 1 and B > 1, static and dynamic, are TestOnsiteRefusesBond of
# tests/test_flex_so_bond_entry.py)
# --------------------------------------------------------------------------

class TestCrossRefusals(_Copies):
    def test_nonso_loader(self):
        from hwave.solver import eliashberg_bond_io
        for name in ("packed_npz", "dense_sidecar"):
            with self.subTest(base=name):
                path = self.copy(name)
                with self.assertRaises(ValueError) as cm:
                    eliashberg_bond_io.load_bond_archive(
                        path, norb=1, nmat_expected=X.FLEX_NMAT,
                        cell_shape_expected=(4, 4, 1), beta_expected=1.0 / X.FLEX_T)
                self.assertEqual(str(cm.exception), T_NONSO_LOADER.format(path=path))

    def test_onsite_dynamic_on_an_archive_run(self):
        """B > 1, frequency = "dynamic": the on-site consumer refuses the
        outputs of a run that wrote the archive, naming chiq.npz."""
        import hwave.sc as sc
        from tests.test_flex_so_bond_entry import _sc_input
        out = os.path.dirname(self.copy("packed_npz"))
        with tempfile.TemporaryDirectory() as d2:
            inp = _sc_input(d2, X.ZEEMAN_EXTERN_ONSITE, out,
                            {"Bond_Channels": False, "Chi0q_Mode": "flex",
                             "Frequency": "dynamic"})
            with self.assertRaises(ValueError) as cm:
                sc.calc_eliashberg(inp)
        self.assertEqual(str(cm.exception),
                         T_ONSITE_BOND.format(path=os.path.join(out, "chiq.npz")))

    def test_capability_backstop(self):
        """The capability check behind the on-site peek: a copy of the
        archive run's chiq checked as an on-site artifact."""
        from hwave.solver import fingerprint
        path = os.path.join(os.path.dirname(self.copy("packed_npz")), "chiq.npz")
        run = X.run_inputs(RUNS["packed_npz"].solver)
        with np.load(path, allow_pickle=False) as z:
            with self.assertRaises(fingerprint.SOArchiveError) as cm:
                fingerprint.check_so_archive(
                    z, path, kind="chiq", run_band=run["run_band"], run_state=run["run_state"],
                    nmat_rule="exact", run_interaction=run["run_interaction"],
                    need="so_onsite")
        self.assertEqual(str(cm.exception), "{} was produced by {} (so_bond), but this run "
                         "needs a so_onsite artifact".format(path, LABEL))


# --------------------------------------------------------------------------
# the wiring
# --------------------------------------------------------------------------

class TestWiring(_Copies):
    def test_heads_equal(self):
        from hwave.solver import so_bond_archive as sa
        from hwave.solver.flex import FLEX
        self.assertEqual(FLEX._SO_BOND_HEAD, sa.SO_BOND_FLEX_HEAD)
        self.assertEqual(FLEX._SO_BOND_HEAD, H)

    def test_no_store_copy(self):
        """The writer receives the store's own chibar slot (detached, no
        copy) and the vertex of the solve."""
        from hwave.solver import flex_so_bond
        real = flex_so_bond.so_bond_archive_pass
        seen = {}

        def wrapped(store, gamma_b, *a, **kw):
            res = real(store, gamma_b, *a, **kw)
            seen["slot"] = store._arr("chibar")
            seen["gamma_b"] = gamma_b
            return res
        for name, extra in (("packed", {}), ("dense", DENSE)):
            with self.subTest(layout=name):
                with mock.patch.object(flex_so_bond, "so_bond_archive_pass", wrapped):
                    run = run_case(os.path.join(self.dir, name), X.ZEEMAN_SZ_ALL,
                                   dict(BASE_PARAM, **extra))
                w = run.written[0]
                self.assertIs(w["gamma_b_obj"], seen["gamma_b"])
                if name == "dense":
                    self.assertIs(w["chi_w_obj"], seen["slot"])
                else:
                    self.assertEqual(len(w["chi_w_obj"].blocks), 3)
                    for got, slot in zip(w["chi_w_obj"].blocks, seen["slot"]):
                        self.assertIs(got, slot)
                self.assertNotIn("so_bond_archive", run.green_info)

    def test_pole_check_once(self):
        from hwave.solver import bond_channels, flex_so_bond
        real_check, real_pass = bond_channels.check_rpa_pole, flex_so_bond.so_bond_archive_pass
        calls, passes = [], []

        def check(*a, **kw):
            calls.append((a, kw))
            return real_check(*a, **kw)

        def wrapped(*a, **kw):
            res = real_pass(*a, **kw)
            passes.append(res)
            return res
        with mock.patch.object(bond_channels, "check_rpa_pole", check), \
                mock.patch.object(flex_so_bond, "so_bond_archive_pass", wrapped):
            run = run_case(self.dir, X.ZEEMAN_SZ_ALL, BASE_PARAM)
        self.assertEqual(len(passes), 1)
        self.assertEqual([kw["where"] for _, kw in calls], [MAP_POLE_WHERE, POLE_WHERE])
        a, kw = calls[1]
        self.assertEqual(kw["consequence"], POLE_CONSEQUENCE)
        self.assertIs(kw["info"], True)
        self.assertEqual(a[0], {"spinful": passes[0].pole[0]})
        self.assertEqual(messages(run.records).count(T_PASS), 1)
        arch = _members(archive_path(run))
        self.assertFalse([k for k in arch if k.startswith("rpa_")])

    def test_pole_members_past_pole(self):
        from hwave.solver import flex_so_bond, so_bond
        real_pass = flex_so_bond.so_bond_archive_pass
        q = (0.25, 0.5, 0.0)

        def wrapped(*a, **kw):
            with mock.patch.object(so_bond, "static_pole_indicator",
                                   return_value=((1.2, q, 1.2), None, None)):
                return real_pass(*a, **kw)
        with mock.patch.object(flex_so_bond, "so_bond_archive_pass", wrapped):
            run = run_case(self.dir, X.ZEEMAN_SZ_ALL, BASE_PARAM)
        warns = [m for m in messages(run.records, "WARNING") if "past its pole" in m]
        self.assertEqual(len(warns), 1, warns)
        self.assertTrue(warns[0].startswith(POLE_WHERE + ": "), warns[0])
        self.assertIn(POLE_CONSEQUENCE, warns[0])
        arch = _members(archive_path(run))
        self.assertIs(bool(arch["rpa_past_pole"]), True)
        self.assertEqual(float(arch["rpa_pole_indicator_spinful"]), 1.2)
        np.testing.assert_array_equal(arch["rpa_pole_q_spinful"], np.array(q))
        self.assertEqual(sorted(k for k in arch if k.startswith("rpa_")),
                         ["rpa_past_pole", "rpa_pole_indicator_spinful", "rpa_pole_q_spinful"])
        # the maps' final check is below the pole: chi0q / chiq carry none
        for n in ("chi0q", "chiq"):
            self.assertFalse([k for k in _members(os.path.join(run.out, n + ".npz"))
                              if k.startswith("rpa_")])
        load(archive_path(run), run.solver)
        # unpatched, below the pole: none
        self.assertFalse([k for k in _members(archive_path(RUNS["packed_npz"]))
                          if k.startswith("rpa_")])

    def test_pass_refusal_text(self):
        """A refusal of the pass's dressing guard propagates with the
        remedy and the note; nothing is published."""
        from hwave.solver import bond_channels, flex_so_bond
        real_pass = flex_so_bond.so_bond_archive_pass

        def strict(*a, **kw):
            kw["cond_tol"] = 0.99
            return real_pass(*a, **kw)
        with mock.patch.object(flex_so_bond, "so_bond_archive_pass", strict):
            with self.assertRaises(bond_channels.BondConditioningError) as cm:
                run_case(self.dir, X.ZEEMAN_SZ_ALL, BASE_PARAM)
        msg = str(cm.exception)
        self.assertTrue(msg.endswith(REMEDY + T_PASS_REFUSAL_NOTE), msg)
        self.assertEqual(msg.count(REMEDY), 1)
        self.assertEqual(cm.exception.cond_tol, 0.99)
        out = os.path.join(self.dir, "out")
        self.assertEqual(os.listdir(out), [])

    def test_output_info_variants(self):
        rows = ((BASE_PARAM, INFO_ARCHIVE),
                (dict(BASE_PARAM, Longitudinal_Bond_Output_Full=False), INFO_NO_ARCHIVE),
                (dict(BASE_PARAM, IterationMax=0), None),
                (dict(BASE_PARAM, IterationMax=0, Longitudinal_Bond_Output_Full=False), None))
        for i, (param, text) in enumerate(rows):
            with self.subTest(row=i):
                run = run_case(os.path.join(self.dir, str(i)), X.ZEEMAN_SZ_ALL, param)
                msgs = messages(run.records)
                for t in (INFO_ARCHIVE, INFO_NO_ARCHIVE):
                    self.assertEqual(msgs.count(t), int(t == text), t)

    def test_equal_paths_refused(self):
        s, _ = X.build_solver(X.write_flex_case(self.dir, X.ZEEMAN_SZ_ALL, param=BASE_PARAM))
        info = dict(OUTPUTS, path_to_output=self.dir, longitudinal_bond="green.npz")
        with self.assertRaises(ValueError) as cm:
            s.validate_output_paths(info, self.dir)
        self.assertEqual(str(cm.exception), EQUAL_PATHS.format(
            "green", "longitudinal_bond", os.path.join(self.dir, "green.npz")))

    def test_validation_and_write_use_one_mapping(self):
        """green.npz and the archive are written exactly to the values of
        the mapping save_results resolves (once)."""
        from hwave.solver import flex_bond
        real = flex_bond.resolve_output_paths
        calls = []

        def moved(info, out, active):
            paths = dict(real(info, out, active))
            for k, new in (("green", "g2.npz"), ("longitudinal_bond", "lb2.npz")):
                if k in paths:
                    paths[k] = os.path.join(os.path.dirname(paths[k]), new)
            calls.append(paths)
            return paths
        inp = X.write_flex_case(self.dir, X.ZEEMAN_SZ_ALL, param=BASE_PARAM, output=OUTPUTS)
        info = inp["file"]["output"]
        out = info["path_to_output"]
        os.makedirs(out)
        solver, gi = X.build_solver(inp)
        solver.validate_output_paths(info, out)
        solver.solve(gi, out)
        with mock.patch.object(flex_bond, "resolve_output_paths", moved):
            solver.save_results(info, gi)
        self.assertEqual(len(calls), 1)
        self.assertEqual(sorted(os.listdir(out)),
                         ["chi0q.npz", "chiq.npz", "g2.npz", "lb2.npz", "sigma.npz"])
        self.assertEqual(str(_members(os.path.join(out, "lb2.npz"))["green_file"]), "g2.npz")
        load(os.path.join(out, "lb2.npz"), solver)

    def test_payload_dropped_by_next_solve(self):
        """The payload of the pass lives in green_info between a solve and
        its save_results only: a failed solve and a solve without the
        archive on the same green_info drop it (the solve resets)."""
        from hwave.solver.flex import FLEX
        inp = X.write_flex_case(self.dir, X.ZEEMAN_SZ_ALL, param=BASE_PARAM, output=OUTPUTS)
        out = inp["file"]["output"]["path_to_output"]
        s1, gi = X.build_solver(inp)
        s1.solve(gi, out)
        self.assertIn("so_bond_archive", gi)                     # positive control
        with mock.patch.object(FLEX, "_so_bond_preflight", side_effect=RuntimeError("x")):
            with self.assertRaises(RuntimeError):
                s1.solve(gi, out)
        self.assertNotIn("so_bond_archive", gi)
        s1.solve(gi, out)
        self.assertIn("so_bond_archive", gi)
        os.makedirs(os.path.join(self.dir, "b"))
        inp2 = X.write_flex_case(os.path.join(self.dir, "b"), X.ZEEMAN_SZ_ALL,
                                 param=dict(BASE_PARAM, Longitudinal_Bond_Output_Full=False),
                                 output=dict(OUTPUTS, longitudinal_bond="lb"))
        info2 = inp2["file"]["output"]
        os.makedirs(info2["path_to_output"])
        s2, _ = X.build_solver(inp2)
        s2.solve(gi, info2["path_to_output"])
        self.assertNotIn("so_bond_archive", gi)
        s2.save_results(info2, gi)
        self.assertEqual(sorted(os.listdir(info2["path_to_output"])),
                         ["chi0q.npz", "chiq.npz", "green.npz", "sigma.npz"])

    def test_one_predicate(self):
        """An admitted archive without the pass's payload: the plain green
        file, no archive, no temporary file and the output notice without
        the archive."""
        from hwave.solver import fingerprint
        inp = X.write_flex_case(self.dir, X.ZEEMAN_SZ_ALL, param=BASE_PARAM, output=OUTPUTS)
        info = inp["file"]["output"]
        out = info["path_to_output"]
        os.makedirs(out)
        s, gi = X.build_solver(inp)
        s.solve(gi, out)
        self.assertIsNotNone(s._so_admission.archive)
        gi.pop("so_bond_archive")
        root, handler = logging.getLogger(), _Records()
        level = root.level
        root.addHandler(handler)
        root.setLevel(logging.INFO)
        try:
            s.save_results(info, gi)
        finally:
            root.removeHandler(handler)
            root.setLevel(level)
        msgs = messages(handler.records)
        self.assertEqual(msgs.count(INFO_NO_ARCHIVE), 1)
        self.assertNotIn(INFO_ARCHIVE, msgs)
        self.assertEqual(sorted(os.listdir(out)),
                         ["chi0q.npz", "chiq.npz", "green.npz", "sigma.npz"])
        g = _members(os.path.join(out, "green.npz"))
        self.assertEqual(str(g["fp_green"]),
                         fingerprint.green_digest(g["green"], fingerprint.GREEN_LAYOUT))

    def test_iter0_output_full(self):
        run = run_case(self.dir, X.ZEEMAN_SZ_ALL, dict(BASE_PARAM, IterationMax=0))
        msgs = messages(run.records)
        self.assertEqual(msgs.count(T_NO_ARCHIVE_ITER0), 1)
        self.assertNotIn(OLD_ITER0, msgs)
        self.assertNotIn(T_PASS, msgs)
        self.assertEqual(run.written, ())
        self.assertEqual(sorted(os.listdir(run.out)), ["green.npz", "sigma.npz"])
        for n in ("green", "sigma"):
            self.assertEqual(str(_members(os.path.join(run.out, n + ".npz"))["block_layout"]),
                             "sz_packed")
        self.assertIsNone(run.solver._so_admission.archive)
        # output_full = false keeps the old notice
        run = run_case(os.path.join(self.dir, "off"), X.ZEEMAN_SZ_ALL,
                       dict(BASE_PARAM, IterationMax=0, Longitudinal_Bond_Output_Full=False))
        self.assertEqual(messages(run.records).count(OLD_ITER0), 1)
        self.assertNotIn(T_NO_ARCHIVE_ITER0, messages(run.records))


# --------------------------------------------------------------------------
# the green commit
# --------------------------------------------------------------------------

def _tmp_greens(d):
    return sorted(n for n in os.listdir(d) if n.startswith(".green.") and n.endswith(".tmp.npz"))


class TestGreenCommit(_Copies):
    def _rerun(self, path, patches, param=None):
        """Re-run the dense archive into the directory of ``path`` (another
        G: one map) under ``patches``; returns the exception and the logs."""
        d = os.path.dirname(os.path.dirname(path))
        records = []
        handler = _Records()
        root = logging.getLogger()
        root.addHandler(handler)
        level = root.level
        root.setLevel(logging.INFO)
        try:
            with contextlib.ExitStack() as st:
                for p in patches:
                    st.enter_context(p)
                with self.assertRaises(BaseException) as cm:
                    run_case(d, X.ZEEMAN_SZ_ALL,
                             dict(BASE_PARAM, **DENSE, **(param or {"IterationMax": 1})))
        finally:
            root.removeHandler(handler)
            root.setLevel(level)
        records.extend(handler.records)
        return cm.exception, records

    @staticmethod
    def _selective(fn, target, exc):
        """``fn`` raising ``exc`` when its first (``os.remove``) or second
        (``os.replace``) path argument is ``target``."""
        real = getattr(os, fn)

        def f(*a, **kw):
            hit = a[0] if fn == "remove" else a[1]
            if os.path.abspath(hit) == os.path.abspath(target):
                raise exc
            return real(*a, **kw)
        return mock.patch.object(os, fn, f)

    def test_failed_rewrite_leaves_consistent_pair(self):
        # the index commit fails after the green commit: the old index is removed
        path = self.copy()
        d = os.path.dirname(path)
        load(path, RUNS["dense_npz"].solver)
        boom = OSError("index replace refused")
        exc, records = self._rerun(path, [self._selective("replace", path, boom)])
        self.assertIs(exc, boom)
        self.assertEqual(_tmp_greens(d), [])
        self.assertFalse(os.path.exists(path))
        green = os.path.join(d, "green.npz")
        self.assertEqual(messages(records, "WARNING").count(T_OLD_INDEX_REMOVED.format(
            path=path, green=green, error=boom)), 1)
        # the old-index removal fails too: the surviving index is refused
        path = self.copy()
        d = os.path.dirname(path)
        old_fp = str(_members(path)["fp_green"])
        rm = PermissionError("old index removal refused")
        exc, records = self._rerun(path, [self._selective("replace", path, boom),
                                          self._selective("remove", path, rm)])
        self.assertIs(exc, boom)
        self.assertEqual(_tmp_greens(d), [])
        green = os.path.join(d, "green.npz")
        self.assertEqual(messages(records, "WARNING").count(T_OLD_INDEX_KEPT.format(
            path=path, green=green, error=rm)), 1)
        new_fp = str(_members(green)["fp_green"])
        self.assertNotEqual(new_fp, old_fp)
        with self.assertRaises(ValueError) as cm:
            load(path, RUNS["dense_npz"].solver)
        self.assertEqual(str(cm.exception), FP_GREEN_MISMATCH.format(
            path=path, label=LABEL, run=new_fp[:12], file=old_fp[:12]))

    def test_green_tmp_removed_on_failure(self):
        """A writer exception before the green commit: the temporary green
        file is removed, the old green.npz is unchanged and the original
        exception propagates; a failed removal is one WARNING."""
        from hwave.solver import so_bond_archive as sa
        boom = RuntimeError("writer failed before the commit")
        path = self.copy()
        d = os.path.dirname(path)
        green = os.path.join(d, "green.npz")
        with open(green, "rb") as fh:
            before = fh.read()
        exc, _ = self._rerun(path, [mock.patch.object(sa, "write_so_bond_archive",
                                                      side_effect=boom)])
        self.assertIs(exc, boom)
        self.assertEqual(_tmp_greens(d), [])
        with open(green, "rb") as fh:
            self.assertEqual(fh.read(), before)
        load(path, RUNS["dense_npz"].solver)
        # the removal of the temporary green file fails as well
        rm = PermissionError("tmp removal refused")
        real_remove = os.remove
        seen = []

        def remove(p, *a, **kw):
            if os.path.basename(p).startswith(".green."):
                seen.append(p)
                raise rm
            return real_remove(p, *a, **kw)
        exc, records = self._rerun(path, [mock.patch.object(sa, "write_so_bond_archive",
                                                            side_effect=boom),
                                          mock.patch.object(os, "remove", remove)])
        self.assertIs(exc, boom)
        self.assertEqual(len(seen), 1)
        self.assertEqual([m for m in messages(records, "WARNING") if "could not remove" in m],
                         [T_CLEANUP_FAILED.format(name=seen[0], error=rm)])
        os.remove(seen[0])


if __name__ == "__main__":
    unittest.main()
