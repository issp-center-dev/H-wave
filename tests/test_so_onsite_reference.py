"""Opt-in bitwise comparison of the spin-orbital on-site outputs (and of the
static spin-orbital hwave_sc bond runs) with the reference revision
:data:`SO_REFERENCE_COMMIT`.

The spin-orbital bond-channel FLEX added after that revision must leave
every spin-orbital path it does not enable unchanged. This module runs the
same inputs on the reference checkout and on this tree, each in its own
subprocess (``PYTHONPATH=<tree>/src``, BLAS/OpenMP threads pinned to one),
and compares the outputs:

* the two output directories list the same file names;
* every member of every ``.npz`` archive is equal (``numpy.array_equal`` with
  the same dtype and shape; strings by value), except exactly the member
  ``producer_version`` -- ``hwave.__version__`` plus ``git describe``, which
  names the producing tree and so differs between the two -- which must be
  present in this tree's archive wherever the reference archive has it;
* every other file (``energy.dat``, ``eigenvalue.dat``, ...) byte for byte,
  except the value of the same member where it is written as a text header
  line (``# producer_version = ...`` in ``eigenvalue.dat``): the line itself
  must be at the same place in both files.

The static bond runs set ``[eliashberg] bond_memory_cap_gb = 8.0``: unset,
the cap is 0.8 x the memory psutil reports available at run time, and
``eigenvalue.dat`` records it (``# bond_memory_cap = ...``), so two runs of
one tree differ there.

The runs (inputs written once, by this tree's test writers, and read by both
trees):

* ``flex_onsite`` -- spin-orbital on-site FLEX (``calc_scheme = "general"``,
  ``flex_hartree_fock = true``, bond gate off), the two-orbital chain of
  ``tests/test_flex_spin_orbital.py`` with a spin-mixing on-site block and
  every on-site interaction type;
* ``rpa_onsite`` -- spin-orbital RPA (``calc_scheme = "general"``) of the
  Rashba-like square lattice of ``tests/test_sc_spin_orbital.py``, and
  ``rpa_onsite_reduced`` -- the same with ``calc_scheme = "reduced"`` on the
  density-diagonal types of that model (CoulombIntra, CoulombInter, Hund,
  Ising; reduced refuses Exchange and PairHop). The ``rpa.py`` changes after
  the reference revision (the ``Interaction`` constructor flag and the
  ``_so_provenance`` signature) sit on the path both schemes share;
  ``rpa_onsite_6x5`` -- the general case on a 6 x 5 grid, where the two
  expressions of the unnormalised inverse FFT round differently (the
  one-body H(k) must keep RPA's ``ifftn * nvol``; ``fp_green`` included);
* ``sc_static_onsite`` -- hwave_sc, ``frequency = "static"``,
  ``chi0q_mode = "calc"``, on that model;
* ``sc_bond_soc_all_types`` and ``sc_bond_zeeman_extern`` -- hwave_sc static
  with ``[eliashberg] bond_channels = true`` on the ``SOC_ALL_TYPES`` and
  ``ZEEMAN_EXTERN`` fixtures of ``tests/so_bond_entry_fixtures.py`` (the
  static pole indicator and the bond-run assembly of that path were
  refactored after the reference revision).

Each subprocess reports the ``hwave`` package it imported; the test
asserts it lies under that tree's ``src`` and that the reference side's
``producer_version`` names the reference commit while this tree's does not,
so no equality row can pass by running one tree twice.

Non-vacuity rows: the ``flex_onsite`` input plus the gate-only key
``longitudinal_bond_max_shells = 2`` with the gate off is inert on both trees
(the only effect is the reference revision's two "ignored" WARNINGs, of the
FLEX and of the RPA parameter readers)
and gives the same outputs; a different ``Nmat`` gives different outputs.

NUMPY >= 2 ONLY. With numpy < 2 (1.26.4 measured) two runs of the SAME tree
on the same input already differ at round-off (5.6e-17 in ``chi0q`` of the
RPA case, and in the FLEX and bond hwave_sc members), so a bitwise
comparison there says nothing about the code; it is a unittest skip with
that reason (after the checkout check, so an unusable checkout still fails
under the require flag).

OPT-IN. The comparison runs only when ``HWAVE_SO_REFERENCE_CHECKOUT`` names
a clean checkout at :data:`SO_REFERENCE_COMMIT` (verified by
:func:`tests.test_flex_second_order_compat.develop_checkout`); otherwise it
is a unittest skip, and with ``HWAVE_REQUIRE_SO_REFERENCE=1`` a failure. CI
does not provision that checkout yet."""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np

#: the reference revision of this comparison
SO_REFERENCE_COMMIT = "a92b533d91246565f2785f44a5d2f91f23054b67"
#: the reference checkout
SO_REFERENCE_PATH_ENV = "HWAVE_SO_REFERENCE_CHECKOUT"
#: set to a true-ish value: an unusable reference checkout is a failure
SO_REFERENCE_REQUIRE_ENV = "HWAVE_REQUIRE_SO_REFERENCE"

#: numpy < 2 is not run-to-run reproducible on these paths (module docstring)
NUMPY_BELOW_2 = int(np.__version__.split(".")[0]) < 2

#: the one archive member allowed to differ (it names the producing tree)
EXCLUDED_MEMBERS = frozenset({"producer_version"})

#: the gate-only key of the inert row, and the WARNING it gets with the gate off
INERT_KEY = ("longitudinal_bond_max_shells", 2)
INERT_WARNINGS = [
    ("WARNING", "hwave.solver.flex",
     "[mode.param] longitudinal_bond_max_shells set but longitudinal_bond_channels is not "
     "true; these bond-only options are ignored (not parsed)."),
    ("WARNING", "hwave.solver.rpa",
     "[mode.param] longitudinal_bond_max_shells set but longitudinal_bond_channels=false; "
     "these options only apply to longitudinal_bond_channels=true and are ignored here."),
]
#: the text form of the excluded member (a header line of ``eigenvalue.dat``)
EXCLUDED_LINE_PREFIXES = tuple("# {} = ".format(m).encode() for m in sorted(EXCLUDED_MEMBERS))
#: the fixed host cap of the static bond runs (unset, it follows psutil)
BOND_MEMORY_CAP_GB = 8.0

_RUN = r'''
import json, os, sys
import hwave
sys.stderr.write("HWAVE_PACKAGE=" + os.path.realpath(os.path.dirname(hwave.__file__)) + "\n")
sys.stderr.flush()
kind, path = sys.argv[1:3]
with open(path) as fh:
    inp = json.load(fh)
if kind == "hwave":
    from hwave import qlms
    qlms.run(input_dict=inp)
elif kind == "sc":
    import hwave.sc as sc
    sc.calc_eliashberg(inp)
else:
    raise SystemExit("unknown kind " + kind)
'''

_LOG_LINE = re.compile(r"^\S+ \S+ (WARNING|ERROR|CRITICAL) (\S+): (.*)$")


def _flex_input(d, out, **par_extra):
    from tests import test_flex_spin_orbital as F
    files = F._write_model(d, True, types=F._SU2_TYPES + ("PairLift",), soc=0.3)
    inter = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
    inter.update(files)
    par = {"T": 0.5, "filling": 0.4, "CellShape": [4, 1, 1], "SubShape": [1, 1, 1],
           "Nmat": 16, "IterationMax": 100, "Mix": 0.5, "EPS": 12,
           "mixing_scheme": "linear", "flex_hartree_fock": True}
    par.update(par_extra)
    return {"log": {"print_level": 1},
            "mode": {"mode": "FLEX", "enable_spin_orbital": True, "calc_scheme": "general",
                     "param": par},
            "file": {"input": {"path_to_input": "", "interaction": inter},
                     "output": {"path_to_output": out, "chi0q": "chi0q", "chiq": "chiq",
                                "sigma": "sigma", "green": "green", "energy": "energy.dat"}}}


def _sc_model(d):
    from tests import test_sc_spin_orbital as S
    return S, S._write_model(d, so=True, soc=0.3)


def _rpa_input(d, out, scheme="general", cell=None):
    if scheme == "reduced":
        from tests import test_sc_spin_orbital as S
        files = S._write_model(d, so=True, soc=0.3,
                               types=("CoulombIntra", "CoulombInter", "Hund", "Ising"))
    else:
        S, files = _sc_model(d)
    inp = S._input(d, files, True, {})
    inp["mode"]["calc_scheme"] = scheme
    if cell is not None:
        inp["mode"]["param"]["CellShape"] = list(cell)
    del inp["eliashberg"]
    inp["log"] = {"print_level": 1}
    inp["file"]["output"] = {"path_to_output": out, "chi0q": "chi0q", "chiq": "chiq"}
    return inp


def _sc_onsite_input(d, out):
    S, files = _sc_model(d)
    inp = S._input(d, files, True, {"frequency": "static", "chi0q_mode": "calc",
                                    "solver_mode": "eigenvalue", "num_eigenvalues": 8})
    inp["file"]["output"]["path_to_output"] = out
    return inp


def _sc_bond_input(fixture_name):
    def build(d, out):
        from tests import so_bond_entry_fixtures as E
        inp = E.write_case(d, getattr(E, fixture_name),
                           eliashberg={"bond_memory_cap_gb": BOND_MEMORY_CAP_GB})
        inp["file"]["output"]["path_to_output"] = out
        return inp
    return build


#: name -> (subprocess kind, input builder (input dir, output dir) -> input dict)
CASES = {
    "flex_onsite": ("hwave", _flex_input),
    "rpa_onsite": ("hwave", _rpa_input),
    "rpa_onsite_reduced": ("hwave", lambda d, out: _rpa_input(d, out, scheme="reduced")),
    # a grid side that is not a power of two: ifftn * nvol and
    # ifftn(norm="forward") round differently there, and the archive's
    # Green-function digest (fp_green) and state members would show it
    "rpa_onsite_6x5": ("hwave", lambda d, out: _rpa_input(d, out, cell=(6, 5, 1))),
    "sc_static_onsite": ("sc", _sc_onsite_input),
    "sc_bond_soc_all_types": ("sc", _sc_bond_input("SOC_ALL_TYPES")),
    "sc_bond_zeeman_extern": ("sc", _sc_bond_input("ZEEMAN_EXTERN")),
}


def _files(root):
    out = []
    for top, _, names in os.walk(root):
        for n in names:
            out.append(os.path.relpath(os.path.join(top, n), root))
    return sorted(out)


def _same_member(a, b):
    if a.dtype != b.dtype or a.shape != b.shape:
        return False
    if a.dtype.kind in "fc":
        return bool(np.array_equal(a, b, equal_nan=True))
    return bool(np.array_equal(a, b))


def _same_text(a, b, prefixes=EXCLUDED_LINE_PREFIXES):
    """Byte for byte, except the value of a line that starts with one of
    ``prefixes`` (at the same line number in both files)."""
    if a == b:
        return True
    la, lb = a.split(b"\n"), b.split(b"\n")
    if len(la) != len(lb):
        return False
    for x, y in zip(la, lb):
        if x == y:
            continue
        if not any(x.startswith(p) and y.startswith(p) for p in prefixes):
            return False
    return True


def output_mismatches(ref_dir, new_dir, excluded=EXCLUDED_MEMBERS):
    """The list of differences between the reference output directory and
    this tree's (empty when they agree): file sets, ``.npz`` member sets and
    values (``excluded`` members only required present where the reference
    has them), other files byte for byte."""
    bad = []
    ref_files, new_files = _files(ref_dir), _files(new_dir)
    if ref_files != new_files:
        bad.append("file sets differ: only reference {}, only this tree {}".format(
            sorted(set(ref_files) - set(new_files)), sorted(set(new_files) - set(ref_files))))
    for fn in sorted(set(ref_files) & set(new_files)):
        a, b = os.path.join(ref_dir, fn), os.path.join(new_dir, fn)
        if not fn.endswith(".npz"):
            with open(a, "rb") as fa, open(b, "rb") as fb:
                if not _same_text(fa.read(), fb.read()):
                    bad.append("{}: bytes differ".format(fn))
            continue
        with np.load(a, allow_pickle=False) as da, np.load(b, allow_pickle=False) as db:
            ra, rb = set(da.files), set(db.files)
            if ra != rb:
                bad.append("{}: member sets differ: only reference {}, only this tree {}"
                           .format(fn, sorted(ra - rb), sorted(rb - ra)))
            for m in sorted(ra & rb):
                if m in excluded:
                    continue
                if not _same_member(da[m], db[m]):
                    bad.append("{}: member {} differs".format(fn, m))
    return bad


def _warnings(stderr):
    """The ``(level, logger, message)`` of every WARNING-or-worse log line."""
    out = []
    for line in stderr.splitlines():
        m = _LOG_LINE.match(line)
        if m:
            out.append(m.groups())
    return out


class TestOutputMismatches(unittest.TestCase):
    """The comparison helper, without a reference checkout."""

    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="so_ref_helper_")
        self.a, self.b = os.path.join(self.d, "a"), os.path.join(self.d, "b")
        for x in (self.a, self.b):
            os.makedirs(x)

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def _write(self, root, version="1 (abc)", value=1.0, text="x\n", extra=None,
               header="# producer_version = 1 (abc)\n"):
        members = {"producer_version": np.str_(version), "v": np.array([value, np.nan]),
                   "tag": np.str_("spin_orbital")}
        members.update(extra or {})
        np.savez(os.path.join(root, "g.npz"), **members)
        with open(os.path.join(root, "energy.dat"), "w") as fh:
            fh.write(text)
        with open(os.path.join(root, "eigenvalue.dat"), "w") as fh:
            fh.write("# a = 1\n" + header + "0.5\n")

    def test_only_producer_version_may_differ(self):
        self._write(self.a)
        self._write(self.b, version="1 (def-dirty)",
                    header="# producer_version = 1 (def-dirty)\n")
        self.assertEqual(output_mismatches(self.a, self.b), [])

    def test_every_other_difference_is_reported(self):
        self._write(self.a)
        for kw, expect in (({"value": 2.0}, "member v differs"),
                           ({"text": "y\n"}, "energy.dat: bytes differ"),
                           ({"extra": {"tag": np.str_("other")}}, "member tag differs"),
                           ({"extra": {"new": np.zeros(1)}}, "member sets differ"),
                           ({"header": ""}, "eigenvalue.dat: bytes differ"),
                           ({"header": "# producer_versions = 1 (abc)\n"},
                            "eigenvalue.dat: bytes differ"),
                           ({"header": "# bond_memory_cap = 1 (abc)\n"},
                            "eigenvalue.dat: bytes differ")):
            with self.subTest(kw=kw):
                self._write(self.b, **kw)
                bad = output_mismatches(self.a, self.b)
                self.assertTrue(any(expect in s for s in bad), bad)

    def test_missing_producer_version_is_reported(self):
        self._write(self.a)
        np.savez(os.path.join(self.b, "g.npz"), v=np.array([1.0, np.nan]),
                 tag=np.str_("spin_orbital"))
        for fn in ("energy.dat", "eigenvalue.dat"):
            shutil.copy(os.path.join(self.a, fn), os.path.join(self.b, fn))
        self.assertTrue(any("producer_version" in s for s in output_mismatches(self.a, self.b)))

    def test_file_sets(self):
        self._write(self.a)
        self._write(self.b)
        with open(os.path.join(self.b, "extra.dat"), "w") as fh:
            fh.write("")
        self.assertTrue(any("file sets differ" in s for s in output_mismatches(self.a, self.b)))


class TestSpinOrbitalOnsiteReference(unittest.TestCase):
    """The spin-orbital on-site and static hwave_sc outputs of this tree
    equal the reference revision's (module docstring)."""

    @classmethod
    def setUpClass(cls):
        from tests.test_flex_second_order_compat import develop_checkout
        ref, why = develop_checkout(expected_commit=SO_REFERENCE_COMMIT,
                                    path_env=SO_REFERENCE_PATH_ENV,
                                    require_env=SO_REFERENCE_REQUIRE_ENV)
        if ref is None:
            raise unittest.SkipTest(why)
        if NUMPY_BELOW_2:
            raise unittest.SkipTest(
                "numpy {} < 2: two runs of one tree differ at round-off, so the bitwise "
                "comparison needs numpy >= 2".format(np.__version__))
        cls.ref = ref
        cls.mine = os.getcwd()
        cls.work = tempfile.mkdtemp(prefix="so_onsite_reference_")
        cls.runs = {}
        try:
            for name, (kind, build) in CASES.items():
                cls._both(name, kind, build)
            flex_kind, flex_build = CASES["flex_onsite"]
            cls._both("flex_onsite_inert_key", flex_kind,
                      lambda d, out: flex_build(d, out, **dict([INERT_KEY])))
            cls.runs["flex_onsite_other_nmat", "mine"] = cls._run(
                cls.mine, "flex_onsite_other_nmat", flex_kind,
                lambda d, out: flex_build(d, out, Nmat=24))
        except BaseException:
            # tearDownClass does not run when setUpClass raises
            shutil.rmtree(cls.work, ignore_errors=True)
            raise

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.work, ignore_errors=True)

    @classmethod
    def _both(cls, name, kind, build):
        for side, tree in (("reference", cls.ref), ("mine", cls.mine)):
            cls.runs[name, side] = cls._run(tree, name + "_" + side, kind, build)

    @classmethod
    def _run(cls, tree, tag, kind, build):
        """Write the inputs into a fresh directory, run them on ``tree`` and
        return ``(output dir, stderr)``."""
        from tests.test_flex_second_order_compat import (own_tree_subprocess_env,
                                                         reference_subprocess_env)
        base = os.path.join(cls.work, tag)
        out = os.path.join(base, "out")
        os.makedirs(base)
        # the inputs live under one fixed name per case, so that every
        # recorded input path is the same on both trees
        shared = os.path.join(cls.work, "input_" + tag.rsplit("_", 1)[0])
        if os.path.exists(shared):
            shutil.rmtree(shared)
        os.makedirs(shared)
        inp = build(shared, out)
        path = os.path.join(base, "input.json")
        with open(path, "w") as fh:
            json.dump(inp, fh)
        env = (own_tree_subprocess_env(tree) if tree == cls.mine
               else reference_subprocess_env(tree))
        env["PYTHONDONTWRITEBYTECODE"] = "1"
        proc = subprocess.run([sys.executable, "-B", "-c", _RUN, kind, path], env=env,
                              cwd=base, capture_output=True, text=True)
        if proc.returncode != 0:
            raise AssertionError("{} failed on {} (exit {}):\n{}".format(
                tag, tree, proc.returncode, proc.stderr[-4000:]))
        return out, proc.stderr

    def _assert_same(self, name, ref_name=None):
        ref_out, _ = self.runs[ref_name or name, "reference"]
        new_out, _ = self.runs[name, "mine"]
        files = _files(ref_out)
        self.assertTrue(files, "{}: the reference wrote no output".format(name))
        self.assertEqual(output_mismatches(ref_out, new_out), [], name)
        return files

    def test_outputs_equal_the_reference(self):
        expected = {
            "flex_onsite": {"chi0q.npz", "chiq.npz", "sigma.npz", "green.npz", "energy.dat"},
            "rpa_onsite": {"chi0q.npz", "chiq.npz"},
            "rpa_onsite_reduced": {"chi0q.npz", "chiq.npz"},
            "rpa_onsite_6x5": {"chi0q.npz", "chiq.npz"},
            "sc_static_onsite": {"eigenvalue.dat", "gap.dat", "gap.npz"},
            "sc_bond_soc_all_types": {"eigenvalue.dat", "gap.dat", "gap.npz"},
            "sc_bond_zeeman_extern": {"eigenvalue.dat", "gap.dat", "gap.npz"},
        }
        for name in CASES:
            with self.subTest(case=name):
                files = self._assert_same(name)
                self.assertTrue(expected[name] <= set(files), (name, files))
                # producer_version is recorded and excluded, never silently absent
                with self.subTest(case=name, check="producer_version"):
                    self._assert_producer_version_kept(name)

    def _assert_producer_version_kept(self, name):
        ref_out, _ = self.runs[name, "reference"]
        new_out, _ = self.runs[name, "mine"]
        for fn in _files(ref_out):
            if not fn.endswith(".npz"):
                continue
            with np.load(os.path.join(ref_out, fn)) as da, \
                    np.load(os.path.join(new_out, fn)) as db:
                if "producer_version" in da.files:
                    self.assertIn("producer_version", db.files, (name, fn))

    def test_each_side_ran_its_own_tree(self):
        """The imported ``hwave`` lies under the side's own ``src``, and the
        recorded ``producer_version`` names the reference commit on the
        reference side only."""
        marker = SO_REFERENCE_COMMIT[:7]
        trees = {"reference": self.ref, "mine": self.mine}
        for (name, side), (out, err) in sorted(self.runs.items()):
            with self.subTest(case=name, side=side):
                pkgs = re.findall(r"^HWAVE_PACKAGE=(.*)$", err, flags=re.M)
                self.assertEqual(pkgs, [os.path.realpath(
                    os.path.join(trees[side], "src", "hwave"))])
                versions = []
                for fn in _files(out):
                    if fn.endswith(".npz"):
                        with np.load(os.path.join(out, fn)) as z:
                            if "producer_version" in z.files:
                                versions.append(str(z["producer_version"]))
                self.assertTrue(versions, (name, side))
                for v in versions:
                    if side == "reference":
                        self.assertIn(marker, v)
                    else:
                        self.assertNotIn(marker, v)

    def test_inert_gate_only_key(self):
        """``longitudinal_bond_max_shells = 2`` with the gate off: the same
        outputs as the reference with and without the key, and the only
        added WARNINGs are the two "ignored" ones, on both trees."""
        self._assert_same("flex_onsite_inert_key")
        self._assert_same("flex_onsite_inert_key", ref_name="flex_onsite")
        for side in ("reference", "mine"):
            with self.subTest(side=side):
                plain = _warnings(self.runs["flex_onsite", side][1])
                keyed = _warnings(self.runs["flex_onsite_inert_key", side][1])
                added = list(keyed)
                for w in plain:
                    if w in added:
                        added.remove(w)
                self.assertEqual(added, INERT_WARNINGS, (side, keyed))

    def test_a_different_nmat_differs(self):
        ref_out, _ = self.runs["flex_onsite", "reference"]
        new_out, _ = self.runs["flex_onsite_other_nmat", "mine"]
        bad = output_mismatches(ref_out, new_out)
        self.assertTrue(any("sigma.npz" in s for s in bad), bad)
        self.assertTrue(any("energy.dat" in s for s in bad), bad)


if __name__ == "__main__":
    unittest.main()
