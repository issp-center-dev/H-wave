"""The ``python -m hwave.initial`` commands prepare, check and collect.

Every test writes under ``tempfile.mkdtemp()``. ``hwave`` and ``python``
are run through small shims on ``PATH`` that start the interpreter running
these tests, with ``PYTHONPATH`` at the package under test, so
``commands.txt`` runs this checkout whatever is installed.
"""
import contextlib
import io
import json
import os
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import unittest

import numpy as np

import hwave
from tests import initial_fixtures
from tests.heavy_tests import heavy
from tests.initial_reference_capture import dump_toml

_SRC = os.path.dirname(os.path.dirname(os.path.abspath(hwave.__file__)))

_INTERNAL = re.compile(r"codex|ai-review|antigravity|probe|superpowers|\bsdd\b|claude lens|"
                       r"ledger|§|design section|\bplan\b|\bO[0-9]\b|\bR[0-9]\b", re.I)


def _read_toml(path):
    import tomli
    with open(path, "rb") as f:
        return tomli.load(f)


def _write_toml(path, doc):
    with open(path, "w") as f:
        f.write(dump_toml(doc))


def _edit(path, fn):
    doc = _read_toml(path)
    fn(doc)
    _write_toml(path, doc)


def _sha(path):
    import hashlib
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def _geometry_uhf(L):
    text = "1 0 0\n0 1 0\n0 0 1\n0 0 0\n{} 0 0\n0 {} 0\n0 0 1\n".format(L, L)
    return text + "".join("{} {} 0 0\n".format(x, y) for y in range(L) for x in range(L))


_ONEBODYG = "=====\nNCisAjs 2\n=====\n=====\n=====\n0 0 0 0\n0 0 1 0\n"



# --------------------------------------------------------------------------
# the staging cases, shared by the candidate runs (TestPrepare) and the
# directories of the preliminary paramagnetic run (test_initial_stage0)
# --------------------------------------------------------------------------
#
# STAGING_CASES is a list of ``(name, build, expected)``:
#
# * ``build(tc, name) -> (toml, info)`` writes a source input under
#   ``tc.path(name)`` and returns its TOML and what the checks need;
# * ``expected`` is ``{"refused": words}`` (prepare exits 1, naming every
#   word, and writes no WORK; for the TOML and for its directory) or
#   ``{"check": fn, "after_run": fn}``, each ``fn(tc, target_dir, info)``
#   asserting on one target directory (a candidate run or a preliminary
#   paramagnetic run): ``check`` after prepare wrote it, ``after_run``
#   after hwave ran in it (``None``: nothing to check).


def _square_source(tc, name):
    src = tc.model(name)
    return os.path.join(src, "input.toml")


def _edited_source(edit):
    def build(tc, name):
        toml = _square_source(tc, name)
        _edit(toml, edit)
        return toml, {}
    return build


def _check_output_paths(tc, rd, info):
    from hwave.results import open_run
    doc = _read_toml(os.path.join(rd, "input.toml"))
    tc.assertEqual(doc["file"]["output"]["path_to_output"], "out/dir")
    # an omitted energy is written, collect reads it
    tc.assertEqual(doc["file"]["output"]["energy"], "energy.dat")
    run = open_run(os.path.join(rd, "input.toml"), base_dir=rd)
    tc.assertEqual(str(run.output_dir), os.path.join(os.path.realpath(rd), "out", "dir"))


def _after_output_paths(tc, rd, info):
    tc.assertTrue(os.path.isfile(os.path.join(rd, "out", "dir", "energy.dat")))


def _build_input_staging(tc, name):
    src = tc.model(name)
    toml = os.path.join(src, "input.toml")
    # absolute, nested and '..' sources; a basename shared by the two tables
    other = tc.path(name + "_other")
    os.makedirs(os.path.join(src, "nested"))
    os.makedirs(other)
    shutil.move(os.path.join(src, "transfer.dat"), os.path.join(src, "nested",
                                                                 "transfer.dat"))
    shutil.move(os.path.join(src, "coulombintra.dat"),
                os.path.join(other, "coulombintra.dat"))
    shutil.move(os.path.join(src, "geom.dat"), os.path.join(other, "geom.dat"))
    os.makedirs(os.path.join(src, "uhf"))
    with open(os.path.join(src, "uhf", "geom.dat"), "w") as f:
        f.write(_geometry_uhf(4))
    with open(os.path.join(src, "uhf", "OneBodyG.dat"), "w") as f:
        f.write(_ONEBODYG)
    rel_other = os.path.join("..", "..", os.path.basename(other), "coulombintra.dat")

    def fn(d):
        fin = d["file"]["input"]
        fin.update(path_to_input="uhf", geometry_uhf="geom.dat",
                   onebodyg_uhf="OneBodyG.dat")
        fin["interaction"] = {"path_to_input": "nested",
                              "Geometry": os.path.join(other, "geom.dat"),
                              "Transfer": "transfer.dat",
                              "CoulombIntra": rel_other}
        d["file"]["output"]["onebodyg"] = "greenone.dat"
    _edit(toml, fn)
    sources = {"Geometry": os.path.join(other, "geom.dat"),
               "Transfer": os.path.join(src, "nested", "transfer.dat"),
               "CoulombIntra": os.path.join(other, "coulombintra.dat"),
               "geometry_uhf": os.path.join(src, "uhf", "geom.dat"),
               "onebodyg_uhf": os.path.join(src, "uhf", "OneBodyG.dat")}
    return toml, {"sources": sources, "roots": [src, other]}


def _check_input_staging(tc, rd, info):
    from hwave.solver import initial_facade
    sources = info["sources"]
    doc = _read_toml(os.path.join(rd, "input.toml"))
    fin = doc["file"]["input"]
    # TOML key order: [file.input] keys come before its interaction table
    tc.assertEqual(fin["geometry_uhf"], "inputs/geom.dat")
    tc.assertEqual(fin["interaction"]["Geometry"], "geom-2.dat")
    files = initial_facade.referenced_input_files(os.path.join(rd, "input.toml"))
    tc.assertEqual(set(files), set(sources))
    for key, path in files.items():
        tc.assertTrue(str(path).startswith(os.path.realpath(rd)) or
                      str(path).startswith(rd), (key, path))
        tc.assertFalse(os.path.islink(path))
        tc.assertTrue(stat.S_ISREG(os.lstat(path).st_mode))
        tc.assertEqual(_sha(path), _sha(sources[key]), key)


def _after_greenone(tc, rd, info):
    tc.assertTrue(os.path.isfile(os.path.join(rd, "output", "greenone.dat")))


def _inputs_dir_source(inter_dir):
    """A source TOML with [file.input] path_to_input = "inputs" and the
    model files under ``inter_dir``."""
    def build(tc, name):
        base = tc.path(name)
        model_toml = initial_fixtures.square(os.path.join(base, inter_dir), 4, U=4.0,
                                             subshape=(2, 2, 1), param={"IterationMax": 1})
        doc = _read_toml(model_toml)
        os.remove(model_toml)
        os.makedirs(os.path.join(base, "inputs"), exist_ok=True)
        with open(os.path.join(base, "inputs", "geometry_uhf.dat"), "w") as f:
            f.write(_geometry_uhf(4))
        with open(os.path.join(base, "inputs", "OneBodyG.dat"), "w") as f:
            f.write(_ONEBODYG)
        doc["file"]["input"].update(path_to_input="inputs", geometry_uhf="geometry_uhf.dat",
                                    onebodyg_uhf="OneBodyG.dat")
        doc["file"]["input"]["interaction"]["path_to_input"] = inter_dir
        doc["file"]["output"]["onebodyg"] = "greenone.dat"
        _write_toml(os.path.join(base, "input.toml"), doc)
        return os.path.join(base, "input.toml"), {}
    return build


def _check_path_to_input(tc, rd, info):
    from hwave.initial.writer import solver_initial_path
    t = os.path.join(rd, "input.toml")
    fin = _read_toml(t)["file"]["input"]
    tc.assertEqual(fin["path_to_input"], "")
    tc.assertEqual(fin["initial"], "initial.npz")
    tc.assertEqual(fin["geometry_uhf"], "inputs/geometry_uhf.dat")
    tc.assertEqual(fin["interaction"]["path_to_input"], "inputs")
    tc.assertTrue(os.path.samefile(solver_initial_path(t), os.path.join(rd, "initial.npz")))


def _after_path_to_input(tc, rd, info):
    with open(os.path.join(rd, "hwave.log")) as f:
        tc.assertIn("read initial from initial.npz", f.read())
    _after_greenone(tc, rd, info)


STAGING_CASES = [
    # test_exact_case_variant_refused: a case variant of an exact-case key
    ("exact_case_Initial", _edited_source(
        lambda d: d["file"]["input"].update(Initial="x.npz")), {"refused": ("'initial'",)}),
    ("exact_case_SCF_Status", _edited_source(
        lambda d: d["file"]["output"].update(SCF_Status="s.dat")),
     {"refused": ("'scf_status'",)}),
    ("exact_case_Mode", _edited_source(lambda d: d.update(Mode=d.pop("mode"))),
     {"refused": ("[mode]",)}),
    # test_prepare_output_paths
    ("output_paths", _edited_source(
        lambda d: d["file"]["output"].update(path_to_output="out/dir")),
     {"check": _check_output_paths, "after_run": _after_output_paths}),
] + [
    ("output_bad_{}_{}".format(key, i), _edited_source(
        lambda d, key=key, value=value: d["file"]["output"].update({key: value})),
     {"refused": (key, word)})
    for i, (key, value, word) in enumerate((
        ("path_to_output", "/abs/out", "absolute"), ("path_to_output", "../out", "'..'"),
        ("green", "", 'green = ""'), ("energy", "", 'energy = ""'),
        ("green", "../g", "'..'"), ("energy", "/e.dat", "absolute")))
] + [
    # test_prepare_input_staging: absolute, nested, '..', a basename
    # collision, geometry_uhf/onebodyg_uhf
    ("input_staging", _build_input_staging,
     {"check": _check_input_staging, "after_run": _after_greenone}),
    # test_prepare_path_to_input_inputs: [file.input] path_to_input = "inputs"
    ("path_to_input_inputs", _inputs_dir_source("inputs"),
     {"check": _check_path_to_input, "after_run": _after_path_to_input}),
    ("path_to_input_data", _inputs_dir_source("data"),
     {"check": _check_path_to_input, "after_run": _after_path_to_input}),
]


def staging_cases(prefix=""):
    """The staging cases whose name starts with ``prefix``."""
    return [c for c in STAGING_CASES if c[0].startswith(prefix)]


def check_refused_case(tc, case, *extra):
    """Run ``prepare`` on the refused ``case`` (its TOML, then its
    directory) with the arguments ``extra``: exit code 1, every expected
    word in the message, and no WORK."""
    name, build, expected = case
    toml, _ = build(tc, name)
    for inp in (toml, os.path.dirname(toml)):
        work = tc.path("W_" + name)
        rc, _, err = tc.main("prepare", inp, "--out", work, *extra)
        tc.assertEqual(rc, 1, (name, inp))
        for word in expected["refused"]:
            tc.assertIn(word, err, name)
        tc.assertFalse(os.path.exists(work), name)


def _pin_tier(argv):
    """``prepare`` argv with ``--tier large`` unless a tier is given: these
    tests check the large tier (the default tier depends on the host
    memory)."""
    if argv and argv[0] == "prepare" and "--tier" not in argv and "--help" not in argv:
        return argv + ["--tier", "large"]
    return argv

class _Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        self.bin = os.path.join(self.tmp, "_bin")
        os.makedirs(self.bin)
        for name, body in (("python", '"$@"'),
                           ("hwave", '-c "import sys; from hwave.qlms import main; '
                                     'sys.argv[0] = \'hwave\'; main()" "$@"')):
            p = os.path.join(self.bin, name)
            with open(p, "w") as f:
                f.write('#!/bin/sh\nexec "{}" {}\n'.format(sys.executable, body))
            os.chmod(p, os.stat(p).st_mode | stat.S_IEXEC)
        self.env = dict(os.environ, PATH=self.bin + os.pathsep + os.environ.get("PATH", ""),
                        PYTHONPATH=_SRC, OPENBLAS_NUM_THREADS="1")

    def path(self, *parts):
        return os.path.join(self.tmp, *parts)

    def model(self, name="src", *, L=4, U=4.0, T=0.0, param=None, subshape=(2, 2, 1)):
        p = {"IterationMax": 1}
        p.update(param or {})
        toml = initial_fixtures.square(self.path(name), L, U=U, T=T, subshape=subshape,
                                       param=p)
        return os.path.dirname(toml)

    def prepare(self, src, work, **kw):
        from hwave.initial import cli
        kw.setdefault("n_random", 0)
        kw.setdefault("max_candidates", 0)
        kw.setdefault("tier", "large")
        out = io.StringIO()
        rc = cli.prepare(src, work, stdout=out, **kw)
        self.assertEqual(rc, 0)
        return out.getvalue()

    def main(self, *argv):
        from hwave.initial import cli
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = cli.main(_pin_tier([str(a) for a in argv]))
        return rc, out.getvalue(), err.getvalue()

    def manifest(self, work):
        with open(os.path.join(work, "manifest.json")) as f:
            return json.load(f)

    def run_ids(self, work):
        return [e["run_id"] for e in self.manifest(work)["runs"]]

    def bash(self, script, cwd=None):
        return subprocess.run(["bash", script], cwd=cwd or self.tmp, env=self.env,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    def run_hwave(self, run_dir, log="hwave.log"):
        with open(os.path.join(run_dir, log), "w") as f:
            return subprocess.run(["hwave", "input.toml"], cwd=run_dir, env=self.env,
                                  stdout=f, stderr=subprocess.STDOUT).returncode

    def collect(self, work):
        from hwave.initial import cli
        out, err = io.StringIO(), io.StringIO()
        rc = cli.collect(work, stdout=out, stderr=err)
        self.assertEqual(rc, 0)
        with open(os.path.join(work, "coverage.json")) as f:
            doc = json.load(f)
        return {r["run_id"]: r for r in doc["rows"] if r["tried"]}, doc, err.getvalue()


class TestPrepare(_Base):
    @heavy
    def test_prepare_layout_and_commands(self):
        from hwave.initial.writer import solver_initial_path
        src = self.model()
        work = self.path("WORK")
        msg = self.prepare(src, work, max_candidates=2)
        self.assertIn("bytes of staged input", msg)
        ids = self.run_ids(work)
        self.assertTrue(1 <= len(ids) <= 4, ids)
        for name in ("coverage.json", "manifest.json", "commands.txt", "candidates.json"):
            self.assertTrue(os.path.isfile(os.path.join(work, name)), name)
        with open(os.path.join(work, "coverage.json")) as f:
            self.assertEqual(json.load(f)["phase_marker"], "pre_run")
        m = self.manifest(work)
        self.assertGreater(m["staged_bytes"], 0)
        with open(os.path.join(work, "commands.txt")) as f:
            lines = f.read().splitlines()
        self.assertEqual(lines[0], 'cd "$(dirname "$0")" || exit 1')
        self.assertTrue(any("run from WORK" in ln for ln in lines if ln.startswith("#")))
        runs = [ln for ln in lines if ln.startswith("(")]
        self.assertEqual(runs, ['(cd "{}" && python -m hwave.initial check '
                                '--require-provenance input.toml && hwave input.toml > '
                                'hwave.log 2>&1)'.format(i) for i in ids])
        for e in m["runs"]:
            rd = os.path.join(work, e["run_id"])
            doc = _read_toml(os.path.join(rd, "input.toml"))
            self.assertEqual(doc["mode"]["param"]["SubShape"], [2, 2, 1])
            self.assertEqual(doc["file"]["input"]["initial"], "initial.npz")
            self.assertTrue(os.path.samefile(solver_initial_path(os.path.join(rd, "input.toml")),
                                             os.path.join(rd, "initial.npz")))
            self.assertEqual(_sha(os.path.join(rd, "initial.npz")), e["initial_sha256"])
            self.assertEqual(_sha(os.path.join(rd, "input.toml")), e["input_toml_sha256"])
        # bash WORK/commands.txt from another directory
        res = self.bash(os.path.join(work, "commands.txt"), cwd=self.bin)
        self.assertEqual(res.returncode, 0, res.stdout)
        for i in ids:
            rd = os.path.join(work, i)
            self.assertTrue(os.path.isfile(os.path.join(rd, "output", "energy.dat")), i)
            rc, out, err = self.main("check", "--require-provenance",
                                     os.path.join(rd, "input.toml"))
            self.assertEqual(rc, 0, err)
            self.assertIn("verified", out)

    def test_prepare_rewrites_existing_initial_keys(self):
        src = self.model(param={"IterationMax": 7, "EPS": 9, "Mix": 0.25})
        with open(os.path.join(src, "geometry_uhf.dat"), "w") as f:
            f.write(_geometry_uhf(4))
        with open(os.path.join(src, "OneBodyG.dat"), "w") as f:
            f.write(_ONEBODYG)
        _edit(os.path.join(src, "input.toml"), lambda d: (
            d["file"]["input"].update(initial="old.npz", initial_uhf="old_uhf.dat",
                                      initial_mode="random", geometry_uhf="geometry_uhf.dat",
                                      onebodyg_uhf="OneBodyG.dat"),
            d["file"]["output"].update(scf_status="status.dat", onebodyg="greenone.dat")))
        work = self.path("WORK")
        self.prepare(src, work)
        with open(os.path.join(work, "coverage.json")) as f:
            labels = {r["family"]: r["label"] for r in json.load(f)["rows"] if r["tried"]}
        self.assertEqual(labels["para"], "para")
        for i in self.run_ids(work):
            doc = _read_toml(os.path.join(work, i, "input.toml"))
            fin = doc["file"]["input"]
            self.assertEqual(fin["initial"], "initial.npz")
            self.assertEqual(fin["path_to_input"], "")
            for k in ("initial_uhf", "initial_mode"):
                self.assertNotIn(k, fin)
            self.assertEqual(fin["geometry_uhf"], "inputs/geometry_uhf.dat")
            self.assertEqual(fin["onebodyg_uhf"], "inputs/OneBodyG.dat")
            p = doc["mode"]["param"]
            self.assertEqual((p["IterationMax"], p["EPS"], p["Mix"]), (7, 9, 0.25))
            self.assertEqual(doc["file"]["output"]["scf_status"], "status.dat")
            with open(os.path.join(work, i, "input.toml")) as f:
                self.assertTrue(f.readline().startswith("# Written by python -m "
                                                        "hwave.initial prepare"))

    def test_prepare_mixed_case_keys(self):
        src = self.model(subshape=None)
        def fn(d):
            p = d["mode"]["param"]
            ncond = p.pop("Ncond")
            p.pop("IterationMax")
            p["subshape"] = [1, 1, 1]
            p["FILLING"] = ncond / 32
            p["iterationmax"] = 3
            inter = d["file"]["input"]["interaction"]
            d["file"]["input"]["interaction"] = {
                ("PATH_TO_INPUT" if k == "path_to_input" else k): v
                for k, v in inter.items()}
        _edit(os.path.join(src, "input.toml"), fn)
        work = self.path("WORK")
        self.prepare(src, work, subshapes=[(2, 2, 1)])
        for i in self.run_ids(work):
            doc = _read_toml(os.path.join(work, i, "input.toml"))
            p = doc["mode"]["param"]
            self.assertNotIn("SubShape", p)
            self.assertEqual(p["subshape"], [2, 2, 1])
            # replaced in place: the key keeps its position
            self.assertEqual(list(p).index("subshape"),
                             list(_read_toml(os.path.join(src, "input.toml"))["mode"]
                                  ["param"]).index("subshape"))
            self.assertEqual((p["FILLING"], p["iterationmax"]), (0.5, 3))
            inter = doc["file"]["input"]["interaction"]
            self.assertEqual(inter["PATH_TO_INPUT"], "inputs")
            self.assertNotIn("path_to_input", inter)
            self.assertEqual(inter["Geometry"], "geom.dat")
        # two spellings of one case-insensitive key are refused
        _edit(os.path.join(src, "input.toml"),
              lambda d: d["mode"]["param"].update(SubShape=[2, 2, 1]))
        rc, _, err = self.main("prepare", src, "--out", self.path("W2"))
        self.assertEqual(rc, 1)
        self.assertIn("'subshape' and 'SubShape'", err)
        self.assertFalse(os.path.exists(self.path("W2")))

    def test_exact_case_variant_refused(self):
        for case in staging_cases("exact_case_"):
            check_refused_case(self, case)

    def test_prepare_output_paths(self):
        for case in staging_cases("output_"):
            name, build, expected = case
            if "refused" in expected:
                check_refused_case(self, case)
                continue
            toml, info = build(self, name)
            work = self.path("W_" + name)
            self.prepare(os.path.dirname(toml), work)
            for i in self.run_ids(work):
                expected["check"](self, os.path.join(work, i), info)

    def _tree(self, root, skip):
        out = {}
        for dp, dns, fns in os.walk(root):
            dns[:] = [d for d in dns if os.path.join(dp, d) != skip]
            for fn in fns:
                p = os.path.join(dp, fn)
                out[p] = _sha(p)
        return out

    def test_prepare_input_staging(self):
        (name, build, expected), = staging_cases("input_staging")
        toml, info = build(self, name)
        work = self.path("WORK")
        before = self._tree(self.tmp, work)
        self.prepare(os.path.dirname(toml), work)
        self.assertEqual(self._tree(self.tmp, work), before)
        ids = self.run_ids(work)
        for i in ids:
            expected["check"](self, os.path.join(work, i), info)
        rd = os.path.join(work, ids[0])
        self.assertEqual(self.run_hwave(rd), 0)
        expected["after_run"](self, rd, info)
        self.assertEqual(self._tree(self.tmp, work), before)

    @heavy
    def test_prepare_path_to_input_inputs(self):
        for name, build, expected in staging_cases("path_to_input_"):
            toml, info = build(self, name)
            work = self.path("W_" + name)
            self.prepare(os.path.dirname(toml), work)
            ids = self.run_ids(work)
            for i in ids:
                expected["check"](self, os.path.join(work, i), info)
            res = self.bash(os.path.join(work, "commands.txt"))
            self.assertEqual(res.returncode, 0, res.stdout)
            for i in ids:
                expected["after_run"](self, os.path.join(work, i), info)

    @heavy
    def test_work_relocatable(self):
        src = self.model(param={"IterationMax": 1})
        work = self.path("a", "WORK")
        self.prepare(src, work)
        for dp, dns, fns in os.walk(work):
            for n in dns + fns:
                self.assertFalse(os.path.islink(os.path.join(dp, n)), n)
            self.assertFalse([n for n in fns if n.startswith("hwave-run.")])

        def check_path(s, where):
            self.assertFalse(os.path.isabs(s), (where, s))
            self.assertNotIn("..", re.split(r"[/\\]", s), (where, s))

        def strings(v):
            if isinstance(v, str):
                yield v
            elif isinstance(v, dict):
                for x in v.values():
                    yield from strings(x)
            elif isinstance(v, list):
                for x in v:
                    yield from strings(x)

        for i in self.run_ids(work):
            for s in strings(_read_toml(os.path.join(work, i, "input.toml"))):
                check_path(s, i)
        for e in self.manifest(work)["runs"]:
            for s in strings(e):
                check_path(s, "manifest")
        with open(os.path.join(work, "commands.txt")) as f:
            text = f.read()
        self.assertNotIn(self.tmp, text)
        self.assertNotIn(os.path.realpath(self.tmp), text)
        self.assertNotIn("..", text)
        moved = self.path("b", "moved")
        os.makedirs(self.path("b"))
        shutil.move(work, moved)
        shutil.rmtree(src)
        res = self.bash(os.path.join(moved, "commands.txt"))
        self.assertEqual(res.returncode, 0, res.stdout)
        rows, _, _ = self.collect(moved)
        self.assertTrue(rows)
        for r in rows.values():
            self.assertIn(r["status"], ("converged", "not_converged"), r["run_id"])
        target = sorted(rows)[0]
        with open(os.path.join(moved, target, "inputs", "transfer.dat"), "a") as f:
            f.write("\n")
        rows, _, _ = self.collect(moved)
        for rid, r in rows.items():
            if rid == target:
                self.assertEqual(r["status"], "seed_unverified")
                self.assertIn("transfer.dat", r["error_summary"])
            else:
                self.assertIn(r["status"], ("converged", "not_converged"))
        self.assertFalse([n for dp, _, fns in os.walk(moved) for n in fns
                          if n.startswith("hwave-run.")])

    def test_exit_codes(self):
        src = self.model()
        work = self.path("WORK")
        self.prepare(src, work)
        rd = os.path.join(work, self.run_ids(work)[0])
        toml = os.path.join(rd, "input.toml")
        self.assertEqual(self.main("check", toml)[0], 0)
        self.assertEqual(self.main("check", "--require-provenance", toml)[0], 0)
        with np.load(os.path.join(rd, "initial.npz"), allow_pickle=False) as data:
            members = {k: data[k] for k in ("green", "green_sublattice")}
        plain = os.path.join(rd, "plain.npz")
        np.savez(plain, **members)
        rc, _, err = self.main("check", toml, plain)
        self.assertEqual(rc, 0)
        self.assertIn("no provenance", err)
        self.assertEqual(self.main("check", "--require-provenance", toml, plain)[0], 4)
        bad = os.path.join(rd, "bad.npz")
        with open(bad, "wb") as f:
            f.write(b"not an npz")
        self.assertEqual(self.main("check", toml, bad)[0], 1)
        # path resolution: initial resolves under [file.input] path_to_input
        _edit(toml, lambda d: d["file"]["input"].update(path_to_input="inputs"))
        rc, _, err = self.main("check", "--require-provenance", toml)
        self.assertEqual(rc, 1)
        self.assertIn("inputs/initial.npz", err)
        # argument errors
        rc, _, err = self.main("prepare", src, "--out", self.path("W2"), "--subshape",
                               "2,2,1", "--subshape", "2,2,1")
        self.assertEqual(rc, 2)
        self.assertIn("more than once", err)
        self.assertEqual(self.main("prepare", src, "--out", self.path("W2"), "--subshape",
                                   "2,2")[0], 2)
        self.assertEqual(self.main("check")[0], 2)
        self.assertEqual(self.main("bogus")[0], 2)
        self.assertFalse(os.path.exists(self.path("W2")))

    def test_prepare_reentry_and_foreign_work(self):
        src = self.model()
        work = self.path("WORK")
        self.prepare(src, work)
        snap = self._tree(work, None)
        msg = self.prepare(src, work)
        self.assertIn("already prepared", msg)
        self.assertEqual(self._tree(work, None), snap)
        # a different command
        rc, _, err = self.main("prepare", src, "--out", work, "--max-candidates", "3")
        self.assertEqual(rc, 1)
        self.assertIn("command", err)
        # an input changed since
        with open(os.path.join(src, "transfer.dat"), "a") as f:
            f.write("\n")
        from hwave.initial import cli
        with self.assertRaisesRegex(ValueError, "input file Transfer changed"):
            cli.prepare(src, work, n_random=0, max_candidates=0, tier="large",
                        stdout=io.StringIO())
        # a foreign directory
        foreign = self.path("foreign")
        os.makedirs(foreign)
        with open(os.path.join(foreign, "notes.txt"), "w") as f:
            f.write("x")
        rc, _, err = self.main("prepare", src, "--out", foreign)
        self.assertEqual(rc, 1)
        self.assertIn("holds no manifest.json", err)
        self.assertEqual(os.listdir(foreign), ["notes.txt"])
        # a moved WORK re-enters (the command is compared without the WORK path)
        src2 = self.model("src2")
        w2 = self.path("W2")
        self.prepare(src2, w2)
        shutil.move(w2, self.path("W3"))
        self.assertIn("already prepared", self.prepare(src2, self.path("W3")))

    def test_prepare_refuses_outputs_on_run_files(self):
        src = self.model()
        toml = os.path.join(src, "input.toml")
        cases = (({"path_to_output": "", "green": "initial"}, "initial.npz"),
                 ({"path_to_output": "inputs"}, "path_to_output"),
                 ({"path_to_output": ".", "energy": "input.toml"}, "input.toml"),
                 ({"path_to_output": "", "eigen": "inputs/transfer.dat"}, "inputs"),
                 ({"path_to_output": "", "onebodyg": "hwave.log"}, "hwave.log"))
        for i, (out, word) in enumerate(cases):
            _edit(toml, lambda d: d["file"].update(output=dict(out)))
            w = self.path("W{}".format(i))
            rc, _, err = self.main("prepare", src, "--out", w)
            self.assertEqual(rc, 1, out)
            self.assertIn(word, err, out)
            self.assertFalse(os.path.exists(w))
        # outputs in the run directory itself, on other names, are fine
        _edit(toml, lambda d: d["file"].update(output={"path_to_output": "",
                                                       "green": "green_out"}))
        self.prepare(src, self.path("W_ok"))

    def test_prepare_toml_syntax_error_names_file(self):
        src = self.model()
        with open(os.path.join(src, "input.toml"), "a") as f:
            f.write("bad = = 1\n")
        for arg in (src, os.path.join(src, "input.toml")):
            rc, _, err = self.main("prepare", arg, "--out", self.path("W"))
            self.assertEqual(rc, 1)
            self.assertIn(os.path.join(src, "input.toml"), err)
            self.assertIn("not a valid TOML file", err)
            self.assertEqual(len(err.strip().splitlines()), 1, err)

    def test_prepare_extra_candidates(self):
        src = self.model()
        pattern = {"patternspec_version": 1, "q": ["1/2", "0", "0"], "subshape": [2, 1, 1],
                   "phase": 0, "sign": 1, "modulation": "cos",
                   "entries": [{"orb_a": 0, "orb_b": 0, "r": [0, 0, 0], "cz": [1.0, 0.0]}]}
        extra = self.path("extra.json")
        with open(extra, "w") as f:
            json.dump([{"family": "user", "pattern": pattern, "label": "stripe"},
                       {"family": "literature", "pattern": pattern, "label": "stripe",
                        "citation": "a paper"}], f)
        work = self.path("WORK")
        rc, _, err = self.main("prepare", src, "--out", work, "--extra-candidates", extra,
                               "--max-candidates", "0")
        self.assertEqual(rc, 0, err)
        with open(os.path.join(work, "coverage.json")) as f:
            rows = json.load(f)["rows"]
        fams = sorted(r["family"] for r in rows if r["family"] in ("user", "literature"))
        self.assertEqual(fams, ["literature", "user"])
        self.assertIn(os.path.realpath(extra), self.manifest(work)["input"]["files_sha256"])
        with open(extra, "w") as f:
            json.dump([{"family": "user", "pattern": pattern, "bogus": 1}], f)
        rc, _, err = self.main("prepare", src, "--out", self.path("W2"),
                               "--extra-candidates", extra)
        self.assertEqual(rc, 1)
        self.assertIn("bogus", err)

    def test_hand_edited_subshape_caught(self):
        src = self.model()
        work = self.path("WORK")
        self.prepare(src, work)
        rid = self.run_ids(work)[-1]
        toml = os.path.join(work, rid, "input.toml")
        for B in ([4, 1, 1], [1, 1, 1]):
            _edit(toml, lambda d: d["mode"]["param"].update(SubShape=B))
            rc, _, err = self.main("check", "--require-provenance", toml)
            self.assertEqual(rc, 1, B)
            self.assertIn("subshape", err)
            rows, _, _ = self.collect(work)
            self.assertEqual(rows[rid]["status"], "seed_unverified")


class TestCollect(_Base):
    def prepared(self, **param):
        # T = 0.3, U = 2: the paramagnet is stable, so the para run converges to it
        p = {"IterationMax": 300, "EPS": 10}
        p.update(param)
        src = self.model(param=p, U=2.0, T=0.3)
        work = self.path("WORK")
        self.prepare(src, work)
        return work, self.run_ids(work)

    def para(self, ids):
        return [i for i in ids if i.startswith("para-")][0]

    def test_collect_status_converged(self):
        work, ids = self.prepared()
        rid = self.para(ids)
        self.assertEqual(self.run_hwave(os.path.join(work, rid)), 0)
        # a run record naming another directory is not consulted: collect
        # reads every run with its run directory as the base directory
        from hwave.results import write_run_record
        write_run_record(os.path.join(work, rid, "input.toml"), cwd=self.tmp, argv=["hwave"])
        rows, doc, _ = self.collect(work)
        r = rows[rid]
        self.assertEqual(r["status"], "converged")
        self.assertEqual(r["run_dir"], rid)
        self.assertGreater(r["iterations"], 0)
        self.assertEqual((r["iteration_max"], r["eps"], r["mix"]), (300, 10, 0.5))
        self.assertEqual(r["scf_eps_tol"], 1e-10)
        self.assertIsNotNone(r["free_energy_per_cell"])
        self.assertIsNotNone(r["order_parameters"]["per_orbital"])
        self.assertEqual(r["final_label"], "para")
        self.assertEqual(r["h_used"], 0)
        self.assertEqual(r["seed_base"], "bare")
        self.assertEqual(doc["phase_marker"], "post_run")
        self.assertEqual(doc["ranking"]["lists"][0]["entries"][0]["run_id"], rid)
        for other in ids:
            if other != rid:
                self.assertEqual(rows[other]["status"], "not_run")
        with open(os.path.join(work, "summary.txt")) as f:
            self.assertIn(rid, f.read())

    def test_collect_status_not_converged(self):
        work, ids = self.prepared(IterationMax=2)
        rid = [i for i in ids if not i.startswith("para-")][0]
        self.run_hwave(os.path.join(work, rid))
        r = self.collect(work)[0][rid]
        self.assertEqual(r["status"], "not_converged")
        self.assertEqual(r["iterations"], 2)

    def test_collect_status_zero_steps(self):
        work, ids = self.prepared(IterationMax=0)
        rid = self.para(ids)
        self.run_hwave(os.path.join(work, rid))
        r = self.collect(work)[0][rid]
        self.assertEqual(r["status"], "failed")
        self.assertEqual(r["failure_kind"], "zero_steps")
        self.assertEqual(r["error_summary"], "no SCF step performed")
        self.assertEqual(r["iterations"], 0)

    def _no_status_file(self):
        src = self.model(param={"IterationMax": 300, "EPS": 10}, U=2.0, T=0.3)
        _edit(os.path.join(src, "input.toml"),
              lambda d: d["file"]["output"].update(scf_status=""))
        work = self.path("WORK")
        self.prepare(src, work)
        return work, self.para(self.run_ids(work))

    def test_collect_status_no_status_file_with_log(self):
        work, rid = self._no_status_file()
        self.run_hwave(os.path.join(work, rid))
        r = self.collect(work)[0][rid]
        self.assertEqual(r["status"], "converged")
        self.assertIsNone(r["iterations"])
        self.assertLessEqual(r["scf_residual"], 1e-10)

    def test_collect_status_no_status_file_without_log(self):
        work, rid = self._no_status_file()
        self.run_hwave(os.path.join(work, rid), log="scheduler.out")
        r = self.collect(work)[0][rid]
        self.assertEqual(r["status"], "convergence_unverified")
        self.assertIsNotNone(r["free_energy_per_cell"])

    def test_collect_status_log_rules(self):
        work, rid = self._no_status_file()
        rd = os.path.join(work, rid)
        self.run_hwave(rd)
        log = os.path.join(rd, "hwave.log")
        with open(log) as f:
            good = f.read()
        cases = (
            (re.sub(r"succeeded: rest=\S+, eps=", "succeeded: rest=0.001, eps=", good),
             "not_converged", "exceeds the tolerance"),
            (good + "Traceback (most recent call last):\n  File x\nKeyError: 'boom'\n",
             "failed", "KeyError: 'boom'"),
            (good + "2026-01-01 00:00:00,0 ERROR qlms.uhfk: abort\n", "failed", "abort"),
            ("".join(ln for ln in good.splitlines(True) if "succeeded" not in ln),
             "not_converged", "no convergence message"),
        )
        for text, status, summary in cases:
            with open(log, "w") as f:
                f.write(text)
            r = self.collect(work)[0][rid]
            self.assertEqual(r["status"], status, text[-200:])
            if summary is not None:
                self.assertIn(summary, r["error_summary"])
            if status == "failed":
                self.assertEqual(r["failure_kind"], "crash")

    def test_collect_status_corrupt_status_file(self):
        work, ids = self.prepared()
        rid = self.para(ids)
        rd = os.path.join(work, rid)
        self.run_hwave(rd)
        with open(os.path.join(rd, "output", "scf_status.dat"), "w") as f:
            f.write("format_version = 2\n")
        rows, doc, err = self.collect(work)
        r = rows[rid]
        self.assertEqual(r["status"], "failed")
        self.assertEqual(r["failure_kind"], "unreadable_output")
        self.assertTrue(r["error_summary"].startswith("scf_status.dat unreadable: "),
                        r["error_summary"])
        self.assertIn("format_version", r["error_summary"])
        self.assertIn("scf_status.dat unreadable", err)
        self.assertIsNone(r["free_energy_per_cell"])
        ranked = [e["run_id"] for lst in doc["ranking"]["lists"] for e in lst["entries"]]
        self.assertNotIn(rid, ranked)
        self.assertIn(rid, [u["run_id"] for u in doc["ranking"]["unranked"]])
        # no absolute path in the collected result
        for text in (r["error_summary"], err):
            self.assertNotIn(self.tmp, text)
            self.assertNotIn(os.path.realpath(self.tmp), text)
        self.assertIn("output/scf_status.dat", r["error_summary"])

    def test_collect_final_label_carries_q(self):
        # U = 8 at half filling: the Neel seed (spin z at q = (1/2, 1/2))
        # converges to the Neel state, the ferro seed stays ferro; the names
        # must tell them apart. The para seed is not used: it is exactly
        # spin-symmetric (G_up = G_dn, no up-down block), a subspace the
        # Hartree-Fock map preserves, so whether it ever reaches the Neel
        # state depends on round-off alone (measured: macOS Accelerate in
        # 35 iterations; Linux OpenBLAS stays spin-symmetric and stalls at a
        # residual of about 1e-2 until IterationMax).
        src = self.model(U=8.0, param={"IterationMax": 2000, "EPS": 8})
        work = self.path("WORK")
        self.prepare(src, work, max_candidates=None)
        ids = self.run_ids(work)
        neel_id = [i for i in ids if i.startswith("baseline-q1_2_1_2_0-spin-z-")]
        ferro_id = [i for i in ids if i.startswith("ferro-")]
        self.assertEqual((len(neel_id), len(ferro_id)), (1, 1), ids)
        for i in neel_id + ferro_id:
            self.run_hwave(os.path.join(work, i))
        rows = self.collect(work)[0]
        neel, ferro = rows[neel_id[0]], rows[ferro_id[0]]
        self.assertEqual((neel["status"], ferro["status"]), ("converged", "converged"))
        self.assertLess(neel["free_energy_per_cell"], ferro["free_energy_per_cell"] - 0.1)
        self.assertIn("q=(1/2,1/2,0)", neel["final_label"])
        self.assertIn("q=(0,0,0)", ferro["final_label"])
        self.assertNotEqual(neel["final_label"], ferro["final_label"])
        # q follows the on-site spin channel it belongs to, never a bond one
        for r in (neel, ferro):
            self.assertRegex(r["final_label"], r"^spin-z \S+-orbital q=\(", r["final_label"])
            self.assertNotRegex(r["final_label"], r"bond-[a-z-]+ q=", r["final_label"])

    def test_dominant_q_near_tie(self):
        from hwave.initial.cli import _dominant_q
        op = {"q": ["0,0,0", "1/2,0,0", "0,1/2,0", "1/2,1/2,0"],
              "spin_amplitude": [0.0, 0.30, 0.29, 0.01],
              "charge_amplitude": [1.0, 0.0, 0.0, 0.0]}
        self.assertEqual(_dominant_q(op), ["1/2,0,0", "0,1/2,0"])
        op["spin_amplitude"] = [0.0, 0.3, 0.1, 0.0]
        self.assertEqual(_dominant_q(op), ["1/2,0,0"])
        op["spin_amplitude"] = [0.0] * 4
        self.assertEqual(_dominant_q(op), [])

    def test_collect_broken_manifest_one_line(self):
        work, ids = self.prepared()
        m = self.manifest(work)
        for breaker, word in ((lambda m: m.pop("runs"), "runs"),
                              (lambda m: m["runs"][0].pop("input_sha256"), "input_sha256"),
                              (lambda m: m["runs"][0].update(run_dir="/abs"), "/abs")):
            bad = json.loads(json.dumps(m))
            breaker(bad)
            with open(os.path.join(work, "manifest.json"), "w") as f:
                json.dump(bad, f)
            rc, _, err = self.main("collect", work)
            self.assertEqual(rc, 1)
            self.assertIn("manifest.json", err)
            self.assertIn(word, err)
            self.assertNotIn("Traceback", err)
            self.assertEqual(len(err.strip().splitlines()), 1, err)

    def test_collect_duplicate_run_id_refused(self):
        work, ids = self.prepared()
        m = self.manifest(work)
        m["runs"].append(dict(m["runs"][0]))
        with open(os.path.join(work, "manifest.json"), "w") as f:
            json.dump(m, f)
        rc, _, err = self.main("collect", work)
        self.assertEqual(rc, 1)
        self.assertIn("duplicate run_id", err)
        self.assertIn(m["runs"][0]["run_id"], err)
        self.assertNotIn("Traceback", err)

    def test_collect_refuses_symlinked_outputs(self):
        # an output name in WORK that is a symbolic link to a file outside
        # WORK: collect refuses it and never writes through the link
        work, ids = self.prepared()
        self.collect(work)
        victim = self.path("victim.txt")
        for name in ("summary.txt", "coverage.json"):
            with open(victim, "w") as f:
                f.write("keep me\n")
            dest = os.path.join(work, name)
            os.remove(dest)
            os.symlink(victim, dest)
            rc, _, err = self.main("collect", work)
            self.assertEqual(rc, 1, name)
            self.assertIn("symbolic link", err)
            self.assertIn(name, err)
            self.assertNotIn("Traceback", err)
            with open(victim) as f:
                self.assertEqual(f.read(), "keep me\n")
            self.assertTrue(os.path.islink(dest))
            os.remove(dest)
        # a link at the old fixed temporary name is never opened either
        os.symlink(victim, os.path.join(work, "coverage.json.tmp"))
        self.collect(work)
        with open(victim) as f:
            self.assertEqual(f.read(), "keep me\n")
        self.assertFalse(os.path.islink(os.path.join(work, "coverage.json")))
        self.assertEqual([n for n in os.listdir(work) if n.startswith(".")], [])

    def test_write_text_replaces_atomically(self):
        from hwave.initial import cli
        dest = self.path("out.txt")
        cli._write_text(dest, "a\n")
        cli._write_text(dest, "b\n")
        with open(dest) as f:
            self.assertEqual(f.read(), "b\n")
        self.assertEqual(sorted(os.listdir(self.tmp)), ["_bin", "out.txt"])
        dangling = self.path("dangling.txt")
        os.symlink(self.path("nowhere"), dangling)
        with self.assertRaisesRegex(ValueError, "symbolic link"):
            cli._write_text(dangling, "x")
        self.assertFalse(os.path.exists(self.path("nowhere")))

    def test_collect_status_input_changed(self):
        work, ids = self.prepared()
        rid = self.para(ids)
        rd = os.path.join(work, rid)
        self.run_hwave(rd)
        with open(os.path.join(rd, "inputs", "coulombintra.dat"), "a") as f:
            f.write("\n")
        rows = self.collect(work)[0]
        self.assertEqual(rows[rid]["status"], "seed_unverified")
        self.assertIn("coulombintra.dat", rows[rid]["error_summary"])
        # the initial file replaced by hand
        rid2 = [i for i in ids if i != rid][0]
        shutil.copyfile(os.path.join(rd, "initial.npz"),
                        os.path.join(work, rid2, "initial.npz"))
        rows = self.collect(work)[0]
        self.assertEqual(rows[rid2]["status"], "seed_unverified")
        self.assertIn("initial.npz", rows[rid2]["error_summary"])

    def test_collect_status_never_run(self):
        work, ids = self.prepared()
        rows = self.collect(work)[0]
        self.assertEqual({r["status"] for r in rows.values()}, {"not_run"})

    def test_collect_idempotent(self):
        work, ids = self.prepared()
        self.run_hwave(os.path.join(work, self.para(ids)))
        self.collect(work)
        a = self._read(work)
        self.collect(work)
        self.assertEqual(self._read(work), a)

    def _read(self, work):
        out = {}
        for n in ("coverage.json", "summary.txt", "manifest.json", "candidates.json"):
            with open(os.path.join(work, n), "rb") as f:
                out[n] = f.read()
        return out

    def test_collect_continues_after_read_error(self):
        work, ids = self.prepared()
        self.assertGreaterEqual(len(ids), 2)
        for i in ids:
            self.run_hwave(os.path.join(work, i))
        bad = self.para(ids)
        with open(os.path.join(work, bad, "output", "green.npz"), "wb") as f:
            f.write(b"corrupt")
        rows = self.collect(work)[0]
        self.assertEqual(rows[bad]["status"], "failed")
        self.assertEqual(rows[bad]["failure_kind"], "unreadable_output")
        self.assertIn("green", rows[bad]["error_summary"])
        for i in ids:
            if i != bad:
                self.assertIn(rows[i]["status"], ("converged", "not_converged"))
                self.assertIsNone(rows[i]["error_summary"])


class TestHelp(_Base):
    def test_cli_help_has_no_internal_terms(self):
        texts = []
        for argv in (["--help"], ["prepare", "--help"], ["check", "--help"],
                     ["collect", "--help"]):
            rc, out, err = self.main(*argv)
            self.assertEqual(rc, 0, argv)
            texts.append(out)
        for t in texts:
            self.assertIsNone(_INTERNAL.search(t), _INTERNAL.search(t))
            for code in ("0", "1", "2", "3", "4"):
                self.assertRegex(t, r"\n  {}  ".format(code))
        res = subprocess.run([sys.executable, "-m", "hwave.initial", "--help"],
                             env=self.env, stdout=subprocess.PIPE, text=True)
        self.assertEqual(res.returncode, 0)
        self.assertIn("prepare", res.stdout)

    def test_no_profile_table_in_command_output(self):
        src = self.model()
        work = self.path("WORK")
        self.prepare(src, work)
        toml = os.path.join(work, self.run_ids(work)[0], "input.toml")
        for argv in (["check", "--require-provenance", toml], ["collect", work]):
            res = subprocess.run([sys.executable, "-m", "hwave.initial"] + argv,
                                 env=self.env, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, text=True)
            self.assertEqual(res.returncode, 0, res.stdout)
            self.assertNotIn("Statistics", res.stdout)
            self.assertNotIn("msec", res.stdout)
        # the process exits normally: exit handlers run
        code = ("import atexit, runpy, sys; atexit.register(lambda: print('exit handler ran')); "
                "sys.argv = ['hwave.initial', 'check', '--require-provenance', {!r}]; "
                "runpy.run_module('hwave.initial', run_name='__main__')").format(toml)
        res = subprocess.run([sys.executable, "-c", code], env=self.env,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        self.assertEqual(res.returncode, 0, res.stdout)
        self.assertIn("exit handler ran", res.stdout)
        self.assertNotIn("Statistics", res.stdout)

    def test_toml_serialiser(self):
        import tomli
        from hwave.initial.cli import dumps_toml
        doc = {"a": 1, "b": {"c": "x\"\\\n\x7fé", "d": [1, 2.5, True, "s"],
                             "2Sz": 0, "e": {"f": -0.0}}, "g": {}}
        self.assertEqual(tomli.loads(dumps_toml(doc)), doc)
        for bad, key in (({"t": {"x": [[1]]}}, "x"), ({"t": {"y": float("nan")}}, "y"),
                         ({"t": {"z": None}}, "z"), ({"t": {"w": [{"a": 1}]}}, "w")):
            with self.assertRaisesRegex(ValueError, key):
                dumps_toml(bad)


if __name__ == "__main__":
    unittest.main()
