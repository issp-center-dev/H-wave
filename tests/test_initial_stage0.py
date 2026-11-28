"""``prepare --stage0``: the paramagnetic Hartree-Fock step of the small tier
run as separate run directories (the preliminary paramagnetic run), and the
re-entry of the same ``prepare`` command after it ran.

Every test writes under ``tempfile.mkdtemp()``. The solver runs as
``bash WORK/stage0_commands.txt`` runs it, through the shims of
``tests/test_initial_cli.py``.
"""
import io
import json
import os
import re
import shutil
import unittest
from unittest import mock

import numpy as np

from tests.heavy_tests import heavy
from tests.test_initial_cli import (_Base, _read_toml, _sha, check_refused_case,
                                    staging_cases)

STAGE0 = ("stage0_para", "stage0_para_check")
PENDING_TEXT = "re-run the same prepare command after these finish"
_KW = {"n_random": 0, "max_candidates": 0, "tier": "small"}


def _strip(doc):
    """``doc`` without ``[mode.param] SubShape``, ``T``, ``EPS`` and
    ``[file.output] green``: the keys the preliminary paramagnetic run sets
    differently from the candidate runs."""
    doc = json.loads(json.dumps(doc))
    for k in list(doc["mode"]["param"]):
        if k.lower() in ("subshape", "t", "eps"):
            del doc["mode"]["param"][k]
    doc["file"]["output"].pop("green", None)
    return doc


class _Stage0(_Base):
    def prep(self, src, work, *, para_check=False, **kw):
        """``cli.prepare(..., stage0=True)``; ``(exit code, stdout)``."""
        from hwave.initial import cli
        args = dict(_KW, **kw)
        out = io.StringIO()
        rc = cli.prepare(src, work, stdout=out, stage0=True, para_check=para_check, **args)
        return rc, out.getvalue()

    def src(self, name="src", **param):
        p = {"IterationMax": 200}
        p.update(param)
        return self.model(name, param=p)

    def coverage(self, work):
        with open(os.path.join(work, "coverage.json")) as f:
            return json.load(f)

    def tree(self, root):
        out = {}
        for dp, _dns, fns in os.walk(root):
            for fn in fns:
                p = os.path.join(dp, fn)
                out[p] = _sha(p)
        return out

    def run_stage0(self, work, names=STAGE0):
        """Run hwave in the named directories, as stage0_commands.txt does."""
        for n in names:
            if os.path.isdir(os.path.join(work, n)):
                self.assertEqual(self.run_hwave(os.path.join(work, n)), 0, n)

    def status_path(self, work, name="stage0_para"):
        return os.path.join(work, name, "output", "scf_status.dat")

    def edit_status(self, work, name, **values):
        path = self.status_path(work, name)
        with open(path) as f:
            lines = f.read().splitlines()
        out = []
        for ln in lines:
            key = ln.split("=")[0].strip()
            out.append("{} = {}".format(key, values[key]) if key in values else ln)
        with open(path, "w") as f:
            f.write("\n".join(out) + "\n")

    def ran_and_failed(self, work, name, text="RuntimeError: the solver stopped"):
        """``name`` as a run that crashed: no outputs, a traceback in its
        log."""
        rd = os.path.join(work, name)
        shutil.rmtree(os.path.join(rd, "output"), ignore_errors=True)
        with open(os.path.join(rd, "hwave.log"), "w") as f:
            f.write("Traceback (most recent call last):\n  File \"x\", line 1\n" + text + "\n")


class TestStage0Pending(_Stage0):
    def test_stage0_writes_pending(self):
        """The first call writes the preliminary paramagnetic run directories
        (SubShape [1,1,1], T = T_hint, EPS = max(EPS, 12), green forced),
        their start files with provenance, stage0_commands.txt and the
        manifest (stage para_pending, SHA-256 records), and exits 3."""
        from hwave.initial import cli, hints
        from hwave.initial.facts import model_facts
        from hwave.initial.writer import check_initial_state
        src = self.src()
        work = self.path("WORK")
        rc, out = self.prep(src, work, para_check=True)
        self.assertEqual(rc, cli.EXIT_STAGE0_PENDING)
        self.assertIn(PENDING_TEXT, out)
        self.assertIn("stage0_commands.txt", out)
        self.assertNotRegex(out, r"\bstage 0\b")
        self.assertEqual(sorted(os.listdir(work)),
                         ["manifest.json", "stage0_commands.txt"] + list(STAGE0))
        facts = model_facts(os.path.join(src, "input.toml"), subshape=(1, 1, 1))
        th = hints.t_hint_bare(facts.problem)
        starts = hints.para_starts(facts, th)
        m = self.manifest(work)
        self.assertEqual(m["stage"], "para_pending")
        self.assertIs(m["command"]["stage0"], True)
        self.assertEqual(m["runs"], [])
        entries = m["stage0"]["runs"]
        self.assertEqual([e["run_dir"] for e in entries], list(STAGE0))
        self.assertEqual(m["stage0"]["T"], th.T_hint)
        for e, start, cid in zip(entries, ("main", "check"),
                                 ("para-hf-start", "para-hf-check-start")):
            rd = os.path.join(work, e["run_dir"])
            self.assertEqual(e["start"], start)
            self.assertEqual(e["candidate_id"], cid)
            self.assertEqual(_sha(os.path.join(rd, "input.toml")), e["input_toml_sha256"])
            self.assertEqual(_sha(os.path.join(rd, "initial.npz")), e["initial_sha256"])
            for s in e["staged"]:
                self.assertEqual(_sha(os.path.join(rd, s["name"])), s["sha256"])
            doc = _read_toml(os.path.join(rd, "input.toml"))
            p = doc["mode"]["param"]
            self.assertEqual(p["SubShape"], [1, 1, 1])
            self.assertEqual(p["T"], th.T_hint)
            self.assertEqual(p["EPS"], 12)
            self.assertEqual(doc["file"]["output"]["green"], "green")
            self.assertEqual(doc["file"]["input"]["initial"], "initial.npz")
            res = check_initial_state(os.path.join(rd, "input.toml"), require_provenance=True)
            self.assertTrue(res.verified)
            with np.load(os.path.join(rd, "initial.npz"), allow_pickle=False) as d:
                self.assertEqual(json.loads(str(d["candidate_id"][()])), cid)
                np.testing.assert_array_equal(d["green"], starts[start])
        with open(os.path.join(work, "stage0_commands.txt")) as f:
            lines = f.read().splitlines()
        self.assertEqual(lines[0], 'cd "$(dirname "$0")" || exit 1')
        self.assertEqual([ln for ln in lines if ln.startswith("(")],
                         ['(cd "{}" && python -m hwave.initial check --require-provenance '
                          'input.toml && hwave input.toml > hwave.log 2>&1)'.format(n)
                          for n in STAGE0])
        # without --para-check: one directory
        w2 = self.path("W2")
        self.assertEqual(self.prep(src, w2)[0], cli.EXIT_STAGE0_PENDING)
        self.assertEqual([e["run_dir"] for e in self.manifest(w2)["stage0"]["runs"]],
                         ["stage0_para"])
        self.assertFalse(os.path.exists(os.path.join(w2, "stage0_para_check")))

    def test_stage0_refused_with_large(self):
        src = self.src()
        rc, _, err = self.main("prepare", src, "--out", self.path("W1"), "--tier", "large",
                               "--stage0")
        self.assertEqual(rc, 1)
        self.assertIn("--stage0", err)
        self.assertFalse(os.path.exists(self.path("W1")))
        # a SubShape the enumeration refuses is refused before anything is written
        rc, _, err = self.main("prepare", src, "--out", self.path("W3"), "--tier", "small",
                               "--stage0", "--subshape", "3,1,1")
        self.assertEqual(rc, 1)
        self.assertFalse(os.path.exists(self.path("W3")))
        # auto choosing the large tier is refused too, before anything is written
        from hwave.initial import enumerate as en
        with mock.patch.object(en, "SMALL_TIER_MAX_PAIRS", 0):
            rc, _, err = self.main("prepare", src, "--out", self.path("W2"), "--tier", "auto",
                                   "--stage0")
        self.assertEqual(rc, 1)
        self.assertIn("large", err)
        self.assertFalse(os.path.exists(self.path("W2")))

    def test_stage0_command_compared(self):
        """--stage0 is part of the command: the same WORK prepared without it
        is refused, naming it; collect refuses a pending WORK."""
        from hwave.initial import cli
        src = self.src()
        work = self.path("WORK")
        self.assertEqual(self.prep(src, work)[0], 3)
        with self.assertRaisesRegex(ValueError, "command stage0"):
            cli.prepare(src, work, stdout=io.StringIO(), **_KW)
        with self.assertRaisesRegex(ValueError, "preliminary paramagnetic run"):
            cli.collect(work, stdout=io.StringIO(), stderr=io.StringIO())

    def test_stage0_states_not_yet_run(self):
        """Not yet run: the commands are printed again, exit 3, nothing
        changes."""
        src = self.src()
        work = self.path("WORK")
        self.assertEqual(self.prep(src, work, para_check=True)[0], 3)
        snap = self.tree(work)
        rc, out = self.prep(src, work, para_check=True)
        self.assertEqual(rc, 3)
        self.assertIn(PENDING_TEXT, out)
        self.assertEqual(self.tree(work), snap)
        # one of the two run: still pending
        self.run_stage0(work, ["stage0_para"])
        snap = self.tree(work)
        self.assertEqual(self.prep(src, work, para_check=True)[0], 3)
        self.assertEqual(self.tree(work), snap)

    def test_stage0_states_sha_mismatch(self):
        """A stage-0 input that differs from the manifest is refused, naming
        the file."""
        from hwave.initial import cli
        src = self.src()
        work = self.path("WORK")
        self.assertEqual(self.prep(src, work)[0], 3)
        for rel in ("input.toml", "initial.npz", os.path.join("inputs", "transfer.dat")):
            with self.subTest(file=rel):
                path = os.path.join(work, "stage0_para", rel)
                with open(path, "rb") as f:
                    keep = f.read()
                with open(path, "ab") as f:
                    f.write(b"\n")
                with self.assertRaises(ValueError) as cm:
                    self.prep(src, work)
                self.assertIn(os.path.join("stage0_para", rel), str(cm.exception))
                with open(path, "wb") as f:
                    f.write(keep)
        os.remove(os.path.join(work, "stage0_para", "initial.npz"))
        with self.assertRaisesRegex(ValueError, "initial.npz"):
            self.prep(src, work)
        self.assertFalse(os.path.exists(os.path.join(work, "commands.txt")))
        self.assertIsInstance(cli.EXIT_STAGE0_PENDING, int)

    def test_stage0_states_ran_and_failed(self):
        """A traceback without a status file: the bare fallback with the
        reason, exit 0; the rows are seeded on the bare H0."""
        src = self.src()
        work = self.path("WORK")
        self.assertEqual(self.prep(src, work)[0], 3)
        self.ran_and_failed(work, "stage0_para")
        # a re-entry that fails is rolled back: WORK stays pending
        snap = self.tree(work)
        from hwave.initial import enumerate as en
        with mock.patch.object(en, "candidate_states", side_effect=RuntimeError("boom")):
            with self.assertRaisesRegex(RuntimeError, "boom"):
                self.prep(src, work)
        self.assertEqual(self.tree(work), snap)
        rc, out = self.prep(src, work)
        self.assertEqual(rc, 0, out)
        cov = self.coverage(work)
        self.assertEqual(cov["hint"]["status"], "skipped")
        self.assertIn("traceback", cov["hint"]["skip_reason"])
        self.assertIn("preliminary paramagnetic run", cov["hint"]["skip_reason"])
        tried = [r for r in cov["rows"] if r["tried"]]
        self.assertTrue(tried)
        self.assertEqual({r["seed_base"] for r in tried}, {"bare"})
        m = self.manifest(work)
        self.assertEqual(m["stage"], "done")
        self.assertEqual(m["stage0"]["runs"][0]["outcome"], "failed")
        self.assertIsNone(m["para"]["green"])
        self.assertTrue(m["runs"])
        self.assertFalse([e for e in m["runs"] if e["run_dir"].startswith("stage0_")])
        # done: the same command now has nothing to do
        rc, out = self.prep(src, work)
        self.assertEqual(rc, 0)
        self.assertIn("nothing to do", out)


class TestStage0Runs(_Stage0):
    """Branches that run the solver on the preliminary paramagnetic run."""

    @heavy
    def test_stage0_states(self):
        """A valid status file decides (converged with a readable green ->
        converged; converged = false, an unreadable status file or green ->
        ran and failed); without it the log decides."""
        from hwave.initial import cli
        src = self.src()
        base = self.path("W0")
        self.assertEqual(self.prep(src, base)[0], 3)
        self.run_stage0(base)
        cases = {}

        def variant(name, fn):
            w = self.path("W_" + name)
            shutil.copytree(base, w)
            fn(w)
            cases[name] = w

        variant("converged", lambda w: None)
        variant("not_converged", lambda w: self.edit_status(w, "stage0_para",
                                                            converged="false"))
        variant("zero_steps", lambda w: self.edit_status(
            w, "stage0_para", converged="false", iterations="0", residual="nan",
            note="no SCF step performed"))
        variant("unreadable_status", lambda w: self.edit_status(
            w, "stage0_para", iterations="many"))
        variant("green_missing", lambda w: os.remove(
            os.path.join(w, "stage0_para", "output", "green.npz")))
        variant("log_only", lambda w: os.remove(self.status_path(w)))
        def scaled_green(w):
            path = os.path.join(w, "stage0_para", "output", "green.npz")
            with np.load(path, allow_pickle=False) as d:
                members = {k: d[k] for k in d.files}
            members["green"] = members["green"] * 1.5
            np.savez(path, **members)
        variant("green_scaled", scaled_green)
        variant("log_ended_badly", lambda w: (os.remove(self.status_path(w)), open(
            os.path.join(w, "stage0_para", "hwave.log"), "a").write(
                "2026-01-01 00:00:00 ERROR qlms: the solver stopped\n")))
        variant("outputs_without_log", lambda w: (
            os.remove(self.status_path(w)),
            os.remove(os.path.join(w, "stage0_para", "hwave.log"))))
        want = {"converged": ("computed", None),
                "not_converged": ("skipped", "did not converge"),
                "zero_steps": ("skipped", "no SCF step performed"),
                "unreadable_status": ("skipped", "unreadable"),
                "green_missing": ("skipped", "cannot be read"),
                "log_only": ("computed", None),
                "green_scaled": ("skipped", "not a state of this input"),
                "log_ended_badly": ("skipped", "solver error"),
                "outputs_without_log": ("skipped", "hwave.log")}
        for name, w in cases.items():
            with self.subTest(case=name):
                rc, out = self.prep(src, w)
                self.assertEqual(rc, cli.EXIT_OK, out)
                hint = self.coverage(w)["hint"]
                status, word = want[name]
                self.assertEqual(hint["status"], status, hint["skip_reason"])
                if word is not None:
                    self.assertIn(word, hint["skip_reason"])
                if name == "zero_steps":
                    # the status file says why; no output is read
                    self.assertNotIn("cannot be read", hint["skip_reason"])
                m = self.manifest(w)
                self.assertEqual(m["stage"], "done")
                outcome = m["stage0"]["runs"][0]["outcome"]
                self.assertEqual(outcome, "converged" if status == "computed" else
                                 "not_converged" if name == "not_converged" else "failed")
                rows = [r for r in self.coverage(w)["rows"] if r["tried"]]
                self.assertEqual({r["seed_base"] for r in rows},
                                 {"para_hf" if status == "computed" else "bare"})

    @heavy
    def test_stage0_states_not_finished(self):
        """No status file and no final line of the solver in hwave.log (a run
        still in progress): not finished, exit 3, nothing changes; once its
        final line is there, it is decided."""
        src = self.src()
        work = self.path("WORK")
        self.assertEqual(self.prep(src, work, para_check=True)[0], 3)
        self.run_stage0(work)
        # stage0_para finished; stage0_para_check is "still running": a
        # partial log, no outputs
        rd = os.path.join(work, "stage0_para_check")
        shutil.rmtree(os.path.join(rd, "output"))
        with open(os.path.join(rd, "hwave.log")) as f:
            head = [ln for ln in f.read().splitlines()
                    if "calculation succeeded" not in ln and "calculation failed" not in ln]
        with open(os.path.join(rd, "hwave.log"), "w") as f:
            f.write("\n".join(head[:5]) + "\n")
        snap = self.tree(work)
        rc, out = self.prep(src, work, para_check=True)
        self.assertEqual(rc, 3, out)
        self.assertIn(PENDING_TEXT, out)
        self.assertEqual(self.tree(work), snap)
        # the same with outputs present but no final line
        w2 = self.path("W2")
        shutil.copytree(work, w2)
        os.remove(self.status_path(w2))
        with open(os.path.join(w2, "stage0_para", "hwave.log"), "w") as f:
            f.write("2026-01-01 00:00:00 INFO qlms: step 1\n")
        self.assertEqual(self.prep(src, w2, para_check=True)[0], 3)
        # the run finishes: decided, exit 0
        self.run_stage0(work, ["stage0_para_check"])
        rc, out = self.prep(src, work, para_check=True)
        self.assertEqual(rc, 0, out)
        self.assertEqual(self.coverage(work)["hint"]["status"], "computed")

    @heavy
    def test_stage0_para_check_split(self):
        """With --para-check, when exactly one of the two converges, re-entry
        continues with it and warns."""
        from hwave.initial import hints
        src = self.src()
        base = self.path("W0")
        self.assertEqual(self.prep(src, base, para_check=True)[0], 3)
        self.run_stage0(base)
        for failed, kept, warning in (
                ("stage0_para_check", "main", hints.PARA_CHECK_FAILED_WARNING),
                ("stage0_para", "check", hints.PARA_CHECK_SWITCHED_WARNING),
                (None, "main", None)):
            with self.subTest(failed=failed):
                w = self.path("W_{}".format(failed))
                shutil.copytree(base, w)
                if failed is not None:
                    self.ran_and_failed(w, failed)
                rc, out = self.prep(src, w, para_check=True)
                self.assertEqual(rc, 0, out)
                cov = self.coverage(w)
                self.assertEqual(cov["hint"]["status"], "computed")
                with open(os.path.join(w, "para_hf", "para_hf.json")) as f:
                    rec = json.load(f)
                if failed is not None:
                    self.assertEqual(rec["selected_start"], kept)
                if warning is not None:
                    prefix = warning.split("{}")[0]
                    self.assertTrue(any(x.startswith(prefix) for x in rec["warnings"]),
                                    rec["warnings"])
                    self.assertTrue(any(x.startswith(prefix)
                                        for x in cov["model"]["warnings"]))
                else:
                    self.assertIsNotNone(rec["free_energy_check"])
                m = self.manifest(w)
                self.assertEqual(m["stage0"]["selected"], rec["selected_start"])
        # neither converges: the bare fallback
        w = self.path("W_both")
        shutil.copytree(base, w)
        for n in STAGE0:
            self.ran_and_failed(w, n)
        rc, _ = self.prep(src, w, para_check=True)
        self.assertEqual(rc, 0)
        self.assertEqual(self.coverage(w)["hint"]["status"], "skipped")

    @heavy
    def test_stage0_reentry_proceeds(self):
        """After the preliminary paramagnetic run converged, the same prepare
        command continues to the candidate runs (exit 0, not 3), with the
        same para reference as the in-process step; collect lists only the
        candidate runs."""
        from hwave.initial import cli
        src = self.src()
        work = self.path("WORK")
        self.assertEqual(self.prep(src, work, para_check=True)[0], 3)
        res = self.bash(os.path.join(work, "stage0_commands.txt"), cwd=self.bin)
        self.assertEqual(res.returncode, 0, res.stdout)
        rc, out = self.prep(src, work, para_check=True)
        self.assertEqual(rc, 0, out)
        m = self.manifest(work)
        self.assertEqual(m["stage"], "done")
        self.assertEqual(m["para"]["status"], "converged")
        self.assertEqual(m["para"]["green"], "para_hf/green.npz")
        self.assertEqual(_sha(os.path.join(work, "para_hf", "green.npz")),
                         m["para"]["green_sha256"])
        self.assertEqual(_sha(os.path.join(work, "para_hf", "green.npz")),
                         _sha(os.path.join(work, "stage0_para", "output", "green.npz")))
        cov = self.coverage(work)
        self.assertEqual(cov["hint"]["status"], "computed")
        tried = [r for r in cov["rows"] if r["tried"]]
        self.assertEqual({r["seed_base"] for r in tried}, {"para_hf"})
        # the in-process para step gives the same reference
        ref = self.path("W_ref")
        self.assertEqual(cli.prepare(src, ref, stdout=io.StringIO(), para_check=True,
                                     **_KW), 0)
        with open(os.path.join(ref, "para_hf", "para_hf.json")) as f:
            a = json.load(f)
        with open(os.path.join(work, "para_hf", "para_hf.json")) as f:
            b = json.load(f)
        self.assertAlmostEqual(a["free_energy"], b["free_energy"], delta=1e-10)
        self.assertEqual((a["status"], a["selected_start"], a["T"], a["eps"]),
                         (b["status"], b["selected_start"], b["T"], b["eps"]))
        with np.load(os.path.join(ref, "para_hf", "green.npz")) as ga, \
                np.load(os.path.join(work, "para_hf", "green.npz")) as gb:
            np.testing.assert_allclose(ga["green"], gb["green"], atol=1e-10)
        self.assertEqual([r["run_id"] for r in self.coverage(ref)["rows"]],
                         [r["run_id"] for r in cov["rows"]])
        # collect ignores the stage-0 directories
        rows, _, _ = self.collect(work)
        self.assertEqual(sorted(rows), sorted(e["run_id"] for e in m["runs"]))
        self.assertEqual({r["status"] for r in rows.values()}, {"not_run"})
        rc, out = self.prep(src, work, para_check=True)
        self.assertEqual(rc, 0)
        self.assertIn("nothing to do", out)

    def test_stage0_existing_para_hf_refused_twice(self):
        """A WORK/para_hf left by an interrupted re-entry: the re-entry is
        refused, the rollback leaves WORK unchanged (para_hf and the
        manifest included), and the next re-entry is refused the same way."""
        src = self.src()
        work = self.path("WORK")
        self.assertEqual(self.prep(src, work)[0], 3)
        self.run_stage0(work)
        left = os.path.join(work, "para_hf")
        os.mkdir(left)
        with open(os.path.join(left, "para_hf.json"), "w") as f:
            f.write("{}\n")
        before = self.tree(work)
        self.assertIn(os.path.join(left, "para_hf.json"), before)
        for attempt in (1, 2):
            with self.subTest(attempt=attempt):
                with self.assertRaisesRegex(ValueError, "exists before the preliminary "
                                                        "paramagnetic run was read"):
                    self.prep(src, work)
                self.assertEqual(self.tree(work), before)
        self.assertEqual(self.manifest(work)["stage"], "para_pending")

    @heavy
    def test_stage0_work_relocatable(self):
        """Move WORK after the first call, run its stage0_commands.txt there;
        the same command on the moved WORK passes the checks and exits 0;
        then, with the sources deleted, commands.txt and collect work: no
        seed_unverified."""
        src = self.src(IterationMax=3)
        work = self.path("a", "WORK")
        self.assertEqual(self.prep(src, work)[0], 3)
        moved = self.path("b", "moved")
        os.makedirs(self.path("b"))
        shutil.move(work, moved)
        res = self.bash(os.path.join(moved, "stage0_commands.txt"))
        self.assertEqual(res.returncode, 0, res.stdout)
        rc, out = self.prep(src, moved)
        self.assertEqual(rc, 0, out)
        for dp, dns, fns in os.walk(moved):
            for n in dns + fns:
                self.assertFalse(os.path.islink(os.path.join(dp, n)), n)
        with open(os.path.join(moved, "stage0_commands.txt")) as f:
            self.assertNotIn(self.tmp, f.read())
        shutil.rmtree(src)
        res = self.bash(os.path.join(moved, "commands.txt"))
        self.assertEqual(res.returncode, 0, res.stdout)
        rows, _, _ = self.collect(moved)
        self.assertTrue(rows)
        for r in rows.values():
            self.assertNotEqual(r["status"], "seed_unverified", r.get("error_summary"))
            self.assertIn(r["status"], ("converged", "not_converged"), r["run_id"])


class TestStage0Staging(_Stage0):
    """The staging rules of the candidate runs, repeated for the
    preliminary paramagnetic run (the shared cases of test_initial_cli)."""

    def test_stage0_staging_refusals(self):
        for case in STAGING_REFUSED:
            check_refused_case(self, case, "--tier", "small", "--stage0", "--para-check")

    @heavy
    def test_stage0_staging_like_runs(self):
        for name, build, expected in STAGING_ACCEPTED:
            with self.subTest(case=name):
                toml, info = build(self, name)
                work = self.path("W_" + name)
                rc, _ = self.prep(os.path.dirname(toml), work, para_check=True)
                self.assertEqual(rc, 3)
                for n in STAGE0:
                    expected["check"](self, os.path.join(work, n), info)
                res = self.bash(os.path.join(work, "stage0_commands.txt"))
                self.assertEqual(res.returncode, 0, res.stdout)
                for n in STAGE0:
                    with open(os.path.join(work, n, "hwave.log")) as f:
                        self.assertIn("read initial from initial.npz", f.read())
                    if expected["after_run"] is not None:
                        expected["after_run"](self, os.path.join(work, n), info)
                rc, out = self.prep(os.path.dirname(toml), work, para_check=True)
                self.assertEqual(rc, 0, out)
                first = os.path.join(work, self.run_ids(work)[0])
                run_doc = _strip(_read_toml(os.path.join(first, "input.toml")))
                run_inputs = sorted(os.listdir(os.path.join(first, "inputs")))
                for n in STAGE0:
                    rd = os.path.join(work, n)
                    self.assertEqual(_strip(_read_toml(os.path.join(rd, "input.toml"))),
                                     run_doc)
                    self.assertEqual(sorted(os.listdir(os.path.join(rd, "inputs"))),
                                     run_inputs)
                    for fn in run_inputs:
                        p = os.path.join(rd, "inputs", fn)
                        self.assertFalse(os.path.islink(p))
                        self.assertEqual(_sha(p), _sha(os.path.join(first, "inputs", fn)))


STAGING_REFUSED = [c for c in staging_cases() if "refused" in c[2]]
STAGING_ACCEPTED = [c for c in staging_cases() if "refused" not in c[2]]


class TestStage0Text(unittest.TestCase):
    def test_help_names_the_option(self):
        from hwave.initial import cli
        text = cli.build_parser().format_help()
        sub = [a for a in cli.build_parser()._subparsers._group_actions][0]
        prep = sub.choices["prepare"].format_help()
        self.assertIn("--stage0", prep)
        self.assertIn("preliminary paramagnetic run", re.sub(r"\s+", " ", prep))
        self.assertNotRegex(prep + text, r"\bstage 0\b")
        codes = dict(cli.EXIT_CODES)
        self.assertNotIn("never returned", codes[cli.EXIT_STAGE0_PENDING])
        self.assertIn("preliminary paramagnetic run", codes[cli.EXIT_STAGE0_PENDING])


if __name__ == "__main__":
    unittest.main()
