"""The coverage report and the ranking of ``hwave.initial`` (design sections
8.2 and 8.3, the summary template of section 12)."""
import copy
import dataclasses
import logging
import os
import shutil
import tempfile
import unittest

import numpy as np

from tests import initial_fixtures
from tests.initial_reference_capture import dump_toml

#: design section 8.2 (formatting marks removed, internal references in plain words)
EXCLUSIONS = (
    "multi-q superpositions (double-Q, triple-Q, vortex crystals)",
    "sliding-phase variants beyond the two kept phases, including an Umklapp phase "
    "offset not pinned by symmetry",
    "combinations within a degenerate eigenspace other than the canonical basis",
    "rotated cells such as √3×√3, and cells outside subshapes",
    "orbital and bond patterns in the large tier (it has spin and charge patterns only), "
    "and bond patterns other than hints",
    "one amplitude per candidate",
    "hints describe the paramagnet only (static, at T_hint)",
    "spirals constructed only under a passed spin-rotation check",
    "Extern fields (unused by UHFk)",
    "q with a denominator above max_denominator",
)
SENTENCE = ("The selected state is the lowest free energy among the candidates tried. "
            "It is not a proof of the ground state.")

_SEED_ALL = ("seed_base", "shell", "moments_h0", "m0", "m_inf", "m_target", "m_used",
             "h_hi", "h_used", "amplitude_used", "cap_reason")

_SHELL = {"h0": [{"group": "all", "open_shell": False, "shell_degeneracy": None,
                  "shell_fill": None}], "h_used": None}
_MOM = {"spin": [[0.0, 0.0, 0.0]] * 4, "charge": [1.0] * 4}


class _Tmp(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def dest(self, name="m"):
        return os.path.join(self.tmp, name)


def _cset(dest, *, two_sz=None):
    """A real candidate set: 4x4 square, SubShapes 2x2 and 4x1, a small
    quota (quota rows, duplicates, a quota-boundary split) and one user
    row that the 2x2 cell cannot hold (unrepresentable)."""
    from hwave.initial import Candidate, candidate_states
    from hwave.initial.patterns import Entry, PatternSpec
    toml = initial_fixtures.square(dest, 4, U=2.0, T=0.1, subshape=(2, 2, 1),
                                   two_sz=two_sz)
    spec = PatternSpec(q=("1/4", "0", "0"), subshape=(4, 1, 1), phase=0.0, sign=1,
                       modulation="cos", entries=(Entry(0, 0, (0, 0, 0), cz=1.0),))
    user = Candidate.user(pattern=spec, label="quarter")
    return candidate_states(toml, subshapes=[(2, 2, 1), (4, 1, 1)], max_candidates=3,
                            extra=[user], tier="large")


def _with_rows(cset, rows):
    return dataclasses.replace(cset, rows=rows)


def _seeded(row, **extra):
    row = dict(row)
    row.update(seed_base="bare", shell=_SHELL, moments_h0=_MOM, m0=0.0, m_inf=1.0,
               m_target=0.5, m_used=0.5, h_hi=8.0, h_used=4.0, amplitude_used=0.5,
               cap_reason=None)
    row.update(extra)
    return row


def _extra_reason_rows(cset):
    """Rows for every not-tried reason the enumerator does not produce on
    the fixture, with the seed fields the design section 7 table allows."""
    from hwave.initial.seed import BLOCK_REASON, CONSTRAINT_REASON, DEGENERATE_REASON
    base = next(r for r in cset.rows if r["tried"] and r["family"] == "baseline")

    def nt(suffix, reason, **seed):
        r = dict(base)
        r.update(run_id=base["run_id"] + suffix, tried=False, not_tried_reason=reason)
        r.update(seed)
        return r
    return [
        nt("-x1", "no response"),
        nt("-x2", CONSTRAINT_REASON, seed_base="bare", m0=0.25),
        nt("-x3", DEGENERATE_REASON, seed_base="bare", shell=_SHELL, moments_h0=_MOM,
           m0=0.0, m_inf=5e-7),
        nt("-x4", BLOCK_REASON + ": transverse spin in 2Sz mode"),
        nt("-x5", "refused: writer check: not Hermitian", seed_base="bare", m0=0.0,
           m_inf=1.0, h_used=2.0),
    ]


def _full_rows(cset):
    """Every row of ``cset`` with the seed fields filled where a seed was
    made, plus :func:`_extra_reason_rows`."""
    rows = []
    for r in cset.rows:
        r = dict(r)
        if r["tried"]:
            if r["family"] == "para":
                r.update(seed_base="bare", shell=_SHELL, moments_h0=_MOM, h_used=0.0)
            else:
                r = _seeded(r)
        rows.append(r)
    return rows + _extra_reason_rows(cset)


class TestSchema(_Tmp):

    def test_schema_keys_and_nulls(self):
        from hwave.initial import report
        cset = _cset(self.dest())
        rows = _full_rows(cset)
        reasons = {r["not_tried_reason"].split(" ")[0] if r["not_tried_reason"] else None
                   for r in rows}
        for want in ("quota:", "duplicate", "unrepresentable", "no", "refused:"):
            self.assertIn(want, reasons)
        doc = report.coverage_document(_with_rows(cset, rows), phase_marker="pre_run")
        report.validate_coverage(doc)
        self.assertEqual(doc["schema_version"], 4)
        self.assertEqual(doc["phase_marker"], "pre_run")
        for key in ("model", "tier", "hint", "cells", "rows", "limits"):
            self.assertIn(key, doc)
        for key in ("input_sha256", "hwave_version", "hwave_initial_version",
                    "spin_orbital", "norb", "norb_phys", "iflag_fock",
                    "interacting_orbitals", "n_c", "site_groups", "Ncond",
                    "Ncond_per_cell", "T", "2Sz", "cellshape", "su2", "warnings"):
            self.assertIn(key, doc["model"])
        self.assertNotIn("hint", doc["model"])
        self.assertEqual(doc["tier"]["name"], "large")
        self.assertEqual(doc["tier"]["chosen_by"], "auto (interim rule)")
        for key in ("estimate", "budget", "constants"):
            self.assertIn(key, doc["tier"])
        self.assertEqual(doc["hint"]["status"], "skipped")
        self.assertEqual(doc["limits"]["T_hint"]["value"], 0.1)
        self.assertIn("constants", doc["limits"]["T_hint"])
        self.assertEqual(doc["limits"]["label_thresholds"]["dominant_weight"], 0.1)
        self.assertEqual(doc["limits"]["label_thresholds"]["tr_minor_weight"], 0.1)
        self.assertEqual(doc["limits"]["degeneracy_tol"], 1e-4)
        self.assertEqual(len(doc["rows"]), len(rows))
        for r in doc["rows"]:
            for key in report.ROW_FIELDS:
                self.assertIn(key, r)
            self.assertEqual(r["phase_marker"], "pre_run")
            for key in report.POST_RUN_FIELDS:
                self.assertIsNone(r[key], (r["run_id"], key))
            reason = r["not_tried_reason"] or ""
            if reason.startswith(("quota", "duplicate", "unrepresentable")) \
                    or reason in ("no response",) or reason.startswith("refused: block"):
                for key in _SEED_ALL:
                    self.assertIsNone(r[key], (r["run_id"], key))
        # the I1.9 fields are kept
        for key in ("pattern_label", "proposal_note", "duplicate_of", "amplitude_requested"):
            self.assertIn(key, doc["rows"][0])
        self.assertTrue(doc["cells"][0]["skipped_patterns"])
        # the not-tried counts keep every reason apart, duplicates merged
        line = next(l for l in report.summary_text(doc).splitlines()
                    if l.startswith("Rows:"))
        for part in ("1 no response,", "1 no response (constraint)",
                     "1 no response (degenerate with para)", "1 unrepresentable",
                     " duplicate,", "1 refused: writer check: not Hermitian"):
            self.assertIn(part, line)
        self.assertNotIn("duplicate of baseline", line)
        self.assertIn("quota: baseline (duplicate of a row cut by its quota)", line)

    def test_validate_refuses(self):
        from hwave.initial import report
        from hwave.initial.seed import CONSTRAINT_REASON, DEGENERATE_REASON
        cset = _cset(self.dest())
        doc = report.coverage_document(_with_rows(cset, _full_rows(cset)),
                                       phase_marker="pre_run")

        def bad(fn, what):
            d = copy.deepcopy(doc)
            fn(d)
            with self.assertRaises(ValueError, msg=what) as cm:
                report.validate_coverage(d)
            return str(cm.exception)

        def row(d, pred):
            return next(r for r in d["rows"] if pred(r))

        quota = lambda r: (r["not_tried_reason"] or "").startswith("quota")
        self.assertIn("m0", bad(lambda d: row(d, quota).update(m0=0.1), "quota m0"))
        bad(lambda d: row(d, lambda r: r["not_tried_reason"] == CONSTRAINT_REASON)
            .update(m_inf=1.0), "constraint m_inf")
        bad(lambda d: row(d, lambda r: r["not_tried_reason"] == DEGENERATE_REASON)
            .update(h_used=1.0), "degenerate h_used")
        bad(lambda d: row(d, lambda r: r["family"] == "para").update(m0=0.0), "para m0")
        bad(lambda d: row(d, lambda r: r["family"] == "para").update(h_used=1.0),
            "para h_used")
        bad(lambda d: row(d, quota).update(status="converged"), "post-run on pre_run")
        bad(lambda d: row(d, quota).update(tried=True), "tried with a reason")
        bad(lambda d: d["rows"][0].pop("cap_reason"), "missing row key")
        bad(lambda d: d["model"].pop("site_groups"), "missing model key")
        bad(lambda d: d["tier"].update(chosen_by="interim"), "chosen_by")
        bad(lambda d: d.update(schema_version=1), "schema version")
        bad(lambda d: d["limits"]["structural_exclusions"].pop(), "exclusions")
        bad(lambda d: d["limits"].update(fixed_sentence=SENTENCE[:-1]), "sentence")
        bad(lambda d: d["hint"].update(q_set=["0,0,0"]), "skipped hint q_set")
        bad(lambda d: d["model"]["su2"].update(sz_residual=0.0), "su2 sz_residual")
        # a seeded para row has h_used = 0 exactly
        bad(lambda d: row(d, lambda r: r["family"] == "para").update(h_used=None),
            "seeded para h_used null")
        bad(lambda d: row(d, quota).update(failure_kind="crash"), "failure_kind not failed")
        for ok in ("auto", "user", "auto (interim rule)"):
            d = copy.deepcopy(doc)
            d["tier"]["chosen_by"] = ok
            report.validate_coverage(d)

    def test_post_run_rules(self):
        from hwave.initial import report
        cset = _cset(self.dest())
        rows = _full_rows(cset)
        for r in rows:
            if r["tried"]:
                r.update(status="converged", free_energy_total=-10.0, energy_total=-9.0,
                         run_dir=r["run_id"])
        doc = report.coverage_document(_with_rows(cset, rows), phase_marker="post_run")
        report.validate_coverage(doc)
        self.assertTrue(all(r["phase_marker"] == "post_run" for r in doc["rows"]))
        d = copy.deepcopy(doc)
        next(r for r in d["rows"] if r["tried"]).update(status=None)
        with self.assertRaises(ValueError):
            report.validate_coverage(d)
        d = copy.deepcopy(doc)
        next(r for r in d["rows"] if r["tried"]).update(status="done")
        with self.assertRaises(ValueError):
            report.validate_coverage(d)
        d = copy.deepcopy(doc)
        next(r for r in d["rows"] if not r["tried"]).update(status="not_run")
        with self.assertRaises(ValueError):
            report.validate_coverage(d)
        with self.assertRaises(ValueError):
            report.coverage_document(cset, phase_marker="after")

    def test_quota_boundary_split_is_reported(self):
        from hwave.initial import report
        cset = _cset(self.dest())
        doc = report.coverage_document(cset, phase_marker="pre_run")
        qb = doc["limits"]["quota_boundary"]
        b4 = [x for x in qb if x["subshape"] == [4, 1, 1]]
        self.assertEqual(len(b4), 1, qb)
        self.assertEqual(b4[0]["order"], 4)
        self.assertIn("1/4,0,0", b4[0]["tried_q"])
        self.assertIn("3/4,0,0", b4[0]["untried_q"])
        text = report.summary_text(doc, None)
        self.assertIn("quota reached inside q-order n = 4", text)

    def test_sz_warning_surfaces(self):
        from hwave.initial import report
        from hwave.initial.symmetry import SZ_WARNING
        cset = _cset(self.dest())
        model = dict(cset.model)
        msg = SZ_WARNING + ": PairLift; the solver drops them under the 2Sz block mask"
        model["warnings"] = list(model["warnings"]) + [msg]
        doc = report.coverage_document(dataclasses.replace(cset, model=model),
                                       phase_marker="pre_run")
        self.assertIn(msg, doc["model"]["warnings"])
        text = report.summary_text(doc, None)
        self.assertIn("WARNING: " + msg, text)

    def test_tier_chosen_by_user(self):
        from hwave.initial import report
        cset = _cset(self.dest())
        doc = report.coverage_document(cset, phase_marker="pre_run", tier_chosen_by="user")
        self.assertEqual(doc["tier"]["chosen_by"], "user")
        with self.assertRaises(ValueError):
            report.coverage_document(cset, phase_marker="pre_run", tier_chosen_by="me")


class TestExclusions(_Tmp):

    def test_exclusions_and_sentence_verbatim(self):
        from hwave.initial import report
        self.assertEqual(tuple(report.STRUCTURAL_EXCLUSIONS), EXCLUSIONS)
        self.assertEqual(report.FIXED_SENTENCE, SENTENCE)
        cset = _cset(self.dest())
        doc = report.coverage_document(cset, phase_marker="pre_run")
        self.assertEqual(doc["limits"]["structural_exclusions"], list(EXCLUSIONS))
        self.assertEqual(doc["limits"]["fixed_sentence"], SENTENCE)
        ranking = report.rank(doc["rows"], facts=doc["model"])
        for text in (report.summary_text(doc, None), report.summary_text(doc, ranking)):
            for i, item in enumerate(EXCLUSIONS, 1):
                self.assertIn("({}) {}".format(i, item), text)
            self.assertIn(SENTENCE, text)
            self.assertTrue(text.rstrip().endswith(SENTENCE))
            self.assertIn("Tier: large (chosen by auto (interim rule))", text)


# --------------------------------------------------------------------------
# ranking
# --------------------------------------------------------------------------

def _facts(T=0.1, *, groups=((0,), (1,)), interacting=(0, 1), ambiguous=(),
           source="clustered", cellshape=(4, 4, 1)):
    return {"T": T, "Ncond_per_cell": "2", "2Sz": None, "cellshape": list(cellshape),
            "interacting_orbitals": list(interacting),
            "site_groups": {"source": source, "groups": [list(g) for g in groups],
                            "ambiguous": list(ambiguous)}}


def _r(run_id, status="converged", F=None, E=None, *, cellshape=(4, 4, 1), cap=None,
       tried=True, **extra):
    row = {"run_id": run_id, "candidate_id": run_id, "family": "baseline",
           "pattern_label": "lab " + run_id, "final_label": "fin " + run_id,
           "tried": tried, "not_tried_reason": None if tried else "quota: baseline",
           "status": status if tried else None, "free_energy_total": F,
           "energy_total": E, "cellshape": list(cellshape), "cap_reason": cap}
    row.update(extra)
    return row


class TestRanking(unittest.TestCase):

    def test_ranking_rules(self):
        from hwave.initial import report
        vol = 16
        rows = [
            _r("a", F=-1.0 * vol, E=-0.9 * vol),
            _r("b", F=-0.99998 * vol),                 # dF per site 1e-5: near-degenerate
            _r("c", F=-0.5 * vol, cap="target not reached at the cap"),
            _r("d", "convergence_unverified", F=-2.0 * vol),
            _r("e", "not_converged", F=-3.0 * vol),
            _r("f", "failed"),
            _r("g", "not_run"),
            _r("h", "seed_unverified", F=-4.0 * vol),
            _r("i", F=-1.0 * 8, cellshape=(4, 2, 1)),
            _r("j", tried=False),
            _r("k", F=-0.7 * vol, T=0.2),
            _r("l", F=-0.6 * vol, **{"2Sz": 0}),
        ]
        rk = report.rank(rows, facts=_facts())
        main = rk.lists[0]
        self.assertEqual(main.kind, "converged")
        self.assertEqual([e.run_id for e in main.entries], ["a", "b", "c"])
        self.assertEqual(main.group["cellshape"], [4, 4, 1])
        a, b, c = main.entries
        self.assertAlmostEqual(a.free_energy_per_cell, -1.0)
        self.assertAlmostEqual(a.free_energy_per_site, -0.5)     # two interacting sites
        self.assertEqual(a.key_source, "FreeEnergy_Total")
        self.assertAlmostEqual(b.delta_per_site, 1e-5)
        self.assertIn("near-degenerate", a.flags)
        self.assertIn("near-degenerate", b.flags)
        self.assertEqual(a.near_degenerate_with, ("b",))
        self.assertNotIn("near-degenerate", c.flags)
        self.assertIn("cap hit", c.flags)
        self.assertEqual(main.degeneracy_scale, "per site")
        # the CellShape and T splits
        groups = [(tuple(l.group["cellshape"]), l.group["T"], l.group["2Sz"])
                  for l in rk.lists]
        self.assertEqual(sorted(groups[1:], key=repr),
                         sorted([((4, 2, 1), 0.1, None), ((4, 4, 1), 0.2, None),
                                 ((4, 4, 1), 0.1, 0)], key=repr))
        self.assertEqual(len(rk.lists), 4)
        self.assertTrue(any("2Sz 0 differs from the model's free" in w
                            for w in rk.warnings), rk.warnings)
        self.assertTrue(any("CellShape [4, 2, 1]" in w and "k-mesh" in w
                            for w in rk.warnings), rk.warnings)
        self.assertTrue(any("T = 0.2" in w for w in rk.warnings), rk.warnings)
        # the flagged list
        self.assertEqual([e.run_id for e in rk.flagged_unverified.entries], ["d"])
        self.assertEqual(rk.flagged_unverified.kind, "convergence_unverified")
        self.assertIn("convergence unverified", rk.flagged_unverified.entries[0].flags)
        # unranked, never dropped
        un = {r["run_id"]: r for r in rk.unranked}
        self.assertEqual(sorted(un), ["e", "f", "g", "h"])
        self.assertIn("not converged", un["e"]["unranked_reason"])
        self.assertNotIn("j", un)
        js = rk.to_json()
        self.assertEqual(js["lists"][0]["entries"][0]["run_id"], "a")

    def test_missing_free_energy(self):
        from hwave.initial import report
        rows = [_r("a", E=-16.0), _r("b", F=-15.0)]
        rk = report.rank(rows, facts=_facts(T=0.0))
        e = {x.run_id: x for x in rk.lists[0].entries}
        self.assertEqual(e["a"].key_source, "Energy_Total (T = 0)")
        self.assertAlmostEqual(e["a"].free_energy_per_cell, -1.0)
        rk = report.rank(rows, facts=_facts(T=0.1))
        self.assertEqual([x.run_id for x in rk.lists[0].entries], ["b"])
        un = {r["run_id"]: r for r in rk.unranked}
        self.assertIn("FreeEnergy_Total", un["a"]["unranked_reason"])
        self.assertIn("T > 0", un["a"]["unranked_reason"])

    def test_ambiguous_grouping_judged_per_cell(self):
        from hwave.initial import report
        facts = _facts(groups=((0, 1),), ambiguous=("group 0 is wider than site_tol",))
        vol = 16
        rows = [_r("a", F=-1.0 * vol), _r("b", F=(-1.0 + 5e-5) * vol),
                _r("c", F=(-1.0 + 4e-4) * vol)]
        rk = report.rank(rows, facts=facts)
        lst = rk.lists[0]
        self.assertEqual(lst.degeneracy_scale, "per cell")
        a, b, c = lst.entries
        self.assertIsNone(a.free_energy_per_site)
        self.assertIn("ambiguous", a.free_energy_per_site_reason)
        self.assertIn("near-degenerate", a.flags)
        self.assertEqual(b.near_degenerate_with, ("a",))
        self.assertNotIn("near-degenerate", c.flags)
        # a user grouping is never ambiguous: per site again
        rk = report.rank(rows, facts=_facts(groups=((0, 1),), source="user"))
        self.assertEqual(rk.lists[0].degeneracy_scale, "per site")
        self.assertAlmostEqual(rk.lists[0].entries[0].free_energy_per_site, -1.0)
        # per site with two interacting groups: 4e-4 per cell is 2e-4 per site
        rk = report.rank(rows, facts=_facts())
        self.assertNotIn("near-degenerate", rk.lists[0].entries[2].flags)

    def test_energy_fields(self):
        from hwave.initial import report
        out = report.energy_fields(_r("a", F=-32.0), facts=_facts(interacting=(1,)))
        self.assertAlmostEqual(out["free_energy_per_cell"], -2.0)
        self.assertAlmostEqual(out["free_energy_per_site"], -2.0)   # one interacting site
        out = report.energy_fields(_r("a", F=-32.0), facts=_facts(interacting=()))
        self.assertIsNone(out["free_energy_per_site"])
        self.assertIn("no interacting site", out["free_energy_per_site_reason"])

    def test_summary_lists_every_status(self):
        from hwave.initial import report
        vol = 16
        rows = [_r("a", F=-1.0 * vol), _r("d", "convergence_unverified", F=-2.0 * vol),
                _r("e", "not_converged", F=-3.0 * vol),
                _r("c", F=-0.5 * vol, cap="target not reached at the cap")]
        rk = report.rank(rows, facts=_facts())
        doc = {"schema_version": 4, "phase_marker": "post_run",
               "model": dict(_facts(), norb_phys=2, n_c=2, spin_orbital=False,
                             iflag_fock=True, Ncond=32, warnings=[],
                             su2={"status": "passed"}),
               "tier": {"name": "large", "chosen_by": "auto (interim rule)",
                        "estimate": {"bytes": None, "seconds": None, "class": None}},
               "hint": {"status": "skipped", "skip_reason": "large tier"},
               "rows": rows,
               "limits": {"T_hint": {"value": 0.1, "rule": "T_hint = T", "constants": {}},
                          "label_thresholds": {"dominant_weight": 0.1,
                                               "tr_minor_weight": 0.1},
                          "quota_boundary": [],
                          "structural_exclusions": list(EXCLUSIONS),
                          "fixed_sentence": SENTENCE}}
        text = report.summary_text(doc, rk)
        for rid in ("fin a", "fin d", "fin e", "fin c"):
            self.assertIn(rid, text)
        self.assertIn("not_converged", text)
        self.assertIn("cap hit", text)
        self.assertIn("lowest among the candidates tried", text.lower())


class TestRankingFixRound(unittest.TestCase):

    def test_failed_reasons_are_specific(self):
        from hwave.initial import report
        rows = [
            _r("z", "failed", iterations=0, error_summary="no SCF step performed"),
            _r("z2", "failed", iterations=0),
            _r("x", "failed", exit_code=1, error_summary="Traceback: ValueError"),
            _r("u", "failed", iterations=12,
               error_summary="green.npz: the 'green' member is not readable"),
            _r("k", "failed", failure_kind="crash", error_summary="solver abort"),
            _r("n", "not_converged", F=-1.0, iterations=2, iteration_max=2,
               scf_residual=3.5e-3, scf_eps_tol=1e-8),
            _r("s", "seed_unverified"),
            _r("g", "not_run"),
        ]
        rk = report.rank(rows, facts=_facts())
        un = {r["run_id"]: r["unranked_reason"] for r in rk.unranked}
        self.assertEqual(un["z"], "failed (zero SCF steps): no SCF step performed")
        self.assertEqual(un["z2"], "failed (zero SCF steps): no SCF step performed")
        self.assertEqual(un["x"], "failed (solver crash or abort), exit code 1: "
                                  "Traceback: ValueError")
        self.assertEqual(un["u"], "failed (outputs missing or unreadable): green.npz: the "
                                  "'green' member is not readable")
        self.assertEqual(un["k"], "failed (solver crash or abort): solver abort")
        self.assertIn("2 iterations of IterationMax = 2", un["n"])
        self.assertIn("residual 3.500e-03", un["n"])
        self.assertIn("tolerance 1.000e-08", un["n"])
        self.assertTrue(un["s"].startswith("seed unverified"))
        self.assertTrue(un["g"].startswith("not run"))
        self.assertEqual(len(set(un[k] for k in ("z", "x", "u"))), 3)
        text = report.summary_text(_summary_doc(rows), rk)
        for k in ("z", "x", "u", "n"):
            self.assertIn(un[k], text)

    def test_group_key_is_canonical(self):
        from hwave.initial import report
        facts = dict(_facts(T=0), **{"Ncond_per_cell": "1", "2Sz": 0})
        rows = [_r("a", F=-16.0, Ncond_per_cell=1), _r("b", F=-15.0, Ncond_per_cell="1"),
                _r("c", F=-14.0, Ncond_per_cell=1.0, **{"2Sz": "0"}),
                _r("d", F=-13.0, T=0, **{"2Sz": 0.0}), _r("e", F=-12.0)]
        rk = report.rank(rows, facts=facts)
        self.assertEqual(len(rk.lists), 1)
        self.assertEqual(rk.warnings, ())
        self.assertEqual([e.run_id for e in rk.lists[0].entries], ["a", "b", "c", "d", "e"])
        self.assertEqual(rk.lists[0].group["Ncond_per_cell"], "1")
        # free spelled two ways is one group; a different value is not
        facts = dict(_facts(), **{"2Sz": None})
        rk = report.rank([_r("a", F=-16.0, **{"2Sz": "free"}), _r("b", F=-15.0),
                          _r("c", F=-14.0, Ncond_per_cell="3/2")], facts=facts)
        self.assertEqual(len(rk.lists), 2)
        self.assertEqual(len(rk.warnings), 1)
        self.assertIn("Ncond per cell = 3/2", rk.warnings[0])

    def test_tolerance_clusters_and_notation(self):
        from hwave.initial import report
        vol = 16
        rows = [_r("a", F=-1.0 * vol), _r("b", F=(-1.0 + 6e-4) * vol),
                _r("c", F=(-1.0 + 1.2e-3) * vol), _r("d", F=-0.5 * vol)]
        rk = report.rank(rows, facts=_facts(groups=((0,),), interacting=(0,)),
                         degeneracy_tol=1e-3)
        self.assertEqual(rk.degeneracy_tol, 1e-3)
        self.assertEqual(rk.to_json()["degeneracy_tol"], 1e-3)
        text = report.summary_text(_summary_doc(rows), rk)
        lines = [l for l in text.splitlines() if l.strip().startswith("near-degenerate:")]
        self.assertEqual(lines, ["  near-degenerate: a, b, c (|dF| per site < 0.001 in "
                                 "input energy units)"])
        self.assertIn("2Sz free", text)
        self.assertNotIn("2Sz None", text)
        self.assertIn("2Sz = free", text)

    def test_unverified_groups_kept_apart(self):
        from hwave.initial import report
        vol = 16
        rows = [_r("a", "convergence_unverified", F=-1.0 * vol),
                _r("b", "convergence_unverified", F=-2.0 * 8, cellshape=(4, 2, 1)),
                _r("c", "convergence_unverified", F=(-1.0 + 1e-6) * vol)]
        rk = report.rank(rows, facts=_facts())
        fl = rk.flagged_unverified
        self.assertIsNone(fl.group)
        e = {x.run_id: x for x in fl.entries}
        self.assertEqual(e["b"].delta_per_cell, 0.0)
        self.assertEqual(e["a"].delta_per_cell, 0.0)
        self.assertEqual(e["b"].near_degenerate_with, ())
        self.assertEqual(e["a"].near_degenerate_with, ("c",))
        self.assertEqual(e["b"].group["cellshape"], [4, 2, 1])
        self.assertTrue(any("within each group" in w for w in rk.warnings))

    def test_non_finite_energy_unranked(self):
        from hwave.initial import report
        rk = report.rank([_r("a", F=float("nan")), _r("b", F=float("inf")),
                          _r("c", F=-1.0)], facts=_facts())
        self.assertEqual([e.run_id for e in rk.lists[0].entries], ["c"])
        un = {r["run_id"]: r["unranked_reason"] for r in rk.unranked}
        self.assertIn("not finite", un["a"])
        self.assertIn("not finite", un["b"])


def _summary_doc(rows):
    return {"schema_version": 4, "phase_marker": "post_run",
            "model": dict(_facts(), norb_phys=2, n_c=2, spin_orbital=False,
                          iflag_fock=True, Ncond=32, warnings=[], su2={"status": "passed"}),
            "tier": {"name": "large", "chosen_by": "auto (interim rule)",
                     "estimate": {"bytes": None, "seconds": None, "class": None}},
            "hint": {"status": "skipped", "skip_reason": "large tier"},
            "rows": rows,
            "limits": {"T_hint": {"value": 0.1, "rule": "T_hint = T", "constants": {}},
                       "label_thresholds": {"dominant_weight": 0.1, "tr_minor_weight": 0.1},
                       "quota_boundary": [], "structural_exclusions": list(EXCLUSIONS),
                       "fixed_sentence": SENTENCE}}


# --------------------------------------------------------------------------
# order parameters
# --------------------------------------------------------------------------

def _run_hwave(toml):
    from hwave import qlms
    logger = logging.getLogger("qlms")
    level = logger.level
    root = logging.getLogger()
    root_level, handlers = root.level, list(root.handlers)
    old = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(toml)))
    try:
        qlms.run(input_file=os.path.basename(toml))
    finally:
        os.chdir(old)
        logger.setLevel(level)
        root.setLevel(root_level)
        for h in list(root.handlers):
            if h not in handlers:
                root.removeHandler(h)


class TestOrderParameters(_Tmp):

    def two_orbital_run(self):
        import tomli
        t = {}
        for a in range(2):
            t[((1, 0, 0), a, a)] = -1.0
            t[((-1, 0, 0), a, a)] = -1.0
        t[((0, 0, 0), 0, 1)] = -0.3
        t[((0, 0, 0), 1, 0)] = -0.3
        centres = [(0.0, 0.0, 0.0), (0.0005, 0.0, 0.0)]
        toml = initial_fixtures.write_input(
            self.dest(), cellshape=(4, 1, 1), norb=2, transfer=t, centers=centres,
            interactions={"CoulombIntra": {((0, 0, 0), 0, 0): 4.0,
                                           ((0, 0, 0), 1, 1): 4.0}},
            param={"Ncond": 8, "T": 0.1, "SubShape": [2, 1, 1], "IterationMax": 20,
                   "Mix": 0.5})
        with open(toml, "rb") as f:
            doc = tomli.load(f)
        doc["file"]["output"]["green"] = "green"
        with open(toml, "w") as f:
            f.write(dump_toml(doc))
        _run_hwave(toml)
        return toml, np.array(centres)

    def test_order_parameters_per_orbital_and_site(self):
        from hwave.initial import report
        toml, centres = self.two_orbital_run()
        groups = {"source": "clustered", "groups": [[0, 1]], "ambiguous": []}
        out = report.order_parameter_records(toml, site_groups=groups,
                                             orbital_centres=centres)
        po, ps = out["per_orbital"], out["per_site"]
        self.assertIsNotNone(po)
        self.assertIsNotNone(ps)
        self.assertIsNone(out["per_site_reason"])
        self.assertEqual(len(po["spin"]), 4)       # 2 orbitals x 2 cell positions
        self.assertEqual(len(ps["spin"]), 2)       # 1 site x 2 cell positions
        po_spin = np.array(po["spin"])
        ps_spin = np.array(ps["spin"])
        np.testing.assert_allclose(ps_spin[0], po_spin[0] + po_spin[1], atol=1e-12)
        np.testing.assert_allclose(ps_spin[1], po_spin[2] + po_spin[3], atol=1e-12)
        np.testing.assert_allclose(np.array(ps["charge"]),
                                   [po["charge"][0] + po["charge"][1],
                                    po["charge"][2] + po["charge"][3]], atol=1e-12)
        self.assertEqual(po["q"], ["0,0,0", "1/2,0,0"])
        self.assertEqual(len(po["spin_fourier"]), 2)
        self.assertEqual(len(po["spin_fourier"][0]), 3)      # full spin vector
        self.assertEqual(len(po["charge_fourier"]), 2)
        np.testing.assert_allclose(po["site_positions"][1], [0.0005, 0.0, 0.0])
        # ambiguous grouping: per orbital only, per site null with the reason
        amb = dict(groups, ambiguous=["site 0 joins orbitals 0, 1 wider than site_tol"])
        out = report.order_parameter_records(toml, site_groups=amb,
                                             orbital_centres=centres)
        self.assertIsNotNone(out["per_orbital"])
        self.assertIsNone(out["per_site"])
        self.assertIn("ambiguous", out["per_site_reason"])
        self.assertIn("wider than site_tol", out["per_site_reason"])
        # a user grouping is used even when the clustering would be ambiguous
        out = report.order_parameter_records(
            toml, site_groups=dict(amb, source="user", ambiguous=[]),
            orbital_centres=centres)
        self.assertIsNotNone(out["per_site"])


if __name__ == "__main__":
    unittest.main()
