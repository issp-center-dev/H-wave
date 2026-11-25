"""The baseline enumerator, deduplication and quotas
(``hwave.initial.enumerate``; design sections 5, 5.1, 5.2, 5.3 and 8.1).

The models are small square lattices and chains (``tests/initial_fixtures``),
one-orbital unless a test needs two sites or an SU(2)-breaking term.
"""
import json
import os
import shutil
import tempfile
import time
import unittest

from tests import initial_fixtures

#: the 120-degree pattern of design section 8.1 (Skill step 2)
_TRI_120 = {"patternspec_version": 1, "q": ["2/3", "1/3", "0"], "subshape": [3, 3, 1],
            "phase": 0, "sign": 1, "modulation": "spiral", "axes": ["x", "y"],
            "entries": [{"orb_a": 0, "orb_b": 0, "r": [0, 0, 0], "s": [1.0, 0.0]}]}

#: a collinear Neel pattern written by hand (literature row)
_NEEL = {"patternspec_version": 1, "q": ["1/2", "1/2", "0"], "subshape": [2, 2, 1],
         "phase": 0, "sign": 1, "modulation": "cos",
         "entries": [{"orb_a": 0, "orb_b": 0, "r": [0, 0, 0], "cz": [1.0, 0.0]}]}

_RESERVED = ("para", "ferro", "random")


def _literature(obj, citation="A. Author, J. Phys. 1 (2000) 1"):
    from hwave.initial import Candidate, PatternSpec
    return Candidate.literature(pattern=PatternSpec.from_json(obj), citation=citation)


class _Tmp(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="initial_enumerate_")
        self._n = 0
        self._t0 = time.perf_counter()

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)
        # each fast test method stays under 5 s (tests/heavy_tests.py)
        self.assertLess(time.perf_counter() - self._t0, 5.0)

    def dest(self):
        self._n += 1
        return os.path.join(self.root, "m{}".format(self._n))

    def square(self, L=4, **kw):
        kw.setdefault("U", 4.0)
        kw.setdefault("T", 0.1)
        kw.setdefault("subshape", (2, 2, 1))
        return initial_fixtures.square(self.dest(), L, **kw)

    def two_site_chain(self, *, delta=0.5, interactions=None, ncond=4, two_sz=None,
                       cellshape=(4, 1, 1)):
        """Two orbitals at different centres of a chain cell (two sites),
        on-site energies +delta/-delta (a crystal field), U on both."""
        transfer = {((0, 0, 0), 0, 1): -1.0, ((0, 0, 0), 1, 0): -1.0,
                    ((1, 0, 0), 1, 0): -1.0, ((-1, 0, 0), 0, 1): -1.0}
        if delta:
            transfer[((0, 0, 0), 0, 0)] = delta
            transfer[((0, 0, 0), 1, 1)] = -delta
        param = {"Ncond": ncond, "T": 0.1, "SubShape": [1, 1, 1]}
        if two_sz is not None:
            param["2Sz"] = two_sz
        inter = interactions if interactions is not None else {
            "CoulombIntra": {((0, 0, 0), 0, 0): 4.0, ((0, 0, 0), 1, 1): 4.0}}
        return initial_fixtures.write_input(
            self.dest(), cellshape=cellshape, norb=2, transfer=transfer,
            centers=[(0.0, 0.0, 0.0), (0.5, 0.0, 0.0)], interactions=inter, param=param)

    @staticmethod
    def rows_of(cset, **match):
        return [r for r in cset.rows if all(r[k] == v for k, v in match.items())]


# --------------------------------------------------------------------------
# sign patterns
# --------------------------------------------------------------------------

class TestSignPatterns(unittest.TestCase):
    def test_sign_patterns_complete(self):
        from hwave.initial.enumerate import sign_patterns
        for n, count in ((3, 4), (4, 8)):
            reps = list(sign_patterns(n, cap=1000))
            self.assertEqual(len(reps), count)
            self.assertEqual(len(set(reps)), count)
            self.assertTrue(all(eta[0] == 1 for eta in reps))
            # every class {eta, -eta} once
            classes = {min(eta, tuple(-x for x in eta)) for eta in reps}
            self.assertEqual(len(classes), count)
            # the order: all +, then by the number of -1, then lexicographic
            # on the -1 positions
            minus = [tuple(i for i, x in enumerate(eta) if x == -1) for eta in reps]
            self.assertEqual(minus, sorted(minus, key=lambda m: (len(m), m)))
            pairs = list(sign_patterns(n, cap=1000, both_signs=True))
            self.assertEqual(len(pairs), 2 * count)
            for k in range(count):
                self.assertEqual(pairs[2 * k], reps[k])
                self.assertEqual(pairs[2 * k + 1], tuple(-x for x in reps[k]))
        self.assertEqual(list(sign_patterns(3, cap=1000)),
                         [(1, 1, 1), (1, -1, 1), (1, 1, -1), (1, -1, -1)])
        self.assertEqual(list(sign_patterns(1)), [(1,)])
        self.assertEqual(list(sign_patterns(1, both_signs=True)), [(1,), (-1,)])

    def test_sign_patterns_capped_lazily(self):
        """Review Focus 4: 24 sites, 16 patterns in the stated order, fast."""
        from hwave.initial.enumerate import sign_patterns
        t0 = time.perf_counter()
        got = list(sign_patterns(24))
        self.assertLess(time.perf_counter() - t0, 0.5)
        want = [tuple([1] * 24)]
        for p in range(1, 16):
            eta = [1] * 24
            eta[p] = -1
            want.append(tuple(eta))
        self.assertEqual(got, want)
        # a pair counts twice toward the cap, and the cap may split it
        both = list(sign_patterns(24, cap=5, both_signs=True))
        self.assertEqual(len(both), 5)
        self.assertEqual(both[4], want[2])
        # lazy: an iterator, not a list
        it = sign_patterns(60)
        self.assertEqual(next(it), tuple([1] * 60))

    def test_sign_patterns_arguments(self):
        from hwave.initial.enumerate import sign_patterns
        for bad in (0, -1, 1.5, True):
            with self.assertRaises(ValueError):
                sign_patterns(bad)
        for bad in (0, -2, 2.0):
            with self.assertRaises(ValueError):
                sign_patterns(3, cap=bad)


# --------------------------------------------------------------------------
# limits and arguments
# --------------------------------------------------------------------------

class TestLimits(_Tmp):
    def test_defaults_are_the_spec_values(self):
        from hwave.initial.enumerate import Limits
        lim = Limits()
        self.assertEqual((lim.top_k, lim.max_candidates, lim.quota_instability,
                          lim.max_site_patterns, lim.n_random, lim.amplitude,
                          lim.max_denominator, lim.base_seed),
                         (3, 32, 12, 16, 4, 0.5, 16, 0))
        with self.assertRaises(ValueError):
            Limits(max_candidates=-1)
        with self.assertRaises(ValueError):
            Limits(amplitude=1.0)

    def test_unknown_limit_keyword(self):
        from hwave.initial.enumerate import Limits, candidate_states
        toml = self.square()
        with self.assertRaises(TypeError) as cm:
            candidate_states(toml, max_cand=3)
        self.assertIn("max_cand", str(cm.exception))
        with self.assertRaises(TypeError):
            Limits.from_keywords(n_randoms=2)
        # every field name is accepted
        import dataclasses
        names = [f.name for f in dataclasses.fields(Limits)]
        Limits.from_keywords(**{n: getattr(Limits(), n) for n in names})

    def test_tier_and_subshape_arguments(self):
        from hwave.initial.enumerate import candidate_states
        from hwave.initial.facts import NO_SUBSHAPE_WARNING
        toml = self.square()
        for tier in ("small", "medium"):
            with self.assertRaises(NotImplementedError) as cm:
                candidate_states(toml, tier=tier)
            self.assertIn("arrives in a later release", str(cm.exception))
        with self.assertRaises(ValueError):
            candidate_states(toml, tier="huge")
        cset = candidate_states(toml, n_random=0)
        self.assertEqual(cset.tier, "large")
        self.assertEqual({tuple(r.subshape) for r in cset.runs}, {(2, 2, 1)})
        # no SubShape in the input and none given
        bare = self.square(subshape=None)
        with self.assertRaises(ValueError) as cm:
            candidate_states(bare)
        self.assertIn(NO_SUBSHAPE_WARNING, str(cm.exception))
        # duplicate subshapes are refused
        with self.assertRaises(ValueError):
            candidate_states(toml, subshapes=[(2, 2, 1), [2, 2, 1]])


# --------------------------------------------------------------------------
# reserved rows and run identity
# --------------------------------------------------------------------------

class TestReserved(_Tmp):
    def test_reserved_counts(self):
        from hwave.initial.enumerate import candidate_states
        subshapes = [(2, 2, 1), (4, 4, 1)]
        # SU(2): one ferro direction, -O merged
        cset = candidate_states(self.square(), subshapes=subshapes)
        self.assertEqual(cset.model["su2"]["status"], "passed")
        reserved = [r for r in cset.rows if r["family"] in _RESERVED]
        self.assertEqual(len(reserved), 2 * (1 + 1) + 4)
        self.assertEqual(len([r for r in cset.runs if r.candidate.family in _RESERVED]),
                         2 * (1 + 1) + 4)
        ferro = [c for c in cset.candidates if c.family == "ferro"]
        self.assertEqual(len(ferro), 1)
        self.assertEqual(ferro[0].subshape, (1, 1, 1))
        para = [c for c in cset.candidates if c.family == "para"]
        self.assertEqual(len(para), 1)
        self.assertEqual(para[0].subshape, (1, 1, 1))
        # random: round-robin over subshapes, default_rng(base_seed + i)
        rnd = [c for c in cset.candidates if c.family == "random"]
        self.assertEqual([c.seed for c in rnd], [0, 1, 2, 3])
        self.assertEqual([c.subshape for c in rnd],
                         [(2, 2, 1), (4, 4, 1), (2, 2, 1), (4, 4, 1)])
        self.assertEqual([c.candidate_id for c in rnd],
                         ["random-s0", "random-s1", "random-s2", "random-s3"])
        for c in rnd:
            runs = [r for r in cset.runs if r.candidate is c]
            self.assertEqual([r.subshape for r in runs], [c.subshape])

        # 2Sz = 2: z only, and no map merges +z and -z
        cset2 = candidate_states(self.square(two_sz=2), subshapes=subshapes)
        reserved = [r for r in cset2.rows if r["family"] in _RESERVED]
        self.assertEqual(len(reserved), 2 * (1 + 2) + 4)
        ferro = [c for c in cset2.candidates if c.family == "ferro"]
        self.assertEqual(len(ferro), 2)
        self.assertEqual(sorted(c.sign for c in ferro), [-1, 1])
        # the uniform S^z cannot respond in 2Sz mode
        for r in cset2.rows:
            if r["family"] == "ferro":
                self.assertFalse(r["tried"])
                self.assertEqual(r["not_tried_reason"], "no response (constraint)")

    def test_ferro_merge_by_tr(self):
        """Ising breaks SU(2) but keeps time reversal: x, y, z ferro, each
        with its -O merged by time reversal."""
        from hwave.initial.enumerate import candidate_states
        toml = self.two_site_chain(delta=0.0, interactions={
            "CoulombIntra": {((0, 0, 0), 0, 0): 4.0, ((0, 0, 0), 1, 1): 4.0},
            "Ising": {((0, 0, 0), 0, 1): 1.0, ((0, 0, 0), 1, 0): 1.0}})
        cset = candidate_states(toml, n_random=0)
        self.assertEqual(cset.model["su2"]["status"], "failed")
        self.assertTrue(cset.model["tr"]["invariant"])
        ferro = [c for c in cset.candidates if c.family == "ferro"]
        self.assertEqual(len(ferro), 3)
        self.assertEqual({c.sign for c in ferro}, {1})
        self.assertEqual([c.candidate_id.split("-")[3] for c in ferro], ["x", "y", "z"])
        for c in ferro:
            self.assertIn("-O merged (time_reversal)", c.sources)
            row = self.rows_of(cset, candidate_id=c.candidate_id)[0]
            self.assertIn("-O merged (time_reversal)", row["sources"])
        # without time reversal (2Sz = 2) the -O ferro is its own candidate
        toml2 = self.two_site_chain(delta=0.0, two_sz=2)
        cset2 = candidate_states(toml2, n_random=0)
        self.assertEqual(sorted(c.sign for c in cset2.candidates if c.family == "ferro"),
                         [-1, 1])

    def test_run_identity_two_subshapes(self):
        """R2: one literature pattern, two SubShapes, two run_ids."""
        from hwave.initial.enumerate import candidate_states
        lit = _literature(_NEEL)
        cset = candidate_states(self.square(), subshapes=[(2, 2, 1), (4, 4, 1)],
                                extra=[lit])
        runs = [r for r in cset.runs if r.candidate.candidate_id == lit.candidate_id]
        self.assertEqual([r.run_id for r in runs],
                         [lit.candidate_id + "-b2x2x1", lit.candidate_id + "-b4x4x1"])
        rows = self.rows_of(cset, candidate_id=lit.candidate_id)
        self.assertEqual([r["run_id"] for r in rows], [r.run_id for r in runs])
        self.assertTrue(all(r["tried"] for r in rows))
        self.assertTrue(all(r["citation"] == lit.citation for r in rows))
        # para and ferro: one candidate, one run per B
        para = [r for r in cset.runs if r.candidate.family == "para"]
        self.assertEqual([r.subshape for r in para], [(2, 2, 1), (4, 4, 1)])
        self.assertEqual(len({r.candidate.candidate_id for r in para}), 1)

    def test_v2a_literature_not_representable(self):
        """V2a: the 120-degree pattern on a 2x2 cell is listed, not dropped."""
        from hwave.initial.enumerate import candidate_states
        toml = initial_fixtures.square(self.dest(), 12, U=2.0, T=0.1, ncond=144,
                                       subshape=(2, 2, 1))
        lit = _literature(_TRI_120)
        cset = candidate_states(toml, extra=[lit], n_random=0)
        rows = self.rows_of(cset, candidate_id=lit.candidate_id)
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertFalse(row["representable"])
        self.assertFalse(row["tried"])
        self.assertEqual(row["not_tried_reason"], "unrepresentable")
        self.assertEqual(row["reason"], "q_x denominator 3 does not divide B_x = 2")
        self.assertEqual(row["proposed_subshape"], (3, 3, 1))
        self.assertIsNone(row["proposed_cellshape"])
        self.assertIn("√3×√3", row["proposal_note"])
        self.assertNotIn(row["run_id"], [r.run_id for r in cset.runs])
        self.assertIn(lit, cset.candidates)
        # no hint fields in the large tier
        hint = cset.model["hint"]
        self.assertEqual(hint["status"], "skipped")
        self.assertIsNotNone(hint["skip_reason"])
        self.assertEqual(hint["q_set"], ())
        self.assertIsNone(hint["gamma_self_adjoint"])
        self.assertIsNone(hint["method"])
        for r in cset.rows:
            self.assertIsNone(r["lambda"])
            self.assertIsNone(r["lambda_imag"])


# --------------------------------------------------------------------------
# baseline and deduplication
# --------------------------------------------------------------------------

class TestBaselineAndDedup(_Tmp):
    def _neel_rows(self, cset, B):
        tag = "-b{}x{}x{}".format(*B)
        return [r for r in cset.rows if r["family"] == "baseline"
                and r["q_string"] == "1/2,1/2,0" and r["channel"] == "spin-z"
                and r["run_id"].endswith(tag)]

    def test_neel_merge_by_translation(self):
        from hwave.initial.enumerate import candidate_states
        cset = candidate_states(self.square(), n_random=0)
        rows = self._neel_rows(cset, (2, 2, 1))
        self.assertEqual([r["sign"] for r in rows], [1, -1])
        kept, neg = rows
        self.assertTrue(kept["tried"])
        self.assertFalse(neg["tried"])
        self.assertEqual(neg["not_tried_reason"],
                         "duplicate of {} (translation)".format(kept["run_id"]))
        self.assertIn("merged {} (translation)".format(neg["run_id"]), kept["sources"])
        # the frozen candidate is not changed
        cand = next(c for c in cset.candidates if c.candidate_id == kept["candidate_id"])
        self.assertEqual(cand.sources, ("baseline",))
        # the candidate is named by its pattern (O8) and the slug is ASCII
        self.assertEqual(cand.label, "spin-z ferro-orbital q=(1/2,1/2,0) φ=0")
        self.assertEqual(cand.candidate_id, "baseline-q1_2_1_2_0-spin-z-ferro-orbital-q-1-2"
                         "-1-2-0-phi-0-r0-phi0")
        # ferro is reserved; the baseline q = 0 spin row is its duplicate
        g = [r for r in cset.rows if r["family"] == "baseline" and r["q_string"] == "0,0,0"
             and r["channel"] == "spin-z"]
        ferro_run = next(r for r in cset.runs if r.candidate.family == "ferro")
        self.assertEqual(g[0]["not_tried_reason"],
                         "duplicate of {} (translation)".format(ferro_run.run_id))
        self.assertTrue(g[1]["not_tried_reason"].startswith(
            "duplicate of {} (".format(ferro_run.run_id)))

    def test_dedup_targets_run_ids(self):
        from hwave.initial.enumerate import candidate_states
        cset = candidate_states(self.square(), subshapes=[(2, 2, 1), (4, 4, 1)])
        run_ids = {r.run_id for r in cset.runs}
        for B in ((2, 2, 1), (4, 4, 1)):
            kept, neg = self._neel_rows(cset, B)
            self.assertIn(kept["run_id"], run_ids)
            self.assertTrue(kept["run_id"].endswith("-b{}x{}x{}".format(*B)))
            self.assertEqual(neg["not_tried_reason"],
                             "duplicate of {} (translation)".format(kept["run_id"]))
        # every duplicate names a run kept by deduplication on its own
        # SubShape (a later quota may still cut it: then its own row says so)
        by_id = {r["run_id"]: r for r in cset.rows}
        for r in cset.rows:
            reason = r["not_tried_reason"] or ""
            if reason.startswith("duplicate of "):
                target = reason[len("duplicate of "):].split(" ")[0]
                self.assertIn(by_id[target]["not_tried_reason"], (None, "quota: baseline"))
                self.assertEqual(target.rsplit("-b", 1)[1], r["run_id"].rsplit("-b", 1)[1])
        # generated rows run only on the SubShape they were written for
        for run in cset.runs:
            if run.candidate.family == "baseline":
                self.assertEqual(run.subshape, run.candidate.subshape)

    def test_crystal_field_emits_both_signs(self):
        """A q = 0 charge transfer between two inequivalent sites: no map
        takes it to -O, so both signs are runs."""
        from hwave.initial.enumerate import candidate_states
        cset = candidate_states(self.two_site_chain(), n_random=0)
        self.assertEqual(cset.model["su2"]["status"], "passed")
        rows = [r for r in cset.rows if r["family"] == "baseline" and r["channel"] == "charge"
                and r["q_string"] == "0,0,0"]
        self.assertEqual([r["sign"] for r in rows], [1, -1])
        self.assertTrue(all(r["tried"] for r in rows))
        # the in-cell spin pattern (+,-) is merged with its -O (a spin map)
        spin = [r for r in cset.rows if r["family"] == "baseline" and r["channel"] == "spin-z"
                and r["q_string"] == "0,0,0" and "η=(+,-)" in r["pattern_label"]]
        self.assertEqual(len(spin), 2)
        self.assertTrue(spin[0]["tried"])
        self.assertTrue(spin[1]["not_tried_reason"].startswith(
            "duplicate of {} (".format(spin[0]["run_id"])))

    def test_directions_and_spirals(self):
        """SU(2): z only, plus one spiral per sign pattern at q not equal to
        -q; no SU(2): x, y, z and no spiral."""
        from hwave.initial.enumerate import candidate_states
        L = 4
        cset = candidate_states(initial_fixtures.square(self.dest(), L, U=4.0, T=0.1,
                                                        subshape=(4, 1, 1)), n_random=0)
        chans = [r["channel"] for r in cset.rows if r["family"] == "baseline"]
        self.assertEqual(set(chans), {"spin-z", "spin-spiral", "charge"})
        q14 = [r for r in cset.rows if r["family"] == "baseline" and r["q_string"] == "1/4,0,0"
               and r["channel"] in ("spin-z", "spin-spiral")]
        self.assertEqual([r["channel"] for r in q14], ["spin-z", "spin-z", "spin-spiral"])
        # even n = 4: phases phi0 and phi0 + pi/4, no -O row
        self.assertEqual([r["sign"] for r in q14], [1, 1, 1])
        self.assertAlmostEqual(q14[1]["phase"] - q14[0]["phase"], 3.141592653589793 / 4)
        cset2 = candidate_states(self.two_site_chain(delta=0.0, interactions={
            "CoulombIntra": {((0, 0, 0), 0, 0): 4.0, ((0, 0, 0), 1, 1): 4.0},
            "Ising": {((0, 0, 0), 0, 1): 1.0, ((0, 0, 0), 1, 0): 1.0}}), n_random=0)
        chans = {r["channel"] for r in cset2.rows if r["family"] == "baseline"}
        self.assertEqual(chans, {"spin-x", "spin-y", "spin-z", "charge"})

    def test_odd_order_rows(self):
        """q = 1/3 on a chain: two phase rows; the charge pattern has no sign
        map, so its second row is the translated -O."""
        from hwave.initial.enumerate import candidate_states
        from hwave.initial.patterns import NEG_TRANSLATED
        toml = initial_fixtures.write_input(
            self.dest(), cellshape=(6, 1, 1),
            transfer={((1, 0, 0), 0, 0): -1.0, ((-1, 0, 0), 0, 0): -1.0},
            interactions={"CoulombIntra": {((0, 0, 0), 0, 0): 4.0}},
            param={"Ncond": 6, "T": 0.1, "SubShape": [3, 1, 1]})
        cset = candidate_states(toml, n_random=0)
        charge = [r for r in cset.rows if r["family"] == "baseline"
                  and r["q_string"] == "1/3,0,0" and r["channel"] == "charge"]
        self.assertEqual(len(charge), 2)
        self.assertIsNone(charge[0]["equivalent_to"])
        self.assertEqual(charge[1]["equivalent_to"], NEG_TRANSLATED)
        self.assertAlmostEqual(charge[1]["phase"] - charge[0]["phase"], 3.141592653589793 / 3)


    def test_dedup_keeps_reserved_first(self):
        """Whatever the input order, a reserved run is the kept one."""
        from hwave.initial import enumerate as en, symmetry
        from hwave.initial.facts import model_facts
        toml = self.square()
        f = model_facts(toml)
        su2, tr = symmetry.su2_check(f.problem), symmetry.tr_check(f.problem)
        res = en.reserved_candidates(f, [(2, 2, 1)], su2=su2, tr=tr, n_random=0)
        base = en.baseline_candidates(f, (2, 2, 1), su2=su2, tr=tr)
        runs, _ = en.expand_runs(base + res, [(2, 2, 1)], cellshape=f.cellshape)
        kept, rows = en.deduplicate(runs, facts={(2, 2, 1): f}, su2=su2, tr=tr)
        ferro = next(r for r in runs if r.candidate.family == "ferro")
        self.assertIn(ferro, kept)
        g = rows[0]
        self.assertEqual(g["channel"], "spin-z")
        self.assertEqual(g["q_string"], "0,0,0")
        self.assertEqual(g["not_tried_reason"],
                         "duplicate of {} (translation)".format(ferro.run_id))

    def test_reserved_runs_never_removed(self):
        """Two reserved runs related by a map are both kept."""
        import dataclasses
        from hwave.initial import enumerate as en, symmetry
        from hwave.initial.facts import model_facts
        f = model_facts(self.square())
        su2, tr = symmetry.su2_check(f.problem), symmetry.tr_check(f.problem)
        ferro = next(c for c in en.reserved_candidates(f, [(2, 2, 1)], su2=su2, tr=tr,
                                                       n_random=0) if c.family == "ferro")
        copy = dataclasses.replace(ferro, candidate_id=ferro.candidate_id + "-copy")
        runs, _ = en.expand_runs([ferro, copy], [(2, 2, 1)], cellshape=f.cellshape)
        kept, rows = en.deduplicate(runs, facts={(2, 2, 1): f}, su2=su2, tr=tr)
        self.assertEqual(kept, runs)
        self.assertTrue(all(r["tried"] for r in rows))

    def test_invariant_prefilter_keeps_every_map(self):
        """The pre-filter of deduplicate never hides a map: the invariants
        agree on translated, time-reversed and rotated images, and
        deduplication without the pre-filter gives the same rows."""
        from unittest import mock
        import numpy as np
        from hwave.initial import enumerate as en, symmetry
        from hwave.initial.facts import model_facts
        from hwave.initial.patterns import Entry, cos_spec, realise
        toml = self.square()
        f = model_facts(toml)
        p = f.problem
        spec = cos_spec((Entry(0, 0, (0, 0, 0), cz=1.0, cx=0.4), Entry(0, 0, (1, 0, 0), c0=0.3j)),
                        ("1/2", "0", "0"), phase=0.3, subshape=(2, 2, 1))
        O = realise(spec, f, (2, 2, 1))
        o = symmetry._to_r(O, p)
        u = symmetry.spin_rotation((0.3, -0.5, 0.8), 1.1, p.norb)
        images = [symmetry._translate(o, p, (1, 0, 0)), symmetry._time_reverse_r(o, p),
                  u @ o @ u.conj().T]
        su2 = symmetry.su2_check(p)
        self.assertEqual(su2.status, "passed")
        ref = symmetry.map_invariants(O, problem=p, su2=su2)
        for img in images:
            ok = np.fft.ifftn(img, axes=(0, 1, 2), norm="forward").reshape(O.shape)
            self.assertTrue(en._same_invariants(
                ref, symmetry.map_invariants(ok, problem=p, su2=su2)))
        # without SU(2) the per-axis norms tell x from z
        ref_off = symmetry.map_invariants(O, problem=p)
        self.assertTrue(en._same_invariants(ref_off, symmetry.map_invariants(
            np.fft.ifftn(images[1], axes=(0, 1, 2), norm="forward").reshape(O.shape),
            problem=p)))
        self.assertFalse(en._same_invariants(ref_off, symmetry.map_invariants(
            np.fft.ifftn(images[2], axes=(0, 1, 2), norm="forward").reshape(O.shape),
            problem=p)))
        for toml in (self.square(), self.two_site_chain(delta=0.0, interactions={
                "CoulombIntra": {((0, 0, 0), 0, 0): 4.0, ((0, 0, 0), 1, 1): 4.0},
                "Ising": {((0, 0, 0), 0, 1): 1.0, ((0, 0, 0), 1, 0): 1.0}})):
            a = en.candidate_states(toml, n_random=2, max_candidates=1000)
            with mock.patch.object(en, "_candidates_by_invariants",
                                   lambda kept, inv, i: list(kept)):
                b = en.candidate_states(toml, n_random=2, max_candidates=1000)
            self.assertEqual(a.to_json()["rows"], b.to_json()["rows"])

    def test_translate_strided_input(self):
        """The cached translation writes a C-ordered result for any input
        layout (a strided or Fortran-ordered O(R) included)."""
        import numpy as np
        from hwave.initial import symmetry
        from hwave.initial.facts import model_facts
        p = model_facts(self.square()).problem
        rng = np.random.default_rng(1)
        o = rng.normal(size=tuple(p.shape) + (p.nd, p.nd)) + 0j
        want = symmetry._translate(np.ascontiguousarray(o), p, (1, 1, 0))
        self.assertGreater(np.linalg.norm(want), 0)
        for variant in (np.asfortranarray(o), np.swapaxes(np.swapaxes(o, 0, 1).copy(), 0, 1)):
            np.testing.assert_array_equal(symmetry._translate(variant, p, (1, 1, 0)), want)

    def test_invariant_tolerance(self):
        import numpy as np
        from hwave.initial import enumerate as en
        ref = (np.array([-0.5, 0.5]), np.array([1.0, -1.0]), np.array([0.0, 1.0]))
        near = tuple(x + 1e-9 for x in ref)
        far = (ref[0], ref[1] + 1e-3, ref[2])
        inv = {0: ref, 1: near, 2: far, 3: ref}
        self.assertEqual(en._candidates_by_invariants([0, 1, 2], inv, 3), [0, 1])
        self.assertTrue(en._same_invariants(ref, near))
        self.assertFalse(en._same_invariants(ref, far))

    def test_site_pattern_cap(self):
        """max_site_patterns bounds the sign patterns of one (q, channel),
        an unmerged -eta counted."""
        from hwave.initial.enumerate import candidate_states
        transfer = {((0, 0, 0), a, b): -1.0 for a in range(3) for b in range(3) if a != b}
        transfer.update({((1, 0, 0), 0, 2): -0.5, ((-1, 0, 0), 2, 0): -0.5})
        for a, e in enumerate((0.3, 0.0, -0.3)):
            transfer[((0, 0, 0), a, a)] = e
        toml = initial_fixtures.write_input(
            self.dest(), cellshape=(4, 1, 1), norb=3, transfer=transfer,
            centers=[(0.0, 0.0, 0.0), (0.3, 0.0, 0.0), (0.6, 0.0, 0.0)],
            interactions={"CoulombIntra": {((0, 0, 0), a, a): 4.0 for a in range(3)}},
            param={"Ncond": 6, "T": 0.1, "SubShape": [1, 1, 1]})
        full = candidate_states(toml, n_random=0, max_candidates=1000)
        capped = candidate_states(toml, n_random=0, max_candidates=1000, max_site_patterns=3)

        def count(cset, channel):
            return len([r for r in cset.rows if r["family"] == "baseline"
                        and r["channel"] == channel])
        # three sites: 4 classes; spin -O merged (listed as duplicates), charge
        # at q = 0 (+,+,+) is the uniform charge (removed), the rest twice
        self.assertEqual(count(full, "spin-z"), 8)
        self.assertEqual(count(full, "charge"), 6)
        self.assertEqual(count(capped, "spin-z"), 6)
        # the uniform charge (+,+,+) is skipped before the cap and takes no
        # slot: (+,-,+), (-,+,-), (+,+,-), the last pair cut by the cap
        self.assertEqual(count(capped, "charge"), 3)
        skipped = full.cells[0]["skipped_patterns"]
        self.assertEqual(len(skipped), 1)
        self.assertEqual(skipped[0]["reason"], "uniform density change, fixed by N")
        self.assertEqual(skipped[0]["eta"], (1, 1, 1))
        # the name follows the sign the pattern was given: (+,-,+) has its
        # largest entry negative after the trace is removed, so it is named
        # (-,+,-); (+,+,-) likewise (-,-,+)
        labels = [r["pattern_label"] for r in capped.rows if r["channel"] == "charge"]
        self.assertEqual(len(labels), 3)
        self.assertIn("η=(-,+,-)", labels[0])
        names = " ".join(r["pattern_label"] for r in full.rows if r["channel"] == "charge")
        self.assertIn("η=(-,-,+)", names)
        self.assertNotIn("η=(+,+,-)", names)

    def test_spin_orbital_model(self):
        """Spin-orbital mode: physical indices, no SU(2) check, so x, y and z,
        each -O merged by time reversal."""
        from hwave.initial.enumerate import candidate_states
        transfer = {((1, 0, 0), a, a): -1.0 for a in range(2)}
        transfer.update({((-1, 0, 0), a, a): -1.0 for a in range(2)})
        toml = initial_fixtures.write_input(
            self.dest(), cellshape=(4, 1, 1), norb=2, transfer=transfer,
            interactions={"CoulombIntra": {((0, 0, 0), 0, 0): 3.0}},
            mode={"enable_spin_orbital": True},
            param={"Ncond": 4, "T": 0.1, "SubShape": [2, 1, 1]})
        cset = candidate_states(toml, n_random=1)
        self.assertEqual(cset.model["su2"]["status"], "not applicable (spin-orbital)")
        self.assertEqual(cset.model["interacting_orbitals"], (0,))
        ferro = [c for c in cset.candidates if c.family == "ferro"]
        self.assertEqual(len(ferro), 3)
        for c in ferro:
            self.assertEqual([(e.orb_a, e.orb_b) for e in c.pattern.entries], [(0, 0)])
        self.assertTrue(any(r.candidate.family == "random" for r in cset.runs))
        neel = [r for r in cset.rows if r["family"] == "baseline" and r["q_string"] == "1/2,0,0"
                and r["channel"] == "spin-x"]
        self.assertEqual(len(neel), 2)
        self.assertTrue(neel[0]["tried"])
        self.assertTrue(neel[1]["not_tried_reason"].startswith("duplicate of"))


# --------------------------------------------------------------------------
# quotas and edge cases
# --------------------------------------------------------------------------

class TestQuotasAndEdges(_Tmp):
    def test_quota_never_cuts_reserved(self):
        from hwave.initial.enumerate import candidate_states
        lit = _literature(_NEEL)
        toml = self.square()
        full = candidate_states(toml, extra=[lit])
        cset = candidate_states(toml, extra=[lit], max_candidates=0)
        for r in cset.rows:
            if r["family"] in _RESERVED + ("literature",):
                self.assertTrue(r["tried"], r["run_id"])
            if r["family"] == "baseline":
                self.assertFalse(r["tried"])
                if not (r["not_tried_reason"] or "").startswith(("duplicate", "no response")):
                    self.assertRegex(r["not_tried_reason"],
                                     r"^quota: baseline( \(duplicate of \S+; referent not run\))?$")
        self.assertEqual([r.candidate.family for r in cset.runs if r.candidate.family == "baseline"],
                         [])
        # every candidate keeps its rows
        self.assertEqual([c.candidate_id for c in cset.candidates],
                         [c.candidate_id for c in full.candidates])
        # a small quota keeps the first generated rows, in order
        three = candidate_states(toml, extra=[lit], max_candidates=3)
        kept = [r.run_id for r in three.runs if r.candidate.family == "baseline"]
        tried_full = [r.run_id for r in full.runs if r.candidate.family == "baseline"]
        self.assertEqual(kept, tried_full[:3])
        self.assertGreater(len(tried_full), 3)

    def test_apply_quotas_instability_cut_first(self):
        from hwave.initial.enumerate import Limits, apply_quotas
        from hwave.initial import Candidate, CandidateRun, PatternSpec, Entry

        def cand(family, i):
            spec = PatternSpec(q=("1/2", "0", "0"), subshape=(2, 1, 1), phase=0.0, sign=1,
                               modulation="cos",
                               entries=(Entry(0, 0, (0, 0, 0), cz=1.0 + i),))
            return Candidate(candidate_id="{}-{}".format(family, i), family=family,
                             sources=(family,), pattern=spec)
        runs = [CandidateRun(cand("baseline", i), (2, 1, 1)) for i in range(3)]
        runs += [CandidateRun(cand("instability", i), (2, 1, 1)) for i in range(3)]
        rows = [{"run_id": r.run_id, "family": r.candidate.family, "tried": True,
                 "not_tried_reason": None, "subshape": [2, 1, 1]} for r in runs]
        kept, out = apply_quotas(runs, rows, Limits(max_candidates=4, quota_instability=2))
        self.assertEqual([r.run_id for r in kept],
                         ["baseline-0-b2x1x1", "baseline-1-b2x1x1", "baseline-2-b2x1x1",
                          "instability-0-b2x1x1"])
        self.assertEqual([r["not_tried_reason"] for r in out],
                         [None, None, None, None, "quota: instability", "quota: instability"])
        kept, out = apply_quotas(runs, rows, Limits(max_candidates=2, quota_instability=2))
        self.assertEqual([r["not_tried_reason"] for r in out],
                         [None, None, "quota: baseline"] + ["quota: instability"] * 3)

    def test_quota_spread_over_q(self):
        """Three sites, SU(2) broken (x, y, z): the baseline is round-robin
        over q by sign-pattern rank, so Neel is tried and the C4 partners
        (0,1/2) and (1/2,0) get the same rows."""
        import collections
        from hwave.initial.enumerate import candidate_states
        L = 4
        tr = {((0, 0, 0), a, b): -1.0 for a in range(3) for b in range(3) if a != b}
        for a in range(3):
            for R in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)):
                tr[(R, a, a)] = -1.0
        toml = initial_fixtures.write_input(
            self.dest(), cellshape=(L, L, 1), norb=3, transfer=tr,
            centers=[(0, 0, 0), (0.3, 0, 0), (0.6, 0, 0)],
            interactions={"CoulombIntra": {((0, 0, 0), a, a): 4.0 for a in range(3)},
                          "Ising": {((0, 0, 0), 0, 1): 0.5, ((0, 0, 0), 1, 0): 0.5}},
            param={"Ncond": 3 * L * L, "T": 0.1, "SubShape": [2, 2, 1]})
        # q = 0 alone has 18 counted rows (3 axes x 4 sign patterns, charge
        # 3 pairs): a q-major order would spend the whole quota there
        cset = candidate_states(toml, n_random=0, max_candidates=16)
        self.assertEqual(cset.model["su2"]["status"], "failed")
        base = [r for r in cset.rows if r["family"] == "baseline"]
        tried = collections.Counter(r["q_string"] for r in base if r["tried"])
        self.assertEqual(sum(tried.values()), 16)
        self.assertGreater(tried["0,1/2,0"], 0)
        self.assertEqual(tried["0,1/2,0"], tried["1/2,0,0"])
        neel = [r for r in base if r["q_string"] == "1/2,1/2,0" and r["sign"] == 1
                and r["channel"] == "spin-z" and "η=(+,+,+)" in r["pattern_label"]]
        self.assertEqual(len(neel), 1)
        self.assertTrue(neel[0]["tried"])
        for q in ("0,1/2,0", "1/2,0,0"):
            self.assertTrue(any(r["tried"] and r["channel"] == "spin-z" and r["q_string"] == q
                                for r in base))

    def test_duplicate_referent_cut_by_quota(self):
        from hwave.initial.enumerate import Limits, apply_quotas
        from hwave.initial import Candidate, CandidateRun, PatternSpec, Entry
        spec = PatternSpec(q=("1/2", "0", "0"), subshape=(2, 1, 1), phase=0.0, sign=1,
                           modulation="cos", entries=(Entry(0, 0, (0, 0, 0), cz=1.0),))
        runs = [CandidateRun(Candidate(candidate_id="baseline-{}".format(i), family="baseline",
                                       sources=("baseline",), pattern=spec), (2, 1, 1))
                for i in range(2)]
        rows = [{"run_id": r.run_id, "family": "baseline", "tried": True,
                 "not_tried_reason": None} for r in runs]
        rows.append({"run_id": "baseline-9-b2x1x1", "family": "baseline", "tried": False,
                     "not_tried_reason": "duplicate of baseline-1-b2x1x1 (translation)"})
        rows.append({"run_id": "baseline-8-b2x1x1", "family": "baseline", "tried": False,
                     "not_tried_reason": "duplicate of baseline-0-b2x1x1 (translation)"})
        kept, out = apply_quotas(runs, rows, Limits(max_candidates=1))
        self.assertEqual([r.run_id for r in kept], ["baseline-0-b2x1x1"])
        self.assertEqual(out[2]["not_tried_reason"],
                         "quota: baseline (duplicate of baseline-1-b2x1x1; referent not run)")
        self.assertEqual(out[3]["not_tried_reason"],
                         "duplicate of baseline-0-b2x1x1 (translation)")

    def test_user_rows_never_removed_by_dedup(self):
        """A user ferro equals the reserved ferro: it stays a run, marked."""
        from hwave.initial import Candidate, Entry, PatternSpec
        from hwave.initial.enumerate import candidate_states
        ferro = PatternSpec(q=("0", "0", "0"), subshape=(1, 1, 1), phase=0.0, sign=1,
                            modulation="cos", entries=(Entry(0, 0, (0, 0, 0), cz=2.0),))
        user = Candidate.user(pattern=ferro, label="my ferro")
        cset = candidate_states(self.square(), extra=[user], n_random=0)
        row = self.rows_of(cset, candidate_id=user.candidate_id)[0]
        self.assertTrue(row["tried"])
        self.assertIn(row["run_id"], [r.run_id for r in cset.runs])
        f_run = next(r for r in cset.runs if r.candidate.family == "ferro")
        self.assertEqual(row["duplicate_of"], "{} (translation)".format(f_run.run_id))
        f_row = self.rows_of(cset, run_id=f_run.run_id)[0]
        self.assertIn("equivalent {} (translation)".format(row["run_id"]), f_row["sources"])

    def test_colliding_extra_ids_refused_up_front(self):
        from unittest import mock
        from hwave.initial import enumerate as en
        a = _literature(_NEEL)
        b = _literature(_NEEL, citation="B. Author")
        self.assertEqual(a.candidate_id, b.candidate_id)
        with mock.patch.object(en, "model_facts", side_effect=AssertionError("loaded")):
            with self.assertRaises(ValueError) as cm:
                en.candidate_states(self.square(), extra=[a, b])
        self.assertIn("taken=", str(cm.exception))

    def test_nc_zero_only_para(self):
        from hwave.initial.enumerate import candidate_states, NC_ZERO_NOTE
        user = _literature(_NEEL)
        cset = candidate_states(self.square(U=0.0), extra=[user])
        self.assertEqual(cset.model["n_c"], 0)
        self.assertEqual({r.candidate.family for r in cset.runs}, {"para"})
        self.assertIn(NC_ZERO_NOTE, cset.model["warnings"])
        for r in cset.rows:
            if r["family"] != "para":
                self.assertFalse(r["tried"])
                self.assertEqual(r["not_tried_reason"], "no response")

    def test_empty_and_full_only_para(self):
        from hwave.initial.enumerate import candidate_states, EMPTY_FULL_NOTE
        L = 4
        # the solver reads Ncond = 0 as "not given", so the empty model is
        # written with Filling = 0
        empty = initial_fixtures.write_input(
            self.dest(), cellshape=(L, L, 1), transfer=initial_fixtures.square_transfer(),
            interactions={"CoulombIntra": {((0, 0, 0), 0, 0): 4.0}},
            param={"Filling": 0.0, "T": 0.1, "SubShape": [2, 2, 1]})
        for toml in (empty, self.square(L, ncond=2 * L * L)):
            cset = candidate_states(toml, subshapes=[(2, 2, 1), (4, 4, 1)])
            self.assertEqual([r.candidate.family for r in cset.runs], ["para", "para"])
            self.assertIn(EMPTY_FULL_NOTE, cset.model["warnings"])
            others = [r for r in cset.rows if r["family"] != "para"]
            self.assertTrue(others)
            for r in others:
                self.assertFalse(r["tried"])
                self.assertEqual(r["not_tried_reason"], "no response (degenerate with para)")

    def test_deterministic(self):
        from hwave.initial.enumerate import candidate_states
        from hwave.initial import CandidateSet
        lit = _literature(_NEEL)
        toml = self.square()
        a = candidate_states(toml, subshapes=[(2, 2, 1), (4, 1, 1)], extra=[lit])
        b = candidate_states(toml, subshapes=[(2, 2, 1), (4, 1, 1)], extra=[lit])
        ja = json.dumps(a.to_json(), sort_keys=True)
        self.assertEqual(ja, json.dumps(b.to_json(), sort_keys=True))
        # the JSON round trip is exact
        self.assertEqual(CandidateSet.from_json(json.loads(ja)).to_json(), a.to_json())
        # ids and run ids are unique
        ids = [c.candidate_id for c in a.candidates]
        self.assertEqual(len(ids), len(set(ids)))
        rids = [r["run_id"] for r in a.rows]
        self.assertEqual(len(rids), len(set(rids)))

    def test_t_hint_never_changes_the_candidate_set(self):
        """T_hint affects only the hint record (design section 6.4 Scope)."""
        from hwave.initial.enumerate import candidate_states
        toml = self.square()
        sets = [candidate_states(toml, subshapes=[(2, 2, 1), (4, 1, 1)], t_hint=t)
                for t in (None, 0.01, 2.0)]
        hints = [s.model["hint"]["T_hint"] for s in sets]
        self.assertEqual(len(set(hints)), 3)
        ref = sets[0].to_json()
        for s in sets[1:]:
            j = s.to_json()
            for key in ("candidates", "runs", "rows", "cells", "tier"):
                self.assertEqual(j[key], ref[key], key)


if __name__ == "__main__":
    unittest.main()
