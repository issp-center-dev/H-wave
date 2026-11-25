"""Tests of ``hwave.initial.patterns`` (PatternSpec v1) and
``hwave.initial.candidates`` (Candidate, CandidateRun, CandidateSet, ids).

A PatternSpec is the canonical description of a seed pattern: bond-slot
entries ``(orb_a, orb_b, r)`` with coefficients in the sigma (x) orbital-pair
basis, a rational q, a phase, a sign and a modulation. ``orb_a``/``orb_b``
are PHYSICAL orbital indices of the primitive cell, also in spin-orbital
mode (the spin lives in the coefficients), never the interleaved
``2*orb + spin`` solver index.
"""
import copy
import dataclasses
import importlib.util
import json
import os
import shutil
import tempfile
import unittest

import numpy as np

from hwave.initial import (
    Candidate, CandidateRun, CandidateSet, Entry, PatternSpec,
    make_candidate_id, validate_patternspec)
from hwave.initial.patterns import SCHEMA_PATH, canonical_json

#: spec section 8.1: the triangular 120-degree spiral
SPIRAL_120 = {
    "patternspec_version": 1, "q": ["2/3", "1/3", "0"],
    "subshape": [3, 3, 1], "phase": 0, "sign": 1,
    "modulation": "spiral", "axes": ["x", "y"],
    "entries": [{"orb_a": 0, "orb_b": 0, "r": [0, 0, 0], "s": [1.0, 0.0]}],
}

#: spec section 8.1: the collinear stripe entry, in a cos pattern
STRIPE = {
    "patternspec_version": 1, "q": ["1/2", "0", "0"],
    "subshape": [2, 1, 1], "phase": 0, "sign": 1, "modulation": "cos",
    "entries": [{"orb_a": 0, "orb_b": 0, "r": [0, 0, 0], "cz": [1.0, 0.0]}],
}

PARA = {
    "patternspec_version": 1, "q": ["0", "0", "0"], "subshape": [1, 1, 1],
    "phase": 0, "sign": 1, "modulation": "none", "entries": [],
}

RANDOM = {
    "patternspec_version": 1, "q": None, "subshape": [2, 2, 1],
    "phase": None, "sign": None, "modulation": "random", "entries": [],
    "seed": 3,
}

NEEL = {
    "patternspec_version": 1, "q": ["1/2", "1/2", "0"],
    "subshape": [2, 2, 1], "phase": 0, "sign": 1, "modulation": "cos",
    "entries": [{"orb_a": 0, "orb_b": 0, "r": [0, 0, 0], "cz": [1.0, 0.0]}],
}


_DROP = object()


def doc(base, **changes):
    """A copy of ``base`` with keys changed (``_DROP`` removes one)."""
    out = copy.deepcopy(base)
    for k, v in changes.items():
        if v is _DROP:
            out.pop(k, None)
        else:
            out[k] = v
    return out


def valid_documents():
    """Twelve valid PatternSpec v1 documents."""
    return [
        SPIRAL_120,
        STRIPE,
        PARA,
        RANDOM,
        # several entries, complex coefficients, an explicit +r/-r pair
        doc(STRIPE, q=["0", "0", "0"], subshape=[1, 1, 1], entries=[
            {"orb_a": 0, "orb_b": 1, "r": [1, 0, 0],
             "c0": [0.5, -0.25], "cx": [0.0, 1.0]},
            {"orb_a": 1, "orb_b": 0, "r": [-1, 0, 0], "c0": [0.5, 0.25]},
        ]),
        doc(STRIPE, sign=-1, phase=0.7853981633974483),
        doc(SPIRAL_120, axes=["y", "z"], entries=[
            {"orb_a": 0, "orb_b": 0, "r": [0, 0, 0], "s": [1.0, 0.0]},
            {"orb_a": 1, "orb_b": 1, "r": [0, 0, 0], "s": [0.0, -1.0]},
        ]),
        doc(STRIPE, entries=[
            {"orb_a": 0, "orb_b": 1, "r": [0, 0, 0], "cy": [0.0, 1.0]}]),
        doc(STRIPE, axes=_DROP, seed=_DROP),
        doc(STRIPE, axes=None, seed=None),
        doc(RANDOM, seed=0, subshape=[1, 1, 1]),
        doc(STRIPE, q=["-1/3", "2/3", "0"], subshape=[3, 3, 1]),
    ]


def invalid_documents():
    """Twelve invalid documents, each breaking one v1 rule that the
    published JSON schema can express as well."""
    entry = {"orb_a": 0, "orb_b": 0, "r": [0, 0, 0]}
    return [
        doc(STRIPE, patternspec_version=2),
        doc(STRIPE, entries=[dict(entry, cz=1.0)]),            # not [re, im]
        doc(STRIPE, entries=[dict(entry, s=[1.0, 0.0])]),      # s with cos
        doc(STRIPE, axes=["x", "y"]),                          # axes with cos
        doc(SPIRAL_120, entries=[dict(entry, s=[1.0, 0.0],
                                      cz=[1.0, 0.0])]),        # c* in spiral
        doc(SPIRAL_120, axes=["x", "x"]),                      # not distinct
        doc(SPIRAL_120, axes=_DROP),                           # spiral, axes?
        doc(PARA, entries=[dict(entry, cz=[1.0, 0.0])]),       # none, entries
        doc(PARA, phase=0.5),                                  # none, phase
        doc(RANDOM, q=["0", "0", "0"]),                        # random, q
        doc(STRIPE, seed=5),                                   # seed with cos
        doc(STRIPE, q=["1/0", "0", "0"]),                      # not a Fraction
    ]


class TestPatternSpec(unittest.TestCase):

    def test_round_trip_examples(self):
        for example in (SPIRAL_120, STRIPE, PARA, RANDOM):
            with self.subTest(modulation=example["modulation"]):
                spec = PatternSpec.from_json(example)
                again = PatternSpec.from_json(spec.to_json())
                self.assertEqual(again, spec)
                self.assertEqual(hash(again), hash(spec))
                # the canonical JSON keeps the user's document (only the
                # always-written axes/seed keys are added as null)
                expected = dict({"axes": None, "seed": None}, **example)
                self.assertEqual(spec.to_json(), expected)
                # byte-stable, both as canonical and as plain json.dumps
                self.assertEqual(canonical_json(spec.to_json()),
                                 canonical_json(again.to_json()))
                self.assertEqual(json.dumps(spec.to_json()),
                                 json.dumps(again.to_json()))
                self.assertEqual(
                    PatternSpec.from_json(json.loads(json.dumps(
                        spec.to_json()))), spec)

    def test_spiral_fields(self):
        spec = PatternSpec.from_json(SPIRAL_120)
        self.assertEqual(spec.q, ("2/3", "1/3", "0"))
        self.assertEqual(spec.subshape, (3, 3, 1))
        self.assertEqual(spec.axes, ("x", "y"))
        self.assertEqual(spec.entries,
                         (Entry(0, 0, (0, 0, 0), s=1 + 0j),))
        self.assertIsInstance(spec.phase, float)

    def test_none_and_random_rules(self):
        violations = {
            "none with entries": doc(PARA, entries=[
                {"orb_a": 0, "orb_b": 0, "r": [0, 0, 0], "c0": [1.0, 0.0]}]),
            "none with q != 0": doc(PARA, q=["1/2", "0", "0"]),
            "none with phase != 0": doc(PARA, phase=0.1),
            "none with sign != 1": doc(PARA, sign=-1),
            "none with a seed": doc(PARA, seed=1),
            "none with axes": doc(PARA, axes=["x", "y"]),
            "random with entries": doc(RANDOM, entries=[
                {"orb_a": 0, "orb_b": 0, "r": [0, 0, 0], "c0": [1.0, 0.0]}]),
            "random without a seed": doc(RANDOM, seed=None),
            "random with a bool seed": doc(RANDOM, seed=True),
            "random with a negative seed": doc(RANDOM, seed=-1),
            "random with q": doc(RANDOM, q=["0", "0", "0"]),
            "random with phase": doc(RANDOM, phase=0.0),
            "random with sign": doc(RANDOM, sign=1),
            "cos with a seed": doc(STRIPE, seed=2),
            "cos without q": doc(STRIPE, q=None),
            "cos with sign 0": doc(STRIPE, sign=0),
            "cos without entries": doc(STRIPE, entries=[]),
            "spiral with a seed": doc(SPIRAL_120, seed=2),
        }
        for name, bad in violations.items():
            with self.subTest(name):
                with self.assertRaises(ValueError) as cm:
                    validate_patternspec(bad)
                self.assertIn("$", str(cm.exception))
                with self.assertRaises(ValueError):
                    PatternSpec.from_json(bad)

    def test_error_names_json_path(self):
        bad = doc(STRIPE, entries=[
            {"orb_a": 0, "orb_b": 0, "r": [0, 0, 0], "cz": [1.0, 0.0]},
            {"orb_a": 0, "orb_b": 1, "r": [0, 0], "cz": [1.0, 0.0]}])
        with self.assertRaisesRegex(ValueError, r"\$\.entries\[1\]\.r"):
            validate_patternspec(bad)
        bad = doc(STRIPE, entries=[
            {"orb_a": 0, "orb_b": 0, "r": [0, 0, 0], "cz": [1.0, "x"]}])
        with self.assertRaisesRegex(ValueError, r"\$\.entries\[0\]\.cz"):
            validate_patternspec(bad)
        with self.assertRaisesRegex(ValueError, r"\$\.q\[0\]"):
            validate_patternspec(doc(STRIPE, q=["0.5", "0", "0"]))
        with self.assertRaisesRegex(ValueError, r"\$\.colour"):
            validate_patternspec(doc(STRIPE, colour="red"))
        with self.assertRaisesRegex(ValueError, r"\$\.entries\[0\]\.cz"):
            validate_patternspec(doc(STRIPE, entries=[
                {"orb_a": 0, "orb_b": 0, "r": [0, 0, 0],
                 "cz": [float("nan"), 0.0]}]))
        with self.assertRaisesRegex(ValueError, r"\$"):
            validate_patternspec([STRIPE])

    def test_duplicate_entry_refused(self):
        entries = [
            {"orb_a": 0, "orb_b": 1, "r": [1, 0, 0], "cz": [1.0, 0.0]},
            {"orb_a": 0, "orb_b": 1, "r": [1, 0, 0], "c0": [1.0, 0.0]},
        ]
        with self.assertRaisesRegex(ValueError, r"\$\.entries\[1\].*duplic"):
            validate_patternspec(doc(STRIPE, entries=entries))
        with self.assertRaises(ValueError):
            PatternSpec(q=("1/2", "0", "0"), subshape=(2, 1, 1), phase=0.0,
                        sign=1, modulation="cos",
                        entries=(Entry(0, 1, (1, 0, 0), cz=1),
                                 Entry(0, 1, (1, 0, 0), c0=1)))

    def test_explicit_minus_r_partner_kept(self):
        entries = [
            {"orb_a": 0, "orb_b": 1, "r": [1, 0, 0], "c0": [0.5, 0.0]},
            {"orb_a": 1, "orb_b": 0, "r": [-1, 0, 0], "c0": [0.5, 0.0]},
        ]
        spec = PatternSpec.from_json(doc(STRIPE, entries=entries))
        # two entries of A; nothing doubled, nothing merged, nothing added
        self.assertEqual(len(spec.entries), 2)
        self.assertEqual(
            sorted(e.r for e in spec.entries), [(-1, 0, 0), (1, 0, 0)])
        self.assertEqual(
            [e.c0 for e in spec.entries], [0.5 + 0j, 0.5 + 0j])
        out = spec.to_json()["entries"]
        self.assertEqual(len(out), 2)
        for e in out:
            self.assertEqual(e["c0"], [0.5, 0.0])

    def test_entries_canonical_order(self):
        """Entry order does not matter: the PatternSpec sorts by
        (orb_a, orb_b, r), so equal patterns are equal objects and give
        equal JSON."""
        entries = [
            {"orb_a": 1, "orb_b": 0, "r": [0, 0, 0], "cz": [1.0, 0.0]},
            {"orb_a": 0, "orb_b": 0, "r": [0, 1, 0], "cz": [1.0, 0.0]},
            {"orb_a": 0, "orb_b": 0, "r": [0, -1, 0], "cz": [1.0, 0.0]},
        ]
        a = PatternSpec.from_json(doc(STRIPE, entries=entries))
        b = PatternSpec.from_json(doc(STRIPE, entries=entries[::-1]))
        self.assertEqual(a, b)
        self.assertEqual(canonical_json(a.to_json()),
                         canonical_json(b.to_json()))
        self.assertEqual([(e.orb_a, e.orb_b, e.r) for e in a.entries],
                         [(0, 0, (0, -1, 0)), (0, 0, (0, 1, 0)),
                          (1, 0, (0, 0, 0))])

    def test_constructor_validates(self):
        spec = PatternSpec(q=["1/2", "0", "0"], subshape=[2, 1, 1],
                           phase=0, sign=1, modulation="cos",
                           entries=[Entry(0, 0, [0, 0, 0], cz=1.0)])
        self.assertEqual(spec, PatternSpec.from_json(STRIPE))
        self.assertEqual(spec.patternspec_version, 1)
        with self.assertRaises(ValueError):          # s in a cos pattern
            PatternSpec(q=("1/2", "0", "0"), subshape=(2, 1, 1), phase=0.0,
                        sign=1, modulation="cos",
                        entries=(Entry(0, 0, (0, 0, 0), s=1),))
        with self.assertRaises(ValueError):          # c* in a spiral
            PatternSpec(q=("2/3", "1/3", "0"), subshape=(3, 3, 1),
                        phase=0.0, sign=1, modulation="spiral",
                        axes=("x", "y"),
                        entries=(Entry(0, 0, (0, 0, 0), s=1, cz=1),))
        with self.assertRaises(ValueError):
            PatternSpec(q=("0", "0", "0"), subshape=(1, 1, 1), phase=0.0,
                        sign=1, modulation="cosine", entries=())
        with self.assertRaises(ValueError):          # bool is not an int
            Entry(True, 0, (0, 0, 0), cz=1)
        with self.assertRaises(ValueError):          # str is not a number
            Entry(0, 0, (0, 0, 0), cz="1")
        with self.assertRaises(dataclasses.FrozenInstanceError):
            spec.phase = 1.0

    def test_numpy_values_normalised(self):
        """numpy scalars become plain Python values, so the JSON is the
        same as for Python inputs (numpy 1.26 and 2 alike)."""
        spec = PatternSpec(
            q=("1/2", "0", "0"), subshape=tuple(np.array([2, 1, 1])),
            phase=np.float64(0.0), sign=np.int64(1), modulation="cos",
            entries=(Entry(np.int64(0), np.int32(0), tuple(np.zeros(3, int)),
                           cz=np.complex128(1.0)),))
        self.assertEqual(canonical_json(spec.to_json()),
                         canonical_json(PatternSpec.from_json(
                             STRIPE).to_json()))
        for v in spec.subshape + spec.entries[0].r:
            self.assertIs(type(v), int)

    def test_so_physical_index_bound(self):
        """In spin-orbital mode the geometry carries 2*norb_phys orbitals,
        but PatternSpec indices are PHYSICAL: with ``norb_phys`` given, an
        interleaved index (here 3 = 2*1 + 1) is refused."""
        good = doc(STRIPE, entries=[
            {"orb_a": 1, "orb_b": 1, "r": [0, 0, 0], "cz": [1.0, 0.0]}])
        validate_patternspec(good, norb_phys=2)
        bad = doc(STRIPE, entries=[
            {"orb_a": 3, "orb_b": 1, "r": [0, 0, 0], "cz": [1.0, 0.0]}])
        validate_patternspec(bad)       # no model: only orb >= 0 is checked
        with self.assertRaisesRegex(ValueError,
                                    r"\$\.entries\[0\]\.orb_a.*physical"):
            validate_patternspec(bad, norb_phys=2)
        with self.assertRaisesRegex(ValueError, r"orb_b"):
            validate_patternspec(doc(STRIPE, entries=[
                {"orb_a": 0, "orb_b": -1, "r": [0, 0, 0],
                 "cz": [1.0, 0.0]}]))

    def test_provenance_member_pickle_free(self):
        """The canonical JSON of a pattern stored as a 0-d unicode npz
        member loads with allow_pickle=False and decodes to the same
        pattern (the writer's provenance encoding)."""
        spec = PatternSpec.from_json(SPIRAL_120)
        tmp = tempfile.mkdtemp()
        try:
            path = os.path.join(tmp, "p.npz")
            np.savez(path, pattern=np.array(canonical_json(spec.to_json())))
            with np.load(path, allow_pickle=False) as z:
                member = z["pattern"]
            self.assertEqual(member.shape, ())
            self.assertEqual(member.dtype.kind, "U")
            self.assertEqual(PatternSpec.from_json(json.loads(str(member))),
                             spec)
        finally:
            shutil.rmtree(tmp)

    def test_schema_file_agrees_with_validator(self):
        with open(SCHEMA_PATH, encoding="utf-8") as f:
            schema = json.load(f)
        self.assertEqual(schema["properties"]["patternspec_version"],
                         {"const": 1})
        valid, invalid = valid_documents(), invalid_documents()
        self.assertEqual((len(valid), len(invalid)), (12, 12))
        for i, d in enumerate(valid):
            with self.subTest(valid=i):
                validate_patternspec(d)
        for i, d in enumerate(invalid):
            with self.subTest(invalid=i):
                with self.assertRaises(ValueError):
                    validate_patternspec(d)
        if importlib.util.find_spec("jsonschema") is not None:
            import jsonschema
            cls = jsonschema.validators.validator_for(schema)
            cls.check_schema(schema)
            checker = cls(schema)
            for i, d in enumerate(valid):
                with self.subTest(schema_valid=i):
                    self.assertEqual(list(checker.iter_errors(d)), [])
            for i, d in enumerate(invalid):
                with self.subTest(schema_invalid=i):
                    self.assertFalse(checker.is_valid(d))

    def test_schema_shipped_in_package(self):
        self.assertEqual(os.path.basename(SCHEMA_PATH),
                         "patternspec.schema.json")
        self.assertEqual(
            os.path.dirname(SCHEMA_PATH),
            os.path.dirname(importlib.util.find_spec(
                "hwave.initial").origin))


def candidate_of(family, pattern, label, **kw):
    taken = kw.pop("taken", set())
    seed = pattern.seed
    cid = make_candidate_id(family, pattern.q, label, 0, 0, False, taken,
                            seed=seed)
    return Candidate(candidate_id=cid, family=family, sources=(family,),
                     pattern=pattern, amplitude_requested=0.5, citation=None,
                     label=label, **kw)


class TestCandidates(unittest.TestCase):

    def test_candidate_derived_properties(self):
        spec = PatternSpec.from_json(NEEL)
        c = Candidate.user(pattern=spec, label="Neel")
        self.assertEqual(c.q, ("1/2", "1/2", "0"))
        self.assertEqual(c.subshape, (2, 2, 1))
        self.assertEqual(c.phase, 0.0)
        self.assertEqual(c.sign, 1)
        self.assertIsNone(c.seed)
        self.assertEqual(c.family, "user")
        self.assertEqual(c.sources, ("user",))
        self.assertEqual(c.amplitude_requested, 0.5)
        self.assertIsNone(c.citation)
        names = {f.name for f in dataclasses.fields(Candidate)}
        self.assertEqual(names, {"candidate_id", "family", "sources",
                                 "pattern", "amplitude_requested",
                                 "citation", "label"})
        for attr in ("q", "subshape", "phase", "sign", "seed"):
            with self.subTest(attr):
                with self.assertRaises(AttributeError):
                    setattr(c, attr, None)
        r = candidate_of("random", PatternSpec.from_json(RANDOM), None)
        self.assertEqual(r.seed, 3)
        self.assertIsNone(r.q)
        self.assertIsNone(r.phase)
        self.assertEqual(r.candidate_id, "random-s3")

    def test_literature_requires_citation(self):
        spec = PatternSpec.from_json(NEEL)
        for citation in (None, "", "   "):
            with self.subTest(citation=citation):
                with self.assertRaises(ValueError):
                    Candidate.literature(pattern=spec, citation=citation)
        with self.assertRaises(TypeError):
            Candidate.literature(pattern=spec)
        c = Candidate.literature(pattern=spec, citation="PRB 1, 1 (2000)",
                                 label="Neel")
        self.assertEqual(c.citation, "PRB 1, 1 (2000)")
        self.assertEqual(c.sources, ("literature",))
        with self.assertRaises(ValueError):      # citation only literature
            dataclasses.replace(Candidate.user(pattern=spec),
                                citation="PRB 1, 1 (2000)")
        with self.assertRaises(ValueError):
            dataclasses.replace(c, citation=None)

    def test_candidate_validation(self):
        spec = PatternSpec.from_json(NEEL)
        c = Candidate.user(pattern=spec)
        for bad in ({"family": "other"}, {"amplitude_requested": 1.0},
                    {"amplitude_requested": 0.0},
                    {"amplitude_requested": float("nan")},
                    {"candidate_id": ""}, {"sources": "user"},
                    {"pattern": NEEL}):
            with self.subTest(bad):
                with self.assertRaises(ValueError):
                    dataclasses.replace(c, **bad)
        # para <-> modulation none, random <-> modulation random
        with self.assertRaises(ValueError):
            dataclasses.replace(c, family="para")
        with self.assertRaises(ValueError):
            dataclasses.replace(c, family="random")
        with self.assertRaises(ValueError):
            Candidate.user(pattern=PatternSpec.from_json(RANDOM))
        with self.assertRaises(ValueError):
            Candidate.user(pattern=PatternSpec.from_json(PARA))
        # list sources are normalised to a tuple
        self.assertEqual(
            dataclasses.replace(c, sources=["user", "x"]).sources,
            ("user", "x"))

    def test_run_id_format(self):
        spec = PatternSpec.from_json(NEEL)
        c = Candidate.literature(pattern=spec, citation="ref", label="Néel")
        self.assertEqual(c.candidate_id, "literature-q1_2_1_2_0-neel-r0-phi0")
        run = CandidateRun(c, (2, 2, 1))
        self.assertEqual(run.run_id,
                         "literature-q1_2_1_2_0-neel-r0-phi0-b2x2x1")
        self.assertNotIn("run_id",
                         {f.name for f in dataclasses.fields(CandidateRun)})
        with self.assertRaises(AttributeError):
            run.run_id = "x"
        obj = run.to_json()
        self.assertEqual(obj["run_id"], run.run_id)
        self.assertEqual(CandidateRun.from_json(obj), run)
        bad = dict(obj, run_id="literature-q1_2_1_2_0-neel-r0-phi0-b4x4x1")
        with self.assertRaisesRegex(ValueError, "run_id"):
            CandidateRun.from_json(bad)
        self.assertEqual(CandidateRun(c, [4, 4, 1]).run_id[-7:], "-b4x4x1")
        with self.assertRaises(ValueError):
            CandidateRun(c, (0, 1, 1))

    def test_make_candidate_id(self):
        taken = set()
        q = ("1/2", "1/2", "0")
        a = make_candidate_id("instability", q, "spin-z", 1, 0, False, taken)
        self.assertEqual(a, "instability-q1_2_1_2_0-spin-z-r1-phi0")
        b = make_candidate_id("instability", q, "spin-z", 1, 0, False, taken)
        c = make_candidate_id("instability", q, "spin-z", 1, 0, False, taken)
        self.assertEqual((b, c), (a + "-2", a + "-3"))
        self.assertEqual(taken, {a, b, c})
        self.assertEqual(
            make_candidate_id("baseline", ("2/3", "1/3", "0"),
                              "Spin z, Néel φ=0", 0, 1, True, set()),
            "baseline-q2_3_1_3_0-spin-z-neel-0-r0-phi1-neg")
        self.assertEqual(
            make_candidate_id("baseline", ("-1/3", "4/6", "0"), None, 2, 0,
                              False, set()),
            "baseline-qm1_3_2_3_0-pattern-r2-phi0")
        self.assertEqual(
            make_candidate_id("random", None, None, 0, 0, False, set(),
                              seed=7), "random-s7")
        t2 = {"random-s7"}
        self.assertEqual(
            make_candidate_id("random", None, None, 0, 0, False, t2, seed=7),
            "random-s7-2")
        with self.assertRaises(ValueError):
            make_candidate_id("random", None, None, 0, 0, False, set())
        with self.assertRaises(ValueError):
            make_candidate_id("bogus", q, "x", 0, 0, False, set())
        with self.assertRaises(ValueError):
            make_candidate_id("baseline", q, "x", -1, 0, False, set())

    def test_candidate_set_json_round_trip(self):
        neel = PatternSpec.from_json(NEEL)
        ferro = PatternSpec.from_json(doc(
            STRIPE, q=["0", "0", "0"], subshape=[1, 1, 1]))
        orbital = PatternSpec.from_json(doc(STRIPE, entries=[
            {"orb_a": 0, "orb_b": 0, "r": [0, 0, 0], "c0": [1.0, 0.0]},
            {"orb_a": 1, "orb_b": 1, "r": [0, 0, 0], "c0": [-1.0, 0.0]}]))
        taken = set()
        cands = [
            candidate_of("para", PatternSpec.from_json(PARA), "para",
                         taken=taken),
            candidate_of("ferro", ferro, "spin-z ferro", taken=taken),
            candidate_of("random", PatternSpec.from_json(RANDOM), None,
                         taken=taken),
            candidate_of("instability", neel, "spin-z", taken=taken),
            candidate_of("baseline", neel, "spin-z", taken=taken),
            candidate_of("orbital", orbital, "orbital-polarization",
                         taken=taken),
            Candidate.user(pattern=neel, label="spin-z", taken=taken),
            Candidate.literature(pattern=PatternSpec.from_json(SPIRAL_120),
                                 citation="ref", label="120", taken=taken),
        ]
        self.assertEqual(
            {c.family for c in cands},
            {"para", "ferro", "random", "instability", "baseline",
             "orbital", "user", "literature"})
        self.assertEqual(len({c.candidate_id for c in cands}), len(cands))
        runs = [CandidateRun(c, b) for c in cands
                for b in ((2, 2, 1), (4, 4, 1))]
        cset = CandidateSet(
            model={"norb": 2, "spin_orbital": False, "input_sha256":
                   {"input.toml": "0" * 64}},
            tier="large",
            cells=[{"subshape": [2, 2, 1], "q": [["0", "0", "0"]]}],
            candidates=cands, runs=runs,
            rows=[{"run_id": r.run_id, "tried": True, "lambda": None}
                  for r in runs],
            limits={"top_k": 3, "amplitude": 0.5})
        obj = cset.to_json()
        text = canonical_json(obj)
        back = CandidateSet.from_json(json.loads(text))
        self.assertEqual(back, cset)
        self.assertEqual(canonical_json(back.to_json()), text)
        self.assertEqual(json.dumps(back.to_json()), json.dumps(obj))
        # runs refer to their candidate by id, the run_id is checked
        self.assertEqual(set(obj["runs"][0]),
                         {"candidate_id", "subshape", "run_id"})
        bad = json.loads(text)
        bad["runs"][0]["run_id"] = "x"
        with self.assertRaisesRegex(ValueError, "run_id"):
            CandidateSet.from_json(bad)
        bad = json.loads(text)
        bad["runs"][0]["candidate_id"] = "nobody"
        with self.assertRaises(ValueError):
            CandidateSet.from_json(bad)

    def test_candidate_set_consistency(self):
        spec = PatternSpec.from_json(NEEL)
        c = Candidate.user(pattern=spec)
        other = Candidate.user(pattern=PatternSpec.from_json(STRIPE),
                               label="stripe")
        kw = dict(model={}, tier="large", cells=[], rows=[], limits={})
        with self.assertRaisesRegex(ValueError, "candidate_id"):
            CandidateSet(candidates=[c, c], runs=[], **kw)
        with self.assertRaisesRegex(ValueError, "run"):
            CandidateSet(candidates=[c], runs=[CandidateRun(other, (2, 2, 1))],
                         **kw)
        with self.assertRaisesRegex(ValueError, "run_id"):
            CandidateSet(candidates=[c], runs=[CandidateRun(c, (2, 2, 1))] * 2,
                         **kw)
        with self.assertRaisesRegex(ValueError, "tier"):
            CandidateSet(candidates=[c], runs=[], **dict(kw, tier="auto"))
        with self.assertRaises(ValueError):       # not JSON
            CandidateSet(candidates=[c], runs=[], **dict(kw, model={"x": {1}}))
        with self.assertRaises(ValueError):       # NaN is not JSON
            CandidateSet(candidates=[c], runs=[],
                         **dict(kw, limits={"x": float("nan")}))
        # tuples in JSON fields are stored in their JSON (list) form, so
        # the round trip is exact
        s = CandidateSet(candidates=[c], runs=[], **dict(kw, cells=[
            {"subshape": (2, 2, 1)}]))
        self.assertEqual(s.to_json()["cells"], [{"subshape": [2, 2, 1]}])
        with self.assertRaises(TypeError):        # read-only copies
            s.model["x"] = 1
        self.assertIsInstance(s.cells, tuple)
        self.assertEqual(CandidateSet.from_json(s.to_json()), s)

    def test_candidate_json_round_trip(self):
        c = Candidate.literature(pattern=PatternSpec.from_json(SPIRAL_120),
                                 citation="ref", label="120 degree")
        obj = c.to_json()
        self.assertNotIn("q", obj)
        self.assertEqual(Candidate.from_json(obj), c)
        with self.assertRaises(ValueError):
            Candidate.from_json(dict(obj, q=["0", "0", "0"]))


if __name__ == "__main__":
    unittest.main()
