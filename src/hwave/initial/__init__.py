"""Candidate initial states for UHF calculations (preview).

This package builds a deterministic, auditable set of candidate initial
states for multi-start mean-field calculations. The selected state is the
lowest free energy among the candidates tried. It is not a proof of the
ground state.

Available so far: the pattern description :class:`PatternSpec` (v1, with
its published JSON schema ``patternspec.schema.json``) and the candidate
types :class:`Candidate`, :class:`CandidateRun` and :class:`CandidateSet`.
"""
from hwave.initial.candidates import (
    Candidate, CandidateRun, CandidateSet, make_candidate_id)
from hwave.initial.patterns import Entry, PatternSpec, validate_patternspec

__all__ = [
    "Candidate", "CandidateRun", "CandidateSet", "make_candidate_id",
    "Entry", "PatternSpec", "validate_patternspec",
]
