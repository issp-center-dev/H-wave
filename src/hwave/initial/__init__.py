"""Candidate initial states for UHF calculations (preview).

This package builds a deterministic, auditable set of candidate initial
states for multi-start mean-field calculations. The selected state is the
lowest free energy among the candidates tried. It is not a proof of the
ground state.

Available so far: the pattern description :class:`PatternSpec` (v1, with
its published JSON schema ``patternspec.schema.json``), the candidate
types :class:`Candidate`, :class:`CandidateRun` and :class:`CandidateSet`,
the enumerator :func:`candidate_states` with its :class:`Limits`, and the
writer :func:`initial_state` with its check :func:`check_initial_state`
(:class:`CheckResult`, :class:`InitialStateError`); the enumerator and the
writer are loaded on first use, since they read the model through the
solver.
"""
from hwave.initial.candidates import (
    Candidate, CandidateRun, CandidateSet, make_candidate_id)
from hwave.initial.patterns import Entry, PatternSpec, validate_patternspec

__all__ = [
    "Candidate", "CandidateRun", "CandidateSet", "make_candidate_id",
    "Entry", "PatternSpec", "validate_patternspec", "candidate_states", "Limits",
    "initial_state", "check_initial_state", "CheckResult", "InitialStateError",
]

_LAZY = {"candidate_states": "hwave.initial.enumerate", "Limits": "hwave.initial.enumerate",
         "initial_state": "hwave.initial.writer",
         "check_initial_state": "hwave.initial.writer",
         "CheckResult": "hwave.initial.writer", "InitialStateError": "hwave.initial.writer"}


def __getattr__(name):
    if name in _LAZY:
        import importlib
        return getattr(importlib.import_module(_LAZY[name]), name)
    raise AttributeError("module {!r} has no attribute {!r}".format(__name__, name))
