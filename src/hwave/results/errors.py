"""Errors raised by :mod:`hwave.results`.

Every error is a :class:`ResultsError` and carries ``path`` (the file or
directory concerned, or ``None``) and ``hint`` (what to do about it, or
``None``). ``str(error)`` gives the message followed by the hint.
"""
from pathlib import Path

__all__ = ["ResultsError", "UnsupportedLayoutError", "MissingOutputError",
           "MissingConfigError", "AmbiguousRunError", "RunRecordError"]


class ResultsError(Exception):
    """Base class of the :mod:`hwave.results` errors."""

    def __init__(self, message, *, path=None, hint=None):
        self.message = message
        self.path = None if path is None else Path(path)
        self.hint = hint
        super().__init__(message if hint is None else "{} ({})".format(message, hint))


class UnsupportedLayoutError(ResultsError):
    """A file exists but its content or layout cannot be read safely."""


class MissingOutputError(ResultsError):
    """An output is not configured in the input, or its file is absent."""


class MissingConfigError(ResultsError):
    """The input TOML file is missing, or a directory was passed instead."""


class AmbiguousRunError(ResultsError):
    """More than one run matches and none was chosen."""


class RunRecordError(ResultsError):
    """A run record cannot be read, or it belongs to another input file."""
