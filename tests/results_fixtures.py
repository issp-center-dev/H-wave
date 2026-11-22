"""Shared fixtures and helpers of the ``hwave.results`` tests.

``fixture(name)`` returns the input TOML of a small H-wave run. The run
happens once per process, in-process (``hwave.qlms.run``) with the working
directory restored afterwards, into a temporary tree that is removed when the
process exits. Call it from ``setUpClass``. Tests that change files work on
``copy_fixture(name)``, a fresh copy of the tree, and never on the cached one.

Fixtures (each later task adds the ones it uses):

* ``uhfr_tutorial``: docs/tutorial/Hubbard/UHFr (Nsite 8, ``2Sz = 0``), with
  greenone.def replaced by the 16 diagonal entries ``(i, s, i, s)``.
* ``uhfk_cdw``: docs/tutorial/Hubbard/UHFk (V = 4).

Helpers:

* ``_freeze`` (re-exported from ``hwave.results._frozen``): a deep private
  copy with every array marked read-only.
* ``assert_deep_frozen(obj)``: fails unless every array reachable from
  ``obj`` is read-only (its base arrays too) and every dataclass is frozen.
"""
import atexit
import dataclasses
import logging
import os
import shutil
import tempfile
from pathlib import Path

import numpy as np

from hwave.results._frozen import _freeze  # noqa: F401  (shared helper)

REPO = Path(__file__).resolve().parent.parent

_ROOT = None
_CACHE = {}


def _session_root():
    global _ROOT
    if _ROOT is None:
        _ROOT = Path(tempfile.mkdtemp(prefix="hwave-results-"))
        atexit.register(shutil.rmtree, _ROOT, ignore_errors=True)
    return _ROOT


def _copy_inputs(source, dest):
    """Copy an example directory without its outputs and logs."""
    shutil.copytree(source, dest, ignore=shutil.ignore_patterns("output", "run.log"))


def _run(toml):
    """Run H-wave on ``toml`` from its own directory, as ``hwave input.toml``
    run in that directory would; the working directory and the level of the
    ``qlms`` logger are restored afterwards."""
    import hwave.qlms
    cwd = os.getcwd()
    qlms_logger = logging.getLogger("qlms")
    level = qlms_logger.level
    qlms_logger.setLevel(logging.WARNING)
    try:
        os.chdir(toml.parent)
        hwave.qlms.run(input_file=toml.name)
    finally:
        os.chdir(cwd)
        qlms_logger.setLevel(level)


def _build_uhfr_tutorial(dest):
    _copy_inputs(REPO / "docs" / "tutorial" / "Hubbard" / "UHFr", dest)
    nsite = 8
    rows = ["{:5d} {:5d} {:5d} {:5d}".format(i, s, i, s)
            for i in range(nsite) for s in (0, 1)]
    header = ["=" * 31, "NCisAjs {:9d}".format(len(rows)), "=" * 31,
              "======== Green functions ======", "=" * 31]
    (dest / "greenone.def").write_text("\n".join(header + rows) + "\n")
    toml = dest / "input.toml"
    _run(toml)
    return toml


def _build_uhfk_cdw(dest):
    _copy_inputs(REPO / "docs" / "tutorial" / "Hubbard" / "UHFk", dest)
    toml = dest / "input.toml"
    _run(toml)
    return toml


_BUILDERS = {
    "uhfr_tutorial": _build_uhfr_tutorial,
    "uhfk_cdw": _build_uhfk_cdw,
}


def fixture(name):
    """The input TOML of the fixture ``name``; built on first use."""
    if name not in _CACHE:
        if name not in _BUILDERS:
            raise KeyError("unknown results fixture {!r}".format(name))
        dest = _session_root() / name
        _CACHE[name] = _BUILDERS[name](dest)
    return _CACHE[name]


def copy_fixture(name):
    """A fresh copy of the fixture tree (inputs and outputs) in its own
    temporary directory; returns the TOML of the copy."""
    toml = fixture(name)
    dest = Path(tempfile.mkdtemp(dir=_session_root(), prefix=name + "-")) / "run"
    shutil.copytree(toml.parent, dest)
    return dest / toml.name


def assert_deep_frozen(obj, where="result"):
    """Raise AssertionError unless ``obj`` is deeply frozen: every ndarray
    reachable through mappings, sequences and dataclass fields is read-only,
    along with every ndarray it is a view of, and every dataclass is
    ``frozen=True``."""
    if isinstance(obj, np.ndarray):
        if obj.flags.writeable:
            raise AssertionError("{}: array is writeable".format(where))
        base = obj.base
        while isinstance(base, np.ndarray):
            if base.flags.writeable:
                raise AssertionError("{}: array is a view of a writeable array"
                                     .format(where))
            base = base.base
        if obj.dtype == object:
            for i, item in enumerate(obj.flat):
                assert_deep_frozen(item, "{}.flat[{}]".format(where, i))
        return
    if dataclasses.is_dataclass(obj) and not isinstance(obj, type):
        if not type(obj).__dataclass_params__.frozen:
            raise AssertionError("{}: dataclass {} is not frozen"
                                 .format(where, type(obj).__name__))
        for f in dataclasses.fields(obj):
            assert_deep_frozen(getattr(obj, f.name), "{}.{}".format(where, f.name))
        return
    if isinstance(obj, dict):
        for key, value in obj.items():
            assert_deep_frozen(value, "{}[{!r}]".format(where, key))
        return
    if isinstance(obj, (list, tuple)):
        for i, value in enumerate(obj):
            assert_deep_frozen(value, "{}[{}]".format(where, i))
