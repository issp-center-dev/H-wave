"""Shared fixtures and helpers of the ``hwave.results`` tests.

``fixture(name)`` returns the input TOML of a small H-wave run. The run
happens once per process, in-process (``hwave.qlms.run``) with the working
directory restored afterwards, into a temporary tree that is removed when the
process exits. Call it from ``setUpClass``. Tests that change files work on
``copy_fixture(name)``, a fresh copy of the tree, and never on the cached one.

Fixtures (each later task adds the ones it uses):

* ``uhfr_tutorial``: docs/tutorial/Hubbard/UHFr (Nsite 8, ``2Sz = 0``), with
  greenone.def replaced by the 16 diagonal entries ``(i, s, i, s)``.
* ``uhfr_u8``: tests/uhfr/coulombintra (Nsite 16 on a 4 x 4 square lattice,
  U = 8, ``2Sz`` unset) from its own start (``output_ref/green_init.dat``),
  with greenone.def extended to every ``(i, s, i, t)``, ``s, t`` in
  ``{0, 1}`` (64 entries); a collinear Neel state with a tilted axis,
  converged after 1 step.
* ``uhfr_tutorial_zero_steps``: ``uhfr_tutorial`` with ``IterationMax = 0``:
  only ``scf_status.dat`` is written.
* ``uhfk_cdw``: docs/tutorial/Hubbard/UHFk (V = 4).
* ``uhfk_cdw_legacy``: a copy of docs/tutorial/Hubbard/UHFk with its
  git-tracked ``output/`` (``green.npz`` and ``eigen.npz`` written before the
  ``momentum_convention`` and ``green_convention`` markers); **not run**.
* ``uhfk_u8``: tests/uhfk/CoulombIntra (U = 8, EPS 12) from its own start
  (``output_ref/green.dat.npz``), with ``[log] print_check = "check.dat"``;
  converged after 1 step.
* ``uhfk_u8_not_converged``: the same with the start removed,
  ``[file.input] initial_mode = "random"`` and ``IterationMax = 5``, as
  tests/test_uhf_scf_status.py::TestConvergence::test_uhfk_not_converged;
  stops after 5 steps without converging.
* ``uhfk_u8_zero_steps``: ``uhfk_u8`` with ``IterationMax = 0``: only
  ``scf_status.dat`` and an empty ``check.dat`` are written.
* ``uhfk_u8_negative_max``: ``uhfk_u8`` with ``IterationMax = -3`` and
  ``trustme_interaction_range = true`` (which lets it through the range
  check): the writer's ``IterationMax -3 treated as 0`` status file and an
  empty ``check.dat``.
* ``uhfk_u8_nan_residual``: ``uhfk_u8_not_converged`` with the solver's
  residual replaced by NaN from the 3rd step on (``UHFk._calc_phys``
  patched during the run, as tests/test_uhf_scf_status.py does): the
  writer's ``residual not finite after step 5`` status file.
* ``uhfk_fm_axis``: tests/uhfk/CoulombIntra with ``Ncond = 5``,
  ``SubShape = [1, 1, 1]`` and a start ``start.npz`` holding a fully
  polarised ferromagnet along the generic axis
  ``u = (0.6, 0.2537, 0.7734)/|.|`` on every cell (charge ``5/16``), as
  tests/test_uhf_transverse_spin.py::TestUHFkRuns::test_uhfk_ferromagnet_generic_axis.
  The Hubbard U keeps the global spin axis of the start, so the converged
  totals are ``(Sx, Sy, Sz) = 2.5 u``.

Helpers:

* ``_freeze`` (re-exported from ``hwave.results._frozen``): a deep private
  copy with every array marked read-only.
* ``write_greenone(path, entries)``: write a UHFr OneBodyG file;
  ``run_input(toml)``: rerun a (copied, edited) fixture input.
* ``assert_deep_frozen(obj)``: fails unless every array reachable from
  ``obj`` is read-only (its base arrays too) and every dataclass is frozen.
"""
import atexit
import dataclasses
import logging
import os
import shutil
import tempfile
from collections.abc import Mapping
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
    ``qlms`` logger are restored afterwards, and so is the root logger
    (``qlms.run`` calls ``logging.basicConfig(level=INFO)``): its level and
    its handlers."""
    import hwave.qlms
    cwd = os.getcwd()
    qlms_logger = logging.getLogger("qlms")
    root_logger = logging.getLogger()
    level = qlms_logger.level
    root_level = root_logger.level
    root_handlers = list(root_logger.handlers)
    qlms_logger.setLevel(logging.WARNING)
    try:
        os.chdir(toml.parent)
        hwave.qlms.run(input_file=toml.name)
    finally:
        os.chdir(cwd)
        qlms_logger.setLevel(level)
        root_logger.setLevel(root_level)
        for handler in list(root_logger.handlers):
            if handler not in root_handlers:
                root_logger.removeHandler(handler)


def _build_uhfr_tutorial(dest):
    _copy_inputs(REPO / "docs" / "tutorial" / "Hubbard" / "UHFr", dest)
    write_greenone(dest / "greenone.def",
                   [(i, s, i, s) for i in range(8) for s in (0, 1)])
    toml = dest / "input.toml"
    _run(toml)
    return toml


def write_greenone(path, entries):
    """Write a UHFr OneBodyG file (five header lines, then ``i s j t`` per
    line) listing ``entries``."""
    rows = ["{:5d} {:5d} {:5d} {:5d}".format(*e) for e in entries]
    header = ["=" * 31, "NCisAjs {:9d}".format(len(rows)), "=" * 31,
              "======== Green functions ======", "=" * 31]
    Path(path).write_text("\n".join(header + rows) + "\n")


def run_input(toml):
    """Run H-wave on ``toml`` as the fixtures do (for tests that rerun a
    copied fixture)."""
    _run(Path(toml))


def _build_uhfr_u8(dest):
    _copy_inputs(REPO / "tests" / "uhfr" / "coulombintra", dest)
    write_greenone(dest / "greenone.def",
                   [(i, s, i, t) for i in range(16) for s in (0, 1) for t in (0, 1)])
    toml = dest / "input.toml"
    _run(toml)
    return toml


def _build_uhfr_tutorial_zero_steps(dest):
    _copy_inputs(REPO / "docs" / "tutorial" / "Hubbard" / "UHFr", dest)
    toml = dest / "input.toml"
    text = toml.read_text()
    if "IterationMax = 1000" not in text:
        raise RuntimeError("the UHFr tutorial input changed: no IterationMax = 1000")
    toml.write_text(text.replace("IterationMax = 1000", "IterationMax = 0"))
    _run(toml)
    return toml


def _build_uhfk_cdw(dest):
    _copy_inputs(REPO / "docs" / "tutorial" / "Hubbard" / "UHFk", dest)
    toml = dest / "input.toml"
    _run(toml)
    return toml


def _build_uhfk_cdw_legacy(dest):
    shutil.copytree(REPO / "docs" / "tutorial" / "Hubbard" / "UHFk", dest,
                    ignore=shutil.ignore_patterns("run.log"))
    return dest / "input.toml"


def _uhfk_u8_input(dest, **param):
    """tests/uhfk/CoulombIntra copied to ``dest`` with ``print_check`` set
    and ``param`` written into ``[mode.param]``; returns (toml, config)."""
    import tomli
    _copy_inputs(REPO / "tests" / "uhfk" / "CoulombIntra", dest)
    toml = dest / "input.toml"
    with open(toml, "rb") as f:
        config = tomli.load(f)
    config["log"]["print_level"] = 0
    config["log"]["print_check"] = "check.dat"
    config["mode"]["param"].update(param)
    return toml, config


def _write_toml(toml, config):
    """Write the two-level tables of an H-wave input (no dependency on a
    TOML writer)."""
    def value(v):
        if isinstance(v, bool):
            return "true" if v else "false"
        if isinstance(v, str):
            return '"{}"'.format(v)
        if isinstance(v, list):
            return "[{}]".format(", ".join(value(x) for x in v))
        return repr(v)

    lines = []

    def table(name, mapping):
        lines.append("[{}]".format(name))
        subtables = []
        for k, v in mapping.items():
            if isinstance(v, dict):
                subtables.append((k, v))
            else:
                lines.append("  {} = {}".format(k, value(v)))
        for k, v in subtables:
            table("{}.{}".format(name, k), v)

    for name, mapping in config.items():
        table(name, mapping)
    toml.write_text("\n".join(lines) + "\n")


def _build_uhfk_u8(dest):
    toml, config = _uhfk_u8_input(dest)
    _write_toml(toml, config)
    _run(toml)
    return toml


def _build_uhfk_u8_not_converged(dest):
    toml, config = _uhfk_u8_input(dest, IterationMax=5)
    finput = config["file"]["input"]
    finput.pop("initial", None)
    finput["initial_mode"] = "random"
    _write_toml(toml, config)
    _run(toml)
    # the start really was random: its first step differs from that of a
    # zero start (what a misplaced initial_mode silently gives) and from the
    # committed start (output_ref), which converges in 1 step at E = -7.38962582
    first = _first_check_line(dest)
    zero_dest = dest.parent / (dest.name + "-zero-start")
    zero_toml, zero_config = _uhfk_u8_input(zero_dest, IterationMax=1)
    zero_config["file"]["input"].pop("initial", None)
    _write_toml(zero_toml, zero_config)
    _run(zero_toml)
    zero = _first_check_line(zero_dest)
    if first[1:] == zero[1:] or abs(float(first[2]) + 7.38962582) < 1e-3:
        raise RuntimeError("uhfk_u8_not_converged did not start from a random "
                           "Green function: first step {} (zero start: {})"
                           .format(first, zero))
    return toml


def _first_check_line(dest):
    """The fields of the first line of ``output/check.dat`` under ``dest``."""
    text = (dest / "output" / "check.dat").read_text()
    return [field.strip() for field in text.splitlines()[0].split(",")]


def _build_uhfk_u8_zero_steps(dest):
    toml, config = _uhfk_u8_input(dest, IterationMax=0)
    _write_toml(toml, config)
    _run(toml)
    return toml


def _build_uhfk_u8_negative_max(dest):
    toml, config = _uhfk_u8_input(dest, IterationMax=-3,
                                  trustme_interaction_range=True)
    _write_toml(toml, config)
    _run(toml)
    return toml


def _build_uhfk_u8_nan_residual(dest):
    from hwave.solver import uhfk
    toml, config = _uhfk_u8_input(dest, IterationMax=5)
    finput = config["file"]["input"]
    finput.pop("initial", None)
    finput["initial_mode"] = "random"
    _write_toml(toml, config)
    orig = uhfk.UHFk._calc_phys
    calls = []

    def patched(this, *args, **kwargs):
        result = orig(this, *args, **kwargs)
        calls.append(1)
        if len(calls) >= 3:
            this.physics["Rest"] = float("nan")
        return result

    uhfk.UHFk._calc_phys = patched
    try:
        _run(toml)
    finally:
        uhfk.UHFk._calc_phys = orig
    return toml


def _build_uhfk_fm_axis(dest):
    toml, config = _uhfk_u8_input(dest, Ncond=5, SubShape=[1, 1, 1])
    n = 5 / 16
    u = np.array([0.6, 0.2537, 0.7734])
    u = u / np.linalg.norm(u)
    sigma = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]),
             np.array([[1, 0], [0, -1]])]
    rho = (n / 2) * (np.eye(2) + sum(c * s for c, s in zip(u, sigma)))
    green = np.zeros((16, 2, 1, 2, 1), dtype=complex)
    green[0, :, 0, :, 0] = rho.T        # G_st = <c^dag_s c_t> = rho_ts
    np.savez(dest / "start.npz", green=green, momentum_convention="e_plus_ikR")
    config["file"]["input"]["initial"] = "start.npz"
    _write_toml(toml, config)
    _run(toml)
    return toml


_BUILDERS = {
    "uhfr_tutorial": _build_uhfr_tutorial,
    "uhfr_u8": _build_uhfr_u8,
    "uhfr_tutorial_zero_steps": _build_uhfr_tutorial_zero_steps,
    "uhfk_cdw": _build_uhfk_cdw,
    "uhfk_cdw_legacy": _build_uhfk_cdw_legacy,
    "uhfk_u8": _build_uhfk_u8,
    "uhfk_u8_not_converged": _build_uhfk_u8_not_converged,
    "uhfk_u8_zero_steps": _build_uhfk_u8_zero_steps,
    "uhfk_u8_negative_max": _build_uhfk_u8_negative_max,
    "uhfk_u8_nan_residual": _build_uhfk_u8_nan_residual,
    "uhfk_fm_axis": _build_uhfk_fm_axis,
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
    if isinstance(obj, Mapping):
        for key, value in obj.items():
            assert_deep_frozen(value, "{}[{!r}]".format(where, key))
        return
    if isinstance(obj, (list, tuple)):
        for i, value in enumerate(obj):
            assert_deep_frozen(value, "{}[{}]".format(where, i))
