"""Open an H-wave run: its input TOML, the base directory, the output files.

:func:`open_run` reads an input TOML file and returns a :class:`Run`, which
knows where the solver wrote each of its outputs. H-wave resolves
``path_to_input`` and ``path_to_output`` against the directory it was started
from (qlms.py:108-114), which the TOML file does not record. The base
directory is therefore, in this order:

1. ``base_dir=`` when it is given (``base_dir_source == "argument"``);
2. the ``cwd`` of the run record ``hwave-run.<stem>.json`` next to the TOML
   file (``"run_record"``);
3. the directory of the TOML file (``"toml_directory"``).

Output names follow the solvers exactly: ``[file.output]`` and
``[log] print_check`` keys are matched with their exact case, as the solvers
read them (uhfk.py:1796-1897, uhfr.py:1271-1375, uhfk.py:816-823,
uhfr.py:575-577); a key that differs only in case is not written by the
solver, so it is not configured here, and it is listed in
``Run.meta["ignored_output_keys"]``.

:attr:`Run.config` reads each layer of the TOML file the way the solver
does:

* exact-case: the top-level tables ``[mode]``, ``[file]``, ``[log]``
  (qlms.py:98-111), the keys of ``[mode]`` (qlms.py:20-53: a key in other
  case makes the solver stop), of ``[file]`` (qlms.py:104-111), of
  ``[file.input]`` (UHFk: read_input_k.py:152-190; UHFr: ``path_to_input``
  and ``interaction``, qlms.py:108, 146), of ``[file.output]`` and of
  ``[log]``;
* case-insensitive: ``[mode.param]`` (``ParamDict``, base.py:16,
  qlms.py:137) and ``[file.input.interaction]`` (qlms.py:146-149;
  read_input_k.py:93). UHFr also matches ``[file.input]`` ``initial`` and
  ``onebodyg`` in any case (qlms.py:155-159); in :attr:`Run.config` that
  table stays exact-case, so look these two up case-insensitively.

A key of an exact-case layer that matches a key the solver reads only when
case is ignored is not read by the solver. It is listed in
``Run.meta["ignored_keys"]`` (which also holds the ``ignored_output_keys``)
and logged once as a WARNING. A ``[Mode]`` table in other case is an
UnsupportedLayoutError, as the solver stops without a ``[mode]`` table.
"""
from __future__ import annotations

import json
import logging
import os
from collections.abc import Mapping
from types import MappingProxyType
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import tomli
from requests.structures import CaseInsensitiveDict

import hwave
from hwave.results.errors import (
    MissingConfigError, MissingOutputError, RunRecordError,
    UnsupportedLayoutError)

logger = logging.getLogger(__name__)

__all__ = ["Run", "open_run", "as_run", "read_run_record", "write_run_record",
           "run_record_path"]

RUN_RECORD_FORMAT = "hwave-run-record"
RUN_RECORD_VERSION = 1

#: The ``[mode.param]`` defaults of the UHF solvers, as set in
#: ``solver_base.__init__`` (src/hwave/solver/base.py:146-156, ``para_init``).
#: ``EPS`` is the exponent: the solver's threshold is ``10**-EPS``
#: (base.py:237-253). ``2Sz = "free"`` means ``None`` (base.py:194-195).
#: tests/test_results_run.py checks these against the solver.
SOLVER_PARAM_DEFAULTS = {
    "Nsite": 0,
    "Ncond": 0,
    "2Sz": None,
    "Mix": 0.5,
    "EPS": 6,
    "IterationMax": 20000,
    "RndSeed": 1234,
    "T": 0.0,
    "ene_cutoff": 1.0e+2,
}

#: marks "no default given" in Run.param
_NO_DEFAULT = object()

#: modes whose outputs are resolved
SUPPORTED_MODES = ("UHFr", "UHFk")

#: ``[file.output]`` keys each solver reads (exact case)
_OUTPUT_KEYS = {
    "UHFr": ("path_to_output", "energy", "eigen", "green", "initial", "fij",
             "scf_status"),
    "UHFk": ("path_to_output", "energy", "eigen", "green", "onebodyg", "initial",
             "rpa", "export_hamiltonian", "scf_status"),
}

_DIRECTORY_HINT = "pass the input TOML file; directory discovery is not part of this API"
_BASE_DIR_HINT = ("if H-wave was started from another directory, pass "
                  "base_dir=<that directory>")


def _npz(name):
    """The file name ``numpy.savez`` writes for ``name`` (base.py:90-93)."""
    return name if name.endswith(".npz") else name + ".npz"


def _case_insensitive(obj):
    """A deep copy of ``obj`` in which every mapping is a CaseInsensitiveDict."""
    if isinstance(obj, Mapping):
        return CaseInsensitiveDict({k: _case_insensitive(v) for k, v in obj.items()})
    if isinstance(obj, list):
        return [_case_insensitive(v) for v in obj]
    return obj


def _exact_copy(obj):
    """A deep copy of ``obj`` with plain dicts for the mappings."""
    if isinstance(obj, Mapping):
        return {k: _exact_copy(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_exact_copy(v) for v in obj]
    return obj


def _solver_view(raw):
    """``raw`` with the solver's case rules: CaseInsensitiveDicts at
    ``[mode.param]`` and ``[file.input.interaction]``, plain dicts
    elsewhere."""
    config = _exact_copy(raw)
    mode = config.get("mode")
    if isinstance(mode, dict) and isinstance(mode.get("param"), dict):
        mode["param"] = _case_insensitive(mode["param"])
    finput = config.get("file", {}).get("input") \
        if isinstance(config.get("file"), dict) else None
    if isinstance(finput, dict) and isinstance(finput.get("interaction"), dict):
        finput["interaction"] = _case_insensitive(finput["interaction"])
    return config


#: the top-level tables qlms.run reads (qlms.py:98-104), exact-case
_TOP_TABLES = ("mode", "file", "log")
#: the keys of [file] qlms.run reads (qlms.py:104-111)
_FILE_KEYS = ("input", "output")
#: the exact-case [file.input] keys of each solver: UHFk read_input_k.py:93,
#: 152-190; UHFr qlms.py:108, 146 (UHFr matches initial and onebodyg in any
#: case, qlms.py:155-159, so they are not listed)
_INPUT_KEYS = {
    "UHFk": ("path_to_input", "interaction", "initial", "initial_uhf",
             "initial_mode", "geometry_uhf", "onebodyg_uhf"),
    "UHFr": ("path_to_input", "interaction"),
}


def _case_variants(table, keys, where, is_table=False):
    """The entries of ``table`` that equal one of ``keys`` only when case
    is ignored."""
    lower = {k.lower() for k in keys}
    found = []
    for key in table:
        if key not in keys and key.lower() in lower:
            found.append("[{}]".format(key) if is_table else "{} {}".format(where, key))
    return found


def _ignored_input_keys(raw, mode):
    """The keys of the exact-case layers above [file.output] and [log]
    that the solver does not read because of their case."""
    from hwave.qlms import MODE_TABLE_KEYS
    ignored = _case_variants(raw, _TOP_TABLES, "", is_table=True)
    info_mode = raw.get("mode")
    if isinstance(info_mode, dict):
        ignored += _case_variants(info_mode, MODE_TABLE_KEYS, "[mode]")
    info_file = raw.get("file")
    if isinstance(info_file, dict):
        ignored += _case_variants(info_file, _FILE_KEYS, "[file]")
        info_input = info_file.get("input")
        if isinstance(info_input, dict) and mode in _INPUT_KEYS:
            ignored += _case_variants(info_input, _INPUT_KEYS[mode], "[file.input]")
    return ignored


@dataclass(frozen=True)
class Run:
    """An H-wave run: where its input is, and where its outputs are.

    Attributes
    ----------
    root : Path
        The absolute directory against which ``path_to_input`` and
        ``path_to_output`` are resolved (the directory H-wave ran in).
    base_dir_source : str
        How ``root`` was chosen: ``"argument"``, ``"run_record"`` or
        ``"toml_directory"``.
    run_record : dict or None
        The run record that gave ``root``, if any.
    input_path : Path
        The canonical path of the input TOML file.
    config : dict
        The TOML content, with the case rules of the solver (module text):
        ``config["mode"]["param"]`` and ``config["file"]["input"]
        ["interaction"]`` are CaseInsensitiveDicts, every other mapping a
        plain dict.
    mode : str
        ``[mode] mode`` (e.g. ``"UHFk"``).
    spin_orbital : bool
        ``[mode] enable_spin_orbital``, parsed as the solvers do.
    output_dir : Path
        ``root / path_to_output``.
    outputs : Mapping[str, Path]
        Read-only. Every configured output and the file the solver writes for it.
        UHFr eigenvectors are per block: ``eigen:sz-free``, or
        ``eigen:spin-up`` and ``eigen:spin-down``. Empty for other modes.
    meta : Mapping
        Read-only. Notes: ``run_record_path``, ``run_record_note`` (a record bypassed by
        ``base_dir=``), ``unsupported_mode``, ``ignored_output_keys`` (keys
        in other case, which the solver does not read), ``ignored_keys``
        (every such key of an exact-case layer, these included) and
        ``invalid_output_keys`` (a value that is not a file name).
    """
    root: Path
    base_dir_source: str
    run_record: Optional[dict]
    input_path: Path
    config: dict
    mode: str
    spin_orbital: bool
    output_dir: Path
    outputs: Mapping
    meta: Mapping

    def path(self, key):
        """The file of output ``key``. Raises MissingOutputError when the
        output is not configured, or when its file does not exist."""
        if key not in self.outputs:
            configured = ", ".join(sorted(self.outputs)) or "none"
            blocks = [k for k in self.outputs if k.startswith("eigen:")]
            if "unsupported_mode" in self.meta:
                hint = ("outputs of mode {!r} are not resolved by this version of "
                        "hwave.results".format(self.mode))
            elif key.startswith("eigen:") and blocks:
                hint = ("the eigenvector blocks follow [mode.param] 2Sz; this run "
                        "has {}".format(", ".join(blocks)))
            else:
                hint = "configure it in [file.output] of the input file"
            raise MissingOutputError(
                "output {!r} is not configured in {} (configured outputs: {})"
                .format(key, self.input_path, configured),
                path=self.input_path, hint=hint)
        path = self.outputs[key]
        if not path.is_file():
            raise MissingOutputError(
                "output {!r}: file {} does not exist; outputs were resolved against "
                "{} (base directory source: {})".format(
                    key, path, self.root, self.base_dir_source),
                path=path, hint=_BASE_DIR_HINT)
        return path

    def param(self, name, default=_NO_DEFAULT):
        """``[mode.param] name``, looked up case-insensitively. When the key
        is absent, ``default`` is returned if the caller gives one, else the
        solver default (:data:`SOLVER_PARAM_DEFAULTS`), else ``None``.
        ``2Sz = "free"`` gives ``None``, as in the solver. ``EPS`` is
        returned as written (the exponent)."""
        params = self.config.get("mode", {}).get("param", {})
        if name in params:
            value = params[name]
        elif default is not _NO_DEFAULT:
            value = default
        else:
            value = CaseInsensitiveDict(SOLVER_PARAM_DEFAULTS).get(name)
        if name.lower() == "2sz" and value == "free":
            value = None
        return value


def run_record_path(toml):
    """The run record of ``toml``: ``hwave-run.<stem>.json`` next to it."""
    toml = Path(toml)
    return toml.parent / "hwave-run.{}.json".format(toml.stem)


def write_run_record(toml, *, cwd, argv):
    """Write the run record of ``toml`` (see :func:`run_record_path`): the
    directory H-wave runs in (``cwd``, made absolute), the canonical TOML
    path, the command line and the H-wave version. Returns its path."""
    toml = Path(toml)
    if not toml.is_file():
        raise MissingConfigError("input TOML file not found: {}".format(toml), path=toml)
    record = {
        "format": RUN_RECORD_FORMAT,
        "format_version": RUN_RECORD_VERSION,
        "cwd": str(Path(cwd).resolve()),
        "input": str(toml.resolve()),
        "argv": [os.fspath(a) for a in argv],
        "hwave_version": hwave.__version__,
    }
    path = run_record_path(toml)
    with open(path, "w") as f:
        json.dump(record, f, indent=2)
        f.write("\n")
    return path


def read_run_record(toml):
    """The run record of ``toml``, or ``None`` when there is none.

    Raises RunRecordError when the record is not a readable
    ``hwave-run-record`` version 1, or when it was written for another
    input file (e.g. after the directory was moved)."""
    toml = Path(toml)
    path = run_record_path(toml)
    if not path.exists():
        return None
    try:
        with open(path, encoding="utf-8") as f:
            record = json.load(f)
    except (OSError, ValueError) as exc:
        raise RunRecordError("cannot read the run record {}: {}".format(path, exc),
                             path=path) from None

    def bad(what):
        return RunRecordError("run record {}: {}".format(path, what), path=path,
                              hint="pass base_dir= to open the run without it")

    if not isinstance(record, dict):
        raise bad("not a JSON object")
    if record.get("format") != RUN_RECORD_FORMAT:
        raise bad("unknown format {!r} (expected {!r})".format(
            record.get("format"), RUN_RECORD_FORMAT))
    version = record.get("format_version")
    if type(version) is not int or version != RUN_RECORD_VERSION:
        raise bad("unsupported format_version {!r} (this reader knows {})".format(
            version, RUN_RECORD_VERSION))
    for key in ("cwd", "input"):
        value = record.get(key)
        if not isinstance(value, str) or not os.path.isabs(value):
            raise bad("{!r} must be an absolute path, got {!r}".format(key, value))
    argv = record.get("argv")
    if not isinstance(argv, list) or not all(isinstance(a, str) for a in argv):
        raise bad("'argv' must be a list of strings")
    if not isinstance(record.get("hwave_version"), str):
        raise bad("'hwave_version' must be a string")
    if not _same_file(Path(record["input"]), toml):
        raise RunRecordError(
            "run record {} was written for {}, not for {}; the files were moved "
            "or copied after the run".format(path, Path(record["input"]).resolve(),
                                             toml.resolve()),
            path=path,
            hint="pass base_dir=<the directory the run's paths are relative to> "
                 "to open the run without the record")
    return dict(record)


def _same_file(a, b):
    """Whether ``a`` and ``b`` name the same file: ``os.path.samefile`` when
    both exist (so names that differ only in case match on a
    case-insensitive file system), else the resolved paths."""
    if a.exists() and b.exists():
        return os.path.samefile(a, b)
    return a.resolve() == b.resolve()


def _table(parent, key, where, toml):
    value = parent.get(key, {})
    if not isinstance(value, dict):
        raise UnsupportedLayoutError("{} is not a table in {}".format(where, toml),
                                     path=toml)
    return value


def _resolve_outputs(mode, info_output, info_log, two_sz, output_dir, meta):
    """``{key: Path}`` of the files the ``mode`` solver writes."""
    known = _OUTPUT_KEYS[mode]
    ignored = []
    for where, table, keys in (("[file.output]", info_output, known),
                               ("[log]", info_log, ("print_check",))):
        lower = {k.lower(): k for k in keys}
        for key in table:
            if key not in keys and key.lower() in lower:
                ignored.append("{} {}".format(where, key))
    for entry in ignored:
        logger.warning("%s is not read by the solver: keys are matched with their "
                       "exact case", entry)
    invalid = []

    def name(table, where, key, default=None):
        value = table.get(key, default)
        if value is None:
            return None
        if not isinstance(value, str):
            invalid.append("{} {}".format(where, key))
            return None
        return value or None          # "" is not configured

    def out(key, default=None):
        return name(info_output, "[file.output]", key, default)

    outputs = {}
    if out("energy"):
        outputs["energy"] = output_dir / out("energy")
    if mode == "UHFr":
        # uhfr.py:1296-1375; the eigen and fij names are per spin block
        # (uhfr.py:553-560)
        blocks = ["sz-free"] if two_sz is None else ["spin-up", "spin-down"]
        if out("eigen"):
            for block in blocks:
                outputs["eigen:" + block] = output_dir / _npz(block + "_" + out("eigen"))
        if out("green"):
            outputs["green"] = output_dir / out("green")
        if out("initial"):
            outputs["initial"] = output_dir / out("initial")
        if out("fij"):
            prefix = "sz-free_" if two_sz is None else ("Sz0_" if two_sz == 0 else None)
            if prefix is not None:
                outputs["fij"] = output_dir / (prefix + out("fij"))
    else:
        # uhfk.py:1839-1895; green defaults to "green" (uhfk.py:1879)
        if out("eigen"):
            outputs["eigen"] = output_dir / _npz(out("eigen"))
        if out("green", "green"):
            outputs["green"] = output_dir / _npz(out("green", "green"))
        if out("onebodyg"):
            outputs["onebodyg"] = output_dir / out("onebodyg")
        if out("rpa"):
            outputs["rpa"] = output_dir / _npz(out("rpa"))
    # base.py:77-87: on by default, "" disables it
    if out("scf_status", "scf_status.dat"):
        outputs["scf_status"] = output_dir / out("scf_status", "scf_status.dat")
    check = name(info_log, "[log]", "print_check")
    if check:
        outputs["print_check"] = output_dir / check
    if ignored:
        meta["ignored_output_keys"] = tuple(ignored)
    if invalid:
        meta["invalid_output_keys"] = tuple(invalid)
    return outputs


def open_run(toml_path, *, base_dir=None):
    """Open the run of the input TOML file ``toml_path``.

    ``base_dir`` is the directory H-wave was started from; see the module
    text for the default. Raises MissingConfigError for a missing file or a
    directory, UnsupportedLayoutError for a file that is not an H-wave input
    (no ``[mode]`` table), and RunRecordError for an unusable run record
    (``base_dir=`` bypasses the record)."""
    toml = Path(os.fspath(toml_path))
    if toml.is_dir():
        raise MissingConfigError("{} is a directory".format(toml), path=toml,
                                 hint=_DIRECTORY_HINT)
    if not toml.is_file():
        raise MissingConfigError("input TOML file not found: {}".format(toml), path=toml)
    try:
        with open(toml, "rb") as f:
            raw = tomli.load(f)
    except (tomli.TOMLDecodeError, UnicodeDecodeError) as exc:
        raise UnsupportedLayoutError("{} is not a valid TOML file: {}".format(toml, exc),
                                     path=toml) from None

    # the tables and keys that place the files are read as qlms.run reads them
    info_mode = raw.get("mode")
    if not isinstance(info_mode, dict):
        variants = [k for k in raw if k.lower() == "mode" and isinstance(raw[k], dict)]
        if variants:
            raise UnsupportedLayoutError(
                "{} has no [mode] table: [{}] is not read, the solver matches the "
                "table names with their exact case (qlms.py:103, 122-124)"
                .format(toml, variants[0]), path=toml)
        raise UnsupportedLayoutError("{} has no [mode] table; it is not an H-wave input "
                                     "file".format(toml), path=toml)
    mode = info_mode.get("mode")
    if not isinstance(mode, str):
        raise UnsupportedLayoutError("[mode] mode must be a string in {}, got {!r}"
                                     .format(toml, mode), path=toml)
    info_param = _table(info_mode, "param", "[mode.param]", toml)
    info_file = _table(raw, "file", "[file]", toml)
    info_output = _table(info_file, "output", "[file.output]", toml)
    info_log = _table(raw, "log", "[log]", toml)
    path_to_output = info_output.get("path_to_output", "output")
    if not isinstance(path_to_output, str):
        raise UnsupportedLayoutError("[file.output] path_to_output must be a string in {}"
                                     .format(toml), path=toml)

    from hwave.solver.spin_orbital import resolve_flag
    try:
        # [mode] keys are exact-case for the solver (qlms.py:26-53), which
        # resolve_flag itself does not enforce
        spin_orbital = resolve_flag(
            {k: v for k, v in info_mode.items() if k == "enable_spin_orbital"},
            allow_float=(mode == "UHFk"))
    except (TypeError, ValueError) as exc:
        raise UnsupportedLayoutError("[mode] enable_spin_orbital in {}: {}".format(toml, exc),
                                     path=toml) from None

    meta = {}
    record = None
    if base_dir is not None:
        root = Path(os.fspath(base_dir)).resolve()
        source = "argument"
        record_file = run_record_path(toml)
        if record_file.exists():
            meta["run_record_note"] = (
                "the run record {} was not used because base_dir= was given"
                .format(record_file))
    else:
        record = read_run_record(toml)
        if record is not None:
            root = Path(record["cwd"])
            source = "run_record"
            meta["run_record_path"] = str(run_record_path(toml))
        else:
            root = toml.parent.resolve()
            source = "toml_directory"
    logger.debug("run %s: base directory %s (source: %s)", toml, root, source)

    output_dir = root / path_to_output
    ignored = _ignored_input_keys(raw, mode)
    for entry in ignored:
        logger.warning("%s is not read by the solver: keys are matched with their "
                       "exact case", entry)
    if mode in SUPPORTED_MODES:
        two_sz = CaseInsensitiveDict(info_param).get("2Sz")
        if two_sz == "free":
            two_sz = None
        outputs = _resolve_outputs(mode, info_output, info_log, two_sz, output_dir, meta)
    else:
        outputs = {}
        meta["unsupported_mode"] = mode
    ignored += list(meta.get("ignored_output_keys", ()))
    if ignored:
        meta["ignored_keys"] = tuple(ignored)

    return Run(root=root, base_dir_source=source, run_record=record,
               input_path=toml.resolve(), config=_solver_view(raw), mode=mode,
               spin_orbital=spin_orbital, output_dir=output_dir,
               outputs=MappingProxyType(outputs), meta=MappingProxyType(meta))


def as_run(run_or_toml, *, base_dir=None):
    """``run_or_toml`` as a :class:`Run`: a Run is returned as is, anything
    else is opened with :func:`open_run`. ``base_dir`` cannot be combined
    with a Run (its base directory was fixed when it was opened): that
    raises MissingConfigError."""
    if isinstance(run_or_toml, Run):
        if base_dir is not None:
            raise MissingConfigError(
                "base_dir= cannot be given with a Run; its base directory was fixed "
                "when it was opened", path=run_or_toml.input_path,
                hint="call open_run(toml, base_dir=...) to open the run with another "
                     "base directory")
        return run_or_toml
    return open_run(run_or_toml, base_dir=base_dir)
