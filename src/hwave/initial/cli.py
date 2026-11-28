"""The ``python -m hwave.initial`` commands ``prepare``, ``check`` and
``collect`` (preview).

Workflow: ``prepare`` writes one run directory per candidate run under
WORK, the user runs ``hwave`` in each (``bash WORK/commands.txt``), and
``collect`` reads the results back into ``WORK/coverage.json`` and
``WORK/summary.txt``.

**Run directories.** :func:`stage_run_dir` writes every run directory (and
nothing else writes an ``input.toml``): the input TOML is rewritten
(``[file.input] initial = "initial.npz"``, ``[mode.param] SubShape``, the
input files copied into ``inputs/``), written by a small serialiser for the
TOML subset H-wave reads, re-parsed and checked with the solver's input
reader. Comments of the original TOML are not kept.

**Relocatable WORK.** WORK holds no symbolic link, and every path that
running or collecting uses is relative to WORK or to its run directory, so
WORK can be moved, or its sources deleted, and still be run and collected.
The original TOML's absolute path is kept in ``manifest.json`` for
``prepare`` re-entry and display only.

**The preliminary paramagnetic run** (``prepare --stage0``, small tier):
the paramagnetic Hartree-Fock step that the small tier otherwise runs
inside ``prepare`` is written as run directories of its own,
``WORK/stage0_para/`` (and ``WORK/stage0_para_check/`` with
``--para-check``), staged as the candidate runs are (``SubShape = [1, 1,
1]``, ``T = T_hint``, ``EPS = max(EPS, 12)``, ``[file.output] green``
set), with ``WORK/stage0_commands.txt``; ``prepare`` exits with code 3.
After those runs, the same ``prepare`` command checks their inputs against
``manifest.json`` and reads their results: when one converged it becomes
the para reference (``WORK/para_hf/``) and ``prepare`` writes the candidate
runs; when it ran and failed the seeds use the bare ``H0`` and the hints
are skipped with the reason; when it has not run yet, ``prepare`` prints
the commands again and exits with code 3.

**Exit codes** (:data:`EXIT_CODES`): 0 done (for ``check``: verified, or
unverified without ``--require-provenance``); 1 an error or a refusal; 2 an
argument error; 3 a preliminary paramagnetic run is pending (``prepare
--stage0``); 4 ``check --require-provenance`` on an initial-state file
without provenance.
"""
import argparse
import copy
import dataclasses
import hashlib
import json
import math
import numbers
import os
import pathlib
import re
import secrets
import shutil
import sys
from dataclasses import dataclass
from typing import Mapping, Optional, Tuple

import numpy as np

__all__ = ["main", "prepare", "check", "collect", "stage_run_dir", "StagedRun",
           "StagedFile", "dumps_toml", "EXIT_CODES", "EXIT_OK", "EXIT_ERROR",
           "EXIT_USAGE", "EXIT_STAGE0_PENDING", "EXIT_NO_PROVENANCE", "MANIFEST_VERSION"]

EXIT_OK = 0
EXIT_ERROR = 1
EXIT_USAGE = 2
EXIT_STAGE0_PENDING = 3
EXIT_NO_PROVENANCE = 4

#: the exit codes, as listed in ``--help``
EXIT_CODES = (
    (EXIT_OK, "done (check: verified, or unverified without --require-provenance)"),
    (EXIT_ERROR, "error, or the check refused the file"),
    (EXIT_USAGE, "argument error"),
    (EXIT_STAGE0_PENDING, "prepare --stage0: a preliminary paramagnetic run is pending; "
                          "run bash WORK/stage0_commands.txt, then re-run the same "
                          "prepare command"),
    (EXIT_NO_PROVENANCE, "check --require-provenance: the initial-state file has no "
                         "provenance"),
)

#: the version of ``manifest.json``
MANIFEST_VERSION = 1

#: the names inside WORK and each run directory
MANIFEST = "manifest.json"
COVERAGE = "coverage.json"
CANDIDATES = "candidates.json"
COMMANDS = "commands.txt"
SUMMARY = "summary.txt"
RUN_TOML = "input.toml"
RUN_INITIAL = "initial.npz"
RUN_INPUTS = "inputs"
RUN_LOG = "hwave.log"
STAGE0_COMMANDS = "stage0_commands.txt"

#: ``manifest.json`` ``stage``: the preliminary paramagnetic run is written
#: and has not been read back yet, or prepare is done
STAGE_PENDING = "para_pending"
STAGE_DONE = "done"
#: the directories of the preliminary paramagnetic run: ``(directory, the
#: start of :func:`hwave.initial.hints.para_starts`, the candidate id of
#: its start file)``; the second only with ``--para-check``
STAGE0_RUNS = (("stage0_para", "main", "para-hf-start"),
               ("stage0_para_check", "check", "para-hf-check-start"))
#: ``[file.output] green`` of a preliminary paramagnetic run when the input
#: names none
STAGE0_GREEN = "green"
#: the hint printed with the commands of the preliminary paramagnetic run
PENDING_HINT = "re-run the same prepare command after these finish"

#: ``[file.input]`` keys whose files are staged (output requests of UHFk)
STAGED_INPUT_KEYS = ("geometry_uhf", "onebodyg_uhf")
#: ``[file.input]`` keys removed from a run's input (an initial state)
REMOVED_INPUT_KEYS = ("initial", "initial_uhf", "initial_mode")
#: ``[file.output]`` keys that name a file under ``path_to_output``
OUTPUT_FILE_KEYS = ("energy", "eigen", "green", "onebodyg", "initial", "rpa",
                    "export_hamiltonian", "scf_status")
#: written when ``[file.output] energy`` is absent: UHFk writes energy.dat
#: only when the key is present, and collect ranks on it
DEFAULT_ENERGY = "energy.dat"

_HEADER = ("# Written by python -m hwave.initial prepare from {}.\n"
           "# Comments of the original input are not kept.\n")

_BARE_KEY = re.compile(r"[A-Za-z0-9_-]+")


class _Usage(Exception):
    """An argument error (exit code 2)."""


# --------------------------------------------------------------------------
# the TOML serialiser
# --------------------------------------------------------------------------

def _toml_string(s):
    # JSON escapes are TOML basic-string escapes; DEL must be escaped too
    return json.dumps(s, ensure_ascii=False).replace("\x7f", "\\u007f")


def _toml_key(k):
    return k if _BARE_KEY.fullmatch(k) else _toml_string(k)


def _toml_scalar(v, where):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, numbers.Integral):
        return str(int(v))
    if isinstance(v, numbers.Real):
        f = float(v)
        if not math.isfinite(f):
            raise ValueError("{}: {!r} is not a finite number; the input writer of "
                             "hwave.initial does not write it".format(where, v))
        return repr(f)
    if isinstance(v, str):
        return _toml_string(v)
    raise ValueError("{}: a value of type {} is outside the TOML subset hwave.initial "
                     "writes (tables, strings, numbers, booleans, flat arrays)".format(
                         where, type(v).__name__))


def _toml_value(v, where):
    if isinstance(v, (list, tuple)):
        items = []
        for i, x in enumerate(v):
            if isinstance(x, (list, tuple, dict)):
                raise ValueError("{}[{}]: nested arrays and arrays of tables are outside "
                                 "the TOML subset hwave.initial writes".format(where, i))
            items.append(_toml_scalar(x, "{}[{}]".format(where, i)))
        return "[" + ", ".join(items) + "]"
    return _toml_scalar(v, where)


def dumps_toml(doc) -> str:
    """``doc`` as TOML text: tables (nested dicts), strings, numbers,
    booleans and flat arrays of them, the subset of TOML H-wave inputs use.
    Anything else raises ``ValueError`` naming the key."""
    out = []

    def table(path, tab):
        scalars = [(k, v) for k, v in tab.items() if not isinstance(v, dict)]
        subs = [(k, v) for k, v in tab.items() if isinstance(v, dict)]
        name = ".".join(_toml_key(k) for k in path)
        if path and (scalars or not subs):
            if out:
                out.append("")
            out.append("[{}]".format(name))
        for k, v in scalars:
            where = "[{}] {}".format(name, k) if path else k
            out.append("{} = {}".format(_toml_key(k), _toml_value(v, where)))
        for k, v in subs:
            table(path + (k,), v)

    if not isinstance(doc, dict):
        raise ValueError("the document must be a table")
    for k in doc:
        if not isinstance(k, str):
            raise ValueError("key {!r} is not a string".format(k))
    table((), doc)
    return "\n".join(out) + "\n"


# --------------------------------------------------------------------------
# run directories
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class StagedFile:
    """One input file copied into a run directory: ``table`` and ``key``
    of the TOML entry, ``source`` the file as the solver would open it
    from the original input's directory, ``name`` its path in the run
    directory (``inputs/<file>``), ``sha256`` and ``nbytes`` of the copy."""
    table: str
    key: str
    source: str
    name: str
    sha256: str
    nbytes: int


@dataclass(frozen=True)
class StagedRun:
    """A written run directory: its ``input.toml``, the ``initial.npz``
    the solver will read (written afterwards by the caller) and the
    staged input files."""
    toml_path: pathlib.Path
    initial_path: pathlib.Path
    staged: Tuple[StagedFile, ...]


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _table(doc, key, where):
    v = doc.get(key)
    if v is None:
        v = {}
        doc[key] = v
    if not isinstance(v, dict):
        raise ValueError("{} must be a table".format(where))
    return v


def _ci_duplicates(tab, where):
    seen = {}
    for k in tab:
        low = k.lower()
        if low in seen:
            raise ValueError("{}: {!r} and {!r} are the same key (the solver reads the keys "
                             "of this table in any case); keep one".format(
                                 where, seen[low], k))
        seen[low] = k


def _ci_key(tab, key):
    for k in tab:
        if k.lower() == key.lower():
            return k
    return None


def _parts(path):
    return pathlib.PurePath(path).parts


def _check_relative(value, where):
    if not isinstance(value, str):
        raise ValueError("{} must be a string, got {!r}".format(where, value))
    if os.path.isabs(value):
        raise ValueError("{} = {!r} is an absolute path: the runs of one WORK would write "
                         "to the same place; give a path relative to the run directory"
                         .format(where, value))
    if ".." in _parts(value):
        raise ValueError("{} = {!r} contains '..': every run writes inside its own run "
                         "directory; give a path below it".format(where, value))


def _check_outputs(foutput, flog):
    pto = foutput.get("path_to_output", "output")
    _check_relative(pto, "[file.output] path_to_output")
    for key in OUTPUT_FILE_KEYS:
        if key in foutput:
            _check_relative(foutput[key], "[file.output] {}".format(key))
    for key in ("green", "energy"):
        if foutput.get(key) == "":
            raise ValueError("[file.output] {} = \"\" disables the {} output, which collect "
                             "reads; remove the key or name a file".format(key, key))
    if "print_check" in flog:
        _check_relative(flog["print_check"], "[log] print_check")
    # nothing the solver writes may land on the run directory's own files
    reserved = {RUN_TOML.lower(), RUN_INITIAL.lower(), RUN_LOG.lower()}

    def clobbers(rel):
        rel = pathlib.PurePosixPath(os.path.normpath(rel).replace(os.sep, "/")).as_posix()
        low = rel.lower()
        return (low in reserved or low == RUN_INPUTS.lower()
                or low.startswith(RUN_INPUTS.lower() + "/"))

    if clobbers(pto) and os.path.normpath(pto) not in (".", ""):
        raise ValueError("[file.output] path_to_output = {!r} lies on the run directory's "
                         "{} or input copies; choose another directory".format(pto, RUN_INPUTS))
    targets = [("[file.output] {}".format(k), foutput[k]) for k in OUTPUT_FILE_KEYS
               if isinstance(foutput.get(k), str) and foutput.get(k)]
    if isinstance(flog.get("print_check"), str) and flog["print_check"]:
        targets.append(("[log] print_check", flog["print_check"]))
    for where, name in targets:
        for cand in (name, name if name.endswith(".npz") else name + ".npz"):
            rel = os.path.join(pto, cand)
            if clobbers(rel):
                raise ValueError("{} = {!r} with path_to_output = {!r} would overwrite the "
                                 "run directory's {}; choose another name".format(
                                     where, name, pto, os.path.normpath(rel)))
    return pto


def _unique(name, taken):
    stem, ext = os.path.splitext(name)
    out, n = name, 1
    while out.lower() in taken:
        n += 1
        out = "{}-{}{}".format(stem, n, ext)
    taken.add(out.lower())
    return out


def _plan_run_doc(doc, *, source_dir, overrides):
    """The rewritten document and the copies to make:
    ``(new_doc, [(table, key, source, name)])``; nothing is written."""
    new = copy.deepcopy(dict(doc))
    mode = _table(new, "mode", "[mode]")
    param = _table(mode, "param", "[mode.param]")
    _ci_duplicates(param, "[mode.param]")
    ovr = dict(overrides or {})
    green = ovr.pop("green", None)
    for key, value in ovr.items():
        k = _ci_key(param, key)
        param[k if k is not None else key] = value

    ffile = _table(new, "file", "[file]")
    finput = _table(ffile, "input", "[file.input]")
    foutput = _table(ffile, "output", "[file.output]")
    flog = new.get("log", {})
    if not isinstance(flog, dict):
        raise ValueError("[log] must be a table")
    if green is not None:
        foutput["green"] = green
    _check_outputs(foutput, flog)
    if "energy" not in foutput:
        foutput["energy"] = DEFAULT_ENERGY

    pti = finput.get("path_to_input", "")
    if not isinstance(pti, str):
        raise ValueError("[file.input] path_to_input must be a string, got {!r}".format(pti))
    inter = finput.get("interaction", {})
    if not isinstance(inter, dict):
        raise ValueError("[file.input.interaction] must be a table")
    _ci_duplicates(inter, "[file.input.interaction]")
    ipti_key = _ci_key(inter, "path_to_input")
    ipti = inter.get(ipti_key, "") if ipti_key is not None else ""
    if not isinstance(ipti, str):
        raise ValueError("[file.input.interaction] {} must be a string, got {!r}".format(
            ipti_key, ipti))

    copies, taken, by_source = [], set(), {}

    def stage(table, key, value, base):
        if not isinstance(value, str) or not value:
            raise ValueError("{} {} must name a file, got {!r}".format(table, key, value))
        # located as the solver locates it, from the directory hwave runs in
        source = os.path.join(str(source_dir), os.path.join(base, value))
        if not os.path.isfile(source):
            raise ValueError("{} {} = {!r} names {}, which is not a file".format(
                table, key, value, source))
        real = os.path.realpath(source)
        if real not in by_source:
            by_source[real] = _unique(os.path.basename(value), taken)
        copies.append((table, key, source, by_source[real]))
        return by_source[real]

    new_input = {}
    for key, value in finput.items():
        if key in REMOVED_INPUT_KEYS:
            continue
        if key == "path_to_input":
            new_input[key] = ""
        elif key in STAGED_INPUT_KEYS:
            new_input[key] = RUN_INPUTS + "/" + stage("[file.input]", key, value, pti)
        elif key == "interaction":
            new_inter = {}
            for ikey, ivalue in inter.items():
                if ikey == ipti_key:
                    new_inter[ikey] = RUN_INPUTS
                else:
                    new_inter[ikey] = stage("[file.input.interaction]", ikey, ivalue, ipti)
            if ipti_key is None:
                new_inter = dict([("path_to_input", RUN_INPUTS)] + list(new_inter.items()))
            new_input[key] = new_inter
        else:
            new_input[key] = value
    if "path_to_input" not in new_input:
        new_input = dict([("path_to_input", "")] + list(new_input.items()))
    new_input["initial"] = RUN_INITIAL
    ffile["input"] = new_input
    return new, copies


def _inside(path, root):
    path, root = os.path.realpath(path), os.path.realpath(root)
    return os.path.commonpath([path, root]) == root


def stage_run_dir(run_dir, doc, *, source_dir, overrides, source_name=None) -> StagedRun:
    """Write the run directory ``run_dir`` (created; it must not exist or be
    empty) for the parsed input TOML ``doc`` of the directory ``source_dir``.

    ``overrides``: the per-directory values, ``[mode.param]`` keys (e.g.
    ``SubShape``, or ``T``; an existing key is replaced in its
    own spelling) and ``"green"`` for ``[file.output] green``.

    The rules: ``[file.input] path_to_input = ""`` and ``initial =
    "initial.npz"`` (the file in ``run_dir``); ``initial_uhf`` and
    ``initial_mode`` removed; every file of ``[file.input]``
    (``geometry_uhf``, ``onebodyg_uhf``) and of ``[file.input.interaction]``
    located as the solver locates it from ``source_dir`` and copied (bytes,
    never a link) to ``inputs/<name>``, names made unique with ``-2``,
    ``-3`` in key order; ``[file.input] geometry_uhf``/``onebodyg_uhf`` =
    ``"inputs/<name>"``, ``[file.input.interaction] path_to_input =
    "inputs"`` and each entry ``<name>``. ``[file.output] path_to_output``
    and the output names must be relative and below the run directory;
    ``green = ""`` and ``energy = ""`` are refused; an absent ``energy`` is
    set to ``"energy.dat"``. The written TOML is re-parsed and checked by
    the solver's input reader. ``source_name`` (the original file name)
    goes into the header comment. ``initial = "initial.npz"`` is always
    written; the caller writes that file (for the preliminary paramagnetic
    run, its starting npz). Output names that would overwrite
    ``input.toml``, ``initial.npz``, ``hwave.log`` or ``inputs/`` are
    refused.

    Raises ``ValueError`` naming the key on any refusal; nothing is written
    outside ``run_dir``."""
    import tomli
    from hwave.solver import initial_facade
    from hwave.initial.writer import solver_initial_path
    run_dir = pathlib.Path(os.path.abspath(os.fspath(run_dir)))
    new, copies = _plan_run_doc(doc, source_dir=os.path.abspath(os.fspath(source_dir)),
                                overrides=overrides)
    text = dumps_toml(new)
    if tomli.loads(text) != new:
        raise ValueError("the rewritten input does not read back as written")
    if run_dir.exists() and any(run_dir.iterdir()):
        raise ValueError("{}: the run directory exists and is not empty".format(run_dir))
    run_dir.mkdir(parents=True, exist_ok=True)
    inputs = run_dir / RUN_INPUTS
    staged, done = [], {}
    if copies:
        inputs.mkdir()
    for table, key, source, name in copies:
        dest = inputs / name
        if name not in done:
            shutil.copyfile(source, str(dest))
            done[name] = (_sha256(dest), dest.stat().st_size)
        sha, nbytes = done[name]
        staged.append(StagedFile(table=table, key=key, source=source,
                                 name=RUN_INPUTS + "/" + name, sha256=sha, nbytes=nbytes))
    toml_path = run_dir / RUN_TOML
    with open(toml_path, "w", encoding="utf-8") as f:
        f.write(_HEADER.format(source_name or "the original input"))
        f.write(text)
    initial_facade.load_problem(toml_path)
    initial = run_dir / RUN_INITIAL
    target = solver_initial_path(toml_path)
    if target is None or os.path.realpath(target) != os.path.realpath(initial):
        raise AssertionError("{}: the solver would read the initial state from {}, not {}"
                             .format(toml_path, target, initial))
    for key, path in initial_facade.referenced_input_files(toml_path).items():
        if not path.is_file() or not _inside(path, run_dir):
            raise AssertionError("{}: {} resolves to {}, which is not a staged file in the "
                                 "run directory".format(toml_path, key, path))
    return StagedRun(toml_path=toml_path, initial_path=initial, staged=tuple(staged))


# --------------------------------------------------------------------------
# prepare
# --------------------------------------------------------------------------

def _plain(v):
    """``v`` as a plain JSON value (numpy scalars and read-only mappings
    included)."""
    if v is None or isinstance(v, (bool, str)):
        return v
    if isinstance(v, numbers.Integral):
        return int(v)
    if isinstance(v, numbers.Real):
        f = float(v)
        return f if math.isfinite(f) else None
    if isinstance(v, Mapping):
        return {str(k): _plain(x) for k, x in v.items()}
    if isinstance(v, (list, tuple, np.ndarray)):
        return [_plain(x) for x in v]
    raise TypeError("cannot record {!r} as JSON".format(v))


def _read_json(path, what):
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except (OSError, ValueError) as exc:
        raise ValueError("{} {}: {}".format(what, path, exc)) from None


def _extra_candidates(path):
    from hwave.initial.candidates import Candidate
    from hwave.initial.patterns import PatternSpec
    data = _read_json(path, "--extra-candidates")
    if not isinstance(data, list):
        raise ValueError("--extra-candidates {}: a JSON list of objects {{family, pattern, "
                         "label, citation}} expected".format(path))
    taken, out = set(), []
    allowed = ("family", "pattern", "label", "citation", "amplitude_requested")
    for i, obj in enumerate(data):
        where = "--extra-candidates {} [{}]".format(path, i)
        if not isinstance(obj, dict):
            raise ValueError("{}: an object expected".format(where))
        unknown = sorted(set(obj) - set(allowed))
        if unknown:
            raise ValueError("{}: unknown key(s) {}; the keys are {}".format(
                where, ", ".join(unknown), ", ".join(allowed)))
        family = obj.get("family")
        try:
            pattern = PatternSpec.from_json(obj.get("pattern"))
            kw = {"pattern": pattern, "label": obj.get("label"), "taken": taken}
            if "amplitude_requested" in obj:
                kw["amplitude_requested"] = obj["amplitude_requested"]
            if family == "user":
                if obj.get("citation") is not None:
                    raise ValueError("citation: only literature candidates carry one")
                out.append(Candidate.user(**kw))
            elif family == "literature":
                out.append(Candidate.literature(citation=obj.get("citation"), **kw))
            else:
                raise ValueError("family must be 'user' or 'literature', got {!r}".format(
                    family))
        except ValueError as exc:
            raise ValueError("{}: {}".format(where, exc)) from None
    return out


def _orbital_site_file(path):
    data = _read_json(path, "--orbital-site")
    if not isinstance(data, list) or not all(
            isinstance(x, numbers.Integral) and not isinstance(x, bool) for x in data):
        raise ValueError("--orbital-site {}: a JSON list of integers (one site label per "
                         "physical orbital) expected".format(path))
    return [int(x) for x in data]


def _realpath(p):
    return None if p is None else os.path.realpath(os.fspath(p))


def _input_record(toml, input_sha256, files):
    return {"toml": str(toml), "toml_sha256": _sha256(toml),
            "input_sha256": dict(input_sha256),
            "files_sha256": {str(k): (_sha256(k) if os.path.isfile(k) else None)
                             for k in files}}


def _differences(old, new):
    out = []
    oc, nc = old.get("command", {}), new["command"]
    for k in sorted(set(oc) | set(nc)):
        if oc.get(k) != nc.get(k):
            out.append("command {}: {} (WORK) vs {} (now)".format(
                k, json.dumps(oc.get(k)), json.dumps(nc.get(k))))
    oi, ni = old.get("input", {}), new["input"]
    if oi.get("toml_sha256") != ni["toml_sha256"]:
        out.append("the input TOML {} changed".format(ni["toml"]))
    for group in ("input_sha256", "files_sha256"):
        a, b = oi.get(group, {}) or {}, ni[group]
        for k in sorted(set(a) | set(b)):
            if a.get(k) != b.get(k):
                out.append("input file {} changed".format(k))
    return out


def _seed_fields(res, cap_ok=True):
    from hwave.initial.report import SEED_FIELDS
    out = {}
    for k in SEED_FIELDS:
        if k == "amplitude_requested":
            continue
        out[k] = _plain(getattr(res, k))
    if not cap_ok:
        out["cap_reason"] = None
    return out


def _rho(G, problem):
    g = np.asarray(G).reshape(problem.nvol, problem.nd, problem.nd)
    return np.swapaxes(g, -1, -2)


#: the para work directory inside WORK (the small tier)
PARA_DIR = "para_hf"


def _para_green(facts, cache, para=None):
    """The reference of the labels on the cell of ``facts``: the h = 0 seed
    on the base the seeds use, the para Hartree-Fock Hamiltonian folded to
    the cell when ``para`` is a converged para reference, else the bare
    ``H0`` paramagnet."""
    from hwave.initial.hints import fold_para, t_hint_bare
    from hwave.initial.seed import seed_pattern
    on_para = para is not None and para.status == "converged"
    key = (tuple(facts.problem.subshape), on_para)
    if key not in cache:
        p = facts.problem
        if on_para:
            base, seed_base = fold_para(para.H_para_k, facts, p.subshape), "para_hf"
        else:
            base, seed_base = p.ham_trans_k, "bare"
        cache[key] = seed_pattern(facts, base, None, T_run=p.T, W0=t_hint_bare(p).W0,
                                  seed_base=seed_base).G
    return cache[key]


def _state_label(G, facts, cache, *, q=None, phase=None, sign=None, spiral=False,
                 ordering=None, para=None):
    """The label of the state ``G`` from ``delta rho`` against the
    reference of its cell (:func:`_para_green`: the seed base at h = 0;
    the rho = G^T convention)."""
    from hwave.initial import labels
    from hwave.initial.symmetry import hermitian_green
    from hwave.initial.writer import HERMITE_TOL
    p = facts.problem
    g5 = np.asarray(G, dtype=np.complex128).reshape(p.nvol, p.ns, p.norb, p.ns, p.norb)
    gh = hermitian_green(g5, p.shape)
    off = float(np.max(np.abs(g5 - gh), initial=0.0))
    if off > HERMITE_TOL:
        raise ValueError("the Green function is not Hermitian (max |G(r)_ab - "
                         "conj(G(-r)_ba)| / 2 = {:.3g})".format(off))
    # Hermitian to round-off: the Hermitian parts, so that a delta rho of
    # round-off size is not judged against its own norm
    ref = hermitian_green(np.asarray(_para_green(facts, cache, para)).reshape(g5.shape),
                          p.shape)
    d = _rho(gh, p) - _rho(ref, p)
    if float(np.max(np.abs(d), initial=0.0)) <= labels.PARA_TOL:
        d = np.zeros_like(d)
    return labels.label(d, facts, q=q, eta=None, phase=phase, sign=sign,
                        spiral=spiral, ordering=ordering).name


#: q vectors whose weight is at least this fraction of the largest one are
#: named together (a near tie)
Q_TIE = 0.9
#: below this squared amplitude a state has no ordering vector
_Q_FLOOR = 1e-16


def _dominant_q(op):
    """The ordering vectors of a converged state, from the order parameters
    ``op`` of :mod:`hwave.results` (the JSON of :func:`report.order_parameter_records`,
    per orbital): the weight of q is ``|S(q)|^2 + |dn(q)|^2`` (the q = 0
    density, fixed by N, left out); the q of the largest weight and every q
    within :data:`Q_TIE` of it, as ``"a,b,c"`` strings in the reader's
    order. ``[]`` when nothing is ordered."""
    qs = op.get("q") or []
    spin = op.get("spin_amplitude") or [0.0] * len(qs)
    charge = op.get("charge_amplitude") or [0.0] * len(qs)
    w = []
    for q, a, b in zip(qs, spin, charge):
        zero = all(x.strip() == "0" for x in q.split(","))
        w.append(float(a) ** 2 + (0.0 if zero else float(b) ** 2))
    top = max(w, default=0.0)
    if top <= _Q_FLOOR:
        return []
    return [q for q, x in zip(qs, w) if x >= Q_TIE * top]


def _command_lines(run_ids):
    return ['(cd "{}" && python -m hwave.initial check --require-provenance '
            '{} && hwave {} > {} 2>&1)'.format(rid, RUN_TOML, RUN_TOML, RUN_LOG)
            for rid in run_ids]


def _commands_text(run_ids, *, stage0=False):
    name = STAGE0_COMMANDS if stage0 else COMMANDS
    if stage0:
        what = ["# The preliminary paramagnetic run written by python -m hwave.initial "
                "prepare --stage0.",
                "# bash {} works from any directory. For GNU parallel and job".format(name)]
        done = ("# Each line checks the start, then runs hwave; then {}, and the "
                "candidate runs are written.".format(PENDING_HINT))
    else:
        what = ["# Candidate runs written by python -m hwave.initial prepare.",
                "# bash {} works from any directory. For GNU parallel and job".format(name)]
        done = ("# Each line checks the initial state, then runs hwave; collect WORK "
                "reads the results.")
    lines = (['cd "$(dirname "$0")" || exit 1'] + what +
             ["# arrays, which run single lines, run from WORK: cd WORK && parallel < "
              + name, done] + _command_lines(run_ids))
    return "\n".join(lines) + "\n"


def _refuse_link(path):
    """``ValueError`` when ``path`` is a symbolic link (dangling or not)."""
    if os.path.islink(os.fspath(path)):
        raise ValueError("{}: is a symbolic link; hwave.initial writes its outputs "
                         "only as regular files inside WORK. Remove the link and "
                         "repeat the command".format(path))


def _write_text(path, text):
    """Write ``text`` to ``path`` atomically: into a fresh temporary file
    next to it (created exclusively, so an existing file or symbolic link of
    that name is never opened), then renamed onto ``path``. A symbolic link
    at ``path`` is refused (``ValueError``): the outputs of WORK are written
    in WORK only, never through a link to a file elsewhere."""
    path = pathlib.Path(os.fspath(path))
    _refuse_link(path)
    while True:
        tmp = path.parent / ".{}.{}.tmp".format(path.name, secrets.token_hex(8))
        try:
            fd = os.open(str(tmp), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o666)
            break
        except FileExistsError:
            continue
    try:
        try:
            f = os.fdopen(fd, "w", encoding="utf-8")
        except BaseException:
            os.close(fd)
            raise
        with f:
            f.write(text)
        os.replace(str(tmp), str(path))
    except BaseException:
        if os.path.lexists(str(tmp)):
            os.remove(str(tmp))
        raise


def _write_json(path, obj):
    _write_text(path, json.dumps(obj, indent=1, allow_nan=False) + "\n")


def _para_of(work, cset, base):
    """``(para, manifest entry)`` of the para step of ``cset`` kept in
    ``WORK/para_hf`` (the small tier): the para reference rebuilt from its
    files (:func:`hwave.initial.hints.para_from_work`), so the seeds,
    the labels and ``collect`` share one reference; ``(None, None)`` when
    no para Green function was written (the large tier, or a skipped or
    failed para step: the seeds use the bare ``H0``)."""
    from hwave.initial.hints import PARA_GREEN, PARA_RECORD, para_from_work
    root = work / PARA_DIR
    if cset.tier != "small" or not (root / PARA_RECORD).is_file():
        return None, None
    para = para_from_work(base, root)
    entry = {"dir": PARA_DIR, "record": PARA_DIR + "/" + PARA_RECORD,
             "record_sha256": _sha256(root / PARA_RECORD), "green": None,
             "green_sha256": None, "status": None if para is None else para.status,
             "T": None if para is None else float(para.T)}
    if para is not None:
        entry.update(green=PARA_DIR + "/" + PARA_GREEN,
                     green_sha256=_sha256(root / PARA_GREEN))
    return para, entry


def _check_para_entry(manifest, work):
    """The para entry of ``manifest`` (``None`` without a para Green
    function), its files checked: paths inside WORK, SHA-256 as prepare
    wrote them (``ValueError`` naming the file otherwise)."""
    entry = manifest.get("para")
    if not entry or not entry.get("green"):
        return None
    for key in ("dir", "record", "green"):
        name = entry.get(key)
        if not isinstance(name, str) or os.path.isabs(name) or ".." in _parts(name):
            raise ValueError("{}: para.{} {!r} is not a path inside WORK".format(
                work / MANIFEST, key, name))
    for key in ("record", "green"):
        path = work / entry[key]
        if not path.is_file() or _sha256(path) != entry.get(key + "_sha256"):
            raise ValueError("{}: the para reference file {} is missing or differs from "
                             "the file prepare wrote; the labels of WORK refer to it"
                             .format(work, entry[key]))
    return entry


def _uses_para(cset):
    """Whether the rows of ``cset`` were seeded on a para reference, or its
    hints were computed from one (the small tier)."""
    hint = (cset.model or {}).get("hint") or {}
    return hint.get("status") == "computed" or any(
        r.get("seed_base") == "para_hf" for r in cset.rows)


def _collect_para(cache, facts):
    """The para reference of the WORK being collected (``cache["para_dir"]``
    set by :func:`collect`), rebuilt once."""
    from hwave.initial.hints import para_from_work
    root = cache.get("para_dir")
    if root is None:
        return None
    if "para" not in cache:
        cache["para"] = para_from_work(facts, root)
    return cache["para"]


# --------------------------------------------------------------------------
# the preliminary paramagnetic run (prepare --stage0)
# --------------------------------------------------------------------------

def _stage0_overrides(doc, toml, thint):
    """The per-directory values of a preliminary paramagnetic run:
    ``SubShape = [1, 1, 1]``, ``T = T_hint``, ``EPS`` of the para step and
    ``[file.output] green`` (the input's name, else :data:`STAGE0_GREEN`)."""
    from hwave.solver import initial_facade
    out = (doc.get("file") or {}).get("output") if isinstance(doc.get("file"), dict) else None
    green = out.get("green") if isinstance(out, dict) else None
    if not isinstance(green, str) or not green:
        green = STAGE0_GREEN
    return {"SubShape": [1, 1, 1], "T": float(thint.T_hint),
            "EPS": initial_facade.para_eps(toml), "green": green}


def _check_stage0_tier(base, tier, lim, subshapes):
    """Refuse ``--stage0`` unless the small tier is chosen (checked before
    anything is written), and refuse ``SubShape`` lists the enumeration
    would refuse."""
    from hwave.initial import cells, enumerate as en
    from hwave.initial.facts import NO_SUBSHAPE_WARNING
    from hwave.initial.symmetry import su2_check
    if tier == "large":
        raise ValueError("--stage0 writes the preliminary paramagnetic run of the small "
                         "tier; the large tier has none. Give --tier small or auto")
    choice = en.choose_tier(base, "auto" if tier is None else tier,
                            su2=su2_check(base.problem),
                            limits=en.Limits.from_keywords(**lim))
    if choice.name != "small":
        raise ValueError("--stage0 writes the preliminary paramagnetic run of the small "
                         "tier, but the tier chosen for this input is {} ({}); give "
                         "--tier small, or prepare without --stage0".format(
                             choice.name, "; ".join(choice.warnings) or choice.chosen_by))
    if subshapes is None:
        if NO_SUBSHAPE_WARNING in base.warnings:
            raise ValueError(NO_SUBSHAPE_WARNING)
    else:
        cells.check_subshapes(en._subshape_list(subshapes), base.cellshape)


def _stage0_entry(name, start, cid, staged, rd):
    return {"run_dir": name, "start": start, "candidate_id": cid,
            "input_toml": RUN_TOML, "input_toml_sha256": _sha256(staged.toml_path),
            "initial": RUN_INITIAL, "initial_sha256": _sha256(staged.initial_path),
            "input_sha256": {s.key: s.sha256 for s in staged.staged},
            "staged": [{"table": s.table, "key": s.key, "name": s.name,
                        "sha256": s.sha256, "bytes": s.nbytes} for s in staged.staged],
            "outcome": None, "reason": None}


def _print_stage0(work, names, stdout):
    print("{}: the preliminary paramagnetic run is written to {}; it runs outside prepare."
          .format(work, " and ".join(names)), file=stdout)
    print("Run: bash {} (or, from WORK, cd {} && parallel < {}):".format(
        work / STAGE0_COMMANDS, work, STAGE0_COMMANDS), file=stdout)
    for line in _command_lines(names):
        print("  " + line, file=stdout)
    print("Then {}: python -m hwave.initial prepare with the same arguments writes the "
          "candidate runs (exit code {} until then).".format(PENDING_HINT,
                                                            EXIT_STAGE0_PENDING),
          file=stdout)


def _write_stage0(work, doc, toml, base, thint, *, para_check, record, stdout):
    """Write the preliminary paramagnetic run into the new or empty WORK
    ``work`` and return :data:`EXIT_STAGE0_PENDING`."""
    import hwave
    from hwave.initial.hints import para_starts
    from hwave.initial.writer import check_initial_state, write_start_state
    overrides = _stage0_overrides(doc, toml, thint)
    # the staging rules, before anything is written
    _plan_run_doc(doc, source_dir=str(toml.parent), overrides=overrides)
    starts = para_starts(base, thint)
    runs = STAGE0_RUNS[:2 if para_check else 1]
    created = not work.exists()
    work.mkdir(parents=True, exist_ok=True)
    made = []
    try:
        entries, staged_bytes = [], 0
        for name, start, cid in runs:
            rd = work / name
            made.append(rd)
            staged = stage_run_dir(rd, doc, source_dir=toml.parent, overrides=overrides,
                                   source_name=toml.name)
            write_start_state(staged.toml_path, starts[start], out=staged.initial_path,
                              candidate_id=cid)
            if not check_initial_state(staged.toml_path, require_provenance=True).verified:
                raise AssertionError("{}: the written start is not verified".format(rd))
            staged_bytes += sum(s.nbytes for s in {s.name: s for s in staged.staged}.values())
            entries.append(_stage0_entry(name, start, cid, staged, rd))
        _write_text(work / STAGE0_COMMANDS, _commands_text([r[0] for r in runs], stage0=True))
        made.append(work / STAGE0_COMMANDS)
        manifest = dict(record, manifest_version=MANIFEST_VERSION,
                        hwave_version=hwave.__version__, work_dir_at_prepare=str(work),
                        stage=STAGE_PENDING,
                        stage0={"T": float(thint.T_hint), "eps": overrides["EPS"],
                                "selected": None, "runs": entries},
                        staged_bytes=staged_bytes, runs=[], para=None)
        _write_json(work / MANIFEST, manifest)
    except BaseException:
        if created:
            shutil.rmtree(work, ignore_errors=True)
        else:
            for p in made:
                if p.is_dir():
                    shutil.rmtree(p, ignore_errors=True)
                elif p.exists():
                    p.unlink()
        raise
    _print_stage0(work, [r[0] for r in runs], stdout)
    return EXIT_STAGE0_PENDING


def _inside_work(name, where):
    if not isinstance(name, str) or not name or os.path.isabs(name) or ".." in _parts(name):
        raise ValueError("{}: {!r} is not a path inside WORK".format(where, name))


def _stage0_runs(manifest, work):
    """The run entries of the preliminary paramagnetic run in
    ``manifest``, their paths checked."""
    st = manifest.get("stage0")
    runs = st.get("runs") if isinstance(st, dict) else None
    if not isinstance(runs, list) or not runs:
        raise ValueError("{}: stage {} without the entries of the preliminary paramagnetic "
                         "run".format(work / MANIFEST, STAGE_PENDING))
    for i, e in enumerate(runs):
        where = "{}: stage0.runs[{}]".format(work / MANIFEST, i)
        if not isinstance(e, dict) or not isinstance(e.get("staged"), list):
            raise ValueError("{}: an entry of prepare --stage0 expected".format(where))
        for key in ("run_dir", "input_toml", "initial"):
            _inside_work(e.get(key), where)
        for s in e["staged"]:
            _inside_work(s.get("name") if isinstance(s, dict) else None, where)
    return runs


def _verify_stage0(work, e):
    """``ValueError`` naming the file when an input of the preliminary
    paramagnetic run ``e`` is missing or differs from what prepare wrote."""
    files = [(e["input_toml"], e.get("input_toml_sha256")),
             (e["initial"], e.get("initial_sha256"))]
    files += [(s["name"], s.get("sha256")) for s in e["staged"]]
    for rel, sha in files:
        path = work / e["run_dir"] / rel
        shown = os.path.join(e["run_dir"], rel)
        if not path.is_file():
            raise ValueError("{}: {} is missing; the preliminary paramagnetic run reads the "
                             "files prepare wrote. Repeat prepare into a new directory"
                             .format(work, shown))
        if _sha256(path) != sha:
            raise ValueError("{}: {} differs from the file prepare wrote (SHA-256 in {}); "
                             "the preliminary paramagnetic run must run on the files "
                             "prepare wrote. Repeat prepare into a new directory".format(
                                 work, shown, MANIFEST))


@dataclass(frozen=True)
class _Stage0Outcome:
    """The state of one preliminary paramagnetic run: ``state`` is
    ``"pending"`` (not run yet), ``"failed"`` (ran and failed: no usable
    Green function), ``"not_converged"`` (the SCF did not converge) or
    ``"converged"``; ``para`` its para reference (``None`` when pending),
    ``green`` its output green file, ``reason`` why it failed."""
    state: str
    para: object = None
    green: Optional[pathlib.Path] = None
    reason: Optional[str] = None


def _stage0_outcome(work, e, base, thint):
    """Classify the preliminary paramagnetic run ``e`` of WORK: a valid
    status file decides; without one, ``hwave.log`` (none and no output:
    not run yet)."""
    from hwave.initial import hints
    from hwave.initial.writer import check_green
    from hwave.results import (ResultsError, UnsupportedLayoutError, open_run, read_energy,
                               read_scf_status, read_uhfk_green)
    rd = work / e["run_dir"]
    toml = rd / e["input_toml"]

    def failed(reason):
        msg = _relative_text("the preliminary paramagnetic run in {} failed: {}".format(
            e["run_dir"], reason), rd, work)
        return _Stage0Outcome("failed", hints.para_failed(base, thint, msg), None, msg)

    try:
        run = open_run(toml, base_dir=rd)
    except ResultsError as exc:
        return failed(str(exc))
    eps = run.param("EPS")
    iteration_max = run.param("IterationMax")
    eps_tol = 10.0 ** (-float(eps))
    iterations = residual = None
    status_path = _status_file(run)
    if status_path is not None:
        try:
            st = read_scf_status(run)
        except UnsupportedLayoutError as exc:
            return failed("{} unreadable: {}".format(status_path.name, exc))
        iterations, residual = st.iterations, _plain(st.residual)
        if st.iterations == 0:
            return failed(st.note or "no SCF step performed")
        converged = bool(st.converged)
    else:
        log = rd / RUN_LOG
        if not log.is_file():
            if not _any_output(run):
                return _Stage0Outcome("pending")
            return failed("its outputs exist but {} does not, so the convergence cannot be "
                          "verified".format(RUN_LOG))
        verdict, residual, summary = _log_verdict(log, eps_tol)
        if verdict != "failed" and not _log_ended(log):
            # no status file and no final line of the solver: still running,
            # or stopped from outside; not decided here
            return _Stage0Outcome("pending")
        if verdict == "failed" or summary is not None:
            return failed(summary or "the run did not converge")
        converged = verdict == "converged"
    try:
        g = read_uhfk_green(run)
        if g.green is None:
            raise UnsupportedLayoutError("{}: {}".format(
                g.source, g.meta.get("green_unavailable", "no green")))
        free = read_energy(run).free_total
        if free is None:
            raise UnsupportedLayoutError("no FreeEnergy_Total in the energy output")
        G = np.asarray(g.green)
        try:
            check_green(G, base.problem, what="{}: green".format(g.source))
        except ValueError as exc:
            return failed("the output green is not a state of this input ({}); it is not "
                          "used as the para reference".format(exc))
        para = hints.para_from_green(base, thint, G, free_energy=free,
                                     converged=converged, iterations=iterations,
                                     scf_residual=residual, eps=eps,
                                     iteration_max=iteration_max)
    except (ResultsError, OSError, ValueError) as exc:
        return failed("the output green or energy cannot be read: {}".format(exc))
    return _Stage0Outcome("converged" if converged else "not_converged", para,
                          pathlib.Path(g.source))


def _stage0_reentry(work, manifest, base, thint):
    """The para reference from the preliminary paramagnetic run of WORK
    (``(para, selected entry or None, outcomes)``), or ``None`` while a run
    has not run yet. The inputs are checked against the manifest first."""
    from hwave.initial.hints import choose_para
    runs = _stage0_runs(manifest, work)
    T = manifest["stage0"].get("T")
    if not isinstance(T, numbers.Real) or abs(float(T) - float(thint.T_hint)) > \
            1e-12 * max(1.0, abs(float(thint.T_hint))):
        raise ValueError("{}: the preliminary paramagnetic run was written at T_hint = {}, "
                         "this input gives {}".format(work / MANIFEST, T, thint.T_hint))
    for e in runs:
        _verify_stage0(work, e)
    outcomes = [_stage0_outcome(work, e, base, thint) for e in runs]
    if any(o.state == "pending" for o in outcomes):
        return None
    if len(outcomes) == 1:
        para, idx = outcomes[0].para, 0
    else:
        para, sel = choose_para(outcomes[0].para, outcomes[1].para)
        idx = 0 if sel == "main" else 1
    return para, (runs[idx], outcomes[idx]), outcomes


def _keep_stage0_para(work, para, selected):
    """Keep the para reference ``para`` of the preliminary paramagnetic run
    in ``WORK/para_hf`` as the in-process step does (its green file copied
    when it has one, and the record); returns the directory."""
    from hwave.initial.hints import PARA_GREEN, write_para_record
    root = work / PARA_DIR
    # Not atomic on purpose: a re-entry that dies between the mkdir below
    # and the manifest update leaves WORK/para_hf behind, and every later
    # re-entry is refused here (prepare's rollback leaves a WORK/para_hf it
    # did not make, and the manifest, as they are). That fails safe (a
    # refusal, never a mixed reference), so the directory is not written
    # through a temporary name.
    if root.exists():
        raise ValueError("{}: exists before the preliminary paramagnetic run was read; "
                         "repeat prepare into a new directory".format(root))
    root.mkdir()
    entry, outcome = selected
    start = None
    if para.G_para is not None and outcome.green is not None:
        shutil.copyfile(str(outcome.green), str(root / PARA_GREEN))
        start = entry["start"]
    write_para_record(para, start, str(root))
    return root


def prepare(input_path, out, *, subshapes=None, tier=None, extra_candidates=None,
            t_hint=None, amplitude=None, max_candidates=None, orbital_site=None,
            para_check=False, stage0=False, stdout=None, **limits) -> int:
    """``python -m hwave.initial prepare`` (returns the exit code).

    Enumerates the candidates of the UHFk input ``input_path`` (a TOML, or
    a directory holding exactly one), and writes WORK = ``out``:
    ``coverage.json`` (``pre_run``), ``manifest.json``, ``candidates.json``,
    ``commands.txt`` and one run directory per run (:func:`stage_run_dir`,
    then ``initial.npz``, checked with ``check_initial_state``).

    In the small tier the paramagnetic Hartree-Fock step runs into
    ``WORK/para_hf`` (``green.npz``, ``para_hf.json`` and its runs; their
    SHA-256 recorded in ``manifest.json`` under ``para``): the seeds are
    built on its Hamiltonian, and the labels of prepare and collect take
    the h = 0 seed on that base as their reference. ``para_check``: the
    para step also runs from the uniform orbital occupation (small tier;
    refused with the large tier given explicitly).

    ``stage0`` (small tier; refused when the large tier is given or
    chosen): the para step runs outside ``prepare``, as the preliminary
    paramagnetic run (module text): the first call writes
    ``WORK/stage0_para/`` (and ``stage0_para_check/`` with ``para_check``),
    ``stage0_commands.txt`` and ``manifest.json`` (``stage =
    "para_pending"``, the SHA-256 of their inputs) and returns
    :data:`EXIT_STAGE0_PENDING`. The same command on that WORK checks those
    inputs against the manifest (a difference is refused, naming the file)
    and classifies each run: not run yet (no ``hwave.log``, status file or
    output) returns :data:`EXIT_STAGE0_PENDING` again and changes nothing;
    otherwise the valid status file, else the log, decides. A converged run
    becomes the para reference in ``WORK/para_hf`` and the candidate runs
    are written; a run that ran and failed, or did not converge, gives the
    bare-``H0`` seeds with the hints skipped for that reason. With
    ``para_check``, when exactly one of the two runs converged it is kept,
    with a warning. ``stage0`` is part of the command: a WORK prepared
    with it re-enters only with it. Without an interacting orbital the
    para step is skipped at once, so ``stage0`` changes nothing.

    ``extra_candidates`` and ``orbital_site`` are JSON file paths;
    ``**limits`` are further :class:`~hwave.initial.enumerate.Limits`
    fields. An existing non-empty WORK written by the same command for the
    same inputs is left as it is; any other non-empty directory is refused,
    naming the difference. Raises ``ValueError`` (exit code 1 in
    :func:`main`) on a refusal."""
    from hwave.initial.enumerate import candidate_states
    from hwave.initial.facts import _resolve_toml, model_facts
    from hwave.initial.report import coverage_document
    from hwave.initial.writer import InitialStateError, check_initial_state, seed_and_write
    from hwave.initial.candidates import _thaw_json
    import hwave
    import tomli

    stdout = stdout or sys.stdout
    if tier == "medium":
        raise NotImplementedError("tier 'medium' arrives in a later release; this release "
                                  "has the large and small tiers")
    if para_check and tier == "large":
        raise ValueError("--para-check checks the para Hartree-Fock step of the small tier; "
                         "the large tier has none")
    if stage0 and tier == "large":
        raise ValueError("--stage0 writes the preliminary paramagnetic run of the small "
                         "tier; the large tier has none. Give --tier small or auto")
    toml = pathlib.Path(os.path.realpath(_resolve_toml(input_path)))
    work = pathlib.Path(os.path.abspath(os.fspath(out)))
    lim = dict(limits)
    if amplitude is not None:
        lim["amplitude"] = amplitude
    if max_candidates is not None:
        lim["max_candidates"] = max_candidates
    extra_path = _realpath(extra_candidates)
    site_path = _realpath(orbital_site)
    command = {
        "input": str(toml),
        "subshape": None if subshapes is None else [list(map(int, b)) for b in subshapes],
        "tier": tier, "extra_candidates": extra_path, "t_hint": t_hint,
        "orbital_site": site_path, "para_check": bool(para_check),
        "stage0": bool(stage0), "limits": {k: lim[k] for k in sorted(lim)},
    }
    with open(toml, "rb") as f:
        try:
            doc = tomli.load(f)
        except tomli.TOMLDecodeError as exc:
            raise ValueError("{}: not a valid TOML file: {}".format(toml, exc)) from None
    site = _orbital_site_file(site_path) if site_path else None
    # the input is read by the solver's rules first (case variants refused)
    base = model_facts(toml, subshape=(1, 1, 1), orbital_site=site)
    files = [p for p in (extra_path, site_path) if p is not None]
    record = {"command": command, "input": _input_record(toml, base.input_sha256, files)}
    from hwave.initial.hints import t_hint_bare
    # the hint temperature of the preliminary paramagnetic run (only then)
    thint = t_hint_bare(base.problem, t_user=t_hint) if stage0 else None
    reentry, stage0_record = None, None

    if work.exists() and not work.is_dir():
        raise ValueError("{}: --out names a file, not a directory".format(work))
    if work.is_dir() and any(work.iterdir()):
        mpath = work / MANIFEST
        if not mpath.is_file():
            raise ValueError("{}: not empty and holds no {} of hwave.initial; give an empty "
                             "or new directory".format(work, MANIFEST))
        old = _read_json(mpath, "manifest")
        diff = _differences(old, record)
        if diff:
            raise ValueError("{}: written by another prepare command or for other inputs: "
                             "{}. Give a new directory, or repeat the command that wrote "
                             "it".format(work, "; ".join(diff)))
        if old.get("stage") == STAGE_PENDING:
            # the preliminary paramagnetic run: classified before anything else
            if thint is None:
                thint = t_hint_bare(base.problem, t_user=t_hint)
            reentry = _stage0_reentry(work, old, base, thint)
            if reentry is None:
                _print_stage0(work, [e["run_dir"] for e in old["stage0"]["runs"]], stdout)
                return EXIT_STAGE0_PENDING
            stage0_record = old["stage0"]
        else:
            print("{}: already prepared by this command for these inputs; nothing to do. "
                  "Run bash {} and then python -m hwave.initial collect {}".format(
                      work, work / COMMANDS, work), file=stdout)
            return EXIT_OK

    # the staging rules are checked before anything is enumerated or written
    _plan_run_doc(doc, source_dir=str(toml.parent), overrides={"SubShape": [1, 1, 1]})
    extras = _extra_candidates(extra_path) if extra_path else ()
    if stage0 and reentry is None and base.n_c > 0:
        _check_stage0_tier(base, tier, lim, subshapes)
        return _write_stage0(work, doc, toml, base, thint, para_check=para_check,
                             record=record, stdout=stdout)

    created = not work.exists()
    work.mkdir(parents=True, exist_ok=True)
    # the rollback removes only what this call made: a WORK/para_hf that
    # exists already (a re-entry interrupted after it was written, or a
    # concurrent prepare) is left as it is, so a refused re-entry changes
    # nothing and the next one is refused again
    made = [] if (work / PARA_DIR).exists() else [work / PARA_DIR]
    try:
        if reentry is not None:
            from hwave.initial.hints import para_from_work
            para0, selected, outcomes = reentry
            root = _keep_stage0_para(work, para0, selected)
            # the seeds and labels use the reference read back from
            # WORK/para_hf, as after the in-process step
            given = {"para": para_from_work(base, root) or para0}
            stage0_record = dict(stage0_record, selected=selected[0]["start"], runs=[
                dict(e, outcome=o.state, reason=o.reason)
                for e, o in zip(stage0_record["runs"], outcomes)])
        else:
            given = {"para_dir": work / PARA_DIR}
        cset = candidate_states(toml, subshapes=subshapes,
                                tier="auto" if tier is None else tier, t_hint=t_hint,
                                extra=extras, orbital_site=site, para_check=para_check,
                                **given, **lim)
        para, para_entry = _para_of(work, cset, base)
        if not (work / PARA_DIR).exists() and work / PARA_DIR in made:
            made.remove(work / PARA_DIR)
        rows = [_thaw_json(r) for r in cset.rows]
        index = {r["run_id"]: i for i, r in enumerate(rows)}
        kept, entries, staged_bytes, cache = [], [], 0, {}
        for run in cset.runs:
            rd = work / run.run_id
            made.append(rd)
            staged = stage_run_dir(rd, doc, source_dir=toml.parent,
                                   overrides={"SubShape": list(run.subshape)},
                                   source_name=toml.name)
            row = rows[index[run.run_id]]
            try:
                path, res, facts = seed_and_write(staged.toml_path, run,
                                                  out=staged.initial_path, para=para,
                                                  t_hint=t_hint)
            except InitialStateError as exc:
                # a writer refusal of a user or literature row becomes a
                # "refused:" row; for a generated row it is an error
                if run.candidate.family not in ("user", "literature"):
                    raise
                shutil.rmtree(rd)
                made.pop()
                row.update(tried=False, not_tried_reason="refused: {}".format(exc))
                continue
            if path is None:
                shutil.rmtree(rd)
                made.pop()
                row.update(_seed_fields(res), tried=False,
                           not_tried_reason=res.not_tried_reason)
                continue
            result = check_initial_state(staged.toml_path, require_provenance=True)
            if not result.verified:
                raise AssertionError("{}: the written seed is not verified".format(rd))
            spec = run.candidate.pattern
            row.update(_seed_fields(res))
            modulated = spec.modulation in ("cos", "spiral")
            row["label"] = _state_label(
                res.G, facts, cache, q=spec.q if modulated else None,
                phase=spec.phase if modulated else None,
                sign=spec.sign if modulated else None, spiral=spec.modulation == "spiral",
                para=para)
            kept.append(run)
            staged_bytes += sum(s.nbytes for s in {s.name: s for s in staged.staged}.values())
            entries.append({
                "run_id": run.run_id, "run_dir": run.run_id, "subshape": list(run.subshape),
                "input_toml": RUN_TOML, "input_toml_sha256": _sha256(staged.toml_path),
                "initial": RUN_INITIAL, "initial_sha256": _sha256(staged.initial_path),
                "input_sha256": {s.key: s.sha256 for s in staged.staged},
                "staged": [{"table": s.table, "key": s.key, "name": s.name,
                            "sha256": s.sha256, "bytes": s.nbytes} for s in staged.staged],
            })
        cset2 = dataclasses.replace(cset, runs=tuple(kept), rows=rows)
        chosen = "auto (interim rule)" if tier in (None, "auto") else "user"
        cov = coverage_document(cset2, phase_marker="pre_run", tier_chosen_by=chosen)
        _write_json(work / CANDIDATES, dict(cset2.to_json(), tier_chosen_by=chosen))
        made.append(work / CANDIDATES)
        _write_json(work / COVERAGE, cov)
        made.append(work / COVERAGE)
        _write_text(work / COMMANDS, _commands_text([e["run_id"] for e in entries]))
        made.append(work / COMMANDS)
        manifest = dict(record, manifest_version=MANIFEST_VERSION,
                        hwave_version=hwave.__version__, work_dir_at_prepare=str(work),
                        stage=STAGE_DONE, stage0=stage0_record,
                        staged_bytes=staged_bytes, runs=entries, para=para_entry)
        _write_json(work / MANIFEST, manifest)
    except BaseException:
        if created:
            shutil.rmtree(work, ignore_errors=True)
        else:
            for p in made:
                if p.is_dir():
                    shutil.rmtree(p, ignore_errors=True)
                elif p.exists():
                    p.unlink()
        raise
    n_not = sum(1 for r in rows if not r["tried"])
    print("{}: {} runs written ({} rows not tried, see {}); {} bytes of staged input "
          "copies.".format(work, len(entries), n_not, COVERAGE, staged_bytes), file=stdout)
    print("Run them with: bash {} (or, from WORK, cd {} && parallel < {}); then: "
          "python -m hwave.initial collect {}".format(work / COMMANDS, work, COMMANDS, work),
          file=stdout)
    return EXIT_OK


# --------------------------------------------------------------------------
# check
# --------------------------------------------------------------------------

def check(config, npz=None, *, require_provenance=False, stdout=None, stderr=None) -> int:
    """``python -m hwave.initial check`` (returns the exit code)."""
    from hwave.initial.writer import (InitialStateError, ProvenanceMissingError,
                                      check_initial_state)
    stdout, stderr = stdout or sys.stdout, stderr or sys.stderr
    try:
        res = check_initial_state(config, npz, require_provenance=require_provenance)
    except ProvenanceMissingError as exc:
        print("refused: {}".format(exc), file=stderr)
        return EXIT_NO_PROVENANCE
    except (InitialStateError, ValueError, OSError) as exc:
        print("refused: {}".format(exc), file=stderr)
        return EXIT_ERROR
    for w in res.warnings:
        print("warning: {}".format(w), file=stderr)
    print("{}: {} (member {})".format(config, "verified" if res.verified else
                                      "checked, not verified", res.key), file=stdout)
    return EXIT_OK


# --------------------------------------------------------------------------
# collect
# --------------------------------------------------------------------------

_SUCCEEDED = re.compile(r"calculation succeeded: rest=(\S+), eps=")
_FAILED = re.compile(r"calculation failed: rest=(\S+), eps=")
_ERROR_LINE = re.compile(r"^\S+ \S+ (ERROR|CRITICAL) (.*)$")


def _float(text):
    try:
        return float(text)
    except ValueError:
        return None


def _log_verdict(log_path, eps_tol):
    """``(status, residual, summary)`` of a run from its ``hwave.log``."""
    with open(log_path, encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()
    errors = [m.group(2) for m in map(_ERROR_LINE.match, lines) if m]
    tb = [i for i, ln in enumerate(lines) if ln.startswith("Traceback (most recent call")]
    if tb:
        tail = [ln for ln in lines[tb[-1]:] if ln.strip()]
        return "failed", None, "traceback in {}: {}".format(RUN_LOG, tail[-1].strip())
    if errors:
        return "failed", None, "solver error in {}: {}".format(RUN_LOG, errors[-1])
    for ln in reversed(lines):
        m = _SUCCEEDED.search(ln)
        if m:
            r = _float(m.group(1).rstrip(","))
            if r is None or (eps_tol is not None and not r <= eps_tol):
                return "not_converged", r, "the final residual exceeds the tolerance"
            return "converged", r, None
        m = _FAILED.search(ln)
        if m:
            return "not_converged", _float(m.group(1).rstrip(",")), None
    return "not_converged", None, "{} has no convergence message".format(RUN_LOG)


def _log_ended(log_path):
    """Whether ``hwave.log`` holds a final line of the solver (its
    convergence message, success or failure, or an error)."""
    with open(log_path, encoding="utf-8", errors="replace") as f:
        for ln in f:
            if _SUCCEEDED.search(ln) or _FAILED.search(ln) or _ERROR_LINE.match(ln) \
                    or ln.startswith("Traceback (most recent call"):
                return True
    return False


def _seed_problem(entry, run_dir, run_toml):
    """Why the run's seed or inputs are not the prepared ones (``None``
    when they are)."""
    from hwave.solver import initial_facade
    from hwave.initial.writer import check_initial_state, solver_initial_path
    try:
        target = solver_initial_path(run_toml)
    except Exception as exc:  # an unreadable input.toml
        return "{} cannot be read: {}".format(RUN_TOML, exc)
    initial = run_dir / entry["initial"]
    if target is None or not initial.is_file() or not os.path.exists(target) \
            or not os.path.samefile(target, initial):
        return "{} no longer names {} as [file.input] initial".format(RUN_TOML, entry["initial"])
    if _sha256(initial) != entry["initial_sha256"]:
        return "{} differs from the file prepare wrote".format(entry["initial"])
    try:
        files = initial_facade.referenced_input_files(run_toml)
    except ValueError as exc:
        return "{}: {}".format(RUN_TOML, exc)
    want = entry["input_sha256"]
    for key in sorted(set(files) | set(want)):
        path = files.get(key)
        if path is None or key not in want:
            return "the input file set changed ({})".format(key)
        if not path.is_file() or not _inside(path, run_dir):
            return "{} = {} is not a staged file of the run directory".format(key, path)
        if _sha256(path) != want[key]:
            return "the staged input {} ({}) changed after prepare".format(
                os.path.relpath(path, run_dir), key)
    try:
        check_initial_state(run_toml, require_provenance=True)
    except ValueError as exc:
        return "the initial-state check refused it: {}".format(exc)
    return None


def _status_file(run):
    path = run.outputs.get("scf_status")
    return path if path is not None and path.is_file() else None


def _any_output(run):
    for key, path in run.outputs.items():
        if key == "print_check":
            continue
        for p in (path, pathlib.Path(str(path) + ".npz")):
            if p.exists():
                return True
    return False


def _read_outputs(row, run, run_toml, model, facts_cache, para_cache):
    """Fill the energies, order parameters and final label of a run whose
    outputs must be readable (``ResultsError`` otherwise). ``run`` is the
    :class:`hwave.results.Run`, opened with the run directory as its base."""
    from hwave.results import UnsupportedLayoutError, read_energy, read_uhfk_green
    from hwave.initial.facts import model_facts
    from hwave.initial.report import energy_fields, order_parameter_records, orbital_centres
    en = read_energy(run)
    row["energy_total"] = _plain(en.total)
    row["free_energy_total"] = _plain(en.free_total)
    row["rest"] = _plain(en.values.get("Rest"))
    e = energy_fields(row, facts=model)
    for k in ("free_energy_per_cell", "free_energy_per_site", "free_energy_per_site_reason"):
        row[k] = e[k]
    key = str(run_toml)
    if key not in facts_cache:
        facts_cache[key] = model_facts(run_toml)
    facts = facts_cache[key]
    row["order_parameters"] = order_parameter_records(
        run, site_groups=model["site_groups"], orbital_centres=orbital_centres(facts))
    g = read_uhfk_green(run)
    G = g.sublattice.green if g.sublattice is not None else g.green
    if G is None:
        raise UnsupportedLayoutError("{}: the output green cannot be read".format(run_toml))
    qs = _dominant_q(row["order_parameters"]["per_orbital"])
    # the ordering vectors come from the on-site amplitudes, so they follow
    # the on-site channels of the name, not a bond channel
    row["final_label"] = _state_label(G, facts, para_cache,
                                      ordering=[q.split(",") for q in qs],
                                      para=_collect_para(para_cache, facts))


_OUTPUT_ROW_FIELDS = ("energy_total", "free_energy_total", "rest", "free_energy_per_cell",
                      "free_energy_per_site", "free_energy_per_site_reason",
                      "order_parameters", "final_label")


def _relative_text(text, run_dir, work):
    """``text`` with the run directory and WORK written relative (so the
    collected result does not depend on where WORK is)."""
    if not isinstance(text, str):
        return text
    roots = []
    for base, rel in ((run_dir, "."), (work, "WORK")):
        for p in {str(base), os.path.realpath(str(base))}:
            roots.append((p, rel))
    for p, rel in sorted(roots, key=lambda x: -len(x[0])):
        text = text.replace(p + os.sep, "" if rel == "." else rel + "/").replace(p, rel)
    return text


def _collect_row(row, entry, work, model, warnings, facts_cache, para_cache):
    n = len(warnings)
    _collect_row_paths(row, entry, work, model, warnings, facts_cache, para_cache)
    run_dir = work / entry["run_dir"]
    row["error_summary"] = _relative_text(row.get("error_summary"), run_dir, work)
    warnings[n:] = [_relative_text(w, run_dir, work) for w in warnings[n:]]


def _collect_row_paths(row, entry, work, model, warnings, facts_cache, para_cache):
    from hwave.results import ResultsError, UnsupportedLayoutError, open_run, read_scf_status
    run_dir = work / entry["run_dir"]
    run_toml = run_dir / entry["input_toml"]
    row["run_dir"] = entry["run_dir"]
    try:
        # every reader gets this Run: base directory = the run directory, so
        # no run record is consulted and a moved WORK reads the same
        run = open_run(run_toml, base_dir=run_dir)
        row["iteration_max"] = _plain(run.param("IterationMax"))
        row["eps"] = _plain(run.param("EPS"))
        row["mix"] = _plain(run.param("Mix"))
    except ResultsError as exc:
        row.update(status="seed_unverified", error_summary=str(exc))
        return
    problem = _seed_problem(entry, run_dir, run_toml)
    if problem is not None:
        row.update(status="seed_unverified", error_summary=problem)
        return
    eps = row["eps"]
    eps_tol = (10.0 ** (-float(eps)) if isinstance(eps, numbers.Real)
               and not isinstance(eps, bool) else None)
    try:
        status_path = _status_file(run)
        st = None
        if status_path is not None:
            try:
                st = read_scf_status(run)
            except UnsupportedLayoutError as exc:
                # every output reader refuses a run whose status file cannot be
                # read (the zero-step guard), so the run's results are not
                # usable: failed, never ranked
                name = status_path.name
                warnings.append("{}: {} unreadable: {}".format(entry["run_id"], name, exc))
                row.update(status="failed", failure_kind="unreadable_output",
                           error_summary="{} unreadable: {}".format(name, exc))
                return
        if st is not None:
            row["iterations"] = st.iterations
            row["scf_residual"] = _plain(st.residual)
            row["scf_eps_tol"] = _plain(st.eps)
            row["scf_note"] = st.note
            if st.iterations == 0:
                row.update(status="failed", failure_kind="zero_steps",
                           error_summary=st.note or "no SCF step performed")
                return
            if st.converged:
                row["status"] = "converged"
                _read_outputs(row, run, run_toml, model, facts_cache, para_cache)
            else:
                row["status"] = "not_converged"
            return
        log = run_dir / RUN_LOG
        has_out = _any_output(run)
        if not log.is_file() and not has_out:
            row["status"] = "not_run"
            return
        if not log.is_file():
            row["status"] = "convergence_unverified"
            _read_outputs(row, run, run_toml, model, facts_cache, para_cache)
            return
        status, residual, summary = _log_verdict(log, eps_tol)
        row["scf_residual"] = residual
        row["scf_eps_tol"] = eps_tol
        if status == "failed":
            row.update(status="failed", failure_kind="crash", error_summary=summary)
            return
        row["status"] = status
        if summary:
            row["error_summary"] = summary
        if status == "converged":
            _read_outputs(row, run, run_toml, model, facts_cache, para_cache)
    except (ResultsError, OSError, ValueError) as exc:
        for k in _OUTPUT_ROW_FIELDS:
            row[k] = None
        row.update(status="failed", failure_kind="unreadable_output",
                   error_summary="{}: {}".format(type(exc).__name__, exc))


_ENTRY_KEYS = {"run_id": str, "run_dir": str, "input_toml": str, "initial": str,
               "initial_sha256": str, "input_sha256": dict}


def _check_manifest(manifest, path):
    """The run entries of ``manifest`` by run id; ``ValueError`` naming the
    file and the key when its shape is not the one prepare writes."""
    if not isinstance(manifest, dict):
        raise ValueError("{}: a JSON object expected".format(path))
    if manifest.get("manifest_version") != MANIFEST_VERSION:
        raise ValueError("{}: manifest_version {!r} is not {}".format(
            path, manifest.get("manifest_version"), MANIFEST_VERSION))
    runs = manifest.get("runs")
    if not isinstance(runs, list):
        raise ValueError("{}: 'runs' must be a list of run entries".format(path))
    out = {}
    for i, e in enumerate(runs):
        if not isinstance(e, dict):
            raise ValueError("{}: runs[{}] must be an object".format(path, i))
        for key, typ in _ENTRY_KEYS.items():
            if not isinstance(e.get(key), typ):
                raise ValueError("{}: runs[{}].{} is missing or not a {}".format(
                    path, i, key, "string" if typ is str else "object"))
        for name in (e["run_dir"], e["input_toml"], e["initial"]):
            if os.path.isabs(name) or ".." in _parts(name):
                raise ValueError("{}: runs[{}]: {!r} is not a path inside WORK".format(
                    path, i, name))
        if e["run_id"] in out:
            raise ValueError("{}: runs[{}]: duplicate run_id {!r}; each run is listed "
                             "once".format(path, i, e["run_id"]))
        out[e["run_id"]] = e
    return out


def collect(work, *, stdout=None, stderr=None) -> int:
    """``python -m hwave.initial collect`` (returns the exit code).

    Reads every run of WORK and writes ``coverage.json`` (``post_run``)
    and ``summary.txt``. Status per run, in this order: ``seed_unverified``
    when ``input.toml`` no longer names the prepared ``initial.npz``, a
    staged file or the npz differs from the manifest SHA-256, or the
    initial-state check refuses it; otherwise a valid status file decides
    (``iterations = 0`` gives ``failed``); otherwise ``hwave.log``
    (``not_run``, ``convergence_unverified``, ``failed``, ``not_converged``,
    ``converged``). The final labels take the reference of prepare's
    labels (the para reference kept in ``WORK/para_hf`` in the small tier,
    refused when its files differ from the manifest). A status file that
    exists but cannot be read makes the
    row ``failed`` (``"<file> unreadable: <message>"``): the output readers
    refuse such a run. Outputs of a converged run must be readable, else
    ``failed``. Every reader gets ``base_dir`` = the run directory.
    Idempotent: the result depends only on WORK."""
    from hwave.initial.candidates import CandidateSet, _thaw_json
    from hwave.initial.report import coverage_document, rank, summary_text
    stdout, stderr = stdout or sys.stdout, stderr or sys.stderr
    work = pathlib.Path(os.path.abspath(os.fspath(work)))
    for name in (COVERAGE, SUMMARY):
        _refuse_link(work / name)          # before anything is written
    manifest = _read_json(work / MANIFEST, "manifest")
    if isinstance(manifest, dict) and manifest.get("stage") == STAGE_PENDING:
        raise ValueError("{}: the preliminary paramagnetic run of WORK is pending; run bash "
                         "{}, then repeat the prepare command, before collect".format(
                             work, work / STAGE0_COMMANDS))
    entries = _check_manifest(manifest, work / MANIFEST)
    data = _read_json(work / CANDIDATES, "candidate set")
    if not isinstance(data, dict):
        raise ValueError("{}: a JSON object expected".format(work / CANDIDATES))
    chosen = data.pop("tier_chosen_by", "auto (interim rule)")
    try:
        cset = CandidateSet.from_json(data)
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("{}: not a candidate set of this version: {}: {}".format(
            work / CANDIDATES, type(exc).__name__, exc)) from None
    model = _thaw_json(cset.model)
    rows, warnings, facts_cache, para_cache = [], [], {}, {}
    para_entry = _check_para_entry(manifest, work)
    if para_entry is None and _uses_para(cset):
        raise ValueError("{}: the candidate set was prepared on a paramagnetic Hartree-Fock "
                         "reference (seed base para_hf), but the manifest names no para "
                         "reference file; the labels cannot be computed against it. Repeat "
                         "prepare into a new directory".format(work / MANIFEST))
    if para_entry is not None:
        para_cache["para_dir"] = work / para_entry["dir"]
    for r in cset.rows:
        row = _thaw_json(r)
        if row.get("tried"):
            entry = entries.get(row["run_id"])
            if entry is None:
                raise ValueError("{}: run {} is not in the manifest".format(
                    work, row["run_id"]))
            _collect_row(row, entry, work, model, warnings, facts_cache, para_cache)
        rows.append(row)
    ranking = rank(rows, facts=model)
    cset2 = dataclasses.replace(cset, rows=rows)
    doc = coverage_document(cset2, phase_marker="post_run", tier_chosen_by=chosen,
                            ranking=ranking)
    _write_json(work / COVERAGE, doc)
    text = summary_text(doc, ranking)
    if warnings:
        text += "\nCollect warnings:\n" + "\n".join("  " + w for w in warnings) + "\n"
    _write_text(work / SUMMARY, text)
    for w in warnings:
        print("warning: {}".format(w), file=stderr)
    stdout.write(text)
    return EXIT_OK


# --------------------------------------------------------------------------
# the command line
# --------------------------------------------------------------------------

def _subshape_arg(text):
    try:
        vals = [int(x) for x in text.split(",")]
    except ValueError:
        vals = []
    if len(vals) != 3 or any(v < 1 for v in vals):
        raise argparse.ArgumentTypeError(
            "{!r}: a SubShape is three positive integers, e.g. 2,2,1".format(text))
    return tuple(vals)


def _positive_float(text):
    try:
        v = float(text)
    except ValueError:
        raise argparse.ArgumentTypeError("{!r} is not a number".format(text)) from None
    if not math.isfinite(v) or v <= 0:
        raise argparse.ArgumentTypeError("{!r}: a finite number > 0 expected".format(text))
    return v


def _nonneg_int(text):
    try:
        v = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError("{!r} is not an integer".format(text)) from None
    if v < 0:
        raise argparse.ArgumentTypeError("{!r}: an integer >= 0 expected".format(text))
    return v


def _epilog():
    return "exit codes:\n" + "\n".join("  {}  {}".format(c, t) for c, t in EXIT_CODES)


class _Parser(argparse.ArgumentParser):
    def error(self, message):
        self.print_usage(sys.stderr)
        print("{}: error: {}".format(self.prog, message), file=sys.stderr)
        raise _Usage(message)


def build_parser():
    """The argument parser of ``python -m hwave.initial``."""
    p = _Parser(
        prog="python -m hwave.initial",
        description="Candidate initial states for UHFk (preview). prepare writes one run "
                    "directory per candidate, you run hwave in each (bash "
                    "WORK/commands.txt), collect reads the results back. The selected "
                    "state is the lowest free energy among the candidates tried. It is not "
                    "a proof of the ground state.",
        epilog=_epilog(), formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", metavar="{prepare,check,collect}",
                           parser_class=_Parser)
    sub.required = True
    pp = sub.add_parser(
        "prepare", help="write the candidate runs into WORK", epilog=_epilog(),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Enumerate the candidate initial states of a UHFk input and write "
                    "WORK: coverage.json, manifest.json, commands.txt and one run "
                    "directory per run (input.toml, initial.npz, inputs/ with copies of "
                    "the input files). Comments of the input TOML are not kept.")
    pp.add_argument("input", metavar="INPUT_DIR",
                    help="a directory holding exactly one UHFk *.toml, or the TOML")
    pp.add_argument("--out", required=True, metavar="WORK", help="the work directory")
    pp.add_argument("--subshape", action="append", type=_subshape_arg, metavar="Bx,By,Bz",
                    help="a SubShape to evaluate every candidate on (repeatable); "
                         "default: the input's SubShape")
    pp.add_argument("--tier", choices=("auto", "large", "medium", "small"),
                    help="the enumeration tier: large (spin and charge patterns), small "
                         "(adds orbital patterns and instability hints from a paramagnetic "
                         "Hartree-Fock run); auto chooses small when the pair space |P| <= "
                         "400, the k points N_q <= 4096 and N_q*|P|^3 <= 1e11 (fixed limits, "
                         "the same on every host), else large; the host memory never "
                         "chooses, it only refuses a small tier whose memory estimate does "
                         "not fit; medium arrives in a later release")
    pp.add_argument("--extra-candidates", metavar="FILE.json",
                    help="a JSON list of {family: user|literature, pattern, label, "
                         "citation}; the pattern follows patternspec.schema.json")
    pp.add_argument("--t-hint", type=_positive_float, metavar="T",
                    help="the hint temperature, in input energy units (default: the input "
                         "T, or a rule from the bands and k mesh at T = 0); it changes only "
                         "the instability rows of the small tier, never the other "
                         "candidates")
    pp.add_argument("--amplitude", type=float, metavar="A",
                    help="the requested seed amplitude, in (0, 1) (default 0.5)")
    pp.add_argument("--max-candidates", type=_nonneg_int, metavar="N",
                    help="the generated rows per SubShape (default 32)")
    pp.add_argument("--orbital-site", metavar="FILE",
                    help="a JSON list: one site label per physical orbital")
    pp.add_argument("--para-check", action="store_true",
                    help="small tier: also run the paramagnetic Hartree-Fock step from "
                         "the uniform orbital occupation and keep the lower free energy")
    pp.add_argument("--stage0", action="store_true",
                    help="small tier: run the paramagnetic Hartree-Fock step outside "
                         "prepare, as the preliminary paramagnetic run: prepare writes "
                         "WORK/stage0_para (and WORK/stage0_para_check with --para-check) "
                         "and WORK/stage0_commands.txt and exits with code 3; run bash "
                         "WORK/stage0_commands.txt, then re-run the same prepare command, "
                         "which writes the candidate runs (refused when the large tier is given "
                         "or chosen by auto)")
    pc = sub.add_parser(
        "check", help="check an initial-state file against the input that reads it",
        epilog=_epilog(), formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Check the initial-state npz the input CONFIG.toml reads (or NPZ) "
                    "before hwave runs: provenance, SubShape, shape, Hermiticity, "
                    "eigenvalue bounds, electron number and blocks.")
    pc.add_argument("config", metavar="CONFIG.toml")
    pc.add_argument("npz", nargs="?", metavar="NPZ")
    pc.add_argument("--require-provenance", action="store_true",
                    help="refuse a file without provenance (exit code 4)")
    pl = sub.add_parser(
        "collect", help="read the results of WORK", epilog=_epilog(),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Read every run of WORK and write WORK/coverage.json (post_run) and "
                    "WORK/summary.txt; safe to repeat.")
    pl.add_argument("work", metavar="WORK")
    return p


def main(argv=None) -> int:
    """The entry point of ``python -m hwave.initial``; returns the exit
    code (:data:`EXIT_CODES`)."""
    from hwave.results import ResultsError
    parser = build_parser()
    try:
        args = parser.parse_args(argv)
    except _Usage:
        return EXIT_USAGE
    except SystemExit as exc:          # --help
        return int(exc.code or 0)
    try:
        if args.command == "prepare":
            subs = args.subshape
            if subs is not None and len(set(subs)) != len(subs):
                dup = sorted({b for b in subs if subs.count(b) > 1})
                try:
                    parser.error("--subshape {} is given more than once".format(
                        ",".join(map(str, dup[0]))))
                except _Usage:
                    return EXIT_USAGE
            return prepare(args.input, args.out, subshapes=subs, tier=args.tier,
                           extra_candidates=args.extra_candidates, t_hint=args.t_hint,
                           amplitude=args.amplitude, max_candidates=args.max_candidates,
                           orbital_site=args.orbital_site, para_check=args.para_check,
                           stage0=args.stage0)
        if args.command == "check":
            return check(args.config, args.npz, require_provenance=args.require_provenance)
        return collect(args.work)
    except (ValueError, NotImplementedError, OSError, ResultsError) as exc:
        print("error: {}".format(exc), file=sys.stderr)
        return EXIT_ERROR
