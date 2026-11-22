"""The convergence status of a UHFr or UHFk run.

:func:`read_scf_status` reads, in this order:

1. ``scf_status.dat`` (``[file.output] scf_status``, on by default, ``""``
   disables it), format version 1, written by
   ``solver_base._save_scf_status`` (src/hwave/solver/base.py:524-586):
   nine ``key = value`` lines, ``format_version``, ``solver``,
   ``converged``, ``residual``, ``eps``, ``iterations``, ``iteration_max``,
   ``residual_definition`` and ``note``. Every key is required, once; a bad
   value, a broken invariant or another version raises
   UnsupportedLayoutError naming the key.
2. When there is no status file: the ``[log] print_check`` trace, if one is
   configured and present (uhfk.py:816-869, uhfr.py:575-626): one line
   ``step, Rest, energy, NCond, Sz`` per SCF step, Rest to 8 significant
   digits.
3. Otherwise the presence of the output files.
"""
from __future__ import annotations

import dataclasses
import math
import numbers
import re
from pathlib import Path
from typing import Optional

from hwave.results._io import read_key_value_text
from hwave.results.errors import MissingOutputError, UnsupportedLayoutError
from hwave.results.run import as_run

__all__ = ["ScfStatus", "read_scf_status", "require_scf_steps"]

FORMAT_VERSION = 1

#: the keys of format version 1, in the order they are written (base.py:564-574)
KEYS = ("format_version", "solver", "converged", "residual", "eps", "iterations",
        "iteration_max", "residual_definition", "note")

#: ``residual_definition`` per solver (uhfk.py:1792-1793, uhfr.py:1267-1268):
#: ``2 * ||G_new - G_prev||_F / G.size`` over the whole Green function
RESIDUAL_DEFINITIONS = {"UHFk": "two_frobenius_over_size:uhfk_all_r",
                        "UHFr": "two_frobenius_over_size:uhfr_full"}

#: the relative width of the band around ``eps`` in which a ``print_check``
#: residual, rounded to 8 significant digits, cannot be compared with ``eps``
PRINT_CHECK_BAND = 1e-7

NO_STEP_NOTE = "no SCF step performed"

NOT_RUN_NOTE = (
    "no output file was found; the run may not have started, or it stopped "
    "before writing any output (qlms.run creates the output directory before "
    "solving)")

PREDATES_NOTE = (
    "this output predates the convergence status file and no print_check file "
    "was configured; the status is only in the run log")

_STATUS_DISABLED = 'the convergence status file is disabled (scf_status = "")'
_PREDATES = "this output predates the convergence status file"
_NO_PRINT_CHECK = "no print_check file was configured"
_PRINT_CHECK_ABSENT = "the configured print_check file is absent"

_NOT_FINITE = re.compile(r"residual not finite after step ([1-9][0-9]*)")
_NEGATIVE_MAX = re.compile(r"IterationMax (-[1-9][0-9]*) treated as 0")
_COUNT = re.compile(r"[0-9]+")

_CONVERGED = {"converged": True, "not_converged": False}


@dataclasses.dataclass(frozen=True)
class ScfStatus:
    """The convergence status of a run.

    Attributes
    ----------
    status : str
        ``"converged"``, ``"not_converged"``, ``"failed"`` (the SCF ran but
        ``energy`` is configured and its file is absent), ``"unknown"`` or
        ``"not_run"`` (no configured output file exists).
    converged : bool or None
        True for ``"converged"``, False for ``"not_converged"``, else None.
    residual : float or None
        The residual of the last step; NaN when no step was performed;
        ``None`` when it was not finite (``meta["residual_not_finite"]``)
        or is not known.
    eps : float or None
        The threshold ``10**-EPS``.
    iterations, iteration_max : int or None
        The steps performed and the step limit (a negative IterationMax is
        reported as 0, as the writer does).
    residual_definition : str or None
        The residual as the status file names it; ``None`` from other
        sources.
    source_kind : str
        ``"scf_status"``, ``"print_check"`` or ``"none"``.
    note : str
        The note of the status file, or what was concluded and why.
    source : Path or None
        The file read; ``None`` for ``source_kind == "none"``.
    meta : dict
        From a status file: ``format_version``, ``solver`` and
        ``note_kind`` (``""``, ``"no_step"``, ``"residual_not_finite"`` with
        ``nonfinite_step``, or ``"iteration_max_negative"`` with
        ``iteration_max_requested``). From ``print_check``:
        ``residual_digits`` (8), ``incomplete_line`` (a cut last line),
        ``iteration_max_requested`` (a negative IterationMax, reported as 0)
        and, for ``"failed"``, ``scf_verdict``.
        ``residual_not_finite`` is True when the last residual was NaN or
        infinite.
    """
    status: str
    converged: Optional[bool]
    residual: Optional[float]
    eps: Optional[float]
    iterations: Optional[int]
    iteration_max: Optional[int]
    residual_definition: Optional[str]
    source_kind: str
    note: str
    source: Optional[Path]
    meta: dict


def read_scf_status(run_or_toml, *, base_dir=None) -> ScfStatus:
    """The convergence status of a run (a :class:`hwave.results.Run` or the
    path of its input TOML file); see the module text for the sources.

    From ``print_check``, the last residual decides: ``converged`` below
    ``eps * (1 - 1e-7)``, ``not_converged`` above ``eps * (1 + 1e-7)``,
    ``unknown`` in between (the residual is rounded to 8 digits). An empty
    file means that no step was performed (``IterationMax = 0``); it cannot
    be told from a run that stopped during its first step. ``failed`` when
    steps were performed but the configured ``energy`` file is absent. A
    last line without its newline was cut while being written (a hard
    kill): it is dropped, the run is ``not_converged`` and
    ``meta["incomplete_line"]`` gives its number. A negative
    ``IterationMax`` is reported as 0, as the writer does
    (``meta["iteration_max_requested"]``).

    Without records: ``not_run`` when no configured output file exists,
    ``failed`` when ``energy`` is configured and absent but another output
    exists, else ``unknown``. A run of zero steps with ``scf_status = ""``
    and no ``print_check`` writes nothing and is therefore reported as
    ``not_run``.

    Raises UnsupportedLayoutError for an unreadable status or ``print_check``
    file, and for a mode other than UHFr and UHFk."""
    run = as_run(run_or_toml, base_dir=base_dir)
    if run.mode not in RESIDUAL_DEFINITIONS:
        raise UnsupportedLayoutError(
            "the convergence status of mode {!r} is not read by this version of "
            "hwave.results".format(run.mode), path=run.input_path)
    status_path = run.outputs.get("scf_status")
    if status_path is not None and status_path.is_file():
        return _from_status_file(run, status_path)
    check_path = run.outputs.get("print_check")
    if check_path is not None and check_path.is_file():
        return _from_print_check(run, check_path)
    return _from_outputs(run)


def require_scf_steps(run, output):
    """Refuse to read ``output`` (a key of ``run.outputs``) of a run whose
    status file records zero SCF steps.

    After zero steps the solver writes ``scf_status.dat`` and returns
    before any other writer (base.py ``_save_scf_status``; uhfk.py:1783-1794,
    uhfr.py:1263-1266). Any other output file found in the output directory
    is therefore left over from an earlier run, whatever its time stamp.
    Every reader except :func:`read_scf_status` calls this before it opens
    its file.

    When the run's ``scf_status`` file exists, it decides: the rule applies
    when it records ``iterations = 0`` (a status file that cannot be read
    raises UnsupportedLayoutError as :func:`read_scf_status` does). Without
    a status file (``scf_status = ""``, or an older output) the
    configuration decides: ``[mode.param] IterationMax <= 0`` (``false``
    counts as 0, as ``int(IterationMax)`` in the solver) means the
    solver performed no step (its ``solve`` loops ``range(IterationMax)``,
    uhfk.py:835, uhfr.py:598), so it wrote no other output in this run.
    Raises MissingOutputError, with a hint naming the status file or
    IterationMax."""
    path = run.outputs.get(output)
    status_path = run.outputs.get("scf_status")
    if status_path is None or not status_path.is_file():
        steps = run.param("IterationMax")
        shown = repr(steps)
        if isinstance(steps, bool):
            # the solver takes int(IterationMax) (base.py:552) and loops
            # range(IterationMax): false is 0 steps, true is 1
            shown = "false" if not steps else "true"
            steps = int(steps)
        if not isinstance(steps, numbers.Real) or not steps <= 0:
            return
        raise MissingOutputError(
            "output {!r} is not from this run: [mode.param] IterationMax = {} "
            "in {}, so the run performed no SCF step and the solver wrote no "
            "output{}".format(output, shown, run.input_path, _leftover(path)),
            path=path if path is not None else run.input_path,
            hint="[mode.param] IterationMax = {} performs no SCF step; rerun "
                 "with IterationMax > 0".format(shown))
    status = _from_status_file(run, status_path)
    if status.iterations != 0:
        return
    raise MissingOutputError(
        "output {!r} is not from this run: the run performed no SCF step, and "
        "the solver then writes no output other than the status file{}".format(
            output, _leftover(path)),
        path=path if path is not None else status_path,
        hint="{} records iterations = 0 ({}); rerun with IterationMax > 0".format(
            status_path, status.note or "no SCF step performed"))


def _leftover(path):
    if path is None or not path.exists():
        return ""
    return " ({} is left over from an earlier run)".format(path)


# ------------------------------------------------------------ status file


def _from_status_file(run, path):
    def bad(key, message):
        return UnsupportedLayoutError("{}: {!r}: {}".format(path, key, message),
                                      path=path)

    rows = read_key_value_text(path)
    table, lines = {}, {}
    for key, value, lineno in rows:
        if key in table:
            raise bad(key, "line {} repeats the key of line {}".format(
                lineno, lines[key]))
        table[key], lines[key] = value, lineno
    if "format_version" not in table:
        raise bad("format_version", "the key is missing; this is not a convergence "
                  "status file this reader knows")
    if table["format_version"] != str(FORMAT_VERSION):
        raise UnsupportedLayoutError(
            "{}: {!r}: version {!r}; this version of hwave.results reads version {}"
            .format(path, "format_version", table["format_version"], FORMAT_VERSION),
            path=path, hint="update hwave.results")
    for key in table:
        if key not in KEYS:
            raise bad(key, "unknown key in format version {}".format(FORMAT_VERSION))
    for key in KEYS:
        if key not in table:
            raise bad(key, "the key is missing")

    solver = table["solver"]
    if solver != run.mode:
        raise bad("solver", "{!r}, but the input is a {} run".format(solver, run.mode))
    converged_text = table["converged"]
    if converged_text not in ("true", "false"):
        raise bad("converged", "{!r} is not 'true' or 'false'".format(converged_text))
    converged = converged_text == "true"

    def number(key):
        try:
            return float(table[key])
        except ValueError:
            raise bad(key, "{!r} is not a number".format(table[key])) from None

    def count(key):
        if not _COUNT.fullmatch(table[key]):
            raise bad(key, "{!r} is not a non-negative integer".format(table[key]))
        return int(table[key])

    residual = number("residual")
    eps = number("eps")
    if not (math.isfinite(eps) and eps > 0):
        raise bad("eps", "{!r} is not a finite number greater than 0".format(table["eps"]))
    iterations = count("iterations")
    iteration_max = count("iteration_max")
    if table["residual_definition"] != RESIDUAL_DEFINITIONS[solver]:
        raise bad("residual_definition", "{!r}; {} writes {!r}".format(
            table["residual_definition"], solver, RESIDUAL_DEFINITIONS[solver]))
    if iterations > iteration_max:
        raise bad("iterations", "{} steps, more than iteration_max = {}".format(
            iterations, iteration_max))

    note = table["note"]
    meta = {"format_version": FORMAT_VERSION, "solver": solver}
    not_finite = _NOT_FINITE.fullmatch(note)
    negative_max = _NEGATIVE_MAX.fullmatch(note)
    if note == "":
        meta["note_kind"] = ""
    elif note == NO_STEP_NOTE:
        meta["note_kind"] = "no_step"
    elif not_finite:
        meta["note_kind"] = "residual_not_finite"
        meta["nonfinite_step"] = int(not_finite.group(1))
    elif negative_max:
        meta["note_kind"] = "iteration_max_negative"
        meta["iteration_max_requested"] = int(negative_max.group(1))
    else:
        raise bad("note", "{!r} is not a note of format version {}".format(
            note, FORMAT_VERSION))
    kind = meta["note_kind"]

    if iterations == 0:
        if converged:
            raise bad("converged", "true, but no SCF step was performed")
        if not math.isnan(residual):
            raise bad("residual", "{!r}, but no SCF step was performed (expected "
                      "nan)".format(table["residual"]))
        if kind not in ("no_step", "iteration_max_negative"):
            raise bad("note", "{!r}, but no SCF step was performed".format(note))
        if kind == "iteration_max_negative" and iteration_max != 0:
            raise bad("iteration_max", "{}, but the note says it was treated as 0"
                      .format(iteration_max))
    elif kind == "residual_not_finite":
        if meta["nonfinite_step"] != iterations or not math.isnan(residual):
            raise bad("note", "{!r} does not match iterations = {} and residual = "
                      "{}".format(note, iterations, table["residual"]))
        if converged:
            raise bad("converged", "true, but the residual was not finite")
        meta["residual_not_finite"] = True
        residual = None
    elif kind != "":
        raise bad("note", "{!r}, but {} SCF steps were performed".format(
            note, iterations))
    elif not (math.isfinite(residual) and residual >= 0):
        raise bad("residual", "{!r} is not a finite number >= 0 (a non-finite "
                  "residual is written as nan with the note 'residual not finite "
                  "after step N')".format(table["residual"]))
    if converged and not residual < eps:
        raise bad("converged", "true, but residual = {} is not below eps = {}"
                  .format(table["residual"], table["eps"]))

    status = "converged" if converged else "not_converged"
    return ScfStatus(status=status, converged=converged, residual=residual, eps=eps,
                     iterations=iterations, iteration_max=iteration_max,
                     residual_definition=table["residual_definition"],
                     source_kind="scf_status", note=note, source=Path(path),
                     meta=meta)


# ------------------------------------------------------------- print_check


def _config_eps(run):
    """``10**-EPS`` from the input, checked as the solver checks it
    (base.py:237-253)."""
    exponent = run.param("EPS")
    usable = False
    eps = None
    if isinstance(exponent, (int, float)) and not isinstance(exponent, bool):
        try:
            eps = float(pow(10, -exponent))
            usable = math.isfinite(float(exponent)) and math.isfinite(eps) and eps > 0
        except (OverflowError, ValueError):
            usable = False
    if not usable:
        raise UnsupportedLayoutError(
            "{}: {!r}: {!r} gives no usable convergence threshold 10^(-EPS)".format(
                run.input_path, "EPS", exponent), path=run.input_path)
    return eps


def _from_print_check(run, path):
    eps = _config_eps(run)
    meta = {"residual_digits": 8}
    iteration_max = run.param("IterationMax")
    if isinstance(iteration_max, bool):
        iteration_max = int(iteration_max)     # as the writer does (base.py:552)
    if not isinstance(iteration_max, int):
        iteration_max = None
    elif iteration_max < 0:
        # the writer's rule (base.py:552-556): a negative IterationMax runs
        # no step and is reported as 0
        meta["iteration_max_requested"] = iteration_max
        iteration_max = 0
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise UnsupportedLayoutError("{} cannot be read as text: {}".format(path, exc),
                                     path=path) from None
    common = dict(eps=eps, iteration_max=iteration_max, residual_definition=None,
                  source_kind="print_check", source=Path(path))
    if "iteration_max_requested" in meta:
        clamp = "; IterationMax {} treated as 0".format(meta["iteration_max_requested"])
    else:
        clamp = ""
    if text == "":
        return ScfStatus(status="not_converged", converged=False,
                         residual=float("nan"), iterations=0,
                         note=NO_STEP_NOTE + clamp, meta=meta, **common)
    lines = text.splitlines()
    incomplete = None
    if not text.endswith("\n"):
        # every line is written with its newline (uhfk.py:850, uhfr.py:609):
        # a last line without one was cut, e.g. by a hard kill
        incomplete = len(lines)
        meta["incomplete_line"] = incomplete
        lines = lines[:-1]
    residual = float("nan")
    rest_text = "nan"
    for lineno, line in enumerate(lines, start=1):
        fields = line.split(",")
        if len(fields) != 5:
            raise UnsupportedLayoutError(
                "{}: line {} has {} comma-separated fields, not 5 (step, Rest, "
                "energy, NCond, Sz): {!r}".format(path, lineno, len(fields), line),
                path=path)
        rest_text = fields[1].strip()
        try:
            residual = float(rest_text)
        except ValueError:
            raise UnsupportedLayoutError(
                "{}: line {}: the residual {!r} is not a number".format(
                    path, lineno, rest_text), path=path) from None
    iterations = len(lines)
    if incomplete is not None:
        note = ("the print_check trace ends with an incomplete line {} (the run was "
                "stopped while writing it); {} complete SCF steps".format(
                    incomplete, iterations))
        if not iterations:
            residual = float("nan")       # no step performed
        else:
            note += ", the last residual {}".format(rest_text)
            if not math.isfinite(residual):
                meta["residual_not_finite"] = True
                residual = None
        return ScfStatus(status="not_converged", converged=False,
                         residual=residual,
                         iterations=iterations, note=note + clamp, meta=meta,
                         **common)
    if not math.isfinite(residual):
        verdict = "not_converged"        # a non-finite residual never converges
        meta["residual_not_finite"] = True
        residual = None
    elif residual < eps * (1 - PRINT_CHECK_BAND):
        verdict = "converged"
    elif residual > eps * (1 + PRINT_CHECK_BAND):
        verdict = "not_converged"
    else:
        verdict = "unknown"
    note = "from the print_check trace: {} SCF steps, the last residual {} " \
        "against eps = {!r}".format(iterations, rest_text, eps)
    if verdict == "unknown":
        note += "; the residual is too close to eps to decide at 8 significant digits"
    note += clamp
    energy = run.outputs.get("energy")
    if energy is not None and not energy.is_file():
        meta["scf_verdict"] = verdict
        return ScfStatus(
            status="failed", converged=None, residual=residual, iterations=iterations,
            note="{}; but the energy file {} is absent, so the run stopped before it "
                 "saved its results".format(note, energy.name),
            meta=meta, **common)
    return ScfStatus(status=verdict, converged=_CONVERGED.get(verdict),
                     residual=residual, iterations=iterations, note=note, meta=meta,
                     **common)


# ------------------------------------------------------------ no records


def _from_outputs(run):
    common = dict(converged=None, residual=None, eps=None, iterations=None,
                  iteration_max=None, residual_definition=None, source_kind="none",
                  source=None)
    present = sorted(key for key, path in run.outputs.items() if path.is_file())
    meta = {"present_outputs": tuple(present)}
    if not present:
        return ScfStatus(status="not_run", note=NOT_RUN_NOTE, meta=meta, **common)
    energy = run.outputs.get("energy")
    if energy is not None and not energy.is_file():
        return ScfStatus(
            status="failed",
            note="the energy file {} is absent although other outputs exist ({}); "
                 "the run stopped before it saved its results".format(
                     energy.name, ", ".join(present)),
            meta=meta, **common)
    first = _PREDATES if "scf_status" in run.outputs else _STATUS_DISABLED
    second = _PRINT_CHECK_ABSENT if "print_check" in run.outputs else _NO_PRINT_CHECK
    note = "{} and {}; the status is only in the run log".format(first, second)
    return ScfStatus(status="unknown", note=note, meta=meta, **common)
