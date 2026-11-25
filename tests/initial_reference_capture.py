"""Reference data for the ``hwave.initial`` tests, captured at a fixed revision.

``--facade`` writes ``tests/initial_ref/facade_golden_d72f35f6.npz``: for each
input of :func:`facade_cases`, the block structure UHFk detects
(``block_info``, ``block_to_group``, ``group_nconds``), the Hamiltonian of
the first iteration (``_make_ham`` on the initial Green function) and the
Green function after the first ``_green``. It must run in a checkout at
:data:`REVISION` whose ``src`` tree is unmodified, with that checkout's
``src`` first on ``PYTHONPATH``; it refuses anything else. To capture, copy
this file into the ``tests`` directory of a clean checkout at the revision
and run, from that checkout::

    PYTHONPATH=$PWD/src:$PWD python -B tests/initial_reference_capture.py \\
        --facade --out <path>.npz

The capture uses only solver code present at the revision (``QLMSkInput``,
``UHFk`` and its methods), so the golden values are those of the solver
before the facade existed. The tests rebuild the same inputs with
:func:`facade_cases`.

``--writer`` writes ``tests/initial_ref/writer_reference.npz``, the
initial-state file of :func:`writer_reference_case` as
``hwave.initial.writer.initial_state`` writes it on the current tree (the
cross-numpy content reference of ``tests/test_initial_writer.py``)::

    PYTHONPATH=$PWD/src:$PWD python -B tests/initial_reference_capture.py --writer

``--validation`` writes ``tests/initial_ref/v1a_i1.json`` and
``tests/initial_ref/v1b.json``, the frozen references of the validation
cases V1a and V1b of ``tests/test_initial_validation_i1.py`` (the V1 model
and solver settings are :data:`V1_MODEL` and :data:`V1_SOLVER`; each file
records them with the H-wave and numpy versions it was captured with)::

    PYTHONPATH=$PWD/src:$PWD python -B tests/initial_reference_capture.py --validation

They are captured once; a later run that contradicts them is a finding to
report, not a reason to re-capture.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

import numpy as np

#: the revision the golden file is captured at (the branch base of PR I1)
REVISION = "d72f35f6"

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: the committed golden file of ``--facade``
FACADE_GOLDEN = os.path.join(_REPO, "tests", "initial_ref",
                             "facade_golden_{}.npz".format(REVISION))


# ---------------------------------------------------------------------------
# A writer for the TOML subset of H-wave inputs (tomli only reads)
# ---------------------------------------------------------------------------

def _toml_value(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, float):
        return repr(v)
    if isinstance(v, str):
        return json.dumps(v)
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(_toml_value(x) for x in v) + "]"
    raise TypeError("cannot write {!r} as TOML".format(v))


def dump_toml(doc):
    """``doc`` (nested dicts of scalars and flat lists) as TOML text."""
    out = []

    def table(name, tab):
        scalars = [(k, v) for k, v in tab.items() if not isinstance(v, dict)]
        subs = [(k, v) for k, v in tab.items() if isinstance(v, dict)]
        if name:
            out.append("[{}]".format(name))
        for k, v in scalars:
            out.append("{} = {}".format(k, _toml_value(v)))
        for k, v in subs:
            table("{}.{}".format(name, k) if name else k, v)

    table("", doc)
    return "\n".join(out) + "\n"


def _read_toml(path):
    import tomli
    with open(path, "rb") as f:
        return tomli.load(f)


# ---------------------------------------------------------------------------
# The inputs
# ---------------------------------------------------------------------------

def _copy_case(src_rel, dest, edit=None):
    shutil.copytree(os.path.join(_REPO, src_rel), dest)
    for name in ("output",):
        shutil.rmtree(os.path.join(dest, name), ignore_errors=True)
    toml = os.path.join(dest, "input.toml")
    if edit is not None:
        doc = _read_toml(toml)
        edit(doc)
        with open(toml, "w") as f:
            f.write(dump_toml(doc))
    return toml


def _write_so_case(dest):
    """The spin-orbital input of tests/test_uhfk_so_make_ham.py (every
    interaction type, complex one-body table), SubShape 2 on a 6-site chain,
    random initial Green function."""
    from tests.test_uhfk_so_make_ham import NSO, _write
    os.makedirs(dest)
    files = _write(dest)
    inter = {"path_to_input": "./", "Geometry": "geom.dat", "Transfer": "transfer.dat"}
    inter.update(files)
    doc = {"log": {"print_level": 0},
           "mode": {"mode": "UHFk", "enable_spin_orbital": True, "flag_fock": True,
                    "param": {"Ncond": NSO, "T": 0.2, "CellShape": [6, 1, 1],
                              "SubShape": [2, 1, 1], "RndSeed": 3,
                              "IterationMax": 1, "EPS": 8, "Mix": 0.5}},
           "file": {"input": {"path_to_input": "", "initial_mode": "random",
                              "interaction": inter},
                    "output": {"path_to_output": "output"}}}
    toml = os.path.join(dest, "input.toml")
    with open(toml, "w") as f:
        f.write(dump_toml(doc))
    return toml


def _write_degenerate_case(dest):
    """4x4 square lattice, nearest-neighbour hopping, U = 0.5, half filling,
    T = 0, zero initial Green function: the first Hamiltonian is the bare
    band, whose Fermi shell (12 states for 6 electrons) is exactly
    degenerate, so the first ``_green`` depends on the tie order of the
    T = 0 occupation sort."""
    os.makedirs(dest)
    with open(os.path.join(dest, "geom.dat"), "w") as f:
        f.write("1 0 0\n0 1 0\n0 0 1\n1\n0 0 0\n")
    rows = [((1, 0, 0), -1.0), ((-1, 0, 0), -1.0), ((0, 1, 0), -1.0), ((0, -1, 0), -1.0)]
    with open(os.path.join(dest, "transfer.dat"), "w") as f:
        f.write("Transfer\n1\n1\n 1\n")
        for (rx, ry, rz), v in rows:
            f.write("{} {} {} 1 1 {:.12f} 0.0\n".format(rx, ry, rz, v))
    with open(os.path.join(dest, "coulombintra.dat"), "w") as f:
        f.write("CoulombIntra\n1\n1\n 1\n0 0 0 1 1 0.5 0.0\n")
    doc = {"log": {"print_level": 0},
           "mode": {"mode": "UHFk",
                    "param": {"Ncond": 16, "T": 0.0, "CellShape": [4, 4, 1],
                              "SubShape": [1, 1, 1], "IterationMax": 1,
                              "EPS": 8, "Mix": 0.5, "RndSeed": 1}},
           "file": {"input": {"path_to_input": "",
                              "interaction": {"path_to_input": "./",
                                              "Geometry": "geom.dat",
                                              "Transfer": "transfer.dat",
                                              "CoulombIntra": "coulombintra.dat"}},
                    "output": {"path_to_output": "output"}}}
    toml = os.path.join(dest, "input.toml")
    with open(toml, "w") as f:
        f.write(dump_toml(doc))
    return toml


#: spin-resolved (9-column) two-orbital chain with spin-flip hopping
_W90SPIN_ROWS = {((0, 0, 0), 1, 1, 0, 0): -0.3, ((0, 0, 0), 1, 1, 1, 1): -0.2,
                 ((0, 0, 0), 2, 2, 0, 0): 0.3, ((0, 0, 0), 2, 2, 1, 1): 0.35,
                 ((1, 0, 0), 1, 1, 0, 0): -1.0, ((1, 0, 0), 1, 1, 1, 1): -1.0,
                 ((1, 0, 0), 2, 2, 0, 0): -0.7, ((1, 0, 0), 2, 2, 1, 1): -0.7,
                 ((1, 0, 0), 1, 2, 0, 0): -0.2 + 0.05j,
                 ((1, 0, 0), 1, 2, 1, 1): -0.2 + 0.05j,
                 ((1, 0, 0), 1, 2, 0, 1): 0.08j, ((1, 0, 0), 2, 1, 1, 0): -0.08j}


def _write_w90spin_case(dest):
    """A spin-resolved 9-column Transfer (#323) with spin-flip hopping,
    CoulombIntra, SubShape 2 on a 4-site chain, random initial."""
    from tests.test_w90_spin_columns_uhfk import (
        _write_coulomb_intra, _write_geom, _write_spin_transfer)
    os.makedirs(dest)
    _write_geom(dest, 2)
    _write_coulomb_intra(dest, 2, 2.0)
    _write_spin_transfer(dest, 2, _W90SPIN_ROWS)
    doc = {"log": {"print_level": 0},
           "mode": {"mode": "UHFk",
                    "param": {"Ncond": 6, "T": 0.2, "CellShape": [4, 1, 1],
                              "SubShape": [2, 1, 1], "IterationMax": 1,
                              "EPS": 8, "Mix": 0.5, "RndSeed": 7}},
           "file": {"input": {"path_to_input": "", "initial_mode": "random",
                              "interaction": {"path_to_input": "./",
                                              "Geometry": "geom.dat",
                                              "Transfer": "transfer.dat",
                                              "CoulombIntra": "coulombintra.dat"}},
                    "output": {"path_to_output": "output"}}}
    toml = os.path.join(dest, "input.toml")
    with open(toml, "w") as f:
        f.write(dump_toml(doc))
    return toml


def _subshape_2x2(doc):
    doc["mode"]["param"]["SubShape"] = [2, 2, 1]
    finput = doc["file"]["input"]
    finput.pop("initial", None)
    finput["initial_mode"] = "random"


def _fock_off(doc):
    doc["mode"]["flag_fock"] = False


#: ``{name: builder(dest) -> input.toml path}``, in a fixed order
_CASES = (
    # the five UHFk cases of tests/test_uhf_outputs_reference.py CASES
    ("uhfk_coulombintra", lambda d: _copy_case("tests/uhfk/CoulombIntra", d)),
    ("uhfk_2sz0", lambda d: _copy_case("tests/uhfk/2sz0", d)),
    ("uhfk_exchange", lambda d: _copy_case("tests/uhfk/Exchange", d)),
    ("uhfk_coulombinter", lambda d: _copy_case("tests/uhfk/CoulombInter", d)),
    ("tutorial_hubbard_uhfk", lambda d: _copy_case("docs/tutorial/Hubbard/UHFk", d)),
    ("uhfk_2sz1", lambda d: _copy_case("tests/uhfk/2sz1", d)),
    ("so_all_types", _write_so_case),
    ("exchange_subshape_2x2", lambda d: _copy_case("tests/uhfk/Exchange", d, _subshape_2x2)),
    ("pairhop_fock_off", lambda d: _copy_case("tests/uhfk/PairHop", d, _fock_off)),
    ("degenerate_shell_t0", _write_degenerate_case),
    ("w90spin_flip", _write_w90spin_case),
)

FACADE_CASE_NAMES = tuple(name for name, _ in _CASES)


def facade_cases(root):
    """Write every input under ``root/<name>/``; ``{name: input.toml path}``."""
    return {name: build(os.path.join(root, name)) for name, build in _CASES}


# ---------------------------------------------------------------------------
# Capture
# ---------------------------------------------------------------------------

def build_solver(toml_path):
    """The UHFk instance of ``toml_path`` as ``hwave`` builds it (cwd = the
    TOML's directory), after the ``solve`` prologue up to the block
    detection; and the green-function input."""
    import hwave.qlmsio as qlmsio
    import hwave.solver.uhfk as sol_uhfk
    from hwave.solver.base import ParamDict
    doc = _read_toml(toml_path)
    cwd = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(toml_path)))
    try:
        info_log = dict(doc.get("log", {}))
        info_log["print_level"] = info_log.get("print_level", 1)
        info_log["print_step"] = info_log.get("print_step", 1)
        info_mode = dict(doc["mode"])
        info_mode["param"] = ParamDict(info_mode.get("param", {}) or {})
        info_inputfile = dict(doc["file"]["input"])
        info_inputfile["path_to_input"] = info_inputfile.get("path_to_input", "")
        read_io = qlmsio.read_input_k.QLMSkInput(info_inputfile)
        ham_info = read_io.get_param("ham")
        green_info = read_io.get_param("green")
        solver = sol_uhfk.UHFk(ham_info, info_log, info_mode)
        solver._make_ham_trans()
        solver._make_ham_inter()
        solver._detect_blocks()
        solver.Green = solver._initial_green(green_info)
    finally:
        os.chdir(cwd)
    return solver, green_info


def capture_case(toml_path):
    """``{key: array}`` of one input (see the module docstring)."""
    solver, _ = build_solver(toml_path)
    blocks = [np.asarray(b, dtype=np.int64) for b in solver.block_info]
    out = {
        "block_sizes": np.array([b.size for b in blocks], dtype=np.int64),
        "block_index": np.concatenate(blocks),
        "block_to_group": np.array(solver.block_to_group, dtype=np.int64),
        "group_nconds": np.array(solver.group_nconds, dtype=np.int64),
    }
    solver._make_ham()
    out["ham"] = np.array(solver.ham, dtype=np.complex128)
    solver._diag()
    solver._green()
    out["green"] = np.array(solver.Green, dtype=np.complex128)
    return out


def _check_checkout():
    """Refuse to capture unless HEAD is :data:`REVISION`, ``src`` is
    unmodified, and ``hwave`` is imported from this checkout."""
    head = subprocess.run(["git", "-C", _REPO, "rev-parse", "HEAD"],
                          capture_output=True, text=True, check=True).stdout.strip()
    want = subprocess.run(["git", "-C", _REPO, "rev-parse", REVISION + "^{commit}"],
                          capture_output=True, text=True, check=True).stdout.strip()
    if head != want:
        sys.exit("refusing to capture: HEAD is {}, not {} ({})".format(head, REVISION, want))
    dirty = subprocess.run(["git", "-C", _REPO, "status", "--porcelain", "--", "src"],
                           capture_output=True, text=True, check=True).stdout
    if dirty.strip():
        sys.exit("refusing to capture: src is modified:\n" + dirty)
    import hwave
    src = os.path.realpath(os.path.join(_REPO, "src"))
    if not os.path.realpath(hwave.__file__).startswith(src + os.sep):
        sys.exit("refusing to capture: hwave is imported from {}, not {}".format(
            hwave.__file__, src))
    return head


def capture_facade(out_path):
    import tempfile
    head = _check_checkout()
    root = tempfile.mkdtemp(prefix="initial_capture_")
    try:
        arrays = {"revision": np.array(head), "numpy_version": np.array(np.__version__),
                  "cases": np.array(list(FACADE_CASE_NAMES))}
        for name, toml in facade_cases(root).items():
            for key, value in capture_case(toml).items():
                arrays["{}__{}".format(name, key)] = value
    finally:
        shutil.rmtree(root, ignore_errors=True)
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    np.savez_compressed(out_path, **arrays)
    print("wrote {} at revision {}".format(out_path, head))


#: the committed reference of ``--writer``
WRITER_REFERENCE = os.path.join(_REPO, "tests", "initial_ref", "writer_reference.npz")


def writer_reference_case(dest):
    """The input and candidate of the writer reference: a 4 x 4 Hubbard
    square lattice (hopping 1, U = 2, T = 0.1, half filling), SubShape
    2 x 2, a Neel seed. ``(toml, candidate)``."""
    from hwave.initial import Candidate
    from hwave.initial.patterns import Entry, PatternSpec
    from tests import initial_fixtures
    toml = initial_fixtures.square(dest, 4, U=2.0, T=0.1, subshape=(2, 2, 1))
    spec = PatternSpec(q=("1/2", "1/2", "0"), subshape=(2, 2, 1), phase=0.0, sign=1,
                       modulation="cos", entries=(Entry(0, 0, (0, 0, 0), cz=1.0),))
    return toml, Candidate.user(pattern=spec, label="Neel")


def capture_writer(out_path):
    import tempfile
    from hwave.initial.writer import initial_state
    root = tempfile.mkdtemp(prefix="initial_capture_")
    try:
        toml, cand = writer_reference_case(os.path.join(root, "case"))
        os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
        initial_state(toml, cand, out=out_path)
    finally:
        shutil.rmtree(root, ignore_errors=True)
    print("wrote {} with numpy {}".format(out_path, np.__version__))


# ---------------------------------------------------------------------------
# The validation cases V1a and V1b (tests/test_initial_validation_i1.py)
# ---------------------------------------------------------------------------

#: the committed references of ``--validation``
V1A_REFERENCE = os.path.join(_REPO, "tests", "initial_ref", "v1a_i1.json")
V1B_REFERENCE = os.path.join(_REPO, "tests", "initial_ref", "v1b.json")

#: the V1 model: square lattice 16 x 16, hopping 1, U = 4, half filling,
#: SubShape 2 x 2, Fock on; the solver settings of every run
V1_MODEL = {"L": 16, "U": 4.0, "subshape": [2, 2, 1]}
V1_SOLVER = {"IterationMax": 3000, "EPS": 10, "Mix": 0.5}
#: V1a temperature; V1b is T = 0
V1A_T = 0.02
#: the T_N bisection: a Neel seed at T, ordered when max |m_i| > the
#: threshold; ``steps`` halvings of [lo, hi]; IterationMax of those runs
TN_BISECTION = {"lo": 0.02, "hi": 1.0, "steps": 12, "threshold": 1e-6,
                "IterationMax": 400000, "EPS": 14}


def v1_input(dest, T, *, solver=None, finput=None):
    """The V1 input at temperature ``T`` under ``dest`` (``solver``
    overrides :data:`V1_SOLVER`); the TOML path."""
    from tests import initial_fixtures
    param = dict(V1_SOLVER, **(solver or {}))
    return initial_fixtures.write_input(
        dest, cellshape=(V1_MODEL["L"], V1_MODEL["L"], 1),
        transfer=initial_fixtures.square_transfer(),
        interactions={"CoulombIntra": {((0, 0, 0), 0, 0): V1_MODEL["U"]}},
        param=dict(param, Ncond=V1_MODEL["L"] ** 2, T=T,
                   SubShape=list(V1_MODEL["subshape"])),
        finput=finput)


def run_solver(toml_path):
    """Run ``hwave`` on ``toml_path`` in this process, in the TOML's
    directory (as ``commands.txt`` does), its standard output into
    ``hwave.log`` there."""
    import contextlib
    import hwave.qlms
    run_dir = os.path.dirname(os.path.abspath(toml_path))
    cwd = os.getcwd()
    os.chdir(run_dir)
    try:
        with open("hwave.log", "w") as log, contextlib.redirect_stdout(log):
            hwave.qlms.run(input_file=os.path.basename(toml_path))
    finally:
        os.chdir(cwd)


def prepare_run_collect(toml_path, work, **prepare_kw):
    """``prepare`` WORK for ``toml_path``, run every run directory of its
    manifest with :func:`run_solver`, ``collect``; the post-run coverage
    document."""
    import io
    from hwave.initial import cli
    rc = cli.prepare(toml_path, work, stdout=io.StringIO(), **prepare_kw)
    if rc != 0:
        raise RuntimeError("prepare returned {}".format(rc))
    with open(os.path.join(work, "manifest.json")) as f:
        manifest = json.load(f)
    for entry in manifest["runs"]:
        run_solver(os.path.join(work, entry["run_dir"], entry["input_toml"]))
    rc = cli.collect(work, stdout=io.StringIO(), stderr=io.StringIO())
    if rc != 0:
        raise RuntimeError("collect returned {}".format(rc))
    with open(os.path.join(work, "coverage.json")) as f:
        return json.load(f)


def max_moment(row):
    """The largest local spin moment |m_i| of a collected row."""
    return max(row["order_parameters"]["per_site"]["moment"])


def neel_candidate():
    """The collinear Neel pattern at q = (1/2, 1/2) on 2 x 2, as a user row."""
    from hwave.initial import Candidate
    from hwave.initial.patterns import Entry, PatternSpec
    spec = PatternSpec(q=("1/2", "1/2", "0"), subshape=(2, 2, 1), phase=0.0, sign=1,
                       modulation="cos", entries=(Entry(0, 0, (0, 0, 0), cz=1.0),))
    return Candidate.user(pattern=spec, label="Neel")


def neel_run(dest, T, *, solver=None):
    """Run the V1 model at ``T`` from the Neel seed; ``(max |m_i|,
    converged)``."""
    from hwave.initial.writer import initial_state
    from hwave.results import read_scf_status, read_uhfk_order_parameters
    toml = v1_input(dest, T, solver=solver, finput={"initial": "initial.npz"})
    initial_state(toml, neel_candidate(), out=os.path.join(dest, "initial.npz"))
    run_solver(toml)
    op = read_uhfk_order_parameters(toml)
    return float(np.max(op.moment)), bool(read_scf_status(toml).converged)


def bisect_t_n(root, *, lo, hi, steps, threshold, IterationMax, EPS):
    """T_N of the V1 model: bisection of [lo, hi] on max |m_i| > threshold
    from the Neel seed. ``(lo, hi, record)``: ordered at lo, not at hi.

    The runs use a tighter ``EPS`` than the candidate runs: near T_N the
    moment decays by a factor close to 1 per iteration, so at EPS = 10 a run
    stops with a remnant moment above 1e-6 well above T_N."""
    record = []

    def ordered(T):
        m, conv = neel_run(os.path.join(root, "T{:d}".format(len(record))), T,
                           solver={"IterationMax": IterationMax, "EPS": EPS})
        record.append({"T": T, "moment": m, "converged": conv})
        return m > threshold

    if not ordered(lo) or ordered(hi):
        raise RuntimeError("T_N is not bracketed by [{}, {}]: {}".format(lo, hi, record))
    for _ in range(steps):
        mid = 0.5 * (lo + hi)
        if ordered(mid):
            lo = mid
        else:
            hi = mid
    return lo, hi, record


def stoner_t_n():
    """An independent T_N of the V1 model: the root of U chi0(Q, T) = 1 at
    Q = (1/2, 1/2), chi0 the bare static Lindhard function per spin on the
    L x L mesh (the para Hartree-Fock bands of one orbital at half filling
    are the bare bands shifted by U/2, and E_F = 0)."""
    from scipy.optimize import brentq
    L, U = V1_MODEL["L"], V1_MODEL["U"]
    k = 2.0 * np.pi * np.arange(L) / L
    kx, ky = np.meshgrid(k, k, indexing="ij")
    e = -2.0 * (np.cos(kx) + np.cos(ky))     # e(k + Q) = -e(k)

    def chi0(T):
        f = 0.5 * (1.0 - np.tanh(e / (2.0 * T)))
        safe = np.where(np.abs(e) < 1e-12, 1.0, e)
        # (f(e) - f(-e)) / (-2 e); the limit at e = 0 is f (1 - f) / T
        v = np.where(np.abs(e) < 1e-12, f * (1.0 - f) / T, (2.0 * f - 1.0) / (-2.0 * safe))
        return float(v.mean())

    return brentq(lambda T: U * chi0(T) - 1.0, 0.05, 2.0, xtol=1e-12)


def _rows_record(cov):
    """``{run_id: {status, final_label, free_energy_per_site, moment}}`` of
    the tried rows of a post-run coverage document."""
    out = {}
    for r in cov["rows"]:
        if not r["tried"]:
            continue
        conv = r["status"] == "converged"
        out[r["run_id"]] = {
            "status": r["status"], "final_label": r["final_label"],
            "free_energy_per_site": r["free_energy_per_site"],
            "moment": max_moment(r) if conv else None}
    return out


def _provenance_record(description):
    import hwave
    return {"description": description, "hwave_version": hwave.__version__,
            "numpy_version": np.__version__,
            "captured_by": "tests/initial_reference_capture.py --validation",
            "model": dict(V1_MODEL, hopping=1.0, filling="half", fock=True),
            "solver": dict(V1_SOLVER)}


def capture_validation(out_dir):
    """Write ``v1a_i1.json`` and ``v1b.json`` under ``out_dir``."""
    import tempfile
    root = tempfile.mkdtemp(prefix="initial_capture_")
    try:
        cov = prepare_run_collect(v1_input(os.path.join(root, "a_src"), V1A_T),
                                  os.path.join(root, "a_work"))
        rows = _rows_record(cov)
        neel = [rid for rid in rows if rid.startswith("baseline-q1_2_1_2_0-spin-z")
                and "-neg-" not in rid]
        if len(neel) != 1 or rows[neel[0]]["status"] != "converged":
            raise RuntimeError("V1a: the Neel row did not converge: {}".format(neel))
        lo, hi, record = bisect_t_n(os.path.join(root, "tn"), **TN_BISECTION)
        v1a = _provenance_record(
            "V1a: T = 0.02, the default large-tier candidate set; the Neel row's "
            "moment max |m_i|; T_N = the upper end of a bisection on max |m_i| > "
            "1e-6 from the Neel seed (ordered at T_N_lower, not at T_N); "
            "T_N_stoner = the root of U chi0(Q, T) = 1 on the same mesh, an "
            "independent check")
        v1a.update(T=V1A_T, neel_run_id=neel[0], neel_moment=rows[neel[0]]["moment"],
                   T_N=hi, T_N_lower=lo, T_N_bisection=dict(TN_BISECTION),
                   T_N_stoner=stoner_t_n(),
                   T_N_runs=record, rows=rows)
        cov = prepare_run_collect(v1_input(os.path.join(root, "b_src"), 0.0),
                                  os.path.join(root, "b_work"))
        v1b = _provenance_record(
            "V1b: T = 0 (open Fermi shell), the default large-tier candidate set; "
            "status, final label, free energy per site and max |m_i| of every "
            "tried row")
        v1b.update(T=0.0, T_hint=cov["hint"]["T_hint"], rows=_rows_record(cov))
    finally:
        shutil.rmtree(root, ignore_errors=True)
    os.makedirs(out_dir, exist_ok=True)
    for name, doc in (("v1a_i1.json", v1a), ("v1b.json", v1b)):
        path = os.path.join(out_dir, name)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(doc, f, indent=1, sort_keys=True, ensure_ascii=False)
            f.write("\n")
        print("wrote {} with numpy {}".format(path, np.__version__))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--facade", action="store_true",
                        help="capture the facade golden file")
    parser.add_argument("--writer", action="store_true",
                        help="write the initial-state writer reference")
    parser.add_argument("--validation", action="store_true",
                        help="write the V1a and V1b validation references")
    parser.add_argument("--out", default=None,
                        help="output .npz path (--validation: output directory)")
    args = parser.parse_args(argv)
    if args.facade + args.writer + args.validation != 1:
        parser.error("pass exactly one of --facade, --writer and --validation")
    if _REPO not in sys.path:
        sys.path.insert(0, _REPO)
    if args.facade:
        capture_facade(args.out or FACADE_GOLDEN)
    elif args.writer:
        capture_writer(args.out or WRITER_REFERENCE)
    else:
        capture_validation(args.out or os.path.dirname(V1A_REFERENCE))


if __name__ == "__main__":
    main()
