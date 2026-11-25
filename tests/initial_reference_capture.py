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


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--facade", action="store_true",
                        help="capture the facade golden file")
    parser.add_argument("--out", default=FACADE_GOLDEN, help="output .npz path")
    args = parser.parse_args(argv)
    if not args.facade:
        parser.error("nothing to do: pass --facade")
    if _REPO not in sys.path:
        sys.path.insert(0, _REPO)
    capture_facade(args.out)


if __name__ == "__main__":
    main()
