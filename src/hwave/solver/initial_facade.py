"""The solver entry point of ``hwave.initial``.

``hwave.initial`` builds candidate initial states for UHFk. It needs the
model as the solver sees it (lattice, orbitals, one-body H(k), interaction
tables, block structure) and the solver's Hartree-Fock potential, without
running a self-consistent calculation. This module is the only solver
module it imports; it does not import ``hwave.initial``.

* :func:`load_problem` reads an input TOML layer by layer as ``hwave`` does
  and returns the frozen :class:`UHFkProblem`.
* :func:`solver_blocks`, :func:`hf_tables`: the block structure and the
  interaction tables of a problem.
* :func:`accumulate_hf_projected`, :func:`hf_potential`: the Hartree-Fock
  potential UHFk adds to ``T(k)``, from the stored Green function.

Every function takes and returns the stored Green function of UHFk,
``G[r, s, a, t, b] = <c^dag_{sa}(r) c_{tb}(0)>`` (the layout of the
``green`` output), and Hamiltonians in the solver's index order.
"""
import copy
import os
import pathlib
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Optional, Tuple

import numpy as np

from . import hamiltonian_k
from .uhfk_blocks import HFTables, SolverBlocks, accumulate_hf_projected

__all__ = ["UHFkProblem", "SolverBlocks", "HFTables", "load_problem",
           "solver_blocks", "hf_tables", "accumulate_hf_projected",
           "hf_potential"]

#: the top-level tables ``hwave`` reads, by exact name (qlms.py:98-104)
_TOP_TABLES = ("mode", "file", "log")
#: the keys of ``[file]`` ``hwave`` reads (qlms.py:104-111)
_FILE_KEYS = ("input", "output")
#: the ``[file.input]`` keys of UHFk, matched exactly
#: (qlms.py:108, read_input_k.py:93, :163-201)
_INPUT_KEYS = ("path_to_input", "interaction", "initial", "initial_uhf",
               "initial_mode", "geometry_uhf", "onebodyg_uhf")
#: the ``[file.output]`` keys of UHFk, matched exactly (qlms.py:112,
#: base.scf_status_file_name, UHFk.save_results)
_OUTPUT_KEYS = ("path_to_output", "energy", "eigen", "green", "onebodyg",
                "initial", "rpa", "export_hamiltonian", "scf_status")


# eq=False: the fields hold arrays (identity comparison and hash)
@dataclass(frozen=True, eq=False)
class UHFkProblem:
    """A UHFk model as the solver sets it up, before any iteration.

    Lattice: ``cellshape`` (the lattice), ``subshape`` (the unit cell the
    solver works in), ``shape`` (the k-grid, ``cellshape / subshape``),
    ``nvol`` (its size), ``norb`` (orbitals of the folded cell; spin
    included in spin-orbital mode), ``norb_orig`` (Geometry ``norb``),
    ``norb_phys`` (physical orbitals of the folded cell), ``ns``, ``nd``
    (``ns * norb``).

    Modes: ``enable_spin_orbital``, ``iflag_fock``, ``sz_free``,
    ``nconds`` (``(Ncond,)`` or ``(N_up, N_down)``), ``T``,
    ``spin_resolved_transfer`` (a 9-column Transfer, normal mode).

    One-body: ``ham_trans_k`` ``(nvol, nd, nd)``, the solver's ``T(k)``
    (the same array UHFk adds), and ``ham_trans_r`` ``(nx, ny, nz, nd,
    nd)``, the folded Transfer table as a dense array in the same index
    order (``hamiltonian_k.to_k(ham_trans_r)`` is ``ham_trans_k``).

    Interactions: ``tables`` (:class:`HFTables`) and ``param_ham_orig``,
    the tables of the original (primitive) cell as the solver normalised
    them, before the sublattice fold.

    Geometry: ``geometry`` (of the original cell), ``orbital_site`` (the
    ``[mode.param] orbital_site`` labels, ``None`` when absent or in
    spin-orbital mode).

    ``blocks`` (:class:`SolverBlocks`) and ``input_files`` (the absolute
    paths of the Geometry, Transfer and interaction files read).
    """
    cellshape: Tuple[int, int, int]
    subshape: Tuple[int, int, int]
    shape: Tuple[int, int, int]
    nvol: int
    norb: int
    norb_orig: int
    norb_phys: int
    ns: int
    nd: int
    enable_spin_orbital: bool
    iflag_fock: bool
    sz_free: bool
    nconds: Tuple[int, ...]
    T: float
    spin_resolved_transfer: bool
    ham_trans_k: np.ndarray
    ham_trans_r: np.ndarray
    tables: HFTables
    param_ham_orig: Mapping
    geometry: Mapping
    orbital_site: Optional[Tuple[int, ...]]
    blocks: SolverBlocks
    input_files: Tuple[pathlib.Path, ...]


def _case_variants(table, keys, where, is_table=False):
    """Messages for the entries of ``table`` that equal one of ``keys``
    only when case is ignored (``hwave`` would not read them)."""
    lower = {k.lower(): k for k in keys}
    out = []
    for key in table:
        if key not in keys and key.lower() in lower:
            if is_table:
                out.append("[{}] is not read: table names are case-sensitive, "
                           "write it as [{}].".format(key, lower[key.lower()]))
            else:
                out.append("{} {!r} is not read: keys of {} are case-sensitive, "
                           "write it as {!r}.".format(where, key, where,
                                                     lower[key.lower()]))
    return out


def _read_config(toml_path):
    """The parsed TOML, refused (``ValueError``) where ``hwave`` would
    ignore a key for its case or exit on the placement of a key."""
    import tomli
    from ..qlms import check_mode_tables
    with open(toml_path, "rb") as f:
        doc = tomli.load(f)
    errors = _case_variants(doc, _TOP_TABLES, "", is_table=True)
    info_mode = doc.get("mode", {})
    if not isinstance(info_mode, Mapping):
        errors.append("[mode] must be a table")
        info_mode = {}
    errors += check_mode_tables(info_mode)
    if errors:
        raise ValueError("{}: {}".format(toml_path, " ".join(errors)))
    if "mode" not in info_mode:
        raise ValueError("{}: mode is not defined in [mode].".format(toml_path))
    if info_mode["mode"] != "UHFk":
        raise ValueError("{}: [mode] mode = {!r}: only UHFk inputs are supported"
                         .format(toml_path, info_mode["mode"]))
    # the [file] layers by the UHFk rules (UHFr reads some keys in any case)
    info_file = doc.get("file", {})
    if isinstance(info_file, Mapping):
        errors += _case_variants(info_file, _FILE_KEYS, "[file]")
        for name, keys in (("input", _INPUT_KEYS), ("output", _OUTPUT_KEYS)):
            tab = info_file.get(name, {})
            if isinstance(tab, Mapping):
                errors += _case_variants(tab, keys, "[file.{}]".format(name))
    if errors:
        raise ValueError("{}: {}".format(toml_path, " ".join(errors)))
    return doc


def _build_solver(toml_path, *, subshape=None, with_green=False):
    """The UHFk instance of ``toml_path`` after the ``solve`` prologue
    (``T(k)``, the interaction tables and the blocks), as ``hwave`` would
    build it with its working directory at the TOML's directory.

    Returns ``(solver, green_info, input_files)``: ``green_info`` is the
    green-function input (``initial`` and the other ``[file.input]`` files,
    read only when ``with_green``; otherwise empty), ``input_files`` the
    absolute paths of the Geometry, Transfer and interaction files read.
    """
    from requests.structures import CaseInsensitiveDict
    from ..qlmsio.read_input_k import QLMSkInput
    from .base import ParamDict
    from .uhfk import UHFk

    toml_path = pathlib.Path(os.path.abspath(os.fspath(toml_path)))
    doc = _read_config(toml_path)
    base_dir = str(toml_path.parent)

    info_log = dict(doc.get("log", {}))
    info_log["print_level"] = info_log.get("print_level", 1)
    info_log["print_step"] = info_log.get("print_step", 1)

    info_mode = dict(doc["mode"])
    info_mode["param"] = ParamDict(info_mode.get("param", {}) or {})
    if subshape is not None:
        info_mode["param"]["SubShape"] = [int(n) for n in subshape]

    finput = dict(doc.get("file", {}).get("input", {}))
    # the solver joins every path with its cwd; hwave runs at the TOML's
    # directory, so the two path_to_input values are taken from there
    # (os.path.join keeps an absolute value as it is)
    interaction = CaseInsensitiveDict(finput.get("interaction", {}))
    interaction["path_to_input"] = os.path.join(
        base_dir, interaction.get("path_to_input", ""))
    info_inputfile = {"path_to_input": os.path.join(base_dir, finput.get("path_to_input", "")),
                      "interaction": interaction}
    if with_green:
        for key in ("initial", "initial_uhf", "initial_mode", "geometry_uhf",
                    "onebodyg_uhf"):
            if key in finput:
                info_inputfile[key] = finput[key]

    try:
        read_io = QLMSkInput(info_inputfile)
        ham_info = read_io.get_param("ham")
        green_info = read_io.get_param("green") if with_green else {}
        solver = UHFk(ham_info, info_log, info_mode)
        solver._make_ham_trans()
        solver._make_ham_inter()
        solver._detect_blocks()
    except SystemExit as exc:
        raise ValueError("{}: UHFk refused the input (exit status {}); the logged "
                         "error says why.".format(toml_path, exc.code)) from None
    return solver, green_info, _interaction_files(interaction)


def _interaction_files(interaction):
    """The Geometry, Transfer and interaction files of
    ``[file.input.interaction]``, as ``QLMSkInput`` joins them, in the TOML
    key order."""
    base = interaction["path_to_input"]
    return tuple(pathlib.Path(os.path.abspath(os.path.join(base, name)))
                 for key, name in interaction.items()
                 if key.lower() != "path_to_input")


def _orbital_site(solver):
    if solver.enable_spin_orbital:
        # spin-orbital mode never parses the key (as the RPA solver)
        return None
    from .sites import parse_orbital_site
    return parse_orbital_site(solver.param_mod.get("orbital_site"), solver.norb_phys_orig)


def _ham_trans_r(solver):
    """The folded Transfer table, dense, in the index order of
    ``solver.ham_trans``."""
    table = solver.param_ham["Transfer"]
    if solver.enable_spin_orbital or solver.spin_resolved_transfer:
        return hamiltonian_k.dense_table(table, solver.shape, solver.nd)
    tab = hamiltonian_k.dense_table(table, solver.shape, solver.norb)
    nx, ny, nz = solver.shape
    norb = solver.norb
    out = np.zeros((nx, ny, nz, 2 * norb, 2 * norb), dtype=np.complex128)
    for s in range(2):
        out[..., s * norb:(s + 1) * norb, s * norb:(s + 1) * norb] = tab
    return out


def _read_only(a):
    a = np.array(a, copy=True)
    a.flags.writeable = False
    return a


def load_problem(toml_path, *, subshape=None):
    """The UHFk model of the input ``toml_path``, without solving.

    The TOML is read layer by layer as ``hwave`` reads it: ``[mode]``,
    ``[file]``, ``[file.input]`` and ``[file.output]`` by exact key (a case
    variant of a known key is refused, naming the solver spelling;
    misplaced ``[mode]`` keys are refused with the messages of
    ``qlms.check_mode_tables``), ``[mode.param]`` and
    ``[file.input.interaction]`` case-insensitively. Relative paths resolve
    from the TOML's directory, where ``hwave`` runs. Each call reads the
    input files anew, so two problems never share a table.

    Parameters
    ----------
    toml_path : path-like
    subshape : tuple of int, optional
        Overrides ``[mode.param] SubShape``.

    Returns
    -------
    UHFkProblem

    Raises
    ------
    ValueError
        A non-UHFk mode, a refused key, or an input the solver refuses.
    """
    solver, _, input_files = _build_solver(toml_path, subshape=subshape)
    if solver.has_sublattice:
        param_ham_orig = solver.param_ham_orig
    else:
        param_ham_orig = copy.deepcopy(solver.param_ham)
    blocks = SolverBlocks(tuple(solver.block_info), tuple(solver.block_to_group),
                          tuple(solver.group_nconds))
    tables = HFTables(inter_table=MappingProxyType(dict(solver.inter_table)),
                      spin_table=MappingProxyType(dict(solver.spin_table)),
                      shape=tuple(solver.shape),
                      enable_spin_orbital=solver.enable_spin_orbital,
                      norb_phys=solver.norb_phys)
    return UHFkProblem(
        cellshape=tuple(solver.cellshape), subshape=tuple(solver.subshape),
        shape=tuple(solver.shape), nvol=solver.nvol, norb=solver.norb,
        norb_orig=solver.norb_orig, norb_phys=solver.norb_phys, ns=solver.ns,
        nd=solver.nd, enable_spin_orbital=solver.enable_spin_orbital,
        iflag_fock=solver.iflag_fock, sz_free=solver.sz_free,
        nconds=tuple(solver.Nconds), T=solver.T,
        spin_resolved_transfer=solver.spin_resolved_transfer,
        ham_trans_k=_read_only(solver.ham_trans), ham_trans_r=_read_only(_ham_trans_r(solver)),
        tables=tables, param_ham_orig=param_ham_orig,
        geometry=param_ham_orig["Geometry"], orbital_site=_orbital_site(solver),
        blocks=blocks, input_files=input_files)


def solver_blocks(problem):
    """The block structure UHFk solves ``problem`` with: the lifted
    detection (:func:`hwave.solver.uhfk_blocks.detect_blocks`) on the
    problem's ``T(k)`` and interaction tables."""
    from .uhfk_blocks import detect_blocks
    return detect_blocks(
        problem.ham_trans_k, problem.tables.inter_table, problem.tables.spin_table,
        nd=problem.nd, norb=problem.norb, ns=problem.ns, iflag_fock=problem.iflag_fock,
        sz_free=problem.sz_free, enable_spin_orbital=problem.enable_spin_orbital,
        nconds=problem.nconds)


def hf_tables(problem):
    """The interaction tables of ``problem`` (:class:`HFTables`)."""
    return problem.tables


def hf_potential(green, tables, *, include_fock, blocks):
    """The Hartree-Fock potential ``Sigma'(k)`` UHFk adds to ``T(k)`` for
    the stored Green function ``green`` ``(nvol, ns, norb, ns, norb)``:
    :func:`accumulate_hf_projected` on zeros. Returns ``(nvol, nd, nd)``
    complex128 in the Hamiltonian convention and the solver's index
    order."""
    green = np.asarray(green)
    nvol, ns, norb = green.shape[0], green.shape[1], green.shape[2]
    nd = ns * norb
    out = np.zeros((nvol, nd, nd), dtype=np.complex128)
    accumulate_hf_projected(out, green, tables, include_fock=include_fock, blocks=blocks)
    return out
