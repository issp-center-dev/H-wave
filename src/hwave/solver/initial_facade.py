"""The solver entry point of ``hwave.initial``.

``hwave.initial`` builds candidate initial states for UHFk. It needs the
model as the solver sees it (lattice, orbitals, one-body H(k), interaction
tables, block structure) and the solver's Hartree-Fock potential, without
running a self-consistent calculation. This module is the only solver
module it imports; it does not import ``hwave.initial``.

* :func:`load_problem` reads an input TOML layer by layer as ``hwave`` does
  and returns the frozen :class:`UHFkProblem`. Relative paths are taken
  from the TOML's directory (``hwave`` takes them from its working
  directory; the two agree when ``hwave`` runs in the TOML's directory).
* :func:`solver_blocks`, :func:`hf_tables`: the block structure and the
  interaction tables of a problem.
* :func:`accumulate_hf_projected`, :func:`hf_potential`: the Hartree-Fock
  potential UHFk adds to ``T(k)``, from the stored Green function.
* :func:`green_function`: the stored Green function UHFk builds from a
  Hamiltonian (its diagonalisation and occupation fill), or with a uniform
  fill of an open T = 0 Fermi shell.
* :func:`referenced_input_files`: the files an input TOML names, by key.
* :func:`primitive_h0_k`: the bare ``T(k)`` of the original (primitive)
  cell on the ``CellShape`` mesh, built by the solver's
  :mod:`hamiltonian_k` (one Fourier convention).
* Re-exported pure helpers: :func:`to_k` (``hamiltonian_k.to_k``, the
  solver's R -> k transform ``T(k) = sum_R e^{+ikR} T(R)``) and
  :func:`interleaved_index` (``spin_orbital.interleaved_index``, the
  spin-orbital index ``2*orb + spin``).

Every function takes and returns the stored Green function of UHFk,
``G[r, s, a, t, b] = <c^dag_{sa}(r) c_{tb}(0)>`` (the layout of the
``green`` output), and Hamiltonians in the solver's index order.
"""
import os
import pathlib
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Optional, Tuple

import numpy as np

from . import hamiltonian_k
from .hamiltonian_k import to_k
from .spin_orbital import interleaved_index
from .uhfk_blocks import HFTables, SolverBlocks, accumulate_hf_projected
from .uhfk_fill import FillInfo, MuNotConvergedError, ShellRecord

__all__ = ["UHFkProblem", "SolverBlocks", "HFTables", "FillInfo", "ShellRecord",
           "MuNotConvergedError", "TWO_BODY_TYPES", "load_problem",
           "solver_blocks", "hf_tables", "accumulate_hf_projected",
           "hf_potential", "green_function", "referenced_input_files",
           "primitive_h0_k", "to_k", "interleaved_index"]

#: the two-body interaction tables UHFk reads (hartree_fock.HF_TYPE_ORDER
#: and the combined ``Coulomb`` table, uhfk.py ``_check_interaction``)
TWO_BODY_TYPES = ("CoulombIntra", "CoulombInter", "Coulomb", "Hund", "Ising",
                  "PairLift", "Exchange", "PairHop")

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
    ``ncond`` (``Ncond`` as the solver resolves it, ``Filling`` included),
    ``two_sz`` (``2Sz``, ``None`` when Sz is free), ``nconds`` (``(Ncond,)``
    or ``(N_up, N_down)``), ``T``, ``ene_cutoff`` (the Fermi-function cut),
    ``spin_resolved_transfer`` (a 9-column Transfer, normal mode).

    Every field is read-only: arrays have ``writeable = False`` and
    mappings are read-only views of copies.

    One-body: ``ham_trans_k`` ``(nvol, nd, nd)``, the solver's ``T(k)``
    (the same array UHFk adds), and ``ham_trans_r`` ``(nx, ny, nz, nd,
    nd)``, the folded Transfer table as a dense array in the same index
    order (``hamiltonian_k.to_k(ham_trans_r)`` is ``ham_trans_k``).

    Interactions: ``tables`` (:class:`HFTables`) and ``param_ham_orig``,
    the tables of the original (primitive) cell as the solver normalised
    them, before the sublattice fold.

    Geometry: ``geometry`` (of the original cell), ``orbital_site`` (the
    ``[mode.param] orbital_site`` labels, one per physical orbital of the
    original cell in every mode; ``None`` when absent).

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
    ncond: int
    two_sz: Optional[int]
    nconds: Tuple[int, ...]
    T: float
    ene_cutoff: float
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
    """``[mode.param] orbital_site``: one site label per physical orbital of
    the original cell, in spin-orbital mode too (``sites.parse_orbital_site``);
    a spin-orbital-length list is refused with a hint."""
    from .sites import parse_orbital_site
    raw = solver.param_mod.get("orbital_site")
    if (solver.enable_spin_orbital and raw is not None and hasattr(raw, "__len__")
            and not isinstance(raw, (str, bytes)) and len(raw) == solver.norb_orig
            and solver.norb_orig != solver.norb_phys_orig):
        raise ValueError(
            "[mode.param] orbital_site has {} entries, the spin-orbital count; in "
            "spin-orbital mode give one label per physical orbital ({} entries, "
            "orbital a covering spin-orbitals 2a and 2a+1)".format(
                len(raw), solver.norb_phys_orig))
    return parse_orbital_site(raw, solver.norb_phys_orig)


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


class FrozenCaseInsensitiveMap(Mapping):
    """A read-only mapping with the case-insensitive string keys of
    ``requests.structures.CaseInsensitiveDict`` (the solver's table
    container)."""

    def __init__(self, items):
        from requests.structures import CaseInsensitiveDict
        self._store = CaseInsensitiveDict(items)

    def __getitem__(self, key):
        return self._store[key]

    def __iter__(self):
        return iter(self._store)

    def __len__(self):
        return len(self._store)

    def __repr__(self):
        return "FrozenCaseInsensitiveMap({!r})".format(dict(self._store.items()))


def _deep_frozen(obj):
    """A read-only deep copy of ``obj``: arrays copied and made
    read-only, mappings as read-only mappings (case-insensitive ones stay
    case-insensitive), lists as tuples."""
    from requests.structures import CaseInsensitiveDict
    if isinstance(obj, np.ndarray):
        return _read_only(obj)
    if isinstance(obj, CaseInsensitiveDict):
        return FrozenCaseInsensitiveMap({k: _deep_frozen(v) for k, v in obj.items()})
    if isinstance(obj, Mapping):
        return MappingProxyType({k: _deep_frozen(v) for k, v in obj.items()})
    if isinstance(obj, (list, tuple)):
        return tuple(_deep_frozen(v) for v in obj)
    return obj


def load_problem(toml_path, *, subshape=None):
    """The UHFk model of the input ``toml_path``, without solving.

    The TOML is read layer by layer as ``hwave`` reads it: ``[mode]``,
    ``[file]``, ``[file.input]`` and ``[file.output]`` by exact key (a case
    variant of a known key is refused, naming the solver spelling;
    misplaced ``[mode]`` keys are refused with the messages of
    ``qlms.check_mode_tables``), ``[mode.param]`` and
    ``[file.input.interaction]`` case-insensitively. Each call reads the
    input files anew, so two problems never share a table.

    **Paths.** Relative ``path_to_input`` values and file names are taken
    relative to the directory of ``toml_path``, whatever the current
    directory is. ``hwave`` itself takes them relative to its working
    directory, so the problem is the model ``hwave`` solves when ``hwave``
    runs in the TOML's directory (as the run directories of
    ``hwave.initial`` are set up to); started from another directory,
    ``hwave`` may read other files.

    Every array and mapping of the returned problem is a read-only copy
    (``param_ham_orig`` is a case-insensitive read-only mapping; the
    reader's ``source_*`` attributes are not kept).

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
    # every container of the problem is a read-only copy (arrays included),
    # so no consumer can change the model another call sees
    param_ham_orig = _deep_frozen(
        solver.param_ham_orig if solver.has_sublattice else solver.param_ham)
    blocks = SolverBlocks(tuple(solver.block_info), tuple(solver.block_to_group),
                          tuple(solver.group_nconds))
    tables = HFTables(inter_table=_deep_frozen(dict(solver.inter_table)),
                      spin_table=_deep_frozen(dict(solver.spin_table)),
                      shape=tuple(solver.shape),
                      enable_spin_orbital=solver.enable_spin_orbital,
                      norb_phys=solver.norb_phys)
    return UHFkProblem(
        cellshape=tuple(solver.cellshape), subshape=tuple(solver.subshape),
        shape=tuple(solver.shape), nvol=solver.nvol, norb=solver.norb,
        norb_orig=solver.norb_orig, norb_phys=solver.norb_phys, ns=solver.ns,
        nd=solver.nd, enable_spin_orbital=solver.enable_spin_orbital,
        iflag_fock=solver.iflag_fock, sz_free=solver.sz_free,
        ncond=solver.param_mod["Ncond"], two_sz=solver.param_mod["2Sz"],
        nconds=tuple(solver.Nconds), T=solver.T, ene_cutoff=solver.ene_cutoff,
        spin_resolved_transfer=solver.spin_resolved_transfer,
        ham_trans_k=_read_only(solver.ham_trans), ham_trans_r=_read_only(_ham_trans_r(solver)),
        tables=tables, param_ham_orig=param_ham_orig,
        geometry=param_ham_orig["Geometry"], orbital_site=_orbital_site(solver),
        blocks=blocks, input_files=input_files)


def primitive_h0_k(problem):
    """The bare one-body ``T(k)`` of the original (primitive) cell of
    ``problem`` on the ``CellShape`` mesh, ``(prod(CellShape), nd_orig,
    nd_orig)``.

    It is :func:`hamiltonian_k.one_body_k` (the builder UHFk uses) on
    ``param_ham_orig["Transfer"]`` with ``shape = CellShape``, in the index
    order UHFk would use for that cell: spin-block ``s*norb_orig + a`` in
    normal mode (a 7-column Transfer embedded over spin, a 9-column one as
    read), interleaved ``2*orb + spin`` in spin-orbital mode. With
    ``SubShape = [1, 1, 1]`` it equals ``problem.ham_trans_k``; with a
    larger ``SubShape`` it has the same spectrum over the mesh.
    """
    table = problem.param_ham_orig["Transfer"]
    norb = problem.norb_orig
    shape = tuple(int(n) for n in problem.cellshape)
    if problem.enable_spin_orbital:
        return hamiltonian_k.one_body_k(table, shape, norb)
    if problem.spin_resolved_transfer:
        return hamiltonian_k.one_body_k(table, shape, 2 * norb)
    return hamiltonian_k.one_body_k(table, shape, norb, spin_embed=True)


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
    if green.ndim != 5:
        raise ValueError("hf_potential: green must be (nvol, ns, norb, ns, norb), "
                         "got shape {}".format(green.shape))
    nvol, ns, norb = green.shape[0], green.shape[1], green.shape[2]
    nd = ns * norb
    out = np.zeros((nvol, nd, nd), dtype=np.complex128)
    # accumulate_hf_projected checks green and the tables against out
    accumulate_hf_projected(out, green, tables, include_fock=include_fock, blocks=blocks)
    return out


def green_function(ham_k, T, blocks, *, shape, ns, norb, fill="solver",
                   shell_tol=None, ene_cutoff=None, return_info=False):
    """The stored Green function UHFk builds from the Hamiltonian ``ham_k``.

    ``ham_k`` ``(nvol, nd, nd)`` is a Hamiltonian in the solver's index
    order (for example ``T(k)`` plus :func:`hf_potential`), ``blocks`` the
    :class:`SolverBlocks` it is block-diagonal in, ``shape`` the k-grid,
    ``ns`` and ``norb`` the stored layout.

    Returns ``G[r, s, a, t, b] = <c^dag_{sa}(r) c_{tb}(0)>`` ``(nvol, ns,
    norb, ns, norb)``, the array UHFk stores and writes as ``green``. Its
    k-space form is ``G_ab(k) = <c^dag_a c_b>``, the orbital transpose of
    the density matrix ``rho(k) = sum_l f_l |v_l><v_l|``: the expectation
    of ``sum O_ab c^dag_a c_b`` is ``sum_ab O_ab G_ab``.

    * ``fill="solver"``: bit for bit UHFk's ``_diag`` and ``_green``, one
      ``mu`` per mu-group; at T = 0 the ``ncond`` lowest levels of a group,
      ties ordered by ``|k|^2`` and then by block and band order.
    * ``fill="uniform"`` (T = 0; needs ``shell_tol``): the levels within
      ``shell_tol`` of the group's Fermi level (its ``ncond``-th lowest
      level) form the shell; the levels below are filled and each shell
      level by ``(n_group - n_below) / g_shell``. The result does not
      depend on the eigenbasis inside the shell. At T > 0 it is the
      solver's fill.

    ``ene_cutoff`` is the solver's Fermi-function cut (``[mode.param]
    ene_cutoff``, :attr:`UHFkProblem.ene_cutoff`); ``None`` takes the
    solvers' default (``hwave.solver.base.PARAM_DEFAULTS``), so pass the
    problem's value for an input that sets it. With ``return_info`` the
    result is ``(G, FillInfo)``: ``mu`` per group (at T = 0 the midpoint of
    the HOMO-LUMO gap, as the solver) and one :class:`ShellRecord` per
    group (at T = 0 the shell of ``shell_tol``, exact ties when it is
    ``None``; closed at T > 0).

    Raises :class:`MuNotConvergedError` when the T > 0 ``mu`` search fails.
    """
    from .base import PARAM_DEFAULTS
    from .uhfk_fill import diag_blocks, fill_green
    if ene_cutoff is None:
        ene_cutoff = PARAM_DEFAULTS["ene_cutoff"]
    ham_k = np.asarray(ham_k)
    nvol = int(np.prod(shape))
    nd = ns * norb
    if ham_k.shape != (nvol, nd, nd):
        raise ValueError("green_function: ham_k must be (nvol, nd, nd) = {}, got {}"
                         .format((nvol, nd, nd), ham_k.shape))
    ws, vs = diag_blocks(ham_k, blocks.block_info)
    out = fill_green(ws, vs, blocks.block_info, blocks.block_to_group,
                     blocks.group_nconds, shape=tuple(shape), ns=ns, norb=norb, T=T,
                     ene_cutoff=ene_cutoff, fill=fill, shell_tol=shell_tol,
                     return_info=return_info)
    if return_info:
        return out[0], out[2]
    return out[0]


#: the ``[file.input]`` keys that name an initial state, not a model input
_INITIAL_KEYS = ("initial", "initial_uhf")


def referenced_input_files(toml_path):
    """The files the input ``toml_path`` names, as ``{TOML key: absolute
    path}``: every ``[file.input.interaction]`` file (Geometry, Transfer and
    each interaction, keyed as written) and the ``[file.input]`` files
    ``geometry_uhf`` and ``onebodyg_uhf``; not ``initial`` or
    ``initial_uhf`` (an initial state). Paths are resolved as
    :func:`load_problem` resolves them (from the TOML's directory). The
    TOML is read by the rules of :func:`load_problem`."""
    toml_path = pathlib.Path(os.path.abspath(os.fspath(toml_path)))
    doc = _read_config(toml_path)
    from requests.structures import CaseInsensitiveDict
    base_dir = str(toml_path.parent)
    finput = dict(doc.get("file", {}).get("input", {}))
    interaction = CaseInsensitiveDict(finput.get("interaction", {}))
    inter_dir = os.path.join(base_dir, interaction.get("path_to_input", ""))
    input_dir = os.path.join(base_dir, finput.get("path_to_input", ""))
    out = {}
    for key, name in interaction.items():
        if key.lower() != "path_to_input":
            out[key] = pathlib.Path(os.path.abspath(os.path.join(inter_dir, name)))
    for key in _INPUT_KEYS:
        if key in ("path_to_input", "interaction", "initial_mode") or key in _INITIAL_KEYS:
            continue
        if key in finput:
            out[key] = pathlib.Path(os.path.abspath(os.path.join(input_dir, finput[key])))
    return out
