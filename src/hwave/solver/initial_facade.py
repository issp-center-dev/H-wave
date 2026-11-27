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
* :func:`bands_of`: the bands of a Hamiltonian with the solver's
  occupations (:class:`Bands`).
* :func:`fold_hamiltonian_k`: a primitive-cell ``H(k)`` folded onto the
  k-mesh of a ``SubShape``.
* :func:`static_lindhard`, :func:`apply_lindhard`: the static Lindhard
  response of :class:`Bands` on pair-space labels, as a dense matrix and
  as a matrix-free product.
* :func:`run_para`: a UHFk self-consistent run of an input on the
  primitive cell (``SubShape = [1, 1, 1]``) at a given T from a given
  Green function (:class:`ParaRun`).
* :func:`referenced_input_files`: the files an input TOML names, by key.
* :func:`primitive_h0_k`: the bare ``T(k)`` of the original (primitive)
  cell on the ``CellShape`` mesh, built by the solver's
  :mod:`hamiltonian_k` (one Fourier convention).
* Re-exported pure helpers: :func:`to_k` (``hamiltonian_k.to_k``, the
  solver's R -> k transform ``T(k) = sum_R e^{+ikR} T(R)``),
  :func:`interleaved_index` (``spin_orbital.interleaved_index``, the
  spin-orbital index ``2*orb + spin``) and :func:`fold_layout` /
  :func:`unfold_layout` (``fold``, the sublattice layout UHFk reads and
  writes its Green function in), :func:`diag_blocks` and
  :func:`assemble_green` (``uhfk_fill``: the per-block eigen-decomposition
  of ``H(k)`` and the stored Green function of given eigenvectors and
  occupations, the two halves of :func:`green_function`).

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
from .fold import fold_layout, unfold_layout
from .hamiltonian_k import to_k
from .spin_orbital import interleaved_index
from .uhfk_blocks import HFTables, SolverBlocks, accumulate_hf_projected
from .uhfk_fill import (FillInfo, MuNotConvergedError, ShellRecord, assemble_green,
                        diag_blocks, fill_occupations)

__all__ = ["UHFkProblem", "SolverBlocks", "HFTables", "FillInfo", "ShellRecord",
           "MuNotConvergedError", "TWO_BODY_TYPES", "REAL_TWO_BODY_TYPES", "load_problem",
           "solver_blocks", "hf_tables", "accumulate_hf_projected",
           "hf_potential", "green_function", "referenced_input_files",
           "primitive_h0_k", "to_k", "interleaved_index", "fold_layout",
           "unfold_layout", "diag_blocks", "assemble_green", "Bands", "bands_of",
           "fold_hamiltonian_k", "LINDHARD_DEGENERACY_TOL", "LINDHARD_ZERO_T_GAP",
           "static_lindhard",
           "apply_lindhard", "ParaRun", "run_para", "PARA_EPS_MIN", "discard_profile_statistics"]

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


#: the two-body tables whose coefficients must be real. For them UHFk's
#: Hartree-Fock potential uses one coefficient ``J_ab(r)`` for both members
#: of a Hermitian pair of terms, so a complex coefficient gives a potential
#: that is not the mean field of the documented Hamiltonian (with or
#: without the Fock term) and a vertex that is not Hermitian. PairHop is
#: complex-safe and is not listed; CoulombIntra is on-site and real by its
#: Hermitian closure.
REAL_TWO_BODY_TYPES = ("Coulomb", "CoulombInter", "Hund", "Ising", "Exchange", "PairLift")


def _refuse_complex_couplings(toml_path, param_ham):
    """``ValueError`` naming the first complex coefficient of a table of
    :data:`REAL_TWO_BODY_TYPES` (term, bond and 1-based orbitals)."""
    for name in REAL_TWO_BODY_TYPES:
        table = param_ham.get(name)
        if not table:
            continue
        for (irvec, orbvec), v in table.items():
            if complex(v).imag != 0:
                raise ValueError(
                    "{}: {} has the complex coefficient {} at R = {}, orbitals ({}, {}) "
                    "(1-based). hwave.initial accepts complex coefficients only for "
                    "PairHop: for {} the solver's mean-field potential is not that of "
                    "the documented Hamiltonian when the coefficient is complex. Use "
                    "real coefficients.".format(
                        toml_path, name, complex(v), tuple(int(r) for r in irvec),
                        int(orbvec[0]) + 1, int(orbvec[1]) + 1, name))


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
        hint = ""
        if info_mode["mode"] == "UHFr":
            hint = ("; a UHFr input has no lattice geometry: give the UHFk input of the "
                    "same lattice, and write the UHFr initial file with "
                    "hwave.initial.writer.initial_state(..., site_map=, uhfr_toml=)")
        raise ValueError("{}: [mode] mode = {!r}: only UHFk inputs are supported{}"
                         .format(toml_path, info_mode["mode"], hint))
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


def _build_solver(toml_path, *, subshape=None, with_green=False, param=None,
                  prologue=True):
    """The UHFk instance of ``toml_path`` after the ``solve`` prologue
    (``T(k)``, the interaction tables and the blocks), as ``hwave`` would
    build it with its working directory at the TOML's directory.

    ``param`` overrides ``[mode.param]`` keys (case-insensitively, as the
    solver reads them); ``prologue=False`` stops after the constructor (the
    solver's ``solve`` runs the prologue itself).

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
    for key, value in (param or {}).items():
        info_mode["param"][key] = value

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
        if prologue:
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
        A non-UHFk mode, a refused key, an input the solver refuses, or a
        complex coefficient in a table of :data:`REAL_TWO_BODY_TYPES`.
    """
    solver, _, input_files = _build_solver(toml_path, subshape=subshape)
    _refuse_complex_couplings(toml_path, solver.param_ham_orig if solver.has_sublattice
                              else solver.param_ham)
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


# eq=False: the fields hold arrays (identity comparison and hash)
@dataclass(frozen=True, eq=False)
class Bands:
    """The bands of a Hamiltonian with the solver's occupations.

    ``eps_k`` ``(nk, nb)``: the eigenvalues, block by block (the solver's
    block order), ascending within a block. ``vec_k`` ``(nk, nd, nb)``: the
    eigenvectors embedded in the full index space (column ``l`` is the
    eigenvector of band ``l``, zero outside its block). ``kgrid``: the
    k-mesh ``(nx, ny, nz)`` (k in C order, ``k_i = 2 pi m_i / n_i``). ``T``:
    the temperature of the occupations. ``mu``: one chemical potential per
    mu-group. ``band_group`` ``(nk, nb)``: the mu-group of each band.
    ``occ_k`` ``(nk, nb)``: the occupation of each band, the solver's fill
    (:func:`hwave.solver.uhfk_fill.fill_occupations`). ``band_block``
    ``(nb,)``: the solver block of each band (its index in
    ``block_info``). ``ns``, ``norb``: the stored layout of the index space
    (``nd = ns * norb``; ``ns = 1`` in spin-orbital mode, where the index
    is ``2*orb + spin``).

    The stored Green function of the bands is ``G(k)_ab = sum_l
    conj(vec_k[k, a, l]) occ_k[k, l] vec_k[k, b, l]``, ``G(r)`` its
    ``fftn(..., norm="forward")``: :func:`green_function` of the same
    Hamiltonian. Every array is read-only. Not the ``Bands`` of
    :mod:`hwave.results`.
    """
    eps_k: np.ndarray
    vec_k: np.ndarray
    kgrid: Tuple[int, int, int]
    T: float
    mu: Tuple[float, ...]
    band_group: np.ndarray
    occ_k: np.ndarray
    band_block: np.ndarray
    ns: int
    norb: int


def bands_of(H_k, blocks, T, *, shape, ns, norb, ene_cutoff=None):
    """The :class:`Bands` of the Hamiltonian ``H_k`` ``(nvol, nd, nd)``
    (solver index order) with the blocks and mu-groups of ``blocks``
    (:class:`SolverBlocks`, ``group_nconds`` included) at temperature ``T``.

    It composes the solver's per-block diagonalisation
    (:func:`diag_blocks`) and occupation fill
    (:func:`hwave.solver.uhfk_fill.fill_occupations`, ``fill="solver"``),
    the two steps of :func:`green_function`, so ``mu`` and the occupations
    are those of the solver. ``shape`` is the k-grid, ``ns`` and ``norb``
    the stored layout, ``ene_cutoff`` the Fermi-function cut (``None``: the
    solvers' default). Raises :class:`MuNotConvergedError` when the T > 0
    ``mu`` search fails.
    """
    from .base import PARAM_DEFAULTS
    if ene_cutoff is None:
        ene_cutoff = PARAM_DEFAULTS["ene_cutoff"]
    H_k = np.asarray(H_k)
    shape = tuple(int(n) for n in shape)
    nvol = int(np.prod(shape))
    nd = int(ns) * int(norb)
    if H_k.shape != (nvol, nd, nd):
        raise ValueError("bands_of: H_k must be (nvol, nd, nd) = {}, got {}"
                         .format((nvol, nd, nd), H_k.shape))
    ws, vs = diag_blocks(H_k, blocks.block_info)
    dists, mus = fill_occupations(ws, vs, blocks.block_info, blocks.block_to_group,
                                  blocks.group_nconds, shape=shape, T=T,
                                  ene_cutoff=ene_cutoff)
    nb = sum(w.shape[1] for w in ws)
    eps = np.empty((nvol, nb))
    occ = np.empty((nvol, nb))
    vec = np.zeros((nvol, nd, nb), dtype=np.complex128)
    group = np.empty((nvol, nb), dtype=int)
    block = np.empty(nb, dtype=np.int64)
    col = 0
    for kb, (blk, g, w, v, d) in enumerate(zip(blocks.block_info, blocks.block_to_group,
                                               ws, vs, dists)):
        n = w.shape[1]
        idx = np.asarray(blk)
        eps[:, col:col + n] = w
        occ[:, col:col + n] = d
        vec[:, idx, col:col + n] = v
        group[:, col:col + n] = int(g)
        block[col:col + n] = kb
        col += n
    return Bands(eps_k=_read_only(eps), vec_k=_read_only(vec), kgrid=shape, T=float(T),
                 mu=tuple(float(m) for m in mus), band_group=_read_only(group),
                 occ_k=_read_only(occ), band_block=_read_only(block), ns=int(ns),
                 norb=int(norb))


def fold_hamiltonian_k(H_k, *, cellshape, subshape, norb_orig, ns):
    """A primitive-cell Hamiltonian ``H_k`` ``(prod(cellshape), nd_orig,
    nd_orig)`` on the ``CellShape`` mesh, folded onto the ``SubShape``
    cell: ``(nvol, nd, nd)`` on the mesh ``cellshape // subshape``, in the
    folded index order UHFk uses (folded orbital ``a + norb_orig * ir``).

    The fold is done in real space with the solver's Fourier convention:
    ``H(R)`` by :func:`hamiltonian_k.from_k` (the inverse of
    :func:`to_k`), :func:`fold_layout` in the Hamiltonian convention
    (``hamiltonian=True``, the layout ``fold.reshape_interaction`` gives the
    Transfer table), then :func:`to_k` on the supercell mesh. Callers never
    apply :func:`fold_layout` to a k-space array.
    """
    cellshape = tuple(int(n) for n in cellshape)
    subshape = tuple(int(n) for n in subshape)
    norb_orig, ns = int(norb_orig), int(ns)
    nd_orig = ns * norb_orig
    H_k = np.asarray(H_k)
    lvol = int(np.prod(cellshape))
    if H_k.shape != (lvol, nd_orig, nd_orig):
        raise ValueError("fold_hamiltonian_k: H_k must be (prod(cellshape), nd_orig, "
                         "nd_orig) = {}, got {}".format((lvol, nd_orig, nd_orig), H_k.shape))
    h_r = hamiltonian_k.from_k(H_k, cellshape).reshape(lvol, ns, norb_orig, ns, norb_orig)
    sub = fold_layout(h_r, cellshape=cellshape, subshape=subshape, norb_orig=norb_orig,
                      ns=ns, hamiltonian=True)
    shape = tuple(c // b for c, b in zip(cellshape, subshape))
    nd = sub.shape[1] * sub.shape[2]
    return to_k(sub.reshape(shape + (nd, nd)))


#: ``|eps_a - eps_b|`` below which a pair of levels takes the derivative
#: limit ``-f'`` in :func:`static_lindhard` (an absolute energy)
LINDHARD_DEGENERACY_TOL = 1e-10
#: at T = 0, a pair of levels of one block closer than this (an absolute
#: energy) with different occupations is refused like an open Fermi shell:
#: its weight 1/|eps_a - eps_b| would exceed 1e6 and grow without bound as
#: the splitting goes to zero
LINDHARD_ZERO_T_GAP = 1e-6
#: :func:`apply_lindhard` works on chunks of k points holding at most this
#: many complex entries per array (its memory stays independent of |P|^2)
_LINDHARD_CHUNK_ENTRIES = 2 ** 17
#: :func:`static_lindhard` builds its factor ``W`` on chunks of k points of
#: at most this many complex entries (16 MB) per array; a chunk holds
#: ``W``, ``F W`` and their reshaped copies at once
_LINDHARD_DENSE_CHUNK_ENTRIES = 2 ** 20


def _lindhard_kernel(ea, oa, eb, ob, T):
    """The static Lindhard weight of a pair of levels with energies ``ea``,
    ``eb`` and occupations ``oa``, ``ob`` (arrays of one shape):
    ``(ob - oa) / (ea - eb)``, and where ``|ea - eb| <``
    :data:`LINDHARD_DEGENERACY_TOL` the derivative limit ``-f'``: at T > 0
    ``(oa (1 - oa) + ob (1 - ob)) / (2 T)`` (``-f'(eps - mu) = f (1 - f) /
    T`` of each level, averaged), at T = 0 zero when both occupations are
    the same 0 or 1 and ``nan`` otherwise (a degenerate pair at the Fermi
    level: the response is unbounded). At T = 0 a pair closer than
    :data:`LINDHARD_ZERO_T_GAP` with different occupations is also ``nan``
    (a nearly open shell). Nonnegative for levels of one mu-group (the
    occupations fall with the energy)."""
    ea, oa, eb, ob = np.broadcast_arrays(*(np.asarray(x, dtype=float)
                                           for x in (ea, oa, eb, ob)))
    d = ea - eb
    deg = np.abs(d) < LINDHARD_DEGENERACY_TOL
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(deg, 0.0, (ob - oa) / np.where(deg, 1.0, d))
    if T > 0:
        lim = (oa * (1.0 - oa) + ob * (1.0 - ob)) / (2.0 * T)
    else:
        closed = (oa == ob) & ((oa == 0.0) | (oa == 1.0))
        lim = np.where(closed, 0.0, np.nan)
        out = np.where((np.abs(d) < LINDHARD_ZERO_T_GAP) & (oa != ob), np.nan, out)
    return np.where(deg, lim, out)


def _q_mesh(q, shape, who):
    """The mesh index ``(3,)`` of ``q`` (three rationals in reduced lattice
    units: ints, Fractions, strings such as ``"1/3"``, or floats within
    1e-9 of the mesh)."""
    from fractions import Fraction
    q = tuple(q)
    if len(q) != 3:
        raise ValueError("{}: q must have three components, got {!r}".format(who, q))
    m = []
    for i, (x, n) in enumerate(zip(q, shape)):
        if isinstance(x, (float, np.floating)):
            y = float(x) * n
            mi = int(round(y))
            ok = np.isfinite(y) and abs(y - mi) <= 1e-9
        else:
            y = Fraction(x) * n
            ok = y.denominator == 1
            mi = int(y) if ok else 0
        if not ok:
            raise ValueError("{}: q = {} is not on the k-mesh {} (component {} times {} "
                             "is not an integer)".format(who, ",".join(str(v) for v in q),
                                                         list(shape), i, n))
        m.append(mi % int(n))
    return np.array(m, dtype=np.int64)


def _lindhard_setup(bands, q, labels, who):
    """``(alpha, beta, r, kq_idx, kq, F)``: the solver indices of the two
    legs of each label, its bond, the mesh index of ``k + q`` per k (and
    its flat index) and the Lindhard weights ``F[k, l, m]`` (band ``l`` at
    ``k + q``, band ``m`` at ``k``; zero across solver blocks)."""
    if not isinstance(bands, Bands):
        raise ValueError("{}: bands must be a Bands (bands_of), got {!r}"
                         .format(who, type(bands).__name__))
    labels = np.asarray(labels)
    if labels.ndim != 2 or labels.shape[1] != 7 or not np.issubdtype(labels.dtype,
                                                                       np.integer):
        raise ValueError("{}: labels must be an int array (|P|, 7) = (a, s, b, t, rx, ry, "
                         "rz), got {} {}".format(who, labels.dtype, labels.shape))
    labels = labels.astype(np.int64)
    ns, norb = int(bands.ns), int(bands.norb)
    so = ns == 1
    nphys = norb // 2 if so else norb
    a, s, b, t = (labels[:, i] for i in range(4))
    if np.any((a < 0) | (a >= nphys) | (b < 0) | (b >= nphys)):
        raise ValueError("{}: a label orbital is outside 0..{} (the physical orbitals of "
                         "the bands)".format(who, nphys - 1))
    if np.any((s < 0) | (s > 1) | (t < 0) | (t > 1)):
        raise ValueError("{}: a label spin is not 0 or 1".format(who))
    if so:
        alpha, beta = interleaved_index(a, s), interleaved_index(b, t)
    else:
        alpha, beta = s * norb + a, t * norb + b
    alpha = np.asarray(alpha, dtype=np.int64)
    beta = np.asarray(beta, dtype=np.int64)
    iblk = _index_blocks(bands)
    cross = iblk[alpha] != iblk[beta]
    if np.any(cross):
        m = int(np.nonzero(cross)[0][0])
        raise ValueError("{}: label {} joins two solver blocks of the bands; the bands "
                         "have another block structure than the labels (for example a "
                         "fixed 2Sz against labels built without it): use the bands of the "
                         "same model".format(who, [int(v) for v in labels[m]]))
    shape = tuple(int(n) for n in bands.kgrid)
    nk = int(np.prod(shape))
    mq = _q_mesh(q, shape, who)
    kidx = np.stack(np.unravel_index(np.arange(nk), shape), axis=1)
    kq_idx = (kidx + mq[None, :]) % np.array(shape)
    kq = np.ravel_multi_index(tuple(kq_idx.T), shape)
    eps, occ = bands.eps_k, bands.occ_k
    F = _lindhard_kernel(eps[kq][:, :, None], occ[kq][:, :, None], eps[:, None, :],
                         occ[:, None, :], float(bands.T))
    blk = np.asarray(bands.band_block)
    F[:, blk[:, None] != blk[None, :]] = 0.0
    flagged = np.argwhere(np.isnan(F))
    vec = bands.vec_k
    for i0 in range(0, flagged.shape[0], 256):
        fk, fl, fm = flagged[i0:i0 + 256].T
        # the weight of the pair on the labels: sum_m |v_alpha,l(k+q)|^2 |v_beta,m(k)|^2
        w = np.sum(np.abs(vec[kq[fk], :, fl][:, alpha]) ** 2
                   * np.abs(vec[fk, :, fm][:, beta]) ** 2, axis=1)
        hit = w > 1e-20
        if np.any(hit):
            k = int(fk[np.nonzero(hit)[0][0]])
            raise ValueError("{}: the bands at T = 0 have an open Fermi shell on the labels "
                             "(levels of one block within {:g} of each other with different "
                             "occupations at k index {} and k + q); the static response is "
                             "unbounded there: use T > 0 (T_hint)".format(
                                 who, LINDHARD_ZERO_T_GAP, list(kidx[k])))
        F[fk, fl, fm] = 0.0
    return alpha, beta, labels[:, 4:7], kq_idx, kq, F


def _index_blocks(bands):
    """The solver block of each index of the bands' index space (from the
    support of the embedded eigenvectors; ``-1`` for an index no band
    holds)."""
    weight = np.sum(np.abs(np.asarray(bands.vec_k)) ** 2, axis=0)      # (nd, nb)
    blk = np.asarray(bands.band_block)
    out = blk[np.argmax(weight, axis=1)]
    return np.where(weight.max(axis=1) > 0, out, -1)


def _lindhard_phases(kq_idx, r, shape):
    """``e^{i(k+q).r}`` ``(len(kq_idx), |P|)``, reduced exactly on the mesh
    before the exponential."""
    L = int(np.lcm.reduce([int(n) for n in shape]))
    num = np.zeros((kq_idx.shape[0], r.shape[0]), dtype=np.int64)
    for i, n in enumerate(shape):
        num += np.outer(kq_idx[:, i], (r[:, i] % int(n)) * (L // int(n)))
    return np.exp(2j * np.pi * (num % L) / L)


def static_lindhard(bands, q, labels):
    """The static Lindhard response ``chi0(q)`` on the pair-space labels,
    dense ``(|P|, |P|)`` complex128, Hermitian and positive semidefinite.

    ``bands``: :class:`Bands` (:func:`bands_of`, the para Hartree-Fock bands
    of the primitive cell at ``bands.T``). ``labels``: ``(|P|, 7)`` int
    ``(a, s, b, t, rx, ry, rz)``, ``a``, ``b`` physical orbitals (0-based)
    and ``s``, ``t`` spins (0 = up) in every mode, ``r`` a bond in
    primitive cells (any representative modulo the mesh). ``q``: three
    rationals in reduced lattice units on the k-mesh of the bands.

    **Coordinates.** A label ``m`` is the slot of the one-body operator
    ``sum_R c^dag_{as}(R + r) c_{bt}(R)``. A density vector holds the
    rho-components ``x_m = <c^dag_{bt}(R) c_{as}(R + r)>`` (``rho = G^T``),
    a field the coefficients ``h_m`` of ``delta H = sum_R sum_m h_m e^{iq.R}
    c^dag_{as}(R + r) c_{bt}(R)``; a q-modulated vector carries ``e^{iq.R}``
    at the cell of the annihilation operator. To linear order ``delta x =
    chi0 (-h)``:

    ``chi0[m, n] = (1/N) sum_k sum_{l, l'} w_m(k, l, l') conj(w_n(k, l,
    l')) F(eps_l(k+q), eps_l'(k))``, ``w_m(k, l, l') = e^{i(k+q).r_m}
    v_{alpha_m, l}(k+q) conj(v_{beta_m, l'}(k))``

    (``alpha``, ``beta`` the solver indices of ``(a, s)`` and ``(b, t)``,
    ``v`` the eigenvectors of :class:`Bands`, all bands), with
    ``F = (f_l' - f_l) / (eps_l - eps_l')`` from the occupations
    ``occ_k`` and, for ``|eps_l - eps_l'| <``
    :data:`LINDHARD_DEGENERACY_TOL`, the derivative limit ``-f'(eps -
    mu_group) = f (1 - f) / T``. Bands of different solver blocks do not
    pair, and a label whose two legs lie in different blocks of the bands
    is refused (``ValueError``). At ``T = 0`` a degenerate pair contributes
    zero when both levels are filled or both empty; a pair of one block
    closer than :data:`LINDHARD_ZERO_T_GAP` with different occupations (an
    open or nearly open Fermi shell) is refused (``ValueError``) when it
    has weight on the labels: the response is unbounded. Such a pair on
    orbitals outside the labels does not matter.

    Cost ``O(N_k nb^2 |P|^2)``; :func:`apply_lindhard` gives the product
    with a vector without forming the matrix.
    """
    alpha, beta, r, kq_idx, kq, F = _lindhard_setup(bands, q, labels, "static_lindhard")
    shape = tuple(int(n) for n in bands.kgrid)
    vec = bands.vec_k
    nk, _nd, nb = vec.shape
    npair = alpha.size
    chi = np.zeros((npair, npair), dtype=np.complex128)
    if npair == 0:
        return chi
    chunk = max(1, _LINDHARD_DENSE_CHUNK_ENTRIES // max(1, npair * nb * nb))
    for k0 in range(0, nk, chunk):
        ks = np.arange(k0, min(nk, k0 + chunk))
        ph = _lindhard_phases(kq_idx[ks], r, shape)                 # (c, P)
        A = vec[kq[ks]][:, alpha, :]                                # (c, P, nb)
        B = np.conj(vec[ks][:, beta, :])                            # (c, P, nb)
        W = ph[:, :, None, None] * A[:, :, :, None] * B[:, :, None, :]
        Wf = W * F[ks][:, None, :, :]
        W = np.moveaxis(W, 1, 0).reshape(npair, -1)
        Wf = np.moveaxis(Wf, 1, 0).reshape(npair, -1)
        chi += Wf @ W.conj().T
    chi /= nk
    return (chi + chi.conj().T) / 2.0


def apply_lindhard(bands, q, labels, v):
    """``chi0(q) v`` of :func:`static_lindhard` without forming ``chi0``:
    ``v`` ``(|P|,)`` or ``(|P|, c)``, the result of the same shape
    (complex128).

    Per k: the field matrix ``M(k) = sum_m e^{-i(k+q).r_m} v_m
    |alpha_m><beta_m|``, its band form ``V(k+q)^dag M(k) V(k)`` times
    ``F``, back to orbitals, and the label entries times ``e^{i(k+q).r_m}``:
    ``O(N_k (nb^3 + |P|))`` time, and memory ``O(N_k nb^2)`` plus chunks of
    k bounded independently of ``|P|^2``."""
    import scipy.sparse
    # the setup (F, the label maps) is rebuilt at every call; a per-q cache
    # for the repeated products of an iterative eigensolver is left to the
    # medium tier
    alpha, beta, r, kq_idx, kq, F = _lindhard_setup(bands, q, labels, "apply_lindhard")
    shape = tuple(int(n) for n in bands.kgrid)
    vec = bands.vec_k
    nk, nd, nb = vec.shape
    npair = alpha.size
    v = np.asarray(v)
    if v.ndim not in (1, 2) or v.shape[0] != npair:
        raise ValueError("apply_lindhard: v must be (|P|,) or (|P|, c) with |P| = {}, got "
                         "{}".format(npair, v.shape))
    one = v.ndim == 1
    if npair == 0:
        return np.zeros(v.shape, dtype=np.complex128)
    V = v.reshape(npair, -1).astype(np.complex128)
    ncol = V.shape[1]
    S = scipy.sparse.csr_array((np.ones(npair), (alpha * nd + beta, np.arange(npair))),
                               shape=(nd * nd, npair))
    y = np.zeros((npair, ncol), dtype=np.complex128)
    per_k = max(npair * ncol, ncol * nd * nd, nb * nb)
    chunk = max(1, _LINDHARD_CHUNK_ENTRIES // per_k)
    for k0 in range(0, nk, chunk):
        ks = np.arange(k0, min(nk, k0 + chunk))
        c = ks.size
        ph = _lindhard_phases(kq_idx[ks], r, shape)                 # (c, P)
        Mv = np.conj(ph).T[:, :, None] * V[:, None, :]              # (P, c, ncol)
        M = (S @ Mv.reshape(npair, c * ncol)).reshape(nd, nd, c, ncol)
        M = np.moveaxis(M, (2, 3), (0, 1))                          # (c, ncol, nd, nd)
        Ukq = vec[kq[ks]][:, None]                                  # (c, 1, nd, nb)
        Uk = vec[ks][:, None]
        X = np.conj(np.swapaxes(Ukq, -1, -2)) @ M @ Uk               # (c, ncol, nb, nb)
        X *= F[ks][:, None]
        Y = Ukq @ X @ np.conj(np.swapaxes(Uk, -1, -2))              # (c, ncol, nd, nd)
        y += np.einsum("cm,cjm->mj", ph, Y[:, :, alpha, beta])
    y /= nk
    return y[:, 0] if one else y


# eq=False: the fields hold arrays
@dataclass(frozen=True, eq=False)
class ParaRun:
    """The result of :func:`run_para`.

    ``G``: the stored Green function at the end of the run ``(nvol, ns,
    norb, ns, norb)`` (the solver's mixed G, the ``green`` it writes).
    ``H_k``: ``T(k) + hf_potential(G)`` ``(nvol, nd, nd)``, the
    Hamiltonian of ``G`` (block-projected). ``free_energy``: the solver's
    ``FreeEnergy_Total`` (of the whole lattice). ``residual``: the last SCF
    residual (``nan`` after zero steps). ``converged``, ``iterations``: the
    solver's convergence record. ``mu``: one chemical potential per
    mu-group of the last step. ``work_dir``: the directory the run wrote
    its results to (``energy.dat``, ``green.npz``, ``scf_status.dat``, and
    the start Green function ``start.npz``). ``eps``: the ``EPS`` exponent
    the run used (:data:`PARA_EPS_MIN` or the input's, whichever is larger);
    ``iteration_max``: its ``IterationMax`` (the input's).
    """
    G: np.ndarray
    H_k: np.ndarray
    free_energy: float
    residual: float
    converged: bool
    iterations: int
    mu: Tuple[float, ...]
    work_dir: pathlib.Path
    eps: int
    iteration_max: int


#: the smallest ``EPS`` exponent of :func:`run_para` (SCF residual below
#: 1e-12): the para free energies are compared to 1e-8 per cell, which a
#: looser SCF does not resolve
PARA_EPS_MIN = 12


def _input_eps(toml_path):
    """The ``[mode.param] EPS`` exponent of the input (case-insensitive, the
    solvers' default when absent)."""
    import math
    import numbers
    from .base import PARAM_DEFAULTS, ParamDict
    param = ParamDict(_read_config(toml_path)["mode"].get("param", {}) or {})
    eps = param.get("EPS", PARAM_DEFAULTS["EPS"])
    if isinstance(eps, bool) or not isinstance(eps, numbers.Real) \
            or not math.isfinite(float(eps)):
        raise ValueError("{}: [mode.param] EPS = {!r} is not a finite number"
                         .format(toml_path, eps))
    return eps


def run_para(toml_path, *, T, start_G, work_dir):
    """Run UHFk for the input ``toml_path`` on its primitive cell.

    The problem is loaded as :func:`load_problem` loads it, with ``SubShape
    = [1, 1, 1]`` (whatever the input's ``SubShape``), the temperature
    ``T`` and the convergence exponent ``EPS = max(EPS_input,
    PARA_EPS_MIN)`` (an internal reference run, converged tighter than the
    candidate runs); ``IterationMax``, ``Mix`` and every other parameter are
    the input's. The start Green function ``start_G`` (stored layout of
    the primitive cell, ``(prod(CellShape), ns, norb_orig, ns,
    norb_orig)``) is written to ``work_dir/start.npz`` and read by the
    solver as its ``initial`` Green function; the solver's public ``solve``
    then runs in ``work_dir`` and its ``save_results`` writes
    ``energy.dat``, ``green.npz`` and ``scf_status.dat`` there.

    Returns
    -------
    ParaRun

    Raises
    ------
    ValueError
        A refused input, a ``start_G`` of the wrong shape, a negative T,
        or a solver exit (for example a failed chemical-potential search).
    """
    T = float(T)
    if not np.isfinite(T) or T < 0.0:
        raise ValueError("run_para: T must be a finite number >= 0, got {!r}".format(T))
    work_dir = pathlib.Path(os.path.abspath(os.fspath(work_dir)))
    work_dir.mkdir(parents=True, exist_ok=True)
    eps_in = _input_eps(toml_path)
    eps = int(eps_in) if eps_in >= PARA_EPS_MIN else PARA_EPS_MIN
    solver, _, _ = _build_solver(toml_path, subshape=(1, 1, 1),
                                 param={"T": T, "EPS": eps}, prologue=False)
    want = (solver.nvol, solver.ns, solver.norb, solver.ns, solver.norb)
    start_G = np.asarray(start_G, dtype=np.complex128)
    if start_G.shape != want:
        raise ValueError("run_para: start_G must be (nvol, ns, norb, ns, norb) = {} on "
                         "the primitive cell, got {}".format(want, start_G.shape))
    start = work_dir / "start.npz"
    np.savez(start, green=start_G)
    info_out = {"path_to_output": str(work_dir), "energy": "energy.dat",
                "green": "green.npz"}
    try:
        with np.load(start, allow_pickle=False) as data:
            green_info = {"initial": data}
            solver.solve(green_info, str(work_dir))
        solver.save_results(info_out, {})
    except SystemExit as exc:
        raise ValueError("{}: the paramagnetic UHFk run stopped (exit status {}); the "
                         "logged error says why.".format(toml_path, exc.code)) from None
    G = np.array(solver.Green, dtype=np.complex128)
    H_k = np.array(solver.ham_trans, dtype=np.complex128)
    accumulate_hf_projected(H_k, G, solver._hf_tables(), include_fock=solver.iflag_fock,
                            blocks=solver.block_info)
    rest = solver.physics.get("Rest", float("nan")) if solver.scf_iterations else float("nan")
    free = solver.physics.get("FreeEne", {}).get("Total", float("nan"))
    mu = np.atleast_1d(solver._green_list.get("mu", ()))
    return ParaRun(G=_read_only(G), H_k=_read_only(H_k), free_energy=float(np.real(free)),
                   residual=float(rest), converged=bool(solver.scf_converged),
                   iterations=int(solver.scf_iterations),
                   mu=tuple(float(m) for m in mu), work_dir=work_dir, eps=eps,
                   iteration_max=int(solver.param_mod["IterationMax"]))


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


def discard_profile_statistics():
    """Drop the timing statistics of the solver's ``@do_profile`` functions,
    so that the table the solver prints when the interpreter exits is not
    printed (for a command that only reads models, such as ``python -m
    hwave.initial``). The solver's own runs are not affected."""
    from hwave.solver import perf
    perf._perf_db_data._db_count.clear()
    perf._perf_db_data._db_value.clear()
