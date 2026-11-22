"""Readers of the UHFr and UHFk outputs.

* :func:`read_energy`: ``energy.dat`` (:class:`Energy`).
* :func:`read_uhfr_eigen`: the UHFr eigenvalues and eigenvectors
  (:class:`UHFrEigen`).
* :func:`read_uhfr_green`: the UHFr Green function text (:class:`UHFrGreen`).
* :func:`read_uhfr_order_parameters`: the order parameters of a UHFr run
  (:class:`hwave.results.OrderParameters`).
* :func:`read_bands`: the UHFk eigenvalues on the k grid (:class:`Bands`).
* :func:`read_uhfk_green`: the UHFk ``green.npz`` (:class:`GreenR`).
* :func:`read_uhfk_order_parameters`: the order parameters of a UHFk run.
* :func:`so_density_from_green`: the spin density of a spin-orbital UHFk
  Green function.
* :func:`read_onebodyg`: refuses the UHFk ``onebodyg`` text (issue #357).

UHFr numbers its spin-orbitals ``a = site + Nsite * spin`` (``spin`` 0 for
up, 1 for down; uhfr.py:1308-1309), and its Green function is
``G[a, b] = <c^dag_a c_b>``.
"""
from __future__ import annotations

import dataclasses
import math
import numbers
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Optional, Union

import numpy as np

from hwave.results._frozen import _freeze
from hwave.results._io import (
    check_markers, input_file_path, member_shapes, npz_member, open_npz,
    read_geometry, read_key_value_text)
from hwave.results.errors import (
    MissingOutputError, UnsupportedLayoutError)
from hwave.results.order import order_parameters_from_density
from hwave.results.run import as_run
from hwave.results.status import require_scf_steps

__all__ = ["Energy", "read_energy", "UHFrEigen", "read_uhfr_eigen",
           "UHFrGreen", "read_uhfr_green", "read_uhfr_order_parameters",
           "Bands", "read_bands", "GreenR", "read_uhfk_green",
           "read_uhfk_order_parameters", "OneBodyG", "read_onebodyg",
           "so_density_from_green"]

#: the ``energy.dat`` keys of the scalar fields of :class:`Energy`, in lower
#: case: UHFr writes ``Energy_total`` (uhfr.py:1271-1292), UHFk
#: ``Energy_Total`` (uhfk.py:1796-1835)
_ENERGY_FIELDS = (
    ("total", "energy_total"),
    ("band", "energy_band"),
    ("free_total", "freeenergy_total"),
    ("free_band", "freeenergy_band"),
    ("ncond", "ncond"),
    ("sz", "sz"),
    ("sx", "sx"),
    ("sy", "sy"),
)

_MU = "chemicalpotential"

#: the order of :attr:`Energy.meta` ``["missing_keys"]``
_FIELD_ORDER = ("total", "band", "free_total", "free_band", "ncond", "sz", "mu",
                "sx", "sy")


@dataclasses.dataclass(frozen=True)
class Energy:
    """The content of ``energy.dat``.

    Attributes
    ----------
    values : dict[str, float]
        Every line of the file, ``{key: value}`` in file order, the keys
        stripped of their padding (e.g. ``"Energy_CoulombIntra"``).
    total, band : float or None
        ``Energy_Total``/``Energy_total`` and ``Energy_Band``/``Energy_band``.
    free_total, free_band : float or None
        The free energy ``F = E - TS`` (``FreeEnergy_*``).
    ncond, sz : float or None
        ``NCond`` and ``Sz``.
    sx, sy : float or None
        ``Sx`` and ``Sy``, the last two lines of the file; ``None`` in files
        written before they were added.
    mu : float, dict[str, float] or None
        ``ChemicalPotential``; a dict ``{suffix: value}`` when the file has
        one ``ChemicalPotential_<suffix>`` per group (UHFk mu groups ``0``,
        ``1``, ...; UHFr blocks ``spin-up``, ``spin-down``).
    solver : str
        ``"UHFr"`` or ``"UHFk"``, from the input.
    source : Path
        The file read.
    meta : dict
        ``missing_keys``: the names of the fields above whose key is absent
        from the file, in field order.

    The keys are matched case-insensitively (the solvers write them in
    different case).
    """
    values: dict
    total: Optional[float]
    band: Optional[float]
    free_total: Optional[float]
    free_band: Optional[float]
    ncond: Optional[float]
    sz: Optional[float]
    sx: Optional[float]
    sy: Optional[float]
    mu: Union[float, dict, None]
    solver: str
    source: Path
    meta: dict


def read_energy(run_or_toml, *, base_dir=None) -> Energy:
    """Read ``energy.dat`` of a run (a :class:`hwave.results.Run` or the
    path of its input TOML file).

    Raises MissingOutputError when the output is not configured or its file
    is absent, or when the status file records zero SCF steps (the file is
    then left over from an earlier run; :func:`hwave.results.status.require_scf_steps`),
    and UnsupportedLayoutError for a line that is not
    ``key = number``, for a key given twice (compared case-insensitively),
    or for a file with both ``ChemicalPotential`` and
    ``ChemicalPotential_<suffix>``."""
    run = as_run(run_or_toml, base_dir=base_dir)
    require_scf_steps(run, "energy")
    path = run.path("energy")
    values = {}
    seen = {}
    for key, text, lineno in read_key_value_text(path):
        lower = key.lower()
        if lower in seen:
            raise UnsupportedLayoutError(
                "{}: {!r}: line {} repeats the key of line {}".format(
                    path, key, lineno, seen[lower]), path=path)
        seen[lower] = lineno
        try:
            values[key] = float(text)
        except ValueError:
            raise UnsupportedLayoutError(
                "{}: {!r}: line {} has the value {!r}, not a number".format(
                    path, key, lineno, text), path=path) from None

    by_lower = {key.lower(): value for key, value in values.items()}
    fields = {name: by_lower.get(lower) for name, lower in _ENERGY_FIELDS}
    groups = {key[len(_MU) + 1:]: value for key, value in values.items()
              if key.lower().startswith(_MU + "_")}
    if _MU in by_lower and groups:
        first = next(key for key in values if key.lower().startswith(_MU + "_"))
        raise UnsupportedLayoutError(
            "{}: {!r}: the file has both ChemicalPotential and per-group "
            "ChemicalPotential_<suffix> lines".format(path, first), path=path)
    mu = by_lower.get(_MU, groups or None)
    fields["mu"] = mu
    missing = tuple(name for name in _FIELD_ORDER if fields[name] is None)
    return Energy(values=_freeze(values), mu=_freeze(mu), solver=run.mode,
                  source=Path(path), meta={"missing_keys": missing},
                  **{name: fields[name] for name, _ in _ENERGY_FIELDS})


# ===================================================================== UHFr

#: the eigenvector blocks of UHFr (uhfr.py:553-560): one ``sz-free`` block
#: of all ``2 Nsite`` spin-orbitals, or one block per spin when 2Sz is fixed
_SZ_FREE_BLOCKS = ("sz-free",)
_FIXED_SZ_BLOCKS = ("spin-up", "spin-down")

#: at most this many missing Green function entries are named in a message
_MAX_LISTED = 10

_INT = re.compile(r"[+-]?[0-9]+")


def _uhfr_run(run_or_toml, base_dir):
    run = as_run(run_or_toml, base_dir=base_dir)
    if run.mode != "UHFr":
        raise UnsupportedLayoutError(
            "{} is a {!r} run, not a UHFr run".format(run.input_path, run.mode),
            path=run.input_path)
    return run


def _nsite(run):
    """``[mode.param] Nsite``, which the solver requires (uhfr.py:526, 547)."""
    value = run.param("Nsite", None)
    if isinstance(value, bool) or not isinstance(value, numbers.Integral) \
            or value < 1:
        raise UnsupportedLayoutError(
            "[mode.param] Nsite must be a positive integer in {}, got {!r}"
            .format(run.input_path, value), path=run.input_path)
    return int(value)


def _two_sz(run):
    """``[mode.param] 2Sz``: ``None`` (unset or ``"free"``) or an integer."""
    value = run.param("2Sz")
    if value is not None and (isinstance(value, bool)
                              or not isinstance(value, numbers.Integral)):
        raise UnsupportedLayoutError(
            "[mode.param] 2Sz must be an integer or \"free\" in {}, got {!r}"
            .format(run.input_path, value), path=run.input_path)
    return value


def _onebodyg_name(run):
    """The ``[file.input] OneBodyG`` file name, or ``None``. UHFr matches the
    key in any case, the last match winning (qlms.py:155-159)."""
    info_file = run.config.get("file", {})
    info_input = info_file.get("input", {}) if isinstance(info_file, Mapping) else {}
    name = None
    if isinstance(info_input, Mapping):
        for key, value in info_input.items():
            if key.lower() == "onebodyg":
                name = value
    return name


def _uhfr_output(run, key):
    """``run.path(key)``, after :func:`require_scf_steps` (outputs of a
    zero-step run are refused). When the file is configured but absent, the
    hint says why the solver may not have written it."""
    require_scf_steps(run, key)
    try:
        return run.path(key)
    except MissingOutputError as exc:
        if key not in run.outputs:
            raise
        hint = exc.hint
        # a zero-step run was refused by require_scf_steps above
        if key == "green" and _onebodyg_name(run) is None:
            hint = ("UHFr writes the Green function only when [file.input] "
                    "OneBodyG lists the entries to write (uhfr.py:1302-1314); "
                    "configure it and rerun")
        raise MissingOutputError(exc.message, path=exc.path, hint=hint) from None


@dataclasses.dataclass(frozen=True)
class UHFrEigen:
    """The eigenvalues (and eigenvectors) of the UHFr mean-field
    Hamiltonian, per block.

    Attributes
    ----------
    energies : dict[str, ndarray]
        ``{block: eigenvalues}`` in ascending order. The blocks are
        ``"sz-free"`` (``2 Nsite`` values) when 2Sz is not fixed, else
        ``"spin-up"`` and ``"spin-down"`` (``Nsite`` values each)
        (uhfr.py:553-560, 1296-1300).
    eigenvectors : dict[str, ndarray] or None
        ``{block: vectors}``, column ``k`` the eigenvector of
        ``energies[block][k]``, in the block's own basis: the
        spin-orbitals ``site + Nsite * spin`` for ``"sz-free"``, the sites
        for ``"spin-up"``/``"spin-down"``. ``None`` unless read with
        ``with_vectors=True``.
    with_vectors : bool
    sources : tuple of Path
        The files read, one per block, in block order.
    meta : dict
        ``nsite``, ``two_sz``.
    """
    energies: dict
    eigenvectors: Optional[dict]
    with_vectors: bool
    sources: tuple
    meta: dict


def read_uhfr_eigen(run_or_toml, *, base_dir=None, with_vectors=False) -> UHFrEigen:
    """Read the ``<block>_<eigen>.npz`` files of a UHFr run (a
    :class:`hwave.results.Run` or the path of its input TOML file).

    Raises MissingOutputError when ``[file.output] eigen`` is not configured
    or a block's file is absent, or when the status file records zero SCF
    steps (:func:`hwave.results.status.require_scf_steps`), and
    UnsupportedLayoutError for a run that
    is not UHFr, or for a file without the ``eigenvalue`` and
    ``eigenvector`` members of the expected shapes."""
    run = _uhfr_run(run_or_toml, base_dir)
    nsite = _nsite(run)
    two_sz = _two_sz(run)
    blocks = _SZ_FREE_BLOCKS if two_sz is None else _FIXED_SZ_BLOCKS
    size = 2 * nsite if two_sz is None else nsite
    energies, vectors, sources = {}, {}, []
    for block in blocks:
        path = _uhfr_output(run, "eigen:" + block)
        shapes = member_shapes(path)
        for name, shape in (("eigenvalue", (size,)), ("eigenvector", (size, size))):
            if name not in shapes:
                raise UnsupportedLayoutError(
                    "{}: no member {!r}".format(path, name), path=path)
            if shapes[name] != shape:
                raise UnsupportedLayoutError(
                    "{}: member {!r} has shape {}, expected {} for block {!r} "
                    "with Nsite = {}".format(path, name, shapes[name], shape,
                                             block, nsite), path=path)
        with open_npz(path) as data:
            energies[block] = npz_member(data, path, "eigenvalue", "f")
            if with_vectors:
                vectors[block] = npz_member(data, path, "eigenvector", "c")
        sources.append(Path(path))
    return UHFrEigen(energies=_freeze(energies),
                     eigenvectors=_freeze(vectors) if with_vectors else None,
                     with_vectors=bool(with_vectors), sources=tuple(sources),
                     meta={"nsite": nsite, "two_sz": two_sz})


@dataclasses.dataclass(frozen=True)
class UHFrGreen:
    """The UHFr Green function entries written to the ``green`` file.

    UHFr writes one line ``i s j t Re Im`` per entry of the OneBodyG input
    file, in its order: ``G[i + N s, j + N t] = <c^dag_{i s} c_{j t}>`` of
    the saved (mixed) Green function, ``N = Nsite`` (uhfr.py:1302-1314).

    Attributes
    ----------
    green : ndarray (2N, 2N) complex
        The entries of the file at ``site + N * spin``; 0 where the file has
        no entry.
    pairs : ndarray (m, 4) int
        ``(i, s, j, t)`` of each line, in file order.
    mask : ndarray (2N, 2N) bool
        True where the file has an entry.
    source : Path
    meta : dict
        ``nsite``, ``two_sz``, ``onebodyg`` (the OneBodyG file name in the
        input, or None) and ``index`` (the spin-orbital numbering).
    """
    green: np.ndarray
    pairs: np.ndarray
    mask: np.ndarray
    source: Path
    meta: dict


def read_uhfr_green(run_or_toml, *, base_dir=None) -> UHFrGreen:
    """Read the ``green`` text file of a UHFr run (a
    :class:`hwave.results.Run` or the path of its input TOML file).

    Raises MissingOutputError when ``[file.output] green`` is not configured
    or its file is absent (UHFr writes it only when ``[file.input]
    OneBodyG`` is set) or when the status file records zero SCF steps
    (:func:`hwave.results.status.require_scf_steps`), and
    UnsupportedLayoutError for a run that is not
    UHFr, for a line that is not ``i s j t Re Im`` with ``0 <= i, j <
    Nsite``, ``s, t`` in ``{0, 1}`` and finite values, and for an entry
    ``(i, s, j, t)`` given on more than one line (naming the lines)."""
    run = _uhfr_run(run_or_toml, base_dir)
    nsite = _nsite(run)
    two_sz = _two_sz(run)
    path = _uhfr_output(run, "green")
    try:
        with open(path, "r", encoding="utf-8") as f:
            lines = f.read().splitlines()
    except (OSError, UnicodeDecodeError) as exc:
        raise UnsupportedLayoutError("{} cannot be read as text: {}".format(path, exc),
                                     path=path) from None

    def bad(lineno, line, why):
        return UnsupportedLayoutError(
            "{}: line {} is not 'i s j t Re Im' ({}): {!r}".format(
                path, lineno, why, line), path=path)

    green = np.zeros((2 * nsite, 2 * nsite), dtype=np.complex128)
    mask = np.zeros((2 * nsite, 2 * nsite), dtype=bool)
    pairs = []
    first_line = {}
    for lineno, line in enumerate(lines, start=1):
        fields = line.split()
        if not fields:
            continue
        if len(fields) != 6:
            raise bad(lineno, line, "{} fields".format(len(fields)))
        if not all(_INT.fullmatch(x) for x in fields[:4]):
            raise bad(lineno, line, "the indices are not integers")
        i, s, j, t = (int(x) for x in fields[:4])
        if not (0 <= i < nsite and 0 <= j < nsite):
            raise bad(lineno, line, "a site is outside 0..{}".format(nsite - 1))
        if s not in (0, 1) or t not in (0, 1):
            raise bad(lineno, line, "a spin is not 0 or 1")
        try:
            re_part, im_part = float(fields[4]), float(fields[5])
        except ValueError:
            raise bad(lineno, line, "the values are not numbers") from None
        if not (math.isfinite(re_part) and math.isfinite(im_part)):
            raise bad(lineno, line, "a value is not finite")
        key = (i, s, j, t)
        if key in first_line:
            raise UnsupportedLayoutError(
                "{}: the entry ({} {} {} {}) is given twice, on lines {} and {}"
                .format(path, i, s, j, t, first_line[key], lineno), path=path,
                hint="list each entry once in the OneBodyG file and rerun")
        first_line[key] = lineno
        a, b = i + nsite * s, j + nsite * t
        green[a, b] = complex(re_part, im_part)
        mask[a, b] = True
        pairs.append(key)
    pairs = np.array(pairs, dtype=np.int64).reshape(-1, 4)
    return UHFrGreen(green=_freeze(green), pairs=_freeze(pairs), mask=_freeze(mask),
                     source=Path(path),
                     meta={"nsite": nsite, "two_sz": two_sz,
                           "onebodyg": _onebodyg_name(run),
                           "index": "site + Nsite * spin"})


def read_uhfr_order_parameters(run_or_toml, *, base_dir=None, site_positions=None,
                               lattice_shape=None):
    """The order parameters of a UHFr run (a :class:`hwave.results.Run` or
    the path of its input TOML file), from its ``green`` file.

    Every UHFr site is one site: ``site_of_orbital = arange(Nsite)`` and
    ``rho[i] = [[G[i, i], G[i, i+N]], [G[i+N, i], G[i+N, i+N]]]``
    (``N = Nsite``). The cluster is one magnetic cell: ``cell_shape =
    lattice_shape`` and ``n_cell = 1``.

    ``site_positions`` (shape ``(Nsite, 3)``) are the site coordinates in
    units of the lattice vectors, and ``lattice_shape`` the extent of the
    cluster in those units (e.g. ``(4, 4, 1)`` for a periodic 4 x 4
    cluster); the wave vectors are then ``q_a = m_a / L_a``. Without either
    of them only the per-site and uniform quantities are given (the Fourier
    and staggered fields are None).

    The green file must hold ``(i 0 i 0)`` and ``(i 1 i 1)`` for every site
    ``i`` and, when 2Sz is not fixed, also ``(i 0 i 1)`` and ``(i 1 i 0)``.
    With a fixed 2Sz the spin blocks are separate (uhfr.py:553-560): the
    spin-off-diagonal entries are zero by construction, and
    ``transverse_available`` is False; on-site spin-off-diagonal lines in the
    file are then not used, and ``meta["spin_offdiagonal_ignored"]`` (only
    set with a fixed 2Sz) counts them.

    Raises MissingOutputError naming the missing entries (and the errors of
    :func:`read_uhfr_green`); ValueError for bad ``site_positions`` or
    ``lattice_shape`` (from
    :func:`hwave.results.order_parameters_from_density`)."""
    run = _uhfr_run(run_or_toml, base_dir)
    result = read_uhfr_green(run)
    nsite, two_sz = result.meta["nsite"], result.meta["two_sz"]
    spins = ((0, 0), (1, 1)) if two_sz is not None else ((0, 0), (1, 1), (0, 1), (1, 0))
    present = {tuple(p) for p in result.pairs.tolist()}
    missing = [(i, s, i, t) for (s, t) in spins for i in range(nsite)
               if (i, s, i, t) not in present]
    if missing:
        missing.sort()
        listed = ", ".join("({} {} {} {})".format(*m) for m in missing[:_MAX_LISTED])
        if len(missing) > _MAX_LISTED:
            listed += " and {} more".format(len(missing) - _MAX_LISTED)
        name = result.meta["onebodyg"]
        target = "the OneBodyG file {!r}".format(name) if name is not None \
            else "the OneBodyG file"
        raise MissingOutputError(
            "{}: the order parameters need the Green function entries {}, "
            "which the file does not have".format(result.source, listed),
            path=result.source,
            hint="add these lines to {} and rerun".format(target))

    g = result.green
    site = np.arange(nsite)
    rho = np.zeros((nsite, 2, 2), dtype=np.complex128)
    rho[:, 0, 0] = g[site, site]
    rho[:, 1, 1] = g[site + nsite, site + nsite]
    if two_sz is None:
        rho[:, 0, 1] = g[site, site + nsite]
        rho[:, 1, 0] = g[site + nsite, site]
    op = order_parameters_from_density(
        rho, site_of_orbital=site, site_positions=site_positions,
        cell_shape=lattice_shape, lattice_shape=lattice_shape, n_cell=1,
        transverse_available=two_sz is None)
    meta = dict(op.meta)
    meta.update(solver="UHFr", nsite=nsite, two_sz=two_sz,
                green_index="site + Nsite * spin")
    if two_sz is not None:
        meta["spin_offdiagonal_ignored"] = sum(
            1 for i, s, j, t in present if i == j and s != t)
    return dataclasses.replace(op, source=result.source, meta=meta)


# ===================================================================== UHFk
#
# UHFk layouts (uhfk.py):
#
# * ``green.npz`` member ``green``: ``G[r, s, a, t, b] = <c^dag_{s a}(r)
#   c_{t b}(0)>`` on the CellShape lattice, ``r = rz + Nz (ry + Ny rx)``
#   (uhfk.py:1443-1444; docs green_uhfk.rst).
# * member ``green_sublattice`` (when ``prod(SubShape) > 1``, SubShape
#   defaulting to CellShape, uhfk.py:234-238): the same on the lattice of
#   magnetic cells ``CellShape / SubShape``, with the folded orbital
#   ``aa = a + norb_orig * ir`` and ``ir = irx + Bx (iry + By irz)``
#   (uhfk.py:498-500, 532-540).
# * spin-orbital mode: ``ns = 1`` and the orbital index is interleaved,
#   ``a = 2 k + spin`` (spin_orbital.py:151-158).

_KNOWN_ISSUE_357 = (
    "the UHFk onebodyg text output uses a displacement convention that is "
    "inconsistent with the Green function (issue #357)")

#: the tolerance on orbital centres that must coincide
CENTRE_ATOL = 1e-8
#: the largest |NCond - round(NCond)| for which HOMO and LUMO are given
NCOND_INTEGER_ATOL = 1e-8

HOMO_LUMO_NOTE = ("spectral neighbours at integer NCond; not a band gap unless "
                  "the bands are separated")


def _uhfk_run(run_or_toml, base_dir):
    run = as_run(run_or_toml, base_dir=base_dir)
    if run.mode != "UHFk":
        raise UnsupportedLayoutError(
            "{} is a {!r} run, not a UHFk run".format(run.input_path, run.mode),
            path=run.input_path)
    return run


def _shape_param(run, name, default=None):
    value = run.param(name, default)
    if value is None:
        raise UnsupportedLayoutError(
            "[mode.param] {} is missing in {}; UHFk requires it".format(
                name, run.input_path), path=run.input_path)
    if not isinstance(value, (list, tuple)) or len(value) != 3 or not all(
            isinstance(v, numbers.Integral) and not isinstance(v, bool) and v >= 1
            for v in value):
        raise UnsupportedLayoutError(
            "[mode.param] {} must be three positive integers in {}, got {!r}"
            .format(name, run.input_path, value), path=run.input_path)
    return tuple(int(v) for v in value)


@dataclasses.dataclass(frozen=True)
class _Lattice:
    cell_shape: tuple        # CellShape, the lattice of the run
    sub_shape: tuple         # SubShape (default CellShape), the magnetic cell
    super_shape: tuple       # CellShape / SubShape, the lattice of magnetic cells

    @property
    def subvol(self):
        return int(np.prod(self.sub_shape))

    @property
    def nvol(self):
        return int(np.prod(self.super_shape))


def _lattice(run):
    """CellShape and SubShape as the solver reads them (uhfk.py:224-251)."""
    cell = _shape_param(run, "CellShape")
    sub = _shape_param(run, "SubShape", list(cell))
    if any(c % b for c, b in zip(cell, sub)):
        raise UnsupportedLayoutError(
            "[mode.param] SubShape {} does not divide CellShape {} in {}".format(
                sub, cell, run.input_path), path=run.input_path)
    return _Lattice(cell, sub, tuple(c // b for c, b in zip(cell, sub)))


def _r_index(shape):
    """``r_index[r] = (rx, ry, rz)`` for ``r = rz + Nz (ry + Ny rx)``."""
    return np.indices(shape).reshape(3, -1).T.astype(np.int64)


def _uhfk_output(run, key):
    """``run.path(key)``, after :func:`require_scf_steps`."""
    require_scf_steps(run, key)
    return run.path(key)


@dataclasses.dataclass(frozen=True)
class GreenR:
    """A UHFk real-space Green function ``G[r, s, a, t, b] = <c^dag_{s
    a}(r) c_{t b}(0)>``.

    Attributes
    ----------
    green : ndarray (nr, ns, norb, ns, norb) complex, or None
        ``ns = 2`` (``ns = 1`` and ``a = 2 k + spin`` in spin-orbital mode).
        ``None`` when the member is not read (a deflated ``green`` member
        without the ``green_convention`` marker, read with ``legacy=True``;
        see ``meta["legacy_notes"]``) or absent.
    r_index : ndarray (nr, 3) int
        ``(rx, ry, rz)`` of each ``r``, with ``r = rz + Nz (ry + Ny rx)``.
    cell_shape : tuple of int
        The lattice ``(Nx, Ny, Nz)`` that ``r`` runs over: CellShape for the
        ``green`` member, ``CellShape / SubShape`` for the sublattice one.
    sublattice : GreenR or None
        The ``green_sublattice`` member: the same on the lattice of magnetic
        cells, with the folded orbital ``aa = a + norb_orig * ir`` and ``ir
        = irx + Bx (iry + By irz)`` within the cell ``SubShape = (Bx, By,
        Bz)``. ``None`` when the file has no such member.
    spin_orbital : bool
    source : Path
    meta : dict
        ``member`` (the npz member), ``sub_shape``, ``legacy_notes`` (the
        assumptions of :func:`hwave.results._io.check_markers`) and, for a
        ``green`` that is not read, ``green_unavailable``.
    """
    green: Optional[np.ndarray]
    r_index: np.ndarray
    cell_shape: tuple
    sublattice: Optional["GreenR"]
    spin_orbital: bool
    source: Path
    meta: dict


def _check_green_shape(path, name, shape, nr, ns):
    if len(shape) != 5 or shape[0] != nr or shape[1] != ns or shape[3] != ns \
            or shape[2] != shape[4] or shape[2] < 1:
        raise UnsupportedLayoutError(
            "{}: member {!r} has shape {}, expected ({}, {}, norb, {}, norb) for "
            "this run".format(path, name, shape, nr, ns, ns), path=path)
    return shape[2]


def read_uhfk_green(run_or_toml, *, legacy=False, base_dir=None) -> GreenR:
    """Read ``green.npz`` of a UHFk run (a :class:`hwave.results.Run` or
    the path of its input TOML file).

    The convention markers are checked with the legacy table of
    :mod:`hwave.results._io`: a run with ``prod(SubShape) > 1`` must have
    ``green_sublattice``; a file without ``momentum_convention`` needs
    ``legacy=True``; a file with ``green_sublattice`` but no
    ``green_convention`` is read from ``green_sublattice`` only (its
    deflated ``green`` may use the fold of the files written before PR #35,
    issue #36), so ``green`` is then ``None``.

    Raises MissingOutputError when ``green`` is not configured (it is by
    default; ``green = ""`` disables it) or its file is absent, or when the
    status file records zero SCF steps
    (:func:`hwave.results.status.require_scf_steps`); UnsupportedLayoutError
    for a run that is not UHFk, for a marker or member the table refuses,
    and for members whose shapes do not fit CellShape and SubShape."""
    run = _uhfk_run(run_or_toml, base_dir)
    lattice = _lattice(run)
    path = _uhfk_output(run, "green")
    kind = "uhfk_green_sublattice" if lattice.subvol > 1 else "uhfk_green"
    notes = check_markers(path, kind, legacy=legacy)
    shapes = member_shapes(path)
    ns = 1 if run.spin_orbital else 2
    dropped = any(note["member"] == "green" and note["assumed"] is None
                  for note in notes)
    meta_common = {"sub_shape": lattice.sub_shape, "legacy_notes": notes}

    norb = None
    if "green" in shapes:
        norb = _check_green_shape(path, "green", shapes["green"],
                                  int(np.prod(lattice.cell_shape)), ns)
    if "green_sublattice" in shapes:
        norb_sub = _check_green_shape(path, "green_sublattice",
                                      shapes["green_sublattice"], lattice.nvol, ns)
        if norb_sub % lattice.subvol or (
                norb is not None and norb_sub != norb * lattice.subvol):
            raise UnsupportedLayoutError(
                "{}: member 'green_sublattice' has {} orbitals, which is not "
                "{} orbitals x {} cell positions of SubShape {}".format(
                    path, norb_sub, norb if norb is not None else "an integer",
                    lattice.subvol, lattice.sub_shape), path=path)

    with open_npz(path) as data:
        green = None
        if "green" in shapes and not dropped:
            green = npz_member(data, path, "green", "c")
        sub = None
        if "green_sublattice" in shapes:
            sub = GreenR(green=_freeze(npz_member(data, path, "green_sublattice", "c")),
                         r_index=_freeze(_r_index(lattice.super_shape)),
                         cell_shape=lattice.super_shape, sublattice=None,
                         spin_orbital=run.spin_orbital, source=Path(path),
                         meta=dict(_freeze(meta_common), member="green_sublattice"))
    meta = dict(_freeze(meta_common), member="green")
    if green is None:
        meta["green_unavailable"] = (
            "not read: no green_convention marker (issue #36)" if dropped
            else "the file has no 'green' member")
    return GreenR(green=_freeze(green), r_index=_freeze(_r_index(lattice.cell_shape)),
                  cell_shape=lattice.cell_shape, sublattice=sub,
                  spin_orbital=run.spin_orbital, source=Path(path), meta=meta)


@dataclasses.dataclass(frozen=True)
class Bands:
    """The UHFk eigenvalues (and eigenvectors) on the folded k grid.

    Attributes
    ----------
    k_index : ndarray (nk, 3) int
        The archive's ``wavevector_index``.
    k_unit : ndarray (3, 3)
        The archive's ``wavevector_unit``.
    k_cart : ndarray (nk, 3)
        ``k_index @ k_unit``.
    energies : ndarray (nk, nd)
        The archive's ``eigenvalue``, in its column order (the blocks of
        the Hamiltonian one after another, uhfk.py:1839-1875); not sorted.
    sort_index : ndarray (nk, nd) int
        ``np.argsort(energies, axis=1, kind="stable")``.
    eigenvectors : ndarray (nk, nd, nd) complex, or None
        Column ``l`` is the eigenvector of ``energies[k, l]``; only with
        ``with_vectors=True``.
    single_reservoir : bool
        True when ``energy.dat`` has a plain ``ChemicalPotential`` key (one
        chemical potential for all electrons).
    mu : float or None
        The chemical potential of a single reservoir.
    mu_groups : dict[int, float] or None
        ``{g: mu}`` of the ``ChemicalPotential_<g>`` keys, in file order;
        None for a single reservoir.
    homo, lumo, gap : float or None
        For a single reservoir and an integer ``NCond = n_occ`` (within
        1e-8) with ``0 < n_occ < nk * nd``: ``E_flat[n_occ - 1]``,
        ``E_flat[n_occ]`` and their difference, with ``E_flat`` the sorted
        eigenvalues of every k. They are spectral neighbours, not a band gap
        unless the bands are separated (``meta["homo_lumo_note"]``).
        Otherwise None, with the reason in ``meta["homo_lumo_unavailable"]``.
    spin_orbital : bool
    source : Path
        The eigen archive.
    meta : dict
        ``legacy_notes``, ``ncond`` (``NCond`` of ``energy.dat``),
        ``energy_source``, ``energy_unavailable`` and the HOMO/LUMO notes.
    """
    k_index: np.ndarray
    k_unit: np.ndarray
    k_cart: np.ndarray
    energies: np.ndarray
    sort_index: np.ndarray
    eigenvectors: Optional[np.ndarray]
    single_reservoir: bool
    mu: Optional[float]
    mu_groups: Optional[dict]
    homo: Optional[float]
    lumo: Optional[float]
    gap: Optional[float]
    spin_orbital: bool
    source: Path
    meta: dict


def read_bands(run_or_toml, *, with_vectors=False, legacy=False,
               base_dir=None) -> Bands:
    """Read the ``eigen`` archive of a UHFk run (a
    :class:`hwave.results.Run` or the path of its input TOML file), with the
    chemical potential and ``NCond`` of ``energy.dat``.

    When ``energy`` is not configured or its file is absent, the bands are
    still returned, with ``single_reservoir = False`` and no ``mu``,
    ``homo``, ``lumo`` or ``gap`` (``meta["energy_unavailable"]``).

    Raises MissingOutputError when ``eigen`` is not configured or its file is
    absent, or when the status file records zero SCF steps;
    UnsupportedLayoutError for a run that is not UHFk, for a marker the
    legacy table refuses (an archive without ``momentum_convention`` needs
    ``legacy=True``), for members of unexpected shapes, and for a
    ``ChemicalPotential_<g>`` key whose suffix is not an integer."""
    run = _uhfk_run(run_or_toml, base_dir)
    lattice = _lattice(run)
    path = _uhfk_output(run, "eigen")
    notes = check_markers(path, "uhfk_eigen", legacy=legacy)
    shapes = member_shapes(path)
    for name in ("eigenvalue", "eigenvector", "wavevector_unit", "wavevector_index"):
        if name not in shapes:
            raise UnsupportedLayoutError("{}: no member {!r}".format(path, name),
                                         path=path)
    nk_nd = shapes["eigenvalue"]
    if len(nk_nd) != 2 or nk_nd[0] != lattice.nvol:
        raise UnsupportedLayoutError(
            "{}: member 'eigenvalue' has shape {}, expected ({}, nd) for the "
            "{} k points of this run".format(path, nk_nd, lattice.nvol,
                                             lattice.nvol), path=path)
    nk, nd = nk_nd
    for name, shape in (("eigenvector", (nk, nd, nd)), ("wavevector_unit", (3, 3)),
                        ("wavevector_index", (nk, 3))):
        if shapes[name] != shape:
            raise UnsupportedLayoutError(
                "{}: member {!r} has shape {}, expected {}".format(
                    path, name, shapes[name], shape), path=path)
    with open_npz(path) as data:
        energies = npz_member(data, path, "eigenvalue", "f")
        k_unit = npz_member(data, path, "wavevector_unit", "f")
        k_index = npz_member(data, path, "wavevector_index", "i")
        vectors = npz_member(data, path, "eigenvector", "c") if with_vectors else None

    meta = {"legacy_notes": notes}
    single, mu, groups, ncond = False, None, None, None
    try:
        energy = read_energy(run)
    except MissingOutputError as exc:
        meta["energy_unavailable"] = str(exc)
    else:
        meta["energy_source"] = energy.source
        ncond = energy.ncond
        if isinstance(energy.mu, dict):
            groups = {}
            for suffix, value in energy.mu.items():
                if not _INT.fullmatch(suffix):
                    raise UnsupportedLayoutError(
                        "{}: ChemicalPotential_{} is not a UHFk mu group (an "
                        "integer suffix)".format(energy.source, suffix),
                        path=energy.source)
                groups[int(suffix)] = value
        elif energy.mu is not None:
            single, mu = True, float(energy.mu)
    meta["ncond"] = ncond

    homo = lumo = gap = None
    if not single:
        meta["homo_lumo_unavailable"] = (
            "energy.dat has no plain ChemicalPotential (one per mu group, or "
            "none): the electrons are not in a single reservoir")
    elif ncond is None or not abs(ncond - round(ncond)) <= NCOND_INTEGER_ATOL:
        meta["homo_lumo_unavailable"] = "NCond {!r} is not an integer".format(ncond)
    elif not 0 < round(ncond) < nk * nd:
        meta["homo_lumo_unavailable"] = (
            "NCond {} leaves no occupied or no empty state among the {} "
            "eigenvalues".format(ncond, nk * nd))
    else:
        n_occ = int(round(ncond))
        flat = np.sort(energies, axis=None)
        homo, lumo = float(flat[n_occ - 1]), float(flat[n_occ])
        gap = lumo - homo
        meta["homo_lumo_note"] = HOMO_LUMO_NOTE

    return Bands(k_index=_freeze(k_index), k_unit=_freeze(k_unit),
                 k_cart=_freeze(k_index @ k_unit), energies=_freeze(energies),
                 sort_index=_freeze(np.argsort(energies, axis=1, kind="stable")),
                 eigenvectors=_freeze(vectors), single_reservoir=single, mu=mu,
                 mu_groups=groups, homo=homo, lumo=lumo, gap=gap,
                 spin_orbital=run.spin_orbital, source=Path(path),
                 meta=_freeze(meta))


def so_density_from_green(G0, norb_phys_orig, subvol):
    """The per-orbital spin density of a spin-orbital UHFk Green function at
    ``r = 0``.

    ``G0`` is ``Green[0, 0, :, 0, :]`` of shape ``(nd, nd)``, ``nd = 2 *
    norb_phys_orig * subvol``, in the interleaved folded index ``aa = 2 j +
    spin`` with ``j = k + norb_phys_orig * ir`` (the folded index of
    uhfk.py:532-540 applied to the interleaved ``a = 2 k + spin``,
    spin_orbital.py:151-158). Returns a new array ``rho`` of shape
    ``(norb_phys_orig * subvol, 2, 2)`` with ``rho[j, s, t] = G0[2 j + s,
    2 j + t]``; entry ``j`` is physical orbital ``j % norb_phys_orig`` at
    cell position ``j // norb_phys_orig``. Raises ValueError for a bad
    shape or count."""
    for name, value in (("norb_phys_orig", norb_phys_orig), ("subvol", subvol)):
        if isinstance(value, bool) or not isinstance(value, numbers.Integral) \
                or value < 1:
            raise ValueError("{} must be a positive integer, got {!r}"
                             .format(name, value))
    g0 = np.asarray(G0)
    n = int(norb_phys_orig) * int(subvol)
    if g0.shape != (2 * n, 2 * n):
        raise ValueError(
            "G0 must have shape ({0}, {0}) for norb_phys_orig = {1} and subvol = "
            "{2}, got {3}".format(2 * n, norb_phys_orig, subvol, g0.shape))
    j = np.arange(n)
    return np.array(g0.reshape(n, 2, n, 2)[j, :, j, :], dtype=np.complex128)


def _site_of_orbital(value, norb_phys):
    """The validated ``site_of_orbital`` argument (length ``norb_phys``,
    integers covering ``0..n_site - 1``)."""
    if value is None:
        return np.arange(norb_phys)
    raw = np.asarray(value)
    if raw.ndim != 1 or len(raw) != norb_phys:
        raise ValueError(
            "site_of_orbital must give one site per physical orbital of the "
            "Geometry ({} values), got shape {}".format(norb_phys, raw.shape))
    if raw.dtype == np.bool_ or not np.issubdtype(raw.dtype, np.number) \
            or np.iscomplexobj(raw) or not all(
                float(v).is_integer() for v in raw.tolist()):
        raise ValueError("site_of_orbital must hold integers, got {!r}"
                         .format(raw.tolist()))
    sites = raw.astype(np.int64)
    n_site = int(sites.max()) + 1 if len(sites) else 0
    missing = sorted(set(range(n_site)) - set(sites.tolist()))
    if np.any(sites < 0) or missing:
        raise ValueError(
            "site_of_orbital values must cover 0..n_site-1 with no gap, got {!r}"
            .format(sites.tolist()))
    return sites


def _orbital_centres(geometry, so):
    """The centre of each physical orbital ``k``: Geometry entry ``k``, or
    entries ``2k`` and ``2k+1`` in spin-orbital mode, which must agree."""
    centres = np.asarray(geometry.centres)
    if not so:
        return centres
    even, odd = centres[0::2], centres[1::2]
    for k in range(len(even)):
        if np.max(np.abs(even[k] - odd[k])) > CENTRE_ATOL:
            raise UnsupportedLayoutError(
                "Geometry file {}: the spin-orbital entries {} and {} of physical "
                "orbital {} have different centres {} and {}".format(
                    geometry.source, 2 * k, 2 * k + 1, k, even[k].tolist(),
                    odd[k].tolist()), path=geometry.source)
    return even


def read_uhfk_order_parameters(run_or_toml, *, site_of_orbital=None,
                               site_positions=None, legacy=False, base_dir=None):
    """The order parameters of a UHFk run (a :class:`hwave.results.Run` or
    the path of its input TOML file), from ``green.npz`` at ``r = 0``.

    The magnetic cell is SubShape (default CellShape): the density is that
    of ``green_sublattice`` when ``prod(SubShape) > 1``, else of ``green``,
    at ``r = 0``. Its entries are ``j = k + norb_phys_orig * ir``, physical
    orbital ``k`` at cell position ``ir = irx + Bx (iry + By irz)``:

    * normal mode: ``rho[j, s, t] = G[0, s, j, t, j]``;
    * spin-orbital mode: :func:`so_density_from_green` of ``G[0, 0, :, 0,
      :]``.

    ``site_of_orbital`` (length ``norb_phys_orig``, default ``arange``)
    groups the physical orbitals of one cell position into sites; entry
    ``j`` is on site ``ir * n_site_orig + site_of_orbital[k]``. The site
    positions are ``R_ir + tau_site`` in fractional coordinates, with
    ``R_ir = (irx, iry, irz)`` and ``tau`` the Geometry centre of the
    site's orbitals (which must coincide within 1e-8), or
    ``site_positions[site]`` (shape ``(n_site_orig, 3)``) when given. Then
    :func:`hwave.results.order_parameters_from_density` is called with
    ``cell_shape = SubShape``, ``lattice_shape = CellShape`` and ``n_cell``
    the number of magnetic cells. ``transverse_available`` is True: UHFk
    stores the spin-off-diagonal entries in every mode (with a fixed 2Sz and
    no spin-mixing term they are zero, and so are S^x and S^y;
    ``meta["transverse_note"]`` then says so).

    Raises MissingOutputError for a missing Geometry or Green function file
    (and the errors of :func:`read_uhfk_green`); UnsupportedLayoutError when
    the Geometry does not fit the Green function, or when the two
    spin-orbital entries of one physical orbital have different centres;
    ValueError for a bad ``site_of_orbital`` or ``site_positions``, or when
    the orbitals of one site have different centres (naming them)."""
    run = _uhfk_run(run_or_toml, base_dir)
    lattice = _lattice(run)
    green = read_uhfk_green(run, legacy=legacy)
    if lattice.subvol > 1:
        member, values = green.sublattice, green.sublattice.green
    else:
        member, values = green, green.green
    if values is None:
        raise UnsupportedLayoutError(
            "{}: the 'green' member is not readable ({})".format(
                green.source, green.meta.get("green_unavailable")),
            path=green.source)

    geometry = read_geometry(input_file_path(run, "Geometry"))
    so = run.spin_orbital
    norb_orig = geometry.norb
    if so and norb_orig % 2:
        raise UnsupportedLayoutError(
            "Geometry file {}: {} entries cannot be the interleaved spin-orbitals "
            "of spin-orbital mode (an even number)".format(geometry.source, norb_orig),
            path=geometry.source)
    norb_phys = norb_orig // 2 if so else norb_orig
    nfold = values.shape[2]
    if nfold != norb_orig * lattice.subvol:
        raise UnsupportedLayoutError(
            "the Geometry file {} has {} orbital entries, but member {!r} of {} has "
            "{} folded orbitals, which is not {} x {} cell positions".format(
                geometry.source, norb_orig, member.meta["member"], green.source,
                nfold, norb_orig, lattice.subvol),
            path=geometry.source)
    g0 = values[0]
    if so:
        rho = so_density_from_green(g0[0, :, 0, :], norb_phys, lattice.subvol)
    else:
        j = np.arange(nfold)
        rho = g0[:, j, :, j]             # rho[j, s, t] = G[0, s, j, t, j]

    groups = _site_of_orbital(site_of_orbital, norb_phys)
    n_site = int(groups.max()) + 1
    if site_positions is not None:
        tau = np.asarray(site_positions)
        if tau.shape != (n_site, 3) or not np.issubdtype(tau.dtype, np.number) \
                or np.iscomplexobj(tau) or not np.all(np.isfinite(tau)):
            raise ValueError(
                "site_positions must be finite real fractional coordinates of shape "
                "({}, 3), one row per site of site_of_orbital, got {}".format(
                    n_site, tau.shape))
        tau = np.array(tau, dtype=float)
        tau_source = "site_positions"
    else:
        centres = _orbital_centres(geometry, so)
        tau = np.zeros((n_site, 3))
        for site in range(n_site):
            orbitals = np.flatnonzero(groups == site)
            spread = np.max(np.abs(centres[orbitals] - centres[orbitals[0]]))
            if spread > CENTRE_ATOL:
                raise ValueError(
                    "site {} groups the physical orbitals {} whose Geometry centres "
                    "{} differ by more than {:g}; pass site_positions or group "
                    "only orbitals at one centre".format(
                        site, orbitals.tolist(), centres[orbitals].tolist(),
                        CENTRE_ATOL))
            tau[site] = centres[orbitals[0]]
        tau_source = "geometry"

    bx, by, _ = lattice.sub_shape
    ir = np.arange(lattice.subvol)
    cell_pos = np.stack([ir % bx, (ir // bx) % by, ir // (bx * by)], axis=1)
    entries = np.arange(norb_phys * lattice.subvol)
    entry_sites = (entries // norb_phys) * n_site + groups[entries % norb_phys]
    positions = (cell_pos[:, None, :] + tau[None, :, :]).reshape(-1, 3)

    op = order_parameters_from_density(
        rho, site_of_orbital=entry_sites, site_positions=positions,
        cell_shape=lattice.sub_shape, lattice_shape=lattice.cell_shape,
        n_cell=lattice.nvol, transverse_available=True)
    meta = dict(op.meta)
    meta.update(solver="UHFk", spin_orbital=so, green_member=member.meta["member"],
                legacy_notes=_freeze(green.meta["legacy_notes"]),
                geometry=geometry.source, norb_phys_orig=norb_phys,
                site_of_orbital_arg=_freeze(groups), tau_source=tau_source,
                entry_index="j = k + norb_phys_orig * ir, ir = irx + Bx (iry + "
                            "By irz)")
    two_sz = _two_sz(run)
    if two_sz is not None:
        meta["two_sz"] = two_sz
        meta["transverse_note"] = "transverse fixed to zero by 2Sz"
    return dataclasses.replace(op, source=green.source, meta=meta)


@dataclasses.dataclass(frozen=True)
class OneBodyG:
    """The UHFk ``onebodyg`` text output, lines ``i s j t Re Im``.

    Not read by this version: :func:`read_onebodyg` refuses the file (issue
    #357). The class is kept so that the result type is part of the API.

    Attributes
    ----------
    entries : ndarray (m, 4) int
    values : ndarray (m,) complex
    source : Path
    meta : dict
    """
    entries: np.ndarray
    values: np.ndarray
    source: Path
    meta: dict


def read_onebodyg(run_or_toml, *, base_dir=None) -> OneBodyG:
    """Refuse the ``onebodyg`` text output of a UHFk run.

    UHFk writes ``onebodyg`` with a displacement convention that does not
    match its Green function (the displacement sign and the order in which
    it is packed, issue #357), so the file is not read. Use
    :func:`read_uhfk_green` (``green.npz``, ``G[r, s, a, t, b] = <c^dag_{s
    a}(r) c_{t b}(0)>``) instead.

    Raises MissingOutputError when ``onebodyg`` is not configured, when its
    file does not exist, or when the run performed no SCF step;
    UnsupportedLayoutError (issue #357) otherwise, and for a run that is not
    UHFk."""
    run = _uhfk_run(run_or_toml, base_dir)
    require_scf_steps(run, "onebodyg")
    path = run.path("onebodyg")      # MissingOutputError: not configured, absent
    raise UnsupportedLayoutError(
        "{}: {}".format(path, _KNOWN_ISSUE_357),
        path=path,
        hint="read green.npz with read_uhfk_green instead; see issue #357")
