"""Readers of the UHFr and UHFk outputs.

* :func:`read_energy`: ``energy.dat`` (:class:`Energy`).
* :func:`read_uhfr_eigen`: the UHFr eigenvalues and eigenvectors
  (:class:`UHFrEigen`).
* :func:`read_uhfr_green`: the UHFr Green function text (:class:`UHFrGreen`).
* :func:`read_uhfr_order_parameters`: the order parameters of a UHFr run
  (:class:`hwave.results.OrderParameters`).

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
from hwave.results._io import member_shapes, open_npz, read_key_value_text
from hwave.results.errors import (
    MissingOutputError, ResultsError, UnsupportedLayoutError)
from hwave.results.order import order_parameters_from_density
from hwave.results.run import as_run

__all__ = ["Energy", "read_energy", "UHFrEigen", "read_uhfr_eigen",
           "UHFrGreen", "read_uhfr_green", "read_uhfr_order_parameters"]

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
    is absent, and UnsupportedLayoutError for a line that is not
    ``key = number``, for a key given twice (compared case-insensitively),
    or for a file with both ``ChemicalPotential`` and
    ``ChemicalPotential_<suffix>``."""
    run = as_run(run_or_toml, base_dir=base_dir)
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


def _uhfr_output(run, key, what):
    """``run.path(key)``; when the file is configured but absent, the hint
    says why the solver may not have written it."""
    try:
        return run.path(key)
    except MissingOutputError as exc:
        if key not in run.outputs:
            raise
        hint = exc.hint
        steps = None
        try:
            from hwave.results.status import read_scf_status
            steps = read_scf_status(run).iterations
        except ResultsError:
            pass
        if steps == 0:
            hint = ("the run performed no SCF step (IterationMax <= 0), so no "
                    "{} was written".format(what))
        elif key == "green" and _onebodyg_name(run) is None:
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
    or a block's file is absent, and UnsupportedLayoutError for a run that
    is not UHFr, or for a file without the ``eigenvalue`` and
    ``eigenvector`` members of the expected shapes."""
    run = _uhfr_run(run_or_toml, base_dir)
    nsite = _nsite(run)
    two_sz = _two_sz(run)
    blocks = _SZ_FREE_BLOCKS if two_sz is None else _FIXED_SZ_BLOCKS
    size = 2 * nsite if two_sz is None else nsite
    energies, vectors, sources = {}, {}, []
    for block in blocks:
        path = _uhfr_output(run, "eigen:" + block, "eigenvector file")
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
            energies[block] = np.asarray(data["eigenvalue"])
            if with_vectors:
                vectors[block] = np.asarray(data["eigenvector"])
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
    OneBodyG`` is set), and UnsupportedLayoutError for a run that is not
    UHFr, for a line that is not ``i s j t Re Im`` with ``0 <= i, j <
    Nsite``, ``s, t`` in ``{0, 1}`` and finite values, and for an entry
    ``(i, s, j, t)`` given on more than one line (naming the lines)."""
    run = _uhfr_run(run_or_toml, base_dir)
    nsite = _nsite(run)
    two_sz = _two_sz(run)
    path = _uhfr_output(run, "green", "Green function")
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
    ``transverse_available`` is False.

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
    return dataclasses.replace(op, source=result.source, meta=meta)
