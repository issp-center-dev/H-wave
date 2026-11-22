"""Low-level access to H-wave files for the readers of :mod:`hwave.results`.

* :func:`open_npz`, :func:`member_shapes`, :func:`scalar_members`: ``.npz``
  archives, without pickles; the last two never load an array member.
* :func:`read_key_value_text`: the ``key = value`` text of ``energy.dat`` and
  ``scf_status.dat``.
* :func:`read_geometry`: the UHFk Geometry file.
* :func:`check_markers`: the per-marker legacy table of the UHFk ``.npz``
  outputs.
* :func:`input_file_path`: where the solver read an input file of
  ``[file.input.interaction]`` from.

Legacy table (``check_markers``), one row per marker:

* UHFk ``eigen.npz`` and ``green.npz``, marker ``momentum_convention``: when
  absent, assumable as ``e_plus_ikR`` with ``legacy=True``; UHFk always used
  this sign (uhfk.py:1969-1975).
* UHFk ``green.npz`` with ``green_sublattice`` (written whenever
  ``prod(SubShape or CellShape) > 1``, uhfk.py:234-238), marker
  ``green_convention``: never assumed. Before PR #35 the deflated ``green`` member used the
  Hamiltonian fold (uhfk.py:1960-1967; rpa.py:3725-3735; issue #36). When
  the marker is absent, ``green_sublattice`` is used and the deflated
  ``green`` is refused; a file without ``green_sublattice`` is refused.
* ``energy.dat``: no marker; always readable (the energy reader lists absent
  keys).
* ``scf_status.dat``, marker ``format_version``: checked by the status
  reader, not here.
"""
import contextlib
import dataclasses
import os
import zipfile
from collections.abc import Mapping
from pathlib import Path

import numpy as np

from hwave.results._frozen import _freeze
from hwave.results.errors import MissingOutputError, UnsupportedLayoutError

__all__ = ["open_npz", "npz_member", "member_shapes", "scalar_members", "read_key_value_text",
           "Geometry", "read_geometry", "MARKER_KINDS", "check_markers",
           "input_file_path"]

#: the sign convention UHFk writes (uhfk.py:1872, 1969-1975)
MOMENTUM_CONVENTION = "e_plus_ikR"
#: the fold convention of the deflated ``green`` member since PR #35
GREEN_CONVENTION = "green_slot_first"

_LEGACY_HINT = "pass legacy=True to read it with the assumption listed in meta"


def _missing(path, what="file"):
    return MissingOutputError("{} not found: {}".format(what, path), path=path)


# --------------------------------------------------------------------- npz


@contextlib.contextmanager
def open_npz(path):
    """``np.load(path, allow_pickle=False)`` as a context manager that closes
    the archive on exit, also on an error. Raises MissingOutputError for a
    missing file and UnsupportedLayoutError for a file that is not an
    ``.npz`` archive."""
    path = Path(path)
    if not path.is_file():
        raise _missing(path)
    try:
        data = np.load(path, allow_pickle=False)
    except (OSError, ValueError, zipfile.BadZipFile) as exc:
        raise UnsupportedLayoutError("{} is not a readable .npz archive: {}"
                                     .format(path, exc), path=path) from None
    if not isinstance(data, np.lib.npyio.NpzFile):
        raise UnsupportedLayoutError("{} is a .npy array, not an .npz archive"
                                     .format(path), path=path)
    try:
        yield data
    finally:
        data.close()


_MEMBER_KINDS = {"i": "integer", "f": "real", "c": "numeric"}


def npz_member(data, path, name, kind):
    """Member ``name`` of the open archive ``data`` (from :func:`open_npz`,
    the file ``path``), checked against ``kind``:

    * ``"i"``: integer within the int64 range; a float member is accepted
      when every value is such an integer, and the result is int64;
    * ``"f"``: real, as float; a complex member is accepted only when every
      imaginary part is zero;
    * ``"c"``: any integer, real or complex member, in its own dtype.

    Raises UnsupportedLayoutError naming the member for a member that
    cannot be loaded (an object member needs pickles, which are never
    loaded) or does not fit ``kind``."""
    try:
        arr = np.asarray(data[name])
    except (ValueError, OSError, EOFError, zipfile.BadZipFile) as exc:
        raise UnsupportedLayoutError("{}: member {!r} cannot be read: {}".format(
            path, name, exc), path=path) from None
    got = arr.dtype.kind
    bad = None
    if got not in "iufc":
        bad = "has dtype {}".format(arr.dtype)
    elif kind == "i":
        if got == "c" or (got == "f" and not (np.all(np.isfinite(arr))
                                               and np.all(arr == np.round(arr)))):
            bad = "has non-integer values (dtype {})".format(arr.dtype)
        elif arr.size and (
                (got == "u" and arr.max() > np.iinfo(np.int64).max)
                or (got == "f" and not (np.all(arr >= -2.0 ** 63)
                                        and np.all(arr < 2.0 ** 63)))):
            bad = "has values outside the int64 range (dtype {})".format(arr.dtype)
        else:
            arr = arr.astype(np.int64)
    elif kind == "f":
        if got == "c":
            if np.any(arr.imag != 0):
                bad = "has non-zero imaginary parts"
            else:
                arr = arr.real
        if bad is None:
            arr = arr.astype(float)
    if bad is not None:
        raise UnsupportedLayoutError("{}: member {!r} {}, expected {} values".format(
            path, name, bad, _MEMBER_KINDS[kind]), path=path)
    return arr


@contextlib.contextmanager
def _members(path):
    """``(zipfile, {member name: zip entry name})`` of the ``.npy`` members."""
    path = Path(path)
    if not path.is_file():
        raise _missing(path)
    try:
        zf = zipfile.ZipFile(path)
    except (OSError, zipfile.BadZipFile) as exc:
        raise UnsupportedLayoutError("{} is not a readable .npz archive: {}"
                                     .format(path, exc), path=path) from None
    try:
        names = {}
        for entry in zf.namelist():
            if entry.endswith(".npy"):
                names[entry[:-len(".npy")]] = entry
        yield zf, names
    finally:
        zf.close()


def _bad_member(path, name, exc):
    return UnsupportedLayoutError("{}: member {!r} cannot be read: {}"
                                  .format(path, name, exc), path=path)


def member_shapes(path):
    """``{member: shape}`` of the ``.npz`` archive ``path``, read from the
    ``.npy`` headers only (:func:`hwave.solver.npy_header.read_npy_header_shape`,
    npy_header.py:205); no array data is read. Zip entries that are not
    ``.npy`` members are skipped."""
    from hwave.solver.npy_header import read_npy_header_shape
    shapes = {}
    with _members(path) as (zf, names):
        for name, entry in names.items():
            try:
                with zf.open(entry) as fh:
                    shapes[name] = tuple(read_npy_header_shape(fh))
            except Exception as exc:          # malformed header, bad zip entry
                raise _bad_member(path, name, exc) from None
    return shapes


def _python_scalar(arr):
    if arr.dtype.kind == "U":
        return str(arr.item())
    return arr.item()


def scalar_members(path):
    """``{member: value}`` of the 0-d members of the ``.npz`` archive ``path``
    (markers such as ``momentum_convention``), as Python scalars: ``str`` for
    unicode text, else ``.item()`` (a bytes member stays ``bytes``). Only 0-d members are read
    (``numpy.lib.format.read_array`` on the zip member, no pickles); the
    others are identified from their headers and never loaded."""
    from hwave.solver.npy_header import read_npy_header_shape
    values = {}
    with _members(path) as (zf, names):
        for name, entry in names.items():
            try:
                with zf.open(entry) as fh:
                    shape = tuple(read_npy_header_shape(fh))
                if shape != ():
                    continue
                with zf.open(entry) as fh:
                    arr = np.lib.format.read_array(fh, allow_pickle=False)
                values[name] = _python_scalar(arr)
            except Exception as exc:          # object dtype, bad header, ...
                raise _bad_member(path, name, exc) from None
    return values


# -------------------------------------------------------------------- text


def read_key_value_text(path):
    """The ``key = value`` lines of the text file ``path``, in file order, as
    ``[(key, value, lineno)]`` (``lineno`` counts from 1). Each line is split
    on its first ``=`` and both sides are stripped; blank lines are skipped.
    Duplicate keys are kept: callers that want a mapping build it
    themselves. A non-blank line without ``=`` or with an empty key raises
    UnsupportedLayoutError."""
    path = Path(path)
    if not path.is_file():
        raise _missing(path)
    try:
        with open(path, "r", encoding="utf-8") as f:
            lines = f.read().splitlines()
    except (OSError, UnicodeDecodeError) as exc:
        raise UnsupportedLayoutError("{} cannot be read as text: {}".format(path, exc),
                                     path=path) from None
    rows = []
    for lineno, line in enumerate(lines, start=1):
        if not line.strip():
            continue
        key, sep, value = line.partition("=")
        key = key.strip()
        if not sep or not key:
            raise UnsupportedLayoutError(
                "{}: line {} is not a 'key = value' line: {!r}".format(path, lineno, line),
                path=path)
        rows.append((key, value.strip(), lineno))
    return rows


# ---------------------------------------------------------------- geometry


@dataclasses.dataclass(frozen=True)
class Geometry:
    """A UHFk Geometry file.

    Attributes
    ----------
    lattice_vectors : ndarray, shape (3, 3)
        Row ``i`` is the lattice vector ``a_i``.
    centres : ndarray, shape (norb, 3)
        The orbital centres, in fractional coordinates of the lattice vectors.
    source : Path
        The file read.
    """
    lattice_vectors: np.ndarray
    centres: np.ndarray
    source: Path

    @property
    def norb(self):
        """The number of orbitals."""
        return self.centres.shape[0]


def read_geometry(path):
    """Read the Geometry file ``path``: lines 1-3 the lattice vectors (the
    rule of dos.py:113-126), line 4 the number of orbitals ``norb``, then one
    centre per orbital (wan90.py:54-72; docs/en/source/uhfk/filespecification/
    inputfiles_uhfk/geometry.rst). As in the solver, each vector line gives
    its first three numbers. Unlike the solver, which leaves missing centres
    at zero, a file with fewer or more centre lines than ``norb`` (blank
    lines aside) raises UnsupportedLayoutError."""
    path = Path(path)
    if not path.is_file():
        raise _missing(path, "Geometry file")
    try:
        with open(path, "r", encoding="utf-8") as f:
            lines = [line.strip() for line in f.read().splitlines()]
    except (OSError, UnicodeDecodeError) as exc:
        raise UnsupportedLayoutError("{} cannot be read as text: {}".format(path, exc),
                                     path=path) from None

    def bad(message):
        return UnsupportedLayoutError("Geometry file {}: {}".format(path, message),
                                      path=path)

    def vector(lineno):
        if lineno > len(lines):
            raise bad("line {} is missing".format(lineno))
        values = lines[lineno - 1].split()
        try:
            if len(values) < 3:
                raise ValueError
            return [float(v) for v in values[:3]]
        except ValueError:
            raise bad("line {} must hold three numbers, got {!r}"
                      .format(lineno, lines[lineno - 1])) from None

    lattice = np.array([vector(i) for i in (1, 2, 3)], dtype=float)
    if len(lines) < 4:
        raise bad("line 4 (the number of orbitals) is missing")
    try:
        norb = int(lines[3])
    except ValueError:
        raise bad("line 4 must be the number of orbitals, got {!r}"
                  .format(lines[3])) from None
    if norb < 1:
        raise bad("the number of orbitals must be positive, got {}".format(norb))
    for offset, line in enumerate(lines[4:4 + norb]):
        if not line:
            raise bad("line {} is blank, but it must hold the centre of orbital "
                      "{} of {}".format(5 + offset, offset + 1, norb))
    found = len([line for line in lines[4:] if line])
    if found != norb:
        raise bad("{} orbital centres are declared on line 4 but {} centre lines "
                  "follow".format(norb, found))
    centres = np.array([vector(5 + i) for i in range(norb)], dtype=float)
    return Geometry(lattice_vectors=_freeze(lattice), centres=_freeze(centres),
                    source=path)


# ---------------------------------------------------------- legacy markers

#: the ``kind`` values of :func:`check_markers`. ``uhfk_green_sublattice`` is
#: the ``green.npz`` of a UHFk run with the sublattice layout, i.e.
#: ``prod(SubShape or CellShape) > 1`` (SubShape defaults to CellShape,
#: uhfk.py:234-238); ``uhfk_green`` one with ``prod == 1``. The kind only
#: decides whether a missing ``green_sublattice`` is refused: the
#: green_convention rule applies to every file that has ``green_sublattice``.
MARKER_KINDS = ("uhfk_eigen", "uhfk_green", "uhfk_green_sublattice")

_PAYLOAD = {"uhfk_eigen": "eigenvalue", "uhfk_green": "green",
            "uhfk_green_sublattice": "green_sublattice"}


def _marker_value(path, name, shapes, scalars, expected):
    """The value of marker ``name`` (None when absent); a marker that is not
    a single text value (a unicode string, as the solvers write it; bytes
    are refused) equal to ``expected`` raises."""
    if name not in shapes:
        return None
    if name not in scalars or not isinstance(scalars[name], str):
        raise UnsupportedLayoutError(
            "{}: marker {!r} must be a single text value".format(path, name),
            path=path, hint="regenerate the file with the current H-wave")
    if scalars[name] != expected:
        raise UnsupportedLayoutError(
            "{} records {} = {!r}; this version of hwave.results reads {!r} only"
            .format(path, name, scalars[name], expected), path=path,
            hint="regenerate the file with the current H-wave")
    return scalars[name]


def check_markers(path, kind, *, legacy):
    """Check the convention markers of the UHFk ``.npz`` output ``path`` of
    ``kind`` (one of :data:`MARKER_KINDS`) against the legacy table of this
    module. No array member is loaded.

    Returns the list of legacy notes, each ``{"member", "assumed",
    "reason"}``: ``assumed`` is the value assumed for an absent marker
    (only with ``legacy=True``), or ``None`` for a member that must not be
    used (the deflated ``green`` of an untagged sublattice file). The list is
    empty for a file with every marker.

    Raises UnsupportedLayoutError for a marker with another value, for a
    missing payload member, for an absent marker without ``legacy=True``
    (the message then mentions ``legacy=True``; only assumable rows do), and
    for kind ``uhfk_green_sublattice`` without ``green_sublattice``. The
    ``green`` note is given for every file with ``green_sublattice``, whatever
    ``kind`` is. Raises
    ValueError for an unknown ``kind``."""
    if kind not in MARKER_KINDS:
        raise ValueError("unknown marker kind {!r}; expected one of {}"
                         .format(kind, ", ".join(MARKER_KINDS)))
    path = Path(path)
    shapes = member_shapes(path)
    scalars = scalar_members(path)
    notes = []

    payload = _PAYLOAD[kind]
    has_payload = payload in shapes or (
        kind == "uhfk_green" and "green_sublattice" in shapes)
    if not has_payload:
        if kind == "uhfk_green_sublattice":
            raise UnsupportedLayoutError(
                "{} has no green_sublattice member; its deflated 'green' member "
                "cannot be folded back safely without the green_convention marker "
                "(issue #36)".format(path), path=path,
                hint="regenerate the file with the current H-wave")
        raise UnsupportedLayoutError("{} has no {!r} member; members: {}".format(
            path, payload, ", ".join(sorted(shapes)) or "none"), path=path)

    if _marker_value(path, "momentum_convention", shapes, scalars,
                     MOMENTUM_CONVENTION) is None:
        if not legacy:
            raise UnsupportedLayoutError(
                "{} has no momentum_convention marker (written before the "
                "marker was added)".format(path), path=path, hint=_LEGACY_HINT)
        notes.append({
            "member": "momentum_convention",
            "assumed": MOMENTUM_CONVENTION,
            "reason": "the file has no momentum_convention marker; UHFk always "
                      "used the e^{+ikR} sign (uhfk.py:1969-1975, issue #133)",
        })

    # the file decides: every file with green_sublattice has the sublattice
    # layout, whatever kind the caller derived from the input
    if kind != "uhfk_eigen" and _marker_value(
            path, "green_convention", shapes, scalars, GREEN_CONVENTION) is None:
        if "green_sublattice" in shapes and "green" in shapes:
            notes.append({
                "member": "green",
                "assumed": None,
                "reason": "the file has no green_convention marker, so its "
                          "deflated 'green' member may use the fold of the files "
                          "written before PR #35 (issue #36); green_sublattice is "
                          "used instead and 'green' is not read",
            })
    return notes


# ------------------------------------------------------------- input files


def _exact(mapping, key, default=None):
    """``mapping[key]`` with the exact case of ``key``, as qlms.run reads the
    tables and ``[file.input] path_to_input`` (qlms.py:104-108); the
    mappings of :attr:`Run.config` keep the case of the TOML keys."""
    for k, v in mapping.items():
        if k == key:
            return v
    return default


def input_file_path(run, key, *, must_exist=True):
    """The path of the input file ``[file.input.interaction] key`` of ``run``
    (a :class:`hwave.results.Run`), resolved as the solver resolves it, from
    ``run.root``:

    * ``key`` is matched case-insensitively, as both solvers do
      (qlms.py:146-149 for UHFr; read_input_k.py:93-118 for UHFk);
    * UHFr joins the name to ``[file.input] path_to_input``, a key the solver
      reads with its exact case (qlms.py:108, 149);
    * UHFk joins it to ``[file.input.interaction] path_to_input``, matched
      case-insensitively (read_input_k.py:94);
    * both default to ``""``; an absolute name is used as is.

    The ``[file]``, ``[file.input]`` and ``[file.input.interaction]`` tables
    are found with their exact case, as qlms.run does.

    Raises MissingOutputError when the key is not configured, or when
    ``must_exist`` is true and the file does not exist;
    UnsupportedLayoutError for a name that is not a string, or for a mode
    other than UHFr and UHFk."""
    if run.mode not in ("UHFr", "UHFk"):
        raise UnsupportedLayoutError(
            "input files of mode {!r} are not resolved by this version of "
            "hwave.results".format(run.mode), path=run.input_path)
    tables = []
    table = run.config
    for where, name in (("[file]", "file"), ("[file.input]", "input"),
                        ("[file.input.interaction]", "interaction")):
        table = _exact(table, name, {})
        if not isinstance(table, Mapping):
            raise UnsupportedLayoutError("{} is not a table in {}".format(
                where, run.input_path), path=run.input_path)
        tables.append(table)
    info_input, interaction = tables[1], tables[2]
    names = {k.lower(): v for k, v in interaction.items()}
    lookup = key.lower()
    if lookup == "path_to_input" or lookup not in names:
        raise MissingOutputError(
            "input file {!r} is not configured in {}".format(key, run.input_path),
            path=run.input_path,
            hint="configure it in [file.input.interaction] of the input file")
    name = names[lookup]
    if run.mode == "UHFr":
        directory = _exact(info_input, "path_to_input", "")
    else:
        directory = names.get("path_to_input", "")
    for what, value in (("[file.input.interaction] {}".format(key), name),
                        ("path_to_input", directory)):
        if not isinstance(value, str):
            raise UnsupportedLayoutError("{} must be a string in {}, got {!r}"
                                         .format(what, run.input_path, value),
                                         path=run.input_path)
    path = run.root / os.path.join(directory, name)
    if must_exist and not path.is_file():
        raise MissingOutputError(
            "input file {!r}: {} does not exist; input paths were resolved against "
            "{} (base directory source: {})".format(key, path, run.root,
                                                    run.base_dir_source),
            path=path,
            hint="if H-wave was started from another directory, pass "
                 "base_dir=<that directory>")
    return path
