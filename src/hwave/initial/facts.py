"""Model facts: what the enumerator needs to know about a UHFk input.

:func:`model_facts` loads the input through the solver facade and records
the orbital counts, the modes, the interacting orbitals, the site groups,
the lattice, the temperature, the electron count, warnings and the SHA-256
of every input file the TOML names.

Every orbital index here is a **physical** orbital of the original
(primitive) cell, in spin-orbital mode too (where the Geometry lists the
spin-orbitals ``2a`` and ``2a + 1`` of each physical orbital ``a``). This
is the index of PatternSpec ``orb_a``/``orb_b``, of ``[mode.param]
orbital_site`` and of the ``site_of_orbital`` labels of
``hwave.results``.
"""
import hashlib
import itertools
import os
import pathlib
from dataclasses import dataclass
from fractions import Fraction
from types import MappingProxyType
from typing import Mapping, Optional, Tuple, Union

import numpy as np

from hwave.solver import initial_facade

__all__ = ["ModelFacts", "SiteGroups", "model_facts", "cluster_sites",
           "SITE_TOL", "EXTERN_WARNING"]

#: centres closer than this (in the length unit of the Geometry lattice
#: vectors, minimum image) belong to one site (owner-approved value)
SITE_TOL = 1e-3

#: the warning for an ``Extern`` table (design section 5, step 1)
EXTERN_WARNING = ("Extern is not used by UHFk; give a Zeeman field as on-site "
                  "Transfer in spin-orbital mode")

#: the warning for an input without ``SubShape`` (design section 5, step 2)
NO_SUBSHAPE_WARNING = (
    "The input sets no SubShape, so UHFk uses SubShape = CellShape (one "
    "supercell, the Gamma point only); give the candidate SubShape list "
    "explicitly.")

#: lattice images searched per axis beyond the nearest one (a lattice
#: vector shorter than the site tolerance is flagged long before this)
_MAX_EXTRA_IMAGES = 3


@dataclass(frozen=True)
class SiteGroups:
    """Which orbitals share a site.

    ``source`` is ``"orbital_site"`` (``[mode.param] orbital_site``),
    ``"user"`` (the ``orbital_site`` argument of :func:`model_facts`) or
    ``"clustered"`` (orbital centres closer than :data:`SITE_TOL`).
    ``groups`` holds the physical orbital indices of each site, sorted, ordered by
    their first orbital. ``ambiguous`` holds one message per clustered
    group that is wider than the tolerance or joins two images of one
    orbital; such a grouping should be given explicitly.
    """
    source: str
    groups: Tuple[Tuple[int, ...], ...]
    ambiguous: Tuple[str, ...]


# eq=False: ``problem`` holds arrays
@dataclass(frozen=True, eq=False)
class ModelFacts:
    """The facts of one UHFk input (design section 5, step 1; section 4.2).

    ``problem``: the :class:`~hwave.solver.initial_facade.UHFkProblem`.
    ``norb``: orbitals of the original cell (Geometry ``norb``; spin
    included in spin-orbital mode). ``norb_phys``: physical orbitals of the
    original cell. ``spin_orbital``, ``iflag_fock``.

    ``interacting_orbitals``: the physical orbitals with a nonzero entry in
    any two-body table of the original cell (before the sublattice fold;
    two-body indices are physical in every mode). ``n_c``: their number.

    ``site_groups`` (:class:`SiteGroups`). ``cellshape``, ``T``, ``ncond``
    (as the solver resolves ``Ncond`` and ``Filling``), ``ncond_per_cell``
    (``ncond`` over the number of primitive cells; a ``Fraction`` when
    ``ncond`` is an integer), ``two_sz`` (``None`` when Sz is free).

    ``warnings``: messages about the input. ``input_sha256``: the SHA-256
    of every file the TOML names (:func:`initial_facade.referenced_input_files`),
    keyed by its TOML key; the TOML itself is not included.
    """
    problem: "initial_facade.UHFkProblem"
    norb: int
    norb_phys: int
    spin_orbital: bool
    iflag_fock: bool
    interacting_orbitals: Tuple[int, ...]
    n_c: int
    site_groups: SiteGroups
    cellshape: Tuple[int, int, int]
    T: float
    ncond: Union[int, float]
    ncond_per_cell: Union[Fraction, float]
    two_sz: Optional[int]
    warnings: Tuple[str, ...]
    input_sha256: Mapping[str, str]


def _read_toml(path):
    import tomli
    with open(path, "rb") as f:
        return tomli.load(f)


def _mode_of(path):
    try:
        mode = _read_toml(path).get("mode", {})
        return repr(mode.get("mode")) if isinstance(mode, Mapping) else "(no [mode] table)"
    except Exception as exc:  # an unreadable TOML is still listed
        return "(unreadable: {})".format(exc)


def _resolve_toml(input_path):
    """The TOML of ``input_path``: the path itself, or the one ``*.toml``
    of a directory whose mode is UHFk."""
    p = pathlib.Path(os.fspath(input_path))
    if p.is_dir():
        found = sorted(p.glob("*.toml"))
        if len(found) != 1:
            listing = ", ".join("{} (mode {})".format(f.name, _mode_of(f)) for f in found)
            raise ValueError("{}: expected exactly one *.toml input, found {}{}".format(
                p, len(found), ": " + listing if found else " (no *.toml file)"))
        toml = found[0]
        mode = _read_toml(toml).get("mode", {})
        if not isinstance(mode, Mapping) or mode.get("mode") != "UHFk":
            hint = ""
            if isinstance(mode, Mapping) and mode.get("mode") == "UHFr":
                hint = ("; a UHFr input has no lattice geometry: give the UHFk input of "
                        "the same lattice, and write the UHFr initial file with "
                        "hwave.initial.writer.initial_state(..., site_map=, uhfr_toml=)")
            raise ValueError("{}: the only *.toml input, {} (mode {}), is not a UHFk "
                             "input{}".format(p, toml.name, _mode_of(toml), hint))
        return toml
    if not p.is_file():
        raise ValueError("{}: no such input file or directory".format(p))
    return p


def _interacting(param_ham_orig, spin_orbital):
    phys = set()
    for name in initial_facade.TWO_BODY_TYPES:
        table = param_ham_orig.get(name)
        if not table:
            continue
        for (_irvec, orbvec), v in table.items():
            if v != 0:
                phys.update(int(a) for a in orbvec)
    phys = tuple(sorted(phys))
    return phys, len(phys)


def _groups_from_labels(labels):
    by_label = {}
    for i, lab in enumerate(labels):
        by_label.setdefault(lab, []).append(i)
    return tuple(sorted(tuple(g) for g in by_label.values()))


def _shift_ranges(df, heights, tol):
    """Integer shifts ``n`` per axis with ``|df + n| * h <= tol`` possible."""
    out = []
    for d, h in zip(df, heights):
        reach = tol / h if h > 0 else np.inf
        reach = min(reach, _MAX_EXTRA_IMAGES + 0.5)
        lo = int(np.ceil(-d - reach))
        hi = int(np.floor(-d + reach))
        out.append(range(lo, hi + 1))
    return out


def cluster_sites(centers, rvec, tol=SITE_TOL):
    """Group orbitals by their centres.

    ``centers`` ``(norb, 3)`` are fractional coordinates, ``rvec`` the
    lattice vectors (rows). Two orbitals are linked when some lattice image
    of one lies within ``tol`` of the other (Cartesian distance); the
    groups are the connected components (a transitive join). A group is
    ambiguous when it is wider than ``tol`` (the largest distance between
    its members, placed by the links) or when its links join two images of
    one orbital.

    Returns ``(groups, ambiguous)`` as in :class:`SiteGroups`.
    """
    centers = np.asarray(centers, dtype=float).reshape(-1, 3)
    rvec = np.asarray(rvec, dtype=float).reshape(3, 3)
    n = centers.shape[0]
    vol = abs(float(np.linalg.det(rvec)))
    heights = [vol / float(np.linalg.norm(np.cross(rvec[(i + 1) % 3], rvec[(i + 2) % 3])))
               for i in range(3)]
    links = []          # (i, j, shift): centre_j + shift lies within tol of centre_i
    for i in range(n):
        for j in range(i, n):
            df = centers[j] - centers[i]
            for shift in itertools.product(*_shift_ranges(df, heights, tol)):
                if i == j and shift == (0, 0, 0):
                    continue
                if np.linalg.norm((df + np.array(shift)) @ rvec) <= tol:
                    links.append((i, j, shift))
    # connected components
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for i, j, _ in links:
        parent[find(i)] = find(j)
    comps = {}
    for i in range(n):
        comps.setdefault(find(i), []).append(i)
    groups = tuple(sorted(tuple(sorted(c)) for c in comps.values()))

    ambiguous = []
    for g in groups:
        members = set(g)
        glinks = [(i, j, s) for i, j, s in links if i in members]
        # place every member by the links (breadth first from the first)
        image = {g[0]: np.zeros(3, dtype=int)}
        adj = {}
        for i, j, s in glinks:
            adj.setdefault(i, []).append((j, np.array(s)))
            adj.setdefault(j, []).append((i, -np.array(s)))
        queue = [g[0]]
        while queue:
            a = queue.pop(0)
            for b, s in adj.get(a, []):
                if b not in image:
                    image[b] = image[a] + s
                    queue.append(b)
        images = any(not np.array_equal(image[j] - image[i], np.array(s))
                     for i, j, s in glinks)
        if images:
            ambiguous.append("orbitals {} join two lattice images of one orbital "
                             "within site_tol = {:g}".format(
                                 ", ".join(map(str, g)), tol))
            continue
        pos = np.array([(centers[a] + image[a]) @ rvec for a in g])
        width = max((float(np.linalg.norm(p - q)) for p in pos for q in pos), default=0.0)
        if width > tol:
            ambiguous.append("orbitals {} form a group {:.3g} wide, wider than "
                             "site_tol = {:g}".format(", ".join(map(str, g)), width, tol))
    return groups, tuple(ambiguous)


def _min_image_distance(f1, f2, rvec):
    d = np.asarray(f2, dtype=float) - np.asarray(f1, dtype=float)
    return float(np.linalg.norm((d - np.round(d)) @ rvec))


def _physical_centers(problem):
    """The centre of each physical orbital of the original cell. In
    spin-orbital mode the Geometry lists ``2a`` and ``2a + 1`` for orbital
    ``a``; their centres must coincide within :data:`SITE_TOL`."""
    geom = problem.geometry
    centers = np.asarray(geom["center"], dtype=float).reshape(-1, 3)
    rvec = np.asarray(geom["rvec"], dtype=float).reshape(3, 3)
    if not problem.enable_spin_orbital:
        return centers, rvec
    bad = []
    for a in range(centers.shape[0] // 2):
        dist = _min_image_distance(centers[2 * a], centers[2 * a + 1], rvec)
        if dist > SITE_TOL:
            bad.append("orbital {} (spin-orbitals {} and {}, {:.3g} apart)".format(
                a, 2 * a, 2 * a + 1, dist))
    if bad:
        raise ValueError(
            "spin-orbital mode: the Geometry centres of the two spin-orbitals 2a, "
            "2a+1 of a physical orbital must coincide (within site_tol = {:g}); "
            "they do not for {}".format(SITE_TOL, "; ".join(bad)))
    return centers[0::2], rvec


def _site_groups(problem, norb_phys, user, centers, rvec):
    if user is not None:
        labels = _user_labels(user, norb_phys, problem)
        return SiteGroups("user", _groups_from_labels(labels), ())
    if problem.orbital_site is not None:
        return SiteGroups("orbital_site", _groups_from_labels(problem.orbital_site), ())
    groups, ambiguous = cluster_sites(centers, rvec, SITE_TOL)
    return SiteGroups("clustered", groups, ambiguous)


def _user_labels(user, norb, problem):
    if isinstance(user, (str, bytes)) or not hasattr(user, "__len__"):
        raise ValueError("orbital_site must be a list with one site label per orbital, "
                         "got {!r}".format(user))
    labels = []
    for v in user:
        if isinstance(v, bool) or not isinstance(v, (int, np.integer)) or v < 0:
            raise ValueError("orbital_site must hold non-negative integers, got {!r}"
                             .format(user))
        labels.append(int(v))
    if len(labels) != norb:
        hint = ""
        if problem.enable_spin_orbital and len(labels) == problem.norb_orig:
            hint = ("; {} is the spin-orbital count: give one label per physical "
                    "orbital (orbital a covers spin-orbitals 2a and 2a+1)"
                    .format(len(labels)))
        raise ValueError("orbital_site has {} entries but the original cell has {} "
                         "physical orbitals (one label per physical orbital){}"
                         .format(len(labels), norb, hint))
    return labels


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _input_sha256(toml):
    out = {}
    for key, path in initial_facade.referenced_input_files(toml).items():
        if not path.is_file():
            raise ValueError("{}: {} = {!r} names a file that does not exist ({})"
                             .format(toml, key, path.name, path))
        out[key] = _sha256(path)
    return MappingProxyType(out)


def _subshape_given(toml):
    param = _read_toml(toml).get("mode", {}).get("param", {}) or {}
    # [mode.param] keys are case-insensitive (ParamDict)
    return any(k.lower() == "subshape" for k in param)


def model_facts(input_path, *, subshape=None, orbital_site=None):
    """The :class:`ModelFacts` of a UHFk input.

    Parameters
    ----------
    input_path : path-like
        A TOML file, or a directory holding exactly one ``*.toml``, whose
        ``[mode] mode`` is UHFk.
    subshape : tuple of int, optional
        Overrides ``[mode.param] SubShape`` for the loaded problem.
    orbital_site : list of int, optional
        A site label per physical orbital of the original cell; overrides
        ``[mode.param] orbital_site`` and the clustering.

    Raises
    ------
    ValueError
        No single UHFk TOML, a refused input (see
        :func:`initial_facade.load_problem`), a named file that is
        missing, an ``orbital_site`` of the wrong length, or (spin-orbital
        mode) two spin-orbitals of one physical orbital at different
        centres.
    """
    toml = _resolve_toml(input_path)
    problem = initial_facade.load_problem(toml, subshape=subshape)
    so = problem.enable_spin_orbital
    norb = problem.norb_orig
    norb_phys = norb // 2 if so else norb
    interacting, n_c = _interacting(problem.param_ham_orig, so)
    centers, rvec = _physical_centers(problem)

    warnings = []
    if "Extern" in problem.param_ham_orig:
        warnings.append(EXTERN_WARNING)
    if not _subshape_given(toml):
        warnings.append(NO_SUBSHAPE_WARNING)

    ncell = int(np.prod(problem.cellshape))
    ncond = problem.ncond
    if float(ncond).is_integer():
        ncond_per_cell = Fraction(int(ncond), ncell)
    else:
        ncond_per_cell = ncond / ncell

    return ModelFacts(
        problem=problem, norb=norb, norb_phys=norb_phys, spin_orbital=so,
        iflag_fock=problem.iflag_fock, interacting_orbitals=interacting, n_c=n_c,
        site_groups=_site_groups(problem, norb_phys, orbital_site, centers, rvec),
        cellshape=tuple(problem.cellshape), T=problem.T, ncond=ncond,
        ncond_per_cell=ncond_per_cell, two_sz=problem.two_sz,
        warnings=tuple(warnings), input_sha256=_input_sha256(toml))
