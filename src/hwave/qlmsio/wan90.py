"""Functions for reading/writing Wannier90 format files.

This module provides functions to read and write geometry and Hamiltonian data
in Wannier90 format.

Functions
---------
read_geom(name_in)
    Read geometry data from file.
read_geometry(name_in) 
    Read extended geometry data from file.
write_geom(name_out, info_geometry)
    Write geometry data to file.
read_w90_header_dim(name_in)
    Orbital count declared on line 2 of a Wannier90 file.
read_w90_header_marker(name_in)
    Format marker (``w90spin/<version>``) on line 1, or None.
read_w90(name_in, allow_spin_columns=False)
    Read Wannier90 Hamiltonian data from file (7-column, or the
    spin-resolved 9-column Transfer format when allowed).
read_w90_header_dim(name_in)
    Read the orbital count declared on line 2 of a Wannier90 file.
write_w90(name_in, info_int, info_geometry, interact_shape)
    Write Wannier90 Hamiltonian data to file.
"""

from __future__ import print_function

import itertools
import re
import numpy as np
import logging

logger = logging.getLogger("qlms").getChild("wan90")

# Spin-resolved ("extended column") Transfer format, issue #323. Line 1 of
# the file -- a free comment in Wannier90 -- carries the token
# ``w90spin/<version>``; every body row then has nine columns
#     Rx Ry Rz  a  b  s  t  Re  Im
# with ``a, b`` 1-based PHYSICAL orbital indices (1..norb of line 2) and
# ``s, t`` 0-based spins (0 = up, 1 = down). The reader folds each row to the
# spin-block composite index ``s*norb + a`` so the returned table keeps the
# 2-tuple key shape of the 7-column format. The marker, not the column
# count, decides the format: an unmarked 9-column file is refused, never
# guessed at.
SPIN_COLUMN_FORMAT = "w90spin"
SPIN_COLUMN_VERSION = 1
_SPIN_COLUMN_MARKER = re.compile(
    r"\b" + re.escape(SPIN_COLUMN_FORMAT) + r"/(\d+)\b")
_NCOL_PLAIN = 7
_NCOL_SPIN = 9

def read_geom(name_in):
    """Read geometry data from file.
    
    Parameters
    ----------
    name_in : str
        Input filename
        
    Returns
    -------
    dict
        Dictionary containing:
        - norb : int
            Number of orbitals
        - rvec : ndarray
            Real space lattice vectors (3x3)
        - center : ndarray
            Orbital center positions (norb x 3)
            
    Raises
    ------
    SystemExit
        If file not found
    """
    try:
        with open(name_in, 'r') as f:
            # skip header
            l_strip = [s.strip() for s in f.readlines()]
        norb = int(l_strip[3])
    except OSError:
        logger.error("read_geom: file {} not found".format(name_in))
        exit(1)

    rvec = np.zeros((3, 3), dtype=float)
    center = np.zeros((norb, 3), dtype=float)

    for count, line in enumerate(l_strip[:3]):
        values = line.split()
        rvec[count] = (float(values[0]), float(values[1]), float(values[2]))

    for count, line in enumerate(l_strip[4:]):
        values = line.split()
        center[count] = (float(values[0]), float(values[1]), float(values[2]))

    data = {"norb": norb, "rvec": rvec, "center": center}
    return data

def read_geometry(name_in):
    """Read extended geometry data from file.
    
    Parameters
    ----------
    name_in : str
        Input filename
        
    Returns
    -------
    dict
        Dictionary containing:
        - unit_vec : ndarray
            Unit cell vectors
        - degree : ndarray
            Degrees of freedom
        - cell_vec : ndarray
            Cell vectors
        - site2vec : dict
            Mapping of sites to vectors
        - n_orb : int
            Number of orbitals
            
    Raises
    ------
    SystemExit
        If file not found
    """
    info_geometry = {}
    unit_vec_line_end = 3
    cell_vec_line_start = 4
    cell_vec_line_end = 7
    n_kpath_line = 2

    try:
        with open(name_in, "r") as fr:
            lines = fr.readlines()
    except OSError:
        logger.error("read_geom: file {} not found".format(name_in))
        exit(1)
        
    info_geometry["unit_vec"] = np.array([data.split() for data in lines[0:unit_vec_line_end]], dtype=np.float)
    info_geometry["degree"] = np.array(lines[unit_vec_line_end].split(), dtype=np.float)
    info_geometry["cell_vec"] = np.array([data.split() for data in lines[cell_vec_line_start:cell_vec_line_end]], dtype=np.float)
    info_geometry["site2vec"] = {}
    for idx, data in enumerate(lines[cell_vec_line_end:]):
        if len(data.split()) == n_kpath_line:
            break
        else:
            vector = np.array(data.split(), dtype=np.int32)
            info_geometry["site2vec"][idx] = vector
    norb = 0
    for idx, vec in info_geometry["site2vec"].items():
        if [vec[0], vec[1], vec[2]] == [0, 0, 0]:
            norb += 1
    info_geometry["n_orb"] = norb
    return info_geometry

def write_geom(name_out, info_geometry):
    """Write geometry data to file.
    
    Parameters
    ----------
    name_out : str
        Output filename
    info_geometry : dict
        Geometry information dictionary
    """
    with open(name_out, "w") as fw:
        #Unit_cell
        for vec in info_geometry["unit_vec"]:
            fw.write("{} {} {}\n".format(vec[0], vec[1], vec[2]))
        #norb
        n_orb = info_geometry["n_orb"]
        fw.write("{}\n".format(n_orb))
        #wannier_center (temporary)
        for idx in range(n_orb):
            pos = idx*1.0/(n_orb+1)
            fw.write("{} {} {}\n".format(pos, pos, pos))

def read_w90_header_dim(name_in):
    """Orbital count declared in the header of a Wannier90 file.

    The first line of the file is a comment; line 2 (0-based line index 1)
    is the number of orbitals (Wannier90 ``num_wann``). :func:`read_w90`
    skips this line, so the declared dimension is otherwise lost; it is
    read here, separately, so callers can check a table against it.

    Parameters
    ----------
    name_in : str
        Input filename

    Returns
    -------
    int
        The declared number of orbitals (a positive integer).

    Raises
    ------
    ValueError
        If line 2 is missing, blank, not an integer, or not positive.
    """
    with open(name_in, 'r') as f:
        lines = f.readlines()
    raw = lines[1].strip() if len(lines) > 1 else None
    try:
        dim = int(raw)
    except (TypeError, ValueError):
        dim = None
    if dim is None or dim <= 0:
        raise ValueError(
            "{}: line 2 must be the number of orbitals (a positive integer), "
            "got {!r}".format(name_in, raw))
    return dim


def read_w90_header_marker(name_in):
    """Format marker on line 1 of a Wannier90-like file.

    Line 1 is a free comment in Wannier90; the spin-resolved Transfer format
    (issue #323) declares itself there with the token ``w90spin/<version>``.
    :func:`read_w90` reads the marker itself; this accessor exists so a
    reader that must record the format beside the table (the k-space input
    reader's ``source_formats``) can do so without re-parsing the body.

    Parameters
    ----------
    name_in : str
        Input filename

    Returns
    -------
    tuple or None
        ``(format, version)``, e.g. ``("w90spin", 1)``, or ``None`` when
        line 1 carries no marker (every 7-column file, genuine hr.dat
        included).
    """
    with open(name_in, 'r') as f:
        first = f.readline()
    m = _SPIN_COLUMN_MARKER.search(first)
    if m is None:
        return None
    return (SPIN_COLUMN_FORMAT, int(m.group(1)))


def _column_count_error(name_in, lineno, found, spin_format):
    if spin_format:
        expected = "{} (spin-resolved format {}/{} declared on line 1)".format(
            _NCOL_SPIN, SPIN_COLUMN_FORMAT, SPIN_COLUMN_VERSION)
    elif found == _NCOL_SPIN:
        return ValueError(
            "{}: line {} has {} columns but line 1 carries no '{}/{}' marker. "
            "A spin-resolved Transfer file must declare its format on line 1 "
            "(issue #323); an unmarked 9-column file is not guessed at."
            .format(name_in, lineno, found, SPIN_COLUMN_FORMAT,
                    SPIN_COLUMN_VERSION))
    else:
        expected = "{}".format(_NCOL_PLAIN)
    return ValueError("{}: line {} has {} columns, expected {}".format(
        name_in, lineno, found, expected))


def _fold_spin_columns(name_in, lineno, orbvec, s_str, t_str, norb_phys):
    """Spin-block key ``(s*norb + a, t*norb + b)`` of one spin-resolved row
    (``orbvec`` already 0-based), validating the ranges."""
    from hwave.solver.spin_orbital import spinblock_index
    for idx in orbvec:
        if not 0 <= idx < norb_phys:
            raise ValueError(
                "{}: line {}: orbital index {} is outside 1..{} (the norb of "
                "line 2; spin-resolved rows use physical orbital indices)"
                .format(name_in, lineno, idx + 1, norb_phys))
    spins = (int(s_str), int(t_str))
    for sp in spins:
        if sp not in (0, 1):
            raise ValueError(
                "{}: line {}: spin index {} is outside {{0, 1}}. H-wave reads "
                "the spin-resolved format for N_f = 2 flavors only "
                "(0 = up, 1 = down)".format(name_in, lineno, sp))
    return (spinblock_index(orbvec[0], spins[0], norb_phys),
            spinblock_index(orbvec[1], spins[1], norb_phys))


def read_w90(name_in, allow_spin_columns=False):
    """Read Wannier90 Hamiltonian data from file.

    Two body layouts are read. The 7-column form ``Rx Ry Rz a b Re Im`` is
    the Wannier90 ``hr.dat`` layout and H-wave's own interaction files. The
    spin-resolved 9-column form ``Rx Ry Rz a b s t Re Im`` (issue #323) is
    accepted only when ``allow_spin_columns`` is set and the file declares
    it with the ``w90spin/1`` token on line 1; its rows are folded to the
    spin-block composite index ``s*norb + a`` (``norb`` from line 2), so
    the returned keys have the same shape in both cases.
    
    Parameters
    ----------
    name_in : str
        Input filename
        
    Returns
    -------
    dict
        Dictionary mapping (irvec, orbvec) tuples to complex values where:
        - irvec : tuple
            Real space lattice vector indices
        - orbvec : tuple 
            Orbital indices
            
    allow_spin_columns : bool, optional
        Accept the spin-resolved 9-column format. Off by default: only the
        k-space input reader's Transfer entry opts in, so every other caller
        (two-body files, Extern, the direct reads in ``hwave_sc`` and
        ``flex_hf``) refuses such a file with an explicit error.

    Raises
    ------
    SystemExit
        If file not found, or on a duplicate / degeneracy inconsistency
        (the long-standing checks).
    ValueError
        On a column layout the format rules refuse: a spin-resolved marker
        without ``allow_spin_columns``, an unknown marker version, a row
        whose column count disagrees with the declared format (an unmarked
        9-column file included), a spin index outside {0, 1}, an orbital
        index outside the header ``norb``, or a non-unit degeneracy in the
        spin-resolved format.
    """
    try:
        with open(name_in, 'r') as f:
            lines = f.readlines()
    except OSError:
        logger.error("read_geom: file {} not found".format(name_in))
        exit(1)

    # line 1 is a free comment except for the format marker (issue #323)
    marker = _SPIN_COLUMN_MARKER.search(lines[0]) if lines else None
    l_strip = [s.strip() for s in lines[1:]]

    if marker is not None:
        version = int(marker.group(1))
        if not allow_spin_columns:
            raise ValueError(
                "{}: this file declares the spin-resolved (9-column) Transfer "
                "format {}/{}, which is accepted only for the Transfer entry "
                "of the k-space input reader (issue #323); it cannot be read "
                "here".format(name_in, SPIN_COLUMN_FORMAT, version))
        if version != SPIN_COLUMN_VERSION:
            raise ValueError(
                "{}: unknown spin-resolved format version {}/{} (this build "
                "reads {}/{})".format(name_in, SPIN_COLUMN_FORMAT, version,
                                      SPIN_COLUMN_FORMAT, SPIN_COLUMN_VERSION))
        # the spin-block fold needs the physical orbital count of line 2
        norb_phys = read_w90_header_dim(name_in)
        ncol = _NCOL_SPIN
    else:
        norb_phys = None
        ncol = _NCOL_PLAIN

    nr = int(l_strip[1])

    nints_per_line = 15
    skip_line = nr // nints_per_line
    if nr % nints_per_line != 0:
        skip_line += 1

    # Degeneracy (ndegen) block: nr integers, 15 per line. In the Wannier90
    # convention H(k) = sum_R exp(i k.R) H(R) / ndegen(R), so each matrix
    # element must be divided by the degeneracy of its R-point. H-wave's own
    # "wannier90-like" format is sparse (only nonzero entries are listed) and
    # always uses degeneracy 1, in which case this division is a no-op; for a
    # dense standard Wannier90 hr.dat the R-points appear in file order, so the
    # i-th distinct R-vector carries ndegen[i].
    ndegen = []
    for line in l_strip[2:2 + skip_line]:
        ndegen.extend(int(x) for x in line.split())
    if len(ndegen) != nr:
        logger.error(
            "read_w90: degeneracy count mismatch in {}: expected={}, found={}"
            .format(name_in, nr, len(ndegen)))
        exit(1)

    # When every degeneracy is 1 the division is a no-op, which covers H-wave's
    # own files; in that case nr may be a placeholder and the number of listed
    # R-points is left unconstrained. Only a genuine Wannier90 file (some
    # ndegen != 1) needs the positional R -> ndegen[i] mapping, which in turn
    # requires a dense listing.
    all_unit = all(d == 1 for d in ndegen)

    # The spin-resolved format is H-wave's own extension; the non-unit
    # degeneracy machinery below exists for genuine Wannier90 hr.dat files
    # and the extension stays out of it.
    if marker is not None and not all_unit:
        raise ValueError(
            "{}: the spin-resolved Transfer format requires every degeneracy "
            "to be 1 (found {})".format(
                name_in, sorted(set(d for d in ndegen if d != 1))))

    data = {}
    deg_of_r = {}  # irvec -> ndegen, by order of first appearance (non-unit only)
    ndup = 0
    body_first_line = 4 + skip_line     # 1-based file line of the first row
    for i, line in enumerate(l_strip[2 + skip_line:]):
        values = line.split()
        # if data is empty, break
        if len(values) == 0:
            break
        if len(values) != ncol:
            raise _column_count_error(name_in, body_first_line + i,
                                      len(values), marker is not None)
        irvec = (int(values[0]), int(values[1]), int(values[2]))
        orbvec = (int(values[3]) - 1, int(values[4]) - 1)
        if marker is not None:
            orbvec = _fold_spin_columns(name_in, body_first_line + i, orbvec,
                                        values[5], values[6], norb_phys)
            re_im = (values[7], values[8])
        else:
            re_im = (values[5], values[6])

        if all_unit:
            deg = 1
        else:
            # The i-th distinct R-vector (in file order) carries ndegen[i].
            if irvec not in deg_of_r:
                idx_r = len(deg_of_r)
                if idx_r >= len(ndegen):
                    logger.error(
                        "read_w90: more distinct R-points than declared (nr={}) in {}"
                        .format(nr, name_in))
                    exit(1)
                deg_of_r[irvec] = ndegen[idx_r]
            deg = deg_of_r[irvec]

        key = (irvec, orbvec)
        if key in data:
            ndup += 1
        data[key] = (float(re_im[0]) + 1J * float(re_im[1])) / deg

    # Reject duplicate (R, orbital) entries, matching the strictness of the
    # real-space reader (read_input.py). In the spin-resolved format the
    # key carries the spins, so the same (R, a, b) with different (s, t)
    # are distinct entries.
    if ndup > 0:
        logger.error("read_w90: duplicate ({}) entries found in {}".format(
            "R, a, b, s, t" if marker is not None else "R, orbital", name_in))
        exit(1)

    # With non-unit degeneracies the positional mapping requires every R-grid
    # point to be listed (dense file); reject a sparse listing rather than guess.
    if not all_unit and len(deg_of_r) != nr:
        logger.error(
            "read_w90: non-unit degeneracies require a dense file listing all "
            "{} R-points (found {}) in {}".format(nr, len(deg_of_r), name_in))
        exit(1)

    return data

def split_coulomb(data):
    """Split a combined 'Coulomb' table into intra and inter parts.

    The aggregate Coulomb input (zvo_ur.dat) provides both parts at once:
    the r=0 orbital-diagonal entries are the on-site CoulombIntra terms,
    everything else is CoulombInter.  This decomposition is the single
    definition shared by all solvers (UHFk and RPA/FLEX).

    Parameters
    ----------
    data : dict
        {((rx,ry,rz), (alpha,beta)): value} as returned by read_w90.

    Returns
    -------
    (dict, dict)
        (CoulombIntra part, CoulombInter part), same key format.
    """
    intra = {}
    inter = {}
    for (irvec, orbvec), v in data.items():
        if irvec == (0, 0, 0) and orbvec[0] == orbvec[1]:
            intra[(irvec, orbvec)] = v
        else:
            inter[(irvec, orbvec)] = v
    return intra, inter

def write_w90(name_in, info_int, info_geometry, interact_shape):
    """Write Wannier90 Hamiltonian data to file.
    
    Parameters
    ----------
    name_in : str
        Output filename
    info_int : dict
        Interaction information dictionary
    info_geometry : dict
        Geometry information dictionary
    interact_shape : tuple
        Shape of interaction array
    """
    # Build output content as a list of strings
    output = []
    
    # Header information
    norb = info_int["n_orb"]
    Nlattice = np.prod(interact_shape)
    n_total = (Nlattice * norb) ** 2
    
    output.append("# wannier90 format")
    output.append(str(norb))
    output.append(str(n_total))

    # Write list of ones, 15 per line
    ones = np.ones(n_total, dtype=int)
    for i in range(0, n_total, 15):
        output.append(" ".join(map(str, ones[i:i+15])))

    # Interaction data
    fmt = "{} {} {} {} {} {} {}"
    for idx, (site_org, int_value) in info_int["interaction"].items():
        site = info_geometry["site2vec"][idx]
        output.append(fmt.format(
            site_org,
            site[0], site[1], 
            site[2]+1, site[3]+1,
            int_value.real, int_value.imag
        ))

    # Write all output at once
    with open(name_in, "w") as fw:
        fw.write("\n".join(output) + "\n")

if __name__ == "__main__":
    path_to_sample = "../../sample/dir-model"
    import os
    data_geom = read_geom(os.path.join(path_to_sample, "zvo_geom.dat"))
    data_hr = read_w90(os.path.join(path_to_sample,"zvo_hr.dat"))
    print(data_geom)
    print(data_hr)
