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
read_w90(name_in)
    Read Wannier90 Hamiltonian data from file.
write_w90(name_in, info_int, info_geometry, interact_shape)
    Write Wannier90 Hamiltonian data to file.
"""

from __future__ import print_function

import itertools
import numpy as np
import logging

logger = logging.getLogger("qlms").getChild("wan90")

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

def read_w90(name_in):
    """Read Wannier90 Hamiltonian data from file.
    
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
            
    Raises
    ------
    SystemExit
        If file not found
    """
    try:
        with open(name_in, 'r') as f:
            # skip header
            l_strip = [s.strip() for s in f.readlines()[1:]]
    except OSError:
        logger.error("read_geom: file {} not found".format(name_in))
        exit(1)

    nr = int(l_strip[1])

    nints_per_line = 15
    skip_line = nr // nints_per_line
    if nr % nints_per_line != 0:
        skip_line += 1

    mat_size = len(l_strip[2 + skip_line:])
    data = {}
    for idx, line in enumerate(l_strip[2 + skip_line:]):
        values = line.split()
        # if data is empty, break
        if len(values) == 0:
            break
        irvec = (int(values[0]), int(values[1]), int(values[2]))
        orbvec = (int(values[3]) - 1, int(values[4]) - 1)
        data[(irvec, orbvec)] = float(values[5]) + 1J * float(values[6])

    return data

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
