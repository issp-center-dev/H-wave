"""Input file reader for k-space calculations.

This module provides functionality to read and parse input files for k-space 
calculations in QLMS. It handles geometry, transfer integrals, and various 
interaction terms.

Classes
-------
QLMSkInput
    Reads and parses input files for k-space calculations.
"""

import logging
import numpy as np
import sys
import os
import re
from requests.structures import CaseInsensitiveDict
from . import wan90

logger = logging.getLogger("qlms").getChild("read_input")


class QLMSkInput():
    """Input file reader for k-space calculations.
    
    Parameters
    ----------
    info_inputfile : dict
        Dictionary containing input file information with keys:
        - path_to_input : str
            Path to input file directory
        - interaction : dict
            Dictionary specifying interaction files
        - initial : str, optional
            Filename for initial Green's function in k-space (NPZ format)
        - initial_uhf : str, optional
            Filename for initial Green's function in real space
        - geometry_uhf : str, optional
            Filename for geometry in real space
        - onebodyg_uhf : str, optional
            Filename for one-body Green's function indices
    solver_type : str, optional
        Type of solver to use (default: "UHFk")
        
    Attributes
    ----------
    valid_namelist : list
        List of valid input file section names
    ham_param : CaseInsensitiveDict
        Dictionary storing Hamiltonian parameters
    green : CaseInsensitiveDict
        Dictionary storing Green's function data
        
    Methods
    -------
    get_param(key)
        Get parameters by key
    _read_data(file_name, value_type)
        Read data from file into dictionary
    _read_green(file_name) 
        Read Green's function indices from file
    """

    valid_namelist = [s.lower() for s in ["path_to_input", "Geometry", "Transfer", "CoulombIntra", "CoulombInter", "Hund", "Ising", "PairLift", "Exchange", "PairHop", "Extern"]]

    def __init__(self, info_inputfile, solver_type="UHFk"):
        logger.debug(">>> QLMSkInput init")

        # [file.input]
        #   path_to_input
        #   initial
        # [file.input.interaction]
        #   path_to_input
        #   Geometry
        #   Transfer
        #   (and other interaction terms)

        self.ham_param = CaseInsensitiveDict()
        self.green = CaseInsensitiveDict()

        # file.input.interaction
        files = info_inputfile.get("interaction", {})
        interaction_file_dir = files.get("path_to_input", "")

        # check if keyword is valid
        err = 0
        for k, v in files.items():
            if k.lower() not in self.valid_namelist:
                logger.error("Unknown keyword {}".format(k))
                err += 1
        if err > 0:
            logger.fatal("Invalid input.")
            exit(1)

        for k, v in files.items():
            if k == "path_to_input":
                pass
            elif k == "Geometry":
                f = os.path.join(interaction_file_dir, v)
                logger.info("QLMSkInput: read Gemoetry from {}".format(f))
                self.ham_param[k] = wan90.read_geom(f)
            else:
                f = os.path.join(interaction_file_dir, v)
                logger.info("QLMSkInput: read interaction {} from {}".format(k, f))
                self.ham_param[k] = wan90.read_w90(f)

        # file.input
        input_file_dir = info_inputfile.get("path_to_input", "")

        #- initial green function in wave number space: npz format
        if "initial" in info_inputfile:
            file_name = os.path.join(input_file_dir, info_inputfile["initial"])
            logger.info("QLMSkInput: read initial from {}".format(file_name))
            try:
                _data = np.load(file_name)
            except OSError:
                logger.error("read_input_k: file {} not found".format(file_name))
                _data = None
                exit(1)
            self.green["initial"] = _data

        #- initial green function in coordinate space
        if "initial_uhf" in info_inputfile:
            if "initial" in info_inputfile:
                logger.error("initial and initial_uhf can not be specified simultaneously.")
                exit(1)
            else:
                file_name = os.path.join(input_file_dir, info_inputfile["initial_uhf"])
                logger.info("QLMSkInput: read initial_uhf from {}".format(file_name))
                self.green["initial_uhf"] = self._read_data(file_name, "complex")

        #- geometry.dat is required to handle data in coordinate space
        if "geometry_uhf" in info_inputfile:
            file_name = os.path.join(input_file_dir, info_inputfile["geometry_uhf"])
            logger.info("QLMSkInput: read geometry_uhf from {}".format(file_name))
            self.green["geometry_uhf"] = wan90.read_geometry(file_name)

        #- index list for output green function in coordinate space
        if "onebodyg_uhf" in info_inputfile:
            file_name = os.path.join(input_file_dir, info_inputfile["onebodyg_uhf"])
            logger.info("QLMSkInput: read greenone_uhf from {}".format(file_name))
            self.green["onebodyg_uhf"] = self._read_green(file_name)

    def _read_data(self, file_name, value_type="real"):
        """Read data from file into dictionary.
        
        Parameters
        ----------
        file_name : str
            Name of file to read
        value_type : str, optional
            Type of values to read ("real" or "complex")
            
        Returns
        -------
        dict
            Dictionary containing data read from file
            
        Raises
        ------
        FileNotFoundError
            If input file is not found
        """
        info = {}
        try:
            data = np.loadtxt(file_name, skiprows = 5)
        except FileNotFoundError:
            logger.error("read_input_k: file not found: {}".format(file_name))
            exit(1)
        if value_type == "real":
            for _data in data:
                list = tuple([int(i) for i in _data[:-1]])
                value = float(_data[-1])
                info[list] = value
        elif value_type == "complex":
            for _data in data:
                list = tuple([int(i) for i in _data[:-2]])
                value = float(_data[-2]) + 1j * float(_data[-1])
                info[list] = value
        return info

    def _read_green(self, file_name):
        """Read Green's function indices from file.
        
        Parameters
        ----------
        file_name : str
            Name of file to read
            
        Returns
        -------
        ndarray
            Array of Green's function indices
            
        Raises
        ------
        FileNotFoundError
            If input file is not found
        """
        try:
            _data = np.loadtxt(file_name, dtype=np.int32, skiprows = 5)
        except FileNotFoundError:
            logger.error("read_input_k: file not found: {}".format(file_name))
            exit(1)
        return _data

    def get_param(self, key):
        """Get parameters by key.
        
        Parameters
        ----------
        key : str
            Key to retrieve parameters for:
            - "mod"/"parameter": Returns None
            - "ham"/"hamiltonian": Returns Hamiltonian parameters
            - "output"/"green": Returns Green's function data
            
        Returns
        -------
        CaseInsensitiveDict or None
            Requested parameters or None if key is invalid
        """
        if key == "mod" or key == "parameter":
            return None
        elif key == "ham" or key == "hamiltonian":
            return self.ham_param
        elif key == "output" or key == "green":
            return self.green
        else:
            logger.error("Get_param: key must be mod or ham or output.")
            return None

if __name__ == '__main__':
    qml_input = QLMSkInput({ "Geometry": "geometry.dat" })
