import logging
import numpy as np
import sys
import os
import re
from requests.structures import CaseInsensitiveDict
from . import wan90

logger = logging.getLogger("qlms").getChild("read_input")


class QLMSkInput():
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
        try:
            _data = np.loadtxt(file_name, dtype=np.int32, skiprows = 5)
        except FileNotFoundError:
            logger.error("read_input_k: file not found: {}".format(file_name))
            exit(1)
        return _data

    def get_param(self, key):
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
