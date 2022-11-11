import logging
import numpy as np
import sys
import os
import re
from requests.structures import CaseInsensitiveDict
from . import wan90

logger = logging.getLogger("qlms").getChild("read_input")


class QMLSkInput():
    def __init__(self, info_inputfile, solver_type="UHFk"):
        logger.info(">>> QMLSkInput init")

        # [file.input]
        #   path_to_input
        #   initial
        # [file.input.interaction]
        #   path_to_input
        #   Geometry
        #   Transfer
        #   (and other interaction terms)

        self.ham_param = CaseInsensitiveDict()
        self.initial = CaseInsensitiveDict()

        files = info_inputfile.get("interaction", {})

        input_dir = files.get("path_to_input", "")

        for k, v in files.items():
            if k == "path_to_input":
                pass
            elif k == "Geometry":
                f = os.path.join(input_dir, v)
                logger.info("QMLSkInput: read Gemoetry from {}".format(f))
                self.ham_param[k] = wan90.read_geom(f)
            else:
                f = os.path.join(input_dir, v)
                logger.info("QMLSkInput: read interaction {} from {}".format(k, f))
                self.ham_param[k] = wan90.read_w90(f)

        # filename for initial green function
        if "initial" in info_inputfile:
            self.initial["uhfk"] = info_inputfile["initial"]

        # initial green function in coordinate space
        if "initial_uhf" in info_inputfile:
            if "initial" in info_inputfile:
                logger.error("initial and initial_uhf are not specified simultaneously.")
                exit(1)
            else:
                f = os.path.join(input_dir, info_inputfile["initial_uhf"])
                logger.info("QMLSkInput: read initial_uhf from {}".format(f))
                self.initial["uhf"] = self._read_data(f, "complex")
        if "geometry_uhf" in info_inputfile:
            f = os.path.join(input_dir, info_inputfile["geometry_uhf"])
            logger.info("QMLSkInput: read geometry_uhf from {}".format(f))
            self.initial["geometry"] = wan90.read_geometry(f)
        if "greenone_uhf" in info_inputfile:
            f = os.path.join(input_dir, info_inputfile["greenone_uhf"])
            logger.info("QMLSkInput: read greenone_uhf from {}".format(f))
            self.initial["greenone"] = self._read_green(f)

        if len(self.initial) > 0:
            self.ham_param["Initial"] = self.initial
        else:
            self.ham_param["Initial"] = None

    def _read_data(self, file_name, value_type="real"):
        info = {}
        data = np.loadtxt(file_name, skiprows = 5)
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
        return np.loadtxt(file_name, dtype=np.int32, skiprows = 5)

    def get_param(self, key):
        if key == "ham":
            return self.ham_param
        else:
            logger.error("Get_param: key must be mod or ham or output.")
            return None

if __name__ == '__main__':
    qml_input = QMLSkInput({ "Geometry": "geometry.dat" })
