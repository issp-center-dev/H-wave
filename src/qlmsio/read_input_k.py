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
        self.ham_param["Initial"] = info_inputfile.get("initial", None)


    def get_param(self, key):
        if key == "ham":
            return self.ham_param
        else:
            # Add error message
            logger.error("Get_param: key must be mod or ham or output.")
            return None

if __name__ == '__main__':
    qml_input = QMLSInput("namelist.def")
