from __future__ import annotations
from typing import Optional

import sys
import os
# import pprint
import logging

import tomli

import hwave.qlmsio as qlmsio
import hwave.solver.uhfr as sol_uhfr
import hwave.solver.uhfk as sol_uhfk
import hwave.solver.rpa as sol_rpa
from requests.structures import CaseInsensitiveDict

def run(*, input_dict: Optional[dict] = None, input_file: Optional[str] = None):
    if input_dict is None:
        if input_file is None:
            raise RuntimeError("Neither input_dict nor input_file are passed")
        with open(input_file, "rb") as f:
            input_dict = tomli.load(f)
    else:
        if input_file is not None:
            raise RuntimeError("Both input_dict and input_file are passed")

    # Initialize information about log
    info_log = input_dict.get("log", {})
    info_log["print_level"] = info_log.get("print_level", 1)
    info_log["print_step"] = info_log.get("print_step", 1)

    # Initialize information about mode
    info_mode = input_dict.get("mode", {})
    info_file = input_dict.get("file", {"input": {}, "output": {}})
    # Initialize information about input files
    info_inputfile = info_file.get("input", {})
    info_inputfile["path_to_input"] = info_inputfile.get("path_to_input", "")

    # Initialize information about output files
    info_outputfile = info_file.get("output", {})
    info_outputfile["path_to_output"] = info_outputfile.get("path_to_output", "output")
    path_to_output = info_outputfile["path_to_output"]
    os.makedirs(path_to_output, exist_ok=True)

    logger = logging.getLogger("qlms")
    fmt = "%(asctime)s %(levelname)s %(name)s: %(message)s"
    # logging.basicConfig(level=logging.DEBUG, format=fmt)
    logging.basicConfig(level=logging.INFO, format=fmt)

    if "mode" not in info_mode:
        logger.error("mode is not defined in [mode].")
        exit(1)
    mode = info_mode["mode"]
    if mode == "UHFr":
        logger.info("Read def files")
        file_list = CaseInsensitiveDict()
        #interaction files
        for key, file_name in info_inputfile["interaction"].items():
            if key.lower() in ["trans", "coulombinter", "coulombintra", "pairhop", "hund", "exchange", "ising", "pairlift", "interall"]:
                file_list[key] = os.path.join(info_inputfile["path_to_input"], file_name)
            else:
                logging.error("Keyword {} is incorrect.".format(key))
                exit(1)
        #initial and green
        for key, file_name in info_inputfile.items():
            if key.lower() == "initial":
                file_list[key] = os.path.join(info_inputfile["path_to_input"], file_name)
            elif key.lower() == "onebodyg":
                file_list[key] = os.path.join(info_inputfile["path_to_input"], file_name)
        read_io = qlmsio.read_input.QLMSInput(file_list)

        logger.info("Get Hamiltonian information")
        ham_info = read_io.get_param("ham")

        logger.info("Get Green function information")
        green_info = read_io.get_param("green")

        # solver = sol_uhf.UHF(ham_info, info_log, info_mode, mod_param_info)
        solver = sol_uhfr.UHFr(ham_info, info_log, info_mode)

    elif mode == "UHFk":
        logger.info("Read definitions from files")
        read_io = qlmsio.read_input_k.QLMSkInput(info_inputfile)

        # logger.info("Get parameter information")
        # mod_param_info = info_mode #read_io.get_param("mod")
        # pprint.pprint(mod_param_info, width = 1)

        logger.info("Get Hamiltonian information")
        ham_info = read_io.get_param("ham")

        logger.info("Get Green function information")
        green_info = read_io.get_param("green")

        # pprint.pprint(info_mode, width=1)

        solver = sol_uhfk.UHFk(ham_info, info_log, info_mode)

    elif mode == "RPA":
        logger.info("RPA mode")

        logger.info("Read interaction definitions from files")
        read_io = qlmsio.read_input_k.QLMSkInput(info_inputfile)
        ham_info = read_io.get_param("ham")

        solver = sol_rpa.RPA(ham_info, info_log, info_mode)

        green_info = {}
        # logger.info("Start RPA calculation")
        # solver.solve(green_info, path_to_output)

        # logger.info("Save calculation results")
        # solver.save_results(info_outputfile, green_info)

        # logger.info("all done.")

    else:
        logger.warning("mode is incorrect: mode={}.".format(mode))
        exit(0)

    logger.info("Start UHF calculation")
    solver.solve(green_info, path_to_output)
    logger.info("Save calculation results.")
    solver.save_results(info_outputfile, green_info)
    logger.info("All procedures are finished.")


def main():
    args = sys.argv
    if len(args) != 2:
        print("Usage: python3 qlms.py input.toml")
        sys.exit(1)
    run(input_file=args[1])

