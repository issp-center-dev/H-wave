from __future__ import annotations
from typing import Optional

import sys
import os
# import pprint
import logging

import tomli

import hwave
import hwave.qlmsio as qlmsio
import hwave.solver.uhfr as sol_uhfr
import hwave.solver.uhfk as sol_uhfk
import hwave.solver.rpa as sol_rpa
from requests.structures import CaseInsensitiveDict

def run(*, input_dict: Optional[dict] = None, input_file: Optional[str] = None):
    """Run Hartree-Fock calculation with given input parameters.
    
    Parameters
    ----------
    input_dict : dict, optional
        Dictionary containing input parameters. Cannot be used with input_file.
    input_file : str, optional
        Path to TOML input file. Cannot be used with input_dict.
        
    Raises
    ------
    RuntimeError
        If neither or both input_dict and input_file are provided
    """
    # Validate input arguments
    if input_dict is None:
        if input_file is None:
            raise RuntimeError("Neither input_dict nor input_file are passed")
        # Load parameters from TOML file
        with open(input_file, "rb") as f:
            input_dict = tomli.load(f)
    else:
        if input_file is not None:
            raise RuntimeError("Both input_dict and input_file are passed")

    # Initialize logging configuration
    info_log = input_dict.get("log", {})
    info_log["print_level"] = info_log.get("print_level", 1)  # Default print level
    info_log["print_step"] = info_log.get("print_step", 1)    # Default print step

    # Get calculation mode and file paths
    info_mode = input_dict.get("mode", {})
    info_file = input_dict.get("file", {"input": {}, "output": {}})
    
    # Setup input file paths
    info_inputfile = info_file.get("input", {})
    info_inputfile["path_to_input"] = info_inputfile.get("path_to_input", "")

    # Setup output directory
    info_outputfile = info_file.get("output", {})
    info_outputfile["path_to_output"] = info_outputfile.get("path_to_output", "output")
    path_to_output = info_outputfile["path_to_output"]
    os.makedirs(path_to_output, exist_ok=True)

    # Configure logging
    logger = logging.getLogger("qlms")
    fmt = "%(asctime)s %(levelname)s %(name)s: %(message)s"
    logging.basicConfig(level=logging.INFO, format=fmt)

    # Validate calculation mode
    if "mode" not in info_mode:
        logger.error("mode is not defined in [mode].")
        exit(1)
    mode = info_mode["mode"]

    # Initialize solver based on calculation mode
    if mode == "UHFr":
        # Real-space unrestricted Hartree-Fock
        logger.info("Read def files")
        file_list = CaseInsensitiveDict()
        
        # Process interaction file paths
        for key, file_name in info_inputfile["interaction"].items():
            if key.lower() in ["trans", "coulombinter", "coulombintra", "pairhop", 
                             "hund", "exchange", "ising", "pairlift", "interall"]:
                file_list[key] = os.path.join(info_inputfile["path_to_input"], file_name)
            else:
                logging.error("Keyword {} is incorrect.".format(key))
                exit(1)
                
        # Process initial state and Green's function files
        for key, file_name in info_inputfile.items():
            if key.lower() == "initial":
                file_list[key] = os.path.join(info_inputfile["path_to_input"], file_name)
            elif key.lower() == "onebodyg":
                file_list[key] = os.path.join(info_inputfile["path_to_input"], file_name)
        read_io = qlmsio.read_input.QLMSInput(file_list)

        # Read Hamiltonian and Green's function parameters
        logger.info("Get Hamiltonian information")
        ham_info = read_io.get_param("ham")

        logger.info("Get Green function information")
        green_info = read_io.get_param("green")

        solver = sol_uhfr.UHFr(ham_info, info_log, info_mode)

    elif mode == "UHFk":
        # k-space unrestricted Hartree-Fock
        logger.info("Read definitions from files")
        read_io = qlmsio.read_input_k.QLMSkInput(info_inputfile)

        logger.info("Get Hamiltonian information")
        ham_info = read_io.get_param("ham")

        logger.info("Get Green function information")
        green_info = read_io.get_param("green")

        solver = sol_uhfk.UHFk(ham_info, info_log, info_mode)

    elif mode == "RPA":
        # Random Phase Approximation
        logger.info("RPA mode")

        logger.info("Read interaction definitions from files")
        read_io = qlmsio.read_input_k.QLMSkInput(info_inputfile)
        ham_info = read_io.get_param("ham")

        solver = sol_rpa.RPA(ham_info, info_log, info_mode)

        green_info = read_io.get_param("green")
        green_info.update(solver.read_init(info_inputfile))

    else:
        logger.warning("mode is incorrect: mode={}.".format(mode))
        exit(0)

    # Execute calculation
    logger.info("Start UHF calculation")
    solver.solve(green_info, path_to_output)
    logger.info("Save calculation results.")
    solver.save_results(info_outputfile, green_info)
    logger.info("All procedures are finished.")


def main():
    """Command-line interface entry point.
    
    Parses command line arguments and runs the calculation.
    """
    import argparse

    # Setup argument parser
    parser = argparse.ArgumentParser(prog='hwave')
    parser.add_argument('input_toml', nargs='?', default=None, help='input parameter file')
    parser.add_argument('--version', action='store_true', help='show version')

    args = parser.parse_args()

    # Handle version request
    if args.version:
        print('hwave', hwave.__version__)
        sys.exit(0)

    # Validate input file
    if args.input_toml is None:
        parser.print_help()
        sys.exit(1)

    # Run calculation
    run(input_file = args.input_toml)
