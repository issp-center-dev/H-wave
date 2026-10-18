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
import hwave.solver.flex as sol_flex
from requests.structures import CaseInsensitiveDict

#: the keys read from the [mode] table itself (every other solver setting
#: belongs in [mode.param]). [mode] is a plain TOML table, so these names
#: are matched exactly.
MODE_TABLE_KEYS = ("mode", "param", "calc_scheme", "calc_type",
                   "enable_spin_orbital", "flag_fock")


def check_mode_tables(info_mode):
    """Messages for keys of ``[mode]`` / ``[mode.param]`` that the solvers
    would silently ignore: a ``[mode]`` key placed in ``[mode.param]``, an
    unknown key in ``[mode]`` (a solver parameter belongs in
    ``[mode.param]``), or a ``[mode]`` key spelled with different case.
    An empty list means the tables are well placed."""
    from collections.abc import Mapping
    errors = []
    by_lower = {k.lower(): k for k in MODE_TABLE_KEYS}
    for key in info_mode:
        if key in MODE_TABLE_KEYS:
            continue
        if key.lower() in by_lower:
            errors.append("[mode] {!r} is not read: keys of [mode] are case-sensitive, "
                          "write it as {!r}.".format(key, by_lower[key.lower()]))
        else:
            errors.append("[mode] {!r} is not a [mode] key and would be ignored; solver "
                          "parameters belong in [mode.param].".format(key))
    param = info_mode.get("param", {})
    if isinstance(param, Mapping):
        for key in param:
            if key.lower() == "param":
                errors.append("[mode.param] {!r}: [mode.param] has no nested param "
                              "table.".format(key))
            elif key.lower() in by_lower:
                errors.append("[mode.param] {!r} is read from [mode], not from [mode.param]; "
                              "move it to [mode] as {!r}.".format(key, by_lower[key.lower()]))
    return errors


def report_unused_params(info_mode, solver, logger):
    """Warn once about the ``[mode.param]`` keys this run never read
    (a misspelled or misplaced parameter has no effect otherwise)."""
    from hwave.solver.base import ParamDict
    param = info_mode.get("param")
    if not isinstance(param, ParamDict):
        return
    unused = [k for k, _ in param.items() if k.lower() not in param.accessed]
    if unused:
        logger.warning("[mode.param] key(s) not read by this {} run: {}. A misspelled or "
                       "misplaced key has no effect; keys that only apply to other "
                       "settings of the solver are listed as well.".format(
                           info_mode.get("mode"), ", ".join(sorted(unused))))


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

    # A key in the wrong table, or a [mode] key spelled with other case, used
    # to be ignored without a word (GitHub issue #236).
    placement_errors = check_mode_tables(info_mode)
    if placement_errors:
        for msg in placement_errors:
            logger.error(msg)
        exit(1)
    # record which [mode.param] keys the solver reads (reported at the end)
    from hwave.solver.base import ParamDict
    info_mode = dict(info_mode)
    info_mode["param"] = ParamDict(info_mode.get("param", {}) or {})

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

    elif mode == "FLEX":
        # Fluctuation Exchange approximation (self-consistent RPA)
        logger.info("FLEX mode")

        logger.info("Read interaction definitions from files")
        read_io = qlmsio.read_input_k.QLMSkInput(info_inputfile)
        ham_info = read_io.get_param("ham")

        # the OPTIONAL [eliashberg] table of the same input file configures the
        # in-process pairing step ([mode.param] longitudinal_bond_pairing);
        # without that key a present table is ignored with an INFO line, so a
        # single input file can serve both hwave and hwave_sc.
        solver = sol_flex.FLEX(ham_info, info_log, info_mode,
                               eliashberg_param=input_dict.get("eliashberg"))

        green_info = read_io.get_param("green")
        green_info.update(solver.read_init(info_inputfile))

    else:
        logger.warning("mode is incorrect: mode={}.".format(mode))
        exit(0)

    # Execute calculation
    if mode == "FLEX":
        # output-path validation before any work (a collision between two
        # artifacts is refused here, at solve entry and again before writing)
        solver.validate_output_paths(info_outputfile, path_to_output)
    logger.info("Start UHF calculation")
    solver.solve(green_info, path_to_output)
    logger.info("Save calculation results.")
    solver.save_results(info_outputfile, green_info)
    report_unused_params(info_mode, solver, logger)
    logger.info("All procedures are finished.")

    if mode == "FLEX":
        return {
            "scf_converged": bool(getattr(solver, "scf_converged", False)),
            "scf_iterations": int(getattr(solver, "scf_iterations", -1)),
        }
    return {}


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
