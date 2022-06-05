import logging
import tomli
import qlmsio
import solver.uhf as sol_uhf
import sys
import os
import pprint

def main():
    args = sys.argv
    if len(args) != 2:
        print("Usage: python3 qlms.py input.toml")
        exit(1)
    with open(args[1], "rb") as f:
        toml_dict = tomli.load(f)

    #Initialize information about log
    info_log = toml_dict.get("log", {})
    info_log["print_level"] = info_log.get("print_level", 1)
    info_log["print_step"] = info_log.get("print_step", 1)

    #Initialize information about mode
    info_mode = toml_dict.get("mode", {})

    info_file = toml_dict.get("file", {"input": {}, "output": {}})
    #Initialize information about input files
    info_inputfile = info_file.get("input", {})
    info_inputfile["path_to_input"] = info_inputfile.get("path_to_input", "")
    info_inputfile["namelist"] = info_inputfile.get("namelist", "namelist.def")
    path_to_namelist = os.path.join(info_inputfile["path_to_input"], info_inputfile["namelist"])

    #Initialize information about output files
    info_outputfile = info_file.get("output",{})
    info_outputfile["path_to_output"] = info_file.get("path_to_output", "output")
    path_to_output = info_outputfile["path_to_output"]

    logger = logging.getLogger("qlms")
    fmt = "%(asctime)s %(levelname)s %(name)s :%(message)s"
    # logging.basicConfig(level=logging.DEBUG, format=fmt)
    logging.basicConfig(level=logging.INFO, format=fmt)

    logger.info("Read def files")
    read_io = qlmsio.read_input.QMLSInput(path_to_namelist)
    logger.info("Get Parameters information")
    param_info = read_io.get_param("param")
    pprint.pprint(param_info, width = 1)
    logger.info("Get Hamiltonian information")
    ham_info = read_io.get_param("ham")
    logger.info("Get Output information")
    green_info = read_io.get_param("output")
    os.makedirs(path_to_output, exist_ok=True)

    logger.info("Start UHF calculation")
    uhf = sol_uhf.UHF(param_info, ham_info, info_log, info_mode)
    uhf.solve(path_to_output)

    logger.info("Save calculation results.")
    uhf.save_results(info_outputfile , green_info)
    logger.info("All procedures are finished.")

if __name__ == "__main__":
    main()
