import logging
import tomli
import qlmsio
import solver.uhf as sol_uhf
import solver.uhfk as sol_uhfk
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
    info_outputfile["path_to_output"] = info_outputfile.get("path_to_output", "output")
    path_to_output = info_outputfile["path_to_output"]

    logger = logging.getLogger("qlms")
    fmt = "%(asctime)s %(levelname)s %(name)s: %(message)s"
    # logging.basicConfig(level=logging.DEBUG, format=fmt)
    logging.basicConfig(level=logging.INFO, format=fmt)

    if "mode" not in info_mode:
        logger.error("mode is not defined in [mode].")
        exit(1)
    mode = info_mode["mode"]
    if mode == "UHF":
        logger.info("Read def files")
        read_io = qlmsio.read_input.QMLSInput(path_to_namelist)

        logger.info("Get Parameters information")
        mod_param_info = read_io.get_param("mod")
        pprint.pprint(mod_param_info, width = 1)

        logger.info("Get Hamiltonian information")
        ham_info = read_io.get_param("ham")

        logger.info("Get Output information")
        green_info = read_io.get_param("output")
        os.makedirs(path_to_output, exist_ok=True)

        solver = sol_uhf.UHF(ham_info, info_log, info_mode, mod_param_info)

    elif mode == "UHFk":
        logger.info("Read definitions from files")
        read_io = qlmsio.read_input_k.QMLSkInput(info_inputfile)

        # logger.info("Get parameter information")
        # mod_param_info = info_mode #read_io.get_param("mod")
        # pprint.pprint(mod_param_info, width = 1)

        logger.info("Get Hamiltonian information")
        ham_info = read_io.get_param("ham")

        logger.info("Get output information")
        green_info = read_io.get_param("output")
        os.makedirs(path_to_output, exist_ok=True)

        pprint.pprint(info_mode, width = 1)
        
        solver = sol_uhfk.UHFk(ham_info, info_log, info_mode)

    else:
        logger.warning("mode is incorrect: mode={}.".format(mode))
        exit(0)

    logger.info("Start UHF calculation")
    solver.solve(path_to_output)
    logger.info("Save calculation results.")
    solver.save_results(info_outputfile , green_info)
    logger.info("All procedures are finished.")

if __name__ == "__main__":
    main()
