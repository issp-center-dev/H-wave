import logging

import qlmsio
import solver.uhf as sol_uhf

def main():
    logger = logging.getLogger("qlms")
    fmt = "%(asctime)s %(levelname)s %(name)s :%(message)s"
    # logging.basicConfig(level=logging.DEBUG, format=fmt)
    logging.basicConfig(level=logging.INFO, format=fmt)

    logger.info("Read def files")
    read_io = qlmsio.read_input.QMLSInput("namelist.def")
    logger.info("Get Prameters information")
    param_info = read_io.get_param("param")
    logger.info("Get Hamiltonian information")
    ham_info = read_io.get_param("ham")
    logger.info("Get Output information")
    output_info = read_io.get_param("output")

    logger.info("Start UHF calculation")
    uhf = sol_uhf.UHF(param_info, ham_info, output_info)
    uhf.solve()


if __name__ == "__main__":
    main()
