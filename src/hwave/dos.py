"""Density of states (DoS) calculation and plotting utilities.

This module provides functionality for calculating, plotting and analyzing the
density of states from Hartree-Fock calculations.

"""

from __future__ import annotations

import itertools
import os

import numpy as np

import hwave


class DoS:
    """Class for storing and manipulating density of states data.
    
    Parameters
    ----------
    ene : np.ndarray
        Energy grid points
    dos : np.ndarray 
        Density of states values for each orbital at each energy point
        
    Attributes
    ----------
    dos : np.ndarray
        Density of states array with shape (norb, nene)
    ene : np.ndarray
        Energy grid points array with shape (nene,)
    ene_num : int
        Number of energy points
    norb : int
        Number of orbitals
    """

    def __init__(self, ene: np.ndarray, dos: np.ndarray):
        assert ene.shape[0] == dos.shape[1]
        self.ene = ene
        self.dos = dos
        self.ene_num = ene.shape[0]
        self.norb = dos.shape[0]

    def plot(self, filename: str = "", verbose: bool = False):
        """Plot the density of states.

        Creates a plot showing the total DoS and orbital-resolved DoS.

        Parameters
        ----------
        filename : str, optional
            If provided, save plot to this file
        verbose : bool, optional
            If True, print additional output
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError(
                "matplotlib is not installed. Please install matplotlib and try again."
            )

        total_dos = np.sum(self.dos, axis=0)

        plt.plot(self.ene, total_dos, label="Total")
        for i in range(self.norb):
            plt.plot(self.ene, self.dos[i], label=str(i))
        plt.xlabel("Energy")
        plt.ylabel("DOS")
        plt.ylim(0)
        plt.legend()
        if filename != "":
            if verbose:
                print("Saving plot to file: ", filename)
            plt.savefig(filename)
        plt.close()

    def write_dos(self, output: str, verbose: bool = False):
        """Write density of states data to file.

        Parameters
        ----------
        output : str
            Output filename
        verbose : bool, optional
            If True, print additional output
        """
        if verbose:
            print("Writing DOS to file: ", output)
        total_dos = np.sum(self.dos, axis=0)
        with open(output, "w") as fw:
            fw.write("# ene total_dos")
            for j in range(self.norb):
                fw.write(f" dos[iorb={j}]")
            fw.write("\n")
            for i in range(self.ene_num):
                fw.write("{:15.8f} ".format(self.ene[i]))
                fw.write("{:15.8f} ".format(total_dos[i]))
                for j in range(self.norb):
                    fw.write("{:15.8f} ".format(self.dos[j, i]))
                fw.write("\n")


def __read_geom(file_name="./dir-model/zvo_geom.dat"):
    """Read geometry data from file.

    Parameters
    ----------
    file_name : str, optional
        Path to geometry file

    Returns
    -------
    np.ndarray
        Unit cell vectors array with shape (3,3)
    """
    with open(file_name, "r") as fr:
        uvec = np.zeros((3, 3))
        for i, line in enumerate(itertools.islice(fr, 3)):  # take first 3 lines
            uvec[i, :] = np.array(line.split())
    return uvec


def calc_dos(
    input_dict: dict,
    ene_window: list | None = None,
    ene_num: int = 101,
    verbose: bool = False,
) -> DoS:
    """Calculate density of states.

    Parameters
    ----------
    input_dict : dict
        Input parameters dictionary
    ene_window : list, optional
        Energy window [min, max] for DoS calculation
    ene_num : int, optional
        Number of energy points
    verbose : bool, optional
        If True, print additional output

    Returns
    -------
    DoS
        Calculated density of states object

    Raises
    ------
    ImportError
        If required libtetrabz package is not installed
    """
    try:
        import libtetrabz
    except ImportError:
        raise ImportError(
            "libtetrabz is not installed. Please install libtetrabz and try again."
        )

    output_info_dict = input_dict["file"]["output"]
    filename = os.path.join(
        output_info_dict["path_to_output"], output_info_dict["eigen"] + ".npz"
    )
    if verbose:
        print("Reading eigenvalues from file: ", filename)
    data = np.load(os.path.join(filename))
    eigenvalues = data["eigenvalue"]
    Lx, Ly, Lz = input_dict["mode"]["param"]["CellShape"]
    Lxsub, Lysub, Lzsub = input_dict["mode"]["param"].get("SubShape", (Lx,Ly,Lz))
    norb = eigenvalues.shape[1]
    if verbose:
        print("Lx, Ly, Lz, norb: ", Lx, Ly, Lz, norb)
        print("Lxsub, Lysub, Lzsub, norb: ", Lxsub, Lysub, Lzsub, norb)

    input_info_dict = input_dict["file"]["input"]["interaction"]
    file_name = os.path.join(
        input_info_dict["path_to_input"], input_info_dict["Geometry"]
    )
    if verbose:
        print("Reading geometry from file: ", file_name)
    uvec = __read_geom(file_name)
    bvec = 2.0 * np.pi * np.linalg.inv(uvec).T

    ene_min = np.min(eigenvalues)
    ene_max = np.max(eigenvalues)
    if verbose:
        print("Calculated energy min, max: ", ene_min, ene_max)

    if ene_window is None:
        ene_min = np.min(eigenvalues) - 0.2
        ene_max = np.max(eigenvalues) + 0.2
    else:
        ene_min = ene_window[0]
        ene_max = ene_window[1]

    ene = np.linspace(ene_min, ene_max, num=ene_num)
    if verbose:
        print("Energy window min, max, num: ", ene_min, ene_max, ene_num)

    eig = eigenvalues.reshape(int(Lx/Lxsub), int(Ly/Lysub), int(Lz/Lzsub), norb)

    if verbose:
        print("Finish calculating DOS")
    wght = libtetrabz.dos(bvec, eig, ene)
    dos = wght.sum(2).sum(1).sum(0)
    return DoS(dos=dos, ene=ene)


def main():
    """Command-line interface for DoS calculation.

    Parses command line arguments and runs DoS calculation.
    """
    import tomli
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("input", type=str, help="input file of hwave")
    parser.add_argument("-o","--output", type=str, default="dos.dat", help="DoS output")
    parser.add_argument(
        "--ene-window",
        default=None,
        type=float,
        nargs=2,
        help="""energy window; [ene_low, ene_high].
If omitted, [ene_min - 0.2, ene_max + 0.2]""",
    )
    parser.add_argument(
        "--ene-num", default=101, type=int, help="number of energy points"
    )
    parser.add_argument("-p", "--plot", type=str, default="", help="plot DOS to file")
    parser.add_argument("-q", "--quiet", action="store_true", help="calculate quietly")
    parser.add_argument(
        "-v", "--version", action="version", version=f"hwave_dos v{hwave.__version__}"
    )
    args = parser.parse_args()
    verbose = not args.quiet
    file_toml = args.input
    if os.path.exists(file_toml):
        if verbose:
            print("Reading input file: ", file_toml)
        with open(file_toml, "rb") as f:
            input_dict = tomli.load(f)
    else:
        raise ValueError("Input file does not exist")

    dos = calc_dos(
        input_dict,
        ene_window=args.ene_window,
        ene_num=args.ene_num,
        verbose=verbose,
    )
    output_dir = input_dict["file"]["output"]["path_to_output"]
    if args.output != "":
        dos.write_dos(os.path.join(output_dir, args.output), verbose=verbose)
    if args.plot:
        dos.plot(os.path.join(output_dir, args.plot), verbose=verbose)
