"""Density of states (DoS) calculation and plotting utilities.

This module provides functionality for calculating, plotting and analyzing the
density of states from Hartree-Fock calculations.

"""

from __future__ import annotations

import itertools
import logging
import os
import re

import numpy as np

import hwave

logger = logging.getLogger("hwave.dos")


class DoS:
    """Class for storing and manipulating density of states data.
    
    Parameters
    ----------
    ene : np.ndarray
        Energy grid points
    dos : np.ndarray
        Density of states values for each band at each energy point.
        This is band-resolved (one row per eigenvalue/band of the diagonalized
        Hamiltonian), not an orbital-projected DOS.

    Attributes
    ----------
    dos : np.ndarray
        Density of states array with shape (nband, nene)
    ene : np.ndarray
        Energy grid points array with shape (nene,)
    ene_num : int
        Number of energy points
    norb : int
        Number of bands (rows of ``dos``)
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


def read_chemical_potential(file_name: str) -> float | None:
    """Parse the chemical potential from an ``energy.dat`` file.

    The file is written by the UHFk solver.  Two formats are recognised:

    * **Single mu-group** (most common): a line ``ChemicalPotential = <value>``
      is present.  The value is returned directly.
    * **Multiple mu-groups**: only ``ChemicalPotential_<g> = <value>`` lines
      are present (one per group).  The value for the *lowest* group number
      present is returned and a warning is logged because the DoS energy axis
      is ambiguous when multiple Fermi levels exist.
    * **Absent**: no ``ChemicalPotential`` line at all → ``None`` is returned.

    Parameters
    ----------
    file_name : str
        Path to the ``energy.dat`` file produced by UHFk.

    Returns
    -------
    float or None
        Chemical potential value, or ``None`` if not found or file is missing.
    """
    try:
        with open(file_name, "r") as f:
            lines = f.readlines()
    except OSError:
        return None

    # Compiled patterns: exact match first, then per-group match.
    pat_single = re.compile(r"^\s*ChemicalPotential\s*=\s*([-+0-9eE.]+)\s*$")
    pat_group = re.compile(r"^\s*ChemicalPotential_(\d+)\s*=\s*([-+0-9eE.]+)\s*$")

    group_vals: dict[int, float] = {}
    for line in lines:
        m = pat_single.match(line)
        if m:
            try:
                return float(m.group(1))
            except ValueError:
                pass
        m = pat_group.match(line)
        if m:
            try:
                group_vals[int(m.group(1))] = float(m.group(2))
            except ValueError:
                pass

    if group_vals:
        g0 = min(group_vals)
        logger.warning(
            "read_chemical_potential: only per-group ChemicalPotential entries "
            "found in '%s'; using the lowest group %d (mu=%.6g) for the DoS "
            "energy axis. The result is ambiguous for multi-group calculations.",
            file_name, g0, group_vals[g0],
        )
        return group_vals[g0]

    return None


def calc_dos(
    input_dict: dict,
    ene_window: list | None = None,
    ene_num: int = 101,
    verbose: bool = False,
    mu: float | None = None,
) -> DoS:
    """Calculate density of states.

    Notes
    -----
    This consumes the k-resolved eigenvalues written by the UHFk solver (a
    single ``<eigen>.npz`` with an ``eigenvalue`` array of shape
    ``(nk_sub, nband)``) and uses the tetrahedron method over the k-grid. It is
    UHFk-only: the per-block ``<block>_<eigen>.npz`` files written by UHFr are
    not consumed here.

    By default (``mu=None``) the energy axis uses the raw eigenvalue scale.
    When ``mu`` is given, eigenvalues are shifted to ``E - mu`` so that the
    Fermi level appears at zero (standard condensed-matter convention).  Use
    :func:`read_chemical_potential` to obtain ``mu`` from the ``energy.dat``
    file written by the UHFk solver.

    Parameters
    ----------
    input_dict : dict
        Input parameters dictionary
    ene_window : list, optional
        Energy window [min, max] for DoS calculation.  Applied *after* the
        ``mu`` shift when ``mu`` is given.
    ene_num : int, optional
        Number of energy points
    verbose : bool, optional
        If True, print additional output
    mu : float or None, optional
        Chemical potential to subtract from eigenvalues.  When ``None``
        (default), eigenvalues are used as-is (backward-compatible behavior).
        When a float is given, the energy axis is shifted to ``E - mu`` so
        that the Fermi level sits at zero.

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
    if mu is not None:
        if verbose:
            print("Subtracting chemical potential mu = ", mu, " from eigenvalues")
        eigenvalues = eigenvalues - mu
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
        "--subtract-mu",
        action="store_true",
        default=False,
        help=(
            "Subtract the chemical potential μ from eigenvalues so that the "
            "Fermi level appears at E=0.  μ is read from the 'energy.dat' file "
            "in the output directory (written by the UHFk solver).  "
            "If no ChemicalPotential line is found, a warning is emitted and the "
            "raw eigenvalue axis is used."
        ),
    )
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

    mu = None
    if args.subtract_mu:
        output_dir_for_mu = input_dict["file"]["output"]["path_to_output"]
        energy_file = os.path.join(
            output_dir_for_mu,
            input_dict["file"]["output"].get("energy", "energy.dat"),
        )
        mu = read_chemical_potential(energy_file)
        if mu is None:
            logger.warning(
                "--subtract-mu requested but no ChemicalPotential was found in "
                "'%s'; proceeding with the raw eigenvalue axis (Fermi level NOT "
                "at 0).", energy_file,
            )
        elif verbose:
            print("--subtract-mu: using mu = {} from '{}'".format(mu, energy_file))

    dos = calc_dos(
        input_dict,
        ene_window=args.ene_window,
        ene_num=args.ene_num,
        verbose=verbose,
        mu=mu,
    )
    output_dir = input_dict["file"]["output"]["path_to_output"]
    if args.output != "":
        dos.write_dos(os.path.join(output_dir, args.output), verbose=verbose)
    if args.plot:
        dos.plot(os.path.join(output_dir, args.plot), verbose=verbose)
