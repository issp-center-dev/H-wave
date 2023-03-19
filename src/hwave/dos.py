from __future__ import annotations

import itertools
import os

import numpy as np


class DoS:
    dos: np.ndarray
    ene: np.ndarray
    ene_num: int
    norb: int

    def __init__(self, ene: np.ndarray, dos: np.ndarray):
        assert ene.shape[0] == dos.shape[1]
        self.ene = ene
        self.dos = dos
        self.ene_num = ene.shape[0]
        self.norb = dos.shape[0]

    def plot(self, filename: str = ""):
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
            plt.savefig(filename)
        plt.close()

    def write_dos(self, output: str):
        with open(output, "w") as fw:
            fw.write("# ene")
            for j in range(self.norb):
                fw.write(f" dos[iorb={j}]")
            fw.write("\n")
            for i in range(self.ene_num):
                fw.write("{:15.8f} ".format(self.ene[i]))
                for j in range(self.norb):
                    fw.write("{:15.8f} ".format(self.dos[j, i]))
                fw.write("\n")


def __read_geom(file_name="./dir-model/zvo_geom.dat"):
    with open(file_name, "r") as fr:
        uvec = np.zeros((3, 3))
        for i, line in enumerate(itertools.islice(fr, 3)):  # take first 3 lines
            uvec[i, :] = np.array(line.split())
    return uvec


def calc_dos(
    input_dict: dict,
    ene_window: list | None = None,
    ene_num: int = 100,
    verbose: bool = False,
) -> DoS:

    try:
        import libtetrabz
    except ImportError:
        raise ImportError(
            "libtetrabz is not installed. Please install libtetrabz and try again."
        )

    if verbose:
        print("Reading eigenvalues")
    output_info_dict = input_dict["file"]["output"]
    data = np.load(
        os.path.join(
            output_info_dict["path_to_output"], output_info_dict["eigen"] + ".npz"
        )
    )
    eigenvalues = data["eigenvalue"]
    Lx, Ly, Lz = input_dict["mode"]["param"]["CellShape"]
    norb = eigenvalues.shape[1]
    if verbose:
        print("Lx, Ly, Lz, norb: ", Lx, Ly, Lz, norb)
    eigenvalues.reshape(Lx, Ly, Lz, norb)

    if verbose:
        print("Reading geometry")
    input_info_dict = input_dict["file"]["input"]["interaction"]
    file_name = os.path.join(
        input_info_dict["path_to_input"], input_info_dict["Geometry"]
    )
    uvec = __read_geom(file_name)
    bvec = 2.0 * np.pi * np.linalg.inv(uvec).T

    if ene_window is None:
        ene_min = np.min(eigenvalues) - 0.2
        ene_max = np.max(eigenvalues) + 0.2
    else:
        ene_min = ene_window[0]
        ene_max = ene_window[1]

    ene = np.linspace(ene_min, ene_max, num=ene_num)
    if verbose:
        print("ene_min, ene_max, ene_num: ", ene_min, ene_max, ene_num)

    eig = eigenvalues.reshape(Lx, Ly, Lz, norb)
    if verbose:
        print("Calculating DOS")
    wght = libtetrabz.dos(bvec, eig, ene)
    dos = wght.sum(2).sum(1).sum(0)
    return DoS(dos=dos, ene=ene)


def main():
    import tomli
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input", type=str, help="input file of hwave"
    )
    parser.add_argument("--output_dos", type=str, default="dos.dat", help="DoS output")
    parser.add_argument(
        "--ene_window",
        default=None,
        type=float,
        help="energy window; [ene_low, ene_high]. If None, ene_low = ene_min - 0.2, ene_high = ene_max + 0.2",
    )
    parser.add_argument("--ene_num", default=100, type=int, help="energy step")
    parser.add_argument("--plot", action="store_true", help="plot DOS")
    parser.add_argument("--verbose", action="store_true", help="print verbosely")

    args = parser.parse_args()

    file_toml = args.input
    if os.path.exists(file_toml):
        if args.verbose:
            print("Reading input file: ", file_toml)
        with open(file_toml, "rb") as f:
            input_dict = tomli.load(f)
    else:
        raise ValueError("Input file does not exist")

    dos = calc_dos(
        input_dict,
        ene_window=args.ene_window,
        ene_num=args.ene_num,
        verbose=args.verbose,
    )
    if args.output_dos != "":
        dos.write_dos(args.output_dos)
    if args.plot:
        dos.plot("dos.png")
