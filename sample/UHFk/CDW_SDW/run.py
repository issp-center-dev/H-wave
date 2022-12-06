import os
import sys
import subprocess
import os.path

import tomli
import numpy as np

if len(sys.argv) == 1:
    print(f"Usage: python3 {sys.argv[0]} [path_to_uhf_dry]")
    sys.exit(1)
uhf_dry_path = sys.argv[1]

import hwave.qlms

t = 1.0
U = 4.0
Vs = np.linspace(0.0, 2.0, num=9)

baseinput = "input.toml"
with open(baseinput, "rb") as f:
    param_dict = tomli.load(f)


def write_stanin(filepath: str, *, t: float = 1.0, U: float = 0.0, V: float = 0.0):
    std_template = """\
model = "Hubbard"
lattice = "square"
W = 4
L = 4
t = {t}
U = {U}
V = {V}
Ncond = 16
eps = 8
calcmode = "uhfk"
exportall = 0
"""
    stanin = std_template.format(t=t, U=U, V=V)
    with open(filepath, "w") as f:
        f.write(stanin)


f = open("res.dat", "w")
for i, V in enumerate(Vs):
    stan_in = "stan.in"
    write_stanin(stan_in, t=t, U=U, V=V)
    cmd = [uhf_dry_path, stan_in]
    subprocess.run(cmd)

    output_dir = f"output_{i}"
    mode_param = param_dict["mode"]["param"]

    file_input = param_dict["file"]["input"]
    file_input_interaction = {
        "Geometry": "geom.dat",
        "Transfer": "transfer.dat",
    }
    if U != 0.0:
        file_input_interaction["CoulombIntra"] = "coulombintra.dat"
    if V != 0.0:
        file_input_interaction["CoulombInter"] = "coulombinter.dat"
    # if i > 0:
    #     file_input["initial"] = f"output_{i-1}/green.dat.npz"
    # elif "initial" in file_input:
    #     file_input.pop("initial")
    file_input["interaction"] = file_input_interaction
    param_dict["file"]["input"] = file_input

    param_dict["file"]["output"]["path_to_output"] = output_dir

    hwave.qlms.run(input_dict=param_dict)

    green_file = os.path.join(output_dir, "green.dat.npz")
    g = np.load(green_file)["green_sublattice"]
    density = np.real(g[0, 0, :, 0, :]) + np.real(g[0, 1, :, 1, :])
    spin = np.real(g[0, 0, :, 0, :]) - np.real(g[0, 1, :, 1, :])
    CDW = 0.25 * (density[0, 0] - density[1, 1] - density[2, 2] + density[3, 3])
    SDW = 0.25 * (spin[0, 0] - spin[1, 1] - spin[2, 2] + spin[3, 3])
    f.write(f"{V} {CDW} {SDW}\n")
f.close()
