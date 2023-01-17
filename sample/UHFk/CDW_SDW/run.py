import os
import os.path
import sys
import subprocess

import tomli
import numpy as np
import scipy.optimize as opt

import hwave.qlms

if len(sys.argv) == 1:
    print(f"Usage: python3 {sys.argv[0]} [path_to_uhf_dry]")
    sys.exit(1)
uhf_dry_path = sys.argv[1]

t = 1.0
U = 4.0
Vs = np.linspace(0.0, 2.0, num=20)

# parameters for finding phase boundary

Vmin_b = 0.8
Vmax_b = 1.2

baseinput = "input.toml"
with open(baseinput, "rb") as f:
    param_dict = tomli.load(f)


def write_initial(filename, phase="cdw"):
    G = np.zeros((64, 2, 4, 2, 4), dtype=np.complex128)
    if phase == "cdw":
        G[0, 0, 0, 0, 0] = 1.0
        G[0, 0, 3, 0, 3] = 1.0
        G[0, 1, 0, 1, 0] = 1.0
        G[0, 1, 3, 1, 3] = 1.0
    elif phase == "sdw":
        G[0, 0, 0, 0, 0] = 1.0
        G[0, 1, 1, 1, 1] = 1.0
        G[0, 1, 2, 1, 2] = 1.0
        G[0, 0, 3, 0, 3] = 1.0
    np.savez(filename, green_sublattice=G)


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


def run_hwave(output_dir: str, U: float, V: float, initial: str = ""):
    stan_in = "stan.in"
    write_stanin(stan_in, t=t, U=U, V=V)
    cmd = [uhf_dry_path, stan_in]
    subprocess.run(cmd)

    file_input = param_dict["file"]["input"]
    file_input_interaction = {
        "Geometry": "geom.dat",
        "Transfer": "transfer.dat",
    }
    if U != 0.0:
        file_input_interaction["CoulombIntra"] = "coulombintra.dat"
    if V != 0.0:
        file_input_interaction["CoulombInter"] = "coulombinter.dat"
    file_input["interaction"] = file_input_interaction
    if initial != "":
        write_initial("initial.npz", phase=initial)
        file_input["initial"] = "initial.npz"
    param_dict["file"]["input"] = file_input
    param_dict["file"]["output"]["path_to_output"] = output_dir
    hwave.qlms.run(input_dict=param_dict)


def energy_diff(V: float) -> float:
    output_dir = "output_find_energy"
    energy = {}
    for initial in ("cdw", "sdw"):
        run_hwave(output_dir, U=U, V=V, initial=initial)
        with open(os.path.join(output_dir, "energy.dat")) as fe:
            for line in fe:
                words = line.strip().split("=")
                if words[0].strip() == "Energy_Total":
                    energy[initial] = float(words[1].strip())
                    break
    diff = energy["cdw"] - energy["sdw"]
    with open("find_crossing.log", "a") as f:
        f.write(f"{V} {diff} {energy['cdw']} {energy['sdw']}\n")
    return diff


# search for transition point
if os.path.exists("find_crossing.log"):
    os.remove("find_crossing.log")
sol = opt.root_scalar(energy_diff, bracket=[Vmin_b, Vmax_b])
phase_boundary = sol.root

# V sweep
f = open("res.dat", "w")
f.write(f"# Vc = {phase_boundary}\n")
f.write("# $1:  V\n")
f.write("# $2:  E        [GS]\n")
f.write("# $3:  n(pi,pi) [GS]\n")
f.write("# $4:  s(pi,pi) [GS]\n")
f.write("# $5:  E        [from CDW]\n")
f.write("# $6:  n(pi,pi) [from CDW]\n")
f.write("# $7:  s(pi,pi) [from CDW]\n")
f.write("# $8:  E        [from SDW]\n")
f.write("# $9:  n(pi,pi) [from SDW]\n")
f.write("# $10: s(pi,pi) [from SDW]\n")
f.write("# $11: name of stable phase\n")

for i, V in enumerate(Vs):
    res = {"cdw": {}, "sdw": {}}
    for initial in ("cdw", "sdw"):
        output_dir = f"output_{i}_{initial}"
        run_hwave(output_dir, U=U, V=V, initial=initial)

        with open(os.path.join(output_dir, "energy.dat")) as fe:
            for line in fe:
                words = line.strip().split("=")
                if words[0].strip() == "Energy_Total":
                    res[initial]["ene"] = float(words[1].strip())
                    break

        green_file = os.path.join(output_dir, "green.dat.npz")
        g = np.load(green_file)["green_sublattice"]
        density = np.real(g[0, 0, :, 0, :]) + np.real(g[0, 1, :, 1, :])
        spin = np.real(g[0, 0, :, 0, :]) - np.real(g[0, 1, :, 1, :])
        res[initial]["cdw"] = 0.25 * (
            density[0, 0] - density[1, 1] - density[2, 2] + density[3, 3]
        )
        res[initial]["sdw"] = 0.25 * (spin[0, 0] - spin[1, 1] - spin[2, 2] + spin[3, 3])
    f.write(f"{V}")
    if V < phase_boundary:
        stable_phase = "sdw"
    else:
        stable_phase = "cdw"
    for initial in (stable_phase, "cdw", "sdw"):
        r = res[initial]
        f.write(f" {r['ene']} {r['cdw']} {np.abs(r['sdw'])}")
    f.write(f" {stable_phase}")
    f.write(f"\n")
f.close()
print(f"Vc = {phase_boundary}")
