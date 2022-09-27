import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), "../src/qlmsio"))
from wan90 import read_geometry
from wan90 import write_geom
from wan90 import read_w90
from wan90 import write_w90
from read_input import QMLSInput

interact_shape = [3, 3, 1]
#Read information of geom file
info_geometry = read_geometry("geometry.dat")
write_geom("zvo_geom.dat", info_geometry)
read_io = QMLSInput("namelist.def")
ham_info = read_io.get_param("ham")

def transfer_org(ham_transfer, n_orb):
    ham_org = {}
    for site, value in ham_transfer.items():
        if site[0] < n_orb:
            ham_org[site] = value
    return ham_org

n_orb = info_geometry["n_orb"]
info_int = {}
info_int["n_orb"] = n_orb

if "Transfer" in ham_info.keys():
    ham = transfer_org(ham_info["Transfer"], n_orb)
    info_int["interaction"] = {}
    for site, value in ham.items():
        # In this version, spin must be same.
        if site[1] != site[3]:
            print("Error: spin components are different.")
            exit(1)
        elif site[1] == 0:
            info_int["interaction"][site[2]] = [site[0], value]
    write_w90("zvo_hr.def", info_int, info_geometry, interact_shape)

def coulombintra_org(ham_coulomb_intra, n_orb):
    ham_org = {}
    for site, value in ham_coulomb_intra.items():
        if site[0] < n_orb:
            ham_org[site] = value
    return ham_org

info_int["interaction"] = {}
if ham_info["Coulombintra"] is not None:
    # i value
    ham = transfer_org(ham_info["CoulombIntra"], n_orb)
    for site, value in ham.items():
        info_int["interaction"][site[0]] = [site[0], value]

for key in ["Coulombinter"]:
    if ham_info[key] is not None:
        # i j value
        ham = transfer_org(ham_info[key], n_orb)
        for site, value in ham.items():
            info_int["interaction"][site[1]] = [site[0], value]
write_w90("zvo_ur.def", info_int, info_geometry, interact_shape)

for key in ["Hund", "PairHop", "Exchange", "Ising", "PairLift"]:
    info_int["interaction"] = {}
    if ham_info[key] is not None:
        # i j value
        ham = transfer_org(ham_info[key], n_orb)
        for site, value in ham.items():
            info_int["interaction"][site[1]] = [site[0], value]
    write_w90("zvo_{}.def".format(key), info_int, info_geometry, interact_shape)

