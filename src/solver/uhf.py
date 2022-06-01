import itertools
import logging
import numpy as np
import itertools
from requests.structures import CaseInsensitiveDict
import os

logger = logging.getLogger("qlms").getChild("uhf")

class Interact_UHF_base():
    def __init__(self, ham_info, Nsize):
        self.Nsize = Nsize
        self.Ham_tmp = np.zeros(tuple([(2 * self.Nsize) for i in range(4)]), dtype=complex)
        self.Ham_trans_tmp = np.zeros(tuple([(2 * self.Nsize) for i in range(2)]), dtype=complex)
        self.param_ham = self._transform_interall(ham_info)

    #Change interaction to interall type
    def _transform_interall(self, ham_info):
        return ham_info

    def _calc_hartree(self):
        site = np.zeros(4, dtype=np.int32)
        for site_info, value in self.param_ham.items():
            for i in range(4):
                site[i] = site_info[2 * i] + site_info[2 * i + 1] * self.Nsize
            # Diagonal Fock term
            self.Ham_tmp[site[0]][site[1]][site[2]][site[3]] += value
            self.Ham_tmp[site[2]][site[3]][site[0]][site[1]] += value
            if site[1] == site[2]:
                self.Ham_trans_tmp[site[1]][site[2]] += value
        pass

    def _calc_fock(self):
        site = np.zeros(4,dtype=np.int32)
        for site_info, value in self.param_ham.items():
            for i in range(4):
                site[i] = site_info[2 * i] + site_info[2 * i + 1] * self.Nsize
            # OffDiagonal Fock term
            self.Ham_tmp[site[0]][site[3]][site[2]][site[1]] -= value
            self.Ham_tmp[site[2]][site[1]][site[0]][site[3]] -= value
        pass

    def get_ham(self, type):
        self._calc_hartree()
        if type == "hartreefock":
            self._calc_fock()
        return self.Ham_tmp, self.Ham_trans_tmp

class CoulombIntra_UHF(Interact_UHF_base):
    def _transform_interall(self, ham_info):
        param_tmp = {}
        for site_info, value in ham_info.items():
            site_info = tuple([site_info[0], 0, site_info[0], 0, site_info[0], 1, site_info[0], 1])
            param_tmp[site_info] = value
        return param_tmp

class CoulombInter_UHF(Interact_UHF_base):
    def _transform_interall(self, ham_info):
        param_tmp = {}
        for site_info, value in ham_info.items():
            for spin_i, spin_j in itertools.product([0,1], repeat=2):
                site_info = tuple([site_info[0], spin_i, site_info[0], spin_i, site_info[1], spin_j, site_info[1], spin_j])
                param_tmp[site_info] = value
        return param_tmp

class Hund_UHF(Interact_UHF_base):
    def _transform_interall(self, ham_info):
        param_tmp = {}
        for site_info, value in ham_info.items():
            for spin_i in range(2):
                site_info = tuple([site_info[0], spin_i, site_info[0], spin_i, site_info[1], spin_i, site_info[1], spin_i])
                param_tmp[site_info] = -value
        return param_tmp

class PairHop_UHF(Interact_UHF_base):
    def _transform_interall(self, ham_info):
        param_tmp = {}
        for site_info, value in ham_info.items():
            site_info = tuple([site_info[0], 0, site_info[1], 0, site_info[0], 1, site_info[1], 1])
            param_tmp[site_info] += value
            site_info = tuple([site_info[1], 1, site_info[0], 1, site_info[1], 0, site_info[0], 0])
            param_tmp[site_info] = value
        return param_tmp

class Exchange_UHF(Interact_UHF_base):
    def _transform_interall(self, ham_info):
        param_tmp = {}
        for site_info, value in ham_info.items():
            site_info = tuple([site_info[0], 0, site_info[1], 0, site_info[1], 1, site_info[0], 1])
            param_tmp[site_info] = value
            site_info = tuple([site_info[0], 1, site_info[1], 1, site_info[1], 0, site_info[0], 0])
            param_tmp[site_info] = value
        return param_tmp

class Ising_UHF(Interact_UHF_base):
    def _transform_interall(self, ham_info):
        param_tmp = {}
        for site_info, value in ham_info.items():
            for spin_i, spin_j in itertools.product([0,1], repeat=2):
                site_info = tuple([site_info[0], spin_i, site_info[0], spin_i, site_info[1], spin_j, site_info[1], spin_j])
                if spin_i != spin_j:
                    value *= -1.0
                param_tmp[site_info] = value * (1-2*spin_i) * (1-2*spin_j)
        return param_tmp

class PairLift_UHF(Interact_UHF_base):
    def _transform_interall(self, ham_info):
        param_tmp = {}
        for site_info, value in ham_info.items():
            site_info = tuple([site_info[0], 0, site_info[0], 1, site_info[1], 0, site_info[1], 1])
            param_tmp[site_info] = value
            site_info = tuple([site_info[1], 1, site_info[1], 0, site_info[0], 1, site_info[0], 0])
            param_tmp[site_info] = value
        return param_tmp


# It is better to define solver class, since UHF-k and RPA mode will be added.

class UHF(object):
    def __init__(self, param, param_ham, info_log):
        self.param = param
        self.param_ham = param_ham
        self.info_log = info_log
        self.physics = {"Ene": 0, "NCond": 0, "Sz": 0, "Rest": 1.0}
        self.output_list = ["energy", "eigen", "green"]
        # initial values
        para_init = CaseInsensitiveDict({
            "Nsite": 0,
            "Ne": 0,
            "Ncond": 0,
            "2Sz": None,
            "Mix": 0.5,
            "EPS": 6,
            "Print": 0,
            "IterationMax": 20000,
            "RndSeed": 1234,
            "EpsSlater": 6
        })

        for k, v in para_init.items():
            self.param["ModPara"].setdefault(k, v)

        # set initial values
        for key in ["CDataFileHead", "CParaFileHead"]:
            if self.param["ModPara"][key] is not None and type(self.param["ModPara"][key]) == type([]):
                self.param["ModPara"][key] = str(self.param["ModPara"][key][0])

        for key in ["nsite", "ne", "2Sz", "ncond", "eps", "IterationMax", "Print", "RndSeed", "EpsSlater"]:
            if self.param["ModPara"][key] is not None and type(self.param["ModPara"][key]) == type([]):
                self.param["ModPara"][key] = int(self.param["ModPara"][key][0])

        #The type of mix is float.
        self.param["ModPara"]["mix"] = float(self.param["ModPara"]["mix"][0])

        self.param["ModPara"]["EPS"] = pow(10, -self.param["ModPara"]["EPS"])

        output_str = "Show input parameters\n"
        for k, v in self.param["ModPara"].items():
            output_str += "{}: {}\n".format(k, v)
        logger.debug(output_str)

        self.Nsize = self.param["ModPara"]["nsite"]
        # Define list
        param_mod = self.param["Modpara"]
        self.Ncond = param_mod["Ncond"]
        # TwoSz = 0 if param_mod["2Sz"] is None else param_mod["2Sz"]
        self.green_list = {"free": {"label": [i for i in range(2 * self.Nsize)], "occupied": self.Ncond}}

        #ToDo: Add flag for considering Fock term
        self.iflag_fock = True


    def solve(self, path_to_output):
        print_level = self.info_log["print_level"]
        print_step = self.info_log["print_step"]
        logger.info("Set Initial Green's functions")
        # Get label for eigenvalues
        label = 0
        for k, v in self.green_list.items():
            self.green_list[k]["eigen_start"] = label
            label += len(self.green_list[k]["label"])
        self.Green = self._initial_G()
        logger.info("Start UHF calculations")
        param_mod = self.param["Modpara"]

        if print_level > 0:
            logger.info("step, rest, energy, NCond, Sz")
        self._makeham_const()
        self._makeham_mat()
        for i_step in range(param_mod["IterationMax"]):
            self._makeham()
            self._diag()
            self._determine_occupied_num()
            self._green()
            self._calc_energy()
            self._calc_phys()
            if i_step % print_step == 0 and print_level > 0:
                logger.info(
                    "{}, {:.8g}, {:.8g}, {:.4g}, {:.4g} ".format(i_step, self.physics["Rest"], self.physics["Ene"]["Total"],
                                                                 self.physics["NCond"], self.physics["Sz"]))
            if self.physics["Rest"] < param_mod["eps"]:
                break

        if self.physics["Rest"] < param_mod["eps"]:
            logger.info("UHF calculation is succeeded: rest={}, eps={}.".format(self.physics["Rest"], param_mod["eps"]))
        else:
            logger.warning("UHF calculation is failed: rest={}, eps={}.".format(self.physics["Rest"], param_mod["eps"]))

    def _initial_G(self):
        _green_list = self.green_list
        green = np.zeros((2 * self.Nsize, 2 * self.Nsize), dtype=complex)
        if self.param_ham["Initial"] is not None:
            g_info = self.param_ham["Initial"]
            for site_info, value in g_info.items():
                site1 = site_info[0] + site_info[1] * self.Nsize
                site2 = site_info[2] + site_info[3] * self.Nsize
                green[site1][site2] = value
        else:
            np.random.seed(self.param["ModPara"]["RndSeed"])
            rand = np.random.rand(2 * self.Nsize * 2 * self.Nsize).reshape(2 * self.Nsize, 2 * self.Nsize)
            for k, info in _green_list.items():
                v = info["label"]
                for site1, site2 in itertools.product(v, v):
                    green[site1][site2] = 0.01 * (rand[site1][site2] - 0.5)
        return green

    def _makeham_const(self):
        self.Ham_trans = np.zeros((2 * self.Nsize, 2 * self.Nsize), dtype=complex)
        # Transfer integrals
        for site_info, value in self.param_ham["Transfer"].items():
            site1 = site_info[0] + site_info[1] * self.Nsize
            site2 = site_info[2] + site_info[3] * self.Nsize
            self.Ham_trans[site1][site2] += -value

    def _makeham(self):
        self.Ham = np.zeros((2 * self.Nsize, 2 * self.Nsize), dtype=complex)
        self.Ham = self.Ham_trans.copy()
        green_local = self.Green.reshape((2 * self.Nsize) ** 2)
        self.Ham += np.dot(self.Ham_local, green_local).reshape((2 * self.Nsize), (2 * self.Nsize))

    def _makeham_mat(self):
        # TODO Add Hund, Exchange, Ising, PairHop, and PairLift
        self.Ham_local = np.zeros(tuple([(2 * self.Nsize) for i in range(4)]), dtype=complex)
        if self.iflag_fock is True:
            type = "hartreefock"
        else:
            type = "hartree"

        for key in ["CoulombIntra", "CoulombInter", "Hund", "Exchange", "Ising", "PairHop", "PairLift", "InterAll"]:
            if self.param_ham[key] is not None:
                param_ham = self.param_ham[key]
                if key == "CoulombIntra":
                    ham_uhf = CoulombIntra_UHF(param_ham, self.Nsize)
                elif key == "CoulombInter":
                    ham_uhf = CoulombInter_UHF(param_ham, self.Nsize)
                elif key == "Hund":
                    ham_uhf = Hund_UHF(param_ham, self.Nsize)
                elif key == "Exchange":
                    ham_uhf = Exchange_UHF(param_ham, self.Nsize)
                elif key == "Ising":
                    ham_uhf = Ising_UHF(param_ham, self.Nsize)
                elif key == "PairHop":
                    ham_uhf = PairHop_UHF(param_ham, self.Nsize)
                elif key == "PairLift":
                    ham_uhf = PairLift_UHF(param_ham, self.Nsize)
                elif key == "InterAll":
                    ham_uhf = Interact_UHF_base(param_ham, self.Nsize)
                else:
                    logger.warning("key {} is wrong!".format(key))
                    exit(1)

                Ham_local_tmp, Ham_trans_tmp = ham_uhf.get_ham(type=type)
                self.Ham_local += Ham_local_tmp
                self.Ham_trans += Ham_trans_tmp

        self.Ham_local = self.Ham_local.reshape((2 * self.Nsize) ** 2, (2 * self.Nsize) ** 2)

    def _diag(self):
        _green_list = self.green_list
        for k, block_g_info in _green_list.items():
            g_label = block_g_info["label"]
            block_size = len(g_label)
            mat = np.zeros((block_size, block_size), dtype=complex)
            for site1, org_site1 in enumerate(g_label):
                for site2, org_site2 in enumerate(g_label):
                    mat[site1][site2] = self.Ham[org_site1][org_site2]
            # w: The eigenvalues in ascending order, each repeated according to its multiplicity.
            # v: The column v[:, i] is the normalized eigenvector corresponding to the eigenvalue w[i].
            # Will return a matrix object if a is a matrix object.
            w, v = np.linalg.eigh(mat)
            self.green_list[k]["eigenvalue"] = w
            self.green_list[k]["eigenvector"] = v

    def _determine_occupied_num(self):
        pass

    def _green(self):
        _green_list = self.green_list
        # R_SLT = U^{*} in _green
        # L_SLT = U^T in _green
        R_SLT = np.zeros((2 * self.Nsize, 2 * self.Nsize), dtype=complex)
        L_SLT = np.zeros((2 * self.Nsize, 2 * self.Nsize), dtype=complex)
        for k, block_g_info in _green_list.items():
            g_label = block_g_info["label"]
            # occupied: return the occupied number
            occupied_number = block_g_info["occupied"]
            eigenvec = self.green_list[k]["eigenvector"]
            eigen_start = self.green_list[k]["eigen_start"]
            for eigen_i in range(occupied_number):
                evec_i = eigenvec[:, eigen_i]
                for idx, org_site in enumerate(g_label):
                    R_SLT[org_site][eigen_start + eigen_i] = np.conjugate(evec_i[idx])
                    L_SLT[eigen_start + eigen_i][org_site] = evec_i[idx]
        # Copy Green_old
        self.Green_old = self.Green.copy()
        RMat = np.dot(R_SLT, L_SLT)
        self.Green = RMat.copy()

    def _calc_energy(self):
        _green_list = self.green_list
        Ene = {}
        Ene["band"] = 0
        for k, block_g_info in _green_list.items():
            eigenvalue = self.green_list[k]["eigenvalue"]
            occupied_number = block_g_info["occupied"]
            Ene["band"] += np.sum(eigenvalue[:occupied_number])

        Ene["InterAll"] = 0
        green_local = self.Green.reshape((2 * self.Nsize) ** 2)
        Ene["InterAll"] -= np.dot(green_local.T, np.dot(self.Ham_local, green_local))/2.0

        ene = 0
        for value in Ene.values():
            ene += value
        Ene["Total"] = ene
        self.physics["Ene"] = Ene
        logger.debug(Ene)

    def _calc_phys(self):
        n = 0
        for site in range(2 * self.Nsize):
            n += self.Green[site][site]
        self.physics["NCond"] = n.real

        sz = 0
        for site in range(self.Nsize):
            sz += self.Green[site][site]
            sz -= self.Green[site + self.Nsize][site + self.Nsize]
        self.physics["Sz"] = (0.5 * sz).real

        rest = 0.0
        mix = self.param["ModPara"]["mix"]
        for site1, site2 in itertools.product(range(2 * self.Nsize), range(2 * self.Nsize)):
            rest += abs(self.Green[site1][site2] - self.Green_old[site1][site2])**2
            self.Green[site1][site2] = self.Green_old[site1][site2] * (1.0 - mix) + mix * self.Green[site1][site2]
        self.physics["Rest"] = np.sqrt(rest) / (2.0 * self.Nsize * self.Nsize)
        self.Green[np.where(abs(self.Green) < 1e-12)] = 0

    def get_results(self):
        return (self.physics, self.Green)

    def save_results(self, info_outputfile, green_info):
        path_to_output = info_outputfile["path_to_output"]

        if "energy" in info_outputfile.keys():
            output_str  = "Energy_total = {}\n".format(self.physics["Ene"]["Total"].real)
            output_str += "Energy_band = {}\n".format(self.physics["Ene"]["band"])
            output_str += "Energy_interall = {}\n".format(self.physics["Ene"]["InterAll"].real)
            output_str += "NCond = {}\n".format(self.physics["NCond"])
            output_str += "Sz = {}\n".format(self.physics["Sz"])
            with open(os.path.join(path_to_output, info_outputfile["energy"]), "w") as fw:
                fw.write(output_str)

        if "eigen" in info_outputfile.keys():
            output_str = ""
            for key, _green_list in self.green_list.items():
                eigenvalue = _green_list["eigenvalue"]
                #output_str += "{}, {}\n".format(key, eigenvalue.shape[0])
                for idx, value in enumerate(eigenvalue):
                    output_str += "{}, {}\n".format(idx, value)
            with open(os.path.join(path_to_output, info_outputfile["eigen"]), "w") as fw:
                fw.write(output_str)

        if "green" in info_outputfile.keys():
            with open(os.path.join(path_to_output, info_outputfile["green"]), "w") as fw:
                pass
            if "OneBodyG" in green_info:
                _green_info = green_info["OneBodyG"]
                _green_info = np.array(_green_info, dtype=np.int32)
                output_str = ""
                for info in _green_info:
                    nsite1 = info[0] + self.Nsize*info[1]
                    nsite2 = info[2] + self.Nsize*info[3]
                    output_str += "{} {} {} {} {} {}\n".format(info[0], info[1], info[2], info[3], self.Green[nsite1][nsite2].real, self.Green[nsite1][nsite2].imag)
                with open(os.path.join(path_to_output, info_outputfile["green"]), "w") as fw:
                    fw.write(output_str)

    def get_Ham(self):
        return self.Ham


class param():
    def __init__(self):
        self.para = CaseInsensitiveDict({
            "Nsite": 0,
            "Ne": 0,
            "Ncond": 0,
            "2Sz": None,
            "Mix": 0.5,
            "EPS": 8,
            "Print": 0,
            "IterationMax": 20000,
            "RndSeed": 0,
            "EpsSlater": 6,
            "NMPTrans": 1,
        })

    def set_param(self, param):
        for k, v in param.items():
            self.para[k] = v
        return self.para

    def print_param(self):
        print("########Input parameters ###########")
        for key in ["tmp_eps", "eps_int", "mix", "print"]:
            print("{} {}".format(key, self.para[key]))
        print("####################################")

    def output_phys(self):
        pass


class UHF_SzFix(UHF):
    def __init__(self, param, param_ham, totalsz):
        super(UHF_SzFix, self).__init__(param, param_ham)
        Ncond = self.Ncond
        TwoSz = totalsz * 2
        self.green_list = {
            "up": {"label": [i for i in range(self.Nsize)], "value": 0.5, "occupied": int((self.Ncond + TwoSz) / 2)},
            "down": {"label": [i for i in range(self.Nsize, 2 * self.Nsize)], "value": -0.5,
                     "occupied": int((self.Ncond - TwoSz) / 2)}}
