import itertools
import logging
import numpy as np
from requests.structures import CaseInsensitiveDict

logger = logging.getLogger("qlms").getChild("uhf")


class UHF(object):
    def __init__(self, param, param_ham, output):
        self.param = param
        self.param_ham = param_ham
        self.output = output
        self.physics = {"Ene": 0, "NCond": 0, "Sz": 0, "Rest": 1.0}
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
            "EpsSlater": 6,
            "NMPTrans": 1,
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
        self.param["ModPara"]["mix"] = float(self.param["ModPara"]["mix"])

        self.param["ModPara"]["EPS"] = pow(10, -self.param["ModPara"]["EPS"])

        output_str = "Show input parameters\n"
        for k, v in self.param["ModPara"].items():
            output_str += "{}: {}\n".format(k, v)
        logger.debug(output_str)

        output_str = "Show numerical parameters\n"
        for key in ["eps", "mix", "print"]:
            output_str += "{} {}\n".format(key, self.param["ModPara"][key])
        logging.info(output_str)
        self.Nsize = len(param["Locspin"])

        # Define list
        param_mod = self.param["Modpara"]
        self.Ncond = param_mod["Ncond"]
        # TwoSz = 0 if param_mod["2Sz"] is None else param_mod["2Sz"]
        self.green_list = {"free": {"label": [i for i in range(2 * self.Nsize)], "occupied": self.Ncond}}

    def solve(self):
        logger.info("Set Initial Green's functions")
        # Get label for eigenvalues
        label = 0
        for k, v in self.green_list.items():
            self.green_list[k]["eigen_start"] = label
            label += len(self.green_list[k]["label"])
        self.Green = self._initial_G()
        logger.info("Start UHF calculations")
        param_mod = self.param["Modpara"]
        logger.info("step, rest, energy, NCond, Sz")
        self._makeham_const()
        self._makeham_mat()
        for i_step in range(param_mod["IterationMax"]):
            self._makeham()
            self._diag()
            self._detemine_occupied_num()
            self._green()
            self._calc_energy()
            self._calc_phys()
            logger.info(
                "{}, {:.8g}, {:.8g}, {:.4g}, {:.4g} ".format(i_step, self.physics["Rest"], self.physics["Ene"]["Total"],
                                                             self.physics["NCond"], self.physics["Sz"]))
            if self.physics["Rest"] < param_mod["eps"]:
                break
        pass

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
        if self.param_ham["InterAll"] is not None:
            site = np.zeros(4, int)
            for site_info, value in self.param_ham["InterAll"].items():
                for i in range(4):
                    site[i] = site_info[2 * i] + site_info[2 * i + 1] * self.Nsize
                if site[1] == site[2]:
                    self.Ham_trans[site[0]][site[3]] += value

    def _makeham(self):
        self.Ham = np.zeros((2 * self.Nsize, 2 * self.Nsize), dtype=complex)
        self.Ham = self.Ham_trans.copy()
        green_local = self.Green.reshape((2 * self.Nsize) ** 2)
        self.Ham += np.dot(self.Ham_local, green_local).reshape((2 * self.Nsize), (2 * self.Nsize))

    def _makeham_mat(self):
        # TODO Add Hund, Exchange and PairHop
        self.Ham_local = np.zeros(tuple([(2 * self.Nsize) for i in range(4)]), dtype=complex)
        # CoulombIntra
        if self.param_ham["CoulombIntra"] is not None:
            for site_info, value in self.param_ham["CoulombIntra"].items():
                u_site = site_info[0] + 0 * self.Nsize
                d_site = site_info[0] + 1 * self.Nsize
                self.Ham_local[u_site][u_site][d_site][d_site] += value
                self.Ham_local[d_site][d_site][u_site][u_site] += value
                self.Ham_local[u_site][d_site][d_site][u_site] -= value
                self.Ham_local[d_site][u_site][u_site][d_site] -= value

        # CoulombInter
        if self.param_ham["CoulombInter"] is not None:
            for site_info, value in self.param_ham["CoulombInter"].items():
                u_site = np.zeros(2, dtype=int)
                d_site = np.zeros(2, dtype=int)
                for i in range(2):
                    u_site[i] = site_info[i] + 0 * self.Nsize
                    d_site[i] = site_info[i] + 1 * self.Nsize

                # Hartree term
                for i in range(2):
                    self.Ham_local[u_site[i]][u_site[i]][u_site[1 - i]][u_site[1 - i]] += value
                    self.Ham_local[u_site[i]][u_site[i]][d_site[1 - i]][d_site[1 - i]] += value
                    self.Ham_local[d_site[i]][d_site[i]][u_site[1 - i]][u_site[1 - i]] += value
                    self.Ham_local[d_site[i]][d_site[i]][d_site[1 - i]][d_site[1 - i]] += value

                # Fock term
                # diagonal
                self.Ham_local[u_site[0]][u_site[1]][u_site[1]][u_site[0]] -= value
                self.Ham_local[u_site[1]][u_site[0]][u_site[0]][u_site[1]] -= value
                self.Ham_local[d_site[0]][d_site[1]][d_site[1]][d_site[0]] -= value
                self.Ham_local[d_site[1]][d_site[0]][d_site[0]][d_site[1]] -= value
                # off-diagonal
                self.Ham_local[u_site[0]][d_site[1]][d_site[1]][u_site[0]] -= value
                self.Ham_local[d_site[1]][u_site[0]][u_site[0]][d_site[1]] -= value
                self.Ham_local[u_site[1]][d_site[0]][d_site[0]][u_site[1]] -= value
                self.Ham_local[d_site[0]][u_site[1]][u_site[1]][d_site[0]] -= value

        # InterAll
        if self.param_ham["InterAll"] is not None:
            site = np.zeros(4, int)
            for site_info, value in self.param_ham["InterAll"].items():
                for i in range(4):
                    site[i] = site_info[2 * i] + site_info[2 * i + 1] * self.Nsize
                # Diagonal Fock term
                self.Ham_local[site[0]][site[1]][site[2]][site[3]] += value
                self.Ham_local[site[2]][site[3]][site[0]][site[1]] += value
                # OffDiagonal Fock term
                self.Ham_local[site[0]][site[3]][site[2]][site[1]] -= value
                self.Ham_local[site[2]][site[1]][site[0]][site[3]] -= value

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

        Ene["InterAll"] -= np.dot(green_local.T, np.dot(self.Ham_local, green_local)) / 2.0

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
            rest += abs(self.Green[site1][site2] - self.Green
            self.Green[site1][site2] = self.Green_old[site1][site2] * (1.0 - mix) + mix * self.Green[site1][site2]
        self.physics["Rest"] = np.sqrt(rest) / (2.0 * self.Nsize * self.Nsize)
        self.Green[np.where(abs(self.Green) < 1e-12)] = 0

    def get_results(self):
        return (self.physics, self.Green)

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


class UHF_SzFix(UHF):
    def __init__(self, param, param_ham, output, totalsz):
        super(UHF_SzFix, self).__init__(param, param_ham, output)
        Ncond = self.Ncond
        TwoSz = totalsz * 2
        self.green_list = {
            "up": {"label": [i for i in range(self.Nsize)], "value": 0.5, "occupied": int((self.Ncond + TwoSz) / 2)},
            "down": {"label": [i for i in range(self.Nsize, 2 * self.Nsize)], "value": -0.5,
                     "occupied": int((self.Ncond - TwoSz) / 2)}}

