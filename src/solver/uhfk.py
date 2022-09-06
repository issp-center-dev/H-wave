import itertools
import logging
import numpy as np
import os
from .base import solver_base

logger = logging.getLogger("qlms").getChild("uhfk")

class UHFk(solver_base):

    def __init__(self, param_ham, info_log, info_mode, param_mod=None):
        super().__init__(param_ham, info_log, info_mode, param_mod)

        self.init_lattice()
        self.init_orbit()
        self.init_param()

        #self.dump_param_ham()

        self.init_interaction()

        #self.dump_param_ham()

        self.show_param()

        # local data
        self.physics = {
            "Ene": { "Total": 0.0, "Band": 0.0 },
            "NCond": 0.0,
            "Sz": 0.0,
            "Rest": 1.0
        }

        nvol = self.nvol
        nd = self.nd
        
        self._green_list = {
            "eigenvalue":  np.zeros((nvol,nd), dtype=np.complex128),
            "eigenvector": np.zeros((nvol,nd,nd), dtype=np.complex128)
        }

        #self.Green = np.zeros((nvol,ns,norb,ns,norb), dtype=np.complex128),

        # work area
        self.ham = np.zeros((nvol,nd,nd), dtype=np.complex128)

    def init_lattice(self):
        # check parameters
        if not "CellShape" in self.param_mod:
            logger.info("CellShape is missing. abort")
            exit(1)
            
        Lx,Ly,Lz = self.param_mod.get("CellShape")
        self.cellshape = (Lx,Ly,Lz)

        Bx,By,Bz = self.param_mod.get("SubShape", [1,1,1])
        self.subshape = (Bx,By,Bz)
        self.subvol = Bx * By * Bz

        # check consistency
        # XXX use reciprocal lattice
        err = 0
        for i in range(3):
            if self.cellshape[i] % self.subshape[i] != 0:
                err += 1
        if err > 0:
            logger.info("SubShape is not compatible with CellShape")
            logger.info("abort.")
            exit(1)

        if any([ _ > 1 for _ in self.subshape]):
            logger.info("enable sublattice")
            self.has_sublattice = True
        else:
            self.has_sublattice = False

        # replace by lattice of supercells
        nx, ny, nz = Lx//Bx, Ly//By, Lz//Bz
        nvol = nx * ny * nz

        self.shape = (nx, ny, nz)
        self.nvol = nvol

    def init_orbit(self):
        norb = self.param_ham["Geometry"]["norb"]
        ns = 2  # spin dof

        # take account of supercell
        self.norb_orig = norb
        self.norb = norb * self.subvol

        self.nd = self.norb * ns
        self.ns = ns

    def init_param(self):
        Ncond = self.param_mod["Ncond"]
        self.Ncond = Ncond

        self.T = self.param_mod.get("T", 0.0)
        self.ene_cutoff = self.param_mod.get("ene_cutoff", 1e+2)

    def show_param(self):
        logger.info("Show parameters")
        logger.info("    Cell Shape  = {}".format(self.cellshape))
        logger.info("    Sub Shape   = {}".format(self.subshape))
        logger.info("    Block       = {}".format(self.shape))
        logger.info("    Block volume= {}".format(self.nvol))
        logger.info("    Num orbit   = {}".format(self.norb_orig))
        logger.info("    Num orbit eff= {}".format(self.norb))
        logger.info("    nspin       = {}".format(self.ns))
        logger.info("    nd          = {}".format(self.nd))

        logger.info("    Ncond       = {}".format(self.Ncond))
        logger.info("    T           = {}".format(self.T))
        logger.info("    E_cutoff    = {}".format(self.ene_cutoff))

        logger.info("    Mix         = {}".format(self.param_mod["Mix"]))
        logger.info("    RndSeed     = {}".format(self.param_mod["RndSeed"]))
        logger.info("    IterationMax= {}".format(self.param_mod["IterationMax"]))
        logger.info("    EPS         = {}".format(self.param_mod["EPS"]))

    def init_interaction(self):
        # reinterpret interaction coefficient on sublattice
        if self.has_sublattice:
            Bx,By,Bz = self.subshape
            
            for type in self.param_ham.keys():
                if type in ["Geometry", "Initial"]:
                    pass
                else:
                    tbl = {}
                    for (irvec,orbvec), v in self.param_ham[type].items():

                        rx,ry,rz = irvec
                        Rx = rx // Bx
                        Ry = ry // By
                        Rz = rz // Bz
                        rx0 = rx % Bx
                        ry0 = ry % By
                        rz0 = rz % Bz
                        rr = rx0 + Bx * (ry0 + By * (rz0))

                        alpha,beta = orbvec

                        irvec2 = (Rx,Ry,Rz)
                        orbvec2 = (alpha, beta + rr * self.norb_orig)

                        tbl[(irvec2,orbvec2)] = v

                    # replace
                    self.param_ham[type] = tbl

    def dump_param_ham(self):
        for type in self.param_ham.keys():
            if type in ["Geometry"]:
                print("type =", type)
                print(self.param_ham[type])
            elif type in ["Initial"]:
                print("type =", type)
            else:
                print("type =", type)

                for (irvec,orbvec), v in self.param_ham[type].items():
                    print("\t",irvec,orbvec," = ",v)
                

    def solve(self, path_to_output):
        print_level = self.info_log["print_level"]

        logger.info("Start UHFk calculations")

        if print_level > 0:
            logger.info("step, rest, energy, NCond, Sz")

        self._make_ham_trans() # T_ab(k)
        self._make_ham_inter() # J_ab(r)

        self.Green = self._initial_green()

        is_converged = False

        for i_step in range(self.param_mod["IterationMax"]):
            self._make_ham()
            self._diag()
            self._green()
            self._calc_energy()
            self._calc_phys()  # do update in calc_phys
            if print_level > 0 and i_step % self.info_log["print_step"] == 0:
                logger.info(
                    "{}, {:.8g}, {:.8g}, {:.4g}, {:.4g} "
                    .format(i_step,
                            self.physics["Rest"],
                            self.physics["Ene"]["Total"],
                            self.physics["NCond"],
                            self.physics["Sz"]))
            if self.physics["Rest"] < self.param_mod["eps"]:
                is_converged = True
                break
            
        if is_converged:
            logger.info("UHFk calculation succeeded: rest={}, eps={}."
                        .format(self.physics["Rest"], self.param_mod["eps"]))
        else:
            logger.info("UHFk calculation failed: rest={}, eps={}."
                        .format(self.physics["Rest"], self.param_mod["eps"]))
            
    def _initial_green(self):
        logger.info(">>> _initial_green")

        nx,ny,nz = self.shape
        nvol     = self.nvol
        norb     = self.norb
        nd       = self.nd
        ns       = self.ns

        # green function G_{ab,st}(r) : gab_r(r,s,a,t,b)

        if self.param_ham["Initial"] is not None:
            file_name = self.param_ham["Initial"]
            logger.info("read green function from file {}".format(file_name))
            green = self._read_green(file_name)

        else:
            logger.info("initialize green function with random numbers")
            
            np.random.seed(self.param_mod["RndSeed"])
            rand = np.random.rand(nvol * nd * nd).reshape(nvol,ns,norb,ns,norb)

            green = 0.01 * (rand - 0.5)

        return green

    def _make_ham_trans(self):
        logger.info(">>> _make_ham_trans")
        
        nx,ny,nz = self.shape
        nvol     = self.nvol
        norb     = self.norb
        nd       = self.nd

        # data structure of T_ab(r): T(rx,ry,rz,a,b)
        tab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)

        for (irvec,orbvec), v in self.param_ham["Transfer"].items():
            tab_r[(*irvec, *orbvec)] = v

        # check hermiticity
        t = np.conjugate(
            np.transpose(
                np.flip(np.roll(tab_r, -1, axis=(0,1,2)), axis=(0,1,2)),
                (0,1,2,4,3)
            )
        )
        if not np.allclose(t, tab_r):
            logger.info("hermiticity check failed: |T_ba(-r)^* - T_ab(r)| = {}".format(np.sum(np.abs(t - tab_r))) )

        # fourier transform
        tab_k = np.fft.ifftn(tab_r, axes=(0,1,2), norm='forward')

        # 2x2 unit matrix to introduce spin index
        spin = np.eye(2)

        # T_{a,si,b,sj}(k)
        ham = np.einsum(
            'kab, st -> ksatb', tab_k.reshape(nvol,norb,norb), spin
        ).reshape(nvol,nd,nd)

        # store
        self.ham_trans = ham

    def _make_ham_inter(self):
        logger.info(">>> _make_ham_inter")

        nx,ny,nz = self.shape
        nvol     = self.nvol
        norb     = self.norb
        nd       = self.nd

        #----------------
        # interaction table
        #----------------
        self.inter_table = {}
        self.spin_table = {}
        
        #----------------
        # Coulomb Intra and Coulomb Inter
        #----------------
        if 'Coulomb' in self.param_ham.keys():
            # assume zvo_ur.dat
            # divide into r=0 (coulomb intra) and r!=0 (coulomb inter)

            # coulomb intra: diagonal part of r=0 cell
            uab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)

            # coulomb inter: off-diagonal part
            vab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)

            for (irvec,orbvec), v in self.param_ham["Coulomb"].items():
                alpha, beta = orbvec
                if irvec == (0,0,0) and alpha == beta:
                    uab_r[(*irvec,alpha,beta)] = v
                else:
                    vab_r[(*irvec,alpha,beta)] = v

            # coulomb intra
            # interaction coeffs
            self.inter_table["CoulombIntra"] = uab_r # r=0 component
            # spin combination
            self.spin_table["CoulombIntra"] = np.zeros((2,2,2,2), dtype=int)
            self.spin_table["CoulombIntra"][0,1,1,0] = 1
            self.spin_table["CoulombIntra"][1,0,0,1] = 1

            #XXX
            # n_up n_up = n_up -> include in one-body term

            # coulomb inter
            # J~ab(r) = Jab(r) + Jba(-r)
            vba = np.conjugate(
                np.transpose(
                    np.flip(np.roll(vab_r, -1, axis=(0,1,2)), axis=(0,1,2)),
                    (0,1,2,4,3)
                )
            )

            # interaction coeffs
            self.inter_table["CoulombInter"] = vab_r + vba
            # spin combination
            self.spin_table["CoulombInter"] = np.zeros((2,2,2,2), dtype=int)
            self.spin_table["CoulombInter"][0,0,0,0] = 1
            self.spin_table["CoulombInter"][1,1,1,1] = 1
            self.spin_table["CoulombInter"][0,1,1,0] = 1
            self.spin_table["CoulombInter"][1,0,0,1] = 1

        else:
            self.inter_table["CoulombIntra"] = None
            self.inter_table["CoulombInter"] = None

            # coulomb intra
            if 'CoulombIntra' in self.param_ham.keys():
                uab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)

                # only r=0 and a=b component
                for (irvec,orbvec), v in self.param_ham["CoulombIntra"].items():
                    alpha, beta = orbvec
                    if irvec == (0,0,0) and alpha == beta:
                        uab_r[(*irvec, *orbvec)] = v

                # interaction coeffs
                self.inter_table["CoulombIntra"] = uab_r
                # spin combination
                self.spin_table["CoulombIntra"] = np.zeros((2,2,2,2), dtype=int)
                self.spin_table["CoulombIntra"][0,1,1,0] = 1
                self.spin_table["CoulombIntra"][1,0,0,1] = 1

            if 'CoulombInter' in self.param_ham.keys():
                vab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)

                for (irvec,orbvec), v in self.param_ham["CoulombInter"].items():
#                    if irvec != (0,0,0):
                        vab_r[(*irvec, *orbvec)] = v

                vba = np.conjugate(
                    np.transpose(
                        np.flip(np.roll(vab_r, -1, axis=(0,1,2)), axis=(0,1,2)),
                        (0,1,2,4,3)
                    )
                )

                # interaction coeffs
                self.inter_table["CoulombInter"] = vab_r + vba
                # spin combination
                self.spin_table["CoulombInter"] = np.zeros((2,2,2,2), dtype=int)
                self.spin_table["CoulombInter"][0,0,0,0] = 1
                self.spin_table["CoulombInter"][1,1,1,1] = 1
                self.spin_table["CoulombInter"][0,1,1,0] = 1
                self.spin_table["CoulombInter"][1,0,0,1] = 1

        #----------------
        # Hund
        #----------------        
        if 'Hund' in self.param_ham.keys():
            jab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)

            for (irvec,orbvec), v in self.param_ham["Hund"].items():
                jab_r[(*irvec, *orbvec)] = v

            # J~ab(r) = Jab(r) + Jba(-r)
            jba = np.conjugate(
                np.transpose(
                    np.flip(np.roll(jab_r, -1, axis=(0,1,2)), axis=(0,1,2)),
                    (0,1,2,4,3)
                )
            )

            # interaction coeffs : -J^{Hund} by convention
            self.inter_table["Hund"] = -(jab_r + jba)
            # spin combination
            self.spin_table["Hund"] = np.zeros((2,2,2,2), dtype=int)
            self.spin_table["Hund"][0,0,0,0] = 1
            self.spin_table["Hund"][1,1,1,1] = 1
        else:
            self.inter_table["Hund"] = None

        #----------------
        # Ising
        #----------------
        if 'Ising' in self.param_ham.keys():
            jab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)
        
            for (irvec,orbvec), v in self.param_ham["Ising"].items():
                jab_r[(*irvec, *orbvec)] = v

            # J~ab(r) = Jab(r) + Jba(-r)
            jba = np.conjugate(
                np.transpose(
                    np.flip(np.roll(jab_r, -1, axis=(0,1,2)), axis=(0,1,2)),
                    (0,1,2,4,3)
                )
            )

            # interaction coeffs
            self.inter_table["Ising"] = jab_r + jba
            # spin combination
            self.spin_table["Ising"] = np.zeros((2,2,2,2), dtype=int)
            self.spin_table["Ising"][0,0,0,0] = 1
            self.spin_table["Ising"][1,1,1,1] = 1
            self.spin_table["Ising"][0,1,1,0] = -1
            self.spin_table["Ising"][1,0,0,1] = -1
        else:
            self.inter_table["Ising"] = None

        #----------------
        # PairLift
        #----------------
        if 'PairLift' in self.param_ham.keys():
            jab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)

            for (irvec,orbvec), v in self.param_ham["PairLift"].items():
                jab_r[(*irvec, *orbvec)] = v

            # J~ab(r) = Jab(r) + Jba(-r)
            jba = np.conjugate(
                np.transpose(
                    np.flip(np.roll(jab_r, -1, axis=(0,1,2)), axis=(0,1,2)),
                    (0,1,2,4,3)
                )
            )

            # interaction coeffs
            self.inter_table["PairLift"] = jab_r + jba
            # spin combination
            self.spin_table["PairLift"] = np.zeros((2,2,2,2), dtype=int)
            self.spin_table["PairLift"][0,0,1,1] = 1
            self.spin_table["PairLift"][1,1,0,0] = 1
        else:
            self.inter_table["PairLift"] = None

        #----------------
        # Exchange
        #----------------
        if 'Exchange' in self.param_ham.keys():
            jab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)

            for (irvec,orbvec), v in self.param_ham["Exchange"].items():
                jab_r[(*irvec, *orbvec)] = v

            # J~ab(r) = Jab(r) + Jba(-r)
            jba = np.conjugate(
                np.transpose(
                    np.flip(np.roll(jab_r, -1, axis=(0,1,2)), axis=(0,1,2)),
                    (0,1,2,4,3)
                )
            )

            # interaction coeffs : -J^{Ex} by convention
            self.inter_table["Exchange"] = -(jab_r + jba)
            # spin combination
            self.spin_table["Exchange"] = np.zeros((2,2,2,2), dtype=int)
            self.spin_table["Exchange"][0,1,0,1] = 1
            self.spin_table["Exchange"][1,0,1,0] = 1
        else:
            self.inter_table["Exchange"] = None

        #----------------
        # PairHop
        #----------------
        if 'PairHop' in self.param_ham.keys():
            jab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)

            for (irvec,orbvec), v in self.param_ham["PairHop"].items():
                jab_r[(*irvec, *orbvec)] = v

            # J~ab(r) = Jab(r) + Jba(-r)
            jba = np.conjugate(
                np.transpose(
                    np.flip(np.roll(jab_r, -1, axis=(0,1,2)), axis=(0,1,2)),
                    (0,1,2,4,3)
                )
            )

            # interaction coeffs
            self.inter_table["PairHop"] = jab_r + jba
            # spin combination
            self.spin_table["PairHop"] = np.zeros((2,2,2,2), dtype=int)
            self.spin_table["PairHop"][0,1,1,0] = 1
            self.spin_table["PairHop"][1,0,0,1] = 1
        else:
            self.inter_table["PairHop"] = None

    def _make_ham(self):
        logger.info(">>> _make_ham")

        nx,ny,nz = self.shape
        nvol     = self.nvol
        nd       = self.nd
        norb     = self.norb
        ns       = self.ns

        # green function G_{ab,st}(r) : gab_r(r,s,a,t,b)
        gab_r = self.Green

        # diagonal part G_{bb,st}(r=0) : gbb(s,t,b)
        gbb = np.diagonal(gab_r, axis1=2, axis2=4)[0,:,:,:]

        # hamiltonian H_{ab,st}(k) : ham(k,(s,a),(t,b))
        ham = np.zeros((nvol,nd,nd), dtype=np.complex128)

        # transfer term  -T_{ab}(k)
        logger.info("Transfer")
        ham -= self.ham_trans
        
        # interaction term: Coulomb type
        for type in ['CoulombIntra', 'CoulombInter', 'Hund', 'Ising', 'PairLift', 'Exchange']:
            if self.inter_table[type] is not None:
                logger.info(type)

                # coefficient of interaction term J_{ab}(r)
                jab_r = self.inter_table[type].reshape(nvol,norb,norb)
                # and its spin combination  Spin(s1,s2,s3,s4)
                spin = self.spin_table[type]

                # non-cross term
                #   sum_r J_{ab}(r) G_{bb,uv}(0) Spin{s,u,v,t}
                hh0 = np.einsum('uvb, suvt -> stb', gbb, spin)
                #hh1 = np.einsum('rab, uvb, suvt -> rsta', jab_r, gbb, spin)
                hh1 = np.einsum('rab, stb -> rsta', jab_r, hh0)
                hh2 = np.einsum('rsta, ab -> rsatb', hh1, np.eye(norb, norb))
                hh3 = np.sum(hh2, axis=0).reshape(nd,nd)

                ham += np.broadcast_to(hh3, (nvol,nd,nd))
                
                # cross term
                #   - sum_r J_{ab}(r) G_{ba,uv}(r) Spin{s,u,t,v} e^{ikr}
                hh4 = np.einsum('rab, rubva, sutv -> rsatb', jab_r, gab_r, spin)

                #   fourier transform: sum_r (*) e^{ikr}
                hh5 = np.fft.ifftn(hh4.reshape(nx,ny,nz,nd,nd), axes=(0,1,2), norm='forward')

                ham -= hh5.reshape(nvol, nd, nd)

        # for type in ['CoulombIntra']:
        #     if self.inter_table[type] is not None:
        #         logger.info(type)

        #         # coefficient U_a
        #         ua = self.inter_table[type]
        #         # spin combination
        #         spin = self.spin_table[type]

        #         # non-cross term
        #         #   sum_r U_{a} G_{aa,uv}(0) Spin{s,u,v,t}
        #         hh = np.einsum('a, uva, suvt -> ast', ua, gbb, spin)
        #         hh0 = np.transpose(
        #             np.broadcast_to(hh, (norb,norb,ns,ns)),
        #             (2,0,3,1)).reshape(nd,nd)

        #         ham += np.broadcast_to(hh0, (nvol,nd,nd))

        #         # cross term
        #         #   - sum_r U_{a} G_{aa,uv}(0) Spin{s,u,t,v}
        #         hh = np.einsum('a, uva, sutv -> ast', ua, gbb, spin)
        #         hh0 = np.transpose(
        #             np.broadcast_to(hh, (norb,norb,ns,ns)),
        #             (2,0,3,1)).reshape(nd,nd)

        #         ham -= np.broadcast_to(hh0, (nvol,nd,nd))

        # interaction term: PairHop type
        for type in ['PairHop']:
            if self.inter_table[type] is not None:
                logger.info(type)

                # coefficient of interaction term J_{ab}(r)
                jab_r = self.inter_table[type].reshape(nvol,norb,norb)
                # and its spin combination  Spin(s1,s2,s3,s4)
                spin = self.spin_table[type]

                # green function G_{ab,st}(-r)
                gab_mr = np.flip(
                    np.roll(
                        gab_r.reshape(nx,ny,nz,ns,norb,ns,norb), -1, axis=(0,1,2)
                    ), axis=(0,1,2)
                ).reshape(nvol,ns,norb,ns,norb)

                # # non-cross term
                # #   sum_r J_{ab}(r) G_{ab,uv}(-r) Spin{s,u,v,t} e^{ikr}
                # hh = np.einsum('rab, ruavb, suvt -> rsatb', jab_r, gab_mr, spin)
                # #   fourier transform: sum_r (*) e^{ikr}
                # hh2 = np.fft.ifftn(hh.reshape(nx,ny,nz,nd,nd), axes=(0,1,2), norm='forward')

                # ham += hh2.reshape(nvol, nd, nd)
                
                # # cross term
                # #   - sum_r J_{ab}(r) G_{ab,uv}(-r) Spin{s,u,t,v} e^{ikr}
                # hh = np.einsum('rab, ruavb, sutv -> rsatb', jab_r, gab_mr, spin)
                # #   fourier transform: sum_r (*) e^{ikr}
                # hh2 = np.fft.ifftn(hh.reshape(nx,ny,nz,nd,nd), axes=(0,1,2), norm='forward')

                # ham -= hh2.reshape(nvol, nd, nd)

                hh1 = np.einsum('ruavb, suvt -> rsatb', gab_mr, spin)
                hh2 = np.einsum('ruavb, sutv -> rsatb', gab_mr, spin)
                hh3 = np.einsum('rab, rsatb -> rsatb', jab_r, (hh1 - hh2))
                hh4 = np.fft.ifftn(hh3.reshape(nx,ny,nz,nd,nd), axes=(0,1,2), norm='forward')

                ham += hh4.reshape(nvol, nd, nd)

        # store
        self.ham = ham

    def _diag(self):
        logger.info(">>> _diag")

        # hamiltonian H_{ab,st}(k) : ham(k,(s,a),(t,b))
        mat = self.ham

        # diagonalize
        w,v = np.linalg.eigh(mat)

        # store
        self._green_list["eigenvalue"] = w
        self._green_list["eigenvector"] = v

    def _green(self):
        logger.info(">>> _green")

        nx,ny,nz = self.shape
        nvol     = self.nvol
        nd       = self.nd
        norb     = self.norb
        ns       = self.ns

        w = self._green_list["eigenvalue"]
        v = self._green_list["eigenvector"]

        if self.T == 0:
            # fill lowest Ncond eigenmodes
            ev = np.sort(w.flatten())
            ev0 = ev[self.Ncond]

            dist = np.where(w < ev0, 1.0, 0.0)
        else:
            from scipy import optimize

            def _fermi(t, mu, ev):
                w = (ev - mu) / t
                v = np.where( w > self.ene_cutoff, 0.0, 1.0 / (1.0 + np.exp(w)) )
                return v

            def _calc_delta_n(mu):
                ff = _fermi(self.T, mu, w)
                nn = np.einsum('kal,kl,kal->', np.conjugate(v), ff, v)
                return nn.real - occupied_number

            ev = np.sort(w.flatten())
            occupied_number = self.Ncond

            # find mu s.t. <n>(mu) = N0
            is_converged = False

            if (_calc_delta_n(ev[0]) * _calc_delta_n(ev[-1])) < 0.0:
                logger.info("+++ find mu: try bisection")
                mu, r = optimize.bisect(_calc_delta_n, ev[0], ev[-1], full_output=True, disp=False)
                is_converged = r.converged
            if not is_converged:
                logger.info("+++ find mu: try newton")
                mu, r = optimize.newton(_calc_delta_n, ev[0], full_output=True)
                is_converged = r.converged
            if not is_converged:
                logger.info("+++ find mu: not converged. abort")
                exit(1)

            self._green_list["mu"] = mu

            logger.info("mu = {}".format(mu))

            dist = _fermi(self.T, mu, w)

        # G_ab(k) = sum_l v_al(k)^* v_bl(k) where ev_l(k) < ev0
        gab_k = np.einsum('kal, kl, kbl -> kab', np.conjugate(v), dist, v)

        gab_k = gab_k.reshape(nx,ny,nz,nd,nd)

        # G_ab(r) = 1/V sum_k G_ab(k) e^{-ikr}
        gab_r = np.fft.fftn(gab_k, axes=(0,1,2), norm='forward')

        # store
        self.Green_prev = self.Green
        self.Green = gab_r.reshape(nvol,ns,norb,ns,norb)

    def _calc_phys(self):
        logger.info(">>> _calc_phys")

        nx,ny,nz = self.shape
        nvol     = self.nvol
        nd       = self.nd
        norb     = self.norb
        ns       = self.ns

        # expectation value of Ncond
        gab_r = self.Green.reshape(nvol,nd,nd)

        n = np.sum(np.diagonal(gab_r, axis1=1, axis2=2))
        self.physics["NCond"] = n.real

        logger.info("ncond = {}".format(n))

        # expectation value of Sz
        gab_r = self.Green

        sz = 0.0
        sz += np.sum(np.diagonal(gab_r[:,0,:,0,:], axis1=1, axis2=2))
        sz -= np.sum(np.diagonal(gab_r[:,1,:,1,:], axis1=1, axis2=2))

        self.physics["Sz"] = 0.5 * sz.real

        logger.info("sz = {}".format(sz))

        # residue
        rest = np.linalg.norm(self.Green - self.Green_prev)
        self.physics["Rest"] = rest / np.size(self.Green) * 2

        logger.info("rest = {}".format(rest))

        # update
        mix = self.param_mod["Mix"]

        g_new = (1.0 - mix) * self.Green_prev + mix * self.Green

        tol = self.param_mod.get("tolerance", 1.0e-12)
        g_new[np.where(abs(g_new) < tol)] = 0.0

        self.Green = g_new

    def _calc_energy(self):
        logger.info(">>> _calc_energy")

        nx,ny,nz = self.shape
        nvol     = self.nvol
        nd       = self.nd
        norb     = self.norb
        ns       = self.ns

        occupied_number = self.Ncond

        energy = {}
        energy_total = 0.0

        if self.T == 0:
            ev = np.sort(self._green_list["eigenvalue"].flatten())
            energy["Band"] = np.sum(ev[:occupied_number])
            logger.info("energy: Band = {}".format(energy["Band"]))
            energy_total += energy["Band"]
        else:
            w = self._green_list["eigenvalue"]
            v = self._green_list["eigenvector"]
            mu = self._green_list["mu"]
            T = self.T

            def _fermi(t, mu, ev):
                w = (ev - mu) / t
                v = np.where( w > self.ene_cutoff, 0.0, 1.0 / (1.0 + np.exp(w)) )
                return v

            ln_e = np.log1p(np.exp(-(w - mu) / T))

            nn = np.einsum('kal,kl,kal->', np.conjugate(v), _fermi(T, mu, w), v)

            e_band = mu * nn - T * np.sum(ln_e)

            energy["Band"] = e_band.real
            energy_total += energy["Band"]
            logger.info("energy: Band = {}".format(e_band))

        for type in self.inter_table.keys():
            if self.inter_table[type] is not None:
                jab_r = self.inter_table[type].reshape(nvol,norb,norb)
                spin  = self.spin_table[type]
                gab_r = self.Green

                if type == "PairHop":
                    w1 = np.einsum('stuv, rsavb, rtaub -> rab', spin, gab_r, gab_r)
                    w2 = np.einsum('stuv, rsaub, rtavb -> rab', spin, gab_r, gab_r)
                    ee = np.einsum('rab, rab ->', jab_r, w1-w2)
                    energy[type] = -ee/2.0

                else:
                    w1 = np.einsum('stuv, rvasa, rubtb -> rab', spin, gab_r, gab_r)
                    w2 = np.einsum('stuv, rubsa, rvatb -> rab', spin, gab_r, gab_r)
                    ee = np.einsum('rab, rab->', jab_r, w1-w2)
                    energy[type] = -ee/2.0

                energy_total += energy[type].real
                logger.info("energy: {} = {}".format(type, energy[type]))
            else:
                # logger.info("energy: {} skip".format(type))
                pass

        energy["Total"] = energy_total
        self.physics["Ene"] = energy.copy()

    def get_results(self):
        return (self.physics, self.Green)

    def save_results(self, info_outputfile, green_info):
        path_to_output = info_outputfile["path_to_output"]

        if "energy" in info_outputfile.keys():
            file_name = os.path.join(path_to_output, info_outputfile["energy"])

            with open(file_name, "w") as fw:
                type_ex = ["Total", "Band"]

                for type in type_ex:
                    fw.write("Energy_{:<12} = {}\n".format(type, self.physics["Ene"][type].real))

                for type in self.physics["Ene"].keys():
                    if type in type_ex:
                        pass
                    else:
                        fw.write("Energy_{:<12s} = {}\n".format(type, self.physics["Ene"][type].real))

                fw.write("NCond   = {}\n".format(self.physics["NCond"]))
                fw.write("Sz      = {}\n".format(self.physics["Sz"]))

            logger.info("save_results: save energy in file {}".format(file_name))

        if "eigen" in info_outputfile.keys():
            file_name = os.path.join(path_to_output, info_outputfile["eigen"])
            np.savez(file_name,
                     eigenvalue  = self._green_list["eigenvalue"],
                     eigenvector = self._green_list["eigenvector"])
            logger.info("save_results: save eigenvalues and eigenvectors in file {}".format(file_name))
                
        if "green" in info_outputfile.keys():
            file_name = os.path.join(path_to_output, info_outputfile["green"])
            self._save_green(file_name)
            logger.info("save_results: save green function in file {}".format(file_name))

        if "initial" in info_outputfile.keys():
            #XXX
            logger.info("save_results: save initial is not supported")
            pass

    def _read_green(self, file_name):
        try:
            v = np.load(file_name)
        except FileNotFoundError:
            logger.info("_read_green: file {} not found".format(file_name))
            return None
            
        if "green" in v.files:
            data = v["green"]
        else:
            data = None
        return data

    def _save_green(self, file_name):
        np.savez(file_name, green = self.Green)


