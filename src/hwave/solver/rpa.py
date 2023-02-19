from __future__ import annotations
from typing import Optional

import sys, os
import numpy as np
import numpy.fft as FFT
import itertools
from requests.structures import CaseInsensitiveDict

from .perf import do_profile

import logging
logger = logging.getLogger(__name__)

#import read_input_k
import hwave.qlmsio.read_input_k as read_input_k

#XXX
#debug = True
debug = False

class Lattice:
    """
    Lattice parameters:

      Lattice.cellshape = (Lx, Ly, Lz)
      Lattice.cellvol
        original box shape and its volume

      Lattice.subshape = (Bx, By, Bz)
      Lattice.subvol
        shape of supercell and its volume

      Lattice.shape = (Nx, Ny, Nz)
      Lattice.nvol
        box shape in units of supercell and its volume

      Lattice.has_sublattice = True | False
        whether sublattice (except equiv to original lattice) is defined

    Constructor:
      Lattice(param) where param contains
        CellShape (mandatory)
        SubShape  (optional: default to same as CellShape)

    Note:
      Bx, By, Bz must be chosen so that Lx, Ly, Lz are divisable 
      by Bx, By, Bz, respectively. otherwise, initialization fails.

    """
    def __init__(self, params):
        self._init_lattice(params)
        self._show_params()

    def _init_lattice(self, params):
        if "CellShape" not in params:
            logger.error("Lattice initialization failed: 'CellShape' not found.")
            sys.exit(1)

        cell = params.get("CellShape")
        if type(cell) is not list:
            cell = [ cell ]
        cell_len = len(cell)
        if cell_len < 1 or cell_len > 3:
            logger.error("dimension of CellShape must be one, two, or three.")
            sys.exit(1)
        if cell_len < 3:
            cell.extend([1] * (3 - cell_len))

        Lx,Ly,Lz = cell
        self.cellshape = (Lx,Ly,Lz)
        self.cellvol = Lx * Ly * Lz
        self.celldim = cell_len

        if self.cellvol == 0:
            logger.error("invalid CellShape.")
            sys.exit(1)

        subcell = params.get("SubShape", [Lx,Ly,Lz])
        if type(subcell) is not list:
            subcell = [ subcell ]
        if len(subcell) != cell_len:
            logger.error("dimension of SubShape does not match with that of CellShape.")
            sys.exit(1)
        if len(subcell) < 3:
            subcell.extend([1] * (3 - len(subcell)))

        Bx,By,Bz = subcell
        self.subshape = (Bx,By,Bz)
        self.subvol = Bx * By * Bz

        if self.subvol == 0:
            logger.error("invalid SubShape.")
            sys.exit(1)

        self.has_sublattice = (self.subvol > 1)

        # check consistency
        # XXX use reciprocal lattice
        if not all([ self.cellshape[i] % self.subshape[i] == 0 for i in range(3) ]):
            logger.error("SubShape is not compatible with CellShape.")
            sys.exit(1)

        # replace by lattice of supercells
        nx, ny, nz = Lx//Bx, Ly//By, Lz//Bz
        nvol = nx * ny * nz

        self.shape = (nx, ny, nz)
        self.nvol = nvol

    def _show_params(self):
        logger.info("Lattice parameters:")
        logger.info("    CellShape       = {}".format(self.cellshape))
        logger.info("    cell volume     = {}".format(self.cellvol))
        logger.info("    cell dimension  = {}".format(self.celldim))
        logger.info("    SubShape        = {}".format(self.subshape))
        logger.info("    subshape volume = {}".format(self.subvol))
        logger.info("    Shape           = {}".format(self.shape))
        logger.info("    shape volume    = {}".format(self.nvol))
        logger.info("    has_sublattice  = {}".format(self.has_sublattice))

class Interaction:
    """
    Construct Hamiltonian from input
    """
    def __init__(self, lattice, param_ham):
        self.lattice = lattice
        self.param_ham = param_ham

        self._init_interaction()

        self.norb = param_ham["Geometry"]["norb"]

        self._make_ham_trans()
        self._make_ham_inter()

        pass

    def _init_interaction(self):
        # reinterpret interaction coefficient on sublattice
        if self.lattice.has_sublattice:
            # backup
            self.param_ham_orig = copy.deepcopy(self.param_ham)

            # replace by sublatticed versions
            for type in self.param_ham.keys():
                if type in ["Initial"]:
                    pass
                elif type in ["Geometry"]:
                    tbl = self._reshape_geometry(self.param_ham[type])
                    self.param_ham[type] = tbl
                else:
                    tbl = self._reshape_interaction(self.param_ham[type])
                    self.param_ham[type] = tbl
        pass

    def _reshape_geometry(self, geom):
        Bx,By,Bz = self.lattice.subshape
        bvol = self.lattice.subvol

        norb = geom['norb']

        geom_new = {}
        geom_new['norb'] = geom['norb'] * bvol
        geom_new['rvec'] = np.matmul(np.diag([Bx, By, Bz]), geom['rvec'])

        sc = np.array([1.0/Bx, 1.0/By, 1.0/Bz])
        cw = [ sc * geom['center'][k] for k in range(norb) ]

        centerv = np.zeros((norb * bvol, 3), dtype=np.double)
        k = 0
        for bz,by,bx in itertools.product(range(Bz),range(By),range(Bx)):
            for i in range(norb):
                centerv[k] = cw[i] + np.array([bx, by, bz]) * sc
                k += 1
        geom_new['center'] = centerv

        return geom_new

    def _reshape_interaction(self, ham):
        Bx,By,Bz = self.lattice.subshape
        nx,ny,nz = self.lattice.shape

        norb_orig = self.param_ham_orig["Geometry"]["norb"]

        def _reshape_orbit(a, x):
            return a + norb_orig * ( x[0] + Bx * (x[1] + By * (x[2])))

        def _round(x, n):
            return x % n if x >= 0 else x % -n

        ham_new = {}
        for (irvec,orbvec), v in ham.items():
            rx,ry,rz = irvec
            alpha,beta = orbvec

            for bz,by,bx in itertools.product(range(Bz),range(By),range(Bx)):

                # original cell index of both endes
                #   x0 -> x1=x0+r
                x0,y0,z0 = bx, by, bz
                x1,y1,z1 = x0 + rx, y0 + ry, z0 + rz

                # decompose into supercell-index and cell-index within supercell
                #   x0 = 0 + x0
                #   x1 = X + xr
                xx1,xr1 = x1 // Bx, x1 % Bx
                yy1,yr1 = y1 // By, y1 % By
                zz1,zr1 = z1 // Bz, z1 % Bz

                # find orbital index within supercell
                aa = _reshape_orbit(alpha,(x0,y0,z0))
                bb = _reshape_orbit(beta, (xr1,yr1,zr1))

                # check wrap-around: maybe overwritten by duplicate entries
                ir = (_round(xx1, nx), _round(yy1, ny), _round(zz1, nz))
                ov = (aa, bb)

                ham_new[(ir, ov)] = v

        return ham_new

    def _export_interaction(self, type, file_name):
        intr = self.param_ham[type]

        min_r = [0,0,0]
        max_r = [0,0,0]
        for (irvec,orbvec), v in self.param_ham[type].items():
            for k in range(3):
                min_r[k] = irvec[k] if irvec[k] < min_r[k] else min_r[k]
                max_r[k] = irvec[k] if irvec[k] > max_r[k] else max_r[k]
        rshape = [ max_r[i]-min_r[i]+1 for i in range(3) ]
        rsize = rshape[0] * rshape[1] * rshape[2]

        with open(file_name, "w") as fw:
            # write header
            fw.write("{} with sublattice for uhfk\n".format(type))
            # write number of orbitals
            fw.write("{}\n".format(self.norb))
            # write number of points of box enclosing transport vectors
            fw.write("{}\n".format(rsize))
            # write multiplicity factors (nominal)
            for i in range(rsize):
                if i > 0 and i % 15 == 0:
                    fw.write("\n")
                fw.write(" 1")
            fw.write("\n")
            # write index and elements
            for (irvec,orbvec), v in self.param_ham[type].items():
                if (abs(v) > 1.0e-12):
                    fw.write("{:3} {:3} {:3} {:3} {:3}  {:.12f} {:.12f}\n".format(
                        *irvec, orbvec[0]+1, orbvec[1]+1, v.real, v.imag
                    ))

    def _make_ham_trans(self):
        nx,ny,nz = self.lattice.shape
        nvol = self.lattice.nvol
        norb = self.norb
        ns = 2
        nd = norb * ns

        if 'Transfer' not in self.param_ham.keys():
            logger.warn("Transfer not found")
            self.ham_trans_r = None
            self.ham_trans_q = None
            return

        tab_r = np.zeros((nx,ny,nz,norb,norb), dtype=complex)

        for (irvec,orbvec), v in self.param_ham["Transfer"].items():
            tab_r[(*irvec,*orbvec)] = v

        # Fourier transform
        tab_q = FFT.ifftn(tab_r, axes=(0,1,2)) * nvol

        # 2x2 unit matrix for spin dof
        spin = np.eye(2)

        # T_{a,s,b,t}(r)
        ham_r = np.einsum('rab,st->rsatb',
                          tab_r.reshape(nvol,norb,norb),
                          spin
                          ).reshape(nx,ny,nz,nd,nd)

        # T_{a,s,b,t}(k)
        ham_q = np.einsum('kab,st->ksatb',
                          tab_q.reshape(nvol,norb,norb),
                          spin
                          ).reshape(nx,ny,nz,nd,nd)

        logger.debug("ham_trans_r shape={}, size={}".format(ham_r.shape, ham_r.size))
        logger.debug("ham_trans_r nonzero count={}".format(ham_r[abs(ham_r) > 1.0e-8].size))
        
        logger.debug("ham_trans_q shape={}, size={}".format(ham_q.shape, ham_q.size))
        logger.debug("ham_trans_q nonzero count={}".format(ham_q[abs(ham_q) > 1.0e-8].size))

        self.ham_trans_r = ham_r
        self.ham_trans_q = ham_q

    def _make_ham_inter(self):
        nx,ny,nz = self.lattice.shape
        nvol = self.lattice.nvol
        norb = self.norb
        ns = 2
        nd = norb * ns

        # Interaction Hamiltonian W[r,b,bp,a,ap]
        #   H = W(r)^{\beta\beta^\prime\alpha\alpah^\prime}
        #        * c_{i\alpha}^\dagger c_{i\alpha^\prime} c_{j\beta^\prime}^\dagger c_{j\beta}
        ham_r = np.zeros((nx,ny,nz,*(ns,norb)*4), dtype=complex)

        # spin(a,ap,bp,b)  0: up, 1: down
        spin_table = {
            'CoulombIntra': { (0,0,1,1): 1, (1,1,0,0): 1 },
            'CoulombInter': { (0,0,0,0): 1, (1,1,1,1): 1, (0,0,1,1): 1, (1,1,0,0): 1 },
            'Hund':         { (0,0,0,0): 1, (1,1,1,1): 1 },
            'Ising':        { (0,0,0,0): 1, (1,1,1,1): 1, (0,0,1,1): -1, (1,1,0,0): -1 },
            'PairLift':     { (0,1,0,1): 1, (1,0,1,0): 1 },
            'Exchange':     { (0,1,1,0): 1, (1,0,0,1): 1 },
            #--
            'PairHop':      { (0,0,1,1): 1, (1,1,0,0): 1 },
        }

        # coulomb-type interactions
        def _append_inter(type):
            logger.debug("_append_inter {}".format(type))
            spins = spin_table[type]
            for (irvec,orbvec), v in self.param_ham[type].items():
                a, b = orbvec
                for spinvec, w in spins.items():
                    s1,s2,s3,s4 = spinvec
                    # beta beta' alpha alpha'
                    orb = (s4, b, s3, b, s1, a, s2, a)
                    ham_r[(*irvec, *orb)] += v * w

        # pairhop-type interaction
        #   H^PH = P^{\alpha\alpha^\prime}
        #        * c_{i\alpha\up}^\dagger c_{j\alpha^\prime\up}
        #            c_{i\alpha\down}^\dagger c_{j\alpha^\prime\down}
        #        + (up <-> down)
        def _append_pairhop(type):
            spins = spin_table[type]
            for (irvec,orbvec), v in self.param_ham[type].items():
                # take account of same-site interaction only
                if (irvec == (0,0,0)):
                    a, b = orbvec
                    for spinvec, w in spins.items():
                        s1,s2,s3,s4 = spinvec
                        # beta beta' alpha alpha'
                        orb = (s4, a, s3, b, s1, a, s2, b)
                        ham_r[(*irvec, *orb)] += v * w

        if 'CoulombIntra' in self.param_ham.keys():
            _append_inter('CoulombIntra')

        if 'CoulombInter' in self.param_ham.keys():
            _append_inter('CoulombInter')

        if 'Hund' in self.param_ham.keys():
            _append_inter('Hund')

        if 'Ising' in self.param_ham.keys():
            _append_inter('Ising')

        if 'PairLift' in self.param_ham.keys():
            _append_inter('PairLift')

        if 'Exchange' in self.param_ham.keys():
            _append_inter('Exchange')

        if 'PairHop' in self.param_ham.keys():
            _append_pairhop('PairHop')

        # reshape to W(r)^{bb'aa'}, a,b=(spin,alpha)
        ham_r = ham_r.reshape(nx,ny,nz,*(nd,)*4)

        logger.debug("ham_inter_r shape={}, size={}".format(ham_r.shape, ham_r.size))
        logger.debug("ham_inter_r nonzero count={}".format(ham_r[abs(ham_r) > 1.0e-8].size))

        # Fourier transform W(q)^{bb'aa'}
        ham_q = FFT.ifftn(ham_r, axes=(0,1,2)) * nvol

        logger.debug("ham_inter_q shape={}, size={}".format(ham_q.shape, ham_q.size))
        logger.debug("ham_inter_q nonzero count={}".format(ham_q[abs(ham_q) > 1.0e-8].size))

        self.ham_inter_r = ham_r
        self.ham_inter_q = ham_q

class RPA:
    """
    RPA calculation
    """
    @do_profile
    def __init__(self, param_ham, info_log, info_mode):
        self.param_ham = param_ham
        self.info_log = info_log
        self.param_mod = CaseInsensitiveDict(info_mode.get("param", {}))

        self.lattice = Lattice(self.param_mod)

        self._init_param()
        self._show_params()

        self.ham_info = Interaction(self.lattice, param_ham)

        pass

    def _init_param(self):
        logger.debug(">>> RPA._init_param")

        self.T = self.param_mod.get("T", 0.0)
        self.ene_cutoff = self.param_mod.get("ene_cutoff", 1.0e+2)

        self.nmat = self.param_mod.get("Nmat", 1024)

        self.norb = self.param_ham["geometry"]["norb"]
        self.ns = 2  # spin dof
        self.nd = self.norb * self.ns

        self.Nstate = self.lattice.nvol * self.nd

        if "Ncond" in self.param_mod:
            self.Ncond = self.param_mod["Ncond"]
            self.filling = 1.0 * self.Ncond / self.Nstate
        elif "Nelec" in self.param_mod:
            self.Ncond = self.param_mod["Nelec"]
            self.filling = 1.0 * self.Ncond / self.Nstate
        elif "filling" in self.param_mod:
            self.filling = self.param_mod['filling']
            self.Ncond = self.Nstate * self.filling
            # force Ncond to be integer when T=0
            if self.T == 0.0:
                round_mode = self.param_mod.get("Ncond_round_mode", "strict")
                self.Ncond = self._round_to_int(self.Ncond, round_mode)
        else:
            logger.error("Ncond or Nelec or filling is missing. abort")
            sys.exit(1)

        self.mu_value = self.param_mod.get("mu", None)

        pass

    def _round_to_int(self, val, mode):
        import math
        mode = mode.lower()  # case-insensitive
        if   mode == "as-is":
            ret = val  # not rounding to int
        elif mode == "round":
            ret = round(val)
        elif mode == "round-up":
            ret = math.ceil(val)
        elif mode == "round-down":
            ret = math.floor(val)
        elif mode == "round-to-zero":
            ret = int(val)
        elif mode == "round-off":
            nn = math.floor(val)
            rr = val - nn
            ret = nn if rr < 0.5 else nn+1
        elif mode == "strict":
            if val != round(val):
                logger.error("value not integer")
                sys.exit(1)
            ret = round(val)
        elif mode == "exact":  # "round" with warning
            if val != round(val):
                logger.warning("value not integer")
            ret = round(val)
        else:
            logger.error("round_to_int: unknown mode {}".format(mode))
            sys.exit(1)
        return ret

    def _show_params(self):
        logger.debug(">>> RPA._show_params")
        logger.info("RPA parameters:")
        logger.info("    norbit          = {}".format(self.norb))
        logger.info("    nspin           = {}".format(self.ns))
        logger.info("    nd              = {}".format(self.nd))
        logger.info("    Nmat            = {}".format(self.nmat))

        logger.info("    Ncond           = {}".format(self.Ncond))
        logger.info("    filling         = {:.3f}".format(self.filling))
        logger.info("    T               = {}".format(self.T))
        logger.info("    E_cutoff        = {:e}".format(self.ene_cutoff))

        if self.mu_value is not None:
            logger.info("    mu              = {}".format(self.mu_value))

        # logger.info("    RndSeed         = {}".format(self.param_mod["RndSeed"]))
        # logger.info("    strict_hermite  = {}".format(self.strict_hermite))
        # logger.info("    hermite_tol     = {}".format(self.hermite_tolerance))
        pass

    @do_profile
    def solve(self, green_info, path_to_output):
        logger.info("Start RPA calculations")

        beta = 1.0/self.T

        if "chi0q" in green_info and green_info["chi0q"] is not None:
            # use chi0q input
            chi0q = green_info["chi0q"]
            if chi0q.shape[0] != self.nmat:
                logger.info("nmat resized from {} to {}".format(self.nmat, chi0q.shape[0]))
                self.nmat = chi0q.shape[0]
        else:
            self._calc_epsilon_k(green_info)

            if self.mu_value is not None:
                mu = self.mu_value
            else:
                dist, mu = self._find_mu(self.Ncond, self.T)

            green0 = self._calc_green(beta, mu)

            chi0q = self._calc_chi0q(green0, beta)

            green_info["chi0q"] = chi0q

        # solve
        sol = self._solve_rpa(chi0q, self.ham_info.ham_inter_q)

        # adhoc store
        green_info["chiq"] = sol

        logger.info("End RPA calculations")
        pass

    @do_profile
    def save_results(self, info_outputfile, green_info):
        logger.info("Save RPA results")
        path_to_output = info_outputfile["path_to_output"]

        freq_range = info_outputfile.get("matsubara_frequency", "all")
        freq_index = self._find_index_range(freq_range)
        logger.debug("freq_index = {}".format(freq_index))

        if "chiq" in info_outputfile.keys():
            file_name = os.path.join(path_to_output, info_outputfile["chiq"])
            np.savez(file_name,
                     chiq = green_info["chiq"][freq_index],
                     freq_index = freq_index,
                     )
            logger.info("save_results: save chiq in file {}".format(file_name))

        if "chi0q" in info_outputfile.keys():
            file_name = os.path.join(path_to_output, info_outputfile["chi0q"])
            np.savez(file_name,
                     chi0q = green_info["chi0q"][freq_index],
                     freq_index = freq_index,
                     )
            logger.info("save_results: save chi0q in file {}".format(file_name))
        
        pass

    def _find_index_range(self, freq_range):
        nmat = self.nmat

        if type(freq_range) == int:
            # e.g. freq_range = 0
            freq_index = [ freq_range ]
        elif type(freq_range) == list:
            if len(freq_range) == 1:
                # e.g. freq_range = [index]
                freq_index = [ i for i in freq_range ]
            elif len(freq_range) == 2:
                # e.g. freq_range = [min,max]
                freq_index = [ i for i in range(freq_range[0], freq_range[1]+1) ]
            elif len(freq_range) >= 3:
                # e.g. freq_range = [min,max,step]
                freq_index = [ i for i in range(freq_range[0], freq_range[1]+1, freq_range[2]) ]
            else:
                raise ValueError("invalid value for matsubara_frequency")
        elif type(freq_range) == str:
            freq_range = freq_range.lower()
            if freq_range == "all":
                freq_index = [ i for i in range(nmat) ]
            elif freq_range == "center":
                freq_index = [ nmat//2 ]
            elif freq_range == "none":
                freq_index = []
            else:
                raise ValueError("invalid value for matsubara_frequency")
        else:
            raise ValueError("invalid value type for matsubara_frequency")

        return freq_index

    @do_profile
    def read_init(self, info_inputfile):
        logger.debug("RPA read initial configs")
        info = {}

        path_to_input = info_inputfile.get("path_to_input", "")

        if "chi0q_init" in info_inputfile.keys():
            file_name = os.path.join(path_to_input, info_inputfile["chi0q_init"])
            info["chi0q"] = self._read_chi0q(file_name)

        return info

    def _read_chi0q(self, file_name):
        logger.debug("RPA._read_chi0q")

        try:
            logger.debug("read chi0q from {}".format(file_name))
            data = np.load(file_name)
            chi0q = data["chi0q"]
            logger.debug("chi0q: shape={}".format(chi0q.shape))
            assert len(chi0q.shape) == 6, "unexpected shape: {}".format(chi0q.shape)
        except Exception as e:
            logger.error("read_chi0q failed: {}".format(e))
            sys.exit(1)

        # check size
        cs = chi0q.shape
        try:
            assert cs[1] == self.lattice.nvol, "lattice volume"
            assert cs[2] == self.nd, "shape[2]"
            assert cs[3] == self.nd, "shape[3]"
            assert cs[4] == self.nd, "shape[4]"
            assert cs[5] == self.nd, "shape[5]"
        except AssertionError as e:
            logger.error("unexpected data size {}".format(e))
            sys.exit(1)

        return chi0q

    @do_profile
    def _calc_epsilon_k(self, green_info):
        logger.debug(">>> RPA._calc_epsilon_k")

        nvol = self.lattice.nvol
        nd = self.nd

        # input green function
        # g0 = _init_green(green_info)

        # H0(k) = T_ab(k) + (H * G0)_ab(k)
        H0 = self.ham_info.ham_trans_q
        # XXX assume g0 = 0

        # diagonalize H0(k)
        w,v = np.linalg.eigh(H0.reshape(nvol,nd,nd))

        self.H0_eigenvalue = w
        self.H0_eigenvector = v

    @do_profile
    def _find_mu(self, Ncond, T):
        from scipy import optimize

        # load eigenvalues and eigenvectors
        w = self.H0_eigenvalue
        v = self.H0_eigenvector
        # fetch parameters
        # Ncond = self.Ncond
        # T = self.T
        ene_cutoff = self.ene_cutoff

        ev = np.sort(w.flatten())
        occupied_number = Ncond

        def _fermi(t, mu, ev):
            w_ = (ev - mu) / t
            mask_ = w_ < ene_cutoff
            w1_ = np.where( mask_, w_, 0.0 )
            v1_ = 1.0 / (1.0 + np.exp(w1_))
            v_ = np.where( mask_, v1_, 0.0 )
            return v_

        def _calc_delta_n(mu):
            ff = _fermi(T, mu, w)
            nn = np.einsum('kal,kl,kal->', np.conjugate(v), ff, v)
            return nn.real - occupied_number

        # find mu s.t. <n>(mu) = N0
        is_converged = False
        if (_calc_delta_n(ev[0]) * _calc_delta_n(ev[-1])) < 0.0:
            logger.debug("RPA._find_mu: try bisection")
            mu, r = optimize.bisect(_calc_delta_n, ev[0], ev[-1], full_output=True, disp=False)
            is_converged = r.converged
        if not is_converged:
            logger.debug("RPA._find_mu: try newton")
            mu, r = optimize.newton(_calc_delta_n, ev[0], full_output=True)
            is_converged = r.converged
        if not is_converged:
            logger.error("RPA._find_mu: not converged. abort")
            sys.exit(1)

        logger.info("RPA._find_mu: mu = {}".format(mu))

        dist = _fermi(T, mu, w)

        return dist, mu

    @do_profile
    def _calc_green(self, beta, mu):
        """
        ew: eigenvalues  ew[k,i]    i-th eigenvalue of wavenumber k
        ev: eigenvectors ev[k,a,i]  i-th eigenvector corresponding to ew[k,i]
        beta: inverse temperature
        mu: chemical potential
        """
        logger.debug(">>> RPA._calc_green")

        # load eigenvalues and eigenvectors
        ew = self.H0_eigenvalue
        ev = self.H0_eigenvector

        nx,ny,nz = self.lattice.shape
        nvol = self.lattice.nvol
        nd = self.nd
        nmat = self.nmat

        iomega = (np.arange(nmat) * 2 + 1 - nmat) * np.pi / beta

        # 1 / (iw_{n} - (e_i(k) - mu)) -> g[n,k,i]
        g = 1.0 / (np.tile(1j * iomega, (nd,nvol,1)).T - np.tile((ew - mu), (nmat,1,1)))

        # G_ab^j(k,iw_n) = d_{a,j} d_{b,j}^* / (iw_{n} - (e_j(k) - mu))
        green = np.einsum('kaj,kbj,lkj->jlkab', ev, np.conj(ev), g)

        return green

    @do_profile
    def _calc_chi0q(self, green_kw, beta):
        # green function
        #   green_kw[j,l,k,a,b]: G_ab^j(k,ie)
        #     l : matsubara freq i\epsilon_l = i\pi(2*l+1-N)/\beta for l=0..N-1
        #     k : wave number kx,ky,kz
        #     a,b : orbital and spin index a=(s,alpha) b=(t,beta)
        #     j : diagnoalized orbital index

        logger.debug(">>> RPA._calc_chi0q")

        nx,ny,nz = self.lattice.shape
        nvol = self.lattice.nvol
        nd = self.nd
        nmat = self.nmat

        # Fourier transform from wave number space to coordinate space
        green_rw = FFT.ifftn(green_kw.reshape(nd,nmat,nx,ny,nz,nd*nd), axes=(2,3,4))

        # Fourier transform from matsubara freq to imaginary time
        omg = np.exp(-1j * np.pi * (1.0/nmat - 1.0) * np.arange(nmat))

        green_rt = np.einsum('jtv,t->jtv',
                             FFT.fft(green_rw.reshape(nd,nmat,nvol*nd*nd), axis=1),
                             omg).reshape(nd,nmat,nx,ny,nz,nd,nd)

        # calculate chi0(r,t)[a,a',b,b'] = G(r,t)[a,b] * G(-r,-t)[b',a']
        green_rev = np.flip(np.roll(green_rt, -1, axis=(1,2,3,4)), axis=(1,2,3,4)).reshape(nd,nmat,nvol,nd,nd)

        sgn = np.full(nmat, -1)
        sgn[0] = 1

        chi0_rt = np.einsum('jlrab,jlrdc,l->lracbd',
                            green_rt.reshape(nd,nmat,nvol,nd,nd),
                            green_rev,
                            sgn)

        # Fourier transform to wave number space
        chi0_qt = FFT.fftn(chi0_rt.reshape(nmat,nx,ny,nz,nd**4), axes=(1,2,3))

        # Fourier transform to matsubara freq
        omg = np.exp(1j * np.pi * (-1) * np.arange(nmat))

        chi0_qw = FFT.ifft(
            np.einsum('tv,t->tv', chi0_qt.reshape(nmat,nvol*nd**4), omg),
            axis=0).reshape(nmat,nvol,nd,nd,nd,nd) * (-1.0/beta)

        return chi0_qw

    @do_profile
    def _solve_rpa(self, chi0q, ham):
        nvol = self.lattice.nvol
        nmat = self.nmat
        nd = self.nd
        ndx = nd**2  # combined index a = (alpha, alpha')

        # 1 - X^0(l,k,aa,bb) W(k,bb,cc)
        mat  = np.tile(np.eye(ndx, dtype=np.complex128), (nmat,nvol,1,1))
        mat += np.einsum('lkab,kbc->lkac',
                         chi0q.reshape(nmat,nvol,ndx,ndx),
                         ham.reshape(nvol,ndx,ndx))
        
        # [ 1 - X^0 W ]^-1
        matrev = np.linalg.inv(mat)

        # [ 1 - X^0 W ]^-1 X^0
        sol = np.einsum('lkab,lkbc->lkac',
                        matrev,
                        chi0q.reshape(nmat,nvol,ndx,ndx))

        return sol.reshape(nmat,nvol,nd,nd,nd,nd)

    @do_profile
    def _calc_chi0q_exact(self, eps_k, beta, mu):
        nx,ny,nz = self.lattice.shape

        fermi = 1.0 / (1.0 + np.exp(beta*(eps_k - mu)))

        eta = 1.0e-12

        chi0q_exact = np.zeros((nx,ny,nz), dtype=complex)
        for iqx,iqy,iqz in itertools.product(range(nx),range(ny),range(nz)):
            for ikx,iky,ikz in itertools.product(range(nx),range(ny),range(nz)):
                kqx = (ikx + iqx) % nx
                kqy = (iky + iqy) % ny
                kqz = (ikz + iqz) % nz

                chi0q_exact[iqx,iqy,iqz] += -( fermi[kqx,kqy,kqz] - fermi[ikx,iky,ikz] ) / (eps_k[kqx,kqy,kqz] - eps_k[ikx,iky,ikz] + 1j * eta)

        return (chi0q_exact/(nx*ny*nz)).real

    @do_profile
    def _get_green(self, eps_k, beta, mu):
        nmat = self.nmat
        nvol = self.nvol

        iomega = (np.arange(nmat)*2+1-nmat) * np.pi / beta
        eps = eps_k.reshape(nvol)

        # g = np.zeros((nmat,nvol), dtype=complex)
        # for i in range(nvol):
        #     for l in range(nmat):
        #         g[l][i] = 1.0 / (1j * iomega[l] - (eps[i] - mu))

        g = 1.0 / (np.tile(1j * iomega, (nvol,1)).T - np.tile((eps-mu), (nmat,1)))

        return g

def run(*, input_dict: Optional[dict] = None, input_file: Optional[str] = None):
    logger = logging.getLogger("run")

    if input_dict is None:
        raise RuntimeError("input_dict not passed")

    # Initialize information about log
    info_log = input_dict.get("log", {})
    info_log["print_level"] = info_log.get("print_level", 1)
    info_log["print_step"] = info_log.get("print_step", 1)

    # Initialize information about mode
    info_mode = input_dict.get("mode", {})
    info_file = input_dict.get("file", {"input": {}, "output": {}})
    # Initialize information about input files
    info_inputfile = info_file.get("input", {})
    info_inputfile["path_to_input"] = info_inputfile.get("path_to_input", "")

    # Initialize information about output files
    info_outputfile = info_file.get("output", {})
    info_outputfile["path_to_output"] = info_outputfile.get("path_to_output", "output")
    path_to_output = info_outputfile["path_to_output"]
    os.makedirs(path_to_output, exist_ok=True)

    if "mode" not in info_mode:
        logger.error("mode is not defined in [mode].")
        sys.exit(1)

    mode = info_mode["mode"]

    if mode == "RPA":
        logger.info("RPA mode")

        logger.info("Read interaction definitions from files")
        read_io = read_input_k.QLMSkInput(info_inputfile)
        ham_info = read_io.get_param("ham")

        solver = RPA(ham_info, info_log, info_mode)

        green_info = read_io.get_param("green")
        green_info.update( solver.read_init(info_input_file) )

        logger.info("Start RPA calculation")
        solver.solve(green_info, path_to_output)

        logger.info("Save calculation results")
        solver.save_results(info_outputfile, green_info)

        logger.info("all done.")

    elif mode == "RPAtest":
        logger.info("RPA test mode")

        beta = 10.0
        mu = 0.0

        nx = 2**5
        ny = 2**5
        nz = 1

        kx_array = np.linspace(0, 2*np.pi, nx, endpoint=False)
        ky_array = np.linspace(0, 2*np.pi, ny, endpoint=False)
        kz_array = np.linspace(0, 2*np.pi, nz, endpoint=False)

        t = 1.0
        t1 = 0.5

        eps_k = np.zeros((nx,ny,nz), dtype=complex)
        for ix,iy,iz in itertools.product(range(nx),range(ny),range(nz)):
            kx=kx_array[ix]
            ky=ky_array[iy]
            kz=kz_array[iz]

            eps_k[ix,iy,iz] = 2.0 * t * (np.cos(kx)+np.cos(ky)+np.cos(kz)) + 2.0 * t1 * np.cos(kx+ky+kz)

        rpa = RPA({}, info_log, info_mode)
        x = rpa._calc_chi0q_exact(eps_k, beta, mu)

        with open("test_chi0q_exact.dat", "w") as fw:
            for iqx in range(nx):
                for iqy in range(ny):
                    fw.write("{} {} {}\n".format(kx_array[iqx],ky_array[iqy],x[iqx,iqy,0]))

        # g = rpa._get_green(eps_k, beta, mu)

        # y = rpa._calc_chi0q(g, beta)

        # with open("test_chi0q.dat", "w") as fw:
        #     for iqx in range(nx):
        #         for iqy in range(ny):
        #             fw.write("{} {} {}\n".format(kx_array[iqx],ky_array[iqy],(y[0,:,0,0,0,0].reshape(nx,ny,nz))[iqx,iqy,0].real))

        #---
        # param_ham = {}
        # solver = RPA(param_ham, info_log, info_mode)

        # logger.info("Start RPA calculation")
        # solver.solve(green_info, path_to_output)

        # logger.info("Save calculation results")
        # solver.save_results(info_outputfile, green_info)

        logger.info("all done.")

    else:
        logger.warning("mode is incorrect: mode={}.".format(mode))
        sys.exit(0)

    pass

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="increase verbosity", action="count", default=0)
    parser.add_argument("-q", "--quiet", help="decrease verbosity", action="count", default=0)
    parser.add_argument("input_file", nargs='?', default="input.toml", help="parameter file in TOML format")
    args = parser.parse_args()

    fmt = "%(asctime)s %(levelname)s %(name)s: %(message)s"
    log_level = 20 - (args.verbose - args.quiet) * 10
    logging.basicConfig(level=log_level, format=fmt)

    try:
        import toml
        params = toml.load(args.input_file)
    except Exception as e:
        print(e)
        sys.exit(1)

    run(input_dict=params)
