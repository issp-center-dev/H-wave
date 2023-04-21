from __future__ import annotations
from typing import Optional

import sys, os
import numpy as np
import numpy.fft as FFT
import itertools
import copy
from requests.structures import CaseInsensitiveDict

try:
    from .perf import do_profile
except ImportError:
    from functools import wraps
    def do_profile(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)
        return wrapper

import logging
logger = logging.getLogger(__name__)

#import read_input_k
import hwave.qlmsio.read_input_k as read_input_k

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
        logger.debug(">>> Lattice.__init__")

        self._init_lattice(params)
        self._show_params()

    def _init_lattice(self, params):
        logger.debug(">>> Lattice._init_lattice")

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
    def __init__(self, lattice, param_ham, info_mode):
        logger.debug(">>> Interaction.__init__")

        self.lattice = lattice
        self.param_ham = param_ham

        self._has_interaction = False
        self._has_interaction_coulomb = False
        self._has_interaction_exchange = False
        self._has_interaction_pairhop = False


        # mode options
        self.enable_spin_orbital = info_mode.get("enable_spin_orbital", False)

        # initialize, and reshape if use sublattice
        self._init_interaction()

        self.norb = param_ham["Geometry"]["norb"]

        # create hamiltonian
        self._make_ham_trans()
        self._make_ham_inter()

        pass

    def has_interaction(self):
        return self._has_interaction

    def has_interaction_coulomb(self):
        return self._has_interaction_coulomb

    def has_interaction_exchange(self):
        return self._has_interaction_exchange

    def has_interaction_pairhop(self):
        return self._has_interaction_pairhop

    def _init_interaction(self):
        logger.debug(">>> Interaction._init_interaction")

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
                elif type in ["Transfer"]:
                    tbl = self._reshape_interaction(self.param_ham[type], self.enable_spin_orbital)
                    self.param_ham[type] = tbl
                else:
                    tbl = self._reshape_interaction(self.param_ham[type], False)
                    self.param_ham[type] = tbl
        pass

    def _reshape_geometry(self, geom):
        logger.debug(">>> Interaction._reshape_geometry")

        Bx,By,Bz = self.lattice.subshape
        bvol = self.lattice.subvol

        norb = geom['norb']

        geom_new = {}
        geom_new['norb'] = geom['norb'] * bvol
        geom_new['rvec'] = np.matmul(np.diag([Bx, By, Bz]), geom['rvec'])

        sc = np.array([1.0/Bx, 1.0/By, 1.0/Bz])
        cw = [ sc * geom['center'][k] for k in range(norb) ]

        centerv = np.zeros((norb * bvol, 3), dtype=np.float64)
        k = 0
        for bz,by,bx in itertools.product(range(Bz),range(By),range(Bx)):
            for i in range(norb):
                centerv[k] = cw[i] + np.array([bx, by, bz]) * sc
                k += 1
        geom_new['center'] = centerv

        return geom_new

    def _reshape_interaction(self, ham, enable_spin_orbital):
        logger.debug(">>> Interaction._reshape_interaction")

        Bx,By,Bz = self.lattice.subshape
        nx,ny,nz = self.lattice.shape

        norb_orig = self.param_ham_orig["Geometry"]["norb"]

        def _reshape_orbit_(a, x):
            return a + norb_orig * ( x[0] + Bx * (x[1] + By * (x[2])))

        def _reshape_orbit_spin(a, x):
            a_, s_ = a%norb_orig, a//norb_orig
            return a_ + norb_orig * ( x[0] + Bx * (x[1] + By * (x[2] + Bz * s_)))

        if enable_spin_orbital:
            _reshape_orbit = _reshape_orbit_spin
        else:
            _reshape_orbit = _reshape_orbit_

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
        logger.debug(">>> Interaction._export_interaction")

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
        logger.debug(">>> Interaction._make_ham_trans")

        nx,ny,nz = self.lattice.shape
        nvol = self.lattice.nvol
        norb = self.norb
        ns = 2
        nd = norb * ns

        if 'Transfer' not in self.param_ham.keys():
            logger.warning("Transfer not found")
            self.ham_trans_r = None
            self.ham_trans_q = None
            self.ham_extern_r = None
            self.ham_extern_q = None
            return

        if self.enable_spin_orbital == True:
            # assume orbital index includes spin index
            tab_r = np.zeros((nx,ny,nz,nd,nd), dtype=np.complex128)

            for (irvec,orbvec), v in self.param_ham["Transfer"].items():
                tab_r[(*irvec,*orbvec)] = v

            # Fourier transform
            tab_q = FFT.ifftn(tab_r, axes=(0,1,2)) * nvol

            self.ham_trans_r = tab_r.reshape(nvol,nd,nd)
            self.ham_trans_q = tab_q.reshape(nvol,nd,nd)

        else:
            tab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)

            for (irvec,orbvec), v in self.param_ham["Transfer"].items():
                if orbvec[0] < norb and orbvec[1] < norb:
                    tab_r[(*irvec,*orbvec)] = v
                else:
                    pass  # skip spin dependence

            # Fourier transform
            tab_q = FFT.fftn(tab_r, axes=(0,1,2))

            # N.B. spin degree of freedom not included
            self.ham_trans_r = tab_r.reshape(nvol,norb,norb)
            self.ham_trans_q = tab_q.reshape(nvol,norb,norb)

        logger.debug("ham_trans_r shape={}, size={}, nonzero_count={}".format(
            self.ham_trans_r.shape,
            self.ham_trans_r.size,
            self.ham_trans_r[abs(self.ham_trans_r) > 1.0e-8].size,
        ))
        logger.debug("ham_trans_q shape={}, size={}, nonzero_count={}".format(
            self.ham_trans_q.shape,
            self.ham_trans_q.size,
            self.ham_trans_q[abs(self.ham_trans_q) > 1.0e-8].size,
        ))

        if 'Extern' in self.param_ham.keys():
            logger.info("read External field")

            hab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)

            for (irvec,orbvec), v in self.param_ham["Extern"].items():
                if orbvec[0] < norb and orbvec[1] < norb:
                    hab_r[(*irvec,*orbvec)] = v
                else:
                    pass  # skip spin dependence

            # Fourier transform
            hab_q = FFT.fftn(hab_r, axes=(0,1,2))

            # N.B. spin degree of freedom not included
            self.ham_extern_r = hab_r.reshape(nvol,norb,norb)
            self.ham_extern_q = hab_q.reshape(nvol,norb,norb)
        else:
            self.ham_extern_r = None
            self.ham_extern_q = None


    def _make_ham_inter(self):
        logger.debug(">>> Interaction._make_ham_inter")

        nx,ny,nz = self.lattice.shape
        nvol = self.lattice.nvol
        norb = self.norb
        ns = 2
        nd = norb * ns

        # Interaction Hamiltonian W[r,b,bp,a,ap]
        #   H = W(r)^{\beta\beta^\prime\alpha\alpah^\prime}
        #        * c_{i\alpha}^\dagger c_{i\alpha^\prime} c_{j\beta^\prime}^\dagger c_{j\beta}
        ham_r = np.zeros((nx,ny,nz,*(ns,norb)*4), dtype=np.complex128)

        # spin(a,ap,bp,b)  0: up, 1: down
        spin_table = {
            'CoulombIntra': { (0,0,1,1): 1, (1,1,0,0): 1 },
            'CoulombInter': { (0,0,0,0): 1, (1,1,1,1): 1, (0,0,1,1): 1, (1,1,0,0): 1 },
            'Hund':         { (0,0,0,0): -1, (1,1,1,1): -1 },
            'Ising':        { (0,0,0,0): 1, (1,1,1,1): 1, (0,0,1,1): -1, (1,1,0,0): -1 },
            'PairLift':     { (0,1,0,1): 1, (1,0,1,0): 1 },
            'Exchange':     { (0,1,1,0): -1, (1,0,0,1): -1 },
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
                        orb = (s4, b, s3, a, s1, a, s2, b)
                        ham_r[(*irvec, *orb)] += v * w

        if 'CoulombIntra' in self.param_ham.keys():
            _append_inter('CoulombIntra')
            self._has_interaction = True
            self._has_interaction_coulomb = True

        if 'CoulombInter' in self.param_ham.keys():
            _append_inter('CoulombInter')
            self._has_interaction = True
            self._has_interaction_coulomb = True

        if 'Hund' in self.param_ham.keys():
            _append_inter('Hund')
            self._has_interaction = True
            self._has_interaction_coulomb = True

        if 'Ising' in self.param_ham.keys():
            _append_inter('Ising')
            self._has_interaction = True
            self._has_interaction_coulomb = True

        if 'PairLift' in self.param_ham.keys():
            _append_inter('PairLift')
            self._has_interaction = True
            self._has_interaction_exchange = True

        if 'Exchange' in self.param_ham.keys():
            _append_inter('Exchange')
            self._has_interaction = True
            self._has_interaction_exchange = True

        if 'PairHop' in self.param_ham.keys():
            _append_pairhop('PairHop')
            self._has_interaction = True
            self._has_interaction_pairhop = True

        # reshape to W(r)^{bb'aa'}, a,b=(spin,alpha)
        ham_r = ham_r.reshape(nx,ny,nz,*(nd,)*4)

        logger.debug("ham_inter_r shape={}, size={}".format(ham_r.shape, ham_r.size))
        logger.debug("ham_inter_r nonzero count={}".format(ham_r[abs(ham_r) > 1.0e-8].size))

        # Fourier transform W(q)^{bb'aa'}
        ham_q = FFT.fftn(ham_r, axes=(0,1,2))

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
        logger.debug(">>> RPA.__init__")

        self.param_ham = param_ham
        self.info_log = info_log
        self.param_mod = CaseInsensitiveDict(info_mode.get("param", {}))

        self.lattice = Lattice(self.param_mod)
        self.ham_info = Interaction(self.lattice, param_ham, info_mode)

        self._set_scheme(info_mode)

        self._init_param()
        self._show_params()

        pass

    def _set_scheme(self, info_mode):
        # handle calc_scheme: must be called after setting up interactions

        self.calc_scheme = info_mode.get("calc_scheme", "auto")

        # auto choose
        if self.calc_scheme == "auto":
            if not self.ham_info.has_interaction():
                logger.error("calc_scheme must be specified for chi0q-only mode.")
                sys.exit(1)
            else:
                if self.ham_info.has_interaction_exchange():
                    self.calc_scheme = "squashed"
                    logger.info("auto mode for calc_scheme: set to squashed")
                else:
                    self.calc_scheme = "reduced"
                    logger.info("auto mode for calc_scheme: set to reduced")

        # consistency check
        if self.calc_scheme == "reduced" and self.ham_info.has_interaction_exchange():
            logger.error("calc_scheme=reduced is not compatible with exchange-type interaction.")
            sys.exit(1)

        # calc chiq if interaction term exists; otherwise chi0q-only mode
        self.calc_chiq = self.ham_info.has_interaction()

        # chi0q in reduced mode if calc_scheme is reduced or squashed
        self.enable_reduced = self.calc_scheme.lower() in ["reduced", "squashed"]
        
    def _init_param(self):
        logger.debug(">>> RPA._init_param")

        self.T = self.param_mod.get("T", 0.0)
        self.ene_cutoff = self.param_mod.get("ene_cutoff", 1.0e+2)

        self.nmat = self.param_mod.get("Nmat", 1024)

        self.norb = self.param_ham["geometry"]["norb"]
        self.ns = 2  # spin dof
        self.nd = self.norb * self.ns

        self.coeff_tail = self.param_mod.get("coeff_tail", 0.0)
        self.ext = self.param_mod.get("coeff_extern", 0.0)

        # exclusive options: mu and Ncond/filling
        have_mu = "mu" in self.param_mod.keys()
        have_Ncond = "Ncond" in self.param_mod.keys() or "Nelec" in self.param_mod.keys()
        have_filling = "filling" in self.param_mod.keys()

        if have_mu and (have_Ncond or have_filling):
            # conflicting parameters
            logger.error("both mu and Ncond or filling are specified")
            sys.exit(1)
        elif not (have_mu or have_Ncond or have_filling):
            # none specified
            logger.error("none of mu, Ncond, nor filling is specified")
            sys.exit(1)

        self.Nstate = self.lattice.nvol * self.nd

        if have_Ncond or have_filling:
            self.calc_mu = True
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

        if have_mu:
            self.calc_mu = False
            self.mu_value = self.param_mod.get("mu", None)

        # range of matsubara frequency in matrix calculation and output
        self.freq_range = self.param_mod.get("matsubara_frequency", "all")
        self.freq_index = self._find_index_range(self.freq_range)
        logger.debug("freq_index = {}".format(self.freq_index))

        # check parameters
        err = 0
        if self.T < 0.0:
            logger.error("T must be greater than or equal to zero: T={}".format(self.T))
            err += 1
        if self.calc_mu and self.Ncond <= 0:
            logger.error("Ncond must be greater than zero: Ncond={}".format(self.Ncond))
            err += 1
        if self.nmat <= 0:
            logger.error("Nmat must be greater than zero: Nmat={}".format(self.nmat))
            err += 1
        if err > 0:
            sys.exit(1)

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
        logger.info("    coeff_tail      = {}".format(self.coeff_tail))

        if self.calc_mu == True:
            logger.info("    Ncond           = {}".format(self.Ncond))
            logger.info("    filling         = {:.3f}".format(self.filling))
        else:
            logger.info("    mu              = {}".format(self.mu_value))

        logger.info("    T               = {}".format(self.T))
        logger.info("    E_cutoff        = {:e}".format(self.ene_cutoff))
        logger.info("    coeff_extern    = {}".format(self.ext))

        # logger.info("    RndSeed         = {}".format(self.param_mod["RndSeed"]))
        # logger.info("    strict_hermite  = {}".format(self.strict_hermite))
        # logger.info("    hermite_tol     = {}".format(self.hermite_tolerance))
        logger.info("    freq_range      = {}".format(self.freq_range))
        logger.info("    calc_chiq       = {}".format(self.calc_chiq))
        logger.info("    spin_orbital    = {}".format(self.ham_info.enable_spin_orbital))
        logger.info("    calc_scheme     = {}".format(self.calc_scheme))
        pass

    @do_profile
    def solve(self, green_info, path_to_output):
        logger.info("Start RPA calculations")

        beta = 1.0/self.T

        if "chi0q" in green_info and green_info["chi0q"] is not None:
            # use chi0q input
            chi0q = green_info["chi0q"]
            if chi0q.shape[0] != self.nmat:
                logger.info("partial range in matsubara frequency: {} in {}".format(chi0q.shape[0], self.nmat))
                #self.nmat = chi0q.shape[0]
        else:
            self._calc_epsilon_k(green_info)

            if self.calc_mu:
                if self.spin_mode == "spin-free":
                    Ncond = self.Ncond/2
                else:
                    Ncond = self.Ncond
                dist, mu = self._find_mu(Ncond, self.T)
            else:
                mu = self.mu_value

            green0, green0_tail = self._calc_green(beta, mu)
            #XXX
            self.green0 = green0
            self.green0_tail = green0_tail

            chi0q = self._calc_chi0q(green0, green0_tail, beta)

            # filter by matsubara freq range
            if len(self.freq_index) < self.nmat:
                chi0q = chi0q[:,self.freq_index]
                logger.info("filter range in matsubara frequency: {} in {}".format(chi0q.shape[0], self.nmat))

            # nblock
            if self.spin_mode in [ "spin-free", "spinful" ]:
                assert chi0q.shape[0] == 1
                chi0q = chi0q[0]
            else:
                assert chi0q.shape[0] == 2
                pass

            green_info["chi0q"] = chi0q

        if self.calc_chiq:

            if self.spin_mode == "spinful":
                chi0q_orig = chi0q
                ham_orig = self.ham_info.ham_inter_q

                if self.calc_scheme == "reduced":
                    # alpha=alpha', beta=beta' case
                    nvol = self.lattice.nvol
                    nd = self.nd
                    ham = np.einsum('kaabb->kab',
                                    ham_orig.reshape(nvol,*(nd,)*4)).reshape(nvol,*(nd,)*2)
                elif self.calc_scheme == "squashed":
                    logger.error("squash is not available with spin-orbital interaction")
                    sys.exit(1)
                else:
                    ham = ham_orig

            elif self.spin_mode == "spin-diag":
                chi0q_orig = chi0q
                ham_orig = self.ham_info.ham_inter_q

                if self.calc_scheme == "reduced":
                    nblock,nfreq,nvol,norb1,norb2 = chi0q_orig.shape

                    norb = self.norb
                    ns = self.ns
                    nd = norb * ns

                    spin_tensor = np.identity(2)
                    chi0q = np.einsum('glkab,gh->lkgahb',
                                      chi0q_orig,
                                      spin_tensor).reshape(nfreq,nvol,nd,nd)

                    ham = np.einsum('kaabb->kab',
                                    ham_orig.reshape(nvol,*(nd,)*4)).reshape(nvol,*(nd,)*2)

                elif self.calc_scheme == "squashed":
                    nblock,nfreq,nvol,norb1,norb2 = chi0q_orig.shape

                    norb = self.norb
                    ns = self.ns
                    nd = norb * ns

                    spin_tensor = np.zeros((2,2,2,2), dtype=np.int32)
                    spin_tensor[0,0,0,0] = 1
                    spin_tensor[1,1,1,1] = 1

                    chi0q = np.einsum('glkab,gtuv->lkgtauvb',
                                      chi0q_orig,
                                      spin_tensor).reshape(nfreq,nvol,ns,ns,norb,ns,ns,norb)

                    ham = np.einsum('ksauatbvb->ksuatvb',
                                    ham_orig.reshape(nvol,*(ns,norb)*4)).reshape(nvol,*(ns,ns,norb)*2)

                else:
                    nblock,nfreq,nvol,norb1,norb2,norb3,norb4 = chi0q_orig.shape

                    norb = self.norb
                    ns = self.ns
                    nd = norb * ns

                    spin_tensor = np.zeros((2,2,2,2), dtype=np.int32)
                    spin_tensor[0,0,0,0] = 1
                    spin_tensor[1,1,1,1] = 1

                    chi0q = np.einsum('glkabcd,gtuv->lkgatbucvd',
                                      chi0q_orig,
                                      spin_tensor).reshape(nfreq,nvol,nd,nd,nd,nd)
                    ham = ham_orig

            elif self.spin_mode == "spin-free":
                # introduce spin degree of freedom
                chi0q_orig = chi0q
                ham_orig = self.ham_info.ham_inter_q

                if self.calc_scheme == "reduced":
                    # alpha=alpha', beta=beta' case
                    nfreq,nvol,norb1,norb2 = chi0q_orig.shape

                    norb = self.norb
                    ns = self.ns
                    nd = norb * ns

                    spin_tensor = np.identity(ns)
                    chi0q = np.einsum('lkab,st->lksatb',
                                      chi0q_orig.reshape(nfreq,nvol,norb,norb),
                                      spin_tensor).reshape(nfreq,nvol,nd,nd)

                    ham = np.einsum('ksasatbtb->ksatb',
                                    ham_orig.reshape(nvol,*(ns,norb)*4)).reshape(nvol,*(nd,)*2)

                elif self.calc_scheme == "squashed":
                    # norb**2 squash
                    nfreq,nvol,norb1,norb2 = chi0q_orig.shape

                    norb = self.norb
                    ns = self.ns
                    nd = norb * ns

                    spin_tensor = np.zeros((2,2,2,2), dtype=np.int32)
                    spin_tensor[0,0,0,0] = 1
                    spin_tensor[1,1,1,1] = 1

                    chi0q = np.einsum('lkab,stuv->lkstauvb',
                                      chi0q_orig.reshape(nfreq,nvol,norb,norb),
                                      spin_tensor).reshape(nfreq,nvol,ns,ns,norb,ns,ns,norb)

                    ham = np.einsum('ksauatbvb->ksuatvb',
                                    ham_orig.reshape(nvol,*(ns,norb)*4)).reshape(nvol,*(ns,ns,norb)*2)

                else:
                    # general nd**4 case
                    nfreq,nvol,norb1,norb2,norb3,norb4 = chi0q_orig.shape

                    norb = self.norb
                    ns = self.ns
                    nd = norb * ns

                    spin_tensor = np.zeros((2,2,2,2), dtype=np.int32)
                    spin_tensor[0,0,0,0] = 1
                    spin_tensor[1,1,1,1] = 1

                    chi0q = np.einsum('lkabcd,stuv->lksatbucvd',
                                      chi0q_orig.reshape(nfreq,nvol,norb,norb,norb,norb),
                                      spin_tensor).reshape(nfreq,nvol,nd,nd,nd,nd)
                    ham = ham_orig

            # solve
            sol = self._solve_rpa(chi0q, ham)

            # adhoc store
            green_info["chiq"] = sol

        logger.info("End RPA calculations")
        pass

    @do_profile
    def save_results(self, info_outputfile, green_info):
        logger.info("Save RPA results")
        path_to_output = info_outputfile["path_to_output"]

        self._init_wavevec()

        if "chiq" in info_outputfile.keys():
            if self.calc_chiq == True:
                file_name = os.path.join(path_to_output, info_outputfile["chiq"])
                np.savez(file_name,
                         chiq = green_info["chiq"],
                         freq_index = self.freq_index,
                         wavevector_unit = self.kvec,
                         wavevector_index = self.wavenum_table,
                )
                logger.info("save_results: save chiq in file {}".format(file_name))
            else:
                logger.info("save_results: chiq not calculated. skip")

        if "chi0q" in info_outputfile.keys():
            file_name = os.path.join(path_to_output, info_outputfile["chi0q"])
            np.savez(file_name,
                     chi0q = green_info["chi0q"],
                     freq_index = self.freq_index,
                     wavevector_unit = self.kvec,
                     wavevector_index = self.wavenum_table,
                     )
            logger.info("save_results: save chi0q in file {}".format(file_name))

        pass

    def _init_wavevec(self):
        # wave vectors on sublatticed geometry
        def _klist(n):
            return np.roll( (np.arange(n)-(n//2)), -(n//2) )

        geom = self.param_ham["Geometry"]
        rvec = geom["rvec"]
        omg = np.dot(rvec[0], np.cross(rvec[1], rvec[2]))
        kvec = np.array([
            np.cross(rvec[(i+1)%3], rvec[(i+2)%3])/omg * 2*np.pi/self.lattice.shape[i]
            for i in range(3) ])

        self.kvec = kvec  # store reciprocal lattice vectors

        nx,ny,nz = self.lattice.shape
        nvol = self.lattice.nvol

        self.wavenum_table = np.array([(i,j,k) for i in _klist(nx) for j in _klist(ny) for k in _klist(nz)])

        wtable = np.zeros((nx,ny,nz,3), dtype=float)
        for ix, kx in enumerate(_klist(nx)):
            for iy, ky in enumerate(_klist(ny)):
                for iz, kz in enumerate(_klist(nz)):
                    v = kvec[0] * kx + kvec[1] * ky + kvec[2] * kz
                    wtable[ix,iy,iz] = v
        self.wave_table = wtable.reshape(nvol,3)

    def _find_index_range(self, freq_range):
        # decode matsubara frequency index list
        #   1. single index
        #   2. min, max, step
        #   3. keyword
        # note index n in [0 .. Nmat-1] corresponds to
        #   w_n = (2*n-Nmat) * pi / beta

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
            elif freq_range == "center" or freq_range == "zero":
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
        logger.info("RPA read initial configs")
        info = {}

        path_to_input = info_inputfile.get("path_to_input", "")

        if "chi0q_init" in info_inputfile.keys():
            file_name = os.path.join(path_to_input, info_inputfile["chi0q_init"])
            info["chi0q"] = self._read_chi0q(file_name)

        if "trans_mod" in info_inputfile.keys():
            file_name = os.path.join(path_to_input, info_inputfile["trans_mod"])
            info["trans_mod"] = self._read_trans_mod(file_name)

        if "green_init" in info_inputfile.keys():
            file_name = os.path.join(path_to_input, info_inputfile["green_init"])
            info["green_init"] = self._read_green(file_name)

        return info

    def _read_chi0q(self, file_name):
        logger.debug(">>> RPA._read_chi0q")

        try:
            logger.debug("read chi0q from {}".format(file_name))
            data = np.load(file_name)
            chi0q = data["chi0q"]
            logger.debug("chi0q: shape={}".format(chi0q.shape))
        except Exception as e:
            logger.error("read_chi0q failed: {}".format(e))
            sys.exit(1)

        # check size
        if self.calc_scheme == "general":
            if len(chi0q.shape) == 6:
                # spin-free or spinful
                #   shape = (nmat,nvol,nd,nd,nd,nd) where nd = norb or norb*nspin
                cs = chi0q.shape
                assert cs[1] == self.lattice.nvol, "lattice volume"
                nd = cs[2]
                assert nd == self.nd or nd == self.norb, "shape[2]"
                assert cs[3] == nd, "shape[3]"
                assert cs[4] == nd, "shape[4]"
                assert cs[5] == nd, "shape[5]"

                if nd == self.nd:
                    self.spin_mode = "spinful"
                else:
                    self.spin_mode = "spin-free"
            elif len(chi0q.shape) == 7:
                # spin-diagonal
                #   shape = (nblock,nmat,nvol,norb,norb,norb,norb)
                cs = chi0q.shape
                assert cs[0] == 2, "spin block"
                assert cs[2] == self.lattice.nvol, "lattice volume"
                nd = cs[3]
                assert nd == self.norb, "shape[3]"
                assert cs[4] == nd, "shape[4]"
                assert cs[5] == nd, "shape[5]"
                assert cs[6] == nd, "shape[6]"

                self.spin_mode = "spin-diag"
            else:
                assert False, "unexpected shape for general scheme: {}".format(chi0q.shape)

        elif self.calc_scheme == "reduced" or self.calc_scheme == "squashed":
            # reduced: shape = (nmat,nvol,nd,nd) where nd = norb or norb*nspin
            if len(chi0q.shape) == 4:
                # spin-free or spinful
                #   shape = (nmat,nvol,nd,nd) where nd = norb or norb*nspin
                cs = chi0q.shape
                assert cs[1] == self.lattice.nvol, "lattice volume"
                nd = cs[2]
                assert nd == self.nd or nd == self.norb, "shape[2]"
                assert cs[3] == nd, "shape[3]"

                if nd == self.nd:
                    self.spin_mode = "spinful"
                else:
                    self.spin_mode = "spin-free"
            elif len(chi0q.shape) == 5:
                # spin-diagonal
                #   shape = (nblock,nmat,nvol,norb,norb)
                cs = chi0q.shape
                assert cs[0] == 2, "spin block"
                assert cs[2] == self.lattice.nvol, "lattice volume"
                nd = cs[3]
                assert nd == self.norb, "shape[3]"
                assert cs[4] == nd, "shape[4]"

                self.spin_mode = "spin-diag"
            else:
                assert False, "unexpected shape for reduced scheme: {}".format(chi0q.shape)
        else:
            assert False, "unknown scheme: {}".format(self.calc_scheme)

        logger.info("read_chi0q: shape={}, spin_mode={}".format(chi0q.shape, self.spin_mode))

        return chi0q

    def _read_trans_mod(self, file_name):
        logger.debug(">>> RPA._read_trans_mod")

        try:
            logger.info("read trans_mod from {}".format(file_name))
            data = np.load(file_name)
            tab_r = data["trans_mod"]
            logger.debug("read_trans_mod: shape={}".format(tab_r.shape))
        except Exception as e:
            logger.error("read_trans_mod failed: {}".format(e))
            sys.exit(1)

        if self.lattice.has_sublattice:
            # use reshape green to convert layout
            tab_r = self.lattice._reshape_green(tab_r)

        nx,ny,nz = self.lattice.shape
        nvol = self.lattice.nvol
        nd = self.nd

        tab_k = np.fft.fftn(tab_r.reshape(nx,ny,nz,nd,nd), axes=(0,1,2)).reshape(nvol,nd,nd)

        return tab_k

    def _read_green(self, file_name):
        logger.debug(">>> RPA._read_green")
        try:
            logger.info("read green from {}".format(file_name))
            data = np.load(file_name)
            green = data["green"]
            logger.debug("read_green: shape={}".format(green.shape))
        except Exception as e:
            logger.error("read_green failed: {}".format(e))
            sys.exit(1)

        nvol = self.lattice.nvol
        nd = self.nd
        return green.reshape(nvol,nd,nd)

    def _reshape_green(self, green_):
        # convert green function into sublattice

        Lx,Ly,Lz = self.lattice.cellshape
        Lvol = Lx * Ly * Lz
        Bx,By,Bz = self.lattice.subshape
        Bvol = Bx * By * Bz
        Nx,Ny,Nz = self.lattice.shape
        Nvol = Nx * Ny * Nz

        norb_orig = self.norb_orig
        norb = self.norb
        ns = self.ns

        # check array size
        #assert(green.shape == (Lvol,ns,norb_orig,ns,norb_orig))
        green = green_.reshape(Lvol,ns,norb_orig,ns,norb_orig)

        def _pack_index(x, n):
            _ix, _iy, _iz = x
            _nx, _ny, _nz = n
            return _ix + _nx * (_iy + _ny * (_iz))

        def _unpack_index(x, n):
            _nx, _ny, _nz = n
            _ix = x % _nx
            _iy = (x // _nx) % _ny
            _iz = (x // (_nx * _ny)) % _nz
            return (_ix, _iy, _iz)

        green_sub = np.zeros((Nvol,ns,norb,ns,norb), dtype=np.complex128)

        for isite in range(Nvol):
            ixx, iyy, izz = _unpack_index(isite, (Nx,Ny,Nz))
            ix0, iy0, iz0 = ixx * Bx, iyy * By, izz * Bz

            for aa, bb in itertools.product(range(norb), range(norb)):
                a, ri = aa % norb_orig, aa // norb_orig
                b, rj = bb % norb_orig, bb // norb_orig

                rix, riy, riz = _unpack_index(ri, (Bx,By,Bz))
                rjx, rjy, rjz = _unpack_index(rj, (Bx,By,Bz))

                ix = (ix0 + rjx - rix) % Lx
                iy = (iy0 + rjy - riy) % Ly
                iz = (iz0 + rjz - riz) % Lz

                jsite = _pack_index((ix,iy,iz), (Lx,Ly,Lz))

                for s, t in itertools.product(range(ns), range(ns)):
                    green_sub[isite, s, aa, t, bb] = green[jsite, s, a, t, b]

        return green_sub

    def _calc_trans_mod(self, g0):
        logger.debug(">>> RPA._calc_trans_mod")

        nx,ny,nz = self.lattice.shape
        nvol = self.lattice.nvol
        nd = self.nd
        norb = self.norb

        gg = g0[0]
        ww = self.ham_info.ham_inter_r.reshape(nvol,nd,nd,nd,nd)

        hh1 = np.einsum('rbacd,cd->rab', ww, gg)
        hh2 = np.einsum('rcdab,dc->rab', ww, gg)
        hh3 = np.sum(hh1+hh2, axis=0)/2

        if self.ham_info.enable_spin_orbital:
            H0r = self.ham_info.ham_trans_r.reshape(nvol,nd,nd)
        else:
            H0r = np.einsum('kab,st->ksatb',
                            self.ham_info.ham_trans_r.reshape(nvol,norb,norb),
                            np.eye(2)).reshape(nvol,nd,nd)
        H0r[0] += hh3

        H0 = np.fft.fftn(H0r.reshape(nx,ny,nz,nd,nd), axes=(0,1,2)).reshape(nvol,nd,nd)
        return H0

    @do_profile
    def _calc_epsilon_k(self, green_info):
        logger.debug(">>> RPA._calc_epsilon_k")

        nvol = self.lattice.nvol

        # Find transfer term H0(k)
        if "trans_mod" in green_info:
            logger.debug("calc_epsilon_k: use trans_mod")
            H0 = green_info['trans_mod']
            do_spin_orbital = True

        elif "green_init" in green_info:
            logger.debug("calc_epsilon_k: use initial green")
            H0 = self._calc_trans_mod(green_info["green_init"])
            do_spin_orbital = True

        else:
            H0 = self.ham_info.ham_trans_q
            do_spin_orbital = self.ham_info.enable_spin_orbital

        if do_spin_orbital:
            # H0(k) = T_{a~b~}(k) + h sigma_z_{ss'} H_{ab}(k)
            #   T_{a~b~} with extended orbital index
            #   H_{ab} with bare orbital index
            #   sigma_z = diag(1,-1) Pauli matrix
            #   h coefficient

            if self.ham_info.ham_extern_q is not None:
                H1 = self.ham_info.ham_extern_q * self.ext
                Sz = np.diag((1,-1))
                H0 += np.einsum('kab,st->ksatb', H1, Sz).reshape(H0.shape)

            # check if diagonal
            ns = self.ns
            norb = self.norb

            Htmp = H0.reshape(nvol,ns,norb,ns,norb)
            if np.allclose(Htmp[:,0,:,1,:], 0) and np.allclose(Htmp[:,1,:,0,:], 0):
                if np.allclose(Htmp[:,0,:,0,:], Htmp[:,1,:,1,:]):
                    logger.info("H is spin-free")
                    H0 = Htmp[:,0,:,0,:].reshape(1,nvol,norb,norb)
                    self.spin_mode = "spin-free"
                else:
                    logger.info("H is spin-diagnoal")
                    Hnew = np.zeros((2,nvol,norb,norb), dtype=np.complex128)
                    Hnew[0] = Htmp[:,0,:,0,:]
                    Hnew[1] = Htmp[:,1,:,1,:]
                    H0 = Hnew
                    self.spin_mode = "spin-diag"
            else:
                logger.debug("H is spinful")
                H0 = H0.reshape(1,*H0.shape)
                self.spin_mode = "spinful"

        else:
            # H0(k) = T_{ab}(k) x 1_{ss'} + h Sz_{ss'} H_{ab}(k)
            #   T_{a~b~} with bare orbital index

            if self.ham_info.ham_extern_q is not None:
                H1 = self.ham_info.ham_extern_q * self.ext

                ns = self.ns
                norb = self.norb

                Hnew = np.zeros((2,nvol,norb,norb), dtype=np.complex128)
                Hnew[0] = H0 + H1
                Hnew[1] = H0 - H1
                H0 = Hnew
                self.spin_mode = "spin-diag"
                logger.info("H is spin-diagnoal")
            else:
                logger.debug("H is spin-free")
                H0 = H0.reshape(1,*H0.shape)
                self.spin_mode = "spin-free"

        # diagonalize H0(k)
        w,v = np.linalg.eigh(H0)

        self.H0_eigenvalue = w
        self.H0_eigenvector = v

    @do_profile
    def _find_mu(self, Ncond, T):
        logger.debug(">>> RPA._find_mu")

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
            nn = np.einsum('gkal,gkl,gkal->', np.conjugate(v), ff, v)
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
        ew: eigenvalues  ew[g,k,i]    i-th eigenvalue of wavenumber k in block g
        ev: eigenvectors ev[g,k,a,i]  i-th eigenvector corresponding to ew[g,k,i]
        beta: inverse temperature
        mu: chemical potential
        """
        logger.debug(">>> RPA._calc_green")

        # load eigenvalues and eigenvectors
        ew = self.H0_eigenvalue
        ev = self.H0_eigenvector

        nx,ny,nz = self.lattice.shape
        #nvol = self.lattice.nvol

        nblock,nvol,nd = ew.shape
        assert nvol == self.lattice.nvol

        nmat = self.nmat

        iomega = (np.arange(nmat) * 2 + 1 - nmat) * np.pi / beta

        # 1 / (iw_{n} - (e_i(k) - mu)) -> g[n,k,i]
        wn = np.transpose(np.tile(1j * iomega, (nblock,nvol,nd,1)), axes=(0,3,1,2))
        ek = np.transpose(np.tile((ew - mu), (nmat,1,1,1)), axes=(1,0,2,3))

        # tail improvement
        aa = self.coeff_tail
        g = 1.0 / (wn - ek) - aa / wn

        # G_ab(k,iw_n) = sum_j d_{a,j} d_{b,j}^* / (iw_{n} - (e_j(k) - mu))
        green = np.einsum('gkaj,gkbj,glkj->glkab', ev, np.conj(ev), g)
        green_tail = np.einsum('gkaj,gkbj,gl->glkab', ev, np.conj(ev), np.tile(np.ones(nmat), (nblock,1))) * aa * 0.5 * beta

        return green, green_tail

    @do_profile
    def _calc_chi0q(self, green_kw, green0_tail, beta):
        # green function
        #   green_kw[g,l,k,a,b]: G_ab(k,ie) in block g
        #     l : matsubara freq i\epsilon_l = i\pi(2*l+1-N)/\beta for l=0..N-1
        #     k : wave number kx,ky,kz
        #     a,b : orbital and spin index a=(s,alpha) b=(t,beta)

        logger.debug(">>> RPA._calc_chi0q")

        nx,ny,nz = self.lattice.shape
        #nvol = self.lattice.nvol

        nblock,nmat,nvol,nd,nd2 = green_kw.shape
        assert nvol == self.lattice.nvol
        assert nmat == self.nmat

        # Fourier transform from matsubara freq to imaginary time
        omg = np.exp(-1j * np.pi * (1.0/nmat - 1.0) * np.arange(nmat))

        green_kt = np.einsum('gtv,t->gtv',
                             FFT.fft(green_kw.reshape(nblock,nmat,nvol*nd*nd), axis=1),
                             omg).reshape(nblock,nmat,nx,ny,nz,nd,nd)
        green_kt -= green0_tail.reshape(nblock,nmat,nx,ny,nz,nd,nd)

        # Fourier transform from wave number space to coordinate space
        green_rt = FFT.ifftn(green_kt.reshape(nblock,nmat,nx,ny,nz,nd*nd), axes=(2,3,4))

        # calculate chi0(r,t)[a,a',b,b'] = G(r,t)[a,b] * G(-r,-t)[b',a']
        green_rev = np.flip(np.roll(green_rt, -1, axis=(1,2,3,4)), axis=(1,2,3,4)).reshape(nblock,nmat,nvol,nd,nd)

        sgn = np.full(nmat, -1)
        sgn[0] = 1

        if self.enable_reduced:
            # reduced index calculation
            chi0_rt = np.einsum('glrab,glrba,l->glrab',
                                green_rt.reshape(nblock,nmat,nvol,nd,nd),
                                green_rev,
                                sgn)
            nd_shape=(nd,nd)
            nds=nd**2
        else:
            chi0_rt = np.einsum('glrab,glrdc,l->glracbd',
                                green_rt.reshape(nblock,nmat,nvol,nd,nd),
                                green_rev,
                                sgn)
            nd_shape=(nd,nd,nd,nd)
            nds=nd**4

        # Fourier transform to wave number space
        chi0_qt = FFT.fftn(chi0_rt.reshape(nblock,nmat,nx,ny,nz,nds), axes=(2,3,4))

        # Fourier transform to matsubara freq
        omg = np.exp(1j * np.pi * (-1) * np.arange(nmat))

        chi0_qw = FFT.ifft(
            np.einsum('gtv,t->gtv', chi0_qt.reshape(nblock,nmat,nvol*nds), omg),
            axis=1).reshape(nblock,nmat,nvol,*nd_shape) * (-1.0/beta)

        return chi0_qw

    @do_profile
    def _solve_rpa(self, chi0q, ham):
        logger.debug(">>> RPA._solve_rpa")

        nvol = self.lattice.nvol
        #nmat = self.nmat
        nmat = chi0q.shape[0]
        #nd = self.nd
        #ndx = nd**2  # combined index a = (alpha, alpha')
        chi_shape = chi0q.shape  # [nmat,nvol,(spin_orbital structure)]
        ndx = np.prod(chi_shape[2:2+(len(chi_shape)-2)//2])

        # 1 + X^0(l,k,aa,bb) W(k,bb,cc)
        mat  = np.tile(np.eye(ndx, dtype=np.complex128), (nmat,nvol,1,1))
        mat += np.einsum('lkab,kbc->lkac',
                         chi0q.reshape(nmat,nvol,ndx,ndx),
                         ham.reshape(nvol,ndx,ndx))

        # [ 1 + X^0 W ]^-1 X^0
        sol = np.linalg.solve(mat, chi0q.reshape(nmat,nvol,ndx,ndx))

        return sol.reshape(chi_shape)


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
        green_info.update( solver.read_init(info_inputfile) )

        logger.info("Start RPA calculation")
        solver.solve(green_info, path_to_output)

        logger.info("Save calculation results")
        solver.save_results(info_outputfile, green_info)

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
