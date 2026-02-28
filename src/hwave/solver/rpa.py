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

        # calc_type: "ring" (default) or "ring+ladder"
        self.calc_type = info_mode.get("calc_type", "ring")
        if self.calc_type not in ["ring", "ring+ladder"]:
            logger.error("calc_type must be 'ring' or 'ring+ladder', got '{}'".format(self.calc_type))
            sys.exit(1)

        # auto choose
        if self.calc_scheme == "auto":
            if not self.ham_info.has_interaction():
                logger.error("calc_scheme must be specified for chi0q-only mode.")
                sys.exit(1)
            else:
                if self.calc_type == "ring+ladder":
                    # ladder diagrams require general scheme (full rank-4 tensor)
                    self.calc_scheme = "general"
                    logger.info("auto mode for calc_scheme: set to general (ring+ladder)")
                elif self.ham_info.has_interaction_exchange():
                    self.calc_scheme = "squashed"
                    logger.info("auto mode for calc_scheme: set to squashed")
                else:
                    self.calc_scheme = "reduced"
                    logger.info("auto mode for calc_scheme: set to reduced")

        # consistency check
        if self.calc_scheme == "reduced" and self.ham_info.has_interaction_exchange():
            logger.error("calc_scheme=reduced is not compatible with exchange-type interaction.")
            sys.exit(1)
        if self.calc_type == "ring+ladder" and self.calc_scheme != "general":
            logger.error("calc_type='ring+ladder' requires calc_scheme='general' or 'auto'.")
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
        elif have_Ncond and have_filling:
            # both Ncond or filling
            logger.error("both Ncond and filling are specified")
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
        logger.info("    calc_type       = {}".format(self.calc_type))
        pass

    @do_profile
    def solve(self, green_info, path_to_output):
        """Solve the RPA equation to calculate susceptibility.

        This is the main method that performs RPA calculations. It either calculates 
        or loads chi0q, transforms interaction Hamiltonians based on spin state, 
        and solves the RPA equation.

        Parameters
        ----------
        green_info : dict
            Dictionary containing Green's function information and calculation parameters.
            May include pre-calculated chi0q.
        path_to_output : str
            Path to directory where output files will be saved.

        Notes
        -----
        The calculation flow is:
        1. Calculate/load chi0q
        2. Transform interactions based on spin state (spin-free/spinful/spin-diag)
        3. Transform tensors based on calculation scheme (reduced/squashed/general)
        4. Solve RPA equation to get chiq
        """
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

                if self.calc_scheme == "reduced" or self.calc_scheme == "squashed":
                    # Treat combined spin-orbital indices as general orbitals.
                    # squashed degenerates to reduced; block structure is
                    # exploited by _find_block_diagonal inside _solve_rpa.
                    nvol = self.lattice.nvol
                    nd = self.nd
                    ham = np.einsum('kaabb->kab',
                                    ham_orig.reshape(nvol,*(nd,)*4)).reshape(nvol,*(nd,)*2)
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

            # solve longitudinal (ring) RPA
            sol = self._solve_rpa(chi0q, ham)

            # adhoc store
            green_info["chiq"] = sol

            # Solve transverse (ladder) RPA if requested
            if self.calc_type == "ring+ladder":
                chi0q_pm, ham_pm = self._build_transverse_channel(
                    chi0q_orig, ham_orig)
                sol_pm = self._solve_rpa(chi0q_pm, ham_pm)
                green_info["chiq_pm"] = sol_pm

        logger.info("End RPA calculations")
        pass

    @do_profile
    def save_results(self, info_outputfile, green_info):
        """Save calculation results to files.

        Parameters
        ----------
        info_outputfile : dict
            Dictionary containing output file configuration.
        green_info : dict
            Dictionary containing calculation information and results.

        Notes
        -----
        Saves the following information:
        - Calculated correlation functions (chi0q/chiq)
        - Matsubara frequency indices
        - Wave vector unit vectors
        - Wave vector indices
        """
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

        kx = _klist(nx)
        ky = _klist(ny)
        kz = _klist(nz)
        # Build (nx,ny,nz) grids for each k-component
        kx_g, ky_g, kz_g = np.meshgrid(kx, ky, kz, indexing='ij')
        # wtable[ix,iy,iz,:] = kvec[0]*kx + kvec[1]*ky + kvec[2]*kz
        wtable = (kx_g[..., np.newaxis] * kvec[0]
                + ky_g[..., np.newaxis] * kvec[1]
                + kz_g[..., np.newaxis] * kvec[2])
        self.wave_table = wtable.reshape(nvol, 3)

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
        """Read chi0q from file and validate its shape.

        Parameters
        ----------
        file_name : str
            Path to the file containing chi0q data.

        Returns
        -------
        ndarray
            The loaded chi0q array.

        Raises
        ------
        AssertionError
            If the loaded data doesn't match expected dimensions or format.

        Notes
        -----
        Performs checks for:
        - Lattice volume consistency
        - Number of orbitals consistency 
        - Spin freedom consistency
        - Calculation scheme compatibility
        """
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
        """Read transfer integral modifications from file.

        Parameters
        ----------
        file_name : str
            Path to the file containing transfer integral modifications.

        Returns
        -------
        ndarray
            The transfer integrals in k-space.

        Raises
        ------
        RuntimeError
            If file reading fails.

        Notes
        -----
        - Converts layout if sublattice exists
        - Performs Fourier transform to k-space
        """
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

        # Build index mapping tables (vectorized)
        # Supercell site indices
        isite_arr = np.arange(Nvol)
        ixx = isite_arr % Nx
        iyy = (isite_arr // Nx) % Ny
        izz = (isite_arr // (Nx * Ny)) % Nz
        ix0 = ixx * Bx  # (Nvol,)
        iy0 = iyy * By
        iz0 = izz * Bz

        # Orbital decomposition
        orb_arr = np.arange(norb)
        a_arr = orb_arr % norb_orig          # original orbital index
        ri_arr = orb_arr // norb_orig         # sublattice index
        rix = ri_arr % Bx
        riy = (ri_arr // Bx) % By
        riz = (ri_arr // (Bx * By)) % Bz

        # Compute jsite for all (isite, aa, bb) combinations
        # drx[aa,bb] = rjx[bb] - rix[aa], etc.
        drx = rix[np.newaxis, :] - rix[:, np.newaxis]  # (norb, norb)
        dry = riy[np.newaxis, :] - riy[:, np.newaxis]
        drz = riz[np.newaxis, :] - riz[:, np.newaxis]

        # ix[isite, aa, bb] = (ix0[isite] + drx[aa, bb]) % Lx
        jx = (ix0[:, np.newaxis, np.newaxis] + drx[np.newaxis, :, :]) % Lx  # (Nvol, norb, norb)
        jy = (iy0[:, np.newaxis, np.newaxis] + dry[np.newaxis, :, :]) % Ly
        jz = (iz0[:, np.newaxis, np.newaxis] + drz[np.newaxis, :, :]) % Lz
        jsite_map = jx + Lx * (jy + Ly * jz)  # (Nvol, norb, norb)

        # Source orbital indices: a_arr[aa], a_arr[bb]
        a_src = a_arr  # (norb,) - maps aa -> original orbital a
        b_src = a_arr  # (norb,) - maps bb -> original orbital b

        # Gather using advanced indexing
        # green[jsite, s, a, t, b] -> green_sub[isite, s, aa, t, bb]
        green_sub = green[jsite_map][:, :, :, :, a_src, :, :][:, :, :, :, :, :, b_src]
        # Shape: (Nvol, norb, norb, ns, norb, ns, norb) - need to select diagonal
        # Use explicit indexing for clarity
        green_sub = green[
            jsite_map[:, :, :, np.newaxis, np.newaxis],   # (Nvol, norb, norb, 1, 1)
            np.arange(ns)[np.newaxis, np.newaxis, np.newaxis, :, np.newaxis],  # s
            a_src[np.newaxis, :, np.newaxis, np.newaxis, np.newaxis],          # a
            np.arange(ns)[np.newaxis, np.newaxis, np.newaxis, np.newaxis, :],  # t
            b_src[np.newaxis, np.newaxis, :, np.newaxis, np.newaxis],          # b
        ]  # shape: (Nvol, norb, norb, ns, ns)

        # Transpose to match (Nvol, ns, norb, ns, norb)
        green_sub = green_sub.transpose(0, 3, 1, 4, 2)

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

        # diagonalize H0(k) with optional block decomposition
        nblock_spin = H0.shape[0]
        nd_block = H0.shape[-1]
        blocks = self._find_block_diagonal(H0.reshape(nblock_spin * nvol, nd_block, nd_block))

        if blocks is not None and len(blocks) > 1:
            logger.info("_calc_epsilon_k: orbital block structure detected, "
                        "nd={} -> {} blocks of sizes {}".format(
                            nd_block, len(blocks), [len(b) for b in blocks]))
            w = np.zeros((nblock_spin, nvol, nd_block), dtype=np.float64)
            v = np.zeros((nblock_spin, nvol, nd_block, nd_block), dtype=np.complex128)
            col_offset = 0
            for blk_idx in blocks:
                idx = np.array(blk_idx)
                ix = np.ix_(idx, idx)
                H0_blk = H0[:, :, ix[0], ix[1]]
                wb, vb = np.linalg.eigh(H0_blk)
                nb = len(idx)
                w[:, :, col_offset:col_offset + nb] = wb
                v[np.ix_(np.arange(nblock_spin), np.arange(nvol), idx,
                  np.arange(col_offset, col_offset + nb))] = vb
                col_offset += nb
        else:
            w,v = np.linalg.eigh(H0)

        self.H0_eigenvalue = w
        self.H0_eigenvector = v

    @do_profile
    def _find_mu(self, Ncond, T):
        logger.debug(">>> RPA._find_mu")

        from scipy import optimize

        # load eigenvalues (eigenvectors not needed thanks to unitarity)
        w = self.H0_eigenvalue
        # fetch parameters
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

        # Exploit unitarity of eigenvectors:
        # Tr[V† diag(f) V] = sum_l f(ε_l)
        # This eliminates the O(nd²) einsum per iteration.
        def _calc_delta_n(mu):
            ff = _fermi(T, mu, w)
            return np.sum(ff).real - occupied_number

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

        nblock,nvol,nd = ew.shape
        assert nvol == self.lattice.nvol

        nmat = self.nmat

        iomega = (np.arange(nmat) * 2 + 1 - nmat) * np.pi / beta

        # 1 / (iw_{n} - (e_i(k) - mu)) -> g[g,l,k,i] via broadcasting
        # iomega: (nmat,), ew: (nblock, nvol, nd)
        wn = 1j * iomega[np.newaxis, :, np.newaxis, np.newaxis]  # (1,nmat,1,1)
        ek = (ew - mu)[:, np.newaxis, :, :]  # (nblock,1,nvol,nd)

        # tail improvement
        aa = self.coeff_tail
        g = 1.0 / (wn - ek) - aa / wn  # (nblock,nmat,nvol,nd)

        # G_ab(k,iw_n) = sum_j V_{a,j} V*_{b,j} * g_j
        # = (V * g) @ V†  -- use matmul (BLAS) instead of einsum
        ev_conj_t = np.conj(ev).swapaxes(-2, -1)  # (nblock,nvol,nd,nd): V†[g,k]

        # Vg = V * g: broadcast g into eigenvector columns
        Vg = ev[:, np.newaxis, :, :, :] * g[:, :, :, np.newaxis, :]  # (nblock,nmat,nvol,nd,nd)
        green = Vg @ ev_conj_t[:, np.newaxis, :, :, :]  # (nblock,nmat,nvol,nd,nd)

        # Tail: G_tail = V @ V† * aa * 0.5 * beta = I * aa * 0.5 * beta (unitarity)
        # But original code retains V V† form for non-complete basis cases
        VVt = ev @ ev_conj_t  # (nblock,nvol,nd,nd)
        green_tail = VVt[:, np.newaxis, :, :, :] * np.ones((1, nmat, 1, 1, 1)) * aa * 0.5 * beta

        return green, green_tail

    @do_profile
    def _calc_chi0q(self, green_kw, green0_tail, beta):
        """Calculate the bare susceptibility chi0q.

        Parameters
        ----------
        green_kw : ndarray
            Green's function in k-space and Matsubara frequency.
            Shape: [g,l,k,a,b] where:
            - g: block index
            - l: Matsubara frequency index
            - k: wave number index
            - a,b: orbital and spin indices
        green0_tail : ndarray
            High-frequency tail correction for Green's function.
        beta : float
            Inverse temperature.

        Returns
        -------
        ndarray
            The calculated chi0q.

        Notes
        -----
        Calculation steps:
        1. Fourier transform from Matsubara freq to imaginary time
        2. Transform from k-space to real space
        3. Calculate chi0 in real space and imaginary time
        4. Transform back to k-space and Matsubara frequency
        """
        logger.debug(">>> RPA._calc_chi0q")

        nx,ny,nz = self.lattice.shape
        #nvol = self.lattice.nvol

        nblock,nmat,nvol,nd,nd2 = green_kw.shape
        assert nvol == self.lattice.nvol
        assert nmat == self.nmat

        # Fourier transform from Matsubara freq to imaginary time
        omg = np.exp(-1j * np.pi * (1.0/nmat - 1.0) * np.arange(nmat))

        # In-place multiply to avoid extra copy
        green_flat = green_kw.reshape(nblock, nmat, nvol * nd * nd)
        green_kt = FFT.fft(green_flat, axis=1)
        green_kt *= omg[np.newaxis, :, np.newaxis]
        green_kt = green_kt.reshape(nblock, nmat, nx, ny, nz, nd, nd)
        green_kt -= green0_tail.reshape(nblock, nmat, nx, ny, nz, nd, nd)

        # Fourier transform from wave number space to coordinate space
        green_rt = FFT.ifftn(green_kt.reshape(nblock, nmat, nx, ny, nz, nd * nd), axes=(2, 3, 4))

        # calculate chi0 in real space and imaginary time
        green_rev = np.flip(np.roll(green_rt, -1, axis=(1, 2, 3, 4)), axis=(1, 2, 3, 4)).reshape(nblock, nmat, nvol, nd, nd)

        sgn = np.full(nmat, -1)
        sgn[0] = 1

        if self.enable_reduced:
            # reduced index calculation
            # chi0[g,l,r,a,b] = G[g,l,r,a,b] * G_rev[g,l,r,b,a] * sgn[l]
            green_rt_5d = green_rt.reshape(nblock, nmat, nvol, nd, nd)
            chi0_rt = (green_rt_5d
                       * green_rev.swapaxes(-2, -1)
                       * sgn[np.newaxis, :, np.newaxis, np.newaxis, np.newaxis])
            nd_shape = (nd, nd)
            nds = nd ** 2
        else:
            # General index: chi0[g,l,r,a,c,b,d] = G[g,l,r,a,b] * G_rev[g,l,r,d,c] * sgn[l]
            # Use outer product via broadcasting instead of einsum
            G_fwd = green_rt.reshape(nblock, nmat, nvol, nd, nd)
            sgn_bc = sgn[np.newaxis, :, np.newaxis, np.newaxis, np.newaxis]
            # G_fwd[:,:,:,a,b] * sgn -> shape (g,l,r,a,b)
            # G_rev[:,:,:,d,c]        -> shape (g,l,r,d,c)
            # result[:,:,:,a,c,b,d] = G_fwd[:,:,:,a,b] * G_rev[:,:,:,d,c] * sgn
            chi0_rt = ((G_fwd * sgn_bc)[:, :, :, :, np.newaxis, :, np.newaxis]
                       * green_rev[:, :, :, np.newaxis, :, np.newaxis, :])
            # shape: (g,l,r,a,d,b,c) -> need (g,l,r,a,c,b,d)
            chi0_rt = chi0_rt.transpose(0, 1, 2, 3, 6, 5, 4)
            nd_shape = (nd, nd, nd, nd)
            nds = nd ** 4

        # Fourier transform to wave number space
        chi0_qt = FFT.fftn(chi0_rt.reshape(nblock, nmat, nx, ny, nz, nds), axes=(2, 3, 4))

        # Fourier transform to matsubara freq
        omg2 = np.exp(1j * np.pi * (-1) * np.arange(nmat))

        chi0_qt_flat = chi0_qt.reshape(nblock, nmat, nvol * nds)
        chi0_qt_flat *= omg2[np.newaxis, :, np.newaxis]
        chi0_qw = FFT.ifft(chi0_qt_flat, axis=1).reshape(nblock, nmat, nvol, *nd_shape) * (-1.0 / beta)

        return chi0_qw

    def _calc_chi0q_transverse(self, green_kw, green0_tail, beta):
        """Calculate the transverse bare susceptibility chi0_+-(q,iω).

        chi0_+-[a,c,b,d](r,τ) = -G_↑[a,b](r,τ) * G_↓[d,c](-r,-τ)

        This crosses spin-up and spin-down Green's functions, unlike the
        longitudinal chi0 which uses same-spin products.

        Parameters
        ----------
        green_kw : ndarray, shape (2, nmat, nvol, norb, norb)
            Green's function with block 0=↑, block 1=↓.
        green0_tail : ndarray
            High-frequency tail correction.
        beta : float
            Inverse temperature.

        Returns
        -------
        ndarray
            Transverse chi0_+- with block dimension removed (shape depends
            on enable_reduced).
        """
        logger.debug(">>> RPA._calc_chi0q_transverse")

        nx, ny, nz = self.lattice.shape
        nblock, nmat, nvol, nd, nd2 = green_kw.shape
        assert nblock == 2, "Transverse chi0 requires spin-diag (nblock=2)"

        # Fourier transform from Matsubara freq to imaginary time
        omg = np.exp(-1j * np.pi * (1.0/nmat - 1.0) * np.arange(nmat))

        green_kt = (FFT.fft(green_kw.reshape(nblock, nmat, nvol*nd*nd), axis=1)
                    * omg[np.newaxis, :, np.newaxis]
                    ).reshape(nblock, nmat, nx, ny, nz, nd, nd)
        green_kt -= green0_tail.reshape(nblock, nmat, nx, ny, nz, nd, nd)

        # Fourier transform from k-space to real space
        green_rt = FFT.ifftn(green_kt.reshape(nblock, nmat, nx, ny, nz, nd*nd),
                             axes=(2, 3, 4)).reshape(nblock, nmat, nvol, nd, nd)

        # G_↓(-r,-τ): flip r and τ, then shift
        green_dn_rev = np.flip(np.roll(green_rt[1:2], -1, axis=(1, 2, 3, 4)),
                               axis=(1, 2, 3, 4)).reshape(nmat, nvol, nd, nd)

        # G_↑(r,τ)
        green_up_rt = green_rt[0].reshape(nmat, nvol, nd, nd)

        sgn = np.full(nmat, -1)
        sgn[0] = 1

        if self.enable_reduced:
            # chi0_+-[l,r,a,b] = G_↑[l,r,a,b] * G_↓_rev[l,r,b,a] * sgn[l]
            # (same contraction as longitudinal but crossing spin blocks)
            chi0_rt = (green_up_rt
                       * green_dn_rev.swapaxes(-2, -1)
                       * sgn[:, np.newaxis, np.newaxis, np.newaxis])
            nd_shape = (nd, nd)
            nds = nd**2
        else:
            # chi0_+-[l,r,a,c,b,d] = G_↑[l,r,a,b] * G_↓_rev[l,r,d,c] * sgn[l]
            chi0_rt = np.einsum('lrab,lrdc,l->lracbd',
                                green_up_rt, green_dn_rev, sgn)
            nd_shape = (nd, nd, nd, nd)
            nds = nd**4

        # Fourier transform to k-space
        chi0_qt = FFT.fftn(chi0_rt.reshape(nmat, nx, ny, nz, nds), axes=(1, 2, 3))

        # Fourier transform to Matsubara frequency
        omg = np.exp(1j * np.pi * (-1) * np.arange(nmat))
        chi0_qw = FFT.ifft(
            chi0_qt.reshape(nmat, nvol*nds) * omg[:, np.newaxis],
            axis=0).reshape(nmat, nvol, *nd_shape) * (-1.0/beta)

        return chi0_qw

    def _build_transverse_channel(self, chi0q_orig, ham_orig):
        """Build transverse (ladder) channel for RPA.

        The transverse susceptibility chi_+-(q) describes spin-flip
        correlations <S^+(q) S^-(-q)>.

        For paramagnetic systems, the transverse bare susceptibility
        has the same numerical values as the longitudinal chi0:
            chi0_+-[a,c,b,d] = -G_↑[a,b] * G_↓[d,c] = chi0_orb[a,c,b,d]

        The transverse vertex W_+- is obtained by crossing the Hartree
        vertex (Fock exchange):
            W_+-[a,c,b,d] = Gamma_H[(↑,d),(↓,c),(↑,a),(↓,b)]

        Parameters
        ----------
        chi0q_orig : ndarray
            Original bare susceptibility (before spin inflation).
        ham_orig : ndarray
            Original interaction Hamiltonian in spin-orbital space.

        Returns
        -------
        chi0q_pm : ndarray
            Transverse bare susceptibility, shape matches general scheme.
        ham_pm : ndarray
            Transverse vertex, shape (nvol, norb, norb, norb, norb).
        """
        norb = self.norb
        ns = self.ns
        nd = norb * ns
        nvol = self.lattice.nvol

        # --- Build chi0_+- ---
        # For paramagnetic and spin-diagonal cases, chi0_+- = chi0_orb
        # (same numerical values, just interpreted as transverse bubble)
        if self.spin_mode == "spin-free":
            # chi0q_orig shape: (nfreq, nvol, norb, norb, norb, norb) for general
            # or (nfreq, nvol, norb, norb) for reduced
            if chi0q_orig.ndim == 4:
                # reduced: expand to general
                # chi0q_pm[:, :, l1, l2, l3, l2] = chi0q_orig[:, :, l1, l3]
                # i.e., delta_{l2,l4} structure
                nfreq, nvol_c, n1, n2 = chi0q_orig.shape
                chi0q_pm = np.zeros((nfreq, nvol_c, norb, norb, norb, norb),
                                    dtype=np.complex128)
                # Vectorized: broadcast chi0q_orig into diagonal l2=l4 positions
                for l2 in range(norb):
                    chi0q_pm[:, :, :, l2, :, l2] = chi0q_orig
            else:
                chi0q_pm = chi0q_orig.copy()

        elif self.spin_mode == "spin-diag":
            # Compute exact chi0_+- from G_↑ and G_↓ Green's functions.
            # chi0_+-[a,c,b,d](r,τ) = -G_↑[a,b](r,τ) * G_↓[d,c](-r,-τ)
            if hasattr(self, 'green0') and self.green0 is not None:
                chi0q_pm_full = self._calc_chi0q_transverse(
                    self.green0, self.green0_tail, 1.0 / self.T)
                # Filter by freq_index if needed
                if len(self.freq_index) < self.nmat:
                    chi0q_pm = chi0q_pm_full[self.freq_index]
                else:
                    chi0q_pm = chi0q_pm_full
                # Expand reduced to general if needed for vertex contraction
                if self.enable_reduced and chi0q_pm.ndim == 4:
                    nfreq, nvol_c, n1, n2 = chi0q_pm.shape
                    chi0q_pm_gen = np.zeros(
                        (nfreq, nvol_c, norb, norb, norb, norb),
                        dtype=np.complex128)
                    for l2 in range(norb):
                        chi0q_pm_gen[:, :, :, l2, :, l2] = chi0q_pm
                    chi0q_pm = chi0q_pm_gen
            else:
                # Fallback: Green's functions not stored (e.g. chi0q provided externally)
                logger.warning(
                    "spin-diag transverse channel: Green's functions not available, "
                    "falling back to spin-↑ block approximation for chi0_+-")
                if chi0q_orig.ndim == 5:
                    chi0q_pm = chi0q_orig[0].copy()
                else:
                    chi0q_pm = chi0q_orig[0].copy()

        elif self.spin_mode == "spinful":
            # Already in spin-orbital space
            # chi0_+- requires re-extraction with proper spin indices
            # For now, extract from the chi0_SO structure
            # chi0_SO[a,c,b,d] -> chi0_+-[a_orb,c_orb,b_orb,d_orb]
            #   = chi0_SO[(↑,a),(↓,c),(↑,b),(↓,d)]  (when s1=↑=s3, s2=↓=s4)
            #   which is just the orbital diagonal of chi0
            if chi0q_orig.ndim == 6:
                nfreq = chi0q_orig.shape[0]
                chi0q_pm = chi0q_orig[:, :,
                                      0:norb, norb:2*norb,
                                      0:norb, norb:2*norb].copy()
            else:
                nfreq = chi0q_orig.shape[0]
                chi0q_pm = chi0q_orig[:, :, 0:norb, 0:norb].copy()

        # --- Build W_+- (transverse vertex) ---
        # The transverse vertex for the +- (spin-flip) channel:
        #   W_+-[a,c,b,d] = ham[↑a,↑c,↑b,↑d] - ham[↓d,↓b,↑c,↑a]
        #
        # First term: same-spin (↑↑↑↑) block of the Hartree vertex
        # Second term: cross-spin (↓↓↑↑) block (subtracted)
        #
        # Verification for each interaction type:
        #   CoulombIntra U: ham[↑↑↑↑]=0, ham[↓↓↑↑]=U  → W_+-= -U  (correct)
        #   CoulombInter V: ham[↑↑↑↑]=V, ham[↓↓↑↑]=V  → W_+-=  0  (SU(2) correct)
        #   Exchange J:     ham[↑↑↑↑]=0, ham[↓↓↑↑]=0   → W_+-=  0  (SU(2) correct)
        #   PairLift J:     ham[↑↑↑↑]=0, ham[↓↓↑↑]=0   → W_+-=  0  (SU(2) correct)
        #   Hund J:         ham[↑↑↑↑]=-J, ham[↓↓↑↑]=0  → W_+-= -J  (not SU(2) alone)
        #   Ising J:        ham[↑↑↑↑]=J, ham[↓↓↑↑]=-J  → W_+-= 2J  (not SU(2) alone)
        #   Full Kanamori:  W_+- = -(U-2J) = W_zz  (SU(2) correct)
        ham_4d = ham_orig.reshape(nvol, nd, nd, nd, nd)

        # Vectorized extraction using index arrays (replaces norb^4 Python loop)
        # Same-spin block: ham[(↑,a),(↑,c),(↑,b),(↑,d)]
        up_block = ham_4d[:, :norb, :norb, :norb, :norb]  # (nvol, norb, norb, norb, norb)
        # Cross-spin block: ham[(↓,d),(↓,b),(↑,c),(↑,a)] with index reordering
        cross_block = ham_4d[:, norb:, norb:, :norb, :norb]  # (nvol, norb, norb, norb, norb)
        # cross_block[:, d, b, c, a] -> need to transpose to match ham_pm[:, a, c, b, d]
        cross_reordered = cross_block.transpose(0, 4, 3, 2, 1)  # (nvol, a, c, b, d)
        ham_pm = up_block - cross_reordered

        logger.info("ring+ladder: built transverse channel "
                    "(chi0_pm shape={}, ham_pm shape={})".format(
                        chi0q_pm.shape, ham_pm.shape))

        return chi0q_pm, ham_pm

    @do_profile
    def _solve_rpa(self, chi0q, ham):
        """Solve the RPA equation.

        Parameters
        ----------
        chi0q : ndarray
            Bare susceptibility.
        ham : ndarray
            Interaction Hamiltonian.

        Returns
        -------
        ndarray
            The RPA susceptibility chiq.

        Notes
        -----
        Solves the equation:
        chiq = [1 + chi0q * W]^(-1) * chi0q
        where W is the interaction vertex.

        When the matrices have block-diagonal structure (e.g. spin-diagonal case),
        the solver automatically detects and exploits this to reduce problem size.
        Block structure is cached per ham shape+content to avoid re-detection.

        Frequency parallelization: since ham is frequency-independent, each
        frequency slice can be solved independently. When nmat is large enough,
        the solve is distributed across threads using concurrent.futures.
        """
        logger.debug(">>> RPA._solve_rpa")

        nvol = self.lattice.nvol
        nmat = chi0q.shape[0]
        chi_shape = chi0q.shape  # [nmat,nvol,(spin_orbital structure)]
        ndx = np.prod(chi_shape[2:2+(len(chi_shape)-2)//2])

        chi0q_2d = chi0q.reshape(nmat, nvol, ndx, ndx)
        ham_2d = ham.reshape(nvol, ndx, ndx)

        # Use cached block structure if available for same ham
        cache_key = (ham_2d.shape, ham_2d.dtype)
        if not hasattr(self, '_block_cache'):
            self._block_cache = {}

        if cache_key in self._block_cache:
            cached_hash, blocks = self._block_cache[cache_key]
            # Verify by hash (cheap compared to re-detection)
            cur_hash = hash(ham_2d.data.tobytes())
            if cur_hash != cached_hash:
                blocks = self._find_block_diagonal(ham_2d)
                self._block_cache[cache_key] = (cur_hash, blocks)
        else:
            blocks = self._find_block_diagonal(ham_2d)
            self._block_cache[cache_key] = (hash(ham_2d.data.tobytes()), blocks)

        # Determine thread-parallel chunking for frequency axis
        # LAPACK releases the GIL, so threading gives real parallelism.
        # Only parallelize when there is enough work per thread.
        # LAPACK's batched zgesv already uses internal threads for large batches,
        # so explicit threading only helps for very large problems where the
        # per-chunk overhead is negligible relative to compute.
        # Users can force threading via HWAVE_RPA_THREADS env var.
        import os as _os
        n_workers = int(_os.environ.get("HWAVE_RPA_THREADS", "0"))
        if n_workers == 0:
            import multiprocessing
            n_workers = min(multiprocessing.cpu_count(), 4)
        # Heuristic: only parallelize for large multi-orbital problems
        # where ndx >= 16 (8+ orbitals spinful) and enough frequency points.
        # For small ndx, batched LAPACK is faster than thread pool overhead.
        use_parallel = (n_workers > 1 and nmat >= 4 * n_workers
                        and ndx >= 16 and nmat * nvol * ndx * ndx >= 1000000)

        if blocks is not None and len(blocks) > 1:
            logger.info("_solve_rpa: block-diagonal structure detected, "
                        "ndx={} -> {} blocks of sizes {}".format(
                            ndx, len(blocks), [len(b) for b in blocks]))
            sol = np.zeros_like(chi0q_2d)
            for block_idx in blocks:
                idx = np.array(block_idx)
                ix = np.ix_(idx, idx)
                chi0q_blk = chi0q_2d[:, :, ix[0], ix[1]]
                ham_blk = ham_2d[:, ix[0], ix[1]]
                nb = len(idx)

                if use_parallel:
                    sol[:, :, ix[0], ix[1]] = self._solve_rpa_parallel(
                        chi0q_blk, ham_blk, nb, n_workers)
                else:
                    mat_blk = (np.eye(nb, dtype=np.complex128)
                               + (chi0q_blk @ ham_blk[np.newaxis, :, :, :]))
                    sol[:, :, ix[0], ix[1]] = np.linalg.solve(mat_blk, chi0q_blk)
        else:
            if use_parallel:
                sol = self._solve_rpa_parallel(chi0q_2d, ham_2d, ndx, n_workers)
            else:
                mat = (np.eye(ndx, dtype=np.complex128)
                       + (chi0q_2d @ ham_2d[np.newaxis, :, :, :]))
                sol = np.linalg.solve(mat, chi0q_2d)

        return sol.reshape(chi_shape)

    @staticmethod
    def _solve_rpa_parallel(chi0q_2d, ham_2d, ndx, n_workers):
        """Solve RPA equation with frequency-axis thread parallelism.

        Parameters
        ----------
        chi0q_2d : ndarray, shape (nmat, nvol, ndx, ndx)
            Bare susceptibility in 2D matrix form.
        ham_2d : ndarray, shape (nvol, ndx, ndx)
            Frequency-independent interaction Hamiltonian.
        ndx : int
            Matrix dimension.
        n_workers : int
            Number of threads.

        Returns
        -------
        sol : ndarray, shape (nmat, nvol, ndx, ndx)
            RPA susceptibility.
        """
        from concurrent.futures import ThreadPoolExecutor
        nmat = chi0q_2d.shape[0]
        sol = np.empty_like(chi0q_2d)
        eye = np.eye(ndx, dtype=np.complex128)
        ham_bc = ham_2d[np.newaxis, :, :, :]  # broadcast-ready

        # Split frequency axis into contiguous chunks
        chunk_size = (nmat + n_workers - 1) // n_workers
        slices = [slice(i, min(i + chunk_size, nmat))
                  for i in range(0, nmat, chunk_size)]

        def _solve_chunk(sl):
            chi0q_chunk = chi0q_2d[sl]
            mat = eye + (chi0q_chunk @ ham_bc)
            sol[sl] = np.linalg.solve(mat, chi0q_chunk)

        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            list(pool.map(_solve_chunk, slices))

        return sol

    def _find_block_diagonal(self, ham_2d):
        """Detect block-diagonal structure from the interaction Hamiltonian.

        Parameters
        ----------
        ham_2d : ndarray, shape (nvol, ndx, ndx)
            Interaction Hamiltonian reshaped to 2D matrices.

        Returns
        -------
        list of list of int, or None
            List of index groups forming independent blocks.
            Returns None if no block structure is found (single block).
        """
        ndx = ham_2d.shape[-1]
        if ndx <= 1:
            return None

        # Sum absolute values over nvol to get connectivity pattern
        connectivity = np.sum(np.abs(ham_2d), axis=0)

        # Build adjacency from non-zero off-diagonal entries
        threshold = 1.0e-12
        adj = (np.abs(connectivity) > threshold) | (np.abs(connectivity.T) > threshold)

        # Union-Find with path halving for connected components
        parent = list(range(ndx))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]  # path halving
                x = parent[x]
            return x

        rows, cols = np.where(np.triu(adj, k=1))
        for i, j in zip(rows, cols):
            ri, rj = find(i), find(j)
            if ri != rj:
                if ri > rj:
                    ri, rj = rj, ri
                parent[rj] = ri

        # Collect components
        components = {}
        for i in range(ndx):
            r = find(i)
            if r not in components:
                components[r] = []
            components[r].append(i)

        if len(components) <= 1:
            return None

        return list(components.values())


def run(*, input_dict: Optional[dict] = None, input_file: Optional[str] = None):
    """Main entry point for running RPA calculations.

    Parameters
    ----------
    input_dict : dict, optional
        Dictionary containing input parameters. Must be provided if input_file is None.
    input_file : str, optional
        Path to input file. Not used if input_dict is provided.

    Raises
    ------
    RuntimeError
        If input_dict is None.

    Notes
    -----
    The input dictionary should contain the following sections:
    - log: Logging configuration
        - print_level: Verbosity level (default: 1)
        - print_step: Print frequency (default: 1)
    - mode: Calculation mode configuration
        - mode: Must be "RPA" for RPA calculations
    - file: File paths configuration
        - input: Input file paths
            - path_to_input: Input directory (default: "")
        - output: Output file paths
            - path_to_output: Output directory (default: "output")

    The function performs the following steps:
    1. Initializes logging and file paths
    2. Reads interaction definitions
    3. Creates RPA solver instance
    4. Performs RPA calculations
    5. Saves results
    """
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
