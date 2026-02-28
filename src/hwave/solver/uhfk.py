import itertools
import logging
import numpy as np
import os
from .base import solver_base
from .perf import do_profile

logger = logging.getLogger("qlms").getChild("uhfk")

class UHFk(solver_base):
    @do_profile
    def __init__(self, param_ham, info_log, info_mode, param_mod=None):
        """Initialize the UHF solver.
        
        Parameters
        ----------
        param_ham : dict
            Hamiltonian parameters
        info_log : dict
            Logging configuration  
        info_mode : dict
            Solver mode parameters
        param_mod : dict, optional
            Additional model parameters
        """
        self.name = "uhfk"
        super().__init__(param_ham, info_log, info_mode, param_mod)

        # Initialize solver components
        self._init_mode(info_mode)  # Set solver modes like Fock term
        self._init_param()          # Initialize parameters
        self._init_lattice()        # Setup lattice geometry
        self._init_orbit()          # Setup orbital structure
        # self._dump_param_ham()
        self._check_interaction()   # Validate interactions
        self._init_interaction()    # Setup interaction terms
        # self._dump_param_ham()
        self._init_wavevec()        # Setup k-vectors

        self._show_param()          # Display parameters

        # Initialize physical quantities
        self.physics = {
            "Ene": { "Total": 0.0, "Band": 0.0 },
            "NCond": 0.0,  # Total particle number
            "Sz": 0.0,     # Total spin
            "Rest": 1.0    # Convergence measure
        }

        # Setup dimensions
        nvol = self.nvol  # Number of k-points
        nd = self.nd      # Total number of basis states

        # Storage for eigenvalues/vectors
        self._green_list = {
#            "eigenvalue":  np.zeros((nvol,nd), dtype=np.complex128),
#            "eigenvector": np.zeros((nvol,nd,nd), dtype=np.complex128)
        }

        #self.Green = np.zeros((nvol,ns,norb,ns,norb), dtype=np.complex128),

        # Initialize Hamiltonian matrix
        self.ham = np.zeros((nvol,nd,nd), dtype=np.complex128)

    def _init_mode(self, param):
        """Initialize solver modes.
        
        Parameters
        ----------
        param : dict
            Mode parameters
        """
        # Enable/disable Fock term
        self.iflag_fock = self._to_bool(param.get('flag_fock', True), 'flag_fock')
        
        # Enable/disable spin-orbital coupling
        self.enable_spin_orbital = self._to_bool(
            param.get('enable_spin_orbital', False), 'enable_spin_orbital'
        )

    def _to_bool(self, value, param_name):
        """Convert a value to boolean with type validation.
        
        Parameters
        ----------
        value : any
            Value to convert to boolean
        param_name : str
            Parameter name for error messages
            
        Returns
        -------
        bool
            Converted boolean value
            
        Raises
        ------
        ValueError
            If value cannot be converted to boolean
        """
        if isinstance(value, bool):
            return value
        if isinstance(value, str):
            if value.lower() in ('true', 'yes', '1', 'on'):
                logger.warning(
                    f"Parameter '{param_name}' should be boolean, not string. "
                    f"Use '{param_name} = true' instead of '{param_name} = \"{value}\"'"
                )
                return True
            elif value.lower() in ('false', 'no', '0', 'off'):
                logger.warning(
                    f"Parameter '{param_name}' should be boolean, not string. "
                    f"Use '{param_name} = false' instead of '{param_name} = \"{value}\"'"
                )
                return False
            else:
                raise ValueError(
                    f"Parameter '{param_name}' has invalid value '{value}'. "
                    f"Expected boolean (true/false)."
                )
        if isinstance(value, (int, float)):
            if value == 0:
                return False
            elif value == 1:
                return True
            else:
                raise ValueError(
                    f"Parameter '{param_name}' has invalid numeric value {value}. "
                    f"Expected 0 or 1 (or boolean true/false)."
                )
        raise ValueError(
            f"Parameter '{param_name}' has invalid type {type(value).__name__}. "
            f"Expected boolean."
        )

    @do_profile
    def _init_param(self):
        """Initialize solver parameters.
        
        Sets up key parameters including:
        - Cell shape
        - Temperature
        - Number of electrons/filling
        - Sz constraints
        - Numerical precision parameters
        """
        # check and store parameters

        # - cell shape
        if not "CellShape" in self.param_mod:
            logger.error("CellShape is missing. abort")
            exit(1)

        # - temperature
        self.T = self.param_mod.get("T", 0.0)

        # - number of electrons
        if "filling" in self.param_mod:
            self.filling = self.param_mod["filling"]
            round_mode = self.param_mod.get("Ncond_round_mode", "strict")

            Lx,Ly,Lz = self.param_mod.get("CellShape")
            norb = self.param_ham["Geometry"]["norb"]
            if self.enable_spin_orbital:
                Nstate = Lx*Ly*Lz*norb  # norb already includes spin
            else:
                Nstate = Lx*Ly*Lz*norb*2

            ncond = self._round_to_int(Nstate * self.filling, round_mode)
            self.param_mod["Ncond"] = ncond  # overwrite
        else:
            self.filling = None
            ncond = self.param_mod["Ncond"]

        # - Sz fixed or free
        if self.param_mod["2Sz"] is None:
            self.sz_free = True
            self.Nconds = [ ncond ]
        else:
            self.sz_free = False
            twosz = self.param_mod["2Sz"]
            self.Nconds = [ (ncond + twosz)//2, (ncond - twosz)//2 ]

        # - calc condition
        self.ene_cutoff = self.param_mod.get("ene_cutoff", 1e+2)

        # - strict hermiticity check
        self.strict_hermite = self.param_mod.get("strict_hermite", False)
        self.hermite_tolerance = self.param_mod.get("hermite_tolerance", 1.0e-8)

        # - check strictness
        if self.relax_checks:
            self.strict_hermite = False

    @do_profile
    def _init_lattice(self):
        """Initialize lattice geometry.
        
        Sets up:
        - Cell shape and volume
        - Sublattice structure
        - Consistency checks for sublattice
        
        Parameters are read from self.param_mod.
        """
        Lx,Ly,Lz = self.param_mod.get("CellShape")
        self.cellshape = (Lx,Ly,Lz)
        self.cellvol = Lx * Ly * Lz

        Bx,By,Bz = self.param_mod.get("SubShape", [Lx,Ly,Lz])
        self.subshape = (Bx,By,Bz)
        self.subvol = Bx * By * Bz

        self.has_sublattice = (self.subvol > 1)

        # check consistency
        # XXX use reciprocal lattice
        err = 0
        for i in range(3):
            if self.cellshape[i] % self.subshape[i] != 0:
                err += 1
        if err > 0:
            logger.error("SubShape is not compatible with CellShape. abort")
            exit(1)

        # replace by lattice of supercells
        nx, ny, nz = Lx//Bx, Ly//By, Lz//Bz
        nvol = nx * ny * nz

        self.shape = (nx, ny, nz)
        self.nvol = nvol

    @do_profile
    def _init_orbit(self):
        """Initialize orbital structure.

        Sets up:
        - Number of orbitals
        - Spin degrees of freedom
        - Total basis dimension
        Takes into account supercell structure if present.

        When enable_spin_orbital is True, the orbital index from input files
        already includes spin (Wannier90 format: index = 2*orb + spin).
        In this case, ns=1 since spin is encoded in the orbital index.
        """
        norb = self.param_ham["Geometry"]["norb"]

        if self.enable_spin_orbital:
            # Spin is already encoded in orbital index (Wannier90 SOI format)
            ns = 1
            if norb % 2 != 0:
                logger.error("In spin-orbital mode, norb must be even (got {})".format(norb))
                exit(1)
            self.norb_phys_orig = norb // 2
        else:
            # Spin is separate degree of freedom
            ns = 2
            self.norb_phys_orig = norb

        # take account of supercell
        self.norb_orig = norb
        self.norb = norb * self.subvol

        self.nd = self.norb * ns
        self.ns = ns

        if self.enable_spin_orbital:
            self.norb_phys = self.norb_phys_orig * self.subvol
        else:
            self.norb_phys = self.norb

    @do_profile
    def _show_param(self):
        """Display solver parameters.
        
        Logs key parameters including:
        - Enabled features (Fock term, spin-orbital coupling)
        - Lattice dimensions
        - Orbital/spin structure
        - Physical parameters (filling, temperature)
        - Numerical parameters
        """
        logger.info("Show parameters")
        logger.info("    Enable Fock    = {}".format(self.iflag_fock))
        logger.info("    Cell Shape     = {}".format(self.cellshape))
        logger.info("    Sub Shape      = {}".format(self.subshape))
        logger.info("    Block          = {}".format(self.shape))
        logger.info("    Block volume   = {}".format(self.nvol))
        logger.info("    Num orbit      = {}".format(self.norb_orig))
        logger.info("    Num orbit eff  = {}".format(self.norb))
        logger.info("    nspin          = {}".format(self.ns))
        logger.info("    nd             = {}".format(self.nd))

        logger.info("    Ncond          = {}".format(self.Nconds))
        if self.filling is not None:
            logger.info("    filling        = {}".format(self.filling))
        logger.info("    T              = {}".format(self.T))
        logger.info("    E_cutoff       = {}".format(self.ene_cutoff))

        logger.info("    Mix            = {}".format(self.param_mod["Mix"]))
        logger.info("    RndSeed        = {}".format(self.param_mod["RndSeed"]))
        logger.info("    IterationMax   = {}".format(self.param_mod["IterationMax"]))
        logger.info("    EPS            = {}".format(self.param_mod["EPS"]))
        logger.info("    2Sz            = {}".format(self.param_mod["2Sz"]))
        logger.info("    strict_hermite = {}".format(self.strict_hermite))
        logger.info("    hermite_tol    = {}".format(self.hermite_tolerance))

        logger.info("    spin_orbital   = {}".format(self.enable_spin_orbital))

    @do_profile
    def _reshape_geometry(self, geom):
        """Reshape geometry for sublattice structure.
        
        Parameters
        ----------
        geom : dict
            Original geometry definition
        
        Returns
        -------
        dict
            Reshaped geometry for sublattice
        """
        Bx,By,Bz = self.subshape
        bvol = Bx * By * Bz

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

    @do_profile
    def _reshape_interaction(self, ham, enable_spin_orbital):
        """Reshape interaction terms for sublattice.
        
        Parameters
        ----------
        ham : dict
            Original interaction terms
        enable_spin_orbital : bool
            Whether spin-orbital coupling is enabled
        
        Returns
        -------
        dict
            Reshaped interaction terms
        """
        Bx,By,Bz = self.subshape
        nx,ny,nz = self.shape

        norb_orig = self.param_ham_orig["Geometry"]["norb"]

        def _reshape_orbit_(a, x):
            return a + self.norb_orig * ( x[0] + Bx * (x[1] + By * (x[2])))

        def _reshape_orbit_spin(a, x):
            a_, s_ = a%norb_orig, a//norb_orig
            return a_ + norb_orig * ( x[0] + Bx * (x[1] + By * (x[2] + Bz * s_)))

        if enable_spin_orbital:
            _reshape_orbit = _reshape_orbit_spin
        else:
            _reshape_orbit = _reshape_orbit_

        def _round(x, n):
            return x % n

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

    @do_profile
    def _reshape_green(self, green):
        """Convert Green function to sublattice basis.
        
        Parameters
        ----------
        green : ndarray
            Green function in original basis
        
        Returns
        -------
        ndarray
            Green function in sublattice basis
        """
        # convert green function into sublattice

        Lx,Ly,Lz = self.cellshape
        Lvol = Lx * Ly * Lz
        Bx,By,Bz = self.subshape
        Bvol = Bx * By * Bz
        Nx,Ny,Nz = self.shape
        Nvol = Nx * Ny * Nz

        norb_orig = self.norb_orig
        norb = self.norb
        ns = self.ns

        # check array size
        assert(green.shape == (Lvol,ns,norb_orig,ns,norb_orig))

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

        def _pack_site(x, n):
            _ix, _iy, _iz = x
            _nx, _ny, _nz = n
            return _iz + _nz * (_iy + _ny * (_ix))

        def _unpack_site(idx, n):
            _nx, _ny, _nz = n
            _iz = idx % _nz
            _iy = (idx // _nz) % _ny
            _ix = (idx // (_nz * _ny)) % _nx
            return (_ix,_iy,_iz)

        green_sub = np.zeros((Nvol,ns,norb,ns,norb), dtype=np.complex128)

        for isite in range(Nvol):
            ixx, iyy, izz = _unpack_site(isite, (Nx,Ny,Nz))
            ix0, iy0, iz0 = ixx * Bx, iyy * By, izz * Bz

            for aa, bb in itertools.product(range(norb), range(norb)):
                a, ri = aa % norb_orig, aa // norb_orig
                b, rj = bb % norb_orig, bb // norb_orig

                rix, riy, riz = _unpack_index(ri, (Bx,By,Bz))
                rjx, rjy, rjz = _unpack_index(rj, (Bx,By,Bz))

                ix = (ix0 + rjx - rix) % Lx
                iy = (iy0 + rjy - riy) % Ly
                iz = (iz0 + rjz - riz) % Lz

                jsite = _pack_site((ix,iy,iz), (Lx,Ly,Lz))

                for s, t in itertools.product(range(ns), range(ns)):
                    green_sub[isite, s, aa, t, bb] = green[jsite, s, a, t, b]

        return green_sub

    @do_profile
    def _deflate_green(self, green_sub):
        """Convert Green function back to original basis.
        
        Parameters
        ----------
        green_sub : ndarray
            Green function in sublattice basis
        
        Returns
        -------
        ndarray
            Green function in original basis
        """
        # convert green function back to original lattice
        Lx,Ly,Lz = self.cellshape
        Lvol = Lx * Ly * Lz
        Bx,By,Bz = self.subshape
        Bvol = Bx * By * Bz
        Nx,Ny,Nz = self.shape
        Nvol = Nx * Ny * Nz

        norb_orig = self.norb_orig
        norb = self.norb
        ns = self.ns

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

        def _pack_site(x, n):
            _ix, _iy, _iz = x
            _nx, _ny, _nz = n
            return _iz + _nz * (_iy + _ny * (_ix))

        def _unpack_site(idx, n):
            _nx, _ny, _nz = n
            _iz = idx % _nz
            _iy = (idx // _nz) % _ny
            _ix = (idx // (_nz * _ny)) % _nx
            return (_ix,_iy,_iz)

        green = np.zeros((Lvol,ns,norb_orig,ns,norb_orig), dtype=np.complex128)

        for jsite in range(Lvol):
            ix, iy, iz = _unpack_site(jsite, (Lx,Ly,Lz))

            ixx, irx = ix // Bx, ix % Bx
            iyy, iry = iy // By, iy % By
            izz, irz = iz // Bz, iz % Bz

            isite = _pack_site((ixx,iyy,izz), (Nx,Ny,Nz))
            ir = _pack_index((irx,iry,irz), (Bx,By,Bz))

            for a, b in itertools.product(range(norb_orig), range(norb_orig)):
                aa = a
                bb = b + norb_orig * ir
                for s, t in itertools.product(range(ns), range(ns)):
                    green[jsite, s, a, t, b] = green_sub[isite, s, aa, t, bb]

        return green

    @do_profile
    def _init_interaction(self):
        """Initialize interaction terms.
        
        Processes interaction terms for sublattice structure if needed.
        Handles different types of interactions (Coulomb, Hund, etc).
        """
        # reinterpret interaction coefficient on sublattice
        if self.has_sublattice:
            # backup
            import copy
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

    def _init_wavevec(self):
        """Initialize wave vectors.
        
        Sets up:
        - k-vectors for sublatticed geometry
        - Reciprocal lattice vectors
        - Wave number tables
        """
        # wave vectors on sublatticed geometry
        def _klist(n):
            return np.roll( (np.arange(n)-(n//2)), -(n//2) )

        geom = self.param_ham["Geometry"]
        rvec = geom["rvec"]
        omg = np.dot(rvec[0], np.cross(rvec[1], rvec[2]))
        kvec = np.array([
            np.cross(rvec[(i+1)%3], rvec[(i+2)%3])/omg * 2*np.pi/self.shape[i]
            for i in range(3) ])

        self.kvec = kvec  # store reciprocal lattice vectors

        nx,ny,nz = self.shape
        nvol = self.nvol

        self.wavenum_table = np.array([(i,j,k) for i in _klist(nx) for j in _klist(ny) for k in _klist(nz)])

        wtable = np.zeros((nx,ny,nz,3), dtype=float)
        for ix, kx in enumerate(_klist(nx)):
            for iy, ky in enumerate(_klist(ny)):
                for iz, kz in enumerate(_klist(nz)):
                    v = kvec[0] * kx + kvec[1] * ky + kvec[2] * kz
                    wtable[ix,iy,iz] = v
        self.wave_table = wtable.reshape(nvol,3)

    @do_profile
    def _dump_param_ham(self):
        if logger.getEffectiveLevel() < logging.INFO:
            return

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

    @do_profile
    def _check_interaction(self):
        """Validate interaction terms.
        
        Performs checks on:
        - Cell size compatibility
        - Orbital index validity
        - Hermiticity of terms
        - Spin-orbital mode compatibility
        """
        self._check_cellsize()
        self._check_orbital_index()
        self._check_hermite()
        self._check_spin_orbital_compatibility()

    def _check_spin_orbital_compatibility(self):
        """Check if interaction terms are compatible with spin-orbital mode.

        When enable_spin_orbital=True, interaction terms are supported by
        converting the spin-orbital Green function to a virtual (ns=2, norb_phys)
        form for the interaction contractions. The interaction file orbital indices
        must refer to physical orbitals [0, norb/2).
        """
        if not self.enable_spin_orbital:
            return

        supported_types = [
            "CoulombIntra", "CoulombInter", "Coulomb",
            "Hund", "Ising", "PairLift", "Exchange", "PairHop"
        ]

        found_interactions = []
        for itype in supported_types:
            if itype in self.param_ham and self.param_ham[itype]:
                found_interactions.append(itype)

        if found_interactions:
            logger.info(
                "Spin-orbital mode: interaction terms {} will be handled "
                "using virtual spin decomposition (norb_phys={})".format(
                    found_interactions, self.norb_phys_orig)
            )

    def _check_cellsize(self):
        err = 0
        for k in self.param_ham.keys():
            if k in ["Geometry", "Initial"]:
                pass
            else:
                min_r, max_r = [0,0,0], [0,0,0]
                for (irvec,orbvec), v in self.param_ham[k].items():
                    for i in range(3):
                        min_r[i] = irvec[i] if irvec[i] < min_r[i] else min_r[i]
                        max_r[i] = irvec[i] if irvec[i] > max_r[i] else max_r[i]
                width_r = [ max_r[i]-min_r[i]+1 for i in range(3) ]

                logger.debug(
                    k+": range=({},{})({},{})({},{}), width={},{},{}, cellshape={},{},{}".format(
                        min_r[0], max_r[0],
                        min_r[1], max_r[1],
                        min_r[2], max_r[2],
                        width_r[0], width_r[1], width_r[2],
                        self.cellshape[0], self.cellshape[1], self.cellshape[2]
                ))

                if all([ width_r[i] <= self.cellshape[i] for i in range(3) ]):
                    logger.debug("range check for {} ok.".format(k))
                else:
                    err += 1
                    msg = "range check for {} failed.".format(k)
                    if self.relax_checks:
                        logger.warning(msg)
                    else:
                        logger.error(msg)
        if err > 0:
            msg = "_check_cellsize failed. interaction range exceeds cell shape."
            if self.relax_checks:
                logger.warning(msg)
            else:
                logger.error(msg)
                exit(1)

    def _check_orbital_index(self):
        norb = self.param_ham["Geometry"]["norb"]
        # In spin-orbital mode, interaction orbital indices refer to physical orbitals
        norb_inter = self.norb_phys_orig if self.enable_spin_orbital else norb
        err = 0
        for k in self.param_ham.keys():
            if k in ["Geometry", "Initial"]:
                pass
            elif k == "Transfer":
                fail = 0
                for (irvec,orbvec), v in self.param_ham[k].items():
                    # allow spin-orbital interaction, i.e. twice norb
                    if not all([ 0 <= orbvec[i] < norb*2 for i in range(2) ]):
                        fail += 1
                if fail > 0:
                    logger.error("orbital index check failed for {}.".format(k))
                    err += 1
                else:
                    logger.debug("orbital index check for {} ok.".format(k))
            elif k == "CoulombIntra":
                fail = 0
                for (irvec,orbvec), v in self.param_ham[k].items():
                    if not all([ 0 <= orbvec[i] < norb_inter for i in range(2) ]):
                        fail += 1
                for (irvec,orbvec), v in self.param_ham[k].items():
                    if not orbvec[0] == orbvec[1]:
                        fail += 1
                if fail > 0:
                    logger.error("orbital index check failed for {}.".format(k))
                    err += 1
                else:
                    logger.debug("orbital index check for {} ok.".format(k))
            else:
                fail = 0
                for (irvec,orbvec), v in self.param_ham[k].items():
                    if not all([ 0 <= orbvec[i] < norb_inter for i in range(2) ]):
                        fail += 1
                if fail > 0:
                    logger.error("orbital index check failed for {}.".format(k))
                    err += 1
                else:
                    logger.debug("orbital index check for {} ok.".format(k))
        if err > 0:
            logger.error("_check_orbital_index failed. invalid orbital index found in interaction definitions.")
            exit(1)

    def _check_hermite(self):
        max_print = 32

        if "Transfer" in self.param_ham:
            nx,ny,nz = self.param_mod.get("CellShape")
            norb = self.norb

            tab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)
            for (irvec,orbvec), v in self.param_ham["Transfer"].items():
                tab_r[(*irvec, *orbvec)] = v

            t = np.conjugate(
                    np.transpose(
                        np.flip(np.roll(tab_r, -1, axis=(0,1,2)), axis=(0,1,2)),
                        (0,1,2,4,3)
                    )
                )

            if not np.allclose(t, tab_r, atol=self.hermite_tolerance, rtol=0.0):
                msg = "Hermiticity check failed: |T_ba(-r)^* - T_ab(r)| = {}".format(np.sum(np.abs(t - tab_r)))
                if self.strict_hermite:
                    logger.error(msg)
                    errlist = np.array(np.nonzero(~np.isclose(t, tab_r, atol=self.hermite_tolerance, rtol=0.0))).transpose()
                    for idx in errlist[:max_print]:
                        logger.error("    index={}, T_ab={}, T_ba^*={}".format(idx, tab_r[tuple(idx)], t[tuple(idx)]))
                    if len(errlist) > max_print:
                        logger.error("    ...")
                    logger.error("{} entries".format(len(errlist)))
                    exit(1)
                else:
                    logger.warning(msg)

    @do_profile
    def solve(self, green_info, path_to_output):
        """Main solver routine.
        
        Performs UHF iteration until convergence:
        1. Constructs Hamiltonian
        2. Diagonalizes
        3. Updates Green function
        4. Calculates physical quantities
        
        Parameters
        ----------
        green_info : dict
            Initial Green function information
        path_to_output : str
            Path for output files
        """
        print_level = self.info_log["print_level"]
        print_check = self.info_log.get("print_check", None)

        logger.info("Start UHFk calculations")

        if print_level > 0:
            logger.info("step, rest, energy, NCond, Sz")
        if print_check is not None:
            fch = open(os.path.join(path_to_output, print_check), "w")

        self._make_ham_trans() # T_ab(k)
        self._make_ham_inter() # J_ab(r)
        self._detect_blocks()  # Detect block-diagonal structure

        self.Green = self._initial_green(green_info)

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
            if print_check is not None:
                fch.write("{}, {:.8g}, {:.8g}, {:.4g}, {:.4g}\n"
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
            logger.info("Total Energy = {}.".format(self.physics["Ene"]["Total"]))
        else:
            logger.info("UHFk calculation failed: rest={}, eps={}."
                        .format(self.physics["Rest"], self.param_mod["eps"]))
        if print_check is not None:
            fch.close()

    @do_profile
    def _initial_green(self, green_info):
        logger.debug(">>> _initial_green")

        _data = None
        if "initial" in green_info and green_info["initial"] is not None:
            logger.info("load initial green function from file")
            _data = self._read_green_from_data(green_info["initial"])
        elif "initial_uhf" in green_info:
            if not "geometry_uhf" in green_info:
                logger.error("initial green function in coord space requires geometry.dat")
                exit(1)
            logger.info("load initial green function in coord space from file")
            _data = self._initial_green_uhf(green_info["initial_uhf"], green_info["geometry_uhf"])
        elif "initial_mode" in green_info:
            if green_info["initial_mode"] == "zero":
                logger.info("initialize green function with zeros")
                _data = self._initial_green_zero()
            elif green_info["initial_mode"] in ["one", "unity"]:
                logger.info("initialize green function with identity")
                _data = self._initial_green_one()
            elif green_info["initial_mode"] == "conventional_random":
                logger.info("initialize green function with random numbers (conventional)")
                _data = self._initial_green_random_conv()
                #_data = self._initial_green_random_conv_reshape()
            elif green_info["initial_mode"] == "random":
                logger.info("initialize green function with random numbers")
                _data = self._initial_green_random_simple()

        if _data is None: # default = zero
            logger.info("initialize green function with zeros")
            _data = self._initial_green_zero()
        return _data

    def _initial_green_zero(self):
        nvol     = self.nvol
        norb     = self.norb
        ns       = self.ns
        nd = ns * norb

        green = np.zeros((nvol,ns,norb,ns,norb), dtype=np.complex128)
        return green

    def _initial_green_one(self):
        nvol     = self.nvol
        norb     = self.norb
        ns       = self.ns
        nd = ns * norb

        green = np.zeros((nvol,ns,norb,ns,norb), dtype=np.complex128)
        green[0] = np.einsum('ab,st->satb', np.identity(norb), np.identity(ns))
        return green

    def _initial_green_random_simple(self):
        nvol     = self.nvol
        norb     = self.norb
        ns       = self.ns
        nd = ns * norb

        np.random.seed(self.param_mod["RndSeed"])
        rand = np.random.rand(nvol * nd * nd).reshape(nvol,ns,norb,ns,norb)
        green = 0.01 * (rand - 0.5)

        return green

    def _initial_green_random_conv_reshape(self):
        lx,ly,lz = self.cellshape
        lvol     = self.cellvol
        norb     = self.norb_orig if self.has_sublattice else self.norb
        ns       = self.ns
        nd = ns * norb

        np.random.seed(self.param_mod["RndSeed"])
        #rand = np.random.rand(nvol * nd * nd).reshape(nvol,ns,norb,ns,norb)
        #green = 0.01 * (rand - 0.5)

        # G[(s,i,a),(t,j,b)] -> G[ij,s,a,t,b]
        x1 = np.random.rand(lvol * nd * lvol * nd).reshape(ns,lvol,norb,ns,lvol,norb)
        x2 = np.transpose(x1,(1,4,0,2,3,5))
        x3 = np.average(x2,axis=0).reshape(lvol,ns,norb,ns,norb)

        green = 0.01 * (x3 - 0.5)

        if self.has_sublattice:
            return self._reshape_green(green)
        else:
            return green

    def _initial_green_random_conv(self):
        nvol     = self.nvol
        norb     = self.norb
        ns       = self.ns
        nd = ns * norb

        np.random.seed(self.param_mod["RndSeed"])
        #rand = np.random.rand(nvol * nd * nd).reshape(nvol,ns,norb,ns,norb)
        #green = 0.01 * (rand - 0.5)

        # G[(s,i,a),(t,j,b)] -> G[ij,s,a,t,b]
        x1 = np.random.rand(nvol * nd * nvol * nd).reshape(ns,nvol,norb,ns,nvol,norb)
        x2 = np.transpose(x1,(1,4,0,2,3,5))
        x3 = np.average(x2,axis=0).reshape(nvol,ns,norb,ns,norb)

        green = 0.01 * (x3 - 0.5)

        return green

    def _initial_green_uhf(self, info, geom):
        lx,ly,lz = self.cellshape
        lvol     = self.cellvol
        norb     = self.norb
        nd       = self.nd
        ns       = self.ns

        def _pack_site(ix,iy,iz):
            return ix + lx * (iy + ly * iz)

        tbl = {}
        for idx, (ix,iy,iz,a) in geom["site2vec"].items():
            site = _pack_site(ix,iy,iz)
            tbl[idx] = (site, a)

        green_uhf = np.zeros((lvol,lvol,ns,norb,ns,norb),dtype=np.complex128)
        for (i,s,j,t),v in info.items():
            isite, a = tbl[i]
            jsite, b = tbl[j]
            green_uhf[isite,jsite,s,a,t,b] = v

        if self.has_sublattice:
            return self._reshape_green(green_uhf[0])
        else:
            return green_uhf[0]

    @do_profile
    def _make_ham_trans(self):
        logger.debug(">>> _make_ham_trans")

        nx,ny,nz = self.shape
        nvol     = self.nvol
        norb     = self.norb
        nd       = self.nd

        if self.enable_spin_orbital == True:
            tab_r = np.zeros((nx,ny,nz,nd,nd), dtype=np.complex128)

            for (irvec,orbvec), v in self.param_ham["Transfer"].items():
                tab_r[(*irvec, *orbvec)] = v

            # fourier transform
            tab_k = np.fft.ifftn(tab_r, axes=(0,1,2), norm='forward')

            # Reshape from (nx,ny,nz,nd,nd) to (nvol,nd,nd)
            ham = tab_k.reshape(nvol, nd, nd)

        else:
            # data structure of T_ab(r): T(rx,ry,rz,a,b)
            tab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)

            has_spin_dep = False
            for (irvec,orbvec), v in self.param_ham["Transfer"].items():
                if orbvec[0] < norb and orbvec[1] < norb:
                    tab_r[(*irvec, *orbvec)] = v
                else:
                    has_spin_dep = True
            if has_spin_dep:
                logger.warning(
                    "Transfer has orbital indices >= norb (spin-dependent terms) "
                    "but enable_spin_orbital is False. "
                    "These terms are ignored. Set enable_spin_orbital = true to use them."
                )

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

    @do_profile
    def _make_ham_inter(self):
        logger.debug(">>> _make_ham_inter")

        nx,ny,nz = self.shape
        nvol     = self.nvol
        # In spin-orbital mode, interaction arrays use physical orbital count
        norb     = self.norb_phys if self.enable_spin_orbital else self.norb
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
            self.inter_table["CoulombInter"] = (vab_r + vba)/2
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
                self.inter_table["CoulombInter"] = (vab_r + vba)/2
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
            self.inter_table["Hund"] = -(jab_r + jba)/2
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

            # interaction coeffs: J_ij Sz_i Sz_j where Sz = 1/2 sigma, sigma=+1,-1
            self.inter_table["Ising"] = (jab_r + jba)/2 / 4
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
            self.inter_table["PairLift"] = (jab_r + jba)/2
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
            self.inter_table["Exchange"] = -(jab_r + jba)/2
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
            self.inter_table["PairHop"] = (jab_r + jba)/2
            # spin combination
            self.spin_table["PairHop"] = np.zeros((2,2,2,2), dtype=int)
            self.spin_table["PairHop"][0,1,1,0] = 1
            self.spin_table["PairHop"][1,0,0,1] = 1
        else:
            self.inter_table["PairHop"] = None

    @do_profile
    def _detect_blocks(self):
        """Detect block-diagonal structure from transfer and interaction terms.

        Analyzes the connectivity in the (spin, orbital) index space by
        combining the non-zero patterns of:
        1. Transfer Hamiltonian ham_trans (summed over k-points)
        2. Interaction terms inter_table x spin_table (expanded to nd x nd)

        The detected blocks are used by _diag() and _green() to decompose
        the problem into smaller independent sub-problems.

        Sets self.block_info: list of arrays, each containing nd-space indices
        for one independent block. Also sets self.group_nconds and
        self.block_to_group for global chemical potential determination.
        """
        logger.debug(">>> _detect_blocks")

        nd = self.nd
        norb = self.norb
        ns = self.ns

        # Build connectivity matrix in nd x nd space
        connectivity = np.zeros((nd, nd), dtype=np.float64)

        # 1. Transfer Hamiltonian: sum |ham_trans(k)| over k
        connectivity += np.sum(np.abs(self.ham_trans), axis=0)

        # 2. Interaction terms: expand (norb_inter, norb_inter) x spin(2,2,2,2) -> (nd, nd)
        norb_inter = self.norb_phys if self.enable_spin_orbital else norb

        for type_name in self.inter_table:
            if self.inter_table[type_name] is not None:
                jab = self.inter_table[type_name]
                spin = self.spin_table[type_name]

                # Sum |J_ab(r)| over r to get orbital connectivity
                jab_sum = np.sum(np.abs(jab.reshape(-1, norb_inter, norb_inter)), axis=0)

                # For Hartree: spin_mat_h[s,t] = any(spin[s,:,:,t] != 0)
                spin_mat_h = np.any(spin != 0, axis=(1, 2)).astype(float)
                orb_connected_h = (jab_sum > 0).any(axis=1).astype(float)

                if self.enable_spin_orbital:
                    # In spin-orbital mode: (s, a) -> so_idx = 2*a + s
                    for s in range(2):
                        for t in range(2):
                            if spin_mat_h[s, t] > 0:
                                for a in range(norb_inter):
                                    if orb_connected_h[a]:
                                        connectivity[2 * a + s, 2 * a + t] += 1.0
                    if self.iflag_fock:
                        spin_fock = np.max(np.abs(spin), axis=(1, 3))
                        for s in range(2):
                            for t in range(2):
                                if spin_fock[s, t] > 0:
                                    for a in range(norb_inter):
                                        for b in range(norb_inter):
                                            if jab_sum[a, b] > 0:
                                                connectivity[2 * a + s, 2 * b + t] += jab_sum[a, b]
                else:
                    # Normal mode: (s, a) -> nd_idx = s * norb + a
                    for s in range(ns):
                        for t in range(ns):
                            if spin_mat_h[s, t] > 0:
                                for a in range(norb):
                                    if orb_connected_h[a]:
                                        connectivity[s * norb + a, t * norb + a] += 1.0

                    if self.iflag_fock:
                        spin_mat_f = np.any(
                            spin.transpose(0, 2, 1, 3).reshape(ns * ns, ns * ns) != 0
                        ).astype(float) if ns > 1 else np.ones((1, 1))
                        spin_fock = np.max(np.abs(spin), axis=(1, 3))
                        for s in range(ns):
                            for t in range(ns):
                                if spin_fock[s, t] > 0:
                                    connectivity[s * norb:(s + 1) * norb,
                                                 t * norb:(t + 1) * norb] += jab_sum

        # Symmetrize
        connectivity = connectivity + connectivity.T

        # Find connected components
        threshold = 1.0e-12
        adj = np.abs(connectivity) > threshold

        # BFS/label propagation
        labels = np.arange(nd)
        changed = True
        while changed:
            changed = False
            for i in range(nd):
                neighbors = np.where(adj[i])[0]
                if len(neighbors) > 0:
                    min_label = min(labels[i], np.min(labels[neighbors]))
                    if labels[i] != min_label:
                        labels[i] = min_label
                        changed = True
                    for n in neighbors:
                        if labels[n] != min_label:
                            labels[n] = min_label
                            changed = True

        unique_labels = np.unique(labels)
        blocks = [np.where(labels == lbl)[0] for lbl in unique_labels]

        # Sort blocks by first index
        blocks.sort(key=lambda b: b[0])

        # Group blocks by spin sector for shared chemical potential.
        # Instead of distributing Ncond per block (which can violate the
        # global energy minimum), we group blocks that share the same
        # electron reservoir and find a single mu per group.
        if self.sz_free:
            # All blocks share one global mu
            # group_nconds[g] = target electron count for group g
            # block_to_group[b] = which group block b belongs to
            block_to_group = [0] * len(blocks)
            group_nconds = [self.Nconds[0]]
        else:
            # Classify blocks by spin content
            ncond_up, ncond_down = self.Nconds[0], self.Nconds[1]
            block_to_group = []
            # group 0 = up-only, group 1 = down-only, group 2 = mixed
            group_nconds_dict = {}  # group_id -> ncond
            for blk in blocks:
                if self.enable_spin_orbital:
                    n_up = np.sum(blk % 2 == 0)
                    n_down = np.sum(blk % 2 == 1)
                else:
                    n_up = np.sum(blk < norb)
                    n_down = np.sum(blk >= norb)
                if n_up > 0 and n_down == 0:
                    block_to_group.append(0)
                    group_nconds_dict[0] = ncond_up
                elif n_up == 0 and n_down > 0:
                    block_to_group.append(1)
                    group_nconds_dict[1] = ncond_down
                else:
                    # Mixed spin block: use total ncond
                    block_to_group.append(2)
                    group_nconds_dict[2] = ncond_up + ncond_down
            # Build ordered list: group_nconds[g] for each unique group
            unique_groups = sorted(set(block_to_group))
            group_remap = {g: i for i, g in enumerate(unique_groups)}
            block_to_group = [group_remap[g] for g in block_to_group]
            group_nconds = [group_nconds_dict[g] for g in unique_groups]

        self.block_info = blocks
        self.block_to_group = block_to_group
        self.group_nconds = group_nconds

        # Log detected structure
        if len(blocks) == 1:
            logger.info("Block detection: single block (nd={})".format(nd))
        else:
            block_desc = ", ".join(
                ["{}".format(len(b)) for b in blocks]
            )
            logger.info("Block detection: {} blocks of sizes [{}], "
                        "{} mu-group(s) with Nconds={}".format(
                            len(blocks), block_desc,
                            len(group_nconds), group_nconds))

    def _so_to_virtual_green(self, G_so):
        """Convert spin-orbital Green function to virtual (ns=2, norb_phys) form.

        In spin-orbital mode, index = 2*orb + spin (interleaved).
        This converts to the standard (nvol, 2, norb_phys, 2, norb_phys) form
        so existing interaction contractions can be reused.

        Parameters
        ----------
        G_so : ndarray, shape (nvol, 1, nd, 1, nd)
            Green function in spin-orbital basis

        Returns
        -------
        ndarray, shape (nvol, 2, norb_phys, 2, norb_phys)
            Green function in virtual spin-separate form
        """
        nvol = G_so.shape[0]
        norb_phys = self.norb_phys
        G_flat = G_so.reshape(nvol, self.nd, self.nd)
        G_virt = np.zeros((nvol, 2, norb_phys, 2, norb_phys), dtype=np.complex128)
        for s in range(2):
            for t in range(2):
                G_virt[:, s, :, t, :] = G_flat[:, s::2, t::2]
        return G_virt

    def _virtual_ham_to_so(self, H_virt):
        """Convert virtual (ns=2, norb_phys) Hamiltonian to spin-orbital form.

        Inverse of _so_to_virtual_green for the Hamiltonian matrix.
        Converts from (nvol, 2, norb_phys, 2, norb_phys) ordering
        where index = s*norb_phys + a, to spin-orbital ordering
        where index = 2*a + s.

        Parameters
        ----------
        H_virt : ndarray, shape (nvol, 2, norb_phys, 2, norb_phys)
            Hamiltonian in virtual spin-separate form

        Returns
        -------
        ndarray, shape (nvol, nd, nd)
            Hamiltonian in spin-orbital basis
        """
        nvol = H_virt.shape[0]
        norb_phys = self.norb_phys
        nd = self.nd
        H_so = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for s in range(2):
            for t in range(2):
                H_so[:, s::2, t::2] = H_virt[:, s, :, t, :]
        return H_so

    @do_profile
    def _make_ham(self):
        logger.debug(">>> _make_ham")

        nx,ny,nz = self.shape
        nvol     = self.nvol
        nd       = self.nd
        ns       = self.ns

        # In spin-orbital mode, use virtual (ns=2, norb_phys) form for interactions
        if self.enable_spin_orbital:
            norb_inter = self.norb_phys
            nd_virt = 2 * norb_inter
            gab_r = self._so_to_virtual_green(self.Green)
        else:
            norb_inter = self.norb
            nd_virt = nd
            gab_r = self.Green

        # green function G_{ab,st}(r) : gab_r(r,s,a,t,b)
        # diagonal part G_{bb,st}(r=0) : gbb(s,t,b)
        gbb = np.diagonal(gab_r, axis1=2, axis2=4)[0,:,:,:]

        # hamiltonian H_{ab,st}(k) : ham(k,(s,a),(t,b))
        ham = np.zeros((nvol,nd,nd), dtype=np.complex128)

        # transfer term  T_{ab}(k) (note convention)
        logger.debug("Transfer")
        ham += self.ham_trans

        # interaction term: Coulomb type
        for type in ['CoulombIntra', 'CoulombInter', 'Hund', 'Ising', 'PairLift', 'Exchange']:
            if self.inter_table[type] is not None:
                logger.debug(type)

                # coefficient of interaction term J_{ab}(r)
                jab_r = self.inter_table[type].reshape(nvol,norb_inter,norb_inter)
                # and its spin combination  Spin(s1,s2,s3,s4)
                spin = self.spin_table[type]

                # non-cross term
                #   sum_r J_{ab}(r) G_{bb,uv}(0) Spin{s,u,v,t}
                hh0 = np.einsum('uvb, suvt -> stb', gbb, spin)
                hh1 = np.einsum('rab, stb -> rsta', jab_r, hh0)
                hh2 = np.einsum('rsta, ab -> rsatb', hh1, np.eye(norb_inter, norb_inter))
                hh3 = np.sum(hh2, axis=0)  # shape: (2, norb_inter, 2, norb_inter)

                if self.enable_spin_orbital:
                    hh3_nd = self._virtual_ham_to_so(
                        np.broadcast_to(hh3.reshape(1, 2, norb_inter, 2, norb_inter),
                                       (nvol, 2, norb_inter, 2, norb_inter)).copy()
                    )
                else:
                    hh3_nd = np.broadcast_to(hh3.reshape(nd, nd), (nvol, nd, nd))

                ham += hh3_nd

                # cross term
                #   - sum_r J_{ab}(r) G_{ba,uv}(r) Spin{s,u,t,v} e^{ikr}
                if self.iflag_fock:
                    hh4 = np.einsum('rab, rubva, sutv -> rsatb', jab_r, gab_r, spin)

                    #   fourier transform: sum_r (*) e^{ikr}
                    hh5 = np.fft.ifftn(hh4.reshape(nx,ny,nz,nd_virt,nd_virt), axes=(0,1,2), norm='forward')

                    if self.enable_spin_orbital:
                        ham -= self._virtual_ham_to_so(
                            hh5.reshape(nvol, 2, norb_inter, 2, norb_inter)
                        )
                    else:
                        ham -= hh5.reshape(nvol, nd, nd)

        # interaction term: PairHop type
        for type in ['PairHop']:
            if self.inter_table[type] is not None:
                logger.debug(type)

                # coefficient of interaction term J_{ab}(r)
                jab_r = self.inter_table[type].reshape(nvol,norb_inter,norb_inter)
                # and its spin combination  Spin(s1,s2,s3,s4)
                spin = self.spin_table[type]

                # non-cross and cross term
                #   + sum_r J_{ab}(r) G_{ab,uv}(-r) Spin{s,u,v,t} e^{ikr}
                #   - sum_r J_{ab}(r) G_{ab,uv}(-r) Spin{s,u,t,v} e^{ikr}

                if self.iflag_fock:
                    hh1 = np.einsum('rvbua, suvt -> rsbta', np.conjugate(gab_r), spin)
                    hh2 = np.einsum('rvbua, sutv -> rsbta', np.conjugate(gab_r), spin)
                    hh3 = np.einsum('rab, rsbta -> rsatb', jab_r, (hh1 - hh2))
                    hh4 = np.fft.ifftn(hh3.reshape(nx,ny,nz,nd_virt,nd_virt), axes=(0,1,2), norm='forward')
                else:
                    hh1 = np.einsum('rvbua, suvt -> rsbta', np.conjugate(gab_r), spin)
                    hh3 = np.einsum('rab, rsbta -> rsatb', jab_r, hh1)
                    hh4 = np.fft.ifftn(hh3.reshape(nx,ny,nz,nd_virt,nd_virt), axes=(0,1,2), norm='forward')

                if self.enable_spin_orbital:
                    ham += self._virtual_ham_to_so(
                        hh4.reshape(nvol, 2, norb_inter, 2, norb_inter)
                    )
                else:
                    ham += hh4.reshape(nvol, nd, nd)

        # Enforce block structure: zero out cross-block entries
        if len(self.block_info) > 1:
            mask = np.zeros((nd, nd), dtype=bool)
            for blk in self.block_info:
                idx = np.array(blk)
                ix = np.ix_(idx, idx)
                mask[ix] = True
            ham[:, ~mask] = 0.0

        # store
        self.ham = ham

    @do_profile
    def _diag(self):
        logger.debug(">>> _diag")

        nvol = self.nvol
        blocks = self.block_info

        # Diagonalize each block independently
        eigenvalues = []
        eigenvectors = []
        for blk in blocks:
            idx = np.array(blk)
            ix = np.ix_(idx, idx)
            mat_blk = self.ham[:, ix[0], ix[1]]  # (nvol, blk_size, blk_size)
            w, v = np.linalg.eigh(mat_blk)
            eigenvalues.append(w)
            eigenvectors.append(v)

        # Store as list of per-block arrays
        self._green_list["eigenvalue"] = eigenvalues
        self._green_list["eigenvector"] = eigenvectors

    @do_profile
    def _green(self):
        logger.debug(">>> _green")

        nx,ny,nz = self.shape
        nvol     = self.nvol
        nd       = self.nd
        norb     = self.norb
        ns       = self.ns

        ws_list = self._green_list["eigenvalue"]
        vs_list = self._green_list["eigenvector"]
        blocks = self.block_info
        nblock = len(blocks)
        block_to_group = self.block_to_group
        group_nconds = self.group_nconds
        ngroup = len(group_nconds)

        # Find a single global mu per group (shared across all blocks in the group)
        self._green_list["mu"] = np.zeros(ngroup, dtype=float)

        # Group blocks by mu-group
        group_blocks = [[] for _ in range(ngroup)]
        for b in range(nblock):
            group_blocks[block_to_group[b]].append(b)

        # Find mu for each group using all eigenvalues in that group
        group_dists = [None] * nblock
        for g in range(ngroup):
            blist = group_blocks[g]
            ws_group = [ws_list[b] for b in blist]
            vs_group = [vs_list[b] for b in blist]
            ncond = group_nconds[g]

            dists, mu = self._find_dist_group(ws_group, vs_group, ncond)

            self._green_list["mu"][g] = mu
            logger.debug("mu[group {}] = {}".format(g, mu))
            for i, b in enumerate(blist):
                group_dists[b] = dists[i]

        # Build full Green's function in k-space by assembling blocks
        gab_k = np.zeros((nvol, nd, nd), dtype=np.complex128)

        for k in range(nblock):
            v = vs_list[k]
            dist = group_dists[k]

            # G_ab(k) for this block
            gg_blk = np.einsum('kal, kl, kbl -> kab', np.conjugate(v), dist, v)

            # Place into full matrix
            idx = np.array(blocks[k])
            ix = np.ix_(idx, idx)
            gab_k[:, ix[0], ix[1]] = gg_blk

        # G_ab(r) = 1/V sum_k G_ab(k) e^{-ikr}
        gab_r = np.fft.fftn(gab_k.reshape(nx, ny, nz, nd, nd),
                            axes=(0, 1, 2), norm='forward')

        # store
        self.Green_prev = self.Green
        self.Green = gab_r.reshape(nvol, ns, norb, ns, norb)

    def _find_dist_group(self, ws_list, vs_list, ncond):
        """Find occupation distribution using a single mu across multiple blocks.

        Parameters
        ----------
        ws_list : list of ndarray
            Eigenvalues for each block in the group, shape (nvol, block_size)
        vs_list : list of ndarray
            Eigenvectors for each block, shape (nvol, block_size, block_size)
        ncond : int
            Total electron count for the entire group

        Returns
        -------
        dists : list of ndarray
            Occupation distribution for each block
        mu : float
            Chemical potential
        """
        if self.T == 0:
            return self._find_dist_group_zero_t(ws_list, vs_list, ncond)
        else:
            return self._find_dist_group_nonzero_t(ws_list, vs_list, ncond)

    def _find_dist_group_zero_t(self, ws_list, vs_list, ncond):
        def _ksq_table(width):
            nx,ny,nz = self.shape
            nvol = self.nvol

            kx = np.roll( (np.arange(nx) - nx//2), -nx//2) ** 2
            ky = np.roll( (np.arange(ny) - ny//2), -ny//2) ** 2
            kz = np.roll( (np.arange(nz) - nz//2), -nz//2) ** 2

            rr = np.zeros((nx,ny,nz))
            rr += np.broadcast_to(kx.reshape(nx,1,1),(nx,ny,nz))
            rr += np.broadcast_to(ky.reshape(1,ny,1),(nx,ny,nz))
            rr += np.broadcast_to(kz.reshape(1,1,nz),(nx,ny,nz))

            return np.broadcast_to(rr.reshape(nvol,1), (nvol,width))

        # Collect all eigenvalues across blocks with block labels
        all_ww = []
        all_ksq = []
        all_block_idx = []  # which block each eigenvalue belongs to
        all_local_idx = []  # flat index within that block's w array
        for b, w in enumerate(ws_list):
            k_sq = _ksq_table(w.shape[1]).flatten()
            ww = w.flatten()
            all_ww.append(ww)
            all_ksq.append(k_sq)
            all_block_idx.append(np.full(ww.size, b, dtype=int))
            all_local_idx.append(np.arange(ww.size))

        all_ww = np.concatenate(all_ww)
        all_ksq = np.concatenate(all_ksq)
        all_block_idx = np.concatenate(all_block_idx)
        all_local_idx = np.concatenate(all_local_idx)

        # Sort globally and pick lowest ncond states
        ev_idx = np.lexsort((all_ksq, all_ww))[0:ncond]

        # Distribute back to per-block arrays
        dists = [np.zeros(w.size) for w in ws_list]
        for idx in ev_idx:
            b = all_block_idx[idx]
            li = all_local_idx[idx]
            dists[b][li] = 1.0
        dists = [d.reshape(w.shape) for d, w in zip(dists, ws_list)]

        return dists, 0.0

    def _find_dist_group_nonzero_t(self, ws_list, vs_list, ncond):
        from scipy import optimize

        # Collect all eigenvalues for bracket search
        all_ev = np.sort(np.concatenate([w.flatten() for w in ws_list]))
        occupied_number = ncond

        def _fermi(t, mu, ev):
            w_ = (ev - mu) / t
            mask_ = w_ < self.ene_cutoff
            w1_ = np.where( mask_, w_, 0.0 )
            v1_ = 1.0 / (1.0 + np.exp(w1_))
            v_ = np.where( mask_, v1_, 0.0 )
            return v_

        def _calc_delta_n(mu):
            nn = 0.0
            for w, v in zip(ws_list, vs_list):
                ff = _fermi(self.T, mu, w)
                nn += np.einsum('kal,kl,kal->', np.conjugate(v), ff, v).real
            return nn - occupied_number

        # find mu s.t. <n>(mu) = N0
        is_converged = False
        if (_calc_delta_n(all_ev[0]) * _calc_delta_n(all_ev[-1])) < 0.0:
            logger.debug("+++ find mu: try bisection")
            mu, r = optimize.bisect(_calc_delta_n, all_ev[0], all_ev[-1], full_output=True, disp=False)
            is_converged = r.converged
        if not is_converged:
            logger.debug("+++ find mu: try newton")
            mu, r = optimize.newton(_calc_delta_n, all_ev[0], full_output=True)
            is_converged = r.converged
        if not is_converged:
            logger.error("find mu: not converged. abort")
            exit(1)

        dists = [_fermi(self.T, mu, w) for w in ws_list]
        return dists, mu

    @do_profile
    def _calc_phys(self):
        logger.debug(">>> _calc_phys")

        nx,ny,nz = self.shape
        nvol     = self.nvol
        nd       = self.nd
        norb     = self.norb
        ns       = self.ns

        # expectation value of Ncond
        gab_r = self.Green.reshape(nvol,nd,nd)

        n = np.sum(np.diagonal(gab_r[0])) * nvol
        self.physics["NCond"] = n.real

        logger.debug("ncond = {}".format(n))

        # expectation value of Sz
        if self.enable_spin_orbital:
            # In spin-orbital basis: index = 2*orb + spin
            # Sz = sum_a (G[2a, 2a] - G[2a+1, 2a+1]) / 2
            g0 = np.diagonal(gab_r[0])  # shape (nd,)
            sz_diag = np.zeros(nd)
            sz_diag[0::2] = 1.0   # up spins
            sz_diag[1::2] = -1.0  # down spins
            sz = np.sum(g0 * sz_diag) * nvol
            self.physics["Sz"] = 0.5 * sz.real
        else:
            gab_r = self.Green
            sigma_z = np.array(np.array([[1,0],[0,-1]]))
            sz = np.sum(np.diagonal(np.einsum('satb,st->ab', gab_r[0], sigma_z))) * nvol
            self.physics["Sz"] = 0.5 * sz.real

        logger.debug("sz = {}".format(sz))

        # residue
        rest = np.linalg.norm(self.Green - self.Green_prev)
        self.physics["Rest"] = rest / np.size(self.Green) * 2

        logger.debug("rest = {}".format(rest))

        # update
        mix = self.param_mod["Mix"]

        g_new = (1.0 - mix) * self.Green_prev + mix * self.Green

        g_new[np.where(abs(g_new) < self.threshold)] = 0.0

        self.Green = g_new

    @do_profile
    def _calc_energy(self):
        logger.debug(">>> _calc_energy")

        nx,ny,nz = self.shape
        nvol     = self.nvol
        nd       = self.nd
        norb     = self.norb
        ns       = self.ns

        energy = {}
        energy_total = 0.0

        ws_list = self._green_list["eigenvalue"]
        vs_list = self._green_list["eigenvector"]
        mus = self._green_list["mu"]
        nblock = len(ws_list)
        block_to_group = self.block_to_group
        group_nconds = self.group_nconds
        ngroup = len(group_nconds)

        if self.T == 0:
            # Compute band energy using global sorting per group
            e_band = 0.0
            for g in range(ngroup):
                # Collect eigenvalues from all blocks in this group
                blist = [b for b in range(nblock) if block_to_group[b] == g]
                all_ev = np.sort(np.concatenate(
                    [ws_list[b].flatten() for b in blist]
                ))
                e_band += np.sum(all_ev[:group_nconds[g]])

            energy["Band"] = e_band
            logger.debug("energy: Band = {}".format(energy["Band"]))
            energy_total += energy["Band"]
        else:
            T = self.T

            def _fermi(t, mu, ev):
                w_ = (ev - mu) / t
                mask_ = w_ < self.ene_cutoff
                w1 = np.where( mask_, w_, 0.0 )
                v1 = 1.0 / (1.0 + np.exp(w1))
                v_ = np.where( mask_, v1, 0.0 )
                return v_

            e_band = 0.0
            for k in range(nblock):
                w = ws_list[k]
                v = vs_list[k]
                mu = mus[block_to_group[k]]

                wt = -(w - mu) / T
                mask_ = wt < self.ene_cutoff
                ln_e = np.where( mask_, np.log1p(np.exp(wt)), wt)

                nn = np.einsum('kal,kl,kal->', np.conjugate(v), _fermi(T, mu, w), v)

                e_band += mu * nn - T * np.sum(ln_e)

            energy["Band"] = e_band.real
            energy_total += energy["Band"]
            logger.debug("energy: Band = {}".format(e_band))

        # In spin-orbital mode, convert Green to virtual form for interaction energy
        if self.enable_spin_orbital:
            norb_inter = self.norb_phys
            gab_r_inter = self._so_to_virtual_green(self.Green)
        else:
            norb_inter = norb
            gab_r_inter = self.Green

        for type in self.inter_table.keys():
            if self.inter_table[type] is not None:
                jab_r = self.inter_table[type].reshape(nvol,norb_inter,norb_inter)
                spin  = self.spin_table[type]
                gab_r = gab_r_inter

                if type == "PairHop":
                    if self.iflag_fock:
                        w1 = np.einsum('stuv, rsavb, rtaub -> rab', spin, gab_r, gab_r)
                        w2 = np.einsum('stuv, rsaub, rtavb -> rab', spin, gab_r, gab_r)
                        ee = np.einsum('rab, rab ->', jab_r, w1-w2)
                    else:
                        w1 = np.einsum('stuv, rsavb, rtaub -> rab', spin, gab_r, gab_r)
                        ee = np.einsum('rab, rab ->', jab_r, w1)
                    energy[type] = -ee/2.0*nvol

                else:
                    if self.iflag_fock:
                        w1 = np.einsum('stuv, vasa, ubtb -> ab', spin, gab_r[0], gab_r[0])
                        w1b = np.broadcast_to(w1, (nvol,norb_inter,norb_inter))
                        w2 = np.einsum('stuv, rubsa, rvatb -> rab', spin, gab_r, gab_r)
                        ee = np.einsum('rab, rab->', jab_r, w1b-w2)
                    else:
                        w1 = np.einsum('stuv, vasa, ubtb -> ab', spin, gab_r[0], gab_r[0])
                        w1b = np.broadcast_to(w1, (nvol,norb_inter,norb_inter))
                        ee = np.einsum('rab, rab->', jab_r, w1b)
                    energy[type] = -ee/2.0*nvol

                energy_total += energy[type].real
                logger.debug("energy: {} = {}".format(type, energy[type]))
            else:
                # logger.info("energy: {} skip".format(type))
                pass

        energy["Total"] = energy_total
        self.physics["Ene"] = energy.copy()

    @do_profile
    def get_results(self):
        return (self.physics, self.Green)

    @do_profile
    def save_results(self, info_outputfile, green_info):
        """Save results to files.
        
        Parameters
        ----------
        info_outputfile : dict
            Output file configuration
        green_info : dict
            Green function information
            
        Saves:
        - Energy and physical quantities
        - Eigenvalues/vectors
        - Green functions
        - One-body quantities
        - RPA data if requested
        """
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

            # Reconstruct full eigenvalue/eigenvector arrays from blocks
            nvol = self.nvol
            nd = self.nd
            blocks = self.block_info
            eg_list = self._green_list["eigenvalue"]
            ev_list = self._green_list["eigenvector"]

            egg = np.zeros((nvol, nd), dtype=np.float64)
            evv = np.zeros((nvol, nd, nd), dtype=np.complex128)

            col_offset = 0
            for i, blk in enumerate(blocks):
                idx = np.array(blk)
                blk_size = len(idx)
                egg[:, col_offset:col_offset + blk_size] = eg_list[i]
                # Place eigenvectors in correct rows
                evv[np.ix_(np.arange(nvol), idx,
                     np.arange(col_offset, col_offset + blk_size))] = ev_list[i]
                col_offset += blk_size

            # # wavevec[k,eigen_index] = \vec(k)
            # wv = self.wave_table
            # wvs = wv.shape
            # wvv = np.transpose(np.broadcast_to(wv, ((egs[0]*egs[2]),wvs[0],wvs[1])), (1,0,2))

            file_name = os.path.join(path_to_output, info_outputfile["eigen"])
            np.savez(file_name,
                     eigenvalue  = egg,
                     eigenvector = evv,
                     wavevector_unit = self.kvec,
                     wavevector_index = self.wavenum_table,
                     )
            logger.info("save_results: save eigenvalues and eigenvectors in file {}".format(file_name))

        # export green function to "green.npz",
        # or to a file of the name specified by 'green' keyword.
        # if the keyword has empty string "", nothing exported.
        fname = info_outputfile.get("green", "green")
        if fname != "":
            file_name = os.path.join(path_to_output, fname)
            self._save_green(file_name)

        if "onebodyg" in info_outputfile.keys():
            file_name = os.path.join(path_to_output, info_outputfile["onebodyg"])
            self._save_greenone(file_name, green_info)

        if "initial" in info_outputfile.keys():
            logger.warning("save_results: save initial is not supported")
            pass

        if "rpa" in info_outputfile.keys():
            file_name = os.path.join(path_to_output, info_outputfile["rpa"])
            self._save_trans_mod(file_name)

        if "export_hamiltonian" in info_outputfile.keys():
            self._export_hamiltonian(path_to_output, info_outputfile["export_hamiltonian"])

    @do_profile
    def _read_green(self, file_name):
        try:
            v = np.load(file_name)
        except FileNotFoundError:
            logger.info("_read_green: file {} not found".format(file_name))
            return None
        return self._read_green_from_data(v)

    def _read_green_from_data(self, ginfo):
        if self.has_sublattice:
            if "green_sublattice" in ginfo.files:
                logger.debug("_read_green: read green_sublattice")
                data = ginfo["green_sublattice"]
            elif "green" in ginfo.files:
                logger.debug("_read_green: read green and reshape")
                data = self._reshape_green(ginfo["green"])
            else:
                logger.debug("_read_green: no data found")
                data = None
        else:
            if "green" in ginfo.files:
                logger.debug("_read_green: read green")
                data = ginfo["green"]
            else:
                logger.debug("_read_green: no data found")
                data = None
        return data

    @do_profile
    def _save_green(self, file_name):
        if self.has_sublattice:
            green_orig = self._deflate_green(self.Green)
            np.savez(file_name, green = green_orig, green_sublattice = self.Green)
        else:
            np.savez(file_name, green = self.Green)
        logger.info("save_results: save green function to file {}".format(file_name))

    @do_profile
    def _save_greenone(self, file_name, green_info):
        if not "onebodyg_uhf" in green_info or not "geometry_uhf" in green_info:
            logger.error("_save_greenone: onebodyg_uhf and geometry_uhf are required")
            return None

        if self.has_sublattice:
            gr = self._deflate_green(self.Green)
        else:
            gr = self.Green

        geom = green_info["geometry_uhf"]
        greenone = green_info["onebodyg_uhf"]

        lx,ly,lz = self.cellshape
        lvol     = self.cellvol
        norb     = self.norb
        nd       = self.nd
        ns       = self.ns

        tbl = geom["site2vec"]

        def _pack_site(ix,iy,iz):
            return ix + lx * (iy + ly * iz)

        with open(file_name, "w") as fw:
            # G(i,a,s; j,b,t) := G(a,s;b,t;j-i)
            for (i,s,j,t) in greenone:
                (ix,iy,iz,a) = tbl[i]
                (jx,jy,jz,b) = tbl[j]

                kx = (jx - ix + lx) % lx
                ky = (jy - iy + ly) % ly
                kz = (jz - iz + lz) % lz
                idx = _pack_site(kx,ky,kz)

                v = gr[idx,s,a,t,b]

                fw.write("{:3} {:3} {:3} {:3}  {:.12e} {:.12e}\n".format(
                    i, s, j, t, v.real, v.imag
                ))

        logger.info("save_results: save greenone to file {}".format(file_name))

    def _save_trans_mod(self, file_name):
        nx,ny,nz = self.shape
        nvol = self.nvol
        nd = self.nd

        tab_k = self.ham
        tab_r = np.fft.fftn(tab_k.reshape(nx,ny,nz,nd,nd), axes=(0,1,2)).reshape(nvol,nd,nd) / nvol

        if self.has_sublattice:
            # use deflate_green to convert to original lattice

            norb = self.norb
            ns = self.ns

            tab_r_defl = self._deflate_green(tab_r.reshape(nvol,ns,norb,ns,norb))

            lvol = self.cellvol
            norb_orig = self.norb_orig
            nd0 = norb_orig * ns

            np.savez(file_name, trans_mod=tab_r_defl.reshape(lvol,nd0,nd0), trans_mod_sublattice=tab_r)
        else:
            np.savez(file_name, trans_mod=tab_r)
        logger.info("save_results: save trans_mod to file {}".format(file_name))

    def _export_geometry(self, file_name):
        geom = self.param_ham["Geometry"]
        with open(file_name, "w") as fw:
            # write unit vector
            rvec = geom["rvec"]
            for k in range(3):
                fw.write("{:.8f} {:.8f} {:.8f}\n".format(rvec[k][0], rvec[k][1], rvec[k][2]))
            # write number of orbitals
            fw.write("{}\n".format(geom["norb"]))
            # write center coordinates in fractional coordinates
            center = geom["center"]
            for k in range(len(center)):
                fw.write("{:.8f} {:.8f} {:.8f}\n".format(center[k][0], center[k][1], center[k][2]))

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
                        *irvec, *orbvec, v.real, v.imag
                    ))

    def _export_hamiltonian(self, path_to_output, prefix):
        for type in self.param_ham.keys():
            if type in ["Geometry"]:
                file_name = os.path.join(path_to_output, prefix + type + ".dat")
                self._export_geometry(file_name)
                logger.info("export Geometry to {}".format(file_name))
            elif type in ["Initial"]:
                pass
            else:
                file_name = os.path.join(path_to_output, prefix + type + ".dat")
                self._export_interaction(type, file_name)
                logger.info("export {} to {}".format(type, file_name))
