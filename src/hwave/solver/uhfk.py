import logging
import numpy as np
import os
from .base import solver_base
from .perf import do_profile
from . import fold
from .kgrid import reverse_fft_axes
from . import hartree_fock
from . import spin_orbital
from . import hamiltonian_k
from . import uhfk_blocks
from . import uhfk_fill
from ..qlmsio import wan90

logger = logging.getLogger("qlms").getChild("uhfk")

class UHFk(solver_base):
    @do_profile
    def __init__(self, param_ham, info_log, info_mode, param_mod=None):
        """Initialize the UHF solver.
        
        Parameters
        ----------
        param_ham : dict
            Hamiltonian parameters. The solver normalises tables IN this
            container: a sublattice fold rewrites them, and in normal mode
            the spin block of a spin-orbital-format Transfer (orbital
            indices above the Geometry norb) is removed from it with one
            warning, also without a sublattice. A container reused for a
            second solver therefore no longer holds the original Transfer;
            build each solver from its own ``QLMSkInput`` result.
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
        # A spin-resolved (9-column) Transfer (issue #323) is consumed in the
        # normal (ns = 2) mode: the reader folds it to spin-block keys
        # s*norb + a, which is the spin-major index of this solver's arrays,
        # so norb, geom.dat and the two-body files keep their physical
        # indices and the 2Sz structure check applies unchanged. It
        # contradicts spin-orbital mode (Geometry norb = 2 * physical).
        self.spin_resolved_transfer = (
            spin_orbital.spin_resolved_transfer_format(param_ham) is not None)
        spin_orbital.check_spin_resolved_clash(param_ham, self.enable_spin_orbital)
        # the reader folded the rows with the line-2 orbital count; it must
        # be the Geometry norb this solver decodes the keys with
        spin_orbital.check_spin_resolved_norb(
            param_ham, param_ham["Geometry"]["norb"], "UHFk")
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
            "FreeEne": { "Total": 0.0, "Band": 0.0 },
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
        self.enable_spin_orbital = spin_orbital.resolve_flag(param, allow_float=True)

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

        # - GPU execution is implemented only for RPA/FLEX and the dynamic
        #   Eliashberg solver; warn instead of silently ignoring the request.
        from hwave.solver.backend import as_bool as _as_bool
        if _as_bool(self.param_mod.get("gpu", False)):
            logger.warning(
                "gpu=true is not supported by the UHF solvers and is "
                "ignored; GPU execution is available for RPA, FLEX, and "
                "the dynamic Eliashberg solver.")

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
            if self.enable_spin_orbital:
                # Sz is not a good quantum number when spin is folded into
                # the orbital index; block detection would silently assign
                # the total Ncond to the spin-mixed blocks and ignore the
                # constraint. Fail fast instead.
                logger.error(
                    "2Sz-fixed mode is not supported with "
                    "enable_spin_orbital; omit 2Sz (Sz-free mode).")
                raise ValueError(
                    "2Sz-fixed mode is not supported with "
                    "enable_spin_orbital; use Sz-free mode.")
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
            self.norb_phys_orig = spin_orbital.physical_norb(
                norb, True, source="Geometry")
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
        logger.info("    spin-resolved Transfer = {}".format(self.spin_resolved_transfer))

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
            Reshaped geometry (see fold.reshape_geometry)
        """
        return fold.reshape_geometry(geom, self.subshape)

    @do_profile
    def _reshape_interaction(self, ham, enable_spin_orbital,
                             drop_spin_block=False):
        """Reshape interaction terms for sublattice.

        Parameters
        ----------
        ham : dict
            Original interaction terms
        enable_spin_orbital : bool
            Whether THIS table uses spin-orbital indices (Transfer in SO mode)
        drop_spin_block : bool
            True for one-body tables only (see fold.reshape_interaction)

        Returns
        -------
        dict
            Reshaped interaction terms (see fold.reshape_interaction)
        """
        return fold.reshape_interaction(
            ham, self.subshape, self.shape,
            norb_so_orig=self.param_ham_orig["Geometry"]["norb"],
            norb_phys_orig=self.norb_phys_orig,
            enable_spin_orbital=enable_spin_orbital,
            drop_spin_block=drop_spin_block)

    @do_profile
    def _reshape_green(self, green):
        """Convert Green function to sublattice basis.

        ``fold.fold_layout`` with the Green convention (the within-cell
        offset on the first orbital slot, ``dr = ri - rj``), the inverse of
        the default ``_deflate_green``.
        
        Parameters
        ----------
        green : ndarray
            Green function in original basis
        
        Returns
        -------
        ndarray
            Green function in sublattice basis

        Raises
        ------
        ValueError
            If ``green`` is not ``(Lvol, ns, norb_orig, ns, norb_orig)``.
        """
        return fold.fold_layout(
            green, cellshape=self.cellshape, subshape=self.subshape,
            norb_orig=self.norb_orig, ns=self.ns)

    @do_profile
    def _deflate_green(self, green_sub, hamiltonian=False):
        """Convert Green function back to original basis.

        Parameters
        ----------
        green_sub : ndarray
            Green function in sublattice basis
        hamiltonian : bool, optional
            Select the cross-sublattice slot convention. With
            ``green[jsite,s,a,t,b] = green_sub[isite,s,aa,t,bb]``, where the
            original site ``jsite`` decomposes per axis into a supercell site
            ``isite`` and the within-cell offset ``ir``, the default (Green)
            carries ``ir`` on the FIRST orbital slot
            (``aa = a + norb_orig*ir``, ``bb = b``); ``hamiltonian=True`` carries
            it on the SECOND slot (``aa = a``, ``bb = b + norb_orig*ir``). The
            two differ because the folded Green ``self.Green`` from ``_green()``
            is the orbital transpose of the Hamiltonian-convention object (it
            stores ``G_ab(k) = sum_l conj(V_al) f_l V_bl``). Use the default for
            a Green function and ``hamiltonian=True`` for Hamiltonian / transfer
            quantities (e.g. ``trans_mod``).
            The layout is ``fold.unfold_layout``.

        Returns
        -------
        ndarray
            Green function in original basis
        """
        return fold.unfold_layout(
            green_sub, cellshape=self.cellshape, subshape=self.subshape,
            norb_orig=self.norb_orig, ns=self.ns, hamiltonian=hamiltonian)

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
                    # The spin block of a normal-mode table was dropped at
                    # the boundary (_drop_spin_block_transfer), so the fold
                    # raises on any leftover instead of dropping. A
                    # spin-resolved table folds with the spin stride like a
                    # spin-orbital one: s*norb + a -> s*(norb*subvol) + a'.
                    tbl = self._reshape_interaction(
                        self.param_ham[type],
                        self.enable_spin_orbital or self.spin_resolved_transfer)
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
        self._drop_spin_block_transfer()
        self._check_hermite()
        self._check_spin_orbital_compatibility()

    def _drop_spin_block_transfer(self):
        """Normalize a normal-mode Transfer once, at the boundary.

        A spin-orbital-format file read with enable_spin_orbital = False
        carries its spin block as orbital indices >= norb; the documented
        behaviour is to ignore those entries with a warning. They are
        dropped here, right after the index check, so every later consumer
        (the hermiticity check, the sublattice fold, the H(k) builder)
        sees an in-range table and raises on anything else instead of
        guarding for them individually.

        The replacement is made in the ``param_ham`` container the caller
        passed in (the same container the sublattice fold rewrites), so a
        second solver built from the same ``QLMSkInput`` result receives
        the normalised table and no further warning. Build each solver
        from its own reader result when the original table matters.
        """
        if self.enable_spin_orbital or self.spin_resolved_transfer:
            return
        if "Transfer" not in self.param_ham:
            return
        table, dropped = hamiltonian_k.drop_spin_block(
            self.param_ham["Transfer"], self.norb_orig, label="Transfer")
        if dropped:
            self.param_ham["Transfer"] = table

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
                # In spin-orbital mode the Geometry norb already includes spin
                # (norb = nd), so valid Transfer indices are [0, norb). In normal
                # mode the file may carry indices up to 2*norb: the spin block
                # of a spin-orbital-format file (dropped at the boundary by
                # _drop_spin_block_transfer), or the spin-block keys s*norb + a
                # of a spin-resolved Transfer (all of them consumed).
                trans_max = norb if self.enable_spin_orbital else norb * 2
                fail = 0
                for (irvec,orbvec), v in self.param_ham[k].items():
                    if not all([ 0 <= orbvec[i] < trans_max for i in range(2) ]):
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
            # this check runs before the sublattice reshape, so Transfer still
            # carries the original-geometry orbital indices
            norb = self.norb_orig
            # a spin-resolved table carries both spins (s*norb + a); the
            # spin block of a normal-mode table was dropped at the boundary,
            # so an out-of-range entry here is a defect and raises
            dim = 2 * norb if self.spin_resolved_transfer else norb
            tab_r = hamiltonian_k.dense_table(
                self.param_ham["Transfer"], (nx, ny, nz), dim)

            t = np.conjugate(
                    np.transpose(
                        reverse_fft_axes(tab_r, (0, 1, 2)),
                        (0,1,2,4,3)
                    )
                )

            if not np.allclose(t, tab_r, atol=self.hermite_tolerance, rtol=0.0):
                msg = "Hermiticity check failed: |T_ba(-r)^* - T_ab(r)| = {}".format(np.sum(np.abs(t - tab_r)))
                if self.strict_hermite:
                    logger.error(msg)
                    errlist = np.array(np.nonzero(~np.isclose(t, tab_r, atol=self.hermite_tolerance, rtol=0.0))).transpose()
                    for idx in errlist[:max_print]:
                        note = ""
                        if self.spin_resolved_transfer:
                            # decode the spin-block composite indices
                            (a, s), (b, t_) = (
                                spin_orbital.spinblock_split(int(i), norb)
                                for i in idx[3:5])
                            note = " [(orbital, spin) = ({}, {}), ({}, {}); 1-based orbitals]".format(
                                a + 1, s, b + 1, t_)
                        logger.error("    index={}{}, T_ab={}, T_ba^*={}".format(
                            idx, note, tab_r[tuple(idx)], t[tuple(idx)]))
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

        # the convergence record (scf_status.dat), defined before the loop
        # so that a run of zero steps has one too
        self._reset_scf_status()

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

            if self._record_scf_step(i_step, self.physics["Rest"], self.param_mod["eps"]):
                break

        if self.scf_converged:
            logger.info("UHFk calculation succeeded: rest={}, eps={}."
                        .format(self.physics["Rest"], self.param_mod["eps"]))
            logger.info("Total Energy = {}.".format(self.physics["Ene"]["Total"]))
        elif self.scf_iterations > 0:
            # (after zero steps save_results logs the one warning there is)
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
            try:
                _data = self._read_green_from_data(green_info["initial"])
            except ValueError as exc:
                # a "green" array of another lattice/orbital layout (e.g. a
                # file from a run with another CellShape) cannot be folded
                logger.error("initial green function: {}".format(exc))
                exit(1)
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
        if self.has_sublattice:
            logger.error("initial_uhf requires SubShape = [1, 1, 1]; give a sublattice-ordered initial state through "
                         "initial (green.npz) instead")
            exit(1)

        lx,ly,lz = self.cellshape
        lvol     = self.cellvol
        norb     = self.norb
        nd       = self.nd
        ns       = self.ns

        def _pack_site(ix,iy,iz):
            return iz + lz * (iy + ly * ix)

        tbl = {}
        for idx, (ix,iy,iz,a) in geom["site2vec"].items():
            site = _pack_site(ix,iy,iz)
            tbl[idx] = (site, a)

        green_uhf = np.zeros((lvol,lvol,ns,norb,ns,norb),dtype=np.complex128)
        for (i,s,j,t),v in info.items():
            isite, a = tbl[i]
            _, b = tbl[j]
            ix,iy,iz,_ = geom["site2vec"][i]
            jx,jy,jz,_ = geom["site2vec"][j]
            # G(d) = <c^dag(R+d)c(R)>, with the same sign as _save_greenone.
            jsite = _pack_site((ix-jx) % lx, (iy-jy) % ly, (iz-jz) % lz)
            green_uhf[isite,jsite,s,a,t,b] = v

        # Primitive-cell translation invariance makes the origin row sufficient
        # for this initial guess. Sublattice-ordered states need green_sublattice
        # input instead: a single origin row cannot preserve their site dependence.
        return green_uhf[0]

    @do_profile
    def _make_ham_trans(self):
        logger.debug(">>> _make_ham_trans")

        table = self.param_ham["Transfer"]
        if self.enable_spin_orbital or self.spin_resolved_transfer:
            # the table's indices already carry the spin: interleaved
            # 2*orb + spin (kept in this solver's order) or the spin-block
            # s*norb + a of a spin-resolved Transfer -- both (nd, nd)
            ham = hamiltonian_k.one_body_k(table, self.shape, self.nd)
        else:
            # physical orbitals, spin-independent: embed over the spin axis.
            # The spin block of a spin-orbital-format file was dropped at
            # the boundary (_drop_spin_block_transfer), so anything out of
            # range here raises.
            ham = hamiltonian_k.one_body_k(table, self.shape, self.norb,
                                           spin_embed=True)

        # store
        self.ham_trans = ham

    @do_profile
    def _make_ham_inter(self):
        """Build the per-type interaction and spin tables.

        The normal (spin-major) mode delegates to the shared kernel
        ``hartree_fock.build_interaction_tables`` (#181 Phase B); the
        ``discarded`` report it returns (declared entries the per-type
        semantics drop, e.g. a CoulombIntra entry that is not r = 0,
        a == b) is ignored here exactly as those entries were silently
        skipped before, and the ambiguous aggregate/explicit Coulomb case
        keeps its legacy logger.error + exit(1).
        """
        logger.debug(">>> _make_ham_inter")
        norb = self.norb_phys if self.enable_spin_orbital else self.norb
        try:
            tables = hartree_fock.build_interaction_tables(
                self.param_ham, norb, self.shape)
        except ValueError as exc:
            logger.error(str(exc))
            exit(1)
        self.inter_table = tables.inter_table
        self.spin_table = tables.spin_table
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
        The detection itself is :func:`hwave.solver.uhfk_blocks.detect_blocks`.
        """
        logger.debug(">>> _detect_blocks")

        nd = self.nd
        found = uhfk_blocks.detect_blocks(
            self.ham_trans, self.inter_table, self.spin_table,
            nd=nd, norb=self.norb, ns=self.ns, iflag_fock=self.iflag_fock,
            sz_free=self.sz_free, enable_spin_orbital=self.enable_spin_orbital,
            nconds=tuple(self.Nconds))
        blocks = list(found.block_info)
        group_nconds = list(found.group_nconds)

        self.block_info = blocks
        self.block_to_group = list(found.block_to_group)
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
        """Spin-orbital Green function ``(nvol, 1, nd, 1, nd)`` (interleaved
        index ``2*orb + spin``) -> the virtual spin-separate form
        ``(nvol, 2, norb_phys, 2, norb_phys)``: the spin-block reorder of
        :mod:`hwave.solver.spin_orbital`, so the normal-mode interaction
        contractions apply unchanged."""
        nvol = G_so.shape[0]
        norb_phys = self.nd // 2
        G_sb = spin_orbital.interleaved_to_spinblock(
            G_so.reshape(nvol, self.nd, self.nd), norb_phys)
        return np.ascontiguousarray(G_sb).reshape(nvol, 2, norb_phys, 2, norb_phys)

    def _virtual_ham_to_so(self, H_virt):
        """Inverse of :meth:`_so_to_virtual_green` for a Hamiltonian
        ``(nvol, 2, norb_phys, 2, norb_phys)`` -> ``(nvol, nd, nd)`` in the
        interleaved order."""
        nvol = H_virt.shape[0]
        return np.ascontiguousarray(spin_orbital.spinblock_to_interleaved(
            H_virt.reshape(nvol, self.nd, self.nd), self.nd // 2))

    @do_profile
    def _make_ham(self):
        logger.debug(">>> _make_ham")

        nvol     = self.nvol
        nd       = self.nd

        # hamiltonian H_{ab,st}(k) : ham(k,(s,a),(t,b))
        ham = np.zeros((nvol,nd,nd), dtype=np.complex128)

        # transfer term  T_{ab}(k) (note convention)
        logger.debug("Transfer")
        ham += self.ham_trans

        # interaction term and the block mask, in place
        # (uhfk_blocks.accumulate_hf_projected: in spin-orbital mode the
        # interactions run in the virtual spin-block frame, #83)
        uhfk_blocks.accumulate_hf_projected(
            ham, self.Green, self._hf_tables(), include_fock=self.iflag_fock,
            blocks=self.block_info, debug=logger.debug)

        # store
        self.ham = ham

    def _hf_tables(self):
        """The interaction tables of this solver as a
        :class:`hwave.solver.uhfk_blocks.HFTables`."""
        return uhfk_blocks.HFTables(
            inter_table=self.inter_table, spin_table=self.spin_table,
            shape=tuple(self.shape), enable_spin_orbital=self.enable_spin_orbital,
            norb_phys=self.norb_phys)

    @do_profile
    def _diag(self):
        logger.debug(">>> _diag")

        # Diagonalize each block independently (uhfk_fill.diag_blocks);
        # stored as lists of per-block arrays
        eigenvalues, eigenvectors = uhfk_fill.diag_blocks(self.ham, self.block_info)
        self._green_list["eigenvalue"] = eigenvalues
        self._green_list["eigenvector"] = eigenvectors

    @do_profile
    def _green(self):
        logger.debug(">>> _green")

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
        group_blocks = uhfk_fill.group_blocks(block_to_group, ngroup)

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

        # G(k) assembled from the blocks, then G(r) (uhfk_fill.assemble_green)
        green = uhfk_fill.assemble_green(vs_list, group_dists, blocks, shape=self.shape,
                                         ns=self.ns, norb=self.norb)

        # store
        self.Green_prev = self.Green
        self.Green = green

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
        """``uhfk_fill.find_dist_group_zero_t`` on this solver's k-grid."""
        return uhfk_fill.find_dist_group_zero_t(ws_list, vs_list, ncond,
                                                shape=self.shape)

    def _find_dist_group_nonzero_t(self, ws_list, vs_list, ncond):
        """``uhfk_fill.find_dist_group_nonzero_t`` at this solver's T."""
        try:
            return uhfk_fill.find_dist_group_nonzero_t(
                ws_list, vs_list, ncond, T=self.T, ene_cutoff=self.ene_cutoff)
        except uhfk_fill.MuNotConvergedError:
            exit(1)

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

        # expectation value of the spin vector (Sx, Sy, Sz)
        #   S_alpha = (1/2) sum_a sum_{st} (sigma_alpha)_{st} <c^dag_{as} c_{at}>,
        # read in the spin-block frame m[s, a, t, b] of the on-site Green
        # function G_{ij} = <c^dag_i c_j> (issue #350). The mean-field state
        # can point along any spin axis, so Sx and Sy can be nonzero also
        # without spin-orbit coupling.
        if self.enable_spin_orbital:
            # the virtual spin-block form of the interleaved index
            m = self._so_to_virtual_green(self.Green)[0]
            g0 = np.diagonal(gab_r[0])           # rho_{ii}, interleaved
            sz_diag = np.where(spin_orbital.spin_of(np.arange(nd)) == 0, 1.0, -1.0)
            sz = np.sum(g0 * sz_diag) * nvol
        else:
            # normal mode: G[0] is already in the (s, a, t, b) frame
            m = self.Green[0]
            sigma_z = np.array(np.array([[1,0],[0,-1]]))
            sz = np.sum(np.diagonal(np.einsum('satb,st->ab', m, sigma_z))) * nvol
        # Sx = (1/2)(ud + du), Sy = -(i/2)(ud - du)
        ud = np.diagonal(m[0, :, 1, :])          # <c^dag_up c_down> = <S^+>
        du = np.diagonal(m[1, :, 0, :])          # <c^dag_down c_up> = <S^->
        sx = np.sum(ud + du) * nvol
        sy = np.sum(-1j * (ud - du)) * nvol
        self.physics["Sz"] = 0.5 * sz.real
        self.physics["Sx"] = 0.5 * sx.real
        self.physics["Sy"] = 0.5 * sy.real

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
        free_energy = {}
        energy_total = 0.0
        free_total = 0.0

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

            # At T=0 the entropy vanishes, so internal energy and free energy
            # coincide.
            e_band_internal = e_band
            e_band_free = e_band
            logger.debug("energy: Band = {}".format(e_band))
        else:
            T = self.T

            def _fermi(t, mu, ev):
                w_ = (ev - mu) / t
                mask_ = w_ < self.ene_cutoff
                w1 = np.where( mask_, w_, 0.0 )
                v1 = 1.0 / (1.0 + np.exp(w1))
                v_ = np.where( mask_, v1, 0.0 )
                return v_

            # Internal band energy  E_band = sum_n eps_n f(eps_n)
            # Free   band energy    F_band = mu*N + Omega_0
            #                              = mu*N - T*sum_n ln(1+exp(-(eps_n-mu)/T))
            e_band_internal = 0.0
            e_band_free = 0.0
            for k in range(nblock):
                w = ws_list[k]
                v = vs_list[k]
                mu = mus[block_to_group[k]]

                fermi = _fermi(T, mu, w)

                wt = -(w - mu) / T
                mask_ = wt < self.ene_cutoff
                ln_e = np.where( mask_, np.log1p(np.exp(wt)), wt)

                nn = np.einsum('kal,kl,kal->', np.conjugate(v), fermi, v)

                e_band_internal += np.sum(w * fermi)
                e_band_free += mu * nn - T * np.sum(ln_e)

            e_band_internal = np.real(e_band_internal)
            e_band_free = np.real(e_band_free)
            logger.debug("energy: Band (internal) = {}".format(e_band_internal))
            logger.debug("free energy: Band = {}".format(e_band_free))

        energy["Band"] = e_band_internal
        free_energy["Band"] = e_band_free
        energy_total += e_band_internal
        free_total += e_band_free

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
                        w1 = np.einsum('stuv, rvbsa, rubta -> rab', spin, np.conjugate(gab_r), np.conjugate(gab_r), optimize=True)
                        w2 = np.einsum('stuv, rubsa, rvbta -> rab', spin, np.conjugate(gab_r), np.conjugate(gab_r), optimize=True)
                        ee = np.einsum('rab, rab ->', jab_r, w1-w2)
                    else:
                        w1 = np.einsum('stuv, rvbsa, rubta -> rab', spin, np.conjugate(gab_r), np.conjugate(gab_r), optimize=True)
                        ee = np.einsum('rab, rab ->', jab_r, w1)
                    energy[type] = -ee/2.0*nvol

                else:
                    if self.iflag_fock:
                        w1 = np.einsum('stuv, vasa, ubtb -> ab', spin, gab_r[0], gab_r[0])
                        w1b = np.broadcast_to(w1, (nvol,norb_inter,norb_inter))
                        w2 = np.einsum('stuv, rubsa, rtbva -> rab', spin, np.conjugate(gab_r), gab_r, optimize=True)
                        ee = np.einsum('rab, rab->', jab_r, w1b-w2)
                    else:
                        w1 = np.einsum('stuv, vasa, ubtb -> ab', spin, gab_r[0], gab_r[0])
                        w1b = np.broadcast_to(w1, (nvol,norb_inter,norb_inter))
                        ee = np.einsum('rab, rab->', jab_r, w1b)
                    energy[type] = -ee/2.0*nvol

                # Interaction (double-counting corrected) energy is identical
                # for the internal energy and the free energy.
                energy_total += energy[type].real
                free_total += energy[type].real
                logger.debug("energy: {} = {}".format(type, energy[type]))
            else:
                # logger.info("energy: {} skip".format(type))
                pass

        energy["Total"] = energy_total
        free_energy["Total"] = free_total
        self.physics["Ene"] = energy.copy()
        self.physics["FreeEne"] = free_energy.copy()

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

        # refuse a status file that would overwrite another output, then
        # write the convergence status first; after zero SCF steps it is the
        # only result there is
        self.check_output_paths(info_outputfile)
        if not self._save_scf_status(info_outputfile, "UHFk",
                                     "two_frobenius_over_size:uhfk_all_r"):
            return

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

                # Free energy F = E - T*S (Helmholtz). At T=0 it equals the
                # internal energy above; at finite T it is the quantity the
                # self-consistent loop minimizes.
                if "FreeEne" in self.physics:
                    for type in type_ex:
                        fw.write("FreeEnergy_{:<8} = {}\n".format(
                            type, self.physics["FreeEne"][type].real))

                fw.write("NCond   = {}\n".format(self.physics["NCond"]))
                fw.write("Sz      = {}\n".format(self.physics["Sz"]))
                # Chemical potential (Fermi level). One value per mu-group; at
                # T=0 it is the midpoint of the HOMO-LUMO gap.
                mu_arr = self._green_list.get("mu", None)
                if mu_arr is not None:
                    mu_arr = np.atleast_1d(mu_arr)
                    if mu_arr.size == 1:
                        fw.write("ChemicalPotential = {}\n".format(float(mu_arr[0])))
                    else:
                        for g, mu in enumerate(mu_arr):
                            fw.write("ChemicalPotential_{} = {}\n".format(g, float(mu)))
                # Transverse spin, in both modes, as the last two lines: the
                # mean-field state can point along any spin axis, so Sx and
                # Sy can be nonzero also without spin-orbit coupling.
                fw.write("Sx      = {}\n".format(self.physics["Sx"]))
                fw.write("Sy      = {}\n".format(self.physics["Sy"]))

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
                     momentum_convention = "e_plus_ikR",
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

    def _output_artifacts(self, info_outputfile):
        """The files :meth:`save_results` writes besides the status file."""
        from hwave.solver.base import npz_path
        path = info_outputfile.get("path_to_output", "output")
        out = []

        def name(key, default=None):
            value = info_outputfile.get(key, default)
            return value if isinstance(value, str) else None

        if name("energy") is not None:
            out.append(("energy", os.path.join(path, name("energy"))))
        if name("eigen") is not None:
            out.append(("eigen", npz_path(os.path.join(path, name("eigen")))))
        green = name("green", "green")
        if green:
            out.append(("green", npz_path(os.path.join(path, green))))
        if name("onebodyg") is not None:
            out.append(("onebodyg", os.path.join(path, name("onebodyg"))))
        if name("rpa") is not None:
            out.append(("rpa", npz_path(os.path.join(path, name("rpa")))))
        prefix = name("export_hamiltonian")
        if prefix is not None:
            for type in self.param_ham.keys():
                if type != "Initial":
                    out.append(("export_hamiltonian",
                                os.path.join(path, prefix + type + ".dat")))
        return out

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
            # Tag the deflated "green" key with the fold convention used since
            # PR #35: the within-cell offset sits on the FIRST orbital slot (the
            # inverse of the Hamiltonian/transfer deflate). Older files lacked
            # this marker and stored "green" with the Hamiltonian convention, so
            # a reader (e.g. RPA green_init) needs the marker to fold "green"
            # back without silently using the wrong sign. See issue #36.
            np.savez(file_name, green = green_orig,
                     green_sublattice = self.Green,
                     green_convention = np.array("green_slot_first"),
                     # issue #133: UHFk always used the documented e^{+ikR}
                     # sign; the stamp makes that machine-readable
                     momentum_convention = "e_plus_ikR")
        else:
            np.savez(file_name, green = self.Green,
                     momentum_convention = "e_plus_ikR")
        logger.info("save_results: save green function to file {}".format(file_name))

    @do_profile
    def _save_greenone(self, file_name, green_info):
        if not "onebodyg_uhf" in green_info or not "geometry_uhf" in green_info:
            logger.error("_save_greenone: onebodyg_uhf and geometry_uhf are required")
            return None

        if self.enable_spin_orbital:
            # In spin-orbital mode spin is folded into the orbital index (ns=1),
            # so the (i, s, j, t) one-body Green output cannot be expressed on
            # the (spin, orbital, spin, orbital) axes. Reject rather than crash.
            logger.error(
                "_save_greenone: one-body Green output (onebodyg) is not "
                "supported in spin-orbital mode")
            return None

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
            return iz + lz * (iy + ly * ix)

        with open(file_name, "w") as fw:
            # Output <c^dag(i,a,s) c(j,b,t)> for each requested pair.
            for (i,s,j,t) in greenone:
                (ix,iy,iz,a) = tbl[i]
                (jx,jy,jz,b) = tbl[j]

                if self.has_sublattice:
                    bx,by,bz = self.subshape
                    nx,ny,nz = self.shape
                    a += self.norb_orig * (ix % bx + bx * (iy % by + by * (iz % bz)))
                    b += self.norb_orig * (jx % bx + bx * (jy % by + by * (jz % bz)))

                    # _green uses exp(-ik.r), so the creation site supplies
                    # the positive displacement: G(r) = <c^dag(R+r)c(R)>.
                    kx = (ix // bx - jx // bx) % nx
                    ky = (iy // by - jy // by) % ny
                    kz = (iz // bz - jz // bz) % nz
                    # Match _green's (nx, ny, nz, ...) C-order reshape.
                    idx = kz + nz * (ky + ny * kx)
                else:
                    # As in the folded branch, d = R_i - R_j.
                    kx = (ix - jx) % lx
                    ky = (iy - jy) % ly
                    kz = (iz - jz) % lz
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

            tab_r_defl = self._deflate_green(tab_r.reshape(nvol,ns,norb,ns,norb), hamiltonian=True)

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
            # write index and elements. Orbital indices are stored 0-based in
            # memory but the wannier90-like format is 1-based (read_w90 subtracts
            # one), so emit orbvec + 1 to keep the file round-trippable.
            for (irvec,orbvec), v in self.param_ham[type].items():
                if (abs(v) > 1.0e-12):
                    fw.write("{:3} {:3} {:3} {:3} {:3}  {:.12f} {:.12f}\n".format(
                        *irvec, orbvec[0] + 1, orbvec[1] + 1, v.real, v.imag
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
