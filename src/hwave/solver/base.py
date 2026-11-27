"""Base solver class for Hartree-Fock calculations.

This module provides the base solver class that implements common functionality
for Hartree-Fock calculations, including parameter handling and validation.

"""
import math
import os
import sys
from collections.abc import ItemsView, ValuesView
from requests.structures import CaseInsensitiveDict
import logging
logger = logging.getLogger("qlms").getChild("solver")


class ParamDict(CaseInsensitiveDict):
    """The solver's ``[mode.param]`` table, recording which keys it read.

    Every read through ``[]``, ``get`` or ``in`` adds the lower-cased key to
    :attr:`accessed`. Iterating the table (``items()``, ``values()``, e.g.
    to log every parameter) does not count as a read. A ``ParamDict`` built
    from another one, or returned by :meth:`copy`, shares its record, so
    every copy a solver makes of the table reports into the same set
    (``copy.deepcopy`` gives an independent record). ``hwave`` uses the
    record to report the keys of the input file that the run never read
    (GitHub issue #236).
    """

    def __init__(self, data=None, **kwargs):
        if isinstance(data, ParamDict):
            # a copy is not a read: take the entries without touching them
            self.accessed = data.accessed
            super().__init__(dict(data._store.values()), **kwargs)
        else:
            self.accessed = set()
            super().__init__(data, **kwargs)

    def __getitem__(self, key):
        self.accessed.add(key.lower())
        return super().__getitem__(key)

    def __contains__(self, key):
        self.accessed.add(key.lower())
        return key.lower() in self._store

    def items(self):
        return _UntrackedItems(self)

    def values(self):
        return _UntrackedValues(self)

    def copy(self):
        return ParamDict(self)


class _UntrackedItems(ItemsView):
    """``items()`` view of a :class:`ParamDict` that does not record reads."""

    def __iter__(self):
        for key, value in self._mapping._store.values():
            yield (key, value)

    def __contains__(self, item):
        key, value = item
        entry = self._mapping._store.get(key.lower())
        return entry is not None and (entry[1] is value or entry[1] == value)


class _UntrackedValues(ValuesView):
    """``values()`` view of a :class:`ParamDict` that does not record reads."""

    def __iter__(self):
        for _, value in self._mapping._store.values():
            yield value


def scf_status_file_name(info_outputfile):
    """The ``[file.output] scf_status`` file name of the UHF solvers:
    ``"scf_status.dat"`` when the key is absent, ``""`` when disabled.
    The key is matched exactly, as the other ``[file.output]`` keys are.
    Raises ValueError when the value is not a string."""
    value = info_outputfile.get("scf_status", "scf_status.dat")
    if not isinstance(value, str):
        raise ValueError(
            "[file.output] scf_status must be a file name (a string), or \"\" "
            "to disable the file; got {!r}".format(value))
    return value


def npz_path(path):
    """The file ``numpy.savez`` writes for ``path``: ``.npz`` is appended
    unless the name already ends with it."""
    return path if path.endswith(".npz") else path + ".npz"


def _same_file_key(path):
    """A comparable form of ``path``: absolute, symbolic links and ``.``/``..``
    resolved, and case-folded on EVERY platform. ``os.path.normcase`` folds
    case only on Windows, but the default file systems of macOS are
    case-insensitive too, so two names that differ only in case are treated
    as the same file everywhere (conservative on a case-sensitive file
    system)."""
    return os.path.realpath(path).casefold()


def check_spin_flip_with_fixed_sz(trans_abs, norb):
    """Reject one-body spin-flip terms in 2Sz-fixed mode.

    Such terms genuinely break Sz conservation, so a fixed 2Sz is
    ill-defined.  Numerical noise must not abort the run: cross-spin
    entries are compared against the overall transfer scale (Wannier
    transfers carry ~1e-14 cross-spin residues, and a k-sum amplifies
    them past any absolute epsilon), so entries below 1e-8 of the scale
    are treated as zero (and masked by the block structure downstream).

    Parameters
    ----------
    trans_abs : ndarray, shape (2*norb, 2*norb)
        Non-negative magnitudes of the one-body terms in spin-block
        order (possibly summed over k).
    norb : int
        Number of orbitals per spin sector.

    Raises
    ------
    ValueError
        If a cross-spin entry exceeds 1e-8 of the largest entry.
    """
    scale = trans_abs.max() if trans_abs.size > 0 else 0.0
    if scale <= 0.0:
        return
    cross = max(trans_abs[:norb, norb:].max(initial=0.0),
                trans_abs[norb:, :norb].max(initial=0.0))
    if cross > 1.0e-8 * scale:
        msg = ("2Sz-fixed mode is incompatible with spin-flip one-body "
               "terms: Sz is not conserved, so a fixed 2Sz is "
               "ill-defined. Use the Sz-free mode (omit 2Sz) for such "
               "systems.")
        logger.error(msg)
        raise ValueError(msg)
    if cross > 0.0:
        # nonzero but below the noise threshold: masked, not fatal --
        # say so, in case it is a deliberately weak physical term
        logger.warning(
            "2Sz-fixed mode: cross-spin one-body entries up to {:.3e} "
            "({:.1e} relative to the transfer scale) are below the "
            "noise threshold (1e-8) and are ignored. If they are "
            "physical (e.g. weak spin-orbit coupling), use the Sz-free "
            "mode.".format(cross, cross / scale))


#: the default ``[mode.param]`` values of every solver (``solver_base``);
#: code outside the solvers that needs one of them reads it from here
PARAM_DEFAULTS = {
    "Nsite":  0,
    "Ncond":  0,
    "2Sz":    None,
    "Mix":    0.5,
    "EPS":    6,
    "IterationMax": 20000,
    "RndSeed": 1234,
    "T":      0.0,
    "ene_cutoff": 1.0e+2,
}

class solver_base():
    """Base solver class for Hartree-Fock calculations.
    
    Parameters
    ----------
    param_ham : dict
        Hamiltonian parameters
    info_log : dict
        Logging configuration
    info_mode : dict
        Calculation mode parameters
    param_mod : dict, optional
        Model parameters to override defaults
        
    Attributes
    ----------
    param_mod : CaseInsensitiveDict
        Model parameters
    param_ham : dict
        Hamiltonian parameters  
    info_log : dict
        Logging configuration
    threshold : float
        Cutoff threshold for Green's function elements
    relax_checks : bool
        Whether to relax parameter validation checks
    """

    def __init__(self, param_ham, info_log, info_mode, param_mod=None):
        logger = logging.getLogger("qlms").getChild(self.name)

        self.param_mod = param_mod
        self.param_ham = param_ham
        self.info_log = info_log
        info_mode_param = info_mode.get("param", CaseInsensitiveDict({}))

        # constants
        self.threshold = 1.0e-12  # cutoff on elements of green function

        # default values
        para_init = CaseInsensitiveDict(PARAM_DEFAULTS)

        range_list = {
            "T":     [ 0, None ],
            "IterationMax": [ 0, None ],
            "Mix":   [ 0.0, 1.0 ],
            "EPS":   [ 0, None ],
        }

        if param_mod is not None:
            for k, v in para_init.items():
                self.param_mod.setdefault(k, v)

            # integer parameters
            for key in ["Nsite", "Ncond", "2Sz", "EPS", "IterationMax", "RndSeed"]:
                if self.param_mod[key] is not None and type(self.param_mod[key]) == type([]):
                    self.param_mod[key] = int(self.param_mod[key][0])
            # float parameters
            for key in ["Mix", "T", "ene_cutoff"]:
                if self.param_mod[key] is not None and type(self.param_mod[key]) == type([]):
                    self.param_mod[key] = float(self.param_mod[key][0])

            # overwrite by info_mode_param
            for key, value in info_mode_param.items():
                self.param_mod[key] = value
        else:
            self.param_mod = ParamDict(para_init)
            if isinstance(info_mode_param, ParamDict):
                self.param_mod.accessed = info_mode_param.accessed
            for key, value in info_mode_param.items():
                self.param_mod[key] = value

        # relax parameter checks
        self.relax_checks = self.param_mod.get("trustme_interaction_range", False)
        if self.relax_checks:
            logger.warning("TRUST-ME mode enabled. parameter checks are relaxed")

        # Sz mode
        if "2Sz" in self.param_mod and self.param_mod["2Sz"] == "free":
            self.param_mod["2Sz"] = None

        # Ncond and filling
        if "Nelec" in self.param_mod:
            if self.param_mod["Ncond"] > 0:
                logger.warning("Both Ncond and Nelec are specified: Ncond used")
            else:
                self.param_mod["Ncond"] = self.param_mod["Nelec"]

        if "filling" in self.param_mod:
            if self.param_mod["Ncond"] > 0:
                logger.error("Both Ncond and filling are specified.")
                sys.exit(1)
            else:
                filling = self.param_mod["filling"]
                if filling < 0.0 or filling > 1.0:
                    logger.error("Parameter range check failed for filling")
                    sys.exit(1)
        else:
            if self.param_mod["Ncond"] == 0:
                logger.error("Neigher Ncond nor filling is specified.")
                sys.exit(1)

        # checks
        if self._check_info_mode(info_mode) != 0:
            logger.error("Parameter check failed for info_mode.")
            sys.exit(1)
        if self._check_param_range(self.param_mod, range_list) != 0:
            msg = "Parameter range check failed for param_mod."
            if self.relax_checks:
                logger.warning(msg)
            else:
                logger.error(msg)
                sys.exit(1)
        if self._check_param_mod(self.param_mod) != 0:
            msg = "Parameter check failed for param_mod."
            if self.relax_checks:
                logger.warning(msg)
            else:
                logger.error(msg)
                sys.exit(1)

        # canonicalize; the threshold must be a finite positive number even
        # when the range checks are relaxed (a NaN EPS passes them, and a
        # large EPS underflows to 0, so that no run could ever converge)
        eps_exponent = self.param_mod["EPS"]
        canonical = None
        try:
            canonical = pow(10, -eps_exponent)
            eps = float(canonical)
            usable = math.isfinite(float(eps_exponent)) and math.isfinite(eps) and eps > 0
        except (OverflowError, TypeError, ValueError):
            usable = False
        if not usable:
            logger.error("EPS = {!r} gives no usable convergence threshold "
                         "10^(-EPS) (it must be finite and greater than 0)"
                         .format(eps_exponent))
            sys.exit(1)
        self.param_mod["EPS"] = canonical

    def _check_info_mode(self, info_mode):
        """Check validity of info_mode parameters.
        
        Parameters
        ----------
        info_mode : dict
            Mode parameters to validate
            
        Returns
        -------
        int
            Number of validation errors found
        """
        logger = logging.getLogger("qlms").getChild(self.name)

        fix_list = {
            "mode": ["UHFr", "UHFk"],
        }

        exit_code = 0
        for key, _fix_list in fix_list.items():
            if key not in info_mode:
                logger.warning("mode.{} is not defined.".format(key))
                exit_code += 1
            elif info_mode[key] not in _fix_list:
                logger.warning("mode.{} in mode section is incorrect: {}.".format(key, _fix_list))
                exit_code += 1
        return exit_code

    def _check_param_mod(self, param_mod):
        """Check validity of model parameters.
        
        Parameters
        ----------
        param_mod : dict
            Model parameters to validate
            
        Returns
        -------
        int
            Number of validation errors found
        """
        logger = logging.getLogger("qlms").getChild(self.name)

        error_code = 0
        if "2Sz" in param_mod and param_mod["2Sz"] is not None:
            ncond_eo = param_mod["Ncond"] % 2
            twosz_eo = param_mod["2Sz"] % 2
            if ((ncond_eo == 0 and twosz_eo == 1) or
                (ncond_eo == 1 and twosz_eo == 0)):
                logger.error("mode.param.2Sz must be even(odd) when Ncond is even(odd).")
                error_code += 1
        return error_code

    def _check_param_range(self, param_mod, range_list):
        """Check if parameters are within valid ranges.
        
        Parameters
        ----------
        param_mod : dict
            Model parameters to validate
        range_list : dict
            Valid ranges for parameters
            
        Returns
        -------
        int
            Number of validation errors found
        """
        logger = logging.getLogger("qlms").getChild(self.name)

        error_code = 0
        for key, [ vmin, vmax ] in range_list.items():
            if key not in param_mod:
                logger.warning("mode.param.{} is not defined.".format(key))
                error_code += 1
            elif param_mod[key] is not None:
                if vmin is not None and param_mod[key] < vmin:
                    logger.warning("mode.param.{} must be greater than {}.".format(key, vmin))
                    error_code += 1
                if vmax is not None and param_mod[key] > vmax:
                    logger.warning("mode.param.{} must be smaller than {}.".format(key, vmax))
                    error_code += 1
        return error_code

    def _round_to_int(self, val, mode):
        """Round a value to integer according to specified mode.
        
        Parameters
        ----------
        val : float
            Value to round
        mode : str
            Rounding mode to use
            
        Returns
        -------
        int
            Rounded integer value
            
        Raises
        ------
        SystemExit
            If rounding fails or mode is invalid
        """
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

    def _check_spin_flip_with_fixed_sz(self, trans_abs, norb):
        """Reject one-body spin-flip terms in 2Sz-fixed mode: see
        :func:`check_spin_flip_with_fixed_sz`."""
        check_spin_flip_with_fixed_sz(trans_abs, norb)

    def solve(self, path_to_output):
        """Solve the Hartree-Fock equations.
        
        Parameters
        ----------
        path_to_output : str
            Path to output file
        """
        pass

    def get_results(self):
        """Get calculation results.
        
        Returns
        -------
        tuple
            (physics, Green's function) results
        """
        return (self.physics, self.Green)

    def save_results(self, info_outputfile, green_info):
        """Save calculation results.

        Parameters
        ----------
        info_outputfile : dict
            Output file configuration
        green_info : dict
            Green's function information
        """
        pass

    def _reset_scf_status(self):
        """Set the convergence record to "no step performed yet".

        Called by the UHF solvers before their SCF loop, so the record is
        defined even when the loop runs no step (``IterationMax = 0``)."""
        self.scf_converged = False
        self.scf_iterations = 0
        self.scf_residual = float("nan")

    def _record_scf_step(self, i_step, residual, eps):
        """Record a completed SCF step ``i_step`` (0-based) with residual
        ``residual``; returns True when it has converged (``residual < eps``).
        A non-finite residual (NaN or infinite) never converges."""
        self.scf_iterations = i_step + 1
        self.scf_residual = float(residual)
        if self.scf_residual < eps:
            self.scf_converged = True
        return self.scf_converged

    def _output_artifacts(self, info_outputfile):
        """``[(label, path)]`` of the result files ``save_results`` writes
        for ``info_outputfile``, besides the status file; the UHF solvers
        override it."""
        return []

    def check_output_paths(self, info_outputfile, print_check=None):
        """Refuse an unusable ``[file.output] scf_status`` without writing
        anything: a value that is not a string, or a status file that would
        be the same file as another result file of this run (as written:
        defaults, the ``.npz`` that ``numpy.savez`` appends, generated
        prefixes; names that differ only in case count as the same file) or
        as the ``[log] print_check`` file. Raises ValueError.

        Every name the solver can write is checked, also outputs that a
        given configuration does not produce (for example after zero SCF
        steps, or a UHFr fij file that an odd Ncond skips).

        Called by ``qlms.run`` before the SCF starts and again at the start
        of ``save_results``, so direct callers of the solver are covered."""
        name = scf_status_file_name(info_outputfile)
        if name == "":
            return
        if print_check is None:
            print_check = getattr(self, "info_log", {}).get("print_check")
        path_to_output = info_outputfile.get("path_to_output", "output")
        target = _same_file_key(os.path.join(path_to_output, name))
        artifacts = list(self._output_artifacts(info_outputfile))
        if isinstance(print_check, str):
            artifacts.append(("[log] print_check",
                              os.path.join(path_to_output, print_check)))
        for label, path in artifacts:
            if _same_file_key(path) == target:
                raise ValueError(
                    "[file.output] scf_status = {!r} is the same file as the "
                    "output of {} ({}); choose another name, or \"\" to "
                    "disable the status file".format(name, label, path))

    def _save_scf_status(self, info_outputfile, solver_name, residual_definition):
        """Write the convergence status file (format version 1) and report
        whether the other result files are to be written.

        The file name is ``[file.output] scf_status`` (default
        ``"scf_status.dat"``; ``""`` disables the file). Each line is
        ``"{key} = {value}\\n"``, in the order written below. Invariants:

        * ``eps`` is finite and > 0;
        * ``0 <= iterations <= iteration_max``;
        * ``iterations > 0`` implies ``residual`` finite and >= 0, OR
          ``residual = nan`` with the note ``residual not finite after
          step <iterations>`` (the last completed step produced a NaN or
          infinite residual; ``converged`` is then ``false``);
        * ``converged = true`` implies ``iterations > 0`` and
          ``residual < eps``;
        * ``iterations = 0`` implies ``converged = false`` and
          ``residual = nan``.

        A negative ``IterationMax`` (only reachable with relaxed parameter
        checks) runs no step; it is written as ``iteration_max = 0`` with
        the note ``IterationMax <value> treated as 0``.

        When no SCF step was performed, one WARNING is logged and False is
        returned: there is no result to write besides the status itself.
        """
        file_name = scf_status_file_name(info_outputfile)
        iterations = int(self.scf_iterations)
        iteration_max = int(self.param_mod["IterationMax"])
        residual = float(self.scf_residual)
        note = ""
        if iteration_max < 0:
            note = "IterationMax {} treated as 0".format(iteration_max)
        elif iterations == 0:
            note = "no SCF step performed"
        elif not math.isfinite(residual):
            residual = float("nan")
            note = "residual not finite after step {}".format(iterations)
        converged = bool(self.scf_converged) and math.isfinite(residual)
        if file_name != "":
            rows = (
                ("format_version", "1"),
                ("solver", solver_name),
                ("converged", "true" if converged else "false"),
                ("residual", repr(residual)),
                ("eps", repr(float(self.param_mod["eps"]))),
                ("iterations", str(iterations)),
                ("iteration_max", str(max(iteration_max, 0))),
                ("residual_definition", residual_definition),
                ("note", note),
            )
            path = os.path.join(info_outputfile["path_to_output"], file_name)
            with open(path, "w") as fw:
                fw.write("".join("{} = {}\n".format(k, v) for k, v in rows))
            logger.info("save_results: save convergence status in file {}".format(path))
        if iterations == 0:
            if file_name != "":
                rest = "no result file other than {} is written".format(file_name)
            else:
                rest = "no result file is written"
            logger.warning("no SCF step was performed (IterationMax = {}); {}"
                           .format(iteration_max, rest))
            return False
        return True
