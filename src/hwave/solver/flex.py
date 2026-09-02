"""FLEX (Fluctuation Exchange Approximation) solver.

This module implements the FLEX approximation for itinerant electron systems.
FLEX extends RPA by using dressed (self-consistent) Green's functions instead
of bare ones, iterating G -> chi0 -> chi -> Sigma -> G until convergence.

The solver inherits from the RPA class to reuse infrastructure for:
- Lattice and interaction setup
- Diagonalization and chemical potential search
- Bare susceptibility calculation
- RPA equation solving
- Block-diagonal structure detection
"""

from __future__ import annotations
from typing import Optional

import os
import numpy as np
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

from .rpa import (RPA, Lattice, Interaction, MOMENTUM_CONVENTION,
                  PAIRLIFT_INERT_WARNING)
from .density_projection import project_density_pairs
from . import backend as _bk
from . import bubble
from . import matsubara as _ms


def _eigvals_small(M):
    """Eigenvalues of a batched ``(..., nd, nd)`` matrix for the mu search.

    For ``nd <= 2`` the eigenvalues are evaluated in closed form ON THE
    INPUT'S ARRAY BACKEND -- the trace/determinant formula for 2x2 blocks and
    the diagonal element itself for 1x1 -- which both avoids the batched
    LAPACK geev cost and, on the GPU, avoids the device-to-host round trip
    (CuPy has no general eigensolver). For ``nd >= 3`` this falls back to
    host ``numpy.linalg.eigvals``.

    The mu search consumes the eigenvalues only as an unordered set
    (``N(mu) = sum_j 1/(lam_j + mu)``), so no eigenvalue ordering is
    guaranteed.

    Parameters
    ----------
    M : ndarray
        Batched square matrices, shape ``(..., nd, nd)`` (numpy or cupy).

    Returns
    -------
    lam : ndarray
        Eigenvalues, shape ``(..., nd)``; same backend as ``M`` for
        ``nd <= 2``, host numpy otherwise.
    """
    nd = M.shape[-1]
    if nd == 1:
        return M[..., 0]
    if nd == 2:
        xp = _bk.array_module_of(M)
        a = M[..., 0, 0]
        b = M[..., 0, 1]
        c = M[..., 1, 0]
        d = M[..., 1, 1]
        tr = a + d
        # Discriminant in the SHIFT-INVARIANT form (a-d)^2 + 4bc, not
        # tr^2 - 4*det: the mu-search operator is M = i*w*I - H0 - Sigma with
        # |w| up to (Nmat-1)*pi*T, and tr^2 - 4*det subtracts two nearly equal
        # O(w^2) terms, losing the O(1) remainder to cancellation (found by
        # review; pinned by test_closed_form_2x2_stable_under_matsubara_shift).
        dm = a - d
        disc = xp.sqrt(dm * dm + 4.0 * (b * c))
        return xp.stack([0.5 * (tr + disc), 0.5 * (tr - disc)], axis=-1)
    return np.linalg.eigvals(_bk.to_host(M))


class _AndersonMixer:
    r"""Anderson (Pulay/DIIS-type) acceleration of the FLEX SCF update.

    The SCF loop is a fixed-point iteration ``sigma -> F(sigma)``. Linear
    mixing takes ``x_{k+1} = x_k + beta r_k`` with the residual
    ``r_k = F(x_k) - x_k``. Anderson acceleration extrapolates over a short
    history of ``(x, r)`` pairs: with difference columns
    ``dR = [r_k - r_{k-1}, ...]`` and ``dX = [x_k - x_{k-1}, ...]`` it solves
    the (regularized) least-squares problem ``min_gamma |r_k - dR gamma|``
    via the normal equations and updates

        x_{k+1} = x_k + beta r_k - (dX + beta dR) gamma.

    With no usable history (first step, or ``depth=1``) this reduces exactly
    to the linear-mixing step. The histories live on the iterate's array
    backend (device arrays under GPU execution; only the tiny Gram system is
    solved on the host), so acceleration adds ``2*depth`` sigma-sized arrays
    of memory but no host round trips. If the extrapolated iterate is not
    finite (near-singular history), the mixer falls back to the plain linear
    step and restarts its history.
    """

    def __init__(self, mix, depth):
        self.mix = float(mix)
        self.depth = max(1, int(depth))
        self._xs = []   # iterate history (flattened, backend arrays)
        self._rs = []   # residual history

    def step(self, sigma, sigma_new):
        """Return the next iterate from ``sigma`` and ``F(sigma)=sigma_new``."""
        xp = _bk.array_module_of(sigma)
        shape = sigma.shape
        x = sigma.reshape(-1)
        r = sigma_new.reshape(-1) - x

        self._xs.append(x)
        self._rs.append(r)
        if len(self._xs) > self.depth:
            self._xs.pop(0)
            self._rs.pop(0)

        m = len(self._xs)
        if m == 1:
            return (x + self.mix * r).reshape(shape)

        # difference columns as (m-1, N) stacks
        dR = xp.stack([self._rs[i + 1] - self._rs[i] for i in range(m - 1)])
        dX = xp.stack([self._xs[i + 1] - self._xs[i] for i in range(m - 1)])

        # normal equations on the host (the system is (m-1) x (m-1))
        G = _bk.to_host(dR.conj() @ dR.T)
        b = _bk.to_host(dR.conj() @ r)
        # Tikhonov guard against a (near-)degenerate history
        lam = 1.0e-10 * max(float(np.trace(G).real) / (m - 1), 1.0e-300)
        try:
            gamma = np.linalg.solve(G + lam * np.eye(m - 1), b)
        except np.linalg.LinAlgError:
            gamma = None

        if gamma is not None:
            x_next = (x + self.mix * r
                      - (dX + self.mix * dR).T @ xp.asarray(gamma))
            if bool(xp.isfinite(x_next).all()):
                return x_next.reshape(shape)

        # fallback: plain linear step, restart the history from this point
        self._xs = [x]
        self._rs = [r]
        return (x + self.mix * r).reshape(shape)


def _auto_general_remediation(solver):
    """#167: the 2.0 default calc_scheme='auto' promotes every declared
    cross-carrying interaction to 'general', so a general-ONLY rejection can
    now stop a configuration that H-wave 1.0.x (== explicit 'reduced') ran
    without complaint. When that is how the run got here, name the
    resolution that chose 'general' and the one-line way back; an explicit
    'general' request was the user's own choice, so it gets no such hint.
    Mirrors the 4-axis chi0q remediation hint in rpa.py.
    """
    token = getattr(solver, "_scheme_resolution", None)
    if not str(token).startswith("auto:"):
        return ""
    return (" calc_scheme='auto' resolved to 'general' ({}) because a "
            "cross-carrying interaction is declared; to run this "
            "configuration as in H-wave 1.0.x request calc_scheme = "
            "'reduced' explicitly (a documented approximation)."
            .format(token))


def _scheme_stamp(solver):
    # #167: describes THIS run (never inside _freq_meta, which passes
    # through the INPUT file's provenance). Plain str -> <U. FLEX's
    # save_results is a distinct function from RPA's; this module-level
    # helper mirrors RPA's closure semantics without importing it.
    # getattr: tests drive save_results on __new__-built stubs (mirrors
    # rpa.py's fallback -- a stub without calc_scheme stamps "unknown"
    # rather than raising AttributeError).
    scheme = getattr(solver, "calc_scheme", "unknown")
    return {
        "calc_scheme": str(scheme),
        "calc_scheme_requested": str(getattr(solver, "calc_scheme_requested",
                                             scheme)),
        "scheme_resolution": str(getattr(solver, "_scheme_resolution", None)
                                 or "explicit"),
    }


class FLEX(RPA):
    """FLEX solver that extends RPA with self-consistent Green's functions.

    The FLEX approximation computes the self-energy from spin and charge
    fluctuations, then updates the Green's function self-consistently.

    The SCF loop is:
        1. Compute G(k, iwn) from H0 and Sigma
        2. Compute chi0(q, ivn) = -T/Nk sum_k G(k+q) G(k)
        3. Inflate chi0 and ham to common spin-orbital space
        4. Decompose ham into spin/charge channels
        5. Compute chi_s, chi_c via RPA equations
        6. Construct V_eff from chi_s, chi_c, chi0
        7. Compute Sigma(k, iwn) = T/Nk sum_q V_eff(q) G(k-q)
        8. Check convergence; if not, go to 1

    Parameters
    ----------
    param_ham : dict
        Hamiltonian parameters.
    info_log : dict
        Logging configuration.
    info_mode : dict
        Calculation mode parameters. FLEX-specific parameters include:
        - IterationMax: Maximum SCF iterations (default: 100)
        - Mix: Mixing parameter for self-energy update (default: 0.2)
        - EPS: Convergence criterion exponent (default: 6, i.e., 1e-6)
    """

    @do_profile
    def __init__(self, param_ham, info_log, info_mode):
        logger.debug(">>> FLEX.__init__")

        # Initialize RPA infrastructure (lattice, interaction, params)
        super().__init__(param_ham, info_log, info_mode)

        # FLEX-specific parameters
        self._init_flex_param()

    def _init_flex_param(self):
        """Initialize FLEX-specific parameters."""
        logger.debug(">>> FLEX._init_flex_param")

        # #167: FLEX has no transverse channel, so ring+ladder is invalid
        # for every scheme that reaches FLEX's own logic -- rejected here as
        # step 0, before any scheme logic runs.
        #
        # Note this guard is not the only stop, and deliberately not the
        # first one: an EXPLICIT non-general scheme with
        # calc_type='ring+ladder' is already rejected upstream by the
        # inherited RPA _set_scheme ("calc_type='ring+ladder' requires
        # calc_scheme='general' or 'auto'", logger.error + sys.exit(1)), so
        # explicit reduced never gets here. That precedence is 1.0.x
        # behaviour and is preserved on purpose (pinned by
        # tests/test_auto_scheme_exactness.py::TestFLEXAutoResolution::
        # test_explicit_reduced_ring_ladder_is_stopped_by_the_inherited_gate).
        # What this step-0 guard covers is the rest: 'auto' (which 1.0.x
        # resolved to general and which would otherwise reach the ring-only
        # general path) and an explicit 'general'.
        if getattr(self, "calc_type", "ring") == "ring+ladder":
            msg = ("FLEX does not support calc_type='ring+ladder' (the "
                   "transverse ladder channel); FLEX general is ring-only. "
                   "Use the default calc_type='ring' with "
                   "calc_scheme='general'.")
            logger.error(msg)
            raise ValueError(msg)

        # #167: 'auto' is resolved HERE, at construction -- FLEX's rule reads
        # only the declared interaction types, never H0(k), so there are no
        # late-arriving inputs to wait for (contrast RPA, whose predicate
        # depends on trans_mod / green_init and so resolves at
        # read_init/solve).
        if self.calc_scheme_requested == "auto":
            self._resolve_flex_auto_scheme()
        else:
            self._emit_flex_reduced_diagnostic()

        # FLEX consumes the reduced-shape (4-dim) chi0q and reduces the
        # interaction via the density-density diagonal ('kaabb->kab') on the
        # reduced path.  The 'general' scheme selects the paramagnetic
        # full-vertex multi-orbital path (v1: spin-free only; the spin_mode
        # guard is enforced in solve(), where spin_mode is determined).
        scheme = self.calc_scheme.lower()
        if scheme == "general":
            # FLEX general is the paramagnetic full-vertex path: it sums the
            # RPA *ring* (bubble) series with the full rank-4 vertex.
            if getattr(self.ham_info, "enable_spin_orbital", False):
                raise ValueError(
                    "calc_scheme='general' FLEX (v1) does not support "
                    "enable_spin_orbital; deferred to the generalized FLEX "
                    "solver.")
            self._flex_general = True
        elif scheme == "reduced":
            if getattr(self.ham_info, "enable_spin_orbital", False):
                raise ValueError(
                    "calc_scheme='reduced' FLEX does not support "
                    "enable_spin_orbital; deferred to the generalized FLEX "
                    "solver.")
            self._flex_general = False
        else:
            # ring+ladder was already rejected as step 0 above, and 'auto'
            # was resolved there too: only a genuinely unsupported scheme
            # name can reach here. The assertion is on the REQUESTED string
            # (the one FLEX's resolver keys off), so a mis-cased 'AUTO' --
            # which the inherited _set_scheme also treats as an explicit,
            # unsupported name -- still reaches the actionable ValueError
            # below instead of an AssertionError.
            assert self.calc_scheme_requested != "auto", \
                "auto must be resolved above"
            msg = ("FLEX requires calc_scheme='reduced' or 'general', "
                   "got '{}'.".format(self.calc_scheme))
            logger.error(msg)
            raise ValueError(msg)
        # Exchange/PairHop under reduced are rejected by the
        # consistency check inherited from RPA (one policy for both
        # solvers, #107): their vertex has no density-diagonal content, so
        # the schemes would drop them entirely -- the former warning here
        # called that an 'approximation', which for those types it is not.

        self.max_iter = int(self.param_mod.get("IterationMax", 100))
        self.mix = float(self.param_mod.get("Mix", 0.2))

        # SCF update scheme: plain linear mixing (default, unchanged), or
        # Anderson acceleration over a short iterate/residual history.
        self.mixing_scheme = str(
            self.param_mod.get("mixing_scheme", "linear")).lower()
        if self.mixing_scheme not in ("linear", "anderson"):
            raise ValueError(
                "mixing_scheme must be 'linear' or 'anderson', got '{}'."
                .format(self.mixing_scheme))
        self.anderson_depth = int(self.param_mod.get("anderson_depth", 5))

        # Matsubara-axis representation (design: docs/design/ir-matsubara.md,
        # Stage 2): the default uniform grid, or the sparse-ir intermediate
        # representation. IR computes chi0/Sigma natively on sparse nodes (no
        # uniform-FFT tau products), so the coeff_tail machinery is bypassed.
        self.matsubara_basis = str(
            self.param_mod.get("matsubara_basis", "uniform")).lower()
        if self.matsubara_basis not in ("uniform", "ir"):
            raise ValueError(
                "[mode.param] matsubara_basis must be 'uniform' or 'ir', "
                "got '{}'.".format(self.matsubara_basis))
        self.use_ir = (self.matsubara_basis == "ir")
        self.ir_tol = float(self.param_mod.get("ir_tol", 1.0e-8))
        self.ir_wmax = self.param_mod.get("ir_wmax")
        self.sigma_init_on_error = str(
            self.param_mod.get("sigma_init_on_error", "warn")).lower()
        if self.sigma_init_on_error not in ("warn", "abort", "zero"):
            raise ValueError(
                "[mode.param] sigma_init_on_error must be 'warn', 'abort' "
                "or 'zero', got '{}'.".format(self.sigma_init_on_error))
        # Stage 3 (design: docs/design/ir-matsubara-stage3.md): keep outputs
        # on the sparse nodes instead of densifying onto the Nmat grid.
        self.write_densified = _bk.as_bool(
            self.param_mod.get("write_densified", True))
        if not self.write_densified and not self.use_ir:
            raise ValueError(
                "[mode.param] write_densified = false requires "
                "[mode.param] matsubara_basis = 'ir'.")
        self._ir_axF = None
        self._ir_axB = None

        eps_exp = self.param_mod.get("EPS", 6)
        if isinstance(eps_exp, float) and eps_exp < 1.0:
            self.eps = eps_exp
        else:
            self.eps = 10.0 ** (-int(eps_exp))

        # Lazy cache for the q-independent MYO S/C matrices (general path only);
        # built once on first use in _inflate_chi0q_and_ham_general and reused
        # across SCF iterations. None until then.
        self._myo_sc_cache = None

        logger.info("FLEX parameters:")
        logger.info("    max_iter        = {}".format(self.max_iter))
        logger.info("    mix             = {}".format(self.mix))
        logger.info("    eps             = {:e}".format(self.eps))
        # The inherited _show_params ran before FLEX resolved 'auto', so it
        # printed the deferred marker; restate the settled value and the
        # token that produced it (#167).
        logger.info("    calc_scheme     = {} ({})".format(
            self.calc_scheme, self._scheme_resolution))

    def _resolve_flex_auto_scheme(self):
        """Constructor-time FLEX auto resolution (#167): general iff any
        declared type has flex_forcing.

        FLEX's rule is H0-INDEPENDENT -- the dressed Green function can
        hybridise during the SCF iteration however diagonal H0(k) starts
        out -- so the decision reads only the declared interaction types and
        can be taken here, at construction. Atomic; no fingerprint is
        recorded (there are no late-arriving inputs to re-check), which also
        makes the inherited RPA resolver a no-op at read_init/solve.
        """
        from hwave.solver import scheme as _scheme
        types = _scheme.declared_types(self._scheme_source_tables())
        chosen, token = _scheme.resolve_flex(types)
        assert token in _scheme.RESOLUTION_TOKENS
        # commit
        self.calc_scheme = chosen
        self.enable_reduced = (chosen == "reduced")
        self._scheme_resolution = token
        self._scheme_fingerprint = None
        # side effects (best-effort: a logging failure must never sink a
        # committed, valid resolution)
        try:
            if chosen == "reduced" and "PairLift" in types:
                logger.warning(PAIRLIFT_INERT_WARNING)
            legacy = _scheme.legacy_1_0_resolution(types, self.calc_type)
            est = self._log_scheme_memory_estimate(
                "auto-resolved calc_scheme (FLEX)")
            if chosen != legacy:
                logger.warning(
                    "FLEX calc_scheme='auto' resolved to '{}' ({}); H-wave 1.0.x "
                    "resolved this input to '{}'. Reason: {} carry cross-family "
                    "vertex content and the dressed Green function can hybridise "
                    "during iteration. chiq output becomes rank-6; memory/solve "
                    "cost grow from O(nd^2)/O(nd^3) to O(nd^4)/O(nd^6) per "
                    "(q, freq) (estimated principal array {:.3f} GB). To keep the "
                    "1.0.x behaviour request calc_scheme = 'reduced' explicitly "
                    "(a documented approximation).".format(
                        chosen, token, legacy,
                        ", ".join(sorted(t for t in types
                                         if _scheme.CAPABILITIES[t].flex_forcing)),
                        est / 1e9))
            else:
                logger.info("auto mode for calc_scheme: set to {} ({})".format(
                    chosen, token))
        except Exception as exc:
            logger.error("scheme resolution logging failed: {!r}".format(exc))

    def _emit_flex_reduced_diagnostic(self):
        """Explicit reduced + a flex_forcing type: an approximation
        UNCONDITIONAL of H0 (#167 §5). Emitted once, at construction.

        ADVISORY ONLY (as RPA's counterpart): fail-closed discovery belongs
        to the AUTO path, which decides on its verdict. Here the user named
        the scheme, so an input the registry cannot judge skips the
        diagnostic at debug level rather than failing a construction that
        1.0.x completed.
        """
        from hwave.solver import scheme as _scheme
        if self.calc_scheme_requested != "reduced":
            return
        try:
            types = _scheme.declared_types(self._scheme_source_tables())
        except ValueError as exc:
            logger.debug("reduced-exactness diagnostic skipped: %s", exc)
            return
        forcing = sorted(t for t in types
                         if _scheme.CAPABILITIES[t].flex_forcing)
        if forcing:
            logger.warning(
                "calc_scheme='reduced' with {} is an approximation in FLEX "
                "regardless of the transfer structure (the dressed Green "
                "function can hybridise during iteration). calc_scheme="
                "'general' or 'auto' carries the full vertex.".format(
                    ", ".join(forcing)))

    def _emit_reduced_exactness_diagnostic(self, *, trans_mod_present,
                                           green_init_present):
        # RPA's explicit-reduced diagnostic is H0-CONDITIONAL and its text
        # cites RPA instabilities and RPA fixture figures; neither applies
        # to FLEX, which warns at construction, unconditionally of H0
        # (_emit_flex_reduced_diagnostic).
        return

    def _check_reduced_rejects_spinful(self):
        """Reject ``self.spin_mode == 'spinful'`` under
        ``calc_scheme='reduced'`` (``self._flex_general`` False) at SOLVE
        time.

        Split out as its own method (rather than inlined in
        ``_solve_impl``) so a unit test can pin this exact check by
        constructing the semantic state directly (setting
        ``_flex_general``/``spin_mode`` on a bare instance) instead of
        driving a full solve.

        This complements the CONSTRUCTOR-time
        ``enable_spin_orbital=True`` guard in ``_init_flex_param``, but is
        NOT merely defense-in-depth for an otherwise-unreachable state:
        ``self.spin_mode`` is determined here from H0(k)'s block
        structure by ``_calc_epsilon_k``, which is driven not only by
        ``enable_spin_orbital`` but also by the PUBLIC ``trans_mod`` /
        ``green_init`` inputs (``RPA.read_init``) -- a spin-mixing
        ``trans_mod`` reaches ``spin_mode == 'spinful'`` here even with
        ``enable_spin_orbital=False`` (see
        ``tests/test_rpa_flex_equivalence_table.py``'s
        ``TestReducedSpinfulGuard.test_public_trans_mod_route_reaches_spinful_and_is_rejected``).
        The reduced scheme projects the interaction onto its
        density-diagonal part (``_inflate_chi0q_and_ham``); that
        projection has never been validated against a genuinely
        spin-mixing H0, so this rejects rather than silently computing
        unvalidated physics. spin-free and spin-diag are unaffected.
        """
        if (not self._flex_general) and self.spin_mode == "spinful":
            raise ValueError(
                "calc_scheme='reduced' FLEX does not support "
                "spin_mode='spinful'; spin-free and spin-diag are "
                "supported. Deferred to the generalized FLEX solver.")

    def read_init(self, info_inputfile):
        """Read initial configs, plus the FLEX-specific ``sigma_init``.

        In addition to the inherited RPA inputs (``green_init`` / ``trans_mod``,
        used to build H0(k)), FLEX can warm-start its SCF loop from a previously
        saved self-energy: pass ``[file.input] sigma_init = "sigma.npz"`` (a
        ``sigma.npz`` written by an earlier FLEX run). Starting from a converged
        neighbouring solution (e.g. the next-higher temperature) instead of
        ``Sigma = 0`` is often decisive for convergence near a magnetic
        instability, where the zero-start transient makes the SCF oscillate.
        """
        info = super().read_init(info_inputfile)
        if "sigma_init" in info_inputfile:
            path_to_input = info_inputfile.get("path_to_input", "")
            file_name = os.path.join(path_to_input,
                                     info_inputfile["sigma_init"])
            sigma, ir_meta = self._read_sigma(file_name)
            if ir_meta is not None and not self.use_ir:
                raise ValueError(
                    "sigma_init file '{}' holds sparse-IR node data "
                    "(frequency_grid=sparse_ir_nodes); an IR-native seed "
                    "requires [mode.param] matsubara_basis = 'ir' in this "
                    "run, or re-run the seeding FLEX with [mode.param] "
                    "write_densified = true.".format(file_name))
            info["sigma_init"] = sigma
            if ir_meta is not None:
                info["sigma_init_ir"] = ir_meta
        return info

    def _read_sigma(self, file_name):
        """Load a saved self-energy array from a FLEX ``sigma.npz``.

        Returns ``(sigma, ir_meta)``: ``ir_meta`` is ``None`` for uniform
        (densified) files, or the node metadata of an IR-native file
        (Stage 3; ``read_init`` decides whether this run can consume it).

        Validates the recorded ``cell_shape`` against this run's lattice:
        ``Nvol = Lx*Ly*Lz`` is a single dimension of the sigma array, so an
        aspect-ratio change (e.g. ``[2,8,1]`` vs ``[4,4,1]``) would pass the
        plain shape check in ``solve()`` while silently mapping the seed onto
        the wrong k-points. Files written before the key existed load with a
        warning instead (the grid cannot be verified).
        """
        from hwave.solver.ir_axis import is_ir_native, ir_native_meta
        logger.info("FLEX: read initial self-energy from {}".format(file_name))
        data = np.load(file_name)
        if "cell_shape" in data:
            saved = tuple(int(x) for x in data["cell_shape"])
            if saved != tuple(self.lattice.shape):
                raise ValueError(
                    "sigma_init was written on CellShape {} but this run uses "
                    "{}; the k-point layout differs even if the total volume "
                    "matches. Regenerate sigma_init at this CellShape.".format(
                        list(saved), list(self.lattice.shape)))
        else:
            logger.warning(
                "sigma_init file '{}' carries no cell_shape metadata (written "
                "by an older H-wave); make sure its CellShape matches this "
                "run's {} -- a mismatched k-point layout cannot be detected "
                "from the array shape alone.".format(
                    file_name, list(self.lattice.shape)))
        meta = ir_native_meta(data) if is_ir_native(data) else None
        # Fourier-sign provenance gate (issue #133): sigma is a k-labeled
        # quantity; a pre-#133 file carries flipped labels. Legacy files
        # are accepted only when the payload is q-even (see the validator).
        from hwave.solver.rpa import validate_momentum_convention
        sig = data["sigma"]
        # layout (nblock, nmat_or_nodes, nvol, nd, nd): nvol is axis 2 --
        # do NOT search by size, nmat == nvol collisions are realistic
        validate_momentum_convention(data, file_name, sig, 2,
                                     self.lattice.shape)
        return sig, meta

    @do_profile
    def solve(self, green_info, path_to_output):
        """Solve the FLEX equations, restoring host-backed public state.

        Thin wrapper around :meth:`_solve_restoring_host_attrs` that
        guarantees the solver's
        public array attributes (``H0_eigenvalue``/``H0_eigenvector``, and the
        stored ``green0``/``green0_tail``) are NumPy-backed after the call --
        on normal completion AND after a GPU-path exception. Under GPU
        execution ``_solve_impl`` converts these to CuPy in place; without the
        ``finally`` a mid-solve error would leave a reused or inspected solver
        object holding device arrays (issue #63).
        """
        return self._solve_restoring_host_attrs(green_info, path_to_output)

    def _solve_impl(self, green_info, path_to_output):
        """Solve the FLEX equations self-consistently.

        Parameters
        ----------
        green_info : dict
            Dictionary containing Green's function information.
        path_to_output : str
            Path to output directory.

        Notes
        -----
        FLEX starts the SCF loop from zero self-energy and recomputes chi0q from
        the dressed Green's function every iteration, so a `chi0q_init` entry
        loaded by the inherited RPA `read_init` is NOT consumed.  `green_init`
        and `trans_mod`, however, ARE consumed: the inherited `_calc_epsilon_k`
        uses them to build the (mean-field-shifted) transfer H0(k), exactly as
        in RPA.
        """
        logger.info("Start FLEX calculations")

        beta = 1.0 / self.T
        nvol = self.lattice.nvol
        nmat = self.nmat
        norb = self.norb
        ns = self.ns
        nd = self.nd

        # Step 1: Compute band structure and chemical potential.
        # _calc_epsilon_k determines self.spin_mode from H0(k); the general
        # full-vertex path (v1) is paramagnetic only, so guard here once
        # spin_mode is known.
        self._calc_epsilon_k(green_info)

        if self._flex_general and self.spin_mode != "spin-free":
            raise ValueError(
                "calc_scheme='general' FLEX (v1) supports spin_mode='spin-free' "
                "only, got '{}'. spin-diag/spinful are deferred to the "
                "generalized FLEX solver.{}".format(
                    self.spin_mode, _auto_general_remediation(self)))

        self._check_reduced_rejects_spinful()

        if self.use_ir:
            self._ir_setup(beta)

        if self.calc_mu:
            # spin-free counts one spin, so the target is halved (as in
            # RPA._find_mu / RPA.solve).  Ncond_target is reused every SCF
            # iteration to re-solve mu from the DRESSED Green's function.
            if self.spin_mode == "spin-free":
                Ncond_target = self.Ncond / 2
            else:
                Ncond_target = self.Ncond
            # initial mu from the non-interacting bands (Sigma = 0)
            dist, mu = self._find_mu(Ncond_target, self.T)
        else:
            Ncond_target = None
            mu = self.mu_value

        self.mu = mu

        # GPU (CuPy) execution: resolve the backend once, then move the SCF
        # loop's working arrays to the device. The eigendecomposition of H0 and
        # the initial mu search above stay on the host; from here on every
        # array-producing method dispatches on its inputs' array module, so the
        # loop below runs end-to-end on the GPU. The chemical-potential search
        # (_find_mu_dressed) is the one exception: it needs a general
        # (non-Hermitian) eigensolver, which CuPy does not provide, so it
        # brings its operator to the host internally.
        xp, gpu_active = _bk.get_backend(self.use_gpu, logger=logger,
                                         required=self.gpu_required)
        if gpu_active:
            logger.info("FLEX: GPU backend active (CuPy); moving H0 "
                        "eigenpairs and interaction to the device.")
            # VRAM preflight: estimate the larger of a dressed-G/sigma tensor
            # and a chi/vertex tensor, then allow for several live arrays.
            # The reduced chi/vertex is inflated to spin-orbital space, while
            # the general path carries a rank-4 orbital block. The 5x factor is
            # a rough advisory estimate; transient FFT/einsum/solve workspace
            # is not included, so actual peak usage may be higher.
            nblk0, _, _, nd0 = self.H0_eigenvector.shape
            nfreq_f = self._ir_axF.n_freq if self.use_ir else nmat
            nfreq_b = self._ir_axB.n_freq if self.use_ir else nmat
            green_elems = nblk0 * nfreq_f * nvol * nd0**2
            if self._flex_general:
                vertex_elems = nfreq_b * nvol * nd0**4
            else:
                inflated_nd = (
                    nd0 if self.spin_mode == "spinful" else self.ns * self.norb
                )
                vertex_elems = nfreq_b * nvol * inflated_nd**2
            resident_bytes = max(green_elems, vertex_elems) * 16
            _bk.warn_if_device_memory_short(
                5 * resident_bytes, logger, label="the FLEX SCF loop")
            self.H0_eigenvalue = xp.asarray(self.H0_eigenvalue)
            self.H0_eigenvector = xp.asarray(self.H0_eigenvector)

        # Step 2: Compute bare Green's function G0(k, iwn)
        if self.use_ir:
            green0, green0_tail = self._calc_green_ir(beta, mu)
        else:
            green0, green0_tail = self._calc_green(beta, mu)

        # Store for reference (as host arrays; the loop below keeps using the
        # backend-local green0_tail)
        self.green0 = _bk.to_host(green0)
        self.green0_tail = (None if green0_tail is None
                            else _bk.to_host(green0_tail))

        # High-frequency tail contract (coeff_tail): RPA's tail acceleration
        # is a two-step pair -- _calc_green subtracts aa/(i w_n) in FREQUENCY
        # space, and _calc_chi0q subtracts the analytic tau-space constant
        # green0_tail (= VV† aa beta/2) after the Matsubara FFT.  The dressed
        # Green's function below comes from _calc_dressed_green, which returns
        # the FULL physical G, so the frequency-space term must be subtracted
        # here before handing G to _calc_chi0q; otherwise G(tau) is uniformly
        # shifted by -aa/2 and chi0q is O(1) wrong.
        #
        # The tail subtraction is applied ONLY to the chi0q transform (the one
        # paired with green0_tail).  The self-energy convolution keeps the
        # full physical G: its FFT pipeline is an exact cyclic frequency
        # convolution and needs no tau-space tail reconstruction -- measured
        # on the 8x8 Hubbard fixture, reconstructing the "true" G(tau) there
        # does not improve the Nmat convergence of Sigma.
        aa = 0.0 if self.use_ir else self.coeff_tail
        if aa != 0.0:
            iomega = (xp.arange(nmat) * 2 + 1 - nmat) * np.pi / beta
            ev = self.H0_eigenvector
            VVt = ev @ xp.conj(ev).swapaxes(-2, -1)  # (nblock, nvol, nd, nd)
            green_tail_w = ((aa / (1j * iomega))[np.newaxis, :, np.newaxis,
                                                 np.newaxis, np.newaxis]
                            * VVt[:, np.newaxis])
        else:
            green_tail_w = None

        # Initialize the self-energy: zero by default, or warm-start from a
        # provided ``sigma_init`` (see read_init). Shape must match this run's
        # (nblock, nmat, nvol, nd_block, nd_block) -- warm-start needs the same
        # k-mesh and Matsubara grid.
        nblock = green0.shape[0]
        nd_block = green0.shape[-1]
        nfreq_axis = self._ir_axF.n_freq if self.use_ir else nmat
        expected_uniform = (nblock, nmat, nvol, nd_block, nd_block)
        expected = (nblock, nfreq_axis, nvol, nd_block, nd_block)
        sigma_init = green_info.get("sigma_init")
        seed_ir_meta = green_info.get("sigma_init_ir")
        if sigma_init is not None:
            if seed_ir_meta is not None:
                # IR-native seed (Stage 3): the file holds node values;
                # read_init already guaranteed this run is IR.
                expected_native = (nblock, seed_ir_meta["freq_n"].size,
                                   nvol, nd_block, nd_block)
                if sigma_init.shape != expected_native:
                    raise ValueError(
                        "IR-native sigma_init shape {} does not match its "
                        "own metadata {} (nblock, n_nodes, Nvol, nd, nd); "
                        "the file is inconsistent or from a different "
                        "k-mesh. Regenerate sigma_init at this CellShape."
                        .format(sigma_init.shape, expected_native))
                seed = self._ir_sigma_init_native(
                    np.array(sigma_init, dtype=np.complex128, copy=True),
                    seed_ir_meta, beta)
            else:
                if sigma_init.shape != expected_uniform:
                    raise ValueError(
                        "sigma_init shape {} does not match this run's {} "
                        "(nblock, Nmat, Nvol, nd, nd); warm-start requires "
                        "the same k-mesh and Matsubara grid. Regenerate "
                        "sigma_init at this Nmat/CellShape.".format(
                            sigma_init.shape, expected_uniform))
                seed = np.array(sigma_init, dtype=np.complex128, copy=True)
                if self.use_ir:
                    seed = self._ir_sigma_init(seed)
            # xp.asarray moves the (host-loaded) seed to the device under GPU
            # execution; on numpy it is a plain cast.
            if seed is None:
                sigma = xp.zeros(expected, dtype=np.complex128)
            else:
                sigma = xp.asarray(seed)
                logger.info("FLEX: warm-starting the SCF loop from sigma_init")
        else:
            sigma = xp.zeros(expected, dtype=np.complex128)

        # Prepare interaction Hamiltonian (full spin-orbital space) as a
        # LOCAL backend copy: ham_info is shared state handed in by the
        # caller, so it must not be mutated to a device array.
        ham_orig = self.ham_info.ham_inter_q
        if gpu_active:
            ham_orig = xp.asarray(ham_orig)

        # Main SCF loop
        diff = float("inf")
        converged = False
        n_iter_done = 0
        mixer = (_AndersonMixer(self.mix, self.anderson_depth)
                 if self.mixing_scheme == "anderson" else None)
        for iteration in range(self.max_iter):
            logger.info("FLEX iteration {}/{}".format(iteration + 1, self.max_iter))

            # Re-solve mu so the DRESSED G reproduces the target particle
            # number: as Sigma grows the frozen non-interacting mu no longer
            # yields Ncond electrons, so the run would otherwise converge to a
            # different filling than requested.  mu is solved for the current
            # sigma BEFORE building green_kw so the stored (mu, green_kw) pair
            # is self-consistent (N(green_kw) == Ncond).  A fixed mu
            # (calc_mu=False) is left untouched.
            if self.calc_mu:
                mu = self._find_mu_dressed(sigma, beta, Ncond_target)
                self.mu = mu

            # Step 3: Compute dressed Green's function G(k, iwn)
            green_kw = self._calc_dressed_green(beta, mu, sigma)

            # Tail-subtracted G for the Matsubara-FFT transforms (see the
            # coeff_tail contract above); identical to green_kw when aa == 0.
            if green_tail_w is not None:
                green_scf = green_kw - green_tail_w
            else:
                green_scf = green_kw

            # Step 4: Compute chi0(q, ivn) from dressed G. Four-way dispatch:
            # (general/reduced) x (uniform/ir). The reduced IR chi0 is a
            # density-density shortcut; general needs the full rank-6 bubble.
            if self.use_ir:
                if self._flex_general:
                    chi0q_raw = self._calc_chi0q_general_ir(green_scf, beta)
                else:
                    chi0q_raw = self._calc_chi0q_ir(green_scf, beta)
            else:
                chi0q_raw = self._calc_chi0q(green_scf, green0_tail, beta)

            # Remove spin block dimension
            if self.spin_mode in ["spin-free", "spinful"]:
                assert chi0q_raw.shape[0] == 1
                chi0q_raw = chi0q_raw[0]
            else:
                assert chi0q_raw.shape[0] == 2

            # Step 5: Inflate chi0q and ham, then compute spin/charge
            # susceptibilities and V_eff.  The general (paramagnetic full-vertex)
            # path keeps the full rank-4 orbital vertex; the reduced
            # path uses the density-density reduction.
            # Step 6: Compute self-energy Sigma(k, iwn).
            if self._flex_general:
                chi0q_out, v_eff, chi_s, chi_c = \
                    self._flex_compute_veff_general(chi0q_raw, ham_orig)
                if self.use_ir:
                    sigma_new = self._calc_self_energy_general_ir(
                        green_kw, v_eff, beta)
                else:
                    sigma_new = self._calc_self_energy_general(
                        green_kw, v_eff, beta)
            else:
                chi0q_out, v_eff, chi_s, chi_c = self._flex_compute_veff(
                    chi0q_raw, ham_orig)
                if self.use_ir:
                    sigma_new = self._calc_self_energy_ir(green_kw, v_eff,
                                                          beta)
                else:
                    sigma_new = self._calc_self_energy(green_kw, v_eff, beta)

            if self.use_ir:
                self._ir_coeff_decay_check(sigma_new)

            # Step 7: Mix and check convergence
            diff = self._calc_convergence(sigma, sigma_new)
            logger.info("  convergence: |dSigma|/|Sigma| = {:.3e}".format(diff))

            n_iter_done = iteration + 1
            if mixer is not None:
                sigma = mixer.step(sigma, sigma_new)
            else:
                sigma = (1.0 - self.mix) * sigma + self.mix * sigma_new

            if diff < self.eps:
                logger.info("FLEX converged after {} iterations".format(
                    iteration + 1))
                converged = True
                break

        if not converged:
            logger.warning("FLEX did not converge after {} iterations "
                           "(diff={:.3e}, eps={:.3e})".format(
                               self.max_iter, diff, self.eps))

        # SCF outcome, for programmatic consumers (sweeps, tests, mixer
        # comparisons)
        self.scf_converged = converged
        self.scf_iterations = n_iter_done

        if self.max_iter == 0:
            # No SCF iteration ran: green_kw / chi_s / chi_c / chi0q_out were
            # never computed.  Warn-and-return instead of dereferencing them.
            if gpu_active:
                self.H0_eigenvalue = _bk.to_host(self.H0_eigenvalue)
                self.H0_eigenvector = _bk.to_host(self.H0_eigenvector)
            logger.warning("FLEX IterationMax=0: no SCF step performed; "
                           "no results stored.")
            return

        # Final-output consistency: during the loop green_kw was built from the
        # PRE-mix sigma, while `sigma` below is the POST-mix estimate, so the
        # stored (mu, green, sigma) triple would not satisfy the Dyson equation
        # green = [G0^{-1} - sigma]^{-1} (noticeable for non-converged or
        # IterationMax=1 runs; negligible once converged).  Rebuild the dressed
        # G from the final stored sigma -- and, for calc_mu, re-solve mu for it
        # -- so the stored triple is mutually consistent and N(green) == Ncond.
        if self.calc_mu:
            mu = self._find_mu_dressed(sigma, beta, Ncond_target)
            self.mu = mu
        green_kw = self._calc_dressed_green(beta, mu, sigma)

        # Physical observables from the final (consistent) dressed G: particle
        # number N and spin Sz.  In fixed-mu mode this N is the mu-N single
        # point; in calc_mu mode it equals the target Ncond.
        physics = self._calc_physics_dressed(green_kw, mu, beta)
        self.physics = physics
        logger.info("FLEX: NCond = {}, Sz = {}, ChemicalPotential = {}".format(
            physics["NCond"], physics["Sz"], physics["mu"]))

        # Restore the solver's public attributes to host arrays so the
        # post-solve object state is backend-independent.
        if gpu_active:
            self.H0_eigenvalue = _bk.to_host(self.H0_eigenvalue)
            self.H0_eigenvector = _bk.to_host(self.H0_eigenvector)

        # Store results (as host arrays: everything downstream -- writers,
        # green_info consumers -- is numpy). On the IR path every stored
        # frequency axis is densified back onto the run's uniform Nmat grid
        # so the output files keep their exact format (design OQ-1) and the
        # Stage-1 dynamic Eliashberg loader works unchanged.
        if self.use_ir and self.write_densified:
            axF, axB = self._ir_axF, self._ir_axB
            sigma = self._ir_densify(sigma, axF, 1)
            green_kw = self._ir_densify(green_kw, axF, 1)
            chi_s = self._ir_densify(chi_s, axB, 0)
            chi_c = self._ir_densify(chi_c, axB, 0)
            chi0q_out = self._ir_densify(chi0q_out, axB, 0)
        else:
            # uniform run, or IR-native (write_densified=false): the arrays
            # stay on their working frequency axis (sparse nodes for IR),
            # host-copied for the numpy-only writers/consumers.
            sigma = _bk.to_host(sigma)
            green_kw = _bk.to_host(green_kw)
            chi_s = _bk.to_host(chi_s)
            chi_c = _bk.to_host(chi_c)
            chi0q_out = _bk.to_host(chi0q_out)
        self.sigma = sigma
        self.green_kw = green_kw
        self.chi_s = chi_s
        self.chi_c = chi_c

        # Store in green_info for output
        green_info["chi0q"] = chi0q_out
        green_info["chiq_s"] = chi_s
        green_info["chiq_c"] = chi_c
        green_info["sigma"] = sigma
        green_info["green"] = green_kw
        green_info["physics"] = physics

        logger.info("End FLEX calculations")

    # ------------------------------------------------------------------
    # IR-basis Matsubara axis (Stage 2, docs/design/ir-matsubara.md 3.3/3.4)
    # ------------------------------------------------------------------

    def _ir_setup(self, beta):
        """Build the fermionic/bosonic IR axes once per solve."""
        if self._ir_axF is not None:
            return
        from hwave.solver.ir_axis import IRAxis
        wmax = self.ir_wmax
        if wmax is None:
            ew = _bk.to_host(self.H0_eigenvalue)
            band = 2.0 * float(np.abs(ew).max())
            u = float(np.abs(_bk.to_host(
                self.ham_info.ham_inter_q)).max())
            wmax = 3.0 * (band + u)
            if not np.isfinite(wmax) or wmax <= 0.0:
                raise ValueError(
                    "ir_wmax auto-estimate is not a positive finite number "
                    "(band bound {} + interaction scale {}); set "
                    "[mode.param] ir_wmax explicitly (a real-frequency "
                    "bandwidth in the same energy units as the "
                    "Hamiltonian).".format(band, u))
            logger.info("IR: auto ir_wmax = %.6g (override with "
                        "[mode.param] ir_wmax)", wmax)
        wmax = float(wmax)
        self._ir_axF = IRAxis(beta=beta, wmax=wmax, eps=self.ir_tol,
                              statistics="F")
        self._ir_axB = IRAxis(beta=beta, wmax=wmax, eps=self.ir_tol,
                              statistics="B")
        logger.info("IR: Lambda=%.3g eps=%.1e -> L_F=%d (nodes %d), "
                    "L_B=%d (nodes %d)", beta * wmax, self.ir_tol,
                    self._ir_axF.L, self._ir_axF.n_freq,
                    self._ir_axB.L, self._ir_axB.n_freq)
        if self.coeff_tail != 0.0:
            logger.info("IR: coeff_tail is not used on the IR path (the "
                        "fermionic basis carries the 1/(i w) tail exactly); "
                        "the option is ignored.")

    def _freq_omegas(self, beta, xp):
        """Matsubara frequencies of the working axis (uniform or IR nodes)."""
        if self._ir_axF is not None:
            return xp.asarray(self._ir_axF.freq_n) * (np.pi / beta)
        nmat = self.nmat
        return (xp.arange(nmat) * 2 + 1 - nmat) * np.pi / beta

    @do_profile
    def _calc_green_ir(self, beta, mu):
        """Bare Green's function on the fermionic IR nodes (no tail terms:
        the basis represents the 1/(i w) asymptotics within ir_tol)."""
        ew = self.H0_eigenvalue
        ev = self.H0_eigenvector
        xp = _bk.array_module_of(ew)
        iomega = xp.asarray(self._ir_axF.freq_n) * (np.pi / beta)
        wn = 1j * iomega[np.newaxis, :, np.newaxis, np.newaxis]
        ek = (ew - mu)[:, np.newaxis, :, :]
        g = 1.0 / (wn - ek)
        ev_conj_t = xp.conj(ev).swapaxes(-2, -1)
        Vg = ev[:, np.newaxis, :, :, :] * g[:, :, :, np.newaxis, :]
        green = Vg @ ev_conj_t[:, np.newaxis, :, :, :]
        return green, None

    @do_profile
    def _calc_chi0q_ir(self, green_kw, beta):
        r"""chi0 on the bosonic IR nodes, computed natively from G on the
        fermionic nodes (reduced scheme).

        Validates the ``enable_reduced`` requirement (the user-facing
        diagnostic below is kept even though the kernel has no reduced-only
        restriction of its own), then delegates to ``bubble.ir_bubble``.
        """
        if not self.enable_reduced:
            raise ValueError(
                "matsubara_basis='ir' chi0 supports the reduced scheme only")
        nx, ny, nz = self.lattice.shape
        workers = getattr(self, "fft_workers", 1)
        return bubble.ir_bubble(
            green_kw, self._ir_axF, self._ir_axB,
            spatial_shape=(nx, ny, nz), scheme="reduced", workers=workers)

    @do_profile
    def _calc_chi0q_general_ir(self, green_kw, beta):
        # General (full-vertex/MYO) counterpart of `_calc_chi0q_ir`; called
        # only on the general path (mirrors that method's `enable_reduced`
        # guard, which is why this one has no analogous check).
        r"""chi0 (full rank-6 orbital bubble) on the bosonic IR nodes, general
        (MYO) scheme.

        Delegates to ``bubble.ir_bubble`` (no ``enable_reduced`` guard, as
        with the legacy body -- this method is only called on the general
        path).
        """
        nx, ny, nz = self.lattice.shape
        workers = getattr(self, "fft_workers", 1)
        return bubble.ir_bubble(
            green_kw, self._ir_axF, self._ir_axB,
            spatial_shape=(nx, ny, nz), scheme="general", workers=workers)

    @do_profile
    def _calc_self_energy_ir(self, green_kw, v_eff, beta):
        """Sigma on the fermionic IR nodes: V (bosonic nodes) and G
        (fermionic nodes) are both evaluated on the FERMIONIC tau nodes
        (the product V*G is anti-periodic), multiplied there, and fitted
        back. Physical transforms -> no explicit 1/beta (pinned by
        test_sigma_gate_one_iteration)."""
        axF, axB = self._ir_axF, self._ir_axB
        nx, ny, nz = self.lattice.shape
        nvol = self.lattice.nvol
        nblock = green_kw.shape[0]
        nd_block = green_kw.shape[-1]
        nd_v = v_eff.shape[-1]
        xp = _bk.array_module_of(green_kw)
        workers = getattr(self, "fft_workers", 1)
        ntF = axF.n_tau

        # G -> (r, tau_F)
        g = xp.moveaxis(
            green_kw.reshape(nblock, axF.n_freq, nvol * nd_block ** 2),
            1, -1)
        g_tau = xp.moveaxis(axF.freq_to_tau(g), -1, 1).reshape(
            nblock, ntF, nx, ny, nz, nd_block ** 2)
        green_rt = _bk.spatial_ifftn(g_tau, axes=(2, 3, 4), workers=workers
                                     ).reshape(nblock, ntF, nvol,
                                               nd_block, nd_block)

        # V_eff -> (r, tau_F) via the bosonic basis evaluated on tau_F
        v = xp.moveaxis(
            v_eff.reshape(axB.n_freq, nvol * nd_v * nd_v), 0, -1)
        v_tau = xp.moveaxis(
            axB.freq_to_tau_points(v, axF.tau), -1, 0).reshape(
            ntF, nx, ny, nz, nd_v * nd_v)
        v_rt = _bk.spatial_ifftn(v_tau, axes=(1, 2, 3), workers=workers
                                 ).reshape(ntF, nvol, nd_v, nd_v)

        # Sigma(r,tau) = V(r,tau) * G(r,tau), spin-slot sliced as uniform
        if nd_block != nd_v:
            norb = self.norb
            sigma_rt = xp.zeros((nblock, ntF, nvol, nd_block, nd_block),
                                dtype=np.complex128)
            for gblk in range(nblock):
                s = gblk if nblock == self.ns else 0
                sl = slice(s * norb, (s + 1) * norb)
                sigma_rt[gblk] = v_rt[:, :, sl, sl] * green_rt[gblk]
        else:
            sigma_rt = v_rt[np.newaxis] * green_rt

        nd_sig = sigma_rt.shape[-1]
        sigma_kt = _bk.spatial_fftn(
            sigma_rt.reshape(nblock, ntF, nx, ny, nz, nd_sig ** 2),
            axes=(2, 3, 4), workers=workers)
        sigma_kw = axF.tau_to_freq(xp.moveaxis(
            sigma_kt.reshape(nblock, ntF, nvol * nd_sig ** 2), 1, -1))
        return xp.moveaxis(sigma_kw, -1, 1).reshape(
            nblock, axF.n_freq, nvol, nd_sig, nd_sig)

    @do_profile
    def _calc_self_energy_general_ir(self, green_kw, v_eff, beta):
        r"""Sigma on the fermionic IR nodes, general (MYO) scheme: the IR-node
        transport of `_calc_self_energy_ir` with the rank-4 orbital contraction
        `_sigma_orbital_contract` of `_calc_self_energy_general` (Hadamard ->
        orbital sum). Pure-orbital (nd_block == nd_v == norb), no spin slice.
        Physical transforms -> no explicit 1/beta.
        """
        axF, axB = self._ir_axF, self._ir_axB
        nx, ny, nz = self.lattice.shape
        nvol = self.lattice.nvol
        nblock = green_kw.shape[0]
        norb = green_kw.shape[-1]
        ndx = v_eff.shape[-1]           # norb ** 2
        xp = _bk.array_module_of(green_kw)
        workers = getattr(self, "fft_workers", 1)
        ntF = axF.n_tau

        # G -> (r, tau_F)
        g = xp.moveaxis(
            green_kw.reshape(nblock, axF.n_freq, nvol * norb * norb), 1, -1)
        g_tau = xp.moveaxis(axF.freq_to_tau(g), -1, 1).reshape(
            nblock, ntF, nx, ny, nz, norb * norb)
        green_rt = _bk.spatial_ifftn(g_tau, axes=(2, 3, 4), workers=workers
                                     ).reshape(nblock, ntF, nvol, norb, norb)

        # V_eff -> (r, tau_F) via the bosonic basis evaluated on tau_F
        v = xp.moveaxis(v_eff.reshape(axB.n_freq, nvol * ndx * ndx), 0, -1)
        v_tau = xp.moveaxis(
            axB.freq_to_tau_points(v, axF.tau), -1, 0).reshape(
            ntF, nx, ny, nz, ndx * ndx)
        v_rt = _bk.spatial_ifftn(v_tau, axes=(1, 2, 3), workers=workers
                                 ).reshape(ntF, nvol, ndx, ndx)

        # Sigma(r,tau) = rank-4 orbital contraction (general)
        sigma_rt = self._sigma_orbital_contract(v_rt, green_rt)

        sigma_kt = _bk.spatial_fftn(
            sigma_rt.reshape(nblock, ntF, nx, ny, nz, norb * norb),
            axes=(2, 3, 4), workers=workers)
        sigma_kw = axF.tau_to_freq(xp.moveaxis(
            sigma_kt.reshape(nblock, ntF, nvol * norb * norb), 1, -1))
        return xp.moveaxis(sigma_kw, -1, 1).reshape(
            nblock, axF.n_freq, nvol, norb, norb)

    def _number_from_eigs_ir(self, lam, mu):
        """N(mu) (and dN/dmu) on the IR path: the k-summed trace of G at the
        fermionic nodes has closed form through the eigenvalues lam of M,
        and the particle number is the beta^- evaluation of its fermionic
        fit (n_total = -Re Tr G(beta^-)); the derivative passes through the
        same LINEAR fit (d/dmu of 1/(lam+mu) per node)."""
        axF = self._ir_axF
        denom = lam + mu
        g = _bk.to_host((1.0 / denom).sum(axis=(0, 2, 3)))
        dg = _bk.to_host((-1.0 / denom ** 2).sum(axis=(0, 2, 3)))
        n = -float(np.real(axF.fit_from_freq(g) @ axF.u_beta_minus))
        dn = -float(np.real(axF.fit_from_freq(dg) @ axF.u_beta_minus))
        return n, dn

    @do_profile
    def _ir_densify(self, arr, ax, freq_axis):
        """Evaluate a node-resolved array back onto the run's uniform grid
        (output compatibility; the frequency axis is ``freq_axis``).

        The uniform-grid OUTPUT is Nmat/L times larger than the node array
        (GB-scale at production sizes), so only the small coefficient array
        crosses the device boundary; the expansion runs as one host GEMM per
        leading index, written straight into the preallocated output buffer,
        so peak extra memory beyond that buffer is the coefficient array
        (design R-4: memory-aware densification)."""
        xp = _bk.array_module_of(arr)
        a = xp.ascontiguousarray(xp.moveaxis(arr, freq_axis, -1))
        coeffs = _bk.to_host(ax.fit_from_freq(a))       # (..., L): small
        _, ev = ax.uniform_matrices(self.nmat)          # (L, nmat)
        evT = np.ascontiguousarray(ev.T)                # (nmat, L)
        out = np.empty(arr.shape[:freq_axis] + (self.nmat,)
                       + arr.shape[freq_axis + 1:], dtype=np.complex128)
        # Write DIRECTLY in the output layout: for each leading index before
        # the frequency axis, out[idx] viewed as (nmat, trailing) is a
        # contiguous matrix filled by one GEMM (the transposed coefficient
        # operand is handled natively by BLAS) -- no GB-scale axis-move copy.
        lead_shape = arr.shape[:freq_axis]
        c = coeffs.reshape(lead_shape + (-1, ax.L))     # (lead..., trail, L)
        for idx in np.ndindex(lead_shape):
            np.matmul(evT, c[idx].T, out=out[idx].reshape(self.nmat, -1))
        return out

    def _ir_coeff_decay_check(self, sigma, label="sigma"):
        """Always-on layer-1 diagnostic (design Sec. 5): warn when the tail
        of the fitted coefficients stops decaying (out-of-basis content)."""
        axF = self._ir_axF
        c = _bk.to_host(axF.fit_from_freq(
            np.moveaxis(_bk.to_host(sigma), 1, -1)))
        mag = np.abs(c).max(axis=tuple(range(c.ndim - 1)))
        tail = float(mag[-max(1, axF.L // 10):].max())
        peak = float(mag.max()) or 1.0
        if tail > 10.0 * axF.eps * peak:
            logger.warning(
                "IR coefficient tail of %s is not decaying (ratio %.2e > "
                "10*ir_tol): the object exceeds the basis bandwidth -- "
                "raise [mode.param] ir_wmax or tighten [mode.param] ir_tol.",
                label, tail / peak)

    def _ir_sigma_init(self, seed):
        """Uniform-grid sigma_init -> fermionic nodes (uniform -> IR
        migration). The fit residual over the full uniform grid is checked
        against 100*ir_tol (relative); above it, behavior follows
        [mode.param] sigma_init_on_error ('warn' default / 'abort' /
        'zero'). Returns None to request the zero start."""
        axF = self._ir_axF
        a = np.moveaxis(seed, 1, -1)
        coeffs = axF.fit_from_uniform(a, self.nmat)
        resid = float(np.abs(axF.eval_to_uniform(coeffs, self.nmat)
                             - a).max())
        scale = float(np.abs(a).max()) or 1.0
        rel = resid / scale
        logger.info("IR sigma_init fit: max uniform residual %.3e (rel "
                    "%.3e)", resid, rel)
        if not self._sigma_seed_residual_ok(rel):
            return None
        return np.ascontiguousarray(
            np.moveaxis(axF.eval_to_freq(coeffs), -1, 1))

    def _sigma_seed_residual_ok(self, rel):
        """Shared sigma_init residual policy: True = use the fitted seed,
        False = fall back to the zero start; 'abort' raises."""
        if rel <= 100.0 * self._ir_axF.eps:
            return True
        msg = ("sigma_init IR fit residual (rel {:.3e}) exceeds "
               "100*ir_tol; the seed may exceed the basis "
               "bandwidth.".format(rel))
        if self.sigma_init_on_error == "abort":
            raise ValueError(
                msg + " Raise [mode.param] ir_wmax, or set [mode.param] "
                "sigma_init_on_error='warn'/'zero'.")
        if self.sigma_init_on_error == "zero":
            logger.warning(
                "%s Falling back to the zero start "
                "(sigma_init_on_error='zero').", msg)
            return False
        logger.warning(
            "%s Using the fitted seed anyway "
            "(sigma_init_on_error='warn').", msg)
        return True

    def _ir_sigma_init_native(self, seed, meta, beta):
        """IR-native sigma_init -> this run's fermionic nodes (Stage 3,
        design Sec. 4.2). Cross-temperature seeding is DELIBERATE (sweep
        chains): the file's integer node indices are interpreted on the
        run's beta (index-matched rescaling, mirroring the uniform warm
        start), logged when the betas differ. Returns None to request the
        zero start (residual policy)."""
        axF = self._ir_axF
        file_beta = float(meta["beta"])
        if not np.isclose(file_beta, beta, rtol=1e-9, atol=1e-9 * beta):
            logger.info(
                "sigma_init: seed ir_beta %.12g differs from this run's "
                "beta %.12g (temperature-sweep warm start); node values "
                "are consumed index-matched on the run's frequencies.",
                file_beta, beta)
        freq_n = np.asarray(meta["freq_n"], dtype=np.int64)
        if np.array_equal(freq_n, axF.freq_n):
            return np.ascontiguousarray(seed)      # already run-node values
        # Fit in the WRITER's basis (reconstructed from the stored params)
        # at its own nodes, then evaluate at THIS run's node indices. A fit
        # onto the run basis would be underdetermined for every
        # high-T -> low-T sweep step (the run L grows with beta while the
        # file only carries its own ~L_file nodes); the writer-basis fit is
        # determined by construction and realizes the index-matched
        # semantics literally (the seed as a function of n, evaluated at
        # the new n's).
        from hwave.solver.ir_axis import IRAxis
        ax_file = IRAxis(beta=file_beta, wmax=meta["wmax"], eps=meta["tol"],
                         statistics="F")
        a = np.moveaxis(seed, 1, -1)
        coeffs = ax_file.fit_from_freq_points(a, freq_n)
        resid = float(np.abs(
            ax_file.eval_to_freq_points(coeffs, freq_n) - a).max())
        scale = float(np.abs(a).max()) or 1.0
        rel = resid / scale
        logger.info("IR-native sigma_init fit (writer basis): max residual "
                    "at file nodes %.3e (rel %.3e)", resid, rel)
        if not self._sigma_seed_residual_ok(rel):
            return None
        out = ax_file.eval_to_freq_points(coeffs, axF.freq_n)
        return np.ascontiguousarray(np.moveaxis(out, -1, 1))

    @do_profile
    def _calc_dressed_green(self, beta, mu, sigma):
        """Compute dressed Green's function G(k, iwn) = [G0^{-1} - Sigma]^{-1}.

        Parameters
        ----------
        beta : float
            Inverse temperature.
        mu : float
            Chemical potential.
        sigma : ndarray
            Self-energy, shape (nblock, nmat, nvol, nd, nd).

        Returns
        -------
        ndarray
            Dressed Green's function, shape (nblock, nmat, nvol, nd, nd).
        """
        logger.debug(">>> FLEX._calc_dressed_green")

        ew = self.H0_eigenvalue
        ev = self.H0_eigenvector
        xp = _bk.array_module_of(sigma)

        nblock, nvol, nd = ew.shape

        # Matsubara frequencies of the working axis (uniform or IR nodes)
        iomega = self._freq_omegas(beta, xp)

        # Reconstruct H0 in orbital basis from eigendecomposition
        # H0 = ev @ diag(ew) @ ev†, using matmul for BLAS efficiency
        H0_k = xp.matmul(ev * ew[:, :, np.newaxis, :], xp.conj(ev).swapaxes(-2, -1))

        # G^{-1}(k, iwn) = (iwn + mu) * I - H0(k) - Sigma(k, iwn)
        eye = xp.eye(nd, dtype=np.complex128)

        # Vectorized construction of G^{-1} for all frequencies
        # iomega shape: (nmat,) -> broadcast to (1, nmat, 1, 1, 1)
        iw = 1j * iomega[np.newaxis, :, np.newaxis, np.newaxis, np.newaxis]
        # H0_k shape: (nblock, nvol, nd, nd) -> (nblock, 1, nvol, nd, nd)
        H0_exp = H0_k[:, np.newaxis, :, :, :]

        green_inv = (iw + mu) * eye - H0_exp - sigma

        # G(k, iwn) = [G^{-1}]^{-1}
        green = xp.linalg.inv(green_inv)

        return green

    @do_profile
    def _calc_number_dressed(self, sigma, mu, beta):
        """Particle number carried by the dressed Green's function at ``mu``.

        The self-energy shifts and renormalizes the dressed G, so the actual
        electron count differs from the non-interacting Fermi count.  The count
        is evaluated with the standard tail-subtracted Matsubara sum: the
        non-interacting reference is summed analytically (Fermi function) and
        only the small, absolutely-convergent difference ``G - G0`` is summed
        over the finite Matsubara grid::

            N(mu) = sum_{block,k,a} f(eps_a(k) - mu)
                    + (1/beta) sum_{k,n} Tr[ G(k,iwn) - G0(k,iwn) ]

        Here ``G0(k,iwn) = [(iwn+mu)I - H0(k)]^{-1}`` (same mu), whose analytic
        occupation is exactly the Fermi term, and ``G - G0 = G0 Sigma G`` decays
        as ``1/(iwn)^2`` so the residual sum needs no e^{iwn 0+} convergence
        factor and is real.  At ``Sigma = 0`` the residual vanishes and ``N``
        reduces exactly to the ``_find_mu`` Fermi count.

        The total is compared against the SAME target convention as
        :meth:`RPA._find_mu`: the sum runs over all spin blocks, so for
        spin-free (one block) it is the one-spin count (target ``Ncond/2``).

        Parameters
        ----------
        sigma : ndarray
            Self-energy, shape (nblock, nmat, nvol, nd_block, nd_block).
        mu : float
            Trial chemical potential.
        beta : float
            Inverse temperature.

        Returns
        -------
        float
            Particle number N(mu).
        """
        nmat = self.nmat
        ew = self.H0_eigenvalue                       # (nblock, nvol, nd_block)
        xp = _bk.array_module_of(sigma)

        # analytic non-interacting reference (Fermi function)
        n_ref = self._fermi_occupation(1.0 / beta, mu, ew).sum()

        # dressed correction Tr[G - G0], summed over the finite Matsubara grid
        green = self._calc_dressed_green(beta, mu, sigma)
        iomega = (xp.arange(nmat) * 2 + 1 - nmat) * np.pi / beta
        trG = xp.einsum('bnkaa->bnk', green)          # (nblock, nmat, nvol)
        trG0 = (1.0 / ((1j * iomega)[np.newaxis, :, np.newaxis, np.newaxis]
                       + (mu - ew)[:, np.newaxis, :, :])).sum(axis=-1)
        corr = (trG - trG0).sum() / beta

        return float(n_ref + corr.real)

    @staticmethod
    def _fermi_occupation(t, mu, ev, ene_cutoff=1.0e2):
        """Fermi function with the same overflow guard as RPA._find_mu.
        Must stay arithmetically identical to rpa.py's module-level
        _masked_fermi_delta_n and to _find_mu's internal _fermi closure
        (both in src/hwave/solver/rpa.py) -- three copies of the same
        masked Fermi factor."""
        xp = _bk.array_module_of(ev)
        w = (ev - mu) / t
        mask = w < ene_cutoff
        w1 = xp.where(mask, w, 0.0)
        v1 = 1.0 / (1.0 + xp.exp(w1))
        return xp.where(mask, v1, 0.0)

    @do_profile
    def _calc_occupation_dressed(self, green_kw, mu, beta):
        r"""k-summed occupation per (spin block, orbital) from the dressed G.

        The per-orbital analogue of :meth:`_calc_number_dressed`: the same
        tail-subtracted Matsubara sum, resolved on each orbital ``a`` instead of
        traced.  The non-interacting reference is projected onto the orbital
        basis with the H0 eigenvectors ``V``::

            n_a(k) = sum_j |V_aj|^2 f(eps_j(k) - mu)
                     + (1/beta) sum_n [ G_aa(k,iwn) - G0_aa(k,iwn) ]
            G0_aa(k,iwn) = sum_j |V_aj|^2 / (iwn + mu - eps_j(k))

        Summing the result over orbitals recovers ``_calc_number_dressed`` (by
        unitarity ``sum_a |V_aj|^2 = 1``); it is split per orbital here so N and
        Sz can be assembled from spin-resolved partial sums.

        Parameters
        ----------
        green_kw : ndarray
            Dressed Green's function, shape (nblock, nmat, nvol, nd_block, nd_block).
        mu : float
            Chemical potential.
        beta : float
            Inverse temperature.

        Returns
        -------
        ndarray
            Occupation summed over k, shape (nblock, nd_block).
        """
        xp = _bk.array_module_of(green_kw)
        if self._ir_axF is not None:
            # IR path: n_a = -Re G_aa(tau = beta^-) directly through the
            # fermionic basis (tail-free by construction).
            axF = self._ir_axF
            g_aa = xp.einsum('bnkaa->bnka', green_kw)   # (nb, nw, nvol, a)
            g_aa = _bk.to_host(xp.moveaxis(g_aa, 1, -1))  # (nb, nvol, a, nw)
            n_k = -np.real(axF.fit_from_freq(g_aa) @ axF.u_beta_minus)
            return n_k.sum(axis=1)                       # (nblock, nd_block)

        nmat = self.nmat
        ew = self.H0_eigenvalue                       # (nblock, nvol, nd_block)
        ev = self.H0_eigenvector                      # (nblock, nvol, a, j)
        vsq = xp.abs(ev) ** 2                          # |V_aj|^2

        # bare per-orbital reference: n0_a = sum_j |V_aj|^2 f(eps_j - mu)
        f = self._fermi_occupation(1.0 / beta, mu, ew)     # (nblock, nvol, j)
        n0 = xp.einsum('bkaj,bkj->bka', vsq, f)            # (nblock, nvol, a)

        # dressed correction (1/beta) sum_n [G_aa - G0_aa]
        iomega = (xp.arange(nmat) * 2 + 1 - nmat) * np.pi / beta
        g_aa = xp.einsum('bnkaa->bnka', green_kw)          # (nblock, nmat, nvol, a)
        denom = ((1j * iomega)[np.newaxis, :, np.newaxis, np.newaxis]
                 + (mu - ew)[:, np.newaxis, :, :])         # (nblock, nmat, nvol, j)
        g0_aa = xp.einsum('bkaj,bnkj->bnka', vsq, 1.0 / denom)
        corr = (g_aa - g0_aa).sum(axis=1) / beta           # sum over n -> (nblock, nvol, a)

        nocc = n0 + corr.real
        return nocc.sum(axis=1)                            # sum over k -> (nblock, nd_block)

    @do_profile
    def _calc_physics_dressed(self, green_kw, mu, beta):
        """Assemble the particle number N and spin Sz from the dressed G.

        Uses :meth:`_calc_occupation_dressed` and the FLEX spin-block orbital
        ordering (``s*norb + a``).  Conventions per spin mode:

        - spin-free (one block computed): ``N = 2 * sum(nocc)``, ``Sz = 0``.
        - spin-diag (block 0 = up, block 1 = down):
          ``N = N_up + N_down``, ``Sz = (N_up - N_down)/2``.
        - spinful (single block, ``nd = 2*norb``, up = orbitals ``[0, norb)``,
          down = ``[norb, 2*norb)``): ``N = sum(nocc)``,
          ``Sz = (N_up - N_down)/2``.

        Parameters
        ----------
        green_kw : ndarray
            Dressed Green's function (nblock, nmat, nvol, nd_block, nd_block).
        mu : float
            Chemical potential.
        beta : float
            Inverse temperature.

        Returns
        -------
        dict
            ``{"NCond": N_total, "Sz": Sz, "mu": mu}`` (plain floats).
        """
        nocc = self._calc_occupation_dressed(green_kw, mu, beta)
        norb = self.norb

        if self.spin_mode == "spin-free":
            n_total = 2.0 * nocc.sum()
            sz = 0.0
        elif self.spin_mode == "spin-diag":
            n_up = nocc[0].sum()
            n_down = nocc[1].sum()
            n_total = n_up + n_down
            sz = 0.5 * (n_up - n_down)
        else:  # spinful: spin folded into the orbital index (spin-block order)
            n_up = nocc[0, :norb].sum()
            n_down = nocc[0, norb:].sum()
            n_total = nocc.sum()
            sz = 0.5 * (n_up - n_down)

        return {"NCond": float(n_total.real if np.iscomplexobj(n_total)
                               else n_total),
                "Sz": float(sz.real if np.iscomplexobj(sz) else sz),
                "mu": float(mu)}

    @do_profile
    def _matsubara_number_operator(self, sigma, beta):
        r"""Eigenvalues of the mu-independent part of ``G^{-1}``, for the mu
        search.

        During the chemical-potential search ``sigma`` is held FIXED, so

            G^{-1}(k, iwn; mu) = (iwn + mu) I - H0(k) - Sigma(k, iwn)
                               = M(k, iwn) + mu I,
            M(k, iwn) = iwn I - H0(k) - Sigma(k, iwn)   (mu-independent).

        Diagonalizing ``M`` ONCE then gives ``Tr[G(mu)] = sum_j 1/(lam_j + mu)``
        for ANY trial ``mu`` -- one eigenvalue decomposition per iteration
        instead of one full matrix inversion per bisection step (the mu search
        was ~50% of solve()).  ``M`` is non-Hermitian, but its eigenvalues sit
        at ``lam_j ~ iwn - (real band+Sigma)``, i.e. ``|Im(lam_j)| ~ |wn| >=
        pi/beta > 0``, so ``lam_j + mu`` (mu real) never touches zero: the
        1/(lam+mu) sum is well conditioned.

        Parameters
        ----------
        sigma : ndarray
            Self-energy, shape (nblock, nmat, nvol, nd_block, nd_block).
        beta : float
            Inverse temperature.

        Returns
        -------
        lam : ndarray
            Eigenvalues of ``M``, shape (nblock, nmat, nvol, nd_block).
        ew : ndarray
            H0 band energies ``self.H0_eigenvalue`` (returned for convenience,
            used by :meth:`_number_from_eigs` for the analytic reference).
        """
        ew = self.H0_eigenvalue
        ev = self.H0_eigenvector
        xp = _bk.array_module_of(sigma)
        nblock, nvol, nd = ew.shape

        iomega = self._freq_omegas(beta, xp)
        H0_k = xp.matmul(ev * ew[:, :, np.newaxis, :],
                         xp.conj(ev).swapaxes(-2, -1))       # (nb, nvol, nd, nd)
        eye = xp.eye(nd, dtype=np.complex128)
        iw = 1j * iomega[np.newaxis, :, np.newaxis, np.newaxis, np.newaxis]
        M = iw * eye - H0_k[:, np.newaxis, :, :, :] - sigma  # (nb,nmat,nvol,nd,nd)

        # M is non-Hermitian, but for nd <= 2 the eigenvalues have a closed
        # form evaluated on M's own backend (_eigvals_small), so the mu search
        # stays on the GPU end-to-end. Only nd >= 3 needs LAPACK geev, which
        # CuPy does not provide -- that case runs on the host, and ew is
        # returned on the matching backend so _number_from_eigs never mixes
        # host and device arrays.
        lam = _eigvals_small(M)                               # (nb,nmat,nvol,nd)
        if _bk.array_module_of(lam) is np:
            return lam, _bk.to_host(ew)
        return lam, ew

    def _number_from_eigs(self, lam, ew, mu, beta, with_deriv=False):
        """Particle number N(mu) (and optionally dN/dmu) from precomputed
        ``M`` eigenvalues.

        Cheap closure evaluated for every trial ``mu`` in the mu search:
        identical formula and result as :meth:`_calc_number_dressed` (verified
        to machine precision in the tests), but reusing the eigenvalues from
        :meth:`_matsubara_number_operator` instead of re-inverting G::

            N(mu) = sum_{block,k,a} f(eps_a - mu)
                    + (1/beta) sum_{k,n} [ sum_j 1/(lam_j + mu)
                                           - sum_a 1/(iwn + mu - eps_a) ]

        Because ``sigma`` (hence ``lam``) is FIXED during the search, the
        derivative is analytic and cheap to evaluate from the same
        eigenvalues, which lets the search use Newton's method::

            dN/dmu = sum_a (1/T) f_a (1 - f_a)
                     + (1/beta) sum_{k,n} [ -sum_j 1/(lam_j + mu)^2
                                            + sum_a 1/(iwn + mu - eps_a)^2 ]

        Parameters
        ----------
        lam : ndarray
            ``M`` eigenvalues from :meth:`_matsubara_number_operator`,
            shape (nblock, nmat, nvol, nd_block).
        ew : ndarray
            H0 band energies, shape (nblock, nvol, nd_block).
        mu : float
            Trial chemical potential.
        beta : float
            Inverse temperature.
        with_deriv : bool, optional
            If True, return the tuple ``(N, dN/dmu)`` instead of ``N``.

        Returns
        -------
        float or tuple of float
            ``N(mu)``, or ``(N(mu), dN/dmu)`` when ``with_deriv`` is True.
        """
        nmat = self.nmat
        T = 1.0 / beta
        xp = _bk.array_module_of(lam)
        iomega = (xp.arange(nmat) * 2 + 1 - nmat) * np.pi / beta

        f = self._fermi_occupation(T, mu, ew)
        n_ref = f.sum()

        denomG = lam + mu                                    # (nb, nmat, nvol, nd)
        denomG0 = ((1j * iomega)[np.newaxis, :, np.newaxis, np.newaxis]
                   + (mu - ew)[:, np.newaxis, :, :])
        trG = (1.0 / denomG).sum(axis=-1)                    # (nb, nmat, nvol)
        trG0 = (1.0 / denomG0).sum(axis=-1)
        # plain Python floats: the safeguarded-Newton loop does host-side
        # comparisons on every trial mu, so device scalars would sync anyway.
        n = float(n_ref + ((trG - trG0).sum() / beta).real)

        if not with_deriv:
            return n

        # dN/dmu: the Fermi part is +f(1-f)/T; each 1/(x+mu) term contributes
        # -1/(x+mu)^2 (same trG - trG0 combination, both differentiated).
        dref = ((1.0 / T) * f * (1.0 - f)).sum()
        dtrG = (-1.0 / denomG ** 2).sum(axis=-1)
        dtrG0 = (-1.0 / denomG0 ** 2).sum(axis=-1)
        dn = float(dref + ((dtrG - dtrG0).sum() / beta).real)
        return n, dn

    @do_profile
    def _find_mu_dressed(self, sigma, beta, Ncond):
        """Re-solve mu so the dressed Green's function carries ``Ncond``.

        Holds ``sigma`` fixed and solves ``N(mu) = Ncond``.  The ``M``
        eigenvalues are computed ONCE via :meth:`_matsubara_number_operator``,
        so both ``N(mu)`` and its analytic derivative ``dN/dmu`` are cheap sums
        over those eigenvalues (:meth:`_number_from_eigs`) -- no per-step matrix
        inversion.  ``N(mu)`` is smooth and monotonically increasing, so the
        root is found with a **safeguarded Newton** iteration (Numerical
        Recipes ``rtsafe``): a Newton step when it stays inside the current
        bracket and makes progress, otherwise a bisection step.  This gets
        Newton's quadratic convergence (a handful of evaluations from the
        previous iteration's mu) while the maintained bracket guarantees
        convergence even if a Newton step misbehaves.

        The initial bracket is the bare eigenvalue range widened by the
        self-energy scale, because Sigma can push spectral weight outside the
        bare band edges.

        Parameters
        ----------
        sigma : ndarray
            Current self-energy (fixed during the mu search).
        beta : float
            Inverse temperature.
        Ncond : float
            Target particle number (already halved for spin-free, matching
            :meth:`RPA._find_mu`).

        Returns
        -------
        float
            Chemical potential mu such that N(mu) = Ncond.
        """
        # Diagonalize the mu-independent part of G^{-1} once; each trial mu is
        # then a cheap sum over eigenvalues (see _matsubara_number_operator).
        # lam/ew come back on a consistent backend (device for nd<=2 under
        # GPU execution, host otherwise); _number_from_eigs dispatches on it
        # and returns plain floats, so the Newton/bisection logic below is
        # backend-independent.
        lam, ew = self._matsubara_number_operator(sigma, beta)
        w = ew

        def _delta_n(mu, with_deriv=False):
            if self._ir_axF is not None:
                n, dn = self._number_from_eigs_ir(lam, mu)
                if with_deriv:
                    return n - Ncond, dn
                return n - Ncond
            r = self._number_from_eigs(lam, ew, mu, beta, with_deriv=with_deriv)
            if with_deriv:
                return r[0] - Ncond, r[1]
            return r - Ncond

        # Initial bracket: the bare band range widened by the self-energy scale.
        # This is only a starting guess -- a general (multi-orbital, off-diagonal
        # or non-Hermitian) Sigma can shift the dressed spectral weight, and
        # hence the root, by more than max|Re Sigma|.  So EXPAND the bracket
        # (doubling the pad) until N(mu) changes sign, instead of failing on the
        # first guess.  N(mu) -> 0 as mu -> -inf and -> Nstate as mu -> +inf, so
        # a finite root is always bracketed after enough expansion.
        #
        # NOTE: do NOT warm-start this search from the previous iteration's mu
        # with a narrow bracket. With a large (not-yet-converged) Sigma the
        # truncated-Matsubara N(mu) is not strictly monotone, so a narrow
        # bracket can select a DIFFERENT root than the wide band-range bracket
        # and silently change the SCF trajectory (observed at 64^2, Nmat=1024:
        # the warm-started run diverged while the cold-bracket run was stable).
        pad = float(_bk.array_module_of(sigma).abs(sigma.real).max()) + 1.0
        lo = float(w.min()) - pad
        hi = float(w.max()) + pad
        f_lo = _delta_n(lo)
        f_hi = _delta_n(hi)
        for _ in range(60):
            if f_lo == 0.0:
                return lo
            if f_hi == 0.0:
                return hi
            if f_lo * f_hi < 0.0:
                break
            span = hi - lo
            lo -= span
            hi += span
            f_lo = _delta_n(lo)
            f_hi = _delta_n(hi)
        else:
            # 60 doublings span ~1e18 * initial width: a real root cannot be
            # this far out. Fail loudly rather than return garbage. Raise a
            # catchable exception (not sys.exit) so parameter sweeps, pipelines,
            # and notebooks can handle/aggregate the failure instead of having
            # the whole interpreter torn down.
            raise RuntimeError(
                "FLEX._find_mu_dressed: chemical-potential root not bracketed "
                "after expansion to [{}, {}] (N-Ncond = {}, {}). Check the "
                "target filling/Ncond and temperature.".format(
                    lo, hi, f_lo, f_hi))

        # orient the bracket so f(lo) < 0 < f(hi) (N increases with mu)
        if f_lo > 0.0:
            lo, hi = hi, lo

        mu = 0.5 * (lo + hi)
        dmu_old = abs(hi - lo)
        dmu = dmu_old
        f, df = _delta_n(mu, with_deriv=True)

        xtol = 1.0e-12
        for _ in range(100):
            # take bisection when the Newton step leaves the bracket, is not
            # shrinking the interval fast enough, or the derivative is flat;
            # otherwise take Newton.  df -> 0 is the CHARGE-GAP regime: N(mu) is
            # flat across the gap, so a Newton step is unreliable (and f/df ill
            # defined).  The rtsafe product test below already routes tiny df to
            # bisection, but guard df == 0 explicitly so f/df is never formed.
            newton_out = ((mu - hi) * df - f) * ((mu - lo) * df - f) > 0.0
            slow = abs(2.0 * f) > abs(dmu_old * df)
            if newton_out or slow or df == 0.0:
                dmu_old = dmu
                dmu = 0.5 * (hi - lo)
                mu = lo + dmu
            else:
                dmu_old = dmu
                dmu = f / df
                mu = mu - dmu

            if abs(dmu) < xtol:
                break

            f, df = _delta_n(mu, with_deriv=True)
            # keep the sign-change bracket [lo, hi] around the root
            if f < 0.0:
                lo = mu
            else:
                hi = mu

        logger.info("FLEX._find_mu_dressed: mu = {}".format(mu))
        return mu

    @do_profile
    def _flex_compute_veff(self, chi0q_raw, ham_orig):
        """Inflate chi0q, decompose into spin/charge channels, and compute V_eff.

        This method handles the spin inflation (same as RPA.solve()) and then
        decomposes the interaction into spin and charge channels for FLEX.

        Parameters
        ----------
        chi0q_raw : ndarray
            Raw chi0q from _calc_chi0q (before spin inflation).
        ham_orig : ndarray
            Original interaction Hamiltonian in full spin-orbital space.

        Returns
        -------
        chi0q_inflated : ndarray
            Chi0q after spin inflation (for output).
        v_eff : ndarray
            Effective FLEX interaction V_eff(q, ivn).
        chi_s : ndarray
            Spin susceptibility.
        chi_c : ndarray
            Charge susceptibility.
        """
        nvol = self.lattice.nvol
        norb = self.norb
        ns = self.ns
        nd = self.nd

        # Inflate chi0q to reduced spin-orbital space (same logic as RPA.solve())
        chi0q, ham = self._inflate_chi0q_and_ham(chi0q_raw, ham_orig)

        # Decompose ham into spin and charge channels
        ham_s, ham_c = self._build_spin_charge_vertices(ham)

        # Solve RPA for spin channel: chi_s = [1 - chi0*V_s]^{-1} chi0
        chi_s = self._solve_rpa(chi0q, ham_s)

        # Solve RPA for charge channel: chi_c = [1 + chi0*V_c]^{-1} chi0
        chi_c = self._solve_rpa(chi0q, ham_c)

        # Compute V_eff
        v_eff = self._calc_veff(chi0q, chi_s, chi_c, ham)

        return chi0q, v_eff, chi_s, chi_c

    @do_profile
    def _flex_compute_veff_general(self, chi0q_raw, ham_orig):
        """General (paramagnetic full-vertex) counterpart of _flex_compute_veff.

        Mirrors the return signature ``(chi0q_out, v_eff, chi_s, chi_c)`` of
        :meth:`_flex_compute_veff`, but uses the general-path methods that keep
        the full rank-4 orbital vertex (MYO-convention S/C matrices) instead of
        the density-density reduction.

        Parameters
        ----------
        chi0q_raw : ndarray
            Raw rank-6 chi0q ``(nmat, nvol, norb, norb, norb, norb)`` from
            :meth:`_calc_chi0q` after the spin block dimension has been stripped.
        ham_orig : ndarray
            Original interaction Hamiltonian in full spin-orbital space.

        Returns
        -------
        chi0q_out : ndarray
            chi0q for public output, in the RPA ``[a,c,b,d]`` orbital convention
            (rank-6), consistent with the reduced path; this is also the layout
            used internally for the S/C math.
        v_eff : ndarray
            Effective FLEX interaction ``(nmat, nvol, norb^2, norb^2)``.
        chi_s : ndarray
            Spin susceptibility for public output, in the RPA ``[a,c,b,d]``
            orbital convention (rank-6), consistent with chi0q_out and the
            reduced path; this is also the layout used internally for the S/C
            math.
        chi_c : ndarray
            Charge susceptibility for public output, in the RPA ``[a,c,b,d]``
            orbital convention (rank-6), consistent with chi0q_out and the
            reduced path; this is also the layout used internally for the S/C
            math.
        """
        chi0q, Us, Uc = self._inflate_chi0q_and_ham_general(chi0q_raw, ham_orig)
        chi_s, chi_c = self._solve_channels_general(chi0q, Us, Uc)
        v_eff = self._calc_veff_general(chi0q, chi_s, chi_c, Us, Uc)
        # chi0q / chi_s / chi_c are already in the native RPA [a,c,b,d]
        # orbital-pair convention (see _inflate_chi0q_and_ham_general), which is
        # what the public output and the Eliashberg loader expect: the loader
        # (_compute_vertices_flex) builds S @ chi @ S with S/C matrices in that
        # native layout (the two builders are one implementation since the #113
        # adjudication), so a transposed chi here would build a transposed
        # pairing vertex and a wrong static lambda (issue #78).  The "myo"
        # chi_convention tag marks these as general-path (orbital-pair shape).
        #
        # This used to transpose the pair axes on the way IN, which is what
        # issue #91 was: the transpose reached V_eff, so Sigma was built from
        # the transposed bubble -- Sigma_mn came out as chi0_nm G_mn instead of
        # chi0_mn G_mn, i.e. right on the orbital diagonal and wrong off it.
        # (Sigma itself is NOT simply transposed: G stays in its own slots.)  It is gone, so there is nothing left
        # to undo here either.
        #
        # Effect on the SAVED output, relative to the previous behaviour:
        #   chi0q  unchanged -- it was already transposed back at this boundary.
        #   chi_s / chi_c  unchanged.  Issue #78 had been fixed at this same
        #     boundary, by solving from the transposed chi0 and transposing the
        #     channels back: transpose(solve(transpose(chi0), U)).  Removing the
        #     transpose at its source gives solve(chi0, U) with no compensation,
        #     and the push-through identity makes the two equal -- so this
        #     retires the compensation rather than changing its result.  The
        #     identity needs S/C symmetric under the orbital-pair transpose,
        #     which holds whenever the on-site interaction parameters are
        #     symmetric in their orbital indices (U'_ab = U'_ba, J_ab = J_ba,
        #     ...): the physical case, and the only one the S/C construction is
        #     meant for.
        #   Sigma  changed off the orbital diagonal.  This is issue #91, and the
        #     part the output-boundary fix could not reach: V_eff is built from
        #     chi0 INSIDE this function, so a transposed chi0 reaches Sigma no
        #     matter what the boundary does.
        #
        # Note that since issue #93 the READERS reject non-Hermitian-closed
        # files; an asymmetric TABLE reaches here only via internal
        # construction. For such a table the
        # two forms differ, and solve(chi0, U) is the consistent value because
        # it applies the S/C matrices in the same index order that RPA and the
        # Eliashberg loader's S @ chi @ S use.
        return chi0q, v_eff, chi_s, chi_c

    @do_profile
    def _inflate_chi0q_and_ham(self, chi0q_raw, ham_orig):
        """Apply spin inflation to chi0q and ham, matching RPA.solve() logic.

        Takes raw chi0q (orbital-only for spin-free, spin-block for spin-diag,
        or spin-orbital for spinful) and inflates to the reduced space used
        for FLEX calculations.

        Parameters
        ----------
        chi0q_raw : ndarray
            Raw chi0q before spin inflation.
        ham_orig : ndarray
            Full interaction Hamiltonian.

        Returns
        -------
        chi0q : ndarray
            Inflated chi0q in reduced spin-orbital space.
        ham : ndarray
            Interaction Hamiltonian in matching reduced space.
        """
        nvol = self.lattice.nvol
        norb = self.norb
        ns = self.ns
        nd = self.nd
        xp = _bk.array_module_of(chi0q_raw)

        if self.spin_mode == "spin-free":
            # chi0q_raw shape: (nmat, nvol, norb, norb) for reduced
            nfreq = chi0q_raw.shape[0]

            # Inflate to spin-orbital reduced space.
            # Equivalent to a Kronecker product with I_ns: scatter chi0q onto
            # the spin-block diagonal (same spin block for row and col),
            # leaving off-diagonal spin blocks exactly zero. Bit-identical to
            # np.einsum('lkab,st->lksatb', chi0q, I_ns), but faster.
            chi0q_src = chi0q_raw.reshape(nfreq, nvol, norb, norb)
            chi0q = xp.zeros((nfreq, nvol, nd, nd), dtype=chi0q_src.dtype)
            for s in range(ns):
                sl = slice(s * norb, (s + 1) * norb)
                chi0q[..., sl, sl] = chi0q_src

            # 'ksasatbtb->ksatb' through the (ns, norb) factorized view
            # picks the same entries as the combined-index pair diagonal
            # and lands them in the same spin-major (nd, nd) layout, so it
            # routes through the shared projection
            ham = project_density_pairs(ham_orig, nvol, nd, xp)

        elif self.spin_mode == "spin-diag":
            # chi0q_raw shape: (nblock=2, nmat, nvol, norb, norb) for reduced
            nblock_s, nfreq, nvol_c, norb1, norb2 = chi0q_raw.shape

            # Inflate per-spin-block chi0q to spin-orbital reduced space.
            # Source block g goes to spin block (g, g) on the diagonal, with
            # off-diagonal spin blocks exactly zero. Bit-identical to
            # np.einsum('glkab,gh->lkgahb', chi0q, I_ns), but faster.
            chi0q = xp.zeros((nfreq, nvol, nd, nd), dtype=chi0q_raw.dtype)
            for s in range(ns):
                sl = slice(s * norb, (s + 1) * norb)
                chi0q[..., sl, sl] = chi0q_raw[s]

            # 'ksasatbtb->ksatb' through the (ns, norb) factorized view
            # picks the same entries as the combined-index pair diagonal
            # and lands them in the same spin-major (nd, nd) layout, so it
            # routes through the shared projection
            ham = project_density_pairs(ham_orig, nvol, nd, xp)

        elif self.spin_mode == "spinful":
            # chi0q_raw shape: (nmat, nvol, nd, nd) for reduced
            chi0q = chi0q_raw

            ham = project_density_pairs(ham_orig, nvol, nd, xp)

        return chi0q, ham

    @do_profile
    def _inflate_chi0q_and_ham_general(self, chi0q_raw, ham_orig):
        """Pass chi0q through; build the full MYO S/C interaction matrices.

        This is the paramagnetic (spin-free) full-vertex general-path analogue
        of :meth:`_inflate_chi0q_and_ham`.  Unlike the reduced path,
        which collapses the vertex onto the density-density diagonal and works
        in the spin-orbital ``nd = norb * ns`` space, the MYO paramagnetic
        formalism (cond-mat/0407094 Eqs. 5/6) is purely in ORBITAL space and
        keeps the full rank-4 vertex.  Hence the working dimension here is
        ``no = self.norb`` (NOT ``self.nd``): chi0q stays orbital-only and the
        S/C matrices are ``norb^2 x norb^2`` MYO Kanamori blocks.

        chi0q is passed through UNCHANGED in its orbital-pair index order.
        ``RPA._calc_chi0q`` labels the rank-6 bare bubble as
        ``chi0[..., a, c, b, d]`` (rpa.py: ``chi0[a,c,b,d] = G[a,b] . G_rev[d,c]``);
        viewed as a ``(norb^2, norb^2)`` matrix with row ``(a,c)`` and column
        ``(b,d)`` that is ``chi0_{(mn),(mu nu)}(q) = -(T/N) sum_k G_{m mu}(k+q)
        G_{nu n}(k)``, which is exactly the order consumed by the downstream
        general-path machinery -- the S/C matrices,
        :meth:`_solve_channels_general`, :meth:`_calc_veff_general` and
        :meth:`_sigma_orbital_contract`.

        Earlier revisions applied an "RPA -> MYO" orbital-pair transpose here
        (and undid it on the way out).  That was wrong (issue #91): the
        round-trip left the public susceptibilities untouched but transposed
        ``V_eff``, so Sigma was built from the transposed bubble and came out
        wrong off the orbital diagonal.  The
        index order is now pinned by a ground truth independent of this codebase
        -- the O(U^2) self-energy of an exactly-diagonalized 3-orbital model, see
        ``tests/test_flex_sopt_index_order.py``.

        Parameters
        ----------
        chi0q_raw : ndarray
            Bare susceptibility, shape (nmat,nvol,norb,norb,norb,norb)
            (6-dim: freq, vol, 4 orbital legs) from ``RPA._calc_chi0q`` for the
            spin-free general scheme, in RPA orbital-pair convention
            ``[a,c,b,d]``, after the (size-1) spin-block dimension has been
            stripped in ``solve()``.
        ham_orig : ignored
            Present for signature parity with ``_inflate_chi0q_and_ham``; the
            general path rebuilds the interaction from the real-space
            ``self.ham_info.param_ham`` rather than the pre-inflated
            ``ham_inter_q``.

        Returns
        -------
        chi0q : ndarray
            ``chi0q_raw`` unchanged (6-dim ``(nmat,nvol,m,n,mu,nu)``, i.e. the
            native RPA ``[a,c,b,d]`` orbital-pair order).
        Us : ndarray
            MYO spin (S) interaction matrices, shape ``(nvol, norb^2, norb^2)``.
        Uc : ndarray
            MYO charge (C) interaction matrices, shape ``(nvol, norb^2, norb^2)``.

        Notes
        -----
        The S/C matrices are built via ``build_sc_matrices_myo`` from an
        ``inter_k`` dict assembled with ``hwave.sc._build_interaction_k``.  For
        on-site Kanamori interactions the S/C matrices are CONSTANT over q, so
        the k-array ordering used to build ``inter_k`` (a plain linspace grid)
        need not match the FFT q-grid for this v1 path; this is reshaped to
        ``(nvol, norb^2, norb^2)`` purely as ``Nx*Ny*Nz`` independent copies.
        The reshape yields the matrix-per-q form that the downstream
        channel-solver (``_solve_rpa``) consumes as ``ham``.  Because the S/C
        matrices are q-independent constants for on-site Kanamori, they are
        cached across SCF iterations; chi0q itself is passed straight through
        and so needs no per-iteration work here.
        """
        logger.debug(">>> FLEX._inflate_chi0q_and_ham_general")

        # RPA's _calc_chi0q labels the rank-6 bare bubble as chi0[..., a, c, b, d]
        # (rpa.py: chi0[a,c,b,d] = G[a,b]·G_rev[d,c]), i.e. as a (norb^2, norb^2)
        # matrix with row=(a,c) and col=(b,d) it is
        #     chi0_{(mn),(μν)}(q) = -(T/N) Σ_k G_{mμ}(k+q) G_{νn}(k),
        # which is ALREADY the orbital-pair order that the downstream general-path
        # machinery (MYO S/C matrices, _solve_channels_general, _calc_veff_general,
        # _sigma_orbital_contract) consumes.  No orbital-pair transpose here.
        #
        # This used to apply a `.transpose(0, 1, 4, 5, 2, 3)` "RPA→MYO" conversion.
        # That transpose was wrong (issue #91): it is invisible on the orbital
        # DIAGONAL and it cancelled out of the public chi_s/chi_c (which were
        # transposed back on the way out -- exactly so for symmetric S/C, see
        # _flex_compute_veff_general), but it transposed V_eff and hence the
        # orbital off-diagonal of Sigma.  For a density-only interaction
        # (CoulombIntra) it made the general path disagree with the reduced path
        # by ~20% of the diagonal scale while both susceptibilities still matched
        # to machine precision.  The correct index order is pinned by an exact
        # ground truth independent of this codebase: the O(U^2) self-energy of a
        # 3-orbital exactly-diagonalized model (tests/test_flex_sopt_index_order.py).
        chi0q = chi0q_raw
        assert chi0q.ndim == 6

        # The MYO S/C matrices depend only on the interaction (q-resolved for
        # off-site entries, constant over q for on-site Kanamori), never on the
        # SCF state, so build them once and cache across iterations. The chi0
        # transpose above stays OUTSIDE the cache (chi0 changes every
        # iteration).
        cache = getattr(self, "_myo_sc_cache", None)
        if cache is None:
            from hwave.sc import _build_interaction_k
            from hwave.solver._sc_matrices_myo import build_sc_matrices_myo

            # PairLift does not contribute to the particle-hole spin/charge
            # (S/C) vertex: its contribution is S=C=0 (verified against the full
            # 4-index RPA; cf. the same treatment in hwave.sc), so the MYO S/C
            # builder correctly omits it. It is therefore physically inert here,
            # but warn so a configured PairLift term is not silently ignored --
            # matching the Eliashberg-path wording.
            if "PairLift" in self.ham_info.param_ham:
                logger.warning(
                    "PairLift is configured but does not contribute to the S/C "
                    "pairing vertex (S=C=0); it is ignored in the general FLEX "
                    "calculation.")

            # OFF-SITE entries are allowed ONLY where FLEX at one iteration
            # is MEASURED equal to the RPA ring, element-complete: CoulombInter
            # with SAME-orbital pairs (a == b), without sublattice folding.
            # For that class the vertex is V(q) on the density slots alone and
            # the equivalence holds to 1e-16 at one and two orbitals, on 4x4,
            # non-cubic 4x6 and 3D 4x4x2 lattices including a z-direction bond.
            #
            # Everything else stays rejected, each for a measured reason:
            #
            #   * CoulombInter with a != b off-site, Hund, Ising -- the MYO
            #     S/C builder applies the full on-site Kanamori slot mapping,
            #     which places q-dependent values into the Fierz (Case 2)
            #     inter-orbital slots; for R != 0 the particle-hole pair
            #     behind those slots is NON-LOCAL and not representable by a
            #     q-only matrix (the locality argument measured for the
            #     transverse channel). Off-site Hund / Ising differ from the
            #     ring by 3e-2 / 7e-2 even at ONE orbital, where no
            #     inter-orbital slot exists to blame -- an unadjudicated
            #     vertex-content disagreement, not a grid issue.
            #   * off-site combined with sublattice folding -- folding turns
            #     part of an a == b bond into intra-cell inter-orbital
            #     coupling, and the two solvers then differ by 2e-2 (the
            #     folded analogue of the #104 content); the equivalence claim
            #     no longer holds, so it is deferred to the #107 unification.
            #   * Exchange, PairHop -- non-local particle-hole pair off-site.
            #   * CoulombIntra -- `uhfk.py` reads only its r = 0 component
            #     (#106).
            def _normalized(tbl_dict):
                # The aggregate 'Coulomb' input carries CoulombIntra (the
                # r = 0 orbital-diagonal entries) and CoulombInter
                # (everything else) in one table; wan90.split_coulomb is the
                # shared decomposition. Without this normalization an
                # aggregate declaration silently produced a ZERO vertex in
                # this path -- neither the guard below nor
                # _build_interaction_k knows the 'Coulomb' key (measured:
                # chiq_s off by 1e-1 against the identical explicit
                # declaration).
                if "Coulomb" not in tbl_dict:
                    return tbl_dict
                if ("CoulombIntra" in tbl_dict
                        or "CoulombInter" in tbl_dict):
                    raise ValueError(
                        "Coulomb cannot be specified together with "
                        "CoulombIntra or CoulombInter")
                from hwave.qlmsio import wan90
                intra, inter = wan90.split_coulomb(tbl_dict["Coulomb"])
                # Rebuild through .copy(), which preserves the container's
                # CLASS. A dict comprehension here returned a plain dict
                # and so dropped the reader's CaseInsensitiveDict: the
                # reader stores each table under the spelling the USER
                # wrote (read_input_k.QLMSkInput), so every lookup below
                # -- the off-site guard AND _build_interaction_k -- then
                # missed a non-canonically-spelled type and silently
                # dropped that interaction. Measured before the fix: with
                # an aggregate Coulomb declaration, 'hund' produced chiq_s
                # identical to omitting Hund entirely, while 'Hund'
                # differed from it by 4.7e-2.
                out = tbl_dict.copy()
                del out["Coulomb"]
                out["CoulombIntra"] = intra
                out["CoulombInter"] = inter
                return out

            has_fold = tuple(getattr(self.lattice, "subshape",
                                     (1, 1, 1))) != (1, 1, 1)
            # Under folding the guard must scan the PRE-fold table:
            # _init_interaction canonicalizes displacements modulo the folded
            # grid, and a folded dimension of size one maps every off-site
            # displacement to (0,0,0) -- e.g. CellShape=[4,4,1] with
            # SubShape=[4,1,1] turns +-x bonds into zero-displacement
            # inter-sublattice entries, which the folded table cannot
            # distinguish from genuinely on-site input.
            scan_ham = self.ham_info.param_ham
            if has_fold:
                # fail CLOSED: validating off-site input against the folded
                # table is exactly the bypass this scan exists to prevent
                scan_ham = getattr(self.ham_info, "param_ham_orig", None)
                if scan_ham is None:
                    raise ValueError(
                        "sublattice folding is active but the pre-fold "
                        "interaction table (param_ham_orig) is missing, so "
                        "off-site input cannot be validated (the folded "
                        "table canonicalizes displacements and can hide "
                        "off-site entries).")
            scan_ham = _normalized(scan_ham)
            # PairLift is deliberately absent: hwave.solver.vertex_table
            # gives it NO particle-hole S/C content, so an off-site
            # PairLift declaration contributes exactly zero here and the
            # answer is right without a guard. Listing it would reject a
            # configuration that computes correctly.
            for itype in ("CoulombIntra", "CoulombInter", "Hund",
                          "Exchange", "PairHop", "Ising"):
                if itype not in scan_ham:
                    continue
                for (irvec, orbvec) in scan_ham[itype]:
                    if tuple(irvec) == (0, 0, 0):
                        continue
                    ok = (itype == "CoulombInter"
                          and orbvec[0] == orbvec[1]
                          and not has_fold)
                    if not ok:
                        raise ValueError(
                            "FLEX calc_scheme='general' accepts off-site "
                            "entries only for CoulombInter with equal "
                            "orbitals (a == b) and without sublattice "
                            "folding; interaction '{}' has an off-site "
                            "entry irvec={}, orbvec={}{}. For that entry "
                            "class the general path is measured equal to "
                            "the RPA ring; other off-site classes are not "
                            "representable by a q-only vertex or carry "
                            "unadjudicated vertex content.{}".format(
                                itype, tuple(irvec), tuple(orbvec),
                                ", with sublattice folding" if has_fold
                                else "",
                                _auto_general_remediation(self)))
                    # (Reader-bypassing internal tables only: file input
                    # rejects one-sided declarations since #93.)
                    # One-sided TABLES are fine: BOTH solvers reduce
                    # a declaration to its reversal-symmetric part (the
                    # physical coefficient of n_a(i) n_a(i+R), even in R by
                    # the site sum) -- the S/C builders since #114, and the
                    # ring's _make_ham_inter since the same reading was
                    # given to rpa.py. Measured: one-sided v and v/2 at
                    # both ends give bit-identical chiq in both solvers.

            no = self.norb
            nx, ny, nz = self.lattice.shape

            # Build k-space interactions from the raw real-space param_ham.
            # The grid contract: linspace(0, 2pi, n, endpoint=False) per axis,
            # C-order flattened -- the same points, order and flattening as
            # chi0's spatial FFT axis, verified by the element-complete
            # equivalence with the RPA ring for off-site (q-dependent)
            # entries. On-site terms are constant over q.
            kx = np.linspace(0, 2.0 * np.pi, nx, endpoint=False)
            ky = np.linspace(0, 2.0 * np.pi, ny, endpoint=False)
            kz = np.linspace(0, 2.0 * np.pi, nz, endpoint=False)
            inter_k = _build_interaction_k(
                kx, ky, kz, _normalized(self.ham_info.param_ham), no)

            # MYO S/C matrices: (nx, ny, nz, norb^2, norb^2).
            Us, Uc = build_sc_matrices_myo(inter_k, no, nx, ny, nz)

            # Reshape to (nvol, norb^2, norb^2) for the downstream channel solver.
            nvol = self.lattice.nvol
            Us = Us.reshape(nvol, no * no, no * no)
            Uc = Uc.reshape(nvol, no * no, no * no)

            self._myo_sc_cache = (Us, Uc)
        else:
            Us, Uc = cache

        # The S/C matrices are built (and cached) on the host; mirror them to
        # chi0q's backend so the downstream channel solve stays on one device.
        xp = _bk.array_module_of(chi0q)
        if xp is not np:
            Us = xp.asarray(Us)
            Uc = xp.asarray(Uc)

        return chi0q, Us, Uc

    @do_profile
    def _solve_channels_general(self, chi0q, Us, Uc):
        """Solve the spin and charge RPA channels for the general MYO path.

        In the MYO / Takimoto-Hotta-Ueda paramagnetic full-vertex convention
        (cond-mat/0407094 Eq. 4) the channel susceptibilities are

            chi_s = [I - chi0 . Us]^{-1} . chi0      (spin   channel)
            chi_c = [I + chi0 . Uc]^{-1} . chi0      (charge channel)

        i.e. the spin channel enters with a MINUS sign (it is the Stoner /
        magnetic-instability channel that diverges as chi0.Us -> I) and the
        charge channel with a PLUS sign.  This mirrors the sign discussion in
        :meth:`_build_spin_charge_vertices` for the reduced path.

        The inherited matrix solver :meth:`RPA._solve_rpa` computes
        ``[I + chi0 . ham]^{-1} . chi0``.  To realize the MYO signs we therefore
        pass ``ham_s = -Us`` (so ``+chi0.(-Us) = -chi0.Us``) and
        ``ham_c = +Uc``.  ``_solve_rpa`` accepts the orbital-space ``chi0q`` --
        a 6-dimensional array of shape ``(nmat, nvol, norb, norb, norb, norb)``
        (four orbital legs ``m,n,mu,nu``) -- and the ``(nvol, norb^2, norb^2)``
        interaction matrices directly (it reshapes chi0q to
        ``(nmat, nvol, norb^2, norb^2)`` internally), so ``Us``/``Uc`` from
        :meth:`_inflate_chi0q_and_ham_general` are used as-is.

        Parameters
        ----------
        chi0q : ndarray
            Bare susceptibility, shape ``(nmat, nvol, norb, norb, norb, norb)``
            (six dimensions: frequency, volume, and four orbital legs).
        Us : ndarray
            MYO spin (S) interaction matrices ``(nvol, norb^2, norb^2)``.
        Uc : ndarray
            MYO charge (C) interaction matrices ``(nvol, norb^2, norb^2)``.

        Returns
        -------
        chi_s : ndarray
            Spin-channel RPA susceptibility, same shape as ``chi0q``.
        chi_c : ndarray
            Charge-channel RPA susceptibility, same shape as ``chi0q``.
        """
        logger.debug(">>> FLEX._solve_channels_general")

        # _solve_rpa computes [I + chi0 ham]^-1 chi0; pass -Us / +Uc to obtain
        # chi_s = [I - chi0 Us]^-1 chi0 and chi_c = [I + chi0 Uc]^-1 chi0.
        chi_s = self._solve_rpa(chi0q, -Us)
        chi_c = self._solve_rpa(chi0q, +Uc)
        return chi_s, chi_c

    @do_profile
    def _build_spin_charge_vertices(self, ham_inflated):
        """Build spin and charge interaction vertices from inflated ham.

        The inflated ham is in reduced (nd x nd) form. We decompose it into
        spin and charge channels.

        For the Hubbard model with U on-site:
            ham_inflated has structure:
                [W_↑↑  W_↑↓]     where W_↑↑ = same-spin, W_↑↓ = cross-spin
                [W_↓↑  W_↓↓]

            Physical susceptibilities:
                chi_s = chi0 * [1 - U_s chi0]^{-1}  (spin, Stoner-enhanced)
                chi_c = chi0 * [1 + U_c chi0]^{-1}  (charge, suppressed)
            where:
                U_s = W_cross - W_same   (for Hubbard: U - 0 = U)
                U_c = W_cross + W_same   (for Hubbard: U + 0 = U)

            _solve_rpa convention: [1 + chi0 * ham]^{-1} chi0
                => ham_s = -U_s  (to get [1 - U_s chi0]^{-1})
                => ham_c = +U_c  (to get [1 + U_c chi0]^{-1})

        Parameters
        ----------
        ham_inflated : ndarray, shape (nvol, nd, nd)
            Interaction in reduced spin-orbital space.

        Returns
        -------
        ham_s : ndarray, shape (nvol, nd, nd)
            Spin channel vertex (with sign for _solve_rpa convention).
        ham_c : ndarray, shape (nvol, nd, nd)
            Charge channel vertex.
        """
        logger.debug(">>> FLEX._build_spin_charge_vertices")

        nvol = self.lattice.nvol
        norb = self.norb
        ns = self.ns
        nd = self.nd

        # ham_inflated is (nvol, nd, nd) where nd = norb * ns
        # Reshape to (nvol, ns, norb, ns, norb)
        ham_so = ham_inflated.reshape(nvol, ns, norb, ns, norb)

        # Same-spin block: (s, a) -> (s, a) same s
        # Cross-spin block: (s, a) -> (t, a) with s != t
        # Average over both spins for symmetry:
        w_same = 0.5 * (ham_so[:, 0, :, 0, :] + ham_so[:, 1, :, 1, :])
        w_cross = 0.5 * (ham_so[:, 1, :, 0, :] + ham_so[:, 0, :, 1, :])

        # U_s = w_cross - w_same (for Hubbard: U - 0 = U)
        # U_c = w_cross + w_same (for Hubbard: U + 0 = U)
        u_s = w_cross - w_same  # (nvol, norb, norb)
        u_c = w_cross + w_same  # (nvol, norb, norb)

        # For _solve_rpa [1 + chi0 * ham]^{-1} chi0:
        #   ham_s = -U_s  -> [1 - U_s chi0]^{-1} chi0  (Stoner enhancement)
        #   ham_c = +U_c  -> [1 + U_c chi0]^{-1} chi0  (charge suppression)
        # Inflate channel vertices to spin-orbital reduced space by scattering
        # onto the spin-block diagonal (Kronecker product with I_ns).
        # Bit-identical to np.einsum('kab,st->ksatb', u, I_ns), but faster.
        xp = _bk.array_module_of(ham_inflated)
        ham_s = xp.zeros((nvol, nd, nd), dtype=u_s.dtype)
        ham_c = xp.zeros((nvol, nd, nd), dtype=u_c.dtype)
        for s in range(ns):
            sl = slice(s * norb, (s + 1) * norb)
            ham_s[..., sl, sl] = -u_s
            ham_c[..., sl, sl] = u_c

        return ham_s, ham_c

    @do_profile
    def _calc_veff(self, chi0q, chi_s, chi_c, ham_inflated):
        """Compute effective FLEX interaction V_eff(q, ivn).

        V_eff = W * [3/2*chi_s + 1/2*chi_c - chi0] * W

        chi0 is subtracted exactly once: at zeroth order chi_s = chi_c = chi0,
        so the kernel reduces to chi0 and V_eff = W chi0 W reproduces the
        second-order (SOPT) bubble.  This single subtraction removes the
        double-counted second-order diagram contained in chi_s and chi_c.

        Parameters
        ----------
        chi0q : ndarray
            Bare susceptibility (inflated), shape (nfreq, nvol, nd, nd).
        chi_s : ndarray
            Spin susceptibility, same shape.
        chi_c : ndarray
            Charge susceptibility, same shape.
        ham_inflated : ndarray
            Interaction in reduced space, shape (nvol, nd, nd).

        Returns
        -------
        ndarray
            Effective interaction V_eff, shape (nfreq, nvol, nd, nd).
        """
        logger.debug(">>> FLEX._calc_veff")

        nfreq = chi0q.shape[0]
        nvol = chi0q.shape[1]
        nd = chi0q.shape[-1]

        # Reshape to matrices
        chi0q_2d = chi0q.reshape(nfreq, nvol, nd, nd)
        chi_s_2d = chi_s.reshape(nfreq, nvol, nd, nd)
        chi_c_2d = chi_c.reshape(nfreq, nvol, nd, nd)
        ham_2d = ham_inflated.reshape(nvol, nd, nd)

        # Fluctuation susceptibility.
        # FLEX kernel: 3/2 chi_s + 1/2 chi_c - chi0.  chi0 must be subtracted
        # exactly ONCE: at zeroth order chi_s = chi_c = chi0, so the kernel is
        # 3/2 chi0 + 1/2 chi0 - chi0 = chi0, giving V_eff = W chi0 W ~ U^2 (the
        # second-order bubble).  Subtracting chi0 in *both* terms would cancel
        # the O(U^0) part and lose the leading SOPT self-energy.
        fluct_chi = 1.5 * chi_s_2d + 0.5 * chi_c_2d - chi0q_2d

        # V_eff = W * fluct_chi * W
        # Use batched matmul instead of einsum for better BLAS utilization
        # W @ fluct_chi: broadcast (nvol, nd, nd) @ (nfreq, nvol, nd, nd)
        xp = _bk.array_module_of(chi0q)
        tmp = xp.matmul(ham_2d, fluct_chi)
        # tmp @ W: (nfreq, nvol, nd, nd) @ (nvol, nd, nd)
        v_eff = xp.matmul(tmp, ham_2d)

        return v_eff

    @do_profile
    def _calc_veff_general(self, chi0q, chi_s, chi_c, Us, Uc):
        r"""Compute the MYO full-vertex effective interaction V_eff(q, ivn).

        Implements the fluctuation part of the paramagnetic full-vertex
        interaction (Takimoto-Hotta-Ueda / MYO convention, cond-mat/0407094
        Eq. 4)::

            V(q) = 3/2 Us.chi_s.Us + 1/2 Uc.chi_c.Uc
                   - 1/4 (Us+Uc).chi0.(Us+Uc)

        The first-order constant terms (``+3/2 Us - 1/2 Uc``) are intentionally
        EXCLUDED so that, like the reduced path :meth:`_calc_veff`, the
        second-order (SOPT) limit is reproduced: at zeroth order
        ``chi_s = chi_c = chi0`` and (with ``Us = Uc = U``) the kernel reduces
        to ``U.chi0.U``, the leading second-order bubble.

        All three susceptibilities are six-dimensional orbital-space arrays of
        shape ``(nmat, nvol, norb, norb, norb, norb)`` (four orbital legs).
        They are flattened to ``(nmat, nvol, ndx, ndx)`` matrices with
        ``ndx = norb^2`` BEFORE the matrix products, using the same flatten
        convention as the MYO S/C builder: the row index is
        ``(l1 * norb + l2)`` (``idx12``) and the column index is
        ``(l3 * norb + l4)`` (``idx34``), matching ``Us`` / ``Uc`` / ``chi0q``.
        ``Us`` / ``Uc`` are ``(nvol, ndx, ndx)`` and are broadcast over the
        frequency axis.

        Parameters
        ----------
        chi0q : ndarray
            Bare susceptibility, shape ``(nmat, nvol, norb, norb, norb, norb)``.
        chi_s : ndarray
            Spin-channel RPA susceptibility, same shape as ``chi0q``.
        chi_c : ndarray
            Charge-channel RPA susceptibility, same shape as ``chi0q``.
        Us : ndarray
            MYO spin (S) interaction matrices ``(nvol, norb^2, norb^2)``.
        Uc : ndarray
            MYO charge (C) interaction matrices ``(nvol, norb^2, norb^2)``.

        Returns
        -------
        ndarray
            Effective interaction V_eff, shape
            ``(nmat, nvol, norb^2, norb^2)``.
        """
        logger.debug(">>> FLEX._calc_veff_general")

        nmat, nvol = chi0q.shape[0], chi0q.shape[1]
        no = chi0q.shape[2]
        ndx = no * no

        # Flatten the four orbital legs to (row, col) matrices BEFORE matmul.
        chi0_2d = chi0q.reshape(nmat, nvol, ndx, ndx)
        chis_2d = chi_s.reshape(nmat, nvol, ndx, ndx)
        chic_2d = chi_c.reshape(nmat, nvol, ndx, ndx)

        # Broadcast the (nvol, ndx, ndx) interaction matrices over frequency.
        UsB = Us[np.newaxis]            # (1, nvol, ndx, ndx)
        UcB = Uc[np.newaxis]
        Uspc = (Us + Uc)[np.newaxis]

        term_s = 1.5 * (UsB @ chis_2d @ UsB)
        term_c = 0.5 * (UcB @ chic_2d @ UcB)
        term_0 = 0.25 * (Uspc @ chi0_2d @ Uspc)

        v_eff = term_s + term_c - term_0    # (nmat, nvol, ndx, ndx)
        return v_eff

    @do_profile
    def _calc_self_energy(self, green_kw, v_eff, beta):
        """Compute self-energy Sigma(k, iwn) via FFT convolution.

        Sigma(k, iwn) = T/Nk * sum_{q,m} V_eff(q, ivm) * G(k-q, iwn-ivm)

        The convolution is computed efficiently using FFT:
        1. Transform G and V_eff to (r, tau) space
        2. Multiply in (r, tau) space: Sigma(r,tau) = V(r,tau) * G(r,tau)
        3. Transform back to (k, iwn) space

        ``green_kw`` is the FULL physical Green's function (NOT the
        tail-subtracted one used for chi0q): the self-energy convolution is an
        exact cyclic frequency convolution and needs no tau-space tail
        reconstruction.

        Parameters
        ----------
        green_kw : ndarray
            Dressed Green's function, shape (nblock, nmat, nvol, nd_block, nd_block).
        v_eff : ndarray
            Effective interaction, shape (nfreq, nvol, nd, nd).
        beta : float
            Inverse temperature.

        Returns
        -------
        ndarray
            Self-energy, shape (nblock, nmat, nvol, nd_block, nd_block).
        """
        logger.debug(">>> FLEX._calc_self_energy")

        nx, ny, nz = self.lattice.shape
        nvol = self.lattice.nvol
        nmat = self.nmat
        nblock = green_kw.shape[0]
        nd_block = green_kw.shape[-1]
        nd_v = v_eff.shape[-1]
        nfreq = v_eff.shape[0]

        # The self-energy convolution Sigma(r,tau) = V_eff(r,tau) * G(r,tau)
        # requires V_eff and G to share the same imaginary-time grid, i.e. the
        # bosonic frequency count must equal the fermionic one.  A smaller
        # nfreq would silently leave tau slices at zero (n_common below), so
        # fail loudly instead of returning a wrong self-energy.
        if nfreq != nmat:
            raise ValueError(
                "FLEX self-energy requires a full bosonic frequency grid: "
                "V_eff has nfreq={} but Nmat={}".format(nfreq, nmat))

        xp = _bk.array_module_of(green_kw)
        workers = getattr(self, "fft_workers", 1)

        # --- Transform Green's function to (r, tau) space ---

        # Matsubara freq -> imaginary time for G (fermionic)
        green_flat = green_kw.reshape(nblock, nmat, nvol * nd_block * nd_block)
        green_kt = _ms.fermion_to_tau(green_flat, axis=1).reshape(
            nblock, nmat, nx, ny, nz, nd_block * nd_block)

        # k-space -> real-space for G
        green_rt = _bk.spatial_ifftn(green_kt, axes=(2, 3, 4), workers=workers
                                     ).reshape(nblock, nmat, nvol, nd_block, nd_block)

        # --- Transform V_eff to (r, tau) space ---

        # Bosonic Matsubara freq -> imaginary time
        v_flat = v_eff.reshape(nfreq, nvol * nd_v * nd_v)
        v_qt = _ms.boson_to_tau(v_flat, axis=0).reshape(
            nfreq, nx, ny, nz, nd_v * nd_v)

        # q-space -> real-space for V_eff
        v_rt = _bk.spatial_ifftn(v_qt, axes=(1, 2, 3), workers=workers
                                 ).reshape(nfreq, nvol, nd_v, nd_v)

        # --- Compute Sigma(r, tau) ---
        n_common = min(nfreq, nmat)

        # If nd_block != nd_v, each Green block only sees its own spin slot
        # of V_eff (spin-free: nd_block = norb, nd_v = norb*ns = nd)
        if nd_block != nd_v:
            # The inflated V_eff carries the (raw) chi0q block g on its spin
            # diagonal slot (g, g) (see _inflate_chi0q_and_ham), so:
            #   spin-diag (nblock == ns): Green block g pairs with slot g
            #   spin-free (nblock == 1):  both slots are identical, use slot 0
            # Sigma stays in block space (nd_block x nd_block); the spin
            # off-diagonal slots of the inflated form are exactly zero.
            norb = self.norb
            sigma_rt = xp.zeros((nblock, nmat, nvol, nd_block, nd_block),
                                dtype=np.complex128)
            for g in range(nblock):
                s = g if nblock == self.ns else 0
                sl = slice(s * norb, (s + 1) * norb)
                sigma_rt[g, :n_common] = (
                    v_rt[:n_common, :, sl, sl]
                    * green_rt[g, :n_common]
                )
        else:
            # Sigma(r,tau) = V_eff(r,tau) * G(r,tau)  (element-wise product)
            sigma_rt = xp.zeros((nblock, nmat, nvol, nd_v, nd_v),
                                dtype=np.complex128)
            sigma_rt[:, :n_common] = v_rt[:n_common] * green_rt[:, :n_common]

        # --- Transform Sigma back to (k, iwn) space ---
        nd_sig = sigma_rt.shape[-1]

        # Real-space -> k-space
        sigma_kt = _bk.spatial_fftn(
            sigma_rt.reshape(nblock, nmat, nx, ny, nz, nd_sig * nd_sig),
            axes=(2, 3, 4), workers=workers
        ).reshape(nblock, nmat, nvol * nd_sig * nd_sig)

        # Imaginary time -> Matsubara freq (fermionic)
        sigma_kw = (_ms.tau_to_fermion(sigma_kt, axis=1)
                    .reshape(nblock, nmat, nvol, nd_sig, nd_sig) * (1.0 / beta))

        return sigma_kw

    @do_profile
    def _sigma_orbital_contract(self, v_rt, green_rt):
        r"""Per-(r,tau) rank-4 orbital contraction.

        Implements the orbital part of MYO (cond-mat/0407094) Eq. 3::

            Sigma_{mn} = sum_{mu,nu} V_{(mu m),(nu n)} G_{mu nu}

        Parameters
        ----------
        v_rt : ndarray, shape (nfreq, nvol, norb^2, norb^2)
            MYO vertex in (r, tau).  The last two axes flatten the orbital
            pairs (mu, m) [row] and (nu, n) [col], i.e.
            ``v_rt.reshape(..., norb, norb, norb, norb)`` has axes
            (freq, vol, mu, m, nu, n).
        green_rt : ndarray, shape (nblock, nfreq, nvol, norb, norb)
            Green's function in (r, tau), last two axes (mu, nu).

        Returns
        -------
        ndarray, shape (nblock, nfreq, nvol, norb, norb)
            Self-energy in (r, tau), last two axes (m, n).
        """
        xp = _bk.array_module_of(green_rt)
        nfreq, nvol = v_rt.shape[0], v_rt.shape[1]
        no = green_rt.shape[-1]
        v6 = v_rt.reshape(nfreq, nvol, no, no, no, no)   # (f, r, mu, m, nu, n)
        # Sigma[g,f,r,m,n] = sum_{mu,nu} v6[f,r,mu,m,nu,n] * green_rt[g,f,r,mu,nu]
        return xp.einsum('frumvn,gfruv->gfrmn', v6, green_rt)

    @do_profile
    def _calc_self_energy_general(self, green_kw, v_eff, beta):
        r"""Compute Sigma(k, iwn) for the general (full-vertex) FLEX path.

        Implements MYO (cond-mat/0407094) Eq. 3::

            Sigma_{mn}(k, iwn)
              = T/Nk * sum_{q,iv} sum_{mu,nu}
                       V_{(mu m),(nu n)}(q, iv) G_{mu nu}(k-q, iwn-iv)

        The frequency/momentum FFT transport (G -> (r,tau), V -> (r,tau),
        per-(r,tau) combine, then back to (k, iwn)) is IDENTICAL to the reduced
        :meth:`_calc_self_energy`; see that method for the derivation.  The ONLY
        difference is the per-(r,tau) combine step: the reduced path uses an
        element-wise Hadamard ``V * G``, whereas here it is the rank-4 orbital
        contraction :meth:`_sigma_orbital_contract`.

        Because the general path is pure-orbital (``nd_block == nd_v == norb``),
        the spin-inflation branch of the reduced method is not needed and is
        dropped.

        Parameters
        ----------
        green_kw : ndarray, shape (nblock, nmat, nvol, norb, norb)
            Dressed Green's function (orbital axes (mu, nu)).
        v_eff : ndarray, shape (nfreq, nvol, norb^2, norb^2)
            MYO effective interaction; last two axes flatten the orbital pairs
            (mu, m) [row] and (nu, n) [col].
        beta : float
            Inverse temperature.

        Returns
        -------
        ndarray, shape (nblock, nmat, nvol, norb, norb)
            Self-energy (orbital axes (m, n)).
        """
        logger.debug(">>> FLEX._calc_self_energy_general")

        nx, ny, nz = self.lattice.shape
        nvol = self.lattice.nvol
        nmat = self.nmat
        nblock = green_kw.shape[0]
        norb = green_kw.shape[-1]
        ndx = v_eff.shape[-1]   # = norb^2
        nfreq = v_eff.shape[0]

        # Shared with _calc_self_energy: V_eff and G must share the imaginary-
        # time grid (bosonic freq count == fermionic), else tau slices would be
        # silently left at zero.  Fail loudly.
        if nfreq != nmat:
            raise ValueError(
                "FLEX self-energy requires a full bosonic frequency grid: "
                "V_eff has nfreq={} but Nmat={}".format(nfreq, nmat))

        workers = getattr(self, "fft_workers", 1)

        # --- Transform Green's function to (r, tau) space ---
        # (transport identical to _calc_self_energy)
        green_flat = green_kw.reshape(nblock, nmat, nvol * norb * norb)
        green_kt = _ms.fermion_to_tau(green_flat, axis=1).reshape(
            nblock, nmat, nx, ny, nz, norb * norb)
        green_rt = _bk.spatial_ifftn(green_kt, axes=(2, 3, 4), workers=workers
                                     ).reshape(nblock, nmat, nvol, norb, norb)

        # --- Transform V_eff to (r, tau) space ---
        # (transport identical to _calc_self_energy)
        v_flat = v_eff.reshape(nfreq, nvol * ndx * ndx)
        v_qt = _ms.boson_to_tau(v_flat, axis=0).reshape(
            nfreq, nx, ny, nz, ndx * ndx)
        v_rt = _bk.spatial_ifftn(v_qt, axes=(1, 2, 3), workers=workers
                                 ).reshape(nfreq, nvol, ndx, ndx)

        # --- Compute Sigma(r, tau): rank-4 orbital contraction (the only new
        # bug-prone step; the rest is reused transport). ---
        sigma_rt = self._sigma_orbital_contract(v_rt, green_rt)

        # --- Transform Sigma back to (k, iwn) space ---
        # (transport identical to _calc_self_energy)
        sigma_kt = _bk.spatial_fftn(
            sigma_rt.reshape(nblock, nmat, nx, ny, nz, norb * norb),
            axes=(2, 3, 4), workers=workers
        ).reshape(nblock, nmat, nvol * norb * norb)

        sigma_kw = (_ms.tau_to_fermion(sigma_kt, axis=1)
                    .reshape(nblock, nmat, nvol, norb, norb) * (1.0 / beta))

        return sigma_kw

    def _calc_convergence(self, sigma_old, sigma_new):
        """Calculate convergence criterion for self-energy.

        Parameters
        ----------
        sigma_old : ndarray
            Previous self-energy.
        sigma_new : ndarray
            New self-energy.

        Returns
        -------
        float
            Relative difference |sigma_new - sigma_old| / |sigma_new|.
        """
        xp = _bk.array_module_of(sigma_new)
        diff = float(xp.linalg.norm(sigma_new - sigma_old))
        norm = float(xp.linalg.norm(sigma_new))
        if norm < 1.0e-30:
            return diff
        return diff / norm

    @do_profile
    def save_results(self, info_outputfile, green_info):
        """Save FLEX calculation results.

        Parameters
        ----------
        info_outputfile : dict
            Output file configuration.
        green_info : dict
            Calculation results.
        """
        logger.info("Save FLEX results")
        path_to_output = info_outputfile["path_to_output"]

        self._init_wavevec()

        # FLEX never applies the matsubara_frequency filter to its outputs
        # (unlike RPA.solve): every stored frequency axis is the full nmat
        # grid, so the metadata must describe that grid -- writing the
        # restricted user option (self.freq_index) would make the strict
        # hwave_sc loader reject a valid file.
        full_freq_index = np.arange(self.nmat)
        if len(self.freq_index) < self.nmat:
            logger.warning(
                "matsubara_frequency is restricted but FLEX outputs always "
                "hold the full frequency grid; the option is ignored.")

        # Stage 3: IR-native files replace the uniform-grid freq_index with
        # the sparse-node schema (design ir-matsubara-stage3.md Sec. 3).
        # freq_index is OMITTED on purpose -- its positional uniform-grid
        # semantics would lie about a node-resolved axis. (getattr: tests
        # drive save_results on __new__-built stubs without __init__.)
        ir_native = (getattr(self, "use_ir", False)
                     and not getattr(self, "write_densified", True))

        def _freq_meta(statistics):
            if not ir_native:
                return {"freq_index": full_freq_index}
            ax = self._ir_axF if statistics == "F" else self._ir_axB
            return {"matsubara_basis": "ir",
                    "frequency_grid": "sparse_ir_nodes",
                    "ir_freq_n": ax.freq_n,
                    "ir_beta": ax.beta,
                    "ir_wmax": ax.wmax,
                    "ir_tol": ax.eps,
                    "ir_L": ax.L,
                    "ir_statistics": ax.statistics}

        # Save physical observables (particle number, spin, chemical potential)
        # as a plain-text energy file, mirroring UHFr/UHFk.  In fixed-mu mode
        # the NCond line is the mu-N single point for the given mu.
        if "energy" in info_outputfile:
            physics = green_info.get("physics", getattr(self, "physics", None))
            if physics is None:
                logger.warning("save_results: no physics data to write to "
                               "'{}'".format(info_outputfile["energy"]))
            else:
                file_name = os.path.join(path_to_output,
                                         info_outputfile["energy"])
                with open(file_name, "w") as fw:
                    fw.write("NCond = {}\n".format(physics["NCond"]))
                    fw.write("Sz = {}\n".format(physics["Sz"]))
                    fw.write("ChemicalPotential = {}\n".format(physics["mu"]))
                logger.info("save_results: save energy in file {}".format(
                    file_name))

        # Save chi0q
        if "chi0q" in info_outputfile:
            file_name = os.path.join(path_to_output, info_outputfile["chi0q"])
            # coeff_tail provenance (issue #80): the tail correction changes
            # chi0q at O(1); record the producing value. On the IR path the
            # uniform-grid tail machinery is BYPASSED (aa = 0.0 if self.use_ir
            # in solve(); the fermionic basis carries the 1/(i w) tail
            # exactly), so the configured value was never applied -- omit the
            # key rather than claim a correction that did not happen. This
            # covers IR-native AND densified-IR output. (getattr: tests drive
            # save_results on __new__-built stubs without __init__.)
            from hwave.solver.rpa import TAIL_ENDPOINT_CONVENTION
            tail_meta = ({} if getattr(self, "use_ir", False)
                         else {"coeff_tail": getattr(self, "coeff_tail", 0.0),
                               # endpoint convention marker (issue #134)
                               "tail_endpoint": TAIL_ENDPOINT_CONVENTION})
            np.savez(file_name,
                     chi0q=green_info["chi0q"],
                     # full grid size: lets consumers locate the zero bosonic
                     # frequency (index nmat//2) unambiguously (run
                     # provenance only on IR-native files)
                     nmat=self.nmat,
                     momentum_convention=MOMENTUM_CONVENTION,
                     wavevector_unit=self.kvec,
                     wavevector_index=self.wavenum_table,
                     # FLEX chi0q comes from the same spin-block-ordered RPA
                     # internals; tag it so the chi0q consumers (RPA read_chi0q,
                     # hwave_sc _load_chi0q) accept the file in SO mode.
                     index_convention="spin_block",
                     **tail_meta,
                     **_freq_meta("B"),
                     **_scheme_stamp(self))
            logger.info("save_results: save chi0q in file {}".format(file_name))

        # Save susceptibilities (spin and charge channels separately for
        # Eliashberg). The chi_convention tag tells the downstream consumer
        # (_load_flex_susceptibilities / _compute_vertices_flex) how to
        # interpret the array layout: "myo" marks the general full-vertex
        # path's orbital-pair shape (nd = norb^2), "kuroki" the reduced path's
        # spin-orbital shape (nd = norb*ns). (The two builders' S/C VALUES
        # were unified by the #113 adjudication; the tag remains a layout /
        # provenance discriminator.)
        #
        # sc_vertex_version records which S/C vertex content the run used, so
        # a downstream Eliashberg run cannot silently pair these
        # susceptibilities with vertices from a different adjudication.
        # Version 2 = the #113 exact-diagonalization values (Hund / Exchange /
        # Ising corrected); files without the field predate #113.
        common_meta = dict(nmat=self.nmat,
                           wavevector_unit=self.kvec,
                           wavevector_index=self.wavenum_table,
                           chi_convention=("myo" if self._flex_general
                                           else "kuroki"),
                           sc_vertex_version=2,
                           # Fourier-sign provenance (issue #133): q labels
                           # follow the documented e^{+iqR} convention
                           momentum_convention=MOMENTUM_CONVENTION,
                           **_freq_meta("B"),
                           **_scheme_stamp(self))
        if self._flex_general:
            # Self-describing index-order marker for the ORBITAL-PAIR files
            # only: their axes are the pairs (a,c) and (b,d), stored in that
            # order. It lets the reader tell these files from pre-fix general
            # outputs, which stored the same arrays orbital-pair TRANSPOSED
            # under the same "myo" tag and are otherwise indistinguishable, and
            # makes any future layout change fail fast instead of being
            # silently misread.
            #
            # Deliberately NOT stamped on reduced files: those axes are
            # spin-orbital (s*norb + a), not four orbital legs -- the loader has
            # to extract a spin block before the array becomes an orbital-pair
            # object at all (see sc._to_orbital_pair). Calling them "acbd" would
            # be false machine-readable metadata. The reader accepts
            # marker-less "kuroki" files; only "myo" requires the marker.
            common_meta["chi_orbital_layout"] = "acbd"

        if "chiq_s" in green_info:
            file_name = os.path.join(path_to_output,
                                     info_outputfile.get("chiq_s", "chiq_s"))
            np.savez(file_name, chiq_s=green_info["chiq_s"], **common_meta)
            logger.info("save_results: save chiq_s in file {}".format(file_name))

        if "chiq_c" in green_info:
            file_name = os.path.join(path_to_output,
                                     info_outputfile.get("chiq_c", "chiq_c"))
            np.savez(file_name, chiq_c=green_info["chiq_c"], **common_meta)
            logger.info("save_results: save chiq_c in file {}".format(file_name))

        if "chiq" in info_outputfile:
            file_name = os.path.join(path_to_output, info_outputfile["chiq"])
            save_dict = dict(**common_meta)
            if "chiq_s" in green_info:
                save_dict["chiq_s"] = green_info["chiq_s"]
            if "chiq_c" in green_info:
                save_dict["chiq_c"] = green_info["chiq_c"]
            np.savez(file_name, **save_dict)
            logger.info("save_results: save chiq in file {}".format(file_name))

        # Save self-energy
        if "sigma" in info_outputfile:
            file_name = os.path.join(path_to_output, info_outputfile["sigma"])
            # #167: no scheme stamp here -- sigma/green are outside the
            # spec's listed npz-stamp scope (deliberate, not an omission).
            np.savez(file_name,
                     sigma=green_info.get("sigma"),
                     wavevector_unit=self.kvec,
                     wavevector_index=self.wavenum_table,
                     cell_shape=np.array(self.lattice.shape),
                     momentum_convention=MOMENTUM_CONVENTION,
                     **_freq_meta("F"))
            logger.info("save_results: save sigma in file {}".format(file_name))

        # Save Green's function
        if "green" in info_outputfile:
            file_name = os.path.join(path_to_output, info_outputfile["green"])
            # #167: no scheme stamp here -- sigma/green are outside the
            # spec's listed npz-stamp scope (deliberate, not an omission).
            green_extra = {}
            if ir_native:
                # native-only provenance key (the densified green.npz key
                # set stays unchanged -- design R-S3-2)
                green_extra["cell_shape"] = np.array(self.lattice.shape)
            np.savez(file_name,
                     green=green_info.get("green"),
                     wavevector_unit=self.kvec,
                     wavevector_index=self.wavenum_table,
                     # temperature provenance (issue #86): the Eliashberg
                     # consumer rebuilds the Matsubara grid from its own
                     # beta and validates it against this field, so a file
                     # from a different temperature fails fast instead of
                     # silently corrupting the pair bubble
                     beta=1.0 / self.T,
                     momentum_convention=MOMENTUM_CONVENTION,
                     **green_extra,
                     **_freq_meta("F"))
            logger.info("save_results: save green in file {}".format(file_name))
