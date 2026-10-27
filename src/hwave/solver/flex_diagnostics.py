"""Diagnostics of the FLEX self-consistency (GitHub issue #282).

Two read-only helpers for the FLEX SCF loop. Neither changes a number the
solver computes; they only add log lines and sharpen one error message.

* The **Stoner factor** of a channel is

      lambda = max_q  max Re eig[ -chi0(q, i nu = 0) . W(q) ]   (real eigenvalues),

  where ``W`` is the vertex handed to ``RPA._solve_rpa`` (which solves
  ``chi = [1 + chi0 W]^-1 chi0``). The RPA-type resummation of that channel
  is statically unstable when ``lambda >= 1``. For the reduced scheme
  ``W = -U_s`` (spin) and ``W = +U_c`` (charge), so the spin value is the
  familiar ``max_q U chi0(q, 0)`` of the one-band Hubbard model; for the
  general scheme ``W = -S`` / ``+C`` with the full orbital vertex; for the
  spin-orbital (spinful) scheme ``W`` is the antisymmetrized vertex of the
  single channel. ``chi0`` is the bubble of the CURRENT dressed Green
  function, so the value is the dressed Stoner factor of that iteration.

* :class:`ScfMonitor` remembers the last completed iteration (its number,
  the relative self-energy step and the Stoner factors), logs the Stoner
  line, and warns once when the iteration looks divergent. The chemical-
  potential search reads it to explain a failure (:func:`mu_failure_message`).
"""

import math

import numpy as np

from . import backend as _bk


def static_boson_index(solver):
    """Index of the bosonic i nu = 0 slice on the solver's working axis.

    Uniform grid: ``nmat // 2`` (the centered even grid). IR: the node with
    ``freq_n == 0`` (present for every basis sparse-ir builds in practice);
    should it be absent, the node closest to zero is used and the second
    return value is False so the log can say so.
    """
    axB = getattr(solver, "_ir_axB", None)
    if getattr(solver, "use_ir", False) and axB is not None:
        freq_n = np.asarray(axB.freq_n)
        idx = int(np.argmin(np.abs(freq_n)))
        return idx, bool(freq_n[idx] == 0)
    return int(solver.nmat) // 2, True


#: An eigenvalue of -chi0 W counts as real when |Im| <= REAL_TOL * (the
#: largest |eigenvalue| at that q). Only a REAL eigenvalue reaching 1 makes
#: 1 + chi0 W singular; a complex pair with Re >= 1 does not. The tolerance
#: is loose on purpose: chi0 W is similar to a Hermitian matrix in the
#: paramagnetic cases, but degenerate eigenvalues (spin multiplets) of a
#: badly conditioned product may split into pairs with a tiny imaginary part.
REAL_TOL = 1.0e-6


def static_slice(chi0q, static_index):
    """The i nu = 0 slice of a bubble as a host ``(nvol, ndx, ndx)`` array
    (one device-to-host copy, shared by every channel of an iteration)."""
    nvol = chi0q.shape[1]
    trail = chi0q.shape[2:]
    ndx = int(np.prod(trail[:len(trail) // 2]))
    return np.asarray(_bk.to_host(chi0q[static_index])).reshape(nvol, ndx, ndx)


def static_instability(c, vertex, lattice_shape):
    """Stoner factor of one channel from the static bubble slice ``c``.

    ``lambda = max_q max { Re e : e eigenvalue of -c(q) W(q), e real }``,
    where ``W`` is the vertex handed to ``RPA._solve_rpa`` (it solves
    ``[1 + chi0 W]^-1 chi0``) and "real" is ``|Im e| <= REAL_TOL max|e|``.

    Parameters
    ----------
    c : ndarray ``(nvol, ndx, ndx)``
        Host static slice from :func:`static_slice`.
    vertex : array
        Reshapeable to ``(nvol, ndx, ndx)``.
    lattice_shape : tuple of int
        ``(nx, ny, nz)``; the flat q index is C-ordered over it.

    Returns
    -------
    (float, tuple, float)
        The Stoner factor, the q vector at which it is attained (units of
        the reciprocal lattice vectors) and ``max_q max Re e`` over ALL
        eigenvalues (equal to the first value unless a complex pair lies
        further right). ``(nan, None, nan)`` when nothing is finite.
    """
    nvol, ndx = c.shape[0], c.shape[-1]
    w = np.asarray(_bk.to_host(vertex)).reshape(nvol, ndx, ndx)
    if not (np.all(np.isfinite(c)) and np.all(np.isfinite(w))):
        return float("nan"), None, float("nan")
    ev = np.linalg.eigvals(-(c @ w))                       # (nvol, ndx)
    scale = np.abs(ev).max(axis=-1, keepdims=True)
    real = np.abs(ev.imag) <= REAL_TOL * np.maximum(scale, 1e-300)
    lam = np.where(real, ev.real, -np.inf).max(axis=-1)
    max_re = float(ev.real.max())
    if not np.any(np.isfinite(lam)):
        return float("nan"), None, max_re
    iq = int(np.argmax(lam))
    idx = np.unravel_index(iq, tuple(int(n) for n in lattice_shape))
    q = tuple(float(i) / float(n) for i, n in zip(idx, lattice_shape))
    return float(lam[iq]), q, max_re


def format_stoner(stoner):
    """One log fragment for a ``{channel: (value, q[, max_re])}`` mapping.
    ``max_re`` (complex eigenvalue pairs further right than every real one)
    is printed only when it differs from the value."""
    parts = []
    for name, entry in stoner.items():
        val, q = entry[0], entry[1]
        if q is None:
            txt = "{} {:.6g}".format(name, val)
        else:
            txt = "{} {:.6f} at q = ({})".format(
                name, val, ", ".join("{:.4g}".format(x) for x in q))
        if len(entry) > 2 and math.isfinite(entry[2]) and not (
                math.isfinite(val) and entry[2] <= val + 1e-12 * max(1.0, abs(val))):
            txt += " (max Re over complex pairs {:.6f})".format(entry[2])
        parts.append(txt)
    return "; ".join(parts)


class ScfMonitor:
    """Per-solve record of the FLEX iteration, for logs and error messages.

    The divergence WARNING fires at most once per solve and never stops the
    run. Its criterion is deliberately narrow: some channel's Stoner factor
    is >= 1 (the bubble of the current Green function is statically
    unstable, so the RPA-type sums of that channel are not physical) AND the
    relative self-energy step has been >= ``step_floor`` for the last
    ``window`` iterations in a row; or the step is not finite. A cold start
    may begin with a Stoner factor slightly above 1 and recover (the
    one-orbital tutorial sample starts at 1.056 for one iteration), which is
    why one unstable iteration is not enough. The criterion was checked not
    to fire on any converging run of the FLEX test suite or the tutorial
    samples; it is a heuristic, which is why it only warns.
    """

    window = 3
    step_floor = 0.5

    def __init__(self, logger):
        self.logger = logger
        self.iteration = None
        self.step = None
        self.stoner = None
        self.warned = False
        self._run = 0

    def record(self, iteration, step, stoner, last):
        """Record a completed iteration and log its Stoner line.

        ``stoner`` is a ``{channel: (value, q)}`` mapping or None (no bubble
        of that form, e.g. the bond-resolved path). ``last`` marks the final
        iteration of the loop (converged or IterationMax reached).
        """
        self.iteration = int(iteration)
        self.step = float(step)
        self.stoner = stoner
        unstable = bool(stoner) and any(
            e[0] >= 1.0 for e in stoner.values() if math.isfinite(e[0]))
        if unstable and self.step >= self.step_floor:
            self._run += 1
        else:
            self._run = 0
        diverging = (not math.isfinite(self.step)) or self._run >= self.window
        if stoner:
            line = "  Stoner factor max_q eig[U chi0(q, i nu = 0)]: " + format_stoner(stoner)
            if iteration == 1 or last or (diverging and not self.warned):
                self.logger.info(line)
            else:
                self.logger.debug(line)
        if diverging and not self.warned:
            self.warned = True
            self.logger.warning(
                "FLEX iteration %d looks divergent: |dSigma|/|Sigma| = %.3e%s. %s",
                self.iteration, self.step,
                "" if not stoner else (
                    " with a Stoner factor >= 1 for {} iterations in a row".format(self._run)
                    if self._run >= self.window else ""),
                ADVICE)

    def context(self):
        """The last recorded iteration, as a sentence (empty before the first)."""
        if self.iteration is None:
            return ""
        s = " The last completed iteration was {} with |dSigma|/|Sigma| = {:.3e}".format(
            self.iteration, self.step)
        if self.stoner:
            s += " and Stoner factor " + format_stoner(self.stoner)
        return s + "."


ADVICE = ("Likely causes: too large a step in a warm-started scan (U or T; "
          "the previous self-energy is far from the new fixed point), too "
          "large a mixing (Mix), or the starting point. Use a smaller U/T "
          "step, a smaller Mix, or a cold start (no sigma_init).")


#: N(mu) counts as flat over the band window when it varies there by less
#: than this fraction of the number of states (see mu_failure_message).
FLAT_FRACTION = 1.0e-2


def band_window(ew, beta):
    """The physical mu window used to classify a failed root search: the
    bare (or Hartree-Fock) band range widened by 1 + 10 T on each side.
    It does NOT depend on the self-energy, so it stays at physical energies
    when Sigma has blown up (the expanded search bracket does not)."""
    pad = 1.0 + 10.0 / float(beta)
    return float(ew.min()) - pad, float(ew.max()) + pad


def mu_failure_message(where, lo, hi, f_lo, f_hi, ncond, nstate,
                       sigma_finite, window, f_window, monitor=None):
    """The error text when the chemical-potential root is not bracketed.

    ``ncond`` is the target and ``nstate`` the number of states it is
    counted over (both per spin for a spin-free run). A HEALTHY N(mu) is
    continuous and runs from 0 (mu -> -inf) to ``nstate`` (mu -> +inf), so it
    crosses every target strictly inside ``(0, nstate)``. Hence:

    * the self-energy is not finite -> diverged;
    * the target is not strictly inside ``(0, nstate)`` -> no filling can
      reach it: the filling/Ncond (and temperature) message;
    * otherwise the target was reachable and the search still failed, so
      N(mu) itself is broken -> diverged.

    The values at the ends of the EXPANDED bracket (``|mu|`` up to ~1e30)
    are not used to decide: on the IR path N(mu) is the basis fit of
    ``sum 1/(lambda + mu)`` and tends to 0 at BOTH ends whatever the
    self-energy is. To say HOW N(mu) is broken, ``f_window = N(mu) - Ncond``
    is read on a few points of the physical ``window`` (:func:`band_window`),
    where both evaluations of N(mu) are valid: not finite, flat (changes by
    less than ``FLAT_FRACTION * nstate``), outside the physical range
    ``[0, nstate]``, or crossing the target non-monotonically.
    """
    head = ("{}: chemical-potential root not bracketed after expansion to "
            "[{}, {}] (N-Ncond = {}, {}).".format(where, lo, hi, f_lo, f_hi))
    ctx = monitor.context() if monitor is not None else ""
    ncond, nstate = float(ncond), float(nstate)
    fw = [float(x) for x in f_window]
    wtxt = "[{:.6g}, {:.6g}]".format(window[0], window[1])
    if sigma_finite and not (0.0 < ncond < nstate):
        return head + (" The target Ncond = {:.6g} is not inside (0, {:.6g}), "
                       "the range of particle numbers these states can hold "
                       "(per spin for a spin-free run). Check the target "
                       "filling/Ncond and temperature.".format(ncond, nstate)) + ctx
    if not sigma_finite:
        why = "the self-energy is not finite"
    elif not all(math.isfinite(x) for x in fw):
        why = "N(mu) is not finite over the band window {}".format(wtxt)
    else:
        n = [x + ncond for x in fw]
        tol = FLAT_FRACTION * max(1.0, nstate)
        if max(n) - min(n) < tol:
            why = ("N(mu) does not depend on mu over the band window {}: it "
                   "changes by {:.3g} of {:.6g} states".format(
                       wtxt, max(n) - min(n), nstate))
        elif min(n) < -tol or max(n) > nstate + tol:
            why = ("N(mu) leaves the physical range [0, {:.6g}] over the band "
                   "window {} (from {:.6g} to {:.6g})".format(
                       nstate, wtxt, min(n), max(n)))
        elif min(n) <= ncond <= max(n):
            why = ("N(mu) over the band window {} runs from {:.6g} to {:.6g} "
                   "but does not cross Ncond = {:.6g} monotonically".format(
                       wtxt, n[0], n[-1], ncond))
        else:
            # a self-energy that moved spectral weight far outside the band
            # window, or (IR path) a basis bandwidth ir_wmax that does not
            # cover the dressed spectrum
            why = ("N(mu) stays within [{:.6g}, {:.6g}] over the band window "
                   "{} and never reaches Ncond = {:.6g}; on the IR path also "
                   "check that ir_wmax covers the spectrum".format(
                       min(n), max(n), wtxt, ncond))
    return head + (" The self-energy has diverged ({}); the target filling is "
                   "reachable, so this is not a problem of the filling or the "
                   "temperature.".format(why)) + ctx + " " + ADVICE
