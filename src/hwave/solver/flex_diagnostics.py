"""Diagnostics of the FLEX self-consistency (GitHub issue #282).

Two read-only helpers for the FLEX SCF loop. Neither changes a number the
solver computes; they only add log lines and sharpen one error message.

* The **Stoner factor** of a channel is

      lambda = max_q  max Re eig[ -chi0(q, i nu = 0) . W(q) ],

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


def static_instability(chi0q, vertex, static_index, lattice_shape):
    """Stoner factor ``max_q max Re eig[-chi0(q, 0) W(q)]`` of one channel.

    Parameters
    ----------
    chi0q : array, shape ``(nfreq, nvol, ...)``
        Bubble in the layout ``_solve_rpa`` takes (the trailing axes split
        into two equal halves, the matrix row and column).
    vertex : array
        The vertex ``W`` passed to ``_solve_rpa``, reshapeable to
        ``(nvol, ndx, ndx)``.
    static_index : int
        Position of i nu = 0 on the frequency axis of ``chi0q``.
    lattice_shape : tuple of int
        ``(nx, ny, nz)``; the flat q index is C-ordered over it.

    Returns
    -------
    (float, tuple)
        The Stoner factor and the q vector at which it is attained, in
        units of the reciprocal lattice vectors. ``(nan, None)`` when the
        bubble or the vertex is not finite.
    """
    nvol = chi0q.shape[1]
    trail = chi0q.shape[2:]
    ndx = int(np.prod(trail[:len(trail) // 2]))
    c = np.asarray(_bk.to_host(chi0q[static_index])).reshape(nvol, ndx, ndx)
    w = np.asarray(_bk.to_host(vertex)).reshape(nvol, ndx, ndx)
    if not (np.all(np.isfinite(c)) and np.all(np.isfinite(w))):
        return float("nan"), None
    lam = np.linalg.eigvals(-(c @ w)).real.max(axis=-1)
    iq = int(np.argmax(lam))
    idx = np.unravel_index(iq, tuple(int(n) for n in lattice_shape))
    q = tuple(float(i) / float(n) for i, n in zip(idx, lattice_shape))
    return float(lam[iq]), q


def format_stoner(stoner):
    """One log fragment for a ``{channel: (value, q)}`` mapping."""
    parts = []
    for name, (val, q) in stoner.items():
        if q is None:
            parts.append("{} {:.6g}".format(name, val))
        else:
            parts.append("{} {:.6f} at q = ({})".format(
                name, val, ", ".join("{:.4g}".format(x) for x in q)))
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
            v >= 1.0 for v, _ in stoner.values() if math.isfinite(v))
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


def mu_failure_message(where, lo, hi, f_lo, f_hi, nstate, sigma_finite,
                       monitor=None):
    """The error text when the chemical-potential root is not bracketed.

    ``N(mu) - Ncond`` is evaluated at both ends of the fully expanded
    bracket. ``N`` runs from 0 to ``nstate`` as mu goes from -inf to +inf,
    so a root that cannot be bracketed has two causes that the two values
    tell apart:

    * N(mu) does not vary over the bracket (``|f_hi - f_lo|`` below
      ``1e-6 * nstate``), or something is not finite: the self-energy has
      diverged and the particle number no longer depends on mu;
    * N(mu) varies but never crosses the target: the target filling (or the
      temperature) is not reachable.
    """
    varies = (math.isfinite(f_lo) and math.isfinite(f_hi)
              and abs(f_hi - f_lo) > 1.0e-6 * max(1.0, float(nstate)))
    head = ("{}: chemical-potential root not bracketed after expansion to "
            "[{}, {}] (N-Ncond = {}, {}).".format(where, lo, hi, f_lo, f_hi))
    ctx = monitor.context() if monitor is not None else ""
    if sigma_finite and varies:
        return head + (" N(mu) varies over the bracket but never reaches the "
                       "target. Check the target filling/Ncond and "
                       "temperature.") + ctx
    if not sigma_finite:
        why = "the self-energy is not finite"
    elif not (math.isfinite(f_lo) and math.isfinite(f_hi)):
        why = "N(mu) is not finite"
    else:
        why = "N(mu) does not depend on mu over the whole bracket"
    return head + (" The self-energy has diverged ({}); this is not a problem "
                   "of the filling or the temperature.".format(why)) + ctx \
        + " " + ADVICE
