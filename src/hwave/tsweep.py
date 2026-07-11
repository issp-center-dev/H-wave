"""hwave_tsweep -- temperature-continuation post-tool.

Runs FLEX(+Eliashberg) across a descending temperature ladder, chaining each
rung's converged self-energy (sigma_init) and dynamic gap (seed_eigenvector)
into the next.  See docs/.../2026-07-11-hwave-tsweep-continuation.md.
"""
import copy
import logging
import math
import os

logger = logging.getLogger("qlms").getChild("tsweep")


def build_ladder(cont):
    """Resolve the temperature ladder from a [continuation] dict."""
    temps = cont.get("temperatures")
    if temps is not None:
        ladder = [float(t) for t in temps]
    else:
        if not all(k in cont for k in ("T_start", "T_stop", "num")):
            raise ValueError(
                "[continuation] needs either `temperatures` or "
                "`T_start`/`T_stop`/`num` (+ optional `spacing`).")
        t0, t1, n = float(cont["T_start"]), float(cont["T_stop"]), int(cont["num"])
        spacing = cont.get("spacing", "linear")
        if n < 1:
            raise ValueError("[continuation] num must be >= 1.")
        if n == 1:
            ladder = [t0]
        elif spacing == "log":
            lg0, lg1 = math.log(t0), math.log(t1)
            ladder = [math.exp(lg0 + (lg1 - lg0) * i / (n - 1)) for i in range(n)]
        elif spacing == "linear":
            ladder = [t0 + (t1 - t0) * i / (n - 1) for i in range(n)]
        else:
            raise ValueError("spacing must be 'linear' or 'log', got %r" % spacing)
    seed_on = cont.get("warm_start", True) or cont.get("seed_gap", True)
    if seed_on and any(ladder[i] <= ladder[i + 1] for i in range(len(ladder) - 1)):
        logger.warning("temperature ladder does not strictly descend; "
                       "warm-start continuation is most effective on a "
                       "descending ladder.")
    return ladder
