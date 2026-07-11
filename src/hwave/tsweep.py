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


def eliashberg_frequency(base):
    return str(base.get("eliashberg", {}).get("frequency", "static"))


def resolve_sigma_name(base):
    name = base.get("file", {}).get("output", {}).get("sigma", "sigma")
    return name if name.endswith(".npz") else name + ".npz"


def resolve_gap_name(base):
    return "gap_dynamic.npz" if eliashberg_frequency(base) == "dynamic" else "gap.npz"


def rung_dir(output_dir, idx, T):
    return os.path.join(output_dir, "%03d_T%g" % (idx, T))


def preflight(base, cont):
    param = base.get("mode", {}).get("param", {})
    for field in ("CellShape", "Nmat"):
        if field not in param:
            raise ValueError("[mode.param] missing required field %r." % field)
    has_filling, has_ncond = "filling" in param, "Ncond" in param
    if has_filling == has_ncond:
        raise ValueError(
            "[mode.param] must set exactly one of `filling`/`Ncond`.")
    if cont.get("run_eliashberg", True) and "eliashberg" not in base:
        raise ValueError(
            "run_eliashberg=true but the base config has no [eliashberg] "
            "section; add [eliashberg] or set "
            "[continuation] run_eliashberg = false.")


def make_rung_dicts(base, T, rung_out, run_eliashberg,
                    sigma_init=None, seed_gap=None):
    canon = copy.deepcopy(base)
    canon.pop("continuation", None)
    canon.setdefault("file", {}).setdefault("input", {})
    canon["file"].setdefault("output", {})
    canon["file"]["input"].pop("sigma_init", None)
    if "eliashberg" in canon:
        canon["eliashberg"].pop("seed_eigenvector", None)
    canon["mode"]["param"]["T"] = T
    canon["file"]["output"]["path_to_output"] = rung_out
    if run_eliashberg:
        canon["file"]["input"]["path_to_flex_output"] = rung_out

    flex = copy.deepcopy(canon)
    if sigma_init is not None:
        flex["file"]["input"]["sigma_init"] = sigma_init

    eli = None
    if run_eliashberg:
        eli = copy.deepcopy(canon)
        if seed_gap is not None:
            eli["eliashberg"]["seed_eigenvector"] = seed_gap
    return flex, eli


def parse_leading_eig(path):
    if not os.path.exists(path):
        raise ValueError("eigenvalue file not found: %s" % path)
    for line in open(path):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        cols = s.split()
        try:
            idx = int(cols[0])
        except (ValueError, IndexError):
            continue
        if idx != 0:
            continue
        re, im = float(cols[1]), float(cols[2])
        match = int(cols[4]) if len(cols) >= 5 else -1
        return re, im, match
    raise ValueError("no leading (index 0) eigenpair row in %s" % path)


_SUMMARY_HEADER = ("# idx  T  status  error_stage  Re_lambda  Im_lambda  "
                   "parity_match  flex_converged  flex_iter\n")


def write_summary(path, rows):
    with open(path, "w") as fw:
        fw.write(_SUMMARY_HEADER)
        for r in rows:
            fw.write("%d %.12g %s %s %.6f %.6f %d %d %d\n" % (
                r["idx"], r["T"], r["status"], r["error_stage"],
                r["re"], r["im"], r["match"],
                r["flex_converged"], r["flex_iter"]))
