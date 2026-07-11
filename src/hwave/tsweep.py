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
    if cont.get("run_eliashberg", True):
        if "eliashberg" not in base:
            raise ValueError(
                "run_eliashberg=true but the base config has no [eliashberg] "
                "section; add [eliashberg] or set "
                "[continuation] run_eliashberg = false.")
        # The per-rung summary reports the parity-selected leading eigenpair
        # (Re/Im/parity), which sc.py writes to eigenvalue.dat only in the
        # "# Eigenvalue analysis" block -- present only for solver_mode
        # "eigenvalue"/"both".  The default "iteration" writes a bare leading
        # value with no parity/Im and no index-0 row, which parse_leading_eig
        # cannot read.  Require an eigenvalue-analysis mode up front rather than
        # failing on rung 0.
        smode = str(base["eliashberg"].get("solver_mode", "iteration"))
        if smode not in ("eigenvalue", "both"):
            raise ValueError(
                "[eliashberg] solver_mode = %r cannot drive a temperature "
                "sweep: hwave_tsweep reads the parity-selected leading "
                "eigenpair from the eigenvalue-analysis block, which is only "
                "written for solver_mode 'eigenvalue' or 'both'. Set "
                "[eliashberg] solver_mode = \"eigenvalue\"." % smode)


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


def _abspath(base_dir, p):
    return p if os.path.isabs(p) else os.path.normpath(os.path.join(base_dir, p))


def run(input_dict, base_dir=".", keep_going=False, dry_run=False):
    from . import qlms, sc                       # lazy: avoids heavy import at module load
    # Work on our own copy: run() resolves input paths to absolute below, and
    # the importable API must not mutate the caller's dict (mirrors the
    # non-mutation guarantee of make_rung_dicts).
    input_dict = copy.deepcopy(input_dict)
    base_dir = os.path.abspath(base_dir)
    cont = input_dict.get("continuation", {})
    preflight(input_dict, cont)
    ladder = build_ladder(cont)
    run_eli = cont.get("run_eliashberg", True)
    warm_start = cont.get("warm_start", True)
    seed_gap_on = cont.get("seed_gap", True) and eliashberg_frequency(input_dict) == "dynamic"
    out_dir = _abspath(base_dir, cont.get("output_dir", "tsweep"))
    summary_file = cont.get("summary_file", "lambda_vs_T.dat")

    # FLEX writes the self-energy only when [file.output] sigma is set (flex.py).
    # warm_start chains that file into the next rung, so ensure it is written;
    # inject the default name when the user omitted it rather than silently
    # cold-starting every rung.  tsweep already owns the per-rung output layout.
    fout = input_dict.setdefault("file", {}).setdefault("output", {})
    if warm_start and not fout.get("sigma"):
        fout["sigma"] = "sigma"
        logger.info("warm_start: [file.output] sigma not set; using 'sigma' so "
                    "the self-energy is written and chained between rungs.")
    sigma_name = resolve_sigma_name(input_dict)
    gap_name = resolve_gap_name(input_dict)

    # resolve base input paths to absolute so CWD does not matter
    fin = input_dict.setdefault("file", {}).setdefault("input", {})
    if "path_to_input" in fin:
        fin["path_to_input"] = _abspath(base_dir, fin["path_to_input"])
    inter = fin.get("interaction", {})
    if "path_to_input" in inter:
        inter["path_to_input"] = _abspath(base_dir, inter["path_to_input"])

    rows, prev = [], None                          # prev = seedable prev rung dir or None
    for idx, T in enumerate(ladder):
        rdir = rung_dir(out_dir, idx, T)
        rout = os.path.join(rdir, "output")
        sigma_init = os.path.join(prev, sigma_name) if (prev and warm_start) else None
        seed = os.path.join(prev, gap_name) if (prev and seed_gap_on) else None
        if dry_run:
            logger.info("[dry-run] rung %d T=%g dir=%s sigma_init=%s seed=%s",
                        idx, T, rdir, sigma_init, seed)
            rows.append(dict(idx=idx, T=T, status="dry", error_stage="none",
                             re=float("nan"), im=float("nan"), match=-1,
                             flex_converged=-1, flex_iter=-1))
            prev = rout
            continue
        flex_dict, eli_dict = make_rung_dicts(
            input_dict, T, rout, run_eli, sigma_init=sigma_init, seed_gap=seed)
        eig_name = eli_dict["eliashberg"].get("output_eigenvalue", "eigenvalue.dat") \
            if eli_dict is not None else "eigenvalue.dat"
        row = dict(idx=idx, T=T, status="ok", error_stage="none",
                   re=float("nan"), im=float("nan"), match=-1,
                   flex_converged=-1, flex_iter=-1)
        try:
            res = qlms.run(input_dict=flex_dict) or {}
            row["flex_converged"] = 1 if res.get("scf_converged") else 0
            row["flex_iter"] = int(res.get("scf_iterations", -1))
            if not res.get("scf_converged", False):
                row["status"] = "not_converged"
            if run_eli:
                sc.calc_eliashberg(eli_dict)
                re, im, match = parse_leading_eig(
                    os.path.join(rout, eig_name))
                row["re"], row["im"], row["match"] = re, im, match
        except Exception as exc:                   # noqa: BLE001 - record & stop/continue
            row["status"] = "error"
            row["error_stage"] = "eliashberg" if row["flex_iter"] >= 0 else "flex"
            logger.error("rung %d T=%g failed at %s: %s",
                         idx, T, row["error_stage"], exc)
        rows.append(row)
        # seedable iff not error and required files exist
        ok_seed = row["status"] != "error" and _seed_files_present(
            rout, warm_start, seed_gap_on, sigma_name, gap_name)
        if row["status"] == "error" and not keep_going:
            break
        prev = rout if ok_seed else None
        if row["status"] == "error":               # keep_going: cold-start next
            logger.warning("continuing after error at T=%g; next rung cold-started", T)

    os.makedirs(out_dir, exist_ok=True)
    write_summary(os.path.join(out_dir, summary_file), rows)
    return rows


def _seed_files_present(rout, warm_start, seed_gap_on, sigma_name, gap_name):
    if warm_start and not os.path.exists(os.path.join(rout, sigma_name)):
        logger.warning("expected sigma %s missing in %s; next rung cold-starts",
                       sigma_name, rout)
        return False
    if seed_gap_on and not os.path.exists(os.path.join(rout, gap_name)):
        logger.warning("expected gap %s missing in %s; next rung cold-starts",
                       gap_name, rout)
        return False
    return True


def main():
    import argparse
    import tomli
    parser = argparse.ArgumentParser(
        prog="hwave_tsweep",
        description="Temperature-continuation driver for FLEX(+Eliashberg).")
    parser.add_argument("input", help="base TOML with a [continuation] section")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the resolved ladder and seed wiring; run nothing")
    parser.add_argument("--keep-going", action="store_true",
                        help="on a rung error, cold-start the next rung instead of stopping")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    with open(args.input, "rb") as f:
        input_dict = tomli.load(f)
    run(input_dict, base_dir=os.path.dirname(os.path.abspath(args.input)),
        keep_going=args.keep_going, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
