"""hwave_tsweep -- temperature-continuation post-tool.

Runs FLEX(+Eliashberg) across a descending temperature ladder, chaining each
rung's converged self-energy (sigma_init) and dynamic gap (seed_eigenvector)
into the next.  See docs/.../2026-07-11-hwave-tsweep-continuation.md.
"""
import copy
import hashlib
import json
import logging
import math
import os

import numpy as np

logger = logging.getLogger("qlms").getChild("tsweep")

MANIFEST_NAME = "tsweep_manifest.json"
_MANIFEST_VERSION = 1


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


def _atomic_write_text(path, text):
    """Write ``text`` to ``path`` atomically: a partially written summary or
    manifest must never be observable (a killed process leaves either the old
    file or the fully new one, never a truncated checkpoint). os.replace is
    atomic within a filesystem, and the temp file is co-located with the
    target so it always is. The temp name carries the pid so two processes
    writing the same directory cannot clobber each other's staging file, and
    the parent directory is fsynced after the rename so the rename itself is
    durable across a machine/filesystem crash (not just a process kill)."""
    directory = os.path.dirname(path) or "."
    tmp = "%s.tmp.%d" % (path, os.getpid())
    with open(tmp, "w") as fw:
        fw.write(text)
        fw.flush()
        os.fsync(fw.fileno())
    os.replace(tmp, path)
    try:
        dfd = os.open(directory, os.O_DIRECTORY)
        try:
            os.fsync(dfd)
        finally:
            os.close(dfd)
    except (OSError, AttributeError):
        pass                                    # best-effort (e.g. non-POSIX)


def write_summary(path, rows):
    lines = [_SUMMARY_HEADER]
    for r in rows:
        lines.append("%d %.12g %s %s %.6f %.6f %d %d %d\n" % (
            r["idx"], r["T"], r["status"], r["error_stage"],
            r["re"], r["im"], r["match"],
            r["flex_converged"], r["flex_iter"]))
    _atomic_write_text(path, "".join(lines))


def read_summary_rows(path):
    """Parse a summary file back into row dicts (inverse of write_summary).

    Returns only the clean, contiguous prefix starting at idx 0: the FIRST
    malformed/short line, or a row whose index is not the next expected one
    (a duplicate or out-of-order index), TRUNCATES the result. A damaged
    checkpoint must shorten the reusable prefix, never be silently patched over
    with last-row-wins semantics that could resurrect a stale/incompatible
    rung."""
    rows = []
    if not os.path.exists(path):
        return rows
    expected = 0
    for line in open(path):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        c = s.split()
        if len(c) < 9:
            break
        try:
            row = dict(
                idx=int(c[0]), T=float(c[1]), status=c[2], error_stage=c[3],
                re=float(c[4]), im=float(c[5]), match=int(c[6]),
                flex_converged=int(c[7]), flex_iter=int(c[8]))
        except ValueError:
            break
        if row["idx"] != expected:              # duplicate / out-of-order
            break
        rows.append(row)
        expected += 1
    return rows


# Every interaction input the k-space reader accepts (QLMSkInput.valid_namelist:
# the combined ``Coulomb`` and the external field ``Extern`` are easy to miss).
# Changing any of their filenames OR contents must invalidate a resume.
_INTERACTION_KEYS = ("Geometry", "Transfer", "Coulomb", "CoulombIntra",
                     "CoulombInter", "Hund", "Exchange", "Ising", "PairLift",
                     "PairHop", "Extern", "Interall", "Initial")

# Eliashberg keys that vary per rung or are purely operational (output paths,
# the continuation seed) and therefore must NOT enter the physics fingerprint.
_ELI_OPERATIONAL_KEYS = ("seed_eigenvector", "output_eigenvalue", "output_gap",
                         "path_to_flex_output")


def _file_digest(path):
    """SHA-256 of a file's bytes, or a stable sentinel if it cannot be read.
    Used so editing an interaction file in place invalidates a resume."""
    try:
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(65536), b""):
                h.update(chunk)
        return h.hexdigest()
    except OSError:
        return "unreadable:" + os.path.basename(str(path))


def config_fingerprint(base, run_eli, base_dir="."):
    """A stable hash of the shape/physics-relevant configuration (everything
    that must NOT change between the original sweep and a resume, except the
    per-rung temperature). File existence alone does not prove a rung is
    reusable -- a resume against a differently-shaped/parameterised config must
    fail fast rather than mix incompatible results.

    The fingerprint covers the FULL ``[mode.param]`` (minus the per-rung ``T``),
    the FULL ``[eliashberg]`` block (minus per-rung/operational keys), the
    continuation seeding flags, the interaction file *names*, and -- crucially
    -- a content digest of each interaction file, so an in-place edit of a
    Geometry/Transfer/Coulomb/... file also invalidates the resume."""
    param = dict(base.get("mode", {}).get("param", {}))
    param.pop("T", None)                       # per-rung; excluded
    eli = {k: v for k, v in base.get("eliashberg", {}).items()
           if k not in _ELI_OPERATIONAL_KEYS}
    cont = base.get("continuation", {})
    fin = base.get("file", {}).get("input", {})
    inter = fin.get("interaction", {})
    ipath = _abspath(base_dir,
                     inter.get("path_to_input", fin.get("path_to_input", ".")))
    digests = {}
    for key in _INTERACTION_KEYS:
        fn = inter.get(key)
        if fn:
            digests[key] = _file_digest(_abspath(ipath, fn))
    src = {
        "mode": base.get("mode", {}).get("mode"),
        "param": param,                        # full mode.param minus T
        "eliashberg": eli,                     # full eliashberg minus operational
        "run_eliashberg": bool(run_eli),
        "warm_start": cont.get("warm_start", True),
        "seed_gap": cont.get("seed_gap", True),
        # names identify the files; digests catch in-place content edits.
        "interaction_files": {k: inter.get(k) for k in _INTERACTION_KEYS},
        "interaction_digests": digests,
    }
    blob = json.dumps(src, sort_keys=True, default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def _manifest_path(out_dir):
    return os.path.join(out_dir, MANIFEST_NAME)


def write_manifest(out_dir, ladder, fingerprint):
    _atomic_write_text(_manifest_path(out_dir), json.dumps(
        {"version": _MANIFEST_VERSION, "ladder": list(ladder),
         "fingerprint": fingerprint}, indent=2) + "\n")


def load_manifest(out_dir):
    path = _manifest_path(out_dir)
    if not os.path.exists(path):
        return None
    try:
        with open(path) as f:
            return json.load(f)
    except (ValueError, OSError):
        return None


def _validate_resume(manifest, ladder, fingerprint):
    """Fail fast unless the existing sweep matches this run's ladder and
    physics configuration; a silent mix of incompatible rungs is worse than
    a clear refusal."""
    if manifest is None:
        raise ValueError(
            "resume requested but no %s was found in the output directory; "
            "cannot prove existing rungs belong to this sweep. Run without "
            "resume to start fresh." % MANIFEST_NAME)
    if not isinstance(manifest, dict) or manifest.get("version") != _MANIFEST_VERSION:
        raise ValueError(
            "resume manifest version mismatch: expected version %d. Start "
            "fresh in a new output_dir or regenerate the sweep." %
            _MANIFEST_VERSION)
    m_ladder = [float(t) for t in manifest.get("ladder", [])]
    # Reject non-finite entries explicitly: NaN comparisons are always False,
    # so a NaN manifest temperature would otherwise slip through the tolerance
    # check and be treated as matching a real ladder value.
    bad = (len(m_ladder) != len(ladder)
           or any(not math.isfinite(a) or not math.isfinite(b)
                  for a, b in zip(m_ladder, ladder))
           or any(abs(a - b) > 1e-12 * max(1.0, abs(b))
                  for a, b in zip(m_ladder, ladder)))
    if bad:
        raise ValueError(
            "resume ladder mismatch: the existing sweep ran %r but this run "
            "resolves to %r. Ladders must be identical (and finite) to "
            "resume." % (m_ladder, list(ladder)))
    if manifest.get("fingerprint") != fingerprint:
        raise ValueError(
            "resume configuration mismatch: the shape/physics fingerprint of "
            "this run differs from the recorded sweep (CellShape/SubShape/"
            "Nmat/filling/Ncond/interaction/eliashberg). Resume only a sweep "
            "with an identical configuration, or start fresh in a new "
            "output_dir.")


def _eig_parseable(rout, eig_name):
    try:
        parse_leading_eig(os.path.join(rout, eig_name))
        return True
    except (ValueError, OSError):
        return False


def _npz_parseable(path, key):
    """Return whether an NPZ seed can be opened and contains readable data."""
    try:
        with np.load(path, allow_pickle=False) as data:
            value = data[key]
            # Force NumPy/zipfile to read and CRC-check the payload instead of
            # merely accepting a valid archive header.
            return value.size > 0
    except (OSError, ValueError, KeyError, EOFError):
        return False


def _resume_scan(out_dir, ladder, old_by_idx, run_eli, warm_start,
                 seed_gap_on, sigma_name, gap_name, eig_name):
    """Longest contiguous prefix of completed rungs.

    A rung qualifies only when its recorded summary row is present and
    non-error, its temperature matches the ladder, its eigenvalue.dat is
    parseable (for an Eliashberg sweep), AND every output the rung was
    configured to produce (sigma when ``warm_start``, gap when ``seed_gap_on``)
    is present and parseable on disk. This on-disk validation applies to EVERY
    completed rung -- including the final ladder rung -- so a fully resumed
    sweep can never report success while retaining a missing/corrupt output.
    Returns ``(prefix_rows, prev, start_idx)`` where ``prev`` is the seed
    source (always a validated rung) for ``ladder[start_idx]``."""
    prefix_rows, prev, start_idx = [], None, 0
    for idx, T in enumerate(ladder):
        row = old_by_idx.get(idx)
        rout = os.path.join(rung_dir(out_dir, idx, T), "output")
        completed = (
            row is not None
            and row.get("status") in ("ok", "not_converged")
            and abs(row["T"] - T) <= 1e-12 * max(1.0, abs(T))
            and (not run_eli or _eig_parseable(rout, eig_name))
            # the rung's own promised seed outputs must exist+parse, or it did
            # not truly complete -> recompute from here (final rung included).
            and _seed_files_present(rout, warm_start, seed_gap_on,
                                    sigma_name, gap_name, quiet=True))
        if not completed:
            break
        prefix_rows.append(row)
        start_idx = idx + 1
        prev = rout
    return prefix_rows, prev, start_idx


def _abspath(base_dir, p):
    return p if os.path.isabs(p) else os.path.normpath(os.path.join(base_dir, p))


def run(input_dict, base_dir=".", keep_going=False, dry_run=False,
        resume=False):
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
    # Gap seeding requires the Eliashberg solver to actually run: a FLEX-only
    # sweep (run_eliashberg=false) never writes gap_dynamic.npz, even if a
    # leftover [eliashberg] block still says frequency="dynamic". Without the
    # run_eli guard, resume would demand a gap file no rung can produce (so no
    # rung is ever reusable) and normal continuation would cold-start every rung.
    seed_gap_on = (run_eli and cont.get("seed_gap", True)
                   and eliashberg_frequency(input_dict) == "dynamic")
    out_dir = _abspath(base_dir, cont.get("output_dir", "tsweep"))
    summary_file = cont.get("summary_file", "lambda_vs_T.dat")
    resume = resume or bool(cont.get("resume", False))
    # eigenvalue.dat name (per-rung override honored via output_eigenvalue); the
    # resume scan needs it before make_rung_dicts is called.
    eig_name = input_dict.get("eliashberg", {}).get(
        "output_eigenvalue", "eigenvalue.dat")
    # Fingerprint the shape/physics config BEFORE resolving input paths to
    # absolute, so it is stable across machines/CWD (see config_fingerprint);
    # base_dir lets it still locate + content-hash the interaction files.
    fingerprint = config_fingerprint(input_dict, run_eli, base_dir=base_dir)

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

    os.makedirs(out_dir, exist_ok=True)
    summary_path = os.path.join(out_dir, summary_file)

    # Resume vs fresh. Resume validates the recorded ladder+fingerprint (fail
    # fast on mismatch), then skips the longest contiguous prefix of completed,
    # seedable rungs and restarts at the first incomplete one -- seeded from the
    # last valid rung. A fresh run records the manifest and warns (but does not
    # refuse) if it is about to overwrite an existing sweep.
    rows, prev, start_idx = [], None, 0            # prev = seedable prev rung dir or None
    if resume:
        _validate_resume(load_manifest(out_dir), ladder, fingerprint)
        old_by_idx = {r["idx"]: r for r in read_summary_rows(summary_path)}
        rows, prev, start_idx = _resume_scan(
            out_dir, ladder, old_by_idx, run_eli, warm_start, seed_gap_on,
            sigma_name, gap_name, eig_name)
        if start_idx >= len(ladder):
            logger.info("resume: all %d rungs already complete; nothing to do.",
                        len(ladder))
        else:
            logger.info("resume: %d completed rung(s) reused; restarting at "
                        "rung %d (T=%g)%s", start_idx, start_idx,
                        ladder[start_idx],
                        "" if prev is None else " seeded from " + prev)
        write_summary(summary_path, rows)          # normalize/refresh checkpoint
    else:
        if load_manifest(out_dir) is not None or os.path.exists(summary_path):
            logger.warning("output_dir %s already contains a sweep; a fresh "
                           "run overwrites it rung-by-rung. Pass resume=true "
                           "(--resume) to continue the existing sweep instead.",
                           out_dir)
        # Invalidate the old checkpoint before installing the new manifest.
        # Otherwise a kill before the first completed rung could leave a
        # new-config manifest beside old-config rows and outputs, which a later
        # resume would incorrectly reuse.
        write_summary(summary_path, rows)
        write_manifest(out_dir, ladder, fingerprint)

    for idx, T in enumerate(ladder):
        if idx < start_idx:
            continue
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
            write_summary(summary_path, rows)
            prev = rout
            continue
        flex_dict, eli_dict = make_rung_dicts(
            input_dict, T, rout, run_eli, sigma_init=sigma_init, seed_gap=seed)
        eig_name = eli_dict["eliashberg"].get("output_eigenvalue", "eigenvalue.dat") \
            if eli_dict is not None else "eigenvalue.dat"
        # Remove this rung's expected checkpoint/seed outputs before running it
        # so a later "file present" check reflects THIS execution -- a solver
        # that returns success without rewriting a file must not let a stale
        # output from a previous sweep be validated and chained as a seed.
        for fn in (sigma_name if warm_start else None,
                   gap_name if seed_gap_on else None,
                   eig_name if run_eli else None):
            if fn:
                try:
                    os.remove(os.path.join(rout, fn))
                except OSError:
                    pass
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
        # Atomically checkpoint after every rung so an interruption leaves a
        # valid (never truncated) summary that a later --resume can read.
        write_summary(summary_path, rows)
        # seedable iff not error and required files exist
        ok_seed = row["status"] != "error" and _seed_files_present(
            rout, warm_start, seed_gap_on, sigma_name, gap_name)
        if row["status"] == "error" and not keep_going:
            break
        prev = rout if ok_seed else None
        if row["status"] == "error":               # keep_going: cold-start next
            logger.warning("continuing after error at T=%g; next rung cold-started", T)

    write_summary(summary_path, rows)
    return rows


def _seed_files_present(rout, warm_start, seed_gap_on, sigma_name, gap_name,
                        quiet=False):
    if warm_start and not _npz_parseable(os.path.join(rout, sigma_name),
                                         "sigma"):
        if not quiet:
            logger.warning("expected sigma %s missing or corrupt in %s; next "
                           "rung cold-starts", sigma_name, rout)
        return False
    if seed_gap_on and not _npz_parseable(os.path.join(rout, gap_name), "gap"):
        if not quiet:
            logger.warning("expected gap %s missing or corrupt in %s; next "
                           "rung cold-starts", gap_name, rout)
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
    parser.add_argument("--resume", action="store_true",
                        help="skip completed rungs of an existing sweep and "
                             "restart at the first incomplete one (validates "
                             "the recorded ladder/configuration)")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    with open(args.input, "rb") as f:
        input_dict = tomli.load(f)
    run(input_dict, base_dir=os.path.dirname(os.path.abspath(args.input)),
        keep_going=args.keep_going, dry_run=args.dry_run, resume=args.resume)


if __name__ == "__main__":
    main()
