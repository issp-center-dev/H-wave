"""``hwave_sigma_split``: convert a total-form FLEX self-energy archive into
the ``"split"`` form (``sigma = sigma_static + sigma_fluct``) that a run with
``flex_hartree_fock = true`` accepts as ``sigma_init``.

    hwave_sigma_split total.npz out.npz --static static.npz [--uhfk-spin-major]
    hwave_sigma_split total.npz out.npz --zero-static
    hwave_sigma_split total.npz out.npz --uhfk-trans-mod trans_mod.npz --bare-transfer transfer.dat

The static correction is the interaction-only mean field in the FLEX
reciprocal-space layout. With ``--uhfk-trans-mod`` it is built from a
native UHFk ``trans_mod`` archive: the archive is Fourier-transformed to
the FLEX C-order momentum layout with the solver's convention
(``e^{+ikR}``), the bare transfer assembled from ``--bare-transfer``
(the same Transfer input the FLEX run reads) is subtracted, the result is
validated spin-diagonal with equal spin blocks and reduced to the
spin-free block. External one-body fields are not supported by
``flex_hartree_fock`` and are therefore not part of the recipe.

The total self-energy is symmetrized to
``(Sigma(k, i w) + Sigma(k, -i w)^dag) / 2`` before the split, and the
relative deviation removed is reported. An archive whose deviation exceeds
``--asymmetry-tol`` (default 1e-3) is refused: a FLEX run without the
Hermitian projection can converge to a symmetry-broken, non-causal branch.

    hwave_sigma_split --check sigma.npz [--asymmetry-tol TOL]

only reports the relative deviation of an existing archive (exit status 1
when it exceeds the tolerance); use it to audit archives written before the
projection existed.
"""
import argparse
import logging
import sys

_USAGE = ("hwave_sigma_split total.npz out.npz "
          "(--static static.npz [--uhfk-spin-major] | --zero-static | "
          "--uhfk-trans-mod trans_mod.npz --bare-transfer transfer.dat) [--asymmetry-tol TOL] "
          "[--force]\n       hwave_sigma_split --check sigma.npz [--asymmetry-tol TOL]")


def _parser():
    p = argparse.ArgumentParser(prog="hwave_sigma_split", usage=_USAGE,
                                description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("total", help="total-form sigma archive (no sigma_convention or \"total\")")
    p.add_argument("out", nargs="?", help="output split archive (not with --check)")
    p.add_argument("--check", action="store_true",
                   help="only report the relative deviation max|Sigma(k,iw) - Sigma(k,-iw)^dag| / "
                        "max|Sigma| of the archive; exit status 1 above --asymmetry-tol")
    p.add_argument("--static", metavar="static.npz",
                   help="archive with the static correction (member sigma_static; sigma as fallback)")
    p.add_argument("--zero-static", action="store_true", help="write sigma_static = 0")
    p.add_argument("--uhfk-trans-mod", metavar="trans_mod.npz",
                   help="native UHFk trans_mod archive (requires --bare-transfer)")
    p.add_argument("--bare-transfer", metavar="transfer.dat",
                   help="the bare Transfer input of the run (Wannier90-style text or .npz)")
    p.add_argument("--uhfk-spin-major", action="store_true",
                   help="the --static array is spin-major (nvol, 2 norb, 2 norb)")
    p.add_argument("--asymmetry-tol", metavar="TOL", type=float, default=None,
                   help="largest relative deviation max|Sigma(k,iw) - Sigma(k,-iw)^dag| / max|Sigma| "
                        "that is symmetrized away (default 1e-3); a larger one is refused")
    p.add_argument("--force", action="store_true", help="overwrite an existing output")
    return p


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    log = logging.getLogger("hwave_sigma_split")
    parser = _parser()
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:
        return 0 if exc.code == 0 else 1
    from hwave.solver.flex_hf import sigma_split_convert, SIGMA_SPLIT_ASYMMETRY_TOL
    tol = SIGMA_SPLIT_ASYMMETRY_TOL if args.asymmetry_tol is None else args.asymmetry_tol
    if not tol >= 0.0:
        log.error("--asymmetry-tol must be a non-negative number, got %r", args.asymmetry_tol)
        parser.print_usage(sys.stderr)
        return 1
    if args.check:
        extra = [o for o, v in (("out", args.out), ("--static", args.static),
                                ("--zero-static", args.zero_static),
                                ("--uhfk-trans-mod", args.uhfk_trans_mod),
                                ("--bare-transfer", args.bare_transfer),
                                ("--uhfk-spin-major", args.uhfk_spin_major),
                                ("--force", args.force)) if v]
        if extra:
            log.error("--check takes the archive and --asymmetry-tol only (got %s)", ", ".join(extra))
            parser.print_usage(sys.stderr)
            return 1
        from hwave.solver.flex_hf import sigma_hermitian_check
        try:
            res = sigma_hermitian_check(args.total, logger=log)
        except (ValueError, OSError, KeyError) as exc:
            log.error("%s", exc)
            return 1
        ok = res["asymmetry"] <= tol
        (log.info if ok else log.warning)(
            "%s: relative deviation max|Sigma(k,iw) - Sigma(k,-iw)^dag| / max|Sigma| = %.3e "
            "(tolerance %.1e, sigma %s, %s grid)%s", args.total, res["asymmetry"], tol,
            tuple(res["shape"]), res["grid"],
            "" if ok else "; not the self-energy of a converged Hermitian run: repeat it with "
            "flex_sigma_hermitian_projection = true")
        return 0 if ok else 1
    if args.out is None:
        log.error("the output archive is required (or use --check)")
        parser.print_usage(sys.stderr)
        return 1
    try:
        summary = sigma_split_convert(
            args.total, args.out, static_path=args.static, zero_static=args.zero_static,
            uhfk_trans_mod=args.uhfk_trans_mod, bare_transfer=args.bare_transfer,
            uhfk_spin_major=args.uhfk_spin_major, force=args.force, logger=log,
            asymmetry_tol=tol)
    except (ValueError, FileExistsError, OSError, KeyError) as exc:
        log.error("%s", exc)
        parser.print_usage(sys.stderr)
        return 1
    log.info("wrote %s (nmat=%d, nvol=%d, norb=%d; max|sigma_static| = %.3e, max|sigma_fluct| = %.3e; "
             "relative Hermitian deviation removed %.3e)",
             summary["out"], summary["nmat"], summary["nvol"], summary["norb"],
             summary["static_max"], summary["fluct_max"], summary["asymmetry"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
