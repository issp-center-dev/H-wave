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
"""
import argparse
import logging
import sys

_USAGE = ("hwave_sigma_split total.npz out.npz "
          "(--static static.npz [--uhfk-spin-major] | --zero-static | "
          "--uhfk-trans-mod trans_mod.npz --bare-transfer transfer.dat) [--force]")


def _parser():
    p = argparse.ArgumentParser(prog="hwave_sigma_split", usage=_USAGE,
                                description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("total", help="total-form sigma archive (no sigma_convention or \"total\")")
    p.add_argument("out", help="output split archive")
    p.add_argument("--static", metavar="static.npz",
                   help="archive with the static correction (member sigma_static; sigma as fallback)")
    p.add_argument("--zero-static", action="store_true", help="write sigma_static = 0")
    p.add_argument("--uhfk-trans-mod", metavar="trans_mod.npz",
                   help="native UHFk trans_mod archive (requires --bare-transfer)")
    p.add_argument("--bare-transfer", metavar="transfer.dat",
                   help="the bare Transfer input of the run (Wannier90-style text or .npz)")
    p.add_argument("--uhfk-spin-major", action="store_true",
                   help="the --static array is spin-major (nvol, 2 norb, 2 norb)")
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
    from hwave.solver.flex_hf import sigma_split_convert
    try:
        summary = sigma_split_convert(
            args.total, args.out, static_path=args.static, zero_static=args.zero_static,
            uhfk_trans_mod=args.uhfk_trans_mod, bare_transfer=args.bare_transfer,
            uhfk_spin_major=args.uhfk_spin_major, force=args.force, logger=log)
    except (ValueError, FileExistsError, OSError, KeyError) as exc:
        log.error("%s", exc)
        parser.print_usage(sys.stderr)
        return 1
    log.info("wrote %s (nmat=%d, nvol=%d, norb=%d; max|sigma_static| = %.3e, max|sigma_fluct| = %.3e)",
             summary["out"], summary["nmat"], summary["nvol"], summary["norb"],
             summary["static_max"], summary["fluct_max"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
