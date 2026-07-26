"""Regenerate the LEGACY (pre-bond-channels) static-Eliashberg golden files.

The golden ``eigenvalue.dat`` / ``gap.dat`` committed next to this script were
produced **at commit 712b1a0** -- the merge commit this branch
(``feature/bond-channels-eliashberg``) was cut from, i.e. the last state of the
solver before the bond-channel work existed.  ``tests/test_sc_legacy_golden.py``
asserts that HEAD reproduces them **byte for byte** when neither
``bond_channels`` nor ``spectral_shift`` is present in the input, which is the
PR's critical invariant: this is research software, so an existing user's output
file must not move.

Regeneration recipe (exactly what was run)::

    git worktree add /tmp/hwave-712b1a0 712b1a0
    cd /tmp/hwave-712b1a0
    PYTHONPATH=/tmp/hwave-712b1a0/src python \
        <this file> --outdir tests/sc/legacy_golden
    git worktree remove /tmp/hwave-712b1a0

``chi0q.npz`` is committed alongside the goldens and is *loaded*
(``chi0q_mode="load"``), never recomputed: ``_calc_chi0q_internal`` is not
bit-reproducible run to run (~1e-14 ULP jitter in the FFT/BLAS kernels; issue
#85), while everything downstream of chi0q is.  Pinning it is what makes a
byte-level golden possible at all.  It was produced once with
``--make-chi0q``; regenerating it invalidates the goldens.
"""
import argparse
import os
import sys

import numpy as np


# The k-grid/Matsubara grid are deliberately tiny so the golden run is a
# fraction of a second, and single-orbital so gap.dat stays small.
RPA_INPUT = "tests/rpa/input"


def build_input(outdir, chi0q_mode="load"):
    """The legacy configuration the goldens were produced from.

    Nothing bond-specific and no ``spectral_shift``: this is exactly what a
    pre-branch user's input looked like.  ``solver_mode="iteration"`` is the
    path that grew the ``eigenvalue_note`` comment lines, so it is the one the
    byte-identity claim has to cover.
    """
    return {
        "mode": {"param": {"T": 1.0, "CellShape": [4, 4, 1],
                           "SubShape": [1, 1, 1], "Nmat": 128,
                           "filling": 0.5}},
        "file": {
            "input": {"interaction": {"path_to_input": RPA_INPUT,
                                      "Geometry": "geom.dat",
                                      "Transfer": "transfer.dat",
                                      "CoulombIntra": "coulombintra.dat"}},
            "output": {"path_to_output": outdir},
        },
        "eliashberg": {"chi0q_mode": chi0q_mode,
                       "solver_mode": "iteration",
                       "max_iter": 40,
                       "convergence_tol": 1.0e-8},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--make-chi0q", action="store_true",
                        help="(re)compute chi0q.npz instead of running the "
                             "solver; INVALIDATES the committed goldens")
    args = parser.parse_args()

    import hwave.sc as sc
    print("hwave.sc from:", sc.__file__, file=sys.stderr)

    os.makedirs(args.outdir, exist_ok=True)
    if args.make_chi0q:
        chi0q = sc._calc_chi0q_internal(build_input(args.outdir,
                                                    chi0q_mode="calc"))
        np.savez(os.path.join(args.outdir, "chi0q.npz"), chi0q=chi0q)
        print("wrote chi0q.npz", chi0q.shape, file=sys.stderr)
        return
    sc.calc_eliashberg(build_input(args.outdir))


if __name__ == "__main__":
    main()
