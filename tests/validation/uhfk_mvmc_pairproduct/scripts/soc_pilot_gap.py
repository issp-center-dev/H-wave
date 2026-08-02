"""Run H-wave UHFk on a case dir and report the HOMO-LUMO gap at each
candidate ``Ncond`` value.

The SCF is driven with whichever ``Ncond`` the fixture's ``input.toml``
carries. The sorted eigenvalues from ``eigen.npz`` are then used to
compute a gap for every requested candidate ``Ncond``: this is a
non-interacting-style band-structure gap check that identifies which
fillings sit inside a spectral gap of the converged Hamiltonian.

Marker ``OK`` requires the candidate ``Ncond`` to be even (Ne must fit
into the AntiParallel /
spin-mixed sector) AND the gap to be at least ``1e-3``. Any other
outcome is ``REJECT``.
"""
from __future__ import annotations

import argparse
import os
import sys
import numpy as np


def _shim_numpy_float():
    """H-wave still uses ``np.float`` in a few legacy paths; provide a
    shim so the SCF loop does not crash on newer NumPy."""
    if not hasattr(np, "float"):
        np.float = float  # noqa: NPY001


def run_and_report(case_dir, ncond_candidates):
    """Run H-wave in ``case_dir``, then print HOMO/LUMO/gap per candidate."""
    case_dir = os.path.abspath(case_dir)
    orig_cwd = os.getcwd()
    os.chdir(case_dir)
    try:
        _shim_numpy_float()
        sys.argv = ["hwave", "input.toml"]
        from hwave.qlms import main
        main()
    finally:
        os.chdir(orig_cwd)

    eig_path = os.path.join(case_dir, "output", "eigen.npz")
    eig = np.load(eig_path)
    ws = np.sort(np.real(np.asarray(eig["eigenvalue"]).reshape(-1)))

    any_ok = False
    for ncond in ncond_candidates:
        if ncond <= 0 or ncond >= ws.size:
            print(f"  Ncond={ncond:3d}: out of range (spectrum has "
                  f"{ws.size} eigenvalues)")
            continue
        homo, lumo = ws[ncond - 1], ws[ncond]
        gap = lumo - homo
        marker = "OK" if (ncond % 2 == 0 and gap >= 1e-3) else "REJECT"
        if marker == "OK":
            any_ok = True
        print(f"  Ncond={ncond:3d}: HOMO={homo: .6e} "
              f"LUMO={lumo: .6e} gap={gap: .3e} [{marker}]")
    return any_ok


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True,
                    help="fixture directory containing input.toml")
    ap.add_argument("--ncond", default="4,6,8,10,12",
                    help="comma-separated candidate Ncond values")
    args = ap.parse_args()
    ok = run_and_report(
        args.case, [int(x) for x in args.ncond.split(",")],
    )
    sys.exit(0 if ok else 2)
