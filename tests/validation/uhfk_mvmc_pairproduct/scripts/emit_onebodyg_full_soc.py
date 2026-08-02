"""Emit OneBodyG.dat with the full Cartesian product of (i, s, j, t).

Format: 5 header lines (skipped by ``_read_green``'s ``np.loadtxt(...,
skiprows=5)``), followed by rows of ``i s j t`` where ``i, j`` are
packed site indices in ``[0, Ns)`` and ``s, t in {0, 1}``. In SOC mode
``_save_greenone`` emits every ``(s, t)`` combination, so requesting the
full Cartesian product is the coverage-maximising choice.

Usage::

    python emit_onebodyg_full_soc.py --cell 4,4,1 --out OneBodyG.dat
"""
from __future__ import annotations

import argparse


def emit(cell, out):
    lx, ly, lz = cell
    Ns = lx * ly * lz
    rows = []
    for i in range(Ns):
        for s in (0, 1):
            for j in range(Ns):
                for t in (0, 1):
                    rows.append((i, s, j, t))

    header = [
        "===============================",
        "NCisAjs         {}".format(len(rows)),
        "===============================",
        "======== Green functions ======",
        "===============================",
    ]
    with open(out, "w") as fw:
        for line in header:
            fw.write(line + "\n")
        for (i, s, j, t) in rows:
            fw.write("    {:>3d}     {:>1d}     {:>3d}     {:>1d}\n".format(
                i, s, j, t))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True, help="lx,ly,lz")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    emit(tuple(int(x) for x in args.cell.split(",")), args.out)
