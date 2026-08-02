"""Emit a real-space Transfer.dat with NN hopping + Rashba SOC.

Format: H-wave wannier90-like (parsed by ``qlmsio.wan90.read_w90``).
Layout::

    <header line>
    <num_wann placeholder>
    <num_entries>
    <ndegen row (15 per line, all 1s)>
    <data rows: irx iry irz iorb1 iorb2 re im>

Orbital indices are 1-based Wannier90 SOI packing ``2 * orb + spin + 1``.
With one physical orbital per unit cell, ``spin = 0`` (up) becomes file
index ``1``; ``spin = 1`` (down) becomes ``2``.

Real-space entries (per plan §Task11, docstring-derived and made
Hermitian; ``alpha`` = fixture Rashba coefficient)::

    NN hopping (spin-diagonal, both spins), for R in {+ex, -ex, +ey, -ey}:
        T_{up,up}(R) = -t
        T_{dn,dn}(R) = -t

    Rashba (spin off-diagonal, from
      H_R = alpha * sum_r [ i c^dag_{r,up} c_{r+ex,dn}
                          - i c^dag_{r,up} c_{r-ex,dn}
                          -   c^dag_{r,up} c_{r+ey,dn}
                          +   c^dag_{r,up} c_{r-ey,dn}
                          + h.c. ]):
        T_{up,dn}(+ex) = +i alpha
        T_{up,dn}(-ex) = -i alpha
        T_{up,dn}(+ey) =   -alpha
        T_{up,dn}(-ey) =   +alpha
        T_{dn,up}(R)   = conj(T_{up,dn}(-R))    (Hermitian conjugate)

Usage::

    python emit_rashba_transfer.py \\
        --cell 4,4,1 --t 1.0 --alpha 0.5 --out Transfer.dat
"""
from __future__ import annotations

import argparse


UP = 0
DN = 1


def _pack_index(orb: int, spin: int) -> int:
    """Wannier90 SOI packing: 1-based file index = 2*orb + spin + 1."""
    return 2 * orb + spin + 1


def emit(cell, t, alpha, out):
    """Write NN hopping + Rashba SOC entries to ``out`` in wannier90-like format.

    ``cell`` is retained in the signature only for interface parity with
    ``emit_onebodyg_full_soc``: the Transfer entries are indexed by the
    hopping vector ``R``, not by site, so the cell shape does not change
    the entry list. A ``cell`` value with any axis < 2 in a direction
    where the Rashba/NN entry lives would fold ``R = +1`` onto
    ``R = -1`` in H-wave's Hermite check; the primary use case (2D
    square with lx, ly >= 2 and lz = 1) is unaffected.
    """
    lx, ly, lz = cell
    entries = []

    def rec(rx, ry, rz, orb_a, spin_a, orb_b, spin_b, val):
        i1 = _pack_index(orb_a, spin_a)
        i2 = _pack_index(orb_b, spin_b)
        entries.append((rx, ry, rz, i1, i2, val.real, val.imag))

    # NN hopping (spin-diagonal), R in {+ex, -ex, +ey, -ey}, orb = 0 only.
    for R in [(+1, 0, 0), (-1, 0, 0), (0, +1, 0), (0, -1, 0)]:
        rec(R[0], R[1], R[2], 0, UP, 0, UP, complex(-t, 0.0))
        rec(R[0], R[1], R[2], 0, DN, 0, DN, complex(-t, 0.0))

    # Rashba SOC, up -> down block.
    a = float(alpha)
    rec(+1, 0, 0, 0, UP, 0, DN, complex(0.0, +a))
    rec(-1, 0, 0, 0, UP, 0, DN, complex(0.0, -a))
    rec(0, +1, 0, 0, UP, 0, DN, complex(-a, 0.0))
    rec(0, -1, 0, 0, UP, 0, DN, complex(+a, 0.0))

    # Rashba SOC, down -> up block. Hermitian conjugate of the up->down
    # block: T_{ba}(R) = conj(T_{ab}(-R)). For real alpha this yields:
    rec(+1, 0, 0, 0, DN, 0, UP, complex(0.0, +a))
    rec(-1, 0, 0, 0, DN, 0, UP, complex(0.0, -a))
    rec(0, +1, 0, 0, DN, 0, UP, complex(+a, 0.0))
    rec(0, -1, 0, 0, DN, 0, UP, complex(-a, 0.0))

    n = len(entries)
    with open(out, "w") as fw:
        fw.write("Transfer in wannier90-like format for uhfk (NN hopping + Rashba SOC)\n")
        fw.write("1\n")                # num_wann placeholder (unused by read_w90)
        fw.write("{}\n".format(n))
        ones = ["1"] * n
        # ndegen block: nr integers, 15 per line
        for i in range(0, n, 15):
            fw.write(" ".join(ones[i:i + 15]) + "\n")
        for (rx, ry, rz, i1, i2, re, im) in entries:
            fw.write("  {:4d} {:4d} {:4d} {:4d} {:4d}  {:.12f}  {:.12f}\n".format(
                rx, ry, rz, i1, i2, re, im))


def _self_check():
    """Sanity check: verify the emitted Transfer is Hermitian.

    H-wave's ``_check_hermite`` requires T_{o1,o2}(R) = conj(T_{o2,o1}(-R)).
    We emit to a temporary file, re-parse the data rows, and assert the
    condition holds exactly.
    """
    import os
    import tempfile

    with tempfile.NamedTemporaryFile(
            mode="w", suffix=".dat", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        emit((4, 4, 1), 1.0, 0.5, tmp_path)
        with open(tmp_path, "r") as fr:
            lines = fr.read().splitlines()
    finally:
        os.unlink(tmp_path)

    # Skip header (1) + num_wann (1) + nr (1) + ndegen rows.
    header_lines = 3
    n = int(lines[2])
    ndegen_rows = (n + 14) // 15
    data_lines = lines[header_lines + ndegen_rows:]
    table = {}
    for line in data_lines:
        toks = line.split()
        if len(toks) < 7:
            continue
        rx, ry, rz, i1, i2 = (int(toks[0]), int(toks[1]), int(toks[2]),
                              int(toks[3]) - 1, int(toks[4]) - 1)
        re, im = float(toks[5]), float(toks[6])
        table[(rx, ry, rz, i1, i2)] = complex(re, im)

    for (rx, ry, rz, o1, o2), v in table.items():
        partner = table.get((-rx, -ry, -rz, o2, o1), 0.0 + 0.0j)
        if abs(v - partner.conjugate()) > 1e-14:
            raise AssertionError(
                "non-Hermitian entry: R=({},{},{}), (o1,o2)=({},{}), "
                "T={}, conj(T_partner)={}".format(
                    rx, ry, rz, o1, o2, v, partner.conjugate()))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", help="lx,ly,lz")
    ap.add_argument("--t", type=float)
    ap.add_argument("--alpha", type=float)
    ap.add_argument("--out")
    ap.add_argument("--self-check", action="store_true",
                    help="run in-memory Hermiticity check and exit")
    args = ap.parse_args()
    if args.self_check:
        _self_check()
        print("self-check passed: emitted Transfer is Hermitian.")
    else:
        missing = [name for name, value in
                   [("--cell", args.cell), ("--t", args.t),
                    ("--alpha", args.alpha), ("--out", args.out)]
                   if value is None]
        if missing:
            ap.error("missing required arguments: " + ", ".join(missing))
        emit(tuple(int(x) for x in args.cell.split(",")),
             args.t, args.alpha, args.out)
