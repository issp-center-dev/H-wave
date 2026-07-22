#!/usr/bin/env python3
"""Compare H-wave and ComplexUHF outputs for an APBC validation case."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _read_energy_hwave(path: Path) -> float:
    for ln in path.read_text().splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        # Looking for "Energy_Total = <value>"
        if "Energy_Total" in ln:
            return float(ln.split("=")[-1].strip())
    raise ValueError(f"Energy_Total not found in {path}")


def _read_energy_complexuhf(path: Path) -> float:
    # ComplexUHF "zvo_result.dat" emits two lines:
    #     " energy <value> "
    #     " num    <value> "
    # (mVMC-1.4.0 src/ComplexUHF/output.c:64-65). Find the "energy"
    # keyword line and return the following whitespace-separated token.
    for ln in path.read_text().splitlines():
        toks = ln.split()
        if len(toks) >= 2 and toks[0] == "energy":
            return float(toks[1])
    raise ValueError(f"no 'energy' line in {path}")


def _read_greenone(path: Path) -> dict[tuple[int, int, int, int], complex]:
    """Parse a 6-column one-body Green output ``(i, s, j, t, re, im)``.

    Handles BOTH the non-SOC (Sz-diagonal, off-block-diagonal entries
    are zero) and the SOC (Rashba cross-spin ``s != t`` non-zero
    entries) layouts. Both mVMC-1.4.0 ComplexUHF's ``zvo_UHF_cisajs.dat``
    (mVMC-1.4.0 ``src/ComplexUHF/output.c:99``, ``print("%d %d %d %d
    %.10lf %.10lf")``) and H-wave UHFk's ``greenone.dat`` write the
    same layout, so a single parser reads both sides of the cross-check.

    Under SOC + APBC + SubShape > [1, 1, 1], the
    ComplexUHF UHF binary emits non-zero cross-spin ``s != t`` rows
    from the Rashba SOC term; the parser preserves them verbatim so
    the caller can compare complex cross-spin entries against the
    bridge's shipping A density at 1e-6 tolerance.
    """
    out: dict[tuple[int, int, int, int], complex] = {}
    for ln in path.read_text().splitlines():
        parts = ln.split()
        if len(parts) < 6:
            continue
        try:
            i, s, j, t = (int(parts[k]) for k in range(4))
            re, im = float(parts[4]), float(parts[5])
        except ValueError:
            continue
        out[(i, s, j, t)] = complex(re, im)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--hwave", required=True, type=Path)
    ap.add_argument("--complexuhf", required=True, type=Path)
    ap.add_argument("--case", required=True)
    ap.add_argument("--etol", type=float, default=1.0e-6, help="energy rel tol")
    ap.add_argument("--gtol", type=float, default=1.0e-6, help="green abs tol")
    args = ap.parse_args()

    # --- Energy
    e_h = _read_energy_hwave(args.hwave / "energy.dat")
    # ComplexUHF output filename: try common candidates. ComplexUHF (UHF
    # binary, not VMC) writes "<CDataFileHead>_result.dat", which with the
    # vmcdry default is "zvo_result.dat".
    for cand in ("zvo_result.dat", "zvo_energy.dat", "uhf_zvo_energy.dat",
                 "energy.dat"):
        e_path = args.complexuhf / cand
        if e_path.exists():
            e_c = _read_energy_complexuhf(e_path)
            break
    else:
        print(f"[{args.case}] complexuhf energy file not found; "
              "tried zvo_result.dat, zvo_energy.dat etc.")
        return 1
    rel = abs(e_h - e_c) / max(abs(e_c), 1.0e-12)
    print(f"[{args.case}] energy: hwave={e_h:.10e}  complexuhf={e_c:.10e}  rel={rel:.2e}  "
          f"{'OK' if rel < args.etol else 'MISMATCH'}")

    # --- One-body Green (optional; skip cleanly if files differ in format)
    h_green_path = args.hwave / "greenone.dat"
    c_green_path = None
    # ComplexUHF writes "<CDataFileHead>_UHF_cisajs.dat" (output.c:99).
    for cand in ("zvo_UHF_cisajs.dat", "zvo_cisajs.dat", "uhf_cisajs.dat"):
        if (args.complexuhf / cand).exists():
            c_green_path = args.complexuhf / cand
            break

    g_ok = True
    g_msg = "not compared (file missing on one side)"
    if h_green_path.exists() and c_green_path is not None:
        g_h = _read_greenone(h_green_path)
        g_c = _read_greenone(c_green_path)
        keys_both = set(g_h) & set(g_c)
        if not keys_both:
            g_msg = "no common (i,s,j,t) keys"
            g_ok = False
        else:
            max_diff = max(abs(g_h[k] - g_c[k]) for k in keys_both)
            g_ok = max_diff <= args.gtol
            g_msg = f"max |delta| = {max_diff:.3e} over {len(keys_both)} pairs"
    print(f"[{args.case}] greenone: {'OK' if g_ok else 'MISMATCH'}  {g_msg}")

    return 0 if (rel < args.etol and g_ok) else 1


if __name__ == "__main__":
    sys.exit(main())
