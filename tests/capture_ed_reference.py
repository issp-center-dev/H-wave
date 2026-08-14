#!/usr/bin/env python3
"""Capture reference values of tests/test_rpa_vs_ed_oracle.py BEFORE the
extraction refactor (#151 harness).  Run once from the repo root:
    python tests/capture_ed_reference.py
The extraction gate (test_ed_oracle_extraction_gate.py) pins the refactored
utility against this file at round-off, so a coordinated utility/map change
cannot slip through compensating errors."""
import numpy as np

import tests.test_rpa_vs_ed_oracle as oracle

out = {}
out["h1"] = oracle.H1
for kind in ("CoulombIntra", "CoulombInter", "Exchange", "PairHop"):
    out["hf_" + kind] = oracle._hf_h1(kind, oracle.V1)
out["base_chi"] = oracle._chi_ed()
for kind in ("CoulombIntra", "CoulombInter"):
    out["first_" + kind] = oracle._ed_first_order(kind)[0]
# kernel sample: rebuild the K matrix exactly as _chi_ed does
H = np.zeros((oracle.DIM, oracle.DIM), dtype=complex)
for p in range(oracle.NMODE):
    for q in range(oracle.NMODE):
        if oracle.H1[p, q] != 0:
            H += oracle.H1[p, q] * (oracle.CD[p] @ oracle.C[q])
N = sum(oracle.CD[p] @ oracle.C[p] for p in range(oracle.NMODE))
ev, _ = np.linalg.eigh(H - oracle.MU * N)
ev = ev - ev.min()
w = np.exp(-oracle.BETA * ev); w /= w.sum()
dE = ev[None, :] - ev[:, None]
dw = w[:, None] - w[None, :]
K = np.where(np.abs(dE) > 1e-10, dw / np.where(dE == 0, 1, dE),
             oracle.BETA * w[:, None])
out["kernel_sample"] = K

# Addition 1: slot-mapped base_chi
out["base_chi_slotmapped"] = oracle._to_solver_slots(out["base_chi"])

# Write to npz and print diagnostics
np.savez_compressed("tests/ed_oracle_reference.npz", **out)

# Print shape, dtype, and size for each array
print("Captured arrays:")
for key in sorted(out.keys()):
    arr = out[key]
    print(f"  {key:24s}: shape={str(arr.shape):15s} dtype={str(arr.dtype):7s} size={arr.nbytes:10d} bytes")

# Print total file size
import os
file_size = os.path.getsize("tests/ed_oracle_reference.npz")
print(f"\nTotal file size: {file_size:,d} bytes ({file_size / 1024 / 1024:.2f} MB)")
print(f"Number of arrays: {len(out)}")
print("wrote tests/ed_oracle_reference.npz:", sorted(out))
