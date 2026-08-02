# mVMC PairProduct E2E validation (Docker, non-CI)

End-to-end check of the H-wave UHFk -> bridge -> mVMC PairProduct
pipeline. Runs inside the Docker sandbox; NOT part of pytest CI.

## Workflow

```bash
# 1. Build mVMC binaries (UHF + vmcdry + vmc.out) once
chmod +x ../apbc_complexuhf/build_complexuhf.sh
../apbc_complexuhf/build_complexuhf.sh

# 2. Run a case
chmod +x ./run.sh
./run.sh case_pbc
./run.sh case_apbc

# 2'. Override the rank-lift noise amplitude (default 1e-8):
EPSILON_NOISE=1e-7 ./run.sh case_apbc
# Or pick a different RNG seed (default 7919):
RNG_SEED=42 ./run.sh case_apbc
```

`run.sh case_X` orchestrates:

1. H-wave UHFk SCF on `case_X/input.toml` -> `output/{eigen,occupation,greenone,energy}.{npz,dat}`.
2. `vmcdry.out` on `case_X/stan.in` -> mVMC `.def` files (orbitalidx, trans, namelist, modpara, ...).
3. `tools/uhfk_to_mvmc.py` -> `zqp_orbital_uhfk.dat` with `--check-density`
   against `greenone.dat` (tol 1e-10).
4. Append `InOrbital ...zqp_orbital_uhfk.dat` to `namelist.def`.
5. `vmc.out namelist.def` with `NVMCCalMode=1, NSROptItrStep=0`
   (initial-WF measurement, no optimization).
6. `compare.py` reads H-wave `Energy_Total` and mVMC `<H>` from
   `zvo_out_001.dat` and reports `delta` plus pass/fail.

## Cases

### `case_pbc` — PBC 1D Hubbard L=8, t=1, U=2, Ncond=2 (1 e per spin)

Drops to Ne=2 because L=8 quarter filling (Ne=4) has a Fermi-level
degeneracy at k=+/-pi/4 that the bridge correctly rejects. With Ne=2
only k=0 is occupied; the (k, -k) pair construction is degenerate with
the naive same-k pairing, so the PairProduct -> mVMC pipeline is the
cleanest possible test case.

```
[case_pbc] H-wave UHFk total = -3.7500000000e+00
[case_pbc] mVMC initial mean = -3.7480000000e+00  (n_bins=1; stderr undefined)
[case_pbc] delta = +2.0000e-03  (0.05% rel)
[case_pbc] OK  (within 1% rel; single-bin VMC sample)
```

### `case_apbc` — APBC 1D Hubbard L=8, t=1, U=2, Ncond=4 (closed shell, 1 (k, -k) pair shell)

APBC k_phys = (2n+1)pi/L, so the two lowest per spin are k=+/-pi/8
(paired, non-degenerate). Ne=4 (per spin = 2) fills the inner shell
exactly (HOMO-LUMO gap 1.083, closed shell).

**Status: FIXED via rank-lift noise in the bridge.**

```
[case_apbc] H-wave UHFk total = -6.3910362601e+00
[case_apbc] mVMC initial mean = -6.3700363544e+00  (n_bins=1; stderr undefined)
[case_apbc] delta = +2.1000e-02  (0.33% rel)
[case_apbc] OK  (within 1% rel; single-bin VMC sample)
```

Before the rank-lift noise was added the bridge's F (built from a
single (k, -k) shell, rank min(Ne_up, Ne_down) = 2 in an 8x8 matrix)
produced ``<H>`` = -5.396 instead of -6.391: mVMC's Pfaffian Slater
evaluation went ill-conditioned and sampled a non-translation-
invariant state (``<n_{0,up}>`` = 5e-4 vs expected 0.25). Sweeping
``NSPGaussLeg``, swapping (k, -k) F for naive same-k F, switching to
tilde basis, and disabling Gutzwiller / Jastrow / TransSym had no
effect; the bridge's own Pfaffian-mimic check
(``tests/test_uhfk_to_mvmc_pfaffian_density.py``) confirmed F itself
was mathematically correct.

The root cause is that mVMC's Pfaffian Slater algorithm becomes
singular when the F matrix is exactly rank-deficient. The bridge now
mirrors ComplexUHF's own workaround (``mVMC-1.4.0/src/ComplexUHF/
output.c:274``), adding a tiny uniform noise to each averaged
``params[idx]`` before writing. Noise amplitude defaults to ``1e-8``
(``--epsilon-noise 0`` disables it). For ``ComplexType 1``
(``orbitalidx.def``) the noise is applied to both real and imag parts;
for ``ComplexType 0`` only the real part is perturbed. The default
amplitude is small enough that ``<H>`` recovers UHF energy to <1% rel
without significantly perturbing the per-site density.

### `case_apbc_halffill` — APBC 1D Hubbard L=8, t=1, U=2, Ncond=8 (half filling, 2 (k, -k) shells)

```
[case_apbc_halffill] H-wave UHFk total = -6.4525037190e+00
[case_apbc_halffill] mVMC initial mean = -6.4775037190e+00  (n_bins=1)
[case_apbc_halffill] delta = -2.5000e-02  (0.39% rel)
[case_apbc_halffill] OK  (within 1% rel; single-bin VMC sample)
```

Half-filling adds the outer (k=+/-3pi/8) shell, lifting F to rank 4.
mVMC's Pfaffian then reconstructs the correct Slater determinant and
the energy agrees with H-wave UHF to <1% relative error. This is the
recommended APBC reference case for the bridge.

## Implementation notes

- `case_*/stan.in` omits `NSROptItrSmp` (vmcdry hangs on the
  "specified but not used" warning when `NSROptItrStep = 0`).
- `case_*/input.toml` keeps `flag_fock = false` in `[mode]` (NOT
  `[mode.param]`): `uhfk.py:82` reads it from `info_mode` directly,
  and Fock on collapses CoulombIntra into a single spin-mixed block
  that the bridge correctly rejects with `column_spin = -1`.
- `run.sh` post-edits `orbitalidx.def` to `ComplexType 1` because
  StdFace defaults to `ComplexType 0` for real-coefficient Hubbard, but
  the bridge always writes complex F values.
- Builds reuse the existing apbc_complexuhf mVMC build at
  `tests/validation/apbc_complexuhf/build/mvmc/build/src/mVMC/`.
