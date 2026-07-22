# APBC validation against ComplexUHF (Docker, non-CI)

Cross-check H-wave UHFk APBC output against mVMC ComplexUHF inside the
Docker sandbox. NOT part of pytest CI.

## Workflow

```
cd tests/validation/apbc_complexuhf
chmod +x build_complexuhf.sh run.sh
./build_complexuhf.sh
./run.sh case_1d_hubbard
```

`build_complexuhf.sh` copies the read-only mVMC-1.4.0/ tree into
build/mvmc/ and runs `cmake` + `make UHF vmcdry.out`. mVMC source is
NEVER modified.

`run.sh` runs both solvers on the same physical case and invokes
`compare.py` to diff total energy and per-site one-body Green
function.

## Case `case_1d_hubbard`

1D Hubbard chain, L = 8, t = 1, U = 2, **quarter filling** (Ncond = 4).
APBC in x. Quarter filling avoids L = 8 half-filling degeneracy
(spec section 7).

H-wave inputs (parsed by `qlmsio.read_input_k.QLMSkInput`):

* `input.toml` — driver config, declares `BoundaryCondition =
  ["antiperiodic", "periodic", "periodic"]`.
* `[file.input.interaction]`: `Geometry.dat` (read by `wan90.read_geom`),
  `Transfer.dat`, `CoulombIntra.dat` (wannier90-like format).
* `[file.input]`: `geometry_uhf.dat` (read by `wan90.read_geometry`,
  different format — supplies `site2vec` used to attach the inverse APBC
  gauge phase to the output Green function) and `OneBodyG.dat` (the
  `(i, s, j, t)` index list to print).

ComplexUHF inputs:

* `stan.in` — StdFace input. APBC along the chain is requested via
  `phase0 = 180.0` (degrees). StdFace converts this to
  `ExpPhase = exp(i pi) = -1` and ultimately writes a `-1` factor on the
  wrap-around hopping in `trans.def`.
* `complexuhf_modpara_override.txt` — also flips `NMPTrans` to `-1` as a
  documentary marker. See the caveat below: this alone does **not**
  enable APBC in the standalone `UHF` binary.

### Caveat: how APBC is actually enabled in ComplexUHF

mVMC documents APBC as "set NMPTrans negative". In `UHF` (the
ComplexUHF binary), that bit is read into `APFlag`
(`mVMC-1.4.0/src/ComplexUHF/readdef.c:382-386`), but `APFlag` is never
consumed anywhere else in the ComplexUHF source tree — the
Hamiltonian build path `xsetmem`, `readdef`, `makeham`, `cal_energy`
contain no other `APFlag` reference. Only the VMC binary
(`src/mVMC/`) actually uses `APFlag` to flip the sign of QPTrans
entries on read.

For UHF, the only working route is to bake APBC into `Trans.def`
directly: a `-1` factor on each wrap-around bond. StdFace produces
that when `phase0 = 180.0` (deg) is set in `stan.in`, because
`StdFace_ModelUtil.c::StdFace_Hopping` multiplies each bond's `t` by
`Cphase`. We therefore drive APBC through `stan.in`, not through the
NMPTrans override.

## Tolerances (spec 7.2)

* energy: rel diff < 1e-6
* per-site Green: max |delta| < 1e-6

## Validation results

case_1d_hubbard (2026-06-26, mVMC 1.4.0, H-wave at branch
feature/uhfk-antiperiodic-bc):

```
[case_1d_hubbard] energy:   hwave=-6.3910362638e+00  complexuhf=-6.3910362601e+00  rel=5.81e-10  OK
[case_1d_hubbard] greenone: max |delta| = 1.187e-09 over 16 (i,s,j,t) pairs    OK
```

Both observables agree well within the 1e-6 tolerance. The H-wave
energy `-6.39103626...` also matches the closed-form HF total for
quarter-filled L=8 APBC Hubbard with U=2:
`E = -8 cos(pi/8) + U L (n/2)^2 = -7.391036... + 1.0 = -6.391036...`.

## Case `case_1d_hubbard_sublattice`

Same physical problem as `case_1d_hubbard` (L=8 1D Hubbard, t=1, U=2,
quarter filling, APBC in x), but the H-wave side uses **SubShape =
[2, 1, 1]**. Confirms that APBC under sublattice fold reproduces the
same physical Green function as without sublattice fold, end-to-end
against ComplexUHF.

```
[case_1d_hubbard_sublattice] energy:   hwave=-6.3910362638e+00  complexuhf=-6.3910362601e+00  rel=5.81e-10  OK
[case_1d_hubbard_sublattice] greenone: max |delta| = 1.187e-09 over 16 (i,s,j,t) pairs    OK
```

Energy and per-site Green agree with both ComplexUHF and the
SubShape=[1,1,1] run to ~1e-9. This validates the v2 sublattice APBC
implementation end-to-end (spectrum, fold, gauge unwinding) when
combined with the green deflate fix on develop (PR #35).

## Runtime workaround

`qlmsio.wan90.read_geometry` uses the deprecated `np.float` alias
(removed in NumPy >= 1.20). `run.sh` injects a `np.float = float`
shim before importing `hwave.qlms`; H-wave source is not modified.
