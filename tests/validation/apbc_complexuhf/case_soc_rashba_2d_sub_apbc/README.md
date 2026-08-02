# case_soc_rashba_2d_sub_apbc (v3.6 apbc_complexuhf fixture)

ComplexUHF cross-check fixture for the v3.6 SOC + APBC + SubShape
combination. Mirrors the H-wave shipping fixture at
`tests/validation/uhfk_mvmc_pairproduct/case_soc_rashba_2d_sub_apbc/`
with identical Hamiltonian (Rashba SOC on 6x4 square lattice, APBC in
x, folded BZ [3, 2, 1], Ncond = 8). ComplexUHF (mVMC-1.4.0's UHF
binary) is run under this fixture to produce a bit-independent
`zvo_UHF_cisajs.dat` for the Phase 6 G2a-in-memory-A gate.

Ships in v3.6 (spec `docs/superpowers/specs/2026-07-09-uhfk-mvmc-
pairproduct-general-v36-design.md`, plan `docs/superpowers/plans/
2026-07-09-uhfk-mvmc-pairproduct-general-v36.md` Task 5).

## Geometry / interaction

Reused verbatim from the H-wave fixture (bit-identical Hamiltonian
matches force the ComplexUHF cross-check to compare equal solvers, not
different discretizations):

- `CellShape` = `[6, 4, 1]`
- `SubShape` = `[2, 2, 1]` (folded BZ `L_folded = [3, 2, 1]`)
- `BoundaryCondition` = `["antiperiodic", "periodic", "periodic"]`
- Rashba SOC: `t = 1`, `alpha = 0.5`, `U = 2`
- Total physical sites: 24 (spin-orbital dim 48)
- `Ncond = 8` (closed-shell selection at theta = (pi, 0, 0) via
  Phase 2 step 2b's `soc_pilot_gap.py`)

## APBC realization

- **H-wave side**: `input.toml` `BoundaryCondition = ["antiperiodic",
  "periodic", "periodic"]` (theta = (pi, 0, 0)) — the same APBC gauge
  the shipping bridge composes on top of the sub_offset phase.
- **ComplexUHF side**: `stan.in` `phase0 = 180.0` -> StdFace bakes the
  (-1) wrap sign into `Trans.def`; ComplexUHF's UHF binary reads the
  signed Trans and the resulting Hamiltonian matches H-wave's exactly.
- **NMPTrans**: `complexuhf_modpara_override.txt` records the legacy
  `NMPTrans = -1` convention; ComplexUHF reads this into the `APFlag`
  bit (mVMC-1.4.0 `src/ComplexUHF/readdef.c:382-386`) but never
  consumes it (no other reference to `APFlag` in the ComplexUHF source
  tree), so the APBC signal actually lives in phase0/Trans.def and the
  NMPTrans override is a no-op-equivalent record.

## v3.6 gate role (Phase 6 preview)

The Phase 6 `run.sh case_soc_rashba_2d_sub_apbc` will produce
`${WORK_DIR}/complexuhf/zvo_UHF_cisajs.dat` from this fixture. The v3.6
seven-gate contract consumes it in:

- **G2a-in-memory-A**: `compare_against_onebodyg_uhf_general(A_ship,
  complexuhf/zvo_UHF_cisajs.dat)` at 1e-6 tolerance. Cross-checks the
  shipping A density directly against a bit-independent SOC-APBC-SCF
  ground truth.
- **G2b**: `compare_against_green_sublattice(gauge_lift-lifted
  green_sublattice, complexuhf/zvo_UHF_cisajs.dat)` at 1e-6 tolerance.
  Cross-checks the density-check's gauge_lift composition against the
  same ComplexUHF ground truth.

Phase 5 is HARNESS-ONLY: the fixture files are committed and
`compare.py` gains SOC 4-index parsing, but no live ComplexUHF SCF is
run yet. Phase 6 first-runs G2a/G2b through `run.sh`.
