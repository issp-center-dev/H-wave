# case_pbc_sz2 — 1D PBC L=8 Hubbard, 2Sz=2 mVMC E2E fixture

Sz-fixed 2Sz=2 (Ne_up=3, Ne_down=1) — the v3 A ケース for the
InOrbitalGeneral pair form. StdFace side generates orbitalidx_general.def
(6-column) via `OrbitalGeneral 1` in `stan.in`.

## Status (v3)

Runnable end-to-end (harness-gated). Precondition: pytest
`tests/test_uhfk_to_mvmc_case_pbc_sz2_occupation.py` and the
`run.sh case_pbc_sz2` occupation gate both call
`scripts/assert_occupation.py case_pbc_sz2 <work>` and must pass.

Target occupation (single source of truth: `_CASE_TARGETS["case_pbc_sz2"]`
in `scripts/assert_occupation.py`):

- `down`: k=(0, 0, 0)
- `up`:   k=(0, 0, 0), (+1, 0, 0)=+π/4, (-1, 0, 0)=-π/4

Same-spin up-up excess consumed by the non-self canonical block
(+π/4, -π/4); self-pair k=0 covers the cross pair; self-pair k=π stays
empty. Section 3.2 pair-closure conditions hold on this table.

## Known limitation (v3 in progress)

Full `run.sh case_pbc_sz2` currently stops at the bridge invocation with
`N_up = 3 != N_down = 1; v1 spec section 7 requires AntiParallel
Sz-fixed sector`. The pytest harness-gate (Step 1 + Step 1.5) is green:

- Step 1 (H-wave UHFk SCF): converges to the target FM state.
- Step 1.5 (`scripts/assert_occupation.py case_pbc_sz2 <work>`): passes.
- Step 2 (`vmcdry.out`): writes `namelist.def`, `orbitalidxgen.def`, ...
- Step 3 (`tools/uhfk_to_mvmc.py --check-density`): fails fast because
  the v1 parallel-Slater dispatch cannot ingest `orbitalidxgen.def`
  (6-column general format). The General dispatch is wired in Task 11
  and the full E2E gate follows in Task 12.

`stan.in` also required `NSPGaussLeg` to be omitted: with `2Sz != 0`
StdFace projects on the fixed-Sz subspace by construction, so vmcdry
warns "NSPGaussLeg is SPECIFIED but will NOT be USED" and halts before
writing `modpara.def` / `namelist.def`. StdFace has no `OrbitalGeneral`
keyword; setting `2Sz = 2` alone triggers the `orbitalidxgen.def` write
and the commented `# OrbitalGeneral orbitalidxgen.def` line in
`namelist.def` (see `StdFace_main.c:1426`).
