# RPA/FLEX equivalence table -- calibration log

Schema (`docs/superpowers/plans/2026-08-18-equivalence-table.md`, Appendix B):
one section per calibration event -- date, commit, stage
(local/chita/CI-run-N), runner descriptor (platform, Python, numpy,
scipy versions), then a table `cell_id | observable | residual |
decision(atol/re-select/diverges) | reason`. Fixture re-selection
attempts each get their own row (candidate params + per-observable
residuals + rejection reason). Task 7 owns the full candidate-atol
calibration pass; this file's first entries (below) are Task 6's own
conditioning-fixture (cell 38 / FC) attempts -- the fixture-selection
step that has to happen before Task 7 can calibrate anything against
cell 38, since the plan requires the conditioning fixture to actually
demonstrate the required amplification before it is frozen into the
registry.

---

## Event 1 -- Task 6, cell 38 (FC) conditioning-fixture selection

- **Date:** 2026-08-18
- **Commit:** a93bf5bdea0359295353f97d93953735cbda80d0 (parent of this
  task's commit -- the working tree these measurements were taken
  against; see the commit this file ships in for the exact source SHA
  once committed)
- **Stage:** local
- **Runner:** macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13,
  numpy 1.26.4, scipy 1.17.1

**Objective (plan Task 6 Step 2):** FC = E1's files (`tests/
equivalence_input/orb1`, off-site same-orbital `coulombinter.dat`,
V=1.0) with `filling=0.5`, `CellShape=(4,4,1)`, `Nmat=256`, starting at
`T=0.2` (Appendix A's stated candidate). The mu/Green divergence
diagnostic (`TestMuGreenDivergenceDiagnostic`'s apparatus,
`_diagnostic_residuals` in `tests/test_rpa_flex_equivalence_table.py`),
run independently on this fixture and on the benign cell-8 fixture
(`general.ring.onsite_u_v_hund.mu`, E2, T=2.0), must show >= 10x
amplification on at least one of (assertion 2 the mu seam, assertion 4
the composed seam, assertion 5 downstream chi0q), with
`max(residual, 1e-15)` denominators. If not met, halve T (up to twice)
and re-measure, logging every attempt.

**van-Hove-shoulder verification:** on FC's actual `CellShape=(4,4,1)`
16-k-point mesh, H0's eigenvalues are
`[-3, -2,-2,-2,-2, -1,-1,-1,-1, 1,1, 2,2,2,2, 5]`. `filling=0.5` targets
8-of-16 one-spin states, which lands the Fermi level EXACTLY inside the
4-fold-degenerate `eps=-1.0` cluster -- a coarse-mesh analogue of a van
Hove shoulder (confirmed by direct `H0_eigenvalue` inspection; a finer
64x64 k-mesh on the identical dispersion shows the true DOS peak near
`eps~-0.93`, consistent with this cluster, and the `T=0.2`/`filling=0.5`
`_find_mu` root on that fine mesh lands at `mu~-0.51`, inside the
shoulder). No re-selection of `filling` was needed -- 0.5 is kept as
Appendix A's stated value.

### Attempt 1 -- T=0.2 (Appendix A's stated starting value)

| cell_id / checkpoint | observable | benign residual (cell 8, T=2.0) | FC residual (T=0.2) | amplification | decision | reason |
|---|---|---|---|---|---|---|
| diagnostic assertion 1 (shared-mu identity control) | mu (scalar) | 0.0 | 0.0 | n/a (control) | accept | both sides bit-identical, as the control predicts -- FLEX inherits `_find_mu` unchanged |
| diagnostic assertion 2 (the mu seam) | mu (scalar) | 0.0 | 1.070254995738651e-12 | >= 1.07e-12 / max(0.0, 1e-15) = **1070x** | **accept** | clears the >=10x bar by ~2-3 orders of magnitude on the FIRST candidate; no halving needed |
| diagnostic assertion 3 (Green-builder seam, bare) | Green (max\|diff\|) | 1.1190705717918148e-16 | 2.2887833992611187e-16 | ~2.0x | accept (not required to amplify -- informational) | not one of the plan's three amplification-triple members |
| diagnostic assertion 4 (composed seam, dressed-zero-Sigma mu) | Green (max\|diff\|) | 1.1190705717918148e-16 | 6.206335383118183e-17 | ~0.55x (no amplification) | accept (criterion already met by assertion 2) | this seam does not amplify on FC at T=0.2; not needed since assertion 2 alone satisfies the >=10x criterion |
| diagnostic assertion 5 (downstream chi0q) | chi0q (max\|diff\|) | 1.110233872252785e-16 | 1.1102640084031525e-16 | ~1.0x (no amplification) | accept (criterion already met by assertion 2) | likewise not needed |
| `general.ring.onsite_coulombinter.conditioning.mu` (cell 38 itself, ORDINARY Equiv row via the full one-shot `solve()` pipeline -- a DIFFERENT computation than the isolated-Sigma=0 diagnostic above) | chi0q | -- | 0.0 (bit-identical) | -- | atol = ceiling (`chi0q_mu`=1e-10), provisional per Task 5's convention | full-pipeline chi0q matches to round-off despite the seam-level amplification found in isolation -- FLEX's post-mix mu differs from the diagnostic's isolated Sigma=0 mu (flex.py:741-743 re-solves mu against the iteration's own self-energy before finishing), so this is not a contradiction |
| `general.ring.onsite_coulombinter.conditioning.mu` (cell 38, chiq) | chiq | -- | 3.3066883022456364e-12 | -- | atol = ceiling (`chiq_mu`=1e-10), provisional per Task 5's convention | comfortably inside ceiling |

**Decision:** ACCEPT T=0.2, Nmat=256, filling=0.5 as FC (Appendix A's
stated starting candidate, unchanged). The amplification criterion is
satisfied by assertion 2 alone (~1070x), so the halve-T re-choice
procedure was not exercised -- this is the only logged attempt.

**Monotonicity note (not a further attempt, recorded for completeness):**
a scan of T=0.2/0.1/0.05 on the SAME diagnostic apparatus (informational
only, not part of the accept/reject decision since T=0.2 already
clears the bar) shows assertion 2's residual growing monotonically as T
shrinks: `1.070254995738651e-12` (T=0.2) ->
`1.5136780717739384e-12` (T=0.1) -> `1.6872059305228504e-12` (T=0.05),
consistent with the derivation in `TestConditioningAmplification`'s
docstring (dN/dmu grows as T shrinks near the near-degenerate eps=-1.0
cluster, so the two independent mu-root-finders' numerical fixed points
drift further apart). Assertions 3-5's Green/chi0q residuals over the
same scan stay within a factor of a few of their benign-fixture values
and do not show the same clean monotonic growth at this coarse mesh --
consistent with assertion 2 (the mu seam) being the dominant, and only
strictly necessary, amplified channel for this fixture.

**Timing (informational, not a runtime concern):** the full `solve()`
pair for cell 38 (RPA + FLEX, Nmat=256) takes ~0.16s locally
(rpa ~0.15s, flex ~0.01s) -- despite Nmat=256, this fixture is
1-orbital/16-k-point (nd=1), so it is not a CI-budget risk. The
isolated diagnostic apparatus (no `solve()`, Sigma=0 only) is faster
still, a few milliseconds per fixture.
