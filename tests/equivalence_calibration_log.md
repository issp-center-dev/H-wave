# RPA/FLEX equivalence table -- calibration log

Schema (`docs/superpowers/plans/2026-08-18-equivalence-table.md`,
Appendix B, as amended by Task 7): one section per calibration event --
date, commit, stage (local/chita/CI-run-N), runner descriptor(s)
(platform, Python, numpy, scipy versions), then a table with the
CANONICAL 7-COLUMN shape (this header is now the authoritative
schema doc; Appendix B's own prose sketch, `cell_id | observable |
residual | decision(atol/re-select/diverges) | reason`, is a 5-column
ABBREVIATION of the same idea and is superseded here -- every event in
this file uses the 7-column shape below, never the 5-column one):

  `cell_id/checkpoint | observable | primary residual | secondary
  residual | derived value | decision | reason`

where "primary residual"/"secondary residual" are event-specific (e.g.
Event 1's benign-fixture/FC-fixture pair; Event 2's local-runner/
chita-runner pair) and "derived value" is whichever single number the
decision turns on (an amplification ratio, a candidate atol, ...) --
named precisely in each event's own table header. Fixture re-selection
attempts each get their own row (candidate params + per-observable
residuals + rejection reason) using the same 7-column shape. Task 6
owns the conditioning-fixture (cell 38 / FC) selection (Event 1, below
-- the fixture-selection step that has to happen before Task 7 can
calibrate anything against cell 38, since the plan requires the
conditioning fixture to actually demonstrate the required amplification
before it is frozen into the registry). Task 7 owns the full
candidate-atol calibration pass (Event 2, below).

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

---

## Event 2 -- Task 7, candidate-atol calibration (local + chita)

- **Date:** 2026-08-18
- **Commit:** b922a1c13b85b2d319bc65ee8c45183dc6ab2a47 (parent of this task's commit -- both measurement runs below were taken against this working tree; see the commit this task ships in for the exact source SHA once committed, per the same convention Event 1 uses).
- **Stage:** local + chita (both advisory pre-freeze proxies; the gating CI runners are Task 9's job)
- **Runner (local):** macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1
- **Runner (chita):** Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1 -- bundle-transferred branch (`git bundle create` -> `scp` -> `git fetch <bundle> feature/equivalence-table:refs/heads/equiv-calib-task7` -> `git checkout` under `~/hwave-bubble-gpu-check`), run with `~/miniconda3/envs/hwave_gpu/bin/python -m tests.equivalence_measure` (CPU only, no CUDA).

**chita transfer finding (STRUCTURAL, fixed before any residual below was trusted):** the `hwave_gpu` conda env has its own installed `hwave` package under `site-packages` (`pip install .`, version 1.1.dev0) that is NOT the bundle-transferred checkout -- `python -c "import hwave; print(hwave.__file__)"` resolved to `~/miniconda3/envs/hwave_gpu/lib/python3.11/site-packages/hwave/__init__.py`, not `~/hwave-bubble-gpu-check/src/hwave`. The first measurement attempt against that stale installed copy produced 2 spurious structural errors (`ValueError not raised` on cell 26, a fragment mismatch on cell 27) that do NOT reproduce locally or against the correct source. Fix: `PYTHONPATH=src` prepended to every chita invocation (verified: `import hwave` then resolves under `~/hwave-bubble-gpu-check/src/hwave`); every number below is from a `PYTHONPATH=src`-correct run.

**chita segfault finding (STRUCTURAL, fixed before the 3-invocation protocol below was run):** with the environment's default BLAS/OMP thread pool, `python -m tests.equivalence_measure` SEGFAULTS on chita partway through a single process -- deterministically after cell 26 (`general.ring.offsite_coulombinter_sameorb.subshape`), while constructing cell 27's (`general.ring.onsite_coulombinter.subshape.mu`) solver. Isolated single-cell and two-cell reproductions of cell 27 alone (or cells 26+27 alone) never crash -- only the full ~27-cell accumulation does, consistent with OpenBLAS/OMP worker-thread accumulation across many repeated small-matrix solver constructions eventually failing `pthread_create` (`free -h` on chita shows 246 GiB available at crash time, ruling out OOM). Fix: `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1` set for every chita invocation (verified: 3/3 clean full-registry runs afterward, 0/3 clean before, including one run with `PYTHONFAULTHANDLER=1` that produced no additional diagnostic -- consistent with a native-thread-pool failure below what `faulthandler`'s SIGSEGV handler can attribute a Python frame to). This finding also shapes `.github/workflows/equivalence-calibration.yml`: the same four env vars are pinned there too, since the gating `ubuntu-latest` runners share the same OpenBLAS-backed numpy/scipy wheels and could plausibly hit the same failure mode at CI's own registry size.

**Protocol:** `python -m tests.equivalence_measure` run 3x on each runner (`find . -name __pycache__ -exec rm -rf {} +` before each invocation); every run: 90 JSON lines, 0 `"error"` records, `source_sha` == `b922a1c13b85b2d319bc65ee8c45183dc6ab2a47` on all 6 invocations. Every quantity below is the runner-local MAX over its 3 invocations (residuals and per-invocation seconds were bitwise identical across the 3 invocations on both runners for every cell measured -- the one-shot IterationMax=1 construction is fully deterministic given a fixed source tree, so MAX-over-3 and any single invocation agree here; the 3-invocation protocol is kept per Appendix B regardless, and Task 9's real gating-runner MIN/MAX will matter once cross-Python-version variation enters).

**Candidate-atol rule applied (Global Constraints, Task-7 local+chita variant):** `raw_max = max(local_max, chita_max)`; `floored = max(raw_max, 1e-15)`; `bound = 10 * floored`, rounded UP to a power of ten; `candidate atol = min(bound, ceiling)`. STOP condition (`raw_max > ceiling`): did not occur -- 0 of 42 comparison-cell/observable pairs. Margin-insufficient condition (`bound > ceiling`, so atol would have to be set EQUAL to ceiling with a recorded coordinator decision): did not occur either -- 0 of 42 pairs strictly exceed; exactly ONE pair lands EXACTLY at its ceiling (`bound == ceiling`, which the validator explicitly allows and is not a margin-insufficient case): cell 38's `chiq` (the conditioning cell -- expected, this is the cell the plan deliberately amplifies).

### Comparison-cell candidate atols (42 rows = 21 Equiv cells x 2 observables; every Equiv cell in the registry)

| cell_id | observable | local residual (max/3) | chita residual (max/3) | candidate atol | decision | reason |
|---|---|---|---|---|---|---|
| `general.ring.onsite_coulombintra.fixedmu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_coulombintra.fixedmu` | chiq | 2.498e-16 | 1.943e-16 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_coulombinter.fixedmu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_coulombinter.fixedmu` | chiq | 1.249e-16 | 1.110e-16 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_hund.fixedmu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_hund.fixedmu` | chiq | 1.110e-16 | 1.110e-16 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_ising.fixedmu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_ising.fixedmu` | chiq | 1.249e-16 | 1.110e-16 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_exchange.fixedmu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_exchange.fixedmu` | chiq | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_pairhop.fixedmu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_pairhop.fixedmu` | chiq | 1.110e-16 | 9.715e-17 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_pairlift.fixedmu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_pairlift.fixedmu` | chiq | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_u_v_hund.mu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `general.ring.onsite_u_v_hund.mu` | chiq | 3.053e-16 | 2.776e-16 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `general.ring.onsite_full_kanamori.mu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `general.ring.onsite_full_kanamori.mu` | chiq | 2.498e-16 | 1.943e-16 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_coulombintra.spinfree.mu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_coulombintra.spinfree.mu` | chiq | 3.886e-16 | 3.053e-16 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_coulombinter.spinfree.mu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_coulombinter.spinfree.mu` | chiq | 1.249e-16 | 1.249e-16 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_hund.spinfree.mu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_hund.spinfree.mu` | chiq | 1.110e-16 | 1.110e-16 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_ising.spinfree.mu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_ising.spinfree.mu` | chiq | 1.249e-16 | 1.249e-16 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_pairlift.spinfree.mu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_pairlift.spinfree.mu` | chiq | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_coulombintra.spindiag.mu` | chi0q | 1.249e-16 | 1.249e-16 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_coulombintra.spindiag.mu` | chiq | 4.164e-16 | 4.719e-16 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_coulombinter.spindiag.mu` | chi0q | 1.249e-16 | 1.249e-16 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_coulombinter.spindiag.mu` | chiq | 1.527e-16 | 1.388e-16 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `general.ring.offsite_coulombinter_sameorb.mu` | chi0q | 4.885e-15 | 4.899e-15 | 1.0e-13 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `general.ring.offsite_coulombinter_sameorb.mu` | chiq | 2.021e-14 | 2.026e-14 | 1.0e-12 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.offsite_coulombinter.mu` | chi0q | 4.885e-15 | 4.899e-15 | 1.0e-13 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.offsite_coulombinter.mu` | chiq | 3.553e-14 | 3.561e-14 | 1.0e-12 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `general.ring.onsite_coulombinter.subshape.mu` | chi0q | 1.665e-16 | 9.715e-17 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `general.ring.onsite_coulombinter.subshape.mu` | chiq | 1.943e-16 | 1.110e-16 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `general.ring.onsite_coulombinter.coefftail.mu` | chi0q | 1.110e-16 | 1.110e-16 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `general.ring.onsite_coulombinter.coefftail.mu` | chiq | 1.388e-16 | 1.388e-16 | 1.0e-14 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `general.ring.onsite_coulombinter.conditioning.mu` | chi0q | 6.432e-13 | 6.432e-13 | 1.0e-11 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `general.ring.onsite_coulombinter.conditioning.mu` | chiq | 3.307e-12 | 3.373e-12 | 1.0e-10 | atol | 10x max(local,chita), floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 (== ceiling) |

### Divergence-diagnostic checkpoints (10 rows = 5 assertions x 2 fixtures; informational -- re-confirms Task 6's Event 1 findings on a second (chita) runner, both `*_diag` ceilings (`mu_diag`= `green_diag`=1e-10) hold with wide margin; no registry atol is gated by this table, only `TestMuGreenDivergenceDiagnostic`'s own hard-coded assertions are)

| cell_id/checkpoint | observable | local residual (max/3) | chita residual (max/3) | derived value | decision | reason |
|---|---|---|---|---|---|---|
| diagnostic assertion 1 | benign | 0.000e+00 | 0.000e+00 | n/a | accept | within `mu_diag` ceiling (1e-10) on both runners |
| diagnostic assertion 1 | fc | 0.000e+00 | 0.000e+00 | n/a | accept | within `mu_diag` ceiling (1e-10) on both runners |
| diagnostic assertion 2 | benign | 0.000e+00 | 0.000e+00 | n/a | accept | within `mu_diag` ceiling (1e-10) on both runners |
| diagnostic assertion 2 | fc | 1.070e-12 | 1.070e-12 | 9.3e+01x under ceiling | accept | within `mu_diag` ceiling (1e-10) on both runners |
| diagnostic assertion 3 | benign | 1.119e-16 | 1.119e-16 | 1.0e+05x under ceiling | accept | within `green_diag` ceiling (1e-10) on both runners |
| diagnostic assertion 3 | fc | 2.289e-16 | 2.289e-16 | 1.0e+05x under ceiling | accept | within `green_diag` ceiling (1e-10) on both runners |
| diagnostic assertion 4 | benign | 1.119e-16 | 1.119e-16 | 1.0e+05x under ceiling | accept | within `green_diag` ceiling (1e-10) on both runners |
| diagnostic assertion 4 | fc | 6.206e-17 | 7.850e-17 | 1.0e+05x under ceiling | accept | within `green_diag` ceiling (1e-10) on both runners |
| diagnostic assertion 5 | benign | 1.110e-16 | 9.714e-17 | 1.0e+05x under ceiling | accept | within `chi0q_mu` ceiling (1e-10) on both runners |
| diagnostic assertion 5 | fc | 1.110e-16 | 5.586e-17 | 1.0e+05x under ceiling | accept | within `chi0q_mu` ceiling (1e-10) on both runners |

**Cross-runner amplification cross-check (assertion 2, the conditioning cell's own seam):** FC/benign = 1.070e-12 / max(0.000e+00, 1e-15) on local; 1.070e-12 / max(0.000e+00, 1e-15) on chita -- both >= 10x, reconfirming Task 6's Event-1 finding on a second, independent platform (Linux x86_64 vs macOS arm64).

### Timing (informational -- the benchmark artifact owns the authoritative table)

- `module_total_seconds` (3 invocations, local): 0.557, 0.545, 0.548
- `module_total_seconds` (3 invocations, chita): 0.632, 0.626, 0.638
- `python -m unittest tests.test_rpa_flex_equivalence_table` process wall time: local 0.85s, chita 0.96s (both far under the 120s freeze-time budget; the 4-gating-runner CI sample lands in Task 9).
