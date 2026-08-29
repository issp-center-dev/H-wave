# RPA/FLEX equivalence table -- calibration log

Schema: one section per calibration event -- date, commit, stage
(which development machine, or which CI run), runner descriptor(s)
(platform, Python, numpy, scipy versions), then a table with the
CANONICAL 7-COLUMN shape (this header is the authoritative schema for
that table; every event in this file uses the 7-column shape below):

  `cell_id/checkpoint | observable | primary residual | secondary
  residual | derived value | decision | reason`

where "primary residual"/"secondary residual" are event-specific (e.g.
Event 1's benign-fixture/FC-fixture pair; Event 2's macOS/Linux runner
pair) and "derived value" is whichever single number the
decision turns on (an amplification ratio, a candidate atol, ...) --
named precisely in each event's own table header. Fixture re-selection
attempts each get their own row (candidate params + per-observable
residuals + rejection reason) using the same 7-column shape. Event 1
below records the conditioning-fixture (cell 38 / FC) selection, which
has to precede any calibration against cell 38: that fixture must be
shown to demonstrate the required amplification before it is frozen
into the registry. Event 2 records the full candidate-atol calibration
pass. Event 3 records the freeze calibration on the four gating
continuous-integration runners, whose primary/secondary residual
columns are the MIN and MAX over those runners.

---

## Event 1 -- cell 38 (FC) conditioning-fixture selection

- **Date:** 2026-08-18
- **Commit:** a93bf5bdea0359295353f97d93953735cbda80d0 (the working
  tree these measurements were taken against; see the commit this file
  ships in for the exact source SHA once committed)
- **Stage:** development machine
- **Runner:** macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13,
  numpy 1.26.4, scipy 1.17.1

**Objective:** FC = the one-orbital fixture files (`tests/
equivalence_input/orb1`, off-site same-orbital `coulombinter.dat`,
V=1.0) with `filling=0.5`, `CellShape=(4,4,1)`, `Nmat=256`, starting at
`T=0.2` (the first candidate tried). The mu/Green divergence
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
the stated starting value.

### Attempt 1 -- T=0.2 (the stated starting value)

| cell_id / checkpoint | observable | benign residual (cell 8, T=2.0) | FC residual (T=0.2) | amplification | decision | reason |
|---|---|---|---|---|---|---|
| diagnostic assertion 1 (shared-mu identity control) | mu (scalar) | 0.0 | 0.0 | n/a (control) | accept | both sides bit-identical, as the control predicts -- FLEX inherits `_find_mu` unchanged |
| diagnostic assertion 2 (the mu seam) | mu (scalar) | 0.0 | 1.070254995738651e-12 | >= 1.07e-12 / max(0.0, 1e-15) = **1070x** | **accept** | clears the >=10x bar by ~2-3 orders of magnitude on the FIRST candidate; no halving needed |
| diagnostic assertion 3 (Green-builder seam, bare) | Green (max\|diff\|) | 1.1190705717918148e-16 | 2.2887833992611187e-16 | ~2.0x | accept (not required to amplify -- informational) | not one of the three members of the amplification triple |
| diagnostic assertion 4 (composed seam, dressed-zero-Sigma mu) | Green (max\|diff\|) | 1.1190705717918148e-16 | 6.206335383118183e-17 | ~0.55x (no amplification) | accept (criterion already met by assertion 2) | this seam does not amplify on FC at T=0.2; not needed since assertion 2 alone satisfies the >=10x criterion |
| diagnostic assertion 5 (downstream chi0q) | chi0q (max\|diff\|) | 1.110233872252785e-16 | 1.1102640084031525e-16 | ~1.0x (no amplification) | accept (criterion already met by assertion 2) | likewise not needed |
| `general.ring.offsite_coulombinter.conditioning.mu` (cell 38 itself, ORDINARY Equiv row via the full one-shot `solve()` pipeline -- a DIFFERENT computation than the isolated-Sigma=0 diagnostic above) | chi0q | -- | 0.0 (bit-identical) | -- | atol = ceiling (`chi0q_mu`=1e-10), provisional per the registry's convention | full-pipeline chi0q matches to round-off despite the seam-level amplification found in isolation -- FLEX's post-mix mu differs from the diagnostic's isolated Sigma=0 mu (flex.py:741-743 re-solves mu against the iteration's own self-energy before finishing), so this is not a contradiction |
| `general.ring.offsite_coulombinter.conditioning.mu` (cell 38, chiq) | chiq | -- | 3.3066883022456364e-12 | -- | atol = ceiling (`chiq_mu`=1e-10), provisional per the registry's convention | comfortably inside ceiling |

**Decision:** ACCEPT T=0.2, Nmat=256, filling=0.5 as FC (the stated
starting candidate, unchanged). The amplification criterion is
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

## Event 2 -- candidate-atol calibration on two development machines

- **Date:** 2026-08-18
- **Commit:** b922a1c13b85b2d319bc65ee8c45183dc6ab2a47 (both measurement runs below were taken against this working tree; see the commit this file ships in for the exact source SHA once committed, per the same convention Event 1 uses).
- **Stage:** two development machines (both advisory pre-freeze proxies; the gating CI runners are measured separately)
- **Runner (macOS arm64):** macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1
- **Runner (Linux x86_64):** Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1 -- the branch was transferred to that machine as a git bundle and checked out there, then run with that machine's own conda Python via `python -m tests.equivalence_measure` (CPU only, no CUDA).

**Transfer finding on the Linux x86_64 machine (STRUCTURAL, fixed before any residual below was trusted):** the conda environment used there has its own installed `hwave` package under `site-packages` (`pip install .`, version 1.1.dev0) that is NOT the transferred checkout -- `python -c "import hwave; print(hwave.__file__)"` resolved to that environment's `site-packages/hwave/__init__.py`, not the checkout's `src/hwave`. The first measurement attempt against that stale installed copy produced 2 spurious structural errors (`ValueError not raised` on cell 26, a fragment mismatch on cell 27) that do NOT reproduce on the other machine or against the correct source. Fix: `PYTHONPATH=src` prepended to every invocation on that machine (verified: `import hwave` then resolves under the checkout's `src/hwave`); every number below is from a `PYTHONPATH=src`-correct run.

**Segfault finding on the Linux x86_64 machine (STRUCTURAL, fixed before the 3-invocation protocol below was run):** with the environment's default BLAS/OMP thread pool, `python -m tests.equivalence_measure` SEGFAULTS on that machine partway through a single process -- deterministically after cell 26 (`general.ring.offsite_coulombinter_sameorb.subshape`), while constructing cell 27's (`general.ring.onsite_coulombinter.subshape.mu`) solver. Isolated single-cell and two-cell reproductions of cell 27 alone (or cells 26+27 alone) never crash -- only the full ~27-cell accumulation does, consistent with OpenBLAS/OMP worker-thread accumulation across many repeated small-matrix solver constructions eventually failing `pthread_create` (`free -h` there shows 246 GiB available at crash time, ruling out OOM). Fix: `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1` set for every invocation on that machine (verified: 3/3 clean full-registry runs afterward, 0/3 clean before, including one run with `PYTHONFAULTHANDLER=1` that produced no additional diagnostic -- consistent with a native-thread-pool failure below what `faulthandler`'s SIGSEGV handler can attribute a Python frame to). This finding also shapes `.github/workflows/equivalence-calibration.yml`: the same four env vars are pinned there too, since the gating `ubuntu-latest` runners share the same OpenBLAS-backed numpy/scipy wheels and could plausibly hit the same failure mode at CI's own registry size.

**Protocol:** `python -m tests.equivalence_measure` run 3x on each runner (`find . -name __pycache__ -exec rm -rf {} +` before each invocation); every run: 90 JSON lines, 0 `"error"` records, `source_sha` == `b922a1c13b85b2d319bc65ee8c45183dc6ab2a47` on all 6 invocations. Every quantity below is the runner-local MAX over its 3 invocations (residuals and per-invocation seconds were bitwise identical across the 3 invocations on both runners for every cell measured -- the one-shot IterationMax=1 construction is fully deterministic given a fixed source tree, so MAX-over-3 and any single invocation agree here; the 3-invocation protocol is kept regardless, and the real gating-runner MIN/MAX will matter once cross-Python-version variation enters).

**Candidate-atol rule applied (the two-development-machine variant):** `raw_max = max(macos_max, linux_max)`; `floored = max(raw_max, 1e-15)`; `bound = 10 * floored`, rounded UP to a power of ten; `candidate atol = min(bound, ceiling)`. STOP condition (`raw_max > ceiling`): did not occur -- 0 of 42 comparison-cell/observable pairs. Margin-insufficient condition (`bound > ceiling`, so atol would have to be set EQUAL to ceiling with a recorded maintainer decision): did not occur either -- 0 of 42 pairs strictly exceed; exactly ONE pair lands EXACTLY at its ceiling (`bound == ceiling`, which the validator explicitly allows and is not a margin-insufficient case): cell 38's `chiq` (the conditioning cell -- expected, this is the cell deliberately chosen to amplify).

### Comparison-cell candidate atols (42 rows = 21 Equiv cells x 2 observables; every Equiv cell in the registry)

| cell_id | observable | macOS arm64 residual (max/3) | Linux x86_64 residual (max/3) | candidate atol | decision | reason |
|---|---|---|---|---|---|---|
| `general.ring.onsite_coulombintra.fixedmu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_coulombintra.fixedmu` | chiq | 2.498e-16 | 1.943e-16 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_coulombinter.fixedmu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_coulombinter.fixedmu` | chiq | 1.249e-16 | 1.110e-16 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_hund.fixedmu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_hund.fixedmu` | chiq | 1.110e-16 | 1.110e-16 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_ising.fixedmu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_ising.fixedmu` | chiq | 1.249e-16 | 1.110e-16 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_exchange.fixedmu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_exchange.fixedmu` | chiq | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_pairhop.fixedmu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_pairhop.fixedmu` | chiq | 1.110e-16 | 9.715e-17 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_pairlift.fixedmu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_pairlift.fixedmu` | chiq | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-12 |
| `general.ring.onsite_u_v_hund.mu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `general.ring.onsite_u_v_hund.mu` | chiq | 3.053e-16 | 2.776e-16 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `general.ring.onsite_full_kanamori.mu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `general.ring.onsite_full_kanamori.mu` | chiq | 2.498e-16 | 1.943e-16 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_coulombintra.spinfree.mu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_coulombintra.spinfree.mu` | chiq | 3.886e-16 | 3.053e-16 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_coulombinter.spinfree.mu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_coulombinter.spinfree.mu` | chiq | 1.249e-16 | 1.249e-16 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_hund.spinfree.mu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_hund.spinfree.mu` | chiq | 1.110e-16 | 1.110e-16 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_ising.spinfree.mu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_ising.spinfree.mu` | chiq | 1.249e-16 | 1.249e-16 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_pairlift.spinfree.mu` | chi0q | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_pairlift.spinfree.mu` | chiq | 1.110e-16 | 9.714e-17 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_coulombintra.spindiag.mu` | chi0q | 1.249e-16 | 1.249e-16 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_coulombintra.spindiag.mu` | chiq | 4.164e-16 | 4.719e-16 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_coulombinter.spindiag.mu` | chi0q | 1.249e-16 | 1.249e-16 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.onsite_coulombinter.spindiag.mu` | chiq | 1.527e-16 | 1.388e-16 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `general.ring.offsite_coulombinter_sameorb.mu` | chi0q | 4.885e-15 | 4.899e-15 | 1.0e-13 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `general.ring.offsite_coulombinter_sameorb.mu` | chiq | 2.021e-14 | 2.026e-14 | 1.0e-12 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.offsite_coulombinter.mu` | chi0q | 4.885e-15 | 4.899e-15 | 1.0e-13 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `reduced.ring.offsite_coulombinter.mu` | chiq | 3.553e-14 | 3.561e-14 | 1.0e-12 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `general.ring.onsite_coulombinter.subshape.mu` | chi0q | 1.665e-16 | 9.715e-17 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `general.ring.onsite_coulombinter.subshape.mu` | chiq | 1.943e-16 | 1.110e-16 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `general.ring.onsite_coulombinter.coefftail.mu` | chi0q | 1.110e-16 | 1.110e-16 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `general.ring.onsite_coulombinter.coefftail.mu` | chiq | 1.388e-16 | 1.388e-16 | 1.0e-14 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `general.ring.offsite_coulombinter.conditioning.mu` | chi0q | 6.432e-13 | 6.432e-13 | 1.0e-11 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 |
| `general.ring.offsite_coulombinter.conditioning.mu` | chiq | 3.307e-12 | 3.373e-12 | 1.0e-10 | atol | 10x max over both runners, floored 1e-15, rounded up to a power of ten, capped at ceiling 1e-10 (== ceiling) |

### Divergence-diagnostic checkpoints (10 rows = 5 assertions x 2 fixtures; informational -- re-confirms Event 1's findings on a second (Linux x86_64) runner, both `*_diag` ceilings (`mu_diag`= `green_diag`=1e-10) hold with wide margin; no registry atol is gated by this table, only `TestMuGreenDivergenceDiagnostic`'s own hard-coded assertions are)

| cell_id/checkpoint | observable | macOS arm64 residual (max/3) | Linux x86_64 residual (max/3) | derived value | decision | reason |
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

**Cross-runner amplification cross-check (assertion 2, the conditioning cell's own seam):** FC/benign = 1.070e-12 / max(0.000e+00, 1e-15) on macOS arm64; 1.070e-12 / max(0.000e+00, 1e-15) on Linux x86_64 -- both >= 10x, reconfirming Event 1's finding on a second, independent platform (Linux x86_64 vs macOS arm64).

### Timing (informational -- the benchmark artifact owns the authoritative table)

- `module_total_seconds` (3 invocations, macOS arm64): 0.557, 0.545, 0.548
- `module_total_seconds` (3 invocations, Linux x86_64): 0.632, 0.626, 0.638
- `python -m unittest tests.test_rpa_flex_equivalence_table` process wall time: macOS arm64 0.85s, Linux x86_64 0.96s (both far under the 120s freeze-time budget; the 4-gating-runner CI sample is measured separately).

---

## Event 3 -- freeze calibration on the four gating continuous-integration runners

- **Date:** 2026-08-19
- **Commit:** 8144bf3f9e9539bad4759a2fbd1b24f52f7bef33 (the measured source commit; every one of the 16 samples below was produced from a checkout of exactly this revision, and each measurement sample's own `source_sha` field was checked equal to it).
- **Stage:** the four gating continuous-integration runners (`ubuntu-latest` x Python 3.9 / 3.10 / 3.11 / 3.12), `.github/workflows/equivalence-calibration.yml`, **workflow run 32204319966, attempt 1** -- the first and only iteration of the calibration loop. All four jobs completed successfully; there are no superseded runs and no superseded attempts.
- **Runner (Python 3.9):** Linux-6.17.0-1022-azure-x86_64-with-glibc2.39, Python 3.9.25, numpy 1.26.4, scipy 1.13.1
- **Runner (Python 3.10):** Linux-6.17.0-1022-azure-x86_64-with-glibc2.39, Python 3.10.21, numpy 1.26.4, scipy 1.15.3
- **Runner (Python 3.11):** Linux-6.17.0-1022-azure-x86_64-with-glibc2.39, Python 3.11.16, numpy 1.26.4, scipy 1.17.1
- **Runner (Python 3.12):** Linux-6.17.0-1022-azure-x86_64-with-glibc2.39, Python 3.12.13, numpy 1.26.4, scipy 1.17.1

**Protocol:** `python -m tests.equivalence_measure` run 3x on each of the four runners (12 measurement samples), plus one timed `python -m unittest tests.test_rpa_flex_equivalence_table` process per runner (4 timing samples) -- 16 artifacts in total. Each job runs from a fresh checkout with `OPENBLAS_NUM_THREADS`/`OMP_NUM_THREADS`/`MKL_NUM_THREADS`/`NUMEXPR_NUM_THREADS` pinned to 1 (the pin Event 2's segfault finding motivated; no crash occurred on any of the 12 invocations here). Every measurement sample: 90 emitted JSON lines plus the workflow's appended `process_seconds` line, 0 `"error"` records, `source_sha` == `8144bf3f9e9539bad4759a2fbd1b24f52f7bef33`. `python3 -m tests.equivalence_freeze_check <artifacts> --source-sha 8144bf3f9e9539bad4759a2fbd1b24f52f7bef33` reports **VALIDATION OK (12 measurement + 4 unittest samples)** over the set.

**Bound rule applied (the gating-runner variant, replacing Event 2's two-development-machine variant):** `raw_max` = the MAX over all 12 measurement samples (4 runners x 3 invocations); `floored = max(raw_max, 1e-15)`; `bound = 10 * floored` rounded UP to a power of ten; the bound must land at or under the mapped policy ceiling. STOP condition (`raw_max > ceiling`): did not occur -- 0 of 42 pairs. Margin-insufficient condition (`bound > ceiling`): did not occur -- 0 of 42 pairs. Exactly one pair lands EXACTLY at its ceiling (`bound == ceiling`, which the rule permits): the conditioning cell's `chiq`, as in Event 2 -- expected, since that cell is deliberately chosen to amplify.

**Outcome: the bounds freeze.** Recomputing all 42 bounds from the gating-runner MAX gives, in every case, exactly the atol already recorded in `tests/equivalence_cells.py` -- **no bound changed**. No cell's residual exceeds its ceiling, so there is no diverging cell and no cell straddling the boundary; no `Diverges` classification changes; the 120s freeze-time budget holds with a MAX `unittest_module_process_seconds` of **1.466 s** (see `tests/equivalence_benchmark.md`, Section 2, for the timing side). Because nothing measurement-affecting therefore had to change, iteration 1 of the calibration loop converged and the registry's `PROVENANCE` record moves to `status = "frozen"`, naming this commit and this run.

### Comparison-cell frozen atols (42 rows = 21 Equiv cells x 2 observables; every Equiv cell in the registry)

Residual columns are the MIN and MAX over all 12 measurement samples; the bound rule uses the MAX column.

| cell_id | observable | gating-runner MIN residual | gating-runner MAX residual | frozen atol | decision | reason |
|---|---|---|---|---|---|---|
| `general.ring.onsite_coulombintra.fixedmu` | chi0q | 9.714e-17 | 9.714e-17 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12 |
| `general.ring.onsite_coulombintra.fixedmu` | chiq | 1.943e-16 | 2.220e-16 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12 |
| `general.ring.onsite_coulombinter.fixedmu` | chi0q | 9.714e-17 | 9.714e-17 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12 |
| `general.ring.onsite_coulombinter.fixedmu` | chiq | 1.110e-16 | 1.110e-16 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12 |
| `general.ring.onsite_hund.fixedmu` | chi0q | 9.714e-17 | 9.714e-17 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12 |
| `general.ring.onsite_hund.fixedmu` | chiq | 1.110e-16 | 1.110e-16 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12 |
| `general.ring.onsite_ising.fixedmu` | chi0q | 9.714e-17 | 9.714e-17 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12 |
| `general.ring.onsite_ising.fixedmu` | chiq | 1.110e-16 | 1.110e-16 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12 |
| `general.ring.onsite_exchange.fixedmu` | chi0q | 9.714e-17 | 9.714e-17 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12 |
| `general.ring.onsite_exchange.fixedmu` | chiq | 9.714e-17 | 9.714e-17 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12 |
| `general.ring.onsite_pairhop.fixedmu` | chi0q | 9.714e-17 | 9.714e-17 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12 |
| `general.ring.onsite_pairhop.fixedmu` | chiq | 9.715e-17 | 9.715e-17 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12 |
| `general.ring.onsite_pairlift.fixedmu` | chi0q | 9.714e-17 | 9.714e-17 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12 |
| `general.ring.onsite_pairlift.fixedmu` | chiq | 9.714e-17 | 9.714e-17 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12 |
| `general.ring.onsite_u_v_hund.mu` | chi0q | 9.714e-17 | 9.714e-17 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `general.ring.onsite_u_v_hund.mu` | chiq | 2.776e-16 | 2.776e-16 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `general.ring.onsite_full_kanamori.mu` | chi0q | 9.714e-17 | 9.714e-17 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `general.ring.onsite_full_kanamori.mu` | chiq | 1.943e-16 | 2.498e-16 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `reduced.ring.onsite_coulombintra.spinfree.mu` | chi0q | 9.714e-17 | 9.714e-17 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `reduced.ring.onsite_coulombintra.spinfree.mu` | chiq | 3.053e-16 | 3.608e-16 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `reduced.ring.onsite_coulombinter.spinfree.mu` | chi0q | 9.714e-17 | 9.714e-17 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `reduced.ring.onsite_coulombinter.spinfree.mu` | chiq | 1.249e-16 | 1.388e-16 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `reduced.ring.onsite_hund.spinfree.mu` | chi0q | 9.714e-17 | 9.714e-17 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `reduced.ring.onsite_hund.spinfree.mu` | chiq | 1.110e-16 | 1.110e-16 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `reduced.ring.onsite_ising.spinfree.mu` | chi0q | 9.714e-17 | 9.714e-17 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `reduced.ring.onsite_ising.spinfree.mu` | chiq | 1.110e-16 | 1.249e-16 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `reduced.ring.onsite_pairlift.spinfree.mu` | chi0q | 9.714e-17 | 9.714e-17 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `reduced.ring.onsite_pairlift.spinfree.mu` | chiq | 9.714e-17 | 9.714e-17 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `reduced.ring.onsite_coulombintra.spindiag.mu` | chi0q | 1.249e-16 | 1.249e-16 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `reduced.ring.onsite_coulombintra.spindiag.mu` | chiq | 4.719e-16 | 4.719e-16 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `reduced.ring.onsite_coulombinter.spindiag.mu` | chi0q | 1.249e-16 | 1.249e-16 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `reduced.ring.onsite_coulombinter.spindiag.mu` | chiq | 1.388e-16 | 1.388e-16 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `general.ring.offsite_coulombinter_sameorb.mu` | chi0q | 4.899e-15 | 4.899e-15 | 1.0e-13 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `general.ring.offsite_coulombinter_sameorb.mu` | chiq | 2.026e-14 | 2.026e-14 | 1.0e-12 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `reduced.ring.offsite_coulombinter.mu` | chi0q | 4.899e-15 | 4.899e-15 | 1.0e-13 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `reduced.ring.offsite_coulombinter.mu` | chiq | 3.561e-14 | 3.561e-14 | 1.0e-12 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `general.ring.onsite_coulombinter.subshape.mu` | chi0q | 9.715e-17 | 9.715e-17 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `general.ring.onsite_coulombinter.subshape.mu` | chiq | 1.110e-16 | 1.110e-16 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `general.ring.onsite_coulombinter.coefftail.mu` | chi0q | 1.110e-16 | 1.110e-16 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `general.ring.onsite_coulombinter.coefftail.mu` | chiq | 1.388e-16 | 1.388e-16 | 1.0e-14 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `general.ring.offsite_coulombinter.conditioning.mu` | chi0q | 6.432e-13 | 6.432e-13 | 1.0e-11 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 |
| `general.ring.offsite_coulombinter.conditioning.mu` | chiq | 3.316e-12 | 3.373e-12 | 1.0e-10 | freeze unchanged | 10x the gating-runner MAX, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10 (== ceiling, which the rule permits) |

The two largest ordinary residuals are the off-site `CoulombInter` pair -- `general.ring.offsite_coulombinter_sameorb.mu` chiq at 2.026157e-14 and `reduced.ring.offsite_coulombinter.mu` chiq at 3.561041e-14 -- which round to the same 1.0e-12 bound. The conditioning cell is larger still by design: chi0q 6.432077e-13 -> 1.0e-11, chiq 3.373345e-12 -> 1.0e-10, the one bound sitting exactly on its ceiling.

### Divergence-diagnostic checkpoints (10 rows = 5 assertions x 2 fixtures; informational -- reconfirms Events 1 and 2 on the four gating runners; no registry atol is gated by this table)

Residual columns are the MIN and MAX of the four runner-local MAXima.

| cell_id/checkpoint | observable | gating-runner MIN residual | gating-runner MAX residual | derived value | decision | reason |
|---|---|---|---|---|---|---|
| diagnostic assertion 1 | benign | 0.000e+00 | 0.000e+00 | n/a | accept | within the `mu_diag` ceiling (1e-10) on all four runners |
| diagnostic assertion 1 | fc | 0.000e+00 | 0.000e+00 | n/a | accept | within the `mu_diag` ceiling (1e-10) on all four runners |
| diagnostic assertion 2 | benign | 0.000e+00 | 0.000e+00 | n/a | accept | within the `mu_diag` ceiling (1e-10) on all four runners |
| diagnostic assertion 2 | fc | 1.070e-12 | 1.070e-12 | n/a | accept | within the `mu_diag` ceiling (1e-10) on all four runners |
| diagnostic assertion 3 | benign | 1.119e-16 | 1.119e-16 | n/a | accept | within the `green_diag` ceiling (1e-10) on all four runners |
| diagnostic assertion 3 | fc | 2.289e-16 | 2.289e-16 | n/a | accept | within the `green_diag` ceiling (1e-10) on all four runners |
| diagnostic assertion 4 | benign | 1.119e-16 | 1.119e-16 | n/a | accept | within the `green_diag` ceiling (1e-10) on all four runners |
| diagnostic assertion 4 | fc | 7.850e-17 | 7.850e-17 | n/a | accept | within the `green_diag` ceiling (1e-10) on all four runners |
| diagnostic assertion 5 | benign | 9.714e-17 | 9.714e-17 | n/a | accept | within the `chi0q_mu` ceiling (1e-10) on all four runners |
| diagnostic assertion 5 | fc | 5.586e-17 | 1.111e-16 | n/a | accept | within the `chi0q_mu` ceiling (1e-10) on all four runners |

**Cross-runner amplification cross-check (assertion 2, the conditioning cell's own seam):** FC/benign = 1.070e-12 / max(0.000e+00, 1e-15) = **1.07e+03x** on every one of the four gating runners, identical across Python 3.9, 3.10, 3.11 and 3.12 -- the >= 10x conditioning criterion Event 1 established on one development machine and Event 2 reconfirmed on a second now holds on all four gating runners as well. The only checkpoint that varies at the precision shown across the four runners is assertion 5 on the FC fixture (5.586e-17 on Python 3.9/3.10/3.12, 1.111e-16 on Python 3.11); both values sit ~1e+05x under the `chi0q_mu` ceiling, so nothing turns on the difference. Assertion 5 on the benign fixture also differs across runners, but far below the precision these tables carry -- 9.714456569e-17 on Python 3.9/3.10/3.12 against 9.714496443e-17 on Python 3.11 -- so the three-digit columns above render it identically.

### Timing (informational -- `tests/equivalence_benchmark.md`, Section 2, owns the authoritative table)

- `module_total_seconds` (MAX of 3 invocations, per runner): Python 3.9 0.668, 3.10 0.840, 3.11 0.966, 3.12 0.986
- `unittest_module_process_seconds` (1 sample per runner): Python 3.9 0.983, 3.10 1.022, 3.11 1.112, 3.12 1.466 -- MAX **1.466 s** against the 120 s freeze-time budget, roughly 82x under it.

---

## Event 4 -- conditioning transfer-gain experiment (cell 38, replacing the retired amplification-ratio criterion)

- **Date:** 2026-08-29
- **Commit (dev-provisional, NOT the frozen PROVENANCE revision --
  see the note below):** 6c973db32287d6e9f22aedd791352b24d8e0e9a0 (the
  working tree these measurements were taken against; see the commit
  this file ships in for the exact source SHA once committed, per the
  same convention Event 1 uses). Deliberately NOT recorded with the
  bare `- **Commit:**` label Events 1-3 use: this section is a
  development-machine measurement appended AFTER Event 3's freeze, and
  `TestBenchmarkRegistryTie.test_frozen_provenance_matches_the_
  calibration_log_commit` reads the file's LAST bare `- **Commit:**`
  line as the calibration log's authoritative frozen revision -- that
  must keep pointing at Event 3's `8144bf3f9e9539bad4759a2fbd1b24f52f7bef33`
  until a future freeze supersedes it.
- **Stage:** development machine (provisional)
- **Runner:** macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13,
  numpy 2.4.6, scipy 1.17.1

**Objective:** `RPA._find_mu` was driven to a round-off particle-number
residual (a preceding commit on this branch), which collapses
diagnostic assertion 2 -- the seam Events 1-3's `>= 10x` amplification
criterion carried -- to ~0 on both the benign and FC fixtures. That
criterion is RETIRED as of this commit: it gated a root-finder
artifact the fix removed, and Events 1-3's amplification narratives
above are historical (accurate at the time, describing a seam that no
longer amplifies). Cell 38's conditioning obligation is replaced by a
deterministic perturbation transfer-gain experiment
(`TestConditioningTransferGain`): `gain = max|chiq(chi0 + eps*E) -
chiq(chi0)| / eps`, with `E` a one-hot real perturbation at the
in-test argmax-cond (freq, k) direction of `cond_2(1 + chi0 W)`, run
independently on the benign and FC fixtures and requiring a `>= 10x`
FC/benign gain contrast.

| cell_id / checkpoint | observable | benign value (cell 8, T=2.0) | FC value (T=0.2) | derived value | decision | reason |
|---|---|---|---|---|---|---|
| `TestConditioningTransferGain` argmax-cond direction | (l0, k, cond) | (16, 10, 3.6249143215242605) -- benign's argmax IS the zero-bosonic-frequency slot | (127, 10, 70.48956521627714), tied with (129, 10, 70.48956521614214) at ~1.9e-12 relative -- FC's argmax is NOT Omega=0 (FC's zero slot, l0=128, is cond 1.9718478092582121, ~35x better conditioned); lexicographic tie-break selects l0=127 | -- | accept | reconfirms the argmax-cond direction the test's own `_perturbation_target` selects, on-branch |
| `TestConditioningTransferGain` transfer gain | gain (max\|d chiq\| / eps) | 1.7551392977877 @ (16, 10) | 1277.6858808294378 @ (127, 10) | ratio = **727.9683626478662x** (~728x) | accept | clears the `>= 10x` FC/benign contrast bar by ~73x; linear across eps 1e-6..1e-10 |
| `TestConditioningTransferGain` propagated builder difference | chiq (max\|rpa-built - flex-built\|) | 2.2204481854326364e-16 | 7.376256836354663e-14 | worst = 7.376256836354663e-14 | accept | far inside the derived `chiq_propagated` ceiling below |

**Derived provisional `POLICY_CEILINGS` keys** (all marked `# provisional (dev) -- finalized by the CI calibration, see calibration log` in `tests/equivalence_cells.py`):

- `chiq_gain` = `10**ceil(log10(10 * 1277.6858808294378))` = **1e5** (upper ceiling on the FC gain).
- `chiq_gain_fc_min` = largest power of ten `<= 1277.6858808294378 / 10` = **1e2** (lower band on the FC gain).
- `chiq_propagated` = `max(1e-15, 10**ceil(log10(10 * 7.376256836354663e-14)))` = **1e-12** (upper ceiling on the worst measured propagated rpa-vs-flex chiq builder difference across both fixtures).

**Decision:** ACCEPT the transfer-gain experiment and its three derived
provisional keys as measured, on this development machine, on the
commit above. Not yet a CI freeze -- these are dev-provisional values
pending the four-gating-runner calibration pass (the same pattern
Events 1-2 preceded Event 3's freeze).

**Note on Events 1-3:** the `>= 10x` amplification criterion those
events measured and froze (assertion 2, the mu seam) is RETIRED by
this commit -- the root-finder artifact it depended on was removed by
the preceding `RPA._find_mu` fix, so assertion 2 now sits at round-off
on both fixtures and the amplification claim no longer holds. Events
1-3's narratives remain accurate historical records of what was true
at the time; they are not amended.

---

## Event 5 -- two-stage CI freeze of every remaining ceiling (Task 7)

- **Date:** 2026-08-30
- **Commit (dev+CI measured, NOT the frozen PROVENANCE revision --
  see the note at the end of this Event):** 6f588ead496b77a077f3d8a642224aa8b3c8a9cf
  (the working tree both stages below were measured against; deliberately
  NOT recorded with the bare `- **Commit:**` label Events 1-3 use, for
  the same reason Event 4 avoided it -- `TestBenchmarkRegistryTie`'s
  provenance-tie tests read the file's LAST bare `- **Commit:**` line
  as authoritative, and that must keep pointing at Event 3's
  `8144bf3f9e9539bad4759a2fbd1b24f52f7bef33` until a future freeze that
  also re-measures and re-freezes the ORIGINAL 42-cell candidate-atol
  table (Event 2/3's own subject) supersedes it. See the closing note.)
- **Stage 1 (dev machine):** macOS-26.5.2-arm64-arm-64bit-Mach-O,
  Python 3.13.13, numpy 2.4.6, scipy 1.17.1, BLAS/LAPACK = Accelerate
  (`numpy.show_config(mode="dicts")`: `blas.name`/`lapack.name` =
  `"accelerate"`, both `detection method = "system"`). Two
  `python -m tests.equivalence_measure` invocations
  (`/tmp/calib-dev-1.json`, `/tmp/calib-dev-2.json`), both 113 lines, 0
  `"error"` records, `source_sha` == `6f588ead496b77a077f3d8a642224aa8b3c8a9cf`
  on both.
- **Stage 2 (the four gating continuous-integration runners):**
  `ubuntu-latest` x Python 3.9/3.10/3.11/3.12,
  `.github/workflows/equivalence-calibration.yml`, **workflow run
  33278664447, attempt 1** -- the first and only iteration of this
  calibration loop; all four jobs completed successfully (16
  artifacts: 12 measurement + 4 unittest-timing), and
  `python3 -m tests.equivalence_freeze_check /tmp/calib-ci --source-sha
  6f588ead496b77a077f3d8a642224aa8b3c8a9cf` reports **VALIDATION OK
  (12 measurement + 4 unittest samples, source_sha
  6f588ead496b77a077f3d8a642224aa8b3c8a9cf)**. Runner descriptors (from
  each sample's own `runner` metadata block; BLAS/LAPACK is NOT
  obtainable from these artifacts -- `equivalence_measure` does not
  call `numpy.show_config()` -- but all four use the same pip-wheel
  numpy 1.26.4, whose manylinux wheels are OpenBLAS-backed):
  - Python 3.9: Linux-6.17.0-1022-azure-x86_64-with-glibc2.39, Python
    3.9.25, numpy 1.26.4, scipy 1.13.1
  - Python 3.10: Linux-6.17.0-1022-azure-x86_64-with-glibc2.39, Python
    3.10.21, numpy 1.26.4, scipy 1.15.3
  - Python 3.11: Linux-6.17.0-1022-azure-x86_64-with-glibc2.39, Python
    3.11.16, numpy 1.26.4, scipy 1.17.1
  - Python 3.12: Linux-6.17.0-1022-azure-x86_64-with-glibc2.39, Python
    3.12.14, numpy 1.26.4, scipy 1.17.1

**Objective:** this branch's `RPA._find_mu` fix (the mu/Green-seam fix
Task 3 landed) drives the mu root-finder to round-off. Task 3's
before/after mu-seam record (transcribed below -- the same mu-seam
diagnostic, `measure: "M1"`, on the FC-family fixtures at descending
T) shows the seam collapsing cleanly:

| fixture | T | `dmu` before (red) | `dmu` after (green) |
|---|---|---|---|
| `benign_cell8` | 2.0 | 0.0 | 0.0 |
| `fc_cell38_T0.2` | 0.2 | 1.070254995738651e-12 | 0.0 |
| `fc_cell38_T0.1` | 0.1 | 1.5136780717739384e-12 | 1.1102230246251565e-16 |
| `fc_cell38_T0.05` | 0.05 | 1.6872059305228504e-12 | 0.0 |
| `fc_cell38_T0.025` | 0.025 | 7.233103005432895e-13 | 1.1102230246251565e-16 |

This is the mechanism behind every ceiling this Event tightens: the
diagnostic checkpoints and the mu-coupled comparison cells that used to
carry the old root-finder residual now sit at round-off, so their
POLICY_CEILINGS entries -- most still at their original dev-only
provisional values, several capped at 1e-10 as a wide placeholder --
can be re-derived tight against the fix. This Event does that in two
stages: Stage 1 re-derives every remaining diagnostic-apparatus ceiling
(`mu_diag`, `green_diag`, `counter_cross_nd_le2`, `counter_cross_geev`,
`mu_number_residual`, `green_dyson`) plus the transfer-gain trio
(`chiq_gain`, `chiq_gain_fc_min`, `chiq_propagated`); Stage 2
re-derives the two remaining mu-coupled comparison ceilings (`chi0q_mu`,
`chiq_mu`) from the worst mu-coupled per-cell residual in the same
report. `chi0q_fixed`/`chiq_fixed` are out of scope (unaffected by the
mu seam; stay 1e-12 per the plan).

**Derivation rule (per key, both stages):** worst = MAX over {the dev
machine, all 4 CI runners} of the measured residual; ceiling = `1e-15
if worst == 0 else max(1e-15, 10**ceil(log10(10 * worst)))`.
`chiq_gain_fc_min` uses the MIN/floor variant instead (largest power of
ten `<= min_measured_fc_gain / 10`); see its own row below.

### Stage 1 -- diagnostic-apparatus checkpoints (33 dev + 12x33 CI records, condensed to one row per (metric, fixture))

| cell_id/checkpoint | observable | dev-machine MAX residual | CI-runner MAX residual (over 4x3) | derived value | decision | reason |
|---|---|---|---|---|---|---|
| assertion1 (`mu_diag`) | benign | 0.000000e+00 | 0.000000e+00 | -- | accept | contributes to `mu_diag`'s worst |
| assertion1 (`mu_diag`) | fc | 0.000000e+00 | 0.000000e+00 | -- | accept | contributes to `mu_diag`'s worst |
| assertion1 (`mu_diag`) | geev | 0.000000e+00 | 0.000000e+00 | -- | accept | contributes to `mu_diag`'s worst |
| assertion2 (`mu_diag`, the mu seam) | benign | 0.000000e+00 | 6.427495e-17 | -- | accept | nonzero on one CI runner (OpenBLAS single-ULP noise, per the WARNING carried into this task -- the identity is structural at Sigma=0) |
| assertion2 (`mu_diag`) | fc | 0.000000e+00 | 0.000000e+00 | -- | accept | contributes to `mu_diag`'s worst |
| assertion2 (`mu_diag`) | geev | 1.765081e-16 | 1.765081e-16 | -- | accept | the WORST value feeding `mu_diag` below |
| assertion3 (`green_diag`, Green seam, bare) | benign | 1.119071e-16 | 1.118863e-16 | -- | accept | contributes to `green_diag`'s worst |
| assertion3 (`green_diag`) | fc | 6.206335e-17 | 7.850462e-17 | -- | accept | contributes to `green_diag`'s worst |
| assertion3 (`green_diag`) | geev | 4.163336e-16 | 3.775166e-16 | -- | accept | contributes to `green_diag`'s worst |
| assertion4 (`green_diag`, composed seam) | benign | 1.119071e-16 | 1.118863e-16 | -- | accept | contributes to `green_diag`'s worst |
| assertion4 (`green_diag`) | fc | 6.206335e-17 | 7.850462e-17 | -- | accept | contributes to `green_diag`'s worst |
| assertion4 (`green_diag`) | geev | 3.775166e-16 | **4.518280e-16** | -- | accept | the WORST value feeding `green_diag` below (CI Python 3.10) |
| assertion5 (`chi0q_mu`, downstream chi0q) | benign | 9.714461e-17 | 9.714477e-17 | -- | accept | contributes to `chi0q_mu`'s worst |
| assertion5 (`chi0q_mu`) | fc | 5.551275e-17 | 5.551320e-17 | -- | accept | contributes to `chi0q_mu`'s worst |
| assertion5 (`chi0q_mu`) | geev | **2.498237e-16** | 1.127572e-16 | -- | accept | the WORST diagnostic value feeding `chi0q_mu` below (dev machine); see Stage 2 for the full `chi0q_mu` derivation |
| counter_cross_at_mu_rpa (`counter_cross_nd_le2`) | benign | 0.000000e+00 | 0.000000e+00 | -- | accept | structural zero at Sigma=0 |
| counter_cross_at_mu_rpa (`counter_cross_nd_le2`) | fc | 0.000000e+00 | 0.000000e+00 | -- | accept | structural zero |
| counter_cross_at_mu_rpa (`counter_cross_geev`) | geev | 0.000000e+00 | 0.000000e+00 | -- | accept | structural zero |
| counter_cross_at_mu_flex (`counter_cross_nd_le2`) | benign | 0.000000e+00 | 0.000000e+00 | -- | accept | structural zero |
| counter_cross_at_mu_flex (`counter_cross_nd_le2`) | fc | 0.000000e+00 | 0.000000e+00 | -- | accept | structural zero |
| counter_cross_at_mu_flex (`counter_cross_geev`) | geev | 0.000000e+00 | 0.000000e+00 | -- | accept | structural zero |
| number_residual_rpa (`mu_number_residual`) | benign/fc/geev | 0.000000e+00 | 0.000000e+00 | -- | accept | structural zero on all three fixtures, both stages |
| number_residual_flex (`mu_number_residual`) | benign/fc/geev | 0.000000e+00 | 0.000000e+00 | -- | accept | structural zero on all three fixtures, both stages |
| dyson_residual_eigenbasis (`green_dyson`) | benign | 8.906697e-16 | 8.906697e-16 | -- | accept | contributes to `green_dyson`'s worst |
| dyson_residual_eigenbasis (`green_dyson`) | fc | 2.286149e-16 | 2.288783e-16 | -- | accept | contributes to `green_dyson`'s worst |
| dyson_residual_eigenbasis (`green_dyson`) | geev | 8.899115e-16 | 8.881871e-16 | -- | accept | contributes to `green_dyson`'s worst |
| dyson_residual_inv (`green_dyson`) | benign | 3.554448e-16 | 3.554448e-16 | -- | accept | contributes to `green_dyson`'s worst |
| dyson_residual_inv (`green_dyson`) | fc | 2.307666e-16 | 2.288783e-16 | -- | accept | contributes to `green_dyson`'s worst |
| dyson_residual_inv (`green_dyson`) | geev | **9.930137e-16** | 9.036561e-16 | -- | accept | the WORST value feeding `green_dyson` below (dev machine) |

**WARNING carried from review, confirmed here:** `counter_cross_*` and
`mu_number_residual` measure exactly 0.0 on the dev machine, both
stages, every fixture -- the identity is structural at Sigma=0. CI DID
show a nonzero value on this one checkpoint (assertion 2, benign
fixture): exactly 6.42749509214523e-17 on Python 3.9 (all 3
invocations) and Python 3.12 (all 3 invocations), while Python 3.10
and Python 3.11 measured exactly 0.0 (all 3 invocations each) --
consistent with OpenBLAS single-ULP noise on an O(N) quantity that
happens to round differently on some interpreter builds, exactly the
possibility flagged before this run. It does not touch
`counter_cross_*`/`mu_number_residual` themselves (those stayed
exactly 0.0 everywhere, all 4 runners), only `mu_diag` (already
dominated by the geev fixture's 1.765081e-16 either way).

**Stage 1 -- per-key derived ceilings:**

| key | worst (MAX over dev+CI) | derivation | ceiling | prior value | change |
|---|---|---|---|---|---|
| `mu_diag` | 1.765081e-16 (assertion2, geev, dev) | `10**ceil(log10(1.765081e-15))` | **1e-14** | 1e-14 | unchanged |
| `green_diag` | 4.518280e-16 (assertion4, geev, CI 3.10) | `10**ceil(log10(4.518280e-15))` | **1e-14** | 1e-10 | **tightened 1e4x** |
| `counter_cross_nd_le2` | 0.0 | worst==0 floor | **1e-15** | 1e-15 | unchanged |
| `counter_cross_geev` | 0.0 | worst==0 floor | **1e-15** | 1e-15 | unchanged |
| `mu_number_residual` | 0.0 | worst==0 floor | **1e-15** | 1e-15 | unchanged |
| `green_dyson` | 9.930137e-16 (dyson_residual_inv, geev, dev) | `10**ceil(log10(9.930137e-15))` | **1e-14** | 1e-14 | unchanged |

### Transfer-gain trio (`chiq_gain`, `chiq_gain_fc_min`, `chiq_propagated`) -- gate-passage freeze

The calibration workflow runs `equivalence_measure` (Stage 1/2 above)
plus a single timed `python -m unittest tests.test_rpa_flex_equivalence_table`
per runner; `TestConditioningTransferGain` (the transfer-gain
experiment) lives in that unittest module, not in
`equivalence_measure`, so none of the 33x(1+12) diagnostic records
above cover it and no raw CI measurement of the gain exists yet. Per
the controller's ruling: since all four gating runners ran
`tests.test_rpa_flex_equivalence_table` GREEN at the current
dev-provisional values (workflow run 33278664447 attempt 1 -- every
"Equivalence table unittest timing" step succeeded, and a nonzero exit
there fails the job per the workflow's `exit $status`), that IS
sufficient CI evidence to freeze these three keys at their provisional
values -- gate-passage rather than a raw cross-runner residual table. A
raw-measurement CI record for the gain (each runner reporting its own
gain/cond/`chiq_propagated`) is deferred to a future calibration pass.

Re-measured on the dev machine for this freeze (`TestConditioningTransferGain.setUpClass`
run standalone, 2026-08-30) to confirm nothing drifted since Event 4:

| cell_id / checkpoint | observable | benign value (cell 8, T=2.0) | FC value (T=0.2) | derived value | decision | reason |
|---|---|---|---|---|---|---|
| `TestConditioningTransferGain` argmax-cond direction | (l0, k, cond) | (16, 10, 3.6249143215242605) | (127, 10, 70.48956521627714) | -- | accept | bit-identical to Event 4's dev measurement -- no drift from the intervening FC-argmax-cond-direction/FLEX-SCF-trajectory fixes (Tasks 5-6), which touched documentation/pinning, not this computation |
| `TestConditioningTransferGain` transfer gain | gain (max\|d chiq\| / eps) | 1.7551392977877 | 1277.6858808294378 | ratio 727.9683626478662x | accept | bit-identical to Event 4; clears the `>=10x` FC/benign contrast bar by ~73x |
| `TestConditioningTransferGain` propagated builder difference | chiq (max\|rpa-built - flex-built\|) | 2.2204481854326364e-16 | 7.376256836354663e-14 | worst = 7.376256836354663e-14 | accept | bit-identical to Event 4 |

| key | worst | derivation | ceiling | prior value | change |
|---|---|---|---|---|---|
| `chiq_gain` | FC gain 1277.6858808294378 | `10**ceil(log10(10 * 1277.6858808294378))` | **1e5** | 1e5 (dev-provisional) | frozen unchanged -- CI evidence is gate-passage (see above) |
| `chiq_gain_fc_min` | FC gain 1277.6858808294378 | largest power of ten `<= 1277.6858808294378 / 10` | **1e2** | 1e2 (dev-provisional) | frozen unchanged -- gate-passage |
| `chiq_propagated` | 7.376256836354663e-14 (FC) | `10**ceil(log10(7.376256836354663e-13))` | **1e-12** | 1e-12 (dev-provisional) | frozen unchanged -- gate-passage |

### Stage 2 -- mu-coupled comparison cells (`chi0q_mu`, `chiq_mu`)

Every `(cell, chi0q|chiq)` pair whose registry entry maps to the
`chi0q_mu`/`chiq_mu` ceiling keys (every mu-coupled comparison cell
except the fixed-mu `_fixedmu` group, which uses `chi0q_fixed`/
`chiq_fixed` and is out of scope), dev-machine MAX (primary residual)
vs. the CI-runner MAX over all 12 measurement samples (secondary
residual), in the canonical 7-column shape:

| cell_id/checkpoint | observable | dev-machine MAX (primary) | CI-runner MAX 4x3 (secondary) | derived value (worst) | decision | reason |
|---|---|---|---|---|---|---|
| `general.ring.offsite_coulombinter.conditioning.mu` (cell 38) | chi0q | 5.551275e-17 | 5.551320e-17 | 5.551320e-17 | recalibrate | stale pre-#160 literal would raise `_candidate_atol` against the tightened ceiling; see the per-cell recalibration table below |
| `general.ring.offsite_coulombinter.conditioning.mu` (cell 38) | chiq | **7.294815e-14** | 6.843859e-14 | **7.294815e-14** | recalibrate | worst value feeding `chiq_mu` below; stale literal recalibrated (see below) |
| `general.ring.offsite_coulombinter_sameorb.mu` | chi0q | 5.551921e-17 | 4.167915e-17 | 5.551921e-17 | recalibrate | stale pre-#160 literal recalibrated (see below) |
| `general.ring.offsite_coulombinter_sameorb.mu` | chiq | 2.220865e-16 | 1.944572e-16 | 2.220865e-16 | recalibrate | stale pre-#160 literal recalibrated (see below) |
| `general.ring.onsite_coulombinter.coefftail.mu` | chi0q | 1.110223e-16 | 1.110245e-16 | 1.110245e-16 | accept | existing literal atol already covers this residual |
| `general.ring.onsite_coulombinter.coefftail.mu` | chiq | 1.387779e-16 | 1.387795e-16 | 1.387795e-16 | accept | existing literal atol already covers this residual |
| `general.ring.onsite_coulombinter.subshape.mu` | chi0q | **1.387794e-16** | 9.714795e-17 | **1.387794e-16** | accept | worst per-cell value feeding `chi0q_mu` below (before the assertion5 diagnostic is folded in); existing literal atol already covers it |
| `general.ring.onsite_coulombinter.subshape.mu` | chiq | 1.665345e-16 | 1.110254e-16 | 1.665345e-16 | accept | existing literal atol already covers this residual |
| `general.ring.onsite_full_kanamori.mu` | chi0q | 9.714461e-17 | 9.714477e-17 | 9.714477e-17 | accept | existing literal atol already covers this residual |
| `general.ring.onsite_full_kanamori.mu` | chiq | 2.498003e-16 | 2.498006e-16 | 2.498006e-16 | accept | existing literal atol already covers this residual |
| `general.ring.onsite_u_v_hund.mu` | chi0q | 9.714461e-17 | 9.714477e-17 | 9.714477e-17 | accept | existing literal atol already covers this residual |
| `general.ring.onsite_u_v_hund.mu` | chiq | 2.498004e-16 | 2.775562e-16 | 2.775562e-16 | accept | existing literal atol already covers this residual |
| `reduced.ring.offsite_coulombinter.mu` | chi0q | 5.551921e-17 | 4.167915e-17 | 5.551921e-17 | recalibrate | stale pre-#160 literal recalibrated (see below) |
| `reduced.ring.offsite_coulombinter.mu` | chiq | 3.609021e-16 | 3.333702e-16 | 3.609021e-16 | recalibrate | stale pre-#160 literal recalibrated (see below) |
| `reduced.ring.onsite_coulombinter.spindiag.mu` | chi0q | 1.249090e-16 | 1.249428e-16 | 1.249428e-16 | accept | existing literal atol already covers this residual |
| `reduced.ring.onsite_coulombinter.spindiag.mu` | chiq | 1.526687e-16 | 1.387822e-16 | 1.526687e-16 | accept | existing literal atol already covers this residual |
| `reduced.ring.onsite_coulombinter.spinfree.mu` | chi0q | 9.714461e-17 | 9.714477e-17 | 9.714477e-17 | accept | existing literal atol already covers this residual |
| `reduced.ring.onsite_coulombinter.spinfree.mu` | chiq | 1.110225e-16 | 1.387780e-16 | 1.387780e-16 | accept | existing literal atol already covers this residual |
| `reduced.ring.onsite_coulombintra.spindiag.mu` | chi0q | 1.249090e-16 | 1.249428e-16 | 1.249428e-16 | accept | existing literal atol already covers this residual |
| `reduced.ring.onsite_coulombintra.spindiag.mu` | chiq | 3.886063e-16 | 4.718469e-16 | 4.718469e-16 | accept | existing literal atol already covers this residual |
| `reduced.ring.onsite_coulombintra.spinfree.mu` | chi0q | 9.714461e-17 | 9.714477e-17 | 9.714477e-17 | accept | existing literal atol already covers this residual |
| `reduced.ring.onsite_coulombintra.spinfree.mu` | chiq | 3.885783e-16 | 3.608233e-16 | 3.885783e-16 | accept | existing literal atol already covers this residual |
| `reduced.ring.onsite_hund.spinfree.mu` | chi0q | 9.714461e-17 | 9.714477e-17 | 9.714477e-17 | accept | existing literal atol already covers this residual |
| `reduced.ring.onsite_hund.spinfree.mu` | chiq | 9.714480e-17 | 1.110230e-16 | 1.110230e-16 | accept | existing literal atol already covers this residual |
| `reduced.ring.onsite_ising.spinfree.mu` | chi0q | 9.714461e-17 | 9.714477e-17 | 9.714477e-17 | accept | existing literal atol already covers this residual |
| `reduced.ring.onsite_ising.spinfree.mu` | chiq | 9.714461e-17 | 1.249006e-16 | 1.249006e-16 | accept | existing literal atol already covers this residual |
| `reduced.ring.onsite_pairlift.spinfree.mu` | chi0q | 9.714461e-17 | 9.714477e-17 | 9.714477e-17 | accept | existing literal atol already covers this residual |
| `reduced.ring.onsite_pairlift.spinfree.mu` | chiq | 9.714461e-17 | 9.714477e-17 | 9.714477e-17 | accept | existing literal atol already covers this residual |

`chi0q_mu`'s overall worst also includes the Stage-1 diagnostic
checkpoint (assertion 5, which is directly gated against `chi0q_mu`):
`max(1.387794e-16 [this table], 2.498237e-16 [assertion5/geev/dev]) =
2.498237e-16`. `chiq_mu` has no diagnostic checkpoint, so its worst is
this table's own maximum, 7.294815e-14 (cell 38's chiq, dev machine).

| key | worst (MAX over dev+CI, both this table and Stage 1's assertion5) | derivation | ceiling | prior value | change |
|---|---|---|---|---|---|
| `chi0q_mu` | 2.498237e-16 | `10**ceil(log10(2.498237e-15))` | **1e-14** | 1e-10 | **tightened 1e4x** |
| `chiq_mu` | 7.294815e-14 | `10**ceil(log10(7.294815e-13))` | **1e-12** | 1e-10 | **tightened 1e2x** |

**Per-cell atol decisions (Stage 2):** every mu-coupled cell's own
`_candidate_atol` bound was re-checked against the new, tighter
`chi0q_mu`=1e-14 / `chiq_mu`=1e-12 ceilings. All but three cleared with
their EXISTING two-development-machine literal residuals unchanged
(those literals, measured back in Event 2 before the #160 fix, were
already at round-off and did not shift enough to matter -- e.g. cell
8's chi0q literal 1.110e-16 still bounds to 1e-14, exactly at the new
ceiling, which the rule permits). Three cells' literals were STALE
enough to make `_candidate_atol` raise against the new ceilings (their
residual dropped by 2-4 orders of magnitude once the #160 fix landed,
but their registry literal still carried the pre-fix value) and were
RECALIBRATED to the dev+CI-worst values in this table, via a new
`_measured_equiv_recalibrated` helper (`tests/equivalence_cells.py`,
mirrors `_measured_equiv`'s bound rule exactly, just sourced from {dev
machine, CI-runner MAX} instead of {macOS, Linux dev machine}):

| cell_id | observable | old literal (pre-#160, Event 2) | new literal (dev, CI-MAX) | old atol | new atol | resists tightening? |
|---|---|---|---|---|---|---|
| `general.ring.offsite_coulombinter_sameorb.mu` | chi0q | 4.884981379996393e-15 / 4.898859387959854e-15 | 5.551921408376704e-17 / 4.1679153779129384e-17 | 1e-13 | **1e-14** | no (lands exactly on the new `chi0q_mu` ceiling -- `_candidate_atol` floors the raw residual at 1e-15 BEFORE the 10x multiply, since 5.55e-17 < 1e-15: floored 1e-15 -> x10 -> 1e-14) |
| `general.ring.offsite_coulombinter_sameorb.mu` | chiq | 2.0206059344954128e-14 / 2.0261571408243727e-14 | 2.220864572139768e-16 / 1.944571711743353e-16 | 1e-12 | **1e-14** | no |
| `reduced.ring.offsite_coulombinter.mu` | chi0q | 4.884981379996393e-15 / 4.898859387959854e-15 | 5.551921408376704e-17 / 4.1679153779129384e-17 | 1e-13 | **1e-14** | no (same 1e-15 pre-multiply floor as the sameorb cell above) |
| `reduced.ring.offsite_coulombinter.mu` | chiq | 3.552713730991191e-14 / 3.5610405641548485e-14 | 3.6090211750999277e-16 / 3.3337017449116085e-16 | 1e-12 | **1e-14** | no |
| `general.ring.offsite_coulombinter.conditioning.mu` (cell 38) | chi0q | 6.431799537609854e-13 / 6.432077093507842e-13 | 5.551275218863831e-17 / 5.551319785776992e-17 | 1e-11 | **1e-14** | no (same 1e-15 pre-multiply floor) |
| `general.ring.offsite_coulombinter.conditioning.mu` (cell 38) | chiq | 3.3066883022456364e-12 / 3.373345155219826e-12 | 7.294815206161733e-14 / 6.843858684819496e-14 | 1e-10 (== old ceiling) | **1e-12** (== new ceiling) | **no** -- lands exactly on the new `chiq_mu` ceiling, same "conditioning cell is the tightest" pattern Events 2/3 recorded against the old ceiling |

No cell in the registry resisted tightening (no cell's bound needed
its CEILING to hold it down after recalibration); all three
recalibrated cells' new chiq bounds land comfortably under the new
`chiq_mu` ceiling, except cell 38's chiq, which lands exactly on
`chiq_mu`=1e-12 by construction (it is the deliberately-chosen
worst-case conditioning cell). All three recalibrated cells' new chi0q
bounds land EXACTLY ON the new `chi0q_mu`=1e-14 ceiling (permitted
equality, not a margin-insufficient case) -- the same pattern already
acknowledged above for cell 8's existing chi0q literal, and for the
same reason: each cell's raw chi0q residual (5.2-5.6e-17) is below the
1e-15 floor `_candidate_atol` applies before the 10x multiply, so
every one of them bounds to the same `10 * 1e-15 = 1e-14`.
`chi0q_fixed`/`chiq_fixed` stay 1e-12 unchanged (out of scope; no
fixed-mu cell's residual moved).

### Summary: every `POLICY_CEILINGS` key, before -> after this Event

| key | prior value | frozen value | change |
|---|---|---|---|
| `mu_diag` | 1e-14 (dev-provisional) | **1e-14** | unchanged, provisional label removed |
| `green_diag` | 1e-10 (placeholder) | **1e-14** | tightened 1e4x |
| `chi0q_mu` | 1e-10 (placeholder) | **1e-14** | tightened 1e4x |
| `chiq_mu` | 1e-10 (placeholder) | **1e-12** | tightened 1e2x |
| `chi0q_fixed` | 1e-12 | **1e-12** | unchanged (out of scope) |
| `chiq_fixed` | 1e-12 | **1e-12** | unchanged (out of scope) |
| `counter_cross_nd_le2` | 1e-15 (dev-provisional) | **1e-15** | unchanged, provisional label removed |
| `counter_cross_geev` | 1e-15 (dev-provisional) | **1e-15** | unchanged, provisional label removed |
| `mu_number_residual` | 1e-15 (dev-provisional) | **1e-15** | unchanged, provisional label removed |
| `green_dyson` | 1e-14 (dev-provisional) | **1e-14** | unchanged, provisional label removed |
| `chiq_gain` | 1e5 (dev-provisional) | **1e5** | unchanged, frozen via CI gate-passage |
| `chiq_gain_fc_min` | 1e2 (dev-provisional) | **1e2** | unchanged, frozen via CI gate-passage |
| `chiq_propagated` | 1e-12 (dev-provisional) | **1e-12** | unchanged, frozen via CI gate-passage |

**Verification after applying every change above:**
`python -m unittest discover tests -v` and `pytest tests -q` both pass
in full (`tests/equivalence_calibration_log.md`'s own commit records
the exact counts; see the commit this Event ships in).

**Note on `PROVENANCE` (why it is untouched by this Event):**
`tests.equivalence_cells.PROVENANCE` names the source revision and
calibration run that CONFIRMED the original 42-cell candidate-atol
table (Events 1-3) -- it is tied, by
`test_frozen_provenance_names_the_most_recent_benchmark_commit` and
`test_frozen_provenance_matches_the_calibration_log_commit`
(`tests/test_rpa_flex_equivalence_table.py`), to BOTH
`tests/equivalence_benchmark.md`'s most recent section's commit line
AND this file's most recent BARE `- **Commit:**` line -- and moving it
would require updating `tests/equivalence_benchmark.md` too (out of
this task's file scope: no cell was added, removed, renamed, or had
its fixture enlarged, so the registry docstring's own maintenance
checklist does not call for a new benchmark section here). This Event
follows Event 4's discipline instead: its own Commit line is
deliberately non-bare, so it is invisible to both provenance-tie
parsers, and `PROVENANCE` keeps naming Event 3's revision/run
(`8144bf3f9e9539bad4759a2fbd1b24f52f7bef33`, `"32204319966 attempt
1"`) as the confirmation for that original table -- which this Event
does not reopen, aside from the three cells recalibrated above (their
own `ObservableSpec.provenance` strings carry this Event's own
commit/run identity directly, so nothing about them is left
unattributed).

### Addendum: the gate margin with the Newton polish disabled

Recorded by the merge-gate review measurement, 2026-08-30, as part of
this branch's final review. The question addressed: how much margin
does the `RPA._find_mu` Newton polish step actually buy over the
underlying `brentq` bracket search, and how much of the gate's
sensitivity comes from the polish versus from `brentq` itself?

Two comparisons were made against the diagnostic and comparison
checkpoints this Event calibrated above:

- **Polish disabled, `brentq` alone (`xtol=1e-14`), on every fixture
  in the diagnostic suite:** every checkpoint stays green except one --
  `number_residual_rpa` (the `mu_number_residual` metric) on the FC
  fixture measures 1.776e-15 against its 1e-15 ceiling, a single-ULP
  overshoot at the O(8) target for that residual. `mu_diag` and both
  `counter_cross_*` metrics, and every other fixture's
  `number_residual_rpa`, remain green with the polish off.
- **The full pre-branch regression (the state this Event's bisect
  identified as broken, i.e. before the mu/Green-seam fix landed) is
  still caught decisively either way:** assertion 2's mu-seam residual
  measures 1.07e-12 against the `mu_diag` ceiling of 1e-14 -- two
  decades of margin, unaffected by whether the polish step runs.

So the gate distinguishes two different things at two different
scales: brentq-vs-pre-branch-regression is a two-decade margin (the
polish is not what catches that), while polish-vs-brentq-alone is a
single-ULP margin on one metric on one fixture. The Newton polish
step's own contribution to the gate is real but narrow -- it is the
difference between `number_residual_rpa` passing and failing by one
ULP on the FC fixture, not the difference between catching or missing
the regression this Event's fix addresses.
