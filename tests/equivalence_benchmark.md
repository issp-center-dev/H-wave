# RPA/FLEX equivalence table -- benchmark

Section schema: header (commit, protocol, runner matrix), then a table `cell_id | max_seconds_over_runs_and_runners` covering EVERY cell (`tests.equivalence_measure` times all proof shapes -- rejections, construction-only rows, and the two multirun exceptions included), then `module_total_seconds_per_runner`, `process_overhead_seconds`, the `unittest_module_process_seconds` sample(s), and the MEASURED SOURCE SHA + the authoritative run identity. **This file grows one section per calibration stage**: Section 1 below is the CANDIDATE pass taken on two development machines (advisory proxy runners, 3 measurement invocations each). Section 2 is the FREEZE pass on the 4 gating continuous-integration runners (3 measurement invocations each plus the 4 unittest-timing samples the 120s budget actually gates), taken once a closed calibration loop produced a freezable measurement; it records the authoritative run identity and the rule used to choose it.

**Registry docstring maintenance checklist:** any future cell ADDITION or ENLARGEMENT (a new fixture, a heavier `Nmat`/`CellShape`, a new proof shape) invalidates this file's per-cell timing table for the changed/added cell(s) and REQUIRES a new benchmark section here (never an in-place edit of a PAST section -- each section is an immutable record of one calibration event, exactly like `tests/equivalence_calibration_log.md`). `tests/equivalence_cells.py`'s own module docstring cross-references this rule; `TestBenchmarkRegistryTie` (in `tests/test_rpa_flex_equivalence_table.py`) enforces that this file's MOST RECENT section's cell_id column is an EXACT multiset match (duplicates detected, not masked by set semantics) against `CELLS` -- a cell added to the registry without a matching new benchmark row fails that test.

---

## Section 1 -- candidate pass on two development machines

- **Commit:** b922a1c13b85b2d319bc65ee8c45183dc6ab2a47 (the revision the measurements were taken at -- see `tests/equivalence_calibration_log.md`, Event 2, for the identical measured-source-commit convention and the full residual side of this same sweep).
- **Protocol:** `python -m tests.equivalence_measure` run 3x per runner (`find . -name __pycache__ -exec rm -rf {} +` before each invocation); per-cell/module timings below are the MAX over each runner's 3 invocations, then the MAX across the 2 runners (the "MAX over 3 measurement-process invocations per gating runner" rule, applied here to the 2 advisory development machines in place of the 4 CI gating runners; residuals were bitwise identical across all 3 invocations on both runners for every cell, so MAX-over-3 and any single invocation agree).
- **Runner matrix:** macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) + Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1, `PYTHONPATH=src`, `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1` -- see the calibration log's Event 2 for why the thread-count pin is required on that machine).

### Per-cell timing (37 rows -- every cell in `CELLS`)

| cell_id | max_seconds_over_runs_and_runners |
|---|---|
| `general.ring.onsite_coulombintra.fixedmu` | 0.343860 |
| `general.ring.onsite_coulombinter.fixedmu` | 0.011977 |
| `general.ring.onsite_hund.fixedmu` | 0.010564 |
| `general.ring.onsite_ising.fixedmu` | 0.011406 |
| `general.ring.onsite_exchange.fixedmu` | 0.011167 |
| `general.ring.onsite_pairhop.fixedmu` | 0.010949 |
| `general.ring.onsite_pairlift.fixedmu` | 0.009918 |
| `general.ring.onsite_u_v_hund.mu` | 0.013040 |
| `general.ring.onsite_full_kanamori.mu` | 0.013822 |
| `reduced.ring.onsite_coulombintra.spinfree.mu` | 0.006978 |
| `reduced.ring.onsite_coulombinter.spinfree.mu` | 0.007061 |
| `reduced.ring.onsite_hund.spinfree.mu` | 0.008863 |
| `reduced.ring.onsite_ising.spinfree.mu` | 0.007015 |
| `reduced.ring.onsite_pairlift.spinfree.mu` | 0.006716 |
| `reduced.ring.onsite_coulombintra.spindiag.mu` | 0.008119 |
| `reduced.ring.onsite_coulombinter.spindiag.mu` | 0.008412 |
| `reduced.ring.onsite_exchange.reject` | 0.001470 |
| `reduced.ring.onsite_pairhop.reject` | 0.001359 |
| `general.ring.offsite_coulombinter_sameorb.mu` | 0.006669 |
| `reduced.ring.offsite_coulombinter.mu` | 0.006491 |
| `general.ring.offsite_coulombinter_interorb.flexreject` | 0.008423 |
| `general.ring.offsite_hund.flexreject` | 0.007674 |
| `general.ring.offsite_ising.flexreject` | 0.008440 |
| `general.ring.offsite_exchange.flexreject` | 0.007598 |
| `general.ring.offsite_pairhop.flexreject` | 0.007679 |
| `general.ring.offsite_coulombintra_literalkey.reject` | 0.000717 |
| `general.ring.offsite_coulombinter_sameorb.subshape` | 0.005541 |
| `general.ring.onsite_coulombinter.subshape.mu` | 0.119363 |
| `general.ring.onsite_coulombinter.coefftail.mu` | 0.013468 |
| `auto.density.resolution` | 0.001703 |
| `auto.exchange.resolution` | 0.005980 |
| `auto.pairhop.resolution` | 0.001612 |
| `chi0q_init.reuse` | 0.012345 |
| `ringladder.general.onsite_coulombintra` | 0.001727 |
| `so.general.construction.reject` | 0.002737 |
| `so.reduced.construction.reject` | 0.002049 |
| `general.ring.offsite_coulombinter.conditioning.mu` | 0.016424 |

### Module totals and process overhead

| quantity | macOS arm64 | Linux x86_64 |
|---|---|---|
| `module_total_seconds_per_runner` (MAX of 3 invocations) | 0.557 | 0.638 |
| `process_seconds` (MAX of 3 invocations, external wall clock via `date +%s.%N`) | 0.781 | 0.819 |
| `process_overhead_seconds` (`process_seconds - module_total_seconds`, import/JSON-emission/interpreter startup-shutdown) | 0.224 | 0.181 |
| `unittest_module_process_seconds` (`python -m unittest tests.test_rpa_flex_equivalence_table`, 1 sample) | 0.853 | 0.962 |

**120s freeze-time budget:** both samples above are informational only here (this pass has 2 advisory proxies, not the 4 gating runners the budget is actually defined over) -- both are ~125x under the 120s budget regardless, so no headroom concern is expected once the real 4-sample MAX is measured in CI.

**Source SHA / run identity:** `b922a1c13b85b2d319bc65ee8c45183dc6ab2a47` on all 6 measurement invocations (3 per machine) and both unittest-timing samples, verified equal via `tests.equivalence_freeze_check`'s own source-SHA check applied to this pass's raw JSON output. No GitHub Actions run ID/attempt exists for this section -- both machines are proxy runs, not CI jobs; the CI section appended later records the authoritative `(run ID, attempt)` pair. **The rule for choosing WHICH `(run ID, attempt)` is authoritative when a calibration produces more than one is not yet written down anywhere.** Whoever appends the first CI section has to decide it and state it there.

---

## Section 2 -- freeze pass on the four gating continuous-integration runners

- **Commit:** 8144bf3f9e9539bad4759a2fbd1b24f52f7bef33 (the revision every measurement below was taken at -- see `tests/equivalence_calibration_log.md`, Event 3, for the residual side of this same sweep and for the freeze decision it supports).
- **Authoritative run identity:** `.github/workflows/equivalence-calibration.yml`, **workflow run 32204319966, attempt 1**, event `push`, head SHA `8144bf3f9e9539bad4759a2fbd1b24f52f7bef33`, completed successfully with all four jobs green. This was the FIRST AND ONLY iteration of the calibration loop: the measurements required no change to any tolerance, to any `Diverges` classification, or to the 120 s budget, so no re-measurement at a later source commit was needed and there are no superseded runs or attempts.
- **Rule for choosing the authoritative `(run ID, attempt)`** (Section 1 recorded that this rule was not written down; it is written down here, as the rule actually exercised): take runs of the calibration workflow whose head SHA is an EXACT match for the measured source commit, whose triggering event is `push`, and whose status is completed and conclusion successful; among those, take the LATEST completed successful attempt. Applying it here selected run 32204319966 attempt 1, which was also the only candidate.
- **Protocol:** `python -m tests.equivalence_measure` run 3x per gating runner (12 measurement samples) plus one timed `python -m unittest tests.test_rpa_flex_equivalence_table` process per runner (4 timing samples). Each job starts from a fresh checkout, so no `__pycache__` clearing step is needed; `OPENBLAS_NUM_THREADS`/`OMP_NUM_THREADS`/`MKL_NUM_THREADS`/`NUMEXPR_NUM_THREADS` are pinned to 1 by the workflow (Section 1's segfault finding; no crash occurred here on any of the 12 invocations). Per-cell timings below are the MAX over each runner's 3 invocations, then the MAX across the 4 runners -- the "MAX over 3 measurement-process invocations per gating runner" rule in its intended form, this time on the gating runners themselves rather than on advisory proxies.
- **Runner matrix:** all four jobs on `ubuntu-latest` (Linux-6.17.0-1022-azure-x86_64-with-glibc2.39): Python 3.9.25 (numpy 1.26.4, scipy 1.13.1), Python 3.10.21 (numpy 1.26.4, scipy 1.15.3), Python 3.11.16 (numpy 1.26.4, scipy 1.17.1), Python 3.12.13 (numpy 1.26.4, scipy 1.17.1).

### Per-cell timing (37 rows -- every cell in `CELLS`)

| cell_id | max_seconds_over_runs_and_runners |
|---|---|
| `general.ring.onsite_coulombintra.fixedmu` | 0.637397 |
| `general.ring.onsite_coulombinter.fixedmu` | 0.014503 |
| `general.ring.onsite_hund.fixedmu` | 0.012623 |
| `general.ring.onsite_ising.fixedmu` | 0.013951 |
| `general.ring.onsite_exchange.fixedmu` | 0.013634 |
| `general.ring.onsite_pairhop.fixedmu` | 0.013624 |
| `general.ring.onsite_pairlift.fixedmu` | 0.012230 |
| `general.ring.onsite_u_v_hund.mu` | 0.015882 |
| `general.ring.onsite_full_kanamori.mu` | 0.017074 |
| `reduced.ring.onsite_coulombintra.spinfree.mu` | 0.008333 |
| `reduced.ring.onsite_coulombinter.spinfree.mu` | 0.008121 |
| `reduced.ring.onsite_hund.spinfree.mu` | 0.008028 |
| `reduced.ring.onsite_ising.spinfree.mu` | 0.008141 |
| `reduced.ring.onsite_pairlift.spinfree.mu` | 0.008735 |
| `reduced.ring.onsite_coulombintra.spindiag.mu` | 0.010343 |
| `reduced.ring.onsite_coulombinter.spindiag.mu` | 0.010411 |
| `reduced.ring.onsite_exchange.reject` | 0.002052 |
| `reduced.ring.onsite_pairhop.reject` | 0.001620 |
| `general.ring.offsite_coulombinter_sameorb.mu` | 0.008270 |
| `reduced.ring.offsite_coulombinter.mu` | 0.007771 |
| `general.ring.offsite_coulombinter_interorb.flexreject` | 0.010640 |
| `general.ring.offsite_hund.flexreject` | 0.009817 |
| `general.ring.offsite_ising.flexreject` | 0.010766 |
| `general.ring.offsite_exchange.flexreject` | 0.009674 |
| `general.ring.offsite_pairhop.flexreject` | 0.011266 |
| `general.ring.offsite_coulombintra_literalkey.reject` | 0.000524 |
| `general.ring.offsite_coulombinter_sameorb.subshape` | 0.007587 |
| `general.ring.onsite_coulombinter.subshape.mu` | 0.182526 |
| `general.ring.onsite_coulombinter.coefftail.mu` | 0.021914 |
| `auto.density.resolution` | 0.002459 |
| `auto.exchange.resolution` | 0.002030 |
| `auto.pairhop.resolution` | 0.001694 |
| `chi0q_init.reuse` | 0.021667 |
| `ringladder.general.onsite_coulombintra` | 0.002732 |
| `so.general.construction.reject` | 0.003745 |
| `so.reduced.construction.reject` | 0.002805 |
| `general.ring.offsite_coulombinter.conditioning.mu` | 0.016737 |

### Module totals and process overhead

| quantity | Python 3.9 | Python 3.10 | Python 3.11 | Python 3.12 |
|---|---|---|---|---|
| `module_total_seconds_per_runner` (MAX of 3 invocations) | 0.668 | 0.840 | 0.966 | 0.986 |
| `process_seconds` (MAX of 3 invocations, external wall clock via `date +%s.%N`) | 0.926 | 1.094 | 1.549 | 1.277 |
| `process_overhead_seconds` (the two rows above differenced, i.e. import/JSON-emission/interpreter startup-shutdown) | 0.258 | 0.254 | 0.584 | 0.291 |
| `unittest_module_process_seconds` (`python -m unittest tests.test_rpa_flex_equivalence_table`, 1 sample) | 0.983 | 1.022 | 1.112 | 1.466 |

**120s freeze-time budget:** the gated quantity is the MAX of the four `unittest_module_process_seconds` samples above, which is **1.466 s** (Python 3.12) against the **120 s** budget -- roughly 82x of headroom. The budget holds, so it needed no change and the freeze was not blocked on timing.

**Conditioning-cell amplification:** the diagnostic's assertion 2 (the chemical-potential seam) amplifies by **1.07e+03x** on every one of the four gating runners -- identical across Python 3.9, 3.10, 3.11 and 3.12, and far above the >= 10x conditioning criterion. See `tests/equivalence_calibration_log.md`, Event 3, for the full checkpoint table.

**Source SHA / run identity:** `8144bf3f9e9539bad4759a2fbd1b24f52f7bef33` on all 12 measurement invocations, verified by `python3 -m tests.equivalence_freeze_check <artifacts> --source-sha 8144bf3f9e9539bad4759a2fbd1b24f52f7bef33`, which reports **VALIDATION OK (12 measurement + 4 unittest samples)**. Recomputing all 42 (cell, observable) bounds from the gating-runner MAX using the registry's own rounding rule reproduces every recorded atol exactly: **no bound changed**, no bound exceeds its policy ceiling, and no cell diverges or straddles its ceiling. On the strength of this section and Event 3, the registry's `PROVENANCE` record moves to `status = "frozen"`, naming this commit and this run.
