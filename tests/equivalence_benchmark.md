# RPA/FLEX equivalence table -- benchmark

Section schema: header (commit, protocol, runner matrix), then a table `cell_id | max_seconds_over_runs_and_runners` covering EVERY cell (`tests.equivalence_measure` times all proof shapes -- rejections, construction-only rows, and the two multirun exceptions included), then `module_total_seconds_per_runner`, `process_overhead_seconds`, the `unittest_module_process_seconds` sample(s), and the MEASURED SOURCE SHA + the authoritative run identity. **This file grows one section per calibration stage**: Section 1 below is the CANDIDATE pass taken on two development machines (advisory proxy runners, 3 measurement invocations each). A continuous-integration section is appended later (the 4 gating runners x 3 invocations + the 4 unittest-timing samples the 120s budget actually gates, with real GitHub Actions run IDs/attempts) once a closed calibration loop produces a freezable measurement.

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
