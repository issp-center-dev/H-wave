# IR general-scheme SDD progress

Branch: feat/ir-general-scheme
Plan: docs/superpowers/plans/2026-07-12-ir-general-scheme.md
Spec: docs/superpowers/specs/2026-07-12-ir-general-scheme-design.md

- Task 1: _calc_chi0q_general_ir (rank-6 IR bubble)         — complete (ab93469..467f2f9, review clean)
- Task 2: _calc_self_energy_general_ir                       — complete (467f2f9..1b1f769, review clean)
- Task 3: SCF dispatch + guard removal + VRAM preflight      — complete (1b1f769..5b248d1, review clean)
- Task 4: dynamic Eliashberg e2e + IR-native round-trip      — complete (no production fix needed; see task-4-report.md)
- Task 5: CPU/GPU parity + dtype + static handling           — pending

## Minor findings (for final review triage)
- T1: _calc_chi0q_general_ir has no scheme guard comment (low-level primitive; fine).
- T1: test _make_general_solver flips use_ir/matsubara_basis post-construction to bypass the guard; Task 3 removes the guard -> simplify to construct with matsubara_basis='ir' directly.
- T4: dynamic Eliashberg already consumes general+IR FLEX outputs (densified and IR-native) cleanly through the existing chi_convention="myo" tag; only test-side output wiring (explicit save_results + interaction-file plumbing) was missing. lam_ir=-0.4948 vs lam_u=-0.4944 (rel diff 0.08%), lam_native=-0.4954, all finite/consistent.
