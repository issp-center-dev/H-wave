# IR Stage 3: IR-native FLEX outputs (design addendum to ir-matsubara.md) — v4

## 1. Motivation (measured)

With `matsubara_basis="ir"` (Stage 2) the FLEX SCF costs ~0.1 s/iteration,
but a fixed cost dominates the run: densifying every output back onto the
uniform Nmat grid. Profiling (32x32, Nmat=2048, 15 iters, CPU): total 5.5 s
= densify 3.7 s (67%; of which 3.4 s is sparse-ir `basis.uhat()` evaluated
at all Nmat uniform frequencies and only 0.3 s the expansion GEMMs)
+ basis construction 1.1 s + SCF 0.7 s. On clavius (64x64, Nmat=4096, GPU)
the same structure leaves a ~10 s floor out of 12 s. Stage 3 removes the
densify + large-file cost by writing outputs natively on the sparse nodes.

## 2. Parameter

- `[mode.param] write_densified` (bool, default `true`).
  - `true` (default): Stage-2 behavior, byte-compatible files (OQ-1
    semantics preserved).
  - `false`: requires `matsubara_basis="ir"`; a `ValueError` otherwise
    (the option is meaningless on the uniform path — fail fast, do not
    warn-and-ignore).
  - With `false`, `solve()` skips `_ir_densify` entirely and
    `save_results()` writes the IR-native schema below.

## 3. IR-native npz schema

Every IR-native file stores node VALUES (not IR coefficients): values at
integer Matsubara indices are self-describing and independent of any
sparse-ir SVE reproducibility across versions (OQ-3 spirit). Coefficients
would require bit-reproducing the writer's basis to interpret.

**Discriminator (D-1)**: a file is IR-native iff `matsubara_basis == "ir"`,
`frequency_grid == "sparse_ir_nodes"`, AND `ir_freq_n` is present. IR
provenance alone is NOT sufficient: the Stage-1 `gap_dynamic.npz` already
carries `matsubara_basis="ir"` on a DENSIFIED array. A sparse-grid tag or
`ir_freq_n` without the other required discriminator fields is malformed and
readers reject it rather than falling through to a legacy uniform-grid path.

Common added metadata (all IR-native files):
- `matsubara_basis="ir"`.
- `frequency_grid="sparse_ir_nodes"`: explicit schema marker used by both
  in-repo readers and external tooling/migration checks.
- `ir_freq_n` (int64, 1-D): the node Matsubara indices n, with
  i*omega = n*pi/beta (odd fermionic / even bosonic), strictly increasing,
  symmetric under negation. One key per file; each file has exactly one
  frequency axis and `ir_statistics` states its parity.
- `ir_beta`, `ir_wmax`, `ir_tol` (floats), `ir_L` (int),
  `ir_statistics` ("F"/"B").
- `freq_index` is OMITTED (its uniform-grid positional semantics would
  lie); `nmat` is kept as run provenance but is NOT the frequency-axis
  length.

Per file:
- `chi0q.npz` / `chiq_s` / `chiq_c` / combined `chiq`: bosonic nodes,
  array shape (n_freq_B, nvol, nd, nd) [reduced]; keeps `wavevector_*`,
  `chi_convention` / `index_convention` unchanged.
- `sigma.npz` / `green.npz`: fermionic nodes, shape
  (nblock, n_freq_F, nvol, nd, nd); `sigma` keeps `cell_shape` (as today).
  `green` gains `cell_shape` in IR-native files ONLY (provenance; the
  densified `green.npz` key set stays unchanged — documented as a
  native-only key, not a legacy-compat key).

## 4. Consumer contracts (from the downstream-reader audit)

Full np.load consumer audit (repo-wide):

| Reader | Reads FLEX chi/sigma/green? | Action |
|---|---|---|
| `sc.py _load_chi0q` (static Eliashberg) | yes (chi0q) | fail fast on IR-native |
| `sc.py _read_flex_chi_raw` -> static `_load_flex_susceptibilities` | yes (chiq_s/c) | fail fast on IR-native (static machinery is uniform-only) |
| `sc.py _load_flex_susceptibilities_full` + `_load_flex_green` (dynamic) | yes | IR-aware: return arrays + IR meta (see below) |
| `RPA._read_chi0q` (`chi0q_init`) | yes (chi0q) | fail fast on IR-native |
| `flex._read_sigma` (`sigma_init`) | yes (sigma) | IR-aware (see below) |
| `dos.py`, `uhfk.py` green/eigen, `read_input_k` green_init | no (different file families) | no change; recorded here for completeness |
| `sample/RPA/view.py`, external scripts | RPA outputs / positional indexing | RPA unaffected; external risk mitigated by opt-in + D-1 + docs (R-S3-1) |

Rule: any reader that assumes uniform-grid semantics MUST apply D-1 FIRST
and raise with an actionable message. Mandatory because several readers
have permissive legacy fallbacks (`_static_freq_position` treats a missing
`freq_index` as "use the center row") that would otherwise silently
misread an IR-native file. The guard is a single shared helper
(`sc._reject_ir_native(data, file_name, hint)`), called before
`_read_freq_meta`/`_static_freq_position` in every uniform-only reader.

### 4.1 Dynamic Eliashberg fast path (the intended consumer)

- **Loader boundary (resolves OQ-S3-A)**: `sc.py` stays the canonical file
  reader but stays transform-free. `_read_flex_chi_raw`,
  `_load_flex_susceptibilities_full` and `_load_flex_green` gain an
  `allow_ir=False` keyword; when True they return `(arrays..., ir_meta)`
  where `ir_meta` is `None` for uniform files or a small dict
  (`freq_n`, `beta`, `wmax`, `tol`, `L`, `statistics`) per file family
  (chi: bosonic; green: fermionic). When False (default, used by every
  static caller) they apply D-1 and raise — the return ARITY changes only
  under `allow_ir=True`, so every existing call site keeps its current
  unpacking untouched; the two modes are covered by the fail-fast and
  chain tests respectively. All validation and refitting happens in
  `eliashberg_dynamic` (which owns the run's axes).
- **Preflight (must-fix M-1)**: `solve_dynamic` applies D-1 BEFORE the
  existing `_npz_freq_size`-vs-Nmat equality check; for IR-native inputs
  that uniform-grid check is skipped, and the memory guard
  (`check_memory`) is sized from the IR node count (post-refit axis
  length), not config Nmat. Mixed inputs (chiq_s native, chiq_c densified,
  etc.) are rejected — all three files must be the same encoding.
- `[eliashberg] matsubara_basis="uniform"` + IR-native input: fail fast
  (switch to `"ir"` or re-run FLEX with `write_densified = true`).
- `[eliashberg] matsubara_basis="ir"` + IR-native input: validate
  `np.isclose(ir_beta, run_beta, rtol=1e-9, atol=1e-9*run_beta)` (error
  message prints both values); warn if run `ir_wmax` < file `ir_wmax` or
  run `ir_tol` > file `ir_tol` (weaker target basis — refit quality is
  empirical, not guaranteed; see R-S3-3); then refit node values onto the
  RUN's axes via `IRAxis.fit_from_freq_points` (Sec. 5), replacing the
  full-grid `_ir_compress` (no `drop_constant` column: node values carry
  no uniform-FFT delta(tau) artifact). Refit residual at the file's own
  nodes is logged with the same >5%-scale warning as `_ir_compress`.
- **Ingestion primitive returns NODE VALUES, not coefficients** (round-2
  must-fix): the dynamic pipeline (`compute_vertices_flex_dynamic`,
  `calc_g2_dynamic`, the IR kernel) consumes values on the RUN axis'
  nodes. The refit is therefore the composition
  `values@file_nodes -> fit_from_freq_points -> coefficients ->
  eval_to_freq -> values@run_nodes`, wrapped in one helper
  `_ir_refit_nodes(arr, file_meta, run_axis, label)` in
  `eliashberg_dynamic`. When the file node set equals the run axis'
  `freq_n` exactly (`np.array_equal`), the helper RETURNS THE INPUT VALUES
  UNCHANGED (pure pass-through, no fit/eval round trip — a projection
  would not be zero-cost and could perturb values). This exact-equality
  case is the common one (same beta, auto wmax, same sparse-ir version).
- `chiq_s` and `chiq_c` each carry their own `ir_freq_n`/meta and are
  validated and refit INDEPENDENTLY (both must be bosonic; green
  fermionic). The same-encoding requirement (all native or all densified)
  is about encoding, not about node-set equality across files.

### 4.2 sigma_init warm start (sweep chains, e.g. bp-ICl2 lambda(P,T))

- `flex._read_sigma` changes to return the pair `(sigma_array, ir_meta)`
  with `ir_meta=None` for uniform files (it is private with exactly one
  call site); `read_init()` unpacks it and stores
  `info["sigma_init"]` (plain ndarray in both encodings) and
  `info["sigma_init_ir"]` (the meta dict or absent). This keeps every
  downstream consumer of `green_info["sigma_init"]` type-stable
  (resolves OQ-S3-B).
- In `solve()`: an IR-native seed requires `matsubara_basis="ir"` for this
  run (fail fast otherwise: "re-run the seeding FLEX with
  write_densified=true to seed a uniform run" — v1 scope, revisit only if
  a real workflow needs it; recorded as OQ-S3-D disposition). The seed is
  brought to the run's nodes via the WRITER's basis: reconstruct
  `IRAxis(file ir_beta/ir_wmax/ir_tol, "F")`, `fit_from_freq_points` at
  the file's own nodes (determined by construction), then
  `eval_to_freq_points` at the RUN's node indices. A direct fit onto the
  run basis would be UNDERDETERMINED on every high-T -> low-T sweep step
  (the run's L grows with beta while the file carries only ~L_file nodes)
  — found in implementation, recorded as OQ-S3-J. The existing
  `sigma_init_on_error` residual policy applies unchanged (the residual is
  evaluated at the file's nodes, in the writer basis). Shape check
  compares against the file's node count, mirroring `expected_uniform`
  today.
- **Cross-temperature seeding semantics (deliberate, OQ-S3-G)**: sweeps
  seed from a NEIGHBOURING temperature, so the seed's `ir_beta` will
  usually differ from the run's beta. This is NOT an error — it mirrors
  the uniform path today, where a different-T sigma.npz is consumed
  index-matched on the Nmat grid (frequencies implicitly rescaled by the
  beta ratio; a warm start needs a good starting point, not exactness).
  The IR-native equivalent: the file's integer `ir_freq_n` are interpreted
  as node indices ON THE RUN'S beta (same index-matched rescaling), then
  refit. When `ir_beta` differs from the run beta, an INFO line logs both
  values and the ratio; the `sigma_init_on_error` residual policy is the
  quality guard. (Contrast Sec. 4.1: for the dynamic Eliashberg the chi
  files are PHYSICS INPUT, not a warm start, so there beta mismatch is a
  hard error.)
- Uniform seed into an IR run: unchanged Stage-2 behavior
  (`_ir_sigma_init`).

### 4.3 Canonical fail-fast messages

Every guard names the file, the fact ("holds sparse-IR node data"), and
the exact TOML fix(es) with their section:

| Case | Canonical message (suffix after "file '<name>' holds sparse-IR node data (frequency_grid=sparse_ir_nodes); ...") |
|---|---|
| static hwave_sc chi0q / static FLEX-chi path | "the static Eliashberg solver requires uniform-grid files. Re-run FLEX with `[mode.param] write_densified = true`, or switch to `frequency = \"dynamic\"` with `[eliashberg] matsubara_basis = \"ir\"`." |
| `RPA chi0q_init` | "chi0q_init requires a uniform-grid file. Re-run FLEX with `[mode.param] write_densified = true`." |
| dynamic Eliashberg, `[eliashberg] matsubara_basis = "uniform"` | "set `[eliashberg] matsubara_basis = \"ir\"` to consume it, or re-run FLEX with `[mode.param] write_densified = true`." |
| `sigma_init` into a uniform FLEX run | "an IR-native seed requires `[mode.param] matsubara_basis = \"ir\"` in this run, or re-run the seeding FLEX with `[mode.param] write_densified = true`." |
| mixed encodings among chiq_s/chiq_c/green | "all dynamic-Eliashberg inputs must share one encoding; re-run FLEX so chiq_s/chiq_c/green are all densified or all IR-native." |
| `write_densified=false` without IR | "`[mode.param] write_densified = false` requires `[mode.param] matsubara_basis = \"ir\"`." |
| beta mismatch (dynamic chi/green input) | "ir_beta {file} differs from this run's beta {run} (rel {x}); the susceptibilities are physics input and must match. Re-run FLEX at this temperature." |

## 5. `IRAxis.fit_from_freq_points(arr, freq_n)`

Mirror of `eval_to_tau_points`: least-squares fit of values given at
ARBITRARY integer Matsubara indices `freq_n` onto this axis' coefficients.
Like every `IRAxis` transform it CONTRACTS THE LAST AXIS (callers move the
frequency axis to the last position first — the existing `_ir_compress`
call pattern). Matrix orientation follows the existing convention:
`design = basis.uhat(freq_n).T` with shape (npts, L), fit matrix
`np.linalg.pinv(design).T` with shape (npts, L) applied as `arr @ fit`
(values (..., npts) -> coefficients (..., L)).
Validation (must-fix M-3; all raise `ValueError`):
- `freq_n` 1-D, integer, strictly increasing (implies unique);
- parity matches the axis (`all odd` for F, `all even` for B);
- `len(freq_n) >= L`.
The design-matrix rank (explicit `numpy.linalg.matrix_rank`) must equal
`L`, else raise (a rank-deficient node set cannot determine the
coefficients). The design-matrix condition number is logged at DEBUG; a
WARNING is emitted above 1e10 (ill-conditioned fits can amplify noise into
lambda shifts without any file-format error — R-S3-4). Matrices are cached
per (`freq_n.tobytes()` hash, backend), like the tau-point cache.

## 6. Non-goals (v1)

- `Nmat` remains a required parameter (defines the output grid when
  densified, the uniform sigma_init seed grid, and the dynamic-Eliashberg
  gap output grid). Making it optional under pure-IR chains is a
  follow-up.
- `gap_dynamic.npz` stays densified (cost is small; can adopt the schema
  later — its provenance-only `matsubara_basis` key is why D-1 requires
  `ir_freq_n` too).
- The `general` calc_scheme and RPA solver are untouched.
- IR-native sigma_init into a uniform FLEX run: unsupported in v1
  (fail fast with the regenerate hint).

## 6b. Documentation checklist (part of the deliverable, not an afterthought)

FLEX tutorial (EN/JA) must state: when to keep `write_densified = true`
(default; any legacy/external analysis, static hwave_sc); when `false` is
safe (pure IR chains: IR FLEX -> dynamic-IR Eliashberg, IR->IR sigma_init
sweeps) and what it buys (removes the fixed densify+write cost, ~10 s and
GB-scale files at production sizes); how to RECOGNIZE a native file
(`frequency_grid="sparse_ir_nodes"` in the npz keys); and the recovery
path (re-run FLEX densified — cheap because `sigma_init` can seed from
the native sigma — or the documented 5-line python snippet that densifies
a native file offline via `ir_axis`; no new CLI utility in v1). Eliashberg
tutorial gains the matching one-paragraph note. A warning box tells users
NOT to point legacy scripts at native outputs.

## 7. Test plan

- Round trip: IR run with `write_densified=false` -> files carry the
  schema (D-1 keys, freq_n negation symmetry, no freq_index); densifying
  the stored node values offline reproduces the `write_densified=true`
  outputs within ~ir_tol.
- End-to-end chain: FLEX(ir, native) -> dynamic Eliashberg(ir) lambda
  equals the FLEX(ir, densified) -> dynamic Eliashberg(ir) lambda within
  eps-scale (the refit replaces the lossy full-grid compress, so agreement
  should be TIGHTER than the Stage-1 uniform gate).
- sigma_init chain: run A (ir, native sigma) seeds run B (ir); B converges
  and matches a B seeded from A's densified sigma.
- Fail-fast matrix: static hwave_sc chi0q path, static FLEX-chi path, RPA
  chi0q_init, uniform-run sigma_init, dynamic-uniform loader, and mixed
  native/densified chi inputs each raise the documented message;
  `write_densified=false` without `matsubara_basis="ir"` raises.
- `fit_from_freq_points`: identity when `freq_n` == the axis' own nodes
  (degenerates to `fit_from_freq`); wrong-parity, non-monotonic,
  duplicate, and underdetermined (len < L) inputs raise; refit across a
  different node set (file basis with different eps) matches direct
  evaluation within the looser eps; rank-deficient set raises.
- Memory guard: IR-native dynamic load sizes `check_memory` from node
  count (assert via a large-Nmat config that would fail the uniform-sized
  guard but passes IR-sized).
- Migration: pre-Stage-3 densified files load unchanged through every
  reader (pin by re-running an existing uniform fixture through the
  updated loaders); uniform-only readers reject native files BEFORE any
  legacy fallback (`_static_freq_position` center-row path unreachable on
  native input); native sigma from run A seeds run B without manual file
  handling (one TOML line).
- Cross-temperature seed: run A at T1 seeds run B at T2 (ir->ir);
  INFO log records both betas; B converges (semantics of OQ-S3-G).

## 8. Risks

- R-S3-1 (silent misread by EXTERNAL analysis scripts that index the
  frequency axis positionally, e.g. `shape[0]//2` static slice): mitigated
  by default `write_densified=true` (opt-in only), D-1, and docs; in-repo
  readers all gain the mandatory guard (Sec. 4).
- R-S3-2 (npz key-set change breaking strict loaders): only NEW opt-in
  files carry new keys; default output stays byte-compatible with Stage 2
  (including `green.npz` NOT gaining `cell_shape` when densified).
- R-S3-3 (refit quality between differing bases): bounded empirically, not
  by theory; guarded by the weaker-basis warning (wmax/eps comparison),
  the logged refit residual, and the >5%-scale warning.
- R-S3-4 (hidden numerical change via ill-conditioned refit): guarded by
  rank check (raise) + condition-number warning in Sec. 5; the
  same-node-set fast path keeps the common case exactly equivalent to
  Stage 2's `fit_from_freq`.

## 9. Resolution record (open questions from review round 1)

- OQ-S3-A (who refits): RESOLVED — sc.py reads and returns raw arrays +
  meta behind `allow_ir=True`; eliashberg_dynamic owns validation and
  refit (it owns the run axes). Static callers keep `allow_ir=False` and
  fail fast.
- OQ-S3-B (`_read_sigma` return type): RESOLVED — the array stays a plain
  ndarray under the existing key; IR meta travels under the new sibling
  key `sigma_init_ir` in the info dict. No return-type fork.
- OQ-S3-C (`ir_freq_n` naming): RESOLVED — one common key; each file has
  exactly one frequency axis, and `ir_statistics` + parity validation in
  `fit_from_freq_points` prevent cross-parity misuse mechanically.
- OQ-S3-D (IR-native sigma_init into uniform runs): RESOLVED — v1
  unsupported with fail-fast + regenerate hint; revisit only on a concrete
  workflow need (sweeps chain ir->ir, the fast case, or uniform->uniform).
- OQ-S3-E (general/full-vertex + IR-native output, round 2): RESOLVED by
  construction — `write_densified=false` requires `matsubara_basis="ir"`,
  and Stage 2 already rejects `matsubara_basis="ir"` with
  `calc_scheme="general"` at `__init__`; native general output therefore
  cannot arise. No additional guard needed (covered by the existing
  init-time test).
- OQ-S3-F (writer-side `ir_freq_n` symmetry validation, round 2):
  RESOLVED — the writer stores the run axis' `freq_n`, whose negation
  symmetry is already ASSERTED at `IRAxis` construction; the writer
  inherits that guarantee (no re-validation). Readers validate only what
  they rely on (parity, monotonicity, count/rank) because refit onto the
  run axis does not require the file set to be symmetric.
- Round-2 must-fix (fast path returned coefficients): RESOLVED in Sec. 4.1
  — the ingestion helper `_ir_refit_nodes` returns run-axis NODE VALUES,
  with exact node-set equality short-circuiting to a pass-through of the
  stored values (no projection).
- OQ-S3-G (cross-temperature sigma_init seeding, user-lens round):
  RESOLVED — supported deliberately with index-matched rescaling semantics
  (Sec. 4.2), mirroring the uniform warm start; INFO log + residual policy
  as guards; dynamic-Eliashberg chi inputs remain a hard beta error
  (physics input, not warm start). Canonical messages in Sec. 4.3.
- OQ-S3-H (write_densified granularity, user-lens round): RESOLVED — v1 is
  one global flag; per-file granularity (e.g. native sigma + densified
  chi) deferred until a concrete workflow needs it. The recovery cost is
  low because a densified re-run can seed from the native sigma.
- OQ-S3-I (native->densified converter utility, user-lens round):
  RESOLVED — no new CLI surface in v1; the docs carry a 5-line python
  snippet (load npz -> IRAxis(file params) -> fit_from_freq_points ->
  eval_to_uniform -> save), and re-running FLEX densified stays the
  supported path. Revisit if users ask.
- OQ-S3-J (implementation round): a run-basis fit of the sigma_init seed
  is underdetermined for every high-T -> low-T sweep step (run L grows
  with beta; the file has only ~L_file nodes) — exactly the intended use
  case. RESOLVED by fitting in the reconstructed WRITER basis at the
  file's own nodes and evaluating at the run's node indices; this also
  realizes the index-matched semantics literally. Note the deliberate
  contrast: sigma_init MAY reconstruct the writer basis because a warm
  start needs only a good starting point (sparse-ir version drift merely
  perturbs the seed, and the residual policy guards it); the chi/green
  ingestion path must NOT depend on basis reconstruction (physics input,
  OQ-3), which is why `_ir_refit_nodes` fits file node VALUES onto the
  run's own basis instead.
