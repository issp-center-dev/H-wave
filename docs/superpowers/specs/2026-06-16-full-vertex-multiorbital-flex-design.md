# Full-vertex multi-orbital FLEX (`calc_scheme="general"`)

**Date:** 2026-06-16 (rev. 2026-06-19 after Codex design review)
**Branch:** `feature/full-vertex-flex`
**Status:** DESIGN — Codex review round 1 incorporated (2 critical + majors
addressed: V_eff gate §4.5, dressed-`G` transverse §4.3, S/C builder location
§4.3, index convention §4.4, spin-mode split §4.6, SOC guard §4.1, validation
equalities §6, memory budget §8). **One open decision for the user: confirm the
V_eff reference (§3/§9).** Awaiting user spec review before writing the plan.

## 1. Problem

H-wave's FLEX solver only runs in the density-density approximation. The whole
self-energy pipeline (`_inflate_chi0q_and_ham` → `_build_spin_charge_vertices` →
`_solve_rpa` → `_calc_veff` → `_calc_self_energy`) operates on the reduced
`(nvol, nd, nd)` interaction obtained by the diagonal contraction `kaabb->kab`
(`flex.py:384,401,409`). This **drops the orbital off-diagonal vertices** of
Exchange, Hund, PairLift and PairHop. `_init_flex_param` (`flex.py:84-115`)
therefore (a) rejects `calc_scheme="general"`/`ring+ladder` outright and (b)
warns whenever Exchange/PairHop are present that they are reduced to their
density-density part.

The density-density treatment is a recognised simplification but is **not** what
full multi-orbital FLEX requires: the transverse (spin-flip) Hund term and
pair-hopping affect spin fluctuations, the spin/charge vertex, and ultimately
Tc. We want full-vertex FLEX available as an option.

## 2. Goal / scope

Implement **full-vertex multi-orbital FLEX** including both the longitudinal
(spin/charge) channel with full rank-4 `(nd², nd²)` vertices **and** the
transverse (spin-flip, ladder-equivalent) channel, selected via
`calc_scheme="general"` (consistent with RPA, where `general` already drives the
`ring+ladder` full-tensor path).

Decisions locked in (from brainstorming):

- **Physics scope:** B — longitudinal **and** transverse full vertex.
- **Selection:** reuse `calc_scheme="general"`; no new FLEX-specific flag.
- **Definition of done (validation):** primary = a brute-force, no-approximation
  reference that matches the optimized self-energy functional in a single shot;
  secondary = degenerate-reduction + symmetry/analytic checks. See §6.

### Out of scope (explicit)

- **SOC / `enable_spin_orbital` with general-FLEX.** The transverse channel
  extraction assumes Sz conservation (`rpa.py:1982-1990` warns for genuine spin
  mixing). General-FLEX + SOC is rejected by a fail-fast guard in
  `_init_flex_param` (§4.1) until the transverse extractor is generalized.
- **`spinful` spin mode in general-FLEX (v1).** v1 supports **`spin-free` and
  `spin-diag`** only (§4.6). `spinful` (and SOC, above) are deferred: their
  transverse-bubble construction differs materially (`rpa.py:1974-2002`) and
  needs its own validation.
- **External numerical reference data.** None available now; the brute-force
  reference (§6) is the implementation-correctness anchor. External
  cross-checks are a later, separate task.
- **Reduced/squashed (density-density) FLEX numerics.** The reduced path's
  numerics and outputs stay unchanged and are pinned by a regression test (§7);
  the shared control flow (`_init_flex_param`, `solve()` dispatch) is edited, so
  the path is regression-guarded rather than literally untouched.

## 3. Reference formulation (what we cross-check the formulas against)

- **S/C interaction matrices (full Kanamori):** the multi-orbital matrix RPA of
  Kuroki et al. (5-band RPA + Eliashberg, [arXiv:0902.3691]) and Takimoto et al.
  (PRB 69, 104504) use the **full Kanamori vertex** (`U, U', J, J'`) in
  `(norb², norb²)` matrix form — **not** density-density. H-wave already builds
  exactly these S (spin) / C (charge) matrices in `sc.py`
  `_build_sc_matrices_all_q` (`sc.py:422-523`), reused by the SC pairing vertex
  `_compute_vertices_general`. **Caveat (per Codex review):**
  `_compute_vertices_general` is a *pairing* vertex
  (`V^s = 3/2 S χ_s S − 1/2 C χ_c C + 1/2 (S + C)`), **not** a normal-state
  self-energy kernel. It validates the S/C matrices and χ construction, but it
  **cannot** independently fix the FLEX self-energy coefficients.
- **Normal-state FLEX self-energy `V_eff` — RESOLVED, dual-sourced.** The exact
  `V_eff` is transcribed in §4.5 from **Mochizuki, Yanase, Ogata
  [cond-mat/0407094]** (user-confirmed primary) and **corroborated** by
  **Takimoto, Hotta, Ueda [cond-mat/0309575] = PRB 69, 104504 (2004)**. Both
  give the **same coefficient structure**:
  `(3/2) Ûˢ χˢ Ûˢ + (1/2) Ûᶜ χᶜ Ûᶜ − (1/4)(Ûˢ+Ûᶜ) χ⁰ (Ûˢ+Ûᶜ) + (3/2)Ûˢ − (1/2)Ûᶜ`
  (THU writes the charge channel as an orbital channel `Û°/χ°` and the
  subtraction with the cross-spin bare bubble `χ̄^(σσ̄)`, but the coefficients
  `3/2, 1/2, −1/4` and the `+(3/2)Ûˢ−(1/2)Û°` constants match). The gate is
  satisfied.
- **Both canonical references are paramagnetic / SU(2)-symmetric.** Neither
  carries a *separate spin transverse (`χ_+-`) channel*: the spin-flip physics
  is bundled into the `3/2 χˢ` coefficient, and neither does spin-dependent
  (`G↑≠G↓`) calculations. THU *does* split the **orbital** channel into
  longitudinal/transverse (`χ°_C/L/T`), but that is orbital, not spin. This
  **reframes the original "A vs B" choice** (§2, §4.6): standard multi-orbital
  FLEX = the paramagnetic full-vertex form; a separate spin-flip self-energy
  (the literal "B") is **not** part of these references and is treated as a
  deferred future extension needing its own reference.
- **S/C interaction matrices (full Kanamori):** Kuroki et al.
  ([arXiv:0902.3691]) and THU use the **full Kanamori vertex** (`U, U', J, J'`)
  in `(norb², norb²)` matrix form — **not** density-density. H-wave already
  builds these S/C matrices in `sc.py:_build_sc_matrices_all_q` (`sc.py:422-523`)
  for the SC pairing vertex. **Caveat:** `_compute_vertices_general` is a
  *pairing* vertex, **not** a self-energy kernel — it validates the S/C matrices
  and χ construction (cross-check against MYO/THU Eq.(10)) but does not by itself
  fix the FLEX self-energy coefficients (those come from §4.5).

Caveat recorded honestly: RPA/FLEX is itself an approximation. It resums the
RPA particle-hole bubble+ladder series but **omits the Aslamazov-Larkin (AL)
and Maki-Thompson (MT) vertex corrections** (two fluctuation propagators joined
by a triangular fermion loop / mode-mode coupling). "Full vertex" here means
*not density-density reduced*, **not** *exact* and **not** AL/MT-inclusive. This
boundary is now documented in the FLEX manual (`docs/{en,ja}/source/flex/
tutorial/tu-index.rst`, "Scope of the approximation"). **Future-work idea
(user, 2026-06-20):** after FLEX converges, optionally evaluate AL and MT once
(a non-self-consistent **1-shot diagnostic**) to gauge whether they are
significant for a given system — recorded in §10.

## 4. Architecture

Keep the existing reduced path untouched; add a **parallel general path**
dispatched from `solve()` on `self.calc_scheme == "general"`.

### 4.1 Dispatch / guard change (`_init_flex_param`)

- Replace the `general`-rejection branch: `general` → enable full-vertex FLEX;
  `reduced`/`squashed` → existing density-density path (unchanged default).
- Move the "exchange/pairhop reduced to density-density" warning so it fires
  **only** for `reduced`/`squashed` (under `general` nothing is dropped).
- Require the general-scheme rank-4 `chi0q` (already produced by RPA in
  `general`); fail loudly if a reduced `chi0q` is supplied in general-FLEX.
- **Fail-fast guard:** reject `calc_scheme="general"` when
  `ham_info.enable_spin_orbital` is true, and when `spin_mode == "spinful"`
  (both deferred from v1, §2). Clear error pointing to the supported modes.

### 4.2 New general methods (mirroring the reduced ones)

| Reduced (keep) | New general (rank-4) |
|---|---|
| `_inflate_chi0q_and_ham` | `_inflate_chi0q_and_ham_general` |
| `_build_spin_charge_vertices` | `_build_spin_charge_vertices_general` |
| `_calc_veff` | `_calc_veff_general` |
| `_calc_self_energy` | `_calc_self_energy_general` |
| (n/a) | `_build_transverse_channel_flex` → χ_+- via `_solve_rpa` (see §4.3) |

Each unit has one clear purpose, a shape-typed interface, and a degenerate
limit that reproduces the reduced method (testable in isolation — see §7).

### 4.3 Reuse from RPA / SC — with the corrections from the Codex review

- **S/C interaction matrices.** There is **no** named S/C builder in `rpa.py`;
  the general-scheme inflation there is **inline inside `solve()`**
  (`rpa.py:875-999`). The actual reusable full `(norb², norb²)` S/C builder is
  **`_build_sc_matrices_all_q` in `sc.py` (`sc.py:422-523`)**. Plan: **extract a
  shared helper** (e.g. to a common module) and call it from FLEX, RPA-general,
  and SC, rather than claiming a non-existent `rpa.py` reuse.
- `_solve_rpa` already inverts general-scheme matrices (flattens rank-4 to
  `(ndx, ndx)`, `rpa.py:2068-2138`); verify on rank-4 χ for both channels.
- **Transverse channel — must NOT reuse `_build_transverse_channel` as-is**
  (Codex critical). In `spin-diag` mode it builds `chi0_+-` from the **bare**
  `self.green0`/`green0_tail` (`rpa.py:1938-1956`), but FLEX needs the
  **dressed** `green_kw` recomputed each SCF iteration. Plan: add
  `_build_transverse_channel_flex` that mirrors the RPA algebra
  (`W_+-` from the per-interaction sign table `rpa.py:2011-2018` is reused — it
  depends only on the interaction, not on `G`) but takes the **current dressed
  `G`** as an explicit argument. The `spin-free` branch (paramagnetic,
  `chi0_+- = chi0_orb`, `rpa.py:1921-1935`) needs no `G` but must still consume
  the iteration's `chi0`.

### 4.4 Core change: the self-energy contraction

The reduced self-energy is an **element-wise (Hadamard)** product
`Σ_{ab}(r,τ) = V_eff,{ab}(r,τ) · G_{ab}(r,τ)` (`flex.py:642-645`). The general
self-energy is the **rank-4 contraction**

```
Σ_{l1 l2}(r,τ) = Σ_{l3 l4} V_{l1 l3 l2 l4}(r,τ) · G_{l3 l4}(r,τ)
```

i.e. `V_eff` is `(nd², nd²)` and contracts the internal orbital pair `(l3,l4)`
of `G` at each `(r,τ)`. **Index convention (fixed now):** orbital pairs are
flattened as `row = (l1,l2)`, `col = (l3,l4)` with `flat = l_outer*nd + l_inner`
(C/row-major), matching the `(ndx,ndx)` flattening `_solve_rpa` already uses for
χ and the S/C matrices (`rpa.py:2068-2138`). Worked example (norb=2, nd=2):
`V` is 4×4 with rows/cols ordered `(00,01,10,11)`; the density-density limit is
`V` nonzero only on entries `(ab),(ab)` (Codex finding A1.1). The brute-force
check (§6) confirms this map end-to-end. The FFT transport (Matsubara↔τ, k↔r) of
`_calc_self_energy` is reused unchanged; only the per-`(r,τ)` product becomes a
batched matmul over the flattened orbital indices.

### 4.5 `V_eff` assembly (`_calc_veff_general`) — formula pinned to MYO

Reference (user-confirmed): **Mochizuki, Yanase, Ogata, [cond-mat/0407094]**,
multi-orbital FLEX (paramagnetic / Sz-conserving). Transcribed verbatim
(via ar5iv) — interaction matrices `Ûˢ` (spin) and `Ûᶜ` (charge) from Eq.(10)
in the orbital-pair basis (intra `U`, inter `U'`, Hund `J_H`, pair-hopping
`J' = J_H`):

```
χ⁰_{mn,μν}(q) = −(T/N) Σ_k G_{μm}(k+q) G_{nν}(k)            (Eq.5)

χˢ(q) = [I − χ⁰(q) Ûˢ]⁻¹ χ⁰(q)                              (Eq.4)
χᶜ(q) = [I + χ⁰(q) Ûᶜ]⁻¹ χ⁰(q)

V(q) =  (3/2) Ûˢ χˢ(q) Ûˢ                                    (Eq.4)
      + (1/2) Ûᶜ χᶜ(q) Ûᶜ
      − (1/4)(Ûˢ+Ûᶜ) χ⁰(q) (Ûˢ+Ûᶜ)        ← double-counting subtraction
      + (3/2) Ûˢ − (1/2) Ûᶜ                ← first-order (Hartree-Fock) constants

Σ_{mn}(k) = (T/N) Σ_q Σ_{μν} V_{μm,νn}(q) G_{μν}(k−q)        (Eq.3)
```

**Key consequence — there is NO separately-computed transverse channel.** In
this paramagnetic matrix FLEX the spin fluctuation enters as the single matrix
`χˢ` with coefficient `3/2` (= 1 longitudinal Sz + 2 transverse S±, equal by
SU(2)); the transverse (spin-flip) physics is **already bundled** in `3/2 Ûˢ χˢ
Ûˢ`. A separate `χ_+-` ladder is only a distinct object when SU(2) is broken
(spin-polarised / `spin-diag` with `G↑≠G↓`), which MYO does **not** cover. This
**reframes the earlier "A vs B" decision** (see §2, §4.6): for the paramagnetic
`spin-free` case the full-vertex self-energy IS the MYO formula and needs no
ladder term. (Codex A1.2/A1.3 resolved: coefficients and the
`−(1/4)(Ûˢ+Ûᶜ)χ⁰(Ûˢ+Ûᶜ)` subtraction are now explicit, not deferred.)

Implementation notes / cross-checks:
- `Ûˢ, Ûᶜ` must equal H-wave's existing `sc.py:_build_sc_matrices_all_q` S/C
  matrices (the shared helper, §4.3) — verify element-wise against MYO Eq.(10).
- Sign mapping to `_solve_rpa` (`[1 + χ⁰·ham]⁻¹χ⁰`): `ham_s = −Ûˢ`, `ham_c = +Ûᶜ`
  (same convention as the reduced `_build_spin_charge_vertices`).
- The `+(3/2)Ûˢ − (1/2)Ûᶜ` constants are the first-order term; the reduced
  `_calc_veff` omits them (`flex.py:491-532`). Reconcile how H-wave handles
  first-order (Hartree) in FLEX — likely via the static term / `μ`, not in
  `V_eff` — and keep `_calc_veff_general` consistent with the reduced treatment
  to preserve the degenerate-reduction check (§6.1). Single-orbital identity:
  MYO `= reduced + U` (the constant), confirming the fluctuation part matches.
- Self-energy index ordering follows Eq.(3): `V_{μm,νn}` contracted with
  `G_{μν}` — this fixes §4.4's flatten map (`μ,ν` are the contracted/"outer"
  pair).

### 4.6 Spin-mode coverage (v1: `spin-free`, `spin-diag`)

The general path is **not** one uniform branch — χ0 inflation and the transverse
bubble differ by spin mode (`flex.py:369-411`, `rpa.py:1919-2002`):

- **`spin-free`** (paramagnetic): `chi0_+- = chi0_orb` (no `G` needed);
  self-energy expands orbital→spin-orbital like the reduced `nd_block != nd_v`
  case (`flex.py:618-665`).
- **`spin-diag`**: `chi0_+-` is built from the **dressed** spin-resolved `G`
  each iteration (§4.3); self-energy tensor is full `(nd, nd)` per spin block.
- **`spinful` / SOC**: deferred (§2), guarded (§4.1).

Each supported mode states its `chi0_+-` source and self-energy tensor shape and
gets its own validation target (§7).

## 5. Data flow (general path, one SCF iteration)

```
H0, Σ_in ─▶ dressed G(k,iω)              [_calc_dressed_green, reused]
        ─▶ χ0(q,iν) rank-4              [_calc_chi0q general, RPA]
        ─▶ inflate χ0, build S,C        [_inflate_*_general, _build_*_general]
        ─▶ χ_s, χ_c (longitudinal)      [_solve_rpa, general]
        ─▶ χ0_+-(dressed G), W_+- ─▶ χ_+-  [_build_transverse_channel_flex + _solve_rpa]
        ─▶ V_eff(q,iν) rank-4           [_calc_veff_general]
        ─▶ Σ_out(k,iω) = ΣV·G           [_calc_self_energy_general, rank-4 contraction]
        ─▶ mix, check convergence       [reused]
```

## 6. Validation = definition of done

**Primary — no-approximation 1-shot match.** A standalone brute-force reference
computes one evaluation of the self-energy functional `Σ[G]` by **direct
summation**, with no FFT convolution and no density-density reduction:

- `χ0_{l1l2l3l4}(q,iν) = −T/N Σ_k G_{l3l1}(k,iω) G_{l2l4}(k+q, iω+iν)` (explicit
  k-sum),
- `χ_s, χ_c, χ_+-` by matrix RPA,
- `Σ_{l1l2}(k,iω) = T/N Σ_{q,iν} Σ_{l3l4} V_{l1l3l2l4}(q,iν) G_{l3l4}(k−q, iω−iν)`
  (explicit q,iν double sum).

The optimized general path (rank-4 FFT convolution) must match this reference to
~1e-10 on a small system (e.g. 2 orbitals, 4×4 k-grid, few Matsubara
frequencies). One shot suffices: self-consistency only iterates `Σ[G]`, so a
correct functional is the core guarantee. "Slow is OK" — the reference is
test-only. The brute-force is written from the FLEX **diagrammatic** definition,
independently of the optimized `V_eff`, so a shared formula error is not masked.

**Scope of the primary check (Codex A2.1, stated honestly):** the 1-shot match
proves **algorithmic equivalence** — the optimized path computes the *same*
functional as the direct-sum reference. It does **not** by itself prove the
chosen FLEX *formula* is physically correct, because both encode the same
diagrammatic content. It is conclusive only **after** §3's reference formula is
transcribed and the brute-force is written from the diagrammatic definition. The
secondary checks below guard the formula choice.

**Secondary — limits & symmetry (guards the formula choice), as precise
equalities (Codex A2.2):**

1. **Degenerate reduction:** general-FLEX with only `U, U'` (no Hund/PairHop/
   Exchange) reproduces reduced-FLEX to ~1e-10 for each of
   `χ0, χ_s, χ_c, Σ, energy` (reduced path keeps `CoulombInter`, `flex.py:384-401`).
2. **Single-orbital SU(2) limit:** `χ_+-(q,iν) == χ_s(q,iν)` to ~1e-10, and
   general-FLEX `Σ`/`energy` == reduced-FLEX (no off-diagonal exists). (Replaces
   the vague "transverse == longitudinal" wording.)
3. **RPA consistency:** first-iteration `χ_s`, `χ_c` from general-FLEX ==
   RPA general-scheme `χ` for the same `χ0`/interaction (~1e-10).
4. **SU(2) vertex identity:** for full Kanamori, `W_+- == W_zz` element-wise
   (`rpa.py:2011-2018`); `χ_s`, `χ_c`, `χ_+-` Hermitian; physical-sign checks.
5. **Dressed-G transverse (Codex A2.3):** a 2-iteration run must show `χ_+-`
   recomputed from the **dressed** `G` (changes between iterations), proving the
   FLEX transverse builder does **not** reuse the stored bare `green0`.

## 7. Testing

- **Unit:** each new general method — output shape; degenerate input ⇒ equals
  the corresponding reduced method.
- **Brute-force equivalence (primary):** small-system `Σ[G]` from optimized path
  == direct-sum reference (~1e-10). New test module, marked slow. Run for
  **both `spin-free` and `spin-diag`** (Codex A2.3/A4.3).
- **Per-spin-mode (Codex A4.3):** the reduction, SU(2), and RPA-consistency
  checks (§6.1-§6.4) each run for `spin-free` **and** `spin-diag`.
- **Dressed-G regression (Codex A2.3):** 2-iteration run asserts `χ_+-` differs
  between iterations (uses dressed `G`, not `green0`).
- **Reduced-path regression (Codex A4.4):** existing reduced/squashed FLEX
  `energy.dat` outputs unchanged after the dispatch/guard edits.
- **Integration — reduction:** 2-orbital, `U, U'` only ⇒ general-FLEX
  `energy.dat` == reduced-FLEX reference.
- **Integration — new physics:** 2-orbital with Hund `J` ⇒ general differs from
  reduced, density-density warning is absent, symmetry checks pass.
- **Guard tests:** `enable_spin_orbital`+general and `spinful`+general fail fast
  with a clear message.
- **Tutorial:** add a `calc_scheme="general"` variant of `iron_2orb`.

Tests run from repo root, `atol≈1e-8`/`1e-10` as appropriate, matching the
existing FLEX test conventions.

## 8. Risks / mitigations

| Risk | Mitigation |
|---|---|
| **V_eff coefficient / double-counting wrong** (no in-tree baseline; **critical**) | §4.5 is a GATE: transcribe the equation from a named reference (§3) + corroborate with a second, **before** coding; brute-force from the diagrammatic definition (§6 primary, equivalence only); SU(2)/single-orbital/RPA-consistency equalities (§6 secondary) guard the formula |
| Transverse bubble uses bare `green0` mid-iteration (**critical**) | Do **not** reuse `_build_transverse_channel`; add `_build_transverse_channel_flex` taking dressed `G`; 2-iteration regression (§6.5/§7) |
| Optimized rank-4 FFT contraction bug | 1-shot match to direct-sum reference (§6 primary) |
| Regression of the reduced path | Reduced numerics unchanged; explicit reduced-path regression test (§7) |
| Performance / memory of `(nfreq, nvol, nd², nd²)` tensors | `O(nd⁴)` memory before matmul; add explicit memory/flop estimate per `norb` in the plan; decide materialized vs streamed `V_eff` (§9) |
| SOC / `spinful` interplay | Out of scope (v1); fail-fast guards + guard tests (§4.1, §7) |
| S/C reuse claim was wrong (helper is in `sc.py`, not `rpa.py`) | Extract shared `_build_sc_matrices_all_q`-based helper (§4.3) |

## 9. Open items for the implementation plan

- **Name the specific multi-orbital FLEX reference** for the `V_eff`
  longitudinal/transverse coefficients + subtraction, and **transcribe the exact
  equation into §4.5** (gate; candidates: Takimoto–Hotta–Ueda PRB 69 104504;
  Mochizuki–Yanase–Ogata cond-mat/0407094; Kubo PRB 75 224509). **User to
  confirm the reference.**
- Extract the shared full S/C builder from `sc.py:_build_sc_matrices_all_q`
  (NOT a non-existent `rpa.py` helper) for FLEX/RPA/SC reuse.
- Design `_build_transverse_channel_flex` (dressed-`G` argument) per spin mode.
- Memory/flop budget for `(nfreq, nvol, nd², nd²)`; materialized vs streamed.
- Decide the small-system parameters for the brute-force test (k-grid, Nmat).
- **Reconcile §6/§7 with the confirmed scope.** §6/§7 still list `spin-diag` +
  separate `χ_+-` targets from the earlier review round; if v1 is locked to
  paramagnetic `spin-free` (per §3/§4.5 reframing), drop the separate-transverse
  targets. **Pending user scope confirmation** (the separate spin-flip channel
  is not in the MYO/THU references).

## 10. Deferred / future work (recorded so it is not lost)

- **Tier 2 — SU(2)-broken (magnetic, `G↑≠G↓`) FLEX** with a *separate* spin-flip
  `χ_+-` self-energy channel. Not covered by MYO/THU (both paramagnetic); needs
  its own reference.
- **Tier 3 — generalized FLEX with spin-orbit coupling** ([arXiv:1503.02544],
  Sr₂IrO₄): full `2m×2m` spin-orbital generalized susceptibility, Hugenholtz
  antisymmetrized vertices, **both** particle-hole and particle-particle RPA
  channels, `Σ = δΦ/δG` (Baym-Kadanoff). A *new solver*, not an extension of the
  s/c-decomposition path; connects to H-wave `enable_spin_orbital`. Paramagnetic
  limit should reduce to v1 (Tier 1) — v1 is its validation台. Single-reference
  (weaker dual-source grounding); higher validation difficulty.
- **AL/MT 1-shot diagnostic (user, 2026-06-20).** After FLEX converges,
  optionally evaluate the Aslamazov-Larkin and Maki-Thompson vertex corrections
  **once** (non-self-consistent) and report their magnitude relative to the FLEX
  self-energy, as a *diagnostic* of whether AL/MT matter for the system — without
  including them in the self-consistency. Cheap-ish (one extra pass), genuinely
  useful (AL drives orbital/charge fluctuations à la Onari-Kontani). Needs the
  AL/MT diagram expressions in the multi-orbital matrix form; a separate small
  design of its own.
