# Full-vertex multi-orbital FLEX (`calc_scheme="general"`)

**Date:** 2026-06-16 (rev. 2026-06-20)
**Branch:** `feature/full-vertex-flex`
**Status:** DESIGN — **scope locked: v1 = paramagnetic (`spin-free`) only.**
Codex review rounds 1–2 incorporated. V_eff formula pinned to Mochizuki–Yanase–
Ogata, corroborated by Takimoto–Hotta–Ueda (§3). Round-2 fixes: §6.1 corrected
(general == reduced only for **single orbital**; multi-orbital is a diagnostic
difference, not equality — `Ûˢ≠Ûᶜ`); §4.4 index wiring **frozen in physical
indices** with the brute-force as independent ground truth; §4.5 first-order
constants **excluded** to match the reduced convention. Implementation-ready once
the user approves; the plan's first task re-confirms the MYO index placement from
the PDF.

## 1. Problem

H-wave's FLEX solver only runs in the density-density approximation. The whole
self-energy pipeline (`_inflate_chi0q_and_ham` → `_build_spin_charge_vertices` →
`_solve_rpa` → `_calc_veff` → `_calc_self_energy`) operates on the reduced
`(nvol, nd, nd)` interaction obtained by the diagonal contraction `kaabb->kab`
(`flex.py:384,401,409`). This **drops the orbital off-diagonal vertices** of
Exchange, Hund, PairLift and PairHop. `_init_flex_param` (`flex.py:84-115`)
therefore (a) rejects `calc_scheme="general"` outright and (b) warns whenever
Exchange/PairHop are present that they are reduced to their density-density part.

The density-density treatment is a recognised simplification but is **not** what
full multi-orbital FLEX requires: the inter-orbital Hund and pair-hopping vertices
enter the spin/charge fluctuations and affect the self-energy. We want full-vertex
FLEX available as an option.

## 2. Goal / scope

Implement **paramagnetic full-vertex multi-orbital FLEX**: the spin/charge
self-energy with full rank-4 `(nd², nd²)` Kanamori vertices (`U, U', J, J'`),
selected via `calc_scheme="general"` (consistent with RPA, where `general`
already means the full-tensor path).

Decisions locked in (from brainstorming):

- **Physics scope (v1):** paramagnetic (SU(2)-symmetric, `G↑=G↓`), `spin-free`
  only. The spin fluctuation enters as the single matrix `χ_s` with coefficient
  `3/2` (= 1 longitudinal Sz + 2 transverse S±, equal by SU(2)); **no separate
  transverse channel is computed** — it is bundled in `3/2 χ_s`. This is exactly
  the standard multi-orbital FLEX of the references in §3.
- **Selection:** reuse `calc_scheme="general"`; no new FLEX-specific flag.
- **Definition of done (validation):** primary = a brute-force, no-approximation
  reference that matches the optimized self-energy functional `Σ[G]` in a single
  shot; secondary = degenerate-reduction + single-orbital + RPA-consistency
  equalities. See §6.
- **Follow-up (separate design):** the SOC-capable **generalized FLEX** (Tier 3,
  §10) is the agreed next phase and subsumes magnetic (`G↑≠G↓`) and the explicit
  transverse channel; v1 is its paramagnetic-limit validation台.

### Out of scope (explicit, v1)

- **Magnetic / SU(2)-broken (`G↑≠G↓`) FLEX and any separate transverse `χ_+-`
  self-energy.** Deferred to the generalized phase (§10). In v1 a fail-fast
  guard (§4.1) rejects `spin-diag` and `spinful` under `calc_scheme="general"`.
- **SOC / `enable_spin_orbital` with general-FLEX.** Deferred to §10; rejected by
  the same guard.
- **AL / MT vertex corrections.** Outside the FLEX class entirely (documented in
  the manual; a 1-shot diagnostic is a future idea, §10). "Full vertex" = *not
  density-density reduced*, **not** *exact*.
- **External numerical reference data.** None available; the brute-force
  reference (§6) is the implementation-correctness anchor.
- **Reduced/squashed (density-density) FLEX numerics.** Unchanged and pinned by a
  regression test (§7); only the shared control flow (`_init_flex_param`,
  `solve()` dispatch) is edited, so the path is regression-guarded rather than
  literally untouched.

## 3. Reference formulation

- **V_eff — pinned, dual-sourced.** Primary (user-confirmed): **Mochizuki,
  Yanase, Ogata [cond-mat/0407094]** (multi-orbital FLEX, paramagnetic).
  Corroborating: **Takimoto, Hotta, Ueda [cond-mat/0309575] = PRB 69, 104504
  (2004)**. Both give the **same coefficient structure** (§4.5). THU writes the
  charge channel as an orbital channel and the subtraction with the cross-spin
  bare bubble `χ̄^(σσ̄)`, but the coefficients `3/2, 1/2, −1/4` and the
  `+(3/2)Ûˢ−(1/2)Ûᶜ` constants match. Both are paramagnetic / SU(2)-symmetric
  and carry **no separate spin transverse channel** — confirming the v1 design.
- **S/C interaction matrices (full Kanamori).** Kuroki et al. [arXiv:0902.3691]
  and THU use the full Kanamori vertex in `(norb², norb²)` matrix form. H-wave
  already builds these S (spin) / C (charge) matrices in
  `sc.py:_build_sc_matrices_all_q` (`sc.py:422-523`) for the SC pairing vertex.
  **Caveat:** `_compute_vertices_general` is a *pairing* vertex, **not** a
  self-energy kernel — it validates the S/C matrices and χ construction but does
  not fix the FLEX self-energy coefficients (those are §4.5).

Honest caveat: RPA/FLEX resums the bubble+ladder series but omits AL/MT vertex
corrections; "full vertex" means *not density-density reduced*, not *exact*.
Documented in `docs/{en,ja}/source/flex/tutorial/tu-index.rst`.

## 4. Architecture

Keep the reduced path; add a **parallel general path** dispatched from `solve()`
on `self.calc_scheme == "general"`. Because v1 is paramagnetic, the path is the
straightforward rank-4 generalisation of the existing reduced methods — **no new
transverse machinery**.

### 4.1 Dispatch / guards (`_init_flex_param`)

- Replace the `general`-rejection branch: `general` → full-vertex path;
  `reduced`/`squashed` → existing density-density path (unchanged default).
- Move the "Exchange/PairHop reduced to density-density" warning so it fires
  **only** for `reduced`/`squashed` (under `general` nothing is dropped).
- Require the general-scheme rank-4 `chi0q`; fail loudly if a reduced `chi0q` is
  supplied in general-FLEX.
- **Fail-fast guard (v1):** reject `calc_scheme="general"` unless
  `spin_mode == "spin-free"` — i.e. reject `spin-diag`, `spinful`, and
  `enable_spin_orbital` with a clear message pointing to the (future) generalized
  solver (§10).

### 4.2 New general methods (mirroring the reduced ones)

| Reduced (keep) | New general (rank-4) |
|---|---|
| `_inflate_chi0q_and_ham` | `_inflate_chi0q_and_ham_general` |
| `_build_spin_charge_vertices` | `_build_spin_charge_vertices_general` |
| `_calc_veff` | `_calc_veff_general` |
| `_calc_self_energy` | `_calc_self_energy_general` |

Each unit has one clear purpose, a shape-typed interface, and a degenerate limit
that reproduces the reduced method (testable in isolation — §7). **No transverse
builder** is needed in v1.

### 4.3 Reuse from SC / RPA (corrected per Codex review)

- **S/C interaction matrices.** There is **no** named S/C builder in `rpa.py`
  (its general inflation is inline in `solve()`, `rpa.py:875-999`). The reusable
  full `(norb², norb²)` S/C builder is **`sc.py:_build_sc_matrices_all_q`
  (`sc.py:422-523`)**. Plan: **extract a shared helper** (common module) used by
  FLEX-general, RPA-general, and SC — rather than a non-existent `rpa.py` reuse.
  Verify the matrices equal MYO/THU Eq.(10) element-wise.
- `_solve_rpa` already inverts general-scheme matrices (flattens rank-4 to
  `(ndx, ndx)`, `rpa.py:2068-2138`); verify on the rank-4 `χ_s`, `χ_c`.
- Reuse the existing general-scheme `chi0q` from `_calc_chi0q` (RPA) and the
  spin inflation; for `spin-free` paramagnetic this is straightforward and needs
  no Green's-function-resolved transverse bubble.

### 4.4 Core change: the self-energy contraction

The reduced self-energy is an **element-wise (Hadamard)** product
`Σ_{ab}(r,τ) = V_eff,{ab}(r,τ) · G_{ab}(r,τ)` (`flex.py:642-645`). The general
self-energy is the **rank-4 contraction**, per MYO Eq.(3):

```
Σ_{mn}(r,τ) = Σ_{μν} V_{μm,νn}(r,τ) · G_{μν}(r,τ)
```

**Index convention — pinned now in PHYSICAL orbital indices (Codex re-review).**
The ground-truth contraction is MYO Eq.(3) in explicit orbital indices:

```
Σ_{mn}(r,τ) = Σ_{μν} V_{μm,νn}(r,τ) · G_{μν}(r,τ)
```

with `V` the rank-4 object built from MYO Eq.(4) and `χ⁰_{mn,μν}(q) =
−(T/N) Σ_k G_{μm}(k+q) G_{nν}(k)` (Eq.5). Note the wiring: the **contracted**
indices are `μ` (first slot of `V`'s first pair, = first slot of `G`) and `ν`
(first slot of `V`'s second pair, = second slot of `G`); the **output** indices
`m, n` are the *second* slot of each `V` pair. This is **not** a plain
`V_matrix · vec(G)` in the `(mn),(μν)` flatten that `χ`/S/C use
(`sc.py:463-468`) — it needs a deliberate relabel/transpose.

**How this is locked (avoids Codex's "shared implicit convention" trap):** the
brute-force reference (§6) implements Eq.(3) **directly in physical indices**
(nested `m,n,μ,ν` loops, no flatten), so it does **not** share the optimized
path's flatten convention — it is an independent ground truth. The optimized
path's flatten/transpose must reproduce it. The plan will record the chosen
flatten and a worked `norb=2` example, but the *physics wiring above is frozen
here* and the brute-force enforces it. The FFT transport (Matsubara↔τ, k↔r) of
`_calc_self_energy` is reused unchanged; only the per-`(r,τ)` product becomes the
batched orbital contraction with this wiring.

**Caveat on the source:** these index placements are transcribed from an ar5iv
rendering of MYO; the plan's first task is to **re-derive Eq.(3)–(5) index
placement from the paper PDF** and confirm before coding, since a mis-transcribed
slot would propagate into both the formula and the brute-force.

### 4.5 `V_eff` assembly (`_calc_veff_general`) — formula (MYO, dual-sourced)

Interaction matrices `Ûˢ` (spin), `Ûᶜ` (charge) from MYO Eq.(10) in the
orbital-pair basis (intra `U`, inter `U'`, Hund `J_H`, pair-hopping `J'=J_H`):

```
χ⁰_{mn,μν}(q) = −(T/N) Σ_k G_{μm}(k+q) G_{nν}(k)            (Eq.5)

χˢ(q) = [I − χ⁰(q) Ûˢ]⁻¹ χ⁰(q)                              (Eq.4)
χᶜ(q) = [I + χ⁰(q) Ûᶜ]⁻¹ χ⁰(q)

V(q) =  (3/2) Ûˢ χˢ(q) Ûˢ
      + (1/2) Ûᶜ χᶜ(q) Ûᶜ
      − (1/4)(Ûˢ+Ûᶜ) χ⁰(q) (Ûˢ+Ûᶜ)        ← double-counting subtraction
      + (3/2) Ûˢ − (1/2) Ûᶜ                ← first-order (Hartree-Fock) constants

Σ_{mn}(k) = (T/N) Σ_q Σ_{μν} V_{μm,νn}(q) G_{μν}(k−q)        (Eq.3)
```

Notes / cross-checks:
- The `3/2 χˢ` term already includes the transverse (spin-flip) contribution by
  SU(2) — no separate `χ_+-` is computed (the reason v1 is clean).
- `Ûˢ, Ûᶜ` must equal `sc.py:_build_sc_matrices_all_q` (the shared helper, §4.3) —
  verify element-wise vs MYO Eq.(10).
- `_solve_rpa` sign mapping (`[1 + χ⁰·ham]⁻¹χ⁰`): `ham_s = −Ûˢ`, `ham_c = +Ûᶜ`
  (same convention as the reduced `_build_spin_charge_vertices`).
- **First-order constants `+(3/2)Ûˢ − (1/2)Ûᶜ` — DECISION (Codex re-review):
  EXCLUDE them from `_calc_veff_general`**, matching the reduced `_calc_veff`
  which omits them (`flex.py:491-532`). H-wave handles the first-order (Hartree-
  Fock) shift outside `V_eff` (static term / `μ`); `V_eff` carries only the
  *fluctuation* part. `_calc_veff_general` therefore implements only the three
  fluctuation terms (`3/2 ÛˢχˢÛˢ + 1/2 ÛᶜχᶜÛᶜ − 1/4(Ûˢ+Ûᶜ)χ⁰(Ûˢ+Ûᶜ)`). This is
  recorded as an **acceptance criterion**: the brute-force reference (§6) also
  omits the constants, so both compare the fluctuation self-energy. (If a later
  check shows H-wave does *not* otherwise capture the multi-orbital first-order
  term, revisit — but v1 matches the established reduced convention.)

## 5. Data flow (general path, one SCF iteration)

```
H0, Σ_in ─▶ dressed G(k,iω)              [_calc_dressed_green, reused]
        ─▶ χ⁰(q,iν) rank-4              [_calc_chi0q general, RPA]
        ─▶ inflate χ⁰, build Ûˢ,Ûᶜ      [_inflate_*_general, _build_*_general]
        ─▶ χ_s, χ_c                     [_solve_rpa, general]
        ─▶ V(q,iν) rank-4              [_calc_veff_general, MYO Eq.4]
        ─▶ Σ_out(k,iω) = Σ V·G          [_calc_self_energy_general, rank-4 contraction]
        ─▶ mix, check convergence       [reused]
```

## 6. Validation = definition of done

**Primary — no-approximation 1-shot match.** A standalone brute-force reference
computes one evaluation of `Σ[G]` by **direct summation** (no FFT convolution, no
density-density reduction): `χ⁰` by explicit k-sum (MYO Eq.5), `χ_s, χ_c` by
matrix RPA, and `Σ` by the explicit `q, iν` double sum (MYO Eq.3). The optimized
general path must match this to ~1e-10 on a small system (e.g. 2 orbitals, 4×4
k-grid, few Matsubara). One shot suffices (self-consistency only iterates
`Σ[G]`). "Slow is OK" — test-only. The brute-force is written from MYO Eqs.(3-5)
**independently** of the optimized `V_eff` code path.

Honest scope (Codex A2.1): the 1-shot match proves **algorithmic equivalence**
(optimized path computes the same functional as the direct sum), not by itself
the physical correctness of the formula. With the formula now dual-sourced (§3)
and the brute-force written straight from MYO's equations, plus the secondary
checks below, the correctness case is complete for v1.

**Secondary — limits (precise equalities). NOTE (Codex re-review): the
brute-force 1-shot is the real correctness anchor; the reduced path is a
*different, cruder* approximation, so general-FLEX only equals it in the
single-orbital limit, not for multi-orbital density-density.**

1. **Single orbital == reduced:** general-FLEX `χ_s, χ_c, Σ, energy` ==
   reduced-FLEX to ~1e-10. Here `Ûˢ = Ûᶜ = U`, so the MYO per-channel form and
   the reduced single-`W` form coincide (fluctuation part). This is the only
   exact general-vs-reduced equality.
2. **Multi-orbital density-density — DIAGNOSTIC, not equality.** general-FLEX with
   only `U, U'` (no Hund/PairHop) will **differ** from reduced-FLEX, because the
   reduced `_calc_veff` sandwiches a single full `W` on both sides
   (`flex.py:521-539`) whereas MYO uses per-channel `Ûˢ ≠ Ûᶜ` (`Ûˢ=U'`, `Ûᶜ=−U'`
   / `2U'` off-diagonal, `sc.py:478-509`). Record the difference as expected;
   do **not** assert equality. (Pre-existing observation: the reduced
   `_build_spin_charge_vertices` builds separate `u_s, u_c` at `flex.py:467-470`
   but `_calc_veff` then uses the full `W` — an internal inconsistency in the
   reduced path; out of scope to fix here, but it is *why* reduced ≠ general.)
3. **RPA consistency:** first-iteration `χ_s`, `χ_c` from general-FLEX == RPA
   general-scheme `χ` for the same `χ⁰`/interaction (~1e-10).
4. **S/C matrix identity:** `Ûˢ, Ûᶜ` == MYO/THU Eq.(10) element-wise (compare to
   `sc.py:_build_sc_matrices_all_q`); `χ_s, χ_c` Hermitian.

## 7. Testing

- **Unit:** each new general method — output shape; degenerate input ⇒ equals the
  reduced method.
- **Brute-force equivalence (primary):** small-system `Σ[G]`, optimized ==
  direct-sum (MYO Eqs.) to ~1e-10. New test module, marked slow.
- **Integration — single-orbital reduction:** 1-orbital ⇒ general-FLEX
  `energy.dat` == reduced-FLEX reference (the only exact general-vs-reduced
  equality, §6.1).
- **Integration — multi-orbital diagnostic (NOT equality):** 2-orbital `U, U'`
  only ⇒ general-FLEX **differs** from reduced-FLEX (record the difference;
  assert it is nonzero and bounded, not that it vanishes).
- **Integration — new physics:** 2-orbital with Hund `J` ⇒ general differs from
  reduced; density-density warning absent; S/C-matrix and RPA-consistency checks
  pass.
- **Reduced-path regression:** existing reduced/squashed FLEX `energy.dat`
  outputs unchanged after the dispatch/guard edits.
- **Guard tests:** `spin-diag`+general, `spinful`+general, and
  `enable_spin_orbital`+general each fail fast with a clear message.
- **Tutorial:** add a `calc_scheme="general"` variant of `iron_2orb`.

Tests run from repo root, `atol≈1e-8`/`1e-10` as appropriate, matching the
existing FLEX test conventions.

## 8. Risks / mitigations

| Risk | Mitigation |
|---|---|
| V_eff formula correctness | Dual-sourced (MYO + THU, §3); brute-force from MYO Eqs. (§6 primary); degenerate/single-orbital/RPA-consistency equalities (§6 secondary) |
| Optimized rank-4 FFT contraction bug | 1-shot match to direct-sum reference |
| Regression of the reduced path | Reduced numerics unchanged; explicit reduced-path regression test (§7) |
| First-order (Hartree) handling differs from reduced | Reconcile per §4.5; degenerate-reduction test catches mismatch |
| Performance / memory of `(nfreq, nvol, nd², nd²)` tensors | `O(nd⁴)` memory before matmul; add explicit memory/flop estimate per `norb` in the plan; decide materialized vs streamed `V_eff` |
| S/C reuse claim was wrong (helper is in `sc.py`) | Extract shared `_build_sc_matrices_all_q`-based helper (§4.3) |
| Misuse on magnetic/SOC systems | Fail-fast guards + guard tests (§4.1, §7) |

## 9. Open items for the implementation plan

- **First task — confirm MYO Eq.(3)–(5) index placement from the paper PDF**
  (not just the ar5iv summary), since a mis-transcribed index slot would corrupt
  both the formula and the brute-force (§4.4 caveat).
- Extract the shared full S/C builder from `sc.py:_build_sc_matrices_all_q` for
  FLEX/RPA/SC reuse; verify against MYO/THU Eq.(10).
- Reconcile first-order (Hartree) treatment between `_calc_veff_general` and the
  reduced `_calc_veff` so the degenerate-reduction check holds.
- **Pin the orbital index wiring (§4.4)** between the `(mn),(μν)` χ/S/C matrix
  convention and MYO Eq.(3) `V_{μm,νn} G_{μν}`; record the worked map and lock it
  with the brute-force check.
- Memory/flop budget for `(nfreq, nvol, nd², nd²)`; materialized vs streamed.
- Small-system parameters for the brute-force test (k-grid, Nmat).

## 10. Deferred / future work (recorded so it is not lost)

- **Generalized FLEX with spin-orbit coupling (agreed next phase)**
  ([arXiv:1503.02544], Sr₂IrO₄): full `2m×2m` spin-orbital generalized
  susceptibility, Hugenholtz antisymmetrized vertices, **both** particle-hole and
  particle-particle RPA channels, `Σ = δΦ/δG` (Baym-Kadanoff). A *new solver*,
  not an extension of the s/c-decomposition path; connects to H-wave
  `enable_spin_orbital`. **Subsumes** magnetic (`G↑≠G↓`) and the explicit
  transverse channel (so a standalone "SU(2)-broken transverse" tier is not
  pursued separately). Paramagnetic, SOC-free limit must reduce to v1 — **v1 is
  its validation台**. Single-reference (weaker dual-source); higher validation
  difficulty; its own design + spec.
- **AL/MT 1-shot diagnostic (user, 2026-06-20).** After FLEX converges,
  optionally evaluate the Aslamazov-Larkin and Maki-Thompson vertex corrections
  **once** (non-self-consistent) and report their magnitude relative to the FLEX
  self-energy, as a *diagnostic* of whether AL/MT matter for the system. Cheap-ish
  (one pass); needs the AL/MT diagram expressions in multi-orbital matrix form;
  its own small design.
  - **I/O already in place.** The diagnostic's inputs — converged `G` (`green`),
    `χ_s`/`χ_c` (`chiq_s`/`chiq_c`), `χ⁰` (`chi0q`), `Σ` (`sigma`) — are already
    saved by `FLEX.save_results` (`flex.py:705-760`) as `.npz` with
    `freq_index`/wavevector metadata, **when the matching `[file.output]` keys are
    set** (the bundled samples set them). So it can be a post-processing tool (like
    `hwave_sc`) needing little/no new I/O. Caveats: (a) saving is key-gated;
    (b) saved `sigma`/`χ` are reduced `(nd,nd)` today and general rank-4 under
    full-vertex FLEX — detect which from metadata.
