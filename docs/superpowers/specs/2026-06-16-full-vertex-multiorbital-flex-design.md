# Full-vertex multi-orbital FLEX (`calc_scheme="general"`)

**Date:** 2026-06-16
**Branch:** `develop` (feature branch to be created)
**Status:** DESIGN — awaiting user review before writing the implementation plan

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
  extraction assumes Sz conservation (`rpa.py:1987-1990` warns for genuine spin
  mixing). General-FLEX + SOC is deferred to a follow-up.
- **External numerical reference data.** None available now; the brute-force
  reference (§6) is the correctness anchor. External cross-checks are a later,
  separate task.
- **No change to the existing reduced/squashed (density-density) FLEX path.** It
  remains the default and must stay bit-for-bit unchanged.

## 3. Reference formulation (what we cross-check the formulas against)

- **Longitudinal + pairing:** the multi-orbital matrix RPA of Kuroki et al.
  (5-band RPA + Eliashberg, [arXiv:0902.3691]) and Takimoto et al.
  (PRB 69, 104504) use the **full Kanamori vertex** (`U, U', J, J'`) in
  `(norb², norb²)` matrix form — **not** density-density. H-wave's SC module
  already implements this exact form in `_compute_vertices_general`
  (`sc.py`: singlet `V^s = 3/2 S χ_s S − 1/2 C χ_c C + 1/2 (S + C)`), so the
  S/C matrix construction and the longitudinal susceptibility are a known,
  cross-checkable formulation.
- **Transverse FLEX self-energy:** must be **pinned to a specific multi-orbital
  FLEX reference** (the V_eff coefficient that combines longitudinal and
  transverse χ without double counting). This is the part with no in-tree
  baseline and is the **main correctness risk**. The pinned reference will be
  named in the implementation plan; the brute-force check (§6) is written from
  the FLEX *diagrammatic* definition, independently of the optimized V_eff.

Caveat recorded honestly: RPA/FLEX is itself an approximation (bubble+ladder
summation, no vertex corrections beyond RPA). "Full vertex" here means *not
density-density reduced*, not *exact*.

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

### 4.2 New general methods (mirroring the reduced ones)

| Reduced (keep) | New general (rank-4) |
|---|---|
| `_inflate_chi0q_and_ham` | `_inflate_chi0q_and_ham_general` |
| `_build_spin_charge_vertices` | `_build_spin_charge_vertices_general` |
| `_calc_veff` | `_calc_veff_general` |
| `_calc_self_energy` | `_calc_self_energy_general` |
| (n/a) | transverse: reuse RPA `_build_transverse_channel` → χ_+- via `_solve_rpa` |

Each unit has one clear purpose, a shape-typed interface, and a degenerate
limit that reproduces the reduced method (testable in isolation — see §7).

### 4.3 Reuse from RPA (do not duplicate)

- General-scheme spin inflation and the full S (spin) / C (charge) interaction
  matrices `(nd², nd²)` — reuse RPA's general-scheme construction.
- `_solve_rpa` already inverts general-scheme matrices; verify on rank-4 χ.
- `_build_transverse_channel` (`rpa.py:1886-2034`) returns `chi0_+-` and the
  transverse vertex `W_+-`; the per-interaction sign table there
  (`rpa.py:2011-2018`) is the SU(2) cross-check anchor.

### 4.4 Core change: the self-energy contraction

The reduced self-energy is an **element-wise (Hadamard)** product
`Σ_{ab}(r,τ) = V_eff,{ab}(r,τ) · G_{ab}(r,τ)` (`flex.py:642-645`). The general
self-energy is the **rank-4 contraction**

```
Σ_{l1 l2}(r,τ) = Σ_{l3 l4} V_{l1 l3 l2 l4}(r,τ) · G_{l3 l4}(r,τ)
```

i.e. `V_eff` is `(nd², nd²)` and contracts the internal orbital pair `(l3,l4)`
of `G` at each `(r,τ)`. The exact index ordering of `V` (shown here as
`l1 l3 l2 l4`) is fixed together with the `V_eff` convention (§3, §4.5) and is
verified by the brute-force check (§6) rather than asserted here. The FFT
transport (Matsubara↔τ, k↔r) of `_calc_self_energy` is reused unchanged; only
the per-`(r,τ)` product becomes a batched matmul over the orbital indices.

### 4.5 `V_eff` assembly (`_calc_veff_general`)

Assemble the rank-4 `V_eff(q, iν)` from:
- longitudinal spin (`3/2 · S χ_s S`), longitudinal charge (`1/2 · C χ_c C`),
- the transverse contribution (`χ_+-` with `W_+-`),
- minus the single double-counting subtraction (the reduced code subtracts `χ0`
  exactly once, `flex.py:526-532`; the general form must preserve that the
  zeroth order reproduces the second-order bubble `W χ0 W`).

The exact longitudinal/transverse coefficients are pinned per §3 and verified by
§6. This subsection is the physics-critical one.

## 5. Data flow (general path, one SCF iteration)

```
H0, Σ_in ─▶ dressed G(k,iω)              [_calc_dressed_green, reused]
        ─▶ χ0(q,iν) rank-4              [_calc_chi0q general, RPA]
        ─▶ inflate χ0, build S,C        [_inflate_*_general, _build_*_general]
        ─▶ χ_s, χ_c (longitudinal)      [_solve_rpa, general]
        ─▶ χ0_+-, W_+- ─▶ χ_+-          [_build_transverse_channel + _solve_rpa]
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

**Secondary — limits & symmetry (guards the formula choice):**

1. **Degenerate reduction:** general-FLEX with only `U, U'` (no Hund/PairHop/
   Exchange) reproduces reduced-FLEX to ~1e-10 (χ0, χ_s, χ_c, Σ, energy).
2. **Single orbital:** general-FLEX == reduced-FLEX (no off-diagonal exists);
   transverse == longitudinal in this limit.
3. **RPA consistency:** first-iteration χ_s/χ_c match RPA general-scheme χ for
   the same χ0/interaction.
4. **SU(2) limit:** for full Kanamori, `W_+- = W_zz` (`rpa.py:2018`); χ
   Hermiticity / physical-sign checks.

## 7. Testing

- **Unit:** each new general method — output shape; degenerate input ⇒ equals
  the corresponding reduced method.
- **Brute-force equivalence (primary):** small-system `Σ[G]` from optimized path
  == direct-sum reference (~1e-10). New test module, marked slow.
- **Integration — reduction:** 2-orbital, `U, U'` only ⇒ general-FLEX
  `energy.dat` == reduced-FLEX reference.
- **Integration — new physics:** 2-orbital with Hund `J` ⇒ general differs from
  reduced, density-density warning is absent, symmetry checks pass.
- **Tutorial:** add a `calc_scheme="general"` variant of `iron_2orb`.

Tests run from repo root, `atol≈1e-8`/`1e-10` as appropriate, matching the
existing FLEX test conventions.

## 8. Risks / mitigations

| Risk | Mitigation |
|---|---|
| Transverse V_eff coefficient / double-counting wrong (no in-tree baseline) | Pin to a named multi-orbital FLEX reference (§3); brute-force from the diagrammatic definition (§6 primary); SU(2)/single-orbital limits (§6 secondary) |
| Optimized rank-4 FFT contraction bug | 1-shot match to direct-sum reference |
| Regression of the reduced path | Reduced path untouched; degenerate-reduction test pins equality |
| Performance of `(nd², nd²)` matrices | Acceptable; reuse RPA batched-matmul patterns; document cost |
| SOC interplay | Out of scope; fail-fast guard for general-FLEX + `enable_spin_orbital` |

## 9. Open items for the implementation plan

- Name the specific transverse-FLEX reference for the V_eff coefficients.
- Confirm RPA's general-scheme S/C builder is directly reusable by FLEX (vs. a
  thin shared helper).
- Decide the small-system parameters for the brute-force test (k-grid, Nmat).
