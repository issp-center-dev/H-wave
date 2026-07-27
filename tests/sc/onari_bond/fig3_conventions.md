# Onari Fig. 3 — frozen convention record (spec P0-2 gate)

**Task:** P0-2 of `docs/superpowers/specs/2026-07-27-dynamic-bond-channels-design.md` §2.
**Status:** **GATE: MISMATCH** — see §5. Phase A (tasks 2–10) may proceed; Phase B
implementation is blocked pending spec amendment + design re-review.
**Recorded:** 2026-07-28.

## 0. Source and provenance

| | |
|---|---|
| Paper | S. Onari, R. Arita, K. Kuroki, H. Aoki, *"Phase diagram of the two-dimensional extended Hubbard model: Phase transitions between different pairing symmetries when charge and spin fluctuations coexist"* |
| Journal | Phys. Rev. B **70**, 094523 (2004) |
| Preprint used | `arXiv:cond-mat/0312314v3 [cond-mat.supr-con] 7 Oct 2004` (6 pages) |
| Retrieved from | `https://arxiv.org/pdf/cond-mat/0312314` |
| Extraction | `pdftotext -layout`, **plus page-image visual verification** of every equation transcribed below (all prefactors, signs, primes and phase arguments were read off the rendered page, not the text layer) |
| Equation inventory | The paper contains **exactly nine numbered equations, (1)–(9)**, all on p. 2 (§II FORMULATION), plus two unnumbered displays (`χ_sp`, `χ_ch`) between Eqs. (4) and (5). There are no equations anywhere else in the paper. |

All equations of the paper are reproduced or accounted for in §1 below; nothing
relevant was omitted.

## 1. What the paper actually contains (verbatim)

Notation of the paper: `k ≡ (k, ω_n)` (fermionic), `q ≡ (q, ε_n)` with
`ε_n = 2nπT` (bosonic). `χ̄` = irreducible susceptibility. `Δr (= 0, ±x̂, ±ŷ)`.

**Eq. (1) — Hamiltonian**

```
H = − Σ^nn_{i,j} Σ_σ t_ij c†_iσ c_jσ  +  U Σ_i n_i↑ n_i↓
    + ½ Σ^nn_{i,j} Σ_{σσ'} V_ij n_iσ n_jσ'
```

**Eq. (2) — Éliashberg's equation**

```
λ φ(k) = − (T/N) Σ_{k'} Γ(k, k') G(k') G(−k') φ(k')
```

"where φ is the gap function, G Green's function, and Γ the pairing interaction
with `k ≡ (k, ω_n)`. The eigenvalue λ, a measure of the pairing, becomes unity at
`T = T_C`."

**Eq. (3) — singlet pairing vertex**

```
Γ_s(k, k') =
  Σ_{Δr,Δr'} {  3/2 [V_m χ_sp V_m](k−k'; Δr; Δr') e^{i(k·Δr + k'·Δr')}
              − 1/2 [V_d χ_ch V_d](k−k'; Δr; Δr') e^{i(k·Δr + k'·Δr')}
              + 1/2  V_s(0; Δr; Δr')              e^{i(k·Δr' − k'·Δr)}  } ,
```

**Eq. (4) — triplet pairing vertex**

```
Γ_t(k, k') =
  Σ_{Δr,Δr'} { − 1/2 [V_m χ_sp V_m](k−k'; Δr; Δr') e^{i(k·Δr + k'·Δr')}
              − 1/2 [V_d χ_ch V_d](k−k'; Δr; Δr') e^{i(k·Δr + k'·Δr')}
              + 1/2  V_t(0; Δr; Δr')              e^{i(k·Δr' − k'·Δr)}  }
```

**Unnumbered display (between Eqs. 4 and 5) — dressed susceptibilities**

```
χ_sp = χ̄ / (1 + V_m χ̄) ,        χ_ch = χ̄ / (1 + V_d χ̄)
```

**Eq. (5) — irreducible (bare) bond susceptibility**

```
χ̄(q; Δr; Δr') = − (T/N) Σ_{k'} e^{i k'·(Δr − Δr')} G(k' + q) G(k')
```

**Eqs. (6), (7) — density / magnetic couplings**

```
V_d(q; Δr; Δr') = { U + 4[V_x cos(q_x) + V_y cos(q_y)] , Δr = Δr' = 0
                  { −V_x                               , Δr = Δr' = ±x̂     (6)
                  { −V_y                               , Δr = Δr' = ±ŷ

V_m(q; Δr; Δr') = { −U   , Δr = Δr' = 0
                  { −V_x , Δr = Δr' = ±x̂                                   (7)
                  { −V_y , Δr = Δr' = ±ŷ
```

"We have found that the q dependence of `V_m` and `V_d` does not in fact affect
`Γ_s` and `Γ_t` significantly."

**Eqs. (8), (9) — constant (first-order, particle–particle) terms**

```
             { 2U            , Δr =  Δr' = 0
             { V_x           , Δr =  Δr' = ±x̂
V_s(q;·;·) = { V_x e^{±i q_x}, Δr = −Δr' = ±x̂                              (8)
             { V_y           , Δr =  Δr' = ±ŷ
             { V_y e^{±i q_y}, Δr = −Δr' = ±ŷ

             {  V_x           , Δr =  Δr' = ±x̂
V_t(q;·;·) = { −V_x e^{±i q_x}, Δr = −Δr' = ±x̂                             (9)
             {  V_y           , Δr =  Δr' = ±ŷ
             { −V_y e^{±i q_y}, Δr = −Δr' = ±ŷ
```

"`V_s(0; Δr; Δr')`, `V_t(0; Δr; Δr')`, appearing in the last lines in eqs.(3,4)
respectively, are constant terms involving U and V". Note the argument used in
Eqs. (3)/(4) is **q = 0**, so the `e^{±iq}` phases evaluate to 1 there.
`V_t` has **no `Δr = Δr' = 0` entry** (no `2U` for triplet).

**Numerical parameters (p. 2, left column, verbatim)**

> "For the calculation we take an `N = 32 × 32` lattice, the temperature
> `T = 0.02`, and the Matsubara frequency for fermions
> `−(2N_c − 1)πT ≤ ω_n ≤ (2N_c − 1)πT` with `N_c = 1024`."

**Matrix dimension**

> "When the off-site interaction V is introduced all the vertices
> (`V_m, V_d, V_s, V_t`) as well as the susceptibilities become
> `(Z+1)×(Z+1)` matrices for the lattice coordination number `Z(= 4` for the
> square lattice)."

**Fig. 3 caption (the acceptance target)**

> "FIG. 3. The maximum eigenvalue, λ, of Éliashberg's equation for the
> spin-singlet (solid line) and triplet (dotted) channel as a function of V for
> `n = 0.7` with the dominant orbital symmetry indicated. The CDW phase is
> identified from the divergence in the charge susceptibility."

Footnote 24 (attached to the Fig. 3 sentence):

> "A kink in λ for singlet pairing at `V ≃ 1.2` corresponds to a change in the
> functional form in the singlet pairing channel."

Related (Fig. 4): the triplet state is the `Γ⁻₅` irrep, "may be called f in that
the gap function `(∼ sin(k_x) + sin(k_y) + const.[sin(2k_x) + sin(2k_y)])`
changes sign as `+ − + − +−` around the Fermi surface", shown for `n = 0.7`,
`V = 1.3`; the `d_{x²−y²}` panel is `n = 0.7`, `V = 0.5`.

**What the paper does NOT contain.** A full-text scan (case-insensitive) returns
**zero** occurrences of: *self-energy*, *self energy*, *Dyson*, *chemical
potential*, *self-consist(ent)*, *Hartree*, *Fock*, *first-order*, *cutoff*,
*tail*, *convergen(ce)*. There is **no self-energy equation, no Dyson equation,
no effective-interaction-for-Σ equation, no μ/density procedure, no
tail/asymptotic treatment, no SCF convergence criterion and no mixing scheme**
anywhere in the paper. The self-energy side of the FLEX loop is delegated
entirely by citation:

> "We adopt the FLEX developed by Bickers *et al.*¹⁷⁻²⁰ …"
> "Esirgen *et al.*¹³,²¹,²² have extended the FLEX method to general lattice
> Hamiltonians including the extended Hubbard model. **Following them we
> introduce the pairing interaction,**" [→ Eqs. (3)–(9)]

i.e. the paper adopts Esirgen–Bickers only for the *pairing interaction* it then
writes out; it never states what it adopted for Σ. Refs. 13, 21, 22 are
G. Esirgen, H. B. Schüttler, N. E. Bickers, PRL **82**, 1217 (1999);
G. Esirgen, N. E. Bickers, PRB **55**, 2122 (1997); PRB **57**, 5376 (1998).

## 2. Verbatim transcription into the spec's notation

Sign mapping (static spec §4.3): **`S = −V_m`, `C = V_d`** — sc.py dresses the
spin channel as `I − Sχ̄` where Onari uses `I + V_m χ̄`.

| Onari | Spec object | Transcription |
|---|---|---|
| Eq. (5) | `bond_bubble_dynamic`, dyn-spec §3.1 | `χ̄_{(m,idx),(m',idx')}(q,iν) = −(T/N) Σ_{k,ñ} e^{i k·(Δr_m − Δr_{m'})} G(k+q, iω_{ñ+j̃}) G(k, iω_ñ)` |
| Eq. (7) → `S` | `S_bond`, static spec §4.3 | `S_00 = +U` ; `S_mm = +V_x (m=±x̂)`, `+V_y (m=±ŷ)` ; bond-diagonal |
| Eq. (6) → `C` | `C_bond`, static spec §4.3 | `C_00 = U + 4[V_x cos q_x + V_y cos q_y] = U + 2V_ab(q)` ; `C_mm = −V_x / −V_y` |
| `χ_sp`, `χ_ch` | dyn-spec §3.2 | `chi_s = (I − χ̄S)^{-1}χ̄ = χ̄(I − Sχ̄)^{-1}` ; `chi_c = (I + χ̄C)^{-1}χ̄ = χ̄(I + Cχ̄)^{-1}` |
| Eq. (3) line 1–2 | `F_s`, dyn-spec §3.3 | `F_s = +1.5·S χ_s S − 0.5·C χ_c C` (since `Sχ_sS = V_m χ_sp V_m`) |
| Eq. (4) line 1–2 | `F_t`, dyn-spec §3.3 | `F_t = −0.5·S χ_s S − 0.5·C χ_c C` |
| Eqs. (8),(9) | `V^pp_η = Q_η L_η Q_η + Q_η(D_off+D_off†)Q_η`, static spec §4.5 | single band, real inversion-symmetric ⇒ `V^pp_η = L_η(m=0) + D_off(1 + ηP)`: singlet `2U` at `m=0`, `+V` diagonal, `+V` on `m ↔ m̄`; triplet `0` at `m=0`, `+V` diagonal, `−V` on `m ↔ m̄` |
| Eq. (2) | dyn-spec §3.3 normative matvec | `Y = −(T/N) Σ_{k',n'}[ F_η(k−k', iω_n−iω_{n'}) e^{i(k·Δr_m + k'·Δr_{m'})} + ½ V^pp_η e^{i(k·Δr_{m'} − k'·Δr_m)} ] X(k', iω_{n'})`, `X = G Δ G` |

Every line of this table is a **MATCH** against the corresponding spec equation
(checked term by term, including prefactors `3/2, −1/2, −1/2, −1/2, +1/2`, the
`−(T/N)` in Eqs. (2) and (5), the two distinct phase structures, the bond-
diagonal vs `q`-dependent placement of `V`, and the absence of a `2U` triplet
entry). Notable confirmations:

- **Per-spin convention.** Eq. (5) carries no spin factor 2, and
  `χ_sp = χ̄/(1 − Uχ̄)` at `V = 0`; this is the single-spin (`[↑↑]`) convention.
  Confirms dyn-spec §1.1 "per-spin components ([↑↑], not the spin trace)".
- **`C_00` Hartree term.** `4V(cos q_x + cos q_y) = 2·V(q)` with
  `V(q) = Σ_R V(R)e^{iq·R}`; matches static spec §4.3's `+2·V_ab(q)` exactly.
- **Fock signs.** Eq. (6)/(7) off-site entries `−V_x` (both channels) map to
  `S^F = +V`, `C^F = −V` — exactly the static spec §4.3 "sign result".
- **`V^pp` parity structure.** `D_off(1 + ηP)` reproduces Eq. (8) `(+V, +V)` and
  Eq. (9) `(+V, −V)` for the `{+R, −R}` block, and `Q_t` annihilates the `m=0`
  block, reproducing Eq. (9)'s missing `2U`. The static spec's derived pp vertex
  is the paper's, not a reinterpretation.
- **No second-order subtraction in Γ.** Onari's Eqs. (3)/(4) contain no
  `−(S+C)χ̄(S+C)/4`-type double-counting subtraction, matching dyn-spec §3.3
  `F_s/F_t`. (The subtraction in dyn-spec §4.3 lives on the **Σ** side, where
  the paper is silent — see §3 item 1.)
- **Pair weight.** Eq. (2)'s `G(k')G(−k')` with `k ≡ (k,ω_n)` is the spec §2
  metric `w(k,iω_n) = G(k,iω_n)G(−k,−iω_n)`.
- **Per-parity maximum.** Fig. 3 plots the maximum λ *within each of the singlet
  and triplet channels*, matching dyn-spec §1.1 "algebraic-largest λ per parity
  sector".

## 3. The four P0-2 answers

### Item 1 — the Σ effective-interaction equation

**Answer: the paper contains no Σ equation. → MISMATCH (unsatisfiable from this
source).**

Spec §4.3 states:

```
W(q,iν) = 1.5 · S_bond χ_s(q,iν) S_bond
        + 0.5 · C_bond χ_c(q,iν) C_bond
        − 0.25 · (S_bond + C_bond) χ̄(q,iν) (S_bond + C_bond)
```

and spec §4.2 states the external-leg phase attachment:

```
Σ(k,iω_n) = (T/N) Σ_{q,ν} Σ_{m,m'} e^{+i k·Δr_m} · W_{mm'}(q,iν)
            · e^{−i (k−q)·Δr_{m'}} · G(k−q, iω_n−iν)
```

There is **nothing in Onari to quote against either**. The paper writes only the
*pairing* vertex (Eqs. 3–4) and its ingredients (Eqs. 5–9); it never writes Σ,
never writes a Dyson equation, and never states which of the three FLEX
`W`-conventions (channel form with `−¼(S+C)χ̄(S+C)`, `W·K·W` reduced form, or the
Bickers `3/2 V_s + 1/2 V_c − ...` form) it used. The Σ convention is delegated by
citation to Esirgen *et al.* refs. 13, 21, 22, which are **not** transcribed
here (out of scope for this task; PRB 55/57 are paywalled and not on arXiv, the
USC mirror `physics.usc.edu/~esirgen/papers/pdf/flex1.pdf` returns HTML, not a
PDF).

What *can* be said, and is the strongest statement the paper supports:

- the **prefactor pattern** `3/2` (spin) and `1/2` (charge) with `S`/`C`
  sandwiching, and the `(Z+1)×(Z+1)` bond-matrix structure, are the paper's
  (Eqs. 3–4, 6–7) and are reproduced by the first two lines of the spec's `W`;
- the **`−0.25·(S+C)χ̄(S+C)` double-counting subtraction is not from Onari** — it
  is inherited from H-wave's own general-path FLEX (`flex.py:2219–2223`);
- the **external-leg phase attachment in Σ is not from Onari either**. It is
  *consistent with* Eq. (3)'s ph-channel attachment `e^{i(k·Δr + k'·Δr')}` under
  the pp→ph leg re-orientation (one leg conjugated because the outgoing leg is
  `k−q` rather than a second incoming pair leg), but the paper fixes only the
  **pp** attachment, so the Σ attachment remains an unvalidated design choice.

Consequence for the spec: dyn-spec §4.2's sentence *"The paper's Σ equation
(P0-2) is cross-checked against this form"* is **false as written** and must be
amended.

### Item 2 — is a k-dependent first-order (Fock) V term included in Σ?

**Answer: the paper does not say. → MISMATCH (unsatisfiable from this source).**

The words *Hartree*, *Fock*, *first-order*, *self-energy* and *chemical
potential* do not appear in the paper. There is therefore no statement to
compare with spec §4.4's default (`k`-dependent first-order Fock from V
**excluded**; `sigma_fock_bond` as the pre-specified opt-in).

Two adjacent facts, recorded so they are not mistaken for an answer:

- The paper **does** include first-order `U` and `V` terms in the **pairing**
  vertex, explicitly and with bond structure: Eqs. (8)/(9) contribute
  `2U`, `±V` including the `Δr = −Δr'` (exchange-orientation) entries. That is
  the pp channel, not Σ; the spec already reproduces it (`V^pp_η`). It is
  **not** evidence about Σ.
- The Esirgen–Bickers "general lattice Hamiltonian" FLEX formulation that Onari
  cites is a Φ-derivable/Baym–Kadanoff construction for the *full* instantaneous
  two-body interaction, in which the first-order Hartree–Fock contribution to Σ
  is conventionally present. This is a **plausible inference from the cited
  formalism, not a claim by Onari, and it has not been verified against the
  Esirgen papers in this task.** It must not be used as the convention freeze.
  If it turns out to be true, the correct H-wave response is exactly the
  pre-specified `sigma_fock_bond = true` branch (spec §4.4) — which is why the
  practical risk to Phase B is bounded even though the gate fails.

### Item 3 — Matsubara cutoff convention for `N_c = 1024`

**Answer: stated explicitly, and the window convention MATCHES the spec's
frequency map — but the spec's numeric value is wrong by a factor 2.**

Paper (verbatim): "the Matsubara frequency for fermions
`−(2N_c − 1)πT ≤ ω_n ≤ (2N_c − 1)πT` with `N_c = 1024`."

With `ω_n = (2n+1)πT`, that window is `−N_c ≤ n ≤ N_c − 1`, i.e.

> **2·N_c = 2048 fermionic Matsubara frequencies**, symmetric about zero, with
> one more negative than positive frequency.

Spec §3.1's uniform map is `ñ = n − Nmat/2`, `ω_ñ = (2ñ+1)πT`,
`ñ ∈ {−Nmat/2, …, Nmat/2 − 1}` ⇒ `ω ∈ [−(Nmat−1)πT, (Nmat−1)πT]`. Equating
windows gives **`Nmat = 2·N_c = 2048`**. The structural convention (symmetric,
`Nmat/2` negative and `Nmat/2` positive, extra slot on the negative side) is a
**MATCH**; `N_c` is the paper's *cutoff parameter*, not the frequency count.

Therefore, in the spec:

- §1.1 "Nc = 1024 Matsubara frequencies" is a **mis-paraphrase** → should read
  "`N_c = 1024` ⇒ `Nmat = 2N_c = 2048` fermionic frequencies".
- §5.3's example TOML `Nmat = 1024` is **half the paper's window** → should be
  `Nmat = 2048`.
- Corroboration: the prior H-wave reduced-FLEX benchmark carried into the
  fixtures was run at `Nmat = 2048` (dyn-spec §1.2), i.e. already at the paper's
  window.

**Tail treatment: the paper says nothing.** No tail/asymptotic/high-frequency
correction, no `coeff_tail`-equivalent, is mentioned. The spec's circular
(wrap-around) uniform convention and its IR `u_zero_plus` instantaneous
convention (§3.3) are therefore **unconstrained by the paper** — neither
validated nor contradicted.

**Bosonic axis**: the paper states only `ε_n = 2nπT` (text under Eq. 7), which
matches spec §3.1's `ν_j̃ = 2πT·j̃`. No bosonic window size is given.

### Item 4 — μ / density procedure

**Answer: the paper does not say. → unconstrained.**

The paper fixes the **band filling** `n` (`n = 0.7` for Fig. 3; `n` is the
vertical axis of the Fig. 2 phase diagram) but never states how the chemical
potential is obtained, never says whether it is re-determined from the dressed
`G` or held at its non-interacting value, and never states a density-convergence
tolerance. "chemical potential" and "self-consistent" do not appear in the text.

The spec §1.1/§4.4 choice — "μ re-determined from the dressed G each FLEX
iteration (existing `_find_mu_dressed`), density fixed to `n = 0.7`" — is
standard FLEX practice and is **not contradicted**, but it is **not confirmed by
the paper** and cannot be recorded as a frozen paper convention. It is hereby
frozen as an **H-wave convention choice**, labelled as such.

## 4. Ambiguities in the paper (recorded, not resolved)

1. **Matrix ordering of the dressing.** `χ_sp = χ̄/(1 + V_m χ̄)` is scalar-style
   notation for `(Z+1)×(Z+1)` matrices and does not fix left vs right division.
   Two readings: (a) `χ̄ (1 + V_m χ̄)^{-1}` (right division — the literal reading
   of the numerator-on-the-left notation), (b) `(1 + V_m χ̄)^{-1} χ̄`. Spec §3.2
   uses `chi_s = solve(I − χ̄S, χ̄) = (I − χ̄S)^{-1}χ̄`, which by the push-through
   identity equals `χ̄(I − Sχ̄)^{-1}` = **reading (a)**. Recorded as a MATCH under
   the literal reading; reading (b) differs in general (`S_bond` is bond-diagonal
   but not proportional to the identity: `U` vs `V_x` vs `V_y`). Not resolvable
   from the paper.
2. **`V_s`/`V_t` `q`-dependence is unused.** Eqs. (8)/(9) define a `q`-dependent
   object, but Eqs. (3)/(4) evaluate it only at `q = 0`, where the `e^{±iq}`
   phases are 1. The spec's frequency- and `q`-independent `V^pp_η` is therefore
   correct for reproducing the paper, but the paper's own `q`-dependent
   definition is left unexplained.
3. **Anisotropy-only labels.** Eqs. (6)–(9) are written for the tetragonal
   `(V_x, V_y)` case; the Fig. 3 run is isotropic (`V_x = V_y = V`,
   `t_x = t_y = 1`), used throughout §III A.
4. **`V_m`, `V_d` `q`-dependence declared unimportant.** "We have found that the
   q dependence of `V_m` and `V_d` does not in fact affect `Γ_s` and `Γ_t`
   significantly." This is an observation, not an approximation they apply — the
   `q`-dependent `C_00` is written in Eq. (6) and should be kept.
5. **Singlet-branch symmetry change at `V ≈ 1.2`** (footnote 24). The singlet
   curve in Fig. 3 has a kink there from a change of the dominant *singlet* form
   factor. This is inside the even-parity branch and does **not** move the
   even/odd crossing; but it means the "dominant orbital symmetry indicated" on
   the singlet branch is not a single irrep across the sweep. Relevant to the
   spec §1.2 branch-identity clause: the clause compares **parity** branches,
   which stays well defined.
6. **Stale equation numbers in the predecessor spec.** The static spec
   (`2026-07-25-…`) cites "Onari Eq. 14" (§4.2) and "Onari Eqs. 12–13" (§4.4).
   The paper has **only nine equations**. The correct citations are **Eq. (5)**
   for the bubble and the **unnumbered display between Eqs. (4) and (5)** for the
   dressed susceptibilities. The dynamic spec's §1.1 citations (bubble Eq. 5,
   vertices Eqs. 6–7, pairing vertices Eqs. 3–4, first-order triplet Eq. 9) are
   all **correct**.

## 5. Gate decision

**GATE: MISMATCH.**

Spec §2 P0-2 makes it a *hard gate for Phase B implementation* to record "the Σ
effective-interaction equation (verbatim, transcribed into this spec's
notation)" and "whether the k-dependent first-order (Fock) V term is included in
Σ". **Neither can be produced from Onari *et al.*: the paper contains no
self-energy equation and no statement about first-order terms in Σ at all.**

This is reported as MISMATCH, not as a pass, because the gate's decision
procedure treats the paper as the authority for §4.3's `W` and §4.4's first-order
ownership, and the paper is not an authority for either. Concretely:

- §4.2's claim "The paper's Σ equation (P0-2) is cross-checked against this form"
  is **unexecutable**; the Σ external-leg phase assignment is pinned by the SOPT
  test and the Σ oracle **only**, with no external reference.
- §4.3's `W` equation is validated by the paper **only in its first two terms'
  prefactors and vertex sandwiching** (which are Onari's Γ prefactors, an
  analogy, not the same object); the `−0.25(S+C)χ̄(S+C)` term is H-wave's own
  convention.
- §4.4's exclusion of the k-dependent Fock term is an H-wave convention with
  **no paper support either way**. The pre-specified `sigma_fock_bond` toggle
  covers the alternative, so the *implementation* risk is bounded — but the
  *convention freeze* the gate demands has not happened.
- Independently, item 3 found a **concrete factual error** in the frozen
  comparison target: `N_c = 1024` means `Nmat = 2048`, while spec §1.1/§5.3 say
  1024. Left uncorrected, the Phase B production run would use half the paper's
  Matsubara window.

**Required actions (spec §2):** stop before Phase B; amend the spec; re-run
design review. Amendments needed:

1. §1.1 / §5.3: `Nmat = 2N_c = 2048` (window `−(2N_c−1)πT ≤ ω_n ≤ (2N_c−1)πT`).
2. §4.2: delete/replace "The paper's Σ equation (P0-2) is cross-checked against
   this form"; state that the Σ phase assignment is pinned solely by the SOPT
   expansion (§4.5) and the tiny-grid Σ oracle (§4.6).
3. §4.3: label the `W` equation as an **H-wave convention** (general-path FLEX,
   `flex.py:2219–2223`), not an Onari transcription; keep the Onari-sourced
   prefactor/vertex structure claim, dropping the rest.
4. §4.4: state that Onari is silent on first-order terms in Σ; either (a) freeze
   `sigma_fock_bond = false` as an explicit H-wave choice with a documented
   caveat on the Fig. 3 comparison, or (b) resolve the convention from Esirgen
   *et al.* PRB **55**, 2122 (1997) / PRB **57**, 5376 (1998) / PRL **82**, 1217
   (1999) and freeze accordingly. Option (b) requires obtaining those papers
   (institutional access; not on arXiv).
5. §1.1: mark the μ/density procedure as an H-wave convention choice, not a
   paper convention.
6. Static spec §4.2/§4.4: fix the "Onari Eq. 14" / "Eqs. 12–13" citations to
   Eq. (5) and the unnumbered `χ_sp`/`χ_ch` display.

**Phase A is unaffected** and may proceed: every equation Phase A depends on
(bubble Eq. 5, vertices Eqs. 6–7, `F_s`/`F_t` from Eqs. 3–4, `V^pp` from
Eqs. 8–9, the `−(T/N)` Eliashberg form Eq. 2, the two phase structures, the
`w = G(k)G(−k)` metric) is present in the paper, transcribed in §1–§2 above, and
**MATCHES the spec term by term**.

## 6. Frozen values for downstream use (Task 9 and Phase B fixtures)

| Item | Frozen value | Source |
|---|---|---|
| Lattice / dispersion | 2D square, `t_ij = 1.0` NN, `a = 1`, no `t′` | Eq. (1) + text |
| `U` | 4 (`U/t = 4`) | Fig. 2 caption, §III A text |
| `V` | isotropic NN `V_x = V_y = V`, swept to the CDW boundary (Fig. 3 shows up to `V ≈ 1.3`) | Eqs. (6)–(9), Fig. 3/4 |
| Bond set | `Δr ∈ {0, ±x̂, ±ŷ}`, `(Z+1)×(Z+1) = 5×5`, `Z = 4` | text after Eq. (9) |
| Filling | `n = 0.7` | Fig. 3 caption |
| `T` | 0.02 | p. 2 |
| k-mesh | `N = 32 × 32` | p. 2 |
| Fermionic Matsubara | `N_c = 1024` ⇒ **`Nmat = 2048`**, window `−(2N_c−1)πT … (2N_c−1)πT` | p. 2 |
| Bosonic Matsubara | `ε_n = 2nπT` (window size unstated) | text under Eq. (7) |
| χ convention | per-spin (no factor 2 in Eq. 5) | Eq. (5) |
| λ definition | max λ **within each of singlet / triplet**, `λ → 1` at `T_c` | Eq. (2), Fig. 3 caption |
| Triplet state | `Γ⁻₅` ("f"), `∼ sin k_x + sin k_y + c[sin 2k_x + sin 2k_y]` | Fig. 4 + text |
| Singlet kink | `V ≃ 1.2`, change of singlet form factor (not a parity crossing) | footnote 24 |
| CDW boundary | identified by divergence of `χ_ch`; `V ≈ 1.3` at `n = 0.7` | Fig. 3/5 captions |
| Σ convention | **NOT STATED IN THE PAPER** | — |
| First-order Fock in Σ | **NOT STATED IN THE PAPER** | — |
| Tail treatment | **NOT STATED IN THE PAPER** | — |
| μ procedure | **NOT STATED IN THE PAPER** (H-wave choice: `_find_mu_dressed` on dressed G, `n` fixed) | — |
