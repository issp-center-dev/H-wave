# RPA / UHFk geometry-`norb` convention unification (spin-orbital)

**Date:** 2026-06-09
**Branch:** `feature/flex-approximation`
**Status:** design approved, spec for review

## 1. Problem

In spin-orbital (SO) mode (`enable_spin_orbital = true`), the two k-space
solvers disagree on what `Geometry["norb"]` in `geom.dat` means:

| Solver | `Geometry["norb"]` interpreted as | internal physical orbitals | `nd` (total basis) | SO transfer file index count |
|---|---|---|---|---|
| **UHFk** | **SO count** (`2·norb_phys`) | `norb//2` (`uhfk.py:256`) | `norb·ns = norb` (ns=1) | `2·norb_phys` |
| **RPA** | **physical count** (`norb_phys`) | `norb` (`rpa.py:197,638`) | `norb·2` | `2·norb_phys` |

Consequence: for an `N`-physical-orbital SO system, UHFk requires
`geom.dat norb = 2N` while RPA requires `geom.dat norb = N`. **The same input
files are not portable between the two solvers.** The final `nd` agrees
(`2N`), so the bug is invisible whenever `norb_phys = 1` (every pre-existing
RPA SO test), which is why it went unnoticed.

Wannier90 SOI output uses `num_wann = total spin-orbitals = 2N`. UHFk follows
this convention; **RPA is the outlier.**

This convention split is also the root cause that currently forces the
fail-fast guard in `rpa.py` `_init_interaction` (commit `bd777a6`) rejecting
SO + sublattice-folding + multi-orbital: the interaction-folding helper
`_reshape_orbit_spin` behaves differently between the two solvers purely
because `norb_orig` means different things (see §3).

## 2. Goal & scope

**Goal:** unify on the **SO-count** convention (UHFk / Wannier90), so the same
`geom.dat` + transfer files work in both solvers, AND remove the SO + sublattice
+ multi-orbital guard once the folded path is verified to match UHFk.

**In scope:**
- Change RPA to interpret `Geometry["norb"]` as SO count in SO mode, deriving
  the physical orbital count `norb_phys = norb // 2` (even-check, fail-fast on
  odd).
- Audit every raw `Geometry["norb"]` / `self.norb` / `self.nd` site in `rpa.py`
  for physical-vs-SO intent (§4).
- Test-first removal of the SO+sublattice+multiorbital guard (§5).
- Migrate RPA SO test fixtures to SO-count geom files (§6).
- Update RPA SO documentation (EN/JA) to state the SO-count convention (§7).

**Out of scope:**
- UHFk and the I/O layer (already SO-count / W90-compliant).
- Non-SO behavior (unchanged).
- Any change to the transfer-file index convention (stays interleaved
  `2*orb+spin`, 1-based in file) — only the `geom.dat norb` *value* changes.

## 3. Key invariant (the design decision)

**RPA internally keeps `self.norb` = physical orbital count and
`self.nd = 2·norb_phys = SO count`, exactly as today.** The only change is the
*derivation* of `self.norb` from `geom.dat` in SO mode:

```
non-SO mode (unchanged):   self.norb = geom_norb            ; ns = 2 ; nd = 2·norb
SO mode (changed):         norb_phys = geom_norb // 2       ; ns = 2 ; nd = 2·norb_phys = geom_norb
                           self.norb = norb_phys            (even-check geom_norb, fail-fast on odd)
```

This is a *narrow derivation change at the init sites*, not a sweeping rewrite,
because RPA already treats `self.norb` as the physical orbital count everywhere
downstream (`_calc_epsilon_k` reshape `(nvol, ns, norb, ns, norb)`, transverse
slices `[:norb]`/`[norb:2·norb]`, `_solve_rpa`). Those sites stay correct
unchanged once `self.norb` carries the physical count.

**Construction-order constraint (D3 + N1) — derive AFTER folding.** This is
subtle: `_init_interaction()` (rpa.py:195) *folds the geometry* when
`has_sublattice` — `_reshape_geometry` replaces `param_ham["Geometry"]["norb"]`
with `orig_norb × subvol` (rpa.py:249,268). Only afterwards does rpa.py:197 read
`self.norb = param_ham["Geometry"]["norb"]` (the **post-fold** value), and then
`_make_ham_trans`/`_make_ham_inter` (rpa.py:200-201) consume `self.norb`.

Therefore the `norb_phys` derivation must operate on the **post-fold** value and
must stay **after** `_init_interaction()` — *not* hoisted before it (an earlier
draft of this spec said "hoist before"; that is wrong because it would drop the
`× subvol` factor). `_init_interaction` and `_reshape_interaction` read raw
`param_ham["Geometry"]["norb"]` / lattice attributes and do **not** reference
`self.norb`, so leaving the assignment after folding introduces no dependency
problem (verified: neither method uses `self.norb`/`self.norb_phys`).

Implementation at rpa.py:197 (and the analogous `RPA._init_param` rpa.py:638 —
verify whether its `param_ham` is folded or raw, and match):

```
geom_norb = param_ham["Geometry"]["norb"]     # post-fold (= SO_count × subvol in SO mode)
if enable_spin_orbital:
    # even-check target (P1): param_ham_orig exists ONLY when has_sublattice
    # (rpa.py:240-242); otherwise use the current (unfolded) geom norb.
    check_norb = param_ham_orig["Geometry"]["norb"] if has_sublattice else geom_norb
    require check_norb even -> else fail-fast ValueError naming geom.dat value
    self.norb = geom_norb // 2                # = norb_phys × subvol
else:
    self.norb = geom_norb
self.ns = 2 ; self.nd = self.norb * 2
```

**`RPA._init_param` fold-state contract (P3).** `RPA.__init__` constructs the
`Interaction` (rpa.py:579) — which mutates the *shared* `param_ham["Geometry"]`
in place when `has_sublattice` — *before* calling `_init_param` (rpa.py:583).
So `_init_param` (rpa.py:638) reads the **post-fold** geom norb when
`has_sublattice`, and the **raw** value otherwise — the same value
`Interaction` saw. Apply the identical SO derivation there
(`self.norb = geom_norb // 2` in SO mode), keeping the two paths consistent by
construction (shared mutated dict). State this explicitly; do not leave it as
"verify and match".

**Two distinct orbital-index spaces (the root of P4).** RPA mixes two index
conventions and unification must respect both:
- **Transfer** is **spin-orbital-indexed** (file index `2*orb+spin`, range
  `[0, SO_count)`); its fold stride is the SO count.
- **Two-body interactions and Extern** are **physical-orbital-indexed**
  (`_make_ham_inter` writes `a,b` into the `norb=self.norb` physical axis with
  spin carried explicitly by `spin_table`, rpa.py:491-497); their fold stride
  is `norb_phys = SO_count/2`.

Before unification (geom norb = physical) these happened to coincide for the
only unguarded SO case (`norb_phys=1`). After unification (geom norb = SO count)
they diverge by a factor of 2, so `_reshape_interaction` must use the right
stride per term — see §4 (P4).

**Why this fixes the Transfer fold path.** `_reshape_interaction`'s
`_reshape_orbit_spin` (rpa.py:295) decodes `s_ = a // norb_orig` where
`norb_orig = Geometry["norb"]`. Under the new convention `norb_orig` = SO count,
and every Transfer spin-orbital index `a` satisfies `a < SO count`, so
`s_ = 0` always — the helper collapses to the opaque-generalized-orbital fold,
**identical to UHFk's verified behavior** (uhfk.py:372, where `norb_orig` is
already the SO count). The divergence disappears for Transfer. The non-Transfer
stride is handled separately (P4).

## 4. Audit table — `rpa.py` `norb` / `nd` / raw-geom sites

`norb_phys` = physical orbital count; `SOcount` = `2·norb_phys` = `geom_norb` in
SO mode.

| Site | Current code | Intended dimension | Action |
|---|---|---|---|
| `Interaction.__init__` `rpa.py:197` | `self.norb = Geometry["norb"]` (reads **post-fold** value, after `_init_interaction` at :195) | physical (× subvol) | SO: `self.norb = post_fold_norb//2`; even-check on pre-fold backup. Keep assignment **after** `_init_interaction` (N1) |
| `RPA._init_param` `rpa.py:638` | `self.norb = geometry["norb"]` | physical | SO: `self.norb = geom_norb//2`; even-check. Verify fold state of this `param_ham` and match |
| `rpa.py:639-640` | `self.ns=2; self.nd=self.norb*2` | `nd = SOcount` | OK once `self.norb=norb_phys` |
| `_init_interaction` guard `rpa.py:232` | `Geometry["norb"] > 1` | — | remove (test-first, §5) |
| `_reshape_geometry` `rpa.py:265,268,274` | `norb = geom['norb']`; `center` sized `norb*bvol`; `geom_new['norb']=geom['norb']*bvol` | SO mode: W90 gives one center per spin-orbital → SO-count centers is correct for the *input*; verify no downstream consumer assumes physical-count centers | audit + test; expected no-op but must confirm |
| `_reshape_interaction` `rpa.py:284-297` (fold stride `norb_orig`, :290) | `norb_orig = param_ham_orig["Geometry"]["norb"]` used as fold stride for **all** terms | **per-term, differs (P4)** | **Blocker P4 — split the stride.** Transfer is SO-indexed → stride = SO count (pre-fold geom norb). Two-body / Extern are **physical-indexed** (`_make_ham_inter` writes `a,b` into the `norb=self.norb` physical axis, rpa.py:491-497) → stride must be `norb_phys = geom_norb//2` in SO mode, else folded indices overrun `self.norb=norb_phys·subvol`. Currently `_reshape_interaction` is called with `enable_spin_orbital=True` only for Transfer (rpa.py:252) and `False` for the rest (rpa.py:255), so distinguish via the solver's `self.enable_spin_orbital`: `stride = geom_norb` for the SO-indexed Transfer call, `geom_norb//2` for non-Transfer when `self.enable_spin_orbital`, `geom_norb` otherwise (non-SO, unchanged). `param_ham_orig` always exists here (only reached inside the `has_sublattice` branch). |
| `_make_ham_trans` (transfer) `rpa.py:377` | `norb = self.norb` | physical | OK |
| `_make_ham_trans` (`Extern` alloc) `rpa.py:439-455` (esp. 442) | `hab_r` sized `(nx,ny,nz,norb,norb)`, `norb=self.norb`; skips spin (comment :453) | physical | OK once `self.norb=norb_phys` (N2: this block is in `_make_ham_trans`, *not* `_make_ham_inter`) |
| `_make_ham_inter` (two-body) `rpa.py:466,473` | `norb=self.norb`, `nd=norb*ns`, tensor `(ns,norb)*4`, explicit `spin_table` | physical, ns=2 | OK once `self.norb=norb_phys` |
| `_read_chi0q` shape checks `rpa.py:1185,1201,1218,1232` | classify `nd==self.nd`→spinful, `nd==self.norb`→spin-free | `self.nd=SOcount`, `self.norb=norb_phys` | OK once derivation fixed; add SO regression |
| `_calc_trans_mod` `rpa.py:1390,1393` | `norb=self.norb` | physical | OK |
| `_calc_epsilon_k` `rpa.py:1415,1448-1449` | `nd=self.nd; norb=self.norb` | nd=SOcount, norb=physical | OK |
| `RPA._reshape_green` `rpa.py:1306-1349` (`self.norb_orig`, 1316) | `green` reshaped `(Lvol,ns,norb_orig,ns,norb_orig)` | physical (ns=2) | **D1 corrected → see B1.** This method is **not dead**: `_read_trans_mod` at rpa.py:1281 calls `self.lattice._reshape_green(tab_r)`, but `Lattice` (rpa.py:71-174) has **no** `_reshape_green` — the intended target is **`self._reshape_green`** (this RPA method). The `self.lattice.` is a mis-wire (pre-existing since 2023, commit 2d7128c). So do **not** delete it; resolve via B1. `self.norb_orig` is still never assigned → must be set to `norb_phys` (physical) if this path is fixed. |

Anything the audit finds to be physical-count-specific while reading raw
`geom['norb']` (now SO count) gets an explicit `norb_phys` substitution.
(Citations re-verified 2026-06-09; earlier draft mis-cited 1393 as
`_calc_epsilon_k` and 466 as `_make_ham_extern`.)

## 5. Guard removal (test-first)

The guard at `rpa.py:229-232` raises `NotImplementedError` for
`enable_spin_orbital and has_sublattice and Geometry["norb"] > 1`.

Order of operations (TDD):
1. **Rewrite** the guard tests in `test_rpa_so_multiorb.py`
   (`TestRPAMultiOrbitalSOSublatticeGuard`, currently asserting the exception)
   to assert *success* — RED until the refactor lands.
2. **Add** equivalence/invariance/validation tests (RED first):
   - **SO vs non-SO**: a spin-independent multi-orbital system encoded in SO
     form yields the same `chi0q` as the equivalent non-SO (`ns=2`) run, *with*
     `SubShape != [1,1,1]` sublattice folding. (Extends the existing
     `test_rpa_so_multiorb` pattern to the folded case.)
   - **Multi-orbital folded vs unfolded invariance**: `chi0q` for a
     `norb_phys ≥ 2` SO system is invariant under choice of `SubShape`
     (physical observable must not depend on the fold), mirroring
     `TestUHFkSublatticeInvarianceSO`. (Existing invariance coverage is
     1-orbital only — D5.)
   - **Odd-`norb` fail-fast (D5)**: SO mode with odd `geom_norb` raises a clear
     `ValueError`/error (not a downstream shape crash).
   - **SO transfer index range validation (D7)**: RPA must reject an SO
     transfer file whose index is outside `[0, geom_norb)` — UHFk has this
     check (exercised in `test_block_matrix.py`); RPA currently has no
     equivalent. The `_reshape_orbit_spin` no-op argument (§3) depends on every
     interaction/transfer index being `< SO count`, so this validation is a
     correctness precondition, not just hygiene.
   - **Extern / two-body in SO mode (D5)**: pin behavior of the `Extern`
     allocation in `_make_ham_trans` (`rpa.py:439-455`) and the two-body
     `_make_ham_inter` (`rpa.py:461+`) under the SO-count convention (both use
     `norb=self.norb`; expected correct once `self.norb=norb_phys`).
   - **Multi-orbital SO + sublattice + two-body interaction fold (P4)**: the
     key newly-unguarded path. A `norb_phys ≥ 2` SO system with a two-body term
     (e.g. CoulombInter) **and** `SubShape != [1,1,1]` must give the same
     `chi0q` folded vs unfolded. This is the test that catches the
     `_reshape_interaction` physical-vs-SO stride bug (P4); it must be RED
     before the stride split and GREEN after. Without it the guard removal is
     unsafe.
3. **Implement** §3/§4 until those go GREEN.
4. **Delete** the guard only after green; keep the odd-`norb` fail-fast.

If an equivalence test cannot be made to pass (folded path genuinely wrong),
stop and re-scope — do not delete the guard.

## 5b. Pre-existing `trans_mod` + sublattice bug (B1)

The audit surfaced a **pre-existing, independent** bug (since 2023, commit
`2d7128c`). `RPA._read_trans_mod` (rpa.py:1281) calls
`self.lattice._reshape_green(tab_r)` inside `if self.lattice.has_sublattice`,
but `Lattice` (rpa.py:71-174) has no `_reshape_green` method — the real method
is `RPA._reshape_green` (rpa.py:1306). So **any** `trans_mod` input combined
with sublattice folding raises `AttributeError`, regardless of SO mode. No test
exercises this path (only `test_flex_analytical.py` uses `trans_mod`, without a
folding `SubShape`), so it has gone unnoticed.

This is **orthogonal to geom-norb unification**: it breaks for non-SO and
`norb_phys=1` SO sublattice runs too, and removing the multi-orbital chi0q guard
(§5) does not newly expose it. But the spec must not let guard removal advertise
"SO + sublattice + multi-orbital supported" while this sub-path crashes.

**Decision (to confirm with user) — recommended: narrow-guard + defer.**
- **Option B1-guard (recommended):** add a focused fail-fast — `trans_mod`
  provided together with `has_sublattice` raises a clear `NotImplementedError`
  pointing at this known gap — and track the real fix as a separate follow-up
  issue. Keeps this spec scoped to the chi0q core path.
- **Option B1-fix:** correct the mis-wire to `self._reshape_green(tab_r)`,
  assign `self.norb_orig = norb_phys` (SO) / `norb` (non-SO), and add a
  `trans_mod` + sublattice (+ SO) regression test. Larger scope; pulls a 2023
  bug into this change.

Either way: **do not delete `RPA._reshape_green`** (corrects D1). Whichever
option, it becomes a gate item — guard removal in §5 must not leave a reachable
`AttributeError` in the SO + sublattice space.

## 6. Fixture migration

`geom.dat` files are shared between SO and non-SO test cases, so the SO `norb`
value cannot be doubled in place. Add **new SO-count geom fixtures**; transfer
files are unchanged (their index count was already SO-count).

| New fixture | `norb` | Replaces (for SO cases) | Used by |
|---|---|---|---|
| `tests/rpa/input/geom_so.dat` | 2 | `geom.dat` (norb=1) | `test_rpa_spin.py` SO cases; `test_rpa_1orb.py::test_U_and_V_spin_orbital` (D4) |
| `tests/rpa/input/geom_so_2orb.dat` | 4 | `geom_2orb.dat` (norb=2) | `test_rpa_so_multiorb.py` SO run |

Update each SO test config to point `Geometry` at the new SO-count file.
Affected SO-using test files (verified 2026-06-09): `test_rpa_spin.py`,
`test_rpa_1orb.py` (the `test_U_and_V_spin_orbital` case at line 278, hardcodes
`geom.dat` at line 181 — **D4**), and `test_rpa_so_multiorb.py`. Non-SO configs
keep the existing physical-count geom files; the non-SO comparison run in
`test_rpa_so_multiorb` keeps `geom_2orb.dat` (norb=2).

**No migration needed** (verified — D6): `test_rpa_ladder.py` has **no SO
full-run case** (no `enable_spin_orbital=True` solver run); its SO-adjacent
tests use `object.__new__` stubs with `norb`/`nd` set directly.
`test_block_matrix.py` likewise constructs `norb`/`nd` via stubs and is
convention-agnostic.

## 7. Documentation

Update RPA SO sections (source `.rst` only; `build/html` + `_sources` are
generated — regenerate, do not hand-edit):
- `docs/{en,ja}/source/algorithm/rpa.rst`
- `docs/{en,ja}/source/rpa/filespecification/config/index_config.rst`
- `docs/{en,ja}/source/rpa/filespecification/inputfiles/interaction.rst`
- `docs/{en,ja}/source/rpa/filespecification/outputfiles/chiq_rpa.rst`

State: in SO mode, `geom.dat`'s `norb` is the **spin-orbital count**
(`= 2·N_phys = Wannier90 num_wann`), matching UHFk; transfer index convention
is interleaved `2*orb+spin`. Add a one-line migration note (RPA SO users must
double `geom.dat norb`).

## 8. Backward compatibility

Breaking change to RPA SO input (`develop`/pre-release, acceptable). The
breaking surface is limited to the `geom.dat norb` *value* in SO mode; transfer
files and the transfer index convention are unchanged. Fail-fast on odd
`geom_norb` in SO mode gives a clear migration signal. Record in
CHANGELOG / migration note.

## 9. Verification

- New equivalence/invariance tests in §5 GREEN.
- Existing `test_rpa_so_multiorb` (unfolded) still GREEN with SO-count fixtures.
- `chi0q.npz` continues to carry `index_convention="spin_block"`.
- Full suite GREEN (currently 395).
- Codex re-review of the implementation diff.

## 10. Risks

- **Hidden physical-vs-SO assumption** in a raw `geom['norb']` site not in the
  audit table (mitigated by §4 audit + equivalence tests).
- **Folded path genuinely incorrect** even after unification (mitigated by
  test-first guard removal in §5 — stop if equivalence fails).
- **`_reshape_geometry` center array** semantics under SO-count: resolved by
  audit — raw-file norb (= SO count) gives one center per spin-orbital, correct
  for W90 SOI; keep as-is with SO-count geom files.
- **`trans_mod` + sublattice `AttributeError`** (B1, corrects D1): `RPA._reshape_green`
  is *not* dead — it is the mis-wired target of `self.lattice._reshape_green`
  (rpa.py:1281). Pre-existing since 2023, orthogonal to geom-norb, untested.
  Resolve per §5b (recommended: narrow `trans_mod`+sublattice guard + defer);
  do **not** delete the method. Bind `self.norb_orig = norb_phys` if fixed.
- **Implementation must re-verify every audit-table line number** before
  writing tasks (D2 showed two stale citations in the first draft).
- **Physical-vs-SO fold stride (P4)** is the highest-risk item: the two-body
  interaction fold stride must be `norb_phys`, not the SO count. Gated by the
  §5 multi-orbital SO + sublattice + two-body invariance test — do not remove
  the guard until it is GREEN.
