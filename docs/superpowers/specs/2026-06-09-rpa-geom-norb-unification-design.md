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

**Construction-order constraint (D3).** In `Interaction.__init__`,
`self._init_interaction()` (rpa.py:195) currently runs *before*
`self.norb = Geometry["norb"]` (rpa.py:197). The unified derivation
(`norb_phys`, even-check) and `self.norb`/`self.nd` assignment must be computed
**before** `_init_interaction()` so any norb-dependent interaction setup and the
(removed) guard logic see the physical count. Implementation: hoist the
geom-norb parsing / `norb_phys` derivation to the top of `__init__`, then call
`_init_interaction()`. Apply the same ordering in `RPA._init_param`
(rpa.py:638) if `self.norb` is consumed before its assignment there.

**Why this fixes the fold path automatically.** `_reshape_interaction`'s
`_reshape_orbit_spin` (rpa.py:295) decodes `s_ = a // norb_orig` where
`norb_orig = Geometry["norb"]`. Under the new convention `norb_orig` = SO count,
and every interaction spin-orbital index `a` satisfies `a < SO count`, so
`s_ = 0` always — the helper collapses to the opaque-generalized-orbital fold,
**identical to UHFk's verified behavior** (uhfk.py:372, where `norb_orig` is
already the SO count). The divergence disappears, which is the precondition for
removing the guard.

## 4. Audit table — `rpa.py` `norb` / `nd` / raw-geom sites

`norb_phys` = physical orbital count; `SOcount` = `2·norb_phys` = `geom_norb` in
SO mode.

| Site | Current code | Intended dimension | Action |
|---|---|---|---|
| `Interaction.__init__` `rpa.py:197` | `self.norb = Geometry["norb"]` | physical | SO: `self.norb = geom_norb//2`; even-check |
| `RPA._init_param` `rpa.py:638` | `self.norb = geometry["norb"]` | physical | SO: `self.norb = geom_norb//2`; even-check |
| `rpa.py:639-640` | `self.ns=2; self.nd=self.norb*2` | `nd = SOcount` | OK once `self.norb=norb_phys` |
| `_init_interaction` guard `rpa.py:232` | `Geometry["norb"] > 1` | — | remove (test-first, §5) |
| `_reshape_geometry` `rpa.py:265,268,274` | `norb = geom['norb']`; `center` sized `norb*bvol`; `geom_new['norb']=geom['norb']*bvol` | SO mode: W90 gives one center per spin-orbital → SO-count centers is correct for the *input*; verify no downstream consumer assumes physical-count centers | audit + test; expected no-op but must confirm |
| `_reshape_interaction` `rpa.py:290` | `norb_orig = Geometry["norb"]` | SO count (for `s_=0` fold) | OK; becomes SO-count → matches UHFk |
| `_make_ham_trans` `rpa.py:377` | `norb = self.norb` | physical | OK |
| `_make_ham_inter` `rpa.py:466` | `norb = self.norb` | physical | OK |
| `_read_chi0q` shape checks `rpa.py:1185,1201,1218,1232` | classify `nd==self.nd`→spinful, `nd==self.norb`→spin-free | `self.nd=SOcount`, `self.norb=norb_phys` | OK once derivation fixed; add SO regression |
| `_calc_trans_mod` `rpa.py:1390,1393` | `norb=self.norb` | physical | OK |
| `_calc_epsilon_k` `rpa.py:1415,1448-1449` | `nd=self.nd; norb=self.norb` | nd=SOcount, norb=physical | OK |
| `RPA._reshape_green` `rpa.py:1306-1349` (`self.norb_orig`, 1316) | `green` reshaped `(Lvol,ns,norb_orig,ns,norb_orig)` | physical (ns=2) | **Resolve D1:** `self.norb_orig` is **never assigned** in rpa.py and this method has **no caller** (only `self.lattice._reshape_green` at rpa.py:1281 is invoked). Confirm dead, then **delete `RPA._reshape_green`**; if instead it is to be revived, assign `self.norb_orig = norb_phys` in SO mode. |
| `_make_ham_inter` `Extern`/two-body `rpa.py:442` | physical-orbital indexed allocation | physical | **Resolve D5:** specify SO-mode behavior; add test (§5) |

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
   - **Extern / two-body in SO mode (D5)**: pin behavior of the
     `_make_ham_inter` allocation (`rpa.py:442`) under the SO-count convention.
3. **Implement** §3/§4 until those go GREEN.
4. **Delete** the guard only after green; keep the odd-`norb` fail-fast.

If an equivalence test cannot be made to pass (folded path genuinely wrong),
stop and re-scope — do not delete the guard.

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
- **Dead `RPA._reshape_green` / unassigned `self.norb_orig`** (D1): delete the
  dead method during the refactor to avoid leaving a latent `AttributeError`;
  if revived, bind `self.norb_orig = norb_phys`.
- **Implementation must re-verify every audit-table line number** before
  writing tasks (D2 showed two stale citations in the first draft).
