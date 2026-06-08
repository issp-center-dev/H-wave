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
| `_make_ham_extern` `rpa.py:466` | `norb = self.norb` | physical | OK |
| `_read_chi0q` shape checks `rpa.py:1185,1201,1218,1232` | classify `nd==self.nd`→spinful, `nd==self.norb`→spin-free | `self.nd=SOcount`, `self.norb=norb_phys` | OK once derivation fixed; add SO regression |
| `_calc_epsilon_k` `rpa.py:1393-1394,1448-1449` | `nd=self.nd; norb=self.norb` | nd=SOcount, norb=physical | OK |
| `self.norb_orig` use `rpa.py:1316-1317` | `norb_orig=self.norb_orig` | verify where `norb_orig` is set in RPA and its dimension | audit |

Anything the audit finds to be physical-count-specific while reading raw
`geom['norb']` (now SO count) gets an explicit `norb_phys` substitution.

## 5. Guard removal (test-first)

The guard at `rpa.py:229-232` raises `NotImplementedError` for
`enable_spin_orbital and has_sublattice and Geometry["norb"] > 1`.

Order of operations (TDD):
1. **Rewrite** the guard tests in `test_rpa_so_multiorb.py`
   (`TestRPAMultiOrbitalSOSublatticeGuard`, currently asserting the exception)
   to assert *success* — RED until the refactor lands.
2. **Add** equivalence/invariance tests (RED first):
   - **SO vs non-SO**: a spin-independent multi-orbital system encoded in SO
     form yields the same `chi0q` as the equivalent non-SO (`ns=2`) run, *with*
     `SubShape != [1,1,1]` sublattice folding. (Extends the existing
     `test_rpa_so_multiorb` pattern to the folded case.)
   - **Folded vs unfolded invariance**: `chi0q` for an SO system is invariant
     under choice of `SubShape` (physical observable must not depend on the
     fold), mirroring `TestUHFkSublatticeInvarianceSO`.
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
| `tests/rpa/input/geom_so.dat` | 2 | `geom.dat` (norb=1) | `test_rpa_spin.py` SO cases, `test_rpa_ladder.py` SO cases (if any) |
| `tests/rpa/input/geom_so_2orb.dat` | 4 | `geom_2orb.dat` (norb=2) | `test_rpa_so_multiorb.py` SO run |

Update each SO test config to point `Geometry` at the new SO-count file.
Non-SO configs keep the existing physical-count geom files. The non-SO
comparison runs in `test_rpa_so_multiorb` keep `geom_2orb.dat` (norb=2).
`test_block_matrix.py` and the stub-based ladder tests construct `norb`/`nd`
directly via `object.__new__` and are convention-agnostic — no change.

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
- **`_reshape_geometry` center array** semantics under SO-count (flagged for
  audit; expected correct for W90 SOI but unverified downstream).
