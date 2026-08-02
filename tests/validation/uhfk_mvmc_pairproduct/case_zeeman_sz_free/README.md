# case_zeeman_sz_free — 1D PBC L=8 Zeeman-Hubbard Sz-free mVMC E2E fixture

Sz-free (no `2Sz` constraint) 1D Hubbard with `flag_fock = false` and a
single shared mu-group; `column_spin` stays in {0, 1} (no mixed block).
v3 B ケース.

## Status (v3)

Runnable end-to-end (harness-gated). See `scripts/assert_occupation.py`
`_CASE_TARGETS["case_zeeman_sz_free"]` for the target occupation.

Precondition: pytest + `run.sh case_zeeman_sz_free` both call
`assert_occupation.py case_zeeman_sz_free` after SCF. Update the
target dict there if the fixture Transfer.dat / physical parameters
change.

Finalized target occupation (single source of truth:
`_CASE_TARGETS["case_zeeman_sz_free"]` in `scripts/assert_occupation.py`):

- `up`:   k=(0, 0, 0), (+1, 0, 0)=+π/4, (-1, 0, 0)=-π/4
- `down`: k=(0, 0, 0)

§3.2 pair-closure walkthrough over canonical blocks:
- self-pair k=(0,0,0): `NN_up=1, NN_down=1` → `n_cross=1`, both
  excesses zero (even) ✓
- self-pair k=(4,0,0)=π: empty ✓
- non-self (+π/4, -π/4): canonical=(-1,0,0), partner=(+1,0,0);
  `NN_up_k = NN_up_p = 1`, `NN_down = 0` → `excess_up_k = excess_up_p
  = 1` (same-spin up-up pair) ✓; down excess balanced 0/0 ✓
- non-self (±π/2), (±3π/4): empty ✓
- Total Ne = 3 + 1 = 4 (even) ✓

## Known limitation — Zeeman entries could not be materialized

The plan intended `Transfer.dat` to carry on-site Zeeman-like entries
`0 0 0 1 1 -0.3 0.0` (up on-site) and `0 0 0 2 2 +0.3 0.0` (down
on-site) with the file's `orb ∈ [norb+1, 2*norb]` convention encoding
the spin block. In H-wave's `uhfk.py`:

- `_check_hermite` (`uhfk.py:877`) allocates `tab_r` with shape
  `(nx, ny, nz, norb, norb)` and hard-errors with `IndexError` when a
  Transfer entry has `orbvec ≥ norb`.
- `_make_ham_trans` (`uhfk.py:1136-1146`) silently drops Transfer
  entries with `orbvec ≥ norb` when `enable_spin_orbital = false`,
  making the Zeeman terms no-ops even if the Hermite check were
  bypassed.
- `enable_spin_orbital = true` would let those entries take effect,
  but `tools/uhfk_to_mvmc.py:91-93` explicitly rejects that mode in v1
  (`ERROR: enable_spin_orbital is not supported in v1`).

Consequently `Transfer.dat` carries only the standard NN hopping. The
"Sz-free B ケース" role of this fixture is preserved because the
`input.toml` still omits `2Sz` (single shared mu-group) and
`flag_fock = false` keeps `column_spin ∈ {0, 1}`; H-wave converges
from its random initial Green function to a symmetry-broken FM state
that happens to coincide with the case_pbc_sz2 occupation but is
reached via the Sz-free (no `2Sz`) branch of the bridge, exactly the
dispatch row the B ケース is meant to cover.

## Known limitation (v3 in progress)

Full `run.sh case_zeeman_sz_free` currently stops at the bridge
invocation with `N_up = 3 != N_down = 1; v1 spec section 7 requires
AntiParallel Sz-fixed sector`. The pytest harness-gate (Step 1 +
Step 1.5) is green:

- Step 1 (H-wave UHFk SCF): converges to the target FM state.
- Step 1.5 (`scripts/assert_occupation.py case_zeeman_sz_free <work>`):
  passes.
- Step 2 (`vmcdry.out`): writes `namelist.def`, `orbitalidxgen.def`
  (Sz-free branch of StdFace also emits general format when the state
  is spin-imbalanced), ...
- Step 3 (`tools/uhfk_to_mvmc.py --check-density`): fails fast because
  the v1 parallel-Slater dispatch cannot ingest the general-format
  orbital file. The General dispatch is wired in Task 11 and the full
  E2E gate follows in Task 12.

`stan.in` has `2Sz` removed (Sz-free), and `NSPGaussLeg` is not
specified (default projection off), matching the input.toml with no
`2Sz` key.
