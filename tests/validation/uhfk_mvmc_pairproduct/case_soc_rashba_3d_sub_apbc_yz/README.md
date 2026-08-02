# case_soc_rashba_3d_sub_apbc_yz (v3.7 shipping fixture)

xy-plane Rashba SOC + general complex z-direction spin-mixing hopping
fixture with antiperiodic boundary condition in y and z (yz plane
APBC, x periodic). Common lattice across all 4 v3.7 shipping fixtures.

## Geometry / interaction

- `CellShape` = `[4, 4, 4]`
- `SubShape` = `[2, 2, 2]` (folded BZ `L_folded = [2, 2, 2]`)
- `BoundaryCondition` = `["periodic", "antiperiodic", "antiperiodic"]` (2D APBC in yz plane)
- 24-row Transfer.dat: xy Rashba (16 rows from v3.6) + z spin-diagonal hopping (4 rows, t_z = -1) + general complex z-direction spin-mixing hopping (4 rows, 0.3+0.4j).
- U = 2 on-site Hubbard
- Total physical sites: 64 (spin-orbital dim 128)

The z-direction spin-mixing block is not Rashba SOC. It gives
`H_z(k_z) = (0.6 cos k_z - 0.8 sin k_z) sigma_x`, whose even-in-`k_z`
term breaks spin-1/2 time-reversal symmetry. The full hopping remains
Hermitian, so this fixture is valid for testing the mapping and exercises
general complex hopping. A time-reversal-symmetric 3D SOC fixture is a v3.8
follow-up.

## Ncond selection

Phase 2b: `soc_pilot_gap.py` under this fixture's own
`BoundaryCondition` (theta = (0, pi, pi)):

```
PYTHONPATH=src python tests/validation/uhfk_mvmc_pairproduct/scripts/soc_pilot_gap.py \
    --case tests/validation/uhfk_mvmc_pairproduct/case_soc_rashba_3d_sub_apbc_yz \
    --ncond 12,14,16,18,20,22,24
```

Measured with `flag_fock = true` (see `input.toml`). The Fock term
changes the converged spectrum, so this table supersedes the earlier
Fock-off one; the pinned filling is unchanged.

- `Ncond= 12`: gap = 3.795e-01 (PASS)
- `Ncond= 14`: gap = 5.002e-02 (PASS, marginal -- just above the 5e-2 floor)
- `Ncond= 16`: gap = 2.933e-01 (PASS)
- `Ncond= 18`: gap = 5.170e-02 (PASS, marginal -- just above the 5e-2 floor)
- `Ncond= 20`: gap = 4.628e-01 (PASS)
- `Ncond= 22`: gap = 1.245e-02 (FAIL; below the 5e-2 floor)
- `Ncond= 24`: gap = 5.737e-01 (PASS)

(The script's own `[OK]`/`[REJECT]` marker applies a laxer `gap >=
1e-3` threshold plus an even-`Ncond` parity check per its docstring;
the `>= 5e-2` cutoff above is the stricter spec §1.2 selection
criterion used to pick the pinned value. Note `Ncond=22` prints
`[OK]` from the script but does not clear `>= 5e-2` here.)

This fixture's gap-passing candidates (`gap >= 5e-2`) are `{12, 14,
16, 18, 20, 24}` -- six of the seven candidates, which is precisely
why the scalar gap check is not a sufficient selection criterion for
this fixture. `Ncond = 20` was INITIALLY selected across the 4 v3.7
siblings on the gap criterion alone (cross-fixture physics
consistency, and inside spec §1.2's "expected range 16-24 for
closed-shell filling"); the partner-balance analysis below is what
moved this fixture off that value.

### Ncond=20 blocked Phase 2d: partner-balance invariant (v3.7 Task 2d fix)

`soc_pilot_gap.py`'s gap check above is a SCALAR HOMO-LUMO check only.
It does not verify `build_pair_list`'s SOC pairing invariant, which
`general_fij_builder.py` requires independently for every non-self
canonical/partner k-row pair: `n_occ(k) == n_occ(partner(k))`. At the
initially-pinned `Ncond = 20`, this fixture VIOLATES that invariant
(canonical k-row 2, partner row 1: `n_occ_k=2 != n_occ_p=4`), which
raises inside `phase2_produce_manifest.py` (via
`soc_apbc_topology_guard._build_A_ship_mutated` ->
`build_pair_list`) and blocked Phase 2d composite-manifest generation
entirely -- the gap check alone cannot see this because a scalar
HOMO/LUMO gap says nothing about how the occupied states distribute
across the two rows of a canonical pair.

A read-only frozen-spectrum re-slice (same methodology as the gap
table above: re-derive occupation at each candidate `Ncond` via aufbau
filling of the Ncond=20-converged spectrum, `step_occupation` +
`find_partner_rows` + `build_pair_list`, without rerunning SCF at each
candidate) against all 7 gap-table candidates gives:

- `Ncond= 12`: gap PASS, partner-balance FAIL (k=2/partner=1: `2 != 1`)
- `Ncond= 14`: gap PASS, partner-balance FAIL (k=6/partner=5: `2 != 1`)
- `Ncond= 16`: gap PASS, partner-balance **PASS**
- `Ncond= 18`: gap PASS, partner-balance FAIL (k=2/partner=1: `2 != 3`)
- `Ncond= 20`: gap PASS, partner-balance FAIL (k=2/partner=1: `2 != 4`)
- `Ncond= 22`: gap FAIL, partner-balance FAIL (k=2/partner=1: `3 != 4`)
- `Ncond= 24`: gap PASS, partner-balance **PASS**

Candidates satisfying BOTH `gap >= 5e-2` AND partner-balance: `{16,
24}`. `Ncond = 24` is pinned (over 16) for more gap headroom
(5.737e-01 vs 2.933e-01) and to stay closer to the xy/xz siblings'
shared `Ncond = 20` for cross-fixture comparability.

Re-verified under `flag_fock = true`: the usable set is still `{16,
24}` and the pinned `Ncond = 24` still satisfies both criteria, so
enabling the Fock term does not move this fixture's filling. (Under
the superseded Fock-off spectrum, `Ncond = 18` and `22` failed the
gap check on Fermi-level degeneracy; with Fock on, 18 clears the gap
but fails partner-balance and 22 fails both. The usable set is
unchanged either way.)

**Caveat**: as with the gap table above, this is a frozen-spectrum
estimate (aufbau re-fill of one self-consistent solution's spectrum),
not independent SCF reruns per candidate -- the authoritative check is
the actual SCF run + `phase2_produce_manifest.py` run at the pinned
`Ncond`, both of which are re-verified below. `soc_pilot_gap.py` itself
was intentionally left unextended with a partner-balance check (v3.7
Task 2d Step 7 was deemed non-essential given this fixture's own
README now documents the caveat explicitly); a future fixture's Ncond
selection should scan for partner-balance as well as gap, not gap
alone.

Pinned: `Ncond = 24` (gap = 5.737e-01 >= 5e-2, partner-balance holds).
Fresh H-wave SCF at `Ncond = 24` with `flag_fock = true` confirms
convergence (45 iterations, rest=6.478e-15,
`Energy_Total=-91.080493088...`) and `phase2_produce_manifest.py`
completes with a self-check-passing 30-entry manifest (composite
`i=7 s=0 j=57 t=1`, `|G_c|=2.8956e-02`). Its active y/z axes contribute
20 distinct, positive-threshold mutation evaluations; the 10 inactive
x-axis entries remain at zero per spec §4b. The original M-4 sign flip
was structurally invisible on `L_folded=2`, leaving 16 effective entries.
An interim fix made M-ship-4 and M-ship-5 both omit `sub_offset`, so only
18 of the 20 active entries were distinct. The current M-ship-4 instead
halves the named-axis `sub_offset` contribution, while M-ship-5 omits it.
The producer rejects structurally degenerate active mutations and
threshold-policy drift before writing a manifest. Measured 2026-07-19.

## v3.7 gate role

Phase 6 seven-gate run: `bash tests/validation/uhfk_mvmc_pairproduct/run.sh case_soc_rashba_3d_sub_apbc_yz` produces the G0-writer-check + G1 + G2a-emitted-F + G2a-in-memory-A + G2b + G3 + G4 records at 1e-6 tolerance for the SOC + APBC(yz) + SubShape[2,2,2] triple validation.

**G4 guard fix (v3.7 Task 2d)**: this fixture is the only one of the 4
whose active APBC axes exclude x (`active_apbc_axes = [y, z]`). Its
spec-compliant composite (selected by `_pick_composite_v37`, which
correctly requires `dr != 0` along every ACTIVE direction only, per
spec §4 addendum) differs in `sub_offset` along y and z but NOT x
(`sub_offset_i=(0,1,0)`, `sub_offset_j=(0,0,1)`). This exposed a latent
bug in `soc_apbc_topology_guard.py`'s semantic re-check: it had a
defensive `sub_offset_x(i) != sub_offset_x(j)` assertion hardcoded to
axis 0, a carryover from the v3.6 single-axis (always-x) design that
was never generalized when the v3.7 per-direction addendum was added
to the producer. It silently under-verified xy/xz/xyz (whose composite
happens to also differ along x, an active axis for all three) and
outright rejected this fixture's fully spec-compliant composite
(checking an inactive axis instead of the true active ones). Fixed by
generalizing the check to iterate over every axis with `theta_d != 0`
(derived from the manifest's own `theta_radians`), matching the
producer's `_pick_composite_v37` criterion exactly. Verified no
regression on xy/xz/xyz or the v3.6 shipping fixture
(`case_soc_rashba_2d_sub_apbc`); regression tests added to
`tests/test_soc_apbc_topology_guard_semantic.py`.

## Phase 2a fixture-format note

`Geometry.dat` and `CoulombIntra.dat` are byte-identical copies of the
v3.6 `case_soc_rashba_2d_sub_apbc` fixture's files, not scaled to 64
sites: both are parsed via `wan90.read_w90` / `wan90.read_geom` in the
primitive-cell-orbital basis (`norb = 2`, one physical site's up/down
spin-orbitals), and `uhfk.py` combines that `norb` with `CellShape`
from `input.toml` (`Nstate = Lx*Ly*Lz*norb`) to build the full 64-site
system -- the on-site `U = 2` term and the 1-site primitive cell are
unchanged from v3.6, so these two files do not depend on `CellShape`
at all. `Transfer.dat`'s ndegen block (declared count 24) wraps at 15
entries per line (15 + 9), matching `read_w90`'s `nints_per_line = 15`
contract and the `emit_rashba_transfer.py` generator convention;
`OneBodyG.dat` uses the standard 5-line header
(`emit_onebodyg_full_soc.py --cell 4,4,4`), matching the `skiprows=5`
contract documented in that script and used by v3.6's fixture. All 5
of these files are byte-identical across the 4 v3.7 fixtures.
