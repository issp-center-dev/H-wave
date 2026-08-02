# case_soc_rashba_3d_sub_apbc_xyz (v3.7 shipping fixture)

xy-plane Rashba SOC + general complex z-direction spin-mixing hopping
fixture with antiperiodic boundary condition in x, y, and z (full APBC
on all three axes). Common lattice across all 4 v3.7 shipping
fixtures.

## Geometry / interaction

- `CellShape` = `[4, 4, 4]`
- `SubShape` = `[2, 2, 2]` (folded BZ `L_folded = [2, 2, 2]`)
- `BoundaryCondition` = `["antiperiodic", "antiperiodic", "antiperiodic"]` (APBC on x, y, and z)
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
`BoundaryCondition` (theta = (pi, pi, pi)):

```
PYTHONPATH=src python tests/validation/uhfk_mvmc_pairproduct/scripts/soc_pilot_gap.py \
    --case tests/validation/uhfk_mvmc_pairproduct/case_soc_rashba_3d_sub_apbc_xyz \
    --ncond 12,14,16,18,20,22,24
```

Measured with `flag_fock = true` (see `input.toml`), anchored to this
fixture's own pinned `Ncond = 12` SCF solution. The Fock term changes
the converged spectrum, so this table supersedes the earlier Fock-off
one; the pinned filling is unchanged.

- `Ncond= 12`: gap = 7.680e-02 (PASS)
- `Ncond= 14`: gap = 8.882e-16 (FAIL; degenerate HOMO/LUMO)
- `Ncond= 16`: gap = 9.062e-03 (FAIL; below the 5e-2 floor)
- `Ncond= 18`: gap = 1.063e-01 (PASS)
- `Ncond= 20`: gap = 6.193e-01 (PASS)
- `Ncond= 22`: gap = 7.542e-03 (FAIL; below the 5e-2 floor)
- `Ncond= 24`: gap = 4.441e-16 (FAIL; degenerate HOMO/LUMO)

(The script's own `[OK]`/`[REJECT]` marker applies a laxer `gap >=
1e-3` threshold plus an even-`Ncond` parity check per its docstring;
the `>= 5e-2` cutoff above is the stricter spec §1.2 selection
criterion used to pick the pinned value. Note `Ncond=16` and `22`
print `[OK]` from the script but do not clear `>= 5e-2` here.)

This fixture's gap-passing candidates (`gap >= 5e-2`) are `{12, 18,
20}`. `Ncond = 20` was INITIALLY selected across the 4 v3.7 siblings
on the gap criterion alone (cross-fixture physics consistency, and
inside spec §1.2's "expected range 16-24 for closed-shell filling");
the partner-balance analysis below is what moved this fixture off that
value.

### Ncond=20 blocked Phase 2d: partner-balance invariant (v3.7 Task 2d fix)

`soc_pilot_gap.py`'s gap check above is a SCALAR HOMO-LUMO check only.
It does not verify `build_pair_list`'s SOC pairing invariant, which
`general_fij_builder.py` requires independently for every non-self
canonical/partner k-row pair: `n_occ(k) == n_occ(partner(k))`. At the
initially-pinned `Ncond = 20`, this fixture VIOLATES that invariant
(canonical k-row 4, partner row 3: `n_occ_k=2 != n_occ_p=4`), which
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

- `Ncond= 12`: gap PASS, partner-balance **PASS**
- `Ncond= 14`: gap FAIL (Fermi-level degeneracy; partner-balance also
  unevaluable for the same reason)
- `Ncond= 16`: gap FAIL, partner-balance FAIL (k=4/partner=3: `2 != 3`)
- `Ncond= 18`: gap PASS, partner-balance FAIL (k=4/partner=3: `2 != 4`)
- `Ncond= 20`: gap PASS, partner-balance FAIL (k=4/partner=3: `2 != 4`)
- `Ncond= 22`: gap FAIL, partner-balance FAIL (k=4/partner=3: `3 != 4`)
- `Ncond= 24`: gap FAIL (Fermi-level degeneracy)

(Both this table and the Phase 2b table above are now re-sliced from
the same anchor -- the on-disk `flag_fock = true`, `Ncond = 12`
converged solution -- so their gap values agree entry for entry. They
did not agree in earlier revisions of this README, which anchored the
two tables to different self-consistent runs; that discrepancy is
resolved rather than merely documented.)

Candidates satisfying BOTH `gap >= 5e-2` AND partner-balance: `{12}`
-- the ONLY value in the scanned range. `Ncond = 12` is pinned.

Re-verified under `flag_fock = true`: `Ncond = 12` remains the unique
usable value, so enabling the Fock term does not move this fixture's
filling.

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

Pinned: `Ncond = 12` (gap = 7.680e-02 >= 5e-2; partner-balance holds).
Fresh H-wave SCF at `Ncond = 12` with `flag_fock = true` confirms
convergence (38 iterations, rest=8.709e-15,
`Energy_Total=-55.673475845...`) and `phase2_produce_manifest.py`
completes with a self-check-passing 30-entry manifest (composite
`i=2 s=0 j=61 t=1`, `|G_c|=2.1529e-02`). All 30 entries are distinct
active-axis, positive-threshold mutation evaluations. The original M-4
sign flip was structurally invisible on `L_folded=2`, leaving 24
effective entries. An interim fix made M-ship-4 and M-ship-5 both omit
`sub_offset`, so only 27 of the 30 active entries were distinct. The
current M-ship-4 instead halves the named-axis `sub_offset` contribution,
while M-ship-5 omits it. The producer rejects structurally degenerate
active mutations and threshold-policy drift before writing a manifest.
Measured 2026-07-19.

## v3.7 gate role

Phase 6 seven-gate run: `bash tests/validation/uhfk_mvmc_pairproduct/run.sh case_soc_rashba_3d_sub_apbc_xyz` produces the G0-writer-check + G1 + G2a-emitted-F + G2a-in-memory-A + G2b + G3 + G4 records at 1e-6 tolerance for the SOC + APBC(xyz) + SubShape[2,2,2] triple validation.

**G4 guard fix (v3.7 Task 2d)**: while re-deriving this fixture's
manifest at the new `Ncond = 12`, `scripts/soc_apbc_topology_guard.py`
was also generalized to check `sub_offset` differs along every ACTIVE
APBC axis (derived from `theta_radians`) rather than a hardcoded x
check -- see `case_soc_rashba_3d_sub_apbc_yz/README.md`'s "G4 guard
fix" note for the full root cause (a y/z-only fixture, not this one,
is what exposed it). This fixture's own composite already differs
along x (an active axis here), so it was unaffected by the bug either
way; re-verified G4 PASS after the fix.

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
