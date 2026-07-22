# case_soc_rashba_3d_sub_apbc_xy (v3.7 shipping fixture)

xy-plane Rashba SOC + general complex z-direction spin-mixing hopping
fixture with antiperiodic boundary condition in x and y (xy plane
APBC, z periodic). Common lattice across all 4 v3.7 shipping fixtures.

## Geometry / interaction

- `CellShape` = `[4, 4, 4]`
- `SubShape` = `[2, 2, 2]` (folded BZ `L_folded = [2, 2, 2]`)
- `BoundaryCondition` = `["antiperiodic", "antiperiodic", "periodic"]` (2D APBC in xy plane)
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
`BoundaryCondition` (theta = (pi, pi, 0)):

```
PYTHONPATH=src python tests/validation/uhfk_mvmc_pairproduct/scripts/soc_pilot_gap.py \
    --case tests/validation/uhfk_mvmc_pairproduct/case_soc_rashba_3d_sub_apbc_xy \
    --ncond 12,14,16,18,20,22,24
```

Measured with `flag_fock = true` (see `input.toml`). The Fock term
changes the converged spectrum, so this table supersedes the earlier
Fock-off one.

Both criteria are listed per candidate: the scalar HOMO-LUMO gap, and
the `build_pair_list` partner-balance invariant `n_occ(k) ==
n_occ(partner(k))` that `general_fij_builder.py` requires for every
non-self canonical/partner k-row pair. The gap check alone cannot see
a partner imbalance (see the yz sibling's README for the case that
forced this second criterion into the selection). Partner-balance is a
frozen-spectrum re-slice of this fixture's own converged spectrum
(`step_occupation` + `find_partner_rows` + `build_pair_list`, no SCF
rerun per candidate).

- `Ncond= 12`: gap = 1.192e-02 (FAIL; below the 5e-2 floor); partner-balance FAIL (k=4/partner=2: `1 != 2`)
- `Ncond= 14`: gap = 1.332e-15 (FAIL; degenerate HOMO/LUMO)
- `Ncond= 16`: gap = 3.966e-03 (FAIL; below the 5e-2 floor); partner-balance FAIL (k=4/partner=2: `2 != 3`)
- `Ncond= 18`: gap = 2.480e-01 (PASS); partner-balance FAIL (k=4/partner=2: `2 != 3`)
- `Ncond= 20`: gap = 5.796e-01 (PASS); partner-balance **PASS**
- `Ncond= 22`: gap = 8.882e-16 (FAIL; degenerate HOMO/LUMO)
- `Ncond= 24`: gap = 4.476e-01 (PASS); partner-balance **PASS**

(The script's own `[OK]`/`[REJECT]` marker applies a laxer `gap >=
1e-3` threshold plus an even-`Ncond` parity check per its docstring;
the `>= 5e-2` cutoff above is the stricter spec §1.2 selection
criterion used to pick the pinned value. Note `Ncond=12` and `16`
print `[OK]` from the script but do not clear `>= 5e-2` here.)

Candidates satisfying BOTH `gap >= 5e-2` AND partner-balance: `{20,
24}`. `Ncond = 20` stays pinned: it has the larger gap (5.796e-01 vs
4.476e-01), it matches the xz sibling's pinned filling, and it sits
inside spec §1.2's "expected range 16-24 for closed-shell filling".
The switch to `flag_fock = true` therefore does not move this
fixture's filling.

Pinned: `Ncond = 20` (gap = 5.796e-01 >= 5e-2, partner-balance holds).

## v3.7 gate role

Phase 6 seven-gate run: `bash tests/validation/uhfk_mvmc_pairproduct/run.sh case_soc_rashba_3d_sub_apbc_xy` produces the G0-writer-check + G1 + G2a-emitted-F + G2a-in-memory-A + G2b + G3 + G4 records at 1e-6 tolerance for the SOC + APBC(xy) + SubShape[2,2,2] triple validation.

The regenerated G4 manifest keeps the 30-entry per-direction schema
and selects composite `i=5 s=0 j=42 t=1`,
`|G_c|=3.7700e-02`. Its active x/y axes contribute 20 distinct,
positive-threshold mutation evaluations; the 10 z-axis entries remain
at zero by the inactive-axis policy in spec §4b. The original M-4
sign flip was structurally invisible on `L_folded=2`, leaving 16
effective entries. An interim fix made M-ship-4 and M-ship-5 both omit
`sub_offset`, so only 18 of the 20 active entries were distinct. The
current M-ship-4 instead halves the named-axis `sub_offset` contribution,
while M-ship-5 omits it. The producer rejects structurally degenerate
active mutations and threshold-policy drift before writing a manifest.

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
contract documented in that script and used by v3.6's fixture.
