# case_2d_hubbard_sublattice — 2D SubShape=[2,2,1] mVMC E2E fixture

2D square Hubbard 4x4 lattice, PBC in both directions, t=1, U=1, Ne=2
closed shell (Ne_per_spin=1, only k=(0,0) occupied to avoid Fermi
degeneracy). Folded lattice 2x2x1 = 4 folded cells with 4-site
sublattice (subvol=4). StdFace `Lsub=2 Wsub=2` on the mVMC side.

## Status (v2.1)

Runnable end-to-end after the positive-Bloch sign-flip in
`unfold_amplitude_columns` (v2.1, spec §3.3 revision).

Reference run inside the sandboxed build environment:

```
bash tests/validation/uhfk_mvmc_pairproduct/run.sh case_2d_hubbard_sublattice
```

Latest reference (NVMCSample=10000, epsilon_noise=1e-8):

- H-wave UHFk total energy: **-7.938**
- mVMC initial-mean energy: **-7.941**
- relative delta: **~0.05 %**  (within the 1 % rel tolerance)
