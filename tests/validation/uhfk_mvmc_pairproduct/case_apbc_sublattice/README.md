# case_apbc_sublattice — 1D SubShape=[2,1,1] mVMC E2E fixture

L=8 Hubbard APBC Ncond=4, folded lattice L_folded=4 with 2-site
sublattice. StdFace `Lsub=2` on the mVMC side imposes matching
sublattice translation symmetry on `orbitalidx.def`.

## Status (v2.1)

Runnable end-to-end after the positive-Bloch sign-flip in
`unfold_amplitude_columns` (v2.1, spec §3.3 revision). The bridge's
density check now agrees element-wise with `greenone.dat` for
SubShape > [1, 1, 1] (verified to 1e-14 in
`test_bridge_F_reconstructs_correct_density_apbc_subshape_2`).

Reference run inside the sandboxed build environment:

```
bash tests/validation/uhfk_mvmc_pairproduct/run.sh case_apbc_sublattice
```

Latest reference (NVMCSample=10000, epsilon_noise=1e-8):

- H-wave UHFk total energy: **-6.391**
- mVMC initial-mean energy: **-6.372**
- relative delta: **~0.29 %**  (within the 1 % rel tolerance)
