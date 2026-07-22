# case_soc_rashba_3d_sub_apbc_xy

ComplexUHF mirror fixture for the v3.7 64-site SOC validation case at
`tests/validation/uhfk_mvmc_pairproduct/case_soc_rashba_3d_sub_apbc_xy/`.
It uses `CellShape = [4, 4, 4]`, `SubShape = [2, 2, 2]`, `U = 2`,
`Ncond = 20`, and APBC in x and y.

The committed files were bootstrapped with mVMC 1.4.0 `vmcdry.out` from
the orthorhombic `stan.in` in a scratch directory. Orthorhombic StdFace
maps `phase0`, `phase1`, and `phase2` to x, y, and z, so this case uses
the phase mask `(180, 180, 0)` degrees. At runtime the harness replaces
the generated hopping baseline with the H-wave bridge's `trans.def` and
supplies the SCF-seeded `initial.def`; neither runtime file is part of
this pre-substitution bundle.

`orbitalidx.def` and `orbitalidxpara.def` are syntactically valid dummy
tables retained only for ComplexUHF's post-SCF `MakeOrbitalFile` path.
They are not a physical parametrization of the substituted SOC
Hamiltonian. The resulting `zqp_orbital.dat` is physically meaningless
and must not be reused as an mVMC initial state.

There is intentionally no `greenone.def`: ComplexUHF `cal_cisajs`
writes all `(2 * Nsite)^2 = 16384` one-body rows unconditionally.
