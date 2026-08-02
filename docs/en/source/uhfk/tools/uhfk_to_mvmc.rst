.. highlight:: none

uhfk_to_mvmc.py: UHFk to mVMC PairProduct bridge
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``tools/uhfk_to_mvmc.py`` converts the SCF result of H-wave UHFk into the
mVMC ``InOrbital`` / ``InOrbitalAntiParallel`` / ``InOrbitalGeneral``
initial wave function file (``zqp_orbital_uhfk.dat``). The PairProduct
state of mVMC can thus be initialized from the Slater matrix of an H-wave
UHF solution.

For the conventions and the construction, see
:doc:`../../algorithm/uhfk_to_mvmc`. This section describes usage.

Supported cases
^^^^^^^^^^^^^^^

Orbitals and lattice
""""""""""""""""""""

- Single orbital (``norb_orig = 1``).
- Sublattice folding supports any ``SubShape`` that divides ``CellShape``
  in every direction. If ``SubShape`` is omitted it falls back to
  ``CellShape``, matching ``_init_lattice`` in ``uhfk.py``.

Conversion paths
""""""""""""""""

There are three paths, selected automatically by the CLI according to the
nature of the SCF solution (see :ref:`uhfk_to_mvmc_dispatch_en`).

AntiParallel path
   Handles Sz-fixed solutions (``2Sz = 0``). Consumes the 3-column or
   4-column form of ``orbitalidx.def``.

General path
   Handles Sz-fixed solutions with ``2Sz != 0``, and Sz-free solutions
   without a mixed block. Consumes the 6-column form of
   ``orbitalidx_general.def``. The user has to generate
   ``orbitalidx_general.def`` with StdFace, setting ``2Sz`` in ``stan.in``
   to a non-zero value for the ``2Sz != 0`` case.

General-SOC path
   Handles Sz-free solutions with a mixed block, that is, solutions
   obtained with spin-orbit coupling enabled
   (``enable_spin_orbital = true``). Supports Zeeman, Rashba, Dresselhaus,
   and general :math:`\sigma_x` / :math:`\sigma_y` type one-body coupling.
   Emits the Sz-non-conserving ``F[up_i, down_j]`` and
   ``F[down_i, up_j]``.

For spin-imbalanced Slater states, same-spin pair components
(``F[up, up]``, ``F[down, down]``) are supported by emitting same-spin
excess pairs from the canonical ``(k, partner(k))`` block.

Boundary conditions
"""""""""""""""""""

Both periodic and antiperiodic boundary conditions are supported
deterministically.

Only when spin-orbit coupling, an antiperiodic direction, and
``SubShape > [1, 1, 1]`` are combined is the set of lattice shapes and
directions restricted to the following.

==============================  ==========================================
Lattice shape                   Supported antiperiodic directions
==============================  ==========================================
``CellShape = [6, 4, 1]``       one direction (x, y, or z)
``SubShape  = [2, 2, 1]``
``CellShape = [4, 4, 4]``       several directions (xy, xz, yz, xyz)
``SubShape  = [2, 2, 2]``
==============================  ==========================================

The restriction does not apply, and the run is supported unconditionally,
in any of these cases:

- spin-orbit coupling is not used,
- ``SubShape = [1, 1, 1]``,
- spin-orbit coupling is used but every direction is periodic.

A combination that falls under none of these and is not in the table is
rejected before dispatch. The error message lists the supported
combinations. Adding a new combination requires a fixture and gate
validation first.

Other properties
""""""""""""""""

- The Slater state is projected to T = 0. A finite-temperature SCF that
  leaves fractional occupation near the Fermi level is rejected, with a
  prompt to recompute at a smaller ``T``.
- A small uniform random number is added to ``params[idx]`` to avoid the
  singularity of the mVMC Pfaffian Slater evaluation for a rank-deficient
  :math:`F`. The default amplitude is ``1e-8`` and can be changed with
  ``--epsilon-noise``. For ``ComplexType 1`` it is added to both the real
  and imaginary parts, for ``ComplexType 0`` to the real part only. This
  is the same technique used by ComplexUHF in mVMC itself
  (``mVMC-1.4.0/src/ComplexUHF/output.c:274``).

Limitations
^^^^^^^^^^^

- Multiple orbitals (``norb_orig > 1``) are out of scope.
- Two-body Sz-non-conserving interactions (spin-flip Coulomb, Hund
  coupling, pair hopping) are out of scope. The only supported two-body
  term is CoulombIntra (on-site :math:`U`).
- The antiperiodic twist is verified only for :math:`\theta` whose
  components are in :math:`\{0, \pi\}`. A general complex twist is out of
  scope.
- Combinations of spin-orbit coupling, antiperiodic boundaries, and
  ``SubShape > [1, 1, 1]`` outside the table above are rejected.

Workflow
^^^^^^^^

1. Generate the mVMC input set (``orbitalidx.def`` and so on) with
   StdFace. Set ``CalcMode = 2`` and ``Lsub`` according to the sublattice
   translational symmetry. For antiperiodic boundaries set
   ``phase0 = 180.0``.
2. Run the H-wave UHFk SCF, requesting ``occupation.npz`` in addition to
   ``eigen.npz`` under ``[file.output]``::

       [file.output]
         path_to_output = "output"
         eigen          = "eigen.npz"
         green          = "green.npz"
         occupation     = "occupation.npz"
         onebodyg       = "greenone.dat"

3. Run the bridge::

       python tools/uhfk_to_mvmc.py \
           --input        input.toml \
           --eigen        output/eigen.npz \
           --occupation   output/occupation.npz \
           --geometry     geometry_uhf.dat \
           --orbitalidx   mvmc_inputs/orbitalidx.def \
           --output       mvmc_inputs/zqp_orbital_uhfk.dat \
           --check-density \
           --onebodyg-uhf output/greenone.dat

4. Register the output file in the mVMC ``namelist.def``. Use the key that
   matches the path taken: ``InOrbital``, ``InOrbitalAntiParallel``, or
   ``InOrbitalGeneral``::

       InOrbitalGeneral zqp_orbital_uhfk.dat

   mVMC then initializes the PairProduct parameters from this file.

.. _uhfk_to_mvmc_dispatch_en:

Path selection
^^^^^^^^^^^^^^

The CLI parses ``--orbitalidx`` and ``input.toml`` and selects the path
from three values.

``is_antiparallel_metadata``
   Whether the metadata of ``orbitalidx.def`` meets the AntiParallel
   requirements.
``orbitalidx_format``
   The column form of ``orbitalidx.def`` (``antiparallel`` or
   ``general``).
``is_soc_mode``
   ``enable_spin_orbital`` in ``input.toml``.

The correspondence is as follows.

==========  ==================  ==========  ===========================
metadata    format              SOC         Path
==========  ==================  ==========  ===========================
``True``    ``antiparallel``    ``False``   AntiParallel
``True``    ``general``         ``False``   forced-General
``False``   ``general``         ``False``   General
``False``   ``antiparallel``    ``False``   rejected
any         ``general``         ``True``    General-SOC
any         ``antiparallel``    ``True``    rejected (6 columns needed)
==========  ==================  ==========  ===========================

The columns correspond, in order, to ``is_antiparallel_metadata``,
``orbitalidx_format``, and ``is_soc_mode``.

``is_antiparallel_metadata`` is ``True`` only when **all** of the
following hold.

- ``2Sz`` in ``input.toml`` is explicitly 0.
- ``N_up == N_down``.
- The values of ``column_spin`` are within ``{0, 1}``.
- ``column_mu_group`` has exactly 2 unique values.
- ``column_spin`` and ``column_mu_group`` are in bijection.

The forced-General path is the branch taken when input carrying
AntiParallel metadata is given in the 6-column form. If the occupation
set satisfies the ``(k_row, local_band)`` pair-closure, ``F[up, down]``
reproduces the :math:`F` of the AntiParallel path to :math:`10^{-12}`.
Otherwise, for instance with an up-up excess arising from spin canting, a
warning is issued, but the output is still a valid ``InOrbitalGeneral``
state. That state cannot be represented on the AntiParallel path.

``(False, antiparallel)`` is rejected because a solution that does not
meet the AntiParallel requirements cannot be represented in the 3-column
or 4-column form. Rerun StdFace with a non-zero ``2Sz`` or a
Zeeman-driven Sz-free setting and regenerate ``orbitalidx_general.def``.

BoundaryCondition contract
^^^^^^^^^^^^^^^^^^^^^^^^^^

The bridge normalizes ``BoundaryCondition`` by delegating to H-wave's
shared helper ``normalize_boundary_condition``. The accepted forms are
below; they are case-insensitive and surrounding whitespace is stripped.

- Periodic: ``"p"``, ``"periodic"``
- Antiperiodic: ``"ap"``, ``"antiperiodic"``

Any other string raises ``ValueError`` before dispatch. There is no
fallback that silently interprets an unknown value. If the key is
omitted, all directions default to periodic.

The ``twist_offset`` stored in ``eigen.npz`` is also checked against
``BoundaryCondition`` in ``input.toml``, so that a stale pairing of input
and eigen data is rejected.

trans.def emission
^^^^^^^^^^^^^^^^^^

When spin-orbit coupling is used, the ``trans.def`` produced by mVMC's
``vmcdry.out`` keeps only the spin-diagonal entries and silently drops the
:math:`s \neq t` terms. To fill that gap, the bridge reads H-wave's
``Transfer.dat`` and emits ``trans.def`` in the ``(i, s, j, t, re, im)``
form. For the mapping convention see
:doc:`../../algorithm/uhfk_to_mvmc`.

The required CLI flags are ``--transfer`` (read ``Transfer.dat``) and
``--emit-trans`` (write ``trans.def``).

Combining spin-orbit coupling with ``SubShape > [1, 1, 1]`` additionally
requires ``--emit-orbitalidx``. Under a folded lattice the
``orbitalidxgen.def`` generated by StdFace cannot represent
Sz-non-conserving pair classes and over-groups them. With this flag the
bridge emits an ``orbitalidx_general.def`` that does not merge classes.

Class consistency check
^^^^^^^^^^^^^^^^^^^^^^^

On the General path, ``aggregate_general_orbital_params`` checks,
**before** averaging, that the signed :math:`F` components assigned to
each class of ``orbitalidx_general.def`` agree within
``class_consistency_tol`` (default :math:`10^{-8}`). If they do not, it
raises ``ClassInconsistencyError`` together with the offending index and
the observed maximum residual.

This check prevents differing values from being averaged silently when the
Slater state does not respect the symmetry assumed by the StdFace-generated
classes, for example when a symmetric Hamiltonian yields a UHF ground state
with spontaneously broken symmetry.

Density matrix check
^^^^^^^^^^^^^^^^^^^^

With ``--check-density``, a one-body density matrix is built from the
converted wave function and compared element-wise against the H-wave
result. The tolerance is :math:`10^{-10}` and exceeding it is fatal. A
mismatch indicates an error in H-wave's boundary-condition handling, in the
bridge, or in the geometry assumptions.

The reference for the comparison depends on the case.

- Normally the density matrix built from the :math:`(k, -k)` pair
  construction is compared against H-wave's ``greenone.dat`` in the
  physical basis.
- When spin-orbit coupling is combined with ``SubShape > [1, 1, 1]``, the
  check switches to ``compare_against_green_sublattice``. H-wave's
  ``green_sublattice`` is lifted to the physical basis with ``gauge_lift``
  and compared element-wise against :math:`\overline{A} A^{T}` built from
  the Slater matrix that is actually emitted. The folded path of
  ``greenone.dat`` is not used as the reference for this combination.

Because the latter takes the emitted output itself as its subject, a
regression on the emission path cannot pass undetected.

Validation
^^^^^^^^^^

The conversion is validated by seven gates. Each gate is an independent
comparison, and the conversion is judged sound only when all of them pass.

G0-writer-check
   The :math:`F` obtained with the rank-lift noise disabled agrees with the
   aggregated ``(mapping, params)``. This checks the soundness of the
   writing step itself. Tolerance :math:`10^{-10}`.
G1
   The density built from the emitted Slater matrix agrees with
   ``green_sublattice`` lifted to the physical basis by ``gauge_lift``.
   Tolerance :math:`10^{-10}`.
G2a-emitted-F
   The density obtained by skew-SVD projection of the emitted :math:`F`
   agrees with the one-body Green function of ComplexUHF. Tolerance
   :math:`10^{-6}`.
G2a-in-memory-A
   The density built from the in-memory Slater matrix agrees with
   ComplexUHF. Tolerance :math:`10^{-6}`.
G2b
   ``green_sublattice`` lifted by ``gauge_lift`` agrees with ComplexUHF.
   Tolerance :math:`10^{-6}`.
G3
   The relative difference between mVMC's :math:`\langle H \rangle` and
   H-wave's ``Energy_Total`` is at most 1 %.
G4
   The composite element is preserved on the current SCF, and applying a
   mutation produces a difference exceeding the threshold
   :math:`T_M = \max(10^{-5},\, 0.10 |G_\mathrm{base}|)`. This confirms
   that the conversion depends on the topology of the lattice.

The three G2 gates require not only agreement but **contraction**.
ComplexUHF is seeded from H-wave's converged density, so if the seed
already met the tolerance, a solver that returned the seed unchanged would
pass. Each G2 therefore requires both that the seed start at least ten
times the tolerance away from the reference and that the converged result
fall inside the tolerance, and it records both values and the contraction
ratio. Non-finite values are rejected before the comparison, since a
comparison against NaN makes both bound tests false and would otherwise
pass silently.

.. note::

   The H-wave run used for validation has to set ``flag_fock = true``. The
   reference solver ``ComplexUHF`` fixes the on-site exchange term at
   compile time (``#define Fock 1`` in ``src/ComplexUHF/include/Def.h``)
   and has no runtime switch. A solution obtained with
   ``flag_fock = false`` is a stationary point of a different mean-field
   functional and agrees with no fixed point of ComplexUHF.

The gates are run with::

    bash tests/validation/uhfk_mvmc_pairproduct/run.sh <case name>

Each gate prints a line beginning with ``<gate name> PASS mode=...``.

Exit codes
^^^^^^^^^^

- ``0``: success.
- ``2``: a fail-fast guard rejected the input. This covers out-of-scope
  modes, an Sz-free SCF, residual fractional occupation at finite
  temperature, and inconsistencies in ``orbitalidx.def``, the geometry, or
  the boundary conditions. The name of the failing check is written to
  stderr.
- ``3``: ``--check-density`` found a mismatch exceeding the tolerance.
