==========================================
Tutorial: FLEX solver
==========================================

This tutorial demonstrates how to use the FLEX
(Fluctuation Exchange Approximation) solver in H-wave.
FLEX extends RPA by using dressed (self-consistent) Green's functions
instead of bare ones, providing a more accurate description of
correlated electron systems.

The sample files for this tutorial are located in
``docs/en/source/flex/sample`` directory.


Overview
----------------------------

The FLEX approximation [1]_ is a self-consistent diagrammatic method
for itinerant electron systems. Unlike RPA, which uses the bare
Green's function :math:`G_0`, FLEX iterates the following self-consistent loop
until convergence:

1. Compute the dressed Green's function :math:`G(\mathbf{k}, i\omega_n)`
   from the Dyson equation.
2. Compute the bare susceptibility :math:`\chi_0(\mathbf{q}, i\nu_m)`
   from the dressed :math:`G`.
3. Decompose the interaction into spin and charge channels.
4. Solve the RPA equations for spin/charge susceptibilities
   :math:`\chi_s` and :math:`\chi_c`.
5. Construct the effective interaction :math:`V_{\mathrm{eff}}`.
6. Compute the self-energy :math:`\Sigma(\mathbf{k}, i\omega_n)` via
   FFT convolution.
7. Check convergence; if not converged, go to step 1.

.. note::

   When the electron number is fixed through ``filling`` / ``Ncond`` (rather
   than a fixed ``mu``), FLEX re-solves the chemical potential :math:`\mu` from
   the *dressed* Green function at every SCF iteration so that the target
   filling is maintained self-consistently as the self-energy grows.  A
   ``FLEX._find_mu_dressed: mu = ...`` line is therefore printed each iteration,
   and the converged :math:`\mu` (and the exact iteration count) differ from a
   calculation that keeps :math:`\mu` fixed at its non-interacting value.  All
   iteration counts and convergence values shown in this tutorial are
   illustrative and may vary slightly with the version and platform.


Theory
----------------------------

Dressed Green's function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The dressed Green's function is obtained from the Dyson equation:

.. math::

   G(\mathbf{k}, i\omega_n)
   = \left[ G_0^{-1}(\mathbf{k}, i\omega_n) - \Sigma(\mathbf{k}, i\omega_n) \right]^{-1}

where :math:`G_0^{-1}(\mathbf{k}, i\omega_n) = i\omega_n + \mu - H_0(\mathbf{k})`
is the inverse bare Green's function.

Spin and charge susceptibilities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The bare susceptibility is computed from the dressed Green's function:

.. math::

   \chi_0(\mathbf{q}, i\nu_m) = -\frac{T}{N_k} \sum_{\mathbf{k}, n}
   G(\mathbf{k}+\mathbf{q}, i\omega_n + i\nu_m)\, G(\mathbf{k}, i\omega_n)

The spin and charge susceptibilities are:

.. math::

   \chi_s = \left[ I - \chi_0 \, U_s \right]^{-1} \chi_0

.. math::

   \chi_c = \left[ I + \chi_0 \, U_c \right]^{-1} \chi_0

where :math:`U_s` and :math:`U_c` are the spin and charge interaction vertices
decomposed from the full interaction Hamiltonian.
For the single-band Hubbard model, :math:`U_s = U_c = U`.

Effective interaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The effective FLEX interaction combines spin and charge fluctuations [1]_:

.. math::

   V_{\mathrm{eff}}(\mathbf{q}, i\nu_m)
   = W \left[ \frac{3}{2}\chi_s
            + \frac{1}{2}\chi_c - \chi_0 \right] W

where :math:`W` is the bare interaction vertex.
:math:`\chi_0` is subtracted **once**: at lowest order
:math:`\chi_s = \chi_c = \chi_0`, so the bracket reduces to :math:`\chi_0` and
:math:`V_{\mathrm{eff}} = W \chi_0 W` (the second-order :math:`U^2` bubble).
This single subtraction removes the double-counted second-order diagram
contained in :math:`\chi_s` and :math:`\chi_c`.

Self-energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The self-energy is computed via convolution in real space and
imaginary time:

.. math::

   \Sigma(\mathbf{r}, \tau) = V_{\mathrm{eff}}(\mathbf{r}, \tau) \cdot G(\mathbf{r}, \tau)

This element-wise (Hadamard) product is efficiently evaluated using FFT.

.. _flex_scope:

Scope of the approximation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

FLEX is a conserving (Baym--Kadanoff) approximation that resums the RPA
particle--hole bubble and ladder series for the self-energy. It does **not**
include the Aslamazov--Larkin (AL) or Maki--Thompson (MT) vertex corrections,
in which the electron couples to *two* fluctuation propagators through a
triangular fermion loop (mode--mode coupling). These higher-order corrections
are outside the FLEX class and are **not** evaluated here. They can matter when
charge/orbital fluctuations driven by two spin fluctuations are important
(e.g. the orbital-fluctuation mechanism of Onari and Kontani [2]_).

In addition, the ``calc_scheme = "reduced"`` scheme decomposes the
interaction via its density--density part for the spin/charge vertices.
Since version 2.0 ``reduced`` is no longer the default: ``calc_scheme =
"auto"`` is, and ``reduced`` is the explicit opt-out that reproduces the
H-wave 1.0.x behaviour (a documented approximation, warned about once at
construction). ``auto`` selects ``general`` whenever ``Exchange`` or
``PairHop`` is present -- these have **no** density--density vertex content
at all, so ``reduced`` cannot even approximate them, and supplying either
under an explicit ``reduced`` raises a ``ValueError`` directing you to
``calc_scheme = "general"`` (in earlier development builds this input was
accepted with a warning while the interaction silently had zero effect).
``auto`` likewise selects ``general`` whenever ``CoulombInter``, ``Hund``,
``Ising``, or the aggregate ``Coulomb`` interaction is declared, because
these carry cross-family vertex content that ``reduced`` cannot represent
exactly.
``PairLift`` is accepted everywhere: its particle-hole vertex is exactly
zero, so omitting it from the susceptibility channels is exact, not an
approximation. Accordingly, in this scheme "FLEX" means *not exact*: it
is the density--density, AL/MT-free fluctuation-exchange level of
approximation.

The alternative ``calc_scheme = "general"`` instead **retains** the full
off-diagonal Kanamori vertices. It is a paramagnetic full-vertex
formulation following Mochizuki, Yanase, and Ogata (MYO) [3]_ (and
corroborated by Takimoto, Hotta, and Ueda (THU) [4]_): the full matrix-form
spin (:math:`\hat{U}^s`) and charge (:math:`\hat{U}^c`) interaction matrices
are built in the MYO convention, the matrix RPA is solved for
:math:`\chi_s`/:math:`\chi_c`, and the fluctuation interaction is assembled
as :math:`V = \tfrac{3}{2}\hat{U}^s[\chi_s - \chi_0]\hat{U}^s
+ \tfrac{1}{2}\hat{U}^c[\chi_c - \chi_0]\hat{U}^c + W^{(2)}`, the ring
series from third order on plus the second-order kernel :math:`W^{(2)}`
selected by ``flex_second_order`` (see below).
Under ``"general"`` the off-diagonal vertices are therefore **not** dropped
and the density--density reduction warning is suppressed. The AL/MT vertex
corrections noted above remain outside the FLEX class even in the
``"general"`` scheme.

.. [1] N. E. Bickers and D. J. Scalapino,
   Ann. Phys. (N.Y.) **193**, 206 (1989).

.. [2] H. Kontani and S. Onari,
   Phys. Rev. Lett. **104**, 157001 (2010).

.. [3] M. Mochizuki, Y. Yanase, and M. Ogata,
   J. Phys. Soc. Jpn. (cond-mat/0407094).

.. [4] T. Takimoto, T. Hotta, and K. Ueda,
   Phys. Rev. B **69**, 104504 (2004); cond-mat/0309575.


.. _flex_second_order_tutorial:

Second order of the effective interaction (``flex_second_order``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Under ``calc_scheme = "general"`` the second-order part :math:`W^{(2)}` of
the effective interaction is built by default from the exact LOCAL
second-order kernel, ``[mode.param] flex_second_order = "local"``. It
carries the complete second order of every on-site interaction term, the
direct skeleton of the off-site density terms and every mixed
on-site/off-site diagram; only the exchange skeleton of two off-site
vertices is left to ``longitudinal_bond_channels = true`` (see
:ref:`flex_second_order_kernel` for the formulas and the classification).

Earlier versions used the Takimoto-Hotta-Ueda expression
:math:`-\tfrac{1}{4}(\hat{U}^s+\hat{U}^c)\chi_0(\hat{U}^s+\hat{U}^c)`,
which is exact at second order only for a single-band Hubbard
interaction. Set

.. code-block:: toml

   [mode.param]
     flex_second_order = "takimoto"

to reproduce results produced with H-wave 2.0.0 and earlier releases;
the option stays available throughout the 2.x series.

**What changes.** Single-band inputs containing only ``CoulombIntra``
agree under both kernel choices (``"local"`` and ``"takimoto"``) to within
numerical round-off, and so do single-band :math:`U + V` inputs run with
the bond-resolved channels (``longitudinal_bond_channels = true``), which
already carried the exact direct :math:`V` second order. Every
multi-orbital on-site interaction (:math:`U'`, ``Hund``, ``Ising``,
``Exchange``, ``PairHop``, ``PairLift``) and every off-site interaction
run without the bond-resolved channels gives a different
``calc_scheme = "general"`` result.

**On-site same-orbital rows (degenerate rows) are now refused.** Under
``"local"`` an on-site same-orbital row of ``CoulombInter``, ``Hund``,
``Ising``, ``Exchange``, ``PairHop`` or ``PairLift`` -- a row with
``rx = ry = rz = 0`` and identical orbital indices -- stops the run at
start-up: such a row is not a two-body term, because it reduces to a
one-body level shift, to an effective ``CoulombIntra``, or to identically
zero. The message gives the equivalent declaration for that type (usually
a ``CoulombIntra`` entry plus a level shift in the transfer file).
Rewrite the interaction file as the message says, or use
``flex_second_order = "takimoto"`` as an immediate workaround.

**Convergence.** The SCF trajectory can change with the kernel. On the
2-orbital self-consistency fixture (on-site :math:`U`, :math:`U'`,
``Hund`` plus off-site :math:`V`, Anderson mixing) both values needed the
same iteration count -- 9 without and 11 with ``flex_hartree_fock =
true`` -- and a full FLEX iteration costs about 1.2x more under
``"local"``. If a run that used to converge now stalls, raise
``IterationMax``, lower ``Mix``, switch to (or keep) Anderson mixing
(``mixing_scheme = "anderson"``; the solver default is ``"linear"``), and
start the loop from :math:`\Sigma = 0` -- that is, omit ``sigma_init``
from the ``[file.input]`` section -- rather than continuing from a seed
written with the other kernel: the solver warns when the seed's recorded
``flex_second_order`` differs from the run's (and notes it when the seed
carries no record, i.e. a reduced-scheme archive or one written before
this key was introduced). An archive written by H-wave 2.0.0 or an earlier
release carries no kernel record at all, so only that informational line is
logged; if warm-starting such a seed stalls under
``"local"``, either restart from :math:`\Sigma = 0` or set
``flex_second_order = "takimoto"`` to match the kernel the seed was
produced with.


.. _flex_near_instability:

Operating near an instability
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A strongly coupled run can approach a magnetic or charge instability:
the Stoner factor climbs above about 0.99, or the ``guards:`` line of
the log reports a bond conditioning minimum (``cond_min``) within an
order of magnitude of ``[mode.param] longitudinal_bond_cond_tol``. In
that regime the self-consistency is stiff, several fixed points can
exist, and the following operating rules apply. Every key named in
this subsection belongs to the ``[mode.param]`` table:

.. code-block:: toml

   [mode.param]
     mixing_scheme = "linear"
     Mix = 0.02
     flex_hartree_fock = true
     flex_guard_policy = "warn"

- **Mix in small linear steps.** Use ``mixing_scheme = "linear"`` with
  ``Mix <= 0.03``. Larger steps overshoot into the instability region,
  where the dressed vertices are enormous and the guard stops the run.
- **Warm-start in small temperature steps.** Walk down in temperature
  and start each run from the converged self-energy of the previous,
  slightly higher temperature (``[file.input] sigma_init``). A cold
  start at the lowest temperature of a sweep can fail where the
  descending ladder converges.
- **Cross-check an Anderson-mixing result by a linear recomputation.**
  Anderson mixing has been observed to land on an unphysical fixed
  point -- one whose pairing eigenvalue saturates instead of growing --
  on the single-band Hubbard model at :math:`U = 8t`. Recompute the
  same point with linear mixing before reporting it.
- **Set** ``flex_guard_policy = "warn"`` **in** ``[mode.param]`` **only
  for a transient excursion.** It lets an iteration whose equal-time
  density symmetry (``flex_hf_density_tol``) or bond conditioning
  (``longitudinal_bond_cond_tol``) violates its tolerance continue
  instead of ending the run. An exactly singular solve, or non-finite
  numbers, still end the run under either policy: the policy tolerates
  a nearly singular denominator, not a singular one. The final state
  is checked regardless of the policy, so the run is refused when the
  equal-time density of the final state deviates by more than
  ``flex_hf_density_tol``, and also when the last map needed the
  policy for a bond-guard violation -- that refusal reports the
  conditioning minima of the last map. ``flex_guard_violations`` in
  the outputs counts the violations that were let through: one per
  iteration for the density, and one per offending (channel,
  frequency, q-point) finding for the bond guard, so a single
  iteration can contribute several.

Read the ``guards:`` line the solver logs for each iteration and
confirm that the violation decays: the density deviation must fall
back below ``flex_hf_density_tol``, ``cond_min`` must rise back above
``longitudinal_bond_cond_tol``, and both must stay there.

.. code-block:: text

      guards: density hermitian 3.2e-10 (tol 1.0e-08)  bond cond_min spin 4.1e-03 charge 2.7e-01 (tol 1.0e-03)

The line reports the relative Frobenius deviation of the Hermitian
symmetry of the equal-time density against ``flex_hf_density_tol``
and, with ``longitudinal_bond_channels = true``, the smallest
conditioning score of the spin and of the charge RPA denominator of
that iteration's map against ``longitudinal_bond_cond_tol``.

In a sweep, keep the default ``flex_guard_policy = "refuse"`` for the
automated runs and turn ``"warn"`` on for a single point you have
already diagnosed; a post-processing script should check that
``flex_guard_violations`` is ``0`` before accepting a converged result.
Since the count is per event, one stiff iteration of the bond guard
can contribute many of them -- read the ``guards:`` lines to see how
many iterations were actually affected.

The tolerances themselves (``flex_hf_density_tol``,
``longitudinal_bond_cond_tol``) are knobs for a deliberately stiff
study, not a way to silence a guard: lowering a floor lets a run
continue with numbers that are dominated by amplified round-off.
Raising ``flex_hf_density_tol`` or lowering
``longitudinal_bond_cond_tol`` is therefore not the remedy for a
final-state refusal -- a smaller ``Mix`` and finer temperature steps
are.

Sample 1: Single-orbital Hubbard model
-----------------------------------------

Model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first sample is a **single-orbital Hubbard model** on a
two-dimensional square lattice at half filling (:math:`n = 1`).

.. math::

   H = -\sum_{\langle i,j \rangle, \sigma} t_{ij}\,
       c^\dagger_{i\sigma} c_{j\sigma}
     + U \sum_i n_{i\uparrow} n_{i\downarrow}

with nearest-neighbor hopping :math:`t = 1.0`,
next-nearest-neighbor hopping :math:`t' = 0.5`,
on-site Coulomb repulsion :math:`U = 4.0`,
and temperature :math:`T = 0.5`.

This model is known to exhibit strong antiferromagnetic (AF) spin
fluctuations with a peak at :math:`\mathbf{Q} = (\pi, \pi)`.

Prepare input files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The sample files are in ``docs/en/source/flex/sample/1orb/``.

**Parameter file** (``input.toml``):

.. literalinclude:: ../sample/1orb/input.toml

Key parameters:

- ``mode = "FLEX"``: Selects the FLEX solver.
- ``T = 0.5``: Temperature.
- ``CellShape = [8, 8, 1]``: 8 x 8 k-point mesh for a 2D system.
- ``Nmat = 64``: Number of Matsubara frequencies.
- ``filling = 0.5``: target electron number per site (half filling). Specifying
  ``filling`` (or ``Ncond``) makes FLEX re-solve the chemical potential
  :math:`\mu` from the dressed Green's function at every SCF iteration so the
  filling is conserved as the self-energy grows; specifying ``mu`` instead holds
  it fixed.
- ``coeff_tail = 1.0`` (optional): high-frequency tail-acceleration coefficient
  for the Matsubara sums. ``coeff_tail = 1`` matches the exact
  :math:`1/(i\omega_n)` coefficient of :math:`G` (unitarity), so it accelerates
  convergence in ``Nmat`` without biasing the result. Also supported by the RPA
  solver. Ignored (unnecessary) when ``matsubara_basis = "ir"``.
- ``IterationMax = 100``: Maximum number of SCF iterations.
- ``Mix = 0.2``: Mixing parameter for self-energy update
  (:math:`\Sigma_{\mathrm{new}} = (1 - \alpha)\Sigma_{\mathrm{old}} + \alpha\Sigma_{\mathrm{calc}}`).
- ``EPS = 6``: Convergence criterion :math:`10^{-6}`.

**Geometry** (``geom.dat``):

.. literalinclude:: ../sample/1orb/geom.dat

A single orbital at the origin.

**Transfer integrals** (``transfer.dat``):

.. literalinclude:: ../sample/1orb/transfer.dat

Nearest-neighbor (:math:`t = 1.0`) and next-nearest-neighbor
(:math:`t' = 0.5`) hopping on the square lattice.

**On-site interaction** (``coulombintra.dat``):

.. literalinclude:: ../sample/1orb/coulombintra.dat

On-site Coulomb repulsion :math:`U = 4.0`.


Run the calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    $ cd docs/en/source/flex/sample/1orb
    $ hwave input.toml

The output log shows the SCF convergence:

.. code-block:: text

    FLEX iteration 1/100
    FLEX._find_mu_dressed: mu = -0.398893
      convergence: |dSigma|/|Sigma| = 1.000e+00
    FLEX iteration 2/100
    FLEX._find_mu_dressed: mu = -0.291966
      convergence: |dSigma|/|Sigma| = 9.876e-01
    ...
    FLEX iteration 64/100
    FLEX._find_mu_dressed: mu = -0.249146
      convergence: |dSigma|/|Sigma| = 9.241e-07
    FLEX converged after 64 iterations


Results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After convergence, the solver produces the following output files
in the ``output`` directory:

- ``chi0q.npz``: Bare susceptibility :math:`\chi_0(\mathbf{q}, i\nu_m)`
- ``chiq_s.npz``: Spin susceptibility :math:`\chi_s(\mathbf{q}, i\nu_m)`
- ``chiq_c.npz``: Charge susceptibility :math:`\chi_c(\mathbf{q}, i\nu_m)`
- ``chiq.npz``: Combined susceptibility file
- ``sigma.npz``: Self-energy :math:`\Sigma(\mathbf{k}, i\omega_n)` (with
  ``flex_hartree_fock = true`` also its two components ``sigma_static`` and
  ``sigma_fluct``, the marker ``sigma_convention = "split"`` and the
  convergence provenance; see :ref:`flex_bond_hf_tutorial`). Under
  ``calc_scheme = "general"`` it also records ``flex_second_order`` and
  ``flex_second_order_schema``.
- ``green.npz``: Dressed Green's function :math:`G(\mathbf{k}, i\omega_n)`
  (with the same two ``flex_second_order`` fields under ``calc_scheme =
  "general"``)
- ``energy.dat``: Text file with the particle number ``NCond``, spin
  ``Sz``, and the converged ``ChemicalPotential`` :math:`\mu`.

.. note::

   The ``energy.dat`` output (enabled by the ``energy`` key in
   ``[file.output]``) is written from the final dressed Green function.  In
   the fixed-:math:`\mu` mode (specify ``mu`` instead of ``filling`` /
   ``Ncond``) its ``NCond`` line gives the particle number at that
   :math:`\mu`, so running the solver at several fixed :math:`\mu` values
   traces the :math:`\mu`-:math:`N` relation.  ``Sz`` is zero for a
   paramagnetic (spin-free) calculation and non-zero only for
   spin-dependent (spin-diagonal / spinful) runs.

.. _flex_sigma_init:

Warm-starting the SCF loop (``sigma_init``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By default FLEX starts the self-consistency loop from :math:`\Sigma = 0`.
Setting ``sigma_init`` in ``[file.input]`` to a ``sigma.npz`` written by an
earlier FLEX run instead seeds the loop from that self-energy:

.. warning::

   A multi-orbital ``sigma.npz`` produced by ``calc_scheme = "general"``
   before the orbital-pair transpose fix is wrong off the orbital diagonal.
   The corrected solver converges to the corrected fixed point regardless of
   the seed, so such a file does not poison the result, but it is a poor warm
   start -- it can slow convergence or steer the iteration towards a different
   solution branch, which defeats the purpose of seeding. Prefer regenerating
   it, and note that the run's other outputs must be regenerated in any case --
   see :ref:`the migration warning <flex_general_transpose_fix>`.

.. code-block:: toml

   [file.input]
     sigma_init = "sigma.npz"

This is often decisive near a magnetic instability (low temperature, strong
spin fluctuations), where the :math:`\Sigma = 0` transient makes the SCF
*oscillate* (the residual ``|dSigma|/|Sigma|`` stalls around 1 instead of
decreasing) and the run hits ``IterationMax`` without converging. Starting from
a converged neighbouring solution -- for example, stepping the temperature down
and feeding each run the previous (higher-:math:`T`) ``sigma.npz`` -- begins the
iteration near the fixed point and avoids the oscillation. The seed must have
the same ``CellShape`` and ``Nmat`` as the current run (both are fail-fast
errors: ``sigma.npz`` records its ``CellShape``, so even a same-volume
aspect-ratio change like ``[2,8,1]`` vs ``[4,4,1]`` is caught), so keep
``Nmat`` and ``CellShape`` fixed across a continuation sweep.

A ``calc_scheme = "general"`` seed also records which second-order kernel
produced it. Seeding across kernels is allowed, and the solver warns when
the seed's ``flex_second_order`` differs from the current run's (and notes
it when the seed carries no such record, i.e. a reduced-scheme archive or
one written before this key was introduced -- H-wave 2.0.0 and earlier). If
such a run stalls, restart it from :math:`\Sigma = 0` -- see
:ref:`flex_second_order_tutorial`.

.. note::

   The ``sigma_init`` path is resolved relative to ``[file.input]
   path_to_input``, while the previous run wrote its ``sigma.npz`` under
   ``[file.output] path_to_output``. In a sweep, either copy the previous
   ``sigma.npz`` into the input directory, or point ``sigma_init`` at the
   previous output directory with a relative path, e.g.::

      [file.input]
        path_to_input = "."
        sigma_init = "run_T0.50/output/sigma.npz"

.. _flex_bond_hf_tutorial:

Hartree-Fock renormalisation and the bond-resolved channel (experimental)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``flex_hartree_fock = true`` adds the self-consistent Hartree-Fock
self-energy of every interaction term to the FLEX loop, and
``longitudinal_bond_channels = true`` (which requires it) builds the
fluctuation part on the bond-resolved pair basis so that the exchange
crossing of an off-site ``CoulombInter`` / ``Hund`` / ``Ising`` enters
self-consistently (see :ref:`flex_bond_hf` for the equations and the domain).
Both are ``calc_scheme = "general"``, spin-free, uniform-grid options, and
both run on the GPU (``gpu = true``, via CuPy) as well as on the CPU: on
its own, ``flex_hartree_fock = true`` uses the ordinary FLEX GPU path (see
the ``gpu`` option below); with ``longitudinal_bond_channels = true``
(which requires ``flex_hartree_fock = true``) the dressing, the effective
interaction and the self-energy transport run on the GPU instead, while
the bond arrays stay in host memory and are transferred one frequency
batch at a time -- the batch is chosen against both the host cap
(``longitudinal_bond_memory_cap_gb``) and the free device memory
(``longitudinal_bond_freq_batch`` overrides both and is refused when it
exceeds either). Results agree with the CPU path to round-off; the
outputs record ``longitudinal_bond_device`` and ``longitudinal_bond_nb``.
Once a run has been checked with ``longitudinal_bond_guard_freqs = "all"``,
production GPU runs can switch to ``"static"``: on a GPU the full guard costs
about 2.5 times the static one, versus about 20% on the CPU.
A complete input for a single-band square lattice with an on-site ``U`` and a
nearest-neighbour ``V`` (the interaction files follow the Wannier90-style
format of the :ref:`interaction input <Ch:Config_rpa>`; ``coulombinter.dat``
lists the four bonds ``(+-1, 0, 0)``, ``(0, +-1, 0)``):

.. code-block:: toml

   [mode]
     mode = "FLEX"
     calc_scheme = "general"
   [mode.param]
     T = 0.02
     filling = 0.35
     CellShape = [16, 16, 1]
     Nmat = 2048
     mixing_scheme = "anderson"
     anderson_depth = 8
     Mix = 0.2
     EPS = 8
     flex_hartree_fock = true
     longitudinal_bond_channels = true
     # longitudinal_bond_memory_cap_gb = 8.0  # refuse before running if the estimate exceeds it
     # longitudinal_bond_freq_batch = 64      # override the automatic frequency batch
     # longitudinal_bond_guard_freqs = "all"  # "static" checks only the zero frequency by SVD
     # longitudinal_bond_output_full = true   # dynamic chi_s_w / chi_c_w archive (doubles the memory)
   [file.input]
     path_to_input = "."
   [file.input.interaction]
     path_to_input = "."
     Geometry = "geom.dat"
     Transfer = "transfer.dat"
     CoulombIntra = "coulombintra.dat"
     CoulombInter = "coulombinter.dat"
   [file.output]
     path_to_output = "output"
     sigma = "sigma"
     green = "green"
     chiq = "chiq"         # the static longitudinal_bond_* keys go here
     # longitudinal_bond = "longitudinal_bond.npz"   # the dynamic archive (output_full only)
     energy = "energy.dat"

Set ``[file.output] chiq`` in such runs: without it the static bond keys of
BOTH channels are written into ``chiq_s.npz``. Output archive names without
a ``.npz`` suffix get one appended, as for every other ``.npz`` output.

The run is accepted only when the three residuals of the split state stay
below ``EPS`` for three consecutive iterations; the log prints them with a
``[pass k/3]`` counter, and every archive records ``scf_converged`` and the
last residuals (see the output reference). The ``hf_density_error`` field
reports how well the Hartree-Fock density closes on the target filling.

*Starting from a UHFk mean field.* Pass the UHFk ``trans_mod`` archive as
today (``[file.input] trans_mod = "trans_mod.npz"``): with
``flex_hartree_fock = true`` the mean field is NOT folded into the band but
becomes the initial static self-energy, so the first iteration starts from
the UHFk solution and the Hartree-Fock term is never counted twice. No
converter is needed for this workflow. The UHFk solution must be
paramagnetic (``2Sz = 0``, equal spin blocks, no spin mixing): a
magnetically ordered ``trans_mod`` is refused, because the Hartree-Fock
FLEX is spin-free.

*Restarting from a previous run.* A ``sigma.npz`` written with
``flex_hartree_fock = true`` carries the split components and seeds a new
run exactly as stored (``[file.input] sigma_init``). Remove (or comment
out) ``trans_mod`` and ``green_init`` from ``[file.input]`` in that case:
the archive already contains the static part, and a mean field given at
the same time is refused. A ``sigma.npz`` from a run without
``flex_hartree_fock`` (or a hand-made total self-energy) is refused as a
seed and must be converted; use ``--zero-static`` when that run had no
mean field (its static part was absorbed by the chemical potential), and
``--uhfk-trans-mod`` when it was started from a UHFk band modification::

   hwave_sigma_split total.npz seed.npz --zero-static
   hwave_sigma_split total.npz seed.npz --static static.npz
   hwave_sigma_split total.npz seed.npz --uhfk-trans-mod trans_mod.npz --bare-transfer transfer.dat

The last form builds the static correction from a native UHFk ``trans_mod``
archive and the run's bare Transfer input (the archive is Fourier-transformed
with the solver's convention, the bare transfer is subtracted and the result
is reduced to the spin-free block). ``--force`` overwrites an existing output.

**Spin susceptibility** :math:`\chi_s(\mathbf{q}, i\nu_0)`:

.. figure:: ../sample/1orb/chi_s.png
   :width: 60%
   :align: center

   Static spin susceptibility :math:`\chi_s(\mathbf{q})` of the
   single-orbital Hubbard model at half filling.
   The peak at :math:`\mathbf{Q} = (\pi, \pi)` indicates strong
   antiferromagnetic spin fluctuations, consistent with the nesting
   of the Fermi surface.

**Self-energy** :math:`\mathrm{Im}\,\Sigma(\mathbf{k}, i\omega_0)`:

.. figure:: ../sample/1orb/sigma_kspace.png
   :width: 60%
   :align: center

   Imaginary part of the self-energy at the lowest Matsubara frequency.
   The k-dependence reflects the scattering of quasiparticles by
   spin fluctuations, with stronger damping near the antiferromagnetic
   hot spots.

**Self-energy vs Matsubara frequency**:

.. figure:: ../sample/1orb/sigma_matsubara.png
   :width: 80%
   :align: center

   Frequency dependence of the self-energy at selected k-points.
   The imaginary part :math:`\mathrm{Im}\,\Sigma(i\omega_n) < 0`
   shows quasiparticle damping, while the
   :math:`1/\omega_n` tail at high frequencies indicates
   Fermi liquid behavior.


.. _flex_bond_pairing_tutorial:

Pairing eigenvalues with the bond-resolved vertex
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once a Hartree-Fock FLEX run with the bond-resolved channel
(``longitudinal_bond_channels = true``, above) has produced a dressed
Green function, the linearized Eliashberg equation can be solved with the
SAME bond-resolved, frequency-resolved pairing vertex the self-energy used
-- rather than the on-site ``chiq_s.npz`` / ``chiq_c.npz`` vertex of
:ref:`the dynamic-frequency section <sc_dynamic_frequency>`. That is all
the in-process entry needs. The dynamic bond archive
``longitudinal_bond.npz`` (schema 2), which is written only with
``longitudinal_bond_output_full = true``, is needed solely by the
post-processing entry through ``hwave_sc``. Two entries reach this kernel
and build it from the same arrays, so they agree to round-off on the
uniform Matsubara grid:

- **In-process**, at the end of the FLEX solve itself: set
  ``[mode.param] longitudinal_bond_pairing = "singlet" | "triplet" |
  "both"`` in the FLEX input. Requires ``longitudinal_bond_channels =
  true`` and ``IterationMax >= 1``. An optional ``[eliashberg]`` table in
  the SAME input file configures the solver (every key defaults exactly as
  it does for ``hwave_sc``); ``pairing_type`` is refused there (the
  channel is chosen by ``longitudinal_bond_pairing``), and ``gpu`` /
  ``gpu_required`` / ``bond_channels`` / ``chi0q_mode`` / ``frequency`` /
  ``flex_bond_archive`` / ``bond_green`` / ``bond_max_shells`` are ignored
  with an INFO line -- the pairing step runs on the FLEX solve's own
  backend and data, not on a re-read archive. Writes
  ``eliashberg_bond_<type>.npz``, ``gap_bond_<type>.dat`` and
  ``eigenvalue_bond_<type>.dat`` per requested channel AFTER every FLEX
  output file, so a failure of the pairing step never costs the FLEX
  results.
- **Post-processing**, through ``hwave_sc``: set ``[eliashberg]
  bond_channels = true`` with ``frequency = "dynamic"`` and
  ``chi0q_mode = "flex"`` in a separate ``hwave_sc`` input, pointing
  ``[file.input] path_to_flex_output`` at the FLEX run's output directory.
  Requires ``longitudinal_bond_output_full = true`` on that FLEX run (the
  archive is written only then) under a version writing archive schema 2
  -- an older schema-1 archive is refused with a message asking for a
  re-run. Solves ONE channel per run (``[eliashberg] pairing_type``), and
  writes the same file names as the on-site dynamic solver
  (``gap_dynamic.npz``, ``gap.dat``, ``eigenvalue.dat``).

See :ref:`the output reference <subsec:eliashberg_bond_outputs>` for the
complete key sets, :ref:`the configuration reference's eliashberg section
<eliashberg_bond_dynamic_config>` for the bond-specific keys and
:ref:`the dynamic-frequency section <sc_dynamic_frequency>` for the general
solver controls.

Memory and the IR basis
""""""""""""""""""""""""""""""""

The archive's two dominant members, ``chi_s_w`` and ``chi_c_w``, are each
``Nmat * nvol * ND**2 * 16`` bytes (complex128; ``ND = B * norb**2`` with
``B`` the number of bond channels the FLEX run kept, channel 0 the on-site
term), so the archive as a whole is about ``2 * Nmat * nvol * ND**2 * 16``
bytes -- the same doubling ``longitudinal_bond_output_full = true`` (above)
already costs the FLEX solve's own persistent memory. For in-process-only
use leave ``longitudinal_bond_output_full`` at ``false``: the archive is
written only for the post-processing entry, and on a multi-orbital model it
is tens of GB on disk. The pairing kernel adds its own working set on top
(the pair bubble, the hoisted vertex blocks, the eigensolver vectors).
Worked numbers for two representative models (CuO2-type: :math:`B=5`, ``norb`` = 3, ``ND`` = 45, ``nvol`` = 1024,
``Nmat`` = 1024, ``num_eigenvalues`` = 10; single band: :math:`B=5`,
``norb`` = 1, a 32x32 lattice, ``Nmat`` = 1024):

.. list-table::
   :header-rows: 1
   :widths: 46 27 27

   * - Quantity
     - CuO2-type (GB)
     - Single band (GB)
   * - Dynamic archive on disk (``chi_s_w`` + ``chi_c_w``)
     - 68
     - < 1
   * - In-process, uniform grid, hoisted vertex blocks (``B**2`` blocks)
     - 34
     - < 1
   * - In-process, uniform grid, ``"stream"`` residency, device memory
       (one A100; one matvec = one FLEX-sized transport, 15-30 s)
     - ~8.3
     - < 1
   * - In-process, IR (:math:`L_B \approx 80`), device memory
       (dominated by the hoisted blocks; everything fits, matvec in
       seconds)
     - ~2.7
     - < 1
   * - Post-processing, uniform grid, host memory before the kernel
       (one archive member + the vertex accumulator)
     - ~68
     - < 1
   * - Post-processing, IR, host memory before ``finish``
       (one archive member + the IR coefficients)
     - ~39
     - < 1

GB here means :math:`10^9` bytes; the cap key ``bond_memory_cap_gb`` is in
binary GiB (:math:`1024^3` bytes, about 7 % larger), so set it from these
numbers with that margin in mind.

For a model whose ``B**2`` hoisted blocks (or, in post-processing, whose
archive members) do not fit in memory, set ``[eliashberg] matsubara_basis
= "ir"``: the IR basis replaces the ``Nmat``-point Matsubara axis by its
sparse sampling nodes (typically 50-100), shrinking every frequency-axis
row above by roughly ``Nmat / L``. It requires the optional
`sparse-ir <https://sparse-ir.readthedocs.io>`_ package
(``pip install sparse-ir``). ``[eliashberg] bond_memory_cap_gb`` caps the
host side of this path explicitly (the existing static-bond key of the
same name, reused here); when even the streaming residency would not fit,
the in-process entry refuses BEFORE the first FLEX map (a long FLEX run is
never lost to a post-hoc allocation failure) and the post-processing entry
refuses with the full memory table. The kernel picks its residency
(``"device"``, ``"host"`` or ``"stream"``: vertex blocks on the device, on
the host, or rebuilt per application) automatically from the free device /
host memory and ``bond_memory_cap_gb``; the choice is recorded in the
outputs as ``bond_residency``.

A single-band example (both entries)
""""""""""""""""""""""""""""""""""""""

On-site :math:`U = 4t` and nearest-neighbour :math:`V = t` on a 32x32
lattice, ``Nmat`` 1024, solving both pairing channels in-process:

.. code-block:: toml

   [mode]
     mode = "FLEX"
     calc_scheme = "general"
   [mode.param]
     T = 0.02
     filling = 0.45
     CellShape = [32, 32, 1]
     Nmat = 1024
     IterationMax = 100
     mixing_scheme = "anderson"
     anderson_depth = 8
     Mix = 0.2
     EPS = 8
     flex_hartree_fock = true
     longitudinal_bond_channels = true
     longitudinal_bond_output_full = true   # only for the hwave_sc entry below; omit for in-process-only runs
     longitudinal_bond_pairing = "both"
   [file.input]
     path_to_input = "."
   [file.input.interaction]
     path_to_input = "."
     Geometry = "geom.dat"
     Transfer = "transfer.dat"
     CoulombIntra = "coulombintra.dat"
     CoulombInter = "coulombinter.dat"
   [file.output]
     path_to_output = "output"
     sigma = "sigma"
     green = "green"
     chiq = "chiq"
     longitudinal_bond = "longitudinal_bond.npz"
     energy = "energy.dat"
   [eliashberg]
     solver_mode = "eigenvalue"
     num_eigenvalues = 6

.. code-block:: bash

   $ hwave input.toml

This writes ``eliashberg_bond_singlet.npz`` / ``gap_bond_singlet.dat`` /
``eigenvalue_bond_singlet.dat`` and the ``triplet`` equivalents into
``output/``. To solve the same archive through the post-processing entry
instead (e.g. to compare bases without re-running FLEX), run the FLEX
input above once (``longitudinal_bond_pairing`` may stay ``"none"`` if
only the post-processing entry is wanted), then a second input through
``hwave_sc``:

.. code-block:: toml

   [mode]
     mode = "SC"
   [mode.param]
     T = 0.02
     CellShape = [32, 32, 1]
     Nmat = 1024
     filling = 0.45
   [file.input]
     path_to_flex_output = "output"
   [file.input.interaction]
     path_to_input = "."
     Geometry = "geom.dat"
     Transfer = "transfer.dat"
     CoulombIntra = "coulombintra.dat"
     CoulombInter = "coulombinter.dat"
   [file.output]
     path_to_output = "output_sc"
   [eliashberg]
     chi0q_mode = "flex"
     frequency = "dynamic"
     bond_channels = true
     pairing_type = "singlet"
     solver_mode = "eigenvalue"
     num_eigenvalues = 6

.. code-block:: bash

   $ hwave_sc input_sc.toml

Both give the leading eigenvalue and gap to round-off on the uniform grid.

.. note::

   ``hwave_sc`` writes ``gap_dynamic.npz`` under that fixed name whatever
   ``pairing_type`` is (``[eliashberg] output_gap`` / ``output_eigenvalue``
   rename only the text files ``gap.dat`` / ``eigenvalue.dat``); solve the
   singlet and the triplet channel into two different ``path_to_output``
   directories, or the second run overwrites the first archive.

A three-band example on the IR basis
""""""""""""""""""""""""""""""""""""""

For a model whose bond-resolved dynamic archive does not fit the uniform
grid's memory (e.g. a CuO2-type three-band model), set ``matsubara_basis
= "ir"`` in the ``[eliashberg]`` table of the FLEX input.

The post-processing entry accepts ``matsubara_basis = "ir"`` as well; the
in-process entry is shown because it never writes the uniform-grid archive.

.. code-block:: toml

   [mode.param]
     ...
     flex_hartree_fock = true
     longitudinal_bond_channels = true
     longitudinal_bond_pairing = "singlet"
   [eliashberg]
     matsubara_basis = "ir"
     ir_fit_tol = 0.1
     solver_mode = "eigenvalue"
     num_eigenvalues = 4

``ir_fit_tol`` (default 0.5) is the componentwise relative residual of
fitting the uniform-grid bond vertex onto the IR basis; real uniform-FFT
archives with the constant retained (``ir_keep_static_chi = true``)
typically give a 0.1-0.2 residual, so the tighter value above is a
reasonable check on a well-resolved run and may need loosening on a
coarser ``Nmat``.

.. note::

   The instantaneous (frequency-independent) part of the vertex enters the
   IR kernel through the exact Matsubara-sum midpoint
   :math:`\tfrac12(F(0^+) - F(\beta^-))`, not :math:`F(0^+)`. On a
   frequency-even pair amplitude the two coincide (its even-:math:`l` IR
   coefficients vanish); the old prescription additionally mapped
   odd-frequency components into the even sector and never the other way,
   so it was block-triangular with respect to frequency parity and had the
   same spectrum -- eigenvalues of parity-projected solves and of the
   eigenvalue solver are unchanged, while eigenvectors and an unprojected
   iteration were affected. What the correction removes is that spurious
   parity leakage (and with it the disabled projection, which in a triplet
   run can converge to the singlet-sector value, and the even-parity
   admixture in odd-channel gap functions). The same correction applies to
   the ordinary on-site dynamic IR solver; see
   :ref:`the change note <sc_dynamic_ir_instantaneous_en>` for who should
   re-run.

Nmat dependence of the IR result
""""""""""""""""""""""""""""""""""""

On the uniform grid the two entries agree to round-off. Against the IR
basis, however, a nonzero instantaneous vertex makes the two differ by
:math:`O(\beta / N_{\rm mat})`: the uniform grid truncates the Matsubara
sum of the frequency-flat term at a finite ``Nmat``, while the IR
representation evaluates it analytically. Measured on the single-band
:math:`U=4`, :math:`V=1`, 4x4 fixture at :math:`T=0.5`:

.. list-table::
   :header-rows: 1
   :widths: 34 33 33

   * - ``Nmat``
     - singlet, uniform vs. IR
     - triplet, uniform vs. IR
   * - 64
     - 5.7 %
     - 3.0 %
   * - 128
     - 1.7 %
     - 0.9 %

The two converge to the same continuum limit as ``Nmat`` grows; for
production use, check this shift at your own model's ``Nmat`` before
trusting the IR result at face value, the same way the on-site dynamic
solver's IR section recommends. The IR representation of a uniform-FFT
archive also carries its own parity asymmetry (unrelated to the above),
decaying as ``Nmat^-2`` (measured 7.6e-3 at Nmat 64, 1.3e-3 at Nmat 128 on
the same fixture) -- the reason ``[eliashberg] parity_leakage_tol``
defaults to a looser ``2e-2`` on the IR basis than the uniform grid's
``1e-8``.

Pairing on a non-converged FLEX solve
""""""""""""""""""""""""""""""""""""""""

The in-process pairing step runs even when the FLEX SCF loop did not
converge (the cost of the run is already sunk): the kernel uses the LAST
map's dressed susceptibility together with the FINAL (post-mix) Green
function, a WARNING is logged, and the outputs record
``scf_converged = false`` and ``state = "mixed: last-map chi, final
green"`` (a converged run instead records
``state = "last_map_chi / final_green"``) -- both in the npz metadata and
as ``#`` header lines of ``eigenvalue_bond_<type>.dat``. Treat such a
result as a diagnostic of the trajectory the SCF loop was following, not
as a self-consistent pairing eigenvalue. A pipeline that reads
``eigenvalue_bond_<type>.dat`` should test the ``# scf_converged=true``
header line (or the ``scf_converged`` key of the npz) before using the
number.


Sample 2: Two-orbital Hubbard model
-----------------------------------------

Model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The second sample is a **two-orbital Hubbard model** with
inter-orbital Coulomb interaction and Hund's coupling:

.. math::

   H = \sum_{\mathbf{k},\alpha,\beta,\sigma}
       \varepsilon_{\alpha\beta}(\mathbf{k})\,
       c^\dagger_{\mathbf{k}\alpha\sigma} c_{\mathbf{k}\beta\sigma}
     + U \sum_{i,\alpha} n_{i\alpha\uparrow} n_{i\alpha\downarrow}
     + V \sum_{i,\alpha\neq\beta} n_{i\alpha} n_{i\beta}
     - 2J \sum_{i,\alpha\neq\beta}
       \mathbf{S}_{i\alpha} \cdot \mathbf{S}_{i\beta}

with :math:`U = 4.0`, :math:`V = 1.0`, :math:`J = 0.5`,
and temperature :math:`T = 1.0`.

Prepare input files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The sample files are in ``docs/en/source/flex/sample/2orb/``.

**Parameter file** (``input.toml``):

.. literalinclude:: ../sample/2orb/input.toml

Key differences from the 1-orbital sample:

- ``T = 1.0``: Higher temperature for stability.
- ``IterationMax = 200``: More iterations for multi-orbital convergence.
- Additional interaction files: ``CoulombInter`` and ``Hund``.

**Geometry** (``geom.dat``):

.. literalinclude:: ../sample/2orb/geom.dat

Two orbitals per unit cell.

**Transfer integrals** (``transfer.dat``):

.. literalinclude:: ../sample/2orb/transfer.dat

Intra-orbital hopping (:math:`t = 1.0` along y) and
inter-orbital hybridization (:math:`t' = 0.5`).

**On-site interaction** (``coulombintra.dat``):

.. literalinclude:: ../sample/2orb/coulombintra.dat

Intra-orbital Coulomb :math:`U = 4.0` on both orbitals.

**Inter-orbital Coulomb** (``coulombinter.dat``):

.. literalinclude:: ../sample/2orb/coulombinter.dat

Inter-orbital Coulomb :math:`V = 1.0`.

**Hund's coupling** (``hund.dat``):

.. literalinclude:: ../sample/2orb/hund.dat

Hund's coupling :math:`J = 0.5`.


Run the calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    $ cd docs/en/source/flex/sample/2orb
    $ hwave input.toml

.. code-block:: text

    FLEX iteration 1/200
    FLEX._find_mu_dressed: mu = 0.000000
      convergence: |dSigma|/|Sigma| = 1.000e+00
    FLEX iteration 2/200
    FLEX._find_mu_dressed: mu = 0.000000
      convergence: |dSigma|/|Sigma| = 3.587e-01
    ...
    FLEX iteration 59/200
    FLEX._find_mu_dressed: mu = 0.000000
      convergence: |dSigma|/|Sigma| = 8.870e-07
    FLEX converged after 59 iterations

(This is a particle-hole symmetric half-filled model, so :math:`\mu = 0`.)


Results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Spin susceptibility** :math:`\chi_s(\mathbf{q}, i\nu_0)`:

.. figure:: ../sample/2orb/chi_s.png
   :width: 60%
   :align: center

   Static spin susceptibility of the two-orbital model.
   The peak at :math:`\mathbf{Q} = (\pi, \pi)` is enhanced by
   Hund's coupling, which promotes ferromagnetic alignment within
   each site while allowing antiferromagnetic inter-site correlations.

**Self-energy** :math:`\mathrm{Im}\,\Sigma(\mathbf{k}, i\omega_0)`:

.. figure:: ../sample/2orb/sigma_kspace.png
   :width: 60%
   :align: center

   Imaginary part of the self-energy for the two-orbital model.
   The orbital-dependent k-structure reflects the different
   scattering channels in the multi-orbital system.

**Self-energy vs Matsubara frequency**:

.. figure:: ../sample/2orb/sigma_matsubara.png
   :width: 80%
   :align: center

   Frequency dependence of the self-energy for the two-orbital model.
   The larger magnitude compared to the single-orbital case reflects
   the enhanced correlations from inter-orbital interactions.


Plotting
----------------------------

The figures above can be reproduced using the plotting script:

.. code-block:: bash

    $ cd docs/en/source/flex/sample
    $ python plot_results.py

Or from within a single sample directory:

.. code-block:: bash

    $ cd docs/en/source/flex/sample/1orb
    $ python ../plot_results.py


Output file format
----------------------------

The FLEX solver produces NumPy ``.npz`` files with the following contents:

.. note::

   Under ``calc_scheme = "general"`` every archive a FLEX run writes --
   ``chi0q.npz``, ``chiq.npz``, ``chiq_s.npz``, ``chiq_c.npz``,
   ``sigma.npz``, ``green.npz`` and the dedicated bond archive -- records
   the two provenance members ``flex_second_order`` and
   ``flex_second_order_schema`` described below, in addition to the
   contents listed for that file. See :ref:`rpa_chiq_provenance`.

``chi0q.npz``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``chi0q``: Bare susceptibility :math:`\chi_0(\mathbf{q}, i\nu_m)`,
  shape ``(nmat, nvol, nd, nd)``.
- ``freq_index``: Matsubara frequency indices.
- ``wavevector_unit``: k-point vectors.
- ``wavevector_index``: Wavenum table.

``chiq_s.npz``, ``chiq_c.npz``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``chiq_s`` / ``chiq_c``: Spin / charge susceptibility,
  same shape as ``chi0q``.
- ``chi_convention``: which spin/charge vertex the susceptibilities are meant
  to be paired with, and which shape family they have: ``"kuroki"`` for the
  reduced scheme (spin-orbital shape, ``nd = norb * ns``) or
  ``"myo"`` for the general full-vertex scheme (orbital-pair shape,
  ``nd = norb^2``; the historical MYO-vs-Kuroki difference in the
  ``C(ab,ab)`` charge vertex was resolved by the exact-diagonalization
  adjudication of the per-type vertex content, so the two builders now
  coincide). The
  Eliashberg loader (``hwave_sc``) uses this tag to interpret the orbital
  indices; it is essential for two-orbital systems, where the spin-orbital and
  orbital-pair dimensions coincide (both ``4``) and shape alone is ambiguous.
- ``chi_orbital_layout``: written by the **general** scheme only, value
  ``"acbd"`` — the four orbital legs are stored as the pairs ``(a,c)`` (row)
  and ``(b,d)`` (column). Reduced-scheme files do not carry it: their axes
  are spin-orbital (``s*norb + a``), not four orbital legs, so the loader has
  to extract a spin block before the array is an orbital-pair object at all.
  The marker exists so that a file written by a pre-fix build of the
  general path — which stored the arrays orbital-pair transposed under the same
  ``"myo"`` tag, and is indistinguishable by tag alone — is rejected on load
  with a regenerate message instead of silently producing a transposed pairing
  vertex, and so that any future layout change fails fast rather than being
  misread.

``sigma.npz``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``sigma``: Self-energy :math:`\Sigma(\mathbf{k}, i\omega_n)`,
  shape ``(nblock, nmat, nvol, nd_block, nd_block)``
  where ``nblock`` is the number of spin blocks (1 for spin-free mode).
- ``flex_second_order`` / ``flex_second_order_schema``: written by
  ``calc_scheme = "general"`` runs only -- the second-order kernel of the
  effective interaction (``local`` | ``takimoto``, a 0-d ``<U8`` string
  array) and the schema version of that record (``1``). Both are absent
  from reduced-scheme and RPA archives and from archives written before
  this key was introduced (H-wave 2.0.0 and earlier); readers never require
  them, and a ``sigma_init`` seed
  recording a different kernel is accepted with a warning. See
  :ref:`rpa_chiq_provenance`.

.. note::

   Multi-orbital ``sigma.npz`` and ``green.npz`` files written by
   ``calc_scheme = "general"`` *before* the orbital-pair transpose fix are wrong
   off the orbital diagonal and must be regenerated -- including any that are
   still being consumed as ``sigma_init`` seeds or fed to ``hwave_sc`` via
   ``bond_green``. The ``chiq_s``/``chiq_c`` of such a run are unaffected. See
   :ref:`the migration warning <flex_general_transpose_fix>`.

``green.npz``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``green``: Dressed Green's function :math:`G(\mathbf{k}, i\omega_n)`,
  same shape as ``sigma``.
- ``flex_second_order`` / ``flex_second_order_schema``: as for
  ``sigma.npz`` above.

These output files can also be used as input for the
Eliashberg equation solver (``hwave_sc``) to analyze
superconducting instabilities. See :doc:`/rpa/tutorial/sc-index` for details.


.. _flex_params:

FLEX-specific parameters
----------------------------

The FLEX solver accepts the following parameters in the
``[mode.param]`` section:

.. list-table::
   :header-rows: 1
   :widths: 20 10 10 60

   * - Parameter
     - Type
     - Default
     - Description
   * - ``IterationMax``
     - int
     - 100
     - Maximum number of SCF iterations.
   * - ``Mix``
     - float
     - 0.2
     - Mixing parameter :math:`\alpha` for self-energy update.
       Smaller values give more stable convergence but slower progress.
   * - ``EPS``
     - int/float
     - 6
     - Convergence criterion. If integer :math:`n`, the threshold is
       :math:`10^{-n}`. If float < 1, used directly as threshold.
   * - ``mixing_scheme``
     - str
     - "linear"
     - Self-energy update scheme. ``"linear"`` is the conventional linear
       mixing :math:`\Sigma \leftarrow (1-\alpha)\Sigma +
       \alpha\Sigma_{\mathrm{new}}`. ``"anderson"`` enables Anderson
       acceleration (Pulay/DIIS-type extrapolation over a short
       iterate/residual history), which reaches the same fixed point in far
       fewer iterations (e.g. 78 -> 13 on the 8x8 Hubbard benchmark at
       :math:`U=3.5`, ``Mix=0.2``). Falls back to a plain linear step
       automatically if the history becomes degenerate.
   * - ``anderson_depth``
     - int
     - 5
     - History depth :math:`m` of the Anderson acceleration. Memory grows by
       :math:`2m` sigma-sized arrays (kept on the device under GPU
       execution).
   * - ``matsubara_basis``
     - str
     - "uniform"
     - Matsubara-axis representation: ``"uniform"`` (default, unchanged) or
       ``"ir"`` (the sparse-ir intermediate representation; chi0 and Sigma
       are computed NATIVELY on sparse nodes, so the uniform-FFT
       :math:`O(\beta/N_{\mathrm{mat}})` discretization artifacts do not
       arise by construction). ``Nmat`` keeps its role as the output grid
       (all output files are densified onto it). Supported with
       ``calc_scheme = "reduced"`` and ``"general"``.
       Requires the optional
       `sparse-ir <https://sparse-ir.readthedocs.io>`_ package. The mu
       search becomes the basis evaluation of
       :math:`n = -\mathrm{Tr}\,G(\tau=\beta^-)`, and ``coeff_tail`` is
       unnecessary (ignored).
   * - ``ir_tol``
     - float
     - 1e-8
     - IR basis cutoff accuracy.
   * - ``ir_wmax``
     - float
     - auto
     - Real-frequency bandwidth of the IR basis (same energy units as the
       Hamiltonian); auto-estimated from the band range and interaction
       scale when omitted (a fail-fast error asks for an explicit value if
       the estimate cannot be formed). An always-on coefficient-decay
       diagnostic warns when the bandwidth is insufficient.
   * - ``sigma_init_on_error``
     - str
     - "warn"
     - Behavior when the IR fit residual of a (uniform-grid) ``sigma_init``
       exceeds 100x ``ir_tol``: ``"warn"`` (use it, warn), ``"abort"``, or
       ``"zero"`` (fall back to the zero start).
   * - ``write_densified``
     - bool
     - true
     - IR runs only. ``true`` (default): all output files are densified
       onto the uniform ``Nmat`` grid (unchanged format, readable by every
       tool). ``false``: outputs stay on the sparse IR nodes — the fixed
       densify+write cost disappears (the dominant remaining cost of an IR
       run), files shrink by ~``Nmat``/L, and downstream IR-aware
       consumers (the dynamic Eliashberg solver with
       ``[eliashberg] matsubara_basis = "ir"``, and ``sigma_init``
       chaining into another IR FLEX run) read them directly. Uniform-only
       readers (static ``hwave_sc``, ``chi0q_init``, legacy scripts)
       reject such files with an explicit error. See the note below.
   * - ``flex_second_order``
     - str
     - "local"
     - ``calc_scheme = "general"`` only. Second-order kernel of the
       effective interaction: ``"local"`` (default) is the exact local
       second order of every accepted interaction term;
       ``"takimoto"`` keeps the legacy Takimoto-Hotta-Ueda expression and
       reproduces the results of H-wave 2.0.0 and earlier releases. See
       :ref:`flex_second_order_tutorial` and the configuration reference.
   * - ``gpu``
     - bool
     - false
     - Set ``true`` to run the SCF loop (dressed G, chi0q, chiq, V_eff, and
       the self-energy) on a GPU via CuPy. When CuPy or a CUDA device is
       unavailable the solver warns and falls back to the CPU (numpy) path
       (identical result). The chemical-potential search also runs on the GPU
       via closed-form eigenvalues when each spin block has at most 2
       components (single-orbital, or e.g. a spin-reduced two-orbital model);
       only larger blocks fall back to a host non-Hermitian
       eigendecomposition. ``flex_hartree_fock = true`` on its own uses
       this ordinary GPU path; with ``longitudinal_bond_channels = true``
       the bond-resolved GPU path described in
       :ref:`flex_bond_hf_tutorial` applies instead.
   * - ``fft_workers``
     - int
     - 1
     - Number of worker threads for the spatial FFTs (parallelized via
       ``scipy.fft``). The default ``1`` keeps the serial numpy path,
       unchanged from previous releases (opt-in); ``-1`` uses all cores.
       Ignored on the GPU. Set a smaller number when running several
       calculations concurrently.

All other parameters (``T``, ``CellShape``, ``Nmat``, ``filling``, etc.)
are shared with the RPA solver. See :ref:`Ch:Config_rpa` for details.

.. note::

   **Running FLEX on the IR basis.** Install the optional dependency once
   (``pip install sparse-ir``), then add a single line under
   ``[mode.param]`` of any existing FLEX input — every other line, including
   ``Nmat``, stays as it is:

   .. code-block:: toml

      [mode]
      mode = "FLEX"
      calc_scheme = "reduced"     # or "general"
      [mode.param]
      CellShape = [64, 64, 1]
      T = 0.05
      Nmat = 4096                 # still required: the output grid
      matsubara_basis = "ir"      # opt in to the sparse-IR axis
      # ir_tol = 1e-8             # optional: basis cutoff accuracy
      # ir_wmax = 30.0            # optional: bandwidth (auto-estimated)

   The SCF then runs on a few dozen sparse nodes instead of ``Nmat``
   frequencies (e.g. 4096 -> 42 at :math:`T=0.05`), while all output files
   are densified back onto the ``Nmat`` grid, so downstream tools (including
   the dynamic Eliashberg solver) work unchanged. ``coeff_tail`` is ignored
   on this path — the IR basis carries the :math:`1/(i\omega)` tail exactly.

   **IR-native outputs.** Adding ``write_densified = false`` keeps the
   outputs on the sparse nodes. Use it for pure IR chains — IR FLEX
   feeding the dynamic Eliashberg solver
   (``[eliashberg] matsubara_basis = "ir"``) or seeding the next IR FLEX
   run of a temperature sweep via ``sigma_init`` (cross-temperature seeds
   are supported) — where it removes the fixed densification and file-size
   cost entirely. You can recognize such a file by the
   ``frequency_grid = "sparse_ir_nodes"`` key in the ``.npz``.

   .. warning::

      Do **not** point legacy analysis scripts (anything that indexes the
      frequency axis positionally, e.g. the static slice at ``Nmat/2``) at
      ``write_densified = false`` outputs — the frequency axis holds
      sparse nodes, not the uniform grid. All H-wave readers detect this
      and stop with an explicit error; external scripts may not. To
      recover a uniform-grid file, re-run FLEX with
      ``write_densified = true`` (cheap: seed it with ``sigma_init`` from
      the native run), or densify offline::

         import numpy as np
         from hwave.solver.ir_axis import IRAxis
         d = np.load("chiq_s.npz")
         ax = IRAxis(float(d["ir_beta"]), float(d["ir_wmax"]),
                     float(d["ir_tol"]), str(d["ir_statistics"]))
         c = ax.fit_from_freq_points(np.moveaxis(d["chiq_s"], 0, -1),
                                     d["ir_freq_n"])
         chi_u = np.moveaxis(ax.eval_to_uniform(c, nmat=4096), -1, 0)

.. note::

   The FLEX solver accepts ``calc_scheme`` in ``"reduced"`` or ``"general"``.
   The ``"reduced"`` scheme consumes the reduced-shape susceptibility and
   solves with the density-density part of the interaction; it **rejects**
   ``Exchange`` and ``PairHop`` (whose
   vertex has no density-density content — the input would silently have
   zero effect). The ``"general"`` scheme is the paramagnetic
   full-vertex path: it keeps the full Kanamori vertices (MYO formula, see
   :ref:`above <flex_scope>`), but it is **spin-free only** — it raises a
   ``ValueError`` for ``spin_mode = "spin-diag"`` or ``"spinful"`` and
   rejects ``enable_spin_orbital``. Off-site input is accepted for
   ``CoulombInter``, ``Hund`` and ``Ising`` (same-orbital or inter-orbital,
   with or without sublattice folding): each enters the RING vertex as its
   Hartree (density) part :math:`V(q)` only — the exchange crossing of an
   off-site term is not representable by a :math:`q`-only vertex and is
   left out, the same approximation the RPA ring makes, and the solver
   logs a warning saying so. (At SECOND order the default kernel
   ``flex_second_order = "local"`` is exact for the direct skeleton and for
   every mixed on-site/off-site diagram of those terms; only the exchange
   skeleton of two off-site vertices is left out there, see
   :ref:`flex_second_order_tutorial`.) For every such class the general path is
   measured element-complete equal to the RPA ring. (The omitted
   crossing is available, statically, in the RPA solver's experimental
   bond-resolved longitudinal channel, ``longitudinal_bond_channels =
   true``, see :ref:`rpa_longitudinal_bond`, and self-consistently in the
   Hartree-Fock FLEX, see :ref:`flex_bond_hf_tutorial`.) Off-site ``Exchange``
   and ``PairHop`` raise a ``ValueError``. An off-site ``Exchange`` has no
   effect a :math:`q`-dependent spin/charge vertex could carry (verified
   by exact diagonalization): its physics is spin-flip (transverse), which
   this spin-free path does not compute, and the small remainder needs a
   bond-resolved vertex. An off-site ``PairHop`` has no local-pair form.
   Remove such entries from the interaction files for a FLEX run; an
   interaction set prepared for the RPA solver (which reads them) is not
   reusable as-is. On-site ``Exchange`` and
   ``PairHop`` off-diagonal vertices **are kept** (the point of the scheme),
   but ``PairLift`` contributes ``S=C=0`` to the particle-hole vertex and is
   **inert** (accepted with a note that it is exactly zero). The general path writes ``chiq_s``/
   ``chiq_c`` in the MYO convention (tagged ``chi_convention="myo"``), which
   ``hwave_sc`` reads back automatically. In all schemes
   ``calc_type = "ring+ladder"`` is **not** supported (the solver raises a
   ``ValueError``).

   .. _flex_general_transpose_fix:

   .. warning::

      **Result change for** ``calc_scheme = "general"``. Until this fix the
      general path applied a spurious orbital-pair transpose when building its
      effective interaction, so the self-energy was built from the transposed
      bubble: ``sigma.npz`` — and everything derived from it (occupations,
      energies, Eliashberg eigenvalues) — was correct on the orbital diagonal
      and wrong **off** it. The scheme has only ever existed in development
      builds, so this affects work done on ``develop``, not any released
      version. Single-orbital runs are unaffected (the transpose is the
      identity there); multi-orbital runs are affected whenever the Green
      function has orbital off-diagonal weight. For a
      density-only interaction (``CoulombIntra`` alone) the general and reduced
      schemes now agree on the self-energy to machine precision, as they must.

      The saved susceptibilities do **not** change. ``chi0q.npz`` was already
      transposed back at the output boundary, and ``chiq_s``/``chiq_c`` were
      already corrected separately (the ``chi_orbital_layout`` marker and the
      ``[a,c,b,d]`` write order). Removing the transpose at its source makes
      those output-boundary corrections unnecessary rather than changing their
      result: measured on a 2-orbital Kanamori model the saved ``chi0q`` is
      bit-identical and ``chiq_s``/``chiq_c`` agree to ~10⁻¹⁴ relative, the
      residual being floating-point error of the matrix products and the
      channel linear solve. That equivalence is exact algebra only while the
      spin/charge interaction matrices are symmetric under the orbital-pair
      transpose, which holds whenever the on-site interaction parameters are
      symmetric in their orbital indices — the physical case. H-wave does not
      currently check interaction files for that symmetry.

      **What to regenerate.** ``sigma.npz`` and ``green.npz`` from any
      multi-orbital ``calc_scheme = "general"`` run made before this fix, plus
      anything derived from them. ``chi0q``/``chiq_s``/``chiq_c`` do not need
      regenerating on this account. Single-orbital runs are unaffected, and
      ``calc_scheme = "reduced"`` was never affected.

   **Memory with** ``"general"`` **+ IR.** ``matsubara_basis = "ir"`` (see
   above) works with ``calc_scheme = "general"``, but IR only compresses the
   *frequency* axis, typically by ``Nmat/L`` (20-40x); it does not change the
   scaling in :math:`N_{\mathrm{orb}}`. Because ``"general"`` keeps the full
   rank-4 orbital vertex, chi0q/chiq/sigma storage already scales as
   :math:`O(N_{\mathrm{orb}}^4)` (and orbital contractions scale worse), so
   ``"general"`` + IR extends the reachable ``Nmat``/:math:`\beta` but not the
   reachable :math:`N_{\mathrm{orb}}`. Watch memory on multi-orbital
   ``"general"`` runs (the FLEX GPU path's VRAM preflight already warns when a
   run looks under-provisioned).


Sample 3: Iron pnictide 2-orbital model
-----------------------------------------

Model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The third sample is a **two-orbital minimal model for iron-based
superconductors** proposed by Raghu et al. [5]_
The model describes the Fe-As plane using :math:`d_{xz}` and
:math:`d_{yz}` orbitals on a square lattice (1-Fe unit cell).

.. math::

   H_0(\mathbf{k}) = \begin{pmatrix}
   \varepsilon_x(\mathbf{k}) & \varepsilon_{xy}(\mathbf{k}) \\
   \varepsilon_{xy}(\mathbf{k}) & \varepsilon_y(\mathbf{k})
   \end{pmatrix}

where

.. math::

   \varepsilon_x(\mathbf{k}) &= -2t_1 \cos k_x - 2t_2 \cos k_y - 4t_3 \cos k_x \cos k_y \\
   \varepsilon_y(\mathbf{k}) &= -2t_2 \cos k_x - 2t_1 \cos k_y - 4t_3 \cos k_x \cos k_y \\
   \varepsilon_{xy}(\mathbf{k}) &= -4t_4 \sin k_x \sin k_y

with :math:`t_1 = -1.0`, :math:`t_2 = 1.3`, :math:`t_3 = t_4 = -0.85`.

The interactions follow the Kanamori parameterization:

.. math::

   H_{\mathrm{int}} = U \sum_{i,\alpha} n_{i\alpha\uparrow} n_{i\alpha\downarrow}
   + U' \sum_{i,\alpha\neq\beta} n_{i\alpha} n_{i\beta}
   - 2J \sum_{i,\alpha\neq\beta} \mathbf{S}_{i\alpha} \cdot \mathbf{S}_{i\beta}
   + J' \sum_{i,\alpha\neq\beta} c^\dagger_{i\alpha\uparrow} c^\dagger_{i\alpha\downarrow}
     c_{i\beta\downarrow} c_{i\beta\uparrow}

with :math:`U = 1.5`, :math:`J = J' = 0.25`, :math:`U' = U - 2J = 1.0`,
at temperature :math:`T = 0.1` and half filling (:math:`n = 2`).

The Fermi surface consists of hole pockets at :math:`\Gamma` and
electron pockets at :math:`M = (\pi, 0)` / :math:`(0, \pi)`.
The nesting between these pockets drives strong spin fluctuations
at :math:`\mathbf{Q} = (\pi, 0)`, which is the hallmark of iron pnictides.

.. [5] S. Raghu, X.-L. Qi, C.-X. Liu, D. J. Scalapino, and S.-C. Zhang,
   Phys. Rev. B **77**, 220503(R) (2008).

Prepare input files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The sample files are in ``docs/en/source/flex/sample/iron_2orb/``.

**Parameter file** (``input.toml``):

.. literalinclude:: ../sample/iron_2orb/input.toml

**Geometry** (``geom.dat``):

.. literalinclude:: ../sample/iron_2orb/geom.dat

Two orbitals (:math:`d_{xz}` and :math:`d_{yz}`) at the same site.

**Transfer integrals** (``transfer.dat``):

.. literalinclude:: ../sample/iron_2orb/transfer.dat

The hopping parameters produce the characteristic two-pocket Fermi surface.

**Interactions**:

.. literalinclude:: ../sample/iron_2orb/coulombintra.dat

.. literalinclude:: ../sample/iron_2orb/coulombinter.dat

.. literalinclude:: ../sample/iron_2orb/hund.dat

.. literalinclude:: ../sample/iron_2orb/exchange.dat


Run the calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    $ cd docs/en/source/flex/sample/iron_2orb
    $ hwave input.toml

.. code-block:: text

    FLEX iteration 1/200
    FLEX._find_mu_dressed: mu = 1.562757
      convergence: |dSigma|/|Sigma| = 1.000e+00
    FLEX iteration 2/200
    FLEX._find_mu_dressed: mu = 1.551623
      convergence: |dSigma|/|Sigma| = 7.139e-01
    ...
    FLEX iteration 62/200
    FLEX._find_mu_dressed: mu = 1.512917
      convergence: |dSigma|/|Sigma| = 8.716e-07
    FLEX converged after 62 iterations


Results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Spin and charge susceptibilities**:

.. figure:: ../sample/iron_2orb/chi_spin_charge.png
   :width: 90%
   :align: center

   Static spin susceptibility :math:`\chi_s(\mathbf{q})` (left) and charge
   susceptibility :math:`\chi_c(\mathbf{q})` (right).
   The spin susceptibility peaks at :math:`\mathbf{Q} = (\pi, 0)` and
   :math:`(0, \pi)`, reflecting the nesting between hole and electron
   Fermi pockets. This is qualitatively different from the single-band
   Hubbard model where :math:`\chi_s` peaks at :math:`(\pi, \pi)`.

**Orbital-resolved self-energy**:

.. figure:: ../sample/iron_2orb/sigma_orbital.png
   :width: 90%
   :align: center

   Imaginary part of the self-energy at the lowest Matsubara frequency,
   resolved by orbital. The :math:`d_{xz}` orbital shows stronger
   scattering along :math:`k_y` direction, while :math:`d_{yz}` shows
   stronger scattering along :math:`k_x`. This orbital anisotropy
   arises from the orbital character of the Fermi surface.

**Self-energy vs Matsubara frequency**:

.. figure:: ../sample/iron_2orb/sigma_matsubara_orbital.png
   :width: 90%
   :align: center

   Frequency dependence of the orbital-resolved self-energy at
   high-symmetry k-points. At :math:`M = (\pi, 0)`, the
   :math:`d_{yz}` orbital is more strongly damped than :math:`d_{xz}`,
   reflecting orbital-selective correlations.

**Plotting script**:

.. code-block:: bash

    $ python plot_results.py


Sample 3b: Full-vertex (general) variant
-----------------------------------------

The iron pnictide model above is also provided as a full-vertex variant
that selects ``calc_scheme = "general"`` explicitly. Since ``Exchange`` is
present in the interaction files, Sample 3's ``calc_scheme = "auto"`` now
resolves to the **same** general scheme automatically — the two samples run
the identical full-vertex path, and this variant differs only in making the
choice explicit. It uses the **same** model and interaction files
(``CoulombIntra``, ``CoulombInter``, ``Hund``, ``Exchange``) and retains
the full off-diagonal Kanamori vertices (the spin-flip Hund and
pair-hopping / exchange terms). This is the paramagnetic full-vertex MYO
formulation [3]_ (corroborated by THU [4]_).

This variant is appropriate for multi-orbital models in which the
Hund/exchange/pair-hopping off-diagonal vertices matter — such as the iron
pnictide model here — and a paramagnetic FLEX is sufficient. Note that the
``"general"`` scheme is **spin-free only** (it raises a ``ValueError`` for
``spin_mode = "spin-diag"``/``"spinful"`` and does not support
``enable_spin_orbital``) and does not support ``calc_type = "ring+ladder"``.

The sample files are in
``docs/en/source/flex/sample/iron_2orb_general/``.

**Parameter file** (``input.toml``):

.. literalinclude:: ../sample/iron_2orb_general/input.toml

The only difference from Sample 3 is the explicit
``calc_scheme = "general"`` in the ``[mode]`` section (Sample 3's
``"auto"`` resolves to the same scheme); the geometry, transfer, and
interaction files are identical, and so are the results.


.. _flex_tips:

Tips
----------------------------

- **Convergence issues**: If the SCF loop does not converge, try
  reducing ``Mix`` (e.g., 0.1 or 0.05) or increasing the temperature.
  Strong correlations near magnetic instabilities can make convergence
  difficult.

- **Matsubara frequencies**: A sufficient number of Matsubara
  frequencies (``Nmat``) is needed for accurate results. A good rule
  of thumb is :math:`N_{\mathrm{mat}} \geq 10 / T` to capture the
  low-frequency structure.

- **k-point mesh**: The mesh size (``CellShape``) should be large enough
  to resolve the momentum structure of susceptibilities and self-energy.
  For 2D systems, 8x8 is sufficient for qualitative results; 32x32 or
  larger is recommended for quantitative calculations.

- **Computational cost**: FLEX is more expensive than RPA due to the
  SCF loop. The cost scales as
  :math:`O(N_{\mathrm{iter}} \times N_k \times N_\omega \times N_d^3)`
  where :math:`N_d = N_{\mathrm{orb}} \times N_{\mathrm{spin}}`.

- **Connection to Eliashberg equation**: The FLEX output files
  (``chiq_s.npz`` and ``chiq_c.npz``) can be used with ``hwave_sc``
  by setting ``chi0q_mode = "flex"`` in the ``[eliashberg]`` section.
  This enables analysis of superconducting instabilities with
  FLEX-level spin and charge fluctuations.

  .. note::

     On the ``reduced`` route only the density-density
     susceptibility :math:`\chi_{(a,a),(b,b)}` is stored, so the pairing vertex
     is FLEX-dressed in full **only** for ``CoulombIntra``-only models (or
     ``norb = 1``). With ``CoulombInter``, ``Hund`` or ``Ising`` the
     off-density channels enter undressed and the solver warns.
     ``Exchange`` and ``PairHop`` are **rejected** with a reduced
     susceptibility: they have no density-diagonal vertex content at all, so
     nothing of them would be dressed -- and the scheme refuses them at the
     FLEX/RPA stage anyway. Use ``calc_scheme = "general"`` for the complete
     vertex. See
     :ref:`the Eliashberg supported-interactions note <sc_supported_inter>`.

  .. note::

     ``chi0q_mode = "flex"`` is **rejected** together with
     ``[eliashberg] bond_channels = true``: the bond-resolved path builds its
     own bond-resolved :math:`\bar\chi` bubble directly from the Green
     function, so it never reads a ``chi0q``/``chiq`` file. To feed a
     FLEX-dressed Green function into the bond path, run FLEX to convergence
     and point ``[eliashberg] bond_green`` at its ``green.npz`` output; see
     the ``bond_green`` parameter in :doc:`/rpa/tutorial/sc-index`.


Implementation details and limitations
-----------------------------------------

Supported interaction types
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The FLEX solver supports the following interaction types:

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Interaction type
     - Support
     - Notes
   * - ``CoulombIntra``
     - Yes
     - Intra-orbital Coulomb repulsion :math:`U`
   * - ``CoulombInter``
     - Yes
     - Inter-orbital Coulomb repulsion :math:`V`
   * - ``Hund``
     - Yes
     - Hund's coupling :math:`J`
   * - ``Exchange``
     - general only
     - Exchange interaction :math:`J'` (no density-density vertex
       content: rejected under ``reduced``; ``auto``
       selects ``general``)
   * - ``Ising``
     - Yes
     - Ising-type interaction
   * - ``PairLift``
     - Inert
     - Pair lifting interaction (particle-hole vertex exactly zero;
       accepted in every scheme, with no effect on the
       susceptibility channels)
   * - ``PairHop``
     - general only
     - Pair hopping interaction (no density-density vertex content:
       rejected under ``reduced``; ``auto`` selects
       ``general``)
   * - ``InterAll``
     - **No**
     - Arbitrary 4-body interaction (UHFr solver only)

.. note::

   The ``InterAll`` format is not available in k-space solvers (RPA/FLEX).
   Interactions described by ``InterAll`` should be decomposed into
   the individual interaction types listed above.


Momentum-dependent interactions (long-range interactions)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Interactions are specified in Wannier90-format input files.
Each line has the format ``rx ry rz a b Re Im``,
where ``(rx, ry, rz)`` is a real-space lattice vector.

**On-site interactions** (only ``rx = ry = rz = 0``):

After FFT, these become momentum-independent interactions
:math:`W(\mathbf{q}) = W_0`, constant across all q-points.
All current sample files use this case.

**Long-range interactions** (with non-zero ``(rx, ry, rz)``):

These are automatically transformed into momentum-dependent
interactions via FFT:

.. math::

   W(\mathbf{q}) = \sum_{\mathbf{r}} W(\mathbf{r})\, e^{+i\mathbf{q}\cdot\mathbf{r}}

For example, to include nearest-neighbor Coulomb interactions,
add entries with lattice vectors ``(1,0,0)``, ``(0,1,0)``, etc.
to the interaction file.

.. note::

   There is no built-in facility to automatically discretize a
   continuous :math:`1/r` Coulomb potential. Users must explicitly
   specify the interaction value at each lattice point in the input file.


Spin-charge channel decomposition constraint
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The FLEX solver decomposes the interaction tensor into spin and charge
channels for the self-energy calculation. This decomposition reduces
the 4-body interaction tensor
:math:`W_{\alpha\sigma,\beta\sigma',\alpha\sigma,\beta\sigma'}`
to a 2-body contracted form :math:`W_{\alpha\sigma,\beta\sigma'}`.

Specifically, same-spin and cross-spin components are separated:

.. math::

   U_s &= W_{\mathrm{cross}} - W_{\mathrm{same}} \\
   U_c &= W_{\mathrm{cross}} + W_{\mathrm{same}}

where :math:`W_{\mathrm{same}}` is the same-spin interaction and
:math:`W_{\mathrm{cross}}` is the cross-spin interaction.

This contraction is exact for **density-density type interactions**.
``CoulombIntra``, ``CoulombInter``, ``Hund``, and ``Ising``
are all density-density type and are handled correctly.
``Exchange`` and ``PairHop`` have **no** density-density vertex content,
so this reduction cannot represent them at all: the solver rejects them
under ``reduced`` with a ``ValueError`` pointing to
``calc_scheme = "general"`` (``auto`` selects it for you).
``PairLift``'s particle-hole vertex is exactly zero, so it is accepted
in every scheme and its absence from the channels is exact.


Spin degrees of freedom (spin-free mode)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The FLEX solver operates in **spin-free mode** by default.
This mode assumes SU(2) spin symmetry to reduce computational cost.

In spin-free mode:

- Green's functions are represented in orbital space only
  (shape: ``(1, nmat, nvol, norb, norb)``).
- Susceptibilities and effective interactions are internally inflated to
  spin-orbital space (``nd = norb × ns``) for computation.
- After self-energy computation, the properties guaranteed by SU(2) symmetry
  — :math:`\Sigma_{\uparrow\uparrow} = \Sigma_{\downarrow\downarrow}` and
  :math:`\Sigma_{\uparrow\downarrow} = 0` — are used to contract back
  to orbital space.

.. note::

   Spin-free mode assumes a paramagnetic state (no magnetic order).
   Describing magnetically ordered phases requires an extension that
   treats spin degrees of freedom explicitly.
