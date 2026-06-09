==========================================
Tutorial: Eliashberg equation solver
==========================================

This tutorial demonstrates how to use ``hwave_sc``,
the linearized Eliashberg equation solver included in H-wave.
The tool analyzes superconducting instabilities by computing
the leading eigenvalue of the linearized Eliashberg equation
using the bare susceptibility :math:`\chi_0(\mathbf{q})` from H-wave's RPA solver.

The sample files for this tutorial are located in
``docs/en/source/rpa/sample_sc`` directory.


Overview of the workflow
----------------------------

The calculation proceeds in two steps:

1. **RPA calculation** (``hwave``): Compute the bare susceptibility
   :math:`\chi_0(\mathbf{q})` and save it to ``chi0q.npz``.
2. **Eliashberg solver** (``hwave_sc``): Load ``chi0q.npz``, reconstruct
   the Green's function, compute RPA vertices, and solve the
   linearized Eliashberg equation.


Model
----------------------------

This tutorial uses a **two-orbital tight-binding model** on a
two-dimensional square lattice at 3/4 filling.
The Hamiltonian is:

.. math::

   H = \sum_{\mathbf{k},\alpha,\beta,\sigma}
       \varepsilon_{\alpha\beta}(\mathbf{k})\,
       c^\dagger_{\mathbf{k}\alpha\sigma} c_{\mathbf{k}\beta\sigma}
     + U \sum_{i,\alpha} n_{i\alpha\uparrow} n_{i\alpha\downarrow}
     + \sum_{i,\alpha\neq\beta} V_{\alpha\beta}\,
       n_{i\alpha} n_{i\beta}

with on-site Coulomb repulsion :math:`U = 0.4` and
inter-orbital Coulomb interaction :math:`V`.

This sample is based on the conduction layer model of the organic conductor
:math:`\beta`\ -(meso-DMBEDT-TTF)\ :math:`_2`\ PF\ :math:`_6`,
with transfer integrals obtained from extended Huckel calculations. [1]_


Theory
----------------------------

Superconducting susceptibility
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The RPA charge susceptibility :math:`\hat{X}^c` and spin susceptibility
:math:`\hat{X}^s` are given by:

.. math::

   \hat{X}^c = (\hat{I} + \hat{X}^{(0)} (\hat{U} + 2\hat{V}))^{-1} \hat{X}^{(0)}

.. math::

   \hat{X}^s = (\hat{I} - \hat{X}^{(0)} \hat{U})^{-1} \hat{X}^{(0)}

where :math:`\hat{X}^{(0)}` is the bare susceptibility,
:math:`\hat{U}` the on-site interaction matrix,
and :math:`\hat{V}` the inter-site interaction matrix.

Linearized Eliashberg equation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The linearized Eliashberg equation for singlet superconductivity reads:

.. math::

   \lambda_S \Sigma^a_{\alpha\sigma;\beta\bar{\sigma}}(\mathbf{k})
   = -\frac{T}{N_L} \sum_{\mathbf{k}',n',\alpha',\beta'}
   P^S_{\alpha\sigma;\beta\bar{\sigma}}(\mathbf{k} - \mathbf{k}')
   G^{(0)}_{\alpha\alpha'}(\mathbf{k}', i\varepsilon_{n'})
   G^{(0)}_{\beta\beta'}(-\mathbf{k}', -i\varepsilon_{n'})
   \Sigma^a_{\alpha'\sigma;\beta'\bar{\sigma}}(\mathbf{k}')

The pairing interaction for the singlet channel is:

.. math::

   \hat{P}^S = \hat{U} + \hat{V}
   + \frac{3}{2} \hat{U} \hat{X}^s \hat{U}
   - \frac{1}{2} (\hat{U} + 2\hat{V}) \hat{X}^c (\hat{U} + 2\hat{V})

For the triplet channel:

.. math::

   \hat{P}^T = \hat{V}
   - \frac{1}{2} \hat{U} \hat{X}^s \hat{U}
   - \frac{1}{2} (\hat{U} + 2\hat{V}) \hat{X}^c (\hat{U} + 2\hat{V})

The superconducting transition corresponds to :math:`\lambda_S = 1`
(:math:`\lambda_T = 1` for triplet). When :math:`\lambda > 1`,
the normal state is unstable toward superconductivity.
Note that only positive eigenvalues indicate SC instability;
negative eigenvalues correspond to sign-changing gap functions
but do not satisfy the self-consistency condition :math:`\Delta = K\Delta`.

``hwave_sc`` implements two numerical methods to solve the Eliashberg equation:
self-consistent power iteration (converges to the mode with the largest eigenvalue)
and Arnoldi eigenvalue analysis.

.. [1] K. Yoshimi, M. Nakamura, and H. Mori,
   J. Phys. Soc. Jpn. **76**, 024706 (2007);
   `arXiv:cond-mat/0608466 <https://arxiv.org/abs/cond-mat/0608466>`_.


Prepare input files
----------------------------

Parameter file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a TOML parameter file ``input.toml``:

.. literalinclude:: ../sample_sc/input.toml

The file contains the following sections:

``[mode.param]`` section
""""""""""""""""""""""""""""""""

- ``T``: Temperature.
- ``CellShape``: Size of the k-point mesh (32 x 32 x 1 for a 2D system).
- ``Nmat``: Number of Matsubara frequencies (512).
- ``filling``: Electron filling per orbital per spin (0.75 = 3/4 filling).

``[file]`` section
""""""""""""""""""""""""""""""""

- ``[file.input.interaction]``: Specifies files for geometry, transfer integrals,
  and interaction parameters. These files are shared with the RPA step.
- ``[file.output]``: Output directory and filenames for :math:`\chi_0(\mathbf{q})`
  and :math:`\chi(\mathbf{q})`.

``[eliashberg]`` section
""""""""""""""""""""""""""""""""

This section controls the Eliashberg solver. Key parameters:

- ``solver_mode``: ``"iteration"`` (self-consistent power iteration),
  ``"eigenvalue"`` (Arnoldi eigenvalue analysis), or ``"both"``.
- ``chi0q_mode``: ``"load"`` reads :math:`\chi_0(\mathbf{q})` from the RPA output
  file; ``"calc"`` computes it internally.
- ``pairing_type``: ``"singlet"`` or ``"triplet"``.
- ``init_gap``: Initial gap symmetry for iteration.
  Options include ``"cos"`` (:math:`\cos k_x + \cos k_y`),
  ``"dx2y2"`` (:math:`\cos k_x - \cos k_y`), ``"random"``, etc.
- ``max_iter``: Maximum number of self-consistent iterations.
- ``alpha``: Mixing parameter (0 = no mixing, 1 = full mixing of old solution).
- ``convergence_tol``: Convergence criterion on the gap function.
- ``num_eigenvalues``: Number of eigenvalues to compute in eigenvalue mode.
- ``eigenvalue_method``: ``"arnoldi"`` (default), ``"subspace"``, or
  ``"shift-invert-gmres"`` / ``"shift-invert-bicgstab"`` / ``"shift-invert-lgmres"``.

Interaction definition files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The interaction definition files use the Wannier90-like format,
shared with the RPA solver. See :ref:`Ch:Config_rpa` for details.

``Geometry`` (``geom.dat``):

.. literalinclude:: ../sample_sc/geom.dat

Defines a unit cell with 2 orbitals.

``Transfer`` (``transfer.dat``):

.. literalinclude:: ../sample_sc/transfer.dat

Defines the hopping integrals for the two-orbital model.

``CoulombIntra`` (``coulombintra.dat``):

.. literalinclude:: ../sample_sc/coulombintra.dat

On-site Coulomb repulsion :math:`U = 0.4` on each orbital.

``CoulombInter`` (``coulombinter.dat``):

.. literalinclude:: ../sample_sc/coulombinter.dat

Inter-orbital and inter-site Coulomb interactions.


Step 1: Run RPA calculation
----------------------------

First, compute the bare susceptibility by running the RPA solver:

.. code-block:: bash

    $ hwave input.toml

This generates ``output/chi0q.npz`` and ``output/chiq.npz``.
The RPA step takes a few seconds for a 32 x 32 mesh.


Step 2: Run Eliashberg solver
---------------------------------

Next, run the Eliashberg equation solver using the same input file:

.. code-block:: bash

    $ hwave_sc input.toml

The solver performs the following steps:

1. Loads :math:`\chi_0(\mathbf{q})` from ``output/chi0q.npz``.
2. Reads the interaction files and builds the Hamiltonian.
3. Constructs the non-interacting Green's function :math:`G(\mathbf{k}, i\omega_n)`.
4. Computes the RPA charge and spin vertices
   :math:`V_c(\mathbf{q})` and :math:`V_s(\mathbf{q})`.
5. Solves the linearized Eliashberg equation by self-consistent
   iteration and/or eigenvalue analysis.

The output log shows the iteration progress:

.. code-block:: text

    hwave_sc: === Self-consistent iteration ===
    hwave_sc: Iteration    0: eigenvalue = 0.924446, diff = 3.544353e-01
    hwave_sc: Iteration    1: eigenvalue = 0.817270, diff = 7.893848e-02
    ...
    hwave_sc: Iteration  192: eigenvalue = 0.959725, diff = 9.900091e-06
    hwave_sc: Converged at iteration 193
    hwave_sc: Iteration result: eigenvalue = 0.959725, converged = True, n_iter = 193

followed by the eigenvalue analysis:

.. code-block:: text

    hwave_sc: === Eigenvalue analysis ===
    hwave_sc: Leading eigenvalues:
    hwave_sc:     0: -1.345390 (|ev| = 1.345390)
    hwave_sc:     1: -1.121834 (|ev| = 1.121834)
    hwave_sc:     2: -1.100188 (|ev| = 1.100188)
    ...
    hwave_sc:     4: 0.922684 (|ev| = 0.922684)

An eigenvalue :math:`\lambda > 1` (positive) indicates a superconducting instability
at the given temperature. Negative eigenvalues correspond to sign-changing gap symmetries
but do not indicate SC instability.


Results
----------------------------

Gap function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following figures show the gap function :math:`\Sigma_{\alpha\beta}(\mathbf{k})`
in momentum space obtained from the self-consistent iteration.

**Singlet channel** (:math:`\lambda \approx 0.96`):

.. figure:: ../sample_sc/gap_singlet.png
   :width: 90%
   :align: center

   Singlet gap function in k-space. Left: intra-orbital component
   :math:`\mathrm{Re}\,\Sigma_{00}(\mathbf{k})`.
   Right: inter-orbital component :math:`\mathrm{Re}\,\Sigma_{01}(\mathbf{k})`.
   The inter-orbital component is about 5 times larger than the intra-orbital one,
   indicating that inter-orbital pairing is dominant.

**Triplet channel** (:math:`\lambda \approx 1.58`):

.. figure:: ../sample_sc/gap_triplet.png
   :width: 90%
   :align: center

   Triplet gap function in k-space. Left: intra-orbital component
   :math:`\mathrm{Re}\,\Sigma_{00}(\mathbf{k})`.
   Right: inter-orbital component :math:`\mathrm{Re}\,\Sigma_{01}(\mathbf{k})`.
   Unlike the singlet case, the intra-orbital component dominates.

Eigenvalue spectrum
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The eigenvalue spectrum of the linearized Eliashberg equation is shown below
for both singlet and triplet channels.

.. figure:: ../sample_sc/eigenvalue_spectrum.png
   :width: 70%
   :align: center

   Positive eigenvalue spectrum :math:`\lambda` of the linearized Eliashberg equation.
   The dashed red line indicates :math:`\lambda = 1` (SC instability criterion).
   All singlet eigenvalues are below 1, while two triplet eigenvalues exceed 1.

The Arnoldi eigenvalue analysis finds multiple eigenvalues.
The figure shows only positive eigenvalues, which are relevant for the
SC instability criterion :math:`\lambda = 1`.
For the singlet channel, the largest positive eigenvalue is
:math:`\lambda \approx 0.92 < 1`, consistent with the
self-consistent iteration result (:math:`\lambda \approx 0.96`),
meaning no singlet SC instability at this temperature.

In the triplet channel, the leading positive eigenvalue is
:math:`\lambda \approx 1.58 > 1`, indicating a triplet SC instability.

Plotting script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The above figures can be reproduced using the plotting script included
in the sample directory:

.. code-block:: bash

    $ python plot_results.py


Output files
----------------------------

The solver produces the following files in the ``output`` directory:

``gap.dat``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The converged gap function :math:`\Delta_{\alpha\beta}(\mathbf{k})` in k-space.
Each line contains:

.. code-block:: text

    kx  ky  kz  Re(Δ_00)  Im(Δ_00)  Re(Δ_01)  Im(Δ_01)  Re(Δ_10)  Im(Δ_10)  Re(Δ_11)  Im(Δ_11)

where :math:`\alpha, \beta` are orbital indices.

``eigenvalue.dat``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The eigenvalues of the linearized Eliashberg equation:

.. code-block:: text

    # Iteration eigenvalue
    9.59724792e-01
    # Eigenvalue analysis
    # index  Re(eigenvalue)  Im(eigenvalue)  |eigenvalue|
       0 -1.34539047e+00  0.00000000e+00  1.34539047e+00
       1 -1.12183387e+00  0.00000000e+00  1.12183387e+00
       ...


Physical interpretation
----------------------------

The leading eigenvalue :math:`\lambda` of the linearized Eliashberg equation
determines whether superconductivity emerges:

- :math:`\lambda > 1`: The normal state is unstable toward superconductivity.
  The corresponding eigenvector gives the symmetry of the gap function.
- :math:`\lambda < 1` (for all positive eigenvalues): The normal state is stable
  at this temperature.

Negative eigenvalues correspond to sign-changing gap functions
(such as :math:`s_\pm`-wave in multi-orbital systems),
but they do **not** indicate SC instability regardless of their magnitude.
The self-consistency condition :math:`\Delta = K\Delta` requires :math:`\lambda = 1`
(not :math:`\lambda = -1`).

By varying the temperature and finding the point where the largest
positive eigenvalue reaches :math:`\lambda = 1`, one can determine
the superconducting transition temperature :math:`T_c`.

In this example, the singlet channel has :math:`\lambda_S \approx 0.96 < 1`
(no SC instability), while the triplet channel has
:math:`\lambda_T \approx 1.58 > 1` (SC instability) at :math:`T = 0.1`.

Singlet vs. triplet comparison
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By changing ``pairing_type`` to ``"triplet"`` in the input file,
one can compare singlet and triplet instabilities.
At :math:`T = 0.1` with the same parameters,
the triplet channel yields a leading eigenvalue
:math:`\lambda_T \approx 1.58`, which is larger than the singlet value
(:math:`\lambda_S \approx 0.96`).
This indicates that the triplet SC state is dominant at this temperature,
consistent with the results of Ref. [1]_ where the triplet
SC state competes with the singlet SC state for :math:`T > 0.05`,
while the singlet SC transition dominates at lower temperatures
(:math:`T < 0.05`) due to the enhancement of spin fluctuations.


Supported interactions
----------------------------

The Eliashberg solver supports all interaction types available in H-wave:

- ``CoulombIntra`` (:math:`U`): Intra-orbital Coulomb repulsion
- ``CoulombInter`` (:math:`V`): Inter-orbital Coulomb repulsion
- ``Hund`` (:math:`J`): Hund's coupling
- ``Exchange`` (:math:`J'`): Pair-hopping (exchange)
- ``Ising`` (:math:`I`): Ising-type spin interaction
- ``PairHop`` (:math:`P`): Pair hopping

For systems with ``Hund``, ``Exchange``, ``Ising``, or ``PairHop``
interactions, the solver automatically uses the generalized
:math:`S`/:math:`C` matrix formulation (Kuroki et al., PRB 79, 224511)
with four-index vertex structure.


Tips
----------------------------

- For large systems, set ``chi0q_mode = "calc"`` to compute
  :math:`\chi_0(\mathbf{q})` internally and avoid loading
  a large file.
- The ``"arnoldi"`` eigenvalue method is fastest for finding
  a few leading eigenvalues. For degenerate eigenvalues,
  ``"subspace"`` may be more robust.
- Different ``init_gap`` symmetries can be used to target
  specific pairing channels in iteration mode.
  The eigenvalue mode finds all leading symmetries automatically.
- The ``pairing_type = "triplet"`` option analyzes triplet
  pairing instabilities using the appropriate vertex.
