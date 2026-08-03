.. highlight:: none

.. _Ch:Config_rpa:

Parameter files
--------------------------------

The parameter file specifies calculation conditions and parameters for H-wave
in TOML format. It is composed of the following three sections.

#. ``mode`` section for specifying calculation conditions,

#. ``log`` section for setting standard outputs,

#. ``file`` section for setting file paths: It contains ``input`` and ``output`` subsections.

An example of the file is shown below:

.. literalinclude:: ../../sample/input.toml


File format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
TOML format


Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``mode`` section
================================

- ``mode``

  **Type :**
  String

  **Description :**
  This parameter specifies the calculation mode.
  Set to ``"RPA"`` for calculations of the Random Phase Approximation,
  or to ``"FLEX"`` for the Fluctuation Exchange Approximation.

- ``enable_spin_orbital`` (default value is ``false``)

  **Type :**
  Boolean

  **Description :**
  This parameter specifies whether to allow spin-orbital interaction.
  If it is set to ``true``, the orbital indices in the Transfer term use the interleaved
  spin-orbital convention: the combined index is :math:`2 \alpha + s`, where
  :math:`\alpha` is the physical-orbital index (0-based) and :math:`s \in \{0, 1\}`
  is the spin index.

  In spin-orbital mode, the ``Norbit`` value in the geometry file (``geom.dat``) is
  the **spin-orbital count** (= 2 × the number of physical orbitals = Wannier90
  ``num_wann``), the same convention as UHFk.

  .. note::

     **Migration (RPA):** the geometry ``norb`` for spin-orbital input is now the
     spin-orbital count; double any pre-existing RPA spin-orbital ``geom.dat``
     ``Norbit``.

- ``calc_scheme`` (default value is ``"auto"``)

  **Type :**
  String

  **Description :**
  This parameter specifies how the spin and orbitals are treated in the calculation. The parameter takes one of the following options.

  - ``general``: Generalized orbitals combining spins and orbitals are considered. The susceptibility matrix takes the most general form, with the size of :math:`N_\text{orb}^4 N_\text{spin}^4 N_k N_\omega`.

  - ``reduced``: Generalized orbitals combining spins and orbitals are considered. The components of the susceptibility matrix with :math:`\alpha=\alpha^\prime` and :math:`\beta=\beta^\prime` are considered. The size of the matrix turns to :math:`N_\text{orb}^2 N_\text{spin}^2 N_k N_\omega`. For the two-body interaction terms, only CoulombIntra, CoulombInter, Ising and Hund are allowed. 

  - ``squashed``: Spins and orbitals are separately treated, and for the orbitals :math:`\alpha=\alpha^\prime` and :math:`\beta=\beta^\prime` are considered. The size of the susceptilibity matrix becomes :math:`N_\text{orb}^2 N_\text{spin}^4 N_k N_\omega`. See :ref:`Ch:Algorithm` for details.

  - ``auto``: scheme is automatically chosen according to the specifications of interaction terms. This option is not available when only ``chi0q`` is to be calculated.

- ``calc_type`` (default value is ``"ring"``)

  **Type :**
  String

  **Description :**
  This parameter specifies which RPA diagrams to include.

  - ``ring``: Standard RPA (ring diagram only). Computes the longitudinal susceptibility.

  - ``ring+ladder``: Includes the transverse (ladder) susceptibility :math:`\chi_{+-}(\mathbf{q})` in addition to the standard RPA. This requires the ``general`` calculation scheme (automatically selected). See :ref:`Ch:Algorithm` for details.


``mode.param`` section
================================

``mode.param`` section contains the parameters for the calculation.

- ``CellShape``

  **Type :**
  Integer array

  **Description :**
  This parameter specifies the shape of the lattice Lx, Ly, Lz.

- ``SubShape``

  **Type :**
  Integer array (default value is ``[`` Lx, Ly, Lz ``]``)

  **Description :**
  This parameter specifies the shape of the sublattice Bx, By, Bz.

- ``T``

  **Type :**
  Float (default value is 0)

  **Description :**
  This parameter specifies the temperature.
  It must be greater than or equal to zero.

- ``mu``

  **Type :**
  Float or None (default value is None)

  **Description :**
  This parameter specifies the chemical potential :math:`\mu`.
  If it is not specified, the value of :math:`\mu` will be calculated so that
  the expectation value of the number of electrons equals to ``Ncond``.
  If both ``mu`` and ``Ncond`` or ``filling`` are specified, the program terminates with error.

- ``Ncond``

  **Type :**
  Integer

  **Description :**
  This parameter specifies the number of conduction electrons.
  It must be greater than or equal to one.

- ``filling``

  **Type :**
  Float

  **Description :**
  This parameter specifies the filling ratio of electrons with respect to the number of states.
  Both ``Ncond`` and ``filling`` are specified, the program will be terminated with error.

- ``Ncond_round_mode``

  **Type :**
  String (default value is ``"strict"``)

  **Description :**
  This parameter specifies how the number of electrons calculated from the ``filling`` parameter is rounded to an integer value when the temperature is zero. The parameter must take one of the following values.

    - ``as-is``:  the value is not rounded to an integer. (returns a floating-point number)
    - ``round-up``:  the value is rounded up.
    - ``round-down``:  the value is rounded down.
    - ``round-off``:  the value is rounded to the closest integer. (0.5 is rounded up.)
    - ``round``:  the value is rounded by ``round`` function. (0.5 is rounded down.)
    - ``strict``:  if the value is not an integer value, the program terminates with error.
    - ``exact``:  if the value is not an integer value, a warning message will be shown and the value is rounded to an integer as ``round``.

- ``Nmat``

  **Type :**
  Integer (default value is 1024)

  **Description :**
  This parameter specifies the cut-off of Matsubara frequency.
  It must be an even number greater than zero. Matsubara frequency is defined as follows:

      - Boson: :math:`\omega_n = \dfrac{2\pi (n-\texttt{Nmat}/2)}{\beta}`
      - Fermion: :math:`\omega_n = \dfrac{\pi (2n+1-\texttt{Nmat})}{\beta}`

  with the indices :math:`n` between 0 and ``Nmat-1``.

- ``coeff_tail``

  **Type :**
  Float (default value is 0.0)

  **Description :**
  This parameter specifies the magnitude of the correction when correcting the tails of the Fourier transformation.
  After Fourier transforming the diagonalized one-body Green function to the imaginary time representation by subtracting :math:`\texttt{coeff\_tail}/(i \omega_n)`, the term :math:`-\beta/2\cdot\texttt{coeff\_tail}` is added to the one-body Green function.
  In the FLEX solver the same tail treatment is applied to the *dressed* Green function before the bare susceptibility :math:`\chi_0(q)` is computed, so that ``coeff_tail`` accelerates the frequency summation without changing the physical result. (The FLEX self-energy convolution keeps the full Green function and is unaffected.)
  Since issue #134 the susceptibility kernels also restore the Green
  function's equal-time discontinuity at the bubble's :math:`\tau = 0`
  sample (the tail piece carries the jump; the sample is the mean of the
  two branches), which makes ``coeff_tail = 1.0`` converge at
  :math:`O(1/N_{\rm mat}^2)`. Earlier versions omitted this endpoint and
  ``coeff_tail`` then *slowed* the convergence by a constant factor
  (still :math:`O(1/N_{\rm mat})`); results produced with a nonzero
  ``coeff_tail`` before this fix are not comparable with current ones
  at the same ``Nmat``.
  The value must be a finite real number; ``NaN`` and infinities are
  rejected. Only ``0.0`` (off) and ``1.0`` (the physical :math:`1/i\omega_n`
  coefficient) are recommended: fractional values cancel only part of the
  equal-time jump, remain :math:`O(1/N_{\rm mat})` and can converge more
  slowly than ``coeff_tail = 0.0``.
  ``chi0q.npz`` files written with a nonzero ``coeff_tail`` carry a
  ``tail_endpoint = "branch_mean_v1"`` marker recording the endpoint
  treatment; ``chi0q_init`` and the ``hwave_sc`` chi0q loader refuse a
  nonzero-tail file without it (produced before the fix), since the
  pre-fix error cannot be detected from the array itself. Recompute such
  bubbles instead of reusing the files.

- ``spinful_vertex_exchange``

  **Type :**
  Boolean (default value is ``true``)

  **Description :**
  Spinful (``enable_spin_orbital``) calculations resum the susceptibility
  with a single vertex tensor. Since issue #137 that tensor is the
  antisymmetrized bare particle-hole vertex: the direct (ring) wiring
  plus the exchange wiring of the on-site interaction terms. The
  exchange part is what corrects the spin-flip pair components of
  :math:`\chi(q)` (in non-spin-orbital calculations the analogous
  content is provided separately by ``calc_type = "ring+ladder"``);
  without it those components are returned as the bare bubble at any
  interaction strength, and, because spin is not conserved, the error
  leaks into every component. The construction was verified against
  exact diagonalization at first order in the coupling for
  CoulombIntra, Exchange and PairLift, and reproduces the established
  transverse (ring+ladder) series in the spin-conserving limit.
  Setting ``spinful_vertex_exchange = false`` restores the previous
  ring-only vertex (results produced before this fix): use it only to
  reproduce old runs. The exchange wiring of OFF-site interaction
  terms depends on both fermionic momenta and is not representable in
  this resummation; it remains excluded (as it is in the
  non-spin-orbital ladder).

- ``matsubara_frequency``

  **Type :**
  Integer, List of Integers, or String (default value is ``"all"``)

  **Description :**
  This parameter specifies the indices of Matsubara frequency for which the susceptibility matrix :math:`\chi(\vec{q})` is calculated.
  The value must be one of the following:

    - *an integer value* : a single index value.
    - ``[`` *min*, *max* (, *step*) ``]`` : every *step* index from *min* to *max*. If *step* is omitted, it is assumed to be 1.
    - all : all indices
    - center : corresponds to ``Nmat/2``.
    - none : nothing will be calculated.

  When the susceptibility matrix :math:`\chi(\vec{q})` or the irreducible susceptibility matrix :math:`\chi_0(\vec{q})` are stored to files, the values at the specified freqneucy are exported.


- ``coeff_extern``

  **Type :**
  Float (default value is 0.0)

  **Description :**
  This parameter specifies the coefficient :math:`h` of the external field given by the form :math:`\pm h H_{\alpha\beta}(r_{ij})`. The definition of the matrix :math:`H_{\alpha\beta}(r_{ij})` will be provided by an input file. The sign :math:`+` and :math:`-` correspond to spin up and down, respectively.
  


- ``RndSeed``

  **Type :**
  Integer (default value is 1234)

  **Description :**
  This parameter specifies the seed of the random number generator.

- ``ene_cutoff``

  **Type :**
  Float (default value is 100.0)

  **Description :**
  This parameter specifies the upper cutoff of the exponent in the Fermi distribution function to avoid overflow during the calculation.

- ``gpu``

  **Type :**
  Boolean (default value is false)

  **Description :**
  When set to true, the main computation (the Green's function, the chi0q
  FFT pair bubble, the spin inflation, and the batched RPA solve; for FLEX
  the whole SCF loop) runs on a GPU via CuPy. When CuPy or a CUDA device is
  unavailable, the solver warns and falls back to the CPU (numpy) path with
  an identical result. Install CuPy as the precompiled binary wheel matching
  your CUDA version (e.g. ``pip install cupy-cuda12x`` for CUDA 12.x); see
  the `CuPy installation guide <https://docs.cupy.dev/en/stable/install.html>`_.

- ``fft_workers``

  **Type :**
  Integer (default value is 1)

  **Description :**
  Number of worker threads for the spatial FFTs (parallelized via
  ``scipy.fft``). The default 1 keeps the serial numpy path, unchanged
  from previous releases (opt-in); -1 uses all cores. Ignored on the
  GPU. Set a smaller number when running several calculations
  concurrently.

- ``mixing_scheme``

  **Type :**
  String (default value is "linear"; FLEX mode only)

  **Description :**
  Self-energy update scheme of the FLEX SCF loop. ``"linear"`` is the
  conventional linear mixing; ``"anderson"`` enables Anderson acceleration
  (Pulay/DIIS-type extrapolation over a short iterate/residual history),
  which reaches the same fixed point in far fewer iterations. Anderson
  acceleration is less sensitive to step-size instability than linear
  mixing, so a somewhat larger ``Mix`` (e.g. 0.3--0.5) can reduce the
  iteration count further. Falls back to a plain linear step automatically
  if the history becomes degenerate.

- ``anderson_depth``

  **Type :**
  Integer (default value is 5; FLEX mode only)

  **Description :**
  History depth of the Anderson acceleration. Memory grows by 2*depth
  sigma-sized arrays (kept on the device under GPU execution).


``log`` section
================================

- ``print_level``

  **Type :**
  Integer (default value is 1)

  **Description :**
  This parameter specifies the verbosity of the standard output.
  If it is set to 1, a detailed information will be shown.


``file`` section
================================

This section consists of ``input`` and ``output`` subsections.
They specify the settings of the input and output files, respectively, on the types of files, the directories to be located or stored, and the names of the files.


``file.input`` section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``path_to_input``

  **Type :**
  String (default value is ``""`` (blank string))

  **Description :**
  This parameter specifies the directory in which the input files are located.

- ``chi0q_init``

  **Type :**
  String

  **Description :**
  This parameter specifies the filename of the pre-calculated irreducible susceptibility
  :math:`\chi_0(\vec{q})` to be used for the calculation of the susceptibility matrix.
  The input file is in NumPy binary format that corresponds to the output format of
  ``chi0q`` in ``file.output`` section.

- ``trans_mod``

  **Type :**
  String

  **Description :**
  This parameter specifies the filename of the initial configuration exported from UHFk by the parameter ``file.output.rpa``. It contains the one-body interaction term involving the approximated two-body interaction terms via UHF method.

- ``green_init``

  **Type :**
  String

  **Description :**
  This parameter specifies the filename of the initial Green's function for RPA calculation. The file format corresponds to the output file of ``green`` of UHFk. When ``trans_mod`` is specified, ``green_init`` is not used.

- ``sigma_init``

  **Type :**
  String (FLEX mode only)

  **Description :**
  This parameter specifies the filename of a ``sigma.npz`` (written by an
  earlier FLEX run) used to seed the FLEX SCF loop instead of Sigma = 0. The
  path is resolved relative to ``path_to_input``. The recorded ``CellShape``
  and the array's ``Nmat`` must match the current run (a mismatch is a
  fail-fast error). See the FLEX tutorial section "Warm-starting the SCF
  loop" for the sweep workflow.


``file.input.interaction`` section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This section describes the relation of the interaction types and geometry information
to the definition files.

- ``path_to_input``

  **Type :**
  String

  **Description :**
  This parameter specifies the directory in which the input files are located.
  It is independent from ``path_to_input`` in ``file.input`` section.

- ``Geometry``

  **Type :**
  String

  **Description :**
  This parameter specifies the filename for the geometry information.

- ``Transfer``, ``CoulombIntra``, ``CoulombInter``, ``Hund``, ``Ising``, ``Exchange``, ``PairLift``, ``PairHop``, ``Extern``

  **Type :**
  String

  **Description :**
  These parameters specify the filenames for the definitions of the corresponding interaction terms. If none of two-body interaction term (CoulombIntra, CoulombInter, Hund, Ising, Exchange, PairLift, or PairHop) is specified, the program only calculates ``chi0q`` and exits.


``file.output`` section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``path_to_output``

  **Type :**
  String (default value is ``"output"``)

  **Description :**
  This parameter specifies the directory in which the output files are stored.

- ``chi0q``

  **Type :**
  String

  **Description :**

  This parameter specifies the name of the file to store the irreducible susceptibility matrix
  :math:`\chi_0(\vec{q})`.
  If it is not set, no output will be generated.

- ``chiq``

  **Type :**
  String

  **Description :**
  This parameter specifies the name of the file to store the susceptibility matrix
  :math:`\chi(\vec{q})`.
  If it is not set, no output will be generated.
