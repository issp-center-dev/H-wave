.. highlight:: none

.. _Ch:Config:

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
  Set to ``"UHFk"`` for calculations of the wave-number space UHF.


- ``flag_fock``

  **Type :**
  Boolean (Default value is ``true``)

  **Description :**
  If ``true``, include the Fock term in the Hamiltonian, otherwise, exclude it (Hartree approximation).

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
  It must be between 0 and 1.
  Both ``Ncond`` and ``filling`` are specified, the program will be terminated with error.

- ``Ncond_round_mode``

  **Type :**
  String (default value is ``"strict"``)

  **Description :**
  This parameter specifies how the number of electrons calculated from the ``filling`` parameter is rounded to an integer value. The parameter must take one of the following values.

    - ``as-is``:  the value is not rounded to an integer. (returns a floating-point number)
    - ``round-up``:  the value is rounded up.
    - ``round-down``:  the value is rounded down.
    - ``round-off``:  the value is rounded to the closest integer. (0.5 is rounded up.)
    - ``round``:  the value is rounded by ``round`` function. (0.5 is rounded down.)
    - ``strict``:  if the value is not an integer value, the program terminates with error.
    - ``exact``:  if the value is not an integer value, a warning message will be shown and the value is rounded to an integer as ``round``.

- ``IterationMax``

  **Type :**
  Integer (default value is 20000)

  **Description :**
  This parameter specifies the maximum number of iterations.
  It must be greater than or equal to zero.

- ``EPS``

  **Type :**
  Integer (default value is 6)

  **Description :**
  This parameter specifies the convergence criterion.
  It is considered convergent when the norm of the difference between the previous and the new Green's function falls below :math:`10^{\rm -EPS}`.
  The residue is defined by
  :math:`R = \sum_{i, j}^{N} \sqrt{ \left| G_{ij}^{\rm new} - G_{ij}^{\rm old} \right|^2} / 2N^2`.
  It must be greater than or equal to zero.

- ``Mix``

  **Type :**
  Float (default value is 0.5)

  **Description :**
  This parameter specifies the ratio :math:`\alpha` of simple-mixing
  when the Green's function is updated by the previous and the new one.
  It must be between 0 and 1.
  If it is set to 1, the previous value will not be mixed.

  See :ref:`Algorithms <algorithm_sec>` section for the simple-mixing algorithm.

- ``RndSeed``

  **Type :**
  Integer (default value is 1234)

  **Description :**
  This parameter specifies the seed of the random number generator.

- ``ene_cutoff``

  **Type :**
  Float (default value is 100.0)

  **Description :**
  This parameter specifies the cutoff to avoid overflow in the calculations of the Fermi distribution function.

- ``strict_hermite``

  **Type :**
  Boolean (default value is false)

  **Description :**
  This parameter specifies the strictness of checking Hermiticity of the interaction definitions when they are read from files.
  If it is true, the program will be terminated with error when there are deviations greater than ``hermite_tolerance``. If it is false, only the warning messages will be shown, and the calculation continues.

- ``hermite_tolerance``

  **Type :**
  Float (default value is :math:`10^{-8}`)

  **Description :**
  This parameter specifies the tolerance for the Hermiticity condition
  :math:`|t_{ij} - t_{ji}^*| < \varepsilon`.

- ``trustme_interaction_range``

  **Type :**
  Boolean (default value is false)

  **Description :**

  Relax the checks that the distances of hopping and interaction are within the halves of ``CellShape``;
  if true, warn and continue execution, otherwise stop.

``log`` section
================================

- ``print_level``

  **Type :**
  Integer (default value is 1)

  **Description :**
  This parameter specifies the verbosity of the standard output.
  If it is set to 1, a detailed information will be shown.

- ``print_step``

  **Type :**
  Integer (default value is 1)

  **Description :**
  This parameter specifies the interval to write calculation logs to the standard output
  during the iteration.
  It must be greater than or equal to one.

- ``print_check``

  **Type :**
  String (default value is None)

  **Description :**
  This parameter specifies the filename to which the calculation logs are written
  during the iteration besides the standard output.
  If it is not set, no output file will be generated.


``file`` section
================================

This section consists of ``input`` and ``output`` subsections.
The former specifies the settings on the input files (such as locations and name of the files),
and the latter on the output files (such as locations to store).

``file.input`` section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``path_to_input``

  **Type :**
  String (default value is ``""``)

  **Description :**
  This parameter specifies the directory in which the input files are located.

- ``initial``

  **Type :**
  String

  **Description :**
  This parameter specifies the filename of the initial configuration
  of the one-body Green's function.
  The input file is in NumPy binary format that corresponds to the output format of
  ``green`` in ``file.output`` section.

- ``initial_mode``

  **Type :**
  String

  **Description :**
  This parameter specifies how the initial value of the Green's function is set when the initial configuration is not given from files. It takes one of the following:

  - ``zero``: set the initial configuration to be zero. (default)

  - ``one`` or ``unity``: set the initial configuration as :math:`G_{\alpha\sigma,\beta\sigma^\prime}(\vec{r}) = \delta_{\vec{r},0}\delta_{\alpha\beta}\delta_{\sigma\sigma^\prime}`.
  
  - ``random``: set the initial configuration by random numbers.

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

- ``Transfer``, ``CoulombIntra``, ``CoulombInter``, ``Hund``, ``Ising``, ``Exchange``, ``PairLift``, ``PairHop``

  **Type :**
  String

  **Description :**
  These parameters specify the filenames for the definitions of the corresponding interaction terms.

``file.output`` section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``path_to_output``

  **Type :**
  String (default value is ``"output"``)

  **Description :**
  This parameter specifies the directory in which the output files are stored.

- ``energy``

  **Type :**
  String

  **Description :**

  This parameter specifies the name of the file to store the energy values.
  If it is not set, no output will be generated.

- ``eigen``

  **Type :**
  String

  **Description :**
  This parameter specifies the name of the file to store the eigenvalues and eigenvectors
  of the Hamiltonian.
  If it is not set, no output will be generated.

- ``green``

  **Type :**
  String (default value is ``"green"``)

  **Description :**
  This parameter specifies the name of the file to store the one-body Green's functions.
  If it is set to an empty string, no output will be generated.

- ``rpa``

  **Type :**
  String

  **Description :**
  This parameter specifies the name of the file to store the one-body term involving the approximated two-body interaction term by UHF method. The output may be used as an initial configuration of the RPA calculation by specifying it in the parameter ``file.input.trans_mod``.
  If it is not set, no output will be generated.
