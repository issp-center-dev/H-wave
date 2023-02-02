.. highlight:: none

.. _Ch:Config_UHFR:

Parameter file
--------------------------------

The parameter file specifies calculation conditions of H-wave in TOML format.
This file consists of these three sections:

#. ``mode`` section to set calculation mode,

#. ``log`` section to set standard log output,

#. ``file`` section to set file and directory paths.
   This section consists of two subsections, ``input`` and ``output``.

A sample file reads as follows:

::

    [log]
    print_level = 1
    print_step = 20
    [mode]
    mode = "UHFr"
    [mode.param]
    Nsite = 8
    2Sz = 0
    Ncond = 8
    IterationMax = 1000
    EPS = 8
    RndSeed = 123456789
    T = 0.0
    [file]
    [file.input]
    path_to_input = ""
    OneBodyG = "greenone.def"
    [file.input.interaction]
    Trans = "trans.def"
    CoulombIntra = "coulombintra.def"
    [file.output]
    path_to_output = "output"
    energy = "energy.dat"
    eigen = "eigen"
    green = "green.dat"

File format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
TOML


Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``mode`` section
================================

- ``mode``

  **Type :**
  String

  **Description :**
  This parameter specifies calculation mode.
  Set to ``UHFr`` when using coordinate-space UHF.

- ``flag_fock``

  **Type :**
  Boolean (default value is ``true``)

  **Description :**
  When this parameter is ``true``, the Fock term is considered.
  If it is ``false``, only the Hartree term is considered.


``mode.param`` section
================================

``mode.param`` section sets parameters of the calculation.

- ``T``

  **Type :**
  Float (default value is ``0``)

  **Description :**
  This parameter specifies temperature. It must be greater than or equal to zero.

- ``2Sz``

  **Type :**
  Integer, String, or None (default value is None)

  **Description :**
  Twice the size of the :math:`z` compoonent of the total spin is specified
  when it is set to a fixed value.
  In this case, the up and down spin components are calculated separately.
  If this parameter is not given, or it is set to ``"free"``, the spin space is
  not separated in the calculation.

  This parameter should not be specified when ``Sz`` is not conserved (e.g. when
  the spin-orbital interaction is present).

  ``2Sz`` takes a value between ``-Nsite`` and ``Nsite``.

- ``Nsite``

  **Type :**
  Integer

  **Description :**
  This parameter specifies the number of sites. It must be greater than or equal to one.

- ``Ncond``

  **Type :**
  Integer

  **Description :**
  This parameter specifies the number of conduction electrons. It must be greater than or equal to one.


- ``IterationMax``

  **Type :**
  Integer (default value is ``20000``)

  **Description :**
  This parameter specifies the maximum number of iterations. It must be greater than or equal to zero.


- ``EPS``

  **Type :**
  Integer (default value is ``6``)
  
  **Description :**
  This parameter specifies the convergence criterion.
  The solver iteration will be terminated when
  the norm of the difference between the previous and new Green's function falls
  below :math:`10^{\rm -EPS}`.
  The residue is defined by
  :math:`R = \sum_{i,j}^{N}\sqrt{ \left| G_{ij}^{\rm new} - G_{ij}^{\rm old} \right|^2} / 2N^2`.
  It must be greater than or equal to zero.

- ``Mix``

  **Type :**
  Float (default value is ``0.5``)

  **Description :**
  This parameter specifies the ratio :math:`\alpha` of simple-mixing
  when the Green's function is updated by the previous and the new one.
  It must be between 0 and 1.
  If it is set to 1, the previous value will not be mixed.
  See :ref:`Algorithms <algorithm_sec>` section for simple-mixing algorithm.

- ``RndSeed``

  **Type :**
  Integer (default value is ``1234``)

  **Description :**
  This parameter specifies the seed of random numbers.

- ``ene_cutoff``

  **Type :**
  Float (default value is ``100.0``)
  
  **Description :**
  This parameter specifies a cut-off to avoid overflow when the Fermi distribution function is calculated.

- ``strict_hermite``

  **Type :**
  Boolean (default value is ``false``)

  **Description :**
  This parameter specifies strictness of Hermiticity checks when the interaction definitions are read from files.
  If it is set to ``true``, the program stops when the deviation larger than ``hermite_tolerance`` is detected.
  If it is set to ``false``, a warning message will be shown and the program execution continues. 

- ``hermite_tolerance``

  **Type :**
  Float (default value is :math:`10^{-8}`)

  **Description :**
  This parameter specifies the tolerance of the deviation from Hermiticity condition
  :math:`|t_{ij} - t_{ji}^*| < \varepsilon`.

``log`` section
================================

- ``print_level``

  **Type :**
  Integer (default value is ``1``)

  **Description :**
  This parameter specifies verbosity of the standard log output.
  When it is set to ``1``, the detailed information will be printed.

- ``print_step``

  **Type :**
  Integer (default value is ``1``)
  
  **Description :**
  This parameter specifies the interval between outputs of calculation logs to the standard output during iterations.
  It must be greater than or equal to one.

- ``print_check``

  **Type :**
  String

  **Description :**
  This parameter specifies the output logfile to which the calculation logs are written during the iterations besides the standard output.
  If it is not given, the logs are not exported to files.

``file`` section
================================

This section consists of ``input`` and ``output`` subsections.
The former specifies settings on input files (e.g. locations and names of files),
while the latter on output files, as described below.

``file.input`` section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``path_to_input``

  **Type :**
  String (default value is ``""``)

  **Description :**
  This parameter specifies the directory in which the input files are located.

- ``Initial``

  **Type :**
  String (default value is ``""``)

  **Description :**
  This parameter specifies the input file for the initial configuration.

- ``OneBodyG``

  **Type :**
  String (default value is ``""``)

  **Description :**
  This parameter specifies the input file that contains a list of indices of one-body Green's function to export.
  

``file.input.interaction`` section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``Trans``

  **Type :**
  String (default value is ``""``)

  **Description :**
  This parameter specifies the input file for general one-body interaction term.

- ``InterAll``

  **Type :**
  String (default value is ``""``)

  **Description :**
  This parameter specifies the input file for generalized two-body interaction term.

- ``CoulombIntra``

  **Type :**
  String (default value is ``""``)

  **Description :**
  This parameter specifies the input file for on-site Coulomb interaction term.

- ``CoulombInter``

  **Type :**
  String (default value is ``""``)

  **Description :**
  This parameter specifies the input file for inter-site Coulomb interaction term.

- ``Hund``

  **Type :**
  String (default value is ``""``)

  **Description :**
  This parameter specifies the input file for Hund interaction term.

- ``PairHop``

  **Type :**
  String (default value is ``""``)

  **Description :**
  This parameter specifies the input file for pair-hopping term.

- ``Exchange``

  **Type :**
  String (default value is ``""``)

  **Description :**
  This parameter specifies the input file for exchange interaction term.

- ``Ising``

  **Type :**
  String (default value is ``""``)

  **Description :**
  This parameter specifies the input file for Ising interaction term.

- ``PairLift``

  **Type :**
  String (default value is ``""``)

  **Description :**
  This parameter specifies the input file for pair-lift interaction term.

``file.output`` section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``path_to_output``

  **Type :**
  String (default value is ``"output"``)

  **Description :**
  This parameter specifies the directory to store the output files.

- ``energy``

  **Type :**
  String

  **Description :**
  This parameter specifies the output file for energies.
  If it is not given, the output is not exported.

- ``eigen``

  **Type :**
  String

  **Description :**
  This parameter specifies the output file for eigenvalues of Hamiltonian.
  If it is not given, the output is not exported.

- ``green``

  **Type :**
  String

  **Description :**
  This parameter specifies the output file for one-body Green's function.
  If it is not given, the output is not exported.

- ``initial``

  **Type :**
  String

  **Description :**
  This parameter specifies the output file for one-body Green's function in a format suitable for the initial configuration.
  If it is not given, the output is not exported.

- ``fij``

  **Type :**
  String

  **Description :**
  This parameter specifies the output file for pair-orbital factor :math:`f_{ij}`.
  If it is not given, the output is not exported.
