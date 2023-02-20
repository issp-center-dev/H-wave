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
  Set to ``"RPA"`` for calculations of the Random Phase Approximation.


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
  When ``Ncond`` is specified, the ``filling`` parameter will not be used.

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
  It must be greater than or equal to zero.

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
  The input file is in NumPy zip format that corresponds to the output format of
  ``green`` of the UHF calculation.
  If the initial configuration is not given, it is assumed to be zero.

- ``chi0q_init``

  **Type :**
  String

  **Description :**
  This parameter specifies the filename of the pre-calculated irreducible susceptibility
  :math:`\chi_0(\vec{q})` to be used for the calculation of the susceptibility matrix.
  The input file is in NumPy binary format that corresponds to the output format of
  ``chi0q`` in ``file.output`` section.

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

- ``matsubara_frequency``

  **Type :**
  Integer, List of Integers, or String (default value is ``"all"``)

  **Description :**
  This parameter specifies the range of Matsubara frequency that the susceptibility matrices :math:`\chi(\vec{q})` and :math:`\chi_0(\vec{q})` will be stored at.
  The value must be one of the following:

    - *an integer value* : a single index value.
    - ``[`` *min*, *max* (, *step*) ``]`` : every *step* index from *min* to *max*. If *step* is omitted, it is assumed to be 1.
    - all : all indices
    - center : corresponds to ``Nmat/2``.
    - none : nothing will be stored.
