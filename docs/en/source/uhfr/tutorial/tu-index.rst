==================
Tutorial
==================

To use H-wave, you need to prepare the input files:

#. parameter file to set calculation conditions, 
#. definition files of the Hamiltonian,
#. specifications to output the results, 

before performing calculations.
In the following, we will provide a tutorial
using a sample in ``docs/tutorial/Hubbard/UHFr`` directory.
   

Create a parameter file
------------------------------------------

The parameter file contains information to control inputs and outputs of the program.
An example is given in the directory ``docs/tutorial/Hubbard/UHFr`` by the name ``input.toml``.
The content of the file will be as follows:

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

The parameter file is described in TOML format.

In the ``log`` section, ``print_level`` specifies the level of the standard output, and ``print_step`` specifies the number of steps between printing logs.

In the ``mode`` section, the calculation mode and the basic parameters are specified.

In the ``file.input`` section, ``path_to_input`` specifies the directory in which input files are located, ``OneBodyG`` specifies the definition file of the one-body Green's function to output, and ``Initial`` specifies the initial configuration.
If ``OneBodyG`` is missing, the Green's function will not be exported.
When ``Initial`` is not specified, a random configuration will be generated for the initial state.

In the ``file.input.interaction`` section, the input files to define Hamiltonian are specified.

In the ``file.output`` section, ``path_to_output`` specifies the directory to which the results will be written. ``energy`` specifies the filename for the energy, ``eigen`` specifies the filename for the eigenvalues and eigenvectors of the Hamiltonian, and ``green`` specifies the filename for the one-body Green's function.
If these keywords are missing, the corresponding results will not be exported.

See File format section for the details.


Create definition files for Hamiltonian
---------------------------------------

Next, we will create input files that defines the Hamiltonian.

Transfer term
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The input file associated with a keyword ``Trans`` (``trans.def`` in this tutorial)
provides definitions of Hamltonian for Transfer term of the electron system:

.. math::

   \mathcal{H} = -\sum_{ij\sigma_1\sigma_2}
   t_{ij\sigma_1\sigma_2}c_{i\sigma_1}^{\dagger}c_{j\sigma_2}^{\phantom\dagger}.
   
The content of the file is as follows:

::

    ========================
    NTransfer      64
    ========================
    ========i_j_s_tijs======
    ========================
        4     0     0     0         1.000000000000000         0.000000000000000
        0     0     4     0         1.000000000000000        -0.000000000000000
        4     1     0     1         1.000000000000000         0.000000000000000
        0     1     4     1         1.000000000000000        -0.000000000000000
        2     0     0     0         1.000000000000000         0.000000000000000
        0     0     2     0         1.000000000000000        -0.000000000000000
        2     1     0     1         1.000000000000000         0.000000000000000
        0     1     2     1         1.000000000000000        -0.000000000000000
    ...


See :ref:`Subsec:Trans` for the details.

Two-body interaction term
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this tutorial, we consider a two-body interaction Hamiltonian of the electron system
of the form:

.. math::

   \mathcal{H} = \sum_{i} U_i n_{i\uparrow}n_{i\downarrow}.

The definition is given in the file associated with the keyword ``CoulombIntra``
(``coulombintra.def`` in the present case). 
The content of the file is as follows:
   
::

    =============================================
    NCoulombIntra          8
    =============================================
    ================== CoulombIntra ================
    =============================================
        0         8.000000000000000
        1         8.000000000000000
        2         8.000000000000000
        3         8.000000000000000
        4         8.000000000000000
     ...


There are a number of keywords provided to concicely describe the Hamiltonian,
besides ``CoulombIntra``.
See sections :ref:`Subsec:interall` - :ref:`Subsec:pairlift` for the details.

Specify output components
----------------------------

Next, we will provide the files that describe the output components.


Setting indices of one-body Green's functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A file associated with the keyword ``OneBodyG`` (``greenone.def`` in this tutorial) specifies
the indices of one-body Green's functions to be calculated
:math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2} \rangle`.
The content of the file will be as follows:

::

    ===============================
    NCisAjs         16
    ===============================
    ======== Green functions ======
    ===============================
        0     0     0     0
        0     0     1     0
        0     0     2     0
        0     0     3     0
        0     0     4     0
     ...

See :ref:`Subsec:onebodyg` for the details of the file format
to specify indices of the one-body Green's functions.

Run
--------------------------

All the input files have been created, and we are ready to run the program.
Type in the command with the parameter file (``input.toml`` in this tutorial) as an argument:

.. code-block:: bash

    $ hwave input.toml

The calculation is launched, and the logs will be shown as follows:

::

    2022-12-01 09:37:30,114 INFO qlms: Read def files
    2022-12-01 09:37:30,116 INFO qlms: Get Hamiltonian information
    2022-12-01 09:37:30,116 INFO qlms: Get Green function information
    2022-12-01 09:37:30,116 INFO qlms.uhfr: Show input parameters
      Nsite               : 8
      Ncond               : 8
      2Sz                 : 0
      Mix                 : 0.5
      EPS                 : 1e-08
      IterationMax        : 1000
      RndSeed             : 123456789
      T                   : 0.0
      ene_cutoff          : 100.0
      threshold           : 1e-12
    2022-12-01 09:37:30,117 INFO qlms: Start UHF calculation
    2022-12-01 09:37:30,117 INFO qlms.uhfr: Set Initial Green's functions
    2022-12-01 09:37:30,117 INFO qlms.uhfr: Initialize green function by random numbers
    2022-12-01 09:37:30,117 INFO qlms.uhfr: Start UHFr calculations
    2022-12-01 09:37:30,117 INFO qlms.uhfr: step, rest, energy, NCond, Sz
    2022-12-01 09:37:30,119 INFO qlms.uhfr: 0, 0.022144468, -27.16081+0j, 8, -7.425e-16
    2022-12-01 09:37:30,134 INFO qlms.uhfr: 20, 1.2083848e-05, -3.399532+0j, 8, -1.055e-15
    2022-12-01 09:37:30,145 INFO qlms.uhfr: UHFr calculation is succeeded: rest=5.7552848630056134e-09, eps=1e-08.
    2022-12-01 09:37:30,145 INFO qlms: Save calculation results.
    2022-12-01 09:37:30,146 INFO qlms: All procedures are finished.
    --------------------------------------------------------------------------------
    Statistics
      function                         :  total elapsed  : average elapsed : ncalls
    --------------------------------------------------------------------------------
      hwave.solver.uhfr.__init__       :      0.357 msec :      0.357 msec :      1
      hwave.solver.uhfr._initial_G     :      0.090 msec :      0.090 msec :      1
      hwave.solver.uhfr._makeham_const :      0.839 msec :      0.839 msec :      1
      hwave.solver.uhfr._makeham_mat   :      0.309 msec :      0.309 msec :      1
      hwave.solver.uhfr._makeham       :      6.001 msec :      0.176 msec :     34
      hwave.solver.uhfr._diag          :      2.468 msec :      0.073 msec :     34
      hwave.solver.uhfr._green         :      3.107 msec :      0.091 msec :     34
      hwave.solver.uhfr._calc_energy   :      1.990 msec :      0.059 msec :     34
      hwave.solver.uhfr._calc_phys     :     12.929 msec :      0.380 msec :     34
      hwave.solver.uhfr.solve          :     28.290 msec :     28.290 msec :      1
      hwave.solver.uhfr.save_results   :      0.852 msec :      0.852 msec :      1
    --------------------------------------------------------------------------------
		
The log messages on reading the input files are presented, followed by the information
on the process of UHF calculations.
The results are written in the ``output`` directory, according to the settings in ``file.output`` section of the input toml file:
``energy.dat`` for the eigenvalues,
``spin-up_eigen.npz`` and ``spin-down_eigen.npz`` for the eigenvectors, and
``green.dat`` for the one-body Green's functions.
See File format section for the details of the output files.

Generate input files using StdFace library
----------------------------------------------

The definition files for the Hamiltonian can be generated easily by using StdFace library.
StdFace is a program to create input files for the predefined lattice models or from the information based on the effective models provided in Wannier90 format.
It is used for the exact-diagonalization solver :math:`{\mathcal H}\Phi` and the many-variable variational Monte Carlo solver mVMC for the genration of input files.
In this section, we describe the procedure to create the input files for UHF using the StdFace library by a concrete example.

First, download the source files of StdFace library with git:

.. code-block:: bash

    $ git clone https://github.com/issp-center-dev/StdFace.git

Then, type in the following commands to compile the sources:

.. code-block:: bash

    $ cd StdFace
    $ mkdir build && cd build
    $ cmake -DUHF=ON ../
    $ make

If the compilation is successful, you will find an executable file ``uhf_dry.out`` in ``src`` directory.

In this tutorial, we use sample files in ``docs/tutorial/Hubbard/UHFr`` directory.
There is a file named ``stan.in`` which is an input to ``uhf_dry.out``.
The content of the file is as follows:

::

    model = "Hubbard"
    lattice = "square"
    a0W = 2
    a0L = 2
    a1W = -2
    a1L = 2
    t = 1.0
    U = 8.0
    ncond = 8
    2Sz = 0

- ``model`` is a keyword that specifies the model to consider. At present, only the value ``Hubbard`` is supported, which corresponds to the Hubbard model with the number of electrons fixed.
- ``lattice`` is a keyword that specifies the crystal structure. In this example, the square lattice ``square`` is chosen.
- ``a0W`` and ``a0L`` are parameters to set x axis as a vector ``(a0W, a0L)``, and ``a1W`` and ``a1L`` to set y axis as a vector ``(a1W, a1L)``.
- ``t`` corresponds to the hopping, and ``U`` to the on-site Coulomb interaction.
- ``ncond`` and ``2Sz`` are given for compatibility with :math:`{\mathcal H}\Phi` and mVMC. 
  It is noted that these parameters must also be specified in the parameter file. 

See, for example, the manual of :math:`{\mathcal H}\Phi` for the details of input parameters.

With the file shown above as an argument, run ``uhf_dry.out`` as follows:

.. code-block:: bash

    $ cd path_to_Hwave/docs/tutorial/Hubbard/UHFr
    $ ln -s path_to_Stdface/build/src/uhf_dry.out .
    $ ./uhf_dry.out stan.in

If the command is completed, a set of Hamiltonian definition files are generated in the current directory.
We use ``trans.def`` and ``coulombintra.def`` for the tutorial.
