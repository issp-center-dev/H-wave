==================
Tutorial
==================

To use H-wave in wave-number space (UHFk),
you need to prepare input files:

#. an input parameter file,
#. interaction definition files,

before running the program.
For the latter, you can use the outputs of external programs such as RESPACK,
or you may create the files using the StdFace library.

In the following, we provide a tutorial based on a sample in
``docs/tutorial/Hubbard/UHFk`` directory.
The interaction definition files are generated using StdFace library.


Create a parameter file
--------------------------------

We prepare an input parameter file that contains basic parameters as well as settings on
inputs and outputs.
A sample file can be found in ``docs/tutorial/Hubbard/UHFk`` directory by a filename
``input.toml``. The content of the file is as follows:

.. literalinclude:: ../sample/input.toml

The file is written in TOML format, and organized by the following sections.

``[log]`` section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section contains settings on log outputs.
``print_level`` sets the verbosity of the standard output,
and ``print_step`` sets the interval of the outputs.

``[mode]`` section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section contains settings on calculation mode and basic parameters.
We choose ``mode`` to be either the coordinate space UHF (``UHFr``) or the wave-number space UHF (``UHFk``).
The calculation parameters are specified in ``[mode.param]`` subsection.

``[file]`` section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``[file.input]`` subsection contains settings
on the directory for the input files by ``path_to_input``,
and the filename for the initial configuration by ``initial``.
If it is not specified, a random configuration will be generated.

``[file.input.interaction]`` subsection contains a list of files
associated with the geometry information and the interactions distinguished by the keywords.

``[file.output]`` subsection contains filenames
for the physical observables such as the energy by ``energy``,
for the eigenvalues and eigenvectors of the Hamiltonian by ``eigen``,
and for the one-body Green's functions by ``green``.
If they are not specified, the corresponding data will not be outputted.

See File format section for the details.


Create interaction definition files
----------------------------------------

You need to prepare data files on the geometry information of the lattice, and
the coefficients of the interactions required to construct the Hamiltonian.
The association between the types of information and the filenames is provided
in ``[file.input.interaction]`` section.

``Geometry``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The file associated with Geometry provides the geometrical information of the lattice.
An example of the file is shown below.

.. literalinclude:: ../sample/geom.dat

It contains the primitive vectors (lines 1--3), the number of orbitals (line 4),
and the Wannier centers of the orbitals (line 5 onwards).

``Transfer``, ``CoulombIntra``, ``CoulombInter``, ``Hund``, etc
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The file associated with Transfer contains the coefficients of Hamiltonian
corresponding to the transfer term of the electron systems.
The coefficients of the two-body interaction terms are stored in the files
associated with the respective keywords.
The defined types include
CoulombIntra, CoulombInter, Hund, Ising, Exchange, PairLift, and PairHop,
in accordance with the coordinate-space UHF.
These files are written in Wannier90(-like) format, as exemplified below.

.. literalinclude:: ../sample/transfer.dat

It contains
a comment (line 1),
the number of orbitals (line 2),
the number of cells ``nrpts`` of the rectangular cuboid that accommodates translation vectors (line 3),
the multiplicity factors (``nrpts`` elements with 15 points a line),
and the elements of the coefficient matrix.

Each element of the matrix consists of
translation vector :math:`r_x, r_y, r_z`,
indices of orbitals :math:`\alpha, \beta`,
and the real and imaginary part of the coefficient.

Run
----------------------------------------

Once you prepare all the input files, you can perform the calculation
by running H-wave with the input parameter file (``input.toml`` in this tutorial)
as an argument.

.. code-block:: bash

    $ hwave input.toml

The calculation starts with the logs as shown below:

.. literalinclude:: ../sample/run.log

The logs on the input files are shown, followed by the logs on the iterations of
the wave-number space UHF calculation.
The program will yield, according to the settings in ``[file.output]`` section,
the output files ``energy.dat``, ``eigen.npz``, and ``green.npz`` in ``output`` directory.

See File format section for the details of the output files.

Compile and run StdFace library
----------------------------------------------------------------

The interaction definition files can be generated easily using StdFace library.
We will provide a short instruction how to use it.

The source package of StdFace library that supports input formats of the
wave-number space UHF is available from the repository as follows.

.. code-block:: bash

    $ git clone https://github.com/issp-center-dev/StdFace.git

Then, the library is to be compiled with the commands:

.. code-block:: bash

    $ cd StdFace
    $ mkdir build && cd build
    $ cmake -DUHF=ON ..
    $ make

If the compilation is successful, you can find the executable module ``uhf_dry.out`` in ``src`` directory.

An input to ``uhf_dry.out`` can be found as ``stan.in`` in the sample directory,
which reads:

.. literalinclude:: ../sample/stan.in

- ``model`` is a keyword to choose the target model.
  Currently, only ``Hubbard`` is supported that denotes Hubbard model with the number of electrons fixed.

- ``lattice`` is a keyword to specify the lattice structure.
  In this example, the square lattice ``square`` is chosen.
  ``W`` and ``L`` denote the size of the lattice.

- ``t`` and ``V`` denote parameters of the hopping and the neighbor-site Coulomb interaction, respectively.

- ``calcmode = "uhfk""`` specifies the output to be in the Wannier90(-like) format.
  If ``exportall = 0`` is given, the outputs are compactified with zero components omitted.

See Section :ref:`Ch:HowToWannier90` for the details of input files.

Then, run ``uhf_dry.out`` with the file above as an input:

.. code-block:: bash

    $ cd path_to_Hwave/docs/tutorial/Hubbard/UHFk
    $ ln -s path_to_Stdface/build/src/uhf_dry.out .
    $ ./uhf_dry.out stan.in

When the program finishes, there are generated
a geometry information file ``geom.dat``
and interaction definition files ``transfer.dat`` and ``coulombinter.dat``,
in the current directory.

