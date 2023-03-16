==================
Tutorial
==================

To use H-wave for the calculation of Random Phase Approximation (RPA),
you need to prepare input files:

#. an input parameter file,
#. interaction definition files,

before running the program.

In the following, we provide a tutorial based on a sample in
``docs/tutorial/Hubbard/RPA`` directory.
The interaction definition files can be generated using StdFace library.
See :ref:`Ch:StdFace` section for the details.

Create a parameter file
--------------------------------

We prepare an input parameter file that contains basic parameters as well as settings on
inputs and outputs.
A sample file can be found in ``docs/tutorial/Hubbard/RPA`` directory by a filename
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
We choose ``mode`` to be ``"RPA"``.
The calculation parameters are specified in ``[mode.param]`` subsection.

``[file]`` section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``[file.input]`` subsection contains settings
on the directory for the input files by ``path_to_input``,
and the filename for the initial configuration by ``initial``.
If the latter is not specified, a random configuration will be generated.

``[file.input.interaction]`` subsection contains a list of files
associated with the geometry information and the interactions distinguished by the keywords.

``[file.output]`` subsection contains filenames
for the physical observables such as the energy by ``energy``,
for the eigenvalues and eigenvectors of the Hamiltonian by ``eigen``,
and for the one-body Green's functions by ``green``.
If they are not specified, the corresponding data will not be outputted.

See :ref:`Ch:Config_rpa` section for the details.


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
in accordance with the UHF calculation.
These files are written in Wannier90(-like) format, as exemplified below.

.. literalinclude:: ../sample/transfer.dat

It contains
a comment (line 1),
the number of orbitals (line 2),
the number of cells ``nrpts`` of the rectangular cuboid that accommodates translation vectors (line 3),
the multiplicity factors (``nrpts`` elements, with 15 elements per line),
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

.. literalinclude:: ../sample/runlog.dat

The logs on the input files are shown, followed by the logs on the process of RPA calculations.
The program will yield, according to the settings in ``[file.output]`` section,
the output files ``chi0q.npz`` and ``chiq.npz`` in ``output`` directory.

See :ref:`Sec:outputfile_rpa` section for the details of the output files.

A tool is prepared in ``sample/RPA/view.py`` for visualizing the calculation results
as a post-process.
Let us copy the script file to the current directory, and run the script as follows:

.. code-block:: bash

    $ python3 view.py

The script reads ``output/chi0q.npz`` and ``output/chiq.npz``, and writes the values of
the charge susceptibility :math:`\chi_c(\vec{q})`
and the spin susceptibility :math:`\chi_s(\vec{q})`
at Matsubara frequency :math:`{\rm i}\omega_m=0`
for each :math:`\vec{q}`
to the standard output.
It also produces the figures of these quantities in PNG format shown as below:

.. tabularcolumns:: CC

.. raw:: latex

	 \begingroup
	 \renewcommand{\hline}{}

.. list-table::

  * - .. figure:: ../sample/chic.png

         :math:`\chi_c(\vec{q})`

    - .. figure:: ../sample/chis.png

         :math:`\chi_s(\vec{q})`

.. raw:: latex

	 \endgroup



