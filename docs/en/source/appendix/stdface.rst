.. _Ch:StdFace:

Generation of interaction files using StdFace library
================================================================

Compile StdFace library
----------------------------------------------------------------

The interaction definition files can be generated easily using StdFace library.
We will provide a short instruction how to use it.

The source package of StdFace library that supports input formats of the Hwave is available from the repository as follows.

.. code-block:: bash

    $ git clone https://github.com/issp-center-dev/StdFace.git

Then, the library is to be compiled with the commands:

.. code-block:: bash

    $ cd StdFace
    $ mkdir build && cd build
    $ cmake -DHWAVE=ON ..
    $ make

If the compilation is successful, you can find the executable module ``hwave_dry.out`` in ``src`` directory.

An input to ``hwave_dry.out`` can be found as ``stan.in`` in the sample directory,
which reads:

.. literalinclude:: ../rpa/sample/stan.in

- ``model`` is a keyword to choose the target model.
  Currently, only ``Hubbard`` is supported that denotes Hubbard model with the number of electrons fixed.

- ``lattice`` is a keyword to specify the lattice structure.
  In this example, the square lattice ``square`` is chosen.
  ``W`` and ``L`` denote the size of the lattice.

- ``t`` and ``V`` denote parameters of the hopping and the neighbor-site Coulomb interaction, respectively.

- ``calcmode = "uhfk"`` and ``calcmode = "rpa"`` specify the output to be in the Wannier90(-like) format.
  ``calcmode = "uhfr"`` specifies the output for input files of UHFr. The default is ``calcmode = "uhfk"``.
  If ``exportall = 0`` is given, the outputs are compactified with zero components omitted.

See Section :ref:`Ch:HowToExpert` , :ref:`Ch:HowToWannier90` , :ref:`Ch:HowToWannier90_rpa` for the details of input files.


Run StdFace library
----------------------------------------------------------------

Then, run ``hwave_dry.out`` with the file above as an input:

.. code-block:: bash

    $ cd path_to_Hwave/docs/tutorial/Hubbard/RPA
    $ ln -s path_to_Stdface/build/src/hwave_dry.out .
    $ ./uhf_dry.out stan.in

When the program finishes, a geometry information file ``geom.dat``
and interaction definition files ``transfer.dat`` and ``coulombinter.dat``,
are generated in the current directory.

