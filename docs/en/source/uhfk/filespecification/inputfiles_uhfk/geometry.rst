.. highlight:: none

Geometry input file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Geometry input file describes the geometry information of the lattice.
An example of the file is shown as follows.

::

   3.7599302871   0.0000000000   0.0000000000
   0.0000000000   3.7599302871   0.0000000000
   0.0000000000   0.0000000000   5.4822004186
        10
   -7.179835091886330E-003 -3.812050198019962E-002  1.639284152926924E-003
    1.346463812592166E-002  6.709778405878775E-003 -6.812442303544219E-003
    0.495705070884200      -0.457955704941170      -4.077818544354700E-003
   -1.577970702078565E-004 -2.999005205319096E-004 -1.190284144276225E-004
   -1.302397074478660E-003 -5.021621895411691E-003 -3.514564279609852E-004
    0.504124376959700       0.457760356450585      -2.634809811615298E-003
    0.499384075989520      -0.494227365093439       6.927730957590197E-003
   -5.164444920392309E-003  3.667887236852975E-002  4.972296517752579E-003
    0.500170586121734       0.499747448247510       2.760670734661295E-003
    0.500734036298328       0.494793997305026      -2.212377045150314E-003


File format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  Lines 1-3: ``[ax_i] [ay_i] [az_i]``

-  Line 4: ``[Norbit]``

-  Lines 5-: ``[vx_i] [vy_i] [vz_i]``

Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  ``[ax_i]``, ``[ay_i]``, ``[az_i]``

   **Type :**
   Float

   **Description :**
   These parameters for ``i`` from 1 to 3 specify the primitive vectors :math:`\vec{a}_1, \vec{a}_2, \vec{a}_3`.

-  ``[Norbit]``

   **Type :**
   Integer

   **Description :**
   This parameter specifies the number of orbitals :math:`N_\text{orbit}` in a unitcell.

-  ``[vx_i]``, ``[vy_i]``, ``[vz_i]``

   **Type :**
   Float

   **Description :**
   These parameters specify the Wannier center :math:`\vec{v}_i` of each orbital
   in the fractional coordinates.


Usage rules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- The indices of the orbitals are implicitly assigned from 1 to :math:`N_\text{orbit}`
  in the order of the Wannier centers.

.. raw:: latex
