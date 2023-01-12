.. highlight:: none

.. _subsec:energy.dat:

Energy
~~~~~~~~~~~

The values of energy, number of particles, and spin obtained by the calculations
of the UHF method are written to this file.
The filename is specified by a keyword ``energy`` in ``file.output`` section of the input parameter file.
An example of the output is presented below.

::

    Energy_total = -5.88984624257707
    Energy_band = -0.9265413257740396
    Energy_interall = -4.963304916803031
    NCond = 8.000000000000007
    Sz = 3.2822430107160017e-07

File format
^^^^^^^^^^^

-  Energy_total = ``[energy_total]``

-  Energy_band = ``[energy_band]``

-  Energy_interall = ``[energy_interall]``

-  NCond = ``[ncond]``

-  Sz = ``[sz]``


Parameters
^^^^^^^^^^

-  ``[energy_total]``

   **Type :**
   Float

   **Description :**
   The value of the total energy which is calculated using the eigenvectors obtained by the UHF method.

-  ``[energy_band]``

   **Type :**
   Float

   **Description :**
   The value of the energy which is derived from
   the eigenvalues of the Hamiltonian obtained by the UHF mothod.

-  ``[energy_interall]``

   **Type :**
   Float

   **Description :**
   The value of the energy of the interaction terms.

-  ``[ncond]``

   **Type :**
   Float

   **Description :**
   The expectation value of the number of particles
   :math:`\sum_{i}\langle n_{i}\rangle` .

-  ``[sz]``

   **Type :**
   Float

   **Description :**
   The expectation value of the :math:`z` component of the total spin
   :math:`S_z = \sum_{i}\langle (n_{i\uparrow}-n_{i\downarrow})\rangle/2` .


.. raw:: latex

   \newpage
