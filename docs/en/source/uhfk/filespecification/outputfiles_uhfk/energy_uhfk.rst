.. highlight:: none

.. _subsec:energy_uhfk.dat:

energy
~~~~~~~~~~

The values of the energy, the number of electrons, and the spin
obtained by the wave-number space UHF method
are outputted.
The filename is specified by the keyword ``energy`` in ``file.output`` section in the parameter file.

An example of the file is shown as follows.

::

    Energy_Total = -5.88984624257707
    Energy_Band = -0.9265413257740396
    Energy_Coulomb = -4.963304916803031
    NCond = 8.000000000000007
    Sz = 3.2822430107160017e-07


File format
^^^^^^^^^^^^

-  ``Energy_Total = [energy_total]``

-  ``Energy_Band = [energy_band]``

-  ``Energy_``\{type} ``= [energy_``\type ``]``

-  ``NCond = [ncond]``

-  ``Sz = [sz]``

Parameters
^^^^^^^^^^

-  ``[energy_total]``

   **Type :**
   Float

   **Description :**
   The value of the energy calculated from the eigenvectors obtained by the UHF mothod.

-  ``[energy_band]``

   **Type :**
   Float

   **Description :**
   The value of the energy derived only from the eigenvalues of the Hamiltonian
   obtained by the UHF method.

-  ``[energy_``\type ``]``

   **Type :**
   Float

   **Description :**
   The value of the energy calculated separately for each interaction type.

-  ``[ncond]``

   **Type :**
   Float

   **Description :**
   The expectation value of the total number of electrons denoted by
   :math:`\sum_{i}\langle n_{i}\rangle`.

-  ``[sz]``

   **Type :**
   Float

   **Description :**
   The expectation value of the :math:`z` component of total spin :math:`S_z`
   denoted by 
   :math:`\sum_{i}\langle (n_{i\uparrow}-n_{i\downarrow})\rangle/2`.


.. raw:: latex
