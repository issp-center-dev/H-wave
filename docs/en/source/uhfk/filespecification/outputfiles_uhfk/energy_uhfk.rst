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
    FreeEnergy_Total = -5.91234017896271
    FreeEnergy_Band = -0.9490352621596976
    NCond = 8.000000000000007
    Sz = 3.2822430107160017e-07
    ChemicalPotential = 0.123456789


File format
^^^^^^^^^^^^

-  ``Energy_Total = [energy_total]``

-  ``Energy_Band = [energy_band]``

-  ``Energy_``\{type} ``= [energy_``\type ``]``

-  ``FreeEnergy_Total = [free_energy_total]``

-  ``FreeEnergy_Band = [free_energy_band]``

-  ``NCond = [ncond]``

-  ``Sz = [sz]``

-  ``ChemicalPotential = [mu]``

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

-  ``[free_energy_total]``

   **Type :**
   Float

   **Description :**
   The Helmholtz free energy :math:`F = E - TS`, i.e. the quantity minimized
   by the finite-temperature self-consistent loop. At :math:`T=0` it coincides
   with the internal energy ``[energy_total]``.

-  ``[free_energy_band]``

   **Type :**
   Float

   **Description :**
   The band contribution to the free energy,
   :math:`\mu N + \Omega_0` with
   :math:`\Omega_0 = -T\sum_n \ln(1+e^{-(\varepsilon_n-\mu)/T})`.
   At :math:`T=0` it coincides with ``[energy_band]``.

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

-  ``[mu]``

   **Type :**
   Float

   **Description :**
   The chemical potential (Fermi level). At :math:`T=0` it is taken as the
   midpoint of the HOMO-LUMO gap (the :math:`T\to 0^+` limit of the
   finite-temperature :math:`\mu`). Note that for an insulator the particle
   number is independent of :math:`\mu` inside the gap, so the
   finite-temperature :math:`\mu` is not uniquely fixed within the gap.
   When several :math:`\mu`-groups exist, the value is written per group as
   ``ChemicalPotential_{g}``.


.. raw:: latex
