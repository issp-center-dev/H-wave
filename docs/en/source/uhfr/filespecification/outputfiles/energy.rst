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
    FreeEnergy_total = -5.91234017896271
    FreeEnergy_band = -0.9490352621596976
    NCond = 8.000000000000007
    Sz = 3.2822430107160017e-07
    ChemicalPotential = 0.123456789
    Sx = 0.0
    Sy = 0.0

.. note::

   **Changed in this version.** UHFr and normal-mode UHFk now append ``Sx`` and ``Sy`` as the
   last two lines of the file; all other lines are unchanged. Files written
   by older versions lack them. A script that reads the last line (for
   example ``tail -n 1`` to obtain ``ChemicalPotential``) now gets ``Sy``;
   read the value by its key instead. The mean-field state can point along
   any spin axis, so Sx and Sy can be nonzero also without spin-orbit
   coupling.

File format
^^^^^^^^^^^

-  Energy_total = ``[energy_total]``

-  Energy_band = ``[energy_band]``

-  Energy_interall = ``[energy_interall]``

-  FreeEnergy_total = ``[free_energy_total]``

-  FreeEnergy_band = ``[free_energy_band]``

-  NCond = ``[ncond]``

-  Sz = ``[sz]``

-  ChemicalPotential = ``[mu]``

-  Sx = ``[sx]``

-  Sy = ``[sy]``


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
   The expectation value of the number of particles,
   :math:`\sum_{i}\langle n_{i}\rangle` .

-  ``[sz]``

   **Type :**
   Float

   **Description :**
   The expectation value of the :math:`z` component of the total spin,
   :math:`S_z = \sum_{i}\langle (n_{i\uparrow}-n_{i\downarrow})\rangle/2` .

-  ``[mu]``

   **Type :**
   Float

   **Description :**
   The chemical potential (Fermi level). At :math:`T=0` it is taken as the
   midpoint of the HOMO-LUMO gap (the :math:`T\to 0^+` limit of the
   finite-temperature :math:`\mu`). Note that for an insulator the particle
   number is independent of :math:`\mu` inside the gap, so the
   finite-temperature :math:`\mu` is not uniquely fixed within the gap.
   When several blocks (e.g. spin up/down) exist, the value is written per
   block as ``ChemicalPotential_{block}``.

-  ``[sx]``, ``[sy]``

   **Type :**
   Float

   **Description :**
   The expectation values of the :math:`x` and :math:`y` components of the
   total spin, in units of :math:`\hbar` and summed over all sites like
   ``[sz]``:
   :math:`S_x = \sum_{i}\mathrm{Re}\,(\langle c^{\dagger}_{i\uparrow}c_{i\downarrow}\rangle + \langle c^{\dagger}_{i\downarrow}c_{i\uparrow}\rangle)/2`
   and
   :math:`S_y = \sum_{i}\mathrm{Re}\,[-i(\langle c^{\dagger}_{i\uparrow}c_{i\downarrow}\rangle - \langle c^{\dagger}_{i\downarrow}c_{i\uparrow}\rangle)]/2`.
   The mean-field state can point along any spin axis, so Sx and Sy can be
   nonzero also without spin-orbit coupling.
   With a fixed ``2Sz`` the two spin blocks are separate, and Sx and Sy
   are zero.


.. raw:: latex

   \newpage
