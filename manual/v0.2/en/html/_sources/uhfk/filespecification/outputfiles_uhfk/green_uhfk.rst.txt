.. highlight:: none

.. _Subsec:green_uhfk:

green
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The one-body Green's function
:math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}^{\phantom{\dagger}}\rangle`
calculated by the wave-number space UHF method is exported in NumPy zip (npz) format. 
Using the string (referred to as *green_str*) specified by the keyword ``green``
in ``file.output`` section in the parameter file,
the filename is chosen as *green_str*\ ``.npz``.

The data is bound to the key ``green``. The data format is a numpy ndarray
with the layout ``ndarray(r, s, a, t, b)``, where

- ``r`` denotes a linearlinzed index of translation vector :math:`[r_x\ r_y\ r_z]`,
  where the indices are packed into ``r`` by ``r`` :math:`= r_z + N_z \cdot (r_y + N_y r_x)`.
- ``a``, ``b`` denote the indices of the orbitals :math:`\alpha, \beta`,
- ``s``, ``t`` denote the indices of the spins :math:`\sigma_1, \sigma_2`.

The output can be used as an initial configuration of the Green's function
specified by the keyword ``initial`` in ``file.input`` section.

When the sublattice is considerd,
the Green's function in unit of the sublattice is also stored with the key ``green_sublattice``.
The indices of the data are regarded as those of the sublattice.

The following code is an example for reading the data from the output file.

.. code-block:: python

    import numpy as np
    data = np.load("green.dat.npz")
    green = data["green"]

.. raw:: latex
