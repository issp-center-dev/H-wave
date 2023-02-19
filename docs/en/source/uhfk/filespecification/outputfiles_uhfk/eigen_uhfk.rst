.. highlight:: none

.. _subsec:eigen_uhfk.dat:

eigen
~~~~~~~~~~

The eigenvalues and eigenvectors of the Hamiltonian at the convergence
are exported in NumPy zip (npz) format.
Using the string (referred to as *eigen_str*) specified by the keyword ``eigen``
in ``file.output`` section in the parameter file,
the filename is chosen as *eigen_str*\ ``.npz``.

The following code is an example for reading the data from the output file.

.. code-block:: python

    import numpy as np
    data = np.load("eigen_str.npz")
    eigenvalue = data["eigenvalue"]
    eigenvector = data["eigenvector"]

    wavevector_unit = data["wavevector_unit"]
    wavevector_index = data["wavevector_index"]

``eigenvalue`` contains the eigenvalues :math:`\lambda_l(\vec{k})` for each wave number.
The wave number is taken in unit of sublattice when the sublattice is considered.
The data format is a numpy ndarray with the layout as ``eigenvalue[k][l]``, where
``k`` refers to the linearlized index of the wave number vector :math:`\vec{k}` (see below),
and ``l`` refers to the index of eigenvalue.
When ``Sz`` is fixed, ``l`` is given by ``l = l' + Norb * s`` where ``l'`` is the
index of the eigenvalue in a cell, and ``s`` refers to the spin index
(0 for up-spin, and 1 for down-spin).

``eigenvector`` contains the corresponding eigenvectors.
The data format is a numpy ndarray with the layout as ``eigenvector[k][j][l]``, where
``k`` and ``l`` refer to the indices of the corresponding wave number and eigenvalue,
and ``j`` refers to the index of the orbital and spin in a cell.

``wavevector_unit`` and ``wavevector_index`` refer to the information of the wave number vectors.
``wavevector_unit`` contains the unit wave number vectors given by
:math:`2\pi\vec{b}_i/N_i` with :math:`\vec{b}_i` being reciprocal lattice vectors.
``wavevector_index`` contains the map from the index ``k``
to the indices of the wave number vector ``(kx, ky, kz)``.
The wave number vector that corresponds to the index ``k`` is obtained by

.. code-block:: python

    k_vec = np.dot(wavevector_index[k], wavevector_unit)


.. raw:: latex
