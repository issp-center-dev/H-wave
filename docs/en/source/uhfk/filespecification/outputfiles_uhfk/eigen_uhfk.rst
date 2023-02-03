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

``eigenvalue`` contains the eigenvalues :math:`\lambda_l(\vec{k})` for each wave number.
The wave number is taken in unit of sublattice when the sublattice is considered.
The data format is a numpy ndarray with the layout as ``eigenvalue[k][l]``, where
``k`` refers to the linearlized index of the wave number vector :math:`\vec{k}`,
and ``l`` refers to the index of eigenvalues in a cell.
The indices of the wave number vector ``(kx, ky, kz)`` is packed into the linearlized index
``k`` by ``k = kz + Nz * (ky + Ny * kx)``.


``eigenvector`` contains the corresponding eigenvectors.
The data format is a numpy ndarray with the layout as ``eigenvector[k][l][j]``, where
``k`` and ``l`` refer to the indices of the corresponding wave number and eigenvalue,
and ``j`` refers to the index of the orbital and spin in a cell.

.. raw:: latex
