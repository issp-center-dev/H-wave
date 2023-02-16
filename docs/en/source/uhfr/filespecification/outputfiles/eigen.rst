.. highlight:: none

.. _subsec:eigen.dat:

eigen
~~~~~~~~~~

The eigenvalues and eigenvectors of the Hamiltonian obtained by the UHF method are exported in npz (numpy zip archive) format.

The filename is chosen, with the string specified for the ``eigen`` keyword in ``file.output`` section (indicated by *eigen_str* hereafter), as ``{key}_``\ *eigen_str*\ ``.npz``, where
``{key}`` turns to be:

- ``sz-free`` if the parameter ``2Sz`` is not specified in ``mode.param`` section,

- ``spin-up`` and ``spin-down`` otherwise. (yields two files)

The code shown below is an example for reading data from the file using Python.

.. code-block:: python

    import numpy as np
    data = np.load("key_eigen_str.npz")
    eigenvalue = data["eigenvalue"]
    eigenvector = data["eigenvector"]

``eigenvalue`` holds a list of eigenvalues in ascending order.
The number of eigenvalues is ``N`` if ``2Sz`` is specified, or ``2N`` otherwise,
for ``N`` being the total number of sites.

``eigenvector`` holds the corresponding eigenvectors as a two-dimensional array: 
The first index refers to the site index ``i_site`` and the spin index ``s_spin`` (0 for up-spin, and 1 for down-spin) by:

- ``i_site`` if `2Sz` is specified, or

- ``i_site + s_spin * N`` otherwise.

The second index refers to the index of the eigenvalues.


.. raw:: latex

   \newpage
