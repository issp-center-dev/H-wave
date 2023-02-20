.. highlight:: none

.. _subsec:chiq_rpa:

chiq and chi0q
~~~~~~~~~~~~~~~~~~~~

The susceptibility matrix and the irreducible susceptibility matrix are exported
in NumPy zip (npz) format.
Using the string (referred to as *chiq_str*) specified by the keyword ``chiq`` or ``chi0q``
in ``file.output`` section in the parameter file,
the filename is chosen as *chiq_str*\ ``.npz``.

The following code is an example for reading the data from the output file.

.. code-block:: python

    import numpy as np
    data = np.load("chiq_str.npz")
    chiq = data["chiq"]
    freq_index = data["freq_index"]


``chiq`` contains the susceptibility matrix :math:`\chi^{\alpha\alpha^\prime\beta\beta^\prime}(\vec{q})` for each wave number and Matsubara frequency.


.. raw:: latex
