.. highlight:: none

.. _subsec:chiq_rpa:

chiq and chi0q
~~~~~~~~~~~~~~~~~~~~

The susceptibility matrix and the irreducible susceptibility matrix are exported
in NumPy zip (npz) format.
Using the string (referred to as *chiq_str*) specified by the keyword ``chiq`` or ``chi0q``
in ``file.output`` section in the parameter file,
the filename is chosen as *chiq_str*\ ``.npz``.

The data are bound to keys ``chiq`` and ``chi0q``, respectively.
They are represented by the array ``ndarray(l,q,a,ap,b,bp)``, with the following indices:

  - ``l``: labels of Matsubara frequency. The map from the label to the index is provided by
    the array associated to the key ``freq_index``.

  - ``q``: linearlized index of wave-number indices :math:`[ q_x\ q_y\ q_z ]`, where
    :math:`q = q_z + N_z\cdot(q_y + N_y\cdot q_x)`.

  - ``a``, ``ap``, ``b``, ``bp``: indices of the generalized orbitals corresponding to
    :math:`\tilde\alpha`, :math:`\tilde\alpha^\prime`,
    :math:`\tilde\beta`, and :math:`\tilde\beta^\prime`, respectively.
    A generalized orbital index :math:`\tilde\alpha` refers to the orbital :math:`\alpha`
    and spin :math:`\sigma` by :math:`\tilde\alpha = \alpha + N_\text{orb}\cdot\sigma`,
    where :math:`N_\text{orb}` denotes the number of orbitals.

The value or the range of Matsubara frequency to be exported is specified by
``matsubara_frequency`` parameter.

The output file of ``chi0q`` can be used as a pre-calculated input of the irreducible
susceptibility by specifying the file to ``chi0q_init`` in ``file.input`` section.

When the sublattice is considered, the indices of the wave numbers and the orbitals are
regarded as those of the sublattice.

The following code is an example for reading the data from the output file.

.. code-block:: python

    import numpy as np
    data = np.load("chiq_str.npz")
    chiq = data["chiq"]
    freq_index = data["freq_index"]


.. raw:: latex
