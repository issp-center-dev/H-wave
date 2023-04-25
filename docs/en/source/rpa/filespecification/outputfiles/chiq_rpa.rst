.. highlight:: none

.. _subsec:chiq_rpa:

chiq and chi0q
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The susceptibility matrix and the irreducible susceptibility matrix are exported
in NumPy zip (npz) format.
Using the string (referred to as *chiq_str*) specified by the keyword ``chiq`` or ``chi0q``
in ``file.output`` section in the parameter file,
the filename is chosen as *chiq_str*\ ``.npz``.

The file contains several arrays bound to the following keys:

- ``chiq`` or ``chi0q``:

  The susceptibility matrix or the irreducible susceptibility matrix. Their data layout is described in the following sections.

- ``freq_index``:

  The value or the range of Matsubara frequency is specified by ``matsubara_frequency`` parameter. The array bound to ``freq_index`` relates the index of the output data and the label of the actual Matsubara frequency.

- ``wavevector_unit`` and ``wavevector_index``:

  These arrays refer to the information of the wave number vectors. See :ref:`Output files of UHFk <Subsec:eigen_uhfk.dat>` for details.

When the sublattice is considered, the indices of the wave numbers and the orbitals are
regarded as those of the sublattice.

The output file of ``chi0q`` can be used as a pre-calculated input of the irreducible
susceptibility by specifying the file to ``chi0q_init`` in ``file.input`` section.


Data format of chi0q
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Data format of ``chi0q`` relies on the presence of spin-orbital interaction and external field, and the value of ``mode.calc_scheme`` parameter, and takes one of the following:

- "spin-free" case:

  If ``enable_spin_orbital`` parameter is set to ``false``, or even if it is set to ``true`` when :math:`T_{\tilde\alpha\tilde\beta}(k)` is diagonal and symmetric with respect to spin degree of freedom, while the external field is not present, the spin-independent irreducible susceptibility matrix is exported. 
  
  - When ``calc_scheme = general``, the array format takes the form of ``ndarray(l,q,a,ap,b,bp)`` whose indices are given as follows:

    - ``l``: label of Matsubara frequency. The map from the label to the index is provided by the aforementioned array ``freq_index``.

    - ``q``: linearlized index of wave-number indices :math:`[ q_x\ q_y\ q_z ]`, where :math:`q = q_z + N_z\cdot(q_y + N_y\cdot q_x)`.

    - ``a``, ``ap``, ``b``, ``bp``: indices of the orbitals not including spin degree of freedom. They correspond to :math:`\alpha`, :math:`\alpha^\prime`, :math:`\beta`, :math:`\beta^\prime`. 

  - When ``calc_scheme = reduced`` or ``squashed``, the array format takes the form of ``ndarray(l,q,a,b)`` whose indices are same as the above.

- "spin-diagonal" case:

  If ``enable_spin_orbital`` parameter is set to ``false`` and the external field is present, or it is set to ``true`` while :math:`T_{\tilde\alpha\tilde\beta}(k)` is diagonal with respect to spin degree of freedom, the spin-up and spin-down components of the irreducible susceptibility matrix are exported.

  - When ``calc_scheme = general``, the array format takes the form of ``ndarray(s,l,q,a,ap,b,bp)``, where ``s = 0`` denotes spin-up component and ``s = 1`` does spin-down component. The other indices are same as the above.

  - When ``calc_scheme = reduced`` or ``squashed``, the array format takes the form of ``ndarray(s,l,q,a,b)``. The indices are same as above.

- "spinful" case:

  If ``enable_spin_orbital`` parameter is set to ``true``, and :math:`T_{\tilde\alpha\tilde\beta}(k)` takes a general form, the irreducible susceptibility matrix with the generalized orbital indices is exported.

  - When ``calc_scheme = general``, the array format takes the form of ``ndarray(l,q,a,ap,b,bp)``, where ``a``, ``ap``, ``b``, and ``bp`` corresponding to the generalized orbital indices including spin degree of freedom denoted by :math:`\tilde\alpha`, :math:`\tilde\alpha^\prime`, :math:`\tilde\beta`, and :math:`\tilde\beta^\prime`, respectively.

  - When ``calc_scheme = reduced``, the array format takes the form of ``ndarray(l,q,a,b)``, where ``a`` and ``b`` corresponding to the generalized orbital indices :math:`\tilde\alpha` and :math:`\tilde\beta`, respectively.


Data format of chiq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Data format of ``chiq`` takes the following form depending on the value of ``calc_scheme`` parameter:

- When ``calc_scheme = general``, the array format takes the form of ``ndarray(l,q,a,ap,b,bp)``, where ``a``, ``ap``, ``b``, and ``bp`` correspond to the generalized orbital indices including spin degree of freedom denoted by :math:`\tilde\alpha`, :math:`\tilde\alpha^\prime`, :math:`\tilde\beta`, and :math:`\tilde\beta^\prime`, respectively.

- When ``calc_scheme = reduced``, the array format takes the form of ``ndarray(l,q,a,b)``, where ``a`` and ``b`` correspond to the generalized orbital indices :math:`\tilde\alpha` and :math:`\tilde\beta`, respectively.

- When ``calc_scheme = squashed``, the array format takes the form of ``ndarray(l,q,s1,s2,a,s3,s4,b)``, where ``a`` and ``b`` correspond to the orbital indices :math:`\alpha` and :math:`\beta`, respectively, and ``s1``, ``s2``, ``s3``, ``s4`` denote spin indices :math:`\sigma`, :math:`\sigma^\prime`, :math:`\sigma_1`, :math:`\sigma_1^\prime`, respectively. See :ref:`Algorithm<Algorithm_sec>` section for the notation.


Example for reading data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following code is an example for reading the data from the output file.

.. code-block:: python

    import numpy as np
    data = np.load("chiq_str.npz")

    chiq = data["chiq"]
    freq_index = data["freq_index"]


.. raw:: latex
