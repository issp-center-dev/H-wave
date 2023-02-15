.. highlight:: none

.. _Subsec:cgcisajs:

green
~~~~~~~~~~

The one-body Green's function
:math:`\langle c_{i\sigma_1}^{\dagger} c_{j\sigma_2}^{\phantom{\dagger}} \rangle`
obtained by the UHF calculation is to be exported.
The indices of the component to be written are specified by ``OneBodyG``.
The filename is associated with a keyword ``green`` in ``file.output`` section of the input parameter file.
An example of the output is presented below.

::

    0 0 0 0  0.9517526553947047  0.0
    0 0 1 0 -0.03971951040016314 0.0
    0 0 2 0  0.09202884754223833 0.0
    0 0 3 0 -0.039719448981075135 0.0
    0 0 4 0  0.09202884754219534 0.0
    0 0 5 0 -0.03971947216145664 0.0
    0 0 6 0  0.09202884753253462 0.0
    0 0 7 0  0.09202884754259735 0.0
    0 1 0 1  0.04824734460529617 0.0
    0 1 1 1  0.03971951040016307 0.0
    …

File format
^^^^^^^^^^^

-  ``[i] [s1] [j] [s2]  [v.real] [v.imag]``


Parameters
^^^^^^^^^^

-  ``[i]``, ``[j]``

   **Type :**
   Integer

   **Description :**
   An integer giving a site index (:math:`0 \le i, j < {\rm Nsite}`).

-  ``[s1]``, ``[s2]``

   **Type :**
   Integer

   **Description :**
   An integer giving a spin index: either 0 (up-spin) or 1 (down-spin).

-  ``[v.real]``, ``[v.imag]``

   **Type :**
   Float

   **Description :**
   Values for the real and imaginary parts of
   :math:`\langle c_{i\sigma_1}^{\dagger} c_{j\sigma_2}^{\phantom{\dagger}} \rangle`,
   respectively.


.. raw:: latex

   \newpage
