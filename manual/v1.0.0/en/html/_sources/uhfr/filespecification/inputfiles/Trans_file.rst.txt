.. highlight:: none

.. raw:: latex

   \newpage

.. _Subsec:Trans:

Trans file
~~~~~~~~~~~~~~~~~

This file determines the coefficients of the transfer integrals
:math:`t_{ij\sigma_1\sigma_2}`
in the Hamiltonian      

  .. math::

     \begin{aligned}
     \mathcal{H} = -\sum_{ij\sigma_1\sigma_2} t_{ij\sigma_1\sigma_2}c_{i\sigma_1}^{\dagger}c_{j\sigma_2}^{\phantom{\dagger}}.
     \end{aligned}

An example of the file format is presented below.

::

    ======================== 
    NTransfer      24  
    ======================== 
    ========i_j_s_tijs====== 
    ======================== 
        0     0     2     0   1.000000  0.000000
        2     0     0     0   1.000000  0.000000
        0     1     2     1   1.000000  0.000000
        2     1     0     1   1.000000  0.000000
        2     0     4     0   1.000000  0.000000
        4     0     2     0   1.000000  0.000000
        2     1     4     1   1.000000  0.000000
        4     1     2     1   1.000000  0.000000
        4     0     6     0   1.000000  0.000000
        6     0     4     0   1.000000  0.000000
        4     1     6     1   1.000000  0.000000
        6     1     4     1   1.000000  0.000000
        6     0     8     0   1.000000  0.000000
        8     0     6     0   1.000000  0.000000
    …

File format
^^^^^^^^^^^

-  Line 1: Header

-  Line 2: ``[ntransfer] [count]``

-  Lines 3-5: Header

-  Lines 6-: ``[i] [s1] [j] [s2]  [v.real] [v.imag]``

Parameters
^^^^^^^^^^

-  ``[ntransfer]``

   **Type :**
   String (blank is not allowed)

   **Description :**
   An arbitrary keyword for the total number of transfer integrals.

-  ``[count]``

   **Type :**
   Integer (blank is not allowed)

   **Description :**
   An integer giving the total number of transfer integrals.

-  ``[i]``, ``[j]``

   **Type :**
   Integer (blank is not allowed)

   **Description :**
   An integer giving a site index (:math:`0 \le i, j < {\rm Nsite}`).

-  ``[s1]``, ``[s2]``

   **Type :**
   Integer (blank is not allowed)

   **Description :**
   An integer giving a spin index: either 0 (up-spin) or 1 (down-spin).

-  ``[v.real]``, ``[v.imag]``

   **Type :**
   Float (blank is not allowed)

   **Description :**
   Values for the real and imaginary parts of :math:`t_{ij\sigma_1\sigma_2}`, respectively.

Usage rules
^^^^^^^^^^^

-  Headers cannot be omitted.

-  Since the Hamiltonian should be hermite, the coefficients must satisfy the relation :math:`t_{ij\sigma_1\sigma_2}=t_{ji\sigma_2\sigma_1}^{\dagger}`.
   Otherwise, the program is terminated or a warning is reported, depending on the ``strict_hermite`` parameter.

-  The program is terminated with error if there are duplicated entries.

-  The program is terminated with error when the number of entries is different from ``[count]``.

-  The program is terminated with error if
   ``[i]``, ``[j]``, ``[s1]``, or ``[s2]``
   are outside the range of the defined values.


.. raw:: latex

   \newpage
