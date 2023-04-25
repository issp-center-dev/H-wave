.. highlight:: none

.. _Subsec:interall:

InterAll file
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This file determines the coefficients of the generalized two-body interaction integrals
:math:`I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}`
of the Hamiltonian

.. math::

   \mathcal{H} = \sum_{i,j,k,l}\sum_{\sigma_1,\sigma_2, \sigma_3, \sigma_4}
   I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}c_{i\sigma_1}^{\dagger}c_{j\sigma_2}^{\phantom{\dagger}}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}^{\phantom{\dagger}}.

An example of the file format is presented below.

::

    ====================== 
    NInterAll      36  
    ====================== 
    ========zInterAll===== 
    ====================== 
    0    0    0    1    1    1    1    0   0.50  0.0
    0    1    0    0    1    0    1    1   0.50  0.0
    0    0    0    0    1    0    1    0   0.25  0.0
    0    0    0    0    1    1    1    1  -0.25  0.0
    0    1    0    1    1    0    1    0  -0.25  0.0
    0    1    0    1    1    1    1    1   0.25  0.0
    2    0    2    1    3    1    3    0   0.50  0.0
    2    1    2    0    3    0    3    1   0.50  0.0
    2    0    2    0    3    0    3    0   0.25  0.0
    2    0    2    0    3    1    3    1  -0.25  0.0
    2    1    2    1    3    0    3    0  -0.25  0.0
    2    1    2    1    3    1    3    1   0.25  0.0
    4    0    4    1    5    1    5    0   0.50  0.0
    4    1    4    0    5    0    5    1   0.50  0.0
    4    0    4    0    5    0    5    0   0.25  0.0
    4    0    4    0    5    1    5    1  -0.25  0.0
    4    1    4    1    5    0    5    0  -0.25  0.0
    4    1    4    1    5    1    5    1   0.25  0.0
    …

File format
^^^^^^^^^^^

-  Line 1: Header

-  Line 2: ``[ninterall] [count]``

-  Lines 3-5: Header

-  Lines 6-: ``[i] [s1] [j] [s2] [k] [s3] [l] [s4] [v.real] [v.imag]``

Parameters
^^^^^^^^^^

-  ``[ninterall]``

   **Type :**
   String (blank is not allowed)

   **Description :**
   An arbitrary keyword for the total number of the two-body interactions.

-  ``[count]``

   **Type :**
   Integer (blank is not allowed)

   **Description :**
   An integer giving the total number of the two-body interactions.

-  ``[i]``, ``[j]``, ``[k]``, ``[l]``

   **Type :**
   Integer (blank is not allowed)

   **Description :**
   An integer giving a site index (:math:`0 \le i, j, k, l < {\rm Nsite}`).

-  ``[s1]``, ``[s2]``, ``[s3]``, ``[s4]``

   **Type :**
   Integer (blank is not allowed)

   **Description :**
   An integer giving a spin index: either 0 (up-spin) or 1 (down-spin).


-  ``[v.real]``, ``[v.imag]``

   **Type :**
   Float (blank is not allowed)

   **Description :**
   Values for the real and imaginary parts of :math:`I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}`, respectively.

Usage ruels
^^^^^^^^^^^

-  Headers cannot be omitted.

-  Since the Hamiltonian should be hermite, the coefficients must satisfy the relation
   :math:`I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}=I_{lkji\sigma_4\sigma_3\sigma_2\sigma_1}^{\dagger}`.
   Otherwise, the program is terminated or a warning is reported, depending on the ``strict_hermite`` parameter.
   It is noted that the term of the Hermite conjugate for 
   :math:`I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4} c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}`
   should be read
   :math:`I_{lkji\sigma_4\sigma_3\sigma_2\sigma_1} c_{l\sigma_4}^{\dagger}c_{k\sigma_3}c_{j\sigma_2}^{\dagger}c_{i\sigma_1}`.

-  The program is terminated with error if there are duplicated entries.

-  The program is terminated with error when the number of entries is different from ``[count]``.

-  The program is terminated with error if
   ``[i]``, ``[j]``, ``[k]``, ``[l]``, ``[s1]``, ``[s2]``, ``[s3]``, or ``[s4]``
   are outside the range of the defined values.

.. raw:: latex

   \newpage
