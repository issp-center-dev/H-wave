.. highlight:: none

.. _Subsec:onebodyg:

OneBodyG file
~~~~~~~~~~~~~~~~~~~~

This file determines the list of indices of the one-body Green's function
:math:`\langle c_{i\sigma_1}^{\dagger} c_{j\sigma_2}^{\phantom{\dagger}} \rangle`
to be exported.
An example of the file format is presented below.

::

    ===============================
    NCisAjs         24
    ===============================
    ======== Green functions ======
    ===============================
        0     0     0     0
        0     1     0     1
        1     0     1     0
        1     1     1     1
        2     0     2     0
        2     1     2     1
        3     0     3     0
        3     1     3     1
        4     0     4     0
        4     1     4     1
        5     0     5     0
        5     1     5     1
        6     0     6     0
        6     1     6     1
        7     0     7     0
        7     1     7     1
        8     0     8     0
        8     1     8     1
        9     0     9     0
        9     1     9     1
       10     0    10     0
       10     1    10     1
       11     0    11     0
       11     1    11     1

File format
^^^^^^^^^^^

-  Line 1: Header

-  Line 2: ``[ncisajs] [count]``

-  Lines 3-5: Header

-  Lines 6-: ``[i]  [s1]  [j]  [s2]``

Parameters
^^^^^^^^^^

-  ``[ncisajs]``

   **Type :**
   String (blank is not allowed)

   **Description :**
   An arbitrary keyword for the total number of components of the one-body Green's functions.

-  ``[count]``

   **Type :**
   Integer (blank is not allowed)

   **Description :**
   An integer giving the total number of components of the one-body Green's functions.

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

Usage rules
^^^^^^^^^^^

-  Headers cannot be omitted.

-  The program is terminated with error if there are duplicated entries.

-  The program is terminated with error when the number of entries is different from ``[count]``.

-  The program is terminated with error if
   ``[i]``, ``[j]``, ``[s1]``, or ``[s2]``
   are outside the range of the defined values.


.. raw:: latex

   \newpage
