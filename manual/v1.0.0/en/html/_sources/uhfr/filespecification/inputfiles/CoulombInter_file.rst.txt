.. highlight:: none

CoulombInter file
~~~~~~~~~~~~~~~~~~~~~~~~

This file determines the coefficients of the off-site Coulomb interactions
given by

.. math:: \mathcal{H} = \sum_{i,j}V_{ij} n_ {i}n_{j} .

An example of the file format is presented below.

::

    ====================== 
    NCoulombInter 6  
    ====================== 
    ========CoulombInter ====== 
    ====================== 
       0     1 -0.125000
       1     2 -0.125000
       2     3 -0.125000
       3     4 -0.125000
       4     5 -0.125000
       5     0 -0.125000

File format
^^^^^^^^^^^

-  Line 1: Header

-  Line 2: ``[ncoulombinter] [count]``

-  Lines 3-5: Header

-  Lines 6-: ``[i] [j] [val]``

Parameters
^^^^^^^^^^

-  ``[ncoulombinter]``

   **Type :**
   String (blank is not allowed)

   **Description :**
   An arbitrary keyword for the total number of the off-site Coulomb interactions.

-  ``[count]``

   **Type :**
   Integer (blank is not allowed)

   **Description :**
   An integer giving the total number of the off-site Coulomb interactions.

-  ``[i]``, ``[j]``

   **Type :**
   Integer (blank is not allowed)

   **Description :**
   An integer giving a site index (:math:`0 \le i, j < {\rm Nsite}`).

-  ``[val]``

   **Type :**
   Float (blank is not allowed)

   **Description :**
   A value for :math:`V_{ij}`.

Usage rules
^^^^^^^^^^^

-  Headers cannot be omitted.

-  The program is terminated with error if there are duplicated entries.

-  The program is terminated with error when the number of entries is different from ``[count]``.

-  The program is terminated with error if
   ``[i]`` or ``[j]``
   are outside the range of the defined values.

.. raw:: latex

   \newpage
