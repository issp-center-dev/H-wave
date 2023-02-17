.. highlight:: none

Hund file
~~~~~~~~~

This file determines the coefficients of the Hund couplings given by

.. math:: \mathcal{H} = -\sum_{i,j}J_{ij}^{\rm Hund} (n_{i\uparrow}n_{j\uparrow}+n_{i\downarrow}n_{j\downarrow}) .

An example of the file format is presented below.

::

    ====================== 
    NHund 6  
    ====================== 
    ========Hund ====== 
    ====================== 
       0     1 -0.250000
       1     2 -0.250000
       2     3 -0.250000
       3     4 -0.250000
       4     5 -0.250000
       5     0 -0.250000

File format
^^^^^^^^^^^

-  Line 1: Header

-  Line 2: ``[nhund] [count]``

-  Lines 3-5: Header

-  Lines 6-: ``[i] [j] [val]``

Parameters
^^^^^^^^^^

-  ``[nhund]``

   **Type :**
   String (blank is not allowed)

   **Description :**
   An arbitrary keyword for the total number of the Hund couplings.

-  ``[count]``

   **Type :**
   Integer (blank is not allowed)

   **Description :**
   An integer giving the total number of the Hund couplings.

-  ``[i]``, ``[j]``

   **Type :**
   Integer (blank is not allowed)

   **Description :**
   An integer giving a site index (:math:`0 \le i, j < {\rm Nsite}`).

-  ``[val]``

   **Type :**
   Float (blank is not allowed)

   **Description :**
   A value for :math:`J_{ij}^{\rm Hund}`.

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
