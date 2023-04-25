.. highlight:: none

CoulombIntra file
~~~~~~~~~~~~~~~~~~~~~~~~

This file determines the coefficients of the on-site Coulomb interactions
given by

.. math:: \mathcal{H} = \sum_{i}U_i n_{i \uparrow}n_{i \downarrow} .

An example of the file format is presented below.

::

    ====================== 
    NCoulombIntra 6  
    ====================== 
    ==== CoulombIntra ====
    ====================== 
       0  4.000000
       1  4.000000
       2  4.000000
       3  4.000000
       4  4.000000
       5  4.000000

File format
^^^^^^^^^^^

-  Line 1: Header

-  Line 2: ``[ncoulombintra] [count]``

-  Lines 3-5: Header

-  Lines 6-: ``[i] [val]``

Parameters
^^^^^^^^^^

-  ``[ncoulombintra]``

   **Type :**
   String (blank is not allowed)

   **Description :**
   An arbitrary keyword for the total number of the on-site Coulomb interactions.

-  ``[count]``

   **Type :**
   Integer (blank is not allowed)

   **Description :**
   An integer giving the total number of the on-site Coulomb interactions.

-  ``[i]``

   **Type :**
   Integer (blank is not allowed)

   **Description :**
   An integer giving a site index (:math:`0 \le i < {\rm Nsite}`).

-  ``[val]``

   **Type :**
   Float (blank is not allowed)

   **Description :**
   A value for :math:`U_i`.

Usage rules
^^^^^^^^^^^

-  Headers cannot be omitted.

-  The program is terminated with error if there are duplicated entries.

-  The program is terminated with error when the number of entries is different from ``[count]``.

-  The program is terminated with error if ``[i]`` is outside the range of the defined value.


.. raw:: latex

   \newpage
