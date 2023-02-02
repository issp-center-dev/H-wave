.. highlight:: none

PairHop file
~~~~~~~~~~~~~~~~~~~

This file determines the coefficients of the pair-hopping couplings
given by

.. math::
   \mathcal{H} = \sum_{i,j}J_{ij}^{\rm Pair}
   (c_ {i \uparrow}^{\dagger} c_{j\uparrow}^{\phantom{\dagger}} c_{i \downarrow}^{\dagger} c_{j  \downarrow}^{\phantom{\dagger}} + \textit{h.c.}) .

An example of the file format is presented below.

::

    ====================== 
    NPairhop 6
    ====================== 
    ========Pairhop ====== 
    ====================== 
       0     1  0.50000
       1     2  0.50000
       2     3  0.50000
       3     4  0.50000
       4     5  0.50000
       5     0  0.50000

File format
^^^^^^^^^^^

-  Line 1: Header

-  Line 2: ``[npairhop] [count]``

-  Lines 3-5: Header

-  Lines 6-: ``[i] [j] [val]``

Parameters
^^^^^^^^^^

-  ``[npairhop]``

   **Type :**
   String (blank is not allowed)

   **Description :**
   An arbitrary keyword for the total number of the pair-hopping couplings.

-  ``[count]``

   **Type :**
   Integer (blank is not allowed)

   **Description :**
   An integer giving the total number of the pair-hopping couplings.

-  ``[i]``, ``[j]``

   **Type :**
   Integer (blank is not allowed)

   **Description :**
   An integer giving a site index (:math:`0 \le i, j < {\rm Nsite}`).

-  ``[val]``

   **Type :**
   Float (blank is not allowed)

   **Description :**
   A value for :math:`J_{ij}^{\rm Pair}`.

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
