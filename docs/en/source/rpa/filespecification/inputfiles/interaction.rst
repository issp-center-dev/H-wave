.. highlight:: none

Interaction definition files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The interaction definition files describe the coefficients
:math:`T_{\alpha\beta}(r_{ij})`, :math:`J_{\alpha\beta}(r_{ij})`, :math:`V_{\alpha\beta}(r_{ij})`, or :math:`U_{\alpha}`
of the one-body and two-body Hamiltonian denoted by the following expressions. 
They are given in Wannier90(-like) format.
It is noted that the generalized two-body interaction term (InterAll) is not supported
in the random phase approximation, 
    
    **Transfer**:
      :math:`\sum_{ij\alpha\beta\sigma} T_{\alpha\beta}(r_{ij})\,c_{i\alpha\sigma}^{\dagger}c_{j\beta\sigma}^{\phantom{\dagger}}`
    **Extern**:
      :math:`\sum_{ij\alpha\beta\sigma_1\sigma_2} H_{\alpha\beta}(r_{ij})\,\sigma^z_{\sigma_1\sigma_2}\,c_{i\alpha\sigma_1}^{\dagger}c_{j\beta\sigma_2}^{\phantom{\dagger}}, \quad \sigma^z = \text{diag}(1,-1)`
    **CoulombIntra**:
      :math:`\sum_{i\alpha} U_\alpha\,n_ {i\alpha\uparrow} n_{i\alpha\downarrow}, \quad n_{i\alpha\sigma}=c_{i\alpha\sigma}^{\dagger}c_{i\alpha\sigma}^{\phantom{\dagger}}`
    **CoulombInter**:
      :math:`\sum_{ij\alpha\beta} V_{\alpha\beta}(r_{ij})\,n_{i\alpha} n_{j\beta}, \quad n_{i\alpha}=n_{i\alpha\uparrow}+n_{i\alpha\downarrow}`
    **Hund**:
      :math:`\sum_{ij\alpha\beta} J_{\alpha\beta}^{\rm Hund}(r_{ij}) \left( n_{i\alpha\uparrow} n_{j\beta\uparrow} + n_{i\alpha\downarrow} n_{j\beta\downarrow} \right)`
    **Ising**:
      :math:`\sum_{ij\alpha\beta} J_{\alpha\beta}^{\rm Ising}(r_{ij}) S^{z}_{i\alpha} S^{z}_{j\beta}, \quad S^{z}_{i\alpha}=\frac{1}{2}(n_{i\alpha\uparrow} - n_{i\alpha\downarrow})`
    **PairHop**:
      :math:`\sum_{ij\alpha\beta} J_{\alpha\beta}^{\rm PH}(r_{ij})\,c_{i\alpha\uparrow}^{\dagger} c_{j\beta\uparrow}^{\phantom{\dagger}} c_{i\alpha\downarrow}^{\dagger} c_{j\beta\downarrow}^{\phantom{\dagger}} + \textit{h.c.}`
    **Exchange**:
      :math:`\sum_{ij\alpha\beta} J_{\alpha\beta}^{\rm Ex}(r_{ij})\,c_{i\alpha\uparrow}^\dagger c_{j\beta\uparrow}^{\phantom{\dagger}} c_{j\beta\downarrow}^\dagger c_{i\alpha\downarrow}^{\phantom{\dagger}}`
    **PairLift**:
      :math:`\sum_{ij\alpha\beta} J_{\alpha\beta}^{\rm PairLift}(r_{ij})\,c_{i\alpha\uparrow}^{\dagger} c_{i\alpha\downarrow}^{\phantom{\dagger}} c_{j\beta\uparrow}^{\dagger} c_{j\beta\downarrow}^{\phantom{\dagger}}`


An example of the file is shown below.

::

   wannier90 format for vmcdry.out or HPhi -sdry
       10
      245
    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1
    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1
   ...
    1    1    1    1    1
   -3   -3   -2    1    1  -0.0000269645  -0.0000000000
   -3   -3   -2    1    2  -0.0000071722  -0.0000018600
   -3   -3   -2    1    3  -0.0000083990   0.0000010972
   -3   -3   -2    1    4  -0.0000000990   0.0000000427
   -3   -3   -2    1    5  -0.0000018628  -0.0000003609
   -3   -3   -2    1    6  -0.0000129504  -0.0000014047
   -3   -3   -2    1    7  -0.0000189169   0.0000024697
   -3   -3   -2    1    8   0.0000238115   0.0000014316
   -3   -3   -2    1    9   0.0000036708  -0.0000003266
   -3   -3   -2    1   10   0.0000361752   0.0000003247
   -3   -3   -2    2    1  -0.0000071722   0.0000018600
   -3   -3   -2    2    2   0.0000105028  -0.0000000000
   ...


File format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  Line 1: Header

-  Line 2: ``[Norbit]``

-  Line 3: ``[Npts]``

-  Lines 4 - :math:`\lceil N_\text{pts} / 15 \rceil + 3`:

      ``[n1] [n2] ...``

-  Line :math:`\lceil N_\text{pts} / 15 \rceil + 4` onwards:

      ``[rx] [ry] [rz] [alpha] [beta]  [J.real] [J.imag]``

Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  ``[Norbit]``

   **Type :**
   Integer

   **Description :**
   This parameter specifies the number of orbitals :math:`N_\text{orbit}` in a unit cell.

-  ``[Npts]``

   **Type :**
   Integer

   **Description :**
   This parameter specifies the number of cells in a rectangular cuboid
   that accommodates entire translation vectors.

-  ``[n1]``, ``[n2]``, ...

   **Type :**
   Integer

   **Description :**
   These parameters specify the multiplicity of cells (ordinary 1), 
   with 15 points in a line.

-  ``[rx]``, ``[ry]``, ``[rz]``

   **Type :**
   Integer

   **Description :**
   These parameters specify the translation vector.
   
-  ``[alpha]``, ``[beta]``

   **Type :**
   Integer

   **Description :**
   These parameters specify the indices of the orbitals.
   ``[alpha]`` corresponds to the orbital :math:`\alpha` in the original cell,
   and ``[beta]`` corresponds to the orbital :math:`\beta` in the cell displaced
   by :math:`\vec{r}`.

-  ``[J.real]``, ``[J.imag]``

   **Type :**
   Float

   **Description :**
   These parameters specify the real and imaginary parts of the coefficient
   :math:`J_{\alpha\beta}(\vec{r})`.


Usage rules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  Header cannot be omitted.

-  The unspecified elements of the coefficient matrix are assumed to be zero.

-  The translation vectors need to be enclosed within the CellShape. If the range of ``r_x``, ``r_y``, or ``r_z`` exceeds the extent of ``x``, ``y``, or ``z`` dimension of CellShape, the program terminates with an error.

-  When ``mode.enable_spin_orbital`` is set to ``true``, the orbital indices of Transfer term are interpreted as the extended orbital indices including spin degree of freedom that ranges from 1 to :math:`2 N_\text{orbital}`. Otherwise, only the entries with the orbital indices from 1 to :math:`N_\text{orbital}` are taken into account.

.. raw:: latex
