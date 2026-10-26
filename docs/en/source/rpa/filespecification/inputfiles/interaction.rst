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
      :math:`\frac{1}{2}\sum_{ij\alpha\beta} V_{\alpha\beta}(r_{ij})\,n_{i\alpha} n_{j\beta}, \quad n_{i\alpha}=n_{i\alpha\uparrow}+n_{i\alpha\downarrow}`
    **Hund**:
      :math:`-\frac{1}{2}\sum_{ij\alpha\beta} J_{\alpha\beta}^{\rm Hund}(r_{ij}) \left( n_{i\alpha\uparrow} n_{j\beta\uparrow} + n_{i\alpha\downarrow} n_{j\beta\downarrow} \right)`
    **Ising**:
      :math:`\frac{1}{2}\sum_{ij\alpha\beta} J_{\alpha\beta}^{\rm Ising}(r_{ij}) (n_{i\alpha\uparrow} - n_{i\alpha\downarrow})(n_{j\beta\uparrow} - n_{j\beta\downarrow})`
    **PairHop**:
      :math:`\frac{1}{2}\sum_{ij\alpha\beta} \left[ J_{\alpha\beta}^{\rm PH}(r_{ij})\,c_{i\alpha\uparrow}^{\dagger} c_{j\beta\uparrow}^{\phantom{\dagger}} c_{i\alpha\downarrow}^{\dagger} c_{j\beta\downarrow}^{\phantom{\dagger}} + \textit{h.c.} \right]`
    **Exchange**:
      :math:`\sum_{ij\alpha\beta} J_{\alpha\beta}^{\rm Ex}(r_{ij})\,c_{i\alpha\uparrow}^\dagger c_{j\beta\uparrow}^{\phantom{\dagger}} c_{j\beta\downarrow}^\dagger c_{i\alpha\downarrow}^{\phantom{\dagger}}`
    **PairLift**:
      :math:`\frac{1}{2}\sum_{ij\alpha\beta} \left[ J_{\alpha\beta}^{\rm PairLift}(r_{ij})\,c_{i\alpha\uparrow}^{\dagger} c_{i\alpha\downarrow}^{\phantom{\dagger}} c_{j\beta\uparrow}^{\dagger} c_{j\beta\downarrow}^{\phantom{\dagger}} + \textit{h.c.} \right]`

.. note::

   Some of the expressions above were corrected to match what the solvers
   have always implemented; **no input file needs to be changed**, and no
   numerical result changes.

   - **Factor 1/2** (CoulombInter, Hund, Ising, PairHop, PairLift): the
     sums run over every ordered pair, i.e. over every declared row. Each
     row :math:`(R, \alpha, \beta)` with value :math:`X` must be
     accompanied by its partner :math:`(-R, \beta, \alpha)` with value
     :math:`X^{*}` (see below), and the two describe the same bond. For
     these types the partner row gives the same term as the row itself
     (for PairHop and PairLift, the same Hermitian-completed bracket), so
     each row carries half of it: a bond declared as a CoulombInter pair
     with value :math:`V` contributes :math:`V n_{i\alpha} n_{j\beta}`
     once. This is the convention of RESPACK ``zvo_ur.dat``, which lists
     both orientations of every bond. The factor was missing from earlier
     editions of this page. **Exchange** needs no factor: its summand has
     no explicit ``h.c.``, and the partner row supplies the Hermitian
     conjugate of the row's term, so the sum as written gives the bond
     term once.

   - **Hund**: the implemented convention carries a MINUS sign, as written
     above. A ferromagnetic (same-spin attractive) Hund coupling is
     therefore declared with a POSITIVE coefficient
     :math:`J^{\rm Hund}_{\alpha\beta}(r_{ij})`.
   - **PairLift**: each declared row enters together with its Hermitian
     conjugate, exactly as ``PairHop`` does; the ``+ h.c.`` above was
     missing from earlier editions of this page.

.. _rpa_interaction_orientation:

.. note::

   **Convention.** In every expression above :math:`r_{ij} = R_j - R_i`: a row
   ``[rx] [ry] [rz] [alpha] [beta] ...`` places the orbital ``[alpha]`` in the
   original cell and ``[beta]`` in the cell displaced by
   :math:`\vec{r} = (r_x, r_y, r_z)`, for the one-body and the two-body files
   alike.

   **Correction of the bond placement (issue #232).** A two-body row and a
   one-body row with the same ``(R, a, b)`` act on the same bond, as the
   Wannier90 ``hr.dat`` and RESPACK ``zvo_*r.dat`` formats define them.
   Earlier versions did not do this everywhere: the RPA solver (released
   versions up to and including v1.0.1 among them), the FLEX vertex
   (including ``flex_second_order = "local"`` and the bond-resolved
   channels) and the Eliashberg vertices put an off-site two-body row on
   the bond OPPOSITE its one-body row. The released ``UHFk`` placed it
   correctly (development versions between issues #193 and #232, which
   were never part of a release, had moved it, and the Hartree-Fock term
   of ``FLEX``, to the opposite bond as well). From the first release that
   includes issue #232, all solvers follow the convention above. Results
   change only for off-site rows that are not symmetric under
   :math:`R \to -R` (:math:`X_{ab}(R) \neq X_{ab}(-R)`), the generic case
   of multi-orbital RESPACK input: a bond between two different orbitals
   :math:`a \neq b` is declared as :math:`X_{ab}(R)` together with its
   partner :math:`X_{ba}(-R)`, and :math:`X_{ab}(-R)` -- the same pair of
   orbitals on the mirrored bond -- is in general a different number;
   every on-site interaction, every real single-orbital input and every
   table with :math:`X_{ab}(R) = X_{ab}(-R)` is unaffected
   (a complex off-site ``PairHop`` is affected even with a single orbital,
   because its coefficient at :math:`-R` is the complex conjugate).
   **No input file needs to be changed.**

   **Reproducing earlier results.** An earlier RPA, FLEX-vertex or
   Eliashberg result is reproduced by negating the displacement of every
   off-site two-body row (``[rx] [ry] [rz]`` -> ``[-rx] [-ry] [-rz]``; not
   the one-body rows such as ``Transfer``, not the on-site rows, and
   without swapping the orbital indices, which is wrong for complex
   couplings) -- that negated file describes the Hamiltonian those versions
   actually solved. The exception is the bond-resolved pairing kernel
   (``longitudinal_bond_pairing`` and the dynamic bond path of
   ``hwave_sc``): with more than one orbital it also coupled the gap to the
   bond channels with a wrong pair geometry, corrected at the same time, so
   multi-orbital results of that kernel can change even for an off-site
   table that is symmetric under :math:`R \to -R`, and cannot be
   reproduced from a modified file. Released ``UHFk`` results are
   unchanged. An earlier
   ``FLEX`` run with ``flex_hartree_fock = true`` evaluated the mean-field
   term and the vertices on different bonds and cannot be reproduced with
   a single file. Keep such reproduction runs in a separate directory:
   their outputs carry the current version stamps although they describe
   the old Hamiltonian. FLEX susceptibility, self-energy and Green-function
   files now carry ``sc_vertex_version = 3`` and the dynamic bond archive
   ``bond_archive_schema = 3``; the Eliashberg solvers and a FLEX
   ``sigma_init`` seed refuse an older file when the interaction declares
   an off-site row that is not symmetric under :math:`R \to -R`, so
   regenerate those files by re-running FLEX. Remove any workaround added
   for the old placement, and restart affected self-consistent runs from
   scratch rather than from earlier seeds.


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

-  The unspecified elements of the coefficient matrix are assumed to be zero. Note that a declared entry whose Hermitian partner :math:`X_{ba}(-R)` is unspecified is rejected at read time (as of version 2.0): both directions of a coupling must be declared. A nonzero on-site CoulombInter row of an orbital with itself (:math:`R = 0`, :math:`\alpha = \beta`) is also rejected: each declared row carries half of :math:`V n n`, so such a row is the Hubbard term :math:`U = V` plus a one-body level shift :math:`V/2`; declare :math:`U = V` in CoulombIntra and add :math:`V/2` to the orbital's on-site energy in the Transfer file (or absorb it in the chemical potential).

-  The translation vectors need to be enclosed within the CellShape. If the range of ``r_x``, ``r_y``, or ``r_z`` exceeds the extent of ``x``, ``y``, or ``z`` dimension of CellShape, the program terminates with an error.

-  When ``mode.enable_spin_orbital`` is set to ``true``, the orbital indices of the
   Transfer term use the interleaved spin-orbital convention and range from 1 to
   :math:`2 N_\text{orbital}`.  The combined index is :math:`2\alpha + s + 1`
   (1-based), where :math:`\alpha = 0, \ldots, N_\text{orbital}-1` is the
   physical-orbital index and :math:`s \in \{0,1\}` is the spin index.
   Otherwise, only the entries with the orbital indices from 1 to
   :math:`N_\text{orbital}` are taken into account.

   .. note::

      **Migration (RPA):** the geometry ``norb`` for spin-orbital input is now the
      spin-orbital count (= 2 × physical orbitals); double any pre-existing RPA
      spin-orbital ``geom.dat`` ``Norbit``.

.. raw:: latex
