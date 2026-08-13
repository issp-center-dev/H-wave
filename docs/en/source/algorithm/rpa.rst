.. highlight:: none

.. _algorithm_sec:

Random Phase Approximation
==========================

The random phase approximation (RPA) is a method to detect the response to the fluctuations of
one-body operators by the effect of electron correlations, starting from the non-interacting
state.
While in the UHF approximation an initial guess of the configuration is required, 
the RPA method enables to infer the ordered phase that emerges from the second-order transition.
H-wave implements RPA method using Matsubara frequency, and allows to compare with the dynamical
observables measured in the experiments by analytical continuation.

In the following, the algorithm is described.
In the RPA mode of H-wave, the Hamiltonian given below will be considered:

.. math::
    \begin{aligned}
     {\cal H}&={\cal H}_0+{\cal H}_{\rm int},\\
     {\cal H}_0&=\sum_{i\alpha;j\beta}
      t_{ij}^{\alpha \beta}c_{i\alpha}^{\dagger}
      c_{j\beta}^{\mathstrut},\\
     {\cal H}_{\rm int}&=\frac{1}{2}\sum_{ij}\sum_{\alpha, \alpha', \beta, \beta'}W_{ij}^{\beta\beta',\alpha\alpha'}
      c_{i\alpha}^{\dagger}c_{i\alpha'}c_{j\beta'}^{\dagger}c_{j\beta}
    \end{aligned}

The input table lists each interaction in both of its orderings: an on-site
Coulomb term, for instance, fills the up-down slot and the down-up slot alike.
Since the sum above runs over all four indices without restriction, it counts
every such term twice, and the :math:`1/2` cancels that.  The transfer table is
given in the same way, satisfying
:math:`t_{\alpha\beta}(R)=t_{\beta\alpha}(-R)^{*}`; both directions being
present already, no separate Hermitian conjugate is added to
:math:`{\cal H}_0`.

Applying the Fourier transformation

.. math::
    \begin{aligned}
    c_{i\alpha}
    =\frac{1}{\sqrt{N_L}}\sum_{\bf{k}}
    e^{i \bf{k}\cdot \bf{r}_{i}}c_{\bf{k},\alpha}^{\mathstrut},
    \end{aligned}

the Hamiltonian is rewritten in the following form

.. note::

   Every real-space-coefficient build (:math:`R \to k/q`) of the RPA
   module follows this :math:`e^{+i\bf{k}\cdot\bf{r}}` convention (the
   Wannier90-style sign, shared with UHFk):
   :math:`\varepsilon({\bf k}) = \sum_{\bf R} t({\bf R})
   e^{+i {\bf k}\cdot{\bf R}}` and :math:`W({\bf q}) = \sum_{\bf R}
   W({\bf R}) e^{+i {\bf q}\cdot{\bf R}}` (the convolution transforms
   inside the susceptibility machinery are self-inverse pairs and are not
   affected). In 1.0.x the non-spin-orbital path used the opposite
   sign throughout -- a self-consistent global :math:`{\bf k} \to
   -{\bf k}` relabeling, so its stored ``chi0q``/``chiq`` carried negated
   momentum labels whenever the tensor is not elementwise even under
   :math:`{\bf q} \to -{\bf q}` on the FFT grid -- while the
   spin-orbital path mixed the two signs between the transfer and the
   interaction, so its ``chiq`` was solved from :math:`\chi_0({\bf q})`
   with :math:`W(-{\bf q})`: for interactions with
   :math:`W({\bf q}) \neq W(-{\bf q})` (directional bonds) the old
   spin-orbital ``chiq`` is wrong, not merely relabeled. Momentum-space
   files written since this fix carry ``momentum_convention =
   "e_plus_ikR"``; loaders reject unmarked legacy files unless their
   content is :math:`{\bf q}`-even (for which the two conventions
   coincide).

.. math::
    \begin{aligned}
     {\cal H}&=\sum_{{\bf k}\alpha\beta}
     \varepsilon_{\alpha\beta}({\bf k})c_{{\bf k}\alpha}^{\dagger}
     c_{{\bf k}\beta}^{\mathstrut} \nonumber\\
    &+\frac{1}{2N_L}\sum_{{\bf k} {\bf k}'{\bf q}}\sum_{\alpha\beta\alpha'\beta'}
     W^{\beta\beta',\alpha\alpha'}_{{\bf q}}
     c_{{\bf k}+{\bf q},\alpha}^{\dagger}
      c_{{\bf k},\alpha'}^{\mathstrut}
      c_{{\bf k}'-{\bf q},\beta'}^{\dagger}
      c_{{\bf k}',\beta}^{\mathstrut}
    \end{aligned}

In the random phase approximation, the density fluctuation by the effect of electron correlation
is detected with respect to :math:`{\cal H}_0`.
The one-body part is diagonalized to build the Green's functions, and the
interaction term is expressed in that basis by the exact change of basis

.. math::
    \begin{aligned}
    &W^{\beta\beta',\alpha\alpha'}_{\bf{q}}c_{\bf{k}+\bf{q},\alpha}^{\dagger}c_{\bf{k},\alpha'}^{\mathstrut}
    c_{\bf{k}'-\bf{q},\beta'}^{\dagger} c_{\bf{k}',\beta}^{\mathstrut}\nonumber\\
    &= W^{\beta\beta',\alpha\alpha'}_{\bf{q}}
    \sum_{\gamma_1 \gamma_2 \gamma_1' \gamma_2'}
    u_{\alpha \gamma_1, \bf{k}+\bf{q}}^* d_{\bf{k}+\bf{q},\gamma_1}^{\dagger}
    u_{\alpha' \gamma_2, \bf{k}} d_{\bf{k},\gamma_2}^{\mathstrut}
    u_{\beta' \gamma_1', \bf{k}'-\bf{q}}^* d_{\bf{k}'-\bf{q},\gamma_1'}^{\dagger}
    u_{\beta  \gamma_2', \bf{k}'}d_{\bf{k}',\gamma_2'}^{\mathstrut}.
    \end{aligned}

Each bilinear carries two independent band indices.

Here, 

.. math::
    \begin{aligned}
    c_{\bf{k},\alpha} = \sum_{\gamma} u_{\alpha \gamma, \bf{k}} d_{\bf{k}, \gamma}
    \end{aligned}

and :math:`d_{\bf{k}, \gamma}` denotes the annihilation operator that diagonalizes :math:`{\cal H}_0`. (:math:`\gamma` refers to the index of the eigenvalue.)
Then, the irreducible one-body Green's function is written as

.. math::
    \begin{aligned}
     G^{(0)\alpha\beta}_{\gamma}({\bf k}, i\epsilon_{n})=
      \frac{u^{\alpha\gamma}({\bf k})u^{*\beta\gamma}({\bf k})}{i\epsilon_{n}-\xi^{\gamma}({\bf k})+\mu}.
    \end{aligned}

The full one-body Green's function is the sum of these band contributions,

.. math::
    \begin{aligned}
     G^{(0)\alpha\beta}({\bf k}, i\epsilon_{n})
      = \sum_{\gamma=1}^{n_{\rm orb}} G^{(0)\alpha\beta}_{\gamma}({\bf k}, i\epsilon_{n}),
    \end{aligned}

H-wave retains both sums.  The band index is contracted away as soon as it
appears, so the interaction is kept in the orbital basis throughout and the
diagonalization is used only to construct the Green's functions.

By using the irreducible susceptibility, the susceptibility matrix from the RPA
is obtained as follows:

.. math::
    \begin{aligned}
    X^{\alpha\alpha', \beta\beta'}(q)&=
    X^{(0)\alpha\alpha', \beta\beta'}(q)
    - \sum_{\alpha_1\alpha_1'\beta_1\beta_1'}
    X^{(0)\alpha\alpha', \beta_1\beta_1'}(q)
    W^{\beta_1'\beta_1, \alpha_1'\alpha_1}_{\bf q}
    X^{\alpha_1 \alpha_1' , \beta \beta'}(q),
    \end{aligned}

Combining indices such as :math:`\alpha\alpha^\prime` into one index, they are expressed
in the matrix form. Then finally it leads to the expression:

.. math::
    \begin{aligned}
     \hat{X}(q)&=\hat{X}^{(0)}(q)-\hat{X}^{(0)}(q)\hat{W}(q)\hat{X}(q)\nonumber\\
     &=\left[\hat{I}+\hat{X}^{(0)}(q)\hat{W}(q)\right]^{-1}\hat{X}^{(0)}(q).
    \end{aligned}

.. note::

   **Index-pair convention between the interaction and the susceptibility.**
   The two objects above label their index pairs with opposite orderings
   inside each pair, and converting between them is part of the RPA
   equation, not an afterthought.

   The two pair slots of the bare susceptibility carry the bilinears
   below -- the left pair with its creation index *second*, the right
   pair with its creation index *first* --

   .. math::
      X^{(0)\alpha\alpha',\beta\beta'}(q) \;\sim\;
      \Big\langle \big(c^{\dagger}_{\alpha'}c^{\mathstrut}_{\alpha}\big)(-q)\;;\;
      \big(c^{\dagger}_{\beta}c^{\mathstrut}_{\beta'}\big)(q)\Big\rangle ,

   which is what the Green-function product :math:`G^{\alpha\beta}(k+q)
   G^{\beta'\alpha'}(k)` above encodes. The interaction
   :math:`W^{\beta\beta',\alpha\alpha'}_{\bf q}` multiplies
   :math:`c^{\dagger}_{\alpha}c^{\mathstrut}_{\alpha'}
   c^{\dagger}_{\beta'}c^{\mathstrut}_{\beta}`, i.e. in *each* pair its
   ordering is the reverse of the susceptibility's.  The vertex entering
   the RPA equation above is therefore the interaction tensor with the
   two indices of *each pair* transposed -- written there as
   :math:`W^{\beta_1'\beta_1,\alpha_1'\alpha_1}_{\bf q}` -- i.e. as a
   matrix,

   .. math::
      \hat{W}(q)_{(\beta\beta'),(\alpha\alpha')}
      = W^{\beta'\beta,\,\alpha'\alpha}_{\bf q},

   with both index pairs reordered.

   The conversion is the identity on density (pair-diagonal) components
   and, more generally, on any real Hermitian-closed declaration, since
   the transposed slot then carries the same value.  It is observable
   only for a *complex* pair-crossing interaction -- a complex
   Hermitian-closed ``PairHop`` -- where omitting it returns the
   susceptibility of the complex-conjugate Hamiltonian.  H-wave applies
   the conversion when it assembles the vertex for the ring solve.  The
   transverse (ladder) assembly is *not* routed through it: that
   assembly re-pairs the interaction tensor itself and therefore
   consumes the Hamiltonian convention directly.  The stored
   ``chiq``/``chi0q`` and the interaction files keep the conventions
   documented in their own sections.

In the above formula, orbitals and spins were treated as unified generalised orbitals.
Of the arrays needed to perform the calculations,
the susceptibility ( :math:`X^{(0)\alpha\alpha^\prime, \beta\beta^\prime}({\bf q},i\omega_m)`, :math:`X^{\alpha\alpha^\prime, \beta\beta^\prime}({\bf q},i\omega_m)`) is the largest multidimensional array,
given by :math:`N_{\rm orb}^4 N_{\rm spin}^4 N_k N_{\omega}`, where the memory cost and computational complexity increase as the size increases.
As explained below, the size of the multidimensional array of susceptibilities can be reduced by separating orbits and spins:
for the two-body interactions handled in H-wave's RPA mode, separating orbits and spins results in

.. math::
    \begin{aligned}
    & W^{\beta\sigma_1\sigma_1',\alpha\sigma\sigma'}_{\bf{q}}c_{\bf{k}+\bf{q},\alpha \sigma}^{\dagger}c_{\bf{k},\alpha \sigma'}^{\mathstrut}
    c_{\bf{k}'-\bf{q},\beta\sigma_1'}^{\dagger} c_{\bf{k}',\beta\sigma_1}^{\mathstrut}.
    \end{aligned}

Since the scattering keeps the orbital index, the irreducible
susceptibility becomes

.. math::
    \begin{aligned}
     X^{(0)\alpha, \beta}_{\sigma\sigma'\sigma_1\sigma_1'}({\bf q},i\omega_m)=
      -\frac{T}{N_L}
      \sum_{{\bf k},n}
      G^{(0)\alpha\beta}_{\sigma\sigma_1'}({\bf k}+{\bf q}, i\omega_m+ i\epsilon_{n})
      G^{(0)\beta\alpha}_{\sigma_1\sigma'}({\bf k}, i\epsilon_{n}).
    \end{aligned}

The array size can be reduced to :math:`N_{\rm orb}^2 N_{\rm spin}^4 N_k N_{\omega}`.
Then susceptibility matrix by RPA is obtained as follows:

.. math::
    \begin{aligned}
    X^{\alpha, \beta}_{\sigma\sigma'\sigma_1\sigma_1'}(q)&=
    X^{(0)\alpha, \beta}_{\sigma\sigma'\sigma_1\sigma_1'}(q)
    - \sum_{\alpha_2,\alpha_3}
      \sum_{\sigma_2\sigma_2'\sigma_3\sigma_3'}
    X^{(0)\alpha, \alpha_2}_{\sigma\sigma'\sigma_2\sigma_2'}(q)
    W^{\alpha_2, \alpha_3}_{\sigma_2'\sigma_2, \sigma_3'\sigma_3}({\bf q})
    X^{\alpha_3, \beta}_{\sigma_3\sigma_3',\sigma_1\sigma_1'}(q).
    \end{aligned}

If :math:`\alpha\sigma\sigma'` is regarded as a single index,
it can be put into matrix form and, as in the case of generalised orbitals, can be used as a

.. math::
    \begin{aligned}
     \hat{X}(q)&=\hat{X}^{(0)}(q)-\hat{X}^{(0)}(q)\hat{W}(q)\hat{X}(q)\nonumber\\
     &=\left[\hat{I}+\hat{X}^{(0)}(q)\hat{W}(q)\right]^{-1}\hat{X}^{(0)}(q).
    \end{aligned}

The above formula is the general formula for the RPA method.

It is noted that the vertex correction may be taken into account as a means to consider
higher order correlations. See, for example, reference [1]_ for the details.


Sector structure and block decomposition
*****************************************

When the susceptibility decomposes into several blocks, the RPA equation can be
solved independently for each of them, which reduces the computational cost
substantially.

If :math:`{\cal H}_0` has a conserved quantity :math:`Q`, the Green's function is
block diagonal in it and the susceptibility separates into blocks carrying the
same label.  The label of each side of the correlator is the change
:math:`\Delta Q` its bilinear produces, the left-hand side being read as
:math:`A^{\dagger}` and labelled by :math:`A`; with that convention the selection
rule reads :math:`\Delta Q_{\rm L} = \Delta Q_{\rm R}`.  Spin is the usual
example: when :math:`[{\cal H}_0, S_z]=0` the density
(:math:`\Delta S_z=0`) and spin-flip (:math:`\Delta S_z=\pm 1`) channels are
strictly decoupled, which is what allows the longitudinal (ring) and transverse
(ladder) channels to be solved separately.  Solving block by block is not an
approximation: it returns the same result as solving the whole.

The decomposition holds only if the **interaction respects the same quantity**.
If :math:`X^{(0)}` is block diagonal but :math:`\hat{W}` is not, neither
:math:`\hat{X}^{(0)}\hat{W}` nor its inverse is, and the sectors mix from
second order in the interaction onwards.  Spin-orbit coupling makes
:math:`S_z` non-conserving and so puts mixing into :math:`X^{(0)}` itself,
while ``PairLift`` breaks the conservation on the interaction side even when
:math:`{\cal H}_0` preserves it.

H-wave therefore does not assume a block structure but detects it numerically,
from the connectivity of **both** the interaction matrix and the bare
susceptibility.

1. Sum the absolute values of the interaction Hamiltonian over all k-points,
   and those of the bare susceptibility over momentum and frequency, to
   obtain a connectivity pattern matrix.
2. Build an adjacency graph from non-zero off-diagonal entries (threshold: :math:`10^{-12}`).
3. Find connected components via label propagation (union-find algorithm).

If the matrix decomposes into :math:`m` blocks of sizes :math:`n_1, n_2, \ldots, n_m`,
the computational cost of solving the RPA equation is reduced from
:math:`O(N^3)` to :math:`O(n_1^3 + n_2^3 + \cdots + n_m^3)`,
where :math:`N = n_1 + n_2 + \cdots + n_m`.

This optimization is applied automatically and is transparent to the user.


Transverse susceptibility (ladder diagram)
*******************************************

In addition to the standard (ring diagram) RPA susceptibility,
H-wave can compute the transverse susceptibility
:math:`\chi_{+-}(\mathbf{q})`, which describes spin-flip correlations.
Writing it with the bilinears themselves avoids the ambiguity of the
:math:`S^{\pm}` labels: the stored array holds

.. math::

   \chi_{+-,\alpha\gamma;\beta\delta}(\mathbf{q}) =
   \Big\langle \big(c^{\dagger}_{\gamma\downarrow}
   c^{\mathstrut}_{\alpha\uparrow}\big)(-\mathbf{q}) \;;\;
   \big(c^{\dagger}_{\beta\uparrow}
   c^{\mathstrut}_{\delta\downarrow}\big)(\mathbf{q}) \Big\rangle ,

with the same index-pair convention as the longitudinal susceptibility.

.. _rpa_which_array:

Which array holds which channel
--------------------------------

:math:`\chi_{+-}` is delivered in its own array, ``chiq_pm``, and only when
``calc_type = "ring+ladder"`` (which requires ``calc_scheme = "general"``).
The array ``chiq`` always holds the longitudinal (ring) result.

Whether ``chiq`` has such slots at all depends on the scheme.  Under
``calc_scheme = "general"`` it is indexed by pairs that each carry a spin, so
it always has slots in which a pair is spin-off-diagonal.  ``"reduced"``
stores density-pair components only and has no such slots at all.

**Those slots are not the transverse susceptibility, and whenever the bubble
is obtained by inflating a spin-free or spin-diagonal one they are not
computed: they are identically zero.**  The inflation builds only the
same-spin slots.  This covers every calculation without
``enable_spin_orbital``, and also an ``enable_spin_orbital`` calculation whose
one-body Hamiltonian the solver detects as spin-free or spin-diagonal -- what
decides the question is the detected spin structure, not the setting.  A zero
in a susceptibility array reads naturally as a computed result; here it is an
absence.  On a two-orbital on-site model the longitudinal slots reach
:math:`\max|\chi| = 1.53` while every slot with a spin-off-diagonal pair is
exactly :math:`0`.

A genuinely spinful calculation under ``calc_scheme = "general"`` is the
exception: with ``enable_spin_orbital`` and a Hamiltonian that actually mixes
the spins, the bubble is built directly on the generalized orbital index, so
these slots are computed and are in general nonzero.  Even then the transverse
susceptibility is what ``chiq_pm`` holds, not what these slots hold.

Adding the ladder does not change the longitudinal answer.  On identical
input the ``chiq`` of ``"ring"`` and of ``"ring+ladder"`` agree bit for bit;
what ``"ring+ladder"`` adds is the separate ``chiq_pm``.  So the choice of
``calc_type`` is a choice of whether the transverse channel is computed at
all, not a choice of approximation for the longitudinal one.

.. note::

   In an SU(2)-symmetric paramagnet the transverse response follows from the
   longitudinal one by symmetry, :math:`\chi_{+-} = 2\chi_{zz}`, with
   :math:`\chi_{zz}` formed from the longitudinal result --

   .. math::

      \chi_{zz} = \frac{1}{4}\sum_{\sigma\sigma'}\sigma\sigma'\,
      \chi^{\sigma\sigma'} ,\qquad \sigma,\sigma' = \pm 1 ,

   where :math:`\chi^{\sigma\sigma'}` are the spin-diagonal-pair slots -- and
   not read off the spin-off-diagonal slots, which are empty.  Where the
   symmetry holds this reproduces ``chiq_pm`` to floating-point round-off
   (measured on the fixture below: largest **absolute** elementwise difference
   :math:`2.2\times 10^{-16}`, against a peak :math:`\chi_{+-}` of
   :math:`1.4`), so the ladder is not needed for it.

   The relation fails as soon as SU(2) is broken, and it fails without a
   small parameter: the discrepancy is not a correction that can be estimated
   and tolerated.  On a single-orbital square lattice (:math:`U = 4`,
   :math:`4\times 4`, :math:`N_{\rm mat} = 64`) with ``coeff_extern = 0.35``
   -- a splitting of :math:`0.7` between the two spin species -- the static
   inferred response deviates from ``chiq_pm`` at the peak wavevector by
   :math:`+2\%` at :math:`T = 0.5`, :math:`-4\%` at :math:`T = 1.0`,
   :math:`-25\%` at :math:`T = 0.3` and :math:`-37\%` at :math:`T = 0.2`.  It
   errs in both directions, it grows as the temperature falls, and at
   :math:`{\bf q} = 0`, :math:`T = 0.2` it is wrong by more than a factor of
   two.  At :math:`T = 0.3` the inferred maximum sits at a different
   wavevector from the true one, so even the position of the peak is not
   safe.  Under a field, a magnetic order parameter, or spin-orbit coupling,
   the transverse channel has to be computed.

.. admonition:: Reproducing the figures quoted above
   :class: note

   Both runs use ``calc_scheme = "general"``, ``mu = 0.0``,
   ``CellShape = [4, 4, 1]``, ``SubShape = [1, 1, 1]`` and ``Nmat = 64``,
   comparing ``calc_type = "ring"`` with ``calc_type = "ring+ladder"`` on
   otherwise identical input.

   The two-orbital figures use ``tests/rpa/input_2orb`` (``geom.dat``,
   ``transfer.dat``, ``coulombintra.dat``) at :math:`T = 0.5`.  Take
   ``CoulombIntra`` only: an off-site ``CoulombInter`` is rejected by
   ``calc_type = "ring+ladder"``.  The field figures use ``tests/rpa/input``
   (``geom.dat``, ``transfer.dat``, ``coulombintra.dat``, ``extern.dat``)
   with ``coeff_extern = 0.35``.

   Read the static response at Matsubara index :math:`N_{\rm mat}/2 = 32`,
   not at index :math:`0`: the frequency axis is
   :math:`\omega_n = (2n - N_{\rm mat})\pi/\beta`, so index :math:`0` is the
   most negative frequency.  The percentages are
   :math:`(\chi^{\rm inferred} - \chi_{+-})/\chi_{+-}` evaluated at the
   wavevector where :math:`\chi_{+-}` itself peaks.

The transverse bare susceptibility is

.. math::

   X^{(0)}_{+-,\alpha\gamma;\beta\delta}(\mathbf{q}, i\omega_m)
   = -\frac{T}{N_L} \sum_{\mathbf{k},n}
     G^{(0)}_{\alpha\beta,\uparrow}(\mathbf{k}+\mathbf{q}, i\omega_m + i\epsilon_{n})\,
     G^{(0)}_{\delta\gamma,\downarrow}(\mathbf{k}, i\epsilon_{n})

The transverse vertex :math:`W_{+-}` is obtained by crossing the
Hartree (Fock exchange) vertex from the longitudinal channel:

.. math::

   W_{+-} = W^{\rm spin-flip} - \left[W^{\rm cross-spin}\right]^{\rm crossed}

The interaction tensor enters this assembly in the HAMILTONIAN convention:
the index-pair conversion the ring solve applies (see the note in the
previous section) must NOT be applied here, because the crossing above
re-pairs the tensor itself.  The exact permutation is the one implemented in
``_assemble_transverse_vertex``.

The transverse vertex is built from the cross-spin and spin-flip blocks of
the interaction tensor only. The same-spin block does not enter: a same-spin
interaction cannot connect the up and down propagators of the transverse loop,
so it contributes self-energy but no vertex.

Each orbital pair is symmetrised with the mean of the two declarations, since
an interaction file may write the same operator either way
(:math:`n_a n_b = n_b n_a`, and :math:`X_{ab} = X_{ba}` for Exchange).
This matches the convention UHFk uses. The partner in the mean depends on the
interaction type (equivalently, on the slot family the type occupies):
density-density types and Exchange average with the plain transpose,
while PairHop averages with the conjugated transpose, because its two
declarations are Hermitian partners (:math:`P_{ba} = P_{ab}^{*}`) rather than
the same coefficient. For a complex Hermitian-closed Exchange the physical
coupling :math:`(J_{01} + J_{10})/2` is therefore real, while a complex
Hermitian-closed PairHop keeps its full complex value.

:math:`W_{+-}` is a matrix in the orbital pairs, and each interaction
contributes only to the slots it occupies.  The vertex draws on just two spin
blocks of the interaction tensor: the **cross-spin** block
(:math:`\uparrow\uparrow\downarrow\downarrow` and
:math:`\downarrow\downarrow\uparrow\uparrow`) and the **spin-flip** block
(:math:`\uparrow\downarrow\uparrow\downarrow` and
:math:`\downarrow\uparrow\downarrow\uparrow`).  The same-spin blocks do not
contribute: a same-spin interaction cannot connect the up and down propagators
of the transverse loop, so it produces self-energy but no vertex.

For on-site interactions the components are

.. list-table::
   :header-rows: 1
   :widths: 24 34 24

   * - interaction
     - block occupied
     - component of :math:`W_{+-}`
   * - ``CoulombIntra`` :math:`U`
     - cross-spin
     - :math:`-U`
   * - ``CoulombInter`` :math:`V`
     - cross-spin and same-spin (the latter unused)
     - :math:`-V`
   * - ``Ising`` :math:`I`
     - as above
     - :math:`+I`
   * - ``PairHop`` :math:`J`
     - cross-spin
     - :math:`-J`
   * - ``Exchange`` :math:`J`
     - spin-flip
     - :math:`-(J + J^{\rm T})/2`
   * - ``Hund`` :math:`J`
     - same-spin only
     - :math:`0`
   * - ``PairLift`` :math:`J`
     - double spin flip (not used by :math:`W_{+-}`)
     - :math:`0`

``Hund`` and ``PairLift`` vanish not because the coupling is weak but because
their components lie outside the two blocks the vertex uses.

On the ``Ising`` normalization: all Wannier90-like k-space solvers read the
Ising file this way; the UHFk factor-1/4 discrepancy was resolved in version
2.0, while the separate real-space UHFr reader keeps its :math:`S^z`
convention.

.. note::

   These values differ from those published before H-wave's transverse channel
   was checked against exact diagonalization; four of the earlier entries were
   incorrect and one type was missing. Transverse susceptibilities produced by
   earlier versions with ``CoulombInter``, ``Hund``, ``Ising`` or ``Exchange``
   should be recomputed. Only ``chiq_pm`` is affected -- it does not feed
   ``chiq``, the self-energy, or the Eliashberg vertex.

``calc_type = "ring+ladder"`` validates, before the longitudinal solve, that
the **assembled transverse vertex is independent of** :math:`q` on the
(sublattice-folded) lattice, to a relative tolerance of :math:`10^{-10}`;
input failing this is rejected. The transverse pair
:math:`c^\dagger_{i a \uparrow} c_{j b \downarrow}` is non-local for an
off-site term, so such a term's vertex is not a function of :math:`q` alone
and cannot be represented. In practice this rejects off-site
``CoulombInter``, ``Ising`` and ``Exchange``, while off-site ``Hund`` and
``PairLift`` are accepted because their transverse vertex vanishes. Note that a set of declarations whose members cancel or disagree is
rejected earlier, at read time (as of version 2.0: declaration files must be
Hermitian-closed); an inter-site pair that folds into the supercell under
``SubShape`` becomes an intra-cell orbital pair and is representable. The
longitudinal (``ring``) channel is unaffected.

.. warning::

   Off-site ``PairHop`` entries are silently discarded when the interaction is
   read, before this check runs, so they are neither rejected nor included.
   Do not rely on off-site ``PairHop`` in any RPA calculation. A DIAGONAL
   PairHop entry (equal orbitals) denotes the density term
   :math:`2P\, n_\uparrow n_\downarrow`; the interaction reader stores it
   with coefficient :math:`P` rather than :math:`2P`, consistently in both the
   longitudinal and transverse channels. Validation of such degenerate entries
   is tracked separately.

The transverse RPA susceptibility is obtained as

.. math::

   \hat{X}_{+-}(\mathbf{q})
   = \left[\hat{I} + \hat{X}^{(0)}_{+-}(\mathbf{q})\, \hat{W}_{+-}\right]^{-1}
     \hat{X}^{(0)}_{+-}(\mathbf{q})

To enable the transverse channel calculation, set ``calc_type = "ring+ladder"``
in the input TOML file. This requires the ``general`` calculation scheme
(automatically selected).

.. note::

   In the spin-orbital mode, when the Hamiltonian genuinely mixes spins
   (e.g. spin-orbit coupling), the transverse channel extracts only the
   :math:`S_z`-conserving block :math:`G_\uparrow G_\downarrow` of the
   bubble; the spin-mixing cross terms are not included, and a warning
   is emitted. The transverse susceptibility of a spin-mixing system is
   therefore an approximation in the current implementation.


Spin-orbital mode
*****************************

H-wave supports a spin-orbital mode, in which the input files index spin
and orbital in interleaved rather than block-separated form.  The solver
remaps that input to its internal spin-block order, which is also the
order the stored susceptibilities carry (they are marked
``index_convention = "spin_block"``).

In the normal mode, the combined index is :math:`i = s \cdot n_{\rm orb} + a`
(spin-block first), where :math:`s = 0, 1` is the spin index and
:math:`a = 0, \ldots, n_{\rm orb}-1` is the orbital index.
In the spin-orbital mode, the index is :math:`i = 2a + s`
(interleaved), which naturally accommodates spin-orbit coupling.

The spin-orbital mode is activated by setting ``enable_spin_orbital = true``
in the input TOML file. In this mode:

- The Hamiltonian is constructed in the full :math:`2n_{\rm orb} \times 2n_{\rm orb}` space
  without assuming spin conservation.
- All interaction types (``CoulombIntra``, ``CoulombInter``, ``Hund``, ``Exchange``,
  ``Ising``, ``PairLift``, ``PairHop``) are supported.
- The block decomposition is applied automatically when possible.
- Both the ``reduced`` and ``general`` calculation schemes are supported for spin-orbital systems.
- The ``Norbit`` value in the geometry file (``geom.dat``) is the **spin-orbital
  count** (= 2 × the number of physical orbitals = Wannier90 ``num_wann``), the same
  convention as UHFk.

.. note::

   **Migration (RPA):** the geometry ``norb`` for spin-orbital input is now the
   spin-orbital count; double any pre-existing RPA spin-orbital ``geom.dat``
   ``Norbit``.


.. [1] `K. Yoshimi, T. Kato, H. Maebashi, J. Phys. Soc. Jpn. 78, 104002 (2009). <https://journals.jps.jp/doi/10.1143/JPSJ.78.104002>`_
