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

The tabulated :math:`W` is complete under the exchange of its two bilinears --
an on-site Coulomb term, for instance, occupies both the up-down and the
down-up slots -- so the unrestricted sum over sites and over all four indices
contains each pair ordering twice, which is what the :math:`1/2` removes.  The
transfer table is likewise required to be Hermitian-closed, i.e. to satisfy
:math:`t_{\alpha\beta}(R)=t_{\beta\alpha}(-R)^{*}`, so :math:`{\cal H}_0`
carries no separate Hermitian conjugate.

Applying the Fourier transformation

.. math::
    \begin{aligned}
    c_{i\alpha}
    =\frac{1}{\sqrt{N_L}}\sum_{\bf{k}}
    e^{i \bf{k}\cdot \bf{r}_{i}}c_{\bf{k},\alpha}^{\mathstrut},
    \end{aligned}

the Hamiltonian is rewritten in the following form

.. note::

   This :math:`e^{+i\bf{k}\cdot\bf{r}}` convention (the Wannier90-style
   sign, shared with UHFk) is what every real-space-coefficient build
   (:math:`R \to k/q`) of the RPA module follows since issue #133:
   :math:`\varepsilon({\bf k}) = \sum_{\bf R} t({\bf R})
   e^{+i {\bf k}\cdot{\bf R}}` and :math:`W({\bf q}) = \sum_{\bf R}
   W({\bf R}) e^{+i {\bf q}\cdot{\bf R}}` (the convolution transforms
   inside the susceptibility machinery are self-inverse pairs and are not
   affected). Before that fix the non-spin-orbital path used the opposite
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

Each bilinear carries two independent band indices.  Keeping only its
band-diagonal part (:math:`\gamma_1=\gamma_2`) would drop the interband
contributions and make the result depend on the choice of basis inside a
degenerate subspace (see the discussion at the irreducible susceptibility
below), so no such restriction is made; H-wave keeps the interaction in the
orbital basis throughout and uses the diagonalization only to construct the
Green's functions.

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

and the irreducible susceptibility is the particle-hole bubble built from it,

.. math::
    \begin{aligned}
     X^{(0)\alpha\alpha', \beta\beta'}({\bf q},i\omega_m)=
      -\frac{T}{N_L}
      \sum_{{\bf k},n}
      G^{(0)\alpha\beta}({\bf k}+{\bf q}, i\omega_m+ i\epsilon_{n})
      G^{(0)\beta'\alpha'}({\bf k}, i\epsilon_{n}),
    \end{aligned}

so that the two band sums it contains are independent: the particle and the
hole may sit on different bands.  Tying them to a common band would drop the
interband particle-hole transitions, which carry most of the response away
from :math:`{\bf q}=0`.  On a six-site two-orbital chain the density response
is unchanged at :math:`{\bf q}=0` -- there the total density operator is
exactly band diagonal -- but loses 29% at :math:`q=\pi/3`, 71% at
:math:`2\pi/3`, and at the zone boundary the band-diagonal bubble is
essentially zero (0.004 against 0.788), because the states at :math:`k` and
:math:`k+q` are then nearly orthogonal within a band.  The restriction is
also not basis independent inside a degenerate subspace, and for
spin-degenerate bands it removes the spin-flip channel, breaking SU(2)
symmetry.  No such restriction is made.

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

   the same reordering that appears between Eqs. (16) and (20) of the
   H-wave paper (`arXiv:2308.00324 [cond-mat.str-el]
   <https://arxiv.org/abs/2308.00324>`_).  Both index pairs are affected.

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


Block-diagonal optimization
*****************************

When the interaction Hamiltonian has a block-diagonal structure
(e.g., due to spin conservation or orbital decoupling),
the RPA equation can be solved independently for each block,
significantly reducing the computational cost.

The block structure is detected automatically from the COMBINED
connectivity of the interaction matrix and the bare susceptibility.  Both
are needed: a block-diagonal interaction does not by itself make
:math:`\hat{I}+\hat{X}^{(0)}\hat{W}` block diagonal if :math:`X^{(0)}`
couples indices across blocks, as happens for spin-mixing bands.

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

The resulting vertex, for on-site interactions:

- ``CoulombIntra`` :math:`U`: :math:`W_{+-} = -U`
- ``CoulombInter`` :math:`V`: :math:`W_{+-} = -V`
- ``Hund`` :math:`J`: :math:`W_{+-} = 0`
- ``Exchange`` :math:`J`: :math:`W_{+-} = -(J + J^{\rm T})/2`
- ``Ising`` :math:`I`: :math:`W_{+-} = +I` (all Wannier90-like k-space solvers now read
  the Ising file in this normalization -- the UHFk factor-1/4 discrepancy
  was resolved with issue #106; the separate real-space UHFr reader keeps
  its S^z convention)
- ``PairLift`` :math:`J`: :math:`W_{+-} = 0`
- ``PairHop`` :math:`J`: :math:`W_{+-} = -J`

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
rejected earlier, at read time (issue #93: declaration files must be
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
- Block-diagonal optimization is applied automatically when possible.
- The ``squashed`` calculation scheme is also supported for spin-orbital systems.
- The ``Norbit`` value in the geometry file (``geom.dat``) is the **spin-orbital
  count** (= 2 × the number of physical orbitals = Wannier90 ``num_wann``), the same
  convention as UHFk.

.. note::

   **Migration (RPA):** the geometry ``norb`` for spin-orbital input is now the
   spin-orbital count; double any pre-existing RPA spin-orbital ``geom.dat``
   ``Norbit``.


.. [1] `K. Yoshimi, T. Kato, H. Maebashi, J. Phys. Soc. Jpn. 78, 104002 (2009). <https://journals.jps.jp/doi/10.1143/JPSJ.78.104002>`_
