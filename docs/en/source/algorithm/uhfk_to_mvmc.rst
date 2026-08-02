.. highlight:: none

UHFk to mVMC PairProduct conversion
***************************************

``tools/uhfk_to_mvmc.py`` converts a self-consistent H-wave UHFk solution
into the initial PairProduct wave function of mVMC. This chapter states the
conventions and the construction. For usage, see
:doc:`../uhfk/tools/uhfk_to_mvmc`.

The conversion has four stages: extraction of the occupied Slater columns,
gauge composition under sublattice folding, construction of the pair
amplitude :math:`F`, and, when spin-orbit coupling is present, emission of
``trans.def``.

Bloch amplitude sign convention
===============================

The bridge adopts the **positive-Bloch amplitude convention**, that is, the
real-space amplitude is defined as

.. math::

   c_R = \frac{1}{\sqrt{N_\mathrm{folded}}} \sum_k c_k e^{+i k R}

This sign is fixed by H-wave's internal conventions and is not a free
choice. H-wave treats antiperiodic boundary conditions with the
negative-sign gauge

.. math::

   \tilde{c}_r = e^{-i \theta r / L_\mathrm{phys}} c_r

and uses ``ifftn(..., norm="forward")`` for the inverse Fourier transform.
Composing the two makes the Bloch transform on the tilde side positive.

For ``SubShape = [1, 1, 1]`` the :math:`(k, -k)` time-reversal sum
symmetrizes the plane-wave factor, so the density does not depend on the
choice of sign. For ``SubShape > [1, 1, 1]`` the sublattice envelope of the
folded eigenvectors is non-trivial, the two conventions are distinguishable,
and only the positive convention agrees element-wise with H-wave's
``greenone.dat``.

Sublattice folding and gauge composition
========================================

When ``SubShape`` divides ``CellShape`` in every direction, the lattice is
folded onto a supercell. The folded orbital index encodes the cell position
``sub_offset`` within the supercell together with the original orbital
index. Physical sites are reconstructed by decomposing that index.

Two phases multiply the Slater matrix under folding, and they are composed:

.. math::

   \exp\left[-i k_\mathrm{folded} \cdot (\mathrm{folded\_cell} + \mathrm{sub\_offset})\right]
   \times
   \exp\left[-i \theta \cdot r_\mathrm{phys} / L_\mathrm{phys}\right]

The first factor is the gauge of the intra-supercell offset, the second is
the antiperiodic boundary twist. ``build_slater_orbitals`` applies the
composed phase to the emitted Slater matrix.

On the verification side, ``gauge_lift`` uses the same composed transform to
lift H-wave's ``green_sublattice``, which is stored in the folded Bloch
basis, into the physical basis. Because both share one transform, the
density check can take the emitted Slater matrix itself as its subject.

Construction of the pair amplitude F
====================================

The one-body density matrix convention is

.. math::

   G_{ij} = \langle c^\dagger_i c_j \rangle = \left(\overline{A} A^{T}\right)_{ij}

where :math:`A` denotes the matrix of occupied Slater columns.

Canonical pairs
---------------

Pairs are emitted from one canonical row per unordered pair
:math:`\{k, \mathrm{partner}(k)\}`. The canonical row is defined as follows.

- For a self-pair (:math:`\mathrm{partner}(k) = k`), :math:`k` itself is
  canonical.
- Otherwise, the row whose wavevector index tuple is lexicographically
  smaller is canonical, with the smaller row index as tie-break.

:math:`F` is a :math:`(2N_\mathrm{site}, 2N_\mathrm{site})` complex
antisymmetric matrix; the lower triangle is filled by
:math:`F_{ji} = -F_{ij}`.

Partner balance
---------------

For the canonical pair emission to close consistently, the two rows of a
pair must carry the same occupation. The bridge rejects, before conversion,
any occupation for which either of the following holds.

- Some pair has :math:`n_\mathrm{occ}(k) \neq n_\mathrm{occ}(\mathrm{partner}(k))`.
- A canonical self-pair :math:`k` has odd :math:`n_\mathrm{occ}(k)`.

``validate_general_prerequisites`` performs this check and reports a failing
occupation as an error. It does not silently emit an approximate result.

When choosing a filling, the spectral gap alone is not sufficient. The gap
says nothing about how the occupation is distributed over the two rows of a
canonical pair, so the partner-balance invariant has to be checked
separately.

Projection to the density
-------------------------

Obtaining the one-body density from :math:`F` uses a skew-SVD projection.
The leading :math:`2 N_\mathrm{pairs}` left singular vectors span the
occupied subspace. This projection is used to verify that the emitted
:math:`F` represents the intended Slater state.

Transfer mapping convention
===========================

When spin-orbit coupling is present, the ``trans.def`` produced by mVMC's
``vmcdry.out`` keeps only the spin-diagonal entries and drops the
:math:`s \neq t` terms. The bridge therefore emits ``trans.def`` directly
from H-wave's ``Transfer.dat``.

The mapping is as follows. The mVMC convention is
:math:`H = -\sum \mathrm{trans}\, c^\dagger c`, and ComplexUHF builds its
bare Hamiltonian as :math:`K = -\mathrm{trans}`. In H-wave's negative-Bloch
physical basis, an entry :math:`(R, s, t, v)` of ``Transfer.dat`` maps at
displacement :math:`+R` as

.. math::

   K[i, t;\, i+R, s] &= \overline{v} \\
   \mathrm{trans}[i, t;\, i+R, s] &= -\overline{v}

The site endpoints stay :math:`i \to i+R`, the **spin endpoints are
swapped**, and the coefficient is conjugated and negated. On a real
spin-diagonal entry the swap is a no-op and the rule reduces to
:math:`\mathrm{trans} = -v`.

This general rule is derived. It follows the path from H-wave's
:math:`e^{+ikR}` transform under ``ifftn(norm="forward")``, through the
Hermitian transpose, to ComplexUHF's :math:`K = -\mathrm{trans}`. As
numerical support, reconstructing H-wave's bare :math:`K` from all stored
eigenpairs gives a maximum absolute deviation of
:math:`1.1 \times 10^{-12}` for this rule, with no entry above
:math:`10^{-10}`.

The older rule, which does not swap the spin endpoints (:math:`-v` for
:math:`s = t` and :math:`+v` for :math:`s \neq t`), was not derived; it was
pinned so that the numbers agreed, and it is a special case rather than an
equivalent. The two rules emit the same matrix only for matrices satisfying
:math:`v[t,s] = -\overline{v[s,t]}` at fixed :math:`R`. A system with
in-plane Rashba coupling alone satisfies that condition, so no difference
appears there, but a spin-symmetric hopping carrying both a real and an
imaginary part does not. In the same reconstruction the older rule gives a
maximum absolute deviation of :math:`6.0 \times 10^{-1}`.

The boundary wrap phase is applied **after** conjugation. A row crossing the
boundary along an antiperiodic direction is multiplied by
:math:`e^{+i\theta_d}` for a crossing at positive :math:`R` and
:math:`e^{-i\theta_d}` for a crossing at negative :math:`R`.

.. note::

   The verified range covers only :math:`\theta` whose components are in
   :math:`\{0, \pi\}`, where the wrap phase is real. A general complex
   twist is out of scope and has not been verified.
