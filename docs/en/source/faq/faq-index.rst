.. _Ch:FAQ:

FAQ
****

Superconductivity / Eliashberg equation usage
=================================================================

Separating the spin and charge contributions to superconductivity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Q. I want to analyze superconductivity by separating the charge and spin
contributions. Is that possible?**

Yes. In the dynamic Eliashberg solver (``[eliashberg] frequency = "dynamic"``,
which reads the FLEX susceptibilities) the pairing vertex is built from the spin
(:math:`\chi_{\mathrm s}`) and charge (:math:`\chi_{\mathrm c}`)
susceptibilities and is **linear** in the two channels. You can therefore
compare each channel's effect on the pairing eigenvalue :math:`\lambda` by
zeroing the other channel in the ``[eliashberg]`` section:

- ``zero_chi_c = true`` keeps the spin-fluctuation term plus the bare term;
- ``zero_chi_s = true`` keeps the charge-fluctuation term plus the bare term.

Running the full calculation together with the two zeroed calculations and
comparing the resulting :math:`\lambda` provides a diagnostic indication of
which fluctuation channel has the stronger effect on pairing. Both flags default
to ``false`` (ordinary runs are unaffected), and the diagnostic works for both
``pairing_type = "singlet"`` and ``"triplet"``.

.. note::

   The instantaneous (bare) interaction term is retained in every case, and the
   leading eigenvalue is not a linear function of the pairing vertex, so the
   eigenvalues do **not** add up: :math:`\lambda_{\mathrm s} +
   \lambda_{\mathrm c} \neq
   \lambda_{\mathrm{full}}` in general, where :math:`\lambda_{\mathrm s}` and
   :math:`\lambda_{\mathrm c}` denote the ``zero_chi_c = true``
   (spin-plus-bare) and ``zero_chi_s = true`` (charge-plus-bare) runs,
   respectively. The bare term is present in both and would also be counted
   twice by a simple sum. Use the decomposition to compare the relative
   strength of the two channels, not as an exact additive split.

See :ref:`sc_channel_decomposition` for the vertex formulas and further details.

Estimating :math:`T_c` when :math:`\lambda` does not reach 1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Q. My leading eigenvalue never reaches 1. How do I estimate** :math:`T_c`\ **?**

:math:`\lambda = 1` is the self-consistency condition for the superconducting
transition, and for a fixed pairing channel :math:`\lambda` grows as the
temperature is lowered, so :math:`\lambda < 1` at a given :math:`T` simply
means that :math:`T` is above :math:`T_c` -- run at lower temperatures rather
than treating it as a failure. ``hwave_tsweep`` automates a descending
temperature ladder, chaining each rung's converged self-energy and (for the
dynamic solver) leading eigenvector into the next; see
:ref:`sc_tsweep`. Because a FLEX+Eliashberg :math:`T_c` is only
semi-quantitative (see :ref:`faq_flex_tc`), comparing the position and shape
of :math:`\lambda(T)` or :math:`\lambda(P)` across a parameter sweep is more
robust than reading off the absolute value where :math:`\lambda` crosses 1.

FLEX does not converge at low temperature (cold-start failure)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Q. My FLEX run at low** :math:`T` **fails to converge, even though a similar
run at higher** :math:`T` **works. What should I do?**

A cold start (:math:`\Sigma = 0`) directly at low :math:`T` often fails to
converge, because the SCF fixed point moves further from :math:`\Sigma = 0`
as correlations strengthen. Instead, start from a higher temperature where
convergence is easy and step down, feeding each temperature's converged
self-energy into the next run as a warm start (``[file.input] sigma_init``,
see :ref:`flex_sigma_init`); ``hwave_tsweep`` automates this whole ladder
(:ref:`sc_tsweep`).

The leading eigenvalue or gap symmetry jumps between temperatures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Q. Stepping through a temperature sweep, the leading eigenvalue's gap
symmetry suddenly switches to a different one. Is my calculation broken?**

Not necessarily. The dynamic (frequency-dependent) Eliashberg kernel is
non-Hermitian and can have *exceptional points* where two real eigenvalues
collide and split into a complex-conjugate pair, so the "leading" branch can
jump discontinuously between neighbouring temperatures even when the FLEX
self-energy varies smoothly. Set ``[eliashberg] seed_eigenvector`` to a
neighbouring run's ``gap_dynamic.npz`` to select the eigenpair whose
eigenvector best overlaps the seed instead of the algebraically largest one,
which follows one physical branch (e.g. the d-wave mode) continuously across
the sweep. See :ref:`sc_seed_eigenvector` for the full mechanism, including
``sigma_shift`` for a masked or complexifying eigenvalue.

How do I identify the gap symmetry I obtained? (the ``match`` column)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Q. I set** ``init_gap`` **to seed a particular symmetry, but I am not sure
the solver actually converged to it. How do I check?**

``init_gap`` (static solver) only seeds the self-consistent iteration; the
reported result is the leading eigenpair of the requested channel *parity*
(even for singlet, odd for triplet), which need not be the symmetry you
seeded, and the eigenvalue analysis also reports opposite-parity solutions
that are mathematically valid but unphysical for that channel. The trailing
``match`` column of ``eigenvalue.dat`` marks this: ``1`` when the eigenvector
lies in the requested channel-parity sector (physical), ``0`` when it lies in
the opposite sector (spurious). Inspect the reported gap function itself
(``gap.dat`` / the figures in the tutorial) to identify which symmetry was
actually obtained. See :ref:`sc_eigenvalue_dat` for the file format and the
parity note above it for why opposite-parity solutions appear at all.

:math:`\lambda(T)` is non-monotonic -- is that physical?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Q. For a fixed pairing channel,** :math:`\lambda` **dips and then rises
again as I lower the temperature. Is that a real effect?**

For a fixed pairing channel, :math:`\lambda` should increase monotonically as
:math:`T \to 0`. A dip usually signals a convergence problem -- a cold start
or an under-converged FLEX self-energy at that rung rather than a warm-started
one (:ref:`flex_sigma_init`) -- or a branch switch at an exceptional point
(:ref:`sc_seed_eigenvector`), not new physics. Check monotonicity across the
full temperature ladder (:ref:`sc_tsweep`) and use ``seed_eigenvector`` to
confirm you are tracking one physical branch before concluding the dip is
real.

Dynamic-mode prerequisites and choosing ``ir_wmax``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Q. What do I need before I can use** ``frequency = "dynamic"``\ **, and how
do I choose** ``ir_wmax``\ **?**

``frequency = "dynamic"`` requires frequency-resolved input that only a
preceding FLEX run (``chi0q_mode = "flex"``) produces -- see
:ref:`sc_dynamic_frequency` for the required files and the ``Nmat``
consistency checks. On the IR path (``matsubara_basis = "ir"``) you also need
the optional `sparse-ir <https://sparse-ir.readthedocs.io>`_ package. Choose
``ir_wmax`` large enough that :math:`\lambda` is converged with respect to it
(check against a larger value or the automatic estimate); making it too large
makes the delta(:math:`\tau`)-derived constant that the IR loader has to
isolate ill-conditioned. See :ref:`sc_dynamic_ir_en` for the automatic
estimate, the isolation diagnostic, and measured eigenvalue sensitivity.

GPU out-of-memory at large ``Nmat``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Q. My dynamic Eliashberg run with** ``gpu = true`` **aborts with an
out-of-memory error at large** ``Nmat``\ **. What can I do?**

The two large resident device tensors (the pair bubble and the pairing
vertex) scale as :math:`N_{\mathrm{orb}}^4 N_k N_{\mathrm{mat}}`, so GPU
memory use grows directly with ``Nmat`` (and with the k-mesh size and orbital
count). Reduce ``Nmat`` or ``CellShape``, or run fewer concurrent solves on
the same device. Before the large allocation, the solver logs an advisory
estimate of the required device memory against the free amount; check that
log line against your device's limit before assuming the run will fit. See
:ref:`sc_dynamic_gpu_en` for the memory formula and a reference speed/memory
point.

Two-orbital ``chi_convention`` (``kuroki`` vs. ``myo``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Q. My two-orbital FLEX+Eliashberg results look wrong even though the file
shapes match what I expect. What should I check?**

The FLEX susceptibility files carry a ``chi_convention`` tag (``"kuroki"`` for
the reduced scheme, ``"myo"`` for the general full-vertex scheme) that the
Eliashberg loader uses to interpret the stored orbital layout. For
:math:`N_{\mathrm{orb}} = 2` specifically, the reduced spin-orbital dimension
and the general orbital-pair dimension coincide (both equal 4), so the array
*shape alone* cannot distinguish them -- only the ``chi_convention`` tag can.
Combining files with the wrong or mismatched convention silently builds an
incorrect pairing vertex rather than raising an error. See
:ref:`sc_dynamic_frequency` for the full explanation and the historical
mislabeling this tag now guards against.


FLEX approximation
=================================================================

Validity range of FLEX
^^^^^^^^^^^^^^^^^^^^^^^^

**Q. What are the limitations of FLEX? When should I not trust it?**

FLEX is a conserving, weak-to-intermediate-coupling approximation built from
the RPA particle-hole bubble and ladder series; it does not include the
Aslamazov-Larkin/Maki-Thompson vertex corrections, and it cannot describe a
Mott insulating state. Near a magnetic or charge ordering instability it
tends to overestimate the fluctuations rather than entering the ordered
phase, so derived transition temperatures are semi-quantitative at best (see
:ref:`faq_flex_tc`). See :ref:`flex_scope` for the precise scope, including
the density-density (``reduced``) vs. full-vertex (``general``) distinction.

My FLEX run diverges or will not converge near an instability
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Q. As I approach what looks like a magnetic or charge instability, FLEX
stops converging (or the susceptibility blows up). What is happening?**

As the system approaches a magnetic or charge instability, the RPA-type
denominator (the Stoner factor) approaches 1, :math:`\chi` grows large, and
the self-consistency loop becomes harder to converge. Approach the region
gradually with a warm-start temperature ladder (:ref:`flex_sigma_init`,
:ref:`sc_tsweep`) and gentler mixing -- see :ref:`flex_tips` for concrete
``Mix``/``mixing_scheme`` suggestions. If the instability is physical
(genuine magnetic or charge order), keep in mind that FLEX itself cannot
enter the ordered or insulating phase: a large susceptibility or pairing
eigenvalue in that region is a sign of the approximation nearing its limits,
not a reliable prediction of a superconducting state there.

.. _faq_flex_tc:

How should I read a FLEX-derived :math:`T_c`?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Q. My FLEX+Eliashberg calculation gives a specific** :math:`T_c` **value.
How much should I trust that number?**

Treat it as semi-quantitative, not quantitative. The absolute :math:`T_c`
depends near-exponentially on the interaction strength, so it is sensitive
both to how the interaction was calibrated and to the approximation itself
(:ref:`flex_scope`). Compare trends -- the position and shape of
:math:`\lambda(P, T)` across a parameter sweep -- rather than relying on the
absolute value of :math:`T_c` from a single calculation.

Choosing ``Nmat``, mixing, and ``IterationMax`` for FLEX convergence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Q. How should I set** ``Nmat``\ **,** ``mixing_scheme``\ **, and**
``IterationMax`` **for a well-converged FLEX run?**

``Nmat`` must be even and large enough for the Matsubara sums to resolve the
low-frequency structure at the target temperature (a good rule of thumb is
:math:`N_{\mathrm{mat}} \gtrsim 10/T`, see :ref:`flex_tips`). Anderson
mixing (``mixing_scheme = "anderson"``) together with a modest ``Mix`` and a
warm start typically reaches self-consistency in far fewer iterations than
plain linear mixing at low :math:`T`, so raise ``IterationMax`` only after
trying it. See :ref:`flex_params` for the full parameter table
(``IterationMax``, ``Mix``, ``EPS``, ``mixing_scheme``, ``anderson_depth``).

``calc_scheme = "auto"`` promotes hybridised/FLEX runs to ``general``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Q. My run got noticeably slower and its output arrays got bigger after
upgrading, even though I did not change my input. Why?**

Since version 2.0, ``calc_scheme = "auto"`` (the default) resolves to a
scheme that is *exact* for the declared input rather than always preferring
the cheaper ``reduced`` scheme. For RPA, ``reduced`` is kept only when every
declared interaction is density-only or the one-body Hamiltonian conserves
orbital flavour, so a hybridised multi-orbital model declaring
``CoulombInter``, ``Hund``, or ``Ising`` is promoted to ``general``. For
FLEX, *every* run declaring one of those interactions (or the aggregate
``Coulomb`` interaction) promotes to ``general`` unconditionally, hybridised
or not, because the self-consistent loop cannot exploit the RPA-only
exactness argument (see :ref:`flex_scope`). ``general`` uses a rank-6
``chiq`` and costs more (memory/solve cost grow from
:math:`O(N_d^2)/O(N_d^3)` to :math:`O(N_d^4)/O(N_d^6)` per
:math:`(\mathbf{q},\omega)`). To keep the pre-2.0 behaviour explicitly --
with a one-time warning that the result is then an approximation for that
input -- request ``calc_scheme = "reduced"``. See the
:ref:`calc_scheme parameter description <rpa_calc_scheme_auto>` for the full
migration note, including how to recover the old 4-axis density-channel
layout from the rank-6 array.
