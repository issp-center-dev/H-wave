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

:math:`\lambda = 1` identifies :math:`T_c` when :math:`\lambda` is the
eigenvalue of a single physical branch tracked continuously across
temperature; :math:`\lambda < 1` at a given :math:`T` only tells you that the
normal state is stable at that particular temperature, not how close
:math:`T` is to :math:`T_c`. Finding the crossing therefore requires running
at lower temperatures. ``hwave_tsweep`` automates a descending temperature
ladder, chaining each rung's final saved self-energy (and, for the dynamic
solver, its leading eigenvector) into the next -- see :ref:`sc_tsweep`, and
check the summary's ``flex_converged`` column, since an unconverged rung's
self-energy is what seeds the next rung regardless. Because a
FLEX+Eliashberg :math:`T_c` is only semi-quantitative (see
:ref:`faq_flex_tc`), comparing the position and shape of :math:`\lambda(T)`
or :math:`\lambda(P)` across a parameter sweep is more robust than reading
off the absolute value where a tracked branch crosses 1.

FLEX does not converge at low temperature (cold-start failure)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Q. My FLEX run at low** :math:`T` **fails to converge, even though a similar
run at higher** :math:`T` **works. What should I do?**

A cold start (:math:`\Sigma = 0`) directly at low :math:`T` often fails to
converge, because the SCF fixed point moves further from :math:`\Sigma = 0`
as correlations strengthen. Instead, start from a higher temperature where
convergence is easy and step down, feeding each temperature's final
self-energy into the next run as a warm start (``[file.input] sigma_init``,
see :ref:`flex_sigma_init`). ``hwave_tsweep`` automates this whole ladder and
records each rung's convergence status, so you can check whether a given
warm start actually came from a converged rung (:ref:`sc_tsweep`).

The leading eigenvalue or gap symmetry jumps between temperatures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Q. Stepping through a temperature sweep, the leading eigenvalue's gap
symmetry suddenly switches to a different one. Is my calculation broken?**

Not necessarily. The dynamic (frequency-dependent) Eliashberg kernel is
non-Hermitian and can have *exceptional points* where two real eigenvalues
collide and split into a complex-conjugate pair, so the "leading" branch can
jump discontinuously between neighbouring temperatures even when the FLEX
self-energy varies smoothly. Setting ``[eliashberg] seed_eigenvector`` to a
neighbouring run's ``gap_dynamic.npz`` helps track one physical branch (e.g.
the d-wave mode) continuously across the sweep: among the eigenpairs actually
*returned* by the eigensolver (governed by ``num_eigenvalues`` and, on the
shift-invert path, ``sigma_shift``), it selects the one whose eigenvector
best overlaps the seed, instead of the algebraically largest one. If the
physical branch is not among the returned eigenpairs at all, increasing
``num_eigenvalues`` or moving ``sigma_shift`` near it is what brings it back
into range. See :ref:`sc_seed_eigenvector` for the full mechanism.

How do I identify the gap symmetry I obtained? (the ``match`` column)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Q. I set** ``init_gap`` **to seed a particular symmetry, but I am not sure
the solver actually converged to it. How do I check?**

``init_gap`` (static solver) only seeds the self-consistent iteration; the
reported result is the leading eigenpair of the requested channel *parity*
(even for singlet, odd for triplet), which need not be the symmetry you
seeded, and the eigenvalue analysis also reports opposite-parity solutions
that are mathematically valid but unphysical for that channel. The trailing
``match`` column of ``eigenvalue.dat`` is a dominance test, not an exact
one: it is ``1`` when projecting the eigenvector onto the requested
channel-parity sector retains at least 90% of its norm (physical), and ``0``
otherwise (spurious). If no computed eigenpair reaches that threshold, the
solver warns and falls back to reporting the unmatched leading eigenpair --
in that case also try increasing ``num_eigenvalues`` or double-check
``pairing_type``. Inspect the reported gap function itself (``gap.dat`` /
the figures in the tutorial) to identify which symmetry was actually
obtained. See :ref:`sc_eigenvalue_dat` for the file format and the parity
note above it for why opposite-parity solutions appear at all.

:math:`\lambda(T)` is non-monotonic -- is that physical?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Q. For a fixed pairing channel,** :math:`\lambda` **dips and then rises
again as I lower the temperature. Is that a real effect?**

Rule out the usual non-physical causes before concluding it is. For a fixed,
continuously-tracked pairing branch, :math:`\lambda` is expected to increase
as :math:`T \to 0`, so an unexplained dip is worth checking against: a
cold-started or under-converged FLEX self-energy at that rung rather than a
warm-started, converged one (check the ``flex_converged`` column in the
``hwave_tsweep`` summary, :ref:`flex_sigma_init`); a branch switch at an
exceptional point, where the "leading" eigenpair silently changes symmetry
between rungs (:ref:`sc_seed_eigenvector`); and, on the IR path, an
under-resolved ``ir_wmax`` or FLEX ``Nmat`` (:ref:`sc_dynamic_ir_en`). Step
through the full temperature ladder (:ref:`sc_tsweep`), checking these one by
one, and use ``seed_eigenvector`` to confirm you are tracking one physical
branch, before treating a dip as a genuine non-monotonic effect.

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
(check against a larger value or the automatic estimate); making it too
large can instead make the delta(:math:`\tau`)-derived constant that the IR
loader has to isolate ill-conditioned. If that happens, the remedies are to
use the automatic ``ir_wmax`` estimate (or a value near
``3*(bandwidth + max interaction)``), raise the FLEX ``Nmat``, set
``ir_keep_static_chi = true`` to retain a genuinely static component instead
of isolating it, or fall back to ``matsubara_basis = "uniform"``. See
:ref:`sc_dynamic_ir_en` for the automatic estimate, the isolation
diagnostic, and measured eigenvalue sensitivity.

GPU out-of-memory at large ``Nmat``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Q. My dynamic Eliashberg run with** ``gpu = true`` **aborts with an
out-of-memory error at large** ``Nmat``\ **. What can I do?**

The two large resident device tensors (the pair bubble and the pairing
vertex) scale as :math:`N_{\mathrm{orb}}^4 N_k N_{\mathrm{mat}}`, so GPU
memory use grows directly with ``Nmat`` (and with the k-mesh size and orbital
count). Reduce ``Nmat`` or ``CellShape``, or run fewer concurrent solves on
the same device. Before the large allocation, the solver always logs the
resident tensor size; only when that estimate exceeds the device's free
memory does it additionally log a warning naming the required and free
amounts (advisory only -- CuPy itself still raises a clear out-of-memory
error on the actual allocation). Watch for that warning rather than
assuming a silent fit. See :ref:`sc_dynamic_gpu_en` for the memory formula
and a reference speed/memory point.

Two-orbital ``chi_convention`` (``kuroki`` vs. ``myo``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Q. My two-orbital FLEX+Eliashberg results look wrong even though the file
shapes match what I expect. What should I check?**

The FLEX susceptibility files carry a ``chi_convention`` tag (``"kuroki"`` for
the reduced scheme, ``"myo"`` for the general full-vertex scheme) that the
Eliashberg loader uses to interpret the stored orbital layout. The loader
raises rather than guesses whenever it can tell something is wrong: a
``chi_s``/``chi_c`` pair tagged with two different conventions, an
unrecognized tag, or -- for :math:`N_{\mathrm{orb}} \neq 2`, where the array
shape unambiguously implies one convention -- a tag that disagrees with the
shape, all raise an explicit error rather than mis-reading the file. The one
case that cannot be caught this way is :math:`N_{\mathrm{orb}} = 2`
specifically: there the reduced spin-orbital dimension and the general
orbital-pair dimension coincide (both equal 4), so the array *shape alone*
cannot disambiguate them. If such a file is internally consistent but simply
carries the wrong tag for what it actually is (e.g. a kuroki-layout file
mistagged ``"myo"``), the loader has no way to detect that, and it silently
builds an incorrect pairing vertex. See :ref:`sc_dynamic_frequency` for the
full explanation and the historical mislabeling this tag now guards against.


FLEX approximation
=================================================================

Validity range of FLEX
^^^^^^^^^^^^^^^^^^^^^^^^

**Q. What are the limitations of FLEX? When should I not trust it?**

FLEX is a conserving, weak-to-intermediate-coupling approximation built from
the RPA particle-hole bubble and ladder series; it does not include the
Aslamazov-Larkin/Maki-Thompson vertex corrections, and it cannot describe a
Mott insulating state. The paramagnetic FLEX solved here also cannot enter a
magnetically or charge-ordered phase, so a large susceptibility or pairing
eigenvalue obtained near such an instability should be treated as unreliable
rather than as a confirmed prediction, and derived transition temperatures
are semi-quantitative at best (see :ref:`faq_flex_tc`). See :ref:`flex_scope`
for the precise scope, including the density-density (``reduced``) vs.
full-vertex (``general``) distinction.

My FLEX run diverges or will not converge near an instability
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Q. As I approach what looks like a magnetic or charge instability, FLEX
stops converging (or the susceptibility blows up). What is happening?**

As the system approaches a magnetic instability, the Stoner factor (the
interaction times :math:`\chi_0`, RPA-summed) approaches 1, so the RPA-type
denominator (:math:`1 -` the Stoner factor) approaches **zero** -- not the
other way around -- and :math:`\chi` grows large, making the self-consistency
loop harder to converge. Approach the region gradually with a warm-start
temperature ladder (:ref:`flex_sigma_init`, :ref:`sc_tsweep`) and gentler
mixing -- see :ref:`flex_tips` and :ref:`flex_params` for concrete
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
