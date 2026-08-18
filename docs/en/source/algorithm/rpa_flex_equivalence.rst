.. This page is GENERATED from tests/equivalence_cells.py by
   docs/tools/render_equivalence_table.py. Do not edit it by hand:
   re-run `python docs/tools/render_equivalence_table.py --write` and
   commit the result. A test in
   tests/test_rpa_flex_equivalence_table.py fails whenever this page
   and the registry disagree.

.. _rpa_flex_equivalence:

RPA and FLEX: support and equivalence
=====================================

H-wave computes the susceptibilities of a given model along two
routes: the RPA solver and the one-shot (first-iteration) FLEX
solver. Where the two routes describe the same approximation they
are expected to produce the same arrays, and this page records, for
a fixed set of small reference inputs, which configurations each
solver accepts and how closely the two agree where both run.

The comparison is one of IMPLEMENTATION, not of physics: a row
stating that the two solvers agree says that the two code paths
produce the same numbers on that input, not that either result is a
better description of the material. Likewise, a configuration
marked as rejected is one the solver refuses to run today; the
entry records the behaviour, not a statement about whether the
quantity is meaningful.

Every claim on this page is checked by
``tests/test_rpa_flex_equivalence_table.py``, which runs both
solvers on each input listed below and asserts the recorded
outcome. The page itself is generated from the same records, so it
cannot drift away from what the tests assert.

Status of these records
***********************

The tolerances below are **provisional**. They were measured on
development machines and have not yet been confirmed on the
project's continuous-integration runners, so the exact numbers
may still move. The per-tolerance notes name the machines and
the source revision each measurement came from.

How to read the table
*********************

The *RPA* and *FLEX* columns give each solver's outcome on the
input, and the *Relationship* column summarises the pair:

``EQUIV``
   Both solvers run and every compared array agrees to the recorded
   tolerance.

``DIVERGES(<issue>)``
   Both solvers run, but at least one compared array differs by more
   than the agreement policy allows. The difference is bounded and
   pinned by the tests, and the named issue tracks unifying the two
   implementations.

``AUTO-RESOLVES(<scheme>)``
   The input leaves ``calc_scheme = "auto"``, and both solvers
   resolve it to the named scheme. These rows check the resolution
   only; no susceptibility is computed.

The remaining labels
   No comparison is possible, because at least one solver does not
   run the input: ``BOTH-REJECT``, ``FLEX-REJECT·RPA-SUPPORTED``,
   ``RPA-REJECT·FLEX-SUPPORTED``, ``RPA-ONLY``, ``FLEX-ONLY``,
   ``RPA-REJECT`` and ``FLEX-REJECT``, naming the side that does
   not. The next section says what happens instead.

Support matrix
**************

.. list-table::
   :header-rows: 1
   :widths: 26 20 10 9 8 8 19

   * - Configuration
     - Interactions
     - Scheme
     - Spin mode
     - RPA
     - FLEX
     - Relationship
   * - ``auto.density.resolution``
     - ``CoulombInter`` (on-site)
     - ``auto`` -> ``reduced``
     - spin-free
     - supported
     - supported
     - ``AUTO-RESOLVES(reduced)``
   * - ``auto.exchange.resolution``
     - ``Exchange`` (on-site)
     - ``auto`` -> ``general``
     - spin-free
     - supported
     - supported
     - ``AUTO-RESOLVES(general)``
   * - ``auto.pairhop.resolution``
     - ``PairHop`` (on-site)
     - ``auto`` -> ``general``
     - spin-free
     - supported
     - supported
     - ``AUTO-RESOLVES(general)``
   * - ``chi0q_init.reuse``
     - ``CoulombInter`` (off-site)
     - ``general``
     - spin-free
     - supported
     - not applicable
     - ``RPA-ONLY``
   * - ``general.ring.offsite_coulombinter_interorb.flexreject``
     - ``CoulombInter`` (off-site)
     - ``general``
     - spin-free
     - supported
     - rejected
     - ``FLEX-REJECT·RPA-SUPPORTED``
   * - ``general.ring.offsite_coulombinter_sameorb.mu``
     - ``CoulombInter`` (off-site)
     - ``general``
     - spin-free
     - supported
     - supported
     - ``EQUIV``
   * - ``general.ring.offsite_coulombinter_sameorb.subshape``
     - ``CoulombInter`` (off-site)
     - ``general``
     - spin-free
     - supported
     - rejected
     - ``FLEX-REJECT·RPA-SUPPORTED``
   * - ``general.ring.offsite_coulombintra_literalkey.reject``
     - ``CoulombIntra`` (off-site)
     - ``general``
     - spin-free
     - rejected
     - rejected
     - ``BOTH-REJECT``
   * - ``general.ring.offsite_exchange.flexreject``
     - ``Exchange`` (off-site)
     - ``general``
     - spin-free
     - supported
     - rejected
     - ``FLEX-REJECT·RPA-SUPPORTED``
   * - ``general.ring.offsite_hund.flexreject``
     - ``Hund`` (off-site)
     - ``general``
     - spin-free
     - supported
     - rejected
     - ``FLEX-REJECT·RPA-SUPPORTED``
   * - ``general.ring.offsite_ising.flexreject``
     - ``Ising`` (off-site)
     - ``general``
     - spin-free
     - supported
     - rejected
     - ``FLEX-REJECT·RPA-SUPPORTED``
   * - ``general.ring.offsite_pairhop.flexreject``
     - ``PairHop`` (off-site)
     - ``general``
     - spin-free
     - supported
     - rejected
     - ``FLEX-REJECT·RPA-SUPPORTED``
   * - ``general.ring.onsite_coulombinter.coefftail.mu``
     - ``CoulombInter`` (on-site)
     - ``general``
     - spin-free
     - supported
     - supported
     - ``EQUIV``
   * - ``general.ring.onsite_coulombinter.conditioning.mu``
     - ``CoulombInter`` (off-site)
     - ``general``
     - spin-free
     - supported
     - supported
     - ``EQUIV``
   * - ``general.ring.onsite_coulombinter.fixedmu``
     - ``CoulombInter`` (on-site)
     - ``general``
     - spin-free
     - supported
     - supported
     - ``EQUIV``
   * - ``general.ring.onsite_coulombinter.subshape.mu``
     - ``CoulombInter`` (on-site)
     - ``general``
     - spin-free
     - supported
     - supported
     - ``EQUIV``
   * - ``general.ring.onsite_coulombintra.fixedmu``
     - ``CoulombIntra`` (on-site)
     - ``general``
     - spin-free
     - supported
     - supported
     - ``EQUIV``
   * - ``general.ring.onsite_exchange.fixedmu``
     - ``Exchange`` (on-site)
     - ``general``
     - spin-free
     - supported
     - supported
     - ``EQUIV``
   * - ``general.ring.onsite_full_kanamori.mu``
     - ``CoulombInter``, ``CoulombIntra``, ``Exchange``, ``Hund``, ``Ising``, ``PairHop`` (on-site)
     - ``general``
     - spin-free
     - supported
     - supported
     - ``EQUIV``
   * - ``general.ring.onsite_hund.fixedmu``
     - ``Hund`` (on-site)
     - ``general``
     - spin-free
     - supported
     - supported
     - ``EQUIV``
   * - ``general.ring.onsite_ising.fixedmu``
     - ``Ising`` (on-site)
     - ``general``
     - spin-free
     - supported
     - supported
     - ``EQUIV``
   * - ``general.ring.onsite_pairhop.fixedmu``
     - ``PairHop`` (on-site)
     - ``general``
     - spin-free
     - supported
     - supported
     - ``EQUIV``
   * - ``general.ring.onsite_pairlift.fixedmu``
     - ``PairLift`` (on-site)
     - ``general``
     - spin-free
     - supported
     - supported
     - ``EQUIV``
   * - ``general.ring.onsite_u_v_hund.mu``
     - ``CoulombInter``, ``CoulombIntra``, ``Hund`` (on-site)
     - ``general``
     - spin-free
     - supported
     - supported
     - ``EQUIV``
   * - ``reduced.ring.offsite_coulombinter.mu``
     - ``CoulombInter`` (off-site)
     - ``reduced``
     - spin-free
     - supported
     - supported
     - ``EQUIV``
   * - ``reduced.ring.onsite_coulombinter.spindiag.mu``
     - ``CoulombInter`` (on-site)
     - ``reduced``
     - spin-diag
     - supported
     - supported
     - ``EQUIV``
   * - ``reduced.ring.onsite_coulombinter.spinfree.mu``
     - ``CoulombInter`` (on-site)
     - ``reduced``
     - spin-free
     - supported
     - supported
     - ``EQUIV``
   * - ``reduced.ring.onsite_coulombintra.spindiag.mu``
     - ``CoulombIntra`` (on-site)
     - ``reduced``
     - spin-diag
     - supported
     - supported
     - ``EQUIV``
   * - ``reduced.ring.onsite_coulombintra.spinfree.mu``
     - ``CoulombIntra`` (on-site)
     - ``reduced``
     - spin-free
     - supported
     - supported
     - ``EQUIV``
   * - ``reduced.ring.onsite_exchange.reject``
     - ``Exchange`` (on-site)
     - ``reduced``
     - spin-free
     - rejected
     - rejected
     - ``BOTH-REJECT``
   * - ``reduced.ring.onsite_hund.spinfree.mu``
     - ``Hund`` (on-site)
     - ``reduced``
     - spin-free
     - supported
     - supported
     - ``EQUIV``
   * - ``reduced.ring.onsite_ising.spinfree.mu``
     - ``Ising`` (on-site)
     - ``reduced``
     - spin-free
     - supported
     - supported
     - ``EQUIV``
   * - ``reduced.ring.onsite_pairhop.reject``
     - ``PairHop`` (on-site)
     - ``reduced``
     - spin-free
     - rejected
     - rejected
     - ``BOTH-REJECT``
   * - ``reduced.ring.onsite_pairlift.spinfree.mu``
     - ``PairLift`` (on-site)
     - ``reduced``
     - spin-free
     - supported
     - supported
     - ``EQUIV``
   * - ``ringladder.general.onsite_coulombintra``
     - ``CoulombIntra`` (on-site)
     - ``general``
     - spin-free
     - supported
     - rejected
     - ``FLEX-REJECT·RPA-SUPPORTED``
   * - ``so.general.construction.reject``
     - ``CoulombIntra`` (on-site)
     - ``general``
     - spinful
     - supported
     - rejected
     - ``FLEX-REJECT·RPA-SUPPORTED``
   * - ``so.reduced.construction.reject``
     - ``CoulombIntra`` (on-site)
     - ``reduced``
     - spinful
     - supported
     - rejected
     - ``FLEX-REJECT·RPA-SUPPORTED``

Configurations one solver does not run
**************************************

For these inputs the two solvers cannot be compared. The entries
record what the unsupported side does instead.

``chi0q_init.reuse``
   * FLEX: nothing to compare -- accepted; no corresponding option semantics (flex.py:408) -- FLEX starts every SCF loop from zero self-energy and recomputes chi0q from the dressed Green's function each iteration (flex.py:446-451's own docstring), so a chi0q_init entry loaded by the inherited RPA read_init is never consumed.

``general.ring.offsite_coulombinter_interorb.flexreject``
   * FLEX refuses the input during the solve and raises ``ValueError`` whose message contains ``interaction 'CoulombInter'``.

``general.ring.offsite_coulombinter_sameorb.subshape``
   * FLEX refuses the input during the solve and raises ``ValueError`` whose message contains ``with sublattice folding``.

``general.ring.offsite_coulombintra_literalkey.reject``
   * RPA refuses the input at construction and raises ``ValueError`` whose message contains ``the documented operator is on-site and same-orbital``.
   * FLEX refuses the input at construction and raises ``ValueError`` whose message contains ``the documented operator is on-site and same-orbital``.

``general.ring.offsite_exchange.flexreject``
   * FLEX refuses the input during the solve and raises ``ValueError`` whose message contains ``interaction 'Exchange'``.

``general.ring.offsite_hund.flexreject``
   * FLEX refuses the input during the solve and raises ``ValueError`` whose message contains ``interaction 'Hund'``.

``general.ring.offsite_ising.flexreject``
   * FLEX refuses the input during the solve and raises ``ValueError`` whose message contains ``interaction 'Ising'``.

``general.ring.offsite_pairhop.flexreject``
   * FLEX refuses the input during the solve and raises ``ValueError`` whose message contains ``interaction 'PairHop'``.

``reduced.ring.onsite_exchange.reject``
   * RPA refuses the input at construction and raises ``ValueError`` whose message contains ``reduced``.
   * FLEX refuses the input at construction and raises ``ValueError`` whose message contains ``reduced``.

``reduced.ring.onsite_pairhop.reject``
   * RPA refuses the input at construction and raises ``ValueError`` whose message contains ``reduced``.
   * FLEX refuses the input at construction and raises ``ValueError`` whose message contains ``reduced``.

``ringladder.general.onsite_coulombintra``
   * FLEX refuses the input at construction and raises ``ValueError`` whose message contains ``ring+ladder``.

``so.general.construction.reject``
   * FLEX refuses the input at construction and raises ``ValueError`` whose message contains ``enable_spin_orbital``.

``so.reduced.construction.reject``
   * FLEX refuses the input at construction and raises ``ValueError`` whose message contains ``enable_spin_orbital``.

Recorded agreement
******************

Arrays are compared elementwise on complex values: the two results
count as agreeing when ``|a - b| <= atol`` holds for every element,
with no relative term. ``chi0q`` is the bare susceptibility and
``chiq`` the interacting one; the comparator names the mapping
applied before the comparison, which reconciles the two solvers'
storage layouts (``identity`` when both already store the same
shape).

Each entry also carries the measurement note recorded with the
tolerance -- the machines, the source revision and the observed
differences it was derived from. Those notes are maintenance
records rather than user guidance.

``general.ring.offsite_coulombinter_sameorb.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-13``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 4.885e-15, chita 4.899e-15 -- candidate atol 1.0e-13 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-12``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 2.021e-14, chita 2.026e-14 -- candidate atol 1.0e-12 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).

``general.ring.onsite_coulombinter.coefftail.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.110e-16, chita 1.110e-16 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.388e-16, chita 1.388e-16 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).

``general.ring.onsite_coulombinter.conditioning.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-11``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 6.432e-13, chita 6.432e-13 -- candidate atol 1.0e-11 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-10``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 3.307e-12, chita 3.373e-12 -- candidate atol 1.0e-10 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).

``general.ring.onsite_coulombinter.fixedmu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.110e-16, chita 9.714e-17 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.249e-16, chita 1.110e-16 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``general.ring.onsite_coulombinter.subshape.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.665e-16, chita 9.715e-17 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.943e-16, chita 1.110e-16 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).

``general.ring.onsite_coulombintra.fixedmu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.110e-16, chita 9.714e-17 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 2.498e-16, chita 1.943e-16 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``general.ring.onsite_exchange.fixedmu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.110e-16, chita 9.714e-17 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.110e-16, chita 9.714e-17 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``general.ring.onsite_full_kanamori.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.110e-16, chita 9.714e-17 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 2.498e-16, chita 1.943e-16 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).

``general.ring.onsite_hund.fixedmu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.110e-16, chita 9.714e-17 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.110e-16, chita 1.110e-16 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``general.ring.onsite_ising.fixedmu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.110e-16, chita 9.714e-17 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.249e-16, chita 1.110e-16 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``general.ring.onsite_pairhop.fixedmu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.110e-16, chita 9.714e-17 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.110e-16, chita 9.715e-17 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``general.ring.onsite_pairlift.fixedmu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.110e-16, chita 9.714e-17 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.110e-16, chita 9.714e-17 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``general.ring.onsite_u_v_hund.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.110e-16, chita 9.714e-17 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 3.053e-16, chita 2.776e-16 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).

``reduced.ring.offsite_coulombinter.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-13``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 4.885e-15, chita 4.899e-15 -- candidate atol 1.0e-13 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).
   * ``chiq`` (comparator ``reduced_blocks``): agree to ``1.0e-12``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 3.553e-14, chita 3.561e-14 -- candidate atol 1.0e-12 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).

``reduced.ring.onsite_coulombinter.spindiag.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.249e-16, chita 1.249e-16 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).
   * ``chiq`` (comparator ``reduced_blocks``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.527e-16, chita 1.388e-16 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).

``reduced.ring.onsite_coulombinter.spinfree.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.110e-16, chita 9.714e-17 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).
   * ``chiq`` (comparator ``reduced_blocks``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.249e-16, chita 1.249e-16 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).

``reduced.ring.onsite_coulombintra.spindiag.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.249e-16, chita 1.249e-16 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).
   * ``chiq`` (comparator ``reduced_blocks``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 4.164e-16, chita 4.719e-16 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).

``reduced.ring.onsite_coulombintra.spinfree.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.110e-16, chita 9.714e-17 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).
   * ``chiq`` (comparator ``reduced_blocks``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 3.886e-16, chita 3.053e-16 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).

``reduced.ring.onsite_hund.spinfree.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.110e-16, chita 9.714e-17 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).
   * ``chiq`` (comparator ``reduced_blocks``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.110e-16, chita 1.110e-16 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).

``reduced.ring.onsite_ising.spinfree.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.110e-16, chita 9.714e-17 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).
   * ``chiq`` (comparator ``reduced_blocks``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.249e-16, chita 1.249e-16 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).

``reduced.ring.onsite_pairlift.spinfree.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.110e-16, chita 9.714e-17 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).
   * ``chiq`` (comparator ``reduced_blocks``): agree to ``1.0e-14``.
     Measured: candidate (Task 7; local macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1; chita Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1; commit b922a1c13b85): max|diff| local 1.110e-16, chita 9.714e-17 -- candidate atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-10).

Related tests
*************

A maintenance index of tests elsewhere in the suite that cover the
same ground, each with the note recorded alongside it. They are
listed for orientation only: the claims on this page are
established by ``tests/test_rpa_flex_equivalence_table.py`` alone.

``general.ring.offsite_coulombinter_interorb.flexreject``
   * FLEX -- ``tests.test_flex_offsite_general::TestOffsiteGeneralFLEX::test_rejected_offsite_classes``: FLEX rejects off-site CoulombInter with a != b (inter-orbital) under calc_scheme='general'.

``general.ring.offsite_coulombinter_sameorb.mu``
   * FLEX -- ``tests.test_flex_offsite_general::TestOffsiteGeneralFLEX::test_one_orbital_offsite_v_matches_the_rpa_ring_exactly``: The same fixture (tests/rpa/input's coulombinter.dat, filling=0.75) element-complete equal between RPA and FLEX for off-site same-orbital CoulombInter.

``general.ring.offsite_exchange.flexreject``
   * FLEX -- ``tests.test_flex_offsite_general::TestOffsiteGeneralFLEX::test_rejected_offsite_classes``: FLEX rejects off-site Exchange under calc_scheme='general'.

``general.ring.offsite_hund.flexreject``
   * FLEX -- ``tests.test_flex_offsite_general::TestOffsiteGeneralFLEX::test_rejected_offsite_classes``: FLEX rejects off-site Hund under calc_scheme='general'.

``general.ring.offsite_ising.flexreject``
   * FLEX -- ``tests.test_flex_offsite_general::TestOffsiteGeneralFLEX::test_rejected_offsite_classes``: FLEX rejects off-site Ising under calc_scheme='general'.

``general.ring.offsite_pairhop.flexreject``
   * FLEX -- ``tests.test_flex_offsite_general::TestOffsiteGeneralFLEX::test_rejected_offsite_classes``: tests/test_flex_offsite_general.py's general rejection-classes coverage (this exact PairHop case is not one of its enumerated subTest rows -- the module's docstring names off-site PairHop as rejected for the same non-local-pair reason as Exchange; linked for the surrounding module context, not a per-type pin).

``reduced.ring.offsite_coulombinter.mu``
   * RPA -- ``tests.test_rpa_flex_oneshot_equivalence::TestReducedOneShot::test_matrix_cells``: The oneshot suite's TestReducedOneShot case (a) -- this exact fixture/filling -- pins the chiq uu-ud/uu+ud block equivalence under calc_scheme='reduced'.

``reduced.ring.onsite_coulombinter.spindiag.mu``
   * RPA -- ``tests.test_rpa_flex_equivalence_table::TestReducedSpinDiagAcceptance::test_reduced_rpa_also_solves_with_extern``: RPA accepts the E4 Extern/Zeeman spin-diag route under calc_scheme='reduced' (Task 4).
   * FLEX -- ``tests.test_rpa_flex_equivalence_table::TestReducedSpinDiagAcceptance::test_reduced_flex_solves_with_extern_coulombinter``: FLEX accepts the E4 Extern/Zeeman spin-diag route under calc_scheme='reduced' with CoulombInter (Task 4 deliverable 1).

``reduced.ring.onsite_coulombintra.spindiag.mu``
   * RPA -- ``tests.test_rpa_flex_equivalence_table::TestReducedSpinDiagAcceptance::test_reduced_rpa_also_solves_with_extern``: RPA accepts the E4 Extern/Zeeman spin-diag route under calc_scheme='reduced' (Task 4).
   * FLEX -- ``tests.test_rpa_flex_equivalence_table::TestReducedSpinDiagAcceptance::test_reduced_flex_solves_with_extern_coulombintra``: FLEX accepts the E4 Extern/Zeeman spin-diag route under calc_scheme='reduced' -- the frozen both-accept expectation Appendix A records for this row (Task 4 deliverable 1).

``ringladder.general.onsite_coulombintra``
   * FLEX -- ``tests.test_spinful_transverse_ed::TestSpinfulVertexExchangeOptOutRingLadderRejection::test_ring_ladder_rejects_the_opt_out_combination``: The transverse ring+ladder pipeline's own rejection contract (a companion instance to FLEX's blanket calc_type='ring+ladder' rejection).
   * FLEX -- ``tests.test_bond_transverse_w::TestGateW1OnsiteReduction::test_collapsed_output_matches_plain_chiq_pm_at_omega_zero``: RPA's ring+ladder transverse channel resolves correctly on the collapsed (non-bond) path -- the capability FLEX does not implement.

Scope
*****

The table certifies the inputs it lists and nothing beyond them.
All of them are deliberately small, so that the whole set runs as
part of the ordinary test suite:

* Lattices ``2x1x1``, ``4x4x1`` with sublattices ``1x1x1``, ``2x1x1``.
* Temperatures ``0.2``, ``2.0``.
* Matsubara frequency counts ``8``, ``32``, ``256``.
* Diagram selections ``ring``, ``ring+ladder``.
* Calculation schemes ``general``, ``reduced``.
* Both a fixed chemical potential and a fixed electron count.

Agreement measured on inputs of this size does not by itself
establish agreement on a production-sized model, and the tolerances
are the measured round-off levels of these inputs, not general
accuracy guarantees.

Accelerated backends
--------------------

The table makes no statement about the GPU backends: every input is
run on the CPU backend of both solvers. That the GPU backend
reproduces the CPU one is a separate, per-solver property, covered
by ``tests/test_gpu_backend.py``, ``tests/test_rpa_gpu.py``,
``tests/test_flex_gpu.py`` and ``tests/test_sc_gpu.py``.

FLEX-only numerical representations
-----------------------------------

The intermediate-representation (IR) route has no counterpart on
the RPA side, so it cannot appear as a row here: comparing the IR
route with the dense one is a comparison of FLEX with itself. That
comparison is covered by ``tests/test_flex_ir.py`` and
``tests/test_flex_ir_general.py``.
