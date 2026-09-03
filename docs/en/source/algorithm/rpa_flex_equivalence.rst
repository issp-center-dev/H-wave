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

The tolerances below are **confirmed**: they were reproduced on the
project's continuous-integration runners at source revision
``8144bf3f9e9539bad4759a2fbd1b24f52f7bef33`` (workflow run
``32204319966 attempt 1``), and hold as of that revision. The notes
below record where each residual was originally measured.

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
     - ``auto`` -> ``general``
     - spin-free
     - supported
     - supported
     - ``AUTO-RESOLVES(general)``
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
   * - ``general.ring.offsite_coulombinter.conditioning.mu``
     - ``CoulombInter`` (off-site)
     - ``general``
     - spin-free
     - supported
     - supported
     - ``EQUIV``
   * - ``general.ring.offsite_coulombinter_interorb.mu``
     - ``CoulombInter`` (off-site)
     - ``general``
     - spin-free
     - supported
     - supported
     - ``EQUIV``
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
     - supported
     - ``EQUIV``
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
   * - ``general.ring.offsite_hund.mu``
     - ``Hund`` (off-site)
     - ``general``
     - spin-free
     - supported
     - supported
     - ``EQUIV``
   * - ``general.ring.offsite_ising.mu``
     - ``Ising`` (off-site)
     - ``general``
     - spin-free
     - supported
     - supported
     - ``EQUIV``
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
   * FLEX: nothing to compare -- accepted; no corresponding option semantics (flex.py:408) -- FLEX starts every SCF loop from zero self-energy and recomputes chi0q from the dressed Green's function each iteration (flex.py:446-451's own docstring), so a chi0q_init entry loaded by the inherited RPA read_init is never consumed. Measured: supplying a chi0q_init DELIBERATELY PERTURBED to twice the value FLEX itself computes leaves every one of FLEX's output arrays bitwise unchanged.

``general.ring.offsite_coulombintra_literalkey.reject``
   * RPA refuses the input at construction and raises ``ValueError`` whose message contains ``the documented operator is on-site and same-orbital``.
   * FLEX refuses the input at construction and raises ``ValueError`` whose message contains ``the documented operator is on-site and same-orbital``.

``general.ring.offsite_exchange.flexreject``
   * FLEX refuses the input during the solve and raises ``ValueError`` whose message contains ``off-site 'Exchange' entry``.

``general.ring.offsite_pairhop.flexreject``
   * FLEX refuses the input during the solve and raises ``ValueError`` whose message contains ``off-site 'PairHop' entry``.

``reduced.ring.onsite_exchange.reject``
   * RPA refuses the input at construction and raises ``ValueError`` whose message contains ``has no density-diagonal content``.
   * FLEX refuses the input at construction and raises ``ValueError`` whose message contains ``has no density-diagonal content``.

``reduced.ring.onsite_pairhop.reject``
   * RPA refuses the input at construction and raises ``ValueError`` whose message contains ``has no density-diagonal content``.
   * FLEX refuses the input at construction and raises ``ValueError`` whose message contains ``has no density-diagonal content``.

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

``general.ring.offsite_coulombinter.conditioning.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: RECALIBRATED for the #160 mu/Green-seam fix (calibration log Event 5, cell general.ring.offsite_coulombinter.conditioning.mu): measured on the macOS arm64 development machine (Python 3.13.13, numpy 2.4.6, scipy 1.17.1) max|diff| 5.551275e-17; MAX over the four gating CI runners (ubuntu-latest x Python 3.9-3.12, workflow run 33278664447 attempt 1) max|diff| 5.551320e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-14). Supersedes this cell's original two-development-machine literal, which dated from before the #160 fix and no longer reflected this cell's actual (much smaller) post-fix residual.
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-12``.
     Provenance: RECALIBRATED for the #160 mu/Green-seam fix (calibration log Event 5, cell general.ring.offsite_coulombinter.conditioning.mu): measured on the macOS arm64 development machine (Python 3.13.13, numpy 2.4.6, scipy 1.17.1) max|diff| 7.294815e-14; MAX over the four gating CI runners (ubuntu-latest x Python 3.9-3.12, workflow run 33278664447 attempt 1) max|diff| 6.843859e-14. Resulting atol 1.0e-12 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12). Supersedes this cell's original two-development-machine literal, which dated from before the #160 fix and no longer reflected this cell's actual (much smaller) post-fix residual.

``general.ring.offsite_coulombinter_interorb.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: NEW comparison cell (#181 Tier 1, calibration log Event 7, cell general.ring.offsite_coulombinter_interorb.mu; FLEX-REJECT.RPA-SUPPORTED before): measured on the macOS arm64 development machine (Python 3.13.13, numpy 2.4.6, scipy 1.17.1) max|diff| 9.714461e-17; MAX over the three gating CI runners (ubuntu-latest x Python 3.10-3.12, workflow run 33706826878 attempt 1 (PROVISIONAL, measured at 76e09354; refrozen in the calibration commit)) max|diff| 9.714755e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-14).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Provenance: NEW comparison cell (#181 Tier 1, calibration log Event 7, cell general.ring.offsite_coulombinter_interorb.mu; FLEX-REJECT.RPA-SUPPORTED before): measured on the macOS arm64 development machine (Python 3.13.13, numpy 2.4.6, scipy 1.17.1) max|diff| 1.110224e-16; MAX over the three gating CI runners (ubuntu-latest x Python 3.10-3.12, workflow run 33706826878 attempt 1 (PROVISIONAL, measured at 76e09354; refrozen in the calibration commit)) max|diff| 1.110249e-16. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``general.ring.offsite_coulombinter_sameorb.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: RECALIBRATED for the #160 mu/Green-seam fix (calibration log Event 5, cell general.ring.offsite_coulombinter_sameorb.mu): measured on the macOS arm64 development machine (Python 3.13.13, numpy 2.4.6, scipy 1.17.1) max|diff| 5.551921e-17; MAX over the four gating CI runners (ubuntu-latest x Python 3.9-3.12, workflow run 33278664447 attempt 1) max|diff| 4.167915e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-14). Supersedes this cell's original two-development-machine literal, which dated from before the #160 fix and no longer reflected this cell's actual (much smaller) post-fix residual.
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Provenance: RECALIBRATED for the #160 mu/Green-seam fix (calibration log Event 5, cell general.ring.offsite_coulombinter_sameorb.mu): measured on the macOS arm64 development machine (Python 3.13.13, numpy 2.4.6, scipy 1.17.1) max|diff| 2.220865e-16; MAX over the four gating CI runners (ubuntu-latest x Python 3.9-3.12, workflow run 33278664447 attempt 1) max|diff| 1.944572e-16. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12). Supersedes this cell's original two-development-machine literal, which dated from before the #160 fix and no longer reflected this cell's actual (much smaller) post-fix residual.

``general.ring.offsite_coulombinter_sameorb.subshape``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: NEW comparison cell (#181 Tier 1, calibration log Event 7, cell general.ring.offsite_coulombinter_sameorb.subshape; FLEX-REJECT.RPA-SUPPORTED before): measured on the macOS arm64 development machine (Python 3.13.13, numpy 2.4.6, scipy 1.17.1) max|diff| 2.776646e-17; MAX over the three gating CI runners (ubuntu-latest x Python 3.10-3.12, workflow run 33706826878 attempt 1 (PROVISIONAL, measured at 76e09354; refrozen in the calibration commit)) max|diff| 4.163442e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Provenance: NEW comparison cell (#181 Tier 1, calibration log Event 7, cell general.ring.offsite_coulombinter_sameorb.subshape; FLEX-REJECT.RPA-SUPPORTED before): measured on the macOS arm64 development machine (Python 3.13.13, numpy 2.4.6, scipy 1.17.1) max|diff| 5.645067e-16; MAX over the three gating CI runners (ubuntu-latest x Python 3.10-3.12, workflow run 33706826878 attempt 1 (PROVISIONAL, measured at 76e09354; refrozen in the calibration commit)) max|diff| 9.267834e-16. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``general.ring.offsite_hund.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: NEW comparison cell (#181 Tier 1, calibration log Event 7, cell general.ring.offsite_hund.mu; FLEX-REJECT.RPA-SUPPORTED before): measured on the macOS arm64 development machine (Python 3.13.13, numpy 2.4.6, scipy 1.17.1) max|diff| 9.714461e-17; MAX over the three gating CI runners (ubuntu-latest x Python 3.10-3.12, workflow run 33706826878 attempt 1 (PROVISIONAL, measured at 76e09354; refrozen in the calibration commit)) max|diff| 9.714755e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-14).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Provenance: NEW comparison cell (#181 Tier 1, calibration log Event 7, cell general.ring.offsite_hund.mu; FLEX-REJECT.RPA-SUPPORTED before): measured on the macOS arm64 development machine (Python 3.13.13, numpy 2.4.6, scipy 1.17.1) max|diff| 9.714462e-17; MAX over the three gating CI runners (ubuntu-latest x Python 3.10-3.12, workflow run 33706826878 attempt 1 (PROVISIONAL, measured at 76e09354; refrozen in the calibration commit)) max|diff| 9.714755e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``general.ring.offsite_ising.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: NEW comparison cell (#181 Tier 1, calibration log Event 7, cell general.ring.offsite_ising.mu; FLEX-REJECT.RPA-SUPPORTED before): measured on the macOS arm64 development machine (Python 3.13.13, numpy 2.4.6, scipy 1.17.1) max|diff| 9.714461e-17; MAX over the three gating CI runners (ubuntu-latest x Python 3.10-3.12, workflow run 33706826878 attempt 1 (PROVISIONAL, measured at 76e09354; refrozen in the calibration commit)) max|diff| 9.714755e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-14).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Provenance: NEW comparison cell (#181 Tier 1, calibration log Event 7, cell general.ring.offsite_ising.mu; FLEX-REJECT.RPA-SUPPORTED before): measured on the macOS arm64 development machine (Python 3.13.13, numpy 2.4.6, scipy 1.17.1) max|diff| 1.110224e-16; MAX over the three gating CI runners (ubuntu-latest x Python 3.10-3.12, workflow run 33706826878 attempt 1 (PROVISIONAL, measured at 76e09354; refrozen in the calibration commit)) max|diff| 1.110249e-16. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``general.ring.onsite_coulombinter.coefftail.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.110e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 1.110e-16. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-14).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.388e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 1.388e-16. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``general.ring.onsite_coulombinter.fixedmu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.110e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 9.714e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.249e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 1.110e-16. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``general.ring.onsite_coulombinter.subshape.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.665e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 9.715e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-14).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.943e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 1.110e-16. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``general.ring.onsite_coulombintra.fixedmu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.110e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 9.714e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 2.498e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 1.943e-16. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``general.ring.onsite_exchange.fixedmu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.110e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 9.714e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.110e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 9.714e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``general.ring.onsite_full_kanamori.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.110e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 9.714e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-14).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 2.498e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 1.943e-16. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``general.ring.onsite_hund.fixedmu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.110e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 9.714e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.110e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 1.110e-16. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``general.ring.onsite_ising.fixedmu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.110e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 9.714e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.249e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 1.110e-16. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``general.ring.onsite_pairhop.fixedmu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.110e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 9.714e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.110e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 9.715e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``general.ring.onsite_pairlift.fixedmu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.110e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 9.714e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.110e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 9.714e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``general.ring.onsite_u_v_hund.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.110e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 9.714e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-14).
   * ``chiq`` (comparator ``general_from_flex_channels``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 3.053e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 2.776e-16. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``reduced.ring.offsite_coulombinter.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: RECALIBRATED for the #160 mu/Green-seam fix (calibration log Event 5, cell reduced.ring.offsite_coulombinter.mu): measured on the macOS arm64 development machine (Python 3.13.13, numpy 2.4.6, scipy 1.17.1) max|diff| 5.551921e-17; MAX over the four gating CI runners (ubuntu-latest x Python 3.9-3.12, workflow run 33278664447 attempt 1) max|diff| 4.167915e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-14). Supersedes this cell's original two-development-machine literal, which dated from before the #160 fix and no longer reflected this cell's actual (much smaller) post-fix residual.
   * ``chiq`` (comparator ``reduced_blocks``): agree to ``1.0e-14``.
     Provenance: RECALIBRATED for the #160 mu/Green-seam fix (calibration log Event 5, cell reduced.ring.offsite_coulombinter.mu): measured on the macOS arm64 development machine (Python 3.13.13, numpy 2.4.6, scipy 1.17.1) max|diff| 3.609021e-16; MAX over the four gating CI runners (ubuntu-latest x Python 3.9-3.12, workflow run 33278664447 attempt 1) max|diff| 3.333702e-16. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12). Supersedes this cell's original two-development-machine literal, which dated from before the #160 fix and no longer reflected this cell's actual (much smaller) post-fix residual.

``reduced.ring.onsite_coulombinter.spindiag.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.249e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 1.249e-16. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-14).
   * ``chiq`` (comparator ``reduced_blocks``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.527e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 1.388e-16. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``reduced.ring.onsite_coulombinter.spinfree.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.110e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 9.714e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-14).
   * ``chiq`` (comparator ``reduced_blocks``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.249e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 1.249e-16. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``reduced.ring.onsite_coulombintra.spindiag.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.249e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 1.249e-16. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-14).
   * ``chiq`` (comparator ``reduced_blocks``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 4.164e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 4.719e-16. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``reduced.ring.onsite_coulombintra.spinfree.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.110e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 9.714e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-14).
   * ``chiq`` (comparator ``reduced_blocks``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 3.886e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 3.053e-16. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``reduced.ring.onsite_hund.spinfree.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.110e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 9.714e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-14).
   * ``chiq`` (comparator ``reduced_blocks``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.110e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 1.110e-16. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``reduced.ring.onsite_ising.spinfree.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.110e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 9.714e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-14).
   * ``chiq`` (comparator ``reduced_blocks``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.249e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 1.249e-16. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

``reduced.ring.onsite_pairlift.spinfree.mu``
   * ``chi0q`` (comparator ``identity``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.110e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 9.714e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-14).
   * ``chiq`` (comparator ``reduced_blocks``): agree to ``1.0e-14``.
     Provenance: measured on two development machines at commit b922a1c13b85: macOS arm64 (macOS-26.5.2-arm64-arm-64bit-Mach-O, Python 3.13.13, numpy 1.26.4, scipy 1.17.1) max|diff| 1.110e-16; Linux x86_64 (Linux-5.4.0-216-generic-x86_64-with-glibc2.31, Python 3.11.15, numpy 1.26.4, scipy 1.13.1) max|diff| 9.714e-17. Resulting atol 1.0e-14 (10x the larger residual, floored at 1e-15, rounded up to a power of ten; policy ceiling 1.0e-12).

Related tests
*************

A maintenance index of tests elsewhere in the suite that cover the
same ground, each with the note recorded alongside it. They are
listed for orientation only: the claims on this page are
established by ``tests/test_rpa_flex_equivalence_table.py`` alone.

``general.ring.offsite_coulombinter_interorb.mu``
   * FLEX -- ``tests.test_flex_offsite_general::TestOffsiteGeneralFLEX::test_two_orbital_offsite_interorbital_classes_match_the_ring``: The same fixture element-complete equal between RPA and FLEX for off-site inter-orbital CoulombInter (#181 Tier 1).

``general.ring.offsite_coulombinter_sameorb.mu``
   * FLEX -- ``tests.test_flex_offsite_general::TestOffsiteGeneralFLEX::test_one_orbital_offsite_v_matches_the_rpa_ring_exactly``: The same fixture (tests/rpa/input's coulombinter.dat, filling=0.75) element-complete equal between RPA and FLEX for off-site same-orbital CoulombInter.

``general.ring.offsite_coulombinter_sameorb.subshape``
   * FLEX -- ``tests.test_flex_offsite_general::TestOffsiteGeneralFLEX::test_offsite_under_sublattice_folding_matches_the_ring``: Same-orbital off-site CoulombInter under SubShape=(2,1,1) and (4,1,1) element-complete equal between RPA and FLEX (#181 Tier 1).

``general.ring.offsite_exchange.flexreject``
   * FLEX -- ``tests.test_flex_offsite_general::TestOffsiteGeneralFLEX::test_rejected_offsite_classes``: FLEX rejects off-site Exchange under calc_scheme='general'.

``general.ring.offsite_hund.mu``
   * FLEX -- ``tests.test_flex_offsite_general::TestOffsiteGeneralFLEX::test_two_orbital_offsite_interorbital_classes_match_the_ring``: The same fixture element-complete equal between RPA and FLEX for off-site Hund (#181 Tier 1).

``general.ring.offsite_ising.mu``
   * FLEX -- ``tests.test_flex_offsite_general::TestOffsiteGeneralFLEX::test_two_orbital_offsite_interorbital_classes_match_the_ring``: The same fixture element-complete equal between RPA and FLEX for off-site Ising (#181 Tier 1).

``general.ring.offsite_pairhop.flexreject``
   * FLEX -- ``tests.test_flex_offsite_general::TestOffsiteGeneralFLEX::test_rejected_offsite_classes``: FLEX rejects off-site PairHop under calc_scheme='general' (this fixture is one of the enumerated subTest rows since #181 Tier 1).

``reduced.ring.offsite_coulombinter.mu``
   * RPA -- ``tests.test_rpa_flex_oneshot_equivalence::TestReducedOneShot::test_matrix_cells``: An existing one-shot comparison over this same fixture and filling pins the chiq uu-ud/uu+ud block equivalence under calc_scheme='reduced'.

``reduced.ring.onsite_coulombinter.spindiag.mu``
   * RPA -- ``tests.test_rpa_flex_equivalence_table::TestReducedSpinDiagAcceptance::test_reduced_rpa_also_solves_with_extern``: RPA accepts the Extern/Zeeman spin-diagonal route (the tests/equivalence_input/spin fixture) under calc_scheme='reduced'.
   * FLEX -- ``tests.test_rpa_flex_equivalence_table::TestReducedSpinDiagAcceptance::test_reduced_flex_solves_with_extern_coulombinter``: FLEX accepts the Extern/Zeeman spin-diagonal route (the tests/equivalence_input/spin fixture) under calc_scheme='reduced' with CoulombInter.

``reduced.ring.onsite_coulombintra.spindiag.mu``
   * RPA -- ``tests.test_rpa_flex_equivalence_table::TestReducedSpinDiagAcceptance::test_reduced_rpa_also_solves_with_extern``: RPA accepts the Extern/Zeeman spin-diagonal route (the tests/equivalence_input/spin fixture) under calc_scheme='reduced'.
   * FLEX -- ``tests.test_rpa_flex_equivalence_table::TestReducedSpinDiagAcceptance::test_reduced_flex_solves_with_extern_coulombintra``: FLEX accepts the Extern/Zeeman spin-diagonal route (the tests/equivalence_input/spin fixture) under calc_scheme='reduced' -- the recorded expectation for this row is that both solvers accept it.

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
