.. _Ch:FAQ:

FAQ
****

Separating the spin and charge contributions to superconductivity
=================================================================

**Q. I want to analyze the superconductivity by separating the charge and spin
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
comparing the resulting :math:`\lambda` shows which fluctuation channel
dominates the pairing. Both flags default to ``false`` (ordinary runs are
unaffected), and the diagnostic works for both ``pairing_type = "singlet"`` and
``"triplet"``.

.. note::

   The instantaneous (bare) interaction term is retained in every case, and the
   linearized-gap eigenvalue problem is nonlinear, so the eigenvalues do **not**
   add up: :math:`\lambda_{\mathrm s} + \lambda_{\mathrm c} \neq
   \lambda_{\mathrm{full}}` in general, where :math:`\lambda_{\mathrm s}` and
   :math:`\lambda_{\mathrm c}` denote the spin-plus-bare and charge-plus-bare
   runs, respectively. Use the decomposition to compare the relative strength
   of the two channels, not as an exact additive split.

See :ref:`sc_channel_decomposition` for the vertex formulas and further details.
