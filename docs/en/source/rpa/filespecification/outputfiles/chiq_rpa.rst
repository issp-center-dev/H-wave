.. highlight:: none

.. _subsec:chiq_rpa:

chiq and chi0q
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The susceptibility matrix and the irreducible susceptibility matrix are exported
in NumPy zip (npz) format.
Using the string (referred to as *chiq_str*) specified by the keyword ``chiq`` or ``chi0q``
in ``file.output`` section in the parameter file,
the filename is chosen as *chiq_str*\ ``.npz``.

The file contains several arrays bound to the following keys:

- ``chiq`` or ``chi0q``:

  The susceptibility matrix or the irreducible susceptibility matrix. Their data layout is described in the following sections.

- ``momentum_convention``:

  The Fourier-sign provenance marker, ``"e_plus_ikR"`` (as of version 2.0): the
  momentum labels follow :math:`M(k) = \sum_R M(R) e^{+ikR}`. Loaders
  reject files recording a different value; files written before this
  field existed are accepted only when their content is elementwise even
  under :math:`q \to -q` (for which the two conventions coincide).

- ``freq_index``:

  The value or the range of Matsubara frequency is specified by ``matsubara_frequency`` parameter. The array bound to ``freq_index`` relates the index of the output data and the label of the actual Matsubara frequency.

- ``wavevector_unit`` and ``wavevector_index``:

  These arrays refer to the information of the wave number vectors. See :ref:`Output files of UHFk <Subsec:eigen_uhfk.dat>` for details.

- ``index_convention``:

  A string marker (always ``"spin_block"``) recording that the spin-orbital axes are ordered as spin*norb+orb (spin-block), distinct from UHFk's interleaved 2*orb+spin output. A ``chi0q`` file lacking this marker (produced before the convention fix) is rejected when reloaded in spin-orbital mode.

- ``coeff_tail``:

  The value of the ``coeff_tail`` parameter (Matsubara high-frequency tail correction) used when the susceptibility was computed. The tail correction changes :math:`\chi_0(\mathbf{q})` at :math:`O(1)` for moderate ``Nmat``, so runs that consume this file (e.g. ``hwave_sc``) warn when their own ``coeff_tail`` setting differs from the recorded value. When a pre-computed ``chi0q_init`` file is passed through, the input file's value is re-saved. The key is omitted when the producing run did not record one (older versions) and for FLEX output on the IR Matsubara basis (``matsubara_basis = "ir"``, including densified output), where the uniform-grid tail correction does not apply.

When the sublattice is considered, the indices of the wave numbers and the orbitals are
regarded as those of the sublattice.

The output file of ``chi0q`` can be used as a pre-calculated input of the irreducible
susceptibility by specifying the file to ``chi0q_init`` in ``file.input`` section.


Scheme provenance
^^^^^^^^^^^^^^^^^

Every ``chiq``/``chi0q`` file written by RPA and FLEX (2.0 and later) carries three plain-string fields, loaded as 0-d arrays (use ``.item()``):

- ``calc_scheme``: the scheme the file was computed with (``reduced`` | ``general``).
- ``calc_scheme_requested``: the value in the input (``auto`` | ``reduced`` | ``general``).
- ``scheme_resolution``: how it was decided. Closed vocabulary (changing it is an output-format change): ``explicit``, ``auto:ring_ladder``, ``auto:general_only``, ``auto:no_discarded_content``, ``auto:exact:diagonal_transfer``, ``auto:exact:folded_diagonal``, ``auto:mixed:transfer``, ``auto:mixed:extern``, ``auto:mixed:trans_mod``, ``auto:mixed:green_init``, ``auto:flex_forcing``.

Files written by 1.0.x lack these fields; readers do not require them.


Data format of chi0q
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Data format of ``chi0q`` relies on the presence of spin-orbital interaction and external field, and the value of ``mode.calc_scheme`` parameter, and takes one of the following:

- "spin-free" case:

  If ``enable_spin_orbital`` parameter is set to ``false``, or even if it is set to ``true`` when :math:`T_{\tilde\alpha\tilde\beta}(k)` is diagonal and symmetric with respect to spin degree of freedom, while the external field is not present, the spin-independent irreducible susceptibility matrix is exported. 
  
  - When ``calc_scheme = general``, the array format takes the form of ``ndarray(l,q,a,ap,b,bp)`` whose indices are given as follows:

    - ``l``: label of Matsubara frequency. The map from the label to the index is provided by the aforementioned array ``freq_index``.

    - ``q``: linearlized index of wave-number indices :math:`[ q_x\ q_y\ q_z ]`, where :math:`q = q_z + N_z\cdot(q_y + N_y\cdot q_x)`.

    - ``a``, ``ap``, ``b``, ``bp``: indices of the orbitals not including spin degree of freedom. They correspond to :math:`\alpha`, :math:`\alpha^\prime`, :math:`\beta`, :math:`\beta^\prime`. 

  - When ``calc_scheme = reduced``, the array format takes the form of ``ndarray(l,q,a,b)`` whose indices are same as the above.

- "spin-diagonal" case:

  If ``enable_spin_orbital`` parameter is set to ``false`` and the external field is present, or it is set to ``true`` while :math:`T_{\tilde\alpha\tilde\beta}(k)` is diagonal with respect to spin degree of freedom, the spin-up and spin-down components of the irreducible susceptibility matrix are exported.

  - When ``calc_scheme = general``, the array format takes the form of ``ndarray(s,l,q,a,ap,b,bp)``, where ``s = 0`` denotes spin-up component and ``s = 1`` does spin-down component. The other indices are same as the above.

  - When ``calc_scheme = reduced``, the array format takes the form of ``ndarray(s,l,q,a,b)``. The indices are same as above.

- "spinful" case:

  If ``enable_spin_orbital`` parameter is set to ``true``, and :math:`T_{\tilde\alpha\tilde\beta}(k)` takes a general form, the irreducible susceptibility matrix with the generalized orbital indices is exported.

  - When ``calc_scheme = general``, the array format takes the form of ``ndarray(l,q,a,ap,b,bp)``, where ``a``, ``ap``, ``b``, and ``bp`` corresponding to the generalized orbital indices including spin degree of freedom denoted by :math:`\tilde\alpha`, :math:`\tilde\alpha^\prime`, :math:`\tilde\beta`, and :math:`\tilde\beta^\prime`, respectively.

  - When ``calc_scheme = reduced``, the array format takes the form of ``ndarray(l,q,a,b)``, where ``a`` and ``b`` corresponding to the generalized orbital indices :math:`\tilde\alpha` and :math:`\tilde\beta`, respectively.


Data format of chiq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Data format of ``chiq`` takes the following form depending on the value of ``calc_scheme`` parameter:

- When ``calc_scheme = general``, the array format takes the form of ``ndarray(l,q,a,ap,b,bp)``, where ``a``, ``ap``, ``b``, and ``bp`` correspond to the generalized orbital indices including spin degree of freedom denoted by :math:`\tilde\alpha`, :math:`\tilde\alpha^\prime`, :math:`\tilde\beta`, and :math:`\tilde\beta^\prime`, respectively.

- When ``calc_scheme = reduced``, the array format takes the form of ``ndarray(l,q,a,b)``, where ``a`` and ``b`` correspond to the generalized orbital indices :math:`\tilde\alpha` and :math:`\tilde\beta`, respectively.

``chiq`` holds the longitudinal (ring) susceptibility. Under ``calc_scheme = general`` it has slots in which a pair of indices is spin-off-diagonal; ``reduced`` stores density-pair components only and has no such slots. Those slots are **not** the transverse susceptibility, and whenever the bubble is obtained by inflating a spin-free or spin-diagonal one they are **not computed** -- the inflation builds only the same-spin slots, leaving the rest identically zero. That covers every calculation without ``enable_spin_orbital``, and also an ``enable_spin_orbital`` calculation whose one-body Hamiltonian is detected as spin-free or spin-diagonal. Do not read those zeros as a computed transverse response. (In a genuinely spinful ``general`` calculation the bubble is built directly on the generalized orbital index, so these slots are computed and are generally nonzero; the transverse susceptibility is still ``chiq_pm``.) See :ref:`Algorithm<rpa_which_array>` for how to obtain the transverse channel.

.. admonition:: Migrating from ``calc_scheme = "squashed"`` (removed in 2.0)
   :class: note

   ``squashed`` computed the same susceptibility as ``reduced`` at several
   times the cost; the slots of its 8-axis output in which a pair of spin
   indices differed were structurally zero. Configurations now fail with an
   error at start-up; use ``calc_scheme = "reduced"``.

   Analysis scripts reading a legacy 8-axis file can convert it::

       chiq8 = data["chiq"]           # (l, q, s1, s2, a, s3, s4, b)
       l, q, _, _, norb, _, _, _ = chiq8.shape
       nd = 2 * norb
       chiq4 = np.zeros(chiq8.shape[:2] + (nd, nd), dtype=chiq8.dtype)
       for s1 in (0, 1):
           for s3 in (0, 1):
               chiq4[:, :, s1*norb:(s1+1)*norb, s3*norb:(s3+1)*norb] = \
                   chiq8[:, :, s1, s1, :, s3, s3, :]

   Only the pair-diagonal slots (``s1 == s2``, ``s3 == s4``) carry data;
   every other slot of the 8-axis array is exactly zero. Susceptibility
   files produced by older ``squashed`` runs remain loadable as
   ``chi0q_init`` under ``reduced``, subject to the same provenance
   validation as any legacy file (e.g. the momentum-convention marker):
   the two schemes always shared one bubble representation.

When ``calc_type = ring+ladder``, the ``chiq`` file additionally contains the array ``chiq_pm``, which holds the transverse susceptibility :math:`\chi_{+-}(q)`. Its layout is ``ndarray(l,q,a,ap,b,bp)`` where ``a``, ``ap``, ``b``, ``bp`` are **orbital** indices :math:`\alpha`, :math:`\gamma`, :math:`\beta`, :math:`\delta` that do *not* include the spin degree of freedom -- the spin structure is fixed by the :math:`+-` labels. It is therefore smaller than ``chiq``, whose corresponding axes run over the generalized (spin-orbital) indices. The longitudinal ``chiq`` is unaffected by the presence of the ladder: on identical input it is bit-identical to the ``calc_type = ring`` result.

When ``longitudinal_bond_channels = true`` (experimental, ``calc_type = ring``, see :ref:`Ch:Algorithm`), the ``chiq`` file additionally contains the bond-resolved longitudinal objects, all under the ``longitudinal_bond_`` prefix; ``chiq`` itself is unchanged. With :math:`B` bond channels (channel 0 is the on-site :math:`R = 0`), :math:`n_d = n_{\rm orb}^2` and :math:`N_D = B\, n_d`: ``longitudinal_bond_chi_s`` and ``longitudinal_bond_chi_c`` are the dressed static spin and charge susceptibilities, ``ndarray(q, I, J)`` with the bond-major index :math:`I = m\, n_d + l_1 n_{\rm orb} + l_2` (``longitudinal_bond_index_order``); ``longitudinal_bond_chiq_s_static`` and ``longitudinal_bond_chiq_c_static`` are their :math:`(m = 0, m' = 0)` blocks, ``ndarray(q, a, b, c, d)`` over orbital indices, i.e. the static spin/charge channels of the standard ring with the off-site exchange crossing included. ``longitudinal_bond_delta_r`` (``(B, 3)`` integer displacements), ``longitudinal_bond_reverse`` (the channel of :math:`-R`), ``longitudinal_bond_spatial_shape``, ``longitudinal_bond_q_convention``, ``longitudinal_bond_spin_mode``, ``longitudinal_bond_normalization``, ``longitudinal_bond_types`` (always the capability tuple ``CoulombInter, Hund, Ising`` -- the types whose bond blocks the channel builds, declared or not; an undeclared type contributes zero blocks), ``longitudinal_bond_max_shells`` (``-1`` when unset), ``longitudinal_bond_cond_min_s`` / ``longitudinal_bond_cond_min_c`` (the smallest conditioning score of the spin / charge RPA denominator over :math:`q`; the run is refused when it reaches the instability floor) and ``longitudinal_bond_schema`` (``1``) describe the layout.


Example for reading data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following code is an example for reading the data from the output file.

.. code-block:: python

    import numpy as np
    data = np.load("chiq_str.npz")

    chiq = data["chiq"]
    freq_index = data["freq_index"]


.. raw:: latex
