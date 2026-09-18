.. highlight:: none

.. _Ch:Config_rpa:

Parameter files
--------------------------------

The parameter file specifies calculation conditions and parameters for H-wave
in TOML format. It is composed of the following three sections.

#. ``mode`` section for specifying calculation conditions,

#. ``log`` section for setting standard outputs,

#. ``file`` section for setting file paths: It contains ``input`` and ``output`` subsections.

An example of the file is shown below:

.. literalinclude:: ../../sample/input.toml


File format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
TOML format


Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``mode`` section
================================

- ``mode``

  **Type :**
  String

  **Description :**
  This parameter specifies the calculation mode.
  Set to ``"RPA"`` for calculations of the Random Phase Approximation,
  or to ``"FLEX"`` for the Fluctuation Exchange Approximation.

- ``enable_spin_orbital`` (default value is ``false``)

  **Type :**
  Boolean

  **Description :**
  This parameter specifies whether to allow spin-orbital interaction.
  If it is set to ``true``, the orbital indices in the Transfer term use the interleaved
  spin-orbital convention: the combined index is :math:`2 \alpha + s`, where
  :math:`\alpha` is the physical-orbital index (0-based) and :math:`s \in \{0, 1\}`
  is the spin index.

  In spin-orbital mode, the ``Norbit`` value in the geometry file (``geom.dat``) is
  the **spin-orbital count** (= 2 × the number of physical orbitals = Wannier90
  ``num_wann``), the same convention as UHFk.

  .. note::

     **Migration (RPA):** the geometry ``norb`` for spin-orbital input is now the
     spin-orbital count; double any pre-existing RPA spin-orbital ``geom.dat``
     ``Norbit``.

.. _rpa_calc_scheme_auto:

- ``calc_scheme`` (default value is ``"auto"``)

  **Type :**
  String

  **Description :**
  This parameter specifies how the spin and orbitals are treated in the calculation. The parameter takes one of the following options.

  - ``general``: Generalized orbitals combining spins and orbitals are considered. The susceptibility matrix takes the most general form, with the size of :math:`N_\text{orb}^4 N_\text{spin}^4 N_k N_\omega`.

  - ``reduced``: Generalized orbitals combining spins and orbitals are considered. The components of the susceptibility matrix with :math:`\alpha=\alpha^\prime` and :math:`\beta=\beta^\prime` are considered. The size of the matrix turns to :math:`N_\text{orb}^2 N_\text{spin}^2 N_k N_\omega`. For the two-body interaction terms, only CoulombIntra, CoulombInter, Ising and Hund are allowed. 

  - ``auto``: the scheme is chosen so that the result is **exact for the declared input**, preferring ``reduced`` where exactness is provable ("conservative exact"). ``reduced`` is kept only when every declared interaction is density-only, or when the one-body Hamiltonian (``Transfer``, an active ``Extern``, ``trans_mod`` / ``green_init``) conserves the orbital flavour so the discarded cross-family vertex sectors are unreachable; otherwise ``general`` is selected. ``ring+ladder``, ``Exchange`` and ``PairHop`` always select ``general``. This option is not available when only ``chi0q`` is to be calculated. The decision is recorded in every output file as ``scheme_resolution`` (see the ``chiq`` output page for the token list).

  .. note::

     **Version 2.0 / 1.0.x:** in 1.0.x ``auto`` selected ``reduced`` for CoulombInter/Hund/Ising regardless of orbital hybridisation, which is an approximation for hybridised multi-orbital models (measured 2.3e-4 / 3.3e-4 / 3.2e-4 relative on the reference 2-orbital fixture; unbounded near an RPA instability). Since 2.0 such inputs are promoted to ``general`` and a warning names the reason and the cost (``chiq`` becomes rank-6; memory/solve cost grow from :math:`O(N_d^2)/O(N_d^3)` to :math:`O(N_d^4)/O(N_d^6)` per :math:`(q,\omega)`). **Migration:** (1) to keep the 1.0.x behaviour for an existing post-processor, request ``calc_scheme = "reduced"`` explicitly -- the solver then warns once that the result is an approximation for that input; (2) to adapt a 4-axis density-channel reader to the rank-6 layout, take the density-pair slots ``chiq[:, :, a, a, b, b]`` (``numpy.einsum("kqaabb->kqab", chiq)``), which reproduces the 4-axis array exactly when ``reduced`` was exact. ``auto`` is resolved when the calculation starts (``read_init``/``solve``), not at construction; ``RPA.preview_scheme()`` returns the decision without running.

- ``calc_type`` (default value is ``"ring"``)

  **Type :**
  String

  **Description :**
  This parameter specifies which RPA diagrams to include.

  - ``ring``: Standard RPA (ring diagram only). Computes the longitudinal susceptibility.

  - ``ring+ladder``: Includes the transverse (ladder) susceptibility :math:`\chi_{+-}(\mathbf{q})` in addition to the standard RPA. This requires the ``general`` calculation scheme (automatically selected). See :ref:`Ch:Algorithm` for details.


``mode.param`` section
================================

``mode.param`` section contains the parameters for the calculation.

- ``CellShape``

  **Type :**
  Integer array

  **Description :**
  This parameter specifies the shape of the lattice Lx, Ly, Lz.

- ``SubShape``

  **Type :**
  Integer array (default value is ``[`` Lx, Ly, Lz ``]``)

  **Description :**
  This parameter specifies the shape of the sublattice Bx, By, Bz.

- ``T``

  **Type :**
  Float (default value is 0)

  **Description :**
  This parameter specifies the temperature.
  It must be greater than or equal to zero.

- ``mu``

  **Type :**
  Float or None (default value is None)

  **Description :**
  This parameter specifies the chemical potential :math:`\mu`.
  If it is not specified, the value of :math:`\mu` will be calculated so that
  the expectation value of the number of electrons equals to ``Ncond``.
  If both ``mu`` and ``Ncond`` or ``filling`` are specified, the program terminates with error.

- ``Ncond``

  **Type :**
  Integer

  **Description :**
  This parameter specifies the number of conduction electrons.
  It must be greater than or equal to one.

- ``filling``

  **Type :**
  Float

  **Description :**
  This parameter specifies the filling ratio of electrons with respect to the number of states.
  Both ``Ncond`` and ``filling`` are specified, the program will be terminated with error.

- ``Ncond_round_mode``

  **Type :**
  String (default value is ``"strict"``)

  **Description :**
  This parameter specifies how the number of electrons calculated from the ``filling`` parameter is rounded to an integer value when the temperature is zero. The parameter must take one of the following values.

    - ``as-is``:  the value is not rounded to an integer. (returns a floating-point number)
    - ``round-up``:  the value is rounded up.
    - ``round-down``:  the value is rounded down.
    - ``round-off``:  the value is rounded to the closest integer. (0.5 is rounded up.)
    - ``round``:  the value is rounded by ``round`` function. (0.5 is rounded down.)
    - ``strict``:  if the value is not an integer value, the program terminates with error.
    - ``exact``:  if the value is not an integer value, a warning message will be shown and the value is rounded to an integer as ``round``.

- ``Nmat``

  **Type :**
  Integer (default value is 1024)

  **Description :**
  This parameter specifies the cut-off of Matsubara frequency.
  It must be an even number greater than zero. Matsubara frequency is defined as follows:

      - Boson: :math:`\omega_n = \dfrac{2\pi (n-\texttt{Nmat}/2)}{\beta}`
      - Fermion: :math:`\omega_n = \dfrac{\pi (2n+1-\texttt{Nmat})}{\beta}`

  with the indices :math:`n` between 0 and ``Nmat-1``.

- ``coeff_tail``

  **Type :**
  Float (default value is 0.0)

  **Description :**
  This parameter specifies the magnitude of the correction when correcting the tails of the Fourier transformation.
  After Fourier transforming the diagonalized one-body Green function to the imaginary time representation by subtracting :math:`\texttt{coeff\_tail}/(i \omega_n)`, the term :math:`-\beta/2\cdot\texttt{coeff\_tail}` is added to the one-body Green function.
  In the FLEX solver the same tail treatment is applied to the *dressed* Green function before the bare susceptibility :math:`\chi_0(q)` is computed, so that ``coeff_tail`` accelerates the frequency summation without changing the physical result. (The FLEX self-energy convolution keeps the full Green function and is unaffected.)
  Since version 2.0 the susceptibility kernels also restore the Green
  function's equal-time discontinuity at the bubble's :math:`\tau = 0`
  sample (the tail piece carries the jump; the sample is the mean of the
  two branches), which makes ``coeff_tail = 1.0`` converge at
  :math:`O(1/N_{\rm mat}^2)`. Earlier versions omitted this endpoint and
  ``coeff_tail`` then *slowed* the convergence by a constant factor
  (still :math:`O(1/N_{\rm mat})`); results produced with a nonzero
  ``coeff_tail`` before this fix are not comparable with current ones
  at the same ``Nmat``.
  The value must be a finite real number; ``NaN`` and infinities are
  rejected. Only ``0.0`` (off) and ``1.0`` (the physical :math:`1/i\omega_n`
  coefficient) are recommended: fractional values cancel only part of the
  equal-time jump, remain :math:`O(1/N_{\rm mat})` and can converge more
  slowly than ``coeff_tail = 0.0``.
  ``chi0q.npz`` files written with a nonzero ``coeff_tail`` carry a
  ``tail_endpoint = "branch_mean_v1"`` marker recording the endpoint
  treatment; ``chi0q_init`` and the ``hwave_sc`` chi0q loader refuse a
  nonzero-tail file without it (produced before the fix), since the
  pre-fix error cannot be detected from the array itself. Recompute such
  bubbles instead of reusing the files.

- ``spinful_vertex_exchange``

  **Type :**
  Boolean (default value is ``true``)

  **Description :**
  Spinful (``enable_spin_orbital``) calculations resum the susceptibility
  with a single vertex tensor. Since version 2.0 that tensor is the
  antisymmetrized bare particle-hole vertex: the direct (ring) wiring
  plus the exchange wiring of the on-site interaction terms. The
  exchange part is what corrects the spin-flip pair components of
  :math:`\chi(q)` (in non-spin-orbital calculations the analogous
  content is provided separately by ``calc_type = "ring+ladder"``);
  without it those components are returned as the bare bubble at any
  interaction strength, and, because spin is not conserved, the error
  leaks into every component. The construction was verified against
  exact diagonalization at first order in the coupling for
  CoulombIntra, Exchange and PairLift, and reproduces the established
  transverse (ring+ladder) series in the spin-conserving limit.
  Setting ``spinful_vertex_exchange = false`` restores the previous
  ring-only vertex (results produced before this fix): use it only to
  reproduce old runs. The exchange wiring of OFF-site interaction
  terms depends on both fermionic momenta and is not representable in
  this resummation; it remains excluded (as it is in the
  non-spin-orbital ladder).

- ``longitudinal_bond_channels``

  **Type :**
  Boolean (default value is ``false``)

  **Description :**
  **Experimental.** Enables the bond-resolved longitudinal (spin/charge)
  channel for ``calc_type = "ring"`` with ``calc_scheme = "general"`` (or
  ``"auto"`` resolving to it) on a spin-free system (see
  :ref:`rpa_longitudinal_bond`). The exchange (Fock)
  crossing of the off-site ``CoulombInter``, ``Hund`` and ``Ising``
  terms -- which the standard ring omits -- is carried on a bond-enlarged
  pair basis, and the dressed static susceptibilities are written to the
  ``chiq`` file under the ``longitudinal_bond_*`` keys; the ``chiq``
  array inside that file is not overwritten and stays the standard ring
  result. The gate requires at least one
  declared off-site ``CoulombInter``/``Hund``/``Ising`` shell (a
  declared-zero coefficient counts), real coefficients, an even
  ``Nmat``, no sublattice, no ``chi0q_init`` and no
  ``enable_spin_orbital``; an off-site ``Exchange`` or ``PairHop``
  declaration is rejected with an error. It cannot be combined with
  ``transverse_bond_channels = true``. In ``mode = "FLEX"`` the same key
  enables the self-consistent bond-resolved channel of the Hartree-Fock
  FLEX (experimental; requires ``flex_hartree_fock = true``,
  ``calc_scheme = "general"``, a spin-free system, the uniform Matsubara
  grid, no sublattice and an even ``Nmat``; see
  :ref:`flex_bond_hf`). With ``longitudinal_bond_channels = true``
  (which requires ``flex_hartree_fock = true``) and ``gpu = true``,
  the dressing, the effective interaction and the self-energy
  transport run on the GPU (CuPy); the bond arrays stay in host memory
  and one frequency batch is transferred at a time; the batch is
  chosen against both the host cap and the free device memory
  (``longitudinal_bond_freq_batch`` overrides both, see there). Results
  agree with the CPU path to round-off; the outputs record
  ``longitudinal_bond_device`` and ``longitudinal_bond_nb``. The
  bond-specific companions below keep their meaning;
  ``longitudinal_bond_output_full`` and ``longitudinal_bond_freq_batch``
  are FLEX-only.

- ``longitudinal_bond_max_shells``

  **Type :**
  Integer (default: keep every declared off-site shell)

  **Description :**
  Keeps only the ``longitudinal_bond_max_shells`` shortest off-site
  shells (must be >= 1). This is NOT an approximation knob: a truncation
  that would drop a shell carrying a declared nonzero coefficient is
  rejected, so it can only remove declared-zero (spectator) shells; to
  exclude a physical interaction, set its coefficient to zero in the
  interaction file. Ignored with a warning unless
  ``longitudinal_bond_channels = true``.

- ``longitudinal_bond_memory_cap_gb``

  **Type :**
  Float (default value is ``8.0``)

  **Description :**
  Cap on the ESTIMATED peak host (CPU) memory of the bond-resolved solve,
  in binary GiB (GPU memory is not covered); the estimate is logged and
  the run is refused before any expensive step when it exceeds the cap.
  The refusal message reports the estimate and the shapes behind it; on
  a machine with enough physical memory, raise the cap to proceed.
  Ignored with a warning unless ``longitudinal_bond_channels = true``.

- ``longitudinal_bond_guard_freqs``

  **Type :**
  String (default value is ``"all"``)

  **Description :**
  Which bosonic Matsubara frequencies the conditioning guard of the
  bond-resolved dressing inspects. ``"all"`` (default) checks every
  ``(nu, q)`` denominator with a singular-value decomposition and
  refuses the run when one is singular or nearly singular (the
  behaviour of previous versions). ``"static"`` checks only the zero
  frequency and validates the other slices by the solve residual
  ``||(1 -/+ chibar V) chi - chibar|| / max(1, ||chibar||) <= 1e-6``; it
  is a reduced diagnostic that removes the dominant CPU cost of the
  guard and is recorded in the outputs (``longitudinal_bond_guard_freqs``
  member). The default does not depend on ``gpu``. Ignored with a
  warning unless ``longitudinal_bond_channels = true``.

- ``longitudinal_bond_pairing``

  **Type :**
  String (``"none"``, ``"singlet"``, ``"triplet"``, ``"both"``;
  case-insensitive; default ``"none"``; FLEX mode only)

  **Description :**
  With ``longitudinal_bond_channels = true`` (required; the key is refused
  otherwise, as it is with ``IterationMax = 0``), solves the linearized
  Eliashberg equation at the end of the FLEX solve with the bond-resolved,
  frequency-resolved pairing vertex built from the same ``S``, ``C`` and
  dressed ``chi_s``, ``chi_c`` the self-energy used (the in-process
  counterpart of ``[eliashberg] bond_channels`` below). The solver controls
  are read from the ``[eliashberg]`` table of the SAME input file (optional;
  every key defaults to the value ``hwave_sc`` uses; ``pairing_type`` is
  refused there -- the channel(s) are selected by this key instead; ``gpu``
  follows the FLEX solve). Results: ``eliashberg_bond_<type>.npz``,
  ``gap_bond_<type>.dat``, ``eigenvalue_bond_<type>.dat`` per requested
  channel (``[file.output]`` keys of the same names, see below). A failure
  of the pairing step is recorded in the log; that channel's files are then
  simply not written (the FLEX outputs are written first and are never lost
  to it). Memory: see :ref:`the FLEX tutorial section on pairing with the
  bond-resolved vertex <flex_bond_pairing_tutorial>`; for large models set
  ``[eliashberg] matsubara_basis = "ir"``.

- ``flex_hartree_fock``

  **Type :**
  Boolean (default value is ``false``; FLEX mode only)

  **Description :**
  **Experimental.** Adds the self-consistent Hartree-Fock self-energy
  :math:`\Sigma_{\rm HF}[G]` of EVERY accepted interaction term (on-site
  and off-site) to the FLEX self-energy, :math:`\Sigma = \Sigma_{\rm HF}
  + \Sigma_{\rm fluct}`, recomputed from the dressed Green function every
  iteration (see :ref:`flex_bond_hf`). Requires ``calc_scheme =
  "general"`` on a spin-free system, the uniform Matsubara grid
  (``matsubara_basis = "ir"`` is refused), no sublattice, no external
  field and an even ``Nmat``. ``gpu = true`` is accepted: on its own
  the ordinary FLEX GPU path applies (see ``gpu`` below); with
  ``longitudinal_bond_channels = true`` the bond-resolved GPU path
  described there applies instead.
  A user mean field (``trans_mod`` / ``green_init``) is taken as the
  initial static self-energy instead of being folded into the band;
  ``sigma.npz`` gains the members ``sigma_static``, ``sigma_fluct``,
  ``sigma_convention = "split"`` and a provenance block, and only a
  ``"split"`` archive is accepted as ``sigma_init`` (see ``sigma_init``
  and ``hwave_sigma_split``). Mandatory (``true``) when
  ``longitudinal_bond_channels = true`` in FLEX mode. With ``false``
  every output is unchanged.

- ``longitudinal_bond_output_full``

  **Type :**
  Boolean (default value is ``false``; FLEX mode only)

  **Description :**
  Writes the full dynamic bond-resolved spin and charge susceptibilities
  ``chi_s_w`` / ``chi_c_w`` (``ndarray(l, q, I, J)`` over the bosonic
  frequencies) of the last SCF map into the dedicated archive named by
  ``[file.output] longitudinal_bond`` (default ``longitudinal_bond.npz``).
  Doubles the persistent memory of the bond-resolved solve; the archive
  size is logged before writing. Ignored with a warning unless
  ``longitudinal_bond_channels = true``.

- ``longitudinal_bond_freq_batch``

  **Type :**
  Integer (default: chosen automatically; FLEX mode only)

  **Description :**
  Number of bosonic frequencies dressed at once by the bond-resolved
  FLEX (in ``[1, Nmat]``). By default the largest batch that keeps the
  estimated peak memory under ``longitudinal_bond_memory_cap_gb`` is
  used (and, with ``gpu = true``, under the free device memory as
  well); an explicit value is validated against the same cap(s) --
  the host cap always, and the free device memory measured when the
  bond solve starts when ``gpu = true`` -- and refused when it exceeds
  either. Ignored with a warning unless
  ``longitudinal_bond_channels = true``.

- ``flex_second_order``

  **Type :**
  String (default value is ``"local"``; FLEX mode with ``calc_scheme =
  "general"`` only)

  **Description :**
  Selects the second-order part of the FLEX effective interaction on the
  general (full-vertex) path; see
  :ref:`flex_second_order_kernel` for the formulas.

  - ``"local"`` (default) builds the exact LOCAL second order: the
    complete second order of the on-site interaction (``CoulombIntra``,
    ``CoulombInter``, ``Hund``, ``Ising``, ``Exchange``, ``PairHop``,
    ``PairLift`` at :math:`R = 0`), the direct (bubble) skeleton of the
    off-site density terms (``CoulombInter``, ``Hund``, ``Ising`` at
    :math:`R \neq 0`) and every mixed on-site/off-site diagram (the
    off-site vertex enters those in its crossed placement; see
    :ref:`flex_second_order_kernel`). The only
    second-order class it does not carry is the exchange skeleton of TWO
    off-site vertices, which is not representable by a
    :math:`q`-only vertex: full second-order accuracy for an off-site
    interaction additionally needs ``longitudinal_bond_channels = true``.
  - ``"takimoto"`` keeps the legacy Takimoto-Hotta-Ueda expression
    verbatim; its second-order (double-counting subtraction) term is
    :math:`-\tfrac{1}{4}(\hat{U}^s + \hat{U}^c)\,\bar\chi\,
    (\hat{U}^s + \hat{U}^c)`, with :math:`\bar\chi` the bare bubble
    (the full expression is given in
    :ref:`flex_second_order_kernel`). It is the reproduction path for
    results produced with H-wave 2.0.0 and remains available throughout
    the 2.x series.

  **Applicability:** The key is meaningful for ``mode = "FLEX"`` with
  ``calc_scheme = "general"`` only. The default is applied AFTER FLEX
  resolves ``calc_scheme = "auto"``, so an absent key never raises
  anywhere. An explicit key with an explicit ``calc_scheme = "reduced"``
  is refused at start-up; an explicit key with ``calc_scheme = "auto"``
  that resolves to ``"reduced"`` for the declared interaction set is
  refused with a message naming that resolution (an explicit key never
  promotes ``auto`` to ``general``). A malformed value is refused before
  applicability is examined. In ``mode = "RPA"`` an explicit key is
  ignored with the usual one-off warning about FLEX-only keys; the UHF
  solvers ignore it silently. The selected value is logged once at
  start-up and recorded in every output archive of a general-scheme FLEX
  run (see :ref:`the output reference <rpa_chiq_provenance>`).

  **Degenerate declarations:** Under ``"local"`` an ON-SITE same-orbital row
  (:math:`R = 0`, :math:`\alpha = \beta`) of ``CoulombInter``, ``Hund``,
  ``Ising``, ``Exchange``, ``PairHop`` or ``PairLift`` is refused at
  start-up -- a row with ``rx = ry = rz = 0`` and identical orbital
  indices :math:`\alpha = \beta`. Such a row is not a two-body term: it
  reduces to a one-body level shift, to an effective ``CoulombIntra``, or
  to identically zero. The error
  message names the row and gives the equivalent declaration for that
  type. Rewrite the interaction file as the message says, or set
  ``flex_second_order = "takimoto"`` as an immediate workaround (the
  legacy expression accepts every row the general path accepted before).

  **Cost and memory:** Measured with ``tests/sc/second_order_cost.py`` (run
  from the repository root; the 2-orbital ``tests/rpa/input_2orb`` geometry
  and transfer, and the synthetic 3-orbital square lattice the script writes,
  :math:`L = 8`, ``Nmat = 128``, minimum over repeated runs, 2026-09):
  assembling :math:`V_{\rm eff}` under ``"local"`` costs about 1.3x to 1.6x
  the legacy expression (2- and 3-orbital inputs; the ratio varies from run
  to run within that band) and a full FLEX iteration about 1.1x to 1.25x.
  Memory goes the other way: ``"local"`` accumulates the second order in
  frequency batches, while the legacy expression materialises three
  full-size terms, so under ``"takimoto"`` the peak of this step is about
  three times that of ``"local"`` (:math:`V_{\rm eff} + 3\,C` versus
  :math:`V_{\rm eff} + C/4 + F`, with :math:`C` the size of
  :math:`V_{\rm eff}` and :math:`F` the small compiled interaction
  factors); its transient buffers alone are twelve times larger, so a
  large reproduction run may need more RAM than the same run under the
  default.

  **Convergence:** The SCF trajectory can change with the kernel. On the
  2-orbital self-consistency fixture (on-site :math:`U`, :math:`U'`,
  ``Hund`` plus off-site :math:`V`, Anderson mixing) both values needed
  the same number of iterations (9 without and 11 with
  ``flex_hartree_fock = true``), but that is one observation, not a
  guarantee; see the :ref:`FLEX tutorial <flex_second_order_tutorial>`
  for what to do when a run that converged before does not.

- ``matsubara_frequency``

  **Type :**
  Integer, List of Integers, or String (default value is ``"all"``)

  **Description :**
  This parameter specifies the indices of Matsubara frequency for which the susceptibility matrix :math:`\chi(\vec{q})` is calculated.
  The value must be one of the following:

    - *an integer value* : a single index value.
    - ``[`` *min*, *max* (, *step*) ``]`` : every *step* index from *min* to *max*. If *step* is omitted, it is assumed to be 1.
    - all : all indices
    - center : corresponds to ``Nmat/2``.
    - none : nothing will be calculated.

  When the susceptibility matrix :math:`\chi(\vec{q})` or the irreducible susceptibility matrix :math:`\chi_0(\vec{q})` are stored to files, the values at the specified freqneucy are exported.


- ``coeff_extern``

  **Type :**
  Float (default value is 0.0)

  **Description :**
  This parameter specifies the coefficient :math:`h` of the external field given by the form :math:`\pm h H_{\alpha\beta}(r_{ij})`. The definition of the matrix :math:`H_{\alpha\beta}(r_{ij})` will be provided by an input file. The sign :math:`+` and :math:`-` correspond to spin up and down, respectively.
  


- ``RndSeed``

  **Type :**
  Integer (default value is 1234)

  **Description :**
  This parameter specifies the seed of the random number generator.

- ``ene_cutoff``

  **Type :**
  Float (default value is 100.0)

  **Description :**
  This parameter specifies the upper cutoff of the exponent in the Fermi distribution function to avoid overflow during the calculation.

- ``gpu``

  **Type :**
  Boolean (default value is false)

  **Description :**
  When set to true, the main computation (the Green's function, the chi0q
  FFT pair bubble, the spin inflation, and the batched RPA solve; for FLEX
  the whole SCF loop) runs on a GPU via CuPy. When CuPy or a CUDA device is
  unavailable, the solver warns and falls back to the CPU (numpy) path with
  an identical result. Install CuPy as the precompiled binary wheel matching
  your CUDA version (e.g. ``pip install cupy-cuda12x`` for CUDA 12.x); see
  the `CuPy installation guide <https://docs.cupy.dev/en/stable/install.html>`_.
  ``flex_hartree_fock = true`` on its own is accepted and uses this
  ordinary GPU path; with ``longitudinal_bond_channels = true`` (which
  requires ``flex_hartree_fock = true``) the bond-resolved GPU path
  described under ``longitudinal_bond_channels`` above applies instead.

- ``fft_workers``

  **Type :**
  Integer (default value is 1)

  **Description :**
  Number of worker threads for the spatial FFTs (parallelized via
  ``scipy.fft``). The default 1 keeps the serial numpy path, unchanged
  from previous releases (opt-in); -1 uses all cores. Ignored on the
  GPU. Set a smaller number when running several calculations
  concurrently.

- ``mixing_scheme``

  **Type :**
  String (default value is "linear"; FLEX mode only)

  **Description :**
  Self-energy update scheme of the FLEX SCF loop. ``"linear"`` is the
  conventional linear mixing; ``"anderson"`` enables Anderson acceleration
  (Pulay/DIIS-type extrapolation over a short iterate/residual history),
  which reaches the same fixed point in far fewer iterations. Anderson
  acceleration is less sensitive to step-size instability than linear
  mixing, so a somewhat larger ``Mix`` (e.g. 0.3--0.5) can reduce the
  iteration count further. Falls back to a plain linear step automatically
  if the history becomes degenerate.

- ``anderson_depth``

  **Type :**
  Integer (default value is 5; FLEX mode only)

  **Description :**
  History depth of the Anderson acceleration. Memory grows by 2*depth
  sigma-sized arrays (kept on the device under GPU execution).


.. _eliashberg_bond_dynamic_config:

``eliashberg`` section
================================

This section is read by ``hwave_sc`` (the linearized Eliashberg equation
solver) and, for the two keys marked "in-process" below, by ``hwave``
itself when ``[mode.param] longitudinal_bond_pairing`` (above) requests it.
Only the keys specific to the bond-resolved DYNAMIC pairing path are
documented here; the full ``[eliashberg]`` table (``frequency``,
``pairing_type``, ``solver_mode``, ``matsubara_basis``, ``ir_tol``,
``ir_wmax``, ...) is documented in
:ref:`the dynamic-frequency section of the Eliashberg solver tutorial
<sc_dynamic_frequency>`.

- ``bond_channels``

  **Type :**
  Boolean (default value is ``false``; post-processing, ``hwave_sc`` only)

  **Description :**
  With ``frequency = "dynamic"`` and ``chi0q_mode = "flex"``, replaces the
  on-site dynamic pairing kernel by the bond-resolved one: the pairing
  vertex is built from a FLEX run's dynamic bond-resolved susceptibilities
  (``S_bond``, ``C_bond`` and the bond-resolved ``chi_s_w`` / ``chi_c_w``)
  instead of ``chiq_s.npz`` / ``chiq_c.npz``. The producing FLEX run must
  have used ``longitudinal_bond_channels = true`` AND
  ``longitudinal_bond_output_full = true`` under a version writing archive
  schema 2 (see the description of ``longitudinal_bond`` in
  ``file.output`` and :ref:`the archive reference <subsec:chiq_rpa>`); an
  older schema-1 archive is refused with a message asking for a re-run.
  Refused together with ``bond_green`` / ``bond_max_shells`` (the Green
  function and the bond topology come from ``path_to_flex_output`` and the
  archive, not from these static-bond-path keys) and ``zero_chi_s`` /
  ``zero_chi_c`` (not implemented on the bond path). Exactly one
  ``pairing_type`` (``"singlet"`` or ``"triplet"``) is solved per
  ``hwave_sc`` run, as on the on-site dynamic path; a second channel is a
  second run with its own ``path_to_output``. See :ref:`the FLEX tutorial
  section on pairing with the bond-resolved vertex
  <flex_bond_pairing_tutorial>` for the memory model, the sparse-ir
  dependency and worked examples.

- ``flex_bond_archive``

  **Type :**
  String (default value is ``"longitudinal_bond.npz"``; post-processing
  only)

  **Description :**
  Filename of the dynamic bond archive read by the ``bond_channels = true``
  dynamic path. A relative name is joined to ``[file.input]
  path_to_flex_output``; an absolute path is used as given (the same rule
  ``bond_green`` follows).

- ``ir_fit_tol``

  **Type :**
  Float (default value is ``0.5``; both entries)

  **Description :**
  With ``matsubara_basis = "ir"``, refusal threshold of the componentwise
  relative residual of fitting the uniform-grid bond vertex onto the IR
  basis (evaluated separately for the spin and the charge contribution of
  every bond-channel block, then maximized). A residual in
  ``[0.1 * ir_fit_tol, ir_fit_tol)`` warns; ``ir_fit_tol = 0`` skips the
  check entirely (the outputs then record the residual as ``nan``). Real
  uniform-FFT archives typically give a 0.1-0.2 componentwise residual with
  the constant retained (``ir_keep_static_chi = true``); independent of
  ``ir_tol``.

- ``parity_leakage_tol``

  **Type :**
  Float (default: unset -- resolves to ``1e-8`` on the uniform grid and
  ``2e-2`` with ``matsubara_basis = "ir"``; both entries)

  **Description :**
  Refusal threshold of the parity-commutation probe the bond kernel runs
  before every solve (the direct-term-only kernel is the physical one only
  on a definite-parity subspace, so a leakage above this threshold is
  refused rather than silently solved). The IR representation of a
  uniform-FFT archive carries a parity asymmetry of its own, decaying as
  ``Nmat^-2`` (measured 7.6e-3 / 1.3e-3 at Nmat 64 / 128 on a single-band
  test fixture), which is why the IR default is looser than the uniform
  one. A leakage above ``0.1 * parity_leakage_tol`` warns; the measured
  value is recorded in the outputs as ``bond_parity_leakage``.

- ``bond_memory_cap_gb``

  **Type :**
  Float (default: 80% of the measured free host memory; both entries)

  **Description :**
  Host memory cap of the dynamic bond-resolved pairing step, in binary
  GiB (the existing static-bond key of the same name, reused here as the
  host cap of this path; GPU memory is not covered). Admission is checked
  before the vertex is built and again before the kernel's largest arrays;
  the run (or, in-process, the requested channel) is refused with the full
  memory table when even the streaming residency does not fit. In-process,
  admission also runs once BEFORE the first FLEX map so that a long FLEX
  run is not lost to a pairing-step allocation failure after the fact.

**In-process pairing** (``[mode.param] longitudinal_bond_pairing``, see the
``mode.param`` section above): when a FLEX input sets
``longitudinal_bond_pairing`` to something other than ``"none"``, an
``[eliashberg]`` table in the SAME input file is optional and configures
the pairing step at the end of the FLEX solve; every key of the general
table not listed as post-processing-only above (``solver_mode``,
``num_eigenvalues``, ``max_iter``, ``alpha``, ``convergence_tol``,
``matsubara_basis``, ``ir_tol``, ``ir_wmax``, ``ir_keep_static_chi``,
``fft_workers``, ...) defaults exactly as it does for ``hwave_sc``.
``pairing_type`` is refused (the channel(s) are selected by
``longitudinal_bond_pairing`` instead); ``gpu`` / ``gpu_required`` are
ignored (the kernel runs on the FLEX solve's own backend); ``bond_channels``,
``chi0q_mode``, ``frequency``, ``flex_bond_archive``, ``bond_green``,
``bond_max_shells`` are ignored too (an INFO line lists any that are
present -- a shared input file for ``hwave`` then ``hwave_sc`` is a normal
workflow, not a mistake). The ``[file.output]`` keys of the resulting files
are documented under ``file.output`` below.


``log`` section
================================

- ``print_level``

  **Type :**
  Integer (default value is 1)

  **Description :**
  This parameter specifies the verbosity of the standard output.
  If it is set to 1, a detailed information will be shown.


``file`` section
================================

This section consists of ``input`` and ``output`` subsections.
They specify the settings of the input and output files, respectively, on the types of files, the directories to be located or stored, and the names of the files.


``file.input`` section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``path_to_input``

  **Type :**
  String (default value is ``""`` (blank string))

  **Description :**
  This parameter specifies the directory in which the input files are located.

- ``chi0q_init``

  **Type :**
  String

  **Description :**
  This parameter specifies the filename of the pre-calculated irreducible susceptibility
  :math:`\chi_0(\vec{q})` to be used for the calculation of the susceptibility matrix.
  The input file is in NumPy binary format that corresponds to the output format of
  ``chi0q`` in ``file.output`` section.

- ``trans_mod``

  **Type :**
  String

  **Description :**
  This parameter specifies the filename of the initial configuration exported from UHFk by the parameter ``file.output.rpa``. It contains the one-body interaction term involving the approximated two-body interaction terms via UHF method.

- ``green_init``

  **Type :**
  String

  **Description :**
  This parameter specifies the filename of the initial Green's function for RPA calculation. The file format corresponds to the output file of ``green`` of UHFk. When ``trans_mod`` is specified, ``green_init`` is not used.

- ``sigma_init``

  **Type :**
  String (FLEX mode only)

  **Description :**
  This parameter specifies the filename of a ``sigma.npz`` (written by an
  earlier FLEX run) used to seed the FLEX SCF loop instead of Sigma = 0. The
  path is resolved relative to ``path_to_input``. The recorded ``CellShape``
  and the array's ``Nmat`` must match the current run (a mismatch is a
  fail-fast error). See the FLEX tutorial section "Warm-starting the SCF
  loop" for the sweep workflow. With ``flex_hartree_fock = true`` only a
  ``sigma.npz`` written by such a run (``sigma_convention = "split"``,
  carrying ``sigma_static`` and ``sigma_fluct``) is accepted, the two
  components seed the loop exactly as stored, and a mean field
  (``trans_mod`` / ``green_init``) given at the same time is refused; a
  legacy or ``"total"`` archive is converted with ``hwave_sigma_split``
  (see :ref:`flex_bond_hf`).


``file.input.interaction`` section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This section describes the relation of the interaction types and geometry information
to the definition files.

- ``path_to_input``

  **Type :**
  String

  **Description :**
  This parameter specifies the directory in which the input files are located.
  It is independent from ``path_to_input`` in ``file.input`` section.

- ``Geometry``

  **Type :**
  String

  **Description :**
  This parameter specifies the filename for the geometry information.

- ``Transfer``, ``CoulombIntra``, ``CoulombInter``, ``Hund``, ``Ising``, ``Exchange``, ``PairLift``, ``PairHop``, ``Extern``

  **Type :**
  String

  **Description :**
  These parameters specify the filenames for the definitions of the corresponding interaction terms. If none of two-body interaction term (CoulombIntra, CoulombInter, Hund, Ising, Exchange, PairLift, or PairHop) is specified, the program only calculates ``chi0q`` and exits.


``file.output`` section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``path_to_output``

  **Type :**
  String (default value is ``"output"``)

  **Description :**
  This parameter specifies the directory in which the output files are stored.

- ``chi0q``

  **Type :**
  String

  **Description :**

  This parameter specifies the name of the file to store the irreducible susceptibility matrix
  :math:`\chi_0(\vec{q})`.
  If it is not set, no output will be generated.

- ``chiq``

  **Type :**
  String

  **Description :**
  This parameter specifies the name of the file to store the susceptibility matrix
  :math:`\chi(\vec{q})`.
  If it is not set, no output will be generated.

- ``longitudinal_bond``

  **Type :**
  String (default value is ``longitudinal_bond.npz``; FLEX mode only)

  **Description :**
  Name of the dedicated archive of the dynamic bond-resolved
  susceptibilities, written only when ``longitudinal_bond_channels =
  true`` and ``longitudinal_bond_output_full = true`` in FLEX mode (see
  :ref:`flex_bond_hf`). A name without the ``.npz`` suffix gets it
  appended, as for the other ``.npz`` outputs. It must not resolve to the
  same file as any other output of the run (the collision is refused
  before the calculation starts).

- ``eliashberg_bond_singlet``, ``eliashberg_bond_triplet``,
  ``eigenvalue_bond_singlet``, ``eigenvalue_bond_triplet``,
  ``gap_bond_singlet``, ``gap_bond_triplet``

  **Type :**
  String (defaults ``eliashberg_bond_singlet.npz``,
  ``eliashberg_bond_triplet.npz``, ``eigenvalue_bond_singlet.dat``,
  ``eigenvalue_bond_triplet.dat``, ``gap_bond_singlet.dat``,
  ``gap_bond_triplet.dat``; FLEX mode only)

  **Description :**
  Output filenames of the in-process bond-resolved pairing step
  (``[mode.param] longitudinal_bond_pairing``, see the ``mode.param``
  section above): the frequency-resolved gap archive, the leading-eigenvalue
  file and the single-frequency gap slice of each requested channel (see
  :ref:`the output reference <subsec:eliashberg_bond_outputs>`). Archive
  names without a ``.npz`` suffix get one appended, as for the other
  ``.npz`` outputs. Each name must not resolve to the same file as any
  other output of the run (the collision is refused before the calculation
  starts).
