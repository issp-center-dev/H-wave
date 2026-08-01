.. highlight:: none

.. _algorithm_eliashberg:

Linearized Eliashberg equation
================================

Overview
*****************************

The linearized Eliashberg equation solver (``hwave_sc``) analyzes
superconducting instabilities by solving the eigenvalue problem
of the linearized gap equation within the RPA framework.
The superconducting transition temperature :math:`T_c` is determined
by the condition that the leading eigenvalue reaches :math:`\lambda = 1`.

The algorithm proceeds in the following steps:

1. Compute or load the bare susceptibility :math:`\hat{X}^{(0)}(\mathbf{q})`.
2. Construct the non-interacting Green's function :math:`G(\mathbf{k}, i\omega_n)`.
3. Build the RPA pairing vertex :math:`V(\mathbf{q})`.
4. Solve the linearized Eliashberg equation for the leading eigenvalue.


Green's function
*****************************

The non-interacting one-body Hamiltonian in k-space is obtained
by Fourier transformation of the hopping integrals:

.. math::

   \varepsilon_{\alpha\beta}(\mathbf{k})
   = \sum_{\mathbf{R}} t^{\alpha\beta}_{\mathbf{R}}\,
     e^{i\mathbf{k}\cdot\mathbf{R}}

Diagonalizing the Hamiltonian at each k-point:

.. math::

   \sum_\beta \varepsilon_{\alpha\beta}(\mathbf{k})\, u_{\beta m}(\mathbf{k})
   = \xi_m(\mathbf{k})\, u_{\alpha m}(\mathbf{k})

the non-interacting Green's function is given by

.. math::

   G_{\alpha\beta}(\mathbf{k}, i\omega_n)
   = \sum_{m} \frac{u_{\alpha m}(\mathbf{k})\, u^*_{\beta m}(\mathbf{k})}
                    {i\omega_n - (\xi_m(\mathbf{k}) - \mu)}

where :math:`\omega_n = \pi(2n+1)/\beta` are the fermionic Matsubara frequencies,
:math:`\mu` is the chemical potential determined by the filling condition,
and :math:`\beta = 1/T` is the inverse temperature.


Pairing vertex
*****************************

Simple mode
-----------------------------

When only ``CoulombIntra`` (:math:`U`) and ``CoulombInter`` (:math:`V`)
are present, the pairing vertex is computed using the spin (:math:`W_s`)
and charge (:math:`W_c`) channels. :math:`U` and :math:`V` enter as the
symmetrised reading of the declaration file -- the mean with the
reversed-bond partner :math:`(R, a, b) \leftrightarrow (-R, b, a)` --
matching every other route in the package: a one-sided off-site
declaration therefore contributes :math:`v \cos(qR)`, not
:math:`v e^{-iqR}` (earlier versions of the simple mode used the raw
one-sided phase; declarations that already contain both directions are
read unchanged, bit for bit):

.. math::

   W_s = -U, \qquad W_c = U + 2V

The RPA susceptibilities are

.. math::

   \hat{X}^s(\mathbf{q}) = \left[\hat{I} - \hat{X}^{(0)}(\mathbf{q})\, \hat{W}_s\right]^{-1} \hat{X}^{(0)}(\mathbf{q})

.. math::

   \hat{X}^c(\mathbf{q}) = \left[\hat{I} + \hat{X}^{(0)}(\mathbf{q})\, \hat{W}_c\right]^{-1} \hat{X}^{(0)}(\mathbf{q})

The singlet pairing vertex is

.. math::

   V^S_{\alpha\beta}(\mathbf{q})
   = \frac{1}{2}(W_c + W_s)_{\alpha\beta}
     + \frac{3}{2} (W_s\, X^s\, W_s)_{\alpha\beta}
     - \frac{1}{2} (W_c\, X^c\, W_c)_{\alpha\beta}

and the triplet pairing vertex is

.. math::

   V^T_{\alpha\beta}(\mathbf{q})
   = \frac{1}{2}(W_c - W_s)_{\alpha\beta}
     - \frac{1}{2} (W_s\, X^s\, W_s)_{\alpha\beta}
     - \frac{1}{2} (W_c\, X^c\, W_c)_{\alpha\beta}


General mode (S/C matrix formulation)
-----------------------------------------

When ``Hund`` (:math:`J`), ``Exchange`` (:math:`J'`),
``Ising`` (:math:`I`), or ``PairHop`` (:math:`P`) interactions
are present, the solver uses the generalized :math:`S`/:math:`C` matrix
formulation [1]_.

The :math:`S` and :math:`C` matrices are defined in the composite index
space :math:`(l_1, l_2)` where :math:`l_1, l_2` are orbital indices.
For a system with :math:`n_{\rm orb}` orbitals, the matrices have
dimension :math:`n_{\rm orb}^2 \times n_{\rm orb}^2`.
The matrix elements are:

.. list-table::
   :header-rows: 1
   :widths: 30 20 25 25

   * - Index condition
     - Type
     - :math:`S` value
     - :math:`C` value
   * - :math:`l_1 = l_2 = l_3 = l_4`
     - Intra-orbital
     - :math:`U`
     - :math:`U`
   * - :math:`l_1 = l_3 \neq l_2 = l_4`
     - Cross
     - :math:`U' + I - J + J'`
     - :math:`-U' - I + J + J'`
   * - :math:`l_1 = l_2 \neq l_3 = l_4`
     - Density
     - :math:`J - 2I`
     - :math:`2U' - J`
   * - :math:`l_1 = l_4 \neq l_2 = l_3`
     - Pair hop
     - :math:`P`
     - :math:`P`

.. note::

   The Cross-row entries for :math:`J`, :math:`J'` and :math:`I`, and the
   placement of :math:`J'` in the Cross row rather than the last row, follow
   the exact-diagonalization adjudication of the interaction vertices: the
   Hund term contributes :math:`S \mathrel{-}= J`, :math:`C \mathrel{+}= J`,
   the Exchange term :math:`S \mathrel{+}= J'`, :math:`C \mathrel{+}= J'`,
   and the Ising term enters :math:`S` with :math:`+I`. For the
   SU(2)-symmetric Kanamori combination (:math:`J' = J`) these reproduce the
   standard literature matrices (no :math:`J` in the Cross :math:`S` entry and
   :math:`-U' + 2J` in the Cross :math:`C` entry). Susceptibility files
   produced before this correction carry no ``sc_vertex_version`` field and
   are rejected when the interaction set contains Hund, Exchange or Ising.
   For general-scheme (``myo``) files, an unversioned file is also rejected
   when CoulombInter, Hund, Ising or Exchange declares an asymmetric
   on-site inter-orbital coupling (:math:`X_{ab} \neq X_{ba}`): the
   orbital orientation in which the interaction enters the vertex and
   the symmetrised reading of such declarations both differ between
   historical builds and the current one, and an unversioned file does
   not record which semantics it was produced with. PairHop is checked
   against its HERMITIAN partner instead of the plain transpose (its two
   declarations are Hermitian partners, so its orientation never changed
   -- but an on-site declaration that is not Hermitian-closed is also
   rejected, because the conjugated-mean reading of PairHop declarations
   arrived together with the version stamp). PairLift is exempt (its
   particle-hole vertex contribution is exactly zero, so neither the
   orientation nor the symmetrised reading can matter).
   The unaffected cases, stated per rule: for the four
   transpose-checked interactions, transpose-symmetric declarations
   (:math:`X_{ab} = X_{ba}`) -- every physically ordinary input; for
   PairHop, Hermitian-closed declarations
   (:math:`P_{ba} = P_{ab}^{*}`) -- the physically valid form, which
   need not be transpose-symmetric. Reduced-scheme files never depended
   on this orientation.

The RPA susceptibilities are

.. math::

   \hat{X}^s(\mathbf{q}) = \left[\hat{I} - \hat{X}^{(0)}(\mathbf{q})\, \hat{S}\right]^{-1} \hat{X}^{(0)}(\mathbf{q})

.. math::

   \hat{X}^c(\mathbf{q}) = \left[\hat{I} + \hat{X}^{(0)}(\mathbf{q})\, \hat{C}\right]^{-1} \hat{X}^{(0)}(\mathbf{q})

The singlet pairing vertex is

.. math::

   \hat{V}^S(\mathbf{q})
   = \frac{3}{2}\, \hat{S}\, \hat{X}^s(\mathbf{q})\, \hat{S}
     - \frac{1}{2}\, \hat{C}\, \hat{X}^c(\mathbf{q})\, \hat{C}
     + \frac{1}{2}(\hat{S} + \hat{C})

and the triplet pairing vertex is

.. math::

   \hat{V}^T(\mathbf{q})
   = -\frac{1}{2}\, \hat{S}\, \hat{X}^s(\mathbf{q})\, \hat{S}
     - \frac{1}{2}\, \hat{C}\, \hat{X}^c(\mathbf{q})\, \hat{C}
     + \frac{1}{2}(\hat{C} - \hat{S})


Linearized Eliashberg equation
***********************************

The linearized Eliashberg equation is formulated as an eigenvalue problem:

.. math::

   \lambda\, \Sigma_{\alpha\beta}(\mathbf{k})
   = -\frac{T}{N_L} \sum_{\mathbf{k}', n', \alpha', \beta'}
     V_{\alpha\alpha';\beta\beta'}(\mathbf{k} - \mathbf{k}')
     \, G_{\alpha\alpha'}(\mathbf{k}', i\omega_{n'})
     \, G_{\beta\beta'}(-\mathbf{k}', -i\omega_{n'})
     \, \Sigma_{\alpha'\beta'}(\mathbf{k}')

where :math:`\Sigma_{\alpha\beta}(\mathbf{k})` is the anomalous self-energy
(gap function), and the right-hand side defines the Eliashberg kernel
:math:`K[\Sigma]`.

The superconducting instability occurs when :math:`\lambda = 1`.
Only positive eigenvalues are physically relevant for the SC transition.


Two-particle Green's function
-----------------------------

The Matsubara frequency summation is performed analytically to obtain
the two-particle Green's function:

.. math::

   G^{(2)}_{\alpha\beta;\gamma\delta}(\mathbf{q})
   = \frac{T}{N_L} \sum_{\mathbf{k}, n}
     G_{\alpha\gamma}(\mathbf{k}, i\omega_n)\,
     G_{\beta\delta}(-\mathbf{k}+\mathbf{q}, -i\omega_n)

This reduces the Eliashberg kernel to a convolution in k-space,
which is efficiently computed using the Fast Fourier Transform (FFT).


FFT-based kernel evaluation
-----------------------------

The kernel operation :math:`\Sigma_{\rm new} = K[\Sigma_{\rm old}]` involves
a convolution of the pairing vertex :math:`V(\mathbf{q})` with the product
of :math:`G^{(2)}` and the gap function.
This convolution is computed efficiently as:

1. Inverse FFT of :math:`V(\mathbf{q})` and :math:`G^{(2)}(\mathbf{q}) \cdot \Sigma(\mathbf{q})` to real space.
2. Pointwise multiplication in real space.
3. FFT back to k-space.

This reduces the computational cost from :math:`O(N_k^2)` to :math:`O(N_k \log N_k)`.


Numerical methods
*****************************

Self-consistent power iteration
---------------------------------

The self-consistent power iteration converges to the eigenmode
with the largest positive eigenvalue:

.. math::

   \Sigma^{(i+1)} = K[\Sigma^{(i)}], \qquad
   \lambda^{(i)} = \|\Sigma^{(i+1)}\|, \qquad
   \Sigma^{(i+1)} \leftarrow \Sigma^{(i+1)} / \lambda^{(i)}

Linear mixing is applied to stabilize convergence:

.. math::

   \Sigma^{(i+1)} \leftarrow (1-\alpha)\, \Sigma^{(i+1)}_{\rm norm}
   + \alpha\, \Sigma^{(i)}

where :math:`\alpha` is the mixing parameter.
The iteration converges when
:math:`\|\Sigma^{(i+1)} - \Sigma^{(i)}\| < \epsilon`.

The initial gap function can be set to various symmetries
(s-wave, :math:`d_{x^2-y^2}`, :math:`\cos k_x + \cos k_y`, random, etc.)
to target specific pairing channels.


Arnoldi eigenvalue analysis
---------------------------------

The Arnoldi method (implicitly restarted, via ARPACK) finds
the leading eigenvalues of the Eliashberg kernel as a linear operator.
This method efficiently computes multiple eigenvalues simultaneously
without requiring explicit construction of the kernel matrix.


Subspace iteration
---------------------------------

The subspace iteration method propagates multiple vectors simultaneously:

1. Apply the kernel to all vectors: :math:`W = K \cdot V`
2. Compute the Rayleigh quotient: :math:`H = V^T K V`
3. Eigendecompose the small matrix :math:`H`
4. Update the subspace with Ritz vectors
5. Re-orthogonalize via QR decomposition

This method is more robust for degenerate eigenvalues.


Shift-invert method
---------------------------------

The shift-invert transformation :math:`(K - \sigma I)^{-1}`
is used to find eigenvalues near a target value :math:`\sigma`.
The linear system is solved iteratively using BiCGSTAB, GMRES,
or LGMRES.


.. [1] K. Kuroki, S. Onari, R. Arita, H. Usui, Y. Tanaka, H. Kontani,
   and H. Aoki, Phys. Rev. Lett. **101**, 087004 (2008);
   K. Kuroki, H. Usui, S. Onari, R. Arita, and H. Aoki,
   Phys. Rev. B **79**, 224511 (2009).
