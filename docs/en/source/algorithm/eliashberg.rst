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
and charge (:math:`W_c`) channels:

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
     - :math:`U' - I`
     - :math:`-U' + J - I`
   * - :math:`l_1 = l_2 \neq l_3 = l_4`
     - Density
     - :math:`J - 2I`
     - :math:`2U' - J`
   * - :math:`l_1 = l_4 \neq l_2 = l_3`
     - Exchange
     - :math:`J' + P`
     - :math:`J' + P`

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
