.. highlight:: none

.. _algorithm_sec:

Unrestricted Hartree-Fock method
================================

Overview
*****************************

The unrestricted Hartree-Fock approximation is a method to approximate
two-body interactions into one-body terms by taking account of the fluctuation
of the one-body operators up to first order.
For a general two-body interactions, it leads to the following approximation:

.. math::
   \begin{aligned}
   c_{i}^{\dagger} c_{j}^{\dagger} c_{k} c_{l} 
   & \sim
   \langle c_{i}^{\dagger} c_l \rangle c_{j}^{\dagger} c_k
   + c_{i}^{\dagger} c_l \langle c_{j}^{\dagger} c_k\rangle
   - \langle c_{i}^{\dagger} c_k \rangle c_{j}^{\dagger} c_l
   - c_{i}^{\dagger} c_k \langle c_{j}^{\dagger} c_l \rangle \nonumber \\
   & \qquad
   - (\langle c_{i}^{\dagger} c_l \rangle \langle c_{j}^{\dagger} c_k \rangle
   - \langle c_{i}^{\dagger} c_k\rangle \langle c_{j}^{\dagger} c_l\rangle) .
   \end{aligned}

In H-wave, the two-body interaction terms are defined as
   
.. math::
   \begin{aligned}
   \mathcal{H}_\text{InterAll}
   &=
   \sum_{ijkl\alpha\beta\gamma\delta} \sum_{\sigma_1 \sigma_2 \sigma_3 \sigma_4}
   I_{ijkl\alpha\beta\gamma\delta}
   c^\dagger_{i\alpha\sigma_1} c_{j\beta\sigma_2} c^\dagger_{k\gamma\sigma_3} c_{l\delta\sigma_4}
   \nonumber\\
   &=
   \sum_{ijkl\alpha\beta\gamma\delta} \sum_{\sigma_1 \sigma_2 \sigma_3 \sigma_4}
   I_{ijkl\alpha\beta\gamma\delta} (
   c^\dagger_{i\alpha\sigma_1} c^\dagger_{k\gamma\sigma_3} c_{l\gamma\sigma_4} c_{j\beta\sigma_2}
   + c^\dagger_{i\alpha\sigma_1} c_{l\delta\sigma_4} \delta_{j,k}\delta_{\beta,\gamma}\delta_{\sigma_2,\sigma_3} ) .
   \end{aligned}

It is noted that there is a one-body term as depicted in the second term of
the above expression.
Then, the Hamiltonian given by the one-body terms is generally denoted as

.. math::
   \begin{aligned}
   \mathcal{H}_\text{UHF} &= \sum_{ij} H_{ij} c^\dagger_{i} c_{j} = \hat{c}^\dagger H \hat{c}
   \end{aligned}

where we adopt a notation :math:`i\equiv(i, \alpha, \sigma_1), j\equiv(j, \beta, \sigma_2)` for brevity,
:math:`H` denotes a matrix whose elements are :math:`H_{ij}`, and 
:math:`\hat{c}` denotes a column vector whose elements are :math:`c_{i}`.

As :math:`H` is an Hermite matrix,
the Hamiltonian can be transformed into :math:`H=U \hat{\xi} U^\dagger`
where
:math:`\hat{\xi}` is a matrix whose diagonal elements are the eigenvalues of :math:`H`, 
and :math:`U` is a matrix composed of the corresponding eigenvectors.

Then, let :math:`\hat{d} = U^\dagger \hat{c}`, and :math:`\mathcal{H}_\text{UHF}` leads to

.. math::
   \begin{aligned}
   \mathcal{H}_\text{UHF} &= \hat{d}^\dagger \hat{\xi} \hat{d} =  \sum_{k} \xi_k d_k^\dagger d_k .
   \end{aligned}
   
Therefore, the energy derived from the one-body interaction term of the UHF approximation
is obtained by

.. math::
   \begin{aligned}
   E_\text{UHF} = \langle \mathcal{H}_\text{UHF} \rangle = \sum_{k} \xi_k \langle d_k^\dagger d_k \rangle .
   \end{aligned}

In the numerical calculation, as :math:`H` depends on the one-body Green's function
:math:`\langle c_{i}^\dagger c_{j}\rangle`
through the UHF approximation, the equation is iteratively solved to satisfy the self-consistency.
Starting from a one-body Green's function given as an initial value, it is updated through the relation

.. math::
   \begin{aligned}
   \langle c_{i}^\dagger c_{j}\rangle = \sum_{l} U_{il}^* U_{jl} \langle d_l^\dagger d_l \rangle = \sum_{l} \frac{U_{il}^* U_{jl}}{1+\exp^{\beta(\xi_l -\mu)}}
   \end{aligned}

until the one-body Green's function converges. 
Here, :math:`\beta` denotes the inverse temperature :math:`1/ k_B T`, and
:math:`\mu` denotes the chemical potential.
In the canonical calculation in which the number of particles is fixed,
:math:`\mu` is determined to satisfy the relation

.. math::
   \begin{aligned}
   N = \sum_{i} \langle c_i^{\dagger} c_i \rangle
   \end{aligned}

for the number of particles :math:`N` at every step.

In H-wave, the simple-mixing algorithm is employed to update the configuration.
If we denote the one-body Green's function at *n*-th step by 
:math:`\langle c_{i}^\dagger c_{j}\rangle^{(n)}`,
the Green's function at *n+1*-th step is chosen by mixing that of *n*-th step
with the new one obtained in *n+1*-th step as

.. math::
   \begin{aligned}
   \langle c_{i}^\dagger c_{j}\rangle^{(n+1)} := (1-\alpha) \langle c_{i}^\dagger c_{j}\rangle^{(n)} +  \alpha \langle c_{i}^\dagger c_{j}\rangle^{(n+1)}, 
   \end{aligned}

where :math:`\alpha` is a parameter between 0 and 1.
There are other update algorithms such as Anderson mixing, though they are not supported
in the present version of H-wave.

In the coordinate-space UHF mode of H-wave, all interactions are mapped to InterAll form.
The free energy at finite temperature is given by

.. math::
   \begin{aligned}
   F = \mu N -\frac{1}{\beta}\sum_k \ln \left[ 1+\exp (-\beta(\xi_k - \mu)) \right]
    - \sum_{ijkl} I_{ijkl} (\langle c_{i}^{\dagger} c_j\rangle \langle c_{k}^{\dagger} c_l\rangle - \langle c_{i}^{\dagger} c_l\rangle \langle c_{k}^{\dagger} c_j\rangle) .
   \end{aligned}


Extension to wave-number space
******************************

The Hamiltonian given by the one-body terms is rewritten in the wave-number representation
by the Fourier transform :math:`c_i = \dfrac{1}{\sqrt{V}} \sum_k e^{ikr_i} c_k` as

.. math::
   \begin{aligned}
   \mathcal{H}_\text{UHF}
   &=
   \sum_{k\alpha\beta\sigma\sigma^\prime}
   h_{\alpha\beta\sigma\sigma^\prime}(k)\,
   c_{k\alpha\sigma}^\dagger c_{k\beta\sigma^\prime}^{\phantom\dagger}
   \end{aligned}

Here, the interaction is assumed to have translational symmetry so that 
the coefficients depend only on the translation vectors :math:`r_{ij}=r_j - r_i`.
It is noted that InterAll type of interaction is not considered in the
wave-number space UHF mode.

As the Hamiltonian is diagonal with respect to the wave number :math:`k`,
the calculation of the eigenvalues and eigenvectors reduces from
diagonalization of a matrix of the size
:math:`N_\text{site}N_\text{orbit} \times N_\text{site}N_\text{orbit}`
to that of :math:`N_\text{site}` matrices of the size
:math:`N_\text{orbit} \times N_\text{orbit}`,
which lowers the calculation costs.
Here, 
:math:`N_\text{site}` denotes the number of sites, and
:math:`N_\text{orbit}` denotes the number of orbitals including the spin degree of freedom.      

