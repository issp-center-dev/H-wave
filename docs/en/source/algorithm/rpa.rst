.. highlight:: none

.. _algorithm_sec:

Random Phase Approximation
==========================

The random phase approximation (RPA) is a method to detect the response to the fluctuations of
one-body operators by the effect of electron correlations, starting from the non-interacting
state.
While in the UHF approximation an initial guess of the configuration is required, 
the RPA method enables to infer the ordered phase that emerges from the second-order transition.
H-wave implements RPA method using Matsubara frequency, and allows to compare with the dynamical
observables measured in the experiments by analytical continuation.

In the following, the algorithm is described.
In the RPA mode of H-wave, the Hamiltonian given below will be considered:

.. math::
    \begin{aligned}
     {\cal H}&={\cal H}_0+{\cal H}_{\rm int},\\
     {\cal H}_0&=\sum_{\langle i\alpha;j\beta \rangle}
      (t_{ij}^{\alpha \beta}c_{i\alpha}^{\dagger}
      c_{j\beta}^{\mathstrut}+\mbox{H.c.}),\\
     {\cal H}_{\rm int}&=\sum_{ij}\sum_{\alpha, \alpha', \beta, \beta'}W_{ij}^{\beta\beta',\alpha\alpha'}\left(
      c_{i\alpha}^{\dagger}c_{i\alpha'}c_{j\beta'}^{\dagger}c_{j\beta}+\mbox{H.c.}\right)
    \end{aligned}

Applying the Fourier transformation

.. math::
    \begin{aligned}
    c_{i\alpha}
    =\frac{1}{\sqrt{N_L}}\sum_{\bf{k}}
    e^{i \bf{k}\cdot \bf{r}_{i}}c_{\bf{k},\alpha}^{\mathstrut},
    \end{aligned}

the Hamiltonian is rewritten in the following form

.. math::
    \begin{aligned}
     {\cal H}&=\sum_{{\bf k}\alpha\beta}
     (\varepsilon_{\alpha\beta}({\bf k})c_{{\bf k}\alpha}^{\dagger}
     c_{{\bf k}\beta}^{\mathstrut}+\mbox{H.c.}) \nonumber\\
    &+\frac{1}{2N_L}\sum_{{\bf k} {\bf k}'{\bf q}}\sum_{\alpha\beta\alpha'\beta'}
     W^{\beta\beta',\alpha\alpha'}_{{\bf q}}
     c_{{\bf k}+{\bf q},\alpha}^{\dagger}
      c_{{\bf k},\alpha'}^{\mathstrut}
      c_{{\bf k}'-{\bf q},\beta'}^{\dagger}
      c_{{\bf k}',\beta}^{\mathstrut}
    \end{aligned}

In the random phase approximation, the density fluctuation by the effect of electron correlation
is detected with respect to :math:`{\cal H}_0`.
The scattering by the interaction must therefore be considered on the basis
where :math:`{\cal H}_0` is diagonalize, and thus the interaction term is approximated as

.. math::
    \begin{aligned}
    &W^{\beta\beta',\alpha\alpha'}_{\bf{q}}c_{\bf{k}+\bf{q},\alpha}^{\dagger}c_{\bf{k},\alpha'}^{\mathstrut}
    c_{\bf{k}'-\bf{q},\beta'}^{\dagger} c_{\bf{k}',\beta}^{\mathstrut}\nonumber\\
    &\sim W^{\beta\beta',\alpha\alpha'}_{\bf{q}} \sum_{\gamma, \gamma'}
    u_{\alpha \gamma, \bf{k}+\bf{q}}^* d_{\bf{k}+\bf{q},\gamma}^{\dagger}
    u_{\alpha' \gamma, \bf{k}} d_{\bf{k},\gamma}^{\mathstrut}
    u_{\beta' \gamma', \bf{k}'-\bf{q}}^* d_{\bf{k}'-\bf{q},\gamma'}^{\dagger}
    u_{\beta  \gamma', \bf{k}'}d_{\bf{k}',\gamma'}^{\mathstrut}.
    \end{aligned}

Here, 

.. math::
    \begin{aligned}
    c_{\bf{k},\alpha} = \sum_{\gamma} u_{\alpha \gamma, \bf{k}} d_{\bf{k}, \gamma}
    \end{aligned}

and :math:`d_{\bf{k}, \gamma}` denotes the annihilation operator that diagonalizes :math:`{\cal H}_0`. (:math:`\gamma` refers to the index of the eigenvalue.)
Then, the irreducible one-body Green's function is written as

.. math::
    \begin{aligned}
     G^{(0)\alpha\beta}_{\gamma}({\bf k}, i\omega_{n})=
      \frac{u^{\alpha\gamma}({\bf k})u^{*\beta\gamma}({\bf k})}{i\epsilon_{n}-\xi^{\gamma}({\bf k})+\mu}.
    \end{aligned}

The irreducible susceptibility is given in the following form, as it must be closed
within the diagnoalized elements:

.. math::
    \begin{aligned}
     X^{(0)\alpha\alpha', \beta\beta'}({\bf q},i\omega_n)=
      -\frac{T}{N_L}
      \sum_{\gamma=1}^{n_{\rm orb}}\sum_{{\bf k},n}
      G^{(0)\alpha\beta}_{\gamma}({\bf k}+{\bf q}, i\omega_m+ i\epsilon_{n})
      G^{(0)\beta'\alpha'}_{\gamma}({\bf k}, i\epsilon_{n}),
    \end{aligned}

By using the irreducible susceptibility, the susceptibility matrix from the RPA
is obtained as follows:

.. math::
    \begin{aligned}
    X^{\alpha\alpha', \beta\beta'}(q)&=
    X^{(0)\alpha\alpha', \beta\beta'}(q) - \sum_{\alpha_1'\beta_1'}
    X^{(0)\alpha\alpha', \beta_1\beta_1'}(q) W^{\beta_1\beta_1', \alpha_1\alpha_1'}_{\bf q}X^{\alpha_1 \alpha_1' , \beta \beta'}(q),
    \end{aligned}

Combining indices such as :math:`\alpha\alpha^\prime` into one index, they are expressed
in the matrix form. Then finally it leads to the expression:

.. math::
    \begin{aligned}
     \hat{X}(q)&=\hat{X}^{(0)}(q)-\hat{X}^{(0)}(q)\hat{W}(q)\hat{X}(q)\nonumber\\
     &=\left[\hat{I}+\hat{X}^{(0)}(q)\hat{W}(q)\right]^{-1}\hat{X}^{(0)}(q).
    \end{aligned}

It is noted that the vertex correction may be taken into account as a means to consider
higher order correlations. See, for example, reference [1]_ for the details.

.. [1] `K. Yoshimi, T. Kato, H. Maebashi, J. Phys. Soc. Jpn. 78, 104002 (2009). <https://journals.jps.jp/doi/10.1143/JPSJ.78.104002>`_
