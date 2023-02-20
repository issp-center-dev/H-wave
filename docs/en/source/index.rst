#####################################################
Welcome to H-wave documentations!
#####################################################

What is H-wave?
----------------------------------------------------------------

H-wave is a program to perform unrestricted Hartree-Fock (UHF) approximation for itinerant electron systems.
The UHF approximation reduces two-body interaction terms into one-body terms by taking account of the fluctuation :math:`\delta A \equiv A-\langle A \rangle` up to first order. 

Let us consider the inter-site Coulomb interaction for example. 

.. math::

   {\cal H}_V = \sum_{i,j, \sigma, \sigma'}V_{ij} n_ {i\sigma}n_{j\sigma'}

We adopt a notation :math:`i\equiv (i, \sigma)` and :math:`j\equiv (j, \sigma')` for brevity.
The interaction term is approximated to the following form by eliminating the second order contribution of the fluctuation:

.. math::

   \begin{aligned}
   n_ {i}n_{j} &=
   (\langle n_{i} \rangle +\delta n_i) (\langle n_{j} \rangle +\delta n_j)
   - \left[ \langle c_{i}^{\dagger}c_j \rangle +\delta (c_{i}^{\dagger}c_j ) \right]
     \left[ \langle c_{j}^{\dagger}c_i \rangle +\delta (c_{j}^{\dagger}c_i )\right]
   \nonumber\\
   &\sim
   \langle n_{i} \rangle n_j+\langle n_{j} \rangle  n_i
   - \langle c_{i}^{\dagger}c_j \rangle  c_{j}^{\dagger}c_i  -  \langle c_{j}^{\dagger}c_i \rangle c_{i}^{\dagger}c_j 
   -\langle n_{i} \rangle \langle n_j \rangle +  \langle c_{j}^{\dagger}c_i \rangle \langle c_{i}^{\dagger}c_j \rangle
   \end{aligned}

Other types of interactions are reduced to one-body problems in a similar manner. 
In the actual calculation, the above equations are solved iteratively until the expectaion values become self-consistent.

License
----------------------------------------------------------------

The distribution of the program package and the source codes for H-wave follow GNU General Public License version 3 (GPL v3) or later.

Contributors
----------------------------------------------------------------

This software was developed by the following contributors.

-  ver.1.0 (released on 2023/xx/xx)

   -  Developers

      -  Kazuyoshi Yoshimi
	 (The Institute for Solid State Physics, The University of Tokyo)

      -  Tatsumi Aoyama
	 (The Institute for Solid State Physics, The University of Tokyo)

      -  Yuichi Motoyama
	 (The Institute for Solid State Physics, The University of Tokyo)

      -  Takahiro Misawa
	 (Beijing Academy of Quantum Information Sciences (BAQIS))

      -  Kota Ido
	 (The Institute for Solid State Physics, The University of Tokyo)

      -  Akito Kobayashi
	 (Department of Physics, Nagoya University)

      -  Taiki Kawamura
	 (Department of Physics, Nagoya University)

   -  Project coordinator

      -  Takeo Kato
	 (The Institute for Solid State Physics, The University of Tokyo)


Copyright
----------------------------------------------------------------

© *2022- The University of Tokyo. All rights reserved.*

This software was developed with the support of "Project for advancement of software usability in materials science" of The Institute for Solid State Physics, The University of Tokyo.


Operating environment
----------------------------------------------------------------

H-wave was tested on the following platforms

- macOS + python3 (brew)

- Ubuntu Linux + python3 (miniconda)


Contents
----------------------------------------------------------------

.. toctree::
   :maxdepth: 3
   :numbered: 3

   howtouse/ho-index
   uhfr/uhfr-index
   uhfk/uhfk-index
   rpa/rpa-index
   algorithm/al-index
   acknowledgement
   

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
