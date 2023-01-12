.. highlight:: none

.. _Ch:HowToExpert:

Input files for UHFr
--------------------------------

In this section, the input files (\*.def) used in H-wave are described.
They are divided into two types.

(1) Hamiltonian
    
    The Hamiltonian is given in the form of interactions of the electron system.
    They are defined in these files.
    
    | **Trans**
      gives the one-body part expressed by
      :math:`c_{i\sigma_1}^{\dagger}c_{j\sigma_2}^{\phantom{\dagger}}`.
    | **InterAll**
      gives the generalized two-body interactions expressed by
      :math:`c_ {i \sigma_1}^{\dagger} c_{j\sigma_2}^{\phantom{\dagger}} c_{k \sigma_3}^{\dagger} c_{l \sigma_4}^{\phantom{\dagger}}`.
	    
    | Besides, we can specify by the following keywords the interactions that are frequently used.

    | **CoulombIntra**
      gives the on-site Coulomb interaction expressed by
      :math:`n_{i \uparrow} n_{i \downarrow}`, where :math:`n_{i \sigma} = c_{i\sigma}^{\dagger} c_{i\sigma}^{\phantom{\dagger}}`.
    | **CoulombInter**
      gives the inter-site Coulomb interaction expressed by
      :math:`n_{i} n_{j}`, where :math:`n_i = n_{i\uparrow} + n_{i\downarrow}`.
    | **Hund**
      gives the Hund interaction expressed by 
      :math:`n_{i\uparrow} n_{j\uparrow} + n_{i\downarrow} n_{j\downarrow}`.
    | **PairHop**
      gives the pair-hop interaction expressed by
      :math:`c_{i \uparrow}^{\dagger} c_{j\uparrow}^{\phantom{\dagger}} c_{i \downarrow}^{\dagger}c_{j\downarrow}^{\phantom{\dagger}}`.
    | **Exchange**
      gives the exchange interaction expressed by
      :math:`S_i^+ S_j^-`.
    | **Ising**
      gives the Ising interaction expressed by
      :math:`S_i^z S_j^z`.
    | **PairLift**
      gives the pair-lift interaction expressed by
      :math:`c_{i\uparrow}^{\dagger} c_{i\downarrow}^{\phantom{\dagger}} c_{j\uparrow}^{\dagger} c_{j \downarrow}^{\phantom{\dagger}}`.
	    
(2) Green's functions

    | **Initial**
      specifies the one-body Green's function for the initial configuration,
      :math:`\langle c^{\dagger}_{i\sigma_1} c_{j\sigma_2}^{\phantom{\dagger}} \rangle`.
    | **OneBodyG**
      specifies the indices of the one-body Green's function
      :math:`\langle c^{\dagger}_{i\sigma_1} c_{j\sigma_2}^{\phantom{\dagger}} \rangle`
      to be exported.

The details of the file format for each input file is described in the following subsections.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   Trans_file
   InterAll_file
   CoulombIntra_file
   CoulombInter_file
   Hund_file
   PairHop_file
   Exchange_file
   Ising_file
   PairLift_file
   Initial_file
   OneBodyG_file
