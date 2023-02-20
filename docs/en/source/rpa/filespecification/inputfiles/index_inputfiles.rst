.. highlight:: none

.. _Ch:HowToWannier90_rpa:

Input files for RPA
--------------------------------

In this section, the input files for the random phase approximation (RPA) are described.
They are classified into two categories, and written in Wannier90 format.

(1) Geometry information
    
    **Geometry**
      defines the geometrical information of the lattice.

(2) Interaction definitions

    These files defines the Hamiltonian for UHF in the form of electron systems.
    They provide the coefficients of the interaction terms associated with the
    specified keywords.

    The following keywords adopted in :math:`{\mathcal H}\Phi` and mVMC in the Expert Mode
    are accepted.
    
    **Transfer**
      corresponds to one-body term denoted by
      :math:`c_{i\sigma_1}^{\dagger} c_{j\sigma_2}^{\phantom{\dagger}}`.

    **CoulombIntra**
      corresponds to the interaction denoted by 
      :math:`n_{i\uparrow} n_{i\downarrow}`, where :math:`n_{i\sigma} = c_{i\sigma}^{\dagger} c_{i\sigma}^{\phantom{\dagger}}`. 

    **CoulombInter**
      corresponds to the interaction denoted by
      :math:`n_{i} n_{j}`, where :math:`n_i = n_{i\uparrow} + n_{i\downarrow}`.

    **Hund**
      corresponds to the interaction denoted by
      :math:`n_{i\uparrow} n_{j\uparrow} + n_{i\downarrow} n_{j\downarrow}`.
    
    **Ising**
      corresponds to the interaction denoted by
      :math:`S_i^z S_j^z`.

    **Exchange**
      corresponds to the interaction denoted by
      :math:`S_i^+ S_j^-`.

    **PairLift**
      corresponds to the interaction denoted by
      :math:`c_{i\uparrow}^{\dagger} c_{i\downarrow}^{\phantom{\dagger}} c_{j\uparrow}^{\dagger} c_{j\downarrow}^{\phantom{\dagger}}`.

    **PairHop**
      corresponds to the interaction denoted by
      :math:`c_{i\uparrow}^{\dagger} c_{j\uparrow}^{\phantom{\dagger}} c_{i\downarrow}^{\dagger} c_{j\downarrow}^{\phantom{\dagger}}`.

The data formats are described in the following sectoins.

.. toctree::
  :maxdepth: 1

  geometry
  interaction

