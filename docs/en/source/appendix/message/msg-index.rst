List of error messages
======================

- ``mode is not defined in [mode].``

  **description :** ``mode`` parameter is missing in ``[mode]`` section of the input parameter file.

  **mode :** main

- ``Get_param: key must be mod or ham or output.``

  **description :** unsupported keyword is given to ``get_param()``.

  **mode :** UHFr (read_input)

- ``duplicate items found in`` *file*

  **description :** the file *file* contains duplicate entries.

  **mode :** UHFr (read_input)

- ``incorrect number of lines in`` *file* ``: expected=`` *N* ``, found=`` *M*

  **description :** number of lines of the input file does not match the description in the file.

  **mode :** UHFr (read_input)

- ``Unknown keyword`` *keyword*

  **description :** unsupported keyword is found in ``[file.input.interaction]``.

  **mode :** UHFk (read_input_k)

- ``initial and initial_uhf can not be specified simultaneously.``

  **description :** ``initial`` and ``initial_uhf`` cannot be specified simultaneously.

  **mode :** UHFk (read_input_k)

- ``read_input_k: file`` *file* ``not found``

  **description :** the file *file* cannot be found.

  **mode :** UHFk (read_input_k)

- ``Get_param: key must be mod or ham or output.``

  **description :** unsupported keyword is given to ``get_param()``.

  **mode :** UHFk (read_input_k)

- ``read_geom: file`` *file* ``not found``

  **description :** the file *file* specified by ``Geometry`` keyword cannot be found.

  **mode :** UHFk (wan90)

- ``mode.param.2Sz must be even(odd) when Ncond is even(odd).``

  **description :** even/odd mismatch between ``2Sz`` and ``Ncond``.

  **mode :** solver base

- ``range check for`` *type* ``failed.``

  **description :** the value of *type* is not appropriate.

  **mode :** UHFk

- ``_check_cellsize failed. interaction range exceeds cell shape.``

  **description :** some of translation vectors of the interaction description do not lie within the CellShape.

  **mode :** UHFk

- ``Hermiticity check failed: |T_ba(-r)^* - T_ab(r)| =`` *val*

  **description :** Transfer term is not Hermite.

  **mode :** UHFk

- ``Parameter range check failed for param_mod.``

  **description :** the parameter value in [mode.param] is out of range.

  **mode :** solver base

- ``Parameter check failed for param_mod.``

  **description :** the parameter value in [mode.param] is inappropriate.

  **mode :** solver base

- ``Hermite check failed for`` *type*

  **description :** *type* is not Hermite.

  **mode :** UHFr

- ``Parameter check failed for info_mode.``

  **description :** the parameter value in [mode] is inappropriate.

  **mode :** solver base

- ``value not integer``

  **description :** the parameter value is not an integer.

  **mode :** RPA

- ``Lattice initialization failed: 'CellShape' not found.``

  **description :** ``CellShape`` is missing in [mode.param].

  **mode :** RPA

- ``Ncond must be greater than zero: Ncond=`` *Ncond*

  **description :** the value of ``Ncond`` is not appropriate.

  **mode :** RPA

- ``Nmat must be greater than zero: Nmat=`` *Nmat*

  **description :** the value of ``Nmat`` is not appropriate.

  **mode :** RPA

- ``RPA._find_mu: not converged. abort``

  **description :** the calculation of ``mu`` does not converge.

  **mode :** RPA

- ``SubShape is not compatible with CellShape.``

  **description :** the value of SubShape does not divide that of CellShape.

  **mode :** RPA

- ``T must be greater than or equal to zero: T=`` *T*

  **description :** the value of ``T`` is not appropriate.

  **mode :** RPA

- ``both mu and Ncond or filling are specified``

  **description :** ``mu` and ``Ncond`` or ``filling`` should not specified simultaneously.

  **mode :** RPA

- ``dimension of CellShape must be one, two, or three.``

  **description :** the dimension of ``CellShape`` is not appropriate.

  **mode :** RPA

- ``dimension of SubShape does not match with that of CellShape.``

  **description :** the dimension of ``SubShape`` is not appropriate.

  **mode :** RPA

- ``invalid CellShape.``

  **description :** the value of ``CellShape`` is not appropriate.

  **mode :** RPA

- ``invalid SubShape.``

  **description :** the value of ``SubShape`` is not appropriate.

  **mode :** RPA

- ``none of mu, Ncond, nor filling is specified``

  **description :** one of ``mu``, ``Ncond``, or ``filling`` should be specified.

  **mode :** RPA

- ``read_chi0q failed:`` *info*

  **description :** reading ``chi0q`` from file was not successful.

  **mode :** RPA

- ``round_to_int: unknown mode`` *mode*

  **description :** unsupported rounding mode is specified.

  **mode :** RPA

- ``unexpected data size`` *error*

  **description :** data size is not as expected.

  **mode :** RPA

- ``mode is not defined in [mode].``

  **description :** the ``mode`` parameter is missing in ``[mode]``.

  **mode :** RPA

- ``orbital index check failed for`` *type*

  **description :** the indices of the orbitals are inappropriate.

  **mode :** UHFk

- ``initial green function in coord space requires geometry.dat``

  **description :** ``geometry.dat`` must also be specified when the coordinate space Green's function.

  **mode :** UHFk

- ``CellShape is missing. abort``

  **description :** ``CellShape`` parameter is missing.

  **mode :** UHFk

- ``Ncond or Nelec is missing. abort``

  **description :** ``Ncond`` or ``Nelec`` parameter is missing.

  **mode :** UHFk

- ``SubShape is not compatible with CellShape. abort``

  **description :** the value of ``SubShape`` does not divide that of ``CellShape``.

  **mode :** UHFk

- ``_check_orbital_index failed. invalid orbital index found in interaction definitions.``

  **description :** the indices of the orbitals in interaction definition files are inappropriate.

  **mode :** UHFk

- ``_save_greenone: onebodyg_uhf and geometry_uhf are required``

  **description :** ``onebodyg_uhf`` and ``geometry_uhf`` are not provided.

  **mode :** UHFk

- ``find mu: not converged. abort``

  **description :** the calculation of ``mu`` does not converge.

  **mode :** UHFk

- ``range check failed for Initial``

  **description :** the values of ``Initial`` are inappropriate.

  **mode :** UHFr

- ``OneBodyG is required to output green function.``

  **description :** ``OneBodyG`` is missing for the output of Green's function.

  **mode :** UHFr

- ``hermite check failed for Initial``

  **description :** ``Initial`` is not Hermite.

  **mode :** UHFr

- ``Range check failed for Transfer``

  **description :** the indices of ``Transfer`` definition file are out of range.

  **mode :** UHFr

- ``Range check failed for`` *type*

  **description :** the indices of *type* definition file are out of range.

  **mode :** UHFr

- ``parameter range check failed.``

  **description :** the value of the parameter is not appropriate.

  **mode :** UHFr

- ``mode is incorrect: mode=`` *mode*

  **description :** ``mode`` parameter is not appropriate.

  **mode :** UHFr

- ``mode.param.`` *key* ``must be greater than`` *value*

  **description :** the value of parameter *key* in [mode.param] is inappropriate.

  **mode :** solver base [warning]

- ``"mode.param.`` *key* ``must be smaller than`` *value*

  **description :** the value of parameter *key* is [mode.param] is inappropriate.

  **mode :** solver base [warning]

- ``mode.param.`` *key* ``is not defined.``

  **description :** parameter *key* is not found in [mode.param].

  **mode :** solver base [warning]

- ``mode.`` *key* ``in mode section is incorrect:`` *values*

  **description :** ``mode`` parameter in ``[mode]`` section is not valid.

  **mode :** solver base [warning]

- ``mode.`` *key* ``is not defined.``

  **description :** ``mode`` parameter is not found in ``[mode]`` section.

  **mode :** solver base [warning]

- ``TRUST-ME mode enabled. parameter checks are relaxed``

  **description :** ``TRUST-ME`` mode is enabled. the parameter checks will be omitted.

  **mode :** solver base [warning]

- ``value not integer``

  **description :** the specified value is not an integer.

  **mode :** RPA [warning]

- ``mode is incorrect: mode=`` *mode*

  **description :** ``mode`` parameter is not valid.

  **mode :** RPA [warning]

- ``FATAL: 2Sz=`` *value* ``. 2Sz should be even for calculating fij``

  **description :** ``2Sz`` must be an even number for the calculation of :math:`f_{ij}`.

  **mode :** UHFr [warning]

- ``FATAL: Ne=`` *value* ``. Ne should be even for calculating fij``

  **description :** ``Ne`` must be an even number for the calculation of :math:`f_{ij}`.

  **mode :** UHFr [warning]

- ``NOT IMPLEMENTED: Sz even and Sz != 0: this case will be implemented in near future``

  **description :** the calculation of :math:`f_{ij}` is not yet supported when ``Sz`` is an even number except zero.

  **mode :** UHFr [warning]

- ``key`` *key* ``is wrong!``

  **description :** the keyword *key* is invalid.

  **mode :** UHFr [warning]

- ``UHFr calculation is failed: rest=`` *residue* ``, eps=`` *eps*

  **description :** the calculation of UHFr does not converge.

  **mode :** UHFr [warning]

