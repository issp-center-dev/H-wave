.. highlight:: none

occupation.npz
~~~~~~~~~~~~~~

The converged SCF occupation, chemical potential, temperature, and
per-column metadata are saved to ``occupation.npz`` (key name
configurable via ``[file.output].occupation``; default
``occupation.npz``). Consumers (notably the ``tools/uhfk_to_mvmc.py``
bridge) need this information together with ``eigen.npz`` to project
finite-temperature occupations to T=0 Slater determinants without
re-solving the SCF.

Keys
^^^^

- ``occupation`` (shape ``(nvol, nd)``, ``float64``)
    Fermi-Dirac weight :math:`f(\epsilon_{k,n})` per (k-point, band)
    used at the final SCF iteration. At T=0 the values are 0 or 1; at
    T>0 they are fractional.

- ``mu`` (shape ``(n_mu_groups,)``, ``float64``)
    Chemical potential per mu-group. With ``2Sz = 0`` in Sz-fixed mode,
    there are 2 groups (up, down); in Sz-free mode there is 1 group.

- ``T`` (scalar ``float64``)
    Temperature used during SCF.

- ``column_spin`` (shape ``(nd,)``, ``int64``)
    Spin character of each eigenvector column. ``0`` = up-only block,
    ``1`` = down-only block, ``-1`` = mixed (Sz-free).

- ``column_mu_group`` (shape ``(nd,)``, ``int64``)
    Mu-group index for each column. ``mu[column_mu_group[n]]`` gives the
    chemical potential associated with column ``n``.

Conventions
^^^^^^^^^^^

The column layout matches the ``eigenvector`` array in ``eigen.npz``:
columns are arranged by block order from ``UHFk._init_block_structure``,
and rows live in the original ``nd`` index space.
