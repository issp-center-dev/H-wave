"""Linearized Eliashberg equation solver for superconducting instability analysis.

``hwave_sc`` is a post-processing tool that uses the bare susceptibility
chi0q from H-wave's RPA solver (or the dressed susceptibilities of a FLEX
run) to solve the linearized Eliashberg equation and analyze superconducting
instabilities. The implementation lives in :mod:`hwave.solver.eliashberg`,
one module per layer (dev/sc-split-plan.md): ``config``, ``band``,
``vertex``, ``kernel``, ``eigen``, ``io``, ``output``, ``bond``,
``so_solve``, the static entry ``static`` and the dynamic solvers
``dynamic``, ``dynamic_so``, ``bond_dynamic``, ``bond_io``.

This module is the entry shell: the console script ``hwave_sc`` is
``hwave.sc:main``, and :func:`calc_eliashberg`, :func:`main`,
:func:`reject_spin_orbital_mode` and the public constants are re-exported
here for the temperature sweep and for external scripts. Everything else is
imported from its module.
"""

from hwave.solver.eliashberg.static import (  # noqa: F401
    DEFAULT_SOLVER_MODE,
    EXIT_EIGEN_NOT_CONVERGED,
    IN_PROCESS_DEFAULT_SOLVER_MODE,
    calc_eliashberg,
    main,
    reject_spin_orbital_mode,
)

__all__ = [
    "DEFAULT_SOLVER_MODE",
    "EXIT_EIGEN_NOT_CONVERGED",
    "IN_PROCESS_DEFAULT_SOLVER_MODE",
    "calc_eliashberg",
    "main",
    "reject_spin_orbital_mode",
]


if __name__ == "__main__":
    main()
