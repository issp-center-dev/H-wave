"""An independent Hartree-Fock oracle for the half-filled square-lattice
Hubbard model.

The Neel state of the L x L square lattice (nearest-neighbour hopping t,
on-site U, half filling, T = 0) is solved from the gap equation directly,
without any part of H-wave:

* ``eps_k = -2 t (cos kx + cos ky)`` with ``kx, ky = 2 pi j / L``;
* the staggered moment ``m`` solves
  ``m = (1/N) sum_k U m / (2 sqrt(eps_k^2 + U^2 m^2))``;
* the total energy is ``E = -sum_k sqrt(eps_k^2 + U^2 m^2) + N U (m^2 + 1/4)``.

The golden values, asserted in tests/test_results_run.py:

* ``neel_gap_equation(4, 8)``  -> m = 0.4474190440247483, E = -7.38962582304832
* ``neel_gap_equation(16, 2)`` -> m = 0.1917429137967485, E = -291.4365477135228
"""
import numpy as np
from scipy.optimize import brentq


def neel_gap_equation(L: int, U: float, t: float = 1.0):
    """Return ``(m, E)`` of the Hartree-Fock Neel state (see the module text).

    ``m`` is the staggered magnetisation per site (at most 1/2) and ``E`` the
    total energy of the N = L*L site lattice."""
    k = 2.0 * np.pi * np.arange(L) / L
    kx, ky = np.meshgrid(k, k, indexing="ij")
    eps = (-2.0 * t * (np.cos(kx) + np.cos(ky))).ravel()
    n = eps.size

    def gap(m):
        # the gap equation divided by m (m = 0 is always a solution)
        return np.sum(U / (2.0 * np.sqrt(eps ** 2 + (U * m) ** 2))) / n - 1.0

    m = brentq(gap, 1e-6, 0.5, xtol=1e-15, rtol=4 * np.finfo(float).eps)
    energy = -np.sum(np.sqrt(eps ** 2 + (U * m) ** 2)) + n * U * (m ** 2 + 0.25)
    return float(m), float(energy)
