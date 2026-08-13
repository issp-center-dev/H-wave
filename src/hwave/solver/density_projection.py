"""Density-diagonal projections of the rank-4 interaction tensor.

The ``reduced`` calculation scheme solves with the density-density
diagonal of the interaction. Where the particle-hole vertex is itself
density diagonal this is exact, not an approximation: the RPA series
closes on precisely those slots, so the projection discards nothing.
Interactions with no density-diagonal content -- Exchange and PairHop
-- are *refused* by the callers, with a pointer to
``calc_scheme='general'``, rather than being silently dropped, and the
auto selection picks ``general`` when they are present; PairLift is
accepted with a warning, its particle-hole vertex being zero in every
scheme. The two solvers used to carry their own copies of these einsum
reductions; per issue #107 the projection lives here once.

Layout convention:

* ``project_density_pairs``: combined-index form. The tensor is read as
  ``(nvol, nd, nd, nd, nd)`` with ``nd`` the (spin-)orbital dimension of
  the run, and the pair diagonal ``[a, a, b, b]`` is extracted --
  ``'kaabb->kab'``. Used by the plain reduced/spinful branches of both
  solvers; FLEX's spin-resolved reduced branch is the same extraction
  through a ``(ns, norb)`` factorized view (``'ksasatbtb->ksatb'``), so
  it routes here after reshaping.

``xp`` is the array backend (NumPy or CuPy); einsum dispatches on it.
"""


def project_density_pairs(ham_q, nvol, nd, xp):
    """Pair-diagonal 'kaabb->kab' of a combined-index rank-4 tensor."""
    return xp.einsum(
        'kaabb->kab',
        ham_q.reshape(nvol, *(nd,) * 4)).reshape(nvol, nd, nd)
