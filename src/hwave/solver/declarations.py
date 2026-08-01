"""Declaration symmetrisation, single-sourced (#108 increment 2).

An interaction file may declare the same physical bond from both ends:
the two declarations of one coupling are ``(R, a, b)`` and
``(-R, b, a)``. For every type except PairHop they multiply the SAME
operator, so the table entering any vertex is the mean

    T~[r, a, b] = (T[r, a, b] + T[-r, b, a]) / 2.

PairHop's partner entry is the HERMITIAN one, so its mean conjugates and
a Hermitian-closed complex declaration keeps its physical phase
(adjudicated in the #100/#105/#106 series; an asymmetric same-operator
declaration means the mean coefficient, and an antisymmetric one is an
identically zero Hamiltonian).

This module owns that content in BOTH representations it is consumed in:

* :func:`symmetrise_dense` -- real space, on a dense
  ``(nx, ny, nz, norb, norb)`` array. Used by the RPA ring's interaction
  assembly (``rpa._make_ham_inter``), which works on declaration tables
  before any Fourier transform. Dense on purpose: table keys may sit in
  a wrapped canonical form (``(n-1, 0, 0)`` for a ``-x`` bond; folded
  tables store canonicalized displacements), where a ``(-R)`` key lookup
  would silently miss the partner and halve the coefficient.
* :func:`symmetrise_k` -- momentum space, on the FFT'd
  ``(norb, norb, nx, ny, nz)`` tables. Used by the pair-space S/C
  builders in :mod:`hwave.sc` (RPA/FLEX/Eliashberg vertex side), whose
  pipeline transforms first and symmetrises second.

The two forms are the SAME reduction: the Fourier transform is linear,
so ``fft(symmetrise_dense(T)) == symmetrise_k(fft(T))`` exactly in
exact arithmetic (the r-space partner at ``(-R, b, a)`` maps to the
orbital transpose AT ``-q``; PairHop's conjugation at fixed ``q`` is the
Fourier image of {conjugate, ``R -> -R``}). That identity is pinned
executably in tests/test_declarations.py; the numerical results differ
only by float summation order, which is why BOTH forms exist here
instead of one delegating to the other through transforms -- each
consumer keeps its historical bit-exact pipeline while the content
lives in one module.

Relation to UHFk: ``uhfk.py`` stores the HERMITIAN mean for every type
(its table is contracted over both ordered orbital pairs, so a complex
same-operator coefficient p appears as p + conj(p) there). The two
conventions give the same physical Hamiltonian; see
``sc._symmetrise_interactions_k``'s original derivation, now in
:func:`symmetrise_k`'s docstring.
"""

import numpy as np

from hwave.solver.kgrid import reverse_fft_axes


def symmetrise_dense(arr, hermitian=False):
    """Mean of a dense real-space table with its reversed-bond partner.

    Parameters
    ----------
    arr : ndarray
        Dense table, shape ``(nx, ny, nz, norb, norb)``, index 0 of each
        lattice axis = the on-site/origin cell (FFT order).
    hermitian : bool
        False (default): the same-operator mean
        ``0.5 * (T[r, a, b] + T[-r, b, a])`` -- every type but PairHop.
        True: the Hermitian-partner mean, which conjugates the reversed
        term -- PairHop.

    Returns
    -------
    ndarray
        The symmetrised table, same shape.
    """
    rev = np.transpose(reverse_fft_axes(arr, (0, 1, 2)),
                       (0, 1, 2, 4, 3))
    if hermitian:
        rev = np.conjugate(rev)
    return 0.5 * (arr + rev)


def symmetrise_k(inter_k):
    """Momentum-space form of the same reduction, per interaction type.

    The same-operator partner of ``M[b, a](q)`` is the orbital-transposed
    entry AT ``-q`` (arrays carry ``e^{-iqR}`` phases) -- averaging with
    the same-q transpose instead corrupted every off-site interaction,
    collapsing a one-direction bond's ``V(q) = v e^{-iqR}`` to
    ``v cos(qR)`` (measured: the S/C entry vanished outright at
    ``q = pi/2``). The mean used here:

    * every type except PairHop: ``0.5 * (M(q) + M(-q)^T)``. For a
      both-ends declaration this mean is an identity; for an on-site
      complex Hermitian-closed declaration it drops the inert imaginary
      part; for an antisymmetric declaration -- an identically zero
      Hamiltonian -- it gives zero.
    * PairHop: ``0.5 * (M(q) + conj(M(q)^T))`` at the SAME q.
      Conjugation at fixed q is the Fourier image of
      {conjugate, ``R -> -R``} (the #105 derivation), so this is an
      identity for Hermitian-closed input on-site and off-site, and
      preserves the physical complex phase.

    Idempotent, so it is safe that both the all-q and per-q S/C builders
    apply it.

    Parameters
    ----------
    inter_k : dict
        ``{type: ndarray}`` with arrays of shape
        ``(norb, norb, Nx, Ny, Nz)``.

    Returns
    -------
    dict
        Symmetrised tables, same shapes.
    """
    out = {}
    for itype, M in inter_k.items():
        if itype == "PairHop":
            out[itype] = 0.5 * (M + np.conj(M.transpose(1, 0, 2, 3, 4)))
            continue
        Mrev = reverse_fft_axes(M, (2, 3, 4))
        out[itype] = 0.5 * (M + Mrev.transpose(1, 0, 2, 3, 4))
    return out
