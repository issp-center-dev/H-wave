#!/usr/bin/env python3

"""Brute-force, physical-index reference implementations for FLEX.

These functions implement the multi-orbital FLEX building blocks from
Mochizuki-Yanase-Ogata (MYO, cond-mat/0407094) by **explicit nested-loop
direct summation in physical orbital indices** -- no FFT, no flattening of
the rank-4 vertex tensors.  They are deliberately literal and slow; they are
the independent ground truth against which the optimized (FFT + flattened)
general-mode FLEX path is validated to ~1e-10.

Conventions (a SIMPLE test layout, not the production memory layout):

* The Green function ``G`` has shape ``(norb, norb, Nk, nmat)`` on a 1-D
  k-grid of ``Nk`` points, with periodic wrap on both the k index and the
  Matsubara index.

* Irreducible susceptibility / bare bubble:
      chi0_{mn,mu nu}(q, i nu) =
          -(T/N) Sum_{k, i omega} G_{m mu}(k+q, i omega + i nu)
                                   G_{nu n}(k, i omega)

* MYO Eq.(3) (self-energy from a full vertex V):
      Sigma_{mn}(k, i omega) =
          (T/N) Sum_{q, i nu} Sum_{mu nu}
                V_{mu m, nu n}(q, i nu) G_{mu nu}(k-q, i omega - i nu)

The index slots below match these formulas exactly; do not "optimize" the
wiring.

Note on the bubble's index slots (issue #91).  These two formulas are a
*pair*: only one relative orbital-pair order of the bubble reproduces the
exact second-order self-energy when the two are composed as V = U chi0 U.
This file previously carried the bubble with both pairs transposed
(``G_{mu m} G_{n nu}``).  Composed with Eq.(3) that gives, for a density-only
U, ``Sigma_{mn} ~ chi0_{nm} G_{mn}`` instead of ``chi0_{mn} G_{mn}`` -- a
transpose that is invisible on the orbital diagonal and wrong by ~40% off it.
The order below is the one pinned against an exactly-diagonalized 3-orbital
model in ``tests/test_flex_sopt_index_order.py``, which is independent of any
convention choice made inside this codebase.
"""

import numpy as np


def chi0_bruteforce(G, T, Nk):
    """chi0_{mn,mu nu}(q,i nu) by explicit direct summation.

    G: (norb, norb, Nk, nmat). Returns (norb,norb,norb,norb, Nk, nmat).
    chi0[m,n,mu,nu,q,iv] =
        -(T/Nk) Sum_{k,iw} G[m,mu,(k+q)%Nk,(iw+iv)%nmat] G[nu,n,k,iw].
    """
    norb = G.shape[0]; nmat = G.shape[-1]
    chi0 = np.zeros((norb, norb, norb, norb, Nk, nmat), dtype=complex)
    for m in range(norb):
     for n in range(norb):
      for mu in range(norb):
       for nu in range(norb):
        for q in range(Nk):
         for iv in range(nmat):
          s = 0j
          for k in range(Nk):
           for iw in range(nmat):
            s += G[m, mu, (k + q) % Nk, (iw + iv) % nmat] * G[nu, n, k, iw]
          chi0[m, n, mu, nu, q, iv] = -(T / Nk) * s
    return chi0


def sigma_bruteforce(G, V, T, Nk):
    """Sigma_{mn}(k,i omega) by explicit direct summation (MYO Eq.3).

    G: (norb,norb,Nk,nmat). V: (norb,norb,norb,norb,Nk,nmat) indexed
    V[mu,m,nu,n,q,iv] (i.e. V_{mu m, nu n}(q,i nu)). Returns
    (norb,norb,Nk,nmat).
    Sigma[m,n,k,iw] =
        (T/Nk) Sum_{q,iv} Sum_{mu nu}
            V[mu,m,nu,n,q,iv] G[mu,nu,(k-q)%Nk,(iw-iv)%nmat].
    """
    norb = G.shape[0]; nmat = G.shape[-1]
    Sig = np.zeros((norb, norb, Nk, nmat), dtype=complex)
    for m in range(norb):
     for n in range(norb):
      for k in range(Nk):
       for iw in range(nmat):
        s = 0j
        for q in range(Nk):
         for iv in range(nmat):
          for mu in range(norb):
           for nu in range(norb):
            s += V[mu, m, nu, n, q, iv] * G[mu, nu, (k - q) % Nk, (iw - iv) % nmat]
        Sig[m, n, k, iw] = (T / Nk) * s
    return Sig
