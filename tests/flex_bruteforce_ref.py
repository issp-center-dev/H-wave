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
*pair*: only one relative orbital-pair order of the bubble reproduces the exact
second-order self-energy when the two are composed as V = U chi0 U.

The bubble is deliberately NOT the verbatim MYO Eq.(5), which reads
``chi0_{mn,mu nu} = -(T/N) sum_k G_{mu m}(k+q) G_{n nu}(k)`` -- both orbital
pairs the other way round from the form below.  This file used to carry that
verbatim form.  Composed with Eq.(3) it gives, for a density-only U,
``Sigma_{mn} ~ chi0_{nm} G_{mn}`` instead of ``chi0_{mn} G_{mn}``: a difference
invisible on the orbital diagonal and ~40% off it, measured against an
exactly-diagonalized 3-orbital model.

That is a mapping problem, not an error in the paper: the paper's Green
function index order is not the one this codebase uses, so transcribing its
equations term by term does not carry the convention with them.  What is pinned
here is the COMPOSITION -- the pair (bubble, Eq.(3)) must reproduce the exact
O(U^2) self-energy -- and that is what ``tests/test_flex_sopt_index_order.py``
checks, independently of any convention choice made inside this codebase.
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
