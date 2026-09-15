#!/usr/bin/env python3

"""Pins for the single-particle Green function of the sector exact
diagonalisation (``ed_oracle_util.SectorED.green`` / ``density_matrix``).

These are the fast checks that fix the cross-sector Lehmann sum -- its
Boltzmann weighting, its grand-canonical denominators and the global
partition function -- before the chain G4 gate of the off-site second
order is built on them in this module.
"""
import unittest

import numpy as np

from tests import ed_oracle_util as edu


def _fx_chain():
    """``L = 4``, one orbital, ``t = -1``: the brief's single-orbital chain."""
    return edu.EDFixture(L=4, norb=1, t={(0, 0): -1.0}, eps=(0.0,), T=0.5, mu=0.2)


def _fx_orbital():
    """``L = 3``, two orbitals, asymmetric intra/inter-orbital hopping and
    distinct on-site energies: the fixture that can see an orbital index
    order (``tests/test_bond_vs_ed_oracle.py``'s case-M shape, at the same
    ``T`` and ``mu``). ``nmode = 12``, Fock dimension 4096, largest
    ``(N_up, N_dn)`` block 400 -- comfortably inside ``SectorED``."""
    return edu.EDFixture(L=3, norb=2,
                         t={(0, 0): -1.0, (1, 1): -0.7, (0, 1): -0.2, (1, 0): -0.2},
                         eps=(0.0, 0.3), T=0.5, mu=0.2)


class TestGreenFunction(unittest.TestCase):
    """``SectorED.green`` / ``density_matrix``: the single-particle
    observables the chain gate is built on, pinned where they are known in
    closed form or by an independent dense evaluation."""

    def test_free_chain_green_equals_the_one_body_resolvent(self):
        """With no interaction the cross-sector Lehmann sum must collapse to
        ``[i w - (h1 - mu)]^{-1}`` on every mode pair -- which pins the
        Boltzmann weighting, the grand-canonical denominators (``mu``
        enters them, the sectors differ by one particle) and the global
        partition function all at once."""
        fx = _fx_chain()
        iws = 1j * (2 * np.arange(8) + 1 - 8) * np.pi / fx.beta
        G = edu.green_function(fx, iws)
        h1 = fx.build_h1()
        eye = np.eye(fx.nmode)
        for i, iw in enumerate(iws):
            np.testing.assert_allclose(
                G[i], np.linalg.inv(iw * eye - (h1 - fx.mu * eye)), atol=1e-12)

    def test_single_site_interacting_green_equals_the_dense_lehmann_sum(self):
        """One site, one orbital, ``U n_up n_dn``: the sector engine against
        the DENSE Lehmann sum over the full Fock space
        (``EDFixture.annihilators``), which knows nothing about sectors.
        The free pin above cannot see a sector-pairing error that the
        interaction makes matter; this one can."""
        fx = edu.EDFixture(L=1, norb=1, t={(0, 0): 0.0}, eps=(0.3,), T=0.5, mu=0.2)
        terms = edu.canonical_density_terms(fx, [(0, 0, 0, 0.7)])
        iws = 1j * (2 * np.arange(6) + 1 - 6) * np.pi / fx.beta
        G = edu.green_function(fx, iws, terms)

        C = fx.annihilators()
        CD = [c.conj().T for c in C]
        h1 = fx.build_h1()
        H = edu.h_int_from_terms(fx, terms)
        for p in range(fx.nmode):
            for q in range(fx.nmode):
                if h1[p, q] != 0:
                    H = H + h1[p, q] * (CD[p] @ C[q])
        N = sum(CD[p] @ C[p] for p in range(fx.nmode))
        E, V = np.linalg.eigh(H - fx.mu * N)
        E = E - E.min()
        w = np.exp(-fx.beta * E)
        cm = [V.conj().T @ C[p] @ V for p in range(fx.nmode)]
        boltz = w[:, None] + w[None, :]
        dE = E[None, :] - E[:, None]
        ref = np.zeros_like(G)
        for p in range(fx.nmode):
            for q in range(fx.nmode):
                num = cm[p] * np.conj(cm[q]) * boltz
                for i, iw in enumerate(iws):
                    ref[i, p, q] = (num / (iw - dE)).sum()
        ref /= w.sum()
        self.assertGreater(np.abs(ref).max(), 1e-3)            # anti-vacuity
        np.testing.assert_allclose(G, ref, atol=1e-12)

    def test_free_density_matrix_equals_the_fermi_occupation(self):
        """``SectorED.density_matrix`` at zero coupling is the free
        ``<c^dag_p c_q> = [W f(ev) W^dagger]^T`` of ``h1 - mu``, the same
        expression ``hf_h1_from_terms`` Wick-contracts against."""
        fx = _fx_orbital()
        h1 = fx.build_h1()
        ev, W = np.linalg.eigh(h1 - fx.mu * np.eye(fx.nmode))
        f = 1.0 / (np.exp(np.clip(fx.beta * ev, -500, 500)) + 1.0)
        np.testing.assert_allclose(edu.SectorED(fx).density_matrix(),
                                   ((W * f) @ W.conj().T).T, atol=1e-12)



if __name__ == "__main__":
    unittest.main()
