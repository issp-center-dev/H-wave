#!/usr/bin/env python3

"""Bit-identity tests for the block-diagonal scatter rewrite of the
Kronecker-with-identity einsums in flex.py.

The optimization replaces

    np.einsum('lkab,st->lksatb', chi0q, I_ns).reshape(nfreq, nvol, nd, nd)   # spin-free
    np.einsum('glkab,gh->lkgahb', chi0q, I_ns).reshape(nfreq, nvol, nd, nd)  # spin-diag
    np.einsum('kab,st->ksatb', u, I_ns).reshape(nvol, nd, nd)                # vertices

with a preallocated-zeros + per-spin-block scatter.  This is a pure data
rearrangement (Kronecker product with the identity), so the result MUST be
bit-identical to the einsum.  These tests paste the original einsums as the
reference and assert exact equality.
"""

import unittest
import numpy as np


def _scatter_spin_free(chi0q, ns):
    """New implementation for the spin-free branch.

    chi0q: (nfreq, nvol, norb, norb) -> (nfreq, nvol, nd, nd)
    """
    nfreq, nvol, norb, _ = chi0q.shape
    nd = ns * norb
    out = np.zeros((nfreq, nvol, nd, nd), dtype=chi0q.dtype)
    for s in range(ns):
        sl = slice(s * norb, (s + 1) * norb)
        out[..., sl, sl] = chi0q
    return out


def _scatter_spin_diag(chi0q, ns):
    """New implementation for the spin-diag branch.

    chi0q: (nblock, nfreq, nvol, norb, norb) -> (nfreq, nvol, nd, nd)
    where nblock == ns; source block g goes to spin-block (g, g).
    """
    nblock, nfreq, nvol, norb, _ = chi0q.shape
    nd = ns * norb
    out = np.zeros((nfreq, nvol, nd, nd), dtype=chi0q.dtype)
    for s in range(ns):
        sl = slice(s * norb, (s + 1) * norb)
        out[..., sl, sl] = chi0q[s]
    return out


def _scatter_vertex(u, ns):
    """New implementation for _build_spin_charge_vertices.

    u: (nvol, norb, norb) -> (nvol, nd, nd)
    """
    nvol, norb, _ = u.shape
    nd = ns * norb
    out = np.zeros((nvol, nd, nd), dtype=u.dtype)
    for s in range(ns):
        sl = slice(s * norb, (s + 1) * norb)
        out[..., sl, sl] = u
    return out


class TestInflateScatterEquivalence(unittest.TestCase):
    """Bit-identical equivalence of scatter vs Kronecker-identity einsum."""

    def _rand_complex(self, shape):
        rng = np.random.default_rng(12345)
        return rng.standard_normal(shape) + 1j * rng.standard_normal(shape)

    def test_spin_free_ns2(self):
        nfreq, nvol, norb, ns = 5, 7, 3, 2
        nd = ns * norb
        chi0q = self._rand_complex((nfreq, nvol, norb, norb))
        I_ns = np.identity(ns)
        ref = np.einsum('lkab,st->lksatb', chi0q, I_ns).reshape(
            nfreq, nvol, nd, nd)
        new = _scatter_spin_free(chi0q, ns)
        np.testing.assert_array_equal(new, ref)

    def test_spin_free_ns1(self):
        nfreq, nvol, norb, ns = 4, 6, 2, 1
        nd = ns * norb
        chi0q = self._rand_complex((nfreq, nvol, norb, norb))
        I_ns = np.identity(ns)
        ref = np.einsum('lkab,st->lksatb', chi0q, I_ns).reshape(
            nfreq, nvol, nd, nd)
        new = _scatter_spin_free(chi0q, ns)
        np.testing.assert_array_equal(new, ref)

    def test_spin_diag_ns2(self):
        nblock, nfreq, nvol, norb, ns = 2, 5, 7, 3, 2
        nd = ns * norb
        chi0q = self._rand_complex((nblock, nfreq, nvol, norb, norb))
        I_ns = np.identity(ns)
        ref = np.einsum('glkab,gh->lkgahb', chi0q, I_ns).reshape(
            nfreq, nvol, nd, nd)
        new = _scatter_spin_diag(chi0q, ns)
        np.testing.assert_array_equal(new, ref)

    def test_vertex_ns2(self):
        nvol, norb, ns = 7, 3, 2
        nd = ns * norb
        u = self._rand_complex((nvol, norb, norb))
        I_ns = np.eye(ns)
        ref = np.einsum('kab,st->ksatb', u, I_ns).reshape(nvol, nd, nd)
        new = _scatter_vertex(u, ns)
        np.testing.assert_array_equal(new, ref)

    def test_vertex_ns1(self):
        nvol, norb, ns = 6, 2, 1
        nd = ns * norb
        u = self._rand_complex((nvol, norb, norb))
        I_ns = np.eye(ns)
        ref = np.einsum('kab,st->ksatb', u, I_ns).reshape(nvol, nd, nd)
        new = _scatter_vertex(u, ns)
        np.testing.assert_array_equal(new, ref)


if __name__ == '__main__':
    unittest.main()
