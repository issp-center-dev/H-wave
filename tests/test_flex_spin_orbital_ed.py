#!/usr/bin/env python3

"""Exact diagonalization of the spin-orbital FLEX second order with
spin-orbit coupling (issue #83).

The model is the single-site three-orbital model of
``tests/test_flex_second_order_ed_onsite.py`` (G3) with a random complex
spin-mixing one-body block ``F`` added (``F_ab c^dag_{a up} c_{b dn} +
h.c.``), so the bare propagator has spin-off-diagonal blocks and no spin
component is conserved. The method is G3's: the O(x^2) coefficient of the
EXACT fluctuation remainder

    Sigma_ED - Sigma_HF[rho_ED]

(``Sigma_HF`` the test-side Wick functional ``wick_hf_sigma`` of the ED
monomials, evaluated on the exact spin-orbital density, off-diagonal spin
blocks included) must equal the coefficient of the production
``Sigma_fluct`` of ONE spin-orbital FLEX map (``enable_spin_orbital =
true``, ``flex_hartree_fock = true``) on the bare Green function at the
fixed ``mu``: bubble ``_calc_chi0q``, effective interaction
``_flex_compute_veff_spinful``, kernel ``_calc_self_energy_general``.

The coefficient is extracted with G3's stencil and refinement
(``_coeff2``, ``_refine``) over the innermost eighth of the Matsubara grid.
Measured worst relative deviation over the seven types: about 3e-3 at
``Nmat = 1024, x = 0.0125`` and 6e-4 to 8e-4 at ``Nmat = 4096, x =
0.00625``; the heavy gate runs at the latter working point with G3's 2e-3
tolerance. A fast smoke variant at ``Nmat = 512`` (deviation about 6e-3)
keeps three types in the default gate.

Both gates also assert that the compared coefficient is not vacuous and
that its spin-off-diagonal block is a sizeable part of it, so the
spin-mixing content is really being compared.
"""
import logging
import os
import shutil
import tempfile
import unittest

import numpy as np

from tests.heavy_tests import heavy
from tests import test_flex_sopt_index_order as ref
from tests import test_flex_second_order_ed_onsite as G3
from tests.test_second_order_factors import _write_wan

C, CD, N, BETA, MU, H0 = ref.C, ref.CD, ref.NORB, ref.BETA, ref.MU, ref.H0
M = 2 * N

#: amplitude of the spin-mixing block
_SOC = 0.25


def _soc_block():
    rng = np.random.default_rng(11)
    return _SOC * (rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N)))


F = _soc_block()
#: the one-body matrix on the spin-block index ``s * N + a``
HG = np.zeros((M, M), complex)
HG[:N, :N] = H0
HG[N:, N:] = H0
HG[:N, N:] = F
HG[N:, :N] = F.conj().T


def _hamiltonian(terms):
    H = G3._hamiltonian(terms)
    with ref._quiet_matmul():
        for a in range(N):
            for b in range(N):
                if F[a, b] != 0:
                    op = CD[ref._so(a, 0)] @ C[ref._so(b, 1)]
                    H = H + F[a, b] * op + np.conj(F[a, b]) * op.conj().T
    return H


def _eig(H):
    E, V = np.linalg.eigh(H)
    E = E - E.min()
    w = np.exp(-BETA * E)
    return E, V, w, w.sum()


def _gen_op(p):
    """The annihilator of spin-block index ``p``."""
    return C[ref._so(p % N, p // N)]


def _green(H, iws):
    E, V, w, Z = _eig(H)
    with ref._quiet_matmul():
        cm = [V.conj().T @ _gen_op(p) @ V for p in range(M)]
    boltz = w[:, None] + w[None, :]
    dE = E[None, :] - E[:, None]
    G = np.zeros((len(iws), M, M), complex)
    for p in range(M):
        for q in range(M):
            num = cm[p] * np.conj(cm[q]) * boltz
            G[:, p, q] = (num[None] / (iws[:, None, None] - dE[None])).sum(axis=(1, 2))
    return G / Z


def _rho(H):
    E, V, w, Z = _eig(H)
    r = np.zeros((M, M), complex)
    with ref._quiet_matmul():
        for p in range(M):
            for q in range(M):
                op = V.conj().T @ (_gen_op(p).conj().T @ _gen_op(q)) @ V
                r[p, q] = (np.diag(op) * w).sum() / Z
    return r


def _fluct_ed(terms, iws):
    """``Sigma_ED - Sigma_HF[rho_ED]`` on the spin-block index."""
    H = _hamiltonian(terms)
    G = _green(H, iws)
    G0 = np.array([np.linalg.inv((iw + MU) * np.eye(M) - HG) for iw in iws])
    sig = np.linalg.inv(G0) - np.linalg.inv(G)
    hf = G3.wick_hf_sigma(M, G3._mode_terms(terms), G3._gen_to_mode(_rho(H)))
    return sig - G3._mode_to_gen(hf)[None]


def _fluct_prod(terms, nmat):
    """Production ``Sigma_fluct`` of one spin-orbital FLEX map on the bare
    Green function at ``mu = MU``."""
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as flex_mod
    d = tempfile.mkdtemp(prefix="flex_so_ed_")
    try:
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("1 0 0\n0 1 0\n0 0 1\n{}\n".format(M) + "0 0 0\n" * M)
        rows = []
        for p in range(M):
            for q in range(M):
                if HG[p, q] != 0:
                    # the transfer file uses the interleaved index 2 a + s
                    ip = 2 * (p % N) + p // N + 1
                    iq = 2 * (q % N) + q // N + 1
                    rows.append((0, 0, 0, ip, iq, HG[p, q].real, HG[p, q].imag))
        _write_wan(os.path.join(d, "transfer.dat"), "Transfer", M, rows)
        idict = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
        for itype, v in terms:
            _write_wan(os.path.join(d, itype + ".dat"), itype, N, G3._rows(itype, v))
            idict[itype] = itype + ".dat"
        r = read_input_k.QLMSkInput({"path_to_input": d, "interaction": idict})
        par = {"T": 1.0 / BETA, "mu": MU, "CellShape": [1, 1, 1], "SubShape": [1, 1, 1],
               "Nmat": nmat, "IterationMax": 1, "Mix": 1.0, "EPS": 1,
               "flex_hartree_fock": True}
        prev = logging.root.manager.disable
        logging.disable(logging.WARNING)
        try:
            s = flex_mod.FLEX(r.get_param("ham"), {}, {
                "mode": "FLEX", "param": par, "enable_spin_orbital": True,
                "calc_scheme": "general"})
            s._calc_epsilon_k({})
            G = s._calc_dressed_green(BETA, MU, np.zeros((1, nmat, 1, M, M), complex))
            chi0 = s._calc_chi0q(G, np.zeros_like(G), BETA)[0]
            _, v, _ = s._flex_compute_veff_spinful(chi0)
            return s._calc_self_energy_general(G, v, BETA)[0, :, 0]
        finally:
            logging.disable(prev)
    finally:
        shutil.rmtree(d, ignore_errors=True)


class _Gate(unittest.TestCase):

    def _check(self, types, nmat, x, tol):
        iws = 1j * np.pi * (2 * np.arange(nmat) + 1 - nmat) / BETA
        win = slice(nmat // 2 - nmat // 16, nmat // 2 + nmat // 16)
        for t in types:
            ce = G3._refine(lambda h: G3._coeff2(lambda y: _fluct_ed([(t, y)], iws[win]), h), x)
            cp = G3._refine(lambda h: G3._coeff2(lambda y: _fluct_prod([(t, y)], nmat)[win], h),
                            x)
            scale = float(np.abs(ce).max())
            offdiag = float(np.abs(ce[:, :N, N:]).max())
            err = float(np.abs(cp - ce).max()) / scale
            with self.subTest(itype=t):
                self.assertGreater(scale, G3._FLOOR)                       # anti-vacuity
                self.assertGreater(offdiag, 0.05 * scale)                  # spin mixing compared
                self.assertLess(err, tol, "{}: relative deviation {:.2e}".format(t, err))


class TestSecondOrderWithSOCSmoke(_Gate):
    """Fast variant: three types (a density-density one, the spin-flip
    Exchange and the non-spin-conserving PairLift) at ``Nmat = 512``,
    against a tolerance set at about three times the measured 6e-3."""

    def test_three_types(self):
        self._check(("CoulombIntra", "Exchange", "PairLift"), 512, 0.0125, 2.0e-2)


class TestSecondOrderWithSOC(_Gate):
    """The gate: all seven types at G3's working point and tolerance."""

    @heavy
    def test_every_type(self):
        self._check(G3._TYPES, G3._NMAT, G3._X, G3._TOL)


if __name__ == "__main__":
    unittest.main()
