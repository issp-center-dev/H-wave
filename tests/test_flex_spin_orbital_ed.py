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

``TestZeemanSecondOrder`` (gate 0-F2, #88) repeats the comparison with a
Zeeman field ``h = 0.15`` (``+h`` on spin up, ``-h`` on spin down), without
the spin-mixing block (S_z conserved) and with it, the field reaching the
solver through the Transfer file or through an Extern file with
``coeff_extern = 1``. Measured worst deviation over the seven types: 7.2e-4
(no SOC) and 8.5e-4 (SOC) at ``Nmat = 4096`` (tolerance 2e-3), both routes
equal; the up and down diagonal blocks of the ED coefficient differ by 0.20
/ 0.35 of its scale (asserted: 1e-2). The fast smoke variant (three types,
no SOC, Extern, ``Nmat = 512``) measures 5.3e-3 against the 2e-2 bound.
"""
import functools
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


def _soc_block(soc=_SOC):
    """The spin-mixing block ``F`` at amplitude ``soc`` (one fixed draw)."""
    rng = np.random.default_rng(11)
    return soc * (rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N)))


def _one_body(zeeman=0.0, soc=_SOC):
    """The one-body matrix on the spin-block index ``s * N + a``: ``H0`` in
    both spin blocks, the spin-mixing block ``F`` (amplitude ``soc``) and a
    Zeeman field (``+zeeman`` on spin up, ``-zeeman`` on spin down)."""
    F = _soc_block(soc)
    hg = np.zeros((M, M), complex)
    hg[:N, :N] = H0 + zeeman * np.eye(N)
    hg[N:, N:] = H0 - zeeman * np.eye(N)
    hg[:N, N:] = F
    hg[N:, :N] = F.conj().T
    return hg


F = _soc_block()
#: the one-body matrix of the SOC gate (no Zeeman field)
HG = _one_body()


def _hamiltonian(terms, hg=HG):
    """``G3._hamiltonian`` (``H0`` in both spin blocks) plus the rest of
    ``hg`` (spin-mixing block, Zeeman field) as a one-body operator."""
    H = G3._hamiltonian(terms)
    rest = hg.copy()
    rest[:N, :N] -= H0
    rest[N:, N:] -= H0
    with ref._quiet_matmul():
        for p in range(M):
            for q in range(M):
                if rest[p, q] != 0:
                    H = H + rest[p, q] * (_gen_op(p).conj().T @ _gen_op(q))
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


def _fluct_ed(terms, iws, hg=HG):
    """``Sigma_ED - Sigma_HF[rho_ED]`` on the spin-block index."""
    H = _hamiltonian(terms, hg)
    G = _green(H, iws)
    G0 = np.array([np.linalg.inv((iw + MU) * np.eye(M) - hg) for iw in iws])
    sig = np.linalg.inv(G0) - np.linalg.inv(G)
    hf = G3.wick_hf_sigma(M, G3._mode_terms(terms), G3._gen_to_mode(_rho(H)))
    return sig - G3._mode_to_gen(hf)[None]


def _fluct_prod(terms, nmat, hg=HG, extern=None):
    """Production ``Sigma_fluct`` of one spin-orbital FLEX map on the bare
    Green function at ``mu = MU``. ``hg`` goes into the Transfer file; with
    ``extern = h`` an Extern file in physical orbital indices carries a
    Zeeman field ``h`` (``coeff_extern = 1``) on top of it."""
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as flex_mod
    d = tempfile.mkdtemp(prefix="flex_so_ed_")
    try:
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("1 0 0\n0 1 0\n0 0 1\n{}\n".format(M) + "0 0 0\n" * M)
        rows = []
        for p in range(M):
            for q in range(M):
                if hg[p, q] != 0:
                    # the transfer file uses the interleaved index 2 a + s
                    ip = 2 * (p % N) + p // N + 1
                    iq = 2 * (q % N) + q // N + 1
                    rows.append((0, 0, 0, ip, iq, hg[p, q].real, hg[p, q].imag))
        _write_wan(os.path.join(d, "transfer.dat"), "Transfer", M, rows)
        idict = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
        if extern is not None:
            _write_wan(os.path.join(d, "extern.dat"), "Extern", N,
                       [(0, 0, 0, a + 1, a + 1, extern, 0.0) for a in range(N)])
            idict["Extern"] = "extern.dat"
        for itype, v in terms:
            _write_wan(os.path.join(d, itype + ".dat"), itype, N, G3._rows(itype, v))
            idict[itype] = itype + ".dat"
        r = read_input_k.QLMSkInput({"path_to_input": d, "interaction": idict})
        par = {"T": 1.0 / BETA, "mu": MU, "CellShape": [1, 1, 1], "SubShape": [1, 1, 1],
               "Nmat": nmat, "IterationMax": 1, "Mix": 1.0, "EPS": 1,
               "flex_hartree_fock": True}
        if extern is not None:
            par["coeff_extern"] = 1.0
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


# ---------------------------------------------------------------------------
# 0-F2: a Zeeman field, without and with the spin-mixing block (#88)
# ---------------------------------------------------------------------------

#: Zeeman field of the 0-F2 gate (+h on spin up, -h on spin down)
_ZEEMAN = 0.15
#: the smoke variant's types (as TestSecondOrderWithSOCSmoke)
_SMOKE = ("CoulombIntra", "Exchange", "PairLift")
#: smoke tolerance: 3x the measured worst deviation, but not below 2e-2
SMOKE_TOL = 2.0e-2


@functools.lru_cache(maxsize=None)
def _coeff2_ed(itype, zeeman, soc, nmat, x):
    """The refined O(x^2) coefficient of the exact fluctuation remainder of
    ``itype`` in the model ``_one_body(zeeman, soc)``, on the inner eighth
    of the ``nmat`` grid (cached per model and working point)."""
    hg = _one_body(zeeman, soc)
    iws = 1j * np.pi * (2 * np.arange(nmat) + 1 - nmat) / BETA
    win = slice(nmat // 2 - nmat // 16, nmat // 2 + nmat // 16)
    return G3._refine(lambda h: G3._coeff2(lambda y: _fluct_ed([(itype, y)], iws[win], hg), h),
                      x)


def _coeff2_prod(itype, zeeman, soc, route, nmat, x):
    """The production counterpart of :func:`_coeff2_ed`; the field reaches
    the solver through the Transfer file (``route = "transfer"``) or an
    Extern file with ``coeff_extern = 1`` (``route = "extern"``)."""
    win = slice(nmat // 2 - nmat // 16, nmat // 2 + nmat // 16)
    if route == "transfer":
        hg, ext = _one_body(zeeman, soc), None
    elif route == "extern":
        hg, ext = _one_body(0.0, soc), zeeman
    else:
        raise ValueError(route)
    return G3._refine(lambda h: G3._coeff2(
        lambda y: _fluct_prod([(itype, y)], nmat, hg, ext)[win], h), x)


def _coeff2_pair(zeeman, soc, route, nmat, x, types=G3._TYPES):
    """``(c_ed, c_prod)``, each ``(len(types), nwin, M, M)``: the second-order
    coefficients of every type, stacked in the order of ``types``."""
    c_ed = np.array([_coeff2_ed(t, zeeman, soc, nmat, x) for t in types])
    c_prod = np.array([_coeff2_prod(t, zeeman, soc, route, nmat, x) for t in types])
    return c_ed, c_prod


def _rel_dev(c_prod, c_ed):
    """Worst over the types of ``max|c_prod - c_ed| / max|c_ed|`` (G3's
    per-type max-norm relative deviation)."""
    return max(float(np.abs(p - e).max() / np.abs(e).max()) for p, e in zip(c_prod, c_ed))


def _spin_diag_blocks(c):
    """The spin-up and spin-down diagonal blocks of a stacked coefficient."""
    return c[..., :N, :N], c[..., N:, N:]


class TestZeemanSecondOrder(unittest.TestCase):
    """Gate 0-F2 (#88): with a Zeeman field, without (``soc = 0``, S_z
    conserved) and with the spin-mixing block, the second order of the
    spin-orbital FLEX map equals the exact one, the field reaching the solver
    through the Transfer file or through Extern. The up and down diagonal
    blocks of the ED coefficient differ by more than 1e-2 of its scale, so
    the field is load bearing."""

    def _asymmetric(self, c_ed):
        self.assertGreater(np.abs(c_ed).max(), G3._FLOOR)               # anti-vacuity
        up, dn = _spin_diag_blocks(c_ed)
        self.assertGreater(np.abs(up - dn).max() / np.abs(c_ed).max(), 1e-2)

    @heavy
    def test_zeeman_without_and_with_soc(self):
        for soc in (0.0, _SOC):
            for route in ("transfer", "extern"):              # production side; ED gets the same H
                with self.subTest(soc=soc, route=route):
                    c_ed, c_prod = _coeff2_pair(zeeman=_ZEEMAN, soc=soc, route=route,
                                                nmat=4096, x=0.00625)
                    dev = _rel_dev(c_prod, c_ed)
                    self.assertLessEqual(dev, 2e-3, "relative deviation {:.2e}".format(dev))
                    self._asymmetric(c_ed)

    def test_smoke(self):
        """Fast: three types, no spin-mixing block, the Extern route, at
        ``Nmat = 512``."""
        c_ed, c_prod = _coeff2_pair(zeeman=_ZEEMAN, soc=0.0, route="extern", nmat=512,
                                    x=0.0125, types=_SMOKE)
        dev = _rel_dev(c_prod, c_ed)
        self.assertLessEqual(dev, SMOKE_TOL, "relative deviation {:.2e}".format(dev))


if __name__ == "__main__":
    unittest.main()
