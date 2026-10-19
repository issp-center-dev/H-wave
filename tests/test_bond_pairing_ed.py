"""Exact diagonalization of the FULL bond-resolved dynamic pairing vertex
(bare + fluctuation) on an inter-orbital, mirror-asymmetric chain (issue #204).

``tests/test_bond_pairing_placement.py`` pins the leg geometry of
``eliashberg_bond.BondPairKernel`` against a static spin-orbital vertex built
by hand. Here the oracle is an exact diagonalization, so nothing on the
reference side is a vertex at all: the observable is the static pair
susceptibility of the up-down pair operators of the chain,

    chi_{(0 c),(d' d); (j e),(j + dd f)} summed over j
        = Xpp[q = 0] of ``SectorED.pair_correlator`` (``tests/ed_oracle_util.py``),

and the production side is the Bethe-Salpeter ladder of the bond pair kernel,

    chi_prod eta = X[eta + K eta + K K eta],   X[phi] = T sum_n G2(n) phi(n),

with ``G`` the EXACT interacting Green function of the same ED
(``SectorED.green``), the bond vertex ``Gamma = c_s S chibar S + c_c C chibar C
+ V_inst`` (the second-order vertex: ``chi_s = chi_c = chibar``, the dynamic
bond bubble of that ``G``), and a static pair field ``eta`` of definite spin
parity. In the skeleton expansion the pp-irreducible vertex through second
order is exactly that, so ``chi_ED - chi_prod = O(lambda^3)`` when every
interaction is scaled by ``lambda``: the self-energy is inside ``G`` on both
sides. The checks compare the ``lambda^1`` and ``lambda^2`` coefficients
(five-point stencils, ``h = 0.05``) in the singlet and the triplet sector.

Model (the one of the placement test, on ``L = 3``): two orbitals, a hopping
that is neither mirror- nor transpose-symmetric (the inter-orbital bond with
orbital 0 at the origin and 1 at ``+x`` is strong, its mirror carries no
hopping, plus an on-site inter-orbital term), an on-site ``U`` on both
orbitals and an inter-orbital CoulombInter table that is not symmetric under
``R -> -R`` (0.4 on one bond, 0.1 on its mirror). ``beta = 1``, ``mu = 0.15``.

Tolerance. The only approximation on the production side is the Matsubara
truncation, an ``O(1/nmat)`` error (measured: ``c1`` 2.1e-3 / ``c2`` 1.5e-3
relative at ``nmat = 256``, halving with each doubling to 2.6e-4 / 1.8e-4 at
2048). Richardson extrapolation over ``nmat = 128, 256`` leaves at most
4.5e-6 (``c1``) and 2.3e-5 (``c2``); the asserted bounds are 5e-5 and 2e-4.
Each comparison is two-sided: the ED of the OPPOSITE placement of the
two-body rows, the ladder without the fluctuation vertex, and the ladder with
the singlet and triplet spin coefficients interchanged each miss by at least
2.7e-2 (asserted: 1e-2) in every sector they touch.
"""
import functools
import unittest

import numpy as np

from hwave.solver import bond_channels as bc
from hwave.solver import bubble
from hwave.solver.eliashberg_bond import (COEFF, ArrayBlockSource, BondPairKernel,
                                          PairVertexUniform, instantaneous_vertex)
from hwave.solver.eliashberg_dynamic import calc_g2_dynamic
from tests.ed_oracle_util import EDFixture, SectorED, canonical_density_terms
from tests.eliashberg_bond_fixtures import _channel0_from_gate, _tables_from_pairs
from tests.heavy_tests import heavy

L, NORB, BETA, MU, U = 3, 2, 1.0, 0.15, 1.0
ND = NORB * NORB
#: Transfer rows (r, a, b): t c+_{j,a} c_{j+r,b}.
TRANSFER = {(0, 0, 0): -0.4, (0, 1, 1): 0.6, (0, 0, 1): -0.3, (0, 1, 0): -0.3,
            (1, 0, 1): -1.0, (-1, 1, 0): -1.0, (1, 0, 0): -0.5, (-1, 0, 0): -0.5}
#: CoulombInter rows (r, a, b), each declared with its partner (-r, b, a).
ROWS = {(1, 0, 1): 0.4, (-1, 1, 0): 0.4, (-1, 0, 1): 0.1, (1, 1, 0): 0.1}

H = 0.05
LAMS = (-2 * H, -H, 0.0, H, 2 * H)
NMATS = (128, 256)
TOL_C1, TOL_C2, MISS = 5e-5, 2e-4, 1e-2

_K = 2 * np.pi * np.arange(L) / L
#: pair-space basis (c, d, d'): orbital c (up) at the origin, d (down) at d'
PAIRS = [(c, d, dp) for c in range(NORB) for d in range(NORB) for dp in range(L)]
_INDEX = {p: i for i, p in enumerate(PAIRS)}


def _wn(nmat):
    return (2 * np.arange(nmat) + 1 - nmat) * np.pi / BETA


def _parity_projectors():
    """Singlet / triplet: even / odd under the exchange of the two pair legs,
    ``(c, d, d') -> (d, c, -d')``."""
    swap = np.zeros((len(PAIRS), len(PAIRS)))
    for i, (c, d, dp) in enumerate(PAIRS):
        swap[_INDEX[(d, c, (-dp) % L)], i] = 1.0
    eye = np.eye(len(PAIRS))
    return {"singlet": 0.5 * (eye + swap), "triplet": 0.5 * (eye - swap)}


# ---------------------------------------------------------------------------
# exact diagonalization
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=None)
def _ed(lam, sign=+1):
    """``(chi, G)``: the pair susceptibility in the ``PAIRS`` basis and the
    up-spin Green function ``G[n, a, (j, b)]`` on the ``max(NMATS)`` grid.
    ``sign = +1`` places a row ``(r, a, b)`` on the bond of the Transfer row
    ``(r, a, b)`` (``a`` at ``j``, ``b`` at ``j + r``), ``sign = -1`` on the
    opposite bond."""
    fx = EDFixture(L, NORB, {}, [0.0] * NORB, 1.0 / BETA, MU)
    h1 = np.zeros((fx.nmode, fx.nmode), complex)
    for (r, a, b), t in TRANSFER.items():
        for j in range(L):
            for s in range(2):
                h1[fx.mode(j, a, s), fx.mode((j + r) % L, b, s)] += t
    # canonical_density_terms puts (a, b, R) on (a at j + R, b at j); a row and
    # its partner name the same density product, each carrying half of it
    pairs = [(a, a, 0, lam * U) for a in range(NORB)]
    pairs += [(a, b, -sign * r, lam * v / 2) for (r, a, b), v in ROWS.items()]
    ed = SectorED(fx, canonical_density_terms(fx, pairs), h1=h1)
    # Delta_{R,a,b} = sum_j c_{j+R,a,dn} c_{j,b,up}: the pair (c=b up at 0, d=a down at R)
    chi = ed.pair_correlator([(dp, d, c) for (c, d, dp) in PAIRS])[0]
    G = ed.green(1j * _wn(max(NMATS)), rows=[fx.mode(0, a, 0) for a in range(NORB)],
                 cols=[fx.mode(j, b, 0) for j in range(L) for b in range(NORB)])
    return chi, G


def _green_sc(G, nmat):
    """``G_ab(k, i w_n) = sum_r e^{ikr} G_{(0 a),(r b)}`` in the sc layout
    ``(norb, norb, L, 1, 1, nmat)``, on the ``nmat`` grid (the centre of the
    ED grid)."""
    off = (G.shape[0] - nmat) // 2
    Gr = G[off:off + nmat].reshape(nmat, NORB, L, NORB)          # (n, a, r, b)
    g = np.einsum("narb,kr->abkn", Gr, np.exp(1j * np.outer(_K, np.arange(L))))
    return g[:, :, :, None, None, :].copy()


# ---------------------------------------------------------------------------
# production ladder
# ---------------------------------------------------------------------------

def _kernel(g, lam, eta, nmat, fluct=True, coeff=None):
    inter = {}
    for (r, a, b), v in ROWS.items():
        inter.setdefault(("CoulombInter", (r, 0, 0)), np.zeros((NORB, NORB)))[a, b] = lam * v
    whole, onsite, offsite = _tables_from_pairs(inter, lam * U, NORB)
    S0, C0 = _channel0_from_gate(whole, onsite, offsite, NORB, (L, 1, 1))
    topo = bc.resolve_bond_topology(whole, np.eye(3), NORB,
                                    active_types=bc._LONGITUDINAL_ACTIVE_TYPES)
    view = bc.BondSetView(topo)
    S = bc.build_sc_bond_channel(topo, S0, "S")
    C = bc.build_sc_bond_channel(topo, C0, "C")
    green = g[:, :, :, 0, 0, :].transpose(3, 2, 0, 1)[None].copy()
    chib = np.zeros((nmat, L) + S.shape[1:], complex)
    for (m, mp), blk in bubble._iter_bond_dynamic(green, None, BETA, view,
                                                  spatial_shape=(L, 1, 1)):
        chib[:, :, m * ND:(m + 1) * ND, mp * ND:(mp + 1) * ND] = blk
    cs, cc = (coeff or COEFF)[eta] if fluct else (0.0, 0.0)
    gam = cs * S[None] @ chib @ S[None] + cc * C[None] @ chib @ C[None]
    vert = PairVertexUniform(source=ArrayBlockSource({"G": gam}, ND), slot="G",
                             B=view.n_channels, nd=ND, nmat=nmat, nvol=L)
    return BondPairKernel(vert, calc_g2_dynamic(g, BETA), view, xp=np,
                          spatial_shape=(L, 1, 1), norb=NORB, beta=BETA, nfreq=nmat,
                          V_inst=instantaneous_vertex(S, C, ND, eta, (L, 1, 1)),
                          residency="host", host_cap=10 ** 12, device_cap=10 ** 12)


def _pair_sum(g, phi_w):
    """``X_cd(k) = T sum_n G_ce(k, w_n) G_df(-k, -w_n) phi_ef(k, n)``."""
    gm = np.roll(g[:, :, ::-1, :, :, ::-1], 1, axis=2)[:, :, :, 0, 0]
    return np.einsum("cekn,dfkn,efkn->cdk", g[:, :, :, 0, 0], gm, phi_w) / BETA


def _production_chi(g, lam, nmat, **kw):
    """The ladder through ``K K`` in the ``PAIRS`` basis: column ``(e, f, dd)``
    is the response ``kappa_{(0 c),(d' d)}`` to the pair field
    ``eta_{(j e),(j + dd f)} = 1``, split into its two spin-parity parts."""
    out = np.zeros((len(PAIRS), len(PAIRS)), complex)
    phase = np.exp(1j * np.outer(_K, np.arange(L)))
    for eta, sgn in (("singlet", 1.0), ("triplet", -1.0)):
        K = _kernel(g, lam, eta, nmat, **kw)
        for col, (e, f, dd) in enumerate(PAIRS):
            Dr = np.zeros((NORB, NORB, L), complex)
            Dr[e, f, dd] += 0.5
            Dr[f, e, (-dd) % L] += 0.5 * sgn
            phi = np.einsum("abd,kd->abk", Dr, phase)                 # e^{+ikd}
            term = np.ascontiguousarray(np.broadcast_to(phi[..., None], phi.shape + (nmat,)))
            total = term.copy()
            for _ in range(2):
                term = K.matvec(term[:, :, :, None, None, :].ravel()).reshape(
                    NORB, NORB, L, nmat)
                total = total + term
            X = _pair_sum(g, total)
            for row, (c, d, dp) in enumerate(PAIRS):
                out[row, col] += (X[c, d] * phase[:, dp].conj()).mean()   # e^{-ikd'}
    return out


@functools.lru_cache(maxsize=None)
def _production(nmat, sign=+1, fluct=True, swapped=False):
    """The production ladder at every ``LAMS`` point, from the Green
    function of the ED with placement ``sign``."""
    coeff = {"singlet": COEFF["triplet"], "triplet": COEFF["singlet"]} if swapped else None
    return tuple(_production_chi(_green_sc(_ed(lam, sign)[1], nmat), lam, nmat,
                                 fluct=fluct, coeff=coeff) for lam in LAMS)


def _richardson():
    lo, hi = (_production(n) for n in NMATS)
    return tuple(2 * b - a for a, b in zip(lo, hi))


def _coefficients(f):
    """``(c1, c2)`` of ``f(lam) = c0 + c1 lam + c2 lam^2 + ...`` on ``LAMS``."""
    fm2, fm1, f0, f1, f2 = f
    c1 = (8 * (f1 - fm1) - (f2 - fm2)) / (12 * H)
    c2 = (-f2 + 16 * f1 - 30 * f0 + 16 * fm1 - fm2) / (24 * H * H)
    return c1, c2


def _ed_series(sign=+1):
    return tuple(_ed(lam, sign)[0] for lam in LAMS)


def _sector_errors(want, got):
    """``{(order, eta): relative deviation}`` of the coefficients."""
    out = {}
    for order, (a, b) in enumerate(zip(_coefficients(want), _coefficients(got)), 1):
        for eta, P in _parity_projectors().items():
            A, B = P @ a @ P, P @ b @ P
            out[(order, eta)] = np.linalg.norm(A - B) / np.linalg.norm(A)
    return out


class TestPairingVertexAgainstED(unittest.TestCase):

    @heavy
    def test_free_pair_susceptibility(self):
        """``lambda = 0``: the pair-space mapping and the bubble ``X``
        reproduce the ED pair susceptibility (fixes the basis and the
        normalization the coefficient checks rely on)."""
        want = _ed(0.0)[0]
        got = _richardson()[LAMS.index(0.0)]
        dev = np.linalg.norm(got - want) / np.linalg.norm(want)
        self.assertLess(dev, 1e-6, "free pair susceptibility: {:.3e}".format(dev))

    @heavy
    def test_first_and_second_order(self):
        errs = _sector_errors(_ed_series(), _richardson())
        for (order, eta), dev in sorted(errs.items()):
            tol = TOL_C1 if order == 1 else TOL_C2
            self.assertLess(dev, tol, "order {} {}: {:.3e}".format(order, eta, dev))

    @heavy
    def test_opposite_placement_is_missed(self):
        """The production ladder (documented placement) against the ED of the
        opposite placement, both from that ED's own Green function."""
        errs = _sector_errors(_ed_series(-1), _production(NMATS[0], sign=-1))
        for (order, eta), dev in sorted(errs.items()):
            self.assertGreater(dev, MISS, "order {} {}: {:.3e}".format(order, eta, dev))

    @heavy
    def test_second_order_needs_the_fluctuation_vertex(self):
        """The ladder with ``V_inst`` alone misses the ED at O(lambda^2)."""
        errs = _sector_errors(_ed_series(), _production(NMATS[0], fluct=False))
        for eta in ("singlet", "triplet"):
            dev = errs[(2, eta)]
            self.assertGreater(dev, MISS, "{}: {:.3e}".format(eta, dev))

    @heavy
    def test_spin_coefficients_are_not_interchangeable(self):
        """Singlet and triplet ``(c_s, c_c)`` swapped miss the ED at O(lambda^2)."""
        errs = _sector_errors(_ed_series(), _production(NMATS[0], swapped=True))
        for eta in ("singlet", "triplet"):
            dev = errs[(2, eta)]
            self.assertGreater(dev, MISS, "{}: {:.3e}".format(eta, dev))


if __name__ == "__main__":
    unittest.main()
