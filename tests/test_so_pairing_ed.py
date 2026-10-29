"""Exact diagonalization of the spin-orbital pairing vertex at second order
(issue #83, ``hwave.solver.so_pairing``).

``tests/test_so_pairing_vertex.py`` pins the bare part against an
independent first-order construction and the whole vertex against the
non-spin-orbital kernels in the spin-rotation invariant limit. Here the
oracle is an exact diagonalization with spin-orbit coupling, so the
fluctuation part is checked where no spin-resolved formula exists to
compare with. The observable is the static pair susceptibility of the
spin-orbital pair operators of a chain,

    chi_{(A, B, d'), (E, F, d)} = <Delta_(A, B, d'); Delta_(E, F, d)^dagger>,
    Delta_(A, B, d') = sum_j c_{j + d', B} c_{j, A}

(``A``, ``B`` spin-orbitals in the spin-block order ``s * norb + a``,
``TotalNED.pair_correlator`` at total momentum zero), and the production
side is the Bethe-Salpeter ladder through second order,

    chi_prod eta = 2 X[eta + P K eta + P K P K eta],   X[phi] = T sum_n G2(n) phi(n),

with ``G`` the EXACT interacting Green function of the same ED
(``TotalNED.green``), ``K`` the ``BondPairKernel`` with one (on-site)
channel and ``nfreq = nmat`` fed with ``V_bare`` (``V_inst``) and
``V_fluct`` at every bosonic frequency, both from
``so_pairing_vertex_parts`` with ``chi`` the dynamic spinful bubble of that
``G`` (the second-order ``chi``), and ``eta`` a static pair field. In the
skeleton expansion the pp-irreducible vertex through second order is the
bare antisymmetrized vertex plus the particle-hole terms, so ``chi_ED -
chi_prod = O(lambda^3)`` when every interaction is scaled by ``lambda``.
The checks compare the ``lambda^1`` and ``lambda^2`` coefficients
(five-point stencils, ``h = 0.05``) as ``tests/test_bond_pairing_ed.py``
does.

Two conventions of the comparison, both fixed at ``lambda = 0`` or by the
antisymmetry of the pair:

* ``P`` projects onto pair functions with fermionic antisymmetry,
  ``phi_AB(k, n) -> (phi_AB(k, n) - phi_BA(-k, -n)) / 2``. ``V_fluct`` is
  the direct particle-hole image only; on an antisymmetric gap the crossed
  image gives the same, so ``P K`` (direct image at weight one) is the
  physical kernel, while ``K`` alone leaves a symmetric part that a second
  ``K`` would feed back.
* The factor 2: the field ``eta`` is antisymmetric (``1/2`` on ``(E, F,
  d)``, ``-1/2`` on ``(F, E, -d)``, i.e. ``sum eta c+ c+ = Delta^dagger``),
  and both Wick contractions of ``<c c c+ c+>`` contribute ``X[eta]``.
  ``test_free_pair_susceptibility`` pins both to ``1e-6`` at
  ``lambda = 0``.

Model: two orbitals on ``L = 3`` sites, the mirror-asymmetric chain of
``tests/test_bond_pairing_ed.py`` plus a complex spin-flip hopping, so the
Green function has spin-off-diagonal blocks (a quarter of its norm) and
no spin component is conserved -- three quarters of the second-order
coefficient couple pairs of different total S_z (both asserted above
1e-1); every on-site interaction
type (CoulombIntra, CoulombInter, Hund, Exchange, PairHop with a complex
value, PairLift, Ising), declared in physical orbital indices as the
interaction files are. ``beta = 2``, ``mu = 0.8``.

Tolerance. The production side is exact up to the Matsubara truncation, an
``O(1/nmat)`` error; Richardson extrapolation over ``nmat = 128, 256``
leaves 5.8e-6 (``c1``) and 4.5e-5 (``c2``), falling to 1.5e-6 / 1.3e-5 over
``256, 512`` (measured), so the residual is the truncation and not a
formula defect. The bounds are 2e-5 (``c1``, 3.4x the measured 5.8e-6;
tightened from the 5e-5 of ``tests/test_bond_pairing_ed.py``) and 2e-4
(``c2``, 4.4x the measured 4.5e-5). The comparison is two-sided: the fluctuation term read in
the opposite orientation (``-(Gamma chi Gamma)[(d, c), (b, a)]`` at
``-q``, the same tensor only at zero bosonic frequency) misses ``c2`` by
1.9e-2, and the ladder without the fluctuation term by 9.9e-2 (both at
``nmat = 128``; asserted: 1e-2).

``TestTotalNEDAgainstSectorED`` pins the two ``TotalNED`` correlators this
gate adds (``green`` and ``pair_correlator`` with arbitrary spin labels)
to ``SectorED`` on an S_z-conserving model.

Gate 0-S2 (#88) adds a Zeeman field ``h = 0.15`` (``+h`` on spin up, ``-h``
on spin down): ``TestZeemanSpinFlipModel`` runs the comparison above on the
spin-mixing chain (TotalNED only), and ``TestZeemanSzConserving`` on the
chain without spin-flip hopping and without ``PairLift`` (``MODELS``), where
S_z is conserved, so TotalNED must also equal SectorED (``L = 2`` in the
fast gate, the 12-mode ``L = 3`` chain in the full suite). Same tolerances.
"""
import functools
import types
import unittest

import numpy as np

from hwave.solver import bubble, spin_orbital
from hwave.solver.eliashberg_bond import (ArrayBlockSource, BondPairKernel, PairVertexUniform,
                                          _vinst_as_block)
from hwave.solver.eliashberg_dynamic import calc_g2_dynamic
from hwave.solver.so_pairing import so_pairing_vertex_parts
import tests.test_sc_bond_pairing_oracle as O
from tests.ed_oracle_util import EDFixture, SectorED, TotalNED
from tests.heavy_tests import heavy
from tests.test_bond_transverse_ed import _terms_pairlift_offsite
from tests.test_so_pairing_vertex import _gamma, _onsite_decl

L, NP, BETA, MU = 3, 2, 2.0, 0.8
ND = 2 * NP                       # spin-orbitals, spin-block order s * NP + a
#: EDFixture's in-site mode index (interleaved ``2 a + s``) of spin-block index A
_MODE = spin_orbital.spinblock_order(NP)
#: Transfer rows (r, a, b): t c+_{j,a,s} c_{j+r,b,s} (the chain of test_bond_pairing_ed)
TRANSFER = {(0, 0, 0): -0.4, (0, 1, 1): 0.6, (0, 0, 1): -0.3, (0, 1, 0): -0.3,
            (1, 0, 1): -1.0, (-1, 1, 0): -1.0, (1, 0, 0): -0.5, (-1, 0, 0): -0.5}
#: spin-flip hopping (r, a, b): v c+_{j,a,up} c_{j+r,b,dn} + h.c.
SOC = {(1, 0, 1): 0.7 + 0.5j, (1, 0, 0): -0.4 + 0.6j, (0, 0, 1): 0.6j,
       (1, 1, 1): 0.3 - 0.6j}
#: the on-site interaction: (type, value) of ``_onsite_decl``
TYPES = (("CoulombIntra", 1.0), ("CoulombInter", 0.6), ("Hund", 0.5), ("Exchange", 0.5),
         ("PairHop", 0.4 + 0.3j), ("PairLift", 0.3), ("Ising", 0.2))
#: the two models of the module (#88, gate 0-S2): ``"spin_flip"`` is the
#: chain above (spin-flip hopping ``SOC``, every type of ``TYPES``), and
#: ``"no_spin_flip"`` the same chain without ``SOC`` and without the
#: S_z-changing ``PairLift``, so S_z is conserved (also with a Zeeman field)
MODELS = ("spin_flip", "no_spin_flip")

H = 0.05
LAMS = (-2 * H, -H, 0.0, H, 2 * H)
NMATS = (128, 256)
TOL_C1, TOL_C2, MISS = 2e-5, 2e-4, 1e-2

_K = 2 * np.pi * np.arange(L) / L
#: pair basis (A, B, d'): spin-orbital A at the origin, B at d'
PAIRS = [(A, B, dp) for A in range(ND) for B in range(ND) for dp in range(L)]
_INDEX = {p: i for i, p in enumerate(PAIRS)}


def _wn(nmat):
    return (2 * np.arange(nmat) + 1 - nmat) * np.pi / BETA


def _mode(fx, j, A):
    """ED mode of spin-block index ``A`` on site ``j``."""
    return fx.nd * j + int(_MODE[A])


def _antisymmetric_projector():
    """``(1 - swap) / 2`` with swap ``(A, B, d') -> (B, A, -d')``: the pair
    operators obey ``Delta_(B, A, -d') = -Delta_(A, B, d')``."""
    swap = np.zeros((len(PAIRS), len(PAIRS)))
    for i, (A, B, dp) in enumerate(PAIRS):
        swap[_INDEX[(B, A, (-dp) % L)], i] = 1.0
    return 0.5 * (np.eye(len(PAIRS)) - swap)


def _sz_change():
    """``(dSz_row, dSz_col)`` mask of the pair basis: True where the two
    pairs carry different total S_z (nonzero only through spin mixing)."""
    sz = np.array([(A // NP) + (B // NP) for (A, B, _) in PAIRS])
    return sz[:, None] != sz[None, :]


# ---------------------------------------------------------------------------
# exact diagonalization
# ---------------------------------------------------------------------------

def _types(model):
    """The interaction ``(type, value)`` entries of ``model``."""
    if model == "spin_flip":
        return TYPES
    if model == "no_spin_flip":
        return tuple((t, v) for t, v in TYPES if t != "PairLift")
    raise ValueError(model)


def _add_zeeman(fx, h1, zeeman):
    """``+zeeman`` on every spin-up mode, ``-zeeman`` on every spin-down one."""
    for j in range(fx.L):
        for a in range(fx.norb):
            h1[fx.mode(j, a, 0), fx.mode(j, a, 0)] += zeeman
            h1[fx.mode(j, a, 1), fx.mode(j, a, 1)] -= zeeman
    return h1


def _h1(fx, model="spin_flip", zeeman=0.0):
    h1 = np.zeros((fx.nmode, fx.nmode), complex)
    for (r, a, b), t in TRANSFER.items():
        for j in range(L):
            for s in range(2):
                h1[fx.mode(j, a, s), fx.mode((j + r) % L, b, s)] += t
    _types(model)                                   # validates ``model``
    _add_zeeman(fx, h1, zeeman)
    for (r, a, b), v in (SOC.items() if model == "spin_flip" else ()):
        for j in range(L):
            p, q = fx.mode(j, a, 0), fx.mode((j + r) % L, b, 1)
            h1[p, q] += v
            h1[q, p] += np.conj(v)
    return h1


def _terms(fx, t, v):
    """Quartic terms of ``_onsite_decl(t, v)`` on the chain -- the builders
    of ``test_so_pairing_vertex._onsite_terms`` (validated there against
    ``Gamma`` at first order), on ``L`` sites."""
    if t == "CoulombIntra":
        return [(fx.mode(j, a, 0), fx.mode(j, a, 0), fx.mode(j, a, 1), fx.mode(j, a, 1),
                 v * (1.0 if a == 0 else 0.8)) for j in range(L) for a in range(NP)]
    if t == "PairLift":
        return _terms_pairlift_offsite(fx, 0, 1, 0, v)
    return O._terms(t, fx, 0, 1, 0, v)


def _model(lam, model, zeeman):
    """``(fx, terms, h1)`` of ``model`` with every interaction scaled by
    ``lam`` and a Zeeman field ``zeeman``."""
    fx = EDFixture(L, NP, {}, [0.0] * NP, 1.0 / BETA, MU)
    terms = []
    for t, v in _types(model):
        terms += [(p, q, r, s, lam * c) for (p, q, r, s, c) in _terms(fx, t, v)]
    return fx, terms, _h1(fx, model, zeeman)


def _ed(lam, model="spin_flip", zeeman=0.0):
    """``(chi, G)``: the pair susceptibility in the ``PAIRS`` basis and
    ``G[n, A, (j, B)] = G_{(0, A), (j, B)}(i w_n)`` on the ``max(NMATS)``
    grid (cached per ``(lam, model, zeeman)``). Uncached itself: it passes
    every argument positionally to :func:`_ed_cached`, so ``_ed(lam)`` and
    ``_ed(lam, "spin_flip", 0.0)`` share one cache entry (``lru_cache``
    keys the two spellings differently)."""
    return _ed_cached(float(lam), str(model), float(zeeman))


@functools.lru_cache(maxsize=None)
def _ed_cached(lam, model, zeeman):
    fx, terms, h1 = _model(lam, model, zeeman)
    ed = TotalNED(fx, terms, h1=h1)
    chi = ed.pair_correlator(
        [(dp, B % NP, B // NP, A % NP, A // NP) for (A, B, dp) in PAIRS], qs=[0])[0]
    G = ed.green(1j * _wn(max(NMATS)), rows=[_mode(fx, 0, A) for A in range(ND)],
                 cols=[_mode(fx, j, B) for j in range(L) for B in range(ND)])
    return chi, G


def _green_sc(G, nmat):
    """``G_AB(k, i w_n) = sum_r e^{ikr} G_{(0 A),(r B)}`` in the sc layout
    ``(ND, ND, L, 1, 1, nmat)`` on the centre ``nmat`` points of the ED grid."""
    off = (G.shape[0] - nmat) // 2
    Gr = G[off:off + nmat].reshape(nmat, ND, L, ND)
    g = np.einsum("narb,kr->abkn", Gr, np.exp(1j * np.outer(_K, np.arange(L))))
    return g[:, :, :, None, None, :].copy()


# ---------------------------------------------------------------------------
# production ladder
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=None)
def _gamma_onsite(model="spin_flip"):
    """``Gamma`` (L, ND^2, ND^2) of every interaction entry of ``model``,
    from the spin-orbital RPA reader (q-independent: on-site only)."""
    decl = {}
    for t, v in _types(model):
        decl.update(_onsite_decl(t, v))
    gam = _gamma(decl)[0]
    return np.broadcast_to(gam, (L,) + gam.shape).copy()


def _opposite(V):
    """``V[d, c, b, a](-q)`` at ``[a, b, c, d](q)``."""
    return np.roll(V[..., ::-1, :, :], 1, axis=4).transpose(3, 2, 1, 0, 4, 5, 6)


def _kernel(g, lam, nmat, fluct="direct", model="spin_flip"):
    """The one-channel ``BondPairKernel`` of ``V_bare + V_fluct(q, i nu_m)``;
    ``fluct``: ``"direct"`` (as implemented), ``"opposite"`` (the opposite
    reading) or ``"none"``."""
    gam = lam * _gamma_onsite(model)
    gk = g[:, :, :, 0, 0, :].transpose(3, 2, 0, 1)[None].copy()
    chi0 = bubble.dense_bubble(gk, None, BETA, spatial_shape=(L, 1, 1), scheme="general")[0]
    nd2 = ND * ND
    blocks = np.zeros((nmat, L, nd2, nd2), complex)
    V_bare = None
    for m in range(nmat):             # bosonic index m: nu_m = (2 m - nmat) pi / beta
        V_bare, V_fl = so_pairing_vertex_parts(gam, chi0[m].reshape(L, nd2, nd2), ND,
                                               (L, 1, 1))
        if fluct == "direct":
            blocks[m] = _vinst_as_block(V_fl, nd2, L)
        elif fluct == "opposite":
            blocks[m] = _vinst_as_block(_opposite(V_fl), nd2, L)
    view = types.SimpleNamespace(n_channels=1, delta_r=[(0, 0, 0)])
    vert = PairVertexUniform(source=ArrayBlockSource({"G": blocks}, nd2), slot="G",
                             B=1, nd=nd2, nmat=nmat, nvol=L)
    return BondPairKernel(vert, calc_g2_dynamic(g, BETA), view, xp=np,
                          spatial_shape=(L, 1, 1), norb=ND, beta=BETA, nfreq=nmat,
                          V_inst=V_bare, residency="host", host_cap=10 ** 12,
                          device_cap=10 ** 12)


def _antisym(phi):
    """``(phi_AB(k, n) - phi_BA(-k, -n)) / 2`` on ``(ND, ND, L, nmat)``."""
    ref = np.roll(phi[:, :, ::-1, ::-1], 1, axis=2).transpose(1, 0, 2, 3)
    return 0.5 * (phi - ref)


def _pair_sum(g, phi_w):
    """``X_AB(k) = T sum_n G_AE(k, w_n) G_BF(-k, -w_n) phi_EF(k, n)``."""
    gm = np.roll(g[:, :, ::-1, :, :, ::-1], 1, axis=2)[:, :, :, 0, 0]
    return np.einsum("cekn,dfkn,efkn->cdk", g[:, :, :, 0, 0], gm, phi_w) / BETA


def _production_chi(g, lam, nmat, **kw):
    """The ladder through ``K K`` in the ``PAIRS`` basis: column ``(E, F,
    d)`` is the response to the antisymmetric static field of
    ``Delta_(E, F, d)^dagger``."""
    out = np.zeros((len(PAIRS), len(PAIRS)), complex)
    phase = np.exp(1j * np.outer(_K, np.arange(L)))
    K = _kernel(g, lam, nmat, **kw)
    for col, (E, F, dd) in enumerate(PAIRS):
        Dr = np.zeros((ND, ND, L), complex)
        Dr[E, F, dd] += 0.5
        Dr[F, E, (-dd) % L] -= 0.5
        phi = np.einsum("abd,kd->abk", Dr, phase)                     # e^{+ikd}
        term = np.ascontiguousarray(np.broadcast_to(phi[..., None], phi.shape + (nmat,)))
        total = term.copy()
        for _ in range(2):
            term = _antisym(K.matvec(term[:, :, :, None, None, :].ravel()).reshape(
                ND, ND, L, nmat))
            total = total + term
        X = 2.0 * _pair_sum(g, total)
        for row, (A, B, dp) in enumerate(PAIRS):
            out[row, col] += (X[A, B] * phase[:, dp].conj()).mean()   # e^{-ikd'}
    return out


@functools.lru_cache(maxsize=None)
def _production(nmat, fluct="direct", lams=LAMS, model="spin_flip", zeeman=0.0):
    return tuple(_production_chi(_green_sc(_ed(lam, model, zeeman)[1], nmat), lam, nmat,
                                 fluct=fluct, model=model)
                 for lam in lams)


def _richardson(fluct="direct", lams=LAMS, model="spin_flip", zeeman=0.0):
    lo, hi = (_production(n, fluct, lams, model, zeeman) for n in NMATS)
    return tuple(2 * b - a for a, b in zip(lo, hi))


def _stencil(f):
    """``(c1, c2)`` of ``f(lam) = c0 + c1 lam + c2 lam^2 + ...`` on ``LAMS``."""
    fm2, fm1, f0, f1, f2 = f
    c1 = (8 * (f1 - fm1) - (f2 - fm2)) / (12 * H)
    c2 = (-f2 + 16 * f1 - 30 * f0 + 16 * fm1 - fm2) / (24 * H * H)
    return c1, c2


def _ed_series(model="spin_flip", zeeman=0.0):
    return tuple(_ed(lam, model, zeeman)[0] for lam in LAMS)


def _errors(want, got):
    """``{order: relative deviation}`` of the coefficients, both sides
    projected onto the antisymmetric pair space (the ED side lives there
    exactly)."""
    P = _antisymmetric_projector()
    out = {}
    for order, (a, b) in enumerate(zip(_stencil(want), _stencil(got)), 1):
        A, B = P @ a @ P, P @ b @ P
        out[order] = np.linalg.norm(A - B) / np.linalg.norm(A)
    return out


def _coefficients(model="spin_flip", zeeman=0.0):
    """``(dev_c1, dev_c2)``: the relative deviations of the Richardson
    production ``lambda^1`` and ``lambda^2`` coefficients from the ED ones,
    for ``model`` with a Zeeman field ``zeeman``."""
    dev = _errors(_ed_series(model, zeeman), _richardson(model=model, zeeman=zeeman))
    return dev[1], dev[2]


class TestSOPairingVertexAgainstED(unittest.TestCase):

    @heavy
    def test_first_and_second_order(self):
        want = _ed_series()
        # the spin-orbit coupling is load bearing: the Green function has
        # sizeable spin-off-diagonal blocks, and the second-order pair
        # susceptibility couples pairs of different total S_z
        g = _green_sc(_ed(0.0)[1], NMATS[0])
        self.assertGreater(np.linalg.norm(g[:NP, NP:]), 1e-1 * np.linalg.norm(g))
        c2 = _stencil(want)[1]
        self.assertGreater(np.linalg.norm(c2[_sz_change()]), 1e-1 * np.linalg.norm(c2))
        c1, c2 = _coefficients()
        self.assertLess(c1, TOL_C1, "order 1: {:.3e}".format(c1))
        self.assertLess(c2, TOL_C2, "order 2: {:.3e}".format(c2))

    @heavy
    def test_free_pair_susceptibility(self):
        """``lambda = 0``: the pair-space mapping, the factor 2 and the
        bubble ``X`` reproduce the ED pair susceptibility."""
        want = _ed(0.0)[0]
        got = _richardson(lams=(0.0,))[0]
        dev = np.linalg.norm(got - want) / np.linalg.norm(want)
        self.assertLess(dev, 1e-6, "free pair susceptibility: {:.3e}".format(dev))

    @heavy
    def test_opposite_reading_is_missed(self):
        """The fluctuation term read in the opposite orientation (equal to
        the implemented one at zero bosonic frequency only) misses the ED at
        O(lambda^2)."""
        dev = _errors(_ed_series(), _production(NMATS[0], "opposite"))[2]
        self.assertGreater(dev, MISS, "{:.3e}".format(dev))

    @heavy
    def test_second_order_needs_the_fluctuation_vertex(self):
        """The ladder with ``V_bare`` alone misses the ED at O(lambda^2)."""
        dev = _errors(_ed_series(), _production(NMATS[0], "none"))[2]
        self.assertGreater(dev, MISS, "{:.3e}".format(dev))


# ---------------------------------------------------------------------------
# TotalNED's Green function and pair correlator against SectorED
# ---------------------------------------------------------------------------

def _sz_fixture_l2():
    """``(fx, terms)``: the S_z-conserving two-site, two-orbital fixture of
    :class:`TestTotalNEDAgainstSectorED`."""
    fx = EDFixture(2, 2, {(0, 0): -1.0, (1, 1): -0.7, (0, 1): -0.3, (1, 0): -0.2},
                   [0.1, -0.2], 1.0 / 1.5, 0.15)
    terms = [(fx.mode(j, a, 0), fx.mode(j, a, 0), fx.mode(j, a, 1), fx.mode(j, a, 1), 1.0)
             for j in range(fx.L) for a in range(fx.norb)]
    for t, v in (("CoulombInter", 0.5), ("Hund", 0.3), ("Exchange", 0.25),
                 ("PairHop", 0.2), ("Ising", 0.1)):
        terms += O._terms(t, fx, 0, 1, 0, v)
    return fx, terms


class TestTotalNEDAgainstSectorED(unittest.TestCase):
    """On an S_z-conserving model (two sites, two orbitals, inter-orbital
    hopping and S_z-conserving interactions) the total-N engine reproduces
    the (N_up, N_dn) engine."""

    @classmethod
    def setUpClass(cls):
        fx, terms = _sz_fixture_l2()
        cls.fx = fx
        cls.sector = SectorED(fx, terms)
        cls.total = TotalNED(fx, terms)

    def test_green_matches_sector_ed(self):
        iws = 1j * np.pi / self.fx.beta * np.array([-7.0, -1.0, 1.0, 3.0, 11.0])
        want = self.sector.green(iws)
        got = self.total.green(iws)
        self.assertGreater(np.abs(want).max(), 1e-2)
        self.assertLess(np.abs(got - want).max(), 1e-12 * np.abs(want).max())
        rows = [self.fx.mode(0, 1, 1), self.fx.mode(1, 0, 0)]
        cols = [self.fx.mode(1, 1, 1), self.fx.mode(0, 0, 0), self.fx.mode(0, 1, 1)]
        sub = self.total.green(iws, rows=rows, cols=cols)
        self.assertLess(np.abs(sub - self.sector.green(iws, rows=rows, cols=cols)).max(),
                        1e-12 * np.abs(want).max())

    def test_pair_correlator_matches_sector_ed(self):
        fx = self.fx
        chans = [(R, a, b) for R in range(fx.L) for a in range(fx.norb) for b in range(fx.norb)]
        want = self.sector.pair_correlator(chans)
        got = self.total.pair_correlator([(R, a, 1, b, 0) for (R, a, b) in chans])
        self.assertGreater(np.abs(want).max(), 1e-2)
        self.assertLess(np.abs(got - want).max(), 1e-12 * np.abs(want).max())
        got1 = self.total.pair_correlator([(R, a, 1, b, 0) for (R, a, b) in chans], qs=[1])
        self.assertLess(np.abs(got1[0] - want[1]).max(), 1e-12 * np.abs(want).max())


# ---------------------------------------------------------------------------
# 0-S2: a Zeeman field (#88)
# ---------------------------------------------------------------------------

#: Zeeman field of the 0-S2 gates (+h on spin up, -h on spin down)
_ZEEMAN = 0.15


def _zeeman_split(g):
    """``||G_upup - G_dndn|| / ||G||`` of an sc-layout Green function."""
    return np.linalg.norm(g[:NP, :NP] - g[NP:, NP:]) / np.linalg.norm(g)


@functools.lru_cache(maxsize=None)
def _pair_corr(engine, zeeman, L):
    """The opposite-spin pair correlator ``(L, n, n)`` (SectorED's channels
    ``(R, a, b)``, TotalNED's ``(R, a, 1, b, 0)``) of an S_z-conserving model
    with a Zeeman field ``zeeman``, from ``engine`` (``"TotalNED"`` or
    ``"SectorED"``). ``L = 2``: the fixture of
    :class:`TestTotalNEDAgainstSectorED`; ``L = 3``: the ``"no_spin_flip"``
    chain at the largest stencil coupling ``lambda = 2 h``. Cached per
    ``(engine, zeeman, L)``."""
    if L == 2:
        fx, terms = _sz_fixture_l2()
        h1 = _add_zeeman(fx, fx.build_h1().copy(), zeeman)
    elif L == 3:
        fx, terms, h1 = _model(LAMS[-1], "no_spin_flip", zeeman)
    else:
        raise ValueError(L)
    chans = [(R, a, b) for R in range(fx.L) for a in range(fx.norb) for b in range(fx.norb)]
    if engine == "SectorED":
        return SectorED(fx, terms, h1=h1).pair_correlator(chans)
    if engine == "TotalNED":
        return TotalNED(fx, terms, h1=h1).pair_correlator([(R, a, 1, b, 0) for (R, a, b) in chans])
    raise ValueError(engine)


class TestZeemanSpinFlipModel(unittest.TestCase):
    """Gate 0-S2 (i): the spin-flip-hopping chain of
    :class:`TestSOPairingVertexAgainstED` with a Zeeman field (TotalNED only:
    S_z is not conserved) -- same comparison and tolerances."""

    @heavy
    def test_first_and_second_order(self):
        g = _green_sc(_ed(0.0, "spin_flip", _ZEEMAN)[1], NMATS[0])
        self.assertGreater(_zeeman_split(g), 1e-2)                  # the field is load bearing
        self.assertGreater(np.linalg.norm(g[:NP, NP:]), 1e-1 * np.linalg.norm(g))
        c1, c2 = _coefficients(model="spin_flip", zeeman=_ZEEMAN)
        self.assertLess(c1, TOL_C1, "order 1: {:.3e}".format(c1))
        self.assertLess(c2, TOL_C2, "order 2: {:.3e}".format(c2))


class TestZeemanSzConserving(unittest.TestCase):
    """Gate 0-S2 (ii): without spin-flip hopping and without ``PairLift``
    the Zeeman-split model conserves S_z, so the total-N engine must agree
    with the (N_up, N_dn) engine, and the production ladder with the ED."""

    def test_totaled_equals_sectored(self):
        want = _pair_corr("SectorED", _ZEEMAN, 2)
        self.assertGreater(np.abs(want).max(), 1e-2)
        # the field is load bearing
        self.assertGreater(np.abs(want - _pair_corr("SectorED", 0.0, 2)).max(), 1e-3)
        np.testing.assert_allclose(_pair_corr("TotalNED", _ZEEMAN, 2), want, atol=1e-12)

    @heavy
    def test_totaled_equals_sectored_L3(self):
        want = _pair_corr("SectorED", _ZEEMAN, 3)
        self.assertGreater(np.abs(want).max(), 1e-2)
        np.testing.assert_allclose(_pair_corr("TotalNED", _ZEEMAN, 3), want, atol=1e-12)

    @heavy
    def test_first_and_second_order(self):
        g = _green_sc(_ed(0.0, "no_spin_flip", _ZEEMAN)[1], NMATS[0])
        self.assertGreater(_zeeman_split(g), 1e-2)                  # the field is load bearing
        self.assertLess(np.linalg.norm(g[:NP, NP:]), 1e-12 * np.linalg.norm(g))   # S_z conserved
        c1, c2 = _coefficients(model="no_spin_flip", zeeman=_ZEEMAN)
        self.assertLess(c1, TOL_C1, "order 1: {:.3e}".format(c1))
        self.assertLess(c2, TOL_C2, "order 2: {:.3e}".format(c2))


if __name__ == "__main__":
    unittest.main()
