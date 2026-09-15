"""Independent real-space second-order oracle (spec 2.1) and gates G1 (ii),
(iii) and G2 (a): the antisymmetrised skeleton enumerated over
site-orbital-spin records on the 4x4 torus with an operator table written
HERE (independent of hwave.solver.second_order._records), the local/dropped
split by the weight rule of spec 2.4, and the exhaustive comparison against
production for every pair of coupling types on an asymmetric bond.

G2 (a) is compared at TWO levels, because a self-energy comparison alone
cannot see everything a kernel comparison can: the transport at one fixed
``G`` maps ``nmat nvol nd^2`` kernel components onto ``nmat nvol norb^2``
self-energy ones, so it has a nullspace.
:func:`oracle_sigma2` is the self-energy level and :func:`oracle_w2` the
KERNEL level (compared element by element with ``dense_w2``); both are
derived from the same records and the same weight rule, and
:meth:`TestOracle.test_direct_comparison_catches_a_transport_nullspace_mutation`
exhibits a kernel perturbation that only the second one catches.

Three conventions of the discrete representation are pinned by the two G1
gates before the oracle is used on the production kernel, and each is
documented where it is implemented below: the mirrored-row half weight
(the reader's reversal closure emits both ordered entries of a bond), the
two factors 1/beta of a second-order diagram, and the direction of the
reversed propagator ``rev[d] = G(-d, -tau)``."""
import itertools
import unittest

import numpy as np

from hwave.solver import matsubara as _ms, backend as _bk
from tests.test_second_order_factors import _split_for

_SHAPE = (4, 4, 1)
_NVOL = 16
_NMAT = 8
_BETA = 0.5
UP, DN = 0, 1


def _site(x, y, z=0, shape=_SHAPE):
    """Flat site index of a lattice point on ``shape`` (the solver's
    ``(nx, ny, nz)`` order, x-major). Defaults to the 4x4 torus the G1/G2
    gates of this module are written on; the chain gate
    (``tests/test_flex_second_order_ed_chain.py``) passes ``(L, 1, 1)``."""
    nx, ny, nz = shape
    return ((x % nx) * ny + (y % ny)) * nz + (z % nz)


def _coords(i, shape=_SHAPE):
    """Inverse of :func:`_site`: the ``(x, y, z)`` of a flat site index."""
    nx, ny, nz = shape
    x, rem = divmod(i, ny * nz)
    y, z = divmod(rem, nz)
    return x, y, z


def _g(i, a, s, norb):
    """site-orbital-spin index: (i, a, s) -> (i * 2 + s) * norb + a."""
    return (i * 2 + s) * norb + a


#: Types whose transposed row (b, a) is the SAME operator as (a, b). The
#: reader's reversal closure emits BOTH ordered entries of every declared
#: bond, so each ordered row is HALF the bond: a bond declared in both
#: directions with v couples with v, not 2 v. ``CoulombIntra`` has no
#: transposed partner and ``PairHop``'s transposed row is a different
#: operator (its Hermitian partner), so neither is halved. This is the
#: Hamiltonian the solver implements (the k-space Hartree-Fock kernel's
#: convention); gate G1 (iii) below pins it against the Phase B
#: enumeration of ``tests/test_flex_bond_sopt.py``, which reads the same
#: declarations independently.
_MIRRORED_ROW_TYPES = ("CoulombInter", "Hund", "Ising", "Exchange", "PairLift")


def oracle_records(rows_by_type, norb, shape=_SHAPE):
    """(p, q, r, s, x, offsite) records of x c+_p c+_q c_s c_r, one per
    monomial of the documented operator per DECLARED ROW (spec 2.2 table,
    written independently: n_i n_j expansions, Exchange both directions,
    PairHop one orientation, PairLift both placements, Hund with the code's
    minus, Ising density-difference form), at half weight for the mirrored
    types (:data:`_MIRRORED_ROW_TYPES`).

    ``shape`` is the lattice the records are enumerated on, defaulting to
    this module's 4x4 torus; the chain exact-diagonalization gate passes
    ``(L, 1, 1)``. It must be the shape of the ``G_kw`` later handed to
    :func:`oracle_sigma2`."""
    recs = []
    for itype, rows in rows_by_type.items():
        for (rx, ry, rz, a1, b1, vr, vi) in rows:
            a, b, v = a1 - 1, b1 - 1, complex(vr, vi)
            if itype in _MIRRORED_ROW_TYPES:
                v = 0.5 * v
            off = (rx, ry, rz) != (0, 0, 0)
            for i in range(shape[0] * shape[1] * shape[2]):
                xi, yi, zi = _coords(i, shape)
                j = _site(xi + rx, yi + ry, zi + rz, shape)
                g = lambda site, o, sp: _g(site, o, sp, norb)
                if itype == "CoulombIntra":
                    recs.append((g(i, a, UP), g(i, a, DN), g(i, a, UP), g(i, a, DN), v, False))
                elif itype == "CoulombInter":
                    for s1, s2 in itertools.product((UP, DN), repeat=2):
                        recs.append((g(i, a, s1), g(j, b, s2), g(i, a, s1), g(j, b, s2), v, off))
                elif itype == "Hund":
                    for s1 in (UP, DN):
                        recs.append((g(i, a, s1), g(j, b, s1), g(i, a, s1), g(j, b, s1), -v, off))
                elif itype == "Ising":
                    for s1, s2 in itertools.product((UP, DN), repeat=2):
                        recs.append((g(i, a, s1), g(j, b, s2), g(i, a, s1), g(j, b, s2),
                                     v if s1 == s2 else -v, off))
                elif itype == "Exchange":
                    for s1 in (UP, DN):
                        # +v c+_{a s} c_{b s} c+_{b -s} c_{a -s} = -v c+_{a s} c+_{b -s} c_{b s} c_{a -s}
                        recs.append((g(i, a, s1), g(j, b, 1 - s1), g(i, a, 1 - s1), g(j, b, s1), -v, off))
                elif itype == "PairHop":
                    # +v c+_{a up} c_{b up} c+_{a dn} c_{b dn} = -v c+_{a up} c+_{a dn} c_{b up} c_{b dn}
                    recs.append((g(i, a, UP), g(i, a, DN), g(j, b, DN), g(j, b, UP), -v, off))
                elif itype == "PairLift":
                    recs.append((g(i, a, UP), g(j, b, UP), g(j, b, DN), g(i, a, DN), -v, off))
                    recs.append((g(j, b, DN), g(i, a, DN), g(i, a, UP), g(j, b, UP), -v, off))
                else:
                    raise ValueError("oracle: unknown interaction type {!r}".format(itype))
    return recs


def _gamma_entries(recs):
    """Antisymmetrised entries as a dict {(p,q,r,s): (value, offsite)} from
    the records (V[p,q,r,s] += x, V[q,p,s,r] += x; Gamma = V - V^x)."""
    V = {}
    for (p, q, r, s, x, off) in recs:
        for key in ((p, q, r, s), (q, p, s, r)):
            val, o = V.get(key, (0.0, off))
            V[key] = (val + x, off or o)
    G = {}
    for (p, q, r, s), (x, off) in V.items():
        for key, sign in (((p, q, r, s), 1.0), ((p, q, s, r), -1.0)):
            val, o = G.get(key, (0.0, off))
            G[key] = (val + sign * x, off or o)
    return {k: v for k, v in G.items() if abs(v[0]) > 0}


def pair_weight(which, off1, crossed1, off2, crossed2):
    """The weight of one record PAIR under the three weightings of spec 2.4.

    ``off*`` says whether the record's vertex is off-site; ``crossed*``
    whether that record's transport leg sits on the site of its own external
    leg (the ``q``-representable, "crossed" placement). The rule, written
    out as a table by
    :meth:`TestWeightRule.test_the_four_documented_classes`:

    ======================  ======  =======  =========
    class                   exact   local    dropped
    ======================  ======  =======  =========
    on / on                 1/2     1/2      0
    on / off (crossed)      1/2     1        -1/2
    off / off both crossed  1/2     1        -1/2
    off / off one direct    1/2     0        1/2
    ======================  ======  =======  =========

    ``dropped`` is ``exact - local`` by construction, which is what makes
    ``local + dropped`` the exact second order in the gates that use it."""
    if which == "exact":
        return 0.5
    if not off1 and not off2:
        w_local = 0.5
    elif (not off1 or crossed1) and (not off2 or crossed2):
        w_local = 1.0
    else:
        w_local = 0.0
    if which == "local":
        return w_local
    if which == "dropped":
        return 0.5 - w_local
    raise ValueError("oracle: unknown weighting {!r}".format(which))


def _site_of(idx, norb):
    return idx // (2 * norb)


def oracle_sigma2(G_kw, beta, recs, norb, which, shape=_SHAPE):
    """Sigma2 on (k, i w) from the skeleton of spec 2.1 with the weight rule
    of 2.4 ('exact': 1/2 every pair; 'local': 1/2 on/on, 1 when every
    off-site entry in the pair is crossed (s on p's site), 0 otherwise;
    'dropped': exact - local).

    ``shape`` is the lattice ``G_kw`` lives on (default: this module's 4x4
    torus) and must be the one ``recs`` was enumerated on."""
    nx, ny, nz = shape
    nvol = nx * ny * nz
    nmat = G_kw.shape[1]
    M = nvol * 2 * norb
    g_rt = _bk.spatial_ifftn(
        _ms.fermion_to_tau(G_kw[0].reshape(nmat, nvol * norb * norb), axis=0).reshape(nmat, nx, ny, nz, norb * norb),
        axes=(1, 2, 3), workers=1).reshape(nmat, nvol, norb, norb)
    rev = np.zeros_like(g_rt)
    for t in range(nmat):
        for i in range(nvol):
            xi, yi, zi = _coords(i, shape)
            rev[t, i] = g_rt[(-t) % nmat, _site(-xi, -yi, -zi, shape)]
    sgn = -np.ones(nmat); sgn[0] = 1.0

    # site-orbital-spin propagators (spin-diagonal, spin-free)
    def Gf(t, P, Q):   # G_{PQ}(tau) forward: site displacement
        i, j = _site_of(P, norb), _site_of(Q, norb)
        if (P // norb) % 2 != (Q // norb) % 2:
            return 0.0
        xi, yi, zi = _coords(i, shape); xj, yj, zj = _coords(j, shape)
        return g_rt[t, _site(xi - xj, yi - yj, zi - zj, shape), P % norb, Q % norb]

    def Gr(t, P, Q):   # G_{PQ}(-tau)
        # rev[d] = G(-d, -tau) by construction, so G_{PQ}(-tau) = G(x_P - x_Q,
        # -tau) is read at the OPPOSITE displacement x_Q - x_P.
        i, j = _site_of(P, norb), _site_of(Q, norb)
        if (P // norb) % 2 != (Q // norb) % 2:
            return 0.0
        xi, yi, zi = _coords(i, shape); xj, yj, zj = _coords(j, shape)
        return rev[t, _site(xj - xi, yj - yi, zj - zi, shape), P % norb, Q % norb]

    Gam = _gamma_entries(recs)
    by_p = {}
    for (p, q, r, s), (x, off) in Gam.items():
        by_p.setdefault(p, []).append((q, r, s, x, off))
    by_pprime = {}
    for (rp, sp, pp, qp), (x, off) in Gam.items():   # Gamma_{r's',p'q'}
        by_pprime.setdefault(pp, []).append((rp, sp, qp, x, off))
    sig = np.zeros((nmat, nvol, norb, norb), complex)
    for p in range(M):
        if (p // norb) % 2 != UP or _site_of(p, norb) != 0:
            continue            # external up spin at site 0 (translation invariance)
        a = p % norb
        for pp in range(M):
            if (pp // norb) % 2 != UP:
                continue
            j2, b = _site_of(pp, norb), pp % norb
            # Sigma_{p p'} is a matrix element between the row site (p, at 0)
            # and the column site (p', at j2); the solver's spatial index is
            # the row-minus-column displacement R = 0 - j2.
            xj2, yj2, zj2 = _coords(j2, shape)
            rsite = _site(-xj2, -yj2, -zj2, shape)
            for (q, r, s, x1, off1) in by_p.get(p, []):
                crossed1 = (_site_of(s, norb) == _site_of(p, norb))
                for (rp, sp, qp, x2, off2) in by_pprime.get(pp, []):
                    crossed2 = (_site_of(sp, norb) == _site_of(pp, norb))
                    w = pair_weight(which, off1, crossed1, off2, crossed2)
                    if w == 0.0:
                        continue
                    for t in range(nmat):
                        val = Gf(t, r, rp) * Gf(t, s, sp) * Gr(t, qp, q)
                        if val != 0.0:
                            sig[t, rsite, a, b] += -w * x1 * x2 * val * sgn[t]
    # sig is Sigma2(R, tau) on the solver's spatial index ; to (k, i w).
    # TWO factors 1/beta:
    # the second-order diagram carries two Matsubara sums in the solver's
    # discrete representation (the internal bubble and the transport), which
    # is how the Phase B enumeration (_sopt_oracle) normalises the same
    # skeleton -- one explicit 1/beta on the (r, tau) product and one on the
    # transform back to (k, i w).
    kt = _bk.spatial_fftn(sig.reshape(nmat, nx, ny, nz, norb * norb), axes=(1, 2, 3), workers=1)
    return (_ms.tau_to_fermion(kt.reshape(nmat, nvol * norb * norb), axis=0)
            .reshape(1, nmat, nvol, norb, norb) / (beta * beta))


def oracle_w2(G_kw, beta, recs, norb, shape=_SHAPE):
    """The oracle's LOCAL second-order KERNEL ``W2(q, i nu)``, laid out in the
    solver's pair flattening, derived from the same records and the same
    weight rule as :func:`oracle_sigma2` -- independently of production.

    :func:`oracle_sigma2` compares a SELF-ENERGY; the transport that turns a
    kernel into a self-energy at one fixed ``G`` is a linear map from
    ``nmat nvol nd^2`` components onto ``nmat nvol norb^2``, so it has a large
    nullspace and a self-energy agreement cannot see a kernel error inside it
    (the mutation test below exhibits one). This function therefore builds the
    kernel itself.

    Derivation. The solver's transport is a per-``(R, tau)`` product,
    ``Sigma_{ab}(R, tau) = sum_{cd} V_{(c a),(d b)}(R, tau) G_{cd}(R, tau)``,
    followed by one explicit ``1/beta``
    (:meth:`hwave.solver.flex.FLEX._calc_self_energy_general`). Reading the
    skeleton of spec 2.1 with the ``s``-leg as the transport leg, the external
    index pair ``(c a)`` sits on the row site and ``(d b)`` on the column site,
    so the kernel is the skeleton with that leg amputated:

        W2_{(c a),(d b)}(R, tau) = -(1/beta) sgn(tau) sum_{sigma_s}
            sum w Gamma_{p q, r s} Gamma_{r' s', p' q'} G_{r r'}(tau) G_{q' q}(-tau)

    with ``p = (0, a, up)``, ``s = (0, c, sigma_s)``, ``p' = (j, b, up)``,
    ``s' = (j, d, sigma_s)``, ``R = -j`` and the remaining ``1/beta`` left to the
    transport -- the two factors :func:`oracle_sigma2` applies together. Because
    the transport leg of each vertex is pinned to the site of that vertex's
    external leg, every record pair reachable here is ``q``-representable
    ("crossed" in the language of 2.4), so the weight rule of spec 2.4 reduces
    to ``1/2`` for an on-site/on-site pair and ``1`` as soon as one record is
    off-site -- the surviving rows of that table, written out here.

    Returns ``(nmat, nvol, nd, nd)``, comparable element by element with
    ``hwave.solver.second_order.dense_w2``.
    """
    nx, ny, nz = shape
    nvol = nx * ny * nz
    nmat = G_kw.shape[1]
    g_rt = _bk.spatial_ifftn(
        _ms.fermion_to_tau(G_kw[0].reshape(nmat, nvol * norb * norb), axis=0).reshape(nmat, nx, ny, nz, norb * norb),
        axes=(1, 2, 3), workers=1).reshape(nmat, nvol, norb, norb)
    rev = np.zeros_like(g_rt)
    for t in range(nmat):
        for i in range(nvol):
            xi, yi, zi = _coords(i, shape)
            rev[t, i] = g_rt[(-t) % nmat, _site(-xi, -yi, -zi, shape)]
    sgn = -np.ones(nmat); sgn[0] = 1.0

    def Gf(P, Q):      # G_{PQ}(tau) for every tau at once, or None if spin-forbidden
        if (P // norb) % 2 != (Q // norb) % 2:
            return None
        i, j = _site_of(P, norb), _site_of(Q, norb)
        xi, yi, zi = _coords(i, shape); xj, yj, zj = _coords(j, shape)
        return g_rt[:, _site(xi - xj, yi - yj, zi - zj, shape), P % norb, Q % norb]

    def Gr(P, Q):      # G_{PQ}(-tau); see oracle_sigma2 for the reversal direction
        if (P // norb) % 2 != (Q // norb) % 2:
            return None
        i, j = _site_of(P, norb), _site_of(Q, norb)
        xi, yi, zi = _coords(i, shape); xj, yj, zj = _coords(j, shape)
        return rev[:, _site(xj - xi, yj - yi, zj - zi, shape), P % norb, Q % norb]

    Gam = _gamma_entries(recs)
    v1 = {}
    for (p, q, r, s), (x, off) in Gam.items():             # Gamma_{pq,rs}, keyed by (p, s)
        v1.setdefault((p, s), []).append((q, r, x, off))
    v2 = {}
    for (rp, sp, pp, qp), (x, off) in Gam.items():         # Gamma_{r's',p'q'}, keyed by (p', s')
        v2.setdefault((pp, sp), []).append((rp, qp, x, off))

    nd = norb * norb
    W = np.zeros((nmat, nvol, nd, nd), complex)
    for j2 in range(nvol):
        xj2, yj2, zj2 = _coords(j2, shape)
        rsite = _site(-xj2, -yj2, -zj2, shape)
        for a, b, c, d in itertools.product(range(norb), repeat=4):
            p = _g(0, a, UP, norb)
            pp = _g(j2, b, UP, norb)
            acc = np.zeros(nmat, complex)
            for ss in (UP, DN):
                s = _g(0, c, ss, norb)
                sp = _g(j2, d, ss, norb)
                for (q, r, x1, off1) in v1.get((p, s), []):
                    for (rp, qp, x2, off2) in v2.get((pp, sp), []):
                        gf = Gf(r, rp)
                        if gf is None:
                            continue
                        gr = Gr(qp, q)
                        if gr is None:
                            continue
                        w = 0.5 if (not off1 and not off2) else 1.0
                        acc += w * x1 * x2 * gf * gr
            W[:, rsite, c * norb + a, d * norb + b] += -acc * sgn / beta
    Wq = _bk.spatial_fftn(W.reshape(nmat, nx, ny, nz, nd * nd), axes=(1, 2, 3), workers=1)
    return _ms.tau_to_boson(Wq.reshape(nmat, nvol * nd * nd), axis=0).reshape(nmat, nvol, nd, nd)


def _bare_green(s, beta):
    s._calc_epsilon_k({})
    nmat, nvol, norb = s.nmat, s.lattice.nvol, s.norb
    return s._calc_dressed_green(beta, 0.1, np.zeros((1, nmat, nvol, norb, norb), complex))


def _general_path_second_order_sigma(s, G, beta):
    """The SECOND ORDER of today's general path on the same G, from
    production pieces only: the Takimoto-Hotta-Ueda assembly
    ``_calc_veff_general`` with the RPA channels replaced by the bare bubble
    (chi_s = chi_c = chibar), i.e. 3/2 Us chibar Us + 1/2 Uc chibar Uc
    - 1/4 (Us + Uc) chibar (Us + Uc) = U^2 chibar for the Hubbard tensor
    (spec 2.3), pushed through the production transport.

    The full ``_flex_compute_veff_general`` cannot serve as the reference
    here: it resums the RPA ladders, so its self-energy also carries third
    and higher orders in U (measured: 2x the second order at U = 1 on this
    fixture), which an exact second-order oracle must not reproduce."""
    chi0q_raw = s._calc_chi0q(G, np.zeros_like(G), beta)[0]
    chi0q, Us, Uc = s._inflate_chi0q_and_ham_general(chi0q_raw, s.ham_info.ham_inter_q)
    v_eff2 = s._calc_veff_general(chi0q, chi0q, chi0q, Us, Uc)
    return s._calc_self_energy_general(G, v_eff2, beta)


_ONSITE = {
    "U": {"CoulombIntra": [(0, 0, 0, 1, 1, 1.0, 0.0), (0, 0, 0, 2, 2, 1.0, 0.0)]},
    "Up": {"CoulombInter": [(0, 0, 0, 1, 2, 1.0, 0.0), (0, 0, 0, 2, 1, 1.0, 0.0)]},
    "J": {"Hund": [(0, 0, 0, 1, 2, 1.0, 0.0), (0, 0, 0, 2, 1, 1.0, 0.0)]},
    "I": {"Ising": [(0, 0, 0, 1, 2, 1.0, 0.0), (0, 0, 0, 2, 1, 1.0, 0.0)]},
    "X": {"Exchange": [(0, 0, 0, 1, 2, 1.0, 0.0), (0, 0, 0, 2, 1, 1.0, 0.0)]},
    "PH": {"PairHop": [(0, 0, 0, 1, 2, 1.0, 0.0), (0, 0, 0, 2, 1, 1.0, 0.0)]},
    "PL": {"PairLift": [(0, 0, 0, 1, 2, 1.0, 0.0), (0, 0, 0, 2, 1, 1.0, 0.0)]},
}
_OFFSITE = {
    "V": {"CoulombInter": [(1, 0, 0, 1, 2, 1.0, 0.0), (-1, 0, 0, 2, 1, 1.0, 0.0),
                           (0, 1, 0, 1, 1, 0.6, 0.0), (0, -1, 0, 1, 1, 0.6, 0.0)]},   # asymmetric bond
    "JV": {"Hund": [(1, 0, 0, 1, 1, 1.0, 0.0), (-1, 0, 0, 1, 1, 1.0, 0.0)]},
    "IV": {"Ising": [(1, 0, 0, 1, 1, 1.0, 0.0), (-1, 0, 0, 1, 1, 1.0, 0.0)]},
}


def _scaled(rows_by_type, x):
    return {t: [r[:5] + (x * r[5], x * r[6]) for r in rows] for t, rows in rows_by_type.items()}


def _merge(*dicts):
    out = {}
    for d in dicts:
        for t, rows in d.items():
            out.setdefault(t, []).extend(rows)
    return out


class TestOracle(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.s, cls.split = _split_for({"CoulombIntra": [(0, 0, 0, 1, 1, 1.0, 0.0), (0, 0, 0, 2, 2, 1.0, 0.0)]})
        cls.G = _bare_green(cls.s, _BETA)

    def test_g1_hubbard_reproduces_the_general_path(self):
        recs = oracle_records(_ONSITE["U"], 2)
        sig = oracle_sigma2(self.G, _BETA, recs, 2, "exact")
        ref = _general_path_second_order_sigma(self.s, self.G, _BETA)
        self.assertGreater(np.abs(ref).max(), 1e-6)                 # anti-vacuity
        np.testing.assert_allclose(sig, ref, rtol=1e-12, atol=1e-12 * np.abs(ref).max())

    def test_g1_density_interaction_reproduces_the_phase_b_enumeration(self):
        from tests.test_flex_bond_sopt import _sopt_oracle, _interaction_matrices
        rows1 = [(1, 0, 0, 1, 1, 0.3, 0.0), (-1, 0, 0, 1, 1, 0.3, 0.0)]
        s1, _ = _split_for({"CoulombInter": rows1}, norb=2)
        G1 = _bare_green(s1, _BETA)
        recs = oracle_records({"CoulombInter": rows1}, 2)
        exact = oracle_sigma2(G1, _BETA, recs, 2, "exact")
        Wuu, Wud = _interaction_matrices(2, 0.0, 0.0, 0.0, [r[:5] + (r[5], r[6]) for r in rows1])
        d, x = _sopt_oracle(G1, _BETA, 2, Wuu, Wud)
        self.assertGreater(np.abs(d).max(), 1e-6)                   # anti-vacuity
        self.assertGreater(np.abs(x).max(), 1e-3 * np.abs(d).max())  # the exchange class is present
        np.testing.assert_allclose(exact[0], d + x, rtol=1e-12, atol=1e-12 * np.abs(d).max())
        local = oracle_sigma2(G1, _BETA, recs, 2, "local")
        np.testing.assert_allclose(local[0], d, rtol=1e-12, atol=1e-12 * np.abs(d).max())

    def test_g2a_every_pair_direct_kernel_comparison(self):
        """production W2 (accumulate_batch on the solver's chibar) vs the
        oracle's local kernel, every pair of coupling types, asymmetric bond."""
        from hwave.solver.second_order import build_factors, dense_w2
        names = list(_ONSITE) + list(_OFFSITE)
        table = dict(_ONSITE, **_OFFSITE)
        expected_zero = set()          # filled from symmetry: pairs whose oracle kernel vanishes
        for x, y in itertools.combinations_with_replacement(names, 2):
            with self.subTest(pair=(x, y)):
                rows = _merge(_scaled(table[x], 0.7), _scaled(table[y], 0.4)) if x != y else _scaled(table[x], 0.7)
                s, split = _split_for(rows)
                G = _bare_green(s, _BETA)
                f = build_factors(split, s.lattice, 2)
                chi0q = s._calc_chi0q(G, np.zeros_like(G), _BETA)[0].reshape(_NMAT, _NVOL, 4, 4)
                W2 = dense_w2(chi0q, f)
                sig_prod = s._calc_self_energy_general(G, W2, _BETA)
                sig_orc = oracle_sigma2(G, _BETA, oracle_records(rows, 2), 2, "local")
                scale = max(np.abs(sig_orc).max(), 1e-300)
                if np.abs(sig_orc).max() < 1e-13:
                    expected_zero.add((x, y))
                    self.assertLess(np.abs(sig_prod).max(), 1e-13)
                else:
                    self.assertGreater(np.abs(sig_orc).max(), 1e-6)         # anti-vacuity
                    np.testing.assert_allclose(sig_prod, sig_orc, rtol=1e-10, atol=1e-12 * scale)
        self.assertEqual(expected_zero, set())      # record the symmetry-zero pairs here if any appear

    def test_g2a_every_pair_direct_w2_kernel_comparison(self):
        """production W2 (``dense_w2``) vs the oracle's local KERNEL
        (:func:`oracle_w2`), element by element, every pair of coupling types
        on the asymmetric bond.

        The companion test above compares self-energies; the transport that
        produces them has a nullspace (see
        :meth:`test_direct_comparison_catches_a_transport_nullspace_mutation`),
        so the kernel is compared here directly."""
        from hwave.solver.second_order import build_factors, dense_w2
        names = list(_ONSITE) + list(_OFFSITE)
        table = dict(_ONSITE, **_OFFSITE)
        expected_zero = set()          # filled from symmetry: pairs whose oracle kernel vanishes
        for x, y in itertools.combinations_with_replacement(names, 2):
            with self.subTest(pair=(x, y)):
                rows = _merge(_scaled(table[x], 0.7), _scaled(table[y], 0.4)) if x != y else _scaled(table[x], 0.7)
                s, split = _split_for(rows)
                G = _bare_green(s, _BETA)
                f = build_factors(split, s.lattice, 2)
                chi0q = s._calc_chi0q(G, np.zeros_like(G), _BETA)[0].reshape(_NMAT, _NVOL, 4, 4)
                W2 = dense_w2(chi0q, f)
                W2_orc = oracle_w2(G, _BETA, oracle_records(rows, 2), 2)
                ref = np.abs(W2_orc).max()
                if ref < 1e-13:
                    expected_zero.add((x, y))
                    self.assertLess(np.abs(W2).max(), 1e-13)
                else:
                    self.assertGreater(ref, 1e-6)                   # anti-vacuity, on W2 itself
                    self.assertLess(np.abs(W2 - W2_orc).max(), 1e-11 * ref)
        self.assertEqual(expected_zero, set())      # record the symmetry-zero pairs here if any appear

    def test_direct_comparison_catches_a_transport_nullspace_mutation(self):
        """The transport at a fixed ``G`` maps ``nmat nvol nd^2`` kernel
        components onto ``nmat nvol norb^2`` self-energy components, so it has
        a nullspace: a kernel error inside it is invisible to a self-energy
        comparison and visible only to the direct one.

        The nullspace direction is found NUMERICALLY (a singular-value
        decomposition of one ``(tau, R)`` block of the transport, which is a
        per-``(tau, R)`` product and therefore block diagonal there), injected
        at that single block, and carried to ``(q, i nu)`` by the inverse of
        the transport's own transforms."""
        from hwave.solver.second_order import build_factors, dense_w2
        rows = _merge(_scaled(_ONSITE["U"], 0.7), _scaled(_OFFSITE["V"], 0.4))
        s, split = _split_for(rows)
        G = _bare_green(s, _BETA)
        f = build_factors(split, s.lattice, 2)
        chi0q = s._calc_chi0q(G, np.zeros_like(G), _BETA)[0].reshape(_NMAT, _NVOL, 4, 4)
        W2 = dense_w2(chi0q, f)
        W2_orc = oracle_w2(G, _BETA, oracle_records(rows, 2), 2)
        self.assertLess(np.abs(W2 - W2_orc).max(), 1e-11 * np.abs(W2_orc).max())

        norb, nd = 2, 4
        nx, ny, nz = _SHAPE
        g_rt = _bk.spatial_ifftn(
            _ms.fermion_to_tau(G[0].reshape(_NMAT, _NVOL * norb * norb), axis=0)
            .reshape(_NMAT, nx, ny, nz, norb * norb),
            axes=(1, 2, 3), workers=1).reshape(_NMAT, _NVOL, norb, norb)
        # Sigma[a, b] = sum_{c, d} W[(c a), (d b)] G[c, d] at the (tau, R) block (0, 0)
        M = np.zeros((norb * norb, nd * nd), complex)
        for a, b, c, d in itertools.product(range(norb), repeat=4):
            M[a * norb + b, (c * norb + a) * nd + (d * norb + b)] = g_rt[0, 0, c, d]
        _, sv, vh = np.linalg.svd(M)
        rank = int(np.sum(sv > 1e-12 * sv.max()))
        self.assertEqual(rank, norb * norb)             # the block map is onto
        self.assertGreater(vh.shape[0] - rank, 0)       # ... and has a nullspace
        dW_rt = np.zeros((_NMAT, _NVOL, nd, nd), complex)
        dW_rt[0, 0] = vh[rank].conj().reshape(nd, nd) * np.abs(W2).max()
        dW = _ms.tau_to_boson(
            _bk.spatial_fftn(dW_rt.reshape(_NMAT, nx, ny, nz, nd * nd),
                             axes=(1, 2, 3), workers=1).reshape(_NMAT, _NVOL * nd * nd),
            axis=0).reshape(_NMAT, _NVOL, nd, nd)

        sig = s._calc_self_energy_general(G, W2, _BETA)
        sig_mut = s._calc_self_energy_general(G, W2 + dW, _BETA)
        self.assertGreater(np.abs(sig).max(), 1e-6)
        # invisible to the self-energy comparison ...
        self.assertLess(np.abs(sig_mut - sig).max(), 1e-13 * np.abs(sig).max())
        # ... and caught by the direct one
        self.assertGreater(np.abs(W2 + dW - W2_orc).max(), 1e-2 * np.abs(W2_orc).max())

    def test_dropped_class_is_load_bearing_for_v_only(self):
        s, split = _split_for(_OFFSITE["V"])
        G = _bare_green(s, _BETA)
        recs = oracle_records(_OFFSITE["V"], 2)
        dropped = oracle_sigma2(G, _BETA, recs, 2, "dropped")
        local = oracle_sigma2(G, _BETA, recs, 2, "local")
        self.assertGreater(np.abs(dropped).max(), 1e-3 * np.abs(local).max())
        for name in ("J", "X"):
            rows = _merge(_ONSITE[name], _OFFSITE["V"])
            s, split = _split_for(rows)
            G = _bare_green(s, _BETA)
            recs = oracle_records(rows, 2)
            cross_dropped = (oracle_sigma2(G, _BETA, recs, 2, "dropped")
                             - oracle_sigma2(G, _BETA, oracle_records(_OFFSITE["V"], 2), 2, "dropped"))
            self.assertLess(np.abs(cross_dropped).max(), 1e-12)     # mixed on/off: nothing dropped



class TestWeightRule(unittest.TestCase):
    """The locality weight rule of spec 2.4, stated as a table and checked
    on SYNTHETIC records rather than only through the self-energies it
    produces."""

    #: The documented per-record-pair LOCAL weight, written here from the
    #: spec's prose rather than from :func:`pair_weight`: a pair of on-site
    #: records carries 1/2 (both of its two antisymmetrisation copies are
    #: representable, and the skeleton factor halves them); a pair
    #: containing an off-site record carries 1 when EVERY off-site record
    #: in it is in its crossed placement (only one copy is representable,
    #: so it carries the whole diagram) and 0 otherwise.
    _LOCAL = {
        ("on", "on"): 0.5,
        ("on", "off crossed"): 1.0,
        ("on", "off direct"): 0.0,
        ("off crossed", "off crossed"): 1.0,
        ("off crossed", "off direct"): 0.0,
        ("off direct", "off direct"): 0.0,
    }
    _FLAGS = {"on": (False, True), "off crossed": (True, True), "off direct": (True, False)}

    def test_the_four_documented_classes(self):
        for (k1, k2), expect in self._LOCAL.items():
            for (a, b) in ((k1, k2), (k2, k1)):
                off1, cr1 = self._FLAGS[a]
                off2, cr2 = self._FLAGS[b]
                with self.subTest(classes=(a, b)):
                    self.assertEqual(pair_weight("local", off1, cr1, off2, cr2), expect)
                    self.assertEqual(pair_weight("exact", off1, cr1, off2, cr2), 0.5)
                    # dropped is exact - local by construction
                    self.assertAlmostEqual(
                        pair_weight("dropped", off1, cr1, off2, cr2), 0.5 - expect, places=15)
        with self.assertRaises(ValueError):
            pair_weight("nonsense", False, True, False, True)

    def test_synthetic_records_realise_every_class(self):
        """A two-site chain with an on-site ``U`` and an off-site ``V`` bond
        produces records of all three placement classes at one external
        index, and each pair of them is weighted as the table says."""
        shape = (2, 1, 1)
        rows = {"CoulombIntra": [(0, 0, 0, 1, 1, 1.0, 0.0)],
                "CoulombInter": [(1, 0, 0, 1, 1, 0.5, 0.0), (-1, 0, 0, 1, 1, 0.5, 0.0)]}
        Gam = _gamma_entries(oracle_records(rows, 1, shape))
        p = _g(0, 0, UP, 1)
        seen = {}
        for (pk, q, r, s), (x, off) in Gam.items():
            if pk != p:
                continue
            crossed = (_site_of(s, 1) == _site_of(pk, 1))
            key = "on" if not off else ("off crossed" if crossed else "off direct")
            seen.setdefault(key, []).append((off, crossed))
        self.assertEqual(set(seen), set(self._FLAGS),
                         "the synthetic fixture does not realise every placement class")
        for k1, entries1 in seen.items():
            for k2, entries2 in seen.items():
                key = (k1, k2) if (k1, k2) in self._LOCAL else (k2, k1)
                for (off1, cr1) in entries1[:1]:
                    for (off2, cr2) in entries2[:1]:
                        with self.subTest(classes=(k1, k2)):
                            self.assertEqual(pair_weight("local", off1, cr1, off2, cr2),
                                             self._LOCAL[key])


if __name__ == "__main__":
    unittest.main()
