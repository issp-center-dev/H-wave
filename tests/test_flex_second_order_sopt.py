"""G2 (b)-(d) of spec 2026-09-08: END-TO-END second-order coefficients of the
production ``Sigma_fluct`` (frozen bare G at a fixed mu, quadratic fit on a
3 x 3 coupling grid with three-level Richardson extrapolation at h = 2e-3)
equal the independent real-space oracle's for a covering set of pure
interaction types and mixed pairs; the dropped (off/off uncrossed) class is
load-bearing for an off-site V; and the bond gate reproduces the FULL
(exact) oracle for an off-site CoulombInter bond, orbital-diagonal and
inter-orbital alike.

What is asserted and what is recorded: the whole covering set of (b) and
the dropped-class ratio are ASSERTED, and so is the bond gate against the
exact oracle on BOTH off-site CoulombInter bonds -- the orbital-diagonal
one and the orbital-off-diagonal one, the latter since the pair
permutation of spec 2026-09-16 R3 closed issue #192 (it used to miss by
1.6e-3, and is adjudicated against exact diagonalisation by the chain
gate, ``tests/test_flex_second_order_ed_chain.py``). The gate's deviation
on off-site Hund/Ising is RECORDED by a print (spec G2 (c)), because that
class is adjudicated by the chain gate rather than here; the two records
are held to a ceiling (:data:`_ROUNDOFF_CEIL`) so they cannot drift
unnoticed.

The EXTRACTION is checked too, not only its outcome
---------------------------------------------------
Production and oracle are fitted on the same grids and alias the same way,
so their agreement survives an extraction that has stopped measuring a
second-order coefficient at all. :meth:`TestG2Heavy._check_extraction`
therefore holds the fit itself to account on every entry: the relative
least-squares residual stays under :data:`_RESIDUAL_BOUND` and improves as
the grid is refined, every fitted quantity is finite, and the three-level
Richardson estimate agrees with the two-level one from the finer pair
(:data:`_LADDER_BOUND`) -- i.e. the ladder has settled rather than merely
been applied.

Relation to the neighbouring gates
----------------------------------
``tests/test_second_order_oracle.py`` (G2 (a)) compares the production
KERNEL (``dense_w2`` on the solver's chibar) with the oracle for all 55
type pairs at one coupling. This module instead drives the production
SOLVER end to end -- ``_flex_compute_veff_general`` (which resums the RPA
ladders, so its self-energy carries third and higher orders too) and the
bond gate's ``dress_and_build_w`` / ``calc_self_energy_bond`` -- and
isolates the second order by fitting. It is therefore the gate that would
catch a second-order kernel that is correct in isolation but mis-wired,
mis-weighted or double-counted inside the production assembly.

Why a frozen G: the coefficient of U^2 is only defined once the propagator
is held fixed. Every map here is built from ``_bare_green`` at mu = 0.1 --
no SCF, no density update -- on BOTH sides of the comparison, so the two
sides differ only in how the second order is assembled.
"""
import contextlib
import logging
import os
import shutil
import tempfile
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k
from tests.heavy_tests import heavy
from tests.test_second_order_oracle import (_ONSITE, _OFFSITE, _merge, _scaled, _bare_green,
                                            _BETA, oracle_records, oracle_sigma2)
from tests.test_flex_bond_sopt import _fit, _richardson, _rel

_IN2 = "tests/rpa/input_2orb"
_SHAPE = (4, 4, 1)
_NMAT = 8
_NORB = 2
_MU = 0.1

#: Coupling scale of the coarsest extraction grid; the Richardson ladder
#: adds the grids at h/2 and h/4 (the same ladder as
#: ``tests/test_flex_bond_sopt.py``, whose ``_fit``/``_richardson`` this
#: module reuses).
_H = 2.0e-3

#: The covering set. Pure types exercise every diagonal of the second-order
#: kernel; the pairs exercise every mixed class that the locality split can
#: produce -- on/on (the first seven), and on/off with the asymmetric
#: off-site bond ``_OFFSITE["V"]`` (the last five).
_PAIRS = [("U", "Up"), ("U", "J"), ("U", "I"), ("U", "X"), ("U", "PH"), ("U", "PL"), ("Up", "J"),
          ("U", "V"), ("J", "V"), ("I", "V"), ("X", "V"), ("PH", "V")]
_PURE = ["U", "Up", "J", "I", "X", "PH", "PL", "V", "JV", "IV"]

#: Anti-vacuity floors on the ORACLE coefficient that each assertion
#: compares against: a relative comparison against a numerically zero
#: reference would pass for any production value at all.
_FLOOR_PURE = 1e-6
_FLOOR_PAIR = 1e-8

#: Covering-set entries whose second-order coefficient is STRUCTURALLY zero
#: in both production and the oracle, with the symmetry that makes it so.
#: An entry listed here is checked for a TWO-SIDED zero instead of a
#: relative agreement (a relative check against a zero reference is
#: vacuous, and its "relative deviation" is pure round-off).
#:
#: Both entries vanish by SPIN algebra alone -- no property of this
#: fixture's geometry or band structure enters, so they are zero for every
#: input, not merely small here. ``CoulombIntra``'s only monomial is
#: ``n_{a up} n_{a dn}``, so at second order the U vertex enters the
#: skeleton with one up and one down leg on EACH side; the partner vertex
#: must therefore supply an (up, dn) pair on each side too. On-site
#: ``Hund`` has only all-same-spin monomials ``n_{a s} n_{b s}``, and
#: on-site ``PairLift`` only spin-paired ones (``c+_{a up} c+_{b up} c_{b
#: dn} c_{a dn}``, both legs up on one side and both down on the other);
#: neither has an (up, dn) entry on both sides, so the U x J and U x PL
#: cross terms are identically zero. Measured: oracle 2.1e-17 at unit
#: couplings (the oracle is exactly quadratic, so that is the exact cross
#: term), production 5.8e-13 after the Richardson ladder -- against pure
#: coefficients of 3.6e-2 and the smallest NONZERO cross term of the
#: covering set at 3.3e-4.
_KNOWN_ZERO = frozenset({("U", "J"), ("U", "PL")})

#: An ORBITAL-DIAGONAL off-site CoulombInter bond (``v_ab(R) = 0`` for
#: ``a != b``), on two different shells so the two orbitals are not
#: symmetry-equivalent. This is the fixture G2 (c) asserts on: the bond
#: gate reproduces the exact oracle here to 3e-11, and the class it has to
#: recover to do so -- the off/off uncrossed one the local weighting drops
#: -- is 2.1e-2 of the coefficient, so the assertion distinguishes the
#: exact from the local weighting by four orders more than the tolerance.
_V_DIAG = {"CoulombInter": [(1, 0, 0, 1, 1, 1.0, 0.0), (-1, 0, 0, 1, 1, 1.0, 0.0),
                            (0, 1, 0, 2, 2, 0.6, 0.0), (0, -1, 0, 2, 2, 0.6, 0.0)]}

#: Bound on the quadratic fit's RELATIVE least-squares residual on any of
#: the three Richardson grids. Measured: 3.3e-5 at the coarsest grid, and
#: exactly halving at each refinement -- the residual is the cubic term of
#: the map, which is what the ladder then extrapolates away. A residual
#: above this bound means the 3 x 3 quadratic fit no longer describes the
#: sampled maps, and the coefficients it returns are not the ones the gate
#: thinks it is comparing.
_RESIDUAL_BOUND = 1.0e-4

#: Bound on the relative difference between the THREE-level Richardson
#: estimate this module uses and the TWO-level one from the two finer grids
#: alone. Measured: 2.5e-6 at worst over the entries checked. It is much
#: larger than the 1e-8 at which production and oracle agree, because the
#: residual aliasing is the SAME on both sides and cancels in that
#: comparison; what this bound checks is that the ladder itself has settled.
_LADDER_BOUND = 1.0e-4

#: Below this the fit residual is round-off, not the map's cubic term, and
#: neither its size nor its scaling says anything about the extraction.
#: ``PairLift`` alone (whose first order vanishes identically, so the fitted
#: map is tiny) lands here, at 5e-11 falling to 2e-10 -- pure noise, and
#: three orders below every other entry's 1e-6 ... 3e-5.
_RESIDUAL_FLOOR = 1.0e-9

#: RECORDED gate-on deviations. ``JV``/``IV`` are off-site Hund/Ising,
#: where the gate reproduces the exact oracle and the recorded number is
#: round-off. (``V``, the inter-orbital off-site CoulombInter bond of issue
#: #192, is no longer recorded: since the pair permutation of spec
#: 2026-09-16 R3 the gate reproduces the exact oracle there too, so it is
#: ASSERTED alongside the orbital-diagonal ``Vd``.)
_PINNED_GATE = {"JV": 3.5e-11, "IV": 3.0e-11}

#: Ceiling for the entries of :data:`_PINNED_GATE`. They sit at the
#: round-off floor of a 1e-1-sized coefficient, so only their CEILING is
#: meaningful (a lower bound there would be asserting the exact bit pattern
#: of a BLAS reduction).
_ROUNDOFF_CEIL = 1.0e-9

_TABLE = dict(_ONSITE, **_OFFSITE)
_TABLE["Vd"] = _V_DIAG      # not in the covering set; used by G2 (c) only


def _write_wan(path, name, norb, rows):
    rvecs = sorted({tuple(r[:3]) for r in rows})
    with open(path, "w") as fw:
        fw.write("{} in wannier90-like format for uhfk\n{}\n{}\n".format(name, norb, len(rvecs)))
        fw.write(" ".join("1" for _ in rvecs) + "\n")
        for r in rows:
            fw.write("{:4d} {:4d} {:4d} {:4d} {:4d} {: .15e} {: .15e}\n".format(*r))


@contextlib.contextmanager
def _quiet():
    """Silence the solver's per-map INFO/WARNING narration.

    Every general-scheme map with an off-site declaration logs one WARNING
    naming what the q-only vertex does and does not carry, and the bond
    gate's preflight logs its memory table at INFO. Both are correct and
    both are pinned by their own tests; here they would emit some two
    thousand lines over the grids of one test method and bury the RECORDED
    lines this module does mean to print."""
    lg = logging.getLogger("hwave")
    old = lg.level
    lg.setLevel(logging.ERROR)
    try:
        yield
    finally:
        lg.setLevel(old)


def _solver(d, rows_by_type, gate):
    """A FLEX general solver on ``input_2orb``'s geometry/transfer with the
    given interaction rows written into the temporary directory ``d``.

    ``flex_second_order = "local"`` is passed as a PARAMETER rather than
    stamped onto the instance afterwards, so the solver compiles its own
    ``_second_order_factors`` through the production code path; likewise the
    gate is switched on through its own two parameters, so the Phase B
    preflight and vertex preparation run exactly as in a real solve."""
    import hwave.solver.flex as flex_mod
    for f in ("geom.dat", "transfer.dat"):
        shutil.copy(os.path.join(_IN2, f), d)
    idict = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
    for t, rows in rows_by_type.items():
        _write_wan(os.path.join(d, t.lower() + ".dat"), t, _NORB, rows)
        idict[t] = t.lower() + ".dat"
    r = read_input_k.QLMSkInput({"path_to_input": d, "interaction": idict})
    par = {"T": 1.0 / _BETA, "filling": 0.5, "CellShape": list(_SHAPE), "SubShape": [1, 1, 1],
           "Nmat": _NMAT, "IterationMax": 1, "Mix": 1.0, "EPS": 1e-12,
           "flex_second_order": "local"}
    if gate:
        par["flex_hartree_fock"] = True
        par["longitudinal_bond_channels"] = True
    info = {"mode": "FLEX", "param": par, "enable_spin_orbital": False, "calc_scheme": "general"}
    return flex_mod.FLEX(r.get_param("ham"), {}, info)


def _one_map_sigma(rows, gate):
    """Production ``Sigma_fluct`` of ONE map from the frozen bare G.

    Returns ``(sigma, G)``. Both branches feed the bubble an explicitly
    ZERO tail (``np.zeros_like(G)``), which the dense kernel treats exactly
    as ``None``: the frozen G here is the full Green function, not a
    tail-deflated one."""
    from hwave.solver import flex_bond
    with tempfile.TemporaryDirectory() as d, _quiet():
        s = _solver(d, rows, gate)
        if not gate:
            G = _bare_green(s, _BETA)
            chi0q_raw = s._calc_chi0q(G, np.zeros_like(G), _BETA)[0]
            _, v_eff, _, _ = s._flex_compute_veff_general(chi0q_raw, s.ham_info.ham_inter_q)
            return s._calc_self_energy_general(G, v_eff, _BETA), G
        green_info = {}
        s._phase_b_reset(green_info)
        s._phase_b_preflight(green_info)
        G = _bare_green(s, _BETA)
        nmat, nvol, norb = s.nmat, s.lattice.nvol, s.norb
        nd = norb * norb
        B = s._bond_view.n_channels
        with flex_bond.BondBlockStore(nmat, nvol, B * nd, nd, ("chibar", "W")) as store:
            s._phase_b_prepare_vertices()
            flex_bond.assemble_bubble(store, G, np.zeros_like(G), _BETA, s._bond_view, _SHAPE, 1)
            flex_bond.dress_and_build_w(store, s._bond_S, s._bond_C, S_on=s._bond_S_on,
                                        C_on=s._bond_C_on, nb=nmat, output_full=False, nmat=nmat,
                                        nvol=nvol, nd=nd, spatial_shape=_SHAPE,
                                        factors=s._second_order_factors, second_order="local")
            return flex_bond.calc_self_energy_bond(store, G, _BETA, s._bond_view, _SHAPE, norb, 1), G


def _grid_points(h):
    return [(x, y) for x in (0.0, h / 2, h) for y in (0.0, h / 2, h)]


def _rows_at(name_x, name_y, x, y):
    if name_x == name_y:
        # one type on a single ray: Sigma is a function of (x + y), so a
        # pure coefficient c appears in the 3 x 3 fit as c20 = c02 = c and
        # c11 = 2 c (the cross term of (x + y)^2). Only c20 is read back.
        return _scaled(_TABLE[name_x], x + y)
    return _merge(_scaled(_TABLE[name_x], x), _scaled(_TABLE[name_y], y))


def _coefficients(name_x, name_y, gate, which):
    """Richardson-extrapolated (c20, c11, c02) of production and oracle.

    Returns ``(out, diag)`` where ``out[k] = (production, oracle)`` and
    ``diag`` carries what the extraction itself has to be held to:

    ``diag["residuals"]``
        the largest RELATIVE least-squares residual (rss over the nine grid
        points divided by the norm of the fitted values) seen on either
        side, ONE PER GRID, coarsest first -- the number that says whether
        the quadratic fit actually describes the sampled maps, and whether
        it improves as the grid is refined.
    ``diag["two_level"][k]``
        the TWO-level extrapolation ``2 c(h/4) - c(h/2)`` from the two
        finer grids alone, as ``(production, oracle)``. Comparing it with
        the three-level ``out[k]`` is what says the ladder has settled
        rather than merely been applied.
    """
    cs_prod, cs_orc, res = [], [], []
    for f in (1.0, 0.5, 0.25):
        pts = _grid_points(_H * f)
        vals_p, vals_o = [], []
        for (x, y) in pts:
            rows = _rows_at(name_x, name_y, x, y)
            sig, G = _one_map_sigma(rows, gate)
            vals_p.append(sig)
            vals_o.append(oracle_sigma2(G, _BETA, oracle_records(rows, _NORB), _NORB, which))
        cp, rp = _fit(pts, vals_p)
        co, ro = _fit(pts, vals_o)
        r = 0.0
        for vals, rss in ((vals_p, rp), (vals_o, ro)):
            nrm = np.linalg.norm(np.asarray(vals).ravel())
            if nrm > 0.0:
                r = max(r, rss / nrm)
        res.append(r)
        cs_prod.append(cp)
        cs_orc.append(co)
    out = {k: (_richardson([c[k] for c in cs_prod]), _richardson([c[k] for c in cs_orc]))
           for k in ("c20", "c11", "c02")}
    two = {k: (2.0 * cs_prod[2][k] - cs_prod[1][k], 2.0 * cs_orc[2][k] - cs_orc[1][k])
           for k in ("c20", "c11", "c02")}
    return out, {"residuals": res, "two_level": two}


class TestG2Heavy(unittest.TestCase):
    """The covering set is ~700 maps of the production solver plus the same
    number of oracle evaluations; every method here is opt-in (see
    ``tests/heavy_tests.py``)."""

    def _check(self, entry, key, coeffs, floor):
        prod, orc = coeffs[key]
        if entry in _KNOWN_ZERO:
            self.assertLess(np.abs(orc).max(), floor, "{}: listed as a symmetry zero".format(entry))
            self.assertLess(np.abs(prod).max(), floor,
                            "{}: oracle is zero but production is not".format(entry))
            return
        self.assertGreater(np.abs(orc).max(), floor,        # anti-vacuity
                           "{}: oracle {} is at the floor; add it to _KNOWN_ZERO with the "
                           "symmetry that makes it vanish".format(entry, key))
        self.assertLess(_rel(prod, orc), 1e-8, "{} {}".format(entry, key))

    def _check_extraction(self, entry, key, coeffs, diag, floor):
        """The extraction itself: the quadratic fit describes the maps, its
        residual improves as the grid is refined, every fitted quantity is
        finite, and the Richardson ladder has SETTLED (the two-level
        estimate from the finer pair agrees with the three-level one).

        Without this, a fit that had stopped describing the sampled maps --
        or a ladder still moving between levels -- would be invisible: the
        production and oracle sides alias the same way, so their AGREEMENT
        survives an extraction that no longer measures a second-order
        coefficient."""
        res = diag["residuals"]
        self.assertTrue(all(np.isfinite(r) for r in res),
                        "{}: the fit residual is not finite".format(entry))
        self.assertLess(max(res), _RESIDUAL_BOUND,
                        "{}: the quadratic fit does not describe the sampled maps "
                        "(relative residuals {})".format(entry, ["{:.2e}".format(r) for r in res]))
        for coarse, fine in zip(res, res[1:]):
            if coarse <= _RESIDUAL_FLOOR:
                continue                # round-off, not the map's cubic term
            self.assertLess(fine, coarse * 0.9,
                            "{}: the fit residual does not improve with the grid "
                            "({:.2e} -> {:.2e})".format(entry, coarse, fine))
        for side, label in ((0, "production"), (1, "oracle")):
            three = coeffs[key][side]
            two = diag["two_level"][key][side]
            self.assertTrue(np.all(np.isfinite(three)) and np.all(np.isfinite(two)),
                            "{}: a fitted coefficient is not finite".format(entry))
            if np.abs(three).max() <= floor:
                # a structural zero (_KNOWN_ZERO): "has the ladder settled"
                # is a question about round-off there, not about the fit
                continue
            self.assertLess(_rel(two, three), _LADDER_BOUND,
                            "{} {}: the Richardson ladder has not settled on the {} "
                            "side".format(entry, key, label))

    @heavy
    def test_b_covering_set_equals_the_local_oracle(self):
        """G2 (b): the standalone general path's second order equals the
        oracle's LOCAL weighting for every pure type and every mixed pair
        of the covering set."""
        for x in _PURE:
            with self.subTest(pure=x):
                c, diag = _coefficients(x, x, False, "local")
                self._check_extraction(x, "c20", c, diag, _FLOOR_PURE)
                self._check(x, "c20", c, _FLOOR_PURE)
        for x, y in _PAIRS:
            with self.subTest(pair=(x, y)):
                c, diag = _coefficients(x, y, False, "local")
                self._check_extraction((x, y), "c11", c, diag, _FLOOR_PAIR)
                self._check((x, y), "c11", c, _FLOOR_PAIR)

    @heavy
    def test_b_dropped_class_load_bearing_for_v(self):
        """The off/off UNCROSSED class that the local weighting drops is not
        numerically negligible for the off-site V of the covering set: were
        it tiny, ``test_b_covering_set_equals_the_local_oracle`` would pass
        against the exact oracle too and would not be pinning the local
        weighting at all."""
        dropped, _ = _coefficients("V", "V", False, "dropped")
        local, _ = _coefficients("V", "V", False, "local")
        self.assertGreater(np.abs(local["c20"][1]).max(), _FLOOR_PURE)       # anti-vacuity
        self.assertGreater(np.abs(dropped["c20"][1]).max(),
                           1e-3 * np.abs(local["c20"][1]).max())

    @heavy
    def test_c_gate_on_reproduces_the_full_oracle(self):
        """G2 (c): with the bond gate on, the off-site exchange crossing the
        local kernel drops is resummed too, so the second order must be the
        FULL (exact) oracle -- asserted here on an ORBITAL-DIAGONAL off-site
        CoulombInter bond (:data:`_V_DIAG`) and on the ORBITAL-OFF-DIAGONAL
        one (``_OFFSITE["V"]``, whose leading rows are ``v_12(+x)``).

        The inter-orbital entry was RECORDED rather than asserted while the
        bond gate's mixed second-order blocks missed it (issue #192); the
        pair permutation of spec 2026-09-16 R3 closed that, and the
        adjudication against exact diagonalisation lives in the chain gate
        ``tests/test_flex_second_order_ed_chain.py``. Off-site Hund/Ising
        stay RECORDED: that class belongs to the chain gate too."""
        c, diag = _coefficients("Vd", "Vd", True, "exact")
        self._check_extraction("Vd", "c20", c, diag, _FLOOR_PURE)
        prod, orc = c["c20"]
        self.assertGreater(np.abs(orc).max(), _FLOOR_PURE)                   # anti-vacuity
        self.assertLess(_rel(prod, orc), 1e-8, "Vd")
        # and the exact-vs-local distinction is load-bearing on this fixture:
        # without recovering the dropped class the assertion above would miss
        # by 2.1e-2, four orders above the tolerance.
        d, _ = _coefficients("Vd", "Vd", True, "dropped")
        l, _ = _coefficients("Vd", "Vd", True, "local")
        self.assertGreater(np.abs(l["c20"][1]).max(), _FLOOR_PURE)           # anti-vacuity
        self.assertGreater(np.abs(d["c20"][1]).max(), 1e-3 * np.abs(l["c20"][1]).max())
        # the ORBITAL-OFF-DIAGONAL off-site bond (_OFFSITE["V"], whose
        # leading rows are v_12(+x)) is asserted the same way since the pair
        # permutation of spec 2026-09-16 R3 closed issue #192: the mixed
        # (channel-0 x bond) second-order blocks now read the bond-side leg
        # and the bond vertex at the transposed orbital pair, which is the
        # exact mixed class on an inter-orbital bond. It used to miss the
        # exact oracle by 1.6e-3 here (4.0e-2 against exact diagonalisation
        # in tests/test_flex_second_order_ed_chain.py, which adjudicates it).
        c, diag = _coefficients("V", "V", True, "exact")
        self._check_extraction("V", "c20", c, diag, _FLOOR_PURE)
        prod, orc = c["c20"]
        self.assertGreater(np.abs(orc).max(), _FLOOR_PURE)                   # anti-vacuity
        self.assertLess(_rel(prod, orc), 1e-8, "V (inter-orbital off-site bond, issue #192)")
        for name in ("JV", "IV"):
            c, diag = _coefficients(name, name, True, "exact")
            self._check_extraction(name, "c20", c, diag, _FLOOR_PURE)
            prod, orc = c["c20"]
            self.assertGreater(np.abs(orc).max(), _FLOOR_PURE)               # anti-vacuity
            rel = _rel(prod, orc)
            expect = _PINNED_GATE[name]
            print("RECORDED gate-on second order for {}: relative deviation {:.3e} "
                  "(recorded {:.3e})".format(name, rel, expect))
            self.assertTrue(np.isfinite(rel), "{}: the recorded deviation is not finite".format(name))
            # at the round-off floor of a 1e-1-sized coefficient: only the
            # ceiling is a statement about the code (see _ROUNDOFF_CEIL)
            self.assertLess(rel, _ROUNDOFF_CEIL,
                            "{}: the bond gate no longer reproduces the exact oracle "
                            "(measured {:.3e}, recorded {:.3e})".format(name, rel, expect))


if __name__ == "__main__":
    unittest.main()
