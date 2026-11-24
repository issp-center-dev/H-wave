"""First-order spinful bond correlators against exact diagonalization.

The longitudinal bond bilinear ``B_{R,a,b,sigma}(q) = sum_j e^{+iq.j}
c+_{j+R,a sigma} c_{j,b sigma}`` and the transverse one
``S+_{R,a,b}(q) = sum_j e^{+iq.j} c+_{j+R,a up} c_{j,b dn}`` are label
vectors of ``rho^{(d)}_{XY}(q) = sum_j e^{+iq.j} c+_{X,j} c_{Y,j+d}``
(``tests/so_bond_oracles.bond_observable``). Their correlators, evaluated on
the table-driven exact diagonalization, are first pinned to the sector
exact diagonalization of ``tests/ed_oracle_util.py`` on an S_z-conserving
fixture. The first-order vertex part of random combinations within each
family (full first order minus the Hartree-Fock-shifted first order) is
then compared with ``u^T [-chibar(q) Gamma(-q) chibar(q)] v`` built from the
production spinful bond vertex ``so_bond.so_bond_vertex``.
"""

import unittest

import numpy as np

from tests.ed_oracle_util import EDFixture, SectorED
from tests.heavy_tests import heavy
from tests.so_bond_fixtures import chibar_static, production_first_order
from tests.so_bond_oracles import (FAMILIES, FIXTURE_I, FIXTURE_I_SZ, FIXTURE_I_ZEEMAN,
                                   FIXTURE_II, FIXTURE_X, NEEDS_RICHARDSON, rel,
                                   bond_observable, build_terms, ed_cached, family_sketches,
                                   fixture_name, label_list, sketch_matrix, sketch_set,
                                   vertex_part_sketches)

#: the first-order exact-diagonalization tolerance (the central-difference floor)
CHI_ED_TOL = 1e-6
#: the sector exact diagonalization agrees with the table-driven one to this
REEXPRESSION_TOL = 1e-10
#: the coupling of the re-expression check
REEXPRESSION_LAM = 0.3
#: every family's first-order part must reach this fraction of its bare bubble
FAMILY_FLOOR = 1e-3
#: the sector-targeted sketches must reach this fraction of the largest
#: random family sketch (every fixture here carries every interaction type)
TARGETED_FLOOR = 1e-2
#: random sketch pairs per family
N_RANDOM = 8
#: the sector-targeted pairs of ``sketch_set`` follow its eight random ones
_SKETCH_SET_TARGETED = slice(8, 11)


def _channels(fixture):
    return [(R,) for R in range(fixture.nvol)]


def _oracle_observables(fixture, *, lam):
    """``(longitudinal[q, s, s', I, J], transverse[q, I, J])``: the static
    correlators ``(1/N) int <O_I(q, tau) O_J(q)^dagger>_conn`` of the bond
    bilinears of every channel ``R`` (norb_phys = 1), as sketches of the
    table-driven exact diagonalization of ``h1 + lam H_int``. The
    observables are built on orbital 0 only, so a fixture with more physical
    orbitals is refused."""
    if fixture.norb_phys != 1:
        raise ValueError("_oracle_observables: the bond observables cover norb_phys = 1 only, "
                         "got {}".format(fixture.norb_phys))
    ed = ed_cached(fixture, lam)
    labels = label_list(fixture.geom, fixture.nd)
    Rs = [c[0] for c in _channels(fixture)]
    n, L = len(Rs), fixture.nvol
    lon = np.empty((L, 2, 2, n, n), dtype=np.complex128)
    tra = np.empty((L, n, n), dtype=np.complex128)
    for q in range(L):
        # rows: B_{R,s} for s = 0, 1 then S+_R; columns: their conjugate legs
        obs = [bond_observable(fixture, q, R, 0, 0, s, s) for s in (0, 1) for R in Rs]
        obs += [bond_observable(fixture, q, R, 0, 0, 0, 1) for R in Rs]
        M = sketch_matrix(ed, [u for u, _ in obs], [v for _, v in obs], q, labels,
                          fixture.L, fixture.norb_phys)
        lon[q] = M[:2 * n, :2 * n].reshape(2, n, 2, n).transpose(0, 2, 1, 3)
        tra[q] = M[2 * n:, 2 * n:]
    return lon, tra


def _sector_ed_observables(fixture, *, lam):
    """The same correlators from ``SectorED.bond_correlator`` /
    ``bond_correlator_transverse`` of ``h1 + lam H_int`` (the fixture must
    conserve S_z)."""
    fx = EDFixture(fixture.nvol, fixture.norb_phys, {}, [0.0] * fixture.norb_phys,
                   1.0 / fixture.beta, fixture.mu)
    terms = [(p, q, r, s, lam * v)
             for (p, q, r, s, v) in build_terms(fixture.tables_dict, fixture.L, fixture.norb_phys)]
    se = SectorED(fx, terms, h1=fixture.h1)
    return se.bond_correlator(_channels(fixture)), se.bond_correlator_transverse(_channels(fixture))


def _chi_first_order(fixture, *, qs, n_random):
    """``(dev, family_ratio, targeted, table)`` of the first-order check.

    ``dev``: per family ``max |ED - pred| / max |ED|`` over its sketches and
    ``qs``, where ``ED`` is the first-order vertex part of the exact bond
    correlators and ``pred = u^T [-chibar(q) Gamma(-q) chibar(q)] v`` with
    the production vertex in the label form; ``family_ratio``: per family
    ``max |ED| / max |u^T chibar v|``; ``targeted``: each sector-targeted
    sketch of ``sketch_set`` (``(sigma, sigma)`` when the fixture has
    PairLift) over the largest random family sketch. The targeted sketch of
    sector 0 counts in the longitudinal deviation, those of sectors +-1 in
    the transverse one; ``table`` holds the per-sketch residuals (shown on
    failure)."""
    name = fixture_name(fixture)
    qs = tuple(int(q) for q in qs)
    fams = {f: family_sketches(fixture, f, n_random) for f in sorted(FAMILIES)}
    targeted = sketch_set(fixture, "fast")[_SKETCH_SET_TARGETED]
    sketches = [p for f in sorted(fams) for p in fams[f]] + list(targeted)
    ref = vertex_part_sketches(fixture, sketches=sketches, qs=qs,
                               richardson=NEEDS_RICHARDSON[name], x=FIXTURE_X[name])
    pred = production_first_order(fixture, fixture, sketches, qs)
    chibar = chibar_static(fixture)
    bare = np.array([[u @ chibar[q] @ v for (u, v) in sketches] for q in qs])
    spans, lo = {}, 0
    for f in sorted(fams):
        spans[f] = slice(lo, lo + len(fams[f]))
        lo += len(fams[f])
    # the targeted sectors 0, +1, -1: sector 0 is longitudinal, +-1 transverse
    spans_t = {"longitudinal": [lo], "transverse": [lo + 1, lo + 2]}
    dev, ratio = {}, {}
    for f, sl in spans.items():
        cols = list(range(sl.start, sl.stop)) + spans_t[f]
        r, p = ref[:, cols], pred[:, cols]
        dev[f] = float(np.abs(p - r).max() / np.abs(r).max())
        ratio[f] = float(np.abs(ref[:, sl]).max() / np.abs(bare[:, sl]).max())
    rnd = np.abs(ref[:, :lo]).max()
    tgt = [float(np.abs(ref[:, i]).max() / rnd) for i in range(lo, ref.shape[1])]
    lines = ["fixture {} qs {}: dev {}, family ratio {}, targeted {}".format(
        name, qs, {f: "{:.3e}".format(d) for f, d in dev.items()},
        {f: "{:.3e}".format(d) for f, d in ratio.items()}, ["{:.3e}".format(t) for t in tgt])]
    kinds = [f for f in sorted(fams) for _ in fams[f]] + ["targeted"] * len(targeted)
    for s in range(ref.shape[1]):
        lines.append("  sketch {:2d} {:12s} |ED| {:.3e}  |bare| {:.3e}  |ED - pred| {:.3e}".format(
            s, kinds[s], np.abs(ref[:, s]).max(), np.abs(bare[:, s]).max(),
            np.abs(pred[:, s] - ref[:, s]).max()))
    return dev, ratio, tgt, "\n".join(lines)


class TestObservableReexpression(unittest.TestCase):
    """A1-chi ED: the bond bilinears ``B`` and ``S+`` re-expressed as label
    vectors and evaluated on the table-driven exact diagonalization equal
    ``SectorED.bond_correlator`` / ``bond_correlator_transverse`` on an
    S_z-conserving fixture (no spin-flip hopping, no PairLift, Zeeman on)."""

    def test_sketches_equal_sector_ed(self):
        lon_o, tra_o = _oracle_observables(FIXTURE_I_SZ, lam=REEXPRESSION_LAM)
        lon_s, tra_s = _sector_ed_observables(FIXTURE_I_SZ, lam=REEXPRESSION_LAM)
        got = np.concatenate([lon_o.ravel(), tra_o.ravel()])
        want = np.concatenate([lon_s.ravel(), tra_s.ravel()])
        dev = rel(got, want)
        if not dev <= REEXPRESSION_TOL:
            self.fail("re-expressed correlators vs the sector ED: {:.3e} (longitudinal {:.3e}, "
                      "transverse {:.3e})".format(dev, rel(lon_o, lon_s), rel(tra_o, tra_s)))
        # the transverse correlator is non-trivial and differs from the longitudinal one
        self.assertGreaterEqual(np.abs(tra_s).max(), 1e-2 * np.abs(lon_s).max())
        self.assertGreaterEqual(rel(tra_s, lon_s[:, 0, 0]), 1e-3)


class TestChiFirstOrderED(unittest.TestCase):
    """A1-chi ED: the first-order vertex part of the longitudinal and
    transverse bond correlators from exact diagonalization against
    ``u^T [-chibar(q) Gamma_B(-q) chibar(q)] v`` with the production spinful
    bond vertex mapped to the label form, by random sketch pairs within each
    family plus the sector-targeted sketches."""

    def _check(self, fixture, *, qs, n_random=N_RANDOM):
        dev, ratio, targeted, table = _chi_first_order(fixture, qs=qs, n_random=n_random)
        ok = (all(d <= CHI_ED_TOL for d in dev.values())
              and all(r >= FAMILY_FLOOR for r in ratio.values())
              and len(targeted) == 3 and all(t >= TARGETED_FLOOR for t in targeted))
        if not ok:
            self.fail(table)
        return max(dev.values())

    def test_fixture_i(self):
        self.assertLessEqual(self._check(FIXTURE_I, qs=(1, 2)), CHI_ED_TOL)

    def test_fixture_i_zeeman(self):
        self.assertLessEqual(self._check(FIXTURE_I_ZEEMAN, qs=(1,)), CHI_ED_TOL)

    @heavy
    def test_fixture_ii(self):
        self.assertLessEqual(self._check(FIXTURE_II, qs=(1, 2)), CHI_ED_TOL)


if __name__ == "__main__":
    unittest.main()
