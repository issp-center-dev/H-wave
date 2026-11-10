"""Exact diagonalization of the spin-orbital bond second order (gate C1a-ED).

The interacting five-site ring of FIXTURE_I (a random complex spin hop at
``R = +-1``, which breaks ``q -> -q`` and S_z, and every interaction type on
site, at ``R = +-1`` and at ``R = +-2``) is diagonalized exactly with the
interaction scaled by a coupling ``x``: ``build_h1(FIXTURE_I) +
build_terms(x tables)``. The Green function comes from ``TotalNED``
(``tests.ed_cache.total``), the density from ``JWED``; neither touches the
solver.

``TestEDEngines`` checks both engines at ``x = 0`` against the free Green
function and density.

``TestFirstOrderED`` (the ED part of the Hartree-Fock check C1a-HF of
``tests.test_so_bond_second_order``) extracts the first-order
coefficient of the exact self-energy ``Sigma_ED = G0^-1 - G^-1`` (central
difference plus one Richardson step) and compares it with the Wick
Hartree-Fock self-energy at the free density.

``TestSecondOrderEDSmoke`` / ``TestSecondOrderED`` (C1a-ED) extract the
second-order coefficient of the remainder ``Sigma_ED - Sigma_HF[rho_ED]``
with the stencils of ``tests.test_flex_second_order_ed_onsite``
(``_coeff2``: ``(f(2x) - 2 f(x)) / (2 x^2)``, then one Richardson step
``_refine``: ``2 c(x/2) - c(x)``) on the inner eighth of the frequency
grid and compare it with the production chain ``p_b(1/2 Gamma_B chibar_B
Gamma_B) -> calc_self_energy_bond`` at the free Green function
(``tests.so_bond_fixtures.c1a_sigma_prod``). The ``"all"`` row is the gate;
the rows of the single interaction types are reported.

The ED side truncates nothing but the comparison does: the production chain
sums a finite frequency grid. The residual ``dev`` of the ``"all"`` row
therefore falls as ``1/Nmat`` (2.60e-3 at Nmat 512, 3.29e-4 at 4096, a
ratio of 7.9 for the expected 8; the smoke point's doubled coupling
adds an extraction error of 1.8e-4 only); ``TestSecondOrderED.test_all_types``
asserts this ratio (window 6 to 10), so the gate's residual is the
truncated frequency sum and not a defect of the formula.

``TestSecondOrderED`` is opt-in (``HWAVE_FULL_TESTS=1``): on numpy 1.26
(the CI numpy) the gate alone takes 8 s on macOS (4.5 s on numpy 2), over
the 5 s rule of ``tests/heavy_tests.py``; the smoke row and its mutation
rows stay in the fast gate.
"""
import functools
import unittest

import numpy as np

import tests.test_flex_second_order_ed_onsite as G3
from tests import ed_cache
from tests.ed_oracle_util import EDFixture
from tests.heavy_tests import heavy
from tests.so_bond_fixtures import c1a_green, c1a_sigma_prod, rel, win
from tests.so_bond_oracles import (FIXTURE_I, JWED, build_terms, free_density, hf_onebody,
                                   k_of_modes, wick_hf_sigma_k)

#: the C1a-ED working point (the gate passed here without refinement rounds)
NMAT_GATE = 4096
X_GATE = 6.25e-3
#: the smoke point
NMAT_SMOKE = 512
X_SMOKE = 0.0125
#: the first-order (C1a-HF, ED part) point: Nmat, coupling, and the centred
#: frequencies on either side of the middle of the window
NMAT_FIRST = 512
X_FIRST = 1e-2
N_FIRST_HALF = 4

#: the seven interaction types of FIXTURE_I
TYPES = tuple(sorted(FIXTURE_I.tables_dict))

_SHAPE = (5, 1, 1)
_NORB = 1


# --------------------------------------------------------------------------
# ED maps
# --------------------------------------------------------------------------

def _ed_tables(row, x):
    """FIXTURE_I's tables scaled by ``x``: every type (``row = "all"``) or the
    type ``row`` alone."""
    tables = FIXTURE_I.tables_dict
    if row != "all" and row not in tables:
        raise ValueError("_ed_tables: unknown row {!r}".format(row))
    return {t: {k: x * v for k, v in tbl.items()}
            for t, tbl in tables.items() if row in ("all", t)}


def _ed_terms(row, x):
    return build_terms(_ed_tables(row, x), 5, _NORB)


def _iws(nmat):
    """The Matsubara frequencies ``i (2 n + 1 - nmat) pi T`` on ``win(nmat)``."""
    n = np.arange(nmat)[win(nmat)]
    return 1j * (2 * n + 1 - nmat) * np.pi / FIXTURE_I.beta


def _g0(nmat):
    """The free Green function ``(nwin, 5, 2, 2)`` on ``win(nmat)``."""
    return c1a_green(FIXTURE_I, "band", nmat)[win(nmat)]


def _ed_green_k(row, x, nmat):
    """The interacting Green function ``(nwin, 5, 2, 2)`` on ``win(nmat)``
    (``TotalNED`` through ``ed_cache.total``)."""
    fx = FIXTURE_I
    efx = EDFixture(5, _NORB, {}, [0.0], 1.0 / fx.beta, fx.mu)
    gm = ed_cache.total(efx, _ed_terms(row, x), fx.h1).green(_iws(nmat))
    return k_of_modes(gm, _SHAPE, _NORB)


def _ed_sigma(row, x, nmat):
    """The exact self-energy ``G0^-1 - G^-1`` ``(nwin, 5, 2, 2)`` on ``win(nmat)``."""
    return np.linalg.inv(_g0(nmat)) - np.linalg.inv(_ed_green_k(row, x, nmat))


@functools.lru_cache(maxsize=None)
def _ed_density(row, x):
    """The exact density ``<c+_p c_q>`` (interleaved modes) of ``JWED``; read-only."""
    fx = FIXTURE_I
    rho = JWED(fx.h1, _ed_terms(row, x), mu=fx.mu, beta=fx.beta).density()
    rho.flags.writeable = False
    return rho


@functools.lru_cache(maxsize=None)
def _ed_hf(row, x, *, hf_density="ed"):
    """The Hartree-Fock self-energy ``(5, 2, 2)`` of the scaled terms at the
    exact density (``"ed"``) or at the free one (``"free"``, a mutation
    hook); read-only."""
    fx = FIXTURE_I
    if hf_density == "ed":
        rho = _ed_density(row, x)
    elif hf_density == "free":
        rho = free_density(fx.h1, fx.mu, fx.beta)
    else:
        raise ValueError("_ed_hf: unknown hf_density {!r}".format(hf_density))
    out = k_of_modes(hf_onebody(_ed_terms(row, x), rho), _SHAPE, _NORB)
    out.flags.writeable = False
    return out


@functools.lru_cache(maxsize=None)
def _ed_remainder(row, x, nmat, *, hf_density="ed"):
    """``Sigma_ED - Sigma_HF`` ``(nwin, 5, 2, 2)`` on ``win(nmat)``; read-only."""
    out = _ed_sigma(row, x, nmat) - _ed_hf(row, x, hf_density=hf_density)[None]
    out.flags.writeable = False
    return out


def _c_ed(row, x, nmat, level, *, hf_density="ed"):
    """The second-order coefficient of the remainder at refinement ``level``:
    ``_refine(_coeff2)`` at the coupling ``x / 2**level``."""
    return G3._refine(
        lambda h: G3._coeff2(
            lambda y: _ed_remainder(row, y, nmat, hf_density=hf_density), h),
        x / 2 ** level)


def _c_array(row, nmat, *, ring_factor=0.5):
    """The production second order of the unit tables on ``win(nmat)``
    (uncached)."""
    g0 = c1a_green(FIXTURE_I, "band", nmat)
    return c1a_sigma_prod(FIXTURE_I, _ed_tables(row, 1.0), g0,
                          ring_factor=ring_factor)[win(nmat)]


def metric_ed(row, nmat, x, *, ring_factor=0.5, hf_density="ed"):
    """The C1a-ED numbers of one row; asserts nothing.

    ``c0``/``c1`` are the ED coefficients at refinement levels 0 and 1,
    ``ca`` the production array. ``scale = max|c1|``; ``dev``,
    ``extraction``, ``offsite`` (the k-dependent part) and ``spinflip``
    (the spin-off-diagonal blocks) are relative to it. When ``finite`` is
    False every ratio is NaN."""
    c0 = _c_ed(row, x, nmat, 0, hf_density=hf_density)
    c1 = _c_ed(row, x, nmat, 1, hf_density=hf_density)
    ca = _c_array(row, nmat, ring_factor=ring_factor)
    scale = float(np.abs(c1).max()) if np.isfinite(c1).all() else float("nan")
    finite = bool(all(np.isfinite(a).all() for a in (c0, c1, ca)) and scale > 0)
    out = {"row": row, "nmat": nmat, "x": x, "finite": finite, "scale": scale}
    if not finite:
        out.update(dev=float("nan"), extraction=float("nan"), offsite=float("nan"),
                   spinflip=float("nan"))
        return out
    out["dev"] = float(np.abs(ca - c1).max()) / scale
    out["extraction"] = float(np.abs(c0 - c1).max()) / scale
    out["offsite"] = float(np.abs(c1 - c1.mean(axis=1, keepdims=True)).max()) / scale
    out["spinflip"] = max(float(np.abs(c1[..., 0, 1]).max()),
                          float(np.abs(c1[..., 1, 0]).max())) / scale
    return out


def first_order_dev(nmat=NMAT_FIRST, x=X_FIRST):
    """The ED part of C1a-HF: the first-order coefficient of ``Sigma_ED`` on the
    ``2 N_FIRST_HALF`` frequencies at the middle of the window (central
    difference plus one Richardson step) against the Wick Hartree-Fock
    self-energy at the free density (max-norm relative); asserts nothing."""
    fx = FIXTURE_I
    nwin = len(_iws(nmat))
    sel = slice(nwin // 2 - N_FIRST_HALF, nwin // 2 + N_FIRST_HALF)

    def D(h):
        return (_ed_sigma("all", h, nmat)[sel] - _ed_sigma("all", -h, nmat)[sel]) / (2 * h)

    c1 = (4 * D(x / 2) - D(x)) / 3
    ref = wick_hf_sigma_k(build_terms(fx.tables_dict, 5, _NORB),
                          free_density(fx.h1, fx.mu, fx.beta), _SHAPE, _NORB)
    return rel(c1, ref[None])


def _assert_finite(tc, m):
    tc.assertTrue(m["finite"] and all(np.isfinite(m[k]) for k in
                                      ("scale", "dev", "extraction", "offsite", "spinflip")),
                  "non-finite metric: {}".format(m))


# --------------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------------

class TestEDEngines(unittest.TestCase):
    """JWED and TotalNED at ``x = 0`` reproduce the free density and Green
    function."""

    def test_free_engines(self):
        fx = FIXTURE_I
        rho = _ed_density("all", 0.0)
        dev_rho = float(np.abs(rho - free_density(fx.h1, fx.mu, fx.beta)).max())
        dev_g = rel(_ed_green_k("all", 0.0, 512), _g0(512))
        self.assertLessEqual(dev_rho, 1e-12, "density {:.3e}".format(dev_rho))
        self.assertLessEqual(dev_g, 1e-12, "Green function {:.3e}".format(dev_g))


class TestFirstOrderED(unittest.TestCase):
    """C1a-HF, ED part (tolerance 1e-6, the extraction's accuracy; the
    array check of C1a-HF is 1e-12): the first-order coefficient of the exact
    self-energy is the Wick Hartree-Fock self-energy."""

    def test_first_order_coefficient(self):
        dev = first_order_dev()
        self.assertTrue(np.isfinite(dev), "non-finite metric")
        self.assertLessEqual(dev, 1e-6, "first-order deviation {:.3e}".format(dev))


class TestSecondOrderEDSmoke(unittest.TestCase):
    """C1a-ED smoke: the ``"all"`` row at Nmat 512, x = 0.0125."""

    def test_all_types_smoke(self):
        m = metric_ed("all", NMAT_SMOKE, X_SMOKE)
        _assert_finite(self, m)
        self.assertGreater(m["scale"], 1e-6, m)
        self.assertLessEqual(m["dev"], 2e-2, m)
        self.assertLessEqual(m["extraction"], 5e-4, m)
        self.assertGreaterEqual(m["offsite"], 0.1, m)
        self.assertGreaterEqual(m["spinflip"], 0.1, m)

    def test_mutation_rows(self):
        # E1: the ring at its full staged factor 1 (instead of 1/2)
        changed1 = rel(_c_array("all", NMAT_SMOKE, ring_factor=1.0), _c_array("all", NMAT_SMOKE))
        m1 = metric_ed("all", NMAT_SMOKE, X_SMOKE, ring_factor=1.0)
        # E2: the Hartree-Fock subtraction at the free density (instead of the exact one)
        changed2 = rel(_ed_hf("all", X_SMOKE, hf_density="free"), _ed_hf("all", X_SMOKE))
        m2 = metric_ed("all", NMAT_SMOKE, X_SMOKE, hf_density="free")
        report = ("E1 ring_factor=1: applied %.3e, dev %.3e; E2 hf_density=free: applied %.3e, "
                  "dev %.3e" % (changed1, m1["dev"], changed2, m2["dev"]))
        print("\nC1a-ED mutation rows (Nmat %d, x %g): %s" % (NMAT_SMOKE, X_SMOKE, report))
        with self.subTest(row="E1"):
            self.assertGreaterEqual(changed1, 1e-3, report)
            self.assertFalse(np.isfinite(m1["dev"]) and m1["dev"] <= 2e-2, report)
        with self.subTest(row="E2"):
            self.assertGreaterEqual(changed2, 1e-6, report)
            self.assertFalse(np.isfinite(m2["dev"]) and m2["dev"] <= 2e-2, report)


class TestSecondOrderED(unittest.TestCase):
    """C1a-ED at the working point (Nmat 4096, x = 6.25e-3): the ``"all"``
    row is the gate, the single-type rows are reported."""

    @heavy
    def test_all_types(self):
        m = metric_ed("all", NMAT_GATE, X_GATE)
        print("C1a-ED gate:", m)
        _assert_finite(self, m)
        self.assertGreater(m["scale"], 1e-6, m)
        self.assertLessEqual(m["dev"], 2e-3, m)
        self.assertLessEqual(m["extraction"], 5e-4, m)
        self.assertGreaterEqual(m["offsite"], 0.1, m)
        self.assertGreaterEqual(m["spinflip"], 0.1, m)
        # the residual is the production chain's truncated frequency sum:
        # it falls as 1/Nmat between the smoke point and the gate (ratio 8)
        smoke = metric_ed("all", NMAT_SMOKE, X_SMOKE)
        _assert_finite(self, smoke)
        ratio = smoke["dev"] / m["dev"]
        print("C1a-ED dev ratio Nmat %d / %d: %.2f" % (NMAT_SMOKE, NMAT_GATE, ratio))
        self.assertTrue(6.0 <= ratio <= 10.0, (ratio, smoke["dev"], m["dev"]))

    @heavy
    def test_each_type(self):
        rows = [metric_ed(t, 4096, 6.25e-3) for t in TYPES]
        for m in rows:
            print("C1a-ED {:>12s}: dev {:.3e} extraction {:.3e} scale {:.3e}".format(
                m["row"], m["dev"], m["extraction"], m["scale"]))
        for m in rows:
            self.assertTrue(m["finite"], "non-finite metric: {}".format(m))
            self.assertGreater(m["scale"], 1e-6, m)


if __name__ == "__main__":
    unittest.main()
