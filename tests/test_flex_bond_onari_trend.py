"""Onari-type trend milestone of the bond-resolved Hartree-Fock FLEX (#181
Phase B, spec 2026-09-06 section 5): the f-like triplet eigenvalue
``lambda_t`` of the EXISTING static bond Eliashberg chain, fed with the
Phase B ``green.npz`` of each V, rises monotonically with the inter-site
V on the single-band Onari case; every SCF converged by the three-residual
rule with the density closure ``hf_density_error < 1e-6``; a cold start
and a warm start from the adjacent V agree at the two points nearest the
d-to-f crossing (no metastable SCF branch).

Setting and generator: ``tests/sc/onari_bond/generate_flex_bond_fixtures.py``
(``U = 4``, ``n = 0.7``, ``T = 0.02``, ``L = 16``, ``Nmat = 2048``, Anderson
depth 8, ``Mix = 0.2``, ``EPS = 8``, ``IterationMax = 1500``).  The Green
functions are NOT committed; the compact observables file
``tests/sc/onari_bond/flex_bond_observables.json`` is (its per-file SHA-256
values are informational: neither FLEX nor ``np.savez_compressed`` is
bit-reproducible, so the physics is pinned by the lambda table, not by hashes).  By default this
module verifies that file against the pinned tables below; with
``HWAVE_RUN_SLOW_FIXTURES=1`` it regenerates every Green function (hours)
and re-derives the observables.
"""
import json
import os
import unittest

import numpy as np

import hwave.sc as sc
from hwave.solver import bond_channels as bc
from tests.heavy_tests import heavy
from tests.test_bond_onari_milestone import (_odd_spectrum, _kgrid, _interactions, DEG_TOL,
                                             BETA, LAMBDA_RTOL, ATOL_MONO)

_SLOW_ENV = "HWAVE_RUN_SLOW_FIXTURES"
FIXTURE_DIR = os.path.join("tests", "sc", "onari_bond")
OBSERVABLES = os.path.join(FIXTURE_DIR, "flex_bond_observables.json")

#: tracked f-like triplet lambda_t of the Phase B greens (cold starts),
#: recorded from the generator run of 2026-09-06 (see the observables file;
#: every point converged in 10-25 iterations from the legacy seed, the warm
#: starts at V = 1.0 and 1.2 agree with the cold ones to 1e-9 relative).
LAMBDA_FLEX_BOND_16 = {
    "0.00": 0.07975703202716758,
    "0.40": 0.08178290472738786,
    "0.80": 0.09940861974554759,
    "1.00": 0.13084753641943114,
    "1.20": 0.20422236554811116,
}
HF_DENSITY_TOL = 1.0e-6


def _point_from_green(green_path, V):
    """One bond-path point from a Phase B green (sc layout), built by the
    production ``sc._build_bond_operator`` exactly as the milestone does."""
    raw = np.load(green_path)["green"]
    nmat = raw.shape[1]
    Lg = int(round(np.sqrt(raw.shape[2])))
    green = raw[0].reshape(nmat, Lg, Lg, 1, 1, 1).transpose(4, 5, 1, 2, 3, 0).copy()
    kx, ky, kz = _kgrid(Lg)
    it = _interactions(V, True)
    inter_k = sc._build_interaction_k(kx, ky, kz, it, 1)
    bond_set = bc.resolve_interactions(it.get("CoulombInter", {}), np.eye(3), 1)
    carrier = sc._BondGreen(full_sc=green, deflated_kw=sc._green_sc_to_canonical(green),
                            green0_tail=None, coeff_tail_requested=0.0, coeff_tail_applied=0.0,
                            source="external_npz")
    (A, vec_size), provenance, attribution = sc._build_bond_operator(
        bond_set, carrier, it, inter_k, {"norb": 1, "rvec": np.eye(3)}, 1, kx, ky, kz, BETA,
        "triplet", bond_max_shells=None, bond_memory_cap_gb=0.0, g2_tail=False)
    evals, vecs = _odd_spectrum(A, vec_size, attribution["weight"], Lg)
    return {"evals": evals, "vecs": vecs, "weight": attribution["weight"], "L": Lg}


def lambda_sweep(green_paths, v_grid):
    """Invariant-subspace-tracked lambda_t along ``v_grid`` (milestone
    tracking: f_x seed, odd harmonic basis)."""
    points = [_point_from_green(p, V) for p, V in zip(green_paths, v_grid)]
    Lg = points[0]["L"]
    kx, ky, kz = _kgrid(Lg)
    seed = sc._initialize_gap("f_x", 1, kx, ky, kz).ravel()
    basis = sc._odd_harmonic_basis(1, kx, ky, kz)
    records = bc.track_subspace([(p["evals"], p["vecs"]) for p in points], seed=seed,
                                weight=points[0]["weight"], basis=basis, deg_tol=DEG_TOL,
                                capture_tol=0.2)
    return {"lambda": [float(r["lambda"]) for r in records],
            "tracking": [{"V": float(V), "dim": int(r["dim"]), "overlap": float(r.get("overlap", 1.0))}
                         for V, r in zip(v_grid, records)]}


def _accept(obs, testcase):
    """The acceptance rules on an observables dict (committed or fresh)."""
    from tests.sc.onari_bond import generate_flex_bond_fixtures as gen
    testcase.assertEqual(obs["settings"], json.loads(json.dumps(gen.SETTINGS)))
    v_grid = [float(v) for v in obs["settings"]["V_grid"]]
    eps = 10.0 ** (-obs["settings"]["EPS"])
    for rec in obs["records"]:
        testcase.assertTrue(rec["scf_converged"], rec)
        for k in ("scf_sigma_residual", "scf_green_residual", "scf_component_residual"):
            testcase.assertLess(rec[k], eps, (rec["V"], rec["seed"], k))
        testcase.assertLess(rec["hf_density_error"], HF_DENSITY_TOL, (rec["V"], rec["seed"]))
        testcase.assertGreater(min(rec["cond_min_s"], rec["cond_min_c"]), 0.0)
    lam = [obs["lambda_t"]["cold"]["{:.2f}".format(V)] for V in v_grid]
    for i in range(1, len(lam)):
        testcase.assertGreaterEqual(lam[i], lam[i - 1] - ATOL_MONO,
                                    "lambda_t not monotone at V={}".format(v_grid[i]))
    testcase.assertGreater(lam[-1], lam[0])
    for tr in obs["tracking"]:
        testcase.assertEqual(tr["dim"], 2)
    for key, lw in obs["lambda_t"]["warm"].items():
        lc = obs["lambda_t"]["cold"][key]
        testcase.assertLess(abs(lw - lc), LAMBDA_RTOL * abs(lc), "two-seed check at V={}".format(key))
    return lam, v_grid


class TestFlexBondOnariTrend(unittest.TestCase):

    def test_committed_observables_pass_the_acceptance(self):
        with open(OBSERVABLES) as f:
            obs = json.load(f)
        lam, v_grid = _accept(obs, self)
        for V, value in zip(v_grid, lam):
            ref = LAMBDA_FLEX_BOND_16["{:.2f}".format(V)]
            self.assertLess(abs(value - ref), LAMBDA_RTOL * abs(ref), "V={}".format(V))

    @heavy
    def test_regenerated_greens_reproduce_the_pinned_lambda(self):
        if os.environ.get(_SLOW_ENV, "").strip().lower() in ("", "0", "false", "no", "off"):
            self.skipTest("set {}=1 to regenerate the Phase B Onari greens (hours)".format(_SLOW_ENV))
        from tests.sc.onari_bond import generate_flex_bond_fixtures as gen
        obs = gen.generate(gen.DEFAULT_DIR)
        lam, v_grid = _accept(obs, self)
        for V, value in zip(v_grid, lam):
            ref = LAMBDA_FLEX_BOND_16["{:.2f}".format(V)]
            self.assertLess(abs(value - ref), LAMBDA_RTOL * abs(ref), "V={}".format(V))


if __name__ == "__main__":
    unittest.main()
