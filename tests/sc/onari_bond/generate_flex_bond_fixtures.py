#!/usr/bin/env python3
"""Regenerate the bond-resolved Hartree-Fock FLEX (#181 Phase B) Green
functions of the Onari-type trend milestone
(``tests/test_flex_bond_onari_trend.py``) and the compact observables file
``flex_bond_observables.json`` that is committed in their place.

Setting (frozen, spec 2026-09-06 section 5 "Onari-type trend milestone")
------------------------------------------------------------------------
Single-band square lattice, ``t = 1`` nearest neighbour only, ``U = 4``,
``V`` on the four nearest-neighbour bonds, filling ``n = 0.7``, ``T = 0.02``,
``L = 16``, ``Nmat = 2048``, ``V in {0, 0.4, 0.8, 1.0, 1.2}``; FLEX
``calc_scheme = "general"`` with ``flex_hartree_fock = true`` and
``longitudinal_bond_channels = true``, Anderson mixing (depth 8), ``Mix =
0.2``, ``EPS = 8`` (the three residuals of the split state, three
consecutive passes), ``IterationMax = 1500``.  Every V point is an
independent start from the converged self-energy of the standard general
FLEX at the same V ("cold" below: at ``U = 4``, ``beta = 50`` the bare Green
function is inside the RPA spin instability and the bond path's conditioning
check refuses a bare start, so the legacy solution -- no Hartree-Fock term,
Hartree-only vertex -- is converted with ``hwave_sigma_split --zero-static``
and used as the seed); the two points nearest the d-to-f crossing
(``V = 1.0``, ``1.2``) are ALSO run warm-started from the converged split
self-energy of the adjacent lower V (the two-seed check against metastable
SCF branches).

Consumer: the EXISTING static bond Eliashberg chain
(``sc._build_bond_operator`` fed with the FLEX-bond ``green.npz``, exactly
as ``[eliashberg] bond_green`` does), tracked through the odd-parity f-like
triplet subspace with the milestone's own helpers.

What is committed
-----------------
Only ``flex_bond_observables.json`` (a few kB): the settings, per run the
convergence record (iterations, the three residuals, the density-closure
error, the chemical potential, the conditioning minima, wall time and peak
RSS) and the tracked ``lambda_t``.  The Green functions (~1 MB each) are
regenerated on demand into ``tests/sc/onari_bond/_regenerated_flexbond/``
(git-ignored; override with ``HWAVE_FLEXBOND_DIR``) by::

    HWAVE_RUN_SLOW_FIXTURES=1 python3 -m unittest tests.test_flex_bond_onari_trend

or directly with ``python3 tests/sc/onari_bond/generate_flex_bond_fixtures.py``.
Neither FLEX nor ``np.savez_compressed`` is bit-reproducible, so regenerated
files are not hash-pinned; the physics is pinned by the lambda table in the
test at ``LAMBDA_RTOL``.
"""
import argparse
import hashlib
import json
import os
import resource
import shutil
import sys
import tempfile
import time

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))
from tests.sc.onari_bond.generate_fixtures import _write_inputs, _symmetrize   # noqa: E402

U = 4.0
T = 0.02
FILLING = 0.7
L = 16
NMAT = 2048
V_GRID = (0.0, 0.4, 0.8, 1.0, 1.2)
WARM_FROM = {1.0: 0.8, 1.2: 1.0}
EPS = 8
ITERATION_MAX = 1500
MIX = 0.2
DEPTH = 8

SETTINGS = dict(U=U, T=T, filling=FILLING, L=L, Nmat=NMAT, V_grid=list(V_GRID),
                warm_from={str(k): v for k, v in WARM_FROM.items()}, EPS=EPS,
                IterationMax=ITERATION_MAX, Mix=MIX, anderson_depth=DEPTH,
                scheme="general", flex_hartree_fock=True, longitudinal_bond_channels=True)

OBSERVABLES = "flex_bond_observables.json"
DEFAULT_DIR = os.environ.get("HWAVE_FLEXBOND_DIR",
                             os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                          "_regenerated_flexbond"))


def green_name(V, seed):
    return "green_flexbond_L{}_V{:.2f}_{}.npz".format(L, V, seed)


def _flex_input(indir, outdir, sigma_init=None, iteration_max=ITERATION_MAX, phase_b=True):
    inp = {"path_to_input": "",
           "interaction": {"path_to_input": indir, "Geometry": "geom.dat",
                           "Transfer": "transfer.dat", "CoulombIntra": "coulombintra.dat",
                           "CoulombInter": "coulombinter.dat"}}
    if sigma_init is not None:
        inp["path_to_input"] = os.path.dirname(sigma_init)
        inp["sigma_init"] = os.path.basename(sigma_init)
    param = {"T": T, "filling": FILLING, "CellShape": [L, L, 1], "SubShape": [1, 1, 1],
             "Nmat": NMAT, "IterationMax": iteration_max, "Mix": MIX, "EPS": EPS,
             "mixing_scheme": "anderson", "anderson_depth": DEPTH}
    if phase_b:
        param.update(flex_hartree_fock=True, longitudinal_bond_channels=True)
    return {
        "log": {"print_level": 1},
        "mode": {"mode": "FLEX", "calc_scheme": "general", "param": param},
        "file": {"input": inp,
                 "output": {"path_to_output": outdir, "green": "green", "sigma": "sigma",
                            "energy": "energy.dat"}},
    }


def legacy_seed(V, outdir, iteration_max=ITERATION_MAX):
    """The 'cold' start of the Phase B run: at U = 4, beta = 50 the BARE
    Green function sits inside the RPA spin instability (the bond path's
    conditioning check refuses the very first map), so every V point is
    seeded from the converged self-energy of the standard general FLEX at
    the same V (no Hartree-Fock term, Hartree-only vertex), converted to
    the split form with a zero static part by the hwave_sigma_split
    machinery.  Returns (seed path, legacy iteration count)."""
    import hwave.qlms as qlms
    from hwave.solver.flex_hf import sigma_split_convert
    work = tempfile.mkdtemp(prefix="flexbond_legacy_")
    try:
        indir = os.path.join(work, "in")
        _write_inputs(indir, V)
        result = qlms.run(input_dict=_flex_input(indir, work, None, iteration_max, phase_b=False))
        assert result.get("scf_converged"), "legacy FLEX did not converge at V={}".format(V)
        os.makedirs(outdir, exist_ok=True)
        seed = os.path.join(outdir, "seed_legacy_L{}_V{:.2f}.npz".format(L, V))
        sigma_split_convert(os.path.join(work, "sigma.npz"), seed, zero_static=True, force=True)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return seed, int(result.get("scf_iterations", -1))


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run_point(V, outdir, seed="cold", sigma_init=None, iteration_max=ITERATION_MAX):
    """One FLEX run; writes the symmetrised green and returns (record, sigma_path)."""
    import hwave.qlms as qlms
    work = tempfile.mkdtemp(prefix="flexbond_")
    t0 = time.time()
    try:
        indir = os.path.join(work, "in")
        _write_inputs(indir, V)
        result = qlms.run(input_dict=_flex_input(indir, work, sigma_init, iteration_max))
        raw = np.load(os.path.join(work, "green.npz"))
        prov = {k: raw[k] for k in ("scf_converged", "scf_iterations", "scf_sigma_residual",
                                    "scf_green_residual", "scf_component_residual",
                                    "hf_density_error")}
        cs = np.load(os.path.join(work, "chiq_s.npz"))
        cond = {"cond_min_s": float(cs["longitudinal_bond_cond_min_s"]),
                "cond_min_c": float(cs["longitudinal_bond_cond_min_c"])}
        mu = None
        with open(os.path.join(work, "energy.dat")) as f:
            for ln in f:
                if ln.startswith("ChemicalPotential"):
                    mu = float(ln.split("=")[1])
        green, delta = _symmetrize(raw["green"], L)
        os.makedirs(outdir, exist_ok=True)
        gpath = os.path.join(outdir, green_name(V, seed))
        np.savez_compressed(gpath, green=green, L=L, V=V, U=U, T=T, filling=FILLING, nmat=NMAT,
                            seed=seed, symmetrization_residual=delta, **prov)
        spath = os.path.join(outdir, "sigma_L{}_V{:.2f}_{}.npz".format(L, V, seed))
        shutil.copy(os.path.join(work, "sigma.npz"), spath)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    wall = time.time() - t0
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform != "darwin":
        rss *= 1024
    rec = dict(V=V, seed=seed, green=os.path.basename(gpath), sha256=_sha256(gpath),
               scf_converged=bool(prov["scf_converged"]), scf_iterations=int(prov["scf_iterations"]),
               scf_sigma_residual=float(prov["scf_sigma_residual"]),
               scf_green_residual=float(prov["scf_green_residual"]),
               scf_component_residual=float(prov["scf_component_residual"]),
               hf_density_error=float(prov["hf_density_error"]), mu=mu,
               symmetrization_residual=float(delta), wall_s=wall, maxrss_mb=rss / 2 ** 20,
               qlms_result={k: (bool(v) if isinstance(v, (bool, np.bool_)) else int(v))
                            for k, v in result.items()}, **cond)
    return rec, spath


def generate(outdir=DEFAULT_DIR, iteration_max=ITERATION_MAX, v_grid=V_GRID, warm=True):
    """Run every point (cold, then the warm starts) and write the observables."""
    from tests.test_flex_bond_onari_trend import lambda_sweep
    records = []
    sigmas = {}
    for V in v_grid:
        seed, legacy_iterations = legacy_seed(V, outdir, iteration_max)
        rec, spath = run_point(V, outdir, "cold", seed, iteration_max)
        rec["legacy_seed_iterations"] = legacy_iterations
        assert rec["scf_converged"], "FLEX did not converge at V={} ({} iterations)".format(
            V, rec["scf_iterations"])
        records.append(rec)
        sigmas[V] = spath
        print("V={:.2f} cold: {} iterations, {:.0f} s, lambda pending".format(
            V, rec["scf_iterations"], rec["wall_s"]), flush=True)
    if warm:
        for V, V0 in WARM_FROM.items():
            if V not in v_grid or V0 not in sigmas:
                continue
            rec, _ = run_point(V, outdir, "warm_from_V{:.2f}".format(V0), sigmas[V0], iteration_max)
            assert rec["scf_converged"], "warm FLEX did not converge at V={}".format(V)
            records.append(rec)
            print("V={:.2f} warm from {:.2f}: {} iterations, {:.0f} s".format(
                V, V0, rec["scf_iterations"], rec["wall_s"]), flush=True)
    lam_cold = lambda_sweep([os.path.join(outdir, green_name(V, "cold")) for V in v_grid], list(v_grid))
    lam = {"cold": {"{:.2f}".format(V): lam_cold["lambda"][i] for i, V in enumerate(v_grid)},
           "warm": {}}
    for rec in records:
        if rec["seed"] != "cold":
            V = rec["V"]
            # track the warm green as the last point of the sweep up to V
            seq = [os.path.join(outdir, green_name(v, "cold")) for v in v_grid if v < V]
            seq.append(os.path.join(outdir, rec["green"]))
            vs = [v for v in v_grid if v < V] + [V]
            lam["warm"]["{:.2f}".format(V)] = lambda_sweep(seq, vs)["lambda"][-1]
    obs = dict(settings=SETTINGS, records=records, lambda_t=lam, tracking=lam_cold["tracking"])
    with open(os.path.join(outdir, OBSERVABLES), "w") as f:
        json.dump(obs, f, indent=1, sort_keys=True)
    print(json.dumps(lam, indent=1))
    return obs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("outdir", nargs="?", default=DEFAULT_DIR)
    parser.add_argument("--iteration-max", type=int, default=ITERATION_MAX)
    parser.add_argument("--only", type=float, nargs="*", default=None, help="restrict the V grid")
    parser.add_argument("--no-warm", action="store_true")
    a = parser.parse_args()
    generate(a.outdir, a.iteration_max, tuple(a.only) if a.only else V_GRID, warm=not a.no_warm)
