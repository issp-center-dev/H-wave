#!/usr/bin/env python3
"""Regenerate the checked-in FLEX ``green.npz`` fixtures of the single-band
Onari bond-channel acceptance milestone (``tests/test_bond_onari_milestone.py``,
spec S7.7).

Which fixtures are committed
---------------------------
Only the **L = 16** set is committed (~1.3 MB); it is content-hash-verified by
``test_fixtures_are_intact_and_carry_the_documented_setup`` and backs every
default-run milestone test.  The **L = 32** set (~3.9 MB) is deliberately NOT
committed -- it backs a single grid-convergence test, and the repository keeps
no other tracked file above ~150 KB.  ``test_grid_convergence_16_to_32``
regenerates it on demand by calling :func:`generate` below, and is skipped by
default.  To run it::

    HWAVE_RUN_SLOW_FIXTURES=1 pytest tests/test_bond_onari_milestone.py \
        -k grid_convergence
    # or, equivalently, select the marker:
    pytest -m slow tests/test_bond_onari_milestone.py

The regenerated files land in ``tests/sc/onari_bond/_regenerated/`` (override
with ``HWAVE_ONARI_L32_DIR``); that directory is git-ignored and is reused on
subsequent runs, so the FLEX cost is paid once.  Because neither FLEX nor
``np.savez_compressed`` is bit-reproducible, regenerated L = 32 files are
metadata- and physics-verified but NOT hash-pinned.

Apart from that, this script is **not** run by CI -- the committed L = 16
fixtures are content-hash-verified by the test.  It exists so the fixtures are
reproducible from a recorded identity rather than from a path (spec S7.7,
"Reproducibility record").

Model (single-band square lattice, full C4v point group)
-------------------------------------------------------
* orthorhombic unit lattice, 1 orbital at the origin;
* nearest-neighbour hopping only, ``t = 1`` (no ``t'``: a diagonal hopping on
  only ONE of the two diagonals would break C4 and lift the odd-parity E
  doublet the milestone tracks);
* ``U = 4`` on site, ``V`` on the four nearest-neighbour bonds;
* ``filling n = 0.7``, ``T = 0.02`` (``beta = 50``).

FLEX run (per V point, independently self-consistent -- Onari's Green function
is self-consistent, so the milestone regenerates it at every V rather than
freezing one green across the sweep)
------------------------------------------------------------------------------
* ``mode = "FLEX"``, ``calc_scheme = "reduced"``;
* ``CellShape = [L, L, 1]``, ``SubShape = [1, 1, 1]``, ``Nmat = 256``
  (``omega_max = 255 pi / beta = 16 = 2 x`` the bandwidth);
* ``IterationMax = 1500``, ``Mix = 0.2``, ``EPS = 10``,
  ``mixing_scheme = "anderson"``, ``anderson_depth = 8``.

  **Anderson mixing is mandatory near the CDW.**  With plain linear mixing the
  ``V = 1.2`` points stall on a nearly-marginal mode: the step-to-step
  ``|dSigma|/|Sigma|`` drops below ``1e-8`` while the iterate is still far from
  the fixed point, and the "converged" green differs from the true fixed point
  by ~26 % (measured).  The Anderson solution is stable: ``EPS = 8`` and
  ``EPS = 12`` agree to ``2.7e-9`` relative.

  ``qlms.run()`` returns ``{"scf_converged", "scf_iterations"}`` for FLEX
  (``src/hwave/qlms.py``).  This script asserts ``scf_converged`` is ``True``
  before writing a fixture -- a regeneration that hits
  ``IterationMax = 1500`` unconverged raises loudly instead of silently
  committing a fixture built from a non-converged Green function -- and
  stores both fields in the ``.npz`` metadata.

Post-processing before writing the fixture
------------------------------------------
The stored array is the FLEX output **projected onto the exact
``(k, w) -> (-k, -w)`` conjugation symmetry** of the model,

    ``G <- 0.5 * (G(k, iw) + conj(G(-k, -iw)))``,

which is what makes the static pair weight ``w = GG`` exactly real, as the v1
Hermitian bond path requires (``check_hermitian_preconditions``).  The
projection moves the data by at most ~1e-8 (below the FLEX convergence
tolerance) for every point except a stray run; it is applied unconditionally so
the committed fixture satisfies the precondition by construction.

Layout: the arrays keep the H-wave FLEX layout ``(nblock, nmat, nvol, norb,
norb)`` written by ``flex.save_results`` and are stored with
``np.savez_compressed`` (~2.7x smaller; the test only reads ``green``).

Usage::

    python3 tests/sc/onari_bond/generate_fixtures.py [outdir] [-L 16|32]

``-L`` restricts the run to one grid size (default: every point of ``POINTS``).
"""
import argparse
import os
import shutil
import sys
import tempfile

import numpy as np

U = 4.0
T = 0.02
FILLING = 0.7
NMAT = 256
NN = [(-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0)]

# (L, V) points that the milestone test consumes. The L=16 subset is the
# committed one; the L=32 subset is regenerated on demand (see the module
# docstring).
POINTS = [(16, 0.0), (16, 0.05), (16, 0.4), (16, 0.8), (16, 1.0), (16, 1.2),
          (32, 0.0), (32, 0.8), (32, 1.2)]

# ---------------------------------------------------------------------------
# Task 9 (spec ``2026-07-27-dynamic-bond-channels-design.md`` S3.6): the
# "Phase A milestone" dynamic ``lambda_t(V)`` sweep
# (``tests/test_bond_onari_milestone_dynamic.py``) uses U=4, n=0.7, T=0.02,
# 16x16, Nmat=256, V in {0, 0.4, 0.8, 1.0, 1.2} -- IDENTICAL to the L=16
# subset of ``POINTS`` above (same reduced-FLEX scheme, "Hartree-only V" per
# spec S3.6/S5.3, Anderson depth=8 / Mix=0.2 / IterationMax=1500).
#
# Generating a SECOND, independently self-consistent green at the exact same
# parameters would only reproduce the same physics from a different
# (non-bit-identical) SCF trajectory -- FLEX/Anderson mixing is not
# bit-reproducible, see the module docstring's L=32 discussion. Worse, using
# a DIFFERENT green for the static and the dynamic milestone would confound
# the static/dynamic ratio indicator (spec S3.6) with independent FLEX-SCF
# noise on top of the retardation effect it is meant to measure -- the ratio
# is only an apples-to-apples comparison when both arms consume the SAME
# green. The dynamic milestone therefore intentionally REUSES these exact
# ``green_L16_V*.npz`` fixtures (loaded via
# ``tests.test_bond_onari_milestone._load_green``); this constant documents
# that decision at the point of fixture generation rather than forking a
# parallel ``green_V*.npz`` set that ``generate()`` would have to keep in
# sync with ``POINTS`` by hand.
DYNAMIC_MILESTONE_POINTS = [(16, V) for V in (0.0, 0.4, 0.8, 1.0, 1.2)]
assert set(DYNAMIC_MILESTONE_POINTS) <= set(POINTS), (
    "DYNAMIC_MILESTONE_POINTS must be a subset of the committed L=16 POINTS "
    "it reuses")


def _write_inputs(dirname, V):
    os.makedirs(dirname, exist_ok=True)
    with open(os.path.join(dirname, "geom.dat"), "w") as f:
        f.write("  1.000000000000   0.000000000000   0.000000000000\n")
        f.write("  0.000000000000   1.000000000000   0.000000000000\n")
        f.write("  0.000000000000   0.000000000000   1.000000000000\n")
        f.write("1\n")
        f.write("    0.000000000000000e+00     0.000000000000000e+00"
                "     0.000000000000000e+00\n")
    with open(os.path.join(dirname, "transfer.dat"), "w") as f:
        f.write("Transfer in wannier90-like format for uhfk\n1\n{}\n"
                .format(len(NN)))
        f.write(" ".join(["1"] * len(NN)) + "\n")
        for r in NN:
            f.write("{:4d} {:4d} {:4d}    1    1   1.000000000000"
                    "   0.000000000000\n".format(*r))
    with open(os.path.join(dirname, "coulombintra.dat"), "w") as f:
        f.write("CoulombIntra in wannier90-like format for uhfk\n1\n1\n 1\n")
        f.write("   0    0    0    1    1   {:.12f}   0.000000000000\n"
                .format(U))
    with open(os.path.join(dirname, "coulombinter.dat"), "w") as f:
        f.write("CoulombInter in wannier90-like format for uhfk\n1\n{}\n"
                .format(len(NN)))
        f.write(" ".join(["1"] * len(NN)) + "\n")
        for r in NN:
            f.write("{:4d} {:4d} {:4d}    1    1   {:.12f}"
                    "   0.000000000000\n".format(r[0], r[1], r[2], V))


def _flex_input(indir, outdir, L):
    return {
        "mode": {"mode": "FLEX", "calc_scheme": "reduced",
                 "param": {"T": T, "filling": FILLING,
                           "CellShape": [L, L, 1], "SubShape": [1, 1, 1],
                           "Nmat": NMAT, "IterationMax": 1500,
                           "Mix": 0.2, "EPS": 10,
                           "mixing_scheme": "anderson", "anderson_depth": 8}},
        "file": {"input": {"path_to_input": "",
                           "interaction": {
                               "path_to_input": indir,
                               "Geometry": "geom.dat",
                               "Transfer": "transfer.dat",
                               "CoulombIntra": "coulombintra.dat",
                               "CoulombInter": "coulombinter.dat"}},
                 "output": {"path_to_output": outdir, "green": "green"}},
    }


def _symmetrize(green_raw, L):
    """Project onto G(k, iw) = conj(G(-k, -iw)); in/out H-wave FLEX layout."""
    nblock, nmat, nvol, norb, _ = green_raw.shape
    g = green_raw[0].reshape(nmat, L, L, 1, norb, norb).transpose(
        4, 5, 1, 2, 3, 0)
    inv = np.roll(g[:, :, ::-1, ::-1, ::-1, ::-1], (1, 1, 1), (2, 3, 4))
    delta = float(np.abs(g - np.conj(inv)).max())
    g = 0.5 * (g + np.conj(inv))
    out = g.transpose(5, 2, 3, 4, 0, 1).reshape(
        1, nmat, nvol, norb, norb).copy()
    return out, delta


def fixture_name(L, V):
    return "green_L{}_V{:.2f}.npz".format(L, V)


def generate(outdir, points=None):
    """Write a ``green_L{L}_V{V}.npz`` fixture for every ``(L, V)`` in
    ``points`` (default: all of :data:`POINTS`) into ``outdir``.

    Each point runs an independent FLEX SCF and raises if it does not converge,
    so a caller can trust every file this returns.  Called both from the CLI
    below and from ``test_grid_convergence_16_to_32``, which regenerates the
    uncommitted L = 32 set on demand.
    """
    import hwave.qlms as qlms

    if points is None:
        points = POINTS
    os.makedirs(outdir, exist_ok=True)
    for L, V in points:
        work = tempfile.mkdtemp(prefix="onari_bond_")
        try:
            indir = os.path.join(work, "in")
            _write_inputs(indir, V)
            result = qlms.run(input_dict=_flex_input(indir, work, L))
            scf_converged = bool(result.get("scf_converged", False))
            scf_iterations = int(result.get("scf_iterations", -1))
            assert scf_converged, (
                "FLEX did not converge for L={} V={} after {} iterations "
                "(IterationMax=1500) -- refusing to write a fixture from a "
                "non-converged SCF. Increase IterationMax or investigate "
                "before regenerating.".format(L, V, scf_iterations))
            raw = np.load(os.path.join(work, "green.npz"))["green"]
        finally:
            shutil.rmtree(work, ignore_errors=True)
        green, delta = _symmetrize(raw, L)
        path = os.path.join(outdir, fixture_name(L, V))
        np.savez_compressed(
            path, green=green,
            L=L, V=V, U=U, T=T, filling=FILLING, nmat=NMAT,
            symmetrization_residual=delta,
            scf_converged=scf_converged, scf_iterations=scf_iterations,
            provenance=("FLEX reduced, anderson(depth=8), Mix=0.2, EPS=10, "
                        "IterationMax=1500; t=1 NN square lattice; "
                        "symmetrized under (k,w)->(-k,-w)"))
        print("{}  (symmetrization residual {:.2e}, scf_iterations {})".format(
            path, delta, scf_iterations))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("outdir", nargs="?",
                        default=os.path.dirname(os.path.abspath(__file__)),
                        help="output directory (default: this directory)")
    parser.add_argument("-L", type=int, choices=sorted({L for L, _ in POINTS}),
                        default=None,
                        help="restrict to one grid size (default: all)")
    args = parser.parse_args()
    selected = (None if args.L is None
                else [(L, V) for L, V in POINTS if L == args.L])
    sys.exit(generate(args.outdir, points=selected))
