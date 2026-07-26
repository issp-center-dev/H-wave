"""Task 8 -- physics acceptance milestone (spec S7.7): the bond-resolved
static Eliashberg path makes the f-like TRIPLET eigenvalue ``lambda_t`` RISE
with the inter-site Coulomb ``V`` on the single-band Onari case, while the
collapsed-bond (scalar) path does not follow the same trend.

Reference: Onari, Arita, Kuroki, Aoki, cond-mat/0312314 / PRB 70, 094523
(2004); spec ``docs/superpowers/specs/
2026-07-25-bond-channels-eliashberg-design.md`` S7.7.

Setup (fixed and documented; see ``tests/sc/onari_bond/generate_fixtures.py``)
-----------------------------------------------------------------------------
Single-band square lattice, nearest-neighbour ``t = 1`` only (full C4v: a
one-diagonal ``t'`` would break C4 and split the odd-parity E doublet this test
tracks), ``U = 4``, filling ``n = 0.7``, ``T = 0.02``, static frequency,
``Nmat = 256``.  ``V`` grid ``{0.0, 0.4, 0.8, 1.0, 1.2}`` (all below the
``d -> f`` crossover and the CDW instability), plus ``V = 0.05`` for the
small-``V`` continuity check.

**Green provenance.**  For every ``V`` an INDEPENDENTLY self-consistent FLEX
``green.npz`` (regenerated at that ``V`` -- freezing one green across the sweep
would test only the vertex, while Onari's green is self-consistent).  The
**L = 16** fixtures are committed under ``tests/sc/onari_bond/`` and
content-hash-verified below, so the default CI run never runs FLEX.

**The L = 32 fixtures are NOT committed.**  They are ~3.9 MB and back exactly
one test (``test_grid_convergence_16_to_32``); the repository keeps no other
tracked file above ~150 KB, so committing them would be a precedent change in
permanent history for a single grid-convergence check.  That test is therefore
marked ``slow``, **skipped by default**, and regenerates the L = 32 greens on
demand via ``tests/sc/onari_bond/generate_fixtures.py``.  To run it::

    HWAVE_RUN_SLOW_FIXTURES=1 pytest tests/test_bond_onari_milestone.py \
        -k grid_convergence
    # or select the marker (pytest.ini registers it):
    pytest -m slow tests/test_bond_onari_milestone.py

The regenerated files land in ``tests/sc/onari_bond/_regenerated/`` (override
with ``HWAVE_ONARI_L32_DIR``); the directory is git-ignored and reused, so the
FLEX cost is paid once per checkout.  Regenerated files are metadata- and
symmetry-verified but NOT hash-pinned: neither FLEX nor ``np.savez_compressed``
is bit-reproducible.  The *physics* is still pinned -- the stored
``LAMBDA_BOND_32`` / ``LAMBDA_BASE_32`` tables below are asserted at
``LAMBDA_RTOL`` exactly as before.

**Which code is under test.**  The kernel is built by the PRODUCTION function
``sc._build_bond_operator`` -- the same call ``sc.calc_eliashberg`` makes on the
``bond_channels=true`` branch -- fed with the fixture green.  The green must be
the FLEX one the spec's "externally supplied Green function" calls for: at
``U = 4``, ``beta = 50`` a BARE ``_calc_green`` is already past the RPA spin
instability (measured ``min_q sigma_min/sigma_max = 1.3e-3`` and ``lambda`` up
to 65).  ``[eliashberg] bond_green`` wires exactly that into the TOML entry
point, and ``test_milestone_lambda_is_reproducible_from_the_toml_entry_point``
below closes the loop: a full ``calc_eliashberg`` run with
``bond_green = <fixture>`` reproduces the ``lambda_t`` of the direct call.

Baseline arm
------------
The "current scalar path" comparison uses the **2-index
``_compute_vertices_simple`` path** (``Wc = U + 2V(q)``, ``Ws = -U``) -- the
genuine collapsed-bond approximation.  It is NOT the general/4-index path:
at ``norb = 1`` ``_build_sc_matrices_all_q`` drops ``V`` entirely (its Cases 2
and 3 require ``a != b``), so a 4-index baseline would be V-independent by
construction and prove nothing.  Both arms consume the SAME fixture green and
the SAME bubble (the baseline's ``chi0q`` is the ``m = m' = 0`` block of the
bond ``chi_bar``), so the only difference between them is the bond structure of
the vertex.

Measured outcome (see the stored reference tables below)
--------------------------------------------------------
* the tracked f-like triplet ``lambda_t`` rises 0.0705 -> 0.3101 (16x16) and
  0.0739 -> 0.3358 (32x32): monotone, ``Delta lambda = 0.24`` >> the required
  0.02, ratio 4.4;
* the rise is entirely in ``lambda_t^fl`` (``lambda_t^pp <= 0`` at every V, as
  spec S8 predicts for the reversal-complete bare pp vertex), and 74 % of it is
  carried by the CHARGE fluctuations -- Onari's mechanism;
* **the scalar baseline is NOT flat** at these parameters.  Spec S7.7 expected
  ``<= 5 %``; the measured static-path spread is 126 %, because the collapsed
  ``Delta r = 0`` density channel still carries the ``-1/2 Wc chi_c Wc`` term
  and ``chi_c`` grows toward the CDW.  (The flat behaviour recorded in the
  project's earlier benchmark was the DYNAMIC path, where retardation
  suppresses it.)  Per spec S7.7 this comparison is "a stored baselined
  observation, not an invariant", so the observation is stored and pinned
  below instead of being asserted as flatness.  What IS asserted is the
  discriminating statement: the bond path is above the scalar path at every
  ``V > 0`` and the gap grows monotonically to 0.15.
"""
import hashlib
import os

import numpy as np
import pytest

import hwave.sc as sc
from hwave.solver import bond_channels as bc


# ---------------------------------------------------------------------------
# fixtures: committed FLEX green.npz, content-hash verified (spec S7.7
# "Reproducibility record"). The hash pins the COMMITTED BYTES; regenerating
# with tests/sc/onari_bond/generate_fixtures.py reproduces the physics, not the
# hash (neither FLEX nor np.savez_compressed is bit-reproducible).
#
# Only the L=16 set is committed. L=32 is regenerated on demand into
# L32_CACHE_DIR by the slow, skipped-by-default grid-convergence test (see the
# module docstring), so it carries no hash -- its metadata and its physics are
# verified instead.
# ---------------------------------------------------------------------------

FIXTURE_DIR = os.path.join("tests", "sc", "onari_bond")

FIXTURE_SHA256 = {
    (16, 0.0):  "18e868b6df51e1a84064de8c38cdf652d7a04eacbf69a6d520e6b472e04d0b33",
    (16, 0.05): "f773d472bac5fc0bfd23b7c6bd822f1e11f7c9c70f8c52a61d942a5183439bec",
    (16, 0.4):  "c5911512eaad8f43b529b1d38f5700514b14c5a483bfc6bedbf728e8a69a41cf",
    (16, 0.8):  "15b035404699c1c5fa04ae64712aa4fd0f0cd83b1c40d45f212a3f1722567b71",
    (16, 1.0):  "3638458e4fa805fad8c5ef53dd848da3767b434dd15b2badc794e2d29160e0db",
    (16, 1.2):  "6df1b5ef9e26d7eb4ba83c0c1e05785ff278d7cec9fc16bd74f1ff90e60accff",
}

# Where the uncommitted L=32 greens are regenerated to / cached.
L32_CACHE_DIR = os.environ.get(
    "HWAVE_ONARI_L32_DIR", os.path.join(FIXTURE_DIR, "_regenerated"))
# Env-var opt-in for the slow regeneration (the marker `-m slow` also enables
# it; see _require_slow_fixtures).
_SLOW_ENV = "HWAVE_RUN_SLOW_FIXTURES"

U = 4.0
T = 0.02
BETA = 1.0 / T
NMAT = 256
FILLING = 0.7
V_GRID = (0.0, 0.4, 0.8, 1.0, 1.2)
V_GRID_32 = (0.0, 0.8, 1.2)
NN_BONDS = ((-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0))

# --- acceptance thresholds (spec S7.7; NOT tunable knobs) ------------------
SIGMA_FLOOR = 1.0e-3        # off-instability guard, complex128
ATOL_MONO = 1.0e-6          # monotonicity tolerance
DELTA_LAMBDA_MIN = 0.02     # minimum absolute rise
RATIO_MIN = 1.3             # minimum ratio (baseline lambda_t(0) > 0.02)
DEG_TOL = 1.0e-3            # degenerate-cluster tolerance

# --- stored reference numbers (this branch, CPU/numpy) ---------------------
# tracked f-like triplet lambda of the BOND path
LAMBDA_BOND_16 = {0.0: 0.0704950, 0.4: 0.0925329, 0.8: 0.1218159,
                  1.0: 0.1644559, 1.2: 0.3101238}
LAMBDA_BOND_32 = {0.0: 0.0739082, 0.8: 0.1570313, 1.2: 0.3358213}
# ... and of the 2-index collapsed-bond SCALAR baseline (stored observation)
LAMBDA_BASE_16 = {0.0: 0.0704950, 0.4: 0.0699692, 0.8: 0.0637392,
                  1.0: 0.0841373, 1.2: 0.1595172}
LAMBDA_BASE_32 = {0.0: 0.0739082, 0.8: 0.0722234, 1.2: 0.1566309}
LAMBDA_RTOL = 1.0e-5


def _fixture_path(L, V):
    """Path of a fixture green.

    L=16 is committed next to ``generate_fixtures.py``; L=32 lives in the
    git-ignored regeneration cache (see the module docstring).
    """
    directory = FIXTURE_DIR if (L, V) in FIXTURE_SHA256 else L32_CACHE_DIR
    return os.path.join(directory, "green_L{}_V{:.2f}.npz".format(L, V))


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _generator_module():
    """Import ``tests/sc/onari_bond/generate_fixtures.py`` as a module.

    The generator is the single source of truth for the model AND for the
    fixture regeneration, so both the TOML-entry-point test and the L=32
    regeneration below go through it.
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "onari_bond_generate_fixtures",
        os.path.join(FIXTURE_DIR, "generate_fixtures.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _require_slow_fixtures(request):
    """Skip unless the caller explicitly asked for the slow path.

    Enabled either by ``HWAVE_RUN_SLOW_FIXTURES=1`` or by selecting the marker
    (``pytest -m slow``). Everything else -- including a plain ``pytest`` and
    an explicit ``-m "not slow"`` -- skips.
    """
    if os.environ.get(_SLOW_ENV, "").lower() not in ("", "0", "false", "no"):
        return
    markexpr = str(request.config.getoption("markexpr", default="") or "")
    tokens = markexpr.replace("(", " ").replace(")", " ").split()
    if "slow" in tokens and "not" not in tokens:
        return
    pytest.skip(
        "L=32 grid-convergence test is skipped by default: its ~3.9 MB FLEX "
        "green fixtures are deliberately NOT committed and must be "
        "regenerated (minutes of FLEX). Run it with {}=1 or with "
        "`pytest -m slow`; see the module docstring.".format(_SLOW_ENV))


def _ensure_l32_fixtures():
    """Regenerate whichever L=32 greens are missing from the cache dir."""
    missing = [(32, V) for V in V_GRID_32
               if not os.path.exists(_fixture_path(32, V))]
    if not missing:
        return
    print("\nregenerating {} uncommitted L=32 FLEX fixture(s) into {} "
          "(this takes minutes)".format(len(missing), L32_CACHE_DIR))
    _generator_module().generate(L32_CACHE_DIR, points=missing)
    for L, V in missing:
        assert os.path.exists(_fixture_path(L, V)), _fixture_path(L, V)


def _check_fixture_metadata(path, L, V):
    """The recorded setup + the (k,w) -> (-k,-w) symmetry of one fixture.

    Everything ``test_fixtures_are_intact_and_carry_the_documented_setup``
    asserts EXCEPT the content hash, so it applies equally to the committed
    L=16 files and to the regenerated (unhashable) L=32 ones.
    """
    data = np.load(path)
    assert int(data["L"]) == L
    assert float(data["V"]) == pytest.approx(V)
    assert float(data["U"]) == U
    assert float(data["T"]) == T
    assert float(data["filling"]) == FILLING
    assert int(data["nmat"]) == NMAT
    assert data["green"].shape == (1, NMAT, L * L, 1, 1)

    nmat = data["green"].shape[1]
    green = data["green"][0].reshape(nmat, L, L, 1, 1, 1).transpose(
        4, 5, 1, 2, 3, 0)
    inv = np.roll(green[:, :, ::-1, ::-1, ::-1, ::-1], (1, 1, 1), (2, 3, 4))
    assert np.abs(green - np.conj(inv)).max() < 1e-13

    assert float(data["symmetrization_residual"]) < 1e-8
    provenance = str(data["provenance"])
    assert provenance, "fixture {} has empty provenance".format(path)
    if "scf_converged" in data.files:
        assert bool(data["scf_converged"]) is True, (
            "fixture {} was written from a non-converged FLEX SCF "
            "(scf_iterations={})".format(path, int(data["scf_iterations"])))


_GREEN_CACHE = {}


def _load_green(L, V):
    """Fixture green in sc.py layout (1,1,L,L,1,NMAT).

    Committed (L=16) fixtures are hash-verified against ``FIXTURE_SHA256``.
    Regenerated (L=32) fixtures cannot be -- FLEX/savez are not bit-
    reproducible -- so their recorded metadata is verified instead.
    """
    key = (L, V)
    if key in _GREEN_CACHE:
        return _GREEN_CACHE[key]
    path = _fixture_path(L, V)
    if key in FIXTURE_SHA256:
        digest = _sha256(path)
        assert digest == FIXTURE_SHA256[key], (
            "fixture {} has hash {} but {} was recorded -- the committed green "
            "changed; regenerate the reference numbers with "
            "tests/sc/onari_bond/generate_fixtures.py".format(
                path, digest, FIXTURE_SHA256[key]))
    else:
        assert os.path.exists(path), (
            "uncommitted fixture {} is missing; it is regenerated on demand by "
            "_ensure_l32_fixtures()".format(path))
        _check_fixture_metadata(path, L, V)
    raw = np.load(path)["green"]
    nmat = raw.shape[1]
    green = raw[0].reshape(nmat, L, L, 1, 1, 1).transpose(
        4, 5, 1, 2, 3, 0).copy()
    _GREEN_CACHE[key] = green
    return green


# ---------------------------------------------------------------------------
# model / solve helpers
# ---------------------------------------------------------------------------

def _kgrid(L):
    return (np.linspace(0, 2.0 * np.pi, L, endpoint=False),
            np.linspace(0, 2.0 * np.pi, L, endpoint=False),
            np.linspace(0, 2.0 * np.pi, 1, endpoint=False))


def _interactions(V, declare_inter=True):
    it = {"CoulombIntra": {((0, 0, 0), (0, 0)): U}}
    if declare_inter:
        it["CoulombInter"] = {(r, (0, 0)): V for r in NN_BONDS}
    return it


def _odd_basis_matrix(L):
    """Orthonormal basis of the odd-parity (``phi(-k) = -phi(k)``) subspace.

    Columns of the returned ``(L*L, d)`` matrix.  The four self-inverse points
    (Gamma, (pi,0), (0,pi), (pi,pi)) carry no odd component and are excluded,
    so ``d = (L*L - 4) / 2``.
    """
    n = L * L
    idx = np.arange(n).reshape(L, L)
    neg = np.roll(idx[::-1, ::-1], (1, 1), (0, 1))
    cols, seen = [], set()
    for i in range(L):
        for j in range(L):
            a, b = int(idx[i, j]), int(neg[i, j])
            if a == b or (a, b) in seen or (b, a) in seen:
                continue
            seen.add((a, b))
            v = np.zeros(n)
            v[a] = 1.0 / np.sqrt(2.0)
            v[b] = -1.0 / np.sqrt(2.0)
            cols.append(v)
    return np.array(cols).T


def _odd_spectrum(A, n, W, L, n_ev=8):
    """FULL odd-parity spectrum of the kernel, exactly (no ARPACK window).

    ``K = -Gamma W`` is symmetrizable by the pair weight: with the (single
    band, hence per-k scalar) metric ``w`` the matrix
    ``Ktilde = sqrt(w) K sqrt(w)^-1`` is Hermitian (this is what
    ``check_hermitian_preconditions`` verifies), so a dense ``eigh`` restricted
    to the odd subspace gives every odd eigenvalue with no convergence
    parameters and no risk of ARPACK returning an even-parity mode.  ``n`` is
    only 256 (L=16) / 1024 (L=32), so the dense build is cheap.

    Returns ``(evals descending, eigenvectors as ROWS in the physical gap
    basis)`` -- ROWS because that is what ``bc.track_subspace`` consumes
    (``scipy.sparse.linalg.eigs`` would return COLUMNS).
    """
    K = np.empty((n, n), dtype=complex)
    e = np.zeros(n, dtype=complex)
    for i in range(n):
        e[i] = 1.0
        K[:, i] = A.matvec(e)
        e[i] = 0.0
    w = np.real(np.asarray(W).reshape(L * L))
    s = np.sqrt(w)
    Ktilde = (s[:, None] * K) / s[None, :]
    B = _odd_basis_matrix(L)
    Kodd = B.T @ Ktilde @ B
    Kodd = 0.5 * (Kodd + Kodd.conj().T)
    ev, evec = np.linalg.eigh(Kodd)
    order = np.argsort(-ev)
    ev, evec = ev[order], evec[:, order]
    phys = (B @ evec) / s[:, None]
    return ev[:n_ev], phys[:, :n_ev].T.copy()


def _conditioning(chi_bar, S_bond, C_bond):
    """Off-instability guard (spec S7.7): smallest singular value of the two
    ACTUAL solve matrices, and its ratio to the largest, minimized over q."""
    ND = S_bond.shape[-1]
    eye = np.eye(ND, dtype=complex)
    out = {}
    for name, M in (("spin", eye - chi_bar @ S_bond),
                    ("charge", eye + chi_bar @ C_bond)):
        sv = np.linalg.svd(M.reshape(-1, ND, ND), compute_uv=False)
        out[name] = {"sigma_min": float(sv[:, -1].min()),
                     "ratio": float((sv[:, -1] / sv[:, 0]).min())}
    return out


_POINT_CACHE = {}


def _bond_point(L, V, declare_inter=True):
    """One bond-path V point, built by the PRODUCTION ``_build_bond_operator``."""
    key = ("bond", L, V, declare_inter)
    if key in _POINT_CACHE:
        return _POINT_CACHE[key]
    kx, ky, kz = _kgrid(L)
    it = _interactions(V, declare_inter)
    inter_k = sc._build_interaction_k(kx, ky, kz, it, 1)
    bond_set = bc.resolve_interactions(it.get("CoulombInter", {}),
                                       np.eye(3), 1)
    green = _load_green(L, V)

    (A, vec_size), provenance, attribution = sc._build_bond_operator(
        bond_set, green, it, inter_k, {"norb": 1, "rvec": np.eye(3)}, 1,
        kx, ky, kz, BETA, "triplet",
        bond_max_shells=None, bond_memory_cap_gb=0.0)

    # the guard quantities (chi_bar / S / C) are diagnostics, rebuilt from the
    # same public chain the operator used
    chi_bar = bc.bond_bubble(green, bond_set, BETA)
    S0, C0 = sc._build_bond_m0_blocks(bond_set, it, inter_k, 1, kx, ky, kz)
    S_bond, C_bond, _, _ = bc.bare_bond_vertices(bond_set, S0, C0, 1)

    evals, vecs = _odd_spectrum(A, vec_size, attribution["weight"], L)
    out = {"evals": evals, "vecs": vecs, "weight": attribution["weight"],
           "A": A, "A_pp": attribution["pp"], "A_fl": attribution["fl"],
           "conditioning": _conditioning(chi_bar, S_bond, C_bond),
           "n_channels": bond_set.n_channels, "provenance": provenance,
           "diagnostics": attribution["diagnostics"]}
    _POINT_CACHE[key] = out
    return out


def _baseline_point(L, V):
    """One 2-index ``_compute_vertices_simple`` (collapsed-bond) V point.

    ``chi0q`` is the ``m = m' = 0`` block of the bond bubble, i.e. the ordinary
    bare susceptibility of the SAME green -- so the two arms differ only in the
    vertex, never in the input data.
    """
    key = ("base", L, V)
    if key in _POINT_CACHE:
        return _POINT_CACHE[key]
    kx, ky, kz = _kgrid(L)
    it = _interactions(V)
    inter_k = sc._build_interaction_k(kx, ky, kz, it, 1)
    bond_set = bc.resolve_interactions(it["CoulombInter"], np.eye(3), 1)
    green = _load_green(L, V)
    chi_bar = bc.bond_bubble(green, bond_set, BETA)
    chi0q = chi_bar[:, :, :, 0, 0].reshape(1, 1, L, L, 1, 1).copy()

    vertex = sc._compute_vertices(chi0q, inter_k, 1, L, L, 1, 1,
                                  pairing_type="triplet", static_index=0)
    assert isinstance(vertex, tuple), (
        "the baseline must be the 2-index simple path (Pc, Ps); a 4-index "
        "vertex means _compute_vertices took the general branch, which drops "
        "V entirely at norb=1")
    Pc_q, Ps_q = vertex
    G2 = sc._calc_g2(green, BETA)
    A, vec_size = sc._make_kernel_operator(Pc_q + Ps_q, G2, 1, L, L, 1)
    W = bc.pair_weight(green, BETA)
    evals, vecs = _odd_spectrum(A, vec_size, W, L)
    out = {"evals": evals, "vecs": vecs, "weight": W, "A": A}
    _POINT_CACHE[key] = out
    return out


_SWEEP_CACHE = {}


def _sweep(L, v_grid, arm):
    """Invariant-subspace-tracked sweep (spec S7.7 (i)-(v))."""
    key = (L, v_grid, arm)
    if key in _SWEEP_CACHE:
        return _SWEEP_CACHE[key]
    kx, ky, kz = _kgrid(L)
    seed = sc._initialize_gap("f_x", 1, kx, ky, kz).ravel()
    basis = sc._odd_harmonic_basis(1, kx, ky, kz)
    point = _bond_point if arm == "bond" else _baseline_point
    points = [point(L, V) for V in v_grid]
    # capture_tol: _odd_spectrum returns the COMPLETE odd spectrum, so the
    # "widen the eigenvalue window" escalation capture_tol gates is vacuous
    # here; it is lowered only so the anchor point (whose overlap with the bare
    # f_x form factor is 0.21 -- the physical state carries higher odd
    # harmonics too) is not flagged. The anchoring itself is unambiguous: the
    # leading doublet's f-seed overlap beats the runner-up's by ~6x
    # (test_anchor_cluster_is_the_most_f_like).
    records = bc.track_subspace(
        [(p["evals"], p["vecs"]) for p in points], seed=seed,
        weight=points[0]["weight"], basis=basis, deg_tol=DEG_TOL,
        capture_tol=0.2)
    _SWEEP_CACHE[key] = (points, records)
    return points, records


def _attribution(point, record):
    """pp/fl split of the tracked subspace's leading eigenvector."""
    lead = record["cluster"][int(np.argmax(np.real(record["eigenvalues"])))]
    return bc.attribute_lambda(point["vecs"][lead], point["weight"],
                               point["A_pp"], point["A_fl"],
                               op_full=point["A"])


# ===========================================================================
# tests
# ===========================================================================

def test_fixtures_are_intact_and_carry_the_documented_setup():
    """Every committed green.npz matches its recorded hash and metadata, and
    satisfies the (k,w) -> (-k,-w) symmetry the v1 Hermitian path needs.

    Spec S7.7 "convergence status" / "Reproducibility record": the fixture
    must also carry a tight symmetrization residual, a non-empty provenance
    string, and (I-1) an explicit converged-FLEX flag.  Fixtures regenerated
    before I-1 landed do not have the ``scf_converged``/``scf_iterations``
    keys yet -- that is tolerated here (no regeneration of the green data
    itself is required for the metadata fix), but any fixture that DOES carry
    the key must say converged=True.

    Only the COMMITTED (L=16) set is covered here; the regenerated L=32 set has
    no committed bytes to hash, and ``_load_green`` runs the same metadata
    checks on it via ``_check_fixture_metadata``.
    """
    assert {L for L, _ in FIXTURE_SHA256} == {16}, (
        "only the L=16 fixtures are committed; an L=32 hash entry means the "
        "3.9 MB set was re-added to git (see the module docstring)")
    for (L, V), digest in FIXTURE_SHA256.items():
        path = _fixture_path(L, V)
        assert os.path.exists(path), path
        assert _sha256(path) == digest, path
        _check_fixture_metadata(path, L, V)


@pytest.mark.parametrize("V", V_GRID)
def test_off_instability_guard_at_every_V(V):
    """Spec S7.7: ``min_q sigma_min`` (and its ratio to ``sigma_max``) of the
    ACTUAL solve matrices ``I - chi_bar S`` and ``I + chi_bar C`` must stay
    above the floor at EVERY V point -- a point that fails FAILS the run, it is
    not dropped from the metrics."""
    cond = _bond_point(16, V)["conditioning"]
    for channel in ("spin", "charge"):
        assert cond[channel]["ratio"] > SIGMA_FLOOR, (
            "V={} {} channel: min_q sigma_min/sigma_max = {:.3e} <= {} -- the "
            "V grid entered the unstable region".format(
                V, channel, cond[channel]["ratio"], SIGMA_FLOOR))
        assert cond[channel]["sigma_min"] > SIGMA_FLOOR


def test_tracked_f_triplet_lambda_rises_with_V():
    """THE MILESTONE (spec S7.7 acceptance): the tracked f-like triplet
    eigenvalue is monotone non-decreasing in V and rises by well over
    ``Delta lambda_min`` across the grid."""
    points, records = _sweep(16, V_GRID, "bond")
    lam = [r["lambda"] for r in records]

    for V, value in zip(V_GRID, lam):
        assert value == pytest.approx(LAMBDA_BOND_16[V], rel=LAMBDA_RTOL)

    for i in range(1, len(lam)):
        assert lam[i] >= lam[i - 1] - ATOL_MONO, (
            "lambda_t is not monotone at V={}: {} -> {}".format(
                V_GRID[i], lam[i - 1], lam[i]))

    assert lam[-1] - lam[0] >= DELTA_LAMBDA_MIN
    assert lam[0] > 0.02                      # positive baseline
    assert lam[-1] / lam[0] >= RATIO_MIN

    # every point is a genuinely tracked, doubly degenerate odd (E) doublet
    for rec in records:
        assert rec["dim"] == 2
    for rec in records[1:]:
        assert rec["overlap"] >= 0.9


def test_the_rise_lives_in_the_fluctuation_part():
    """Spec S4.5/S7.10b: ``lambda = lambda^pp + lambda^fl`` exactly, the bare
    (reversal-complete) particle-particle vertex is REPULSIVE
    (``lambda^pp <= 0``), and the whole V-rise is carried by ``lambda^fl``."""
    points, records = _sweep(16, V_GRID, "bond")
    attrs = [_attribution(p, r) for p, r in zip(points, records)]

    for V, rec, at in zip(V_GRID, records, attrs):
        assert at["sum_residual"] < 1e-10
        assert at["lambda"] == pytest.approx(rec["lambda"], rel=1e-8)
        assert at["lambda_pp"] <= 1e-12, (
            "V={}: lambda^pp = {:.3e} > 0 -- the bare pp vertex must be "
            "repulsive in the triplet channel".format(V, at["lambda_pp"]))
        assert abs(at["imag"]) < 1e-8 * max(1.0, abs(at["lambda"]))

    fl = [a["lambda_fl"] for a in attrs]
    assert fl[-1] - fl[0] >= DELTA_LAMBDA_MIN
    for i in range(1, len(fl)):
        assert fl[i] >= fl[i - 1] - ATOL_MONO


def test_charge_fluctuations_carry_most_of_the_rise():
    """Onari's mechanism: the fluctuation kernel splits EXACTLY into its spin
    and charge halves, and it is the CHARGE (bond) fluctuations that grow with
    V and dominate the rise."""
    kx, ky, kz = _kgrid(16)
    _, records = _sweep(16, V_GRID, "bond")
    spin_part, charge_part = [], []
    for V, rec in zip(V_GRID, records):
        it = _interactions(V)
        inter_k = sc._build_interaction_k(kx, ky, kz, it, 1)
        bond_set = bc.resolve_interactions(it["CoulombInter"], np.eye(3), 1)
        green = _load_green(16, V)
        chi_bar = bc.bond_bubble(green, bond_set, BETA)
        S0, C0 = sc._build_bond_m0_blocks(bond_set, it, inter_k, 1,
                                          kx, ky, kz)
        S_b, C_b, Vpp_s, Vpp_t = bc.bare_bond_vertices(bond_set, S0, C0, 1)
        chi_s, chi_c = bc.dress_bond(chi_bar, S_b, C_b)
        zero = np.zeros_like(chi_s)
        args = (S_b, C_b, Vpp_s, Vpp_t, green, bond_set, "triplet", BETA)
        _, fl_all, _, _ = bc.make_bond_kernel_parts(chi_s, chi_c, *args)
        _, fl_spin, _, _ = bc.make_bond_kernel_parts(chi_s, zero, *args)
        _, fl_chg, _, _ = bc.make_bond_kernel_parts(zero, chi_c, *args)
        _, fl_none, _, _ = bc.make_bond_kernel_parts(zero, zero, *args)

        point = _bond_point(16, V)
        lead = rec["cluster"][int(np.argmax(np.real(rec["eigenvalues"])))]
        v, W = point["vecs"][lead], point["weight"]

        def rq(op):
            return bc.attribute_lambda(v, W, op, op, op_full=op)["lambda"]

        base = rq(fl_none)
        assert base == pytest.approx(0.0, abs=1e-12)     # no chi => no fl
        s, c, tot = rq(fl_spin) - base, rq(fl_chg) - base, rq(fl_all)
        assert s + c + base == pytest.approx(tot, abs=1e-10)
        spin_part.append(s)
        charge_part.append(c)

    # the charge channel is what switches on with V
    assert charge_part[-1] - charge_part[0] >= DELTA_LAMBDA_MIN
    assert charge_part[-1] > 10.0 * charge_part[0]
    assert (charge_part[-1] - charge_part[0]) > (spin_part[-1] - spin_part[0])


def test_tracked_state_is_f_like_and_not_p_like():
    """Spec S7.7: report/assert the harmonic decomposition of the TRACKED
    subspace on the odd basis.  No point-group separation is claimed -- p- and
    f-like harmonics span the same odd 2D rep -- but the tracked state must be
    overwhelmingly f-like, and its f content must GROW with V."""
    _, records = _sweep(16, V_GRID, "bond")
    f_content, p_content = [], []
    for rec in records:
        h = rec["harmonics"]
        f = h["f_x"] + h["f_y"]
        p = h["p_x"] + h["p_y"]
        assert 0.0 <= f <= rec["dim"] + 1e-8
        f_content.append(f)
        p_content.append(p)
    # the absolute f-capture at V_max (measured 1.36115, growing from 0.42395
    # at V=0): a real pin, unlike "f > 10*max(p, 1e-12)" which the p floor
    # made pass with a ~4000x margin regardless of the actual f value. This is
    # a DIAGNOSTIC pin (not the headline lambda claim, which stays at
    # LAMBDA_RTOL = 1e-5 above), on a quantity derived from eigenvector
    # overlaps in a near-degenerate subspace -- far more sensitive to a
    # BLAS/numpy/ARPACK version bump than a bare eigenvalue. Loosened from
    # rel=1e-4 (review fix: too brittle for a diagnostic) to rel=1e-2, still
    # tight enough to catch a real regression (the qualitative checks below
    # already cover "still f-like, still growing").
    assert f_content[-1] == pytest.approx(1.361154, rel=1e-2)
    assert f_content[-1] > 2.0 * f_content[0]
    assert max(p_content) < 0.01


def test_anchor_cluster_is_the_most_f_like_one():
    """The V=0 anchor (spec S7.7 (v)) is unambiguous: among ALL near-degenerate
    clusters the tracked one has by far the largest f-seed overlap."""
    kx, ky, kz = _kgrid(16)
    seed = sc._initialize_gap("f_x", 1, kx, ky, kz).ravel()[np.newaxis, :]
    point = _bond_point(16, 0.0)
    clusters = bc.cluster_eigenvalues(point["evals"], deg_tol=DEG_TOL)
    scores = [bc.subspace_similarity(point["vecs"][c], seed,
                                     point["weight"])[0] for c in clusters]
    assert np.argmax(scores) == 0                  # the LEADING cluster
    assert scores[0] > 4.0 * max(scores[1:])


def test_bond_path_reduces_to_the_scalar_path_at_V_zero():
    """Spec S7.4: at V=0 (topology declared, values zero -> B=5) the bond
    kernel and the collapsed-bond scalar kernel give the same lambda."""
    _, bond_rec = _sweep(16, V_GRID, "bond")
    _, base_rec = _sweep(16, V_GRID, "base")
    assert bond_rec[0]["lambda"] == pytest.approx(base_rec[0]["lambda"],
                                                  abs=1e-12)


def test_scalar_baseline_is_a_stored_observation_not_flat():
    """Spec S7.7 baseline arm -- the 2-index ``_compute_vertices_simple``
    (collapsed-bond) path on the SAME greens.

    S7.7 anticipated a flat (<= 5 %) scalar baseline.  That expectation does
    NOT hold for the STATIC path at these parameters: the retained Delta r = 0
    density channel still carries -1/2 Wc chi_c Wc, and chi_c grows toward the
    CDW, so the baseline spread is ~126 %.  Per S7.7 this comparison is "a
    stored baselined observation, not an invariant", so the observed numbers
    are PINNED here (a regression guard) and the discriminating statement --
    bond above scalar at every V > 0, with a monotonically growing gap -- is
    what is asserted."""
    _, bond_rec = _sweep(16, V_GRID, "bond")
    _, base_rec = _sweep(16, V_GRID, "base")
    lam_bond = [r["lambda"] for r in bond_rec]
    lam_base = [r["lambda"] for r in base_rec]

    for V, value in zip(V_GRID, lam_base):
        assert value == pytest.approx(LAMBDA_BASE_16[V], rel=LAMBDA_RTOL)

    # the recorded (non-flat) spread of the scalar path. A DIAGNOSTIC pin (the
    # per-V lam_base values above are already pinned individually at the
    # tight LAMBDA_RTOL); loosened from rel=1e-3 (review fix: too brittle for
    # a diagnostic derived quantity) to rel=1e-2.
    spread = max(abs(x - lam_base[0]) for x in lam_base) / lam_base[0]
    assert spread == pytest.approx(1.2628, rel=1e-2)

    gap = [b - s for b, s in zip(lam_bond, lam_base)]
    assert gap[0] == pytest.approx(0.0, abs=1e-12)
    for i in range(1, len(gap)):
        assert gap[i] > 0.0
        assert gap[i] >= gap[i - 1] - ATOL_MONO
    # measured gap[-1] = 0.1506; 0.10 leaves a defensible margin below the
    # exact value pinned by the rel=1e-5 lambda approxes above, rather than
    # the previous 0.15 threshold that was only a 0.4% margin from failing.
    assert gap[-1] >= 0.10


@pytest.mark.slow
def test_grid_convergence_16_to_32(request):
    """Spec S7.7: the 16x16 trend must be REPRODUCED on 32x32 -- same rising,
    monotone, f-like doublet, and the same qualitative bond-vs-scalar contrast.
    A discrepant trend is a failure, not a pass.

    SKIPPED BY DEFAULT.  The L=32 greens are the only fixtures this repository
    would have to commit above ~150 KB (3.9 MB for this one test), so they are
    regenerated on demand instead -- run with ``HWAVE_RUN_SLOW_FIXTURES=1`` or
    ``pytest -m slow`` (module docstring has the details).  The stored
    ``LAMBDA_*_32`` tables are asserted at the unchanged ``LAMBDA_RTOL``, so
    the physics claim is exactly the one it always was; only the *bytes* of the
    green are no longer pinned.
    """
    _require_slow_fixtures(request)
    _ensure_l32_fixtures()
    points32, rec32 = _sweep(32, V_GRID_32, "bond")
    lam32 = [r["lambda"] for r in rec32]
    for V, value in zip(V_GRID_32, lam32):
        assert value == pytest.approx(LAMBDA_BOND_32[V], rel=LAMBDA_RTOL)

    for i in range(1, len(lam32)):
        assert lam32[i] >= lam32[i - 1] - ATOL_MONO
    assert lam32[-1] - lam32[0] >= DELTA_LAMBDA_MIN
    assert lam32[-1] / lam32[0] >= RATIO_MIN

    for point, rec in zip(points32, rec32):
        assert rec["dim"] == 2
        h = rec["harmonics"]
        assert h["f_x"] + h["f_y"] > 10.0 * (h["p_x"] + h["p_y"])
        assert point["conditioning"]["spin"]["ratio"] > SIGMA_FLOOR
        assert point["conditioning"]["charge"]["ratio"] > SIGMA_FLOOR
        at = _attribution(point, rec)
        assert at["lambda_pp"] <= 1e-12

    # same contrast against the collapsed-bond scalar path
    _, base32 = _sweep(32, V_GRID_32, "base")
    lam_base32 = [r["lambda"] for r in base32]
    for V, value in zip(V_GRID_32, lam_base32):
        assert value == pytest.approx(LAMBDA_BASE_32[V], rel=LAMBDA_RTOL)
    assert lam32[0] == pytest.approx(lam_base32[0], abs=1e-12)
    assert lam32[-1] - lam_base32[-1] >= 0.15

    # the two grids agree on the size of the effect to better than 20 %
    _, rec16 = _sweep(16, V_GRID, "bond")
    d16 = rec16[-1]["lambda"] - rec16[0]["lambda"]
    d32 = lam32[-1] - lam32[0]
    assert abs(d32 - d16) / d16 < 0.2


def test_small_V_continuity_across_the_B1_to_B5_topology_change():
    """Spec S7.7: the ``B = 1 -> B = 5`` topology change introduces no
    discontinuity -- V=0 with no CoulombInter declared (B=1) reproduces V=0
    with a declared-but-zero NN topology (B=5) exactly, and V=0.05 (B=5) is
    continuously close to it."""
    p_b1 = _bond_point(16, 0.0, declare_inter=False)
    p_b5 = _bond_point(16, 0.0, declare_inter=True)
    p_005 = _bond_point(16, 0.05, declare_inter=True)
    assert p_b1["n_channels"] == 1
    assert p_b5["n_channels"] == 5
    assert p_005["n_channels"] == 5
    assert p_b1["provenance"].get("collapsed_to_pure_hubbard") is True

    kx, ky, kz = _kgrid(16)
    seed = sc._initialize_gap("f_x", 1, kx, ky, kz).ravel()
    basis = sc._odd_harmonic_basis(1, kx, ky, kz)
    records = bc.track_subspace(
        [(p["evals"], p["vecs"]) for p in (p_b1, p_b5, p_005)],
        seed=seed, weight=p_b1["weight"], basis=basis, deg_tol=DEG_TOL,
        capture_tol=0.2)
    lam_b1, lam_b5, lam_005 = [r["lambda"] for r in records]

    assert lam_b1 == pytest.approx(lam_b5, abs=1e-12)
    # No jump at the topology change: the V = 0.05 response must not exceed the
    # linear extrapolation of the first (smooth, B=5) grid step,
    # (0.05/0.4) * [lambda(0.4) - lambda(0)] = 2.76e-3. Measured: 2.45e-3, i.e.
    # slightly SUB-linear -- the branch is continuous and smooth through
    # B = 1 -> B = 5.
    step = (0.05 / 0.4) * (LAMBDA_BOND_16[0.4] - LAMBDA_BOND_16[0.0])
    assert 0.0 < lam_005 - lam_b5 < step
    assert records[2]["overlap"] > 0.9


def test_runtime_preconditions_hold_at_every_V():
    """The v1 Hermitian-path preconditions (spec S4.5) are satisfied by every
    acceptance point -- the reported lambda really is a real eigenvalue."""
    for V in V_GRID:
        diag = _bond_point(16, V)["diagnostics"]
        assert diag["weight_min_eigenvalue"] > 0.0
        assert diag["weight_hermiticity_residual"] < 1e-10
        assert diag["kernel_hermiticity_relative"] < 1e-8


# ---------------------------------------------------------------------------
# the milestone through the TOML entry point ([eliashberg] bond_green)
# ---------------------------------------------------------------------------

def _write_model_inputs(dirname, V):
    """The exact model files ``generate_fixtures.py`` fed to FLEX.

    Imported from the generator itself so the TOML run below cannot drift away
    from the model the committed green was produced with.
    """
    mod = _generator_module()
    assert mod.U == U and mod.T == T and mod.NMAT == NMAT
    assert mod.FILLING == FILLING
    mod._write_inputs(dirname, V)


def test_milestone_lambda_is_reproducible_from_the_toml_entry_point(tmp_path):
    """Spec Goal ("fed by an externally supplied Green function"):
    ``calc_eliashberg`` with ``[eliashberg] bond_channels = true`` and
    ``bond_green = <FLEX green.npz>`` reproduces the milestone ``lambda_t``
    that ``_build_bond_operator`` gives when called directly with the same
    green -- i.e. the acceptance number is reachable from a TOML input, not
    only from Python.
    """
    L, V = 16, 1.2
    indir = str(tmp_path / "input")
    outdir = str(tmp_path / "output")
    _write_model_inputs(indir, V)

    input_dict = {
        "mode": {"param": {"T": T, "CellShape": [L, L, 1],
                           "SubShape": [1, 1, 1], "Nmat": NMAT,
                           "filling": FILLING}},
        "file": {"input": {"interaction": {
                     "path_to_input": indir,
                     "Geometry": "geom.dat",
                     "Transfer": "transfer.dat",
                     "CoulombIntra": "coulombintra.dat",
                     "CoulombInter": "coulombinter.dat"}},
                 "output": {"path_to_output": outdir}},
        "eliashberg": {"chi0q_mode": "calc",
                       "bond_channels": True,
                       "bond_green": os.path.abspath(_fixture_path(L, V)),
                       "pairing_type": "triplet",
                       "init_gap": "f_x",
                       "solver_mode": "eigenvalue",
                       "num_eigenvalues": 6},
    }
    sc.calc_eliashberg(input_dict)

    vals, rows = {}, []
    with open(os.path.join(outdir, "eigenvalue.dat")) as fh:
        for line in fh:
            if line.startswith("#"):
                if "=" in line:
                    k, _, v = line.lstrip("# ").partition("=")
                    vals[k.strip()] = v.strip()
                continue
            rows.append([float(x) for x in line.split()])

    # the leading (odd-parity, f-like) eigenvalue of the file == the milestone
    lam_file = rows[0][0]
    assert lam_file == pytest.approx(LAMBDA_BOND_16[V], rel=1e-6)
    assert float(vals["lambda_rayleigh"]) == pytest.approx(
        LAMBDA_BOND_16[V], rel=1e-6)

    # ... and the provenance says WHICH green produced it
    assert vals["bond_green"].endswith("green_L16_V1.20.npz")
    assert "external" in vals["approximation"].lower()
    assert "bare" not in vals["approximation"].lower()
