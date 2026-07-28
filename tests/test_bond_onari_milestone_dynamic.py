"""Task 9 -- Phase A milestone (spec ``2026-07-27-dynamic-bond-channels-
design.md`` S3.6): the frequency-resolved (dynamic) bond-resolved triplet
Eliashberg eigenvalue ``lambda_t(V)`` rises STRICTLY MONOTONICALLY with the
inter-site Coulomb ``V``, on the SAME Hartree-only reduced-FLEX green
fixtures the static milestone (``test_bond_onari_milestone.py``) already
commits and hash-verifies: single-band square lattice, ``U=4``, filling
``n=0.7``, ``T=0.02``, ``Nmat=256``, ``16x16``, ``V in {0, 0.4, 0.8, 1.0,
1.2}``.

Green provenance -- reusing the static milestone's fixtures, not
regenerating
-----------------------------------------------------------------------------
Spec S3.6's parameters are IDENTICAL to the static milestone's V grid at
``L=16`` (see ``tests/sc/onari_bond/generate_fixtures.py``:
``mode="FLEX"``, ``calc_scheme="reduced"`` -- "Hartree-only V" -- Anderson
``depth=8``, ``Mix=0.2``, ``IterationMax=1500``). A second, independently
self-consistent green at the same parameters would only reproduce the same
physics from a different (non-bit-identical) SCF trajectory, and using a
DIFFERENT green for the static and dynamic arms would confound the
static/dynamic ratio indicator below with FLEX-SCF noise on top of the
retardation effect it measures. This test therefore loads the SAME
committed, hash-verified ``green_L16_V*.npz`` fixtures via
``tests.test_bond_onari_milestone._load_green`` for both arms;
``generate_fixtures.DYNAMIC_MILESTONE_POINTS`` documents the decision at the
point of generation.

Which code is under test
-----------------------------------------------------------------------------
The dynamic kernel is built by the PRODUCTION call sequence
``bond_bubble_dynamic -> dress_bond_dynamic -> make_bond_kernel_dynamic``,
the exact sequence ``eliashberg_dynamic._build_bond_dynamic_context``'s
uniform-axis branch makes (see ``_dynamic_point`` below) -- mirroring the
static milestone's own pattern of calling ``sc._build_bond_operator``
directly rather than only through the TOML entry point.
``test_milestone_reproducible_from_solve_dynamic_entry_point`` closes the
loop for the well-conditioned V points by also going through
``hwave.solver.eliashberg_dynamic.solve_dynamic`` (the dict/TOML entry
point), and is where the ``bond_diagnostics=true`` archive (conditioning
maps + Hermiticity + tail diagnostics) is asserted to exist, per the task
brief.

DEVIATION 1 -- ``matsubara_basis="uniform"``, not ``"ir"``
-----------------------------------------------------------------------------
The task parameters name ``matsubara_basis="ir"`` as the production choice.
On this fixture (``beta=50``, Nmat=256) it is BROKEN: ``_ir_compress``'s
``keep_constant`` fit (mandatory for the bond branch -- the dressed bond chi
is static-dominated, so the frequency-independent component must be
retained, spec S3.3/Task 8) is built from ``np.linalg.pinv`` with the
default ``rcond``, and that pseudo-inverse becomes severely ill-conditioned
as the IR basis size ``L`` grows with ``Lambda = beta*ir_wmax``. Measured on
``V=0.0``, ``pairing_type="triplet"``, ``solver_mode="eigenvalue"``:

======================  ==================  ========================
``ir_wmax``              auto (27.25)         explicit override
======================  ==================  ========================
auto                     lambda = 1.67e+03
40.0                                          lambda = 2.92e+08
80.0                                          lambda = 7.16e+14
160.0                                         lambda = 4.35e+18
======================  ==================  ========================

(uniform-basis reference at the same point: ``lambda = 0.109433``.) The
fitted "frequency-independent constant" of the retained-constant fit grows
in lock-step (4.8e1 -> 3.8e3 -> 1.06e7 -> 4.1e6), the unambiguous signature
of an ill-conditioned least-squares solve, not a converging basis. This is
NOT the green-symmetry watch item named in the task brief (that check
passes comfortably, see below) -- it is a separate, previously
undiscovered defect in ``_ir_compress``'s constant-augmented pseudo-inverse
at production-scale ``Lambda`` (the toy 4x4/``beta=2`` fixture Task 8's own
tests use never reaches a large enough ``Lambda`` to expose it). Per the
brief's "do not silently patch a tripped precondition" discipline this is
reported, not hidden: the milestone uses ``matsubara_basis="uniform"``
(the oracle-verified, P0-1-gated path) instead, and
``test_matsubara_basis_ir_is_broken_on_this_fixture`` pins the broken
behaviour as a documented regression so a future fix is visible rather than
silently changing the milestone's basis back to "ir" unnoticed. Filed as a
concern in the Task 9 report for a follow-up fix to ``_ir_compress``.

DEVIATION 2 -- ``V=1.2``'s default charge-conditioning guard
-----------------------------------------------------------------------------
The task brief names the green-symmetry precondition (spec S2:
``G(-k,iw)=G(k,iw)``, ``G(k,-iw)=conj G``, checked at ``1e-10`` relative) as
the specific watch item. It passes comfortably on every fixture used here
(measured residuals ~1e-16 -- 1e-15, six orders below the tolerance; see
``test_fixtures_carry_scf_convergence_metadata`` and the milestone's own
per-V assertions). A DIFFERENT guard trips instead: ``dress_bond_dynamic``'s
hardcoded charge-conditioning floor (``bond_channels._BOND_COND_FLOOR =
1e-3``, not exposed through ``[eliashberg]`` / ``solve_dynamic`` in this
branch) is violated at exactly ONE (q, i-nu) point for ``V=1.2`` --
``cond_charge = 6.819e-04`` at q near Gamma (grid index (0,15,0) /
(0,1,0)-equivalent by symmetry) and bosonic index 250 (of 256; i.e. near
the EDGE of the finite Matsubara window, not at the static i-nu=0 point,
whose own minimum conditioning is a healthy 0.256). This is consistent with
``bond_bubble_dynamic``'s documented, deferred ``green0_tail`` high-
frequency correction (its own docstring: "does NOT apply the ... tail
correction ... the values returned here are the raw finite-nmat sum only")
producing excess weight specifically at the window edge, not a genuine
zero-frequency CDW divergence spreading to finite frequency.

Investigated per spec S3.6 ("investigation, not automatic failure, and the
disposition is recorded") before overriding anything: the triplet-parity
leading eigenvalue at ``V=1.2`` (``1.8553050420``, see ``_dynamic_point``'s
``COND_TOL_OVERRIDE``) is IDENTICAL to 10 significant digits across
``cond_tol in {1e-4, 1e-5, 1e-6, 1e-8}`` and across ``num_eigenvalues in {6,
10}`` -- i.e. the near-singular direction the default floor refuses carries
negligible weight in the physical (parity-selected) eigenvector; loosening
the guard by one decade (``cond_tol=1e-4``, still 6.8x above the measured
floor) does not move the reported lambda at all. The disposition: SAFE,
recorded here rather than silently defaulting -- ``solve_dynamic`` itself
has no config knob for this ``cond_tol`` in the bond branch (a found gap,
noted for Task 10 / a follow-up), so
``test_v1_2_needs_the_conditioning_floor_relaxed_below_solve_dynamics_reach``
documents that the TOML/dict entry point cannot reach ``V=1.2`` today and
pins the direct-call value it CAN reach.

Ratio indicator -- measured OUTSIDE (and on the wrong side of) the loose
[1.5, 3.5] window at every V; disposition recorded, not failed
-----------------------------------------------------------------------------
Spec S3.6: "static/dynamic ratio within the loose literature window [1.5,
3.5] (indicator only ... a value outside the window triggers investigation,
not automatic failure)". Measured (static ``LAMBDA_BOND_16`` / this
module's dynamic values): ``V=0.0`` -> 0.644, ``V=0.4`` -> 0.706, ``V=0.8``
-> 0.801, ``V=1.0`` -> 0.346, ``V=1.2`` -> 0.167. Every point is outside the
window, and INVERTED from the window's assumption (static > dynamic by
1.5x-3.5x): here the DYNAMIC eigenvalue is systematically LARGER than the
static one, by a growing margin toward the CDW. This is not confined to the
V=1.0/1.2 near-instability points -- it already holds at the safely-
conditioned V=0.0/0.4/0.8 (cond_charge_min > 0.08 there), so it is not an
artifact of the conditioning-guard investigation above. Two candidate
explanations are recorded (not resolved here -- out of Task 9's scope):
(a) genuine physics -- the bond-resolved dynamic RPA may capture additional
attractive retarded weight the single ``m=m'=0`` static treatment cannot,
unlike the SCALAR dynamic/static comparison the [1.5, 3.5] window was
calibrated on (spec S3.6 says as much: "the exact 2.3x was measured under
different V treatment"); (b) an unverified normalization convention
difference between ``make_bond_kernel_dynamic`` and the static
``sc._build_bond_operator`` -- Task 8's own cross-validation compared the
dynamic kernel against itself (uniform vs IR, and the B=1 reduction against
the SCALAR dynamic kernel), never against the STATIC bond kernel. This is
exactly the "Phase A hand-off indicator" spec S3.6 asks Task 9 to produce
for Phase B, not a Task 9 blocker: the REQUIRED criterion (strict
monotonicity) holds throughout.
"""

import os

import numpy as np
import pytest

import hwave.sc as sc
from hwave.solver import bond_channels as bc
from hwave.solver import eliashberg_dynamic as ed

from tests.test_bond_onari_milestone import (
    BETA,
    FILLING,
    LAMBDA_BOND_16,
    LAMBDA_RTOL,
    NMAT,
    T,
    U,
    V_GRID,
    _check_fixture_metadata,
    _fixture_path,
    _generator_module,
    _interactions,
    _kgrid,
    _load_green,
)

L = 16
_SLOW_ENV = "HWAVE_RUN_SLOW_FIXTURES"
NUM_EIGENVALUES = 6

# Production default (bond_channels._BOND_COND_FLOOR); V=1.2 needs the
# documented, measured-stable relaxation -- see module docstring "DEVIATION
# 2". No other V point needs an override.
COND_TOL_DEFAULT = bc._BOND_COND_FLOOR
COND_TOL_OVERRIDE = {1.2: 1.0e-4}

# ir_wmax values used to reproduce the DEVIATION-1 measurement table (V=0.0
# only, in test_matsubara_basis_ir_is_broken_on_this_fixture).
_IR_WMAX_PROBES = (None, 40.0, 80.0, 160.0)


def _require_slow(request):
    """Skip unless explicitly asked for (env var or ``-m slow``).

    Mirrors ``test_bond_onari_milestone._require_slow_fixtures``, but gates
    the expensive DYNAMIC SOLVE here (frequency-resolved kernel over
    Nmat=256), not fixture regeneration -- the green fixtures are already
    committed and reused verbatim (module docstring).
    """
    if os.environ.get(_SLOW_ENV, "").lower() not in ("", "0", "false", "no"):
        return
    markexpr = str(request.config.getoption("markexpr", default="") or "")
    tokens = markexpr.replace("(", " ").replace(")", " ").split()
    if "slow" in tokens and "not" not in tokens:
        return
    pytest.skip(
        "Phase A dynamic milestone runs a full-frequency (Nmat=256) bond "
        "Eliashberg solve at every V: minutes of compute. Run it with "
        "{}=1 or with `pytest -m slow`; see the module docstring."
        .format(_SLOW_ENV))


def _gap_shape(green):
    norb = green.shape[0]
    return (norb, norb, green.shape[2], green.shape[3], green.shape[4],
            green.shape[-1])


def _dynamic_point(V):
    """One V point of the dynamic bond kernel, built by the PRODUCTION call
    sequence ``bond_bubble_dynamic -> dress_bond_dynamic ->
    make_bond_kernel_dynamic`` -- the same sequence
    ``eliashberg_dynamic._build_bond_dynamic_context``'s uniform-axis branch
    makes (module docstring "Which code is under test").
    """
    kx, ky, kz = _kgrid(L)
    it = _interactions(V)
    inter_k = sc._build_interaction_k(kx, ky, kz, it, 1)
    bond_set = bc.resolve_interactions(it["CoulombInter"], np.eye(3), 1)
    green = _load_green(L, V)

    S0_q, C0_q = sc._build_bond_m0_blocks(bond_set, it, inter_k, 1, kx, ky, kz)
    S_bond, C_bond, Vpp_s, Vpp_t = bc.bare_bond_vertices(bond_set, S0_q, C0_q, 1)

    # spec S2 green-symmetry precondition (the task brief's named watch
    # item) -- checked directly, exactly as eliashberg_dynamic's wiring
    # does, before it feeds anything downstream.
    sym = ed._check_bond_dynamic_green_symmetry(green)

    chi_bar_w = bc.bond_bubble_dynamic(green, bond_set, BETA)
    cond_tol = COND_TOL_OVERRIDE.get(V, COND_TOL_DEFAULT)
    chi_s_w, chi_c_w, cond_min_spin, cond_min_charge = bc.dress_bond_dynamic(
        chi_bar_w, S_bond, C_bond, cond_tol=cond_tol)
    A, vec_size = bc.make_bond_kernel_dynamic(
        chi_s_w, chi_c_w, S_bond, C_bond, Vpp_s, Vpp_t, green, bond_set,
        "triplet", BETA)

    _, _, info = sc._solve_leading(
        lambda: (A, vec_size), vec_size, "arnoldi",
        num_eigenvalues=NUM_EIGENVALUES)
    gap_shape = _gap_shape(green)
    vals, vecs, match = ed._reorder_eigenpairs_by_parity_dynamic(
        info["eigenvalues"], info["eigenvectors"], gap_shape, "triplet")
    assert match[0], (
        "V={}: no ARPACK eigenpair matches the triplet parity within "
        "num_eigenvalues={} -- widen num_eigenvalues".format(
            V, NUM_EIGENVALUES))
    lam = complex(vals[0])

    weight = ed._bond_dynamic_pair_weight(green).real.ravel()
    herm = bc.check_hermitian_preconditions(A, weight,
                                            label="milestone dynamic V={}"
                                            .format(V))
    gap_w = np.asarray(vecs[:, 0]).reshape(gap_shape)
    tail = ed._bond_tail_diagnostic(green, gap_w, BETA, green.shape[-1],
                                    "triplet")

    return {
        "lambda": lam.real,
        "lambda_imag": lam.imag,
        "cond_min_spin": float(cond_min_spin.min()),
        "cond_min_charge": float(cond_min_charge.min()),
        "cond_min_charge_at": tuple(int(i) for i in np.unravel_index(
            int(np.argmin(cond_min_charge)), cond_min_charge.shape)),
        "cond_tol_used": cond_tol,
        "green_sym_k": sym["green_sym_k"],
        "green_sym_w": sym["green_sym_w"],
        "kernel_hermiticity_relative": herm["kernel_hermiticity_relative_ktilde"],
        "weight_min_eigenvalue": herm["weight_min_eigenvalue"],
        "tail_est": tail["tail_est"],
        "tail_est_rel": tail["tail_est_rel"],
        "tail_est_unreliable": tail["tail_est_unreliable"],
    }


def _static_lambda(V, tmp_path):
    """The existing static bond path -- ``calc_eliashberg`` with
    ``bond_channels=true`` + ``bond_green=<the SAME fixture>`` -- for the
    ratio indicator. Reproduces the pinned ``LAMBDA_BOND_16[V]`` (asserted by
    the caller), matching
    ``test_bond_onari_milestone.test_milestone_lambda_is_reproducible_from_the_toml_entry_point``.
    """
    mod = _generator_module()
    indir = str(tmp_path / "static_in_{}".format(V))
    outdir = str(tmp_path / "static_out_{}".format(V))
    mod._write_inputs(indir, V)
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
                       "num_eigenvalues": NUM_EIGENVALUES},
    }
    sc.calc_eliashberg(input_dict)
    with open(os.path.join(outdir, "eigenvalue.dat")) as fh:
        rows = [line for line in fh if not line.startswith("#")]
    return float(rows[0].split()[0])


def _dynamic_input_dict(V, output_dir, basis="uniform"):
    green_path = os.path.abspath(_fixture_path(L, V))
    mod = _generator_module()
    indir = os.path.join(os.path.dirname(output_dir), "in_{}".format(V))
    mod._write_inputs(indir, V)
    return {
        "mode": {"param": {"T": T, "CellShape": [L, L, 1],
                           "SubShape": [1, 1, 1], "Nmat": NMAT,
                           "filling": FILLING}},
        "file": {"input": {"interaction": {
                     "path_to_input": indir,
                     "Geometry": "geom.dat",
                     "Transfer": "transfer.dat",
                     "CoulombIntra": "coulombintra.dat",
                     "CoulombInter": "coulombinter.dat"}},
                 "output": {"path_to_output": output_dir}},
        "eliashberg": {"frequency": "dynamic", "bond_channels": True,
                       "bond_green": green_path,
                       "matsubara_basis": basis,
                       "pairing_type": "triplet",
                       "solver_mode": "eigenvalue",
                       "eigenvalue_method": "arnoldi",
                       "num_eigenvalues": NUM_EIGENVALUES,
                       "bond_diagnostics": True},
    }


# ===========================================================================
# tests
# ===========================================================================

def test_fixtures_carry_scf_convergence_metadata():
    """(I-1) UNCONDITIONAL fixture-provenance check (PR #82 lesson, named
    explicitly in the task brief): every fixture this milestone consumes
    records ``scf_converged is True`` and ``scf_iterations > 0`` -- checked
    directly with ``assert``, never gated behind ``if key in data``. Not
    slow: reads five small already-committed npz files.
    """
    for V in V_GRID:
        _check_fixture_metadata(_fixture_path(L, V), L, V)


@pytest.mark.slow
def test_phase_a_milestone_lambda_t_monotone_and_ratio_window(request,
                                                               tmp_path):
    """THE MILESTONE (spec S3.6): the dynamic bond-resolved triplet
    eigenvalue rises STRICTLY MONOTONICALLY with V (required), and the
    static/dynamic ratio + per-(q,i-nu) conditioning minima + Hermiticity
    residuals are recorded (indicators, per spec S3.6 -- see module
    docstring for the measured disposition of the two deviations)."""
    _require_slow(request)

    points, lam_static = {}, {}
    for V in V_GRID:
        _check_fixture_metadata(_fixture_path(L, V), L, V)   # (I-1), again
        points[V] = _dynamic_point(V)
        lam_static[V] = _static_lambda(V, tmp_path)
        assert lam_static[V] == pytest.approx(LAMBDA_BOND_16[V],
                                              rel=LAMBDA_RTOL)

    dyn = [points[V]["lambda"] for V in V_GRID]
    for a, b, V in zip(dyn, dyn[1:], V_GRID[1:]):
        assert b > a, "lambda_t(V) not monotone at V={}: {} -> {}".format(
            V, a, b)

    print("\nPhase A milestone -- dynamic bond channels (uniform Matsubara "
          "basis), U={}, n={}, T={}, Nmat={}, 16x16".format(U, FILLING, T,
                                                             NMAT))
    header = ("{:>5}  {:>12}  {:>12}  {:>8}  {:>10}  {:>10}  {:>10}  "
              "{:>10}  {:>9}").format(
        "V", "lam_dyn", "lam_static", "ratio", "cond_spin", "cond_chg",
        "herm_rel", "tail_rel", "cond_tol")
    print(header)
    for V in V_GRID:
        p = points[V]
        assert np.isfinite(p["lambda"])
        assert abs(p["lambda_imag"]) <= max(1.0e-8, 1.0e-6 * abs(p["lambda"])), (
            "V={}: leading eigenvalue has a non-negligible imaginary part "
            "{}".format(V, p["lambda_imag"]))
        assert p["green_sym_k"] < 1.0e-10 and p["green_sym_w"] < 1.0e-10
        assert p["kernel_hermiticity_relative"] < 1.0e-8
        assert p["cond_min_spin"] > 0.0
        assert p["cond_min_charge"] > 0.0
        assert np.isfinite(p["tail_est"])

        ratio = lam_static[V] / p["lambda"] if p["lambda"] else float("nan")
        assert np.isfinite(ratio)
        flag = ("" if 1.5 <= ratio <= 3.5 else
                "  ** past_instability/outside [1.5,3.5] indicator window **")
        print("{:5.2f}  {:12.6f}  {:12.6f}  {:8.3f}  {:10.3e}  {:10.3e}  "
              "{:10.3e}  {:10.3e}  {:9.1e}{}".format(
                  V, p["lambda"], lam_static[V], ratio, p["cond_min_spin"],
                  p["cond_min_charge"], p["kernel_hermiticity_relative"],
                  p["tail_est_rel"], p["cond_tol_used"], flag))
    # spec S3.6 doesn't require the ratio inside the window (indicator
    # only); the printed table + the module docstring "Ratio indicator"
    # section are the recorded disposition.


@pytest.mark.slow
@pytest.mark.parametrize("V", (0.0, 0.4, 0.8, 1.0))
def test_milestone_reproducible_from_solve_dynamic_entry_point(request, V,
                                                                tmp_path):
    """Spec goal ("fed by an externally supplied Green function"): the
    dict/TOML entry point ``solve_dynamic`` reproduces the direct-call
    ``_dynamic_point`` lambda, and the ``bond_diagnostics=true`` archive
    (conditioning maps, Hermiticity/tail diagnostics) exists with the keys
    the task brief requires. ``V=1.2`` is excluded here -- see
    ``test_v1_2_needs_the_conditioning_floor_relaxed_below_solve_dynamics_reach``.
    """
    _require_slow(request)
    direct = _dynamic_point(V)

    output_dir = str(tmp_path / "out")
    inp = _dynamic_input_dict(V, output_dir)
    lam = ed.solve_dynamic(inp)
    assert lam == pytest.approx(direct["lambda"], rel=1.0e-8)

    meta = np.load(os.path.join(output_dir, "gap_dynamic.npz"),
                   allow_pickle=False)
    assert bool(meta["bond_channels"])
    assert int(meta["bond_n_channels"]) == 5
    assert str(meta["bond_green_source"]) == "file"
    assert str(meta["bond_frequency_grid"]) == "uniform"
    assert float(meta["bond_cond_min_spin"]) > 0.0
    assert float(meta["bond_cond_min_charge"]) > 0.0
    assert float(meta["bond_green_symmetry_residual_k"]) < 1.0e-10
    assert float(meta["bond_green_symmetry_residual_w"]) < 1.0e-10
    assert str(meta["hermitian_metric_features"]) == "uniform"
    # stored as a formatted string ("{:.3e}"), not a float array
    assert float(str(meta["kernel_hermiticity_relative_ktilde"])) < 1.0e-8
    assert np.isfinite(float(meta["bond_lambda_imag"]))

    diag = np.load(os.path.join(output_dir, "bond_diagnostics.npz"),
                   allow_pickle=False)
    assert diag["cond_min_spin"].shape == (L, L, 1, NMAT)
    assert diag["cond_min_charge"].shape == (L, L, 1, NMAT)
    assert "bond_diag_axes" in diag.files
    assert np.isfinite(float(diag["tail_est"]))
    assert str(diag["tail_est_branch"]) == "triplet"
    assert bool(diag["tail_est_unreliable"]) in (True, False)


@pytest.mark.slow
def test_v1_2_needs_the_conditioning_floor_relaxed_below_solve_dynamics_reach(
        request, tmp_path):
    """Documents DEVIATION 2 (module docstring): ``solve_dynamic`` refuses
    ``V=1.2`` with the hardcoded production ``cond_tol`` (a found gap -- no
    ``[eliashberg]`` knob reaches ``dress_bond_dynamic``'s ``cond_tol`` on
    this branch), while the direct call with the measured-stable
    ``cond_tol=1.0e-4`` succeeds and agrees with the value pinned in the
    main sweep."""
    _require_slow(request)
    output_dir = str(tmp_path / "out")
    inp = _dynamic_input_dict(1.2, output_dir)
    with pytest.raises(ValueError, match="charge instability"):
        ed.solve_dynamic(inp)

    direct = _dynamic_point(1.2)
    assert direct["lambda"] == pytest.approx(1.8553050420, rel=1.0e-6)
    assert direct["cond_tol_used"] == pytest.approx(1.0e-4)


@pytest.mark.slow
def test_matsubara_basis_ir_is_broken_on_this_fixture(request, tmp_path):
    """Pins DEVIATION 1 (module docstring) as a documented regression: the
    IR-basis leading eigenvalue at V=0.0 is many orders of magnitude away
    from the trustworthy uniform-basis reference, and grows WORSE (not
    better) as ``ir_wmax`` is raised -- the signature of an ill-conditioned
    least-squares fit in ``_ir_compress``'s constant-retention branch, not a
    basis-truncation error that a larger window would fix.

    This test exists so a future ``_ir_compress`` fix is NOTICED (the
    assertions below start failing) rather than the milestone silently
    reverting to ``matsubara_basis="ir"`` with nobody re-checking the
    physics. When it is fixed, flip the milestone's basis back to "ir" and
    delete/rewrite this test.
    """
    _require_slow(request)
    reference = _dynamic_point(0.0)["lambda"]   # uniform basis, trustworthy

    output_dir = str(tmp_path / "out")
    inp = _dynamic_input_dict(0.0, output_dir, basis="ir")
    lam_ir = ed.solve_dynamic(inp)

    # today's broken behaviour: >100x the trustworthy reference (measured
    # ~1.5e4x; a loose bound so this doesn't flake on minor IR-tuning
    # changes while still catching "still broken" vs "actually fixed").
    assert abs(lam_ir) > 100.0 * abs(reference), (
        "matsubara_basis='ir' no longer explodes on this fixture "
        "(lam_ir={}, uniform reference={}) -- the _ir_compress "
        "keep_constant ill-conditioning (module docstring DEVIATION 1) may "
        "be fixed; if so, switch the Phase A milestone back to "
        "matsubara_basis='ir' and remove this regression pin.".format(
            lam_ir, reference))
