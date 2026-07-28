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
loop for every V (including 1.2, via ``[eliashberg] bond_cond_tol`` -- Task 9
review IMPORTANT-3) by also going through
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
of an ill-conditioned least-squares solve, not a converging basis.

**Root-cause measurement (Task 9 review MINOR):** the fit design matrix
itself, ``ir_axis.IRAxis.uniform_matrices(nmat, with_constant=True)``, has
``cond() = 3.11e7`` at the production scale (``Lambda = beta*ir_wmax =
1362.7``, 47 bosonic nodes), vs. ``cond() = 184`` at the toy scale
(``Lambda = 26``, 18 nodes) that Task 8's own tests exercise -- a
~170,000x increase, entirely explaining why the toy fixture never exposed
this. This is NOT the green-symmetry watch item named in the task brief
(that check passes comfortably, see below) -- it is a separate,
previously undiscovered defect in ``_ir_compress``'s constant-augmented
pseudo-inverse at production-scale ``Lambda``. Per the brief's "do not
silently patch a tripped precondition" discipline this is reported, not
hidden: the milestone uses ``matsubara_basis="uniform"`` (the
oracle-verified, P0-1-gated path) instead, and
``test_matsubara_basis_ir_is_broken_on_this_fixture`` (slow, this module)
plus ``test_dynamic_bond_ir_is_catastrophically_broken_at_nmat_32``
(FAST/default-CI, ``tests/test_sc_bond_dynamic.py``, same defect at toy
scale: ``Lambda=26`` uniform=1.075495 vs. IR=-3.081435 -- wrong SIGN) pin
the broken behaviour as a documented regression so a future fix is
noticed rather than the milestone silently reverting to ``"ir"``
unnoticed. Filed as a concern in the Task 9 report for a follow-up fix to
``_ir_compress``/``ir_axis.uniform_matrices``.

DEVIATION 2 + Ratio indicator -- SINGLE ROOT CAUSE (Task 9 review
IMPORTANT-1/4/5): the deferred ``green0_tail`` correction, not two
unrelated puzzles
-----------------------------------------------------------------------------
**The first round's "V=1.2's near-singular point is safe to override,
verified stable across cond_tol" claim was WRONG -- refuted by measurement,
not defended.** ``dress_bond_dynamic``'s ``cond_tol`` argument only gates the
``ValueError`` check; it never touches ``np.linalg.solve``, so varying
``cond_tol`` and re-solving cannot show anything about whether the
near-singular block actually MATTERS -- the computation is bit-identical at
every ``cond_tol`` that lets it run. That is why it "was stable": the test
could not have failed. Replaced with a measurement that CAN fail: with the
production ``chi_bar_w`` for ``V=1.2``, replacing (or, equivalently,
zeroing) every ``(q, i-nu)`` block with charge-conditioning score ``< 0.01``
by the SAME-q STATIC (``i-nu=0``) value, then re-dressing and re-solving,
moves the leading triplet eigenvalue from ``1.855305`` to ``0.254`` -- an
**86% shift**. The near-singular region does not carry negligible weight; it
DOMINATES the reported ``V=1.2`` eigenvalue.

Three further measurements (this module's ``test_...`` functions and the
Task 9 report) converge on a single explanation -- ``bond_bubble_dynamic``'s
own documented, deferred ``green0_tail`` high-frequency correction (its
docstring: "does NOT apply the ... tail correction ... the values returned
here are the raw finite-``nmat`` sum only"):

1. **A non-decaying, V-independent edge floor in** ``||chi_bar(q, i-nu)||``.
   The bubble's amplitude (max over q of the Frobenius norm of the
   ``ND x ND`` block) should DECAY toward the edges of the Matsubara window
   (the physical pair bubble falls off ``~1/nu^2`` at large frequency); it
   does not. At every V it plateaus at ``43-45%`` of its zero-frequency peak
   at the exact window edge (``|n~| = nmat/2``), and the SAME ``43-45%``
   floor appears at every V in ``{0.0, 0.4, 0.8, 1.0, 1.2}`` -- consistent
   with a fixed discretization artifact of the missing tail correction, not
   a V-dependent physical feature.
2. **An unphysical conditioning ordering.** The GLOBAL minimum of the
   charge-channel RPA-denominator conditioning score sits at the window
   EDGE, not at the physical static (``i-nu=0``) point, for every V >= 0.4 --
   and the gap between them grows with V (V=0.4: edge 0.347 vs. static
   0.434; V=0.8: edge 0.087 vs. static 0.341; V=1.0: edge 0.0013 vs. static
   0.298; V=1.2: edge 0.00068 vs. static 0.256). A genuine approach to a
   zero-frequency CDW instability would make the STATIC point the worst one,
   not a point 122 bosonic indices away from it.
3. **The ablation shift itself is V-dependent and changes SIGN.** Flattening
   the outer quarter-shell (``|n~| > 3*nmat/8`` -- ``bond_channels.
   tail_estimate``'s own outer-shell definition) to the static value at
   EVERY V (not just 1.2) shifts lambda by +17% (V=0.0), +19% (V=0.4), +17%
   (V=0.8), **-59%** (V=1.0), **-86%** (V=1.2, see
   ``LAMBDA_SHELL_FLAT_16``). The artifact is small and roughly V-independent
   at weak coupling, then becomes the DOMINANT contribution once the charge
   channel approaches its RPA instability (V >= 1.0) -- exactly what an
   uncorrected high-frequency tail feeding back through a near-singular
   denominator would do.

**Disposition (recorded, not silently patched):** ``[eliashberg]
bond_cond_tol=1e-4`` (spec S3.6 hand-off, review IMPORTANT-3) is still used
to let ``V=1.2`` COMPUTE at all -- the guard's job is only to refuse
division by an exact/near-exact zero, and 6.819e-04 is not that -- but the
resulting number must NOT be read as "the physically validated triplet
lambda at V=1.2". **Phase B must not treat ``lambda_dyn``'s absolute scale
as physical at ANY V until ``green0_tail`` lands**; V=1.0 and V=1.2 are the
least trustworthy (largest ablation shift), V=0.0/0.4/0.8 the most (+17-19%,
roughly constant, plausibly closer to the true small-``O(1/Nmat)``
correction the static milestone's own fixtures already absorb via a
different mechanism). This also resolves the ratio-indicator puzzle from the
first round: the flat-chi collapse test (``test_dynamic_kernel_with_
frequency_flat_chi_collapses_onto_static``, ``tests/test_sc_bond_dynamic.
py``) proves ``make_bond_kernel_dynamic`` reduces to the static
``make_bond_kernel_parts`` to machine precision (relative difference
``1.6e-15``) when fed frequency-flat chi -- candidate (b) "a normalization
bug in the dynamic kernel builder" is DEAD. The surviving explanation for
the measured ``static/dynamic`` ratio sitting outside (and inverted from)
the loose ``[1.5, 3.5]`` window is candidate (a), the frequency dependence
of ``chi_bar(q, i-nu)`` -- but per the single-root-cause finding above, a
substantial and V-growing fraction of that frequency dependence is now
understood to be the SAME uncorrected tail artifact, not (only, or even
mostly) genuine retardation physics. The REQUIRED spec S3.6 criterion
(strict monotonicity of the reported, reproducible ``lambda_dyn(V)``
sequence) still holds, and -- reassuringly -- the artifact-suppressed
(outer-shell-flattened) sequence is ALSO monotone, just far gentler in
slope and different in absolute scale (``LAMBDA_SHELL_FLAT_16``). This is
exactly the kind of Phase-A-to-Phase-B hand-off finding spec S3.6 asks
Task 9 to produce, not a Task 9 blocker.

Green-symmetry precondition (the task brief's named watch item): passes
comfortably at every V (measured residuals ~1e-16, six orders below the
1e-10 tolerance) -- see ``test_fixtures_carry_scf_convergence_metadata`` and
the per-V assertions below.
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

# Production default (bond_channels._BOND_COND_FLOOR); V=1.2 needs
# [eliashberg] bond_cond_tol relaxed (spec S3.6 hand-off, review
# IMPORTANT-3) -- see the module docstring's single-root-cause section for
# why this is safe to COMPUTE but not to read as physically validated.
COND_TOL_DEFAULT = bc._BOND_COND_FLOOR
COND_TOL_OVERRIDE = {1.2: 1.0e-4}

# Pinned dynamic lambda_t(V) (uniform basis, this module's own
# _dynamic_point) -- Task 9 review IMPORTANT-2. Full printed precision from
# a cold rerun; reproduced to the asserted rtol on an independent rerun
# during the review-fix round.
LAMBDA_BOND_DYN_16 = {
    0.0: 0.10943297029791513,
    0.4: 0.13102477861872136,
    0.8: 0.15211304475749632,
    1.0: 0.47571721381791665,
    1.2: 1.8553050420058088,
}
LAMBDA_DYN_RTOL = 1.0e-6

# Outer-shell-flattened diagnostic (module docstring point 3 / review
# IMPORTANT-1): chi_bar_w's outer quarter-shell (|n~| > 3*nmat/8, matching
# bond_channels.tail_estimate's own outer-shell definition) replaced by the
# SAME-q static (i-nu=0) value before dressing+solving. A REAL, falsifiable
# sensitivity measurement (unlike the retracted cond_tol-variation claim):
# it moves lambda by 17-86% depending on V. Looser rtol than the primary
# pin -- it is a diagnostic on a hand-modified input, not the production
# path.
LAMBDA_SHELL_FLAT_16 = {
    0.0: 0.12763298481459254,
    0.4: 0.15531005744268794,
    0.8: 0.17871576498360847,
    1.0: 0.1961645311051929,
    1.2: 0.2686665272119033,
}
LAMBDA_SHELL_FLAT_RTOL = 1.0e-3


def _require_slow(request):
    """Skip unless explicitly asked for (env var or ``-m slow``).

    Mirrors ``test_bond_onari_milestone._require_slow_fixtures``, but gates
    the expensive DYNAMIC SOLVE here (frequency-resolved kernel over
    Nmat=256, run TWICE per V once the outer-shell ablation is included),
    not fixture regeneration -- the green fixtures are already committed
    and reused verbatim (module docstring).
    """
    if os.environ.get(_SLOW_ENV, "").lower() not in ("", "0", "false", "no"):
        return
    markexpr = str(request.config.getoption("markexpr", default="") or "")
    tokens = markexpr.replace("(", " ").replace(")", " ").split()
    if "slow" in tokens and "not" not in tokens:
        return
    pytest.skip(
        "Phase A dynamic milestone runs full-frequency (Nmat=256) bond "
        "Eliashberg solves at every V: minutes of compute. Run it with "
        "{}=1 or with `pytest -m slow`; see the module docstring."
        .format(_SLOW_ENV))


def _gap_shape(green):
    norb = green.shape[0]
    return (norb, norb, green.shape[2], green.shape[3], green.shape[4],
            green.shape[-1])


def _solve_bond_dynamic(chi_bar_w, S_bond, C_bond, Vpp_s, Vpp_t, green,
                        bond_set, cond_tol):
    """Dress + build + solve the dynamic bond kernel from an already-built
    ``chi_bar_w`` -- the shared tail of ``_dynamic_point``, reused for both
    the production (raw) computation and the outer-shell-flattened
    diagnostic ablation (module docstring, review IMPORTANT-1) so the two
    only differ in ``chi_bar_w`` itself.

    Returns ``(lam: complex, cond_min_spin, cond_min_charge, A, vec_size)``
    -- the operator is returned too so a caller (the Hermiticity check in
    ``_dynamic_point``) can reuse it instead of rebuilding a third time.
    """
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
        "no ARPACK eigenpair matches the triplet parity within "
        "num_eigenvalues={} -- widen num_eigenvalues".format(
            NUM_EIGENVALUES))
    return complex(vals[0]), cond_min_spin, cond_min_charge, A, vec_size


def _flatten_outer_shell(chi_bar_w, nmat):
    """Replace the outer quarter-shell (``|n~| > 3*nmat/8``) with the
    SAME-q static (``i-nu = nmat//2``) value (module docstring point 3).
    """
    nb2 = nmat // 2
    shell_cut = 3 * nmat // 8
    out = chi_bar_w.copy()
    out[:, :, :, :shell_cut] = chi_bar_w[:, :, :, nb2:nb2 + 1]
    out[:, :, :, nmat - shell_cut:] = chi_bar_w[:, :, :, nb2:nb2 + 1]
    return out


def _dynamic_point(V):
    """One V point of the dynamic bond kernel, built by the PRODUCTION call
    sequence ``bond_bubble_dynamic -> dress_bond_dynamic ->
    make_bond_kernel_dynamic`` -- the same sequence
    ``eliashberg_dynamic._build_bond_dynamic_context``'s uniform-axis branch
    makes (module docstring "Which code is under test"). Also computes the
    outer-shell-flattened diagnostic ablation (module docstring point 3).
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
    lam, cond_min_spin, cond_min_charge, A, _ = _solve_bond_dynamic(
        chi_bar_w, S_bond, C_bond, Vpp_s, Vpp_t, green, bond_set, cond_tol)

    # outer-shell-flattened ablation (module docstring point 3 / review
    # IMPORTANT-1): a REAL sensitivity measurement, not a vacuous one --
    # cond_tol only gates the error check, so it can never move chi_bar_w
    # itself; this actually replaces the suspect data.
    nmat = chi_bar_w.shape[3]
    chi_bar_flat = _flatten_outer_shell(chi_bar_w, nmat)
    lam_flat, _, _, _, _ = _solve_bond_dynamic(
        chi_bar_flat, S_bond, C_bond, Vpp_s, Vpp_t, green, bond_set,
        cond_tol)

    weight = ed._bond_dynamic_pair_weight(green).real.ravel()
    herm = bc.check_hermitian_preconditions(A, weight,
                                            label="milestone dynamic V={}"
                                            .format(V))

    return {
        "lambda": lam.real,
        "lambda_imag": lam.imag,
        "lambda_shell_flat": lam_flat.real,
        "cond_min_spin": float(cond_min_spin.min()),
        "cond_min_charge": float(cond_min_charge.min()),
        "cond_min_charge_at": tuple(int(i) for i in np.unravel_index(
            int(np.argmin(cond_min_charge)), cond_min_charge.shape)),
        "cond_tol_used": cond_tol,
        "green_sym_k": sym["green_sym_k"],
        "green_sym_w": sym["green_sym_w"],
        "kernel_hermiticity_relative": herm["kernel_hermiticity_relative_ktilde"],
        "weight_min_eigenvalue": herm["weight_min_eigenvalue"],
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


def _dynamic_input_dict(V, output_dir, basis="uniform", cond_tol=None):
    green_path = os.path.abspath(_fixture_path(L, V))
    mod = _generator_module()
    indir = os.path.join(os.path.dirname(output_dir), "in_{}".format(V))
    mod._write_inputs(indir, V)
    eli = {"frequency": "dynamic", "bond_channels": True,
          "bond_green": green_path,
          "matsubara_basis": basis,
          "pairing_type": "triplet",
          "solver_mode": "eigenvalue",
          "eigenvalue_method": "arnoldi",
          "num_eigenvalues": NUM_EIGENVALUES,
          "bond_diagnostics": True}
    if cond_tol is not None:
        eli["bond_cond_tol"] = cond_tol
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
        "eliashberg": eli,
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
    eigenvalue rises STRICTLY MONOTONICALLY with V (required), pinned
    exactly (review IMPORTANT-2); the static/dynamic ratio, the
    outer-shell-flattened sensitivity ablation (review IMPORTANT-1), and
    per-(q,i-nu) conditioning minima + Hermiticity residuals are recorded
    (indicators, per spec S3.6 -- see module docstring for the single-
    root-cause disposition)."""
    _require_slow(request)

    points, lam_static = {}, {}
    for V in V_GRID:
        _check_fixture_metadata(_fixture_path(L, V), L, V)   # (I-1), again
        points[V] = _dynamic_point(V)
        lam_static[V] = _static_lambda(V, tmp_path)
        assert lam_static[V] == pytest.approx(LAMBDA_BOND_16[V],
                                              rel=LAMBDA_RTOL)
        assert points[V]["lambda"] == pytest.approx(LAMBDA_BOND_DYN_16[V],
                                                     rel=LAMBDA_DYN_RTOL)
        assert points[V]["lambda_shell_flat"] == pytest.approx(
            LAMBDA_SHELL_FLAT_16[V], rel=LAMBDA_SHELL_FLAT_RTOL)

    dyn = [points[V]["lambda"] for V in V_GRID]
    for a, b, V in zip(dyn, dyn[1:], V_GRID[1:]):
        assert b > a, "lambda_t(V) not monotone at V={}: {} -> {}".format(
            V, a, b)

    # The artifact-suppressed sequence is ALSO monotone -- a real,
    # falsifiable secondary check (module docstring): if a future
    # green0_tail fix or an unrelated regression broke this, it would fail
    # here even though the primary (raw) sequence above might still pass.
    dyn_flat = [points[V]["lambda_shell_flat"] for V in V_GRID]
    for a, b, V in zip(dyn_flat, dyn_flat[1:], V_GRID[1:]):
        assert b > a, (
            "outer-shell-flattened lambda_t(V) not monotone at V={}: "
            "{} -> {}".format(V, a, b))

    print("\nPhase A milestone -- dynamic bond channels (uniform Matsubara "
          "basis), U={}, n={}, T={}, Nmat={}, 16x16".format(U, FILLING, T,
                                                             NMAT))
    print("NOTE: per the module docstring's single-root-cause finding, "
          "lambda_dyn's ABSOLUTE SCALE is not yet physically validated "
          "(pending the deferred green0_tail correction) -- treat the "
          "numbers below as a pinned, reproducible COMPUTATIONAL result, "
          "not a finished physics claim.")
    header = ("{:>5}  {:>12}  {:>12}  {:>12}  {:>8}  {:>9}  {:>10}  "
              "{:>10}  {:>10}  {:>9}").format(
        "V", "lam_dyn", "lam_flat", "lam_static", "ratio", "flat_shift",
        "cond_spin", "cond_chg", "herm_rel", "cond_tol")
    print(header)
    for V in V_GRID:
        p = points[V]
        assert np.isfinite(p["lambda"])
        assert abs(p["lambda_imag"]) <= max(1.0e-8, 1.0e-6 * abs(p["lambda"])), (
            "V={}: leading eigenvalue has a non-negligible imaginary part "
            "{}".format(V, p["lambda_imag"]))
        assert p["green_sym_k"] < 1.0e-10 and p["green_sym_w"] < 1.0e-10
        assert p["kernel_hermiticity_relative"] < 1.0e-8

        # COND_TOL_DEFAULT-referenced (not tautological ">0") checks (review
        # MINOR): V<=1.0 must clear the PRODUCTION floor outright (that is
        # exactly why they need no override); V=1.2 is documented to sit
        # BELOW it (that is why it needs bond_cond_tol relaxed).
        assert p["cond_min_spin"] > COND_TOL_DEFAULT
        if V in COND_TOL_OVERRIDE:
            assert 0.0 < p["cond_min_charge"] < COND_TOL_DEFAULT, (
                "V={} was expected to sit BELOW the default conditioning "
                "floor (documented exception); it no longer does -- "
                "re-check whether the bond_cond_tol override is still "
                "needed.".format(V))
        else:
            assert p["cond_min_charge"] > COND_TOL_DEFAULT

        ratio = lam_static[V] / p["lambda"] if p["lambda"] else float("nan")
        assert np.isfinite(ratio)
        flat_shift = ((p["lambda_shell_flat"] - p["lambda"]) / p["lambda"]
                     if p["lambda"] else float("nan"))
        window_flag = ("" if 1.5 <= ratio <= 3.5 else
                       "  ** outside [1.5,3.5] window (indicator only) **")
        # near_instability is about the CONDITIONING score, not the ratio
        # window -- V=0.0 is nowhere near instability even though its ratio
        # is also outside the window (review MINOR: fixes the previous,
        # misleading "past_instability" label that fired at every V).
        if p["cond_min_charge"] < 0.01:
            window_flag += "  ** near_instability (cond_min_charge<0.01) **"
        print("{:5.2f}  {:12.6f}  {:12.6f}  {:12.6f}  {:8.3f}  {:+9.1%}  "
              "{:10.3e}  {:10.3e}  {:10.3e}  {:9.1e}{}".format(
                  V, p["lambda"], p["lambda_shell_flat"], lam_static[V],
                  ratio, flat_shift, p["cond_min_spin"], p["cond_min_charge"],
                  p["kernel_hermiticity_relative"], p["cond_tol_used"],
                  window_flag))
    # spec S3.6 doesn't require the ratio inside the window (indicator
    # only); the printed table + the module docstring's single-root-cause
    # section are the recorded disposition.


@pytest.mark.slow
@pytest.mark.parametrize("V", V_GRID)
def test_milestone_reproducible_from_solve_dynamic_entry_point(request, V,
                                                                tmp_path):
    """Spec goal ("fed by an externally supplied Green function"): the
    dict/TOML entry point ``solve_dynamic`` reproduces the direct-call
    ``_dynamic_point`` lambda at EVERY V (including 1.2, via ``[eliashberg]
    bond_cond_tol`` -- review IMPORTANT-3, no more "below solve_dynamic's
    reach" carve-out), and the ``bond_diagnostics=true`` archive
    (conditioning maps, Hermiticity/tail diagnostics) exists with the keys
    the task brief requires.
    """
    _require_slow(request)
    kx, ky, kz = _kgrid(L)
    it = _interactions(V)
    inter_k = sc._build_interaction_k(kx, ky, kz, it, 1)
    bond_set = bc.resolve_interactions(it["CoulombInter"], np.eye(3), 1)
    green = _load_green(L, V)
    S0_q, C0_q = sc._build_bond_m0_blocks(bond_set, it, inter_k, 1, kx, ky, kz)
    S_bond, C_bond, Vpp_s, Vpp_t = bc.bare_bond_vertices(bond_set, S0_q, C0_q, 1)
    chi_bar_w = bc.bond_bubble_dynamic(green, bond_set, BETA)
    cond_tol = COND_TOL_OVERRIDE.get(V, COND_TOL_DEFAULT)
    direct_lam, _, _, _, _ = _solve_bond_dynamic(
        chi_bar_w, S_bond, C_bond, Vpp_s, Vpp_t, green, bond_set, cond_tol)

    output_dir = str(tmp_path / "out")
    inp = _dynamic_input_dict(V, output_dir,
                              cond_tol=COND_TOL_OVERRIDE.get(V))
    lam = ed.solve_dynamic(inp)
    assert lam == pytest.approx(direct_lam.real, rel=1.0e-8)

    meta = np.load(os.path.join(output_dir, "gap_dynamic.npz"),
                   allow_pickle=False)
    assert bool(meta["bond_channels"])
    assert int(meta["bond_n_channels"]) == 5
    assert str(meta["bond_green_source"]) == "file"
    assert str(meta["bond_frequency_grid"]) == "uniform"
    assert float(meta["bond_cond_min_spin"]) > 0.0
    assert float(meta["bond_cond_min_charge"]) > 0.0
    assert float(meta["bond_cond_tol"]) == pytest.approx(cond_tol)
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
def test_v1_2_default_bond_cond_tol_still_refuses_via_solve_dynamic(
        request, tmp_path):
    """Without an explicit ``[eliashberg] bond_cond_tol``, ``solve_dynamic``
    still refuses ``V=1.2`` at the production default floor -- the override
    in ``test_milestone_reproducible_from_solve_dynamic_entry_point`` is an
    explicit, documented opt-in (module docstring), not a silently loosened
    default."""
    _require_slow(request)
    output_dir = str(tmp_path / "out")
    inp = _dynamic_input_dict(1.2, output_dir)   # no bond_cond_tol override
    with pytest.raises(ValueError, match="charge instability"):
        ed.solve_dynamic(inp)


@pytest.mark.slow
def test_matsubara_basis_ir_is_broken_on_this_fixture(request, tmp_path):
    """Pins DEVIATION 1 (module docstring) as a documented regression: the
    IR-basis leading eigenvalue at V=0.0 is many orders of magnitude away
    from the trustworthy uniform-basis reference, and grows WORSE (not
    better) as ``ir_wmax`` is raised -- the signature of an ill-conditioned
    least-squares fit in ``_ir_compress``'s constant-retention branch, not a
    basis-truncation error that a larger window would fix. The FAST sibling
    at toy scale (``tests/test_sc_bond_dynamic.py::
    test_dynamic_bond_ir_is_catastrophically_broken_at_nmat_32``) catches
    the same defect on every default CI run; this one pins it at the
    milestone's own production scale.

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
