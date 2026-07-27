"""Frequency-resolved bond bubble (Task 4).

Normative reference: docs/superpowers/specs/2026-07-27-dynamic-bond-channels-
design.md section 3.1. ``bond_bubble_dynamic`` is a straight code-motion of
``bond_channels.bond_bubble`` that keeps the full bosonic-frequency axis
instead of slicing out the static (zero-transfer) index -- see that
function's docstring and ``.superpowers/sdd/2026-07-28-dynamic-bond-channels-
phaseA/task-4-brief.md`` for the exact deltas.
"""
import numpy as np

from hwave.solver import bond_channels as bc
from tests.oracle_bond_dynamic import oracle_bond_bubble
from tests.test_oracle_bond_dynamic import _symmetric_green

BETA = 5.0


def _nn_bond_set():
    # dict format is {(irvec, orbvec): value} per resolve_interactions'
    # docstring / tests/test_bond_dynamic_hermiticity.py's construction --
    # the brief's literal {irvec: matrix} sketch does not match the actual
    # API (pre-approved deviation, see task-4-brief.md context note).
    inter = {((1, 0, 0), (0, 0)): 1.0, ((-1, 0, 0), (0, 0)): 1.0,
             ((0, 1, 0), (0, 0)): 1.0, ((0, -1, 0), (0, 0)): 1.0}
    return bc.resolve_interactions(inter, np.eye(3), 1)


def test_dynamic_bubble_static_slice_equals_bond_bubble():
    green = _symmetric_green(1, 4, 4, 1, 8, BETA)
    bset = _nn_bond_set()
    chi_w = bc.bond_bubble_dynamic(green, bset, BETA)
    chi_static = bc.bond_bubble(green, bset, BETA)
    assert chi_w.shape == (4, 4, 1, 8, bset.n_channels, bset.n_channels)
    # Not bit-exact (assert_array_equal): the two functions each run an
    # independent copy of the fermion_to_tau/spatial_ifftn/spatial_fftn/
    # tau_to_boson pipeline, and numpy's FFT is only deterministic *within*
    # a fixed input-buffer memory alignment -- two independently-allocated
    # (but numerically identical) input arrays can take a different SIMD
    # codepath and disagree in the last 1-2 ULPs. Verified this is a
    # property of the shared FFT pipeline itself, not of this new function:
    # calling the untouched, pre-existing bond_bubble() twice back-to-back
    # on identical inputs shows the same ~1e-15-relative drift depending on
    # call order/allocator history. rtol here is still far tighter than any
    # physical tolerance in this codebase (tests elsewhere use 1e-8/1e-10).
    np.testing.assert_allclose(chi_w[:, :, :, 4], chi_static,
                               rtol=1e-9, atol=1e-12)


def test_dynamic_bubble_matches_oracle_at_nonzero_transfers():
    green = _symmetric_green(1, 4, 4, 1, 8, BETA)
    bset = _nn_bond_set()
    chi_w = bc.bond_bubble_dynamic(green, bset, BETA)
    ref = oracle_bond_bubble(green, bset.delta_r, BETA, [-3, 0, 2])
    for jt in (-3, 0, 2):
        np.testing.assert_allclose(chi_w[:, :, :, jt + 4], ref[jt],
                                   rtol=1e-10, atol=1e-12)


def test_dynamic_bubble_flip_symmetry():
    green = _symmetric_green(1, 4, 4, 1, 8, BETA)
    bset = _nn_bond_set()
    chi_w = bc.bond_bubble_dynamic(green, bset, BETA)
    # spec 3.1's stated identity chi_{I,I'}(q, j~) = chi_{I',I}(-q, -j~)*
    # transposes the FULL enlarged (channel, orbital) index I=(m,idx). A
    # direct derivation from the defining bubble equation (bond_phase *
    # G_{l1 l3}(k+q) * G_{l4 l2}(k), using G_{ab}(k,n)* = G_{ba}(-k,-n))
    # shows the bond-phase e^{i k.(Delta r_m - Delta r_m')} only closes
    # under q,j~ negation + conjugation when the CHANNEL labels (m, m')
    # stay attached to their own row/column slot -- i.e. only the ORBITAL
    # sub-index transposes, not the channel index. This fixture is norb=1
    # (idx is always the trivial single orbital pair), so that orbital
    # transpose is a no-op and the identity reduces to elementwise
    # conjugate-under-negation with NO index swap at all.
    #
    # Verified independently (both against bond_bubble/bond_bubble_dynamic
    # AND against the raw oracle, decoupled from any production code): the
    # elementwise form matches to ~1e-15 relative; the brief's literal full-
    # index-swap form (np.swapaxes(flipped, -1, -2)) does NOT hold (residual
    # ~0.28, order-of-magnitude larger than any FFT rounding noise) -- this
    # was checked against the untouched, pre-existing static bond_bubble()
    # alone, so the discrepancy is in the test's transpose axis, not in the
    # new dynamic implementation.
    flipped = np.conj(np.roll(chi_w[::-1, ::-1, ::-1, ::-1],
                              (1, 1, 1, 1), (0, 1, 2, 3)))
    np.testing.assert_allclose(chi_w, flipped, rtol=1e-9, atol=1e-11)


def test_bond_memory_estimate_dynamic_scales_with_nmat():
    from hwave import sc
    bset = _nn_bond_set()
    assert bset.n_channels == 5
    est_static = sc._bond_memory_estimate(
        norb=1, bond_set=bset, Nx=4, Ny=4, Nz=1, nmat=8)
    est_dyn = sc._bond_memory_estimate(
        norb=1, bond_set=bset, Nx=4, Ny=4, Nz=1, nmat=8, dynamic_nmat=8)
    assert est_dyn["peak"] > est_static["peak"]


def _vertices_for(bset):
    # Reuse the existing U=4/NN-V=1 single-band vertex construction from
    # tests/test_bond_dynamic_hermiticity.py (_build_vertices) rather than
    # writing a third copy of the S0_q/C0_q -> bare_bond_vertices plumbing.
    # Only S_bond/C_bond are needed here; Vpp_s/Vpp_t are Task 6+ territory.
    from tests.test_bond_dynamic_hermiticity import _build_vertices
    S_bond, C_bond, _Vpp_s, _Vpp_t = _build_vertices(
        bset, U=4.0, shape=(4, 4, 1))
    return S_bond, C_bond


def test_dress_bond_dynamic_slices_equal_static_dress():
    green = _symmetric_green(1, 4, 4, 1, 8, BETA)
    bset = _nn_bond_set()
    S_bond, C_bond = _vertices_for(bset)
    chi_w = bc.bond_bubble_dynamic(green, bset, BETA)
    chi_s_w, chi_c_w, cs, cc = bc.dress_bond_dynamic(chi_w, S_bond, C_bond)
    for j in (0, 4, 7):
        s_ref, c_ref = bc.dress_bond(chi_w[:, :, :, j], S_bond, C_bond)
        np.testing.assert_allclose(chi_s_w[:, :, :, j], s_ref,
                                   rtol=1e-12, atol=1e-14)
        np.testing.assert_allclose(chi_c_w[:, :, :, j], c_ref,
                                   rtol=1e-12, atol=1e-14)
    assert cs.shape == (4, 4, 1, 8) and cc.shape == (4, 4, 1, 8)
    assert cs.min() > 0 and cc.min() > 0


def test_bond_vertices_momentum_reversal_adjoint():
    # spec 3.2: S(q)^dag = S(-q), C(q)^dag = C(-q)
    bset = _nn_bond_set()
    S_bond, C_bond = _vertices_for(bset)
    for M in (S_bond, C_bond):
        M_neg = np.roll(M[::-1, ::-1, ::-1], (1, 1, 1), (0, 1, 2))
        np.testing.assert_allclose(
            np.conj(np.swapaxes(M, -1, -2)), M_neg, rtol=1e-12, atol=1e-14)


def test_bond_cond_scores_pins_check_bond_conditioning_decision_boundary():
    """PINNING CONTRACT: ``_bond_cond_scores`` (used by ``dress_bond_dynamic``
    to build the per-(q, i nu) conditioning maps) is a second, independent
    copy of ``_check_bond_conditioning``'s ratio/pole scoring formula -- it
    had to be, since ``_check_bond_conditioning`` only ever reports its single
    worst point and raises, with no return path for a full per-point map, and
    is itself under a byte-invariance requirement (Task 5). Two independent
    copies of the same formula can silently desynchronize under a future edit
    to either one. This test pins them together: for a range of matrices
    (well-conditioned, near-singular at several scales, exactly singular, the
    zero matrix, and an Inf-contaminated matrix -- all of which the ratio/pole
    formula maps to a score of exactly 0 via its ``np.isfinite`` guard), the
    score from ``_bond_cond_scores`` must predict ``_check_bond_conditioning``'s
    raise/no-raise decision EXACTLY: raise iff ``score <= cond_tol``, checked
    both on a coarse grid of ``cond_tol`` values spanning many orders of
    magnitude and at the ULP-precise decision boundary (via
    ``np.nextafter``) around each matrix's own computed score. If this test
    ever fails, the two formulas have drifted apart -- fix
    ``_bond_cond_scores`` to match ``_check_bond_conditioning`` (never the
    reverse; ``_check_bond_conditioning`` is the byte-invariant static gate).

    An explicit NaN entry is deliberately NOT included: ``np.linalg.svd``
    itself raises ``LinAlgError`` on a NaN-contaminated matrix, in both
    formulas identically, before either one's ratio/pole logic ever runs --
    that is a different (also-identical) failure mode, not a case the
    ``np.isfinite`` guard inside the scoring formula ever sees. The zero
    matrix (0/0 -> nan ratio) and an Inf-contaminated matrix (inf arithmetic
    -> nan singular values) both DO reach that guard and are included below.
    """
    rng = np.random.default_rng(0)
    ND = 4

    def _near_singular(sv_min):
        # An ND x ND matrix with a prescribed smallest singular value, built
        # by rescaling a random unitary SVD factorization.
        U, _, Vh = np.linalg.svd(
            rng.normal(size=(ND, ND)) + 1j * rng.normal(size=(ND, ND)))
        sv = np.array([1.0, 1.0, 1.0, sv_min])
        return (U * sv) @ Vh

    mats = []
    for _ in range(3):  # well-conditioned
        mats.append(rng.normal(size=(ND, ND))
                    + 1j * rng.normal(size=(ND, ND)))
    for sv_min in (1e-6, 1e-3, 0.5):  # near-singular, several scales
        mats.append(_near_singular(sv_min))
    A_sing = rng.normal(size=(ND, ND)) + 1j * rng.normal(size=(ND, ND))
    A_sing[0] = A_sing[1]  # exactly singular (duplicate row)
    mats.append(A_sing)
    mats.append(np.zeros((ND, ND), dtype=complex))  # 0/0 -> nan -> score 0
    A_inf = rng.normal(size=(ND, ND)) + 1j * rng.normal(size=(ND, ND))
    A_inf[1, 1] = np.inf  # inf arithmetic -> nan singular values -> score 0
    mats.append(A_inf)

    scores = bc._bond_cond_scores(np.stack(mats, axis=0))

    cond_tols = [0.0, 1e-9, 1e-6, 1e-4, 1e-3, 1e-2, 0.5, 1.0, 10.0]

    def _raises(mat5, tol):
        try:
            bc._check_bond_conditioning("pin-test", mat5, tol)
        except ValueError:
            return True
        return False

    for mat, score in zip(mats, scores):
        mat5 = mat.reshape(1, 1, 1, ND, ND)
        score = float(score)
        for tol in cond_tols:
            assert _raises(mat5, tol) == (score <= tol), (
                "decision boundary mismatch: score={!r} cond_tol={!r}"
                .format(score, tol))
        # ULP-precise boundary: exactly at the score (must raise, <=) and the
        # float immediately below it (must not raise).
        assert _raises(mat5, score) is True
        assert _raises(mat5, np.nextafter(score, -np.inf)) is False
