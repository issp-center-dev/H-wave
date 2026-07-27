"""Task 7: the IR-axis dynamic bond kernel, the tail estimator, and the
frozen IR round-trip contract.

Normative reference: docs/superpowers/specs/2026-07-27-dynamic-bond-channels-
design.md, section 2 (IR-coordinate Hermiticity threshold), section 3.3
(instantaneous per-basis definitions + the tail-estimator pipeline) and
section 3.5 (the IR round-trip contract, steps 1-5, with its two frozen
inequalities). Every number asserted below is defined there; none of them is
tuned here.

The IR path is only exercised when ``sparse-ir`` is installed (the same
import-safe skip guard as tests/test_eliashberg_ir.py: a module-level
``pytest.importorskip`` would turn a missing optional dependency into an
ERROR under ``unittest discover``, which CI runs).
"""
import numpy as np
import pytest

from hwave.solver import bond_channels as bc

try:
    import sparse_ir  # noqa: F401
    _HAVE_SPARSE_IR = True
except ImportError:
    _HAVE_SPARSE_IR = False

BETA = 5.0


# ---------------------------------------------------------------------------
# Step 1/2: the tail estimator (spec S3.3), on an analytic profile
# ---------------------------------------------------------------------------
def test_tail_estimate_on_exact_inverse_square():
    """The estimator's pipeline is fixed by the spec; on an exact ``1/w^2``
    profile the fit is exact, so the returned ``tail_est`` must reproduce the
    literal window-complement sum it models."""
    beta, nmat = 5.0, 64
    n_t = np.arange(nmat) - nmat // 2
    w = (2 * n_t + 1) * np.pi / beta
    X = (1.0 / w ** 2)[None, None, None, None, :]   # (1,1,1,1,nmat) toy
    est, rel, bad = bc.tail_estimate(X, beta, nmat)
    exact = (1.0 / beta) * sum(
        1.0 / (((2 * m + 1) * np.pi / beta) ** 2)
        for m in range(nmat // 2, 64 * nmat))       # positive tail
    exact *= 2                                      # both signs
    assert not bad
    assert abs(est - exact) / exact < 0.05
    # tail_est_rel is the same quantity over the stored-window sum of the
    # SAME profile (spec S3.3 "Relative form"), so it is checkable in closed
    # form too -- and it must be small here (the window holds most of the
    # 1/w^2 weight).
    denom = (1.0 / beta) * float(np.sum(1.0 / w ** 2))
    assert abs(rel - est / denom) <= 1e-12 * abs(rel)
    assert 0.0 < rel < 0.05


def test_tail_estimate_flags_an_unfittable_profile():
    """The reliability flag (relative outer-shell fit residual > 0.2) is the
    spec's guard against reporting a tail number for data the ``|a|/w^2 +
    |b|/w^3`` model does not describe. A profile that oscillates between
    neighbouring frequencies (an X nowhere near its asymptotic regime) is the
    canonical case; a merely noisy but smooth-on-average profile is NOT
    flagged, which is the estimator behaving as specified (the flag is a
    model-adequacy test at the 0.2 relative-residual level, not a
    noise detector)."""
    beta, nmat = 5.0, 64
    n_t = np.arange(nmat) - nmat // 2
    X = (1.0 + 0.9 * (-1.0) ** n_t)[None, :]
    _, _, bad = bc.tail_estimate(X, beta, nmat)
    assert bad
    rng = np.random.default_rng(11)
    _, _, mild = bc.tail_estimate(rng.random((2, 3, nmat)) + 0.5, beta, nmat)
    assert not mild


def test_tail_estimate_reduces_over_k_and_orbital_pairs():
    """The fitted object is the single scalar profile ``x(n) = max_{k,pairs}
    |X|`` (spec S3.3 "Data fitted"): one fit, never per-k. Broadcasting the
    same 1/w^2 profile over extra leading axes with a k-dependent (<=1)
    modulation must therefore give exactly the max-reduced answer."""
    beta, nmat = 5.0, 32
    n_t = np.arange(nmat) - nmat // 2
    w = (2 * n_t + 1) * np.pi / beta
    prof = 1.0 / w ** 2
    scale = np.array([0.25, 1.0, 0.5])[:, None]
    est_multi, rel_multi, _ = bc.tail_estimate(
        (scale * prof[None, :])[None, :, :], beta, nmat)
    est_one, rel_one, _ = bc.tail_estimate(prof[None, :], beta, nmat)
    assert abs(est_multi - est_one) <= 1e-12 * est_one
    # the denominator is the max-reduced profile's window sum as well
    assert abs(rel_multi - rel_one) <= 1e-12 * rel_one
