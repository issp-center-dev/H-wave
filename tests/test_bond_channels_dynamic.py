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
