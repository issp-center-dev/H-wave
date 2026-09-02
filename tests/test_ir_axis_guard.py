"""Regression tests for the IRAxis ill-conditioning guard (issue #153).

At small ``Lambda = beta * wmax`` the sparse-IR sampling node sets stop
supporting a well-posed least-squares fit onto the L basis coefficients, and
``IRAxis`` used to react in one of two unacceptable ways: a bare
``LinAlgError`` out of ``np.linalg.pinv`` when the node set fell below L
points, or -- far worse -- a SILENT success whose transform matrices were
wrong by factors of order one. A FLEX run with an unusually low ``ir_wmax``
at high temperature would then have produced plausible-looking wrong
susceptibilities. Both are now rejected at construction with a ``ValueError``
naming beta, wmax and eps.

WHY THESE TESTS ARE WRITTEN AGAINST THE CONTRACT, NOT AGAINST FIXED
PARAMETER POINTS
------------------------------------------------------------------
The first version of this file pinned specific (beta, wmax, eps) points to
specific regimes. That was wrong: **which regime a given parameter point
lands in is sparse-ir-version and platform dependent, in both directions.**
Measured:

* local -- macOS arm64 / Accelerate, numpy 2.4.6, sparse-ir 2.1.1:
  beta=0.1, wmax=5 (Lambda=0.5) is healthy (n_tau=6 >= L=6); beta=0.001,
  wmax=5 is ill-conditioned (cond 4.5e7) and rejected.
* CI -- Linux / OpenBLAS, sparse-ir 1.1.7:
  beta=0.1, wmax=5 returns only TWO tau sampling nodes, so the structural
  check correctly fires and the point is NOT healthy; while beta=0.001
  constructs acceptably and is NOT rejected.
* reproduced locally against sparse-ir 1.1.7 (pinned in a throwaway venv;
  that version also needs scipy < 1.15 for scipy.linalg.interpolative.seed):
  beta=0.1, wmax=5 is REJECTED for both statistics, confirming the CI
  report; and -- in the opposite direction -- beta=0.02, wmax=5, BOSONIC
  BUILDS cleanly (L=4, n_freq=5) where sparse-ir 2.1.1 returns a single
  Matsubara node and the axis is rejected. The mapping moves both ways, so
  no fixed point can be pinned to a regime. beta=50, wmax=8 is healthy in
  BOTH (L=36, n_freq=36/37, n_tau=36, residuals 2e-15 .. 6e-15), which is
  why the fixture lives there.

Note for anyone reading warnings from an old sparse-ir: 1.1.7 emits
divide-by-zero / overflow / invalid RuntimeWarnings from its OWN
``poly.py`` while evaluating the basis, on every point including ones the
guard then rejects. That noise is library-internal and predates this guard.

The guard behaved correctly in both environments -- only the mapping from
parameters to regimes moved. So the fixed-point tests were replaced by:

1. a healthy fixture at comfortably large Lambda -- the same regime the
   shipped IR tests use (``tests/test_ir_axis.py`` samples beta=50, wmax=8,
   Lambda=400; all 25 shipped IR tests pass on both environments) -- which
   additionally ASSERTS it is healthy in the current environment before
   relying on it, so a future sampling change yields a diagnosis rather than
   a mystery;
2. an environment-adaptive ladder scan that pins the issue's actual
   contract -- a clean error somewhere on the descending-beta ladder, and
   NEVER a silent success with a bad residual anywhere on it;
3. mock- and direct-call-based tests for the raise machinery and message
   content, which depend on no sampling behaviour at all.
"""
import unittest
from unittest import mock

import numpy as np

# Import-safe under BOTH pytest and unittest discovery (the CI gate runs
# `unittest discover`, where a module-level pytest.importorskip would turn a
# missing optional dependency into an import ERROR instead of a skip).
try:
    import sparse_ir  # noqa: F401
    _HAVE_SPARSE_IR = True
except ImportError:
    _HAVE_SPARSE_IR = False

# Healthy fixture: the regime tests/test_ir_axis.py itself uses (Lambda=400).
# Every supported sparse-ir version samples this adequately.
HEALTHY_BETA = 50.0
HEALTHY_WMAX = 8.0
EPS = 1e-8

# The descending-beta ladder for the adaptive scan, at wmax=5, eps=1e-8.
# Provenance -- measured locally (macOS, sparse-ir 2.1.1), where 6 of the 14
# (beta, statistics) points raise and none escape:
#   beta    F                       B
#   0.05    built  resid 2.6e-16    built  resid 4.4e-16
#   0.02    built  resid 5.9e-16    RAISE  structural (1 Matsubara node, L=4)
#   0.01    built  resid 1.6e-15    RAISE  structural
#   0.005   RAISE  ill-cond 9.7e-11 built  resid 8.9e-16
#   0.002   RAISE  ill-cond 3.7e-10 built  resid 2.8e-15
#   0.001   RAISE  ill-cond 9.9e-09 built  resid 1.2e-15
#   0.0005  RAISE  ill-cond         built  resid 3.8e-16
# sparse-ir 1.1.7 distributes RAISE/built differently across the same ladder
# (measured: 4 rejections / 10 builds, versus 6 / 8 on 2.1.1) and inverts
# several individual points -- but in BOTH versions at least one point is
# rejected and NO built point escapes the residual threshold, which is
# exactly what the two assertions below require.
LADDER_BETAS = (0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005)
LADDER_WMAX = 5.0


def _roundtrip_residuals(ax):
    """The two residuals the guard itself checks, recomputed from the stored
    transform matrices: max|eval @ fit - I| in coefficient space."""
    eye = np.eye(ax.L)
    return {tag: float(np.max(np.abs(
        ax._m["eval_" + tag] @ ax._m["fit_" + tag] - eye)))
        for tag in ("freq", "tau")}


@unittest.skipUnless(_HAVE_SPARSE_IR, "sparse-ir not installed")
class TestGuardContractOnBetaLadder(unittest.TestCase):
    """The issue's contract, asserted without assuming WHERE the regime
    boundary sits in this environment: descending beta at fixed wmax/eps must
    produce a clean, contextual error somewhere, and must NEVER produce a
    silently-wrong axis anywhere."""

    def _scan(self):
        from hwave.solver.ir_axis import IRAxis
        rejections, built = [], []
        for beta in LADDER_BETAS:
            for statistics in ("F", "B"):
                try:
                    ax = IRAxis(beta=beta, wmax=LADDER_WMAX, eps=EPS,
                                statistics=statistics)
                except ValueError as exc:
                    rejections.append((beta, statistics, exc))
                else:
                    built.append((beta, statistics, ax,
                                  _roundtrip_residuals(ax)))
        return rejections, built

    def test_ladder_produces_at_least_one_contextual_rejection(self):
        """(a) The guard must actually engage somewhere on the ladder, and
        every rejection it issues must be the contextual ValueError -- not a
        bare LinAlgError (which, note, IS a ValueError subclass)."""
        rejections, built = self._scan()
        self.assertTrue(
            rejections,
            "no point on the descending-beta ladder {} (wmax={}, eps={}) was "
            "rejected; the guard never engages in this environment"
            .format(LADDER_BETAS, LADDER_WMAX, EPS))
        for beta, statistics, exc in rejections:
            where = "beta={}, statistics={}".format(beta, statistics)
            self.assertNotIsInstance(exc, np.linalg.LinAlgError, where)
            msg = str(exc)
            for token in ("beta", "wmax", "eps", "ir_wmax"):
                self.assertIn(token, msg,
                              "{}: message must name {!r}: {}"
                              .format(where, token, msg))
            self.assertIn(repr(beta), msg, where)
            self.assertIn(repr(LADDER_WMAX), msg, where)
            self.assertIn(repr(EPS), msg, where)

    def test_no_ladder_point_constructs_silently_wrong(self):
        """(b) THE central pin of #153. Every point that DOES construct must
        have transform matrices whose measured roundtrip residual is within
        the guard's own threshold -- i.e. the guard never lets a
        silently-wrong axis through, wherever this environment happens to put
        the boundary."""
        from hwave.solver import ir_axis as _mod
        rejections, built = self._scan()
        self.assertTrue(built, "the whole ladder was rejected; the scan "
                               "proves nothing about silent wrongness")
        # NaN-SAFE, mirroring the guard's own predicate: `nan > thr` is
        # False, so a plain `>` comparison would wave a NaN residual through
        # -- the worst silent wrongness of all.
        escapes = [(beta, statistics, res)
                   for beta, statistics, _ax, res in built
                   if not (max(res.values()) <= _mod.IR_FIT_RESIDUAL_MAX)]
        self.assertEqual(
            escapes, [],
            "IRAxis constructed successfully with a roundtrip residual above "
            "IR_FIT_RESIDUAL_MAX={:.1e} (or NaN) -- this is exactly the "
            "silent wrongness #153 exists to prevent: {}"
            .format(_mod.IR_FIT_RESIDUAL_MAX, escapes))
        # and no accepted axis may carry inf/NaN in the matrices the
        # transforms actually use, composites included
        nonfinite = [(beta, statistics,
                      sorted(k for k, v in ax._m.items()
                             if not np.all(np.isfinite(v))))
                     for beta, statistics, ax, _res in built]
        nonfinite = [row for row in nonfinite if row[2]]
        self.assertEqual(nonfinite, [],
                         "accepted axes carry non-finite transform matrices: "
                         "{}".format(nonfinite))


@unittest.skipUnless(_HAVE_SPARSE_IR, "sparse-ir not installed")
class TestGuardFailureRoutes(unittest.TestCase):
    """The raise machinery and message content, pinned WITHOUT depending on
    sparse-ir's sampling behaviour: the failures are forced at a healthy
    parameter point, or by calling the guard directly.

    ``np.linalg.pinv`` is patched on the numpy module itself, which is what
    ``ir_axis`` resolves at call time. That is safe around the whole
    construction: sparse-ir's basis build makes no pinv call -- the side
    effects below record every call they see, and each test asserts it
    intercepted the matrix it meant to.
    """

    def _msg_is_actionable(self, msg):
        for token in ("beta", "wmax", "eps", "ir_wmax"):
            self.assertIn(token, msg, "message must name {!r}: {}"
                          .format(token, msg))
        self.assertIn(repr(HEALTHY_BETA), msg)
        self.assertIn(repr(HEALTHY_WMAX), msg)
        self.assertIn(repr(EPS), msg)

    def test_pinv_linalgerror_is_reraised_with_context(self):
        """A LinAlgError out of pinv itself (e.g. SVD non-convergence on a
        well-shaped matrix) must surface as the contextual ValueError, never
        raw. Both statistics."""
        from hwave.solver.ir_axis import IRAxis
        for statistics in ("F", "B"):
            seen = []

            def _boom(a, *args, **kwargs):
                seen.append(np.shape(a))
                raise np.linalg.LinAlgError("SVD did not converge")

            with mock.patch("numpy.linalg.pinv", side_effect=_boom):
                with self.assertRaises(ValueError) as ctx:
                    IRAxis(beta=HEALTHY_BETA, wmax=HEALTHY_WMAX, eps=EPS,
                           statistics=statistics)
            self.assertTrue(seen, "pinv was never reached")
            self.assertNotIsInstance(ctx.exception, np.linalg.LinAlgError)
            msg = str(ctx.exception)
            self._msg_is_actionable(msg)
            # the underlying failure is quoted, not swallowed
            self.assertIn("SVD did not converge", msg)
            self.assertIn("pseudo-inverting", msg)

    def test_bad_residual_on_the_tau_matrix_is_reported_as_tau(self):
        """Force a bad roundtrip residual on the TAU design matrix only, at
        an otherwise healthy parameter point, and check the error names the
        tau axis (not the Matsubara one) and reports the residual. The two
        matrices are told apart by dtype: uhat() is complex, u() is real."""
        from hwave.solver.ir_axis import IRAxis
        from hwave.solver import ir_axis as _mod
        real_pinv = np.linalg.pinv
        for statistics in ("F", "B"):
            perturbed = []

            def _spoil(a, *args, **kwargs):
                out = real_pinv(a, *args, **kwargs)
                if not np.iscomplexobj(a):      # the tau matrix
                    perturbed.append(np.shape(a))
                    out = out + 1.0e-6          # >> IR_FIT_RESIDUAL_MAX
                return out

            with mock.patch("numpy.linalg.pinv", side_effect=_spoil):
                with self.assertRaises(ValueError) as ctx:
                    IRAxis(beta=HEALTHY_BETA, wmax=HEALTHY_WMAX, eps=EPS,
                           statistics=statistics)
            self.assertTrue(perturbed, "the tau matrix was never perturbed")
            self.assertNotIsInstance(ctx.exception, np.linalg.LinAlgError)
            msg = str(ctx.exception)
            self._msg_is_actionable(msg)
            self.assertIn("ill-conditioned", msg)
            self.assertIn("tau sampling matrix", msg)
            self.assertNotIn("Matsubara sampling matrix", msg)
            self.assertIn("roundtrip residual", msg)
            self.assertIn("{:.1e}".format(_mod.IR_FIT_RESIDUAL_MAX), msg)

    def test_structural_check_rejects_short_matrix_for_either_label(self):
        """The n_nodes < L branch, exercised by calling the guard directly
        with a deliberately short design matrix. Deterministic: it depends on
        no sparse-ir sampling behaviour, and covers BOTH axis labels (no real
        parameter point can be relied on to make tau the failing axis)."""
        from hwave.solver.ir_axis import IRAxis
        ax = IRAxis(beta=HEALTHY_BETA, wmax=HEALTHY_WMAX, eps=EPS,
                    statistics="F")
        for axis_name in ("Matsubara", "tau"):
            short = np.zeros((ax.L - 1, ax.L))
            with self.assertRaises(ValueError) as ctx:
                ax._pinv_checked(short, axis_name)
            self.assertNotIsInstance(ctx.exception, np.linalg.LinAlgError)
            msg = str(ctx.exception)
            self.assertIn("underdetermined", msg)
            self.assertIn("{} sampling node".format(axis_name), msg)
            self.assertIn("L={}".format(ax.L), msg)
            self._msg_is_actionable(msg)
        # a matrix with exactly L rows is NOT short -- boundary is >=, not >
        square = np.eye(ax.L)
        self.assertIsNone(np.testing.assert_allclose(
            ax._pinv_checked(square, "tau"), np.eye(ax.L), atol=1e-12))

    def test_condition_number_diagnostic_cannot_mask_the_guard_error(self):
        """The condition number in the ill-conditioning message runs a SECOND
        SVD while an error is already being raised. If that diagnostic SVD
        fails it must degrade to 'unavailable', not replace the contextual
        ValueError with a bare LinAlgError."""
        from hwave.solver.ir_axis import IRAxis
        real_pinv = np.linalg.pinv

        def _spoil(a, *args, **kwargs):
            return real_pinv(a, *args, **kwargs) + 1.0e-6

        def _cond_boom(a, *args, **kwargs):
            raise np.linalg.LinAlgError("SVD did not converge")

        with mock.patch("numpy.linalg.pinv", side_effect=_spoil), \
                mock.patch("numpy.linalg.cond", side_effect=_cond_boom):
            with self.assertRaises(ValueError) as ctx:
                IRAxis(beta=HEALTHY_BETA, wmax=HEALTHY_WMAX, eps=EPS,
                       statistics="F")
        self.assertNotIsInstance(ctx.exception, np.linalg.LinAlgError)
        msg = str(ctx.exception)
        self._msg_is_actionable(msg)
        self.assertIn("condition number unavailable", msg)
        self.assertIn("roundtrip residual", msg)

    def test_cond_str_helper_degrades_on_nonfinite_input(self):
        """Direct unit check of the diagnostic helper."""
        from hwave.solver.ir_axis import _cond_str
        self.assertIn("e+", _cond_str(np.eye(3) * np.array([1.0, 2.0, 4.0])))
        bad = np.array([[np.nan, 0.0], [0.0, 1.0]])
        self.assertEqual(_cond_str(bad), "unavailable")


@unittest.skipUnless(_HAVE_SPARSE_IR, "sparse-ir not installed")
class TestHealthyPathUnchanged(unittest.TestCase):
    """The guard must be a pure check: healthy axes construct as before, with
    bit-identical transform matrices."""

    def _healthy_axis(self, statistics):
        """Construct the fixture axis AND assert it really is healthy in this
        environment, so a sparse-ir sampling change is diagnosed rather than
        showing up as a baffling failure in the assertions below."""
        from hwave.solver.ir_axis import IRAxis
        try:
            ax = IRAxis(beta=HEALTHY_BETA, wmax=HEALTHY_WMAX, eps=EPS,
                        statistics=statistics)
        except ValueError as exc:
            self.fail(
                "the healthy fixture (beta={}, wmax={}, eps={}, Lambda={}, "
                "statistics={!r}) was REJECTED by the guard in this "
                "environment -- this sparse-ir build samples it differently "
                "than expected, so the fixture must move to a larger "
                "Lambda, not the guard be loosened. Guard said: {}"
                .format(HEALTHY_BETA, HEALTHY_WMAX, EPS,
                        HEALTHY_BETA * HEALTHY_WMAX, statistics, exc))
        self.assertGreaterEqual(
            ax.n_freq, ax.L,
            "fixture has fewer Matsubara nodes ({}) than basis functions "
            "({}) in this environment".format(ax.n_freq, ax.L))
        self.assertGreaterEqual(
            ax.n_tau, ax.L,
            "fixture has fewer tau nodes ({}) than basis functions ({}) in "
            "this environment".format(ax.n_tau, ax.L))
        # all four base matrices carry the node counts they should
        self.assertEqual(ax._m["fit_freq"].shape, (ax.n_freq, ax.L))
        self.assertEqual(ax._m["eval_freq"].shape, (ax.L, ax.n_freq))
        self.assertEqual(ax._m["fit_tau"].shape, (ax.n_tau, ax.L))
        self.assertEqual(ax._m["eval_tau"].shape, (ax.L, ax.n_tau))
        return ax

    def test_healthy_axes_construct_and_roundtrip(self):
        from hwave.solver import ir_axis as _mod
        for statistics in ("F", "B"):
            ax = self._healthy_axis(statistics)
            for tag, resid in _roundtrip_residuals(ax).items():
                # comfortably under the guard threshold in any environment
                self.assertLess(
                    resid, _mod.IR_FIT_RESIDUAL_MAX * 1e-2,
                    "{} roundtrip residual {:.3e} is within 2 decades of the "
                    "guard threshold; the healthy fixture is too marginal"
                    .format(tag, resid))

    def test_transform_matrices_bit_identical_to_unguarded_recipe(self):
        """Pin that the guard changed no arithmetic: rebuild the four base
        matrices with the pre-change recipe (pinv of the raw design matrices)
        and require BITWISE equality with what IRAxis stores."""
        for statistics in ("F", "B"):
            ax = self._healthy_axis(statistics)
            basis = ax._basis
            eval_freq = basis.uhat(ax.freq_n).T
            eval_tau = basis.u(ax.tau).T
            ref = {
                "fit_freq": np.ascontiguousarray(
                    np.linalg.pinv(eval_freq).T),
                "eval_freq": np.ascontiguousarray(eval_freq.T),
                "fit_tau": np.ascontiguousarray(np.linalg.pinv(eval_tau).T),
                "eval_tau": np.ascontiguousarray(eval_tau.T),
            }
            ref["freq_to_tau"] = np.ascontiguousarray(
                ref["fit_freq"] @ ref["eval_tau"])
            ref["tau_to_freq"] = np.ascontiguousarray(
                ref["fit_tau"] @ ref["eval_freq"])
            for name, want in ref.items():
                got = ax._m[name]
                self.assertEqual(got.shape, want.shape, name)
                self.assertTrue(np.array_equal(got, want),
                                "{} ({}) is not bit-identical to the "
                                "unguarded construction".format(
                                    name, statistics))


class TestGuardThresholdMargins(unittest.TestCase):
    """The threshold is a measured constant; pin the margins it was chosen
    for so a later retune cannot silently reopen either failure mode. Needs
    no sparse-ir (it only reads the constant)."""

    def test_threshold_value_and_separation(self):
        from hwave.solver import ir_axis as _mod
        thr = _mod.IR_FIT_RESIDUAL_MAX
        # Largest residual measured over healthy axes (Lambda 0.5 .. 5e4,
        # eps 1e-6 .. 1e-12) was 6.3e-14; smallest residual of a MEASURABLY
        # WRONG axis (beta=0.005, wmax=5, eps=1e-8, macOS) was 9.7e-11.
        self.assertGreater(thr, 6.3e-14 * 10.0,
                           "threshold must clear the healthy maximum by >10x")
        self.assertLess(thr, 9.7e-11,
                        "threshold must still reject the beta=0.005 axis")


class TestSkipsWithoutSparseIR(unittest.TestCase):
    """Mirrors the skip discipline of tests/test_ir_axis.py: with sparse-ir
    absent, IRAxis raises the actionable ImportError rather than anything the
    guard produces."""

    def test_missing_sparse_ir_still_raises_importerror(self):
        from hwave.solver import ir_axis as _mod

        def _boom():
            raise ImportError("no sparse_ir")

        saved = _mod._import_sparse_ir
        _mod._import_sparse_ir = _boom
        try:
            with self.assertRaises(ImportError) as ctx:
                _mod.IRAxis(beta=HEALTHY_BETA, wmax=HEALTHY_WMAX, eps=EPS,
                            statistics="F")
            self.assertIn("sparse-ir", str(ctx.exception))
        finally:
            _mod._import_sparse_ir = saved


if __name__ == "__main__":
    unittest.main()
