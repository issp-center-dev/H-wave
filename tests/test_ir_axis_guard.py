"""Regression tests for the IRAxis ill-conditioning guard (issue #153).

At small ``Lambda = beta * wmax`` the sparse-IR sampling node sets stop
supporting a well-posed least-squares fit onto the L basis coefficients, and
``IRAxis`` used to react in one of two ways, neither acceptable:

* ``beta = 0.02``, ``wmax = 5``, ``eps = 1e-8``, bosonic: sparse-ir returns a
  SINGLE Matsubara sampling point (n = 0) for a basis of size L = 4, so the
  design matrix is underdetermined and collapses to 1-D; ``np.linalg.pinv``
  then raised a bare ``LinAlgError`` ("1-dimensional array given") from deep
  inside the transform-matrix build, naming none of the parameters that
  caused it.
* ``beta <= 0.005``, ``wmax = 5``, ``eps = 1e-8``, fermionic: construction
  SUCCEEDED, but the Matsubara design matrix has condition number 1e6-1e7, so
  the pinv-fitted transform matrices are silently wrong -- fitting a plain
  pole Green's function through them and evaluating on tau lands ~4x off the
  exact answer with NO exception raised. That is the dangerous mode: a FLEX
  run with an unusually low ``ir_wmax`` at high temperature would produce
  plausible-looking wrong susceptibilities.

Both are now rejected at construction time by a single coefficient-space
roundtrip self-check on the fitted matrices (see
``IRAxis._check_fit_conditioning``), raising ``ValueError`` naming beta, wmax
and eps. These tests pin (a) the clean error at the crash point, (b) the same
clean error at the silent-wrongness point, (c) that the healthy regime the
shipped fixtures use still constructs AND that its transform matrices are
bit-identical to the unguarded recipe, and (d) that the module still skips
cleanly where sparse-ir is absent.
"""
import unittest

import numpy as np

# Import-safe under BOTH pytest and unittest discovery (the CI gate runs
# `unittest discover`, where a module-level pytest.importorskip would turn a
# missing optional dependency into an import ERROR instead of a skip).
try:
    import sparse_ir  # noqa: F401
    _HAVE_SPARSE_IR = True
except ImportError:
    _HAVE_SPARSE_IR = False


# The two reproduced failure points of issue #153, and the healthy reference.
WMAX = 5.0
EPS = 1e-8
BETA_CRASH = 0.02        # bosonic: n_freq (1) < L (4) -> used to LinAlgError
BETA_SILENT = 0.001      # fermionic: cond(design) = 4.5e7 -> silently wrong
BETA_HEALTHY = 0.1       # cond 2.5-3.1, roundtrip residual ~1e-15


@unittest.skipUnless(_HAVE_SPARSE_IR, "sparse-ir not installed")
class TestIRAxisConditioningGuard(unittest.TestCase):
    """The two #153 regimes must both raise one clear, parameterised error."""

    def _construct(self, beta, statistics):
        from hwave.solver.ir_axis import IRAxis
        return IRAxis(beta=beta, wmax=WMAX, eps=EPS, statistics=statistics)

    def _assert_actionable(self, ctx, beta):
        msg = str(ctx.exception)
        # The message must name the three parameters that caused it, so the
        # user can act without reading the traceback.
        for token in ("beta", "wmax", "eps"):
            self.assertIn(token, msg,
                          "guard message must name {!r}: {}".format(token, msg))
        self.assertIn(repr(beta), msg.replace("=", "=").replace(",", ","),
                      "guard message must carry the offending beta value: "
                      + msg)
        self.assertIn(repr(WMAX), msg)
        self.assertIn(repr(EPS), msg)
        # ... and the remedy.
        self.assertIn("ir_wmax", msg)

    def test_crash_regime_raises_valueerror_not_linalgerror(self):
        """beta=0.02 (bosonic): underdetermined fit -> ValueError, and
        specifically NOT a bare LinAlgError leaking out of pinv."""
        with self.assertRaises(ValueError) as ctx:
            self._construct(BETA_CRASH, "B")
        self.assertNotIsInstance(ctx.exception, np.linalg.LinAlgError)
        self._assert_actionable(ctx, BETA_CRASH)
        # This regime is structurally underdetermined; say so.
        self.assertIn("4", str(ctx.exception))  # L = 4 appears in the message

    def test_silent_wrongness_regime_raises_valueerror(self):
        """beta=0.001 (fermionic): construction used to SUCCEED with ~4x-wrong
        transforms. It must now fail loudly with the same error class."""
        with self.assertRaises(ValueError) as ctx:
            self._construct(BETA_SILENT, "F")
        self.assertNotIsInstance(ctx.exception, np.linalg.LinAlgError)
        self._assert_actionable(ctx, BETA_SILENT)

    def test_silent_wrongness_boundary_still_rejected(self):
        """beta=0.005 is the LEAST ill-conditioned point that still fits a
        pole Green's function ~110% wrong; it must be rejected too (the guard
        threshold is not tuned so tightly that it only catches beta=0.001)."""
        with self.assertRaises(ValueError):
            self._construct(0.005, "F")

    def test_guard_error_is_not_the_missing_dependency_error(self):
        """A conditioning failure must not be reported as ImportError."""
        from hwave.solver.ir_axis import IRAxis
        with self.assertRaises(ValueError):
            IRAxis(beta=BETA_SILENT, wmax=WMAX, eps=EPS, statistics="F")


@unittest.skipUnless(_HAVE_SPARSE_IR, "sparse-ir not installed")
class TestHealthyPathUnchanged(unittest.TestCase):
    """The guard must be a pure check: healthy axes construct as before, with
    bit-identical transform matrices."""

    def test_healthy_axes_construct_and_roundtrip(self):
        from hwave.solver.ir_axis import IRAxis
        from hwave.solver import ir_axis as _mod
        for statistics in ("F", "B"):
            ax = IRAxis(beta=BETA_HEALTHY, wmax=WMAX, eps=EPS,
                        statistics=statistics)
            eye = np.eye(ax.L)
            for tag in ("freq", "tau"):
                resid = np.max(np.abs(
                    ax._m["eval_" + tag] @ ax._m["fit_" + tag] - eye))
                # measured at this point: 6.7e-16 .. 9.6e-16
                self.assertLess(resid, 1e-13,
                                "{} roundtrip residual {:.3e}"
                                .format(tag, resid))
                # ... and comfortably under the guard threshold (>= 3 decades)
                self.assertLess(resid * 1e3, _mod.IR_FIT_RESIDUAL_MAX)

    def test_shipped_fixture_regime_constructs(self):
        """The parameters the existing IR unit tests use (beta=50, wmax=8)
        must be unaffected by the guard."""
        from hwave.solver.ir_axis import IRAxis
        for statistics in ("F", "B"):
            ax = IRAxis(beta=50.0, wmax=8.0, eps=EPS, statistics=statistics)
            self.assertGreaterEqual(ax.n_freq, ax.L)
            self.assertGreaterEqual(ax.n_tau, ax.L)

    def test_transform_matrices_bit_identical_to_unguarded_recipe(self):
        """Pin that the guard changed no arithmetic: rebuild the four base
        matrices with the pre-change recipe (pinv of the raw design matrices)
        and require BITWISE equality with what IRAxis stores."""
        from hwave.solver.ir_axis import IRAxis
        for statistics in ("F", "B"):
            ax = IRAxis(beta=BETA_HEALTHY, wmax=WMAX, eps=EPS,
                        statistics=statistics)
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


@unittest.skipUnless(_HAVE_SPARSE_IR, "sparse-ir not installed")
class TestGuardThresholdMargins(unittest.TestCase):
    """The threshold is a measured constant; pin the margins it was chosen
    for so a later retune cannot silently reopen either failure mode."""

    def test_threshold_value_and_separation(self):
        from hwave.solver import ir_axis as _mod
        thr = _mod.IR_FIT_RESIDUAL_MAX
        # Largest residual measured over healthy axes (Lambda 0.5 .. 5e4,
        # eps 1e-6 .. 1e-12) was 6.3e-14; smallest residual of a MEASURABLY
        # WRONG axis (beta=0.005, wmax=5, eps=1e-8) was 9.7e-11.
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
                _mod.IRAxis(beta=BETA_HEALTHY, wmax=WMAX, eps=EPS,
                            statistics="F")
            self.assertIn("sparse-ir", str(ctx.exception))
        finally:
            _mod._import_sparse_ir = saved


if __name__ == "__main__":
    unittest.main()
