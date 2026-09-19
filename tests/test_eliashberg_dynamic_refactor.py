"""Golden files for the eigen-driver refactor: the on-site dynamic solver's
outputs must not change (spec 5.3). The golden values are recorded ONCE from
the pre-refactor code and stored under tests/eliashberg_dynamic_golden/."""
import os
import pathlib
import shutil
import tempfile
import unittest

import numpy as np

from tests.test_eliashberg_dynamic import _write_flex_fixture

_GOLDEN = os.path.join("tests", "eliashberg_dynamic_golden")
_CASES = (("iteration", None), ("iteration", 0.5), ("eigenvalue", None), ("both", None))


def _run(tmp, solver_mode, spectral_shift):
    import hwave.sc as sc
    # scipy's eigs()/eigsh() seed their default ARPACK start vector from a
    # fresh numpy.random.Generator built from OS entropy on every call
    # (its "rng" kwarg, default None) -- legacy np.random.seed() has no
    # effect on it. Without pinning this, two runs of the same
    # "eigenvalue"/"both" case diverge at the last bit (verified
    # empirically) and the golden comparison below is never
    # self-consistent. Route any entropy-seeded default_rng() call through
    # a fixed seed for the duration of this run only; an explicit seed
    # (e.g. this module's own fixture RNG) passes through unchanged.
    _real_default_rng = np.random.default_rng

    def _pinned_default_rng(seed=None):
        return _real_default_rng(20260918 if seed is None else seed)

    np.random.default_rng = _pinned_default_rng
    try:
        m = _write_flex_fixture(pathlib.Path(tmp), nmat=8, norb=1, Nx=2, Ny=2, Nz=1)
        inp = {"mode": {"param": {"T": 0.5, "CellShape": [2, 2, 1], "SubShape": [1, 1, 1],
                                  "Nmat": 8, "filling": 0.5}},
               "file": {"input": {"interaction": {"path_to_input": "tests/rpa/input",
                                                  "Geometry": "geom.dat",
                                                  "Transfer": "transfer.dat",
                                                  "CoulombIntra": "coulombintra.dat"}},
                        "output": {"path_to_output": str(tmp)}},
               "eliashberg": {"chi0q_mode": "flex", "frequency": "dynamic",
                              "pairing_type": "singlet", "solver_mode": solver_mode,
                              "max_iter": 50, "num_eigenvalues": 4}}
        if spectral_shift is not None:
            inp["eliashberg"]["spectral_shift"] = spectral_shift
        sc.calc_eliashberg(inp)
    finally:
        np.random.default_rng = _real_default_rng


def _npz_members(path):
    with np.load(path, allow_pickle=False) as d:
        return {k: np.asarray(d[k]) for k in d.files}


def _assert_text_close(case, got, want, tag, rtol=1e-9):
    """Compare two text outputs token by token: numeric tokens to a relative
    tolerance, every other token exactly. The golden files were recorded on
    one platform; another BLAS / FFT build differs in the last bits, which
    the ``%.8e`` formatting can flip in the last printed digit."""
    punct = "()[]{},;:"

    def _strip(text):
        # the sector-weight header line is new in this version and is not part
        # of the recorded golden files; compare everything else verbatim
        return "\n".join(ln for ln in text.splitlines()
                         if not ln.startswith("# gap_sector_weights"))

    ga, gb = _strip(got).split(), _strip(want).split()
    case.assertEqual(len(ga), len(gb), "{}: token count".format(tag))
    for x, y in zip(ga, gb):
        # a number may be wrapped in punctuation ("(spectral_shift=0.038),"):
        # compare the wrapping exactly and the numeric core to the tolerance
        cx, cy = x.strip(punct), y.strip(punct)
        case.assertEqual(x.replace(cx, "", 1), y.replace(cy, "", 1), tag)
        if "=" in cx or "=" in cy:                  # "spectral_shift=0.5"
            kx, _, cx = cx.partition("=")
            ky, _, cy = cy.partition("=")
            case.assertEqual(kx, ky, tag)
        try:
            fx, fy = float(cx), float(cy)
        except ValueError:
            case.assertEqual(x, y, tag)
            continue
        case.assertTrue(abs(fx - fy) <= rtol * max(1.0, abs(fy)),
                        "{}: {} != {}".format(tag, x, y))


def _assert_member_close(a, b, tag):
    """Bitwise for non-float members, ``allclose`` for float / complex ones
    (round-off of the platform's BLAS / FFT, not of the refactor)."""
    if a.dtype.kind in "fc":
        np.testing.assert_allclose(a, b, rtol=1e-10, atol=1e-13, err_msg=tag)
    else:
        np.testing.assert_array_equal(a, b, err_msg=tag)


class TestDynamicGolden(unittest.TestCase):
    def test_outputs_match_golden(self):
        if not os.path.isdir(_GOLDEN):
            self.skipTest("golden directory missing; record it with RECORD_GOLDEN=1")
        for mode, shift in _CASES:
            tag = "{}_{}".format(mode, "noshift" if shift is None else "shift")
            tmp = tempfile.mkdtemp()
            try:
                _run(tmp, mode, shift)
                for fn in ("eigenvalue.dat", "gap.dat"):
                    with open(os.path.join(tmp, fn)) as fa, \
                            open(os.path.join(_GOLDEN, tag, fn)) as fb:
                        _assert_text_close(self, fa.read(), fb.read(),
                                           "{} {}".format(tag, fn))
                a = _npz_members(os.path.join(tmp, "gap_dynamic.npz"))
                b = _npz_members(os.path.join(_GOLDEN, tag, "gap_dynamic.npz"))
                # the current run may carry keys the golden predates (the
                # sector weights); every golden key must still be there and
                # unchanged
                self.assertTrue(set(b) <= set(a), tag)
                for k in b:
                    _assert_member_close(a[k], b[k], "{} {}".format(tag, k))
                # the new keys, and the singlet gap's purity in the
                # conventional (even k, even w) sector
                self.assertIn("gap_sector_weights", a, tag)
                self.assertIn("gap_sector_labels", a, tag)
                labels = [str(x) for x in a["gap_sector_labels"]]
                self.assertEqual(labels, ["even_k_even_w", "odd_k_even_w",
                                          "even_k_odd_w", "odd_k_odd_w"], tag)
                w = np.asarray(a["gap_sector_weights"], dtype=float)
                self.assertAlmostEqual(float(w.sum()), 1.0, places=9, msg=tag)
                # a 2x2x1 grid has no odd-k function at all, so both odd-k
                # sectors are empty. These fixtures are a synthetic random
                # susceptibility archive whose kernel does not commute with
                # the frequency reversal, so the gap keeps a sizeable
                # odd-frequency admixture -- exactly what the new diagnostic
                # is there to make visible, and the reason the purity of a
                # returned gap is asserted on the physical fixtures instead.
                self.assertEqual(float(w[1]), 0.0, tag)
                self.assertEqual(float(w[3]), 0.0, tag)
                # the recorded weights describe the gap stored next to them
                from hwave.solver import eliashberg_dynamic as ed
                recomputed = ed.gap_sector_weights(a["gap"])
                for i, label in enumerate(labels):
                    self.assertAlmostEqual(float(w[i]), recomputed[label],
                                           places=10, msg=tag)
            finally:
                shutil.rmtree(tmp, ignore_errors=True)


class TestEigenDriverUnits(unittest.TestCase):
    """Direct unit coverage for build_seed, run_leading_eigenproblem's
    parity_leakage_policy branches, and write_eigenvalue_file -- exercised
    only indirectly (and never for a leaky, non-commuting kernel) by the
    golden end-to-end runs above."""

    def setUp(self):
        from hwave.solver import eliashberg_dynamic as ed
        self.ed = ed
        self.gap_shape = (1, 1, 2, 2, 1, 4)
        vec_size = int(np.prod(self.gap_shape))
        # A fixed, non-symmetric complex matrix: generic enough that it does
        # NOT commute with the channel's combined (k, orbital, frequency)
        # parity, so the leakage probe inside run_leading_eigenproblem is
        # exercised on a genuinely leaky kernel (not the physical, roughly-
        # parity-preserving FLEX kernel the golden runs use).
        rng = np.random.default_rng(11)
        M = rng.standard_normal((vec_size, vec_size)) + 1j * rng.standard_normal((vec_size, vec_size))

        def matvec(x):
            return M @ x

        self.matvec = matvec
        # Exercise build_seed with an empty eli_param (defaults: init_gap
        # resolved from pairing_type -> "cos" for singlet).
        kx = np.array([0.0, np.pi])
        ky = np.array([0.0, np.pi])
        kz = np.array([0.0])
        self.phi0, self.seed_vec = ed.build_seed(
            {}, "singlet", 1, kx, ky, kz, self.gap_shape, False, None, 4)

    def test_build_seed_shape_and_normalization(self):
        self.assertEqual(self.phi0.shape, self.gap_shape)
        self.assertTrue(np.iscomplexobj(self.phi0))
        self.assertAlmostEqual(np.linalg.norm(self.phi0), 1.0)
        self.assertIsNone(self.seed_vec)

    def test_refuse_raises_before_any_solve_for_iteration(self):
        eli_param = {"solver_mode": "iteration", "max_iter": 5}
        with self.assertRaisesRegex(
                ValueError, "does not commute with the combined parity"):
            self.ed.run_leading_eigenproblem(
                self.matvec, self.gap_shape, eli_param, "singlet",
                phi0=self.phi0, seed_vec=self.seed_vec, use_ir=False, axF=None,
                nmat=4, parity_leakage_policy="refuse")

    def test_refuse_raises_before_any_solve_for_eigenvalue(self):
        eli_param = {"solver_mode": "eigenvalue", "num_eigenvalues": 2}
        with self.assertRaisesRegex(
                ValueError, "does not commute with the combined parity"):
            self.ed.run_leading_eigenproblem(
                self.matvec, self.gap_shape, eli_param, "singlet",
                phi0=self.phi0, seed_vec=self.seed_vec, use_ir=False, axF=None,
                nmat=4, parity_leakage_policy="refuse")

    def test_warn_completes_and_logs_existing_message(self):
        eli_param = {"solver_mode": "iteration", "max_iter": 3}
        with self.assertLogs("qlms.eliashberg_dynamic", level="WARNING") as cm:
            lam, gap_w, eigenvalues_all, eigenvalue_match, note, leakage, weights = \
                self.ed.run_leading_eigenproblem(
                    self.matvec, self.gap_shape, eli_param, "singlet",
                    phi0=self.phi0, seed_vec=self.seed_vec, use_ir=False,
                    axF=None, nmat=4, parity_leakage_policy="warn")
        self.assertTrue(
            any("does not commute with parity" in msg for msg in cm.output),
            cm.output)
        self.assertIsInstance(lam, float)
        self.assertEqual(gap_w.shape, self.gap_shape)
        self.assertIsNone(eigenvalues_all)
        self.assertIsNone(eigenvalue_match)
        # the driver reports the probe value it measured (the "warn" policy
        # probes on the iteration path)
        self.assertIsInstance(leakage, float)
        self.assertGreater(leakage, 1.0e-8)
        # the driver also reports the sector composition of the gap it returns
        self.assertEqual(set(weights), {"even_k_even_w", "odd_k_even_w",
                                        "even_k_odd_w", "odd_k_odd_w"})
        self.assertAlmostEqual(sum(weights.values()), 1.0, places=9)

    def test_parity_probe_runs_at_most_once_and_only_where_it_is_needed(self):
        """The probe is a full matvec pair, so a duplicate one silently
        doubles the cost of every dynamic solve. It runs exactly once on the
        paths that consume it -- the iteration path (which decides whether to
        project) and the "refuse" policy -- and not at all on the eigenvalue
        family under the default policy."""
        from unittest import mock

        real = self.ed._parity_leakage

        def run(eli_param, **kw):
            calls = []

            def counting(*a, **k):
                calls.append(1)
                return real(*a, **k)

            with mock.patch.object(self.ed, "_parity_leakage", side_effect=counting):
                try:
                    with self.assertLogs("qlms.eliashberg_dynamic", level="WARNING"):
                        self.ed.run_leading_eigenproblem(
                            self.matvec, self.gap_shape, eli_param, "singlet",
                            phi0=self.phi0, seed_vec=self.seed_vec, use_ir=False,
                            axF=None, nmat=4, **kw)
                except ValueError:
                    pass                       # the "refuse" policy, as designed
            return len(calls)

        self.assertEqual(run({"solver_mode": "iteration", "max_iter": 3}), 1)
        self.assertEqual(run({"solver_mode": "eigenvalue", "num_eigenvalues": 2}), 0)
        self.assertEqual(run({"solver_mode": "both", "num_eigenvalues": 2}), 0)
        self.assertEqual(run({"solver_mode": "eigenvalue", "num_eigenvalues": 2},
                             parity_leakage_policy="refuse"), 1)
        self.assertEqual(run({"solver_mode": "iteration", "max_iter": 3},
                             parity_leakage_policy="refuse"), 1)

    def test_invalid_policy_raises(self):
        eli_param = {"solver_mode": "iteration"}
        with self.assertRaisesRegex(ValueError, "parity_leakage_policy"):
            self.ed.run_leading_eigenproblem(
                self.matvec, self.gap_shape, eli_param, "singlet",
                phi0=self.phi0, seed_vec=self.seed_vec, use_ir=False, axF=None,
                nmat=4, parity_leakage_policy="bogus")

    def test_write_eigenvalue_file_header_lines(self):
        tmp = tempfile.mkdtemp()
        try:
            path = os.path.join(tmp, "eigenvalue.dat")
            self.ed.write_eigenvalue_file(
                path, 1.5, None, None, None, header_lines=["x=1"])
            with open(path) as f:
                lines = f.read().splitlines()
            self.assertEqual(lines[0], "# Dynamic Eliashberg leading eigenvalue")
            self.assertEqual(lines[1], "# x=1")
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    if os.environ.get("RECORD_GOLDEN") == "1":
        for mode, shift in _CASES:
            tag = "{}_{}".format(mode, "noshift" if shift is None else "shift")
            dst = os.path.join(_GOLDEN, tag)
            os.makedirs(dst, exist_ok=True)
            tmp = tempfile.mkdtemp()
            _run(tmp, mode, shift)
            for fn in ("eigenvalue.dat", "gap.dat", "gap_dynamic.npz"):
                shutil.copy(os.path.join(tmp, fn), dst)
            shutil.rmtree(tmp, ignore_errors=True)
    unittest.main()
