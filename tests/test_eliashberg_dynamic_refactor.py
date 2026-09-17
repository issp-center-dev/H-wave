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
                        self.assertEqual(fa.read(), fb.read(), "{} {}".format(tag, fn))
                a = _npz_members(os.path.join(tmp, "gap_dynamic.npz"))
                b = _npz_members(os.path.join(_GOLDEN, tag, "gap_dynamic.npz"))
                self.assertEqual(set(a), set(b), tag)
                for k in a:
                    np.testing.assert_array_equal(a[k], b[k], err_msg="{} {}".format(tag, k))
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
