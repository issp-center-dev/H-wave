"""The momentum-convention provenance gate (issue #133).

k/q-space NPZ files written since the Fourier-sign alignment carry
momentum_convention = "e_plus_ikR". Loaders reject a mismatching tag
outright, and accept unmarked legacy files ONLY when the stored payload
is elementwise even under q -> -q on the FFT grid (then both
conventions coincide bit-for-bit; every centrosymmetric fixture
qualifies). This mirrors the sc_vertex_version content-conditioned
fail-closed precedent.
"""

import os
import tempfile
import unittest

import numpy as np

import hwave.sc as sc
from hwave.solver.rpa import (MOMENTUM_CONVENTION,
                              validate_momentum_convention)


def _q_even(nx):
    """(nfreq, nvol, 1, 1) payload even under q -> -q."""
    q = np.arange(nx)
    vals = np.cos(2 * np.pi * q / nx)
    return np.tile(vals[None, :, None, None], (2, 1, 1, 1))


def _q_odd(nx):
    q = np.arange(nx)
    vals = np.cos(2 * np.pi * q / nx) + 0.3 * np.sin(2 * np.pi * q / nx)
    return np.tile(vals[None, :, None, None], (2, 1, 1, 1))


class TestValidatorUnit(unittest.TestCase):

    def _npz(self, payload, **extra):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = os.path.join(tmp.name, "x.npz")
        np.savez(path, chi0q=payload, **extra)
        return np.load(path), path

    def test_matching_tag_passes(self):
        data, path = self._npz(_q_odd(4),
                               momentum_convention=MOMENTUM_CONVENTION)
        validate_momentum_convention(data, path, data["chi0q"], 1, (4, 1, 1))

    def test_mismatching_tag_is_rejected(self):
        data, path = self._npz(_q_even(4), momentum_convention="e_minus_ikR")
        with self.assertRaises(ValueError) as cm:
            validate_momentum_convention(data, path, data["chi0q"], 1,
                                         (4, 1, 1))
        self.assertIn("momentum_convention", str(cm.exception))

    def test_unmarked_q_even_payload_is_accepted(self):
        data, path = self._npz(_q_even(4))
        validate_momentum_convention(data, path, data["chi0q"], 1, (4, 1, 1))

    def test_unmarked_q_odd_payload_is_rejected(self):
        data, path = self._npz(_q_odd(4))
        with self.assertRaises(ValueError) as cm:
            validate_momentum_convention(data, path, data["chi0q"], 1,
                                         (4, 1, 1))
        self.assertIn("#133", str(cm.exception))


class TestLoaderGates(unittest.TestCase):
    """Each production loader must reach the gate."""

    def test_flex_green_loader_rejects_unmarked_q_odd_file(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        nx = 4
        q = np.arange(nx)
        prof = 1.0 + 0.3 * np.sin(2 * np.pi * q / nx)
        green = np.ones((1, 4, nx, 1, 1), dtype=complex) * \
            prof[None, None, :, None, None]
        np.savez(os.path.join(d, "green.npz"), green=green, beta=0.5)
        inp = {"mode": {"param": {"T": 2.0}},
               "file": {"output": {"path_to_output": d}},
               "eliashberg": {}}
        with self.assertRaises(ValueError) as cm:
            sc._load_flex_green(inp, 1, nx, 1, 1)
        self.assertIn("#133", str(cm.exception))

    def test_flex_green_loader_accepts_marked_file(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        nx = 4
        q = np.arange(nx)
        prof = 1.0 + 0.3 * np.sin(2 * np.pi * q / nx)
        green = np.ones((1, 4, nx, 1, 1), dtype=complex) * \
            prof[None, None, :, None, None]
        np.savez(os.path.join(d, "green.npz"), green=green, beta=0.5,
                 momentum_convention=MOMENTUM_CONVENTION)
        inp = {"mode": {"param": {"T": 2.0}},
               "file": {"output": {"path_to_output": d}},
               "eliashberg": {}}
        green_loaded = sc._load_flex_green(inp, 1, nx, 1, 1)
        self.assertIsNotNone(green_loaded)

    def test_sc_chi0q_loader_rejects_unmarked_q_odd_file(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        nx = 4
        chi0q = _q_odd(nx)
        np.savez(os.path.join(d, "chi0q.npz"), chi0q=chi0q)
        inp = {"mode": {"param": {"T": 1.0, "Nmat": 2,
                                  "CellShape": [nx, 1, 1]}},
               "file": {"input": {"path_to_flex_output": d},
                        "output": {"path_to_output": d}},
               "eliashberg": {}}
        with self.assertRaises(ValueError) as cm:
            sc._load_chi0q(inp)
        self.assertIn("#133", str(cm.exception))


class TestExternFourierSign(unittest.TestCase):
    """Extern's R -> k build shares the documented sign."""

    def test_extern_q_is_documented_sign(self):
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.rpa as solver_rpa
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        nx = 4
        t = 0.7 + 0.3j
        h = 0.2 + 0.4j
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n1\n0.0 0.0 0.0\n")
        with open(os.path.join(d, "transfer.dat"), "w") as f:
            f.write("hdr\n1\n2\n1 1\n")
            f.write(" 1 0 0 1 1 %.12f %.12f\n" % (t.real, t.imag))
            f.write("-1 0 0 1 1 %.12f %.12f\n"
                    % (np.conj(t).real, np.conj(t).imag))
        with open(os.path.join(d, "extern.dat"), "w") as f:
            f.write("hdr\n1\n2\n1 1\n")
            f.write(" 1 0 0 1 1 %.12f %.12f\n" % (h.real, h.imag))
            f.write("-1 0 0 1 1 %.12f %.12f\n"
                    % (np.conj(h).real, np.conj(h).imag))
        io = read_input_k.QLMSkInput(
            {"path_to_input": d,
             "interaction": {"path_to_input": d, "Geometry": "geom.dat",
                             "Transfer": "transfer.dat",
                             "Extern": "extern.dat"}})
        info_mode = {"mode": "RPA",
                     "param": {"T": 1.0, "mu": 0.0, "CellShape": [nx, 1, 1],
                               "SubShape": [1, 1, 1], "Nmat": 8},
                     "calc_scheme": "reduced"}
        solver = solver_rpa.RPA(io.get_param("ham"), {}, info_mode)
        got = np.asarray(solver.ham_info.ham_extern_q).reshape(nx)
        want = np.array([h * np.exp(2j * np.pi * k / nx)
                         + np.conj(h) * np.exp(-2j * np.pi * k / nx)
                         for k in range(nx)])
        np.testing.assert_allclose(got, want, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
