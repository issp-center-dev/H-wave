"""Input-shape guards of the longitudinal chi0q kernel (issue #125).

RPA._calc_chi0q validated its Green's function's volume and frequency
axes with bare asserts, which vanish under ``python -O`` -- a malformed
array would then proceed into the kernel and produce plausible-looking
wrong output. This is the longitudinal analogue of the transverse
block-count fix; both guards are now ValueErrors naming the expected
and observed values, pinned here under the normal interpreter and in a
``-O`` subprocess.
"""

import os
import subprocess
import sys
import types
import unittest

import numpy as np


class TestChi0qShapeGuards(unittest.TestCase):

    @staticmethod
    def _stub(nx=2, ny=2, nz=1, nmat=4):
        import hwave.solver.rpa as rpa_module

        stub = object.__new__(rpa_module.RPA)
        stub.lattice = types.SimpleNamespace(shape=(nx, ny, nz),
                                             nvol=nx * ny * nz)
        stub.nmat = nmat
        return stub

    def test_wrong_volume_rejected(self):
        stub = self._stub()
        g = np.zeros((1, 4, 5, 2, 2), dtype=complex)  # nvol 5 != 4
        with self.assertRaises(ValueError) as cm:
            stub._calc_chi0q(g, np.zeros_like(g), 1.0)
        msg = str(cm.exception)
        self.assertIn("volume", msg)
        self.assertIn("5", msg)
        self.assertIn("4", msg)

    def test_wrong_frequency_count_rejected(self):
        stub = self._stub()
        g = np.zeros((1, 6, 4, 2, 2), dtype=complex)  # nmat 6 != 4
        with self.assertRaises(ValueError) as cm:
            stub._calc_chi0q(g, np.zeros_like(g), 1.0)
        msg = str(cm.exception)
        self.assertIn("frequency", msg)
        self.assertIn("6", msg)

    def test_guards_survive_python_O(self):
        """The gating CI never runs optimized Python; without this pin a
        regression back to bare asserts would pass CI while vanishing in
        production -O runs (same rationale as the transverse guard)."""
        code = (
            "import types, numpy as np\n"
            "import hwave.solver.rpa as rpa_module\n"
            "stub = object.__new__(rpa_module.RPA)\n"
            "stub.lattice = types.SimpleNamespace(shape=(2, 2, 1), nvol=4)\n"
            "stub.nmat = 4\n"
            "g = np.zeros((1, 4, 5, 2, 2), dtype=complex)\n"
            "try:\n"
            "    stub._calc_chi0q(g, np.zeros_like(g), 1.0)\n"
            "except ValueError as e:\n"
            "    assert 'volume' in str(e), str(e)\n"
            "else:\n"
            "    raise SystemExit('volume guard vanished under -O')\n"
            "g = np.zeros((1, 6, 4, 2, 2), dtype=complex)\n"
            "try:\n"
            "    stub._calc_chi0q(g, np.zeros_like(g), 1.0)\n"
            "except ValueError as e:\n"
            "    assert 'frequency' in str(e), str(e)\n"
            "else:\n"
            "    raise SystemExit('frequency guard vanished under -O')\n"
        )
        env = dict(os.environ)
        src = os.path.join(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))), "src")
        if os.path.isdir(src):
            env["PYTHONPATH"] = src + os.pathsep + env.get("PYTHONPATH", "")
        proc = subprocess.run([sys.executable, "-O", "-c", code],
                              capture_output=True, text=True, env=env)
        self.assertEqual(
            proc.returncode, 0,
            "guards must survive python -O: {}{}".format(
                proc.stdout, proc.stderr))


if __name__ == "__main__":
    unittest.main()
