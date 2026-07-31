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
        # role-BOUND fragments (round-1 review: value-anywhere checks let
        # a message swapping the roles, or an unrelated ValueError
        # containing the digits, pass)
        msg = str(cm.exception)
        self.assertIn("volume axis (5)", msg)
        self.assertIn("lattice (4)", msg)

    def test_wrong_frequency_count_rejected(self):
        stub = self._stub()
        g = np.zeros((1, 6, 4, 2, 2), dtype=complex)  # nmat 6 != 4
        with self.assertRaises(ValueError) as cm:
            stub._calc_chi0q(g, np.zeros_like(g), 1.0)
        msg = str(cm.exception)
        self.assertIn("frequency axis (6)", msg)
        self.assertIn("Nmat (4)", msg)

    def _run_O_subprocess(self, bad_shape, want_fragment):
        """One guard per subprocess (round-1 review: a single combined
        script let a vanished guard hide behind the other's failure),
        with explicit if/raise message validation -- an assert INSIDE the
        -O subprocess would itself be optimized away, which made the
        original version of this test accept an unrelated reshape error."""
        code = (
            "import types, numpy as np\n"
            "import hwave.solver.rpa as rpa_module\n"
            "stub = object.__new__(rpa_module.RPA)\n"
            "stub.lattice = types.SimpleNamespace(shape=(2, 2, 1), nvol=4)\n"
            "stub.nmat = 4\n"
            "g = np.zeros({shape}, dtype=complex)\n"
            "try:\n"
            "    stub._calc_chi0q(g, np.zeros_like(g), 1.0)\n"
            "except ValueError as e:\n"
            "    if {frag!r} not in str(e):\n"
            "        raise SystemExit('wrong ValueError: ' + str(e))\n"
            "else:\n"
            "    raise SystemExit('guard vanished under -O')\n"
        ).format(shape=bad_shape, frag=want_fragment)
        env = dict(os.environ)
        src = os.path.join(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))), "src")
        if os.path.isdir(src):
            env["PYTHONPATH"] = src + os.pathsep + env.get("PYTHONPATH", "")
        proc = subprocess.run([sys.executable, "-O", "-c", code],
                              capture_output=True, text=True, env=env)
        self.assertEqual(
            proc.returncode, 0,
            "guard must survive python -O: {}{}".format(
                proc.stdout, proc.stderr))

    def test_volume_guard_survives_python_O(self):
        """The gating CI never runs optimized Python; without this pin a
        regression back to a bare assert would pass CI while vanishing
        in production -O runs."""
        self._run_O_subprocess((1, 4, 5, 2, 2), "volume axis (5)")

    def test_frequency_guard_survives_python_O(self):
        self._run_O_subprocess((1, 6, 4, 2, 2), "frequency axis (6)")

if __name__ == "__main__":
    unittest.main()
