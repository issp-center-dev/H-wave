"""Input-shape guards of the longitudinal chi0q kernel (issue #125).

RPA._calc_chi0q validated its Green's function's volume and frequency
axes with bare asserts, which vanish under ``python -O``. A wrong
frequency count then produces plausible-looking wrong output; a wrong
volume fails later at the lattice reshape with an unhelpful diagnostic.
This is the longitudinal analogue of the transverse block-count fix;
the guards are now ValueErrors naming the expected and observed values,
pinned here under the normal interpreter and in ``-O`` subprocesses.
The round-3 review added the paired-tail invariant: a green0_tail whose
shape differs from the Green's function -- even at the SAME element
count, where the kernel's reshape succeeds -- must be rejected, in both
the longitudinal and transverse kernels.
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

    def test_mismatched_tail_rejected(self):
        """Same element count, different shape: the reshape succeeds and
        the kernel returns finite, plausible-looking garbage (round-3
        review measured it) -- the paired-tail guard must fire first."""
        stub = self._stub()
        g = np.zeros((1, 4, 4, 2, 2), dtype=complex)
        tail = np.zeros((1, 8, 2, 2, 2), dtype=complex)  # same 64 elements
        with self.assertRaises(ValueError) as cm:
            stub._calc_chi0q(g, tail, 1.0)
        msg = str(cm.exception)
        self.assertIn("green0_tail", msg)
        self.assertIn("(1, 8, 2, 2, 2)", msg)
        self.assertIn("(1, 4, 4, 2, 2)", msg)

    def test_transverse_mismatched_tail_rejected(self):
        """The transverse kernel shares the paired-tail invariant."""
        stub = self._stub()
        g = np.zeros((2, 4, 4, 2, 2), dtype=complex)
        tail = np.zeros((2, 8, 2, 2, 2), dtype=complex)
        with self.assertRaises(ValueError) as cm:
            stub._calc_chi0q_transverse(g, tail, 1.0)
        self.assertIn("green0_tail", str(cm.exception))

    def test_malformed_block_and_orbital_axes_rejected(self):
        """Round 5 reproduced fail-open structural inputs: nblock=0 gave
        a finite EMPTY result, nblock=3 a finite three-block one (which
        the caller's block handling then truncated under -O), nd=0 was
        accepted, and rectangular orbital axes failed only through an
        incidental reshape. All are loud ValueErrors now."""
        stub = self._stub()
        cases = [
            ((0, 4, 4, 2, 2), "block axis (0)"),
            ((3, 4, 4, 2, 2), "block axis (3)"),
            ((1, 4, 4, 0, 0), "square and nonempty"),
            ((1, 4, 4, 2, 3), "square and nonempty"),
        ]
        for shape, frag in cases:
            with self.subTest(shape=shape):
                g = np.zeros(shape, dtype=complex)
                with self.assertRaises(ValueError) as cm:
                    stub._calc_chi0q(g, np.zeros_like(g), 1.0)
                self.assertIn(frag, str(cm.exception))

    def test_transverse_axis_checks(self):
        """The transverse kernel now validates volume, frequency and the
        orbital square too (a wrong transverse frequency count was
        processed normally before round 5)."""
        stub = self._stub()
        cases = [
            ((2, 4, 5, 2, 2), "volume axis (5)"),
            ((2, 6, 4, 2, 2), "frequency axis (6)"),
            ((2, 4, 4, 2, 3), "square and nonempty"),
        ]
        for shape, frag in cases:
            with self.subTest(shape=shape):
                g = np.zeros(shape, dtype=complex)
                with self.assertRaises(ValueError) as cm:
                    stub._calc_chi0q_transverse(g, np.zeros_like(g), 1.0)
                self.assertIn(frag, str(cm.exception))

    def test_guard_precedence_on_combined_malformations(self):
        """When several things are wrong at once, the axis checks fire
        before the tail check (and the transverse block count before its
        tail check), so the FIRST reported error names the most
        fundamental problem. Behavior was correct; pin it (round 4)."""
        stub = self._stub()
        bad_tail = np.zeros((1, 8, 2, 2, 2), dtype=complex)
        g_bad_vol = np.zeros((1, 4, 5, 2, 2), dtype=complex)
        with self.assertRaises(ValueError) as cm:
            stub._calc_chi0q(g_bad_vol, bad_tail, 1.0)
        self.assertIn("volume axis", str(cm.exception))
        g_bad_freq = np.zeros((1, 6, 4, 2, 2), dtype=complex)
        with self.assertRaises(ValueError) as cm:
            stub._calc_chi0q(g_bad_freq, bad_tail, 1.0)
        self.assertIn("frequency axis", str(cm.exception))
        g3 = np.zeros((3, 4, 4, 2, 2), dtype=complex)
        with self.assertRaises(ValueError) as cm:
            stub._calc_chi0q_transverse(
                g3, np.zeros((3, 8, 2, 2, 2), dtype=complex), 1.0)
        self.assertIn("nblock=3", str(cm.exception))

    def test_volume_guard_survives_python_O(self):
        """The gating CI never runs optimized Python; without this pin a
        regression back to a bare assert would pass CI while vanishing
        in production -O runs."""
        self._run_O_subprocess((1, 4, 5, 2, 2), "volume axis (5)")

    def test_frequency_guard_survives_python_O(self):
        self._run_O_subprocess((1, 6, 4, 2, 2), "frequency axis (6)")

if __name__ == "__main__":
    unittest.main()
