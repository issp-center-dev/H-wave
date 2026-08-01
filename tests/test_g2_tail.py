"""The pair-bubble Matsubara tail correction (#86).

sc._calc_g2's truncated sum (1/beta) sum_n G(k,iw_n) G(-k,-iw_n) misses
the leading identity tail of the exact pair bubble (an O(1/Nmat) error),
which can leave G2 slightly indefinite and inject spurious imaginary parts
into the static Eliashberg eigenvalues.
The correction adds the analytic tail of the exact summand model
delta_ij delta_lm / wn^2, i.e. c * identity on the (i,l) gap space with
c = beta/4 - (1/beta) sum_{n in window} 1/wn^2.

These tests pin the correction against an INDEPENDENT band-basis evaluation
of the exact Matsubara sum (the analytic pair bubble
(1 - f(e_a) - f(e_b)) / (e_a + e_b)), the positivity restoration, the exact
shape of the correction term, the config plumbing through the public entry
(including a differently-cased key, per the case-insensitivity defect
class), and the indefiniteness warning.
"""

import logging
import os
import tempfile
import unittest

import numpy as np

import hwave.sc as sc
from hwave.sc import _calc_g2, _warn_if_g2_indefinite


def _model(norb=2, Nx=4, Ny=4, Nz=1, seed=86):
    """Random Hermitian real-space hoppings -> H(k) on the FFT grid."""
    rng = np.random.RandomState(seed)
    hops = {}
    for R in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)]:
        m = rng.randn(norb, norb) + 1j * rng.randn(norb, norb)
        if R == (0, 0, 0):
            m = 0.5 * (m + m.conj().T)
        hops[R] = m
    H = np.zeros((Nx, Ny, Nz, norb, norb), dtype=complex)
    for kx in range(Nx):
        for ky in range(Ny):
            for kz in range(Nz):
                h = np.zeros((norb, norb), dtype=complex)
                for R, m in hops.items():
                    ph = np.exp(2j * np.pi * (kx * R[0] / Nx
                                              + ky * R[1] / Ny
                                              + kz * R[2] / Nz))
                    h += m * ph + m.conj().T * np.conj(ph)
                H[kx, ky, kz] = h
    return H


def _green(H, beta, nmat):
    Nx, Ny, Nz, norb, _ = H.shape
    iw = 1j * (2 * np.arange(nmat) + 1 - nmat) * np.pi / beta
    g = np.zeros((norb, norb, Nx, Ny, Nz, nmat), dtype=complex)
    eye = np.eye(norb)
    for kx in range(Nx):
        for ky in range(Ny):
            for kz in range(Nz):
                for i, w in enumerate(iw):
                    g[:, :, kx, ky, kz, i] = np.linalg.inv(
                        w * eye - H[kx, ky, kz])
    return g


def _exact_pair_bubble(H, beta):
    """Band-basis analytic Matsubara sum -- independent of _calc_g2."""
    Nx, Ny, Nz, norb, _ = H.shape

    def fermi(e):
        return 1.0 / (np.exp(np.clip(beta * e, -500, 500)) + 1.0)

    W = np.zeros((norb, norb, norb, norb, Nx, Ny, Nz), dtype=complex)
    for kx in range(Nx):
        for ky in range(Ny):
            for kz in range(Nz):
                ew, U = np.linalg.eigh(H[kx, ky, kz])
                mk = ((Nx - kx) % Nx, (Ny - ky) % Ny, (Nz - kz) % Nz)
                ew2, V = np.linalg.eigh(H[mk])
                for a in range(norb):
                    for b in range(norb):
                        s = ((1.0 - fermi(ew[a]) - fermi(ew2[b]))
                             / (ew[a] + ew2[b]))
                        W[:, :, :, :, kx, ky, kz] += s * np.einsum(
                            "i,j,l,m->ijlm", U[:, a], U[:, a].conj(),
                            V[:, b], V[:, b].conj())
    return W


def _min_eigenvalue(G2):
    """Min eigenvalue of G2 as the (il),(jm) matrix on the gap space."""
    norb = G2.shape[0]
    nvol = int(np.prod(G2.shape[4:]))
    M = G2.reshape(norb, norb, norb, norb, nvol)
    M = M.transpose(4, 0, 2, 1, 3).reshape(nvol, norb * norb, norb * norb)
    M = 0.5 * (M + np.conj(M.transpose(0, 2, 1)))
    return np.linalg.eigvalsh(M).min()


class TestG2TailCorrection(unittest.TestCase):
    beta = 10.0
    nmat = 64

    def setUp(self):
        self.H = _model()
        self.green = _green(self.H, self.beta, self.nmat)
        self.exact = _exact_pair_bubble(self.H, self.beta)

    def test_correction_approaches_the_exact_band_basis_sum(self):
        bare = _calc_g2(self.green, self.beta, tail=False)
        corr = _calc_g2(self.green, self.beta)  # default tail=True
        err_bare = np.abs(bare - self.exact).max()
        err_corr = np.abs(corr - self.exact).max()
        # At Nmat=64, beta=10 the bare truncation error is ~1.6e-2 and the
        # corrected one ~1.6e-3 (next order); require a solid margin rather
        # than the measured ratio so the pin is not brittle.
        self.assertGreater(err_bare, 1.0e-2)
        self.assertLess(err_corr, err_bare / 5.0)

    def test_correction_restores_positive_semidefiniteness(self):
        bare = _calc_g2(self.green, self.beta, tail=False)
        corr = _calc_g2(self.green, self.beta)
        self.assertLess(_min_eigenvalue(bare), -1.0e-3)
        self.assertGreater(_min_eigenvalue(corr), -1.0e-10)

    def test_correction_is_the_analytic_identity_shift(self):
        """Locks the implemented term to c * delta_ij delta_lm, uniform in k,
        with c the exact window-complement of (1/beta) sum 1/wn^2."""
        bare = _calc_g2(self.green, self.beta, tail=False)
        corr = _calc_g2(self.green, self.beta)
        wn = (2.0 * np.arange(self.nmat) + 1.0 - self.nmat) * np.pi / self.beta
        c = self.beta / 4.0 - np.sum(1.0 / wn**2) / self.beta
        self.assertGreater(c, 0.0)
        norb = bare.shape[0]
        diff = corr - bare
        mask = np.zeros(bare.shape, dtype=bool)
        for i in range(norb):
            for l in range(norb):
                mask[i, i, l, l] = True
        # Off the (i,i,l,l) positions the two calls share every float op, so
        # the difference is exactly zero; on them (base + c) - base rounds,
        # so compare to c with a tight tolerance instead of exact equality.
        np.testing.assert_array_equal(diff[~mask],
                                      np.zeros_like(diff[~mask]))
        np.testing.assert_allclose(diff[mask], c, rtol=0, atol=1.0e-12)

    def test_warning_fires_on_indefinite_g2_and_not_after_correction(self):
        bare = _calc_g2(self.green, self.beta, tail=False)
        corr = _calc_g2(self.green, self.beta)
        with self.assertLogs("hwave_sc", level="WARNING") as cm:
            _warn_if_g2_indefinite(bare, bare.shape[0], tail_enabled=False)
        self.assertTrue(any("positive semi-definite" in m
                            for m in cm.output))
        self.assertTrue(any("g2_tail" in m for m in cm.output))
        # No-warning branch, Python 3.9-compatible (no assertNoLogs): collect
        # records with a handler and require zero warnings.
        records = []

        class _Grab(logging.Handler):
            def emit(self, record):
                records.append(record)

        grab = _Grab(level=logging.WARNING)
        target = logging.getLogger("hwave_sc")
        target.addHandler(grab)
        try:
            _warn_if_g2_indefinite(corr, corr.shape[0], tail_enabled=True)
        finally:
            target.removeHandler(grab)
        self.assertEqual([r.getMessage() for r in records], [])


class TestG2TailGuards(unittest.TestCase):
    """Grid validation, off-switch fidelity, and the diagnostic's branches."""

    beta = 10.0

    def test_odd_or_nonpositive_nmat_is_rejected_with_tail(self):
        H = _model(Nx=2, Ny=2)
        for nmat in (1, 3, 7):
            green = _green(H, self.beta, nmat)
            with self.assertRaises(ValueError) as cm:
                with np.errstate(divide="ignore", invalid="ignore"):
                    _calc_g2(green, self.beta)
            self.assertIn("g2_tail", str(cm.exception))
        empty = np.zeros((2, 2, 2, 2, 1, 0), dtype=complex)
        with self.assertRaises(ValueError):
            _calc_g2(empty, self.beta)

    def test_odd_nmat_guard_survives_python_O(self):
        """The guard is an explicit raise, not an assert: it must still fire
        under python -O (bare-assert defect class)."""
        import subprocess
        import sys
        code = (
            "import numpy as np\n"
            "from hwave.sc import _calc_g2\n"
            "g = np.ones((1, 1, 2, 1, 1, 3), dtype=complex)\n"
            "try:\n"
            "    _calc_g2(g, 10.0)\n"
            "except ValueError:\n"
            "    print('RAISED')\n")
        env = dict(os.environ)
        src = os.path.join(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))), "src")
        if os.path.isdir(src):
            env["PYTHONPATH"] = src + os.pathsep + env.get("PYTHONPATH", "")
        out = subprocess.run([sys.executable, "-O", "-c", code],
                             capture_output=True, text=True, env=env)
        self.assertIn("RAISED", out.stdout)

    def test_tail_off_is_bit_identical_to_the_pre_change_sum(self):
        """The off-switch must reproduce the historical implementation
        exactly, pinned here as a frozen copy of the pre-change matmul."""
        from hwave.solver.kgrid import reverse_fft_axes
        H = _model()
        green = _green(H, self.beta, 32)
        norb = green.shape[0]
        Nx, Ny, Nz, nmat = green.shape[2:]
        nvol = Nx * Ny * Nz
        inv = reverse_fft_axes(green[..., ::-1], (2, 3, 4))
        A = green.reshape(norb * norb, nvol, nmat)
        B = inv.reshape(norb * norb, nvol, nmat)
        As = np.moveaxis(A, 1, 0)
        Bs = np.moveaxis(B, 1, 0)
        ref = np.moveaxis(As @ Bs.transpose(0, 2, 1), 0, 2)
        ref = ref.reshape(norb, norb, norb, norb, Nx, Ny, Nz) / self.beta
        self.assertTrue(np.array_equal(
            _calc_g2(green, self.beta, tail=False), ref))

    def test_non_hermitian_g2_gets_its_own_warning(self):
        """Hermitization must not silently hide a malformed loaded Green
        function: a significant anti-Hermitian part is reported first."""
        norb = 2
        g2 = np.zeros((norb, norb, norb, norb, 1, 1, 1), dtype=complex)
        for i in range(norb):
            for l in range(norb):
                g2[i, i, l, l] = 1.0
        g2[0, 1, 0, 0] = 1.0j  # M[(0,0),(1,0)] = i with no conjugate partner
        with self.assertLogs("hwave_sc", level="WARNING") as cm:
            _warn_if_g2_indefinite(g2, norb, tail_enabled=True)
        self.assertTrue(any("non-Hermitian" in m for m in cm.output))

    def test_size_guard_skips_the_check(self):
        """Both guard dimensions independently: the work bound and the byte
        bound must each trigger the skip on their own."""
        H = _model(Nx=2, Ny=2)
        g2 = _calc_g2(_green(H, self.beta, 32), self.beta)
        for attr in ("_G2_CHECK_MAX_WORK", "_G2_CHECK_MAX_BYTES"):
            with self.subTest(threshold=attr):
                saved = getattr(sc, attr)
                setattr(sc, attr, 1)
                try:
                    with self.assertLogs("hwave_sc", level="INFO") as cm:
                        result = _warn_if_g2_indefinite(
                            g2, 2, tail_enabled=True)
                finally:
                    setattr(sc, attr, saved)
                self.assertIsNone(result)
                self.assertTrue(any("skipped" in m for m in cm.output))

    def test_flex_loaded_green_with_odd_grid_is_rejected(self):
        """The previously-unvalidated production route: _load_flex_green
        accepts a file with an odd frequency axis (it validates shape and
        finiteness, not grid parity), and the tail correction is then the
        first place the invalid grid can corrupt results -- pin that the
        production loader output hits the _calc_g2 guard."""
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        norb, Nx, Ny, Nz, nmat = 1, 2, 1, 1, 5
        green_file = np.ones((1, nmat, Nx * Ny * Nz, norb, norb),
                             dtype=complex)
        np.savez(os.path.join(d, "green.npz"), green=green_file)
        inp = {"file": {"output": {"path_to_output": d}},
               "eliashberg": {}}
        green = sc._load_flex_green(inp, norb, Nx, Ny, Nz)
        self.assertIsNotNone(green)
        self.assertEqual(green.shape[-1], nmat)
        with self.assertRaises(ValueError) as cm:
            _calc_g2(green, self.beta)
        self.assertIn("g2_tail", str(cm.exception))


class TestG2TailDressedGreen(unittest.TestCase):
    """The identity 1/(i wn) coefficient claim, exercised on a DRESSED Green
    function: G = (i wn - H - Sigma(i wn))^{-1} with a one-pole self-energy.
    The self-energy only enters at O(1/wn^2), so the same correction must
    accelerate convergence toward the (numerically converged) exact sum."""

    def test_correction_accelerates_convergence_for_dyson_green(self):
        beta = 10.0
        H = _model(Nx=2, Ny=2, seed=87)
        Nx, Ny, Nz, norb, _ = H.shape
        v2, eps = 0.5, 0.3

        def dressed(nmat):
            iw = 1j * (2 * np.arange(nmat) + 1 - nmat) * np.pi / beta
            g = np.zeros((norb, norb, Nx, Ny, Nz, nmat), dtype=complex)
            eye = np.eye(norb)
            for kx in range(Nx):
                for ky in range(Ny):
                    for kz in range(Nz):
                        for i, w in enumerate(iw):
                            sig = v2 / (w - eps) * eye
                            g[:, :, kx, ky, kz, i] = np.linalg.inv(
                                w * eye - H[kx, ky, kz] - sig)
            return g

        # Converged reference: tail=False at a large Nmat plus its own
        # analytic tail (the same identity-shift formula), so the reference
        # does not share the finite window under test.
        nref = 8192
        ref = _calc_g2(dressed(nref), beta)
        g64 = dressed(64)
        err_bare = np.abs(_calc_g2(g64, beta, tail=False) - ref).max()
        err_corr = np.abs(_calc_g2(g64, beta) - ref).max()
        self.assertGreater(err_bare, 1.0e-2)
        self.assertLess(err_corr, err_bare / 5.0)


class TestG2TailConfigPlumbing(unittest.TestCase):
    """g2_tail must act through the public entry, case-insensitively."""

    def _run(self, eliashberg_extra):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("  1.0 0.0 0.0\n  0.0 1.0 0.0\n  0.0 0.0 1.0\n"
                    "1\n 0.0 0.0 0.0\n")
        with open(os.path.join(d, "transfer.dat"), "w") as f:
            f.write("Transfer\n1\n2\n 1 1\n"
                    "   1    0    0    1    1   1.000000   0.0\n"
                    "  -1    0    0    1    1   1.000000   0.0\n")
        with open(os.path.join(d, "coulombintra.dat"), "w") as f:
            f.write("CoulombIntra\n1\n1\n 1\n"
                    "   0    0    0    1    1   2.000000   0.0\n")
        out = os.path.join(d, "out")
        os.makedirs(out, exist_ok=True)
        eli = {"chi0q_mode": "calc", "frequency": "static",
               "pairing_type": "singlet", "init_gap": "cos",
               "solver_mode": "eigenvalue",
               "eigenvalue_method": "arnoldi", "num_eigenvalues": 1,
               "output_eigenvalue": "eig.dat", "output_gap": "gap.dat"}
        eli.update(eliashberg_extra)
        nmat = getattr(self, "_nmat_override", None) or 8
        inp = {
            "mode": {"param": {"T": 2.0, "filling": 0.5,
                               "CellShape": [4, 4, 1],
                               "SubShape": [1, 1, 1], "Nmat": nmat}},
            "file": {"input": {"path_to_input": d,
                               "interaction": {"path_to_input": d,
                                               "Geometry": "geom.dat",
                                               "Transfer": "transfer.dat",
                                               "CoulombIntra":
                                                   "coulombintra.dat"}},
                     "output": {"path_to_output": out}},
            "eliashberg": eli,
        }
        sc.calc_eliashberg(inp)
        with open(os.path.join(out, "eig.dat")) as f:
            lines = [ln for ln in f if ln.strip()
                     and not ln.lstrip().startswith("#")]
        return float(lines[0].split()[0])

    def test_default_on_off_switch_and_case_insensitive_key(self):
        lam_default = self._run({})
        lam_off = self._run({"g2_tail": False})
        lam_off_cased = self._run({"G2_Tail": False})
        lam_off_string = self._run({"g2_tail": "off"})
        lam_on_explicit = self._run({"g2_tail": True})
        # The switch must change the result at Nmat=8, the default must be
        # ON, the differently-cased key must behave identically to the
        # lowercase one (case-insensitivity defect class, PR #128), and a
        # string boolean ("off") must not be read as truthy.
        self.assertNotAlmostEqual(lam_default, lam_off, places=6)
        self.assertEqual(lam_default, lam_on_explicit)
        self.assertEqual(lam_off, lam_off_cased)
        self.assertEqual(lam_off, lam_off_string)

    def test_odd_nmat_is_rejected_through_the_public_entry(self):
        """An odd grid must fail loudly somewhere on the public path, never
        reach a silent inf-laden G2. On the chi0q_mode='calc' route the RPA
        layer's own even-Nmat check fires first (SystemExit); the _calc_g2
        ValueError guard covers the load/flex routes, where the Green
        function's frequency axis arrives from a file with no earlier
        validation (pinned directly in TestG2TailGuards)."""
        with np.errstate(all="ignore"):
            with self.assertRaises((ValueError, SystemExit)):
                self._run_with_nmat(7)

    def _run_with_nmat(self, nmat):
        # thin wrapper so the odd-Nmat rejection exercises the same public
        # config path as the plumbing test
        self._nmat_override = nmat
        try:
            return self._run({})
        finally:
            self._nmat_override = None


if __name__ == "__main__":
    unittest.main()
