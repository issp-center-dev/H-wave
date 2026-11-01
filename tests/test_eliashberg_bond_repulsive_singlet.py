"""Singlet pairing of the repulsive one-band Hubbard model through the
in-process FLEX bond pairing step, end to end.

The model: square lattice, nearest-neighbour t = 1, U = 3, n = 0.85, T = 0.1
on a 6x6 grid with Nmat = 64, and a nearest-neighbour CoulombInter declared
with value ZERO -- so the bond path builds its displacement channels while
the Hamiltonian is the plain Hubbard model, and the bond and the on-site
(non-bond) kernels must agree.

The dominant (largest-magnitude) eigenvalue of this singlet kernel is a
NEGATIVE on-site s-wave mode (about -2.21 here, -1.93 on 8x8); the physical
leading eigenvalue is the largest real part (+0.506 here). A power iteration
without a spectral shift falls into a sign-alternating two-cycle there; the
in-process step used to report its iterate norm ``|lambda_-|/3`` as "the
leading eigenvalue" (0.643 on 8x8).

Pinned here:

(a) with the defaults (in-process ``solver_mode = "eigenvalue"``) the
    in-process bond pairing reproduces ``hwave_sc``'s non-bond dynamic solver
    on the same FLEX state: lambda, the gap and the top eigenvalues;
(b) the converged gap satisfies the high-frequency identity of the
    instantaneous on-site U, ``lambda Delta(R=0, iw_n) = -U F_loc``;
(c) a constant on-site gap through the bond kernel gives the sign of -U F_loc
    (F_loc > 0) at high frequency, for U > 0 and U < 0;
(d) ``solver_mode = "iteration"`` without a shift reports NaN with
    ``eigen_converged = false`` in every output (the FLEX run itself
    completes), and with the default shift ("auto") converges to (a)'s
    lambda;
and the ``hwave_sc`` command's exit status on a non-converged run.
"""
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

L = 6
NMAT = 64
T = 0.1
FILLING = 0.425


def _write_model(d, U):
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n1\n0.0 0.0 0.0\n")
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("Transfer\n1\n5\n1 1 1 1 1\n")
        for r in ((-1, 0, 0, -1.0), (0, -1, 0, -1.0), (0, 0, 0, 0.0),
                  (0, 1, 0, -1.0), (1, 0, 0, -1.0)):
            f.write("%4d %4d %4d    1    1  %18.12f  0.0\n" % r)
    with open(os.path.join(d, "coulombintra.dat"), "w") as f:
        f.write("CoulombIntra\n1\n1\n1\n   0    0    0    1    1  %18.12f 0.0\n" % U)
    with open(os.path.join(d, "coulombinter.dat"), "w") as f:
        f.write("CoulombInter\n1\n4\n1 1 1 1\n")
        for r in ((-1, 0, 0), (0, -1, 0), (0, 1, 0), (1, 0, 0)):
            f.write("%4d %4d %4d    1    1  0.0  0.0\n" % r)


def _interaction(d):
    return {"path_to_input": d, "Geometry": "geom.dat",
            "Transfer": "transfer.dat", "CoulombIntra": "coulombintra.dat",
            "CoulombInter": "coulombinter.dat"}


def _flex_input(d, out, bond, eli=None):
    p = {"T": T, "filling": FILLING, "CellShape": [L, L, 1],
         "SubShape": [1, 1, 1], "orbital_site": [0], "Nmat": NMAT,
         "IterationMax": 500, "Mix": 0.2, "mixing_scheme": "anderson",
         "anderson_depth": 8, "EPS": 10, "flex_hartree_fock": True,
         "flex_guard_policy": "warn", "longitudinal_bond_channels": bond}
    if bond:
        p["longitudinal_bond_pairing"] = "singlet"
    inp = {"log": {"print_level": 0},
           "mode": {"mode": "FLEX", "calc_scheme": "general", "param": p},
           "file": {"input": {"path_to_input": d, "interaction": _interaction(d)},
                    "output": {"path_to_output": out, "sigma": "sigma",
                               "green": "green", "chiq": "chiq"}}}
    if eli is not None:
        inp["eliashberg"] = eli
    return inp


def _sc_input(d, flex_out, out, eli):
    e = {"chi0q_mode": "flex", "frequency": "dynamic", "bond_channels": False,
         "pairing_type": "singlet"}
    e.update(eli)
    return {"mode": {"mode": "SC", "param": {
                "T": T, "filling": FILLING, "CellShape": [L, L, 1],
                "SubShape": [1, 1, 1], "orbital_site": [0], "Nmat": NMAT}},
            "file": {"input": {"path_to_flex_output": flex_out,
                               "interaction": _interaction(d)},
                     "output": {"path_to_output": out}},
            "eliashberg": e}


def _run_flex(inp):
    """Run the FLEX solve; return the bond kernel's matvec (when the
    in-process pairing step ran) and the green_info of the run."""
    import hwave.qlms as qlms
    from hwave.solver import eliashberg_dynamic as ed
    cap = {}
    real = ed.run_leading_eigenproblem

    def spy(matvec, gap_shape, *a, **k):
        cap["matvec"], cap["shape"] = matvec, gap_shape
        return real(matvec, gap_shape, *a, **k)

    with mock.patch.object(ed, "run_leading_eigenproblem", spy):
        qlms.run(input_dict=inp)
    return cap


def _green(path):
    """G(k, iw_n) of a one-orbital green.npz as (Nx, Ny, Nz, Nmat), and beta."""
    g = np.load(path)
    G = np.zeros((L, L, 1, NMAT), complex)
    for k, (i, j, l) in enumerate(g["wavevector_index"]):
        G[i, j, l, :] = g["green"][0, :, k, 0, 0]
    return G, float(g["beta"])


def _f_loc(G, beta, phi):
    """F_loc = (1/(N beta)) sum_{k,n} G(k, iw_n) phi(k, iw_n) G(-k, -iw_n),
    through the solver's own pair bubble (calc_g2_dynamic)."""
    from hwave.solver import eliashberg_dynamic as ed
    G2 = ed.calc_g2_dynamic(G[np.newaxis, np.newaxis], beta)[0, 0, 0, 0]
    return np.sum(G2 * phi) / phi[..., 0].size


def _onsite(phi):
    """Delta(R = 0, iw_n) = (1/N) sum_k Delta(k, iw_n)."""
    return np.fft.ifftn(phi, axes=(0, 1, 2))[0, 0, 0]


def _spectrum(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if line.startswith("#") or len(p) < 4:
                continue
            rows.append(complex(float(p[1]), float(p[2])))
    return np.array(rows)


class TestRepulsiveHubbardSinglet(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        d = cls.model = os.path.join(cls.tmp, "u3")
        _write_model(d, 3.0)
        cls.bond_out = os.path.join(d, "out_bond")
        cls.nb_out = os.path.join(d, "out_nb")
        cls.sc_out = os.path.join(d, "out_sc")
        cls.cap = _run_flex(_flex_input(d, cls.bond_out, True))
        _run_flex(_flex_input(d, cls.nb_out, False))
        import hwave.sc as sc
        cls.lam_sc = sc.calc_eliashberg(_sc_input(
            d, cls.nb_out, cls.sc_out,
            {"solver_mode": "eigenvalue", "num_eigenvalues": 4}))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def _bond(self, out=None):
        return np.load(os.path.join(out or self.bond_out,
                                    "eliashberg_bond_singlet.npz"),
                       allow_pickle=True)

    # (a) ---------------------------------------------------------------
    def test_a_default_in_process_matches_non_bond_hwave_sc(self):
        zb = self._bond()
        zn = np.load(os.path.join(self.sc_out, "gap_dynamic.npz"))
        lam_b, lam_n = float(zb["eigenvalue"]), float(zn["eigenvalue"])
        self.assertIs(bool(zb["eigen_converged"]), True)
        self.assertEqual(str(zb["eigenvalue_selection"]), "LR_projected")
        self.assertAlmostEqual(lam_b, lam_n, delta=1.0e-8)
        self.assertAlmostEqual(lam_n, self.lam_sc, delta=1.0e-12)
        self.assertGreater(lam_b, 0.0)
        gb, gn = zb["gap"], zn["gap"]
        overlap = abs(np.vdot(gb, gn)) / (np.linalg.norm(gb) * np.linalg.norm(gn))
        self.assertAlmostEqual(overlap, 1.0, delta=1.0e-8)
        sb = _spectrum(os.path.join(self.bond_out, "eigenvalue_bond_singlet.dat"))
        sn = _spectrum(os.path.join(self.sc_out, "eigenvalue.dat"))
        # in-process default num_eigenvalues = 10, hwave_sc asked for 4: the
        # largest-real (projected) solve makes the top 4 identical
        self.assertEqual((len(sb), len(sn)), (10, 4))
        np.testing.assert_allclose(sb[:4].real, sn.real, rtol=0, atol=1.0e-8)

    def test_a_the_dominant_eigenvalue_is_the_negative_onsite_mode(self):
        """The premise of the defect: |lambda_-| > lambda_leading with
        lambda_- < 0 for this kernel."""
        from scipy.sparse.linalg import LinearOperator, eigs
        n = int(np.prod(self.cap["shape"]))
        A = LinearOperator((n, n), matvec=self.cap["matvec"], dtype=complex)
        v = eigs(A, k=1, which="LM", return_eigenvectors=False, tol=1e-10)[0]
        self.assertLess(v.real, -1.5)
        self.assertGreater(abs(v), 3.0 * float(self._bond()["eigenvalue"]))

    def test_a_tsweep_parser_reads_the_bond_leader(self):
        """hwave_tsweep's parser on eigenvalue_bond_singlet.dat: the parsed
        value and row 0 equal the returned leader."""
        from hwave import tsweep
        path = os.path.join(self.bond_out, "eigenvalue_bond_singlet.dat")
        lam = float(self._bond()["eigenvalue"])
        re, _im, match = tsweep.parse_leading_eig(path)
        self.assertAlmostEqual(re, lam, delta=1e-7)
        self.assertEqual(match, 1)
        rows = [l.split() for l in open(path)
                if l.strip() and not l.startswith("#") and len(l.split()) >= 4]
        self.assertEqual(int(rows[0][0]), 0)
        self.assertAlmostEqual(float(rows[0][1]), lam, delta=1e-7)

    # (b) ---------------------------------------------------------------
    def test_b_high_frequency_identity_on_the_converged_gap(self):
        """lambda Delta(R=0, iw_n) -> -U F_loc at large |w_n|: only the
        instantaneous on-site U survives there. Evaluated at the interior
        high frequency n = 3 Nmat/4 (w_n ~ 10.4 t), where the measured ratio
        is 0.924; the fluctuation part still decays, and the top of the grid
        (w ~ 19.8 t) sits at 0.82 -- so the window [0.85, 1.05] passes the
        correct sign and magnitude and fails any sign error (the reported
        defect had ratio -0.31)."""
        zb = self._bond()
        lam = float(zb["eigenvalue"])
        phi = zb["gap"][0, 0]
        G, beta = _green(os.path.join(self.bond_out, "green.npz"))
        F = _f_loc(G, beta, phi)
        self.assertGreater(abs(F), 1.0e-4 * np.max(np.abs(phi)))
        n = 3 * NMAT // 4
        ratio = (lam * _onsite(phi)[n] / (-3.0 * F)).real
        self.assertTrue(0.85 < ratio < 1.05, ratio)
        top = (lam * _onsite(phi)[NMAT - 1] / (-3.0 * F)).real
        self.assertTrue(0.75 < top < 1.05, top)

    # (c) ---------------------------------------------------------------
    def _constant_gap_ratio(self, cap, green_path, U):
        phi = np.ones(cap["shape"], dtype=complex)
        out = cap["matvec"](phi.ravel()).reshape(cap["shape"])[0, 0]
        G, beta = _green(green_path)
        F = _f_loc(G, beta, phi[0, 0])
        self.assertGreater(F.real, 0.0)
        n = 3 * NMAT // 4
        return _onsite(out)[n].real, (-U * F).real

    def test_c_constant_onsite_gap_sign_follows_minus_u(self):
        """Sign-convention pin of the bond kernel's instantaneous on-site
        term, NOT a regression test of the reporting defect: it applies the
        kernel to a constant on-site gap directly, without any eigen-solver,
        so it passes on the code before the fix as well (the kernel was
        right; the iteration's report was wrong). It guards the convention
        lambda Delta(R=0, iw -> inf) = -U F_loc, F_loc > 0, for U > 0 and
        U < 0, on which test (b) relies."""
        for U in (3.0, -1.0):
            with self.subTest(U=U):
                if U == 3.0:
                    cap, green = self.cap, os.path.join(self.bond_out, "green.npz")
                else:
                    d = os.path.join(self.tmp, "um1")
                    _write_model(d, U)
                    cap = _run_flex(_flex_input(d, os.path.join(d, "out"), True))
                    green = os.path.join(d, "out", "green.npz")
                got, want = self._constant_gap_ratio(cap, green, U)
                self.assertEqual(np.sign(got), np.sign(want))
                self.assertEqual(np.sign(got), -np.sign(U))
                # measured 0.977 (U = 3) and 0.987 (U = -1) at n = 3 Nmat/4
                self.assertTrue(0.85 < got / want < 1.15, got / want)

    # (d) ---------------------------------------------------------------
    def test_d_unshifted_iteration_reports_nan_everywhere(self):
        d = self.model
        out = os.path.join(d, "out_iter_none")
        import hwave.qlms as qlms
        from hwave.solver import flex
        captured = {}
        real_write = flex.FLEX.save_results

        def save(self_, info_outputfile, green_info):
            captured.update({k: v for k, v in green_info.items()
                             if str(k).startswith("pairing_")})
            return real_write(self_, info_outputfile, green_info)

        with mock.patch.object(flex.FLEX, "save_results", save):
            qlms.run(input_dict=_flex_input(
                d, out, True,
                eli={"solver_mode": "iteration", "spectral_shift": "none"}))
        self.assertTrue(np.isnan(captured["pairing_singlet_eigenvalue"]))
        self.assertIs(captured["pairing_singlet_eigen_converged"], False)
        z = self._bond(out)
        self.assertTrue(np.isnan(float(z["eigenvalue"])))
        self.assertIs(bool(z["eigen_converged"]), False)
        with open(os.path.join(out, "eigenvalue_bond_singlet.dat")) as fh:
            text = fh.read()
        self.assertIn("# eigen_converged=false", text)
        self.assertIn("two-cycle", text)
        vals = [l for l in text.splitlines() if not l.startswith("#")]
        self.assertTrue(np.isnan(float(vals[0])))
        with open(os.path.join(out, "gap_bond_singlet.dat")) as fh:
            self.assertIn("eigen_converged=false", fh.readline())
        # the FLEX outputs themselves were written
        self.assertTrue(os.path.exists(os.path.join(out, "green.npz")))

    def test_d_default_shift_iteration_converges_to_the_leading_eigenvalue(self):
        d = self.model
        out = os.path.join(d, "out_iter_auto")
        _run_flex(_flex_input(
            d, out, True, eli={"solver_mode": "iteration", "max_iter": 5000,
                               "convergence_tol": 1.0e-11}))
        z = self._bond(out)
        self.assertIs(bool(z["eigen_converged"]), True)
        self.assertAlmostEqual(float(z["eigenvalue"]),
                               float(self._bond()["eigenvalue"]), delta=1.0e-8)
        with open(os.path.join(out, "eigenvalue_bond_singlet.dat")) as fh:
            self.assertIn("spectral_shift='auto'", fh.read())

    # hwave_sc exit status (1d) ----------------------------------------------
    def _hwave_sc(self, out, solver_lines):
        d = self.model
        toml = os.path.join(self.tmp, "sc_{}.toml".format(os.path.basename(out)))
        with open(toml, "w") as fh:
            fh.write(
                "[mode]\nmode = \"SC\"\n[mode.param]\nT = {T}\nfilling = {f}\n"
                "CellShape = [{L}, {L}, 1]\nSubShape = [1, 1, 1]\n"
                "orbital_site = [0]\nNmat = {N}\n"
                "[file.input]\npath_to_flex_output = \"{nb}\"\n"
                "[file.input.interaction]\npath_to_input = \"{d}\"\n"
                "Geometry = \"geom.dat\"\nTransfer = \"transfer.dat\"\n"
                "CoulombIntra = \"coulombintra.dat\"\n"
                "CoulombInter = \"coulombinter.dat\"\n"
                "[file.output]\npath_to_output = \"{out}\"\n"
                "[eliashberg]\nchi0q_mode = \"flex\"\nfrequency = \"dynamic\"\n"
                "bond_channels = false\npairing_type = \"singlet\"\n{extra}"
                .format(T=T, f=FILLING, L=L, N=NMAT, nb=self.nb_out, d=d,
                        out=out, extra=solver_lines))
        src = os.path.join(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))), "src")
        env = dict(os.environ)
        env["PYTHONPATH"] = src + os.pathsep + env.get("PYTHONPATH", "")
        return subprocess.run(
            [sys.executable, "-B", "-c",
             "import sys; from hwave.sc import main; sys.argv[0] = 'hwave_sc'; main()",
             "-q", toml], env=env, capture_output=True, text=True)

    def test_hwave_sc_exit_status(self):
        import hwave.sc as sc
        out = os.path.join(self.tmp, "sc_none")
        r = self._hwave_sc(out, "solver_mode = \"iteration\"\n"
                                "spectral_shift = \"none\"\n")
        self.assertEqual(r.returncode, sc.EXIT_EIGEN_NOT_CONVERGED, r.stderr)
        self.assertIn("did not converge", r.stderr)
        with open(os.path.join(out, "eigenvalue.dat")) as fh:
            text = fh.read()
        self.assertIn("# eigen_converged=false", text)
        self.assertIs(bool(np.load(os.path.join(out, "gap_dynamic.npz"))
                           ["eigen_converged"]), False)
        ok = os.path.join(self.tmp, "sc_eig")
        r = self._hwave_sc(ok, "solver_mode = \"eigenvalue\"\n"
                               "num_eigenvalues = 4\n")
        self.assertEqual(r.returncode, 0, r.stderr)


if __name__ == "__main__":
    unittest.main()
