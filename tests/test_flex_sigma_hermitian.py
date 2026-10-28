"""The frequency-reversal Hermitian symmetry of the FLEX self-energy
(GitHub issue #281):

    Sigma(k, i w) = Sigma(k, -i w)^dag      (at fixed k)

holds for every self-energy of a Hermitian problem. The FLEX map keeps it,
but it does not damp an asymmetric part that is already there: a cold
start past the bare Stoner instability amplifies the round-off asymmetry
with the transient, and linear mixing removes it only slowly. The solver
therefore projects every update (``[mode.param]
flex_sigma_hermitian_projection``, default true), and ``hwave_sigma_split``
symmetrizes the total archive it converts.

What this module pins:

* :class:`TestProjectionHelpers` -- the projection against an explicit
  element-wise formula (the dagger is a conjugate AND a transpose, and the
  frequency index is reversed), idempotence, the deviation measure, and the
  split state (the singleton static axis gives the Hermitian part).
* :class:`TestMapPreservesSymmetry` -- with the projection OFF, a few
  iterations from a cold start keep the relation at round-off on every
  path: reduced and general vertices on the uniform grid and the IR basis,
  a trans_mod mean field, the Hartree-Fock (Phase B) loop, the spin-orbital
  loop with spin-orbit coupling, and the bond gate. This is why the
  projection cannot move a fixed point.
* :class:`TestColdStartReproduction` -- the issue: a two-orbital dimer
  (``orbital_site = [0, 1]``, no inversion) past the bare Stoner
  instability, cold start, linear mixing. Without the projection the
  archive violates the relation by far more than 1e-8 (and
  ``hwave_sigma_split`` refuses it); with the projection it holds exactly,
  the run converges, and the archive is accepted as ``sigma_init`` by a
  ``flex_hartree_fock = true`` run.
* :class:`TestFixedPointUnchanged` -- a converged symmetric run gives the
  same self-energy with and without the projection.
* :class:`TestSeedProjection` -- a seed with an asymmetric remnant: the
  projection removes it at load time (Phase B included); without it the
  density guard refuses the seed, as reported in the issue.
* :class:`TestSigmaSplitSymmetrization` -- the conversion symmetrizes,
  reports the deviation, and refuses above ``--asymmetry-tol``.
* :class:`TestConfigKey` -- strict boolean parsing of the key.
"""
import contextlib
import logging
import os
import shutil
import tempfile
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k
from hwave.solver import flex_mixing as fm
from tests.heavy_tests import heavy

try:
    import sparse_ir  # noqa: F401
    _HAVE_SPARSE_IR = True
except ImportError:
    _HAVE_SPARSE_IR = False


@contextlib.contextmanager
def _quiet():
    prev = logging.root.manager.disable
    logging.disable(logging.WARNING)
    try:
        yield
    finally:
        logging.disable(prev)


def _rel_dev(sigma):
    """Reference deviation measure, written out independently of the
    production helper: max_{k,n,a,b} |S_ab(k, w_n) - conj(S_ba(k, w_{-n}))|
    relative to max |S|."""
    s = np.asarray(sigma)
    nf = s.shape[1]
    worst = 0.0
    for n in range(nf):
        m = nf - 1 - n
        d = s[:, n] - np.conj(np.swapaxes(s[:, m], -1, -2))
        worst = max(worst, float(np.abs(d).max()))
    scale = float(np.abs(s).max())
    return worst / scale if scale > 0 else worst


def _max_im_eig(sigma):
    """Largest eigenvalue of the anti-Hermitian part
    ``(Sigma - Sigma^dag) / 2i`` of ``Sigma(k, i w_n)`` over ``w_n > 0`` (the
    upper half of the centred grid) and every k. Causality of the
    self-energy requires it to be <= 0."""
    s = np.asarray(sigma)
    pos = s[:, s.shape[1] // 2:]
    a = (pos - np.conj(np.swapaxes(pos, -1, -2))) / 2j
    return float(np.linalg.eigvalsh(a).max())


def _antisymmetric(shape, amp, seed=1):
    """A random array ``a`` with ``S a = -a`` (``S: x -> x(-i w)^dag``) and
    ``max |a| = amp``."""
    rng = np.random.default_rng(seed)
    x = rng.normal(size=shape) + 1j * rng.normal(size=shape)
    a = 0.5 * (x - fm.frequency_reversed_dagger(x))
    return a * (amp / float(np.abs(a).max()))


# --------------------------------------------------------------------------
# the two-orbital dimer of the issue (two inequivalent sites, no inversion,
# complex inter-site hopping so that no time-reversal relation helps)
# --------------------------------------------------------------------------

def _w90(path, title, norb, rows):
    rv = sorted(set(r[0] for r in rows))
    with open(path, "w") as f:
        f.write("{}\n{}\n{}\n".format(title, norb, len(rv)) + " ".join(["1"] * len(rv)) + "\n")
        for R, (a, b), v in rows:
            v = complex(v)
            f.write("%d %d %d %d %d %.15f %.15f\n" % (R[0], R[1], R[2], a, b, v.real, v.imag))


def _write_dimer(d, U, V=None):
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1 0 0\n0 1 0\n0 0 1\n2\n0 0 0\n0.3 0.4 0\n")
    t1 = -0.6 * np.exp(0.3j)
    t2 = -0.3 * np.exp(-0.5j)
    rows = []
    for R in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)):
        rows.append((R, (1, 1), -1.0))
        rows.append((R, (2, 2), -0.8))
    rows += [((0, 0, 0), (1, 2), t1), ((0, 0, 0), (2, 1), np.conj(t1)),
             ((-1, 0, 0), (1, 2), t2), ((1, 0, 0), (2, 1), np.conj(t2)),
             ((0, -1, 0), (1, 2), -0.2), ((0, 1, 0), (2, 1), -0.2),
             ((0, 0, 0), (1, 1), 0.1), ((0, 0, 0), (2, 2), -0.1)]
    _w90(os.path.join(d, "transfer.dat"), "Transfer", 2, rows)
    _w90(os.path.join(d, "coulombintra.dat"), "CoulombIntra", 2,
         [((0, 0, 0), (1, 1), U), ((0, 0, 0), (2, 2), U)])
    files = {"CoulombIntra": "coulombintra.dat"}
    if V is not None:
        _w90(os.path.join(d, "coulombinter.dat"), "CoulombInter", 2,
             [((0, 0, 0), (1, 2), V), ((0, 0, 0), (2, 1), V)])
        files["CoulombInter"] = "coulombinter.dat"
    return files


def _dimer_solver(d, files, scheme="reduced", **par):
    import hwave.solver.flex as flex_mod
    idict = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
    idict.update(files)
    r = read_input_k.QLMSkInput({"path_to_input": d, "interaction": idict})
    p = {"T": 0.1, "filling": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1],
         "Nmat": 16, "IterationMax": 4, "Mix": 0.5, "EPS": 14, "orbital_site": [0, 1]}
    p.update(par)
    with _quiet():
        s = flex_mod.FLEX(r.get_param("ham"), {}, {
            "mode": "FLEX", "param": p, "enable_spin_orbital": False, "calc_scheme": scheme})
    return s, r.get_param("green")


def _solve(s, gi, d):
    with _quiet():
        s.solve(gi, d)
    return s


class _TmpDir(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="flex_herm_")
        self.addCleanup(shutil.rmtree, self.d, True)


class TestProjectionHelpers(unittest.TestCase):

    def setUp(self):
        rng = np.random.default_rng(7)
        shape = (2, 6, 5, 3, 3)
        self.s = rng.normal(size=shape) + 1j * rng.normal(size=shape)

    def test_reversed_dagger_elementwise(self):
        r = fm.frequency_reversed_dagger(self.s)
        nf = self.s.shape[1]
        for b in range(2):
            for n in range(nf):
                for k in range(5):
                    np.testing.assert_array_equal(
                        r[b, n, k], self.s[b, nf - 1 - n, k].conj().T)

    def test_projection_is_the_symmetric_part(self):
        p, (a, rel) = fm.project_frequency_hermitian(self.s, with_deviation=True)
        nf = self.s.shape[1]
        for n in range(nf):
            ref = 0.5 * (self.s[:, n] + np.conj(np.swapaxes(self.s[:, nf - 1 - n], -1, -2)))
            np.testing.assert_allclose(p[:, n], ref, rtol=0, atol=1e-15)
        # exactly symmetric, and idempotent bit for bit
        self.assertEqual(_rel_dev(p), 0.0)
        np.testing.assert_array_equal(fm.project_frequency_hermitian(p), p)
        # the deviation measure of the input
        self.assertAlmostEqual(rel, _rel_dev(self.s), places=14)
        self.assertAlmostEqual(a, rel * float(np.abs(self.s).max()), places=12)
        self.assertEqual(fm.frequency_hermitian_deviation(p), (0.0, 0.0))
        z = np.zeros_like(self.s)
        self.assertEqual(fm.frequency_hermitian_deviation(z), (0.0, 0.0))

    def test_inplace_pairwise_equals_the_formula(self):
        for nf in (6, 5, 1):                     # even, odd (IR-like), static
            s = self.s[:, :nf].copy()
            ref = 0.5 * (s + fm.frequency_reversed_dagger(s))
            out, (a, rel) = fm.project_frequency_hermitian_inplace(s.copy(), with_deviation=True)
            np.testing.assert_array_equal(out, ref)          # bit for bit
            full = float(np.abs(s - fm.frequency_reversed_dagger(s)).max())
            self.assertEqual(a, full)
            self.assertEqual(fm.frequency_hermitian_deviation(s)[0], full)
            st = self.s[:, :1].copy()
            tot = s + st
            self.assertEqual(fm.frequency_hermitian_deviation(s, static=st)[0],
                             float(np.abs(tot - fm.frequency_reversed_dagger(tot)).max()))

    def test_split_state(self):
        st = fm.SplitState(static=self.s[:1, :1].copy(), fluct=self.s[:1].copy())
        out = fm.project_split_state(st)
        np.testing.assert_allclose(out.static[0, 0],
                                   0.5 * (st.static[0, 0] + np.conj(np.swapaxes(st.static[0, 0], -1, -2))),
                                   rtol=0, atol=1e-15)
        np.testing.assert_array_equal(out.fluct, fm.project_frequency_hermitian(st.fluct))
        self.assertEqual(_rel_dev(out.total()), 0.0)


class TestMapPreservesSymmetry(_TmpDir):
    """Projection OFF: the map itself keeps the relation at round-off."""

    TOL = 1e-12

    def _check(self, s, gi):
        _solve(s, gi, self.d)
        self.assertGreater(float(np.abs(s.sigma).max()), 1e-3)     # a nontrivial sigma
        self.assertLess(_rel_dev(s.sigma), self.TOL)
        return s

    def test_reduced_uniform(self):
        files = _write_dimer(self.d, 1.5)
        self._check(*_dimer_solver(self.d, files, flex_sigma_hermitian_projection=False))

    def test_general_uniform(self):
        files = _write_dimer(self.d, 1.5, V=0.4)
        self._check(*_dimer_solver(self.d, files, scheme="general",
                                   flex_sigma_hermitian_projection=False))

    def test_trans_mod_mean_field(self):
        files = _write_dimer(self.d, 1.5)
        s, gi = _dimer_solver(self.d, files, flex_sigma_hermitian_projection=False)
        nvol, norb = s.lattice.nvol, s.norb
        shift = np.array([[0.2, 0.1 + 0.15j], [0.1 - 0.15j, -0.1]])
        H = np.zeros((nvol, 2, norb, 2, norb), complex)
        base = np.asarray(s.ham_info.ham_trans_q).reshape(nvol, norb, norb)
        H[:, 0, :, 0, :] = base + shift
        H[:, 1, :, 1, :] = base + shift
        gi["trans_mod"] = H.reshape(nvol, 2 * norb, 2 * norb)
        self._check(s, gi)

    def test_phase_b_hartree_fock(self):
        files = _write_dimer(self.d, 1.5, V=0.4)
        self._check(*_dimer_solver(self.d, files, scheme="general", flex_hartree_fock=True,
                                   flex_sigma_hermitian_projection=False))

    @unittest.skipUnless(_HAVE_SPARSE_IR, "sparse-ir not installed")
    def test_reduced_ir(self):
        files = _write_dimer(self.d, 1.5)
        s, gi = _dimer_solver(self.d, files, matsubara_basis="ir", ir_wmax=20.0,
                              flex_sigma_hermitian_projection=False, write_densified=False)
        self._check(s, gi)

    @unittest.skipUnless(_HAVE_SPARSE_IR, "sparse-ir not installed")
    def test_projection_on_ir_nodes(self):
        """Projection ON on the IR path: the node values are symmetric
        exactly after every update (the map alone leaves round-off), and
        the densified output inherits it."""
        files = _write_dimer(self.d, 1.5)
        s, gi = _dimer_solver(self.d, files, matsubara_basis="ir", ir_wmax=20.0,
                              write_densified=False)
        _solve(s, gi, self.d)
        self.assertEqual(s.sigma.shape[1], s._ir_axF.n_freq)
        self.assertGreater(s.sigma_hermitian_removed_max, 0.0)   # round-off was removed
        self.assertEqual(_rel_dev(s.sigma), 0.0)
        s, gi = _dimer_solver(self.d, files, matsubara_basis="ir", ir_wmax=20.0, Nmat=256)
        _solve(s, gi, self.d)
        self.assertEqual(s.sigma.shape[1], 256)
        self.assertLess(_rel_dev(s.sigma), 1e-13)

    @unittest.skipUnless(_HAVE_SPARSE_IR, "sparse-ir not installed")
    def test_general_ir_densified(self):
        files = _write_dimer(self.d, 1.5, V=0.4)
        s, gi = _dimer_solver(self.d, files, scheme="general", matsubara_basis="ir",
                              ir_wmax=20.0, Nmat=256, flex_sigma_hermitian_projection=False)
        self._check(s, gi)

    def test_spin_orbital_with_soc(self):
        from tests import test_flex_spin_orbital as so
        files = so._write_model(self.d, True, soc=0.2)
        s, gi = so._build(self.d, files, IterationMax=4, flex_sigma_hermitian_projection=False)
        self._check(s, gi)

    def test_bond_gate(self):
        from tests import test_flex_bond_gate as bg
        s, r = bg._flex({"IterationMax": 4, "flex_sigma_hermitian_projection": False})
        self._check(s, r.get_param("green"))


#: the cold-start configuration of the issue on the dimer: U past the bare
#: Stoner instability, small linear mixing. The projected cold start crosses
#: the pole 1 - U chi0 = 0 of the spin channel; with Mix = 0.05 an unlucky
#: step there left it unconverged after 1000 iterations in about one run in
#: ten under round-off-level changes of U (GitHub issue #299). With
#: Mix = 0.06 both runs take the same number of iterations under relative
#: changes of U up to 1e-6, and the unprojected run still reaches the broken
#: branch (larger Mix does not: the unprojected run then stays symmetric).
_COLD = dict(T=0.1, CellShape=[6, 6, 1], Nmat=32, Mix=0.06, EPS=8, IterationMax=1000)
_COLD_U = 2.8


class TestColdStartReproduction(_TmpDir):

    @heavy
    def test_cold_start_remnant_and_split_seed(self):
        from hwave.sigma_split import main as split_main
        files = _write_dimer(self.d, _COLD_U)
        # without the projection the run CONVERGES -- to a symmetry-broken,
        # non-causal branch (relative violation ~0.3), with only the
        # end-of-solve diagnostic to tell
        s_off, gi = _dimer_solver(self.d, files, flex_sigma_hermitian_projection=False, **_COLD)
        with self.assertLogs("hwave.solver.flex", level="WARNING") as cm:
            s_off.solve(gi, self.d)
        self.assertTrue(s_off.scf_converged)
        dev_off = _rel_dev(s_off.sigma)
        self.assertGreater(dev_off, 0.1)
        self.assertGreater(_max_im_eig(s_off.sigma), 0.05)
        self.assertTrue(any("non-causal" in m for m in cm.output))
        s_off.save_results({"path_to_output": self.d, "sigma": "off.npz"}, gi)
        # with the projection (the default): exact, converged and causal
        s_on, gi = _dimer_solver(self.d, files, **_COLD)
        _solve(s_on, gi, self.d)
        self.assertTrue(s_on.scf_converged)
        self.assertLessEqual(_rel_dev(s_on.sigma), 1e-14)
        self.assertLess(_max_im_eig(s_on.sigma), 0.0)
        s_on.save_results({"path_to_output": self.d, "sigma": "on.npz"}, gi)
        with _quiet():
            # the unprojected archive is refused by the converter and flagged
            # by the audit; the projected one passes both
            self.assertEqual(split_main([os.path.join(self.d, "off.npz"),
                                         os.path.join(self.d, "off_split.npz"),
                                         "--zero-static"]), 1)
            self.assertEqual(split_main(["--check", os.path.join(self.d, "off.npz")]), 1)
            self.assertEqual(split_main(["--check", os.path.join(self.d, "on.npz")]), 0)
            self.assertEqual(split_main([os.path.join(self.d, "on.npz"),
                                         os.path.join(self.d, "seed.npz"), "--zero-static"]), 0)
        # the split seed is accepted by a flex_hartree_fock run, even with the
        # consumer's own seed projection switched off
        s_hf, gi = _dimer_solver(self.d, files, scheme="general", flex_hartree_fock=True,
                                 flex_sigma_hermitian_projection=False,
                                 **dict(_COLD, IterationMax=2))
        gi.update(s_hf.read_init({"path_to_input": self.d, "sigma_init": "seed.npz"}))
        _solve(s_hf, gi, self.d)
        self.assertEqual(s_hf.scf_iterations, 2)


class TestBrokenBranchLight(_TmpDir):
    """The fast-gate twin of the cold-start case: a smaller grid and a seed
    carrying a 1e-2 anti-symmetric part. Without the projection the run
    converges to the symmetry-broken, non-causal branch and the end-of-solve
    diagnostic warns; with it the same start stays symmetric and causal.

    Both runs start past the bare Stoner instability, so the symmetric
    iteration has to cross the pole of the spin channel, 1 - U chi0 = 0, on
    its way to the fixed point: a step that lands close to it produces a
    huge transient self-energy, after which linear mixing wanders
    chaotically, and whether it settles within IterationMax is decided by
    rounding (GitHub issue #299: from the bare anti-symmetric seed with
    Mix = 0.1, relative changes of U of 1e-15 left the projected run
    unconverged in about a third of the cases and non-causal in a few
    percent). The seed therefore also carries a small causal damping
    ``-i gamma sign(w_n)``, removed by neither projection, and the mixing is
    chosen so that the projected run crosses the pole in one mild step: it
    then converges in the same number of iterations under relative changes
    of U up to 1e-6, on Accelerate and OpenBLAS alike, while the unprojected
    run still reaches the broken branch."""

    _PAR = dict(T=0.1, CellShape=[6, 6, 1], Nmat=16, Mix=0.1075, EPS=8, IterationMax=400)
    _GAMMA = 0.0875

    def _run(self, proj):
        files = _write_dimer(self.d, _COLD_U)
        s, gi = _dimer_solver(self.d, files, flex_sigma_hermitian_projection=proj, **self._PAR)
        shape = (1, 16, 36, 2, 2)
        # centred grid: the first half of the frequency axis is w_n < 0
        sign = np.where(np.arange(shape[1]) < shape[1] // 2, -1.0, 1.0)
        damping = (-1j * self._GAMMA * sign)[None, :, None, None, None] * np.eye(2)
        gi["sigma_init"] = damping + _antisymmetric(shape, 1e-2)
        with self.assertLogs("hwave.solver.flex", level="INFO") as cm:
            s.solve(gi, self.d)
        return s, [m for m in cm.output if m.startswith("WARNING")]

    def test_unprojected_converges_to_the_broken_branch(self):
        s, warn = self._run(False)
        self.assertTrue(s.scf_converged)
        self.assertGreater(_rel_dev(s.sigma), 0.1)
        self.assertGreater(s.sigma_final_hermitian_deviation, 0.1)
        self.assertGreater(_max_im_eig(s.sigma), 0.01)
        # the seed warning and the end-of-solve warning, both naming the branch
        self.assertTrue(any("sigma_init violates" in m and "non-causal" in m for m in warn))
        self.assertTrue(any("final self-energy violates" in m and "non-causal" in m
                            for m in warn))

    def test_projected_stays_symmetric_and_causal(self):
        s, warn = self._run(True)
        self.assertTrue(s.scf_converged)
        self.assertEqual(_rel_dev(s.sigma), 0.0)
        self.assertLess(_max_im_eig(s.sigma), 0.0)
        self.assertFalse(any("final self-energy violates" in m for m in warn))


class TestDiagnosticsNeverFailTheSolve(_TmpDir):
    """The deviation diagnostics (end of solve with the projection off,
    seed check) must not fail a solve: a failure inside them -- e.g. device
    memory exhausted after convergence -- is logged as a skipped check and
    recorded as NaN, and the results are kept."""

    def _boom(self, *a, **k):
        raise MemoryError("simulated out-of-memory in the diagnostic")

    def _solve_patched(self, s, gi):
        from unittest import mock
        with mock.patch.object(fm, "frequency_hermitian_deviation", self._boom), \
                self.assertLogs("hwave.solver.flex", level="WARNING") as cm:
            s.solve(gi, self.d)
        return cm.output

    def test_end_of_solve_legacy_and_phase_b(self):
        files = _write_dimer(self.d, 1.5, V=0.4)
        for extra in ({}, {"scheme": "general", "flex_hartree_fock": True}):
            with self.subTest(**{k: str(v) for k, v in extra.items()}):
                s, gi = _dimer_solver(self.d, files, flex_sigma_hermitian_projection=False,
                                      IterationMax=3, **extra)
                out = self._solve_patched(s, gi)
                self.assertTrue(any("check of the final self-energy was skipped" in m
                                    and "MemoryError" in m for m in out))
                self.assertTrue(np.isnan(s.sigma_final_hermitian_deviation))
                self.assertEqual(s.scf_iterations, 3)
                self.assertIn("sigma", gi)
                self.assertGreater(float(np.abs(s.sigma).max()), 1e-3)

    def test_seed_check(self):
        files = _write_dimer(self.d, 1.5)
        s, gi = _dimer_solver(self.d, files, IterationMax=1)
        gi["sigma_init"] = _antisymmetric((1, 16, 16, 2, 2), 1e-2)
        out = self._solve_patched(s, gi)
        self.assertTrue(any("check of the seed was skipped" in m for m in out))
        self.assertTrue(np.isnan(s.sigma_seed_hermitian_deviation))
        self.assertEqual(_rel_dev(s.sigma), 0.0)     # the projection itself still ran


class TestIRNativeArchiveCheck(_TmpDir):

    @unittest.skipUnless(_HAVE_SPARSE_IR, "sparse-ir not installed")
    def test_check_and_mixed_metadata(self):
        from hwave.sigma_split import main
        from hwave.solver.flex_hf import sigma_hermitian_check
        files = _write_dimer(self.d, 1.5)
        s, gi = _dimer_solver(self.d, files, matsubara_basis="ir", ir_wmax=20.0,
                              write_densified=False)
        _solve(s, gi, self.d)
        s.save_results({"path_to_output": self.d, "sigma": "ir.npz"}, gi)
        path = os.path.join(self.d, "ir.npz")
        data = dict(np.load(path))
        self.assertNotIn("freq_index", data)
        self.assertEqual(sigma_hermitian_check(path)["grid"], "ir")
        with _quiet():
            self.assertEqual(main(["--check", path]), 0)
        # mixed metadata: IR-native nodes AND a uniform freq_index
        mixed = os.path.join(self.d, "mixed.npz")
        np.savez(mixed, **dict(data, freq_index=np.arange(data["sigma"].shape[1])))
        with self.assertRaisesRegex(ValueError, "both IR-native node metadata and a uniform"):
            sigma_hermitian_check(mixed)
        # inconsistent: a node set that is not reversal-symmetric
        bad = os.path.join(self.d, "bad.npz")
        n = np.array(data["ir_freq_n"]).copy()
        n[-1] += 2
        np.savez(bad, **dict(data, ir_freq_n=n))
        with self.assertRaisesRegex(ValueError, "not reversal-symmetric"):
            sigma_hermitian_check(bad)
        with _quiet():
            self.assertEqual(main(["--check", mixed]), 1)
            self.assertEqual(main(["--check", bad]), 1)


class TestEquivariance(_TmpDir):
    """F(S Sigma) = S F(Sigma) for one unprojected map step from a random
    asymmetric seed (S: Sigma -> Sigma(k, -i w)^dag), on the dimer without
    inversion and with complex inter-orbital hopping: the symmetric
    subspace is invariant, and an asymmetric fixed point comes with its S
    image."""

    def _map(self, files, scheme, seed):
        s, gi = _dimer_solver(self.d, files, scheme=scheme, IterationMax=1, Mix=1.0,
                              flex_sigma_hermitian_projection=False)
        gi["sigma_init"] = seed
        _solve(s, gi, self.d)
        return np.asarray(s.sigma)

    def test_one_map_step(self):
        for scheme, V in (("reduced", None), ("general", 0.4)):
            with self.subTest(scheme=scheme):
                files = _write_dimer(self.d, 1.5, V=V)
                rng = np.random.default_rng(5)
                shape = (1, 16, 16, 2, 2)
                sym = rng.normal(size=shape) + 1j * rng.normal(size=shape)
                sym = 0.05 * fm.project_frequency_hermitian(sym) / np.abs(sym).max()
                seed = sym + _antisymmetric(shape, 0.02, seed=9)
                f = self._map(files, scheme, seed)
                fs = self._map(files, scheme, fm.frequency_reversed_dagger(seed))
                self.assertGreater(_rel_dev(f), 1e-3)          # not vacuous
                scale = float(np.abs(f).max())
                np.testing.assert_allclose(fs, fm.frequency_reversed_dagger(f), rtol=0,
                                           atol=1e-12 * scale)


class TestFixedPointUnchanged(_TmpDir):

    def test_same_sigma_with_and_without_projection(self):
        files = _write_dimer(self.d, 1.5)
        out = {}
        for proj in (False, True):
            s, gi = _dimer_solver(self.d, files, flex_sigma_hermitian_projection=proj,
                                  IterationMax=200, EPS=1e-13, mixing_scheme="anderson")
            _solve(s, gi, self.d)
            self.assertTrue(s.scf_converged)
            out[proj] = (s.sigma, s.mu)
        scale = float(np.abs(out[True][0]).max())
        np.testing.assert_allclose(out[True][0], out[False][0], rtol=0, atol=1e-11 * scale)
        self.assertAlmostEqual(out[True][1], out[False][1], places=11)


def _asymmetric_seed(sigma, rel, seed=11):
    """``sigma`` plus an anti-symmetric perturbation of relative size ``rel``
    (``P(x) = (x - x(-i w)^dag) / 2``, so the projection removes it exactly)."""
    rng = np.random.default_rng(seed)
    x = rng.normal(size=sigma.shape) + 1j * rng.normal(size=sigma.shape)
    a = 0.5 * (x - fm.frequency_reversed_dagger(x))
    a *= rel * float(np.abs(sigma).max()) / float(np.abs(a).max())
    return sigma + a


class TestSeedProjection(_TmpDir):
    """A total archive with a 1e-5 asymmetric remnant (the size reported in
    the issue), converted to a split seed without symmetrization."""

    def _split_seed(self, files, rel):
        s, gi = _dimer_solver(self.d, files, scheme="general", IterationMax=20)
        _solve(s, gi, self.d)
        s.save_results({"path_to_output": self.d, "sigma": "total.npz"}, gi)
        data = dict(np.load(os.path.join(self.d, "total.npz")))
        sig = _asymmetric_seed(data["sigma"], rel)
        self.assertGreater(_rel_dev(sig), 0.5 * rel)
        data.update(sigma=sig, sigma_convention=np.str_("split"),
                    sigma_static=np.zeros_like(sig[:, :1]), sigma_fluct=sig)
        np.savez(os.path.join(self.d, "seed.npz"), **data)

    def _consume(self, files, proj):
        s, gi = _dimer_solver(self.d, files, scheme="general", flex_hartree_fock=True,
                              IterationMax=1, flex_sigma_hermitian_projection=proj)
        gi.update(s.read_init({"path_to_input": self.d, "sigma_init": "seed.npz"}))
        _solve(s, gi, self.d)
        return s

    def test_phase_b_seed(self):
        files = _write_dimer(self.d, 1.5)
        self._split_seed(files, 1e-5)
        with self.assertRaisesRegex(ValueError, r"violates rho_ab\(r\) = conj\(rho_ba\(-r\)\)"):
            self._consume(files, False)
        s = self._consume(files, True)
        self.assertGreater(s.sigma_seed_hermitian_deviation, 1e-6)
        self.assertLess(s.sigma_seed_hermitian_deviation, 1e-3)
        self.assertEqual(_rel_dev(s.sigma), 0.0)

    def test_legacy_seed_and_large_asymmetry_warns(self):
        files = _write_dimer(self.d, 1.5)
        self._split_seed(files, 1e-2)
        data = dict(np.load(os.path.join(self.d, "seed.npz")))
        for k in ("sigma_convention", "sigma_static", "sigma_fluct"):
            data.pop(k)
        np.savez(os.path.join(self.d, "seed.npz"), **data)
        s, gi = _dimer_solver(self.d, files, IterationMax=1)
        gi.update(s.read_init({"path_to_input": self.d, "sigma_init": "seed.npz"}))
        with self.assertLogs("hwave.solver.flex", level="WARNING") as cm:
            s.solve(gi, self.d)
        self.assertTrue(any("not the self-energy of this Hermitian problem" in m
                            and "non-causal" in m for m in cm.output))
        self.assertEqual(_rel_dev(s.sigma), 0.0)


class TestSigmaSplitSymmetrization(_TmpDir):

    def setUp(self):
        super().setUp()
        files = _write_dimer(self.d, 1.5)
        s, gi = _dimer_solver(self.d, files, IterationMax=20)
        _solve(s, gi, self.d)
        s.save_results({"path_to_output": self.d, "sigma": "sym.npz"}, gi)
        data = dict(np.load(os.path.join(self.d, "sym.npz")))
        self.sym = data["sigma"]
        for name, rel in (("small", 1e-5), ("large", 1e-2)):
            np.savez(os.path.join(self.d, name + ".npz"),
                     **dict(data, sigma=_asymmetric_seed(self.sym, rel)))

    def _convert(self, name, **kw):
        from hwave.solver.flex_hf import sigma_split_convert
        with _quiet():
            return sigma_split_convert(os.path.join(self.d, name + ".npz"),
                                       os.path.join(self.d, name + "_out.npz"),
                                       zero_static=True, force=True, **kw)

    def test_symmetrizes_and_reports(self):
        summ = self._convert("small")
        self.assertGreater(summ["asymmetry"], 0.5e-5)
        self.assertLess(summ["asymmetry"], 1e-3)
        out = np.load(os.path.join(self.d, "small_out.npz"))
        self.assertEqual(_rel_dev(out["sigma"]), 0.0)
        # the remnant is removed exactly: back to the symmetric archive
        np.testing.assert_allclose(out["sigma"], self.sym, rtol=0,
                                   atol=1e-14 * float(np.abs(self.sym).max()))
        np.testing.assert_array_equal(out["sigma_fluct"], out["sigma"])
        # a symmetric archive reports round-off only
        self.assertLess(self._convert("sym")["asymmetry"], 1e-13)

    def test_threshold(self):
        with self.assertRaisesRegex(ValueError, r"relative deviation of .* > 1\.0e-03"):
            self._convert("large")
        self.assertFalse(os.path.exists(os.path.join(self.d, "large_out.npz")))
        summ = self._convert("large", asymmetry_tol=0.1)
        self.assertGreater(summ["asymmetry"], 0.5e-2)

    def test_check_mode(self):
        from hwave.sigma_split import main
        with self.assertLogs("hwave_sigma_split", level="INFO") as cm:
            self.assertEqual(main(["--check", os.path.join(self.d, "small.npz")]), 0)
        self.assertTrue(any("relative deviation" in m and "uniform grid" in m for m in cm.output))
        with _quiet():
            self.assertEqual(main(["--check", os.path.join(self.d, "large.npz")]), 1)
            self.assertEqual(main(["--check", os.path.join(self.d, "large.npz"),
                                   "--asymmetry-tol", "0.1"]), 0)
            # --check writes nothing and takes no conversion options
            self.assertEqual(main(["--check", os.path.join(self.d, "small.npz"),
                                   os.path.join(self.d, "x.npz")]), 1)
            self.assertEqual(main(["--check", os.path.join(self.d, "small.npz"),
                                   "--zero-static"]), 1)
            self.assertFalse(os.path.exists(os.path.join(self.d, "x.npz")))

    def test_frequency_metadata_is_checked(self):
        data = dict(np.load(os.path.join(self.d, "small.npz")))
        nf = data["sigma"].shape[1]
        self.assertTrue(np.array_equal(data["freq_index"], np.arange(nf)))
        # a restricted uniform grid
        np.savez(os.path.join(self.d, "restricted.npz"), **dict(data, freq_index=np.arange(nf - 1)))
        with self.assertRaisesRegex(ValueError, "reversing the frequency axis"):
            self._convert("restricted")
        # an incomplete set of sparse-node markers (canonical discriminator)
        np.savez(os.path.join(self.d, "incomplete.npz"),
                 **dict(data, ir_freq_n=np.arange(-nf + 1, nf + 1, 2)))
        with self.assertRaisesRegex(ValueError, "incomplete IR-native discriminator"):
            self._convert("incomplete")
        # no metadata: taken as the full centred grid
        nometa = {k: v for k, v in data.items() if k != "freq_index"}
        np.savez(os.path.join(self.d, "nometa.npz"), **nometa)
        self.assertGreater(self._convert("nometa")["asymmetry"], 0.5e-5)

    def test_cli_unreadable_inputs(self):
        """Missing, empty and corrupt archives: exit status 1, a one-line
        error and no traceback, for --check and for the conversion; no
        output file is written."""
        import subprocess
        import sys
        empty = os.path.join(self.d, "empty.npz")
        open(empty, "wb").close()
        corrupt = os.path.join(self.d, "corrupt.npz")
        with open(os.path.join(self.d, "small.npz"), "rb") as fh:
            head = fh.read(200)
        with open(corrupt, "wb") as fh:
            fh.write(head)                       # a truncated zip
        text = os.path.join(self.d, "text.npz")
        with open(text, "w") as fh:
            fh.write("not an archive\n")
        missing = os.path.join(self.d, "missing.npz")
        env = dict(os.environ, PYTHONPATH=os.pathsep.join(
            (os.path.join(os.getcwd(), "src"), os.getcwd())))
        drive = "import sys; from hwave.sigma_split import main; sys.exit(main(sys.argv[1:]))"
        out = os.path.join(self.d, "never.npz")
        for src in (missing, empty, corrupt, text):
            for argv in (["--check", src], [src, out, "--zero-static"]):
                with self.subTest(src=os.path.basename(src), mode=argv[0]):
                    r = subprocess.run([sys.executable, "-B", "-c", drive] + argv, env=env,
                                       capture_output=True, text=True)
                    self.assertEqual(r.returncode, 1, r.stderr)
                    self.assertNotIn("Traceback", r.stderr)
                    self.assertIn("ERROR hwave_sigma_split", r.stderr)
                    self.assertFalse(os.path.exists(out))

    def test_cli(self):
        from hwave.sigma_split import main
        src, out = os.path.join(self.d, "large.npz"), os.path.join(self.d, "cli.npz")
        with _quiet():
            self.assertEqual(main([src, out, "--zero-static"]), 1)
            self.assertFalse(os.path.exists(out))
            self.assertEqual(main([src, out, "--zero-static", "--asymmetry-tol", "-1"]), 1)
        with self.assertLogs("hwave_sigma_split", level="INFO") as cm:
            self.assertEqual(main([os.path.join(self.d, "small.npz"), out, "--zero-static"]), 0)
        self.assertTrue(any("relative Hermitian deviation removed" in m for m in cm.output))
        with _quiet():
            self.assertEqual(main([src, out, "--zero-static", "--asymmetry-tol", "0.1",
                                   "--force"]), 0)


class TestConfigKey(_TmpDir):

    def test_default_and_strict_parsing(self):
        files = _write_dimer(self.d, 1.5)
        s, _ = _dimer_solver(self.d, files)
        self.assertTrue(s.flex_sigma_hermitian_projection)
        s, _ = _dimer_solver(self.d, files, FLEX_Sigma_Hermitian_Projection="off")
        self.assertFalse(s.flex_sigma_hermitian_projection)
        with self.assertRaisesRegex(ValueError, "flex_sigma_hermitian_projection"):
            _dimer_solver(self.d, files, flex_sigma_hermitian_projection="ture")


if __name__ == "__main__":
    unittest.main()
