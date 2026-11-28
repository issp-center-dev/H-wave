"""The paramagnetic Hartree-Fock reference of ``hwave.initial`` (the para
step of the hints) and the facade functions it uses: ``run_para``,
``Bands``/``bands_of``, ``fold_hamiltonian_k`` and ``hamiltonian_k.from_k``.

The para runs are real UHFk runs of a few iterations on small lattices.
They start from a spin-symmetric (or sector-uniform) density, so the tests
never rely on a symmetric start breaking a symmetry: a converged state
stays in the symmetric subspace, up to rounding.
"""
import dataclasses
import logging
import os
import shutil
import tempfile
import unittest
from unittest import mock

import numpy as np

from tests import initial_fixtures


def _square_v(dest, L=6, *, U=2.0, V=0.5, T=0.1, subshape=None, param=None):
    """Square lattice with on-site U and a nearest-neighbour CoulombInter V
    (the Fock term renormalises the bonds)."""
    p = {"Ncond": L * L, "T": T, "IterationMax": 500, "Mix": 0.5}
    if subshape is not None:
        p["SubShape"] = list(subshape)
    p.update(param or {})
    v = {}
    for r in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)):
        v[(r, 0, 0)] = V
    inter = {"CoulombIntra": {((0, 0, 0), 0, 0): U}, "CoulombInter": v}
    return initial_fixtures.write_input(dest, cellshape=(L, L, 1),
                                        transfer=initial_fixtures.square_transfer(),
                                        interactions=inter, param=p)


def _zeeman(dest, L=4, *, h=0.2, U=2.0, T=0.1):
    """Square lattice with a 9-column Transfer carrying an on-site Zeeman
    field ``-h sigma^z`` (commutes with S_z only)."""
    entries = initial_fixtures.spin_independent_entries(initial_fixtures.square_transfer())
    entries[((0, 0, 0), 0, 0, 0, 0)] = -h
    entries[((0, 0, 0), 0, 0, 1, 1)] = h
    text = initial_fixtures.spin_transfer_text(1, entries)
    return initial_fixtures.write_input(
        dest, cellshape=(L, L, 1), transfer={}, interactions={
            "CoulombIntra": {((0, 0, 0), 0, 0): U}},
        param={"Ncond": L * L, "T": T, "IterationMax": 500},
        extra_files={"transfer.dat": text})


def _metastable(dest, *, iteration_max=500):
    """Two orbitals per site, a crystal field 0.5, U_0 = 4 on orbital 0 and
    an on-site inter-orbital V = 3, one electron per site: the free density
    fills orbital 0 and stays there (a metastable orbital-ordered para),
    the state with the electron in orbital 1 is lower by 0.5 per site."""
    return initial_fixtures.two_orbital(dest, 4, delta=0.5, t=0.1, U=(4.0, 0.0), V=3.0,
                                        T=0.02, param={"IterationMax": iteration_max})


def _so_pairhop(dest, ph):
    """The spin-orbital chain of :func:`initial_fixtures.so_chain` (TR
    invariant H0, SOC 0.3) with a bond PairHop ``ph`` along x (Fock on);
    a complex ``ph`` breaks time reversal of the Hartree-Fock map only."""
    transfer = {}
    for s in range(2):
        transfer[((1, 0, 0), s, s)] = -1.0 + 0.3j * (1 - 2 * s)
        transfer[((-1, 0, 0), s, s)] = -1.0 - 0.3j * (1 - 2 * s)
    inter = {"CoulombIntra": {((0, 0, 0), 0, 0): 2.0},
             "PairHop": {((1, 0, 0), 0, 0): ph, ((-1, 0, 0), 0, 0): np.conj(ph)}}
    return initial_fixtures.write_input(
        dest, cellshape=(6, 1, 1), norb=2, transfer=transfer, interactions=inter,
        centers=[(0.0, 0.0, 0.0)] * 2, mode={"enable_spin_orbital": True},
        param={"Ncond": 4, "T": 0.1, "IterationMax": 1000})


def _neel():
    from hwave.initial import Candidate
    from hwave.initial.patterns import Entry, PatternSpec
    spec = PatternSpec(q=("1/2", "1/2", "0"), subshape=(2, 2, 1), phase=0.0, sign=1,
                       modulation="cos", entries=(Entry(0, 0, (0, 0, 0), cz=1.0),))
    return Candidate.user(pattern=spec, label="Neel")


class _Tmp(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self._n = 0
        # the solver logs its iterations at INFO
        logging.getLogger("qlms").setLevel(logging.WARNING)
        self.addCleanup(logging.getLogger("qlms").setLevel, logging.NOTSET)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def dest(self):
        self._n += 1
        return os.path.join(self.tmp, "case{}".format(self._n))

    def facts(self, toml, subshape=None):
        from hwave.initial.facts import model_facts
        return model_facts(toml, subshape=subshape)

    def para(self, toml, **kw):
        from hwave.initial.hints import para_hf, t_hint_bare
        f = self.facts(toml)
        return f, para_hf(f, t_hint_bare(f.problem), **kw)


def _hermitian_green(p, seed):
    from hwave.initial.symmetry import hermitian_green
    rng = np.random.default_rng(seed)
    shape = (p.nvol, p.ns, p.norb, p.ns, p.norb)
    return hermitian_green(rng.normal(size=shape) + 1j * rng.normal(size=shape), p.shape)


# --------------------------------------------------------------------------
# facade: Bands, from_k, the fold
# --------------------------------------------------------------------------

class TestBands(_Tmp):

    def _ham(self, p, seed=1):
        from hwave.solver import initial_facade as fa
        g = 0.05 * _hermitian_green(p, seed)
        return np.asarray(p.ham_trans_k) + fa.hf_potential(
            g, fa.hf_tables(p), include_fock=p.iflag_fock, blocks=p.blocks)

    def test_bands_of_matches_green_function(self):
        from hwave.solver import initial_facade as fa
        cases = {
            "sz free": initial_fixtures.square(self.dest(), 4, U=2.0, T=0.1,
                                               subshape=(1, 1, 1)),
            "2Sz (two mu)": initial_fixtures.square(self.dest(), 4, U=2.0, T=0.1, two_sz=2,
                                                    subshape=(1, 1, 1)),
            "T = 0": initial_fixtures.square(self.dest(), 4, U=2.0, T=0.0, ncond=10,
                                             subshape=(1, 1, 1)),
            "sublattice 2x2": initial_fixtures.square(self.dest(), 4, U=2.0, T=0.1,
                                                      subshape=(2, 2, 1)),
            "spin-orbital": initial_fixtures.so_chain(self.dest(), 6, U=2.0, T=0.1, soc=0.3,
                                                      param={"SubShape": [1, 1, 1]}),
        }
        for name, toml in cases.items():
            with self.subTest(name):
                p = self.facts(toml).problem
                H = self._ham(p)
                b = fa.bands_of(H, p.blocks, p.T, shape=p.shape, ns=p.ns, norb=p.norb,
                                ene_cutoff=p.ene_cutoff)
                G, info = fa.green_function(H, p.T, p.blocks, shape=p.shape, ns=p.ns,
                                            norb=p.norb, ene_cutoff=p.ene_cutoff,
                                            return_info=True)
                gk = np.einsum("kal,kl,kbl->kab", np.conj(b.vec_k), b.occ_k, b.vec_k)
                gr = np.fft.fftn(gk.reshape(p.shape + (p.nd, p.nd)), axes=(0, 1, 2),
                                 norm="forward").reshape(G.shape)
                np.testing.assert_allclose(gr, G, rtol=0, atol=1e-14)
                self.assertEqual(b.mu, info.mu)
                self.assertEqual(len(b.mu), len(p.blocks.group_nconds))
                self.assertEqual(b.kgrid, p.shape)
                self.assertEqual(b.eps_k.shape, (p.nvol, p.nd))
                self.assertEqual(b.vec_k.shape, (p.nvol, p.nd, p.nd))
                # the eigenpairs: H v = eps v
                np.testing.assert_allclose(H @ b.vec_k, b.vec_k * b.eps_k[:, None, :],
                                           atol=1e-12)
                for g, n in enumerate(p.blocks.group_nconds):
                    self.assertAlmostEqual(float(b.occ_k[b.band_group == g].sum()), n,
                                           places=8)
                self.assertFalse(b.occ_k.flags.writeable)
                if name == "2Sz (two mu)":
                    self.assertEqual(len(b.mu), 2)
                    self.assertNotEqual(b.mu[0], b.mu[1])

    def test_fill_green_is_occupations_then_assembly(self):
        from hwave.solver import uhfk_fill
        p = self.facts(initial_fixtures.square(self.dest(), 4, U=2.0, T=0.1, two_sz=2,
                                               subshape=(1, 1, 1))).problem
        ws, vs = uhfk_fill.diag_blocks(self._ham(p), p.blocks.block_info)
        args = (ws, vs, p.blocks.block_info, p.blocks.block_to_group, p.blocks.group_nconds)
        G, mus = uhfk_fill.fill_green(*args, shape=p.shape, ns=p.ns, norb=p.norb, T=p.T,
                                      ene_cutoff=p.ene_cutoff)
        dists, mus2 = uhfk_fill.fill_occupations(*args, shape=p.shape, T=p.T,
                                                 ene_cutoff=p.ene_cutoff)
        np.testing.assert_array_equal(mus, mus2)
        np.testing.assert_array_equal(
            G, uhfk_fill.assemble_green(vs, dists, p.blocks.block_info, shape=p.shape,
                                        ns=p.ns, norb=p.norb))


class TestFold(_Tmp):

    def test_from_k_inverts_to_k(self):
        from hwave.solver import hamiltonian_k
        rng = np.random.default_rng(5)
        x = rng.normal(size=(3, 4, 2, 2, 3)) + 1j * rng.normal(size=(3, 4, 2, 2, 3))
        np.testing.assert_allclose(hamiltonian_k.from_k(hamiltonian_k.to_k(x), (3, 4, 2)),
                                   x, atol=1e-14)
        k = hamiltonian_k.to_k(x)
        np.testing.assert_allclose(hamiltonian_k.to_k(hamiltonian_k.from_k(k, (3, 4, 2))),
                                   k, atol=1e-14)
        with self.assertRaises(ValueError):
            hamiltonian_k.from_k(k, (3, 4, 1))

    def test_fold_bare_transfer_like_solver(self):
        """The fold of the primitive T(k) is the T(k) the solver builds on
        that SubShape (its Transfer folded by reshape_interaction)."""
        from hwave.solver import initial_facade as fa
        nnn = dict(initial_fixtures.square_transfer())
        for r in ((1, 1, 0), (-1, -1, 0)):
            nnn[(r, 0, 0)] = -0.3
        for r in ((1, -1, 0), (-1, 1, 0)):
            nnn[(r, 0, 0)] = 0.2 + 0.1j if r[0] > 0 else 0.2 - 0.1j
        cases = {
            "square 2x2": (initial_fixtures.write_input(
                self.dest(), cellshape=(4, 4, 1), transfer=nnn,
                interactions={"CoulombIntra": {((0, 0, 0), 0, 0): 1.0}},
                param={"Ncond": 16}), (2, 2, 1)),
            "two orbitals 2x1": (initial_fixtures.two_orbital(
                self.dest(), 4, delta=0.5, U=(1.0, 1.0)), (2, 1, 1)),
            "spin-orbital 3x1": (initial_fixtures.so_chain(
                self.dest(), 6, U=1.0, soc=0.3), (3, 1, 1)),
        }
        for name, (toml, B) in cases.items():
            with self.subTest(name):
                pB = self.facts(toml, subshape=B).problem
                folded = fa.fold_hamiltonian_k(fa.primitive_h0_k(pB),
                                               cellshape=pB.cellshape, subshape=B,
                                               norb_orig=pB.norb_orig, ns=pB.ns)
                np.testing.assert_allclose(folded, pB.ham_trans_k, rtol=0, atol=1e-12)

    def test_fold_para_contract(self):
        """fold_para(H_para_k) on SubShape 2x2 equals the Hamiltonian a UHFk
        loaded with that SubShape builds (_make_ham) from the folded G_para."""
        from hwave.initial.hints import fold_para
        from hwave.solver import initial_facade as fa
        toml = _square_v(self.dest(), subshape=(2, 2, 1))
        f, para = self.para(toml)
        self.assertEqual(para.status, "converged", para.skip_reason)
        B = (2, 2, 1)
        fB = self.facts(toml, subshape=B)
        solver = fa._build_solver(toml, subshape=B)[0]
        solver.Green = fa.fold_layout(para.G_para, cellshape=fB.cellshape, subshape=B,
                                      norb_orig=fB.problem.norb_orig, ns=fB.problem.ns)
        solver._make_ham()
        mine = fold_para(para.H_para_k, fB, B)
        self.assertEqual(mine.shape, solver.ham.shape)
        np.testing.assert_allclose(mine, solver.ham, rtol=0, atol=1e-12)
        # the Fock bond renormalisation is in it: not the bare T(k) plus a constant
        diff = mine - np.asarray(fB.problem.ham_trans_k)
        self.assertGreater(float(np.ptp(np.abs(diff))), 1e-3)

    def test_fold_refuses_wrong_shape(self):
        from hwave.solver import initial_facade as fa
        with self.assertRaises(ValueError):
            fa.fold_hamiltonian_k(np.zeros((8, 2, 2)), cellshape=(4, 4, 1),
                                  subshape=(2, 2, 1), norb_orig=1, ns=2)


# --------------------------------------------------------------------------
# the para step
# --------------------------------------------------------------------------

class TestRunPara(_Tmp):

    def test_input_eps_above_the_minimum_is_kept(self):
        from hwave.solver import initial_facade as fa
        toml = _square_v(self.dest(), param={"EPS": 13})
        self.assertEqual(fa._input_eps(toml), 13)
        p1 = self.facts(toml, subshape=(1, 1, 1)).problem
        start = fa.green_function(p1.ham_trans_k, 0.1, p1.blocks, shape=p1.shape,
                                  ns=p1.ns, norb=p1.norb)
        with mock.patch.object(fa, "_build_solver", wraps=fa._build_solver) as bs:
            fa.run_para(toml, T=0.1, start_G=start, work_dir=os.path.join(self.tmp, "r"))
        self.assertEqual(bs.call_args.kwargs["param"]["EPS"], 13)

    def test_run_para_drives_the_solver_on_the_primitive_cell(self):
        from hwave.solver import initial_facade as fa
        toml = _square_v(self.dest(), subshape=(2, 2, 1))
        p1 = self.facts(toml, subshape=(1, 1, 1)).problem
        start = fa.green_function(p1.ham_trans_k, 0.1, p1.blocks, shape=p1.shape,
                                  ns=p1.ns, norb=p1.norb)
        wd = os.path.join(self.tmp, "run")
        run = fa.run_para(toml, T=0.1, start_G=start, work_dir=wd)
        self.assertTrue(run.converged)
        self.assertGreater(run.iterations, 0)
        self.assertEqual(run.eps, 12)                     # the input has EPS = 8
        self.assertLess(run.residual, 1e-12)
        self.assertEqual(run.G.shape, (36, 2, 1, 2, 1))
        for name in ("start.npz", "green.npz", "energy.dat", "scf_status.dat"):
            self.assertTrue(os.path.isfile(os.path.join(wd, name)), name)
        with np.load(os.path.join(wd, "green.npz")) as data:
            np.testing.assert_array_equal(data["green"], run.G)
        with open(os.path.join(wd, "energy.dat")) as fh:
            line = [ln for ln in fh if ln.startswith("FreeEnergy_Total")][0]
        self.assertAlmostEqual(float(line.split("=")[1]), run.free_energy, places=10)
        np.testing.assert_allclose(
            run.H_k, np.asarray(p1.ham_trans_k) + fa.hf_potential(
                run.G, fa.hf_tables(p1), include_fock=p1.iflag_fock, blocks=p1.blocks),
            atol=1e-14)
        with self.assertRaises(ValueError):
            fa.run_para(toml, T=0.1, start_G=np.zeros((4, 2, 4, 2, 4)), work_dir=wd)


class TestParaHF(_Tmp):

    def test_para_preserves_symmetry(self):
        from hwave.solver import initial_facade as fa
        cases = {
            "square, U and V": (_square_v(self.dest()), ("S_x", "S_y", "S_z")),
            "spin-orbital, TR invariant": (initial_fixtures.so_chain(
                self.dest(), 6, U=2.0, T=0.1, soc=0.3, param={"IterationMax": 500}),
                ("time reversal",)),
        }
        for name, (toml, checked) in cases.items():
            with self.subTest(name):
                f, para = self.para(toml)
                self.assertEqual(para.status, "converged", para.skip_reason)
                self.assertIsNone(para.skip_reason)
                self.assertLessEqual(para.residual, 1e-8)
                self.assertEqual(para.checked, checked)
                p = self.facts(toml, subshape=(1, 1, 1)).problem
                self.assertAlmostEqual(para.T, f.problem.T)       # T > 0: T_hint = T
                H = fa.primitive_h0_k(p) + fa.hf_potential(
                    para.G_para, fa.hf_tables(p), include_fock=p.iflag_fock,
                    blocks=p.blocks)
                np.testing.assert_allclose(para.H_para_k, H, atol=1e-14)
                # the bands are those of H_para at T_hint, and self-consistent
                b = para.bands
                G = fa.green_function(para.H_para_k, para.T, p.blocks, shape=p.shape,
                                      ns=p.ns, norb=p.norb, ene_cutoff=p.ene_cutoff)
                np.testing.assert_allclose(G, para.G_para, atol=1e-6)
                self.assertEqual(b.T, para.T)
                self.assertAlmostEqual(para.free_energy_per_cell * 1.0,
                                       para.free_energy / int(np.prod(f.cellshape)))

    def test_2sz_para_polarised_recorded(self):
        from hwave.initial import report
        from hwave.initial.hints import hint_record, t_hint_bare
        toml = _square_v(self.dest(), param={"2Sz": 2})
        f, para = self.para(toml)
        self.assertGreater(para.iterations, 1)
        self.assertEqual(para.status, "converged", para.skip_reason)
        self.assertEqual(para.checked, ("transverse spin", "sector counts"))
        self.assertLessEqual(para.residual, 1e-8)
        counts = para.polarisation["sector_counts"]
        np.testing.assert_allclose(counts, [19.0, 17.0], atol=1e-8)
        # S_z per site (N_up - N_down) / (2 N_cell), recorded, not a residual
        self.assertAlmostEqual(para.polarisation["spin"][0][2], 2 / 72, places=8)
        self.assertEqual(list(para.polarisation["allowed"]), ["sectors"])
        rec = hint_record(t_hint_bare(f.problem), para, vertex={
            "gamma_self_adjoint": True, "method": "dense", "q_status": {}, "q_set": []})
        self.assertEqual(rec["para_status"], "converged")
        self.assertEqual(rec["para_polarisation"]["sector_counts"], list(counts))
        errs = []
        report._check_hint(rec, errs)
        self.assertEqual(errs, [])

    def test_zeeman_transfer_drops_spin_residual(self):
        from hwave.initial import symmetry
        toml = _zeeman(self.dest())
        f, para = self.para(toml)
        h0 = symmetry.h0_spin_commutators(f.problem)
        self.assertGreater(h0[0], 1e-3)
        self.assertGreater(h0[1], 1e-3)
        self.assertLessEqual(h0[2], 1e-14)
        self.assertEqual(para.status, "converged", para.skip_reason)
        # S_z is H0-induced: recorded, not checked; S_x, S_y still checked
        self.assertEqual(para.checked, ("S_x", "S_y"))
        self.assertEqual(list(para.polarisation["allowed"]), ["z"])
        self.assertGreater(abs(para.polarisation["spin"][0][2]), 1e-3)
        self.assertLessEqual(para.residual, 1e-8)

    def test_metastable_orbital_para_caught_by_check(self):
        toml = _metastable(self.dest())
        f, plain = self.para(toml)
        self.assertEqual(plain.status, "converged", plain.skip_reason)
        self.assertIsNone(plain.free_energy_check)
        self.assertEqual(plain.warnings, ())
        # the free start stays in orbital 0 (the metastable branch)
        g0 = np.real(np.diagonal(plain.G_para[0].reshape(4, 4)))
        self.assertGreater(g0[0] + g0[2], 0.9)
        with self.assertLogs("qlms.initial.hints", level="WARNING") as cm:
            f, checked = self.para(toml, para_check=True)
        self.assertEqual(checked.status, "converged")
        self.assertEqual(len(checked.warnings), 1)
        self.assertIn("para check", checked.warnings[0])
        self.assertIn("metastable", cm.output[0])
        self.assertLess(checked.free_energy_per_cell, plain.free_energy_per_cell - 0.1)
        self.assertEqual(checked.free_energy_check, checked.free_energy_per_cell)
        g1 = np.real(np.diagonal(checked.G_para[0].reshape(4, 4)))
        self.assertGreater(g1[1] + g1[3], 0.9)

    def test_para_check_agrees_without_warning(self):
        """An input with EPS = 8 (whose two starts differed by ~7e-8 per cell
        at that tolerance) or the default 6 gives no para-check warning: the
        para runs use EPS = max(EPS_input, 12)."""
        from hwave.initial.hints import hint_record, t_hint_bare
        for eps, want in ((8, 12), (6, 12)):
            with self.subTest(EPS=eps):
                toml = _square_v(self.dest(), param={"EPS": eps})
                with self.assertNoLogs("qlms.initial.hints", level="WARNING"):
                    f, para = self.para(toml, para_check=True)
                self.assertEqual(para.status, "converged")
                self.assertEqual(para.warnings, ())
                self.assertEqual(para.eps, want)
                self.assertLessEqual(abs(para.free_energy_check - para.free_energy_per_cell),
                                     1e-8)
                self.assertEqual(hint_record(t_hint_bare(f.problem), para,
                                             skip_reason="x")["para_eps"], None)

    def test_failure_falls_back_to_bare(self):
        from hwave.initial import report
        from hwave.initial.hints import hint_record, t_hint_bare
        from hwave.initial.writer import initial_state
        # 6 x 6: on 4 x 4 the T_hint occupations are pinned by the level gaps
        # and the first step already converges
        toml = _square_v(self.dest(), param={"IterationMax": 1})
        f, para = self.para(toml)
        self.assertEqual(para.status, "not_converged")
        self.assertIn("did not converge", para.skip_reason)
        rec = hint_record(t_hint_bare(f.problem), para)
        self.assertEqual(rec["status"], "skipped")
        self.assertEqual(rec["skip_reason"], para.skip_reason)
        self.assertEqual(rec["para_status"], "not_converged")
        errs = []
        report._check_hint(rec, errs)
        self.assertEqual(errs, [])
        cand = _neel()
        out = os.path.join(os.path.dirname(toml), "initial.npz")
        initial_state(toml, cand, out=out, para=para)
        with np.load(out, allow_pickle=False) as data:
            self.assertEqual(str(data["seed_base"][()]), '"bare"')

    def test_no_interaction_is_skipped(self):
        toml = initial_fixtures.square(self.dest(), 4, U=0.0, T=0.1)
        f, para = self.para(toml)
        self.assertEqual(para.status, "skipped")
        self.assertIn("n_c = 0", para.skip_reason)
        self.assertIsNone(para.G_para)

    def test_symmetry_breaking_state_is_not_accepted(self):
        """A converged state with a spin density where the model keeps the
        spin symmetry fails the residual (the reference is not a para)."""
        from hwave.initial.hints import para_from_green, t_hint_bare
        toml = _square_v(self.dest())
        f, para = self.para(toml)
        G = np.array(para.G_para)
        G[0, 0, 0, 0, 0] += 0.01
        G[0, 1, 0, 1, 0] -= 0.01
        res = para_from_green(f, t_hint_bare(f.problem), G, free_energy=para.free_energy,
                              converged=True)
        self.assertEqual(res.status, "not_converged")
        self.assertIn("breaks a symmetry", res.skip_reason)
        self.assertAlmostEqual(res.residual, 0.01, places=10)

    def test_para_from_green_rebuilds_the_result(self):
        from hwave.initial.hints import para_from_green, t_hint_bare
        toml = _square_v(self.dest(), subshape=(2, 2, 1))
        f, para = self.para(toml)
        again = para_from_green(f, t_hint_bare(f.problem), para.G_para,
                                free_energy=para.free_energy, converged=True,
                                iterations=para.iterations)
        np.testing.assert_array_equal(again.H_para_k, para.H_para_k)
        np.testing.assert_array_equal(again.bands.eps_k, para.bands.eps_k)
        self.assertEqual(again.status, para.status)
        self.assertEqual(again.residual, para.residual)

    def test_work_dir_keeps_the_para_green(self):
        import json
        toml = _metastable(self.dest())
        wd = os.path.join(self.tmp, "WORK", "para_hf")
        f, para = self.para(toml, para_check=True, work_dir=wd)
        with np.load(os.path.join(wd, "green.npz")) as data:
            np.testing.assert_array_equal(data["green"], para.G_para)
        with open(os.path.join(wd, "para_hf.json")) as fh:
            rec = json.load(fh)
        self.assertEqual(rec["status"], "converged")
        self.assertEqual(rec["selected_start"], "check")
        self.assertEqual(rec["eps"], 12)
        self.assertAlmostEqual(rec["free_energy_per_cell"], para.free_energy_per_cell)
        self.assertEqual(rec["input_sha256"], dict(f.input_sha256))
        self.assertTrue(os.path.isdir(os.path.join(wd, "main")))
        self.assertTrue(os.path.isdir(os.path.join(wd, "check")))

    def test_hint_record_keys(self):
        from hwave.initial import report
        from hwave.initial.hints import hint_record, t_hint_bare
        f, para = self.para(_square_v(self.dest()))
        th = t_hint_bare(f.problem)
        vertex = {"gamma_self_adjoint": True, "method": "dense", "q_status": {"0,0,0": "ok"},
                  "q_set": ["0,0,0"]}
        for rec in (hint_record(th, para, vertex=vertex),
                    hint_record(th, para, skip_reason="large tier"),
                    hint_record(th, skip_reason="large tier")):
            self.assertEqual(set(rec), set(report._HINT_KEYS))
            errs = []
            report._check_hint(rec, errs)
            self.assertEqual(errs, [])
        rec = hint_record(th, para, vertex=vertex)
        self.assertEqual(rec["status"], "computed")
        self.assertEqual(rec["q_set"], ["0,0,0"])
        self.assertEqual(rec["para_eps"], 12)
        # a computed record never goes out without the vertex fields
        with self.assertRaises(ValueError):
            hint_record(th, para)
        bad = dict(rec, method=None)
        errs = []
        report._check_hint(bad, errs)
        self.assertTrue(any("hint.method" in e for e in errs), errs)
        self.assertEqual(rec["para_free_energy"], para.free_energy_per_cell)
        self.assertEqual(rec["para_residual"], para.residual)
        self.assertEqual(rec["para_polarisation"]["allowed"], [])
        with self.assertRaises(ValueError):
            hint_record(th)


class TestWriterPara(_Tmp):

    def _neel(self):
        return _neel()

    def test_initial_state_reuses_para(self):
        from hwave.initial import hints
        from hwave.initial.writer import initial_state
        from hwave.solver import initial_facade
        toml = _square_v(self.dest(), subshape=(2, 2, 1))
        f, para = self.para(toml)
        self.assertEqual(para.status, "converged")
        out = os.path.join(os.path.dirname(toml), "initial.npz")
        with mock.patch.object(initial_facade, "run_para",
                               wraps=initial_facade.run_para) as rp, \
                mock.patch.object(hints, "para_hf", wraps=hints.para_hf) as ph:
            initial_state(toml, self._neel(), out=out, para=para)
        self.assertEqual(rp.call_count, 0)
        self.assertEqual(ph.call_count, 0)
        with np.load(out, allow_pickle=False) as data:
            self.assertEqual(str(data["seed_base"][()]), '"para_hf"')
            g_para = np.array(data["green_sublattice"])
        bare = os.path.join(os.path.dirname(toml), "bare.npz")
        initial_state(toml, self._neel(), out=bare)
        with np.load(bare, allow_pickle=False) as data:
            self.assertEqual(str(data["seed_base"][()]), '"bare"')
            g_bare = np.array(data["green_sublattice"])
        # the Fock-renormalised base gives another seed
        self.assertGreater(float(np.max(np.abs(g_para - g_bare))), 1e-6)

    def test_para_of_another_model_refused(self):
        from hwave.initial.writer import initial_state
        f, para = self.para(_square_v(self.dest(), V=0.5))
        other = _square_v(self.dest(), V=0.7, subshape=(2, 2, 1))
        with self.assertRaises(ValueError):
            initial_state(other, self._neel(), out=os.path.join(self.tmp, "x.npz"),
                          para=para)
        # same input files, another Ncond in the TOML
        same_files = _square_v(self.dest(), V=0.5, subshape=(2, 2, 1),
                               param={"Ncond": 34})
        with self.assertRaisesRegex(ValueError, "ncond"):
            initial_state(same_files, self._neel(), out=os.path.join(self.tmp, "y.npz"),
                          para=para)
        with self.assertRaises(TypeError):
            initial_state(other, self._neel(), out=os.path.join(self.tmp, "x.npz"),
                          para="para")


class TestMapSymmetry(_Tmp):
    """The map leg of the residual rule: a symmetry enters only when H0 AND
    the Hartree-Fock map keep it."""

    def test_complex_pairhop_takes_tr_out_of_the_residual(self):
        from hwave.initial import symmetry
        cases = {"real": 0.4, "complex": 0.4 + 0.3j}
        out = {}
        for name, ph in cases.items():
            f, para = self.para(_so_pairhop(self.dest(), ph))
            tr = symmetry.tr_check(f.problem)
            self.assertLessEqual(tr.h0_residual, 1e-14)      # H0 is TR invariant
            out[name] = (tr, para)
        tr, para = out["real"]
        self.assertTrue(tr.invariant)
        self.assertEqual(para.status, "converged", para.skip_reason)
        self.assertEqual(para.checked, ("time reversal",))
        self.assertLessEqual(para.residual, 1e-8)
        tr, para = out["complex"]
        self.assertFalse(tr.invariant)
        self.assertGreater(tr.hf_residual, 1e-3)              # the map breaks TR
        self.assertEqual(para.status, "converged", para.skip_reason)
        self.assertEqual(para.checked, ())
        self.assertNotIn("time reversal", para.residual_parts)
        self.assertGreater(para.polarisation["time_reversal_odd"], 1e-4)
        self.assertEqual(len(para.polarisation["allowed"]), 1)

    def test_axis_projection(self):
        from hwave.initial import symmetry
        p = self.facts(initial_fixtures.square(self.dest(), 4, U=1.0,
                                               subshape=(1, 1, 1))).problem
        g = _hermitian_green(p, 3)
        for mu in range(3):
            with self.subTest(axis=mu):
                P = symmetry.axis_projection(g, mu, p.norb)
                np.testing.assert_allclose(symmetry.axis_projection(P, mu, p.norb), P,
                                           atol=1e-15)                    # idempotent
                Q = g - P
                np.testing.assert_allclose(symmetry.axis_projection(Q, mu, p.norb), 0,
                                           atol=1e-15)                    # complement
                # P is invariant under the rotations about mu, Q is not
                axis = np.eye(3)[mu]
                u = symmetry.spin_rotation(axis, 0.7, p.norb)
                np.testing.assert_allclose(symmetry.rotate_green(P, u), P, atol=1e-13)
                self.assertGreater(float(np.max(np.abs(symmetry.rotate_green(Q, u) - Q))),
                                   1e-3)
        # a spin-free G is in every axis subspace
        spin_free = symmetry.axis_projection(symmetry.axis_projection(g, 0, p.norb), 2,
                                             p.norb)
        for mu in range(3):
            np.testing.assert_allclose(symmetry.axis_projection(spin_free, mu, p.norb),
                                       spin_free, atol=1e-15)

    def test_hf_axis_residuals(self):
        from hwave.initial import symmetry
        p = self.facts(initial_fixtures.square(self.dest(), 4, U=1.0,
                                               subshape=(1, 1, 1))).problem
        self.assertEqual(symmetry.hf_axis_residuals(p), (0.0, 0.0, 0.0))
        # a map that rotates every output about z keeps only the z subspace
        u = symmetry.spin_rotation((0.0, 0.0, 1.0), 0.9, p.norb)
        nd = p.nd

        def rotating(problem):
            def sigma(green):
                m = np.asarray(green).reshape(problem.nvol, nd, nd)
                return u @ m @ u.conj().T
            return sigma
        with mock.patch.object(symmetry, "_sigma_map", rotating):
            res = symmetry.hf_axis_residuals(p)
        self.assertGreater(res[0], 1e-2)
        self.assertGreater(res[1], 1e-2)
        self.assertLessEqual(res[2], 1e-14)
        with self.assertRaises(ValueError):
            symmetry.hf_axis_residuals(self.facts(initial_fixtures.so_chain(
                self.dest(), 6, U=1.0)).problem)


class TestParaRecord(_Tmp):

    def test_check_side_switch_warns(self):
        from hwave.initial import hints
        real = hints._one_run
        calls = []

        def first_fails(facts1, thint, start, work_dir):
            run, res = real(facts1, thint, start, work_dir)
            calls.append(work_dir)
            if len(calls) == 1:
                res = dataclasses.replace(res, status="not_converged",
                                          skip_reason="forced failure")
            return run, res
        toml = _square_v(self.dest(), param={"EPS": 12})
        with mock.patch.object(hints, "_one_run", first_fails), \
                self.assertLogs("qlms.initial.hints", level="WARNING") as cm:
            f, para = self.para(toml, para_check=True)
        self.assertEqual(para.status, "converged")
        self.assertEqual(len(para.warnings), 1)
        self.assertIn("forced failure", para.warnings[0])
        self.assertIn("uniform orbital occupation is kept", cm.output[0])

    def test_record_written_when_skipped(self):
        import json
        from hwave.solver import initial_facade
        wd = os.path.join(self.tmp, "skip0")
        f, para = self.para(initial_fixtures.square(self.dest(), 4, U=0.0, T=0.1),
                            work_dir=wd)
        with open(os.path.join(wd, "para_hf.json")) as fh:
            rec = json.load(fh)
        self.assertEqual(rec["status"], "skipped")
        self.assertIsNone(rec["selected_start"])
        self.assertFalse(os.path.exists(os.path.join(wd, "green.npz")))
        wd = os.path.join(self.tmp, "skip1")
        with mock.patch.object(initial_facade, "run_para",
                               side_effect=ValueError("refused for the test")):
            f, para = self.para(_square_v(self.dest()), work_dir=wd)
        self.assertEqual(para.status, "skipped")
        with open(os.path.join(wd, "para_hf.json")) as fh:
            rec = json.load(fh)
        self.assertIn("refused for the test", rec["skip_reason"])
        self.assertEqual(rec["input_digest"], __import__(
            "hwave.initial.hints", fromlist=["x"]).input_digest(f.input_sha256))

    def test_non_numeric_eps_is_a_value_error(self):
        from hwave.solver import initial_facade as fa
        toml = _square_v(self.dest())
        with open(toml) as fh:
            text = fh.read()
        with open(toml, "w") as fh:
            fh.write(text.replace("EPS = 8", 'EPS = "8"'))
        with self.assertRaisesRegex(ValueError, "EPS"):
            fa._input_eps(toml)

    def test_not_converged_reason_names_eps_and_iteration_max(self):
        f, para = self.para(_square_v(self.dest(), param={"IterationMax": 3}))
        self.assertEqual(para.status, "not_converged")
        self.assertIn("EPS = 12", para.skip_reason)
        self.assertIn("IterationMax = 3", para.skip_reason)
        self.assertEqual(para.iteration_max, 3)


class TestParaFromWork(_Tmp):
    """A saved para reference is reused only for the model it belongs to:
    the input files AND the TOML settings they do not fix."""

    VARIANTS = {
        "ncond": ("Ncond = 16", "Ncond = 14"),
        "two_sz": ("Ncond = 16", "Ncond = 16\n2Sz = 0"),
        "iflag_fock": ('mode = "UHFk"', 'mode = "UHFk"\nflag_fock = false'),
        "cellshape": ("CellShape = [4, 4, 1]", "CellShape = [6, 4, 1]"),
        "ene_cutoff": ("Ncond = 16", "Ncond = 16\nene_cutoff = 50.0"),
    }

    def variant(self, src, key):
        """A copy of the model directory of ``src`` (the same referenced
        files) with one TOML setting changed."""
        old, new = self.VARIANTS[key]
        dst = os.path.join(self.tmp, "variant-" + key)
        shutil.copytree(os.path.dirname(src), dst)
        toml = os.path.join(dst, os.path.basename(src))
        with open(toml) as fh:
            text = fh.read()
        self.assertEqual(text.count(old), 1, key)
        with open(toml, "w") as fh:
            fh.write(text.replace(old, new))
        return toml

    def test_other_model_refused(self):
        """Ncond, 2Sz, flag_fock, CellShape or ene_cutoff changed with the same input
        files: para_from_work refuses with the differing setting named, for
        a record with a Green function and for a skipped one; the same model
        is reused through candidate_states(para=...)."""
        from hwave.initial import hints
        from hwave.initial.enumerate import candidate_states
        toml = _square_v(self.dest(), L=4, subshape=(1, 1, 1), param={"Ncond": 16})
        wd = os.path.join(self.tmp, "para")
        f, para = self.para(toml, work_dir=wd)
        self.assertEqual(para.status, "converged", para.skip_reason)
        again = hints.para_from_work(f, wd)
        self.assertEqual(again.status, "converged")
        np.testing.assert_array_equal(again.G_para, para.G_para)
        cs = candidate_states(toml, tier="small", n_random=0, para=again)
        self.assertEqual(cs.model["hint"]["status"], "computed")
        skipped = os.path.join(self.tmp, "skipped")
        os.makedirs(skipped)
        hints.write_para_record(hints.para_failed(f, hints.t_hint_bare(f.problem),
                                                  "failed for the test"), None, skipped)
        self.assertIsNone(hints.para_from_work(f, skipped))
        for key in self.VARIANTS:
            with self.subTest(key):
                other = self.variant(toml, key)
                fo = self.facts(other)
                self.assertEqual(fo.input_sha256, f.input_sha256)
                for work in (wd, skipped):
                    with self.assertRaisesRegex(ValueError, "another model: .*" + key):
                        candidate_states(other, tier="small", n_random=0,
                                         para=hints.para_from_work(fo, work))

    def test_in_memory_para_checked(self):
        """A ParaResult passed directly (no work directory): one computed at
        the default ene_cutoff is refused for the same files with
        ene_cutoff = 50, and one without a model is refused outright, by
        candidate_states and by the writer."""
        from hwave.initial.enumerate import candidate_states
        from hwave.initial.writer import initial_state
        toml = _square_v(self.dest(), L=4, subshape=(1, 1, 1), param={"Ncond": 16})
        f, para = self.para(toml)
        self.assertEqual(para.status, "converged", para.skip_reason)
        self.assertEqual(dict(para.model)["ene_cutoff"], 100.0)
        other = self.variant(toml, "ene_cutoff")
        with self.assertRaisesRegex(ValueError, r"another model \(ene_cutoff differ\)"):
            candidate_states(other, tier="small", n_random=0, para=para)
        bare = dataclasses.replace(para, model=None)
        out = os.path.join(self.tmp, "x.npz")
        with self.assertRaisesRegex(ValueError, "names no model"):
            candidate_states(toml, tier="small", n_random=0, para=bare)
        with self.assertRaisesRegex(ValueError, "names no model"):
            initial_state(toml, _neel(), out=out, para=bare)
        self.assertFalse(os.path.exists(out))
        # the reference with its model is accepted for its own input
        initial_state(toml, _neel(), out=out, para=para)
        self.assertTrue(os.path.exists(out))

    def test_other_input_files_refused(self):
        """Other referenced files: refused (the input digest), also for a
        skipped record."""
        from hwave.initial import hints
        toml = _square_v(self.dest(), L=4, param={"Ncond": 16})
        wd = os.path.join(self.tmp, "skipped")
        os.makedirs(wd)
        f = self.facts(toml)
        hints.write_para_record(hints.para_failed(f, hints.t_hint_bare(f.problem), "x"),
                                None, wd)
        other = _square_v(self.dest(), L=4, U=3.0, param={"Ncond": 16})
        with self.assertRaisesRegex(ValueError, "other input files"):
            hints.para_from_work(self.facts(other), wd)


class TestParaProvenance(_Tmp):

    def test_para_T_and_digest_in_provenance(self):
        import json
        from hwave.initial.hints import input_digest
        from hwave.initial.writer import check_initial_state, initial_state
        toml = _square_v(self.dest(), subshape=(2, 2, 1))
        f, para = self.para(toml)
        out = os.path.join(os.path.dirname(toml), "initial.npz")
        initial_state(toml, _neel(), out=out, para=para)
        with np.load(out, allow_pickle=False) as data:
            self.assertEqual(json.loads(str(data["seed_base_para_T"][()])), para.T)
            self.assertEqual(json.loads(str(data["seed_base_para_digest"][()])),
                             input_digest(f.input_sha256))
        self.assertTrue(check_initial_state(toml, out, require_provenance=True).verified)
        bare = os.path.join(os.path.dirname(toml), "bare.npz")
        initial_state(toml, _neel(), out=bare)
        with np.load(bare, allow_pickle=False) as data:
            self.assertIsNone(json.loads(str(data["seed_base_para_T"][()])))
            self.assertIsNone(json.loads(str(data["seed_base_para_digest"][()])))

    def test_para_digest_checked_against_input(self):
        """require_provenance: a para seed's seed_base_para_digest must be a
        lowercase SHA-256 and the digest of the file's own input_sha256;
        either defect is an InitialStateError (not a missing provenance)."""
        import json
        from hwave.initial.hints import input_digest
        from hwave.initial.writer import (InitialStateError, ProvenanceMissingError,
                                          check_initial_state, initial_state)
        from tests.test_initial_writer import _load, _save
        toml = _square_v(self.dest(), subshape=(2, 2, 1))
        f, para = self.para(toml)
        out = os.path.join(os.path.dirname(toml), "initial.npz")
        initial_state(toml, _neel(), out=out, para=para)
        members = _load(out)
        good = input_digest(f.input_sha256)
        self.assertEqual(json.loads(str(members["seed_base_para_digest"][()])), good)
        other = "0" * 64
        self.assertNotEqual(other, good)
        for bad, msg in (("Z" * 64, "SHA-256"), (good.upper(), "SHA-256"),
                         (other, "not the digest of the file's own input_sha256")):
            with self.subTest(bad=bad):
                m = dict(members)
                m["seed_base_para_digest"] = np.array(json.dumps(bad))
                path = os.path.join(self.tmp, "bad.npz")
                _save(path, m)
                self.assertEqual(json.loads(str(_load(path)["seed_base_para_digest"][()])),
                                 bad)
                with self.assertRaises(InitialStateError) as cm:
                    check_initial_state(toml, path, require_provenance=True)
                self.assertNotIsInstance(cm.exception, ProvenanceMissingError)
                self.assertIn(msg, str(cm.exception))
        self.assertTrue(check_initial_state(toml, out, require_provenance=True).verified)

    def test_para_of_another_T_refused(self):
        from hwave.initial.writer import initial_state
        f, para = self.para(_square_v(self.dest(), T=0.1))
        hot = _square_v(self.dest(), T=0.3, subshape=(2, 2, 1))
        out = os.path.join(self.tmp, "x.npz")
        with self.assertRaisesRegex(ValueError, r"T_hint = 0\.1.*T_hint = 0\.3"):
            initial_state(hot, _neel(), out=out, para=para)
        self.assertFalse(os.path.exists(out))
        # the run whose hint was computed at a user T_hint = 0.1 takes it
        initial_state(hot, _neel(), out=out, para=para, t_hint=0.1)
        with np.load(out, allow_pickle=False) as data:
            self.assertEqual(str(data["seed_base"][()]), '"para_hf"')


if __name__ == "__main__":
    unittest.main()
