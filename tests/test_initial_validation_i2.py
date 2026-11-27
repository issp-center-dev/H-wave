"""Validation cases of ``hwave.initial`` for the small tier: V4, V5 and V6.

* **V4** a Wannier-like multi-orbital model
  (:func:`tests.initial_fixtures.wannier_like`, the settings of
  :data:`tests.initial_reference_capture.V4_MODEL`): two interacting orbitals
  with the on-site Kanamori couplings U = 4 (CoulombIntra), U' = U - 2J = 2.8
  (CoulombInter) and J = 0.6 (Hund, Exchange, PairHop), a ligand that does
  not interact, Wannier90-format files, 4 x 4, T = 0.2, Ncond = 40. A complex
  inter-orbital hopping makes the model complex, so the orbital-pair layouts
  below are pinned (with a real model a layout and its pair-swapped reading
  would agree). Checked against independent code paths of the package:

  - ``n_c = 2`` (the ligand is not interacting);
  - ``static_lindhard`` on the para Hartree-Fock bands against the label-form
    oracle ``tests/so_bond_oracles.lindhard_static`` (1e-10) and against the
    static block of the RPA ``chi0q`` on the same bands (``Nmat = 8192``,
    ``coeff_tail = 0``) to 1e-4, which covers its leading Matsubara
    truncation offset ``1/(pi^2 T Nmat) = 6.2e-5``; every on-site orbital
    pair, ligand included, at every q;
  - the leading eigenvalue at each q equals the Stoner factor of
    ``flex_diagnostics.static_instability`` fed the same chi0 on the
    spin-orbital RPA vertex: the antisymmetrised vertex (direct plus
    exchange crossing) with the Fock term, the direct (ring) vertex without
    it (1e-10);
  - the spin and charge eigenvalues equal those of the multi-orbital S and C
    matrices of ``hwave_sc`` (Kuroki et al.; ``chi0 S`` and ``-chi0 C``) on a
    density-only subcase (CoulombIntra and CoulombInter) and on the full
    Kanamori model (1e-10);
  - the reduced charge and S^z sectors hold the nonzero spectrum of the full
    pair space (S^z three times; 1e-12);
  - seed shape: the density change of the field seed, ``(rho(h) - rho(0))/h``
    at ``h = 1e-6`` along the conjugate field ``O = Gamma x`` of the leading
    pair, equals ``chi0 Gamma x`` on every on-site slot, ligand included
    (1e-5; the field enters as ``H - h O``);
  - multi-start (heavy): ``prepare`` in the small tier on SubShape 1 x 1 and
    2 x 2, every run, ``collect``, against ``tests/initial_ref/v4.json``;
  - the active-window span ``W_act`` near a band touching the window edge
    jumps (0.4 -> 4.6 at delta = -2.4, an accepted property of energy-rank
    bands), deterministically and as recorded in ``v4.json``.

  Multi-start comparison: statuses and free energies per site (1e-6) of every
  spin-polarised row; final labels of the rows that end below the para row
  by more than 1e-3 per site (an ordered state), random rows excepted. Rows
  seeded without a spin component (para, charge, orbital) are not compared:
  their seed keeps the spin symmetry, so whether a run leaves it depends on
  the round-off of the numpy and BLAS build. Rows that reach the same state
  differ by up to 3e-8 per site (EPS = 10); which of them is the lowest is
  therefore not compared, only its final label and free energy.

* **V5** spin-orbital mode (:func:`tests.initial_fixtures.so_two_orbital`,
  :data:`tests.initial_reference_capture.V5_MODEL`): two orbitals with the
  inter-orbital spin-orbit coupling ``i lambda sigma^z`` only, U on both,
  4 x 4 at half filling, T = 0.2. The runs from the Neel seeds with the moment
  along x, y and z: their free energies against ``tests/initial_ref/v5.json``
  (1e-8); E_x = E_y (sigma^z keeps the rotations about z; 1e-8); invariance
  under lambda -> -lambda (the gauge b -> -b of orbital 1; 1e-8); the three
  agree at lambda = 0 (1e-8). Which axis is lower is the frozen value, not a
  claim of the test. The small tier on the same model: time reversal with
  the spin-orbit coupling passes the check, the labels carry no gauge
  warning (the coupling is a spin part), the leading hint is a time-reversal
  odd spin pattern, and a UHFr target is refused.

* **V6** the Fock channels.

  - CoulombInter only (square 8 x 8, V = 1.2): the small tier lists bond-order
    instability rows with the Fock term and none without it; a bond-order
    seed at q = (1/2, 1/2) written by ``prepare`` survives the write -> run ->
    read round trip (``IterationMax = 1``, ``Mix = 0``; 1e-12).
  - PairHop (complex, on-site and on the x bond), with and without the Fock
    term, at the commensurate q = 1/3 and 2/3 along x (3 x 3, SubShape
    3 x 1): the finite-difference Jacobian J of the whole SCF map ``rho ->
    green_function(H0 + hf_potential(rho^T))^T`` at the para state (T = 0.1,
    central differences, step 1e-5, on every slot of the folded cell)
    against ``K(q) = -chi0(q) Gamma(q)`` on the pair space. J holds every q of
    the cell; at q = 0 the fixed electron number removes the uniform density
    change, so there chi0 is replaced by its projection ``chi0 - chi0 n n^T
    chi0 / (n^T chi0 n)`` (n the number operator on the pair space). Each fill
    holds the electron number to the tolerance of its chemical-potential
    bisection, which puts about 1e-7 into J along ``dG/dmu``, enough to move
    the eigenvectors of the small eigenvalues by more than 1e-6 (up to 7e-6
    measured); that multiple is removed (the exact J keeps the electron
    number). The nonzero spectra
    (|lambda| > 1e-8) agree to 1e-6, and the pair-space restriction of the J
    eigenvectors spans the K eigenspace of each eigenvalue at each q (1e-6;
    up to a phase for a single eigenvalue).

Green functions are never compared across numpy versions; the references
hold energies, labels and statuses only.
"""
import contextlib
import copy
import io
import itertools
import json
import logging
import os
import shutil
import tempfile
import unittest
from fractions import Fraction

import numpy as np

from tests import initial_fixtures
from tests import initial_reference_capture as capture
from tests.heavy_tests import heavy

#: the free-energy tolerance per site of the V4 multi-start
ENERGY_TOL = 1e-6
#: an ordered state lies below the para row by more than this per site
MARGIN = 1e-3
#: the V5 tolerance (free energies of the whole lattice)
V5_TOL = 1e-8
#: rows seeded without a spin component (module docstring)
SPINLESS_CHANNELS = ("para", "charge")


def _load(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


@contextlib.contextmanager
def _cwd(path):
    old = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


@contextlib.contextmanager
def _quiet():
    logger = logging.getLogger("qlms")
    level = logger.level
    logger.setLevel(logging.ERROR)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            yield
    finally:
        logger.setLevel(level)


def _read_toml(path):
    import tomli
    with open(path, "rb") as f:
        return tomli.load(f)


def _para(toml, subshape=None):
    """``(facts, para)``: the model facts and the converged para reference."""
    from hwave.initial import hints
    from hwave.initial.facts import model_facts
    f = model_facts(toml, subshape=subshape)
    with _quiet():
        para = hints.para_hf(f, hints.t_hint_bare(f.problem))
    if para.status != "converged":
        raise AssertionError("para step: {} ({})".format(para.status, para.skip_reason))
    return f, para


def _q_index(q, shape):
    return int(np.ravel_multi_index(tuple(int(Fraction(x) * n) for x, n in zip(q, shape)),
                                    tuple(shape)))


def _mesh(shape):
    from hwave.initial import cells
    return [cells.to_qfrac(q) for q in cells.commensurate_q(shape)]


def _onsite_labels(norb, spins=(0,)):
    """Every on-site orbital pair ``(a, s, b, t, 0, 0, 0)`` with ``s, t`` in
    ``spins``, in the order ``(s, a)`` major, ``(t, b)`` minor (the
    spin-block index ``s * norb + a`` of the RPA pair layout)."""
    idx = [(a, s) for s in spins for a in range(norb)]
    return np.array([(a, s, b, t, 0, 0, 0) for (a, s) in idx for (b, t) in idx],
                    dtype=np.int64)


def _lindhard_mesh(bands, labels, shape):
    """``chi0(q)`` on ``labels`` at every mesh q, ``(nvol, n, n)`` in the
    C-ordered flat q index of the mesh."""
    from hwave.solver import initial_facade as fa
    out = np.empty((int(np.prod(shape)), len(labels), len(labels)), dtype=np.complex128)
    for q in _mesh(shape):
        out[_q_index(q, shape)] = fa.static_lindhard(bands, q, labels)
    return out


def _leading(res):
    """The largest eigenvalue of a ``SolveResult`` (None without a pair)."""
    return max((lam for lam, _x in res.pairs), default=None)


def _sector_leading(res, sector):
    return max((lam for (lam, _x), s in zip(res.pairs, res.sectors) if s == sector),
               default=None)


class _Tmp(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="initial_validation_i2_")

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)

    def path(self, *parts):
        return os.path.join(self.root, *parts)


# ---------------------------------------------------------------------------
# V4
# ---------------------------------------------------------------------------

def _rpa_static_chi0(toml, H_up, mu, *, nmat=8192):
    """The static block of the RPA ``chi0q`` (general scheme, spin-free,
    ``coeff_tail = 0``) on the band ``H_up`` ``(nvol, norb, norb)`` at the
    chemical potential ``mu``: ``(nvol, norb, norb, norb, norb)``. The band is
    handed to the RPA solver as its one-body term, which it must reproduce
    exactly (``hwave_sc``'s mean-field band path)."""
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.rpa as sol_rpa
    from hwave.solver.eliashberg.io import (_chi0q_from_prepared_solver,
                                            _load_band_into_solver)
    doc = _read_toml(toml)
    par = doc["mode"]["param"]
    info_mode = {"mode": "RPA", "calc_scheme": "general",
                 "param": {"CellShape": par["CellShape"], "SubShape": [1, 1, 1],
                           "T": par["T"], "Ncond": par["Ncond"], "Nmat": nmat,
                           "coeff_tail": 0.0}}
    with _cwd(os.path.dirname(toml)), _quiet():
        rio = read_input_k.QLMSkInput(copy.deepcopy(doc["file"]["input"]))
        solver = sol_rpa.RPA(rio.get_param("ham"), {"print_level": 0}, info_mode)
        _load_band_into_solver(solver, np.ascontiguousarray(H_up))
        chi = _chi0q_from_prepared_solver(solver, mu, static_only=True)
    return np.asarray(chi)[0]


def _rpa_spinful_vertex(toml_so, *, exchange):
    """The spinful vertex of the spin-orbital RPA ``(nvol, nd^2, nd^2)``:
    the direct tensor plus, with ``exchange``, the on-site exchange crossing
    (:func:`hwave.solver.rpa.spinful_bubble_vertex`)."""
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.rpa as sol_rpa
    from hwave.solver import backend
    doc = _read_toml(toml_so)
    info_mode = dict(doc["mode"])
    info_mode["calc_scheme"] = "general"
    info_mode["_orbital_site_reported"] = True
    with _cwd(os.path.dirname(toml_so)), _quiet():
        rio = read_input_k.QLMSkInput(copy.deepcopy(doc["file"]["input"]))
        solver = sol_rpa.RPA(rio.get_param("ham"), {"print_level": 0}, info_mode)
        nvol, nd = solver.lattice.nvol, 2 * solver.norb
        ham = np.asarray(backend.to_host(solver.ham_info.ham_inter_q)).reshape(
            nvol, nd, nd, nd, nd)
        cross = solver.ham_info.spinful_exchange() if exchange else None
        return sol_rpa.spinful_bubble_vertex(ham, cross)


class TestV4(_Tmp):

    def model(self, name="v4", **kw):
        return capture.v4_input(self.path(name), **kw)

    def test_v4_n_c_and_oracle(self):
        """n_c = 2; static_lindhard on every on-site pair (ligand included)
        and on bond labels, against the label-form oracle (1e-10)."""
        from tests.test_initial_lindhard import _oracle_on_labels
        from hwave.solver import initial_facade as fa
        f, para = _para(self.model())
        self.assertEqual((f.n_c, f.interacting_orbitals), (2, (0, 1)))
        p = f.problem
        b = para.bands
        self.assertEqual(len(b.mu), 1)
        lab = _onsite_labels(3, spins=(0, 1))
        bond = np.array([(0, 0, 2, 0, 1, 0, 0), (2, 1, 1, 1, 0, -1, 0), (1, 0, 0, 0, 1, 1, 0),
                         (0, 0, 1, 1, -1, 0, 0)], dtype=np.int64)
        lab = np.concatenate([lab, bond])
        qs = [(0, 0, 0), (Fraction(1, 4), Fraction(1, 2), 0), (Fraction(3, 4), Fraction(1, 4), 0),
              (Fraction(1, 2), Fraction(1, 2), 0)]
        want = _oracle_on_labels(p, lab, np.asarray(para.H_para_k), b.mu[0], b.T, qs)
        for q in qs:
            got = fa.static_lindhard(b, q, lab)
            scale = max(1.0, float(np.max(np.abs(want[q]))))
            self.assertLessEqual(float(np.max(np.abs(got - want[q]))), 1e-10 * scale,
                                 msg=str(q))
            self.assertGreater(float(np.max(np.abs(got))), 1e-2)

    def test_v4_rpa_chi0q(self):
        """The RPA chi0q on the same para bands, every on-site pair with the
        ligand, every q: chi0[(a, b), (c, d)] = chi0q[a, b, c, d] to 1e-4
        (the truncation offset 6.2e-5 of Nmat = 8192 at T = 0.2)."""
        toml = self.model()
        f, para = _para(toml)
        p = f.problem
        H = np.asarray(para.H_para_k)
        n = p.norb_orig
        # the para bands are spin-free: the RPA receives the up block
        np.testing.assert_allclose(H[:, :n, :n], H[:, n:, n:], rtol=0, atol=1e-12)
        self.assertEqual(float(np.max(np.abs(H[:, :n, n:]))), 0.0)
        rpa = _rpa_static_chi0(toml, H[:, :n, :n], para.bands.mu[0])
        ours = _lindhard_mesh(para.bands, _onsite_labels(n), p.shape)
        err = float(np.max(np.abs(rpa.reshape(p.nvol, n * n, n * n) - ours)))
        offset = 1.0 / (np.pi ** 2 * para.T * 8192)
        self.assertLessEqual(err, 1e-4)
        self.assertLess(err, 1.5 * offset)
        # the pair-swapped reading does not agree: the layout is pinned
        swapped = rpa.transpose(0, 2, 1, 4, 3).reshape(p.nvol, n * n, n * n)
        self.assertGreater(float(np.max(np.abs(swapped - ours))), 1e-3)

    def test_v4_vertex_vs_rpa_stoner(self):
        """At every q the leading eigenvalue equals static_instability on the
        spin-orbital RPA vertex fed the same chi0: antisymmetrised with the
        Fock term, the direct (ring) vertex without it (1e-10)."""
        from hwave.initial import hints
        from hwave.solver import flex_diagnostics
        for fock in (True, False):
            with self.subTest(fock=fock):
                toml = self.model("n{}".format(fock), fock=fock)
                toml_so = self.model("so{}".format(fock), fock=fock, spin_orbital=True)
                f, para = _para(toml)
                p = f.problem
                P = hints.pair_space(f)
                qs = _mesh(p.shape)
                eig = hints.hint_eigenpairs(f, para, P, qs, top_k=64)
                C = _lindhard_mesh(para.bands, _onsite_labels(3, spins=(0, 1)), p.shape)
                W = _rpa_spinful_vertex(toml_so, exchange=fock)
                self.assertEqual(W.shape, (p.nvol, 36, 36))
                best = flex_diagnostics.static_instability(C, W, p.shape)
                lam_max = max(_leading(r) for r in eig.values() if r.pairs)
                self.assertAlmostEqual(best[0], lam_max, delta=1e-10)
                self.assertGreater(lam_max, 1.0)
                for q in qs:
                    iq = _q_index(q, p.shape)
                    lam_q = flex_diagnostics.static_instability(
                        C[iq:iq + 1], W[iq:iq + 1], (1, 1, 1))[0]
                    ours = _leading(eig[q])
                    if ours is None:
                        self.assertLessEqual(lam_q, 1e-9, msg=str(q))
                    else:
                        self.assertAlmostEqual(lam_q, ours, delta=1e-10, msg=str(q))

    def test_v4_kuroki_sc(self):
        """The S^z and charge eigenvalues at every q equal max eig(chi0 S) and
        max eig(-chi0 C) of hwave_sc's S and C matrices (density-only
        subcase and the full Kanamori model; 1e-10)."""
        from hwave.initial import hints
        from hwave.solver import initial_facade as fa
        from hwave.solver.eliashberg.io import _read_interaction_files
        from hwave.solver.eliashberg.vertex import _build_interaction_k, _build_sc_matrices
        cases = {"density": ("CoulombIntra", "CoulombInter"),
                 "kanamori": initial_fixtures.WANNIER_LIKE_TYPES}
        for name, types in cases.items():
            with self.subTest(name):
                toml = self.model(name, types=types)
                f, para = _para(toml)
                p = f.problem
                P = hints.pair_space(f)
                self.assertEqual(P.sectors, (hints.SECTOR_CHARGE, hints.SECTOR_SZ))
                qs = _mesh(p.shape)
                eig = hints.hint_eigenpairs(f, para, P, qs, top_k=64)
                with _cwd(os.path.dirname(toml)), _quiet():
                    _geom, _hr, inter = _read_interaction_files(_read_toml(toml))
                k = [np.linspace(0, 2 * np.pi, n, endpoint=False) for n in p.shape]
                inter_k = _build_interaction_k(*k, inter, 3)
                lab = _onsite_labels(3)
                compared = 0
                for q in qs:
                    S, C = (np.asarray(m).reshape(9, 9) for m in _build_sc_matrices(
                        inter_k, 3, *(int(x * n) for x, n in zip(q, p.shape))))
                    chi = fa.static_lindhard(para.bands, q, lab)
                    want = {hints.SECTOR_SZ: float(np.linalg.eigvals(chi @ S).real.max()),
                            hints.SECTOR_CHARGE: float(np.linalg.eigvals(-chi @ C).real.max())}
                    for sector, value in want.items():
                        got = _sector_leading(eig[q], sector)
                        if got is None:
                            self.assertLessEqual(value, 1e-9, msg=(sector, q))
                        else:
                            compared += 1
                            self.assertAlmostEqual(got, value, delta=1e-10, msg=(sector, q))
                self.assertGreaterEqual(compared, len(qs))

    def test_v4_su2_sectors_equal_full(self):
        """The charge and S^z sectors hold the nonzero spectrum of K on the
        full pair space, S^z three times (1e-12)."""
        from hwave.initial import hints
        from hwave.initial.symmetry import Su2Result
        f, para = _para(self.model())
        P = hints.pair_space(f)
        self.assertEqual(P.sectors, (hints.SECTOR_CHARGE, hints.SECTOR_SZ))
        full = hints.pair_space(f, su2=Su2Result(status="failed", max_residual=1.0,
                                                 method=None, sz_residual=None))
        self.assertEqual(full.sectors, (hints.SECTOR_FULL,))
        for q in ((0, 0, 0), (Fraction(1, 4), Fraction(1, 2), 0), (Fraction(1, 2), Fraction(1, 2), 0)):
            K = -hints.lindhard(para.bands, q, full) @ hints.vertex(f, q, full)
            w = np.sort(np.linalg.eigvals(K).real)
            parts = []
            for s, mult in ((hints.SECTOR_CHARGE, 1), (hints.SECTOR_SZ, 3)):
                Ks = -hints.lindhard(para.bands, q, P, sector=s) @ hints.vertex(f, q, P, sector=s)
                parts.extend(list(np.linalg.eigvals(Ks).real) * mult)
            parts = np.sort(np.array(parts))
            nz_w, nz_p = w[np.abs(w) > 1e-9], parts[np.abs(parts) > 1e-9]
            self.assertEqual(nz_w.size, nz_p.size)
            np.testing.assert_allclose(nz_w, nz_p, rtol=0, atol=1e-12)

    def test_v4_seed_shape(self):
        """(rho(h) - rho(0))/h along the conjugate field O = Gamma x of the
        leading pair (q = 0, S^z) equals chi0 Gamma x on every on-site slot,
        ligand included (h = 1e-6; 1e-5)."""
        from hwave.initial import hints, seed
        from hwave.solver import initial_facade as fa
        f, para = _para(self.model())
        p = f.problem
        P = hints.pair_space(f)
        q = (Fraction(0), Fraction(0), Fraction(0))
        res = hints.hint_eigenpairs(f, para, P, [q])[q]
        self.assertEqual(res.sectors[0], hints.SECTOR_SZ)
        lam_top = res.pairs[0][0]
        self.assertGreater(lam_top, 1.0)
        x = np.asarray(res.x_full[0])
        # the Hermitian phase: x_mbar = conj(x_m) (q = 0, on-site labels)
        lab = [tuple(int(v) for v in r) for r in P.labels]
        bar = np.array([lab.index((b, t, a, s, 0, 0, 0)) for (a, s, b, t, _rx, _ry, _rz) in lab])
        k = int(np.argmax(np.abs(x)))
        z = np.conj(x[bar[k]]) / x[k]
        x = x * np.sqrt(z)
        np.testing.assert_allclose(x[bar], np.conj(x), rtol=0, atol=1e-12)
        O = hints.vertex(f, q, P) @ x
        # O_k: the potential component of label (a, s, b, t, 0) is H(r = 0)[i, j]
        i = np.asarray(hints._nd_index(p, P.labels[:, 0], P.labels[:, 1]))
        j = np.asarray(hints._nd_index(p, P.labels[:, 2], P.labels[:, 3]))
        O_r = np.zeros((p.nvol, p.nd, p.nd), dtype=np.complex128)
        O_r[0, i, j] = O
        O_k = np.fft.ifftn(O_r.reshape(*p.shape, p.nd, p.nd), axes=(0, 1, 2),
                           norm="forward").reshape(p.nvol, p.nd, p.nd)
        H = np.asarray(para.H_para_k)
        th = hints.t_hint_bare(p)
        G0, _ = seed.density_at(f, H, None, 0.0, T_run=p.T, W0=th.W0)
        h = 1e-6
        Gh, _ = seed.density_at(f, H, O_k, h, T_run=p.T, W0=th.W0)
        ext = _onsite_labels(3, spins=(0, 1))
        ie = np.asarray(hints._nd_index(p, ext[:, 0], ext[:, 1]))
        je = np.asarray(hints._nd_index(p, ext[:, 2], ext[:, 3]))
        dG = (Gh - G0).reshape(p.nvol, p.nd, p.nd)
        # the rho-component of label m is G(-r)[j, i] (r = 0 here)
        got = dG[0, je, ie] / h
        cols = np.array([[tuple(v) for v in ext.tolist()].index(m) for m in lab])
        chi = fa.static_lindhard(para.bands, q, ext)[:, cols]
        want = chi @ O
        np.testing.assert_allclose(got, want, rtol=0, atol=1e-5)
        ligand = (ext[:, 0] == 2) | (ext[:, 2] == 2)
        self.assertGreater(float(np.max(np.abs(want[ligand]))), 1e-2)
        # on the pair space K x = -chi0 Gamma x = lambda x
        np.testing.assert_allclose(want[cols], -lam_top * x, rtol=0, atol=1e-10)

    def test_w_act_near_crossing(self):
        """The near-crossing fixture: W_act jumps from 0.4 to 4.6 where the
        wide band reaches the window edge; deterministic, invariant under a
        constant shift, and as recorded."""
        from hwave.initial import hints
        from hwave.solver import initial_facade
        ref = _load(capture.V4_REFERENCE)["w_act_near_crossing"]
        self.assertEqual([r["delta"] for r in ref], list(capture.W_ACT_DELTAS))
        for k, want in enumerate(ref):
            toml = capture.w_act_chain(self.path("w{}".format(k)), want["delta"])
            p = initial_facade.load_problem(toml)
            a, b = hints.t_hint_bare(p), hints.t_hint_bare(p)
            self.assertEqual((a.W_act, a.T_hint), (b.W_act, b.T_hint))
            for key in ("W_act", "T_hint", "D"):
                self.assertAlmostEqual(getattr(a, key) if key != "D" else a.D[0], want[key],
                                       delta=1e-12, msg=(key, want["delta"]))
            self.assertAlmostEqual(a.E_F[0], want["E_F"], delta=1e-12)
        self.assertAlmostEqual(ref[0]["W_act"], 0.4, delta=1e-12)
        self.assertAlmostEqual(ref[1]["W_act"], 4.6, delta=1e-12)
        self.assertAlmostEqual(ref[2]["W_act"], 4.6 - 1e-6, delta=1e-12)

    @heavy
    def test_v4_multistart(self):
        """prepare (small tier, SubShape 1 x 1 and 2 x 2), every run, collect:
        the lowest state and every spin-polarised row as frozen."""
        ref = _load(capture.V4_REFERENCE)
        cov = capture.prepare_run_collect(
            capture.v4_input(self.path("src")), self.path("work"), tier="small",
            subshapes=[list(b) for b in capture.V4_SUBSHAPES])
        self.assertEqual(cov["tier"]["name"], "small")
        self.assertEqual(cov["hint"]["status"], "computed")
        self.assertEqual(cov["hint"]["T_hint"], ref["T_hint"])
        got = {r["run_id"]: r for r in cov["rows"] if r["tried"]}
        self.assertEqual(sorted(got), sorted(ref["rows"]))
        families = {r["family"] for r in got.values()}
        self.assertTrue({"instability", "orbital", "baseline"} <= families, families)
        # the lowest converged row: its state and free energy
        conv = [r for r in got.values() if r["status"] == "converged"]
        lowest = min(conv, key=lambda r: r["free_energy_per_site"])
        self.assertEqual(lowest["final_label"], ref["lowest"]["final_label"])
        self.assertAlmostEqual(lowest["free_energy_per_site"],
                               ref["lowest"]["free_energy_per_site"], delta=ENERGY_TOL)
        para = [r for r in conv if r["family"] == "para"]
        self.assertTrue(para)
        for r in para:
            self.assertEqual(r["final_label"], "para")
            self.assertGreater(r["free_energy_per_site"] - lowest["free_energy_per_site"],
                               MARGIN, r["run_id"])
        para_f = min(r["free_energy_per_site"] for r in para)
        compared = labelled = 0
        for rid, want in ref["rows"].items():
            r = got[rid]
            self.assertEqual(r["channel"], want["channel"], rid)
            if r["channel"] in SPINLESS_CHANNELS or r["family"] == "orbital":
                self.assertIn(r["status"], ("converged", "not_converged"), rid)
                continue
            compared += 1
            self.assertEqual(r["status"], want["status"], rid)
            if want["free_energy_per_site"] is None:
                self.assertIsNone(r["free_energy_per_site"], rid)
                continue
            self.assertAlmostEqual(r["free_energy_per_site"], want["free_energy_per_site"],
                                   delta=ENERGY_TOL, msg=rid)
            if r["family"] != "random" and want["free_energy_per_site"] < para_f - MARGIN:
                labelled += 1
                self.assertEqual(r["final_label"], want["final_label"], rid)
        self.assertGreaterEqual(compared, 10)
        self.assertGreaterEqual(labelled, 3)


# ---------------------------------------------------------------------------
# V5
# ---------------------------------------------------------------------------

class TestV5(_Tmp):

    def test_v5_moments(self):
        """The Neel runs along x, y and z at lambda = +-0.3 and 0: as frozen,
        E_x = E_y, lambda -> -lambda invariant, equal at lambda = 0; each
        moment stays along its seed axis."""
        ref = _load(capture.V5_REFERENCE)["runs"]
        got = capture.v5_energies(self.path("v5"))
        self.assertEqual(sorted(got), sorted(ref))
        for soc, runs in got.items():
            for axis, r in runs.items():
                self.assertTrue(r["converged"], (soc, axis))
                self.assertAlmostEqual(r["free_energy"], ref[soc][axis]["free_energy"],
                                       delta=V5_TOL, msg=(soc, axis))
                m = np.abs(np.array(r["staggered"]))
                k = "xyz".index(axis)
                self.assertGreater(m[k], 0.1, (soc, axis))
                self.assertLess(float(np.max(np.delete(m, k))), 1e-8, (soc, axis))
        lam = capture.v5_soc_key(capture.V5_SOC)
        neg = capture.v5_soc_key(-capture.V5_SOC)
        zero = capture.v5_soc_key(0.0)

        def energy(soc, axis):
            return got[soc][axis]["free_energy"]

        for soc in (lam, neg):
            self.assertAlmostEqual(energy(soc, "x"), energy(soc, "y"), delta=V5_TOL)
        for axis in "xyz":
            self.assertAlmostEqual(energy(lam, axis), energy(neg, axis), delta=V5_TOL)
        self.assertAlmostEqual(energy(zero, "x"), energy(zero, "z"), delta=V5_TOL)
        self.assertAlmostEqual(energy(zero, "y"), energy(zero, "z"), delta=V5_TOL)
        # the coupling splits in-plane and out-of-plane moments (the frozen
        # values say which is lower; the test does not)
        self.assertGreater(abs(energy(lam, "x") - energy(lam, "z")), 1e-3)

    def test_v5_hints_and_refusals(self):
        """Small tier on the V5 model: time reversal with the coupling passes,
        no gauge warning, the leading hint is a TR-odd spin pattern; a UHFr
        target is refused."""
        from hwave.initial import labels, patterns
        from hwave.initial.enumerate import candidate_states
        from hwave.initial.facts import model_facts
        from hwave.initial.writer import InitialStateError, initial_state
        toml = initial_fixtures.so_two_orbital(self.path("so"), capture.V5_SOC,
                                               **capture.V5_MODEL)
        with _quiet():
            cs = candidate_states(toml, tier="small", n_random=0)
        self.assertEqual(cs.tier, "small")
        self.assertIs(cs.model["tr"]["invariant"], True)
        self.assertEqual(cs.model["su2"]["status"], "not applicable (spin-orbital)")
        self.assertEqual(cs.model["hint"]["status"], "computed")
        self.assertFalse([w for w in cs.model.get("warnings", ()) if "gauge" in w.lower()])
        inst = sorted((r for r in cs.rows if r["family"] == "instability"),
                      key=lambda r: r["rank"])
        self.assertTrue(inst)
        top = inst[0]
        self.assertEqual(top["rank"], 0)
        self.assertGreater(top["lambda"], 1.0)
        run = {r.run_id: r for r in cs.runs}[top["run_id"]]
        fB = model_facts(toml, subshape=run.subshape)
        O = patterns.realise(run.candidate.pattern, fB, run.subshape)
        lab = labels.pattern_label(O, fB, q=None, eta=None, phase=None, sign=None,
                                   spiral=False)
        self.assertEqual(lab.tr, "odd")
        self.assertIs(lab.gauge_dependent, False)
        self.assertTrue(lab.name.startswith("spin-"), lab.name)
        # UHFr has no spin-orbital mode
        uhfr = self.path("uhfr.toml")
        with open(uhfr, "w") as fh:
            fh.write(capture.dump_toml({"mode": {"mode": "UHFr", "param": {
                "Nsite": 16, "Ncond": 32, "T": 0.2}},
                "file": {"input": {"path_to_input": "", "interaction": {}},
                         "output": {"path_to_output": "output"}}}))
        out = self.path("x.def")
        with self.assertRaises(InitialStateError) as cm:
            initial_state(toml, capture.v5_neel_candidate("z"), out=out,
                          site_map=[(x, y, 0, a) for a in range(2) for y in range(4)
                                    for x in range(4)][:16], uhfr_toml=uhfr)
        self.assertIn("UHFr has no spin-orbital mode", str(cm.exception))
        self.assertFalse(os.path.exists(out))


# ---------------------------------------------------------------------------
# V6
# ---------------------------------------------------------------------------

def _v6_bond_model(dest, *, fock, subshape=(2, 2, 1)):
    """Square 8 x 8, hopping 1, the nearest-neighbour CoulombInter V = 1.2
    only, Ncond = 48, T = 0.1."""
    V = {((1, 0, 0), 0, 0): 1.2, ((-1, 0, 0), 0, 0): 1.2,
         ((0, 1, 0), 0, 0): 1.2, ((0, -1, 0), 0, 0): 1.2}
    return initial_fixtures.write_input(
        dest, cellshape=(8, 8, 1), transfer=initial_fixtures.square_transfer(),
        interactions={"CoulombInter": V},
        param={"Ncond": 48, "T": 0.1, "SubShape": list(subshape), "IterationMax": 2000},
        mode=None if fock else {"flag_fock": False})


def _has_bond(run):
    return any(tuple(e.r) != (0, 0, 0) for e in run.candidate.pattern.entries)


class TestV6(_Tmp):

    def test_v6_coulombinter(self):
        """Bond-order instability rows with the Fock term, none without it;
        a q = (1/2, 1/2) bond-order seed of prepare survives the round trip."""
        from hwave.initial import cli
        from hwave.initial.enumerate import candidate_states
        from hwave.results import read_uhfk_green
        from tests.test_initial_writer import _load as load_npz, _run_hwave, _set_run
        rows = {}
        for fock in (True, False):
            toml = _v6_bond_model(self.path("f{}".format(fock)), fock=fock)
            with _quiet():
                cs = candidate_states(toml, tier="small", n_random=0)
            inst = [r for r in cs.rows if r["family"] == "instability"]
            self.assertTrue(inst, fock)
            runs = {r.run_id: r for r in cs.runs}
            bond = [r for r in inst if r["run_id"] in runs and _has_bond(runs[r["run_id"]])]
            listed_bond = [r for r in inst if "bond-" in (r.get("label") or "")
                           or "bond-" in (r.get("pattern_label") or "")]
            if fock:
                self.assertGreater(len(bond), 0)
                rows[fock] = (toml, bond, runs)
            else:
                self.assertEqual(bond, [])
                self.assertEqual(listed_bond, [])
                self.assertTrue(all(not _has_bond(r) for r in cs.runs))
        toml, bond, runs = rows[True]
        target = [r for r in bond if r["q"] == ["1/2", "1/2", "0"] or
                  tuple(r["q"]) == ("1/2", "1/2", "0")]
        self.assertTrue(target, [r["q"] for r in bond])
        rid = target[0]["run_id"]
        work = self.path("work")
        with _quiet():
            rc = cli.prepare(toml, work, tier="small", stdout=io.StringIO())
        self.assertEqual(rc, 0)
        manifest = _load(os.path.join(work, "manifest.json"))
        entry = {e["run_id"]: e for e in manifest["runs"]}[rid]
        run_toml = os.path.join(work, entry["run_dir"], entry["input_toml"])
        written = load_npz(os.path.join(work, entry["run_dir"], "initial.npz"))
        _set_run(run_toml)
        with _quiet():
            _run_hwave(run_toml)
        green = read_uhfk_green(run_toml)
        np.testing.assert_allclose(green.green, written["green"], rtol=0, atol=1e-12)
        np.testing.assert_allclose(green.sublattice.green, written["green_sublattice"],
                                   rtol=0, atol=1e-12)
        # a q = (1/2, 1/2) bond modulation: the x bonds of the two sublattices
        # of a row of the 2 x 2 cell differ
        G = written["green_sublattice"]
        self.assertGreater(abs(G[0, 0, 0, 0, 1] - G[0, 0, 1, 0, 0]), 1e-4)

    def test_v6_pairhop_jacobian(self):
        """The finite-difference Jacobian of the SCF map against K(q) at
        q = 0 and (1/2, 0, 0): nonzero spectra (1e-6) and eigenvectors on
        the pair space (1e-6), with and without the Fock term."""
        for fock in (True, False):
            with self.subTest(fock=fock):
                self.check_jacobian(fock)

    def check_jacobian(self, fock):
        from hwave.initial import cells, hints
        from hwave.solver import initial_facade as fa
        tables = {"CoulombIntra": {((0, 0, 0), 0, 0): 2.0, ((0, 0, 0), 1, 1): 1.5},
                  "PairHop": {((0, 0, 0), 0, 1): 0.3 + 0.1j, ((1, 0, 0), 0, 1): 0.05 - 0.02j}}
        toml = initial_fixtures.golden(self.path("ph{}".format(fock)), cellshape=(3, 3, 1),
                                       fock=fock, tables=tables, T=0.1,
                                       param={"IterationMax": 3000})
        f, para = _para(toml)
        self.assertEqual(para.T, 0.1)
        B = (3, 1, 1)
        pB = fa.load_problem(toml, subshape=B)
        gs = (pB.nvol, pB.ns, pB.norb, pB.ns, pB.norb)
        N, nd = pB.nvol, pB.nd

        def scf_map(G):
            H = np.asarray(pB.ham_trans_k) + fa.hf_potential(
                G.reshape(gs), pB.tables, include_fock=pB.iflag_fock, blocks=pB.blocks)
            return fa.green_function(H, pB.T, pB.blocks, shape=pB.shape, ns=pB.ns,
                                     norb=pB.norb, ene_cutoff=pB.ene_cutoff).reshape(N, nd, nd)

        G0 = fa.green_function(hints.fold_para(para.H_para_k, f, B), pB.T, pB.blocks,
                               shape=pB.shape, ns=pB.ns, norb=pB.norb,
                               ene_cutoff=pB.ene_cutoff).reshape(N, nd, nd)
        self.assertLess(float(np.max(np.abs(scf_map(G0) - G0))), 1e-10)
        # the slots of the folded cell as primitive labels (a, s, b, t, r) at
        # an anchor offset (the placement of tests/test_initial_vertex.py)
        nphys, norb_orig = f.problem.norb_phys, pB.norb_orig
        offsets = list(itertools.product(range(B[0]), range(B[1]), range(B[2])))

        def idx(s, a, o):
            ir = o[0] + B[0] * (o[1] + B[1] * o[2])
            return s * pB.norb + a + norb_orig * ir

        def slots(labels, q):
            """``(flat slot, conj phase)`` per (label, anchor offset)."""
            out = []
            for m, (a, s, b, t, *r) in enumerate(labels):
                for oa in offsets:
                    R = [oa[i] - r[i] for i in range(3)]
                    D = tuple((R[i] // B[i]) % pB.shape[i] for i in range(3))
                    oR = tuple(R[i] % B[i] for i in range(3))
                    ph = np.exp(2j * np.pi * float(sum(Fraction(x) * y for x, y in zip(q, R))))
                    flat = np.ravel_multi_index(D + (idx(t, b, oR), idx(s, a, oa)),
                                                tuple(pB.shape) + (nd, nd))
                    out.append((m, int(flat), np.conj(ph)))
            return out

        def restrict(vec, labels, q):
            out = np.zeros(len(labels), dtype=np.complex128)
            for m, flat, cph in slots(labels, q):
                out[m] += cph * vec[flat]
            return out / len(offsets)

        # J in the stored-G coordinates (rho = G^T is a relabelling): the
        # complex-linear extension of the derivative along Hermitian directions
        neg = hints._neg_r(pB.shape)
        n = N * nd * nd
        J = np.zeros((n, n), dtype=np.complex128)
        step = 1e-5
        # slot k and its partner kbar (G(r)_ij <-> G(-r)_ji) share the two
        # Hermitian directions E + E^dag and i(E - E^dag): each derivative is
        # taken once per pair
        r_k, i_k, j_k = np.unravel_index(np.arange(n), (N, nd, nd))
        partner = np.ravel_multi_index((neg[r_k], j_k, i_k), (N, nd, nd))

        def derivative(Hm):
            return ((scf_map(G0 + step * Hm) - scf_map(G0 - step * Hm)) / (2 * step)).reshape(-1)

        for k in range(n):
            kb = int(partner[k])
            if kb < k:
                continue
            E = np.zeros(n, dtype=np.complex128)
            E[k] = 1.0
            E = E.reshape(N, nd, nd)
            Ed = np.conj(np.swapaxes(E[neg], -1, -2))
            d1 = derivative((E + Ed) / 2)
            if kb == k:
                # a self-partner slot (r = -r, i = j) is real: E = E^dag
                J[:, k] = d1
                continue
            d2 = derivative((E - Ed) / 2j)
            # E_k = H1 + i H2 and E_kbar = H1 - i H2
            J[:, k] = d1 + 1j * d2
            J[:, kb] = d1 - 1j * d2
        # Every fill holds the electron number to the tolerance of the mu
        # bisection (xtol 2e-12), which enters each column of J as a multiple
        # of u = dG/dmu, about 1e-7 at this step. The exact J keeps N
        # (tr J[:, k] at r = 0 vanishes), so that multiple is removed: u is
        # the response chi0 n to a uniform potential on every slot of the cell.
        everything = np.array(
            [(a, s, b, t) + tuple(int(v) - (c if v > c // 2 else 0)
                                  for v, c in zip(rr, f.cellshape))
             for rr in itertools.product(*(range(c) for c in f.cellshape))
             for a in range(nphys) for s in range(2) for b in range(nphys) for t in range(2)],
            dtype=np.int64)
        # (only slots inside one solver block: without the Fock term the
        # spins are separate blocks)
        p1 = f.problem
        block = np.full(p1.nd, -1)
        for kb, blk in enumerate(p1.blocks.block_info):
            block[np.asarray(blk)] = kb
        same = (block[np.asarray(hints._nd_index(p1, everything[:, 0], everything[:, 1]))]
                == block[np.asarray(hints._nd_index(p1, everything[:, 2], everything[:, 3]))])
        everything = everything[same]
        zero = (Fraction(0), Fraction(0), Fraction(0))
        dens_all = ((everything[:, 0] == everything[:, 2]) & (everything[:, 1] == everything[:, 3])
                    & np.all(everything[:, 4:] == 0, axis=1)).astype(float)
        u_rho = fa.static_lindhard(para.bands, zero, everything) @ dens_all
        u = np.zeros(n, dtype=np.complex128)
        for m, flat, _cph in slots(everything, zero):
            self.assertEqual(u[flat], 0)
            u[flat] = u_rho[m]
        trace0 = np.zeros(n)
        trace0[[np.ravel_multi_index((0, 0, 0, i, i), tuple(pB.shape) + (nd, nd))
                for i in range(nd)]] = 1.0
        J = J - np.outer(u, trace0 @ J) / (trace0 @ u)
        wJ, vJ = np.linalg.eig(J)
        keep = np.abs(wJ) > 1e-8
        wJ, vJ = wJ[keep], vJ[:, keep]
        self.assertLess(float(np.max(np.abs(wJ.imag))), 1e-6)
        # K(q) per q of the cell; at q = 0 with the electron number fixed
        P = hints.pair_space(f)
        dens = ((P.labels[:, 0] == P.labels[:, 2]) & (P.labels[:, 1] == P.labels[:, 3])
                & np.all(P.labels[:, 4:] == 0, axis=1)).astype(float)
        spectra = []
        for q in cells.commensurate_q(B):
            qf = cells.to_qfrac(q)
            C = hints.lindhard(para.bands, qf, P)
            if all(x == 0 for x in qf):
                C = C - np.outer(C @ dens, dens @ C) / (dens @ C @ dens)
            K = -C @ hints.vertex(f, qf, P)
            w, v = np.linalg.eig(K)
            nz = np.abs(w) > 1e-8
            spectra.append((qf, w[nz], v[:, nz]))
        wK = np.concatenate([s[1] for s in spectra])
        self.assertEqual(wK.size, wJ.size)
        np.testing.assert_allclose(np.sort(wJ.real), np.sort(wK.real), rtol=0, atol=1e-6)
        self.assertGreater(float(np.max(np.abs(wK))), 0.1)

        def basis(M):
            # a set degenerate across q (K(q) and K(-q)) comes back from eig as
            # any mixture of its q components, so the restriction of such a
            # set to one q has the rank of that q's part only, up to the
            # noise of the difference quotients (below 1e-8 relative)
            u, sv, _ = np.linalg.svd(M, full_matrices=False)
            return u[:, sv > 1e-5 * sv[0]]

        # per eigenvalue (degenerate across q: K(q) and K(-q) share their
        # spectrum), the restriction of J's eigenspace to each q spans K(q)'s
        checked = 0
        for qf, w, v in spectra:
            for lam in w:
                kq = np.abs(w - lam) < 1e-6
                kj = np.abs(wJ - lam) < 1e-6
                total = sum(int(np.sum(np.abs(w2 - lam) < 1e-6)) for _q, w2, _v in spectra)
                self.assertEqual(total, int(kj.sum()), (qf, lam))
                QK = basis(v[:, kq])
                QJ = basis(np.stack([restrict(vJ[:, i], P.labels, qf)
                                     for i in np.nonzero(kj)[0]], axis=1))
                self.assertEqual(QK.shape[1], QJ.shape[1], (qf, lam))
                resid = QJ - QK @ (QK.conj().T @ QJ)
                self.assertLessEqual(float(np.linalg.norm(resid, 2)), 1e-6, (qf, lam))
                checked += 1
        self.assertEqual(checked, wK.size)


if __name__ == "__main__":
    unittest.main()
