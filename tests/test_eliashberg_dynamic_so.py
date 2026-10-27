"""The spin-orbital dynamic Eliashberg solver (``[mode] enable_spin_orbital``
with ``[eliashberg] frequency = "dynamic"``, issue #284).

The gap ``Delta_AB(k, i w_n)`` is solved on the fermion-antisymmetric space
``Delta_AB(k, i w_n) = -Delta_BA(-k, -i w_n)`` (``-i w_n`` is the index
``Nmat - 1 - n``) as ``E^dagger K E``, with ``K`` the frequency-resolved
kernel of the non-spin-orbital dynamic solver at ``nd = 2 N`` spin-orbitals
and the spin-orbital pairing vertex at every bosonic frequency.

* ``TestKernelAgainstValidatedKernel``: the production operator equals the
  ``BondPairKernel`` application that ``tests/test_so_pairing_ed.py``
  validates against exact diagonalization with spin-orbit coupling, on the
  ``L = 3`` model of that test (free Green function in the fast gate; the
  exact interacting one in the heavy test).
* ``TestSU2Limit``: with spin-independent hopping and a rotation-invariant
  interaction the spectrum on the antisymmetric space is every non-spin-orbital
  dynamic singlet eigenvalue once and every triplet eigenvalue three times
  (both in their combined-parity sectors, i.e. even and odd frequency), and
  nothing else.
* ``TestStaticLimit``: a frequency-independent vertex reproduces the
  spectrum of the static spin-orbital solver (#269, ``g2_tail = false``), with
  frequency-flat gaps.
* ``TestFlexMode`` / ``TestCalcMode``: the public entry reads the archives
  of a spin-orbital FLEX run (every bosonic frequency, the vertex
  consistency checked at every one of them) or computes the spinful bubble,
  and its eigenvalues equal an independent dense construction.
* ``TestOutputs``: every output gap is fermion-antisymmetric to round-off;
  the eigenvalue.dat / gap.dat / gap_dynamic.npz schema and weights.
* ``TestRefusals``: each unsupported combination is refused through the
  public entry, with the key named, before any output.
* ``TestGPU``: host and device operators agree (skipped without CuPy).
"""
import contextlib
import logging
import os
import shutil
import tempfile
import types
import unittest
from unittest import mock

import numpy as np

import hwave.sc as sc
from hwave.solver import backend, bubble
from hwave.solver import eliashberg_dynamic as dyn
from hwave.solver import eliashberg_dynamic_so as dso
from hwave.solver.so_pairing import so_pairing_vertex, so_pairing_vertex_parts
from tests.heavy_tests import heavy
import tests.test_sc_spin_orbital as S
import tests.test_so_pairing_ed as E
import tests.test_so_pairing_vertex as T


@contextlib.contextmanager
def _quiet():
    prev = logging.root.manager.disable
    logging.disable(logging.WARNING)
    try:
        yield
    finally:
        logging.disable(prev)


def _dense(A, n):
    return np.asarray(A.matmat(np.eye(n, dtype=complex)))


def _reduced_dense(A, basis):
    return _dense(basis.reduce(A), basis.dim)


def _multiset_match(want, got, atol):
    """Largest distance of a greedy one-to-one matching of ``got`` onto
    ``want`` (``inf`` when the sizes differ)."""
    want, got = list(np.asarray(want)), list(np.asarray(got))
    if len(want) != len(got):
        return float("inf")
    worst = 0.0
    for w in sorted(want, key=lambda z: (-z.real, z.imag)):
        j = int(np.argmin([abs(w - g) for g in got]))
        worst = max(worst, abs(w - got[j]))
        got.pop(j)
    return worst


def _bubble_w(g, beta):
    """The general bubble at every bosonic frequency, ``(nmat, nvol, nd^2,
    nd^2)``, of a Green function in the sc layout."""
    nd = g.shape[0]
    shape = g.shape[2:5]
    nvol = int(np.prod(shape))
    gk = g.reshape(nd, nd, nvol, g.shape[-1]).transpose(3, 2, 0, 1)[None].copy()
    c = bubble.dense_bubble(gk, None, beta, spatial_shape=shape, scheme="general")[0]
    return c.reshape(c.shape[0], nvol, nd * nd, nd * nd)


def _rpa(chi0_w, gam):
    nd2 = gam.shape[-1]
    eye = np.eye(nd2)[None]
    return np.stack([np.linalg.solve(eye + c @ gam, c) for c in chi0_w])


# ---------------------------------------------------------------------------
# the production operator against the ED-validated kernel
# ---------------------------------------------------------------------------

def _free_green_ed_model(nmat):
    """The free Green function of the ``L = 3`` spin-mixing chain of
    ``tests/test_so_pairing_ed.py`` (its ``TRANSFER`` and ``SOC`` rows), in
    the sc layout ``(ND, ND, L, 1, 1, nmat)``."""
    L, NP_, ND_ = E.L, E.NP, E.ND
    k = 2 * np.pi * np.arange(L) / L
    h = np.zeros((L, ND_, ND_), complex)
    for (r, a, b), t in E.TRANSFER.items():
        for s in range(2):
            h[:, s * NP_ + a, s * NP_ + b] += t * np.exp(1j * k * r)
    for (r, a, b), v in E.SOC.items():
        h[:, a, NP_ + b] += v * np.exp(1j * k * r)
        h[:, NP_ + b, a] += np.conj(v) * np.exp(-1j * k * r)
    wn = (2 * np.arange(nmat) + 1 - nmat) * np.pi / E.BETA
    g = np.empty((ND_, ND_, L, 1, 1, nmat), complex)
    for i, w in enumerate(wn):
        g[:, :, :, 0, 0, i] = np.linalg.inv(
            (1j * w + E.MU) * np.eye(ND_)[None] - h).transpose(1, 2, 0)
    return g


class TestKernelAgainstValidatedKernel(unittest.TestCase):
    """``tests/test_so_pairing_ed.py`` validates ``_kernel`` (the
    ``BondPairKernel`` fed with ``so_pairing_vertex_parts`` at every bosonic
    frequency, read at ``+q`` and ``+nu``) against exact diagonalization with
    spin-orbit coupling at first and second order. The production operator of
    the spin-orbital dynamic solver applied to the same inputs is that kernel
    to round-off -- on general gaps and on the antisymmetric space the
    solver iterates."""

    NMAT = 16

    def _compare(self, g, lam):
        nmat = g.shape[-1]
        gam = lam * E._gamma_onsite()
        chi_w = _bubble_w(g, E.BETA)
        A, basis = dso.build_operator(gam, chi_w, g, E.BETA)
        K = E._kernel(g, lam, nmat)
        # non-vacuity: spin mixing in G and a sizeable fluctuation term
        self.assertGreater(np.linalg.norm(g[:E.NP, E.NP:]), 1e-1 * np.linalg.norm(g))
        _, V_fl = so_pairing_vertex_parts(gam, chi_w[nmat // 2 + 1], E.ND,
                                          (E.L, 1, 1))
        V_b, _ = so_pairing_vertex_parts(gam, chi_w[0], E.ND, (E.L, 1, 1))
        # V_fluct / V_bare grows as lam (measured 0.14 for the free G at lam = 1)
        self.assertGreater(np.abs(V_fl).max(), 1e-3 * lam * np.abs(V_b).max())
        rng = np.random.default_rng(5)
        for _ in range(3):
            x = rng.standard_normal(basis.n) + 1j * rng.standard_normal(basis.n)
            want = K.matvec(x)
            got = A.matvec(x)
            self.assertLess(np.linalg.norm(got - want), 1e-13 * np.linalg.norm(want))
            # what the solver iterates: E^dagger K E, mapped back = P K on range(P)
            c = basis.to_reduced(x)
            red = basis.to_full(basis.reduce(A).matvec(c)).reshape(basis.shape)
            y = K.matvec(basis.to_full(c)).reshape(basis.shape)
            ref = dso.antisymmetrize(y)
            self.assertLess(np.linalg.norm(red - ref), 1e-13 * np.linalg.norm(ref))
            # the ED test's own projection is the same map
            ref_ed = E._antisym(y[:, :, :, 0, 0, :])
            np.testing.assert_allclose(ref[:, :, :, 0, 0, :], ref_ed, rtol=0,
                                       atol=1e-15 * np.abs(ref_ed).max())

    def test_free_green_function(self):
        self._compare(_free_green_ed_model(self.NMAT), 1.0)

    @heavy
    def test_exact_interacting_green_function(self):
        g = E._green_sc(E._ed(E.H)[1], self.NMAT)
        self._compare(g, E.H)


# ---------------------------------------------------------------------------
# SU(2) limit and static limit (dense, on the chain of test_so_pairing_vertex)
# ---------------------------------------------------------------------------

_NMAT_DENSE = 16
_WN_DENSE = (2 * np.arange(_NMAT_DENSE) + 1 - _NMAT_DENSE) * np.pi / T.BETA
_SU2_DECL_TYPES = (("CoulombIntra", 1.2), ("CoulombInter", 0.6), ("Hund", 0.25),
                   ("Exchange", 0.25), ("PairHop", 0.2))


def _green_chain(h):
    n = h.shape[1]
    g = np.empty((n, n, T.NX, 1, 1, _NMAT_DENSE), complex)
    for i, w in enumerate(_WN_DENSE):
        g[:, :, :, 0, 0, i] = np.linalg.inv(
            (1j * w + T.MU) * np.eye(n)[None] - h).transpose(1, 2, 0)
    return g


def _su2_decl():
    decl = {}
    for t, v in _SU2_DECL_TYPES:
        decl.update(T._onsite_decl(t, v))
    return decl


class TestSU2Limit(unittest.TestCase):
    """Spin-independent hopping, a rotation-invariant interaction set: the
    spectrum of ``E^dagger K E`` is the non-spin-orbital dynamic singlet
    spectrum once and the triplet spectrum three times, each restricted to its
    combined-parity sector ``Delta_ab(k, w_n) = +/- Delta_ba(-k, -w_n)``
    (which contains both frequency parities), and no extra mode."""

    def test_spectrum_is_singlet_once_and_triplet_three_times(self):
        decl = _su2_decl()
        gam = T._gamma(decl)
        h = T._h_orbital()
        g = _green_chain(h)
        np_, nd = T.NP, T.ND
        hso = np.zeros((T.NX, nd, nd), complex)
        hso[:, :np_, :np_] = h
        hso[:, np_:, np_:] = h
        gso = _green_chain(hso)
        # non-spin-orbital side: the general RPA vertex per bosonic index
        chi0 = _bubble_w(g, T.BETA).reshape(_NMAT_DENSE, T.NX, np_, np_, np_, np_)
        chi0q = np.ascontiguousarray(chi0.transpose(2, 3, 4, 5, 1, 0)).reshape(
            np_, np_, np_, np_, T.NX, 1, 1, _NMAT_DENSE)
        kx = 2 * np.pi * np.arange(T.NX) / T.NX
        inter_k = sc._build_interaction_k(kx, np.zeros(1), np.zeros(1), decl, np_)
        G2n = dyn.calc_g2_dynamic(g, T.BETA)
        nshape = (np_, np_, T.NX, 1, 1, _NMAT_DENSE)
        nn = int(np.prod(nshape))
        want = []
        for eta, s, mult in (("singlet", 1.0, 1), ("triplet", -1.0, 3)):
            with _quiet():
                Vn = np.stack([sc._compute_vertices_general(
                    chi0q, inter_k, np_, T.NX, 1, 1, _NMAT_DENSE, pairing_type=eta,
                    static_index=m) for m in range(_NMAT_DENSE)], axis=-1)
            I = np.eye(nn, dtype=complex)
            Kn = np.stack([dyn.eliashberg_kernel_dynamic(
                Vn, G2n, I[:, j].reshape(nshape), np_, T.BETA).ravel()
                for j in range(nn)], axis=1)
            Pn = np.stack([dyn._project_combined_parity_dynamic(
                I[:, j].reshape(nshape), s).ravel() for j in range(nn)], axis=1)
            w, U = np.linalg.eigh(0.5 * (Pn + Pn.conj().T))
            Q = U[:, w > 0.5]
            self.assertLess(np.linalg.norm(Kn @ Pn - Pn @ Kn), 1e-12 * np.linalg.norm(Kn))
            want += [np.linalg.eigvals(Q.conj().T @ Kn @ Q)] * mult
        want = np.concatenate(want)
        # spin-orbital side: the production operator
        chi_w = _rpa(_bubble_w(gso, T.BETA), gam[0])
        A, basis = dso.build_operator(gam, chi_w, gso, T.BETA)
        self.assertEqual(basis.n_self, 0)
        got = np.linalg.eigvals(_reduced_dense(A, basis))
        scale = np.abs(want).max()
        self.assertGreater(scale, 1e-2)
        self.assertEqual(len(got), len(want))
        self.assertLess(_multiset_match(want, got, 0), 1e-11 * scale)
        # the leading SO eigenvalue is a (threefold) triplet here
        lead = got[np.argmax(got.real)]
        self.assertEqual(int(np.sum(np.abs(got - lead) < 1e-10 * scale)), 3)


class TestStaticLimit(unittest.TestCase):
    """A frequency-independent susceptibility (the zero-frequency one at every
    bosonic index) gives a frequency-flat vertex, ``V(tau) = V delta(tau)``:
    the nonzero spectrum on the dynamic antisymmetric space is the spectrum of
    the static spin-orbital solver (#269) with ``g2_tail = false``, and those
    eigenvectors are frequency-flat (pure even frequency)."""

    def test_static_spectrum(self):
        decl = _su2_decl()
        gam = T._gamma(decl)
        gso = _green_chain(T._h_spin_orbital(0.3))
        nd = T.ND
        chi_s = _rpa(_bubble_w(gso, T.BETA)[_NMAT_DENSE // 2][None], gam[0])[0]
        chi_flat = np.broadcast_to(chi_s, (_NMAT_DENSE,) + chi_s.shape).copy()
        A, basis = dso.build_operator(gam, chi_flat, gso, T.BETA)
        vals, vecs = np.linalg.eig(_reduced_dense(A, basis))
        # static spin-orbital solver (production operator, #269)
        Vs = so_pairing_vertex(gam, chi_s, nd, (T.NX, 1, 1))
        G2s = sc._calc_g2(gso, T.BETA, tail=False)
        Ks, _ = sc._make_kernel_operator(Vs, G2s, nd, T.NX, 1, 1)
        sbasis = sc._SOAntisymmetricBasis((nd, nd, T.NX, 1, 1))
        with _quiet():
            svals, _ = sc._so_solve_eigenvalue(Ks, sbasis, sbasis.dim, "arnoldi")
        scale = np.abs(svals).max()
        snz = svals[np.abs(svals) > 1e-10 * scale]
        dnz = np.abs(vals) > 1e-10 * scale
        self.assertGreater(len(snz), 10)
        self.assertLess(_multiset_match(snz, vals[dnz], 0), 1e-12 * scale)
        for i in np.where(dnz)[0]:
            gap = basis.to_full(vecs[:, i]).reshape(basis.shape)
            w = dso.gap_weights(gap, T.NP)
            self.assertLess(w["parity"][1], 1e-18)
            self.assertLess(np.abs(gap - gap[..., :1]).max(), 1e-10)


# ---------------------------------------------------------------------------
# the public entry
# ---------------------------------------------------------------------------

_NMAT_E2E = 8


def _run(inp, spy_solve=True):
    """``calc_eliashberg`` with the operator, basis and eigenpairs captured."""
    cap = {}
    real_solve = sc._so_solve_eigenvalue
    real_green = sc._build_bond_green

    def solve(A, basis, *a, **kw):
        cap["A"], cap["basis"] = A, basis
        out = real_solve(A, basis, *a, **kw)
        cap["vals"], cap["vecs"] = out
        return out

    def green(*a, **kw):
        out = real_green(*a, **kw)
        cap["green"] = np.array(out.full_sc)
        return out

    with _quiet(), mock.patch.object(sc, "_so_solve_eigenvalue", side_effect=solve), \
            mock.patch.object(sc, "_build_bond_green", side_effect=green):
        cap["lam"] = sc.calc_eliashberg(inp)
    return cap


def _dyn_eli(**kw):
    eli = {"frequency": "dynamic", "solver_mode": "eigenvalue", "num_eigenvalues": 6}
    eli.update(kw)
    return eli


def _read_table(path):
    rows = [ln.split() for ln in open(path) if ln.strip() and not ln.startswith("#")]
    return rows[0], np.array([[float(x) for x in r] for r in rows[1:]])


class _FlexFixture(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.d = tempfile.mkdtemp(prefix="sc_so_dyn_")
        cls.files = S._write_model(cls.d, so=True, soc=0.3)
        cls.flex = S._run_flex(cls.d, cls.files, "flex", nmat=_NMAT_E2E)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.d, True)

    def _out(self, name):
        out = os.path.join(self.d, name)
        self.addCleanup(shutil.rmtree, out, True)
        return out

    def _flex_input(self, out, flex=None, files=None, **eli):
        inp = S._flex_input(self.d, files or self.files, flex or self.flex,
                            _dyn_eli(**eli), out=out)
        inp["mode"]["param"]["Nmat"] = _NMAT_E2E
        return inp


class TestFlexMode(_FlexFixture):
    """``chi0q_mode = "flex"``: the stored chiq and chi0q at every bosonic
    frequency and the stored green of a spin-orbital FLEX run."""

    def _stored(self):
        with np.load(os.path.join(self.flex, "chiq.npz")) as z:
            chiq = z["chiq"]
        with np.load(os.path.join(self.flex, "chi0q.npz")) as z:
            chi0q = z["chi0q"]
        with np.load(os.path.join(self.flex, "green.npz")) as z:
            green = z["green"][0]
        nmat, nvol = chiq.shape[:2]
        nd2 = S.ND * S.ND
        Nx, Ny, Nz = S.CELL
        g = green.reshape(nmat, Nx, Ny, Nz, S.ND, S.ND).transpose(4, 5, 1, 2, 3, 0)
        return (chiq.reshape(nmat, nvol, nd2, nd2), chi0q.reshape(nmat, nvol, nd2, nd2),
                np.ascontiguousarray(g))

    def test_uses_stored_chi_at_every_frequency(self):
        inp = self._flex_input(self._out("o1"))
        cap = _run(inp)
        chi, chi0, g = self._stored()
        # the dressed chi is far from the bubble at finite frequency too
        m = _NMAT_E2E // 2 + 2
        self.assertGreater(np.linalg.norm(chi[m] - chi0[m]) / np.linalg.norm(chi0[m]), 1e-2)
        with _quiet():
            gam = sc._so_rpa_vertex(inp)
        A, basis = dso.build_operator(gam, chi, g, 4.0)
        dense = np.linalg.eigvals(_reduced_dense(A, basis))
        for lam in cap["vals"]:
            k = np.argmin(np.abs(dense - lam))
            self.assertLess(abs(dense[k] - lam), 1e-10 * np.abs(dense).max())
        # the frequency dependence of chi is load bearing: chi(nu = 0) at every
        # frequency gives a different leading eigenvalue
        flat = np.broadcast_to(chi[_NMAT_E2E // 2], chi.shape).copy()
        A0, _ = dso.build_operator(gam, flat, g, 4.0)
        d0 = np.linalg.eigvals(_reduced_dense(A0, basis))
        lead = cap["vals"][0]
        self.assertGreater(np.min(np.abs(d0 - lead)), 1e-4 * abs(lead))
        prov = open(os.path.join(inp["file"]["output"]["path_to_output"],
                                 "eigenvalue.dat")).read()
        self.assertIn("flex_vertex_backward_error_max", prov)

    def test_zero_frequency_slice_equals_static_reader(self):
        with _quiet():
            full = sc._load_flex_so(self._flex_input("x"), S.ND, *S.CELL, full=True)
            stat = sc._load_flex_so(self._flex_input("x"), S.ND, *S.CELL)
        m0 = _NMAT_E2E // 2
        np.testing.assert_array_equal(full[0][m0], stat[0])
        np.testing.assert_array_equal(full[2][m0], stat[2])
        np.testing.assert_array_equal(full[1], stat[1])

    def test_vertex_mismatch_is_refused(self):
        """Interaction files that differ from the FLEX run's are refused."""
        d = tempfile.mkdtemp(prefix="sc_so_dyn_mm_")
        self.addCleanup(shutil.rmtree, d, True)
        for f in os.listdir(self.d):
            if f.endswith(".dat"):
                shutil.copy(os.path.join(self.d, f), d)
        with open(os.path.join(d, "coulombinter.dat")) as f:
            txt = f.read()
        with open(os.path.join(d, "coulombinter.dat"), "w") as f:
            f.write(txt.replace("0.500000000000", "0.550000000000"))
        inp = S._flex_input(d, self.files, self.flex, _dyn_eli(), out="o")
        inp["mode"]["param"]["Nmat"] = _NMAT_E2E
        with self.assertRaisesRegex(ValueError, "do not reproduce the vertex"):
            with _quiet():
                sc.calc_eliashberg(inp)
        self.assertFalse(os.path.exists(os.path.join(d, "o")))

    def test_mismatch_at_one_finite_frequency_is_refused(self):
        """A chiq whose zero-frequency slice is right but one finite-frequency
        slice is not (the static check alone would pass): refused, naming
        that bosonic index."""
        d = tempfile.mkdtemp(prefix="sc_so_dyn_f_")
        self.addCleanup(shutil.rmtree, d, True)
        for f in os.listdir(self.flex):
            shutil.copy(os.path.join(self.flex, f), d)
        with np.load(os.path.join(self.flex, "chiq.npz")) as z:
            chiq = np.array(z["chiq"])
        m = _NMAT_E2E // 2 + 2
        chiq[m] = chiq[m + 1]
        S._resave(os.path.join(self.flex, "chiq.npz"), os.path.join(d, "chiq.npz"),
                  chiq=chiq)
        out = self._out("o_f")
        with self.assertRaisesRegex(ValueError, r"bosonic index {}\b".format(m)):
            with _quiet():
                sc.calc_eliashberg(self._flex_input(out, flex=d))
        self.assertFalse(os.path.exists(out) and os.listdir(out))
        # the static entry on the same archive does not see it
        with _quiet():
            sc.calc_eliashberg(S._flex_input(self.d, self.files, d,
                                             {"num_eigenvalues": 2},
                                             out=self._out("o_static")))

    def test_frequency_grid_mismatch_is_refused(self):
        out = self._out("o_g")
        inp = self._flex_input(out)
        inp["mode"]["param"]["Nmat"] = 2 * _NMAT_E2E
        with self.assertRaisesRegex(ValueError, "Nmat"):
            with _quiet(), mock.patch.object(
                    sc, "_load_flex_so", side_effect=AssertionError("loaded")):
                sc.calc_eliashberg(inp)
        self.assertFalse(os.path.exists(out) and os.listdir(out))


class TestCalcMode(_FlexFixture):
    """``chi0q_mode = "calc"``: the spinful bubble of the bare Green function
    at every bosonic frequency, resummed with Gamma."""

    def test_matches_an_independent_bubble(self):
        inp = S._input(self.d, self.files, True, _dyn_eli(chi0q_mode="calc"),
                       nmat=_NMAT_E2E, out=self._out("c1"))
        cap = _run(inp)
        g = cap["green"]
        with _quiet():
            gam = sc._so_rpa_vertex(inp)
        chi = _rpa(_bubble_w(g, 4.0), gam[0])
        A, basis = dso.build_operator(gam, chi, g, 4.0)
        dense = np.linalg.eigvals(_reduced_dense(A, basis))
        scale = np.abs(dense).max()
        for lam in cap["vals"]:
            self.assertLess(np.min(np.abs(dense - lam)), 1e-10 * scale)
        lead = dense[np.argmax(np.abs(dense))]
        self.assertLess(np.min(np.abs(cap["vals"] - lead)), 1e-10 * scale)

    def test_iteration_default_seed_is_random_over_frequency(self):
        out = self._out("c2")
        inp = S._input(self.d, self.files, True,
                       _dyn_eli(chi0q_mode="calc", solver_mode="iteration",
                                spectral_shift="auto", max_iter=3000,
                                convergence_tol=1e-10),
                       nmat=_NMAT_E2E, out=out)
        seen = {}
        real = sc._so_solve_iteration

        def spy(A, basis, sigma_init, *a, **kw):
            seen["init"] = np.array(sigma_init)
            seen["A"], seen["basis"] = A, basis
            return real(A, basis, sigma_init, *a, **kw)

        with _quiet(), mock.patch.object(sc, "_so_solve_iteration", side_effect=spy):
            lam = sc.calc_eliashberg(inp)
        w = dso.gap_weights(seen["init"], S.NP)
        self.assertGreater(min(w["parity"]), 0.1)        # both frequency parities
        self.assertLess(dso.antisymmetry_residual(seen["init"]), 1e-14)
        dense = np.linalg.eigvals(_reduced_dense(seen["A"], seen["basis"]))
        self.assertLess(abs(lam - dense.real.max()), 1e-6 * np.abs(dense).max())
        txt = open(os.path.join(out, "eigenvalue.dat")).read()
        self.assertIn("# init_gap = random", txt)
        self.assertIn("# eigenvalue_selection: iteration", txt)


class TestLMRetry(_FlexFixture):
    """which='LM' without a shift can miss the physical leading eigenvalue
    when large repulsive eigenvalues (which grow with Nmat) fill the
    num_eigenvalues largest-|lambda| set. On this model at Nmat = 32 the four
    largest |lambda| are all negative; as the non-spin-orbital dynamic solver
    does, the run re-solves once with spectral_shift = "auto" and reports the
    positive leading eigenvalue."""

    NMAT = 32

    def _run32(self, out):
        inp = S._input(self.d, self.files, True,
                       _dyn_eli(chi0q_mode="calc", num_eigenvalues=4),
                       nmat=self.NMAT, out=out)
        return _run(inp)

    def _check(self, cap, out, ref):
        self.assertGreater(ref, 0.0)
        self.assertGreater(cap["lam"], 0.0)
        self.assertLess(abs(cap["lam"] - ref), 1e-8 * abs(ref))
        txt = open(os.path.join(out, "eigenvalue.dat")).read()
        self.assertIn("# eigenvalue_selection: LR_retry", txt)
        self.assertIn("second solve with spectral_shift='auto'", txt)

    def test_retry_finds_the_positive_leading_eigenvalue(self):
        from scipy.sparse.linalg import eigs
        out = self._out("lm32")
        cap = self._run32(out)
        red = cap["basis"].reduce(cap["A"])
        # independent of the solver path: the 4 largest |lambda| are negative,
        # the largest real part is positive
        lm = eigs(red, k=4, which="LM", return_eigenvectors=False, tol=1e-12)
        self.assertTrue(np.all(lm.real < 0))
        lr = eigs(red, k=6, which="LR", return_eigenvectors=False, tol=1e-12)
        self._check(cap, out, float(lr.real.max()))

    @heavy
    def test_retry_matches_the_dense_reduced_spectrum(self):
        out = self._out("lm32d")
        cap = self._run32(out)
        dense = np.linalg.eigvals(_reduced_dense(cap["A"], cap["basis"]))
        lead = dense[np.argmax(dense.real)]
        # not among the num_eigenvalues = 4 largest in magnitude
        self.assertGreaterEqual(int(np.sum(np.abs(dense) > abs(lead) * (1 + 1e-9))), 4)
        self._check(cap, out, float(lead.real))


class TestOutputs(_FlexFixture):

    def test_schema_weights_and_antisymmetry(self):
        out = self._out("w")
        inp = self._flex_input(out, num_eigenvalues=8)
        gaps = []
        real = dso.gap_weights

        def spy(gap, n):
            gaps.append(np.array(gap))
            return real(gap, n)

        with mock.patch.object(dso, "gap_weights", side_effect=spy):
            cap = _run(inp)
        # every gap the solver produced (8 eigenvectors + the written one)
        self.assertGreaterEqual(len(gaps), 9)
        for g in gaps:
            R = dyn._reverse_kw_and_orbital(g)
            self.assertLess(np.linalg.norm(g + R), 1e-13 * np.linalg.norm(g))
        z = np.load(os.path.join(out, "gap_dynamic.npz"))
        gap = z["gap"]
        nd, (Nx, Ny, Nz) = S.ND, S.CELL
        self.assertEqual(gap.shape, (nd, nd, Nx, Ny, Nz, _NMAT_E2E))
        self.assertLess(np.linalg.norm(gap + dyn._reverse_kw_and_orbital(gap)),
                        1e-13 * np.linalg.norm(gap))
        self.assertAlmostEqual(float(np.linalg.norm(gap)), 1.0, places=12)
        for key in ("spin_orbital", "index_convention", "norb_phys", "spin_components",
                    "gap_spin_weights", "gap_spin_weights_w",
                    "gap_frequency_parity_weights", "gap_spin_frequency_weights",
                    "eigenvalues", "eigenvector_spin_weights",
                    "eigenvector_frequency_parity_weights", "eigenvalue_match",
                    "iomega", "momentum_convention", "pauli_leakage"):
            self.assertIn(key, z.files)
        for key in ("pairing_type", "eigenvectors", "gap_sector_weights"):
            self.assertNotIn(key, z.files)
        self.assertTrue(bool(z["spin_orbital"]))
        self.assertEqual(str(z["index_convention"]), "spin_block")
        self.assertEqual(z["gap_spin_weights_w"].shape, (_NMAT_E2E, 4))
        self.assertAlmostEqual(float(z["gap_spin_weights_w"].sum()), 1.0, places=12)
        np.testing.assert_allclose(z["gap_spin_weights_w"].sum(0), z["gap_spin_weights"],
                                   atol=1e-14)
        self.assertAlmostEqual(float(z["gap_frequency_parity_weights"].sum()), 1.0,
                               places=12)
        tab = z["gap_spin_frequency_weights"]
        self.assertEqual(tab.shape, (4, 2))
        np.testing.assert_allclose(tab.sum(1), z["gap_spin_weights"], atol=1e-12)
        np.testing.assert_allclose(tab.sum(0), z["gap_frequency_parity_weights"],
                                   atol=1e-12)
        self.assertLess(float(z["pauli_leakage"]), 1e-12)
        # eigenvalue.dat: tsweep's reader works; columns
        from hwave.tsweep import parse_leading_eig
        re_, im_, match = parse_leading_eig(os.path.join(out, "eigenvalue.dat"))
        self.assertEqual(match, 1)
        self.assertAlmostEqual(re_, float(np.real(cap["vals"][0])), places=7)
        txt = open(os.path.join(out, "eigenvalue.dat")).read()
        self.assertTrue(txt.startswith("# Dynamic Eliashberg leading eigenvalue\n"))
        self.assertIn("# spin_orbital = true", txt)
        self.assertIn("match(1=fermion-antisymmetric)  w_psi  w_dx  w_dy  w_dz  "
                      "w_even_w  w_odd_w", txt)
        lam_row, tabl = _read_table(os.path.join(out, "eigenvalue.dat"))
        self.assertEqual(tabl.shape, (8, 11))
        np.testing.assert_allclose(tabl[:, 5:9].sum(1), 1.0, atol=1e-7)
        np.testing.assert_allclose(tabl[:, 9:11].sum(1), 1.0, atol=1e-7)
        self.assertTrue(np.all(tabl[:, 4] == 1))
        np.testing.assert_allclose(tabl[0, 5:9], z["gap_spin_weights"], atol=1e-8)
        # gap.dat: the n0 = Nmat // 2 slice of the archive
        lines = open(os.path.join(out, "gap.dat")).read().splitlines()
        self.assertTrue(lines[0].startswith("# frequency=dynamic"))
        self.assertIn("# spin_orbital = true", lines)
        data = np.loadtxt(os.path.join(out, "gap.dat"))
        sl = gap[..., _NMAT_E2E // 2].reshape(nd * nd, Nx * Ny * Nz).T
        np.testing.assert_allclose(data[:, 3::2], sl.real, atol=1e-8)
        np.testing.assert_allclose(data[:, 4::2], sl.imag, atol=1e-8)

    def test_seed_continuation(self):
        """A spin-orbital gap_dynamic.npz seeds a spin-orbital run (the
        leading eigenpair is kept); it is refused by a non-spin-orbital run."""
        # largest real part in both runs: a seeded run is never re-solved
        # (as on the non-spin-orbital path), so an unshifted which='LM' set
        # need not contain the seeded state here
        out1 = self._out("s1")
        cap1 = _run(self._flex_input(out1, spectral_shift="auto"))
        seed = os.path.join(out1, "gap_dynamic.npz")
        cap2 = _run(self._flex_input(self._out("s2"), seed_eigenvector=seed,
                                     spectral_shift="auto"))
        self.assertAlmostEqual(complex(cap2["vals"][0]), complex(cap1["vals"][0]),
                               places=10)
        with self.assertRaisesRegex(ValueError, "spin_orbital"):
            dyn._load_seed_gap({"seed_eigenvector": seed},
                               (S.ND, S.ND) + tuple(S.CELL) + (_NMAT_E2E,), False,
                               None, _NMAT_E2E)

    def test_unmarked_seed_is_refused(self):
        out1 = self._out("s3")
        _run(self._flex_input(out1))
        bad = os.path.join(self._out("s4"), "seed.npz")
        os.makedirs(os.path.dirname(bad))
        S._resave(os.path.join(out1, "gap_dynamic.npz"), bad, drop=("spin_orbital",))
        out = self._out("s5")
        with self.assertRaisesRegex(ValueError, "spin_orbital = True marker"):
            with _quiet():
                sc.calc_eliashberg(self._flex_input(out, seed_eigenvector=bad))
        self.assertFalse(os.path.exists(out) and os.listdir(out))


class TestRefusals(unittest.TestCase):
    """Every unsupported combination is refused through the public entry,
    before any output, with the key named."""

    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="sc_so_dyn_r_")
        self.addCleanup(shutil.rmtree, self.d, True)

    def _refused(self, pattern, eli, files=None, exc=ValueError, **kw):
        files = files if files is not None else S._write_model(self.d, so=True, soc=0.3)
        base = {"frequency": "dynamic", "chi0q_mode": "calc",
                "solver_mode": "eigenvalue"}
        base.update(eli)
        base = {k: v for k, v in base.items() if v is not None}
        inp = S._input(self.d, files, True, base, nmat=_NMAT_E2E, **kw)
        with self.assertRaisesRegex(exc, pattern):
            with _quiet(), mock.patch.object(
                    sc, "_calc_chi0q_internal",
                    side_effect=AssertionError("chi0q computed before the check")):
                sc.calc_eliashberg(inp)
        out = inp["file"]["output"]["path_to_output"]
        self.assertFalse(os.path.exists(out) and os.listdir(out))

    def test_bond_channels(self):
        self._refused("bond_channels", {"bond_channels": True})

    def test_offsite_row(self):
        files = S._write_model(self.d, so=True, extra_rows={
            "CoulombInter": [(1, 0, 0, 0.2), (-1, 0, 0, 0.2)]})
        self._refused(r"off-site row .* in the CoulombInter file", {}, files=files)

    def test_cross_site_row_under_orbital_site(self):
        files = S._write_model(self.d, so=True, types=("CoulombIntra", "CoulombInter"))
        self._refused("orbital_site", {}, files=files, orbital_site=[0, 1])

    def test_ir_basis(self):
        self._refused("matsubara_basis", {"matsubara_basis": "ir"})
        self._refused("matsubara_basis", {"Matsubara_Basis": "IR"})

    def test_load_mode(self):
        self._refused("chi0q_mode", {"chi0q_mode": "load"})
        self._refused("chi0q_mode", {"chi0q_mode": None})         # default "load"
        self._refused("chi0q_mode", {"chi0q_mode": None, "Chi0q_Mode": "load"})

    def test_zero_chi(self):
        for key, val in (("zero_chi_s", True), ("zero_chi_c", "true"),
                         ("Zero_Chi_S", True), ("zero_chi_s", "ture")):
            with self.subTest(key=key, value=val):
                self._refused("(?i)" + key, {key: val})

    def test_zero_chi_false_is_accepted(self):
        """``false`` is not refused: the run proceeds past the configuration
        checks (stopped here at the chi0q computation)."""
        files = S._write_model(self.d, so=True, soc=0.3)
        for key in ("zero_chi_s", "Zero_Chi_C"):
            inp = S._input(self.d, files, True,
                           {"frequency": "dynamic", "chi0q_mode": "calc",
                            "solver_mode": "eigenvalue", key: False},
                           nmat=_NMAT_E2E)
            with self.subTest(key=key):
                with self.assertRaisesRegex(AssertionError, "reached"):
                    with _quiet(), mock.patch.object(
                            sc, "_calc_chi0q_internal",
                            side_effect=AssertionError("reached")):
                        sc.calc_eliashberg(inp)

    def test_pairing_type(self):
        self._refused("pairing_type", {"pairing_type": "singlet"})

    def test_eigenvalue_method(self):
        for m in ("subspace", "shift-invert-gmres"):
            with self.subTest(method=m):
                self._refused("eigenvalue_method", {"eigenvalue_method": m})

    def test_bond_path_keys(self):
        self._refused("parity_leakage_policy", {"parity_leakage_policy": "refuse"})

    def test_odd_nmat(self):
        self._refused("Nmat", {}, Nmat=15)

    def test_memory_cap(self):
        self._refused("mem_limit_gb", {"mem_limit_gb": 1e-6}, exc=MemoryError)

    def test_case_insensitive_frequency_dispatch(self):
        """``Frequency = "dynamic"`` reaches the spin-orbital dynamic solver
        (the configuration checks read it case-insensitively)."""
        self._refused("matsubara_basis", {"Frequency": "dynamic",
                                          "matsubara_basis": "ir"})

    def test_tsweep(self):
        import hwave.tsweep as ts
        base = {"mode": {"enable_spin_orbital": True,
                         "param": {"CellShape": [2, 2, 1], "Nmat": 8,
                                   "filling": 0.5}},
                "eliashberg": {"solver_mode": "eigenvalue", "frequency": "dynamic"},
                "file": {"input": {}, "output": {}}}
        with self.assertRaisesRegex(ValueError, "temperature sweep"):
            ts.preflight(base, {"run_eliashberg": True, "warm_start": False})

    def test_direct_solve_dynamic_call_is_checked_first(self):
        with self.assertRaisesRegex(ValueError, "chi0q_mode"):
            dyn.solve_dynamic({"mode": {"enable_spin_orbital": True},
                               "eliashberg": {"frequency": "dynamic"}})


class TestPairBubbleDtype(unittest.TestCase):
    """``calc_g2_dynamic`` scales in place only when that keeps the dtype of
    ``G2 / beta``; the result is bit-identical to ``G2 / beta`` for every
    input dtype."""

    def test_bit_identical_to_out_of_place_division(self):
        rng = np.random.default_rng(4)
        shape = (2, 2, 3, 1, 1, 4)
        beta = 3.7
        raw = rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
        for dt in (np.complex64, np.complex128, np.float32, np.float64,
                   np.int64, np.bool_):
            g = (raw.real > 0 if dt is np.bool_ else
                 (np.round(3 * raw.real) if dt is np.int64 else
                  (raw if np.dtype(dt).kind == "c" else raw.real))).astype(dt)
            with self.subTest(dtype=np.dtype(dt).name):
                gi = dyn.reverse_fft_axes(g[..., ::-1], (2, 3, 4))
                want = (g.reshape(4, 3, 4)[:, None] * gi.reshape(4, 3, 4)[None]
                        ).reshape(2, 2, 2, 2, 3, 1, 1, 4) / beta
                got = dyn.calc_g2_dynamic(g, beta)
                self.assertEqual(got.dtype, want.dtype)
                self.assertTrue(np.array_equal(got, want))


class TestMemoryEstimate(unittest.TestCase):

    def test_counts_chi_chi0_and_the_arpack_basis(self):
        nd, Nk, nmat = 4, 32 * 32, 512
        T_ = 16 * nd ** 4 * Nk * nmat
        need = dso.estimate_memory_bytes(nd, Nk, nmat, "eigenvalue", 10, "flex")
        self.assertGreaterEqual(need, 3 * T_)
        n = nd * nd * Nk * nmat
        self.assertGreater(need - dso.estimate_memory_bytes(
            nd, Nk, nmat, "iteration", 10, "flex"), 0)
        self.assertGreaterEqual(
            dso.estimate_memory_bytes(nd, Nk, nmat, "eigenvalue", 10, "flex"),
            2 * T_ + 16 * 25 * (n // 2))
        # larger than the non-spin-orbital dynamic estimate at the same nd
        self.assertGreater(need, dyn.estimate_memory_bytes(nd, Nk, nmat))
        with self.assertRaises(MemoryError):
            dso.check_memory(need, 1.0)
        dso.check_memory(need, 0)


@unittest.skipUnless(backend.gpu_available(), "needs CuPy with a CUDA device")
class TestGPU(unittest.TestCase):
    """The device operator equals the host one (matvec to 1e-12, leading
    eigenvalue), on the complex-SOC chain."""

    def test_host_equals_device(self):
        import cupy
        decl = _su2_decl()
        gam = T._gamma(decl)
        g = _green_chain(T._h_spin_orbital(0.3))
        chi = _rpa(_bubble_w(g, T.BETA), gam[0])
        A, basis = dso.build_operator(gam, chi, g, T.BETA)
        Ad, _ = dso.build_operator(gam, chi, g, T.BETA, xp=cupy)
        rng = np.random.default_rng(2)
        for _ in range(3):
            x = rng.standard_normal(basis.n) + 1j * rng.standard_normal(basis.n)
            a, b = A.matvec(x), Ad.matvec(x)
            self.assertIsInstance(b, np.ndarray)
            self.assertLess(np.linalg.norm(a - b), 1e-12 * np.linalg.norm(a))
        with _quiet():
            va, _ = sc._so_solve_eigenvalue(A, basis, 4, "arnoldi", map_back=False)
            vb, _ = sc._so_solve_eigenvalue(Ad, basis, 4, "arnoldi", map_back=False)
        self.assertLess(abs(va[0] - vb[0]), 1e-10 * abs(va[0]))

    def test_entry_gpu_equals_cpu(self):
        d = tempfile.mkdtemp(prefix="sc_so_dyn_gpu_")
        self.addCleanup(shutil.rmtree, d, True)
        files = S._write_model(d, so=True, soc=0.3)
        lam = {}
        for gpu in (False, True):
            inp = S._input(d, files, True,
                           _dyn_eli(chi0q_mode="calc", gpu=gpu, gpu_required=gpu),
                           nmat=_NMAT_E2E, out="o_{}".format(gpu))
            lam[gpu] = _run(inp)["vals"]
        np.testing.assert_allclose(lam[True], lam[False], rtol=1e-10)


if __name__ == "__main__":
    unittest.main()
