"""The second-order self-energy of the spin-orbital bond FLEX.

The production chain ``Gamma_B chibar_B Gamma_B -> so_bond.p_b ->
flex_bond.calc_self_energy_bond`` is compared, at second order in the
interaction, with an explicit Wick contraction of the documented operators
(``tests.so_bond_oracles``, which never imports the solver).

``TestConventions`` pins the conventions the comparison rests on:

* the imaginary-time grid ``tau_j = j beta / N`` of the oracle DFT
  matrices against ``hwave.solver.matsubara``, and the circular frequency
  sum of the tau route (the discretization control);
* the bridge between the k-space spin-block arrays and the interleaved
  exact-diagonalization modes (``modes_of_k`` / ``k_of_modes``) against the
  one-body table (``build_h1`` / ``band_k``);
* the equal-time sample ``G(tau = 0)`` and its density reading
  ``rho_pq = G_qp``;
* the absence of records whose operators 2 and 3 coincide.

``TestFirstOrderControl`` (C1a-C1) checks that the first-order functional
``Phi_1 = (N/2) sum_IJ Gamma_IJ(0) n_I n_J`` of the label vertex (direct
plus crossed filing) has the Wick Hartree-Fock self-energy as its gradient,
``Sigma_AB = beta dPhi_1/dG_BA``, with the one scalar ``beta``.
"""
import ast
import os
import unittest

import numpy as np

from hwave.solver import matsubara
from tests.so_bond_fixtures import (C1A_NMAT_FAST, C1A_ONSITE_TABLES, c1a_gamma_label,
                                    c1a_green, rel)
from tests.so_bond_oracles import (FIXTURE_I, FIXTURE_II, band_k, build_h1, build_terms,
                                   delta23_monomials, ed_mode, equal_time_density,
                                   equal_time_sample, expand_tables, geom_of, green_tau,
                                   k_of_modes, mesh_grid, mesh_momenta, minus_tau,
                                   modes_of_k, phi1_gradient, tau_to_matsubara,
                                   wick_hf_sigma_k)
from tests.test_so_bond_oracles import _imported_modules

_ORACLE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "so_bond_oracles.py")


def metric_t(nmat):
    """The C1a-T deviation on FIXTURE_I ``"rand"``: the tau route
    ``B [g_ab g_ab g~_ba]`` against the circular frequency sum
    ``T^2 sum_{n,n'} G_ab(n) G_ab(n') G_ba((n + n' - m) mod nmat)``, both in
    the mode space (max-norm relative). Asserts nothing."""
    fx = FIXTURE_I
    beta, shape, norb = fx.beta, fx.spatial_shape, fx.norb_phys
    G = c1a_green(fx, "rand", nmat)
    Gm = modes_of_k(G, shape, norb)                                     # (nmat, M, M)
    g = modes_of_k(green_tau(G, beta), shape, norb)
    tau_route = tau_to_matsubara(g * g * np.swapaxes(minus_tau(g), -1, -2), beta)
    T = 1.0 / beta
    GmT = np.swapaxes(Gm, -1, -2)
    n = np.arange(nmat)
    ref = np.zeros_like(Gm)
    for m in range(nmat):
        for n1 in range(nmat):
            ref[m] += (Gm[n1][None] * Gm * GmT[(n1 + n - m) % nmat]).sum(axis=0)
    return rel(tau_route, T * T * ref)


def fit_scalar(g, sigma):
    """The R6 one-scalar fit ``sigma = s g`` over the whole array:
    ``s = vdot(g, sigma) / vdot(g, g)`` and the residual
    ``max|s g - sigma| / max|sigma|``. Returns ``(residual, s)``."""
    s = np.vdot(g, sigma) / np.vdot(g, g)
    return float(np.abs(s * g - sigma).max() / np.abs(sigma).max()), complex(s)


def metric_c1(kind, nmat, **filing):
    """The C1a-C1 fit on FIXTURE_I: the Wick Hartree-Fock self-energy
    (``hf_onebody`` of the documented terms at the equal-time density, in k,
    the same at every frequency) against the gradient of ``Phi_1`` of the
    label vertex (``filing`` selects a wrong filing of ``Gamma``). Returns
    ``(residual, s)``; asserts nothing."""
    fx = FIXTURE_I
    beta, shape, norb = fx.beta, fx.spatial_shape, fx.norb_phys
    tables = fx.tables_dict
    G = c1a_green(fx, kind, nmat)
    rho = equal_time_density(G, beta, shape, norb)
    sigma = wick_hf_sigma_k(build_terms(tables, shape, norb), rho, shape, norb)
    sigma = np.broadcast_to(sigma[None], G.shape)
    geom = geom_of(fx)
    g = phi1_gradient(G, beta, c1a_gamma_label(fx, tables, geom, **filing), geom)
    return fit_scalar(g, sigma)


class TestConventions(unittest.TestCase):
    def test_oracle_imports_no_hwave(self):
        with open(_ORACLE_PATH, encoding="utf-8") as fh:
            names = _imported_modules(ast.parse(fh.read()))
        self.assertEqual({n for n in names if n.startswith("hwave")}, set())
        # non-vacuity: the scan sees a solver import
        seen = _imported_modules(ast.parse("from hwave.solver import matsubara\n"))
        self.assertTrue({n for n in seen if n.startswith("hwave")})

    def test_tau_dft_is_matsubara(self):
        fx = FIXTURE_I
        beta = fx.beta
        G = c1a_green(fx, "rand", C1A_NMAT_FAST)
        self.assertLessEqual(rel(green_tau(G, beta), matsubara.fermion_to_tau(G, axis=0) / beta),
                             1e-14)
        s = np.random.default_rng(2702).normal(size=G.shape) + 0j
        self.assertLessEqual(rel(tau_to_matsubara(s, beta), beta * matsubara.tau_to_fermion(s, axis=0)),
                             1e-14)

    def test_matsubara_round_trips(self):
        G = c1a_green(FIXTURE_I, "rand", C1A_NMAT_FAST)
        back = matsubara.tau_to_fermion(matsubara.fermion_to_tau(G, axis=0), axis=0)
        self.assertLessEqual(rel(back, G), 1e-14)
        rng = np.random.default_rng(2703)
        W = rng.normal(size=(8, 5, 4, 4)) + 1j * rng.normal(size=(8, 5, 4, 4))
        back = matsubara.tau_to_boson(matsubara.boson_to_tau(W, axis=0), axis=0)
        self.assertLessEqual(rel(back, W), 1e-14)

    def test_tau_route_is_the_circular_sum(self):
        dev = metric_t(C1A_NMAT_FAST)
        self.assertTrue(np.isfinite(dev), "non-finite metric")
        self.assertLessEqual(dev, 1e-13)

    def test_mode_bridge_is_band_k(self):
        for fx in (FIXTURE_I, FIXTURE_II):
            with self.subTest(fixture=fx.L):
                h1 = build_h1(fx.transfer_dict, fx.L, fx.norb_phys, zeeman=fx.zeeman)
                hk = band_k(fx.transfer_dict, fx.spatial_shape, fx.norb_phys, zeeman=fx.zeeman)
                self.assertLessEqual(rel(k_of_modes(h1, fx.spatial_shape, fx.norb_phys), hk),
                                     1e-14)
                self.assertLessEqual(rel(modes_of_k(hk, fx.spatial_shape, fx.norb_phys), h1),
                                     1e-14)
                rng = np.random.default_rng(2704)
                A = rng.normal(size=(3, fx.nvol, fx.nd, fx.nd)) \
                    + 1j * rng.normal(size=(3, fx.nvol, fx.nd, fx.nd))
                back = k_of_modes(modes_of_k(A, fx.spatial_shape, fx.norb_phys),
                                  fx.spatial_shape, fx.norb_phys)
                self.assertEqual(back.shape, A.shape)
                self.assertLessEqual(rel(back, A), 1e-14)
                M = fx.nvol * fx.nd
                bad = rng.normal(size=(M, M)) + 0j
                with self.assertRaisesRegex(ValueError,
                                            r"^k_of_modes: not translation invariant \(deviation "):
                    k_of_modes(bad, fx.spatial_shape, fx.norb_phys)

    def test_density_is_the_transposed_sample(self):
        for fx in (FIXTURE_I, FIXTURE_II):
            shape, norb, nd = fx.spatial_shape, fx.norb_phys, fx.nd
            R = mesh_grid(shape)
            phase = np.exp(1j * (R @ mesh_momenta(shape).T)) / fx.nvol      # (R, k)
            for kind in ("band", "rand", "herm"):
                with self.subTest(fixture=fx.L, kind=kind):
                    G = c1a_green(fx, kind, C1A_NMAT_FAST)
                    S = equal_time_sample(G, fx.beta)
                    np.testing.assert_array_equal(S, green_tau(G, fx.beta)[0])
                    rho = equal_time_density(G, fx.beta, shape, norb)
                    ref = np.einsum("rk,kyx->rxy", phase, S)                # [R, X, Y]
                    X = np.arange(nd)
                    rows = ed_mode(0, X, norb)[None, :, None]                  # (1, X, 1)
                    cols = ed_mode(np.arange(fx.nvol)[:, None], X[None, :], norb)[:, None, :]
                    got = rho[rows, cols]                                       # [R, X, Y]
                    self.assertLessEqual(rel(got, ref), 1e-14)
                    if kind != "rand":
                        self.assertLessEqual(rel(rho.conj().T, rho), 1e-14)

    def test_no_record_has_mode2_equal_mode3(self):
        cases = ((FIXTURE_I.tables_dict, FIXTURE_I.norb_phys, FIXTURE_I.spatial_shape),
                 (FIXTURE_II.tables_dict, FIXTURE_II.norb_phys, FIXTURE_II.spatial_shape),
                 (C1A_ONSITE_TABLES, 1, (5, 1, 1)))
        for tables, norb, shape in cases:
            with self.subTest(tables=sorted(tables)):
                mons = expand_tables(tables, norb)
                self.assertTrue(mons)
                self.assertEqual(delta23_monomials(mons, shape), [])
        # non-vacuity: on-site same-orbital Hund has mode 2 = mode 3
        hund = expand_tables({"Hund": {((0, 0, 0), (0, 0)): 1.0}}, 1)
        self.assertTrue(delta23_monomials(hund, (5, 1, 1)))


class TestFirstOrderControl(unittest.TestCase):
    def test_one_scalar_is_beta(self):
        beta = FIXTURE_I.beta
        for kind in ("band", "rand"):
            with self.subTest(kind=kind):
                residual, s = metric_c1(kind, C1A_NMAT_FAST)
                self.assertTrue(np.isfinite(residual) and np.isfinite(s), "non-finite metric")
                self.assertLessEqual(residual, 1e-12)
                self.assertLessEqual(abs(s - beta), 1e-12 * beta)

    def test_fock_is_carried_by_the_crossed_filing(self):
        for kind in ("band", "rand"):
            with self.subTest(kind=kind):
                residual, s = metric_c1(kind, C1A_NMAT_FAST, crossed=False)
                self.assertTrue(np.isfinite(residual) and np.isfinite(s), "non-finite metric")
                self.assertGreaterEqual(residual, 1e-3)


if __name__ == "__main__":
    unittest.main()
