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

``TestAdjoint`` checks the adjoint gradient of the second-order functional
``Phi_2 = -(T/8) sum_{q,l} Tr[(chibar(q, l) Gamma(-q))**2]``: against the
exhaustive circle stencil of every entry at Nmat 4, and by circle stencils
along random directions (Ruling R7) at Nmat 8.

``TestOnsiteControl`` (C1a-C2) checks, for the on-site CoulombIntra control
on FIXTURE_I's ring, that the explicit Wick second-order self-energy is the
gradient of ``Psi = sum_{q,l} Tr[(chibar Gamma(-q))**2]`` with the one
scalar ``-T/8``, and reports the production second order of one
spin-orbital FLEX map (``chi := chibar``) against the same Wick contraction.
"""
import ast
import os
import unittest
from unittest import mock

import numpy as np

from hwave.solver import matsubara
from tests.so_bond_fixtures import (C1A_NMAT_FAST, C1A_ONSITE_GEOM, C1A_ONSITE_TABLES,
                                    c1a_gamma_label, c1a_green, c1a_onsite_flex_sigma, rel, win)
from tests.so_bond_oracles import (FIXTURE_I, FIXTURE_II, antisym_tensor, band_k, build_h1,
                                   build_terms, delta23_monomials, ed_mode,
                                   equal_time_density, equal_time_sample, expand_tables,
                                   geom_of, green_tau, k_of_modes, mesh_grid, mesh_momenta,
                                   minus_tau, modes_of_k, phi1_gradient, phi2_gradient,
                                   random_green, ring2_gradient, ring2_value,
                                   tau_to_matsubara, wick_hf_sigma_k, wick_sigma2)
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


#: R7: stencil nodes, radii, seed of the first direction and the tolerance
R7_NODES = 8
R7_RADII = (1.0, 0.25)
R7_SEED0 = 2710
R7_TOL = 1e-12
#: R7 non-vacuity floor: a direction with ``|a| < R7_FLOOR max|Sigma|`` is redrawn
R7_FLOOR = 1e-6


def directional_rows(value_fn, sigma, green, beta, n_dir, seed0=R7_SEED0):
    """The R7 directional rows of a gradient ``sigma[n, k, A, B] = beta
    dPhi/dG_BA(k, n)`` of the functional ``value_fn`` (a stack ``(S, nmat,
    nvol, nd, nd)`` of Green functions -> ``(S,)`` values of ``Phi``).

    Direction ``i`` is ``dG = rng.normal + i rng.normal`` of
    ``np.random.default_rng(seed)``, scaled to ``max|dG| = max|G|``; the
    first ``n_dir`` seeds are ``seed0 + i``. For each direction ``a(s) =
    beta (1/(M s)) sum_j Phi(G + s w^j dG) w^{-j}`` (``w = e^{2 pi i/M}``,
    ``M = R7_NODES``) at the radii ``s = 1`` and ``s = 1/4``, ``a = a(1)``,
    and ``b = sum sigma[n, k, A, B] dG[n, k, B, A]`` (bilinear). A direction
    with ``|a| < R7_FLOOR max|sigma|`` is not accepted and is replaced by the
    next unused seed ``seed0 + n_dir + j`` (``j = 0, 1, ...``), with at most
    ``4 n_dir`` draws in all.

    Each row records ``seed``, ``accepted``, ``a``, ``a_quarter``, ``b``,
    ``r_ab = |a - b| / max(max|sigma|, |a|)``, ``r_s = |a(1) - a(1/4)| /
    max(|a(1)|, max|sigma|)``, and the failure report of R7: ``phi`` =
    ``|Phi(G)|``, ``max_g``, ``max_dg``, ``s`` = the radii and ``fourier`` =
    the ``M`` coefficients ``(1/M) sum_j Phi(G + s w^j dG) w^{-jp}`` of each
    radius. Asserts nothing."""
    green = np.asarray(green, dtype=np.complex128)
    sigma = np.asarray(sigma, dtype=np.complex128)
    M = R7_NODES
    w = np.exp(2j * np.pi * np.arange(M) / M)
    smax = float(np.abs(sigma).max())
    gmax = float(np.abs(green).max())
    phi0 = abs(complex(np.asarray(value_fn(green[None])).reshape(-1)[0]))
    radii = np.array(R7_RADII)
    nodes = (radii[:, None] * w[None, :]).reshape(-1)                   # (len(radii) M,)
    rows, accepted = [], 0
    seeds = list(range(seed0, seed0 + n_dir))
    spare = seed0 + n_dir
    while seeds and len(rows) < 4 * n_dir:
        seed = seeds.pop(0)
        rng = np.random.default_rng(seed)
        dg = rng.normal(size=green.shape) + 1j * rng.normal(size=green.shape)
        dg *= gmax / float(np.abs(dg).max())
        vals = np.asarray(value_fn(green[None] + nodes[:, None, None, None, None] * dg[None]))
        vals = vals.reshape(len(radii), M)
        fourier = np.array([[np.sum(v * w ** (-p)) / M for p in range(M)] for v in vals])
        a_s = beta * fourier[:, 1] / radii
        a = complex(a_s[0])
        b = complex(np.sum(sigma * np.swapaxes(dg, -1, -2)))
        ok = abs(a) >= R7_FLOOR * smax
        rows.append({"seed": seed, "accepted": bool(ok), "a": a, "a_quarter": complex(a_s[1]),
                     "b": b, "r_ab": abs(a - b) / max(smax, abs(a)),
                     "r_s": abs(a_s[0] - a_s[1]) / max(abs(a_s[0]), smax),
                     "phi": phi0, "max_g": gmax, "max_dg": float(np.abs(dg).max()),
                     "s": tuple(R7_RADII), "fourier": fourier})
        if ok:
            accepted += 1
        else:
            seeds.append(spare)
            spare += 1
    return rows


def metric_adjoint(fx, kind, nmat, n_dir):
    """The R7 check of ``Sigma_F = -ring2_gradient/8`` (``beta dPhi_2/dG``)
    against ``Phi_2 = -(T/8) ring2_value`` on ``fx`` (its tables on its label
    geometry, the Green function ``c1a_green(fx, kind, nmat)``).

    Returns ``(worst, rows)``: ``rows`` from :func:`directional_rows` and
    ``worst`` the largest ``r_ab`` / ``r_s`` of the accepted rows; when fewer
    than ``n_dir`` directions are accepted the shortfall counts as a relative
    miss of 1. R7 passes iff ``worst <= R7_TOL``. Asserts nothing."""
    beta = fx.beta
    geom = geom_of(fx)
    gamma = c1a_gamma_label(fx, fx.tables_dict, geom)
    G = c1a_green(fx, kind, nmat)
    sigma = -ring2_gradient(G, beta, gamma, geom) / 8.0
    rows = directional_rows(lambda st: -(1.0 / beta) / 8.0 * ring2_value(st, beta, gamma, geom),
                            sigma, G, beta, n_dir)
    good = [r for r in rows if r["accepted"]]
    worst = max([max(r["r_ab"], r["r_s"]) for r in good], default=0.0)
    if len(good) < n_dir:
        worst = max(worst, 1.0)
    return float(worst), rows


def _onsite_wick(kind, nmat):
    """``(G, sigma)``: FIXTURE_I's Green function and the Wick second order
    of the on-site control tables on the 5-ring."""
    fx = FIXTURE_I
    G = c1a_green(fx, kind, nmat)
    terms = build_terms(C1A_ONSITE_TABLES, fx.spatial_shape, fx.norb_phys)
    return G, wick_sigma2(G, fx.beta, terms, fx.spatial_shape, fx.norb_phys)


def metric_c2(kind, nmat):
    """The C1a-C2 fit on FIXTURE_I with the on-site control tables: the Wick
    second-order self-energy against ``beta ring2_gradient`` of the B = 1
    label vertex (expected scalar ``-T/8``). Returns ``(residual, s)``;
    asserts nothing."""
    fx = FIXTURE_I
    G, sigma = _onsite_wick(kind, nmat)
    gamma = c1a_gamma_label(fx, C1A_ONSITE_TABLES, C1A_ONSITE_GEOM)
    g = fx.beta * ring2_gradient(G, fx.beta, gamma, C1A_ONSITE_GEOM)
    return fit_scalar(g, sigma)


def prod_dev(kind, nmat, *, rpa_return):
    """The production second order of one spin-orbital FLEX map with the
    on-site control tables (``rpa_return`` the hook in place of the RPA
    solve) against the Wick second order, max-norm relative. Asserts
    nothing."""
    G, sigma = _onsite_wick(kind, nmat)
    return rel(c1a_onsite_flex_sigma(G, rpa_return=rpa_return), sigma)


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

    def test_wick_sigma2_hubbard_atom(self):
        U = 0.83
        beta = FIXTURE_I.beta
        tables = {"CoulombIntra": {((0, 0, 0), (0, 0)): U}}
        shape = (1, 1, 1)
        terms = build_terms(tables, shape, 1)
        # (a) spin-diagonal random G: Sigma_up = FT[-U^2 g_up g_dn G~_dn], and the mirror
        rng = np.random.default_rng(2705)
        nmat = C1A_NMAT_FAST
        G = np.zeros((nmat, 1, 2, 2), dtype=np.complex128)
        for s in (0, 1):
            G[:, 0, s, s] = rng.normal(size=nmat) + 1j * rng.normal(size=nmat)
        sig = wick_sigma2(G, beta, terms, shape, 1)
        self.assertEqual(sig.shape, G.shape)
        g = green_tau(G, beta)[:, 0]
        gt = minus_tau(g)
        for s, o in ((0, 1), (1, 0)):
            with self.subTest(spin=s):
                ref = tau_to_matsubara(-U * U * g[:, s, s] * g[:, o, o] * gt[:, o, o], beta)
                self.assertLessEqual(rel(sig[:, 0, s, s], ref), 1e-14)
        # (b) analytic half-filled atom G = 1/(i w): Sigma = U^2/(4 i w), first order in 1/N
        errs = {}
        for N in (512, 1024):
            w = (2 * np.arange(N) + 1 - N) * np.pi / beta
            Ga = np.zeros((N, 1, 2, 2), dtype=np.complex128)
            Ga[:, 0, 0, 0] = Ga[:, 0, 1, 1] = 1.0 / (1j * w)
            sa = wick_sigma2(Ga, beta, terms, shape, 1)
            ref = (U * U / (4j * w))[:, None, None] * np.eye(2)
            errs[N] = rel(sa[win(N), 0], ref[win(N)])
            self.assertTrue(np.isfinite(errs[N]), "non-finite metric")
            self.assertLessEqual(errs[N], 10.0 / N, errs)
        ratio = errs[512] / errs[1024]
        self.assertTrue(1.5 <= ratio <= 2.5, (ratio, errs))
        # (c) spin-mixing random G: the spin-off-diagonal entries are populated
        Gm = random_green(2706, nmat, 1, 2)
        sm = wick_sigma2(Gm, beta, terms, shape, 1)
        smax = float(np.abs(sm).max())
        for (a, b) in ((0, 1), (1, 0)):
            self.assertGreater(float(np.abs(sm[..., a, b]).max()), 1e-3 * smax)

    def test_antisym_tensor_refuses_mode2_equal_mode3(self):
        U = antisym_tensor([(0, 0, 1, 1, 0.5)], 2)
        H = {(0, 1, 0, 1): 0.5, (1, 0, 0, 1): -0.5, (0, 1, 1, 0): -0.5, (1, 0, 1, 0): 0.5}
        for idx in np.ndindex(U.shape):
            self.assertEqual(U[idx], H.get(idx, 0.0))
        with self.assertRaisesRegex(ValueError, r"^antisym_tensor: term \(0, 1, 1, 0, 0\.5\) "
                                                r"has mode 2 = mode 3$"):
            antisym_tensor([(0, 1, 1, 0, 0.5)], 2)

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


class TestOnsiteControl(unittest.TestCase):
    def test_one_scalar_is_minus_t_over_8(self):
        T = 1.0 / FIXTURE_I.beta
        for kind in ("band", "rand"):
            with self.subTest(kind=kind):
                residual, s = metric_c2(kind, C1A_NMAT_FAST)
                self.assertTrue(np.isfinite(residual) and np.isfinite(s), "non-finite metric")
                self.assertLessEqual(residual, 1e-12)
                self.assertLessEqual(abs(s + T / 8.0), 1e-12 * T / 8.0)

    def test_production_onsite_second_order(self):
        # reported (Ruling R12): not part of the control verdict
        for kind in ("band", "rand"):
            with self.subTest(kind=kind):
                spy = mock.MagicMock(side_effect=lambda c: c)
                dev = prod_dev(kind, C1A_NMAT_FAST, rpa_return=spy)
                self.assertEqual(spy.call_count, 1)
                self.assertTrue(np.isfinite(dev), "non-finite metric")
                self.assertLessEqual(dev, 1e-12)


class TestAdjoint(unittest.TestCase):
    def test_exhaustive_stencil_nmat4(self):
        for fx in (FIXTURE_I, FIXTURE_II):
            with self.subTest(fixture=fx.L):
                geom = geom_of(fx)
                gamma = c1a_gamma_label(fx, fx.tables_dict, geom)
                G = c1a_green(fx, "rand", 4)
                adj = -ring2_gradient(G, fx.beta, gamma, geom) / 8.0
                ref = phi2_gradient(G, fx.beta, gamma, geom)
                self.assertEqual(adj.shape, G.shape)
                dev = rel(adj, ref)
                self.assertTrue(np.isfinite(dev), "non-finite metric")
                self.assertLessEqual(dev, 1e-12)

    def test_directional_fast(self):
        for fx in (FIXTURE_I, FIXTURE_II):
            for kind in ("band", "rand"):
                with self.subTest(fixture=fx.L, kind=kind):
                    worst, rows = metric_adjoint(fx, kind, C1A_NMAT_FAST, 8)
                    self.assertTrue(np.isfinite(worst), "non-finite metric")
                    self.assertEqual(sum(r["accepted"] for r in rows), 8)
                    self.assertLessEqual(worst, R7_TOL, rows)


if __name__ == "__main__":
    unittest.main()
