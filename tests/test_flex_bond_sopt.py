"""Gate G2 (spec 2026-09-06 section 5): perturbative ownership of the
bond-resolved FLEX vertex. With a FROZEN bare Green function (fixed mu, no
SCF) the fluctuation self-energy of one map is fitted to a quadratic in
(U, V) on the 3 x 3 grid {0, u0/2, u0} x {0, v0/2, v0}; (a) its constant
and linear parts vanish and the fit residual scales as the cube of the
grid; (b) the Hartree-Fock map is linear in every coefficient and the
PairLift map is identically zero on the paramagnetic density; (c) the
quadratic coefficients equal an INDEPENDENT hand-enumerated second-order
evaluation (the two skeleton diagrams, direct and exchange, enumerated in
real space on the torus with explicit spin bookkeeping -- no production
vertex code); (d) with every off-site coefficient zero the gate reproduces
the general path's coefficients."""
import os
import shutil
import tempfile
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k

_IN1 = "tests/rpa/input"
_IN2 = "tests/rpa/input_2orb"
_NMAT = 8
_SHAPE = (4, 4, 1)
_T = 1.0
_MU = 0.1


def _write_wan(path, name, norb, rows):
    rvecs = sorted({tuple(r[:3]) for r in rows})
    with open(path, "w") as fw:
        fw.write("{} in wannier90-like format for uhfk\n{}\n{}\n".format(name, norb, len(rvecs)))
        fw.write(" ".join("1" for _ in rvecs) + "\n")
        for r in rows:
            fw.write("{:4d} {:4d} {:4d} {:4d} {:4d} {: .15e} {: .15e}\n".format(*r))


def _make_inputs(d, norb, U, V, uprime=0.0, pairlift=None):
    src = _IN1 if norb == 1 else _IN2
    for f in ("geom.dat", "transfer.dat"):
        shutil.copy(os.path.join(src, f), d)
    _write_wan(os.path.join(d, "intra.dat"), "CoulombIntra", norb,
               [(0, 0, 0, a, a, U, 0.0) for a in range(1, norb + 1)])
    rows = []
    if norb == 2:
        rows += [(0, 0, 0, 1, 2, uprime, 0.0), (0, 0, 0, 2, 1, uprime, 0.0)]
        rows += [(1, 0, 0, 1, 1, V, 0.0), (-1, 0, 0, 1, 1, V, 0.0),
                 (0, 1, 0, 2, 2, V, 0.0), (0, -1, 0, 2, 2, V, 0.0),
                 (1, 0, 0, 1, 2, 0.6 * V, 0.0), (-1, 0, 0, 2, 1, 0.6 * V, 0.0)]
    else:
        rows += [(1, 0, 0, 1, 1, V, 0.0), (-1, 0, 0, 1, 1, V, 0.0),
                 (0, 1, 0, 1, 1, 0.5 * V, 0.0), (0, -1, 0, 1, 1, 0.5 * V, 0.0)]
    _write_wan(os.path.join(d, "inter.dat"), "CoulombInter", norb, rows)
    inter = {"CoulombIntra": "intra.dat", "CoulombInter": "inter.dat"}
    if pairlift is not None:
        _write_wan(os.path.join(d, "pairlift.dat"), "PairLift", norb, pairlift)
        inter["PairLift"] = "pairlift.dat"
    return inter


def _solver(d, inter, gate):
    import hwave.solver.flex as flex_mod
    idict = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
    idict.update(inter)
    r = read_input_k.QLMSkInput({"path_to_input": d, "interaction": idict})
    par = {"T": _T, "mu": _MU, "CellShape": list(_SHAPE), "SubShape": [1, 1, 1], "Nmat": _NMAT,
           "IterationMax": 1, "Mix": 1.0, "EPS": 1e-12, "flex_hartree_fock": True,
           "longitudinal_bond_channels": gate}
    info = {"mode": "FLEX", "param": par, "enable_spin_orbital": False, "calc_scheme": "general"}
    return flex_mod.FLEX(r.get_param("ham"), {}, info), r.get_param("green")


def _one_map(s, gi, gate):
    """(Sigma_HF (nvol,norb,norb), Sigma_fluct (nmat,nvol,norb,norb), G) of ONE
    map from the frozen bare G at the fixed mu."""
    from hwave.solver import flex_bond
    s._phase_b_reset(gi)
    s._phase_b_preflight(gi)
    s._calc_epsilon_k({})
    beta = 1.0 / s.T
    mu = s.mu_value
    green0, green0_tail = s._calc_green(beta, mu)
    nmat, nvol, norb = s.nmat, s.lattice.nvol, s.norb
    nd = norb * norb
    G = s._calc_dressed_green(beta, mu, np.zeros((1, nmat, nvol, norb, norb), complex))
    static0 = np.zeros((1, 1, nvol, norb, norb), complex)
    _, sigma_hf, _ = s._phase_b_density_and_hf(G, static0, mu, beta, None)
    if gate:
        B = s._bond_view.n_channels
        with flex_bond.BondBlockStore(nmat, nvol, B * nd, nd, ("chibar", "W")) as store:
            s._phase_b_prepare_vertices()
            flex_bond.assemble_bubble(store, G, green0_tail, beta, s._bond_view, _SHAPE, 1)
            flex_bond.dress_and_build_w(store, s._bond_S, s._bond_C, S_on=s._bond_S_on,
                                        C_on=s._bond_C_on, nb=nmat, output_full=False,
                                        nmat=nmat, nvol=nvol, nd=nd, spatial_shape=_SHAPE)
            sig = flex_bond.calc_self_energy_bond(store, G, beta, s._bond_view, _SHAPE, norb, 1)
    else:
        chi0q_raw = s._calc_chi0q(G, green0_tail, beta)[0]
        _, v_eff, _, _ = s._flex_compute_veff_general(chi0q_raw, s.ham_info.ham_inter_q)
        sig = s._calc_self_energy_general(G, v_eff, beta)
    return sigma_hf[0, 0], sig[0], G


# ---------------------------------------------------------------------------
# the independent second-order oracle (real-space enumeration on the torus)
# ---------------------------------------------------------------------------

def _site(x, y):
    nx, ny, _ = _SHAPE
    return (x % nx) * ny + (y % ny)


def _coords(site):
    nx, ny, _ = _SHAPE
    return site // ny, site % ny


def _interaction_matrices(norb, U, V, uprime, V_rows):
    """W^{up up} and W^{up down} on the site-orbital torus (M x M), from the
    same declarations the solver reads: H = sum over DECLARED rows of
    v n_i n_j (each ordered pair once); the symmetric matrix W_ij with
    W_ij = W_ji = v then enters 1/2 sum_ij W_ij n_i n_j with W = 2 v for a
    pair declared in both orders (the reader's Hermitian closure)."""
    nvol = _SHAPE[0] * _SHAPE[1]
    M = nvol * norb
    Wuu = np.zeros((M, M))
    Wud = np.zeros((M, M))
    for s in range(nvol):
        x, y = _coords(s)
        for a in range(norb):
            i = s * norb + a
            Wud[i, i] += U                      # Hubbard: opposite spins only
            for (rx, ry, rz, a1, b1, v, _im) in V_rows:
                if a1 - 1 != a:
                    continue
                j = _site(x + rx, y + ry) * norb + (b1 - 1)
                Wuu[i, j] += v
                Wud[i, j] += v
    return Wuu, Wud


def _sopt_oracle(G, beta, norb, Wuu, Wud):
    """Direct and exchange second-order skeletons with the frozen G, in the
    same discrete (r, tau) representation the solver's transforms realise
    (fermionic tau grid, FFT-grid reversal for G(-r, -tau) with the
    antiperiodic sign, spatial FFTs); returns (Sigma_direct,
    Sigma_exchange) as (nmat, nvol, norb, norb) on (k, i w_n)."""
    from hwave.solver import matsubara as _ms, backend as _bk
    nx, ny, nz = _SHAPE
    nvol = nx * ny * nz
    nmat = G.shape[1]
    M = nvol * norb
    g_rt = _bk.spatial_ifftn(
        _ms.fermion_to_tau(G[0].reshape(nmat, nvol * norb * norb), axis=0).reshape(nmat, nx, ny, nz, norb * norb),
        axes=(1, 2, 3), workers=1).reshape(nmat, nx, ny, nz, norb, norb)
    rev = np.zeros_like(g_rt)
    for t in range(nmat):
        for ix in range(nx):
            for iy in range(ny):
                rev[t, ix, iy, 0] = g_rt[(-t) % nmat, (-ix) % nx, (-iy) % ny, 0]
    g_rt = g_rt.reshape(nmat, nvol, norb, norb)
    rev = rev.reshape(nmat, nvol, norb, norb)
    sgn = -np.ones(nmat)
    sgn[0] = 1.0
    # site-orbital matrices
    Gf = np.zeros((nmat, M, M), complex)
    Gr = np.zeros((nmat, M, M), complex)
    for si in range(nvol):
        xi, yi = _coords(si)
        for sj in range(nvol):
            xj, yj = _coords(sj)
            d = _site(xi - xj, yi - yj)
            for a in range(norb):
                for b in range(norb):
                    Gf[:, si * norb + a, sj * norb + b] = g_rt[:, d, a, b]
                    Gr[:, si * norb + a, sj * norb + b] = rev[:, d, a, b]
    L = Gf * Gr.transpose(0, 2, 1)                          # L_kl = G_kl(tau) G_lk(-tau)
    K = Wuu @ L @ Wuu + Wud @ L @ Wud.T
    sig_d = -(1.0 / beta) * sgn[:, None, None] * Gf * K
    sig_x = (1.0 / beta) * sgn[:, None, None] * np.einsum('ik,tkj,til,lj,tlk->tij', Wuu, Gf, Gf, Wuu, Gr,
                                                          optimize=True)

    def _to_kw(sig):
        out = np.zeros((nmat, nvol, norb, norb), complex)
        for s in range(nvol):
            for a in range(norb):
                for b in range(norb):
                    out[:, s, a, b] = sig[:, s * norb + a, b]      # j = (site 0, b)
        kt = _bk.spatial_fftn(out.reshape(nmat, nx, ny, nz, norb * norb), axes=(1, 2, 3), workers=1)
        return _ms.tau_to_fermion(kt.reshape(nmat, nvol * norb * norb), axis=0).reshape(
            nmat, nvol, norb, norb) / beta
    return _to_kw(sig_d), _to_kw(sig_x)


def _fit(points, values):
    """Unweighted least squares of c20 U^2 + c11 U V + c02 V^2 + c10 U + c01 V + c00
    per element; returns (coeffs dict, residual (rss over the points))."""
    A = np.array([[u * u, u * v, v * v, u, v, 1.0] for (u, v) in points])
    Y = np.array([np.asarray(val).ravel() for val in values])
    C, _, _, _ = np.linalg.lstsq(A, Y, rcond=None)
    fit = A @ C
    rss = np.sqrt(np.sum(np.abs(Y - fit) ** 2))
    shape = np.asarray(values[0]).shape
    names = ("c20", "c11", "c02", "c10", "c01", "c00")
    return {n: C[i].reshape(shape) for i, n in enumerate(names)}, rss


def _rel(a, b):
    return np.linalg.norm(np.asarray(a).ravel() - np.asarray(b).ravel()) / max(
        np.linalg.norm(np.asarray(b).ravel()), 1e-300)


class _Grid:
    """One map per grid point on a temporary input directory."""

    def __init__(self, norb, gate, u0=0.5, v0=0.5, scale=1.0, uprime=0.0):
        self.norb, self.gate = norb, gate
        self.points = [(u * scale, v * scale) for u in (0.0, u0 / 2, u0) for v in (0.0, v0 / 2, v0)]
        self.hf, self.fluct = [], []
        with tempfile.TemporaryDirectory() as d:
            for (U, V) in self.points:
                inter = _make_inputs(d, norb, U, V, uprime=uprime * (U / u0 if u0 else 0.0))
                s, gi = _solver(d, inter, gate)
                hf, fl, G = _one_map(s, gi, gate)
                self.hf.append(hf)
                self.fluct.append(fl)
                self.G = G
                self.rows = _read_rows(os.path.join(d, "inter.dat"))
        self.coeffs, self.rss = _fit(self.points, self.fluct)


def _read_rows(path):
    rows = []
    for ln in open(path).read().splitlines()[4:]:
        p = ln.split()
        rows.append((int(p[0]), int(p[1]), int(p[2]), int(p[3]), int(p[4]), float(p[5]), float(p[6])))
    return rows


_H = 2.0e-3          # coupling scale of the coefficient-extraction grids
_U0 = _V0 = 0.5      # the spec's grid; used for the residual-scaling test


def _richardson(cs):
    """Three-level extrapolation of a fitted coefficient c(h) = a + b1 h +
    b2 h^2 + O(h^3) from the grids at h, h/2, h/4 (the 3 x 3 quadratic fit
    interpolates, so the cubic and quartic terms of Sigma alias into the
    fitted coefficients at O(h) and O(h^2))."""
    return (8.0 * cs[2] - 6.0 * cs[1] + cs[0]) / 3.0


def _coefficients(norb, gate, uprime=0.0):
    grids = [_Grid(norb, gate=gate, scale=_H / _U0 * f, uprime=uprime) for f in (1.0, 0.5, 0.25)]
    out = {k: _richardson([g.coeffs[k] for g in grids]) for k in ("c20", "c11", "c02")}
    return out, grids[0]


class TestG2(unittest.TestCase):
    """Coefficient extraction at h = 2e-3 with three-level Richardson
    extrapolation (a3/a2 is 0.1-1 on this fixture, so the spec's u0 = v0 =
    0.5 grid aliases the cubic terms at 1e-2 into every fitted coefficient;
    the residual-scaling test keeps the spec's grid)."""

    @classmethod
    def setUpClass(cls):
        cls.on, cls.g_on = _coefficients(1, gate=True)
        cls.off, cls.g_off = _coefficients(1, gate=False)
        beta = 1.0 / _T
        G = cls.g_on.G
        rows = cls.g_on.rows
        vscale = _V0 * _H / _U0            # the V of the largest grid point
        Wuu_v, Wud_v = _interaction_matrices(1, 0.0, 0.0, 0.0,
                                             [r[:5] + (r[5] / vscale, r[6]) for r in rows])
        Wuu_u, Wud_u = _interaction_matrices(1, 1.0, 0.0, 0.0, [])
        d_u, x_u = _sopt_oracle(G, beta, 1, Wuu_u, Wud_u)
        d_v, x_v = _sopt_oracle(G, beta, 1, Wuu_v, Wud_v)
        d_uv, x_uv = _sopt_oracle(G, beta, 1, Wuu_u + Wuu_v, Wud_u + Wud_v)
        cls.oracle = dict(c20=d_u + x_u, c02=d_v + x_v, c11=(d_uv + x_uv) - (d_u + x_u) - (d_v + x_v),
                          d_v=d_v, x_v=x_v, x_u=x_u, x_uv=x_uv, x_v_only=x_v)

    def test_a_constant_and_linear_parts_vanish(self):
        c = self.g_on.coeffs
        quad = max(np.abs(c[k]).max() for k in ("c20", "c11", "c02"))
        low = max(np.abs(c[k]).max() for k in ("c10", "c01", "c00"))
        self.assertLess(low, 1e-6 * quad)

    def test_a_residual_scales_as_the_cube_of_the_grid(self):
        grids = [_Grid(1, gate=True, scale=sc) for sc in (1.0, 0.5, 0.25)]
        ref = np.linalg.norm(np.asarray(grids[0].fluct).ravel())
        floor = 1e3 * np.finfo(float).eps * ref
        pts = [(sc, g.rss) for sc, g in zip((1.0, 0.5, 0.25), grids) if g.rss >= floor]
        if len(pts) < 2:
            self.skipTest("fit residual at the round-off floor; exponent passes vacuously")
        sv = np.array([x[0] for x in pts]); rv = np.array([x[1] for x in pts])
        slope = np.polyfit(np.log(sv), np.log(rv), 1)[0]
        self.assertTrue(2.5 <= slope <= 3.5, "slope {}".format(slope))

    def test_b_hartree_fock_map_is_linear(self):
        g = self.g_on
        hf = {p: h for p, h in zip(g.points, g.hf)}
        u0, v0 = max(p[0] for p in g.points), max(p[1] for p in g.points)
        h00 = hf[(0.0, 0.0)]
        du = (hf[(u0, 0.0)] - h00) / u0
        dv = (hf[(0.0, v0)] - h00) / v0
        ref = max(np.linalg.norm(du.ravel()), np.linalg.norm(dv.ravel()))
        for (U, V), h in hf.items():
            pred = h00 + U * du + V * dv
            self.assertLess(np.linalg.norm((h - pred).ravel()), 1e-12 * ref)
        self.assertGreater(np.linalg.norm(dv.ravel()), 1e-3)     # V has a nonzero (Fock) map

    def test_b_pairlift_map_is_zero_on_the_paramagnetic_density(self):
        from hwave.solver import flex_hf
        with tempfile.TemporaryDirectory() as d:
            rows = [(1, 0, 0, 1, 2, 0.3, 0.0), (-1, 0, 0, 2, 1, 0.3, 0.0),
                    (0, 0, 0, 1, 2, 0.2, 0.0), (0, 0, 0, 2, 1, 0.2, 0.0)]
            inter = _make_inputs(d, 2, 0.4, 0.3, uprime=0.2, pairlift=rows)
            s, gi = _solver(d, inter, gate=True)
            hf, _, G = _one_map(s, gi, gate=True)
            self.assertGreater(np.abs(hf).max(), 1e-3)
            from hwave.solver import hartree_fock as _hf
            s._calc_epsilon_k({})
            heff = _hf.heff_eigenpairs(np.asarray(s.H0_k)[0], np.zeros((s.lattice.nvol, 2, 2), complex))
            dens = _hf.equal_time_density(G, heff, s.mu_value, 1.0 / s.T, _SHAPE)
            self.assertGreater(np.abs(dens.rho_r[:, 0, 1]).max(), 1e-6)   # orbital coherences present
            tabs = flex_hf.build_flex_hf_tables({"PairLift": s.ham_info.param_ham["PairLift"]}, 2, _SHAPE)
            self.assertLess(np.abs(flex_hf.hf_map(dens.rho_r, tabs, _SHAPE, 2)).max(), 1e-14)

    def test_c_quadratic_coefficients_equal_the_hand_enumerated_sopt(self):
        for k in ("c20", "c11", "c02"):
            self.assertLess(_rel(self.on[k], self.oracle[k]), 1e-8, k)
        # the exchange skeleton is load-bearing (nonzero) and the Hubbard-only
        # coefficients carry none
        self.assertGreater(np.linalg.norm(self.oracle["x_v"].ravel()), 1e-2 * np.linalg.norm(self.oracle["d_v"].ravel()))
        self.assertLess(np.abs(self.oracle["x_u"]).max(), 1e-14)
        self.assertLess(np.abs(self.oracle["x_uv"] - self.oracle["x_v_only"]).max(), 1e-14)

    def test_c_general_path_second_order_recorded(self):
        """Recorded finding (not a Phase B defect): the general path's
        Takimoto-Hotta-Ueda subtraction -1/4 (S+C) chi (S+C) with the
        density-slot placement (S, C) = (0, 2 V(q)) of a spin-independent
        off-site V removes half of the direct V^2 skeleton and the whole
        U V cross term at second order, and carries no exchange skeleton."""
        self.assertLess(_rel(self.off["c02"], 0.5 * self.oracle["d_v"]), 1e-8)
        self.assertLess(np.abs(self.off["c11"]).max(), 1e-8 * np.abs(self.oracle["c11"]).max())
        self.assertLess(_rel(self.off["c20"], self.oracle["c20"]), 1e-8)

    def test_d_declared_zero_reproduces_the_general_path(self):
        g_on, g_off = self.g_on, self.g_off
        for (p, a), (q, b) in zip(zip(g_on.points, g_on.fluct), zip(g_off.points, g_off.fluct)):
            if p[1] == 0.0:            # V = 0 points: the two paths coincide exactly
                self.assertEqual(p, q)
                np.testing.assert_allclose(a, b, rtol=0, atol=1e-12 * np.abs(b).max())


if __name__ == "__main__":
    unittest.main()
