"""The bond-aware self-energy transport (spec 2026-09-06 section 3.4) and
gate G1: the B = 1 reduction to the general path and the direct-sum
oracle in (k, q, tau) space with explicit external-leg phases."""
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k

_IN2 = "tests/rpa/input_2orb"


def _flex():
    import hwave.solver.flex as flex_mod
    idict = {"path_to_input": _IN2, "Geometry": "geom.dat", "Transfer": "transfer.dat",
             "CoulombInter": "coulombinter.dat"}
    r = read_input_k.QLMSkInput({"path_to_input": _IN2, "interaction": idict})
    par = {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1],
           "Nmat": 8, "IterationMax": 1, "Mix": 1.0, "EPS": 1}
    info = {"mode": "FLEX", "param": par, "enable_spin_orbital": False, "calc_scheme": "general"}
    return flex_mod.FLEX(r.get_param("ham"), {}, info), r


def _green(s, beta):
    nmat, nvol, norb = s.nmat, s.lattice.nvol, s.norb
    return s._calc_dressed_green(beta, 0.1, np.zeros((1, nmat, nvol, norb, norb), complex))


class TestTransport(unittest.TestCase):

    def test_b1_reduces_to_calc_self_energy_general(self):
        from hwave.solver import bond_channels as bc
        from hwave.solver.flex_bond import BondBlockStore, calc_self_energy_bond
        s, _ = _flex()
        import os; os.makedirs("tests/flex/output", exist_ok=True)
        nmat, nvol, norb, nd = s.nmat, s.lattice.nvol, s.norb, s.norb ** 2
        beta = 0.5
        rng = np.random.default_rng(0)
        W = rng.normal(size=(nmat, nvol, nd, nd)) + 1j * rng.normal(size=(nmat, nvol, nd, nd))
        topo = bc.resolve_bond_topology({}, np.eye(3), norb, active_types=bc._LONGITUDINAL_ACTIVE_TYPES)
        view = bc.BondSetView(topo)
        s._calc_epsilon_k({})
        G = _green(s, beta)
        with BondBlockStore(nmat, nvol, nd, nd, ("W",)) as store:
            store.put_freq_batch("W", 0, nmat, W)
            sig = calc_self_energy_bond(store, G, beta, view, (4, 4, 1), norb, 1)
        ref = s._calc_self_energy_general(G, W, beta)
        np.testing.assert_array_equal(sig, ref)

    def test_hermitian_consistency_of_the_transport(self):
        """Sigma(k, iw)^dagger == Sigma(k, -iw) to round-off for a W with the
        bubble's conjugation symmetry W(q, i nu)^dagger = W(q, -i nu); the
        mixed-leg phase assignment breaks this at O(1e-3)."""
        from hwave.solver import bond_channels as bc
        from hwave.solver.flex_bond import BondBlockStore, calc_self_energy_bond
        s, _ = _flex()
        s._calc_epsilon_k({})
        nmat, nvol, norb, nd = s.nmat, s.lattice.nvol, s.norb, s.norb ** 2
        decl = {"CoulombInter": {((1, 0, 0), (0, 1)): 0.3, ((-1, 0, 0), (1, 0)): 0.3,
                                 ((0, 1, 0), (0, 0)): 0.2, ((0, -1, 0), (0, 0)): 0.2}}
        topo = bc.resolve_bond_topology(decl, np.eye(3), norb, active_types=bc._LONGITUDINAL_ACTIVE_TYPES)
        view = bc.BondSetView(topo)
        ND = view.n_channels * nd
        rng = np.random.default_rng(7)
        X = rng.normal(size=(nmat, nvol, ND, ND)) + 1j * rng.normal(size=(nmat, nvol, ND, ND))
        # impose W(q, l)^dagger = W(q, nmat - l) (l -> -l on the centred grid)
        W = np.empty_like(X)
        for l in range(nmat):
            lm = (nmat - l) % nmat
            W[l] = 0.5 * (X[l] + X[lm].conj().swapaxes(-1, -2))
        G = _green(s, 0.5)
        with BondBlockStore(nmat, nvol, ND, nd, ("W",)) as store:
            store.put_freq_batch("W", 0, nmat, W)
            sig = calc_self_energy_bond(store, G, 0.5, view, (4, 4, 1), norb, 1)[0]
        dev = np.max(np.abs(sig.conj().swapaxes(-1, -2) - sig[::-1])) / np.max(np.abs(sig))
        self.assertLess(dev, 1e-12)

    def test_g1_direct_sum_oracle(self):
        """4x4, norb = 2, nmat = 8, B = 5 (declared +-x, +-y), complex hopping,
        a random non-symmetric W: the (k, q, tau) oracle with the explicit
        bond form factor e^{+i(k-q).(R_alpha - R_beta)} on the internal leg
        and k - q reduced modulo the mesh equals the production transport."""
        from hwave.solver import bond_channels as bc
        from hwave.solver.flex_bond import BondBlockStore, calc_self_energy_bond
        from hwave.solver import matsubara as _ms
        s, _ = _flex()
        s._calc_epsilon_k({})
        nmat, nvol, norb, nd = s.nmat, s.lattice.nvol, s.norb, s.norb ** 2
        nx, ny, nz = 4, 4, 1
        beta = 0.5
        decl = {"CoulombInter": {((1, 0, 0), (0, 1)): 0.3, ((-1, 0, 0), (1, 0)): 0.3,
                                 ((0, 1, 0), (0, 0)): 0.2, ((0, -1, 0), (0, 0)): 0.2}}
        topo = bc.resolve_bond_topology(decl, np.eye(3), norb, active_types=bc._LONGITUDINAL_ACTIVE_TYPES)
        view = bc.BondSetView(topo)
        B = view.n_channels; ND = B * nd
        rng = np.random.default_rng(1)
        W = rng.normal(size=(nmat, nvol, ND, ND)) + 1j * rng.normal(size=(nmat, nvol, ND, ND))
        G = _green(s, beta)                                    # (1, nmat, nvol, norb, norb)
        with BondBlockStore(nmat, nvol, ND, nd, ("W",)) as store:
            store.put_freq_batch("W", 0, nmat, W)
            sig = calc_self_energy_bond(store, G, beta, view, (nx, ny, nz), norb, 1)
        # ---- oracle: tau-space product, explicit k, q sums and phases ----
        kx = 2 * np.pi * np.arange(nx) / nx; ky = 2 * np.pi * np.arange(ny) / ny
        kvec = np.array([(kx[i], ky[j], 0.0) for i in range(nx) for j in range(ny) for _ in range(nz)])
        G_t = _ms.fermion_to_tau(G[0].reshape(nmat, nvol * norb * norb), axis=0).reshape(nmat, nvol, norb, norb)
        W_t = _ms.boson_to_tau(W.reshape(nmat, nvol * ND * ND), axis=0).reshape(nmat, nvol, B, nd, B, nd)
        idx = np.arange(nvol).reshape(nx, ny, nz)
        sig_t = np.zeros((nmat, nvol, norb, norb), complex)
        R = np.array(view.delta_r)                             # (B, 3)
        for k in range(nvol):
            ik = np.unravel_index(k, (nx, ny, nz))
            for q in range(nvol):
                iq = np.unravel_index(q, (nx, ny, nz))
                kmq = idx[(ik[0] - iq[0]) % nx, (ik[1] - iq[1]) % ny, (ik[2] - iq[2]) % nz]
                kk = kvec[k]; kq = kvec[kmq]
                for a_ in range(B):
                    for b_ in range(B):
                        ph = np.exp(1j * kq @ (R[a_] - R[b_]))
                        Wblk = W_t[:, q, a_, :, b_, :].reshape(nmat, norb, norb, norb, norb)   # (t, c, a, d, b)
                        sig_t[:, k] += ph * np.einsum('tcadb,tcd->tab', Wblk, G_t[:, kmq])
        sig_t /= nvol
        ref = _ms.tau_to_fermion(sig_t.reshape(nmat, nvol * norb * norb), axis=0).reshape(1, nmat, nvol, norb, norb) / beta
        np.testing.assert_allclose(sig, ref, rtol=1e-10, atol=1e-12)


if __name__ == "__main__":
    unittest.main()


class TestFrequencyOracleSanity(unittest.TestCase):

    def test_large_nmat_frequency_oracle_sanity(self):
        """Pure-frequency oracle with the naive index difference n - l at
        nmat = 256 (wrap terms negligible; loose tolerance), B = 2, 2x2 mesh."""
        from hwave.solver import bond_channels as bc
        from hwave.solver.flex_bond import BondBlockStore, calc_self_energy_bond
        nmat, nx, ny, nz, norb = 256, 2, 2, 1, 1
        nvol, nd = nx * ny * nz, 1
        beta = 4.0
        decl = {"CoulombInter": {((1, 0, 0), (0, 0)): 0.3, ((-1, 0, 0), (0, 0)): 0.3}}
        topo = bc.resolve_bond_topology(decl, np.eye(3), norb, active_types=bc._LONGITUDINAL_ACTIVE_TYPES)
        view = bc.BondSetView(topo)
        B = view.n_channels; ND = B * nd
        wn = (2 * np.arange(nmat) + 1 - nmat) * np.pi / beta
        nu = (2 * np.arange(nmat) - nmat) * np.pi / beta
        ek = np.array([-0.5, 0.2, 0.2, 0.9])
        G = 1.0 / (1j * wn[:, None] - ek[None, :])                       # (nmat, nvol)
        G5 = G.reshape(1, nmat, nvol, 1, 1)
        rng = np.random.default_rng(3)
        amp = rng.normal(size=(nvol, ND, ND)) + 1j * rng.normal(size=(nvol, ND, ND))
        W = amp[None] / (1.0 + nu[:, None, None, None] ** 2)             # decays ~ 1/nu^2
        with BondBlockStore(nmat, nvol, ND, nd, ("W",)) as store:
            store.put_freq_batch("W", 0, nmat, W)
            sig = calc_self_energy_bond(store, G5, beta, view, (nx, ny, nz), norb, 1)[0, :, :, 0, 0]
        kx = 2 * np.pi * np.arange(nx) / nx; ky = 2 * np.pi * np.arange(ny) / ny
        kvec = np.array([(kx[i], ky[j], 0.0) for i in range(nx) for j in range(ny)])
        idx = np.arange(nvol).reshape(nx, ny)
        R = np.array(view.delta_r)
        ref = np.zeros((nmat, nvol), complex)
        for k in range(nvol):
            ik = np.unravel_index(k, (nx, ny))
            for q in range(nvol):
                iq = np.unravel_index(q, (nx, ny))
                kmq = idx[(ik[0] - iq[0]) % nx, (ik[1] - iq[1]) % ny]
                for a_ in range(B):
                    for b_ in range(B):
                        ph = np.exp(1j * kvec[kmq] @ (R[a_] - R[b_]))
                        for n in range(nmat):
                            for l in range(nmat):
                                m = n - l + nmat // 2          # w_n - nu_l -> index n - l + nmat/2
                                if 0 <= m < nmat:
                                    ref[n, k] += ph * W[l, q, a_, b_] * G[m, kmq]
        ref /= beta * nvol
        mid = slice(nmat // 2 - 8, nmat // 2 + 8)
        np.testing.assert_allclose(sig[mid], ref[mid], rtol=1e-3, atol=1e-6)
