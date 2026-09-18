"""Bond-resolved dynamic pairing kernel (spec 10.1)."""
import unittest
from unittest import mock

import numpy as np

from tests.eliashberg_bond_fixtures import physical_fixture


def _dev(fx):
    from hwave.solver.flex_bond import BondDeviceContext
    return BondDeviceContext(np, fx["S"], fx["C"])


def _axes(beta, wmax=12.0, eps=1.0e-8):
    """(axF, axB); skips cleanly when sparse-ir is not installed (the guard
    tests/test_bubble_kernel.py uses)."""
    try:
        import sparse_ir  # noqa: F401
    except ImportError:
        raise unittest.SkipTest("sparse-ir not installed")
    from hwave.solver.ir_axis import IRAxis
    return (IRAxis(beta=beta, wmax=wmax, eps=eps, statistics="F"),
            IRAxis(beta=beta, wmax=wmax, eps=eps, statistics="B"))


class TestAccumulator(unittest.TestCase):
    def test_pairing_coefficients_are_pinned(self):
        """The normative coefficients of spec 4.2, pinned as literals so that a
        mutated COEFF cannot pass by moving the reference with it."""
        from hwave.solver.eliashberg_bond import COEFF
        self.assertEqual(COEFF["singlet"], (1.5, -0.5))
        self.assertEqual(COEFF["triplet"], (-0.5, -0.5))
        self.assertEqual(set(COEFF), {"singlet", "triplet"})

    def test_accumulator_stages_agree(self):                      # 10.1.7
        from hwave.solver.eliashberg_bond import PairVertexAccumulator, COEFF
        fx = physical_fixture(norb=2, shape=(4, 2, 1), nmat=4)
        for eta in ("singlet", "triplet"):
            with _dev(fx) as dev:
                acc = PairVertexAccumulator(dev, pairing_types=(eta,), nb=3, nmat=fx["nmat"],
                                            nvol=fx["nvol"], nd=fx["nd"],
                                            spatial_shape=fx["spatial_shape"])
                acc.add_dressed(fx["store"], cond_tol=1e-6)
                v1 = acc.finish()[eta]
                acc2 = PairVertexAccumulator(dev, pairing_types=(eta,), nb=3, nmat=fx["nmat"],
                                             nvol=fx["nvol"], nd=fx["nd"],
                                             spatial_shape=fx["spatial_shape"])
                acc2.add_channel("spin", fx["store"], "chi_s_w")
                acc2.add_channel("charge", fx["store"], "chi_c_w")
                v2 = acc2.finish()[eta]
                g1 = v1.source.get_freq_batch(v1.slot, 0, fx["nmat"])
                g2 = v2.source.get_freq_batch(v2.slot, 0, fx["nmat"])
                np.testing.assert_array_equal(g1, g2)
                cs, cc = COEFF[eta]
                S, C = fx["S"], fx["C"]
                ref = cs * (S[None] @ fx["chi_s"] @ S[None]) + cc * (C[None] @ fx["chi_c"] @ C[None])
                np.testing.assert_allclose(g1, ref, rtol=0, atol=1e-13)
                # a missing stage is refused
                acc3 = PairVertexAccumulator(dev, pairing_types=(eta,), nb=3, nmat=fx["nmat"],
                                             nvol=fx["nvol"], nd=fx["nd"],
                                             spatial_shape=fx["spatial_shape"])
                acc3.add_channel("spin", fx["store"], "chi_s_w")
                with self.assertRaisesRegex(ValueError, "charge stage"):
                    acc3.finish()

    def test_uniform_refuses_two_types_and_non_finite(self):
        from hwave.solver.eliashberg_bond import PairVertexAccumulator, ArrayBlockSource
        fx = physical_fixture(norb=1, shape=(2, 2, 1), nmat=4)
        with _dev(fx) as dev:
            with self.assertRaisesRegex(ValueError, "exactly one pairing type"):
                PairVertexAccumulator(dev, pairing_types=("singlet", "triplet"), nb=2,
                                      nmat=4, nvol=4, nd=1, spatial_shape=(2, 2, 1))
            bad = fx["chi_s"].copy()
            bad[1, 0, 0, 0] = np.nan
            src = ArrayBlockSource({"chi_s_w": bad}, 1)
            acc = PairVertexAccumulator(dev, pairing_types=("singlet",), nb=2, nmat=4, nvol=4,
                                        nd=1, spatial_shape=(2, 2, 1))
            with self.assertRaisesRegex(ValueError, r"non-finite .* batch \[0, 2\)"):
                acc.add_channel("spin", src, "chi_s_w")

    def test_ir_finish_refuses_a_stage_callable_that_replays_nothing(self):
        """A callable that replays no stage leaves r = a = 0, which would read as
        a perfect fit AND silently disable the constant-vs-scale refusal."""
        from hwave.solver.eliashberg_bond import PairVertexAccumulator
        fx = physical_fixture(norb=1, shape=(4, 4, 1), nmat=8, beta=2.0)
        axF, axB = _axes(fx["beta"])

        def stages(acc):
            acc.add_channel("spin", fx["store"], "chi_s_w")
            acc.add_channel("charge", fx["store"], "chi_c_w")

        with _dev(fx) as dev:
            def build():
                a = PairVertexAccumulator(dev, pairing_types=("singlet",), nb=3, nmat=fx["nmat"],
                                          nvol=fx["nvol"], nd=fx["nd"],
                                          spatial_shape=fx["spatial_shape"], ir=(axF, axB))
                stages(a)
                return a
            with self.assertRaisesRegex(ValueError, "replayed no stage"):
                build().finish(ir_fit_tol=1.0e-4, stage_callable=lambda a: None)
            # one stage only is refused too
            with self.assertRaisesRegex(ValueError, "did not replay the charge stage"):
                build().finish(ir_fit_tol=1.0e-4,
                               stage_callable=lambda a: a.add_channel("spin", fx["store"],
                                                                      "chi_s_w"))
            # and the full replay is accepted (a loose tolerance: this test is
            # about the guard, not about the fit quality). The module logger is
            # qlms.eliashberg_bond; capturing it also keeps the run's output clean.
            with self.assertLogs("qlms.eliashberg_bond", level="INFO") as cm:
                out = build().finish(ir_fit_tol=1.0, stage_callable=stages)
            self.assertTrue(any("componentwise relative residual" in m for m in cm.output))
            self.assertTrue(np.all(np.isfinite(out["singlet"].fit_residual_rel)))
            self.assertGreater(float(out["singlet"].fit_residual_rel.max()), 0.0)


class TestAdmission(unittest.TestCase):
    def test_admission_table(self):                              # 10.1.11
        from hwave.solver.eliashberg_bond import estimate_pair_memory
        tab = estimate_pair_memory(nmat=8, ntau=0, nfreq=8, nvol=16, norb=2, B=3,
                                   num_eigenvalues=2, residency="auto", ir=False, L_B=0,
                                   n_channels=1, in_process=False)
        nd, ND, vec = 4, 12, 4 * 16 * 8 * 16
        blk = 8 * 16 * 16 * 16
        self.assertEqual(tab.rows["G2"], 16 * 16 * 8 * 16)
        self.assertEqual(tab.rows["gap_work"], (3 + 3) * vec)
        self.assertEqual(tab.rows["hoisted_blocks"], 9 * blk)
        self.assertEqual(tab.rows["stream_workspace"], 3 * blk)
        self.assertEqual(tab.rows["eigen_vectors"], (2 * 2 + 1 + 2) * vec)
        self.assertEqual(tab.rows["vertex_slot"], 8 * 16 * ND * ND * 16)   # post-processing uniform
        self.assertEqual(tab.rows["source_member"], 8 * 16 * ND * ND * 16)
        self.assertEqual(tab.rows["dressing_workspace"], 0)                # nb omitted
        h, d = tab.need("device")
        self.assertEqual(d, int(1.25 * (tab.rows["G2"] + tab.rows["gap_work"] + 9 * blk + 3 * blk
                                        + tab.rows["dressing_workspace"])))
        self.assertEqual(h, int(1.25 * (tab.rows["vertex_slot"] + tab.rows["source_member"]
                                        + tab.rows["eigen_vectors"])))
        big = 10 ** 12
        self.assertEqual(tab.choose("auto", host_cap=big, device_cap=big), "device")
        self.assertEqual(tab.choose("auto", host_cap=big, device_cap=tab.need("host")[1]), "host")
        self.assertEqual(tab.choose("auto", host_cap=tab.need("stream")[0],
                                    device_cap=tab.need("stream")[1]), "stream")
        with self.assertRaisesRegex(MemoryError, "stream"):
            tab.choose("auto", host_cap=1, device_cap=1)
        with self.assertRaisesRegex(MemoryError, "device"):
            tab.choose("device", host_cap=big, device_cap=1)
        tab_ip = estimate_pair_memory(nmat=8, ntau=0, nfreq=8, nvol=16, norb=2, B=3,
                                      num_eigenvalues=2, residency="auto", ir=False, L_B=0,
                                      n_channels=1, in_process=True)
        self.assertEqual(tab_ip.rows["vertex_slot"], 0)
        self.assertEqual(tab_ip.rows["source_member"], 0)
        # the dressing batch is a build-phase DEVICE peak: it enters every
        # residency's device need, the streaming one included
        tab_nb = estimate_pair_memory(nmat=8, ntau=0, nfreq=8, nvol=16, norb=2, B=3,
                                      num_eigenvalues=2, residency="auto", ir=False, L_B=0,
                                      n_channels=1, in_process=False, nb=2)
        self.assertEqual(tab_nb.rows["dressing_workspace"], 7 * 2 * 16 * ND * ND * 16)
        self.assertEqual(tab_nb.need("stream")[1],
                         int(1.25 * (tab_nb.rows["G2"] + tab_nb.rows["gap_work"]
                                     + tab_nb.rows["stream_workspace"]
                                     + tab_nb.rows["dressing_workspace"])))
        self.assertEqual(tab_nb.need("device")[1],
                         int(1.25 * (tab_nb.rows["G2"] + tab_nb.rows["gap_work"]
                                     + tab_nb.rows["stream_workspace"]
                                     + tab_nb.rows["dressing_workspace"]
                                     + tab_nb.rows["hoisted_blocks"])))
        self.assertEqual(tab_nb.need("host")[1], tab_nb.need("stream")[1])
        tab_ir = estimate_pair_memory(nmat=8, ntau=5, nfreq=6, nvol=16, norb=2, B=3,
                                      num_eigenvalues=2, residency="auto", ir=True, L_B=4,
                                      n_channels=2, in_process=True)
        self.assertEqual(tab_ir.rows["coefficients_build"], 2 * 2 * 5 * 16 * ND * ND * 16)
        self.assertEqual(tab_ir.rows["coefficients"], 2 * 5 * 16 * ND * ND * 16)
        self.assertEqual(tab_ir.rows["hoisted_blocks"], 9 * 5 * 16 * 16 * 16)


def _b1_topology(norb):
    """A ``B = 1`` (on-site only) topology.

    ``resolve_bond_topology`` always refuses a declaration set without an
    off-site shell, so the degenerate single-channel topology the on-site
    equivalence test (10.1.1) needs is built through the ``BondTopology``
    constructor, which re-validates every invariant (channel 0 is R = 0, the
    reversal involution, the exactly-zero channel-0 coefficients).
    """
    from hwave.solver import bond_channels
    return bond_channels.BondTopology(
        delta_r=[[0, 0, 0]], reverse=[0],
        coeffs={"CoulombInter": np.zeros((1, norb, norb), complex)})


def _g2(fx):
    from hwave.solver.eliashberg_dynamic import calc_g2_dynamic
    return calc_g2_dynamic(fx["green_sc"], fx["beta"])


def _kernel(fx, vertex, residency="host", V_inst=None):
    from hwave.solver.eliashberg_bond import BondPairKernel
    return BondPairKernel(vertex, _g2(fx), fx["view"], xp=np,
                          spatial_shape=fx["spatial_shape"], norb=fx["norb"],
                          beta=fx["beta"], nfreq=fx["nmat"], V_inst=V_inst,
                          residency=residency, host_cap=10 ** 12, device_cap=10 ** 12)


def _uniform_vertex(fx, eta):
    from hwave.solver.eliashberg_bond import PairVertexAccumulator
    with _dev(fx) as dev:
        acc = PairVertexAccumulator(dev, pairing_types=(eta,), nb=4, nmat=fx["nmat"],
                                    nvol=fx["nvol"], nd=fx["nd"],
                                    spatial_shape=fx["spatial_shape"])
        acc.add_channel("spin", fx["store"], "chi_s_w")
        acc.add_channel("charge", fx["store"], "chi_c_w")
        return acc.finish()[eta]


class TestKernelUniform(unittest.TestCase):
    def test_b1_matches_onsite_dynamic_kernel(self):              # 10.1.1
        from hwave.solver import eliashberg_dynamic as ed
        from hwave.solver.eliashberg_bond import (PairVertexUniform, ArrayBlockSource,
                                                  BondPairKernel, COEFF, instantaneous_vertex)
        from hwave.solver.bond_channels import BondSetView
        rng = np.random.default_rng(1)
        for norb in (1, 2):
            nd = norb * norb
            shape, nmat, beta = (2, 2, 1), 4, 1.5
            nvol = 4

            def rc(*s):
                return rng.standard_normal(s) + 1j * rng.standard_normal(s)

            chi_s, chi_c = rc(nmat, nvol, nd, nd), rc(nmat, nvol, nd, nd)
            S0, C0 = rc(nvol, nd, nd), rc(nvol, nd, nd)
            G2 = rc(norb, norb, norb, norb, 2, 2, 1, nmat)
            for eta in ("singlet", "triplet"):
                cs, cc = COEFF[eta]
                Gam = cs * (S0[None] @ chi_s @ S0[None]) + cc * (C0[None] @ chi_c @ C0[None])
                V_inst = instantaneous_vertex(S0, C0, nd, eta, shape)
                # on-site reference: V(q, l) in the sc.py layout, the
                # instantaneous part added at every l
                V_ref = np.empty((norb, norb, norb, norb, 2, 2, 1, nmat), complex)
                for l in range(nmat):
                    V_ref[..., l] = Gam[l].reshape(2, 2, 1, norb, norb, norb, norb).transpose(
                        3, 4, 5, 6, 0, 1, 2) + V_inst
                phi = rc(norb, norb, 2, 2, 1, nmat)
                ref = ed.eliashberg_kernel_dynamic(V_ref, G2, phi, norb, beta)
                view = BondSetView(_b1_topology(norb))
                src = ArrayBlockSource({"Gamma": Gam}, nd)
                vert = PairVertexUniform(source=src, slot="Gamma", B=1, nd=nd, nmat=nmat,
                                         nvol=nvol)
                for residency in ("host", "stream"):
                    K = BondPairKernel(vert, G2, view, xp=np, spatial_shape=shape, norb=norb,
                                       beta=beta, nfreq=nmat, V_inst=V_inst,
                                       residency=residency, host_cap=10 ** 12,
                                       device_cap=10 ** 12)
                    out = K.matvec(phi.ravel()).reshape(phi.shape)
                    np.testing.assert_allclose(out, ref, rtol=0, atol=1e-12,
                                               err_msg="norb {} {} {}".format(norb, eta,
                                                                              residency))

    def test_flat_frequency_matches_static_bond_kernel(self):     # 10.1.2
        from hwave.solver import bond_channels as bc
        from hwave.solver.eliashberg_bond import (PairVertexUniform, ArrayBlockSource,
                                                  BondPairKernel, COEFF)
        fx = physical_fixture(norb=1, shape=(4, 4, 1), nmat=8, beta=2.0,
                              delta_r=((0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)))
        nx, ny, nz = fx["spatial_shape"]
        nmat, nvol, nd, ND = fx["nmat"], fx["nvol"], fx["nd"], fx["ND"]
        # static side: the ED-adjudicated bond kernel on the SAME topology and vertices
        coul = {(tuple(R), (0, 0)): 0.5 for R in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0))}
        bond_set = bc.resolve_interactions(coul, np.eye(3), 1)
        S0q = fx["S"][:, :nd, :nd].reshape(nx, ny, nz, nd, nd)
        C0q = fx["C"][:, :nd, :nd].reshape(nx, ny, nz, nd, nd)
        S_b, C_b, Vpp_s, Vpp_t = bc.bare_bond_vertices(bond_set, S0q, C0q, 1)
        np.testing.assert_allclose(S_b.reshape(nvol, ND, ND), fx["S"], atol=1e-13)
        np.testing.assert_allclose(C_b.reshape(nvol, ND, ND), fx["C"], atol=1e-13)
        l0 = nmat // 2
        chi_s_st = fx["chi_s"][l0].reshape(nx, ny, nz, ND, ND)
        chi_c_st = fx["chi_c"][l0].reshape(nx, ny, nz, ND, ND)
        rng = np.random.default_rng(5)
        phi_k = rng.standard_normal((1, 1, nx, ny, nz)) + 1j * rng.standard_normal(
            (1, 1, nx, ny, nz))
        for eta in ("singlet", "triplet"):
            A_st, _ = bc.make_bond_kernel(chi_s_st, chi_c_st, S_b, C_b, Vpp_s, Vpp_t,
                                          fx["green_sc"], bond_set, eta, fx["beta"],
                                          part="fluctuation")
            ref = A_st.matvec(phi_k.ravel()).reshape(phi_k.shape)
            cs, cc = COEFF[eta]
            Gam_flat = np.broadcast_to(
                (cs * fx["S"] @ chi_s_st.reshape(nvol, ND, ND) @ fx["S"]
                 + cc * fx["C"] @ chi_c_st.reshape(nvol, ND, ND) @ fx["C"])[None],
                (nmat, nvol, ND, ND)).copy()
            vert = PairVertexUniform(source=ArrayBlockSource({"G": Gam_flat}, nd), slot="G",
                                     B=fx["B"], nd=nd, nmat=nmat, nvol=nvol)
            K = BondPairKernel(vert, _g2(fx), fx["view"], xp=np,
                               spatial_shape=fx["spatial_shape"], norb=1, beta=fx["beta"],
                               nfreq=nmat, residency="host", host_cap=10 ** 12,
                               device_cap=10 ** 12)
            phi_w = np.broadcast_to(phi_k[..., None], phi_k.shape + (nmat,)).copy()
            out = K.matvec(phi_w.ravel()).reshape(phi_w.shape)
            # a flat-frequency vertex acting on a flat gap gives a flat output
            np.testing.assert_allclose(out, np.broadcast_to(out[..., :1], out.shape),
                                       atol=1e-10)
            np.testing.assert_allclose(out[..., 0], ref, rtol=0, atol=1e-10, err_msg=eta)

    def test_bruteforce_oracle_norb2_bonds(self):                 # 10.1.3
        fx = physical_fixture(norb=2, shape=(4, 2, 1), nmat=4, beta=1.0,
                              coeffs={("CoulombInter", (1, 0, 0)):
                                      np.array([[0.3, 0.2], [0.2, 0.1]]),
                                      ("CoulombInter", (-1, 0, 0)):
                                      np.array([[0.3, 0.2], [0.2, 0.1]])})
        nx, ny, nz = fx["spatial_shape"]
        nmat, nvol, nd, B, norb = fx["nmat"], fx["nvol"], fx["nd"], fx["B"], 2
        G2 = _g2(fx)
        rng = np.random.default_rng(7)
        phi = (rng.standard_normal((norb, norb, nx, ny, nz, nmat))
               + 1j * rng.standard_normal((norb, norb, nx, ny, nz, nmat)))
        F = np.einsum("iljmxyzn,lmxyzn->ijxyzn", G2, phi)          # (l2, l3, k', w')
        kx = 2 * np.pi * np.arange(nx) / nx
        ky = 2 * np.pi * np.arange(ny) / ny
        dr = np.asarray(fx["view"].delta_r)
        for eta in ("singlet", "triplet"):
            vert = _uniform_vertex(fx, eta)
            Gam = vert.source.get_freq_batch(vert.slot, 0, nmat)   # (l, q, ND, ND)
            ref = np.zeros_like(phi)
            for ix in range(nx):
                for iy in range(ny):
                    for n in range(nmat):
                        for jx in range(nx):
                            for jy in range(ny):
                                qx, qy = (ix - jx) % nx, (iy - jy) % ny
                                q = (qx * ny + qy) * nz
                                for m in range(nmat):
                                    # The uniform-grid tau product is a CIRCULAR
                                    # convolution on the bosonic axis: the phase
                                    # bookkeeping of fermion_to_tau / boson_to_tau /
                                    # tau_to_fermion pins the bosonic index to
                                    # l = n - m + nmat//2 taken MODULO nmat.
                                    l = (n - m + nmat // 2) % nmat
                                    for al in range(B):
                                        pa = np.exp(1j * (kx[ix] * dr[al, 0]
                                                          + ky[iy] * dr[al, 1]))
                                        for be in range(B):
                                            pb = np.exp(1j * (kx[jx] * dr[be, 0]
                                                              + ky[jy] * dr[be, 1]))
                                            blk = Gam[l, q, al * nd:(al + 1) * nd,
                                                      be * nd:(be + 1) * nd]
                                            blk4 = blk.reshape(norb, norb, norb, norb)
                                            ref[:, :, ix, iy, 0, n] += -(1.0 / nvol) * pa * pb * (
                                                np.einsum("abcd,bc->ad", blk4,
                                                          F[:, :, jx, jy, 0, m]))
            K = _kernel(fx, vert)
            out = K.matvec(phi.ravel()).reshape(phi.shape)
            np.testing.assert_allclose(out, ref, rtol=0, atol=1e-12, err_msg=eta)

    def test_parity_commutation(self):                            # 10.1.4
        from hwave.solver.eliashberg_dynamic import _parity_leakage
        from scipy.sparse.linalg import LinearOperator
        fx = physical_fixture(norb=2, shape=(4, 2, 1), nmat=4, beta=1.0)
        for eta in ("singlet", "triplet"):
            K = _kernel(fx, _uniform_vertex(fx, eta))
            n = int(np.prod(K.gap_shape))
            A = LinearOperator((n, n), matvec=K.matvec, dtype=complex)
            self.assertLessEqual(_parity_leakage(A, K.gap_shape, eta), 1e-10, eta)

    def test_parity_leakage_refused(self):                        # 10.1.4b
        from hwave.solver.eliashberg_dynamic import run_leading_eigenproblem, build_seed
        from hwave.solver.eliashberg_bond import ArrayBlockSource, PairVertexUniform
        fx = physical_fixture(norb=1, shape=(4, 2, 1), nmat=4, beta=1.0)
        vert = _uniform_vertex(fx, "singlet")
        Gam = vert.source.get_freq_batch(vert.slot, 0, fx["nmat"]).copy()
        Gam[1, 2, 0, 1] += 0.7          # break the reversal symmetry of one block
        vert2 = PairVertexUniform(ArrayBlockSource({"G": Gam}, 1), "G", fx["B"], 1,
                                  fx["nmat"], fx["nvol"])
        K = _kernel(fx, vert2)
        nx, ny, nz = fx["spatial_shape"]
        kx, ky, kz = (2 * np.pi * np.arange(n) / n for n in (nx, ny, nz))
        phi0, seed = build_seed({}, "singlet", 1, kx, ky, kz, K.gap_shape, False, None,
                                fx["nmat"])
        for mode in ("iteration", "eigenvalue"):
            eli = {"solver_mode": mode, "eigenvalue_method": "arnoldi",
                   "num_eigenvalues": 2, "max_iter": 5}
            with self.assertRaisesRegex(ValueError, "does not commute with the combined parity"):
                run_leading_eigenproblem(K.matvec, K.gap_shape, eli, "singlet", phi0=phi0,
                                         seed_vec=seed, use_ir=False, axF=None,
                                         nmat=fx["nmat"], parity_leakage_policy="refuse")

    def test_dense_oracle_leading_eigenvalue(self):               # 10.1.5
        import scipy.linalg
        import hwave.sc as sc
        from scipy.sparse.linalg import LinearOperator
        fx = physical_fixture(norb=1, shape=(2, 2, 1), nmat=4, beta=1.0)
        K = _kernel(fx, _uniform_vertex(fx, "singlet"))
        n = int(np.prod(K.gap_shape))
        M = np.empty((n, n), complex)
        for j in range(n):
            e = np.zeros(n, complex)
            e[j] = 1.0
            M[:, j] = K.matvec(e)
        ev = scipy.linalg.eigvals(M)
        lead = ev[np.argmax(ev.real)]

        def make():
            return LinearOperator((n, n), matvec=K.matvec, dtype=complex), n

        lam, _, _ = sc._solve_leading(make, n, "arnoldi", num_eigenvalues=3)
        self.assertAlmostEqual(float(np.real(lam)), float(lead.real),
                               delta=1e-8 * max(1.0, abs(lead)))

    def test_residency_modes_agree(self):                         # 10.1.6
        from hwave.solver import backend
        from hwave.solver.eliashberg_bond import BondPairKernel
        fx = physical_fixture(norb=2, shape=(4, 2, 1), nmat=4, beta=1.0)
        vert = _uniform_vertex(fx, "singlet")
        rng = np.random.default_rng(3)
        Kh = _kernel(fx, vert, residency="host")
        phi = rng.standard_normal(int(np.prod(Kh.gap_shape))) + 0j
        ref = Kh.matvec(phi)
        Ks = _kernel(fx, vert, residency="stream")
        np.testing.assert_allclose(Ks.matvec(phi), ref, atol=1e-13)
        Ka = _kernel(fx, vert, residency="auto")
        self.assertEqual(Ka.residency, "host")
        np.testing.assert_allclose(Ka.matvec(phi), ref, atol=1e-13)
        if not backend.gpu_available():
            raise unittest.SkipTest("CUDA device required for the device residency")
        import cupy
        Kd = BondPairKernel(vert, cupy.asarray(_g2(fx)), fx["view"], xp=cupy,
                            spatial_shape=fx["spatial_shape"], norb=2, beta=fx["beta"],
                            nfreq=4, residency="device", host_cap=10 ** 12,
                            device_cap=10 ** 12)
        np.testing.assert_allclose(Kd.matvec(phi), ref, atol=1e-13)

    def test_gap_bond_projection(self):                           # 10.1.10
        from hwave.solver.eliashberg_bond import gap_bond_projection
        fx = physical_fixture(norb=1, shape=(4, 4, 1), nmat=2)
        nx, ny, nz = fx["spatial_shape"]
        kx = 2 * np.pi * np.arange(nx) / nx
        dr = np.asarray(fx["view"].delta_r)
        m = 1
        gap = (np.exp(1j * kx[:, None, None] * dr[m, 0])[None, None, :, :, :, None]
               * np.ones((1, 1, nx, ny, nz, 2)))
        psi = gap_bond_projection(gap, fx["view"], fx["spatial_shape"])
        self.assertEqual(psi.shape, (fx["B"], 1, 1, 2))
        for mm in range(fx["B"]):
            np.testing.assert_allclose(psi[mm], 1.0 if mm == m else 0.0, atol=1e-12)

    def test_instantaneous_vertex_single_band_equals_onsite(self):   # 10.1.12
        import hwave.sc as sc
        from hwave.solver import eliashberg_dynamic as ed
        from hwave.solver.eliashberg_bond import instantaneous_vertex
        fx = physical_fixture(norb=1, shape=(4, 4, 1), nmat=2, U=2.0,
                              delta_r=((0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)))
        nx, ny, nz = fx["spatial_shape"]
        kx, ky, kz = (2 * np.pi * np.arange(n) / n for n in (nx, ny, nz))
        interactions = {"CoulombIntra": {((0, 0, 0), (0, 0)): 2.0},
                        "CoulombInter": {(tuple(R), (0, 0)): 0.5
                                         for R in ((1, 0, 0), (-1, 0, 0), (0, 1, 0),
                                                   (0, -1, 0))}}
        inter_k = sc._build_interaction_k(kx, ky, kz, interactions, 1)
        for eta in ("singlet", "triplet"):
            ref = ed._instantaneous_vertex(inter_k, 1, nx, ny, nz, eta, "myo")
            got = instantaneous_vertex(fx["S"], fx["C"], fx["nd"], eta, fx["spatial_shape"])
            np.testing.assert_allclose(got, ref, atol=1e-14, err_msg=eta)
        fx2 = physical_fixture(norb=2, shape=(4, 2, 1), nmat=2, U=2.0,
                               coeffs={("CoulombInter", (1, 0, 0)):
                                       np.array([[0.0, 0.4], [0.4, 0.0]]),
                                       ("CoulombInter", (-1, 0, 0)):
                                       np.array([[0.0, 0.4], [0.4, 0.0]])})
        interactions2 = {"CoulombIntra": {((0, 0, 0), (a, a)): 2.0 for a in range(2)},
                         "CoulombInter": {((1, 0, 0), (0, 1)): 0.4, ((1, 0, 0), (1, 0)): 0.4,
                                          ((-1, 0, 0), (0, 1)): 0.4, ((-1, 0, 0), (1, 0)): 0.4}}
        inter_k2 = sc._build_interaction_k(kx, 2 * np.pi * np.arange(2) / 2, kz,
                                           interactions2, 2)
        ref2 = ed._instantaneous_vertex(inter_k2, 2, 4, 2, 1, "singlet", "myo")
        got2 = instantaneous_vertex(fx2["S"], fx2["C"], 4, "singlet", (4, 2, 1))
        diff = np.abs(got2 - ref2) > 1e-12
        # only the cross slots (ab,ab)/(ba,ba) may differ: V[a,b,c,d] with (a,b)=(c,d), a != b
        allowed = np.zeros(diff.shape, bool)
        allowed[0, 1, 0, 1] = allowed[1, 0, 1, 0] = True
        self.assertFalse(np.any(diff & ~allowed))

    def test_instantaneous_operator_matches_static_bond_kernel(self):   # 10.1.13 (ARBITER)
        from hwave.solver import bond_channels as bc
        from hwave.solver.eliashberg_dynamic import _project_parity_dynamic
        from hwave.solver.eliashberg_bond import (PairVertexUniform, ArrayBlockSource,
                                                  BondPairKernel, instantaneous_vertex)
        fx = physical_fixture(norb=1, shape=(4, 4, 1), nmat=8, beta=2.0, U=2.0,
                              delta_r=((0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)))
        nx, ny, nz = fx["spatial_shape"]
        nmat, nvol, nd, ND = fx["nmat"], fx["nvol"], fx["nd"], fx["ND"]
        coul = {(tuple(R), (0, 0)): 0.5 for R in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0))}
        bond_set = bc.resolve_interactions(coul, np.eye(3), 1)
        S0q = fx["S"][:, :nd, :nd].reshape(nx, ny, nz, nd, nd)
        C0q = fx["C"][:, :nd, :nd].reshape(nx, ny, nz, nd, nd)
        S_b, C_b, Vpp_s, Vpp_t = bc.bare_bond_vertices(bond_set, S0q, C0q, 1)
        zero = np.zeros((nx, ny, nz, ND, ND), complex)
        G2 = _g2(fx)
        rng = np.random.default_rng(11)
        for eta in ("singlet", "triplet"):
            A_st, _ = bc.make_bond_kernel(zero, zero, S_b, C_b, Vpp_s, Vpp_t, fx["green_sc"],
                                          bond_set, eta, fx["beta"], part="instantaneous")
            V_inst = instantaneous_vertex(fx["S"], fx["C"], nd, eta, fx["spatial_shape"])
            Gam0 = np.zeros((nmat, nvol, ND, ND), complex)
            vert = PairVertexUniform(ArrayBlockSource({"G": Gam0}, nd), "G", fx["B"], nd,
                                     nmat, nvol)
            K = BondPairKernel(vert, G2, fx["view"], xp=np, spatial_shape=fx["spatial_shape"],
                               norb=1, beta=fx["beta"], nfreq=nmat, V_inst=V_inst,
                               residency="host", host_cap=10 ** 12, device_cap=10 ** 12)
            flat = np.ones((1, 1, nx, ny, nz, nmat), complex)
            rnd = (rng.standard_normal((1, 1, nx, ny, nz, nmat))
                   + 1j * rng.standard_normal((1, 1, nx, ny, nz, nmat)))
            odd = rnd - rnd[..., ::-1]
            for phi in (flat, rnd, odd):
                phi = _project_parity_dynamic(phi, eta)
                if np.linalg.norm(phi) == 0:
                    continue
                out = K.matvec(phi.ravel()).reshape(phi.shape)
                # frequency-flat output
                np.testing.assert_allclose(out, np.broadcast_to(out[..., :1], out.shape),
                                           atol=1e-10)
                # static side on the gap whose X equals the dynamic frequency-summed F:
                # X_static(phi_s) = G2_static phi_s with G2_static = T sum_n G G
                # == sum_n G2(n) phi(n)
                F_dyn = np.einsum("iljmxyzn,lmxyzn->ijxyzn", G2, phi).sum(-1)
                G2_st = bc._g2_from_green(fx["green_sc"], fx["beta"])
                phi_s = F_dyn / G2_st[0, 0, 0, 0]          # norb = 1: scalar per k
                ref = A_st.matvec(phi_s.ravel()).reshape(phi.shape[:-1])
                np.testing.assert_allclose(out[..., 0], ref, rtol=0, atol=1e-10, err_msg=eta)



#: A residual band wide enough that the deliberately non-smooth physical
#: fixture is accepted without even the warning decade of spec 4.5.
_WIDE_TOL = 100.0


def _ir_fixture():
    """The physical fixture the IR accumulator tests run on. Its uniform-FFT
    susceptibilities are deliberately NOT IR-representable (they carry the
    finite-Nmat delta(tau) constant and the aliasing images of the discrete
    tau transform), so these tests pin the ACCUMULATION and the residual
    bookkeeping, never the fit quality; the operator equivalence is gated on
    an IR-representable vertex in 10.1.14."""
    return physical_fixture(norb=1, shape=(4, 4, 1), nmat=64, beta=4.0, U=1.0)


def _stage_pair(source):
    """The two post-processing stages of the accumulator over one source."""
    def stages(acc):
        acc.add_channel("spin", source, "chi_s_w")
        acc.add_channel("charge", source, "chi_c_w")
    return stages


def _direct_residual(fx, axB, eta, keep):
    """The spec-4.5 componentwise relative residual computed directly from the
    resident contributions: per stage, ``r = max |eval(fit(G)) - G|`` and
    ``a = max |G|`` over frequency, momentum and the elements of each
    ``(alpha, beta)`` block, reported as ``max_stage r / a``."""
    from hwave.solver.eliashberg_bond import COEFF
    nd, B, nmat = fx["nd"], fx["B"], fx["nmat"]
    fit, ev = axB.uniform_matrices(nmat, with_constant=True)
    cs, cc = COEFF[eta]
    stages = (cs * (fx["S"][None] @ fx["chi_s"] @ fx["S"][None]),
              cc * (fx["C"][None] @ fx["chi_c"] @ fx["C"][None]))
    rel = np.zeros((B, B))
    for G in stages:
        sol = np.einsum("lvij,lc->cvij", G, fit)
        rec = np.einsum("cvij,cl->lvij", sol[:axB.L], ev)
        if keep:
            rec = rec + sol[axB.L][None]
        for al in range(B):
            for be in range(B):
                sl = (slice(None), slice(None), slice(al * nd, (al + 1) * nd),
                      slice(be * nd, (be + 1) * nd))
                a = float(np.abs(G[sl]).max())
                r = float(np.abs(rec[sl] - G[sl]).max())
                rel[al, be] = max(rel[al, be], r / a if a > 0 else 0.0)
    return rel


def _lorentzian_source(fx, g=1.3, kappa=0.0):
    """An ArrayBlockSource whose susceptibilities are EXACTLY representable in
    the bosonic IR basis: the fixture's own static (reversal-closed) bond
    blocks times a Lorentzian bosonic profile ``g^2 / (nu^2 + g^2)``, plus an
    optional frequency-FLAT offset ``kappa`` -- the delta(tau) component the
    smooth basis cannot represent and ``ir_keep_static`` retains."""
    from hwave.solver.eliashberg_bond import ArrayBlockSource
    nmat, beta = fx["nmat"], fx["beta"]
    nu = (2 * np.arange(nmat) - nmat) * np.pi / beta
    prof = (g * g / (nu ** 2 + g * g) + kappa)[:, None, None, None]
    return ArrayBlockSource({"chi_s_w": fx["chi_s"][nmat // 2][None] * prof,
                             "chi_c_w": fx["chi_c"][nmat // 2][None] * prof}, fx["nd"])


def _ir_vs_uniform_at_nodes(fx, axF, axB, eta, *, source, keep, with_vinst):
    """Build the uniform and the IR bond kernel from the SAME source and return
    the max-norm relative difference of their matvecs at the Matsubara nodes
    the two grids share (the comparison :func:`_smooth_vertex_gate` of
    ``tests/test_eliashberg_ir.py`` makes on the on-site path). Densifying the
    IR output instead would be meaningless here: a frequency-FLAT output term
    (the bare vertex, the retained constant) is by construction outside the
    fermionic IR basis."""
    from hwave.solver.eliashberg_bond import (PairVertexAccumulator, BondPairKernel,
                                              instantaneous_vertex)
    from hwave.solver.eliashberg_dynamic import calc_g2_dynamic, _ir_compress
    nmat, beta = fx["nmat"], fx["beta"]
    nx, ny, nz = fx["spatial_shape"]
    stages = _stage_pair(source)
    # both kernels see the same Green function: the IR nodes and the SAME
    # interpolant densified back onto the uniform grid
    green_ir = _ir_compress(fx["green_sc"], axF, nmat, "green")
    green_u = axF.eval_to_uniform(axF.fit_from_freq(green_ir), nmat)
    V_inst = (instantaneous_vertex(fx["S"], fx["C"], fx["nd"], eta, fx["spatial_shape"])
              if with_vinst else None)
    kx = 2 * np.pi * np.arange(nx) / nx
    ky = 2 * np.pi * np.arange(ny) / ny
    form = (np.cos(kx)[:, None] - np.cos(ky)[None, :])[None, None, :, :, None, None]
    with _dev(fx) as dev:
        au = PairVertexAccumulator(dev, pairing_types=(eta,), nb=16, nmat=nmat,
                                   nvol=fx["nvol"], nd=fx["nd"],
                                   spatial_shape=fx["spatial_shape"])
        stages(au)
        ai = PairVertexAccumulator(dev, pairing_types=(eta,), nb=16, nmat=nmat,
                                   nvol=fx["nvol"], nd=fx["nd"],
                                   spatial_shape=fx["spatial_shape"], ir=(axF, axB),
                                   ir_keep_static=keep)
        stages(ai)
        vu = au.finish()[eta]
        vi = ai.finish(ir_fit_tol=1e-3, stage_callable=stages)[eta]
    assert (vi.const is not None) == keep, "ir_keep_static did not reach the vertex"
    Ku = BondPairKernel(vu, calc_g2_dynamic(green_u, beta), fx["view"], xp=np,
                        spatial_shape=fx["spatial_shape"], norb=1, beta=beta, nfreq=nmat,
                        V_inst=V_inst, residency="host", host_cap=10 ** 12,
                        device_cap=10 ** 12)
    Ki = BondPairKernel(vi, calc_g2_dynamic(green_ir, beta), fx["view"], xp=np,
                        spatial_shape=fx["spatial_shape"], norb=1, beta=beta,
                        nfreq=axF.n_freq, V_inst=V_inst, axF=axF, residency="host",
                        host_cap=10 ** 12, device_cap=10 ** 12)
    phi_u = np.ascontiguousarray(np.broadcast_to(
        form, (1, 1, nx, ny, nz, nmat))).astype(complex)
    phi_i = np.ascontiguousarray(np.broadcast_to(
        form, (1, 1, nx, ny, nz, axF.n_freq))).astype(complex)
    out_u = Ku.matvec(phi_u.ravel()).reshape(phi_u.shape)
    out_i = Ki.matvec(phi_i.ravel()).reshape((1, 1, nx, ny, nz, axF.n_freq))
    idx = (axF.freq_n - 1 + nmat) // 2
    inside = (idx >= 0) & (idx < nmat)
    return float(np.abs(out_i[..., inside] - out_u[..., idx[inside]]).max()
                 / np.abs(out_u).max())


class TestKernelIR(unittest.TestCase):
    def test_ir_coefficient_accumulation_equals_direct_fit(self):     # 10.1.8
        from hwave.solver.eliashberg_bond import PairVertexAccumulator, COEFF
        fx = _ir_fixture()
        axF, axB = _axes(fx["beta"])
        nmat = fx["nmat"]
        calls = []

        def stages(acc):
            calls.append(1)
            _stage_pair(fx["store"])(acc)

        def direct(eta):
            cs, cc = COEFF[eta]
            Gam = (cs * (fx["S"][None] @ fx["chi_s"] @ fx["S"][None])
                   + cc * (fx["C"][None] @ fx["chi_c"] @ fx["C"][None]))
            return np.einsum("lvij,lc->cvij", Gam,
                             axB.uniform_matrices(nmat, with_constant=True)[0])

        with _dev(fx) as dev:
            acc = PairVertexAccumulator(dev, pairing_types=("singlet", "triplet"), nb=3,
                                        nmat=nmat, nvol=fx["nvol"], nd=fx["nd"],
                                        spatial_shape=fx["spatial_shape"], ir=(axF, axB))
            stages(acc)
            calls.clear()
            # the uniform-FFT vertex of the physical fixture is NOT smooth (see
            # _ir_fixture); this test is about the bookkeeping, so the residual
            # band is wide open (and well clear of the warning decade) and the
            # residual itself is pinned against a direct computation below
            both = acc.finish(ir_fit_tol=_WIDE_TOL, stage_callable=stages)
            self.assertEqual(len(calls), 1)     # the residual pass replays once
            for eta in ("singlet", "triplet"):
                # batched accumulation over nb = 3 == the one-shot fit
                np.testing.assert_allclose(both[eta].coeffs, direct(eta)[:axB.L], atol=1e-12)
                self.assertIsNone(both[eta].const)
                self.assertEqual(both[eta].fit_residual_rel.shape, (fx["B"], fx["B"]))
                np.testing.assert_allclose(both[eta].fit_residual_rel,
                                           _direct_residual(fx, axB, eta, keep=False),
                                           rtol=1e-12, err_msg=eta)
            # (the two types share one residual array by construction -- the
            # residual pass runs on the FIRST type only, the types differing by
            # scalar coefficients that cancel in r / a; the per-type equality
            # against _direct_residual above is what pins that)
            # "both" in one pass == two single passes
            for eta in ("singlet", "triplet"):
                a1 = PairVertexAccumulator(dev, pairing_types=(eta,), nb=3, nmat=nmat,
                                           nvol=fx["nvol"], nd=fx["nd"],
                                           spatial_shape=fx["spatial_shape"], ir=(axF, axB))
                stages(a1)
                one = a1.finish(ir_fit_tol=_WIDE_TOL, stage_callable=stages)[eta]
                np.testing.assert_array_equal(one.coeffs, both[eta].coeffs)
                # r and a carry the pairing coefficients, which cancel in the
                # ratio up to the floating-point order of the two reductions
                np.testing.assert_allclose(one.fit_residual_rel,
                                           both[eta].fit_residual_rel, rtol=1e-12)
            # ir_fit_tol = 0 skips the residual pass: the callable is NOT invoked
            a0 = PairVertexAccumulator(dev, pairing_types=("singlet",), nb=3, nmat=nmat,
                                       nvol=fx["nvol"], nd=fx["nd"],
                                       spatial_shape=fx["spatial_shape"], ir=(axF, axB))
            stages(a0)
            calls.clear()
            with self.assertLogs("qlms.eliashberg_bond", level="WARNING"):
                v0 = a0.finish(ir_fit_tol=0, stage_callable=stages)["singlet"]
            self.assertEqual(calls, [])
            self.assertTrue(np.all(np.isnan(v0.fit_residual_rel)))
            np.testing.assert_array_equal(v0.coeffs, both["singlet"].coeffs)
            # ir_keep_static retains the constant column of the augmented fit,
            # and the residual then credits it instead of counting it as misfit
            ak = PairVertexAccumulator(dev, pairing_types=("singlet",), nb=3, nmat=nmat,
                                       nvol=fx["nvol"], nd=fx["nd"],
                                       spatial_shape=fx["spatial_shape"], ir=(axF, axB),
                                       ir_keep_static=True)
            stages(ak)
            vk = ak.finish(ir_fit_tol=_WIDE_TOL, stage_callable=stages)["singlet"]
            np.testing.assert_allclose(vk.const, direct("singlet")[axB.L], atol=1e-12)
            np.testing.assert_array_equal(vk.coeffs, both["singlet"].coeffs)
            np.testing.assert_allclose(vk.fit_residual_rel,
                                       _direct_residual(fx, axB, "singlet", keep=True),
                                       rtol=1e-12)
            self.assertGreater(float(np.abs(vk.const).max()), 0.0)
            self.assertLess(float(vk.fit_residual_rel.max()),
                            float(both["singlet"].fit_residual_rel.max()))

    def test_ir_fit_tol_refusal(self):                                # 10.1.9
        from hwave.solver.eliashberg_bond import PairVertexAccumulator, ArrayBlockSource
        fx = _ir_fixture()
        axF, axB = _axes(fx["beta"])
        rng = np.random.default_rng(2)
        noise = (rng.standard_normal(fx["chi_s"].shape)
                 + 1j * rng.standard_normal(fx["chi_s"].shape))
        src = ArrayBlockSource({"chi_s_w": noise, "chi_c_w": noise}, fx["nd"])
        stages = _stage_pair(src)

        def build(dev):
            a = PairVertexAccumulator(dev, pairing_types=("singlet",), nb=8, nmat=fx["nmat"],
                                      nvol=fx["nvol"], nd=fx["nd"],
                                      spatial_shape=fx["spatial_shape"], ir=(axF, axB))
            stages(a)
            return a

        with _dev(fx) as dev:
            # white noise is not representable in a smooth basis: refused
            with self.assertRaisesRegex(ValueError, "raise ir_wmax, increase Nmat"):
                build(dev).finish(ir_fit_tol=1e-4, stage_callable=stages)
            # ir_fit_tol = 0 only warns that the check was skipped
            with self.assertLogs("qlms.eliashberg_bond", level="WARNING") as cm:
                build(dev).finish(ir_fit_tol=0, stage_callable=stages)
            self.assertTrue(any("SKIPPED" in m for m in cm.output))
            # a residual inside [0.1 * tol, tol) is accepted with a warning
            worst = float(build(dev).finish(
                ir_fit_tol=1e9, stage_callable=stages)["singlet"].fit_residual_rel.max())
            self.assertGreater(worst, 0.0)
            with self.assertLogs("qlms.eliashberg_bond", level="WARNING") as cm2:
                build(dev).finish(ir_fit_tol=1.5 * worst, stage_callable=stages)
            self.assertTrue(any("within a decade" in m for m in cm2.output))

    def test_ir_matvec_matches_uniform_matvec(self):                  # 10.1.14
        """The IR and the uniform bond matvec are the same operator.

        Part A pins the DYNAMIC part exactly: on a fully IR-representable
        (Lorentzian) bond vertex carried by the fixture's own vertices,
        topology and Green function, the two agree at the Matsubara nodes the
        grids share to well under the on-site path's ``1e-6``.

        Part B measures the FREQUENCY-FLAT terms -- the bare vertex ``V_inst``
        and the retained ``ir_keep_static`` constant. A frequency-flat vertex
        multiplies the UNREGULARIZED Matsubara sum ``(1/beta) sum_n F(i w_n)``,
        which on IR is the midpoint of the ``tau = 0`` jump,
        ``sum_l F_l * 0.5 * (u_l(0+) - u_l(beta-))`` -- the
        ``IRAxis.u_matsubara_sum`` weights, NOT the one-sided equal-time value
        ``F(0+)``. The uniform grid evaluates the same sum TRUNCATED at Nmat,
        so the two differ by O(beta / Nmat) by construction and the assertion
        is that difference's convergence: doubling Nmat halves it. Contracting
        with ``u_zero_plus`` instead (the half-jump anti-commutes with
        ``i w -> -i w``) gives an O(1) difference that does not converge, and
        breaks the parity commutation the kernel is projected on -- see
        ``TestKernelIR.test_ir_parity_commutation``.
        """
        beta = 2.0

        def fixture(nmat):
            # the axes first, so a machine without sparse-ir skips before
            # paying for the fixture
            axes = _axes(beta, wmax=20.0)
            return physical_fixture(norb=1, shape=(4, 4, 1), nmat=nmat, beta=beta,
                                    U=1.0), axes

        # -- Part A: the dynamic part, no frequency-flat term ---------------
        fx, (axF, axB) = fixture(256)
        for eta in ("singlet", "triplet"):
            d = _ir_vs_uniform_at_nodes(fx, axF, axB, eta,
                                        source=_lorentzian_source(fx), keep=False,
                                        with_vinst=False)
            self.assertLess(d, 1e-6, "{}: IR vs uniform {:.3e}".format(eta, d))

        # -- Part B: the frequency-flat terms converge as O(beta / Nmat) ----
        for eta in ("singlet", "triplet"):
            for keep, with_vinst in ((False, True), (True, False), (True, True)):
                diffs = []
                for nmat in (128, 256):
                    fx, (axF, axB) = fixture(nmat)
                    src = _lorentzian_source(fx, kappa=0.3 if keep else 0.0)
                    diffs.append(_ir_vs_uniform_at_nodes(fx, axF, axB, eta, source=src,
                                                         keep=keep, with_vinst=with_vinst))
                tag = "{} keep_static={} V_inst={}".format(eta, keep, with_vinst)
                # measured: 3.7e-3 -> 1.8e-3 (V_inst) and 1.5e-3 -> 7.5e-4
                # (retained constant); the equal-time truncation error itself
                # is 3.2e-3 -> 1.6e-3 on this fixture
                self.assertLess(diffs[0], 1e-2, tag)
                self.assertLess(diffs[1], 0.6 * diffs[0],
                                "{}: no Nmat convergence {:.3e} -> {:.3e}"
                                .format(tag, diffs[0], diffs[1]))

    def test_ir_parity_commutation(self):                             # 10.1.4 (IR)
        """The IR bond kernel must commute with the combined parity operator
        with BOTH frequency-flat terms live -- the bare vertex ``V_inst`` and
        the ``ir_keep_static`` constant.

        The uniform half (10.1.4) cannot see this: its dense tau grid
        represents the flat term as a single bin, which IS the (truncated)
        Matsubara sum. On IR the sum is evaluated analytically, and it must be
        the UNREGULARIZED sum, i.e. the midpoint of the tau jump. The one-sided
        ``F(0^+)`` carries half the jump on top, and that half ANTI-commutes
        with the frequency reversal: measured leakage 4.30e-01 before the fix
        on this fixture (an IR fit residual of 2.9e-10, so not a data-quality
        effect), 1.3e-12 after."""
        from hwave.solver.eliashberg_bond import (PairVertexAccumulator, BondPairKernel,
                                                  instantaneous_vertex)
        from hwave.solver.eliashberg_dynamic import (calc_g2_dynamic, _ir_compress,
                                                     _parity_leakage)
        from scipy.sparse.linalg import LinearOperator
        beta = 2.0
        axF, axB = _axes(beta, wmax=20.0)
        fx = physical_fixture(norb=1, shape=(4, 4, 1), nmat=128, beta=beta, U=1.0)
        stages = _stage_pair(_lorentzian_source(fx, kappa=0.3))
        with _dev(fx) as dev:
            acc = PairVertexAccumulator(dev, pairing_types=("singlet",), nb=16,
                                        nmat=fx["nmat"], nvol=fx["nvol"], nd=fx["nd"],
                                        spatial_shape=fx["spatial_shape"], ir=(axF, axB),
                                        ir_keep_static=True)
            stages(acc)
            vert = acc.finish(ir_fit_tol=1e-3, stage_callable=stages)["singlet"]
        self.assertIsNotNone(vert.const)            # the retained flat term is live
        G2 = calc_g2_dynamic(_ir_compress(fx["green_sc"], axF, fx["nmat"], "green"), beta)
        for eta in ("singlet", "triplet"):
            V_inst = instantaneous_vertex(fx["S"], fx["C"], fx["nd"], eta,
                                          fx["spatial_shape"])
            self.assertGreater(float(np.abs(V_inst).max()), 0.0)
            K = BondPairKernel(vert, G2, fx["view"], xp=np,
                               spatial_shape=fx["spatial_shape"], norb=fx["norb"],
                               beta=beta, nfreq=axF.n_freq, V_inst=V_inst, axF=axF,
                               residency="host", host_cap=10 ** 12, device_cap=10 ** 12)
            n = int(np.prod(K.gap_shape))
            A = LinearOperator((n, n), matvec=K.matvec, dtype=complex)
            self.assertLessEqual(_parity_leakage(A, K.gap_shape, eta), 1e-10, eta)

    def test_ir_residency_modes_agree(self):                          # 10.1.6 (IR)
        """The IR half of the residency contract: the hoisted blocks and the
        per-matvec rebuild of ``_block_rtau`` must give the same operator.
        Run with the bare vertex AND a retained constant, so the two
        frequency-flat terms (which live outside the hoisted blocks) are part
        of the comparison."""
        from hwave.solver import backend
        from hwave.solver.eliashberg_bond import (PairVertexAccumulator, BondPairKernel,
                                                  instantaneous_vertex)
        from hwave.solver.eliashberg_dynamic import calc_g2_dynamic, _ir_compress
        axF, axB = _axes(2.0, wmax=20.0)
        fx = physical_fixture(norb=2, shape=(4, 2, 1), nmat=64, beta=2.0, U=1.0)
        stages = _stage_pair(_lorentzian_source(fx, kappa=0.3))
        with _dev(fx) as dev:
            acc = PairVertexAccumulator(dev, pairing_types=("singlet",), nb=16,
                                        nmat=fx["nmat"], nvol=fx["nvol"], nd=fx["nd"],
                                        spatial_shape=fx["spatial_shape"], ir=(axF, axB),
                                        ir_keep_static=True)
            stages(acc)
            vert = acc.finish(ir_fit_tol=1e-3, stage_callable=stages)["singlet"]
        self.assertIsNotNone(vert.const)
        V_inst = instantaneous_vertex(fx["S"], fx["C"], fx["nd"], "singlet",
                                      fx["spatial_shape"])
        G2 = calc_g2_dynamic(_ir_compress(fx["green_sc"], axF, fx["nmat"], "green"),
                             fx["beta"])

        def kernel(xp, g2, residency):
            return BondPairKernel(vert, g2, fx["view"], xp=xp,
                                  spatial_shape=fx["spatial_shape"], norb=fx["norb"],
                                  beta=fx["beta"], nfreq=axF.n_freq, V_inst=V_inst,
                                  axF=axF, residency=residency, host_cap=10 ** 12,
                                  device_cap=10 ** 12)

        Kh = kernel(np, G2, "host")
        self.assertEqual(Kh.gap_shape[-1], axF.n_freq)
        rng = np.random.default_rng(13)
        n = int(np.prod(Kh.gap_shape))
        phi = rng.standard_normal(n) + 1j * rng.standard_normal(n)
        ref = Kh.matvec(phi)
        self.assertGreater(np.abs(ref).max(), 0.0)
        Ks = kernel(np, G2, "stream")
        self.assertIsNone(Ks._blocks)                 # rebuilt inside every matvec
        np.testing.assert_allclose(Ks.matvec(phi), ref, atol=1e-13)
        Ka = kernel(np, G2, "auto")
        self.assertEqual(Ka.residency, "host")        # numpy has no separate device
        np.testing.assert_allclose(Ka.matvec(phi), ref, atol=1e-13)
        if not backend.gpu_available():
            raise unittest.SkipTest("CUDA device required for the device residency")
        import cupy
        Kd = kernel(cupy, cupy.asarray(G2), "device")
        self.assertEqual(Kd.residency, "device")
        np.testing.assert_allclose(Kd.matvec(phi), ref, atol=1e-13)


class _NotNumpy(object):
    """An array module that is NOT ``numpy`` and whose every operation IS
    numpy's: it makes the kernel take its device-backend branch on a machine
    with no device, which is the only way to reach the ``device_cap is None``
    path from a CPU test."""

    def __getattr__(self, name):
        return getattr(np, name)


class TestDeviceProbeFallback(unittest.TestCase):
    """``backend.device_available_bytes()`` returns ``None`` when there is no
    device AND when the query itself fails (a busy driver, a container
    without the management interface). A caller that multiplies it loses the
    run to a ``TypeError`` about ``NoneType``; the admission must instead
    fall back to the host cap and say so."""

    def test_helper_falls_back_to_the_host_cap(self):
        from hwave.solver import backend as bk, eliashberg_bond as eb
        with mock.patch.object(bk, "device_available_bytes", return_value=None):
            with self.assertLogs("qlms.eliashberg_bond", "WARNING") as log:
                cap = eb.device_cap_or_host(7.0 * eb._GIB, "probe test")
        self.assertEqual(cap, 7.0 * eb._GIB)
        self.assertIn("device memory probe", "\n".join(log.output))
        self.assertIn("probe test", "\n".join(log.output))
        # a probe that answers is used as before, with the 0.9 headroom
        with mock.patch.object(bk, "device_available_bytes", return_value=1000):
            self.assertEqual(eb.device_cap_or_host(7.0 * eb._GIB, "probe test"), 900.0)

    def test_kernel_admission_survives_a_failed_probe(self):
        """``BondPairKernel(device_cap=None)`` on a device backend: the
        admission runs with the host cap on both sides instead of raising."""
        from hwave.solver import backend as bk, eliashberg_bond as eb
        fx = physical_fixture(norb=1, shape=(4, 2, 1), nmat=4, beta=1.0)
        vert = _uniform_vertex(fx, "singlet")
        seen = {}
        real_choose = eb.AdmissionTable.choose

        def spy(table, requested, host_cap, device_cap):
            seen.update(requested=requested, host=host_cap, device=device_cap)
            return real_choose(table, requested, host_cap, device_cap)

        with mock.patch.object(bk, "device_available_bytes", return_value=None), \
                mock.patch.object(eb.AdmissionTable, "choose", spy):
            with self.assertLogs("qlms.eliashberg_bond", "WARNING") as log:
                K = eb.BondPairKernel(vert, _g2(fx), fx["view"], xp=_NotNumpy(),
                                      spatial_shape=fx["spatial_shape"], norb=fx["norb"],
                                      beta=fx["beta"], nfreq=fx["nmat"], residency="stream",
                                      host_cap=8.0 * eb._GIB, device_cap=None)
        self.assertEqual(K.residency, "stream")
        self.assertEqual(seen["device"], seen["host"])
        self.assertEqual(seen["host"], 8.0 * eb._GIB)
        self.assertIn("device memory probe", "\n".join(log.output))


if __name__ == "__main__":
    unittest.main()
