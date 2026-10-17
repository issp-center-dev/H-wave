"""The documented Fourier sign, pinned per path (issue #133).

The documented convention (docs/en/source/algorithm/rpa.rst, Wannier90
style) is M(k) = sum_R M(R) e^{+ikR}. Before #133 the non-spin-orbital
RPA transfer, Extern, trans_mod, and _calc_trans_mod paths -- and the
interaction/Fierz FFTs shared by BOTH modes -- used e^{-ikR}: the
non-SO mode was a self-consistent global k -> -k relabeling of the
documentation, while the SO mode combined a e^{+ikR} transfer with a
e^{-iqR} interaction, i.e. chiq was solved with chi0(q) against W(-q).
For R-symmetric real interactions W(q) = W(-q) and nothing was visible;
these tests use non-centrosymmetric complex fixtures so every sign is
load-bearing, and the chiq test below FAILS against the pre-#133 code.
"""

import os
import tempfile
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.rpa as solver_rpa


LX = 4
T_HOP = 0.7 + 0.3j       # non-centrosymmetric complex hopping
V_DIR = 0.7              # DIRECTIONAL inter-orbital bond V_{01}(R=+1):
                         # after the adjudicated closure/symmetrisation the
                         # table keeps X_{01}(+1) = X_{10}(-1) = V_DIR with
                         # no X_{01}(-1), so W(q) has e^{+iq} structure and
                         # W(-q) = W(q)^T != W(q) -- a legal file input for
                         # which the Fourier sign is load-bearing. (A
                         # same-orbital complex pair is NOT usable here: the
                         # symmetrised reading folds its phase to a real
                         # coefficient and W(q) degenerates to a cosine.)


def _write_chain(d, so, norb_phys):
    """norb_phys physical orbitals; in SO mode nd = 2*norb_phys with an
    on-site spin flip so the Hamiltonian classifies as spinful. The
    interaction file always uses PHYSICAL orbital indices."""
    nd = 2 * norb_phys if so else norb_phys
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n%d\n" % nd)
        for _ in range(nd):
            f.write("0.0 0.0 0.0\n")
    ents = []
    for i in range(nd):
        ents.append((1, i, i, T_HOP * (1 + 0.2 * i)))
        ents.append((-1, i, i, np.conj(T_HOP * (1 + 0.2 * i))))
    if nd > 1:
        # complex inter-index hopping: q-asymmetric chi0
        ents.append((0, 0, 1, 0.3 + 0.1j))
        ents.append((0, 1, 0, 0.3 - 0.1j))
        ents.append((1, 0, 1, 0.15 + 0.25j))
        ents.append((-1, 1, 0, 0.15 - 0.25j))
    rv = sorted({r for r, _, _, _ in ents})
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("hdr\n%d\n%d\n%s\n" % (nd, len(rv), " ".join("1" * len(rv))))
        for r, i, j, v in ents:
            f.write("%3d 0 0 %d %d %.12f %.12f\n"
                    % (r, i + 1, j + 1, v.real, v.imag))
    with open(os.path.join(d, "coulombinter.dat"), "w") as f:
        f.write("hdr\n%d\n2\n1 1\n" % norb_phys)
        f.write(" 1 0 0 1 2 %.12f 0.0\n" % V_DIR)
        f.write("-1 0 0 2 1 %.12f 0.0\n" % V_DIR)


def _solver(d, so, scheme="reduced", with_v=True, nmat=64):
    inter = {"path_to_input": d, "Geometry": "geom.dat",
             "Transfer": "transfer.dat"}
    if with_v:
        inter["CoulombInter"] = "coulombinter.dat"
    info_mode = {"mode": "RPA",
                 "param": {"T": 0.5, "mu": 0.1, "CellShape": [LX, 1, 1],
                           "SubShape": [1, 1, 1], "Nmat": nmat},
                 "enable_spin_orbital": so, "calc_scheme": scheme}
    io = read_input_k.QLMSkInput({"path_to_input": d, "interaction": inter})
    solver = solver_rpa.RPA(io.get_param("ham"), {}, info_mode)
    return solver, io.get_param("green")


def _eps_doc(k, i):
    t = T_HOP * (1 + 0.2 * i)
    return t * np.exp(2j * np.pi * k / LX) \
        + np.conj(t) * np.exp(-2j * np.pi * k / LX)


class TestTransferSign(unittest.TestCase):

    def test_non_so_transfer_is_documented_sign(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        _write_chain(tmp.name, so=False, norb_phys=1)
        solver, _ = _solver(tmp.name, so=False, with_v=False)
        eps = np.asarray(solver.ham_info.ham_trans_q).reshape(LX)
        want = np.array([_eps_doc(k, 0) for k in range(LX)])
        np.testing.assert_allclose(eps, want, atol=1e-12)

    def test_so_transfer_is_documented_sign(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        _write_chain(tmp.name, so=True, norb_phys=1)
        solver, _ = _solver(tmp.name, so=True, with_v=False)
        eps = np.asarray(solver.ham_info.ham_trans_q).reshape(LX, 2, 2)
        for i in (0, 1):
            want = np.array([_eps_doc(k, i) for k in range(LX)])
            np.testing.assert_allclose(eps[:, i, i], want, atol=1e-12)


class TestInteractionSign(unittest.TestCase):

    def test_w_q_is_documented_sign(self):
        """The density slots of the assembled vertex must carry
        V_ab(q) = sum_R X_ab(R) e^{+iqR}, not V_ab(-q). The directional
        bond gives V_01(q) = V e^{+iq} and V_10(q) = V e^{-iq}, which a
        flipped sign would swap."""
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        _write_chain(tmp.name, so=False, norb_phys=2)
        solver, _ = _solver(tmp.name, so=False, scheme="reduced")
        ham = np.asarray(solver.ham_info.ham_inter_q).reshape(
            LX, 2, 2, 2, 2, 2, 2, 2, 2)
        # _append_inter places the entry's orbital pair as
        # (s4,b, s3,b, s1,a, s2,a): the FIRST pair carries b, the second a
        q = 2.0 * np.pi * np.arange(LX) / LX
        w01 = ham[:, 0, 1, 0, 1, 1, 0, 1, 0]     # a=0, b=1: V_01(q)
        w10 = ham[:, 0, 0, 0, 0, 1, 1, 1, 1]     # a=1, b=0: V_10(q)
        np.testing.assert_allclose(w01, V_DIR * np.exp(+1j * q), atol=1e-12)
        np.testing.assert_allclose(w10, V_DIR * np.exp(-1j * q), atol=1e-12)


class TestSpinOrbitalChiqConsistency(unittest.TestCase):

    def test_so_chiq_combines_chi0_and_w_at_the_same_q(self):
        """The mixed-convention catcher: before #133 the SO mode solved
        [1 + chi0(q) W(-q)]^-1 chi0(q). Build both candidates from the
        solver's OWN chi0q and the analytic W and pin the consistent one;
        this test fails against the pre-#133 code."""
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        _write_chain(tmp.name, so=True, norb_phys=2)
        solver, gi = _solver(tmp.name, so=True, scheme="reduced")
        out = os.path.join(tmp.name, "out")
        os.makedirs(out, exist_ok=True)
        solver.solve(gi, out)
        chi0 = np.asarray(gi["chi0q"])     # (nfreq, nvol, nd, nd), spinful
        chiq = np.asarray(gi["chiq"])
        nfreq = chi0.shape[0]
        nd = chi0.shape[-1]
        # non-vacuous premises: q-asymmetric chi0 AND W(q) != W(-q)
        il0 = nfreq // 2
        self.assertGreater(
            max(np.abs(chi0[il0, q] - chi0[il0, (-q) % LX]).max()
                for q in range(LX)), 1e-4)
        # ANALYTIC W in the documented sign, independent of the solver
        # (using the solver's own ham_inter_q would be self-consistent in
        # any version and catch nothing): density projection of the
        # CoulombInter vertex in spin-block space is
        # w[(s,alpha),(t,beta)] = V(q)^T[alpha, beta] for every spin pair
        # (the paper's Eq.(21) pair transpose; weight 1 for all spins).
        norb_phys = nd // 2
        w_analytic = np.zeros((LX, nd, nd), dtype=complex)
        qs = 2.0 * np.pi * np.arange(LX) / LX
        for iq, q in enumerate(qs):
            vq = np.zeros((norb_phys, norb_phys), dtype=complex)
            vq[0, 1] = V_DIR * np.exp(+1j * q)
            vq[1, 0] = V_DIR * np.exp(-1j * q)
            for st in range(2):
                for tt in range(2):
                    w_analytic[iq,
                               st * norb_phys:(st + 1) * norb_phys,
                               tt * norb_phys:(tt + 1) * norb_phys] = vq.T
        self.assertGreater(
            max(np.abs(w_analytic[q] - w_analytic[(-q) % LX]).max()
                for q in range(LX)), 0.1)
        matches = {}
        for sign in (+1, -1):
            ok = True
            for il in range(nfreq):
                for q in range(LX):
                    qq = q if sign > 0 else (-q) % LX
                    m = np.eye(nd) + chi0[il, q] @ w_analytic[qq]
                    manual = np.linalg.solve(m, chi0[il, q])
                    if not np.allclose(manual, chiq[il, q], atol=1e-10):
                        ok = False
            matches[sign] = ok
        self.assertTrue(matches[+1],
                        "chiq must match the documented-sign W(q) assembly")
        self.assertFalse(matches[-1],
                         "chiq must NOT match the W(-q) assembly")


if __name__ == "__main__":
    unittest.main()
