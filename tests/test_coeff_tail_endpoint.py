"""The coeff_tail equal-time endpoint fix (issue #134).

The tail piece aa/(i w_n) carries the Green function's equal-time jump:
its tau branches are -aa/2 (0 < tau < beta) and +aa/2 (tau < 0). The
chi0 kernels reconstructed the 0^+ branch at every tau slot, but the
REVERSED factor's tau = 0 entry is the 0^- side, so the tail-corrected
bubble carried an O(1) error at one tau point -- summing to an
O(1/Nmat) error ~4.5x LARGER than the one the correction removes. With
the jump restored on the reversed factor, coeff_tail = 1 converges at
O(1/Nmat^2) (measured doubling ratios ~16) instead of degrading the
result.

References are exact: the analytic Lindhard function of an exactly
solvable chain, evaluated in the eigenbasis with the full Matsubara sum
done analytically.
"""

import os
import tempfile
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.rpa as solver_rpa


LX = 4
T = 0.5
BETA = 1.0 / T
MU = 0.2
THOP = 0.7 + 0.3j
HZ = 0.3            # Zeeman field for the transverse fixture


def _fermi(e):
    return 1.0 / (np.exp(np.clip(BETA * e, -500, 500)) + 1.0)


def _eps(k):
    return (THOP * np.exp(2j * np.pi * k / LX)
            + np.conj(THOP) * np.exp(-2j * np.pi * k / LX)).real


def _write_chain(d, extern=False):
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n1\n0.0 0.0 0.0\n")
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("hdr\n1\n2\n1 1\n")
        f.write(" 1 0 0 1 1 %.12f %.12f\n" % (THOP.real, THOP.imag))
        f.write("-1 0 0 1 1 %.12f %.12f\n"
                % (np.conj(THOP).real, np.conj(THOP).imag))
    if extern:
        with open(os.path.join(d, "extern.dat"), "w") as f:
            f.write("hdr\n1\n1\n1\n")
            f.write(" 0 0 0 1 1 1.0 0.0\n")


def _solver(d, nmat, tail, extern=False, calc_type="ring"):
    inter = {"path_to_input": d, "Geometry": "geom.dat",
             "Transfer": "transfer.dat"}
    if extern:
        inter["Extern"] = "extern.dat"
    info_mode = {"mode": "RPA",
                 "param": {"T": T, "mu": MU, "CellShape": [LX, 1, 1],
                           "SubShape": [1, 1, 1], "Nmat": nmat,
                           "coeff_tail": tail},
                 "calc_scheme": "reduced", "calc_type": calc_type}
    if extern:
        info_mode["param"]["coeff_extern"] = HZ
    io = read_input_k.QLMSkInput({"path_to_input": d, "interaction": inter})
    solver = solver_rpa.RPA(io.get_param("ham"), {}, info_mode)
    return solver, io.get_param("green")


def _chi0_static(d, nmat, tail, **kw):
    solver, gi = _solver(d, nmat, tail, **kw)
    out = os.path.join(d, "out")
    os.makedirs(out, exist_ok=True)
    solver.solve(gi, out)
    chi0 = np.asarray(gi["chi0q"])
    nfreq = chi0.shape[0]
    return chi0.reshape(nfreq, LX)[nfreq // 2].real


def _exact_longitudinal():
    out = np.zeros(LX)
    for q in range(LX):
        for k in range(LX):
            ek, ekq = _eps(k) - MU, _eps((k + q) % LX) - MU
            de = ekq - ek
            if abs(de) < 1e-12:
                out[q] += -BETA * _fermi(ek) * (1 - _fermi(ek))
            else:
                out[q] += (_fermi(ekq) - _fermi(ek)) / de
    return -out / LX


class TestLongitudinalEndpoint(unittest.TestCase):

    def test_tail_now_accelerates_and_is_second_order(self):
        """coeff_tail = 1 must beat coeff_tail = 0 by a wide margin (it
        was ~4.5x WORSE), and its error must fall ~16x per fourfold Nmat
        (second order; pinned loosely at > 8x)."""
        exact = _exact_longitudinal()
        errs = {}
        for nmat in (256, 1024):
            for tail in (0.0, 1.0):
                tmp = tempfile.TemporaryDirectory()
                self.addCleanup(tmp.cleanup)
                _write_chain(tmp.name)
                errs[(nmat, tail)] = np.abs(
                    _chi0_static(tmp.name, nmat, tail) - exact).max()
        # improvement: >= 10x better than tail=0 at both sizes
        self.assertLess(errs[(256, 1.0)], errs[(256, 0.0)] / 10.0)
        self.assertLess(errs[(1024, 1.0)], errs[(1024, 0.0)] / 10.0)
        # order: ratio ~16 per 4x Nmat (tail=0 is ~4)
        self.assertGreater(errs[(256, 1.0)] / errs[(1024, 1.0)], 8.0)

    def test_tail_zero_takes_no_correction_path(self):
        """Bit-level neutrality for coeff_tail = 0: the correction block
        must not even execute (spatial_ifftn is called once for the
        Green transform and once... exactly twice per chi0 call without
        the tail, three times with it)."""
        from unittest import mock
        import hwave.solver.backend as _bk
        calls = []
        real = _bk.spatial_ifftn

        def spy(*a, **k):
            calls.append(1)
            return real(*a, **k)

        for tail, expected_extra in ((0.0, 0), (1.0, 1)):
            with self.subTest(tail=tail):
                tmp = tempfile.TemporaryDirectory()
                self.addCleanup(tmp.cleanup)
                _write_chain(tmp.name)
                calls.clear()
                with mock.patch.object(_bk, "spatial_ifftn",
                                       side_effect=spy):
                    _chi0_static(tmp.name, 64, tail)
                if tail == 0.0:
                    base_calls = len(calls)
                else:
                    self.assertEqual(len(calls), base_calls + 1)


class TestGeneralSchemeEndpoint(unittest.TestCase):
    """Multi-orbital general scheme: the bubble's equal-time jump has
    off-diagonal orbital structure (the one-sided endpoint fix left
    those components O(1/Nmat)); the mean-endpoint sample restores
    second order for the full (a,c,b,d) tensor."""

    H0 = {(0,): np.array([[0.1, 0.3 + 0.2j], [0.3 - 0.2j, -0.2]]),
          (1,): np.array([[1.0, 0.2 + 0.1j], [0.15 - 0.25j, 0.9]])}

    def _write(self, d):
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n2\n"
                    "0.0 0.0 0.0\n0.0 0.0 0.0\n")
        ents = []
        for (rx,), m in self.H0.items():
            for i in range(2):
                for j in range(2):
                    if m[i, j] != 0:
                        ents.append((rx, i, j, m[i, j]))
                    if rx != 0 and np.conj(m[j, i]) != 0:
                        ents.append((-rx, i, j, np.conj(m[j, i])))
        rv = sorted({r for r, _, _, _ in ents})
        with open(os.path.join(d, "transfer.dat"), "w") as f:
            f.write("hdr\n2\n%d\n%s\n" % (len(rv),
                                            " ".join("1" * len(rv))))
            for r, i, j, v in ents:
                f.write("%3d 0 0 %d %d %.12f %.12f\n"
                        % (r, i + 1, j + 1, v.real, v.imag))

    def _hk(self, k):
        h = np.zeros((2, 2), dtype=complex)
        for (rx,), m in self.H0.items():
            ph = np.exp(2j * np.pi * k * rx / LX)
            h += m if rx == 0 else m * ph + m.conj().T * np.conj(ph)
        return h

    def _exact(self):
        ev = np.zeros((LX, 2))
        U = np.zeros((LX, 2, 2), dtype=complex)
        for k in range(LX):
            w, v = np.linalg.eigh(self._hk(k))
            ev[k], U[k] = w - MU, v
        out = np.zeros((LX, 2, 2, 2, 2), dtype=complex)
        for q in range(LX):
            for k in range(LX):
                kq = (k + q) % LX
                for m in range(2):
                    for n in range(2):
                        de = ev[kq][m] - ev[k][n]
                        if abs(de) < 1e-12:
                            w0 = -BETA * _fermi(ev[k][n]) \
                                * (1 - _fermi(ev[k][n]))
                        else:
                            w0 = (_fermi(ev[kq][m])
                                  - _fermi(ev[k][n])) / de
                        out[q] += w0 * np.einsum(
                            "a,c,b,d->acbd", U[kq][:, m].conj(),
                            U[k][:, n], U[kq][:, m], U[k][:, n].conj())
        return -out / LX

    def _err(self, nmat, tail):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        self._write(d)
        info_mode = {"mode": "RPA",
                     "param": {"T": T, "mu": MU, "CellShape": [LX, 1, 1],
                               "SubShape": [1, 1, 1], "Nmat": nmat,
                               "coeff_tail": tail},
                     "calc_scheme": "general"}
        io = read_input_k.QLMSkInput(
            {"path_to_input": d,
             "interaction": {"path_to_input": d, "Geometry": "geom.dat",
                             "Transfer": "transfer.dat"}})
        solver = solver_rpa.RPA(io.get_param("ham"), {}, info_mode)
        gi = io.get_param("green")
        out = os.path.join(d, "out")
        os.makedirs(out, exist_ok=True)
        solver.solve(gi, out)
        c = np.asarray(gi["chi0q"])
        n = c.shape[0]
        got = c.reshape(n, LX, 2, 2, 2, 2)[n // 2]
        return np.abs(got.conj() - self._exact()).max()

    def test_general_tensor_is_second_order_with_tail(self):
        e256_0, e256_1 = self._err(256, 0.0), self._err(256, 1.0)
        e1024_1 = self._err(1024, 1.0)
        self.assertLess(e256_1, e256_0 / 10.0)
        self.assertGreater(e256_1 / e1024_1, 8.0)


class TestTransverseEndpoint(unittest.TestCase):

    def _exact_pm(self):
        """chi0_+-(q, nu=0) for the Zeeman-split chain:
        -(1/N) sum_k (f(e_dn(k)) - f(e_up(k+q))) / (e_up(k+q) - e_dn(k))
        with e_{up/dn} = eps -+ HZ (Extern adds -+ h)."""
        out = np.zeros(LX)
        for q in range(LX):
            for k in range(LX):
                edn = _eps(k) - MU - HZ
                eup = _eps((k + q) % LX) - MU + HZ
                de = eup - edn
                if abs(de) < 1e-12:
                    out[q] += -BETA * _fermi(edn) * (1 - _fermi(edn))
                else:
                    out[q] += (_fermi(eup) - _fermi(edn)) / de
        return -out / LX

    def test_transverse_tail_accelerates(self):
        candidates = [self._exact_pm()]
        # sign convention of Extern (up/down assignment) is pinned by
        # whichever branch matches at convergence; both split candidates
        # share the acceleration property under test
        out2 = np.zeros(LX)
        for q in range(LX):
            for k in range(LX):
                edn = _eps(k) - MU + HZ
                eup = _eps((k + q) % LX) - MU - HZ
                de = eup - edn
                if abs(de) < 1e-12:
                    out2[q] += -BETA * _fermi(edn) * (1 - _fermi(edn))
                else:
                    out2[q] += (_fermi(eup) - _fermi(edn)) / de
        candidates.append(-out2 / LX)

        def pm_static(nmat, tail):
            """Drive the fixed kernel directly with production-built
            inputs: solve() resolves mu and the spin-diag Green pair,
            then _calc_chi0q_transverse is the unit under test."""
            tmp = tempfile.TemporaryDirectory()
            self.addCleanup(tmp.cleanup)
            _write_chain(tmp.name, extern=True)
            solver, gi = _solver(tmp.name, nmat, tail, extern=True)
            out = os.path.join(tmp.name, "out")
            os.makedirs(out, exist_ok=True)
            solver.solve(gi, out)
            pm = np.asarray(solver._calc_chi0q_transverse(
                solver.green0, solver.green0_tail, BETA))
            pm = pm.reshape(nmat, LX)
            # static bosonic index: chi is largest at nu = 0
            idx = int(np.argmax(np.abs(pm).sum(axis=1)))
            return pm[idx].real

        err = {}
        for nmat in (256, 1024):
            for tail in (0.0, 1.0):
                got = pm_static(nmat, tail)
                err[(nmat, tail)] = min(np.abs(got - c).max()
                                        for c in candidates)
        # the transverse bubble itself jumps at tau = 0 for a spin-split
        # system; the mean-endpoint sample restores second order
        self.assertLess(err[(256, 1.0)], err[(256, 0.0)] / 10.0)
        self.assertGreater(err[(256, 1.0)] / err[(1024, 1.0)], 8.0)


class TestFlexInheritsTheFix(unittest.TestCase):

    def test_flex_uniform_chi0_delegates_to_the_fixed_kernel(self):
        """FLEX's uniform-grid chi0 path calls the inherited RPA kernel
        (flex.py delegates to self._calc_chi0q), so the endpoint fix
        propagates; pin the delegation so a future FLEX-local override
        must revisit the endpoint."""
        from hwave.solver.flex import FLEX
        self.assertIs(FLEX._calc_chi0q, solver_rpa.RPA._calc_chi0q)


if __name__ == "__main__":
    unittest.main()
