"""The coeff_tail equal-time endpoint fix (issue #134).

The tail piece aa/(i w_n) carries the Green function's equal-time jump:
its tau branches are -aa/2 (0 < tau < beta) and +aa/2 (tau < 0). The
chi0 kernels reconstructed the 0^+ branch at every tau slot, but the
REVERSED factor's tau = 0 entry is the 0^- side, so the tail-corrected
bubble carried an O(1) error at one tau point -- summing to an
O(1/Nmat) error ~4.5x LARGER than the one the correction removes. With
the jump restored on the reversed factor, coeff_tail = 1 converges at
O(1/Nmat^2) (measured fourfold-Nmat ratios ~16) instead of degrading the
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
        # Extern convention resolved (round-1 review): coeff_extern h
        # gives H_up = H0 + h, H_dn = H0 - h, so e_dn = eps - MU - HZ
        # and e_up = eps - MU + HZ; the fixed static index and the single
        # candidate together also pin the bosonic index convention.
        candidates = [self._exact_pm()]

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
            # static bosonic index is nmat//2 on this kernel's centered
            # output; a drifting index convention must fail here, not be
            # absorbed by an argmax search (round-1 review)
            idx = int(np.argmax(np.abs(pm).sum(axis=1)))
            assert idx == nmat // 2, idx
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


class TestSyntheticTailReversal(unittest.TestCase):
    """Missing-spatial-reversal detector (round-1 review): production
    VV^dag is the identity (complete eigenbasis), so the jump is
    momentum-independent and the committed physical fixtures cannot see
    a dropped reverse_fft_axes on the reversed factor's jump. A
    synthetic momentum-DEPENDENT tail on a length-3 axis discriminates:
    the kernel's tau = 0 slice must equal an independently constructed
    mean of the two equal-time branches."""

    def test_kernel_matches_branch_construction_with_k_dependent_tail(self):
        import hwave.solver.matsubara as _ms
        import hwave.solver.backend as _bk

        NX, NMAT, nd = 3, 4, 1
        nx, nmat = NX, NMAT
        rng = np.random.RandomState(134)
        green = (rng.randn(1, nmat, nx, nd, nd)
                 + 1j * rng.randn(1, nmat, nx, nd, nd))
        # frequency-constant, momentum-dependent synthetic tail
        tailk = (rng.randn(1, 1, nx, nd, nd)
                 + 1j * rng.randn(1, 1, nx, nd, nd))
        tail = np.tile(tailk, (1, nmat, 1, 1, 1))

        # independent branch construction: build G(tau) on both branches
        # explicitly, form chi(tau) with the mean at tau = 0, transform.
        # A singleton FFT axis is the identity, so the flat reference is
        # the same for either orientation of the length-3 axis below.
        gkt = _ms.fermion_to_tau(green.reshape(1, nmat, -1), axis=1
                                 ).reshape(1, nmat, nx, nd, nd) - tail
        grt = _bk.spatial_ifftn(gkt.reshape(1, nmat, nx, 1, 1, nd * nd),
                                axes=(2, 3, 4), workers=1
                                ).reshape(1, nmat, nx, nd, nd)
        jr = _bk.spatial_ifftn(
            (2.0 * tailk).reshape(1, 1, nx, 1, 1, nd * nd),
            axes=(2, 3, 4), workers=1).reshape(1, nx, nd, nd)
        grev = np.stack([grt[:, (-l) % nmat] for l in range(nmat)], axis=1)
        grev = np.stack([grev[:, l][:, [( -r) % nx for r in range(nx)]]
                         for l in range(nmat)], axis=1)
        jrev = jr[:, [(-r) % nx for r in range(nx)]]
        chi = np.empty((1, nmat, nx, nd, nd), dtype=complex)
        for l in range(nmat):
            sgn = 1.0 if l == 0 else -1.0
            chi[:, l] = sgn * grt[:, l] * grev[:, l].swapaxes(-2, -1)
        chi[:, 0] = 0.5 * (grt[:, 0] * (grev[:, 0] + jrev
                                        ).swapaxes(-2, -1)
                           + (grt[:, 0] + jr)
                           * grev[:, 0].swapaxes(-2, -1))
        cqt = _bk.spatial_fftn(chi.reshape(1, nmat, nx, 1, 1, nd * nd),
                               axes=(2, 3, 4), workers=1)
        want = _ms.tau_to_boson(cqt.reshape(1, nmat, -1), axis=1
                                ).reshape(1, nmat, nx, nd, nd) * (-1.0)

        # both orientations of the nontrivial axis: an axis tuple that is
        # off by one (round-2's transverse defect was exactly that, on
        # the z slot) is invisible when only x varies (round-3 review)
        import hwave.solver.rpa as R
        for shape in ((NX, 1, 1), (1, 1, NX)):
            with self.subTest(shape=shape):
                class _Lat:
                    pass
                _Lat.shape = shape
                _Lat.nvol = NX

                class _S:
                    lattice = _Lat()
                    nmat = NMAT
                    enable_reduced = True
                    fft_workers = 1

                got = R.RPA._calc_chi0q(_S(), green, tail, 1.0)
                np.testing.assert_allclose(np.asarray(got), want,
                                           atol=1e-12)


class TestNblockActiveTailReducedGeneral(unittest.TestCase):
    """nblock=2 longitudinal consistency with per-block ACTIVE tails
    (round-1 review): reduced[g, l, q, a, b] must equal the general
    scheme's diagonal slots general[g, l, q, a, a, b, b], including the
    endpoint correction, with each spin block carrying a DIFFERENT
    momentum-dependent tail so a block mix-up or a diagonal-only jump
    cannot cancel."""

    def test_reduced_matches_general_diagonal(self):
        NX, NMAT, ND = 3, 4, 2
        rng = np.random.RandomState(1134)
        green = (rng.randn(2, NMAT, NX, ND, ND)
                 + 1j * rng.randn(2, NMAT, NX, ND, ND))
        tailk = (rng.randn(2, 1, NX, ND, ND)
                 + 1j * rng.randn(2, 1, NX, ND, ND))
        tail = np.tile(tailk, (1, NMAT, 1, 1, 1))

        class _Lat:
            shape = (NX, 1, 1)
            nvol = NX

        def _stub(reduced):
            class _S:
                lattice = _Lat()
                nmat = NMAT
                enable_reduced = reduced
                fft_workers = 1
            return _S()

        import hwave.solver.rpa as R
        red = np.asarray(R.RPA._calc_chi0q(_stub(True), green, tail, 1.0))
        gen = np.asarray(R.RPA._calc_chi0q(_stub(False), green, tail, 1.0))
        self.assertEqual(red.shape, (2, NMAT, NX, ND, ND))
        self.assertEqual(gen.shape, (2, NMAT, NX, ND, ND, ND, ND))
        for a in range(ND):
            for b in range(ND):
                np.testing.assert_allclose(
                    red[..., a, b], gen[..., a, a, b, b], atol=1e-13)


class TestTransverseSyntheticTail(unittest.TestCase):
    """Transverse analogues of the synthetic-tail checks (round-2
    review): every committed active-tail transverse fixture has nz = 1,
    so a reversal applied to the wrong axis tuple -- e.g. one that
    still carries the spin axis and therefore reverses (tau, x, y)
    instead of (x, y, z) -- was invisible. A momentum-dependent
    down-spin tail on a (1, 1, 3) lattice pins the exact modular
    spatial reversal in both schemes, against a reference that indexes
    (-r) mod n explicitly."""

    NX, NMAT, ND = 3, 4, 2

    def _fixture(self):
        rng = np.random.RandomState(2134)
        green = (rng.randn(2, self.NMAT, self.NX, self.ND, self.ND)
                 + 1j * rng.randn(2, self.NMAT, self.NX, self.ND, self.ND))
        tailk = (rng.randn(2, 1, self.NX, self.ND, self.ND)
                 + 1j * rng.randn(2, 1, self.NX, self.ND, self.ND))
        tail = np.tile(tailk, (1, self.NMAT, 1, 1, 1))
        return green, tailk, tail

    def _stub(self, reduced):
        NX, NMAT = self.NX, self.NMAT

        class _Lat:
            shape = (1, 1, NX)
            nvol = NX

        class _S:
            lattice = _Lat()
            nmat = NMAT
            enable_reduced = reduced
            fft_workers = 1
        return _S()

    def _reference(self, green, tailk, reduced):
        import hwave.solver.matsubara as _ms
        import hwave.solver.backend as _bk
        NX, NMAT, ND = self.NX, self.NMAT, self.ND
        rev = [(-r) % NX for r in range(NX)]
        gkt = _ms.fermion_to_tau(green.reshape(2, NMAT, -1), axis=1
                                 ).reshape(2, NMAT, NX, ND, ND) \
            - np.tile(tailk, (1, NMAT, 1, 1, 1))
        grt = _bk.spatial_ifftn(
            gkt.reshape(2, NMAT, 1, 1, NX, ND * ND),
            axes=(2, 3, 4), workers=1).reshape(2, NMAT, NX, ND, ND)
        up = grt[0]
        dn_rev = np.stack([grt[1][(-l) % NMAT][rev] for l in range(NMAT)])
        jr = _bk.spatial_ifftn(
            (2.0 * tailk).reshape(2, 1, 1, 1, NX, ND * ND),
            axes=(2, 3, 4), workers=1).reshape(2, NX, ND, ND)
        jump_up, jump_dn_rev = jr[0], jr[1][rev]
        if reduced:
            chi = np.stack([(1.0 if l == 0 else -1.0)
                            * up[l] * dn_rev[l].swapaxes(-2, -1)
                            for l in range(NMAT)])
            chi[0] = 0.5 * (
                up[0] * (dn_rev[0] + jump_dn_rev).swapaxes(-2, -1)
                + (up[0] + jump_up) * dn_rev[0].swapaxes(-2, -1))
            nds = ND * ND
            tshape = (NMAT, NX, ND, ND)
        else:
            chi = np.stack([(1.0 if l == 0 else -1.0)
                            * np.einsum('rab,rdc->racbd', up[l], dn_rev[l])
                            for l in range(NMAT)])
            chi[0] = 0.5 * (
                np.einsum('rab,rdc->racbd', up[0],
                          dn_rev[0] + jump_dn_rev)
                + np.einsum('rab,rdc->racbd', up[0] + jump_up, dn_rev[0]))
            nds = ND ** 4
            tshape = (NMAT, NX, ND, ND, ND, ND)
        cqt = _bk.spatial_fftn(chi.reshape(NMAT, 1, 1, NX, nds),
                               axes=(1, 2, 3), workers=1)
        return _ms.tau_to_boson(
            cqt.reshape(1, NMAT, -1), axis=1).reshape(*tshape) * (-1.0)

    def test_kernel_matches_branch_construction(self):
        import hwave.solver.rpa as R
        green, tailk, tail = self._fixture()
        for reduced in (True, False):
            with self.subTest(reduced=reduced):
                got = np.asarray(R.RPA._calc_chi0q_transverse(
                    self._stub(reduced), green, tail, 1.0))
                want = self._reference(green, tailk, reduced)
                np.testing.assert_allclose(got, want, atol=1e-12)

    def test_reduced_matches_general_diagonal(self):
        """Distinct active up/down tails: reduced[l, q, a, b] must equal
        general[l, q, a, a, b, b] (the longitudinal analogue already
        exists; the transverse kernels are separate code)."""
        import hwave.solver.rpa as R
        green, _, tail = self._fixture()
        red = np.asarray(R.RPA._calc_chi0q_transverse(
            self._stub(True), green, tail, 1.0))
        gen = np.asarray(R.RPA._calc_chi0q_transverse(
            self._stub(False), green, tail, 1.0))
        for a in range(self.ND):
            for b in range(self.ND):
                np.testing.assert_allclose(
                    red[..., a, b], gen[..., a, a, b, b], atol=1e-13)


class TestZeroTailBitwise(unittest.TestCase):
    """coeff_tail = 0 must be BITWISE unchanged by the endpoint fix
    (round-3 review): the spy test proves the correction FFT is
    skipped, this one asserts the stated output guarantee directly.
    The reference replicates the pre-fix zero-tail operation sequence
    verbatim (same ops, same order, same process), so equality must be
    exact -- assert_array_equal, no tolerance."""

    NX, NMAT, ND = 3, 4, 2
    BETA = 0.7

    def _lat_stub(self, reduced):
        NX, NMAT = self.NX, self.NMAT

        class _Lat:
            shape = (1, 1, NX)
            nvol = NX

        class _S:
            lattice = _Lat()
            nmat = NMAT
            enable_reduced = reduced
            fft_workers = 1
        return _S()

    def _green(self, nblock, seed):
        rng = np.random.RandomState(seed)
        return (rng.randn(nblock, self.NMAT, self.NX, self.ND, self.ND)
                + 1j * rng.randn(nblock, self.NMAT, self.NX, self.ND,
                                 self.ND))

    def test_longitudinal_forms(self):
        import hwave.solver.rpa as R
        import hwave.solver.matsubara as _ms
        import hwave.solver.backend as _bk
        from hwave.solver.kgrid import reverse_fft_axes
        NX, NMAT, ND, BETA = self.NX, self.NMAT, self.ND, self.BETA
        green = self._green(2, 34)
        zeros = np.zeros_like(green)
        for reduced in (True, False):
            with self.subTest(reduced=reduced):
                got = np.asarray(R.RPA._calc_chi0q(
                    self._lat_stub(reduced), green, zeros, BETA))
                gkt = _ms.fermion_to_tau(
                    green.reshape(2, NMAT, -1), axis=1).reshape(
                    2, NMAT, 1, 1, NX, ND, ND)
                gkt -= zeros.reshape(2, NMAT, 1, 1, NX, ND, ND)
                grt = _bk.spatial_ifftn(
                    gkt.reshape(2, NMAT, 1, 1, NX, ND * ND),
                    axes=(2, 3, 4), workers=1)
                grev = reverse_fft_axes(grt, (1, 2, 3, 4)).reshape(
                    2, NMAT, NX, ND, ND)
                sgn = np.full(NMAT, -1)
                sgn[0] = 1
                if reduced:
                    g5 = grt.reshape(2, NMAT, NX, ND, ND)
                    chi = (g5 * grev.swapaxes(-2, -1)
                           * sgn[np.newaxis, :, np.newaxis, np.newaxis,
                                 np.newaxis])
                    nd_shape, nds = (ND, ND), ND ** 2
                else:
                    gf = grt.reshape(2, NMAT, NX, ND, ND)
                    sgn_bc = sgn[np.newaxis, :, np.newaxis, np.newaxis,
                                 np.newaxis]
                    chi = ((gf * sgn_bc)[:, :, :, :, np.newaxis, :,
                                         np.newaxis]
                           * grev[:, :, :, np.newaxis, :, np.newaxis, :])
                    chi = chi.transpose(0, 1, 2, 3, 6, 5, 4)
                    nd_shape, nds = (ND, ND, ND, ND), ND ** 4
                cqt = _bk.spatial_fftn(
                    chi.reshape(2, NMAT, 1, 1, NX, nds),
                    axes=(2, 3, 4), workers=1)
                want = _ms.tau_to_boson(
                    cqt.reshape(2, NMAT, NX * nds), axis=1).reshape(
                    2, NMAT, NX, *nd_shape) * (-1.0 / BETA)
                np.testing.assert_array_equal(got, want)

    def test_transverse_forms(self):
        import hwave.solver.rpa as R
        import hwave.solver.backend as _bk
        from hwave.solver.kgrid import reverse_fft_axes
        NX, NMAT, ND, BETA = self.NX, self.NMAT, self.ND, self.BETA
        green = self._green(2, 35)
        zeros = np.zeros_like(green)
        for reduced in (True, False):
            with self.subTest(reduced=reduced):
                got = np.asarray(R.RPA._calc_chi0q_transverse(
                    self._lat_stub(reduced), green, zeros, BETA))
                omg = np.exp(-1j * np.pi * (1.0 / NMAT - 1.0)
                             * np.arange(NMAT))
                gkt = (np.fft.fft(green.reshape(2, NMAT, -1), axis=1)
                       * omg[np.newaxis, :, np.newaxis]).reshape(
                    2, NMAT, 1, 1, NX, ND, ND)
                gkt -= zeros.reshape(2, NMAT, 1, 1, NX, ND, ND)
                grt = _bk.spatial_ifftn(
                    gkt.reshape(2, NMAT, 1, 1, NX, ND * ND),
                    axes=(2, 3, 4), workers=1)
                dn_rev = reverse_fft_axes(grt[1], (0, 1, 2, 3)).reshape(
                    NMAT, NX, ND, ND)
                up = grt[0].reshape(NMAT, NX, ND, ND)
                sgn = np.full(NMAT, -1)
                sgn[0] = 1
                if reduced:
                    chi = (up * dn_rev.swapaxes(-2, -1)
                           * sgn[:, np.newaxis, np.newaxis, np.newaxis])
                    nd_shape, nds = (ND, ND), ND ** 2
                else:
                    chi = np.einsum('lrab,lrdc,l->lracbd', up, dn_rev,
                                    sgn)
                    nd_shape, nds = (ND, ND, ND, ND), ND ** 4
                cqt = _bk.spatial_fftn(
                    chi.reshape(NMAT, 1, 1, NX, nds),
                    axes=(1, 2, 3), workers=1)
                omg2 = np.exp(1j * np.pi * (-1) * np.arange(NMAT))
                want = np.fft.ifft(
                    cqt.reshape(NMAT, NX * nds) * omg2[:, np.newaxis],
                    axis=0).reshape(NMAT, NX, *nd_shape) * (-1.0 / BETA)
                np.testing.assert_array_equal(got, want)


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
