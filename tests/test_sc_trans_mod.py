"""hwave_sc on a mean-field one-body band (issue #273).

``[file.input] trans_mod`` (the T + Sigma_HF that UHFk writes), ``green_init``
(T + Sigma_HF[G], UHFk's Hartree + Fock mean field of G, assembled by RPA; issue
#278) and an ``Extern`` field replace the bare Transfer band in
the RPA H0 assembly that the internal chi0q is built on. hwave_sc used to
take the chemical potential and every kernel Green function from the BARE
Transfer band anyway, so chi0q was evaluated at the wrong filling and the
pairing bubble on the wrong band: silently wrong eigenvalues, or unrelated
shape errors for spin-dependent fields.

Oracle: the band a trans_mod defines is an ordinary one-body Hamiltonian, so
a run whose Transfer file IS that band must give the same eigenvalues.

Tests must be run from the repository root.
"""

import logging
import os
import tempfile
import unittest

import numpy as np

import hwave.sc as sc

NX, NY = 8, 8
T_RUN = 0.1
NMAT = 64
FILLING = 0.4375


# ---------------------------------------------------------------------------
# model files
# ---------------------------------------------------------------------------

def _fold(i, n):
    return i if i <= n // 2 - (1 if n % 2 == 0 else 0) else i - n


def _bare_hr(norb):
    """Real-space one-body table h[ix, iy, a, b] of the bare model."""
    h = np.zeros((NX, NY, norb, norb), dtype=complex)
    for a in range(norb):
        for (rx, ry) in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            h[rx % NX, ry % NY, a, a] += -1.0
        for (rx, ry) in ((1, 1), (-1, -1), (1, -1), (-1, 1)):
            h[rx % NX, ry % NY, a, a] += 0.2
    if norb == 2:
        h[0, 0, 0, 1] += 0.3 + 0.2j
        h[0, 0, 1, 0] += 0.3 - 0.2j
        h[0, 0, 1, 1] += 0.4
    return h


def _nematic_hr(norb):
    """Bare table plus a spin-free nematic mean field: x/y anisotropic
    nearest-neighbour shift, an orbital-dependent on-site shift and (two
    orbitals) a complex inter-orbital bond term -- complex so a transposed
    orbital layout would change the band's eigenvectors and G."""
    h = _bare_hr(norb).copy()
    for a in range(norb):
        for (rx, ry) in ((1, 0), (-1, 0)):
            h[rx % NX, ry % NY, a, a] += -0.15
        for (rx, ry) in ((0, 1), (0, -1)):
            h[rx % NX, ry % NY, a, a] += 0.15
        h[0, 0, a, a] += 0.1 * (1 - 1.5 * a)
        # complex next-nearest-neighbour bond: eps(k) != eps(-k), so the
        # one-orbital oracle is not blind to a k -> -k (Fourier-sign) error
        h[1 % NX, 1 % NY, a, a] += 0.1j
        h[(-1) % NX, (-1) % NY, a, a] += -0.1j
    if norb == 2:
        v = 0.05 + 0.08j
        h[1, 0, 0, 1] += v
        h[NX - 1, 0, 1, 0] += np.conj(v)
    # Hermiticity of the real-space table: h(R)_{ab} = conj(h(-R)_{ba})
    hm = np.conj(np.transpose(h[(-np.arange(NX)) % NX][:, (-np.arange(NY)) % NY],
                              (0, 1, 3, 2)))
    assert np.allclose(h, hm)
    return h


def _trans_mod_from_hr(h_up, h_dn=None, mixing=None):
    """UHFk's trans_mod layout: (cellvol, 2*norb, 2*norb), spin-block
    orbital index s*norb + a, cell index C-ordered over (Nx, Ny, Nz)."""
    norb = h_up.shape[-1]
    if h_dn is None:
        h_dn = h_up
    tm = np.zeros((NX * NY, 2 * norb, 2 * norb), dtype=complex)
    tm[:, :norb, :norb] = h_up.reshape(NX * NY, norb, norb)
    tm[:, norb:, norb:] = h_dn.reshape(NX * NY, norb, norb)
    if mixing is not None:
        tm[0, :norb, norb:] += mixing
        tm[0, norb:, :norb] += np.conj(mixing).T
    return tm


def _write_w90(path, title, h):
    """Write a real-space table as an H-wave wannier90-like file."""
    norb = h.shape[-1]
    rows = []
    for ix in range(NX):
        for iy in range(NY):
            for a in range(norb):
                for b in range(norb):
                    v = h[ix, iy, a, b]
                    if v != 0:
                        rows.append((_fold(ix, NX), _fold(iy, NY), 0, a + 1,
                                     b + 1, float(v.real), float(v.imag)))
    rvecs = sorted({r[:3] for r in rows})
    with open(path, "w") as f:
        f.write(title + "\n{}\n{}\n".format(norb, len(rvecs)))
        for i in range(0, len(rvecs), 15):
            f.write(" ".join("1" for _ in rvecs[i:i + 15]) + "\n")
        for r in rows:
            f.write("{} {} {} {} {} {!r} {!r}\n".format(*r))


def _write_model(d, norb, transfer_hr=None, extern=False, inter=True):
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n{}\n".format(norb))
        for a in range(norb):
            f.write("0.0 0.0 0.0\n" if a == 0 else "0.5 0.5 0.0\n")
    _write_w90(os.path.join(d, "transfer.dat"), "transfer",
               _bare_hr(norb) if transfer_hr is None else transfer_hr)
    u = np.zeros((NX, NY, norb, norb), dtype=complex)
    for a in range(norb):
        u[0, 0, a, a] = 1.4
    _write_w90(os.path.join(d, "coulombintra.dat"), "CoulombIntra", u)
    v = np.zeros((NX, NY, norb, norb), dtype=complex)
    for a in range(norb):
        for (rx, ry) in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            v[rx % NX, ry % NY, a, a] = 0.25
    if inter:
        _write_w90(os.path.join(d, "coulombinter.dat"), "CoulombInter", v)
    if extern:
        e = np.zeros((NX, NY, norb, norb), dtype=complex)
        for a in range(norb):
            e[0, 0, a, a] = 1.0
        _write_w90(os.path.join(d, "extern.dat"), "Extern", e)


def _input(d, out, *, chi0q_tensor="auto", chi0q_mode="calc", extern=False,
           extra_file=None, extra_param=None, eliashberg=None,
           coulomb_inter=True):
    inter = {"path_to_input": d, "Geometry": "geom.dat",
             "Transfer": "transfer.dat", "CoulombIntra": "coulombintra.dat"}
    if coulomb_inter:
        inter["CoulombInter"] = "coulombinter.dat"
    if extern:
        inter["Extern"] = "extern.dat"
    param = {"T": T_RUN, "CellShape": [NX, NY, 1], "SubShape": [1, 1, 1],
             "Nmat": NMAT, "filling": FILLING, "coeff_tail": 1.0}
    param.update(extra_param or {})
    fin = {"path_to_input": d, "interaction": inter}
    fin.update(extra_file or {})
    eli = {"solver_mode": "eigenvalue", "chi0q_mode": chi0q_mode,
           "chi0q_tensor": chi0q_tensor, "pairing_type": "singlet",
           "init_gap": "d_x2y2", "num_eigenvalues": 4,
           "eigenvalue_method": "arnoldi"}
    eli.update(eliashberg or {})
    return {"log": {"print_level": 0},
            "mode": {"mode": "RPA", "param": param},
            "file": {"input": fin,
                     "output": {"path_to_output": out}},
            "eliashberg": eli}


class _Capture:
    """Record lambda, the mu handed to chi0q, the kernel Green function and
    the prepared band, by wrapping the module functions for one run."""

    def __init__(self):
        self.data = {}

    def __enter__(self):
        self._saved = {}

        def wrap(name, hook):
            orig = getattr(sc, name)
            self._saved[name] = orig

            def f(*a, **k):
                return hook(orig, *a, **k)
            setattr(sc, name, f)

        def solve_eig(orig, *a, **k):
            r = orig(*a, **k)
            self.data["eig"] = np.array(r[0])
            return r

        def g2(orig, green_kw, *a, **k):
            self.data["green"] = np.array(green_kw)
            return orig(green_kw, *a, **k)

        def chi0(orig, input_dict, chi0q_tensor="auto", precomputed_mu=None,
                 prepared=None, static_only=False):
            self.data["mu_chi0"] = precomputed_mu
            self.data["prepared"] = prepared
            self.data["static_only"] = static_only
            r = orig(input_dict, chi0q_tensor=chi0q_tensor,
                     precomputed_mu=precomputed_mu, prepared=prepared,
                     static_only=static_only)
            self.data["chi0q"] = np.array(r)
            return r

        def bond_green(orig, *a, **k):
            r = orig(*a, **k)
            self.data["bond_green"] = np.array(r.full_sc)
            return r

        wrap("_solve_eigenvalue", solve_eig)
        wrap("_build_bond_green", bond_green)
        wrap("_calc_g2", g2)
        wrap("_calc_chi0q_internal", chi0)
        return self

    def __exit__(self, *exc):
        for name, orig in self._saved.items():
            setattr(sc, name, orig)
        return False


def _run(inp):
    with _Capture() as cap:
        sc.calc_eliashberg(inp)
    return cap.data


def _leading(data):
    ev = data["eig"]
    return ev[np.argmax(ev.real)].real


def _band_k(h):
    """eps_{ab}(k) = sum_R h_{ab}(R) e^{+ikR}, independent of both code
    paths (direct sum, sc.py k mesh)."""
    norb = h.shape[-1]
    kx = 2 * np.pi * np.arange(NX) / NX
    ky = 2 * np.pi * np.arange(NY) / NY
    eps = np.zeros((NX, NY, norb, norb), dtype=complex)
    for ix in range(NX):
        for iy in range(NY):
            ph = np.exp(1j * (kx[:, None] * _fold(ix, NX)
                              + ky[None, :] * _fold(iy, NY)))
            eps += ph[:, :, None, None] * h[ix, iy][None, None, :, :]
    return eps


# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------

class _ForbidDownstream:
    """Make every chi0q / Green / pair-bubble / kernel builder raise, so a
    refusal test proves the refusal precedes all of them."""

    def __enter__(self):
        import hwave.solver.rpa as rpa
        self._targets = [
            (sc, "_calc_chi0q_internal"), (sc, "_load_chi0q"),
            (sc, "_chi0q_from_prepared_solver"), (sc, "_calc_g2"),
            (sc, "_build_bond_green"), (sc, "_make_kernel_operator"),
            (rpa.RPA, "_calc_green"), (rpa.RPA, "_calc_chi0q")]
        self._saved = [(obj, name, getattr(obj, name))
                       for obj, name in self._targets]

        def forbidden(*a, **k):
            raise AssertionError("chi0q / G / kernel reached before the "
                                 "refusal")
        for obj, name in self._targets:
            setattr(obj, name, forbidden)
        return self

    def __exit__(self, *exc):
        for obj, name, orig in self._saved:
            setattr(obj, name, orig)
        return False


class TestTransModOracle(unittest.TestCase):
    """A spin-free trans_mod run equals a run whose Transfer IS that band."""

    def _oracle(self, norb, chi0q_tensor):
        h = _nematic_hr(norb)
        with tempfile.TemporaryDirectory() as d, \
                tempfile.TemporaryDirectory() as d_or:
            _write_model(d, norb)
            np.savez(os.path.join(d, "tm.npz"),
                     trans_mod=_trans_mod_from_hr(h))
            res = _run(_input(d, os.path.join(d, "out"),
                              chi0q_tensor=chi0q_tensor,
                              extra_file={"trans_mod": "tm.npz"}))
            _write_model(d_or, norb, transfer_hr=h)
            ref = _run(_input(d_or, os.path.join(d_or, "out"),
                              chi0q_tensor=chi0q_tensor))
            bare = _run(_input(d, os.path.join(d, "out_b"),
                               chi0q_tensor=chi0q_tensor))

        lam, lam_ref, lam_bare = _leading(res), _leading(ref), _leading(bare)
        self.assertAlmostEqual(lam, lam_ref, delta=1e-10)
        # the field must matter, or the oracle proves nothing
        self.assertGreater(abs(lam - lam_bare), 1e-2)

        band = res["prepared"]
        self.assertIsNotNone(band)
        # the chi0q solver runs the requested scheme
        self.assertEqual(band.solver.calc_scheme, chi0q_tensor)
        # the static solve computes only the zero bosonic frequency on the
        # band path too (issue #272)
        self.assertTrue(res["static_only"])
        self.assertEqual(res["chi0q"].shape[0], 1)

        # mu handed to chi0q == the root of sum_k,n f(e_n(k) - mu) =
        # filling * norb * Nk, by plain bisection on an independent band
        ev = np.linalg.eigvalsh(_band_k(h)).ravel()
        target = FILLING * norb * NX * NY

        def count(mu):
            x = np.clip((ev - mu) / T_RUN, -700.0, 700.0)
            return np.sum(1.0 / (np.exp(x) + 1.0)) - target
        lo, hi = ev.min() - 10.0, ev.max() + 10.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if count(mid) < 0.0:
                lo = mid
            else:
                hi = mid
        self.assertLessEqual(abs(res["mu_chi0"] - 0.5 * (lo + hi)), 1e-10)

        # kernel G == 1/(i w_n + mu - eps_tm(k)), eps_tm by direct sum
        eps = _band_k(h)                               # (Nx, Ny, a, b)
        g = res["green"]                               # (a, b, Nx, Ny, 1, nmat)
        nmat = g.shape[-1]
        mu = res["mu_chi0"]
        for n in (nmat // 2, nmat // 2 + 3, 0):
            wn = (2 * (n - nmat // 2) + 1) * np.pi * T_RUN
            ginv = (1j * wn + mu) * np.eye(norb)[None, None] - eps
            gref = np.linalg.inv(ginv)                 # (Nx, Ny, a, b)
            gsc = np.transpose(g[:, :, :, :, 0, n], (2, 3, 0, 1))
            self.assertLessEqual(np.abs(gsc - gref).max(), 1e-12,
                                 "frequency index {}".format(n))
        return lam, lam_ref, lam_bare

    def test_one_orbital_reduced(self):
        self._oracle(1, "reduced")

    def test_one_orbital_general(self):
        self._oracle(1, "general")

    def test_two_orbital_reduced(self):
        self._oracle(2, "reduced")

    def test_two_orbital_general(self):
        self._oracle(2, "general")

    def test_auto_scheme_follows_the_frozen_rule(self):
        # chi0q_tensor = "auto" with CoulombInter declared: the prepared
        # solver must run the scheme hwave_sc's frozen auto rule selects
        with tempfile.TemporaryDirectory() as d:
            _write_model(d, 1)
            np.savez(os.path.join(d, "tm.npz"),
                     trans_mod=_trans_mod_from_hr(_nematic_hr(1)))
            inp = _input(d, os.path.join(d, "o"), chi0q_tensor="auto",
                         extra_file={"trans_mod": "tm.npz"})
            expected = sc._auto_chi0q_tensor_scheme(
                inp["file"]["input"]["interaction"])
            res = _run(inp)
        self.assertEqual(expected, "general")
        self.assertEqual(res["prepared"].solver.calc_scheme, expected)


class TestGreenInit(unittest.TestCase):
    """A UHFk solution supplied as green_init defines the same band as the
    trans_mod UHFk saved with it.

    The RPA reader turns green_init into UHFk's T + Hartree + Fock of G
    (issue #278; the per-type comparison is TestGreenInitIsUHFkMeanField
    and tests/test_rpa_green_init_fock.py). CoulombIntra here. Two orbitals,
    so the orbital-dependent Hartree shift changes the band beyond a
    constant that mu would absorb.
    """

    def test_green_init_matches_saved_trans_mod(self):
        import hwave.qlms as qlms
        norb = 2
        with tempfile.TemporaryDirectory() as d:
            _write_model(d, norb, inter=False)
            out = os.path.join(d, "uhfk")
            qlms.run(input_dict={
                "log": {"print_level": 0},
                "mode": {"mode": "UHFk",
                         "param": {"Ncond": int(round(FILLING * 2 * norb
                                                      * NX * NY)),
                                   "IterationMax": 3000, "EPS": 12,
                                   "Mix": 0.5, "RndSeed": 1, "T": T_RUN,
                                   "CellShape": [NX, NY, 1],
                                   "SubShape": [1, 1, 1]}},
                "file": {"input": {"path_to_input": "",
                                   "interaction": {
                                       "path_to_input": d,
                                       "Geometry": "geom.dat",
                                       "Transfer": "transfer.dat",
                                       "CoulombIntra": "coulombintra.dat"}},
                         "output": {"path_to_output": out, "green": "green",
                                    "rpa": "trans_mod"}}})
            tm = np.load(os.path.join(out, "trans_mod.npz"))["trans_mod"]
            # paramagnetic solution (else both inputs would be refused)
            self.assertLess(np.abs(tm[:, :norb, :norb]
                                   - tm[:, norb:, norb:]).max(), 1e-10)

            def lam(extra):
                return _leading(_run(_input(
                    d, os.path.join(d, "o"), coulomb_inter=False,
                    extra_file=extra)))
            lam_tm = lam({"trans_mod": os.path.join(out, "trans_mod.npz")})
            lam_gi = lam({"green_init": os.path.join(out, "green.npz")})
            lam_bare = lam(None)
        self.assertAlmostEqual(lam_gi, lam_tm, delta=1e-8)
        self.assertGreater(abs(lam_tm - lam_bare), 1e-3)


class TestSpinDependentRefused(unittest.TestCase):
    """A spin-dependent band cannot enter the spinless pairing kernel: it is
    refused, naming issue #88, before any chi0q is computed or loaded."""

    def _assert_refused(self, chi0q_mode, extra_file=None, extern=False,
                        extra_param=None, tm=None):
        with tempfile.TemporaryDirectory() as d:
            _write_model(d, 1, extern=extern)
            if tm is not None:
                np.savez(os.path.join(d, "tm.npz"), trans_mod=tm)
            inp = _input(d, os.path.join(d, "out"), chi0q_mode=chi0q_mode,
                         extern=extern, extra_file=extra_file,
                         extra_param=extra_param)
            with _ForbidDownstream(), \
                    self.assertRaises(ValueError) as cm:
                sc.calc_eliashberg(inp)
        self.assertIn("#88", str(cm.exception))
        return str(cm.exception)

    def _ferro(self):
        h = _bare_hr(1)
        m = np.zeros_like(h)
        m[0, 0, 0, 0] = 0.3
        return _trans_mod_from_hr(h + m, h - m)

    def _mixing(self):
        return _trans_mod_from_hr(_bare_hr(1), mixing=np.array([[0.2]]))

    def test_ferro_trans_mod(self):
        for mode in ("calc", "load"):
            with self.subTest(chi0q_mode=mode):
                msg = self._assert_refused(
                    mode, extra_file={"trans_mod": "tm.npz"}, tm=self._ferro())
                self.assertIn("spin-diag", msg)

    def test_spin_mixing_trans_mod(self):
        for mode in ("calc", "load"):
            with self.subTest(chi0q_mode=mode):
                msg = self._assert_refused(
                    mode, extra_file={"trans_mod": "tm.npz"},
                    tm=self._mixing())
                self.assertIn("spin-mixing", msg)

    def test_active_extern(self):
        for mode in ("calc", "load"):
            with self.subTest(chi0q_mode=mode):
                msg = self._assert_refused(
                    mode, extern=True, extra_param={"coeff_extern": 0.2})
                self.assertIn("Extern", msg)

    def test_mixed_case_key_is_honoured(self):
        # 'Trans_Mod' must not be skipped: a ferro field is refused
        self._assert_refused("calc", extra_file={"Trans_Mod": "tm.npz"},
                             tm=self._ferro())

    def test_bond_path_refuses_too(self):
        with tempfile.TemporaryDirectory() as d:
            _write_model(d, 1)
            np.savez(os.path.join(d, "tm.npz"), trans_mod=self._ferro())
            inp = _input(d, os.path.join(d, "out"),
                         extra_file={"trans_mod": "tm.npz"},
                         eliashberg={"bond_channels": True})
            with _ForbidDownstream(), self.assertRaises(ValueError) as cm:
                sc.calc_eliashberg(inp)
        self.assertIn("#88", str(cm.exception))


class TestMixedCaseAndControl(unittest.TestCase):

    def test_mixed_case_spelling_gives_same_lambda(self):
        with tempfile.TemporaryDirectory() as d:
            _write_model(d, 1)
            np.savez(os.path.join(d, "tm.npz"),
                     trans_mod=_trans_mod_from_hr(_nematic_hr(1)))
            a = _leading(_run(_input(d, os.path.join(d, "o1"),
                                     extra_file={"trans_mod": "tm.npz"})))
            b = _leading(_run(_input(d, os.path.join(d, "o2"),
                                     extra_file={"Trans_Mod": "tm.npz"})))
        self.assertEqual(a, b)

    def test_duplicate_spellings_refused(self):
        with tempfile.TemporaryDirectory() as d:
            _write_model(d, 1)
            np.savez(os.path.join(d, "tm.npz"),
                     trans_mod=_trans_mod_from_hr(_bare_hr(1)))
            inp = _input(d, os.path.join(d, "o"),
                         extra_file={"trans_mod": "tm.npz",
                                     "Trans_Mod": "tm.npz"})
            with self.assertRaises(ValueError) as cm:
                sc.calc_eliashberg(inp)
        self.assertIn("more than once", str(cm.exception))

    def test_bare_trans_mod_is_bit_identical(self):
        for tensor in ("reduced", "general"):
            with self.subTest(chi0q_tensor=tensor), \
                    tempfile.TemporaryDirectory() as d:
                _write_model(d, 1)
                np.savez(os.path.join(d, "tm.npz"),
                         trans_mod=_trans_mod_from_hr(_bare_hr(1)))
                with_tm = _run(_input(d, os.path.join(d, "o1"),
                                      chi0q_tensor=tensor,
                                      extra_file={"trans_mod": "tm.npz"}))
                without = _run(_input(d, os.path.join(d, "o2"),
                                      chi0q_tensor=tensor))
                self.assertTrue(np.array_equal(with_tm["eig"],
                                               without["eig"]))
                self.assertTrue(np.array_equal(with_tm["green"],
                                               without["green"]))

    def test_extern_is_a_source_whatever_the_coefficient(self):
        for coeff in (0.0, 0.5):
            self.assertEqual(sc._one_body_band_sources(
                {"file": {"input": {"interaction": {"extern": "e.dat"}}},
                 "mode": {"param": {"Coeff_Extern": coeff}}}), ("Extern",))

    def test_zero_coeff_extern_equals_run_without_extern(self):
        # used to fail with "Unexpected chi0q shape (2, ...)": RPA builds a
        # two-block H0 for any declared Extern
        for tensor in ("reduced", "general"):
            with self.subTest(chi0q_tensor=tensor), \
                    tempfile.TemporaryDirectory() as d:
                _write_model(d, 1, extern=True)
                with_ext = _run(_input(d, os.path.join(d, "o1"),
                                       chi0q_tensor=tensor, extern=True,
                                       extra_param={"coeff_extern": 0.0}))
                without = _run(_input(d, os.path.join(d, "o2"),
                                      chi0q_tensor=tensor))
                for key in ("eig", "green", "chi0q"):
                    self.assertTrue(np.array_equal(with_ext[key],
                                                   without[key]), key)


class TestLoadAndBond(unittest.TestCase):

    def test_load_mode_uses_the_band(self):
        """chi0q written by an RPA run on the trans_mod band, loaded by
        hwave_sc with the same trans_mod: same lambda as the calc run."""
        import hwave.qlms as qlms
        h = _nematic_hr(1)
        with tempfile.TemporaryDirectory() as d:
            _write_model(d, 1)
            np.savez(os.path.join(d, "tm.npz"),
                     trans_mod=_trans_mod_from_hr(h))
            out = os.path.join(d, "rpa")
            rpa_in = _input(d, out, extra_file={"trans_mod": "tm.npz"})
            rpa_in["mode"]["calc_scheme"] = "general"
            rpa_in["file"]["output"]["chi0q"] = "chi0q"
            del rpa_in["eliashberg"]
            qlms.run(input_dict=rpa_in)
            calc = _leading(_run(_input(d, os.path.join(d, "o_calc"),
                                        chi0q_tensor="general",
                                        extra_file={"trans_mod": "tm.npz"})))
            load = _leading(_run(_input(d, out, chi0q_mode="load",
                                        extra_file={"trans_mod": "tm.npz"})))
            load_bare = _leading(_run(_input(d, out, chi0q_mode="load")))
        self.assertAlmostEqual(load, calc, delta=1e-9)
        self.assertGreater(abs(load - load_bare), 1e-3)

    def test_bond_path_uses_the_band(self):
        h = _nematic_hr(1)
        bond = {"bond_channels": True, "bond_coeff_tail": 1.0}
        with tempfile.TemporaryDirectory() as d, \
                tempfile.TemporaryDirectory() as d_or:
            _write_model(d, 1)
            np.savez(os.path.join(d, "tm.npz"),
                     trans_mod=_trans_mod_from_hr(h))
            lam = _leading(_run(_input(d, os.path.join(d, "o"),
                                       extra_file={"trans_mod": "tm.npz"},
                                       eliashberg=bond)))
            _write_model(d_or, 1, transfer_hr=h)
            lam_ref = _leading(_run(_input(d_or, os.path.join(d_or, "o"),
                                           eliashberg=bond)))
            lam_bare = _leading(_run(_input(d, os.path.join(d, "ob"),
                                            eliashberg=bond)))
        self.assertAlmostEqual(lam, lam_ref, delta=1e-10)
        self.assertGreater(abs(lam - lam_bare), 1e-3)


class TestFlexPathsUnchanged(unittest.TestCase):
    """chi and G of the FLEX path come from the FLEX run: trans_mod in the
    hwave_sc input changes nothing and is reported."""

    def test_flex_lambda_unchanged_and_warned(self):
        import hwave.qlms as qlms
        with tempfile.TemporaryDirectory() as d:
            _write_model(d, 1)
            np.savez(os.path.join(d, "tm.npz"),
                     trans_mod=_trans_mod_from_hr(_nematic_hr(1)))
            out = os.path.join(d, "flex")
            qlms.run(input_dict={
                "log": {"print_level": 0},
                "mode": {"mode": "FLEX", "calc_scheme": "reduced",
                         "param": {"T": 0.25, "filling": FILLING,
                                   "CellShape": [4, 4, 1],
                                   "SubShape": [1, 1, 1], "Nmat": 32,
                                   "IterationMax": 200, "Mix": 0.3,
                                   "EPS": 8}},
                "file": {"input": {"path_to_input": "",
                                   "interaction": {
                                       "path_to_input": d,
                                       "Geometry": "geom.dat",
                                       "Transfer": "transfer.dat",
                                       "CoulombIntra": "coulombintra.dat"}},
                         "output": {"path_to_output": out, "chiq_s": "chiq_s",
                                    "chiq_c": "chiq_c", "green": "green"}}})

            def eli(extra):
                inp = {"log": {"print_level": 0},
                       "mode": {"mode": "FLEX",
                                "param": {"T": 0.25, "filling": FILLING,
                                          "CellShape": [4, 4, 1],
                                          "SubShape": [1, 1, 1],
                                          "Nmat": 32}},
                       "file": {"input": {"path_to_input": d,
                                          "path_to_flex_output": out,
                                          "interaction": {
                                              "path_to_input": d,
                                              "Geometry": "geom.dat",
                                              "Transfer": "transfer.dat",
                                              "CoulombIntra":
                                                  "coulombintra.dat"}},
                                "output": {"path_to_output":
                                           os.path.join(d, "eli")}},
                       "eliashberg": {"solver_mode": "eigenvalue",
                                      "chi0q_mode": "flex",
                                      "pairing_type": "singlet",
                                      "num_eigenvalues": 4}}
                inp["file"]["input"].update(extra)
                return inp

            plain = _run(eli({}))
            with self.assertLogs("hwave_sc", level="WARNING") as logs:
                with_tm = _run(eli({"trans_mod": "tm.npz"}))
        self.assertTrue(np.array_equal(plain["eig"], with_tm["eig"]))
        self.assertTrue(any("trans_mod is ignored" in m for m in logs.output),
                        logs.output)

    def test_dynamic_warns(self):
        from hwave.solver import eliashberg_dynamic
        saved_solve = eliashberg_dynamic.solve_dynamic
        saved_prereq = sc._validate_dynamic_prereqs
        eliashberg_dynamic.solve_dynamic = lambda inp: "solved"
        sc._validate_dynamic_prereqs = lambda inp: None
        try:
            inp = {"mode": {"param": {"T": 0.1, "CellShape": [4, 4, 1],
                                      "SubShape": [1, 1, 1],
                                      "filling": 0.5}},
                   "file": {"input": {"green_init": "g.npz"},
                            "output": {"path_to_output": "unused"}},
                   "eliashberg": {"frequency": "dynamic",
                                  "chi0q_mode": "flex"}}
            with self.assertLogs("hwave_sc", level="WARNING") as logs:
                self.assertEqual(sc.calc_eliashberg(inp), "solved")
        finally:
            eliashberg_dynamic.solve_dynamic = saved_solve
            sc._validate_dynamic_prereqs = saved_prereq
        self.assertTrue(any("green_init is ignored" in m
                            for m in logs.output), logs.output)



def _ferro_green(norb=1, m=0.05):
    """A UHFk-layout green (cellvol, 2, norb, 2, norb) with unequal spins."""
    g = np.zeros((NX * NY, 2, norb, 2, norb), dtype=complex)
    for a in range(norb):
        g[0, 0, a, 0, a] = FILLING + m
        g[0, 1, a, 1, a] = FILLING - m
    return g


class TestRefusalOrdering(unittest.TestCase):
    """The refusal is made on the ASSEMBLED band, before any chi0q, Green
    function, pair bubble or kernel is built."""

    def _refused(self, extra_file, files, extern=False, extra_param=None,
                 inter_key="Extern", norb=1):
        with tempfile.TemporaryDirectory() as d:
            _write_model(d, norb, extern=extern)
            for name, payload in files.items():
                np.savez(os.path.join(d, name), **payload)
            inp = _input(d, os.path.join(d, "out"), extra_file=extra_file,
                         extra_param=extra_param)
            if extern:
                inp["file"]["input"]["interaction"][inter_key] = "extern.dat"
            with _ForbidDownstream(), self.assertRaises(ValueError) as cm:
                sc.calc_eliashberg(inp)
        self.assertIn("#88", str(cm.exception))
        return str(cm.exception)

    def test_trans_mod(self):
        h = _bare_hr(1)
        m = np.zeros_like(h)
        m[0, 0, 0, 0] = 0.3
        self._refused(
            {"trans_mod": "tm.npz"},
            {"tm.npz": {"trans_mod": _trans_mod_from_hr(h + m, h - m)}})

    def test_green_init(self):
        self._refused({"green_init": "g.npz"},
                      {"g.npz": {"green": _ferro_green()}})

    def test_green_init_mixed_case(self):
        self._refused({"Green_Init": "g.npz"},
                      {"g.npz": {"green": _ferro_green()}})

    def test_extern_mixed_case(self):
        # 'extern' table key and 'Coeff_Extern' are both honoured
        msg = self._refused({}, {}, extern=True, inter_key="extern",
                            extra_param={"Coeff_Extern": 0.2})
        self.assertIn("Extern", msg)


class TestAssembledBandCheck(unittest.TestCase):
    """The criterion is the assembled band, not the raw inputs."""

    def _band(self, norb=1, extra_file=None, files=None, extern_rows=None,
              coeff=None, inter=None):
        """_mean_field_band on a small model; returns None or raises."""
        with tempfile.TemporaryDirectory() as d:
            _write_model(d, norb)
            for name, payload in (files or {}).items():
                np.savez(os.path.join(d, name), **payload)
            inp = _input(d, os.path.join(d, "out"), extra_file=extra_file,
                         extra_param=({"coeff_extern": coeff}
                                      if coeff is not None else None))
            if extern_rows is not None:
                with open(os.path.join(d, "ext.dat"), "w") as f:
                    f.write("Extern\n{}\n{}\n".format(norb, len(extern_rows)))
                    f.write(" ".join("1" for _ in extern_rows) + "\n")
                    for r in extern_rows:
                        f.write("{} {} {} {} {} {!r} 0.0\n".format(*r))
                inp["file"]["input"]["interaction"]["Extern"] = "ext.dat"
            if inter is not None:
                files_ = inp["file"]["input"]["interaction"]
                for k in ("CoulombIntra", "CoulombInter"):
                    files_.pop(k, None)
                tname, fname = inter
                files_[tname] = fname
            sources = sc._one_body_band_sources(inp)
            return sc._mean_field_band(inp, sources, "general",
                                       (NX, NY, 1), norb)

    def test_green_init_band_is_the_uhf_band(self):
        """Spin-free G with nonzero cells R = (0, +-1), a CoulombInter bond
        of this model. RPA's green_init assembly is UHFk's T + Hartree +
        Fock of G (issue #278): the bond's Fock term -V G(R) is in the band,
        and no warning is emitted."""
        g = np.zeros((NX * NY, 2, 1, 2, 1), dtype=complex)
        g[0, 0, 0, 0, 0] = g[0, 1, 0, 1, 0] = FILLING
        g[1, 0, 0, 0, 0] = g[1, 1, 0, 1, 0] = 0.2              # R = (0, 1)
        g[NY - 1, 0, 0, 0, 0] = g[NY - 1, 1, 0, 1, 0] = 0.2     # R = (0, -1)
        records = []
        handler = logging.Handler()
        handler.emit = records.append
        lg = logging.getLogger("hwave_sc")
        lg.addHandler(handler)
        try:
            band = self._band(extra_file={"green_init": "g.npz"},
                              files={"g.npz": {"green": g}})
        finally:
            lg.removeHandler(handler)
        self.assertFalse(any("Fock" in r.getMessage() for r in records))
        # T + Hartree + Fock: U n + z V (2n), z = 4 bonds, n per spin; the
        # Fock term of the R = (0, +-1) bonds is -V G(R) on those bonds
        eps = np.transpose(band.epsilon_k[:, :, :, :, 0], (2, 3, 0, 1))
        hartree = 1.4 * FILLING + 4 * 0.25 * 2 * FILLING
        fock = np.zeros((NX, NY, 1, 1), dtype=complex)
        fock[0, 1, 0, 0] = fock[0, NY - 1, 0, 0] = -0.25 * 0.2
        self.assertLess(np.abs(eps - (_band_k(_bare_hr(1)) + hartree
                                      + _band_k(fock))).max(), 1e-12)

    def test_green_init_spin_dependent_in_another_cell_is_refused(self):
        # G-level check: a moment in cell 1 is refused even though RPA's
        # Hartree band (g0[0] only) would be spin-free
        g = np.zeros((NX * NY, 2, 1, 2, 1), dtype=complex)
        g[0, 0, 0, 0, 0] = g[0, 1, 0, 1, 0] = FILLING
        g[1, 0, 0, 0, 0] = 0.2
        g[1, 1, 0, 1, 0] = -0.2
        with self.assertRaises(ValueError) as cm:
            self._band(extra_file={"green_init": "g.npz"},
                       files={"g.npz": {"green": g}})
        self.assertIn("#88", str(cm.exception))
        self.assertIn("green_init describes a spin-dependent state",
                      str(cm.exception))

    def test_extern_boundary(self):
        # on-site Extern 1.0: H_up - H_down = 2 * coeff_extern
        rows = [(0, 0, 0, 1, 1, 1.0)]
        band = self._band(extern_rows=rows, coeff=0.4e-10)
        # accepted, and the band used is the checked up-up block: Extern
        # included (bare + coeff), not dropped
        eps = np.transpose(band.epsilon_k[:, :, :, :, 0], (2, 3, 0, 1))
        bare = _band_k(_bare_hr(1))
        self.assertLess(np.abs(eps - (bare + 0.4e-10)).max(), 1e-13)
        self.assertFalse(band.same_as_bare)
        with self.assertRaises(ValueError) as cm:
            self._band(extern_rows=rows, coeff=0.6e-10)
        self.assertIn("#88", str(cm.exception))

    def test_extern_cancelling_on_the_grid(self):
        # R = 0 and R = -NX wrap onto the same slot and cancel exactly
        rows = [(0, 0, 0, 1, 1, 1.0), (-NX, 0, 0, 1, 1, -1.0)]
        band = self._band(extern_rows=rows, coeff=0.2)
        self.assertTrue(band.same_as_bare)

    def test_each_interaction_type(self):
        """Per interaction type (two orbitals, on-site terms): a spin-free G
        gives an accepted band, a ferro G is refused -- by the G-level
        check, for every type (including those whose Hartree term would
        carry no spin dependence)."""
        norb = 2
        m = 0.05

        def g_of(diag, off):
            g = np.zeros((NX * NY, 2, norb, 2, norb), dtype=complex)
            for s_, sign in ((0, 1.0), (1, -1.0)):
                g[0, s_, 0, s_, 0] = g[0, s_, 1, s_, 1] = 0.4 + sign * diag
                g[0, s_, 0, s_, 1] = g[0, s_, 1, s_, 0] = 0.1 + sign * off
            return g
        g_free = g_of(0.0, 0.0)
        g_ferro = g_of(m, m)
        for tname in ("CoulombIntra", "Hund", "Ising", "PairHop",
                      "CoulombInter", "Exchange", "PairLift"):
            with self.subTest(type=tname), tempfile.TemporaryDirectory() as d:
                if tname == "CoulombIntra":
                    rows = [(0, 0, 0, 1, 1, 1.0), (0, 0, 0, 2, 2, 1.0)]
                else:
                    rows = [(0, 0, 0, 1, 2, 0.5), (0, 0, 0, 2, 1, 0.5)]
                path = os.path.join(d, "x.dat")
                with open(path, "w") as f:
                    f.write("{}\n{}\n1\n1\n".format(tname, norb))
                    for r in rows:
                        f.write("{} {} {} {} {} {!r} 0.0\n".format(*r))
                for label, g, ok in (("free", g_free, True),
                                     ("ferro", g_ferro, False)):
                    np.savez(os.path.join(d, "g.npz"), green=g)

                    def call():
                        return self._band(
                            norb=norb,
                            extra_file={"green_init":
                                        os.path.join(d, "g.npz")},
                            inter=(tname, path))
                    if ok:
                        self.assertIsNotNone(call(), label)
                    else:
                        with self.assertRaises(ValueError) as cm:
                            call()
                        # |G_up - G_down| = 2 m = 0.1
                        self.assertIn("#88", str(cm.exception))
                        self.assertIn("max |G_up - G_down| = 1.000e-01",
                                      str(cm.exception))


class TestSpinFreeTolerance(unittest.TestCase):
    """The spin-independence test is ABSOLUTE (SPIN_FREE_ATOL)."""

    def _lam_or_error(self, tm):
        with tempfile.TemporaryDirectory() as d:
            _write_model(d, 1)
            np.savez(os.path.join(d, "tm.npz"), trans_mod=tm)
            try:
                return _leading(_run(_input(
                    d, os.path.join(d, "out"),
                    extra_file={"trans_mod": "tm.npz"}))), None
            except ValueError as exc:
                return None, str(exc)

    def test_large_offset_does_not_hide_a_splitting(self):
        h = _nematic_hr(1)
        shift = np.zeros_like(h)
        shift[0, 0, 0, 0] = 1.0e3
        split = np.zeros_like(h)
        split[0, 0, 0, 0] = 0.5e-6
        tm = _trans_mod_from_hr(h + shift + split, h + shift - split)
        # a relative comparison would call these blocks equal
        self.assertTrue(np.allclose(tm[:, 0, 0], tm[:, 1, 1]))
        lam, err = self._lam_or_error(tm)
        self.assertIsNone(lam)
        self.assertIn("#88", err)

    def test_spin_mixing_threshold(self):
        h = _nematic_hr(1)
        self.assertEqual(sc.SPIN_FREE_ATOL, 1.0e-10)
        lam_ref, _ = self._lam_or_error(_trans_mod_from_hr(h))
        above = _trans_mod_from_hr(h, mixing=np.array([[1.0e-9]]))
        lam, err = self._lam_or_error(above)
        self.assertIsNone(lam)
        self.assertIn("#88", err)
        below = _trans_mod_from_hr(h, mixing=np.array([[1.0e-12]]))
        lam, err = self._lam_or_error(below)
        self.assertIsNone(err)
        self.assertAlmostEqual(lam, lam_ref, delta=1e-9)


class TestMixedCaseGreenInit(unittest.TestCase):

    def test_green_init_spelling_gives_same_lambda(self):
        # two orbitals with unequal occupations: the Hartree shift is
        # orbital-dependent, so the green_init band differs from the bare one
        g = np.zeros((NX * NY, 2, 2, 2, 2), dtype=complex)
        for s_ in range(2):
            g[0, s_, 0, s_, 0] = 0.6
            g[0, s_, 1, s_, 1] = 0.3
        with tempfile.TemporaryDirectory() as d:
            _write_model(d, 2)
            np.savez(os.path.join(d, "g.npz"), green=g)
            a = _leading(_run(_input(d, os.path.join(d, "o1"),
                                     extra_file={"green_init": "g.npz"})))
            b = _leading(_run(_input(d, os.path.join(d, "o2"),
                                     extra_file={"Green_Init": "g.npz"})))
            bare = _leading(_run(_input(d, os.path.join(d, "o3"))))
        self.assertEqual(a, b)
        self.assertGreater(abs(a - bare), 1e-3)


class TestExternalBondGreen(unittest.TestCase):
    """With an external bond_green the band plays no role: trans_mod /
    green_init / Extern are not read, not validated, and reported."""

    def test_external_bond_green_ignores_band_inputs(self):
        bond = {"bond_channels": True}
        with tempfile.TemporaryDirectory() as d:
            _write_model(d, 1)
            # a valid external Green function: the bare model's own
            ref = _run(_input(d, os.path.join(d, "o0"), eliashberg=bond))
            np.savez(os.path.join(d, "green.npz"),
                     green=sc._green_sc_to_canonical(ref["bond_green"]))
            bond_ext = dict(bond, bond_green=os.path.join(d, "green.npz"))
            h = _bare_hr(1)
            m = np.zeros_like(h)
            m[0, 0, 0, 0] = 0.3
            np.savez(os.path.join(d, "ferro.npz"),
                     trans_mod=_trans_mod_from_hr(h + m, h - m))
            plain = _run(_input(d, os.path.join(d, "o1"),
                                eliashberg=bond_ext))
            cases = [
                ("trans_mod", {"trans_mod": "ferro.npz"}, None, None),
                ("trans_mod", {"trans_mod": "does_not_exist.npz"}, None,
                 None),
                ("green_init", {"Green_Init": "does_not_exist.npz"}, None,
                 None),
                ("Extern", {}, "does_not_exist.dat", {"Coeff_Extern": 0.2}),
            ]
            for name, extra_file, extern_file, param in cases:
                with self.subTest(name=name, extra=extra_file):
                    inp = _input(d, os.path.join(d, "o2"),
                                 eliashberg=bond_ext, extra_file=extra_file,
                                 extra_param=param)
                    if extern_file is not None:
                        inp["file"]["input"]["interaction"]["extern"] = \
                            extern_file
                    with self.assertLogs("hwave_sc", level="WARNING") as logs:
                        res = _run(inp)
                    self.assertTrue(np.array_equal(res["eig"], plain["eig"]))
                    self.assertTrue(any(
                        "{} is ignored".format(name) in msg
                        and "bond_green" in msg
                        for msg in logs.output), logs.output)


class TestCheckedBandIsUsedBand(unittest.TestCase):
    """The band that is used is exactly the band that was checked."""

    def test_check_pass_restores_mutated_hamiltonian_arrays(self):
        # a stand-in solver whose assembly modifies ham_trans_r in place,
        # as RPA's spin-orbital green_init branch does
        class _Ham:
            pass

        class _Solver:
            norb = 1

            def __init__(self):
                self.ham_info = _Ham()
                self.ham_info.ham_trans_r = np.ones((4, 2, 2), complex)
                self.ham_info.ham_trans_q = np.ones((4, 2, 2), complex)
                self.ham_info.ham_extern_r = None
                self.ham_info.ham_extern_q = None

            def _calc_epsilon_k(self, gi):
                self.ham_info.ham_trans_r[0] += 7.0     # in-place, as RPA
                self.H0_k = np.zeros((1, 4, 2, 2), complex)
                self._find_block_diagonal(None)

        solver = _Solver()
        before = solver.ham_info.ham_trans_r.copy()
        up = sc._checked_band(solver, {}, ("green_init",))
        self.assertTrue(np.array_equal(solver.ham_info.ham_trans_r, before))
        self.assertEqual(up.shape, (4, 1, 1))
        self.assertFalse(hasattr(solver, "_find_block_diagonal"))

    def test_green_init_band_equals_single_rpa_assembly(self):
        import hwave.qlmsio.read_input_k as rik
        import hwave.solver.rpa as rpa
        norb = 2
        g = np.zeros((NX * NY, 2, norb, 2, norb), dtype=complex)
        for s_ in range(2):
            g[0, s_, 0, s_, 0] = 0.6
            g[0, s_, 1, s_, 1] = 0.3
            g[0, s_, 0, s_, 1] = g[0, s_, 1, s_, 0] = 0.05
        with tempfile.TemporaryDirectory() as d:
            _write_model(d, norb)
            np.savez(os.path.join(d, "g.npz"), green=g)
            inp = _input(d, os.path.join(d, "o"),
                         extra_file={"green_init": "g.npz"})
            band = sc._mean_field_band(inp, ("green_init",), "general",
                                       (NX, NY, 1), norb)
            # a fresh RPA solver, one assembly
            info = {"path_to_input": d,
                    "interaction": inp["file"]["input"]["interaction"],
                    "green_init": "g.npz"}
            rio = rik.QLMSkInput(info)
            fresh = rpa.RPA(rio.get_param("ham"), {"print_level": 0},
                            {"mode": "RPA", "calc_scheme": "general",
                             "param": dict(inp["mode"]["param"])})
            gi = rio.get_param("green")
            gi.update(fresh.read_init(info))
            fresh._calc_epsilon_k(gi)
        self.assertTrue(np.array_equal(band.solver.ham_info.ham_trans_r,
                                       fresh.ham_info.ham_trans_r))
        eps = np.transpose(band.epsilon_k[:, :, :, :, 0], (2, 3, 0, 1))
        ref = np.asarray(fresh.H0_k)[0].reshape(NX, NY, norb, norb)
        self.assertLess(np.abs(eps - ref).max(), 1e-12)
        # and the chi0q solver holds exactly that band
        self.assertTrue(np.array_equal(
            np.asarray(band.solver.H0_k)[0],
            eps.reshape(NX * NY, norb, norb)))

    def test_spin_orbital_green_init_is_refused_upstream(self):
        with tempfile.TemporaryDirectory() as d:
            _write_model(d, 1)
            inp = _input(d, os.path.join(d, "o"),
                         extra_file={"green_init": "g.npz"})
            inp["mode"]["enable_spin_orbital"] = True
            with self.assertRaises(ValueError) as cm:
                sc.calc_eliashberg(inp)
        self.assertIn("enable_spin_orbital", str(cm.exception))

    def test_extern_cancelling_a_split_with_large_offset(self):
        """trans_mod = h + 1e3 with a +-0.3 spin split, cancelled by an
        on-site Extern with coeff_extern = -0.3: the COMBINED band is
        spin-free and accepted, and the band used is that combined band."""
        h = _nematic_hr(1)
        off = np.zeros_like(h)
        off[0, 0, 0, 0] = 1.0e3
        sp = np.zeros_like(h)
        sp[0, 0, 0, 0] = 0.3
        tm = _trans_mod_from_hr(h + off + sp, h + off - sp)
        with tempfile.TemporaryDirectory() as d:
            _write_model(d, 1, extern=True)
            np.savez(os.path.join(d, "tm.npz"), trans_mod=tm)
            inp = _input(d, os.path.join(d, "o"), extern=True,
                         extra_file={"trans_mod": "tm.npz"},
                         extra_param={"coeff_extern": -0.3})
            band = sc._mean_field_band(inp, sc._one_body_band_sources(inp),
                                       "general", (NX, NY, 1), 1)
        eps = np.transpose(band.epsilon_k[:, :, :, :, 0], (2, 3, 0, 1))
        self.assertLess(np.abs(eps - (_band_k(h) + 1.0e3)).max(), 1e-10)


class TestNonFiniteRefused(unittest.TestCase):

    def _refused(self, mode, coeff, nan_row=False, bond=False):
        with tempfile.TemporaryDirectory() as d:
            _write_model(d, 1, extern=True)
            if nan_row:
                with open(os.path.join(d, "extern.dat"), "w") as f:
                    f.write("Extern\n1\n1\n1\n0 0 0 1 1 nan 0.0\n")
            inp = _input(d, os.path.join(d, "o"), chi0q_mode=mode,
                         extern=True, extra_param={"coeff_extern": coeff},
                         eliashberg={"bond_channels": True} if bond else None)
            with _ForbidDownstream(), self.assertRaises(ValueError) as cm:
                sc.calc_eliashberg(inp)
        self.assertIn("non-finite", str(cm.exception))

    def test_nan_inf_coeff_and_nan_table(self):
        for mode, bond in (("calc", False), ("load", False), ("calc", True)):
            for coeff, nan_row in ((float("nan"), False),
                                   (float("inf"), False), (0.2, True)):
                with self.subTest(mode=mode, bond=bond, coeff=coeff,
                                  nan_row=nan_row):
                    self._refused(mode, coeff, nan_row=nan_row, bond=bond)


class TestUHFkLayout(unittest.TestCase):
    """The trans_mod layout the tests write is the one UHFk writes: at
    U = 0 UHFk's trans_mod is the Transfer table itself."""

    def test_uhfk_trans_mod_at_zero_u(self):
        import hwave.qlms as qlms
        norb = 2
        h = _nematic_hr(norb)
        with tempfile.TemporaryDirectory() as d:
            _write_model(d, norb, transfer_hr=h, inter=False)
            with open(os.path.join(d, "coulombintra.dat"), "w") as f:
                f.write("CoulombIntra\n2\n1\n1\n"
                        "0 0 0 1 1 0.0 0.0\n0 0 0 2 2 0.0 0.0\n")
            out = os.path.join(d, "uhfk")
            qlms.run(input_dict={
                "log": {"print_level": 0},
                "mode": {"mode": "UHFk",
                         "param": {"Ncond": 112, "IterationMax": 50,
                                   "EPS": 10, "Mix": 0.5, "T": T_RUN,
                                   "CellShape": [NX, NY, 1],
                                   "SubShape": [1, 1, 1]}},
                "file": {"input": {"path_to_input": "",
                                   "interaction": {
                                       "path_to_input": d,
                                       "Geometry": "geom.dat",
                                       "Transfer": "transfer.dat",
                                       "CoulombIntra": "coulombintra.dat"}},
                         "output": {"path_to_output": out,
                                    "rpa": "trans_mod"}}})
            tm = np.load(os.path.join(out, "trans_mod.npz"))["trans_mod"]
        self.assertEqual(tm.shape, (NX * NY, 2 * norb, 2 * norb))
        self.assertLess(np.abs(tm - _trans_mod_from_hr(h)).max(), 1e-12)


# 4x4 random-start UHFk fixture: converged to a paramagnetic state at
# EPS = 6, whose trans_mod keeps a ~1e-5 spin splitting and mixing.
_RAND_NX = 4


def _run_random_uhfk(d):
    import hwave.qlms as qlms
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n1\n0.0 0.0 0.0\n")
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("Transfer\n1\n4\n 1 1 1 1\n"
                "-1 0 0 1 1 -1.0 0.0\n0 -1 0 1 1 -1.0 0.0\n"
                "0 1 0 1 1 -1.0 0.0\n1 0 0 1 1 -1.0 0.0\n")
    with open(os.path.join(d, "coulombintra.dat"), "w") as f:
        f.write("CoulombIntra\n1\n1\n 1\n0 0 0 1 1 1.5 0.0\n")
    inter = {"path_to_input": d, "Geometry": "geom.dat",
             "Transfer": "transfer.dat", "CoulombIntra": "coulombintra.dat"}
    out = os.path.join(d, "uhfk")
    qlms.run(input_dict={
        "log": {"print_level": 0},
        "mode": {"mode": "UHFk",
                 "param": {"Ncond": 12, "IterationMax": 2000, "EPS": 6,
                           "Mix": 0.5, "RndSeed": 1, "T": 0.1,
                           "CellShape": [_RAND_NX, _RAND_NX, 1],
                           "SubShape": [1, 1, 1]}},
        "file": {"input": {"path_to_input": "", "initial_mode": "random",
                           "interaction": dict(inter)},
                 "output": {"path_to_output": out, "rpa": "trans_mod"}}})
    return inter, os.path.join(out, "trans_mod.npz")


def _random_sc_input(d, inter, tm_path, eli_extra=None):
    eli = {"solver_mode": "eigenvalue", "chi0q_mode": "calc",
           "chi0q_tensor": "reduced", "pairing_type": "singlet",
           "init_gap": "d_x2y2", "num_eigenvalues": 2}
    eli.update(eli_extra or {})
    return {"log": {"print_level": 0},
            "mode": {"mode": "RPA",
                     "param": {"T": 0.1, "CellShape": [_RAND_NX, _RAND_NX, 1],
                               "SubShape": [1, 1, 1], "Nmat": 32,
                               "filling": 0.375, "coeff_tail": 1.0}},
            "file": {"input": {"path_to_input": d, "trans_mod": tm_path,
                               "interaction": dict(inter)},
                     "output": {"path_to_output": os.path.join(d, "o")}},
            "eliashberg": eli}


class TestSpinFreeTol(unittest.TestCase):
    """[eliashberg] spin_free_tol: a converged-but-noisy paramagnetic UHF
    output is refused by default and spin-averaged when allowed."""

    def test_random_start_uhfk_output(self):
        with tempfile.TemporaryDirectory() as d:
            inter, tm_path = _run_random_uhfk(d)
            tm = np.load(tm_path)["trans_mod"]
            k = (np.fft.ifftn(tm.reshape(_RAND_NX, _RAND_NX, 1, 2, 2),
                              axes=(0, 1, 2)) * _RAND_NX ** 2).reshape(-1, 2, 2)
            split = np.abs(k[:, 0, 0] - k[:, 1, 1]).max()
            mix = np.abs(k[:, 0, 1]).max()
            # fixture precondition: noisy, but far below a real moment
            self.assertGreater(max(split, mix), 1e-8)
            self.assertLess(max(split, mix), 1e-3)

            with _ForbidDownstream(), self.assertRaises(ValueError) as cm:
                sc.calc_eliashberg(_random_sc_input(d, inter, tm_path))
            msg = str(cm.exception)
            for word in ("#88", "spin_free_tol", "EPS",
                         "{:.3e}".format(split)):
                self.assertIn(word, msg)

            for key in ("spin_free_tol", "Spin_Free_Tol"):
                with self.subTest(key=key):
                    with self.assertLogs("hwave_sc", level="WARNING") as logs:
                        res = _run(_random_sc_input(d, inter, tm_path,
                                                    {key: 1e-3}))
                    self.assertTrue(any("spin-averaged" in m_
                                        for m_ in logs.output), logs.output)
                    self.assertTrue(any("{:.3e}".format(mix) in m_
                                        for m_ in logs.output), logs.output)
                    band = res["prepared"]
                    eps = band.epsilon_k[0, 0].reshape(-1)
                    avg = 0.5 * (k[:, 0, 0] + k[:, 1, 1])
                    self.assertLess(np.abs(eps - avg).max(), 1e-13)
                    # the chi0q solver holds exactly the averaged band
                    self.assertTrue(np.array_equal(
                        np.asarray(band.solver.H0_k)[0, :, 0, 0], eps))

    def test_bad_values_refused(self):
        for val in (-1.0, float("nan"), "abc", True, np.bool_(True), "1e-3",
                    np.str_("1e-3")):
            with self.subTest(val=val):
                with self.assertRaises(ValueError) as cm:
                    sc._read_spin_free_tol({"spin_free_tol": val})
                self.assertIn("spin_free_tol", str(cm.exception))
        self.assertEqual(sc._read_spin_free_tol({}), sc.SPIN_FREE_ATOL)
        self.assertEqual(sc._read_spin_free_tol({"SPIN_FREE_TOL": 2e-3}),
                         2e-3)


class TestBandGuards(unittest.TestCase):

    def test_load_band_into_solver_detects_a_different_band(self):
        class _Ham:
            ham_extern_q = None

        class _Solver:
            ham_info = _Ham()

            def _calc_epsilon_k(self, gi):
                self.spin_mode = "spin-free"
                tm = gi["trans_mod"]
                self.H0_k = tm[np.newaxis, :, :1, :1] + 1e-9   # not `up`

        up = np.zeros((4, 1, 1), complex)
        with self.assertRaises(RuntimeError):
            sc._load_band_into_solver(_Solver(), up)

    def test_wrong_cellshape_trans_mod(self):
        with tempfile.TemporaryDirectory() as d:
            _write_model(d, 1)
            np.savez(os.path.join(d, "tm.npz"),
                     trans_mod=np.zeros((16, 2, 2), complex))
            inp = _input(d, os.path.join(d, "o"),
                         extra_file={"trans_mod": "tm.npz"})
            with _ForbidDownstream(), self.assertRaises(ValueError) as cm:
                sc.calc_eliashberg(inp)
        self.assertIn("does not match", str(cm.exception))

    def test_non_hermitian_band_refused(self):
        h = _nematic_hr(1)
        tm = _trans_mod_from_hr(h)
        tm[1, 0, 0] += 0.01          # R = (0, 1) without its -R partner
        tm[1, 1, 1] += 0.01
        with tempfile.TemporaryDirectory() as d:
            _write_model(d, 1)
            np.savez(os.path.join(d, "tm.npz"), trans_mod=tm)
            inp = _input(d, os.path.join(d, "o"),
                         extra_file={"trans_mod": "tm.npz"})
            with _ForbidDownstream(), self.assertRaises(ValueError) as cm:
                sc.calc_eliashberg(inp)
        self.assertIn("not Hermitian", str(cm.exception))


class TestGreenInitLoadOracle(unittest.TestCase):
    """green_init on the LOAD path: chi0q written by an RPA run with
    green_init, loaded by hwave_sc with the same green_init, gives the
    lambda of a run whose Transfer IS T + Sigma_HF[G] (U only: the on-site
    Hartree shift U n_a per orbital; no Fock term for a spin-free G)."""

    def test_load_green_init_equals_transfer_oracle(self):
        import hwave.qlms as qlms
        norb = 2
        occ = (0.6, 0.3)
        g = np.zeros((NX * NY, 2, norb, 2, norb), dtype=complex)
        for s_ in range(2):
            for a in range(norb):
                g[0, s_, a, s_, a] = occ[a]
        h_hf = _bare_hr(norb).copy()
        for a in range(norb):
            h_hf[0, 0, a, a] += 1.4 * occ[a]
        with tempfile.TemporaryDirectory() as d, \
                tempfile.TemporaryDirectory() as d_or:
            _write_model(d, norb, inter=False)
            np.savez(os.path.join(d, "g.npz"), green=g)
            out = os.path.join(d, "rpa")
            rpa_in = _input(d, out, coulomb_inter=False,
                            extra_file={"green_init": "g.npz"})
            rpa_in["mode"]["calc_scheme"] = "reduced"
            rpa_in["file"]["output"]["chi0q"] = "chi0q"
            del rpa_in["eliashberg"]
            qlms.run(input_dict=rpa_in)
            res = _run(_input(d, out, chi0q_mode="load", coulomb_inter=False,
                              extra_file={"green_init": "g.npz"}))
            _write_model(d_or, norb, transfer_hr=h_hf, inter=False)
            ref = _run(_input(d_or, os.path.join(d_or, "o"),
                              coulomb_inter=False, chi0q_tensor="reduced"))
            bare = _run(_input(d, os.path.join(d, "ob"),
                               coulomb_inter=False, chi0q_tensor="reduced"))
        self.assertAlmostEqual(_leading(res), _leading(ref), delta=1e-9)
        self.assertGreater(abs(_leading(res) - _leading(bare)), 1e-3)


class TestRound5(unittest.TestCase):

    def test_load_mode_refuses_a_band_that_needs_averaging(self):
        # reviewer probe: a trans_mod split of 1e-6 / 1e-4 with
        # spin_free_tol = 1e-3. calc averages; load must refuse before the
        # chi0q file is read (it was built on a different band).
        for delta in (1e-6, 1e-4):
            with self.subTest(delta=delta), \
                    tempfile.TemporaryDirectory() as d:
                h = _nematic_hr(1)
                hd = h.copy()
                hd[0, 0, 0, 0] += delta
                _write_model(d, 1)
                np.savez(os.path.join(d, "tm.npz"),
                         trans_mod=_trans_mod_from_hr(h, hd))
                tol = {"spin_free_tol": 1e-3}
                with _ForbidDownstream(), \
                        self.assertRaises(ValueError) as cm:
                    sc.calc_eliashberg(_input(
                        d, os.path.join(d, "rpa"), chi0q_mode="load",
                        coulomb_inter=False, eliashberg=tol,
                        extra_file={"trans_mod": "tm.npz"}))
                msg = str(cm.exception)
                self.assertIn("chi0q_mode = 'load'", msg)
                self.assertIn("chi0q_mode = 'calc'", msg)
                # calc mode accepts it with the averaged band
                with self.assertLogs("hwave_sc", level="WARNING") as logs:
                    _run(_input(d, os.path.join(d, "c"), coulomb_inter=False,
                                chi0q_tensor="reduced", eliashberg=tol,
                                extra_file={"trans_mod": "tm.npz"}))
                self.assertTrue(any("spin-averaged" in m_
                                    for m_ in logs.output))

    def _transfer_with_asymmetry(self, d, eps):
        h = _bare_hr(1).copy()
        h[1, 0, 0, 0] += eps          # R = (1, 0) without its -R partner
        _write_model(d, 1, transfer_hr=h, extern=True)

    def test_slightly_non_hermitian_transfer_is_hermitized(self):
        with tempfile.TemporaryDirectory() as d:
            self._transfer_with_asymmetry(d, 1e-7)
            with self.assertLogs("hwave_sc", level="WARNING") as logs:
                res = _run(_input(d, os.path.join(d, "o"), extern=True,
                                  extra_param={"coeff_extern": 0.0}))
        self.assertTrue(any("not exactly Hermitian" in m_
                            for m_ in logs.output), logs.output)
        eps = res["prepared"].epsilon_k
        herm = np.abs(eps - np.conj(np.swapaxes(eps, 0, 1))).max()
        self.assertLess(herm, 1e-15)

    def test_grossly_non_hermitian_transfer_is_refused(self):
        with tempfile.TemporaryDirectory() as d:
            self._transfer_with_asymmetry(d, 1e-2)
            with _ForbidDownstream(), self.assertRaises(ValueError) as cm:
                sc.calc_eliashberg(_input(d, os.path.join(d, "o"),
                                          extern=True,
                                          extra_param={"coeff_extern": 0.0}))
        self.assertIn("not Hermitian", str(cm.exception))


def _green_init_vs_uhfk_difference(typ, rows, G, orbital_site=None):
    """max |H_UHFk - H_green_init| over the real-space one-body table:
    UHFk's mean-field Hamiltonian (Hartree + Fock) of G against the band RPA
    assembles from green_init (the production ``RPA._calc_trans_mod`` on
    copies of the solver's arrays); two orbitals. Returns (difference,
    ham_param, parsed orbital_site, grid)."""
    from hwave.qlmsio import read_input_k
    import hwave.solver.uhfk as uhfk
    import hwave.solver.rpa as rpa
    norb, nd = 2, 4
    with tempfile.TemporaryDirectory() as d:
        _write_model(d, norb, inter=False)
        w = np.zeros((NX, NY, norb, norb), complex)
        for (rx, ry, a, b, v) in rows:
            w[rx % NX, ry % NY, a, b] += v
        _write_w90(os.path.join(d, "x.dat"), typ, w)
        fin = {"path_to_input": d,
               "interaction": {"path_to_input": d, "Geometry": "geom.dat",
                               "Transfer": "transfer.dat", typ: "x.dat"}}
        u = uhfk.UHFk(read_input_k.QLMSkInput(fin).get_param("ham"),
                      {"print_level": 0, "print_step": 1},
                      {"mode": "UHFk", "flag_fock": True,
                       "param": {"T": 0.1, "CellShape": [NX, NY, 1],
                                 "SubShape": [1, 1, 1], "Ncond": 2,
                                 "IterationMax": 1}})
        u._make_ham_trans()
        u._make_ham_inter()
        u._detect_blocks()
        u.Green = G.copy()
        u._make_ham()
        param = {"T": 0.1, "CellShape": [NX, NY, 1], "SubShape": [1, 1, 1],
                 "Nmat": 8, "filling": 0.5}
        if orbital_site is not None:
            param["orbital_site"] = list(orbital_site)
        r = rpa.RPA(read_input_k.QLMSkInput(fin).get_param("ham"),
                    {"print_level": 0, "print_step": 1},
                    {"mode": "RPA", "calc_scheme": "general",
                     "param": param})
        hi = r.ham_info
        saved = {n: np.array(getattr(hi, n), copy=True)
                 for n in ("ham_trans_r", "ham_inter_r")}
        try:
            h_rpa = np.array(r._calc_trans_mod(
                np.array(G.reshape(NX * NY, nd, nd), copy=True)))
        finally:
            for n, v in saved.items():
                setattr(hi, n, v)
        ham = read_input_k.QLMSkInput(fin).get_param("ham")
        site = getattr(hi, "orbital_site", None)
        grid = tuple(r.lattice.shape)

    def to_r(hk):
        return (np.fft.fftn(np.asarray(hk).reshape(NX, NY, 1, nd, nd),
                            axes=(0, 1, 2)) / (NX * NY)).reshape(-1, nd, nd)
    return float(np.abs(to_r(u.ham) - to_r(h_rpa)).max()), ham, site, grid


class TestGreenInitIsUHFkMeanField(unittest.TestCase):
    """Issue #278: for every interaction type and row kind the reader
    accepts (with and without orbital_site), the band RPA assembles from
    green_init equals UHFk's T + Hartree + Fock band of the same G. (Before
    #278 the Fock term was left out and hwave_sc warned about it.)"""

    ROW_KINDS = {
        "onsite_intra": [(0, 0, 0, 0, 0.7), (0, 0, 1, 1, 0.7)],
        "onsite_inter": [(0, 0, 0, 1, 0.3), (0, 0, 1, 0, 0.3)],
        "offsite_intra": [(1, 0, 0, 0, 0.3), (-1, 0, 0, 0, 0.3)],
        "offsite_inter": [(1, 0, 0, 1, 0.3), (-1, 0, 1, 0, 0.3)],
    }
    TYPES = ("CoulombIntra", "CoulombInter", "Hund", "Ising", "Exchange",
             "PairLift", "PairHop", "Coulomb")

    @staticmethod
    def _green():
        rng = np.random.default_rng(1)
        G = np.zeros((NX * NY, 2, 2, 2, 2), complex)
        g0 = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        g0 = (g0 + g0.conj().T) / 4 + 0.4 * np.eye(2)
        x = (rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))) * 0.1
        for s_ in range(2):
            G[0, s_, :, s_, :] = g0
            G[1 * NY, s_, :, s_, :] = x                 # R = (1, 0)
            G[(NX - 1) * NY, s_, :, s_, :] = x.conj().T  # R = (-1, 0)
        return G

    def test_band_equals_uhfk_mean_field(self):
        G = self._green()
        checked = 0
        for typ in self.TYPES:
            for kind, rows in self.ROW_KINDS.items():
                for site in (None, (0, 1)):
                    try:
                        diff, _ham, _psite, _grid = \
                            _green_init_vs_uhfk_difference(typ, rows, G, site)
                    except ValueError:
                        continue        # declaration refused by the reader
                    with self.subTest(type=typ, rows=kind, site=site):
                        self.assertLess(diff, 1e-12)
                    checked += 1
        self.assertGreater(checked, 50)

    def test_spin_flip_green_is_refused_for_every_type(self):
        """A spin-flip (in-plane) G is no longer a warning case: the
        G-level check refuses it for every type and row kind."""
        G = self._green()
        G[0, 0, 0, 1, 0] = G[0, 1, 0, 0, 0] = 0.1
        G[0, 0, 1, 1, 1] = G[0, 1, 1, 0, 1] = 0.05
        for typ in self.TYPES:
            with self.subTest(type=typ), tempfile.TemporaryDirectory() as d:
                rows = (self.ROW_KINDS["onsite_intra"]
                        if typ in ("CoulombIntra", "Coulomb")
                        else self.ROW_KINDS["offsite_inter"])
                _write_model(d, 2, inter=False)
                w = np.zeros((NX, NY, 2, 2), complex)
                for (rx, ry, a, b, v) in rows:
                    w[rx % NX, ry % NY, a, b] += v
                _write_w90(os.path.join(d, "x.dat"), typ, w)
                np.savez(os.path.join(d, "g.npz"), green=G)
                inp = _input(d, os.path.join(d, "o"), coulomb_inter=False,
                             extra_file={"green_init": "g.npz"})
                files_ = inp["file"]["input"]["interaction"]
                files_.pop("CoulombIntra")
                files_[typ] = "x.dat"
                with self.assertRaises(ValueError) as cm:
                    sc._mean_field_band(inp, ("green_init",), "general",
                                        (NX, NY, 1), 2)
                self.assertIn("green_init describes a spin-dependent state",
                              str(cm.exception))

    def test_no_fock_warning_through_the_full_path(self):
        """The green_init band carries the Fock term, so hwave_sc no longer
        warns about Fock-carrying interactions (Hund, off-site PairHop)."""
        g = np.zeros((NX * NY, 2, 2, 2, 2), dtype=complex)
        for s_ in range(2):
            g[0, s_, 0, s_, 0] = 0.6
            g[0, s_, 1, s_, 1] = 0.3
            g[0, s_, 0, s_, 1] = g[0, s_, 1, s_, 0] = 0.1
        for typ, rows, site in (
                ("Hund", [(0, 0, 0, 1, 0.3), (0, 0, 1, 0, 0.3)], None),
                ("PairHop", [(0, 0, 0, 1, 0.3), (0, 0, 1, 0, 0.3)], [0, 1])):
            with self.subTest(type=typ), tempfile.TemporaryDirectory() as d:
                _write_model(d, 2, inter=False)
                w = np.zeros((NX, NY, 2, 2), complex)
                for (rx, ry, a, b, v) in rows:
                    w[rx % NX, ry % NY, a, b] += v
                _write_w90(os.path.join(d, "x.dat"), typ, w)
                np.savez(os.path.join(d, "g.npz"), green=g)
                inp = _input(d, os.path.join(d, "o"), coulomb_inter=False,
                             extra_file={"green_init": "g.npz"})
                inp["file"]["input"]["interaction"][typ] = "x.dat"
                if site is not None:
                    inp["mode"]["param"]["orbital_site"] = site
                records = []
                handler = logging.Handler()
                handler.emit = records.append
                lg = logging.getLogger("hwave_sc")
                lg.addHandler(handler)
                try:
                    sc._mean_field_band(inp, ("green_init",), "general",
                                        (NX, NY, 1), 2)
                finally:
                    lg.removeHandler(handler)
                self.assertFalse(any("Fock" in r.getMessage()
                                     for r in records))


class TestHermiticityRound6(unittest.TestCase):

    def _model(self, d, asym, offset=0.0, norb=1, a=0, b=0):
        h = _bare_hr(norb).copy()
        h[1, 0, a, b] += asym         # R = (1, 0) without its -R partner
        for o in range(norb):
            h[0, 0, o, o] += offset
        _write_model(d, norb, transfer_hr=h, extern=True)

    def test_load_mode_does_not_hermitize(self):
        with tempfile.TemporaryDirectory() as d:
            self._model(d, 1e-7)
            with _ForbidDownstream(), self.assertRaises(ValueError) as cm:
                sc.calc_eliashberg(_input(d, os.path.join(d, "o"),
                                          chi0q_mode="load", extern=True,
                                          extra_param={"coeff_extern": 0.0}))
        msg = str(cm.exception)
        self.assertIn("not exactly Hermitian", msg)
        self.assertIn("chi0q_mode = 'calc'", msg)

    def test_bound_is_relative(self):
        # asymmetry 1e-5 (|H - H^dagger| up to 2e-5): accepted next to a
        # 1e3 on-site offset (bound 1e-6 * max|H| ~ 1e-3), refused without
        # it (bound ~ 5e-6)
        with tempfile.TemporaryDirectory() as d:
            self._model(d, 1e-5, offset=1.0e3)
            with self.assertLogs("hwave_sc", level="WARNING") as logs:
                _run(_input(d, os.path.join(d, "o"), extern=True,
                            extra_param={"coeff_extern": 0.0}))
        self.assertTrue(any("not exactly Hermitian" in m_
                            for m_ in logs.output))
        with tempfile.TemporaryDirectory() as d:
            self._model(d, 1e-5)
            with _ForbidDownstream(), self.assertRaises(ValueError) as cm:
                sc.calc_eliashberg(_input(d, os.path.join(d, "o"),
                                          extern=True,
                                          extra_param={"coeff_extern": 0.0}))
        self.assertIn("not Hermitian", str(cm.exception))

    def test_hermitized_band_is_exactly_the_hermitian_part(self):
        # two orbitals, one-sided inter-orbital hopping at R = (1, 0)
        import hwave.qlmsio.read_input_k as rik
        import hwave.solver.rpa as rpa
        norb = 2
        with tempfile.TemporaryDirectory() as d:
            self._model(d, 1e-7, norb=norb, a=0, b=1)
            inp = _input(d, os.path.join(d, "o"), extern=True,
                         extra_param={"coeff_extern": 0.0})
            band = sc._mean_field_band(inp, ("Extern",), "general",
                                       (NX, NY, 1), norb)
            info = {"path_to_input": d,
                    "interaction": inp["file"]["input"]["interaction"]}
            fresh = rpa.RPA(rik.QLMSkInput(info).get_param("ham"),
                            {"print_level": 0},
                            {"mode": "RPA", "calc_scheme": "general",
                             "param": dict(inp["mode"]["param"])})
        H = np.asarray(fresh.ham_info.ham_trans_q)       # (nvol, a, b)
        expected = 0.5 * (H + np.conj(np.swapaxes(H, -1, -2)))
        eps = np.transpose(band.epsilon_k[:, :, :, :, 0],
                           (2, 3, 0, 1)).reshape(-1, norb, norb)
        self.assertTrue(np.array_equal(eps, expected))
        self.assertGreater(np.abs(H - expected).max(), 1e-8)


class TestGreenInitSpinCheck(unittest.TestCase):
    """A magnetic green_init is refused whatever the spin direction."""

    def _call(self, g, tol=None):
        with tempfile.TemporaryDirectory() as d:
            _write_model(d, 1, inter=False)
            np.savez(os.path.join(d, "g.npz"), green=g)
            inp = _input(d, os.path.join(d, "o"), coulomb_inter=False,
                         extra_file={"green_init": "g.npz"})
            kw = {} if tol is None else {"spin_free_tol": tol}
            return sc._mean_field_band(inp, ("green_init",), "general",
                                       (NX, NY, 1), 1, **kw)

    @staticmethod
    def _g(mz=0.0, mx=0.0):
        g = np.zeros((NX * NY, 2, 1, 2, 1), dtype=complex)
        g[0, 0, 0, 0, 0] = FILLING + mz
        g[0, 1, 0, 1, 0] = FILLING - mz
        g[0, 0, 0, 1, 0] = g[0, 1, 0, 0, 0] = mx
        return g

    def test_in_plane_and_z_moments_refused_paramagnet_accepted(self):
        for label, g in (("in-plane", self._g(mx=0.2)),
                         ("z", self._g(mz=0.2))):
            with self.subTest(moment=label):
                with self.assertRaises(ValueError) as cm:
                    self._call(g)
                msg = str(cm.exception)
                for word in ("#88", "trans_mod", "spin_free_tol",
                             "dimensionless"):
                    self.assertIn(word, msg)
        band = self._call(self._g())
        eps = np.transpose(band.epsilon_k[:, :, :, :, 0], (2, 3, 0, 1))
        self.assertLess(np.abs(eps - (_band_k(_bare_hr(1))
                                      + 1.4 * FILLING)).max(), 1e-12)

    def test_loosened_tolerance_is_consistent(self):
        g = self._g(mz=0.5e-6)          # |G_up - G_down| = 1e-6
        with self.assertRaises(ValueError):
            self._call(g)
        with self.assertLogs("hwave_sc", level="WARNING") as logs:
            band = self._call(g, tol=1e-3)
        # accepted, and the band check then spin-averages the Hartree band
        self.assertTrue(any("spin-averaged" in m_ for m_ in logs.output))
        eps = np.transpose(band.epsilon_k[:, :, :, :, 0], (2, 3, 0, 1))
        self.assertLess(np.abs(eps - (_band_k(_bare_hr(1))
                                      + 1.4 * FILLING)).max(), 1e-12)


class TestHermitianRtolPinned(unittest.TestCase):
    """The bound is exactly HERMITIAN_RTOL * max|H| next to a large
    offset: an asymmetry just below it is accepted, just above refused."""

    def _assembled(self, d, asym):
        import hwave.qlmsio.read_input_k as rik
        import hwave.solver.rpa as rpa
        h = _bare_hr(1).copy()
        h[1, 0, 0, 0] += asym
        h[0, 0, 0, 0] += 1.0e3
        _write_model(d, 1, transfer_hr=h, extern=True)
        inp = _input(d, os.path.join(d, "o"), extern=True,
                     extra_param={"coeff_extern": 0.0})
        info = {"path_to_input": d,
                "interaction": inp["file"]["input"]["interaction"]}
        r = rpa.RPA(rik.QLMSkInput(info).get_param("ham"), {"print_level": 0},
                    {"mode": "RPA", "calc_scheme": "general",
                     "param": dict(inp["mode"]["param"])})
        H = np.asarray(r.ham_info.ham_trans_q)
        herm = np.abs(H - np.conj(np.swapaxes(H, -1, -2))).max()
        bound = 1.0e-6 * np.abs(H).max()     # the documented relative bound
        return inp, herm, bound

    def test_just_below_and_just_above(self):
        with tempfile.TemporaryDirectory() as d:
            _inp, herm1, bound1 = self._assembled(d, 1e-4)
        # herm is linear in the asymmetry: herm = c * asym
        c = herm1 / 1e-4
        for factor, accepted in ((0.99, True), (1.01, False)):
            with self.subTest(factor=factor), \
                    tempfile.TemporaryDirectory() as d:
                inp, herm, bound = self._assembled(d, factor * bound1 / c)
                self.assertEqual(herm < bound, accepted)   # precondition
                if accepted:
                    _run(inp)
                else:
                    with _ForbidDownstream(), \
                            self.assertRaises(ValueError) as cm:
                        sc.calc_eliashberg(inp)
                    self.assertIn("not Hermitian", str(cm.exception))


class TestGreenInitGuardsPinned(unittest.TestCase):

    def _call(self, g, tol=None):
        with tempfile.TemporaryDirectory() as d:
            _write_model(d, 1, inter=False)
            np.savez(os.path.join(d, "g.npz"), green=g)
            inp = _input(d, os.path.join(d, "o"), coulomb_inter=False,
                         extra_file={"green_init": "g.npz"})
            kw = {} if tol is None else {"spin_free_tol": tol}
            return sc._mean_field_band(inp, ("green_init",), "general",
                                       (NX, NY, 1), 1, **kw)

    @staticmethod
    def _g():
        g = np.zeros((NX * NY, 2, 1, 2, 1), dtype=complex)
        g[0, 0, 0, 0, 0] = g[0, 1, 0, 1, 0] = FILLING
        return g

    def test_nan_in_a_cell_rpa_does_not_read_is_refused(self):
        g = self._g()
        g[1, 0, 0, 0, 0] = np.nan          # R != 0: not in the Hartree band
        with self.assertRaises(ValueError) as cm:
            self._call(g)
        self.assertIn("non-finite", str(cm.exception))

    def test_g_tolerance_is_floored_at_1e_10(self):
        # spin_free_tol = 0 tightens the band check, not the G check: G
        # noise of 5e-11 (in a cell the Hartree band does not read, so the
        # band itself is exactly spin-free) is still accepted
        g = self._g()
        g[1, 0, 0, 0, 0] = 0.1 + 5e-11
        g[1, 1, 0, 1, 0] = 0.1
        self.assertIsNotNone(self._call(g, tol=0.0))
        g[1, 0, 0, 0, 0] = 0.1 + 5e-10   # above the floor: refused
        with self.assertRaises(ValueError) as cm:
            self._call(g, tol=0.0)
        self.assertIn("green_init describes a spin-dependent state",
                      str(cm.exception))

if __name__ == "__main__":
    unittest.main()
