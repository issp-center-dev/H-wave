"""hwave_sc on a mean-field one-body band (issue #273).

``[file.input] trans_mod`` (the T + Sigma_HF that UHFk writes), ``green_init``
(T + HF[G]) and an active ``Extern`` field replace the bare Transfer band in
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
                 prepared=None):
            self.data["mu_chi0"] = precomputed_mu
            self.data["prepared"] = prepared
            return orig(input_dict, chi0q_tensor=chi0q_tensor,
                        precomputed_mu=precomputed_mu, prepared=prepared)

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

        # mu handed to chi0q == RPA._find_mu on the trans_mod band
        band = res["prepared"]
        self.assertIsNotNone(band)
        solver = band.solver
        _dist, mu_rpa = solver._find_mu(solver.Ncond / 2, solver.T)
        self.assertLessEqual(abs(res["mu_chi0"] - mu_rpa), 1e-12)

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


class TestGreenInit(unittest.TestCase):
    """A UHFk solution supplied as green_init defines the same band as the
    trans_mod UHFk saved with it.

    On-site interactions only: the RPA reader turns green_init into
    T + HF[G] from the on-site (R = 0) Green function alone, so an off-site
    Fock term of a CoulombInter run is present in UHFk's trans_mod but not
    in the green_init band (the two inputs then define different bands, in
    RPA as well). Two orbitals, so the orbital-dependent Hartree shift
    changes the band beyond a constant that mu would absorb.
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
            with self.assertRaises(ValueError) as cm:
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

    def test_inactive_extern_keeps_bare_path(self):
        # coeff_extern = 0 adds nothing to the band: no RPA band is built
        self.assertEqual(sc._one_body_band_sources(
            {"file": {"input": {"interaction": {"extern": "e.dat"}}},
             "mode": {"param": {"coeff_extern": 0.0}}}), ())
        self.assertEqual(sc._one_body_band_sources(
            {"file": {"input": {"interaction": {"Extern": "e.dat"}}},
             "mode": {"param": {"Coeff_Extern": 0.5}}}), ("Extern",))


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

    def test_green_init_only_onsite_cell_enters(self):
        # spin-free g0[0], spin-dependent cell 1: the RPA band uses g0[0]
        # only, so the band is spin-free and accepted
        g = np.zeros((NX * NY, 2, 1, 2, 1), dtype=complex)
        g[0, 0, 0, 0, 0] = g[0, 1, 0, 1, 0] = FILLING
        g[1, 0, 0, 0, 0] = 0.2
        g[1, 1, 0, 1, 0] = -0.2
        band = self._band(extra_file={"green_init": "g.npz"},
                          files={"g.npz": {"green": g}})
        self.assertIsNotNone(band)

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
        """Per interaction type: a spin-free G gives an accepted band; a ferro
        G gives a refused band whenever RPA's T + HF[G] carries its spin
        dependence (two orbitals, on-site terms).

        RPA's green_init band is T plus the on-site HF term its
        ``_calc_trans_mod`` builds. For CoulombInter, Exchange and PairLift
        that term carries no spin dependence for the ferro patterns below
        (nor for an on-site spin-flip G), so those bands ARE spin-free and
        are accepted -- the same band RPA's own chi0q uses. This pins the
        current RPA behaviour; it is not a statement about the UHF band of
        such a G.
        """
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
        refused = {"CoulombIntra": True, "Hund": True, "Ising": True,
                   "PairHop": True, "CoulombInter": False, "Exchange": False,
                   "PairLift": False}
        for tname, ferro_refused in refused.items():
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
                                     ("ferro", g_ferro, not ferro_refused)):
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
                        self.assertIn("#88", str(cm.exception))


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
        self.assertGreater(abs(a - bare), 1e-6)


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

if __name__ == "__main__":
    unittest.main()
