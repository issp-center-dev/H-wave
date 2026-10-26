"""The static spin-orbital entry of hwave_sc (``[mode] enable_spin_orbital``,
issue #83): ``chi0q_mode = "calc"``, ``"load"`` and ``"flex"``.

The gap is one ``2N x 2N`` spin-orbital matrix ``Delta_AB(k)`` (spin-block
order ``A = s * N + a``); the kernel is ``P K`` with the pairing vertex of
``hwave.solver.so_pairing`` and ``P`` the projection onto fermion-
antisymmetric gaps ``Delta_AB(k) = -Delta_BA(-k)``, applied after every
kernel application.

* ``TestStaticEntry``: the SU(2) limit reproduces the non-spin-orbital
  singlet and (threefold) triplet eigenvalues end to end; a model with
  complex Rashba-like spin-orbit hopping runs and returns an antisymmetric
  gap; ``chi0q_mode = "load"`` of a spin-orbital RPA chi0q equals
  ``chi0q_mode = "calc"`` (input agreement).
* ``TestFlexMode``: ``chi0q_mode = "flex"`` uses the stored chi (zero bosonic
  frequency) and dressed green of a spin-orbital FLEX run (input
  agreement); the archive is identified by its tags (nd = 2 has the shape of a two-orbital
  non-spin-orbital file); non-spin-orbital, reduced and untagged outputs and
  a vertex the FLEX run did not use are refused.
* ``TestRefusals``: every unsupported combination is refused through the
  public entry, with the key named.
* ``TestFlexVertexCheck``: the FLEX vertex-consistency check is a per-q
  backward error (independent of the lattice size; zero-scale blocks).
* ``TestAntisymmetricOperator``: ``matmat`` of the projected kernel equals
  stacked ``matvec`` and returns antisymmetric gaps.
* ``TestArchiveMembers``: only spin-orbital runs write the ``spin_orbital``
  marker; non-spin-orbital RPA / FLEX archives keep their member sets.
* ``TestTransferUnderSOC``: with spin-orbit hopping, the epsilon(k) of
  hwave_sc equals the H0(k) of the spin-orbital RPA.
* ``TestInitialGap``: the spin structure of the spin-orbital seed.
* ``TestSpinDecomposition``: ``so_pairing.so_spin_components`` /
  ``so_spin_decomposition`` on hand-built gaps (pure psi, d_x, d_y, d_z and a
  mixture), written out as explicit 2 x 2 spin matrices.
* ``TestOutputs``: the spin-orbital output schema (eigenvalue.dat header and
  weight columns, gap.dat comments, the gap.npz archive) end to end in the
  SU(2) limit.
"""
import contextlib
import logging
import os
import shutil
import tempfile
import unittest
from unittest import mock

import numpy as np

import hwave.sc as sc
from hwave.solver import so_pairing

#: physical orbitals
NP = 2
ND = 2 * NP
#: the lattice
CELL = [4, 4, 1]

#: the rotation-invariant on-site set (Kanamori form, as in
#: tests/test_flex_spin_orbital.py): J_Exchange = J_Hund - 2 J_Ising
_J_HUND, _J_ISING = 0.12, 0.04
_SU2_VALUES = {"CoulombIntra": None, "CoulombInter": 0.5, "Hund": _J_HUND,
               "Exchange": _J_HUND - 2 * _J_ISING, "Ising": _J_ISING,
               "PairHop": 0.05}
_SU2_TYPES = tuple(_SU2_VALUES)

_IY = np.array([[0, 1], [-1, 0]], complex)


@contextlib.contextmanager
def _quiet():
    prev = logging.root.manager.disable
    logging.disable(logging.WARNING)
    try:
        yield
    finally:
        logging.disable(prev)


def _wfile(path, head, n, rows):
    with open(path, "w") as f:
        f.write("{}\n{}\n1\n 1\n".format(head, n))
        for r in rows:
            f.write("  %d %d %d %d %d  %.12f  %.12f\n" % r)


def _write_model(d, so, soc=0.0, types=_SU2_TYPES, extra_rows=None,
                 n_geom=None):
    """A centrosymmetric two-orbital square lattice (real, symmetric
    hopping) plus, in the spin-orbital encoding, the Rashba-like hopping
    ``soc (sin ky sigma_x - sin kx sigma_y)`` on both orbitals.

    ``so``: geometry of ``2 NP`` spin-orbitals and the transfer on the
    interleaved index ``2 a + s``; the interaction files use the PHYSICAL
    orbital index either way. ``extra_rows``: ``{type: [(R, a, b, v)]}``
    appended to the interaction files. Returns the interaction file map."""
    n = 2 * NP if so else NP
    with open(os.path.join(d, "geom.dat"), "w") as f:
        ng = n if n_geom is None else n_geom
        f.write("1 0 0\n0 1 0\n0 0 1\n%d\n" % ng)
        f.write("0 0 0\n" * ng)
    t = {(0, 0): -1.0, (1, 1): -0.7, (0, 1): 0.15, (1, 0): 0.15}
    eps = {0: 0.0, 1: 0.3}
    hop = []   # (R, a, b, sa, sb, value)
    for (a, b), v in t.items():
        for R in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            for s in range(2):
                hop.append((R, a, b, s, s, v))
    for a, e in eps.items():
        for s in range(2):
            hop.append(((0, 0), a, a, s, s, e))
    if so and soc:
        sx = np.array([[0, 1], [1, 0]], complex)
        sy = np.array([[0, -1j], [1j, 0]], complex)
        # t(+x) = i soc sigma_y, t(+y) = -i soc sigma_x and their Hermitian
        # partners: eps(k) = soc (sin ky sigma_x - sin kx sigma_y)
        for R, m in (((1, 0), 1j * soc * sy), ((-1, 0), -1j * soc * sy),
                     ((0, 1), -1j * soc * sx), ((0, -1), 1j * soc * sx)):
            for a in range(NP):
                for s in range(2):
                    for sp in range(2):
                        if m[s, sp] != 0:
                            hop.append((R, a, a, s, sp, m[s, sp]))
    rows = []
    for R, a, b, s, sp, v in hop:
        if so:
            i, j = 2 * a + s + 1, 2 * b + sp + 1
        elif s == 0 and sp == 0:
            i, j = a + 1, b + 1
        else:
            continue
        v = complex(v)
        rows.append((R[0], R[1], 0, i, j, v.real, v.imag))
    _wfile(os.path.join(d, "transfer.dat"), "Transfer", n, rows)
    offd = [(a, b) for a in range(NP) for b in range(NP) if a != b]
    files = {}
    for nm in types:
        if nm == "CoulombIntra":
            spec = [(0, a, a, 1.6 + 0.3 * a) for a in range(NP)]
        else:
            spec = [(0, a, b, _SU2_VALUES[nm]) for a, b in offd]
        spec = spec + list((extra_rows or {}).get(nm, ()))
        fn = nm.lower() + ".dat"
        _wfile(os.path.join(d, fn), nm, NP,
               [(R, 0, 0, a + 1, b + 1, v, 0.0) for R, a, b, v in spec])
        files[nm] = fn
    return files


def _input(d, files, so, eli, T=0.25, nmat=64, out="out", **param):
    par = {"T": T, "CellShape": list(CELL), "SubShape": [1, 1, 1],
           "Nmat": nmat, "filling": 0.4}
    par.update(param)
    inter = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
    inter.update(files)
    mode = {"mode": "RPA", "param": par}
    if so:
        mode["enable_spin_orbital"] = True
    return {"mode": mode,
            "file": {"input": {"path_to_input": "", "interaction": inter},
                     "output": {"path_to_output": os.path.join(d, out),
                                "chi0q": "chi0q"}},
            "eliashberg": dict(eli)}


def _run(inp):
    """``calc_eliashberg`` with the arguments of ``_save_results`` captured at
    full precision: ``(gap, eigenvalue, eigenvalues, match, provenance)``."""
    seen = {}
    real = sc._save_results
    real_reorder = sc._reorder_eigenpairs_by_parity

    def spy_reorder(*a, **kw):
        out = real_reorder(*a, **kw)
        seen["vecs"] = out[1]
        return out

    def spy(output_dir, sigma, eigenvalue, eigenvalues_eig, *a, **kw):
        seen.update(gap=sigma, ev=eigenvalue, evs=eigenvalues_eig,
                    match=kw.get("eigenvalue_match"),
                    provenance=kw.get("provenance"))
        return real(output_dir, sigma, eigenvalue, eigenvalues_eig, *a, **kw)

    with _quiet(), mock.patch.object(sc, "_save_results", side_effect=spy), \
            mock.patch.object(sc, "_reorder_eigenpairs_by_parity",
                              side_effect=spy_reorder):
        sc.calc_eliashberg(inp)
    return seen


def _spin_weights(gap, n):
    """``(|psi|^2, |d_x|^2, |d_y|^2, |d_z|^2) / |Delta|^2`` of a spin-block
    gap, ``Delta = (psi + d . sigma) i sigma_y`` per orbital pair."""
    g = np.asarray(gap).reshape((2, n, 2, n) + gap.shape[2:])
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    sz = np.diag([1.0, -1.0]).astype(complex)
    basis = [_IY, sx @ _IY, sy @ _IY, sz @ _IY]
    tot = np.linalg.norm(g) ** 2
    return [float(np.sum(np.abs(np.einsum("st,satb...->ab...", m.conj(), g)) ** 2)
                  / 2.0 / tot) for m in basis]


def _reverse(gap):
    return sc._reverse_k_and_orbital(gap)


def _clusters(vals, vecs, rtol=1e-8):
    """Group eigenvalues equal to ``rtol`` and give each cluster's even
    and odd dimension (k -> -k with the orbital transpose): for a
    degenerate cluster the solver may return parity-mixed eigenvectors,
    so the dimensions are the ranks of the projected span."""
    order = np.argsort(vals.real)
    out = []
    for i in order:
        if out and abs(vals[i] - out[-1][0]) <= rtol * abs(out[-1][0]):
            out[-1][1].append(i)
        else:
            out.append([vals[i], [i]])
    res = []
    for v, idx in out:
        V = np.stack([vecs[i].ravel() for i in idx], axis=1)
        dims = []
        for ch in ("singlet", "triplet"):
            P = np.stack([sc._project_gap_parity(vecs[i], ch).ravel()
                          for i in idx], axis=1)
            sv = np.linalg.svd(P, compute_uv=False)
            dims.append(int(np.sum(sv > 1e-6 * np.linalg.norm(V, 2))))
        res.append((v, len(idx), dims[0], dims[1]))
    return res


def _assert_su2_multiplicities(tc, res, so, rtol=1e-8, forbidden_floor=0.5):
    """SU(2) oracle: every allowed singlet eigenvalue of the non-SO kernels
    ``res["singlet"]`` (even sector) appears once and every allowed triplet
    eigenvalue of ``res["triplet"]`` (odd sector) three times in the
    spin-orbital spectrum ``so``; the large Pauli-forbidden triplet modes
    (above ``forbidden_floor``) are absent. ``res[ch]`` / ``so`` are :func:`_run` captures; the non-SO
    sets must cover all but the smallest-magnitude eigenvalues."""
    evs = so["evs"]
    scale = np.abs(evs).max()
    nonzero = np.abs(evs) > 1e-10 * scale
    tc.assertTrue(np.all(so["match"][nonzero]))
    # expected multiplicity of each value: even singlet x 1 + odd triplet x 3
    expected, forbidden = {}, []
    floor = {ch: np.abs(res[ch]["evs"]).min() for ch in res}
    for ch, allowed_even in (("singlet", True), ("triplet", False)):
        for v, n, d_even, d_odd in _clusters(res[ch]["evs"], res[ch]["vecs"],
                                             rtol):
            tc.assertEqual(d_even + d_odd, n)
            if abs(v) <= (1 + 1e-6) * floor[ch]:
                continue        # a cluster cut by the truncation
            d_ok, d_bad = (d_even, d_odd) if allowed_even else (d_odd, d_even)
            key = next((k for k in expected
                        if abs(k - v) <= rtol * abs(v)), v)
            expected[key] = expected.get(key, 0) + d_ok * (
                1 if ch == "singlet" else 3)
            if d_bad:
                forbidden.append((ch, v))
    tc.assertGreater(len(expected), 20)
    for v, count in expected.items():
        got = int(np.sum(np.abs(evs - v) <= rtol * abs(v)))
        with tc.subTest(value=v):
            tc.assertEqual(got, count)
    # the large Pauli-forbidden triplet modes (even k-parity) are absent
    big = sorted((v for ch, v in forbidden
                  if ch == "triplet" and v.real > forbidden_floor),
                 key=lambda x: -x.real)
    tc.assertGreaterEqual(len(big), 2)
    for v in big:
        with tc.subTest(forbidden=v):
            if any(abs(k - v) <= rtol * abs(v) for k in expected):
                continue
            tc.assertFalse(np.any(np.abs(evs - v) <= 1e-6 * abs(v)))
    return expected


class _TmpDir(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="sc_so_")
        self.addCleanup(shutil.rmtree, self.d, True)


class TestStaticEntry(_TmpDir):

    def _leading_physical(self, seen):
        evs, match = seen["evs"], seen["match"]
        self.assertTrue(match[0])
        return evs[0]

    def test_su2_model_reproduces_singlet_and_triplet(self):
        """SU(2) limit: with spin-independent hopping and a rotation-invariant
        interaction, every singlet eigenvalue (even sector of the non-SO
        singlet kernel) appears once in the spin-orbital spectrum and every
        triplet eigenvalue (odd sector of the non-SO triplet kernel) three
        times; the Pauli-forbidden modes (odd singlet, even triplet) are
        absent. Both non-SO sets cover 62 of the 64 eigenvalues, the SO set
        254 of 256 (every nonzero eigenvalue: the symmetric sector sits at
        zero)."""
        eli = {"chi0q_mode": "calc", "solver_mode": "eigenvalue"}
        res = {}
        for ch in ("singlet", "triplet"):
            files = _write_model(self.d, so=False)
            res[ch] = _run(_input(self.d, files, False,
                                  dict(eli, pairing_type=ch, num_eigenvalues=62),
                                  out="out_" + ch))
        files = _write_model(self.d, so=True)
        so = _run(_input(self.d, files, True, dict(eli, num_eigenvalues=254),
                         out="out_so"))
        _assert_su2_multiplicities(self, res, so)

    def test_complex_soc_model_runs(self):
        files = _write_model(self.d, so=True, soc=0.3)
        seen = _run(_input(self.d, files, True,
                           {"chi0q_mode": "calc", "solver_mode": "eigenvalue",
                            "num_eigenvalues": 8}))
        # every computed eigenpair is fermion-antisymmetric: the symmetric
        # (Pauli-forbidden) sector is projected out of the kernel
        self.assertTrue(np.all(seen["match"]))
        lam = seen["evs"][0]
        self.assertTrue(np.isfinite(lam))
        self.assertGreater(abs(lam.real), 1e-3)
        self.assertLess(abs(lam.imag), 1e-10)
        gap = seen["gap"]
        self.assertEqual(gap.shape, (ND, ND) + tuple(CELL))
        np.testing.assert_allclose(gap, -_reverse(gap), rtol=0,
                                   atol=1e-12 * np.abs(gap).max())
        # spin-orbit coupling mixes the spin structure of the gap
        w = _spin_weights(gap, NP)
        self.assertAlmostEqual(sum(w), 1.0, places=12)
        self.assertGreater(min(w[0], max(w[1:])), 1e-6)
        self.assertEqual(seen["provenance"]["spin_orbital"], "true")
        self.assertEqual(seen["provenance"]["index_convention"], "spin_block")

    def test_iteration_matches_eigenvalue(self):
        """The shifted power iteration (with the projection in every step)
        converges to the leading (largest real part) eigenpair of the Arnoldi
        solve."""
        files = _write_model(self.d, so=True, soc=0.3)
        seen = _run(_input(self.d, files, True,
                           {"chi0q_mode": "calc", "solver_mode": "both",
                            "num_eigenvalues": 8, "convergence_tol": 1e-10,
                            "max_iter": 20000, "spectral_shift": "auto"}))
        lam = seen["evs"][0]
        self.assertGreater(lam.real, 0.0)
        self.assertAlmostEqual(seen["ev"], lam.real, delta=1e-7 * abs(lam))
        gap = seen["gap"]
        np.testing.assert_allclose(gap, -_reverse(gap), rtol=0,
                                   atol=1e-12 * np.abs(gap).max())

    def test_calc_and_load_agree(self):
        """Input agreement: the general chi0q of a spin-orbital RPA run, loaded,
        gives the eigenvalues of chi0q_mode = "calc"."""
        from hwave import qlms
        files = _write_model(self.d, so=True, soc=0.3)
        eli = {"solver_mode": "eigenvalue", "num_eigenvalues": 8}
        calc = _run(_input(self.d, files, True, dict(eli, chi0q_mode="calc"),
                           out="out_calc"))
        rpa_in = _input(self.d, files, True, {}, out="out_load")
        rpa_in["mode"]["calc_scheme"] = "general"
        rpa_in["file"]["input"]["path_to_input"] = self.d
        del rpa_in["eliashberg"]
        rpa_in["log"] = {"print_level": 0}
        with _quiet():
            qlms.run(input_dict=rpa_in)
        load = _run(_input(self.d, files, True, dict(eli, chi0q_mode="load"),
                           out="out_load"))
        np.testing.assert_allclose(load["evs"], calc["evs"], rtol=1e-10, atol=0)


class TestRefusals(_TmpDir):
    """Refusals: each unsupported combination is refused through the public entry
    with the key named."""

    def _refused(self, pattern, eli=None, so_files=None, **kw):
        files = so_files if so_files is not None else _write_model(self.d, so=True)
        inp = _input(self.d, files, True,
                     dict({"chi0q_mode": "calc", "solver_mode": "eigenvalue"},
                          **(eli or {})), **kw)
        with self.assertRaisesRegex(ValueError, pattern):
            with _quiet():
                sc.calc_eliashberg(inp)
        out = inp["file"]["output"]["path_to_output"]
        self.assertFalse(os.path.exists(out) and os.listdir(out))

    def test_dynamic(self):
        self._refused(r"frequency", {"frequency": "dynamic", "chi0q_mode": "flex"})

    def test_dynamic_case_insensitive_key(self):
        self._refused(r"frequency", {"Frequency": "dynamic"})

    def test_bond_channels(self):
        self._refused(r"bond_channels", {"bond_channels": True})

    def test_offsite_coulomb_inter(self):
        files = _write_model(self.d, so=True, extra_rows={
            "CoulombInter": [(1, 0, 0, 0.2), (-1, 0, 0, 0.2)]})
        self._refused(r"off-site row .* in the CoulombInter file", so_files=files)

    def test_offsite_row_in_coulomb_aggregate(self):
        _write_model(self.d, so=True, types=())
        _wfile(os.path.join(self.d, "coulomb.dat"), "Coulomb", NP,
               [(0, 0, 0, 1, 1, 1.5, 0.0), (0, 0, 0, 2, 2, 1.5, 0.0),
                (1, 0, 0, 1, 1, 0.3, 0.0), (-1, 0, 0, 1, 1, 0.3, 0.0)])
        self._refused(r"off-site row .* in the Coulomb file", so_files={"Coulomb": "coulomb.dat"})

    def test_cross_site_row_under_orbital_site(self):
        files = _write_model(self.d, so=True, types=("CoulombIntra", "CoulombInter"))
        self._refused(r"orbital_site", so_files=files, orbital_site=[0, 1])

    def test_reduced_chi0q_tensor(self):
        self._refused(r"chi0q_tensor", {"chi0q_tensor": "reduced"})

    def test_reduced_chi0q_file(self):
        out = os.path.join(self.d, "out")
        os.makedirs(out)
        nvol = int(np.prod(CELL))
        np.savez(os.path.join(out, "chi0q.npz"),
                 chi0q=np.zeros((64, nvol, ND, ND), complex),
                 index_convention="spin_block", momentum_convention="e_plus_ikR",
                 spin_orbital=True)
        files = _write_model(self.d, so=True)
        inp = _input(self.d, files, True, {"chi0q_mode": "load",
                                           "solver_mode": "eigenvalue"})
        # refused right after loading, before the Green function is built
        with self.assertRaisesRegex(ValueError, r"general \(4-index\)"):
            with _quiet(), mock.patch.object(
                    sc, "_build_bond_green",
                    side_effect=AssertionError("Green built before the check")):
                sc.calc_eliashberg(inp)
        self.assertEqual(os.listdir(out), ["chi0q.npz"])

    def test_non_so_chi0q_file(self):
        """A REAL non-spin-orbital RPA chi0q over ND orbitals, carrying every
        tag the RPA writes (index_convention = "spin_block" included), has
        exactly the shape of a spin-orbital chi0q over NP physical orbitals;
        only the spin_orbital = True marker, which ONLY a spin-orbital run
        writes, tells them apart: the non-spin-orbital file has no marker."""
        from hwave import qlms
        src = os.path.join(self.d, "nonso")
        os.makedirs(src)
        # the spin-orbital model files read as a plain ND-orbital model
        files = _write_model(src, so=True)
        rpa_in = _input(src, files, False, {}, out=os.path.join(self.d, "out"))
        rpa_in["mode"]["calc_scheme"] = "general"
        rpa_in["file"]["input"]["path_to_input"] = src
        del rpa_in["eliashberg"]
        rpa_in["log"] = {"print_level": 0}
        with _quiet():
            qlms.run(input_dict=rpa_in)
        path = os.path.join(self.d, "out", "chi0q.npz")
        nvol = int(np.prod(CELL))
        with np.load(path) as z:
            self.assertEqual(z["chi0q"].shape, (64, nvol) + (ND,) * 4)
            self.assertEqual(str(z["index_convention"]), "spin_block")
            self.assertEqual(str(z["momentum_convention"]), "e_plus_ikR")
            self.assertNotIn("spin_orbital", z.files)
        before = os.listdir(os.path.join(self.d, "out"))
        files = _write_model(self.d, so=True)
        inp = _input(self.d, files, True, {"chi0q_mode": "load",
                                           "solver_mode": "eigenvalue"})
        # refused at load as unmarked, before the Green function is built
        with self.assertRaisesRegex(
                ValueError, r"(?s)chi0q\.npz' has no spin_orbital marker.*"
                            r"enable_spin_orbital = true"):
            with _quiet(), mock.patch.object(
                    sc, "_build_bond_green",
                    side_effect=AssertionError("Green built before the check")):
                sc.calc_eliashberg(inp)
        self.assertEqual(sorted(os.listdir(os.path.join(self.d, "out"))),
                         sorted(before))
        # the same file with an explicit False marker: refused as such
        with np.load(path) as z:
            items = {k: z[k] for k in z.files}
        np.savez(path, spin_orbital=False, **items)
        with self.assertRaisesRegex(ValueError, r"records spin_orbital = False"):
            with _quiet():
                sc.calc_eliashberg(inp)

    def test_untagged_chi0q_file(self):
        out = os.path.join(self.d, "out")
        os.makedirs(out)
        nvol = int(np.prod(CELL))
        np.savez(os.path.join(out, "chi0q.npz"),
                 chi0q=np.zeros((64, nvol) + (ND,) * 4, complex))
        files = _write_model(self.d, so=True)
        inp = _input(self.d, files, True, {"chi0q_mode": "load",
                                           "solver_mode": "eigenvalue"})
        with self.assertRaisesRegex(ValueError, r"index_convention"):
            with _quiet():
                sc.calc_eliashberg(inp)

    def test_eigenvalue_method(self):
        self._refused(r"eigenvalue_method", {"eigenvalue_method": "shift-invert-gmres"})

    def test_explicit_pairing_type(self):
        for key in ("pairing_type", "Pairing_Type"):
            with self.subTest(key=key):
                self._refused(r"pairing_type", {key: "singlet"})

    def test_flex_spin_charge_file_keys(self):
        for key in ("flex_chi_s", "Flex_Chi_C"):
            with self.subTest(key=key):
                self._refused(key.lower(), {"chi0q_mode": "flex", key: "chiq_s.npz"})

    def test_flex_mode_without_archive(self):
        self._refused(r"chiq\.npz' not found", {"chi0q_mode": "flex"})

    def test_spinful_vertex_exchange_false(self):
        self._refused(r"spinful_vertex_exchange", spinful_vertex_exchange=False)

    def test_odd_spin_orbital_count(self):
        files = _write_model(self.d, so=True, n_geom=3)
        self._refused(r"even Geometry norb", so_files=files)

    def test_interaction_in_spin_orbital_indices(self):
        files = _write_model(self.d, so=True, types=("CoulombIntra",),
                             extra_rows={"CoulombIntra": [(0, 3, 3, 1.0)]})
        self._refused(r"PHYSICAL orbital indices", so_files=files)

    def test_interaction_orbital_index_zero(self):
        """A 1-based orbital index 0 (internally -1) would wrap to the last
        orbital in the vertex; it is refused in the early validation, before
        any Green, susceptibility or vertex work."""
        boom = AssertionError("heavy work before the orbital-index check")
        for nm, rows in (("CoulombIntra", [(0, -1, -1, 1.0)]),
                         ("CoulombInter", [(0, -1, 1, 0.3), (0, 1, -1, 0.3)]),
                         ("Hund", [(0, 0, -1, 0.1), (0, -1, 0, 0.1)])):
            with self.subTest(file=nm):
                files = _write_model(self.d, so=True, extra_rows={nm: rows})
                inp = _input(self.d, files, True, {"chi0q_mode": "calc",
                                                   "solver_mode": "eigenvalue"})
                with self.assertRaisesRegex(
                        ValueError,
                        r"the {} file '{}' has a row on orbital 0, outside the "
                        r"physical orbitals 1\.\.{}".format(
                            nm, files[nm].replace(".", r"\."), NP)), \
                        _quiet(), \
                        mock.patch.object(sc, "_build_bond_green", side_effect=boom), \
                        mock.patch.object(sc, "_calc_chi0q_internal", side_effect=boom), \
                        mock.patch.object(sc, "_load_chi0q", side_effect=boom), \
                        mock.patch.object(sc, "_load_flex_so", side_effect=boom), \
                        mock.patch.object(sc, "_so_static_vertex", side_effect=boom), \
                        mock.patch.object(sc, "_so_rpa_vertex", side_effect=boom):
                    sc.calc_eliashberg(inp)
                out = inp["file"]["output"]["path_to_output"]
                self.assertFalse(os.path.exists(out) and os.listdir(out))

    def test_extern_file(self):
        files = _write_model(self.d, so=True)
        _wfile(os.path.join(self.d, "extern.dat"), "Extern", NP,
               [(0, 0, 0, a + 1, a + 1, 0.1, 0.0) for a in range(NP)])
        for key in ("Extern", "EXTERN"):
            with self.subTest(key=key):
                self._refused(r"Extern", so_files=dict(files, **{key: "extern.dat"}))

    def _mean_field_refused(self, name, spellings):
        files = _write_model(self.d, so=True)
        with open(os.path.join(self.d, "mf.npz"), "wb") as f:
            f.write(b"placeholder")
        for key in spellings:
            with self.subTest(key=key):
                inp = _input(self.d, files, True,
                             {"chi0q_mode": "calc", "solver_mode": "eigenvalue"})
                inp["file"]["input"][key] = "mf.npz"
                boom = mock.Mock(side_effect=AssertionError("heavy work reached"))
                read = mock.Mock(side_effect=AssertionError("input files read"))
                with mock.patch.object(sc, "_so_rpa_vertex", boom), \
                        mock.patch("hwave.solver.rpa.RPA.read_init", boom), \
                        mock.patch.object(sc, "_read_interaction_files", read):
                    with self.assertRaisesRegex(ValueError, name + r".*#273"):
                        with _quiet():
                            sc.calc_eliashberg(inp)
                boom.assert_not_called()
                read.assert_not_called()
                out = inp["file"]["output"]["path_to_output"]
                self.assertFalse(os.path.exists(out) and os.listdir(out))

    def test_trans_mod_file(self):
        self._mean_field_refused("trans_mod", ("trans_mod", "Trans_Mod"))

    def test_green_init_file(self):
        self._mean_field_refused("green_init", ("green_init", "GREEN_INIT"))

    def test_init_gap_without_antisymmetric_component(self):
        """An odd seed on a single k point vanishes (sin k = 0): refused
        before the kernel is built; an even seed on the same grid runs. (The
        model keeps only the on-site transfer, which a 1 x 1 x 1 cell can
        hold.)"""
        files = _write_model(self.d, so=True)
        _wfile(os.path.join(self.d, "transfer.dat"), "Transfer", ND,
               [(0, 0, 0, i + 1, i + 1, 0.3 * (i // 2), 0.0)
                for i in range(ND)])
        boom = mock.Mock(side_effect=AssertionError("kernel built"))
        for seed, refused in (("p_x", True), ("s", False)):
            out = "out_" + seed
            inp = _input(self.d, files, True,
                         {"chi0q_mode": "calc", "solver_mode": "eigenvalue",
                          "num_eigenvalues": 2, "init_gap": seed},
                         CellShape=[1, 1, 1], out=out)
            with self.subTest(init_gap=seed), _quiet():
                if refused:
                    with mock.patch.object(sc, "_make_kernel_operator", boom), \
                            self.assertRaisesRegex(
                                ValueError, r"init_gap = 'p_x' has no fermion-"
                                            r"antisymmetric component"):
                        sc.calc_eliashberg(inp)
                    path = os.path.join(self.d, out)
                    self.assertFalse(os.path.exists(path) and os.listdir(path))
                else:
                    # control: the same grid with an even seed runs
                    sc.calc_eliashberg(inp)
        boom.assert_not_called()

    def test_load_layout_error_is_rewrapped(self):
        """A chi0q file ``_load_chi0q`` itself rejects (no supported layout)
        is refused in load mode with the spin-orbital requirement appended."""
        out = os.path.join(self.d, "out")
        os.makedirs(out)
        nvol = int(np.prod(CELL))
        np.savez(os.path.join(out, "chi0q.npz"),
                 chi0q=np.zeros((64, nvol) + (3,) * 4, complex),
                 index_convention="spin_block", momentum_convention="e_plus_ikR",
                 spin_orbital=True)
        files = _write_model(self.d, so=True)
        inp = _input(self.d, files, True, {"chi0q_mode": "load",
                                           "solver_mode": "eigenvalue"})
        with self.assertRaisesRegex(
                ValueError, r"(?s)matches no supported layout.*\(\[mode\] "
                            r"enable_spin_orbital = true with chi0q_mode = "
                            r"\"load\" needs the general \(4-index\) spinful "
                            r"chi0q of a spin-orbital RPA run over nd = 4"):
            with _quiet(), mock.patch.object(
                    sc, "_build_bond_green",
                    side_effect=AssertionError("Green built before the check")):
                sc.calc_eliashberg(inp)
        self.assertEqual(os.listdir(out), ["chi0q.npz"])

    def test_dynamic_entry_still_refuses(self):
        from hwave.solver import eliashberg_dynamic
        with self.assertRaisesRegex(ValueError, "enable_spin_orbital"):
            eliashberg_dynamic.solve_dynamic({"mode": {"enable_spin_orbital": True}})


class TestStaticVertexLayouts(_TmpDir):

    def test_raw_and_reference_chi0q_layouts_agree(self):
        """``_so_static_vertex`` reads the same static bubble from the raw
        ``(nmat, nvol, nd^4)`` and the reference ``(nd^4, Nx, Ny, Nz, nmat)``
        layouts."""
        files = _write_model(self.d, so=True, soc=0.3)
        inp = _input(self.d, files, True, {})
        Nx, Ny, Nz = CELL
        nvol, nmat = Nx * Ny * Nz, 6
        rng = np.random.default_rng(7)
        raw = 0.05 * (rng.standard_normal((nmat, nvol) + (ND,) * 4)
                      + 1j * rng.standard_normal((nmat, nvol) + (ND,) * 4))
        ref = raw.reshape((nmat, Nx, Ny, Nz) + (ND,) * 4).transpose(
            4, 5, 6, 7, 1, 2, 3, 0)
        for si in (None, 1):
            with self.subTest(static_index=si), _quiet():
                v_raw = sc._so_static_vertex(inp, raw, NP, Nx, Ny, Nz, si)
                v_ref = sc._so_static_vertex(inp, ref, NP, Nx, Ny, Nz, si)
                self.assertGreater(np.abs(v_raw).max(), 0.0)
                np.testing.assert_array_equal(v_raw, v_ref)
        with _quiet():
            v0 = sc._so_static_vertex(inp, raw, NP, Nx, Ny, Nz, None)
            v1 = sc._so_static_vertex(inp, raw, NP, Nx, Ny, Nz, 1)
        self.assertGreater(np.abs(v0 - v1).max(), 1e-6)


class TestFlexVertexCheck(unittest.TestCase):
    """``_check_so_flex_vertex`` measures the backward error PER q block:
    the same mismatch must be refused at every lattice size."""

    M = 16

    def _block(self):
        rng = np.random.default_rng(3)
        m = self.M
        c0 = rng.standard_normal((m, m)) + 1j * rng.standard_normal((m, m))
        c0 = 0.1 * (c0 + c0.conj().T)
        gam = rng.standard_normal((m, m))
        gam = 0.5 * (gam + gam.T)
        chi = np.linalg.solve(np.eye(m) + c0 @ gam, c0)
        return gam, c0, chi, rng.standard_normal((m, m))

    @staticmethod
    def _global_error(gam, chi, chi0):
        """The global-norm measure (for the non-vacuity check only)."""
        c0g = chi0 @ gam
        return (np.linalg.norm(chi + c0g @ chi - chi0)
                / (np.linalg.norm(chi0) + np.linalg.norm(c0g) * np.linalg.norm(chi)))

    def test_error_is_independent_of_nvol(self):
        gam, c0, chi, d = self._block()
        eps = 3.0e-8
        errs = []
        for nvol in (1, 4, 16, 64, 256):
            C0 = np.broadcast_to(c0, (nvol,) + c0.shape).copy()
            C = np.broadcast_to(chi + eps * d, (nvol,) + c0.shape).copy()
            with self.subTest(nvol=nvol):
                with _quiet(), self.assertRaisesRegex(
                        ValueError, r"do not reproduce the vertex .*\(backward "
                                    r"error 3\.00\de-08 >"):
                    sc._check_so_flex_vertex(gam, C, C0, "chiq.npz")
                # the error itself, through the tolerance-free path
                with _quiet(), mock.patch.object(sc, "_SO_FLEX_GAMMA_TOL", np.inf):
                    errs.append(sc._check_so_flex_vertex(gam, C, C0, "chiq.npz"))
                if nvol >= 16:
                    # non-vacuity: a global norm would let this pass
                    self.assertLess(self._global_error(gam, C, C0),
                                    sc._SO_FLEX_GAMMA_TOL)
        self.assertGreater(errs[0], sc._SO_FLEX_GAMMA_TOL)
        np.testing.assert_allclose(errs, errs[0], rtol=1e-10)

    def test_one_bad_block_among_many(self):
        gam, c0, chi, d = self._block()
        nvol = 256
        C0 = np.broadcast_to(c0, (nvol,) + c0.shape).copy()
        C = np.broadcast_to(chi, (nvol,) + c0.shape).copy()
        C[nvol // 3] += 3.0e-8 * d
        with _quiet(), self.assertRaisesRegex(ValueError, r"do not reproduce"):
            sc._check_so_flex_vertex(gam, C, C0, "chiq.npz")
        with _quiet():
            C[nvol // 3] = chi
            self.assertLess(sc._check_so_flex_vertex(gam, C, C0, "chiq.npz"),
                            sc._SO_FLEX_GAMMA_TOL)

    def test_zero_scale_block(self):
        """A block with ``chi0_q = 0``: passes when ``chi_q = 0`` too, and is
        refused (error infinity, no division by zero) when ``chi_q != 0``."""
        gam, c0, chi, d = self._block()
        C0 = np.stack([c0, np.zeros_like(c0)])
        C = np.stack([chi, np.zeros_like(chi)])
        with _quiet(), np.errstate(all="raise"):
            self.assertLess(sc._check_so_flex_vertex(gam, C, C0, "chiq.npz"),
                            sc._SO_FLEX_GAMMA_TOL)
            C[1] = 1e-3 * d
            with self.assertRaisesRegex(ValueError, r"backward error inf"):
                sc._check_so_flex_vertex(gam, C, C0, "chiq.npz")

    def test_nan_is_refused(self):
        gam, c0, chi, d = self._block()
        C0 = np.stack([c0, c0])
        C = np.stack([chi, chi])
        C[1, 0, 0] = np.nan
        with _quiet(), self.assertRaisesRegex(ValueError, r"backward error nan"):
            sc._check_so_flex_vertex(gam, C, C0, "chiq.npz")


class TestAntisymmetricOperator(unittest.TestCase):
    """``_so_antisymmetric_operator(A).matmat(B)`` equals the column-stacked
    ``matvec`` results, and every column is a fermion-antisymmetric gap
    ``Delta_AB(k) = -Delta_BA(-k)``, on even and odd k-grid dimensions."""

    @staticmethod
    def _partner(g):
        """``Delta_BA(-k)`` by explicit index arithmetic (k_i -> (-k_i) mod n)."""
        nd, _, nx, ny, nz = g.shape
        out = np.empty_like(g)
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    out[:, :, ix, iy, iz] = g[:, :, (-ix) % nx, (-iy) % ny,
                                              (-iz) % nz].T
        return out

    def test_matmat_equals_matvec_and_is_antisymmetric(self):
        from scipy.sparse.linalg import aslinearoperator
        nd = 4
        for grid in ((4, 2, 1), (3, 5, 1), (2, 3, 3)):
            with self.subTest(grid=grid):
                shape = (nd, nd) + grid
                n = int(np.prod(shape))
                rng = np.random.default_rng(sum(grid))
                M = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
                PK = sc._so_antisymmetric_operator(aslinearoperator(M), n, shape)
                B = rng.standard_normal((n, 5)) + 1j * rng.standard_normal((n, 5))
                got = PK.matmat(B)
                self.assertEqual(got.shape, (n, 5))
                cols = np.stack([PK.matvec(B[:, j]) for j in range(5)], axis=1)
                np.testing.assert_allclose(got, cols, rtol=0, atol=1e-12)
                for j in range(5):
                    g = got[:, j].reshape(shape)
                    self.assertGreater(np.linalg.norm(g), 1.0)
                    np.testing.assert_allclose(g, -self._partner(g),
                                               rtol=0, atol=1e-10)
                # non-vacuity: the raw kernel output is not antisymmetric
                raw = (M @ B[:, 0]).reshape(shape)
                self.assertGreater(np.linalg.norm(raw + self._partner(raw)), 1.0)


class TestRepulsiveSpectralShift(_TmpDir):
    """A repulsive-dominated spin-orbital kernel under ``spectral_shift``.

    ``P K`` maps the Pauli-forbidden (symmetric) sector -- about half the
    space -- to eigenvalue 0. Under the plain shift ``P K + sigma I`` that
    sector sits at ``+sigma``, above every physical eigenvalue when all of
    them are negative, and the largest-real selection returned it. The model
    is the complex-SOC model with the kernel lowered by one (``K - 1``, so
    ``P (K - 1) = P K - P``: every antisymmetric eigenvalue moves down by 1
    and the leading one is negative); the operator is still built by
    ``_so_antisymmetric_operator``. The solver must return the leading
    antisymmetric eigenvalues, equal to the dense diagonalization of the
    kernel restricted to the antisymmetric subspace."""

    def _lowered_run(self, eli):
        from scipy.sparse.linalg import LinearOperator
        files = _write_model(self.d, so=True, soc=0.4)
        real = sc._so_antisymmetric_operator
        built = {}

        def lowered(A, n, shape):
            Am = LinearOperator((n, n), matvec=lambda v: A.matvec(v) - v,
                                dtype=complex)
            built["op"] = real(Am, n, shape)
            built["shape"] = shape
            return built["op"]

        with mock.patch.object(sc, "_so_antisymmetric_operator",
                               side_effect=lowered):
            seen = _run(_input(self.d, files, True, eli))
        return seen, built

    @staticmethod
    def _dense_antisymmetric(built):
        op, shape = built["op"], built["shape"]
        n = int(np.prod(shape))
        eye = np.eye(n, dtype=complex)
        M = np.stack([op.matvec(eye[:, j]) for j in range(n)], axis=1)
        P = np.stack([sc._project_gap_parity(eye[:, j].reshape(shape),
                                             "triplet").ravel()
                      for j in range(n)], axis=1)
        w, U = np.linalg.eigh(P)
        Q = U[:, w > 0.5]
        vals = np.linalg.eigvals(Q.conj().T @ M @ Q)
        return vals[np.argsort(-vals.real)], Q.shape[1], n

    def test_shifted_arnoldi_returns_the_antisymmetric_leading(self):
        seen, built = self._lowered_run(
            {"chi0q_mode": "calc", "solver_mode": "eigenvalue",
             "num_eigenvalues": 4, "spectral_shift": 3.0})
        dense, n_anti, n = self._dense_antisymmetric(built)
        self.assertLess(n_anti, n)
        # repulsive-dominated: every physical eigenvalue is negative, so the
        # forbidden sector (0 of P K) would outrank all of them
        self.assertLess(dense[0].real, -0.1)
        self.assertTrue(np.all(seen["match"]))
        np.testing.assert_allclose(seen["evs"], dense[:4], rtol=0, atol=1e-9)
        gap = seen["gap"]
        np.testing.assert_allclose(gap, -_reverse(gap), rtol=0,
                                   atol=1e-12 * np.abs(gap).max())

    def test_shifted_auto_and_iteration(self):
        """``spectral_shift = "auto"`` (Arnoldi) and the shifted power
        iteration reach the same leading antisymmetric eigenvalue."""
        seen, built = self._lowered_run(
            {"chi0q_mode": "calc", "solver_mode": "both",
             "num_eigenvalues": 4, "spectral_shift": "auto",
             "convergence_tol": 1e-10, "max_iter": 20000})
        dense, _, _ = self._dense_antisymmetric(built)
        self.assertLess(dense[0].real, -0.1)
        self.assertTrue(np.all(seen["match"]))
        np.testing.assert_allclose(seen["evs"], dense[:4], rtol=0, atol=1e-9)
        self.assertAlmostEqual(seen["ev"], dense[0].real,
                               delta=1e-7 * abs(dense[0]))

    def test_driver_confines_unprojected_seeds(self):
        """``_solve_leading`` on a sector-confined operator: a seed (or a
        Krylov start) with a forbidden component converges to the physical
        leading eigenpair -- in the shifted power iteration the forbidden
        sector is annihilated (``P K + sigma P``), in the shifted Arnoldi it
        lies below the spectrum."""
        from scipy.sparse.linalg import aslinearoperator
        shape = (4, 4, 2, 2, 1)
        n = int(np.prod(shape))
        rng = np.random.default_rng(7)
        B = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
        M = -(B @ B.conj().T) / n - 0.5 * np.eye(n)
        PK = sc._so_antisymmetric_operator(aslinearoperator(M), n, shape)
        eye = np.eye(n)
        P = np.stack([PK.sector_projector(eye[:, j]) for j in range(n)], axis=1)
        w, U = np.linalg.eigh(P)
        Q = U[:, w > 0.5]
        dense = np.sort(np.linalg.eigvalsh(Q.conj().T @ M @ Q))[::-1]
        self.assertLess(dense[0], -0.4)
        seed = rng.standard_normal(n) + 1j * rng.standard_normal(n)
        with _quiet():
            lam, vec, info = sc._solve_leading(
                lambda: (PK, n), n, "iteration", init_vec=seed,
                max_iter=5000, convergence_tol=1e-11, spectral_shift=10.0)
        self.assertTrue(info["converged"])
        self.assertAlmostEqual(lam, dense[0], delta=1e-8)
        g = vec.reshape(shape)
        np.testing.assert_allclose(g, -_reverse(g), rtol=0, atol=1e-12)
        with _quiet():
            lam, vec, info = sc._solve_leading(
                lambda: (PK, n), n, "arnoldi", num_eigenvalues=3,
                spectral_shift=10.0)
        np.testing.assert_allclose(info["eigenvalues"], dense[:3], rtol=0,
                                   atol=1e-9)

    def test_forbidden_leading_is_refused(self):
        """Defensive check: an operator whose computed eigenpairs are all
        Pauli-forbidden (the plain shift, without the sector deflation) is
        refused instead of written under the antisymmetric header."""
        from scipy.sparse.linalg import LinearOperator
        files = _write_model(self.d, so=True, soc=0.4)
        real = sc._so_antisymmetric_operator

        def plain(A, n, shape):
            PK = real(A, n, shape)
            P = PK.sector_projector
            return LinearOperator((n, n), matvec=lambda v: PK.matvec(v) - P(v),
                                  dtype=complex)

        inp = _input(self.d, files, True,
                     {"chi0q_mode": "calc", "solver_mode": "eigenvalue",
                      "num_eigenvalues": 4, "spectral_shift": 3.0})
        with _quiet(), mock.patch.object(sc, "_so_antisymmetric_operator",
                                         side_effect=plain):
            with self.assertRaisesRegex(ValueError, "fermion-antisymmetric"):
                sc.calc_eliashberg(inp)
        self.assertFalse(os.path.exists(os.path.join(self.d, "out", "gap.dat")))


class _Stop(Exception):
    """Ends a run once the array under test has been captured."""


class TestTransferUnderSOC(_TmpDir):

    def test_epsilon_k_equals_spin_orbital_rpa(self):
        """With spin-orbit hopping (soc != 0), the epsilon(k) hwave_sc
        diagonalizes for its Green function (the interleaved Transfer file
        reordered to spin-block) equals the H0(k) the spin-orbital RPA
        assembles, to 1e-12. The complex-conjugate, the k -> -k and the
        un-reordered (interleaved) readings all miss by O(soc), so the
        comparison is not vacuous."""
        from hwave import qlms
        from hwave.solver import rpa as rpa_mod
        from hwave.solver import spin_orbital
        files = _write_model(self.d, so=True, soc=0.35)
        Nx, Ny, Nz = CELL
        nvol = Nx * Ny * Nz
        seen = {}

        def spy_eig(eps, *a, **kw):
            seen["sc"] = np.array(eps, copy=True)
            raise _Stop()

        inp = _input(self.d, files, True, {"chi0q_mode": "calc",
                                           "solver_mode": "eigenvalue"})
        with _quiet(), mock.patch.object(sc, "_calc_eigenvalues",
                                         side_effect=spy_eig):
            with self.assertRaises(_Stop):
                sc.calc_eliashberg(inp)

        real_eps = rpa_mod.RPA._calc_epsilon_k

        def spy_rpa(self_, green_info):
            real_eps(self_, green_info)
            seen["rpa"] = np.array(self_.H0_k, copy=True)
            seen["spin_mode"] = self_.spin_mode
            raise _Stop()

        rpa_in = _input(self.d, files, True, {}, out="out_rpa")
        rpa_in["mode"]["calc_scheme"] = "general"
        rpa_in["file"]["input"]["path_to_input"] = self.d
        del rpa_in["eliashberg"]
        rpa_in["log"] = {"print_level": 0}
        with _quiet(), mock.patch.object(rpa_mod.RPA, "_calc_epsilon_k",
                                         autospec=True, side_effect=spy_rpa):
            with self.assertRaises(_Stop):
                qlms.run(input_dict=rpa_in)

        self.assertEqual(seen["spin_mode"], "spinful")
        H = seen["rpa"].reshape(nvol, ND, ND)
        eps = seen["sc"]
        self.assertEqual(eps.shape, (ND, ND, Nx, Ny, Nz))
        E = np.moveaxis(eps.reshape(ND, ND, nvol), -1, 0)
        self.assertGreater(np.abs(H).max(), 1.0)
        # spin-flip blocks are present, so the reorder is exercised
        self.assertGreater(np.abs(H[:, :NP, NP:]).max(), 0.1)
        np.testing.assert_allclose(E, H, rtol=0, atol=1e-12)
        # the readings a sign or ordering slip would produce all miss
        ix, iy, iz = np.meshgrid(np.arange(Nx), np.arange(Ny), np.arange(Nz),
                                 indexing="ij")
        minus_k = (((-ix) % Nx) * Ny * Nz + ((-iy) % Ny) * Nz
                   + (-iz) % Nz).ravel()
        interleaved = spin_orbital.spinblock_to_interleaved(H, NP)
        for name, other in (("conjugate", H.conj()),
                            ("k -> -k", H[minus_k]),
                            ("interleaved", interleaved)):
            with self.subTest(reading=name):
                self.assertGreater(np.abs(E - other).max(), 0.1)


def _run_flex(d, files, out, so=True, scheme="general", nmat=32,
              cell=None, T=0.25):
    """A FLEX run of the model in ``d`` through the public entry, writing
    ``chiq``, ``chi0q`` and ``green`` into ``d/out``."""
    from hwave import qlms
    inp = _input(d, files, so, {}, T=T, nmat=nmat, out=out,
                 **({"CellShape": list(cell)} if cell else {}))
    inp["mode"]["mode"] = "FLEX"
    inp["mode"]["calc_scheme"] = scheme
    inp["mode"]["param"].update({
        "IterationMax": 200, "Mix": 0.5, "EPS": 10, "mixing_scheme": "linear",
        "flex_hartree_fock": scheme == "general"})
    inp["file"]["input"]["path_to_input"] = d
    inp["file"]["output"].update({"chiq": "chiq", "green": "green",
                                  "chi0q": "chi0q"})
    del inp["eliashberg"]
    inp["log"] = {"print_level": 0}
    with _quiet():
        qlms.run(input_dict=inp)
    return os.path.join(d, out)


def _flex_input(d, files, flex_dir, eli=None, out="out_sc", **param):
    inp = _input(d, files, True,
                 dict({"chi0q_mode": "flex", "solver_mode": "eigenvalue",
                       "num_eigenvalues": 8}, **(eli or {})),
                 nmat=32, out=out, **param)
    inp["file"]["input"]["path_to_flex_output"] = flex_dir
    return inp


def _write_one_orbital(d, soc=0.3, U=2.0):
    """One physical orbital (nd = 2) on the square lattice: nearest-neighbour
    hopping plus the Rashba-like ``soc (sin ky sigma_x - sin kx sigma_y)``;
    CoulombIntra ``U``. The transfer uses the interleaved index ``s + 1``."""
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1 0 0\n0 1 0\n0 0 1\n2\n0 0 0\n0 0 0\n")
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    rows = []
    for R, m in (((1, 0), -np.eye(2) + 1j * soc * sy),
                 ((-1, 0), -np.eye(2) - 1j * soc * sy),
                 ((0, 1), -np.eye(2) - 1j * soc * sx),
                 ((0, -1), -np.eye(2) + 1j * soc * sx)):
        for s_ in range(2):
            for sp in range(2):
                if m[s_, sp] != 0:
                    rows.append((R[0], R[1], 0, s_ + 1, sp + 1,
                                 m[s_, sp].real, m[s_, sp].imag))
    _wfile(os.path.join(d, "transfer.dat"), "Transfer", 2, rows)
    _wfile(os.path.join(d, "coulombintra.dat"), "CoulombIntra", 1,
           [(0, 0, 0, 1, 1, U, 0.0)])
    return {"CoulombIntra": "coulombintra.dat"}


def _resave(src, dst, drop=(), **replace):
    """Copy the npz ``src`` to ``dst`` without the keys ``drop``."""
    with np.load(src) as z:
        items = {k: z[k] for k in z.files if k not in drop}
    items.update(replace)
    np.savez(dst, **items)


class TestFlexMode(unittest.TestCase):
    """``chi0q_mode = "flex"`` reads the stored ``chiq`` (the single
    antisymmetrized spinful channel) and the dressed ``green`` of a
    spin-orbital FLEX run, both tagged ``index_convention = "spin_block"``;
    the archive is identified by its tags, never by its shape."""

    @classmethod
    def setUpClass(cls):
        cls.d = tempfile.mkdtemp(prefix="sc_so_flex_")
        cls.files = _write_model(cls.d, so=True, soc=0.3)
        cls.flex = _run_flex(cls.d, cls.files, "flex")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.d, True)

    def _out(self, name):
        out = os.path.join(self.d, name)
        self.addCleanup(shutil.rmtree, out, True)
        return out

    def _refused(self, pattern, inp):
        with self.assertRaisesRegex(ValueError, pattern):
            with _quiet():
                sc.calc_eliashberg(inp)
        out = inp["file"]["output"]["path_to_output"]
        self.assertFalse(os.path.exists(out) and os.listdir(out))

    def _dense_spectrum(self, inp, chi_static, green):
        """Eigenvalues of ``P K`` built independently of the flex loader from
        a static chi ``(nvol, nd^2, nd^2)`` and a Green function in the
        stored ``(nmat, nvol, nd, nd)`` layout."""
        Nx, Ny, Nz = CELL
        nmat = green.shape[0]
        g = green.reshape(nmat, Nx, Ny, Nz, ND, ND).transpose(4, 5, 1, 2, 3, 0)
        with _quiet():
            gam = sc._so_rpa_vertex(inp)
        Vs = sc._so_static_vertex_from_chi(gam, chi_static, ND, (Nx, Ny, Nz))
        G2 = sc._calc_g2(np.ascontiguousarray(g), 1.0 / 0.25, tail=True)
        n = ND * ND * Nx * Ny * Nz
        A, _ = sc._make_kernel_operator(Vs, G2, ND, Nx, Ny, Nz)
        PK = sc._so_antisymmetric_operator(A, n, (ND, ND, Nx, Ny, Nz))
        return np.linalg.eigvals(PK.matmat(np.eye(n, dtype=complex)))

    def test_uses_stored_chi_and_green(self):
        """Input agreement: the eigenvalues equal those of ``P K`` built from the stored
        chi at the zero bosonic frequency and the stored green."""
        inp = _flex_input(self.d, self.files, self.flex, out=self._out("out_g4"))
        seen = _run(inp)
        with np.load(os.path.join(self.flex, "chiq.npz")) as z:
            chiq = z["chiq"]
        with np.load(os.path.join(self.flex, "chi0q.npz")) as z:
            chi0q = z["chi0q"]
        with np.load(os.path.join(self.flex, "green.npz")) as z:
            green = z["green"][0]
        nmat, nvol = chiq.shape[:2]
        si = nmat // 2
        chi = chiq[si].reshape(nvol, ND * ND, ND * ND)
        chi0 = chi0q[si].reshape(nvol, ND * ND, ND * ND)
        # non-vacuity: the dressed chi is far from the bare bubble
        self.assertGreater(np.linalg.norm(chi - chi0) / np.linalg.norm(chi0), 0.01)
        dense = self._dense_spectrum(inp, chi, green)
        evs = np.asarray(seen["evs"])
        self.assertTrue(np.all(seen["match"]))
        for lam in evs:
            k = np.argmin(np.abs(dense - lam))
            self.assertLess(abs(dense[k] - lam), 1e-10 * abs(lam), msg=str(lam))
        lead = dense[np.argmax(dense.real)]
        self.assertLess(abs(evs[0] - lead), 1e-10 * abs(lead))
        # the leading eigenvalue depends on WHICH chi is read: the bare
        # bubble moves it far beyond the bound above
        dense0 = self._dense_spectrum(inp, chi0, green)
        self.assertGreater(abs(dense0[np.argmax(dense0.real)] - lead), 1e-4 * abs(lead))

    def test_off_centre_static_frequency(self):
        """The zero bosonic frequency is located through the archive's
        ``freq_index`` / ``nmat``, not the centre of the stored axis: chiq and
        chi0q restricted to the original indices 14..21 of nmat = 32 (zero
        frequency, original index 16, at position 2; the centre, position 4,
        is a finite frequency) give the eigenvalues of the full archives --
        also when the configured Nmat (64) differs from the file's nmat."""
        bad = os.path.join(self.d, "flex_off_centre")
        self.addCleanup(shutil.rmtree, bad, True)
        shutil.copytree(self.flex, bad)
        sel = np.arange(14, 22)
        for name, key in (("chiq.npz", "chiq"), ("chi0q.npz", "chi0q")):
            with np.load(os.path.join(self.flex, name)) as z:
                fi, arr = z["freq_index"], z[key]
                self.assertEqual(int(z["nmat"]), 32)
            np.testing.assert_array_equal(fi, np.arange(32))
            _resave(os.path.join(self.flex, name), os.path.join(bad, name),
                    **{key: arr[sel], "freq_index": fi[sel]})
            # non-vacuity: the centre of the restricted axis is another slice
            self.assertGreater(np.abs(arr[sel][len(sel) // 2] - arr[16]).max(),
                               1e-3 * np.abs(arr[16]).max())
        ref = _run(_flex_input(self.d, self.files, self.flex,
                               out=self._out("out_full")))
        for nmat in (32, 64):
            with self.subTest(config_nmat=nmat):
                got = _run(_flex_input(self.d, self.files, bad, Nmat=nmat,
                                       out=self._out("out_off_%d" % nmat)))
                np.testing.assert_allclose(got["evs"], ref["evs"], rtol=1e-12,
                                           atol=0)

    def test_su2_limit_against_non_so_flex(self):
        """Independent oracle at the FLEX level: soc = 0 spin-orbital FLEX
        followed by the spin-orbital hwave_sc reproduces the non-SO general
        FLEX followed by the non-SO singlet / triplet hwave_sc -- every
        allowed singlet eigenvalue once and every allowed triplet eigenvalue
        three times (the multiplicities of the static SU(2) test). The two
        FLEX runs reach the same fixed point (their dressed Green functions
        agree to round-off), checked first."""
        root = tempfile.mkdtemp(prefix="sc_so_flex_su2_")
        self.addCleanup(shutil.rmtree, root, True)
        d_so, d_n = os.path.join(root, "so"), os.path.join(root, "nonso")
        os.makedirs(d_so)
        os.makedirs(d_n)
        f_so = _write_model(d_so, so=True)
        f_n = _write_model(d_n, so=False)
        flex_so = _run_flex(d_so, f_so, "flex")
        flex_n = _run_flex(d_n, f_n, "flex", so=False)
        with np.load(os.path.join(flex_so, "green.npz")) as a, \
                np.load(os.path.join(flex_n, "green.npz")) as b:
            gs, gn = a["green"][0], b["green"][0]
        gs = gs.reshape(gs.shape[:2] + (2, NP, 2, NP))
        for s_ in range(2):
            np.testing.assert_allclose(gs[:, :, s_, :, s_, :], gn, rtol=0,
                                       atol=1e-12)
        self.assertEqual(np.abs(gs[:, :, 0, :, 1, :]).max(), 0.0)
        res = {}
        for ch in ("singlet", "triplet"):
            inp = _input(d_n, f_n, False,
                         {"chi0q_mode": "flex", "solver_mode": "eigenvalue",
                          "num_eigenvalues": 62, "pairing_type": ch},
                         nmat=32, out="sc_" + ch)
            inp["file"]["input"]["path_to_flex_output"] = flex_n
            res[ch] = _run(inp)
        so = _run(_flex_input(d_so, f_so, flex_so, {"num_eigenvalues": 254}))
        _assert_su2_multiplicities(self, res, so, forbidden_floor=0.25)

    def test_archive_content_refusals(self):
        """``_load_flex_so`` refuses a chiq archive whose orbital layout is
        not ``acbd``, that has no ``chiq`` array, whose ``index_convention``
        has another VALUE, or a non-finite chiq / chi0q static slice -- each
        through the public entry, with nothing written."""
        with np.load(os.path.join(self.flex, "chiq.npz")) as z:
            chiq = np.array(z["chiq"])
        with np.load(os.path.join(self.flex, "chi0q.npz")) as z:
            chi0q = np.array(z["chi0q"])
        si = chiq.shape[0] // 2
        nan_chi, nan_chi0 = chiq.copy(), chi0q.copy()
        nan_chi[si, 3, 0, 1, 0, 1] = np.nan
        nan_chi0[si, 5, 1, 0, 1, 0] = np.inf
        cases = (
            ("layout", "chiq.npz", (), {"chi_orbital_layout": "abcd"},
             r"chi_orbital_layout = 'abcd' but this reader only supports "
             r"'acbd'"),
            ("no_chiq", "chiq.npz", ("chiq",), {}, r"has no 'chiq' array"),
            ("index_value", "chiq.npz", (), {"index_convention": "interleaved"},
             r"index_convention = 'interleaved', expected 'spin_block'"),
            # a non-finite chiq is caught earlier, by the whole-array scan
            # of the momentum-convention check
            ("nan_chiq", "chiq.npz", (), {"chiq": nan_chi},
             r"chiq\.npz' contains non-finite values"),
            ("nan_chi0q", "chi0q.npz", (), {"chi0q": nan_chi0},
             r"FLEX chi0q in .*non-finite values at the zero bosonic"),
        )
        for tag, name, drop, replace, pattern in cases:
            bad = os.path.join(self.d, "flex_bad_" + tag)
            self.addCleanup(shutil.rmtree, bad, True)
            shutil.copytree(self.flex, bad)
            _resave(os.path.join(self.flex, name), os.path.join(bad, name),
                    drop, **replace)
            with self.subTest(case=tag), mock.patch.object(
                    sc, "_so_rpa_vertex",
                    side_effect=AssertionError("vertex built before the check")):
                self._refused(pattern, _flex_input(
                    self.d, self.files, bad, out=self._out("out_bad_" + tag)))

    def test_unconverged_flex_warns(self):
        """A chiq archive recording ``scf_converged = false`` is used, with a
        WARNING naming the archive."""
        bad = os.path.join(self.d, "flex_unconverged")
        self.addCleanup(shutil.rmtree, bad, True)
        shutil.copytree(self.flex, bad)
        _resave(os.path.join(self.flex, "chiq.npz"),
                os.path.join(bad, "chiq.npz"), scf_converged=False)
        inp = _flex_input(self.d, self.files, bad, out=self._out("out_unconv"))
        with self.assertLogs("hwave_sc", level="WARNING") as cm:
            sc.calc_eliashberg(inp)
        hits = [m for m in cm.output if "scf_converged = false" in m]
        self.assertEqual(len(hits), 1, cm.output)
        self.assertIn("chiq.npz", hits[0])
        # the converged archive does not warn about convergence
        with self.assertLogs("hwave_sc", level="WARNING") as cm:
            sc.calc_eliashberg(_flex_input(self.d, self.files, self.flex,
                                           out=self._out("out_conv")))
            logging.getLogger("hwave_sc").warning("sentinel")
        self.assertFalse([m for m in cm.output if "scf_converged" in m])

    def test_requires_chi0q_archive(self):
        """Without the run's chi0q the vertex cannot be checked: refused
        before the green or any vertex is built."""
        bad = os.path.join(self.d, "flex_no_chi0q")
        self.addCleanup(shutil.rmtree, bad, True)
        shutil.copytree(self.flex, bad)
        os.remove(os.path.join(bad, "chi0q.npz"))
        with mock.patch.object(
                sc, "_load_flex_green",
                side_effect=AssertionError("green read before the check")), \
                mock.patch.object(
                    sc, "_so_rpa_vertex",
                    side_effect=AssertionError("vertex built before the check")):
            self._refused(r"chi0q\.npz' not found.*\[file\.output\] chi0q.*"
                          r"\[eliashberg\] flex_chi0q",
                          _flex_input(self.d, self.files, bad,
                                      out=self._out("out_no_chi0q")))

    def test_file_keys_are_case_insensitive(self):
        """``flex_chiq`` / ``flex_chi0q`` / ``flex_green`` under mixed-case
        keys select renamed copies; the result equals the default names'."""
        ren = os.path.join(self.d, "flex_renamed")
        self.addCleanup(shutil.rmtree, ren, True)
        os.makedirs(ren)
        for src, dst in (("chiq.npz", "renamed.npz"), ("chi0q.npz", "bare.npz"),
                         ("green.npz", "g2.npz")):
            shutil.copy(os.path.join(self.flex, src), os.path.join(ren, dst))
        ref = _run(_flex_input(self.d, self.files, self.flex,
                               out=self._out("out_default")))
        keys = {"Flex_ChiQ": "renamed.npz", "Flex_Chi0Q": "bare.npz",
                "FLEX_GREEN": "g2.npz"}
        with self.subTest(keys=keys):
            got = _run(_flex_input(self.d, self.files, ren, keys,
                                   out=self._out("out_renamed")))
            np.testing.assert_array_equal(got["evs"], ref["evs"])
        # a lower-cased name the files do not have is not silently defaulted
        for k in keys:
            with self.subTest(missing=k):
                self._refused(r"not found", _flex_input(
                    self.d, self.files, ren, dict(keys, **{k: "absent.npz"}),
                    out=self._out("out_absent")))

    def test_vertex_mismatch_is_refused(self):
        """Interaction files that do not reproduce the FLEX run's vertex are
        refused through the stored chi0q."""
        d = tempfile.mkdtemp(prefix="sc_so_flex_u_")
        self.addCleanup(shutil.rmtree, d, True)
        files = _write_model(d, so=True, soc=0.3)
        with open(os.path.join(d, "coulombintra.dat")) as f:
            txt = f.read()
        with open(os.path.join(d, "coulombintra.dat"), "w") as f:
            f.write(txt.replace("1.600000000000", "1.700000000000"))
        self._refused(r"do not reproduce the vertex",
                      _flex_input(d, files, self.flex, out=self._out("out_mis")))

    def test_norb1_shape_ambiguity_uses_tags(self):
        """N = 1 (nd = 2, the shape of a two-orbital non-spin-orbital file)
        loads by its tags; a copy without a tag is refused."""
        d = tempfile.mkdtemp(prefix="sc_so_flex_n1_")
        self.addCleanup(shutil.rmtree, d, True)
        files = _write_one_orbital(d)
        flex = _run_flex(d, files, "flex")
        inp = _flex_input(d, files, flex, out="out_ok")
        seen = _run(inp)
        self.assertEqual(seen["gap"].shape, (2, 2) + tuple(CELL))
        self.assertTrue(np.isfinite(seen["evs"][0]))
        cases = (("chiq.npz", ("chi_convention",), r"no chi_convention tag"),
                 ("chiq.npz", ("index_convention",), r"no index_convention tag"),
                 ("green.npz", ("index_convention",), r"index_convention = None"))
        for name, drop, pattern in cases:
            with self.subTest(file=name, drop=drop):
                bad = os.path.join(d, "bad_" + "_".join(drop) + name[:-4])
                shutil.copytree(flex, bad)
                _resave(os.path.join(flex, name), os.path.join(bad, name), drop)
                self._refused(pattern, _flex_input(d, files, bad,
                                                   out="out_" + os.path.basename(bad)))

    def test_so_outputs_carry_the_marker(self):
        """The spin-orbital FLEX run writes ``spin_orbital = True`` (0-d
        bool) into its chi0q and chiq archives."""
        for name in ("chi0q.npz", "chiq.npz"):
            with self.subTest(file=name), \
                    np.load(os.path.join(self.flex, name)) as z:
                self.assertIn("spin_orbital", z.files)
                m = np.asarray(z["spin_orbital"])
                self.assertEqual((m.ndim, m.dtype), (0, np.bool_))
                self.assertTrue(bool(m))

    def _assert_chi0q_marker_refused(self, d, files, flex, label):
        """The chi0q companion of ``flex`` replaced by unmarked / False /
        malformed copies of the same shape: refused before the green or any
        vertex is built; the marked original gets past the check (the green
        mock fires), so the marker is what refuses."""
        green_boom = AssertionError("green read before the check")
        for tag, replace, pattern in (
                ("unmarked", None, r"chi0q\.npz' has no spin_orbital marker"),
                ("false", False, r"records spin_orbital = False: it was not "
                                 r"written by a spin-orbital FLEX run"),
                ("malformed", np.array([True]), r"malformed spin_orbital")):
            bad = os.path.join(d, "{}_{}".format(label, tag))
            self.addCleanup(shutil.rmtree, bad, True)
            shutil.copytree(flex, bad)
            src = os.path.join(flex, "chi0q.npz")
            dst = os.path.join(bad, "chi0q.npz")
            if replace is None:
                _resave(src, dst, ("spin_orbital",))
            else:
                _resave(src, dst, spin_orbital=replace)
            with np.load(src) as a, np.load(dst) as b:
                self.assertEqual(a["chi0q"].shape, b["chi0q"].shape)
            with self.subTest(case=label, marker=tag), mock.patch.object(
                    sc, "_load_flex_green", side_effect=green_boom), \
                    mock.patch.object(
                        sc, "_so_flex_static_slice",
                        side_effect=AssertionError("array sliced before the check")), \
                    mock.patch.object(
                        sc, "_so_rpa_vertex",
                        side_effect=AssertionError("vertex built before the check")):
                self._refused(pattern, _flex_input(
                    d, files, bad, out=self._out("out_{}_{}".format(label, tag))))
        with mock.patch.object(sc, "_load_flex_green", side_effect=green_boom), \
                self.assertRaisesRegex(AssertionError, "green read"), _quiet():
            sc.calc_eliashberg(_flex_input(d, files, flex,
                                           out=self._out("out_" + label + "_ok")))

    def test_chiq_metadata_is_checked_before_the_chi0q_marker(self):
        """A chiq with a malformed momentum_convention AND a chi0q with an
        invalid marker: the chiq metadata error is the one raised."""
        bad = os.path.join(self.d, "flex_bad_momentum")
        self.addCleanup(shutil.rmtree, bad, True)
        shutil.copytree(self.flex, bad)
        _resave(os.path.join(self.flex, "chiq.npz"), os.path.join(bad, "chiq.npz"),
                momentum_convention=np.array(["e_plus_ikR", "e_plus_ikR"]))
        _resave(os.path.join(self.flex, "chi0q.npz"),
                os.path.join(bad, "chi0q.npz"), spin_orbital=False)
        with mock.patch.object(
                sc, "_so_flex_static_slice",
                side_effect=AssertionError("array sliced before the check")):
            self._refused(r"chiq\.npz': momentum_convention must be a single "
                          r"value", _flex_input(self.d, self.files, bad,
                                                out=self._out("out_bad_mom")))

    def test_unmarked_chi0q_companion_is_refused(self):
        """A same-shaped chi0q companion without ``spin_orbital = True`` (the
        member set of a non-spin-orbital archive) is refused."""
        self._assert_chi0q_marker_refused(self.d, self.files, self.flex, "u")

    def test_unmarked_chi0q_companion_is_refused_at_zero_vertex(self):
        """Gamma = 0 (U = 0): the stored chi equals chi0, so the vertex
        consistency check alone would accept an unmarked chi0q of the same
        values; the marker refuses it."""
        d = tempfile.mkdtemp(prefix="sc_so_flex_g0_")
        self.addCleanup(shutil.rmtree, d, True)
        files = _write_one_orbital(d, U=0.0)
        flex = _run_flex(d, files, "flex")
        with np.load(os.path.join(flex, "chiq.npz")) as z, \
                np.load(os.path.join(flex, "chi0q.npz")) as z0:
            si = z["chiq"].shape[0] // 2
            nvol = z["chiq"].shape[1]
            chi = z["chiq"][si].reshape(nvol, 4, 4)
            chi0 = z0["chi0q"][si].reshape(nvol, 4, 4)
        # the vertex check passes vacuously here
        with _quiet():
            self.assertLess(sc._check_so_flex_vertex(
                np.zeros((4, 4)), chi, chi0, "chiq.npz"), sc._SO_FLEX_GAMMA_TOL)
        self._assert_chi0q_marker_refused(d, files, flex, "g0")

    def test_refuses_non_so_and_reduced(self):
        """A non-spin-orbital general FLEX output ("myo") and a reduced one
        ("kuroki") are refused, before the green or any vertex is built."""
        d = tempfile.mkdtemp(prefix="sc_so_flex_nonso_")
        self.addCleanup(shutil.rmtree, d, True)
        nfiles = _write_model(d, so=False, types=("CoulombIntra", "CoulombInter"))
        for scheme, tag in (("general", "myo"), ("reduced", "kuroki")):
            flex = _run_flex(d, nfiles, "flex_" + scheme, so=False, scheme=scheme)
            with np.load(os.path.join(flex, "chiq.npz")) as z:
                self.assertEqual(str(z["chi_convention"]), tag)
            with self.subTest(scheme=scheme), mock.patch.object(
                    sc, "_load_flex_green",
                    side_effect=AssertionError("green read before the check")), \
                    mock.patch.object(
                        sc, "_so_rpa_vertex",
                        side_effect=AssertionError("vertex built before the check")):
                self._refused(r"chi_convention = '{}'".format(tag),
                              _flex_input(self.d, self.files, flex,
                                          out=self._out("out_" + scheme)))


#: npz member sets of the NON-spin-orbital chi0q / chiq archives, recorded
#: from the revision before the spin-orbital marker (796738aa) by running the
#: same inputs as TestArchiveMembers; the marker must not change them
_RPA_MEMBERS = ["calc_scheme", "calc_scheme_requested", "coeff_tail",
                "freq_index", "index_convention", "momentum_convention", "nmat",
                "scheme_resolution", "tail_endpoint", "wavevector_index",
                "wavevector_unit"]
_FLEX_PHASE_B = ["density_target_enforced", "flex_guard_policy",
                 "flex_guard_violations", "flex_second_order",
                 "flex_second_order_schema", "hf_density_error",
                 "hf_density_source", "map_iteration", "payload_kind",
                 "scf_component_residual", "scf_converged",
                 "scf_green_residual", "scf_iterations", "scf_sigma_residual",
                 "state_iteration"]
_NON_SO_MEMBERS = {
    ("RPA", "general", "chiq"): sorted(_RPA_MEMBERS + ["chiq"]),
    ("RPA", "general", "chi0q"): sorted(_RPA_MEMBERS + ["chi0q"]),
    ("RPA", "reduced", "chiq"): sorted(_RPA_MEMBERS + ["chiq"]),
    ("RPA", "reduced", "chi0q"): sorted(_RPA_MEMBERS + ["chi0q"]),
    ("FLEX", "general", "chi0q"): sorted(_RPA_MEMBERS + _FLEX_PHASE_B + ["chi0q"]),
    ("FLEX", "general", "chiq"): sorted(
        ["calc_scheme", "calc_scheme_requested", "chi_convention",
         "chi_orbital_layout", "chiq_c", "chiq_s", "freq_index",
         "momentum_convention", "nmat", "sc_vertex_version",
         "scheme_resolution", "wavevector_index", "wavevector_unit"]
        + _FLEX_PHASE_B),
    ("FLEX", "reduced", "chi0q"): sorted(_RPA_MEMBERS + ["chi0q"]),
    ("FLEX", "reduced", "chiq"): sorted(
        ["calc_scheme", "calc_scheme_requested", "chi_convention", "chiq_c",
         "chiq_s", "freq_index", "momentum_convention", "nmat",
         "sc_vertex_version", "scheme_resolution", "wavevector_index",
         "wavevector_unit"]),
}


class TestChi0qInitProvenance(_TmpDir):
    """A spin-orbital RPA run PROPAGATES the ``spin_orbital`` marker of a
    chi0q it did not compute (``chi0q_init`` file or in-memory green_info
    reuse): marked source -> marked output; unmarked source -> the run
    proceeds, the output has no marker (hwave_sc SO load refuses it);
    explicit ``False`` / malformed -> refused before solve."""

    def _rpa(self, d, files, so, out, chi0q_init=None, capture=None,
             inject=None):
        from hwave import qlms
        from hwave.solver.rpa import RPA
        inp = _input(d, files, so, {}, out=out)
        inp["mode"]["calc_scheme"] = "general"
        inp["file"]["input"]["path_to_input"] = d
        inp["file"]["output"].update({"chiq": "chiq", "chi0q": "chi0q"})
        if chi0q_init is not None:
            inp["file"]["input"]["chi0q_init"] = chi0q_init
        del inp["eliashberg"]
        inp["log"] = {"print_level": 0}
        real = RPA.save_results

        def spy(slf, info_outputfile, green_info):
            if capture is not None:
                capture.update(solver=slf, green_info=green_info,
                               outputs=info_outputfile)
            if inject is not None:
                inject(slf, green_info)
            return real(slf, info_outputfile, green_info)

        with mock.patch.object(RPA, "save_results", autospec=True,
                               side_effect=spy):
            qlms.run(input_dict=inp)
        return os.path.join(d, out)

    def _non_so_chi0q(self):
        """A real non-spin-orbital general chi0q over ND orbitals (no
        marker), the shape of a spin-orbital one over NP orbitals."""
        src = os.path.join(self.d, "nonso")
        os.makedirs(src)
        files = _write_model(src, so=True)   # read as a plain ND-orbital model
        with _quiet():
            path = self._rpa(src, files, False, "out")
        path = os.path.join(path, "chi0q.npz")
        with np.load(path) as z:
            self.assertNotIn("spin_orbital", z.files)
        return path

    def _marker(self, out):
        res = {}
        for name in ("chiq", "chi0q"):
            with np.load(os.path.join(out, name + ".npz")) as z:
                res[name] = (bool(z["spin_orbital"]) if "spin_orbital" in z.files
                             else None)
        return res

    def _so_load_refused(self, out, pattern):
        files = _write_model(self.d, so=True)
        inp = _input(self.d, files, True, {"chi0q_mode": "load",
                                           "solver_mode": "eigenvalue"})
        inp["file"]["output"]["path_to_output"] = out
        with self.assertRaisesRegex(ValueError, pattern), _quiet(), \
                mock.patch.object(sc, "_build_bond_green",
                                  side_effect=AssertionError("Green built")):
            sc.calc_eliashberg(inp)

    def test_unmarked_chi0q_init_gives_unmarked_output(self):
        src = self._non_so_chi0q()
        files = _write_model(self.d, so=True)
        with self.assertLogs("hwave.solver.rpa", level="INFO") as logs:
            out = self._rpa(self.d, files, True, "so_out", chi0q_init=src)
        self.assertTrue(any("no spin_orbital = True marker" in m
                            and "hwave_sc cannot load them" in m
                            for m in logs.output))
        self.assertEqual(self._marker(out), {"chiq": None, "chi0q": None})
        self._so_load_refused(out, r"has no spin_orbital marker")

    def test_false_or_malformed_chi0q_init_refused_before_solve(self):
        from hwave.solver.rpa import RPA
        src = self._non_so_chi0q()
        files = _write_model(self.d, so=True)
        for tag, value, pattern in (
                ("false", False, r"records spin_orbital = False"),
                ("str", "true", r"malformed spin_orbital"),
                ("array", np.array([True]), r"malformed spin_orbital")):
            with self.subTest(marker=tag):
                bad = os.path.join(self.d, "bad_{}.npz".format(tag))
                _resave(src, bad, spin_orbital=value)
                with self.assertRaisesRegex(ValueError, pattern), _quiet(), \
                        mock.patch.object(
                            RPA, "solve",
                            side_effect=AssertionError("solve reached")):
                    self._rpa(self.d, files, True, "so_" + tag, chi0q_init=bad)

    def test_marked_chi0q_init_keeps_marker(self):
        files = _write_model(self.d, so=True, soc=0.3)
        with _quiet():
            first = self._rpa(self.d, files, True, "first")
        self.assertEqual(self._marker(first), {"chiq": True, "chi0q": True})
        with _quiet():
            second = self._rpa(self.d, files, True, "second",
                               chi0q_init=os.path.join(first, "chi0q.npz"))
        self.assertEqual(self._marker(second), {"chiq": True, "chi0q": True})
        with np.load(os.path.join(first, "chi0q.npz")) as a, \
                np.load(os.path.join(second, "chi0q.npz")) as b:
            np.testing.assert_array_equal(a["chi0q"], b["chi0q"])

    def test_save_without_solve_writes_no_marker(self):
        """``save_results`` on an UNSOLVED spin-orbital solver, given an
        arbitrary chi0q/chiq (also with forged provenance that is not bound
        to the array), must not emit the marker: only solve() establishes
        it."""
        from hwave.solver.rpa import RPA
        files = _write_model(self.d, so=True, soc=0.3)
        nvol = int(np.prod(CELL))
        rng = np.random.default_rng(5)
        arr = rng.standard_normal((64, nvol) + (ND,) * 4) + 0j
        for tag, meta in (("plain", None),
                          ("forged", {"spin_orbital": True}),
                          ("stale", {"spin_orbital": True,
                                     "fingerprint": "0" * 32})):
            def inject(slf, gi, meta=meta):
                self.assertIsNone(getattr(slf, "_chi0q_init_meta", None))
                gi["chi0q"] = arr
                gi["chiq"] = arr.copy()
                if meta is not None:
                    gi["chi0q_freq_meta"] = dict(meta)
            with self.subTest(case=tag), _quiet(), mock.patch.object(
                    RPA, "solve", autospec=True, return_value=None):
                out = self._rpa(self.d, files, True, "nosolve_" + tag,
                                inject=inject)
            self.assertEqual(self._marker(out), {"chiq": None, "chi0q": None})

    def test_in_memory_reuse(self):
        files = _write_model(self.d, so=True, soc=0.3)
        cap = {}
        with _quiet():
            self._rpa(self.d, files, True, "mem0", capture=cap)
        solver, gi, outputs = cap["solver"], cap["green_info"], cap["outputs"]
        self.assertIs(gi["chi0q_freq_meta"]["spin_orbital"], True)

        def reuse(out, meta_update=None, drop=False):
            g = dict(gi)
            meta = dict(gi["chi0q_freq_meta"])
            if drop:
                meta.pop("spin_orbital")
            meta.update(meta_update or {})
            g["chi0q_freq_meta"] = meta
            path = os.path.join(self.d, out)
            os.makedirs(path)
            o = dict(outputs, path_to_output=path)
            solver.solve(g, path)
            solver.save_results(o, g)
            return path, g

        with _quiet():
            out, _ = reuse("mem_marked")
        self.assertEqual(self._marker(out), {"chiq": True, "chi0q": True})
        with self.assertLogs("hwave.solver.rpa", level="INFO") as logs:
            out, g = reuse("mem_unmarked", drop=True)
        self.assertTrue(any("no spin_orbital = True marker" in m
                            for m in logs.output))
        self.assertEqual(self._marker(out), {"chiq": None, "chi0q": None})
        # a further reuse of the unmarked bubble stays unmarked
        self.assertNotIn("spin_orbital", g["chi0q_freq_meta"])
        with _quiet():
            path = os.path.join(self.d, "mem_chain")
            os.makedirs(path)
            solver.solve(g, path)
            solver.save_results(dict(outputs, path_to_output=path), g)
        self.assertEqual(self._marker(path), {"chiq": None, "chi0q": None})
        for tag, value, pattern in (
                ("false", False, r"records spin_orbital = False"),
                ("str", "yes", r"malformed spin_orbital"),
                ("none", None, r"malformed spin_orbital")):
            with self.subTest(marker=tag), \
                    self.assertRaisesRegex(ValueError, pattern), _quiet(), \
                    mock.patch.object(solver, "_solve_rpa",
                                      side_effect=AssertionError("solved")):
                reuse("mem_" + tag, {"spin_orbital": value})


class TestArchiveMembers(_TmpDir):
    """The ``spin_orbital`` marker is written ONLY by spin-orbital runs: a
    non-spin-orbital RPA / FLEX archive keeps exactly its previous member
    set; a spin-orbital RPA run writes ``spin_orbital = True``."""

    def _rpa(self, files, so, scheme, out):
        from hwave import qlms
        inp = _input(self.d, files, so, {}, out=out)
        inp["mode"]["calc_scheme"] = scheme
        inp["file"]["input"]["path_to_input"] = self.d
        inp["file"]["output"].update({"chiq": "chiq", "chi0q": "chi0q"})
        del inp["eliashberg"]
        inp["log"] = {"print_level": 0}
        with _quiet():
            qlms.run(input_dict=inp)
        return os.path.join(self.d, out)

    def test_non_so_member_sets_unchanged(self):
        files = _write_model(self.d, so=False, types=("CoulombIntra", "CoulombInter"))
        for scheme in ("general", "reduced"):
            dirs = {"RPA": self._rpa(files, False, scheme, "rpa_" + scheme),
                    "FLEX": _run_flex(self.d, files, "flex_" + scheme, so=False,
                                      scheme=scheme)}
            for solver, path in dirs.items():
                for name in ("chiq", "chi0q"):
                    with self.subTest(solver=solver, scheme=scheme, file=name), \
                            np.load(os.path.join(path, name + ".npz")) as z:
                        self.assertEqual(sorted(z.files),
                                         _NON_SO_MEMBERS[(solver, scheme, name)])

    def test_so_rpa_writes_the_marker(self):
        files = _write_model(self.d, so=True, soc=0.3)
        path = self._rpa(files, True, "general", "rpa_so")
        for name in ("chiq", "chi0q"):
            with self.subTest(file=name), \
                    np.load(os.path.join(path, name + ".npz")) as z:
                m = np.asarray(z["spin_orbital"])
                self.assertEqual((m.ndim, m.dtype), (0, np.bool_))
                self.assertTrue(bool(m))
                self.assertEqual(sorted(set(z.files) - {"spin_orbital"}),
                                 _NON_SO_MEMBERS[("RPA", "general", name)])


class TestInitialGap(unittest.TestCase):

    K = [np.linspace(0, 2 * np.pi, n, endpoint=False) for n in CELL]

    def test_even_seed_is_singlet_structured(self):
        gap = so_pairing.so_initial_gap("cos", NP, *self.K)
        proj = sc._project_gap_parity(gap, "triplet")
        self.assertGreater(np.linalg.norm(proj), 0.99)
        np.testing.assert_allclose(proj, gap, atol=1e-14)
        w = _spin_weights(gap, NP)
        self.assertAlmostEqual(w[0], 1.0, places=12)

    def test_odd_seed_is_sz0_triplet(self):
        gap = so_pairing.so_initial_gap("p_x", NP, *self.K)
        np.testing.assert_allclose(sc._project_gap_parity(gap, "triplet"), gap,
                                   atol=1e-14)
        w = _spin_weights(gap, NP)
        self.assertAlmostEqual(w[3], 1.0, places=12)

    def test_random_seed_spans_all_spin_blocks(self):
        gap = so_pairing.so_initial_gap("random", NP, *self.K)
        self.assertEqual(gap.shape, (ND, ND) + tuple(CELL))
        self.assertTrue(np.all(np.abs(gap.reshape(2, NP, 2, NP, -1)).sum(
            axis=(1, 3, 4)) > 0))


def _hand_gap(psi, dx, dy, dz):
    """``(2N, 2N, ...)`` spin-block gap from components ``(N, N, ...)``,
    with the spin blocks written out explicitly (rows/columns up, down):
    ``(psi + d . sigma) i sigma_y = [[-dx + i dy, psi + dz],
    [-psi + dz, dx + i dy]]``."""
    psi, dx, dy, dz = (np.asarray(x, complex) for x in (psi, dx, dy, dz))
    n = psi.shape[0]
    g = np.zeros((2, n, 2, n) + psi.shape[2:], complex)
    g[0, :, 0] = -dx + 1j * dy
    g[0, :, 1] = psi + dz
    g[1, :, 0] = -psi + dz
    g[1, :, 1] = dx + 1j * dy
    return g.reshape((2 * n, 2 * n) + psi.shape[2:])


class TestSpinDecomposition(unittest.TestCase):

    SHAPE = (NP, NP, 3, 2, 1)

    def setUp(self):
        rng = np.random.default_rng(83)
        self.f = [rng.normal(size=self.SHAPE) + 1j * rng.normal(size=self.SHAPE)
                  for _ in range(4)]
        self.zero = np.zeros(self.SHAPE, complex)

    def _check_pure(self, slot):
        parts = [self.f[i] if i == slot else self.zero for i in range(4)]
        gap = _hand_gap(*parts)
        comps = so_pairing.so_spin_components(gap, NP)
        for i, c in enumerate(comps):
            np.testing.assert_allclose(c, parts[i], rtol=0, atol=1e-14)
        w = so_pairing.so_spin_decomposition(gap, NP)
        self.assertEqual(tuple(w), so_pairing.SPIN_COMPONENTS)
        for i, key in enumerate(so_pairing.SPIN_COMPONENTS):
            self.assertAlmostEqual(w[key], 1.0 if i == slot else 0.0, places=14)

    def test_pure_psi(self):
        self._check_pure(0)

    def test_pure_dx(self):
        self._check_pure(1)

    def test_pure_dy(self):
        self._check_pure(2)

    def test_pure_dz(self):
        self._check_pure(3)

    def test_mixture(self):
        amps = (0.5, 1.0, np.sqrt(2.0), 1.5)
        parts = [a * fi / np.linalg.norm(fi) for a, fi in zip(amps, self.f)]
        gap = _hand_gap(*parts)
        comps = so_pairing.so_spin_components(gap, NP)
        for c, ref in zip(comps, parts):
            np.testing.assert_allclose(c, ref, rtol=0, atol=1e-14)
        w = so_pairing.so_spin_decomposition(gap, NP)
        tot = sum(a * a for a in amps)
        for key, a in zip(so_pairing.SPIN_COMPONENTS, amps):
            self.assertAlmostEqual(w[key], a * a / tot, places=14)
        self.assertAlmostEqual(sum(w.values()), 1.0, places=14)
        # |Delta|_F^2 = 2 (|psi|^2 + |d|^2)
        self.assertAlmostEqual(np.linalg.norm(gap) ** 2, 2 * tot, places=12)
        # invariant under the overall phase and normalization of the gap
        w2 = so_pairing.so_spin_decomposition(-3.0j * gap, NP)
        for key in w:
            self.assertAlmostEqual(w2[key], w[key], places=14)
        # the independent test helper agrees
        np.testing.assert_allclose(_spin_weights(gap, NP),
                                   [w[k] for k in so_pairing.SPIN_COMPONENTS],
                                   rtol=0, atol=1e-14)

    def test_bad_shape_and_zero_gap(self):
        with self.assertRaisesRegex(ValueError, "2 N, 2 N"):
            so_pairing.so_spin_decomposition(np.zeros((3, 3, 2)), 2)
        w = so_pairing.so_spin_decomposition(np.zeros((ND, ND, 2)), NP)
        self.assertTrue(all(np.isnan(v) for v in w.values()))


class TestOutputs(_TmpDir):

    def test_header_and_weights(self):
        """SU(2) limit: the leading (singlet) eigenvector has psi weight 1, a
        triplet one d weight 1; eigenvalue.dat carries the spin-orbital header
        lines and the four weight columns; gap.npz holds the
        ``(nd, nd, Nx, Ny, Nz)`` gap with its index_convention tag."""
        files = _write_model(self.d, so=True)
        inp = _input(self.d, files, True,
                     {"chi0q_mode": "calc", "solver_mode": "eigenvalue",
                      "num_eigenvalues": 24})
        seen = _run(inp)
        evs, vecs = seen["evs"], seen["vecs"]
        weights = [so_pairing.so_spin_decomposition(v, NP) for v in vecs]
        self.assertAlmostEqual(weights[0]["psi"], 1.0, delta=1e-10)
        trip = [i for i, w in enumerate(weights)
                if abs(w["dx"] + w["dy"] + w["dz"] - 1.0) <= 1e-10
                and abs(evs[i]) > 1e-8 * abs(evs[0])]
        self.assertTrue(trip)
        # every physical eigenvector is either a pure singlet or a pure triplet
        for i, w in enumerate(weights):
            if abs(evs[i]) > 1e-8 * abs(evs[0]):
                with self.subTest(index=i):
                    self.assertLess(min(w["psi"], 1.0 - w["psi"]), 1e-8)
        out = inp["file"]["output"]["path_to_output"]
        with open(os.path.join(out, "eigenvalue.dat")) as f:
            lines = f.read().splitlines()
        self.assertIn("# spin_orbital = true", lines)
        self.assertIn("# index_convention = spin_block", lines)
        head = [ln for ln in lines if ln.startswith("# index  Re")]
        self.assertEqual(len(head), 1)
        self.assertTrue(head[0].split()[-4:] == ["w_psi", "w_dx", "w_dy", "w_dz"])
        start = lines.index("# Eigenvalue analysis")
        rows = np.loadtxt(lines[start:])
        self.assertEqual(rows.shape, (len(evs), 9))
        np.testing.assert_allclose(rows[:, 1], evs.real, rtol=1e-8,
                                   atol=1e-8 * abs(evs[0]))
        ref = np.array([[w[k] for k in so_pairing.SPIN_COMPONENTS]
                        for w in weights])
        np.testing.assert_allclose(rows[:, 5:], ref, rtol=1e-8, atol=1e-15)
        with open(os.path.join(out, "gap.dat")) as f:
            gap_head = [next(f) for _ in range(4)]
        self.assertEqual(gap_head[0], "# spin_orbital = true\n")
        self.assertEqual(gap_head[1], "# index_convention = spin_block\n")
        self.assertTrue(gap_head[3].startswith("# kx ky kz"))
        with np.load(os.path.join(out, "gap.npz")) as z:
            self.assertEqual(str(z["index_convention"]), "spin_block")
            self.assertTrue(bool(z["spin_orbital"]))
            self.assertEqual(int(z["norb_phys"]), NP)
            self.assertEqual(z["gap"].shape, (ND, ND) + tuple(CELL))
            np.testing.assert_array_equal(z["gap"], seen["gap"])
            np.testing.assert_array_equal(z["eigenvalues"], evs)
            self.assertEqual(z["eigenvectors"].shape,
                             (len(evs), ND, ND) + tuple(CELL))
            np.testing.assert_allclose(z["eigenvector_spin_weights"], ref,
                                       rtol=0, atol=1e-15)
            np.testing.assert_allclose(z["gap_spin_weights"], ref[0],
                                       rtol=0, atol=1e-15)
            self.assertEqual(list(z["spin_components"]),
                             ["psi", "dx", "dy", "dz"])

    def test_iteration_only_writes_gap_weights(self):
        files = _write_model(self.d, so=True, soc=0.3)
        inp = _input(self.d, files, True,
                     {"chi0q_mode": "calc", "solver_mode": "iteration",
                      "max_iter": 50, "output_gap": "g.txt"})
        seen = _run(inp)
        out = inp["file"]["output"]["path_to_output"]
        prov = seen["provenance"]
        w = so_pairing.so_spin_decomposition(seen["gap"], NP)
        self.assertEqual(prov["gap_npz"], "g.npz")
        self.assertEqual(prov["gap_spin_weights"], " ".join(
            "{}={:.8e}".format(k, w[k]) for k in so_pairing.SPIN_COMPONENTS))
        with np.load(os.path.join(out, "g.npz")) as z:
            np.testing.assert_array_equal(z["gap"], seen["gap"])
            self.assertNotIn("eigenvalues", z.files)

    def test_gap_archive_colliding_with_eigenvalue_file(self):
        files = _write_model(self.d, so=True)
        inp = _input(self.d, files, True,
                     {"chi0q_mode": "calc", "solver_mode": "eigenvalue",
                      "output_gap": "res.dat", "output_eigenvalue": "res.npz"})
        with _quiet(), mock.patch.object(
                sc, "_read_interaction_files",
                side_effect=AssertionError("input read before the check")):
            with self.assertRaisesRegex(ValueError,
                                        "output_gap.*output_eigenvalue"):
                sc.calc_eliashberg(inp)
        self.assertFalse(os.path.exists(inp["file"]["output"]["path_to_output"]))

    def test_npz_name(self):
        self.assertEqual(sc._so_gap_npz_name("gap.dat"), "gap.npz")
        self.assertEqual(sc._so_gap_npz_name("gap"), "gap.npz")
        self.assertEqual(sc._so_gap_npz_name("gap.npz"), "gap.npz.npz")


if __name__ == "__main__":
    unittest.main()
