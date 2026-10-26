"""The static spin-orbital entry of hwave_sc (``[mode] enable_spin_orbital``,
issue #83): ``chi0q_mode = "calc"``, ``"load"`` and ``"flex"``.

The gap is one ``2N x 2N`` spin-orbital matrix ``Delta_AB(k)`` (spin-block
order ``A = s * N + a``); the kernel is ``P K`` with the pairing vertex of
``hwave.solver.so_pairing`` and ``P`` the projection onto fermion-
antisymmetric gaps ``Delta_AB(k) = -Delta_BA(-k)``, applied after every
kernel application.

* ``TestStaticEntry``: the SU(2) limit reproduces the non-spin-orbital
  singlet and (threefold) triplet eigenvalues (G1, end to end); a model with
  complex Rashba-like spin-orbit hopping runs and returns an antisymmetric
  gap; ``chi0q_mode = "load"`` of a spin-orbital RPA chi0q equals
  ``chi0q_mode = "calc"`` (G4, part).
* ``TestFlexMode``: ``chi0q_mode = "flex"`` uses the stored chi (zero bosonic
  frequency) and dressed green of a spin-orbital FLEX run (G4, part); the
  archive is identified by its tags (nd = 2 has the shape of a two-orbital
  non-spin-orbital file); non-spin-orbital, reduced and untagged outputs and
  a vertex the FLEX run did not use are refused.
* ``TestRefusals``: every unsupported combination is refused through the
  public entry, with the key named (G5).
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


class _TmpDir(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="sc_so_")
        self.addCleanup(shutil.rmtree, self.d, True)


class TestStaticEntry(_TmpDir):

    def _leading_physical(self, seen):
        evs, match = seen["evs"], seen["match"]
        self.assertTrue(match[0])
        return evs[0]

    @staticmethod
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

    def test_su2_model_reproduces_singlet_and_triplet(self):
        """G1: with spin-independent hopping and a rotation-invariant
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
        evs = so["evs"]
        scale = np.abs(evs).max()
        nonzero = np.abs(evs) > 1e-10 * scale
        self.assertTrue(np.all(so["match"][nonzero]))
        # expected multiplicity of each value: even singlet x 1 + odd triplet x 3
        expected, forbidden = {}, []
        floor = {ch: np.abs(res[ch]["evs"]).min() for ch in res}
        for ch, allowed_even in (("singlet", True), ("triplet", False)):
            for v, n, d_even, d_odd in self._clusters(res[ch]["evs"],
                                                      res[ch]["vecs"]):
                self.assertEqual(d_even + d_odd, n)
                if abs(v) <= (1 + 1e-6) * floor[ch]:
                    continue        # a cluster cut by the truncation
                d_ok, d_bad = (d_even, d_odd) if allowed_even else (d_odd, d_even)
                key = next((k for k in expected
                            if abs(k - v) <= 1e-8 * abs(v)), v)
                expected[key] = expected.get(key, 0) + d_ok * (
                    1 if ch == "singlet" else 3)
                if d_bad:
                    forbidden.append((ch, v))
        self.assertGreater(len(expected), 20)
        for v, count in expected.items():
            got = int(np.sum(np.abs(evs - v) <= 1e-8 * abs(v)))
            with self.subTest(value=v):
                self.assertEqual(got, count)
        # the large Pauli-forbidden triplet modes (even k-parity) are absent
        big = sorted((v for ch, v in forbidden
                      if ch == "triplet" and v.real > 0.5), key=lambda x: -x.real)
        self.assertGreaterEqual(len(big), 2)
        for v in big:
            with self.subTest(forbidden=v):
                if any(abs(k - v) <= 1e-8 * abs(v) for k in expected):
                    continue
                self.assertFalse(np.any(np.abs(evs - v) <= 1e-6 * abs(v)))

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
        """G4 (part): the general chi0q of a spin-orbital RPA run, loaded,
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
    """G5: each unsupported combination is refused through the public entry
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
                 index_convention="spin_block", momentum_convention="e_plus_ikR")
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

    def test_extern_file(self):
        files = _write_model(self.d, so=True)
        _wfile(os.path.join(self.d, "extern.dat"), "Extern", NP,
               [(0, 0, 0, a + 1, a + 1, 0.1, 0.0) for a in range(NP)])
        for key in ("Extern", "EXTERN"):
            with self.subTest(key=key):
                self._refused(r"Extern", so_files=dict(files, **{key: "extern.dat"}))

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
        """G4: the eigenvalues equal those of ``P K`` built from the stored
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
