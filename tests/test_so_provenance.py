"""Version-1 provenance of the spin-orbital RPA and FLEX archives (#88, PR 0).

A spin-orbital RPA run writes the fingerprint and state members into chi0q
and chiq; a spin-orbital Hartree-Fock FLEX run writes them into chiq, green,
sigma and chi0q, with one Green digest (that of the G in green.npz) in all
four. What this module pins:

* :class:`TestMembers` -- ``fingerprint.so_provenance_members``: the member
  schema (dtype, shape), the digests recomputable from the stored ``prov_*``
  tables, the ``calc_mu`` tagged union of the state, the refusal of a
  non-finite member, ``producer_version`` (truncation, git describe cached
  and never failing).
* :class:`TestProducersWrite` -- the RPA and FLEX producers through
  ``hwave.qlms.run``: the band digest equals the one computed from the input
  files, the state is the realized one (closure density in electrons per
  cell over all spin-orbitals, a fixed-mu density against a closed-form
  Fermi sum), the PRE-FOLD tables are hashed under a sublattice, the members
  are built once per run, mixed-case keys, and the non-spin-orbital archives
  are unchanged against the snapshot in ``tests/so_provenance_ref`` (made by
  ``tests/so_provenance_ref/make.py`` on the tree before this change).
"""
import contextlib
import glob
import logging
import os
import shutil
import tempfile
import unittest
from unittest import mock

import numpy as np

import hwave
import hwave.qlmsio.wan90 as wan90
from hwave.solver import fingerprint as fp
import tests.test_so_zeeman_route as Z
from tests.so_provenance_ref import make as ref_make

_REF_DIR = os.path.dirname(os.path.abspath(ref_make.__file__))

_SCHEMA = {
    "fingerprint_version": ("i", 0),
    "producer_stage": ("U", 0), "producer_label": ("U", 0), "producer_version": ("U", 0),
    "fp_transfer": ("U", 0), "fp_extern": ("U", 0), "fp_meta": ("U", 0),
    "fp_band": ("U", 0), "fp_green": ("U", 0),
    "prov_transfer_keys": ("i", 2), "prov_extern_keys": ("i", 2),
    "prov_transfer_vals": ("c", 1), "prov_extern_vals": ("c", 1),
    "prov_extern_present": ("b", 0), "prov_enable_spin_orbital": ("b", 0),
    "state_calc_mu": ("b", 0),
    "prov_coeff_extern": ("f", 0), "state_beta": ("f", 0), "state_mu": ("f", 0),
    "state_density": ("f", 0), "state_target": ("f", 0), "state_mu_tol": ("f", 0),
    "prov_norb": ("i", 0), "prov_norb_phys": ("i", 0), "state_nmat": ("i", 0),
    "prov_index_convention": ("U", 0),
    "prov_cell_shape": ("i", 1), "prov_sub_shape": ("i", 1),
}
_ITEMSIZE = {"i": 8, "f": 8, "c": 16, "b": 1}


@contextlib.contextmanager
def _quiet():
    prev = logging.root.manager.disable
    logging.disable(logging.WARNING)
    try:
        yield
    finally:
        logging.disable(prev)


def _table_from_prov(keys, vals):
    return {(tuple(int(x) for x in k[:3]), tuple(int(x) for x in k[3:])): complex(v)
            for k, v in zip(keys, vals)}


def _meta_from_prov(z):
    return fp.BandMeta(float(z["prov_coeff_extern"]), bool(z["prov_enable_spin_orbital"]),
                       str(z["prov_index_convention"]), int(z["prov_norb"]),
                       int(z["prov_norb_phys"]), tuple(int(x) for x in z["prov_cell_shape"]),
                       tuple(int(x) for x in z["prov_sub_shape"]))


# --------------------------------------------------------------- runners
_RPA_EXTRA = {"calc_scheme": "general", "calc_type": "ring"}


def _so_input(d, *, cell=(4, 1, 1), sub_shape=(1, 1, 1), mu=None, extern_key="Extern",
              coeff_key="coeff_extern"):
    """The 0-E1 spin-orbital fixture (one physical orbital, complex SOC,
    Extern field ``H_11(0) = 0.1``) written into ``d``; returns the input
    description used by the runners and the test-side oracles."""
    inter = Z._so_1orb_fixture(d, "extern", extern_key=extern_key)
    param = dict(Z._E1_PARAM, CellShape=list(cell), SubShape=list(sub_shape))
    param.pop("coeff_extern")
    param[coeff_key] = 1.0
    if mu is not None:
        param.pop("filling")
        param["mu"] = mu
    return {"dir": d, "inter": inter, "param": param, "cell": tuple(cell),
            "sub_shape": tuple(sub_shape)}


def _run(inp, mode, outputs, extra):
    import hwave.qlms as qlms
    d = inp["dir"]
    out_dir = os.path.join(d, "output_" + mode)
    mode_tbl = dict({"mode": mode, "enable_spin_orbital": True,
                     "param": dict(inp["param"])}, **extra)
    run_inp = {"log": {"print_level": 0}, "mode": mode_tbl,
               "file": {"input": {"path_to_input": d,
                                  "interaction": dict(inp["inter"], path_to_input=d)},
                        "output": dict(outputs, path_to_output=out_dir)}}
    with _quiet():
        qlms.run(input_dict=run_inp)
    res = {}
    for k, f in outputs.items():
        with np.load(os.path.join(out_dir, f)) as z:
            res[k] = {m: np.array(z[m]) for m in z.files}
    return res


def _run_rpa_so(inp):
    return _run(inp, "RPA", {"chi0q": "chi0q.npz", "chiq": "chiq.npz"}, _RPA_EXTRA)


def _run_flex_so(inp, iteration_max=3):
    param = dict(Z._F1_PARAM, IterationMax=iteration_max, coeff_extern=1.0,
                 CellShape=list(inp["cell"]))
    return _run(dict(inp, param=param), "FLEX",
                {"chiq": "chiq.npz", "green": "green.npz", "sigma": "sigma.npz",
                 "chi0q": "chi0q.npz"}, {"calc_scheme": "general"})


def _tables_and_meta_from_input(inp):
    """``read_w90`` of the input Transfer and Extern files, and the band
    metadata of the fixture (2 spin-orbitals = 1 physical orbital)."""
    d, inter = inp["dir"], inp["inter"]
    ext_key = [k for k in inter if k.lower() == "extern"]
    T = wan90.read_w90(os.path.join(d, inter["Transfer"]))
    E = wan90.read_w90(os.path.join(d, inter[ext_key[0]])) if ext_key else None
    meta = fp.BandMeta(1.0, True, "spin_block", 2, 1, inp["cell"], inp["sub_shape"])
    return T, E, meta


def _closed_form_fermi_sum(inp, mu):
    """Total electron count of the band at ``mu`` from an independent
    eigendecomposition: H(k) = T(k) + h sigma^z (interleaved 2a+s)."""
    T, E, _ = _tables_and_meta_from_input(inp)
    beta = 1.0 / inp["param"]["T"]
    L = inp["cell"]
    total = 0.0
    for kx in range(L[0]):
        for ky in range(L[1]):
            for kz in range(L[2]):
                k = 2 * np.pi * np.array([kx / L[0], ky / L[1], kz / L[2]])
                H = np.zeros((2, 2), dtype=complex)
                for (R, (a, b)), v in T.items():
                    H[a, b] += v * np.exp(1j * np.dot(k, R))
                for (R, (a, b)), v in E.items():
                    ph = v * np.exp(1j * np.dot(k, R))
                    H[2 * a, 2 * b] += ph
                    H[2 * a + 1, 2 * b + 1] -= ph
                w = np.linalg.eigvalsh(H)
                total += np.sum(1.0 / (1.0 + np.exp(beta * (w - mu))))
    return total


def _reference_outputs(name):
    """The committed snapshot of run ``name``: ``{archive file: NpzFile}``."""
    d = os.path.join(_REF_DIR, name.lstrip("_"))
    return {os.path.basename(f): np.load(f) for f in sorted(glob.glob(os.path.join(d, "*.npz")))}


def _ref_platform():
    with open(os.path.join(_REF_DIR, "platform.txt")) as f:
        return f.read().strip()


# ----------------------------------------------------------------- tests
class _Dir(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="so_prov_")
        self.addCleanup(shutil.rmtree, self.d, True)

    def new_dir(self):
        d = tempfile.mkdtemp(dir=self.d)
        return d


class TestMembers(_Dir):
    def setUp(self):
        super().setUp()
        inp = _so_input(self.new_dir())
        T, E, meta = _tables_and_meta_from_input(inp)
        self.kw = dict(producer_stage="rpa_so", producer_label="test producer",
                       transfer=T, extern=E, meta=meta, fp_green="0" * 64,
                       state=fp.so_state(beta=2.0, nmat=16, mu=-0.3, density=0.8,
                                         calc_mu=True, target=0.8 + 1e-14))

    def test_schema(self):
        m = fp.so_provenance_members(**self.kw)
        self.assertEqual(set(m), set(_SCHEMA))
        self.assertEqual(tuple(m), fp.SO_PROVENANCE_MEMBERS)
        for k, (kind, ndim) in _SCHEMA.items():
            with self.subTest(member=k):
                self.assertEqual(m[k].dtype.kind, kind)
                self.assertEqual(m[k].ndim, ndim)
                if kind in _ITEMSIZE:
                    self.assertEqual(m[k].dtype.itemsize, _ITEMSIZE[kind])
        for k in ("prov_transfer_keys", "prov_extern_keys"):
            self.assertEqual(m[k].shape[1], 5)
        self.assertEqual(m["prov_transfer_keys"].shape[0], m["prov_transfer_vals"].shape[0])
        self.assertEqual(m["prov_cell_shape"].shape, (3,))
        for k in ("fp_transfer", "fp_extern", "fp_meta", "fp_band"):
            self.assertRegex(str(m[k]), r"^[0-9a-f]{64}$")
        self.assertEqual(int(m["fingerprint_version"]), 1)

    def test_digests_recomputable_from_prov(self):
        m = fp.so_provenance_members(**self.kw)
        T = _table_from_prov(m["prov_transfer_keys"], m["prov_transfer_vals"])
        E = (_table_from_prov(m["prov_extern_keys"], m["prov_extern_vals"])
             if bool(m["prov_extern_present"]) else None)
        got = fp.band_digests(T, E, _meta_from_prov(m))
        for k in ("fp_transfer", "fp_extern", "fp_meta", "fp_band"):
            self.assertEqual(str(m[k]), got[k])
        # absent Extern: present False, empty tables, the absent digest
        m2 = fp.so_provenance_members(**dict(self.kw, extern=None))
        self.assertFalse(bool(m2["prov_extern_present"]))
        self.assertEqual(m2["prov_extern_keys"].shape, (0, 5))
        self.assertEqual(m2["prov_extern_vals"].shape, (0,))
        self.assertEqual(str(m2["fp_extern"]), fp.table_digest("extern", None))
        self.assertNotEqual(str(m2["fp_band"]), str(m["fp_band"]))

    def test_state_tagged_union(self):
        s = fp.so_state(beta=2.0, nmat=8, mu=0.3, density=0.75, calc_mu=True, target=0.8)
        self.assertAlmostEqual(s.mu_tol, 0.05, delta=1e-15)
        s = fp.so_state(beta=2.0, nmat=8, mu=0.3, density=0.75, calc_mu=False, target=0.3)
        self.assertEqual(s.mu_tol, 0.0)
        self.assertEqual(s.target, 0.3)
        with self.assertRaisesRegex(ValueError, "state_mu_tol"):
            fp.so_provenance_members(**dict(self.kw, state=s._replace(mu_tol=1e-9)))

    def test_nonfinite_refused(self):
        for field in ("mu", "density", "target", "beta"):
            with self.subTest(field=field):
                bad = self.kw["state"]._replace(**{field: float("nan")})
                with self.assertRaisesRegex(ValueError, "state_" + field):
                    fp.so_provenance_members(**dict(self.kw, state=bad))
        with self.assertRaisesRegex(ValueError, "prov_coeff_extern"):
            fp.so_provenance_members(**dict(self.kw, meta=self.kw["meta"]._replace(
                coeff_extern=float("inf"))))
        with self.assertRaisesRegex(ValueError, "producer_stage"):
            fp.so_provenance_members(**dict(self.kw, producer_stage="rpa"))

    def test_long_producer_version_truncated(self):
        with mock.patch.object(hwave, "__version__", "x" * 400):
            with self.assertLogs("hwave.solver.fingerprint", level="INFO") as cm:
                v = str(fp.so_provenance_members(**self.kw)["producer_version"])
        self.assertEqual(len(v), 256)
        self.assertTrue(v.endswith("…"))
        self.assertTrue(any("producer_version truncated to 256 characters" in r
                            for r in cm.output))
        # read at call time, stored at its exact length
        with mock.patch.object(hwave, "__version__", "9.9.9-test"):
            v = str(fp.so_provenance_members(**self.kw)["producer_version"])
        self.assertTrue(v.startswith("9.9.9-test"), v)

    def test_git_describe_cached_and_never_fails(self):
        fp._git_describe.cache_clear()
        self.addCleanup(fp._git_describe.cache_clear)
        with mock.patch.object(fp.subprocess, "run", side_effect=OSError("no git")) as run:
            with mock.patch.object(hwave, "__version__", "1.2.3"):
                self.assertEqual(fp._producer_version(), "1.2.3")
                self.assertEqual(fp._producer_version(), "1.2.3")
            self.assertEqual(run.call_count, 1)
        fp._git_describe.cache_clear()
        ok = mock.Mock(returncode=0, stdout="v2.0.0-5-gabc\n")
        with mock.patch.object(fp.subprocess, "run", return_value=ok) as run:
            with mock.patch.object(hwave, "__version__", "2.0.0"):
                self.assertEqual(fp._producer_version(), "2.0.0 (v2.0.0-5-gabc)")
            self.assertEqual(run.call_args.kwargs["timeout"], 2.0)


class TestProducersWrite(_Dir):
    def setUp(self):
        super().setUp()
        self.inp_soc = _so_input(self.new_dir())
        p = self.inp_soc["param"]
        self.nvol = int(np.prod(self.inp_soc["cell"]))
        self.ncond = p["filling"] * 2 * self.nvol                 # Nstate = nd * nvol

    def _check(self, npz, stage, inp):
        self.assertEqual(int(npz["fingerprint_version"]), 1)
        self.assertEqual(str(npz["producer_stage"]), stage)
        T, E, meta = _tables_and_meta_from_input(inp)
        self.assertEqual(str(npz["fp_band"]), fp.band_digests(T, E, meta)["fp_band"])
        self.assertEqual(float(npz["state_beta"]), 1.0 / inp["param"]["T"])
        self.assertEqual(int(npz["state_nmat"]), inp["param"]["Nmat"])
        for k, (kind, ndim) in _SCHEMA.items():
            self.assertEqual((npz[k].dtype.kind, npz[k].ndim), (kind, ndim), k)

    def test_rpa_so_calc_mu(self):
        out = _run_rpa_so(self.inp_soc)
        for k in ("chi0q", "chiq"):
            self._check(out[k], "rpa_so", self.inp_soc)
            self.assertEqual(int(out[k]["nmat"]), int(out[k]["state_nmat"]))
        z = out["chi0q"]
        self.assertTrue(bool(z["state_calc_mu"]))
        self.assertAlmostEqual(float(z["state_target"]), self.ncond / self.nvol, delta=1e-15)
        self.assertLessEqual(abs(float(z["state_density"]) - self.ncond / self.nvol),
                             float(z["state_mu_tol"]) + 1e-12)
        self.assertLess(float(z["state_mu_tol"]), 1e-10)
        # the realized mu reproduces the density in the closed form
        self.assertAlmostEqual(
            _closed_form_fermi_sum(self.inp_soc, float(z["state_mu"])) / self.nvol,
            float(z["state_density"]), delta=1e-12)
        # one build for both archives
        for k in _SCHEMA:
            np.testing.assert_array_equal(out["chi0q"][k], out["chiq"][k])

    def test_rpa_so_fixed_mu_density(self):
        inp = _so_input(self.new_dir(), mu=0.1)
        out = _run_rpa_so(inp)
        z = out["chi0q"]
        self._check(z, "rpa_so", inp)
        self.assertFalse(bool(z["state_calc_mu"]))
        self.assertEqual(float(z["state_mu_tol"]), 0.0)
        self.assertEqual(float(z["state_target"]), 0.1)
        self.assertEqual(float(z["state_mu"]), 0.1)
        self.assertAlmostEqual(float(z["state_density"]),
                               _closed_form_fermi_sum(inp, mu=0.1) / self.nvol, delta=1e-12)

    def test_rpa_so_density_counts_every_spin_orbital(self):
        """A spin-free (no field, no SOC) and a spin-diagonal band count both
        spins: the closure density equals the target in electrons per cell."""
        for coeff, expect in ((0.0, "spin-free"), (1.0, "spin-diag")):
            with self.subTest(band=expect):
                d = self.new_dir()
                inter = Z._so_2orb_fixture(d, "extern")
                param = dict(Z._E2_PARAM, coeff_extern=coeff)
                inp = {"dir": d, "inter": inter, "param": param,
                       "cell": (4, 4, 1), "sub_shape": (1, 1, 1)}
                import hwave.solver.rpa as rpa_mod
                modes = []
                orig = rpa_mod.RPA._so_bare_state

                def spy(self_, mu, _orig=orig):
                    modes.append(self_.spin_mode)
                    return _orig(self_, mu)
                with mock.patch.object(rpa_mod.RPA, "_so_bare_state", spy):
                    z = _run_rpa_so(inp)["chi0q"]
                self.assertEqual(modes, [expect])
                target = param["filling"] * 4                      # nd = 4 per cell
                self.assertAlmostEqual(float(z["state_target"]), target, delta=1e-14)
                self.assertLess(abs(float(z["state_density"]) - target), 1e-10)
                self.assertAlmostEqual(float(z["state_mu_tol"]),
                                       abs(float(z["state_density"]) - target), delta=1e-15)

    def test_rpa_sublattice_hashes_prefold(self):
        inp = _so_input(self.new_dir(), sub_shape=(2, 1, 1), cell=(4, 1, 1))
        a = _run_rpa_so(inp)["chi0q"]
        T, E, meta = _tables_and_meta_from_input(inp)
        self.assertEqual(str(a["fp_transfer"]), fp.band_digests(T, E, meta)["fp_transfer"])
        self.assertEqual(str(a["fp_extern"]), fp.band_digests(T, E, meta)["fp_extern"])
        self.assertEqual(str(a["fp_band"]), fp.band_digests(T, E, meta)["fp_band"])
        self.assertEqual(int(a["prov_norb"]), 2)
        self.assertEqual(tuple(a["prov_sub_shape"]), (2, 1, 1))
        # the sublattice run folds (nvol 2): density is per (super)cell of the run
        self.assertLess(float(a["state_mu_tol"]), 1e-10)

    def test_flex_so_bindings_agree(self):
        out = _run_flex_so(self.inp_soc)
        d = fp.green_digest(out["green"]["green"], "b_f_k_o_o_v1")
        for k in ("chiq", "green", "sigma", "chi0q"):
            with self.subTest(archive=k):
                self._check(out[k], "flex_so_onsite", self.inp_soc)
                self.assertEqual(str(out[k]["fp_green"]), d)
                for m in _SCHEMA:
                    np.testing.assert_array_equal(out[k][m], out["green"][m])
        g = out["green"]
        self.assertEqual(float(g["beta"]), float(g["state_beta"]))
        self.assertTrue(bool(g["state_calc_mu"]))
        self.assertLessEqual(float(g["state_mu_tol"]), 1e-10)
        self.assertAlmostEqual(float(g["state_target"]), self.ncond / self.nvol, delta=1e-15)

    def test_flex_so_fixed_mu(self):
        inp = self.inp_soc
        param = dict(Z._F1_PARAM, IterationMax=2, coeff_extern=1.0, mu=0.1)
        param.pop("filling")
        out = _run(dict(inp, param=param), "FLEX", {"green": "green.npz", "sigma": "sigma.npz"},
                   {"calc_scheme": "general"})
        g = out["green"]
        self.assertFalse(bool(g["state_calc_mu"]))
        self.assertEqual(float(g["state_mu"]), 0.1)
        self.assertEqual(float(g["state_target"]), 0.1)
        self.assertEqual(float(g["state_mu_tol"]), 0.0)
        self.assertTrue(0.0 < float(g["state_density"]) < 2.0)

    def test_built_once_per_run(self):
        real = fp.so_provenance_members
        for run, n_archives in ((_run_rpa_so, 2), (_run_flex_so, 4)):
            with self.subTest(run=run.__name__):
                with mock.patch.object(fp, "so_provenance_members", side_effect=real) as m:
                    out = run(_so_input(self.new_dir()))
                self.assertEqual(m.call_count, 1)
                self.assertEqual(len(out), n_archives)

    def test_mixed_case_keys(self):
        base = _run_rpa_so(_so_input(self.new_dir()))["chi0q"]
        mixed = _run_rpa_so(_so_input(self.new_dir(), extern_key="EXTERN",
                                      coeff_key="Coeff_Extern"))["chi0q"]
        self.assertTrue(bool(mixed["prov_extern_present"]))
        self.assertEqual(float(mixed["prov_coeff_extern"]), 1.0)
        self.assertEqual(str(mixed["fp_band"]), str(base["fp_band"]))

    def test_non_so_archives_bit_identical(self):
        """Member set, dtype, shape and values of the non-spin-orbital
        archives equal the pre-change snapshot (exactly on the platform it
        was made on; elsewhere floats within 1e-12 of their scale, since a
        different LAPACK may round differently)."""
        exact = ref_make.platform_tag() == _ref_platform()
        for name, run in ref_make.RUNS.items():
            ref = _reference_outputs(name)
            self.assertTrue(ref)
            out = os.path.join(self.new_dir(), "out")
            os.makedirs(out)
            run(out)
            self.assertEqual(sorted(os.path.basename(f) for f in
                                    glob.glob(os.path.join(out, "*.npz"))), sorted(ref))
            for fname, r in ref.items():
                with np.load(os.path.join(out, fname)) as new:
                    self.assertEqual(set(new.files), set(r.files), (name, fname))
                    for k in r.files:
                        a, b = new[k], r[k]
                        with self.subTest(run=name, archive=fname, member=k):
                            self.assertEqual((a.dtype, a.shape), (b.dtype, b.shape))
                            if exact or a.dtype.kind not in "fc":
                                np.testing.assert_array_equal(a, b)
                            else:
                                scale = max(1.0, float(np.max(np.abs(b)))) if b.size else 1.0
                                np.testing.assert_allclose(a, b, rtol=0, atol=1e-12 * scale)


if __name__ == "__main__":
    unittest.main()
