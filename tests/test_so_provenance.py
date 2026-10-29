"""Version-1 provenance of the spin-orbital RPA and FLEX archives (#88, PR 0).

A spin-orbital RPA run writes the fingerprint and state members into chi0q
and chiq; a spin-orbital Hartree-Fock FLEX run writes them into chiq, green,
sigma and chi0q, with one Green digest (that of the G in green.npz) in all
four. What this module pins:

* :class:`TestMembers` -- ``fingerprint.so_provenance_members``: the member
  schema (dtype, shape), the digests recomputable from the stored ``prov_*``
  tables, the ``calc_mu`` tagged union of the state, the refusal of a
  non-finite member, ``producer_version`` (truncation, git describe cached,
  never failing, and taken only from the package's own repository).
* :class:`TestProducersWrite` -- the RPA and FLEX producers through
  ``hwave.qlms.run``: the band digest equals the one computed from the input
  files, the state is the realized one (closure density in electrons per
  cell over all spin-orbitals, a fixed-mu density against a closed-form
  Fermi sum), the PRE-FOLD tables are hashed under a sublattice, the members
  are built once per run, mixed-case keys, no provenance when trans_mod
  changes the band, and the non-spin-orbital archives against the snapshot
  in ``tests/so_provenance_ref`` (made by ``tests/so_provenance_ref/make.py``
  on the tree before this change): member set, dtype, shape and non-float
  values exact everywhere; float values exact on the recorded platform,
  within 1e-12 of scale elsewhere (CI cannot reproduce the recorded
  platform's bits).
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

    def _git(self, top):
        """A fake ``subprocess.run``: ``rev-parse`` answers ``top``,
        ``describe`` answers a fixed description."""
        def run(cmd, **kw):
            self.assertEqual(kw["timeout"], 2.0)
            if cmd[1] == "rev-parse":
                return mock.Mock(returncode=0, stdout=top + "\n")
            return mock.Mock(returncode=0, stdout="v2.0.0-5-gabc\n")
        return run

    def test_git_describe_only_from_own_repository(self):
        pkg = os.path.dirname(os.path.dirname(os.path.abspath(fp.__file__)))
        own = os.path.dirname(os.path.dirname(pkg))                 # <top>/src/hwave
        self.addCleanup(fp._git_describe.cache_clear)
        for top, expect in ((own, "2.0.0 (v2.0.0-5-gabc)"),       # the package's repo
                            (os.path.dirname(own), "2.0.0"),       # an enclosing repo
                            (self.d, "2.0.0")):                    # unrelated repo
            with self.subTest(top=top):
                fp._git_describe.cache_clear()
                with mock.patch.object(fp.subprocess, "run", side_effect=self._git(top)):
                    with mock.patch.object(hwave, "__version__", "2.0.0"):
                        self.assertEqual(fp._producer_version(), expect)
        self.assertTrue(os.path.isfile(os.path.join(own, "pyproject.toml")))


class TestProducersWrite(_Dir):
    def setUp(self):
        super().setUp()
        self.inp_soc = _so_input(self.new_dir())
        p = self.inp_soc["param"]
        self.ncell = int(np.prod(self.inp_soc["cell"]))           # ORIGINAL cells
        self.ncond = p["filling"] * 2 * self.ncell                # Nstate = nd * ncell

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
        self.assertAlmostEqual(float(z["state_target"]), self.ncond / self.ncell, delta=1e-15)
        self.assertLessEqual(abs(float(z["state_density"]) - self.ncond / self.ncell),
                             float(z["state_mu_tol"]) + 1e-12)
        self.assertLess(float(z["state_mu_tol"]), 1e-10)
        # the realized mu reproduces the density in the closed form
        self.assertAlmostEqual(
            _closed_form_fermi_sum(self.inp_soc, float(z["state_mu"])) / self.ncell,
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
                               _closed_form_fermi_sum(inp, mu=0.1) / self.ncell, delta=1e-12)

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
        self.assertLess(float(a["state_mu_tol"]), 1e-10)

    def test_rpa_sublattice_state_per_original_cell(self):
        """Under a sublattice (2 supercells of 2 cells) the state is per
        ORIGINAL unit cell: the target is filling * nd (nd = 2 spin-orbitals
        per cell), not per folded supercell, and the density is the
        closed-form count of the unfolded band over the 4 original cells."""
        for mu in (None, 0.1):
            with self.subTest(mu=mu):
                inp = _so_input(self.new_dir(), sub_shape=(2, 1, 1), cell=(4, 1, 1), mu=mu)
                a = _run_rpa_so(inp)["chi0q"]
                self.assertEqual(a["chi0q"].shape[1], 2)            # folded: nvol = 2
                ref = _closed_form_fermi_sum(inp, float(a["state_mu"])) / 4
                self.assertAlmostEqual(float(a["state_density"]), ref, delta=1e-12)
                if mu is None:
                    self.assertAlmostEqual(float(a["state_target"]),
                                           inp["param"]["filling"] * 2, delta=1e-14)
                    self.assertLess(float(a["state_mu_tol"]), 1e-10)
                else:
                    self.assertEqual(float(a["state_target"]), 0.1)

    def test_rpa_trans_mod_writes_no_provenance(self):
        """A one-body band changed by trans_mod is not covered by the
        version-1 band digest: the archives keep the spin_orbital marker but
        carry no provenance member, and one INFO line says why."""
        import tests.test_rpa_so_trans_mod as TM
        tc = TM.TestRPASOTransModMultiOrbital()
        tc._tmpdir = self.new_dir()
        tm_path = tc._write_fixture()
        inter = {"Geometry": "geom_so_2orb.dat", "Transfer": "transfer_so_2orb.dat",
                 "CoulombIntra": "coulombintra_2orb.dat"}
        prov = set(fp.SO_PROVENANCE_MEMBERS)
        for trans_mod in (None, tm_path):
            with self.subTest(trans_mod=trans_mod is not None):
                out = os.path.join(self.new_dir(), "out")
                with _quiet():
                    solver, green = TM._run({"enable_spin_orbital": True}, inter, out,
                                            trans_mod_file=trans_mod)
                files = {"path_to_output": out, "chi0q": "chi0q.npz", "chiq": "chiq.npz"}
                with self.assertLogs("hwave.solver.rpa", level="INFO") as cm:
                    solver.save_results(files, green)
                note = [r for r in cm.output if "provenance not written" in r]
                for f in ("chi0q.npz", "chiq.npz"):
                    with np.load(os.path.join(out, f)) as z:
                        self.assertTrue(bool(z["spin_orbital"]))
                        if trans_mod is None:                       # anti-vacuity
                            self.assertLessEqual(prov, set(z.files))
                        else:
                            self.assertFalse(prov & set(z.files))
                if trans_mod is None:
                    self.assertEqual(note, [])
                else:
                    self.assertEqual(len(note), 1)
                    self.assertIn("trans_mod/green_init", note[0])

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
        self.assertAlmostEqual(float(g["state_target"]), self.ncond / self.ncell, delta=1e-15)

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
        """The non-spin-orbital archives against the pre-change snapshot:
        member set, dtype, shape and non-float values exact everywhere; float
        values exact on the recorded platform, within 1e-12 of scale
        elsewhere (CI cannot reproduce the recorded platform's bits).

        A failure here means a non-spin-orbital archive changed. If that
        change is intended (a deliberate new or changed non-spin-orbital
        member), the snapshot must be regenerated, never loosened: run, from
        the repository root on the tree WITH the change,

            PYTHONPATH=$PWD/src:$PWD python3 -B tests/so_provenance_ref/make.py

        which rewrites ``tests/so_provenance_ref/<run>/*.npz`` and
        ``platform.txt``; commit them together with the change and say in
        the commit message which member changed and why."""
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


def _traced_peak(fn):
    """Peak bytes traced by ``tracemalloc`` (NumPy reports its buffers to
    it) while ``fn()`` runs, above the level at its start."""
    import tracemalloc
    tracemalloc.start()
    try:
        base = tracemalloc.get_traced_memory()[0]
        fn()
        return tracemalloc.get_traced_memory()[1] - base
    finally:
        tracemalloc.stop()


class TestRpaGreenDigest(unittest.TestCase):
    """``RPA._so_green_digest`` (the ``fp_green`` of a spin-orbital RPA run)
    streams the canonical G: its digest is byte-identical to the
    whole-array definition (deflated G plus tail, a spin-free/spin-diagonal
    band expanded to spin blocks) for every band kind, and no Green-sized
    array is allocated while it is computed."""

    @staticmethod
    def _whole_array(spin_mode, green, tail):
        # the whole-array definition fp_green had before it was streamed
        g = green + tail
        nb, nmat, nvol, n, _ = g.shape
        if spin_mode == "spinful":
            return g
        full = np.zeros((1, nmat, nvol, 2 * n, 2 * n), dtype=np.complex128)
        for s in (0, 1):
            full[0, :, :, s * n:(s + 1) * n, s * n:(s + 1) * n] = g[s if nb == 2 else 0]
        return full

    @staticmethod
    def _cases(nmat=8):
        rng = np.random.default_rng(7)

        def c(*shape):
            return rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
        # (spin_mode, nblock, n): spin-free (1 block, norb), spin-diagonal
        # (2 blocks, norb), spinful (1 block, 2*norb)
        for mode, nb, n in (("spin-free", 1, 3), ("spin-diag", 2, 3), ("spinful", 1, 6)):
            green = c(nb, nmat, 5, n, n)
            # the tail is constant over frequency (V V^dagger coeff_tail beta/2)
            tail = np.broadcast_to(c(nb, 1, 5, n, n), green.shape).copy()
            yield mode, green, tail

    def _digest(self, mode, green, tail):
        from hwave.solver.rpa import RPA
        from types import SimpleNamespace
        return RPA._so_green_digest(SimpleNamespace(spin_mode=mode), green, tail)

    def test_streaming_equals_whole_array_digest(self):
        for mode, green, tail in self._cases():
            ref = fp.green_digest(self._whole_array(mode, green, tail), "b_f_k_o_o_v1")
            for chunk in (fp.CHUNK_BYTES, 16, 5 * 36 * 16, 3 * 5 * 36 * 16 + 16):
                with self.subTest(mode=mode, chunk=chunk), \
                        mock.patch.object(fp, "CHUNK_BYTES", chunk):
                    self.assertEqual(self._digest(mode, green, tail), ref)
            with self.subTest(mode=mode, tail="zero"):
                z = np.zeros_like(green)
                self.assertEqual(self._digest(mode, green, z), fp.green_digest(
                    self._whole_array(mode, green, z), "b_f_k_o_o_v1"))

    def test_no_green_sized_allocation(self):
        """Peak traced allocation while digesting stays below half of the
        deflated G (one frequency slice per block); the whole-array
        definition, measured the same way, exceeds it (non-vacuity)."""
        for mode, green, tail in self._cases(nmat=64):
            bound = green.nbytes // 2
            nd = 6
            with self.subTest(mode=mode), \
                    mock.patch.object(fp, "CHUNK_BYTES", 5 * nd * nd * 16):
                peak = _traced_peak(lambda: self._digest(mode, green, tail))
                peak_whole = _traced_peak(lambda: fp.green_digest(
                    self._whole_array(mode, green, tail), "b_f_k_o_o_v1"))
                self.assertLess(peak, bound)
                self.assertGreater(peak_whole, bound)


# ------------------------------------------------- consumers (hwave_sc, #88)
_SC_ELI = {"solver_mode": "eigenvalue", "num_eigenvalues": 2}


def _soc_input(d, extern=False):
    """The 0-E1 spin-mixing fixture (complex SOC 0.35+0.15i, one physical
    orbital) for the hwave_sc consumers; ``extern = False`` drops the Extern
    key (the Transfer file is the same either way)."""
    inp = _so_input(d)
    if not extern:
        inp["inter"] = {k: v for k, v in inp["inter"].items() if k.lower() != "extern"}
    return inp


def _zeeman_only_input(d):
    """The 0-E2 fixture (two orbitals, Extern Zeeman field, no SOC): the RPA
    reduces it to the spin-diagonal form."""
    inter = Z._so_2orb_fixture(d, "extern")
    return {"dir": d, "inter": inter, "param": dict(Z._E2_PARAM),
            "cell": (4, 4, 1), "sub_shape": (1, 1, 1)}


def _with(param, over):
    """``param`` with the overrides ``coeff_extern``, ``mu`` (drops
    ``filling``), ``nmat``, ``T`` and ``filling``."""
    param = dict(param)
    for k, v in over.items():
        if k == "mu":
            param.pop("filling", None)
            param["mu"] = v
        elif k == "nmat":
            param["Nmat"] = v
        else:
            param[k] = v
    return param


def _rpa_for_sc(inp, **over):
    """The spin-orbital RPA of ``inp`` (with the parameter overrides);
    returns ``{"dir": its output directory, "path": its chi0q.npz}``."""
    _run(dict(inp, param=_with(inp["param"], over)), "RPA",
         {"chi0q": "chi0q.npz", "chiq": "chiq.npz"}, _RPA_EXTRA)
    out = os.path.join(inp["dir"], "output_RPA")
    return {"dir": out, "path": os.path.join(out, "chi0q.npz")}


def _flex_for_sc(inp, **over):
    """The spin-orbital Hartree-Fock FLEX of ``inp`` (3 maps), writing chiq,
    green, sigma and chi0q; returns the output directory and the paths."""
    param = _with(dict(Z._F1_PARAM, IterationMax=3, coeff_extern=1.0,
                       CellShape=list(inp["cell"])), over)
    _run(dict(inp, param=param), "FLEX",
         {"chiq": "chiq.npz", "green": "green.npz", "sigma": "sigma.npz",
          "chi0q": "chi0q.npz"}, {"calc_scheme": "general"})
    out = os.path.join(inp["dir"], "output_FLEX")
    return {"dir": out, **{k: os.path.join(out, k + ".npz")
                           for k in ("chiq", "green", "sigma", "chi0q")}}


def _resave(path, drop=(), **replace):
    with np.load(path) as z:
        items = {k: np.array(z[k]) for k in z.files if k not in drop}
    items.update(replace)
    np.savez(path, **items)


def _drop_member(path, k):
    _resave(path, drop=(k,))


def _set_member(path, k, v):
    _resave(path, **{k: v})


def _strip_fingerprint_members(path):
    """An archive as written before version 1: no provenance member (the
    spin_orbital marker stays)."""
    _resave(path, drop=fp.SO_PROVENANCE_MEMBERS)


def _overwrite_green_entry(green_path, value):
    with np.load(green_path) as z:
        g = np.array(z["green"])
    g[0, 0, 0, 0, 0] += value
    _set_member(green_path, "green", g)


def _edit_w90(path, R, ab, dv):
    """Add ``dv`` to the row ``(R, a, b)`` (0-based orbitals) of a W90 file
    (appended when absent)."""
    dim, _nr, _deg, rows = Z._read_w90_rows(path)
    key = tuple(R) + (ab[0] + 1, ab[1] + 1)
    for i, r in enumerate(rows):
        if r[:5] == key:
            rows[i] = r[:5] + (r[5] + dv, r[6])
            break
    else:
        rows.append(key + (dv, 0.0))
    Z._write_w90_file(path, dim, rows)


class _ScOut:
    def __init__(self, npz, kernel_green):
        self.npz = npz
        self.kernel_green = kernel_green


def _sc(mode, inp, src=None, frequency="static", transfer_edit=None,
        extern_edit=None, quiet=True, **over):
    """A spin-orbital hwave_sc run on ``inp`` (parameter overrides
    ``over``) in ``chi0q_mode = mode`` reading the archives of ``src`` (the
    output directory of :func:`_rpa_for_sc` / :func:`_flex_for_sc`; default:
    that of ``inp``). ``transfer_edit`` / ``extern_edit`` ``(R, (a, b), dv)``
    change the run's input files first. Returns the gap archive and the G
    the pair bubble was built from."""
    import hwave.sc as sc
    d = inp["dir"]
    if transfer_edit is not None:
        _edit_w90(os.path.join(d, inp["inter"]["Transfer"]), *transfer_edit)
    if extern_edit is not None:
        _edit_w90(os.path.join(d, inp["inter"]["Extern"]), *extern_edit)
    if src is None:
        src = os.path.join(d, "output_RPA" if mode == "load" else "output_FLEX")
    out = src if mode == "load" else os.path.join(d, "output_sc_" + mode)
    eli = dict(_SC_ELI, chi0q_mode=mode)
    if frequency == "dynamic":
        eli.update(frequency="dynamic", spectral_shift="auto")
    run_inp = {
        "mode": {"mode": "RPA", "enable_spin_orbital": True,
                 "param": _with(inp["param"], over)},
        "file": {"input": {"path_to_input": "", "path_to_flex_output": src,
                           "interaction": dict(inp["inter"], path_to_input=d)},
                 "output": {"path_to_output": out, "chi0q": "chi0q"}},
        "eliashberg": eli}
    seen = {}
    real_g2 = sc._calc_g2

    def spy_g2(green, *a, **kw):
        seen["green"] = np.array(green, copy=True)
        return real_g2(green, *a, **kw)

    import hwave.solver.eliashberg_dynamic_so as dso
    real_g2d = dso.calc_g2_dynamic

    def spy_g2d(green, *a, **kw):
        seen["green"] = np.array(green, copy=True)
        return real_g2d(green, *a, **kw)

    ctx = _quiet() if quiet else Z._quiet_info()
    with ctx, mock.patch.object(sc, "_calc_g2", side_effect=spy_g2), \
            mock.patch.object(dso, "calc_g2_dynamic", side_effect=spy_g2d):
        sc.calc_eliashberg(run_inp)
    name = "gap.npz" if frequency == "static" else "gap_dynamic.npz"
    with np.load(os.path.join(out, name)) as z:
        npz = {k: np.array(z[k]) for k in z.files}
    return _ScOut(npz, seen.get("green"))


def _write_members(d, members):
    path = os.path.join(tempfile.mkdtemp(dir=d), "a.npz")
    np.savez(path, **members)
    return path


def _write_raw(d, lines):
    path = os.path.join(tempfile.mkdtemp(dir=d), "raw.dat")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return path


class TestConsumerChecks(_Dir):
    """``fingerprint.check_so_archive`` through the hwave_sc consumers: the
    structural, band, state and green refusals, the legacy policy and the
    exact state-point rule."""

    def setUp(self):
        super().setUp()
        self.inp_soc = _soc_input(self.new_dir())
        self.inp_soc_extern = _soc_input(self.new_dir(), extern=True)
        self.inp_zeeman_only = _zeeman_only_input(self.new_dir())
        T, E, meta = _tables_and_meta_from_input(self.inp_soc_extern)
        self.band = fp.band_inputs(T, E, meta)
        self.kw = dict(producer_stage="rpa_so", producer_label="test producer",
                       transfer=T, extern=E, meta=meta, fp_green="0" * 64,
                       state=fp.SOState(10.0, 64, 0.1, 2.0, True, 2.0, 0.0))
        nan = float("nan")
        self.run_state = fp.SOState(10.0, 64, nan, nan, True, 2.0, nan)

    def test_run_inputs_mixed_case_keys(self):
        """``sc._so_run_inputs`` reads ``CellShape``/``SubShape`` (and
        ``coeff_extern``) case-insensitively: mixed-case keys give the same
        band metadata and digests as the canonical spelling, and an omitted
        ``SubShape`` defaults to a lowercase ``cellshape``."""
        from hwave import sc
        T, E, _ = _tables_and_meta_from_input(self.inp_soc_extern)

        def run(param):
            return sc._so_run_inputs({"mode": {"param": param}}, T, E, 1, 10.0, 64, 0.5)
        canon = run({"CellShape": [4, 2], "SubShape": [1, 1, 1], "coeff_extern": 0.5})
        mixed = run({"cellshape": [4, 2], "SUBSHAPE": [1, 1, 1], "Coeff_Extern": 0.5})
        self.assertEqual(canon.band.meta.cell_shape, (4, 2, 1))
        self.assertEqual(mixed.band.meta, canon.band.meta)
        self.assertEqual(mixed.band.digests, canon.band.digests)
        dflt = run({"cellShape": [4, 2]})
        self.assertEqual(dflt.band.meta.sub_shape, (4, 2, 1))

    def test_field_sweep_cannot_borrow(self):
        _rpa_for_sc(self.inp_soc, coeff_extern=0.0)
        with self.assertRaisesRegex(ValueError, r"coeff_extern: file 0\.0, run 0\.1.*"
                                                r"identical state point \(Nmat excepted in "
                                                r"chi0q_mode = \"load\"\)"):
            _sc("load", self.inp_soc, coeff_extern=0.1)

    def test_state_message_lists_all_fields(self):
        _flex_for_sc(self.inp_soc, filling=0.5, T=0.1)
        with self.assertRaisesRegex(ValueError, r"state differs: .*density .* \(file\) "
                                                r"vs target .*; beta "):
            _sc("flex", self.inp_soc, filling=0.4, T=0.2)

    def test_beta_mismatch(self):
        _flex_for_sc(self.inp_soc, T=0.1)
        with self.assertRaisesRegex(ValueError, r"state differs: beta"):
            _sc("flex", self.inp_soc, T=0.1 * (1 + 1e-9))

    def test_green_digest_mismatch(self):
        out = _flex_for_sc(self.inp_soc)
        _overwrite_green_entry(out["green"], value=1e-3)
        with self.assertRaisesRegex(ValueError, r"fp_green mismatch"):
            _sc("flex", self.inp_soc)

    def test_unversioned_archive_rules(self):
        _strip_fingerprint_members(_rpa_for_sc(self.inp_soc)["path"])
        with self.assertLogs("qlms", "WARNING") as cm:
            _sc("load", self.inp_soc, quiet=False)
        self.assertTrue(any("has no fingerprint_version" in m for m in cm.output))
        _strip_fingerprint_members(_rpa_for_sc(self.inp_soc_extern)["path"])
        with self.assertRaisesRegex(ValueError, r"has an Extern field.*Regenerate"):
            _sc("load", self.inp_soc_extern)

    def test_missing_member_and_bad_stage_and_dup_field(self):
        for mutate, rx in ((lambda p: _drop_member(p, "fp_band"),
                            r"malformed or truncated: member fp_band is missing"),
                           (lambda p: _set_member(p, "producer_stage", "rpa"),
                            r"malformed.*producer_stage"),
                           (lambda p: _set_member(p, "state_beta", 1.0), r"beta")):
            with self.subTest(rx=rx):
                mutate(_rpa_for_sc(self.inp_soc)["path"])
                with self.assertRaisesRegex(ValueError, rx):
                    _sc("load", self.inp_soc)

    def test_duplicated_nmat_refused(self):
        p = _rpa_for_sc(self.inp_soc)["path"]
        with np.load(p) as z:
            self.assertIn("nmat", z.files)                       # anti-vacuity
        _set_member(p, "state_nmat", np.int64(17))
        with self.assertRaisesRegex(ValueError, r"malformed.*nmat differs from state_nmat"):
            _sc("load", self.inp_soc)

    def test_transfer_row_diff(self):
        _rpa_for_sc(self.inp_soc)
        with self.assertRaisesRegex(ValueError, r"Transfer: file \d+ rows, run \d+ rows; "
                                                r"first differing row \(R=\(1, 0, 0\), "
                                                r"a=1, b=2\)"):
            _sc("load", self.inp_soc, transfer_edit=((1, 0, 0), (0, 1), 0.01))

    def test_extern_presence_and_row_diff(self):
        src = _rpa_for_sc(self.inp_soc)["dir"]
        with self.assertRaisesRegex(ValueError, r"Extern: file absent, run present"):
            _sc("load", self.inp_soc_extern, src=src)
        _rpa_for_sc(self.inp_soc_extern)
        with self.assertRaisesRegex(ValueError, r"Extern: .*first differing row "
                                                r"\(R=\(0, 0, 0\), a=1, b=1\)"):
            _sc("load", self.inp_soc_extern, extern_edit=((0, 0, 0), (0, 0), 0.2))

    def test_tampered_provenance_and_inflated_tolerance_refused(self):
        for mutate, rx in ((lambda p: _set_member(p, "prov_coeff_extern", 0.5),
                            r"malformed.*fp_meta"),
                           (lambda p: _set_member(p, "state_mu_tol", 1.0),
                            r"malformed.*state_mu_tol")):
            with self.subTest(rx=rx):
                mutate(_rpa_for_sc(self.inp_soc)["path"])
                with self.assertRaisesRegex(ValueError, rx):
                    _sc("load", self.inp_soc)

    def test_tagged_union_round_trip(self):
        for calc_mu, state in ((True, fp.SOState(10.0, 64, 0.1, 2.0 + 1e-11, True, 2.0, 1e-11)),
                               (False, fp.SOState(10.0, 64, 0.1, 1.97, False, 0.1, 0.0))):
            with self.subTest(calc_mu=calc_mu):
                path = _write_members(self.d, fp.so_provenance_members(
                    **dict(self.kw, state=state)))
                run = self.run_state._replace(
                    calc_mu=calc_mu, target=state.target,
                    mu=state.mu if not calc_mu else float("nan"))
                info = fp.check_so_archive(np.load(path), path, run_band=self.band,
                                           run_state=run, run_has_extern=False,
                                           nmat_rule="exact")
                self.assertFalse(info.legacy)
                self.assertEqual(info.producer_state.calc_mu, calc_mu)
                self.assertEqual(info.producer_state, state)
                self.assertEqual(info.producer_label, "test producer")
                # a fixed mu that differs is refused, in energy units
                if not calc_mu:
                    with self.assertRaisesRegex(ValueError, r"fixed mu 0\.1 \(file\) "
                                                            r"vs 0\.2 \(run\)"):
                        fp.check_so_archive(np.load(path), path, run_band=self.band,
                                            run_state=run._replace(mu=0.2),
                                            run_has_extern=False, nmat_rule="exact")

    def test_nmat_rule(self):
        path = _write_members(self.d, fp.so_provenance_members(**self.kw))
        run = self.run_state._replace(nmat=128)
        with self.assertRaisesRegex(ValueError, r"nmat 64 \(file\) vs 128 \(run\)"):
            fp.check_so_archive(np.load(path), path, run_band=self.band, run_state=run,
                                run_has_extern=False, nmat_rule="exact")
        fp.check_so_archive(np.load(path), path, run_band=self.band, run_state=run,
                            run_has_extern=False, nmat_rule="load")
        with self.assertRaisesRegex(ValueError, r"nmat_rule"):
            fp.check_so_archive(np.load(path), path, run_band=self.band, run_state=run,
                                run_has_extern=False, nmat_rule="loose")

    def test_legacy_checked_first_and_unknown_label(self):
        """No fingerprint_version: the legacy policy, even with other members
        missing or broken; a malformed versioned archive without a readable
        label names <unknown producer>."""
        path = _write_members(self.d, {"spin_orbital": np.bool_(True),
                                       "state_beta": np.str_("broken")})
        with self.assertLogs("qlms", "WARNING"):
            info = fp.check_so_archive(np.load(path), path, run_band=self.band,
                                       run_state=self.run_state, run_has_extern=False,
                                       nmat_rule="exact")
        self.assertEqual(info, fp.SOArchiveInfo(None, None, None, legacy=True))
        m = fp.so_provenance_members(**self.kw)
        m.pop("producer_label")
        path = _write_members(self.d, m)
        with self.assertRaisesRegex(ValueError, r"Regenerate it with <unknown producer>\."):
            fp.check_so_archive(np.load(path), path, run_band=self.band,
                                run_state=self.run_state, run_has_extern=False,
                                nmat_rule="exact")

    def test_green_verified_or_asserted(self):
        g = np.zeros((1, 4, 4, 2, 2), dtype=complex)
        g[0, :, :, 0, 0] = 1.0
        m = fp.so_provenance_members(**dict(self.kw, fp_green=fp.green_digest(
            g, "b_f_k_o_o_v1")))
        path = _write_members(self.d, m)
        call = lambda **kw: fp.check_so_archive(
            np.load(path), path, run_band=self.band, run_state=self.run_state,
            run_has_extern=False, nmat_rule="exact", **kw)
        call(green=g, green_layout="b_f_k_o_o_v1")
        with self.assertLogs("qlms", "INFO") as cm:
            call()
        self.assertTrue(any("fp_green: producer assertion" in r for r in cm.output))
        g2 = g.copy()
        g2[0, 1, 2, 1, 0] = 1e-12
        with self.assertRaisesRegex(ValueError, r"fp_green mismatch"):
            call(green=g2, green_layout="b_f_k_o_o_v1")

    def test_fixed_mu_archive_refused_by_hwave_sc(self):
        _rpa_for_sc(self.inp_soc, mu=0.1)
        with self.assertRaisesRegex(ValueError, r"calc_mu: file False, run True"):
            _sc("load", self.inp_soc)

    def test_residual_cap_boundary(self):
        t = 2.0
        cap = 1e-8 * max(1.0, t)
        for r, ok in ((0.999 * cap, True), (1.001 * cap, False)):
            with self.subTest(r=r):
                path = _write_members(self.d, fp.so_provenance_members(
                    **dict(self.kw, state=fp.so_state(beta=10.0, nmat=64, mu=0.1,
                                                      density=t + r, calc_mu=True,
                                                      target=t))))
                call = lambda: fp.check_so_archive(
                    np.load(path), path, run_band=self.band, run_state=self.run_state,
                    run_has_extern=False, nmat_rule="exact")
                if ok:
                    call()
                else:
                    with self.assertRaisesRegex(ValueError, r"did not reach its filling"):
                        call()

    def test_non_finite_member_refused_at_write(self):
        with self.assertRaisesRegex(ValueError, r"state_density"):
            fp.so_provenance_members(**dict(self.kw, state=self.kw["state"]._replace(
                density=float("nan"))))

    def test_header_dim_bad_line(self):
        for bad in ("", "x", "0"):
            with self.subTest(bad=bad), self.assertRaisesRegex(
                    ValueError, r"line 2 must be the number of orbitals"):
                wan90.read_w90_header_dim(_write_raw(self.d, ["comment", bad, "1", "1",
                                                               "0 0 0 1 1 0.1 0.0"]))

    def test_sigma_green_binding(self):
        out = _flex_for_sc(self.inp_soc)
        _set_member(out["sigma"], "fp_green", "0" * 64)
        with self.assertRaisesRegex(ValueError,
                                    r"sigma\.npz and green\.npz come from different runs"):
            _sc("flex", self.inp_soc)

    def test_flex_chi0q_companion_checked(self):
        """The FLEX chi0q is provenance-checked like chiq: one from another
        coeff_extern is refused by its band, an unversioned one beside a
        versioned chiq is refused."""
        out = _flex_for_sc(self.inp_soc)
        other = _flex_for_sc(_soc_input(self.new_dir()), coeff_extern=0.0)
        shutil.copy(other["chi0q"], out["chi0q"])
        with self.assertRaisesRegex(ValueError, r"chi0q\.npz' \(written by .*band differs: "
                                                r"coeff_extern: file 0\.0, run 1\.0"):
            _sc("flex", self.inp_soc)
        out = _flex_for_sc(self.inp_soc)
        _strip_fingerprint_members(out["chi0q"])
        with self.assertRaisesRegex(ValueError, r"chiq\.npz' \(versioned\) and '.*chi0q\.npz' "
                                                r"\(unversioned\) disagree on fingerprint_version"):
            _sc("flex", self.inp_soc)

    def test_flex_legacy_chiq_with_versioned_chi0q_refused(self):
        """The reversed mix (unversioned chiq, versioned chi0q) is refused
        too, before the Green function is read: otherwise the Green-digest
        and self-energy checks would be skipped."""
        out = _flex_for_sc(self.inp_soc)
        _strip_fingerprint_members(out["chiq"])
        os.remove(out["sigma"])
        os.remove(out["green"])              # a Green read would fail differently
        with self.assertRaisesRegex(ValueError, r"chiq\.npz' \(unversioned\) and '.*chi0q\.npz' "
                                                r"\(versioned\) disagree on fingerprint_version; "
                                                r"regenerate both from one FLEX run"):
            _sc("flex", self.inp_soc)

    def test_nmat_only_difference_names_the_key(self):
        _flex_for_sc(self.inp_soc)
        with self.assertRaisesRegex(ValueError, r"state differs: nmat 16 \(file\) vs 32 "
                                                r"\(run\)\. .*set \[mode\.param\] Nmat = 16"):
            _sc("flex", self.inp_soc, nmat=32)

    def test_legacy_sigma_with_versioned_chiq(self):
        out = _flex_for_sc(self.inp_soc)
        _strip_fingerprint_members(out["sigma"])
        with self.assertRaisesRegex(ValueError, r"sigma\.npz' has no fingerprint_version; "
                                                r"regenerate the FLEX run's sigma"):
            _sc("flex", self.inp_soc)

    def test_wrong_dtype_member_refused(self):
        good = fp.so_provenance_members(**self.kw)
        for k, v in (("state_beta", np.str_("10.0")), ("fp_band", np.int64(3)),
                     ("state_nmat", np.int32(good["state_nmat"])),
                     ("state_beta", np.float32(good["state_beta"])),
                     ("prov_transfer_vals", good["prov_transfer_vals"].astype(np.complex64)),
                     ("prov_cell_shape", good["prov_cell_shape"].astype(np.int32)),
                     ("state_calc_mu", np.int64(good["state_calc_mu"]))):
            with self.subTest(member=k, dtype=str(np.asarray(v).dtype)):
                path = _write_members(self.d, dict(good, **{k: v}))
                with self.assertRaisesRegex(ValueError, r"malformed or truncated: member {} "
                                                        r"has the wrong dtype or rank"
                                            .format(k)):
                    fp.check_so_archive(np.load(path), path, run_band=self.band,
                                        run_state=self.run_state, run_has_extern=False,
                                        nmat_rule="exact")

    def test_big_endian_member_accepted(self):
        """Every value is read through a Python scalar, so the byte order of
        a member does not matter; only its width and kind do."""
        good = fp.so_provenance_members(**self.kw)
        swapped = dict(good, state_nmat=np.asarray(good["state_nmat"]).astype(">i8"),
                       prov_transfer_vals=good["prov_transfer_vals"].astype(">c16"))
        path = _write_members(self.d, swapped)
        with np.load(path) as z:
            self.assertEqual(z["state_nmat"].dtype.byteorder, ">")
            info = fp.check_so_archive(z, path, run_band=self.band,
                                       run_state=self.run_state, run_has_extern=False,
                                       nmat_rule="exact")
        self.assertFalse(info.legacy)

    def test_dynamic_requires_filling(self):
        _flex_for_sc(self.inp_soc)
        with self.assertRaisesRegex(ValueError, r"\[mode\.param\] filling must be specified "
                                                r"\(the spin-orbital archive checks"):
            _sc("flex", self.inp_soc, frequency="dynamic", filling=None)

    def test_missing_sigma_refused(self):
        out = _flex_for_sc(self.inp_soc)
        os.remove(out["sigma"])
        with self.assertRaisesRegex(ValueError, r"self-energy archive .* not found"):
            _sc("flex", self.inp_soc)

    def test_load_mode_keeps_existing_nmat_rule(self):
        _rpa_for_sc(self.inp_soc, nmat=64)
        with self.assertLogs("qlms.solver.fingerprint", "WARNING") as cm:
            out = _sc("load", self.inp_soc, quiet=False, nmat=128)
        self.assertEqual(int(out.npz["state_nmat"]), 128)        # its own G, not the RPA's
        msgs = [m for m in cm.output if "another Matsubara grid" in m]
        self.assertEqual(len(msgs), 1, cm.output)
        self.assertRegex(msgs[0], r"chi0q\.npz' \(written by .*\) was summed on another "
                                  r"Matsubara grid: state_nmat 64 \(file\) vs Nmat 128 "
                                  r"\(run\); .*Set \[mode\.param\] Nmat = 64 to reproduce")

    def test_load_mode_same_nmat_does_not_warn(self):
        _rpa_for_sc(self.inp_soc, nmat=64)
        with self.assertLogs("qlms", "WARNING") as cm:
            logging.getLogger("qlms").warning("sentinel")     # assertNoLogs needs 3.10
            _sc("load", self.inp_soc, quiet=False, nmat=64)
        self.assertFalse(any("another Matsubara grid" in m for m in cm.output), cm.output)

    def test_zeeman_only_rpa_load_refused(self):                     # documented limitation, pin
        _rpa_for_sc(self.inp_zeeman_only)
        with self.assertRaisesRegex(ValueError, r"spin-diag|layout"):
            _sc("load", self.inp_zeeman_only)


class TestScWriters(_Dir):
    """The hwave_sc gap archives record this run's band and the Green
    function and state its kernel used."""

    def setUp(self):
        super().setUp()
        self.inp_soc = _soc_input(self.new_dir())

    def _common(self, npz, stage):
        self.assertEqual(int(npz["fingerprint_version"]), 1)
        self.assertEqual(str(npz["producer_stage"]), stage)
        T, E, meta = _tables_and_meta_from_input(self.inp_soc)
        self.assertEqual(str(npz["fp_band"]), fp.band_digests(T, None, meta)["fp_band"])
        self.assertTrue(bool(npz["state_calc_mu"]))
        self.assertAlmostEqual(float(npz["state_target"]),
                               self.inp_soc["param"]["filling"] * 2, delta=1e-15)
        self.assertEqual(float(npz["state_beta"]), 1.0 / self.inp_soc["param"]["T"])
        for k, (kind, ndim) in _SCHEMA.items():
            self.assertEqual((npz[k].dtype.kind, npz[k].ndim), (kind, ndim), k)

    def test_calc_mode_own_state(self):
        out = _sc("calc", self.inp_soc)
        self._common(out.npz, "sc_so_static")
        self.assertEqual(str(out.npz["fp_green"]),
                         fp.green_digest(out.kernel_green, "sc_o_o_xyz_f"))
        self.assertLessEqual(abs(float(out.npz["state_density"])
                                 - float(out.npz["state_target"])), 1e-10)

    def test_load_records_own_green_and_state(self):
        _rpa_for_sc(self.inp_soc)
        out = _sc("load", self.inp_soc)
        self._common(out.npz, "sc_so_static")
        self.assertEqual(str(out.npz["fp_green"]),
                         fp.green_digest(out.kernel_green, "sc_o_o_xyz_f"))
        self.assertLessEqual(abs(float(out.npz["state_density"])
                                 - float(out.npz["state_target"])), 1e-10)
        self.assertAlmostEqual(float(out.npz["state_mu_tol"]),
                               abs(float(out.npz["state_density"])
                                   - float(out.npz["state_target"])), delta=1e-15)

    def test_flex_records_borrowed_green_and_used_mu(self):
        src = {}
        with np.load(_flex_for_sc(self.inp_soc)["green"]) as z:
            src = {k: np.array(z[k]) for k in z.files}
        out = _sc("flex", self.inp_soc)
        self._common(out.npz, "sc_so_static")
        self.assertEqual(str(out.npz["fp_green"]), fp.green_digest(src["green"], "b_f_k_o_o_v1"))
        self.assertEqual(str(out.npz["fp_green"]),
                         fp.green_digest(out.kernel_green, "sc_o_o_xyz_f"))
        self.assertEqual(float(out.npz["state_mu"]), float(src["state_mu"]))
        self.assertAlmostEqual(float(out.npz["state_density"]), float(src["state_density"]),
                               delta=1e-12 * max(1.0, abs(float(src["state_density"]))))
        self.assertAlmostEqual(float(out.npz["state_mu_tol"]),
                               abs(float(out.npz["state_density"])
                                   - float(out.npz["state_target"])), delta=1e-15)

    def test_flex_density_is_measured_on_the_consumed_green(self):
        """The closure density is FLEX's own routine applied to the G the
        kernel consumed, at the archive mu, with the stored static sigma."""
        from hwave.solver import hartree_fock as hf
        out = _flex_for_sc(self.inp_soc)
        with np.load(out["sigma"]) as z:
            sig = np.array(z["sigma_static"])[0, 0]
        with np.load(out["green"]) as z:
            g, mu = np.array(z["green"]), float(z["state_mu"])
        calls = []
        real = hf.equal_time_density

        def spy(green_kw, heff, mu_, beta, shape, **kw):
            calls.append((np.array(green_kw), mu_, kw))
            return real(green_kw, heff, mu_, beta, shape, **kw)

        heffs = []
        real_heff = hf.heff_eigenpairs

        def spy_heff(H0, s, *a, **kw):
            heffs.append(np.array(s))
            return real_heff(H0, s, *a, **kw)

        with mock.patch.object(hf, "equal_time_density", side_effect=spy), \
                mock.patch.object(hf, "heff_eigenpairs", side_effect=spy_heff):
            _sc("flex", self.inp_soc)
        self.assertEqual(len(calls), 1)
        np.testing.assert_array_equal(calls[0][0], g)
        self.assertEqual(calls[0][1], mu)
        self.assertEqual(calls[0][2], {"sym_tol": None})
        np.testing.assert_array_equal(heffs[0], sig)

    def test_dynamic_writers(self):
        for mode in ("calc", "flex"):
            with self.subTest(mode=mode):
                if mode == "flex":
                    with np.load(_flex_for_sc(self.inp_soc)["green"]) as z:
                        src_fp, src_mu = str(z["fp_green"]), float(z["state_mu"])
                out = _sc(mode, self.inp_soc, frequency="dynamic")
                self._common(out.npz, "sc_so_dynamic")
                self.assertEqual(str(out.npz["fp_green"]),
                                 fp.green_digest(out.kernel_green, "sc_o_o_xyz_f"))
                if mode == "flex":
                    self.assertEqual(str(out.npz["fp_green"]), src_fp)
                    self.assertEqual(float(out.npz["state_mu"]), src_mu)
                else:
                    self.assertLessEqual(abs(float(out.npz["state_density"])
                                             - float(out.npz["state_target"])), 1e-10)

    def test_eigenvalue_dat_carries_the_scalars(self):
        out = _sc("calc", self.inp_soc)
        path = os.path.join(self.inp_soc["dir"], "output_sc_calc", "eigenvalue.dat")
        with open(path) as f:
            text = f.read()
        self.assertIn("# producer_stage = sc_so_static", text)
        self.assertIn("# fp_green = {}".format(str(out.npz["fp_green"])), text)
        self.assertIn("# state_calc_mu = true", text)
        self.assertNotIn("prov_transfer_keys", text)

    def test_legacy_input_gives_truthful_outputs(self):
        _strip_fingerprint_members(_rpa_for_sc(self.inp_soc)["path"])
        out = _sc("load", self.inp_soc).npz
        self.assertEqual(int(out["fingerprint_version"]), 1)
        self.assertTrue(np.isfinite(float(out["state_density"])))

    def test_legacy_flex_writes_no_provenance(self):
        out = _flex_for_sc(self.inp_soc)
        for k in ("chiq", "green", "sigma", "chi0q"):
            _strip_fingerprint_members(out[k])
        npz = _sc("flex", self.inp_soc).npz
        self.assertFalse(set(fp.SO_PROVENANCE_MEMBERS) & set(npz))
        self.assertTrue(bool(npz["spin_orbital"]))


if __name__ == "__main__":
    unittest.main()
