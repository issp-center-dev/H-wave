"""The bond-resolved longitudinal channel inside the FLEX loop (spec
2026-09-06, gates G0, G3, G4): the declared-zero reduction to standalone
HF FLEX, the first-map static slice against the Phase A RPA gate, the
refusal precedence tripwires, the output-path validation, the dedicated
archive, IterationMax = 0, reuse/failure atomicity and the cond_min
provenance."""
import os
import shutil
import tempfile
import unittest
from unittest import mock

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k

_IN2 = "tests/rpa/input_2orb"
_NMAT = 8
_PAR = {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1],
        "Nmat": _NMAT, "IterationMax": 3, "Mix": 0.5, "EPS": 1e-12}


def _reader(path, inter):
    idict = {"path_to_input": path, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
    idict.update(inter)
    return read_input_k.QLMSkInput({"path_to_input": path, "interaction": idict})


def _flex(param_extra=None, path=_IN2, inter=None, gate=True, hf=True, mixing="linear"):
    import hwave.solver.flex as flex_mod
    r = _reader(path, {"CoulombInter": "coulombinter.dat"} if inter is None else inter)
    par = dict(_PAR, flex_hartree_fock=hf, longitudinal_bond_channels=gate, mixing_scheme=mixing)
    par.update(param_extra or {})
    par = {k: v for k, v in par.items() if v is not None}
    info = {"mode": "FLEX", "param": par, "enable_spin_orbital": False, "calc_scheme": "general"}
    return flex_mod.FLEX(r.get_param("ham"), {}, info), r


def _rpa(param_extra=None, path=_IN2, inter=None):
    import hwave.solver.rpa as rpa_mod
    r = _reader(path, {"CoulombInter": "coulombinter.dat"} if inter is None else inter)
    par = {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1],
           "Nmat": _NMAT, "longitudinal_bond_channels": True}
    par.update(param_extra or {})
    par = {k: v for k, v in par.items() if v is not None}
    info = {"mode": "RPA", "param": par, "enable_spin_orbital": False,
            "calc_scheme": "general", "calc_type": "ring"}
    return rpa_mod.RPA(r.get_param("ham"), {}, info), r


def _zeroed_offsite_copy(dst):
    """input_2orb with every off-site CoulombInter coefficient set to zero
    (the shells stay declared)."""
    for f in ("geom.dat", "transfer.dat"):
        shutil.copy(os.path.join(_IN2, f), dst)
    lines = open(os.path.join(_IN2, "coulombinter.dat")).read().splitlines()
    out = lines[:4]
    for ln in lines[4:]:
        parts = ln.split()
        if len(parts) >= 7 and any(int(x) != 0 for x in parts[:3]):
            parts[5], parts[6] = "0.0", "0.0"
        out.append("  ".join(parts))
    with open(os.path.join(dst, "coulombinter.dat"), "w") as fw:
        fw.write("\n".join(out) + "\n")


def _collect(s, gi, out):
    rec = []
    s._iteration_hook = lambda d: rec.append(d)
    s.solve(gi, out)
    return rec


class TestG0DeclaredZero(unittest.TestCase):

    def test_g0_declared_zero_equals_standalone_hf(self):
        with tempfile.TemporaryDirectory() as inp, tempfile.TemporaryDirectory() as out:
            _zeroed_offsite_copy(inp)
            recs = {}
            gis = {}
            for gate in (True, False):
                s, r = _flex(path=inp, gate=gate)
                gi = r.get_param("green")
                recs[gate] = _collect(s, gi, out)
                gis[gate] = (s, gi)
            self.assertEqual(len(recs[True]), len(recs[False]))
            self.assertGreaterEqual(len(recs[True]), 2)
            for a, b in zip(recs[True], recs[False]):
                self.assertEqual(a["iteration"], b["iteration"])
                self.assertAlmostEqual(a["mu"], b["mu"], delta=1e-12)
                for k in ("static_new", "fluct_new", "chi0q", "chiq_s", "chiq_c"):
                    np.testing.assert_allclose(a[k], b[k], rtol=0, atol=1e-12, err_msg=k)
            sg, gg = gis[True]
            sh, gh = gis[False]
            for k in ("sigma", "sigma_static", "sigma_fluct", "green"):
                np.testing.assert_allclose(gg[k], gh[k], rtol=0, atol=1e-12, err_msg=k)
            self.assertAlmostEqual(gg["physics"]["mu"], gh["physics"]["mu"], delta=1e-12)
            # gate-owned additions present only with the gate
            self.assertIn("longitudinal_bond_chi_s", gg)
            self.assertEqual(str(gg["longitudinal_bond_source"]), "last_map")
            self.assertFalse(any(str(k).startswith("longitudinal_bond_") for k in gh))


class TestG3StaticSlice(unittest.TestCase):

    def test_g3_static_slice_equals_phase_a(self):
        with tempfile.TemporaryDirectory() as out:
            rs, rr = _rpa({"mu": 0.3, "filling": None})
            rgi = rr.get_param("green")
            rs.solve(rgi, out)
            fs, fr = _flex({"mu": 0.3, "filling": None, "IterationMax": 1})
            fgi = fr.get_param("green")
            rec = _collect(fs, fgi, out)
        self.assertEqual(len(rec), 1)
        d = rec[0]
        np.testing.assert_allclose(d["bond_static_s"], rgi["longitudinal_bond_chi_s"], rtol=0, atol=1e-12)
        np.testing.assert_allclose(d["bond_static_c"], rgi["longitudinal_bond_chi_c"], rtol=0, atol=1e-12)
        np.testing.assert_allclose(d["chiq_s"][_NMAT // 2], rgi["longitudinal_bond_chiq_s_static"],
                                   rtol=0, atol=1e-12)
        np.testing.assert_allclose(d["chiq_c"][_NMAT // 2], rgi["longitudinal_bond_chiq_c_static"],
                                   rtol=0, atol=1e-12)


class TestRefusals(unittest.TestCase):

    def _tripwired(self, s):
        import hwave.solver.backend as bk
        boom = mock.Mock(side_effect=RuntimeError("expensive step reached"))
        patches = [mock.patch.object(type(s), n, boom) for n in
                   ("_calc_epsilon_k", "_find_mu_dressed", "_calc_green", "_calc_chi0q")]
        patches.append(mock.patch.object(bk, "get_backend", boom))
        return patches

    def _assert_refused(self, s, gi, needle):
        patches = self._tripwired(s)
        for p in patches:
            p.start()
        try:
            with tempfile.TemporaryDirectory() as out:
                with self.assertRaises(ValueError) as cm:
                    s.solve(gi, out)
        finally:
            for p in patches:
                p.stop()
        self.assertIn(needle, str(cm.exception))
        self.assertFalse(any(str(k).startswith("longitudinal_bond_") for k in gi))

    def test_refusal_precedence_tripwires(self):
        # (6) off-site Exchange
        s, r = _flex(inter={"CoulombInter": "coulombinter.dat", "Exchange": "coulombinter.dat"})
        self._assert_refused(s, r.get_param("green"), "Exchange")
        # (6) no declared off-site shell
        s, r = _flex(inter={"CoulombIntra": "coulombintra.dat"})
        self._assert_refused(s, r.get_param("green"), "shell")
        # (7) memory cap
        s, r = _flex({"longitudinal_bond_memory_cap_gb": 1e-9})
        self._assert_refused(s, r.get_param("green"), "GiB")
        # (4) spin-mixing mean field
        s, r = _flex()
        gi = r.get_param("green")
        nvol, norb = s.lattice.nvol, s.norb
        tm = np.zeros((nvol, 2, norb, 2, norb), complex)
        tm[:, 0, :, 1, :] = 0.1
        gi["trans_mod"] = tm
        self._assert_refused(s, gi, "spin")
        # (5) split seed with a mean field (D15)
        with tempfile.TemporaryDirectory() as out:
            s0, r0 = _flex({"IterationMax": 1})
            gi0 = r0.get_param("green")
            s0.solve(gi0, out)
            s0.save_results({"path_to_output": out, "sigma": "sigma.npz"}, gi0)
            s, r = _flex()
            gi = r.get_param("green")
            gi.update(s.read_init({"path_to_input": out, "sigma_init": "sigma.npz"}))
            gi["trans_mod"] = np.zeros((nvol, 2, norb, 2, norb), complex)
            self._assert_refused(s, gi, "split")

    def test_complex_offsite_coefficient_refused(self):
        with tempfile.TemporaryDirectory() as inp:
            for f in ("geom.dat", "transfer.dat"):
                shutil.copy(os.path.join(_IN2, f), inp)
            lines = open(os.path.join(_IN2, "coulombinter.dat")).read().splitlines()
            out = lines[:4]
            for ln in lines[4:]:
                p = ln.split()
                if p[:3] == ["0", "1", "0"] and p[3:5] == ["1", "1"]:
                    p[6] = "0.2"
                if p[:3] == ["0", "-1", "0"] and p[3:5] == ["1", "1"]:
                    p[6] = "-0.2"
                out.append("  ".join(p))
            open(os.path.join(inp, "coulombinter.dat"), "w").write("\n".join(out) + "\n")
            s, r = _flex(path=inp)
            self._assert_refused(s, r.get_param("green"), "real")


class TestOutputs(unittest.TestCase):

    def test_output_paths_collision_refused_at_solve_entry_and_qlms(self):
        from hwave.solver import flex_bond
        with tempfile.TemporaryDirectory() as out:
            paths = flex_bond.resolve_output_paths(
                {"path_to_output": out, "chiq_s": "a", "sigma": "b.npz"}, out,
                ("chiq_s", "chiq_c", "sigma"))
            self.assertEqual(paths["chiq_s"], os.path.abspath(os.path.join(out, "a.npz")))
            self.assertEqual(paths["chiq_c"], os.path.abspath(os.path.join(out, "chiq_c.npz")))
            with self.assertRaises(ValueError) as cm:
                flex_bond.resolve_output_paths({"chiq_s": "x", "sigma": "x.npz"}, out,
                                               ("chiq_s", "sigma"))
            self.assertIn("chiq_s", str(cm.exception))
            self.assertIn("sigma", str(cm.exception))
            s, r = _flex()
            with self.assertRaises(ValueError):
                s.validate_output_paths()          # no mapping stored, no argument
            s.validate_output_paths({"path_to_output": out, "chiq_s": "x", "sigma": "y"})
            with self.assertRaises(ValueError):
                s.validate_output_paths({"path_to_output": out, "chiq_s": "x", "sigma": "x"})
            # stored mapping -> solve entry refuses before any work
            with mock.patch.object(type(s), "_calc_epsilon_k",
                                   mock.Mock(side_effect=RuntimeError("reached"))):
                with self.assertRaises(ValueError):
                    s.solve(r.get_param("green"), out)
            # the dedicated archive counts only when it will be written
            s2, _ = _flex({"longitudinal_bond_output_full": True})
            with self.assertRaises(ValueError):
                s2.validate_output_paths({"path_to_output": out, "chiq_s": "longitudinal_bond"})
            s3, _ = _flex()
            s3.validate_output_paths({"path_to_output": out, "chiq_s": "longitudinal_bond"})
            # qlms.run
            import hwave.qlms as qlms
            d = {"log": {"print_level": 1},
                 "mode": {"mode": "FLEX", "calc_scheme": "general", "enable_spin_orbital": False,
                          "param": dict(_PAR, flex_hartree_fock=True, longitudinal_bond_channels=True)},
                 "file": {"input": {"path_to_input": _IN2,
                                    "interaction": {"path_to_input": _IN2, "Geometry": "geom.dat",
                                                    "Transfer": "transfer.dat",
                                                    "CoulombInter": "coulombinter.dat"}},
                          "output": {"path_to_output": out, "green": "same.npz", "sigma": "same"}}}
            with mock.patch("hwave.solver.flex.FLEX.solve", side_effect=RuntimeError("reached")):
                with self.assertRaises(ValueError):
                    qlms.run(input_dict=d)

    def test_dedicated_archive_schema_and_cross_output_consistency(self):
        with tempfile.TemporaryDirectory() as out:
            s, r = _flex({"longitudinal_bond_output_full": True, "IterationMax": 2})
            gi = r.get_param("green")
            s.solve(gi, out)
            s.save_results({"path_to_output": out, "sigma": "sigma.npz", "green": "green.npz",
                            "chi0q": "chi0q.npz"}, gi)
            for f in ("longitudinal_bond.npz", "chiq_s.npz", "chiq_c.npz", "sigma.npz",
                      "green.npz", "chi0q.npz"):
                self.assertTrue(os.path.exists(os.path.join(out, f)), f)
            b = np.load(os.path.join(out, "longitudinal_bond.npz"))
            self.assertEqual(int(b["bond_archive_schema"]), 1)
            self.assertEqual(str(b["freq_axis"]), "bosonic l -> 2l - nmat")
            nvol, nd = s.lattice.nvol, s.norb ** 2
            B = b["delta_r"].shape[0]
            self.assertEqual(b["chi_s_w"].shape, (_NMAT, nvol, B * nd, B * nd))
            self.assertEqual(b["chi_c_w"].shape, (_NMAT, nvol, B * nd, B * nd))
            for k in ("beta", "T", "nmat", "cell_shape", "momentum_convention", "index_order",
                      "reverse", "types", "longitudinal_bond_schema", "scf_converged",
                      "payload_kind", "hf_density_error", "density_target_enforced",
                      "longitudinal_bond_cond_min_s"):
                self.assertIn(k, b.files, k)
            self.assertEqual(str(b["payload_kind"]), "last_map")
            np.testing.assert_array_equal(b["chi_s_w"][_NMAT // 2], b["longitudinal_bond_chi_s"])
            np.testing.assert_array_equal(b["chi_c_w"][_NMAT // 2], b["longitudinal_bond_chi_c"])
            cs = np.load(os.path.join(out, "chiq_s.npz"))
            self.assertIn("longitudinal_bond_chi_s", cs.files)
            self.assertIn("longitudinal_bond_chi_c", cs.files)
            self.assertEqual(str(cs["longitudinal_bond_source"]), "last_map")
            self.assertNotIn("chi_s_w", cs.files)
            np.testing.assert_array_equal(cs["longitudinal_bond_chi_s"], b["longitudinal_bond_chi_s"])
            np.testing.assert_array_equal(cs["chiq_s"][_NMAT // 2].reshape(nvol, nd, nd),
                                          b["longitudinal_bond_chi_s"][:, :nd, :nd])
            sg = np.load(os.path.join(out, "sigma.npz"))
            self.assertEqual(str(sg["sigma_convention"]), "split")
            np.testing.assert_array_equal(sg["sigma"], sg["sigma_static"] + sg["sigma_fluct"])
            self.assertEqual(str(sg["payload_kind"]), "final_state")
            c0 = np.load(os.path.join(out, "chi0q.npz"))
            self.assertEqual(str(c0["payload_kind"]), "last_map")
            # combined chiq routing: static keys go to chiq, not chiq_s
            s.save_results({"path_to_output": out, "chiq": "chiq.npz"}, gi)
            cq = np.load(os.path.join(out, "chiq.npz"))
            self.assertIn("longitudinal_bond_chi_s", cq.files)
            cs2 = np.load(os.path.join(out, "chiq_s.npz"))
            self.assertNotIn("longitudinal_bond_chi_s", cs2.files)

    def test_iteration_max_zero_omits_last_map_archives(self):
        with tempfile.TemporaryDirectory() as out:
            s, r = _flex({"IterationMax": 0, "longitudinal_bond_output_full": True})
            gi = r.get_param("green")
            s.solve(gi, out)
            self.assertFalse(any(str(k).startswith("longitudinal_bond_") for k in gi))
            s.validate_output_paths({"path_to_output": out, "chiq_s": "longitudinal_bond"})
            s.save_results({"path_to_output": out, "sigma": "sigma.npz", "green": "green.npz",
                            "chi0q": "chi0q.npz", "chiq": "chiq.npz"}, gi)
            for f in ("chi0q.npz", "chiq.npz", "chiq_s.npz", "chiq_c.npz", "longitudinal_bond.npz"):
                self.assertFalse(os.path.exists(os.path.join(out, f)), f)
            sg = np.load(os.path.join(out, "sigma.npz"))
            self.assertEqual(int(sg["scf_iterations"]), 0)
            self.assertTrue(np.isnan(float(sg["scf_sigma_residual"])))

    def test_reused_green_info_and_failure_atomicity(self):
        from hwave.solver import flex_bond
        with tempfile.TemporaryDirectory() as out:
            s, r = _flex({"IterationMax": 1})
            gi = r.get_param("green")
            s.solve(gi, out)
            self.assertIn("longitudinal_bond_chi_s", gi)
            stores = []
            real = flex_bond.BondBlockStore

            def _spy(*a, **k):
                st = real(*a, **k)
                stores.append(st)
                return st
            with mock.patch.object(flex_bond, "BondBlockStore", _spy), \
                    mock.patch.object(flex_bond, "calc_self_energy_bond",
                                      side_effect=RuntimeError("injected")):
                with self.assertRaises(RuntimeError):
                    s.solve(gi, out)
            self.assertEqual(len(stores), 1)
            self.assertTrue(stores[0].released)
            self.assertFalse(any(str(k).startswith("longitudinal_bond_") for k in gi))
            for k in ("sigma", "green", "chiq_s", "physics"):
                self.assertNotIn(k, gi)
            self.assertFalse(hasattr(s, "_bond_last"))
            s.solve(gi, out)
            self.assertIn("longitudinal_bond_chi_s", gi)
            self.assertIn("sigma", gi)

    def test_cond_min_provenance(self):
        with tempfile.TemporaryDirectory() as out:
            s, r = _flex({"IterationMax": 2})
            gi = r.get_param("green")
            s.solve(gi, out)
        for ch in ("s", "c"):
            v = float(gi["longitudinal_bond_cond_min_" + ch])
            self.assertTrue(np.isfinite(v) and v > 0.0)
            self.assertEqual(v, getattr(s._bond_last, "cond_min_" + ch))


class TestFailureClearing(unittest.TestCase):

    def test_caches_and_seed_are_dropped_on_every_failure(self):
        from hwave.solver import flex_bond
        attrs = ("_phase_b_seed", "_hf_tables", "_bond_topo", "_bond_split", "_bond_view",
                 "_bond_S", "_bond_C", "_bond_S_on", "_bond_C_on", "_bond_est", "_bond_last")
        with tempfile.TemporaryDirectory() as out:
            # preflight failure (no declared shell)
            s, r = _flex(inter={"CoulombIntra": "coulombintra.dat"})
            with self.assertRaises(ValueError):
                s.solve(r.get_param("green"), out)
            for a in attrs:
                self.assertFalse(hasattr(s, a), a)
            # mid-loop failure in the transport
            s, r = _flex({"IterationMax": 2})
            gi = r.get_param("green")
            with mock.patch.object(flex_bond, "calc_self_energy_bond",
                                   side_effect=RuntimeError("injected")):
                with self.assertRaises(RuntimeError):
                    s.solve(gi, out)
            for a in attrs:
                self.assertFalse(hasattr(s, a), a)
            self.assertFalse(any(str(k).startswith("longitudinal_bond_") for k in gi))
            # the seed is gone once the loop runs (a successful run)
            s, r = _flex({"IterationMax": 1})
            gi = r.get_param("green")
            s.solve(gi, out)
            self.assertFalse(hasattr(s, "_phase_b_seed"))
            # a reused solver whose stored output mapping now collides: the
            # solve-entry refusal clears the previous results too
            s.validate_output_paths({"path_to_output": out, "sigma": "s"})
            s._info_outputfile = {"path_to_output": out, "chiq_s": "same", "sigma": "same"}
            with self.assertRaises(ValueError):
                s.solve(gi, out)
            for k in ("sigma", "green", "chiq_s", "physics"):
                self.assertNotIn(k, gi)
            for a in attrs + ("sigma", "green_kw", "chi_s", "physics"):
                self.assertFalse(hasattr(s, a), a)


class TestStandaloneHFAdmissibility(unittest.TestCase):

    def test_complex_offsite_coefficient_and_offsite_exchange_refused_at_preflight(self):
        with tempfile.TemporaryDirectory() as inp:
            for f in ("geom.dat", "transfer.dat"):
                shutil.copy(os.path.join(_IN2, f), inp)
            lines = open(os.path.join(_IN2, "coulombinter.dat")).read().splitlines()
            out = lines[:4]
            for ln in lines[4:]:
                p = ln.split()
                if p[:3] == ["0", "1", "0"] and p[3:5] == ["1", "1"]:
                    p[6] = "0.2"
                if p[:3] == ["0", "-1", "0"] and p[3:5] == ["1", "1"]:
                    p[6] = "-0.2"
                out.append("  ".join(p))
            open(os.path.join(inp, "pairlift.dat"), "w").write("\n".join(out).replace("CoulombInter", "PairLift") + "\n")
            open(os.path.join(inp, "coulombinter.dat"), "w").write("\n".join(out) + "\n")
            shutil.copy(os.path.join(_IN2, "coulombinter.dat"), os.path.join(inp, "coulombinter_real.dat"))
            boom = mock.Mock(side_effect=RuntimeError("expensive step reached"))
            for inter in ({"CoulombInter": "coulombinter.dat"},
                          {"CoulombInter": "coulombinter_real.dat", "PairLift": "pairlift.dat"},
                          {"Exchange": "coulombinter_real.dat"}):
                s, r = _flex(path=inp, inter=inter, gate=False)
                with mock.patch.object(type(s), "_calc_epsilon_k", boom), \
                        tempfile.TemporaryDirectory() as o:
                    with self.assertRaises(ValueError) as cm:
                        s.solve(r.get_param("green"), o)
                self.assertIn("flex_hartree_fock", str(cm.exception))


if __name__ == "__main__":
    unittest.main()
