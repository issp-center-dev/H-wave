"""The experimental bond-resolved longitudinal channel of the RPA solver
(``[mode.param] longitudinal_bond_channels = true``; GitHub issue #181,
Tier 3 Phase A): configuration, prerequisites, memory preflight, wiring,
and the ``longitudinal_bond_*`` outputs (spec gates G0 and G4).
"""
import logging
import os
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.rpa as rpa_mod

_IN = "tests/rpa/input_2orb"


def _build(param_extra=None, interactions=None, calc_type="ring",
           mode="RPA", nmat=32, path=_IN, cell=(4, 4, 1)):
    idict = {"path_to_input": path, "Geometry": "geom.dat",
             "Transfer": "transfer.dat"}
    idict.update({"CoulombInter": "coulombinter.dat"}
                 if interactions is None else interactions)
    r = read_input_k.QLMSkInput({"path_to_input": path, "interaction": idict})
    par = {"T": 2.0, "filling": 0.5, "CellShape": list(cell),
           "SubShape": [1, 1, 1], "Nmat": nmat}
    par.update(param_extra or {})
    info = {"mode": mode, "param": par, "enable_spin_orbital": False,
            "calc_scheme": "general", "calc_type": calc_type}
    if mode == "FLEX":
        import hwave.solver.flex as flex_mod
        par.update({"IterationMax": 1, "Mix": 1.0, "EPS": 1})
        del info["calc_type"]
        return flex_mod.FLEX(r.get_param("ham"), {}, info), r
    return rpa_mod.RPA(r.get_param("ham"), {}, info), r


class TestGateConfig(unittest.TestCase):

    def test_defaults_off(self):
        s, _ = _build()
        self.assertFalse(s.longitudinal_bond_channels)
        self.assertIsNone(s.longitudinal_bond_max_shells)
        self.assertEqual(s.longitudinal_bond_memory_cap_gb,
                         rpa_mod.LONGITUDINAL_BOND_MEMORY_CAP_GB_DEFAULT)
        self.assertEqual(rpa_mod.LONGITUDINAL_BOND_MEMORY_CAP_GB_DEFAULT, 8.0)

    def test_on_and_options(self):
        s, _ = _build({"longitudinal_bond_channels": True,
                       "longitudinal_bond_max_shells": 2,
                       "longitudinal_bond_memory_cap_gb": 1.5})
        self.assertTrue(s.longitudinal_bond_channels)
        self.assertEqual(s.longitudinal_bond_max_shells, 2)
        self.assertEqual(s.longitudinal_bond_memory_cap_gb, 1.5)
        s, _ = _build({"Longitudinal_Bond_Channels": True})   # case-insensitive
        self.assertTrue(s.longitudinal_bond_channels)

    def test_bounds_and_types(self):
        for bad in ({"longitudinal_bond_channels": "yes"},
                    {"longitudinal_bond_channels": 1},
                    {"longitudinal_bond_channels": True,
                     "longitudinal_bond_max_shells": 0},
                    {"longitudinal_bond_channels": True,
                     "longitudinal_bond_max_shells": 1.5},
                    {"longitudinal_bond_channels": True,
                     "longitudinal_bond_max_shells": True},
                    {"longitudinal_bond_channels": True,
                     "longitudinal_bond_memory_cap_gb": 0.0},
                    {"longitudinal_bond_channels": True,
                     "longitudinal_bond_memory_cap_gb": float("nan")},
                    {"longitudinal_bond_channels": True,
                     "longitudinal_bond_memory_cap_gb": "8"}):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError) as cm:
                    _build(bad)
                self.assertIn("longitudinal_bond", str(cm.exception))

    def test_dependent_keys_warn_when_off_or_wrong_calc_type(self):
        with self.assertLogs("hwave.solver.rpa", level="WARNING") as cm:
            _build({"longitudinal_bond_max_shells": 1})
        self.assertTrue(any("longitudinal_bond_max_shells" in m
                            for m in cm.output), cm.output)
        with self.assertLogs("hwave.solver.rpa", level="WARNING") as cm:
            s, _ = _build({"longitudinal_bond_channels": True},
                          calc_type="ring+ladder")
        self.assertFalse(s.longitudinal_bond_channels)
        self.assertTrue(any("calc_type" in m for m in cm.output), cm.output)

    def test_both_gates_true_is_rejected_at_construction(self):
        for ct in ("ring", "ring+ladder"):
            with self.subTest(calc_type=ct):
                with self.assertRaises(ValueError) as cm:
                    _build({"longitudinal_bond_channels": True,
                            "transverse_bond_channels": True}, calc_type=ct)
                self.assertIn("calc_type", str(cm.exception))
        # an explicit false next to a true is fine
        s, _ = _build({"longitudinal_bond_channels": True,
                       "transverse_bond_channels": False})
        self.assertTrue(s.longitudinal_bond_channels)

    def test_info_line_logged(self):
        with self.assertLogs("hwave.solver.rpa", level="INFO") as cm:
            _build({"longitudinal_bond_channels": True})
        self.assertTrue(any("longitudinal_bond_channels = True" in m
                            for m in cm.output))

    def test_flex_rejects_only_a_true_key(self):
        with self.assertRaises(ValueError) as cm:
            _build({"longitudinal_bond_channels": True}, mode="FLEX")
        self.assertIn("Phase B", str(cm.exception))
        _build({"longitudinal_bond_channels": False}, mode="FLEX")
        with self.assertRaises(ValueError):
            _build({"longitudinal_bond_channels": "true"}, mode="FLEX")

    def test_flex_rejects_a_case_varied_true_key_under_ring_ladder(self):
        """The RAW value is inspected before the RPA base parser can turn
        it stale under calc_type='ring+ladder'."""
        import hwave.solver.flex as flex_mod
        idict = {"path_to_input": _IN, "Geometry": "geom.dat",
                 "Transfer": "transfer.dat", "CoulombInter": "coulombinter.dat"}
        r = read_input_k.QLMSkInput({"path_to_input": _IN, "interaction": idict})
        par = {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1],
               "SubShape": [1, 1, 1], "Nmat": 32, "IterationMax": 1,
               "Mix": 1.0, "EPS": 1, "Longitudinal_Bond_CHANNELS": True}
        info = {"mode": "FLEX", "param": par, "enable_spin_orbital": False,
                "calc_scheme": "general", "calc_type": "ring+ladder"}
        with self.assertRaises(ValueError) as cm:
            flex_mod.FLEX(r.get_param("ham"), {}, info)
        self.assertIn("Phase A", str(cm.exception))


_EQ2 = "tests/equivalence_input/orb2"


class TestGatePrerequisites(unittest.TestCase):

    def _gate(self, param_extra=None, interactions=None, path=_IN):
        s, _ = _build(dict({"longitudinal_bond_channels": True},
                           **(param_extra or {})), interactions, path=path)
        s.spin_mode = "spin-free"          # what solve() sets before the check
        return s

    def test_refusals(self):
        cases = [
            ({}, {"CoulombInter": "onsite_inter.dat"}, _IN, "declared off-site"),
            ({}, {"PairHop": "offsite_pairhop.dat"}, _EQ2, "PairHop"),
            ({}, {"Exchange": "offsite_exchange.dat"}, _EQ2, "Exchange"),
            ({}, {"CoulombInter": "offsite_coulombinter_interorb.dat",
                  "Exchange": "offsite_exchange.dat"}, _EQ2, "Exchange"),
            ({"SubShape": [2, 1, 1]}, None, _IN, "sublattice"),
        ]
        for extra, inter, path, fragment in cases:
            with self.subTest(fragment=fragment):
                s = self._gate(extra, inter, path)
                with self.assertRaises(ValueError) as cm:
                    s._validate_longitudinal_bond_prereqs()
                self.assertIn("longitudinal_bond_channels", str(cm.exception))
                self.assertIn(fragment, str(cm.exception))

    def test_no_interaction_is_refused_not_a_silent_chi0_run(self):
        """A gate-on run without any interaction has calc_chiq=False; it
        must still reach the 'no declared off-site shell' refusal."""
        s, r = _build({"longitudinal_bond_channels": True}, interactions={})
        self.assertFalse(s.calc_chiq)
        gi = r.get_param("green")
        os.makedirs("tests/rpa/output", exist_ok=True)
        with self.assertRaises(ValueError) as cm:
            s.solve(gi, "tests/rpa/output")
        self.assertIn("declared off-site", str(cm.exception))
        self.assertNotIn("chi0q", gi)

    def test_refusals_fire_before_the_green_function(self):
        """Nothing expensive runs before a gate refusal: the chemical
        potential search, the Green's function and the bubble all trip."""
        for extra, inter, fragment in (
                ({"longitudinal_bond_memory_cap_gb": 1e-9}, None, "memory"),
                ({}, {"CoulombInter": "onsite_inter.dat"}, "declared off-site")):
            with self.subTest(fragment=fragment):
                s, r = _build(dict({"longitudinal_bond_channels": True}, **extra), inter)

                def _boom(*a, **k):
                    raise AssertionError("expensive work ran before the gate check")
                s._find_mu = _boom
                s._calc_green = _boom
                s._calc_chi0q = _boom
                gi = r.get_param("green")
                with self.assertRaises(ValueError) as cm:
                    s.solve(gi, "tests/rpa/output")
                self.assertIn(fragment, str(cm.exception))

    def test_external_chi0q_is_refused_before_inspection(self):
        """A supplied chi0q is refused by the gate before its dtype, shape
        or provenance are inspected (a string tensor still gets the
        gate-owned message)."""
        s, r = _build({"longitudinal_bond_channels": True})
        gi = r.get_param("green")
        gi["chi0q"] = np.array(["not", "a", "bubble"])
        with self.assertRaises(ValueError) as cm:
            s.solve(gi, "tests/rpa/output")
        self.assertIn("longitudinal_bond_channels", str(cm.exception))
        self.assertIn("chi0q", str(cm.exception))

    def test_resolved_backend_state_reaches_the_preflight(self):
        s, r = _build({"longitudinal_bond_channels": True})
        seen = {}
        real = s._longitudinal_bond_resource_preflight

        def _spy(topo, gpu_active=None):
            seen["gpu_active"] = gpu_active
            return real(topo, gpu_active=gpu_active)
        s._longitudinal_bond_resource_preflight = _spy
        gi = r.get_param("green")
        os.makedirs("tests/rpa/output", exist_ok=True)
        s.solve(gi, "tests/rpa/output")
        self.assertIs(seen["gpu_active"], False)      # CPU run: resolved, not None

    def test_spin_orbital_flag_is_read_from_the_interaction(self):
        s = self._gate()
        s.ham_info.enable_spin_orbital = True
        with self.assertRaises(ValueError) as cm:
            s._validate_longitudinal_bond_prereqs()
        self.assertIn("enable_spin_orbital", str(cm.exception))

    def test_odd_nmat_is_the_constructor_exit(self):
        """The documented exception to 'every refusal is a ValueError':
        an odd Nmat is refused by the constructor (sys.exit) before the
        gate is ever reached; the gate's own even-Nmat guard is a
        belt-and-braces check for programmatic callers."""
        with self.assertRaises(SystemExit):
            _build({"longitudinal_bond_channels": True, "Nmat": 31})
        s = self._gate()
        s.nmat = 31
        with self.assertRaises(ValueError) as cm:
            s._validate_longitudinal_bond_prereqs()
        self.assertIn("even", str(cm.exception))

    def test_spin_mode_and_external_chi0q_refused(self):
        s = self._gate()
        s.spin_mode = "spinful"
        with self.assertRaises(ValueError) as cm:
            s._validate_longitudinal_bond_prereqs()
        self.assertIn("spin", str(cm.exception))
        s = self._gate()
        s._chi0q_external = True
        with self.assertRaises(ValueError) as cm:
            s._validate_longitudinal_bond_prereqs()
        self.assertIn("chi0q", str(cm.exception))

    def test_accepts_and_returns_topology_and_split(self):
        s = self._gate()
        topo, split = s._validate_longitudinal_bond_prereqs()
        self.assertEqual(tuple(topo.coeffs), ("CoulombInter", "Hund", "Ising"))
        self.assertGreater(topo.delta_r.shape[0], 1)
        self.assertEqual(split.offsite_types, ("CoulombInter",))
        self.assertFalse(split.has_fold)

    def test_declared_but_zero_shell_is_active(self):
        s = self._gate()
        topo, _ = s._validate_longitudinal_bond_prereqs()
        s2 = self._gate()
        for k in list(s2.ham_info.param_ham["CoulombInter"]):
            s2.ham_info.param_ham["CoulombInter"][k] = 0.0
        topo2, _ = s2._validate_longitudinal_bond_prereqs()
        self.assertEqual(topo2.delta_r.shape, topo.delta_r.shape)
        self.assertTrue(all(np.all(a == 0) for a in topo2.coeffs.values()))

    def test_hund_and_ising_shells_count_and_pairlift_is_inert(self):
        s = self._gate({}, {"Hund": "offsite_hund.dat"})
        topo, split = s._validate_longitudinal_bond_prereqs()
        self.assertEqual(topo.delta_r.shape[0], 3)
        self.assertEqual(split.offsite_types, ("Hund",))
        s = self._gate({}, {"Hund": "offsite_hund.dat"})
        s.ham_info.param_ham["PairLift"] = {((1, 0, 0), (0, 1)): 0.1,
                                            ((-1, 0, 0), (1, 0)): 0.1}
        with self.assertLogs("hwave.solver.rpa", level="INFO") as cm:
            topo, split = s._validate_longitudinal_bond_prereqs()
        self.assertTrue(any("PairLift" in m for m in cm.output))
        self.assertEqual(topo.delta_r.shape[0], 3)

    def test_max_shells_truncation_of_a_declared_nonzero_shell_refused(self):
        s = self._gate({"longitudinal_bond_max_shells": 1})
        s.ham_info.param_ham["CoulombInter"][((2, 0, 0), (0, 0))] = 0.05
        s.ham_info.param_ham["CoulombInter"][((-2, 0, 0), (0, 0))] = 0.05
        with self.assertRaises(ValueError):
            s._validate_longitudinal_bond_prereqs()
        s.ham_info.param_ham["CoulombInter"][((2, 0, 0), (0, 0))] = 0.0
        s.ham_info.param_ham["CoulombInter"][((-2, 0, 0), (0, 0))] = 0.0
        topo, _ = s._validate_longitudinal_bond_prereqs()   # zero shell dropped
        self.assertEqual(topo.delta_r.shape[0], 5)          # 4x4: +-x, +-y

    def test_complex_coefficient_refused(self):
        s = self._gate()
        s.ham_info.param_ham["CoulombInter"][((1, 0, 0), (0, 1))] = 0.1 + 0.2j
        s.ham_info.param_ham["CoulombInter"][((-1, 0, 0), (1, 0))] = 0.1 - 0.2j
        with self.assertRaises(ValueError) as cm:
            s._validate_longitudinal_bond_prereqs()
        self.assertIn("complex", str(cm.exception).lower())

    def test_aggregate_coulomb_feeds_the_topology(self):
        s = self._gate({}, {"CoulombInter": "onsite_inter.dat"})
        del s.ham_info.param_ham["CoulombInter"]
        s.ham_info.param_ham["Coulomb"] = {((0, 0, 0), (0, 0)): 1.0,
                                           ((1, 0, 0), (0, 0)): 0.2,
                                           ((-1, 0, 0), (0, 0)): 0.2}
        topo, split = s._validate_longitudinal_bond_prereqs()
        self.assertEqual(topo.delta_r.shape[0], 3)
        self.assertEqual(split.offsite_types, ("CoulombInter",))


class TestGatePreflight(unittest.TestCase):

    def test_estimate_families_and_cap(self):
        s, _ = _build({"longitudinal_bond_channels": True,
                       "longitudinal_bond_memory_cap_gb": 1e6})
        s.spin_mode = "spin-free"
        topo, _ = s._validate_longitudinal_bond_prereqs()
        est = s._longitudinal_bond_resource_preflight(topo)
        nvol, P, nmat = s.lattice.nvol, s.norb ** 2, s.nmat
        ND = topo.delta_r.shape[0] * P
        U = nvol * ND * ND * 16
        self.assertEqual(est["ND"], ND)
        self.assertAlmostEqual(est["peak_solve"], 1.25 * (5 + 2) * U)
        prep = 16 * nvol * (ND ** 2 + 4 * P * (nmat + 2))
        pair = 16 * nvol * (2 * ND ** 2 + 3 * nmat * P ** 2 + 2 * nmat * P + 8 * P)
        self.assertAlmostEqual(est["bubble_bytes"], 1.25 * max(prep, pair))
        self.assertAlmostEqual(est["peak"], max(est["bubble_bytes"], est["peak_solve"]))
        self.assertEqual(est["cap_bytes"], 1e6 * 1024 ** 3)

    def test_cap_refusal_names_the_shapes(self):
        s, _ = _build({"longitudinal_bond_channels": True,
                       "longitudinal_bond_memory_cap_gb": 1e-9})
        s.spin_mode = "spin-free"
        topo, _ = s._validate_longitudinal_bond_prereqs()
        with self.assertRaises(ValueError) as cm:
            s._longitudinal_bond_resource_preflight(topo)
        msg = str(cm.exception)
        self.assertIn("longitudinal_bond_memory_cap_gb", msg)
        self.assertIn("ND", msg)
        self.assertIn("GiB", msg)

    def test_gpu_fallback_uses_the_resolved_backend_state(self):
        s, _ = _build({"longitudinal_bond_channels": True,
                       "longitudinal_bond_memory_cap_gb": 1e6})
        s.spin_mode = "spin-free"
        topo, _ = s._validate_longitudinal_bond_prereqs()
        s.use_gpu = True
        U = s._longitudinal_bond_resource_preflight(topo)["U"]
        self.assertAlmostEqual(
            s._longitudinal_bond_resource_preflight(topo, gpu_active=False)["peak_solve"],
            1.25 * (5 + 2) * U)
        self.assertAlmostEqual(
            s._longitudinal_bond_resource_preflight(topo, gpu_active=True)["peak_solve"],
            1.25 * (5 + 2 + 1) * U)
        self.assertAlmostEqual(
            s._longitudinal_bond_resource_preflight(topo)["peak_solve"],
            1.25 * (5 + 2 + 1) * U)                      # falls back to use_gpu

    def test_op_count_warning(self):
        s, _ = _build({"longitudinal_bond_channels": True,
                       "longitudinal_bond_memory_cap_gb": 1e9})
        s.spin_mode = "spin-free"
        topo, _ = s._validate_longitudinal_bond_prereqs()
        s.lattice.nvol = 10 ** 9                    # synthetic: force the warning
        with self.assertLogs("hwave.solver.rpa", level="WARNING") as cm:
            s._longitudinal_bond_resource_preflight(topo)
        self.assertTrue(any("SVD" in m for m in cm.output), cm.output)


def _lb_keys(gi):
    return sorted(k for k in gi if str(k).startswith("longitudinal_bond_"))


def _inject_zero_shell(solver):
    """A declared-but-ZERO off-site shell (+-x, orbital pair (0,1)) on top
    of an on-site-only declaration: identical physics to the plain solve,
    but the gate's topology has B=3 (spec G0)."""
    tbl = solver.ham_info.param_ham.setdefault("CoulombInter", {})
    tbl[((1, 0, 0), (0, 1))] = 0.0
    tbl[((-1, 0, 0), (1, 0))] = 0.0


class TestGatePipelineAndOutputs(unittest.TestCase):

    _OUT = "tests/rpa/output"

    def _run(self, param_extra=None, interactions=None):
        s, r = _build(dict({"longitudinal_bond_channels": True},
                           **(param_extra or {})), interactions)
        gi = r.get_param("green")
        os.makedirs(self._OUT, exist_ok=True)
        s.solve(gi, self._OUT)
        return s, gi

    def test_g0_parity_declared_but_zero_topology(self):
        """G0: with every off-site coefficient zero the collapse equals the
        plain general path's static S/C channels (FLEX general at one
        iteration, the primary reference) and the ring's same +- diff
        combination (cross-check), the bond-enlarged bubble's channel-0
        block is the plain static chi0, and the bare vertex has no
        content outside its channel-0 block."""
        s, r = _build({"longitudinal_bond_channels": True},
                      {"CoulombInter": "onsite_inter.dat"})
        _inject_zero_shell(s)
        gi = r.get_param("green")
        os.makedirs(self._OUT, exist_ok=True)
        s.solve(gi, self._OUT)
        cs = np.asarray(gi["longitudinal_bond_chiq_s_static"])
        cc = np.asarray(gi["longitudinal_bond_chiq_c_static"])
        w0 = int(s.nmat) // 2
        norb = s.norb
        chiq = np.asarray(gi["chiq"])[w0]
        same = chiq[:, :norb, :norb, :norb, :norb]
        diff = chiq[:, :norb, :norb, norb:, norb:]
        np.testing.assert_allclose(cc, same + diff, rtol=0, atol=1e-12)
        np.testing.assert_allclose(cs, same - diff, rtol=0, atol=1e-12)

        f, r2 = _build({}, {"CoulombInter": "onsite_inter.dat"}, mode="FLEX")
        _inject_zero_shell(f)
        gf = r2.get_param("green")
        f.solve(gf, self._OUT)
        np.testing.assert_allclose(cs, np.asarray(gf["chiq_s"])[w0], rtol=0, atol=1e-12)
        np.testing.assert_allclose(cc, np.asarray(gf["chiq_c"])[w0], rtol=0, atol=1e-12)

        from hwave.solver import bond_channels as bc, bubble
        topo, _ = s._validate_longitudinal_bond_prereqs()
        chi_bar = np.asarray(bubble.bond_bubble_static(
            s.green0, s.green0_tail, 1.0 / 2.0, bc.BondSetView(topo),
            spatial_shape=tuple(s.lattice.shape)))
        nd = norb * norb
        chi_s = np.asarray(gi["longitudinal_bond_chi_s"])
        self.assertEqual(chi_s.shape, chi_bar.shape)
        self.assertEqual(chi_s.shape[1], 3 * nd)
        # the plain chi0 static slice is chi_bar's channel-0 block
        np.testing.assert_allclose(
            chi_bar[:, :nd, :nd].reshape(-1, norb, norb, norb, norb),
            np.asarray(gi["chi0q"])[w0], rtol=0, atol=1e-13)
        # and the vertex has NO content outside the channel-0 block
        from hwave.solver.offsite import sc_matrices_from_split
        _, split = s._validate_longitudinal_bond_prereqs()
        S0, C0 = sc_matrices_from_split(split, bc._LONGITUDINAL_ACTIVE_TYPES, norb,
                                        *s.lattice.shape)
        for ch, W0 in (("S", S0), ("C", C0)):
            W = bc.build_sc_bond_channel(topo, W0.reshape(-1, nd, nd), ch)
            mask = np.ones(W.shape[1:], bool)
            mask[:nd, :nd] = False
            self.assertEqual(np.abs(W[:, mask]).max(), 0.0)

    def test_keys_present_when_on_absent_when_off_and_saved(self):
        s, gi = self._run()
        keys = _lb_keys(gi)
        self.assertEqual(len(keys), 16, keys)
        B = int(gi["longitudinal_bond_delta_r"].shape[0])
        nd = s.norb ** 2
        nvol = s.lattice.nvol
        self.assertEqual(gi["longitudinal_bond_chi_s"].shape, (nvol, B * nd, B * nd))
        self.assertEqual(gi["longitudinal_bond_chi_c"].shape, (nvol, B * nd, B * nd))
        self.assertEqual(gi["longitudinal_bond_chiq_s_static"].shape, (nvol,) + (s.norb,) * 4)
        self.assertEqual(gi["longitudinal_bond_reverse"].shape, (B,))
        self.assertEqual(list(gi["longitudinal_bond_types"]), ["CoulombInter", "Hund", "Ising"])
        self.assertEqual(int(gi["longitudinal_bond_schema"]), 1)
        self.assertEqual(int(gi["longitudinal_bond_max_shells"]), -1)
        self.assertEqual(tuple(gi["longitudinal_bond_spatial_shape"]), (4, 4, 1))
        self.assertEqual(str(gi["longitudinal_bond_spin_mode"]), "spin-free")
        exact = {
            "longitudinal_bond_index_order": "I = m*norb**2 + l1*norb + l2",
            "longitudinal_bond_q_convention":
                "q = 2*pi*(n_x/N_x, n_y/N_y, n_z/N_z), C-order flattened",
            "longitudinal_bond_normalization": "chi_bar = -(T/N) sum_k G G, per site",
        }
        for k, v in exact.items():
            self.assertEqual(str(gi[k]), v)
            self.assertEqual(np.asarray(gi[k]).dtype.kind, "U")
        for k in ("longitudinal_bond_delta_r", "longitudinal_bond_reverse",
                  "longitudinal_bond_spatial_shape", "longitudinal_bond_max_shells",
                  "longitudinal_bond_schema"):
            self.assertEqual(np.asarray(gi[k]).dtype, np.int64, k)
        for k in ("longitudinal_bond_cond_min_s", "longitudinal_bond_cond_min_c"):
            self.assertEqual(np.asarray(gi[k]).dtype, np.float64, k)
        for k in ("longitudinal_bond_chi_s", "longitudinal_bond_chi_c",
                  "longitudinal_bond_chiq_s_static", "longitudinal_bond_chiq_c_static"):
            self.assertEqual(np.asarray(gi[k]).dtype, np.complex128, k)
        self.assertGreater(float(gi["longitudinal_bond_cond_min_s"]), 1e-3)
        self.assertGreater(float(gi["longitudinal_bond_cond_min_c"]), 1e-3)
        s.save_results({"path_to_output": self._OUT, "chiq": "chiq_lb.npz"}, gi)
        data = np.load(os.path.join(self._OUT, "chiq_lb.npz"), allow_pickle=True)
        self.assertEqual(sorted(k for k in data.files
                                if k.startswith("longitudinal_bond_")), keys)
        self.assertIn("chiq", data.files)
        for k, v in exact.items():
            self.assertEqual(str(data[k]), v)
        for k in keys:
            np.testing.assert_array_equal(np.asarray(data[k]), np.asarray(gi[k]))

        s_off, r = _build({})
        gi2 = r.get_param("green")
        s_off.solve(gi2, self._OUT)
        self.assertEqual(_lb_keys(gi2), [])
        s_off.save_results({"path_to_output": self._OUT, "chiq": "chiq_off.npz"}, gi2)
        data = np.load(os.path.join(self._OUT, "chiq_off.npz"), allow_pickle=True)
        self.assertFalse([k for k in data.files if k.startswith("longitudinal_bond_")])

    def test_reused_green_info_carries_no_stale_keys(self):
        s, gi = self._run()
        self.assertEqual(len(_lb_keys(gi)), 16)
        s_off, _ = _build({})
        s_off.solve(gi, self._OUT)                    # gate off, same container
        self.assertEqual(_lb_keys(gi), [])
        # a container carrying a previous chi0q is the in-memory reuse
        # route (an EXTERNAL chi0q), which the gate refuses like chi0q_init
        s2, _ = _build({"longitudinal_bond_channels": True})
        with self.assertRaises(ValueError) as cm:
            s2.solve(gi, self._OUT)
        self.assertIn("chi0q", str(cm.exception))
        self.assertEqual(_lb_keys(gi), [])
        gi.pop("chi0q")
        s2.solve(gi, self._OUT)
        self.assertEqual(len(_lb_keys(gi)), 16)
        gi.pop("chi0q")
        s_bad, _ = _build({"longitudinal_bond_channels": True},
                          {"CoulombInter": "onsite_inter.dat"})
        with self.assertRaises(ValueError):           # prerequisite refusal
            s_bad.solve(gi, self._OUT)
        self.assertEqual(_lb_keys(gi), [])
        self.assertNotIn("chiq", gi)                  # dropped at entry too

    def test_nonzero_end_to_end_matches_direct_evaluation(self):
        from hwave.solver import bond_channels as bc, bubble
        from hwave.solver.offsite import sc_matrices_from_split
        s, gi = self._run({}, {"CoulombInter": "coulombinter.dat",
                               "Hund": "offsite_hund.dat"})
        topo, split = s._validate_longitudinal_bond_prereqs()
        nx, ny, nz = s.lattice.shape
        nvol, nd = s.lattice.nvol, s.norb ** 2
        S0, C0 = sc_matrices_from_split(split, bc._LONGITUDINAL_ACTIVE_TYPES,
                                        s.norb, nx, ny, nz)
        S0 = S0.reshape(nvol, nd, nd)
        C0 = C0.reshape(nvol, nd, nd)
        chi_bar = np.asarray(bubble.bond_bubble_static(
            s.green0, s.green0_tail, 1.0 / 2.0, bc.BondSetView(topo),
            spatial_shape=(nx, ny, nz)))
        chi_s, cond_s = bc.dress_channel(
            chi_bar, bc.build_sc_bond_channel(topo, S0, "S"), "spin",
            spatial_shape=(nx, ny, nz))
        chi_c, cond_c = bc.dress_channel(
            chi_bar, bc.build_sc_bond_channel(topo, C0, "C"), "charge",
            spatial_shape=(nx, ny, nz))
        np.testing.assert_allclose(np.asarray(gi["longitudinal_bond_chi_s"]), chi_s, rtol=0, atol=1e-12)
        np.testing.assert_allclose(np.asarray(gi["longitudinal_bond_chi_c"]), chi_c, rtol=0, atol=1e-12)
        np.testing.assert_allclose(
            np.asarray(gi["longitudinal_bond_chiq_s_static"]),
            chi_s[:, :nd, :nd].reshape(nvol, s.norb, s.norb, s.norb, s.norb), rtol=0, atol=1e-13)
        np.testing.assert_allclose(
            np.asarray(gi["longitudinal_bond_chiq_c_static"]),
            chi_c[:, :nd, :nd].reshape(nvol, s.norb, s.norb, s.norb, s.norb), rtol=0, atol=1e-13)
        self.assertEqual(float(gi["longitudinal_bond_cond_min_s"]), cond_s)
        self.assertEqual(float(gi["longitudinal_bond_cond_min_c"]), cond_c)
        # the Hund bond blocks make the collapse DIFFER from the plain ring
        w0 = int(s.nmat) // 2
        norb = s.norb
        chiq = np.asarray(gi["chiq"])[w0]
        same = chiq[:, :norb, :norb, :norb, :norb]
        diff = chiq[:, :norb, :norb, norb:, norb:]
        self.assertGreater(np.abs(np.asarray(gi["longitudinal_bond_chiq_s_static"]) - (same - diff)).max(), 1e-6)
        # chiq itself is the plain ring result (untouched by the gate)
        s_off, r = _build({}, {"CoulombInter": "coulombinter.dat", "Hund": "offsite_hund.dat"})
        gi2 = r.get_param("green")
        s_off.solve(gi2, self._OUT)
        np.testing.assert_array_equal(np.asarray(gi["chiq"]), np.asarray(gi2["chiq"]))

    def test_charge_conditioning_failure_publishes_nothing(self):
        """Atomic publication: a refusal in the SECOND dressing leaves no
        longitudinal_bond_* key behind, while chiq (stored earlier) stays."""
        from hwave.solver import bond_channels as bc
        s, r = _build({"longitudinal_bond_channels": True})
        gi = r.get_param("green")
        os.makedirs(self._OUT, exist_ok=True)
        real = bc.dress_channel

        def _fail_charge(chi_bar, W, channel, **kw):
            if channel == "charge":
                raise ValueError("dress_bond: the charge RPA denominator is singular (forced)")
            return real(chi_bar, W, channel, **kw)
        bc.dress_channel = _fail_charge
        try:
            with self.assertRaises(ValueError) as cm:
                s.solve(gi, self._OUT)
        finally:
            bc.dress_channel = real
        self.assertIn("charge", str(cm.exception))
        self.assertEqual(_lb_keys(gi), [])
        self.assertIn("chiq", gi)

    def test_info_line_names_the_gate_outputs(self):
        with self.assertLogs("hwave.solver.rpa", level="INFO") as cm:
            self._run()
        self.assertTrue(any("chiq is the plain" in m for m in cm.output), cm.output[-5:])

    def test_measured_peak_below_preflight_estimate(self):
        import tracemalloc
        from hwave.solver import bond_channels as bc
        from hwave.solver.offsite import sc_matrices_from_split
        s, gi = self._run()
        topo, split = s._validate_longitudinal_bond_prereqs()
        est = s._longitudinal_bond_resource_preflight(topo)
        nx, ny, nz = s.lattice.shape
        nvol, nd = s.lattice.nvol, s.norb ** 2
        S0, C0 = sc_matrices_from_split(split, bc._LONGITUDINAL_ACTIVE_TYPES,
                                        s.norb, nx, ny, nz)
        S0 = S0.reshape(nvol, nd, nd)
        C0 = C0.reshape(nvol, nd, nd)
        tracemalloc.start()
        res = s._run_longitudinal_bond_pipeline(s.green0, s.green0_tail, 0.5, topo, S0, C0)
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        self.assertLess(peak, est["peak"])
        self.assertIsInstance(res, rpa_mod.LongitudinalBondResults)


if __name__ == "__main__":
    unittest.main()
