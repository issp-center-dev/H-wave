"""The spin-orbital bond-resolved FLEX at its public entry (issue #270):
``[mode] enable_spin_orbital = true`` with ``[mode.param]
longitudinal_bond_channels = true``, one map and the converged state.

* :class:`TestEntryOneMap` -- the first map at a fixed
  chemical potential from a zero seed. Its Hartree-Fock self-energy is the
  Wick oracle of the documented operators at the exact free density; its
  fluctuation self-energy is :func:`hwave.solver.flex_so_bond.so_bond_flex_map`
  on an independent free Green function.
* :class:`TestDyson` -- the per-sector Dyson step of an S_z-conserving run
  (both layouts) and its internal-consistency refusal.
* :class:`TestCrossPathB1` -- B = 1: with on-site interactions only
  (B = 1) the bond machinery equals the on-site spin-orbital FLEX.
* :class:`TestOutputs` -- the written archives and their provenance
  (``producer_stage = "flex_so_bond"``; chi0q / chiq carry the last map's
  Green digest, green / sigma the final state's), and the hwave_sc readers
  that still refuse them.
* :class:`TestPole`, :class:`TestMapWiring`, :class:`TestNotices` -- the
  static pole check, the rank-6 layout, the resident buffers, the dressing
  guard's remedy, case-insensitive keys, and the PairHop discard notice.

Every chemical potential is the fixture's fixed ``mu``
(:func:`tests.so_bond_flex_fixtures.fixture_mu`), ``coeff_tail = 0`` unless
stated. Mutations patch module attributes looked up at call time; every
reference is computed before its patch is entered.
"""
import contextlib
import logging
import os
import tempfile
import unittest
from unittest import mock

import numpy as np

from tests import so_bond_entry_fixtures as F
from tests import so_bond_flex_fixtures as X
from tests.so_bond_oracles import build_h1, build_terms, free_density, rel, wick_hf_sigma_k

#: the array tolerance of the one-map checks (the entry's first map, B = 1)
TOL = 1e-12
#: the miss floor of a deliberate break (and of the non-vacuity checks)
MISS = 1e-3
#: the content floor of the Hartree-Fock non-vacuity check
CONTENT = 0.1
#: the floor of the spin-off-diagonal Hartree-Fock content of SOC_ALL_TYPES:
#: the fixture is fixed and its spin-off-diagonal blocks carry 1.77e-2 of the
#: max (measured), so 1e-2 keeps a 1.8x margin; CONTENT = 0.1 would need a
#: fixture with a stronger SOC
SPIN_OFF_CONTENT = 1e-2
#: run-to-run reproducibility of two identical solves (the local numpy 1.26
#: OpenBLAS build is not bitwise reproducible; numpy 2 is)
REPRO_TOL = 1e-14
#: the SCF threshold of the converged comparison: max(10 eps, 1e-8), EPS = 12
CONV_TOL = max(10 * 1e-12, 1e-8)

BETA = 1.0 / X.FLEX_T
REMEDY = (" (remedy: lower [mode.param] longitudinal_bond_cond_tol, a finite number in (0, 1), "
          "default 1.0e-3; the dressing guard is not relaxed by flex_guard_policy)")
POLE_WHERE = "FLEX longitudinal_bond_channels (spin-orbital, final map)"
POLE_CONSEQUENCE = ("The final FLEX state (its effective interaction W and self-energy) is "
                    "built on that susceptibility; it describes no physical state and is "
                    "still reported.")
OUTPUT_INFO = ("spin-orbital bond FLEX: chi0q and chiq are the channel-0 blocks of chibar_B "
               "and chi_B of the last SCF iteration; the bond-resolved arrays are written only "
               "with longitudinal_bond_output_full = true.")
SO_BOND_FLEX_REFUSAL = (
    "[mode] enable_spin_orbital = true (hwave_sc, issue #83): [eliashberg] chi0q_mode = "
    "\"flex\" with bond_channels = true is not available yet: the spin-orbital bond path "
    "builds its bond bubble from the band's Green function and reads no FLEX output. Remove "
    "chi0q_mode (or leave it unset); it is not used on this path.")
PAIRHOP_NOTICE = "PairHop declares 4 off-site term(s)"
STONER_LABEL = "RPA pole indicator max_q eig[-chibar_B Gamma_B] (spinful) (q, i nu = 0)"

#: the fixtures of the entry's first-map checks and their layouts
ENTRY_FIXTURES = ((F.SOC_ALL_TYPES, "dense"), (X.SOC_ALL_TYPES_CPH, "dense"),
                  (X.ZEEMAN_SZ_ALL, "sz_packed"))
#: the B = 1 fixtures and their layouts
B1_FIXTURES = ((F.SOC_ONSITE, "dense"), (X.ZEEMAN_EXTERN_ONSITE, "sz_packed"))


# --------------------------------------------------------------------------
# test-side references
# --------------------------------------------------------------------------

def _zeeman(fixture):
    """The Zeeman field of an Extern fixture (``H1(R = 0) = 1`` on orbital
    0 times ``coeff_extern``; ``0`` without Extern)."""
    if not fixture.extern:
        return 0.0
    if fixture.extern != (((0, 0, 0), (0, 0), 1.0),) or fixture.norb_phys != 1:
        raise ValueError("_zeeman: only the on-site unit Extern of one orbital")
    return float(fixture.coeff_extern)


def wick_hf(fixture, mu):
    """The Wick Hartree-Fock self-energy ``(nvol, nd, nd)`` (spin-block
    index) of ``fixture``'s documented operators at the exact Fermi density
    of its one-body matrix at ``mu`` (interleaved ED modes)."""
    transfer = {(tuple(R), tuple(AB)): v for R, AB, v in F._transfer_rows(fixture)}
    shape, n = fixture.shape, fixture.norb_phys
    rho = free_density(build_h1(transfer, shape, n, zeeman=_zeeman(fixture)), mu, BETA)
    return wick_hf_sigma_k(build_terms(X.interaction_tables(fixture), shape, n), rho, shape, n)


def hf_content(sigma, norb):
    """``(off_site, spin_off_diagonal)`` of a k-space ``(nvol, nd, nd)``
    self-energy, each relative to its max."""
    scale = float(np.abs(sigma).max())
    off_site = float(np.abs(sigma - sigma.mean(axis=0)[None]).max()) / scale
    spin_off = max(float(np.abs(sigma[:, :norb, norb:]).max()),
                   float(np.abs(sigma[:, norb:, :norb]).max())) / scale
    return off_site, spin_off


def independent_green(fixture):
    """``(G, mu)``: :func:`tests.so_bond_entry_fixtures.independent_green` at
    the entry's grid, ``G`` as ``(1, nmat, nvol, nd, nd)``."""
    G, mu = F.independent_green(fixture, nmat=X.FLEX_NMAT, beta=BETA)
    nd = fixture.nd
    G = np.ascontiguousarray(np.transpose(G, (5, 2, 3, 4, 0, 1))).reshape(
        X.FLEX_NMAT, fixture.nvol, nd, nd)
    return G[None], mu


def array_map(solver, G):
    """``so_bond_flex_map`` of the admission of ``solver`` on ``G`` (bubble
    and transport), with a zero tail and the solver's ``cond_tol``."""
    from hwave.solver import flex_so_bond, so_bond
    adm = solver._so_admission
    shape = tuple(int(x) for x in solver.lattice.shape)
    gamma = so_bond.so_bond_vertex(solver.ham_info.so_context.compiled, adm.topo, shape,
                                   plan=adm.plan)
    with so_bond.open_store(adm.plan, adm.sector_map) as store:
        return flex_so_bond.so_bond_flex_map(
            store, gamma, G, G, np.zeros_like(G), BETA, adm.geom, plan=adm.plan,
            cond_tol=solver.longitudinal_bond_cond_tol)


def fockless_hf(rho_r, tables, shape, norb, *, herm_tol=1e-10):
    """``flex_hf.hf_map_spin_orbital`` without the Fock term (a deliberate
    break; same signature)."""
    from hwave.solver import hartree_fock as _hf
    nvol = rho_r.shape[0]
    rho_so = np.ascontiguousarray(rho_r, dtype=np.complex128).reshape(nvol, 2, norb, 2, norb)
    out = np.zeros((nvol, 2 * norb, 2 * norb), dtype=np.complex128)
    _hf.accumulate_hf(out, rho_so, tables.inter_table, tables.spin_table, tuple(shape),
                      include_fock=False)
    return out


def p_b_identity(M, nd, *, out=None):
    """``so_bond.p_b`` replaced by the identity (a deliberate break; same
    signature, dense or packed)."""
    from hwave.solver import so_bond
    if isinstance(M, so_bond.PackedArrays):
        if out is None:
            return so_bond.PackedArrays(M.sector_map, tuple(b.copy() for b in M.blocks))
        for o, b in zip(out.blocks, M.blocks):
            o[...] = b
        return out
    if out is None:
        return np.array(M, copy=True)
    out[...] = M
    return out


def _spy(target):
    return mock.Mock(wraps=target)


def _zero_seed_spy(seen):
    """A wrapper of ``FLEX._solve_phase_b`` recording whether the seed of
    the first map is zero."""
    from hwave.solver.flex import FLEX
    orig = FLEX._solve_phase_b

    def spy(self, *args, **kwargs):
        static, fluct = self._phase_b_seed
        seen.append(bool(np.all(static == 0) and np.all(fluct == 0)))
        return orig(self, *args, **kwargs)
    return mock.patch.object(FLEX, "_solve_phase_b", spy)


def _messages(records, level=None):
    return [r.getMessage() for r in records if level is None or r.levelno == level]


# --------------------------------------------------------------------------
# The entry's first map
# --------------------------------------------------------------------------

class TestEntryOneMap(unittest.TestCase):
    """The entry's first map against independent references
    (dense SOC_ALL_TYPES and SOC_ALL_TYPES_CPH, packed ZEEMAN_SZ_ALL)."""

    def test_hf_matches_wick(self):
        static = {}
        for fx, layout in ENTRY_FIXTURES:
            with self.subTest(fixture=fx.name):
                ref = wick_hf(fx, X.fixture_mu(fx))
                seen = []
                with _zero_seed_spy(seen):
                    payload, solver = X.one_map(fx)
                self.assertEqual(seen, [True])                # zero seed
                self.assertEqual(solver._so_admission.layout, layout)
                self.assertEqual(payload["mu"], X.fixture_mu(fx))
                got = payload["static_new"][0, 0]
                self.assertTrue(np.all(np.isfinite(got)))
                dev = rel(got, ref)
                print("entry first map HF %-18s %.3e" % (fx.name, dev))
                self.assertLessEqual(dev, TOL)
                static[fx.name] = got
        off_site, spin_off = hf_content(static["SOC_ALL_TYPES"], 1)
        print("entry first map HF content: off-site %.3e spin-off-diagonal %.3e" % (off_site, spin_off))
        self.assertGreaterEqual(off_site, CONTENT)
        # spin-off-diagonal content: floor SPIN_OFF_CONTENT (see there), and
        # the two blocks are also compared on their own scale
        self.assertGreaterEqual(spin_off, SPIN_OFF_CONTENT)
        ref = wick_hf(F.SOC_ALL_TYPES, X.fixture_mu(F.SOC_ALL_TYPES))
        got = static["SOC_ALL_TYPES"]
        for blk in ((slice(0, 1), slice(1, 2)), (slice(1, 2), slice(0, 1))):
            dev = rel(got[(Ellipsis,) + blk], ref[(Ellipsis,) + blk])
            print("entry first map HF SOC_ALL_TYPES spin-off-diagonal block %.3e" % dev)
            self.assertLessEqual(dev, TOL)
        cph = rel(static["SOC_ALL_TYPES_CPH"], static["SOC_ALL_TYPES"])
        print("entry first map HF complex PairHop effect %.3e" % cph)
        self.assertGreaterEqual(cph, MISS)

    def test_fluct_matches_array_map(self):
        for fx, _ in ENTRY_FIXTURES:
            with self.subTest(fixture=fx.name):
                payload, solver = X.one_map(fx)
                G, mu = independent_green(fx)
                self.assertEqual(mu, payload["mu"])
                ref = array_map(solver, G).sigma_fluct
                got = payload["fluct_new"]
                self.assertEqual(got.shape, ref.shape)
                self.assertTrue(np.all(np.isfinite(got)))
                dev = rel(got, ref)
                print("entry first map fluct %-18s %.3e" % (fx.name, dev))
                self.assertLessEqual(dev, TOL)

    def test_mutations(self):
        from hwave.solver import flex_hf, so_bond
        fx = F.SOC_ALL_TYPES
        # Fock-less Hartree-Fock: the Wick reference first
        ref_hf = wick_hf(fx, X.fixture_mu(fx))
        spy = _spy(fockless_hf)
        self.assertEqual(spy.call_count, 0)
        with mock.patch.object(flex_hf, "hf_map_spin_orbital", spy):
            payload, _ = X.one_map(fx)
        self.assertGreater(spy.call_count, 0)
        miss = rel(payload["static_new"][0, 0], ref_hf)
        print("entry first map mutation Fock-less HF miss %.3e" % miss)
        self.assertTrue(np.isfinite(miss))
        self.assertGreaterEqual(miss, MISS)
        # P_B = identity: the array map of the unpatched code first
        payload0, solver0 = X.one_map(fx)
        G, _ = independent_green(fx)
        ref = array_map(solver0, G).sigma_fluct
        spy = _spy(p_b_identity)
        self.assertEqual(spy.call_count, 0)
        with mock.patch.object(so_bond, "p_b", spy):
            payload, _ = X.one_map(fx)
        self.assertGreater(spy.call_count, 0)
        miss = rel(payload["fluct_new"], ref)
        print("entry first map mutation P_B identity miss %.3e" % miss)
        self.assertTrue(np.isfinite(miss))
        self.assertGreaterEqual(miss, MISS)


# --------------------------------------------------------------------------
# Per-sector Dyson step
# --------------------------------------------------------------------------

def _off_sector(A, n):
    return max(float(np.abs(A[..., :n, n:]).max()), float(np.abs(A[..., n:, :n]).max()))


class TestDyson(unittest.TestCase):
    """The per-sector Dyson step: only on an S_z-conserving spin-orbital
    bond run (both layouts), with exact off-sector zeros."""

    def _spied(self, fx, **kw):
        from hwave.solver import so_bond
        calls = []
        orig = so_bond.inverse_per_sector

        def spy(A, norb_phys):
            out = orig(A, norb_phys)
            calls.append((np.array(A, copy=True), out))
            return out
        with mock.patch.object(so_bond, "inverse_per_sector", spy):
            res = X.run_flex(fx, **kw)
        return res, calls

    def test_sector_zeros(self):
        fx = X.ZEEMAN_SZ_ALL
        for layout, param in (("sz_packed", {}), ("dense", {"so_bond_layout": "dense"})):
            with self.subTest(layout=layout):
                param = dict(param, IterationMax=3)
                res, calls = self._spied(fx, param=param)
                s = res.solver
                self.assertEqual(s._so_admission.layout, layout)
                self.assertTrue(s._so_admission.sector_dyson)
                # two calls per map (G and the G of the new self-energy) + the final state
                self.assertEqual(len(calls), 2 * 3 + 1)
                n = s.norb
                ev, ew = np.asarray(s.H0_eigenvector), np.asarray(s.H0_eigenvalue)
                H0 = np.matmul(ev * ew[:, :, None, :], np.conj(ev).swapaxes(-2, -1))
                self.assertEqual(_off_sector(H0, n), 0.0)
                for A, G in calls:
                    self.assertEqual(_off_sector(A, n), 0.0)
                    self.assertEqual(_off_sector(G, n), 0.0)
                self.assertEqual(_off_sector(res.green_info["green"], n), 0.0)
        for fx, gate in ((F.SOC_ALL_TYPES, True), (X.ZEEMAN_EXTERN_ONSITE, False)):
            with self.subTest(fixture=fx.name, gate=gate):
                res, calls = self._spied(fx, gate=gate, param={"IterationMax": 1})
                self.assertEqual(calls, [])

    def test_leak_refused(self):
        _, solver = X.one_map(X.ZEEMAN_SZ_ALL)
        sigma = np.zeros((1, X.FLEX_NMAT, X.ZEEMAN_SZ_ALL.nvol, 2, 2), dtype=np.complex128)
        sigma[0, 3, 5, 0, 1] = 1e-3
        solver._so_dyson_where = "iteration 7"
        with self.assertRaises(ValueError) as cm:
            solver._calc_dressed_green(BETA, X.fixture_mu(X.ZEEMAN_SZ_ALL), sigma)
        text = str(cm.exception)
        self.assertRegex(text, r"^spin-orbital bond FLEX: the Green function inverse has an "
                               r"S_z-off-sector entry at iteration 7 \(leak \d\.\d{3}e[-+]\d+\) "
                               r"although the input conserves S_z; this is an internal "
                               r"inconsistency, please report it$")


# --------------------------------------------------------------------------
# B = 1 against the on-site path
# --------------------------------------------------------------------------

def _b1_rows(on, off):
    return {"Sigma_HF": rel(on["static_new"], off["static_new"]),
            "Sigma_fluct": rel(on["fluct_new"], off["fluct_new"]),
            "chi0q": rel(on["chi0q"], off["chi0q"]),
            "chiq": rel(on["chiq_s"], off["chiq_s"])}


class TestCrossPathB1(unittest.TestCase):
    """B = 1: the bond machinery (gate on) against the
    on-site spin-orbital FLEX (gate off)."""

    def _pair(self, fx, param=None):
        on, son = X.one_map(fx, gate=True, param=param)
        off, _ = X.one_map(fx, gate=False, param=param)
        return on, off, son

    def _check(self, fx, layout, param=None, label=""):
        on, off, son = self._pair(fx, param)
        self.assertEqual(int(son._so_admission.geom.B), 1)
        self.assertEqual(son._so_admission.layout, layout)
        for k in ("static_new", "fluct_new", "chi0q", "chiq_s"):
            self.assertEqual(on[k].shape, off[k].shape, k)
            self.assertTrue(np.all(np.isfinite(on[k])), k)
        rows = _b1_rows(on, off)
        print("B = 1 %s%-22s " % (label, fx.name)
              + "  ".join("%s %.3e" % kv for kv in rows.items()))
        for k, v in rows.items():
            self.assertLessEqual(v, TOL, k)
        return on, off

    def test_one_map(self):
        for fx, layout in B1_FIXTURES:
            with self.subTest(fixture=fx.name):
                self._check(fx, layout)

    def test_one_map_with_tail(self):
        for fx, layout in B1_FIXTURES:
            with self.subTest(fixture=fx.name):
                _, off1 = self._check(fx, layout, {"coeff_tail": 1.0}, "tail ")
                off0, _ = X.one_map(fx, gate=False)
                effect = rel(off1["chi0q"], off0["chi0q"])
                print("B = 1 tail effect %-18s %.3e" % (fx.name, effect))
                self.assertGreaterEqual(effect, MISS)

    def test_mutation(self):
        from hwave.solver import so_bond
        for fx, _ in B1_FIXTURES:
            with self.subTest(fixture=fx.name):
                off, _ = X.one_map(fx, gate=False)          # the reference first
                spy = _spy(p_b_identity)
                with mock.patch.object(so_bond, "p_b", spy):
                    on, _ = X.one_map(fx, gate=True)
                self.assertGreater(spy.call_count, 0)
                miss = rel(on["fluct_new"], off["fluct_new"])
                print("B = 1 mutation P_B identity %-18s %.3e" % (fx.name, miss))
                self.assertTrue(np.isfinite(miss))
                self.assertGreaterEqual(miss, MISS)

    def test_converged(self):
        for fx, _ in B1_FIXTURES:
            with self.subTest(fixture=fx.name):
                on = X.converge(fx, gate=True)
                off = X.converge(fx, gate=False)
                a, b = on.green_info, off.green_info
                rows = {"Sigma": rel(a["sigma"], b["sigma"]),
                        "G": rel(a["green"], b["green"]),
                        "mu": abs(a["physics"]["mu"] - b["physics"]["mu"]),
                        "NCond": abs(a["physics"]["NCond"] - b["physics"]["NCond"])}
                print("B = 1 converged %-18s (%d / %d iterations) " % (
                    fx.name, on.solver.scf_iterations, off.solver.scf_iterations)
                    + "  ".join("%s %.3e" % kv for kv in rows.items()))
                for k, v in rows.items():
                    self.assertTrue(np.isfinite(v), k)
                    self.assertLessEqual(v, CONV_TOL, k)


# --------------------------------------------------------------------------
# Outputs and provenance
# --------------------------------------------------------------------------

def _write_run(d, fixture, param, gate=True):
    """Solve and save the FLEX run of ``fixture`` into ``<d>/out`` (green,
    sigma, chi0q, chiq), as ``hwave.qlms.run`` does; returns ``(solver,
    out, records)``."""
    inp = X.write_flex_case(d, fixture, param=param, gate=gate)
    info = inp["file"]["output"]
    info.update({"green": "green", "sigma": "sigma", "chi0q": "chi0q", "chiq": "chiq"})
    out = info["path_to_output"]
    os.makedirs(out, exist_ok=True)
    solver, green_info = X.build_solver(inp)
    root = logging.getLogger()
    handler = X._Records()
    level = root.level
    root.addHandler(handler)
    root.setLevel(logging.INFO)
    try:
        solver.validate_output_paths(info, out)
        solver.solve(green_info, out)
        solver.save_results(info, green_info)
    finally:
        root.removeHandler(handler)
        root.setLevel(level)
    return solver, out, handler.records


def _load(out, name):
    with np.load(os.path.join(out, name + ".npz"), allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def _sc_input(d, fixture, out, eliashberg):
    inp = F.write_case(d, fixture, eliashberg=eliashberg,
                       param={"T": X.FLEX_T, "Nmat": X.FLEX_NMAT})
    inp["file"]["input"]["path_to_flex_output"] = out
    return inp


class TestOutputs(unittest.TestCase):
    """The archives of a spin-orbital bond FLEX run (ZEEMAN_SZ_ALL,
    IterationMax 3, Mix 0.5) and the hwave_sc readers that refuse them."""

    @staticmethod
    def _spy_maps():
        """Patch the map and its density step: every call records a copy
        of the map's Green function, and the mu and the density of the
        map's density step."""
        from hwave.solver.flex import FLEX
        maps, dens = [], []
        real_map = FLEX._phase_b_so_bond_map
        real_dens = FLEX._phase_b_density_and_hf

        def so_map(self, store, gamma_b, green_kw, *args, **kwargs):
            maps.append(np.array(green_kw, dtype=np.complex128, copy=True))
            return real_map(self, store, gamma_b, green_kw, *args, **kwargs)

        def density(self, green_kw, static, mu, *args, **kwargs):
            out = real_dens(self, green_kw, static, mu, *args, **kwargs)
            dens.append((float(mu), float(out[0].n_per_spin)))
            return out
        patches = contextlib.ExitStack()
        patches.enter_context(mock.patch.object(FLEX, "_phase_b_so_bond_map", so_map))
        patches.enter_context(mock.patch.object(FLEX, "_phase_b_density_and_hf", density))
        return patches, maps, dens

    def test_files_and_stage(self):
        from hwave.solver import fingerprint, so_bond
        fx = X.ZEEMAN_SZ_ALL
        patches, maps, _ = self._spy_maps()
        with tempfile.TemporaryDirectory() as d:
            with patches:
                solver, out, records = _write_run(d, fx, {"IterationMax": 3, "Mix": 0.5})
            self.assertEqual(len(maps), 3)
            # the digests of the Green functions the maps were built from,
            # computed here from copies of the map's argument
            orig = fingerprint.green_digest
            loop_digests = [orig(g, fingerprint.GREEN_LAYOUT) for g in maps]
            self.assertFalse(solver.scf_converged)
            files = {n: _load(out, n) for n in ("green", "sigma", "chi0q", "chiq")}
            compiled = solver.ham_info.so_context.compiled
            topo = so_bond.so_bond_topology(compiled, np.eye(3))
            topo_digest = fingerprint.resolved_topology_digest(
                np.asarray(topo.delta_r), np.asarray(topo.reverse), so_bond.SLOT_CONVENTION)
            green_digest = orig(files["green"]["green"], fingerprint.GREEN_LAYOUT)
            for name, z in files.items():
                with self.subTest(archive=name):
                    self.assertEqual(str(z["producer_stage"]), "flex_so_bond")
                    self.assertEqual(str(z["fp_capability"]), "so_bond")
                    self.assertEqual(str(z["fp_resolved_topology"]), topo_digest)
                    self.assertEqual(str(z["producer_label"]),
                                     "hwave FLEX (spin-orbital, bond channels)")
                    self.assertFalse([k for k in z if k.startswith("longitudinal_bond_")])
                    last_map = name in ("chi0q", "chiq")
                    want = loop_digests[-1] if last_map else green_digest
                    self.assertEqual(str(z["fp_green"]), want)
                    self.assertEqual(str(z["payload_kind"]),
                                     "last_map" if last_map else "final_state")
            # Mix < 1: the final (post-mixing) G is not the last map's G
            self.assertNotEqual(loop_digests[-1], green_digest)
            self.assertEqual(files["chi0q"]["chi0q"].ndim, 6)
            self.assertEqual(files["chiq"]["chiq"].ndim, 6)
            self.assertEqual(_messages(records).count(OUTPUT_INFO), 1)

            # the hwave_sc readers: the bond path refuses chi0q_mode = "flex" by
            # configuration; its bond_green reads this green.npz and refuses it
            # at the state step: this FLEX run fixed mu (T and Nmat are equal)
            import hwave.sc as sc
            green = os.path.join(out, "green.npz")
            for eli, text in (
                    ({"Bond_Channels": True, "Chi0q_Mode": "flex"}, SO_BOND_FLEX_REFUSAL),
                    ({"Bond_Channels": True, "Bond_Green": green},
                     "'{}' (written by hwave FLEX (spin-orbital, bond channels)) does not "
                     "match this run: state differs: calc_mu: file False, run True. Archives "
                     "are reused only at the identical state point; re-run hwave FLEX "
                     "(spin-orbital, bond channels) at this filling and temperature."
                     .format(green))):
                with self.subTest(hwave_sc=eli):
                    with tempfile.TemporaryDirectory() as d2:
                        with self.assertRaises(ValueError) as cm:
                            sc.calc_eliashberg(_sc_input(d2, fx, out, eli))
                    self.assertEqual(str(cm.exception), text)

    def test_block_layout_member(self):
        """Every spin-orbital bond FLEX output records the bond block layout
        of the run as a 0-d '<U16' member and nothing else of the layout
        (no sz_sector_* member); IterationMax = 0 writes it on green and
        sigma too."""
        names = ("green", "sigma", "chi0q", "chiq")
        for fx, layout in ((X.ZEEMAN_SZ_ALL, "sz_packed"), (F.SOC_ALL_TYPES, "dense")):
            for it, written in ((3, names), (0, names[:2])):
                with self.subTest(layout=layout, iteration_max=it):
                    with tempfile.TemporaryDirectory() as d:
                        _, out, _ = _write_run(d, fx, {"IterationMax": it, "Mix": 0.5})
                        files = {n: _load(out, n) for n in written}
                        for n in names[len(written):]:
                            self.assertFalse(os.path.exists(os.path.join(out, n + ".npz")))
                    for n, z in files.items():
                        self.assertEqual(str(z["block_layout"]), layout)
                        self.assertEqual(z["block_layout"].shape, ())
                        self.assertEqual(z["block_layout"].dtype, np.dtype("<U16"))
                        self.assertFalse([k for k in z if k.startswith("sz_sector_")])
                        self.assertEqual(len(z), {"green": 67, "sigma": 69,
                                                  "chi0q": 70, "chiq": 71}[n])

    def test_block_layout_absent_on_onsite(self):
        """The spin-orbital on-site outputs do not carry block_layout."""
        with tempfile.TemporaryDirectory() as d:
            solver, out, _ = _write_run(d, X.ZEEMAN_EXTERN_ONSITE, {"IterationMax": 1}, gate=False)
            self.assertFalse(getattr(solver, "_so_bond", False))
            for n in ("green", "sigma", "chi0q", "chiq"):
                self.assertNotIn("block_layout", _load(out, n))

    def test_no_map_no_output_info(self):
        """IterationMax = 0 with chi0q / chiq requested: no map ran, only
        green / sigma are written, and the channel-0 notice is not logged."""
        fx = X.ZEEMAN_SZ_ALL
        with tempfile.TemporaryDirectory() as d:
            _, out, records = _write_run(d, fx, {"IterationMax": 0})
            written = sorted(n for n in ("green", "sigma", "chi0q", "chiq")
                             if os.path.exists(os.path.join(out, n + ".npz")))
        self.assertEqual(written, ["green", "sigma"])
        self.assertNotIn(OUTPUT_INFO, _messages(records))

    def test_last_map_state_point(self):
        """chi0q / chiq carry the state point (mu, density) of the last
        map's Green function, like its digest; green / sigma carry the
        final state's. A chemical-potential search (filling), three maps,
        Mix 0.5: the last map's mu is not the final state's (both densities
        meet the filling target)."""
        from hwave.solver import fingerprint
        fx = X.ZEEMAN_SZ_ALL
        patches, maps, dens = self._spy_maps()
        with tempfile.TemporaryDirectory() as d:
            inp_param = {"IterationMax": 3, "Mix": 0.5}
            inp = X.write_flex_case(d, fx, param=inp_param)
            par = inp["mode"]["param"]
            del par["mu"]
            par["filling"] = fx.filling
            info = inp["file"]["output"]
            info.update({"green": "green", "sigma": "sigma", "chi0q": "chi0q", "chiq": "chiq"})
            out = info["path_to_output"]
            os.makedirs(out, exist_ok=True)
            with patches:
                solver, green_info = X.build_solver(inp)
                solver.validate_output_paths(info, out)
                solver.solve(green_info, out)
                solver.save_results(info, green_info)
            self.assertTrue(solver.calc_mu)
            self.assertEqual((len(maps), len(dens)), (3, 3))
            files = {n: _load(out, n) for n in ("green", "sigma", "chi0q", "chiq")}
        cellvol = solver.lattice.cellvol
        mu_map, n_map = dens[-1]
        target = solver.Ncond / cellvol
        final = files["green"]
        self.assertEqual(float(final["state_mu"]), float(solver.mu))
        self.assertNotEqual(mu_map, float(final["state_mu"]))
        for name in ("chi0q", "chiq"):
            with self.subTest(archive=name):
                z = files[name]
                self.assertEqual(str(z["payload_kind"]), "last_map")
                self.assertEqual(str(z["fp_green"]),
                                 fingerprint.green_digest(maps[-1], fingerprint.GREEN_LAYOUT))
                self.assertEqual(float(z["state_mu"]), mu_map)
                self.assertEqual(float(z["state_density"]), n_map / cellvol)
                self.assertEqual(float(z["state_target"]), target)
                self.assertEqual(float(z["state_mu_tol"]), abs(n_map / cellvol - target))
                self.assertTrue(bool(z["state_calc_mu"]))
                for k in ("state_beta", "state_nmat"):
                    self.assertEqual(z[k], final[k])
        self.assertEqual(float(files["sigma"]["state_mu"]), float(final["state_mu"]))

    def test_bond_stage_is_read_by_the_fingerprint_check(self):
        """The bond stage is no longer deferred (only the dynamic bond stage
        is), and every output of a bond FLEX run, in both block layouts,
        passes the structure check of its kind."""
        from hwave.solver import fingerprint
        self.assertEqual(set(fingerprint.DEFERRED_STAGES), {"sc_so_bond_dynamic"})
        for fx in (X.ZEEMAN_SZ_ALL, F.SOC_ALL_TYPES):
            with tempfile.TemporaryDirectory() as d:
                _, out, _ = _write_run(d, fx, {"IterationMax": 1})
                for kind in ("green", "sigma", "chi0q", "chiq"):
                    with self.subTest(fixture=fx.name, kind=kind):
                        path = os.path.join(out, kind + ".npz")
                        with np.load(path, allow_pickle=False) as z:
                            self.assertEqual(
                                fingerprint.validate_structure(z, path, kind=kind),
                                (kind, "flex_so_bond", 2))


ONSITE_BOND_TEXT = (
    "[mode] enable_spin_orbital = true (hwave_sc, issue #83): '{path}' was written by the "
    "spin-orbital FLEX with bond channels ([mode.param] longitudinal_bond_channels = true): "
    "the chi0q, chiq and green of that run describe the bond-resolved run, which [eliashberg] "
    "chi0q_mode = \"flex\" cannot use. Re-run FLEX with longitudinal_bond_channels = false "
    "for this mode.")


class TestOnsiteRefusesBond(unittest.TestCase):
    """The on-site spin-orbital hwave_sc (chi0q_mode = "flex") refuses the
    outputs of the bond FLEX with a message naming the file, whichever of the
    four archives is the bond one."""

    NAMES = ("chiq", "chi0q", "green", "sigma")

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        root = cls._tmp.name
        cls.dirs = {}
        for key, fx, gate in (("bond_b1", X.ZEEMAN_EXTERN_ONSITE, True),
                              ("bond_b2", X.ZEEMAN_SZ_ALL, True),
                              ("onsite", X.ZEEMAN_EXTERN_ONSITE, False)):
            d = os.path.join(root, key + "_case")
            os.makedirs(d)
            _, out, _ = _write_run(d, fx, {"IterationMax": 1}, gate=gate)
            cls.dirs[key] = out

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def _mixed(self, bond_names, bond="bond_b1"):
        """A directory of the on-site archives with ``bond_names`` taken from
        the bond run."""
        import shutil
        d = tempfile.mkdtemp(dir=self._tmp.name)
        for n in self.NAMES:
            src = self.dirs[bond if n in bond_names else "onsite"]
            shutil.copy(os.path.join(src, n + ".npz"), os.path.join(d, n + ".npz"))
        return d

    def _refusal(self, out, frequency):
        import hwave.sc as sc
        with tempfile.TemporaryDirectory() as d2:
            inp = _sc_input(d2, X.ZEEMAN_EXTERN_ONSITE, out,
                            {"Bond_Channels": False, "Chi0q_Mode": "flex",
                             "Frequency": frequency})
            with self.assertRaises(ValueError) as cm:
                sc.calc_eliashberg(inp)
        return str(cm.exception)

    def test_onsite_hwave_sc_refuses_bond_outputs(self):
        for key in ("bond_b1", "bond_b2"):
            for frequency in ("static", "dynamic"):
                with self.subTest(run=key, frequency=frequency):
                    out = self.dirs[key]
                    self.assertEqual(
                        self._refusal(out, frequency),
                        ONSITE_BOND_TEXT.format(path=os.path.join(out, "chiq.npz")))

    def test_each_bond_archive_is_named(self):
        for name in self.NAMES:
            out = self._mixed((name,))
            for frequency in ("static", "dynamic"):
                with self.subTest(archive=name, frequency=frequency):
                    self.assertEqual(
                        self._refusal(out, frequency),
                        ONSITE_BOND_TEXT.format(path=os.path.join(out, name + ".npz")))

    def test_onsite_outputs_pass_the_peek(self):
        """A directory without any bond archive is not refused by the bond
        text: the existing state-point check decides (the FLEX run fixed mu,
        the hwave_sc input asks for a filling)."""
        import hwave.sc as sc
        from hwave.solver import fingerprint
        with tempfile.TemporaryDirectory() as d2:
            inp = _sc_input(d2, X.ZEEMAN_EXTERN_ONSITE, self.dirs["onsite"],
                            {"Bond_Channels": False, "Chi0q_Mode": "flex",
                             "Frequency": "static"})
            with self.assertRaises(fingerprint.SOArchiveError) as cm:
                sc.calc_eliashberg(inp)
        self.assertIn("state differs", str(cm.exception))
        self.assertNotIn("bond channels", str(cm.exception))

    def test_archive_without_producer_stage_falls_through(self):
        """An archive with no producer_stage member is not the bond refusal:
        the existing checks decide."""
        out = self._mixed(())
        for name in self.NAMES:
            with np.load(os.path.join(out, name + ".npz"), allow_pickle=False) as z:
                members = {k: z[k] for k in z.files if k != "producer_stage"}
            np.savez(os.path.join(out, name + ".npz"), **members)
        msg = self._refusal(out, "static")
        self.assertNotIn("bond channels", msg)
        self.assertIn("member producer_stage is missing", msg)


# --------------------------------------------------------------------------
# Static pole, layout, buffers, guard, keys, notices
# --------------------------------------------------------------------------

class TestFilling(unittest.TestCase):
    """A chemical-potential search on the spin-orbital bond path: the
    fixed-mu fixtures with ``filling`` instead of ``mu``, three maps."""

    def test_mu_search_smoke(self):
        for fx in (X.ZEEMAN_SZ_ALL, F.SOC_ALL_TYPES):
            with self.subTest(fixture=fx.name):
                with tempfile.TemporaryDirectory() as d:
                    inp = X.write_flex_case(d, fx, param={"IterationMax": 3})
                    par = inp["mode"]["param"]
                    del par["mu"]
                    par["filling"] = fx.filling
                    res = X.run_input(inp)
                s = res.solver
                self.assertTrue(s.calc_mu)
                self.assertEqual(len(res.payloads), 3)
                phys = res.green_info["physics"]
                target = fx.filling * fx.nd * fx.nvol
                print("filling smoke %-14s NCond %.12g (target %.12g) mu %.6f" % (
                    fx.name, phys["NCond"], target, phys["mu"]))
                self.assertTrue(np.isfinite(phys["mu"]))
                self.assertLessEqual(abs(phys["NCond"] - target), 1e-8 * target)
                for p in res.payloads:
                    self.assertTrue(np.all(np.isfinite(p["fluct_new"])))


class TestPole(unittest.TestCase):
    """The static RPA pole indicator of every map, the check of
    the final map."""

    def test_final_check_once(self):
        from hwave.solver import bond_channels
        fx = X.ZEEMAN_SZ_ALL
        with tempfile.TemporaryDirectory() as d:
            inp = X.write_flex_case(d, fx, param={"IterationMax": 3})
            solver, green_info = X.build_solver(inp)
            seen = []
            solver._iteration_hook = lambda p: seen.append(
                (p["iteration"], solver._stoner_last, solver._bond_last.pole[0]))
            spy = _spy(bond_channels.check_rpa_pole)
            with mock.patch.object(bond_channels, "check_rpa_pole", spy), \
                    self.assertLogs(level="DEBUG") as cm:
                solver.solve(green_info, inp["file"]["output"]["path_to_output"])
        self.assertEqual(spy.call_count, 1)
        args, kwargs = spy.call_args
        self.assertEqual(kwargs["where"], POLE_WHERE)
        self.assertEqual(kwargs["consequence"], POLE_CONSEQUENCE)
        self.assertEqual(args[0], {"spinful": seen[-1][2]})
        self.assertEqual([s[0] for s in seen], [1, 2, 3])
        for _, stoner, entry in seen:
            self.assertEqual(stoner, {"spinful": entry})
        text = "\n".join(cm.output)
        self.assertNotIn("Eliashberg", text)
        self.assertNotIn("pairing", text)
        # one line per map (INFO on the first and the last, DEBUG between)
        self.assertEqual(sum(STONER_LABEL in line for line in cm.output), 3)


class TestMapWiring(unittest.TestCase):

    def test_chi6_layout(self):
        """CX-S: the rank-6 reading of the channel-0 blocks."""
        from hwave.solver import flex_so_bond
        _, solver = X.one_map(F.SOC_ALL_TYPES)
        nd, nmat, nvol = 2, X.FLEX_NMAT, F.SOC_ALL_TYPES.nvol
        P = nd * nd
        l, q, I, J = np.meshgrid(np.arange(nmat), np.arange(nvol), np.arange(P), np.arange(P),
                                 indexing="ij")
        chi4 = (1e6 * l + 1e4 * q + 100 * I + J).astype(np.complex128)
        fake = flex_so_bond.SOBondFlexMap(np.zeros((1, nmat, nvol, nd, nd), complex),
                                          chi4, chi4 + 1, 0.5, None, None,
                                          ((0.1, (0.0, 0.0, 0.0), 0.1), None, None))
        with mock.patch.object(flex_so_bond, "so_bond_flex_map", return_value=fake) as m:
            chi0q, chiq, _ = solver._phase_b_so_bond_map(None, None, None, None, BETA)
        self.assertEqual(m.call_count, 1)
        self.assertEqual(chi0q.shape, (nmat, nvol, nd, nd, nd, nd))
        self.assertTrue(np.shares_memory(chi0q, chi4))
        for a in range(nd):
            for c in range(nd):
                for b in range(nd):
                    for dd in range(nd):
                        np.testing.assert_array_equal(chi0q[:, :, a, c, b, dd],
                                                      chi4[:, :, a * nd + c, b * nd + dd])
                        np.testing.assert_array_equal(chiq[:, :, a, c, b, dd],
                                                      chi4[:, :, a * nd + c, b * nd + dd] + 1)
        self.assertEqual(solver._stoner_last, {"spinful": (0.1, (0.0, 0.0, 0.0), 0.1)})

    def test_no_resident_copies(self):
        """G1: the published chi0q / chiq are views of the two resident
        buffers of the last map."""
        res = X.run_flex(X.ZEEMAN_SZ_ALL, param={"IterationMax": 2})
        last = res.solver._bond_last
        self.assertTrue(np.shares_memory(res.green_info["chi0q"], last.chi0_c0))
        self.assertTrue(np.shares_memory(res.green_info["chiq"], last.chi_c0))
        self.assertIs(res.green_info["chi0q"].base, last.chi0_c0)

    def test_cond_tol_threaded(self):
        from hwave.solver import bond_channels
        with self.assertRaises(bond_channels.BondConditioningError) as cm:
            X.one_map(F.SOC_ALL_TYPES, param={"longitudinal_bond_cond_tol": 0.99})
        self.assertTrue(str(cm.exception).endswith(REMEDY), str(cm.exception))
        self.assertEqual(cm.exception.cond_tol, 0.99)
        self.assertEqual(cm.exception.channel, "spinful")

    def test_case_insensitive_entry(self):
        def run(spell):
            with tempfile.TemporaryDirectory() as d:
                inp = X.write_flex_case(d, X.ZEEMAN_SZ_ALL)
                par = inp["mode"]["param"]
                del par["longitudinal_bond_channels"]
                keys = (("longitudinal_bond_channels", True), ("so_bond_layout", "dense"),
                        ("longitudinal_bond_memory_cap_gb", 4))
                for k, v in keys:
                    if spell:
                        k = {"longitudinal_bond_channels": "Longitudinal_Bond_Channels",
                             "so_bond_layout": "SO_Bond_Layout",
                             "longitudinal_bond_memory_cap_gb":
                                 "Longitudinal_Bond_Memory_Cap_GB"}[k]
                        v = "Dense" if v == "dense" else v
                    par[k] = v
                res = X.run_input(inp)
            self.assertEqual(res.solver._so_admission.layout, "dense")
            return res.payloads[0]
        lower, mixed = run(False), run(True)
        for k in ("static_new", "fluct_new", "chi0q", "chiq_s"):
            self.assertLessEqual(float(np.abs(lower[k] - mixed[k]).max()),
                                 REPRO_TOL * float(np.abs(lower[k]).max()), k)


class TestNotices(unittest.TestCase):
    """RPA's off-site PairHop discard WARNING: absent where the spin-orbital
    bond vertex includes the rows, present where it was before."""

    def _construct(self, fx, gate, so=True):
        with tempfile.TemporaryDirectory() as d:
            inp = X.write_flex_case(d, fx, gate=gate, so=so)
            with self.assertLogs(level="INFO") as cm:
                X.build_solver(inp)
        return [r.getMessage() for r in cm.records if r.levelno == logging.WARNING]

    def test_pairhop_discard_notice(self):
        for fx in (X.ZEEMAN_SZ_ALL, X.SOC_ALL_TYPES_CPH, F.SOC_ALL_TYPES):
            with self.subTest(fixture=fx.name, gate=True):
                self.assertFalse([m for m in self._construct(fx, True)
                                  if m.startswith(PAIRHOP_NOTICE)])
            with self.subTest(fixture=fx.name, gate=False):
                self.assertEqual(len([m for m in self._construct(fx, False)
                                      if m.startswith(PAIRHOP_NOTICE)]), 1)
        # non-spin-orbital FLEX (gate on and off): unchanged
        for gate in (True, False):
            with self.subTest(fixture="SU2_4X4 non-SO", gate=gate):
                self.assertEqual(len([m for m in self._construct(F.SU2_4X4, gate, so=False)
                                      if m.startswith(PAIRHOP_NOTICE)]), 1)


if __name__ == "__main__":
    unittest.main()
