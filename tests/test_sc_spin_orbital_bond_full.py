"""The heavy (opt-in) tests of hwave_sc in spin-orbital mode with
bond-resolved interaction channels (issue #270): public-entry runs whose
cost does not fit the fast gate. The fast module is
``tests/test_sc_spin_orbital_bond.py``."""

import contextlib
import tempfile
import unittest
from unittest import mock

import numpy as np

import hwave.sc as sc

from tests import so_bond_entry_fixtures as F
from tests import so_run_compare
# the module object only: its test classes stay out of this module's namespace
from tests import test_sc_spin_orbital_bond as fast
from tests.heavy_tests import heavy

#: every solve on its dense branch: the dense layout diagonalizes its
#: antisymmetric space completely (num_eigenvalues >= dim - 1) and every
#: pair-S_z sector is small enough for the dense sector solve, so the two
#: layouts compare their leaders without an ARPACK selection in between
_FULL_SPECTRUM = {"Num_Eigenvalues": 40}


class TestAdmission(unittest.TestCase):
    """The traced phases of the S_z-packed power iteration at Nmat 32 within
    their rows (the traced and derived rows are those of the fast module's
    ``TestAdmission``; its ``test_phase_rows_iteration_smoke`` makes the
    same check at Nmat 8)."""

    @heavy
    def test_phase_rows_iteration(self):
        p = fast._prepare_admission(F.ZEEMAN_EXTERN, eliashberg={"Solver_Mode": "iteration",
                                                                 "Max_Iter": 1000})
        # unset: the iteration leg resolves it to "auto", the rows take None
        self.assertIsNone(p.spectral_shift)
        adm = fast.check_phase_rows(self, p, "iteration")
        self.assertEqual((adm.layout, adm.plan.nmat), ("sz_packed", 32))


class TestLayoutSelectionThroughEntry(unittest.TestCase):
    """An S_z-conserving input through the public entry: the automatic
    layout is the S_z-packed one, a forced dense layout gives the same
    leading eigenvalue and gap."""

    @heavy
    def test_layout_selection_through_entry(self):
        auto = F.run_case(F.ZEEMAN_EXTERN, eliashberg=_FULL_SPECTRUM)
        dense = F.run_case(F.ZEEMAN_EXTERN, layout="dense", eliashberg=_FULL_SPECTRUM)
        a, d = auto.gap_npz, dense.gap_npz
        self.assertEqual(str(a["block_layout"]), "sz_packed")
        self.assertEqual(str(d["block_layout"]), "dense")
        self.assertEqual((str(a["gap_status"]), str(d["gap_status"])),
                         ("certified", "certified"))
        # the solves behind the two leaders: dense branches only
        sel = {r: [l.split(" = ", 1)[1] for l in run.eigenvalue_dat.splitlines()
                   if l.startswith("# selection = ")]
               for r, run in (("auto", auto), ("dense", dense))}
        self.assertEqual(sel["dense"], ["dense-LR"])
        for part in sel["auto"][0].split(";"):
            self.assertIn(part.split(":")[1], ("dense-LR", "null-operator", "empty"), part)
        self.assertTrue(np.isfinite(auto.value))
        self.assertLessEqual(abs(auto.value - dense.value), 1e-10 * max(1.0, abs(dense.value)))
        # a degenerate leader is compared through the leading spans (both
        # listings are ordered by descending real part)
        self.assertTrue(so_run_compare.same_leading_vector(
            d["eigenvalues"], a["gap"], d["gap"], vecs_a=a["eigenvectors"],
            vecs_b=d["eigenvectors"]))


_NAN = float("nan")
#: the pair-S_z 0 sector of the Zeeman fixture (the leader's) is the only
#: sector of this antisymmetric dimension; the +1 and -1 sectors have 6
_DIM_PP0 = fast.ZEEMAN_DIM_PP0
_REAL_DENSE_EIGEN = sc._so_dense_reduced_eigen
_eigen_fails = fast.eigen_fails
_REAL_SOLVE_LEADING = sc._solve_leading


def _last_iterate(dim):
    rng = np.random.default_rng(dim)
    v = rng.standard_normal(dim) + 1j * rng.standard_normal(dim)
    return v / np.linalg.norm(v)


def _residual_fails(dim):
    """The second pair of the ``dim``-dimensional sector moved by 1e-6: its
    residual check fails (the sector is failed, "residual-failed")."""
    def fake(red, d, k, seed_vec=None):
        vals, vecs, M = _REAL_DENSE_EIGEN(red, d, k, seed_vec)
        if d == dim:
            vals = vals.copy()
            vals[1] += 1e-6
        return vals, vecs, M
    return mock.patch.object(sc, "_so_dense_reduced_eigen", fake)


def _leading_fails(dim, modes):
    """``_solve_leading`` on a ``dim``-dimensional operator: the ARPACK solve
    returns the not-converged sentinel ("arnoldi" in ``modes``), the power
    iteration ends unconverged at a fixed vector ("iteration" in
    ``modes``)."""
    def fake(make_operator, vec_size, solver_mode, *a, **kw):
        if vec_size == dim and solver_mode in modes:
            if solver_mode == "arnoldi":
                k = kw.get("num_eigenvalues", 10)
                w = np.eye(vec_size, k, dtype=complex)
                v = np.full(k, np.nan + 0j)
                return v[0], w[:, 0], {"eigenvalues": v, "eigenvectors": w,
                                       "selection": "not-converged"}
            return _NAN, _last_iterate(vec_size), {"converged": False, "n_iter": 7}
        return _REAL_SOLVE_LEADING(make_operator, vec_size, solver_mode, *a, **kw)
    return mock.patch.object(sc, "_solve_leading", fake)


def _dense_dim():
    """The antisymmetric dimension of the Zeeman fixture's full gap."""
    fx = F.ZEEMAN_EXTERN
    return sc._SOAntisymmetricBasis((fx.nd, fx.nd) + tuple(fx.shape)).dim


#: the eigen leg's selection of each forced row
_EIGEN_SELECTION = {
    ("2", "auto"): ("+1:null-operator;0:null-operator;-1:null-operator",),
    ("3", "auto"): ("+1:dense-LR;0:not-converged;-1:dense-LR",),
    ("4", "auto"): ("+1:dense-LR;0:residual-failed;-1:dense-LR",),
    ("5", "auto"): ("+1:empty;0:empty;-1:empty",),
    ("6", "auto"): ("+1:empty;0:residual-failed;-1:empty",),
    ("7", "auto"): ("+1:null-operator;0:dense-LR;-1:dense-LR",),
    ("3", "dense"): ("not-converged",),
    ("7", "dense"): ("LM", "LR", "LR_retry"),
}


class TestReportingThroughEntry(unittest.TestCase):
    """Solver modes and failures through the public entry (the fast module
    keeps one cell). The rows of the sector solvers' state table are forced
    on the Zeeman fixture by patching its sector operators or solvers:

    * 2 -- every sector's kernel exactly zero;
    * 3 -- the leader's sector (pair S_z 0) not converged, in both legs;
    * 4 -- a residual failure in that sector (eigen leg only);
    * 5 -- every sector's gap space empty;
    * 6 -- no certified pair: that sector residual-failed (and its iteration
      unconverged), the other two empty;
    * 7 -- a null sector (+1) next to converged ones.

    The dense layout (forced) has no sectors: its row 3 is the
    not-converged sentinel of the whole-gap ARPACK solve and an unconverged
    iteration, its row 7 the unpatched converged run.

    Each cell runs ``"eigenvalue"``, ``"iteration"`` and ``"both"`` and
    checks ``reported``, ``gap_status``, ``pp_sector``, ``eigen_converged``,
    ``iteration_converged`` and the exit status."""

    @heavy
    def test_state_rows(self):
        lam = F.run_case(F.ZEEMAN_EXTERN).value
        self.assertTrue(np.isfinite(lam))
        nd = _dense_dim()
        # (row, layout, patches(mode)) -> {mode: (reported, gap_status, pp_sector,
        #  eigen_converged, iteration_converged, exit status)}; "pm" = the
        #  pair-S_z sector (+1 or -1) of the first listed pair, checked apart
        cells = (
            ("2", "auto", lambda: (fast.patched_sectors(null=(1, 0, -1)),), {
                "eigenvalue": (0.0, "certified", 1, True, None, 0),
                "iteration": (0.0, "certified", 1, True, None, 0),
                "both": (0.0, "certified", 1, True, True, 0)}),
            ("3", "auto", lambda: (_eigen_fails(_DIM_PP0),
                                   _leading_fails(_DIM_PP0, ("iteration",))), {
                "eigenvalue": (_NAN, "uncertified", "pm", False, None, 3),
                "iteration": (_NAN, "uncertified", 0, False, None, 3),
                "both": (_NAN, "uncertified", 0, False, False, 3)}),
            ("4", "auto", lambda: (_residual_fails(_DIM_PP0),), {
                "eigenvalue": (_NAN, "uncertified", "pm", False, None, 3),
                "iteration": (lam, "certified", 0, True, None, 0),
                "both": (_NAN, "uncertified", 0, False, True, 3)}),
            ("5", "auto", lambda: (fast.patched_sectors(empty=(1, 0, -1)),), {
                "eigenvalue": (_NAN, "none", 2, False, None, 3),
                "iteration": (_NAN, "none", 2, False, None, 3),
                "both": (_NAN, "none", 2, False, False, 3)}),
            ("6", "auto", lambda: (fast.patched_sectors(empty=(1, -1)),
                                   _residual_fails(_DIM_PP0),
                                   _leading_fails(_DIM_PP0, ("iteration",))), {
                "eigenvalue": (_NAN, "none", 2, False, None, 3),
                "iteration": (_NAN, "uncertified", 0, False, None, 3),
                "both": (_NAN, "uncertified", 0, False, False, 3)}),
            ("7", "auto", lambda: (fast.patched_sectors(null=(1,)),), {
                "eigenvalue": (lam, "certified", 0, True, None, 0),
                "iteration": (lam, "certified", 0, True, None, 0),
                "both": (lam, "certified", 0, True, True, 0)}),
            ("3", "dense", lambda: (_leading_fails(nd, ("arnoldi", "iteration")),), {
                "eigenvalue": (_NAN, "none", 2, False, None, 3),
                "iteration": (_NAN, "uncertified", 2, False, None, 3),
                "both": (_NAN, "uncertified", 2, False, False, 3)}),
            ("7", "dense", lambda: (), {
                "eigenvalue": (lam, "certified", 2, True, None, 0),
                "iteration": (lam, "certified", 2, True, None, 0),
                "both": (lam, "certified", 2, True, True, 0)}),
        )
        for row, layout, patches, wants in cells:
            for mode, want in wants.items():
                with self.subTest(row=row, layout=layout, mode=mode):
                    run, keys, status = fast.reporting_cell(
                        self, F.ZEEMAN_EXTERN, mode, layout=layout, patches=patches())
                    self.assertEqual(str(run.gap_npz["block_layout"]),
                                     "sz_packed" if layout == "auto" else "dense")
                    if want[2] == "pm":
                        self.assertIn(int(run.gap_npz["pp_sector"]), (1, -1))
                        want = want[:2] + (int(run.gap_npz["pp_sector"]),) + want[3:]
                    fast.check_reporting_cell(self, run, keys, status, want)
                    # the forced state reached the solve
                    if mode != "iteration":
                        self.assertIn(keys["selection"], _EIGEN_SELECTION[(row, layout)])
                    listed = mode != "iteration" and row not in ("5", "6")
                    # the empty gap space WARNING: exactly when every sector is empty
                    if want[1] == "none":
                        self.assertEqual(fast.empty_gap_warnings(run), int(row == "5"))
                    self.assertEqual("eigenvalues" in run.gap_npz, listed)


class TestZeroFieldThroughEntry(unittest.TestCase):
    """Gate A2-0S3b through the public entry (the fast module's
    ``TestZeroField`` runs it at operator level on a 4x2 mesh): 4x4, on-site
    rows only (B = 1) and nearest-neighbour bonds (B = 5), one and two
    physical orbitals, both layouts, five fields each, complete listings;
    and one field through the Extern and through the Transfer route, which
    give the same listing and eigenvectors."""

    @heavy
    def test_zero_field_limit_through_entry(self):
        from hwave.solver import so_bond
        cases = ((fast.ZF_4X4_B1, 1), (fast.ZF_4X4, 5), (fast.ZF_4X4_NORB2_B1, 1),
                 (fast.ZF_4X4_NORB2, 5))
        for fx, B in cases:
            with self.subTest(fixture=fx.name):
                d = tempfile.TemporaryDirectory(prefix="sc_so_bond_zf_")
                self.addCleanup(d.cleanup)
                with fast._quiet_logs():
                    p = F.prepare(fx, d.name)
                self.assertEqual(so_bond.slot_geometry(p.so_run.bond.topo, p.shape).B, B)
                fast.entry_zero_field_gate(self, fx)
        # the field through the two routes, one field, both layouts
        h = fast.ZF_H[3]
        on_transfer = fast.ZF_4X4._replace(name="ZF_4X4_TRANSFER", extern=())
        dim = sc._SOAntisymmetricBasis((2, 2) + tuple(fast.ZF_4X4.shape)).dim
        for layout in ("dense", "auto"):
            with self.subTest(route_layout=layout):
                runs = [F.run_case(fast._at_field(route, h), layout=layout,
                                   eliashberg={"Num_Eigenvalues": dim + 2})
                        for route in (fast.ZF_4X4, on_transfer)]
                (ga, gb) = (r.gap_npz for r in runs)
                self.assertTrue(bool(ga["eigen_converged"]) and bool(gb["eigen_converged"]))
                self.assertEqual(float(ga["prov_coeff_extern"]), h)
                self.assertFalse(bool(gb["prov_extern_present"]))
                self.assertEqual(len(ga["eigenvalues"]), dim)
                self.assertTrue(so_run_compare.eigenvalues_close(ga["eigenvalues"],
                                                                 gb["eigenvalues"]))
                self.assertTrue(so_run_compare.same_eigenvectors(
                    ga["eigenvalues"], ga["eigenvectors"], gb["eigenvectors"]))
                self.assertTrue(so_run_compare.same_leading_vector(
                    ga["eigenvalues"], ga["gap"], gb["gap"], vecs_a=ga["eigenvectors"],
                    vecs_b=gb["eigenvectors"]))


if __name__ == "__main__":
    unittest.main()
