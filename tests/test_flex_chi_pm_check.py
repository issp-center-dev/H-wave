#!/usr/bin/env python3

"""The undocumented ``output_chi_pm`` verification switch.

FLEX solves the LONGITUDINAL spin channel from the Kuroki :math:`S` matrix. With
this switch it additionally solves the TRANSVERSE (spin-flip) channel from the
crossed (Fock-exchange) Hartree vertex. For a paramagnetic system SU(2) symmetry
forces :math:`\\chi^{+-} = \\chi^{zz}`, so the two must agree -- and because they
are assembled from different vertices, that agreement is a real check on the
vertex construction rather than a restatement of one computation.

The switch is deliberately absent from the manual: it produces a diagnostic
output, not physics, and is redundant by symmetry.

Tests must be run from the repository root.
"""

import os
import tempfile
import unittest

import numpy as np


def _fixture(dirpath):
    os.makedirs(dirpath, exist_ok=True)
    geom = ("  1.000000000000   0.000000000000   0.000000000000\n"
            "  0.000000000000   1.000000000000   0.000000000000\n"
            "  0.000000000000   0.000000000000   1.000000000000\n"
            "2\n"
            "    0.000000000000000e+00     0.000000000000000e+00     0.000000000000000e+00\n"
            "    0.500000000000000e+00     0.500000000000000e+00     0.000000000000000e+00\n")
    transfer = ("Transfer\n2\n9\n 1 1 1 1 1 1 1 1 1\n"
                "   0    1    0    1    1  1.0 0.0\n"
                "   0   -1    0    1    1  1.0 0.0\n"
                "   0    1    0    2    2  1.0 0.0\n"
                "   0   -1    0    2    2  1.0 0.0\n"
                "   0    0    0    1    2  0.5 0.0\n"
                "  -1    0    0    1    2  0.5 0.0\n"
                "   0    0    0    2    1  0.5 0.0\n"
                "   1    0    0    2    1  0.5 0.0\n")
    intra = ("CoulombIntra\n2\n1\n 1\n"
             "   0    0    0    1    1   2.000000000000   0.000000000000\n"
             "   0    0    0    2    2   2.000000000000   0.000000000000\n")
    for name, body in (("geom.dat", geom), ("transfer.dat", transfer),
                       ("coulombintra.dat", intra)):
        with open(os.path.join(dirpath, name), "w") as f:
            f.write(body)


def _solve_with(dirpath, scheme, switch, extra_inter=None, extra_param=None,
                green_info=None, write_fixture=True):
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as solver_flex

    if write_fixture:
        _fixture(dirpath)
    inter = {'path_to_input': dirpath, 'Geometry': 'geom.dat',
             'Transfer': 'transfer.dat',
             'CoulombIntra': 'coulombintra.dat'}
    inter.update(extra_inter or {})
    info_input = {'path_to_input': dirpath, 'interaction': inter}
    ham = read_input_k.QLMSkInput(info_input).get_param("ham")
    param = {'T': 0.05, 'CellShape': [4, 4, 1],
             'SubShape': [1, 1, 1], 'Nmat': 16, 'filling': 0.75,
             'IterationMax': 3, 'Mix': 0.5, 'EPS': 6,
             'output_chi_pm': switch}
    param.update(extra_param or {})
    solver = solver_flex.FLEX(ham, {}, {'mode': 'FLEX', 'param': param,
                                        'calc_scheme': scheme})
    if green_info is None:
        green_info = read_input_k.QLMSkInput(info_input).get_param("green")
    solver.solve(green_info, dirpath)
    return green_info


def _solve(dirpath, scheme, switch, extra_param=None, green_info=None):
    return _solve_with(dirpath, scheme, switch, None, extra_param, green_info)


class TestChiPmCheck(unittest.TestCase):

    def test_off_by_default_produces_nothing(self):
        with tempfile.TemporaryDirectory() as d:
            gi = _solve(d, "general", False)
        self.assertNotIn("chiq_pm", gi)

    def test_transverse_equals_longitudinal_for_paramagnetic(self):
        """The point of the switch: two different vertex constructions must
        give the same spin susceptibility when SU(2) holds."""
        with tempfile.TemporaryDirectory() as d:
            gi = _solve(d, "general", True)
        self.assertIn("chiq_pm", gi)
        pm = np.asarray(gi["chiq_pm"])
        cs = np.asarray(gi["chiq_s"])
        self.assertEqual(pm.shape, cs.shape)
        scale = float(np.abs(cs).max())
        self.assertGreater(scale, 0.0, "fixture produced a null chi_s")
        np.testing.assert_allclose(
            pm, cs, rtol=1e-9, atol=1e-12,
            err_msg="chi^{+-} must equal chi^{zz} for a paramagnetic run; a "
                    "mismatch means the crossed transverse vertex and the "
                    "Kuroki S matrix disagree")

    def test_reduced_scheme_is_skipped_with_a_reason(self):
        """The transverse channel needs the rank-4 orbital chi0q that only the
        general scheme produces; on the reduced path chi_s already IS the
        isotropic spin susceptibility."""
        import logging
        records = []
        handler = logging.Handler()
        handler.emit = records.append
        lg = logging.getLogger("hwave.solver.flex")
        lg.addHandler(handler)
        try:
            with tempfile.TemporaryDirectory() as d:
                gi = _solve(d, "reduced", True)
        finally:
            lg.removeHandler(handler)
        self.assertNotIn("chiq_pm", gi)
        msgs = [r.getMessage() for r in records
                if r.levelno >= logging.WARNING
                and "output_chi_pm" in r.getMessage()]
        self.assertEqual(len(msgs), 1)
        self.assertIn("calc_scheme='general'", msgs[0])



    def test_ir_axis_agrees_too(self):
        """The transverse solve must happen on the same frequency axis chi_s was
        dressed on. Densifying first and dressing after would not commute, and
        would show up here as a spurious mismatch."""
        try:
            import sparse_ir  # noqa: F401
        except ImportError:
            self.skipTest("sparse_ir not installed")
        with tempfile.TemporaryDirectory() as d:
            gi = _solve(d, "general", True,
                        {'matsubara_basis': 'ir', 'ir_wmax': 20.0})
        self.assertIn("chiq_pm", gi)
        np.testing.assert_allclose(
            np.asarray(gi["chiq_pm"]), np.asarray(gi["chiq_s"]),
            rtol=1e-9, atol=1e-12,
            err_msg="chi^{+-} must equal chi^{zz} on the IR axis as well")

    def test_non_su2_interaction_still_produces_output(self):
        """No interaction-dependent filter: bare Ising is not SU(2)-invariant,
        so chi^{+-} differs from the longitudinal channel legitimately -- and
        that difference is the interesting quantity, not a reason to withhold
        the output. Deciding applicability in the loader would mean encoding the
        rotational-invariance relations between the anisotropic terms, which got
        it wrong in both directions when tried."""
        with tempfile.TemporaryDirectory() as d:
            _fixture(d)
            with open(os.path.join(d, "ising.dat"), "w") as f:
                f.write("Ising\n2\n1\n 1\n"
                        "   0    0    0    1    2   0.500000000000   0.000000000000\n"
                        "   0    0    0    2    1   0.500000000000   0.000000000000\n")
            gi = _solve_with(d, "general", True, {"Ising": "ising.dat"},
                             write_fixture=False)
        self.assertIn("chiq_pm", gi)
        pm = np.asarray(gi["chiq_pm"])
        cs = np.asarray(gi["chiq_s"])
        self.assertFalse(
            np.allclose(pm, cs),
            "bare Ising is not SU(2)-invariant, so the transverse and "
            "longitudinal channels must NOT coincide -- if they do, the "
            "transverse vertex is not being built from the crossed vertex")

    def test_stale_key_cleared_on_the_iterationmax_zero_path(self):
        """The early return must not leave a previous run's chiq_pm behind."""
        with tempfile.TemporaryDirectory() as d:
            gi = _solve(d, "general", True)
            self.assertIn("chiq_pm", gi)
            gi2 = _solve(d, "general", True, {'IterationMax': 0},
                         green_info=gi)
        self.assertNotIn("chiq_pm", gi2)

    def test_no_stale_key_when_the_switch_is_off(self):
        """A reused green_info must not keep a chiq_pm from an earlier solve."""
        with tempfile.TemporaryDirectory() as d:
            gi = _solve(d, "general", True)
            self.assertIn("chiq_pm", gi)
            gi2 = _solve(d, "general", False, green_info=gi)
        self.assertNotIn("chiq_pm", gi2)


if __name__ == "__main__":
    unittest.main()
