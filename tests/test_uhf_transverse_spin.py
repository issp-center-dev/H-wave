"""The transverse spin (Sx, Sy) of UHFk normal mode and of UHFr.

The mean-field state can point along any spin axis, so Sx and Sy can be
nonzero also without spin-orbit coupling. Both solvers report the totals,
as for Sz (in units of hbar):

    Sx = 1/2 Re sum_i (ud_i + du_i),   Sy = 1/2 Re sum_i -i (ud_i - du_i),

with ud_i = <c^dag_{i up} c_{i down}> and du_i = <c^dag_{i down} c_{i up}>
on site, and append them as the last two lines of ``energy.dat``.
"""
import copy
import logging
import os
import shutil
import tempfile
import unittest

import numpy as np
import tomli

from hwave import qlms

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _energy(path):
    """``energy.dat`` as an ordered list of ``(key, float)``."""
    rows = []
    with open(path) as f:
        for line in f:
            key, _, value = line.partition("=")
            rows.append((key.strip(), float(value)))
    return rows


def _stub_calc_phys(green, enable_spin_orbital, nd, norb, ns):
    """``UHFk._calc_phys`` on a one-cell stub, as tests/test_block_matrix.py
    does for the spin-orbital mode."""
    import hwave.solver.uhfk as uhfk_module
    stub = object.__new__(uhfk_module.UHFk)
    stub.shape = (1, 1, 1)
    stub.nvol, stub.nd, stub.norb, stub.ns = 1, nd, norb, ns
    stub.enable_spin_orbital = enable_spin_orbital
    stub.Green = np.asarray(green, dtype=complex)
    stub.Green_prev = stub.Green.copy()
    stub.param_mod = {"Mix": 1.0}
    stub.threshold = 1e-12
    stub.physics = {}
    stub._calc_phys()
    return stub.physics


class _Quiet(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self._level = logging.getLogger().level
        logging.getLogger().setLevel(logging.ERROR)

    def tearDown(self):
        logging.getLogger().setLevel(self._level)
        shutil.rmtree(self.d, ignore_errors=True)

    def _run(self, inp):
        qlms.run(input_dict=copy.deepcopy(inp))


class TestUHFkNormalFormula(unittest.TestCase):

    def test_uhfk_normal_stub(self):
        # G[s, a, t, b] = <c^dag_{a s} c_{b t}>: <c^dag_up c_down> = 0.1 + 0.2i
        rho = np.array([[0.6, 0.1 + 0.2j], [0.1 - 0.2j, 0.4]])
        green = np.zeros((1, 2, 1, 2, 1), dtype=complex)
        green[0][:, 0, :, 0] = rho
        phys = _stub_calc_phys(green, False, nd=2, norb=1, ns=2)
        self.assertAlmostEqual(phys["NCond"], 1.0)
        self.assertAlmostEqual(phys["Sz"], 0.1)
        self.assertAlmostEqual(phys["Sx"], 0.1)
        self.assertAlmostEqual(phys["Sy"], 0.2)

    def test_uhfk_normal_matches_so_formula(self):
        """Normal mode on a spin-block on-site Green == the spin-orbital
        branch on the same Green reordered to the interleaved index."""
        from hwave.solver import spin_orbital
        norb = 3
        rng = np.random.default_rng(20261118)
        a = rng.normal(size=(2 * norb, 2 * norb)) + 1j * rng.normal(size=(2 * norb, 2 * norb))
        mat = (a + a.conj().T) / 2          # spin-block index s * norb + orb
        normal = _stub_calc_phys(mat.reshape(1, 2, norb, 2, norb), False,
                                 nd=2 * norb, norb=norb, ns=2)
        so_green = spin_orbital.spinblock_to_interleaved(mat, norb).reshape(
            1, 2 * norb, 2 * norb)
        so = _stub_calc_phys(so_green, True, nd=2 * norb, norb=2 * norb, ns=1)
        for key in ("Sx", "Sy", "Sz", "NCond"):
            np.testing.assert_allclose(normal[key], so[key], rtol=0, atol=1e-14,
                                       err_msg=key)
        # not vacuous: the random state has a transverse spin
        self.assertGreater(abs(normal["Sx"]), 1e-3)
        self.assertGreater(abs(normal["Sy"]), 1e-3)


class TestUHFkRuns(_Quiet):

    def _coulombintra(self, out):
        case = os.path.join(_REPO, "tests", "uhfk", "CoulombIntra")
        with open(os.path.join(case, "input.toml"), "rb") as f:
            inp = tomli.load(f)
        inp["log"]["print_level"] = 0
        inp["file"]["input"]["path_to_input"] = case
        inp["file"]["input"]["interaction"]["path_to_input"] = case
        inp["file"]["output"]["path_to_output"] = out
        return inp

    def test_uhfk_ferromagnet_generic_axis(self):
        """A saturated ferromagnet along a generic axis u: the total spin is
        |S| = Ncond / 2 = 2.5 and points along u, which checks Sx, Sy and
        their signs together.

        The test relies on the SCF keeping the start's global spin axis: the
        Hubbard U is spin-rotation invariant and the start is a fully
        polarised state along u, so the totals (Sx, Sy, Sz) are fixed from
        step 0 and the SCF only redistributes the charge.
        """
        n = 5 / 16
        u = np.array([0.6, 0.2537, 0.7734])
        u = u / np.linalg.norm(u)
        sigma = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]),
                 np.array([[1, 0], [0, -1]])]
        rho = (n / 2) * (np.eye(2) + sum(c * s for c, s in zip(u, sigma)))
        green = np.zeros((16, 2, 1, 2, 1), dtype=complex)
        green[0, :, 0, :, 0] = rho.T       # G_st = <c^dag_s c_t> = rho_ts
        start = os.path.join(self.d, "start.npz")
        np.savez(start, green=green, momentum_convention="e_plus_ikR")

        out = os.path.join(self.d, "out")
        inp = self._coulombintra(out)
        inp["mode"]["param"]["Ncond"] = 5
        inp["mode"]["param"]["SubShape"] = [1, 1, 1]
        inp["file"]["input"]["initial"] = start
        self._run(inp)
        e = dict(_energy(os.path.join(out, "energy.dat")))
        s = np.array([e["Sx"], e["Sy"], e["Sz"]])
        self.assertAlmostEqual(np.linalg.norm(s), 2.5, delta=1e-9)
        np.testing.assert_allclose(s / 2.5, u, rtol=0, atol=1e-8)
        with open(os.path.join(out, "scf_status.dat")) as f:
            status = dict(line.rstrip("\n").split(" = ", 1) for line in f)
        self.assertEqual(status["converged"], "true")
        self.assertLessEqual(int(status["iterations"]), 60)

    def test_uhfk_tilted_neel_totals(self):
        """The golden U = 8 state is a collinear Neel state with a tilted
        axis: the sublattice moments have transverse parts, but the totals
        cancel."""
        out = os.path.join(self.d, "out")
        self._run(self._coulombintra(out))
        rows = _energy(os.path.join(out, "energy.dat"))
        self.assertEqual([k for k, _ in rows[-2:]], ["Sx", "Sy"])
        e = dict(rows)
        self.assertLess(abs(e["Sx"]), 1e-10)
        self.assertLess(abs(e["Sy"]), 1e-10)


class TestUHFr(_Quiet):

    FIELDS = {
        # rows (i s j t re im) of H = -sum t c^dag_{is} c_{jt}, giving H = -sigma_alpha
        "Sx": [(0, 0, 0, 1, 1, 0), (0, 1, 0, 0, 1, 0)],
        "Sy": [(0, 0, 0, 1, 0, -1), (0, 1, 0, 0, 0, 1)],
        "Sz": [(0, 0, 0, 0, 1, 0), (0, 1, 0, 1, -1, 0)],
    }

    def _field_case(self, name, rows):
        d = os.path.join(self.d, name)
        os.makedirs(d)
        body = "".join("{} {} {} {} {} {}\n".format(*r) for r in rows)
        with open(os.path.join(d, "trans.def"), "w") as f:
            f.write("========================\nNTransfer {}\n"
                    "========================\n========i_j_s_tijs======\n"
                    "========================\n".format(len(rows)) + body)
        inp = {
            "log": {"print_level": 0},
            "mode": {"mode": "UHFr", "param": {
                "Nsite": 1, "Ncond": 1, "T": 0.0, "EPS": 14, "IterationMax": 100}},
            "file": {"input": {"path_to_input": d, "interaction": {"Trans": "trans.def"}},
                     "output": {"path_to_output": os.path.join(d, "out"),
                                "energy": "energy.dat"}},
        }
        self._run(inp)
        return _energy(os.path.join(d, "out", "energy.dat"))

    def test_uhfr_field_along_each_axis(self):
        for comp, rows in self.FIELDS.items():
            with self.subTest(component=comp):
                energy = self._field_case(comp, rows)
                self.assertEqual([k for k, _ in energy[-2:]], ["Sx", "Sy"])
                e = dict(energy)
                for other in ("Sx", "Sy", "Sz"):
                    want = 0.5 if other == comp else 0.0
                    self.assertAlmostEqual(e[other], want, places=8, msg=other)

    def test_uhfr_fixed_2sz_transverse_zero(self):
        """With 2Sz fixed the two spin blocks are separate (no spin mixing),
        so Sx and Sy are exactly zero."""
        case = os.path.join(_REPO, "docs", "tutorial", "Hubbard", "UHFr")
        with open(os.path.join(case, "input.toml"), "rb") as f:
            inp = tomli.load(f)
        self.assertEqual(inp["mode"]["param"]["2Sz"], 0)
        inp["log"]["print_level"] = 0
        inp["file"]["input"]["path_to_input"] = case
        out = os.path.join(self.d, "out")
        inp["file"]["output"]["path_to_output"] = out
        self._run(inp)
        rows = _energy(os.path.join(out, "energy.dat"))
        self.assertEqual([k for k, _ in rows[-2:]], ["Sx", "Sy"])
        e = dict(rows)
        self.assertEqual(e["Sx"], 0.0)
        self.assertEqual(e["Sy"], 0.0)


if __name__ == "__main__":
    unittest.main()
