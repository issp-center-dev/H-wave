"""End-to-end self-consistency of the general path under the local second
order (#181 follow-up, spec 2026-09-08 section 4).

A converged multi-orbital SCF under ``flex_second_order = "local"`` -- on-site
U, U' and Hund plus an off-site V, Anderson mixing, with
``flex_hartree_fock`` both false and true -- must

* converge, and within the ``"takimoto"`` iteration count + 50 % (the exact
  kernel must not cost self-consistency: it is a different W, not a harder
  fixed-point problem), and
* have a fixed point that IS a fixed point: re-entering the solver from the
  converged self-energy (written as an archive and read back through
  ``[file.input] sigma_init``, the production warm-start path) with
  ``IterationMax = 1`` must reproduce it to 1e-8 relative.

The second leg is what an "it converged" assertion alone cannot give: a
mixing scheme can park on a state that the ITERATION no longer moves while
the MAP still does (Anderson extrapolation, or a residual measured on a
different quantity than the one that is iterated). Reading the archive back
also pins that the kernel provenance of Task 7 survives the round trip.

Both stay in the fast gate: the six solves of :class:`TestSCF` measure 0.7 s
together at Nmat = 64 on a 4x4 cell (T = 1 converges in ~10 Anderson
iterations), well under the 5 s opt-in rule of ``tests/heavy_tests.py``.
:class:`TestSecondOrderSCFSmoke` covers what the SCF case does not: that the
kernel provenance of Task 7 reaches the written archive.
"""
import os
import shutil
import tempfile
import unittest

import numpy as np

from tests.test_second_order_factors import _write_wan

_IN2 = "tests/rpa/input_2orb"

#: on-site U = 0.8, on-site U' = 0.4 and Hund J = 0.2 on the (1,2) pair, plus
#: an off-site V = 0.3 on the +/- x bond: every on-site channel the local
#: kernel weights separately, plus the off-site class it carries as vpair.
_ROWS = {
    "CoulombIntra": [(0, 0, 0, 1, 1, 0.8, 0.0), (0, 0, 0, 2, 2, 0.8, 0.0)],
    "CoulombInter": [(0, 0, 0, 1, 2, 0.4, 0.0), (0, 0, 0, 2, 1, 0.4, 0.0),
                     (1, 0, 0, 1, 1, 0.3, 0.0), (-1, 0, 0, 1, 1, 0.3, 0.0)],
    "Hund": [(0, 0, 0, 1, 2, 0.2, 0.0), (0, 0, 0, 2, 1, 0.2, 0.0)],
}


def _write_inputs(d):
    """The fixture geometry/transfer plus the interaction files; returns the
    ``[file.input.interaction]`` dict."""
    for f in ("geom.dat", "transfer.dat"):
        shutil.copy(os.path.join(_IN2, f), d)
    idict = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
    for name, rows in _ROWS.items():
        fname = name.lower() + ".dat"
        _write_wan(os.path.join(d, fname), name, 2, rows)
        idict[name] = fname
    return idict


def _run(so, hf, itmax=400, nmat=64, sigma_init=None, save_sigma=None):
    """One general-scheme FLEX solve; returns ``(solver, green_info)``.

    ``sigma_init`` is a path to a self-energy archive, consumed through the
    solver's own ``read_init`` exactly as ``[file.input] sigma_init`` does;
    ``save_sigma`` is a path the converged self-energy is written to.
    """
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as flex_mod
    d = tempfile.mkdtemp(prefix="so_scf_in_")
    out = tempfile.mkdtemp(prefix="so_scf_out_")
    try:
        idict = _write_inputs(d)
        r = read_input_k.QLMSkInput({"path_to_input": d, "interaction": idict})
        par = {"T": 1.0, "filling": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1],
               "Nmat": nmat, "IterationMax": itmax, "Mix": 0.3, "EPS": 8,
               "mixing_scheme": "anderson", "anderson_depth": 5,
               "flex_second_order": so, "flex_hartree_fock": hf}
        s = flex_mod.FLEX(r.get_param("ham"), {}, {"mode": "FLEX", "param": par,
                                                   "enable_spin_orbital": False,
                                                   "calc_scheme": "general"})
        gi = r.get_param("green")
        gi.update(s.read_init(
            {"path_to_input": os.path.dirname(sigma_init) if sigma_init else d,
             **({"sigma_init": os.path.basename(sigma_init)} if sigma_init else {})}))
        s.solve(gi, out)
        if save_sigma is not None:
            s.save_results({"path_to_output": os.path.dirname(save_sigma),
                            "sigma": os.path.basename(save_sigma)}, gi)
    finally:
        shutil.rmtree(d, ignore_errors=True)
        shutil.rmtree(out, ignore_errors=True)
    return s, gi


class TestSecondOrderSCFSmoke(unittest.TestCase):
    """Fast-gate subset: the interaction set of the heavy case runs under
    ``"local"`` and the archive carries the kernel provenance."""

    def test_local_runs_this_interaction_set_and_stamps_the_kernel(self):
        with tempfile.TemporaryDirectory() as tmp:
            seed = os.path.join(tmp, "sigma.npz")
            s, _ = _run("local", False, itmax=2, nmat=16, save_sigma=seed)
            self.assertEqual(s.flex_second_order, "local")
            self.assertIsNotNone(s._second_order_factors)
            self.assertIsNotNone(s._second_order_factors.vpair)       # the off-site V is carried
            z = np.load(seed, allow_pickle=True)
            self.assertEqual(str(z["flex_second_order"]), "local")
            self.assertEqual(int(z["flex_second_order_schema"]), 1)
            self.assertGreater(np.abs(z["sigma"]).max(), 1e-6)        # anti-vacuity


class TestSCF(unittest.TestCase):

    def test_converges_and_is_a_fixed_point(self):
        for hf in (False, True):
            with self.subTest(hf=hf), tempfile.TemporaryDirectory() as tmp:
                seed = os.path.join(tmp, "sigma.npz")
                st, _ = _run("takimoto", hf)
                sl, _ = _run("local", hf, save_sigma=seed)
                self.assertTrue(st.scf_converged)
                self.assertTrue(sl.scf_converged)
                self.assertLessEqual(sl.scf_iterations, int(1.5 * st.scf_iterations) + 5)
                print("RECORDED iterations: takimoto {} local {} (hf={})".format(
                    st.scf_iterations, sl.scf_iterations, hf))
                # anti-vacuity: on THIS interaction set the two kernels really
                # are two different fixed points, so "local converged" is a
                # statement about the local kernel and not about a shared one.
                st_scale = np.abs(np.array(st.sigma)).max()
                kernel_gap = np.abs(np.array(sl.sigma) - np.array(st.sigma)).max() / st_scale
                print("RECORDED local-vs-takimoto sigma gap: {:.3e} (hf={})".format(
                    kernel_gap, hf))
                self.assertGreater(kernel_gap, 1.0e-3)
                # The fixed point IS a fixed point: the converged sigma, read
                # back through the production sigma_init seed path, survives one
                # more application of the map.
                sigma0 = np.array(sl.sigma)
                scale = np.abs(sigma0).max()
                self.assertGreater(scale, 1e-6)                   # non-vacuous |sigma0| floor
                s1, _ = _run("local", hf, itmax=1, sigma_init=seed)
                change = np.abs(np.array(s1.sigma) - sigma0).max() / scale
                print("RECORDED fixed-point relative change: {:.3e} (hf={})".format(change, hf))
                self.assertLess(change, 1.0e-8)


if __name__ == "__main__":
    unittest.main()
