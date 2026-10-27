"""hwave_sc's static ``chi0q_mode = "calc"`` computes chi0 at the zero bosonic
frequency only (issue #272).

* ``TestStaticBubbleKernel``: ``bubble.dense_bubble_static`` equals the
  ``nmat // 2`` slice of ``bubble.dense_bubble`` on random Green functions
  (reduced / general scheme, one and two blocks, tail on / off), also when
  the imaginary-time sum runs over several chunks; the row weights equal the
  ``tau_to_boson`` row.
* ``TestStaticSliceEquality``: ``sc._calc_chi0q_internal(static_only=True)``
  equals the static slice of the full-frequency RPA chi0q for the
  non-spin-orbital reduced and general schemes and the spin-orbital (spinful
  general) bubble, with ``coeff_tail`` 0 and 1 and several even ``Nmat``.
* ``TestEigenvaluesUnchanged``: ``calc_eliashberg`` with the static-only
  chi0q gives the eigenvalues of the full-frequency path (non-spin-orbital
  reduced, general, spin-orbital).
* ``TestPeakMemory``: the peak traced memory of the static-only chi0q no
  longer grows with ``Nmat`` like the full-frequency bubble does (measured in
  a subprocess with ``tracemalloc``).
"""
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest
from unittest import mock

import numpy as np

import hwave.sc as sc
from hwave.solver import bubble
from hwave.solver import matsubara as ms
from tests.test_sc_spin_orbital import _input, _quiet, _run, _write_model

#: relative bound of the static-slice equality (round-off of a reordered sum)
_SLICE_RTOL = 1e-13
#: relative bound of the eigenvalue equality
_EIG_RTOL = 1e-12


def _rel(a, b):
    return np.abs(a - b).max() / np.abs(b).max()


class TestStaticBubbleKernel(unittest.TestCase):

    def _green(self, rng, nblock, nmat, nvol, nd, tail):
        g = (rng.standard_normal((nblock, nmat, nvol, nd, nd))
             + 1j * rng.standard_normal((nblock, nmat, nvol, nd, nd)))
        if not tail:
            return g, None
        # a frequency-constant tail, as build_green produces
        t = (rng.standard_normal((nblock, 1, nvol, nd, nd))
             + 1j * rng.standard_normal((nblock, 1, nvol, nd, nd)))
        return g, np.repeat(t, nmat, axis=1)

    def _check(self, nblock, scheme, tail, nmat=16, chunk_bytes=None):
        rng = np.random.default_rng(nblock * 100 + nmat + 7 * tail)
        shape, nd = (2, 3, 1), 2
        g, t = self._green(rng, nblock, nmat, 6, nd, tail)
        full = bubble.dense_bubble(g, t, 3.0, spatial_shape=shape,
                                   scheme=scheme)
        if chunk_bytes is None:
            st = bubble.dense_bubble_static(g, t, 3.0, spatial_shape=shape,
                                            scheme=scheme)
        else:
            with mock.patch.object(bubble, "STATIC_BUBBLE_CHUNK_BYTES",
                                   chunk_bytes):
                st = bubble.dense_bubble_static(g, t, 3.0,
                                                spatial_shape=shape,
                                                scheme=scheme)
        ref = full[:, nmat // 2:nmat // 2 + 1]
        self.assertEqual(st.shape, ref.shape)
        self.assertLess(_rel(st, ref), _SLICE_RTOL)
        return st, ref

    def test_slice_equality(self):
        for nblock in (1, 2):
            for scheme in ("reduced", "general"):
                for tail in (False, True):
                    with self.subTest(nblock=nblock, scheme=scheme, tail=tail):
                        self._check(nblock, scheme, tail)

    def test_several_chunks(self):
        """A chunk budget of three tau slices (the last chunk partial): the
        chunked sum still equals the slice, the tail endpoint included."""
        for scheme, nds in (("reduced", 4), ("general", 16)):
            slice_bytes = 16 * 2 * 6 * nds
            with self.subTest(scheme=scheme):
                self._check(2, scheme, True, nmat=16,
                            chunk_bytes=3 * slice_bytes)

    def test_tail_is_exercised(self):
        """The tail fixture changes the static slice at O(1) (so the
        equality above covers the endpoint correction)."""
        rng = np.random.default_rng(1)
        g, t = self._green(rng, 1, 16, 6, 2, True)
        on = bubble.dense_bubble_static(g, t, 3.0, spatial_shape=(2, 3, 1),
                                        scheme="general")
        off = bubble.dense_bubble_static(g, np.zeros_like(t), 3.0,
                                         spatial_shape=(2, 3, 1),
                                         scheme="general")
        self.assertGreater(_rel(on, off), 1e-3)

    def test_row_weights(self):
        rng = np.random.default_rng(3)
        for n in (2, 8, 64, 7):
            a = rng.standard_normal(n) + 1j * rng.standard_normal(n)
            full = ms.tau_to_boson(a, axis=0)
            for idx in {n // 2, 0, n - 1}:
                w = ms.tau_to_boson_row_weights(n, idx)
                self.assertLess(abs(w @ a - full[idx]),
                                1e-14 * np.abs(a).sum())


class _Tmp(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.d, ignore_errors=True)


class TestStaticSliceEquality(_Tmp):
    """``static_only`` equals the full-frequency slice ``nmat // 2``."""

    CASES = (
        # (spin-orbital, interaction types, chi0q_tensor)
        (False, ("CoulombIntra",), "reduced"),
        (False, ("CoulombIntra",), "general"),
        (False, None, "general"),
        (True, None, "general"),
    )

    def test_static_slice(self):
        for so, types, tensor in self.CASES:
            kw = {} if types is None else {"types": types}
            files = _write_model(self.d, so=so, soc=0.3 if so else 0.0, **kw)
            for ct in (0.0, 1.0):
                by_tail = []
                for nmat in (2, 16, 64):
                    with self.subTest(so=so, types=types, tensor=tensor,
                                      coeff_tail=ct, nmat=nmat):
                        inp = _input(self.d, files, so,
                                     {"chi0q_mode": "calc"}, nmat=nmat,
                                     coeff_tail=ct)
                        with _quiet():
                            full = sc._calc_chi0q_internal(
                                inp, chi0q_tensor=tensor, precomputed_mu=0.1)
                            st = sc._calc_chi0q_internal(
                                inp, chi0q_tensor=tensor, precomputed_mu=0.1,
                                static_only=True)
                        ref = full[nmat // 2:nmat // 2 + 1]
                        self.assertEqual(st.shape, (1,) + full.shape[1:])
                        self.assertLess(_rel(st, ref), _SLICE_RTOL)
                        by_tail.append(st)
                if ct == 0.0:
                    no_tail = by_tail
                else:
                    # the tail matters at moderate Nmat, so the comparison
                    # above really covers the tail treatment
                    self.assertGreater(_rel(by_tail[1], no_tail[1]), 1e-3)


class TestEigenvaluesUnchanged(_Tmp):
    """``calc_eliashberg`` on the static-only chi0q gives the eigenvalues of
    the full-frequency chi0q (the pre-#272 path, reproduced by computing the
    full grid and handing its static slice on)."""

    def _compare(self, so, types, eli, **param):
        kw = {} if types is None else {"types": types}
        files = _write_model(self.d, so=so, soc=0.3 if so else 0.0, **kw)
        eli = dict({"chi0q_mode": "calc", "solver_mode": "eigenvalue",
                    "num_eigenvalues": 6}, **eli)
        real = sc._calc_chi0q_internal
        calls = []

        def full_path(input_dict, chi0q_tensor="auto", precomputed_mu=None,
                      static_only=False):
            calls.append(static_only)
            chi = real(input_dict, chi0q_tensor=chi0q_tensor,
                       precomputed_mu=precomputed_mu)
            nmat = chi.shape[0]
            return chi[nmat // 2:nmat // 2 + 1]

        new = _run(_input(self.d, files, so, eli, out="new", **param))
        with mock.patch.object(sc, "_calc_chi0q_internal",
                               side_effect=full_path):
            old = _run(_input(self.d, files, so, eli, out="old", **param))
        self.assertEqual(calls, [True])
        new_ev, old_ev = np.asarray(new["evs"]), np.asarray(old["evs"])
        self.assertEqual(new_ev.shape, old_ev.shape)
        np.testing.assert_allclose(new_ev, old_ev, rtol=_EIG_RTOL, atol=0)

    def test_reduced(self):
        for ct in (0.0, 1.0):
            with self.subTest(coeff_tail=ct):
                self._compare(False, ("CoulombIntra",),
                              {"chi0q_tensor": "reduced"}, coeff_tail=ct)

    def test_general(self):
        for ct in (0.0, 1.0):
            with self.subTest(coeff_tail=ct):
                self._compare(False, None, {"chi0q_tensor": "general"},
                              coeff_tail=ct)

    def test_spin_orbital(self):
        for ct in (0.0, 1.0):
            with self.subTest(coeff_tail=ct):
                self._compare(True, None, {"num_eigenvalues": 8},
                              coeff_tail=ct)


_MEMORY_SCRIPT = textwrap.dedent("""
    import json, sys, tempfile, tracemalloc
    import hwave.sc as sc
    from tests.test_sc_spin_orbital import _input, _quiet, _write_model
    nmat, static_only = int(sys.argv[1]), sys.argv[2] == "1"
    d = tempfile.mkdtemp()
    files = _write_model(d, so=True, soc=0.3)
    inp = _input(d, files, True, {"chi0q_mode": "calc"}, nmat=nmat,
                 coeff_tail=1.0, CellShape=[8, 8, 1])
    with _quiet():
        tracemalloc.start()
        chi = sc._calc_chi0q_internal(inp, chi0q_tensor="general",
                                      precomputed_mu=0.1,
                                      static_only=static_only)
        peak = tracemalloc.get_traced_memory()[1]
        tracemalloc.stop()
    print(json.dumps({"peak": peak, "shape": list(chi.shape)}))
""")


class TestPeakMemory(unittest.TestCase):
    """Spin-orbital, nd = 4, 8 x 8 k points: the full bubble is
    ``16 * Nmat * 64 * 4**4`` bytes (1 MiB per frequency)."""

    NMATS = (32, 128)

    def _peak(self, nmat, static_only):
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        env = dict(os.environ)
        env["PYTHONPATH"] = os.pathsep.join(
            [os.path.join(root, "src"), root, env.get("PYTHONPATH", "")])
        out = subprocess.run(
            [sys.executable, "-B", "-c", _MEMORY_SCRIPT, str(nmat),
             "1" if static_only else "0"],
            cwd=root, env=env, capture_output=True, text=True, check=True)
        import json
        # the profiler may print its statistics after the record
        line = [ln for ln in out.stdout.splitlines() if ln.startswith("{")][-1]
        rec = json.loads(line)
        return rec["peak"], rec["shape"]

    def test_peak_does_not_scale_with_nmat(self):
        new = {n: self._peak(n, True) for n in self.NMATS}
        old = {n: self._peak(n, False) for n in self.NMATS}
        lo, hi = self.NMATS
        self.assertEqual(new[hi][1][0], 1)
        self.assertEqual(old[hi][1][0], hi)
        d_new = new[hi][0] - new[lo][0]
        d_old = old[hi][0] - old[lo][0]
        # the full bubble alone grows by (hi - lo) MiB
        self.assertGreater(d_old, (hi - lo) * (1 << 20))
        # what still grows with Nmat is the Green function (nd**2, not nd**4)
        self.assertLess(d_new, 0.2 * d_old)
        self.assertLess(new[hi][0], 0.3 * old[hi][0])


if __name__ == "__main__":
    unittest.main()
