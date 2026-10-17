#!/usr/bin/env python3
"""Timing of the FLEX general path under the two second-order kernels.

Measures, on the CPU, (i) one assembly of the effective interaction
``V_eff`` (``_flex_compute_veff_general``) and (ii) one full FLEX iteration
(``solve`` with ``IterationMax = 1``) under ``flex_second_order = "local"``
and ``"takimoto"``, and prints the ratios quoted in the configuration
reference.  Two inputs: the 2-orbital fixture ``tests/rpa/input_2orb``
(geometry and transfer) and a synthetic 3-orbital square lattice written by
this script (nearest-neighbour hopping ``-1 + 0.2 a`` per orbital ``a``,
inter-orbital on-site hopping ``0.3``); both get on-site ``CoulombIntra``
``0.8``, on-site inter-orbital ``CoulombInter`` ``0.4``, ``Hund`` ``0.2`` and
an off-site ``CoulombInter`` ``0.3`` on the ``+x`` bond of orbital 1.

Run from the repository root::

    PYTHONPATH=src:. python3 tests/sc/second_order_cost.py [L] [Nmat] [repeats]

Defaults: ``L = 8``, ``Nmat = 128``, ``repeats = 5`` (the minimum over the
repeats is reported).  Not a test: nothing is asserted.
"""
import os
import shutil
import sys
import tempfile
import time

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from tests.test_second_order_factors import _write_wan  # noqa: E402

_IN2 = "tests/rpa/input_2orb"


def build(norb, second_order, L, nmat):
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as flex_mod
    d = tempfile.mkdtemp()
    if norb == 2:
        for f in ("geom.dat", "transfer.dat"):
            shutil.copy(os.path.join(_IN2, f), d)
    else:
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("1 0 0\n0 1 0\n0 0 1\n{}\n".format(norb) + "".join("0 0 0\n" for _ in range(norb)))
        rows = []
        for a in range(norb):
            for R in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)):
                rows.append(R + (a + 1, a + 1, -1.0 + 0.2 * a, 0.0))
            for b in range(norb):
                if a != b:
                    rows.append((0, 0, 0, a + 1, b + 1, 0.3, 0.0))
        _write_wan(os.path.join(d, "transfer.dat"), "Transfer", norb, rows)
    _write_wan(os.path.join(d, "intra.dat"), "CoulombIntra", norb,
               [(0, 0, 0, a + 1, a + 1, 0.8, 0.0) for a in range(norb)])
    inter = [(0, 0, 0, a + 1, b + 1, 0.4, 0.0) for a in range(norb) for b in range(norb) if a != b]
    inter += [(1, 0, 0, 1, 1, 0.3, 0.0), (-1, 0, 0, 1, 1, 0.3, 0.0)]
    _write_wan(os.path.join(d, "inter.dat"), "CoulombInter", norb, inter)
    _write_wan(os.path.join(d, "hund.dat"), "Hund", norb,
               [(0, 0, 0, a + 1, b + 1, 0.2, 0.0) for a in range(norb) for b in range(norb) if a != b])
    idict = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat",
             "CoulombIntra": "intra.dat", "CoulombInter": "inter.dat", "Hund": "hund.dat"}
    r = read_input_k.QLMSkInput({"path_to_input": d, "interaction": idict})
    par = {"T": 1.0, "filling": 0.5, "CellShape": [L, L, 1], "SubShape": [1, 1, 1], "Nmat": nmat,
           "IterationMax": 1, "Mix": 1.0, "EPS": 1, "flex_second_order": second_order}
    s = flex_mod.FLEX(r.get_param("ham"), {}, {"mode": "FLEX", "param": par,
                                                "enable_spin_orbital": False, "calc_scheme": "general"})
    shutil.rmtree(d, ignore_errors=True)
    return s


def main(L=8, nmat=128, repeats=5):
    beta = 1.0
    for norb in (2, 3):
        veff, full = {}, {}
        for so in ("local", "takimoto"):
            s = build(norb, so, L, nmat)
            s._calc_epsilon_k({})
            nvol = s.lattice.nvol
            G = s._calc_dressed_green(beta, 0.1, np.zeros((1, nmat, nvol, norb, norb), complex))
            chi0q_raw = s._calc_chi0q(G, np.zeros_like(G), beta)[0]
            best = min(_timed(lambda: s._flex_compute_veff_general(chi0q_raw, s.ham_info.ham_inter_q))
                       for _ in range(repeats))
            veff[so] = best
            best = 1e9
            for _ in range(repeats):
                s = build(norb, so, L, nmat)
                out = tempfile.mkdtemp()
                best = min(best, _timed(lambda: s.solve({"path_to_input": out}, out)))
                shutil.rmtree(out, ignore_errors=True)
            full[so] = best
        print("norb {} (nd {}, L {}, Nmat {}): V_eff assembly local {:.4f} s / takimoto {:.4f} s = {:.2f}; "
              "one iteration local {:.4f} s / takimoto {:.4f} s = {:.2f}".format(
                  norb, norb * norb, L, nmat, veff["local"], veff["takimoto"], veff["local"] / veff["takimoto"],
                  full["local"], full["takimoto"], full["local"] / full["takimoto"]))


def _timed(fn):
    t0 = time.perf_counter()
    fn()
    return time.perf_counter() - t0


if __name__ == "__main__":
    args = [int(a) for a in sys.argv[1:4]]
    main(*args)
