"""``resolve_bond_topology(..., displacements=)``: the coefficient-free channel
set of an explicit displacement list, and the unchanged default path."""
import unittest
from collections import namedtuple

import numpy as np

from hwave.solver import bond_channels as bc
from tests.so_bond_oracles import CELL, T1
from tests.so_bond_oracles import merge_tables as _merge

_Fixture = namedtuple("_Fixture", "args kw")


def _decl(t, R, a, b, v):
    R = tuple(R)
    return {t: {(R, (a, b)): v, (tuple(-x for x in R), (b, a)): np.conj(v)}}


def _two_shell_square(v1=0.3, v2=0.1):
    ci = {}
    for R in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)):
        ci[(R, (0, 0))] = v1
    for R in ((1, 1, 0), (1, -1, 0), (-1, 1, 0), (-1, -1, 0)):
        ci[(R, (0, 0))] = v2
    return ci


_TRIANGULAR = np.array([[1.0, 0.0, 0.0], [0.5, np.sqrt(3) / 2, 0.0], [0.0, 0.0, 1.0]])

#: resolver calls of the non-spin-orbital paths (longitudinal, transverse,
#: PairHop, truncated, two-orbital, non-orthogonal cell)
NON_SO_FIXTURES = (
    _Fixture(({"Hund": T1}, CELL, 2), dict(active_types=("Hund",))),
    _Fixture((_decl("Hund", (1, 0, 0), 0, 0, 0.3), np.eye(3), 1),
             dict(active_types=bc._LONGITUDINAL_ACTIVE_TYPES)),
    _Fixture((_decl("Exchange", (1, 0, 0), 0, 1, 0.2), np.eye(3), 2),
             dict(active_types=bc._TRANSVERSE_ACTIVE_TYPES)),
    _Fixture(({"PairHop": {((1, 0, 0), (0, 0)): 0.3 + 0.1j, ((-1, 0, 0), (0, 0)): 0.3 - 0.1j}},
              np.eye(3), 1), dict(active_types=("PairHop",))),
    _Fixture(({"CoulombInter": _two_shell_square()}, np.eye(3), 1),
             dict(active_types=("CoulombInter",))),
    _Fixture(({"CoulombInter": _two_shell_square(v2=0.0)}, np.eye(3), 1),
             dict(max_shells=1, active_types=("CoulombInter",))),
    _Fixture((_merge(_decl("Hund", (2, 0, 0), 0, 0, 0.5), _decl("CoulombInter", (1, 0, 0), 0, 0, 0.1)),
              np.eye(3), 1), dict(max_shells=2, active_types=bc._LONGITUDINAL_ACTIVE_TYPES)),
    _Fixture(({"CoulombInter": {((1, 0, 0), (0, 1)): 0.5 + 0.1j,
                                ((-1, 0, 0), (1, 0)): 0.5 - 0.1j + 5e-9}}, np.eye(3), 2),
             dict(active_types=("CoulombInter",))),
    _Fixture((_merge(_decl("CoulombInter", (1, 0, 0), 0, 1, 0.4), _decl("Ising", (0, 1, 0), 1, 1, 0.2),
                     _decl("Exchange", (1, -1, 0), 1, 0, 0.15)), _TRIANGULAR, 2),
             dict(active_types=("CoulombInter", "Ising", "Exchange"))),
)

_DISP = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0], [2, 0, 0], [-2, 0, 0]], dtype=np.int64)

_NONEMPTY_MESSAGE = (
    "resolve_bond_topology: displacements= builds a coefficient-free channel set; it requires "
    "empty interactions and active_types=(), got interaction type(s) {} and active_types {}")


def snapshot(fixtures=NON_SO_FIXTURES):
    """``{name: bytes-exact array}`` of every fixture's resolved topology."""
    out = {}
    for i, fx in enumerate(fixtures):
        topo = bc.resolve_bond_topology(*fx.args, **fx.kw)
        out["f{}_delta_r".format(i)] = np.array(topo.delta_r)
        out["f{}_reverse".format(i)] = np.array(topo.reverse)
        for t in topo.coeffs:
            out["f{}_coeffs_{}".format(i, t)] = np.array(topo.coeffs[t])
    return out


class TestNonSOUnchanged(unittest.TestCase):
    """The default ``displacements=None`` is the existing resolver."""

    def test_default_path_is_today(self):
        """Passing ``displacements=None`` explicitly gives the same topology
        as omitting it. This guards the keyword plumbing only: both calls run
        the same resolver code, so it does not pin the resolver's output."""
        for i, fx in enumerate(NON_SO_FIXTURES):
            with self.subTest(fixture=i):
                a = bc.resolve_bond_topology(*fx.args, **fx.kw)
                b = bc.resolve_bond_topology(*fx.args, displacements=None, **fx.kw)
                self.assertEqual(np.asarray(a.delta_r).tobytes(), np.asarray(b.delta_r).tobytes())
                self.assertEqual(np.asarray(a.reverse).tobytes(), np.asarray(b.reverse).tobytes())
                self.assertEqual(tuple(a.coeffs), tuple(b.coeffs))
                self.assertGreater(len(a.coeffs), 0)
                for t in a.coeffs:
                    self.assertEqual(a.coeffs[t].tobytes(), b.coeffs[t].tobytes())

    def test_displacements_requires_empty_inputs(self):
        with self.assertRaises(ValueError) as cm:
            bc.resolve_bond_topology({"Hund": T1}, CELL, 2, active_types=("Hund",),
                                     displacements=_DISP)
        self.assertEqual(str(cm.exception), _NONEMPTY_MESSAGE.format(["Hund"], ("Hund",)))
        with self.assertRaises(ValueError) as cm:
            bc.resolve_bond_topology({}, CELL, 2, active_types=("Hund",), displacements=_DISP)
        self.assertEqual(str(cm.exception), _NONEMPTY_MESSAGE.format([], ("Hund",)))
        with self.assertRaises(ValueError) as cm:
            bc.resolve_bond_topology({"Hund": T1}, CELL, 2, active_types=(), displacements=_DISP)
        self.assertEqual(str(cm.exception), _NONEMPTY_MESSAGE.format(["Hund"], ()))


class TestDisplacementSet(unittest.TestCase):
    """The channels of a displacement list: channel 0 on-site, shells by
    length then lexicographically, the reversal, no coefficients."""

    def test_channels_and_order(self):
        topo = bc.resolve_bond_topology({}, CELL, 2, active_types=(), displacements=_DISP[::-1])
        self.assertIsInstance(topo, bc.BondTopology)
        self.assertEqual(topo.coeffs, {})
        self.assertEqual([tuple(r) for r in topo.delta_r],
                         [(0, 0, 0), (-1, 0, 0), (1, 0, 0), (-2, 0, 0), (2, 0, 0)])
        self.assertEqual(topo.reverse.tolist(), [0, 2, 1, 4, 3])
        self.assertEqual(topo.delta_r.dtype, np.int64)

    def test_same_channels_as_the_declared_shells(self):
        """The channels equal those of the same displacements declared as
        zero-valued CoulombInter rows, for every cutoff."""
        disp = [(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (1, 1, 0), (-1, -1, 0),
                (1, -1, 0), (-1, 1, 0), (2, 0, 0), (-2, 0, 0)]
        rows = {(d, (0, 0)): 0.0 for d in disp if any(d)}
        for max_shells in (None, 0, 1, 2, 3, 7):
            with self.subTest(max_shells=max_shells):
                want = bc.resolve_bond_topology({"CoulombInter": rows}, _TRIANGULAR, 1,
                                                max_shells=max_shells, active_types=("CoulombInter",))
                got = bc.resolve_bond_topology({}, _TRIANGULAR, 1, max_shells=max_shells,
                                               active_types=(), displacements=disp)
                self.assertEqual(got.delta_r.tolist(), want.delta_r.tolist())
                self.assertEqual(got.reverse.tolist(), want.reverse.tolist())

    def test_malformed_displacements_refused(self):
        cases = (
            ([(1, 0)], "resolve_bond_topology: displacements must be an (n, 3) array of integer "
                       "displacements, got shape (1, 2)"),
            ([(0.5, 0, 0), (-0.5, 0, 0)],
             "resolve_bond_topology: displacements must be integral, got [[0.5, 0.0, 0.0], "
             "[-0.5, 0.0, 0.0]]"),
            ([(0, 0, 0), (1, 0, 0)],
             "resolve_bond_topology: displacements must be closed under negation; (1, 0, 0) "
             "has no partner (-1, 0, 0)"),
        )
        for disp, msg in cases:
            with self.subTest(disp=disp):
                with self.assertRaises(ValueError) as cm:
                    bc.resolve_bond_topology({}, CELL, 1, active_types=(), displacements=disp)
                self.assertEqual(str(cm.exception), msg)


if __name__ == "__main__":
    unittest.main()
