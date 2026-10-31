"""The bond topology of the spin-orbital vertex (``so_bond.so_bond_topology``)
and its mesh check (``so_bond.slot_geometry``)."""
import unittest

import numpy as np

from hwave.solver import so_bond
from hwave.solver.so_bond import so_bond_topology
from hwave.solver.so_records import compile_interaction
from tests.so_bond_oracles import COMPILER_TABLES, CELL, closed_row

_PM12 = {(0, 0, 0), (1, 0, 0), (-1, 0, 0), (2, 0, 0), (-2, 0, 0)}


def _zero_outer_shell():
    """CoulombInter with nonzero rows on R = +-1 and declared-zero rows on
    R = +-2 (norb_phys = 2)."""
    tbl = {}
    for row in (closed_row(1, 0, 1, 0.37), closed_row(-1, 0, 1, 0.22),
                closed_row(2, 0, 1, 0.0), closed_row(-2, 1, 1, 0.0)):
        tbl.update(row)
    return {"CoulombInter": tbl}


class TestSOTopology(unittest.TestCase):
    """A1-C (topology part): the channels come from the declared topology of
    the compiled records; a cutoff may drop declared-zero shells only."""

    def test_channels_from_declared_topology(self):
        topo = so_bond_topology(compile_interaction(COMPILER_TABLES, 2), CELL)
        self.assertEqual(tuple(topo.delta_r[0]), (0, 0, 0))
        self.assertEqual({tuple(r) for r in topo.delta_r}, _PM12)
        self.assertEqual([tuple(r) for r in topo.delta_r],
                         [(0, 0, 0), (-1, 0, 0), (1, 0, 0), (-2, 0, 0), (2, 0, 0)])
        rev = np.asarray(topo.reverse)
        self.assertTrue(np.all(rev[rev] == np.arange(len(rev))))
        self.assertTrue(np.all(topo.delta_r[rev] == -topo.delta_r))
        self.assertEqual(topo.coeffs, {})

    def test_declared_zero_shell_is_a_channel(self):
        c = compile_interaction(_zero_outer_shell(), 2)
        self.assertEqual({tuple(r) for r in so_bond_topology(c, CELL).delta_r}, _PM12)

    def test_cutoff_rules(self):
        c = compile_interaction(COMPILER_TABLES, 2)
        with self.assertRaises(ValueError) as cm:
            so_bond_topology(c, CELL, max_shells=1)
        self.assertEqual(
            str(cm.exception),
            "so_bond_topology: max_shells=1 would drop bond shell(s) used by nonzero "
            "interaction records: shell 2 with displacements [(-2, 0, 0), (2, 0, 0)]. "
            "Use max_shells >= 2, or remove those declarations.")
        with self.assertRaises(ValueError) as cm:
            so_bond_topology(c, CELL, max_shells=0)
        self.assertEqual(
            str(cm.exception),
            "so_bond_topology: max_shells=0 would drop bond shell(s) used by nonzero "
            "interaction records: shell 1 with displacements [(-1, 0, 0), (1, 0, 0)]; "
            "shell 2 with displacements [(-2, 0, 0), (2, 0, 0)]. "
            "Use max_shells >= 2, or remove those declarations.")
        self.assertEqual(len(so_bond_topology(c, CELL, max_shells=2).delta_r), 5)
        cut = so_bond_topology(compile_interaction(_zero_outer_shell(), 2), CELL, max_shells=1)
        self.assertEqual(len(cut.delta_r), 3)
        self.assertEqual([tuple(r) for r in cut.delta_r], [(0, 0, 0), (-1, 0, 0), (1, 0, 0)])
        self.assertEqual(cut.reverse.tolist(), [0, 2, 1])

    def test_onsite_only(self):
        c = compile_interaction({"CoulombIntra": {((0, 0, 0), (0, 0)): 1.0}}, 1)
        topo = so_bond_topology(c, CELL, max_shells=0)
        self.assertEqual(topo.delta_r.tolist(), [[0, 0, 0]])
        self.assertEqual(topo.reverse.tolist(), [0])

    def test_refuses_an_uncompiled_interaction(self):
        with self.assertRaises(ValueError) as cm:
            so_bond_topology(COMPILER_TABLES, CELL)
        self.assertEqual(str(cm.exception),
                         "so_bond_topology: compiled must be a CompiledInteraction "
                         "(so_records.compile_interaction), got dict")

    def test_wrapped_displacement_refused(self):
        """Channels distinct as lattice vectors but equal modulo the mesh
        (R = 4 and R = -1 on L = 5) are refused at the mesh check."""
        c = compile_interaction({"CoulombInter": {**closed_row(4, 0, 1, 0.5),
                                                  **closed_row(-1, 0, 1, 0.2)}}, 2)
        topo = so_bond_topology(c, CELL)
        with self.assertRaises(ValueError) as cm:
            so_bond.slot_geometry(topo, (5, 1, 1))
        self.assertEqual(
            str(cm.exception),
            "validate_topology_against_mesh: channels 2 (delta_r=(1, 0, 0)) and 3 "
            "(delta_r=(-4, 0, 0)) collide on the mesh [5, 1, 1] -- both wrap to the same "
            "displacement (1, 0, 0) modulo the mesh (an alias such as R vs R+k*L on some "
            "axis, or the even-mesh self-reversal +L/2 == -L/2); a topology whose channels "
            "are not INJECTIVE modulo the mesh cannot be used at this mesh resolution.")
        geom = so_bond.slot_geometry(topo, (9, 1, 1))
        self.assertEqual(geom.B, 5)


if __name__ == "__main__":
    unittest.main()
