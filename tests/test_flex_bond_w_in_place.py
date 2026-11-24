"""``W`` written in place over ``chibar`` in the bond block store (issue #222).

Without in-process pairing the FLEX map stores the effective interaction
``W`` in the slot of the bond bubble ``chibar`` (``BondBlockStore`` alias
``"W" -> "chibar"``): each dressing batch reads its ``chibar`` batch before
it writes the ``W`` batch over it, and the next map rebuilds ``chibar``
before reading it. The reorder does not change any arithmetic, so the
converged state must be BIT-identical to the two-slot layout; the memory
table's ``chibar_W`` row drops by one ``(nmat, nvol, ND, ND)`` buffer. With
in-process pairing the last map's ``chibar`` is read after the solve, so
that run keeps the two slots.
"""
import tempfile
import unittest
from unittest import mock

import numpy as np

from hwave.solver import flex_bond
from tests.test_flex_bond_gate import _flex


class TestStoreAlias(unittest.TestCase):

    def _store(self):
        return flex_bond.BondBlockStore(4, 2, 6, 3, ("chibar",), {"W": "chibar"})

    def test_alias_shares_the_array(self):
        with self._store() as st:
            st.put_freq_batch("chibar", 0, 4, np.full((4, 2, 6, 6), 1.0 + 2.0j))
            st.put_freq_batch("W", 1, 3, np.full((2, 2, 6, 6), 5.0))
            got = st.get_freq_batch("chibar", 0, 4)
            np.testing.assert_array_equal(got[[0, 3]], 1.0 + 2.0j)
            np.testing.assert_array_equal(got[1:3], 5.0)
            np.testing.assert_array_equal(st.get_pair("W", 1, 0),
                                          st.get_pair("chibar", 1, 0))

    def test_alias_cannot_be_detached_or_released(self):
        with self._store() as st:
            with self.assertRaisesRegex(ValueError, "alias"):
                st.detach("W")
            with self.assertRaisesRegex(ValueError, "alias"):
                st.release_slot("W")
            st.get_freq_batch("chibar", 0, 1)            # still allocated

    def test_slot_names(self):
        with flex_bond.BondBlockStore(2, 1, 2, 1, ("chibar", "chi"), {"W": "chibar"}) as st:
            self.assertEqual(st.slot_names, ("chibar", "chi", "W"))
            st.release_slot("chi")
            self.assertEqual(st.slot_names, ("chibar", "W"))
            with self.assertRaises(AttributeError):
                st.slot_names = ()

    def test_slot_names_drop_aliases_of_freed_targets(self):
        with flex_bond.BondBlockStore(2, 1, 2, 1, ("chibar", "chi"), {"W": "chibar"}) as st:
            st.release_slot("chibar")
            self.assertEqual(st.slot_names, ("chi",))
        with flex_bond.BondBlockStore(2, 1, 2, 1, ("chibar", "chi"), {"W": "chibar"}) as st:
            st.detach("chibar")
            self.assertEqual(st.slot_names, ("chi",))
        st = flex_bond.BondBlockStore(2, 1, 2, 1, ("chibar",), {"W": "chibar"})
        st.release()
        self.assertEqual(st.slot_names, ())

    def test_bad_alias_refused(self):
        with self.assertRaisesRegex(ValueError, "alias"):
            flex_bond.BondBlockStore(2, 1, 2, 1, ("chibar",), {"W": "missing"})
        with self.assertRaisesRegex(ValueError, "alias"):
            flex_bond.BondBlockStore(2, 1, 2, 1, ("chibar", "W"), {"W": "chibar"})


class TestMemoryRow(unittest.TestCase):

    def test_chibar_w_row_drops_by_one_buffer(self):
        kw = dict(nmat=64, nvol=16, norb=2, B=5, depth=0, output_full=False,
                  split_seed=False, n_types=1, freq_batch=None, cap_gb=1e6, mixing="linear")
        two = flex_bond.estimate_bond_memory(**kw)
        one = flex_bond.estimate_bond_memory(w_in_place=True, **kw)
        U = 64 * 16 * (5 * 4) ** 2 * 16
        self.assertEqual(two["persistent_rows"]["chibar_W"], 2 * U)
        self.assertEqual(one["persistent_rows"]["chibar_W"], U)
        self.assertEqual(two["persistent"] - one["persistent"], U)


class TestSolverLayout(unittest.TestCase):

    def test_layout_follows_the_pairing_switch(self):
        s, _ = _flex()
        self.assertEqual(s._bond_store_layout(), (("chibar",), {"W": "chibar"}))
        s, _ = _flex({"longitudinal_bond_output_full": True})
        self.assertEqual(s._bond_store_layout(),
                         (("chibar", "chi_s_w", "chi_c_w"), {"W": "chibar"}))
        s, _ = _flex({"longitudinal_bond_pairing": "singlet"})
        self.assertFalse(s._bond_w_in_place())
        self.assertEqual(s._bond_store_layout(), (("chibar", "W"), {}))


def _solve(extra, in_place):
    """Run the gated FLEX solve with the given layout; return the state and
    the aliases the store was built with."""
    built = []
    real = flex_bond.BondBlockStore

    def spy(*a, **k):
        st = real(*a, **k)
        built.append(dict(st._aliases))
        return st

    s, r = _flex(dict({"IterationMax": 4}, **extra))
    gi = r.get_param("green")
    with tempfile.TemporaryDirectory() as out, \
            mock.patch.object(flex_bond, "BondBlockStore", side_effect=spy), \
            mock.patch.object(type(s), "_bond_w_in_place", return_value=in_place):
        s.solve(gi, out)
    state = {k: np.array(v) for k, v in gi.items()
             if isinstance(v, np.ndarray) or k.startswith("longitudinal_bond_")}
    state["sigma"] = np.array(s.sigma) if getattr(s, "sigma", None) is not None else None
    return state, built, s


class TestBitIdenticalState(unittest.TestCase):
    """Two-orbital ``input_2orb`` (inter-orbital CoulombInter, several bond
    channels, the mixed-block pair permutation active)."""

    def _compare(self, extra):
        a, built_a, sa = _solve(extra, True)
        b, built_b, _ = _solve(extra, False)
        self.assertEqual(built_a, [{"W": "chibar"}])
        self.assertEqual(built_b, [{}])
        self.assertTrue(sa._bond_est_kwargs["w_in_place"])
        self.assertEqual(sorted(a), sorted(b))
        for k in a:
            if a[k] is None:
                self.assertIsNone(b[k], k)
                continue
            np.testing.assert_array_equal(a[k], b[k], err_msg=k)
        self.assertIn("green", a)
        self.assertIn("longitudinal_bond_chi_s", a)

    def test_default(self):
        self._compare({})

    def test_output_full(self):
        self._compare({"longitudinal_bond_output_full": True})


if __name__ == "__main__":
    unittest.main()
