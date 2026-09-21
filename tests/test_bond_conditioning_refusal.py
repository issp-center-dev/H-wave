"""The conditioning guard's refusal carries its location and criteria as
ATTRIBUTES (issue #198): ``bond_channels.BondConditioningError`` (a
``ValueError``, so every existing ``except ValueError`` still catches it)
names the channel, the q-point and the two criteria structurally, and
``dress_batch`` builds its refusal from the score it has already computed
instead of scoring the batch a second time. Both message texts are pinned
in full, unchanged from before."""

import unittest
from unittest import mock

import numpy as np

from hwave.solver import bond_channels as bc


def _refusing_static(sigma_min=1.0e-9, iq=1, Nx=3):
    """A (Nx, 1, 1, 2, 2) stack of I - chibar V blocks that is the identity
    everywhere except q = iq, where it is ``sigma_min * I``: both singular
    values equal, so the RELATIVE criterion is 1 and only the absolute pole
    distance ``sigma_min / max(1, sigma_max) = sigma_min`` trips."""
    mat = np.broadcast_to(np.eye(2, dtype=complex), (Nx, 1, 1, 2, 2)).copy()
    mat[iq, 0, 0] = sigma_min * np.eye(2)
    return mat


def _dress_bond_message(name, q, ratio, pole, cond_tol, smin, smax):
    """The pre-#198 text of the guard's refusal, verbatim."""
    return ("dress_bond: the {} RPA denominator is singular or nearly singular at "
            "the q-point index ({}, {}, {}): sigma_min/sigma_max = {:.3e}, "
            "sigma_min/max(1, sigma_max) = {:.3e}; the smaller of the two is "
            "<= cond_tol = {:.3e} (sigma_min = {:.3e}, sigma_max = {:.3e}). The "
            "bond path has entered the {} instability region, where the dressed "
            "vertices are enormous and numerically meaningless. Reduce the "
            "interaction strength, raise the temperature, refine/reduce the "
            "q-grid, or -- if you deliberately want to study the stiff regime -- "
            "lower cond_tol.".format(name, q[0], q[1], q[2], ratio, pole, cond_tol,
                                     smin, smax, name))


def _dress_batch_message(channel, nu, l, q):
    """The pre-#198 text of dress_batch's refusal, verbatim."""
    return ("dress_batch: the {} RPA denominator is singular or nearly singular at "
            "bosonic Matsubara index {} (grid index l={}) and q-point index ({}, {}, {}); "
            "the bond path has entered the instability region. Reduce the interaction "
            "strength, raise the temperature, refine or reduce the q grid, or lower "
            "cond_tol deliberately.".format(channel, nu, l, q[0], q[1], q[2]))


class TestCheckBondConditioning(unittest.TestCase):

    def test_refusal_is_structured_and_the_message_is_unchanged(self):
        mat = _refusing_static()
        score = bc.bond_conditioning_score(mat)
        worst, iq, ratio, pole, smin, smax = score
        with self.assertRaises(ValueError) as cm:          # still a ValueError
            bc._check_bond_conditioning("spin", mat, 1.0e-3)
        exc = cm.exception
        self.assertIsInstance(exc, bc.BondConditioningError)
        self.assertEqual(exc.channel, "spin")
        self.assertEqual(exc.iq, iq)
        self.assertEqual(exc.q, (1, 0, 0))
        self.assertIsNone(exc.l)
        self.assertEqual(exc.worst, worst)
        self.assertEqual(exc.ratio, ratio)
        self.assertEqual(exc.pole, pole)
        self.assertEqual(exc.smin, smin)
        self.assertEqual(exc.smax, smax)
        self.assertEqual(exc.cond_tol, 1.0e-3)
        self.assertEqual(str(exc), _dress_bond_message(
            "spin", (1, 0, 0), ratio, pole, 1.0e-3, smin, smax))

    def test_passing_the_score_skips_the_second_decomposition(self):
        mat = _refusing_static()
        score = bc.bond_conditioning_score(mat)

        def _no_svd(_mat):
            raise AssertionError("bond_conditioning_score must not run again")

        with mock.patch.object(bc, "bond_conditioning_score", _no_svd):
            with self.assertRaises(bc.BondConditioningError) as cm:
                bc._check_bond_conditioning("charge", mat, 1.0e-3, score=score)
            # and a passing score is returned as before
            ok = (0.5, 0, 0.5, 0.5, 0.5, 1.0)
            self.assertEqual(bc._check_bond_conditioning("charge", mat, 1.0e-3, score=ok), 0.5)
        self.assertEqual(cm.exception.iq, score[1])
        self.assertEqual(cm.exception.channel, "charge")

    def test_q_index_is_decoded_on_the_three_dimensional_grid(self):
        # (Nx, Ny, Nz) = (2, 3, 2): flattened iq = 7 -> (1, 0, 1)
        mat = np.broadcast_to(np.eye(2, dtype=complex), (2, 3, 2, 2, 2)).copy()
        mat[1, 0, 1] = 1.0e-9 * np.eye(2)
        with self.assertRaises(bc.BondConditioningError) as cm:
            bc._check_bond_conditioning("spin", mat, 1.0e-3)
        self.assertEqual(cm.exception.iq, 7)
        self.assertEqual(cm.exception.q, (1, 0, 1))
        self.assertIn("q-point index (1, 0, 1)", str(cm.exception))


class TestDressBatchRefusal(unittest.TestCase):

    def _batch(self, sigma_min=1.0e-9, l_bad=1, q_bad=1, nb=3, nvol=2):
        # chibar = c * I with W = I gives 1 - chibar W = (1 - c) I; put the
        # near-singular block at (l_bad, q_bad)
        cb = np.zeros((nb, nvol, 2, 2), dtype=complex)
        for l in range(nb):
            for q in range(nvol):
                c = 0.3 if (l, q) != (l_bad, q_bad) else 1.0 - sigma_min
                cb[l, q] = c * np.eye(2)
        W = np.broadcast_to(np.eye(2, dtype=complex), (nvol, 2, 2)).copy()
        return cb, W

    def test_refusal_names_the_location_structurally_from_one_score(self):
        cb, W = self._batch()                     # bad block at batch row 1, q = 1
        with mock.patch.object(bc, "bond_conditioning_score",
                               wraps=bc.bond_conditioning_score) as score:
            with self.assertRaises(ValueError) as cm:
                bc.dress_batch(cb, W, "spin", l0=4, nmat=16, spatial_shape=(2, 1, 1))
        self.assertEqual(score.call_count, 1)             # no second decomposition
        exc = cm.exception
        self.assertIsInstance(exc, bc.BondConditioningError)
        self.assertEqual(exc.channel, "spin")
        self.assertEqual(exc.l, 4 + 1)
        self.assertEqual(exc.iq, 1)
        self.assertEqual(exc.q, (1, 0, 0))
        self.assertEqual(exc.cond_tol, bc._BOND_COND_FLOOR)
        self.assertAlmostEqual(exc.pole, 1.0e-9, delta=1.0e-12)
        self.assertEqual(exc.worst, min(exc.ratio, exc.pole))
        self.assertEqual(str(exc), _dress_batch_message("spin", 2 * 5 - 16, 5, (1, 0, 0)))
        # the chained cause is the guard's own refusal from the same score,
        # located on the (n, 1, 1, ND, ND) blocks of THIS call: its q index
        # is the flattened batch position (row 1, q 1 of nvol 2 -> 3), not
        # a grid index, and it carries no l
        cause = exc.__cause__
        self.assertIsInstance(cause, bc.BondConditioningError)
        self.assertEqual(cause.channel, "spin")
        self.assertEqual(cause.iq, 3)
        self.assertEqual(cause.q, (3, 0, 0))
        self.assertIsNone(cause.l)
        self.assertEqual(cause.worst, exc.worst)
        self.assertEqual(cause.smin, exc.smin)
        self.assertEqual(str(cause), _dress_bond_message(
            "spin", (3, 0, 0), exc.ratio, exc.pole, exc.cond_tol, exc.smin, exc.smax))

    def test_static_guard_locates_the_slice_in_batch_coordinates(self):
        # guard_freqs = "static": only the l = nmat // 2 slice is scored;
        # with l0 = 6, nb = 3 and nmat = 16 that is batch row 2, scored on
        # its own (nvol, ND, ND) blocks with the row offset 2 * nvol = 4.
        # The outer refusal reports the grid location; the cause the LOCAL
        # index within the scored slice (q = 1, not 4 + 1).
        cb, W = self._batch(l_bad=2, q_bad=1)
        with self.assertRaises(bc.BondConditioningError) as cm:
            bc.dress_batch(cb, W, "spin", l0=6, nmat=16, spatial_shape=(2, 1, 1),
                           guard_freqs="static")
        exc = cm.exception
        self.assertEqual(exc.l, 8)
        self.assertEqual(exc.iq, 1)
        self.assertEqual(exc.q, (1, 0, 0))
        self.assertEqual(str(exc), _dress_batch_message("spin", 2 * 8 - 16, 8, (1, 0, 0)))
        self.assertEqual(exc.__cause__.iq, 1)
        self.assertEqual(exc.__cause__.q, (1, 0, 0))
        self.assertIsNone(exc.__cause__.l)

    def test_passing_batches_are_unaffected(self):
        cb, W = self._batch(sigma_min=0.5)
        chi, cond = bc.dress_batch(cb, W, "spin", l0=0, nmat=16, spatial_shape=(2, 1, 1))
        self.assertEqual(chi.shape, cb.shape)
        self.assertGreater(cond, bc._BOND_COND_FLOOR)

    def test_production_wrapper_keeps_the_structured_refusal(self):
        # flex_bond._dress appends the SCF iteration to a refusal on the
        # production path; the exception it raises must still be the
        # structured one, with the inner refusal chained
        from hwave.solver import flex_bond as fb
        cb, W = self._batch()
        with self.assertRaises(bc.BondConditioningError) as cm:
            fb._dress(cb, W, "spin", 4, 16, (2, 1, 1), bc._BOND_COND_FLOOR, 3)
        exc = cm.exception
        inner = exc.__cause__
        self.assertIsInstance(inner, bc.BondConditioningError)
        self.assertEqual(str(exc), str(inner) + fb._at(3))
        self.assertIn("SCF iteration 3", str(exc))
        for name in ("channel", "iq", "q", "l", "worst", "ratio", "pole", "smin", "smax", "cond_tol"):
            self.assertEqual(getattr(exc, name), getattr(inner, name), name)
        self.assertEqual((exc.l, exc.q), (5, (1, 0, 0)))


if __name__ == "__main__":
    unittest.main()
