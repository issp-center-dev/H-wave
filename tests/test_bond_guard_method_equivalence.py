"""dress_batch(guard_method="interval") must reproduce guard_method="svd"
on every dressing fixture the suite already has (issue #197, spec 6.2):
the dressed batch, the published minimum and the refusal / warn outputs."""
import unittest

import numpy as np

from hwave.solver import bond_channels as bc
from tests.test_flex_bond_dressing import _problem


def _run(cb, W, channel, method, **kw):
    v = []
    try:
        chi, cond = bc.dress_batch(cb, W, channel, guard_method=method, violations=v, **kw)
        return ("ok", chi, cond, tuple(sorted((k, r[k]) for r in v for k in sorted(r))))
    except bc.BondConditioningError as e:
        return ("refused", str(e), e.channel, e.l, e.iq, e.q, e.worst, e.ratio, e.pole, e.smin, e.smax)


class TestGuardMethodEquivalence(unittest.TestCase):

    def _check(self, cb, W, channel, **kw):
        a = _run(cb, W, channel, "svd", **kw)
        b = _run(cb, W, channel, "interval", **kw)
        self.assertEqual(a[0], b[0])
        if a[0] == "ok":
            np.testing.assert_array_equal(a[1], b[1])
            self.assertEqual(a[2], b[2])
            self.assertEqual(a[3], b[3])
        else:
            self.assertEqual(a[1:], b[1:])

    def test_problem_fixture_both_channels_both_modes(self):
        chi_bar, S, C = _problem()
        nmat, nvol = chi_bar.shape[:2]
        for channel, W in (("spin", S), ("charge", C)):
            for guard_freqs in ("all", "static"):
                for l0, nb in ((0, nmat), (1, 3), (nmat // 2, 1)):
                    with self.subTest(channel=channel, guard_freqs=guard_freqs, l0=l0, nb=nb):
                        self._check(chi_bar[l0:l0 + nb], W, channel, l0=l0, nmat=nmat,
                                    spatial_shape=(nvol, 1, 1), guard_freqs=guard_freqs)

    def test_scaled_up_interaction_refuses_identically(self):
        chi_bar, S, C = _problem()
        nmat, nvol = chi_bar.shape[:2]
        for scale in (4.0, 8.0, 16.0):
            for channel, W in (("spin", S), ("charge", C)):
                for policy in ("refuse", "warn"):
                    with self.subTest(scale=scale, channel=channel, policy=policy):
                        self._check(chi_bar, scale * W, channel, l0=0, nmat=nmat,
                                    spatial_shape=(nvol, 1, 1), guard_policy=policy)

    def test_ties_across_l_and_q(self):
        chi_bar, S, _C = _problem(nmat=4, nvol=3)
        nmat, nvol = chi_bar.shape[:2]
        chi_bar[1, 2] = chi_bar[3, 0]                  # two identical blocks: tie in l and q
        self._check(chi_bar, S, "spin", l0=0, nmat=nmat, spatial_shape=(nvol, 1, 1))
        chi_bar[1, 0] = chi_bar[1, 2]                  # tie in q only
        self._check(chi_bar, S, "spin", l0=0, nmat=nmat, spatial_shape=(nvol, 1, 1))
