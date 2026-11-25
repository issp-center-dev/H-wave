"""Labels and frames of ``hwave.initial`` (``hwave.initial.labels``; design
section 6.5).

A local pattern is decomposed ``O_{(as),(bt)} = 1/2 sum_mu c^mu_ab
sigma^mu_st`` (mu = 0 charge/orbital, x, y, z spin; the spin axes are those
of the input's spinor basis). ``label`` reads an induced ``delta rho`` in
the rho convention (``rho = G^T``) on the folded cell of the problem.
"""
import math
import os
import shutil
import tempfile
import time
import unittest

import numpy as np

from tests import initial_fixtures

_S0 = np.eye(2, dtype=complex)
_SX = np.array([[0, 1], [1, 0]], dtype=complex)
_SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
_SZ = np.array([[1, 0], [0, -1]], dtype=complex)
_SIG = {"c0": _S0, "cx": _SX, "cy": _SY, "cz": _SZ}

#: two orbitals on one site, spin-independent hopping on a chain
_TRANSFER2 = {((1, 0, 0), 0, 0): -1.0, ((-1, 0, 0), 0, 0): -1.0,
              ((1, 0, 0), 1, 1): -0.5, ((-1, 0, 0), 1, 1): -0.5,
              ((0, 0, 0), 0, 1): -0.2, ((0, 0, 0), 1, 0): -0.2}


def _so_table(transfer):
    """A 7-column table spread over spin in the interleaved ``2*orb + spin``
    order of spin-orbital mode."""
    return {(r, 2 * a + s, 2 * b + s): v for (r, a, b), v in transfer.items()
            for s in range(2)}


def _block(c):
    """The 2x2 spin block ``1/2 sum_mu c^mu sigma^mu`` of scalar ``c``."""
    return 0.5 * sum(complex(c.get(k, 0)) * _SIG[k] for k in _SIG)


def _rho(problem, entries):
    """A rho-convention array ``(nvol, ns, norb, ns, norb)`` of ``problem``
    from ``{(r, A, B): {"c0": .., "cz": ..}}`` (``A``, ``B`` folded PHYSICAL
    orbitals; the spin block of each entry is ``1/2 sum c^mu sigma^mu``)."""
    so = problem.enable_spin_orbital
    nphys = problem.nd // 2
    out = np.zeros((problem.nvol, nphys, 2, nphys, 2), dtype=complex)
    for (r, a, b), c in entries.items():
        out[r, a, :, b, :] += _block(c)
    if so:   # interleaved 2*A + s
        return out.reshape(problem.nvol, 1, 2 * nphys, 1, 2 * nphys)
    # spin-block s*norb + A
    return np.transpose(out, (0, 2, 1, 4, 3)).copy()


class _Tmp(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="initial_labels_")
        self._n = 0
        self._t0 = time.perf_counter()

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)
        # each fast test method stays under 5 s (tests/heavy_tests.py)
        self.assertLess(time.perf_counter() - self._t0, 5.0)

    def dest(self):
        self._n += 1
        return os.path.join(self.root, "m{}".format(self._n))

    def facts(self, toml):
        from hwave.initial import facts
        return facts.model_facts(toml)

    def two_orbital(self, *, transfer=None, spin_transfer=None, so=False, cellshape=(3, 1, 1)):
        """The facts of the two-orbital, one-site chain (``SubShape = [1,1,1]``)."""
        extra = None
        if spin_transfer is not None:
            extra = {"transfer.dat": initial_fixtures.spin_transfer_text(2, spin_transfer)}
        if so:
            toml = initial_fixtures.write_input(
                self.dest(), cellshape=cellshape, norb=4,
                transfer=_so_table(_TRANSFER2) if transfer is None else transfer,
                param={"Ncond": 2, "T": 0.1, "SubShape": [1, 1, 1]},
                mode={"enable_spin_orbital": True})
        else:
            toml = initial_fixtures.write_input(
                self.dest(), cellshape=cellshape, norb=2,
                transfer=_TRANSFER2 if transfer is None else transfer,
                param={"Ncond": 2, "T": 0.1, "SubShape": [1, 1, 1]}, extra_files=extra)
        return self.facts(toml)


# --------------------------------------------------------------------------
# decompose and tr_parity
# --------------------------------------------------------------------------

class TestDecompose(unittest.TestCase):
    def test_interleaved_round_trip(self):
        from hwave.initial import labels
        rng = np.random.default_rng(7)
        c = {k: rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)) for k in _SIG}
        block = 0.5 * sum(np.kron(c[k], _SIG[k]) for k in _SIG)   # interleaved (a, s)
        got = labels.decompose(block)
        self.assertEqual(sorted(got), ["c0", "cx", "cy", "cz"])
        for k in _SIG:
            np.testing.assert_allclose(got[k], c[k], atol=1e-14)
        np.testing.assert_allclose(labels.compose(got), block, atol=1e-14)

    def test_spin_block_layout(self):
        from hwave.initial import labels
        rng = np.random.default_rng(8)
        c = {k: rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)) for k in _SIG}
        block = 0.5 * sum(np.kron(_SIG[k], c[k]) for k in _SIG)   # spin-block (s, a)
        got = labels.decompose(block, layout="spin_block")
        for k in _SIG:
            np.testing.assert_allclose(got[k], c[k], atol=1e-14)
        np.testing.assert_allclose(labels.compose(got, layout="spin_block"), block,
                                   atol=1e-14)

    def test_moment_reading(self):
        """``c^mu = tr_s(sigma^mu rho)`` is twice the spin moment: a
        polarised up electron has ``c^0 = c^z = 1``."""
        from hwave.initial import labels
        rho = np.diag([1.0, 0.0]).astype(complex)          # one orbital, up
        got = labels.decompose(rho)
        self.assertAlmostEqual(complex(got["c0"][0, 0]), 1.0)
        self.assertAlmostEqual(complex(got["cz"][0, 0]), 1.0)
        self.assertAlmostEqual(abs(got["cx"][0, 0]) + abs(got["cy"][0, 0]), 0.0)

    def test_bad_shape(self):
        from hwave.initial import labels
        with self.assertRaises(ValueError):
            labels.decompose(np.zeros((3, 3)))
        with self.assertRaises(ValueError):
            labels.decompose(np.zeros((4, 2)))
        with self.assertRaises(ValueError):
            labels.decompose(np.zeros((4, 4)), layout="interleave")


class TestTrParity(unittest.TestCase):
    def test_tr_parity_table(self):
        """Design section 6.5: TR-even = real-symmetric c^0 and
        imaginary-antisymmetric c^i; TR-odd = the converse; mixed
        otherwise."""
        from hwave.initial import labels
        sym = np.array([[1.0, 0.3], [0.3, -0.5]])
        anti = np.array([[0.0, 0.4], [-0.4, 0.0]])
        z = np.zeros((2, 2))
        cases = [
            ({"c0": sym}, "even"),
            ({"c0": 1j * anti}, "odd"),                    # orbital current
            ({"cz": 1j * anti}, "even"),                   # SOC-like i lambda sigma^z
            ({"cx": 1j * anti, "cy": 1j * anti}, "even"),
            ({"cz": sym}, "odd"),                          # a moment
            ({"cy": sym}, "odd"),
            ({"c0": sym, "cz": 1j * anti}, "even"),
            ({"c0": sym, "cz": sym}, "mixed"),
            ({"c0": 1j * anti, "cz": 1j * anti}, "mixed"),
        ]
        for c, want in cases:
            full = {k: c.get(k, z) for k in _SIG}
            with self.subTest(c=sorted(c), want=want):
                self.assertEqual(labels.tr_parity(full), want)
                # missing keys are zero
                self.assertEqual(labels.tr_parity(c), want)
        # a block whose coefficients are complex in a TR-mixed way
        self.assertEqual(labels.tr_parity({"c0": np.array([[0.0, 1 + 1j], [1 - 1j, 0.0]])}),
                         "mixed")

    def test_tr_minor_weight_rule(self):
        """O8: even (odd) when the other parity carries less than
        TR_MINOR_WEIGHT = 0.1 of the squared norm, mixed otherwise."""
        from hwave.initial import labels
        self.assertEqual(labels.TR_MINOR_WEIGHT, 0.1)
        one = np.ones((1, 1))
        for w_odd, want in ((0.0, "even"), (0.09, "even"), (0.11, "mixed"), (0.89, "mixed"),
                            (0.91, "odd"), (1.0, "odd")):
            c = {"c0": math.sqrt(1 - w_odd) * one, "cz": math.sqrt(w_odd) * one}
            with self.subTest(w_odd=w_odd):
                self.assertEqual(labels.tr_parity(c), want)

    def test_tr_parity_is_time_reversal(self):
        """The parity is that of the operator under ``Theta = i sigma^y K``:
        ``Y conj(O) Y^T = +O`` (even) or ``-O`` (odd)."""
        from hwave.initial import labels
        y = np.kron(np.eye(2), np.array([[0, 1], [-1, 0]]))   # interleaved
        rng = np.random.default_rng(3)
        for _ in range(4):
            a = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
            o = a + a.conj().T
            for part, want in (((o + y @ o.conj() @ y.T) / 2, "even"),
                               ((o - y @ o.conj() @ y.T) / 2, "odd")):
                self.assertEqual(labels.tr_parity(labels.decompose(part)), want)
            self.assertEqual(labels.tr_parity(labels.decompose(o)), "mixed")


# --------------------------------------------------------------------------
# gauge
# --------------------------------------------------------------------------

class TestGauge(_Tmp):
    def test_gauge_is_cached_per_problem(self):
        """The Transfer table is decomposed once per problem, not once per
        label."""
        from unittest import mock
        from hwave.initial import labels
        f = self.two_orbital()
        calls = []
        real = labels.decompose

        def counting(*a, **kw):
            calls.append(np.shape(a[0]))
            return real(*a, **kw)
        with mock.patch.object(labels, "decompose", counting):
            for _ in range(3):
                self.assertFalse(labels.gauge_dependent(f.problem))
                self.assertFalse(labels.gauge_dependent(f))
        self.assertEqual(calls.count(f.problem.ham_trans_r.shape), 1)

    def test_real_hopping_no_gauge_warning(self):
        from hwave.initial import labels
        f = self.two_orbital()
        self.assertFalse(labels.gauge_dependent(f.problem))
        self.assertFalse(labels.gauge_dependent(f))

    def test_soc_no_gauge_warning(self):
        """V5: inter-orbital SOC ``i lambda sigma^z`` is an imaginary c^z,
        a spin part, exempt from the gauge warning, in spin-orbital mode
        and as a 9-column table in normal mode."""
        from hwave.initial import labels
        lam = 0.25
        so = _so_table(_TRANSFER2)
        for s, sign in ((0, 1), (1, -1)):
            so[((0, 0, 0), 0 + s, 2 + s)] = 1j * lam * sign
            so[((0, 0, 0), 2 + s, 0 + s)] = -1j * lam * sign
        f = self.two_orbital(transfer=so, so=True)
        self.assertTrue(f.problem.enable_spin_orbital)
        self.assertFalse(labels.gauge_dependent(f.problem))
        # an imaginary spin flip on a bond (orbital 0 up -> 0 down, the
        # interleaved slots 0 and 1) is a spin part too; read in the
        # spin-block order these two slots would be one spin
        so = _so_table(_TRANSFER2)
        so[((1, 0, 0), 0, 1)] = 1j * lam
        so[((-1, 0, 0), 1, 0)] = -1j * lam
        f = self.two_orbital(transfer=so, so=True)
        self.assertFalse(labels.gauge_dependent(f.problem))
        entries = initial_fixtures.spin_independent_entries(_TRANSFER2)
        entries[((0, 0, 0), 0, 1, 0, 0)] = 1j * lam
        entries[((0, 0, 0), 0, 1, 1, 1)] = -1j * lam
        f = self.two_orbital(spin_transfer=entries)
        self.assertTrue(f.problem.spin_resolved_transfer)
        self.assertFalse(labels.gauge_dependent(f.problem))
        # a spin-flip SOC i lambda sigma^x on a bond is exempt as well
        entries = initial_fixtures.spin_independent_entries(_TRANSFER2)
        entries[((1, 0, 0), 0, 1, 0, 1)] = 1j * lam
        entries[((1, 0, 0), 0, 1, 1, 0)] = 1j * lam
        f = self.two_orbital(spin_transfer=entries, cellshape=(3, 1, 1))
        self.assertFalse(labels.gauge_dependent(f.problem))

    def test_complex_hopping_gauge_warning(self):
        """A complex spin-independent hopping (``e^{0.3i}`` on a bond) marks
        the labels gauge-dependent: in a 7-column table, a 9-column table
        and in spin-orbital mode."""
        from hwave.initial import labels
        cplx = dict(_TRANSFER2)
        cplx[((1, 0, 0), 0, 0)] = -np.exp(0.3j)
        cplx[((-1, 0, 0), 0, 0)] = -np.exp(-0.3j)
        f = self.two_orbital(transfer=cplx, cellshape=(3, 1, 1))
        self.assertTrue(labels.gauge_dependent(f.problem))
        self.assertTrue(labels.gauge_dependent(f))
        f = self.two_orbital(transfer=_so_table(cplx), so=True, cellshape=(3, 1, 1))
        self.assertTrue(labels.gauge_dependent(f.problem))
        entries = initial_fixtures.spin_independent_entries(cplx)
        f = self.two_orbital(spin_transfer=entries, cellshape=(3, 1, 1))
        self.assertTrue(f.problem.spin_resolved_transfer)
        self.assertTrue(labels.gauge_dependent(f.problem))

    def test_gauge_threshold_is_the_spin_average(self):
        """R5: the spin-independent part of a 9-column entry is
        ``(T_upup + T_dndn)/2``. An imaginary part only on one spin (a
        mixture of c^0 and c^z) counts by half of it; ``1e-8`` is the
        threshold on that half."""
        from hwave.initial import labels
        for im, want in ((3e-8, True), (1.5e-8, False)):
            entries = initial_fixtures.spin_independent_entries(_TRANSFER2)
            entries[((0, 0, 0), 0, 1, 0, 0)] = -0.2 + 1j * im
            with self.subTest(im=im):
                f = self.two_orbital(spin_transfer=entries)
                self.assertIs(labels.gauge_dependent(f.problem), want)

    def test_label_carries_the_flag(self):
        from hwave.initial import labels
        cplx = dict(_TRANSFER2)
        cplx[((1, 0, 0), 0, 0)] = -np.exp(0.3j)
        cplx[((-1, 0, 0), 0, 0)] = -np.exp(-0.3j)
        for transfer, want in ((None, False), (cplx, True)):
            f = self.two_orbital(transfer=transfer, cellshape=(3, 1, 1))
            rho = _rho(f.problem, {(0, 0, 0): {"cz": 1.0}, (0, 1, 1): {"cz": 1.0}})
            lab = labels.label(rho, f, q=("0", "0", "0"), eta=None, phase=0.0, sign=1,
                               spiral=False)
            self.assertIs(lab.gauge_dependent, want)


# --------------------------------------------------------------------------
# label
# --------------------------------------------------------------------------

class TestLabel(_Tmp):
    def square_facts(self):
        toml = initial_fixtures.square(self.dest(), 4, U=4.0, T=0.1, subshape=(2, 2, 1))
        return self.facts(toml)

    def neel_rho(self, problem, m=0.4, axis="cz"):
        """A Neel delta rho on the 2x2 folded cell of the square lattice:
        folded orbital ``ir = ix + 2 iy`` carries ``(-1)^(ix+iy) m``."""
        entries = {}
        for iy in range(2):
            for ix in range(2):
                entries[(0, ix + 2 * iy, ix + 2 * iy)] = {axis: m * (-1) ** (ix + iy)}
        return _rho(problem, entries)

    def test_example_name(self):
        """The example of design section 6.5."""
        from hwave.initial import labels
        f = self.square_facts()
        self.assertEqual(f.problem.nvol, 4)
        lab = labels.label(self.neel_rho(f.problem), f, q=("1/2", "1/2", "0"), eta="Néel",
                           phase=0.0, sign=1, spiral=False)
        self.assertEqual(lab.name, "spin-z ferro-orbital q=(1/2,1/2,0) Néel φ=0")
        self.assertAlmostEqual(lab.weights["spin-z"], 1.0, places=14)
        self.assertEqual(lab.tr, "odd")
        self.assertFalse(lab.gauge_dependent)
        # every weight is listed; the channel weights sum to 1
        for key in labels.CHANNELS:
            self.assertIn(key, lab.weights)
        self.assertAlmostEqual(sum(lab.weights[k] for k in labels.CHANNELS), 1.0, places=14)
        self.assertAlmostEqual(lab.weights["TR-odd"], 1.0, places=14)
        self.assertAlmostEqual(lab.weights["TR-even"], 0.0, places=14)

    def test_spin_axes_and_name_parts(self):
        from hwave.initial import labels
        from fractions import Fraction
        f = self.square_facts()
        for axis in ("cx", "cy", "cz"):
            lab = labels.label(self.neel_rho(f.problem, axis=axis), f,
                               q=(Fraction(1, 2), Fraction(1, 2), 0), eta=None,
                               phase=math.pi / 4, sign=-1, spiral=False)
            self.assertEqual(lab.name, "spin-{} ferro-orbital q=(1/2,1/2,0) φ=π/4 sign=-1"
                             .format(axis[1]))
        lab = labels.label(self.neel_rho(f.problem), f, q=None, eta=(1,), phase=None,
                           sign=None, spiral=True)
        self.assertEqual(lab.name, "spin-z ferro-orbital spiral")
        lab = labels.label(self.neel_rho(f.problem), f, q=("1/2", "1/2", "0"), eta=(1, -1),
                           phase=-math.pi / 2, sign=1, spiral=False)
        self.assertEqual(lab.name, "spin-z ferro-orbital q=(1/2,1/2,0) η=(+,-) φ=-π/2")
        lab = labels.label(self.neel_rho(f.problem), f, q=("0", "0", "0"), eta=None,
                           phase=0.1234567, sign=1, spiral=False)
        self.assertEqual(lab.name, "spin-z ferro-orbital q=(0,0,0) φ=0.123457")
        lab = labels.label(self.neel_rho(f.problem), f, q=("0", "0", "0"), eta=None,
                           phase=math.pi, sign=1, spiral=False)
        self.assertTrue(lab.name.endswith(" φ=π"))

    def test_rho_convention(self):
        """``label`` reads rho = G^T: a stored G slot (the transpose) of a
        sigma^y state is the -sigma^y state; the weights agree, and the
        transposed bond reads the reverse direction."""
        from hwave.initial import labels
        f = self.square_facts()
        rho = self.neel_rho(f.problem, axis="cy")
        nd = f.problem.nd
        g = np.swapaxes(rho.reshape(4, nd, nd), -1, -2).reshape(rho.shape)
        a = labels.label(rho, f, q=None, eta=None, phase=None, sign=None, spiral=False)
        b = labels.label(g, f, q=None, eta=None, phase=None, sign=None, spiral=False)
        self.assertEqual(a.name, b.name)
        # (nvol, nd, nd) is accepted as well
        c = labels.label(rho.reshape(4, nd, nd), f, q=None, eta=None, phase=None,
                         sign=None, spiral=False)
        self.assertEqual(a.weights, c.weights)

    def test_bond_label_prefix(self):
        """Entries between two sites (in the cell or across cells) are
        labelled as bonds: ``bond-...``."""
        from hwave.initial import labels
        f = self.square_facts()
        # a real bond order between folded sites 0 and 1 (in the cell)
        rho = _rho(f.problem, {(0, 0, 1): {"c0": 0.3}, (0, 1, 0): {"c0": 0.3}})
        lab = labels.label(rho, f, q=("0", "0", "0"), eta=None, phase=0.0, sign=1,
                           spiral=False)
        self.assertTrue(lab.name.startswith("bond-"), lab.name)
        self.assertEqual(lab.name, "bond-charge q=(0,0,0) φ=0")
        self.assertEqual(lab.tr, "even")
        # an imaginary (current) bond across supercells (r != 0; on the 2x2
        # supercell mesh -r is r, so the Hermitian partner sits at the same r)
        rho = _rho(f.problem, {(1, 0, 1): {"c0": 0.3j}, (1, 1, 0): {"c0": -0.3j}})
        lab = labels.label(rho, f, q=None, eta=None, phase=None, sign=None, spiral=False)
        self.assertEqual(lab.name, "bond-current")
        self.assertEqual(lab.tr, "odd")
        # a same-orbital entry at r != 0 is a bond too
        rho = _rho(f.problem, {(1, 0, 0): {"cz": 0.2}})
        lab = labels.label(rho, f, q=None, eta=None, phase=None, sign=None, spiral=False)
        self.assertEqual(lab.name, "bond-spin-z")
        rho = _rho(f.problem, {(1, 0, 1): {"cx": 0.2j}, (1, 1, 0): {"cx": -0.2j}})
        lab = labels.label(rho, f, q=None, eta=None, phase=None, sign=None, spiral=False)
        self.assertEqual(lab.name, "bond-spin-orbital-even")
        self.assertEqual(lab.tr, "even")

    def test_on_site_channel_names(self):
        """The on-site channels of section 6.5 on a two-orbital site."""
        from hwave.initial import labels
        f = self.two_orbital()
        cases = [
            ({(0, 0, 0): {"c0": 0.2}, (0, 1, 1): {"c0": 0.2}}, "charge", "even"),
            ({(0, 0, 0): {"c0": 0.2}, (0, 1, 1): {"c0": -0.2}}, "orbital-polarization",
             "even"),
            ({(0, 0, 1): {"c0": 0.2}, (0, 1, 0): {"c0": 0.2}}, "orbital-hybridization",
             "even"),
            ({(0, 0, 1): {"c0": 0.2j}, (0, 1, 0): {"c0": -0.2j}}, "orbital-current", "odd"),
            ({(0, 0, 0): {"cz": 0.2}, (0, 1, 1): {"cz": 0.2}}, "spin-z ferro-orbital", "odd"),
            ({(0, 0, 0): {"cx": 0.2}, (0, 1, 1): {"cx": -0.2}}, "spin-x antiferro-orbital",
             "odd"),
            ({(0, 0, 1): {"cy": 0.2}, (0, 1, 0): {"cy": 0.2}}, "inter-orbital-spin-y", "odd"),
            ({(0, 0, 1): {"cz": 0.2j}, (0, 1, 0): {"cz": -0.2j}}, "spin-orbital-even", "even"),
        ]
        for entries, name, tr in cases:
            with self.subTest(name=name):
                lab = labels.label(_rho(f.problem, entries), f, q=None, eta=None, phase=None,
                                   sign=None, spiral=False)
                self.assertEqual(lab.name, name)
                self.assertEqual(lab.tr, tr)

    def test_dominant_weights(self):
        """Channels with weight >= 0.1 name the label, heaviest first; the
        rest are only listed. A TR-mixed state says so."""
        from hwave.initial import labels
        f = self.two_orbital()
        # spin-z weight 2 * 0.3^2, orbital-polarization 2 * 0.13^2 (about
        # 0.15, just above the threshold), charge 2 * 0.05^2
        entries = {(0, 0, 0): {"cz": 0.3, "c0": 0.18}, (0, 1, 1): {"cz": 0.3, "c0": -0.08}}
        lab = labels.label(_rho(f.problem, entries), f, q=None, eta=None, phase=None,
                           sign=None, spiral=False)
        total = 2 * 0.3 ** 2 + 2 * 0.13 ** 2 + 2 * 0.05 ** 2
        self.assertAlmostEqual(lab.weights["spin-z"], 2 * 0.09 / total, places=14)
        self.assertAlmostEqual(lab.weights["orbital-polarization"], 2 * 0.0169 / total,
                               places=14)
        self.assertLess(lab.weights["orbital-polarization"], 0.16)
        self.assertAlmostEqual(lab.weights["charge"], 2 * 0.0025 / total, places=14)
        self.assertLess(lab.weights["charge"], labels.DOMINANT_WEIGHT)
        self.assertEqual(lab.name, "spin-z ferro-orbital + orbital-polarization")
        self.assertEqual(lab.tr, "mixed")
        self.assertAlmostEqual(lab.weights["TR-odd"], 2 * 0.09 / total, places=14)
        self.assertAlmostEqual(lab.weights["TR-even"] + lab.weights["TR-odd"], 1.0,
                               places=14)

    def test_spin_orbital_mode_reads_physical_orbitals(self):
        """R8: in spin-orbital mode the label decomposes the interleaved
        ``2*orb + spin`` layout into physical orbitals and spin; the same
        physical delta rho gives the same label as in normal mode."""
        from hwave.initial import labels
        fn = self.two_orbital()
        fs = self.two_orbital(so=True)
        self.assertEqual(fs.problem.nd, fn.problem.nd)
        for entries in ({(0, 0, 1): {"cz": 0.2j}, (0, 1, 0): {"cz": -0.2j}},
                        {(0, 0, 0): {"cx": 0.2}, (0, 1, 1): {"cx": -0.2, "c0": 0.05}},
                        {(1, 0, 1): {"cy": 0.1}, (2, 1, 0): {"cy": 0.1}}):
            a = labels.label(_rho(fn.problem, entries), fn, q=None, eta=None, phase=None,
                             sign=None, spiral=False)
            b = labels.label(_rho(fs.problem, entries), fs, q=None, eta=None, phase=None,
                             sign=None, spiral=False)
            self.assertEqual(a.name, b.name)
            self.assertEqual(a.tr, b.tr)
            for k in a.weights:
                self.assertAlmostEqual(a.weights[k], b.weights[k], places=14)

    def test_zero_delta_rho(self):
        """A zero or tolerance-small delta rho (the rounding noise of a
        converged para run) is ``para``; above PARA_TOL it is named."""
        from hwave.initial import labels
        f = self.two_orbital()
        noise = {(0, 0, 0): {"cz": 3e-12, "c0": 1e-12}, (0, 1, 1): {"cx": -2e-12}}
        for entries in ({}, noise):
            lab = labels.label(_rho(f.problem, entries), f, q=None, eta=None, phase=None,
                               sign=None, spiral=False)
            self.assertEqual(lab.name, "para")
            self.assertEqual(lab.tr, "even")
            self.assertTrue(all(v == 0.0 for v in lab.weights.values()))
        self.assertEqual(labels.PARA_TOL, 1e-8)
        lab = labels.label(_rho(f.problem, {(0, 0, 0): {"cz": 1e-6}}), f, q=None, eta=None,
                           phase=None, sign=None, spiral=False)
        self.assertEqual(lab.name, "spin-z ferro-orbital")

    def test_wrong_shape_refused(self):
        from hwave.initial import labels
        f = self.two_orbital()
        with self.assertRaises(ValueError):
            labels.label(np.zeros((2, 4, 4)), f, q=None, eta=None, phase=None, sign=None,
                         spiral=False)
        with self.assertRaises(ValueError):
            labels.label(np.zeros((3, 4, 4)), f, q=None, eta=(1, 0), phase=None, sign=None,
                         spiral=False)
        with self.assertRaises(ValueError):
            labels.label(np.zeros((3, 4, 4)), f, q=None, eta=None, phase=None, sign=2,
                         spiral=False)
        with self.assertRaises(ValueError):
            labels.label(np.zeros((3, 4, 4)), f, q=(0.5, 0.5, 0), eta=None, phase=None,
                         sign=None, spiral=False)
        # a non-Hermitian delta rho (a bond without its partner) is refused
        with self.assertRaises(ValueError):
            labels.label(_rho(f.problem, {(1, 0, 1): {"c0": 0.3}}), f, q=None, eta=None,
                         phase=None, sign=None, spiral=False)


class TestPatternAndInducedNames(_Tmp):
    """Owner decision O8: a candidate is named by its seeding pattern O; an
    induced or converged state by its dominant components (>= 10 %), on-site
    first, with TR even/odd unless the other parity holds >= 10 %."""

    def setUp(self):
        super().setUp()
        toml = initial_fixtures.square(self.dest(), 4, T=0.1, subshape=(2, 2, 1))
        from hwave.initial import facts
        self.f = facts.model_facts(toml)
        p = self.f.problem
        eta = np.array([(-1) ** (ix + iy) for iy in range(2) for ix in range(2)])
        o = np.kron(np.diag([1.0, -1.0]), np.diag(eta)).astype(complex)  # spin-block
        self.o_k = np.broadcast_to(o, (p.nvol, p.nd, p.nd)).copy()

    def induced(self, h):
        """delta rho = rho(T - h O) - rho(T) by the solver's own fill
        (rho = G^T per displacement)."""
        from hwave.solver import initial_facade as fac
        p = self.f.problem

        def rho(hh):
            g = fac.green_function(p.ham_trans_k - hh * self.o_k, p.T, p.blocks,
                                   shape=p.shape, ns=p.ns, norb=p.norb,
                                   ene_cutoff=p.ene_cutoff)
            return np.swapaxes(g.reshape(p.nvol, p.nd, p.nd), -1, -2)
        return rho(h) - rho(0.0)

    def kw(self):
        return dict(q=("1/2", "1/2", "0"), eta="Néel", phase=0.0, sign=1, spiral=False)

    def test_pattern_gives_the_example_name(self):
        from hwave.initial import labels
        lab = labels.pattern_label(self.o_k, self.f, **self.kw())
        self.assertEqual(lab.name, "spin-z ferro-orbital q=(1/2,1/2,0) Néel φ=0")
        self.assertEqual(lab.ascii_name, "spin-z ferro-orbital q=(1/2,1/2,0) Neel phi=0")
        self.assertEqual(lab.source, "pattern")
        self.assertEqual(lab.tr, "odd")
        self.assertAlmostEqual(lab.weights["spin-z"], 1.0, places=14)
        self.assertEqual(lab.thresholds, labels.THRESHOLDS)
        self.assertEqual(lab.thresholds["dominant_weight"], 0.1)
        self.assertEqual(lab.thresholds["tr_minor_weight"], 0.1)
        # a non-Hermitian pattern is refused
        bad = self.o_k.copy()
        bad[:, 0, 1] = 0.3
        with self.assertRaises(ValueError):
            labels.pattern_label(bad, self.f, **self.kw())

    def test_induced_neel_probe(self):
        """The review's probe: the induced Neel delta rho has a bond tail of
        about half the weight. The name leads with the on-site channel;
        bond channels follow only at >= 10 %; at small and at finite h the
        parity is odd, not mixed."""
        from hwave.initial import labels
        for h in (1e-6, 1.0):
            with self.subTest(h=h):
                lab = labels.label(self.induced(h), self.f, **self.kw())
                self.assertEqual(lab.source, "state")
                self.assertTrue(lab.name.startswith("spin-z ferro-orbital"), lab.name)
                self.assertTrue(lab.name.endswith(" q=(1/2,1/2,0) Néel φ=0"), lab.name)
                head = lab.name[:-len(" q=(1/2,1/2,0) Néel φ=0")]
                listed = head.split(" + ")[1:]
                want = [ch for ch in labels.CHANNELS
                        if ch.startswith("bond-") and lab.weights[ch] >= 0.1]
                self.assertEqual(sorted(listed), sorted(want))
                for ch in labels.CHANNELS:
                    if ch.startswith("bond-") and ch not in want:
                        self.assertLess(lab.weights[ch], 0.1)
                self.assertGreater(lab.weights["spin-z"], 0.4)
                self.assertGreater(lab.weights["bond-spin-z"], 0.1)
                self.assertEqual(lab.tr, "odd")
                self.assertLess(lab.weights["TR-even"], labels.TR_MINOR_WEIGHT)
        # h = 1: the TR-even bond-charge part (about 1 %) is below the
        # threshold, so it is neither listed nor makes the parity mixed
        lab = labels.label(self.induced(1.0), self.f, **self.kw())
        self.assertGreater(lab.weights["bond-charge"], 1e-3)
        self.assertNotIn("bond-charge", lab.name)
        # h = 4 (of order the band width): the even part reaches 10 %
        lab = labels.label(self.induced(4.0), self.f, **self.kw())
        self.assertGreaterEqual(lab.weights["TR-even"], labels.TR_MINOR_WEIGHT)
        self.assertEqual(lab.tr, "mixed")
        self.assertTrue(lab.name.startswith("spin-z ferro-orbital + bond-charge"), lab.name)

    def test_bond_only_state_is_a_bond_label(self):
        """With no on-site channel at 10 % the name is the bond channels."""
        from hwave.initial import labels
        rho = _rho(self.f.problem, {(0, 0, 1): {"c0": 0.3}, (0, 1, 0): {"c0": 0.3},
                                    (0, 0, 0): {"cz": 0.02}})
        lab = labels.label(rho, self.f, q=None, eta=None, phase=None, sign=None,
                           spiral=False)
        self.assertEqual(lab.name, "bond-charge")

    def test_ascii_name(self):
        from hwave.initial import labels
        self.assertEqual(labels.ascii_name("spin-z q=(1/2,0,0) η=(+,-) φ=-π/2 Néel"),
                         "spin-z q=(1/2,0,0) eta=(+,-) phi=-pi/2 Neel")
        self.assertEqual(labels.ascii_name("a☃b"), "ab")
        lab = labels.label(self.induced(1e-3), self.f, q=None, eta=(1, -1),
                           phase=math.pi / 3, sign=None, spiral=False)
        self.assertTrue(lab.ascii_name.isascii())
        self.assertIn("eta=(+,-) phi=pi/3", lab.ascii_name)


if __name__ == "__main__":
    unittest.main()
