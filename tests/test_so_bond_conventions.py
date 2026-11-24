"""Convention pins of the spin-orbital bond-resolved vertex.

The production bond bubble (``bubble.bond_bubble_static`` /
``bubble._iter_bond_dynamic``) and the bond pairing kernel
(``eliashberg_bond.BondPairKernel``) are compared with brute-force
label-form evaluations (``tests/so_bond_oracles.py``). The bubble is scanned
over every candidate slot / momentum / frequency / phase reading; exactly one
must match, and it must be the convention object ``so_bond.SLOT_CONVENTION``.
The production spinful bond vertex (``so_bond.so_bond_vertex``) is checked
against the first-order vertex part of exact-diagonalization bond
correlators.
"""

import functools
import itertools
import re
import unittest
from types import SimpleNamespace

import numpy as np

from hwave.solver import bubble
from hwave.solver.bond_channels import BondSetView, BondTopology
from hwave.solver.eliashberg_bond import ArrayBlockSource, BondPairKernel, PairVertexUniform
from hwave.solver.so_bond import SLOT_CONVENTION, BondSlotGeometry, SlotConvention, slot_geometry
from tests.heavy_tests import heavy
from tests.so_bond_fixtures import production_first_order
from tests.so_bond_oracles import (CONVENTION_BETA, CONVENTION_DR, CONVENTION_GREEN_SEEDS,
                                   CONVENTION_ND, CONVENTION_NMAT_DYNAMIC,
                                   CONVENTION_NMAT_STATIC, CONVENTION_REV,
                                   CONVENTION_SAMPLE_SEED, CONVENTION_SAMPLE_SIZE,
                                   CONVENTION_SHAPE, FIXTURE_I, FIXTURE_II, FIXTURE_X,
                                   NEEDS_RICHARDSON, TYPES, candidate_map, fixture_name,
                                   lindhard_discrete, mesh_add, mesh_momenta, mesh_negate,
                                   random_green, rel, sketch_set, vertex_part_sketches)

DR, REV = CONVENTION_DR, CONVENTION_REV
NVOL = int(np.prod(CONVENTION_SHAPE))
B = DR.shape[0]
NC = B * CONVENTION_ND * CONVENTION_ND
GEOM = SimpleNamespace(delta_r=DR, reverse=REV, spatial_shape=CONVENTION_SHAPE)
FIXTURES = CONVENTION_GREEN_SEEDS
#: two candidate readings are one class when their scaled predictions agree to this
#: (round-off of identical sums taken in a different order; distinct readings differ by O(1))
SIGNATURE_ATOL = 1e-13


# ---------------------------------------------------------------------------
# Production bubble and label bubble per fixture (process-local caches)
# ---------------------------------------------------------------------------

def _nmat(static):
    return CONVENTION_NMAT_STATIC if static else CONVENTION_NMAT_DYNAMIC


@functools.lru_cache(maxsize=None)
def _green(seed, static):
    return random_green(seed, _nmat(static), NVOL, CONVENTION_ND)


@functools.lru_cache(maxsize=None)
def _label(seed, static):
    return lindhard_discrete(_green(seed, static), CONVENTION_BETA, GEOM)


@functools.lru_cache(maxsize=None)
def _production(seed, static):
    """The production code bubble: ``(nvol, NC, NC)`` static, or
    ``(nmat, nvol, NC, NC)`` on every bosonic slice."""
    view = BondSetView(BondTopology(DR, REV, {}))
    g = _green(seed, static)[None].copy()
    if static:
        return bubble.bond_bubble_static(g, None, CONVENTION_BETA, view,
                                         spatial_shape=CONVENTION_SHAPE)
    npair = CONVENTION_ND * CONVENTION_ND
    out = np.zeros((g.shape[1], NVOL, NC, NC), dtype=np.complex128)
    for (m, mp), blk in bubble._iter_bond_dynamic(g, None, CONVENTION_BETA, view,
                                                   spatial_shape=CONVENTION_SHAPE):
        out[:, :, m * npair:(m + 1) * npair, mp * npair:(mp + 1) * npair] = blk
    return out


@functools.lru_cache(maxsize=None)
def _sample(static):
    """A fixed random sample of code entries ``(l, q, r, c)``."""
    rng = np.random.default_rng(CONVENTION_SAMPLE_SEED)
    nmat = _nmat(static)
    full = (NVOL, NC, NC) if static else (nmat, NVOL, NC, NC)
    flat = rng.choice(int(np.prod(full)), size=CONVENTION_SAMPLE_SIZE, replace=False)
    idx = np.unravel_index(flat, full)
    if static:
        return (np.full(CONVENTION_SAMPLE_SIZE, nmat // 2, dtype=np.int64),) + tuple(idx)
    return tuple(idx)


@functools.lru_cache(maxsize=None)
def _scale(seed, static):
    return float(np.abs(_production(seed, static)).max())


def _production_entries(seed, static, entries):
    code = _production(seed, static)
    l, q, r, c = entries
    return code[q, r, c] if static else code[l, q, r, c]


# ---------------------------------------------------------------------------
# Candidate families
# ---------------------------------------------------------------------------

def _flags(static):
    """Static: ``(lambda_r, lambda_c, s_q, rowrev, colrev, p_r, p_c)`` (288).
    Dynamic: ``(lambda_r, lambda_c, s_q, s_nu, rowrev, colrev, t, p_r, p_c)``
    (1152)."""
    pm, bit, ph = (1, -1), (0, 1), (-1, 0, 1)
    if static:
        return list(itertools.product(pm, pm, pm, bit, bit, ph, ph))
    return list(itertools.product(pm, pm, pm, pm, bit, bit, bit, ph, ph))


def _full_flags(f, static):
    if static:
        lam_r, lam_c, s_q, rr, cr, p_r, p_c = f
        return (lam_r, lam_c, s_q, 1, rr, cr, 0, p_r, p_c)
    return f


def _static_part(view):
    lam_r, lam_c, s_q, _s_nu, rr, cr, _t, p_r, p_c = view
    return (lam_r, lam_c, s_q, rr, cr, p_r, p_c)


def _exchange_partner(f):
    lam_r, lam_c, s_q, s_nu, rr, cr, t, p_r, p_c = f
    return (lam_r, lam_c, -s_q, -s_nu, rr, cr, 1 - t, p_r, p_c)


def _prediction(f, seed, static, entries):
    return candidate_map(_label(seed, static), GEOM, CONVENTION_ND, _full_flags(f, static),
                         entries=entries)


def _residuals(flags, fixtures, *, static, entries):
    """``{flags: (residual per fixture)}``, relative to ``max|code|``."""
    res = {}
    for f in flags:
        r = []
        for seed in fixtures:
            got = _production_entries(seed, static, entries)
            r.append(float(np.abs(_prediction(f, seed, static, entries) - got).max()
                           / _scale(seed, static)))
        res[f] = tuple(r)
    return res


def _signature_classes(flags, fixtures, *, atol, entries):
    """Group candidates whose predictions agree on the sampled entries of
    every fixture (relative to ``max|code|``) to ``atol``."""
    sig = np.array([np.concatenate([
        _prediction(f, seed, False, entries) / _scale(seed, False)
        for seed in fixtures]) for f in flags])
    w = np.random.default_rng(0).normal(size=sig.shape[1])
    key = (sig @ w).real
    order = np.argsort(key)
    parent = list(range(len(flags)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    bound = atol * np.abs(w).sum()
    for a_pos, a in enumerate(order):
        for b in order[a_pos + 1:]:
            if key[b] - key[a] > bound:
                break
            if np.abs(sig[a] - sig[b]).max() <= atol:
                parent[root(a)] = root(b)
    groups = {}
    for i in range(len(flags)):
        groups.setdefault(root(i), []).append(flags[i])
    return [tuple(sorted(g, key=lambda f: f[6])) for g in groups.values()]


def _class_residuals(classes, fixtures, *, entries):
    out = {}
    for cls in classes:
        rep = [f for f in cls if f[6] == 0][0]
        out[cls] = _residuals([rep], fixtures, static=False, entries=entries)[rep]
    return out


def _full_residual(f, *, static):
    worst = 0.0
    for seed in FIXTURES:
        pred = candidate_map(_label(seed, static), GEOM, CONVENTION_ND, _full_flags(f, static))
        if static:
            pred = pred[_nmat(static) // 2]
        worst = max(worst, rel(pred, _production(seed, static)))
    return worst


def _table(res):
    lines = ["candidate residuals per fixture {} (sorted):".format(FIXTURES)]
    for key, r in sorted(res.items(), key=lambda kv: max(kv[1])):
        lines.append("  {}  {}".format(key, "  ".join("{:.3e}".format(x) for x in r)))
    return "\n".join(lines)


@functools.lru_cache(maxsize=None)
def _classes():
    return tuple(_signature_classes(_flags(static=False), FIXTURES, atol=SIGNATURE_ATOL,
                                    entries=_sample(False)))


class TestStaticBubbleReading(unittest.TestCase):
    """P0-a, static bubble: of the 288 candidate
    readings exactly one matches ``bond_bubble_static`` and it is the static
    part of ``SLOT_CONVENTION.label_view()``."""

    def test_static_winner(self):
        res = _residuals(_flags(static=True), FIXTURES, static=True, entries=_sample(True))
        win = [f for f, r in res.items() if max(r) <= 1e-12]
        ok = (len(res) == 288 and len(win) == 1
              and all(min(r) >= 1e-3 for f, r in res.items() if f != win[0])
              and _full_residual(win[0], static=True) <= 1e-12)
        if not ok:
            self.fail(_table(res))
        self.assertEqual(_static_part(SLOT_CONVENTION.label_view()), win[0])


class TestDynamicBubbleReading(unittest.TestCase):
    """P0-a, dynamic bubble: the 1152 candidates
    form 576 particle-exchange classes; exactly one class matches every
    bosonic slice of ``_iter_bond_dynamic`` and its ``t = 0`` member is
    ``SLOT_CONVENTION.label_view()``."""

    def test_classes_are_the_particle_exchange_pairs(self):
        classes = _classes()
        self.assertEqual(len(classes), 576)
        for cls in classes:
            self.assertEqual(len(cls), 2)
            a, b = cls
            self.assertEqual(b, _exchange_partner(a))

    def test_dynamic_winner_is_the_convention(self):
        res = _class_residuals(_classes(), FIXTURES, entries=_sample(False))
        win = [c for c, r in res.items() if max(r) <= 1e-12]
        ok = len(win) == 1 and all(min(r) >= 1e-3 for c, r in res.items() if c != win[0])
        if not ok:
            self.fail(_table(res))
        rep = [f for f in win[0] if f[6] == 0][0]
        self.assertLessEqual(_full_residual(rep, static=False), 1e-12)
        self.assertEqual(rep, SLOT_CONVENTION.label_view())


GRAMMAR_MESSAGE = ("slot convention '{}' does not match sr<+|->:sc<+|->:q<+|->:nu<+|->:"
                   "rowrev<0|1>:colrev<0|1>:phr<-|0|+>:phc<-|0|+>")


class TestSlotConvention(unittest.TestCase):
    def test_convention_views(self):
        self.assertTrue(SLOT_CONVENTION.is_operator_reading())
        self.assertEqual(SLOT_CONVENTION.label_view(), (1, -1, -1, 1, 1, 0, 0, -1, 0))
        self.assertEqual(SLOT_CONVENTION.to_string(), "sr-:sc-:q-:nu+:rowrev1:colrev0:phr-:phc0")
        self.assertEqual(len(SLOT_CONVENTION.to_string()), 40)
        self.assertEqual(SlotConvention.from_string(SLOT_CONVENTION.to_string()), SLOT_CONVENTION)

    def test_round_trip_every_convention(self):
        pm, bit, ph = (1, -1), (0, 1), (-1, 0, 1)
        n = 0
        for q, nu, rr, cr, pr, pc, sr, sc in itertools.product(pm, pm, bit, bit, ph, ph, pm, pm):
            conv = SlotConvention(q_sign=q, nu_sign=nu, row_rev=rr, col_rev=cr,
                                  phase_row=pr, phase_col=pc, slot_sign=(sr, sc))
            s = conv.to_string()
            self.assertEqual(s, s.encode("ascii").decode("ascii"))
            self.assertEqual(SlotConvention.from_string(s), conv)
            n += 1
        self.assertEqual(n, 576)

    def test_operator_reading_rejects_other_readings(self):
        base = dict(q_sign=-1, nu_sign=1, row_rev=1, col_rev=0, phase_row=-1, phase_col=0,
                    slot_sign=(-1, -1))
        for field, value in (("q_sign", 1), ("nu_sign", -1), ("row_rev", 0), ("col_rev", 1),
                             ("phase_row", 0), ("phase_col", -1)):
            self.assertFalse(SlotConvention(**dict(base, **{field: value})).is_operator_reading())
        self.assertTrue(SlotConvention(**dict(base, slot_sign=(1, -1), phase_row=1))
                        .is_operator_reading())

    def test_grammar_refusals(self):
        s = "sr-:sc-:q-:nu+:rowrev1:colrev0:phr-:phc0"
        for bad in (s[:-5], s.replace("rowrev1", "rowrev2"), "sc-:sr-" + s[7:], s + " ",
                    s + "\n", s.replace("-", "−"), s.upper(), 5, None):
            with self.assertRaisesRegex(ValueError,
                                        "^" + re.escape(GRAMMAR_MESSAGE.format(bad)) + "$"):
                SlotConvention.from_string(bad)

    def test_field_refusals(self):
        base = dict(q_sign=-1, nu_sign=1, row_rev=1, col_rev=0, phase_row=-1, phase_col=0,
                    slot_sign=(-1, -1))
        for field, value, phrase in (
                ("q_sign", 0, "SlotConvention: q_sign must be +1 or -1, got 0"),
                ("nu_sign", 2, "SlotConvention: nu_sign must be +1 or -1, got 2"),
                ("row_rev", -1, "SlotConvention: row_rev must be 0 or 1, got -1"),
                ("col_rev", 2, "SlotConvention: col_rev must be 0 or 1, got 2"),
                ("phase_row", 2, "SlotConvention: phase_row must be -1, 0 or +1, got 2"),
                ("phase_col", True, "SlotConvention: phase_col must be -1, 0 or +1, got True"),
                ("slot_sign", (1,), "SlotConvention: slot_sign must be a pair of +1/-1, got (1,)"),
                ("slot_sign", (1, 0), "SlotConvention: slot_sign must be a pair of +1/-1, "
                                      "got (1, 0)")):
            with self.assertRaisesRegex(ValueError, "^" + re.escape(phrase) + "$"):
                SlotConvention(**dict(base, **{field: value}))


class TestGeometry(unittest.TestCase):
    def test_phases_and_frozen(self):
        g = BondSlotGeometry.from_arrays(DR, REV, (5, 3, 1))
        self.assertEqual((g.nvol, g.B), (15, 7))
        self.assertEqual(g.phase.shape, (15, 7))
        self.assertEqual(g.q_frac.shape, (15, 3))
        self.assertEqual((g.delta_r.dtype, g.reverse.dtype, g.q_frac.dtype, g.phase.dtype),
                         (np.int64, np.int64, np.float64, np.complex128))
        for arr in (g.delta_r, g.reverse, g.q_frac, g.phase):
            self.assertFalse(arr.flags.writeable)
        self.assertIs(g.convention, SLOT_CONVENTION)
        self.assertEqual(g.spatial_shape, (5, 3, 1))
        np.testing.assert_allclose(g.phase[:, 1], np.exp(2j * np.pi * g.q_frac[:, 0]), atol=1e-15)
        # q_frac in C order: grid index / shape, consistent with the oracle mesh
        np.testing.assert_array_equal(2 * np.pi * g.q_frac, mesh_momenta(CONVENTION_SHAPE))
        np.testing.assert_allclose(g.phase, np.exp(1j * mesh_momenta(CONVENTION_SHAPE) @ DR.T),
                                   atol=1e-15)

    def test_geometry_copies_inputs(self):
        dr, rev = DR.copy(), REV.copy()
        g = BondSlotGeometry.from_arrays(dr, rev, (5, 3, 1))
        self.assertTrue(dr.flags.writeable and rev.flags.writeable)
        self.assertFalse(np.shares_memory(g.delta_r, dr) or np.shares_memory(g.reverse, rev))
        dr[1] = (7, 7, 7)
        rev[1] = 1
        np.testing.assert_array_equal(g.delta_r, DR)
        np.testing.assert_array_equal(g.reverse, REV)

    def test_geometry_mesh_mismatch(self):
        msg = ("BondSlotGeometry: channels 1 (delta_r=(1, 0, 0)) and 2 (delta_r=(-1, 0, 0)) "
               "coincide modulo the mesh (2, 3, 1)")
        with self.assertRaisesRegex(ValueError, "^" + re.escape(msg) + "$"):
            BondSlotGeometry.from_arrays(DR, REV, (2, 3, 1))     # +x and -x coincide on Nx = 2
        with self.assertRaisesRegex(
                ValueError, "^" + re.escape("BondSlotGeometry: reverse must be an involution "
                                            "with delta_r[reverse[m]] == -delta_r[m]; fails at "
                                            "channel 1 (reverse[1] = 1)") + "$"):
            BondSlotGeometry.from_arrays(DR, np.arange(7), (5, 3, 1))
        with self.assertRaisesRegex(
                ValueError, "^" + re.escape("BondSlotGeometry: reverse must have shape (7,) "
                                            "with entries in [0, 7), got [0, 2, 1]") + "$"):
            BondSlotGeometry.from_arrays(DR, [0, 2, 1], (5, 3, 1))
        with self.assertRaisesRegex(
                ValueError, "^" + re.escape("BondSlotGeometry: spatial_shape must be three "
                                            "positive integers, got (5, 3)") + "$"):
            BondSlotGeometry.from_arrays(DR, REV, (5, 3))
        with self.assertRaisesRegex(
                ValueError, "^" + re.escape("BondSlotGeometry: delta_r must be an integer "
                                            "(B, 3) array, got shape (7, 2)") + "$"):
            BondSlotGeometry.from_arrays(DR[:, :2], REV, (5, 3, 1))
        with self.assertRaisesRegex(
                ValueError, "^" + re.escape("BondSlotGeometry: convention must be a "
                                            "SlotConvention, got 'sr-'") + "$"):
            BondSlotGeometry.from_arrays(DR, REV, (5, 3, 1), convention="sr-")

    def test_non_numeric_reverse_refused(self):
        for rev, shown in ((["a"] * 7, "['a', 'a', 'a', 'a', 'a', 'a', 'a']"),
                           ([0, [2], 1, 4, 3, 6, 5], "[0, [2], 1, 4, 3, 6, 5]"),
                           ([0, None, 1, 4, 3, 6, 5], "[0, None, 1, 4, 3, 6, 5]")):
            with self.subTest(rev=rev):
                with self.assertRaisesRegex(
                        ValueError, "^" + re.escape("BondSlotGeometry: reverse must have shape "
                                                    "(7,) with entries in [0, 7), got " + shown)
                        + "$"):
                    BondSlotGeometry.from_arrays(DR, rev, (5, 3, 1))

    def test_slot_geometry_accepts_any_bond_set(self):
        for src in (SimpleNamespace(delta_r=[tuple(d) for d in DR], reverse=list(REV)),
                    BondTopology(DR, REV, {}),
                    BondSetView(BondTopology(DR, REV, {}))):
            g = slot_geometry(src, CONVENTION_SHAPE)
            np.testing.assert_array_equal(g.delta_r, DR)
            np.testing.assert_array_equal(g.reverse, REV)
            self.assertEqual(g.spatial_shape, CONVENTION_SHAPE)
            self.assertIs(g.convention, SLOT_CONVENTION)
        with self.assertRaisesRegex(ValueError, "both wrap to the same displacement "
                                                r"\(1, 0, 0\) modulo the mesh"):
            slot_geometry(SimpleNamespace(delta_r=DR, reverse=REV), (2, 3, 1))


# ---------------------------------------------------------------------------
# The bond pairing kernel against its documented formula
# ---------------------------------------------------------------------------

def _kernel_brute(V, phi, *, sign, nu_flip):
    """``out_AD(k, n) = -(1/N) sum_{alpha, beta} e^{sign ik.R_beta} sum_{k', n'}
    V_{alpha beta; ABCD}(k - k', l) e^{sign ik'.R_alpha} phi_BC(k', n')`` with
    ``l = (n - n' + nmat/2) mod nmat`` (reversed to ``(nmat - l) mod nmat``
    with ``nu_flip``)."""
    nfreq = V.shape[0]
    norb = phi.shape[0]
    V8 = V.reshape(nfreq, NVOL, B, norb, norb, B, norb, norb)
    kph = np.exp(sign * 1j * mesh_momenta(CONVENTION_SHAPE) @ DR.T)          # (k, m)
    weighted = np.einsum("pa,bcpn->apnbc", kph, phi)                   # (alpha, k', n', B, C)
    kminus = mesh_add(CONVENTION_SHAPE)[:, mesh_negate(CONVENTION_SHAPE)]    # [k, k'] -> k - k'
    out = np.zeros_like(phi)
    for k in range(NVOL):
        for n in range(nfreq):
            l = (n - np.arange(nfreq) + nfreq // 2) % nfreq
            if nu_flip:
                l = (nfreq - l) % nfreq
            blk = V8[l[:, None], kminus[k][None, :]]                   # (n', k', a, A, B, b, C, D)
            acc = np.einsum("npaABbCD,apnBC->bAD", blk, weighted)
            out[:, :, k, n] = -(1.0 / NVOL) * np.einsum("b,bAD->AD", kph[k], acc)
    return out


def _kernel_vs_formula(nfreq, seed):
    """Production ``BondPairKernel.matvec`` (``G2 = delta``, so ``F = phi``)
    on a random vertex, its documented formula, the competing reading
    (``e^{+ik.R}`` for a one-point axis, the reversed bosonic index
    otherwise) and the ``e^{+ik.R}`` reading."""
    rng = np.random.default_rng(seed)
    norb = CONVENTION_ND
    npair = norb * norb
    ND = B * npair
    Vq = rng.normal(size=(NVOL, ND, ND)) + 1j * rng.normal(size=(NVOL, ND, ND))
    if nfreq == 1:
        V = Vq[None].copy()
    else:
        h = rng.normal(size=nfreq) + 1j * rng.normal(size=nfreq)
        V = Vq[None] * h[:, None, None, None]
    src = ArrayBlockSource({"V": V.copy()}, npair)
    vert = PairVertexUniform(src, "V", B, npair, nfreq, NVOL)
    G2 = np.zeros((norb,) * 4 + CONVENTION_SHAPE + (nfreq,), dtype=np.complex128)
    for i in range(norb):
        for j in range(norb):
            G2[i, i, j, j] = 1.0
    view = SimpleNamespace(n_channels=B, delta_r=tuple(tuple(int(x) for x in d) for d in DR))
    kern = BondPairKernel(vert, G2, view, xp=np, spatial_shape=CONVENTION_SHAPE, norb=norb,
                          beta=CONVENTION_BETA, nfreq=nfreq, residency="host", host_cap=np.inf)
    gap = (norb, norb, NVOL, nfreq)
    phi = rng.normal(size=gap) + 1j * rng.normal(size=gap)
    got = np.asarray(kern.matvec(phi.reshape((norb, norb) + CONVENTION_SHAPE + (nfreq,)).ravel()))
    got = got.reshape(norb, norb, NVOL, nfreq)
    want = _kernel_brute(V, phi, sign=-1, nu_flip=False)
    plus_r = _kernel_brute(V, phi, sign=+1, nu_flip=False)
    if nfreq == 1:
        return got, want, plus_r, plus_r
    return got, want, _kernel_brute(V, phi, sign=-1, nu_flip=True), plus_r


class TestPairKernelFormula(unittest.TestCase):
    """P0-b: ``BondPairKernel`` applies
    ``out_AD(k) = -(1/N) sum e^{-ik.R_beta} V(k - k') e^{-ik'.R_alpha} F(k')``
    with ``l = (n - n' + nmat/2) mod nmat`` and unit weight."""

    def test_static_one_point_axis(self):
        got, want, flipped, _ = _kernel_vs_formula(nfreq=1, seed=2705)
        self.assertLessEqual(rel(got, want), 1e-12)
        self.assertGreaterEqual(rel(flipped, want), 1e-3)

    def test_factorized_dynamic_vertex(self):
        got, want, nu_flipped, plus_r = _kernel_vs_formula(nfreq=8, seed=2706)
        self.assertLessEqual(rel(got, want), 1e-12)
        self.assertGreaterEqual(rel(nu_flipped, want), 1e-3)
        self.assertGreaterEqual(rel(plus_r, want), 1e-3)


# ---------------------------------------------------------------------------
# Production Gamma_B against the exact first-order vertex part
# ---------------------------------------------------------------------------

#: the number of random sketch pairs at the head of every ``sketch_set``
N_RANDOM = 8
FILING_TOL = 1e-6
#: targeted-sketch floor (relative to the largest random sketch) of a run with every type
TARGETED_FLOOR = 1e-2
#: the floor of a single-type run: one type leaves some sector or bond-orientation sketches
#: weakly populated; at 1e-3 the 1e-6 deviation still compares three significant digits
TARGETED_FLOOR_PER_TYPE = 1e-3


def _restrict(fixture, types):
    """The fixture with only the interaction tables of ``types`` ("all"
    keeps the fixture itself, so its exact-diagonalization cache is
    shared)."""
    if types == "all":
        return fixture
    return fixture._replace(tables=tuple((t, tbl) for t, tbl in fixture.tables if t in types))


def _p0c(fixture, *, types, qs, sketches, richardson=None):
    """``(dev, targeted, table)``: ``dev = max |ED - pred| / max |ED|`` over
    the sketches and ``qs``, where ``ED`` is the first-order vertex part of
    the exact bond correlators and ``pred = u^T [-chibar(q) Gamma(-q)
    chibar(q)] v`` with the production vertex in the label form;
    ``targeted`` is each targeted sketch's ``max_q |ED|`` over the largest
    random one; ``table`` the per-sketch residuals (shown on failure)."""
    name = fixture_name(fixture)
    sub = _restrict(fixture, types)
    qs = tuple(int(q) for q in qs)
    if richardson is None:
        richardson = NEEDS_RICHARDSON[name]
    ref = vertex_part_sketches(sub, sketches=sketches, qs=qs, richardson=richardson,
                               x=FIXTURE_X[name])
    pred = production_first_order(sub, sub, sketch_set(sub, sketches), qs)
    scale = np.abs(ref).max()
    dev = float(np.abs(pred - ref).max() / scale)
    rnd = np.abs(ref[:, :N_RANDOM]).max()
    targeted = [float(np.abs(ref[:, i]).max() / rnd) for i in range(N_RANDOM, ref.shape[1])]
    lines = ["fixture {} types {} qs {}: dev {:.3e}, targeted {}".format(
        name, types, qs, dev, ["{:.3e}".format(t) for t in targeted])]
    for s in range(ref.shape[1]):
        lines.append("  sketch {:2d}  |ED| {:.3e}  residual/scale {:.3e}".format(
            s, np.abs(ref[:, s]).max(), np.abs(pred[:, s] - ref[:, s]).max() / scale))
    return dev, targeted, "\n".join(lines)


def _pred_max(fixture, types, *, qs):
    """``max_sketch |pred|`` of the production vertex of ``types`` alone on
    the fast sketches of the full fixture (no exact diagonalization)."""
    sub = _restrict(fixture, tuple(types))
    return float(np.abs(production_first_order(sub, fixture, sketch_set(fixture, "fast"), tuple(qs))).max())


class TestFilingAgainstED(unittest.TestCase):
    """P0-c: the production ``Gamma_B``
    (``so_bond.so_bond_vertex``), mapped to the label form by the test-side
    ``gamma_code_to_label``, against the first-order vertex part of the
    exact bond correlators (full first order minus the Hartree-Fock-shifted
    first order) through ``u^T [-chibar(q) Gamma(-q) chibar(q)] v`` with the
    analytic static Lindhard ``chibar``.

    This validates the filing of every interaction type into the vertex.
    The code-slot rule itself is not decided here: the label map shares the
    slot rule with production; the dressing-algebra ring check and the
    pairing first-order oracle decide it."""

    def _check(self, fixture, types, *, qs, sketches, floor=TARGETED_FLOOR, richardson=None):
        dev, targeted, table = _p0c(fixture, types=types, qs=qs, sketches=sketches,
                                    richardson=richardson)
        if not (dev <= FILING_TOL and all(t >= floor for t in targeted)):
            self.fail(table)
        return dev, targeted

    def test_fixture_i_combined(self):
        _, targeted = self._check(FIXTURE_I, "all", qs=(1, 2), sketches="fast")
        self.assertEqual(len(targeted), 3)
        whole = _pred_max(FIXTURE_I, TYPES, qs=(1, 2))
        for t in TYPES:
            self.assertGreaterEqual(_pred_max(FIXTURE_I, (t,), qs=(1, 2)),
                                    TARGETED_FLOOR * whole, t)

    @heavy
    def test_fixture_i_per_type(self):
        for t in TYPES:
            with self.subTest(type=t):
                _, targeted = self._check(FIXTURE_I, (t,), qs=range(FIXTURE_I.nvol),
                                          sketches="full", floor=TARGETED_FLOOR_PER_TYPE)
                self.assertEqual(len(targeted), 6)

    @heavy
    def test_fixture_ii_combined(self):
        _, targeted = self._check(FIXTURE_II, "all", qs=(1, 2), sketches="full")
        self.assertEqual(len(targeted), 6)


if __name__ == "__main__":
    unittest.main()
