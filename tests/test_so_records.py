"""The spin-orbital quartic-record compiler (``hwave.solver.so_records``).

The operator built from the compiled records is compared with the operator
of the independent table-driven builder (``tests.so_bond_oracles``), applied
to random vectors of the N = 1 and N = 2 particle sectors of a 5-cell,
two-orbital chain (20 modes); full many-body operators are never built.
"""

import logging
import unittest

import numpy as np
from requests.structures import CaseInsensitiveDict

from hwave.solver.declarations import symmetrise_dense
from hwave.solver.so_records import (INSIGNIFICANT_RTOL, RECORD_DTYPE, REMAINDER_DTYPE, TYPES,
                                     compile_interaction, symmetrised_table)
from hwave.solver.declarations import KNOWN_TYPES
from tests.so_bond_oracles import (COMPILER_TABLES, INTRA, T1, T2, closed_row, rel, apply_terms,
                                   build_terms, fock_states, translate_monomials)

L = 5
NORB = 2
NMODE = 2 * NORB * L


def _vecs(N, seed=3101, k=4):
    """``(states, vecs)``: ``k`` random complex vectors in the N-particle sector."""
    states = fock_states(NMODE, N)
    rng = np.random.default_rng(seed + N)
    vecs = rng.standard_normal((len(states), k)) + 1j * rng.standard_normal((len(states), k))
    return states, vecs


def _record_terms(records, L, norb_phys, order=(0, 1, 2, 3), sign=1.0):
    out = []
    for r in records:
        modes = translate_monomials(r["modes"], r["cells"], L, norb_phys)
        for row in modes:
            out.append(tuple(int(row[i]) for i in order) + (sign * complex(r["v"]),))
    return out


def _apply_records(records, L, norb_phys, N):
    return apply_terms(_record_terms(records, L, norb_phys), *_vecs(N))


def _forms(c, L, norb_phys, N):
    """``(direct, crossed, remainders)`` applied to the sector's random vectors:
    ``sum v c+1 c2 c+3 c4``, ``sum -v c+1 c4 c+3 c2`` and the one-body
    remainders ``sum v c+1 c2``."""
    states, vecs = _vecs(N)
    direct = apply_terms(_record_terms(c.records, L, norb_phys), states, vecs)
    crossed = apply_terms(_record_terms(c.records, L, norb_phys, (0, 3, 2, 1), -1.0), states, vecs)
    rem = []
    for r in c.remainders:
        for p, q in translate_monomials(r["modes"], r["cells"], L, norb_phys):
            rem.append((int(p), int(q), complex(r["v"])))
    return direct, crossed, apply_terms(rem, states, vecs)


def _permuted(tables, seed):
    """The same tables with the type order and every table's entry order
    shuffled."""
    rng = np.random.default_rng(seed)
    names = list(tables)
    out = {}
    for i in rng.permutation(len(names)):
        items = list(tables[names[i]].items())
        out[names[i]] = dict(items[j] for j in rng.permutation(len(items)))
    return out


def _with_zero_shell_and_tiny():
    """A declared-zero CoulombInter shell at R = +-2 and an Exchange bond of
    1e-300 (no 1/2 in the Exchange operator, so its records carry 1e-300)."""
    return {"CoulombInter": {**T1, **closed_row(2, 0, 1, 0.0)},
            "Exchange": closed_row(1, 0, 1, 1e-300)}


_TINY = 2.0 ** -60          # exact power of two: scaling commutes with every rounding


def _tiny_outer_shell(scale=_TINY):
    """CoulombInter and Hund on R = +-1 with their R = +-2 rows scaled by
    ``scale``."""
    return {"CoulombInter": {**T1, **closed_row(2, 1, 0, 0.4 * scale)},
            "Hund": {**closed_row(1, 0, 1, 0.53), **closed_row(-2, 0, 1, 0.25 * scale)}}


def _sparse_mean_equals_symmetrise_dense(table, grid, hermitian=False, itype="Hund"):
    """The sparse partner-key mean equals ``symmetrise_dense`` on a grid on
    which no declared displacement wraps onto another."""
    dense = np.zeros(tuple(grid) + (NORB, NORB), dtype=np.complex128)
    for (R, (a, b)), v in table.items():
        dense[tuple(np.mod(R, grid)) + (a, b)] += v
    want = symmetrise_dense(dense, hermitian=hermitian)
    got = np.zeros_like(want)
    for (R, (a, b)), v in symmetrised_table(itype, table).items():
        got[tuple(np.mod(R, grid)) + (a, b)] += v
    return bool(np.array_equal(got, want))


_SAME_OPERATOR = ("{t} coefficient at R={R} orbitals ({a}, {b}) has an imaginary part ({im:.6g}). "
                  "In H-wave's declaration convention the entry and its partner (−R,b,a) "
                  "multiply the same operator, so the imaginary part has no operator content; "
                  "only PairHop carries a physical complex phase (complex couplings beyond "
                  "PairHop: issue #302). Remove the imaginary part.")
_EXCHANGE = ("Exchange coefficient at R={R} orbitals ({a}, {b}) has an imaginary part ({im:.6g}). "
             "The partner entry (−R,b,a) multiplies the Hermitian conjugate, so an imaginary "
             "part would be a Dzyaloshinskii-Moriya-type term; H-wave's declaration convention "
             "folds it, and that extension is tracked in issue #302. Remove the imaginary part.")


def _expected_message(t, R, a, b, im):
    return (_EXCHANGE if t == "Exchange" else _SAME_OPERATOR).format(t=t, R=R, a=a, b=b, im=im)


class TestOperatorIdentity(unittest.TestCase):
    """Gate A1-C: the records reproduce the documented
    operator, and the crossed re-pairing plus the one-body remainders
    reproduce the direct form."""

    def test_records_equal_the_builder(self):
        # the documented monomials are not normal ordered, but none of them
        # repeats a mode in one cell, so the N = 1 sector is annihilated: it
        # is compared on the scale of the N = 2 one
        c = compile_interaction(COMPILER_TABLES, NORB)
        want_terms = build_terms(COMPILER_TABLES, L, NORB)
        got = {N: _apply_records(c.records, L, NORB, N) for N in (1, 2)}
        want = {N: apply_terms(want_terms, *_vecs(N)) for N in (1, 2)}
        self.assertGreater(np.abs(want[2]).max(), 1e-2)
        self.assertLessEqual(rel(got[2], want[2]), 1e-13)
        self.assertLessEqual(np.abs(got[1] - want[1]).max(), 1e-13 * np.abs(want[2]).max())

    def test_repairing_identity_with_remainders(self):
        c = compile_interaction(COMPILER_TABLES, NORB)
        forms = {N: _forms(c, L, NORB, N) for N in (1, 2)}
        scale = np.abs(forms[2][0]).max()
        self.assertGreater(scale, 1e-2)
        for N in (1, 2):
            direct, crossed, rem = forms[N]
            self.assertGreater(np.abs(rem).max(), 1e-2 * scale)              # non-vacuous
            self.assertLessEqual(np.abs(direct - crossed - rem).max(), 1e-13 * scale)


class TestRecords(unittest.TestCase):
    """Gate A1-C: canonical records, declared topology,
    exact zeros, type keys and the insignificant-shell diagnostic."""

    def test_layout(self):
        c = compile_interaction(COMPILER_TABLES, NORB)
        self.assertEqual(c.norb_phys, NORB)
        self.assertEqual(c.records.dtype, RECORD_DTYPE)
        self.assertEqual(c.records_all.dtype, RECORD_DTYPE)
        self.assertEqual(c.remainders.dtype, REMAINDER_DTYPE)
        self.assertEqual(c.declared_topology.dtype, np.int64)
        self.assertTrue(np.all(c.records["cells"][:, 0] == 0))
        for f in ("records", "records_all", "declared_topology", "remainders"):
            self.assertFalse(getattr(c, f).flags.writeable, f)
        topo = [tuple(r) for r in c.declared_topology]
        self.assertEqual(topo, sorted(set(topo)))
        self.assertIn((0, 0, 0), topo)
        self.assertEqual(set(topo), {tuple(-x for x in d) for d in topo})
        self.assertEqual(set(topo), {(r, 0, 0) for r in (-2, -1, 0, 1, 2)})
        keys = [(r["modes"].tolist(), r["cells"].tolist()) for r in c.remainders]
        self.assertEqual(keys, sorted(keys))
        self.assertEqual(len(keys), len({repr(k) for k in keys}))
        rk = [(str(r["type"]), r["modes"].tolist(), r["cells"].tolist()) for r in c.records_all]
        self.assertEqual(len(rk), len({repr(k) for k in rk}))

    def test_permutation_byte_identical(self):
        ref = compile_interaction(COMPILER_TABLES, NORB)
        for seed in range(88, 93):
            c = compile_interaction(_permuted(COMPILER_TABLES, seed), NORB)
            for f in ("records", "records_all", "declared_topology", "remainders"):
                self.assertEqual(getattr(c, f).tobytes(), getattr(ref, f).tobytes(), f)

    def test_declared_zero_and_tiny(self):
        c = compile_interaction(_with_zero_shell_and_tiny(), NORB)
        self.assertFalse(np.any(c.records["v"] == 0))
        self.assertTrue(np.any(c.records_all["v"] == 0))
        self.assertEqual(len(c.records), int(np.count_nonzero(c.records_all["v"])))
        self.assertIn((2, 0, 0), map(tuple, c.declared_topology))
        self.assertIn((-2, 0, 0), map(tuple, c.declared_topology))
        self.assertTrue(np.any(np.abs(c.records["v"]) == 1e-300))

    def test_non_finite_refused_before_sorting(self):
        for bad, text in ((np.nan, "nan"), (np.inf, "inf")):
            with self.assertRaisesRegex(
                    ValueError, r"^Hund coefficient at R=\(-1, 0, 0\) orbitals \(2, 1\) "
                                r"is not finite \({}\)$".format(text)):
                compile_interaction({"Hund": closed_row(1, 0, 1, bad)}, NORB)

    def test_symmetrisation_matches_dense_and_does_not_wrap(self):
        t = COMPILER_TABLES["Hund"]
        self.assertTrue(_sparse_mean_equals_symmetrise_dense(t, grid=(9, 1, 1)))
        # an asymmetric (unclosed) table exercises the mean itself
        rng = np.random.default_rng(5)
        raw = {((r, 0, 0), (a, b)): complex(*rng.standard_normal(2))
               for r in (-2, -1, 0, 1, 3) for a in range(NORB) for b in range(NORB)}
        self.assertTrue(_sparse_mean_equals_symmetrise_dense(raw, grid=(9, 1, 1)))
        self.assertTrue(_sparse_mean_equals_symmetrise_dense(raw, (9, 1, 1), True, "PairHop"))
        c = compile_interaction({"CoulombInter": {**closed_row(4, 0, 1, 0.5),
                                                  **closed_row(-1, 0, 1, 0.2)}}, NORB)
        self.assertTrue({(4, 0, 0), (-1, 0, 0)} <= {tuple(r) for r in c.declared_topology})
        self.assertEqual(sorted({abs(v) for v in c.records["v"]}), [0.1, 0.25])

    def test_mixed_case_duplicate_keys_refused(self):
        with self.assertRaises(ValueError) as cm:
            compile_interaction({"Hund": T2, "HUND": T2}, NORB)
        self.assertEqual(str(cm.exception),
                         "interaction type given twice with different spelling: 'HUND' and 'Hund'")

    def test_type_keys_coulomb_and_closure(self):
        ref = compile_interaction({"CoulombInter": T1, "Hund": T2}, NORB).records.tobytes()
        self.assertEqual(compile_interaction({"coulombinter": T1, "HUND": T2}, NORB)
                         .records.tobytes(), ref)
        self.assertEqual(compile_interaction(CaseInsensitiveDict({"coulombINTER": T1, "hund": T2}),
                                             NORB).records.tobytes(), ref)
        self.assertEqual(compile_interaction({"Coulomb": {**T1, **INTRA}}, NORB).records.tobytes(),
                         compile_interaction({"CoulombInter": T1, "CoulombIntra": INTRA},
                                             NORB).records.tobytes())
        ignored = {"Geometry": object(), "transfer": object(), "EXTERN": object()}
        self.assertEqual(compile_interaction({**ignored, "CoulombInter": T1, "Hund": T2},
                                             NORB).records.tobytes(), ref)
        with self.assertRaisesRegex(ValueError, "^Coulomb cannot be specified together with "
                                                "CoulombIntra or CoulombInter$"):
            compile_interaction({"Coulomb": T1, "CoulombInter": T1}, NORB)
        with self.assertRaisesRegex(ValueError, r"^unknown interaction type 'Hubbard'"):
            compile_interaction({"Hubbard": T1}, NORB)
        with self.assertRaisesRegex(ValueError, r"^Hund declares R=\(1, 0, 0\) orbitals \(1, 2\) "
                                                r"= 0\.3 with NO partner entry"):
            compile_interaction({"Hund": {((1, 0, 0), (0, 1)): 0.3}}, NORB)
        with self.assertRaisesRegex(ValueError, r"^CoulombInter: the on-site row of orbital 1 "
                                                r"with itself"):
            compile_interaction({"CoulombInter": {((0, 0, 0), (0, 0)): 0.3}}, NORB)
        with self.assertRaisesRegex(ValueError, r"^Hund coefficient at R=\(-1, 0, 0\) orbitals "
                                                r"\(3, 1\) names an orbital outside 1\.\.2$"):
            compile_interaction({"Hund": closed_row(1, 0, 2, 0.3)}, NORB)

    def test_types_are_the_declaration_types(self):
        self.assertEqual(frozenset(TYPES), KNOWN_TYPES)
        self.assertEqual(len(TYPES), len(KNOWN_TYPES))

    def test_non_integral_keys_refused(self):
        # the first offending key in the table's own order is named
        cases = (
            ({((1.7, 0, 0), (0, 1)): 0.3, ((-1.7, 0, 0), (1, 0)): 0.3}, "((1.7, 0, 0), (0, 1))"),
            ({((1, 0, 0), (0, 1.5)): 0.3, ((-1, 0, 0), (1, 0)): 0.3}, "((1, 0, 0), (0, 1.5))"),
            ({((1, 0, 0), (0, True)): 0.3, ((-1, 0, 0), (1, 0)): 0.3}, "((1, 0, 0), (0, True))"),
            ({((1, 0, "0"), (0, 1)): 0.3, ((-1, 0, 0), (1, 0)): 0.3}, "((1, 0, '0'), (0, 1))"),
        )
        for table, shown in cases:
            want = ("CoulombInter: entry key {} has a component that is not an exact integer; "
                    "displacements and orbitals must be integers".format(shown))
            with self.subTest(shown=shown):
                with self.assertRaises(ValueError) as cm:
                    compile_interaction({"CoulombInter": table}, NORB)
                self.assertEqual(str(cm.exception), want)
        # integral floats are exact integers and are read as such
        flt = {((1.0, 0.0, 0.0), (0.0, 1.0)): 0.37, ((-1.0, 0.0, 0.0), (1.0, 0.0)): 0.37}
        self.assertEqual(compile_interaction({"CoulombInter": flt}, NORB).records.tobytes(),
                         compile_interaction({"CoulombInter": closed_row(1, 0, 1, 0.37)},
                                             NORB).records.tobytes())
        with self.assertRaisesRegex(ValueError, r"^PairHop: entry key \(\(0\.5, 0, 0\), \(0, 1\)\) "
                                                r"has a component that is not an exact integer"):
            symmetrised_table("PairHop", {((0.5, 0, 0), (0, 1)): 0.3})

    def test_malformed_keys_refused(self):
        for key in (((0, 0, 0), (0, 1), 5), "R01", ((0, 0), (0, 1)), ((0, 0, 0), (0, 1, 2)), 7):
            want = ("Hund: an entry key must be ((Rx, Ry, Rz), (a, b)), got {!r}".format(key))
            with self.subTest(key=key):
                with self.assertRaises(ValueError) as cm:
                    compile_interaction({"Hund": {key: 0.3}}, NORB)
                self.assertEqual(str(cm.exception), want)

    def test_insignificant_shell_info_only(self):
        with self.assertLogs("hwave.solver.so_records", "INFO") as log:
            c = compile_interaction(_tiny_outer_shell(), NORB)
        lines = [r for r in log.output if "bond shell" in r]
        self.assertEqual(lines, ["INFO:hwave.solver.so_records:so_records: bond shell(s) "
                                 "[((-2, 0, 0), (2, 0, 0))] carry only couplings with "
                                 "|v| <= 1e-14 * max|v|"])
        # the INFO line is the only effect: the tiny shell's records are the
        # unit-scale ones times the (power-of-two) scale, bit for bit
        ref = compile_interaction(_tiny_outer_shell(1.0), NORB)
        outer = np.any(np.abs(ref.records["cells"][:, :, 0]) == 2, axis=1)
        self.assertTrue(outer.any() and not outer.all())
        want = ref.records.copy()
        want["v"][outer] *= _TINY
        self.assertEqual(c.records.tobytes(), want.tobytes())
        self.assertEqual(c.records_all.tobytes(), want.tobytes())
        self.assertEqual(c.declared_topology.tobytes(), ref.declared_topology.tobytes())
        self.assertLessEqual(np.abs(c.records["v"][outer]).max(),
                             INSIGNIFICANT_RTOL * np.abs(c.records["v"]).max())

    def test_all_zero_lists_every_offsite_shell(self):
        zero = {"Hund": {**closed_row(1, 0, 1, 0.0), **closed_row(2, 0, 1, 0.0),
                         **closed_row(0, 0, 1, 0.0)}}
        with self.assertLogs("hwave.solver.so_records", "INFO") as log:
            c = compile_interaction(zero, NORB)
        self.assertEqual(len(c.records), 0)
        self.assertGreater(len(c.remainders), 0)                 # exact zeros are kept
        self.assertTrue(np.all(c.remainders["v"] == 0))
        self.assertEqual(log.output, ["INFO:hwave.solver.so_records:so_records: bond shell(s) "
                                      "[((-2, 0, 0), (2, 0, 0)), ((-1, 0, 0), (1, 0, 0))] carry "
                                      "only couplings with |v| <= 1e-14 * max|v|"])

    def test_significant_tables_log_nothing(self):
        logger = logging.getLogger("hwave.solver.so_records")
        seen = []
        handler = logging.Handler(logging.DEBUG)
        handler.emit = seen.append
        logger.addHandler(handler)
        old = logger.level
        logger.setLevel(logging.DEBUG)
        try:
            compile_interaction(COMPILER_TABLES, NORB)
        finally:
            logger.removeHandler(handler)
            logger.setLevel(old)
        self.assertEqual(seen, [])


class TestComplexContract(unittest.TestCase):
    """Gate A1-C: only PairHop may carry an imaginary part."""

    def test_messages_onsite_every_type(self):
        for t in ("CoulombInter", "Hund", "Ising", "PairLift", "Exchange"):
            with self.assertRaises(ValueError) as cm:
                compile_interaction({t: closed_row(0, 0, 1, 0.3 + 1e-6j)}, NORB)
            self.assertEqual(str(cm.exception), _expected_message(t, (0, 0, 0), 1, 2, 1e-6))

    def test_offsite_reports_first_sorted_entry(self):
        with self.assertRaises(ValueError) as cm:
            compile_interaction({"Hund": closed_row(1, 0, 1, 0.3 + 1e-6j)}, NORB)
        self.assertEqual(str(cm.exception), _expected_message("Hund", (-1, 0, 0), 2, 1, -1e-6))

    def test_threshold(self):
        tol = 1e-10 + 1e-8 * abs(0.3 + 3.1e-9j)                       # 3.1e-9
        compile_interaction({"Hund": closed_row(0, 0, 1, 0.3 + 0.99 * tol * 1j)}, NORB)
        with self.assertRaises(ValueError) as cm:
            compile_interaction({"Hund": closed_row(0, 0, 1, 0.3 + 1.01 * tol * 1j)}, NORB)
        self.assertEqual(str(cm.exception),
                         _expected_message("Hund", (0, 0, 0), 1, 2, 1.01 * tol))
        # above the absolute 1e-10: only the relative term admits it
        compile_interaction({"Hund": closed_row(0, 0, 1, 0.3 + 2e-9j)}, NORB)

    def test_pairhop_accepted(self):
        c = compile_interaction({"PairHop": closed_row(1, 0, 1, 0.4 + 0.3j)}, NORB)
        self.assertGreater(np.abs(c.records["v"].imag).max(), 0.1)


if __name__ == "__main__":
    unittest.main()
