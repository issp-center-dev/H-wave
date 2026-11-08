"""The shared spin-orbital entry rules (``hwave.solver.spin_orbital``).

A public spin-orbital entry builds one :class:`SOEntryContext`: the
``orbital_site`` key is read case-insensitively and never parsed (one
WARNING when present), and the interaction is compiled once, so the
compiler's refusals (complex coefficients, Hermitian closure) apply at the
entry. A nested solver reuses the stored context.
"""

import inspect
import logging
import os
import unittest
from unittest import mock

import numpy as np
from requests.structures import CaseInsensitiveDict

from hwave.solver import so_records, spin_orbital
from hwave.solver.spin_orbital import (ORBITAL_SITE_WARNING, SO_ENTRY_CONTEXT_KEY,
                                       SOEntryContext, build_so_entry_context,
                                       context_for_nested, interaction_tables_digest,
                                       offsite_types)
from tests.so_bond_oracles import (EXPECTED_FLEX_PRODUCER, EXPECTED_RPA_PRODUCER,
                                   UNVERSIONED_REFUSAL_TEXT, V1_MEMBERS,
                                   V1_REFUSAL_TEXT, V2_ONLY_MEMBERS, closed_row)
from tests.so_run_compare import eigenvalues_close, numeric_mismatches

LOG = logging.getLogger("qlms.test_so_entry_rules")
NORB = 2


def W(text):
    return "WARNING:{}:{}".format(LOG.name, text)


def _same_records(test, got, want):
    test.assertEqual(got.records.shape, want.records.shape)
    for field in got.records.dtype.names:
        test.assertTrue(np.array_equal(got.records[field], want.records[field]), field)
    test.assertTrue(np.array_equal(got.declared_topology, want.declared_topology))


#: on-site intra-orbital U, an R = 0 cross-orbital Hund row and one
#: nearest-neighbour CoulombInter row
TABLES = {
    "CoulombIntra": {((0, 0, 0), (0, 0)): 4.0, ((0, 0, 0), (1, 1)): 4.0},
    "Hund": closed_row((0, 0, 0), 0, 1, 0.5),
    "CoulombInter": closed_row((1, 0, 0), 0, 0, 0.7),
}


class TestOrbitalSite(unittest.TestCase):

    def test_warning_once_and_value_not_parsed(self):
        for param in ({"Orbital_Site": "garbage"}, {"orbital_site": [0]},
                      CaseInsensitiveDict({"ORBITAL_SITE": [0, 1]})):
            with self.subTest(param=dict(param)):
                with self.assertLogs(LOG, "WARNING") as cm:
                    ctx = build_so_entry_context(param, TABLES, NORB, LOG)
                self.assertEqual(cm.output, [W(ORBITAL_SITE_WARNING)])
                self.assertTrue(ctx.orbital_site_warned)
                self.assertIsInstance(ctx.compiled, so_records.CompiledInteraction)

    def test_warning_text(self):
        self.assertEqual(
            ORBITAL_SITE_WARNING,
            "orbital_site is ignored in spin-orbital mode; all R=0 rows are treated "
            "exactly (see migration note)")

    def test_param_dict_records_only_the_orbital_site_read(self):
        # the hwave entry reports the [mode.param] keys a run never read
        # (ParamDict.accessed); the lookup must not mark the other keys read
        from hwave.solver.base import ParamDict
        for raw in ({"T": 0.1, "Orbital_Site": [0, 1]}, {"T": 0.1}):
            with self.subTest(raw=raw):
                param = ParamDict(raw)
                build_so_entry_context(param, TABLES, NORB,
                                       logging.getLogger("qlms.test_so_entry_rules.quiet"))
                self.assertEqual(param.accessed, {"orbital_site"})

    def test_absent_key(self):
        for param in ({}, None, {"orbital": [0, 1]}):
            with self.subTest(param=param):
                with self.assertNoLogs(LOG, "WARNING"):
                    ctx = build_so_entry_context(param, TABLES, NORB, LOG)
                self.assertFalse(ctx.orbital_site_warned)


class TestCompile(unittest.TestCase):

    def test_compiled_once_and_equal_to_the_compiler(self):
        real = so_records.compile_interaction
        with mock.patch.object(spin_orbital.so_records, "compile_interaction",
                               side_effect=real) as spy:
            ctx = build_so_entry_context({}, TABLES, NORB, LOG)
        self.assertEqual(spy.call_count, 1)
        _same_records(self, ctx.compiled, real(TABLES, NORB))

    def test_only_interaction_types_reach_the_compiler(self):
        # an Initial table (and Geometry/Transfer) is not an interaction type
        tables = dict(TABLES, Initial={((0, 0, 0), (0, 0)): 1.0},
                      Transfer={((0, 0, 0), (0, 0)): 1.0})
        ctx = build_so_entry_context({}, tables, NORB, LOG)
        _same_records(self, ctx.compiled, so_records.compile_interaction(TABLES, NORB))

    def test_misspelled_type_is_dropped_not_refused(self):
        # the readers own the interaction-key validation; this entry passes
        # only the interaction types (and Coulomb) to the compiler
        tables = dict(TABLES, Hunds=closed_row((1, 0, 0), 0, 1, 0.5))
        ctx = build_so_entry_context({}, tables, NORB, LOG)
        _same_records(self, ctx.compiled, so_records.compile_interaction(TABLES, NORB))

    def test_mixed_case_type_keys(self):
        tables = CaseInsensitiveDict({"coulombintra": TABLES["CoulombIntra"],
                                      "HUND": TABLES["Hund"],
                                      "Coulombinter": TABLES["CoulombInter"]})
        ctx = build_so_entry_context({}, tables, NORB, LOG)
        _same_records(self, ctx.compiled, so_records.compile_interaction(TABLES, NORB))

    def test_complex_refusal_propagates(self):
        # one family, exact text; every family is pinned in test_so_records
        tables = {"Hund": {((0, 0, 0), (0, 1)): 0.5 + 0.25j, ((0, 0, 0), (1, 0)): 0.5 - 0.25j}}
        with self.assertRaises(ValueError) as cm:
            build_so_entry_context({"orbital_site": [0, 1]}, tables, NORB,
                                   logging.getLogger("qlms.test_so_entry_rules.quiet"))
        self.assertEqual(
            str(cm.exception),
            "Hund coefficient at R=(0, 0, 0) orbitals (1, 2) has an imaginary part (0.25). "
            "In H-wave's declaration convention the entry and its partner (−R,b,a) multiply "
            "the same operator, so the imaginary part has no operator content; only PairHop "
            "carries a physical complex phase (complex couplings beyond PairHop: issue #302). "
            "Remove the imaginary part.")

    def test_complex_pairhop_accepted(self):
        tables = {"PairHop": {((0, 0, 0), (0, 1)): 0.3 + 0.2j, ((0, 0, 0), (1, 0)): 0.3 - 0.2j}}
        ctx = build_so_entry_context({}, tables, NORB, LOG)
        self.assertTrue(len(ctx.compiled.records) > 0)


class TestNested(unittest.TestCase):

    def test_context_for_nested(self):
        ctx = build_so_entry_context({}, TABLES, NORB, LOG)
        info_mode = {"mode": "RPA", "param": {"T": 0.1}}
        self.assertIsNone(context_for_nested(info_mode, ctx))
        self.assertEqual(SO_ENTRY_CONTEXT_KEY, "_so_entry_context")
        self.assertIs(info_mode["_so_entry_context"], ctx)
        self.assertNotIn("_so_entry_context", info_mode["param"])

    def test_context_type(self):
        self.assertIsInstance(build_so_entry_context({}, TABLES, NORB, LOG), SOEntryContext)

    def test_context_carries_the_interaction_fingerprints(self):
        from hwave.solver import fingerprint
        ctx = build_so_entry_context({}, TABLES, NORB, LOG)
        want = fingerprint.interaction_inputs(so_records.compile_interaction(TABLES, NORB))
        self.assertIsInstance(ctx.interaction, fingerprint.InteractionInputs)
        self.assertIs(ctx.interaction.records, ctx.compiled.records)
        self.assertIs(ctx.interaction.declared_topology, ctx.compiled.declared_topology)
        self.assertEqual((ctx.interaction.fp_interaction, ctx.interaction.fp_topology),
                         (want.fp_interaction, want.fp_topology))


class TestTablesDigest(unittest.TestCase):
    """The context records a digest of the raw interaction tables it
    compiled, so a nested solver can tell whether it builds its vertex
    from the same tables."""

    def test_context_records_the_tables_digest(self):
        ctx = build_so_entry_context({}, TABLES, NORB, LOG)
        self.assertEqual(ctx.tables_digest, interaction_tables_digest(TABLES))
        self.assertRegex(ctx.tables_digest, "^[0-9a-f]{64}$")

    def test_value_change_changes_the_digest(self):
        other = dict(TABLES, Hund=closed_row((0, 0, 0), 0, 1, 0.5 + 1e-12))
        self.assertNotEqual(interaction_tables_digest(other), interaction_tables_digest(TABLES))

    def test_type_change_changes_the_digest(self):
        # the same entries under another type name
        other = {"CoulombIntra": TABLES["CoulombIntra"], "Exchange": TABLES["Hund"],
                 "CoulombInter": TABLES["CoulombInter"]}
        self.assertNotEqual(interaction_tables_digest(other), interaction_tables_digest(TABLES))

    def test_non_interaction_tables_ignored(self):
        tables = dict(TABLES, Transfer={((0, 0, 0), (0, 0)): 1.0},
                      Extern={((0, 0, 0), (0, 0)): 1.0}, Hunds={((0, 0, 0), (0, 1)): 1.0})
        self.assertEqual(interaction_tables_digest(tables), interaction_tables_digest(TABLES))

    def test_mixed_case_type_keys(self):
        tables = CaseInsensitiveDict({"coulombintra": TABLES["CoulombIntra"],
                                      "HUND": TABLES["Hund"],
                                      "Coulombinter": TABLES["CoulombInter"]})
        self.assertEqual(interaction_tables_digest(tables), interaction_tables_digest(TABLES))

    def test_coulomb_aggregate_equals_its_split(self):
        # the RPA reader keeps the aggregate Coulomb table; the hwave_sc
        # reader splits it (R = 0 diagonal -> CoulombIntra, rest -> CoulombInter)
        inter = dict(TABLES["CoulombInter"])
        inter.update(closed_row((0, 0, 0), 0, 1, 0.9))
        coulomb = dict(TABLES["CoulombIntra"])
        coulomb.update(inter)
        split = {"CoulombIntra": TABLES["CoulombIntra"], "CoulombInter": inter,
                 "Hund": TABLES["Hund"]}
        self.assertEqual(interaction_tables_digest({"coulomb": coulomb, "Hund": TABLES["Hund"]}),
                         interaction_tables_digest(split))

    def test_readers_agree_on_a_coulomb_file(self):
        # the hwave_sc reader and the RPA reader of the same files give one digest
        import shutil
        import tempfile
        import hwave.sc as sc
        import tests.test_sc_spin_orbital as S
        from hwave.qlmsio.read_input_k import QLMSkInput
        d = tempfile.mkdtemp(prefix="so_entry_readers_")
        try:
            S._write_model(d, so=True, soc=0.3, types=("Hund",))
            rows = [(0, 0, 0, 1, 1, 1.6, 0.0), (0, 0, 0, 2, 2, 1.9, 0.0),
                    (0, 0, 0, 1, 2, 0.6, 0.0), (0, 0, 0, 2, 1, 0.6, 0.0),
                    (1, 0, 0, 1, 1, 0.2, 0.0), (-1, 0, 0, 1, 1, 0.2, 0.0),
                    (1, 0, 0, 1, 2, 0.1, 0.0), (-1, 0, 0, 2, 1, 0.1, 0.0)]
            S._wfile(os.path.join(d, "coulomb.dat"), "Coulomb", S.NP, rows)
            inter = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat",
                     "coulomb": "coulomb.dat", "HUND": "hund.dat"}
            ham = QLMSkInput({"path_to_input": d, "interaction": dict(inter)}).get_param("ham")
            _, _, sc_tables = sc._read_interaction_files(
                {"file": {"input": {"interaction": dict(inter)}}})
        finally:
            shutil.rmtree(d, ignore_errors=True)
        self.assertIn("CoulombInter", sc_tables)
        self.assertEqual(interaction_tables_digest(ham), interaction_tables_digest(sc_tables))


class TestOffsiteTypes(unittest.TestCase):

    def test_offsite_types(self):
        # NN CoulombInter + R = 0 cross-orbital Hund -> CoulombInter only
        self.assertEqual(offsite_types(so_records.compile_interaction(TABLES, NORB)),
                         ("CoulombInter",))

    def test_sorted_and_pairhop_excluded(self):
        # an off-site PairHop is dropped by the RPA vertex with its own
        # WARNING, so it is never named as entering through V(q)
        tables = {"PairHop": closed_row((1, 0, 0), 0, 1, 0.2),
                  "Ising": closed_row((0, 1, 0), 1, 1, 0.3),
                  "Exchange": closed_row((1, 0, 0), 0, 1, 0.1),
                  "CoulombInter": closed_row((0, 0, 0), 0, 1, 0.4)}
        self.assertEqual(offsite_types(so_records.compile_interaction(tables, NORB)),
                         ("Exchange", "Ising"))

    def test_zero_valued_offsite_row_not_named(self):
        tables = {"Hund": closed_row((1, 0, 0), 0, 1, 0.0),
                  "CoulombIntra": {((0, 0, 0), (0, 0)): 1.0}}
        self.assertEqual(offsite_types(so_records.compile_interaction(tables, NORB)), ())


class TestSCSupported(unittest.TestCase):
    """``so_pairing.check_so_supported`` (the hwave_sc spin-orbital entries):
    every ``R = 0`` row is accepted, an off-site row is refused with the
    exact text, and the ``site`` argument is ignored."""

    def _input(self, files=None):
        import tests.test_sc_spin_orbital as S
        inter = dict({"CoulombIntra": "u.dat", "CoulombInter": "v.dat"}, **(files or {}))
        return S._input("/nonexistent", inter, True,
                        {"chi0q_mode": "calc", "solver_mode": "eigenvalue"})

    #: R = 0 rows between the two orbitals (cross-site under orbital_site = [0, 1])
    CROSS = {"CoulombIntra": {((0, 0, 0), (0, 0)): 2.0, ((0, 0, 0), (1, 1)): 2.0},
             "CoulombInter": closed_row((0, 0, 0), 0, 1, 0.5),
             "Hund": closed_row((0, 0, 0), 0, 1, 0.1)}

    def _check(self, interactions, **kw):
        from hwave.solver import so_pairing
        return so_pairing.check_so_supported(self._input(), interactions, norb_phys=NORB, **kw)

    def test_cross_site_row_accepted(self):
        for site in (None, [0, 1], (0, 1)):
            with self.subTest(site=site):
                self.assertIsNone(self._check(self.CROSS, site=site))

    def test_site_argument_ignored(self):
        offsite = dict(self.CROSS, CoulombInter={**closed_row((0, 0, 0), 0, 1, 0.5),
                                                 **closed_row((1, 0, 0), 0, 1, 0.2)})
        for tables, refused in ((self.CROSS, False), (offsite, True)):
            seen = []
            for site in ([0, 1], None, "garbage", [0]):
                with self.subTest(refused=refused, site=site):
                    try:
                        self._check(tables, site=site)
                        seen.append(None)
                    except ValueError as exc:
                        seen.append(str(exc))
            self.assertEqual(len(set(seen)), 1, seen)
            self.assertEqual(seen[0] is not None, refused)
        from hwave.solver import so_pairing
        self.assertIsNone(inspect.signature(so_pairing.check_so_supported)
                          .parameters["site"].default)

    def test_offsite_text(self):
        tables = dict(self.CROSS, Hund=closed_row((0, 1, 0), 1, 0, 0.1))
        with self.assertRaises(ValueError) as cm:
            self._check(tables)
        self.assertEqual(
            str(cm.exception),
            "[mode] enable_spin_orbital = true (hwave_sc, issue #83): an off-site row "
            "(R = (0, 1, 0), orbitals (2, 1)) in the Hund file needs [eliashberg] "
            "bond_channels = true: the exchange crossing of an off-site term depends on "
            "k + k', which only the bond-resolved pairing path treats. Set bond_channels = "
            "true (with frequency = \"static\"), or remove the off-site rows.")
        import tests.test_sc_spin_orbital as S
        self.assertEqual(str(cm.exception), S.OFFSITE_TEXT.format(
            R=(0, 1, 0), orbs=(2, 1), src="the Hund file"))
        # bond_channels = false said explicitly: the same text
        from hwave.solver import so_pairing
        inp = self._input()
        inp["eliashberg"]["Bond_Channels"] = False
        with self.assertRaises(ValueError) as cm:
            so_pairing.check_so_supported(inp, tables, norb_phys=NORB)
        self.assertEqual(str(cm.exception), S.OFFSITE_TEXT.format(
            R=(0, 1, 0), orbs=(2, 1), src="the Hund file"))

    def test_offsite_dynamic_text_unchanged(self):
        from hwave.solver import so_pairing
        tables = dict(self.CROSS, Hund=closed_row((0, 1, 0), 1, 0, 0.1))
        inp = self._input()
        inp["eliashberg"]["Frequency"] = "dynamic"
        with self.assertRaises(ValueError) as cm:
            so_pairing.check_so_supported(inp, tables, norb_phys=NORB)
        self.assertEqual(
            str(cm.exception),
            "[mode] enable_spin_orbital = true (hwave_sc, issue #83): an off-site row "
            "(R = (0, 1, 0), orbitals (2, 1)) in the Hund file is not supported; the "
            "spin-orbital pairing vertex supports on-site and in-cell (R = 0) interactions "
            "only (the exchange crossing of an off-site term depends on k + k'). Remove the "
            "off-site rows.")
        import tests.test_sc_spin_orbital as S
        self.assertEqual(str(cm.exception), S.OFFSITE_DYNAMIC_TEXT.format(
            R=(0, 1, 0), orbs=(2, 1), src="the Hund file"))


# ---------------------------------------------------------------------------
# The spin-orbital RPA and FLEX public entries
# ---------------------------------------------------------------------------

#: the two-orbital square-lattice model of tests.test_sc_spin_orbital with
#: Rashba-like spin-orbit hopping on a 2 x 2 lattice (Nmat 16): on-site
#: CoulombIntra plus R = 0 cross-orbital CoulombInter and Hund rows
_TYPES = ("CoulombIntra", "CoulombInter", "Hund")
#: no row beyond the R = 0 ones (all on-site in spin-orbital mode)
CROSS_SITE = {}
#: CROSS_SITE plus one nearest-neighbour CoulombInter row (and its partner)
OFFSITE = {"CoulombInter": [(1, 0, 0, 0.2), (-1, 0, 0, 0.2)]}
#: on-site CoulombIntra only
ONSITE = None

#: pre-switch reference values of this model (chiq at frequency index 8,
#: q index 3): RPA element (0, 1, 0, 1) and FLEX element (0, 3, 0, 3), with
#: orbital_site = [0, 1] (the R = 0 cross-orbital rows treated as inter-site,
#: exchange crossing left out) and without the key
_PRE_SWITCH = {
    "rpa": {"index": (8, 3, 0, 1, 0, 1), "orbital_site": 0.452385798801902,
            "none": 0.5463003298961958},
    "flex": {"index": (8, 3, 0, 3, 0, 3), "orbital_site": 0.31578499105300845,
             "none": 0.3721658659310463},
}


class _Records(logging.Handler):
    def __init__(self):
        super().__init__(logging.DEBUG)
        self.records = []

    def emit(self, record):
        self.records.append((record.levelname, record.getMessage()))


def _warning(text):
    """``(level, message)`` of a captured WARNING record."""
    return ("WARNING", text)


def _info(text):
    """``(level, message)`` of a captured INFO record."""
    return ("INFO", text)


_CACHE = {}


def _numeric_members(root):
    out = {}
    for dp, _, fns in os.walk(root):
        for fn in sorted(fns):
            if fn.endswith(".npz"):
                with np.load(os.path.join(dp, fn), allow_pickle=False) as z:
                    for k in z.files:
                        a = z[k]
                        if a.dtype.kind in "biufc":
                            out[(os.path.relpath(os.path.join(dp, fn), root), k)] = a
    return out


def _all_numeric_equal(a, b):
    """The numeric members of two runs agree (:func:`numeric_mismatches`:
    eigenvalues to 1e-12, eigenvectors and gaps through their spans, other
    members to 1e-13)."""
    return bool(a) and numeric_mismatches(a, b) == []


#: the hwave_sc entries of :func:`_run` (hwave.sc.calc_eliashberg), besides
#: the spin-orbital RPA and FLEX entries "rpa" and "flex" (hwave.qlms.run)
SC_ENTRIES = ("sc_calc", "sc_load", "sc_flex", "sc_dynamic_calc", "sc_dynamic_flex")

#: compile_interaction calls of each cached run (same keys as _CACHE)
_COMPILES = {}
#: files left in the output directory by each cached run (same keys)
_LEFT = {}
#: {(archive, member): text} of the 0-d string members of the archives each
#: cached run wrote (same keys)
_STRINGS = {}


def _string_members(root):
    """``{(archive, member): text}`` of the 0-d string members of the
    archives in ``root``."""
    from hwave.solver import npy_header
    out = {}
    if os.path.isdir(root):
        for fn in sorted(os.listdir(root)):
            if fn.endswith(".npz"):
                with np.load(os.path.join(root, fn)) as z:
                    for k in z.files:
                        # decided from the npy header: no payload is read
                        with z.zip.open(k + ".npy") as fh:
                            h = npy_header.read_npy_header(fh)
                        if h.shape == () and np.dtype(h.dtype).kind == "U":
                            out[(fn, k)] = str(z[k][()])
    return out


def _write(d, model):
    import tests.test_sc_spin_orbital as S
    if isinstance(model, tuple) and model[0] == "complex":
        return _write_complex(d, *model[1:])
    types = ("CoulombIntra",) if model is None else _TYPES
    return S._write_model(d, so=True, soc=0.3, types=types, extra_rows=model)


def _write_complex(d, family, value):
    """The model of :func:`_write` with on-site CoulombIntra and one
    ``family`` file holding the R = 0 row (orbitals 1, 2) = ``value`` and
    its partner (orbitals 2, 1) = conj(``value``)."""
    import tests.test_sc_spin_orbital as S
    files = S._write_model(d, so=True, soc=0.3, types=("CoulombIntra",))
    value = complex(value)
    files[family] = family.lower() + ".dat"
    S._wfile(os.path.join(d, files[family]), family, S.NP,
             [(0, 0, 0, 1, 2, value.real, value.imag),
              (0, 0, 0, 2, 1, value.real, -value.imag)])
    return files


#: the hwave_sc runs use a 3 x 2 lattice: on 2 x 2 every sin k vanishes, the
#: spin-orbit hopping drops out and the RPA producer of the load entry would
#: write a physical-orbital (spin-diagonal) chi0q
SC_CELL = [3, 2, 1]
#: Matsubara frequencies of the hwave_sc runs (and their producers)
SC_NMAT = 8


def _qlms_input(d, files, entry, out, cell=(2, 2, 1), nmat=16, **param):
    import tests.test_sc_spin_orbital as S
    inp = S._input(d, files, True, {}, nmat=nmat, out=out, CellShape=list(cell), **param)
    inp["mode"]["mode"] = {"rpa": "RPA", "flex": "FLEX"}[entry]
    inp["mode"]["calc_scheme"] = "general"
    inp["file"]["input"]["path_to_input"] = d
    inp["file"]["output"]["chiq"] = "chiq"
    if entry == "flex":
        inp["mode"]["param"].update({"IterationMax": 200, "Mix": 0.5, "EPS": 10,
                                     "mixing_scheme": "linear",
                                     "flex_hartree_fock": True})
        inp["file"]["output"].update({"green": "green", "sigma": "sigma"})
    del inp["eliashberg"]
    inp["log"] = {"print_level": 0}
    return inp


#: producer runs of the hwave_sc load / flex entries, by (kind, model)
_PRODUCERS = {}


def _producer(kind, model):
    """The output directory of the spin-orbital ``kind`` (``rpa`` or
    ``flex``) run of ``model`` on :data:`SC_CELL` with :data:`SC_NMAT`
    frequencies, without any extra [mode.param] key; run once per process,
    removed at exit."""
    key = (kind, repr(model))
    if key not in _PRODUCERS:
        import atexit
        import contextlib
        import io
        import shutil
        import tempfile
        from hwave import qlms
        import tests.test_sc_spin_orbital as S
        d = tempfile.mkdtemp(prefix="so_entry_prod_")
        atexit.register(shutil.rmtree, d, True)
        prod = _qlms_input(d, _write(d, model), kind, "prod", cell=SC_CELL, nmat=SC_NMAT)
        if kind == "flex":
            prod["file"]["output"]["chi0q"] = "chi0q"
        with contextlib.redirect_stdout(io.StringIO()), S._quiet():
            qlms.run(input_dict=prod)
        _PRODUCERS[key] = os.path.join(d, "prod")
    return _PRODUCERS[key]


def _sc_input(d, files, entry, model, **param):
    """The hwave_sc input of ``entry`` (one of :data:`SC_ENTRIES`) in ``d``;
    a load entry reads a copy of the RPA producer's chi0q (the producer run
    of ``model``) from its own output directory, a flex entry the FLEX
    producer's outputs."""
    import shutil
    import tests.test_sc_spin_orbital as S
    chi0q_mode = entry.rsplit("_", 1)[1]
    eli = {"chi0q_mode": chi0q_mode, "solver_mode": "eigenvalue", "num_eigenvalues": 2}
    if "dynamic" in entry:
        eli["frequency"] = "dynamic"
    inp = S._input(d, files, True, eli, nmat=SC_NMAT, out="out", CellShape=list(SC_CELL),
                   **param)
    if chi0q_mode == "load":
        # the chi0q archive is read from path_to_output
        shutil.copytree(_producer("rpa", model), os.path.join(d, "out"))
    elif chi0q_mode == "flex":
        inp["file"]["input"]["path_to_flex_output"] = _producer("flex", model)
    return inp


def _cache_key(entry, model, refused, param):
    return (entry, repr(model), refused, repr(sorted(param.items())))


def _run(entry, model, refused=False, **param):
    """One run through a public entry (``hwave.qlms.run`` for ``rpa`` and
    ``flex``, ``hwave.sc.calc_eliashberg`` for :data:`SC_ENTRIES`):
    ``(numeric members of every archive in the output directory -- for
    hwave_sc also the returned eigenvalue --, [(level, message)] of every
    record logged)``; with ``refused`` the run must raise ValueError and the
    first member is that exception. Only the entry's own run is logged and
    counted (not the producer run of an hwave_sc load / flex entry, see
    :func:`_producer`). Cached on the arguments (each test also runs
    alone)."""
    key = _cache_key(entry, model, refused, param)
    if key in _CACHE:
        return _CACHE[key]
    import contextlib
    import io
    import shutil
    import tempfile
    from hwave import qlms
    import hwave.sc as sc
    d = tempfile.mkdtemp(prefix="so_entry_")
    try:
        files = _write(d, model)
        if entry in SC_ENTRIES:
            inp = _sc_input(d, files, entry, model, **param)
        else:
            inp = _qlms_input(d, files, entry, "out", **param)
        root = logging.getLogger()
        handler, level = _Records(), root.level
        root.addHandler(handler)
        root.setLevel(logging.INFO)
        error, lam, compiles = None, None, None
        real = so_records.compile_interaction
        try:
            with contextlib.redirect_stdout(io.StringIO()), \
                    mock.patch.object(spin_orbital.so_records, "compile_interaction",
                                      side_effect=real) as spy:
                try:
                    if entry in SC_ENTRIES:
                        lam = sc.calc_eliashberg(inp)
                    else:
                        qlms.run(input_dict=inp)
                finally:
                    compiles = spy.call_count
        except ValueError as exc:
            if not refused:
                raise
            error = exc
        finally:
            root.removeHandler(handler)
            root.setLevel(level)
        if refused and error is None:
            raise AssertionError("the run was not refused")
        out_dir = inp["file"]["output"]["path_to_output"]
        _LEFT[key] = sorted(os.path.relpath(os.path.join(dp, fn), out_dir)
                            for dp, _, fns in os.walk(out_dir) for fn in fns)
        if refused:
            members = error
        else:
            members = _numeric_members(inp["file"]["output"]["path_to_output"])
            if entry == "sc_load":
                # the producer's archives copied into the output directory
                # are inputs, not outputs of this run
                copied = set(_numeric_members(_producer("rpa", model)))
                members = {k: v for k, v in members.items() if k not in copied}
            if entry in SC_ENTRIES:
                members[("return", "lambda")] = np.asarray(lam)
        _STRINGS[key] = _string_members(inp["file"]["output"]["path_to_output"])
        result = (members, handler.records)
    finally:
        shutil.rmtree(d, ignore_errors=True)
    _CACHE[key] = result
    _COMPILES[key] = compiles
    return result


def _compiles(entry, model, refused=False, **param):
    """compile_interaction calls of the (cached) run :func:`_run`."""
    _run(entry, model, refused, **param)
    return _COMPILES[_cache_key(entry, model, refused, param)]


class TestRPAFlexEntry(unittest.TestCase):
    """The spin-orbital RPA and FLEX entries, and the hwave_sc entries
    (static and dynamic): ``orbital_site`` is ignored (one WARNING per
    public run, results equal to the run without the key), the
    interaction is compiled once at the entry (nested RPA runs included),
    and the off-site INFO line is logged by the RPA entry only."""

    #: the orbital_site values run per entry: valid, invalid type, wrong
    #: length (hwave_sc: the valid one here, the invalid ones in
    #: test_invalid_orbital_site_ignored)
    _SITE_VALUES = dict({"rpa": ([0, 1], "garbage", [0]), "flex": ([0, 1], "garbage", [0])},
                        **{entry: ([0, 1],) for entry in SC_ENTRIES})

    def test_orbital_site_ignored_once_and_unchanged(self):
        # one WARNING per public run (a nested RPA included), the
        # interaction compiled once, and every number equal to the run
        # without the key (tests.so_run_compare)
        for entry, values in self._SITE_VALUES.items():
            for value in values:
                with self.subTest(entry=entry, value=value):
                    ref, l0 = _run(entry, CROSS_SITE)
                    got, l1 = _run(entry, CROSS_SITE, Orbital_Site=value)
                    self.assertEqual(numeric_mismatches(ref, got), [])
                    self.assertEqual(l1.count(_warning(ORBITAL_SITE_WARNING)), 1)
                    self.assertEqual(l0.count(_warning(ORBITAL_SITE_WARNING)), 0)
                    self.assertEqual(_compiles(entry, CROSS_SITE, Orbital_Site=value), 1)
                    self.assertEqual(_compiles(entry, CROSS_SITE), 1)

    def test_invalid_orbital_site_ignored(self):
        # the hwave_sc entries (static and dynamic) never parse the value
        for entry in ("sc_calc", "sc_dynamic_calc"):
            for value in ("garbage", [0]):
                with self.subTest(entry=entry, value=value):
                    got, logs = _run(entry, CROSS_SITE, orbital_site=value)
                    self.assertEqual(numeric_mismatches(_run(entry, CROSS_SITE)[0], got), [])
                    self.assertEqual(logs.count(_warning(ORBITAL_SITE_WARNING)), 1)
                    self.assertEqual(_compiles(entry, CROSS_SITE, orbital_site=value), 1)

    def test_unused_key_still_reported(self):
        # the RPA and FLEX entries report a [mode.param] key the run never
        # read; orbital_site counts as read
        for entry in ("rpa", "flex"):
            for extra in ({}, {"Orbital_Site": [0, 1]}):
                with self.subTest(entry=entry, extra=extra):
                    logs = _run(entry, CROSS_SITE, Foo_Bar=1, **extra)[1]
                    lines = [m for _, m in logs if "not read by this" in m]
                    self.assertEqual(lines, [
                        "[mode.param] key(s) not read by this {} run: Foo_Bar. A "
                        "misspelled or misplaced key has no effect; keys that only "
                        "apply to other settings of the solver are listed as "
                        "well.".format(entry.upper())])

    def test_offsite_info_text(self):
        from hwave.solver.rpa import SO_OFFSITE_INFO
        self.assertEqual(
            SO_OFFSITE_INFO.format(types="CoulombInter, Hund"),
            "enable_spin_orbital: the off-site rows of CoulombInter, Hund enter the "
            "spinful RPA vertex through their direct (Hartree-like) part V(q) only; "
            "the exchange (crossing) part of these off-site rows is not included in "
            "the ring sum. The on-site and in-cell exchange is included. For an exact "
            "treatment of the off-site exchange use hwave_sc with [eliashberg] "
            "bond_channels = true (frequency = \"static\").")

    def test_offsite_info_counts(self):
        from hwave.solver.rpa import SO_OFFSITE_INFO
        logs = _run("rpa", OFFSITE)[1]
        self.assertEqual(logs.count(_info(SO_OFFSITE_INFO.format(types="CoulombInter"))), 1)
        marker = "enable_spin_orbital: the off-site rows"
        self.assertEqual(sum(marker in m for _, m in logs), 1)
        for entry, model in (("flex", ONSITE), ("flex", CROSS_SITE), ("rpa", ONSITE),
                             ("rpa", CROSS_SITE)):
            with self.subTest(entry=entry, model=model):
                self.assertFalse(any(marker in m for _, m in _run(entry, model)[1]))
        # FLEX refuses an off-site row (after its RPA set-up) without the line
        error, logs = _run("flex", OFFSITE, refused=True)
        self.assertIn("supports on-site and in-cell interactions only", str(error))
        self.assertFalse(any(marker in m for _, m in logs))

    def test_flex_class_gate(self):
        from hwave.solver.flex import FLEX
        from hwave.solver.rpa import RPA
        self.assertTrue(RPA._SO_OFFSITE_INFO)
        self.assertFalse(FLEX._SO_OFFSITE_INFO)

    def test_cross_site_crossing_now_included(self):
        # with orbital_site the R = 0 cross-orbital rows used to lose their
        # exchange crossing; now they are on-site: the result moves away
        # from the pre-switch value and equals the run without the key
        for entry in ("rpa", "flex"):
            with self.subTest(entry=entry):
                pre = _PRE_SWITCH[entry]
                got = _run(entry, CROSS_SITE, orbital_site=[0, 1])[0]
                ref = _run(entry, CROSS_SITE)[0]
                chiq = got[("chiq.npz", "chiq")]
                self.assertGreaterEqual(abs(chiq[pre["index"]].real - pre["orbital_site"]), 1e-3)
                key = ("chiq.npz", "chiq")
                self.assertEqual(numeric_mismatches({key: chiq}, {key: ref[key]}), [])

    def test_rpa_chiq_unchanged_without_orbital_site(self):
        for entry in ("rpa", "flex"):
            with self.subTest(entry=entry):
                pre = _PRE_SWITCH[entry]
                chiq = _run(entry, CROSS_SITE)[0][("chiq.npz", "chiq")]
                self.assertAlmostEqual(chiq[pre["index"]].real, pre["none"], delta=1e-10)


#: the hwave_sc entries that consume a producer's archives
_CONSUMERS = ("sc_load", "sc_flex", "sc_dynamic_flex")
#: the archives of the producers of the consumers
_PRODUCED = {"rpa": ("chi0q.npz", "chiq.npz"),
             "flex": ("chi0q.npz", "chiq.npz", "green.npz", "sigma.npz")}
#: the producer stage of each producer
_STAGE = {"rpa": "rpa_so", "flex": "flex_so_onsite"}
#: the payload members of the spin-orbital producer archives
_PAYLOADS = frozenset({"chi0q", "chiq", "green", "sigma", "sigma_static", "sigma_fluct"})
#: the fingerprint and state members of a version-1 archive (an unversioned
#: archive carries none of them)
#: the producer an unversioned archive names, by the consumer's chi0q_mode
_EXPECTED_PRODUCER = {"load": EXPECTED_RPA_PRODUCER, "flex": EXPECTED_FLEX_PRODUCER}


class _PayloadReads:
    """While active, records every payload member (:data:`_PAYLOADS`) read
    from any npz archive (``reads``: ``[(archive file name, member)]``); a
    npy header read for the structure check is not a payload read."""

    def __enter__(self):
        real = np.lib.npyio.NpzFile.__getitem__
        self.reads = []

        def spy(npz, key):
            if key in _PAYLOADS:
                name = getattr(getattr(npz, "zip", None), "filename", None)
                self.reads.append((os.path.basename(str(name)), key))
            return real(npz, key)
        self._patch = mock.patch.object(np.lib.npyio.NpzFile, "__getitem__", spy)
        self._patch.__enter__()
        return self

    def __exit__(self, *exc):
        return self._patch.__exit__(*exc)


def _make_old(path, variant):
    """Rewrite the archive at ``path`` as one written before this version:
    ``"v1"`` (fingerprint_version 1, without the interaction fingerprints)
    or ``"unversioned"`` (no fingerprint or state member)."""
    with np.load(path) as z:
        members = {k: z[k] for k in z.files}
    drop = set(V2_ONLY_MEMBERS) | (set(V1_MEMBERS) if variant == "unversioned" else set())
    members = {k: v for k, v in members.items() if k not in drop}
    if variant == "v1":
        members["fingerprint_version"] = np.int64(1)
    np.savez(path, **members)


def _consumer_dir(test, entry, consumer_model=CROSS_SITE):
    """``(input, root)``: the hwave_sc input of ``entry`` in a fresh
    directory, its interaction files those of ``consumer_model`` and its
    producer's archives (the producer run of :data:`CROSS_SITE`) copied to
    ``root``, which the input reads."""
    import shutil
    import tempfile
    d = tempfile.mkdtemp(prefix="so_entry_consumer_")
    test.addCleanup(shutil.rmtree, d, True)
    inp = _sc_input(d, _write(d, consumer_model), entry, CROSS_SITE)
    if entry == "sc_load":
        return inp, os.path.join(d, "out")
    root = os.path.join(d, "flex")
    shutil.copytree(inp["file"]["input"]["path_to_flex_output"], root)
    inp["file"]["input"]["path_to_flex_output"] = root
    return inp, root


def _refusal(inp):
    """``(SOArchiveError, payload reads)`` of the hwave_sc run of ``inp``;
    any other outcome fails."""
    import contextlib
    import io
    import hwave.sc as sc
    import tests.test_sc_spin_orbital as S
    from hwave.solver import fingerprint
    with contextlib.redirect_stdout(io.StringIO()), S._quiet(), _PayloadReads() as spy:
        try:
            sc.calc_eliashberg(inp)
        except fingerprint.SOArchiveError as exc:
            return exc, spy.reads
    raise AssertionError("the run was not refused")


class TestVersionTwo(unittest.TestCase):
    """Every spin-orbital archive is written at fingerprint_version 2 and
    carries the interaction fingerprints; the consumers refuse a version-1
    or unversioned archive, and an archive of another interaction, before
    any payload is read."""

    def _old_archive_refused(self, entry, archive, variant):
        inp, root = _consumer_dir(self, entry)
        path = os.path.join(root, archive)
        with np.load(path) as z:
            label = str(z["producer_label"])
        _make_old(path, variant)
        error, reads = _refusal(inp)
        producer = (label if variant == "v1"
                    else _EXPECTED_PRODUCER[entry.rsplit("_", 1)[1]])
        text = V1_REFUSAL_TEXT if variant == "v1" else UNVERSIONED_REFUSAL_TEXT
        self.assertEqual(str(error), text.format(path=path, producer=producer))
        self.assertEqual(reads, [])

    def test_old_archives_refused(self):
        cells = [("sc_load", "chi0q.npz")] + [
            (entry, name) for entry in ("sc_flex", "sc_dynamic_flex")
            for name in _PRODUCED["flex"]]
        for entry, name in cells:
            for variant in ("v1", "unversioned"):
                with self.subTest(entry=entry, archive=name, variant=variant):
                    self._old_archive_refused(entry, name, variant)

    def test_round_trips(self):
        from hwave.solver import fingerprint
        produced = {}
        for kind, names in _PRODUCED.items():
            root = _producer(kind, CROSS_SITE)
            for name in names:
                with self.subTest(producer=kind, archive=name):
                    path = os.path.join(root, name)
                    with np.load(path) as z:
                        self.assertEqual(
                            fingerprint.validate_structure(z, path, kind=name[:-4],
                                                           archive=z),
                            (name[:-4], _STAGE[kind], 2))
                        produced[(kind, name)] = (str(z["fp_interaction"]),
                                                  str(z["fp_green"]))
        for entry in _CONSUMERS:
            with self.subTest(entry=entry):
                members = _run(entry, CROSS_SITE)[0]
                strings = _STRINGS[_cache_key(entry, CROSS_SITE, False, {})]
                gaps = {fn for (fn, k), v in strings.items()
                        if k == "producer_stage" and v.startswith("sc_so_")}
                self.assertEqual(len(gaps), 1)
                (gap,) = gaps
                self.assertEqual(int(members[(gap, "fingerprint_version")]), 2)
                kind = "rpa" if entry == "sc_load" else "flex"
                fp_interaction, fp_green = produced[(kind, "chi0q.npz")]
                self.assertEqual(strings[(gap, "fp_interaction")], fp_interaction)
                if kind == "flex":
                    # the borrowed Green function, digested as consumed
                    self.assertEqual(strings[(gap, "fp_green")], fp_green)

    def test_interaction_mismatch_refused(self):
        # the consumer's interaction files (on-site CoulombIntra only) are
        # not the producer's (plus the R = 0 CoulombInter and Hund rows)
        for entry in _CONSUMERS:
            with self.subTest(entry=entry):
                inp, root = _consumer_dir(self, entry, consumer_model=ONSITE)
                path = os.path.join(root, "chi0q.npz" if entry == "sc_load" else "chiq.npz")
                with np.load(path) as z:
                    label = str(z["producer_label"])
                error, reads = _refusal(inp)
                self.assertTrue(str(error).startswith(
                    "'{}' (written by {}) does not match this run: interaction differs: "
                    "file ".format(path, label)), str(error))
                self.assertIn("(the same interaction definitions)", str(error))
                self.assertEqual(reads, [])

    def test_interaction_digest_reader_independent(self):
        # the RPA and FLEX readers and the hwave_sc reader of the same
        # interaction files record one digest, that of the files compiled here
        import shutil
        import tempfile
        import hwave.qlmsio.wan90 as wan90
        import tests.test_sc_spin_orbital as S
        from hwave.solver import fingerprint
        fps = {}
        for kind, names in _PRODUCED.items():
            for name in names:
                with np.load(os.path.join(_producer(kind, CROSS_SITE), name)) as z:
                    fps[(kind, name)] = str(z["fp_interaction"])
        for entry in ("sc_calc", "sc_dynamic_calc") + _CONSUMERS:
            _run(entry, CROSS_SITE)
            for (fn, k), v in _STRINGS[_cache_key(entry, CROSS_SITE, False, {})].items():
                if k == "fp_interaction":
                    fps[(entry, fn)] = v
        d = tempfile.mkdtemp(prefix="so_entry_tables_")
        try:
            tables = {t: wan90.read_w90(os.path.join(d, f))
                      for t, f in _write(d, CROSS_SITE).items()}
        finally:
            shutil.rmtree(d, ignore_errors=True)
        want = fingerprint.interaction_digest(so_records.compile_interaction(tables, S.NP).records)
        # every producer archive and every hwave_sc gap archive is compared
        self.assertEqual({w for w, _ in fps}, {"rpa", "flex", "sc_calc", "sc_dynamic_calc"}
                         | set(_CONSUMERS))
        self.assertEqual(len([w for w in fps if w[0] in ("rpa", "flex")]), 6)
        self.assertEqual(set(fps.values()), {want})

    def test_prefold_tables_under_sublattice(self):
        # the interaction digest is that of the pre-fold tables: a run with a
        # sublattice records the same records as the run without one
        plain_members = _run("rpa", CROSS_SITE)[0]
        plain = _STRINGS[_cache_key("rpa", CROSS_SITE, False, {})]
        members = _run("rpa", CROSS_SITE, SubShape=[2, 1, 1])[0]
        folded = _STRINGS[_cache_key("rpa", CROSS_SITE, False, {"SubShape": [2, 1, 1]})]
        self.assertEqual(tuple(members[("chi0q.npz", "prov_sub_shape")]), (2, 1, 1))
        for name in ("chi0q.npz", "chiq.npz"):
            with self.subTest(archive=name):
                self.assertEqual(folded[(name, "fp_interaction")],
                                 plain[(name, "fp_interaction")])
                for k in ("prov_rec_modes", "prov_rec_cells", "prov_rec_v", "prov_topology"):
                    self.assertEqual(members[(name, k)].tobytes(),
                                     plain_members[(name, k)].tobytes(), k)
                # physical orbitals of the original cell: 2 * NP spin-orbitals
                self.assertLess(int(members[(name, "prov_rec_modes")].max()), 2 * 2)

    def test_onsite_stages_write_the_sentinel(self):
        bond = ("delta_r", "reverse", "slot_convention", "max_shells")
        archives = {}
        for kind, names in _PRODUCED.items():
            for name in names:
                with np.load(os.path.join(_producer(kind, CROSS_SITE), name)) as z:
                    archives[(kind, name)] = (str(z["producer_stage"]),
                                              {k: z[k] for k in z.files})
        for entry in ("sc_calc", "sc_dynamic_calc"):
            members = _run(entry, CROSS_SITE)[0]
            strings = _STRINGS[_cache_key(entry, CROSS_SITE, False, {})]
            for (fn, k), v in strings.items():
                if k == "producer_stage":
                    m = {kk: vv for (f, kk), vv in members.items() if f == fn}
                    m.update({kk: np.str_(vv) for (f, kk), vv in strings.items() if f == fn})
                    archives[(entry, fn)] = (v, m)
        self.assertEqual(sorted(stage for stage, _ in archives.values()),
                         sorted(["rpa_so"] * 2 + ["flex_so_onsite"] * 4
                                + ["sc_so_static", "sc_so_dynamic"]))
        for where, (stage, m) in archives.items():
            with self.subTest(archive=where):
                self.assertEqual(str(m["fp_capability"]), "so_onsite")
                self.assertEqual(str(m["fp_resolved_topology"]), "none:onsite")
                self.assertFalse(set(bond) & set(m))

    def test_flex_green_and_sigma_markers(self):
        from hwave.solver import fingerprint
        root = _producer("flex", CROSS_SITE)
        with np.load(os.path.join(root, "green.npz")) as g:
            self.assertEqual(str(g["green_layout"]), fingerprint.GREEN_LAYOUT)
            self.assertEqual(str(g["green_layout"]), "b_f_k_o_o_v1")
            self.assertEqual(g["spin_orbital"].dtype, np.bool_)
            self.assertTrue(bool(g["spin_orbital"]))
        with np.load(os.path.join(root, "sigma.npz")) as s:
            self.assertEqual(s["spin_orbital"].dtype, np.bool_)
            self.assertTrue(bool(s["spin_orbital"]))
            self.assertNotIn("green_layout", s.files)
        # a run without enable_spin_orbital refuses this green archive by
        # its spin_orbital marker, which is checked first
        import hwave.sc as sc
        with self.assertRaises(ValueError) as cm:
            sc._refuse_so_green_file(os.path.join(root, "green.npz"))
        self.assertTrue(str(cm.exception).startswith(
            "Green-function file '{}' was written by a spin-orbital run (spin_orbital = True)"
            .format(os.path.join(root, "green.npz"))), str(cm.exception))

    def test_green_and_sigma_checked_before_any_payload(self):
        # a malformed green or sigma archive is refused by its structure
        # check, before any payload (chiq, chi0q, green, sigma) is read
        cases = (
            ("green.npz", {"green_layout": np.str_("sc_o_o_xyz_f")}, None,
             "green_layout is not the canonical layout"),
            ("green.npz", {}, "spin_orbital", "marker spin_orbital is missing"),
            ("sigma.npz", {}, "sigma_static", "sigma_static"),
        )
        for entry in ("sc_flex", "sc_dynamic_flex"):
            for name, update, drop, text in cases:
                with self.subTest(entry=entry, archive=name, drop=drop):
                    inp, root = _consumer_dir(self, entry)
                    path = os.path.join(root, name)
                    with np.load(path) as z:
                        members = {k: z[k] for k in z.files if k != drop}
                    members.update(update)
                    np.savez(path, **members)
                    error, reads = _refusal(inp)
                    self.assertIn("'{}'".format(path), str(error))
                    self.assertIn(text, str(error))
                    self.assertEqual(reads, [])

    def test_green_and_sigma_checked_against_the_run(self):
        # the green and sigma archives are compared with this run (band,
        # interaction) before any payload is read
        fake = "0" * 64
        cases = (("fp_band", "does not match this run: band differs: fp_band: file 000000000000"),
                 ("fp_interaction", "does not match this run: interaction differs"))
        for entry in ("sc_flex", "sc_dynamic_flex"):
            for name in ("green.npz", "sigma.npz"):
                for member, text in cases:
                    with self.subTest(entry=entry, archive=name, member=member):
                        inp, root = _consumer_dir(self, entry)
                        path = os.path.join(root, name)
                        with np.load(path) as z:
                            members = {k: z[k] for k in z.files}
                        members[member] = np.str_(fake)
                        np.savez(path, **members)
                        error, reads = _refusal(inp)
                        self.assertTrue(str(error).startswith("'{}' (written by ".format(path)),
                                        str(error))
                        self.assertIn(text, str(error))
                        self.assertEqual(reads, [])

    def test_green_state_checked_against_the_run(self):
        # a green archive of another state point (here: another filling) is
        # refused by the state check before any payload is read
        inp, root = _consumer_dir(self, "sc_flex")
        path = os.path.join(root, "green.npz")
        with np.load(path) as z:
            members = {k: z[k] for k in z.files}
        target = float(members["state_target"])
        members["state_target"] = np.float64(target + 0.1)
        members["state_density"] = np.float64(float(members["state_density"]) + 0.1)
        np.savez(path, **members)
        label = str(members["producer_label"])
        density = float(members["state_density"])
        error, reads = _refusal(inp)
        # the run's target: filling 0.4 times 2 N = 4 spin-orbitals
        self.assertEqual(str(error), (
            "'{p}' (written by {l}) does not match this run: state differs: target {t!r} "
            "(file) vs 1.6 (run) electrons per cell; density {d!r} (file) vs target 1.6 "
            "(run, tol 1.6e-08) electrons per cell. Archives are reused only at the "
            "identical state point; re-run {l} at this filling and temperature.").format(
                p=path, l=label, t=target + 0.1, d=density))
        self.assertEqual(reads, [])

    def test_green_archive_bound_to_its_green_function(self):
        # the green archive's own fp_green is that of the Green function the
        # run loads from it
        inp, root = _consumer_dir(self, "sc_flex")
        path = os.path.join(root, "green.npz")
        with np.load(path) as z:
            members = {k: z[k] for k in z.files}
        members["fp_green"] = np.str_("0" * 64)
        np.savez(path, **members)
        error, _ = _refusal(inp)
        self.assertTrue(str(error).startswith("'{}' (written by ".format(path)), str(error))
        self.assertIn("fp_green mismatch", str(error))
        self.assertIn("the archive records 000000000000...", str(error))

    def test_big_endian_sigma_accepted(self):
        # a sigma_static stored as big-endian complex128 passes the archive
        # check and the loader (which returns it in native order), and gives
        # the eigenvalue of the native file
        import contextlib
        import io
        import hwave.sc as sc
        import tests.test_sc_spin_orbital as S
        loader = sc._load_flex_so_sigma
        lam, loaded = {}, {}
        for order in ("<", ">"):
            inp, root = _consumer_dir(self, "sc_flex")
            path = os.path.join(root, "sigma.npz")
            with np.load(path) as z:
                members = {k: z[k] for k in z.files}
            members["sigma_static"] = members["sigma_static"].astype(order + "c16")
            np.savez(path, **members)
            with np.load(path) as z:
                self.assertEqual(z["sigma_static"].dtype.isnative, order == "<")

            def spy(*a, **k):
                loaded[order] = loader(*a, **k)
                return loaded[order]
            with contextlib.redirect_stdout(io.StringIO()), S._quiet(), \
                    mock.patch.object(sc, "_load_flex_so_sigma", side_effect=spy):
                lam[order] = sc.calc_eliashberg(inp)
        self.assertTrue(eigenvalues_close(lam[">"], lam["<"]))
        self.assertEqual(loaded[">"].dtype, np.dtype(np.complex128))
        self.assertTrue(loaded[">"].dtype.isnative)
        self.assertEqual(loaded[">"].tobytes(), loaded["<"].tobytes())


#: the public entries of the complex-refusal matrix (keys of _run)
_COMPLEX_ENTRIES = {"rpa": "rpa", "flex": "flex", "sc_static": "sc_calc",
                    "sc_dynamic": "sc_dynamic_calc"}
#: the complex coefficient of the R = 0 row (orbitals 1, 2)
_COMPLEX_VALUE = 0.5 + 0.25j


class TestComplexRefusalsPublic(unittest.TestCase):
    """Each non-PairHop family, written with an imaginary part in a real
    interaction file, is refused at each spin-orbital entry with the
    compiler's exact text and leaves no output; a complex PairHop runs."""

    FAMILIES = ("CoulombInter", "Hund", "Ising", "Exchange", "PairLift")
    ENTRIES = ("rpa", "flex", "sc_static", "sc_dynamic")

    def _refused(self, family, entry):
        from tests.test_so_records import _expected_message
        model = ("complex", family, _COMPLEX_VALUE)
        error, _ = _run(_COMPLEX_ENTRIES[entry], model, refused=True)
        self.assertEqual(str(error), _expected_message(family, (0, 0, 0), 1, 2, 0.25))
        self.assertEqual(_LEFT[_cache_key(_COMPLEX_ENTRIES[entry], model, True, {})], [])

    def test_full_matrix(self):
        for family in self.FAMILIES:
            for entry in self.ENTRIES:
                with self.subTest(family=family, entry=entry):
                    self._refused(family, entry)

    def test_complex_pairhop_runs(self):
        for entry in ("rpa", "sc_static"):
            with self.subTest(entry=entry):
                complex_run = _run(_COMPLEX_ENTRIES[entry],
                                   ("complex", "PairHop", _COMPLEX_VALUE))[0]
                real_run = _run(_COMPLEX_ENTRIES[entry], ("complex", "PairHop", 0.5))[0]
                self.assertTrue(complex_run)
                self.assertTrue(all(np.all(np.isfinite(v)) for v in complex_run.values()
                                    if v.dtype.kind in "fc"))
                # the imaginary part is read (it changes the result)
                self.assertFalse(_all_numeric_equal(complex_run, real_run))


class TestNestedRPAEntry(unittest.TestCase):
    """A solver built with a supplied context (a nested run) reuses it: no
    second compile, no WARNING, no off-site INFO."""

    def _ham(self, d, model):
        import tests.test_sc_spin_orbital as S
        from hwave.qlmsio.read_input_k import QLMSkInput
        files = S._write_model(d, so=True, soc=0.3, types=_TYPES, extra_rows=model)
        inter = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
        inter.update(files)
        return QLMSkInput({"path_to_input": d, "interaction": inter}).get_param("ham")

    def _mode(self, **extra):
        param = {"T": 0.25, "CellShape": [2, 2, 1], "SubShape": [1, 1, 1], "Nmat": 16,
                 "filling": 0.4, "orbital_site": [0, 1]}
        return dict({"mode": "RPA", "calc_scheme": "general", "enable_spin_orbital": True,
                     "param": param}, **extra)

    def _build(self, model, info_mode, context_tables=None):
        """``context_tables``: a callable of the reader tables returning the
        tables a context is built from and supplied (None: no context)."""
        import shutil
        import tempfile
        from hwave.solver import rpa
        d = tempfile.mkdtemp(prefix="so_entry_nested_")
        try:
            ham = self._ham(d, model)
            if context_tables is not None:
                ctx = build_so_entry_context({}, context_tables(ham), NORB, LOG)
                context_for_nested(info_mode, ctx)
            root = logging.getLogger()
            handler, level = _Records(), root.level
            root.addHandler(handler)
            root.setLevel(logging.INFO)
            real = so_records.compile_interaction
            try:
                with mock.patch.object(spin_orbital.so_records, "compile_interaction",
                                       side_effect=real) as spy:
                    solver = rpa.RPA(ham, {}, info_mode)
            finally:
                root.removeHandler(handler)
                root.setLevel(level)
        finally:
            shutil.rmtree(d, ignore_errors=True)
        return solver, spy.call_count, handler.records

    def test_public_entry_compiles_once(self):
        solver, calls, logs = self._build(OFFSITE, self._mode())
        self.assertEqual(calls, 1)
        self.assertIsNone(solver.ham_info.orbital_site)
        self.assertIsInstance(solver.ham_info.so_context, SOEntryContext)
        self.assertEqual(logs.count(_warning(ORBITAL_SITE_WARNING)), 1)
        self.assertEqual(sum("enable_spin_orbital: the off-site rows" in m for _, m in logs), 1)

    def test_supplied_context_is_reused(self):
        # the context is compiled from the reader tables the solver reads
        info_mode = self._mode()
        solver, calls, logs = self._build(OFFSITE, info_mode, context_tables=lambda ham: ham)
        ctx = info_mode[SO_ENTRY_CONTEXT_KEY]
        self.assertEqual(calls, 0)
        self.assertIs(solver.ham_info.so_context, ctx)
        self.assertIsNone(solver.ham_info.orbital_site)
        self.assertEqual(logs.count(_warning(ORBITAL_SITE_WARNING)), 0)
        self.assertFalse(any("enable_spin_orbital: the off-site rows" in m for _, m in logs))

    def test_context_from_other_tables_refused(self):
        # a context compiled from other tables would put fingerprints on the
        # archives that do not describe the vertex the solver builds
        other = {"CoulombInter": closed_row((1, 0, 0), 0, 0, 0.2)}
        with self.assertRaises(ValueError) as cm:
            self._build(OFFSITE, self._mode(), context_tables=lambda ham: other)
        self.assertEqual(
            str(cm.exception),
            "enable_spin_orbital: the interaction tables of this solver differ from the ones "
            "its entry context was compiled from, so the recorded interaction fingerprint "
            "would not describe the vertex. Build the context and the solver from the same "
            "interaction files.")

    def test_context_differing_in_one_value_refused(self):
        def nudged(ham):
            tables = {k: dict(v) for k, v in ham.items()
                      if k.lower() in ("coulombintra", "coulombinter", "hund")}
            key = next(iter(tables["CoulombInter"]))
            tables["CoulombInter"][key] += 1e-12
            return tables
        with self.assertRaises(ValueError) as cm:
            self._build(OFFSITE, self._mode(), context_tables=nudged)
        self.assertEqual(str(cm.exception), spin_orbital.TABLES_MISMATCH_TEXT)

    def test_non_so_still_parses_orbital_site(self):
        info_mode = self._mode()
        del info_mode["enable_spin_orbital"]
        import tests.test_orbital_site as O
        from hwave.qlmsio.read_input_k import QLMSkInput
        from hwave.solver.rpa import Interaction, Lattice
        import shutil
        import tempfile
        d = tempfile.mkdtemp(prefix="so_entry_nonso_")
        try:
            ham = QLMSkInput({"path_to_input": d,
                              "interaction": O._write_chain(d, 2)}).get_param("ham")
        finally:
            shutil.rmtree(d, ignore_errors=True)
        param = dict(info_mode["param"], CellShape=[O.NSITE // 2, 1, 1])
        inter = Interaction(Lattice(CaseInsensitiveDict(param)), ham,
                            dict(info_mode, param=param))
        self.assertEqual(inter.orbital_site, (0, 1))
        self.assertIsNone(getattr(inter, "so_context", None))


if __name__ == "__main__":
    unittest.main()
