"""The solver facade of ``hwave.initial``: problem loading, block detection
and the Hartree-Fock potential (``hwave.solver.initial_facade``).

The block detection and the interaction part of ``UHFk._make_ham`` were
lifted into side-effect-free functions that UHFk now calls. Each lifted
result is compared, in the same process and bit for bit, with the frozen
d72f35f6 method bodies of ``tests/uhfk_legacy_d72f35f6.py`` (never with the
delegated UHFk method, which would be circular), and with the golden file
``tests/initial_ref/facade_golden_d72f35f6.npz`` captured at d72f35f6 by
``tests/initial_reference_capture.py --facade`` (integers exactly, floats to
1e-12, so the check holds across numpy and BLAS builds).
"""
import dataclasses
import os
import pathlib
import shutil
import tempfile
import unittest

import numpy as np

from tests import initial_reference_capture as capture
from tests import uhfk_legacy_d72f35f6 as legacy

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _assert_bitwise(tc, got, want, msg=""):
    got = np.asarray(got)
    want = np.asarray(want)
    tc.assertEqual(got.shape, want.shape, msg)
    tc.assertEqual(got.dtype, want.dtype, msg)
    tc.assertEqual(got.tobytes(), want.tobytes(), msg)


def _blocks_equal(tc, got, want, msg):
    tc.assertEqual(len(got), len(want), msg)
    for a, b in zip(got, want):
        _assert_bitwise(tc, a, b, msg)


class _Cases(unittest.TestCase):
    """The golden inputs, written once per class under a temporary dir."""

    @classmethod
    def setUpClass(cls):
        cls.root = tempfile.mkdtemp(prefix="initial_facade_")
        cls.cases = capture.facade_cases(cls.root)
        cls.golden = np.load(capture.FACADE_GOLDEN, allow_pickle=False)

    @classmethod
    def tearDownClass(cls):
        cls.golden.close()
        shutil.rmtree(cls.root, ignore_errors=True)

    def golden_case(self, name):
        g = self.golden
        sizes = g[name + "__block_sizes"]
        index = g[name + "__block_index"]
        blocks = np.split(index, np.cumsum(sizes)[:-1])
        return {"blocks": blocks,
                "block_to_group": g[name + "__block_to_group"].tolist(),
                "group_nconds": g[name + "__group_nconds"].tolist(),
                "ham": g[name + "__ham"]}

    def _variant(self, edit, name="uhfk_2sz0"):
        """A copy of input ``name`` whose parsed TOML ``edit`` changed."""
        src = self.cases[name]
        d = tempfile.mkdtemp(dir=self.root)
        dest = os.path.join(d, "case")
        shutil.copytree(os.path.dirname(src), dest)
        toml = os.path.join(dest, "input.toml")
        doc = capture._read_toml(toml)
        edit(doc)
        with open(toml, "w") as f:
            f.write(capture.dump_toml(doc))
        return toml


class TestGoldenFile(_Cases):
    def test_golden_covers_every_case_at_the_revision(self):
        self.assertEqual(tuple(self.golden["cases"].tolist()), capture.FACADE_CASE_NAMES)
        self.assertTrue(str(self.golden["revision"]).startswith(capture.REVISION))
        self.assertEqual(set(self.cases), set(capture.FACADE_CASE_NAMES))


class TestBlocks(_Cases):
    def test_blocks_identical_on_golden_inputs(self):
        from hwave.solver import initial_facade as fac
        for name, toml in self.cases.items():
            with self.subTest(case=name):
                problem = fac.load_problem(toml)
                sb = fac.solver_blocks(problem)
                solver, _, _ = fac._build_solver(toml)
                view = legacy.LegacyView(solver)
                view._detect_blocks()
                gold = self.golden_case(name)
                # facade == frozen legacy == golden, exactly
                _blocks_equal(self, sb.block_info, view.block_info, name)
                _blocks_equal(self, problem.blocks.block_info, view.block_info, name)
                self.assertEqual(list(sb.block_to_group), view.block_to_group)
                self.assertEqual(list(sb.group_nconds), view.group_nconds)
                self.assertEqual(len(sb.block_info), len(gold["blocks"]))
                for a, b in zip(sb.block_info, gold["blocks"]):
                    np.testing.assert_array_equal(a, b)
                self.assertEqual(list(sb.block_to_group), gold["block_to_group"])
                self.assertEqual(list(sb.group_nconds), gold["group_nconds"])
                # the delegated UHFk attributes: the same values and types
                _blocks_equal(self, solver.block_info, view.block_info, name)
                self.assertIsInstance(solver.block_info, list)
                self.assertEqual(solver.block_to_group, view.block_to_group)
                self.assertIsInstance(solver.block_to_group, list)
                self.assertEqual(solver.group_nconds, view.group_nconds)
                self.assertIsInstance(solver.group_nconds, list)

    def test_solver_blocks_is_frozen_and_read_only(self):
        from hwave.solver import initial_facade as fac
        sb = fac.load_problem(self.cases["uhfk_2sz0"]).blocks
        with self.assertRaises(dataclasses.FrozenInstanceError):
            sb.block_info = ()
        with self.assertRaises(ValueError):
            sb.block_info[0][0] = 99
        self.assertIsInstance(sb.block_info, tuple)
        self.assertIsInstance(sb.block_to_group, tuple)
        self.assertIsInstance(sb.group_nconds, tuple)

    def test_mixed_block_refusal_kept(self):
        """2Sz fixed with a spin-flip one-body term: the lifted detection
        refuses it with the solver's message."""
        from hwave.solver import uhfk_blocks
        nd, norb = 4, 2
        ham_trans = np.zeros((1, nd, nd), dtype=np.complex128)
        ham_trans[0] = np.eye(nd)
        ham_trans[0, 0, 2] = ham_trans[0, 2, 0] = 0.5
        with self.assertRaisesRegex(ValueError, "spin-flip one-body"):
            uhfk_blocks.detect_blocks(
                ham_trans, {}, {}, nd=nd, norb=norb, ns=2, iflag_fock=True,
                sz_free=False, enable_spin_orbital=False, nconds=(1, 1))


class TestMakeHam(_Cases):
    def test_make_ham_bitwise(self):
        from hwave.solver import initial_facade as fac
        for name, toml in self.cases.items():
            with self.subTest(case=name):
                solver, green_info, _ = fac._build_solver(toml, with_green=True)
                solver.Green = solver._initial_green(green_info)
                view = legacy.LegacyView(solver)
                view._make_ham()
                solver._make_ham()
                _assert_bitwise(self, solver.ham, view.ham, name)
                gold = self.golden_case(name)["ham"]
                scale = max(1.0, float(np.abs(gold).max()))
                self.assertEqual(solver.ham.shape, gold.shape)
                self.assertLessEqual(float(np.abs(solver.ham - gold).max()),
                                     1e-12 * scale, name)

    def test_make_ham_bitwise_random_green(self):
        """The frozen and the delegated ``_make_ham`` on a random G of wide
        dynamic range, Fock on and off."""
        from hwave.solver import initial_facade as fac
        rng = np.random.default_rng(364)
        for name in ("uhfk_exchange", "so_all_types", "pairhop_fock_off",
                     "w90spin_flip", "uhfk_2sz0"):
            solver, _, _ = fac._build_solver(self.cases[name])
            shape = (solver.nvol, solver.ns, solver.norb, solver.ns, solver.norb)
            g = rng.normal(size=shape) + 1j * rng.normal(size=shape)
            solver.Green = g * 10.0 ** rng.uniform(-6, 2, size=shape)
            for fock in (True, False):
                with self.subTest(case=name, fock=fock):
                    solver.iflag_fock = fock
                    view = legacy.LegacyView(solver)
                    view._make_ham()
                    solver._make_ham()
                    _assert_bitwise(self, solver.ham, view.ham, name)


def _random_hermitian_green(solver, rng):
    """A G(r) whose G(k) is Hermitian at every k, in the stored layout."""
    nx, ny, nz = solver.shape
    nd = solver.nd
    a = rng.normal(size=(solver.nvol, nd, nd)) + 1j * rng.normal(size=(solver.nvol, nd, nd))
    gk = 0.5 * (a + np.conj(np.swapaxes(a, 1, 2)))
    gr = np.fft.fftn(gk.reshape(nx, ny, nz, nd, nd), axes=(0, 1, 2), norm="forward")
    return gr.reshape(solver.nvol, solver.ns, solver.norb, solver.ns, solver.norb)


class TestHFPotential(_Cases):
    def test_hf_potential_equals_make_ham_minus_transfer(self):
        from hwave.solver import initial_facade as fac
        rng = np.random.default_rng(5)
        for name in ("uhfk_exchange", "so_all_types", "pairhop_fock_off",
                     "w90spin_flip", "uhfk_2sz1"):
            problem = fac.load_problem(self.cases[name])
            solver, _, _ = fac._build_solver(self.cases[name])
            for fock in (True, False):
                with self.subTest(case=name, fock=fock):
                    g = _random_hermitian_green(solver, rng)
                    solver.Green = g
                    solver.iflag_fock = fock
                    solver._make_ham()
                    sigma = fac.hf_potential(g, fac.hf_tables(problem),
                                             include_fock=fock, blocks=problem.blocks)
                    want = solver.ham - problem.ham_trans_k
                    self.assertEqual(sigma.shape, (problem.nvol, problem.nd, problem.nd))
                    self.assertEqual(sigma.dtype, np.complex128)
                    tol = 1e-14 * max(1.0, float(np.abs(solver.ham).max()))
                    self.assertLessEqual(float(np.abs(sigma - want).max()), tol)
                    # the input G is not modified
                    _assert_bitwise(self, g, solver.Green)

    def test_accumulate_in_place_matches_make_ham(self):
        from hwave.solver import initial_facade as fac
        for name in ("uhfk_coulombintra", "so_all_types"):
            with self.subTest(case=name):
                solver, green_info, _ = fac._build_solver(self.cases[name], with_green=True)
                solver.Green = solver._initial_green(green_info)
                problem = fac.load_problem(self.cases[name])
                out = np.zeros((problem.nvol, problem.nd, problem.nd), dtype=np.complex128)
                out += problem.ham_trans_k
                ret = fac.accumulate_hf_projected(
                    out, solver.Green, problem.tables,
                    include_fock=problem.iflag_fock, blocks=problem.blocks)
                self.assertIsNone(ret)
                solver._make_ham()
                _assert_bitwise(self, out, solver.ham, name)


class TestDelegation(_Cases):
    """UHFk reaches its blocks and its Hartree-Fock potential through the
    lifted functions (a re-inlined copy in UHFk would still match the
    frozen oracle; these sentinels would not)."""

    def test_uhfk_calls_the_lifted_functions(self):
        from unittest import mock
        from hwave.solver import initial_facade as fac
        from hwave.solver import uhfk_blocks
        for name in ("uhfk_2sz0", "so_all_types"):
            with self.subTest(case=name):
                with mock.patch.object(uhfk_blocks, "detect_blocks",
                                       wraps=uhfk_blocks.detect_blocks) as det:
                    solver, green_info, _ = fac._build_solver(self.cases[name],
                                                              with_green=True)
                self.assertEqual(det.call_count, 1)
                self.assertIs(det.call_args.args[0], solver.ham_trans)
                solver.Green = solver._initial_green(green_info)
                with mock.patch.object(uhfk_blocks, "accumulate_hf_projected",
                                       wraps=uhfk_blocks.accumulate_hf_projected) as acc:
                    solver._make_ham()
                self.assertEqual(acc.call_count, 1)
                self.assertIs(acc.call_args.args[0], solver.ham)
                self.assertIs(acc.call_args.args[1], solver.Green)
                self.assertIs(acc.call_args.kwargs["blocks"], solver.block_info)

    def test_sentinels_reach_uhfk(self):
        from unittest import mock
        from hwave.solver import initial_facade as fac
        from hwave.solver import uhfk_blocks
        solver, green_info, _ = fac._build_solver(self.cases["uhfk_exchange"],
                                                  with_green=True)
        solver.Green = solver._initial_green(green_info)
        boom = RuntimeError("sentinel")
        with mock.patch.object(uhfk_blocks, "detect_blocks", side_effect=boom):
            with self.assertRaisesRegex(RuntimeError, "sentinel"):
                solver._detect_blocks()
        with mock.patch.object(uhfk_blocks, "accumulate_hf_projected", side_effect=boom):
            with self.assertRaisesRegex(RuntimeError, "sentinel"):
                solver._make_ham()


class TestShapeChecks(_Cases):
    def test_mismatched_tables_raise(self):
        from hwave.solver import initial_facade as fac
        toml = self.cases["uhfk_coulombinter"]
        p1 = fac.load_problem(toml, subshape=(1, 1, 1))
        p2 = fac.load_problem(toml, subshape=(2, 2, 1))
        so = fac.load_problem(self.cases["so_all_types"])
        g1 = np.zeros((p1.nvol, p1.ns, p1.norb, p1.ns, p1.norb), dtype=np.complex128)
        g2 = np.zeros((p2.nvol, p2.ns, p2.norb, p2.ns, p2.norb), dtype=np.complex128)
        gso = np.zeros((so.nvol, so.ns, so.norb, so.ns, so.norb), dtype=np.complex128)
        # the matching calls work
        fac.hf_potential(g1, p1.tables, include_fock=True, blocks=p1.blocks)
        fac.hf_potential(gso, so.tables, include_fock=True, blocks=so.blocks)
        bad = [
            (g1, p2), (g2, p1),                       # another SubShape
            (gso, p1),                                # SO G, normal tables
            (g1[:, :1, :, :1, :], p1),                # wrong spin axis
            (g2[..., :1], p2),                        # wrong orbital count
            (g1[0], p1),                              # not 5-d
        ]
        for i, (g, p) in enumerate(bad):
            with self.subTest(i=i):
                with self.assertRaises(ValueError):
                    fac.hf_potential(g, p.tables, include_fock=True, blocks=p.blocks)
        # an SO problem's tables with a normal-layout G of the same nd
        g_norm = np.zeros((so.nvol, 2, so.nd // 2, 2, so.nd // 2), dtype=np.complex128)
        with self.assertRaisesRegex(ValueError, "spin-orbital"):
            fac.hf_potential(g_norm, so.tables, include_fock=True, blocks=so.blocks)
        # accumulate_hf_projected checks the Hamiltonian too
        out = np.zeros((p2.nvol, p2.nd, p2.nd), dtype=np.complex128)
        with self.assertRaisesRegex(ValueError, "Hamiltonian"):
            fac.accumulate_hf_projected(out, g1, p1.tables, include_fock=True,
                                        blocks=p1.blocks)


class TestProblem(_Cases):
    def test_fields(self):
        from hwave.solver import initial_facade as fac
        from hwave.solver import hamiltonian_k
        for name, toml in self.cases.items():
            with self.subTest(case=name):
                p = fac.load_problem(toml)
                solver, _, _ = fac._build_solver(toml)
                for field in ("cellshape", "subshape", "shape", "nvol", "norb",
                              "norb_orig", "norb_phys", "ns", "nd",
                              "enable_spin_orbital", "iflag_fock", "sz_free", "T",
                              "spin_resolved_transfer"):
                    self.assertEqual(getattr(p, field), getattr(solver, field), field)
                self.assertEqual(p.nconds, tuple(solver.Nconds))
                self.assertIsInstance(p.nconds, tuple)
                _assert_bitwise(self, p.ham_trans_k, solver.ham_trans)
                self.assertEqual(p.ham_trans_r.shape, p.shape + (p.nd, p.nd))
                np.testing.assert_array_equal(
                    hamiltonian_k.to_k(p.ham_trans_r), p.ham_trans_k)
                self.assertIs(p.tables.enable_spin_orbital, p.enable_spin_orbital)
                self.assertEqual(p.tables.shape, p.shape)
                self.assertEqual(p.tables.norb_phys, p.norb_phys)
                self.assertEqual(p.geometry["norb"], p.norb_orig)
                self.assertEqual(p.param_ham_orig["Geometry"]["norb"], p.norb_orig)
                self.assertTrue(p.input_files)
                for f in p.input_files:
                    self.assertTrue(f.is_absolute() and f.is_file(), f)
                with self.assertRaises(dataclasses.FrozenInstanceError):
                    p.T = 1.0
                self._assert_deeply_read_only(p)
        self.assertTrue(fac.load_problem(self.cases["w90spin_flip"]).spin_resolved_transfer)
        self.assertTrue(fac.load_problem(self.cases["so_all_types"]).enable_spin_orbital)
        self.assertFalse(fac.load_problem(self.cases["pairhop_fock_off"]).iflag_fock)
        self.assertEqual(fac.load_problem(self.cases["uhfk_2sz0"]).nconds, (8, 8))

    def _assert_deeply_read_only(self, p):
        from collections.abc import Mapping

        def walk(obj, where):
            if isinstance(obj, np.ndarray):
                self.assertFalse(obj.flags.writeable, where)
            elif isinstance(obj, Mapping):
                with self.assertRaises(TypeError, msg=where):
                    obj["__new__"] = 1
                for k, v in obj.items():
                    walk(v, "{}[{!r}]".format(where, k))
            elif isinstance(obj, (list, dict, set)):
                self.fail("mutable {} at {}".format(type(obj).__name__, where))
            elif isinstance(obj, tuple):
                for i, v in enumerate(obj):
                    walk(v, "{}[{}]".format(where, i))

        for field in dataclasses.fields(p):
            value = getattr(p, field.name)
            if field.name in ("tables", "blocks"):
                for sub in dataclasses.fields(value):
                    walk(getattr(value, sub.name), "{}.{}".format(field.name, sub.name))
            else:
                walk(value, field.name)

    def test_problem_tables_are_copies(self):
        """Writing is refused, and the problem's tables are not the
        solver's (a later solver step cannot change them either)."""
        from hwave.solver import initial_facade as fac
        p = fac.load_problem(self.cases["uhfk_exchange"])
        with self.assertRaises(ValueError):
            p.tables.inter_table["Exchange"][0, 0, 0, 0, 0] = 1.0
        with self.assertRaises(ValueError):
            p.geometry["center"][0, 0] = 1.0
        with self.assertRaises(TypeError):
            p.param_ham_orig["Transfer"][((0, 0, 0), (0, 0))] = 1.0
        # case-insensitive like the solver's container
        self.assertIs(p.param_ham_orig["transfer"], p.param_ham_orig["Transfer"])

    def test_hf_tables_is_the_problem_tables(self):
        from hwave.solver import initial_facade as fac
        p = fac.load_problem(self.cases["uhfk_exchange"])
        t = fac.hf_tables(p)
        self.assertIs(t, p.tables)
        self.assertIsNotNone(t.inter_table["Exchange"])
        with self.assertRaises(TypeError):
            t.inter_table["Exchange"] = None

    def test_orbital_site(self):
        from hwave.solver import initial_facade as fac
        self.assertIsNone(fac.load_problem(self.cases["uhfk_coulombintra"]).orbital_site)

        def labelled(doc):
            doc["mode"]["param"]["orbital_site"] = [0]

        p = fac.load_problem(self._variant(labelled, "uhfk_coulombintra"))
        self.assertEqual(p.orbital_site, (0,))

        def labelled_so(doc):
            doc["mode"]["param"]["orbital_site"] = [0, 0, 1]

        # spin-orbital mode never parses the key (as the RPA solver)
        p = fac.load_problem(self._variant(labelled_so, "so_all_types"))
        self.assertIsNone(p.orbital_site)

    def test_subshape_override_independent(self):
        """Two loads of one TOML with different SubShape do not share
        state (the solver normalises its tables in place)."""
        from hwave.solver import initial_facade as fac
        toml = self.cases["uhfk_coulombinter"]
        p1 = fac.load_problem(toml, subshape=(1, 1, 1))
        t1 = p1.ham_trans_k.copy()
        p2 = fac.load_problem(toml, subshape=(2, 2, 1))
        p3 = fac.load_problem(toml)
        self.assertEqual(p1.subshape, (1, 1, 1))
        self.assertEqual(p1.shape, (4, 4, 1))
        self.assertEqual(p2.subshape, (2, 2, 1))
        self.assertEqual(p2.shape, (2, 2, 1))
        self.assertEqual(p2.norb, 4 * p1.norb)
        _assert_bitwise(self, p1.ham_trans_k, t1)
        self.assertEqual(p3.subshape, (2, 2, 1))      # the TOML's own value
        _assert_bitwise(self, p3.ham_trans_k, p2.ham_trans_k)
        self.assertEqual(sorted(p1.param_ham_orig["Transfer"].items()),
                         sorted(p2.param_ham_orig["Transfer"].items()))

    def test_relative_paths_resolve_from_the_toml_directory(self):
        from hwave.solver import initial_facade as fac
        toml = self.cases["uhfk_2sz1"]
        cwd = os.getcwd()
        os.chdir(tempfile.gettempdir())
        try:
            p = fac.load_problem(toml)
        finally:
            os.chdir(cwd)
        self.assertEqual(p.input_files[0].parent, pathlib.Path(toml).parent)
        self.assertEqual(p.input_files[0].name, "geom.dat")


class TestConfigLayers(_Cases):

    def test_config_layers_like_solver(self):
        from hwave.qlms import check_mode_tables
        from hwave.solver import initial_facade as fac

        def respell(key):
            def edit(doc):
                param = doc["mode"]["param"]
                param[key] = param.pop("SubShape")
            return edit

        base = fac.load_problem(self.cases["uhfk_2sz0"])
        for key in ("subshape", "SUBSHAPE"):
            with self.subTest(key=key):
                p = fac.load_problem(self._variant(respell(key)))
                self.assertEqual(p.subshape, base.subshape)
                _assert_bitwise(self, p.ham_trans_k, base.ham_trans_k)
                _blocks_equal(self, p.blocks.block_info, base.blocks.block_info, key)

        def mode_table_case(doc):
            doc["Mode"] = doc.pop("mode")

        def param_case(doc):
            doc["mode"]["Param"] = doc["mode"].pop("param")

        def mode_key_case(doc):
            doc["mode"]["Mode"] = doc["mode"].pop("mode")

        for edit in (param_case, mode_key_case):
            with self.subTest(edit=edit.__name__):
                toml = self._variant(edit)
                doc = capture._read_toml(toml)
                msgs = check_mode_tables(doc["mode"])
                self.assertTrue(msgs)
                with self.assertRaises(ValueError) as cm:
                    fac.load_problem(toml)
                for m in msgs:
                    self.assertIn(m, str(cm.exception))
        with self.assertRaisesRegex(ValueError, r"\[Mode\].*\[mode\]"):
            fac.load_problem(self._variant(mode_table_case))

        def initial_case(doc):
            doc["file"]["input"]["Initial"] = doc["file"]["input"].pop("initial")

        with self.assertRaisesRegex(ValueError, r"Initial.*'initial'"):
            fac.load_problem(self._variant(initial_case))

        def output_case(doc):
            doc["file"]["output"]["Energy"] = doc["file"]["output"].pop("energy")

        with self.assertRaisesRegex(ValueError, r"Energy.*'energy'"):
            fac.load_problem(self._variant(output_case))

        uhfr = os.path.join(_REPO, "tests", "uhfr", "coulombintra", "input.toml")
        with self.assertRaisesRegex(ValueError, "UHFr"):
            fac.load_problem(uhfr)

        def no_mode(doc):
            del doc["mode"]["mode"]

        with self.assertRaisesRegex(ValueError, "mode"):
            fac.load_problem(self._variant(no_mode))

    def test_solver_exit_becomes_value_error(self):
        from hwave.solver import initial_facade as fac

        def bad_subshape(doc):
            doc["mode"]["param"]["SubShape"] = [3, 1, 1]

        with self.assertRaisesRegex(ValueError, "UHFk refused"):
            fac.load_problem(self._variant(bad_subshape))


if __name__ == "__main__":
    unittest.main()
