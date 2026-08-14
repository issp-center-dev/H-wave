"""Task 6: ``sc.py`` wiring of the bond-resolved Eliashberg path.

Covers, per
``docs/superpowers/specs/2026-07-25-bond-channels-eliashberg-design.md``:

* S5  -- config keys (``bond_channels``, ``bond_max_shells``,
  ``bond_memory_cap_gb``), the ignore-with-warning rule, and every guard.
* S4.3 star -- the Case-2 correction (only ``V_ab(R=0)`` stays in the m=0 Fock
  block; the Hartree ``2 V_ab(q)`` stays in the m=0 charge block).
* S3.2 -- the resource preflight.
* S7.3 -- regression: with the flag off the default path is **bit-identical**.
* S7.4 -- ``V=0`` (declared-but-zero NN topology) reduces to the default path.

Tests must be run from the repository root (fixture paths are relative).
"""
import logging
import os
import shutil
import tempfile
import unittest
import zipfile

import numpy as np
from scipy.sparse.linalg import LinearOperator

import hwave.sc as sc
from tests.approx_util import ApproxTestCase, _caplog


# ---------------------------------------------------------------------------
# fixtures / helpers
# ---------------------------------------------------------------------------

RPA_INPUT = "tests/rpa/input"


def _write_w90(path, norb, entries, header="generated for test"):
    """Write an H-wave "wannier90-like" interaction file.

    ``entries`` is a list of ``(irvec, (orb1_0based, orb2_0based), value)``.
    """
    with open(path, "w") as fw:
        fw.write(header + "\n")
        fw.write("{}\n".format(norb))
        fw.write("{}\n".format(len(entries)))
        fw.write(" ".join(["1"] * len(entries)) + "\n")
        for irvec, orbvec, value in entries:
            value = complex(value)
            fw.write("{:4d} {:4d} {:4d} {:4d} {:4d} {:.12f} {:.12f}\n".format(
                irvec[0], irvec[1], irvec[2], orbvec[0] + 1, orbvec[1] + 1,
                value.real, value.imag))


def _nn_entries(value):
    return [((-1, 0, 0), (0, 0), value),
            ((0, -1, 0), (0, 0), value),
            ((0, 1, 0), (0, 0), value),
            ((1, 0, 0), (0, 0), value)]


def _base_input(outdir, coulomb_inter=None, **eli):
    """Small single-band static Eliashberg input (4x4, chi0q computed
    internally) built on the existing ``tests/rpa/input`` fixture."""
    interaction = {"path_to_input": RPA_INPUT,
                   "Geometry": "geom.dat",
                   "Transfer": "transfer.dat",
                   "CoulombIntra": "coulombintra.dat"}
    if coulomb_inter is not None:
        # absolute paths survive os.path.join(path_to_input, ...)
        interaction["CoulombInter"] = coulomb_inter
    eli_param = {"chi0q_mode": "calc", "solver_mode": "iteration",
                 "max_iter": 40, "convergence_tol": 1e-8}
    eli_param.update(eli)
    return {
        "mode": {"param": {"T": 1.0, "CellShape": [4, 4, 1],
                           "SubShape": [1, 1, 1], "Nmat": 128,
                           "filling": 0.5}},
        "file": {"input": {"interaction": interaction},
                 "output": {"path_to_output": outdir}},
        "eliashberg": eli_param,
    }


_CHI0Q_CACHE = {}


def _deterministic_input(outdir, **eli):
    """Same fixture as ``_base_input`` but with chi0q *loaded* from a file
    computed once.

    ``_calc_chi0q_internal`` (the RPA path) is not bit-reproducible run to run
    -- repeated calls with identical inputs differ at the ~1e-14 level (ULP
    jitter in the FFT/BLAS kernels, pre-existing and unrelated to this work).
    Everything downstream of chi0q *is* bit-reproducible, so a bit-identity
    regression test must pin chi0q by loading the same array both times.
    """
    inp = _base_input(outdir, chi0q_mode="load", **eli)
    if "chi0q" not in _CHI0Q_CACHE:
        _CHI0Q_CACHE["chi0q"] = sc._calc_chi0q_internal(
            _base_input(outdir, chi0q_mode="calc"))
    os.makedirs(outdir, exist_ok=True)
    np.savez(os.path.join(outdir, "chi0q.npz"), chi0q=_CHI0Q_CACHE["chi0q"])
    return inp


def _two_orbital_input(outdir, coulomb_inter, **eli):
    interaction = {"path_to_input": RPA_INPUT,
                   "Geometry": "geom_2orb.dat",
                   "Transfer": "transfer_nonso_2orb.dat",
                   "CoulombIntra": "coulombintra_2orb.dat",
                   "CoulombInter": coulomb_inter}
    eli_param = {"chi0q_mode": "calc", "solver_mode": "iteration",
                 "max_iter": 5}
    eli_param.update(eli)
    return {
        "mode": {"param": {"T": 1.0, "CellShape": [4, 4, 1],
                           "SubShape": [1, 1, 1], "Nmat": 64,
                           "filling": 0.5}},
        "file": {"input": {"interaction": interaction},
                 "output": {"path_to_output": outdir}},
        "eliashberg": eli_param,
    }


# ``_caplog`` and the assertApprox-providing base class now live in
# ``tests.approx_util`` (shared with the other bond/spectral test modules);
# ``_ApproxTestCase`` is kept as a local alias since every test class in this
# module still refers to it by that name.
_ApproxTestCase = ApproxTestCase


def _read_outputs(outdir):
    gap = np.loadtxt(os.path.join(outdir, "gap.dat"))
    with open(os.path.join(outdir, "eigenvalue.dat"), "rb") as f:
        eig_bytes = f.read()
    with open(os.path.join(outdir, "gap.dat"), "rb") as f:
        gap_bytes = f.read()
    return gap, eig_bytes, gap_bytes


def _leading_eigenvalue(outdir):
    with open(os.path.join(outdir, "eigenvalue.dat")) as f:
        lines = [l for l in f if not l.strip().startswith("#") and l.strip()]
    return float(lines[0].split()[0])


# ---------------------------------------------------------------------------
# S7.3 regression: the default path is untouched (bit-identical)
# ---------------------------------------------------------------------------

class TestDefaultPathRegression(_ApproxTestCase):

    def setUp(self):
        self.tmp_path = tempfile.mkdtemp(prefix="hwave_sc_bond_")
        self.addCleanup(shutil.rmtree, self.tmp_path, ignore_errors=True)

    def test_bond_channels_absent_and_false_are_bit_identical(self):
        """The flag defaults to false, and setting it explicitly false must run the
        *same* code path -- byte-for-byte identical gap.dat/eigenvalue.dat."""
        d1 = os.path.join(self.tmp_path, "a")
        d2 = os.path.join(self.tmp_path, "b")
        sc.calc_eliashberg(_deterministic_input(d1))
        sc.calc_eliashberg(_deterministic_input(d2, bond_channels=False))

        gap1, eig1, raw1 = _read_outputs(d1)
        gap2, eig2, raw2 = _read_outputs(d2)
        self.assertTrue(np.array_equal(gap1, gap2))
        self.assertEqual(eig1, eig2)
        self.assertEqual(raw1, raw2)

    def test_bond_options_ignored_with_warning_when_flag_off(self):
        """bond_max_shells / bond_memory_cap_gb set while bond_channels=false are
        ignored WITH A WARNING (spec S5) and change nothing bit-wise."""
        d1 = os.path.join(self.tmp_path, "a")
        d2 = os.path.join(self.tmp_path, "b")
        sc.calc_eliashberg(_deterministic_input(d1))
        with _caplog(logging.WARNING) as records:
            sc.calc_eliashberg(_deterministic_input(d2, bond_channels=False,
                                                    bond_max_shells=1,
                                                    bond_memory_cap_gb=2.0))
        msgs = " ".join(r.getMessage() for r in records)
        self.assertIn("bond_max_shells", msgs)
        self.assertIn("bond_memory_cap_gb", msgs)
        self.assertIn("bond_channels", msgs)

        _, eig1, raw1 = _read_outputs(d1)
        _, eig2, raw2 = _read_outputs(d2)
        self.assertTrue(eig1 == eig2 and raw1 == raw2)

    def test_compute_vertices_general_dressing_is_bit_identical(self):
        """The Task-4 refactor factors the inline dressing of
        ``_compute_vertices_general`` into a shared helper. Its output must stay
        bit-identical to the historical inline code (verbatim copy below)."""
        rng = np.random.default_rng(5)
        norb, Nx, Ny, Nz, nmat = 2, 3, 2, 1, 4
        nd = norb * norb
        chi0q = (rng.normal(size=(norb, norb, norb, norb, Nx, Ny, Nz, nmat))
                 + 1j * rng.normal(size=(norb, norb, norb, norb, Nx, Ny, Nz, nmat)))
        inter_k = {
            "CoulombIntra": rng.normal(size=(norb, norb, Nx, Ny, Nz)).astype(complex),
            "CoulombInter": rng.normal(size=(norb, norb, Nx, Ny, Nz)).astype(complex),
            "Hund": rng.normal(size=(norb, norb, Nx, Ny, Nz)).astype(complex),
        }
        got = sc._compute_vertices_general(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                                           pairing_type="singlet")

        # --- historical inline body (sc.py:1095-1118 before the refactor) ---
        S_all, C_all = sc._build_sc_matrices_all_q(inter_k, norb, Nx, Ny, Nz)
        si = nmat // 2
        chi0_static = chi0q[:, :, :, :, :, :, :, si].reshape(
            nd, nd, Nx, Ny, Nz).transpose(2, 3, 4, 0, 1).copy()
        I_mat = np.broadcast_to(np.eye(nd, dtype=complex),
                                (Nx, Ny, Nz, nd, nd)).copy()
        mat_s = I_mat - chi0_static @ S_all
        mat_c = I_mat + chi0_static @ C_all
        chis = np.linalg.solve(mat_s, chi0_static)
        chic = np.linalg.solve(mat_c, chi0_static)
        Vs_all = (1.5 * (S_all @ chis @ S_all) - 0.5 * (C_all @ chic @ C_all)
                  + 0.5 * (S_all + C_all))
        ref = Vs_all.reshape(Nx, Ny, Nz, norb, norb, norb, norb).transpose(
            3, 4, 5, 6, 0, 1, 2)

        self.assertTrue(np.array_equal(got, ref))


def _ill_conditioned_general_case(target=1.0e-3):
    """Build a DEFAULT-path ``_compute_vertices_general`` input whose spin RPA
    denominator ``I - chi0 @ S`` is invertible but badly conditioned
    (``sigma_min/sigma_max <= 1e-3``).

    ``chi0q`` is scaled so that the smallest eigenvalue of ``I - chi0 @ S`` at
    one q-point is ``target``: the solve is perfectly well defined (the
    dressed chi comes back of order 1e2) but the matrix sits in the region the
    BOND guard refuses.
    """
    rng = np.random.default_rng(11)
    norb, Nx, Ny, Nz, nmat = 2, 2, 1, 1, 2
    nd = norb * norb
    chi0q = rng.normal(
        size=(norb, norb, norb, norb, Nx, Ny, Nz, nmat)).astype(complex)
    inter_k = {
        "CoulombIntra": rng.normal(size=(norb, norb, Nx, Ny, Nz)).astype(complex),
        "CoulombInter": rng.normal(size=(norb, norb, Nx, Ny, Nz)).astype(complex),
        "Hund": rng.normal(size=(norb, norb, Nx, Ny, Nz)).astype(complex),
    }
    S_all, C_all = sc._build_sc_matrices_all_q(inter_k, norb, Nx, Ny, Nz)
    si = nmat // 2

    def _static(arr):
        return arr[:, :, :, :, :, :, :, si].reshape(
            nd, nd, Nx, Ny, Nz).transpose(2, 3, 4, 0, 1).copy()

    eig = np.linalg.eigvals(_static(chi0q)[0, 0, 0] @ S_all[0, 0, 0])
    lam = float(sorted(x.real for x in eig if abs(x.imag) < 1e-12)[-1])
    chi0q = chi0q * ((1.0 - target) / lam)
    return chi0q, inter_k, norb, Nx, Ny, Nz, nmat, _static(chi0q), S_all, C_all


class TestBackwardCompatibilityOfConditioningGuard(_ApproxTestCase):

    def test_general_vertex_ill_conditioned_default_path_still_runs(self):
        """BACKWARD COMPATIBILITY: the bond conditioning guard must NOT leak onto
        the legacy (default, non-bond) ``_compute_vertices_general`` path.

        Before 712b1a0 that path simply called ``np.linalg.solve``; an invertible
        but poorly conditioned RPA denominator returned a (large but finite)
        result. Sharing ``dress_bond`` must not turn such a configuration into a
        ``ValueError``.
        """
        (chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
         chi0_static, S_all, C_all) = _ill_conditioned_general_case()

        # the case really is inside the bond guard's refusal region ...
        I_mat = np.broadcast_to(np.eye(norb * norb, dtype=complex),
                                (Nx, Ny, Nz, norb * norb, norb * norb)).copy()
        mat_s = I_mat - chi0_static @ S_all
        sv = np.linalg.svd(mat_s.reshape(-1, norb * norb, norb * norb),
                           compute_uv=False)
        self.assertLessEqual(float((sv[:, -1] / sv[:, 0]).min()), 1.0e-3)
        # ... and still perfectly invertible (no vanishing singular value)
        self.assertGreater(float(sv[:, -1].min()), 1.0e-8)

        got = sc._compute_vertices_general(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                                           pairing_type="singlet")

        # historical inline computation (712b1a0 and earlier)
        mat_c = I_mat + chi0_static @ C_all
        chis = np.linalg.solve(mat_s, chi0_static)
        chic = np.linalg.solve(mat_c, chi0_static)
        Vs_all = (1.5 * (S_all @ chis @ S_all) - 0.5 * (C_all @ chic @ C_all)
                  + 0.5 * (S_all + C_all))
        ref = Vs_all.reshape(Nx, Ny, Nz, norb, norb, norb, norb).transpose(
            3, 4, 5, 6, 0, 1, 2)
        self.assertTrue(np.all(np.isfinite(got)))
        self.assertTrue(np.array_equal(got, ref))

    def test_ill_conditioned_case_is_still_refused_on_the_bond_path(self):
        """The other half of the contract: the guard stays where it belongs. The
        very same conditioning must still raise through the opt-in bond entry
        point (default ``cond_tol``)."""
        from hwave.solver import bond_channels

        (_, _, _, _, _, _, _, chi0_static, S_all,
         C_all) = _ill_conditioned_general_case()
        with self.assertRaisesRegex(ValueError, r"(?i)singular|cond_tol"):
            bond_channels.dress_bond(chi0_static, S_all, C_all)
        # and cond_tol=None genuinely disables it (what the legacy path passes)
        chis, chic = bond_channels.dress_bond(chi0_static, S_all, C_all,
                                              cond_tol=None)
        self.assertTrue(np.all(np.isfinite(chis)) and np.all(np.isfinite(chic)))


# ---------------------------------------------------------------------------
# S5 guards
# ---------------------------------------------------------------------------

class TestS5Guards(_ApproxTestCase):

    def setUp(self):
        self.tmp_path = tempfile.mkdtemp(prefix="hwave_sc_bond_")
        self.addCleanup(shutil.rmtree, self.tmp_path, ignore_errors=True)

    def test_guard_flex_mode(self):
        """(a) bond_channels=true + chi0q_mode='flex' -> error naming the deferred
        bond-FLEX work."""
        ci = os.path.join(self.tmp_path, "ci.dat")
        _write_w90(ci, 1, _nn_entries(1.0))
        inp = _base_input(self.tmp_path, coulomb_inter=ci, bond_channels=True,
                          chi0q_mode="flex")
        with self.assertRaisesRegex(ValueError, r"(?i)flex"):
            sc.calc_eliashberg(inp)

    def test_guard_multi_orbital(self):
        """(b) bond_channels=true + norb>1 -> error at the top level only."""
        ci = os.path.join(self.tmp_path, "ci2.dat")
        _write_w90(ci, 2, [((1, 0, 0), (0, 1), 0.25), ((-1, 0, 0), (1, 0), 0.25)])
        inp = _two_orbital_input(self.tmp_path, ci, bond_channels=True)
        with self.assertRaisesRegex(ValueError, r"(?i)multi-orbital|norb"):
            sc.calc_eliashberg(inp)

    def test_guard_complex_coulomb_inter(self):
        """(c) bond_channels=true + a non-real CoulombInter entry -> error."""
        ci = os.path.join(self.tmp_path, "ci_cplx.dat")
        _write_w90(ci, 1, [((1, 0, 0), (0, 0), 1.0 + 0.3j),
                           ((-1, 0, 0), (0, 0), 1.0 - 0.3j)])
        inp = _base_input(self.tmp_path, coulomb_inter=ci, bond_channels=True)
        with self.assertRaisesRegex(ValueError, r"(?i)complex|real"):
            sc.calc_eliashberg(inp)

    def test_guard_no_coulomb_inter(self):
        """(d) bond_channels=true with no CoulombInter term at all -> error."""
        inp = _base_input(self.tmp_path, bond_channels=True)
        with self.assertRaisesRegex(ValueError, r"(?i)CoulombInter"):
            sc.calc_eliashberg(inp)

    def test_guard_dynamic_entry_point(self):
        """(e) bond_channels=true on the dynamic entry point -> error, both via
        calc_eliashberg's dispatch and via solve_dynamic directly."""
        ci = os.path.join(self.tmp_path, "ci.dat")
        _write_w90(ci, 1, _nn_entries(1.0))
        inp = _base_input(self.tmp_path, coulomb_inter=ci, bond_channels=True,
                          frequency="dynamic", chi0q_mode="flex")
        with self.assertRaisesRegex(ValueError, r"(?i)dynamic"):
            sc.calc_eliashberg(inp)

        from hwave.solver import eliashberg_dynamic
        with self.assertRaisesRegex(ValueError, r"(?i)dynamic"):
            eliashberg_dynamic.solve_dynamic(inp)

    def test_resource_preflight_errors_and_names_channels(self):
        """S3.2: the preflight refuses to allocate beyond bond_memory_cap_gb and
        names the offending channel count."""
        ci = os.path.join(self.tmp_path, "ci.dat")
        _write_w90(ci, 1, _nn_entries(1.0))
        inp = _base_input(self.tmp_path, coulomb_inter=ci, bond_channels=True,
                          bond_memory_cap_gb=1.0e-9)
        with self.assertRaisesRegex(ValueError, r"(?i)bond_memory_cap_gb"):
            sc.calc_eliashberg(inp)

    def test_resource_preflight_raises_before_green_allocation(self):
        """I-2 review fix: the resource preflight (S3.2) must run BEFORE
        ``green_kw`` is allocated -- previously ``_calc_green`` ran first in the
        bond-channel dispatch branch, so a large-Nmat run could OOM before the
        cap was ever consulted, defeating "never a silent runaway allocation".
        Patch ``_calc_green`` to blow up if called; the preflight's ValueError
        must surface instead of the patched crash."""
        ci = os.path.join(self.tmp_path, "ci.dat")
        _write_w90(ci, 1, _nn_entries(1.0))
        inp = _base_input(self.tmp_path, coulomb_inter=ci, bond_channels=True,
                          bond_memory_cap_gb=1.0e-9)

        def _boom(*args, **kwargs):
            raise AssertionError(
                "_calc_green must not be called before the bond-channel "
                "resource preflight (I-2)")

        orig_calc_green = sc._calc_green
        self.addCleanup(setattr, sc, "_calc_green", orig_calc_green)
        sc._calc_green = _boom
        with self.assertRaisesRegex(ValueError, r"(?i)bond_memory_cap_gb"):
            sc.calc_eliashberg(inp)


def _measure_peak_bytes(fn):
    """Peak heap bytes (numpy arrays included) allocated while ``fn`` runs.

    NumPy registers its data allocations with ``tracemalloc`` (the
    ``numpy.core.multiarray`` domain), so ``get_traced_memory()`` sees the
    big buffers.  Returns ``None`` when that tracking is unavailable, so the
    caller can skip rather than assert on a meaningless zero.
    """
    import tracemalloc

    tracemalloc.start()
    try:
        tracemalloc.reset_peak()
        probe = np.empty(1 << 20, dtype=np.complex128)   # 16 MB
        _, probe_peak = tracemalloc.get_traced_memory()
        del probe
        if probe_peak < (1 << 20) * 16 * 0.9:
            return None                                   # numpy not tracked
        tracemalloc.reset_peak()
        fn()
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    return peak


class TestPreflightMemoryBudget(_ApproxTestCase):

    def test_preflight_bubble_budget_covers_the_measured_bond_bubble_peak(self):
        """The preflight must not UNDERCOUNT ``bond_bubble``'s real high-water
        mark: a run it estimates to fit under ``bond_memory_cap_gb`` would then
        OOM anyway, defeating the whole point of the guard.

        Measure the actual peak of a small ``bond_bubble`` call and require the
        preflight's own bubble+green+chi_bar budget for exactly that case to be
        at least as large.  If ``bond_bubble`` grows a buffer without the
        preflight's documented buffer list growing with it, this fails.
        """
        from hwave.solver import bond_channels as bc

        norb, Nx, Ny, Nz, nmat, beta = 2, 4, 4, 4, 32, 10.0
        bond_set = bc.resolve_interactions(
            {((1, 0, 0), (0, 0)): 0.25, ((-1, 0, 0), (0, 0)): 0.25},
            [(0, 0, 0), (1, 0, 0), (-1, 0, 0)], norb, bond_max_shells=1)

        rng = np.random.default_rng(0)
        shape = (norb, norb, Nx, Ny, Nz, nmat)
        green = (rng.normal(size=shape) + 1j * rng.normal(size=shape))

        bc.bond_bubble(green, bond_set, beta)     # warm up numpy's FFT caches
        peak = _measure_peak_bytes(lambda: bc.bond_bubble(green, bond_set, beta))
        if peak is None:
            self.skipTest("tracemalloc does not track numpy array data here")

        est = sc._bond_memory_estimate(norb, bond_set, Nx, Ny, Nz, nmat)
        # `green` itself is allocated outside the measured region and is budgeted
        # separately (est["green_bytes"]), so the budget bond_bubble's own peak
        # must fit into is its working set plus the chi_bar it returns.
        budget = est["bubble_bytes"] + est["chi_bar_bytes"]
        self.assertGreaterEqual(budget, peak, (
            "the preflight budgets {:.3f} MB for bond_bubble's work buffers but "
            "the measured peak is {:.3f} MB -- the buffer list in "
            "_bond_memory_estimate has desynced from bond_bubble".format(
                budget / 1e6, peak / 1e6)))
        # ... and it must not be wildly conservative either, or the cap becomes
        # useless in the other direction.
        self.assertLessEqual(budget, 3.0 * peak)


def _large_B_bond_set(norb, B):
    """A reversal-closed ``ResolvedInteractionSet`` with ``B`` channels.

    Built directly (bypassing ``resolve_interactions``' shell-cutoff/topology
    logic) purely to exercise ``bare_bond_vertices`` at a large ``ND = nd*B``
    -- the function only reads ``delta_r``/``v_bond``/``reverse``/
    ``n_channels``, so a synthetic 1-D chain of channel pairs ``(+i,0,0)``/
    ``(-i,0,0)`` satisfies its reversal-closure guard without needing a real
    interaction fixture.
    """
    from hwave.solver import bond_channels as bc

    rng = np.random.default_rng(1)
    delta_r = [(0, 0, 0)]
    v_bond = [np.zeros((norb, norb), dtype=complex)]
    reverse = [0]
    n_pairs = (B - 1) // 2
    for i in range(1, n_pairs + 1):
        v = (rng.normal(size=(norb, norb))
             + 1j * rng.normal(size=(norb, norb)))
        pos = len(delta_r)
        delta_r.append((i, 0, 0))
        v_bond.append(v)
        neg = len(delta_r)
        delta_r.append((-i, 0, 0))
        v_bond.append(v.conj().T)
        reverse.append(neg)
        reverse.append(pos)
    return bc.ResolvedInteractionSet(
        delta_r=tuple(delta_r), v_bond=tuple(v_bond),
        reverse=tuple(reverse), n_channels=len(delta_r))


class TestPreflightVertexBudget(_ApproxTestCase):

    def test_preflight_vertex_budget_covers_the_measured_bare_bond_vertices_peak(self):
        """The preflight must not UNDERCOUNT ``bare_bond_vertices``'s real
        high-water mark either: the non-q-resolved ``ND x ND`` temporaries it
        keeps alive (``P``, ``D``, ``Dh``, ``Id``, ``Q_s``, ``Q_t``, ``B_s``,
        ``B_t``, ``Vpp_s``, ``Vpp_t``) were previously omitted from
        ``_bond_memory_estimate`` entirely, so at ``N_q = 1`` with a large ``B``
        a run the preflight waved through could still OOM building the vertices.

        ``bare_bond_vertices`` itself also allocates the q-resolved ``S_bond``/
        ``C_bond`` outputs (at ``N_q = 1`` the same per-array size as an ``ND x
        ND`` temporary); those two are already budgeted elsewhere as 2 of the 9
        ``sc._BOND_N_Q_ARRAYS`` (``est["chi_bar_bytes"]`` is the per-array size),
        so the budget a measurement of THIS call alone must fit into is
        ``vertex_bytes + 2 * chi_bar_bytes`` -- mirroring how the ``bond_bubble``
        peak test above adds ``chi_bar_bytes`` to ``bubble_bytes`` for the same
        reason.

        Mirrors ``test_preflight_bubble_budget_covers_the_measured_bond_bubble_peak``
        for the vertex-construction half of the pipeline.
        """
        from hwave.solver import bond_channels as bc

        norb = 1
        Nx, Ny, Nz = 1, 1, 1
        nmat = 4  # unused by bare_bond_vertices; only shapes _bond_memory_estimate
        B = 1001  # large enough that fixed per-array overhead is negligible

        bond_set = _large_B_bond_set(norb, B)
        nd = norb * norb
        S0_q = np.zeros((Nx, Ny, Nz, nd, nd), dtype=complex)
        C0_q = np.zeros((Nx, Ny, Nz, nd, nd), dtype=complex)
        S0_q[0, 0, 0] = 1.0
        C0_q[0, 0, 0] = 0.5

        bc.bare_bond_vertices(bond_set, S0_q, C0_q, norb)  # warm up

        peak = _measure_peak_bytes(
            lambda: bc.bare_bond_vertices(bond_set, S0_q, C0_q, norb))
        if peak is None:
            self.skipTest("tracemalloc does not track numpy array data here")

        est = sc._bond_memory_estimate(norb, bond_set, Nx, Ny, Nz, nmat)
        budget = est["vertex_bytes"] + 2 * est["chi_bar_bytes"]
        # tracemalloc/numpy also attribute a small, ND-independent amount of
        # measurement-harness bookkeeping (observed up to ~70 KB here, e.g. from
        # ``_measure_peak_bytes``'s own probe allocation and numpy's internal
        # ufunc buffer cache) that is not one of the counted ND x ND buffers; a
        # fixed 1 MB slack comfortably absorbs that noise without masking an
        # actual missing multi-MB buffer.
        overhead_slack = 1 << 20  # 1 MB
        self.assertGreaterEqual(budget + overhead_slack, peak, (
            "the preflight budgets {:.3f} MB for bare_bond_vertices's non-q-"
            "resolved ND x ND temporaries (plus its S_bond/C_bond outputs) but "
            "the measured peak is {:.3f} MB -- _bond_memory_estimate has "
            "desynced from bare_bond_vertices".format(
                budget / 1e6, peak / 1e6)))
        # ... and it must not be wildly conservative either, or the cap becomes
        # useless in the other direction.
        self.assertLessEqual(budget, 3.0 * peak)


# ---------------------------------------------------------------------------
# S4.3 star -- the Case-2 correction, and the V=0 reduction
# ---------------------------------------------------------------------------

class TestCase2CorrectionAndVZeroReduction(_ApproxTestCase):

    def setUp(self):
        self.tmp_path = tempfile.mkdtemp(prefix="hwave_sc_bond_")
        self.addCleanup(shutil.rmtree, self.tmp_path, ignore_errors=True)

    def test_m0_blocks_single_band_match_simple_path(self):
        """Single band: the Case-2-corrected m=0 blocks must be S=U and
        C=U+2V(q) -- exactly the ``_compute_vertices_simple`` Ws=-U/Wc=U+2V pair
        (spec S4.5 "Concretely, from the derivation")."""
        from hwave.solver.bond_channels import resolve_interactions

        norb, Nx, Ny, Nz = 1, 4, 4, 1
        U, V = 4.0, 0.7
        kx = np.linspace(0, 2 * np.pi, Nx, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, Ny, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, Nz, endpoint=False)
        interactions = {
            "CoulombIntra": {((0, 0, 0), (0, 0)): U},
            "CoulombInter": {ir: V for ir in
                             [((-1, 0, 0), (0, 0)), ((1, 0, 0), (0, 0)),
                              ((0, -1, 0), (0, 0)), ((0, 1, 0), (0, 0))]},
        }
        inter_k = sc._build_interaction_k(kx, ky, kz, interactions, norb)
        bond_set = resolve_interactions(interactions["CoulombInter"], np.eye(3), norb)
        S0, C0 = sc._build_bond_m0_blocks(bond_set, interactions, inter_k, norb,
                                          kx, ky, kz)
        self.assertEqual(S0.shape, (Nx, Ny, Nz, 1, 1))
        Vq = inter_k["CoulombInter"][0, 0]
        np.testing.assert_allclose(S0[:, :, :, 0, 0], U, atol=1e-13)
        np.testing.assert_allclose(C0[:, :, :, 0, 0], U + 2.0 * Vq, atol=1e-13)

    def test_case2_correction_keeps_only_rzero_fock(self):
        """S4.3 star (multi-orbital, exercised below the top-level norb guard):
        the m=0 Fock (ab,ab) element carries ONLY V_ab(R=0), while the charge
        (aa,bb) Hartree keeps the full 2 V_ab(q)."""
        from hwave.solver.bond_channels import resolve_interactions

        norb, Nx, Ny, Nz = 2, 4, 1, 1
        Uprime, Vbond = 1.0, 0.25
        kx = np.linspace(0, 2 * np.pi, Nx, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, Ny, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, Nz, endpoint=False)
        interactions = {
            "CoulombInter": {
                ((0, 0, 0), (0, 1)): Uprime,
                ((0, 0, 0), (1, 0)): Uprime,
                ((1, 0, 0), (0, 1)): Vbond,
                ((-1, 0, 0), (1, 0)): Vbond,
            },
        }
        inter_k = sc._build_interaction_k(kx, ky, kz, interactions, norb)
        bond_set = resolve_interactions(interactions["CoulombInter"], np.eye(3), norb)
        S0, C0 = sc._build_bond_m0_blocks(bond_set, interactions, inter_k, norb,
                                          kx, ky, kz)

        i01 = 0 * norb + 1                     # (l1,l2) = (0,1)
        # Fock: ONLY the R=0 component (U'), NOT the full V_01(q).
        np.testing.assert_allclose(S0[:, :, :, i01, i01], Uprime, atol=1e-13)
        np.testing.assert_allclose(C0[:, :, :, i01, i01], -Uprime, atol=1e-13)
        # and it is genuinely a correction: the un-corrected value is q-dependent
        V01_q = inter_k["CoulombInter"][0, 1]
        self.assertFalse(np.allclose(V01_q, Uprime))
        # Hartree (aa,bb) keeps the full q-dependent 2 V_ab(q).
        i00, i11 = 0, norb + 1
        np.testing.assert_allclose(C0[:, :, :, i00, i11], 2.0 * V01_q, atol=1e-13)

    def test_bond_kernel_reduces_to_standard_kernel_when_V_is_zero(self):
        """S7.4: with a DECLARED but zero NN topology (B=5) the whole bond path
        (m=0 blocks -> dressing -> two-part operator) must reproduce the standard
        single-band kernel exactly, when both are fed the SAME bubble."""
        from hwave.solver.bond_channels import (resolve_interactions, bond_bubble,
                                                bare_bond_vertices, dress_bond,
                                                make_bond_kernel)

        norb, Nx, Ny, Nz, nmat = 1, 4, 4, 1, 64
        beta = 2.0
        U = 3.0
        kx = np.linspace(0, 2 * np.pi, Nx, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, Ny, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, Nz, endpoint=False)
        interactions = {
            "CoulombIntra": {((0, 0, 0), (0, 0)): U},
            "CoulombInter": {ir: 0.0 for ir in
                             [((-1, 0, 0), (0, 0)), ((1, 0, 0), (0, 0)),
                              ((0, -1, 0), (0, 0)), ((0, 1, 0), (0, 0))]},
        }
        inter_k = sc._build_interaction_k(kx, ky, kz, interactions, norb)
        hr = {((1, 0, 0), (0, 0)): 1.0, ((-1, 0, 0), (0, 0)): 1.0,
              ((0, 1, 0), (0, 0)): 1.0, ((0, -1, 0), (0, 0)): 1.0}
        eps = sc._build_hamiltonian_k(kx, ky, kz, hr, norb)
        evals, evecs = sc._calc_eigenvalues(eps)
        green = sc._calc_green(evals, evecs, 0.1, beta, nmat)

        bond_set = resolve_interactions(interactions["CoulombInter"], np.eye(3), norb)
        self.assertEqual(bond_set.n_channels, 5)  # declared topology survives V=0

        chi_bar = bond_bubble(green, bond_set, beta)
        S0, C0 = sc._build_bond_m0_blocks(bond_set, interactions, inter_k, norb,
                                          kx, ky, kz)
        S_bond, C_bond, Vpp_s, Vpp_t = bare_bond_vertices(bond_set, S0, C0, norb)
        chi_s, chi_c = dress_bond(chi_bar, S_bond, C_bond)
        A_bond, n_bond = make_bond_kernel(chi_s, chi_c, S_bond, C_bond,
                                          Vpp_s, Vpp_t, green, bond_set,
                                          "singlet", beta, g2_tail=True)

        # standard path fed the SAME bubble (the m=m'=0 block of chi_bar)
        chi0q = chi_bar[:, :, :, 0, 0].transpose(0, 1, 2).reshape(
            1, 1, Nx, Ny, Nz, 1)
        Pc, Ps = sc._compute_vertices_simple(chi0q, inter_k, norb, Nx, Ny, Nz, 1,
                                             pairing_type="singlet",
                                             static_index=0)
        G2 = sc._calc_g2(green, beta)
        A_std, n_std = sc._make_kernel_operator(Pc + Ps, G2, norb, Nx, Ny, Nz)
        self.assertEqual(n_bond, n_std)

        rng = np.random.default_rng(3)
        v = rng.normal(size=n_std) + 1j * rng.normal(size=n_std)
        np.testing.assert_allclose(A_bond.matvec(v), A_std.matvec(v),
                                   rtol=1e-9, atol=1e-11)

    def test_bond_channels_true_with_zero_V_matches_default_run(self):
        """End-to-end S7.4: bond_channels=true with a declared-but-zero V (B=5
        declared topology, zero Fock/Cooper) is numerically equivalent to
        bond_channels=false. Numerical (not bit-wise) per the spec's tolerance
        policy -- the enabled path is a different computation of the same
        quantity."""
        ci = os.path.join(self.tmp_path, "ci0.dat")
        _write_w90(ci, 1, _nn_entries(0.0))
        d_off = os.path.join(self.tmp_path, "off")
        d_on = os.path.join(self.tmp_path, "on")
        sc.calc_eliashberg(_base_input(d_off, coulomb_inter=ci,
                                       bond_channels=False))
        sc.calc_eliashberg(_base_input(d_on, coulomb_inter=ci,
                                       bond_channels=True))
        lam_off = _leading_eigenvalue(d_off)
        lam_on = _leading_eigenvalue(d_on)
        self.assertApprox(lam_off, lam_on, rel=1.0e-7)

    def test_bond_channels_true_is_not_inert_for_nonzero_V(self):
        """The flag must actually route through the bond path: at V != 0 the
        bond-resolved result differs from the default (which, for a single band
        with a 4-index chi0q, carries no inter-site V at all -- the scope
        limitation this feature addresses)."""
        ci = os.path.join(self.tmp_path, "ci1.dat")
        _write_w90(ci, 1, _nn_entries(1.0))
        d_off = os.path.join(self.tmp_path, "off")
        d_on = os.path.join(self.tmp_path, "on")
        sc.calc_eliashberg(_base_input(d_off, coulomb_inter=ci,
                                       bond_channels=False))
        sc.calc_eliashberg(_base_input(d_on, coulomb_inter=ci,
                                       bond_channels=True))
        self.assertGreater(
            abs(_leading_eigenvalue(d_on) - _leading_eigenvalue(d_off)), 1e-3)

    def test_bond_run_writes_provenance_without_breaking_readers(self):
        """S5 output provenance is additive: new comment keys in eigenvalue.dat,
        existing numeric rows (and the tsweep reader) unaffected."""
        from hwave import tsweep

        ci = os.path.join(self.tmp_path, "ci.dat")
        _write_w90(ci, 1, _nn_entries(0.5))
        d_on = os.path.join(self.tmp_path, "on")
        sc.calc_eliashberg(_base_input(d_on, coulomb_inter=ci, bond_channels=True,
                                       solver_mode="both", num_eigenvalues=3))
        text = open(os.path.join(d_on, "eigenvalue.dat")).read()
        self.assertIn("bond_channels", text)
        self.assertIn("bond_max_shells", text)
        self.assertTrue("delta_r" in text or "channels" in text)
        self.assertIn("RPA-ladder", text)
        # a downstream reader still parses the numeric rows
        re, im, match = tsweep.parse_leading_eig(
            os.path.join(d_on, "eigenvalue.dat"))
        self.assertTrue(np.isfinite(re) and np.isfinite(im))

    def test_collapsed_to_pure_hubbard_flag(self):
        """A B=1 run (bond_max_shells=0 with no inter-site V declared beyond
        on-site) is recorded as collapsed_to_pure_hubbard so a single-channel run
        is distinguishable from a genuine multi-channel one (spec S5)."""
        ci = os.path.join(self.tmp_path, "ci_onsite.dat")
        # on-site-only CoulombInter: no bond channels at all -> B = 1
        _write_w90(ci, 1, [((0, 0, 0), (0, 0), 0.0)])
        d_on = os.path.join(self.tmp_path, "on")
        sc.calc_eliashberg(_base_input(d_on, coulomb_inter=ci, bond_channels=True))
        text = open(os.path.join(d_on, "eigenvalue.dat")).read()
        self.assertIn("collapsed_to_pure_hubbard", text)


# ---------------------------------------------------------------------------
# Task 7 -- S4.5 lambda attribution on the symmetrized kernel, the runtime
# Hermiticity/positivity preconditions, the odd f-wave seed, and the S7.7
# invariant-subspace tracker.
# ---------------------------------------------------------------------------

def _analytic_pp_fixture():
    """Spec S7.9: single-band 1D 4-point grid, U=4, V(R)=1, chi_s=chi_c=0.

    green == 1 with nmat=1 and beta=1 gives G2(k) = 1 for every k, so the pair
    weight is w = 1 and the (T/N) G G fold is 1/4.  Returns everything
    make_bond_kernel_parts needs plus the odd gap Delta = sin k.
    """
    from hwave.solver.bond_channels import (resolve_interactions,
                                            bare_bond_vertices)
    norb, Nx, Ny, Nz, nmat = 1, 4, 1, 1, 1
    beta = 1.0
    green = np.ones((norb, norb, Nx, Ny, Nz, nmat), dtype=complex)
    ci = {((1, 0, 0), (0, 0)): 1.0, ((-1, 0, 0), (0, 0)): 1.0}
    bset = resolve_interactions(ci, np.eye(3), norb=norb)
    U = 4.0
    S0 = np.full((1, 1, 1, 1, 1), U, dtype=complex)
    C0 = np.full((1, 1, 1, 1, 1), 8.0, dtype=complex)     # U + 2 V(q=0)
    S_bond, C_bond, Vpp_s, Vpp_t = bare_bond_vertices(bset, S0, C0, norb)
    ND = S_bond.shape[-1]
    chi_zero = np.zeros((Nx, Ny, Nz, ND, ND), dtype=complex)
    kx = 2.0 * np.pi * np.arange(Nx) / Nx
    phi = np.sin(kx).reshape(norb, norb, Nx, Ny, Nz).astype(complex)
    return dict(green=green, beta=beta, bond_set=bset, S_bond=S_bond,
                C_bond=C_bond, Vpp_s=Vpp_s, Vpp_t=Vpp_t, chi=chi_zero,
                phi=phi, norb=norb, shape=(Nx, Ny, Nz))


def _physical_single_band(V=1.0, U=3.0, Nx=4, Ny=4, Nz=1, nmat=64, beta=2.0):
    """A physical (centrosymmetric, spin-degenerate) single-band bond setup:
    real >= 0 pair weight and a Hermitian symmetrized kernel."""
    from hwave.solver.bond_channels import (resolve_interactions, bond_bubble,
                                            bare_bond_vertices, dress_bond)
    norb = 1
    kx = np.linspace(0, 2 * np.pi, Nx, endpoint=False)
    ky = np.linspace(0, 2 * np.pi, Ny, endpoint=False)
    kz = np.linspace(0, 2 * np.pi, Nz, endpoint=False)
    interactions = {
        "CoulombIntra": {((0, 0, 0), (0, 0)): U},
        "CoulombInter": {ir: V for ir in
                         [((-1, 0, 0), (0, 0)), ((1, 0, 0), (0, 0)),
                          ((0, -1, 0), (0, 0)), ((0, 1, 0), (0, 0))]},
    }
    inter_k = sc._build_interaction_k(kx, ky, kz, interactions, norb)
    hr = {((1, 0, 0), (0, 0)): 1.0, ((-1, 0, 0), (0, 0)): 1.0,
          ((0, 1, 0), (0, 0)): 1.0, ((0, -1, 0), (0, 0)): 1.0}
    eps = sc._build_hamiltonian_k(kx, ky, kz, hr, norb)
    evals, evecs = sc._calc_eigenvalues(eps)
    green = sc._calc_green(evals, evecs, 0.1, beta, nmat)
    bset = resolve_interactions(interactions["CoulombInter"], np.eye(3), norb)
    chi_bar = bond_bubble(green, bset, beta)
    S0, C0 = sc._build_bond_m0_blocks(bset, interactions, inter_k, norb,
                                      kx, ky, kz)
    S_bond, C_bond, Vpp_s, Vpp_t = bare_bond_vertices(bset, S0, C0, norb)
    chi_s, chi_c = dress_bond(chi_bar, S_bond, C_bond)
    return dict(green=green, beta=beta, bond_set=bset, S_bond=S_bond,
                C_bond=C_bond, Vpp_s=Vpp_s, Vpp_t=Vpp_t, chi_s=chi_s,
                chi_c=chi_c, norb=norb, shape=(Nx, Ny, Nz),
                kgrid=(kx, ky, kz))


def _gap_parity_residual(sigma):
    """max |Delta(-k) + Delta(k)| (0 for a perfectly ODD gap)."""
    flipped = np.roll(sigma[:, :, ::-1, ::-1, ::-1], (1, 1, 1), (2, 3, 4))
    return np.abs(flipped + sigma).max()


def _gap_even_residual(sigma):
    flipped = np.roll(sigma[:, :, ::-1, ::-1, ::-1], (1, 1, 1), (2, 3, 4))
    return np.abs(flipped - sigma).max()


# --- (a) lambda attribution on the analytic S7.9 fixture -------------------

class TestLambdaAttributionAndOddSeed(_ApproxTestCase):

    def setUp(self):
        self.tmp_path = tempfile.mkdtemp(prefix="hwave_sc_bond_")
        self.addCleanup(shutil.rmtree, self.tmp_path, ignore_errors=True)

    def test_lambda_attribution_analytic_fixture(self):
        """S7.9 + S4.5: lambda_t = -1 splits as lambda_t^pp = -1, lambda_t^fl = 0,
        as REAL Rayleigh quotients on the symmetrized kernel."""
        from hwave.solver import bond_channels as bc

        f = _analytic_pp_fixture()
        A, A_fl, A_pp, vec_size = bc.make_bond_kernel_parts(
            f["chi"], f["chi"], f["S_bond"], f["C_bond"], f["Vpp_s"], f["Vpp_t"],
            f["green"], f["bond_set"], "triplet", f["beta"])
        W = bc.pair_weight(f["green"], f["beta"])

        v = f["phi"].ravel()
        # the fixture's gap IS an eigenvector with lambda = -1
        np.testing.assert_allclose(A.matvec(v), -v, atol=1e-12)

        attr = bc.attribute_lambda(v, W, A_pp, A_fl, op_full=A)
        self.assertApprox(attr["lambda_pp"], -1.0, abs=1e-12)
        self.assertApprox(attr["lambda_fl"], 0.0, abs=1e-12)
        self.assertApprox(attr["lambda"], -1.0, abs=1e-12)
        # sum identity
        self.assertApprox(attr["lambda_pp"] + attr["lambda_fl"],
                          attr["lambda"], abs=1e-12)
        self.assertApprox(attr["sum_residual"], 0.0, abs=1e-12)
        # real Rayleigh quotients (the Hermitian path)
        self.assertIsInstance(attr["lambda_pp"], float)
        self.assertTrue(abs(attr["imag_pp"]) < 1e-12 and abs(attr["imag_fl"]) < 1e-12)

    def test_lambda_attribution_sum_identity_physical(self):
        """On a physical bond setup the two parts add up to the full Rayleigh
        quotient, and (for the converged eigenvector) to the eigenvalue."""
        from hwave.solver import bond_channels as bc
        from scipy.sparse.linalg import eigs

        f = _physical_single_band(V=1.0)
        A, A_fl, A_pp, n = bc.make_bond_kernel_parts(
            f["chi_s"], f["chi_c"], f["S_bond"], f["C_bond"], f["Vpp_s"],
            f["Vpp_t"], f["green"], f["bond_set"], "singlet", f["beta"])
        W = bc.pair_weight(f["green"], f["beta"])

        vals, vecs = eigs(A, k=1, which="LR", maxiter=5000, tol=1e-12)
        v = vecs[:, 0]
        attr = bc.attribute_lambda(v, W, A_pp, A_fl, op_full=A)
        # the Rayleigh quotient of the converged eigenvector IS the eigenvalue,
        # and it is real (the Hermitian path)
        self.assertApprox(attr["lambda"], vals[0].real, rel=1e-8)
        self.assertLess(abs(attr["imag"]), 1e-8 * max(1.0, abs(attr["lambda"])))
        self.assertApprox(attr["lambda_pp"] + attr["lambda_fl"],
                          attr["lambda"], abs=1e-10)

        # the sum identity is an operator identity, so it also holds off the
        # eigenvector -- checked on a d_x2y2 gap, for which BOTH parts are
        # non-zero (the leading eigenvector of this fixture happens to be
        # orthogonal to every bond harmonic, so its lambda^pp vanishes exactly).
        kx, ky, kz = f["kgrid"]
        d = sc._initialize_gap("d_x2y2", 1, kx, ky, kz).ravel().astype(complex)
        attr_d = bc.attribute_lambda(d, W, A_pp, A_fl, op_full=A)
        self.assertGreater(abs(attr_d["lambda_fl"]), 1e-6)
        self.assertGreater(abs(attr_d["lambda_pp"]), 1e-6)
        self.assertApprox(attr_d["lambda_pp"] + attr_d["lambda_fl"],
                          attr_d["lambda"], abs=1e-10)

    def test_bond_run_reports_lambda_attribution(self):
        """End-to-end: a bond run records the pp/fl attribution in the eigenvalue
        file provenance (S4.5/S7.10b)."""
        ci = os.path.join(self.tmp_path, "ci.dat")
        _write_w90(ci, 1, _nn_entries(0.5))
        d_on = os.path.join(self.tmp_path, "on")
        sc.calc_eliashberg(_base_input(d_on, coulomb_inter=ci, bond_channels=True))
        text = open(os.path.join(d_on, "eigenvalue.dat")).read()
        self.assertTrue("lambda_pp" in text and "lambda_fl" in text)
        vals = {}
        for line in text.splitlines():
            if line.startswith("#") and "=" in line:
                k, _, v = line.lstrip("# ").partition("=")
                vals[k.strip()] = v.strip()
        lam = float(vals["lambda_rayleigh"])
        self.assertApprox(float(vals["lambda_pp"]) + float(vals["lambda_fl"]),
                          lam, abs=1e-8)

    # --- (c) the odd f-wave seed --------------------------------------------

    def test_initialize_gap_f_wave_is_odd(self):
        """S7.7: the f-like seed sin kx (cos kx - cos ky) is ODD under k -> -k."""
        Nx = Ny = 8
        kx = np.linspace(0, 2 * np.pi, Nx, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, Ny, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, 1, endpoint=False)
        for mode in ("f_x", "f_y"):
            sigma = sc._initialize_gap(mode, 1, kx, ky, kz)
            self.assertApprox(np.linalg.norm(sigma), 1.0)
            self.assertLess(_gap_parity_residual(sigma), 1e-12)
            self.assertFalse(sc._is_gap_parity(sigma, "singlet"))
            self.assertTrue(sc._is_gap_parity(sigma, "triplet"))

        # f_x is the specified form factor, not a p-wave in disguise
        KX, KY, _ = np.meshgrid(kx, ky, kz, indexing='ij')
        ref = np.sin(KX) * (np.cos(KX) - np.cos(KY))
        ref = ref / np.linalg.norm(ref)
        got = sc._initialize_gap("f_x", 1, kx, ky, kz)[0, 0]
        np.testing.assert_allclose(got, ref, atol=1e-12)
        # and it is linearly independent of p_x
        px = sc._initialize_gap("p_x", 1, kx, ky, kz)[0, 0]
        self.assertLess(abs(np.vdot(px.ravel(), ref.ravel())), 1e-12)

    def test_initialize_gap_existing_seeds_unchanged(self):
        """No regression: the pre-existing seeds keep their parity/values."""
        Nx = Ny = 8
        kx = np.linspace(0, 2 * np.pi, Nx, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, Ny, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, 1, endpoint=False)
        KX, KY, _ = np.meshgrid(kx, ky, kz, indexing='ij')
        cases = {
            "p_x": np.sin(KX),
            "d_x2y2": np.cos(KX) - np.cos(KY),
            "cos": np.cos(KX + KY),
            "s_ext_2d": np.cos(KX) * np.cos(KY),
        }
        for mode, ref in cases.items():
            ref = ref / np.linalg.norm(ref)
            got = sc._initialize_gap(mode, 1, kx, ky, kz)[0, 0]
            np.testing.assert_allclose(got, ref, atol=1e-12)
        self.assertLess(
            _gap_parity_residual(sc._initialize_gap("p_x", 1, kx, ky, kz)), 1e-12)
        self.assertLess(
            _gap_even_residual(sc._initialize_gap("d_x2y2", 1, kx, ky, kz)), 1e-12)
        with self.assertRaisesRegex(ValueError, r"(?i)unknown init_gap"):
            sc._initialize_gap("not_a_mode", 1, kx, ky, kz)


# --- (b) runtime Hermiticity / positivity preconditions (S4.5) -------------

class TestHermitianPreconditions(_ApproxTestCase):

    def test_preconditions_pass_on_physical_green(self):
        from hwave.solver import bond_channels as bc

        f = _physical_single_band(V=1.0)
        A, _, _, n = bc.make_bond_kernel_parts(
            f["chi_s"], f["chi_c"], f["S_bond"], f["C_bond"], f["Vpp_s"],
            f["Vpp_t"], f["green"], f["bond_set"], "singlet", f["beta"])
        W = bc.pair_weight(f["green"], f["beta"])
        diag = bc.check_hermitian_preconditions(A, W)
        self.assertGreater(diag["weight_min_eigenvalue"], 0.0)
        self.assertLess(diag["weight_hermiticity_residual"], 1e-10)
        self.assertLess(diag["kernel_hermiticity_residual"], 1e-8)

    def test_preconditions_error_on_unphysical_green(self):
        """A green that violates the inversion/time-reversal symmetry makes the
        static pair weight w = GG complex/indefinite -- v1 must ERROR (no
        non-Hermitian biorthogonal fallback)."""
        from hwave.solver import bond_channels as bc

        rng = np.random.default_rng(7)
        norb, Nx, Ny, Nz, nmat = 1, 4, 4, 1, 4
        beta = 2.0
        green = (rng.normal(size=(norb, norb, Nx, Ny, Nz, nmat))
                 + 1j * rng.normal(size=(norb, norb, Nx, Ny, Nz, nmat)))
        W = bc.pair_weight(green, beta)
        ident = LinearOperator((norb * norb * Nx * Ny * Nz,) * 2,
                               matvec=lambda v: v, dtype=complex)
        with self.assertRaisesRegex(ValueError, r"(?i)weight|real|positive|Hermitian"):
            bc.check_hermitian_preconditions(ident, W)

    def test_preconditions_error_on_non_hermitian_kernel(self):
        """Even with a perfectly physical weight, a non-Hermitian symmetrized
        kernel must ERROR and name the residual."""
        from hwave.solver import bond_channels as bc

        n = 6
        rng = np.random.default_rng(11)
        M = rng.normal(size=(n, n))
        M_ns = M.copy()
        M_ns[0, 1] += 3.0                      # deliberately non-symmetric
        W = np.ones((n, 1, 1, 1, 1))           # identity metric, w = 1 >= 0
        op_ns = LinearOperator((n, n), matvec=lambda v: M_ns @ v, dtype=complex)
        with self.assertRaisesRegex(ValueError, r"(?i)hermit"):
            bc.check_hermitian_preconditions(op_ns, W)

        M_h = 0.5 * (M + M.T)
        op_h = LinearOperator((n, n), matvec=lambda v: M_h @ v, dtype=complex)
        diag = bc.check_hermitian_preconditions(op_h, W)
        self.assertLess(diag["kernel_hermiticity_residual"], 1e-10)

    def test_precondition_probe_path_matches_dense_verdict(self):
        """Above ``dense_limit`` the Hermiticity residual is estimated
        stochastically; it must still PASS a Hermitian kernel and ERROR on a
        non-Hermitian one (deterministic given the seed)."""
        from hwave.solver import bond_channels as bc

        f = _physical_single_band(V=1.0)
        A, _, _, n = bc.make_bond_kernel_parts(
            f["chi_s"], f["chi_c"], f["S_bond"], f["C_bond"], f["Vpp_s"],
            f["Vpp_t"], f["green"], f["bond_set"], "singlet", f["beta"])
        W = bc.pair_weight(f["green"], f["beta"])
        diag = bc.check_hermitian_preconditions(A, W, dense_limit=0)
        self.assertEqual(diag["method"], "probe")
        self.assertLess(diag["kernel_hermiticity_relative"], 1e-8)

        broken = LinearOperator(
            (n, n), matvec=lambda v: A.matvec(v) + 0.5 * np.roll(v, 1),
            dtype=complex)
        with self.assertRaisesRegex(ValueError, r"(?i)hermit"):
            bc.check_hermitian_preconditions(broken, W, dense_limit=0)

    def test_precondition_catches_a_violation_on_a_small_weight_index(self):
        """I2 (kept; re-referred by the I1 fix): a Hermiticity violation sitting
        on a badly-conditioned (small-``w``) index must still be caught.

        The residual is measured on ``M = W K`` but the criterion is about
        ``K~ = W^-1/2 M W^-1/2``, and the two weight such a violation differently:
        ``M`` suppresses it by ``w_i w_j`` (which HIDES it -- the original I2
        finding), ``K~`` by ``sqrt(w_i w_j)``, which is exactly its effect on
        ``Im lambda``.  The check therefore compares the EXACT ``K~`` residual
        against an equally ``K~``-referred scale (the I1 fix), and the violation
        below is caught while the same kernel without it passes.
        """
        from hwave.solver import bond_channels as bc

        n = 8
        w = np.ones(n)
        w[n - 1] = 1.0e-6                     # the badly-conditioned index
        gamma = np.eye(n, dtype=complex)
        gamma[n - 1, 0] += 1.0j               # a NON-Hermitian vertex entry
        # the kernel always factorizes as K = -Gamma W (make_bond_kernel_parts)
        K = -gamma * w[np.newaxis, :]
        op = LinearOperator((n, n), matvec=lambda v: K @ v, dtype=complex)

        M = w[:, np.newaxis] * K
        resid_M = np.linalg.norm(M - np.conj(M.T))
        scale_M = np.linalg.norm(M)
        atol = rtol = 1.0e-6
        # on M the violation is doubly suppressed (w_{n-1} w_0) and hides ...
        self.assertLess(resid_M, max(atol, rtol * scale_M))

        # ... on K~ it is only suppressed by sqrt(w_{n-1} w_0) and is caught
        s = np.sqrt(w)
        Kt = -(s[:, np.newaxis] * gamma * s[np.newaxis, :])
        self.assertGreater(
            np.linalg.norm(Kt - np.conj(Kt.T)), max(atol, rtol * np.linalg.norm(Kt)))
        with self.assertRaisesRegex(ValueError, r"(?i)hermit") as cm:
            bc.check_hermitian_preconditions(op, w, atol=atol, rtol=rtol)
        # the message quotes the EXACT K~ residual, not a /w_min bound
        self.assertIn("{:.3e}".format(
            np.linalg.norm(Kt - np.conj(Kt.T))), str(cm.exception))

        # the same kernel with a Hermitian vertex passes on the same metric
        op_ok = LinearOperator(
            (n, n), matvec=lambda v: (-np.eye(n) * w[np.newaxis, :]) @ v,
            dtype=complex)
        diag = bc.check_hermitian_preconditions(op_ok, w, atol=atol, rtol=rtol)
        self.assertLess(diag["kernel_hermiticity_residual_ktilde"], 1e-12)

    def test_precondition_criterion_is_ktilde_referred_on_both_sides(self):
        """I1 (the real bug): the residual was referred to ``K~`` (divided by
        ``w_min``) but the tolerance scale stayed ``||M||_F``, so the relative
        criterion was ``w_min`` times stricter than requested -- i.e. it tightened
        with beta.  A kernel that is Hermitian to 1e-12 RELATIVE in the ``K~``
        metric must be accepted no matter how small ``min w`` is."""
        from hwave.solver import bond_channels as bc

        n = 12
        rng = np.random.default_rng(20260726)
        # a pair weight with the wide dynamic range w = GG acquires at large beta
        w = np.geomspace(1.0, 1.0e-6, n)
        G = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        gamma = 0.5 * (G + np.conj(G.T))                    # Hermitian vertex
        A = 0.5 * (G - np.conj(G.T))                        # anti-Hermitian
        eps_rel = 1.0e-12                                   # negligible asymmetry
        gamma_noisy = gamma + eps_rel * A
        # K = -Gamma W (the factorization the bond kernel has by construction)
        K = -gamma_noisy * w[np.newaxis, :]
        op = LinearOperator((n, n), matvec=lambda v: K @ v, dtype=complex)

        # the K~ asymmetry is ~1e-12 relative -- physically fine
        s = np.sqrt(w)
        Kt = -(s[:, np.newaxis] * gamma_noisy * s[np.newaxis, :])
        self.assertLess(
            np.linalg.norm(Kt - np.conj(Kt.T)), 1.0e-10 * np.linalg.norm(Kt))

        # the OLD criterion (K~-referred residual vs an M-referred scale) refuses
        M = w[:, np.newaxis] * K
        resid_M = np.linalg.norm(M - np.conj(M.T))
        self.assertGreater(resid_M / w.min(), max(1e-8, 1e-8 * np.linalg.norm(M)))

        diag = bc.check_hermitian_preconditions(op, w)      # must NOT raise
        self.assertLess(diag["kernel_hermiticity_relative_ktilde"], 1.0e-10)

        # ... while a genuine violation on the SAME badly-conditioned metric is
        # still caught (the fix must not just disable the check)
        bad = gamma + 1.0e-4 * A
        K_bad = -bad * w[np.newaxis, :]
        op_bad = LinearOperator((n, n), matvec=lambda v: K_bad @ v, dtype=complex)
        with self.assertRaisesRegex(ValueError, r"(?i)hermit"):
            bc.check_hermitian_preconditions(op_bad, w)

    def test_precondition_scaled_tolerance_still_passes_the_physical_kernel(self):
        """I2 regression: tightening the criterion must not turn the physical
        bond kernel into a false alarm."""
        from hwave.solver import bond_channels as bc

        f = _physical_single_band(V=1.0)
        A, _, _, n = bc.make_bond_kernel_parts(
            f["chi_s"], f["chi_c"], f["S_bond"], f["C_bond"], f["Vpp_s"],
            f["Vpp_t"], f["green"], f["bond_set"], "triplet", f["beta"])
        W = bc.pair_weight(f["green"], f["beta"])
        diag = bc.check_hermitian_preconditions(A, W)
        self.assertGreater(diag["weight_min_eigenvalue"], 0.0)
        self.assertLess(diag["kernel_hermiticity_residual_ktilde"], 1e-9)

    def test_precondition_dense_branch_checks_the_kernel_parts(self):
        """M8: the S4.5 attribution needs EACH part Hermitian, not only the sum;
        the dense branch reports the part-wise residual and warns when a
        violation cancels between the parts."""
        from hwave.solver import bond_channels as bc

        f = _physical_single_band(V=1.0)
        A, A_fl, A_pp, n = bc.make_bond_kernel_parts(
            f["chi_s"], f["chi_c"], f["S_bond"], f["C_bond"], f["Vpp_s"],
            f["Vpp_t"], f["green"], f["bond_set"], "triplet", f["beta"])
        W = bc.pair_weight(f["green"], f["beta"])

        diag = bc.check_hermitian_preconditions(A, W,
                                                parts={"pp": A_pp, "fl": A_fl})
        self.assertLess(diag["part_hermiticity_residual_pp"], 1e-9)
        self.assertLess(diag["part_hermiticity_residual_fl"], 1e-9)

        # a violation that CANCELS between the parts leaves the full kernel
        # Hermitian -- only the part-wise check can see it
        def _bad(op, sgn):
            return LinearOperator(
                (n, n), matvec=lambda v, op=op, sgn=sgn: (
                    op.matvec(v) + sgn * 0.5 * np.roll(v, 1)),
                dtype=complex)

        with _caplog(logging.WARNING) as records:
            diag2 = bc.check_hermitian_preconditions(
                A, W, parts={"pp": _bad(A_pp, +1.0), "fl": _bad(A_fl, -1.0)})
        self.assertGreater(diag2["part_hermiticity_residual_pp"], 1e-3)
        msgs = " ".join(r.getMessage() for r in records)
        self.assertTrue("pp" in msgs and "hermit" in msgs.lower())

        # the probe branch does not build the parts densely
        diag3 = bc.check_hermitian_preconditions(A, W, dense_limit=0,
                                                 parts={"pp": A_pp, "fl": A_fl})
        self.assertNotIn("part_hermiticity_residual_pp", diag3)


# --- (d) invariant-subspace tracking (S7.7) -------------------------------

def _basis_vec(n, idx):
    v = np.zeros(n, dtype=complex)
    v[idx] = 1.0
    return v


class TestTrackSubspace(_ApproxTestCase):

    def test_track_subspace_follows_degenerate_branch_not_single_vector(self):
        """S7.7 (iii)+(iv): near-degenerate eigenpairs are clustered into an
        invariant subspace and tracked by PRINCIPAL-ANGLE overlap of the whole
        subspace.  This test is discriminating: a naive "max overlap of one
        eigenvector" rule follows the WRONG branch here, because the degenerate
        pair comes back in a rotated basis (each rotated vector only has 1/sqrt(2)
        overlap with the previous single vector) while a nearby non-degenerate
        competitor has 0.8."""
        from hwave.solver import bond_channels as bc

        n = 5
        e1, e2, e3 = (_basis_vec(n, 0), _basis_vec(n, 1), _basis_vec(n, 2))
        e4 = _basis_vec(n, 3)
        p0 = (np.array([1.0, 1.0, 0.99, 0.5]),
              np.stack([e1, e2, e3, e4]))
        u1 = (e1 + e2) / np.sqrt(2.0)
        u2 = (e1 - e2) / np.sqrt(2.0)
        competitor = 0.8 * e1 + 0.6 * e3
        p1 = (np.array([1.02, 1.02, 1.01, 0.4]),
              np.stack([u1, u2, competitor, e4]))

        out = bc.track_subspace([p0, p1], seed=e1, deg_tol=1e-3)
        self.assertEqual(len(out), 2)
        self.assertApprox(out[0]["lambda"], 1.0)
        self.assertEqual(out[0]["dim"], 2)
        # the tracked branch is the 2D degenerate one (1.02), NOT the 1.01
        # competitor a single-vector rule would pick
        self.assertApprox(out[1]["lambda"], 1.02)
        self.assertEqual(out[1]["dim"], 2)
        self.assertApprox(out[1]["overlap"], 1.0, abs=1e-10)
        self.assertTrue(out[1]["captured"])

        # the naive rule really does fail here (this is what makes the test
        # discriminating, not just passing)
        naive = max(range(4), key=lambda i: abs(np.vdot(e1, p1[1][i])))
        self.assertApprox(p1[0][naive], 1.01)

    def test_track_subspace_anchors_on_f_seed_not_leading_lambda(self):
        """S7.7 (v): at V = 0 the tracked subspace is the one with maximal overlap
        with the odd f-seed -- not simply the largest eigenvalue."""
        from hwave.solver import bond_channels as bc

        n = 4
        e1, e2 = _basis_vec(n, 0), _basis_vec(n, 1)
        pts = [(np.array([2.0, 0.5]), np.stack([e1, e2]))]
        out = bc.track_subspace(pts, seed=e2, deg_tol=1e-3)
        self.assertApprox(out[0]["lambda"], 0.5)
        out2 = bc.track_subspace(pts, seed=e1, deg_tol=1e-3)
        self.assertApprox(out2[0]["lambda"], 2.0)

    def test_track_subspace_adaptive_n_ev(self):
        """S7.7: when the tracked subspace is not captured by the current n_ev the
        tracker re-solves with a larger n_ev (points given as callables)."""
        from hwave.solver import bond_channels as bc

        n = 6
        e = [_basis_vec(n, i) for i in range(n)]
        calls = []

        def point0(n_ev):
            calls.append(("p0", n_ev))
            return np.array([1.0]), np.stack([e[0]])

        def point1(n_ev):
            calls.append(("p1", n_ev))
            # the tracked state is only the 4th eigenpair at this V
            vals = np.array([3.0, 2.0, 1.5, 1.1, 0.2, 0.1])[:n_ev]
            vecs = np.stack([e[1], e[2], e[3], e[0], e[4], e[5]])[:n_ev]
            return vals, vecs

        out = bc.track_subspace([point0, point1], seed=e[0], deg_tol=1e-3,
                                n_ev=2, max_n_ev=8, capture_tol=0.5)
        self.assertApprox(out[1]["lambda"], 1.1)
        self.assertTrue(out[1]["captured"])
        self.assertGreater(out[1]["n_ev"], 2)               # it escalated
        self.assertTrue(any(name == "p1" and k > 2 for name, k in calls))

    def test_track_subspace_harmonic_decomposition(self):
        """S7.7: the tracked subspace is reported as a decomposition onto the
        normalized odd basis {sin kx, sin ky, sin kx(cos kx - cos ky),
        sin ky(cos ky - cos kx)}."""
        from hwave.solver import bond_channels as bc

        Nx = Ny = 8
        kx = np.linspace(0, 2 * np.pi, Nx, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, Ny, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, 1, endpoint=False)
        basis = sc._odd_harmonic_basis(1, kx, ky, kz)
        self.assertEqual(set(basis), {"p_x", "p_y", "f_x", "f_y"})

        f_x = basis["f_x"]
        p_x = basis["p_x"]
        mix = (f_x + p_x) / np.linalg.norm(f_x + p_x)

        pts = [(np.array([1.0]), np.stack([f_x])),
               (np.array([1.1]), np.stack([mix]))]
        out = bc.track_subspace(pts, seed=f_x, deg_tol=1e-3, basis=basis)
        self.assertApprox(out[0]["harmonics"]["f_x"], 1.0, abs=1e-8)
        self.assertApprox(out[0]["harmonics"]["p_x"], 0.0, abs=1e-8)
        self.assertApprox(out[1]["harmonics"]["f_x"], 0.5, abs=1e-8)
        self.assertApprox(out[1]["harmonics"]["p_x"], 0.5, abs=1e-8)
        self.assertApprox(out[1]["harmonics"]["p_y"], 0.0, abs=1e-8)

    def test_track_subspace_uses_the_gg_metric(self):
        """S7.7 (i): the inner product is the sqrt(GG)-weighted Hermitian form, so
        a non-uniform weight changes which branch is tracked."""
        from hwave.solver import bond_channels as bc

        n = 3
        e0, e1, _ = (_basis_vec(n, i) for i in range(3))
        seed = 0.8 * e0 + 0.6 * e1
        pts = [(np.array([1.0, 0.9]), np.stack([e0, e1]))]
        # identity metric: the seed leans on e0, so the lambda=1.0 branch wins.
        out = bc.track_subspace(pts, seed=seed, weight=np.ones(n), deg_tol=1e-3)
        self.assertApprox(out[0]["lambda"], 1.0)
        # a metric that suppresses component 0 (as a vanishing GG weight would)
        # flips the anchor onto the lambda=0.9 branch.
        out2 = bc.track_subspace(pts, seed=seed, deg_tol=1e-3,
                                 weight=np.array([1e-6, 1.0, 1.0]))
        self.assertApprox(out2[0]["lambda"], 0.9)

    def test_track_subspace_flags_a_multiplet_cut_by_the_eigenvalue_window(self):
        """REGRESSION: a tracked cluster that reaches the low end of the supplied
        eigenvalue window is only PART of the multiplet, so its basis (and hence
        overlap / dim / harmonics) is an arbitrary slice.  ``track_subspace`` must
        say so instead of reporting it as a clean capture."""
        from hwave.solver import bond_channels as bc

        n = 4
        e = [_basis_vec(n, i) for i in range(n)]
        # window cut inside a degenerate pair: only e[1] of the {e[1], e[2]}
        # multiplet at lambda = 1.0 is inside
        pts = [(np.array([2.0, 1.0]), np.stack([e[0], e[1]]))]
        out = bc.track_subspace(pts, seed=e[1], deg_tol=1e-3)
        self.assertIs(out[0]["truncated"], True)
        # widening the window to hold the whole multiplet clears the flag
        pts2 = [(np.array([2.0, 1.0, 1.0, 0.1]),
                 np.stack([e[0], e[1], e[2], e[3]]))]
        out2 = bc.track_subspace(pts2, seed=e[1], deg_tol=1e-3)
        self.assertIs(out2[0]["truncated"], False)
        self.assertEqual(out2[0]["dim"], 2)

    def test_track_subspace_escalates_n_ev_on_a_window_truncated_cluster(self):
        """A truncated cluster escalates ``n_ev`` for callable points even when its
        overlap already clears ``capture_tol`` -- otherwise a multiplet cut in half
        is reported as captured while its basis is solver-state dependent."""
        from hwave.solver import bond_channels as bc

        n = 5
        e = [_basis_vec(n, i) for i in range(n)]
        calls = []

        def point(n_ev):
            calls.append(n_ev)
            vals = np.array([2.0, 1.0, 1.0, 0.1])[:n_ev]
            vecs = np.stack([e[0], e[1], e[2], e[3]])[:n_ev]
            return vals, vecs

        out = bc.track_subspace([point], seed=e[1], deg_tol=1e-3, n_ev=2,
                                max_n_ev=8, capture_tol=0.5)
        # at n_ev = 2 the overlap with e[1] is a perfect 1.0 -- capture_tol alone
        # would have stopped there and returned the truncated half-multiplet
        self.assertTrue(calls[0] == 2 and max(calls) > 2)
        self.assertEqual(out[0]["dim"], 2)                  # the FULL multiplet
        self.assertIs(out[0]["truncated"], False)
        self.assertTrue(out[0]["captured"])

    def test_track_subspace_is_documented_as_the_sweep_driver_api(self):
        """``track_subspace`` is inherently MULTI-POINT: a single
        ``calc_eliashberg`` run has one point, so it stays a documented public
        helper for sweep drivers rather than being called from the solver."""
        from hwave.solver import bond_channels as bc

        doc = bc.track_subspace.__doc__
        self.assertIsNotNone(doc)
        lowered = doc.lower()
        self.assertIn("sweep", lowered)
        self.assertIn("bond_diagnostics", doc)
        # the ROW convention (scipy eigs returns COLUMNS) must stay documented
        self.assertTrue("column" in lowered and "row" in lowered)


# --- covering tests: kernel split, probe path, metric guards --------------

class TestKernelPartsAndMetricGuards(_ApproxTestCase):

    def test_kernel_parts_sum_to_the_full_operator(self):
        """K = K^fl + K^pp as an OPERATOR identity (not just on the eigenvector),
        and ``make_bond_kernel(part=...)`` selects the same pieces."""
        from hwave.solver import bond_channels as bc

        f = _physical_single_band(V=0.7)
        args = (f["chi_s"], f["chi_c"], f["S_bond"], f["C_bond"], f["Vpp_s"],
                f["Vpp_t"], f["green"], f["bond_set"], "triplet", f["beta"])
        A, A_fl, A_pp, n = bc.make_bond_kernel_parts(*args)
        A_only, n2 = bc.make_bond_kernel(*args)
        A_fl2, _ = bc.make_bond_kernel(*args, part="fluctuation")
        A_pp2, _ = bc.make_bond_kernel(*args, part="instantaneous")
        self.assertEqual(n, n2)

        rng = np.random.default_rng(19)
        for _ in range(3):
            v = rng.normal(size=n) + 1j * rng.normal(size=n)
            full = A.matvec(v)
            np.testing.assert_allclose(full, A_only.matvec(v), atol=1e-12)
            np.testing.assert_allclose(A_fl.matvec(v) + A_pp.matvec(v), full,
                                       atol=1e-12)
            np.testing.assert_allclose(A_fl2.matvec(v), A_fl.matvec(v), atol=1e-12)
            np.testing.assert_allclose(A_pp2.matvec(v), A_pp.matvec(v), atol=1e-12)

        with self.assertRaisesRegex(ValueError, r"(?i)unknown part"):
            bc.make_bond_kernel(*args, part="nonsense")

    def test_attribute_lambda_rejects_a_null_metric_vector(self):
        from hwave.solver import bond_channels as bc

        n = 4
        zero = LinearOperator((n, n), matvec=lambda v: 0.0 * v, dtype=complex)
        v = np.zeros(n, dtype=complex)
        v[0] = 1.0
        with self.assertRaisesRegex(ValueError, r"(?i)rayleigh|positive|W-norm"):
            bc.attribute_lambda(v, np.array([0.0, 1.0, 1.0, 1.0]), zero, zero)

    def test_cluster_eigenvalues_and_subspace_similarity_basics(self):
        from hwave.solver import bond_channels as bc

        cl = bc.cluster_eigenvalues(np.array([1.0, 0.5, 1.0005, 0.4995]),
                                    deg_tol=1e-3)
        self.assertEqual([sorted(c) for c in cl], [[0, 2], [1, 3]])
        # a tighter tolerance splits them
        cl2 = bc.cluster_eigenvalues(np.array([1.0, 0.5, 1.0005, 0.4995]),
                                     deg_tol=1e-6)
        self.assertEqual([len(c) for c in cl2], [1, 1, 1, 1])

        # subspace similarity is invariant under a rotation inside the subspace
        rng = np.random.default_rng(4)
        Q = np.linalg.qr(rng.normal(size=(5, 5)))[0]
        S = Q[:2]
        theta = 0.7
        R = np.array([[np.cos(theta), -np.sin(theta)],
                      [np.sin(theta), np.cos(theta)]])
        s1, _ = bc.subspace_similarity(S, R @ S)
        self.assertApprox(s1, 1.0, abs=1e-12)
        s2, cos = bc.subspace_similarity(S, Q[2:4])
        self.assertApprox(s2, 0.0, abs=1e-12)
        self.assertEqual(len(cos), 2)

    def test_cluster_eigenvalues_uses_the_complex_distance(self):
        """M6: two eigenvalues with the same real part but different imaginary
        parts are NOT degenerate; clustering must use |lambda_i - lambda_j|."""
        from hwave.solver import bond_channels as bc

        ev = np.array([1.0 + 0.0j, 1.0 + 0.5j])
        self.assertEqual(
            [len(c) for c in bc.cluster_eigenvalues(ev, deg_tol=1e-3)], [1, 1])
        # genuinely degenerate complex pairs still group
        ev2 = np.array([1.0 + 0.5j, 1.0 + 0.5j + 1e-6])
        self.assertEqual(
            [len(c) for c in bc.cluster_eigenvalues(ev2, deg_tol=1e-3)], [2])

    def test_subspace_similarity_tie_resolves_to_the_leading_cluster(self):
        """M5: the documented tiebreak -- equal mean-cos^2 scores go to the FIRST
        cluster in lambda-descending order."""
        from hwave.solver import bond_channels as bc

        n = 4
        e0, e1 = _basis_vec(n, 0), _basis_vec(n, 1)
        seed = (e0 + e1) / np.sqrt(2.0)          # exactly 0.5 with both clusters
        pts = [(np.array([2.0, 1.0]), np.stack([e0, e1]))]
        out = bc.track_subspace(pts, seed=seed, deg_tol=1e-3)
        self.assertApprox(out[0]["overlap"], 0.5, abs=1e-12)
        self.assertApprox(out[0]["lambda"], 2.0)


# ---------------------------------------------------------------------------
# Task 7 review fixes: C1 (iteration-mode lambda vs lambda_rayleigh),
# I2 (precondition tolerance scaled by the weight conditioning), I3 (subspace
# tracking on a REAL kernel with the 5-D pair-weight metric), M5-M8.
# ---------------------------------------------------------------------------

def _provenance_keys(path):
    """``# key = value`` comment lines of an eigenvalue file, as a dict."""
    vals = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") and "=" in line:
                k, _, v = line.lstrip("# ").partition("=")
                vals[k.strip()] = v.strip()
    return vals


class TestIterationModeLambdaProvenance(_ApproxTestCase):

    def setUp(self):
        self.tmp_path = tempfile.mkdtemp(prefix="hwave_sc_bond_")
        self.addCleanup(shutil.rmtree, self.tmp_path, ignore_errors=True)

    def test_iteration_mode_lambda_is_flagged_against_lambda_rayleigh(self):
        """C1: on a repulsive-dominant bond kernel the power iterate's norm
        (positive, unconverged) and the signed Rayleigh quotient
        ``lambda_rayleigh`` (negative) land in the SAME file. The artifact must
        record which solver mode produced the number and whether it converged, the
        run must warn, and the file must say in words that the iteration value is
        an unsigned iterate norm."""
        ci = os.path.join(self.tmp_path, "ci.dat")
        _write_w90(ci, 1, _nn_entries(0.5))
        d_on = os.path.join(self.tmp_path, "on")
        with _caplog(logging.WARNING) as records:
            sc.calc_eliashberg(_base_input(d_on, coulomb_inter=ci,
                                           bond_channels=True,
                                           solver_mode="iteration"))
        path = os.path.join(d_on, "eigenvalue.dat")
        vals = _provenance_keys(path)

        # (a) provenance records the mode and the convergence flag
        self.assertEqual(vals["lambda_rayleigh_solver_mode"], "iteration")
        self.assertEqual(vals["lambda_rayleigh_converged"], "False")

        # the situation this guards against really occurs on this kernel: the two
        # numbers in the file have OPPOSITE signs
        lam_rayleigh = float(vals["lambda_rayleigh"])
        lam_iter = _leading_eigenvalue(d_on)
        self.assertTrue(lam_rayleigh < 0.0 < lam_iter)

        # ... and an INDEPENDENT computation of the same kernel says which of the
        # two is the physical eigenvalue: run the identical input through
        # solver_mode="eigenvalue" and take its dominant (largest-|.|) eigenvalue.
        d_ev = os.path.join(self.tmp_path, "ev")
        sc.calc_eliashberg(_base_input(d_ev, coulomb_inter=ci,
                                       bond_channels=True,
                                       solver_mode="eigenvalue",
                                       num_eigenvalues=10))
        # (the eigenvalue-mode file carries BOTH a single-column scalar row and the
        # multi-column spectrum, so it is parsed row by row rather than by loadtxt)
        spectrum = [[float(x) for x in l.split()]
                    for l in open(os.path.join(d_ev, "eigenvalue.dat"))
                    if l.strip() and not l.startswith("#")]
        spectrum = np.array([r for r in spectrum if len(r) >= 4])
        dominant = spectrum[np.argmax(spectrum[:, 3]), 1]      # Re at max |ev|
        self.assertLess(dominant, 0.0)
        # the SIGNED Rayleigh quotient reproduces it (the power iterate's DIRECTION
        # does converge; only its reported norm is unsigned) ...
        self.assertApprox(lam_rayleigh, dominant, rel=1e-6)
        # ... while the iteration-mode number in the file does not, in magnitude or
        # in sign -- exactly the confusion the warning and the note exist for.
        self.assertGreater(abs(lam_iter - abs(dominant)), 1.0)

        # (b) a warning names both quantities and recommends the eigenvalue mode
        msgs = " ".join(r.getMessage() for r in records)
        self.assertIn("lambda_rayleigh", msgs)
        self.assertIn("unsigned", msgs.lower())
        self.assertTrue("solver_mode" in msgs and "eigenvalue" in msgs)

        # (c) the artifact itself distinguishes the two numbers
        lines = open(path).read().splitlines()
        idx = lines.index("# Iteration eigenvalue")
        note = " ".join(l for l in lines[idx + 1:] if l.startswith("#")).lower()
        self.assertIn("unsigned", note)
        self.assertIn("lambda_rayleigh", note)
        self.assertIn("did not converge", note)
        # ... and the numeric row keeps its documented place and shape: the file
        #     holds EXACTLY ONE non-comment line, it sits inside the
        #     "# Iteration eigenvalue" block (after the note, with no further
        #     numeric row), it has a single column, and that column is the value
        #     `_leading_eigenvalue` (and every downstream reader) picks up.
        numeric = [(i, l) for i, l in enumerate(lines)
                   if l.strip() and not l.startswith("#")]
        self.assertEqual(len(numeric), 1)
        row_idx, row = numeric[0]
        self.assertGreater(row_idx, idx)           # after "# Iteration eigenvalue"
        self.assertTrue(all(lines[j].startswith("#") for j in range(idx + 1, row_idx)))
        self.assertEqual(len(row.split()), 1)
        self.assertEqual(float(row), lam_iter)

    def test_iteration_with_spectral_shift_converges_and_matches_eigenvalue_mode(self):
        """The payoff of the power-iteration spectral shift: the repulsive-dominant
        bond kernel (dominant eigenvalue ~ -2.75) makes the plain power iteration
        flip sign every step and never converge. Iterating on K + sigma*I and
        subtracting sigma back converges to the SIGNED physical eigenvalue, which
        must agree with ``solver_mode="eigenvalue"`` and with the Task-7
        ``lambda_rayleigh`` attribution."""
        ci = os.path.join(self.tmp_path, "ci.dat")
        _write_w90(ci, 1, _nn_entries(0.5))
        # d_x2y2 seed: the leading in-sector eigenvector of this kernel. A seed
        # orthogonal to it (the default "cos") converges to a sub-leading mode --
        # a property of power iteration itself, unrelated to the shift.
        common = dict(coulomb_inter=ci, bond_channels=True, init_gap="d_x2y2")

        d_it = os.path.join(self.tmp_path, "iter")
        with _caplog(logging.WARNING) as records:
            sc.calc_eliashberg(_base_input(d_it, solver_mode="iteration",
                                           spectral_shift="auto", alpha=0.0,
                                           max_iter=3000, convergence_tol=1e-9,
                                           **common))
        d_ev = os.path.join(self.tmp_path, "eig")
        sc.calc_eliashberg(_base_input(d_ev, solver_mode="eigenvalue",
                                       num_eigenvalues=10, **common))

        vals_it = _provenance_keys(os.path.join(d_it, "eigenvalue.dat"))
        vals_ev = _provenance_keys(os.path.join(d_ev, "eigenvalue.dat"))

        # (a) it converges now (it structurally could not before)
        self.assertEqual(vals_it["lambda_rayleigh_converged"], "True")

        lam_it = _leading_eigenvalue(d_it)
        lam_ev = _leading_eigenvalue(d_ev)
        # (b) the reported number is the un-shifted, signed eigenvalue and agrees
        #     with the eigenvalue solver ...
        self.assertApprox(lam_it, lam_ev, rel=1e-6)
        # (c) ... and with the signed Rayleigh quotient of the same run
        self.assertApprox(lam_it, float(vals_it["lambda_rayleigh"]), rel=1e-6)
        self.assertApprox(float(vals_it["lambda_rayleigh"]),
                          float(vals_ev["lambda_rayleigh"]), rel=1e-6)

        # (d) the "unsigned iterate norm" warning must NOT fire for a shifted,
        #     converged run -- the value really is a signed eigenvalue
        msgs = " ".join(r.getMessage() for r in records)
        self.assertNotIn("unsigned", msgs.lower())

    def test_eigenvalue_mode_provenance_marks_convergence_not_applicable(self):
        """C1: with ``solver_mode='eigenvalue'`` there is no power iterate, so the
        convergence flag is ``n/a`` and no unsigned-norm warning is emitted."""
        ci = os.path.join(self.tmp_path, "ci.dat")
        _write_w90(ci, 1, _nn_entries(0.5))
        d_on = os.path.join(self.tmp_path, "on")
        with _caplog(logging.WARNING) as records:
            sc.calc_eliashberg(_base_input(d_on, coulomb_inter=ci,
                                           bond_channels=True,
                                           solver_mode="eigenvalue",
                                           num_eigenvalues=3))
        vals = _provenance_keys(os.path.join(d_on, "eigenvalue.dat"))
        self.assertEqual(vals["lambda_rayleigh_solver_mode"], "eigenvalue")
        self.assertEqual(vals["lambda_rayleigh_converged"], "n/a")
        msgs = " ".join(r.getMessage() for r in records)
        self.assertNotIn("unsigned iterate norm", msgs.lower())


class TestPreconditionOptionsThreading(_ApproxTestCase):

    def setUp(self):
        self.tmp_path = tempfile.mkdtemp(prefix="hwave_sc_bond_")
        self.addCleanup(shutil.rmtree, self.tmp_path, ignore_errors=True)

    def test_precondition_options_are_threaded_from_the_eliashberg_config(self):
        """M7: ``bond_precondition_atol/rtol/dense_limit`` reach
        ``check_hermitian_preconditions`` so the acceptance run can tune them."""
        ci = os.path.join(self.tmp_path, "ci.dat")
        _write_w90(ci, 1, _nn_entries(0.5))

        # dense_limit = 0 forces the randomized probe estimator
        d_probe = os.path.join(self.tmp_path, "probe")
        with _caplog(logging.INFO, logger="hwave_sc") as records:
            sc.calc_eliashberg(_base_input(d_probe, coulomb_inter=ci,
                                           bond_channels=True,
                                           bond_precondition_dense_limit=0))
        msgs = " ".join(r.getMessage() for r in records)
        self.assertIn("probe check", msgs)

        # an impossible tolerance must actually reach the check and fail the run
        d_tight = os.path.join(self.tmp_path, "tight")
        with self.assertRaisesRegex(ValueError, r"(?i)hermit"):
            sc.calc_eliashberg(_base_input(d_tight, coulomb_inter=ci,
                                           bond_channels=True,
                                           bond_precondition_atol=0.0,
                                           bond_precondition_rtol=1.0e-30))

    def test_precondition_options_ignored_with_warning_when_flag_off(self):
        """M7: like the other bond options, the precondition knobs are ignored
        with a warning when ``bond_channels=false``."""
        d = os.path.join(self.tmp_path, "off")
        with _caplog(logging.WARNING) as records:
            sc.calc_eliashberg(_base_input(d, bond_channels=False,
                                           bond_precondition_atol=1e-3))
        msgs = " ".join(r.getMessage() for r in records)
        self.assertIn("bond_precondition_atol", msgs)


# number of ARPACK eigenpairs used by the physical-kernel tracking tests.
#
# The kernel is 16 x 16 and its odd sector at V = 0 carries a THREE-fold
# degenerate level at lambda = 0.4265276.  With k = 6 the window ends exactly
# inside that multiplet, so ARPACK returns an arbitrary TWO-dimensional slice
# of a three-dimensional eigenspace -- and which slice depends on the state of
# ARPACK's process-global saved start-vector seed, i.e. on how many eigensolves
# ran earlier in the same test process.  k = 8 puts the whole multiplet
# inside the window at both V points (the next level, 0.2374 / 0.2444, is well
# separated), so the tracked subspace is the true invariant subspace and the
# result is reproducible.
_PHYS_TRACK_NEV = 8

# ARPACK's start vector is process-global state when ``v0`` is omitted, so the
# tests pin it -- exactly as ``sc.py`` does for its eigenvector-producing
# solves (``eigs(..., v0=seed_vec)``).  A pseudo-random vector (rather than
# e.g. ones) is used so it is not accidentally orthogonal to the odd sector.
_PHYS_TRACK_V0_SEED = 20260727

# deg_tol for the physical two-point sweep.  Switching on V = 0.5 splits the
# odd multiplet by only 4.3e-3 (0.4428830 / 0.4386173) while the nearest
# unrelated gap is 1.8e-2, so the tracker must use a tolerance between the two
# to follow the multiplet as ONE invariant subspace instead of picking an
# arbitrary sub-branch of it.
_PHYS_TRACK_DEG_TOL = 1.0e-2


def _physical_tracked_records():
    """Build the two-point physical bond sweep and track it.

    Returns ``(records, weight, shape)``.  Every input is constructed from
    scratch on each call and the eigensolver start vector is pinned, so two
    calls must agree bit-for-bit regardless of what else ran in the process.

    ``scipy.sparse.linalg.eigs`` returns eigenvectors as COLUMNS while
    ``track_subspace`` wants one per ROW -- hence the transpose."""
    from hwave.solver import bond_channels as bc
    from scipy.sparse.linalg import eigs

    points = []
    weight = seed = basis = None
    for V in (0.0, 0.5):
        f = _physical_single_band(V=V)
        A, _, _, n = bc.make_bond_kernel_parts(
            f["chi_s"], f["chi_c"], f["S_bond"], f["C_bond"], f["Vpp_s"],
            f["Vpp_t"], f["green"], f["bond_set"], "triplet", f["beta"])
        v0 = np.random.default_rng(_PHYS_TRACK_V0_SEED).standard_normal(n)
        vals, vecs = eigs(A, k=_PHYS_TRACK_NEV, which="LR", maxiter=20000,
                          tol=0, v0=v0)
        points.append((vals, vecs.T))          # ROWS, not columns
        if weight is None:
            kx, ky, kz = f["kgrid"]
            weight = bc.pair_weight(f["green"], f["beta"])
            basis = sc._odd_harmonic_basis(f["norb"], kx, ky, kz)
            seed = sc._initialize_gap("f_x", f["norb"], kx, ky, kz).ravel()

    recs = bc.track_subspace(points, seed=seed, weight=weight, basis=basis,
                             deg_tol=_PHYS_TRACK_DEG_TOL)
    return recs, weight, f["shape"]


class TestTrackSubspacePhysicalKernel(_ApproxTestCase):

    def test_track_subspace_on_the_physical_kernel_with_the_pair_weight_metric(self):
        """I3: exercise the tracker end-to-end on a REAL bond kernel with the 5-D
        ``pair_weight`` metric (the Task-8 interface), at two ``V`` points."""
        from hwave.solver import bond_channels as bc

        recs, weight, shape = _physical_tracked_records()

        # the metric really is the 5-D per-k orbital-pair form
        self.assertTrue(weight.ndim == 5 and weight.shape[:3] == shape)

        self.assertEqual(len(recs), 2)
        # the eigenvalue window contains the WHOLE tracked multiplet -- nothing is
        # cut off by n_ev (which would make the tracked basis arbitrary)
        self.assertFalse(any(rec["truncated"] for rec in recs))

        for rec in recs:
            self.assertTrue(rec["captured"])
            # well-formed harmonic decomposition: each fraction in [0, 1], and the
            # total cannot exceed the dimension of the tracked subspace
            for name, frac in rec["harmonics"].items():
                self.assertTrue(-1e-10 <= frac <= 1.0 + 1e-10, name)
            self.assertLessEqual(sum(rec["harmonics"].values()), rec["dim"] + 1e-8)
            # W-orthonormal basis of the tracked subspace (5-D metric applied)
            gram = bc._metric_inner(weight, rec["vectors"], rec["vectors"])
            np.testing.assert_allclose(gram, np.eye(rec["dim"]), atol=1e-8)

        # the tracked cluster is STABLE across the two V points: same dimension,
        # dominantly f-wave at both, and a high principal-angle overlap
        self.assertEqual(recs[0]["dim"], recs[1]["dim"])
        self.assertGreater(recs[1]["overlap"], 0.5)
        for rec in recs:
            f_wave = rec["harmonics"]["f_x"] + rec["harmonics"]["f_y"]
            p_wave = rec["harmonics"]["p_x"] + rec["harmonics"]["p_y"]
            self.assertGreater(f_wave, 10.0 * max(p_wave, 1e-12))

    def test_physical_tracking_is_independent_of_prior_eigensolves(self):
        """REGRESSION: the physical tracked sweep must not depend on how many
        eigensolves ran earlier in the process.

        ``scipy.sparse.linalg.eigs`` keeps its Arnoldi start-vector seed in
        process-global (Fortran ``SAVE``) state that EVERY call advances, so a
        solve with no ``v0`` silently depends on test execution order.  Combined
        with an ``n_ev`` window that cuts a degenerate multiplet, that used to move
        the tracked subspace's seed overlap between 0.08 and 0.89 -- the same test
        passed alone and failed after other files had run.  Here the whole
        computation is rebuilt from scratch twice with unrelated ARPACK solves
        churning the global seed in between; the tracked subspace and every number
        derived from it must come back the same (ARPACK still jitters the returned
        eigenpairs at the 1e-14 level, which only permutes members INSIDE a
        degenerate cluster -- hence the sorted comparison of ``cluster``)."""
        from scipy.sparse.linalg import eigs

        first, _, _ = _physical_tracked_records()

        # churn ARPACK's process-global start-vector state
        for _ in range(3):
            eigs(np.diag(np.arange(1.0, 12.0)), k=3, which="LR")

        second, _, _ = _physical_tracked_records()

        self.assertEqual(len(first), len(second))
        for a, b in zip(first, second):
            self.assertEqual(a["dim"], b["dim"])
            self.assertEqual(sorted(a["cluster"]), sorted(b["cluster"]))
            self.assertEqual(a["captured"], b["captured"])
            self.assertEqual(a["truncated"], b["truncated"])
            self.assertApprox(a["lambda"], b["lambda"], rel=1e-10)
            # THE regression: this used to swing over 0.08 .. 0.89
            self.assertApprox(a["overlap"], b["overlap"], rel=1e-8)
            np.testing.assert_allclose(np.sort_complex(a["eigenvalues"]),
                                       np.sort_complex(b["eigenvalues"]),
                                       rtol=1e-10, atol=1e-12)
            self.assertEqual(a["harmonics"].keys(), b["harmonics"].keys())
            for name in a["harmonics"]:
                self.assertApprox(a["harmonics"][name], b["harmonics"][name],
                                  rel=1e-8, abs=1e-12, msg=name)


# ---------------------------------------------------------------------------
# [eliashberg] bond_green -- externally supplied Green function (spec Goal)
# ---------------------------------------------------------------------------

def _bare_green_npz(path, inp):
    """Write the SAME bare Green function ``calc_eliashberg`` would build for
    ``inp``, in the H-wave ``(nblock, nfreq, nvol, norb, norb)`` npz layout."""
    mode_param = inp["mode"]["param"]
    beta = 1.0 / mode_param["T"]
    nmat = mode_param["Nmat"]
    Lx, Ly, Lz = mode_param["CellShape"]
    geom_info, hr, _ = sc._read_interaction_files(inp)
    norb = geom_info["norb"]
    kx = np.linspace(0, 2.0 * np.pi, Lx, endpoint=False)
    ky = np.linspace(0, 2.0 * np.pi, Ly, endpoint=False)
    kz = np.linspace(0, 2.0 * np.pi, Lz, endpoint=False)
    eps = sc._build_hamiltonian_k(kx, ky, kz, hr, norb)
    evals, evecs = sc._calc_eigenvalues(eps)
    mu = sc._determine_mu(evals, beta, mode_param["filling"], norb)
    green = sc._calc_green(evals, evecs, mu, beta, nmat)
    raw = green.transpose(5, 2, 3, 4, 0, 1).reshape(
        1, nmat, Lx * Ly * Lz, norb, norb)
    np.savez(path, green=raw)
    return green


def _numeric_rows(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                rows.append([float(x) for x in line.split()])
    return rows


class TestBondGreenExternalGreenFunction(_ApproxTestCase):

    def setUp(self):
        self.tmp_path = tempfile.mkdtemp(prefix="hwave_sc_bond_")
        self.addCleanup(shutil.rmtree, self.tmp_path, ignore_errors=True)

    def test_bond_green_is_actually_consumed(self):
        """``bond_green`` replaces the bare ``_calc_green``: feeding back exactly
        the bare green reproduces the no-``bond_green`` run number for number
        (proving the layout round-trip is right), while a different green gives a
        different lambda (proving the file is really used)."""
        ci = os.path.join(self.tmp_path, "ci.dat")
        _write_w90(ci, 1, _nn_entries(0.5))

        d_bare = os.path.join(self.tmp_path, "bare")
        inp_bare = _base_input(d_bare, coulomb_inter=ci, bond_channels=True,
                               solver_mode="eigenvalue", num_eigenvalues=3)
        sc.calc_eliashberg(inp_bare)

        gpath = os.path.join(self.tmp_path, "green.npz")
        green = _bare_green_npz(gpath, inp_bare)

        d_ext = os.path.join(self.tmp_path, "ext")
        sc.calc_eliashberg(_base_input(d_ext, coulomb_inter=ci,
                                       bond_channels=True, bond_green=gpath,
                                       solver_mode="eigenvalue",
                                       num_eigenvalues=3))
        got = _numeric_rows(os.path.join(d_ext, "eigenvalue.dat"))
        ref = _numeric_rows(os.path.join(d_bare, "eigenvalue.dat"))
        self.assertEqual(len(got), len(ref))
        for row_got, row_ref in zip(got, ref):
            # ARPACK is not bit-reproducible; everything upstream of it is, so an
            # agreement at 1e-10 pins the green (a different green moves lambda by
            # O(1), see below)
            np.testing.assert_allclose(row_got, row_ref, rtol=0, atol=1e-10)

        # a DIFFERENT green must move lambda (and must not be silently ignored)
        gpath2 = os.path.join(self.tmp_path, "green_scaled.npz")
        raw = np.load(gpath)["green"]
        np.savez(gpath2, green=0.8 * raw)
        d_ext2 = os.path.join(self.tmp_path, "ext2")
        sc.calc_eliashberg(_base_input(d_ext2, coulomb_inter=ci,
                                       bond_channels=True, bond_green=gpath2,
                                       solver_mode="eigenvalue",
                                       num_eigenvalues=3))
        lam_ref = _numeric_rows(os.path.join(d_bare, "eigenvalue.dat"))[0][0]
        lam_alt = _numeric_rows(os.path.join(d_ext2, "eigenvalue.dat"))[0][0]
        self.assertGreater(abs(lam_alt - lam_ref), 1e-6 * max(1.0, abs(lam_ref)))
        self.assertEqual(green.shape[-1], inp_bare["mode"]["param"]["Nmat"])

    def _green_file_with_nmat(self, ci, nmat, name="green.npz"):
        """Write a bare Green npz carrying exactly ``nmat`` Matsubara points."""
        inp = _base_input(os.path.join(self.tmp_path, "gen"), coulomb_inter=ci,
                          bond_channels=True)
        inp["mode"]["param"]["Nmat"] = nmat
        gpath = os.path.join(self.tmp_path, name)
        _bare_green_npz(gpath, inp)
        self.assertEqual(np.load(gpath)["green"].shape[1], nmat)
        return gpath

    def test_bond_green_preflight_runs_once_with_the_files_nmat(self):
        """The FILE's Nmat wins (it defines the bubble's frequency grid), so the
        resource preflight must be run with it -- ONCE -- whether the configured
        Nmat is larger or smaller. Previously the first preflight ran with the
        CONFIGURED Nmat before the file was ever inspected, so an over-large
        configured Nmat could refuse a run that fits comfortably."""
        ci = os.path.join(self.tmp_path, "ci.dat")
        _write_w90(ci, 1, _nn_entries(0.5))
        gpath = self._green_file_with_nmat(ci, 32)

        orig = sc._bond_resource_preflight
        self.addCleanup(setattr, sc, "_bond_resource_preflight", orig)
        for configured in (128, 16):
            seen = []

            def _record(norb, bond_set, Nx, Ny, Nz, nmat, cap_gb, _o=orig):
                seen.append(int(nmat))
                return _o(norb, bond_set, Nx, Ny, Nz, nmat, cap_gb)

            sc._bond_resource_preflight = _record
            inp = _base_input(os.path.join(self.tmp_path, "run%d" % configured),
                              coulomb_inter=ci, bond_channels=True,
                              bond_green=gpath)
            inp["mode"]["param"]["Nmat"] = configured
            sc.calc_eliashberg(inp)
            self.assertEqual(seen, [32], (
                "configured Nmat={}: preflight nmat sequence {} (expected exactly "
                "one call with the file's 32)".format(configured, seen)))

    def _rewrite_green_npz_at_npy_version(self, gpath, version):
        """Re-encode an existing bare-green npz's ``green`` member at the given
        NPY format ``version`` (e.g. ``(3, 0)``).

        ``np.savez`` always picks the lowest header version the data fits, so it
        can never itself produce e.g. format 3.0 for a small plain array; this
        bypasses that to exercise the header-peek's version handling directly, on
        a file ``np.load`` still reads perfectly normally."""
        green_raw = np.load(gpath)["green"]
        with zipfile.ZipFile(gpath, mode="w") as zf:
            with zf.open("green.npy", "w") as fh:
                np.lib.format.write_array(fh, green_raw, version=version)
        return green_raw

    def test_peek_green_npz_nfreq_reads_npy_format_3_0(self):
        """FIX 1 (review round 3): ``_peek_green_npz_nfreq`` supported only NPY
        header versions 1.0/2.0. A ``green`` member written at format 3.0 -- a
        version ``np.load`` accepts fine -- must still make the peek return the
        real ``nfreq``, not ``None`` (which would silently reopen the double
        preflight the header peek exists to close)."""
        ci = os.path.join(self.tmp_path, "ci.dat")
        _write_w90(ci, 1, _nn_entries(0.5))
        gpath = self._green_file_with_nmat(ci, 24, "g24_npy3.npz")
        self._rewrite_green_npz_at_npy_version(gpath, (3, 0))

        # sanity: the file really is format 3.0, and np.load still reads it fine
        with zipfile.ZipFile(gpath) as zf:
            with zf.open("green.npy") as fh:
                self.assertEqual(np.lib.format.read_magic(fh), (3, 0))
        self.assertEqual(np.load(gpath)["green"].shape[1], 24)

        self.assertEqual(sc._peek_green_npz_nfreq(gpath), 24)

    def test_bond_green_preflight_runs_once_with_npy_format_3_0(self):
        """FIX 1 (review round 3), end-to-end: with a format-3.0 ``green``
        member, the resource preflight must still run exactly ONCE, with the
        FILE's Nmat -- not the configured one. ``_load_green_npz`` is patched to
        explode so that if the (buggy) code ever preflights with the wrong
        (configured) Nmat first and only discovers the true Nmat by materialising
        the array, that materialisation-before-the-correct-preflight is caught
        red-handed rather than silently succeeding via a second preflight call."""
        ci = os.path.join(self.tmp_path, "ci.dat")
        _write_w90(ci, 1, _nn_entries(0.5))
        gpath = self._green_file_with_nmat(ci, 24, "g24_npy3.npz")
        self._rewrite_green_npz_at_npy_version(gpath, (3, 0))

        orig_preflight = sc._bond_resource_preflight
        seen = []

        def _record(norb, bond_set, Nx, Ny, Nz, nmat, cap_gb, _o=orig_preflight):
            seen.append(int(nmat))
            return _o(norb, bond_set, Nx, Ny, Nz, nmat, cap_gb)

        def _boom(*args, **kwargs):
            raise AssertionError(
                "_load_green_npz called -- checking that the preflight seen so "
                "far already used the file's Nmat, not the configured one")

        orig_load = sc._load_green_npz
        self.addCleanup(setattr, sc, "_bond_resource_preflight", orig_preflight)
        self.addCleanup(setattr, sc, "_load_green_npz", orig_load)
        sc._bond_resource_preflight = _record
        sc._load_green_npz = _boom

        inp = _base_input(os.path.join(self.tmp_path, "npy3_run"), coulomb_inter=ci,
                          bond_channels=True, bond_green=gpath)
        inp["mode"]["param"]["Nmat"] = 128
        with self.assertRaisesRegex(AssertionError, r"_load_green_npz called"):
            sc.calc_eliashberg(inp)

        self.assertEqual(seen, [24], (
            "preflight nmat sequence {} by the time _load_green_npz was reached "
            "(expected exactly one call, with the file's 24)".format(seen)))

    def test_bond_green_preflight_cap_follows_the_file_not_the_config(self):
        """The memory cap must be judged against the FILE's Nmat alone.

        (a) configured Nmat >> file's: a cap that fits the file must be accepted.
        (b) configured Nmat << file's: a cap that fits only the configured Nmat
            must be REFUSED, and refused BEFORE the Green array is materialised.
        """
        ci = os.path.join(self.tmp_path, "ci.dat")
        _write_w90(ci, 1, _nn_entries(0.5))

        # (a) file nmat=32 (peak 9.86e-5 GB) vs configured 128 (peak 2.21e-4 GB)
        g32 = self._green_file_with_nmat(ci, 32, "g32.npz")
        inp = _base_input(os.path.join(self.tmp_path, "big_cfg"), coulomb_inter=ci,
                          bond_channels=True, bond_green=g32,
                          bond_memory_cap_gb=1.5e-4)
        inp["mode"]["param"]["Nmat"] = 128
        sc.calc_eliashberg(inp)  # must NOT raise
        self.assertTrue(os.path.exists(
            os.path.join(self.tmp_path, "big_cfg", "eigenvalue.dat")))

        # (b) file nmat=128 vs configured 16 (peak 7.81e-5 GB)
        g128 = self._green_file_with_nmat(ci, 128, "g128.npz")
        inp = _base_input(os.path.join(self.tmp_path, "small_cfg"), coulomb_inter=ci,
                          bond_channels=True, bond_green=g128,
                          bond_memory_cap_gb=1.0e-4)
        inp["mode"]["param"]["Nmat"] = 16

        def _boom(*args, **kwargs):
            raise AssertionError(
                "the bond_green array must not be materialised before the "
                "resource preflight has passed")

        orig_load = sc._load_green_npz
        self.addCleanup(setattr, sc, "_load_green_npz", orig_load)
        sc._load_green_npz = _boom
        with self.assertRaisesRegex(ValueError, r"(?i)bond_memory_cap_gb"):
            sc.calc_eliashberg(inp)

    def test_bond_green_provenance_says_which_green_ran(self):
        """The recorded approximation level must distinguish the bare green from
        an externally supplied one -- they are different approximations."""
        ci = os.path.join(self.tmp_path, "ci.dat")
        _write_w90(ci, 1, _nn_entries(0.5))

        d_bare = os.path.join(self.tmp_path, "bare")
        inp_bare = _base_input(d_bare, coulomb_inter=ci, bond_channels=True,
                               solver_mode="eigenvalue", num_eigenvalues=3)
        sc.calc_eliashberg(inp_bare)
        vals = _provenance_keys(os.path.join(d_bare, "eigenvalue.dat"))
        self.assertNotIn("bond_green", vals)
        self.assertIn("bare", vals["approximation"].lower())

        gpath = os.path.join(self.tmp_path, "green.npz")
        _bare_green_npz(gpath, inp_bare)
        d_ext = os.path.join(self.tmp_path, "ext")
        sc.calc_eliashberg(_base_input(d_ext, coulomb_inter=ci,
                                       bond_channels=True, bond_green=gpath,
                                       solver_mode="eigenvalue",
                                       num_eigenvalues=3))
        vals = _provenance_keys(os.path.join(d_ext, "eigenvalue.dat"))
        self.assertEqual(vals["bond_green"], gpath)
        self.assertIn("external", vals["approximation"].lower())
        self.assertNotIn("bare", vals["approximation"].lower())

    def test_bond_green_missing_or_mismatched_file_fails_fast(self):
        """A missing file, a file without a 'green' array, and a green whose grid
        does not match the model all raise instead of being reshaped into silent
        nonsense."""
        ci = os.path.join(self.tmp_path, "ci.dat")
        _write_w90(ci, 1, _nn_entries(0.5))
        inp = _base_input(os.path.join(self.tmp_path, "o"), coulomb_inter=ci,
                          bond_channels=True, solver_mode="eigenvalue",
                          num_eigenvalues=3)

        missing = os.path.join(self.tmp_path, "nope.npz")
        inp["eliashberg"]["bond_green"] = missing
        with self.assertRaisesRegex(FileNotFoundError, r"bond_green"):
            sc.calc_eliashberg(inp)

        empty = os.path.join(self.tmp_path, "empty.npz")
        np.savez(empty, something_else=np.zeros(3))
        inp["eliashberg"]["bond_green"] = empty
        with self.assertRaisesRegex(ValueError, r"(?i)green"):
            sc.calc_eliashberg(inp)

        wrong = os.path.join(self.tmp_path, "wrong.npz")
        np.savez(wrong, green=np.zeros((1, 8, 9, 1, 1), dtype=complex))
        inp["eliashberg"]["bond_green"] = wrong
        with self.assertRaisesRegex(ValueError, r"(?i)does not match the model"):
            sc.calc_eliashberg(inp)

    def test_bond_green_rejects_non_positive_or_odd_nmat(self):
        """The FILE's Nmat overrides the configured one, so it must be validated
        like the configured one: POSITIVE and EVEN.

        ``bond_bubble`` builds its bubble on a CENTERED fermionic grid
        (``iomega_n = (2n+1-nmat) pi/beta``); an odd count puts a fermionic
        frequency exactly at zero, which the bubble then processes without
        complaint -- a silent wrong susceptibility. Reject it up front, naming the
        file and the offending value, and do so BEFORE the resource preflight so
        nothing large is sized (let alone allocated) from a bad grid.
        """
        def _boom(*args, **kwargs):
            raise AssertionError(
                "the bond_green Nmat must be validated BEFORE the resource "
                "preflight / any large allocation")

        orig_preflight = sc._bond_resource_preflight
        orig_load = sc._load_green_npz
        self.addCleanup(setattr, sc, "_bond_resource_preflight", orig_preflight)
        self.addCleanup(setattr, sc, "_load_green_npz", orig_load)

        for nfreq in [0, 1, 7, 129]:
            with self.subTest(nfreq=nfreq):
                ci = os.path.join(self.tmp_path, "ci.dat")
                _write_w90(ci, 1, _nn_entries(0.5))
                bad = os.path.join(self.tmp_path, "bad_nmat_{}.npz".format(nfreq))
                np.savez(bad, green=np.zeros((1, nfreq, 16, 1, 1), dtype=complex))

                sc._bond_resource_preflight = _boom
                sc._load_green_npz = _boom

                inp = _base_input(os.path.join(self.tmp_path, "o_{}".format(nfreq)),
                                  coulomb_inter=ci, bond_channels=True,
                                  solver_mode="eigenvalue", num_eigenvalues=3,
                                  bond_green=bad)
                with self.assertRaises(ValueError) as cm:
                    sc.calc_eliashberg(inp)
                msg = str(cm.exception)
                self.assertIn("bond_green", msg)
                self.assertIn(os.path.basename(bad), msg)
                self.assertIn(str(nfreq), msg)
                self.assertIn("even", msg.lower())

    def test_bond_green_ignored_with_warning_when_flag_off(self):
        """Like every other bond option, ``bond_green`` is ignored with a warning
        when ``bond_channels=false`` (it must not silently swap the green of a
        default-path run)."""
        d = os.path.join(self.tmp_path, "off")
        with _caplog(logging.WARNING) as records:
            sc.calc_eliashberg(_base_input(d, bond_channels=False,
                                           bond_green="/nonexistent/green.npz"))
        msgs = " ".join(r.getMessage() for r in records)
        self.assertIn("bond_green", msgs)

    def test_flex_chi_guard_points_at_bond_green(self):
        """The ``chi0q_mode='flex'`` refusal is about FLEX *chi* ingestion; the
        message must say so and point at ``bond_green`` for a FLEX *green*."""
        ci = os.path.join(self.tmp_path, "ci.dat")
        _write_w90(ci, 1, _nn_entries(1.0))
        inp = _base_input(self.tmp_path, coulomb_inter=ci, bond_channels=True,
                          chi0q_mode="flex", bond_green="green.npz")
        with self.assertRaisesRegex(ValueError, r"bond_green") as cm:
            sc.calc_eliashberg(inp)
        self.assertIn("chi", str(cm.exception).lower())


# ---------------------------------------------------------------------------
# [eliashberg] bond_diagnostics -- opt-in character analysis of the leading
# state (odd-harmonic decomposition + degeneracy clustering + the lambda^pp /
# lambda^fl attribution), emitted into the additive provenance channel.
# ---------------------------------------------------------------------------

def _read_provenance(outdir, name="eigenvalue.dat"):
    """``# key = value`` comment lines of an eigenvalue file, as a dict."""
    out = {}
    with open(os.path.join(outdir, name)) as fh:
        for line in fh:
            if not line.startswith("#"):
                continue
            key, sep, value = line.lstrip("# ").partition("=")
            if sep:
                out[key.strip()] = value.strip()
    return out


def _diagnostics_input(outdir, ci, **eli):
    """Bond run in the odd (triplet) channel -- the sector the odd-harmonic
    basis of ``sc._odd_harmonic_basis`` describes."""
    return _base_input(outdir, coulomb_inter=ci, bond_channels=True,
                       pairing_type="triplet", init_gap="f_x",
                       solver_mode="eigenvalue", num_eigenvalues=6, **eli)


class TestBondDiagnostics(_ApproxTestCase):

    def setUp(self):
        self.tmp_path = tempfile.mkdtemp(prefix="hwave_sc_bond_")
        self.addCleanup(shutil.rmtree, self.tmp_path, ignore_errors=True)

    def test_bond_diagnostics_emits_harmonics_and_clusters(self):
        """``bond_diagnostics = true`` reports the odd-harmonic content of the
        leading (tracked) eigenvector and the degeneracy clustering of the computed
        spectrum, alongside the existing lambda^pp / lambda^fl attribution."""
        ci = os.path.join(self.tmp_path, "ci.dat")
        _write_w90(ci, 1, _nn_entries(0.5))
        d = os.path.join(self.tmp_path, "diag")
        sc.calc_eliashberg(_diagnostics_input(d, ci, bond_diagnostics=True))

        prov = _read_provenance(d)
        self.assertEqual(prov["bond_diagnostics"], "True")

        # (i) harmonic decomposition on the odd basis, one captured fraction per
        #     harmonic, each in [0, 1]
        harm = prov["bond_diagnostics_harmonics"]
        fractions = {}
        for item in harm.split(","):
            name, _, value = item.partition("=")
            fractions[name.strip()] = float(value)
        self.assertEqual(set(fractions), {"p_x", "p_y", "f_x", "f_y"})
        for name, value in fractions.items():
            self.assertTrue(0.0 <= value <= 1.0 + 1e-8, (name, value))
        # a genuinely odd leading state carries non-negligible odd-harmonic weight
        self.assertGreater(sum(fractions.values()), 1.0e-3)

        # (ii) degeneracy clustering of the computed eigenvalues; every eigenpair
        #      is in exactly one cluster
        clusters = eval(prov["bond_diagnostics_eigenvalue_clusters"])  # noqa: S307
        self.assertEqual(sorted(i for c in clusters for i in c), list(range(6)))
        # the reported "leading" cluster is the one containing ROW 0 -- the
        # eigenpair the run actually returns after the channel-parity promotion --
        # not simply the largest-Re-lambda cluster
        leading = eval(prov["bond_diagnostics_leading_cluster"])        # noqa: S307
        self.assertIn(0, leading)
        self.assertIn(leading, clusters)
        self.assertEqual(int(prov["bond_diagnostics_leading_cluster_size"]),
                         len(leading))
        self.assertGreater(float(prov["bond_diagnostics_deg_tol"]), 0.0)

        # (iii) the existing attribution is still there, next to the new keys
        self.assertTrue("lambda_pp" in prov and "lambda_fl" in prov)

        # ... and the numeric rows are UNCHANGED: only comment lines were added, so
        # every downstream reader still sees exactly the same data block
        from hwave import tsweep
        re, im, match = tsweep.parse_leading_eig(
            os.path.join(d, "eigenvalue.dat"))
        self.assertTrue(np.isfinite(re) and np.isfinite(im))

        d_off = os.path.join(self.tmp_path, "off")
        sc.calc_eliashberg(_diagnostics_input(d_off, ci))
        rows_on = _numeric_rows(os.path.join(d, "eigenvalue.dat"))
        rows_off = _numeric_rows(os.path.join(d_off, "eigenvalue.dat"))
        self.assertEqual(len(rows_on), len(rows_off))
        self.assertEqual([len(r) for r in rows_on], [len(r) for r in rows_off])
        self.assertEqual(len(rows_on[-1]), 5)              # index Re Im |ev| match

    def test_bond_diagnostics_leading_cluster_follows_row_zero_not_max_lambda(self):
        """The reported degeneracy is that of the RETURNED state.

        ``_reorder_eigenpairs_by_parity`` promotes the channel-parity eigenpair to
        row 0, so an opposite-parity mode with a LARGER ``Re lambda`` can sit
        further down the file while still being clustered. Reporting the
        largest-``Re lambda`` cluster as "the leading state's degeneracy" would
        then describe a different state than the one the run returns.
        """
        prov = {}
        # rows: 0/1 are the returned (degenerate) pair at 0.02; row 2 is a larger
        # eigenvalue that the parity promotion demoted below them
        evals = np.array([0.02, 0.02 + 1e-9, 0.87])
        kx = np.linspace(0, 2.0 * np.pi, 4, endpoint=False)
        kz = np.linspace(0, 2.0 * np.pi, 1, endpoint=False)
        sc._bond_diagnostics_record(prov, None, evals, None, 1, kx, kx, kz)

        # the largest-Re-lambda cluster is row 2 alone ...
        self.assertEqual(eval(prov["bond_diagnostics_eigenvalue_clusters"])[0], [2])
        # ... but the reported one is the doublet containing row 0 (indices are
        # listed in the library's descending-Re-lambda order within the cluster)
        self.assertEqual(
            sorted(eval(prov["bond_diagnostics_leading_cluster"])), [0, 1])
        self.assertEqual(prov["bond_diagnostics_leading_cluster_size"], 2)

    def test_bond_diagnostics_false_is_byte_identical(self):
        """With the flag off (explicitly false or absent) the bond run's output is
        byte-for-byte what it was before the diagnostics existed.

        Run on the (deterministic) shifted power iteration: ARPACK seeds itself
        from a random start vector, so ``solver_mode='eigenvalue'`` is not
        bit-reproducible run to run and cannot carry a byte-identity claim.
        """
        ci = os.path.join(self.tmp_path, "ci.dat")
        _write_w90(ci, 1, _nn_entries(0.5))
        common = dict(coulomb_inter=ci, bond_channels=True,
                      pairing_type="triplet", init_gap="f_x",
                      solver_mode="iteration", max_iter=20,
                      spectral_shift="auto")
        d_absent = os.path.join(self.tmp_path, "absent")
        d_false = os.path.join(self.tmp_path, "false")
        sc.calc_eliashberg(_base_input(d_absent, **common))
        sc.calc_eliashberg(_base_input(d_false, bond_diagnostics=False, **common))

        _, eig_a, gap_a = _read_outputs(d_absent)
        _, eig_f, gap_f = _read_outputs(d_false)
        self.assertEqual(eig_a, eig_f)
        self.assertEqual(gap_a, gap_f)
        self.assertNotIn(b"bond_diagnostics", eig_a)

    def test_bond_diagnostics_ignored_with_warning_when_flag_off(self):
        """Like every other bond option, ``bond_diagnostics`` is ignored with a
        warning when ``bond_channels=false``."""
        d = os.path.join(self.tmp_path, "off")
        with _caplog(logging.WARNING) as records:
            sc.calc_eliashberg(_base_input(d, bond_channels=False,
                                           bond_diagnostics=True))
        msgs = " ".join(r.getMessage() for r in records)
        self.assertIn("bond_diagnostics", msgs)
        text = open(os.path.join(d, "eigenvalue.dat")).read()
        self.assertNotIn("bond_diagnostics", text)

    def test_bond_diagnostics_on_the_iteration_path_reports_harmonics_only(self):
        """``solver_mode='iteration'`` computes no eigenvalue spectrum, so the
        cluster report is honestly absent (not faked from one number) while the
        harmonic decomposition of the returned gap is still emitted."""
        ci = os.path.join(self.tmp_path, "ci.dat")
        _write_w90(ci, 1, _nn_entries(0.5))
        d = os.path.join(self.tmp_path, "iter")
        sc.calc_eliashberg(_base_input(d, coulomb_inter=ci, bond_channels=True,
                                       pairing_type="triplet", init_gap="f_x",
                                       solver_mode="iteration", max_iter=5,
                                       spectral_shift="auto",
                                       bond_diagnostics=True))
        prov = _read_provenance(d)
        self.assertIn("bond_diagnostics_harmonics", prov)
        self.assertNotIn("bond_diagnostics_eigenvalue_clusters", prov)

    def test_bond_diagnostics_harmonics_match_the_library_helper(self):
        """The emitted numbers ARE ``bond_channels.harmonic_decomposition`` of the
        leading eigenvector in the ``sqrt(GG)`` metric -- the diagnostics wire the
        library helper, they do not re-derive it."""
        from hwave.solver import bond_channels as bc

        ci = os.path.join(self.tmp_path, "ci.dat")
        _write_w90(ci, 1, _nn_entries(0.5))
        d = os.path.join(self.tmp_path, "diag")
        inp = _diagnostics_input(d, ci, bond_diagnostics=True)
        sc.calc_eliashberg(inp)
        prov = _read_provenance(d)
        emitted = {}
        for item in prov["bond_diagnostics_harmonics"].split(","):
            name, _, value = item.partition("=")
            emitted[name.strip()] = float(value)

        # rebuild the SAME quantity: the written gap, the written k-grid, and the
        # pair weight of the same bare green the run built for itself
        Nx = Ny = 4
        gap = np.loadtxt(os.path.join(d, "gap.dat"))
        vec = (gap[:, 3] + 1j * gap[:, 4]).reshape(1, 1, Nx, Ny, 1)
        kx = np.linspace(0, 2.0 * np.pi, Nx, endpoint=False)
        ky = np.linspace(0, 2.0 * np.pi, Ny, endpoint=False)
        kz = np.linspace(0, 2.0 * np.pi, 1, endpoint=False)
        basis = sc._odd_harmonic_basis(1, kx, ky, kz)

        geom_info, hr, interactions = sc._read_interaction_files(inp)
        eps = sc._build_hamiltonian_k(kx, ky, kz, hr, 1)
        ev, evec = sc._calc_eigenvalues(eps)
        beta = 1.0 / inp["mode"]["param"]["T"]
        mu = sc._determine_mu(ev, beta, inp["mode"]["param"]["filling"], 1)
        green = sc._calc_green(ev, evec, mu, beta, inp["mode"]["param"]["Nmat"])
        weight = bc.pair_weight(green, beta, g2_tail=True)

        ref = bc.harmonic_decomposition([vec.ravel()], basis, weight)
        self.assertEqual(set(emitted), set(ref))
        for name, value in ref.items():
            self.assertApprox(emitted[name], value, abs=1e-6)


# ---------------------------------------------------------------------------
# S5 config validation: a malformed value must be REFUSED, never silently
# reinterpreted. A typo'd boolean used to be coerced to false -- silently
# disabling the requested feature AND the five safety guards that come with
# it -- and a malformed integer used to be truncated (int(1.5) -> 1).
# ---------------------------------------------------------------------------

_BOOL_KEYS = ("bond_channels", "bond_diagnostics")
_INT_KEYS = ("bond_max_shells", "bond_precondition_dense_limit")


def _bond_cfg(**kw):
    cfg = {"bond_channels": True}
    cfg.update(kw)
    return cfg


class TestBondConfigValidation(_ApproxTestCase):

    def test_bond_boolean_typo_is_rejected_by_key(self):
        """``bond_channels = "ture"`` must NOT quietly become false."""
        for key in _BOOL_KEYS:
            for bad in ["ture", "flase", "maybe", "", "y e s", 2.5, []]:
                with self.subTest(key=key, bad=bad):
                    cfg = _bond_cfg(**{key: bad})
                    with self.assertRaises(ValueError) as cm:
                        sc._read_bond_config(cfg)
                    self.assertIn(key, str(cm.exception))
                    self.assertTrue(repr(bad) in str(cm.exception)
                                    or str(bad) in str(cm.exception))

    def test_bond_boolean_supported_spellings_still_work(self):
        """Every spelling ``backend.as_bool`` documents keeps working."""
        for spelling, expected in [(True, True), (False, False),
                                   ("true", True), ("True", True), ("yes", True),
                                   ("on", True), ("1", True),
                                   ("false", False), ("False", False), ("no", False),
                                   ("off", False), ("0", False)]:
            with self.subTest(spelling=spelling, expected=expected):
                use_bond = sc._read_bond_config({"bond_channels": spelling})[0]
                self.assertIs(use_bond, expected)
                if expected:
                    diagnostics = sc._read_bond_config(
                        {"bond_channels": True, "bond_diagnostics": spelling})[5]
                    self.assertIs(diagnostics, expected)

    def test_bond_channels_typo_is_rejected_on_the_dynamic_entry_point(self):
        """The dynamic refusal reads the same key, so a typo there must not slip
        a bond_channels request past the guard either."""
        with self.assertRaises(ValueError) as cm:
            sc._reject_bond_channels_dynamic({"bond_channels": "ture"})
        self.assertIn("bond_channels", str(cm.exception))

    def test_bond_integer_options_reject_malformed_values(self):
        """``int(1.5)`` used to silently truncate to 1 and TOML booleans were
        accepted as 0/1; both CHANGE THE MEANING of a malformed config."""
        for key in _INT_KEYS:
            for bad in [1.5, True, False, float("nan"),
                       float("inf"), -1, -2.0, "two", "1.5", []]:
                with self.subTest(key=key, bad=bad):
                    with self.assertRaises(ValueError) as cm:
                        sc._read_bond_config(_bond_cfg(**{key: bad}))
                    self.assertIn(key, str(cm.exception))

    def test_bond_integer_options_accept_integral_values(self):
        for good, expected in [(0, 0), (2, 2), (3.0, 3)]:
            with self.subTest(good=good, expected=expected):
                cfg = _bond_cfg(bond_max_shells=good,
                                bond_precondition_dense_limit=good)
                _, _, max_shells, _, pre_opts, _ = sc._read_bond_config(cfg)
                self.assertTrue(max_shells == expected and isinstance(max_shells, int))
                self.assertEqual(pre_opts["dense_limit"], expected)

    def test_bond_float_options_name_the_key_on_a_bad_value(self):
        """A bare ``float()`` failure never said WHICH [eliashberg] key was
        wrong."""
        keys = ["bond_memory_cap_gb", "bond_precondition_atol",
               "bond_precondition_rtol"]
        for key in keys:
            for bad in ["abc", [], True]:
                with self.subTest(key=key, bad=bad):
                    with self.assertRaises(ValueError) as cm:
                        sc._read_bond_config(_bond_cfg(**{key: bad}))
                    self.assertIn(key, str(cm.exception))

    def test_bond_option_set_to_none_means_unset(self):
        """TOML cannot express a null, so a programmatic ``None`` keeps meaning
        "not specified" (it must not become a validation error)."""
        keys = _BOOL_KEYS + _INT_KEYS + (
            "bond_memory_cap_gb", "bond_precondition_atol",
            "bond_precondition_rtol", "bond_green")
        for key in keys:
            with self.subTest(key=key):
                cfg = _bond_cfg(**{key: None})
                if key == "bond_channels":
                    self.assertIs(sc._read_bond_config(cfg)[0], False)
                    continue
                use_bond, green, max_shells, cap, pre_opts, diag = (
                    sc._read_bond_config(cfg))
                self.assertIs(use_bond, True)
                self.assertTrue(green is None and max_shells is None and diag is False)
                self.assertTrue(
                    cap == sc._BOND_MEMORY_CAP_GB and pre_opts == {})


if __name__ == "__main__":
    unittest.main()
