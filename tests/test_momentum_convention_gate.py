"""The momentum-convention provenance gate (issue #133).

k/q-space NPZ files written since the Fourier-sign alignment carry
momentum_convention = "e_plus_ikR". Loaders reject a mismatching tag
outright, and accept unmarked legacy files ONLY when the stored payload
is elementwise even under q -> -q on the FFT grid (then both
conventions coincide bit-for-bit; every centrosymmetric fixture
qualifies). This mirrors the sc_vertex_version content-conditioned
fail-closed precedent.
"""

import os
import tempfile
import unittest

import numpy as np

import hwave.sc as sc
from hwave.solver.rpa import (MOMENTUM_CONVENTION,
                              validate_momentum_convention)


def _q_even(nx):
    """(nfreq, nvol, 1, 1) payload even under q -> -q."""
    q = np.arange(nx)
    vals = np.cos(2 * np.pi * q / nx)
    return np.tile(vals[None, :, None, None], (2, 1, 1, 1))


def _q_odd(nx):
    q = np.arange(nx)
    vals = np.cos(2 * np.pi * q / nx) + 0.3 * np.sin(2 * np.pi * q / nx)
    return np.tile(vals[None, :, None, None], (2, 1, 1, 1))


class TestValidatorUnit(unittest.TestCase):

    def _npz(self, payload, **extra):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = os.path.join(tmp.name, "x.npz")
        np.savez(path, chi0q=payload, **extra)
        return np.load(path), path

    def test_matching_tag_passes(self):
        data, path = self._npz(_q_odd(4),
                               momentum_convention=MOMENTUM_CONVENTION)
        validate_momentum_convention(data, path, data["chi0q"], 1, (4, 1, 1))

    def test_mismatching_tag_is_rejected(self):
        data, path = self._npz(_q_even(4), momentum_convention="e_minus_ikR")
        with self.assertRaises(ValueError) as cm:
            validate_momentum_convention(data, path, data["chi0q"], 1,
                                         (4, 1, 1))
        self.assertIn("momentum_convention", str(cm.exception))

    def test_unmarked_q_even_payload_is_accepted(self):
        data, path = self._npz(_q_even(4))
        validate_momentum_convention(data, path, data["chi0q"], 1, (4, 1, 1))

    def test_unmarked_q_odd_payload_is_rejected(self):
        data, path = self._npz(_q_odd(4))
        with self.assertRaises(ValueError) as cm:
            validate_momentum_convention(data, path, data["chi0q"], 1,
                                         (4, 1, 1))
        self.assertIn("#133", str(cm.exception))


class TestValidatorHardening(unittest.TestCase):
    """Round-2 attack surface: pair-local tolerance, non-finite content,
    three-axis layouts, and the spin-diag axis collision."""

    def _npz(self, **arrays):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = os.path.join(tmp.name, "x.npz")
        np.savez(path, **arrays)
        return np.load(path), path

    def test_small_odd_channel_is_not_masked_by_a_large_even_one(self):
        """Pair-local scaling: a channel that is order-one q-odd relative
        to ITSELF must reject even when its amplitude sits below the old
        1e-8 * global_max threshold (here 1e-10 of the large channel; the
        residual absolute floor is machine-epsilon-sized, so this is well
        above it)."""
        nx = 4
        q = np.arange(nx)
        big_even = 1.0e6 * np.cos(2 * np.pi * q / nx)
        tiny_odd = 1.0e-4 * np.sin(2 * np.pi * q / nx)
        payload = np.stack([big_even, tiny_odd], axis=-1)[None, :, :]
        data, path = self._npz(chi0q=payload)
        with self.assertRaises(ValueError) as cm:
            validate_momentum_convention(data, path, data["chi0q"], 1,
                                         (nx, 1, 1))
        self.assertIn("#133", str(cm.exception))

    def test_non_finite_unmarked_payload_is_rejected(self):
        payload = _q_even(4)
        payload = payload.astype(float)
        payload[0, 1] = np.nan
        data, path = self._npz(chi0q=payload)
        with self.assertRaises(ValueError) as cm:
            validate_momentum_convention(data, path, data["chi0q"], 1,
                                         (4, 1, 1))
        self.assertIn("non-finite", str(cm.exception))

    def test_float_noise_evenness_is_accepted(self):
        """A genuinely even payload with one element perturbed by
        roundoff-scale ASYMMETRIC noise must still pass."""
        payload = _q_even(4).astype(float)
        payload[0, 1, 0, 0] *= (1.0 + 1.0e-12)   # breaks q-evenness at 1e-12
        data, path = self._npz(chi0q=payload)
        validate_momentum_convention(data, path, data["chi0q"], 1,
                                     (4, 1, 1))

    def test_huge_finite_odd_payload_is_rejected_without_overflow(self):
        """Opposite values near float64.max made |a - rev| and the
        tolerance both infinite (inf > inf is False) before the
        normalized comparison."""
        nx = 4
        payload = np.zeros((1, nx, 1, 1))
        payload[0, 1] = 1.7e308
        payload[0, 3] = -1.7e308
        data, path = self._npz(chi0q=payload)
        with self.assertRaises(ValueError) as cm:
            validate_momentum_convention(data, path, data["chi0q"], 1,
                                         (nx, 1, 1))
        self.assertIn("#133", str(cm.exception))

    def test_zero_paired_with_denormal_is_accepted(self):
        """Deliberate: an exact zero against a minimum denormal sits far
        below the machine-epsilon-scaled absolute floor."""
        nx = 4
        payload = np.ones((1, nx, 1, 1))
        payload[0, 1] = 0.0
        payload[0, 3] = 5e-324
        data, path = self._npz(chi0q=payload)
        validate_momentum_convention(data, path, data["chi0q"], 1,
                                     (nx, 1, 1))

    def test_6d_ref_layout_with_norb_equal_nvol_is_gated(self):
        """Round-3 reproduction: ref (norb, norb, Nx, Ny, Nz, nfreq) with
        norb == nvol was misrouted onto the orbital axis and bypassed."""
        nx = 4
        q = np.arange(nx)
        odd = np.cos(2 * np.pi * q / nx) + 0.3 * np.sin(2 * np.pi * q / nx)
        payload = np.zeros((nx, nx, nx, 1, 1, 2))
        payload[:, :, :, 0, 0, :] = odd[None, None, :, None]
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        np.savez(os.path.join(d, "chi0q.npz"), chi0q=payload)
        inp = {"mode": {"param": {"T": 1.0, "Nmat": 2,
                                  "CellShape": [nx, 1, 1]}},
               "file": {"input": {"path_to_flex_output": d},
                        "output": {"path_to_output": d}},
               "eliashberg": {}}
        with self.assertRaises(ValueError) as cm:
            sc._load_chi0q(inp)
        self.assertIn("#133", str(cm.exception))

    def test_ambiguous_unmarked_6d_layout_fails_closed(self):
        """A shape matching BOTH raw and ref patterns cannot have its
        momentum axes identified; unmarked such files must be refused,
        marked ones accepted on the tag."""
        # grid (2,2,2): shape (8, 8, 2, 2, 2, 2) is ref(norb=8) AND
        # raw(nfreq=8, nvol=8, norb=2)
        payload = np.zeros((8, 8, 2, 2, 2, 2))
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        np.savez(os.path.join(d, "chi0q.npz"), chi0q=payload)
        inp = {"mode": {"param": {"T": 1.0, "Nmat": 8,
                                  "CellShape": [2, 2, 2]}},
               "file": {"input": {"path_to_flex_output": d},
                        "output": {"path_to_output": d}},
               "eliashberg": {}}
        with self.assertRaises(ValueError) as cm:
            sc._load_chi0q(inp)
        self.assertIn("BOTH", str(cm.exception))
        np.savez(os.path.join(d, "chi0q.npz"), chi0q=payload,
                 momentum_convention=MOMENTUM_CONVENTION)
        sc._load_chi0q(inp)   # tag decides; must not raise here

    def test_non_cubic_permuted_and_negative_axes(self):
        """The 3-axis path must respect axis identity on a non-cubic grid
        (4, 2, 1), in permuted order and with negative axis indices."""
        nx, ny = 4, 2
        qx = np.arange(nx)
        odd_x = (np.cos(2 * np.pi * qx / nx)
                 + 0.3 * np.sin(2 * np.pi * qx / nx))
        payload = np.zeros((1, 1, nx, ny, 1, 2))
        payload[0, 0, :, :, 0, :] = odd_x[:, None, None]
        data, path = self._npz(chi0q=payload)
        for axes in ((2, 3, 4), (-4, -3, -2)):
            with self.subTest(axes=axes):
                with self.assertRaises(ValueError):
                    validate_momentum_convention(
                        data, path, data["chi0q"], axes, (nx, ny, 1))

    def test_three_axis_layout_is_gated(self):
        """Already-expanded reference layouts carry (qx, qy, qz) on three
        separate axes."""
        nx = 4
        q = np.arange(nx)
        odd = (np.cos(2 * np.pi * q / nx)
               + 0.3 * np.sin(2 * np.pi * q / nx))
        payload = np.zeros((1, 1, nx, 1, 1, 2))
        payload[0, 0, :, 0, 0, :] = odd[:, None]
        data, path = self._npz(chi0q=payload)
        with self.assertRaises(ValueError):
            validate_momentum_convention(data, path, data["chi0q"],
                                         (2, 3, 4), (nx, 1, 1))

    def test_spin_diag_collision_nfreq_equals_nvol(self):
        """(2, nfreq, nvol, ...) with nfreq == nvol: the loader must probe
        the q axis (2), not the frequency axis -- a size search picked
        axis 1 and slipped a q-odd payload through."""
        nx = 4
        q = np.arange(nx)
        odd = (np.cos(2 * np.pi * q / nx)
               + 0.3 * np.sin(2 * np.pi * q / nx))
        # q-odd along axis 2, constant along the (same-sized) axis 1
        payload = np.tile(odd[None, None, :, None, None], (2, nx, 1, 1, 1))
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        np.savez(os.path.join(d, "chi0q.npz"), chi0q=payload)
        inp = {"mode": {"param": {"T": 1.0, "Nmat": nx,
                                  "CellShape": [nx, 1, 1]}},
               "file": {"input": {"path_to_flex_output": d},
                        "output": {"path_to_output": d}},
               "eliashberg": {}}
        with self.assertRaises(ValueError) as cm:
            sc._load_chi0q(inp)
        self.assertIn("#133", str(cm.exception))


class TestRound4Hardening(unittest.TestCase):
    """Round-4 attack surface: the removed 6D fallback, denormal-scale
    payloads, and the IR-native bypass."""

    def _sc_load(self, payload, cell, **extra):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        np.savez(os.path.join(d, "chi0q.npz"), chi0q=payload, **extra)
        inp = {"mode": {"param": {"T": 1.0, "Nmat": 8, "CellShape": cell}},
               "file": {"input": {"path_to_flex_output": d},
                        "output": {"path_to_output": d}},
               "eliashberg": {}}
        return sc._load_chi0q(inp)

    def test_mismatched_grid_6d_layout_fails_closed(self):
        """Neither the raw nor the ref pattern matches: the removed
        partial fallback previously gated the wrong axis here."""
        payload = np.zeros((8, 8, 4, 2, 1, 2))
        q2 = np.arange(4)
        payload[0, 0, :, 0, 0, 0] = (np.cos(2 * np.pi * q2 / 4)
                                     + 0.3 * np.sin(2 * np.pi * q2 / 4))
        with self.assertRaises(ValueError) as cm:
            self._sc_load(payload, [2, 2, 2])
        self.assertIn("neither", str(cm.exception))

    def test_whole_payload_at_denormal_scale_is_accepted(self):
        """Global maximum itself denormal: clamping the normalization
        scale keeps such payloads below the machine floor instead of
        blowing them up to O(1) and rejecting."""
        nx = 4
        payload = np.zeros((1, nx, 1, 1))
        payload[0, 1] = np.nextafter(0.0, 1.0)
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = os.path.join(tmp.name, "x.npz")
        np.savez(path, chi0q=payload)
        data = np.load(path)
        validate_momentum_convention(data, path, data["chi0q"], 1,
                                     (nx, 1, 1))

    def test_marked_unknown_6d_layout_still_fails_closed(self):
        """A marker establishes the Fourier sign, never the layout: a
        marked file whose 6D shape matches neither pattern was previously
        accepted and silently reshaped downstream."""
        payload = np.zeros((4, 8, 1, 2, 2, 4))
        with self.assertRaises(ValueError) as cm:
            self._sc_load(payload, [2, 2, 2],
                          momentum_convention=MOMENTUM_CONVENTION)
        self.assertIn("neither", str(cm.exception))

    def test_malformed_marker_arrays_are_rejected(self):
        """Exactly one marker value: a multi-element array previously
        authorized via its first element; an empty one crashed with
        IndexError."""
        nx = 4
        q = np.arange(nx)
        odd = (np.cos(2 * np.pi * q / nx)
               + 0.3 * np.sin(2 * np.pi * q / nx))
        payload = np.tile(odd[None, :, None, None], (2, 1, 1, 1))
        for label, marker in (
                ("multi", np.array(["e_plus_ikR", "e_minus_ikR"])),
                ("empty", np.array([], dtype=str))):
            with self.subTest(marker=label):
                tmp = tempfile.TemporaryDirectory()
                self.addCleanup(tmp.cleanup)
                path = os.path.join(tmp.name, "x.npz")
                np.savez(path, chi0q=payload, momentum_convention=marker)
                data = np.load(path)
                with self.assertRaises(ValueError) as cm:
                    validate_momentum_convention(data, path,
                                                 data["chi0q"], 1,
                                                 (nx, 1, 1))
                self.assertIn("single value", str(cm.exception))

    def test_norb_layout_gate_rejects_malformed_shapes(self):
        """With the orbital count supplied (as the production entry does),
        malformed 4D/6D/8D shapes are rejected outright, marked or not
        (round-6 review: a malformed 8D file previously passed both)."""
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        for label, shape in (("8d", (4, 1, 2, 2, 2, 2, 2, 4)),
                             ("4d", (4, 8, 2, 3)),
                             ("6d", (4, 8, 2, 2, 2, 3))):
            for marked in (False, True):
                with self.subTest(shape=label, marked=marked):
                    d = tempfile.mkdtemp(dir=tmp.name)
                    extra = ({"momentum_convention": MOMENTUM_CONVENTION}
                             if marked else {})
                    np.savez(os.path.join(d, "chi0q.npz"),
                             chi0q=np.zeros(shape), **extra)
                    inp = {"mode": {"param": {"T": 1.0, "Nmat": 4,
                                              "CellShape": [2, 2, 2]}},
                           "file": {"input": {"path_to_flex_output": d},
                                    "output": {"path_to_output": d}},
                           "eliashberg": {}}
                    with self.assertRaises(ValueError) as cm:
                        sc._load_chi0q(inp, norb=2)
                    self.assertIn("no supported layout",
                                  str(cm.exception))

    def test_norb_gate_rejects_malformed_spin_diag_shapes(self):
        """5D/7D spin-diag layouts are checked over their FULL trailing
        shapes (round-7: partial slices accepted malformed arrays)."""
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        for label, shape in (("5d", (2, 4, 8, 2, 3)),
                             ("7d", (2, 4, 8, 2, 2, 2, 3))):
            for marked in (False, True):
                with self.subTest(shape=label, marked=marked):
                    d = tempfile.mkdtemp(dir=tmp.name)
                    extra = ({"momentum_convention": MOMENTUM_CONVENTION}
                             if marked else {})
                    np.savez(os.path.join(d, "chi0q.npz"),
                             chi0q=np.zeros(shape), **extra)
                    inp = {"mode": {"param": {"T": 1.0, "Nmat": 4,
                                              "CellShape": [2, 2, 2]}},
                           "file": {"input": {"path_to_flex_output": d},
                                    "output": {"path_to_output": d}},
                           "eliashberg": {}}
                    with self.assertRaises(ValueError) as cm:
                        sc._load_chi0q(inp, norb=2)
                    self.assertIn("no supported layout",
                                  str(cm.exception))

    def test_norb_resolves_the_structural_ambiguity(self):
        """(8, 8, 2, 2, 2, 2) on a 2x2x2 grid is raw for norb = 2 and ref
        for norb = 8; with norb supplied the routing must resolve it (the
        no-norb path keeps failing closed, pinned elsewhere)."""
        payload = np.zeros((8, 8, 2, 2, 2, 2))
        for no in (2, 8):
            with self.subTest(norb=no):
                tmp = tempfile.TemporaryDirectory()
                self.addCleanup(tmp.cleanup)
                d = tmp.name
                np.savez(os.path.join(d, "chi0q.npz"), chi0q=payload,
                         momentum_convention=MOMENTUM_CONVENTION)
                inp = {"mode": {"param": {"T": 1.0, "Nmat": 8,
                                          "CellShape": [2, 2, 2]}},
                       "file": {"input": {"path_to_flex_output": d},
                                "output": {"path_to_output": d}},
                       "eliashberg": {}}
                sc._load_chi0q(inp, norb=no)   # must not raise

    def test_norb_routing_selects_the_right_axes_on_a_3cubed_grid(self):
        """Axis-selection regression on a 3x3x3 grid (on a 2-point grid
        every momentum is self-inverse, so evenness cannot distinguish
        the axes): an UNMARKED payload q-odd along the raw volume axis
        must reject under the raw interpretation and accept under the
        ref interpretation of the same ambiguous shape."""
        n = 3
        nvol = n ** 3
        # ambiguous shape (27, 27, 3, 3, 3, 2): raw for norb=3? no --
        # raw needs shape (nfreq, nvol, no x4) = (nfreq, 27, 3, 3, 3, 3);
        # use ref (norb, norb, 3, 3, 3, nfreq) with norb=27 vs raw with
        # norb=3: shape (27, 27, 3, 3, 3, 3) satisfies BOTH.
        shape = (nvol, nvol, n, n, n, n)
        q = np.arange(n)
        odd = np.cos(2 * np.pi * q / n) + 0.3 * np.sin(2 * np.pi * q / n)
        payload = np.zeros(shape)
        # odd along axis 2 (the ref qx axis); constant along axis 1 (the
        # raw flattened volume axis) -> even under the raw interpretation
        payload += odd[None, None, :, None, None, None]
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        np.savez(os.path.join(d, "chi0q.npz"), chi0q=payload)
        inp = {"mode": {"param": {"T": 1.0, "Nmat": nvol,
                                  "CellShape": [n, n, n]}},
               "file": {"input": {"path_to_flex_output": d},
                        "output": {"path_to_output": d}},
               "eliashberg": {}}
        # ref interpretation (norb = 27): gate probes axes (2,3,4),
        # payload is q-odd there -> reject
        with self.assertRaises(ValueError) as cm:
            sc._load_chi0q(inp, norb=nvol)
        self.assertIn("#133", str(cm.exception))
        # raw interpretation (norb = 3): gate probes axis 1, payload is
        # constant there -> accept
        sc._load_chi0q(inp, norb=n)

    def test_huge_complex_odd_payload_is_rejected(self):
        """Finite complex values whose MAGNITUDE overflows: |x+iy| = inf
        for x, y ~ 1.3e308, and an infinite scale previously normalized
        the asymmetry away."""
        nx = 4
        z = 1.3e308 + 1.3e308j
        payload = np.zeros((1, nx, 1, 1), dtype=complex)
        payload[0, 1] = z
        payload[0, 3] = -z
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = os.path.join(tmp.name, "x.npz")
        np.savez(path, chi0q=payload)
        data = np.load(path)
        with self.assertRaises(ValueError) as cm:
            validate_momentum_convention(data, path, data["chi0q"], 1,
                                         (nx, 1, 1))
        self.assertIn("#133", str(cm.exception))

    def _ir_meta(self):
        return dict(matsubara_basis="ir",
                    frequency_grid="sparse_ir_nodes",
                    ir_freq_n=np.arange(3))

    def test_ir_native_q_odd_chi_is_rejected(self):
        """IR-native chi files were wrongly exempted from the gate; their
        q volume also sits on axis 1."""
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        nx = 4
        q = np.arange(nx)
        odd = (np.cos(2 * np.pi * q / nx)
               + 0.3 * np.sin(2 * np.pi * q / nx))
        chi = np.tile(odd[None, :, None, None], (3, 1, 2, 2))
        np.savez(os.path.join(d, "chiq_s.npz"), chiq_s=chi,
                 **self._ir_meta())
        np.savez(os.path.join(d, "chiq_c.npz"), chiq_c=chi,
                 **self._ir_meta())
        inp = {"mode": {"param": {"T": 1.0, "CellShape": [nx, 1, 1]}},
               "file": {"output": {"path_to_output": d}},
               "eliashberg": {}}
        with self.assertRaises(ValueError) as cm:
            sc._read_flex_chi_raw(inp, allow_ir=True)
        self.assertIn("#133", str(cm.exception))

    def test_ir_native_q_odd_green_is_rejected(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        nx = 4
        q = np.arange(nx)
        prof = 1.0 + 0.3 * np.sin(2 * np.pi * q / nx)
        green = np.ones((1, 3, nx, 1, 1), dtype=complex) *             prof[None, None, :, None, None]
        np.savez(os.path.join(d, "green.npz"), green=green,
                 **self._ir_meta())
        inp = {"mode": {"param": {"T": 2.0}},
               "file": {"output": {"path_to_output": d}},
               "eliashberg": {}}
        with self.assertRaises(ValueError) as cm:
            sc._load_flex_green(inp, 1, nx, 1, 1, allow_ir=True)
        self.assertIn("#133", str(cm.exception))


class TestDynamicGapSeedGate(unittest.TestCase):

    def test_unmarked_k_odd_seed_is_rejected(self):
        from hwave.solver import eliashberg_dynamic as ed
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        nx, nmat = 4, 4
        q = np.arange(nx)
        prof = 1.0 + 0.3 * np.sin(2 * np.pi * q / nx)
        gap = np.ones((1, 1, nx, 1, 1, nmat), dtype=complex) *             prof[None, None, :, None, None, None]
        path = os.path.join(tmp.name, "seed.npz")
        np.savez(path, gap=gap)
        with self.assertRaises(ValueError) as cm:
            ed._load_seed_gap({"seed_eigenvector": path}, gap.shape,
                              False, None, nmat)
        self.assertIn("#133", str(cm.exception))


class TestLoaderGates(unittest.TestCase):
    """Each production loader must reach the gate."""

    def test_flex_green_loader_rejects_unmarked_q_odd_file(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        nx = 4
        q = np.arange(nx)
        prof = 1.0 + 0.3 * np.sin(2 * np.pi * q / nx)
        green = np.ones((1, 4, nx, 1, 1), dtype=complex) * \
            prof[None, None, :, None, None]
        np.savez(os.path.join(d, "green.npz"), green=green, beta=0.5)
        inp = {"mode": {"param": {"T": 2.0}},
               "file": {"output": {"path_to_output": d}},
               "eliashberg": {}}
        with self.assertRaises(ValueError) as cm:
            sc._load_flex_green(inp, 1, nx, 1, 1)
        self.assertIn("#133", str(cm.exception))

    def test_flex_green_loader_accepts_marked_file(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        nx = 4
        q = np.arange(nx)
        prof = 1.0 + 0.3 * np.sin(2 * np.pi * q / nx)
        green = np.ones((1, 4, nx, 1, 1), dtype=complex) * \
            prof[None, None, :, None, None]
        np.savez(os.path.join(d, "green.npz"), green=green, beta=0.5,
                 momentum_convention=MOMENTUM_CONVENTION)
        inp = {"mode": {"param": {"T": 2.0}},
               "file": {"output": {"path_to_output": d}},
               "eliashberg": {}}
        green_loaded = sc._load_flex_green(inp, 1, nx, 1, 1)
        self.assertIsNotNone(green_loaded)

    def test_sc_chi0q_loader_rejects_unmarked_q_odd_file(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        nx = 4
        chi0q = _q_odd(nx)
        np.savez(os.path.join(d, "chi0q.npz"), chi0q=chi0q)
        inp = {"mode": {"param": {"T": 1.0, "Nmat": 2,
                                  "CellShape": [nx, 1, 1]}},
               "file": {"input": {"path_to_flex_output": d},
                        "output": {"path_to_output": d}},
               "eliashberg": {}}
        with self.assertRaises(ValueError) as cm:
            sc._load_chi0q(inp)
        self.assertIn("#133", str(cm.exception))


class TestRound9Hardening(unittest.TestCase):
    """Round-9 audit: tolerance boundary pins, marked-non-finite
    rejection, call-site integration, independent chi gates, 8D routing,
    writer marker values, bounded memory, manifest v2 refusal."""

    def _npz(self, **arrays):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = os.path.join(tmp.name, "x.npz")
        np.savez(path, **arrays)
        return np.load(path), path

    def test_pair_tolerance_boundary_is_pinned(self):
        """1e-8 pair-relative: 3e-9 asymmetry accepted, 3e-8 rejected
        (pins the coefficient itself; zero or 1e-2 would flip one)."""
        nx = 4
        base = np.cos(2 * np.pi * np.arange(nx) / nx) + 2.0
        for eps, ok in ((3.0e-9, True), (3.0e-8, False)):
            with self.subTest(eps=eps):
                payload = np.tile(base[None, :, None, None], (2, 1, 1, 1))
                payload = payload.copy()
                payload[0, 1] *= (1.0 + eps)
                data, path = self._npz(chi0q=payload)
                if ok:
                    validate_momentum_convention(data, path,
                                                 data["chi0q"], 1,
                                                 (nx, 1, 1))
                else:
                    with self.assertRaises(ValueError):
                        validate_momentum_convention(data, path,
                                                     data["chi0q"], 1,
                                                     (nx, 1, 1))

    def test_marked_non_finite_payload_is_rejected(self):
        """A valid marker must not authorize NaN/Inf content (round 9:
        the tag early-return skipped the finiteness scan)."""
        payload = _q_even(4).astype(float)
        payload[0, 1] = np.nan
        data, path = self._npz(chi0q=payload,
                               momentum_convention=MOMENTUM_CONVENTION)
        with self.assertRaises(ValueError) as cm:
            validate_momentum_convention(data, path, data["chi0q"], 1,
                                         (4, 1, 1))
        self.assertIn("non-finite", str(cm.exception))

    def test_rpa_chi0q_init_call_site_is_gated(self):
        """Integration through the RPA solver's own read_init: removing
        the _read_chi0q gate call must fail this test."""
        import sys as _sys
        _sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import test_rpa_fourier_sign as fs
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        fs._write_chain(tmp.name, so=False, norb_phys=1)
        solver, gi = fs._solver(tmp.name, so=False, with_v=False, nmat=8)
        out = os.path.join(tmp.name, "out")
        os.makedirs(out, exist_ok=True)
        solver.solve(gi, out)
        # overwrite the produced chi0q with an unmarked q-odd file
        nx = fs.LX
        q = np.arange(nx)
        odd = (np.cos(2 * np.pi * q / nx)
               + 0.3 * np.sin(2 * np.pi * q / nx))
        chi = np.tile(odd[None, :, None, None], (8, 1, 1, 1))
        np.savez(os.path.join(out, "chi0q.npz"), chi0q=chi)
        with self.assertRaises(ValueError) as cm:
            solver.read_init({"path_to_input": out,
                              "chi0q_init": "chi0q.npz"})
        self.assertIn("#133", str(cm.exception))

    def test_chi_s_and_chi_c_gates_fire_independently(self):
        """One valid file must not shield the other: reject when only
        chi_c is unmarked q-odd."""
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        nx = 4
        q = np.arange(nx)
        even = np.tile(np.cos(2 * np.pi * q / nx)[None, :, None, None],
                       (2, 1, 2, 2))
        odd = np.tile((np.cos(2 * np.pi * q / nx)
                       + 0.3 * np.sin(2 * np.pi * q / nx))
                      [None, :, None, None], (2, 1, 2, 2))
        np.savez(os.path.join(d, "chiq_s.npz"), chiq_s=even,
                 momentum_convention=MOMENTUM_CONVENTION)
        np.savez(os.path.join(d, "chiq_c.npz"), chiq_c=odd)
        inp = {"mode": {"param": {"T": 1.0, "CellShape": [nx, 1, 1]}},
               "file": {"output": {"path_to_output": d}},
               "eliashberg": {}}
        with self.assertRaises(ValueError) as cm:
            sc._read_flex_chi_raw(inp)
        self.assertIn("chiq_c", str(cm.exception))

    def test_valid_8d_ref_layout_with_q_odd_content_is_rejected(self):
        n = 3
        q = np.arange(n)
        odd = np.cos(2 * np.pi * q / n) + 0.3 * np.sin(2 * np.pi * q / n)
        payload = np.zeros((2, 2, 2, 2, n, 1, 1, 2))
        payload += odd[None, None, None, None, :, None, None, None]
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        np.savez(os.path.join(d, "chi0q.npz"), chi0q=payload)
        inp = {"mode": {"param": {"T": 1.0, "Nmat": 2,
                                  "CellShape": [n, 1, 1]}},
               "file": {"input": {"path_to_flex_output": d},
                        "output": {"path_to_output": d}},
               "eliashberg": {}}
        with self.assertRaises(ValueError) as cm:
            sc._load_chi0q(inp, norb=2)
        self.assertIn("#133", str(cm.exception))

    def test_rpa_writer_emits_the_exact_marker_value(self):
        import sys as _sys
        _sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import test_rpa_fourier_sign as fs
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        fs._write_chain(tmp.name, so=False, norb_phys=1)
        solver, gi = fs._solver(tmp.name, so=False, with_v=False, nmat=8)
        out = os.path.join(tmp.name, "out")
        os.makedirs(out, exist_ok=True)
        solver.solve(gi, out)
        solver.save_results({"path_to_output": out, "chi0q": "c.npz"}, gi)
        with np.load(os.path.join(out, "c.npz")) as data:
            self.assertEqual(
                str(np.asarray(data["momentum_convention"]).ravel()[0]),
                MOMENTUM_CONVENTION)

    def test_evenness_scan_memory_is_bounded(self):
        """The legacy scan must not materialize whole-tensor
        temporaries (round 9: ~5.5x the payload). Peak extra allocation
        for a 16 MiB q-even payload must stay well below the payload
        size."""
        import tracemalloc
        nvol = 8192
        base = np.cos(2 * np.pi * np.arange(nvol) / nvol)
        payload = np.tile(base[None, :, None, None],
                          (128, 1, 2, 2)).astype(complex)   # ~67 MB
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = os.path.join(tmp.name, "x.npz")
        np.savez(path, chi0q=payload)
        data = np.load(path)
        arr = data["chi0q"]
        tracemalloc.start()
        tracemalloc.reset_peak()
        validate_momentum_convention(data, path, arr, 1, (nvol, 1, 1))
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        # blocks are absolute-bounded (~30 MB here); the pre-fix
        # whole-tensor pipeline peaked at ~5.5x the payload (~370 MB)
        self.assertLess(peak, int(arr.nbytes * 0.6))

    def test_manifest_v2_upgrade_is_refused(self):
        import hwave.tsweep as ts
        old_manifest = {"version": 2, "ladder": [1.0], "fingerprint": "x"}
        with self.assertRaises(ValueError) as cm:
            ts._validate_resume(old_manifest, [1.0], "x")
        self.assertIn("version", str(cm.exception))


class TestExternFourierSign(unittest.TestCase):
    """Extern's R -> k build shares the documented sign."""

    def test_extern_q_is_documented_sign(self):
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.rpa as solver_rpa
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        nx = 4
        t = 0.7 + 0.3j
        h = 0.2 + 0.4j
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n1\n0.0 0.0 0.0\n")
        with open(os.path.join(d, "transfer.dat"), "w") as f:
            f.write("hdr\n1\n2\n1 1\n")
            f.write(" 1 0 0 1 1 %.12f %.12f\n" % (t.real, t.imag))
            f.write("-1 0 0 1 1 %.12f %.12f\n"
                    % (np.conj(t).real, np.conj(t).imag))
        with open(os.path.join(d, "extern.dat"), "w") as f:
            f.write("hdr\n1\n2\n1 1\n")
            f.write(" 1 0 0 1 1 %.12f %.12f\n" % (h.real, h.imag))
            f.write("-1 0 0 1 1 %.12f %.12f\n"
                    % (np.conj(h).real, np.conj(h).imag))
        io = read_input_k.QLMSkInput(
            {"path_to_input": d,
             "interaction": {"path_to_input": d, "Geometry": "geom.dat",
                             "Transfer": "transfer.dat",
                             "Extern": "extern.dat"}})
        info_mode = {"mode": "RPA",
                     "param": {"T": 1.0, "mu": 0.0, "CellShape": [nx, 1, 1],
                               "SubShape": [1, 1, 1], "Nmat": 8},
                     "calc_scheme": "reduced"}
        solver = solver_rpa.RPA(io.get_param("ham"), {}, info_mode)
        got = np.asarray(solver.ham_info.ham_extern_q).reshape(nx)
        want = np.array([h * np.exp(2j * np.pi * k / nx)
                         + np.conj(h) * np.exp(-2j * np.pi * k / nx)
                         for k in range(nx)])
        np.testing.assert_allclose(got, want, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
