"""Regression tests: hwave_sc must locate the static (zero bosonic frequency)
chi0q slice via the freq_index metadata, not blindly at nmat//2.

RPA's matsubara_frequency option can restrict the frequency axis of
chi0q.npz; freq_index records the original indices (zero frequency sits at
original index Nmat//2).  Loading such a file and slicing at nfreq//2
silently uses a finite-frequency chi0 as the static limit.
"""
import os
import tempfile
import unittest

import numpy as np

from hwave.sc import (
    _static_freq_position,
    _load_chi0q,
    _compute_vertices_simple,
)


class TestStaticFreqPosition(unittest.TestCase):
    def test_full_grid(self):
        freq_index = np.arange(8)
        self.assertEqual(_static_freq_position(freq_index, 8, 8, "f"), 4)

    def test_full_grid_with_file_nmat_ignores_mismatched_config(self):
        # a full-grid file with nmat metadata resolves regardless of the
        # hwave_sc config Nmat (e.g. RPA Nmat=2048, Eliashberg Nmat=1024);
        # without nmat metadata this case is indistinguishable from a
        # [0,K] restriction and raises (see
        # test_indices_exceeding_config_raises)
        freq_index = np.arange(2048)
        self.assertEqual(
            _static_freq_position(freq_index, 2048, 1024, "f",
                                  file_nmat=2048), 1024)

    def test_restricted_range_containing_zero(self):
        # original Nmat=8, restricted to indices 2..6: zero (index 4) at pos 2
        freq_index = np.arange(2, 7)
        self.assertEqual(_static_freq_position(freq_index, 5, 8, "f"), 2)

    def test_file_nmat_takes_precedence_over_config(self):
        # file records its own nmat=8; the (wrong) config Nmat is ignored
        freq_index = np.arange(2, 7)
        self.assertEqual(
            _static_freq_position(freq_index, 5, 1024, "f", file_nmat=8), 2)

    def test_restricted_range_missing_zero_raises(self):
        freq_index = np.arange(0, 4)  # zero frequency (index 4) not included
        with self.assertRaises(ValueError):
            _static_freq_position(freq_index, 4, 8, "f", file_nmat=8)

    def test_indices_exceeding_config_raises(self):
        # freq_index = 0..2047 with config Nmat=1024: not a restriction of
        # the config grid, but ALSO not provably a full grid (it could be a
        # [0,2047] restriction of an even larger run) -- the zero-frequency
        # position is unknowable without nmat metadata, so it must fail
        # with instructions instead of silently assuming a full grid
        freq_index = np.arange(0, 2048)
        with self.assertRaises(ValueError):
            _static_freq_position(freq_index, 2048, 1024, "f")

    def test_non_zero_based_range_exceeding_config_raises(self):
        # legacy restricted range [512,1536] with config Nmat=1024: the
        # indices reach past the config grid, so the config-based lookup is
        # disproven; silently returning the position of 512 (=1024//2, at
        # position 0!) would slice a finite frequency
        freq_index = np.arange(512, 1537)
        with self.assertRaises(ValueError):
            _static_freq_position(freq_index, 1025, 1024, "f")

    def test_zero_based_range_without_provable_reading_raises(self):
        # legacy file (no nmat metadata) with freq_index = 0..3 and config
        # Nmat=8: EITHER a full 4-grid (zero at 2) OR a [0,3] restriction of
        # the 8-grid (which contains NO zero frequency) -- silently picking
        # the full-grid reading would return a finite frequency for the
        # restricted file, so it must fail with instructions
        freq_index = np.arange(0, 4)
        with self.assertRaises(ValueError):
            _static_freq_position(freq_index, 4, 8, "f")

    def test_ambiguous_zero_based_range_raises(self):
        # legacy file with freq_index = 0..600 and config Nmat=1024: this is
        # EITHER a full 601-grid (zero at 300) OR a [0,600] restriction of
        # the 1024-grid (zero at 512) -- silently picking either one is
        # wrong for the other, so it must fail with instructions
        freq_index = np.arange(0, 601)
        with self.assertRaises(ValueError):
            _static_freq_position(freq_index, 601, 1024, "f")

    def test_empty_freq_index_raises_value_error(self):
        # matsubara_frequency="none" output: must be a clean ValueError,
        # not an IndexError from the error-message formatting
        with self.assertRaises(ValueError):
            _static_freq_position(np.array([], dtype=int), 0, 8, "f")

    def test_empty_freq_index_with_data_raises(self):
        # empty metadata beats the length-mismatch fallback: it must fail
        # cleanly even when a data axis exists
        with self.assertRaises(ValueError):
            _static_freq_position(np.array([], dtype=int), 8, 8, "f")

    def test_legacy_file_without_freq_index_delegates_to_data_axis(self):
        # no metadata: the caller must slice the center of the ACTUAL data
        # axis (which only the caller can identify, e.g. 6D ref format)
        self.assertIsNone(_static_freq_position(None, 8, 8, "f"))

    def test_inconsistent_length_delegates_to_data_axis(self):
        # legacy FLEX chi0q files store the FULL grid but a restricted
        # freq_index; the DATA axis is authoritative, so the caller slices
        # its center (with a warning) exactly as before the metadata
        # handling was introduced
        self.assertIsNone(_static_freq_position(np.arange(5), 8, 8, "f"))


class TestLoadChi0qStaticIndex(unittest.TestCase):
    def _write_and_load(self, nmat_config, freq_index, nfreq, file_nmat=None,
                        shape=None):
        with tempfile.TemporaryDirectory() as tmp:
            if shape is None:
                shape = (nfreq, 4, 1, 1)
            chi0q = np.zeros(shape, dtype=np.complex128)
            kwargs = {"chi0q": chi0q, "index_convention": "spin_block"}
            if freq_index is not None:
                kwargs["freq_index"] = np.asarray(freq_index)
            if file_nmat is not None:
                kwargs["nmat"] = file_nmat
            np.savez(os.path.join(tmp, "chi0q.npz"), **kwargs)
            input_dict = {
                "mode": {"param": {"Nmat": nmat_config}},
                "file": {"output": {"path_to_output": tmp}},
            }
            return _load_chi0q(input_dict)

    def test_restricted_file_returns_true_static_index(self):
        chi0q, static_index = self._write_and_load(8, np.arange(2, 7), 5)
        self.assertEqual(static_index, 2)

    def test_full_file_returns_center(self):
        chi0q, static_index = self._write_and_load(8, np.arange(8), 8)
        self.assertEqual(static_index, 4)

    def test_full_file_with_larger_grid_than_config_raises_without_nmat(self):
        # a legacy full-grid file bigger than the config Nmat cannot be
        # proven full (could be a [0,K] restriction of a still larger run):
        # fail with instructions; with nmat metadata it resolves cleanly
        # (see test_restricted_file_with_nmat_metadata_ignores_config)
        with self.assertRaises(ValueError):
            self._write_and_load(1024, np.arange(2048), 2048)

    def test_full_file_with_larger_grid_and_nmat_metadata_works(self):
        chi0q, static_index = self._write_and_load(
            1024, np.arange(2048), 2048, file_nmat=2048)
        self.assertEqual(static_index, 1024)

    def test_legacy_flex_4d_restricted_freq_index_delegates(self):
        # legacy FLEX 4D file: FULL 16-point grid on axis 0 but a restricted
        # 2-entry freq_index that happens to match the LAST axis length
        # (norb=2); the frequency axis of the 4D raw layout is ALWAYS axis 0,
        # so this is a length mismatch -> delegate (None), never a binding
        # of freq_index to the orbital axis
        chi0q, static_index = self._write_and_load(
            1024, np.array([7, 8]), 16, shape=(16, 4, 2, 2))
        self.assertIsNone(static_index)

    def test_small_legacy_file_with_default_config_raises_with_guidance(self):
        # a legacy 8-point file (no nmat key) under the default config
        # Nmat=1024 cannot be disambiguated from a [0,7] restriction that
        # lacks the zero frequency; the error tells the user to set
        # Nmat=8 if the file holds a full grid
        with self.assertRaises(ValueError):
            self._write_and_load(1024, np.arange(8), 8)

    def test_restricted_file_with_nmat_metadata_ignores_config(self):
        chi0q, static_index = self._write_and_load(
            1024, np.arange(2, 7), 5, file_nmat=8)
        self.assertEqual(static_index, 2)

    def test_file_without_zero_frequency_raises(self):
        with self.assertRaises(ValueError):
            self._write_and_load(8, np.arange(0, 4), 4, file_nmat=8)

    def test_ref_format_6d_freq_axis_detected(self):
        # reference-format file (norb, norb, Nx, Ny, Nz, nmat): the frequency
        # axis is LAST; freq_index length must resolve the right axis
        chi0q, static_index = self._write_and_load(
            8, np.arange(8), 8, shape=(2, 2, 2, 1, 1, 8))
        self.assertEqual(static_index, 4)

    def test_ref_format_6d_without_metadata_delegates(self):
        # a metadata-less reference-format file must NOT get a static index
        # guessed from the wrong axis (shape[0] = norb); returning None lets
        # the caller slice the center of the axis it actually uses
        chi0q, static_index = self._write_and_load(
            1024, None, 8, shape=(2, 2, 2, 1, 1, 8))
        self.assertIsNone(static_index)

    def test_legacy_file_without_metadata_delegates(self):
        chi0q, static_index = self._write_and_load(1024, None, 8)
        self.assertIsNone(static_index)


class TestLoadFlexSusceptibilitiesStaticSlice(unittest.TestCase):
    def test_restricted_rpa_chiq_uses_metadata_for_static_slice(self):
        # RPA chiq.npz with matsubara_frequency=[400,600] of Nmat=1024:
        # the zero bosonic frequency (original index 512) sits at stored
        # position 112, NOT at the axis center (100)
        from hwave.sc import _load_flex_susceptibilities

        with tempfile.TemporaryDirectory() as tmp:
            # norb=1 "kuroki" (reduced/spin-orbital) data is physically
            # shaped nd_so = norb*2 = 2, not the orbital-pair nd = norb^2 = 1
            # -- _expand_flex_chi now enforces that the stored shape agrees
            # with the chi_convention tag.
            nfreq, nvol, nd = 201, 2, 2
            freq_index = np.arange(400, 601)
            # encode the stored position in the value so the selected
            # slice is observable
            chiq = np.zeros((nfreq, nvol, nd, nd), dtype=np.complex128)
            for i in range(nfreq):
                chiq[i] = float(i)
            for name in ("chiq_s.npz", "chiq_c.npz"):
                np.savez(os.path.join(tmp, name),
                         chiq=chiq, freq_index=freq_index, nmat=1024,
                         chi_convention="kuroki")
            input_dict = {
                "mode": {"param": {"Nmat": 1024}},
                "file": {"input": {"path_to_flex_output": tmp},
                         "output": {"path_to_output": tmp}},
            }
            chis, chic, green_dressed, conv = _load_flex_susceptibilities(
                input_dict, norb=1, Nx=2, Ny=1, Nz=1)
            self.assertEqual(chis.flat[0].real, 112.0,
                             "static slice must come from freq_index/nmat "
                             "metadata, not the axis center")
            self.assertEqual(chic.flat[0].real, 112.0)

    def test_legacy_chiq_without_metadata_uses_data_axis_center(self):
        from hwave.sc import _load_flex_susceptibilities

        with tempfile.TemporaryDirectory() as tmp:
            # norb=1 "kuroki" (reduced/spin-orbital) data is physically
            # shaped nd_so = norb*2 = 2, not the orbital-pair nd = norb^2 = 1
            # -- _expand_flex_chi now enforces that the stored shape agrees
            # with the chi_convention tag.
            nfreq, nvol, nd = 8, 2, 2
            chiq = np.zeros((nfreq, nvol, nd, nd), dtype=np.complex128)
            for i in range(nfreq):
                chiq[i] = float(i)
            for name in ("chiq_s.npz", "chiq_c.npz"):
                np.savez(os.path.join(tmp, name),
                         chiq=chiq, chi_convention="kuroki")
            input_dict = {
                "mode": {"param": {"Nmat": 1024}},
                "file": {"input": {"path_to_flex_output": tmp},
                         "output": {"path_to_output": tmp}},
            }
            chis, chic, green_dressed, conv = _load_flex_susceptibilities(
                input_dict, norb=1, Nx=2, Ny=1, Nz=1)
            self.assertEqual(chis.flat[0].real, 4.0,
                             "metadata-less file: center of the data axis")


class TestComputeVerticesStaticIndex(unittest.TestCase):
    def test_simple_vertices_use_static_index(self):
        norb, Nx, Ny, Nz, nmat = 1, 2, 1, 1, 8
        rng = np.random.default_rng(3)
        chi0q = 0.1 * (rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat))
                       + 1j * rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat)))
        inter_k = {"CoulombIntra": np.full((norb, norb, Nx, Ny, Nz), 2.0,
                                           dtype=complex)}

        ref = _compute_vertices_simple(chi0q, inter_k, norb, Nx, Ny, Nz, nmat)

        # restricted file: keep frequency slices 2..6 (zero frequency at pos 2)
        chi0q_res = chi0q[..., 2:7]
        got = _compute_vertices_simple(chi0q_res, inter_k, norb, Nx, Ny, Nz,
                                       chi0q_res.shape[-1], static_index=2)

        np.testing.assert_allclose(got[0], ref[0], rtol=0.0, atol=1e-12)
        np.testing.assert_allclose(got[1], ref[1], rtol=0.0, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
