"""Tests for the UHFk-to-mVMC density checks.

See ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
"""
from __future__ import annotations

import sys, os, tempfile
from pathlib import Path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import numpy as np
import pytest

from tools._uhfk_to_mvmc.density_check import (
    density_from_amplitudes,
    compare_against_onebodyg_uhf,
    compare_against_onebodyg_uhf_general,
    DensityMismatchError,
    parse_uhf_cisajs_dat,
)


_DATA_DIR = Path(__file__).parent / "data"

# Tracked snapshots of H-wave outputs for case_soc_rashba_2d_sub, so gauge
# tests do not depend on the gitignored tests/validation/**/work/ tree.
_V35_GREEN_PATH = str(_DATA_DIR / "v35_case_soc_rashba_2d_sub_green.npz")
_V35_EIGEN_PATH = str(_DATA_DIR / "v35_case_soc_rashba_2d_sub_eigen.npz")
_V35_OCC_PATH = str(_DATA_DIR / "v35_case_soc_rashba_2d_sub_occupation.npz")

# case_soc_rashba_2d_sub_apbc snapshots. Same geometry as
# case_soc_rashba_2d_sub but with BoundaryCondition = AP-P-P (theta =
# (pi, 0, 0)); Ncond = 8 to match the physical filling.
_V36_APBC_GREEN_PATH = str(
    _DATA_DIR / "v36_case_soc_rashba_2d_sub_apbc_green.npz"
)
_V36_APBC_EIGEN_PATH = str(
    _DATA_DIR / "v36_case_soc_rashba_2d_sub_apbc_eigen.npz"
)
_V36_APBC_OCC_PATH = str(
    _DATA_DIR / "v36_case_soc_rashba_2d_sub_apbc_occupation.npz"
)
_V36_APBC_COMPOSITE_PATH = str(
    _DATA_DIR / "v36_case_soc_rashba_2d_sub_apbc_composite_element.json"
)


def test_density_from_amplitudes_pbc_k0_only():
    """Single occupied k=0 plane wave on L=4 → G^sigma_{ij} = 1/L."""
    L = 4
    A = np.full((L, 1), 1.0 / np.sqrt(L), dtype=np.complex128)
    G = density_from_amplitudes(A)
    expected = np.full((L, L), 1.0 / L, dtype=np.complex128)
    np.testing.assert_allclose(G, expected, atol=1e-12)


def test_density_check_soc_convention_asymmetric():
    """Under ``is_soc_mode=True``, density_check reads greenone.dat rows
    on mVMC spin-block order (``all_i = i + s * Ns``) and passes when the
    bridge output matches. If the bridge used interleaved order, the
    asymmetric golden values would force a numerical mismatch above
    tolerance."""
    Ns = 2
    G_all = np.zeros((2 * Ns, 2 * Ns), dtype=complex)
    G_all[0, 3] = 0.3
    G_all[1, 2] = 0.2 + 0.1j
    G_all[2, 1] = 0.2 - 0.1j
    G_all[3, 0] = 0.3

    golden = str(_DATA_DIR / "soc_greenone_asymmetric_golden.dat")

    # Should pass without raising.
    compare_against_onebodyg_uhf_general(
        G_all, golden, tol=1e-10, is_soc_mode=True,
    )


def test_density_check_soc_rejects_beyond_tol():
    """SOC mode with an s != t entry mismatched beyond tolerance -> raise
    DensityMismatchError."""
    Ns = 2
    G_all = np.zeros((2 * Ns, 2 * Ns), dtype=complex)
    G_all[0, 3] = 0.3
    G_all[1, 2] = 0.2 + 0.1j + 3e-10  # deliberate error above 1e-10 tol
    G_all[2, 1] = 0.2 - 0.1j
    G_all[3, 0] = 0.3

    golden = str(_DATA_DIR / "soc_greenone_asymmetric_golden.dat")

    with pytest.raises(DensityMismatchError, match="differ"):
        compare_against_onebodyg_uhf_general(
            G_all, golden, tol=1e-10, is_soc_mode=True,
        )


def test_density_check_v3_still_rejects_s_neq_t_when_not_soc_mode():
    """Under ``is_soc_mode=False``, if the reference
    greenone.dat has a non-zero s != t entry while the bridge's G_all is
    zero in that off-diagonal block (which is the natural Sz-diagonal
    output), the comparison must still flag the mismatch.

    The bridge builds a spin-block-diagonal
    G_all, so any non-zero s != t entry in the reference will differ from
    the bridge's zero and trigger DensityMismatchError.
    """
    Ns = 2
    G_all = np.zeros((2 * Ns, 2 * Ns), dtype=complex)
    G_all[0, 0] = 0.5
    G_all[1, 1] = 0.5
    G_all[2, 2] = 0.5
    G_all[3, 3] = 0.5  # spin-block-diagonal only, off-diagonal blocks all 0

    # Reference has a non-zero s != t entry, incompatible with this scope.
    with tempfile.NamedTemporaryFile(
        "w", suffix=".dat", delete=False
    ) as tmp:
        tmp.write("0 0 0 0  0.5 0.0\n")
        tmp.write("1 0 1 0  0.5 0.0\n")
        tmp.write("0 1 0 1  0.5 0.0\n")
        tmp.write("1 1 1 1  0.5 0.0\n")
        tmp.write("0 0 1 1  0.3 0.0\n")  # off-diagonal block: not zero
        path = tmp.name

    try:
        with pytest.raises(DensityMismatchError, match="differ"):
            compare_against_onebodyg_uhf_general(
                G_all, path, tol=1e-10, is_soc_mode=False,
            )
    finally:
        os.unlink(path)


def test_gauge_lift_fft_round_trip():
    """fftn(ifftn(...)) must recover
    the original green_sublattice at 1e-14 element-wise, matching H-wave's
    forward-normalized FFT convention (_save_green at uhfk.py:2500-2510).
    """
    import numpy as np
    # Load the tracked H-wave green snapshot for case_soc_rashba_2d_sub.
    gp = np.load(_V35_GREEN_PATH)
    green_sublattice = gp["green_sublattice"]
    assert green_sublattice.ndim == 5, (
        f"expected raw H-wave 5D green shape (nvol, ns, norb, ns, norb); "
        f"got {green_sublattice.shape}"
    )
    # Extract the SOC block: under enable_spin_orbital=True, ns=1 and
    # norb=2*norb_phys; the physical density lives in
    # green_sub[:, 0, :, 0, :].
    gs_soc = green_sublattice[:, 0, :, 0, :]
    # Reshape onto the 3 spatial axes for 3D FFT.
    # case_soc_rashba_2d_sub: CellShape=[6,4,1], SubShape=[2,2,1] -> L_folded=[3,2,1]
    L_folded = np.array([3, 2, 1], dtype=np.int64)
    gs_reshaped = gs_soc.reshape(
        L_folded[0], L_folded[1], L_folded[2],
        gs_soc.shape[1], gs_soc.shape[2],
    )
    G_k = np.fft.ifftn(gs_reshaped, axes=(0, 1, 2), norm="forward")
    gs_round = np.fft.fftn(G_k, axes=(0, 1, 2), norm="forward")
    assert np.allclose(gs_round, gs_reshaped, atol=1e-14), (
        "FFT round-trip residual exceeds 1e-14: "
        f"max abs diff = {np.max(np.abs(gs_round - gs_reshaped)):.3e}"
    )


def test_gauge_lift_sub_offset_sign_audit():
    """Swapping the sign of sub_offset
    in gauge_lift's phase must make the reconstructed G disagree with the
    correct shipping A density by a Rashba-scale amount. See
    ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    import numpy as np
    from tools._uhfk_to_mvmc.density_check import gauge_lift

    # Load tracked case_soc_rashba_2d_sub H-wave snapshot.
    gp = np.load(_V35_GREEN_PATH)
    green_sublattice = gp["green_sublattice"]

    # Fixed inputs derived from input.toml.
    cell_shape = np.array([6, 4, 1], dtype=np.int64)
    subshape = np.array([2, 2, 1], dtype=np.int64)
    Ns = int(np.prod(cell_shape))  # 24
    site_positions = np.array(
        [[x, y, z] for z in range(1) for y in range(4) for x in range(6)],
        dtype=np.int64,
    )
    assert site_positions.shape == (Ns, 3)
    folded_cell_of = lambda r_phys: r_phys // subshape

    # Sanity: gauge_lift on a site-diagonal, spin-diagonal element must be
    # a real positive occupation number close to Ncond/2Ns = 8/48 = 1/6.
    diag = gauge_lift(
        green_sublattice, site_i=0, spin_i=0, site_j=0, spin_j=0,
        subshape=subshape, cell_shape=cell_shape,
        site_positions=site_positions, folded_cell_of=folded_cell_of,
        boundary_theta=(0.0, 0.0, 0.0),
    )
    assert abs(diag.imag) < 1e-10, f"diagonal G[0,0] must be real; got {diag}"
    assert 0.05 < diag.real < 0.3, (
        f"expected diagonal occupation ~1/6; got {diag.real}"
    )

    # Mutation: gauge_lift with sub_offset sign flipped. The mutation is
    # only visible on off-diagonal-in-sub_offset pairs. Pick two sites
    # with distinct sub_offset.
    # site 0 -> r=(0,0,0), fc=(0,0,0), so=(0,0,0)
    # site 1 -> r=(1,0,0), fc=(0,0,0), so=(1,0,0)  <-- distinct sub_offset
    def gauge_lift_flipped(site_i, spin_i, site_j, spin_j):
        # Same as gauge_lift but so_i and so_j are negated in the phase.
        gs_soc = green_sublattice[:, 0, :, 0, :]
        L_folded = cell_shape // subshape
        gs = gs_soc.reshape(
            L_folded[0], L_folded[1], L_folded[2],
            gs_soc.shape[1], gs_soc.shape[2],
        )
        G_k = np.fft.ifftn(gs, axes=(0, 1, 2), norm="forward")
        r_i = site_positions[site_i]
        r_j = site_positions[site_j]
        fc_i = r_i // subshape
        fc_j = r_j // subshape
        so_i = r_i - fc_i * subshape
        so_j = r_j - fc_j * subshape
        folded_orb_i = int(so_i[0] + subshape[0] * (so_i[1] + subshape[1] * so_i[2]))
        folded_orb_j = int(so_j[0] + subshape[0] * (so_j[1] + subshape[1] * so_j[2]))
        aa = 2 * folded_orb_i + spin_i
        bb = 2 * folded_orb_j + spin_j
        accum = 0.0j
        for kx in range(L_folded[0]):
            for ky in range(L_folded[1]):
                for kz in range(L_folded[2]):
                    k_vec = 2.0 * np.pi * np.array([kx, ky, kz]) / L_folded
                    # FLIPPED sub_offset sign:
                    r_diff = (fc_j.astype(np.float64) - so_j.astype(np.float64)) \
                             - (fc_i.astype(np.float64) - so_i.astype(np.float64))
                    phase = np.exp(-1j * np.dot(k_vec, r_diff))
                    accum += G_k[kx, ky, kz, aa, bb] * phase
        return accum / float(np.prod(L_folded))

    # Compare on site (0, up) vs site (1, up) which have different sub_offset.
    correct = gauge_lift(
        green_sublattice, site_i=0, spin_i=0, site_j=1, spin_j=0,
        subshape=subshape, cell_shape=cell_shape,
        site_positions=site_positions, folded_cell_of=folded_cell_of,
        boundary_theta=(0.0, 0.0, 0.0),
    )
    mutated = gauge_lift_flipped(0, 0, 1, 0)
    delta = abs(correct - mutated)
    assert delta > 1e-3, (
        f"sub_offset sign flip should perturb the G element by >1e-3 "
        f"(Rashba scale). Got delta = {delta:.3e}, correct = {correct}, "
        f"mutated = {mutated}"
    )


def test_gauge_lift_absolute_scale_and_full_2Ns_2Ns_coverage():
    """Build the full 2Ns × 2Ns
    G_ref matrix from gauge_lift and compare to G_direct =
    np.conj(A_ship) @ A_ship.T over every element, including s != t
    cross-spin blocks. Asserts 5e-13 element-wise agreement (FP roundoff
    floor for the 24-site 3D FFT + 8-column A @ A.T sum; still 200x
    tighter than the 1e-10 ship gate in
    compare_against_green_sublattice), 1% absolute scale, and > 0.01
    cross-spin magnitude (proves the test actually exercises Rashba SOC
    off-diagonals). See ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    import numpy as np
    from tools._uhfk_to_mvmc.density_check import gauge_lift
    from tools._uhfk_to_mvmc.general_fij_builder import (
        build_slater_orbitals, build_pair_list, compute_canonical_reps,
    )
    from tools._uhfk_to_mvmc.occupation_step import step_occupation
    from tools._uhfk_to_mvmc.partner_index import find_partner_rows

    gp = np.load(_V35_GREEN_PATH)
    eigen = np.load(_V35_EIGEN_PATH)
    occ = np.load(_V35_OCC_PATH)

    green_sublattice = gp["green_sublattice"]
    eigenvector = eigen["eigenvector"]
    eigenvalue = eigen["eigenvalue"]
    wavevector_index = eigen["wavevector_index"]

    cell_shape = np.array([6, 4, 1], dtype=np.int64)
    subshape = np.array([2, 2, 1], dtype=np.int64)
    Ns = int(np.prod(cell_shape))  # 24
    Ncond = 8
    site_positions = np.array(
        [[x, y, z] for z in range(1) for y in range(4) for x in range(6)],
        dtype=np.int64,
    )
    folded_cell_of = lambda r_phys: r_phys // subshape

    # Step 1: build the shipping A via build_slater_orbitals SOC branch.
    stepped, _ = step_occupation(
        occ["occupation"], eigenvalue, occ["column_spin"],
        occ["column_mu_group"], float(occ["T"]),
        ncond_per_group=[Ncond],
        is_soc_mode=True,
    )
    partner_rows, _ = find_partner_rows(
        wavevector_index, np.zeros(3, dtype=np.float64),
        cell_shape // subshape,
    )
    canonical, _ = compute_canonical_reps(partner_rows, wavevector_index)
    pair_list = build_pair_list(
        stepped, occ["column_spin"], canonical, partner_rows, is_soc_mode=True,
    )
    A_ship = build_slater_orbitals(
        wavevector_index=wavevector_index,
        eigenvector=eigenvector,
        column_spin=occ["column_spin"],
        site_positions=site_positions.astype(np.int64),
        cell_shape=cell_shape,
        subshape=subshape,
        theta=np.zeros(3, dtype=np.float64),
        pair_list=pair_list,
        is_soc_mode=True,
    )
    G_direct = np.conj(A_ship) @ A_ship.T  # (2Ns, 2Ns)

    # Step 2: build G_ref via gauge_lift for all (site_i, spin_i, site_j, spin_j).
    G_ref = np.zeros((2 * Ns, 2 * Ns), dtype=complex)
    for site_i in range(Ns):
        for spin_i in (0, 1):
            for site_j in range(Ns):
                for spin_j in (0, 1):
                    all_i = site_i + spin_i * Ns
                    all_j = site_j + spin_j * Ns
                    G_ref[all_i, all_j] = gauge_lift(
                        green_sublattice, site_i, spin_i, site_j, spin_j,
                        subshape=subshape, cell_shape=cell_shape,
                        site_positions=site_positions,
                        folded_cell_of=folded_cell_of,
                        boundary_theta=(0.0, 0.0, 0.0),
                    )

    max_diff = np.max(np.abs(G_ref - G_direct))
    assert max_diff < 5e-13, (
        f"|G_ref - G_direct|_max = {max_diff:.3e} > 5e-13; the gauge_lift "
        "does not match the shipping-A density."
    )

    scale_ref = np.max(np.abs(G_ref))
    scale_direct = np.max(np.abs(G_direct))
    assert abs(scale_ref - scale_direct) / max(scale_direct, 1e-30) < 0.01, (
        f"|G_ref|_max = {scale_ref}, |G_direct|_max = {scale_direct}; "
        "absolute-scale mismatch > 1% suggests a 1/N_folded scale bug."
    )

    # Cross-spin coverage: pick site 0 (up) x site 0 (dn) etc., cross-spin
    # block G[0..Ns, Ns..2Ns].
    cross_spin_block = G_ref[:Ns, Ns:]
    cross_max = np.max(np.abs(cross_spin_block))
    assert cross_max > 0.01, (
        f"|G_ref[s!=t]|_max = {cross_max:.3e} <= 0.01; the fixture does not "
        "exercise Rashba SOC off-diagonals as expected. Verify fixture "
        "parameters (Ncond, alpha_R) or the pack aa/bb formula."
    )


def test_gauge_lift_catches_orientation_swap_mutation():
    """Build G with the WRONG
    orientation A_ship @ A_ship.conj().T (transposed/conjugated) and
    verify gauge_lift disagrees by > 1e-3 element-wise. See
    ``docs/en/source/algorithm/uhfk_to_mvmc.rst`` for the density orientation.
    """
    import numpy as np
    from tools._uhfk_to_mvmc.density_check import gauge_lift
    from tools._uhfk_to_mvmc.general_fij_builder import (
        build_slater_orbitals, build_pair_list, compute_canonical_reps,
    )
    from tools._uhfk_to_mvmc.occupation_step import step_occupation
    from tools._uhfk_to_mvmc.partner_index import find_partner_rows

    gp = np.load(_V35_GREEN_PATH)
    eigen = np.load(_V35_EIGEN_PATH)
    occ = np.load(_V35_OCC_PATH)

    green_sublattice = gp["green_sublattice"]

    cell_shape = np.array([6, 4, 1], dtype=np.int64)
    subshape = np.array([2, 2, 1], dtype=np.int64)
    Ns = int(np.prod(cell_shape))
    Ncond = 8
    site_positions = np.array(
        [[x, y, z] for z in range(1) for y in range(4) for x in range(6)],
        dtype=np.int64,
    )
    folded_cell_of = lambda r_phys: r_phys // subshape

    stepped, _ = step_occupation(
        occ["occupation"], eigen["eigenvalue"], occ["column_spin"],
        occ["column_mu_group"], float(occ["T"]),
        ncond_per_group=[Ncond],
        is_soc_mode=True,
    )
    partner_rows, _ = find_partner_rows(
        eigen["wavevector_index"], np.zeros(3, dtype=np.float64),
        cell_shape // subshape,
    )
    canonical, _ = compute_canonical_reps(partner_rows, eigen["wavevector_index"])
    pair_list = build_pair_list(
        stepped, occ["column_spin"], canonical, partner_rows, is_soc_mode=True,
    )
    A_ship = build_slater_orbitals(
        wavevector_index=eigen["wavevector_index"],
        eigenvector=eigen["eigenvector"],
        column_spin=occ["column_spin"],
        site_positions=site_positions.astype(np.int64),
        cell_shape=cell_shape,
        subshape=subshape,
        theta=np.zeros(3, dtype=np.float64),
        pair_list=pair_list,
        is_soc_mode=True,
    )
    G_wrong = A_ship @ np.conj(A_ship).T  # WRONG orientation

    # Compare G_wrong to gauge_lift on a couple of cross-spin elements
    # where the orientation matters.
    max_disagree = 0.0
    for site_i in range(0, Ns, 4):  # sample every 4th site to keep it fast
        for site_j in range(0, Ns, 4):
            for spin_i, spin_j in [(0, 1), (1, 0)]:
                all_i = site_i + spin_i * Ns
                all_j = site_j + spin_j * Ns
                gref = gauge_lift(
                    green_sublattice, site_i, spin_i, site_j, spin_j,
                    subshape=subshape, cell_shape=cell_shape,
                    site_positions=site_positions,
                    folded_cell_of=folded_cell_of,
                    boundary_theta=(0.0, 0.0, 0.0),
                )
                disagree = abs(G_wrong[all_i, all_j] - gref)
                if disagree > max_disagree:
                    max_disagree = disagree
    assert max_disagree > 1e-3, (
        f"orientation-swap mutation must disagree with gauge_lift by "
        f">1e-3 on some cross-spin element; got max disagreement "
        f"{max_disagree:.3e}. If this passes trivially, either A_ship "
        "or gauge_lift is real-only for the fixture (Rashba should give "
        "complex off-diagonals)."
    )


def test_compare_against_green_sublattice_soc_mode_passes_on_correct_A_ship():
    """compare_against_green_sublattice with
    is_soc_sublattice_mode=True must PASS at 1e-10 tolerance when the
    input G is the correct shipping-A density.
    """
    import numpy as np
    from tools._uhfk_to_mvmc.density_check import (
        compare_against_green_sublattice,
    )
    from tools._uhfk_to_mvmc.general_fij_builder import (
        build_slater_orbitals, build_pair_list, compute_canonical_reps,
    )
    from tools._uhfk_to_mvmc.occupation_step import step_occupation
    from tools._uhfk_to_mvmc.partner_index import find_partner_rows

    eigen = np.load(_V35_EIGEN_PATH)
    occ = np.load(_V35_OCC_PATH)

    cell_shape = np.array([6, 4, 1], dtype=np.int64)
    subshape = np.array([2, 2, 1], dtype=np.int64)
    Ns = int(np.prod(cell_shape))
    Ncond = 8
    site_positions = np.array(
        [[x, y, z] for z in range(1) for y in range(4) for x in range(6)],
        dtype=np.int64,
    )

    stepped, _ = step_occupation(
        occ["occupation"], eigen["eigenvalue"], occ["column_spin"],
        occ["column_mu_group"], float(occ["T"]),
        ncond_per_group=[Ncond],
        is_soc_mode=True,
    )
    partner_rows, _ = find_partner_rows(
        eigen["wavevector_index"], np.zeros(3, dtype=np.float64),
        cell_shape // subshape,
    )
    canonical, _ = compute_canonical_reps(partner_rows, eigen["wavevector_index"])
    pair_list = build_pair_list(
        stepped, occ["column_spin"], canonical, partner_rows, is_soc_mode=True,
    )
    A_ship = build_slater_orbitals(
        wavevector_index=eigen["wavevector_index"],
        eigenvector=eigen["eigenvector"],
        column_spin=occ["column_spin"],
        site_positions=site_positions.astype(np.int64),
        cell_shape=cell_shape,
        subshape=subshape,
        theta=np.zeros(3, dtype=np.float64),
        pair_list=pair_list,
        is_soc_mode=True,
    )
    G_ship = np.conj(A_ship) @ A_ship.T

    # Must pass at 1e-10 without raising.
    compare_against_green_sublattice(
        G_ship,
        green_path=_V35_GREEN_PATH,
        Ns=Ns,
        tol=1e-10,
        is_soc_sublattice_mode=True,
        cell_shape=cell_shape,
        subshape=subshape,
        site_positions=site_positions,
    )


def test_compare_against_green_sublattice_soc_mode_fires_on_orientation_swap():
    """compare_against_green_sublattice with
    is_soc_sublattice_mode=True must RAISE DensityMismatchError when
    the input G is transposed/conjugated (orientation swap mutation).
    """
    import numpy as np
    from tools._uhfk_to_mvmc.density_check import (
        compare_against_green_sublattice, DensityMismatchError,
    )
    from tools._uhfk_to_mvmc.general_fij_builder import (
        build_slater_orbitals, build_pair_list, compute_canonical_reps,
    )
    from tools._uhfk_to_mvmc.occupation_step import step_occupation
    from tools._uhfk_to_mvmc.partner_index import find_partner_rows

    eigen = np.load(_V35_EIGEN_PATH)
    occ = np.load(_V35_OCC_PATH)

    cell_shape = np.array([6, 4, 1], dtype=np.int64)
    subshape = np.array([2, 2, 1], dtype=np.int64)
    Ns = int(np.prod(cell_shape))
    Ncond = 8
    site_positions = np.array(
        [[x, y, z] for z in range(1) for y in range(4) for x in range(6)],
        dtype=np.int64,
    )

    stepped, _ = step_occupation(
        occ["occupation"], eigen["eigenvalue"], occ["column_spin"],
        occ["column_mu_group"], float(occ["T"]),
        ncond_per_group=[Ncond],
        is_soc_mode=True,
    )
    partner_rows, _ = find_partner_rows(
        eigen["wavevector_index"], np.zeros(3, dtype=np.float64),
        cell_shape // subshape,
    )
    canonical, _ = compute_canonical_reps(partner_rows, eigen["wavevector_index"])
    pair_list = build_pair_list(
        stepped, occ["column_spin"], canonical, partner_rows, is_soc_mode=True,
    )
    A_ship = build_slater_orbitals(
        wavevector_index=eigen["wavevector_index"],
        eigenvector=eigen["eigenvector"],
        column_spin=occ["column_spin"],
        site_positions=site_positions.astype(np.int64),
        cell_shape=cell_shape,
        subshape=subshape,
        theta=np.zeros(3, dtype=np.float64),
        pair_list=pair_list,
        is_soc_mode=True,
    )
    G_wrong = A_ship @ np.conj(A_ship).T  # WRONG orientation.
    with pytest.raises(DensityMismatchError):
        compare_against_green_sublattice(
            G_wrong,
            green_path=_V35_GREEN_PATH,
            Ns=Ns,
            tol=1e-10,
            is_soc_sublattice_mode=True,
            cell_shape=cell_shape,
            subshape=subshape,
            site_positions=site_positions,
        )


# ---------------------------------------------------------------------
# APBC pins on the case_soc_rashba_2d_sub_apbc snapshot.
# ---------------------------------------------------------------------


def _v36_apbc_build_shipping_A():
    """Load the APBC snapshots and build the shipping A matrix under
    ``boundary_theta = (pi, 0, 0)`` (AP-P-P). Returns (A_ship, geometry
    inputs) usable by both the positive pin and the negative control.
    """
    import numpy as np
    from tools._uhfk_to_mvmc.general_fij_builder import (
        build_slater_orbitals, build_pair_list, compute_canonical_reps,
    )
    from tools._uhfk_to_mvmc.occupation_step import step_occupation
    from tools._uhfk_to_mvmc.partner_index import find_partner_rows

    eigen = np.load(_V36_APBC_EIGEN_PATH)
    occ = np.load(_V36_APBC_OCC_PATH)

    cell_shape = np.array([6, 4, 1], dtype=np.int64)
    subshape = np.array([2, 2, 1], dtype=np.int64)
    Ns = int(np.prod(cell_shape))
    Ncond = 8
    site_positions = np.array(
        [[x, y, z] for z in range(1) for y in range(4) for x in range(6)],
        dtype=np.int64,
    )
    theta_apbc = np.array([np.pi, 0.0, 0.0], dtype=np.float64)

    stepped, _ = step_occupation(
        occ["occupation"], eigen["eigenvalue"], occ["column_spin"],
        occ["column_mu_group"], float(occ["T"]),
        ncond_per_group=[Ncond],
        is_soc_mode=True,
    )
    partner_rows, _ = find_partner_rows(
        eigen["wavevector_index"], theta_apbc, cell_shape // subshape,
    )
    canonical, _ = compute_canonical_reps(partner_rows, eigen["wavevector_index"])
    pair_list = build_pair_list(
        stepped, occ["column_spin"], canonical, partner_rows, is_soc_mode=True,
    )
    A_ship = build_slater_orbitals(
        wavevector_index=eigen["wavevector_index"],
        eigenvector=eigen["eigenvector"],
        column_spin=occ["column_spin"],
        site_positions=site_positions.astype(np.int64),
        cell_shape=cell_shape,
        subshape=subshape,
        theta=theta_apbc,
        pair_list=pair_list,
        is_soc_mode=True,
    )
    return A_ship, {
        "cell_shape": cell_shape, "subshape": subshape, "Ns": Ns,
        "site_positions": site_positions, "theta_apbc": theta_apbc,
    }


def test_gauge_lift_apbc_matches_shipping_A_on_real_run():
    """On the case_soc_rashba_2d_sub_apbc snapshot with BoundaryCondition = AP-P-P
    (theta = (pi, 0, 0)), the shipping A density ``conj(A) @ A.T`` must
    equal ``gauge_lift(green_sublattice, ..., boundary_theta=(pi, 0, 0))``
    over every (all_i, all_j) at max_abs_delta < 1e-10.

    Under APBC this proves the two-step gauge composition
    ``exp(-i k_folded . dr_folded) * exp(-i theta . dr_full / L_full)``
    matches the shipping A's own APBC composition. See
    ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    import numpy as np
    from tools._uhfk_to_mvmc.density_check import gauge_lift

    gp = np.load(_V36_APBC_GREEN_PATH)
    green_sublattice = gp["green_sublattice"]

    A_ship, geom = _v36_apbc_build_shipping_A()
    Ns = geom["Ns"]
    G_direct = np.conj(A_ship) @ A_ship.T

    folded_cell_of = lambda r_phys: r_phys // geom["subshape"]

    G_ref = np.zeros((2 * Ns, 2 * Ns), dtype=complex)
    for site_i in range(Ns):
        for spin_i in (0, 1):
            for site_j in range(Ns):
                for spin_j in (0, 1):
                    all_i = site_i + spin_i * Ns
                    all_j = site_j + spin_j * Ns
                    G_ref[all_i, all_j] = gauge_lift(
                        green_sublattice, site_i, spin_i, site_j, spin_j,
                        subshape=geom["subshape"],
                        cell_shape=geom["cell_shape"],
                        site_positions=geom["site_positions"],
                        folded_cell_of=folded_cell_of,
                        boundary_theta=geom["theta_apbc"],
                    )

    max_diff = np.max(np.abs(G_ref - G_direct))
    assert max_diff < 1e-10, (
        f"|G_ref - G_direct|_max = {max_diff:.3e} > 1e-10 under APBC; the "
        "gauge_lift APBC composition does not match the shipping-A APBC "
        "composition on case_soc_rashba_2d_sub_apbc."
    )

    # Cross-spin coverage sanity: the fixture must actually exercise
    # Rashba SOC off-diagonals (otherwise the pin degenerates to a PBC
    # check trivially).
    cross_spin_block = G_ref[:Ns, Ns:]
    cross_max = np.max(np.abs(cross_spin_block))
    assert cross_max > 0.01, (
        f"|G_ref[s!=t]|_max = {cross_max:.3e} <= 0.01; the APBC fixture "
        "does not exercise Rashba SOC off-diagonals."
    )


def test_gauge_lift_apbc_negative_control_theta_zero():
    """At the fixture's selected composite element (i_c, s_c, j_c, t_c),
    compute:

    - ``G_A_c``  = shipping A density at the composite under theta = (pi, 0, 0)
    - ``G_neg_c`` = ``gauge_lift(green_sublattice, ..., boundary_theta = (0, 0, 0))``
      at the same composite (twist gauge intentionally dropped)

    Assert ``|G_neg_c - G_A_c| >= 0.3 * |G_A_c|``. This proves the twist
    gauge in ``gauge_lift`` is what makes the APBC composition work; if
    the twist term were vestigial or already carried elsewhere, dropping
    it would leave ``G_neg_c`` numerically close to ``G_A_c`` at machine
    precision. See ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    import json
    import numpy as np
    from tools._uhfk_to_mvmc.density_check import gauge_lift

    with open(_V36_APBC_COMPOSITE_PATH) as fp:
        manifest = json.load(fp)
    i_c, s_c = int(manifest["i_c"]), int(manifest["s_c"])
    j_c, t_c = int(manifest["j_c"]), int(manifest["t_c"])

    A_ship, geom = _v36_apbc_build_shipping_A()
    Ns = geom["Ns"]

    G_direct = np.conj(A_ship) @ A_ship.T
    all_i = i_c + s_c * Ns
    all_j = j_c + t_c * Ns
    G_A_c = G_direct[all_i, all_j]

    gp = np.load(_V36_APBC_GREEN_PATH)
    green_sublattice = gp["green_sublattice"]
    folded_cell_of = lambda r_phys: r_phys // geom["subshape"]

    G_neg_c = gauge_lift(
        green_sublattice, i_c, s_c, j_c, t_c,
        subshape=geom["subshape"], cell_shape=geom["cell_shape"],
        site_positions=geom["site_positions"],
        folded_cell_of=folded_cell_of,
        boundary_theta=(0.0, 0.0, 0.0),  # twist gauge intentionally dropped
    )

    G_A_c_abs = abs(G_A_c)
    assert G_A_c_abs > 1e-3, (
        f"G_A_c magnitude {G_A_c_abs:.3e} <= 1e-3 at composite "
        f"(i={i_c}, s={s_c}, j={j_c}, t={t_c}); the composite fell "
        "off the manifest; check the snapshot / manifest pairing."
    )
    delta = abs(G_neg_c - G_A_c)
    threshold = 0.3 * G_A_c_abs
    assert delta >= threshold, (
        f"|G_neg_c - G_A_c| = {delta:.3e} < 0.3 * |G_A_c| = "
        f"{threshold:.3e}; dropping the twist gauge did NOT perturb the "
        "composite element by the expected Rashba-scale amount. Either "
        "the twist gauge is vestigial in gauge_lift, or the composite "
        "was accidentally chosen with a sub_offset_x invariant to APBC."
    )


# ---------------------------------------------------------------------
# Adversarial mutation matrix. The validation contract is documented in
# ``docs/en/source/uhfk/tools/uhfk_to_mvmc.rst``.
# ---------------------------------------------------------------------


def _axis_of_mutator(mutator_id):
    """Return the axis index for a per-direction mutator id
    (M-gauge-1-x -> 0, M-gauge-1-y -> 1, M-gauge-1-z -> 2) or None for
    whole-vector mutators / baseline.

    Raises ValueError if the id looks like a per-direction mutator but
    the axis suffix is invalid (e.g., M-gauge-1-w).
    """
    if not mutator_id.startswith(("M-gauge-", "M-ship-")):
        return None  # baseline or unknown
    parts = mutator_id.split("-")
    # Whole-vector: "M-gauge-1" -> 3 parts
    # Per-direction: "M-gauge-1-x" -> 4 parts
    if len(parts) == 3:
        return None
    if len(parts) != 4:
        raise ValueError(f"invalid mutator id shape: {mutator_id!r}")
    axis_char = parts[3]
    if axis_char not in ("x", "y", "z"):
        raise ValueError(
            f"invalid axis in mutator id {mutator_id!r}: "
            f"got {axis_char!r}, expected x, y, or z"
        )
    return "xyz".index(axis_char)


def test_axis_of_mutator_v36_ids_return_none():
    assert _axis_of_mutator("baseline") is None
    assert _axis_of_mutator("M-gauge-1") is None
    assert _axis_of_mutator("M-ship-5") is None


def test_axis_of_mutator_v37_ids_return_axis_index():
    assert _axis_of_mutator("M-gauge-1-x") == 0
    assert _axis_of_mutator("M-gauge-1-y") == 1
    assert _axis_of_mutator("M-gauge-1-z") == 2
    assert _axis_of_mutator("M-ship-5-z") == 2


def test_axis_of_mutator_rejects_invalid_axis():
    with pytest.raises(ValueError):
        _axis_of_mutator("M-gauge-1-w")


def _v36_gauge_lift_at_composite_mutated(
    green_sublattice, i, s, j, t, subshape, cell_shape, site_positions,
    boundary_theta, mutator,
):
    """Shadow copy of tools/_uhfk_to_mvmc/density_check.py::gauge_lift at a
    single (i, s, j, t) element with an M-gauge mutation applied to the
    twist / dr_folded phase. Preserves every other step so any pass /
    fail is attributable to the specific mutation.

    Supported mutators:
      * ``baseline``   -- exact reproduction of the shipping gauge_lift
      * ``M-gauge-1``  -- twist sign: exp(-i theta.dr/L) -> exp(+i theta.dr/L)
      * ``M-gauge-2``  -- twist unit: pass theta / (2*pi) as theta
      * ``M-gauge-3``  -- twist L: divide by L_folded instead of L_full
      * ``M-gauge-4``  -- whole-vector: sub_offset sign in dr_folded;
                          per-direction: omit sub_offset on that axis
      * ``M-gauge-5``  -- dr_folded in twist phase instead of dr_full

    Per-direction ids (``M-gauge-{1..5}-{x,y,z}``) apply the mutation
    to only the named axis component
    of theta / L / dr, while the other two axes keep the baseline
    behavior. M-gauge-4 omits sub_offset on that axis because the
    sign flip is degenerate for ``L_folded=2``. See
    ``_axis_of_mutator``. See
    ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    gs_soc = green_sublattice[:, 0, :, 0, :]
    L_folded = cell_shape // subshape
    gs = gs_soc.reshape(
        L_folded[0], L_folded[1], L_folded[2],
        gs_soc.shape[1], gs_soc.shape[2],
    )
    G_k = np.fft.ifftn(gs, axes=(0, 1, 2), norm="forward")

    r_phys_i = site_positions[i].astype(np.int64)
    r_phys_j = site_positions[j].astype(np.int64)
    fc_i = r_phys_i // subshape
    fc_j = r_phys_j // subshape
    so_i = r_phys_i - fc_i * subshape
    so_j = r_phys_j - fc_j * subshape
    folded_orb_i = int(
        so_i[0] + subshape[0] * (so_i[1] + subshape[1] * so_i[2])
    )
    folded_orb_j = int(
        so_j[0] + subshape[0] * (so_j[1] + subshape[1] * so_j[2])
    )
    aa = 2 * folded_orb_i + int(s)
    bb = 2 * folded_orb_j + int(t)

    # Per-direction mutators (M-gauge-N-{x,y,z}) restrict the base
    # mutation N to a single axis; whole-vector ids (axis is None)
    # keep their original code paths below untouched.
    axis = _axis_of_mutator(mutator)
    base_mutator = (
        mutator if axis is None else "-".join(mutator.split("-")[:-1])
    )

    if base_mutator == "M-gauge-4":
        if axis is None:
            # sub_offset(j) - sub_offset(i) -> -(sub_offset(j) - sub_offset(i))
            dr_folded = (
                (fc_j.astype(np.float64) + (-so_j).astype(np.float64))
                - (fc_i.astype(np.float64) + (-so_i).astype(np.float64))
            )
        else:
            # Omit sub_offset on the named axis only;
            # other axes keep the baseline (+so_i, +so_j) contribution.
            dr_folded = (
                (fc_j.astype(np.float64) + so_j.astype(np.float64))
                - (fc_i.astype(np.float64) + so_i.astype(np.float64))
            )
            dr_folded[axis] = (
                fc_j.astype(np.float64)[axis]
                - fc_i.astype(np.float64)[axis]
            )
    else:
        dr_folded = (fc_j.astype(np.float64) + so_j.astype(np.float64)) \
                    - (fc_i.astype(np.float64) + so_i.astype(np.float64))
    dr_full = r_phys_j.astype(np.float64) - r_phys_i.astype(np.float64)
    L_full = (subshape * L_folded).astype(np.float64)

    theta = np.asarray(boundary_theta, dtype=np.float64)
    if base_mutator == "M-gauge-1":
        if axis is None:
            phase_twist = np.exp(+1j * np.dot(theta, dr_full / L_full))
        else:
            # flip sign of theta on the named axis only.
            theta_mut = theta.copy()
            theta_mut[axis] = -theta_mut[axis]
            phase_twist = np.exp(-1j * np.dot(theta_mut, dr_full / L_full))
    elif base_mutator == "M-gauge-2":
        if axis is None:
            # twist_offset = theta / (2*pi) passed as theta
            theta_wrong = theta / (2.0 * np.pi)
            phase_twist = np.exp(-1j * np.dot(theta_wrong, dr_full / L_full))
        else:
            # theta / (2*pi) on the named axis only.
            theta_wrong = theta.copy()
            theta_wrong[axis] = theta_wrong[axis] / (2.0 * np.pi)
            phase_twist = np.exp(-1j * np.dot(theta_wrong, dr_full / L_full))
    elif base_mutator == "M-gauge-3":
        if axis is None:
            # divide by L_folded instead of L_full
            phase_twist = np.exp(
                -1j * np.dot(theta, dr_full / L_folded.astype(np.float64))
            )
        else:
            # divide by L_folded on the named axis only.
            L_denom = L_full.copy()
            L_denom[axis] = L_folded.astype(np.float64)[axis]
            phase_twist = np.exp(-1j * np.dot(theta, dr_full / L_denom))
    elif base_mutator == "M-gauge-5":
        if axis is None:
            # use dr_folded in twist instead of dr_full
            phase_twist = np.exp(-1j * np.dot(theta, dr_folded / L_full))
        else:
            # use dr_folded on the named axis only.
            dr_denom = dr_full.copy()
            dr_denom[axis] = dr_folded[axis]
            phase_twist = np.exp(-1j * np.dot(theta, dr_denom / L_full))
    else:
        phase_twist = np.exp(-1j * np.dot(theta, dr_full / L_full))

    accum = 0.0j
    for kx in range(int(L_folded[0])):
        for ky in range(int(L_folded[1])):
            for kz in range(int(L_folded[2])):
                k_vec = (
                    2.0 * np.pi
                    * np.array([kx, ky, kz], dtype=np.float64) / L_folded
                )
                phase_folded = np.exp(-1j * np.dot(k_vec, dr_folded))
                accum += G_k[kx, ky, kz, aa, bb] * phase_folded
    accum *= phase_twist
    return accum / float(np.prod(L_folded))


def _v36_build_shipping_A_mutated(mutator):
    """Shadow copy of build_slater_orbitals SOC branch with an M-ship
    mutation applied to the shipping A. Returns (A, geom); the caller
    computes G_direct = conj(A) @ A.T at the composite element.

    Supported mutators:
      * ``baseline``  -- reference shipping A (unmutated)
      * ``M-ship-1``  -- phys_dn sign: exp(-i theta.r/L) -> exp(+i theta.r/L)
      * ``M-ship-2``  -- phys_dn unit: pass theta/(2*pi) as theta
      * ``M-ship-3``  -- phys_dn L: divide by L_folded instead of L_phys
      * ``M-ship-4``  -- whole-vector: kf_dot_r offset sign;
                         per-direction: halve offset on that axis
      * ``M-ship-5``  -- kf_dot_r drops sub_offset entirely

    Per-direction ids (``M-ship-{1..5}-{x,y,z}``) apply the mutation
    to only the named axis component
    of theta / L / sub_offset, while the other two axes keep baseline
    behavior. M-ship-4 halves sub_offset on that axis because the
    sign flip is degenerate for ``L_folded=2``; M-ship-5 omits it. See
    ``_axis_of_mutator``. See
    ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    from tools._uhfk_to_mvmc.general_fij_builder import (
        build_pair_list, compute_canonical_reps,
    )
    from tools._uhfk_to_mvmc.occupation_step import step_occupation
    from tools._uhfk_to_mvmc.partner_index import find_partner_rows
    from tools._uhfk_to_mvmc.sublattice_unfold import decode_physical_site

    eigen = np.load(_V36_APBC_EIGEN_PATH)
    occ = np.load(_V36_APBC_OCC_PATH)

    cell_shape = np.array([6, 4, 1], dtype=np.int64)
    subshape = np.array([2, 2, 1], dtype=np.int64)
    Ns_phys = int(np.prod(cell_shape))
    Ncond = 8
    site_positions = np.array(
        [[x, y, z] for z in range(1) for y in range(4) for x in range(6)],
        dtype=np.int64,
    )
    theta_apbc = np.array([np.pi, 0.0, 0.0], dtype=np.float64)

    L_folded = (cell_shape // subshape).astype(np.float64)
    L_phys = cell_shape.astype(np.float64)

    stepped, _ = step_occupation(
        occ["occupation"], eigen["eigenvalue"], occ["column_spin"],
        occ["column_mu_group"], float(occ["T"]),
        ncond_per_group=[Ncond], is_soc_mode=True,
    )
    partner_rows, _ = find_partner_rows(
        eigen["wavevector_index"], theta_apbc, cell_shape // subshape,
    )
    canonical, _ = compute_canonical_reps(
        partner_rows, eigen["wavevector_index"],
    )
    pair_list = build_pair_list(
        stepped, occ["column_spin"], canonical, partner_rows,
        is_soc_mode=True,
    )

    wavevector_index = eigen["wavevector_index"].astype(np.int64)
    eigenvector = eigen["eigenvector"].astype(np.complex128)
    nvol_folded = wavevector_index.shape[0]
    subvol = int(np.prod(subshape))
    inv_sqrt = 1.0 / np.sqrt(float(nvol_folded))

    folded_cell = np.empty((Ns_phys, 3), dtype=np.int64)
    sub_offset_per_site = np.empty((Ns_phys, 3), dtype=np.int64)
    folded_orb_per_site = np.empty(Ns_phys, dtype=np.int64)
    for i in range(Ns_phys):
        fc, so = decode_physical_site(site_positions[i], subshape)
        folded_cell[i] = fc
        sub_offset_per_site[i] = so
        folded_orb_per_site[i] = (
            so[0] + subshape[0] * (so[1] + subshape[1] * so[2])
        )

    # Per-direction mutators (M-ship-N-{x,y,z}) restrict the base
    # mutation N to a single axis; whole-vector ids (axis is None)
    # keep their original code paths below untouched.
    axis = _axis_of_mutator(mutator)
    base_mutator = (
        mutator if axis is None else "-".join(mutator.split("-")[:-1])
    )

    # phys_dn variants for M-ship-1..3.
    if base_mutator == "M-ship-1":
        if axis is None:
            phys_arg = np.einsum(
                "d,id->i", theta_apbc / L_phys,
                site_positions.astype(np.float64),
            )
            phys_dn = np.exp(+1j * phys_arg)  # sign flipped
        else:
            # flip sign of theta on the named axis only.
            theta_mut = theta_apbc.copy()
            theta_mut[axis] = -theta_mut[axis]
            phys_arg = np.einsum(
                "d,id->i", theta_mut / L_phys,
                site_positions.astype(np.float64),
            )
            phys_dn = np.exp(-1j * phys_arg)
    elif base_mutator == "M-ship-2":
        if axis is None:
            theta_wrong = theta_apbc / (2.0 * np.pi)
            phys_arg = np.einsum(
                "d,id->i", theta_wrong / L_phys,
                site_positions.astype(np.float64),
            )
            phys_dn = np.exp(-1j * phys_arg)
        else:
            # theta / (2*pi) on the named axis only.
            theta_wrong = theta_apbc.copy()
            theta_wrong[axis] = theta_wrong[axis] / (2.0 * np.pi)
            phys_arg = np.einsum(
                "d,id->i", theta_wrong / L_phys,
                site_positions.astype(np.float64),
            )
            phys_dn = np.exp(-1j * phys_arg)
    elif base_mutator == "M-ship-3":
        if axis is None:
            phys_arg = np.einsum(
                "d,id->i", theta_apbc / L_folded,
                site_positions.astype(np.float64),
            )
            phys_dn = np.exp(-1j * phys_arg)
        else:
            # divide by L_folded on the named axis only.
            L_denom = L_phys.copy()
            L_denom[axis] = L_folded[axis]
            phys_arg = np.einsum(
                "d,id->i", theta_apbc / L_denom,
                site_positions.astype(np.float64),
            )
            phys_dn = np.exp(-1j * phys_arg)
    else:
        phys_arg = np.einsum(
            "d,id->i", theta_apbc / L_phys, site_positions.astype(np.float64),
        )
        phys_dn = np.exp(-1j * phys_arg)

    # kf_dot_r variants for M-ship-4/5.
    k_folded_all = (
        2.0 * np.pi * wavevector_index.astype(np.float64) / L_folded
    )
    if base_mutator == "M-ship-4":
        if axis is None:
            kf_dot_r = np.einsum(
                "kd,id->ki", k_folded_all,
                (folded_cell - sub_offset_per_site).astype(np.float64),
            )
        else:
            # Halve sub_offset on the named axis only;
            # other axes keep the baseline contribution.
            pos_arg = (folded_cell + sub_offset_per_site).astype(np.float64)
            pos_arg[:, axis] = (
                folded_cell[:, axis].astype(np.float64)
                + 0.5 * sub_offset_per_site[:, axis].astype(np.float64)
            )
            kf_dot_r = np.einsum("kd,id->ki", k_folded_all, pos_arg)
    elif base_mutator == "M-ship-5":
        if axis is None:
            kf_dot_r = np.einsum(
                "kd,id->ki", k_folded_all, folded_cell.astype(np.float64),
            )
        else:
            # drop sub_offset on the named axis only.
            pos_arg = (folded_cell + sub_offset_per_site).astype(np.float64)
            pos_arg[:, axis] = folded_cell[:, axis].astype(np.float64)
            kf_dot_r = np.einsum("kd,id->ki", k_folded_all, pos_arg)
    else:
        kf_dot_r = np.einsum(
            "kd,id->ki", k_folded_all,
            (folded_cell + sub_offset_per_site).astype(np.float64),
        )

    plane_wave_soc = (
        np.exp(-1j * kf_dot_r) * inv_sqrt * phys_dn[np.newaxis, :]
    )

    A_soc = np.zeros(
        (2 * Ns_phys, 2 * len(pair_list)), dtype=np.complex128
    )
    for p, pair in enumerate(pair_list):
        for column_idx, member in enumerate(pair):
            k_row_member, alpha = int(member[0]), int(member[1])
            slater_col = 2 * p + column_idx
            plane_wave_member = plane_wave_soc[k_row_member]
            for spin in (0, 1):
                hwave_rows = (
                    2 * folded_orb_per_site + spin
                ).astype(np.int64)
                v_vals = np.array(
                    [
                        eigenvector[k_row_member, int(hwave_rows[i]), alpha]
                        for i in range(Ns_phys)
                    ],
                    dtype=np.complex128,
                )
                amp = plane_wave_member * v_vals
                row_start = spin * Ns_phys
                A_soc[row_start:row_start + Ns_phys, slater_col] = amp
    return A_soc, {"Ns": Ns_phys}


_M_GAUGE_IDS = tuple(f"M-gauge-{k}" for k in range(1, 6))
_M_SHIP_IDS = tuple(f"M-ship-{k}" for k in range(1, 6))


@pytest.mark.parametrize("mutator", _M_GAUGE_IDS + _M_SHIP_IDS)
def test_gauge_lift_apbc_mutation_regression_matrix(mutator):
    """For each of the 10 named
    mutations, apply to the shipping ``gauge_lift`` (M-gauge) or
    ``build_slater_orbitals`` (M-ship), evaluated at the fixture's
    composite element ``(i_c, s_c, j_c, t_c)`` from the manifest,
    and assert ``|G_mut - G_base| >= max(1e-5, 0.10 * |G_base|)``.

    Runs against the tracked snapshot (tests/data/v36_case_soc_rashba_2d
    _sub_apbc_*.npz) so the pin remains reproducible without a live SCF
    workspace. This is a unit test on stable inputs; the G4
    workspace gate runs the same mutations against a fresh SCF via
    soc_apbc_topology_guard.py. See
    ``docs/en/source/uhfk/tools/uhfk_to_mvmc.rst``.
    """
    import json

    with open(_V36_APBC_COMPOSITE_PATH) as fp:
        manifest = json.load(fp)
    i_c = int(manifest["i_c"])
    s_c = int(manifest["s_c"])
    j_c = int(manifest["j_c"])
    t_c = int(manifest["t_c"])

    gp = np.load(_V36_APBC_GREEN_PATH)
    green_sublattice = gp["green_sublattice"]

    cell_shape = np.array([6, 4, 1], dtype=np.int64)
    subshape = np.array([2, 2, 1], dtype=np.int64)
    Ns = int(np.prod(cell_shape))
    site_positions = np.array(
        [[x, y, z] for z in range(1) for y in range(4) for x in range(6)],
        dtype=np.int64,
    )
    theta_apbc = np.array([np.pi, 0.0, 0.0], dtype=np.float64)
    all_i = i_c + s_c * Ns
    all_j = j_c + t_c * Ns

    if mutator.startswith("M-gauge-"):
        G_base = _v36_gauge_lift_at_composite_mutated(
            green_sublattice, i_c, s_c, j_c, t_c,
            subshape=subshape, cell_shape=cell_shape,
            site_positions=site_positions,
            boundary_theta=theta_apbc, mutator="baseline",
        )
        G_mut = _v36_gauge_lift_at_composite_mutated(
            green_sublattice, i_c, s_c, j_c, t_c,
            subshape=subshape, cell_shape=cell_shape,
            site_positions=site_positions,
            boundary_theta=theta_apbc, mutator=mutator,
        )
    else:  # M-ship-*
        A_base, _ = _v36_build_shipping_A_mutated("baseline")
        A_mut, _ = _v36_build_shipping_A_mutated(mutator)
        G_direct_base = np.conj(A_base) @ A_base.T
        G_direct_mut = np.conj(A_mut) @ A_mut.T
        G_base = G_direct_base[all_i, all_j]
        G_mut = G_direct_mut[all_i, all_j]

    G_base_abs = abs(G_base)
    delta_M = abs(G_mut - G_base)
    T_M = max(1e-5, 0.10 * G_base_abs)
    assert delta_M >= T_M, (
        f"{mutator}: delta_M = {delta_M:.3e} < T_M = {T_M:.3e} "
        f"(|G_base| = {G_base_abs:.3e}); the mutation did NOT perturb "
        "the composite element above the sensitivity floor. The "
        "shipping code path this mutation targets is therefore not being "
        "exercised at the composite; either regenerate the composite via "
        "producer with a stronger candidate, or investigate why "
        "the mutation is orthogonal to the composite's dr_folded / "
        "sub_offset / twist coordinates."
    )


# ---------------------------------------------------------------------
# A2 positive pin over the four multi-direction APBC fixtures.
# ---------------------------------------------------------------------


@pytest.mark.parametrize("fixture_name,expected_theta,expected_ncond", [
    ("xy",  (np.pi, np.pi, 0.0), 20),
    ("xz",  (np.pi, 0.0, np.pi), 20),
    ("yz",  (0.0, np.pi, np.pi), 24),
    ("xyz", (np.pi, np.pi, np.pi), 12),
])
def test_gauge_lift_apbc_matches_shipping_A_on_real_run_v37(
    fixture_name, expected_theta, expected_ncond,
):
    """Parametrized over the four shipping fixtures under
    CellShape=[4,4,4]/SubShape=[2,2,2]. For
    each fixture, verify the shipping A density matches the gauge-
    lifted green_sublattice at max_abs_delta < 1e-10 under the
    fixture's own boundary_theta -- the G1 in-process self-consistency
    gate extended to multi-direction APBC.

    The test above pins the single-direction (AP-P-P) APBC case at
    1e-10. This test covers multi-direction APBC:
    xy (AP-AP-P), xz (AP-P-AP), yz (P-AP-AP), xyz (AP-AP-AP).
    See ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    from tools._uhfk_to_mvmc.density_check import gauge_lift
    from tools._uhfk_to_mvmc.general_fij_builder import (
        build_slater_orbitals, build_pair_list, compute_canonical_reps,
    )
    from tools._uhfk_to_mvmc.occupation_step import step_occupation
    from tools._uhfk_to_mvmc.partner_index import find_partner_rows

    prefix = f"v37_case_soc_rashba_3d_sub_apbc_{fixture_name}"
    gp = np.load(str(_DATA_DIR / f"{prefix}_green.npz"))
    eigen = np.load(str(_DATA_DIR / f"{prefix}_eigen.npz"))
    occ = np.load(str(_DATA_DIR / f"{prefix}_occupation.npz"))

    green_sublattice = gp["green_sublattice"]
    cell_shape = np.array([4, 4, 4], dtype=np.int64)
    subshape = np.array([2, 2, 2], dtype=np.int64)
    Ns = int(np.prod(cell_shape))  # 64
    site_positions = np.array(
        [[x, y, z] for z in range(4) for y in range(4) for x in range(4)],
        dtype=np.int64,
    )
    theta = np.array(expected_theta, dtype=np.float64)

    # Build shipping A.
    stepped, _ = step_occupation(
        occ["occupation"], eigen["eigenvalue"], occ["column_spin"],
        occ["column_mu_group"], float(occ["T"]),
        ncond_per_group=[expected_ncond], is_soc_mode=True,
    )
    partner_rows, _ = find_partner_rows(
        eigen["wavevector_index"], theta, cell_shape // subshape,
    )
    canonical, _ = compute_canonical_reps(
        partner_rows, eigen["wavevector_index"],
    )
    pair_list = build_pair_list(
        stepped, occ["column_spin"], canonical, partner_rows,
        is_soc_mode=True,
    )
    A_ship = build_slater_orbitals(
        wavevector_index=eigen["wavevector_index"],
        eigenvector=eigen["eigenvector"],
        column_spin=occ["column_spin"],
        site_positions=site_positions,
        cell_shape=cell_shape,
        subshape=subshape,
        theta=theta,
        pair_list=pair_list,
        is_soc_mode=True,
    )
    G_ship = np.conj(A_ship) @ A_ship.T

    # Lift green_sublattice for every (i, s, j, t) element.
    folded_cell_of = lambda r_phys: r_phys // subshape
    G_lift = np.zeros((2 * Ns, 2 * Ns), dtype=np.complex128)
    for site_i in range(Ns):
        for spin_i in (0, 1):
            for site_j in range(Ns):
                for spin_j in (0, 1):
                    all_i = site_i + spin_i * Ns
                    all_j = site_j + spin_j * Ns
                    G_lift[all_i, all_j] = gauge_lift(
                        green_sublattice, site_i, spin_i, site_j, spin_j,
                        subshape=subshape, cell_shape=cell_shape,
                        site_positions=site_positions,
                        folded_cell_of=folded_cell_of,
                        boundary_theta=theta,
                    )

    max_diff = float(np.max(np.abs(G_ship - G_lift)))
    cross_spin_block = G_lift[:Ns, Ns:]
    cross_max = float(np.max(np.abs(cross_spin_block)))
    print(
        f"[{fixture_name}] max_abs_delta={max_diff:.3e} "
        f"cross_spin_max={cross_max:.3e}"
    )

    assert max_diff < 1e-10, (
        f"[{fixture_name}] |G_ship - G_lift|_max = {max_diff:.3e} > 1e-10 "
        "under multi-direction APBC; gauge_lift + build_slater_orbitals "
        "do not agree on this fixture."
    )

    # Cross-spin sanity: fixture MUST exercise Rashba SOC off-diagonals.
    assert cross_max > 0.01, (
        f"[{fixture_name}] |G_lift[s!=t]|_max = {cross_max:.3e} <= 0.01; "
        "the fixture does not exercise Rashba SOC off-diagonals as "
        "expected."
    )


# ---------------------------------------------------------------------
# Negative control over the four multi-direction APBC fixtures.
# ---------------------------------------------------------------------


@pytest.mark.parametrize("fixture_name,expected_theta,expected_ncond", [
    ("xy",  (np.pi, np.pi, 0.0), 20),
    ("xz",  (np.pi, 0.0, np.pi), 20),
    ("yz",  (0.0, np.pi, np.pi), 24),
    ("xyz", (np.pi, np.pi, np.pi), 12),
])
def test_gauge_lift_apbc_negative_control_theta_zero_v37(
    fixture_name, expected_theta, expected_ncond,
):
    """Parametrized negative control over the four shipping fixtures.
    At each fixture's composite element
    (i_c, s_c, j_c, t_c) from composite_element.json:

    - G_A_c    = shipping A density under theta = expected_theta
    - G_neg_c  = gauge_lift with boundary_theta = (0, 0, 0) at same composite

    Assert |G_neg_c - G_A_c| >= 0.3 * |G_A_c|. Proves the multi-
    direction twist gauge is functionally load-bearing on this fixture;
    if the gauge were vestigial or already carried elsewhere, dropping
    it would leave G_neg_c ~ G_A_c at machine precision. See
    ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    import json
    from tools._uhfk_to_mvmc.density_check import gauge_lift
    from tools._uhfk_to_mvmc.general_fij_builder import (
        build_slater_orbitals, build_pair_list, compute_canonical_reps,
    )
    from tools._uhfk_to_mvmc.occupation_step import step_occupation
    from tools._uhfk_to_mvmc.partner_index import find_partner_rows

    prefix = f"v37_case_soc_rashba_3d_sub_apbc_{fixture_name}"

    with open(str(_DATA_DIR / f"{prefix}_composite_element.json")) as fp:
        manifest = json.load(fp)
    i_c = int(manifest["i_c"])
    s_c = int(manifest["s_c"])
    j_c = int(manifest["j_c"])
    t_c = int(manifest["t_c"])

    gp = np.load(str(_DATA_DIR / f"{prefix}_green.npz"))
    eigen = np.load(str(_DATA_DIR / f"{prefix}_eigen.npz"))
    occ = np.load(str(_DATA_DIR / f"{prefix}_occupation.npz"))

    green_sublattice = gp["green_sublattice"]
    cell_shape = np.array([4, 4, 4], dtype=np.int64)
    subshape = np.array([2, 2, 2], dtype=np.int64)
    Ns = int(np.prod(cell_shape))  # 64
    site_positions = np.array(
        [[x, y, z] for z in range(4) for y in range(4) for x in range(4)],
        dtype=np.int64,
    )
    theta = np.array(expected_theta, dtype=np.float64)

    # Build shipping A + G_A_c at composite.
    stepped, _ = step_occupation(
        occ["occupation"], eigen["eigenvalue"], occ["column_spin"],
        occ["column_mu_group"], float(occ["T"]),
        ncond_per_group=[expected_ncond], is_soc_mode=True,
    )
    partner_rows, _ = find_partner_rows(
        eigen["wavevector_index"], theta, cell_shape // subshape,
    )
    canonical, _ = compute_canonical_reps(
        partner_rows, eigen["wavevector_index"],
    )
    pair_list = build_pair_list(
        stepped, occ["column_spin"], canonical, partner_rows,
        is_soc_mode=True,
    )
    A_ship = build_slater_orbitals(
        wavevector_index=eigen["wavevector_index"],
        eigenvector=eigen["eigenvector"],
        column_spin=occ["column_spin"],
        site_positions=site_positions,
        cell_shape=cell_shape,
        subshape=subshape,
        theta=theta,
        pair_list=pair_list,
        is_soc_mode=True,
    )
    G_direct = np.conj(A_ship) @ A_ship.T
    all_i = i_c + s_c * Ns
    all_j = j_c + t_c * Ns
    G_A_c = G_direct[all_i, all_j]

    folded_cell_of = lambda r_phys: r_phys // subshape
    G_neg_c = gauge_lift(
        green_sublattice, i_c, s_c, j_c, t_c,
        subshape=subshape, cell_shape=cell_shape,
        site_positions=site_positions,
        folded_cell_of=folded_cell_of,
        boundary_theta=(0.0, 0.0, 0.0),  # twist gauge intentionally dropped
    )

    G_A_c_abs = abs(G_A_c)
    assert G_A_c_abs > 1e-3, (
        f"[{fixture_name}] G_A_c magnitude {G_A_c_abs:.3e} <= 1e-3 at "
        f"composite (i={i_c}, s={s_c}, j={j_c}, t={t_c}); the composite "
        "fell off the manifest; check the snapshot / manifest pairing."
    )
    delta = abs(G_neg_c - G_A_c)
    threshold = 0.3 * G_A_c_abs
    assert delta >= threshold, (
        f"[{fixture_name}] |G_neg_c - G_A_c| = {delta:.3e} < "
        f"0.3 * |G_A_c| = {threshold:.3e}; dropping the twist gauge did "
        "NOT perturb the composite element by the expected Rashba-scale "
        "amount. Either the twist gauge is vestigial in gauge_lift, or "
        "the composite was accidentally chosen with a coordinate "
        "invariant to APBC."
    )


# ---------------------------------------------------------------------
# Per-direction mutation matrix over the four multi-direction APBC fixtures.
# It parallels ``test_gauge_lift_apbc_mutation_regression_matrix`` above,
# reuses the same mutation semantics, and reads the composite element and
# per-mutation floor T_M from each fixture's own composite_element.json
# instead of recomputing T_M = max(1e-5, 0.10 * |G_c|) locally: the
# producer's floor preserves the one documented relaxation on top of
# that base formula: an inactive axis (theta_axis == 0) has T_M = 0 and
# is handled here via pytest.skip. Every active-axis entry
# has a positive threshold; the producer rejects a structurally
# degenerate active-axis mutation set instead of emitting a vacuous gate.
# ---------------------------------------------------------------------


_V37_MUTATION_KEYS = tuple(
    f"{base_id}-{axis_char}"
    for base_id in (_M_GAUGE_IDS + _M_SHIP_IDS)
    for axis_char in ("x", "y", "z")
)


def _v37_build_shipping_A_mutated(
    mutator, *, eigen, occ, cell_shape, subshape, site_positions, theta,
    ncond,
):
    """Parametrized sibling of ``_v36_build_shipping_A_mutated`` for the
    four multi-direction APBC fixtures. Same mutation
    semantics and the same shadow-copy-of-``build_slater_orbitals``
    approach; ``cell_shape`` / ``subshape`` / ``site_positions`` /
    ``theta`` / ``ncond`` and the eigen / occupation snapshot arrays
    (already-loaded ``np.load`` objects) are parameters instead of the
    single-direction fixture's hardcoded constants, so this helper is reusable
    across xy/xz/yz/xyz without touching ``_v36_build_shipping_A_
    mutated``. See that function's docstring for the mutator list;
    ``mutator`` accepts both whole-vector ids and
    per-direction ids per ``_axis_of_mutator``.
    """
    from tools._uhfk_to_mvmc.general_fij_builder import (
        build_pair_list, compute_canonical_reps,
    )
    from tools._uhfk_to_mvmc.occupation_step import step_occupation
    from tools._uhfk_to_mvmc.partner_index import find_partner_rows
    from tools._uhfk_to_mvmc.sublattice_unfold import decode_physical_site

    Ns_phys = int(np.prod(cell_shape))
    theta = np.asarray(theta, dtype=np.float64)
    L_folded = (cell_shape // subshape).astype(np.float64)
    L_phys = cell_shape.astype(np.float64)

    stepped, _ = step_occupation(
        occ["occupation"], eigen["eigenvalue"], occ["column_spin"],
        occ["column_mu_group"], float(occ["T"]),
        ncond_per_group=[ncond], is_soc_mode=True,
    )
    partner_rows, _ = find_partner_rows(
        eigen["wavevector_index"], theta, cell_shape // subshape,
    )
    canonical, _ = compute_canonical_reps(
        partner_rows, eigen["wavevector_index"],
    )
    pair_list = build_pair_list(
        stepped, occ["column_spin"], canonical, partner_rows,
        is_soc_mode=True,
    )

    wavevector_index = eigen["wavevector_index"].astype(np.int64)
    eigenvector = eigen["eigenvector"].astype(np.complex128)
    nvol_folded = wavevector_index.shape[0]
    inv_sqrt = 1.0 / np.sqrt(float(nvol_folded))

    folded_cell = np.empty((Ns_phys, 3), dtype=np.int64)
    sub_offset_per_site = np.empty((Ns_phys, 3), dtype=np.int64)
    folded_orb_per_site = np.empty(Ns_phys, dtype=np.int64)
    for i in range(Ns_phys):
        fc, so = decode_physical_site(site_positions[i], subshape)
        folded_cell[i] = fc
        sub_offset_per_site[i] = so
        folded_orb_per_site[i] = (
            so[0] + subshape[0] * (so[1] + subshape[1] * so[2])
        )

    # Per-direction mutators (M-ship-N-{x,y,z}) restrict the base
    # mutation N to a single axis; whole-vector ids (axis is
    # None) keep their original code paths below untouched.
    axis = _axis_of_mutator(mutator)
    base_mutator = (
        mutator if axis is None else "-".join(mutator.split("-")[:-1])
    )

    # phys_dn variants for M-ship-1..3.
    if base_mutator == "M-ship-1":
        if axis is None:
            phys_arg = np.einsum(
                "d,id->i", theta / L_phys,
                site_positions.astype(np.float64),
            )
            phys_dn = np.exp(+1j * phys_arg)  # sign flipped
        else:
            theta_mut = theta.copy()
            theta_mut[axis] = -theta_mut[axis]
            phys_arg = np.einsum(
                "d,id->i", theta_mut / L_phys,
                site_positions.astype(np.float64),
            )
            phys_dn = np.exp(-1j * phys_arg)
    elif base_mutator == "M-ship-2":
        if axis is None:
            theta_wrong = theta / (2.0 * np.pi)
            phys_arg = np.einsum(
                "d,id->i", theta_wrong / L_phys,
                site_positions.astype(np.float64),
            )
            phys_dn = np.exp(-1j * phys_arg)
        else:
            theta_wrong = theta.copy()
            theta_wrong[axis] = theta_wrong[axis] / (2.0 * np.pi)
            phys_arg = np.einsum(
                "d,id->i", theta_wrong / L_phys,
                site_positions.astype(np.float64),
            )
            phys_dn = np.exp(-1j * phys_arg)
    elif base_mutator == "M-ship-3":
        if axis is None:
            phys_arg = np.einsum(
                "d,id->i", theta / L_folded,
                site_positions.astype(np.float64),
            )
            phys_dn = np.exp(-1j * phys_arg)
        else:
            L_denom = L_phys.copy()
            L_denom[axis] = L_folded[axis]
            phys_arg = np.einsum(
                "d,id->i", theta / L_denom,
                site_positions.astype(np.float64),
            )
            phys_dn = np.exp(-1j * phys_arg)
    else:
        phys_arg = np.einsum(
            "d,id->i", theta / L_phys, site_positions.astype(np.float64),
        )
        phys_dn = np.exp(-1j * phys_arg)

    # kf_dot_r variants for M-ship-4/5.
    k_folded_all = (
        2.0 * np.pi * wavevector_index.astype(np.float64) / L_folded
    )
    if base_mutator == "M-ship-4":
        if axis is None:
            kf_dot_r = np.einsum(
                "kd,id->ki", k_folded_all,
                (folded_cell - sub_offset_per_site).astype(np.float64),
            )
        else:
            pos_arg = (
                folded_cell + sub_offset_per_site
            ).astype(np.float64)
            pos_arg[:, axis] = (
                folded_cell[:, axis].astype(np.float64)
                + 0.5 * sub_offset_per_site[:, axis].astype(np.float64)
            )
            kf_dot_r = np.einsum("kd,id->ki", k_folded_all, pos_arg)
    elif base_mutator == "M-ship-5":
        if axis is None:
            kf_dot_r = np.einsum(
                "kd,id->ki", k_folded_all, folded_cell.astype(np.float64),
            )
        else:
            pos_arg = (
                folded_cell + sub_offset_per_site
            ).astype(np.float64)
            pos_arg[:, axis] = folded_cell[:, axis].astype(np.float64)
            kf_dot_r = np.einsum("kd,id->ki", k_folded_all, pos_arg)
    else:
        kf_dot_r = np.einsum(
            "kd,id->ki", k_folded_all,
            (folded_cell + sub_offset_per_site).astype(np.float64),
        )

    plane_wave_soc = (
        np.exp(-1j * kf_dot_r) * inv_sqrt * phys_dn[np.newaxis, :]
    )

    A_soc = np.zeros(
        (2 * Ns_phys, 2 * len(pair_list)), dtype=np.complex128
    )
    for p, pair in enumerate(pair_list):
        for column_idx, member in enumerate(pair):
            k_row_member, alpha = int(member[0]), int(member[1])
            slater_col = 2 * p + column_idx
            plane_wave_member = plane_wave_soc[k_row_member]
            for spin in (0, 1):
                hwave_rows = (
                    2 * folded_orb_per_site + spin
                ).astype(np.int64)
                v_vals = np.array(
                    [
                        eigenvector[k_row_member, int(hwave_rows[i]), alpha]
                        for i in range(Ns_phys)
                    ],
                    dtype=np.complex128,
                )
                amp = plane_wave_member * v_vals
                row_start = spin * Ns_phys
                A_soc[row_start:row_start + Ns_phys, slater_col] = amp
    return A_soc, {"Ns": Ns_phys}


@pytest.mark.parametrize("mutation_key", _V37_MUTATION_KEYS)
@pytest.mark.parametrize("fixture_name,expected_theta,expected_ncond", [
    ("xy",  (np.pi, np.pi, 0.0), 20),
    ("xz",  (np.pi, 0.0, np.pi), 20),
    ("yz",  (0.0, np.pi, np.pi), 24),
    ("xyz", (np.pi, np.pi, np.pi), 12),
])
def test_gauge_lift_apbc_mutation_regression_matrix_v37(
    fixture_name, expected_theta, expected_ncond, mutation_key,
):
    """Parallel to ``test_gauge_lift_apbc_mutation_regression_matrix``:
    the same 30-id per-direction mutation matrix
    (M-{gauge,ship}-{1..5}-{x,y,z}) is parametrized over the four
    multi-direction APBC fixtures (xy/xz/yz/xyz on
    CellShape=[4,4,4]/SubShape=[2,2,2]) for 4 x 30 = 120 cases.

    A mutation whose axis carries no APBC twist on this fixture
    (theta_axis == 0) is a no-op by construction, so the
    case is skipped rather than asserted. Otherwise the case loads the
    fixture's composite_element.json manifest, applies the
    mutation to the shipping gauge_lift (M-gauge) or
    build_slater_orbitals (M-ship) at the manifest's composite element
    (i_c, s_c, j_c, t_c), and asserts the resulting delta clears the
    manifest's own T_M_per_mutation floor for that mutation id (see
    the module-level banner above for why T_M is read from the
    manifest rather than recomputed).

    Runs against the tracked snapshot (tests/data/v37_case_soc_rashba
    _3d_sub_apbc_{fixture}_*.npz) so the pin remains reproducible
    without a live SCF workspace. See
    ``docs/en/source/uhfk/tools/uhfk_to_mvmc.rst``.
    """
    import json

    axis = _axis_of_mutator(mutation_key)
    axis_char = "xyz"[axis]
    theta = np.array(expected_theta, dtype=np.float64)

    if theta[axis] == 0.0:
        pytest.skip(
            f"{mutation_key}: theta_{axis_char} = 0 on fixture "
            f"{fixture_name!r}; this axis carries no APBC twist, so "
            "the mutation is a no-op by construction."
        )

    prefix = f"v37_case_soc_rashba_3d_sub_apbc_{fixture_name}"
    with open(str(_DATA_DIR / f"{prefix}_composite_element.json")) as fp:
        manifest = json.load(fp)
    i_c = int(manifest["i_c"])
    s_c = int(manifest["s_c"])
    j_c = int(manifest["j_c"])
    t_c = int(manifest["t_c"])
    T_M = float(manifest["T_M_per_mutation"][mutation_key])

    cell_shape = np.array([4, 4, 4], dtype=np.int64)
    subshape = np.array([2, 2, 2], dtype=np.int64)
    Ns = int(np.prod(cell_shape))
    site_positions = np.array(
        [[x, y, z] for z in range(4) for y in range(4) for x in range(4)],
        dtype=np.int64,
    )
    all_i = i_c + s_c * Ns
    all_j = j_c + t_c * Ns

    if mutation_key.startswith("M-gauge-"):
        gp = np.load(str(_DATA_DIR / f"{prefix}_green.npz"))
        green_sublattice = gp["green_sublattice"]
        G_base = _v36_gauge_lift_at_composite_mutated(
            green_sublattice, i_c, s_c, j_c, t_c,
            subshape=subshape, cell_shape=cell_shape,
            site_positions=site_positions,
            boundary_theta=theta, mutator="baseline",
        )
        G_mut = _v36_gauge_lift_at_composite_mutated(
            green_sublattice, i_c, s_c, j_c, t_c,
            subshape=subshape, cell_shape=cell_shape,
            site_positions=site_positions,
            boundary_theta=theta, mutator=mutation_key,
        )
    else:  # M-ship-*
        eigen = np.load(str(_DATA_DIR / f"{prefix}_eigen.npz"))
        occ = np.load(str(_DATA_DIR / f"{prefix}_occupation.npz"))
        A_base, _ = _v37_build_shipping_A_mutated(
            "baseline", eigen=eigen, occ=occ, cell_shape=cell_shape,
            subshape=subshape, site_positions=site_positions,
            theta=theta, ncond=expected_ncond,
        )
        A_mut, _ = _v37_build_shipping_A_mutated(
            mutation_key, eigen=eigen, occ=occ, cell_shape=cell_shape,
            subshape=subshape, site_positions=site_positions,
            theta=theta, ncond=expected_ncond,
        )
        G_direct_base = np.conj(A_base) @ A_base.T
        G_direct_mut = np.conj(A_mut) @ A_mut.T
        G_base = G_direct_base[all_i, all_j]
        G_mut = G_direct_mut[all_i, all_j]

    G_base_abs = abs(G_base)
    delta_M = abs(G_mut - G_base)
    print(
        f"[{fixture_name}] {mutation_key}: delta_M={delta_M:.3e} "
        f"T_M={T_M:.3e} |G_base|={G_base_abs:.3e}"
    )
    assert delta_M >= T_M, (
        f"[{fixture_name}] {mutation_key}: delta_M = {delta_M:.3e} < "
        f"T_M = {T_M:.3e} (|G_base| = {G_base_abs:.3e}); the mutation "
        "did NOT perturb the composite element above the manifest's "
        "producer-time sensitivity floor. The shipping code path this "
        "mutation targets is therefore not being exercised at the "
        "composite; either regenerate the composite via the "
        "producer, or investigate why the mutation is orthogonal to "
        "the composite's dr_folded / sub_offset / twist coordinates."
    )
