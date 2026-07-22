"""Shared target-occupation assertion helper.

Both pytest and the mVMC E2E validation harness (`run.sh`) invoke this
helper immediately after H-wave SCF, BEFORE the bridge runs. Target
values live only in this module (module-level `_CASE_TARGETS` constant)
so a single edit propagates everywhere.

CLI usage:
    python3 scripts/assert_occupation.py <case_dir> <work_dir>

Programmatic usage:
    from assert_occupation import assert_case_occupation
    assert_case_occupation(case_dir_path, work_dir_path)
"""
from __future__ import annotations

import os
import sys

import numpy as np


# Single source of truth for every case's target. Two supported target
# formats:
#
#   1. Sz-fixed / Sz-free without SOC (column_spin in {0, 1}):
#      {"up": [tuple(wv), ...], "down": [tuple(wv), ...]}
#
#   2. SOC / enable_spin_orbital=true (column_spin == -1, spin packed
#      into orbital label): {"n_per_k": {tuple(wv): int_count, ...}}.
#      No up/down bookkeeping under SOC because Rashba mixes spin.
#
# Keys must match the fixture directory names under
# tests/validation/uhfk_mvmc_pairproduct/.
_CASE_TARGETS: dict[str, dict] = {
    "case_pbc_sz2": {
        # Target FM UHF ground state at PBC L=8, U small enough to
        # remain paramagnetic-like: 3 up (at k=0, ±2π/8=±π/4) + 1 down
        # (at k=0). The target has a cross pair at self-pair k=0
        # + same-spin up-up excess at non-self (+π/4, -π/4) block.
        "up": [(0, 0, 0), (1, 0, 0), (-1, 0, 0)],
        "down": [(0, 0, 0)],
    },
    "case_zeeman_sz_free": {
        # 1D L=8 PBC Hubbard with Ne=4, no 2Sz constraint, flag_fock=false.
        # H-wave's uhfk drops spin-block Transfer
        # indices when enable_spin_orbital=false (uhfk.py:1136-1146) and
        # the Hermite check hard-errors on them (uhfk.py:879). The
        # bridge in turn rejects enable_spin_orbital=true, so the Zeeman
        # entries could not survive in the fixture; the file carries only
        # the standard NN hopping. H-wave then converges to a
        # symmetry-broken FM state selected by the initial random Green
        # seed: Ne_up=3 at k=0 and ±π/4, Ne_down=1 at k=0 — the same
        # target occupation as case_pbc_sz2 but reached via the Sz-free
        # (single-mu-group, no 2Sz) branch. Canonical blocks:
        #   self k=0:    NN_up=1, NN_down=1 → n_cross=1, excess=0/0 ✓
        #   self k=π:    empty ✓
        #   (+1,-1):     NN_up_k=NN_up_p=1, NN_down=0 → excess_up_k=
        #                 excess_up_p=1 (same-spin up-up pair) ✓
        #   (+2,-2):     empty ✓
        #   (+3,-3):     empty ✓
        # Ne = 4 (even). See
        # ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
        "up": [(0, 0, 0), (1, 0, 0), (-1, 0, 0)],
        "down": [(0, 0, 0)],
    },
    "case_soc_rashba_2d_nosub": {
        # 2D 4x4 Rashba + Hubbard SOC (no sublattice
        # folding). CellShape = [4, 4, 1], SubShape = [1, 1, 1], Ncond = 6.
        # Under SOC column_spin is packed (== -1), so target is the total
        # spin-summed occupation per folded k-row; the SCF at t=1,
        # alpha=0.5, U=2 fills k=(0,0,0) with two carriers and the four
        # nearest shells {(0,+1,0),(0,-1,0),(+1,0,0),(-1,0,0)} with one
        # each. Sum = 6 = Ncond.
        "n_per_k": {
            (0, 0, 0): 2,
            (0, 1, 0): 1,
            (0, -1, 0): 1,
            (1, 0, 0): 1,
            (-1, 0, 0): 1,
        },
    },
    "case_soc_rashba_2d_nosub_apbc": {
        # 2D 4x4 Rashba + Hubbard SOC + APBC in x (no
        # sublattice folding). CellShape = [4, 4, 1], SubShape = [1, 1, 1],
        # Ncond = 8, BoundaryCondition = ["antiperiodic", "periodic",
        # "periodic"]. The APBC-in-x twist shifts the k-mesh so the
        # Ncond = 6 gap closes at the Fermi level (see
        # scripts/soc_pilot_gap.py); Ncond = 8 sits inside the next
        # HOMO-LUMO gap (~1.35). Under the twisted mesh the SCF fills
        # k_x in {0, -1} × k_y in {-1, 0, +1} with (0, 0, 0) and
        # (-1, 0, 0) each carrying two carriers and the four flanking
        # k rows carrying one each. Sum = 8 = Ncond.
        "n_per_k": {
            (0, 0, 0): 2,
            (0, 1, 0): 1,
            (0, -1, 0): 1,
            (-1, 0, 0): 2,
            (-1, 1, 0): 1,
            (-1, -1, 0): 1,
        },
    },
    "case_soc_rashba_2d_sub": {
        # 2D 6x4 Rashba + Hubbard SOC with SubShape = [2, 2, 1] (folded
        # BZ [3, 2, 1], Ncond = 8). column_spin packed under SOC.
        # Fills the two self-pair k rows (k_x = 0) with two carriers each
        # and each of the four non-self k rows (k_x = ±1) with one carrier.
        # Sum = 2*2 + 4*1 = 8 = Ncond.
        "n_per_k": {
            (0, 0, 0): 2,
            (0, -1, 0): 2,
            (1, 0, 0): 1,
            (1, -1, 0): 1,
            (-1, 0, 0): 1,
            (-1, -1, 0): 1,
        },
    },
    "case_soc_rashba_2d_sub_apbc": {
        # Shipping fixture. Same geometry (CellShape [6,4,1] /
        # SubShape [2,2,1]) as case_soc_rashba_2d_sub but with
        # BoundaryCondition = ["antiperiodic", "periodic", "periodic"]
        # (theta = (pi, 0, 0)). The twist shifts the folded k-mesh so
        # the SCF fills three k rows with 2 carriers each and two k rows
        # with 1 carrier each; the (1, -1, 0) k row lands above the
        # Fermi level and is empty (dropped from the target because
        # _occupied_per_k_soc omits zero-carrier k rows). Sum = 2*3 +
        # 1*2 = 8 = Ncond.
        "n_per_k": {
            (0, 0, 0): 2,
            (0, -1, 0): 1,
            (1, 0, 0): 2,
            (-1, 0, 0): 2,
            (-1, -1, 0): 1,
        },
    },
    "case_soc_rashba_3d_sub_apbc_xy": {
        # Shipping fixture. 3D 4x4x4 lattice, SubShape = [2, 2, 2]
        # (folded BZ [2, 2, 2], 8 k rows), BoundaryCondition =
        # ["antiperiodic", "antiperiodic", "periodic"] (theta = (pi, pi,
        # 0), xy-plane APBC). Ncond = 20 across all four siblings (see
        # README.md). Unlike the 2D fixtures, every one
        # of the 8 folded k rows is occupied here: the 4 rows with
        # k_z = 0 carry 3 carriers each and the 4 rows with k_z = -1 carry
        # 2 each. Sum = 4*3 + 4*2 = 20 = Ncond.
        "n_per_k": {
            (0, 0, 0): 3,
            (0, 0, -1): 2,
            (0, -1, 0): 3,
            (0, -1, -1): 2,
            (-1, 0, 0): 3,
            (-1, 0, -1): 2,
            (-1, -1, 0): 3,
            (-1, -1, -1): 2,
        },
    },
    "case_soc_rashba_3d_sub_apbc_xz": {
        # Shipping fixture. Same lattice as
        # case_soc_rashba_3d_sub_apbc_xy but BoundaryCondition =
        # ["antiperiodic", "periodic", "antiperiodic"] (theta = (pi, 0,
        # pi), xz-plane APBC). Ncond = 20 (same cross-fixture pin). All 8
        # folded k rows occupied: 2 rows carry 4 carriers each
        # ((0, 0, -1) and (-1, 0, 0)) and the remaining 6 rows carry 2
        # each. Sum = 2*4 + 6*2 = 20 = Ncond.
        "n_per_k": {
            (0, 0, 0): 2,
            (0, 0, -1): 4,
            (0, -1, 0): 2,
            (0, -1, -1): 2,
            (-1, 0, 0): 4,
            (-1, 0, -1): 2,
            (-1, -1, 0): 2,
            (-1, -1, -1): 2,
        },
    },
    "case_soc_rashba_3d_sub_apbc_yz": {
        # Shipping fixture. Same lattice as
        # case_soc_rashba_3d_sub_apbc_xy but BoundaryCondition =
        # ["periodic", "antiperiodic", "antiperiodic"] (theta = (0, pi,
        # pi), yz-plane APBC). Ncond = 24 -- NOT the cross-fixture
        # Ncond=20 pin used by xy/xz: at Ncond=20 this fixture violates
        # build_pair_list's SOC pairing invariant n_occ(k) ==
        # n_occ(partner(k)) (k=2 vs partner=1: 2 != 4), which blocks
        # composite manifest generation. Ncond=24 is the value
        # that satisfies both the gap>=5e-2 floor and the
        # partner-balance invariant (see README.md "Ncond selection").
        # All 8 folded k rows occupied: 4 rows carry 4 carriers each
        # ((0,0,0), (0,0,-1), (0,-1,0), (0,-1,-1)) and the remaining 4
        # rows carry 2 each. Sum = 4*4 + 4*2 = 24 = Ncond. See
        # ``docs/en/source/algorithm/uhfk_to_mvmc.rst`` for partner balance.
        "n_per_k": {
            (0, 0, 0): 4,
            (0, 0, -1): 4,
            (0, -1, 0): 4,
            (0, -1, -1): 4,
            (-1, 0, 0): 2,
            (-1, 0, -1): 2,
            (-1, -1, 0): 2,
            (-1, -1, -1): 2,
        },
    },
    "case_soc_rashba_3d_sub_apbc_xyz": {
        # Shipping fixture. Same lattice as
        # case_soc_rashba_3d_sub_apbc_xy but BoundaryCondition =
        # ["antiperiodic", "antiperiodic", "antiperiodic"] (theta = (pi,
        # pi, pi), full 3D APBC). Ncond = 12 -- NOT the cross-fixture
        # Ncond=20 pin used by xy/xz: at Ncond=20 this fixture violates
        # build_pair_list's SOC pairing invariant n_occ(k) ==
        # n_occ(partner(k)) (k=4 vs partner=3: 2 != 4), which blocks
        # composite manifest generation. Ncond=12 is the ONLY
        # value in {12,...,24} that satisfies both the gap>=5e-2
        # floor and the partner-balance invariant (see README.md "Ncond
        # selection"). All 8 folded k rows occupied: 4 rows carry 2
        # carriers each ((0,0,-1), (0,-1,-1), (-1,0,0), (-1,-1,0)) and
        # the remaining 4 rows carry 1 each. Sum = 4*2 + 4*1 = 12 =
        # Ncond. See ``docs/en/source/algorithm/uhfk_to_mvmc.rst`` for
        # partner balance.
        "n_per_k": {
            (0, 0, 0): 1,
            (0, 0, -1): 2,
            (0, -1, 0): 1,
            (0, -1, -1): 2,
            (-1, 0, 0): 2,
            (-1, 0, -1): 1,
            (-1, -1, 0): 2,
            (-1, -1, -1): 1,
        },
    },
}


def _occupied_by_spin(occupation, column_spin, wavevector_index):
    """Return dict[spin] = set of tuple(wavevector_index) for occupied
    columns.

    Only column_spin entries in {0, 1} contribute. SOC-packed rows
    (column_spin == -1) are ignored here and must be scored via
    :func:`_occupied_per_k_soc` instead.
    """
    out = {"up": [], "down": []}
    nvol, nd = occupation.shape
    for k_row in range(nvol):
        for col in range(nd):
            if occupation[k_row, col] < 0.5:
                continue
            spn = int(column_spin[col])
            if spn == 0:
                out["up"].append(tuple(int(v) for v in wavevector_index[k_row]))
            elif spn == 1:
                out["down"].append(tuple(int(v) for v in wavevector_index[k_row]))
    return out


def _occupied_per_k_soc(occupation, wavevector_index):
    """Return dict[tuple(wv)] = spin-summed int count of occupied columns.

    Used under SOC (enable_spin_orbital=true) where column_spin == -1
    packs spin into the orbital label so a per-spin split is not
    meaningful. Rounds each per-k spin-sum to nearest int (SCF at low
    T yields integer fillings inside the fixture's HOMO-LUMO gap; a
    non-integer sum is a numerical anomaly worth surfacing).
    """
    nvol = occupation.shape[0]
    out = {}
    for k_row in range(nvol):
        # sum over all packed spin-orbital columns at this k-row
        n_sum = float(occupation[k_row].sum())
        # Reject non-integer fillings: at T=0 inside a gap the SCF
        # gives integer occupancy. A 0.5 that persists means either
        # T > 0 or the fixture parameters left a degeneracy at mu, both
        # of which invalidate the target lookup.
        n_int = int(round(n_sum))
        if abs(n_sum - n_int) > 1e-6:
            raise AssertionError(
                f"SOC occupation at k_row={k_row} "
                f"wv={tuple(int(v) for v in wavevector_index[k_row])} is "
                f"non-integer ({n_sum:.6f}); fixture parameters may leave "
                f"a partial filling at the Fermi level"
            )
        if n_int != 0:
            out[tuple(int(v) for v in wavevector_index[k_row])] = n_int
    return out


def assert_case_occupation(case_dir, work_dir):
    """Compare H-wave SCF occupation in ``<work_dir>/output/`` against the
    helper's target for ``case_dir.name``. Raises KeyError on unknown
    case, AssertionError on divergence with observed vs target dump."""
    case_dir = os.fspath(case_dir)
    work_dir = os.fspath(work_dir)
    case_name = os.path.basename(os.path.normpath(case_dir))
    if case_name not in _CASE_TARGETS:
        raise KeyError(
            f"unknown case {case_name!r}; add its target to _CASE_TARGETS "
            f"in {__file__} before running the validation harness"
        )
    target = _CASE_TARGETS[case_name]

    occ_path = os.path.join(work_dir, "output", "occupation.npz")
    eig_path = os.path.join(work_dir, "output", "eigen.npz")
    if not os.path.isfile(occ_path):
        raise AssertionError(
            f"occupation.npz missing at {occ_path}; run H-wave SCF first "
            f"(case {case_name})"
        )
    if not os.path.isfile(eig_path):
        raise AssertionError(
            f"eigen.npz missing at {eig_path}; ensure input.toml requests "
            f"eigen output (case {case_name})"
        )
    occ_data = np.load(occ_path, allow_pickle=False)
    eig_data = np.load(eig_path, allow_pickle=False)
    occupation = occ_data["occupation"]
    column_spin = occ_data["column_spin"]
    wavevector_index = eig_data["wavevector_index"]

    if "n_per_k" in target:
        # SOC path: spin-summed per k. Reject if the fixture's
        # column_spin unexpectedly still contains up/down entries.
        if not np.all(np.asarray(column_spin) == -1):
            raise AssertionError(
                f"case {case_name} declares an n_per_k (SOC) target but "
                f"column_spin contains non-SOC entries "
                f"({np.unique(column_spin).tolist()}); rerun H-wave with "
                f"enable_spin_orbital = true or fix the target"
            )
        observed = _occupied_per_k_soc(occupation, wavevector_index)
        obs_sorted = dict(sorted(observed.items()))
        tgt_sorted = dict(sorted(target["n_per_k"].items()))
        if obs_sorted != tgt_sorted:
            raise AssertionError(
                f"occupation mismatch for case {case_name} (SOC):\n"
                f"  observed n_per_k = {obs_sorted}\n"
                f"  target   n_per_k = {tgt_sorted}\n"
                f"If Transfer.dat / Ncond / EPS changed, update "
                f"_CASE_TARGETS[{case_name!r}]['n_per_k'] in {__file__}."
            )
    else:
        observed = _occupied_by_spin(occupation, column_spin, wavevector_index)
        obs_sets = {k: sorted(v) for k, v in observed.items()}
        tgt_sets = {k: sorted(v) for k, v in target.items()}
        if obs_sets != tgt_sets:
            raise AssertionError(
                f"occupation mismatch for case {case_name}:\n"
                f"  observed = {obs_sets}\n"
                f"  target   = {tgt_sets}\n"
                f"If the fixture Transfer.dat was intentionally changed "
                f"(e.g. Zeeman stabilization), update "
                f"_CASE_TARGETS[{case_name!r}] in {__file__} so pytest and "
                f"run.sh stay in sync."
            )


def _cli(argv):
    if len(argv) != 2:
        print(
            "usage: python3 assert_occupation.py <case_dir> <work_dir>",
            file=sys.stderr,
        )
        return 2
    try:
        assert_case_occupation(argv[0], argv[1])
    except (KeyError, AssertionError) as e:
        print(f"OCCUPATION ASSERTION FAILED: {e}", file=sys.stderr)
        return 1
    print(f"occupation check OK for {os.path.basename(argv[0])}")
    return 0


if __name__ == "__main__":
    sys.exit(_cli(sys.argv[1:]))
