# RPA geom-`norb` Convention Unification (spin-orbital) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make RPA interpret `geom.dat`'s `norb` as the spin-orbital count (UHFk / Wannier90 convention) in spin-orbital mode, then remove the SO+sublattice+multi-orbital fail-fast guard once the folded interaction path is verified.

**Architecture:** RPA internally keeps `self.norb` = physical orbital count and `self.nd = 2·norb_phys` = SO count, exactly as today; only the *derivation* of `self.norb` from `geom.dat` changes (SO mode: `norb_phys = geom_norb // 2`). The interaction sublattice-fold stride is split: Transfer is spin-orbital-indexed (stride = SO count), two-body/Extern are physical-indexed (stride = `norb_phys`). A pre-existing, orthogonal `trans_mod`+sublattice crash is fenced off with a narrow guard.

**Tech Stack:** Python 3.x, numpy, pytest/unittest. All work in `src/hwave/solver/rpa.py`, `src/hwave/sc.py` (none here), tests under `tests/`, fixtures under `tests/rpa/input/`, docs under `docs/{en,ja}/source/`.

**Spec:** `docs/superpowers/specs/2026-06-09-rpa-geom-norb-unification-design.md` (Codex-approved, 5th pass GO).

**Run tests from repo root.** Full suite baseline: `python -m pytest tests/ -q` → currently `395 passed`.

---

## File Structure

- `src/hwave/solver/rpa.py` — all production changes:
  - new module helper `_so_physical_norb(...)`
  - `Interaction.__init__` (rpa.py:197) — SO derivation, post-fold, even-check
  - `RPA._init_param` (rpa.py:638) — derive from `self.ham_info.norb` (DRY)
  - `Interaction._reshape_interaction` (rpa.py:290) — P4 stride split
  - `Interaction._init_interaction` (rpa.py:231) — remove multi-orbital guard
  - `RPA._read_trans_mod` (rpa.py:1245) — B1 narrow guard
  - SO transfer index range validation (in `_make_ham_trans`)
- `tests/rpa/input/geom_so.dat` — new fixture, norb=2 (SO count for 1 physical orbital)
- `tests/rpa/input/geom_so_2orb.dat` — new fixture, norb=4 (SO count for 2 physical orbitals)
- `tests/test_rpa_so_multiorb.py` — guard tests rewritten, equivalence/invariance tests added
- `tests/test_rpa_geom_norb.py` — new: odd-norb fail-fast, index-range, P4 unit test, B1 guard test
- `tests/test_rpa_spin.py`, `tests/test_rpa_1orb.py` — SO configs point at new fixtures
- `docs/{en,ja}/source/...rpa...` — SO-count convention documented

---

## Task 1: Create SO-count geom fixtures

**Files:**
- Create: `tests/rpa/input/geom_so.dat`
- Create: `tests/rpa/input/geom_so_2orb.dat`

These are pure data files (one center line per spin-orbital; positions are
arbitrary for these synthetic tests, mirror the existing fixtures). They change
no behavior yet — Task 2 wires them in.

- [ ] **Step 1: Create `geom_so.dat`** (norb = 2 = SO count for 1 physical orbital)

```
  1.000000000000   0.000000000000   0.000000000000
  0.000000000000   1.000000000000   0.000000000000
  0.000000000000   0.000000000000   1.000000000000
2
    0.000000000000000e+00     0.000000000000000e+00     0.000000000000000e+00
    0.000000000000000e+00     0.000000000000000e+00     0.000000000000000e+00
```

- [ ] **Step 2: Create `geom_so_2orb.dat`** (norb = 4 = SO count for 2 physical orbitals)

```
  1.000000000000   0.000000000000   0.000000000000
  0.000000000000   1.000000000000   0.000000000000
  0.000000000000   0.000000000000   1.000000000000
4
  0.0 0.0 0.0
  0.0 0.0 0.0
  0.0 0.0 0.0
  0.0 0.0 0.0
```

- [ ] **Step 3: Verify both parse**

Run:
```bash
python -c "import hwave.qlmsio.wan90 as w; print(w.read_geom('tests/rpa/input/geom_so.dat')['norb'], w.read_geom('tests/rpa/input/geom_so_2orb.dat')['norb'])"
```
Expected: `2 4`
(If `read_geom` is not the exact accessor, fall back to confirming line 4 is `2` / `4` respectively; the loader is exercised end-to-end in Task 2.)

- [ ] **Step 4: Commit**

```bash
git add tests/rpa/input/geom_so.dat tests/rpa/input/geom_so_2orb.dat
git commit -m "test(rpa): add SO-count geom fixtures (norb=2, norb=4)"
```

---

## Task 2: Convention flip — derive physical norb from SO-count geom

This is the atomic core change. Changing the derivation alone would set
`norb_phys = 1//2 = 0` for the existing norb=1 SO fixtures, so the derivation
**and** the SO test-config migration land together; the gate is that the full
suite stays green (norb_phys=1 and unfolded multi-orbital SO are physics-
preserving: new geom norb = 2× old, cancels the `//2`).

**Files:**
- Modify: `src/hwave/solver/rpa.py` (new helper; `Interaction.__init__` :197; `RPA._init_param` :638)
- Modify: `tests/test_rpa_spin.py:38`, `tests/test_rpa_1orb.py:181`, `tests/test_rpa_so_multiorb.py` (`_run` calls + guard `_construct`)
- Test: `tests/test_rpa_geom_norb.py` (new — odd-norb fail-fast)

- [ ] **Step 1: Write the failing odd-norb fail-fast test**

Create `tests/test_rpa_geom_norb.py`:

```python
#!/usr/bin/env python3
"""geom-norb SO-count convention: derivation, validation, fold stride, guards."""

import os
import unittest

import numpy as np

import hwave.solver.rpa as solver_rpa
from hwave.solver.rpa import _so_physical_norb


class TestSoPhysicalNorb(unittest.TestCase):
    def test_non_so_passthrough(self):
        self.assertEqual(_so_physical_norb(3, False), 3)

    def test_so_halves_even(self):
        self.assertEqual(_so_physical_norb(4, True), 2)

    def test_so_odd_raises(self):
        with self.assertRaises(ValueError):
            _so_physical_norb(3, True)

    def test_so_check_norb_targets_prefold(self):
        # divide the post-fold value (8) but validate evenness on pre-fold (3)
        with self.assertRaises(ValueError):
            _so_physical_norb(8, True, check_norb=3)


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_rpa_geom_norb.py::TestSoPhysicalNorb -v`
Expected: FAIL — `ImportError: cannot import name '_so_physical_norb'`

- [ ] **Step 3: Add the helper** to `src/hwave/solver/rpa.py`

Insert immediately after the existing `validate_chi0q_index_convention` function (before `class Lattice`):

```python
def _so_physical_norb(geom_norb, enable_spin_orbital, *, check_norb=None,
                      source="geom.dat"):
    """Physical orbital count from a geometry ``norb``.

    In spin-orbital mode ``geom.dat``'s ``norb`` is the spin-orbital count
    (= 2 * physical orbitals = Wannier90 num_wann), matching UHFk, so halve it.
    Evenness is validated on ``check_norb`` (the *original*, pre-sublattice-fold
    value) so the error names the user's actual ``geom.dat`` entry, while the
    returned count is derived from ``geom_norb`` (which may be the post-fold
    value ``orig * subvol``).
    """
    if not enable_spin_orbital:
        return geom_norb
    cn = check_norb if check_norb is not None else geom_norb
    if cn % 2 != 0:
        raise ValueError(
            "spin-orbital mode requires an even Geometry norb (the spin-orbital "
            "count = 2 * physical orbitals); got {} in {}".format(cn, source))
    return geom_norb // 2
```

- [ ] **Step 4: Run to verify the helper tests pass**

Run: `python -m pytest tests/test_rpa_geom_norb.py::TestSoPhysicalNorb -v`
Expected: PASS (4 tests)

- [ ] **Step 5: Wire the derivation into `Interaction.__init__`**

In `src/hwave/solver/rpa.py`, replace line 197 (`self.norb = param_ham["Geometry"]["norb"]`) with:

```python
        # geom norb is the spin-orbital count in SO mode (UHFk/W90 convention).
        # _init_interaction may have folded the geometry (norb -> norb*subvol),
        # so read the POST-fold value here; validate evenness on the pre-fold
        # original (param_ham_orig exists only when has_sublattice).
        post_fold_norb = param_ham["Geometry"]["norb"]
        if self.lattice.has_sublattice:
            orig_norb = self.param_ham_orig["Geometry"]["norb"]
        else:
            orig_norb = post_fold_norb
        self.norb = _so_physical_norb(post_fold_norb, self.enable_spin_orbital,
                                      check_norb=orig_norb, source="Geometry")
```

- [ ] **Step 6: Wire the derivation into `RPA._init_param` (DRY — reuse Interaction's result)**

In `src/hwave/solver/rpa.py`, replace line 638 (`self.norb = self.param_ham["geometry"]["norb"]`) with:

```python
        # Stay consistent with the Interaction's physical-orbital count
        # (already SO-halved and validated); avoids re-deriving / re-checking.
        self.norb = self.ham_info.norb
```

(Leave `self.ns = 2` and `self.nd = self.norb * self.ns` unchanged on the following lines.)

- [ ] **Step 7: Migrate `test_rpa_spin.py` SO config**

In `tests/test_rpa_spin.py`, replace line 38 (`'Geometry': 'geom.dat',`) with:

```python
                    'Geometry': 'geom_so.dat' if spin_orbital else 'geom.dat',
```

- [ ] **Step 8: Migrate `test_rpa_1orb.py` SO config**

In `tests/test_rpa_1orb.py`, the `run_test(self, params, params_ref, spin_orbital=False)` helper hardcodes `'Geometry': 'geom.dat'` at line 181. Replace it with:

```python
                    'Geometry': 'geom_so.dat' if spin_orbital else 'geom.dat',
```

- [ ] **Step 9: Migrate `test_rpa_so_multiorb.py` SO runs to the norb=4 fixture**

In `tests/test_rpa_so_multiorb.py`, change every SO (`spin_orbital=True`) usage of `geom_2orb.dat` to `geom_so_2orb.dat`; keep the non-SO comparison on `geom_2orb.dat`:
- `_run("transfer_so_2orb.dat", True, "geom_2orb.dat")` → `_run("transfer_so_2orb.dat", True, "geom_so_2orb.dat")` (lines 102 and 108)
- the guard helper `_construct` (line 76) `"Geometry": "geom_2orb.dat"` → `"Geometry": "geom_so_2orb.dat"`
- leave `_run("transfer_nonso_2orb.dat", False, "geom_2orb.dat")` (line 107) unchanged.

- [ ] **Step 10: Run the SO test files**

Run: `python -m pytest tests/test_rpa_spin.py tests/test_rpa_1orb.py tests/test_rpa_so_multiorb.py -q`
Expected: PASS for `test_rpa_spin.py`, `test_rpa_1orb.py`, and the unfolded `TestRPAMultiOrbitalSO`. The two `TestRPAMultiOrbitalSOSublatticeGuard` tests still PASS (guard still present; `geom_so_2orb.dat` norb=4 → physical 2 > 1 → guard fires). The odd-norb helper tests still pass.

- [ ] **Step 11: Run the full suite (regression gate)**

Run: `python -m pytest tests/ -q`
Expected: `399 passed` (395 baseline + 4 new `TestSoPhysicalNorb` tests; no tests removed). The gate is **zero failures**. If anything in `test_rpa_*` fails, the convention flip lost or doubled a factor — debug before continuing.

> Transient guard note: with the guard still present, its condition
> `Geometry["norb"] > 1` now reads the SO-count value, so it fires for *any* SO
> + sublattice run (subvol > 1), including `norb_phys=1`. No existing RPA test
> combines SO with `SubShape != [1,1,1]` except the multi-orbital guard tests
> (which still expect rejection here), so the suite stays green. Task 5 removes
> the guard entirely.

- [ ] **Step 12: Commit**

```bash
git add src/hwave/solver/rpa.py tests/test_rpa_geom_norb.py tests/test_rpa_spin.py tests/test_rpa_1orb.py tests/test_rpa_so_multiorb.py
git commit -m "feat(rpa): interpret geom norb as spin-orbital count in SO mode

Derive physical norb_phys = geom_norb//2 (post-fold) with even-check on the
pre-fold original; RPA._init_param reuses Interaction.norb. Migrate SO test
fixtures to SO-count geom files. Physics-preserving for norb_phys=1 and
unfolded multi-orbital SO."
```

---

## Task 3: SO transfer index range validation (D7)

The `_reshape_orbit_spin` no-op argument and the SO remap both assume every SO
transfer index is in `[0, geom_norb)` (i.e. `< nd`). Add an explicit check so an
out-of-range index fails loudly instead of silently mis-folding.

**Files:**
- Modify: `src/hwave/solver/rpa.py` (`Interaction._make_ham_trans`, SO branch ~rpa.py:401)
- Test: `tests/test_rpa_geom_norb.py`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_rpa_geom_norb.py`:

```python
class TestSoTransferIndexRange(unittest.TestCase):
    def _construct_with_transfer(self, transfer_name):
        import hwave.qlmsio.read_input_k as read_input_k
        info_mode = {
            "mode": "RPA",
            "param": {"T": 2.0, "filling": 0.5,
                      "CellShape": [8, 1, 1], "SubShape": [1, 1, 1], "Nmat": 32},
            "enable_spin_orbital": True,
            "calc_scheme": "general",
        }
        info_file = {"input": {"path_to_input": "tests/rpa/input",
                               "interaction": {"path_to_input": "tests/rpa/input",
                                               "Geometry": "geom_so_2orb.dat",
                                               "Transfer": transfer_name}},
                     "output": {"path_to_output": "tests/rpa/output"}}
        os.makedirs(info_file["input"]["interaction"]["path_to_input"], exist_ok=True)
        read_io = read_input_k.QLMSkInput(info_file["input"])
        ham = read_io.get_param("ham")
        return solver_rpa.RPA(ham, {}, info_mode)

    def test_out_of_range_so_index_raises(self):
        # geom_so_2orb.dat declares nd=4 (indices 0..3); an index of 4 is illegal.
        with self.assertRaises((ValueError, IndexError)):
            self._construct_with_transfer("transfer_so_oob.dat")
```

- [ ] **Step 2: Create the out-of-range fixture** `tests/rpa/input/transfer_so_oob.dat`

```
SO 2orb with out-of-range index (5 -> 0-based 4, nd=4)
1
2
1 1
-1 0 0 5 5 1.0 0.0
```

- [ ] **Step 3: Run to verify it fails (currently no validation)**

Run: `python -m pytest tests/test_rpa_geom_norb.py::TestSoTransferIndexRange -v`
Expected: FAIL — either it does not raise, or raises a different/late error. (If it already raises `IndexError` from the array write, the test passes; in that case still add the explicit check in Step 4 for a clear message and convert the assertion to `ValueError` only.)

- [ ] **Step 4: Add the explicit range check** in `_make_ham_trans` SO branch

In `src/hwave/solver/rpa.py`, inside the `if self.enable_spin_orbital == True:` branch of `_make_ham_trans`, in the loop over `self.param_ham["Transfer"].items()` (around the `a = _so_interleaved_to_spinblock(orbvec[0])` lines), add a guard before the assignment:

```python
            for (irvec, orbvec), v in self.param_ham["Transfer"].items():
                if not (0 <= orbvec[0] < nd and 0 <= orbvec[1] < nd):
                    raise ValueError(
                        "spin-orbital Transfer index {} out of range [0,{}); "
                        "geom norb (SO count) must cover all transfer indices"
                        .format(orbvec, nd))
                a = _so_interleaved_to_spinblock(orbvec[0])
                b = _so_interleaved_to_spinblock(orbvec[1])
                tab_r[(*irvec, a, b)] = v
```

- [ ] **Step 5: Run to verify it passes**

Run: `python -m pytest tests/test_rpa_geom_norb.py::TestSoTransferIndexRange -v`
Expected: PASS (raises `ValueError`)

- [ ] **Step 6: Run full suite**

Run: `python -m pytest tests/ -q`
Expected: zero failures.

- [ ] **Step 7: Commit**

```bash
git add src/hwave/solver/rpa.py tests/test_rpa_geom_norb.py tests/rpa/input/transfer_so_oob.dat
git commit -m "feat(rpa): validate SO transfer index range in _make_ham_trans"
```

---

## Task 4: P4 — split the sublattice-fold stride (Transfer vs physical-indexed)

`_reshape_interaction` uses one `norb_orig` (pre-fold geom norb = SO count) as
the fold stride for all terms. Transfer is SO-indexed (correct), but two-body /
Extern are physical-indexed and need stride `norb_phys = SO_count // 2`, or
folded indices overrun `self.norb = norb_phys * subvol`. Unit-test the stride in
isolation (stub-based) so the guard need not be removed yet.

**Files:**
- Modify: `src/hwave/solver/rpa.py` (`Interaction._reshape_interaction`, :290)
- Test: `tests/test_rpa_geom_norb.py`

- [ ] **Step 1: Write the failing unit test** (stub bypasses the `__init__` guard)

Add to `tests/test_rpa_geom_norb.py`:

```python
class _FakeLattice:
    def __init__(self, shape, subshape):
        self.shape = shape
        self.subshape = subshape


class TestReshapeInteractionStride(unittest.TestCase):
    def _make_interaction(self, geom_norb_orig, enable_spin_orbital, subshape):
        obj = object.__new__(solver_rpa.Interaction)
        obj.enable_spin_orbital = enable_spin_orbital
        obj.lattice = _FakeLattice(shape=(2, 1, 1), subshape=subshape)
        obj.param_ham_orig = {"Geometry": {"norb": geom_norb_orig}}
        return obj

    def test_two_body_stride_is_physical_in_so_mode(self):
        # SO mode, geom norb (SO count) = 4 -> physical = 2.
        # A two-body term (enable_spin_orbital=False call) on a single physical
        # orbital index 1, folded over subshape (2,1,1), must land within
        # [0, norb_phys*subvol) = [0, 4), not stride by the SO count (-> 5+).
        obj = self._make_interaction(4, True, subshape=(2, 1, 1))
        ham = {((0, 0, 0), (1, 1)): 1.0}
        out = obj._reshape_interaction(ham, enable_spin_orbital=False)
        folded_indices = [a for (_ir, (a, b)) in out.keys() for a in (a, b)]
        self.assertTrue(all(0 <= a < 4 for a in folded_indices),
                        "physical-indexed fold must stride by norb_phys")

    def test_transfer_stride_is_so_count_in_so_mode(self):
        # Transfer (enable_spin_orbital=True call) keeps the SO-count stride.
        obj = self._make_interaction(4, True, subshape=(1, 1, 1))
        ham = {((0, 0, 0), (3, 3)): 1.0}
        out = obj._reshape_interaction(ham, enable_spin_orbital=True)
        # subshape (1,1,1) => identity fold; index 3 preserved (< SO count 4).
        self.assertIn(((0, 0, 0), (3, 3)), out)
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_rpa_geom_norb.py::TestReshapeInteractionStride -v`
Expected: `test_two_body_stride_is_physical_in_so_mode` FAILS (folded index ≥ 4, stride used SO count 4 instead of physical 2). `test_transfer_stride_is_so_count_in_so_mode` may already pass.

- [ ] **Step 3: Implement the stride split**

In `src/hwave/solver/rpa.py`, in `_reshape_interaction`, replace line 290
(`norb_orig = self.param_ham_orig["Geometry"]["norb"]`) with:

```python
        geom_norb_orig = self.param_ham_orig["Geometry"]["norb"]
        # P4: Transfer is spin-orbital-indexed (this call passes
        # enable_spin_orbital=True) -> stride by the SO count. Two-body / Extern
        # are physical-orbital-indexed (called with enable_spin_orbital=False)
        # -> in SO mode stride by the physical count = SO_count // 2.
        if self.enable_spin_orbital and not enable_spin_orbital:
            norb_orig = geom_norb_orig // 2
        else:
            norb_orig = geom_norb_orig
```

(The inner `_reshape_orbit_` / `_reshape_orbit_spin` closures already reference
`norb_orig`; no further change there.)

- [ ] **Step 4: Run to verify it passes**

Run: `python -m pytest tests/test_rpa_geom_norb.py::TestReshapeInteractionStride -v`
Expected: PASS (both)

- [ ] **Step 5: Run full suite**

Run: `python -m pytest tests/ -q`
Expected: zero failures (non-SO unchanged: `self.enable_spin_orbital` False → `norb_orig = geom_norb_orig`).

- [ ] **Step 6: Commit**

```bash
git add src/hwave/solver/rpa.py tests/test_rpa_geom_norb.py
git commit -m "fix(rpa): physical fold stride for non-Transfer interactions in SO mode (P4)"
```

---

## Task 5: Remove the SO+sublattice+multi-orbital guard (test-first)

With P4 in place, the folded multi-orbital SO path is correct. Rewrite the
guard-expectation tests to assert success, add folded/unfolded equivalence and
invariance tests, then delete the guard.

**Files:**
- Modify: `tests/test_rpa_so_multiorb.py` (rewrite `TestRPAMultiOrbitalSOSublatticeGuard`; add equivalence/invariance tests)
- Modify: `src/hwave/solver/rpa.py` (`_init_interaction` :231-237 — remove guard)

- [ ] **Step 1: Rewrite the guard tests + add equivalence/invariance (RED)**

In `tests/test_rpa_so_multiorb.py`, replace the entire
`TestRPAMultiOrbitalSOSublatticeGuard` class with success-expecting +
equivalence tests. The `_run` helper already returns `(solver, green)`; reuse it
with `SubShape` by adding a parameter. First extend `_run` to accept `subshape`:

```python
def _run(transfer, spin_orbital, geom, subshape=(1, 1, 1)):
    info_mode = {
        "mode": "RPA",
        "param": {
            "T": 2.0,
            "filling": 0.5,
            "CellShape": [8, 1, 1],
            "SubShape": list(subshape),
            "Nmat": 32,
        },
        "enable_spin_orbital": spin_orbital,
        "calc_scheme": "general",
    }
    info_file = {
        "input": {
            "path_to_input": "tests/rpa/input",
            "interaction": {
                "path_to_input": "tests/rpa/input",
                "Geometry": geom,
                "Transfer": transfer,
            },
        },
        "output": {"path_to_output": "tests/rpa/output"},
    }
    os.makedirs(info_file["output"]["path_to_output"], exist_ok=True)
    read_io = read_input_k.QLMSkInput(info_file["input"])
    ham = read_io.get_param("ham")
    solver = solver_rpa.RPA(ham, {}, info_mode)
    green = read_io.get_param("green")
    solver.solve(green, info_file["output"]["path_to_output"])
    return solver, green
```

Then replace `TestRPAMultiOrbitalSOSublatticeGuard` with:

```python
class TestRPAMultiOrbitalSOSublatticeFold(unittest.TestCase):
    """Folding (SubShape vol > 1) with multi-orbital SO is now supported."""

    def test_subcell_folding_runs(self):
        solver, _ = _run("transfer_so_2orb.dat", True, "geom_so_2orb.dat",
                         subshape=(2, 1, 1))
        self.assertEqual(solver.spin_mode, "spin-free")

    def test_folded_chi0q_matches_unfolded(self):
        # Physical observable must not depend on the (artificial) sublattice fold.
        s_unfold, g_unfold = _run("transfer_so_2orb.dat", True,
                                  "geom_so_2orb.dat", subshape=(1, 1, 1))
        s_fold, g_fold = _run("transfer_so_2orb.dat", True,
                              "geom_so_2orb.dat", subshape=(2, 1, 1))
        # chi0q layouts differ by the fold; compare a fold-invariant scalar:
        # the sum of eigen-spectrum weight is preserved. Use total |chi0q|.
        self.assertAlmostEqual(float(np.abs(g_unfold["chi0q"]).sum()),
                               float(np.abs(g_fold["chi0q"]).sum()), places=6)
```

> **Note on the invariance assertion:** if `np.abs(chi0q).sum()` proves not
> fold-invariant for this transfer (folding reindexes q-points), fall back to
> comparing the static `chi0q` at `q=0` summed over orbitals, which is fold
> invariant. Decide empirically in Step 2 — the goal is a scalar that must match.

- [ ] **Step 2: Run to verify the new tests fail (guard still present)**

Run: `python -m pytest tests/test_rpa_so_multiorb.py -v`
Expected: the two new `...Fold` tests FAIL with `NotImplementedError` (guard fires for subvol>1). If `test_folded_chi0q_matches_unfolded`'s scalar is not invariant, adjust per the note before proceeding.

- [ ] **Step 3: Remove the guard**

In `src/hwave/solver/rpa.py` `_init_interaction`, delete the guard block
(rpa.py:220-237 — the explanatory comment through the `raise NotImplementedError(msg)`):

```python
        # The spin-orbital -> spin-block remap in _make_ham_trans handles the
        # ... (entire comment block) ...
        if (self.enable_spin_orbital and self.lattice.has_sublattice
                and self.param_ham["Geometry"]["norb"] > 1):
            msg = (...)
            logger.error(msg)
            raise NotImplementedError(msg)
```

Leave the following `if self.lattice.has_sublattice:` folding block intact.

- [ ] **Step 4: Run to verify the new tests pass**

Run: `python -m pytest tests/test_rpa_so_multiorb.py -v`
Expected: PASS (all, including `...Fold` and the existing unfolded `TestRPAMultiOrbitalSO`).

- [ ] **Step 5: Run full suite**

Run: `python -m pytest tests/ -q`
Expected: zero failures.

- [ ] **Step 6: Commit**

```bash
git add src/hwave/solver/rpa.py tests/test_rpa_so_multiorb.py
git commit -m "feat(rpa): support SO + sublattice + multi-orbital; remove guard

P4 stride fix makes the folded multi-orbital SO interaction path correct, so
the bd777a6 fail-fast guard is removed. Folded/unfolded chi0q invariance and
sublattice-fold equivalence tests added."
```

---

## Task 6: B1 — narrow guard for the pre-existing `trans_mod` + sublattice crash

`RPA._read_trans_mod` calls `self.lattice._reshape_green` (rpa.py:1281), which
does not exist on `Lattice` — a pre-existing crash for any `trans_mod` +
sublattice run, orthogonal to SO. Fence it off with a narrow fail-fast so guard
removal in Task 5 leaves no reachable `AttributeError`. Do **not** delete
`RPA._reshape_green`; the real fix is a separate follow-up.

**Files:**
- Modify: `src/hwave/solver/rpa.py` (`RPA._read_trans_mod`, top, ~rpa.py:1245)
- Test: `tests/test_rpa_geom_norb.py`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_rpa_geom_norb.py`:

```python
class TestTransModSublatticeGuard(unittest.TestCase):
    def _solver(self, subshape):
        import hwave.qlmsio.read_input_k as read_input_k
        info_mode = {
            "mode": "RPA",
            "param": {"T": 2.0, "filling": 0.5,
                      "CellShape": [8, 1, 1], "SubShape": list(subshape), "Nmat": 32},
            "enable_spin_orbital": False,
            "calc_scheme": "general",
        }
        info_file = {"input": {"path_to_input": "tests/rpa/input",
                               "interaction": {"path_to_input": "tests/rpa/input",
                                               "Geometry": "geom.dat",
                                               "Transfer": "transfer.dat"}},
                     "output": {"path_to_output": "tests/rpa/output"}}
        read_io = read_input_k.QLMSkInput(info_file["input"])
        ham = read_io.get_param("ham")
        return solver_rpa.RPA(ham, {}, info_mode)

    def test_trans_mod_with_sublattice_is_rejected(self):
        solver = self._solver(subshape=(2, 1, 1))
        with self.assertRaises(NotImplementedError):
            solver._read_trans_mod("tests/rpa/input/transfer.dat")
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_rpa_geom_norb.py::TestTransModSublatticeGuard -v`
Expected: FAIL — raises `AttributeError` (the mis-wire) or reads successfully, not `NotImplementedError`.

- [ ] **Step 3: Add the narrow guard** at the top of `_read_trans_mod`

In `src/hwave/solver/rpa.py`, immediately after the `logger.debug(">>> RPA._read_trans_mod")` line in `_read_trans_mod`, add:

```python
        if self.lattice.has_sublattice:
            # Known pre-existing gap (since 2023): the sublattice reshape path
            # calls a non-existent Lattice._reshape_green. Fail fast instead of
            # raising AttributeError mid-read. Tracked as a separate follow-up.
            raise NotImplementedError(
                "trans_mod input combined with sublattice folding "
                "(SubShape volume > 1) is not yet supported in RPA.")
```

- [ ] **Step 4: Run to verify it passes**

Run: `python -m pytest tests/test_rpa_geom_norb.py::TestTransModSublatticeGuard -v`
Expected: PASS

- [ ] **Step 5: Run full suite**

Run: `python -m pytest tests/ -q`
Expected: zero failures. (No existing test exercises trans_mod + sublattice.)

- [ ] **Step 6: Commit**

```bash
git add src/hwave/solver/rpa.py tests/test_rpa_geom_norb.py
git commit -m "fix(rpa): fail-fast on trans_mod + sublattice (pre-existing gap, B1)"
```

---

## Task 7: Documentation — RPA SO uses the SO-count convention

**Files:**
- Modify: `docs/en/source/rpa/filespecification/config/index_config.rst`
- Modify: `docs/en/source/rpa/filespecification/inputfiles/interaction.rst`
- Modify: `docs/en/source/algorithm/rpa.rst`
- Modify: `docs/ja/source/rpa/filespecification/config/index_config.rst`
- Modify: `docs/ja/source/rpa/filespecification/inputfiles/interaction.rst`
- Modify: `docs/ja/source/algorithm/rpa.rst`

- [ ] **Step 1: Locate the SO / `enable_spin_orbital` / geometry-norb prose**

Run:
```bash
grep -rn "spin_orbital\|spin-orbital\|num_wann\|norb" docs/en/source/rpa docs/ja/source/rpa docs/en/source/algorithm/rpa.rst docs/ja/source/algorithm/rpa.rst
```

- [ ] **Step 2: Edit each file** to state the convention

For each RPA SO mention, add/clarify (EN wording; mirror in JA):

> In spin-orbital mode (`enable_spin_orbital = true`), the `norb` value in the
> geometry file is the **spin-orbital count** (= 2 × number of physical
> orbitals = Wannier90 `num_wann`), matching UHFk. The transfer-file orbital
> index uses the interleaved convention `2*orb + spin`.

Add a one-line migration note where SO input is described:

> **Migration (RPA):** geometry `norb` for spin-orbital input is now the
> spin-orbital count; double any pre-existing RPA SO `geom.dat` `norb`.

- [ ] **Step 3: Verify no stale "physical orbital" claim remains for RPA SO**

Run: `grep -rn "physical" docs/en/source/rpa docs/ja/source/rpa` and confirm none contradicts the SO-count convention.

- [ ] **Step 4: Commit**

```bash
git add docs/en/source docs/ja/source
git commit -m "docs(rpa): geom norb is the spin-orbital count in SO mode (+migration note)"
```

> Generated HTML under `docs/{en,ja}/build/html` and `_sources` is regenerated
> by the Sphinx build, not hand-edited here.

---

## Task 8: Final verification & spec close-out

- [ ] **Step 1: Full suite green**

Run: `python -m pytest tests/ -q`
Expected: zero failures (≈ 395 baseline + new tests).

- [ ] **Step 2: Confirm `chi0q.npz` still carries the convention marker**

Run: `python -m pytest tests/test_chi0q_index_convention.py tests/test_rpa_output.py -q`
Expected: PASS (the `index_convention="spin_block"` output is unchanged).

- [ ] **Step 3: Update the spec status & memory**

- Mark `docs/superpowers/specs/2026-06-09-rpa-geom-norb-unification-design.md`
  status → `IMPLEMENTED`.
- Record the B1 follow-up (fix `self.lattice._reshape_green` mis-wire +
  `self.norb_orig`) as a tracked next-step in `next_steps.md` / memory.

- [ ] **Step 4: Codex review of the full implementation diff**

Dispatch a Codex review (`codex:codex-rescue`, `--fresh`) over the cumulative
diff `git diff origin/feature/flex-approximation...HEAD` for `src/hwave/solver/rpa.py`
and the new/changed tests, focused on the P4 stride split, guard removal, and
the convention derivation. Address findings (verify before fixing).

- [ ] **Step 5: Commit any review fixes, then stop for user review before push.**
