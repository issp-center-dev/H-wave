# Paramagnetic Full-Vertex Multi-Orbital FLEX Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a paramagnetic (`spin-free`) full-vertex multi-orbital FLEX path to H-wave, selected by `calc_scheme="general"`, that keeps the orbital off-diagonal (Hund, pair-hopping) vertices instead of reducing them to density-density.

**Architecture:** Add a parallel "general" code path to `FLEX` (subclass of `RPA`) dispatched from `solve()` on `calc_scheme=="general"`; the existing reduced/squashed density-density path is untouched. The general path builds the full Kanamori `Ûˢ`/`Ûᶜ` matrices (a NEW MYO-convention builder — `sc.py` follows Kuroki and differs at charge `(ab,ab)`; left untouched, see spec §3), solves the matrix RPA for `χ_s`/`χ_c`, assembles the Mochizuki–Yanase–Ogata (MYO) fluctuation interaction, and computes the self-energy by a rank-4 orbital contraction. Correctness is anchored by a brute-force, physical-index, direct-summation reference that must match the optimized path in one shot.

**Tech Stack:** Python, NumPy (FFT via `numpy.fft`), pytest/unittest. Design spec: `docs/superpowers/specs/2026-06-16-full-vertex-multiorbital-flex-design.md`.

**Reference formula (MYO, [cond-mat/0407094], corroborated by THU [cond-mat/0309575]):**
```
χ⁰_{mn,μν}(q) = −(T/N) Σ_k G_{μm}(k+q) G_{nν}(k)                    (Eq.5)
χˢ(q) = [I − χ⁰(q) Ûˢ]⁻¹ χ⁰(q),   χᶜ(q) = [I + χ⁰(q) Ûᶜ]⁻¹ χ⁰(q)   (Eq.4)
V(q) = 3/2 Ûˢ χˢ Ûˢ + 1/2 Ûᶜ χᶜ Ûᶜ − 1/4 (Ûˢ+Ûᶜ) χ⁰ (Ûˢ+Ûᶜ)        (Eq.4, fluctuation part; first-order constants EXCLUDED — match reduced)
Σ_{mn}(k) = (T/N) Σ_q Σ_{μν} V_{μm,νn}(q) G_{μν}(k−q)               (Eq.3)
```

---

## File Structure

- `src/hwave/solver/flex.py` — add general-path methods + dispatch + guards (existing reduced methods untouched).
- `src/hwave/solver/_sc_matrices_myo.py` *(new)* — full `(norb², norb²)` S/C builder in the **MYO Eq.(6) convention** (charge `(ab,ab) = −U'+2J`), used by FLEX-general. NOT a reuse of `sc.py` (which is Kuroki, `−U'+J`); `sc.py` is left untouched. See §3/§4.3 of the spec.
- `tests/test_flex_general.py` *(new)* — unit + integration tests for the general path, guards, and the limit checks.
- `tests/flex_bruteforce_ref.py` *(new)* — standalone physical-index direct-summation reference for `Σ[G]` (imported by the equivalence test).
- `docs/en/source/flex/sample/iron_2orb_general/` and `docs/tutorial/Hubbard/FLEX/iron_2orb_general/` *(new)* — `calc_scheme="general"` tutorial variant.

---

## Task 0: Confirm the MYO index placement from the paper PDF — ✅ DONE (2026-06-20)

Read `cond-mat/0407094` directly from the PDF (`pdftotext`). The transcription in
this plan's reference block and spec §4.4/§4.5 **matches the paper** — no
discrepancy. Confirmed:

- [x] **Eq.(5):** `χ⁰_{mn,μν}(q) = −(T/N) Σ_k G_{μm}(k+q) G_{nν}(k)`.
- [x] **Eq.(4):** `V_{mn,μν}(q) = [ 3/2 Ûˢχ̂ˢÛˢ + 1/2 Ûᶜχ̂ᶜÛᶜ − 1/4 (Ûˢ+Ûᶜ)χ̂⁰(Ûˢ+Ûᶜ) + 3/2 Ûˢ − 1/2 Ûᶜ ]_{mn,μν}`; `χ̂ˢ=[Î−χ̂⁰Ûˢ]⁻¹χ̂⁰`, `χ̂ᶜ=[Î+χ̂⁰Ûᶜ]⁻¹χ̂⁰`.
- [x] **Eq.(3):** `Σ_{mn}(k) = (T/N) Σ_q Σ_{μν} V_{μm,νn}(q) G_{μν}(k−q)` — V matrix is in `(mn),(μν)` convention (Eq.4); Eq.3 reads element `[row=(μ,m), col=(ν,n)]`.
- [x] **Eq.(6) interaction matrices** (`Ûˢ=diag(Û₁ˢ,Û₂ˢ)`, `Ûᶜ=diag(Û₁ᶜ,Û₂ᶜ)`):
  `[U₁ˢ]_{aa,bb}`=U(a=b)/J_H(a≠b); `[U₁ᶜ]_{aa,bb}`=U/(2U'−J_H);
  `[U₂ˢ]_{ab,cd}`(a≠b,c≠d)=U' if (a,b)=(c,d), J' if (a,b)=(d,c), else 0;
  `[U₂ᶜ]_{ab,cd}`=−U'+2J_H if (a,b)=(c,d), J' if (a,b)=(d,c), else 0.

> Reminder for downstream tasks: the brute-force (Task 3) and optimized path
> (Tasks 6–7) use exactly these slots. **RESOLVED (2026-06-20):** both PDFs read —
> `sc.py` follows **Kuroki** (charge `(ab,ab) = −U'+J`); **MYO = −U'+2J** (the only
> differing element). FLEX uses **MYO** (its own builder, Task 2); `sc.py` is left
> untouched. See spec §3/§4.3.

---

## Task 1: Reject general-FLEX for unsupported modes; accept it for spin-free

**Files:**
- Modify: `src/hwave/solver/flex.py` (`_init_flex_param`, ~`flex.py:79-115`)
- Test: `tests/test_flex_general.py`

- [ ] **Step 1: Write the failing guard tests.**
```python
import unittest
import numpy as np
from requests.structures import CaseInsensitiveDict
from hwave.solver.flex import FLEX

def _ham_info(norb=2, enable_so=False):
    # Minimal ham_info stub sufficient for FLEX.__init__ up to _init_flex_param.
    # (Mirror the construction used by the existing FLEX tests; see
    #  tests/test_flex.py for the established fixture pattern.)
    ...

class TestGeneralFlexGuards(unittest.TestCase):
    def _mode(self, scheme, spin="spin-free", enable_so=False):
        return CaseInsensitiveDict({
            "mode": "FLEX",
            "calc_scheme": scheme,
            # ... other required mode params from tests/test_flex.py fixtures
        })

    def test_general_rejected_for_spin_diag(self):
        with self.assertRaises(ValueError):
            FLEX(_ham_info_spin_diag(), {}, self._mode("general", "spin-diag"))

    def test_general_rejected_for_spinful(self):
        with self.assertRaises(ValueError):
            FLEX(_ham_info_spinful(), {}, self._mode("general", "spinful"))

    def test_general_rejected_for_enable_spin_orbital(self):
        with self.assertRaises(ValueError):
            FLEX(_ham_info(enable_so=True), {}, self._mode("general", enable_so=True))

    def test_general_accepted_for_spin_free(self):
        # Should NOT raise.
        FLEX(_ham_info(), {}, self._mode("general", "spin-free"))
```
> Build `_ham_info*` fixtures by copying the established FLEX test setup in `tests/test_flex.py` (do not invent a new fixture shape).

- [ ] **Step 2: Run to verify it fails.**
Run: `PYTHONPATH=src python -m pytest tests/test_flex_general.py::TestGeneralFlexGuards -v`
Expected: FAIL (general currently raises for *all* schemes, so `test_general_accepted_for_spin_free` fails; the reject tests may pass for the wrong reason).

> **IMPORTANT (verified against the code):** `self.spin_mode` is **NOT set at
> `__init__` time** — it is assigned later in `RPA._read_chi0q` (via `read_init`,
> `rpa.py:1220+`) and in `RPA.solve` (`rpa.py:1544+`). So the `spin_mode` guard
> **cannot** live in `_init_flex_param` (it would raise `AttributeError`). Split
> the guards: the scheme acceptance + `enable_spin_orbital` rejection + the
> `_flex_general` flag go in `_init_flex_param` (available at init); the
> `spin_mode != "spin-free"` rejection goes at the **start of `solve()`** (after
> `spin_mode` is determined). The guard tests therefore split accordingly (see
> Step 1 update below).

- [ ] **Step 3a: Update `_init_flex_param`** — accept `general`, reject
`enable_spin_orbital`, set the flag (NO `spin_mode` here):
```python
scheme = self.calc_scheme.lower()
if scheme == "general":
    # Paramagnetic full-vertex FLEX (v1). spin_mode is checked in solve()
    # (not yet set here). enable_spin_orbital IS available now.
    if getattr(self.ham_info, "enable_spin_orbital", False):
        raise ValueError(
            "calc_scheme='general' FLEX (v1) does not support "
            "enable_spin_orbital; deferred to the generalized FLEX solver.")
    self._flex_general = True
elif scheme in ("reduced", "squashed"):
    self._flex_general = False
else:
    if getattr(self, "calc_type", "ring") == "ring+ladder":
        msg = "FLEX does not support calc_type='ring+ladder'."
    else:
        msg = ("FLEX requires calc_scheme='reduced', 'squashed', or "
               "'general', got '{}'.".format(self.calc_scheme))
    logger.error(msg)
    raise ValueError(msg)
```

- [ ] **Step 3b: Add the `spin_mode` guard at the start of `FLEX.solve()`** (after
`super()`/chi0q has set `self.spin_mode`; place it right after `spin_mode` is
known, before building vertices):
```python
if self._flex_general and self.spin_mode != "spin-free":
    raise ValueError(
        "calc_scheme='general' FLEX (v1) supports spin_mode='spin-free' "
        "only, got '{}'. spin-diag/spinful are deferred to the generalized "
        "FLEX solver.".format(self.spin_mode))
```
> Update Step 1's tests: `test_general_rejected_for_enable_spin_orbital` and
> `test_general_accepted_for_spin_free` (no raise) stay at construction time;
> move the `spin-diag`/`spinful` rejection tests to call `solve()` (use the
> `_make_solver` fixture with a spin-split input, or assert the guard via a small
> `solve()` invocation). If constructing a spin-diag input is heavy, assert the
> guard logic directly by setting `solver._flex_general=True; solver.spin_mode=
> "spin-diag"` and calling the extracted guard helper.

Then guard the density-density warning so it only fires on the reduced/squashed path:
```python
if not self._flex_general and (
        self.ham_info.has_interaction_exchange()
        or self.ham_info.has_interaction_pairhop()):
    logger.warning(
        "FLEX uses the density-density reduction; exchange- and "
        "pair-hopping-type interactions ... (off-diagonal vertices are "
        "dropped). Use calc_scheme='general' to keep them.")
```

- [ ] **Step 4: Run to verify it passes.**
Run: `PYTHONPATH=src python -m pytest tests/test_flex_general.py::TestGeneralFlexGuards -v`
Expected: PASS (4 tests).

- [ ] **Step 5: Verify reduced path unchanged.**
Run: `PYTHONPATH=src python -m pytest tests/test_flex.py -q`
Expected: PASS (no regressions).

- [ ] **Step 6: Commit.**
```bash
git add src/hwave/solver/flex.py tests/test_flex_general.py
git commit -m "feat(flex): accept calc_scheme=general for spin-free, guard others"
```

---

## Task 2: New MYO-convention S/C matrix builder for FLEX

**Why not reuse `sc.py`:** `sc.py:_build_sc_matrices_all_q` follows **Kuroki**
(charge `(ab,ab) = −U'+J`). FLEX is pinned to **MYO** (charge `(ab,ab) = −U'+2J`;
the only differing element — all others agree). MYO's V_eff (§4.5) is only
self-consistent with MYO's S/C, so FLEX needs its own MYO-convention builder.
`sc.py` is left UNTOUCHED. Both PDFs were read to confirm this (spec §3).

MYO Eq.(6) elements (orbitals `a`, `b`, `a≠b`; `U`=intra, `U'`=inter, `J`=Hund,
`J'`=pair-hopping):

| case (row pair, col pair) | `Ûˢ` | `Ûᶜ` |
|---|---|---|
| `(aa,aa)` | `U` | `U` |
| `(ab,ab)` (`l1=l3, l2=l4`) | `U'` | `−U'+2J` |
| `(aa,bb)` (`l1=l2, l3=l4`) | `J` | `2U'−J` |
| `(ab,ba)` (`l1=l4, l2=l3`) | `J'` | `J'` |

(Kuroki differs ONLY at `Ûᶜ_{(ab,ab)}`: `−U'+J`.)

**Files:**
- Create: `src/hwave/solver/_sc_matrices_myo.py`
- Test: `tests/test_flex_general.py`

- [ ] **Step 1: Write the failing tests** — the MYO builder produces the table
above, and intentionally diverges from `sc.py` (Kuroki) only at charge `(ab,ab)`:
```python
import numpy as np
from hwave.solver._sc_matrices_myo import build_sc_matrices_myo

def _kanamori_inter_k(norb=2, U=4.0, Up=2.0, J=0.5, Jp=0.5, Nx=2, Ny=2, Nz=1):
    """Build the inter_k dict the builder consumes (same shape sc.py expects:
    each (norb,norb,Nx,Ny,Nz)). Mirror how sc.py:_build_interaction_k /
    _build_sc_matrices_all_q read CoulombIntra/CoulombInter/Hund/Exchange/PairHop
    — inspect sc.py to match the exact keys and (norb,norb,Nx,Ny,Nz) layout."""
    ...

class TestMYOSCMatrices(unittest.TestCase):
    def test_myo_elements(self):
        U, Up, J, Jp = 4.0, 2.0, 0.5, 0.5
        S, C = build_sc_matrices_myo(_kanamori_inter_k(2, U, Up, J, Jp),
                                     norb=2, Nx=2, Ny=2, Nz=1)
        # flatten idx = l1*norb+l2 (row), l3*norb+l4 (col)
        def el(M, l1, l2, l3, l4):
            return M[0, 0, 0, l1*2+l2, l3*2+l4]
        # (ab,ab): a=0,b=1
        self.assertAlmostEqual(el(S, 0,1,0,1), Up)        # U^s = U'
        self.assertAlmostEqual(el(C, 0,1,0,1), -Up + 2*J) # U^c = -U'+2J  (MYO)
        # (aa,bb)
        self.assertAlmostEqual(el(S, 0,0,1,1), J)
        self.assertAlmostEqual(el(C, 0,0,1,1), 2*Up - J)
        # (ab,ba)
        self.assertAlmostEqual(el(S, 0,1,1,0), Jp)
        self.assertAlmostEqual(el(C, 0,1,1,0), Jp)
        # (aaaa)
        self.assertAlmostEqual(el(S, 0,0,0,0), U)
        self.assertAlmostEqual(el(C, 0,0,0,0), U)

    def test_diverges_from_kuroki_only_at_charge_abab(self):
        import hwave.sc as scmod
        ik = _kanamori_inter_k(2, 4.0, 2.0, 0.5, 0.5)
        Sm, Cm = build_sc_matrices_myo(ik, 2, 2, 2, 1)
        Sk, Ck = scmod._build_sc_matrices_all_q(ik, 2, 2, 2, 1)
        np.testing.assert_allclose(Sm, Sk)                # spin identical
        # charge differs ONLY at (ab,ab) = (0,1,0,1) and (1,0,1,0), by +J each
        diff = Cm - Ck
        # zero everywhere except the two (ab,ab) diagonal-of-offdiag entries
        nonzero = np.abs(diff) > 1e-12
        self.assertEqual(int(nonzero.sum()), 2 * 2 * 1)   # 2 entries x (Nx*Ny*Nz=4)? adjust to your grid
```
> Adjust the exact nonzero-count assertion to your chosen tiny grid; the point is: charge differs ONLY at the `(ab,ab)` entries, by `+J`.

- [ ] **Step 2: Run to verify it fails.**
Run: `PYTHONPATH=src python -m pytest tests/test_flex_general.py::TestMYOSCMatrices -v`
Expected: FAIL (`_sc_matrices_myo` does not exist).

- [ ] **Step 3: Implement `build_sc_matrices_myo`.** Use `sc.py:_build_sc_matrices_all_q` (`sc.py:444-523`) as a STRUCTURAL template (same `inter_k` reading, same index-case masks, same `(Nx,Ny,Nz,nd,nd)` output), but set the MYO values. The ONLY change vs the Kuroki code is Case 2 (`l1==l3, l2==l4, l1!=l2`) charge: MYO adds `+2*J_mat` (not `+1*J_mat`). Concretely, in the Case-2 block use `c_q += 2.0 * J_mat[_l1,_l2]` instead of Kuroki's `c_q += J_mat[...]`. Keep the Ising handling identical to sc.py (H-wave separates Hund/Ising; MYO has no Ising so it only affects models that set Ising). Document the one-line difference with a comment citing MYO Eq.(6).

- [ ] **Step 4: Run to verify it passes.**
Run: `PYTHONPATH=src python -m pytest tests/test_flex_general.py::TestMYOSCMatrices -v`
Expected: PASS.

- [ ] **Step 5: Confirm `sc.py` untouched.**
Run: `PYTHONPATH=src python -m pytest tests/test_sc.py -q` and `git status` (only the new file + test changed).
Expected: PASS; `sc.py` unmodified.

- [ ] **Step 6: Commit.**
```bash
git add src/hwave/solver/_sc_matrices_myo.py tests/test_flex_general.py
git commit -m "feat(flex): MYO-convention S/C matrix builder (charge ab,ab = -U'+2J)"
```

---

## Task 3: Brute-force physical-index reference for Σ[G]

Write this BEFORE the optimized path so it is the independent ground truth (no flatten convention shared with the optimized code).

**Files:**
- Create: `tests/flex_bruteforce_ref.py`
- Test: `tests/test_flex_general.py`

- [ ] **Step 1: Implement the reference** straight from MYO Eqs.(3)–(5) with explicit orbital loops (no FFT, no flatten):
```python
import numpy as np

def chi0_bruteforce(G, T, Nk):
    # G: (norb, norb, Nk, nmat) on a 1D k-grid for simplicity (extendable).
    # χ⁰_{mn,μν}(q,iν) = −(T/Nk) Σ_{k,iω} G_{μm}(k+q,iω+iν) G_{nν}(k,iω)
    norb = G.shape[0]; nmat = G.shape[-1]
    chi0 = np.zeros((norb, norb, norb, norb, Nk, nmat), dtype=complex)
    for m in range(norb):
     for n in range(norb):
      for mu in range(norb):
       for nu in range(norb):
        for q in range(Nk):
         for iv in range(nmat):
          s = 0j
          for k in range(Nk):
           for iw in range(nmat):
            kq = (k + q) % Nk; iwv = (iw + iv) % nmat
            s += G[mu, m, kq, iwv] * G[n, nu, k, iw]
          chi0[m, n, mu, nu, q, iv] = -(T / Nk) * s
    return chi0

def sigma_bruteforce(G, V, T, Nk):
    # Σ_{mn}(k,iω) = (T/Nk) Σ_{q,iν} Σ_{μν} V_{μm,νn}(q,iν) G_{μν}(k−q, iω−iν)
    norb = G.shape[0]; nmat = G.shape[-1]
    Sig = np.zeros((norb, norb, Nk, nmat), dtype=complex)
    for m in range(norb):
     for n in range(norb):
      for k in range(Nk):
       for iw in range(nmat):
        s = 0j
        for q in range(Nk):
         for iv in range(nmat):
          kq = (k - q) % Nk; iwv = (iw - iv) % nmat
          for mu in range(norb):
           for nu in range(norb):
            s += V[mu, m, nu, n, q, iv] * G[mu, nu, kq, iwv]
          Sig[m, n, k, iw] = (T / Nk) * s
    return Sig
```
> Use the index slots CONFIRMED in Task 0. Keep it tiny (norb≤2, Nk≤4, nmat≤4) — correctness, not speed.

- [ ] **Step 2: Add a self-test of the reference** (sanity: χ0 Hermiticity / shapes) and run it.
Run: `PYTHONPATH=src python -m pytest tests/test_flex_general.py::TestBruteForceRef -v`
Expected: PASS.

- [ ] **Step 3: Commit.**
```bash
git add tests/flex_bruteforce_ref.py tests/test_flex_general.py
git commit -m "test(flex): physical-index brute-force reference for Sigma[G]"
```

---

## Task 4: `_inflate_chi0q_and_ham_general` (rank-4 χ⁰ + full interaction)

**Files:**
- Modify: `src/hwave/solver/flex.py`
- Test: `tests/test_flex_general.py`

- [ ] **Step 1: Write the failing test** — for `spin-free`, the general inflation returns rank-4 `χ⁰` of shape `(nfreq, nvol, nd, nd, nd, nd)` and reuses the RPA general-scheme `chi0q` unchanged:
```python
class TestInflateGeneral(unittest.TestCase):
    def test_shape_spin_free(self):
        flex = _make_general_flex(norb=2)             # helper: spin-free general FLEX
        chi0_raw = _fake_general_chi0q(flex)          # RPA general-scheme chi0q
        chi0, Us, Uc = flex._inflate_chi0q_and_ham_general(chi0_raw, flex.ham_info_array)
        nd = flex.nd
        self.assertEqual(chi0.shape[-4:], (nd, nd, nd, nd))
```

- [ ] **Step 2: Run to verify it fails.**
Run: `PYTHONPATH=src python -m pytest tests/test_flex_general.py::TestInflateGeneral -v`
Expected: FAIL (method missing).

- [ ] **Step 3: Implement `_inflate_chi0q_and_ham_general`.** For `spin-free`, the RPA general scheme already produces rank-4 `chi0q`; return it (spin-inflated as needed to `nd`) together with the full `Ûˢ`, `Ûᶜ` from Task 2's helper:
```python
def _inflate_chi0q_and_ham_general(self, chi0q_raw, ham_orig):
    """Paramagnetic general path: rank-4 chi0q + full Kanamori S/C matrices."""
    from hwave.solver._sc_matrices_myo import build_sc_matrices_myo
    nx, ny, nz = self.lattice.shape
    # chi0q_raw is the RPA general-scheme tensor for spin-free.
    chi0q = chi0q_raw  # already (nfreq, nvol, nd, nd, nd, nd) for general/spin-free
    Us, Uc = build_sc_matrices_myo(self.inter_k, self.norb, nx, ny, nz)
    return chi0q, Us, Uc
```
> Confirm the exact RPA general-scheme `chi0q` shape/axis order against `rpa.py` `_calc_chi0q` for `spin-free` and adapt the comment/return to match; add an `assert chi0q.ndim == 6`.

- [ ] **Step 4: Run to verify it passes.**
Run: `PYTHONPATH=src python -m pytest tests/test_flex_general.py::TestInflateGeneral -v`
Expected: PASS.

- [ ] **Step 5: Commit.**
```bash
git add src/hwave/solver/flex.py tests/test_flex_general.py
git commit -m "feat(flex): rank-4 chi0q + full S/C inflation for general path"
```

---

## Task 5: `χ_s`, `χ_c` via matrix RPA + RPA-consistency check

**Files:**
- Modify: `src/hwave/solver/flex.py`
- Test: `tests/test_flex_general.py`

- [ ] **Step 1: Write the failing test** — first-iteration `χ_s`, `χ_c` equal the RPA general-scheme susceptibilities for the same `χ⁰`/interaction:
```python
class TestChiGeneralConsistency(unittest.TestCase):
    def test_matches_rpa_general(self):
        flex = _make_general_flex(norb=2)
        chi0, Us, Uc = flex._inflate_chi0q_and_ham_general(_fake_general_chi0q(flex), None)
        chis, chic = flex._solve_channels_general(chi0, Us, Uc)
        chis_ref = _rpa_general_chi(chi0, Us, sign=-1)   # [I - chi0 Us]^-1 chi0
        chic_ref = _rpa_general_chi(chi0, Uc, sign=+1)   # [I + chi0 Uc]^-1 chi0
        np.testing.assert_allclose(chis, chis_ref, atol=1e-10)
        np.testing.assert_allclose(chic, chic_ref, atol=1e-10)
```

- [ ] **Step 2: Run to verify it fails.** Expected: FAIL (`_solve_channels_general` missing).

- [ ] **Step 3: Implement `_solve_channels_general`** by calling the inherited `_solve_rpa` with `ham_s = −Ûˢ`, `ham_c = +Ûᶜ` (MYO/THU convention; same sign mapping as the reduced `_build_spin_charge_vertices`), flattening rank-4 to `(ndx, ndx)` as `_solve_rpa` expects:
```python
def _solve_channels_general(self, chi0q, Us, Uc):
    chi_s = self._solve_rpa(chi0q, -Us)   # [1 - chi0 Us]^-1 chi0
    chi_c = self._solve_rpa(chi0q, +Uc)   # [1 + chi0 Uc]^-1 chi0
    return chi_s, chi_c
```
> Verify `_solve_rpa` accepts the general rank-4 `chi0q` and the `(nvol, nd², nd²)` vertex layout; adapt the reshape if needed (see `rpa.py:2068-2138`).

- [ ] **Step 4: Run to verify it passes.** Expected: PASS.

- [ ] **Step 5: Commit.**
```bash
git add src/hwave/solver/flex.py tests/test_flex_general.py
git commit -m "feat(flex): general-path chi_s/chi_c via matrix RPA"
```

---

## Task 6: `_calc_veff_general` (MYO fluctuation interaction)

**Files:**
- Modify: `src/hwave/solver/flex.py`
- Test: `tests/test_flex_general.py`

- [ ] **Step 1: Write the failing test** — `V` matches a direct construction from the MYO formula on a small input, and reduces correctly for a single orbital (`Ûˢ=Ûᶜ=U` ⇒ matches the reduced `1.5χ_s+0.5χ_c−χ0` sandwiched by `U`):
```python
class TestVeffGeneral(unittest.TestCase):
    def test_single_orbital_matches_reduced_kernel(self):
        flex = _make_general_flex(norb=1, U=3.0)
        chi0, Us, Uc = ...   # norb=1: Us=Uc=[[U]]
        chis, chic = flex._solve_channels_general(chi0, Us, Uc)
        V = flex._calc_veff_general(chi0, chis, chic, Us, Uc)
        Vref = _U_sandwich(3.0, 1.5*chis + 0.5*chic - chi0)   # reduced fluctuation kernel
        np.testing.assert_allclose(V.reshape(Vref.shape), Vref, atol=1e-10)
```

- [ ] **Step 2: Run to verify it fails.** Expected: FAIL (method missing).

- [ ] **Step 3: Implement `_calc_veff_general`** (fluctuation part only; first-order constants excluded per spec §4.5), batched matmul over the `(nd², nd²)` flatten:
```python
def _calc_veff_general(self, chi0q, chi_s, chi_c, Us, Uc):
    # All in the (nvol, nd^2, nd^2) flatten used by _solve_rpa.
    UsS = Us[None, ...]                      # broadcast over frequency
    UcC = Uc[None, ...]
    UspUc = Us + Uc
    term_s = 1.5 * (UsS @ chi_s @ UsS)
    term_c = 0.5 * (UcC @ chi_c @ UcC)
    term_0 = 0.25 * (UspUc[None] @ chi0q @ UspUc[None])
    V = term_s + term_c - term_0
    return V
```
> Match array shapes to whatever `_solve_channels_general` returns; insert the frequency/volume broadcasting accordingly. NO `+(3/2)Ûˢ−(1/2)Ûᶜ` constants (spec §4.5).

- [ ] **Step 4: Run to verify it passes.** Expected: PASS.

- [ ] **Step 5: Commit.**
```bash
git add src/hwave/solver/flex.py tests/test_flex_general.py
git commit -m "feat(flex): MYO fluctuation V_eff for general path"
```

---

## Task 7: `_calc_self_energy_general` (rank-4 contraction)

**Files:**
- Modify: `src/hwave/solver/flex.py`
- Test: `tests/test_flex_general.py`

- [ ] **Step 1: Write the failing test** — the optimized self-energy equals the brute-force reference (Task 3) for the SAME small `G`, `V`, to ~1e-10. This is the PRIMARY correctness check and it locks the index wiring:
```python
from tests.flex_bruteforce_ref import sigma_bruteforce, chi0_bruteforce

class TestSelfEnergyGeneral(unittest.TestCase):
    def test_matches_bruteforce(self):
        flex, G, V = _tiny_general_setup(norb=2, Nk=4, nmat=4)
        sig_fast = flex._calc_self_energy_general(G, V, beta=flex.beta)
        sig_ref = sigma_bruteforce(G, V, T=flex.T, Nk=4)
        np.testing.assert_allclose(sig_fast, sig_ref, atol=1e-10)
```

- [ ] **Step 2: Run to verify it fails.** Expected: FAIL (method missing).

- [ ] **Step 3: Implement `_calc_self_energy_general`** by reusing the FFT transport of `_calc_self_energy` but replacing the per-`(r,τ)` Hadamard product with the rank-4 contraction `Σ_{mn} = Σ_{μν} V_{μm,νn} G_{μν}` (wiring frozen in spec §4.4; the test above is the safety net):
```python
def _calc_self_energy_general(self, green_kw, v_eff, beta):
    # 1. FFT G(k,iω) -> G(r,τ); FFT V(q,iν) -> V(r,τ)  [reuse _calc_self_energy transport]
    # 2. per (r,τ): Sigma_{mn} = sum_{μν} V_{μm,νn} G_{μν}
    #    implement as a batched einsum/matmul with the Task-0-confirmed slots, e.g.
    #    sigma_rt = np.einsum('rmuν...?', v_rt, g_rt)  # EXACT subscripts per confirmed wiring
    # 3. FFT back Sigma(r,τ) -> Sigma(k,iω)
    ...
```
> Write the contraction to match `sigma_bruteforce` exactly. Start from a literal `np.einsum` mirroring the brute-force subscripts, then optimize to batched `matmul` only after the test passes (keep both behind the same test).

- [ ] **Step 4: Run to verify it passes.** Expected: PASS (~1e-10).

- [ ] **Step 5: Commit.**
```bash
git add src/hwave/solver/flex.py tests/test_flex_general.py
git commit -m "feat(flex): rank-4 self-energy contraction (matches brute-force)"
```

---

## Task 8: Wire the general path into `solve()`

**Files:**
- Modify: `src/hwave/solver/flex.py` (`solve()` SCF loop, ~`flex.py:129-244`; `_flex_compute_veff`)
- Test: `tests/test_flex_general.py`

- [ ] **Step 1: Write the failing end-to-end test** — a tiny 2-orbital general FLEX run completes and writes `sigma`/`green`/`chiq`:
```python
class TestGeneralSolveEndToEnd(unittest.TestCase):
    def test_runs_and_saves(self):
        # tiny system; assert convergence flag set and output arrays present
        info = _tiny_general_run_info(norb=2)   # input dict, calc_scheme="general"
        ... run FLEX.solve ...
        self.assertIn("sigma", green_info)
```

- [ ] **Step 2: Run to verify it fails.** Expected: FAIL (dispatch not wired).

- [ ] **Step 3: Dispatch on `self._flex_general`** inside the per-iteration vertex/self-energy computation (analogous to `_flex_compute_veff` → `_calc_self_energy`):
```python
if self._flex_general:
    chi0q, Us, Uc = self._inflate_chi0q_and_ham_general(chi0q_raw, ham_orig)
    chi_s, chi_c = self._solve_channels_general(chi0q, Us, Uc)
    v_eff = self._calc_veff_general(chi0q, chi_s, chi_c, Us, Uc)
    sigma = self._calc_self_energy_general(green_kw, v_eff, beta)
else:
    # existing reduced path (unchanged)
    ...
```

- [ ] **Step 4: Run to verify it passes.** Expected: PASS.

- [ ] **Step 5: Commit.**
```bash
git add src/hwave/solver/flex.py tests/test_flex_general.py
git commit -m "feat(flex): dispatch general path in solve()"
```

---

## Task 9: Limit & regression tests

**Files:**
- Test: `tests/test_flex_general.py`

- [ ] **Step 1: Single-orbital == reduced (exact).**
```python
def test_single_orbital_equals_reduced(self):
    e_gen = _run_flex_energy(norb=1, scheme="general")
    e_red = _run_flex_energy(norb=1, scheme="reduced")
    self.assertAlmostEqual(e_gen, e_red, places=8)
```

- [ ] **Step 2: Multi-orbital U,U' — DIAGNOSTIC difference (not equality).**
```python
def test_multiorbital_density_density_differs(self):
    e_gen = _run_flex_energy(norb=2, scheme="general", U=4, Up=2, J=0)
    e_red = _run_flex_energy(norb=2, scheme="reduced", U=4, Up=2, J=0)
    self.assertGreater(abs(e_gen - e_red), 1e-8)   # they differ (U^s != U^c)
    self.assertLess(abs(e_gen - e_red), 10.0)      # but bounded/sane
```

- [ ] **Step 3: Hund new-physics + no warning.**
```python
def test_hund_no_density_density_warning(self):
    import logging
    with self.assertLogs("qlms.solver.flex", level="WARNING") as cm:
        logging.getLogger("qlms.solver.flex").warning("sentinel")
        _run_flex_energy(norb=2, scheme="general", U=4, Up=2, J=0.5)
    self.assertFalse(any("density-density" in m for m in cm.output))
```

- [ ] **Step 4: Reduced-path regression.**
Run: `PYTHONPATH=src python -m pytest tests/test_flex.py -q`
Expected: PASS (unchanged).

- [ ] **Step 5: Run the whole new suite.**
Run: `PYTHONPATH=src python -m pytest tests/test_flex_general.py -q`
Expected: PASS.

- [ ] **Step 6: Commit.**
```bash
git add tests/test_flex_general.py
git commit -m "test(flex): general-path limit and regression tests"
```

---

## Task 10: Tutorial variant (`iron_2orb`, general scheme)

**Files:**
- Create: `docs/en/source/flex/sample/iron_2orb_general/input.toml` (+ symlink/copy of the iron_2orb `.dat` files)
- Create: `docs/tutorial/Hubbard/FLEX/iron_2orb_general/` (input + interaction files; no bulky `output/`)
- Modify: `README.md` Examples table (add the general FLEX row)

- [ ] **Step 1: Create the input.** Copy `docs/en/source/flex/sample/iron_2orb/input.toml`, set `calc_scheme = "general"`, keep the same `.dat` interaction files.

- [ ] **Step 2: Run it to confirm it works.**
Run: `cd docs/en/source/flex/sample/iron_2orb_general && PYTHONPATH=<repo>/src python -m hwave.qlms input.toml`
Expected: converges; writes `output/` (then remove the bulky `output/` before committing, matching the existing lightweight tutorials).

- [ ] **Step 3: Copy the runnable inputs into the tutorial tree** (inputs + result figures only; exclude `output/`).

- [ ] **Step 4: Add the README Examples row** for `FLEX` + `calc_scheme="general"`.

- [ ] **Step 5: Commit.**
```bash
git add docs/en/source/flex/sample/iron_2orb_general docs/tutorial/Hubbard/FLEX/iron_2orb_general README.md
git commit -m "docs(flex): add general-scheme iron_2orb tutorial sample"
```

---

## Task 11: Manual update + full test sweep

**Files:**
- Modify: `docs/{en,ja}/source/flex/tutorial/tu-index.rst`

- [ ] **Step 1: Document `calc_scheme="general"`** in the FLEX manual: what it does (keeps full Kanamori vertices, MYO formula), when to use it, the paramagnetic-only limitation, and that the density-density warning is suppressed under `general`. Reference the MYO/THU papers.

- [ ] **Step 2: Full test sweep.**
Run: `PYTHONPATH=src python -m pytest tests/test_flex.py tests/test_flex_general.py tests/test_sc.py tests/test_rpa_1orb.py tests/test_rpa_2orb.py -q`
Expected: PASS (pre-existing `np.trapz` failures in `test_dos.py`/`test_rpa_reshape_green_2d.py` are unrelated).

- [ ] **Step 3: Commit.**
```bash
git add docs/en/source/flex/tutorial/tu-index.rst docs/ja/source/flex/tutorial/tu-index.rst
git commit -m "docs(flex): document calc_scheme=general (paramagnetic full-vertex)"
```

---

## Self-Review (planner)

- **Spec coverage:** §4.1 guards → Task 1; §4.3 shared S/C → Task 2; §4.4 self-energy wiring + brute-force → Tasks 3,7; §4.5 V_eff (constants excluded) → Task 6; §5 data flow → Task 8; §6 primary 1-shot + secondary limits → Tasks 3,7,9; §7 tests → Tasks 1,9; manual/tutorial → Tasks 10,11; MYO PDF confirm → Task 0. SOC/spinful guards → Task 1. Performance budget (§9) is deferred (optimize-after-correct, Task 7 step 3) — acceptable for v1, no separate task.
- **Placeholder scan:** test fixtures reference `tests/test_flex.py` patterns deliberately (the executor copies the established fixture) rather than inventing shapes; the only intentional "fill from Task 0" points are index slots gated on the PDF confirmation, which is the correct sequencing.
- **Type consistency:** method names used consistently — `_flex_general` (flag), `_inflate_chi0q_and_ham_general`, `_solve_channels_general`, `_calc_veff_general`, `_calc_self_energy_general`, `build_sc_matrices_myo` (new MYO-convention builder). Tasks 3/7 share `sigma_bruteforce`/`chi0_bruteforce` signatures.
