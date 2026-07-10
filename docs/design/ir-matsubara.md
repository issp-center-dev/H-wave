# Design: IR-basis Matsubara axis for H-wave (FLEX + dynamic Eliashberg)

Status: DRAFT v3 for Phase-0 design review (v1 Codex + v2 Antigravity must_fix addressed).

## 1. Problem and goals

H-wave's FLEX and dynamic-Eliashberg solvers represent the Matsubara axis on a
uniform grid of `Nmat` frequencies and move between iw_n and tau with
phase-twisted FFTs (`solver/matsubara.py`). `Nmat` must grow like
O(beta*omega_max), so low temperatures force Nmat = 2048+ and the dynamic
Eliashberg tensors (norb^4 x Nk x Nmat) dominate memory (GB-scale) and time.

The intermediate representation (IR) of Green's functions (Shinaoka et al.,
PRB 96, 035147 (2017)) compresses the frequency axis to L = O(log(Lambda))
coefficients (Lambda = beta*omega_max), with sparse sampling (Li et al., PRB
101, 035144 (2020)) providing ~L well-conditioned sampling nodes in both tau
and iw_n. FLEX on sparse sampling is established practice (Witt et al., PRB
103, 205148 (2021)).

Goals:
- G1. Reduce the frequency-axis size from Nmat to L (~50-100 at production
  Lambda ~ 1e4, eps ~ 1e-8): ~20-40x memory/compute reduction on that axis.
- G2. Strictly opt-in: default behavior (uniform grid) bit-unchanged.
- G3. Compose with the existing backend abstraction (numpy/cupy), Anderson
  mixing, and sigma_init warm start.
- G4. Controlled, testable accuracy with diagnostics that can detect
  out-of-basis (product-generated) spectral weight, not only node residuals.

Non-goals (v1): RPA solver and UHF stay uniform-grid; the k-axis treatment
(uniform mesh + spatial FFT, including the single 1/N fold from ifftn) is
unchanged; no change to default output file formats.

## 2. Staging

- Stage 1: dynamic Eliashberg on IR. FLEX stays uniform-grid; the dynamic
  loader compresses the full-grid FLEX outputs to IR **before** any large
  tensor is built (see 3.2), and the Eliashberg kernel runs entirely on
  sampling nodes. Acceptance gate: matvec outputs compared against the
  current full-grid operator on small fixtures BEFORE any eigensolver runs.
- Stage 2: FLEX SCF on IR ([mode.param] switch).
- Stage 3 (follow-up, out of scope here): IR-native npz outputs after the
  downstream readers (hwave_sc loaders, analysis scripts) are audited.

## 3. Architecture

### 3.1 `solver/ir_axis.py` (wraps the optional `sparse-ir` package)

`IRAxis(beta, wmax, eps, statistics)` holds a FiniteTempBasis plus sampling
objects, with the following **explicit conventions**:

- Frequencies are identified by INTEGER Matsubara indices n (fermionic
  iw_n = (2n+1) pi T, bosonic i nu_m = 2m pi T), never by float frequency
  matching.
- Transforms are pure basis changes with NO temperature factors inside; the
  1/beta factors stay at the product exits, exactly as in the current
  `matsubara.py` contract.
- All transform applications contract the LAST axis (frequency), matching
  the existing (..., nmat)/(nblock, nmat, ...) layouts via moveaxis at the
  call site; matrices are built on the host once and mirrored to the device
  on demand (xp-dispatch as in backend.py).
- API: `freq_nodes` (integer n, guaranteed symmetric: see 3.4),
  `tau_nodes` (guaranteed symmetric under tau -> beta - tau: see 3.4),
  `fit_from_freq` / `eval_to_freq` / `fit_from_tau` / `eval_to_tau`
  (fits are precomputed pseudo-inverses of the sampling matrices),
  `eval_to_uniform(nmat)` (densification for I/O),
  `eval_at_beta_minus()` (row vector u_l(beta^-), for particle numbers),
  and `check_nodes` (oversampled diagnostic nodes, see Sec. 5).
- Large-batch applications are CHUNKED along the k axis (memory-aware, same
  pattern as the existing chi0q chunking) so transient copies never exceed
  a configurable fraction of the array size.

### 3.2 Stage-1 data flow (dynamic Eliashberg)

1. Load full-grid npz **in k-chunks**; each chunk's uniform frequency axis
   is fitted to IR by least squares over ALL stored uniform frequencies
   (design matrix U_l(iw_n_uniform); one QR shared by all chunks). The
   uniform grid thus serves as an oversampled residual check at load time:
   the max residual over the FULL grid is logged per object. Full uniform
   tensors of G2/vertex are NEVER built: `green` is compressed first
   (norb^2 x Nk x L_F), then G2 is formed directly on fermionic nodes;
   chi_s/chi_c are compressed chunk-wise, then the vertex is assembled on
   bosonic nodes (norb^4 x Nk x L_B).
2. Kernel application per matvec (tau-product, mirroring the current
   kernel's algebra): F = G2 * phi at fermionic nodes -> fit_from_freq ->
   eval_to_tau -> spatial ifftn -> multiply with V(r, tau nodes)
   (precomputed once, the IR analogue of the hoisted Vs_rt) -> spatial
   fftn -> fit_from_tau -> eval_to_freq. Temperature factors keep their
   current placement (inside G2, none in the kernel).
3. Eigenproblem size drops from norb^2*Nk*Nmat to norb^2*Nk*L_F.
4. Gap output: densified to the run's uniform grid, reproducing the EXACT
   current metadata semantics (centered grid, `iomega` array, the
   `nmat//2` smallest-positive-frequency slice in gap.dat); loader
   round-trip is a test. Densification (here and in Stage 2 outputs) is
   ALWAYS executed in memory-aware k-chunks -- the full uniform-grid tensor
   exists only inside the npz writer's chunk loop, never as a whole in
   memory -- so a large output Nmat cannot OOM an otherwise-small IR run.
   Densified fermionic objects carry the exact IR tail (the basis
   represents the 1/(iw) asymptotics within eps), so downstream Pade /
   MaxEnt consumers see smooth eps-accurate tails; the npz gains a
   `matsubara_basis="ir"` + (eps, wmax, L) metadata record documenting the
   provenance.

### 3.3 Stage-2 data flow (FLEX)

- G(k) on fermionic nodes; chi0(r, tau nodes) via the tau-product rules of
  3.4; bosonic fit -> chi0(q, i nu nodes); RPA solve batched over L_B x Nk;
  V_eff on bosonic nodes; Sigma = V*G product at tau nodes; fermionic fit
  -> Sigma at freq nodes.
- `coeff_tail` is NOT used on the IR path: the fermionic basis represents
  the 1/(iw) tail within eps by construction, so the tail-subtraction /
  green0_tail machinery is bypassed (documented contract; the uniform path
  keeps it unchanged).
- Convergence measure and Anderson mixing operate on the Sigma freq-node
  vector; additionally the COEFFICIENT-space change |d sigma_l| is logged
  each iteration as a diagnostic (a node-converged but coefficient-drifting
  state indicates an inadequate wmax; see Sec. 5).
- sigma_init: stores IR coefficients + metadata (beta, wmax, eps,
  statistics, L, cell_shape, index conventions). Loading onto a different
  basis (e.g. the next temperature) is explicitly an INTERPOLATION:
  evaluate the stored coefficients on a dense intermediate grid (the more
  stable route), fit to the new basis, and validate the residual on the new
  check nodes against a threshold (default 100*ir_tol, configurable).
  Behavior above threshold is configurable via
  `sigma_init_on_error = "warn" (default) | "abort" | "zero"`:
  automated sweeps keep running on "warn" (seed used, residual logged) or
  degrade to the zero start on "zero"; "abort" fail-fasts for interactive
  use. MIGRATION: a uniform-grid sigma.npz (no IR metadata) is accepted by
  fitting its full grid to the current basis with the same residual policy,
  so uniform -> IR transitions need no file conversion.

### 3.4 Particle number and dressed mu search (Stage 2) -- REVISED

Sparse nodes are interpolation nodes, NOT quadrature nodes, so the current
tail-subtracted Matsubara sums cannot be reused. The IR formulation is:

- For a trial mu, G at the fermionic nodes has closed form through the
  already-available eigenvalues lam of M (per-node, per-k traces; the
  closed-form `_eigvals_small` machinery is reused as a per-node evaluator
  only -- it does NOT define a sum).
- The k-summed trace g(iw_node) is a length-L_F vector; fit_from_freq gives
  g_l and the particle number is n(mu) = -sum_l g_l u_l(beta^-) via
  `eval_at_beta_minus` (the standard IR route to equal-time observables;
  sign/normalization pinned against the uniform-grid result by an
  equivalence test at Sigma = 0 AND at a random complex Sigma).
- Root finding: the SAME safeguarded bracket + expansion logic as the
  uniform path, but with derivative-free safeguarded iteration (Brent-style)
  initially; the analytic-derivative Newton acceleration is deferred until
  the IR n(mu) is proven equivalent (per-trial cost is an L_F-length fit --
  negligible either way).

### 3.5 Tau-product conventions on sparse nodes -- REVISED

There is no index-reversal trick on sparse nodes; the rules are defined in
continuous variables and realized as exact node permutations:

- Node symmetry is REQUIRED and enforced at construction: the fermionic
  frequency node set must be integer-symmetric under n -> -n-1 and the tau
  node set symmetric under tau -> beta - tau. If the library returns an
  asymmetric set, it is symmetrized (union with its mirror) and the fit
  pseudo-inverses are rebuilt; symmetry is asserted in INTEGER arithmetic
  (frequencies); tau nodes are paired with their mirrors within an
  absolute tolerance of 1e-14*beta and replaced by the exactly
  symmetrized pair means, after which the fit pseudo-inverses are rebuilt
  from the snapped nodes (so conditioning reflects the actual nodes used).
- Fermionic reflection: G(-tau) = -G(beta - tau) (anti-periodicity) -- a
  node permutation with a global -1 sign. Bosonic reflection:
  chi(-tau) = +chi(beta - tau) -- the same permutation, sign +1.
- r -> -r stays the existing reverse+roll index map on the uniform k grid
  (unchanged).
- chi0(r, tau) = -G(r, tau) * G(-r, -tau) is then well-defined at every tau
  node; no endpoint (tau = 0, beta) evaluation occurs inside products
  (nodes are interior); beta^- evaluation happens only through basis
  functions (3.4).
- The current centered-grid reversal/roll/sign-vector code remains the
  uniform path's implementation; equivalence of the two implementations is
  pinned by tests on both chi0 and the Eliashberg kernel matvec.

## 4. Configuration surface

- `[mode.param] matsubara_basis = "uniform" (default) | "ir"` (Stage 2) and
  `[eliashberg] matsubara_basis` (Stage 1).
- `ir_tol` (default 1e-8): basis cutoff eps.
- `ir_wmax` (default: auto = 3 * (single-particle band range + max
  interaction scale), ALWAYS logged, always overridable). The auto rule is
  a heuristic starting point only; its adequacy for product-generated
  spectra is enforced by the runtime diagnostics and the wmax-sweep
  validation below, not assumed. If the auto estimate cannot be formed
  (missing/degenerate inputs, non-finite band range) or falls below a
  sanity floor, the run aborts with a message stating the computed pieces
  and asking for an explicit `ir_wmax` -- never a silent default.
- `ir_check_interval` (default 10): SCF period of the two-basis check when
  `ir_check = true`.
- `matsubara_basis` consistency: when a FLEX run and its consuming
  Eliashberg step disagree (one "ir", one "uniform"), an INFO message
  states the combination and that it is supported through the densified
  interface, so a forgotten flag is visible in the log.
- `ir_check` (default false): two-basis validation mode, Sec. 5.
- Fail-fast with an actionable install hint if `sparse-ir` is missing.

## 5. Accuracy management -- REVISED

Node-only fit residuals cannot see out-of-basis weight, so three layers:

1. ALWAYS ON, cheap: coefficient-decay diagnostic. For every fitted object
   the tail ratio max(|g_l|, l in last 10%) / max(|g_l|) is checked; above
   10*eps a warning names the object and remedy (raise ir_wmax / ir_tol).
2. ALWAYS ON at Stage-1 load: the full uniform grid acts as an oversampled
   residual check (3.2), reported per object.
3. OPT-IN `ir_check = true`: two-basis validation. Every
   `ir_check_interval`-th SCF iteration (and once per Eliashberg solve)
   the products chi0, Sigma (Stage 2) or the kernel matvec on a probe
   vector (Stage 1) are re-evaluated on a check basis (eps/100 and wmax*2,
   with oversampled check nodes); the max deviation is reported and a
   threshold aborts with guidance. The check obeys the SAME memory-aware
   chunking as production transforms (3.1); its transform matrices are
   built once per run and cached (host + device) alongside the production
   ones, so no per-check reallocation/transfer churn.

Equivalence testing (all vs the uniform path on the existing fixtures):
- chi0, Sigma, mu, NCond, lambda within eps-scaled tolerances;
- MONOTONIC convergence toward the uniform result as eps is tightened and
  wmax is raised (no fixed halving-rate claim);
- a wmax-sweep test demonstrating insensitivity of lambda / chi_s peak to
  wmax doubling at the default auto choice;
- small-Nmat fixtures assert accuracy only, never a compression ratio
  (at Lambda ~ 10, L ~ Nmat is expected).

## 6. Interplay with existing features

- GPU: transform matrices as device arrays; products already xp; chunked
  matmuls per 3.1.
- fft_workers: spatial FFTs unchanged.
- Anderson mixing: unchanged; coefficient-norm diagnostic added (3.3).
- Parity projection (dynamic Eliashberg): with integer-symmetric nodes the
  combined parity operator is an exact node permutation composed with the
  existing orbital/k transposes. Tests: (a) node-permutation correctness
  under n -> -n-1, (b) the full parity-classification suite (singlet/
  triplet/odd-frequency) re-run on the IR path, (c) permutation composed
  with k/orbital transposes matches the uniform operator on densified
  vectors.

## 7. Open questions and dispositions

- OQ-1 (output compatibility): RESOLVED by policy -- IR runs densify all
  outputs to the uniform grid by default with byte-identical metadata
  semantics (Nmat remains a required parameter); IR-native output deferred
  to Stage 3 after a downstream-reader audit. Loader round-trip tests pin
  the semantics (including the nmat//2 slice convention).
- OQ-2 (wmax auto-choice adequacy for product spectra): RESOLVED by
  policy + validation -- heuristic default, layered runtime diagnostics
  (Sec. 5) as the actual guard, wmax-sweep in the test plan, always
  overridable and logged.
- OQ-3 (node symmetry guarantees across sparse-ir versions): RESOLVED by
  construction -- symmetrize-and-rebuild at IRAxis construction with an
  integer-arithmetic assertion, so no dependence on library guarantees.
- OQ-4 (IR-vs-uniform tolerance expectations): RESOLVED by documentation --
  IR introduces a controlled systematic error at scale eps, distinct from
  backend round-off; all IR equivalence tolerances are eps-scaled and
  documented as such.
- R-1 (sparse-ir API stability): pin `sparse-ir >= 1.0`; all API
  touchpoints confined to ir_axis.py.
- R-2 (small-fixture compression): tests assert accuracy only (Sec. 5).
- R-3 (Anderson node-vs-coefficient convergence): coefficient-space change
  logged each iteration; documented as the wmax-inadequacy indicator.
- R-4 (transform batch memory): chunked application required by 3.1.

## 8. Validation & rollout

- Stage 1 PR: ir_axis.py + dynamic-Eliashberg IR path + matvec-vs-full-grid
  gate + equivalence tests (CPU/GPU) + docs; default OFF.
- Stage 2 PR: FLEX IR path + equivalence tests (incl. mu/NCond and the
  parity suite) + docs; default OFF.
- Production shakedown: one bp-ICl2 (P, T) point, IR vs uniform, lambda
  compared at eps-scaled tolerance, before enabling in sweeps.
