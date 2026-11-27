.. _initial_states:

***********************************
Candidate initial states (preview)
***********************************

A UHF calculation converges to a state that depends on its initial Green
function. ``python -m hwave.initial`` prepares a set of candidate initial
states for a UHFk input, one run directory each; you run ``hwave`` in each
directory, and the tool collects the results and ranks them by free energy.

.. note::

   This is a preview. The commands, options and file formats may still
   change. This release provides two tiers. The large tier is the default
   candidate set: spin and charge patterns at every q the cell can
   represent, plus the reserved rows ``para``, ``ferro`` and ``random``,
   which no quota cuts. The small tier adds orbital patterns and
   instability rows from a static analysis of the paramagnet (see
   :ref:`initial_states_small_tier`). A medium tier for larger models, and a
   tier choice from a calibrated cost model, arrive in a later release.

The selected state is the lowest free energy among the candidates tried. It is
not a proof of the ground state.

Workflow
========

.. code-block:: bash

   $ python -m hwave.initial prepare INPUT_DIR --out WORK
   $ bash WORK/commands.txt          # or: cd WORK && parallel < commands.txt
   $ python -m hwave.initial collect WORK

``python -m hwave.initial --help`` and ``python -m hwave.initial <command>
--help`` list every option and the exit codes.

``prepare``
-----------

``INPUT_DIR`` is a directory holding exactly one UHFk input TOML file, or the
TOML file itself. Options:

- ``--out WORK`` (required): the work directory.
- ``--subshape Bx,By,Bz`` (repeatable): evaluate every candidate on each of
  these ``SubShape`` cells. Default: the input's ``SubShape``.
- ``--extra-candidates FILE.json``: further candidates, a JSON list of
  ``{"family": "user" | "literature", "pattern": ..., "label": ...,
  "citation": ...}``. ``pattern`` follows the JSON schema installed with the
  package, ``hwave/initial/patternspec.schema.json`` (in Python:
  ``importlib.resources.files("hwave.initial") / "patternspec.schema.json"``);
  its orbital indices are physical orbital indices (also in spin-orbital
  mode).
- ``--t-hint T``: the hint temperature, in the energy units of the input.
  The small tier evaluates its hints at this temperature (see
  :ref:`initial_states_small_tier`); it changes only the instability rows
  and their order, never the other candidates, and the candidate runs stay
  at the input ``T``. The large tier only records it. Default: the input
  ``T`` when ``T > 0``; at ``T = 0`` a rule from the bands and the k mesh,
  ``max(W/(2 L), 1e-4 W)``, with ``W`` the energy span of the bands around
  the Fermi level (at least 1e-3) and ``L`` the shortest ``CellShape`` axis
  longer than 2. The rule used is shown in the summary.
- ``--amplitude A`` (the requested seed amplitude, in (0, 1), default 0.5),
  ``--max-candidates N`` (generated rows per ``SubShape``, default 32),
  ``--orbital-site FILE`` (a JSON list with one site label per physical
  orbital).
  At ``T = 0`` the seed field can reach the requested amplitude exactly
  where two levels cross at the Fermi level. The seed is then taken at
  that crossing, with the crossing levels filled uniformly when that meets
  the request, or else with the occupation just past the crossing.
  ``amplitude_used`` in ``coverage.json`` can therefore exceed the request.
  On a spin-symmetric model, a pattern without a spin component (para,
  charge, orbital) gives a seed with identical spin-up and spin-down parts.
- ``--tier auto|large|small``: the candidate set (``medium`` is reserved
  for a later release). ``auto``, the default, chooses the small tier when
  the model is small enough (see :ref:`initial_states_auto_tier`), else the
  large tier.
- ``--para-check`` (small tier): the preliminary paramagnetic run also
  starts from the uniform orbital occupation, and the lower free energy is
  kept (see :ref:`initial_states_small_tier`). Refused with ``--tier
  large``; with ``auto`` choosing the large tier it is ignored with a
  warning.
- ``--stage0`` (small tier): the preliminary paramagnetic run is not run
  inside ``prepare`` but written as run directories of its own, for a
  batch system (see :ref:`initial_states_preliminary_run`). Refused when the
  large tier is given or chosen.

``prepare`` writes into WORK:

- ``coverage.json``: the coverage report before the runs (see below);
- ``manifest.json``, ``candidates.json``: what was prepared, with the
  SHA-256 of each run's ``input.toml``, ``initial.npz`` and copied inputs;
- ``commands.txt``: a first line that changes to WORK, a few comment lines
  starting with ``#``, then one line per run, relative to WORK (a job array
  that picks lines by number must skip the first lines);
- one directory per run, holding ``input.toml``, ``initial.npz`` and
  ``inputs/``.

Each run's ``input.toml`` is the input rewritten as follows:
``[file.input] path_to_input = ""`` and ``initial = "initial.npz"``;
``[mode.param] SubShape`` set (an existing key is replaced in its own
spelling); ``initial_uhf`` and ``initial_mode`` removed; every input file
(the interaction files, ``geometry_uhf``, ``onebodyg_uhf``) **copied** into
``inputs/``, with ``[file.input.interaction] path_to_input = "inputs"``;
``[file.output] energy = "energy.dat"`` added when absent. Comments of the
original TOML are not kept. Refused: absolute output paths or paths with
``..``, ``green = ""`` or ``energy = ""``, and output names that would
overwrite the run directory's own files. Because the inputs are copied, not
linked, WORK can be moved, or its sources deleted, and still be run and
collected.

Keys are read with the solver's case rules: the tables ``[mode]``,
``[file]``, ``[file.input]``, ``[file.output]`` and their keys are
case-sensitive, ``[mode.param]`` and ``[file.input.interaction]`` keys are
not. A case variant of a known key in a case-sensitive table (for example
``[file.input] Initial``) is refused with the correct spelling, because the
solver would ignore it.

Running ``prepare`` again with the same command and inputs on its own WORK
does nothing; a different command, changed inputs or a non-empty directory
not written by ``prepare`` are refused, naming the difference.

Complex coefficients are accepted only for ``PairHop``. ``prepare`` refuses
a complex coefficient of ``CoulombInter``, ``Hund``, ``Ising``,
``Exchange``, ``PairLift`` or the combined ``Coulomb`` table, naming the
term, the bond and the orbitals: for these terms the solver's mean-field
potential with a complex coefficient is not that of the documented
Hamiltonian, so no candidate built on it would be meaningful.

With a fixed ``2Sz``, the solver keeps only the spin-conserving parts of the
interaction terms. When an interaction term breaks spin conservation about z
(for example ``PairLift``), ``prepare`` records a warning naming the terms,
and the summary shows it.

.. _initial_states_small_tier:

The small tier
--------------

The small tier adds two families to the large-tier set:

- ``orbital``: the orbital polarisation ``n_a - n_b`` of consecutive
  interacting orbitals ``a``, ``b`` of a site, uniform and with the site
  sign patterns, at every q the cell can represent (counted in
  ``--max-candidates`` like the baseline rows);
- ``instability``: patterns from a static linear-response analysis of the
  paramagnet, the hints. They order the instability rows only: no hint
  removes a candidate of another family.

The hints are computed in four steps.

1. **The preliminary paramagnetic run.** UHFk runs on the primitive cell
   (``SubShape = [1, 1, 1]``) at the hint temperature ``T_hint``, with the
   input's ``Ncond``, ``2Sz``, ``flag_fock`` and interaction tables, from the
   free density of the bare one-body Hamiltonian, with ``EPS`` raised to at
   least 12 (recorded as ``para_eps``). The result is the reference of the
   hints only when the run converged and keeps the symmetries of the model
   within 1e-8: the spin rotations that the one-body Hamiltonian and the
   interactions both keep, the electron count of each sector with a fixed
   ``2Sz``, and time reversal in spin-orbital mode when the model has it.
   Otherwise the hints are skipped with the reason (``hint.skip_reason``)
   and every seed is built on the bare one-body Hamiltonian. With
   ``--para-check`` a second run starts from the uniform orbital
   occupation; the lower free energy is kept, and a difference above 1e-8
   per cell is reported as a warning (the other start reached a
   metastable, orbital-ordered state). The reference is kept in
   ``WORK/para_hf/`` (``green.npz``, ``para_hf.json`` and the runs), so
   ``collect`` needs no recomputation. The seeds of every family are then
   built on its Hamiltonian ``H0 + Sigma'`` folded to each cell
   (``seed_base = "para_hf"``), and the labels of ``prepare`` and
   ``collect`` are taken relative to it.
2. **The vertex** ``Gamma(q)``: the derivative of the solver's Hartree-Fock
   potential with respect to the density, on the one-body slots the solver
   reads and writes (the pair space ``P``, restricted to the interacting
   orbitals). It follows the input's ``flag_fock``: the Hartree terms
   always, the Fock terms only when the runs have them, so bond-order
   patterns from a ``CoulombInter`` bond appear only with the Fock term.
3. **The bare susceptibility** ``chi0(q)``: the static Lindhard response of
   the paramagnetic bands at ``T_hint``, on the same pair space. At
   ``q = 0`` it is taken at the fixed electron number of each sector (with
   a fixed ``2Sz`` each spin, otherwise the whole system, non-interacting
   orbitals included), as in the solver, which refits the chemical
   potential at every iteration. A uniform change of that number is then
   no instability (``lambda = 0``); with a fixed ``2Sz`` this includes the
   uniform ``S^z``.
4. **The eigenvalues** ``lambda`` of ``K(q) = -chi0(q) Gamma(q)`` at every q
   of the ``CellShape`` mesh, the ``top_k = 3`` largest positive ones per q
   (when the spin-rotation check passes, the charge and ``S^z`` sectors are
   solved separately and merged before the cut). Each eigenvector ``x``
   gives an instability row whose seeding pattern is the conjugate field
   ``O = Gamma x``.

``lambda`` is the static Stoner-type factor of that pattern: at ``lambda >=
1`` the paramagnet at ``T_hint`` is unstable towards it within Hartree-Fock,
and a larger ``lambda`` means a stronger tendency. ``lambda < 1`` at
``T_hint`` is normal and does not mean that the pattern cannot order at the
input ``T``. Each instability row records ``rank`` (its place in the
ordering), ``lambda`` and ``lambda_imag`` (0 unless the vertex is not
Hermitian, which is reported). Rows are ordered by decreasing ``lambda``;
eigenvalues within a relative 1e-8 of each other form one degenerate set,
ordered by a fixed rule (q in a fixed order, then the position of the
eigenvector in that q's eigenvalue order), and a degenerate set at one q
is given a canonical basis, so the rows and their ids do not depend on the
round-off of the eigensolver. Every hint row is listed; at most 12 rows
(``quota_instability``) are run per ``SubShape``, and the others say why
(``quota: instability``, or ``unrepresentable`` with the cell that would
represent their q).

**The hint temperature smears the hints.** Each transition enters
``chi0`` with the weight ``(f(e) - f(e'))/(e' - e) <= 1/(4 T_hint)``, so
structure within about ``T_hint`` of the Fermi level or from energy
differences below ``T_hint`` (small gaps, imperfect nesting, split van
Hove peaks) is suppressed: ``lambda`` drops, peaks broaden, the leading q
can shift and close rows can swap. The summary repeats this sentence with
the value used.

The hint record of ``coverage.json`` (``hint``) holds the hint temperature
and its rule, the status of the paramagnetic run (``para_status``,
``para_free_energy``, ``para_residual``, ``para_polarisation``,
``para_eps``), whether the vertex is Hermitian, the method, the q set and
the status of each q (``ok``, ``fewer than top_k positive eigenvalues`` or
``no positive eigenvalue``).

.. _initial_states_auto_tier:

Choosing the tier
-----------------

In this release ``--tier auto`` depends on the model only, so an input gets
the same tier on every host: the small tier when the pair space has
``|P| <= 400`` slots, the mesh has ``N_q <= 4096`` points and ``N_q |P|^3 <=
1e11``; otherwise the large tier, with a warning. ``--tier small`` given
explicitly runs even beyond these limits, with a warning. The host memory
never chooses the tier: it only refuses a small tier whose memory estimate,
``16 (|P|^2 + N_q top_k |P| + N_k nb^2 + N_k nb |P|)`` bytes (``N_k`` k points,
``nb`` bands), exceeds the available memory. The ``tier`` record of
``coverage.json`` gives the sizes, the estimate and the rule.

.. _initial_states_preliminary_run:

The preliminary paramagnetic run on a batch system
--------------------------------------------------

The paramagnetic run of the small tier normally runs inside ``prepare``.
For a large model it can run as a job of its own:

.. code-block:: bash

   $ python -m hwave.initial prepare INPUT_DIR --out WORK --tier small --stage0
   $ bash WORK/stage0_commands.txt      # exit code 3 above: run this first
   $ python -m hwave.initial prepare INPUT_DIR --out WORK --tier small --stage0

The first ``prepare`` writes ``WORK/stage0_para/`` (and
``WORK/stage0_para_check/`` with ``--para-check``), staged like the candidate
runs (``SubShape = [1, 1, 1]``, ``T = T_hint``, ``EPS`` raised to at least
12), and ``WORK/stage0_commands.txt``, and exits with code 3. The same
``prepare`` command then checks those inputs against ``manifest.json`` (a
changed file is refused, naming it) and reads the run: when it converged
and keeps the symmetries of the model (step 1 of
:ref:`initial_states_small_tier`) it becomes the reference and the
candidate runs are written; when it ran and failed, did not converge or
broke one of those symmetries, the seeds use the bare one-body Hamiltonian
and the hints are skipped with that reason; when it has not run yet, or has not
finished, ``prepare`` prints the commands again, changes nothing and exits
with code 3. ``--stage0`` is part of the command: a WORK prepared with it
re-enters only with it.

The initial-state files
-----------------------

Each ``initial.npz`` holds the Green function UHFk reads, and its
provenance: the candidate and run ids, the label, q, phase and sign of the
pattern, ``SubShape``, ``CellShape``, ``Ncond``, ``2Sz``, the field
strength used (``h_used``), the SHA-256 of the inputs, and the base of the
seed:
``seed_base`` is ``"para_hf"`` (the paramagnetic reference of the small
tier) or ``"bare"``; on the reference it also records the temperature of
that reference (``seed_base_para_T``) and the digest of its inputs
(``seed_base_para_digest``). The library writer refuses a reference that
belongs to another input or another hint temperature.

Running the candidates
----------------------

``bash WORK/commands.txt`` works from any directory. Each line enters one
run directory, checks its initial-state file and runs ``hwave input.toml >
hwave.log 2>&1``. For GNU parallel or job arrays, run the lines from WORK.

``check``
---------

.. code-block:: bash

   $ python -m hwave.initial check CONFIG.toml [NPZ] [--require-provenance]

checks the initial-state file that ``CONFIG.toml`` reads (or ``NPZ``) before
``hwave`` runs: provenance, ``SubShape``, shape, Hermiticity, eigenvalue
bounds, electron number and the solver's block structure. A file without
provenance is accepted as unverified; with ``--require-provenance`` it is
refused with exit code 4.

``collect``
-----------

``collect WORK`` reads every run and writes ``WORK/coverage.json`` (after the
runs) and ``WORK/summary.txt``; it can be repeated. Each run gets a status:

- ``seed_unverified``: ``input.toml`` no longer names the prepared
  ``initial.npz``, or a copied input or the npz differs from the file
  ``prepare`` wrote, or the check refuses it;
- otherwise the run's ``scf_status.dat`` decides when it exists:
  ``converged``, ``not_converged``, or ``failed`` when no SCF step was
  performed. A status file that exists but cannot be read makes the run
  ``failed`` (``"scf_status.dat unreadable: ..."``); it is not ranked;
- without a status file: ``not_run`` (no ``hwave.log`` and no output),
  ``convergence_unverified`` (outputs without ``hwave.log``, ranked in a
  separate flagged list), ``failed`` (a traceback or a solver error),
  ``not_converged`` or ``converged`` from ``hwave.log``. A run that was
  killed (no status file, no traceback) is ``not_converged`` with the
  reason ``hwave.log has no convergence message``.

The outputs of a converged run must be readable, otherwise it is
``failed``. ``collect`` exits with 0 also when some runs failed or did not
converge; read their status in the report.

``summary.txt`` holds the model (orbitals, filling, ``T``, ``2Sz``,
``CellShape``), the tier and the hint temperature with its rule, the label
thresholds, the number of rows and why the others were not tried, the
ranked table (id, label, free energy per cell, difference per site from the
lowest row, status, final label), the flagged rows and warnings, the limits
below, and the sentence quoted at the top of this page.

Exit codes
==========

====  ======================================================================
code  meaning
====  ======================================================================
0     done (``check``: verified, or unverified without
      ``--require-provenance``)
1     an error, or the check refused the file
2     an argument error
3     ``prepare --stage0``: a preliminary paramagnetic run is pending; run
      ``bash WORK/stage0_commands.txt``, then the same ``prepare`` command
4     ``check --require-provenance``: the file has no provenance
====  ======================================================================

The coverage report
===================

``coverage.json`` lists one row per candidate and ``SubShape``. Rows are
listed, never dropped: a row that is not tried says why (``duplicate of
<run_id> (<map>)`` for a state equal to another one by a symmetry,
``quota: <family>``, or ``unrepresentable`` with the reason and a proposed
cell, e.g. ``q_x denominator 3 does not divide B_x = 2`` and a 3 × 3
proposal).

- The families are ``para`` and ``ferro`` (reserved, one run per
  ``SubShape``), ``random``, ``baseline`` (spin and charge patterns at the q
  vectors the cell can represent), ``orbital`` and ``instability`` (the
  small tier), and ``user`` and ``literature`` (from
  ``--extra-candidates``; never cut by a quota).
- A candidate is named by its seeding pattern, e.g. ``spin-z ferro-orbital
  q=(1/2,1/2,0) φ=0``. The state a run ends in has a ``final_label``,
  naming the components that hold at least 10% of the weight; a part is
  called time-reversal (TR) even or odd when the other parity holds less
  than 10% of the weight, otherwise mixed. The ordering vector follows the
  on-site components it describes; bond components come after it, e.g.
  ``spin-z ferro-orbital q=(1/2,1/2,0) + bond-charge`` for the Néel state
  with its uniform change of the hopping. The thresholds are recorded in the
  report.
- ``collect`` ranks the converged runs by free energy per site (per cell
  when the number of sites is not defined) within a comparison group of
  equal ``T``, ``Ncond`` per cell, ``CellShape`` and ``2Sz``. Rows within
  1e-4 per site of each other (in the energy units of the input) are
  flagged ``near-degenerate``. Unconverged
  and failed rows are listed with their reason, not ranked.

Limits
======

UHF overestimates order and ordering temperatures, and the k mesh and the
cells are finite. No candidate covers: multi-q superpositions (double-Q,
triple-Q, vortex crystals); sliding-phase variants beyond the two kept
phases; combinations within a degenerate eigenspace other than the canonical
basis; rotated cells such as √3 × √3, and cells outside the given
``SubShape`` list; orbital and bond patterns in the large tier, and bond
patterns other than the hints; more than one amplitude per candidate;
spirals unless the spin-rotation check passes; q vectors with a
denominator above 16. The summary repeats these limits.

The hints describe the paramagnet only, statically and at ``T_hint``. A
state reached by a first-order transition or at strong coupling may have
no hint at all; this is why the reserved and baseline rows always run and
a hint never removes a row. A small ``lambda`` is not evidence against a
state, and the largest ``lambda`` does not select the result: only the
free energies of the converged runs do, within the candidates tried.

Only UHFk inputs are supported by the command line. The library function
``hwave.initial.initial_state`` also writes a UHFr ``Initial`` file from a
UHFk input of the same lattice and a site map. Spin-orbital mode is
supported by UHFk only.
