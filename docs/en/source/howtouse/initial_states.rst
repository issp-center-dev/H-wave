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
   change, and this release provides the large tier only. The large tier is
   the default candidate set: spin and charge patterns at every q the cell
   can represent, plus the reserved rows ``para``, ``ferro`` and ``random``,
   which no quota cuts. Later releases add smaller candidate sets chosen
   from an instability analysis of the paramagnet.

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
  It is recorded in the report only and never changes the candidate set or
  the seeds. Default: the input ``T`` when ``T > 0``; at ``T = 0`` a rule
  from the bands and the k mesh, ``max(W/(2 L), 1e-4 W)``, with ``W`` the
  energy span of the bands around the Fermi level (at least 1e-3) and ``L``
  the shortest ``CellShape`` axis longer than 2. The rule used is shown in the summary.
- ``--amplitude A`` (the requested seed amplitude, in (0, 1), default 0.5),
  ``--max-candidates N`` (generated rows per ``SubShape``, default 32),
  ``--orbital-site FILE`` (a JSON list with one site label per physical
  orbital), ``--tier`` (``auto`` resolves to ``large`` in this release).

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

With a fixed ``2Sz``, the solver keeps only the spin-conserving parts of the
interaction terms. When an interaction term breaks spin conservation about z
(for example ``PairLift``), ``prepare`` records a warning naming the terms,
and the summary shows it.

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
3     reserved for a later release (a preliminary paramagnetic run);
      never returned now
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
  vectors the cell can represent), and ``user`` and ``literature`` (from
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
``SubShape`` list; orbital and bond patterns in the large tier; more than
one amplitude per candidate; spirals unless the spin-rotation check passes;
q vectors with a denominator above 16. The summary repeats these limits.

Only UHFk inputs are supported by the command line. The library function
``hwave.initial.initial_state`` also writes a UHFr ``Initial`` file from a
UHFk input of the same lattice and a site map. Spin-orbital mode is
supported by UHFk only.
