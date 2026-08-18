#!/usr/bin/env python3
"""Render the RPA/FLEX support and equivalence page from the cell registry.

This is an EXPLICIT generator, run by a developer:

.. code-block:: bash

    python docs/tools/render_equivalence_table.py            # print to stdout
    python docs/tools/render_equivalence_table.py --write    # update the page

``render()`` is a pure function of ``tests/equivalence_cells.py``'s
``CELLS`` and ``PROVENANCE`` records: cells are emitted in ``cell_id``
order, floats through fixed format strings, and nothing is read from the
clock, the environment, or a live solver run. Two calls therefore produce
byte-identical text, and ``TestRenderedEquivalenceTableDrift`` (in
``tests/test_rpa_flex_equivalence_table.py``) fails whenever the committed
page and the registry disagree.

Two constraints shape the design:

  * The Sphinx build NEVER imports the registry -- the committed RST is
    static, so a source distribution without the test tree still builds
    the documentation. Only this script and the drift test import it.
  * The registry is never mutated here. This module only reads records
    (they are frozen dataclasses with read-only mapping views anyway).

The registry's ``PROVENANCE`` record is rendered into the status section
and is the ONLY thing the freeze commit changes: flipping ``status`` from
"candidate" to "frozen" (and filling in the measured source commit and
the workflow run IDs) is a metadata-only regeneration of this page.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import List, Mapping, Sequence

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from tests.equivalence_cells import (  # noqa: E402  (needs the path bootstrap above)
    CELLS,
    PROVENANCE,
    Cell,
    Diverges,
    Equiv,
    ExecuteReject,
    Site,
    Status,
    relationship,
)

OUTPUT_PATH = os.path.join(
    _REPO_ROOT, "docs", "en", "source", "algorithm", "rpa_flex_equivalence.rst"
)

_STATUS_LABEL = {
    Status.SUPPORTED: "supported",
    Status.REJECT: "rejected",
    Status.NOT_APPLICABLE: "not applicable",
}

_INTERACTION_CLASS_LABEL = {
    "onsite": "on-site",
    "offsite": "off-site",
}

_SITE_LABEL = {
    Site.CONSTRUCTOR: "at construction",
    Site.SOLVE: "during the solve",
}


# ---------------------------------------------------------------------------
# Small deterministic formatting helpers
# ---------------------------------------------------------------------------


def _heading(text: str, underline: str) -> List[str]:
    return [text, underline * len(text), ""]


def _literal(text: str) -> str:
    return "``{}``".format(text)


def _shape(shape) -> str:
    return "x".join(str(n) for n in shape)


def _number(value: float) -> str:
    return repr(float(value))


def _tolerance(value: float) -> str:
    return "{:.1e}".format(value)


def _unique_sorted(values) -> list:
    """Sorted distinct values -- sorting the VALUES, not their rendered
    text, so numeric envelope entries come out in numeric order.
    """

    return sorted(set(values))


def _interactions(cell) -> str:
    types = ", ".join(_literal(name) for name in sorted(cell.fixture.interactions))
    label = _INTERACTION_CLASS_LABEL.get(cell.interaction_class, cell.interaction_class)
    return "{} ({})".format(types, label)


def _scheme(cell) -> str:
    requested = cell.fixture.requested_scheme
    if requested == cell.resolved_scheme:
        return _literal(cell.resolved_scheme)
    return "{} -> {}".format(_literal(requested), _literal(cell.resolved_scheme))


def _table_row(values: Sequence[str]) -> List[str]:
    lines = ["   * - {}".format(values[0])]
    lines.extend("     - {}".format(value) for value in values[1:])
    return lines


def _definition(term: str, body: Sequence[str]) -> List[str]:
    lines = [term]
    lines.extend("   {}".format(line) if line else "" for line in body)
    lines.append("")
    return lines


# ---------------------------------------------------------------------------
# Section renderers
# ---------------------------------------------------------------------------


def _header_comment() -> List[str]:
    return [
        ".. This page is GENERATED from tests/equivalence_cells.py by",
        "   docs/tools/render_equivalence_table.py. Do not edit it by hand:",
        "   re-run `python docs/tools/render_equivalence_table.py --write` and",
        "   commit the result. A test in",
        "   tests/test_rpa_flex_equivalence_table.py fails whenever this page",
        "   and the registry disagree.",
        "",
        ".. _rpa_flex_equivalence:",
        "",
    ]


def _introduction() -> List[str]:
    return [
        "H-wave computes the susceptibilities of a given model along two",
        "routes: the RPA solver and the one-shot (first-iteration) FLEX",
        "solver. Where the two routes describe the same approximation they",
        "are expected to produce the same arrays, and this page records, for",
        "a fixed set of small reference inputs, which configurations each",
        "solver accepts and how closely the two agree where both run.",
        "",
        "The comparison is one of IMPLEMENTATION, not of physics: a row",
        "stating that the two solvers agree says that the two code paths",
        "produce the same numbers on that input, not that either result is a",
        "better description of the material. Likewise, a configuration",
        "marked as rejected is one the solver refuses to run today; the",
        "entry records the behaviour, not a statement about whether the",
        "quantity is meaningful.",
        "",
        "Every claim on this page is checked by",
        "``tests/test_rpa_flex_equivalence_table.py``, which runs both",
        "solvers on each input listed below and asserts the recorded",
        "outcome. The page itself is generated from the same records, so it",
        "cannot drift away from what the tests assert.",
        "",
    ]


def _status_section(provenance: Mapping) -> List[str]:
    lines = _heading("Status of these records", "*")
    status = provenance["status"]
    if status == "candidate":
        lines.extend([
            "The tolerances below are **provisional**. They were measured on",
            "development machines and have not yet been confirmed on the",
            "project's continuous-integration runners, so the exact numbers",
            "may still move. The per-tolerance notes name the machines and",
            "the source revision each measurement came from.",
            "",
        ])
    elif status == "frozen":
        source_sha = provenance["source_sha"]
        if not source_sha:
            raise ValueError(
                "a frozen provenance record must name the measured source "
                "revision, got source_sha={!r}".format(source_sha)
            )
        run_ids = tuple(provenance["run_ids"])
        sentence = (
            "The tolerances below are **confirmed**: they were reproduced "
            "on the project's continuous-integration runners at source "
            "revision {}".format(_literal(str(source_sha)))
        )
        if run_ids:
            sentence += " (workflow runs {})".format(
                ", ".join(_literal(str(run_id)) for run_id in run_ids)
            )
        sentence += (
            ", and hold as of that revision. The notes below record where "
            "each residual was originally measured."
        )
        lines.extend([sentence, ""])
    else:
        raise ValueError(
            "unknown provenance status {!r}: expected 'candidate' or "
            "'frozen'".format(status)
        )
    return lines


def _reading_section() -> List[str]:
    lines = _heading("How to read the table", "*")
    lines.extend([
        "The *RPA* and *FLEX* columns give each solver's outcome on the",
        "input, and the *Relationship* column summarises the pair:",
        "",
        "``EQUIV``",
        "   Both solvers run and every compared array agrees to the recorded",
        "   tolerance.",
        "",
        "``DIVERGES(<issue>)``",
        "   Both solvers run, but at least one compared array differs by more",
        "   than the agreement policy allows. The difference is bounded and",
        "   pinned by the tests, and the named issue tracks unifying the two",
        "   implementations.",
        "",
        "``AUTO-RESOLVES(<scheme>)``",
        "   The input leaves ``calc_scheme = \"auto\"``, and both solvers",
        "   resolve it to the named scheme. These rows check the resolution",
        "   only; no susceptibility is computed.",
        "",
        "The remaining labels",
        "   No comparison is possible, because at least one solver does not",
        "   run the input: ``BOTH-REJECT``, ``FLEX-REJECT·RPA-SUPPORTED``,",
        "   ``RPA-REJECT·FLEX-SUPPORTED``, ``RPA-ONLY``, ``FLEX-ONLY``,",
        "   ``RPA-REJECT`` and ``FLEX-REJECT``, naming the side that does",
        "   not. The next section says what happens instead.",
        "",
    ])
    return lines


def _support_matrix(cells) -> List[str]:
    lines = _heading("Support matrix", "*")
    lines.extend([
        ".. list-table::",
        "   :header-rows: 1",
        "   :widths: 26 20 10 9 8 8 19",
        "",
    ])
    lines.extend(_table_row([
        "Configuration", "Interactions", "Scheme", "Spin mode", "RPA", "FLEX",
        "Relationship",
    ]))
    for cell in cells:
        lines.extend(_table_row([
            _literal(cell.cell_id),
            _interactions(cell),
            _scheme(cell),
            cell.expected_spin_mode,
            _STATUS_LABEL[cell.rpa.status],
            _STATUS_LABEL[cell.flex.status],
            _literal(relationship(cell)),
        ]))
    lines.append("")
    return lines


def _reject_step(proof) -> ExecuteReject:
    """The proof's single ``ExecuteReject`` step. The registry schema
    admits exactly one for a REJECT side; a record that violates that is
    a generation-time error rather than a silently mis-rendered row.
    """

    for step in proof.steps:
        if isinstance(step, ExecuteReject):
            return step
    raise ValueError("a rejected side carries no ExecuteReject step")


def _side_outcome(solver: str, proof) -> List[str]:
    if proof.status is Status.REJECT:
        step = _reject_step(proof)
        return ["* {} refuses the input {} and raises {} whose message "
                "contains {}.".format(
                    solver,
                    _SITE_LABEL[step.site],
                    _literal(step.exc_type),
                    _literal(step.fragment),
                )]
    if proof.status is Status.NOT_APPLICABLE:
        return ["* {}: nothing to compare -- {}".format(solver, proof.reason)]
    return []


def _unsupported_section(cells) -> List[str]:
    lines = _heading("Configurations one solver does not run", "*")
    lines.extend([
        "For these inputs the two solvers cannot be compared. The entries",
        "record what the unsupported side does instead.",
        "",
    ])
    for cell in cells:
        body: List[str] = []
        body.extend(_side_outcome("RPA", cell.rpa))
        body.extend(_side_outcome("FLEX", cell.flex))
        if body:
            lines.extend(_definition(_literal(cell.cell_id), body))
    return lines


def _comparison_entries(comparison):
    """(observable, comparator, bound text, provenance) for each observable
    of a comparison proof, in a deterministic order (diverging observables
    first, each group sorted by name).
    """

    entries = []
    if isinstance(comparison, Equiv):
        for name in sorted(comparison.observables):
            spec = comparison.observables[name]
            entries.append((
                name,
                spec.comparator,
                "agree to {}".format(_literal(_tolerance(spec.atol))),
                spec.provenance,
            ))
    elif isinstance(comparison, Diverges):
        for name in sorted(comparison.diverging):
            spec = comparison.diverging[name]
            entries.append((
                name,
                spec.comparator,
                "differ by more than the agreement policy allows; the "
                "difference is pinned to the bracket ({}, {}], and {} tracks "
                "unifying the two implementations".format(
                    _literal(_tolerance(spec.ceiling)),
                    _literal(_tolerance(spec.regression_bound)),
                    comparison.step5_issue,
                ),
                spec.provenance,
            ))
        for name in sorted(comparison.others):
            spec = comparison.others[name]
            entries.append((
                name,
                spec.comparator,
                "agree to {}".format(_literal(_tolerance(spec.atol))),
                spec.provenance,
            ))
    return entries


def _tolerances_section(cells) -> List[str]:
    lines = _heading("Recorded agreement", "*")
    lines.extend([
        "Arrays are compared elementwise on complex values: the two results",
        "count as agreeing when ``|a - b| <= atol`` holds for every element,",
        "with no relative term. ``chi0q`` is the bare susceptibility and",
        "``chiq`` the interacting one; the comparator names the mapping",
        "applied before the comparison, which reconciles the two solvers'",
        "storage layouts (``identity`` when both already store the same",
        "shape).",
        "",
        "Each entry also carries the measurement note recorded with the",
        "tolerance -- the machines, the source revision and the observed",
        "differences it was derived from. Those notes are maintenance",
        "records rather than user guidance.",
        "",
    ])
    for cell in cells:
        entries = _comparison_entries(cell.comparison)
        if not entries:
            continue
        body: List[str] = []
        for name, comparator, bound, provenance in entries:
            body.append("* {} (comparator {}): {}.".format(
                _literal(name), _literal(comparator), bound
            ))
            body.append("  Provenance: {}".format(provenance))
        lines.extend(_definition(_literal(cell.cell_id), body))
    return lines


def _links_section(cells) -> List[str]:
    lines = _heading("Related tests", "*")
    lines.extend([
        "A maintenance index of tests elsewhere in the suite that cover the",
        "same ground, each with the note recorded alongside it. They are",
        "listed for orientation only: the claims on this page are",
        "established by ``tests/test_rpa_flex_equivalence_table.py`` alone.",
        "",
    ])
    for cell in cells:
        body: List[str] = []
        for solver, proof in (("RPA", cell.rpa), ("FLEX", cell.flex)):
            for link in proof.links:
                body.append("* {} -- {}: {}".format(
                    solver, _literal(link.test_id), link.claim
                ))
        if body:
            lines.extend(_definition(_literal(cell.cell_id), body))
    return lines


def _envelope(cells) -> List[str]:
    return [
        "* Lattices {} with sublattices {}.".format(
            ", ".join(_literal(_shape(s)) for s in _unique_sorted(
                cell.fixture.CellShape for cell in cells)),
            ", ".join(_literal(_shape(s)) for s in _unique_sorted(
                cell.fixture.SubShape for cell in cells)),
        ),
        "* Temperatures {}.".format(
            ", ".join(_literal(_number(t)) for t in _unique_sorted(
                cell.fixture.T for cell in cells)),
        ),
        "* Matsubara frequency counts {}.".format(
            ", ".join(_literal(str(n)) for n in _unique_sorted(
                cell.fixture.Nmat for cell in cells)),
        ),
        "* Diagram selections {}.".format(
            ", ".join(_literal(t) for t in _unique_sorted(
                cell.fixture.calc_type for cell in cells)),
        ),
        "* Calculation schemes {}.".format(
            ", ".join(_literal(scheme) for scheme in _unique_sorted(
                cell.resolved_scheme for cell in cells)),
        ),
        "* Both a fixed chemical potential and a fixed electron count.",
    ]


def _scope_section(cells) -> List[str]:
    lines = _heading("Scope", "*")
    lines.extend([
        "The table certifies the inputs it lists and nothing beyond them.",
        "All of them are deliberately small, so that the whole set runs as",
        "part of the ordinary test suite:",
        "",
    ])
    lines.extend(_envelope(cells))
    lines.extend([
        "",
        "Agreement measured on inputs of this size does not by itself",
        "establish agreement on a production-sized model, and the tolerances",
        "are the measured round-off levels of these inputs, not general",
        "accuracy guarantees.",
        "",
    ])
    lines.extend(_heading("Accelerated backends", "-"))
    lines.extend([
        "The table makes no statement about the GPU backends: every input is",
        "run on the CPU backend of both solvers. That the GPU backend",
        "reproduces the CPU one is a separate, per-solver property, covered",
        "by ``tests/test_gpu_backend.py``, ``tests/test_rpa_gpu.py``,",
        "``tests/test_flex_gpu.py`` and ``tests/test_sc_gpu.py``.",
        "",
    ])
    lines.extend(_heading("FLEX-only numerical representations", "-"))
    lines.extend([
        "The intermediate-representation (IR) route has no counterpart on",
        "the RPA side, so it cannot appear as a row here: comparing the IR",
        "route with the dense one is a comparison of FLEX with itself. That",
        "comparison is covered by ``tests/test_flex_ir.py`` and",
        "``tests/test_flex_ir_general.py``.",
        "",
    ])
    return lines


# ---------------------------------------------------------------------------
# render
# ---------------------------------------------------------------------------


def render(cells: Sequence[Cell] = CELLS, provenance: Mapping = PROVENANCE) -> str:
    """Render the equivalence page as RST text (LF line endings, trailing
    newline). Pure: the same records always render the same bytes.
    """

    ordered = sorted(cells, key=lambda cell: cell.cell_id)

    lines: List[str] = []
    lines.extend(_header_comment())
    lines.extend(_heading("RPA and FLEX: support and equivalence", "="))
    lines.extend(_introduction())
    lines.extend(_status_section(provenance))
    lines.extend(_reading_section())
    lines.extend(_support_matrix(ordered))
    lines.extend(_unsupported_section(ordered))
    lines.extend(_tolerances_section(ordered))
    lines.extend(_links_section(ordered))
    lines.extend(_scope_section(ordered))

    while lines and lines[-1] == "":
        lines.pop()
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--write",
        action="store_true",
        help="write docs/en/source/algorithm/rpa_flex_equivalence.rst "
             "instead of printing to standard output",
    )
    args = parser.parse_args(argv)

    text = render()
    if args.write:
        with open(OUTPUT_PATH, "w", encoding="utf-8", newline="") as stream:
            stream.write(text)
        sys.stderr.write(
            "wrote {}\n".format(os.path.relpath(OUTPUT_PATH, _REPO_ROOT))
        )
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
