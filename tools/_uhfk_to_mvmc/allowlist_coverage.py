"""Static allowlist coverage checker.

Iterates every case_*/input.toml under a fixture root, calls the CLI
allowlist predicate on the parsed (theta, sub_shape, cell_shape,
is_soc_mode) tuple, and enforces:

  * Every case_*/input.toml is parsed via
    tools._uhfk_to_mvmc.input_loader.load_input_toml -- the same
    shallow [mode] + [mode.param] merge the production CLI
    (tools/uhfk_to_mvmc.py) uses to source is_soc_mode. A fixture
    that is malformed (no [mode]/[mode.param], TOML parse error, or
    missing CellShape / BoundaryCondition) fails the whole check
    loudly (AllowlistCoverageError); it is never silently
    skipped.

  * Every discovered fixture must be accepted by the CLI predicate
    (tools._uhfk_to_mvmc.allowlist_predicate.is_supported_triple) --
    the exact function the CLI itself calls, so the checker and the
    CLI cannot drift.

  * If require_all_expected=True, every key in the passed expected_map
    must correspond to exactly one discovered fixture (missing -> fail),
    and its recorded (apbc_mask, sub_shape, cell_shape) must match.

  * No duplicate fixtures on the same *allowed* (gated) triple. A
    triple is "gated" when is_supported_triple can only accept it by
    matching the curated frozensets in
    allowlist_predicate.py -- i.e. none of that function's trivial
    early-return branches (non-SOC, SubShape=[1,1,1], or SOC+PBC)
    apply. Only gated triples participate in duplicate-inventory
    checking: on the real fixture tree, several *non-gated* fixtures
    legitimately share a (mask, sub_shape, cell_shape) triple (e.g.
    apbc_complexuhf/case_1d_hubbard, uhfk_mvmc_pairproduct/case_apbc
    and .../case_apbc_halffill are three distinct non-SOC fixtures --
    different Ncond/interaction -- all on CellShape=[8,1,1],
    SubShape=[1,1,1], APBC-x). Flagging those as "duplicates" would
    reject real, already-shipped fixture state, so non-gated triples
    are excluded from this check entirely; only the curated allowlist
    slots (one fixture per gated triple, by design) are guarded.

    Within the gated bucket, a triple claimed by 2+ fixtures is
    tolerated -- not a duplicate error -- when either:
      (a) EXACTLY two fixtures share that triple and both have the
          identical basename (the apbc_complexuhf/ <-> uhfk_mvmc_pairproduct/
          mirror convention: same case_* directory name, different
          fixture family root). The mirror convention is strictly a
          2-fixture pattern by design, so 3+ fixtures sharing both the
          triple and the basename do NOT qualify here and must instead
          satisfy (b) or fail closed, or
      (b) every fixture sharing that triple is a key in expected_map
          (explicit inventory declaration via _V36_EXPECTED_FIXTURE_MAP /
          _V37_EXPECTED_FIXTURE_MAP).
    Otherwise it fails closed as a duplicate.

Fixture roots are canonicalized via Path.resolve(strict=True); repo-
relative resolved paths are used as inventory keys so mirror fixtures
under different roots do not collide by basename, and symlink-escape /
case-mismatched aliases are rejected.

"""
from __future__ import annotations

from pathlib import Path

from tools._uhfk_to_mvmc.allowlist_predicate import (
    is_supported_triple, apbc_mask_of, REJECT_MESSAGE,
)
from tools._uhfk_to_mvmc.boundary_input import (
    normalize_boundary_condition_list,
)
from tools._uhfk_to_mvmc.input_loader import load_input_toml


class AllowlistCoverageError(RuntimeError):
    """Raised on any allowlist / inventory violation."""


def _load_fixture_triple(input_toml_path):
    """Parse the (apbc_mask, sub_shape, cell_shape, is_soc_mode, theta)
    tuple out of a case_*/input.toml file.

    Sources every field -- including is_soc_mode -- from
    input_loader.load_input_toml, the exact shallow [mode] +
    [mode.param] merge the production CLI (tools/uhfk_to_mvmc.py) uses
    ([mode.param] wins on key collision; see that loader's docstring).

    Raises AllowlistCoverageError if input.toml cannot be loaded (no
    [mode]/[mode.param] section, TOML parse error) or is missing
    CellShape / BoundaryCondition. A malformed fixture is a real
    problem the coverage checker must surface loudly -- fail the
    whole check, never silently skip the fixture or leak an
    uncontextualized KeyError/ValueError past this function.
    """
    try:
        merged = load_input_toml(str(input_toml_path))
    except (KeyError, ValueError) as exc:
        raise AllowlistCoverageError(
            f"{input_toml_path}: cannot load input.toml: {exc}"
        ) from exc
    if "CellShape" not in merged:
        raise AllowlistCoverageError(
            f"{input_toml_path}: missing CellShape in [mode.param]"
        )
    if "BoundaryCondition" not in merged:
        raise AllowlistCoverageError(
            f"{input_toml_path}: missing BoundaryCondition in "
            f"[mode.param]"
        )
    is_soc_mode = bool(merged.get("enable_spin_orbital", False))
    cell_shape = tuple(int(c) for c in merged["CellShape"])
    sub_shape = tuple(
        int(s) for s in merged.get("SubShape", list(cell_shape))
    )
    theta = normalize_boundary_condition_list(merged["BoundaryCondition"])
    apbc_mask = apbc_mask_of(theta)
    return apbc_mask, sub_shape, cell_shape, is_soc_mode, theta


def _relpath_anchor(root):
    """Return the anchor directory for computing fixture inventory
    keys: the nearest ancestor of the already-resolved `root` that
    contains both tests/ and tools/ directories (the repo root), so
    real usage against tests/validation/ produces keys matching the
    repo-relative expected_map format (e.g.
    "tests/validation/uhfk_mvmc_pairproduct/case_...").

    Falls back to `root` itself when no such ancestor exists -- this
    is expected and correct for standalone fixture roots with no
    repo context, such as pytest's tmp_path in
    tests/test_uhfk_to_mvmc_allowlist_coverage.py.
    """
    ancestor = root
    while ancestor != ancestor.parent:
        if (ancestor / "tests").is_dir() and (ancestor / "tools").is_dir():
            return ancestor
        ancestor = ancestor.parent
    return root


def _is_gated_triple(apbc_mask, sub_shape, is_soc_mode):
    """True iff (is_soc_mode, sub_shape, apbc_mask) can only be
    accepted by is_supported_triple via the curated
    allowlist frozensets -- i.e. none of that function's trivial
    early-return branches apply.

    Mirrors only the early-return *conditions* of
    allowlist_predicate.is_supported_triple (non-SOC / SubShape=1 /
    zero APBC mask); the frozenset policy data itself is never
    duplicated here and stays the single source of truth in
    allowlist_predicate.py.
    """
    if not is_soc_mode:
        return False
    if not any(s != 1 for s in sub_shape):
        return False
    if apbc_mask == (0, 0, 0):
        return False
    return True


def check_fixture_allowlist_coverage(
    fixture_root, expected_map, require_all_expected=True,
):
    """Iterate case_*/input.toml under fixture_root and enforce the
    declared allowlist inventory, duplicate, and fail-loud checks.

    Parameters
    ----------
    fixture_root : path
        Directory to scan for case_*/input.toml. Canonicalized via
        Path.resolve(strict=True) before use.
    expected_map : dict
        {repo-relative fixture path (str): (apbc_mask, sub_shape,
        cell_shape)}. Each key MUST correspond to an actual on-disk
        fixture when require_all_expected=True, and also exempts
        gated-triple collisions from the duplicate check (see module
        docstring).
    require_all_expected : bool
        If True, every key in expected_map MUST be present in the
        discovered fixture inventory (missing -> fail), and its
        recorded triple must match what's on disk.

    Raises
    ------
    AllowlistCoverageError
        On any violation: a malformed case_*/input.toml (fails loud,
        see _load_fixture_triple), rejected by the CLI predicate,
        missing required fixture, expected-triple mismatch, duplicate
        fixtures on the same gated triple, or a case_* directory /
        input.toml file that is a symlink resolving outside the repo
        root.
    """
    root = Path(fixture_root).resolve(strict=True)
    repo_root = _relpath_anchor(root)

    discovered = {}  # relpath -> triple
    by_gated_triple = {}  # gated triple -> [relpath, ...] (sorted order)

    input_tomls = sorted(
        root.rglob("case_*/input.toml"), key=lambda p: str(p)
    )
    for input_toml in input_tomls:
        # Reject a case_* directory that is a symlink resolving outside
        # repo_root, instead of letting Path.relative_to's bare
        # ValueError escape this function uncontextualized.
        try:
            fixture_dir = input_toml.parent.resolve(strict=True)
        except OSError as exc:
            raise AllowlistCoverageError(
                f"{input_toml}: cannot resolve fixture directory: {exc}"
            ) from exc
        try:
            relpath = str(fixture_dir.relative_to(repo_root))
        except ValueError as exc:
            raise AllowlistCoverageError(
                f"{input_toml}: fixture directory symlink escapes "
                f"repo root"
            ) from exc

        # Reject an input.toml that is itself a symlink resolving
        # outside repo_root -- previously silently admitted, since
        # only the parent directory (not the file) was ever resolved.
        try:
            resolved_input_toml = input_toml.resolve(strict=True)
        except OSError as exc:
            raise AllowlistCoverageError(
                f"{input_toml}: cannot resolve input.toml: {exc}"
            ) from exc
        try:
            resolved_input_toml.relative_to(repo_root)
        except ValueError as exc:
            raise AllowlistCoverageError(
                f"{input_toml}: input.toml symlink escapes repo root"
            ) from exc

        try:
            apbc_mask, sub_shape, cell_shape, is_soc_mode, theta = (
                _load_fixture_triple(input_toml)
            )
        except AllowlistCoverageError:
            # _load_fixture_triple already raises AllowlistCoverageError
            # with file context on any malformed fixture; re-raise as-is
            # (fail loud, never skip) rather than let a bare
            # KeyError/ValueError from a future change escape uncaught.
            raise

        if not is_supported_triple(
            theta, sub_shape, cell_shape, is_soc_mode
        ):
            raise AllowlistCoverageError(
                f"{relpath}: not in the supported allowlist. "
                f"{REJECT_MESSAGE}"
            )

        triple = (apbc_mask, sub_shape, cell_shape)
        discovered[relpath] = triple
        if _is_gated_triple(apbc_mask, sub_shape, is_soc_mode):
            by_gated_triple.setdefault(triple, []).append(relpath)

    for triple, relpaths in by_gated_triple.items():
        if len(relpaths) < 2:
            continue
        basenames = {Path(p).name for p in relpaths}
        # Mirror convention is strictly a 2-fixture pattern (one under
        # uhfk_mvmc_pairproduct/, one under apbc_complexuhf/); 3+
        # fixtures sharing both the triple and the basename are NOT a
        # mirror pair and must still be caught below as a duplicate.
        is_mirror_set = len(basenames) == 1 and len(relpaths) == 2
        all_expected = all(p in expected_map for p in relpaths)
        if not (is_mirror_set or all_expected):
            raise AllowlistCoverageError(
                f"duplicate fixtures on allowed triple {triple}: "
                f"{relpaths!r}"
            )

    if require_all_expected:
        for expected_key, expected_triple in expected_map.items():
            if expected_key not in discovered:
                raise AllowlistCoverageError(
                    f"missing required fixture: {expected_key}"
                )
            if discovered[expected_key] != expected_triple:
                raise AllowlistCoverageError(
                    f"{expected_key}: expected triple {expected_triple}, "
                    f"got {discovered[expected_key]}"
                )
