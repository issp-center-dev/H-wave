# Enforce Pytest CI Design

## Context

The `CI Python 3.9+` workflow runs pytest with `continue-on-error: true`.
Consequently, the visible pull-request checks can remain green after pytest
fails. The preceding `unittest discover` step covers the legacy
`unittest.TestCase` population but does not guarantee execution of the
function-style pytest regressions added by PR #82.

Additionally, `pytest.ini` uses the setup.cfg-style `[tool:pytest]` section
header. Pytest therefore ignores its marker declarations and default options.
`tests/conftest.py` currently compensates for marker registration only.

## Goals

1. Make pytest failures fail the GitHub Actions job on every tested Python
   version.
2. Make `pytest.ini` a valid pytest configuration.
3. Keep the slow L=32 fixture test skipped unless explicitly selected.
4. Apply one clear coverage policy without silently relying on ignored
   configuration.

## Non-goals

- Refactoring the test suite from `unittest` to pytest.
- Running the expensive L=32 fixture regeneration on every pull request.
- Changing numerical tolerances or scientific acceptance criteria.
- Bundling this repository-level CI change into PR #82.

## Workflow Design

This work is tracked in a separate issue and implemented in a separate pull
request based on `develop`.

The `CI Python 3.9+` workflow will:

- retain `unittest discover` for compatibility;
- run pytest without any step- or job-level `continue-on-error`;
- explicitly exclude slow tests in the ordinary PR matrix command;
- propagate pytest's nonzero status to the job.

The repository's authoritative tested matrix is Python 3.9, 3.10, 3.11, and
3.12, matching all three existing test/sample workflows. The existing
workflow and README badge display name remains `CI Python 3.9+` so the current
required-check identity is not invalidated. A workflow comment states that
the certified matrix is currently 3.9--3.12. The broader Poetry
installation constraint `python = "^3.9"` is not treated as a promise that
every future 3.x release has already been CI-certified. Python 3.13 support is
tracked separately: the 2026-07-27 local run exposed two existing
byte-identical golden failures on 3.13, so silently adding it here would mix a
numerical-compatibility change into CI enforcement.

The separate `Run tests` workflow is not treated as a substitute for the
required `CI Python 3.9+` checks; the primary workflow must be independently
correct.

## Pytest Configuration

`pytest.ini` will use the valid `[pytest]` header. Marker declarations move
back to the configuration file and the temporary registration in
`tests/conftest.py` is removed once redundant.

The complete resulting `pytest.ini` is:

```ini
[pytest]
testpaths = tests
python_files = test_*.py
python_classes = Test*
python_functions = test_*
addopts =
    -v
    --tb=short
    --strict-markers
    --disable-warnings
markers =
    unit: Unit tests
    integration: Integration tests
    slow: Slow running tests
    solver: Solver-specific tests
    io: Input/output tests
```

Thus every existing key has been audited: discovery patterns, verbosity,
traceback mode, strict markers, warning display, and all five marker
declarations are retained; the three coverage-report arguments,
`--cov=src/hwave`, and `--cov-fail-under` are removed from global `addopts`.
There is no global marker filter. `tests/conftest.py` contains only the
temporary marker-registration hook, so that file is removed; no unrelated
fixture or hook is deleted.
`--disable-warnings` is intentionally retained to keep the already large
numerical-warning summary out of four CI logs; warning filters are not changed
and targeted local runs can override `addopts`.

The ordinary Python 3.9--3.12 matrix uses the exact command:

```bash
pytest tests/ -m "not slow"
```

Slow exclusion is workflow-local so the explicit developer commands below
can collect slow tests. The existing `unittest discover` step remains
temporarily duplicated by pytest for compatibility; deduplicating it is a
separate cleanup.

A dedicated Python 3.11 coverage job, in addition to the ordinary 3.11 matrix
entry, uses:

```bash
pytest tests/ -m "not slow" \
  --cov=src/hwave --cov-branch \
  --cov-report=term-missing --cov-report=xml \
  --cov-fail-under=80
```

This defines the coverage target (`src/hwave`), branch coverage, reports, and
single enforcement environment. `pytest-cov` is already a development
dependency. A 2026-07-27 baseline on the available Python 3.13 environment
measured 85% statement coverage and 83% combined statement/branch coverage
with the exact selection. The Python 3.11 job must independently measure at
least 80% before the CI pull request is opened. If it does not, tests are
added; the threshold and `--cov-branch` are not silently weakened. The
generated `coverage.xml` is uploaded through the existing Codecov step as a
review aid, but Codecov transport remains non-blocking; the local
`--cov-fail-under=80` result is the enforcement.

The workflow keeps `name: CI Python 3.9+`. The matrix job gains the stable
display name `Tests (Python ${{ matrix.python-version }})`, yielding the four
intended checks `CI Python 3.9+ / Tests (Python 3.9)` through
`... (Python 3.12)`. The coverage job is named `Coverage (Python 3.11)`.
All five must be configured as required branch-protection checks. Because
branch protection is external state, the issue records the current and
intended check names and the PR
cannot claim merge protection is complete until a maintainer verifies that
configuration. Both jobs install `poetry install --no-interaction --with dev`,
which supplies pytest and pytest-cov, followed by the existing optional
`sparse-ir` install.

## Slow-test Selection

The marker classifies test cost; the environment variable is the only
execution authorization. `_require_slow_fixtures` is simplified to inspect
only `HWAVE_RUN_SLOW_FIXTURES`, not pytest's marker-expression string:

| Selection | Environment | Result |
|---|---|---|
| `-m "not slow"` | unset or `1` | deselected |
| plain pytest / direct node | unset | collected, then skipped |
| plain pytest / direct node | `1` | executed |
| `-m slow` or a compound positive expression | unset | selected, then skipped |
| `-m slow` or a compound positive expression | `1` | executed |

No new global hook is added. Developers run the L=32 fixture test with:

```bash
HWAVE_RUN_SLOW_FIXTURES=1 pytest -m slow
```

The environment variable is required; `pytest -m slow` alone only selects the
test and it remains skipped. A direct-node command with the same environment
variable is equivalent.

## Testing

Before publication:

1. On the CI-supported versions, run the exact pytest command used by CI and
   confirm zero failures. On Python 3.11, also run the exact coverage command
   and confirm combined branch coverage is at least 80%.
2. Run `python -m unittest discover tests/ -q`.
3. Run `pytest --collect-only -m slow` and
   `pytest --collect-only -m "not slow"` and assert the exact L=32 node ID is
   present only in the first output.
4. Confirm `pytest --markers` lists every project marker exactly once.
5. Confirm a deliberately failing temporary pytest test makes the exact
   matrix command return nonzero; the temporary test is never committed.
6. Add a repository test that parses every workflow YAML and rejects
   `continue-on-error: true` on a step whose command invokes pytest or
   unittest, or on its containing job. Validate YAML with `actionlint` where
   available and assert each test step's `run` value equals its declared exact
   command, with no pipeline, `||`, wrapper, or trailing status overwrite.
   Add PyYAML 6 as an explicit development dependency and use
   `yaml.BaseLoader`, which preserves GitHub Actions' `on` key as text instead
   of applying YAML 1.1 boolean coercion.
7. With `HWAVE_RUN_SLOW_FIXTURES` absent, confirm plain pytest reports the
   exact L=32 node as skipped and the matrix command reports it as deselected;
   use a lightweight monkeypatched guard test to prove both positive opt-ins
   without regenerating L=32 data.
8. Run `git diff --check`.
9. From the repository root, run `pytest --trace-config --collect-only` and
   assert the active config path is the repository's `pytest.ini`.

The first pull-request run supplies workflow-level positive evidence across
all tested Python versions. A deliberately failing commit is not pushed;
the static workflow assertions plus local nonzero-exit test establish the
negative path without polluting shared history.

Branch-protection migration is ordered: first push the workflow and let all
five intended check names complete once; then an administrator adds those
exact names as required while retaining any old required identity; finally
verify a PR is blocked by all five before removing an obsolete required name
and before merging.

## Delivery

The issue explains that this is a pre-existing repository-level merge-safety
problem surfaced by PR #82. The pull request remains separate so the
scientific implementation diff and CI-policy diff can be reviewed and reverted
independently.
