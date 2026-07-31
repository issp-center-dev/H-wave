"""Detect the pytest.ini header regression (issue #84).

pytest.ini carried the [tool:pytest] header -- valid only in setup.cfg --
for a long time, and pytest silently ignored the whole file: collection
still succeeded, so nothing noticed. This test fails if that ever
happens again, by asserting the configuration the file declares was
actually APPLIED, not merely that collection worked.

pytest-style on purpose: it needs the pytest config object. unittest
discovery imports this module but collects no function-style tests;
since the CI pytest step gates (issue #84 follow-through), a failure
here still fails CI.
"""


def test_pytest_ini_is_actually_loaded(request):
    config = request.config

    # the config file pytest resolved must be our pytest.ini
    assert config.inipath is not None, (
        "pytest resolved NO ini file: pytest.ini is missing or unreadable")
    assert config.inipath == config.rootpath / "pytest.ini", (
        "resolved ini is not the repository's own pytest.ini: {}".format(
            config.inipath))

    # the marker registry from the ini must be present (a wrong header
    # silently drops it; --strict-markers would then fail elsewhere, but
    # only if a marked test is collected -- assert it directly)
    markers = config.getini("markers")
    assert any(str(m).startswith("slow:") for m in markers), (
        "pytest.ini's marker registry was not applied -- most likely the "
        "section header regressed from [pytest] to [tool:pytest], which "
        "pytest ignores in a standalone pytest.ini (issue #84): {}".format(
            markers))

    # addopts must have been applied too
    assert config.getoption("--strict-markers"), (
        "pytest.ini addopts were not applied (issue #84 header regression)")


def test_coverage_floor_is_exact():
    """The CI floor is enforced by tests/_coverage_floor.py on the
    integer line counts from coverage.xml -- NOT by the rounded
    --cov-fail-under comparison, which at any finite precision passes a
    band just below the floor (79.50% passed at zero decimals; 79.995%
    would at two, and the round-5 review showed 6151/7689 = 79.9974% is
    three added statements away in this tree). Pin the exact predicate
    on those boundary cases."""
    from tests._coverage_floor import passes

    assert passes(4, 5, 80)             # exactly 80.00%
    assert passes(80, 100, 80)
    assert not passes(7999, 10000, 80)  # 79.99%
    assert not passes(6151, 7689, 80)   # 79.9974% -- rounds to 80.00
    assert not passes(0, 0, 80)         # nothing measured is not a pass


def test_coverage_report_precision_is_two_decimals():
    """.coveragerc sets two-decimal report precision (display and the
    soft --cov-fail-under comparison). Asserted through coverage's
    PUBLIC get_option API only -- a round-5 review found the earlier
    version of this test bound to undocumented internals, which a
    coverage upgrade could rename (silently skipping the test under the
    broad ImportError guard)."""
    try:
        from coverage import Coverage
    except ImportError:
        import pytest
        pytest.skip("coverage not installed (bare-pytest environments)")

    cov = Coverage(data_file=None)
    assert cov.get_option("report:precision") == 2, (
        ".coveragerc precision=2 was not loaded (cwd must be the repo "
        "root)")


def test_coverage_floor_cli(tmp_path):
    """CLI-level pin of the exact-floor checker: passing, failing,
    malformed and missing XML must produce the right exit codes, and
    the failure modes must be loud (nonzero), never a silent pass."""
    import subprocess
    import sys

    script = "tests/_coverage_floor.py"

    def run(*args):
        return subprocess.run([sys.executable, script, *args],
                              capture_output=True, text=True)

    def xml(covered, valid):
        p = tmp_path / "cov_{}_{}.xml".format(covered, valid)
        p.write_text('<coverage lines-covered="{}" lines-valid="{}">'
                     '</coverage>'.format(covered, valid))
        return str(p)

    assert run(xml(4, 5), "80").returncode == 0
    assert run(xml(7999, 10000), "80").returncode == 1
    assert run(xml(6151, 7689), "80").returncode == 1

    bad = tmp_path / "bad.xml"
    bad.write_text("<coverage")
    assert run(str(bad), "80").returncode != 0
    assert run(str(tmp_path / "missing.xml"), "80").returncode != 0
