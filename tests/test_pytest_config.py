"""Detect the pytest.ini header regression (issue #84).

pytest.ini carried the [tool:pytest] header -- valid only in setup.cfg --
for a long time, and pytest silently ignored the whole file: collection
still succeeded, so nothing noticed. This test fails if that ever
happens again, by asserting the configuration the file declares was
actually APPLIED, not merely that collection worked.

pytest-style on purpose: it needs the pytest config object, and since
the CI pytest step gates (issue #84 follow-through), a failure here
fails CI even though unittest discovery cannot see this module.
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


def test_coverage_gate_boundary_semantics():
    """The CI coverage floor compares the ROUNDED total (PR #126 round-3
    review measured 79.50% passing an 80 floor at coverage's default
    zero-decimal precision). Pin the boundary at the two-decimal
    precision .coveragerc sets, and that the config is actually loaded."""
    try:
        from coverage import Coverage
        from coverage.results import should_fail_under
    except ImportError:
        import pytest
        pytest.skip("coverage not installed (bare-pytest environments)")

    assert should_fail_under(79.99, 80, 2)
    assert not should_fail_under(80.00, 80, 2)
    # at the DEFAULT precision the leak this guards against reappears
    assert not should_fail_under(79.50, 80, 0)

    cov = Coverage()
    assert cov.config.precision == 2, (
        ".coveragerc precision=2 was not loaded (cwd must be the repo "
        "root; got precision={})".format(cov.config.precision))
