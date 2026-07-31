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
    assert config.inipath.name == "pytest.ini", str(config.inipath)

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
