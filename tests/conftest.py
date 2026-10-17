"""Shared pytest configuration for the H-wave test suite.

Marker registration only. ``pytest.ini`` declares the project's markers under a
``[tool:pytest]`` header, which is the ``setup.cfg`` section name -- pytest does
not read it out of a ``pytest.ini`` file, so the declarations there never reach
the runtime and any use of a marker raises ``PytestUnknownMarkWarning``.
Registering here keeps the markers honest without changing how the suite is
invoked (fixing the ini header would also switch on ``addopts``, i.e. coverage
gating, which is a separate decision).
"""


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "slow: long-running test, skipped by default. Select with "
        "`pytest -m slow`. Currently used by "
        "tests/test_tsweep_e2e.py and tests/test_eliashberg_dynamic.py. "
        "tests/test_bond_onari_milestone.py::TestBondOnariMilestone::"
        "test_grid_convergence_16_to_32 (which regenerates the uncommitted "
        "L=32 FLEX green fixtures) is a separate, unittest-only case: it does "
        "NOT carry this marker (that module is unittest-discoverable, so "
        "pytest markers do not gate it under `python -m unittest`) and is "
        "opted into solely via HWAVE_RUN_SLOW_FIXTURES=1, e.g.:\n"
        "    HWAVE_RUN_SLOW_FIXTURES=1 python -m unittest "
        "tests.test_bond_onari_milestone.TestBondOnariMilestone."
        "test_grid_convergence_16_to_32")
