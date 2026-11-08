"""The shard map of the FULL CI job (``test-full`` in
``.github/workflows/ci-python39.yml``).

``test-full`` runs the registered heavy tests (``tests/heavy_tests.py``) as
:data:`N_SHARDS` parallel jobs. Each job takes the registered tests of the
FILES listed for it below -- file granularity, so the module- and
class-level caches of a test file (and the shared ED cache,
``tests/ed_cache.py``, between files of one job) stay inside one process --
plus, in shard 1, the registry's own consistency gate.

The split is a longest-processing-time-first (LPT) partition of the
measured per-file time of the registered tests (:func:`lpt_partition`):
files in decreasing time, each to the currently lightest shard. MEASURED
input, see :data:`MEASUREMENT`; re-measure and re-partition when a file's
cost moves a lot, with::

    python -m tests.heavy_shards plan JUNIT.xml

on the JUnit report of a full run of the registered tests.

``tests/test_heavy_tests_registry.py`` fails the FAST gate when this map
and the registry disagree: a file holding a registered test that is in no
shard, or in two, a shard naming a file that holds none, or a shard count
that differs from the workflow's matrix.
"""

#: Number of ``test-full`` jobs. Must match the ``shard`` matrix of the
#: workflow (checked by the registry test).
N_SHARDS = 3

#: Where the per-file times behind :data:`SHARDS` come from.
MEASUREMENT = (
    "summed test time of the registered tests per file, projected for the "
    "CI runner: CI run 42418446187 (2026-11-07, Python 3.12, numpy 1.26.4) "
    "times the after/before ratio of the ED speed-ups of 2026-11-08 "
    "measured on a 32-core Linux host (numpy 1.26.4, "
    "OPENBLAS_NUM_THREADS=4); projected shard loads 460 / 461 / 465 s")

#: Files that must share a shard: they diagonalize identical Hamiltonians,
#: which ``tests/ed_cache.py`` shares only within one process (the case-M
#: off-site CoulombInter ray). :func:`lpt_partition` callers treat each
#: group as one item.
TOGETHER = (
    ("test_bond_longitudinal_ed", "test_bond_vs_ed_oracle"),
)

#: ``SHARDS[i]`` = the bare module names (as in ``HEAVY_TESTS``) of shard
#: ``i + 1``, with their projected seconds (see :data:`MEASUREMENT`).
SHARDS = (
    (
        "test_so_pairing_ed",                      # 193
        "test_offsite_exchange_ed_longitudinal",   # 96
        "test_eliashberg_dynamic_so",              # 79
        "test_bond_pairing_ed",                    # 37
        "test_flex_second_order_sopt",             # 32
        "test_flex_sigma_hermitian",               # 13
        "test_bond_supercell",                     # 11
        "test_flex_bond_onari_trend",              # 0 (skips without its fixtures)
        "test_flex_bond_gpu",                      # 0 (skips without a CUDA device)
    ),
    (
        "test_spinful_transverse_ed",              # 188
        "test_flex_spin_orbital_ed",               # 111
        "test_flex_second_order_ed_chain",         # 59
        "test_so_bond_chi_ed",                     # 53
        "test_flex_bond_device_store",             # 18
        "test_uhfk_orientation",                   # 16
        "test_flex_second_order_ed_onsite",        # 13
        "test_flex_bond_memory",                   # 2
    ),
    (
        "test_bond_longitudinal_ed",               # 109 } one item (TOGETHER)
        "test_bond_vs_ed_oracle",                  # 61  }
        "test_so_bond_pairing_ed",                 # 135
        "test_so_bond_conventions",                # 56
        "test_bond_transverse_ed",                 # 44
        "test_bond_half_grid",                     # 30
        "test_rpa_vs_ed_oracle",                   # 15
        "test_so_bond_pairing_oracle",             # 8
        "test_eliashberg_bond_entry",              # 6
    ),
)


def shard_of_module():
    """``{module: shard number (1-based)}``; raises ``ValueError`` when a
    module is listed twice."""
    out = {}
    for i, mods in enumerate(SHARDS, 1):
        for mod in mods:
            if mod in out:
                raise ValueError("{} is in shards {} and {}".format(mod, out[mod], i))
            out[mod] = i
    return out


def lpt_partition(weights, n):
    """Longest-processing-time-first partition of ``{item: weight}`` into
    ``n`` bins: ``(bins, loads)``, items sorted by decreasing weight (ties
    by name, so the result is deterministic), each placed in the currently
    lightest bin (ties: the lowest index)."""
    bins = [[] for _ in range(n)]
    loads = [0.0] * n
    for item, w in sorted(weights.items(), key=lambda kv: (-kv[1], kv[0])):
        i = loads.index(min(loads))
        bins[i].append(item)
        loads[i] += w
    return bins, loads


def _module_times_from_junit(path):
    """``{module: seconds}`` summed over the testcases of a JUnit report."""
    import xml.etree.ElementTree as ET

    out = {}
    for case in ET.parse(path).getroot().iter("testcase"):
        mod = case.get("classname", "").split(".")
        mod = mod[1] if mod and mod[0] == "tests" and len(mod) > 1 else mod[0]
        out[mod] = out.get(mod, 0.0) + float(case.get("time", 0.0))
    return out


def _main(argv=None):
    import argparse

    from tests import heavy_tests

    parser = argparse.ArgumentParser(prog="python -m tests.heavy_shards")
    sub = parser.add_subparsers(dest="cmd", required=True)
    plan = sub.add_parser("plan", help="LPT partition of a measured run")
    plan.add_argument("junit_xml")
    args = parser.parse_args(argv)

    registered = {mod for mod, _c, _m, _r in heavy_tests.HEAVY_TESTS}
    times = _module_times_from_junit(args.junit_xml)
    weights = {mod: times.get(mod, 0.0) for mod in registered}
    # a TOGETHER group is one item of the partition
    items = {(mod,): w for mod, w in weights.items()
             if not any(mod in grp for grp in TOGETHER)}
    for grp in TOGETHER:
        items[tuple(grp)] = sum(weights.get(mod, 0.0) for mod in grp)
    bins, loads = lpt_partition(items, N_SHARDS)
    for i, (groups, load) in enumerate(zip(bins, loads), 1):
        print("shard {}: {:.0f} s".format(i, load))
        for grp in groups:
            for mod in grp:
                print("    {!r},  # {:.1f} s".format(mod, weights[mod]))
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
