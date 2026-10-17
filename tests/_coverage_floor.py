"""Exact coverage-floor check for CI (issue #84 / PR #126).

``--cov-fail-under`` compares the floor against the ROUNDED total, so at
any finite precision some band below the floor still passes (at the
default zero decimals, 79.50% passed an 80 floor; at two decimals,
79.995% still would -- and the review showed that value is attainable
in this tree). This checker instead compares the integer line counts
from coverage's XML report exactly: covered/valid >= floor/100 with no
floating point and no rounding.

CI invokes it right after the pytest+coverage step:

    python tests/_coverage_floor.py coverage.xml 80

The pure predicate is unit-tested in test_pytest_config.py. The file is
named with a leading underscore so neither pytest (test_*.py) nor
unittest discovery (test*.py) collects it as a test module.
"""

import sys
import xml.etree.ElementTree as ET


def passes(covered, valid, floor_percent):
    """True iff covered/valid >= floor_percent/100, exactly.

    Integer arithmetic only: covered * 100 >= valid * floor_percent is
    the cross-multiplied form, exact for any int inputs. A report with
    no valid lines fails (nothing was measured -- that is never a pass).
    """
    covered = int(covered)
    valid = int(valid)
    floor_percent = int(floor_percent)
    if valid <= 0:
        return False
    return covered * 100 >= valid * floor_percent


def main(argv):
    xml_path, floor_percent = argv[1], int(argv[2])
    root = ET.parse(xml_path).getroot()
    covered = int(root.get("lines-covered"))
    valid = int(root.get("lines-valid"))
    ok = passes(covered, valid, floor_percent)
    pct = 100.0 * covered / valid if valid else 0.0
    print("exact coverage floor: {}/{} lines = {:.4f}% vs {}%: {}".format(
        covered, valid, pct, floor_percent, "OK" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
