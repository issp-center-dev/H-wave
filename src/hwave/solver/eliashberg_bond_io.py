"""Compatibility shell: this module moved to ``hwave.solver.eliashberg.bond_io``
(dev/sc-split-plan.md, step 6). It stays for one release so that
``hwave.solver.eliashberg_bond_io`` keeps importing; the name resolves to the moved
module itself (not a copy of its attributes), so patching or reading
``hwave.solver.eliashberg_bond_io.<name>`` reaches the code that runs.
"""
import sys

from hwave.solver.eliashberg import bond_io as _moved

sys.modules[__name__] = _moved
