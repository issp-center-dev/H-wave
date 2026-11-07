"""Compatibility shell: this module moved to ``hwave.solver.eliashberg.bond_dynamic``
(dev/sc-split-plan.md, step 6). It stays for one release so that
``hwave.solver.eliashberg_bond`` keeps importing; the name resolves to the moved
module itself (not a copy of its attributes), so patching or reading
``hwave.solver.eliashberg_bond.<name>`` reaches the code that runs.
"""
import sys

from hwave.solver.eliashberg import bond_dynamic as _moved

sys.modules[__name__] = _moved
