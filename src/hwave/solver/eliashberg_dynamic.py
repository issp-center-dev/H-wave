"""Compatibility shell: this module moved to ``hwave.solver.eliashberg.dynamic``
(dev/sc-split-plan.md, step 6). It stays for one release so that
``hwave.solver.eliashberg_dynamic`` keeps importing; the name resolves to the moved
module itself (not a copy of its attributes), so patching or reading
``hwave.solver.eliashberg_dynamic.<name>`` reaches the code that runs.
"""
import sys

from hwave.solver.eliashberg import dynamic as _moved

sys.modules[__name__] = _moved
