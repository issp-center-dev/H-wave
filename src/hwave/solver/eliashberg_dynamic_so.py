"""Compatibility shell: this module moved to ``hwave.solver.eliashberg.dynamic_so``
(dev/sc-split-plan.md, step 6). It stays for one release so that
``hwave.solver.eliashberg_dynamic_so`` keeps importing; the name resolves to the moved
module itself (not a copy of its attributes), so patching or reading
``hwave.solver.eliashberg_dynamic_so.<name>`` reaches the code that runs.
"""
import sys

from hwave.solver.eliashberg import dynamic_so as _moved

sys.modules[__name__] = _moved
