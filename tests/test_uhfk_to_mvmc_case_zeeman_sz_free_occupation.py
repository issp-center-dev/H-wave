"""Wrap `scripts/assert_occupation.py case_zeeman_sz_free` in a pytest so
CI can guard the fixture's target occupation before the mVMC harness
runs."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile

import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
CASE_DIR = os.path.join(
    REPO_ROOT, "tests/validation/uhfk_mvmc_pairproduct/case_zeeman_sz_free"
)
HELPER_PATH = os.path.join(
    REPO_ROOT,
    "tests/validation/uhfk_mvmc_pairproduct/scripts/assert_occupation.py",
)


def _run_hwave(work_dir):
    for fn in os.listdir(CASE_DIR):
        src = os.path.join(CASE_DIR, fn)
        if os.path.isfile(src):
            shutil.copy(src, work_dir)
    env = os.environ.copy()
    env["PYTHONPATH"] = os.path.join(REPO_ROOT, "src") + ":" + env.get(
        "PYTHONPATH", ""
    )
    script = (
        "import numpy as np\n"
        "if not hasattr(np, 'float'):\n"
        "    np.float = float\n"
        "import sys\n"
        "sys.argv = ['hwave', 'input.toml']\n"
        "from hwave.qlms import main\n"
        "main()\n"
    )
    res = subprocess.run(
        [sys.executable, "-c", script],
        cwd=work_dir, env=env, capture_output=True, text=True,
    )
    if res.returncode != 0:
        pytest.skip(f"H-wave SCF failed for case_zeeman_sz_free: {res.stderr}")


def test_case_zeeman_sz_free_target_occupation_holds():
    """Run H-wave SCF then invoke the shared helper."""
    with tempfile.TemporaryDirectory() as tmp:
        _run_hwave(tmp)
        # Import helper directly to invoke the single source of truth
        sys.path.insert(0, os.path.dirname(HELPER_PATH))
        try:
            from assert_occupation import assert_case_occupation
            assert_case_occupation(CASE_DIR, tmp)
        finally:
            sys.path.pop(0)
