#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
CASE="${1:?usage: run.sh <case-dir>}"

UHF="${HERE}/build/mvmc/build/src/ComplexUHF/UHF"
VMCDRY="${HERE}/build/mvmc/build/src/mVMC/vmcdry.out"
if [[ ! -x "${UHF}" || ! -x "${VMCDRY}" ]]; then
  echo "Binaries missing. Run build_complexuhf.sh first." >&2
  exit 1
fi

WORK="${HERE}/${CASE}/work"
rm -rf "${WORK}"
mkdir -p "${WORK}/hwave" "${WORK}/complexuhf"

# ---- H-wave ----
cp "${HERE}/${CASE}/input.toml"        "${WORK}/hwave/"
cp "${HERE}/${CASE}/Geometry.dat"      "${WORK}/hwave/"
cp "${HERE}/${CASE}/Transfer.dat"      "${WORK}/hwave/"
cp "${HERE}/${CASE}/CoulombIntra.dat"  "${WORK}/hwave/"
cp "${HERE}/${CASE}/geometry_uhf.dat"  "${WORK}/hwave/"
cp "${HERE}/${CASE}/OneBodyG.dat"      "${WORK}/hwave/"

(
  cd "${WORK}/hwave"
  # hwave.qlms has no __main__ block; invoke the CLI entry point directly.
  # The numpy>=1.20 shim restores the legacy ``np.float`` alias that
  # qlmsio.wan90.read_geometry still uses (3 call sites). The shim only
  # affects this child process and does not touch H-wave source.
  PYTHONPATH=/workspace/H-wave-dev/src python3 -c "
import numpy as np
if not hasattr(np, 'float'):
    np.float = float  # numpy>=1.20 compat shim for wan90.read_geometry
import sys
sys.argv = ['hwave', 'input.toml']
from hwave.qlms import main
main()
"
)

# ---- ComplexUHF ----
cp "${HERE}/${CASE}/stan.in" "${WORK}/complexuhf/"
(
  cd "${WORK}/complexuhf"
  "${VMCDRY}" stan.in
  # APBC is encoded in stan.in via `phase0 = 180.0`, which makes StdFace
  # emit Trans.def with -1 on the wrap-around bond. The override below sets
  # NMPTrans = -1 (ComplexUHF's APFlag bit) only as a legacy marker; see
  # complexuhf_modpara_override.txt for why this alone is insufficient.
  bash "${HERE}/${CASE}/complexuhf_modpara_override.txt"
  "${UHF}" namelist.def
)

# ---- Compare ----
python3 "${HERE}/compare.py" \
  --hwave      "${WORK}/hwave/output" \
  --complexuhf "${WORK}/complexuhf" \
  --case       "${CASE}"
