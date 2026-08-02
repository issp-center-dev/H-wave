#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="${HERE}/../../../../mVMC-1.4.0"
BUILD="${HERE}/build/mvmc"

if [[ ! -d "${SRC}" ]]; then
  echo "mVMC-1.4.0 not found at ${SRC}" >&2
  exit 1
fi

rm -rf "${BUILD}"
mkdir -p "${BUILD}"
cp -r "${SRC}"/* "${BUILD}/"
rm -rf "${BUILD}/build"   # remove stale CMake cache from host
cd "${BUILD}" && mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release > "${HERE}/build/cmake.log" 2>&1
make UHF vmcdry.out vmc.out -j2 > "${HERE}/build/make.log" 2>&1
echo "Built UHF:     ${BUILD}/build/src/ComplexUHF/UHF"
echo "Built vmcdry:  ${BUILD}/build/src/mVMC/vmcdry.out"
echo "Built vmc.out: ${BUILD}/build/src/mVMC/vmc.out"
