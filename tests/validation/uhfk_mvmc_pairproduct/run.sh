#!/usr/bin/env bash
# Orchestrate the H-wave UHFk -> bridge -> mVMC PairProduct E2E for one case.
# Usage: ./run.sh case_pbc   (or case_apbc)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="${HERE}/../../.."

_validate_transfer_hermiticity() {
  local transfer_path="$1"
  PYTHONPATH="${ROOT}/src:${ROOT}" python3 - "${transfer_path}" <<'PYEOF'
import sys

from tools._uhfk_to_mvmc.transfer_hermiticity import (
    check_transfer_dat_hermiticity,
    TransferHermiticityError,
)

path = sys.argv[1]
try:
    check_transfer_dat_hermiticity(path)
except (FileNotFoundError, TransferHermiticityError) as exc:
    print(f"run.sh: Transfer.dat Hermiticity check failed: {exc}", file=sys.stderr)
    raise SystemExit(1) from exc
print(f"Transfer.dat Hermiticity check: PASS ({path})")
PYEOF
}

if [[ "${1:-}" == "--all-v37" ]]; then
  set -eux
  for name in xy xz yz xyz; do
    _validate_transfer_hermiticity \
      "${HERE}/case_soc_rashba_3d_sub_apbc_${name}/Transfer.dat"
  done
  for name in xy xz yz xyz; do
    bash "$0" case_soc_rashba_3d_sub_apbc_${name}
  done
  echo "--all-v37: SEVEN-GATE PASS x 4 fixtures"
  exit 0
fi

# vmc.out / UHF are built against conda-provided libopenblas / libmpi that
# usually live under a Python env's lib dir (e.g. miniconda3/envs/py39/lib).
# The system loader would otherwise fail with "libopenblas.so.0 not found"
# on machines where the shared libs are not on the default ld path.
#
# Codex adversarial-review 2026-07-12 hardening: MVMC_LD_LIBRARY_PATH is
# not accepted as-is. If set, it MUST be:
#   1. A single absolute path (no colons).
#   2. An existing directory.
#   3. Owned by the current user AND not world-writable.
# Otherwise run.sh drops it and falls back to the auto-detected conda
# path (or nothing, in which case vmc.out is expected to fail at load
# time and the harness reports that failure).
_validate_mvmc_ld_path() {
  local p="$1"
  [[ -n "$p" ]] || return 1
  # No colon-separated lists (would allow prepending malicious dirs).
  if [[ "$p" == *:* ]]; then
    echo "MVMC_LD_LIBRARY_PATH must be a single directory; got ${p!q}" >&2
    return 1
  fi
  # Absolute path only.
  if [[ "$p" != /* ]]; then
    echo "MVMC_LD_LIBRARY_PATH must be absolute; got ${p!q}" >&2
    return 1
  fi
  # Directory exists.
  if [[ ! -d "$p" ]]; then
    echo "MVMC_LD_LIBRARY_PATH is not a directory: ${p!q}" >&2
    return 1
  fi
  # Owner check: must be the current user.
  local owner
  owner="$(stat -c '%u' "$p" 2>/dev/null)"
  if [[ "$owner" != "$(id -u)" ]]; then
    echo "MVMC_LD_LIBRARY_PATH not owned by current user (uid=${owner}): ${p!q}" >&2
    return 1
  fi
  # World-writable check.
  local perm
  perm="$(stat -c '%a' "$p" 2>/dev/null)"
  if [[ "${perm: -1}" =~ [2367] ]]; then
    echo "MVMC_LD_LIBRARY_PATH is world-writable (perm=${perm}): ${p!q}" >&2
    return 1
  fi
  return 0
}

MVMC_LD_LIBRARY_PATH_ACCEPTED=""
if [[ -n "${MVMC_LD_LIBRARY_PATH:-}" ]]; then
  if _validate_mvmc_ld_path "${MVMC_LD_LIBRARY_PATH}"; then
    MVMC_LD_LIBRARY_PATH_ACCEPTED="${MVMC_LD_LIBRARY_PATH}"
  else
    echo "run.sh: dropping unsafe MVMC_LD_LIBRARY_PATH; using auto-detect" >&2
  fi
fi
if [[ -z "${MVMC_LD_LIBRARY_PATH_ACCEPTED}" ]]; then
  # Auto-detect a trusted conda env under the current user's home.
  _autodetect="$(readlink -f "$HOME/miniconda3/envs/py39/lib" 2>/dev/null || true)"
  if [[ -n "${_autodetect}" ]] && _validate_mvmc_ld_path "${_autodetect}" \
      && [[ -f "${_autodetect}/libopenblas.so.0" ]]; then
    MVMC_LD_LIBRARY_PATH_ACCEPTED="${_autodetect}"
  fi
fi
if [[ -n "${MVMC_LD_LIBRARY_PATH_ACCEPTED}" ]]; then
  # Codex re-review 2026-07-12 finding 3 fix: CLEAR the caller-supplied
  # LD_LIBRARY_PATH before prepending our validated path. Otherwise the
  # attacker sets LD_LIBRARY_PATH=/malicious in their env and run.sh
  # appends /malicious unchanged to the accepted path. The vmc.out /
  # UHF child processes then load shared libs from /malicious first.
  # After the reset LD_LIBRARY_PATH contains ONLY the validated
  # MVMC_LD_LIBRARY_PATH_ACCEPTED entry; nothing inherited from the
  # caller's environment survives.
  export LD_LIBRARY_PATH="${MVMC_LD_LIBRARY_PATH_ACCEPTED}"
else
  # No validated path -> caller's LD_LIBRARY_PATH still poses a risk;
  # unset it for the child processes so unvalidated inherited dirs
  # cannot be searched.
  unset LD_LIBRARY_PATH
fi
# Codex 3rd review 2026-07-12 finding: LD_LIBRARY_PATH sanitization
# alone leaves LD_PRELOAD and LD_AUDIT as inherited injection paths.
# glibc treats both as separate loader-time hooks that can execute
# attacker-supplied code inside vmcdry.out / vmc.out / UHF regardless
# of LD_LIBRARY_PATH. Drop them here alongside LD_LIBRARY_PATH so the
# native solver environment is clean of loader-injection controls.
# (Other LD_* vars like LD_BIND_NOW change behavior but do not inject
# code; they are left alone. If a future glibc adds a new injection
# vector, extend this list.)
unset LD_PRELOAD
unset LD_AUDIT
# Codex 4th review 2026-07-12 finding: a caller-controlled BASH_ENV +
# DEBUG trap can re-set LD_PRELOAD / LD_AUDIT between our unset above
# and the child-process exec, defeating the top-level scrub. Also
# reject BASH_ENV / ENV so no attacker-supplied script auto-sources
# on shell startup. Every native solver invocation is additionally
# wrapped in `run_native` (defined below) which uses `env -u` to
# strip these vars per-command; the top-level unsets are kept as a
# belt-and-braces defense so subshells that forget to call
# run_native still start clean.
unset BASH_ENV
unset ENV
# Neutralize any inherited traps (DEBUG, ERR, RETURN, etc.) that could
# re-poison the loader env between commands.
trap - DEBUG ERR RETURN EXIT

# Per-command sanitizer for all native solver executions. Ensures that
# no in-shell trap or startup hook can re-inject LD_PRELOAD / LD_AUDIT
# between our top-level unsets and the actual exec: env -u LD_PRELOAD
# -u LD_AUDIT explicitly REMOVES them at exec time, and LD_LIBRARY_PATH
# is passed only if it was validated above (otherwise -u LD_LIBRARY_PATH
# explicitly removes it too).
#
# Codex 5th review 2026-07-12 finding: an unqualified `env` can be
# shadowed by a caller-exported bash function (BASH_FUNC_env%%=...) or
# by PATH poisoning to a fake env binary. Both defeat the sanitize.
# Fix: call env via the absolute /usr/bin/env path, and additionally
# unset the exported function form (`unset -f env` catches a function
# shadow that already made it into this shell's scope). Fail closed
# if /usr/bin/env is missing.
if [[ ! -x "/usr/bin/env" ]]; then
  echo "run.sh: /usr/bin/env not executable; refusing to run without" >&2
  echo "  a trusted env(1) for per-command loader sanitization." >&2
  exit 1
fi
readonly _RUN_NATIVE_ENV="/usr/bin/env"
unset -f env 2>/dev/null || true
run_native() {
  if [[ -n "${LD_LIBRARY_PATH:-}" ]]; then
    "${_RUN_NATIVE_ENV}" -u LD_PRELOAD -u LD_AUDIT \
        LD_LIBRARY_PATH="${LD_LIBRARY_PATH}" "$@"
  else
    "${_RUN_NATIVE_ENV}" -u LD_PRELOAD -u LD_AUDIT \
        -u LD_LIBRARY_PATH "$@"
  fi
}

CASE="${1:-case_pbc}"

# Codex v3.4 Rev.4 finding: reject path-traversing CASE values so a
# `../whatever` or absolute path cannot cause the later `rm -rf` on
# ${CASE_DIR}/work to escape the fixture directory and delete
# unrelated work dirs. Only accept plain fixture directory names
# (no `/`, no `..`, no leading `.`) that resolve inside ${HERE}.
if [[ "${CASE}" =~ [/] ]] || [[ "${CASE}" == ".."* ]] || [[ "${CASE}" == "."* ]]; then
    echo "ERROR: CASE argument '${CASE}' must be a plain fixture name" >&2
    echo "(no '/', no leading '.'; use e.g. case_pbc, case_soc_rashba_2d_nosub)." >&2
    exit 1
fi
CASE_DIR="${HERE}/${CASE}"
[[ -d "${CASE_DIR}" ]] || { echo "case not found: ${CASE_DIR}" >&2; exit 1; }
# Codex Rev.5 finding: reject symlink CASE directories so
# `case_alias -> case_pbc` cannot cause rm -rf on case_pbc's
# work dir when the user typed `case_alias`.
if [[ -L "${CASE_DIR}" ]]; then
    echo "ERROR: CASE_DIR ${CASE_DIR} is a symlink; refuse to risk" >&2
    echo "cross-fixture rm -rf. Use the target fixture name directly." >&2
    exit 1
fi
# Belt-and-braces: after directory resolution, confirm CASE_DIR is
# actually a direct child of HERE AND matches the requested name
# byte-for-byte (guards against nested symlink shenanigans).
HERE_REAL="$(cd "${HERE}" && pwd -P)"
CASE_DIR_REAL="$(cd "${CASE_DIR}" && pwd -P)"
if [[ "${CASE_DIR_REAL}" != "${HERE_REAL}/${CASE}" ]]; then
    echo "ERROR: resolved CASE_DIR ${CASE_DIR_REAL} does not match" >&2
    echo "expected ${HERE_REAL}/${CASE}; refuse cross-fixture operation." >&2
    exit 1
fi
# From this point on, use the resolved canonical path for all
# destructive operations.
CASE_DIR="${CASE_DIR_REAL}"

_validate_transfer_hermiticity "${CASE_DIR}/Transfer.dat"

# We reuse the mVMC build under apbc_complexuhf/build/mvmc to avoid two copies.
MVMC_BUILD="${HERE}/../apbc_complexuhf/build/mvmc/build"
VMCDRY="${MVMC_BUILD}/src/mVMC/vmcdry.out"
VMC="${MVMC_BUILD}/src/mVMC/vmc.out"
[[ -x "${VMCDRY}" ]] || { echo "vmcdry.out missing; run ../apbc_complexuhf/build_complexuhf.sh first" >&2; exit 1; }
[[ -x "${VMC}" ]] || { echo "vmc.out missing; rerun build_complexuhf.sh after pulling the updated script" >&2; exit 1; }

WORK="${CASE_DIR}/work"
rm -rf "${WORK}"
mkdir -p "${WORK}"

# ---- step 1: H-wave UHFk SCF ----
HWAVE_WORK="${WORK}/hwave"
mkdir -p "${HWAVE_WORK}"
for f in CoulombIntra.dat Geometry.dat OneBodyG.dat Transfer.dat geometry_uhf.dat input.toml; do
  cp "${CASE_DIR}/${f}" "${HWAVE_WORK}/${f}"
done

(
  cd "${HWAVE_WORK}"
  # Shim the deprecated np.float alias before importing hwave (matches
  # apbc_complexuhf/run.sh). Use PYTHONPATH so this works whether the
  # package is editable-installed or not.
  PYTHONPATH="${ROOT}/src" python3 -c "
import numpy as np
if not hasattr(np, 'float'):
    np.float = float
import sys
sys.argv = ['hwave', 'input.toml']
from hwave.qlms import main
main()
"
)

[[ -f "${HWAVE_WORK}/output/eigen.npz" ]]
[[ -f "${HWAVE_WORK}/output/occupation.npz" ]]
[[ -f "${HWAVE_WORK}/output/greenone.dat" ]]
echo "  H-wave UHFk SCF: ${HWAVE_WORK}/output/*"

# ---- step 1.5: harness-gate — assert target occupation per case ----
if grep -qE "case_(pbc_sz2|zeeman_sz_free|soc_rashba_2d_nosub|soc_rashba_2d_nosub_apbc|soc_rashba_2d_sub|soc_rashba_2d_sub_apbc|soc_rashba_3d_sub_apbc_(xy|xz|yz|xyz))" <<< "${CASE}"; then
  python3 "${HERE}/scripts/assert_occupation.py" "${CASE_DIR}" "${HWAVE_WORK}" || {
    echo "harness-gate: occupation assertion failed for ${CASE}" >&2
    exit 1
  }
fi

# ---- step 1.6: SOC canonical-mix gate for case_soc_rashba_2d_sub ----
# Spec §4.2.3 criterion 5: SubShape [2,2,1] on CellShape [6,4,1] must
# yield a folded BZ with at least one non-self canonical block. If a
# future refactor of find_partner_rows / compute_canonical_reps or of
# the fixture's SubShape choice silently degenerates every folded k to a
# self-pair, the block-mix contract is broken and this fixture no longer
# exercises the non-self canonical path — fail loud.
if [[ "${CASE}" == "case_soc_rashba_2d_sub" ]]; then
  HWAVE_WORK="${HWAVE_WORK}" PYTHONPATH="${ROOT}/src:${ROOT}" python3 - <<'PYEOF' || exit 1
import os, sys
import numpy as np
from tools._uhfk_to_mvmc.general_fij_builder import compute_canonical_reps
from tools._uhfk_to_mvmc.partner_index import find_partner_rows

hwave_work = os.environ["HWAVE_WORK"]
eig = np.load(os.path.join(hwave_work, "output", "eigen.npz"), allow_pickle=False)
# Fixture pins BoundaryCondition = periodic in all directions -> theta = 0.
theta = np.zeros(3)
# CellShape=[6,4,1] / SubShape=[2,2,1] -> folded BZ [3,2,1].
L_folded = np.array([3, 2, 1], dtype=np.int64)
partner_rows, _ = find_partner_rows(eig["wavevector_index"], theta, L_folded)
canonical, self_pairs = compute_canonical_reps(partner_rows, eig["wavevector_index"])
n_non_self = len(canonical) - len(self_pairs)
if n_non_self < 1:
    print(
        f"harness-gate: case_soc_rashba_2d_sub reported {n_non_self} non-self "
        f"canonical blocks; SubShape choice degenerated",
        file=sys.stderr,
    )
    sys.exit(1)
print(f"  canonical mix: {n_non_self} non-self, {len(self_pairs)} self")
PYEOF
fi

# ---- step 2: vmcdry to produce mVMC .def files ----
MVMC_WORK="${WORK}/mvmc"
mkdir -p "${MVMC_WORK}"
cp "${CASE_DIR}/stan.in" "${MVMC_WORK}/stan.in"
cd "${MVMC_WORK}"
run_native "${VMCDRY}" stan.in </dev/null > vmcdry.log 2>&1
[[ -f "namelist.def" ]]
# StdFace's FermionHubbard generator writes ``ComplexType 0`` (real
# orbitalidx.def). The bridge produces complex F values (APBC carries a
# nontrivial phase, even under PBC the (k,-k) construction is complex
# off-diagonally), so flip the header to ``ComplexType 1`` here. The
# numerical content of the value rows remains real-zero in PBC and
# real+imag in APBC; mVMC accepts both interchangeably under
# ComplexType=1.
sed -i -E 's/^ComplexType +0/ComplexType 1/' orbitalidx.def
# v3 A/B cases: StdFace's FermionHubbard also emits orbitalidxgen.def
# (6-column General) when 2Sz != 0. Flip its ComplexType too.
if [[ -f orbitalidxgen.def ]]; then
  sed -i -E 's/^ComplexType +0/ComplexType 1/' orbitalidxgen.def
fi
echo "  vmcdry.out: ${MVMC_WORK}/{namelist,modpara,orbitalidx,trans,coulombintra,...}.def"

# ---- step 3: bridge ----
# Rank-lift noise amplitude. Default 1e-8 matches the bridge's own CLI
# default (stable plateau where mVMC <H> reaches UHF within VMC stderr
# at NVMCSample = 10000). Override per run via:
#     EPSILON_NOISE=1e-7 ./run.sh case_apbc
EPSILON_NOISE="${EPSILON_NOISE:-1.0e-8}"
RNG_SEED="${RNG_SEED:-7919}"

# v3 A/B routing: if StdFace produced orbitalidxgen.def (fires when
# 2Sz != 0 or lGC=1), prefer it and route the bridge through the
# General path. Otherwise stay on the v1/v2 AntiParallel-only path.
if [[ -f "${MVMC_WORK}/orbitalidxgen.def" ]]; then
  BRIDGE_ORBITALIDX="${MVMC_WORK}/orbitalidxgen.def"
else
  BRIDGE_ORBITALIDX="${MVMC_WORK}/orbitalidx.def"
fi

cd "${ROOT}"
# v3.1 spec §3.8: SOC path additionally emits mVMC trans.def because
# StdFace's FermionHubbardGC generator drops Rashba s != t transfer
# entries. Non-SOC fixtures leave vmcdry's trans.def untouched.
BRIDGE_SOC_ARGS=()
if [[ "${CASE}" == case_soc_* ]]; then
  BRIDGE_SOC_ARGS=(
    --transfer   "${HWAVE_WORK}/Transfer.dat"
    --emit-trans "${MVMC_WORK}/trans.def"
  )
fi
# v3.2 spec §1: SOC + SubShape > [1, 1, 1] needs a bridge-emitted
# orbitalidxgen.def with all-unique classes (StdFace's over-grouping
# under SOC + sublattice folding trips the class-consistency guard).
# The bridge writes to a sidecar path here; the E2E harness then moves
# it over StdFace's orbitalidxgen.def before mVMC runs so the mVMC
# InOrbitalGeneral consumer sees the same classes the bridge used.
BRIDGE_ORBITALIDX_ARGS=()
BRIDGE_EMIT_ORBITALIDX_PATH=""
if [[ "${CASE}" == "case_soc_rashba_2d_sub" || \
      "${CASE}" == "case_soc_rashba_2d_sub_apbc" || \
      "${CASE}" == case_soc_rashba_3d_sub_apbc_* ]]; then
  BRIDGE_EMIT_ORBITALIDX_PATH="${MVMC_WORK}/orbitalidxgen.def.bridge"
  BRIDGE_ORBITALIDX_ARGS=(
    --emit-orbitalidx "${BRIDGE_EMIT_ORBITALIDX_PATH}"
  )
fi
# Seven-gate cases inject --debug-writer so F_pre/F_post are dumped
# into ${MVMC_WORK} and later frozen into ${WORK}/bridge/ for G0/G2a.
BRIDGE_DEBUG_ARGS=()
if [[ "${CASE}" == "case_soc_rashba_2d_sub_apbc" || \
      "${CASE}" == case_soc_rashba_3d_sub_apbc_* ]]; then
  BRIDGE_DEBUG_ARGS=(--debug-writer)
fi
python3 tools/uhfk_to_mvmc.py \
    --input        "${HWAVE_WORK}/input.toml" \
    --eigen        "${HWAVE_WORK}/output/eigen.npz" \
    --occupation   "${HWAVE_WORK}/output/occupation.npz" \
    --geometry     "${HWAVE_WORK}/geometry_uhf.dat" \
    --orbitalidx   "${BRIDGE_ORBITALIDX}" \
    --output       "${MVMC_WORK}/zqp_orbital_uhfk.dat" \
    --check-density \
    --onebodyg-uhf "${HWAVE_WORK}/output/greenone.dat" \
    --epsilon-noise "${EPSILON_NOISE}" \
    --rng-seed "${RNG_SEED}" \
    "${BRIDGE_SOC_ARGS[@]}" \
    "${BRIDGE_ORBITALIDX_ARGS[@]}" \
    "${BRIDGE_DEBUG_ARGS[@]}"
echo "  bridge wrote: ${MVMC_WORK}/zqp_orbital_uhfk.dat (density check OK, epsilon=${EPSILON_NOISE})"
if [[ "${CASE}" == case_soc_* ]]; then
  echo "  bridge wrote: ${MVMC_WORK}/trans.def (SOC: Rashba entries preserved)"
fi
# SOC + SubShape > 1 override: bridge's all-unique-classes
# orbitalidxgen.def replaces StdFace's over-grouped version so mVMC
# consumes the same class layout the bridge used to write
# zqp_orbital_uhfk.dat.
if [[ -n "${BRIDGE_EMIT_ORBITALIDX_PATH}" ]]; then
  mv "${BRIDGE_EMIT_ORBITALIDX_PATH}" "${MVMC_WORK}/orbitalidxgen.def"
  echo "  bridge wrote: ${MVMC_WORK}/orbitalidxgen.def (SOC+SubShape: all-unique classes)"
fi

# ---- step 4: register InOrbital(General) in namelist.def ----
NAMELIST="${MVMC_WORK}/namelist.def"
# Remove any prior entry for InOrbital* lines, then append the bridge file
# under the appropriate keyword.
sed -i '/^InOrbital/d' "${NAMELIST}"
sed -i '/^InOrbitalAntiParallel/d' "${NAMELIST}"
sed -i '/^InOrbitalGeneral/d' "${NAMELIST}"
if [[ "${BRIDGE_ORBITALIDX}" == *orbitalidxgen.def ]]; then
  # v3 General path: tell mVMC to consume the 6-column class table and
  # bridge params via InOrbitalGeneral. Uncomment the commented
  # "OrbitalGeneral" line StdFace emitted and comment the AntiParallel
  # "Orbital" line so mVMC does NOT try to parse orbitalidx.def as the
  # main pair table.
  sed -i -E 's|^# OrbitalGeneral|  OrbitalGeneral|' "${NAMELIST}"
  sed -i -E 's|^         Orbital  |#        Orbital  |' "${NAMELIST}"
  # StdFace also emits OrbitalParallel orbitalidxpara.def for 2Sz != 0
  # (grand-canonical-style triplet pair). mVMC rejects when more than one
  # OrbitalX* keyword is active alongside OrbitalGeneral, so comment it.
  sed -i -E 's|^ OrbitalParallel|# OrbitalParallel|' "${NAMELIST}"
  echo "InOrbitalGeneral ${MVMC_WORK}/zqp_orbital_uhfk.dat" >> "${NAMELIST}"
else
  echo "InOrbital ${MVMC_WORK}/zqp_orbital_uhfk.dat" >> "${NAMELIST}"
fi
# Disable Gutzwiller / Jastrow / GeneralRBM projection lines: the bridge
# produces a pure Slater initial WF, and StdFace defaults to nonzero
# projection parameters that would multiply the Slater by a non-trivial
# correlator. Without InGutzwiller / InJastrow files mVMC starts those
# parameters from random values, giving an energy that does NOT match
# the UHF Slater. Commenting these lines makes mVMC ignore the
# correlators entirely.
sed -i 's/^Gutzwiller/#Gutzwiller/' "${NAMELIST}"
sed -i 's/^Jastrow/#Jastrow/' "${NAMELIST}"
echo "  namelist.def updated to read InOrbital ${MVMC_WORK}/zqp_orbital_uhfk.dat"
echo "  Gutzwiller / Jastrow projection disabled (pure Slater test)"

# ---- step 5: vmc.out (NVMCCalMode=1, NSROptItrStep=0 -> measurement only) ----
cd "${MVMC_WORK}"
run_native "${VMC}" namelist.def > vmc.log 2>&1 \
    || { echo "vmc.out failed; see ${MVMC_WORK}/vmc.log" >&2; tail -50 vmc.log >&2; exit 1; }

# zvo_out_xxx.dat holds bin-by-bin <H>, <H^2>, ... lines.
OUT_FILE="$(ls output/zvo_out_*.dat 2>/dev/null | head -n1)"
[[ -n "${OUT_FILE}" ]] || { echo "mVMC produced no zvo_out output" >&2; exit 1; }
echo "  vmc.out -> ${MVMC_WORK}/${OUT_FILE}"

# ---- v3.6/v3.7 SOC + SubShape APBC: seven-gate branch ----
# Other cases fall through to the legacy 3-arg energy compare at the bottom.
if [[ "${CASE}" == "case_soc_rashba_2d_sub_apbc" ]] || \
   [[ "${CASE}" == case_soc_rashba_3d_sub_apbc_* ]]; then

  # Static/live gates apply to the 3D APBC cases.
  PHASE5_EXPECTED_MASK=""
  case "${CASE}" in
    case_soc_rashba_3d_sub_apbc_xy)
      PHASE5_EXPECTED_MASK="1, 1, 0"
      ;;
    case_soc_rashba_3d_sub_apbc_xz)
      PHASE5_EXPECTED_MASK="1, 0, 1"
      ;;
    case_soc_rashba_3d_sub_apbc_yz)
      PHASE5_EXPECTED_MASK="0, 1, 1"
      ;;
    case_soc_rashba_3d_sub_apbc_xyz)
      PHASE5_EXPECTED_MASK="1, 1, 1"
      ;;
    case_soc_rashba_3d_sub_apbc_*)
      echo "ERROR: unknown case name '${CASE}'; define PHASE5_EXPECTED_MASK" >&2
      echo "for this case before enabling it." >&2
      exit 1
      ;;
  esac

  # normalize_mvmc_output per spec §5.6: pick the unique zvo_out_001.dat
  # under the two v3.6-supported layouts (${WORK_DIR}/mvmc/output/,
  # ${WORK_DIR}/mvmc/) and copy it to ${WORK_DIR}/mvmc/zvo_out_selected.dat.
  # Bash tests in tests/test_normalize_mvmc_output_bash_v36.py pin the
  # 0/multiple/nested cases.
  normalize_mvmc_output() {
    local -a candidate_dirs=(
      "${WORK}/mvmc/output"
      "${WORK}/mvmc"
    )
    local selector="${ZVO_OUT_SELECTOR:-zvo_out_001.dat}"
    local -a found=()
    for dir in "${candidate_dirs[@]}"; do
      if [[ -d "$dir" ]]; then
        while IFS= read -r -d '' f; do
          found+=("$f")
        done < <(find "$dir" -maxdepth 1 -type f -name "$selector" -print0)
      fi
    done
    if (( ${#found[@]} == 0 )); then
      echo "G3-normalize: no ${selector} in ${candidate_dirs[*]}" >&2
      return 1
    fi
    if (( ${#found[@]} > 1 )); then
      echo "G3-normalize: multiple ${selector}: ${found[*]}" >&2
      return 1
    fi
    cp "${found[0]}" "${WORK}/mvmc/zvo_out_selected.dat"
  }

  # ---- step 6.1: bridge zero-noise producer for G0-writer-check ----
  # Re-runs the bridge with --epsilon-noise 0 --debug-writer into
  # ${WORK}/bridge_zeronoise so parse_emitted_F(bridge_zeronoise) has
  # rank exactly 2*N_pairs (no rank-lift noise applied). G0-writer-check
  # asserts that this zero-noise F equals the F reconstructed from the
  # aggregated (mapping, params) at 1e-10 (see spec §1.4).
  BRIDGE_ZERO="${WORK}/bridge_zeronoise"
  mkdir -p "${BRIDGE_ZERO}"
  cd "${ROOT}"
  BRIDGE_ZN_ORBIDX_ARGS=()
  if [[ -n "${BRIDGE_EMIT_ORBITALIDX_PATH}" ]]; then
    BRIDGE_ZN_ORBIDX_ARGS=(
      --emit-orbitalidx "${BRIDGE_ZERO}/orbitalidxgen.def.bridge"
    )
  fi
  python3 tools/uhfk_to_mvmc.py \
      --input        "${HWAVE_WORK}/input.toml" \
      --eigen        "${HWAVE_WORK}/output/eigen.npz" \
      --occupation   "${HWAVE_WORK}/output/occupation.npz" \
      --geometry     "${HWAVE_WORK}/geometry_uhf.dat" \
      --orbitalidx   "${BRIDGE_ORBITALIDX}" \
      --output       "${BRIDGE_ZERO}/zqp_orbital_uhfk.dat" \
      --no-check-density \
      --epsilon-noise 0 \
      --rng-seed "${RNG_SEED}" \
      --debug-writer \
      "${BRIDGE_SOC_ARGS[@]}" \
      "${BRIDGE_ZN_ORBIDX_ARGS[@]}"
  # namelist.def / modpara.def / orbitalidx_general.def all needed by
  # parse_emitted_F. Copy from the shipping bridge dir (mVMC) into the
  # zero-noise dir so both workspaces are self-contained.
  cp "${MVMC_WORK}/namelist.def" "${BRIDGE_ZERO}/namelist.def"
  cp "${MVMC_WORK}/modpara.def"  "${BRIDGE_ZERO}/modpara.def"
  cp "${MVMC_WORK}/orbitalidxgen.def" "${BRIDGE_ZERO}/orbitalidx_general.def"
  echo "  bridge zero-noise -> ${BRIDGE_ZERO}/"

  # ---- step 6.2: freeze shipping bridge outputs into ${WORK}/bridge ----
  # The compare.py §5.3b dispatcher expects gate inputs under
  # ${WORK}/bridge/{...}; the shipping bridge wrote to ${MVMC_WORK}/ so
  # we materialize the required subset into ${WORK}/bridge/.
  BRIDGE_SHIP="${WORK}/bridge"
  mkdir -p "${BRIDGE_SHIP}"
  # F_pre_noise / F_post_aggregate: v3.6 --debug-writer must have been
  # set on step 3's bridge invocation. Re-issue if missing.
  if [[ ! -f "${MVMC_WORK}/F_pre_noise.npz" ]]; then
    echo "  step 6.2: re-running shipping bridge with --debug-writer to dump F snapshots"
    python3 tools/uhfk_to_mvmc.py \
        --input        "${HWAVE_WORK}/input.toml" \
        --eigen        "${HWAVE_WORK}/output/eigen.npz" \
        --occupation   "${HWAVE_WORK}/output/occupation.npz" \
        --geometry     "${HWAVE_WORK}/geometry_uhf.dat" \
        --orbitalidx   "${BRIDGE_ORBITALIDX}" \
        --output       "${BRIDGE_SHIP}/zqp_orbital_uhfk.dat" \
        --no-check-density \
        --epsilon-noise "${EPSILON_NOISE}" \
        --rng-seed "${RNG_SEED}" \
        --debug-writer \
        "${BRIDGE_SOC_ARGS[@]}"
  else
    cp "${MVMC_WORK}/F_pre_noise.npz"       "${BRIDGE_SHIP}/F_pre_noise.npz"
    cp "${MVMC_WORK}/F_post_aggregate.npz"  "${BRIDGE_SHIP}/F_post_aggregate.npz"
    cp "${MVMC_WORK}/zqp_orbital_uhfk.dat"  "${BRIDGE_SHIP}/zqp_orbital_uhfk.dat"
  fi
  cp "${MVMC_WORK}/namelist.def"       "${BRIDGE_SHIP}/namelist.def"
  cp "${MVMC_WORK}/modpara.def"        "${BRIDGE_SHIP}/modpara.def"
  cp "${MVMC_WORK}/orbitalidxgen.def"  "${BRIDGE_SHIP}/orbitalidx_general.def"
  echo "  bridge shipping snapshot -> ${BRIDGE_SHIP}/"

  # ---- step 6.3: ComplexUHF SCF for G2a/G2b ----
  if [[ "${CASE}" == "case_soc_rashba_2d_sub_apbc" ]]; then
    COMPLEXUHF_CASE="${HERE}/../apbc_complexuhf/case_soc_rashba_2d_sub_apbc"
  else
    COMPLEXUHF_CASE="${HERE}/../apbc_complexuhf/${CASE}"
  fi
  COMPLEXUHF_WORK="${WORK}/complexuhf"
  UHF="${HERE}/../apbc_complexuhf/build/mvmc/build/src/ComplexUHF/UHF"
  if [[ ! -x "${UHF}" ]]; then
    echo "ERROR: ComplexUHF UHF binary missing at ${UHF}; run apbc_complexuhf/build_complexuhf.sh" >&2
    exit 1
  fi
  mkdir -p "${COMPLEXUHF_WORK}"
  cp "${COMPLEXUHF_CASE}/stan.in" "${COMPLEXUHF_WORK}/stan.in"
  (
    cd "${COMPLEXUHF_WORK}"
    complexuhf_files_before_vmcdry=()
    case "${CASE}" in
      case_soc_rashba_3d_sub_apbc_*)
        mapfile -d '' -t complexuhf_files_before_vmcdry < <(
          find . -mindepth 1 -maxdepth 1 -type f -printf '%f\0' | sort -z
        )
        ;;
      case_soc_rashba_2d_sub_apbc)
        ;;
      *)
        echo "ERROR: unsupported ComplexUHF snapshot case: ${CASE}" >&2
        exit 1
        ;;
    esac
    run_native "${VMCDRY}" stan.in </dev/null > vmcdry.log 2>&1
    complexuhf_vmcdry_files=()
    case "${CASE}" in
      case_soc_rashba_3d_sub_apbc_*)
        mapfile -d '' -t complexuhf_files_after_vmcdry < <(
          find . -mindepth 1 -maxdepth 1 -type f -printf '%f\0' | sort -z
        )
        mapfile -d '' -t complexuhf_vmcdry_files < <(
          comm -z -13 \
            <(printf '%s\0' "${complexuhf_files_before_vmcdry[@]}") \
            <(printf '%s\0' "${complexuhf_files_after_vmcdry[@]}")
        )
        ;;
      case_soc_rashba_2d_sub_apbc)
        ;;
      *)
        echo "ERROR: unsupported ComplexUHF snapshot case: ${CASE}" >&2
        exit 1
        ;;
    esac
    # BEGIN task6a2-authoritative-bundle-staging
    case "${CASE}" in
      case_soc_rashba_3d_sub_apbc_*)
        complexuhf_bundle_files=(namelist.def modpara.def locspn.def geometry.dat coulombintra.def orbitalidx.def orbitalidxpara.def)
        for bundle_file in "${complexuhf_bundle_files[@]}"; do
          cp "${COMPLEXUHF_CASE}/${bundle_file}" "./${bundle_file}" || {
            echo "ERROR: ${COMPLEXUHF_CASE}: failed to stage ${bundle_file}" >&2
            exit 1
          }
          if ! cmp -s "${COMPLEXUHF_CASE}/${bundle_file}" "./${bundle_file}"; then
            echo "ERROR: ${COMPLEXUHF_CASE}: staged bundle mismatch: ${bundle_file}" >&2
            exit 1
          fi
        done

        # StdFace output is only a discard-and-verify baseline for v3.7.
        # Remove every newly-created vmcdry file outside the committed bundle
        # and the files supplied or retained by adjacent harness steps.
        for workspace_file in "${complexuhf_vmcdry_files[@]}"; do
          case "${workspace_file}" in
            namelist.def|modpara.def|locspn.def|geometry.dat|coulombintra.def|orbitalidx.def|orbitalidxpara.def|trans.def|initial.def|vmcdry.log)
              ;;
            *)
              rm -f -- "./${workspace_file}" || {
                echo "ERROR: ${COMPLEXUHF_CASE}: failed to delete unexpected vmcdry file: ${workspace_file}" >&2
                exit 1
              }
              ;;
          esac
        done
        for workspace_file in "${complexuhf_vmcdry_files[@]}"; do
          case "${workspace_file}" in
            namelist.def|modpara.def|locspn.def|geometry.dat|coulombintra.def|orbitalidx.def|orbitalidxpara.def|trans.def|initial.def|vmcdry.log)
              ;;
            *)
              if [[ -e "./${workspace_file}" ]]; then
                echo "ERROR: ${COMPLEXUHF_CASE}: unexpected vmcdry file survived staging: ${workspace_file}" >&2
                exit 1
              fi
              ;;
          esac
        done
        ;;
      case_soc_rashba_2d_sub_apbc)
        ;;
      *)
        echo "ERROR: unsupported ComplexUHF staging case: ${CASE}" >&2
        exit 1
        ;;
    esac
    # END task6a2-authoritative-bundle-staging
    # StdFace's FermionHubbardGC drops Rashba s != t entries. Replace
    # StdFace's trans.def with the H-wave bridge's SOC-preserving
    # trans.def (built from Transfer.dat) so ComplexUHF's Hamiltonian
    # matches H-wave's. Same substitution as v3.1 SOC (§3.8).
    cp "${MVMC_WORK}/trans.def" "./trans.def"
    # Apply the case-specific ComplexUHF modpara override (NMPTrans etc.)
    case "${CASE}" in
      case_soc_rashba_2d_sub_apbc)
        bash "${COMPLEXUHF_CASE}/complexuhf_modpara_override.txt"
        ;;
      case_soc_rashba_3d_sub_apbc_*)
        # The committed modpara.def is authoritative.
        ;;
      *)
        echo "ERROR: unsupported ComplexUHF override case: ${CASE}" >&2
        exit 1
        ;;
    esac
  )
  # Codex adversarial-review 2026-07-12 hardening: seed ComplexUHF with
  # H-wave's converged density in the physical basis. Without seeding,
  # ComplexUHF's random init lands in a different broken-symmetry
  # minimum (same energy at 0.07%, different density at 4.76e-2)
  # rendering G2a/G2b element-level comparisons meaningless. With
  # seeding, both solvers converge to the SAME minimum and G2a/G2b at
  # 1e-6 becomes achievable.
  # The seed perturbation must be large enough that G2 starts OUTSIDE
  # its own tolerance -- otherwise a no-op solver passes by handing the
  # seed straight back -- and small enough to stay in H-wave's basin.
  # All fixtures in this seven-gate branch now use flag_fock=true and
  # share ComplexUHF's functional. At 1e-3 the measured initial density
  # delta is hundreds of times the 1e-6 tolerance, while remaining well
  # inside the common attraction basin and converging to O(1e-8).
  PERTURB_SCALE="1e-3"
  echo "  seed perturb-scale = ${PERTURB_SCALE} (shared Fock-consistent value)"
  python3 "${HERE}/scripts/seed_complexuhf_from_hwave.py" \
      --hwave-workspace "${HWAVE_WORK}" \
      --complexuhf-workspace "${COMPLEXUHF_WORK}" \
      --initial-def initial.def \
      --perturb-scale "${PERTURB_SCALE}"

  (
    cd "${COMPLEXUHF_WORK}"
    if [[ -n "${PHASE5_EXPECTED_MASK}" ]]; then
      if ! python3 - "${ROOT}" "${COMPLEXUHF_CASE}" "${WORK}" "${UHF}" \
          "${PHASE5_EXPECTED_MASK}" > uhf.log 2>&1 <<'PY'
import sys
from pathlib import Path

sys.path.insert(0, sys.argv[1])
from tools._uhfk_to_mvmc.phase5_gate import run_phase5_gate

fixture_root = Path(sys.argv[2])
try:
    result = run_phase5_gate(
        fixture_root=fixture_root,
        workspace=Path(sys.argv[3]),
        expected_phase_mask=tuple(
            int(value.strip()) for value in sys.argv[5].split(",")
        ),
        uhf_binary=sys.argv[4],
    )
except Exception as exc:
    raise SystemExit(
        f"Phase 5 gate failed for {fixture_root.name}: {exc}"
    ) from exc
sys.stdout.write(result.stdout)
sys.stderr.write(result.stderr)
PY
      then
        echo "Gated ComplexUHF run failed (gate check or solver); tail uhf.log:" >&2
        tail -50 uhf.log >&2
        exit 1
      fi
    else
      # Preserve the established launch path byte-for-byte.
      run_native "${UHF}" namelist.def > uhf.log 2>&1 || {
        echo "ComplexUHF UHF failed; tail uhf.log:" >&2
        tail -50 uhf.log >&2
        exit 1
      }
    fi
    # Residual SCF-log sanity check. The primary non-vacuity invariant
    # lives in compare.py: every G2 gate requires the seeded density to
    # start at least 10*tol outside H-wave's reference and the converged
    # density to finish inside tol. With the restored shipped
    # perturb_scale=1e-3, ComplexUHF should also report non-zero steps.
    scf_steps=$(awk '/finished at/ {print $(NF-1); exit}' uhf.log 2>/dev/null)
    if [[ -z "${scf_steps}" ]]; then
      echo "ComplexUHF: could not parse SCF step count from uhf.log" >&2
      tail -20 uhf.log >&2
      exit 1
    fi
    if (( scf_steps < 1 )); then
      echo "ComplexUHF terminated at ${scf_steps} SCF steps; the seed" >&2
      echo "was already at the SCF fixed point (circular G2a/G2b)." >&2
      echo "seed_complexuhf_from_hwave must apply a non-zero" >&2
      echo "--perturb-scale ${PERTURB_SCALE} so the solver actually iterates." >&2
      exit 1
    fi
  )
  echo "  ComplexUHF SCF -> ${COMPLEXUHF_WORK}/zvo_UHF_cisajs.dat"

  # ---- step 6.4: normalize mVMC output + surface hwave/energy.dat ----
  normalize_mvmc_output || { echo "G3-normalize failed" >&2; exit 1; }
  echo "  normalize_mvmc_output -> ${WORK}/mvmc/zvo_out_selected.dat"
  # Spec §5.6 canonical path for G3 is ${WORK_DIR}/hwave/energy.dat;
  # H-wave writes to hwave/output/energy.dat, so surface a copy at the
  # canonical location.
  cp "${HWAVE_WORK}/output/energy.dat" "${HWAVE_WORK}/energy.dat"
  # G4 guard reads ${WORK}/hwave/green.npz + eigen.npz + occupation.npz;
  # H-wave writes them under hwave/output/. Surface canonical copies.
  cp "${HWAVE_WORK}/output/green.npz"      "${HWAVE_WORK}/green.npz"
  cp "${HWAVE_WORK}/output/eigen.npz"      "${HWAVE_WORK}/eigen.npz"
  cp "${HWAVE_WORK}/output/occupation.npz" "${HWAVE_WORK}/occupation.npz"

  # ---- step 7: seven-gate check ----
  # Each gate command emits ONE anchored PASS line per §5.3b metadata
  # table; awk index() verifies exactly one line starts at column 1 with
  # the expected `<GATE> PASS mode=X artifact_source=Y helper=Z ` prefix.
  GATE_LOGS="${WORK}/gate_logs"
  mkdir -p "${GATE_LOGS}"
  # Manifest for G4.
  COMPOSITE_MANIFEST="${CASE_DIR}/composite_element.json"

  run_gate() {
    # Usage: run_gate GATE_NAME EXPECTED_PREFIX CMD...
    local gate="$1" prefix="$2" ; shift 2
    local log="${GATE_LOGS}/${gate}.log"
    "$@" > "${log}" 2>&1 || {
      echo "  ${gate} FAILED (non-zero exit): ${log}" >&2
      tail -20 "${log}" >&2 || true
      return 1
    }
    awk -v p="${prefix}" 'index($0, p) == 1 { n++ } END { exit (n == 1 ? 0 : 1) }' "${log}" || {
      echo "  ${gate} FAILED: no anchored PASS matching '${prefix}' in ${log}" >&2
      tail -20 "${log}" >&2 || true
      return 1
    }
    echo "  ${gate} PASS"
  }

  # Absolute path to compare.py (workspace-mode dispatcher) and G4 guard.
  COMPARE_PY="${HERE}/compare.py"
  G4_GUARD="${HERE}/scripts/soc_apbc_topology_guard.py"

  cd "${ROOT}"
  # compare.py + G4 guard import `tools._uhfk_to_mvmc.*` which sits under
  # ${ROOT}; PYTHONPATH must include ${ROOT} for both to resolve. The v3.6
  # dispatch tests already pin this same import path.
  export PYTHONPATH="${ROOT}/src:${ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
  run_gate G0-writer-check \
    "G0-writer-check PASS mode=g0-writer-check artifact_source=bridge_zeronoise helper=parse_emitted_F " \
    python3 "${COMPARE_PY}" --workspace "${WORK}" --mode g0-writer-check --gtol_writer 1e-10
  run_gate G1 \
    "G1 PASS mode=g1 artifact_source=bridge+hwave helper=build_slater_orbitals+gauge_lift " \
    python3 "${COMPARE_PY}" --workspace "${WORK}" --mode g1 --gtol_g1 1e-10
  run_gate G2a-emitted-F \
    "G2a-emitted-F PASS mode=g2a-emitted-F artifact_source=bridge+complexuhf helper=parse_emitted_F+pair_product_density_from_F " \
    python3 "${COMPARE_PY}" --workspace "${WORK}" --mode g2a-emitted-F --gtol_ship 1e-6
  run_gate G2a-in-memory-A \
    "G2a-in-memory-A PASS mode=g2a-in-memory artifact_source=hwave+complexuhf helper=build_slater_orbitals " \
    python3 "${COMPARE_PY}" --workspace "${WORK}" --mode g2a-in-memory --gtol_ship 1e-6
  run_gate G2b \
    "G2b PASS mode=g2b artifact_source=hwave+complexuhf helper=gauge_lift " \
    python3 "${COMPARE_PY}" --workspace "${WORK}" --mode g2b --gtol_gauge 1e-6
  run_gate G3 \
    "G3 PASS mode=g3 artifact_source=hwave+mvmc helper=energy_relative_delta " \
    python3 "${COMPARE_PY}" --workspace "${WORK}" --mode g3 --etol 0.01
  run_gate G4 \
    "G4 PASS mode=g4 artifact_source=hwave+bridge+composite-manifest helper=soc_apbc_topology_guard " \
    python3 "${G4_GUARD}" --workspace "${WORK}" --mode g4 --composite-manifest "${COMPOSITE_MANIFEST}"

  echo "SEVEN-GATE PASS on ${CASE}"
  exit 0
fi

# ---- step 6: compare energies (legacy 3-arg, non-v3.6 cases) ----
cd "${HERE}"
python3 compare.py "${HWAVE_WORK}/output/energy.dat" "${MVMC_WORK}/${OUT_FILE}" "${CASE}"
