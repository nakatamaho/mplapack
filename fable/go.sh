#!/usr/bin/env bash
# go.sh: multi-pass Fortran->C++ conversion with signature convergence detection
set -euo pipefail
shopt -s nullglob

PASSES="${1:-2}"

ROOT="/home/docker/mplapack"
# Ensure Python can import the in-tree 'fable' package.
export PYTHONPATH="${ROOT}:${PYTHONPATH:-}"

LAPACK_ROOT="${ROOT}/external/lapack/work/internal/lapack-3.9.1"
BLAS_SRC="${LAPACK_ROOT}/BLAS/SRC"
LAPACK_SRC="${LAPACK_ROOT}/SRC"

FABLE="${ROOT}/fable"
MPBLAS_REF="${ROOT}/mpblas/reference"
MPLAPACK_REF="${ROOT}/mplapack/reference"

MPBLAS_HDR="${MPBLAS_REF}/mpblas_generic.h"
MPLAPACK_HDR="${MPLAPACK_REF}/mplapack_generic.h"
SIG_PY="${FABLE}/mplapack_signatures.py"

move_cpp_or_die() {
  local src_dir="$1"
  local dst_dir="$2"
  local files=("${src_dir}"/*.cpp)

  # Fail fast if conversion produced no .cpp outputs.
  if (( ${#files[@]} == 0 )); then
    echo "ERROR: No .cpp files produced in: ${src_dir}" >&2
    exit 1
  fi

  /bin/mv -f "${files[@]}" "${dst_dir}/"
}

# Return the path to the per-pass signature snapshot file.
sig_snapshot_path() {
  local pass="$1"
  echo "${SIG_PY%.py}.pass${pass}.py"
}

# Remove any existing signature snapshots from older runs.
cleanup_signature_snapshots() {
  local snaps=( "${SIG_PY%.py}.pass"*.py )
  if (( ${#snaps[@]} )); then
    /bin/rm -f "${snaps[@]}"
  fi
}

# Save a per-pass snapshot of mplapack_signatures.py for later diffing.
snapshot_signatures() {
  local pass="$1"
  local dst
  dst="$(sig_snapshot_path "${pass}")"
  /bin/cp -f "${SIG_PY}" "${dst}"
  echo "[SIG] snapshot: ${dst}"
}

# Compare two signature snapshots and print a human-readable summary.
# Exit code: 0 if converged (no semantic changes), 1 if changed, 2 on error.
report_signature_diff_and_check_converged() {
  local prev_py="$1"
  local cur_py="$2"
  python3 - "${prev_py}" "${cur_py}" <<'PY'
import sys
import importlib.util
import uuid

prev_path, cur_path = sys.argv[1], sys.argv[2]

def load(path: str):
    name = f"sig_{uuid.uuid4().hex}"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load signatures from: {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    fs = getattr(mod, "FUNCTION_SIGNATURES", {}) or {}
    fr = getattr(mod, "FUNCTION_RETURNS", {}) or {}
    fs = {str(k).lower(): v for (k, v) in dict(fs).items()}
    fr = {str(k).lower(): v for (k, v) in dict(fr).items()}
    return fs, fr

def delta(d0, d1):
    k0, k1 = set(d0), set(d1)
    added = sorted(k1 - k0)
    removed = sorted(k0 - k1)
    common = k0 & k1
    changed = sorted(k for k in common if d0[k] != d1[k])
    return added, removed, changed

try:
    fs0, fr0 = load(prev_path)
    fs1, fr1 = load(cur_path)
except Exception as e:
    print(f"[SIGDIFF] ERROR: {e}", file=sys.stderr)
    sys.exit(2)

added, removed, changed = delta(fs0, fs1)
r_added, r_removed, r_changed = delta(fr0, fr1)

print(f"[SIGDIFF] {prev_path} -> {cur_path}")
print(f"[SIGDIFF] FUNCTION_SIGNATURES: prev={len(fs0)} cur={len(fs1)} added={len(added)} removed={len(removed)} changed={len(changed)}")
if added:
    print("[SIGDIFF] added:")
    for k in added:
        print(f"  + {k}")
if removed:
    print("[SIGDIFF] removed:")
    for k in removed:
        print(f"  - {k}")
if changed:
    print("[SIGDIFF] changed:")
    for k in changed:
        print(f"  * {k}: {fs0[k]} -> {fs1[k]}")

print(f"[SIGDIFF] FUNCTION_RETURNS: prev={len(fr0)} cur={len(fr1)} added={len(r_added)} removed={len(r_removed)} changed={len(r_changed)}")
if r_added:
    print("[SIGDIFF] returns added:")
    for k in r_added:
        print(f"  + {k}")
if r_removed:
    print("[SIGDIFF] returns removed:")
    for k in r_removed:
        print(f"  - {k}")
if r_changed:
    print("[SIGDIFF] returns changed:")
    for k in r_changed:
        print(f"  * {k}: {fr0[k]} -> {fr1[k]}")

converged = (not added and not removed and not changed and not r_added and not r_removed and not r_changed)
print(f"[SIGDIFF] converged={converged}")
sys.exit(0 if converged else 1)
PY
}

# ------------------------------------------------------------
# Internal seed pass to bootstrap signatures (reduces outer passes)
# ------------------------------------------------------------
LAPACK_SEED_SRCS1=(
  dgeesx.f
  dgeevx.f
  dggevx.f
  dggsvd3.f
  dlaed1.f
  dlaed7.f
  dlags2.f
  dlarre.f
  dlartgs.f
  dlasd6.f
  dlasq3.f
  dlatdf.f
  dstegr.f
  dtgsy2.f
  zgeesx.f
  zgeevx.f
  zggevx.f
  zggsvd3.f
  zlaed7.f
  zlaev2.f
  zlags2.f
  zlatdf.f
  zstegr.f
  ztgsy2.f
)

LAPACK_SEED_SRCS2=(
  dtgsy2.f
  ztgsy2.f
)

run_one_pass() {
  local pass="$1"
  echo "=== PASS ${pass}/${PASSES} ==="

  # Remove leftovers from a previous aborted run to avoid mixing outputs.
  /bin/rm -f "${BLAS_SRC}"/*.cpp "${LAPACK_SRC}"/*.cpp

  ( cd "${BLAS_SRC}"  && bash "${FABLE}/convert_blas_all.sh" )
  move_cpp_or_die "${BLAS_SRC}" "${MPBLAS_REF}"

  ( cd "${LAPACK_SRC}" && bash "${FABLE}/convert_lapack_all.sh" )
  move_cpp_or_die "${LAPACK_SRC}" "${MPLAPACK_REF}"

  # Generate prototype headers.
  bash "${FABLE}/gen_include_mpblas.sh"
  bash "${FABLE}/gen_include_mplapack.sh"

  # Generate signatures so the next pass can use them.
  python3 "${FABLE}/gen_mplapack_signatures.py" "${MPBLAS_HDR}" "${MPLAPACK_HDR}" > "${SIG_PY}"

  ### internal pass to reduce iterations ###
  ( cd "${LAPACK_SRC}" && parallel -j "${JOBS:-$(nproc)}" 'echo "Converting {}"; bash "'"${FABLE}/convert_lapack.sh"'" "{}"' ::: "${LAPACK_SEED_SRCS1[@]}" )
  move_cpp_or_die "${LAPACK_SRC}" "${MPLAPACK_REF}"

  # Generate prototype headers.
  bash "${FABLE}/gen_include_mpblas.sh"
  bash "${FABLE}/gen_include_mplapack.sh"

  # Generate signatures so the next pass can use them.
  python3 "${FABLE}/gen_mplapack_signatures.py" "${MPBLAS_HDR}" "${MPLAPACK_HDR}" > "${SIG_PY}"

  ### internal pass to reduce iterations ###
  ( cd "${LAPACK_SRC}" && parallel -j "${JOBS:-$(nproc)}" 'echo "Converting {}"; bash "'"${FABLE}/convert_lapack.sh"'" "{}"' ::: "${LAPACK_SEED_SRCS2[@]}" )
  move_cpp_or_die "${LAPACK_SRC}" "${MPLAPACK_REF}"

  # Generate prototype headers.
  bash "${FABLE}/gen_include_mpblas.sh"
  bash "${FABLE}/gen_include_mplapack.sh"

  # Generate signatures so the next pass can use them.
  python3 "${FABLE}/gen_mplapack_signatures.py" "${MPBLAS_HDR}" "${MPLAPACK_HDR}" > "${SIG_PY}"

  echo "[SEED] signatures updated: ${SIG_PY}"
  ### internal pass to reduce iterations ###

  echo "PASS ${pass} done. signatures: ${SIG_PY}"
}

rm -f "${FABLE}/mplapack_name_map.txt"
bash "${FABLE}/gen_mplapack_name_map.sh"

rm -f "${SIG_PY}"
cleanup_signature_snapshots
rm -f "${MPBLAS_HDR}"
rm -f "${MPLAPACK_HDR}"

for pass in $(seq 1 "${PASSES}"); do
  if [ "${pass}" -ge 2 ]; then
    export MPLAPACK_SIGNATURES_PY="${SIG_PY}"
    export FABLE_DEBUG_SIGNATURES=1
  fi

  run_one_pass "${pass}"
  snapshot_signatures "${pass}"

  if [ "${pass}" -ge 2 ]; then
    prev_sig="$(sig_snapshot_path "$((pass - 1))")"
    cur_sig="$(sig_snapshot_path "${pass}")"
    if report_signature_diff_and_check_converged "${prev_sig}" "${cur_sig}"; then
      echo "[CONVERGE] signatures converged at pass ${pass}; stopping early."
      break
    else
      status=$?
      if (( status == 2 )); then
        echo "ERROR: signature diff failed." >&2
        exit 2
      fi
    fi
  fi
done

bash "${FABLE}/patch_lapack.sh"

echo "ALL DONE"
