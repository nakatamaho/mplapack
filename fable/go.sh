#!/usr/bin/env bash
# two_pass_convert.sh
set -euo pipefail
shopt -s nullglob

PASSES="${1:-3}"

ROOT="/home/docker/mplapack"
# Ensure 'import fable.cout' works from this script.
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

debug_show_fable_cout() {
  echo "[DEBUG] python3=$(command -v python3)"
  python3 - <<'PY'
import fable.cout
print("[DEBUG] fable.cout.__file__ =", fable.cout.__file__)
print("[DEBUG] has _load_mplapack_signatures =", hasattr(fable.cout, "_load_mplapack_signatures"))
print("[DEBUG] has _lookup_routine_signature =", hasattr(fable.cout, "_lookup_routine_signature"))
print("[DEBUG] has _propagate_out_inout_through_calls =", hasattr(fable.cout, "_propagate_out_inout_through_calls"))
PY
}

debug_check_signatures_loaded() {
  python3 - <<'PY'
import fable.cout as c
import mplapack_signatures as s
print("[DEBUG] mplapack_signatures.__file__ =", getattr(s, "__file__", None))
fs = getattr(c, "FUNCTION_SIGNATURES", {}) or {}
print("[DEBUG] FUNCTION_SIGNATURES size =", len(fs))
print("[DEBUG] clatdf sig =", fs.get("clatdf"))
PY
}

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

  # Generate signatures so the next pass can use them if the pipeline supports it.
  python3 "${FABLE}/gen_mplapack_signatures.py" "${MPBLAS_HDR}" "${MPLAPACK_HDR}" > "${SIG_PY}"

  echo "PASS ${pass} done. signatures: ${SIG_PY}"
}

rm -f "${FABLE}/mplapack_name_map.txt"
bash "${FABLE}/gen_mplapack_name_map.sh"
rm -f "${SIG_PY}"
rm -f "${MPBLAS_REF}/mpblas_generic.h"
rm -f "${MPLAPACK_REF}/mplapack_generic.h"
for pass in $(seq 1 "${PASSES}"); do
  if [ "${pass}" -ge 2 ]; then
    export MPLAPACK_SIGNATURES_PY="${SIG_PY}"
    export FABLE_DEBUG_SIGNATURES=1
    # (2) Confirm which fable/cout.py is actually imported.
    echo "[DEBUG] MPLAPACK_SIGNATURES_PY=${MPLAPACK_SIGNATURES_PY}"
    if [[ ! -s "${SIG_PY}" ]]; then
      echo "ERROR: signatures file missing or empty: ${SIG_PY}" >&2
      exit 1
    fi
    debug_show_fable_cout

    # (3) Confirm signatures are loaded and contain CLATDF.
    debug_check_signatures_loaded
  fi
  run_one_pass "${pass}"
done

bash "${FABLE}/patch_lapack.sh"

echo "ALL DONE"
