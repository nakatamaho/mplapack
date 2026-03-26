#!/usr/bin/env bash
# go_testing.sh: multi-pass Fortran->C++ conversion for LAPACK TESTING sources.
#
# Scope:
#   - Convert LAPACK TESTING/{EIG,LIN,MATGEN} Fortran sources to C++.
#   - Generate include headers using:
#       * gen_include_mplapack_eig.sh
#       * gen_include_mplapack_lin.sh
#       * gen_include_mplapack_matgen.sh
#   - BLAS conversion is intentionally skipped.

set -euo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${ROOT:-$(cd -- "${SCRIPT_DIR}/.." && pwd)}"

# ------------------------------------------------------------
# Exception lists: sources that should NOT be converted in TESTING pass
# ------------------------------------------------------------
# Exclude any file whose basename starts with these prefixes (case-insensitive).
# We skip single precision real/complex families: s* and c*.
EXCLUDE_PREFIXES=( s c )
# Exclude specific files by name (case-insensitive).
EXCLUDE_FILES=( dlaran xerbla ilaenv xlaenv zchkab dchkab zdrvab zdrvac zdrvgbx zdrvgex zdrvhex zdrvpox zdrvsyx zebchvxx zerrab zerrac zdrvgbx zdrvgex zdrvhex zdrvpox zerrgex zerrhex zerrpox zerrsyx zerrvxx ddrvab ddrvac ddrvgbx ddrvgex ddrvpox ddrvsyx debchvxx derrab derrac derrgex derrpox derrsyx derrvxx )

# Force-include specific utility routines even if they match excluded prefixes.
# Example: chkxer is required by many LAPACK tests but starts with "c".
FORCE_INCLUDE_STEMS=( chkxer )

# Hand-maintained files to keep when cleaning old generated outputs.
KEEP_HAND_WRITTEN_FILES=(
  Mxerbla.cpp
  Mxlaenv.cpp
  iMlaenv.cpp
  Rlaran.cpp
  mplapack_eig.h.in
  mplapack_eig_binary128.h.in
  mplapack_eig_binary80.h.in
  mplapack_eig_dd.h.in
  mplapack_eig_double.h.in
  mplapack_eig_gmp.h.in
  mplapack_eig_mpfr.h.in
  mplapack_eig_qd.h.in
  mplapack_lin.h.in
  mplapack_lin_binary128.h.in
  mplapack_lin_binary80.h.in
  mplapack_lin_dd.h.in
  mplapack_lin_double.h.in
  mplapack_lin_gmp.h.in
  mplapack_lin_mpfr.h.in
  mplapack_lin_qd.h.in
  mplapack_matgen.h.in
  mplapack_matgen_binary128.h.in
  mplapack_matgen_binary80.h.in
  mplapack_matgen_dd.h.in
  mplapack_matgen_double.h.in
  mplapack_matgen_gmp.h.in
  mplapack_matgen_mpfr.h.in
  mplapack_matgen_qd.h.in
  Makefile.am
  Makefile.in
  Makefile
)

GENERATED_CLEAN_DIRS=(
  "${ROOT}/mplapack/test/lin/common"
  "${ROOT}/mplapack/test/eig/common"
  "${ROOT}/mplapack/test/matgen"
)

for target_dir in "${GENERATED_CLEAN_DIRS[@]}"; do
  [ -d "${target_dir}" ] || continue

  find_args=(
    "${target_dir}"
    -maxdepth 1
    -type f
    "("
      -name '*'
    ")"
  )

  for keep in "${KEEP_HAND_WRITTEN_FILES[@]}"; do
    find_args+=( ! -name "${keep}" )
  done

  find "${find_args[@]}" -delete
done

echo "MPLAPACK testdirs have been cleaned up"
ls "${ROOT}/mplapack/test/lin/common"
ls "${ROOT}/mplapack/test/eig/common"
ls "${ROOT}/mplapack/test/matgen"

PASSES="${1:-2}"

# Default LAPACK version
LAPACK_VERSION="${LAPACK_VERSION:-3.12.1}"

export PYTHONPATH="${ROOT}:${PYTHONPATH:-}"

FABLE="${ROOT}/fable"
LAPACK_ROOT="${ROOT}/external/lapack/work/internal/lapack-${LAPACK_VERSION}"
TESTING_ROOT="${LAPACK_ROOT}/TESTING"

EIG_SRC="${TESTING_ROOT}/EIG"
LIN_SRC="${TESTING_ROOT}/LIN"
MATGEN_SRC="${TESTING_ROOT}/MATGEN"

# Where converted C++ files should be placed in the MPLAPACK tree.
# Adjust these if your repository uses a different layout.
MPLAPACK_TEST_ROOT="${MPLAPACK_TEST_ROOT:-${ROOT}/mplapack/test}"
EIG_DST="${MPLAPACK_TEST_ROOT}/eig/common/"
LIN_DST="${MPLAPACK_TEST_ROOT}/lin/common/"
MATGEN_DST="${MPLAPACK_TEST_ROOT}/matgen/"

# Generated TESTING headers (produced by gen_include_mplapack_{eig,lin,matgen}.sh).
# Allow override via environment variables, but provide safe defaults.
EIG_HDR="${EIG_HDR:-${EIG_DST%/}/mplapack_eig_generic.h}"
LIN_HDR="${LIN_HDR:-${LIN_DST%/}/mplapack_lin_generic.h}"
MATGEN_HDR="${MATGEN_HDR:-${MATGEN_DST%/}/mplapack_matgen_generic.h}"

# Optional: existing headers from your main BLAS/LAPACK conversion.
# These are used only to enrich the generated signatures for the next pass.
MPBLAS_HDR_DEFAULT="${MPBLAS_HDR_DEFAULT:-${ROOT}/mpblas/reference/mpblas_generic.h}"
MPLAPACK_HDR_DEFAULT="${MPLAPACK_HDR_DEFAULT:-${ROOT}/mplapack/reference/mplapack_generic.h}"

SIG_PY="${FABLE}/mplapack_signatures.py"
JOBS="${JOBS:-$(nproc)}"

die() {
  echo "ERROR: $*" >&2
  exit 1
}

ensure_dir() {
  local d="$1"
  /bin/mkdir -p "${d}"
}

move_cpp_or_die() {
  local src_dir="$1"
  local dst_dir="$2"
  local files=("${src_dir}"/*.cpp)

  if (( ${#files[@]} == 0 )); then
    die "No .cpp files produced in: ${src_dir}"
  fi

  ensure_dir "${dst_dir}"
  /bin/mv -f "${files[@]}" "${dst_dir}/"
}

sig_snapshot_path() {
  local pass="$1"
  echo "${SIG_PY%.py}.pass${pass}.py"
}

cleanup_signature_snapshots() {
  local snaps=( "${SIG_PY%.py}.pass"*.py )
  if (( ${#snaps[@]} )); then
    /bin/rm -f "${snaps[@]}"
  fi
}

snapshot_signatures() {
  local pass="$1"
  local dst
  dst="$(sig_snapshot_path "${pass}")"
  /bin/cp -f "${SIG_PY}" "${dst}"
  echo "[SIG] snapshot: ${dst}"
}

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

convert_dir() {
  local src_dir="$1"
  local dst_dir="$2"
  local mode="${3-}"

  ensure_dir "${dst_dir}"
  /bin/rm -f "${src_dir}"/*.cpp

  (
    cd "${src_dir}"
    # Collect Fortran sources, then filter out unnecessary precisions.
    local all=( *.f *.F *.f90 )
    if (( ${#all[@]} == 0 )); then
      die "No Fortran sources found in: ${src_dir}"
    fi

    local files=()
    local src base base_lower stem stem_lower skip force inc inc_lower pfx ex ex_lower ex_stem
    for src in "${all[@]}"; do
      base="${src##*/}"      # e.g., cbdt01.f
      base_lower="${base,,}"
      stem="${base%%.*}"     # e.g., cbdt01
      stem_lower="${stem,,}" # case-insensitive compare

      force=false
      for inc in "${FORCE_INCLUDE_STEMS[@]}"; do
        inc_lower="${inc,,}"
        if [[ "${stem_lower}" == "${inc_lower}" ]]; then
          force=true
          break
        fi
      done

      skip=false
      if ! $force; then
        for pfx in "${EXCLUDE_PREFIXES[@]}"; do
          if [[ "${stem_lower}" == "${pfx}"* ]]; then
            skip=true
            break
          fi
        done
      fi

      # Check exact-file exceptions
      if ! $skip && ! $force && (( ${#EXCLUDE_FILES[@]} > 0 )); then
        for ex in "${EXCLUDE_FILES[@]}"; do
          ex_lower="${ex,,}"
          ex_stem="${ex_lower%%.*}"  # allow both "foo.f" and "foo"
          if [[ "${base_lower}" == "${ex_lower}" || "${stem_lower}" == "${ex_stem}" ]]; then
            skip=true
            break
          fi
        done
      fi

      $skip && continue

      files+=( "${src}" )
    done

    if (( ${#files[@]} == 0 )); then
      die "After filtering prefixes (${EXCLUDE_PREFIXES[*]}), no Fortran sources remain in: ${src_dir}"
    fi
    if command -v parallel >/dev/null 2>&1; then
      parallel -j "${JOBS}" --halt soon,fail=1 \
        "echo 'Converting {}'; bash '${FABLE}/convert_lapack.sh' '{}' '${mode}'" \
        ::: "${files[@]}"
    else
      # Fallback without GNU parallel.
      for f in "${files[@]}"; do
        echo "Converting ${f}"
        bash "${FABLE}/convert_lapack.sh" "${f}" "${mode}"
      done
    fi
  )

  move_cpp_or_die "${src_dir}" "${dst_dir}"
}

generate_includes() {
  bash "${FABLE}/gen_include_mplapack_matgen.sh"
  bash "${FABLE}/gen_include_mplapack_lin.sh"
  bash "${FABLE}/gen_include_mplapack_eig.sh"
}

generate_signatures() {
  # Nounset-safe defaults (EIG_HDR etc may be unset under "set -u").
  local eig_hdr="${EIG_HDR:-${EIG_DST%/}/mplapack_eig_generic.h}"
  local lin_hdr="${LIN_HDR:-${LIN_DST%/}/mplapack_lin_generic.h}"
  local matgen_hdr="${MATGEN_HDR:-${MATGEN_DST%/}/mplapack_matgen_generic.h}"

  # Precedence: earlier < later. Later headers override earlier ones on name collisions.
  local -a inputs=(
    "${MPBLAS_HDR_DEFAULT}"
    "${MPLAPACK_HDR_DEFAULT}"
    "${eig_hdr}"
    "${lin_hdr}"
    "${matgen_hdr}"
  )

  # De-duplicate by path while preserving order.
  local -A seen=()
  local -a headers=()
  local h
  for h in "${inputs[@]}"; do
    [[ -n "${h}" ]] || die "Header path is empty (check *_HDR variables)"
    [[ -f "${h}" ]] || die "Missing header for signature generation: ${h}"
    if [[ -z "${seen[${h}]+x}" ]]; then
      seen["${h}"]=1
      headers+=( "${h}" )
    fi
  done

  local tmp
  if command -v mktemp >/dev/null 2>&1; then
    tmp="$(mktemp "${SIG_PY}.tmp.XXXXXX")"
  else
    tmp="${SIG_PY}.tmp.$$"
  fi

  python3 "${FABLE}/gen_mplapack_signatures.py" "${headers[@]}" > "${tmp}"
  /bin/mv -f "${tmp}" "${SIG_PY}"

  echo "[SIG] updated: ${SIG_PY} (inputs: ${#headers[@]} headers)"
}

run_one_pass() {
  local pass="$1"
  echo "=== PASS ${pass}/${PASSES} ==="

  convert_dir "${MATGEN_SRC}" "${MATGEN_DST}" "matgen"
  convert_dir "${LIN_SRC}" "${LIN_DST}" "lin"
  convert_dir "${EIG_SRC}" "${EIG_DST}" "eig"

  generate_includes
  generate_signatures

  echo "PASS ${pass} done. signatures: ${SIG_PY}"
}

[[ -d "${EIG_SRC}" ]] || die "Missing: ${EIG_SRC}"
[[ -d "${LIN_SRC}" ]] || die "Missing: ${LIN_SRC}"
[[ -d "${MATGEN_SRC}" ]] || die "Missing: ${MATGEN_SRC}"

# Name map affects symbol/name rewriting during conversion.
rm -f "${FABLE}/mplapack_name_map.txt"
bash "${FABLE}/gen_mplapack_name_map.sh"

cleanup_signature_snapshots

for pass in $(seq 1 "${PASSES}"); do
  if (( pass >= 2 )); then
    export MPLAPACK_SIGNATURES_PY="${SIG_PY}"
    export FABLE_DEBUG_SIGNATURES=1
  fi

  run_one_pass "${pass}"
  snapshot_signatures "${pass}"

  if (( pass >= 2 )); then
    prev_sig="$(sig_snapshot_path "$((pass - 1))")"
    cur_sig="$(sig_snapshot_path "${pass}")"
    if report_signature_diff_and_check_converged "${prev_sig}" "${cur_sig}"; then
      echo "[CONVERGE] signatures converged at pass ${pass}; stopping early."
      break
    else
      status=$?
      if (( status == 2 )); then
        die "Signature diff failed."
      fi
    fi
  fi
done

python3 "${FABLE}/rename_routine_literals.py" --in-place --map "${FABLE}/mplapack_name_map.txt" --map "${FABLE}/mplapack_testing_name_map.txt" "${ROOT}/mplapack/test/eig"
python3 "${FABLE}/rename_routine_literals.py" --in-place --map "${FABLE}/mplapack_name_map.txt" --map "${FABLE}/mplapack_testing_name_map.txt" "${ROOT}/mplapack/test/lin"
python3 "${FABLE}/rename_routine_literals.py" --in-place --map "${FABLE}/mplapack_name_map.txt" --map "${FABLE}/mplapack_testing_name_map.txt" "${ROOT}/mplapack/test/matgen"

bash "${FABLE}/sync_test_inputs.sh"
bash "${FABLE}/patch_lapack_test_${LAPACK_VERSION}.sh"

# Regenerate explicit testing source lists in Makefile.am from the
# current set of *.cpp files in the corresponding source directory.
regen_testing_makefile_sources() {
  local makefile_am="$1"
  local source_dir="$2"
  local variable_name="$3"
  local path_prefix="$4"

  python3 - "$makefile_am" "$source_dir" "$variable_name" "$path_prefix" <<'PY'
from pathlib import Path
import sys

makefile_am = Path(sys.argv[1])
source_dir = Path(sys.argv[2])
variable_name = sys.argv[3]
path_prefix = sys.argv[4]

cpp_files = sorted(path.name for path in source_dir.glob("*.cpp"))
if not cpp_files:
    raise SystemExit(f"No .cpp files found in {source_dir}")

lines = makefile_am.read_text().splitlines(keepends=True)

start = None
prefix = f"{variable_name} = "
for i, line in enumerate(lines):
    if line.startswith(prefix):
        start = i
        break

if start is None:
    raise SystemExit(f"{variable_name} block not found in {makefile_am}")

end = start + 1
while end <= len(lines) and lines[end - 1].rstrip().endswith("\\"):
    if end == len(lines):
        raise SystemExit(f"Unterminated {variable_name} block in {makefile_am}")
    end += 1

replacement = [f"{variable_name} = \\\n"]
for i, name in enumerate(cpp_files):
    item = f"{path_prefix}{name}"
    current_suffix = " \\\n" if i != len(cpp_files) - 1 else "\n"
    replacement.append(f"{item}{current_suffix}")

new_lines = lines[:]
new_lines[start:end] = replacement

old_text = "".join(lines)
new_text = "".join(new_lines)
if new_text != old_text:
    makefile_am.write_text(new_text)
PY
}

regen_testing_makefile_sources \
  "${ROOT}/mplapack/test/eig/Makefile.am" \
  "${ROOT}/mplapack/test/eig/common" \
  "EIG_SOURCES" \
  "common/"

sed -i \
  -e '/^[[:space:]]*common\/Cchkee\.cpp \\$/d' \
  -e '/^[[:space:]]*common\/Cchkdmd\.cpp \\$/d' \
  -e '/^[[:space:]]*common\/Rchkee\.cpp \\$/d' \
  -e '/^[[:space:]]*common\/Rchkdmd\.cpp \\$/d' \
  -e '/^[[:space:]]*common\/Mxerbla\.cpp \\$/d' \
  -e '/^[[:space:]]*common\/Mxlaenv\.cpp \\$/d' \
  "${ROOT}/mplapack/test/eig/Makefile.am"

regen_testing_makefile_sources \
  "${ROOT}/mplapack/test/lin/Makefile.am" \
  "${ROOT}/mplapack/test/lin/common" \
  "LIN_SOURCES" \
  "common/"

sed -i \
  -e '/^[[:space:]]*common\/Cchkaa\.cpp \\$/d' \
  -e '/^[[:space:]]*common\/Cchkrfp\.cpp \\$/d' \
  -e '/^[[:space:]]*common\/Rchkaa\.cpp \\$/d' \
  -e '/^[[:space:]]*common\/Rchkrfp\.cpp \\$/d' \
  -e '/^[[:space:]]*common\/Mxerbla\.cpp \\$/d' \
  -e '/^[[:space:]]*common\/Mxlaenv\.cpp \\$/d' \
  "${ROOT}/mplapack/test/lin/Makefile.am"

regen_testing_makefile_sources \
  "${ROOT}/mplapack/test/matgen/Makefile.am" \
  "${ROOT}/mplapack/test/matgen" \
  "MATGEN_SOURCES" \
  ""

echo "ALL DONE"
