#!/usr/bin/env bash
set -euo pipefail

# Copy LAPACK TESTING .in files into MPLAPACK test tree.
# - LIN: ztest(.in/.rfp) -> Ctest(.in/.rfp), dtest(.in/.rfp) -> Rtest(.in/.rfp)
# - EIG: parse Makefile to collect .in files used by z*.out and d*.out (EIG/*),
#        then rename z*.in -> C*.in, d*.in -> R*.in, and copy.
#
# Default source: directory where this script lives (expects Makefile and .in files there).
# Default destination root: $HOME/mplapack/mplapack

usage() {
  cat <<'USAGE' >&2
Usage:
  sync_test_inputs.sh [--src DIR] [--dest-root DIR] [--dry-run]

Options:
  --src DIR        Source directory containing Makefile and *.in (default: script dir)
  --dest-root DIR  MPLAPACK root directory (default: $HOME/mplapack/mplapack)
  --dry-run        Print actions only; do not copy
USAGE
}

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
root_dir="$(cd -- "${script_dir}/.." && pwd -P)"

lapack_version="${LAPACK_VERSION:-3.12.1}"

src_dir="${root_dir}/external/lapack/work/internal/lapack-${lapack_version}/TESTING"
dest_root="${root_dir}/mplapack"
dry_run=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --src)
      [[ $# -ge 2 ]] || { usage; exit 2; }
      src_dir="$2"
      shift 2
      ;;
    --dest-root)
      [[ $# -ge 2 ]] || { usage; exit 2; }
      dest_root="$2"
      shift 2
      ;;
    --dry-run)
      dry_run=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Error: unknown argument: $1" >&2
      usage
      exit 2
      ;;
  esac
done

makefile="${src_dir}/Makefile"
if [[ ! -f "${makefile}" ]]; then
  echo "Error: Makefile not found: ${makefile}" >&2
  exit 1
fi

lin_dest="${dest_root}/test/lin"
eig_dest="${dest_root}/test/eig"

mkdir -p "${lin_dest}" "${eig_dest}"

do_copy() {
  local src="$1"
  local dst="$2"

  if [[ ! -f "${src}" ]]; then
    echo "Error: source file not found: ${src}" >&2
    exit 1
  fi

  if [[ "${dry_run}" -eq 1 ]]; then
    echo "cp -f '${src}' '${dst}'"
  else
    cp -f "${src}" "${dst}"
    echo "copied: ${src} -> ${dst}"
  fi
}

echo "Source dir : ${src_dir}"
echo "Dest LIN   : ${lin_dest}"
echo "Dest EIG   : ${eig_dest}"
echo

# ---- 1) LIN: explicit mapping
declare -A lin_map=(
  ["ztest.in"]="Ctest.in"
  ["ztest_rfp.in"]="Ctest_rfp.in"
  ["dtest.in"]="Rtest.in"
  ["dtest_rfp.in"]="Rtest_rfp.in"
)

echo "[LIN] Copying required inputs..."
for k in "${!lin_map[@]}"; do
  do_copy "${src_dir}/${k}" "${lin_dest}/${lin_map[$k]}"
done
echo

# ---- 2) EIG: parse Makefile for z*.out / d*.out rules that run EIG/*
# Extract any token ending with .in from those rules.
echo "[EIG] Collecting *.in dependencies for z*.out and d*.out from Makefile..."

mapfile -t eig_in_files < <(
  awk '
    # Match rules like: zxxx.out: something.in EIG/xeigtstz
    #                   dxxx.out: something.in EIG/xeigtstd
    ($1 ~ /^[zd][^[:space:]]*\.out:$/) && ($0 ~ /EIG\//) {
      for (i=2; i<=NF; i++) if ($i ~ /\.in$/) print $i
    }
  ' "${makefile}" | sort -u
)

if [[ "${#eig_in_files[@]}" -eq 0 ]]; then
  echo "Error: no EIG .in dependencies found in ${makefile}" >&2
  exit 1
fi

echo "[EIG] Copying ${#eig_in_files[@]} input files..."
for in_name in "${eig_in_files[@]}"; do
  base="$(basename "${in_name}")"
  case "${base}" in
    z*.in)
      dst_name="C${base:1}"
      ;;
    d*.in)
      dst_name="R${base:1}"
      ;;
    *)
      dst_name="${base}"
      ;;
  esac
  do_copy "${src_dir}/${base}" "${eig_dest}/${dst_name}"
done

echo
echo "Done."
