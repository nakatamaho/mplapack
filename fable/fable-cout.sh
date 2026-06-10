#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "${script_dir}/.." && pwd)"

lapack_version="${LAPACK_VERSION:-3.12.1}"
default_file="${repo_root}/external/lapack/work/internal/lapack-${lapack_version}/BLAS/SRC/zrotg.f90"

show_help() {
    cat <<'EOF'
Usage:
  fable-cout.sh [FORTRAN_FILE]

Examples:
  fable-cout.sh
  fable-cout.sh external/lapack/work/internal/lapack-3.12.1/BLAS/SRC/drotg.f90
  LAPACK_VERSION=3.9.1 fable-cout.sh
EOF
}

case "${1:-}" in
    -h|--help)
        show_help
        exit 0
        ;;
esac

target="${1:-$default_file}"

if [[ ! "$target" = /* ]]; then
    target="${repo_root}/${target}"
fi

export PYTHONPATH="${repo_root}${PYTHONPATH:+:${PYTHONPATH}}"
exec python -m fable.command_line.cout "$target"
