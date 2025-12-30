#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Locate MPLAPACK root and LAPACK TESTING directories
# ------------------------------------------------------------

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mplapack_root="$(cd "${script_dir}/.." && pwd)"
lapack_version="${LAPACK_VERSION:-3.9.1}"

# Path to the FABLE Fortran->C++ converter (override via env if needed)
FABLE_CONVERT="${FABLE_CONVERT:-${mplapack_root}/fable/convert_lapack.sh}"

testing_root="${mplapack_root}/external/lapack/work/internal/lapack-${lapack_version}/TESTING"
default_dirs=(
  "${testing_root}/EIG"
  "${testing_root}/LIN"
  "${testing_root}/MATGEN"
)

# Allow overriding target directories by command-line arguments
if [[ "$#" -gt 0 ]]; then
  target_dirs=( "$@" )
else
  target_dirs=( "${default_dirs[@]}" )
fi

for d in "${target_dirs[@]}"; do
  if [[ ! -d "$d" ]]; then
    echo "Error: directory not found: $d" >&2
    exit 1
  fi
done

# ------------------------------------------------------------
# Exception lists: files that should NOT be converted
# ------------------------------------------------------------

# Exclude any file whose basename starts with these prefixes
EXCLUDE_PREFIXES=( s c )

# Basenames that will be converted manually
EXCLUDE_BASENAMES_MANUAL=( )

# Basenames that are not needed on the C++ side
EXCLUDE_BASENAMES_UNUSED=( la_constants )

# Other helper/utility routines to be excluded from this pass
EXCLUDE_BASENAMES_MISC=( dlaruv )

# Combined basename exception list (for iteration)
EXCLUDE_BASENAMES=(
  "${EXCLUDE_BASENAMES_MANUAL[@]}"
  "${EXCLUDE_BASENAMES_UNUSED[@]}"
  "${EXCLUDE_BASENAMES_MISC[@]}"
)

# ------------------------------------------------------------
# Collect all Fortran files (*.f, *.F, *.f90, *.F90) except those in exceptions
# ------------------------------------------------------------

files=()

mapfile -t candidates < <(
  find "${target_dirs[@]}" -maxdepth 1 -type f \( -iname '*.f' -o -iname '*.f90' \) | sort
)

for src in "${candidates[@]}"; do
    base="$(basename "$src")"     # e.g., zlatmr.f or ZLATMR.F
    stem="${base%%.*}"            # e.g., zlatmr or ZLATMR
    stem_lc="${stem,,}"           # normalize for comparisons

    skip=false

    # Check prefix exceptions (e.g., s*, c*)
    for pfx in "${EXCLUDE_PREFIXES[@]}"; do
        pfx_lc="${pfx,,}"
        if [[ "$stem_lc" == "$pfx_lc"* ]]; then
            skip=true
            break
        fi
    done

    # Check exact-basename exceptions
    if ! $skip; then
        for ex in "${EXCLUDE_BASENAMES[@]}"; do
            ex_lc="${ex,,}"
            if [[ "$stem_lc" == "$ex_lc" ]]; then
                skip=true
                break
            fi
        done
    fi

    $skip && continue

    files+=( "$src" )
done

#for src in "${files[@]}"; do
#    echo "Converting $src"
#    bash "$FABLE_CONVERT" "$src"
#done

export FABLE_CONVERT
parallel -j "${JOBS:-$(nproc)}" '
     echo "Converting {}"
     # Decide mode from the parent directory name (EIG/LIN/MATGEN)
     parent_dir="$(basename "$(dirname "{}")")"
     mode=""
     case "${parent_dir,,}" in
         lin)    mode="lin" ;;
         eig)    mode="eig" ;;
         matgen) mode="matgen" ;;
         *)      mode="" ;;
     esac

     if [ -n "$mode" ]; then
         bash "$FABLE_CONVERT" "{}" "$mode"
     else
         bash "$FABLE_CONVERT" "{}"
     fi
  ' ::: "${files[@]}"