#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
LAPACK_VERSION="${LAPACK_VERSION:-3.12.1}"

# Path to the FABLE Fortran->C++ converter
FABLE_CONVERT="${SCRIPT_DIR}/convert_blas.sh"

# Enable nullglob to avoid literal "*.f*" when no files exist
shopt -s nullglob

# ------------------------------------------------------------
# Exception lists: files that should NOT be converted
# ------------------------------------------------------------

# Exclude any file whose basename starts with these prefixes
EXCLUDE_PREFIXES=( s c )

# Basenames that will be converted manually
EXCLUDE_BASENAMES_MANUAL=( )

# Basenames that are not needed on the C++ side
EXCLUDE_BASENAMES_UNUSED=( icamax isamax )

# Other helper/utility routines to be excluded from this pass
EXCLUDE_BASENAMES_MISC=( xerbla lsame xerbla_array dsdot )

# Combined basename exception list (for iteration)
EXCLUDE_BASENAMES=(
  "${EXCLUDE_BASENAMES_MANUAL[@]}"
  "${EXCLUDE_BASENAMES_UNUSED[@]}"
  "${EXCLUDE_BASENAMES_MISC[@]}"
)

# Hand-maintained files to keep when cleaning old generated outputs.
KEEP_HAND_WRITTEN_FILES=(
  Mxerbla.cpp
  Mlsame.cpp
  mplapackinit.cpp
  mpblas.h.in
  mpblas_binary128.h.in
  mpblas_binary80.h.in
  mpblas_dd.h.in
  mpblas_double.h.in
  mpblas_gmp.h.in
  mpblas_mpfr.h.in
  mpblas_qd.h.in
  Makefile.am
  Makefile.in
  Makefile
)

# ------------------------------------------------------------
# Collect all Fortran files (*.f, *.f90) except those in exceptions
# ------------------------------------------------------------

files=()

for src in *.f *.f90 *.F; do
    base="${src##*/}"       # e.g., zgemv.f
    stem="${base%%.*}"      # e.g., zgemv

    skip=false

    # Check prefix exceptions (e.g., s*, c*)
    for pfx in "${EXCLUDE_PREFIXES[@]}"; do
        if [[ "$stem" == "$pfx"* ]]; then
            skip=true
            break
        fi
    done

    # Check exact-basename exceptions
    if ! $skip; then
        for ex in "${EXCLUDE_BASENAMES[@]}"; do
            if [[ "$stem" == "$ex" ]]; then
                skip=true
                break
            fi
        done
    fi

    $skip && continue

    files+=( "$src" )
done

# ------------------------------------------------------------
# Convert each selected Fortran file using convert_blas.sh
# (convert_blas.sh / cout.py decide the .cpp name and do all postprocessing)
# ------------------------------------------------------------
MPBLAS_REF="${ROOT}/mpblas/reference"

find_args=(
  "${MPBLAS_REF}"
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

echo "MPBLAS dir has been cleaned up"
ls "${MPBLAS_REF}"

export FABLE_CONVERT
parallel -j "${JOBS:-$(nproc)}" '
     echo "Converting {}"
     bash "$FABLE_CONVERT" "{}"
 ' ::: "${files[@]}"


case "${LAPACK_VERSION}" in
    3.9.1)
        patch < "${SCRIPT_DIR}/${LAPACK_VERSION}/patch-blas"
        ;;
    3.12.1)
        patch < "${SCRIPT_DIR}/${LAPACK_VERSION}/blas/patch-Crotg.cpp"
        patch < "${SCRIPT_DIR}/${LAPACK_VERSION}/blas/patch-RCnrm2.cpp"
        patch < "${SCRIPT_DIR}/${LAPACK_VERSION}/blas/patch-Rnrm2.cpp"
        patch < "${SCRIPT_DIR}/${LAPACK_VERSION}/blas/patch-Rgemm.cpp"
        patch < "${SCRIPT_DIR}/${LAPACK_VERSION}/blas/patch-Rrotg.cpp"
        patch < "${SCRIPT_DIR}/${LAPACK_VERSION}/blas/patch-Rrotmg.cpp"
        ;;
    *)
        echo "Unsupported LAPACK_VERSION: ${LAPACK_VERSION}" >&2
        exit 1
        ;;
esac
