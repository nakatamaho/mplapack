#!/usr/bin/env bash
set -euo pipefail

# Path to the FABLE Fortran->C++ converter
FABLE_CONVERT="$HOME/mplapack/fable/convert_lapack.sh"

# Enable nullglob to avoid literal "*.f*" when no files exist
shopt -s nullglob

# ------------------------------------------------------------
# Exception lists: files that should NOT be converted
# ------------------------------------------------------------

# Exclude any file whose basename starts with these prefixes
EXCLUDE_PREFIXES=( s c )

# Basenames that will be converted manually
EXCLUDE_BASENAMES_MANUAL=( dgbcon dgbrfs dgbtrf dgecon dgedmd dgedmdq dgees dgedmdq)

# Basenames that are not needed on the C++ side
EXCLUDE_BASENAMES_UNUSED=( )

# Other helper/utility routines to be excluded from this pass
EXCLUDE_BASENAMES_MISC=( )

# Combined basename exception list (for iteration)
EXCLUDE_BASENAMES=(
  "${EXCLUDE_BASENAMES_MANUAL[@]}"
  "${EXCLUDE_BASENAMES_UNUSED[@]}"
  "${EXCLUDE_BASENAMES_MISC[@]}"
)

# ------------------------------------------------------------
# Collect all Fortran files (*.f, *.f90) except those in exceptions
# ------------------------------------------------------------

files=()

for src in *.f*; do
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

: "${FABLE_CONVERT:?FABLE_CONVERT is not set}"

for src in "${files[@]}"; do
    echo "Converting $src"
    if ! bash "$FABLE_CONVERT" "$src"; then
        echo "Warning: conversion failed for $src (skipped)" >&2
        continue
    fi
done

