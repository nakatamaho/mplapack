#!/usr/bin/env bash
set -euo pipefail

# Path to the FABLE Fortran->C++ converter
FABLE_CONVERT="$HOME/mplapack/fable/convert_blas.sh"

# Enable nullglob to avoid literal "*.f*" when no files exist
shopt -s nullglob

# ------------------------------------------------------------
# Exception lists: files that should NOT be converted
# ------------------------------------------------------------

# Exclude any file whose basename starts with these prefixes
EXCLUDE_PREFIXES=( s c )

# Basenames that will be converted manually (hand-written C++)
EXCLUDE_BASENAMES_MANUAL=( dznrm2 dnrm2 drotg )

# Basenames that are not needed on the C++ side
EXCLUDE_BASENAMES_UNUSED=( icamax isamax )

# Other helper/utility routines to be excluded from this pass
EXCLUDE_BASENAMES_MISC=( lsame xerbla )

# Combined basename exception list (for iteration)
EXCLUDE_BASENAMES=(
  "${EXCLUDE_BASENAMES_MANUAL[@]}"
  "${EXCLUDE_BASENAMES_UNUSED[@]}"
  "${EXCLUDE_BASENAMES_MISC[@]}"
)

# ------------------------------------------------------------
# Special output name mappings
#   key   = Fortran basename (without extension)
#   value = C++ output filename
# ------------------------------------------------------------
declare -A SPECIAL_NAME_MAP=(
  [idamax]="iRamax.cpp"
  [izamax]="iCamax.cpp"
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
# ------------------------------------------------------------

for src in "${files[@]}"; do
    base="${src##*/}"
    stem="${base%%.*}"

    # Special naming cases first
    if [[ -n "${SPECIAL_NAME_MAP[$stem]+set}" ]]; then
        out="${SPECIAL_NAME_MAP[$stem]}"
        echo "Converting $src -> $out (special naming)"
        bash "$FABLE_CONVERT" "$src" > "$out"
        continue
    fi

    first="${stem:0:1}"
    rest="${stem:1}"

    # Determine output prefix based on BLAS naming convention:
    #   d/s -> R   (real)
    #   z/c -> C   (complex)
    #   others -> capitalize first letter
    case "$first" in
        d|D|s|S)
            prefix="R"
            ;;
        z|Z|c|C)
            prefix="C"
            ;;
        *)
            prefix=$(printf '%s' "$first" | tr '[:lower:]' '[:upper:]')
            ;;
    esac

    out="${prefix}${rest}.cpp"

    echo "Converting $src -> $out"
    bash "$FABLE_CONVERT" "$src" > "$out"
done
