#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

FABLE_CONVERT="${SCRIPT_DIR}/convert_lapack.sh"

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
EXCLUDE_BASENAMES_UNUSED=( la_constants dsgesv zcgesv dsposv zcposv dla_gerfsx_extended dla_gbrfsx_extended dla_porfsx_extended dla_syrfsx_extended zla_gerfsx_extended zla_gbrfsx_extended zla_porfsx_extended zla_syrfsx_extended zla_herfsx_extended lsamen ilaslr ilaslc ilaclc ilaclr iladiag icmax1 xerbla xerbla_array zgbrfsx zgbsvxx zgerfsx zherfsx zporfsx zsyrfsx dgerfsx dgbrfsx dporfsx dsyrfsx zgesvxx zhesvxx zposvxx zsysvxx dgbsvxx dgesvxx dposvxx dsysvxx zla_gbrcond_c zla_gbrcond_x zla_gbrpvgrw zla_gercond_c zla_gercond_x zla_gerpvgrw zla_hercond_c zla_hercond_x zla_herpvgrw zla_porcond_c zla_porcond_x zla_porpvgrw zla_syrcond_c zla_syrcond_x zla_syrpvgrw zla_gbamv zla_geamv zla_heamv zla_lin_berr zla_syamv zla_wwaddw dla_gbrcond dla_gbrpvgrw dla_gercond dla_gerpvgrw dla_porcond dla_porpvgrw dla_syrcond dla_syrpvgrw dla_gbamv dla_geamv dla_lin_berr dla_syamv dla_wwaddw zlag2c zlat2c dlag2s dlat2s zlag2z zlat2c disnan dlaisnan \
dlacon zlacon dlabad
)
# lsamen is not used?
# ilaslr is only used in slarf.
# iladiag is not used?
# ilaslc is not used?
# ilaclc is only used in clarf.
# ilaclr is only used in clarf.
# The following ones use xBLAS so we do not support them.
# zgbsvxx
# zgbrfsx
# zgerfsx
# zherfsx
# zporfsx
# zsyrfsx
# dgerfsx
# dgbrfsx
# dporfsx
# dsyrfsx
# zgesvxx
# zhesvxx
# zposvxx
# zsysvxx
# dgbsvxx
# dgesvxx
# dposvxx
# dsysvxx
# zla_gbrcond_c
# zla_gbrcond_x
# zla_gbrpvgrw
# zla_gercond_c
# zla_gercond_x
# zla_gerpvgrw
# zla_hercond_c
# zla_hercond_x
# zla_herpvgrw
# zla_porcond_c
# zla_porcond_x
# zla_porpvgrw
# zla_syrcond_c
# zla_syrcond_x
# zla_syrpvgrw
# zla_gbamv
# zla_geamv
# zla_heamv
# zla_lin_berr
# zla_syamv
# zla_wwaddw
# dla_gbrcond
# dla_gbrpvgrw
# dla_gercond
# dla_gerpvgrw
# dla_porcond
# dla_porpvgrw
# dla_syrcond
# dla_syrpvgrw
# dla_gbamv
# dla_geamv
# dla_lin_berr
# dla_syamv
# dla_wwaddw
# The following ones converts to a lower precision, so we dont need them.
# clag2c
# clat2c
# dlag2s
# dlat2s
# deprecated
# dlacon
# zlacon

# Other helper/utility routines to be excluded from this pass
EXCLUDE_BASENAMES_MISC=( dlaruv )

# Combined basename exception list (for iteration)
EXCLUDE_BASENAMES=(
  "${EXCLUDE_BASENAMES_MANUAL[@]}"
  "${EXCLUDE_BASENAMES_UNUSED[@]}"
  "${EXCLUDE_BASENAMES_MISC[@]}"
)

# Hand-maintained files to keep when cleaning old generated outputs.
KEEP_HAND_WRITTEN_FILES=(
  Mexponent.cpp
  Mmaxloc.cpp
  Mlsamen.cpp
  Mmaxval.cpp
  Mminval.cpp
  Mxerbla.cpp
  Mxlaenv.cpp
  Risinf.cpp
  Risnan.cpp
  Rlabad.cpp
  Rlamch.cpp
  Rlaruv.cpp
  Rlaran.cpp
  Rroundup_lwork.cpp
  iMlaenv.cpp
  iMladiag.cpp
  iMlaver.cpp
  iMlaver.cpp.in
  Makefile
  Makefile.am
  Makefile.in
  mplapack.h.in
  mplapack_binary128.h.in
  mplapack_binary80.h.in
  mplapack_dd.h.in
  mplapack_double.h.in
  mplapack_gmp.h.in
  mplapack_mpfr.h.in
  mplapack_qd.h.in
)

# ------------------------------------------------------------
# Collect all Fortran files (*.f, *.f90) except those in exceptions
# ------------------------------------------------------------

files=()

for src in *.f *.F *.f90; do
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

#for src in "${files[@]}"; do
#    echo "Converting $src"
#    bash "$FABLE_CONVERT" "$src"
#done

GENERATED_CLEAN_DIRS=(
  "${ROOT}/mplapack/reference"
)

# Keep the hand-maintained iMlaenv.cpp in common test outputs, but remove the
# stale generated copy left in mplapack/reference from the previous output layout.
rm -f "${ROOT}/mplapack/reference/iMlaenv.cpp"

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

echo "MPLAPACK dir has been cleaned up"
ls "${ROOT}/mplapack/reference"

export FABLE_CONVERT
parallel -j "${JOBS:-$(nproc)}" '
     echo "Converting {}"
     bash "$FABLE_CONVERT" "{}"
 ' ::: "${files[@]}"
