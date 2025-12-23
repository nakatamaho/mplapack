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
EXCLUDE_BASENAMES_MANUAL=( )

# Basenames that are not needed on the C++ side
EXCLUDE_BASENAMES_UNUSED=( la_constants dsgesv zcgesv dsposv zcposv dla_gerfsx_extended dla_gbrfsx_extended dla_porfsx_extended dla_syrfsx_extended zla_gerfsx_extended zla_gbrfsx_extended zla_porfsx_extended zla_syrfsx_extended zla_herfsx_extended lsamen ilaslr ilaslc ilaclc ilaclr iladiag icmax1 xerbla xerbla_array)

# lsamen is not used?
# ilaslr is only used in slarf.
# iladiag is not used?
# ilaslc is not used?
# ilaclc is only used in clarf.
# ilaclr is only used in clarf.

# Other helper/utility routines to be excluded from this pass
EXCLUDE_BASENAMES_MISC=( dlaruv )

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

for src in *.f; do
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

export FABLE_CONVERT
parallel -j "${JOBS:-$(nproc)}" '
     echo "Converting {}"
     bash "$FABLE_CONVERT" "{}"
 ' ::: "${files[@]}"


echo "Post-processing mplapack/reference/iMlaenv.cpp with fix_iMlaenv.py"
python $HOME/mplapack/fable/fix_iMlaenv.py iMlaenv.cpp
python $HOME/mplapack/fable/fix_iMlaenv.py iMparmq.cpp

patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Cbbcsd.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Cpptrf.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Cpotf2.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Clarnv.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Ctrexc.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Clartg.cpp
patch -p3    < ~/mplapack/fable/3.9.1/lapack/patch-Rbbcsd.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Rlarrb.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Rlarnv.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Rstebz.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Rlarrk.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Claic1.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Csyequb.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Cgesvj.cpp
patch -p3    < ~/mplapack/fable/3.9.1/lapack/patch-Rgesvj.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Rlacon.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Clacon.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Rlasy2.cpp
patch -p3    < ~/mplapack/fable/3.9.1/lapack/patch-Cgejsv.cpp
patch -p3    < ~/mplapack/fable/3.9.1/lapack/patch-Rgejsv.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Rlaln2.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-iMieeeck.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-Rladiv.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-iMlaenv.cpp
patch -p3 -R < ~/mplapack/fable/3.9.1/lapack/patch-iMparmq.cpp
