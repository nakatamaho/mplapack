#!/bin/bash
fable_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${fable_script_dir}/clang_format_common.sh"

cd ~/mplapack/mplapack/test/matgen

if [ `uname` = "Linux" ]; then
    SED=sed
else
    SED=gsed
fi

FILES=`ls *cpp | grep -v Rlamch | grep -v iMlaenv | grep -v Mutils`
for filename in $FILES; do
ctags -x --c++-kinds=pf --language-force=c++ --_xformat='%{typeref} %{name} %{signature};' "${filename}" | sed -E 's/^typename[[:space:]]*:[[:space:]]*//' > "${filename%.*}.hpp"
done

printf "REAL Rlamch(const char *cmach);" > Rlamch.hpp
printf "INTEGER iMlaenv2stage(INTEGER const ispec, const char *name, const char *opts, INTEGER const n1, INTEGER const n2, INTEGER const n3, INTEGER const n4);" > iMlaenv2stage.hpp
printf "INTEGER iMlaenv(INTEGER const ispec, const char *name, const char *opts, INTEGER const n1, INTEGER const n2, INTEGER const n3, INTEGER const n4);" > iMlaenv.hpp

cat *hpp \
  | grep -v abs1 \
  | grep -vE '^[[:space:]]*-[[:space:]]+' \
  | grep -v main \
  | LC_ALL=C sort | uniq > header_all

rm *hpp

MPLIBS="gmp mpfr binary128 dd qd double binary80"
for mplib in $MPLIBS; do
    if [ x"$mplib" = x"gmp" ]; then
        cp header_all mplapack_matgen_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_matgen_${mplib}.h 
        sed -i -e 's/COMPLEX/mpc_class/g' mplapack_matgen_${mplib}.h 
        sed -i -e 's/REAL/mpf_class/g' mplapack_matgen_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_matgen_${mplib}.h 
        sed -i -e "s/iMlaenv/iMlaenv_${mplib}/g" mplapack_matgen_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_matgen_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_matgen_${mplib}.h 
    fi

    if [ x"$mplib" = x"mpfr" ]; then
        cp header_all mplapack_matgen_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_matgen_${mplib}.h 
        sed -i -e 's/COMPLEX/mpcomplex/g' mplapack_matgen_${mplib}.h 
        sed -i -e 's/REAL/mpreal/g' mplapack_matgen_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_matgen_${mplib}.h 
        sed -i -e "s/iMlaenv/iMlaenv_${mplib}/g" mplapack_matgen_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_matgen_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_matgen_${mplib}.h 
    fi

    if [ x"$mplib" = x"double" ]; then
        cp header_all mplapack_matgen_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_matgen_${mplib}.h 
        sed -i -e 's/COMPLEX/std::complex<double>/g' mplapack_matgen_${mplib}.h 
        sed -i -e 's/REAL/double/g' mplapack_matgen_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_matgen_${mplib}.h 
        sed -i -e "s/iMlaenv/iMlaenv_${mplib}/g" mplapack_matgen_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_matgen_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_matgen_${mplib}.h 
    fi

    if [ x"$mplib" = x"dd" ]; then
        cp header_all mplapack_matgen_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_matgen_${mplib}.h 
        sed -i -e 's/COMPLEX/dd_complex/g' mplapack_matgen_${mplib}.h 
        sed -i -e 's/REAL/dd_real/g' mplapack_matgen_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_matgen_${mplib}.h 
        sed -i -e "s/iMlaenv/iMlaenv_${mplib}/g" mplapack_matgen_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_matgen_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_matgen_${mplib}.h 
    fi

    if [ x"$mplib" = x"qd" ]; then
        cp header_all mplapack_matgen_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_matgen_${mplib}.h 
        sed -i -e 's/COMPLEX/qd_complex/g' mplapack_matgen_${mplib}.h 
        sed -i -e 's/REAL/qd_real/g' mplapack_matgen_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_matgen_${mplib}.h 
        sed -i -e "s/iMlaenv/iMlaenv_${mplib}/g" mplapack_matgen_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_matgen_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_matgen_${mplib}.h 
    fi

    if [ x"$mplib" = x"binary128" ]; then
        cp header_all mplapack_matgen_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_matgen_${mplib}.h 
        sed -i -e 's/COMPLEX/std::complex<mplapack_binary128_t>/g' mplapack_matgen_${mplib}.h
        sed -i -e 's/REAL/mplapack_binary128_t/g' mplapack_matgen_${mplib}.h
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_matgen_${mplib}.h 
        sed -i -e "s/iMlaenv/iMlaenv_${mplib}/g" mplapack_matgen_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_matgen_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_matgen_${mplib}.h 
    fi

    if [ x"$mplib" = x"binary80" ]; then
        cp header_all mplapack_matgen_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_matgen_${mplib}.h 
        sed -i -e 's/COMPLEX/std::complex<mplapack_binary80_t>/g' mplapack_matgen_${mplib}.h
        sed -i -e 's/REAL/mplapack_binary80_t/g' mplapack_matgen_${mplib}.h
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_matgen_${mplib}.h 
        sed -i -e "s/iMlaenv/iMlaenv_${mplib}/g" mplapack_matgen_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_matgen_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_matgen_${mplib}.h 
    fi

    fable_clang_format_stdout mplapack_matgen_${mplib}.h | LC_ALL=C sort > l ; mv l mplapack_matgen_${mplib}.h 
    cat ~/mplapack/mplapack/test/matgen/mplapack_matgen_${mplib}.h.in mplapack_matgen_${mplib}.h > ~/mplapack/include/mplapack_matgen_${mplib}.h
    rm mplapack_matgen_${mplib}.h
    echo "#endif" >> ~/mplapack/include/mplapack_matgen_${mplib}.h

done

mv header_all mplapack_matgen_generic.h

for f in mplapack_matgen_generic.h; do
fable_clang_format_inplace "$f"
done
