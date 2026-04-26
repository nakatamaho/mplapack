#!/bin/bash
fable_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${fable_script_dir}/clang_format_common.sh"

cd ~/mplapack/mplapack/test/eig/common

if [ `uname` = "Linux" ]; then
    SED=sed
else
    SED=gsed
fi

FILES=`ls *cpp | grep -v Rlamch | grep -v Mlaenv | grep -v Mutils`
for filename in $FILES; do
ctags -x --c++-kinds=pf --language-force=c++ --_xformat='%{typeref} %{name} %{signature};' "${filename}" | sed -E 's/^typename[[:space:]]*:[[:space:]]*//' > "${filename%.*}.hpp"
done

printf "REAL Rlamch(const char *cmach);\n" > Rlamch.hpp
printf "INTEGER iMlaenv2stage(INTEGER const ispec, const char *name, const char *opts, INTEGER const n1, INTEGER const n2, INTEGER const n3, INTEGER const n4);\n" > iMlaenv2stage.hpp
printf "INTEGER iMlaenv(INTEGER const ispec, const char *name, const char *opts, INTEGER const n1, INTEGER const n2, INTEGER const n3, INTEGER const n4);\n" > iMlaenv.hpp

awk '1' ./*.hpp \
  | grep -v abs1 \
  | grep -vE '^[[:space:]]*-' \
  | grep -vF 'common(int argc, char const *argv[]);' \
  | grep -vF 'common(int argc, const char *argv[]);' \
  | grep -vF 'program_' \
  | grep -v main \
  | LC_ALL=C sort -u > header_all

rm ./*.hpp

MPLIBS="gmp mpfr binary128 dd qd double binary80"
for mplib in $MPLIBS; do
    if [ x"$mplib" = x"gmp" ]; then
        cat header_all | grep -v mpfr > mplapack_eig_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_eig_${mplib}.h 
        sed -i -e 's/COMPLEX/mpc_class/g' mplapack_eig_${mplib}.h 
        sed -i -e 's/REAL/mpf_class/g' mplapack_eig_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/Mxerbla/Mxerbla_${mplib}/g" mplapack_eig_${mplib}.h
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_eig_${mplib}.h
        sed -i -e "s/\<iMlaenv2stage\>/iMlaenv2stage_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/\<iMlaenv\>/iMlaenv_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_eig_${mplib}.h 
    fi

    if [ x"$mplib" = x"mpfr" ]; then
        cat header_all | grep -v gmp > mplapack_eig_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_eig_${mplib}.h 
        sed -i -e 's/COMPLEX/mpcomplex/g' mplapack_eig_${mplib}.h 
        sed -i -e 's/REAL/mpreal/g' mplapack_eig_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/Mxerbla/Mxerbla_${mplib}/g" mplapack_eig_${mplib}.h
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_eig_${mplib}.h
        sed -i -e "s/\<iMlaenv2stage\>/iMlaenv2stage_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/\<iMlaenv\>/iMlaenv_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_eig_${mplib}.h 
    fi

    if [ x"$mplib" = x"double" ]; then
        cat header_all | grep -v gmp | grep -v mpfr > mplapack_eig_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_eig_${mplib}.h 
        sed -i -e 's/COMPLEX/std::complex<double>/g' mplapack_eig_${mplib}.h 
        sed -i -e 's/REAL/double/g' mplapack_eig_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/Mxerbla/Mxerbla_${mplib}/g" mplapack_eig_${mplib}.h
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_eig_${mplib}.h
        sed -i -e "s/\<iMlaenv2stage\>/iMlaenv2stage_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/\<iMlaenv\>/iMlaenv_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_eig_${mplib}.h 
    fi

    if [ x"$mplib" = x"dd" ]; then
        cat header_all | grep -v gmp | grep -v mpfr > mplapack_eig_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_eig_${mplib}.h 
        sed -i -e 's/COMPLEX/dd_complex/g' mplapack_eig_${mplib}.h 
        sed -i -e 's/REAL/dd_real/g' mplapack_eig_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/Mxerbla/Mxerbla_${mplib}/g" mplapack_eig_${mplib}.h
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_eig_${mplib}.h
        sed -i -e "s/\<iMlaenv2stage\>/iMlaenv2stage_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/\<iMlaenv\>/iMlaenv_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_eig_${mplib}.h 
    fi

    if [ x"$mplib" = x"qd" ]; then
        cat header_all | grep -v gmp | grep -v mpfr > mplapack_eig_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_eig_${mplib}.h 
        sed -i -e 's/COMPLEX/qd_complex/g' mplapack_eig_${mplib}.h 
        sed -i -e 's/REAL/qd_real/g' mplapack_eig_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/Mxerbla/Mxerbla_${mplib}/g" mplapack_eig_${mplib}.h
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_eig_${mplib}.h
        sed -i -e "s/\<iMlaenv2stage\>/iMlaenv2stage_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/\<iMlaenv\>/iMlaenv_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_eig_${mplib}.h 
    fi

    if [ x"$mplib" = x"binary128" ]; then
        cat header_all | grep -v gmp | grep -v mpfr > mplapack_eig_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_eig_${mplib}.h 
        sed -i -e 's/COMPLEX/std::complex<mplapack_binary128_t>/g' mplapack_eig_${mplib}.h
        sed -i -e 's/REAL/mplapack_binary128_t/g' mplapack_eig_${mplib}.h
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/Mxerbla/Mxerbla_${mplib}/g" mplapack_eig_${mplib}.h
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_eig_${mplib}.h
        sed -i -e "s/\<iMlaenv2stage\>/iMlaenv2stage_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/\<iMlaenv\>/iMlaenv_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_eig_${mplib}.h 
    fi

    if [ x"$mplib" = x"binary80" ]; then
        cat header_all | grep -v gmp | grep -v mpfr > mplapack_eig_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_eig_${mplib}.h 
        sed -i -e 's/COMPLEX/std::complex<mplapack_binary80_t>/g' mplapack_eig_${mplib}.h
        sed -i -e 's/REAL/mplapack_binary80_t/g' mplapack_eig_${mplib}.h
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/Mxerbla/Mxerbla_${mplib}/g" mplapack_eig_${mplib}.h
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_eig_${mplib}.h
        sed -i -e "s/\<iMlaenv2stage\>/iMlaenv2stage_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/\<iMlaenv\>/iMlaenv_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_eig_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_eig_${mplib}.h 
    fi

    fable_clang_format_stdout mplapack_eig_${mplib}.h | LC_ALL=C sort > l ; mv l mplapack_eig_${mplib}.h 
    cat ~/mplapack/mplapack/test/eig/common/mplapack_eig_${mplib}.h.in mplapack_eig_${mplib}.h > ~/mplapack/include/mplapack_eig_${mplib}.h
    rm mplapack_eig_${mplib}.h
    echo "#endif" >> ~/mplapack/include/mplapack_eig_${mplib}.h


done

mv header_all mplapack_eig_generic.h

for f in mplapack_eig_generic.h; do
fable_clang_format_inplace "$f"
done
