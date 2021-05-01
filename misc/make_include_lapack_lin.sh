#!/bin/bash

cd ~/mplapack/mplapack/test/lin/common

if [ `uname` = "Linux" ]; then
    SED=sed
else
    SED=gsed
fi

FILES=`ls *cpp | grep -v Rlamch | grep -v Mlaenv | grep -v Mutils`
for filename in $FILES; do
/usr/local/bin/ctags -x --c++-kinds=pf --language-force=c++ --_xformat='%{typeref} %{name} %{signature};' ${filename} |  tr ':' ' ' | sed -e 's/^typename //' >  ${filename%.*}.hpp
done

cat hand/*hpp *hpp ~/mplapack/misc/special.hpp | grep -v abs1 | grep -v abs2 | grep -v abssq | grep -v arr_c | grep -v common | grep -v UNHANDLED_function | grep -v arr_ref | grep -v random_mplapack_gmp | grep -v mplapack_Rlaruv_gmp_initialize | grep -v mplapack_Rlaruv_mpfr_initialize | sort |uniq > header_all

rm *hpp

MPLIBS="gmp mpfr _Float128 dd qd double _Float64x"
for mplib in $MPLIBS; do
    if [ x"$mplib" = x"gmp" ]; then
        cat header_all | grep -v mpfr > mplapack_lin_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_lin_${mplib}.h 
        sed -i -e 's/COMPLEX/mpc_class/g' mplapack_lin_${mplib}.h 
        sed -i -e 's/REAL/mpf_class/g' mplapack_lin_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_lin_${mplib}.h
        sed -i -e "s/iMlaenv2stage/iMlaenv2stage_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMlaenv(/iMlaenv_${mplib}(/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_lin_${mplib}.h 
    fi

    if [ x"$mplib" = x"mpfr" ]; then
        cat header_all | grep -v gmp > mplapack_lin_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_lin_${mplib}.h 
        sed -i -e 's/COMPLEX/mpcomplex/g' mplapack_lin_${mplib}.h 
        sed -i -e 's/REAL/mpreal/g' mplapack_lin_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_lin_${mplib}.h
        sed -i -e "s/iMlaenv2stage/iMlaenv2stage_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMlaenv(/iMlaenv_${mplib}(/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_lin_${mplib}.h 
    fi

    if [ x"$mplib" = x"double" ]; then
        cat header_all | grep -v gmp | grep -v mpfr > mplapack_lin_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_lin_${mplib}.h 
        sed -i -e 's/COMPLEX/std::complex<double>/g' mplapack_lin_${mplib}.h 
        sed -i -e 's/REAL/double/g' mplapack_lin_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_lin_${mplib}.h
        sed -i -e "s/iMlaenv2stage/iMlaenv2stage_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMlaenv(/iMlaenv_${mplib}(/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_lin_${mplib}.h 
    fi

    if [ x"$mplib" = x"dd" ]; then
        cat header_all | grep -v gmp | grep -v mpfr > mplapack_lin_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_lin_${mplib}.h 
        sed -i -e 's/COMPLEX/dd_complex/g' mplapack_lin_${mplib}.h 
        sed -i -e 's/REAL/dd_real/g' mplapack_lin_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_lin_${mplib}.h
        sed -i -e "s/iMlaenv2stage/iMlaenv2stage_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMlaenv(/iMlaenv_${mplib}(/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_lin_${mplib}.h 
    fi

    if [ x"$mplib" = x"qd" ]; then
        cat header_all | grep -v gmp | grep -v mpfr > mplapack_lin_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_lin_${mplib}.h 
        sed -i -e 's/COMPLEX/qd_complex/g' mplapack_lin_${mplib}.h 
        sed -i -e 's/REAL/qd_real/g' mplapack_lin_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_lin_${mplib}.h
        sed -i -e "s/iMlaenv2stage/iMlaenv2stage_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMlaenv(/iMlaenv_${mplib}(/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_lin_${mplib}.h 
    fi

    if [ x"$mplib" = x"_Float128" ]; then
        cat header_all | grep -v gmp | grep -v mpfr > mplapack_lin_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_lin_${mplib}.h 
        sed -i -e 's/COMPLEX/std::complex<_Float128>/g' mplapack_lin_${mplib}.h 
        sed -i -e 's/REAL/_Float128/g' mplapack_lin_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_lin_${mplib}.h
        sed -i -e "s/iMlaenv2stage/iMlaenv2stage_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMlaenv(/iMlaenv_${mplib}(/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_lin_${mplib}.h 
    fi

    if [ x"$mplib" = x"_Float64x" ]; then
        cat header_all | grep -v gmp | grep -v mpfr > mplapack_lin_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_lin_${mplib}.h 
        sed -i -e 's/COMPLEX/std::complex<_Float64x>/g' mplapack_lin_${mplib}.h 
        sed -i -e 's/REAL/_Float64x/g' mplapack_lin_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_lin_${mplib}.h
        sed -i -e "s/iMlaenv2stage/iMlaenv2stage_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMlaenv(/iMlaenv_${mplib}(/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_lin_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_lin_${mplib}.h 
    fi

    clang-format -style="{BasedOnStyle: llvm, IndentWidth: 4, ColumnLimit: 10000 }" mplapack_lin_${mplib}.h | sort > l ; mv l mplapack_lin_${mplib}.h 
    cat ~/mplapack/misc/mplapack_lin_${mplib}.h.in mplapack_lin_${mplib}.h > ~/mplapack/include/mplapack_lin_${mplib}.h
    rm mplapack_lin_${mplib}.h
    echo "#endif" >> ~/mplapack/include/mplapack_lin_${mplib}.h

done
#rm header_all mplapack.h *hpp
#rm mplapack.h
