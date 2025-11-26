#!/bin/bash

cd ~/mplapack/mplapack/reference

if [ `uname` = "Linux" ]; then
    SED=sed
else
    SED=gsed
fi

FILES=`ls *cpp`
for filename in $FILES; do
/usr/local/bin/ctags -x --c++-kinds=pf --language-force=c++ --_xformat='%{typeref} %{name} %{signature};' ${filename} |  tr ':' ' ' | sed -e 's/^typename //' >  ${filename%.*}.hpp
done

cat *hpp | sort | grep -v abs1 | grep -v abs2 | grep -v Rlam | grep -v ___mplapack_ | grep -v abssq | grep -v ___random_mplapack_gmp > header_all

rm *hpp

#MPLIBS="gmp mpfr _Float128 dd qd double _Float64x"
for mplib in $MPLIBS; do
    if [ x"$mplib" = x"gmp" ]; then
        cat header_all | grep -v mpfr > mplapack_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_${mplib}.h 
        sed -i -e 's/COMPLEX/mpc_class/g' mplapack_${mplib}.h 
        sed -i -e 's/REAL/mpf_class/g' mplapack_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_${mplib}.h
        sed -i -e "s/iMlaenv2stage/iMlaenv2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMlaenv(/iMlaenv_${mplib}(/g" mplapack_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparam2stage/iMparam2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_${mplib}.h 
    fi

    if [ x"$mplib" = x"mpfr" ]; then
        cat header_all | grep -v gmp > mplapack_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_${mplib}.h 
        sed -i -e 's/COMPLEX/mpcomplex/g' mplapack_${mplib}.h 
        sed -i -e 's/REAL/mpreal/g' mplapack_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_${mplib}.h
        sed -i -e "s/iMlaenv2stage/iMlaenv2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMlaenv(/iMlaenv_${mplib}(/g" mplapack_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparam2stage/iMparam2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_${mplib}.h 
    fi

    if [ x"$mplib" = x"double" ]; then
        cat header_all | grep -v gmp | grep -v mpfr > mplapack_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_${mplib}.h 
        sed -i -e 's/COMPLEX/std::complex<double>/g' mplapack_${mplib}.h 
        sed -i -e 's/REAL/double/g' mplapack_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_${mplib}.h
        sed -i -e "s/iMlaenv2stage/iMlaenv2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMlaenv(/iMlaenv_${mplib}(/g" mplapack_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparam2stage/iMparam2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_${mplib}.h 
    fi

    if [ x"$mplib" = x"dd" ]; then
        cat header_all | grep -v gmp | grep -v mpfr > mplapack_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_${mplib}.h 
        sed -i -e 's/COMPLEX/dd_complex/g' mplapack_${mplib}.h 
        sed -i -e 's/REAL/dd_real/g' mplapack_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_${mplib}.h
        sed -i -e "s/iMlaenv2stage/iMlaenv2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMlaenv(/iMlaenv_${mplib}(/g" mplapack_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparam2stage/iMparam2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_${mplib}.h 
    fi

    if [ x"$mplib" = x"qd" ]; then
        cat header_all | grep -v gmp | grep -v mpfr > mplapack_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_${mplib}.h 
        sed -i -e 's/COMPLEX/qd_complex/g' mplapack_${mplib}.h 
        sed -i -e 's/REAL/qd_real/g' mplapack_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_${mplib}.h
        sed -i -e "s/iMlaenv2stage/iMlaenv2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMlaenv(/iMlaenv_${mplib}(/g" mplapack_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparam2stage/iMparam2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_${mplib}.h 
    fi

    if [ x"$mplib" = x"_Float128" ]; then
        cat header_all | grep -v gmp | grep -v mpfr > mplapack_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_${mplib}.h 
        sed -i -e 's/COMPLEX/std::complex<_Float128>/g' mplapack_${mplib}.h 
        sed -i -e 's/REAL/_Float128/g' mplapack_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_${mplib}.h
        sed -i -e "s/iMlaenv2stage/iMlaenv2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMlaenv(/iMlaenv_${mplib}(/g" mplapack_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparam2stage/iMparam2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_${mplib}.h 
    fi

    if [ x"$mplib" = x"_Float64x" ]; then
        cat header_all | grep -v gmp | grep -v mpfr > mplapack_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_${mplib}.h 
        sed -i -e 's/COMPLEX/std::complex<_Float64x>/g' mplapack_${mplib}.h 
        sed -i -e 's/REAL/_Float64x/g' mplapack_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_${mplib}.h
        sed -i -e "s/iMlaenv2stage/iMlaenv2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMlaenv(/iMlaenv_${mplib}(/g" mplapack_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparam2stage/iMparam2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_${mplib}.h 
    fi

    clang-format -style="{BasedOnStyle: llvm, IndentWidth: 4, ColumnLimit: 10000 }" mplapack_${mplib}.h | sort > l ; mv l mplapack_${mplib}.h 
    cat ~/mplapack/misc/mplapack_${mplib}.h.in mplapack_${mplib}.h > ~/mplapack/include/mplapack_${mplib}.h
    rm mplapack_${mplib}.h
    echo "#endif" >> ~/mplapack/include/mplapack_${mplib}.h

done

mv header_all mplapack_generic.h

for f in mplapack_generic.h; do
clang-format-19 -i -style '{
    BasedOnStyle: llvm,
    IndentWidth: 4,
    ColumnLimit: 10000,
    SortIncludes: false,
    AlignEscapedNewlines: LeftWithLastLine,
    SpaceBeforeRangeBasedForLoopColon: false,
    PointerAlignment: Right,
    NamespaceIndentation: Inner,
    AlwaysBreakTemplateDeclarations: No,
    BreakBeforeConceptDeclarations: Never,
  }' "$f"
done
