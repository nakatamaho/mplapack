#!/bin/bash

cd ~/mplapack/mplapack/reference

if [ `uname` = "Linux" ]; then
    SED=sed
else
    SED=gsed
fi

FILES=`ls *cpp | grep -v Rlamch`
for filename in $FILES; do
ctags -x --c++-kinds=pf --language-force=c++ --_xformat='%{typeref} %{name} %{signature};' ${filename} |  tr ':' ' ' | sed -e 's/^typename //' >  ${filename%.*}.hpp
done

printf "REAL Rlamch(const char *cmach);\nREAL Rlamc3(REAL a, REAL b);\n" > Rlamch.hpp

cat *hpp | sort | grep -v iseed_is_all_minus_one | grep -v rlaruv_print_nondet_banner_once | grep -v fixed_point | grep -v Mxerbla | grep -v abs1 | grep -v abs2 | grep -v ___mplapack_ | grep -v nondeterministic | grep -v advance_iseed | grep -v iseed_to_seed64 | grep -v abssq | grep -v ___random_mplapack_gmp > header_all

rm *hpp

MPLIBS="gmp mpfr binary128 dd qd double binary80"
for mplib in $MPLIBS; do
    if [ x"$mplib" = x"gmp" ]; then
        cat header_all | grep -v mpfr > mplapack_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_${mplib}.h 
        sed -i -e 's/COMPLEX/mpc_class/g' mplapack_${mplib}.h 
        sed -i -e 's/REAL/mpf_class/g' mplapack_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_${mplib}.h
        sed -i -e "s/\<iMlaenv2stage\>/iMlaenv2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/\<iMlaenv\>/iMlaenv_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparam2stage/iMparam2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/\<Rroundup_lwork\>/Rroundup_lwork_${mplib}/g" mplapack_${mplib}.h
        printf "void mplapack_gmp_initialize(void);" >> mplapack_${mplib}.h
    fi

    if [ x"$mplib" = x"mpfr" ]; then
        cat header_all | grep -v gmp > mplapack_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_${mplib}.h 
        sed -i -e 's/COMPLEX/mpcomplex/g' mplapack_${mplib}.h 
        sed -i -e 's/REAL/mpreal/g' mplapack_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_${mplib}.h
        sed -i -e "s/\<iMlaenv2stage\>/iMlaenv2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/\<iMlaenv\>/iMlaenv_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparam2stage/iMparam2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/\<Rroundup_lwork\>/Rroundup_lwork_${mplib}/g" mplapack_${mplib}.h
        printf "void ___mplapack_mpfr_initialize(void);" >> mplapack_${mplib}.h
        printf "void mplapack_mpfr_finalize(void);" >> mplapack_${mplib}.h
    fi

    if [ x"$mplib" = x"double" ]; then
        cat header_all | grep -v gmp | grep -v mpfr > mplapack_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_${mplib}.h 
        sed -i -e 's/COMPLEX/std::complex<double>/g' mplapack_${mplib}.h 
        sed -i -e 's/REAL/double/g' mplapack_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_${mplib}.h
        sed -i -e "s/\<iMlaenv2stage\>/iMlaenv2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/\<iMlaenv\>/iMlaenv_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparam2stage/iMparam2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/\<Rroundup_lwork\>/Rroundup_lwork_${mplib}/g" mplapack_${mplib}.h
    fi

    if [ x"$mplib" = x"dd" ]; then
        cat header_all | grep -v gmp | grep -v mpfr > mplapack_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_${mplib}.h 
        sed -i -e 's/COMPLEX/dd_complex/g' mplapack_${mplib}.h 
        sed -i -e 's/REAL/dd_real/g' mplapack_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_${mplib}.h
        sed -i -e "s/\<iMlaenv2stage\>/iMlaenv2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/\<iMlaenv\>/iMlaenv_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparam2stage/iMparam2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/\<Rroundup_lwork\>/Rroundup_lwork_${mplib}/g" mplapack_${mplib}.h
    fi

    if [ x"$mplib" = x"qd" ]; then
        cat header_all | grep -v gmp | grep -v mpfr > mplapack_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_${mplib}.h 
        sed -i -e 's/COMPLEX/qd_complex/g' mplapack_${mplib}.h 
        sed -i -e 's/REAL/qd_real/g' mplapack_${mplib}.h 
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_${mplib}.h
        sed -i -e "s/\<iMlaenv2stage\>/iMlaenv2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/\<iMlaenv\>/iMlaenv_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparam2stage/iMparam2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/\<Rroundup_lwork\>/Rroundup_lwork_${mplib}/g" mplapack_${mplib}.h
    fi

    if [ x"$mplib" = x"binary128" ]; then
        cat header_all | grep -v gmp | grep -v mpfr > mplapack_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_${mplib}.h 
        sed -i -e 's/COMPLEX/std::complex<mplapack_binary128_t>/g' mplapack_${mplib}.h
        sed -i -e 's/REAL/mplapack_binary128_t/g' mplapack_${mplib}.h
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_${mplib}.h
        sed -i -e "s/\<iMlaenv2stage\>/iMlaenv2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/\<iMlaenv\>/iMlaenv_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparam2stage/iMparam2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/\<Rroundup_lwork\>/Rroundup_lwork_${mplib}/g" mplapack_${mplib}.h
    fi

    if [ x"$mplib" = x"binary80" ]; then
        cat header_all | grep -v gmp | grep -v mpfr > mplapack_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mplapack_${mplib}.h
        sed -i -e 's/COMPLEX/std::complex<mplapack_binary80_t>/g' mplapack_${mplib}.h
        sed -i -e 's/REAL/mplapack_binary80_t/g' mplapack_${mplib}.h
        sed -i -e "s/Rlamch/Rlamch_${mplib}/g" mplapack_${mplib}.h
        sed -i -e "s/Mlsamen/Mlsamen_${mplib}/g" mplapack_${mplib}.h
        sed -i -e "s/\<iMlaenv2stage\>/iMlaenv2stage_${mplib}/g" mplapack_${mplib}.h
        sed -i -e "s/\<iMlaenv\>/iMlaenv_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMlaver/iMlaver_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMieeeck/iMieeeck_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparam2stage/iMparam2stage_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/iMparmq/iMparmq_${mplib}/g" mplapack_${mplib}.h 
        sed -i -e "s/\<Rroundup_lwork\>/Rroundup_lwork_${mplib}/g" mplapack_${mplib}.h
    fi

    clang-format-19 -style="{BasedOnStyle: llvm, IndentWidth: 4, ColumnLimit: 10000 }" mplapack_${mplib}.h | sort > l ; mv l mplapack_${mplib}.h 
    cat ~/mplapack/mplapack/reference/mplapack_${mplib}.h.in mplapack_${mplib}.h > ~/mplapack/include/mplapack_${mplib}.h
    rm mplapack_${mplib}.h
    echo "#endif" >> ~/mplapack/include/mplapack_${mplib}.h

done

mv header_all ~/mplapack/mplapack/reference/mplapack_generic.h

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
