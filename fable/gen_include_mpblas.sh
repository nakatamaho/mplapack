#!/bin/bash
fable_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${fable_script_dir}/clang_format_common.sh"

if [ `uname` = "Linux" ]; then
    SED=sed
else
    SED=gsed
fi

cd ~/mplapack/mpblas/reference

FILES=`ls *cpp | grep -v mplapackinit.cpp`

for _file in $FILES; do
ctags -x --c++-kinds=pf --language-force=c++ --_xformat='%{typeref} %{name} %{signature};' ${_file} |  tr ':' ' ' | sed -e 's/^typename //' > ${_file%.*}.hpp
done
ctags -x --c++-kinds=pf --language-force=c++ --_xformat='%{typeref} %{name} %{signature};' Mxerbla.cpp |  tr ':' ' ' | sed -e 's/^typename //' > Mxerbla.hpp
ctags -x --c++-kinds=pf --language-force=c++ --_xformat='%{typeref} %{name} %{signature};' Mlsame.cpp |  tr ':' ' ' | sed -e 's/^typename //' > Mlsame.hpp

cat *hpp | grep -v abssq > header_all
rm *hpp

MPLIBS="gmp mpfr binary128 dd qd double binary80"
for mplib in $MPLIBS; do
    if [ x"$mplib" = x"gmp" ]; then
        cp header_all mpblas_${mplib}.h
        sed -i -e 's/INTEGER/mplapackint/g' mpblas_${mplib}.h
        sed -i -e 's/COMPLEX/mpfc_class/g' mpblas_${mplib}.h
        sed -i -e 's/REAL/mpf_class/g' mpblas_${mplib}.h
        sed -i -e "s/Mlsame/Mlsame_${mplib}/g" mpblas_${mplib}.h
        sed -i -e "s/Mxerbla/Mxerbla_${mplib}/g" mpblas_${mplib}.h
    fi

    if [ x"$mplib" = x"mpfr" ]; then
        cp header_all mpblas_${mplib}.h
        sed -i -e 's/INTEGER/mplapackint/g' mpblas_${mplib}.h
        sed -i -e 's/COMPLEX/mpc_class/g' mpblas_${mplib}.h
        sed -i -e 's/REAL/mpfr_class/g' mpblas_${mplib}.h
        sed -i -e "s/Mlsame/Mlsame_${mplib}/g" mpblas_${mplib}.h
        sed -i -e "s/Mxerbla/Mxerbla_${mplib}/g" mpblas_${mplib}.h
    fi

    if [ x"$mplib" = x"double" ]; then
        cp header_all mpblas_${mplib}.h
        sed -i -e 's/INTEGER/mplapackint/g' mpblas_${mplib}.h
        sed -i -e 's/COMPLEX/std::complex<double>/g' mpblas_${mplib}.h
        sed -i -e 's/REAL/double/g' mpblas_${mplib}.h
        sed -i -e "s/Mlsame/Mlsame_${mplib}/g" mpblas_${mplib}.h
        sed -i -e "s/Mxerbla/Mxerbla_${mplib}/g" mpblas_${mplib}.h
    fi

    if [ x"$mplib" = x"dd" ]; then
        cp header_all mpblas_${mplib}.h
        sed -i -e 's/INTEGER/mplapackint/g' mpblas_${mplib}.h
        sed -i -e 's/COMPLEX/dd_complex/g' mpblas_${mplib}.h
        sed -i -e 's/REAL/dd_real/g' mpblas_${mplib}.h
        sed -i -e "s/Mlsame/Mlsame_${mplib}/g" mpblas_${mplib}.h
        sed -i -e "s/Mxerbla/Mxerbla_${mplib}/g" mpblas_${mplib}.h
    fi


    if [ x"$mplib" = x"qd" ]; then
        cp header_all mpblas_${mplib}.h
        sed -i -e 's/INTEGER/mplapackint/g' mpblas_${mplib}.h
        sed -i -e 's/COMPLEX/qd_complex/g' mpblas_${mplib}.h
        sed -i -e 's/REAL/qd_real/g' mpblas_${mplib}.h
        sed -i -e "s/Mlsame/Mlsame_${mplib}/g" mpblas_${mplib}.h
        sed -i -e "s/Mxerbla/Mxerbla_${mplib}/g" mpblas_${mplib}.h
    fi

    if [ x"$mplib" = x"binary128" ]; then
        cp header_all mpblas_${mplib}.h
        sed -i -e 's/INTEGER/mplapackint/g' mpblas_${mplib}.h
        sed -i -e 's/COMPLEX/std::complex<mplapack_binary128_t>/g' mpblas_${mplib}.h
        sed -i -e 's/REAL/mplapack_binary128_t/g' mpblas_${mplib}.h
        sed -i -e "s/Mlsame/Mlsame_${mplib}/g" mpblas_${mplib}.h
        sed -i -e "s/Mxerbla/Mxerbla_${mplib}/g" mpblas_${mplib}.h
    fi

    if [ x"$mplib" = x"binary80" ]; then
        cp header_all mpblas_${mplib}.h
        sed -i -e 's/INTEGER/mplapackint/g' mpblas_${mplib}.h
        sed -i -e 's/COMPLEX/std::complex<mplapack_binary80_t>/g' mpblas_${mplib}.h
        sed -i -e 's/REAL/mplapack_binary80_t/g' mpblas_${mplib}.h
        sed -i -e "s/Mlsame/Mlsame_${mplib}/g" mpblas_${mplib}.h
        sed -i -e "s/Mxerbla/Mxerbla_${mplib}/g" mpblas_${mplib}.h
    fi

    fable_clang_format_stdout mpblas_${mplib}.h | LC_ALL=C sort > l ; mv l mpblas_${mplib}.h
    cat ~/mplapack/mpblas/reference/mpblas_${mplib}.h.in mpblas_${mplib}.h > ~/mplapack/include/mpblas_${mplib}.h
    rm mpblas_${mplib}.h
    echo "#endif" >> ~/mplapack/include/mpblas_${mplib}.h

done
mv header_all mpblas_generic.h

for f in mpblas_generic.h; do
fable_clang_format_inplace "$f"
done
