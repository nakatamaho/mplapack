cd ~/mplapack/mpblas/reference

#!/bin/bash

if [ `uname` = "Linux" ]; then
    SED=sed
else
    SED=gsed
fi

MPLIBS="gmp mpfr _Float128 dd qd double _Float64x"
cat *hpp > header_all

for mplib in $MPLIBS; do
    if [ x"$mplib" = x"gmp" ]; then
        cp header_all mpblas_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mpblas_${mplib}.h 
        sed -i -e 's/COMPLEX/mpc_class/g' mpblas_${mplib}.h 
        sed -i -e 's/REAL/mpc_class/g' mpblas_${mplib}.h 
        sed -i -e "s/Mlsame/Mlsame_${mplib}/g" mpblas_${mplib}.h 
        sed -i -e "s/Mxerbla/Mxerbla_${mplib}/g" mpblas_${mplib}.h 
    fi

    if [ x"$mplib" = x"mpfr" ]; then
        cp header_all mpblas_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mpblas_${mplib}.h 
        sed -i -e 's/COMPLEX/mpcomplex/g' mpblas_${mplib}.h 
        sed -i -e 's/REAL/mpreal/g' mpblas_${mplib}.h 
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

    if [ x"$mplib" = x"_Float128" ]; then
        cp header_all mpblas_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mpblas_${mplib}.h 
        sed -i -e 's/COMPLEX/std::complex<_Float128>/g' mpblas_${mplib}.h 
        sed -i -e 's/REAL/_Float128/g' mpblas_${mplib}.h 
        sed -i -e "s/Mlsame/Mlsame_${mplib}/g" mpblas_${mplib}.h 
        sed -i -e "s/Mxerbla/Mxerbla_${mplib}/g" mpblas_${mplib}.h 
    fi

    if [ x"$mplib" = x"_Float64x" ]; then
        cp header_all mpblas_${mplib}.h 
        sed -i -e 's/INTEGER/mplapackint/g' mpblas_${mplib}.h 
        sed -i -e 's/COMPLEX/std::complex<_Float64x>/g' mpblas_${mplib}.h 
        sed -i -e 's/REAL/_Float64x/g' mpblas_${mplib}.h 
        sed -i -e "s/Mlsame/Mlsame_${mplib}/g" mpblas_${mplib}.h 
        sed -i -e "s/Mxerbla/Mxerbla_${mplib}/g" mpblas_${mplib}.h 
    fi

    clang-format -style="{BasedOnStyle: llvm, IndentWidth: 4, ColumnLimit: 10000 }" mpblas_${mplib}.h > l ; mv l mpblas_${mplib}.h 
    cat ~/mplapack/misc/mpblas_${mplib}.h.in mpblas_${mplib}.h > ~/mplapack/include/mpblas_${mplib}.h
    echo "#endif" >> ~/mplapack/include/mpblas_${mplib}.h

done
