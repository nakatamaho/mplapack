#!/bin/bash

#this script only works on GNU/Linux

MPLIBS="gmp mpfr _Float128 dd qd double _Float64x"
for mplib in $MPLIBS; do
        cd ~/mplapack/mpblas/reference/.libs
        nm --defined-only libmpblas_${mplib}.so | awk '$2=="T"' | awk '{print $3}' | c++filt | sed 's/(/ /g' | awk '{print $1}' | tr a-z A-Z | awk '{print "#define SYMBOL_" "GCC_" $1 }' > ~/mplapack/include/mplapack_benchmark_${mplib}.h_1
        nm --defined-only libmpblas_${mplib}.so | awk '$2=="T"' | awk '{print "\""$3"\""}'  > ~/mplapack/include/mplapack_benchmark_${mplib}.h_2
        paste ~/mplapack/include/mplapack_benchmark_${mplib}.h_1 ~/mplapack/include/mplapack_benchmark_${mplib}.h_2 > ~/mplapack/include/mplapack_benchmark_${mplib}.h
        rm ~/mplapack/include/mplapack_benchmark_${mplib}.h_1 ~/mplapack/include/mplapack_benchmark_${mplib}.h_2
        echo "//mpblas symbols" > ll
        cat ~/mplapack/misc/header_vanilla ll ~/mplapack/include/mplapack_benchmark_${mplib}.h > ll2 ; mv ll2 ~/mplapack/include/mplapack_benchmark_${mplib}.h
        cd ~/mplapack/mplapack/reference/.libs
        nm --defined-only libmplapack_${mplib}.so | awk '$2=="T"' | awk '{print $3}' | c++filt | sed 's/(/ /g' | awk '{print $1}' | tr a-z A-Z | awk '{print "#define SYMBOL_" "GCC_" $1 }' > ~/mplapack/include/mplapack_benchmark_${mplib}.h_1
        nm --defined-only libmplapack_${mplib}.so | awk '$2=="T"' | awk '{print "\""$3"\""}' > ~/mplapack/include/mplapack_benchmark_${mplib}.h_2
        paste ~/mplapack/include/mplapack_benchmark_${mplib}.h_1 ~/mplapack/include/mplapack_benchmark_${mplib}.h_2 >~/mplapack/include/mplapack_benchmark_${mplib}.h_3
        echo "" > ll
        echo "//mplapack symbols" >> ll
        cat ll ~/mplapack/include/mplapack_benchmark_${mplib}.h_3 > ~/mplapack/include/mplapack_benchmark_${mplib}.h_2 
        cat ~/mplapack/include/mplapack_benchmark_${mplib}.h ~/mplapack/include/mplapack_benchmark_${mplib}.h_2 > ~/mplapack/include/mplapack_benchmark_${mplib}.h_3
clang-format -style="{BasedOnStyle: llvm, IndentWidth: 4, ColumnLimit: 10000, SortIncludes: false}" ~/mplapack/include/mplapack_benchmark_${mplib}.h_3 > ~/mplapack/include/mplapack_benchmark_${mplib}.h
        rm -f ll ll2 ~/mplapack/include/mplapack_benchmark_${mplib}.h_*
done
