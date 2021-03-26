#!/bin/bash
#this script generates include/mplapack_benchmark_gmp.h, include/mplapack_benchmark__Float128.h
#include/mplapack_benchmark_mpfr.h, include/mplapack_benchmark_dd.h, include/mplapack_benchmark_qd.h
#include/mplapack_benchmark_double.h

for p in gmp mpfr _Float128 double dd qd
do
rm $p.h
for n in libmpblas_"$p"_ref*.o
do
m=`echo $n | sed s/libmplapack_"$p"_ref_la.// | sed s/\\\.o//`
mupper=`echo $m | tr '[a-z]' '[A-Z]' `
k=`nm $n | grep $m | awk '{print $3}' | grep ^_Z`
echo -e "#define SYMBOL_GCC_$mupper\t\t\"$k\"" >> $p.h
done

done
