#!/bin/bash
#this script generates include/mpack_benchmark_gmp.h, include/mpack_benchmark___float128.h
#include/mpack_benchmark_mpfr.h, include/mpack_benchmark_dd.h, include/mpack_benchmark_qd.h
#include/mpack_benchmark_double.h

for p in gmp mpfr __float128 double dd qd
do
rm $p.h
for n in libmblas_"$p"_ref*.o
do
m=`echo $n | sed s/libmlapack_"$p"_ref_la.// | sed s/\\\.o//`
mupper=`echo $m | tr '[a-z]' '[A-Z]' `
k=`nm $n | grep $m | awk '{print $3}' | grep ^_Z`
echo -e "#define SYMBOL_GCC_$mupper\t\t\"$k\"" >> $p.h
done

done
