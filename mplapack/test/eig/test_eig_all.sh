#!/bin/bash

EIGREALS=`ls xeigtstR_*`
EIGCOMPLEXES=`ls xeigtstC_*`

TESTREALS=`ls R*.in [a-z]*.in  | grep -v double`
TESTCOMPLEXES=`ls C*.in [a-z]*.in | grep -v double`

./xeigtstR_double < Rbal_double.in > log.xeigtstR_double.Rbal.in
./xeigtstC_double < Cbal_double.in > log.xeigtstC_double.Cbal.in

for eigreal in $EIGREALS; do
    for testreal in $TESTREALS; do
        ./$eigreal < $testreal > log.$eigreal.$testreal
    done
done

for eigcomplex in $EIGCOMPLEXS; do
    for testcomplex in $TESTCOMPLEXS; do
        ./$eigcomplex < $testcomplex > log.$eigcomplex.$testcomplex
    done
done
