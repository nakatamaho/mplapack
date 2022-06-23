#!/bin/bash

EIGREALS=`ls xeigtstR_*`
EIGCOMPLEXES=`ls xeigtstC_*`

TESTREALS=`ls R*.in [a-z]*.in  | grep -v double | grep -v ^log`
TESTCOMPLEXES=`ls C*.in [a-z]*.in | grep -v double | grep -v ^log`

/usr/bin/time ./xeigtstR_double < Rbal_double.in > log.xeigtstR_double.Rbal.in
/usr/bin/time ./xeigtstC_double < Cbal_double.in > log.xeigtstC_double.Cbal.in

for eigreal in $EIGREALS; do
    for testreal in $TESTREALS; do
        echo testing $eigreal $testreal
        /usr/bin/time ./$eigreal < $testreal > log.$eigreal.$testreal
    done
done

for eigcomplex in $EIGCOMPLEXES; do
    for testcomplex in $TESTCOMPLEXES; do
        echo testing $eigcomplex $testcomplex
        /usr/bin/time ./$eigcomplex < $testcomplex > log.$eigcomplex.$testcomplex
    done
done
