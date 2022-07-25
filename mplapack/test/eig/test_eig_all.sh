#!/bin/bash

JOBS=4
EIGREALS=`ls xeigtstR_*`
EIGCOMPLEXES=`ls xeigtstC_*`

TESTREALS=`ls R*.in [a-z]*.in  | grep -v double | grep -v ^log`
TESTCOMPLEXES=`ls C*.in [a-z]*.in | grep -v double | grep -v ^log`

rm -f .parallel.test_eig_all.sh

echo "/usr/bin/time ./xeigtstR_double < Rbal_double.in >& log.xeigtstR_double.Rbal.in" >> .parallel.test_eig_all.sh
echo "/usr/bin/time ./xeigtstC_double < Cbal_double.in >& log.xeigtstC_double.Cbal.in" >> .parallel.test_eig_all.sh

for eigreal in $EIGREALS; do
    for testreal in $TESTREALS; do
        echo testing $eigreal $testreal
        echo "/usr/bin/time ./$eigreal < $testreal >& log.$eigreal.$testreal" >> .parallel.test_eig_all.sh
    done
done

for eigcomplex in $EIGCOMPLEXES; do
    for testcomplex in $TESTCOMPLEXES; do
        echo testing $eigcomplex $testcomplex
        echo "/usr/bin/time ./$eigcomplex < $testcomplex >& log.$eigcomplex.$testcomplex" >> .parallel.test_eig_all.sh
    done
done

cat .parallel.test_eig_all.sh | parallel --jobs $JOBS

rm .parallel.test_eig_all.sh
