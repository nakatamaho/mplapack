#!/bin/bash

export WINEPATH="/usr/x86_64-w64-mingw32/lib/;/usr/lib/gcc/x86_64-w64-mingw32/9.3-win32/;/usr/lib/gcc/x86_64-w64-mingw32/9.3-posix;/home/docker/MPLAPACK_MINGW/bin"
export WINEDEBUG="-all"

EIGREALS=`ls *xeigtstR_* | grep -v log`
EIGCOMPLEXES=`ls *xeigtstC_* | grep -v log`

TESTREALS=`ls R*.in [a-z]*.in  | grep -v double | grep -v ^log`
TESTCOMPLEXES=`ls C*.in [a-z]*.in | grep -v double | grep -v ^log`

/usr/bin/time wine64 ./*xeigtstR_double.exe < Rbal_double.in > log.xeigtstR_double.Rbal.in
/usr/bin/time wine64 ./*xeigtstC_double.exe < Cbal_double.in > log.xeigtstC_double.Cbal.in

for eigreal in $EIGREALS; do
    for testreal in $TESTREALS; do
        echo testing $eigreal $testreal
        /usr/bin/time wine64 ./$eigreal < $testreal > log.$eigreal.$testreal
    done
done

for eigcomplex in $EIGCOMPLEXES; do
    for testcomplex in $TESTCOMPLEXES; do
        echo testing $eigcomplex $testcomplex
        /usr/bin/time wine64 ./$eigcomplex < $testcomplex > log.$eigcomplex.$testcomplex
    done
done
