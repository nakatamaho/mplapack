#!/bin/bash

JOBS=4

export WINEPATH="/usr/x86_64-w64-mingw32/lib/;/usr/lib/gcc/x86_64-w64-mingw32/9.3-win32/;/usr/lib/gcc/x86_64-w64-mingw32/9.3-posix;/home/docker/MPLAPACK_MINGW/bin"
export WINEDEBUG="-all"

LINREALS=`ls *xlintstR_*exe | grep -v log`
LINCOMPLEXES=`ls *xlintstC_*exe | grep -v log`

LINREAL_RFPS=`ls *xlintstrfR_*exe | grep -v log`
LINCOMPLEX_RFPS=`ls *xlintstrfC_*exe | grep -v log`

rm -f .parallel.test_lin_all_mingw.sh

for linreal in $LINREALS; do
    echo "wine64 $linreal < ./Rtest.in >& log.$linreal" >> .parallel.test_lin_all_mingw.sh
done

for lincomplex in $LINCOMPLEXES; do
    echo "wine64 $lincomplex < ./Ctest.in >& log.$lincomplex" >> .parallel.test_lin_all_mingw.sh
done

for linreal_rfp in $LINREAL_RFPS; do
    echo "wine64 $linreal_rfp < ./Rtest_rfp.in >& log.$linreal_rfp" >> .parallel.test_lin_all_mingw.sh
done

for lincomplex_rfp in $LINCOMPLEX_RFPS; do
    echo "wine64 $lincomplex_rfp < ./Ctest_rfp.in >& log.$lincomplex_rfp" >> .parallel.test_lin_all_mingw.sh
done

cat .parallel.test_lin_all_mingw.sh | parallel --jobs $JOBS

rm .parallel.test_lin_all_mingw.sh
