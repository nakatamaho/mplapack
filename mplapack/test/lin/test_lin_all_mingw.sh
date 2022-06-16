#!/bin/bash

export WINEPATH="/usr/x86_64-w64-mingw32/lib/;/usr/lib/gcc/x86_64-w64-mingw32/9.3-win32/;/usr/lib/gcc/x86_64-w64-mingw32/9.3-posix;/home/docker/MPLAPACK_MINGW/bin"
export WINEDEBUG="-all"

LINREALS=`ls *xlintstR_*exe | grep -v log`
LINCOMPLEXES=`ls *xlintstC_*exe | grep -v log`

LINREAL_RFPS=`ls *xlintstrfR_*exe | grep -v log`
LINCOMPLEX_RFPS=`ls *xlintstrfC_*exe | grep -v log`

for linreal in $LINREALS; do
    wine64 $linreal < ./Rtest.in >& log.$linreal
done

for lincomplex in $LINCOMPLEXES; do
    wine64 $lincomplex < ./Ctest.in >& log.$lincomplex
done

for linreal_rfp in $LINREAL_RFPS; do
    wine64 $linreal_rfp < ./Rtest_rfp.in >& log.$linreal_rfp
done

for lincomplex_rfp in $LINCOMPLEX_RFPS; do
    wine64 $lincomplex_rfp < ./Ctest_rfp.in >& log.$lincomplex_rfp
done
