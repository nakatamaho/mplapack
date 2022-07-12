#!/bin/bash

LINREALS=`ls xlintstR_*`
LINCOMPLEXES=`ls xlintstC_*`

LINREAL_RFPS=`ls xlintstrfR_*`
LINCOMPLEX_RFPS=`ls xlintstrfC_*`

for linreal in $LINREALS; do
    echo testing $linreal Rtest.in
    ./$linreal < ./Rtest.in >& log.$linreal
done

for lincomplex in $LINCOMPLEXES; do
    echo testing $lincomplex Ctest.in
    /usr/bin/time ./$lincomplex < ./Ctest.in >& log.$lincomplex
done

for linreal_rfp in $LINREAL_RFPS; do
    echo testing $linreal_rfp Rtest_rfp.in
    /usr/bin/time ./$linreal_rfp < ./Rtest_rfp.in >& log.$linreal_rfp
done

for lincomplex_rfp in $LINCOMPLEX_RFPS; do
    echo testing $lincomplex_rfp Ctest_rfp.in
    /usr/bin/time ./$lincomplex_rfp < ./Ctest_rfp.in >& log.$lincomplex_rfp
done
