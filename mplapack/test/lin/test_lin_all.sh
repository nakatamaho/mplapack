#!/bin/bash

LINREALS=`ls xlintstR_*`
LINCOMPLEXS=`ls xlintstC_*`

LINREAL_RFPS=`ls xlintstrfR_*`
LINCOMPLEX_RFPS=`ls xlintstrfC_*`

for linreal in $LINREALS; do
    ./$linreal < ./Rtest.in >& log.$linreal
done

for lincomplex in $LINCOMPLEXS; do
    ./$lincomplex < ./Ctest.in >& log.$lincomplex
done

for linreal_rfp in $LINREAL_RFPS; do
    ./$linreal_rfp < ./Rtest_rfp.in >& log.$linreal_rfp
done

for lincomplex_rfp in $LINCOMPLEX_RFPS; do
    ./$lincomplex_rfp < ./Ctest_rfp.in >& log.$lincomplex_rfp
done
