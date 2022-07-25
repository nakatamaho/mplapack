#!/bin/bash

JOBS=4

LINREALS=`ls xlintstR_*`
LINCOMPLEXES=`ls xlintstC_*`

LINREAL_RFPS=`ls xlintstrfR_*`
LINCOMPLEX_RFPS=`ls xlintstrfC_*`

rm -f .parallel.test_lin_all.sh

for linreal in $LINREALS; do
    echo "./$linreal < ./Rtest.in >& log.$linreal" >> .parallel.test_lin_all.sh
done

for lincomplex in $LINCOMPLEXES; do
    echo "/usr/bin/time ./$lincomplex < ./Ctest.in >& log.$lincomplex" >> .parallel.test_lin_all.sh
done

for linreal_rfp in $LINREAL_RFPS; do
    echo "/usr/bin/time ./$linreal_rfp < ./Rtest_rfp.in >& log.$linreal_rfp" >> .parallel.test_lin_all.sh
done

for lincomplex_rfp in $LINCOMPLEX_RFPS; do
    echo "/usr/bin/time ./$lincomplex_rfp < ./Ctest_rfp.in >& log.$lincomplex_rfp" >> .parallel.test_lin_all.sh
done

cat .parallel.test_lin_all.sh | parallel --jobs $JOBS

rm .parallel.test_lin_all.sh
