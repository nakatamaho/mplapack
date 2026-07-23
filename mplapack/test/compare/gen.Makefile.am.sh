#!/bin/bash

set -eu

if [ `uname` = "Linux" ]; then
    SED=sed
else
    SED=gsed
fi
rm -f _tmpfilelist
pushd common ; ls *test.cpp > ../_tmpfilelist_test ; popd

MPLIBS="gmp mpfr binary128 dd qd double binary80"

for _mplib in $MPLIBS; do
    rm -f $_mplib/Makefile.am
    rm -f ${_mplib}/_filelist_test
    if [ "$_mplib" = "mpfr" ]; then
        _mplapack_mpfr_archive=""
        _mplapack_mpfr_lib=""
    else
        _mplapack_mpfr_archive=',$(top_builddir)/mplapack/reference/.libs/libmplapack_mpfr.a'
        _mplapack_mpfr_lib=' -lmplapack_mpfr'
    fi
    $SED \
        -e "s|%%MPLIB%%|${_mplib}|g" \
        -e "s|%%MPLAPACK_MPFR_ARCHIVE%%|${_mplapack_mpfr_archive}|g" \
        -e "s|%%MPLAPACK_MPFR_LIB%%|${_mplapack_mpfr_lib}|g" \
        $_mplib/Makefile.am.part > ${_mplib}/Makefile.am.part.${_mplib}
    FILES=`cat _tmpfilelist_test`
    for _file in $FILES; do
        echo "${_file}_${_mplib} \\" | $SED 's/\.cpp//g' >> ${_mplib}/_filelist_test
    done

    $SED -i -e "$ s/\\\//g" ${_mplib}/_filelist_test
    $SED -e "/%%insert here1%%/e cat ${_mplib}/_filelist_test" ${_mplib}/Makefile.am.part.${_mplib} > ${_mplib}/Makefile.am.part.2nd
    $SED -i -e "s/%%insert here1%%//g" ${_mplib}/Makefile.am.part.2nd

    rm -f ${_mplib}/Makefile.am.part.3rd
    FILES=`cat _tmpfilelist_test | grep -v Mutils`
    for _file in $FILES; do
        _function=`echo $_file | $SED -e "s/\./ /g" | awk '{print $1}'`
        grep %%FUNCTION%% ${_mplib}/Makefile.am.part.2nd | $SED "s/%%FUNCTION%%/$_function/g" >> ${_mplib}/Makefile.am.part.3rd 
        echo "" >> ${_mplib}/Makefile.am.part.3rd
    done

    $SED -e "/%%insert here2%%/e cat ${_mplib}/Makefile.am.part.3rd" ${_mplib}/Makefile.am.part.2nd | grep -v %%FUNCTION%% | grep -v %%TEST%% > ${_mplib}/Makefile.am
    $SED -i -e "s/%%insert here2%%//g" ${_mplib}/Makefile.am
    if grep -Eq -- '(-lmpblas_|libmpblas_|%%(MPLIB|MPLAPACK_MPFR_(ARCHIVE|LIB))%%)' ${_mplib}/Makefile.am; then
        echo "error: stale library name or unresolved placeholder in ${_mplib}/Makefile.am" >&2
        exit 1
    fi
    rm ${_mplib}/Makefile.am.part.* ${_mplib}/_filelist_test
done

rm -f _tmpfilelist_test
