#!/bin/bash

rm -f _tmpfilelist
pushd common ; ls *debug.cpp > ../_tmpfilelist_debug ; popd
pushd common ; ls *test.cpp  > ../_tmpfilelist_test ; popd

MPLIBS="gmp mpfr __float128 dd qd double longdouble"

for _mplib in $MPLIBS; do
    rm -f $_mplib/Makefile.am
    rm -f ${_mplib}/_filelist_debug ${_mplib}/_filelist_test
    cp $_mplib/Makefile.am.part ${_mplib}/Makefile.am.part.${_mplib}
    FILES=`cat _tmpfilelist_debug`
    for _file in $FILES; do
        echo "${_file}_${_mplib}_ref ${_file}_${_mplib} \\" | sed 's/\.cpp//g' >> ${_mplib}/_filelist_debug
    done
    FILES=`cat _tmpfilelist_test`
    for _file in $FILES; do
        echo "${_file}_${_mplib}_ref ${_file}_${_mplib} \\" | sed 's/\.cpp//g' >> ${_mplib}/_filelist_debug
    done

    sed -i -e "$ s/\\\//g" ${_mplib}/_filelist_debug
    sed -e "/%%insert here1%%/e cat ${_mplib}/_filelist_debug" ${_mplib}/Makefile.am.part.${_mplib} > ${_mplib}/Makefile.am.part.2nd
    sed -i -e "s/%%insert here1%%//g" ${_mplib}/Makefile.am.part.2nd

    rm -f ${_mplib}/Makefile.am.part.3rd
    FILES=`cat _tmpfilelist_debug | grep -v Mutils`
    for _file in $FILES; do
        _function=`echo $_file | sed -e "s/\./ /g" | awk '{print $1}'`
        grep %%FUNCTION%% ${_mplib}/Makefile.am.part.2nd | sed "s/%%FUNCTION%%/$_function/g" >> ${_mplib}/Makefile.am.part.3rd 
        echo "" >> ${_mplib}/Makefile.am.part.3rd
    done

    FILES=`cat _tmpfilelist_test | grep -v Mutils`
    for _file in $FILES; do
        _test=`echo $_file | grep -v Mutils | sed -e "s/\./ /g" | awk '{print $1}'`
        grep %%TEST%% ${_mplib}/Makefile.am.part.2nd | sed "s/%%TEST%%/$_test/g" >> ${_mplib}/Makefile.am.part.3rd
        echo "" >> ${_mplib}/Makefile.am.part.3rd
    done

    sed -e "/%%insert here2%%/e cat ${_mplib}/Makefile.am.part.3rd" ${_mplib}/Makefile.am.part.2nd | grep -v %%FUNCTION%% | grep -v %%TEST%% > ${_mplib}/Makefile.am
    sed -i -e "s/%%insert here2%%//g" ${_mplib}/Makefile.am
    rm ${_mplib}/Makefile.am.part.* ${_mplib}/_filelist_debug
done

rm -f _tmpfilelist_debug _tmpfilelist_test 
