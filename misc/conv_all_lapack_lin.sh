#!/bin/bash

CACHED=no
if [ x"$CACHED" = x"no" ]; then
rm -f LAPACK_LIN_LIST LAPACK_LIN_LIST_  LAPACK_LIN_LIST__
echo "sed \\" > LAPACK_LIN_LIST
echo "-e 's///g'" >> LAPACK_LIN_LIST__

FILES_SUBSET=`ls ~/mplapack/external/lapack/work/internal/lapack-3.9.1/TESTING/LIN/d*.f ~/mplapack/external/lapack/work/internal/lapack-3.9.1/TESTING/LIN/z*.f ~/mplapack/external/lapack/work/internal/lapack-3.9.1/TESTING/LIN/a*.f`

for _file in $FILES_SUBSET; do
oldfilename=`basename $_file | sed -e 's/\.f$//'` 
oldfilenameUP=`basename $_file | sed -e 's/\.f$//' | tr a-z A-Z`
newfilename=`basename $_file | sed -e 's/^dzsum1/RCsum1/g' -e 's/^zdscal/CRscal/g' -e 's/^zdrot/CRrot/g' -e 's/^dcabs/RCabs/g' -e 's/^dzasum/RCasum/g' -e 's/^dznrm2/RCnrm2/g' | sed -e 's/^d/R/' -e 's/^a/A/' -e 's/^z/C/' -e 's/^id/iR/' -e 's/^iz/iC/' -e 's/^ila/iMla/' -e 's/\.f$//'`
echo "-e 's/$oldfilename/$newfilename/g' \\" >> LAPACK_LIN_LIST_
echo "-e 's/$oldfilenameUP/$newfilename/g' \\" >> LAPACK_LIN_LIST_
done
cat LAPACK_LIN_LIST_ | sort -r > LAPACK_LIN_LIST___
mv LAPACK_LIN_LIST___  LAPACK_LIN_LIST_
cat LAPACK_LIN_LIST LAPACK_LIN_LIST_  LAPACK_LIN_LIST__ > ll
mv ll LAPACK_LIN_LIST
rm LAPACK_LIN_LIST_*
else
cp ~/mplapack/misc/LAPACK_LIN_LIST .
fi

cp ~/mplapack/misc/LAPACK_LIST .
cp ~/mplapack/misc/BLAS_LIST .

for _file in $FILES_SUBSET ; do
perl -i.bk -pe 's/[^[:ascii:]]//g;' $_file
bash ~/mplapack/misc/fable_convert_lapack.sh $_file
oldfilename=`basename $_file | sed -e 's/\.f$//'`
newfilename=`basename $_file | sed -e 's/^dzsum1/RCsum1/g' -e 's/^zdscal/CRscal/g' -e 's/^zdrot/CRrot/g' -e 's/^dcabs/RCabs/g' -e 's/^dzasum/RCasum/g' -e 's/^dznrm2/RCnrm2/g' | sed -e 's/^d/R/' -e 's/^z/C/' -e 's/^a/A/' -e 's/^id/iR/' -e 's/^iz/iC/' -e 's/^ila/iMla/' -e 's/\.f$//'`

if [ ! -e $newfilename ]; then
cat ${oldfilename}.cpp | bash BLAS_LIST | bash LAPACK_LIN_LIST | bash LAPACK_LIST | sed 's/dlamch/Rlamch/g' > ${newfilename}.cpp_
mv ${newfilename}.cpp_  ${newfilename}.cpp
sed -i -e 's/const &/const /g' ${newfilename}.cpp
sed -i -e 's/, a\[/, \&a\[/g' ${newfilename}.cpp
sed -i -e 's/, b\[/, \&b\[/g' ${newfilename}.cpp
sed -i -e 's/, c\[/, \&c\[/g' ${newfilename}.cpp
sed -i -e 's/, d\[/, \&d\[/g' ${newfilename}.cpp
sed -i -e 's/, e\[/, \&e\[/g' ${newfilename}.cpp
sed -i -e 's/, h\[/, \&h\[/g' ${newfilename}.cpp
sed -i -e 's/, v\[/, \&v\[/g' ${newfilename}.cpp
sed -i -e 's/, u\[/, \&u\[/g' ${newfilename}.cpp
sed -i -e 's/, q\[/, \&q\[/g' ${newfilename}.cpp
sed -i -e 's/, t\[/, \&t\[/g' ${newfilename}.cpp
sed -i -e 's/, w\[/, \&w\[/g' ${newfilename}.cpp
sed -i -e 's/, x\[/, \&x\[/g' ${newfilename}.cpp
sed -i -e 's/, z\[/, \&z\[/g' ${newfilename}.cpp
sed -i -e 's/, ab\[/, \&ab\[/g' ${newfilename}.cpp
sed -i -e 's/, afac\[/, \&afac\[/g' ${newfilename}.cpp
sed -i -e 's/, ap\[/, \&ap\[/g' ${newfilename}.cpp
sed -i -e 's/, vt\[/, \&vt\[/g' ${newfilename}.cpp
sed -i -e 's/, tau\[/, \&tau\[/g' ${newfilename}.cpp
sed -i -e 's/, dum\[/, \&dum\[/g' ${newfilename}.cpp
sed -i -e 's/, cdum\[/, \&cdum\[/g' ${newfilename}.cpp
sed -i -e 's/, ipiv\[/, \&ipiv\[/g' ${newfilename}.cpp
sed -i -e 's/, work\[/, \&work\[/g' ${newfilename}.cpp
sed -i -e 's/, rwork\[/, \&rwork\[/g' ${newfilename}.cpp
sed -i -e 's/, lwork\[/, \&lwork\[/g' ${newfilename}.cpp
sed -i -e 's/, iwork\[/, \&iwork\[/g' ${newfilename}.cpp
fi
rm -f ${oldfilename}.cpp
done
