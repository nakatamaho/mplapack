cd ~/mplapack/mpblas/reference
FILES=`ls ~/mplapack/external/lapack/work/internal/lapack-3.9.1/BLAS/SRC/d*.f ~/mplapack/external/lapack/work/internal/lapack-3.9.1/BLAS/SRC/z*.f ~/mplapack/external/lapack/work/internal/lapack-3.9.1/BLAS/SRC/id*.f ~/mplapack/external/lapack/work/internal/lapack-3.9.1/BLAS/SRC/iz*.f | grep -v dsdot`

rm -f BLAS_LIST BLAS_LIST_  BLAS_LIST__
echo "sed \\" > BLAS_LIST
echo "-e 's///g'" >> BLAS_LIST__

for _file in $FILES; do
oldfilename=`basename $_file | sed -e 's/\.f$//'` 
oldfilenameUP=`basename $_file | sed -e 's/\.f$//' | tr a-z A-Z`
newfilename=`basename $_file | sed -e 's/^zdscal/CRscal/g' -e 's/^zdrot/CRrot/g' -e 's/^dcabs/RCabs/g' -e 's/^dzasum/RCasum/g' -e 's/^dznrm2/RCnrm2/g' | sed -e 's/^d/R/' -e 's/^z/C/' -e 's/^id/iR/' -e 's/^iz/iC/' -e 's/\.f$//'`
echo "-e 's/$oldfilename/$newfilename/g' \\" >> BLAS_LIST_
echo "-e 's/$oldfilenameUP/$newfilename/g' \\" >> BLAS_LIST_
done
cat BLAS_LIST_ | sort -r > BLAS_LIST___
mv BLAS_LIST___  BLAS_LIST_
cat BLAS_LIST BLAS_LIST_  BLAS_LIST__ > ll
mv ll BLAS_LIST
rm BLAS_LIST_*

for _file in $FILES; do
bash ~/mplapack/misc/fable_convert_blas.sh $_file
oldfilename=`basename $_file | sed -e 's/\.f$//'`
newfilename=`basename $_file | sed -e 's/^zdscal/CRscal/g' -e 's/^zdrot/CRrot/g' -e 's/^dcabs/RCabs/g' -e 's/^dzasum/RCasum/g' -e 's/^dznrm2/RCnrm2/g' | sed -e 's/^d/R/' -e 's/^z/C/' -e 's/^id/iR/' -e 's/^iz/iC/' -e 's/\.f$//'`
cat ${oldfilename}.cpp | bash BLAS_LIST > ${newfilename}.cpp_
mv ${newfilename}.cpp_  ${newfilename}.cpp
sed -i -e 's/const &/const /g' ${newfilename}.cpp
grep -v star]\) ${newfilename}.cpp > ${newfilename}.cpp_ ; mv ${newfilename}.cpp_ ${newfilename}.cpp
rm ${oldfilename}.cpp
done
mv BLAS_LIST ~/mplapack/misc/
