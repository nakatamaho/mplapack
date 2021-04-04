cd ~/mplapack/mpblas/reference
FILES=`ls ~/mplapack/external/lapack/work/internal/lapack-3.9.1/BLAS/SRC/d*.f ~/mplapack/external/lapack/work/internal/lapack-3.9.1/BLAS/SRC/z*.f ~/mplapack/external/lapack/work/internal/lapack-3.9.1/BLAS/SRC/id*.f ~/mplapack/external/lapack/work/internal/lapack-3.9.1/BLAS/SRC/iz*.f`

rm -f BLAS_LIST
echo "sed \\" >> BLAS_LIST
for _file in $FILES; do
oldfilename=`basename $_file | sed -e 's/\.f$//'`
newfilename=`basename $_file | sed -e 's/^d/R/' -e 's/^z/C/' -e 's/^id/iR/' -e 's/^iz/iC/' -e 's/\.f$//'`
echo "-e 's/$oldfilename/$newfilename/g' \\" >> BLAS_LIST
done
echo "-e 's///g'" >> BLAS_LIST

i=0
for _file in $FILES; do
bash ~/mplapack/misc/fem_convert_blas.sh $_file
newfilename=`basename $_file | sed -e 's/^d/R/' -e 's/^z/C/' -e 's/^id/iR/' -e 's/^iz/iC/' -e 's/\.f$//'`
cat ${newfilename}.cpp | bash BLAS_LIST > ${newfilename}.cpp_ 
mv ${newfilename}.cpp_  ${newfilename}.cpp
if [ $i -ge 10 ]; then
exit
fi
i=$(( i+1 ))
done