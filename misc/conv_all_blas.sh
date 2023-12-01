cd ~/mplapack/mpblas/reference
LAPACKVERSION=3.12.0

# dsdot and sdsdot are mixed precsion version
FILES=`ls ~/mplapack/external/lapack/work/internal/lapack-$LAPACKVERSION/BLAS/SRC/d*.f*  \
          ~/mplapack/external/lapack/work/internal/lapack-$LAPACKVERSION/BLAS/SRC/z*.f*  \
          ~/mplapack/external/lapack/work/internal/lapack-$LAPACKVERSION/BLAS/SRC/id*.f* \ 
          ~/mplapack/external/lapack/work/internal/lapack-$LAPACKVERSION/BLAS/SRC/iz*.f* \
           | grep -v dsdot -v sdsdot`

# filename (subroutine) name conversion
# dgemm.f -> DGEMM -> Rgemm
# dnrm2.f90 -> DNRM2 -> Rnrm2
# etc.

rm BLAS_TO_MPBLAS_FUNCTIONLIST.txt
echo "sed \\"      > BLAS_TO_MPBLAS_FUNCTIONLIST.txt

for _file in $FILES; do
    oldfilename=`basename $_file | sed -e 's/\.f$//' | sed -e 's/\.f90$//' ` 
    newfilename=`basename $_file | sed -e 's/^zdscal/CRscal/g' -e 's/^zdrot/CRrot/g' -e 's/^dcabs/RCabs/g' -e 's/^dzasum/RCasum/g' -e 's/^dznrm2/RCnrm2/g' | sed -e 's/^d/R/' -e 's/^z/C/' -e 's/^id/iR/' -e 's/^iz/iC/' -e 's/\.f$//' -e 's/\.f90$//' `
     echo "-e 's/$oldfilename/$newfilename/g' \\"   >> BLAS_TO_MPBLAS_FUNCTIONLIST.txt
done
sed -i '$ s/\\$//' BLAS_TO_MPBLAS_FUNCTIONLIST.txt
exit

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
