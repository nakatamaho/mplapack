cd ~/mplapack/mplapack/reference
FILES=`ls ~/mplapack/external/lapack/work/internal/lapack-3.9.1/SRC/d*.f ~/mplapack/external/lapack/work/internal/lapack-3.9.1/SRC/z*.f ~/mplapack/external/lapack/work/internal/lapack-3.9.1/SRC/id*.f ~/mplapack/external/lapack/work/internal/lapack-3.9.1/SRC/iz*.f | grep -v dsdot | grep dpot` #sed -n '50,60p'`

rm -f LAPACK_LIST LAPACK_LIST_  LAPACK_LIST__
echo "sed \\" > LAPACK_LIST
echo "-e 's///g'" >> LAPACK_LIST__

for _file in $FILES; do
oldfilename=`basename $_file | sed -e 's/\.f$//'` 
oldfilenameUP=`basename $_file | sed -e 's/\.f$//' | tr a-z A-Z`
newfilename=`basename $_file | sed -e 's/^zdscal/CRscal/g' -e 's/^zdrot/CRrot/g' -e 's/^dcabs/RCabs/g' -e 's/^dzasum/RCasum/g' -e 's/^dznrm2/RCnrm2/g' | sed -e 's/^d/R/' -e 's/^z/C/' -e 's/^id/iR/' -e 's/^iz/iC/' -e 's/\.f$//'`
echo "-e 's/$oldfilename/$newfilename/g' \\" >> LAPACK_LIST_
echo "-e 's/$oldfilenameUP/$newfilename/g' \\" >> LAPACK_LIST_
done
cat LAPACK_LIST_ | sort -r > LAPACK_LIST___
mv LAPACK_LIST___  LAPACK_LIST_
cat LAPACK_LIST LAPACK_LIST_  LAPACK_LIST__ > ll
mv ll LAPACK_LIST
rm LAPACK_LIST_*

cp ~/mplapack/misc/BLAS_LIST .

rm -f FILELIST
for _file in $FILES; do
bash ~/mplapack/misc/fem_convert_blas.sh $_file
oldfilename=`basename $_file | sed -e 's/\.f$//'`
newfilename=`basename $_file | sed -e 's/^zdscal/CRscal/g' -e 's/^zdrot/CRrot/g' -e 's/^dcabs/RCabs/g' -e 's/^dzasum/RCasum/g' -e 's/^dznrm2/RCnrm2/g' | sed -e 's/^d/R/' -e 's/^z/C/' -e 's/^id/iR/' -e 's/^iz/iC/' -e 's/\.f$//'`
cat ${oldfilename}.cpp | bash BLAS_LIST | bash LAPACK_LIST > ${newfilename}.cpp_
mv ${newfilename}.cpp_  ${newfilename}.cpp
sed -i -e 's/const &/const /g' ${newfilename}.cpp
/usr/local/bin/ctags -x --c++-kinds=pf --language-force=c++ --_xformat='%{typeref} %{name} %{signature};' ${newfilename}.cpp |  tr ':' ' ' | sed -e 's/^typename //' > ${newfilename}.hpp
rm ${oldfilename}.cpp
echo "${newfilename}.cpp \\" >> FILELIST
done

sed -e "/%%insert here%%/e cat FILELIST" Makefile.am.in > Makefile.am
sed -i -e "s/%%insert here%%//g" Makefile.am
head -c -1 Makefile.am > Makefile.am_
mv Makefile.am_ Makefile.am
sed -i -e '$s/\\//' Makefile.am

rm BLAS_LIST

mv LAPACK_LIST ~/mplapack/misc/
cat *hpp |sort > mplapack.h
cat ~/mplapack/misc/header mplapack.h > mplapack.h_
clang-format -style="{BasedOnStyle: llvm, IndentWidth: 4, ColumnLimit: 10000 }" mplapack.h_ > mplapack.h
rm mplapack.h_

