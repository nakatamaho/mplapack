FILES=`ls ~/mplapack/external/lapack/work/internal/lapack-3.9.1/SRC/d*.f ~/mplapack/external/lapack/work/internal/lapack-3.9.1/SRC/z*.f ~/mplapack/external/lapack/work/internal/lapack-3.9.1/SRC/id*.f ~/mplapack/external/lapack/work/internal/lapack-3.9.1/SRC/iz*.f | grep -v dsdot`

#FILES_SUBSET=`ls ~/mplapack/external/lapack/work/internal/lapack-3.9.1/SRC/d*.f ~/mplapack/external/lapack/work/internal/lapack-3.9.1/SRC/z*.f ~/mplapack/external/lapack/work/internal/lapack-3.9.1/SRC/id*.f ~/mplapack/external/lapack/work/internal/lapack-3.9.1/SRC/iz*.f | grep -v dsdot | grep -e dpot -e disn -e isna -e uum -e lauu2 -e trt`

FILES_SUBSET=`ls ~/mplapack/external/lapack/work/internal/lapack-3.9.1/SRC/d*.f ~/mplapack/external/lapack/work/internal/lapack-3.9.1/SRC/z*.f ~/mplapack/external/lapack/work/internal/lapack-3.9.1/SRC/i{d,l}*.f ~/mplapack/external/lapack/work/internal/lapack-3.9.1/SRC/iz*.f | \
grep -v dsg | grep -v zcg | grep -v svd | grep -v gesvj | grep -v x.f | grep -v la_syr | \
grep -e lae2 -e laev2 -e laev2 -e lassq -e lanst -e lansy -e lanhe -e lapy2 -e larfg -e lapy3 -e ladiv -e larfg -e lartg -e lart \
-e laset -e laset -e Rlasr -e lasr -e potf2  \
-e lacgv -e potf2 -e lascl -e lascl -e lasrt  \
-e steqr -e steqr -e sterf -e sytd2 -e hetd2  \
-e latrd -e latrd -e sytrd -e hetrd -e larf  \
-e larf -e  org2l -e ung2l -e org2r -e ung2r  \
-e larft -e larft -e larfb -e larfb -e orgqr  \
-e ungqr -e orgql -e ungql -e orgtr -e ungtr  \
-e syev -e  heev -e  potrf -e potrf -e lacrm  \
-e trti2 -e trti2 -e trtri -e trtri -e getf2  \
-e getf2 -e laswp -e laswp -e getrf -e getrf  \
-e getri -e getri -e getrs -e getrs -e gesv  \
-e gesv -e  trtrs -e trtrs -e lasyf -e lasyf \
-e lahef -e lacrt -e laesy -e rot -e   spmv   \
-e spr -e   symv -e  syr -e   max1 -e dsum1 \
-e nan -e dcombssq -e iladl`

CACHED=no
if [ x"$CACHED" = x"no" ]; then
rm -f LAPACK_LIST LAPACK_LIST_  LAPACK_LIST__
echo "sed \\" > LAPACK_LIST
echo "-e 's///g'" >> LAPACK_LIST__

FILES_SUBSET=`ls ~/mplapack/external/lapack/work/internal/lapack-3.9.1/SRC/d*.f ~/mplapack/external/lapack/work/internal/lapack-3.9.1/SRC/z*.f ~/mplapack/external/lapack/work/internal/lapack-3.9.1/SRC/i*.f | grep -v dsdot | grep -v zggsvp3 | grep il`

for _file in $FILES_SUBSET; do
oldfilename=`basename $_file | sed -e 's/\.f$//'` 
oldfilenameUP=`basename $_file | sed -e 's/\.f$//' | tr a-z A-Z`
newfilename=`basename $_file | sed -e 's/^dzsum1/RCsum1/g' -e 's/^zdscal/CRscal/g' -e 's/^zdrot/CRrot/g' -e 's/^dcabs/RCabs/g' -e 's/^dzasum/RCasum/g' -e 's/^dznrm2/RCnrm2/g' | sed -e 's/^d/R/' -e 's/^z/C/' -e 's/^id/iR/' -e 's/^iz/iC/' -e 's/^ila/iMla/' -e 's/\.f$//'`
echo "-e 's/$oldfilename/$newfilename/g' \\" >> LAPACK_LIST_
echo "-e 's/$oldfilenameUP/$newfilename/g' \\" >> LAPACK_LIST_
done
cat LAPACK_LIST_ | sort -r > LAPACK_LIST___
mv LAPACK_LIST___  LAPACK_LIST_
cat LAPACK_LIST LAPACK_LIST_  LAPACK_LIST__ > ll
mv ll LAPACK_LIST
rm LAPACK_LIST_*
else
cp ~/mplapack/misc/LAPACK_LIST .
fi

cp ~/mplapack/misc/BLAS_LIST .

for _file in $FILES_SUBSET ; do
bash ~/mplapack/misc/fem_convert_lapack.sh $_file
oldfilename=`basename $_file | sed -e 's/\.f$//'`
newfilename=`basename $_file | sed -e 's/^dzsum1/RCsum1/g' -e 's/^zdscal/CRscal/g' -e 's/^zdrot/CRrot/g' -e 's/^dcabs/RCabs/g' -e 's/^dzasum/RCasum/g' -e 's/^dznrm2/RCnrm2/g' | sed -e 's/^d/R/' -e 's/^z/C/' -e 's/^id/iR/' -e 's/^iz/iC/' -e 's/^ila/iMla/' -e 's/\.f$//'`

if [ ! -e $newfilename ]; then
cat ${oldfilename}.cpp | bash BLAS_LIST | bash LAPACK_LIST | sed 's/dlamch/Rlamch/g' > ${newfilename}.cpp_
mv ${newfilename}.cpp_  ${newfilename}.cpp
sed -i -e 's/const &/const /g' ${newfilename}.cpp
sed -i -e 's/, a\[/, \&a\[/g' ${newfilename}.cpp
sed -i -e 's/, b\[/, \&b\[/g' ${newfilename}.cpp
sed -i -e 's/, c\[/, \&c\[/g' ${newfilename}.cpp
sed -i -e 's/, d\[/, \&d\[/g' ${newfilename}.cpp
sed -i -e 's/, e\[/, \&e\[/g' ${newfilename}.cpp
sed -i -e 's/, h\[/, \&h\[/g' ${newfilename}.cpp
sed -i -e 's/, v\[/, \&v\[/g' ${newfilename}.cpp
sed -i -e 's/, t\[/, \&t\[/g' ${newfilename}.cpp
sed -i -e 's/, tau\[/, \&tau\[/g' ${newfilename}.cpp
sed -i -e 's/, w\[/, \&w\[/g' ${newfilename}.cpp
sed -i -e 's/, z\[/, \&z\[/g' ${newfilename}.cpp
sed -i -e 's/, ipiv\[/, \&ipiv\[/g' ${newfilename}.cpp
sed -i -e 's/, work\[/, \&work\[/g' ${newfilename}.cpp
fi
/usr/local/bin/ctags -x --c++-kinds=pf --language-force=c++ --_xformat='%{typeref} %{name} %{signature};' ${newfilename}.cpp |  tr ':' ' ' | sed -e 's/^typename //' > ${newfilename}.hpp
rm -f ${oldfilename}.cpp
done
