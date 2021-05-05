
if [ $# -ne 1 ]; then
  echo "need a fortran file to convert"
  exit 1
fi

cp /home/docker/mplapack/misc/cout_mplapack.py /home/docker/modules/cctbx_project/fable/
cp /home/docker/mplapack/misc/intrinsics.py  /home/docker/modules/cctbx_project/fable/
cp /home/docker/mplapack/misc/cout_mplapack.commandline.py /home/docker/modules/cctbx_project/fable/command_line/cout_mplapack.py
cp /home/docker/mplapack/misc/fable_mplapack.cout /home/docker/build36/bin/

__output=${1%.*}.cpp
output=`basename $__output`
cp $1 $output
sed -i \
-e 's/A( 1, 1 )/mplp_a11/g' \
-e 's/A( 1, 2 )/mplp_a12/g' \
-e 's/A( 1, 3 )/mplp_a13/g' \
-e 's/A( 1, 4 )/mplp_a14/g' \
-e 's/A( 1, 5 )/mplp_a15/g' \
-e 's/A( 1, 6 )/mplp_a16/g' \
-e 's/A( 1, 7 )/mplp_a17/g' \
-e 's/A( 1, 8 )/mplp_a18/g' \
-e 's/A( 1, 9 )/mplp_a19/g' \
-e 's/A( 1, 10 )/mplp_a1_10/g' \
-e 's/A( 1, 11 )/mplp_a1_11/g' \
-e 's/A( 1, 12 )/mplp_a1_12/g' \
-e 's/A( 1, 13 )/mplp_a1_13/g' \
-e 's/A( 1, 14 )/mplp_a1_14/g' \
-e 's/B( 1, 1 )/mplp_b11/g' \
-e 's/B( 1, 2 )/mplp_b12/g' \
-e 's/B( 1, 3 )/mplp_b13/g' \
-e 's/B( 1, 4 )/mplp_b14/g' \
-e 's/B( 1, 5 )/mplp_b15/g' \
-e 's/B( 1, 6 )/mplp_b16/g' \
-e 's/B( 1, 7 )/mplp_b17/g' \
-e 's/B( 1, 8 )/mplp_b18/g' \
-e 's/B( 1, 9 )/mplp_b19/g' \
-e 's/B( 1, 10 )/mplp_b1_10/g' \
-e 's/B( 1, 11 )/mplp_b1_11/g' \
-e 's/B( 1, 12 )/mplp_b1_12/g' \
-e 's/B( 1, 13 )/mplp_b1_13/g' \
-e 's/B( 1, 14 )/mplp_b1_14/g' \
-e 's/C( 1, 1 )/mplp_c11/g' \
-e 's/C( 1, 2 )/mplp_c12/g' \
-e 's/C( 1, 3 )/mplp_c13/g' \
-e 's/C( 1, 4 )/mplp_c14/g' \
-e 's/C( 1, 5 )/mplp_c15/g' \
-e 's/C( 1, 6 )/mplp_c16/g' \
-e 's/C( 1, 7 )/mplp_c17/g' \
-e 's/C( 1, 8 )/mplp_c18/g' \
-e 's/C( 1, 9 )/mplp_c19/g' \
-e 's/C( 1, 10 )/mplp_c1_10/g' \
-e 's/C( 1, 11 )/mplp_c1_11/g' \
-e 's/C( 1, 12 )/mplp_c1_12/g' \
-e 's/C( 1, 13 )/mplp_c1_13/g' \
-e 's/C( 1, 14 )/mplp_c1_14/g' $output
unifdef -U_OPENMP -m $output
fable_mplapack.cout $output | sed -e 's/max(1/max((INTEGER)1/g' -e 's/max(0/max((INTEGER)0/g' -e 's/max(2/max((INTEGER)2/g' -e 's/max(7/max((INTEGER)7/g' -e 's/max({1/max({(INTEGER)1/g' -e 's/max({0/max({(INTEGER)0/g' -e 's/max({2/max({(INTEGER)2/g' -e 's/max({7/max({(INTEGER)2/g'  -e 's/fem::double0/0\.0/g' | sed 's/fem::bool0/false/g' | sed 's/fem::int0/0/g' | sed 's/fem::fint/int/g' | sed 's/int /INTEGER /g' | sed -e "/dimension/d" | sed 's/fem:://g' | sed 's/double/REAL/g' | sed 's|//C|//|g' | sed -e 's/str_cref/const char */g' -e 's/arr_cref<REAL>/REAL \*/g' -e 's/arr_ref<REAL>/REAL \*/g'| sed -e 's/arr_cref<INTEGER>/INTEGER \*/g' -e 's/arr_ref<INTEGER>/INTEGER \*/g' -e 's/arr_cref<bool>/bool \*/g' -e 's/arr_ref<bool>/bool \*/g' | sed -e 's/arr_ref<REAL, 2> /REAL \*/g' -e 's/0\.0e0/0\.0/g' -e 's/1\.0e0/1\.0/g' -e 's/0\.0e+0/0\.0/g' -e 's/1\.0e+0/1\.0/g' -e 's/dabs/abs/g' -e 's/arr_cref<REAL, 2> /REAL \*/g' -e 's/std::complex<REAL>/COMPLEX/g' -e 's/_MPLAPACK_REPLACE_dble//g' -e 's/_MPLAPACK_REPLACE_dimag//g' -e 's/_MPLAPACK_REPLACE_dsign/sign/g' -e 's/_MPLAPACK_REPLACE_dconjg/conj/g' -e 's/_MPLAPACK_REPLACE_dsqrt/sqrt/g' -e 's/1\.e0/1\.0/g' -e 's/2\.e0/2\.0/g' -e 's/4096\.e0/4096\.0/g' -e 's/16777216\.e0/16777216/g' -e 's/0\.e0/0\.0/g' -e 's/\.0e0/.0/g' -e 's/arr_cref<COMPLEX, 2> /COMPLEX \*/g' -e 's/COMPLEX0/ \(0\.0, 0\.0\)/g'  | sed -e 's/arr_ref<COMPLEX, 2> /COMPLEX \*/g' -e 's/arr_cref<COMPLEX> /COMPLEX \*/g' > ${output}_
clang-format -style="{BasedOnStyle: llvm, IndentWidth: 4, ColumnLimit: 10000 }" ${output}_ > ${output}
cat $output | sed -e 's/_MPLAPACK_REPLACE_dcmplx/COMPLEX/g' -e 's/str_ref /char \*/g'   -e 's/arr_ref<COMPLEX> /COMPLEX \*/g' -e 's/arr_cref<COMPLEX> /COMPLEX \*/g' -e 's/arr_ref<COMPLEX, 2> /COMPLEX \*/g' -e 's/cabs/abs/g' -e 's/dabs1/RCabs1/g' -e 's/int /INTEGER /g' -e 's/arr_ref<int>/INTEGER \*/g' -e 's/arr_cref<int>/INTEGER \*/g'   -e 's/arr_cref<int, 2>/INTEGER \*/g' -e 's/str<1>/char/g'  -e 's/conjg/conj/g' -e 's/mINTEGER/mint/g' -e 's/poINTEGER/point/g' -e 's/PrINTEGER/Print/g'  -e 's/arr_1d<3, int> isave(fill0)/INTEGER isave[3]/g' -e 's/arr_1d<1, int> idum(fill0)/INTEGER idum[1]/g' -e 's/arr_1d<2, REAL> colssq(fill0)/REAL colssq[2]/g' -e 's/arr_1d<2, REAL> ssq(fill0)/REAL ssq[2]/g' -e 's/arr_1d<1, REAL> dum(fill0)/REAL dum[1]/g' -e 's/arr_1d<1, COMPLEX> dum(fill0)/COMPLEX dum[1]/g' -e 's/arr_1d<1, COMPLEX> zdum(fill0)/COMPLEX zdum[1]/g' -e 's/= char0//g' -e 's/arr_1d<1, COMPLEX> cdummy(fill0)/COMPLEX cdummy[1]/g' -e 's/arr_1d<1, REAL> rdummy(fill0)/REAL rdummy[1]/g' -e 's/arr_1d<1, COMPLEX> cdum(fill0)/COMPLEX cdum[1]/g'  -e 's/arr_1d<1, COMPLEX> dummy(fill0)/COMPLEX dummy[1]/g' -e 's/arr_1d<1, COMPLEX> dummy1(fill0)/COMPLEX dummy1[1]/g' -e 's/common &cmn,//g'  > ${output}_

LINE=`cat -n ${output}_ | grep 'side + trans' | head -1 | awk '{print $1}'` 
if [ x$LINE != x"" ]; then
   sed -i -e "$((LINE))i char side_trans[3];" -e "$((LINE))i side_trans[0]=side[0];" -e "$((LINE))i side_trans[1]=trans[0];" -e "$((LINE))i side_trans[2]='\\\0';" ${output}_
fi
sed -i -e 's/side + trans/side_trans/g' ${output}_

sed -i -e 's/mplp_a11/\&a[(1-1)+(1-1)*lda]/g' \
       -e 's/mplp_a12/\&a[(1-1)+(2-1)*lda]/g' \
       -e 's/mplp_a13/\&a[(1-1)+(3-1)*lda]/g' \
       -e 's/mplp_a14/\&a[(1-1)+(4-1)*lda]/g' \
       -e 's/mplp_a15/\&a[(1-1)+(5-1)*lda]/g' \
       -e 's/mplp_a16/\&a[(1-1)+(6-1)*lda]/g' \
       -e 's/mplp_a17/\&a[(1-1)+(7-1)*lda]/g' \
       -e 's/mplp_a18/\&a[(1-1)+(8-1)*lda]/g' \
       -e 's/mplp_a19/\&a[(1-1)+(9-1)*lda]/g' \
       -e 's/mplp_a1_10/\&a[(1-1)+(10-1)*lda]/g' \
       -e 's/mplp_a1_11/\&a[(1-1)+(11-1)*lda]/g' \
       -e 's/mplp_a1_12/\&a[(1-1)+(12-1)*lda]/g' \
       -e 's/mplp_a1_13/\&a[(1-1)+(13-1)*lda]/g' \
       -e 's/mplp_a1_14/\&a[(1-1)+(14-1)*lda]/g' \
       -e 's/mplp_b11/\&b[(1-1)+(1-1)*ldb]/g' \
       -e 's/mplp_b12/\&b[(1-1)+(2-1)*ldb]/g' \
       -e 's/mplp_b13/\&b[(1-1)+(3-1)*ldb]/g' \
       -e 's/mplp_b14/\&b[(1-1)+(4-1)*ldb]/g' \
       -e 's/mplp_b15/\&b[(1-1)+(5-1)*ldb]/g' \
       -e 's/mplp_b16/\&b[(1-1)+(6-1)*ldb]/g' \
       -e 's/mplp_b17/\&b[(1-1)+(7-1)*ldb]/g' \
       -e 's/mplp_b18/\&b[(1-1)+(8-1)*ldb]/g' \
       -e 's/mplp_b19/\&b[(1-1)+(9-1)*ldb]/g' \
       -e 's/mplp_b1_10/\&b[(1-1)+(10-1)*ldb]/g' \
       -e 's/mplp_b1_11/\&b[(1-1)+(11-1)*ldb]/g' \
       -e 's/mplp_b1_12/\&b[(1-1)+(12-1)*ldb]/g' \
       -e 's/mplp_b1_13/\&b[(1-1)+(13-1)*ldb]/g' \
       -e 's/mplp_b1_14/\&b[(1-1)+(14-1)*ldb]/g' \
       -e 's/mplp_c11/\&c[(1-1)+(1-1)*ldc]/g' \
       -e 's/mplp_c12/\&c[(1-1)+(2-1)*ldc]/g' \
       -e 's/mplp_c13/\&c[(1-1)+(3-1)*ldc]/g' \
       -e 's/mplp_c14/\&c[(1-1)+(4-1)*ldc]/g' \
       -e 's/mplp_c15/\&c[(1-1)+(5-1)*ldc]/g' \
       -e 's/mplp_c16/\&c[(1-1)+(6-1)*ldc]/g' \
       -e 's/mplp_c17/\&c[(1-1)+(7-1)*ldc]/g' \
       -e 's/mplp_c18/\&c[(1-1)+(8-1)*ldc]/g' \
       -e 's/mplp_c19/\&c[(1-1)+(9-1)*ldc]/g' \
       -e 's/mplp_c1_10/\&c[(1-1)+(10-1)*ldc]/g' \
       -e 's/mplp_c1_11/\&c[(1-1)+(11-1)*ldc]/g' \
       -e 's/mplp_c1_12/\&c[(1-1)+(12-1)*ldc]/g' \
       -e 's/mplp_c1_13/\&c[(1-1)+(13-1)*ldc]/g' \
       -e 's/mplp_c1_14/\&c[(1-1)+(14-1)*ldc]/g' \
       ${output}_
cat ${output}_ | sed -e '/#include <fem.hpp> \/\/ Fortran EMulation library of fable module/,/\/\/  =====================================================================/d' > ${output}
cat ${output} | sed -e 's/lsame/Mlsame/g' -e 's/xerbla/Mxerbla/g' -e 's/ilaenv/iMlaenv/g' | grep -v '(\[star\])' >  ${output}_
cat /home/docker/mplapack/misc/header_lapack ${output}_ | sed '/namespace placeholder_please_replace/d' > ${output}
clang-format -style="{BasedOnStyle: llvm, IndentWidth: 4, ColumnLimit: 10000 }" ${output} > ${output}_
mv ${output}_ ${output}


