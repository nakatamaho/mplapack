
if [ $# -ne 1 ]; then
  echo "need a fortran file to convert"
  exit 1
fi

cp /home/docker/mplapack/misc/cout_mplapack.py /home/docker/modules/cctbx_project/fable/
cp /home/docker/mplapack/misc/cout_mplapack.commandline.py /home/docker/modules/cctbx_project/fable/command_line/cout_mplapack.py
cp /home/docker/mplapack/misc/fable_mplapack.cout /home/docker/build36/bin/

__output=${1%.*}.cpp
output=`basename $__output`

fable_mplapack.cout $1 | sed -e 's/max(1/max((INTEGER)1/g' -e 's/max(0/max((INTEGER)0/g' -e 's/max(2/max((INTEGER)2/g' -e 's/max(7/max((INTEGER)7/g' -e 's/max({1/max({(INTEGER)1/g' -e 's/max({0/max({(INTEGER)0/g' -e 's/max({2/max({(INTEGER)2/g' -e 's/max({7/max({(INTEGER)2/g'  -e 's/fem::double0/0\.0/g' | sed 's/fem::bool0/false/g' | sed 's/fem::int0/0/g' | sed 's/fem::fint/int/g' | sed 's/int /INTEGER /g' | sed -e "/dimension/d" | sed 's/fem:://g' | sed 's/double/REAL/g' | sed 's|//C|//|g' | sed -e 's/str_cref/const char */g' -e 's/arr_cref<REAL>/REAL \*/g' -e 's/arr_ref<REAL>/REAL \*/g'| sed -e 's/arr_cref<INTEGER>/INTEGER \*/g' -e 's/arr_ref<INTEGER>/INTEGER \*/g' -e 's/arr_cref<bool>/bool \*/g' -e 's/arr_ref<bool>/bool \*/g' | sed -e 's/arr_ref<REAL, 2> /REAL \*/g' -e 's/0\.0e0/0\.0/g' -e 's/1\.0e0/1\.0/g' -e 's/0\.0e+0/0\.0/g' -e 's/1\.0e+0/1\.0/g' -e 's/dabs/abs/g' -e 's/arr_cref<REAL, 2> /REAL \*/g' -e 's/std::complex<REAL>/COMPLEX/g' -e 's/_MPLAPACK_REPLACE_dble//g' -e 's/_MPLAPACK_REPLACE_dimag//g' -e 's/_MPLAPACK_REPLACE_dsign/sign/g' -e 's/_MPLAPACK_REPLACE_dconjg/conj/g' -e 's/_MPLAPACK_REPLACE_dsqrt/sqrt/g' -e 's/1\.e0/1\.0/g' -e 's/2\.e0/2\.0/g' -e 's/4096\.e0/4096\.0/g' -e 's/16777216\.e0/16777216/g' -e 's/0\.e0/0\.0/g' -e 's/\.0e0/.0/g' -e 's/arr_cref<COMPLEX, 2> /COMPLEX \*/g' -e 's/COMPLEX0/ \(0\.0, 0\.0\)/g'  | sed -e 's/arr_ref<COMPLEX, 2> /COMPLEX \*/g' -e 's/arr_cref<COMPLEX> /COMPLEX \*/g' > ${output}_
clang-format -style="{BasedOnStyle: llvm, IndentWidth: 4, ColumnLimit: 10000 }" ${output}_ > ${output}
cat $output | sed -e 's/_MPLAPACK_REPLACE_dcmplx/COMPLEX/g' -e 's/str_ref /char \*/g'   -e 's/arr_ref<COMPLEX> /COMPLEX \*/g' -e 's/arr_cref<COMPLEX> /COMPLEX \*/g' -e 's/arr_ref<COMPLEX, 2> /COMPLEX \*/g' -e 's/cabs/abs/g' -e 's/dabs1/RCabs1/g' -e 's/int /INTEGER /g' -e 's/arr_ref<int>/INTEGER \*/g' -e 's/arr_cref<int>/INTEGER \*/g'   -e 's/arr_cref<int, 2>/INTEGER \*/g' -e 's/str<1>/char/g'  -e 's/conjg/conj/g' -e 's/mINTEGER/mint/g' -e 's/poINTEGER/point/g' -e 's/arr_1d<3, int> isave(fill0)/INTEGER isave[3]/g' -e 's/arr_1d<1, int> idum(fill0)/INTEGER idum[1]/g' -e 's/arr_1d<2, REAL> colssq(fill0)/REAL colssq[2]/g' -e 's/arr_1d<2, REAL> ssq(fill0)/REAL ssq[2]/g' -e 's/arr_1d<1, REAL> dum(fill0)/REAL dum[1]/g' -e 's/arr_1d<1, COMPLEX> dum(fill0)/COMPLEX dum[1]/g' -e 's/arr_1d<1, COMPLEX> zdum(fill0)/COMPLEX zdum[1]/g' -e 's/= char0//g' -e 's/arr_1d<1, COMPLEX> cdummy(fill0)/COMPLEX cdummy[1]/g' -e 's/arr_1d<1, REAL> rdummy(fill0)/REAL rdummy[1]/g' -e 's/arr_1d<1, COMPLEX> cdum(fill0)/COMPLEX cdum[1]/g'  -e 's/arr_1d<1, COMPLEX> dummy(fill0)/COMPLEX dummy[1]/g' -e 's/arr_1d<1, COMPLEX> dummy1(fill0)/COMPLEX dummy1[1]/g' > ${output}_

LINE=`cat -n ${output}_ | grep 'side + trans' | head -1 | awk '{print $1}'` 
if [ x$LINE != x"" ]; then
   sed -i -e "$((LINE))i char side_trans[3];" -e "$((LINE))i side_trans[0]=side[0];" -e "$((LINE))i side_trans[1]=trans[0];" -e "$((LINE))i side_trans[2]='\\\0';" ${output}_
fi
sed -i -e 's/side + trans/side_trans/g' ${output}_
cat ${output}_ | sed -e '/#include <fem.hpp> \/\/ Fortran EMulation library of fable module/,/\/\/  =====================================================================/d' > ${output}
cat ${output} | sed -e 's/lsame/Mlsame/g' -e 's/xerbla/Mxerbla/g' -e 's/ilaenv/iMlaenv/g' > ${output}_
cat /home/docker/mplapack/misc/header_lapack ${output}_ | sed '/namespace placeholder_please_replace/d' > ${output}
clang-format -style="{BasedOnStyle: llvm, IndentWidth: 4, ColumnLimit: 10000 }" ${output} > ${output}_
mv ${output}_ ${output}


