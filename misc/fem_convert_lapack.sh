if [ $# -ne 1 ]; then
  echo "need a fortran file to convert"
  exit 1
fi

cp /home/docker/mplapack/misc/cout.py /home/docker/modules/cctbx_project/fable/cout.py
__output=${1%.*}.cpp
output=`basename $__output`

fable.cout $1 | sed -e 's/max(1/max((INTEGER)1/g'  -e 's/fem::double0/0\.0/g' | sed 's/fem::bool0/false/g' | sed 's/fem::int0/0/g' | sed 's/fem::fint/int/g' | sed 's/int/INTEGER/g' | sed -e "/dimension/d" | sed 's/fem:://g' | sed 's/double/REAL/g' | sed 's|//C|//|g' | sed -e 's/str_cref/const char */g' -e 's/arr_cref<REAL>/REAL \*/g' -e 's/arr_ref<REAL>/REAL \*/g'| sed -e 's/arr_cref<INTEGER>/INTEGER \*/g' -e 's/arr_ref<INTEGER>/INTEGER \*/g' -e 's/arr_ref<bool>/bool \*/g' | sed -e 's/arr_ref<REAL, 2> /REAL \*/g' -e 's/0\.0e0/0\.0/g' -e 's/1\.0e0/1\.0/g' -e 's/0\.0e+0/0\.0/g' -e 's/1\.0e+0/1\.0/g' -e 's/dabs/abs/g' -e 's/arr_cref<REAL, 2> /REAL \*/g' -e 's/std::complex<REAL>/COMPLEX/g' -e 's/_MPLAPACK_REPLACE_dble//g' -e 's/_MPLAPACK_REPLACE_dimag//g' -e 's/_MPLAPACK_REPLACE_dsign/sign/g' -e 's/_MPLAPACK_REPLACE_dconjg/conj/g' -e 's/_MPLAPACK_REPLACE_dsqrt/sqrt/g' -e 's/1\.e0/1\.0/g' -e 's/2\.e0/2\.0/g' -e 's/4096\.e0/4096\.0/g' -e 's/16777216\.e0/16777216/g' -e 's/0\.e0/0\.0/g' -e 's/\.0e0/.0/g' -e 's/arr_cref<COMPLEX, 2> /COMPLEX \*/g' -e 's/COMPLEX0/ \(0\.0, 0\.0\)/g'  > ${output}_
clang-format -style="{BasedOnStyle: llvm, IndentWidth: 4, ColumnLimit: 10000 }" ${output}_ > ${output}
cat $output | sed -e 's/_MPLAPACK_REPLACE_dcmplx/COMPLEX/g' -e 's/arr_ref<COMPLEX> /COMPLEX \*/g' -e 's/arr_cref<COMPLEX> /COMPLEX \*/g' -e 's/arr_ref<COMPLEX, 2> /COMPLEX \*/g' -e 's/cabs/abs/g' -e 's/dabs1/RCabs1/g' > ${output}_
cat ${output}_ | sed -e '/#include <fem.hpp> \/\/ Fortran EMulation library of fable module/,/\/\/  =====================================================================/d' > ${output}
cat ${output} | sed -e 's/lsame/Mlsame/g' -e 's/xerbla/Mxerbla/g' > ${output}_
cat /home/docker/mplapack/misc/header_blas ${output}_ | sed '/namespace placeholder_please_replace/d' > ${output}
clang-format -style="{BasedOnStyle: llvm, IndentWidth: 4, ColumnLimit: 10000 }" ${output} > ${output}_
mv ${output}_ ${output}
