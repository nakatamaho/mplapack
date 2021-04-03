if [ $# -ne 1 ]; then
  echo "need a fortran file to convert"
  exit 1
fi

cp /home/docker/mplapack/misc/cout.py /home/docker/modules/cctbx_project/fable/cout.py
__output=${1%.*}.cpp
output=`basename $__output | sed 's/^d/R/' `
_output=${output}_

fable.cout $1 | sed 's/fem::double0/0\.0/g' | sed 's/fem::bool0/false/g' | sed 's/fem::int0/0/g' | sed 's/fem::fint/int/g' | sed 's/int/INTEGER/g' | sed -e "/dimension/d" | sed 's/fem:://g' | sed 's/double/REAL/g' | sed 's|//C|//|g' | sed -e 's/str_cref/const char */g'  -e 's/arr_cref<REAL>/REAL \*/g'| sed -e 's/arr_cref<INTEGER>/INTEGER \*/g' | sed 's/arr_ref<REAL, 2> /REAL \*/g'  > $_output
cat $_output | sed -e '/#include <fem.hpp> \/\/ Fortran EMulation library of fable module/,/\/\/  =====================================================================/d' > ${_output}__
cat ${_output}__ | sed -e 's/lsame/Mlsame/g' -e 's/xerbra/Mxerbra/g' -e 's/disnan/Risnan/g'  -e 's/dgemv/Rgemv/g' -e 's/ddot/Rdot/g' -e 's/dscal/Rscal/g' -e 's/DPOTF2/Rpotf2/g' -e 's/dpotf2/Rpotf2/g'> ${_output}___
cat ${_output}___ | sed '/namespace placeholder_please_replace/d' > ${_output}__
clang-format -style="{BasedOnStyle: llvm, IndentWidth: 4, ColumnLimit: 10000 }" ${_output}__ > $output ; rm $_output ${_output}__  ${_output}___ 

