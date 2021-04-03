if [ $# -ne 1 ]; then
  echo "need a fortran file to convert"
  exit 1
fi

cp /home/docker/mplapack/misc/cout.py /home/docker/modules/cctbx_project/fable/cout.py
__output=${1%.*}.cpp
output=`basename $__output`
_output=${output}_

fable.cout $1 #| sed 's/fem::double0/0\.0/g' | sed 's/fem::int0/0/g' | sed 's/fem::fint/int/g' | sed 's/int/INTEGER/g' | sed -e "/dimension/d" | sed 's/fem:://g' | sed 's/double/REAL/g' | sed 's|//C|//|g' | sed -e 's/str_cref/const char */g'  -e 's/arr_cref<REAL>/REAL \*/g'  -e 's/arr_cref<REAL, 2>/REAL \*/g' -e 's/arr_cref<INTEGER>/INTEGER \*/g'  > $_output
clang-format -style="{BasedOnStyle: llvm, IndentWidth: 4, ColumnLimit: 10000 }" $_output > $output ; rm $_output
