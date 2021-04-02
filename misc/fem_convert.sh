if [ $# -ne 1 ]; then
  echo "need a fortran file to convert"
  exit 1
fi

__output=${1%.*}.cpp
output=`basename $__output`
_output=${output}_
fable.cout $1 | sed 's/fem::double0/0\.0/g' | sed 's/fem::int0/0/g'  | sed 's/int/INTEGER/g' | sed -e "/dimension/d" | sed 's/fem:://g' | sed 's/double/REAL/g' | sed 's|//C|//|g' | sed 's/arr_cref<REAL>/REAL \*/g' > $_output
clang-format -style="{BasedOnStyle: llvm, IndentWidth: 4, ColumnLimit: 10000 }" $_output > $output ; rm $_output
