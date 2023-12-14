if [ $# -ne 1 ]; then
  echo "need a fortran file to convert"
  exit 1
fi

cp /home/docker/mplapack/misc/cout_mplapack.py /home/docker/modules/cctbx_project/fable/
cp /home/docker/mplapack/misc/intrinsics.py  /home/docker/modules/cctbx_project/fable/
cp /home/docker/mplapack/misc/cout_mplapack.commandline.py /home/docker/modules/cctbx_project/fable/command_line/cout_mplapack.py
cp /home/docker/mplapack/misc/fable_mplapack.cout /home/docker/build38/bin/

__output=${1%.*}
output=`basename $__output`

# convert to C++
/home/docker/build38/bin/fable_mplapack.cout ${output}.f > ${output}.cpp
/home/docker/build38/bin/fable.cout ${output}.f > ${output}.cpp_vanilla
#rm ${output}.f ${output}.cpp_vanilla

# remove some stuffs
sed -i '/#include <fem.hpp>/,/using fem::common;/d' ${output}.cpp
sed -i 's|//C|//|g' ${output}.cpp
sed -i '/.*(\[ld.* \* star\])/d' ${output}.cpp
sed -i '/.*(\[star\]);/d' ${output}.cpp
sed -i '/namespace placeholder_please_replace/d' ${output}.cpp

clang-format -style="{BasedOnStyle: llvm, IndentWidth: 4, ColumnLimit: 10000 }" ${output}.cpp > ${output}.cpp_ ; mv ${output}.cpp_ ${output}.cpp
cat ${output}.cpp | sed -e 's/lsame/Mlsame/g' -e 's/xerbla/Mxerbla/g' > ${output}.cpp_        ; mv ${output}.cpp_ ${output}.cpp
cat /home/docker/mplapack/misc/header_blas ${output}.cpp              > ${output}.cpp_        ; mv ${output}.cpp_ ${output}.cpp
bash -c "$(cat /home/docker/mplapack/misc/BLAS_TO_MPBLAS_FUNCTIONLIST.sh; echo "${output}.cpp")"
