bash -x ~/mplapack/misc/conv_all_lapack.sh 
bash -x ~/mplapack/misc/make_include_lapack.sh
#FILES=`ls *cpp`
#for _file in $FILES; do
#    clang-format -style="{BasedOnStyle: llvm, IndentWidth: 4, ColumnLimit: 10000 }" $_file > l ; mv l $_file
#done
#patch -p0 < patch