FILES=`ls *cpp *h`

for _file in $FILES; do
clang-format -style="{BasedOnStyle: llvm, IndentWidth: 4, ColumnLimit: 10000}" $_file > ${_file}__ ; mv ${_file}__ ${_file}
done



