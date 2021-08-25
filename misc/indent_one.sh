#!/bin/bash

if [ $# -lt 1 ]; then
    echo "specify filename"
    exit
fi
_file=$1
clang-format -style="{BasedOnStyle: llvm, IndentWidth: 4, ColumnLimit: 10000, SortIncludes: false}" $_file > ${_file}__ ; mv ${_file}__ ${_file}




