#!/bin/bash

if [ $# -le 2 ]; then
  echo "add_ldlibpath requies <prefix> and atleast one file"
  exit 2
fi

if [ `uname` = "Linux" ]; then
    SED=sed
elif [ `uname` = "Darwin" ]; then
    SED=gsed
else
    SED=sed
fi

argc=$#
argv=("$@")

for (( j=1; j<argc; j++ )); do
    $SED -i -e "s|%%PREFIX%%|${argv[0]}|g" ${argv[j]}
done
