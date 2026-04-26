#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 file1.cpp [file2.cpp ...]"
    exit 1
fi

. "$(dirname "$0")/clang_format_common.sh"

for file in "$@"; do
    echo "Formatting: $file"
    run_clang_format "$file"
done
