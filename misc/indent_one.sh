#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 file1.cpp [file2.cpp ...]"
    exit 1
fi

for file in "$@"; do
    echo "Formatting: $file"
    clang-format-19 -i -style '{
        BasedOnStyle: llvm,
        IndentWidth: 4,
        ColumnLimit: 10000,
        SortIncludes: false,
        AlignEscapedNewlines: LeftWithLastLine,
        SpaceBeforeRangeBasedForLoopColon: false,
        PointerAlignment: Right,
        NamespaceIndentation: Inner,
        AlwaysBreakTemplateDeclarations: No,
        BreakBeforeConceptDeclarations: Never,
        ReflowComments: true,
        SpacesInLineCommentPrefix: { Minimum: 1, Maximum: 1 },
    }' "$file"
done
