for f in *.cpp *.cc *.h *.hpp *.h.in; do
  clang-format-19 -i -style '{
    BasedOnStyle: llvm,
    IndentWidth: 4,
    ColumnLimit: 10000,
    SortIncludes: false,
    AlignEscapedNewlines: LeftWithLastLine,
    SpaceBeforeRangeBasedForLoopColon: false,
    PointerAlignment: Left,
    NamespaceIndentation: Inner,
    AlwaysBreakTemplateDeclarations: No,
    BreakBeforeConceptDeclarations: Never,
    MaxEmptyLinesToKeep: 0
  }' "$f"
done




