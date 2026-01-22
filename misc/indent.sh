for f in *.cpp *.cc *.h *.hpp *.h.in; do

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
}' "$f"

#  MaxEmptyLinesToKeep: 0
done
