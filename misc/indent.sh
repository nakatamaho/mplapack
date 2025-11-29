for f in *.cpp *.cc *.h *.hpp *.h.in; do

python3 ~/mplapack/misc/strip_lapack_comments.py "$f" 
python3 ~/mplapack/misc/normalize_comment_prefix.py "$f"

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
done
