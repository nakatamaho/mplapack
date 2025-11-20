#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <fortran_file>" >&2
    exit 1
fi

src="$1"

# Temporary file for generated C++ code
tmp_cpp="$(mktemp)"

# Run fable cout, normalize comments, drop leading comment block,
# then format with clang-format
PYTHONPATH="$HOME/mplapack" \
  python -m fable.command_line.cout "$src" \
  | sed 's|//C|//|g' \
  | python -c 'import sys
started = False
for line in sys.stdin:
    if not started:
        # Skip leading blank lines
        if line.strip() == "":
            continue
        # Skip leading comment lines (// ...)
        if line.lstrip().startswith("//"):
            continue
        # First non-comment, non-blank line: start output from here
        started = True
    sys.stdout.write(line)
' > "$tmp_cpp"

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
  }' "$tmp_cpp"

# Print formatted code to stdout
cat "$tmp_cpp"

# Clean up
rm -f "$tmp_cpp"
