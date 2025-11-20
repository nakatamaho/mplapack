#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <fortran_file>" >&2
    exit 1
fi

src="$1"

# Directory of this script (to find header_blas.txt)
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
header="${script_dir}/header_blas.txt"

if [ ! -f "$header" ]; then
    echo "Error: header file not found: $header" >&2
    exit 1
fi

# Temporary files for generated C++ code
tmp_body="$(mktemp)"
tmp_cpp="$(mktemp)"

# Run fable cout, normalize comments, drop leading comment block
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
' > "$tmp_body"

# Prepend header_blas.txt
cat "$header" "$tmp_body" > "$tmp_cpp"

# Format with clang-format
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
  }' "$tmp_cpp"

# Print formatted code to stdout
cat "$tmp_cpp"

# Clean up
rm -f "$tmp_body" "$tmp_cpp"
