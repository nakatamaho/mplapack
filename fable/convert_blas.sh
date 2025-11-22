#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <fortran_file>" >&2
    exit 1
fi

src="$1"

# Directory of this script (used to find header_blas.txt and mplapack_name_map.txt)
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
header="${script_dir}/header_blas.txt"
name_map="${script_dir}/mplapack_name_map.txt"

# MPLAPACK repository root (one level above this script)
mplapack_root="$(cd "${script_dir}/.." && pwd)"

# Make sure Python can import the 'fable' package from MPLAPACK root.
# This assumes the 'fable' package is located at ${mplapack_root}/fable.
export PYTHONPATH="${mplapack_root}:${PYTHONPATH:-}"

if [ ! -f "$header" ]; then
    echo "Error: header file not found: $header" >&2
    exit 1
fi

# Enforce that mplapack_name_map.txt exists next to this script.
if [ ! -f "$name_map" ]; then
    echo "Error: MPLAPACK name map not found: $name_map" >&2
    exit 1
fi

# Resolve absolute path of the source Fortran file.
if [[ "$src" = /* ]]; then
    src_abs="$src"
else
    src_dir="$(cd "$(dirname "$src")" && pwd)"
    src_base="$(basename "$src")"
    src_abs="${src_dir}/${src_base}"
fi

# Temporary files
tmp_body="$(mktemp)"
tmp_cpp="$(mktemp)"

# Run fable cout in the script directory so that cout.py can load
# mplapack_name_map.txt from the same directory as this script.
# Then strip leading blank lines and leading // comments using the
# Python filter you specified.
(
    cd "$script_dir"
    python -m fable.command_line.cout "$src_abs"
) | python -c 'import sys
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

# Prepend MPLAPACK BLAS header
cat "$header" "$tmp_body" > "$tmp_cpp"

# Format with clang-format (C++ indentation and style)
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
