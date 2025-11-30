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
    src_dir_tmp="$(cd "$(dirname "$src")" && pwd)"
    src_base_tmp="$(basename "$src")"
    src_abs="${src_dir_tmp}/${src_base_tmp}"
fi

# Derive source directory and base name.
src_dir="$(cd "$(dirname "$src_abs")" && pwd)"
src_base="$(basename "$src_abs")"

# Fortran routine name is assumed to be the basename without extension.
fortran_name="${src_base%.*}"
lower_fortran="${fortran_name,,}"

# Map Fortran name to MPLAPACK C++ name using mplapack_name_map.txt first.
mapped_name=""
while IFS= read -r line; do
    # Strip comments starting with #
    line="${line%%#*}"
    # Trim leading/trailing whitespace
    line="$(echo "$line" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
    [ -z "$line" ] && continue
    set -- $line
    src_name="$1"
    dst_name="$2"
    if [ "${src_name,,}" = "$lower_fortran" ]; then
        mapped_name="$dst_name"
        break
    fi
done < "$name_map"

# If not found in the map, apply the default MPLAPACK rule: s/d -> R, c/z -> C.
if [ -z "$mapped_name" ]; then
    first_char="${lower_fortran:0:1}"
    tail="${lower_fortran:1}"
    case "$first_char" in
        s|d)
            mapped_name="R${tail}"
            ;;
        c|z)
            mapped_name="C${tail}"
            ;;
        *)
            mapped_name="$fortran_name"
            ;;
    esac
fi

cpp_generated="${src_dir}/${mapped_name}.cpp"

# Temporary files
tmp_body="$(mktemp)"
tmp_cpp="$(mktemp)"

# Run fable cout in the script directory so that cout.py can load
# mplapack_name_map.txt from the same directory as this script.
(
    cd "$script_dir"
    python -m fable.command_line.cout "$src_abs" > /dev/null
)

# Ensure the expected generated C++ file exists.
if [ ! -f "$cpp_generated" ]; then
    echo "Error: expected generated C++ file not found: $cpp_generated" >&2
    exit 1
fi

# Strip leading blank lines and leading // comments from the generated C++ file.
python -c 'import sys
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
' < "$cpp_generated" > "$tmp_body"

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
    PointerAlignment: Right,
    NamespaceIndentation: Inner,
    AlwaysBreakTemplateDeclarations: No,
    BreakBeforeConceptDeclarations: Never,
  }' "$tmp_cpp"


python3 "${script_dir}/strip_lapack_comments.py" "$tmp_cpp"

#python3 "${script_dir}/normalize_comment_prefix.py" "$tmp_cpp"

# Overwrite the generated C++ file with the formatted version
cp "$tmp_cpp" "$cpp_generated"

# Clean up
rm -f "$tmp_body" "$tmp_cpp"
