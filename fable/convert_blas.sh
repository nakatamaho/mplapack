#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <fortran_file>" >&2
    exit 1
fi

src="$1"

# Directory of this script (used to find mplapack/mpblas/reference/mpblas.h.in and mplapack_name_map.txt)
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${script_dir}/clang_format_common.sh"
header="${script_dir}/../mpblas/reference/mpblas.h.in"
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

if [ -z "$mapped_name" ]; then
   echo "Error: not in the $mapped_name" 
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

python3 "${script_dir}/strip_boilerplate_comments.py" --inplace "$cpp_generated"
cat "$header" "$cpp_generated" > "$tmp_cpp"
python3 "${script_dir}/add_attribution.py" "$tmp_cpp" "$src"

# Trim Fortran-style right-padded routine names in Mxerbla("NAME ", ...).
python3 - "$tmp_cpp" << 'EOF_MXERBLA'
import re
import sys
from pathlib import Path

path = Path(sys.argv[1])
text = path.read_text()

MXERBLA_RE = re.compile(r'(Mxerbla\(\s*")([^"]*)("\s*,)')

def repl(m: re.Match) -> str:
    prefix, name, suffix = m.groups()
    return prefix + name.rstrip() + suffix

text2 = MXERBLA_RE.sub(repl, text)
if text2 != text:
    path.write_text(text2)
EOF_MXERBLA

# Format with clang-format (C++ indentation and style)
fable_clang_format_inplace "$tmp_cpp"

# Overwrite the generated C++ file with the formatted version
cp "$tmp_cpp" "$cpp_generated"

# Clean up
rm -f "$tmp_body" "$tmp_cpp"
