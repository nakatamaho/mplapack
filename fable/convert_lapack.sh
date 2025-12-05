#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <fortran_file>" >&2
    exit 1
fi

src="$1"

# Directory of this script (used to find header_blas.txt and mplapack_name_map.txt)
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
header="${script_dir}/header_lapack.txt"
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

python3 "${script_dir}/strip_boilerplate_comments.py" "$tmp_cpp"
#python3 "${script_dir}/add_attribution.py" "$tmp_cpp" "$src"

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

# Simplify trivial (1 - 1) index arithmetic in the formatted C++ file
python3 - "$tmp_cpp" << 'EOF'
import sys
import re
from pathlib import Path

path = Path(sys.argv[1])
text = path.read_text()

# ---------------------------------------------------------------------------
# 1) Line-based simplification inside [...] while respecting comments
#    (DO NOT touch comments or block comments).
# ---------------------------------------------------------------------------

lines = text.splitlines(keepends=True)

BRACKET_RE = re.compile(r"\[(.*?)\]")

def simplify_expr(expr: str) -> str:
    """Simplify an index expression that may contain (1 - 1), *0, +0, etc."""
    e = expr
    
    # -------------------------------------------------------------------------
    # STEP 1: Handle multiplications FIRST (before touching + (1-1) patterns)
    #         This prevents "+ (1 - 1) * lda" from becoming "* lda"
    # -------------------------------------------------------------------------
    
    # (1 - 1) * NAME or NAME * (1 - 1) -> 0
    e = re.sub(
        r"\(\s*1\s*-\s*1\s*\)\s*\*\s*[A-Za-z_][A-Za-z0-9_]*",
        "0",
        e,
    )
    e = re.sub(
        r"[A-Za-z_][A-Za-z0-9_]*\s*\*\s*\(\s*1\s*-\s*1\s*\)",
        "0",
        e,
    )
    
    # (1 - 1) -> 0 (standalone, after multiplication patterns are handled)
    e = re.sub(r"\(\s*1\s*-\s*1\s*\)", "0", e)
    
    # 0 * NAME -> 0
    e = re.sub(r"0\s*\*\s*[A-Za-z_][A-Za-z0-9_]*", "0", e)
    # NAME * 0 -> 0
    e = re.sub(r"[A-Za-z_][A-Za-z0-9_]*\s*\*\s*0", "0", e)
    
    # -------------------------------------------------------------------------
    # STEP 2: Now handle addition/subtraction of zero
    # -------------------------------------------------------------------------
    
    # Remove "0 + " at the beginning or after opening paren
    e = re.sub(r"^\s*0\s*\+\s*", "", e)
    e = re.sub(r"\(\s*0\s*\+\s*", "(", e)
    
    # Remove "+ 0" at the end or before closing paren
    e = re.sub(r"\s*\+\s*0\s*$", "", e)
    e = re.sub(r"\s*\+\s*0\s*\)", ")", e)
    
    # Remove "- 0" patterns
    e = re.sub(r"\s*\-\s*0\s*$", "", e)
    e = re.sub(r"\s*\-\s*0\s*\)", ")", e)
    
    # Clean up leading/trailing operators and spaces
    e = re.sub(r"^\s*\+\s*", "", e)  # Leading +
    e = re.sub(r"\s*\+\s*$", "", e)  # Trailing +
    
    # Collapse spaces
    e = re.sub(r"\s+", " ", e).strip()
    
    return e

def simplify_line(line: str, in_block_comment: bool):
    """Simplify indices on a single line, respecting comments.

    - If inside a block comment (/* ... */), return line unchanged.
    - If the line is a line comment (//...), return unchanged.
    - Otherwise, simplify expressions inside [...] on this line.
    """
    stripped = line.lstrip()

    # Handle block comment state
    if in_block_comment:
        if "*/" in stripped:
            in_block_comment = False
        return line, in_block_comment

    if stripped.startswith("/*"):
        if "*/" not in stripped:
            in_block_comment = True
        return line, in_block_comment

    # Single-line comment: leave unchanged
    if stripped.startswith("//"):
        return line, in_block_comment

    # Apply simplification only to this line
    def repl(m: re.Match) -> str:
        inner = m.group(1)
        return "[" + simplify_expr(inner) + "]"

    new_line = BRACKET_RE.sub(repl, line)
    return new_line, in_block_comment

out_lines = []
in_block_comment = False
for line in lines:
    new_line, in_block_comment = simplify_line(line, in_block_comment)
    out_lines.append(new_line)

text = "".join(out_lines)

# ---------------------------------------------------------------------------
# 2) Global simplification of double-parenthesized MIN/MAX across newlines:
#      ((min(i + 2, n)) - 1)   -> (min(i + 2, n) - 1)
#      ((MAX(i, j)) - 1)       -> (MAX(i, j) - 1)
#    This does NOT touch comments, because it matches only MIN/MAX(...) calls.
# ---------------------------------------------------------------------------

pattern_minmax = re.compile(
    r"\(\s*\(\s*("
    r"(?:[Mm][Ii][Nn]|[Mm][Aa][Xx])"   # MIN, min, MAX, max
    r"\("
    r"(?:[^()]*|\([^()]*\))*"          # allow one level of nested parentheses
    r"\)"
    r")\s*\)\s*-\s*1\s*\)",
    re.DOTALL,
)

text = pattern_minmax.sub(r"(\1 - 1)", text)

# ---------------------------------------------------------------------------
# 3) Replace "norm == '1'" with Mlsame(norm, "1") in code (not in comments)
# ---------------------------------------------------------------------------

NORM_EQ_ONE_RE = re.compile(r"\bnorm\s*==\s*'1'")

def _rewrite_norm_eq_one(line: str) -> str:
    # Do not touch C++ line comments.
    idx = line.find("//")
    if idx >= 0:
        code, comment = line[:idx], line[idx:]
    else:
        code, comment = line, ""
    code = NORM_EQ_ONE_RE.sub('Mlsame(norm, "1")', code)
    return code + comment

lines = text.splitlines(keepends=True)
lines = [_rewrite_norm_eq_one(line) for line in lines]

# ---------------------------------------------------------------------------
# 4) Drop bogus "cabs1(v) = abs(v.real()) + abs(v.imag());" line
# ---------------------------------------------------------------------------

CABS1_STMT_RE = re.compile(
    r"^\s*cabs1\s*\(\s*([A-Za-z_][A-Za-z0-9_]*)\s*\)\s*=\s*"
    r"abs\s*\(\s*\1\s*\.\s*real\s*\(\s*\)\s*\)\s*\+\s*"
    r"abs\s*\(\s*\1\s*\.\s*imag\s*\(\s*\)\s*\)\s*;\s*$"
)

lines = [line for line in lines if not CABS1_STMT_RE.match(line)]

# ---------------------------------------------------------------------------
# 5) Rewrite iMlaenv / iMlaenv2stage(..., "Txxxx", ...) where
#    T is a type letter:
#      Z -> C   (complex double  -> complex)
#      D -> R   (double real     -> real)
#    The rest of the routine name is lowercased:
#      "ZGEQR "        -> "Cgeqr "
#      "ZHETRD_2STAGE" -> "Chetrd_2stage"
#      "DSYTRD_2STAGE" -> "Rsytrd_2stage"
#    C++ line comments are not touched.
# ---------------------------------------------------------------------------

IMLAENV_T_RE = re.compile(
    r'(iMlaenv(?:2stage)?\s*\(\s*\d+\s*,\s*")([ZDSC])([^"]*)(\s*")'
)

def _rewrite_imlaenv_typed(line: str) -> str:
    # Skip C++ line comments
    idx = line.find("//")
    if idx >= 0:
        code, comment = line[:idx], line[idx:]
    else:
        code, comment = line, ""

    def repl(m: re.Match) -> str:
        prefix, tchar, name, suffix = m.groups()
        # Map type letter: Z->C, D->R, others unchanged
        tmap = {"Z": "C", "D": "R"}
        new_t = tmap.get(tchar, tchar)
        # name contains the rest, e.g. "GEQR " or "SYTRD_2STAGE"
        stripped = name.rstrip()
        trailing = name[len(stripped):]  # keep trailing spaces
        return prefix + new_t + stripped.lower() + trailing + suffix

    code = IMLAENV_T_RE.sub(repl, code)
    return code + comment

lines = [_rewrite_imlaenv_typed(line) for line in lines]

# ---------------------------------------------------------------------------
# 6) Rewrite lwork_* / lwrk_* variables with type letter:
#    lwork_zxxxx -> lwork_Cxxxx
#    lwork_dxxxx -> lwork_Rxxxx
#    lwrk_zxxxx  -> lwrk_Cxxxx
#    lwrk_dxxxx  -> lwrk_Rxxxx
#    (Other type letters are left unchanged.)
# ---------------------------------------------------------------------------

LWORK_T_RE = re.compile(r'\b(lwork|lwrk)_([ZzDd])([A-Za-z0-9_]+)')

def _rewrite_lwork_typed(line: str) -> str:
    # Do not touch C++ line comments
    idx = line.find("//")
    if idx >= 0:
        code, comment = line[:idx], line[idx:]
    else:
        code, comment = line, ""

    def repl(m: re.Match) -> str:
        prefix, tchar, rest = m.groups()
        # Map type letter Z/z -> C, D/d -> R
        tmap = {"Z": "C", "z": "C", "D": "R", "d": "R"}
        new_t = tmap.get(tchar, tchar)
        return f"{prefix}_{new_t}{rest}"

    code = LWORK_T_RE.sub(repl, code)
    return code + comment

lines = [_rewrite_lwork_typed(line) for line in lines]

# ---------------------------------------------------------------------------
# Fix broken DO loop translation:
#   for (i = i1; i <= i2; i = i + i3)
# should be
#   for (i = i1; i3 > 0 ? i <= i2 : i >= i2; i = i + i3)
# ---------------------------------------------------------------------------

FOR_I_I1_I2_I3_RE = re.compile(
    r'\bfor\s*\(\s*i\s*=\s*i1\s*;\s*i\s*<=\s*i2\s*;\s*i\s*=\s*i\s*\+\s*i3\s*\)'
)

def _fix_do_i1_i2_i3(line: str) -> str:
    # Do not touch C++ line comments
    idx = line.find("//")
    if idx >= 0:
        code, comment = line[:idx], line[idx:]
    else:
        code, comment = line, ""

    code = FOR_I_I1_I2_I3_RE.sub(
        'for (i = i1; i3 > 0 ? i <= i2 : i >= i2; i = i + i3)',
        code,
    )
    return code + comment

lines = [_fix_do_i1_i2_i3(line) for line in lines]

path.write_text("".join(lines))
EOF

# Overwrite the generated C++ file with the formatted and postprocessed version
cp "$tmp_cpp" "$cpp_generated"
