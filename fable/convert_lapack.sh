#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    echo "Usage: $0 <fortran_file> [lin|eig|matgen]" >&2
    exit 1
fi

src="$1"
#echo "DBG: src='$1' mode='${2-}' script_dir='${script_dir-}'" >&2
mode="${2:-}"

# Directory of this script (used to find headers and name maps)
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${script_dir}/clang_format_common.sh"
core_name_map="${script_dir}/mplapack_name_map.txt"
testing_name_map="${script_dir}/mplapack_testing_name_map.txt"

# Expose both maps to the converter (harmless if ignored).
# This is a future-proof hook for fable.command_line.cout / cout.py.
name_map_files_env="${core_name_map}:${testing_name_map}"


: "${FABLE_SMALL_CHAR:=12}"

fable_cout_env=()
case "$mode" in
  lin|eig)
    fable_cout_env=(env FABLE_SMALL_ARRAY_SIZE="690" FABLE_SMALL_CHAR="0" FABLE_SUPPRESS_COMMON=True FABLE_SUPPRESS_SAVE=1 FABLE_NAME_MAP_FILES="$name_map_files_env" FABLE_COMMON_AS_GLOBALS=1)
    ;;
  matgen)
    fable_cout_env=(env FABLE_SMALL_ARRAY_SIZE="1000" FABLE_SMALL_CHAR="0" FABLE_SUPPRESS_COMMON=True FABLE_SUPPRESS_SAVE=1 FABLE_NAME_MAP_FILES="$name_map_files_env" FABLE_COMMON_AS_GLOBALS=1)
    ;;
  *)
    fable_cout_env=(env FABLE_SMALL_ARRAY_SIZE="10000" FABLE_SMALL_CHAR="$FABLE_SMALL_CHAR" FABLE_NAME_MAP_FILES="$name_map_files_env" FABLE_COMMON_AS_GLOBALS=1)
    ;;
esac

case "$mode" in
    "")
        header="${script_dir}/../mplapack/reference/mplapack.h.in"
        ;;
    lin)
        header="${script_dir}/../mplapack/test/lin/common/mplapack_lin.h.in"
        ;;
    eig)
        header="${script_dir}/../mplapack/test/eig/common/mplapack_eig.h.in"
        ;;
    matgen)
        header="${script_dir}/../mplapack/test/matgen/mplapack_matgen.h.in"
        ;;
    *)
        echo "Error: unknown mode '$mode' (allowed: lin, eig, matgen)" >&2
        exit 2
        ;;
esac

# Optional but recommended: fail early if the selected header is missing
if [ ! -f "$header" ]; then
    echo "Error: header not found: $header" >&2
    exit 2
fi

# MPLAPACK repository root (one level above this script)
mplapack_root="$(cd "${script_dir}/.." && pwd)"

# Make sure Python can import the 'fable' package from MPLAPACK root.
# This assumes the 'fable' package is located at ${mplapack_root}/fable.
export PYTHONPATH="${mplapack_root}:${PYTHONPATH:-}"

if [ ! -f "$header" ]; then
    echo "Error: header file not found: $header" >&2
    exit 1
fi

# Enforce that both name maps exist next to this script.
if [ ! -f "$core_name_map" ]; then
    echo "Error: MPLAPACK core name map not found: $core_name_map" >&2
    exit 1
fi
if [ ! -f "$testing_name_map" ]; then
    echo "Error: MPLAPACK testing name map not found: $testing_name_map" >&2
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

# Decide which name map should take precedence.
# - For testing conversions (lin/eig/matgen or sources under .../TESTING/...),
#   prefer the testing map first.
# - Otherwise prefer the core map first.
is_testing=false
case "$mode" in
    lin|eig|matgen) is_testing=true ;;
esac
if [[ "$src_abs" == *"/TESTING/"* ]]; then
    is_testing=true
fi

map_files=()
if $is_testing; then
    map_files=("$testing_name_map" "$core_name_map")
else
    map_files=("$core_name_map" "$testing_name_map")
fi

# Map Fortran name to MPLAPACK C++ name using BOTH name maps.
# The first match wins (based on the precedence above).
mapped_name=""
for map_file in "${map_files[@]}"; do
    while IFS= read -r line; do
        # Strip comments starting with #
        line="${line%%#*}"
        # Trim leading/trailing whitespace
        line="$(echo "$line" | sed -e "s/^[[:space:]]*//" -e "s/[[:space:]]*$//")"
        [ -z "$line" ] && continue
        set -- $line
        src_name="$1"
        dst_name="$2"
        if [ "${src_name,,}" = "$lower_fortran" ]; then
            mapped_name="$dst_name"
            break
        fi
    done < "$map_file"
    if [ -n "$mapped_name" ]; then
        break
    fi
done

# If not found in the map, apply the default MPLAPACK rule:
#   s/d -> R, c/z -> C, a -> A
# Other leading letters must be listed in either name map file.
if [ -z "$mapped_name" ]; then
    first_char="${lower_fortran:0:1}"
    rest_name="${lower_fortran:1}"
    case "$first_char" in
        s|d) mapped_name="R${rest_name}" ;;
        c|z) mapped_name="C${rest_name}" ;;
        a)   mapped_name="A${rest_name}" ;;
        *)
            echo "Error: routine '${fortran_name}' is not in ${core_name_map} or ${testing_name_map} and has no default mapping rule" >&2
            exit 1
            ;;
    esac
fi

cpp_generated="${src_dir}/${mapped_name}.cpp"

# Temporary files
tmp_body="$(mktemp)"
tmp_cpp="$(mktemp)"

# Run fable cout in the script directory so that cout.py can load
# mplapack_name_map.txt and mplapack_testing_name_map.txt from the same directory as this script.
(
    cd "$script_dir"
    "${fable_cout_env[@]}" python -m fable.command_line.cout "$src_abs" > /dev/null
)

# Ensure the expected generated C++ file exists.
if [ ! -f "$cpp_generated" ]; then
    echo "Error: expected generated C++ file not found: $cpp_generated" >&2
    exit 1
fi

# Prepend MPLAPACK LAPACK header
cat "$header" "$cpp_generated" > "$tmp_cpp"

python3 "${script_dir}/strip_boilerplate_comments.py" --inplace "$tmp_cpp"
python3 "${script_dir}/add_attribution.py" --inplace "$tmp_cpp" "$src"

# Format with clang-format (C++ indentation and style)
fable_clang_format_inplace "$tmp_cpp"

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

    # (i) -> i, (0) -> 0 when the whole index is just a single
    # parenthesized identifier or integer literal.
    m_simple = re.fullmatch(
              r"\(\s*([A-Za-z_][A-Za-z0-9_]*|\d+)\s*\)", e)
    if m_simple:
        e = m_simple.group(1)

    
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
    e = re.sub(r"\b0\s*\*\s*[A-Za-z_][A-Za-z0-9_]*", "0", e)
    # NAME * 0 -> 0
    
    e = re.sub(r"[A-Za-z_][A-Za-z0-9_]*\s*\*\s*0\b", "0", e)
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
# 3) Replace trans/norm comparisons with Mlsame(...)
#
# Target patterns (before // comment):
#   (trans == 'N' || trans == 'n')
#   (trans == "T" || trans == "t")
#   trans == 'C' || trans == 'c'
#   norm == "1" or norm == '1'
# and the same with single or double quotes.
# ---------------------------------------------------------------------------

TRANS_EQ_PATTERNS = [
    # Norm: norm == "1" or norm == '1'
    (
        re.compile(r'\bnorm\s*==\s*(?:"1"|\'1\')'),
        'Mlsame(norm, "1")'
    ),

    # Trans: parenthesized forms, e.g.:
    #   bool notran = (trans == 'N' || trans == 'n');
    (
        re.compile(r'\(\s*trans\s*==\s*[\'"]N[\'"]\s*\|\|\s*trans\s*==\s*[\'"]n[\'"]\s*\)'),
        'Mlsame(trans, "N")'
    ),
    (
        re.compile(r'\(\s*trans\s*==\s*[\'"]T[\'"]\s*\|\|\s*trans\s*==\s*[\'"]t[\'"]\s*\)'),
        'Mlsame(trans, "T")'
    ),
    (
        re.compile(r'\(\s*trans\s*==\s*[\'"]C[\'"]\s*\|\|\s*trans\s*==\s*[\'"]c[\'"]\s*\)'),
        'Mlsame(trans, "C")'
    ),

    # Trans: plain forms, e.g.:
    #   } else if (trans == "T" || trans == "t") {
    (
        re.compile(r'trans\s*==\s*[\'"]N[\'"]\s*\|\|\s*trans\s*==\s*[\'"]n[\'"]'),
        'Mlsame(trans, "N")'
    ),
    (
        re.compile(r'trans\s*==\s*[\'"]T[\'"]\s*\|\|\s*trans\s*==\s*[\'"]t[\'"]'),
        'Mlsame(trans, "T")'
    ),
    (
        re.compile(r'trans\s*==\s*[\'"]C[\'"]\s*\|\|\s*trans\s*==\s*[\'"]c[\'"]'),
        'Mlsame(trans, "C")'
    ),

    # Uplo: parenthesized OR forms, e.g.:
    #   bool upper = (uplo == 'U' || uplo == 'u');
    (
        re.compile(r'\(\s*uplo\s*==\s*[\'"]U[\'"]\s*\|\|\s*uplo\s*==\s*[\'"]u[\'"]\s*\)'),
        'Mlsame(uplo, "U")'
    ),
    (
        re.compile(r'\(\s*uplo\s*==\s*[\'"]L[\'"]\s*\|\|\s*uplo\s*==\s*[\'"]l[\'"]\s*\)'),
        'Mlsame(uplo, "L")'
    ),

    # Uplo: plain OR forms, just in case:
    #   if (uplo == "U" || uplo == "u")
    (
        re.compile(r'uplo\s*==\s*[\'"]U[\'"]\s*\|\|\s*uplo\s*==\s*[\'"]u[\'"]'),
        'Mlsame(uplo, "U")'
    ),
    (
        re.compile(r'uplo\s*==\s*[\'"]L[\'"]\s*\|\|\s*uplo\s*==\s*[\'"]l[\'"]'),
        'Mlsame(uplo, "L")'
    ),

    # Uplo: simple fallback comparisons
    (
        re.compile(r'\buplo\s*==\s*[\'"]U[\'"]'),
        'Mlsame(uplo, "U")'
    ),
    (
        re.compile(r'\buplo\s*==\s*[\'"]L[\'"]'),
        'Mlsame(uplo, "L")'
    ),
]

def _rewrite_trans_eq(line: str) -> str:
    # Do not touch C++ line comments.
    idx = line.find("//")
    if idx >= 0:
        code, comment = line[:idx], line[idx:]
    else:
        code, comment = line, ""
    for pattern, replacement in TRANS_EQ_PATTERNS:
        code = pattern.sub(replacement, code)
    return code + comment


lines = text.splitlines(keepends=True)
lines = [_rewrite_trans_eq(line) for line in lines]

# ---------------------------------------------------------------------------
# 4) Drop bogus "cabs1(v) = abs(v.real()) + abs(v.imag());" line
#    and similarly "abs1(v) = abs(v.real()) + abs(v.imag());"
# ---------------------------------------------------------------------------

CABS1_STMT_RE = re.compile(
    r"^\s*cabs1\s*\(\s*([A-Za-z_][A-Za-z0-9_]*)\s*\)\s*=\s*"
    r"abs\s*\(\s*\1\s*\.\s*real\s*\(\s*\)\s*\)\s*\+\s*"
    r"abs\s*\(\s*\1\s*\.\s*imag\s*\(\s*\)\s*\)\s*;\s*$"
)

ABS1_STMT_RE = re.compile(
    r"^\s*abs1\s*\(\s*([A-Za-z_][A-Za-z0-9_]*)\s*\)\s*=\s*"
    r"abs\s*\(\s*\1\s*\.\s*real\s*\(\s*\)\s*\)\s*\+\s*"
    r"abs\s*\(\s*\1\s*\.\s*imag\s*\(\s*\)\s*\)\s*;\s*$"
)

ABSSQ_STMT_RE = re.compile(
    r'^\s*abssq\s*\([^)]*\)\s*='
)

CABS2_STMT_RE = re.compile(
    r'^\s*cabs2\s*\([^)]*\)\s*='
)

lines = [
    line
    for line in lines
    if not CABS1_STMT_RE.match(line)
    and not CABS2_STMT_RE.match(line)
    and not ABS1_STMT_RE.match(line)
    and not ABSSQ_STMT_RE.match(line)
]

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
# Fix only conservative, obvious broken Fortran DO translations:
#
#   for (i = i1; i <= i2; i = i + i3)
#   for (j = jfirst; j <= jlast; j = j + jinc)
#   for (i = i1; i <= i2; i = i + inc)
#   for (i = ifirst; i <= ilast; i = i + iinc)
#
# but do not rewrite every variable-step loop such as:
#
#   for (krow = iloz; krow <= ihiz; krow = krow + nv)
#
# because that causes unnecessary churn in handwritten or already-correct code.
# ---------------------------------------------------------------------------

FOR_VARIABLE_STEP_RE = re.compile(
    r'\bfor\s*\(\s*'
    r'(?P<idx>[A-Za-z_]\w*)\s*=\s*(?P<first>[A-Za-z_]\w*)\s*;\s*'
    r'(?P=idx)\s*<=\s*(?P<last>[A-Za-z_]\w*)\s*;\s*'
    r'(?P=idx)\s*=\s*(?P=idx)\s*\+\s*(?P<step>[A-Za-z_]\w*)\s*'
    r'\)'
)

LIKELY_FORTRAN_DO_BOUND_RE = re.compile(
    r'^(?:'
    r'[ijklmn][12]|'
    r'[ijklmn]first|[ijklmn]last|'
    r'ifirst|ilast|jfirst|jlast|kfirst|klast|'
    r'lfirst|llast|mfirst|mlast|nfirst|nlast'
    r')$'
)

LIKELY_FORTRAN_DO_STEP_RE = re.compile(
    r'^(?:'
    r'inc|'
    r'i3|j3|k3|l3|m3|n3|'
    r'iinc|jinc|kinc|linc|minc|ninc|'
    r'istep|jstep|kstep|lstep|mstep|nstep'
    r')$'
)


def _should_fix_translated_do(first: str, last: str, step: str) -> bool:
    return (
        LIKELY_FORTRAN_DO_BOUND_RE.match(first) is not None
        and LIKELY_FORTRAN_DO_BOUND_RE.match(last) is not None
        and LIKELY_FORTRAN_DO_STEP_RE.match(step) is not None
    )


def _fix_translated_do_loop(line: str) -> str:
    # Do not touch C++ line comments.
    idx = line.find("//")
    if idx >= 0:
        code, comment = line[:idx], line[idx:]
    else:
        code, comment = line, ""

    def _repl(m: re.Match) -> str:
        idx_name = m.group("idx")
        first = m.group("first")
        last = m.group("last")
        step = m.group("step")

        if not _should_fix_translated_do(first, last, step):
            return m.group(0)

        return (
            "for ({} = {}; {} > 0 ? {} <= {} : {} >= {}; {} = {} + {})".format(
                idx_name,
                first,
                step,
                idx_name,
                last,
                idx_name,
                last,
                idx_name,
                idx_name,
                step,
            )
        )

    code = FOR_VARIABLE_STEP_RE.sub(_repl, code)
    return code + comment


lines = [_fix_translated_do_loop(line) for line in lines]

path.write_text("".join(lines))

EOF

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

# Overwrite the generated C++ file with the formatted and postprocessed version
cp "$tmp_cpp" "$cpp_generated"
