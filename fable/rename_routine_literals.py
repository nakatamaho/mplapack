#!/usr/bin/env python3
import argparse
import os
import re
from pathlib import Path
from typing import Dict, Tuple, List, Iterable


def load_name_map(map_path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with map_path.open("r", encoding="utf-8", errors="surrogateescape") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 2:
                continue
            src, dst = parts
            mapping[src.lower()] = dst
    return mapping


def _split_map_arg(map_arg: str) -> List[str]:
    """
    Allow either:
      --map A --map B
    or legacy convenience:
      --map A,B
    """
    parts: List[str] = []
    for x in map_arg.split(","):
        x = x.strip()
        if x:
            parts.append(x)
    return parts


def load_name_maps(map_args: Iterable[str]) -> Dict[str, str]:
    """
    Load and merge multiple mapping files.
    Precedence rule: later files override earlier ones.
    """
    merged: Dict[str, str] = {}
    for arg in map_args:
        for one in _split_map_arg(arg):
            map_path = Path(os.path.expanduser(one))
            if not map_path.exists():
                raise FileNotFoundError(f"map file not found: {map_path}")
            m = load_name_map(map_path)
            merged.update(m)
    return merged


def map_token(token_with_spaces: str, mapping: Dict[str, str]) -> str | None:
    """
    Map only by mapping files.
    - Trim trailing spaces before lookup (e.g., "ZGGES " -> "zgges").
    - No fallback rules (D/S->R, C/Z->C are disabled).
    """
    trimmed = token_with_spaces.rstrip(" ")
    if not trimmed:
        return None

    key = trimmed.lower()
    return mapping.get(key, None)


# Token finder inside string literals:
# - begins with D/S/C/Z
# - followed by 4..20 chars of [A-Z0-9_]
# - optional trailing spaces (to support "ZGGES " patterns)
# - bounded by non-word char or end
TOKEN_RE = re.compile(
    r'(?<![A-Za-z0-9_])([DSCZ][A-Z0-9_]{4,20})([ ]*)(?=[^A-Za-z0-9_]|$)'
)

# Whole-literal routine tag (optionally padded with trailing spaces).
# Used for srnamt / Chkxer("XXXX  ", ...) style strings: padding is meaningless.
TOKEN_ONLY_RE = re.compile(r'^([DSCZ][A-Za-z0-9_]{4,20})[ ]*$')

# C/C++ normal string literal matcher (handles \" and \\).
# NOTE: This is intentionally *not* a full C++ lexer; it is sufficient for
# fable-generated code where string literals are emitted as normal literals.
CPP_STRING_RE = re.compile(r'"((?:\\.|[^"\\])*)"')

# Match format_#### assignments that may use adjacent string literals:
#   format_9999 = "...."
#   format_9999 = "...."
#                "....";
FORMAT_ASSIGN_RE = re.compile(
    r'(\bformat_\d+\s*=\s*)(?P<lits>(?:"(?:\\.|[^"\\])*"\s*)+);',
    flags=re.DOTALL,
)

# Match fem-style write statements using the comma operator:
#   write(unit, fmt), "Left", "DGGEV1", ...;
# Capture everything after the first comma (after the write(...) call)
# up to the terminating semicolon.
WRITE_STMT_RE = re.compile(
    r'(\bwrite\s*\(\s*[^)]*?\)\s*,)(?P<tail>.*?);',
    flags=re.DOTALL,
)


def rewrite_literal_content(content: str, mapping: Dict[str, str]) -> Tuple[str, int]:
    """
    Replace routine tokens found inside a single C++ string literal content,
    but only if the token exists in mapping.
    """
    # Special case: the literal is exactly one routine token (possibly padded).
    # Examples:
    #   Chkxer("Cgesvdq ", ...)
    #   srnamt = "Cgghrd ";
    # In these cases, trailing spaces are Fortran fixed-length leftovers and
    # should be removed.
    m0 = TOKEN_ONLY_RE.fullmatch(content)
    if m0:
        tok = m0.group(1)
        mapped = mapping.get(tok.lower())
        if mapped is not None:
            return mapped, 1
        # Even if already renamed / not found in mapping, still drop padding.
        return tok, 0

    hits = 0

    def repl(m: re.Match) -> str:
        nonlocal hits
        base = m.group(1)
        spaces = m.group(2)

        mapped = map_token(base + spaces, mapping)
        if mapped is None:
            return m.group(0)

        hits += 1
        # Preserve original padding spaces.
        # LAPACK tests often rely on fixed-width routine tags and aligned banners
        # (e.g., "XXXXXX " or "ZLANHF         ***" inside FORMAT strings).
        return mapped + spaces

    new_content = TOKEN_RE.sub(repl, content)
    return new_content, hits


def apply_patterns(text: str, mapping: Dict[str, str]) -> Tuple[str, int]:
    """
    Targeted replacements only (avoid rewriting unrelated strings):
      - srnamt = "...."
      - write(unit, fmt), "...."
      - write(unit, fmt), ...;   (rewrite routine tokens inside any string literals)
      - format_#### = "...."  (possibly multiple adjacent literals)
      - chkxer/Chkxer/chkxer2/... (case-insensitive)
      - Mxerbla/Mxerbla_array
      - Alaerh/Alahd/Alaesm ... , "...."
    """
    total_hits = 0

    patterns: List[re.Pattern] = [
        re.compile(r'(\bsrnamt\s*=\s*")([^"]*)(")'),

        re.compile(r'(\bwrite\s*\(\s*[^,]*,\s*")([^"]*)(")'),

        # chkxer vs Chkxer: handle both
        re.compile(r'(\bchkxer\d*\s*\(\s*")([^"]*)(")', flags=re.IGNORECASE),

        re.compile(r'(\bMxerbla(?:_array)?\s*\(\s*")([^"]*)(")'),

        re.compile(r'(\bAlaerh\d*\s*\(\s*[^,]*,\s*")([^"]*)(")'),
        re.compile(r'(\bAlahd\s*\(\s*[^,]*,\s*")([^"]*)(")'),
        re.compile(r'(\bAlaesm\s*\(\s*[^,]*,\s*")([^"]*)(")'),
    ]

    for pat in patterns:

        def outer_repl(m: re.Match) -> str:
            nonlocal total_hits
            content = m.group(2)
            new_content, hits = rewrite_literal_content(content, mapping)
            if new_content == content:
                return m.group(0)
            total_hits += hits
            return m.group(1) + new_content + m.group(3)

        text = pat.sub(outer_repl, text)

    # format_#### assignments: rewrite routine tokens inside the format strings.
    def format_repl(m: re.Match) -> str:
        nonlocal total_hits
        prefix = m.group(1)
        lits = m.group("lits")

        def one_lit(sm: re.Match) -> str:
            nonlocal total_hits
            content = sm.group(1)
            new_content, hits = rewrite_literal_content(content, mapping)
            total_hits += hits
            return '"' + new_content + '"'

        new_lits = CPP_STRING_RE.sub(one_lit, lits)
        return prefix + new_lits + ";"

    text = FORMAT_ASSIGN_RE.sub(format_repl, text)

    # write(...), ...; statements: rewrite routine tokens inside any string literals
    # in the comma-operator argument list.
    def write_stmt_repl(m: re.Match) -> str:
        nonlocal total_hits
        prefix = m.group(1)
        tail = m.group("tail")

        def one_lit(sm: re.Match) -> str:
            nonlocal total_hits
            content = sm.group(1)
            new_content, hits = rewrite_literal_content(content, mapping)
            total_hits += hits
            return '"' + new_content + '"'

        new_tail = CPP_STRING_RE.sub(one_lit, tail)
        return prefix + new_tail + ";"

    text = WRITE_STMT_RE.sub(write_stmt_repl, text)
    return text, total_hits


def iter_files(paths: List[Path], suffixes: Tuple[str, ...]) -> List[Path]:
    out: List[Path] = []
    for p in paths:
        if p.is_file():
            if p.suffix in suffixes:
                out.append(p)
            continue
        if p.is_dir():
            for root, _, files in os.walk(p):
                for fn in files:
                    fp = Path(root) / fn
                    if fp.suffix in suffixes:
                        out.append(fp)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Rename LAPACK routine tokens inside string literals using mapping txt files (map-only)"
    )
    ap.add_argument(
        "--map",
        action="append",
        required=True,
        help="Path to map txt (repeatable). Later maps override earlier ones.",
    )
    ap.add_argument("--in-place", action="store_true",
                    help="Rewrite files in place")
    ap.add_argument(
        "--suffix",
        action="append",
        default=[".cpp", ".hpp", ".h"],
        help="File suffix to process (repeatable)",
    )
    ap.add_argument("paths", nargs="*",
                    default=["."], help="Files/directories to scan")
    args = ap.parse_args()

    mapping = load_name_maps(args.map)

    targets = [Path(os.path.expanduser(x)) for x in args.paths]
    suffixes = tuple(args.suffix)

    files = iter_files(targets, suffixes)
    changed_files = 0
    changed_hits = 0

    for fp in files:
        old = fp.read_text(encoding="utf-8", errors="surrogateescape")
        new, hits = apply_patterns(old, mapping)
        if hits and new != old:
            changed_files += 1
            changed_hits += hits
            if args.in_place:
                fp.write_text(new, encoding="utf-8", errors="surrogateescape")

    mode = "in-place" if args.in_place else "dry-run"
    print(f"[{mode}] changed_files={changed_files} replaced_tokens={changed_hits}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
