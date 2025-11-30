#!/usr/bin/env python3
"""
Simplify index expressions inside [...] in C++ source files.

If no arguments are given, all *.cpp files in the current directory are
processed in-place. If file names are given as arguments, only those
files are processed.

Rules applied inside brackets (on non-comment lines only):
  - (1 - 1)       -> 0
  - 0 * NAME      -> 0
  - NAME * 0      -> 0
  - "+ 0", "0 +", "- 0" removed
  - whitespace is collapsed

Example:
    sva[(iwork[p - 1]) - 1] -> sva[iwork[p - 1] - 1]
"""

import sys
import re
from pathlib import Path
from typing import List

BRACKET_RE = re.compile(r"\[(.*?)\]")


def simplify_expr(expr: str) -> str:
    """Simplify an index expression that may contain (1 - 1), *0, +0, etc."""
    e = expr
    # (1 - 1) -> 0
    e = re.sub(r"\(\s*1\s*-\s*1\s*\)", "0", e)
    # 0 * NAME -> 0
    e = re.sub(r"\b0\s*\*\s*[A-Za-z_][A-Za-z0-9_]*", "0", e)
    # NAME * 0 -> 0
    e = re.sub(r"[A-Za-z_][A-Za-z0-9_]*\s*\*\s*0\b", "0", e)
    # "+ 0", "0 +", "- 0" -> remove
    e = re.sub(r"\+\s*0\b", "", e)
    e = re.sub(r"\b0\s*\+", "", e)
    e = re.sub(r"\-\s*0\b", "", e)
    # Collapse spaces
    e = re.sub(r"\s+", " ", e).strip()
    return e


def simplify_line(line: str, in_block_comment: bool) -> (str, bool):
    """Simplify indices on a single line, respecting comments.

    - If inside a block comment (/* ... */), the line is returned unchanged.
    - If the line starts with '//' after indentation, it is returned unchanged.
    - Otherwise, simplify expressions inside [...] on this line only.
    """
    stripped = line.lstrip()

    # Handle block comment start/end
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

    # Apply simplification only to this line (no DOTALL)
    def repl(m: re.Match) -> str:
        inner = m.group(1)
        return "[" + simplify_expr(inner) + "]"

    new_line = BRACKET_RE.sub(repl, line)
    return new_line, in_block_comment


def process_file(path: Path) -> None:
    """Read a file, simplify index expressions, and write back in-place."""
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True)
    out: List[str] = []

    in_block_comment = False
    for line in lines:
        new_line, in_block_comment = simplify_line(line, in_block_comment)
        out.append(new_line)

    path.write_text("".join(out), encoding="utf-8")


def main(argv: List[str]) -> None:
    if len(argv) > 1:
        files = [Path(p) for p in argv[1:]]
    else:
        files = list(Path(".").glob("*.cpp"))

    if not files:
        print("simplify_index_expr: no files to process", file=sys.stderr)
        return

    for p in files:
        if not p.is_file():
            print(f"simplify_index_expr: {p} is not a file, skipping", file=sys.stderr)
            continue
        process_file(p)


if __name__ == "__main__":
    main(sys.argv)
