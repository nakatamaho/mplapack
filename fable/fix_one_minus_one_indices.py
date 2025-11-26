#!/usr/bin/env python3
import sys
import re
from pathlib import Path

def simplify_index_expr(expr: str) -> str:
    """Simplify an index expression that may contain (1 - 1)."""
    e = expr

    # (1 - 1) -> 0
    e = re.sub(r"\(\s*1\s*-\s*1\s*\)", "0", e)

    # 0 * name -> 0
    e = re.sub(r"\b0\s*\*\s*[A-Za-z_][A-Za-z0-9_]*", "0", e)
    # name * 0 -> 0
    e = re.sub(r"[A-Za-z_][A-Za-z0-9_]*\s*\*\s*0\b", "0", e)

    # + 0 / 0 + / - 0 -> remove
    e = re.sub(r"\+\s*0\b", "", e)
    e = re.sub(r"\b0\s*\+", "", e)
    e = re.sub(r"\-\s*0\b", "", e)

    # Collapse multiple spaces
    e = re.sub(r"\s+", " ", e).strip()
    return e

def fix_file(path: Path) -> None:
    text = path.read_text()

    # Simplify inside [...] (indexes), allowing newlines inside.
    pat_bracket = re.compile(r"\[(.*?)\]", re.DOTALL)

    def repl(m):
        inner = m.group(1)
        simplified = simplify_index_expr(inner)
        return "[" + simplified + "]"

    new_text = pat_bracket.sub(repl, text)

    # Optional: global cleanup of leftover (1 - 1), 0*X, X*0, +0, 0+
    new_text = re.sub(r"\(\s*1\s*-\s*1\s*\)", "0", new_text)
    new_text = re.sub(r"\b0\s*\*\s*[A-Za-z_][A-Za-z0-9_]*", "0", new_text)
    new_text = re.sub(r"[A-Za-z_][A-Za-z0-9_]*\s*\*\s*0\b", "0", new_text)
    new_text = re.sub(r"\+\s*0\b", "", new_text)
    new_text = re.sub(r"\b0\s*\+", "", new_text)
    new_text = re.sub(r"\-\s*0\b", "", new_text)

    path.write_text(new_text)

def main(argv):
    if len(argv) < 2:
        print("Usage: fix_one_minus_one_indices.py file1.cpp [file2.cpp ...]", file=sys.stderr)
        sys.exit(1)

    for name in argv[1:]:
        p = Path(name)
        if not p.is_file():
            print(f"Skipping {name} (not a file)", file=sys.stderr)
            continue
        fix_file(p)

if __name__ == "__main__":
    main(sys.argv)
