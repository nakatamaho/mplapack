#!/usr/bin/env python3
"""
Strip netlib LAPACK boilerplate comment blocks from C++ / header files.

Usage:
    strip_lapack_comments.py FILE [FILE ...]
"""

import sys
from pathlib import Path

# Keywords that identify LAPACK boilerplate header lines.
LAPACK_HEADER_KEYWORDS = [
    "-- LAPACK auxiliary routine --",
    "LAPACK is a software package provided by",
]

# Exact section labels from netlib LAPACK sources.
LAPACK_SECTION_PATTERNS = [
    ".. Scalar Arguments ..",
    ".. Array Arguments ..",
    ".. Parameters ..",
    ".. Local Scalars ..",
    ".. External Functions ..",
    ".. External Subroutines ..",
    ".. Intrinsic Functions ..",
    ".. Executable Statements ..",
]

def is_comment_line(line: str) -> bool:
    """Return True if the line is a C++ line comment starting with '//'."""
    stripped = line.lstrip()
    return stripped.startswith("//")


def is_lapack_boilerplate_block(block_lines) -> bool:
    """Return True if this comment block looks like netlib LAPACK boilerplate.

    We only delete blocks that clearly match LAPACK headers or section labels,
    not ordinary explanatory comments like 'Test the input parameters.'.
    """
    has_header = False
    has_section = False
    has_separator = False

    for line in block_lines:
        # Remove leading indentation and the '//' prefix
        stripped = line.lstrip()
        if not stripped.startswith("//"):
            continue
        text = stripped[2:].strip()

        # Header lines (top of file)
        for kw in LAPACK_HEADER_KEYWORDS:
            if kw in text:
                has_header = True
                break

        # Section markers like '.. Scalar Arguments ..'
        for kw in LAPACK_SECTION_PATTERNS:
            if text == kw:
                has_section = True
                break

        # Separator line "====..."
        if "====" in text:
            has_separator = True

    # Delete only clearly boilerplate sections, not generic comments.
    return has_header or has_section or has_separator


def strip_lapack_comments(path: Path) -> None:
    """Remove LAPACK boilerplate comment blocks from a file in-place."""
    try:
        src = path.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True)
    except OSError as e:
        print(f"strip_lapack_comments: failed to read {path}: {e}", file=sys.stderr)
        return

    out_lines = []
    i = 0
    n = len(src)

    while i < n:
        line = src[i]

        if not is_comment_line(line):
            # Not a comment line: keep as-is.
            out_lines.append(line)
            i += 1
            continue

        # Collect a contiguous block of line comments.
        block_start = i
        block = []
        while i < n and is_comment_line(src[i]):
            block.append(src[i])
            i += 1

        # Decide whether this block is LAPACK boilerplate.
        if is_lapack_boilerplate_block(block):
            # Skip entire block.
            continue

        # Keep the whole block.
        out_lines.extend(block)

    try:
        path.write_text("".join(out_lines), encoding="utf-8")
    except OSError as e:
        print(f"strip_lapack_comments: failed to write {path}: {e}", file=sys.stderr)


def main(argv):
    if len(argv) < 2:
        print("Usage: strip_lapack_comments.py FILE [FILE...]", file=sys.stderr)
        sys.exit(1)

    for arg in argv[1:]:
        p = Path(arg)
        if not p.is_file():
            print(f"strip_lapack_comments: {p} is not a file, skipping", file=sys.stderr)
            continue
        strip_lapack_comments(p)


if __name__ == "__main__":
    main(sys.argv)
