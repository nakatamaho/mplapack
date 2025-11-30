"""
Strip netlib LAPACK boilerplate comment *lines* from C++ / header files.

Usage:
    strip_lapack_comments.py FILE [FILE ...]
"""

import sys
from pathlib import Path

# Keywords that identify LAPACK boilerplate header lines.
LAPACK_HEADER_KEYWORDS = [
    "-- LAPACK computational routine --",
    "-- LAPACK auxiliary routine --",
    "LAPACK is a software package provided by",
    "LAPACK is a software package provided by",
    "Univ. of Tennessee",
    "Univ. of California Berkeley",
    "Univ. of Colorado Denver",
    "NAG Ltd",
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


def is_lapack_boilerplate_text(text: str) -> bool:
    """Return True if this comment text is LAPACK boilerplate.

    We delete only lines that clearly match LAPACK headers or section labels,
    not ordinary explanatory comments like 'Test the input parameters.'.
    """
    # Header lines
    for kw in LAPACK_HEADER_KEYWORDS:
        if kw in text:
            return True

    # Section markers like '.. Scalar Arguments ..'
    for kw in LAPACK_SECTION_PATTERNS:
        if text == kw:
            return True

    # Separator line "====..."
    if "====" in text:
        return True

    # Standalone '..' lines that belong to these sections
    if text == "..":
        return True

    return False


def strip_lapack_comments(path: Path) -> None:
    """Remove LAPACK boilerplate comment *lines* from a file in-place."""
    try:
        src = path.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True)
    except OSError as e:
        print(f"strip_lapack_comments: failed to read {path}: {e}", file=sys.stderr)
        return

    out_lines = []

    for line in src:
        if not is_comment_line(line):
            out_lines.append(line)
            continue

        # Extract comment text after '//'
        stripped = line.lstrip()
        text = stripped[2:].strip()

        # If this line is LAPACK boilerplate, skip it.
        if is_lapack_boilerplate_text(text):
            continue

        # Otherwise keep the comment line.
        out_lines.append(line)

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
