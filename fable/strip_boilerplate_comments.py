#!/usr/bin/env python3
"""
Strip netlib BLAS/LAPACK boilerplate comment *blocks* from C++ / header files.

Usage:
    strip_boilerplate_comments.py FILE [FILE ...]
"""

import sys
from pathlib import Path

# ---------------------------------------------------------------------------
#  Boilerplate detection (ported / unified from cout.py)
# ---------------------------------------------------------------------------

# Section keywords that appear inside ".. ... .." markers in netlib sources.
SECTION_KEYWORDS = [
    "Scalar Argument",
    "Array Argument",
    "Parameter",
    "Local Scalar",
    "External Function",
    "External Subroutine",
    "Intrinsic Function",
    "Local Array",
    "Statement Function",
    "Executable Statement",
    "Local Parameter",
    "Common Block",
    "Data Statement",
    "Equivalence",
    "Save statement",
]

# Phrases that clearly indicate LAPACK/BLAS boilerplate headers.
HEADER_KEYWORDS = [
    "-- Reference BLAS",
    "-- Reference LAPACK",
    "-- LAPACK is a software package",
    "LAPACK is a software package provided by",
    "Univ. of Tennessee",
    "Univ. of California Berkeley",
    "Univ. of Colorado Denver",
    "NAG Ltd",
]

# LAPACK routine type headers.
LAPACK_ROUTINE_HEADERS = [
    "LAPACK computational routine",
    "LAPACK driver routine",
    "LAPACK auxiliary routine",
]

# Sentinel phrases: after we see this in the header comment block, we no longer
# treat following lines as boilerplate. We also never delete these lines.
SENTINEL_PHRASES = [
    "test the input arguments",
    "test the input parameters",
]

# Lines that must *never* be considered boilerplate.
NEVER_BOILERPLATE_PHRASES = [
    "test the input arguments",
    "test the input parameters",
    "quick return",
]


def is_comment_line(line: str) -> bool:
    """Return True if the line is a C++ line comment starting with '//'."""
    stripped = line.lstrip()
    return stripped.startswith("//")


def extract_comment_text(line: str) -> str:
    """
    Extract the meaningful comment text from a C++ line comment.

    Handles both:
      // ......
      //C .....
    by stripping leading whitespace, the '//' marker, and an optional leading
    'C' used by cout.py to tag Fortran comments.
    """
    stripped = line.lstrip()
    if not stripped.startswith("//"):
        return ""

    # Remove leading // and following whitespace.
    text = stripped[2:].lstrip()

    # Normalize Fortran-style prefix: //C... -> ...
    if text.startswith("C"):
        # If it is exactly 'C' or 'C ' / 'C\t', treat it as a prefix.
        if len(text) == 1 or text[1].isspace():
            text = text[1:].lstrip()

    return text.rstrip("\n\r")


def is_blas_boilerplate_text(text: str) -> bool:
    """
    Return True if this comment *content* is BLAS/LAPACK boilerplate.

    This is the C++-file analogue of cout.py:is_blas_boilerplate_comment.
    We take plain comment text (no //, no leading 'C').
    """
    if text is None:
        return False

    t_stripped = text.strip()
    lower = t_stripped.lower()

    # 1) Never treat these lines as boilerplate.
    for phrase in NEVER_BOILERPLATE_PHRASES:
        if phrase in lower:
            return False

    # 2) Simple ".." / "..." lines.
    if t_stripped in ("..", "..."):
        return True

    # 3) BLAS/LAPACK header phrases.
    for kw in HEADER_KEYWORDS:
        if kw in t_stripped:
            return True

    # 4) LAPACK routine-type headers.
    for phrase in LAPACK_ROUTINE_HEADERS:
        if phrase in t_stripped:
            return True

    # 5) Section markers with ".." and well-known section keywords.
    if ".." in t_stripped:
        for keyword in SECTION_KEYWORDS:
            if keyword in t_stripped:
                return True

    # 6) Long separator lines made of '=' or '-'.
    if t_stripped.startswith("====") and len(t_stripped) > 20:
        return True
    if t_stripped.startswith("----") and len(t_stripped) > 20:
        return True

    # 7) Simple headings "Purpose", "Arguments".
    if t_stripped in ("Purpose", "Arguments"):
        return True

    # 8) "=== ... ===" style headings.
    if t_stripped.startswith("===") and t_stripped.endswith("==="):
        return True

    return False


def _is_empty_or_dotdot(text: str) -> bool:
    """
    Return True if the comment text is "empty padding" around a boilerplate block:
      - empty string
      - ".."
      - ".. foo .." style
    """
    t = text.strip()
    if t == "":
        return True
    if t == "..":
        return True
    if t.startswith("..") and t.endswith(".."):
        return True
    return False


def _compute_boilerplate_comment_indices(comment_texts):
    """
    Given a list of comment texts (in order), compute indices that should be
    removed as boilerplate, including surrounding padding lines.

    This mirrors produce_comments() / comment_manager._identify_boilerplate_blocks
    in cout.py, with the sentinel "Test the input arguments/parameters".
    """
    n = len(comment_texts)
    if n == 0:
        return set()

    # 1) Find sentinel index (first occurrence of "test the input arguments/parameters").
    sentinel_idx = None
    for i, text in enumerate(comment_texts):
        tl = text.strip().lower()
        for phrase in SENTINEL_PHRASES:
            if phrase in tl:
                sentinel_idx = i
                break
        if sentinel_idx is not None:
            break

    # 2) Mark obvious boilerplate lines, only *before* sentinel.
    boilerplate_indices = []
    skip_flags = [False] * n
    for i, text in enumerate(comment_texts):
        if sentinel_idx is not None and i >= sentinel_idx:
            continue
        if is_blas_boilerplate_text(text):
            boilerplate_indices.append(i)
            skip_flags[i] = True

    if not boilerplate_indices:
        return set()

    # 3) Group boilerplate indices into blocks (distance <= 10 lines).
    blocks = []
    current_block = [boilerplate_indices[0]]
    for idx in boilerplate_indices[1:]:
        if idx - current_block[-1] <= 10:
            current_block.append(idx)
        else:
            blocks.append((min(current_block), max(current_block)))
            current_block = [idx]
    blocks.append((min(current_block), max(current_block)))

    # 4) For each block, mark everything inside, plus leading/trailing padding.
    to_skip = set()

    for start, end in blocks:
        # Mark block itself.
        for i in range(start, end + 1):
            to_skip.add(i)

        # Preceding padding comments.
        j = start - 1
        while j >= 0:
            if _is_empty_or_dotdot(comment_texts[j]):
                to_skip.add(j)
                j -= 1
            else:
                break

        # Following padding comments.
        j = end + 1
        while j < n:
            if _is_empty_or_dotdot(comment_texts[j]):
                to_skip.add(j)
                j += 1
            else:
                break

    return to_skip


# ---------------------------------------------------------------------------
#  File processing
# ---------------------------------------------------------------------------

def strip_boilerplate_comments(path: Path) -> None:
    """Remove BLAS/LAPACK boilerplate comment *blocks* from a file in-place."""
    try:
        src_lines = path.read_text(
            encoding="utf-8", errors="ignore"
        ).splitlines(keepends=True)
    except OSError as e:
        print(f"strip_boilerplate_comments: failed to read {path}: {e}", file=sys.stderr)
        return

    # Collect comment lines and their texts (in order of appearance).
    comment_indices = []   # indices into src_lines
    comment_texts = []     # extracted comment text (without //, without C)

    for idx, line in enumerate(src_lines):
        if not is_comment_line(line):
            continue
        text = extract_comment_text(line)
        comment_indices.append(idx)
        comment_texts.append(text)

    # Decide which comment *indices* (in comment_texts) are boilerplate.
    boilerplate_comment_indices = _compute_boilerplate_comment_indices(comment_texts)

    # Map them back to original src line indices.
    skip_src_indices = set()
    for local_idx in boilerplate_comment_indices:
        src_idx = comment_indices[local_idx]
        skip_src_indices.add(src_idx)

    # Write back file without skipped lines.
    out_lines = []
    for idx, line in enumerate(src_lines):
        if idx in skip_src_indices:
            continue
        out_lines.append(line)

    try:
        path.write_text("".join(out_lines), encoding="utf-8")
    except OSError as e:
        print(f"strip_boilerplate_comments: failed to write {path}: {e}", file=sys.stderr)


# ---------------------------------------------------------------------------
#  CLI
# ---------------------------------------------------------------------------

def main(argv):
    if len(argv) < 2:
        print("Usage: strip_boilerplate_comments.py FILE [FILE...]", file=sys.stderr)
        sys.exit(1)

    for arg in argv[1:]:
        p = Path(arg)
        if not p.is_file():
            print(f"strip_boilerplate_comments: {p} is not a file, skipping", file=sys.stderr)
            continue
        strip_boilerplate_comments(p)


if __name__ == "__main__":
    main(sys.argv)
