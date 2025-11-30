#!/usr/bin/env python3
"""
Add or normalize LAPACK attribution header comments in MPLAPACK sources.

Usage:
    add_attribution.py CPP_FILE FORTRAN_FILE

For the given C++ file, if there is an existing LAPACK attribution block
(contains 'Derived from LAPACK routine', 'Original LAPACK authors',
or 'LAPACK is a software package provided'), replace that whole comment
block with a canonical form:

// Derived from LAPACK routine <ROUTINE>.
// Original LAPACK authors:
//   <author1>
//   <author2>
//   ...

Authors are extracted from the given Fortran source file, by parsing the
"Authors" section:

*  Authors:
*  ========
*
*> \\author Univ. of Tennessee
*> \\author Univ. of California Berkeley
*> \\author Univ. of Colorado Denver
*> \\author NAG Ltd.
"""

import sys
from pathlib import Path
from typing import List

# Markers that indicate we are inside an old LAPACK attribution block.
LAPACK_ATTR_MARKERS = [
    "Derived from LAPACK routine",
    "Original LAPACK authors",
    "LAPACK is a software package provided",
    "Univ. of Tennessee",
    "Univ. of California Berkeley",
    "Univ. of Colorado Denver",
    "NAG Ltd",
]


def debug(msg: str) -> None:
    """Print debug message to stderr."""
    print(f"[add_attribution] {msg}", file=sys.stderr)


def is_comment_line(line: str) -> bool:
    """Return True if the line is a C++ line comment starting with '//'."""
    stripped = line.lstrip()
    return stripped.startswith("//")


def line_contains_lapack_attr(line: str) -> bool:
    """Return True if this line looks like part of a LAPACK attribution block."""
    stripped = line.lstrip()
    if not stripped.startswith("//"):
        return False
    text = stripped[2:].strip()
    for kw in LAPACK_ATTR_MARKERS:
        if kw in text:
            return True
    return False


def extract_authors_from_fortran(path: Path) -> List[str]:
    """Extract author lines from a LAPACK Fortran file.

    Look for a section like:

    *  Authors:
    *  ========
    *
    *> \\author Univ. of Tennessee
    *> \\author Univ. of California Berkeley
    *> \\author Univ. of Colorado Denver
    *> \\author NAG Ltd.

    Return a list of author strings without the '\\author' prefix.
    """
    debug(f"Reading Fortran file for authors: {path}")
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    authors: List[str] = []
    in_authors_section = False

    for idx, line in enumerate(lines, start=1):
        stripped = line.lstrip()

        # Detect start of Authors section
        if not in_authors_section:
            if "Authors:" in stripped:
                debug(f"Found 'Authors:' section at Fortran line {idx}")
                in_authors_section = True
            continue

        # Once in Authors section, collect \author lines.
        if "\\author" in stripped:
            # Typical form: "*> \\author Univ. of Tennessee"
            after = stripped.split("\\author", 1)[1]
            author = after.strip()
            if author:
                debug(f"Found author at Fortran line {idx}: {author}")
                authors.append(author)

        # Heuristic: stop when we leave the comment block (non-* line).
        if stripped and not stripped.startswith("*"):
            debug(f"Leaving Authors section at Fortran line {idx}")
            break

    debug(f"Total authors extracted: {len(authors)}")
    return authors


def build_canonical_header(
    routine_name: str, authors: List[str], indent: str = ""
) -> List[str]:
    """Build the canonical LAPACK attribution header block."""
    lines: List[str] = []
    lines.append(f"{indent}// Derived from LAPACK routine {routine_name}.")
    lines.append(f"{indent}// Original LAPACK authors:")

    if not authors:
        # Fallback generic authors block if extraction failed.
        lines.append(f"{indent}//   Univ. of Tennessee,")
        lines.append(f"{indent}//   Univ. of California Berkeley,")
        lines.append(f"{indent}//   Univ. of Colorado Denver,")
        lines.append(f"{indent}//   NAG Ltd.")
    else:
        for a in authors:
            lines.append(f"{indent}//   {a}")

    lines.append(f"{indent}//")
    return [ln + "\n" for ln in lines]


def normalize_file(cpp_path: Path, fortran_path: Path) -> None:
    """Normalize LAPACK attribution header in a single C++ file."""
    debug(f"Normalizing attribution for C++ file: {cpp_path}")
    src = cpp_path.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True)
    out: List[str] = []

    # Routine name from C++ file stem, uppercased (e.g. Cbdsqr.cpp -> CBDSQR)
    routine = cpp_path.stem.upper()
    debug(f"Detected routine name from C++ stem: {routine}")

    # Extract authors from the given Fortran file (if it exists).
    authors: List[str] = []
    if fortran_path.is_file():
        authors = extract_authors_from_fortran(fortran_path)
    else:
        debug(f"Fortran file does not exist: {fortran_path}, authors list will be empty")

    i = 0
    n = len(src)
    replaced_block = False
    found_any_attr = False

    while i < n:
        line = src[i]

        if not is_comment_line(line) or not line_contains_lapack_attr(line):
            out.append(line)
            i += 1
            continue

        found_any_attr = True
        debug(f"Found LAPACK attribution marker at C++ line {i + 1}: {line.rstrip()}")

        # Determine indentation from this line.
        leading = line[: len(line) - len(line.lstrip())]

        # Skip all contiguous comment lines that look like LAPACK attr
        # or are empty comment lines immediately around them.
        start_i = i
        while i < n and is_comment_line(src[i]):
            text = src[i].lstrip()[2:].strip()
            if text == "" or line_contains_lapack_attr(src[i]):
                i += 1
            else:
                break
        debug(f"Skipping old LAPACK attribution block from C++ line {start_i + 1} to {i}")

        if not replaced_block:
            # Insert canonical header once, at the position of the first removed line.
            debug(f"Inserting canonical attribution header at C++ line {start_i + 1}")
            out.extend(build_canonical_header(routine, authors, indent=leading))
            replaced_block = True
        # Old block is dropped entirely.

    if not found_any_attr:
        debug("No LAPACK attribution block found in this C++ file")
    elif not replaced_block:
        debug("Attribution markers were found but no replacement was inserted (logic bug?)")
    else:
        debug("Attribution block successfully replaced")

    cpp_path.write_text("".join(out), encoding="utf-8")


def main(argv):
    if len(argv) != 3:
        print("Usage: add_attribution.py CPP_FILE FORTRAN_FILE", file=sys.stderr)
        sys.exit(1)

    cpp_path = Path(argv[1])
    ft_path = Path(argv[2])

    debug(f"CPP_FILE     = {cpp_path}")
    debug(f"FORTRAN_FILE = {ft_path}")

    if not cpp_path.is_file():
        print(f"add_attribution: {cpp_path} is not a file", file=sys.stderr)
        sys.exit(1)
    if not ft_path.is_file():
        print(f"add_attribution: WARNING: {ft_path} is not a file, authors will be empty", file=sys.stderr)

    normalize_file(cpp_path, ft_path)


if __name__ == "__main__":
    main(sys.argv)
