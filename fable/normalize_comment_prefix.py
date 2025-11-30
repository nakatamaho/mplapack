#!/usr/bin/env python3
"""
Normalize LAPACK attribution header comments in MPLAPACK sources.

For each C++/header file, if there is an existing LAPACK attribution block
(contains 'Derived from LAPACK routine', 'Original LAPACK authors',
or 'LAPACK is a software package provided'), replace that whole comment
block with a canonical form:

// Derived from LAPACK routine <ROUTINE>.
// Original LAPACK authors:
//   <author1>
//   <author2>
//   ...

Authors are extracted from the corresponding Fortran source file in the
same directory, by parsing the "Authors" section:

*  Authors:
*  ========
*
*> \\author Univ. of Tennessee
*> \\author Univ. of California Berkeley
*> \\author Univ. of Colorado Denver
*> \\author NAG Ltd.

The Fortran source is expected to be in the same directory as the C++ file,
with a matching stem name, e.g.:

Cbdsqr.cpp  <->  cbdsqr.f   or cbdsqr.f90
"""

import sys
from pathlib import Path
from typing import List, Optional

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


def find_fortran_for_cpp(cpp_path: Path) -> Optional[Path]:
    """Locate the Fortran source corresponding to a given C++ file.

    Search patterns in the same directory as cpp_path:
        <stem>.f
        <stem>.f90
        <stem>.F
        <stem>.F90
    """
    stem = cpp_path.stem.lower()
    directory = cpp_path.parent

    for ext in (".f", ".f90", ".F", ".F90"):
        candidate = directory / (stem + ext)
        if candidate.is_file():
            return candidate

    return None


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
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    authors: List[str] = []
    in_authors_section = False

    for line in lines:
        stripped = line.lstrip()

        # Detect start of Authors section
        if not in_authors_section:
            if "Authors:" in stripped:
                in_authors_section = True
            continue

        # Once in Authors section, collect \author lines.
        if "\\author" in stripped:
            # Typical form: "*> \\author Univ. of Tennessee"
            after = stripped.split("\\author", 1)[1]
            author = after.strip()
            if author:
                authors.append(author)

        # Heuristic: stop when we leave the comment block (non-* line).
        if stripped and not stripped.startswith("*"):
            break

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


def normalize_file(path: Path) -> None:
    """Normalize LAPACK attribution header in a single C++/header file."""
    src = path.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True)
    out: List[str] = []

    # Routine name from file stem, uppercased (e.g. Cbdsqr.cpp -> CBDSQR)
    routine = path.stem.upper()

    # Try to find and parse the corresponding Fortran source in the same directory.
    authors: List[str] = []
    fpath = find_fortran_for_cpp(path)
    if fpath is not None:
        authors = extract_authors_from_fortran(fpath)

    i = 0
    n = len(src)
    replaced_block = False

    while i < n:
        line = src[i]

        if not is_comment_line(line) or not line_contains_lapack_attr(line):
            out.append(line)
            i += 1
            continue

        # We found a line that is part of a LAPACK attribution block.
        # Determine indentation from this line.
        leading = line[: len(line) - len(line.lstrip())]

        # Skip all contiguous comment lines that look like LAPACK attr
        # or are empty comment lines immediately around them.
        while i < n and is_comment_line(src[i]):
            text = src[i].lstrip()[2:].strip()
            if text == "" or line_contains_lapack_attr(src[i]):
                i += 1
            else:
                # Stop at the first comment line that is not clearly part of the LAPACK attr.
                break

        if not replaced_block:
            # Insert canonical header once, at the position of the first removed line.
            out.extend(build_canonical_header(routine, authors, indent=leading))
            replaced_block = True
        # Old block is dropped entirely.

    path.write_text("".join(out), encoding="utf-8")


def main(argv):
    if len(argv) < 2:
        print("Usage: normalize_lapack_header.py FILE [FILE...]", file=sys.stderr)
        sys.exit(1)

    for arg in argv[1:]:
        p = Path(arg)
        if not p.is_file():
            print(f"normalize_lapack_header: {p} is not a file, skipping", file=sys.stderr)
            continue
        normalize_file(p)


if __name__ == "__main__":
    main(sys.argv)
