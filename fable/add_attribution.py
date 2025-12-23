#!/usr/bin/env python3
"""Add or normalize BLAS/LAPACK attribution header comments in MPLAPACK sources.

This script can either:

  * output the normalized content to stdout ("stdout mode"), or
  * edit files in-place ("in-place mode").

In-place mode will *not* rewrite a file when the content is unchanged.

Usage:
  add_attribution.py [--inplace|-i] CPP_FILE FORTRAN_FILE [CPP_FILE FORTRAN_FILE ...]
  add_attribution.py --stdout CPP_FILE FORTRAN_FILE

Notes:
  * If neither --inplace nor --stdout is given, the default is in-place mode
    (backward compatible with older versions of this script).
  * In stdout mode, exactly one CPP/Fortran pair is supported.

For each C++ file, this script will either:

  * replace an existing BLAS/LAPACK attribution block (lines containing
    'Derived from ... routine', 'Original ... authors', or known author
    markers), or

  * insert a new canonical attribution block after the leading license
    comment (if no attribution is present).

The canonical form is:

  // Derived from <LIB> routine <ROUTINE>.
  // Original <LIB> authors:
  //   <author1>
  //   <author2>
  //   ...

where <LIB> is "BLAS" or "LAPACK" depending on the Fortran file path, and
<ROUTINE> is taken from the Fortran file stem (e.g. zsymm.f -> ZSYMM).

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

from __future__ import annotations

import sys
from pathlib import Path
from typing import List, Tuple


# Markers that indicate we are inside an old attribution block.
ATTR_MARKERS = [
    "Derived from LAPACK routine",
    "Derived from BLAS routine",
    "Original LAPACK authors",
    "Original BLAS authors",
    "LAPACK is a software package provided",
    "Univ. of Tennessee",
    "Univ. of California Berkeley",
    "Univ. of Colorado Denver",
    "NAG Ltd",
]


def _print_usage(out) -> None:
    prog = Path(sys.argv[0]).name
    print(
        "Usage:\n"
        f"  {prog} [--inplace|-i] CPP_FILE FORTRAN_FILE [CPP_FILE FORTRAN_FILE ...]\n"
        f"  {prog} --stdout CPP_FILE FORTRAN_FILE\n\n"
        "Modes:\n"
        "  --inplace, -i   Edit file(s) in-place (default).\n"
        "  --stdout        Print normalized content to stdout.\n",
        file=out,
    )


def is_comment_line(line: str) -> bool:
    """Return True if the line is a C++ line comment starting with '//'"""
    stripped = line.lstrip()
    return stripped.startswith("//")


def line_contains_attr_marker(line: str) -> bool:
    """Return True if this line looks like part of an attribution block."""
    stripped = line.lstrip()
    if not stripped.startswith("//"):
        return False
    text = stripped[2:].strip()
    for kw in ATTR_MARKERS:
        if kw in text:
            return True
    return False


def extract_authors_from_fortran(path: Path) -> List[str]:
    """Extract author lines from a Fortran file."""
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


def detect_library_name(fortran_path: Path) -> str:
    """Detect whether the routine comes from BLAS or LAPACK based on the path.

    If the resolved path contains 'blas' (case-insensitive), we treat it as BLAS.
    Otherwise we fall back to LAPACK.
    """
    lib = "LAPACK"
    try:
        path_str = str(fortran_path.resolve()).lower()
    except Exception:
        path_str = str(fortran_path).lower()

    if "blas" in path_str:
        lib = "BLAS"

    return lib


def build_canonical_header(
    routine_name: str, lib_name: str, authors: List[str], indent: str = ""
) -> List[str]:
    """Build the canonical attribution header block."""
    lines: List[str] = []
    lines.append(f"{indent}// Derived from {lib_name} routine {routine_name}.")
    lines.append(f"{indent}// Original {lib_name} authors:")

    if not authors:
        # Fallback generic authors block if extraction failed.
        lines.append(f"{indent}//   Univ. of Tennessee,")
        lines.append(f"{indent}//   Univ. of California Berkeley,")
        lines.append(f"{indent}//   Univ. of Colorado Denver,")
        lines.append(f"{indent}//   NAG Ltd.")
    else:
        for a in authors:
            lines.append(f"{indent}//   {a}")

    # Add a blank line (no '//') after the header.
    lines.append("")
    return [ln + "\n" for ln in lines]


def insert_header_if_missing(
    src: List[str], routine: str, lib_name: str, authors: List[str]
) -> List[str]:
    """Insert canonical header if no existing attribution block is found.

    Insert after the leading license block (/* ... */) if present,
    otherwise at the top of the file after initial blank lines.
    """
    n = len(src)
    i = 0

    # Skip initial blank lines
    while i < n and src[i].strip() == "":
        i += 1

    # If there is a leading C-style comment block (license), skip it
    if i < n and src[i].lstrip().startswith("/*"):
        while i < n and "*/" not in src[i]:
            i += 1
        if i < n:
            i += 1  # skip the line with */
        # Skip blank lines after license block
        while i < n and src[i].strip() == "":
            i += 1

    insert_at = i
    header_lines = build_canonical_header(
        routine, lib_name, authors, indent="")
    return src[:insert_at] + header_lines + src[insert_at:]


def replace_existing_header(
    src: List[str], routine: str, lib_name: str, authors: List[str]
) -> List[str]:
    """Replace an existing attribution block with the canonical header."""
    out: List[str] = []
    n = len(src)
    i = 0
    replaced_block = False

    while i < n:
        line = src[i]

        if not is_comment_line(line) or not line_contains_attr_marker(line):
            out.append(line)
            i += 1
            continue

        # Determine indentation from this line.
        leading = line[: len(line) - len(line.lstrip())]

        # Skip all contiguous comment lines that look like attribution
        # or are empty comment lines immediately around them.
        while i < n and is_comment_line(src[i]):
            text = src[i].lstrip()[2:].strip()
            if text == "" or line_contains_attr_marker(src[i]):
                i += 1
            else:
                break

        # Drop any real blank lines immediately following the attribution
        # block so the transformation is idempotent (the canonical header
        # already ends with a single blank line).
        while i < n and src[i].strip() == "":
            i += 1

        if not replaced_block:
            out.extend(build_canonical_header(
                routine, lib_name, authors, indent=leading))
            replaced_block = True
        # Old block is dropped entirely.

    return out if replaced_block else src


def normalize_text(cpp_path: Path, fortran_path: Path, src_text: str) -> str:
    """Return normalized C++ content as text."""
    src_lines = src_text.splitlines(keepends=True)

    # Routine name: prefer the Fortran stem, uppercased (e.g. zsymm.f -> ZSYMM).
    if fortran_path.is_file():
        routine = fortran_path.stem.upper()
    else:
        routine = cpp_path.stem.upper()

    # Library name: BLAS if the Fortran path contains 'blas', otherwise LAPACK.
    lib_name = detect_library_name(fortran_path)

    # Extract authors from the given Fortran file (if it exists).
    authors: List[str] = []
    if fortran_path.is_file():
        authors = extract_authors_from_fortran(fortran_path)

    # Check if there is any existing attribution block.
    has_attr = any(
        is_comment_line(line) and line_contains_attr_marker(line) for line in src_lines
    )

    if has_attr:
        new_src = replace_existing_header(
            src_lines, routine, lib_name, authors)
    else:
        new_src = insert_header_if_missing(
            src_lines, routine, lib_name, authors)

    return "".join(new_src)


def normalize_file(cpp_path: Path, fortran_path: Path) -> Tuple[str, bool]:
    """Normalize attribution header in a single C++ file.

    Returns (output_text, changed).
    """
    src_text = cpp_path.read_text(encoding="utf-8", errors="ignore")
    out_text = normalize_text(cpp_path, fortran_path, src_text)
    return out_text, (out_text != src_text)


def main(argv: List[str]) -> int:
    inplace = False
    stdout_mode = False
    args: List[str] = []

    it = iter(argv[1:])
    for a in it:
        if a in ("-h", "--help"):
            _print_usage(sys.stdout)
            return 0
        if a in ("--inplace", "-i"):
            inplace = True
            continue
        if a == "--stdout":
            stdout_mode = True
            continue
        if a == "--":
            args.extend(list(it))
            break
        args.append(a)

    if inplace and stdout_mode:
        print("add_attribution: choose either --inplace or --stdout", file=sys.stderr)
        return 2

    # Backward compatible default: in-place if no explicit mode is given.
    if not inplace and not stdout_mode:
        inplace = True

    if stdout_mode:
        if len(args) != 2:
            _print_usage(sys.stderr)
            return 1
        pairs = [(args[0], args[1])]
    else:
        if len(args) < 2 or (len(args) % 2) != 0:
            _print_usage(sys.stderr)
            return 1
        pairs = [(args[i], args[i + 1]) for i in range(0, len(args), 2)]

    for idx, (cpp_arg, ft_arg) in enumerate(pairs):
        cpp_path = Path(cpp_arg)
        ft_path = Path(ft_arg)

        if not cpp_path.is_file():
            print(
                f"add_attribution: {cpp_path} is not a file", file=sys.stderr)
            return 3
        if not ft_path.is_file():
            print(
                f"add_attribution: WARNING: {ft_path} is not a file, authors will be empty",
                file=sys.stderr,
            )

        out_text, changed = normalize_file(cpp_path, ft_path)

        if stdout_mode:
            sys.stdout.write(out_text)
            return 0

        # In-place mode
        if changed:
            cpp_path.write_text(out_text, encoding="utf-8")

        # If unchanged: do not rewrite.
        # Continue to next pair.
        _ = idx

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
