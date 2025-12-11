#!/usr/bin/env python3
import sys
import os
import re


def _split_code_comment(line: str) -> tuple[str, str]:
    """Split a C++ line into (code, comment) by '//'. Very simple heuristic."""
    idx = line.find("//")
    if idx < 0:
        return line, ""
    return line[:idx], line[idx:]


def _rewrite_c3_slice(code: str) -> str:
    """
    Rewrite 'c3(1, 1) == "G"' style comparisons to '*c3 == 'G''.
    """
    pattern = re.compile(
        r'c3\s*\(\s*1\s*,\s*1\s*\)\s*==\s*"([^"])"'
    )
    return pattern.sub(lambda m: f"*c3 == '{m.group(1)}'", code)


def _rewrite_c1_eq_char(code: str) -> str:
    """
    Rewrite 'c1 == "S"' to 'c1 == 'S'' for single-character comparisons.
    """
    pattern = re.compile(r'\bc1\s*==\s*"([^"])"')
    return pattern.sub(lambda m: f"c1 == '{m.group(1)}'", code)


def _rewrite_eq_to_strncmp(code: str, var: str, length: int) -> str:
    """
    Rewrite comparisons like 'c2 == "GE"' to 'strncmp(c2, "GE", 2) == 0'.
    This version allows any non-quote characters inside the literal,
    so strings like "QR " are also handled correctly.
    """
    pattern = re.compile(
        rf'\b{var}\s*==\s*"([^"]{{{length}}})"'
    )

    def repl(m: re.Match) -> str:
        lit = m.group(1)
        return f'strncmp({var}, "{lit}", {length}) == 0'

    return pattern.sub(repl, code)


def fix_text(text: str) -> str:
    """Apply all iMlaenv-specific fixes to the given source text."""
    # Safety guard: only touch files that clearly look like iMlaenv.cpp
    if "iMlaenv(" not in text:
        return text

    lines = text.splitlines(keepends=True)
    fixed_lines: list[str] = []

    for line in lines:
        code, comment = _split_code_comment(line)

        # 1) Fix c3(1, 1) == "G" style slicing.
        code = _rewrite_c3_slice(code)

        # 2) Fix c1 == "S" → c1 == 'S'.
        code = _rewrite_c1_eq_char(code)

        # 3) Rewrite c2/c3/c4 == "ABC" → strncmp(cX, "ABC", n) == 0.
        #    c2 and c4 are length-2; c3 is length-3 in iMlaenv.
        code = _rewrite_eq_to_strncmp(code, "c2", 2)
        code = _rewrite_eq_to_strncmp(code, "c4", 2)
        code = _rewrite_eq_to_strncmp(code, "c3", 3)

        fixed_lines.append(code + comment)

    return "".join(fixed_lines)


def process_file(path: str) -> None:
    """Read a file, apply fixes, and overwrite it in-place if changed."""
    with open(path, "r", encoding="utf-8") as f:
        original = f.read()

    fixed = fix_text(original)

    if fixed == original:
        return

    tmp_path = path + ".tmp"
    with open(tmp_path, "w", encoding="utf-8", newline="") as f:
        f.write(fixed)

    # Atomic replace on POSIX systems.
    os.replace(tmp_path, path)


def main() -> None:
    if len(sys.argv) < 2:
        print("Usage: fix_iMlaenv.py path/to/iMlaenv.cpp [...]",
              file=sys.stderr)
        sys.exit(1)

    for p in sys.argv[1:]:
        process_file(p)


if __name__ == "__main__":
    main()
