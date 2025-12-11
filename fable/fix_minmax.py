#!/usr/bin/env python3
import re
import sys
from pathlib import Path

# Convert max({a, b, c}) / min({a, b, c})  --> max(a, b, c) / min(a, b, c)
# Works also with spaces and newlines inside the braces.

pattern = re.compile(r'\b(max|min)\s*\(\s*\{([^{}]*?)\}\s*\)', re.DOTALL)


def fix_text(text: str) -> str:
    def repl(m: re.Match) -> str:
        func = m.group(1)      # "max" or "min"
        args = m.group(2)      # "a, b, c"
        return f"{func}({args})"
    return pattern.sub(repl, text)


def main() -> None:
    # No arguments: read from stdin, write to stdout
    if len(sys.argv) == 1:
        sys.stdout.write(fix_text(sys.stdin.read()))
        return

    # With filenames: rewrite files in-place
    for path_str in sys.argv[1:]:
        path = Path(path_str)
        src = path.read_text(encoding="utf-8")
        fixed = fix_text(src)
        path.write_text(fixed, encoding="utf-8")


if __name__ == "__main__":
    main()
