#!/usr/bin/env python3
import sys
from pathlib import Path
import re

def normalize_comment_prefix(path: Path) -> None:
    src = path.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True)
    out = []
    for line in src:
        stripped = line.lstrip()
        if stripped.startswith("//"):
            indent_len = len(line) - len(stripped)
            indent = line[:indent_len]
            body = stripped[2:]  # after //
            # If there is any non-space text after //, normalize to one space
            if body.strip():
                body = " " + body.lstrip()
            line = indent + "//" + body
        out.append(line)
    path.write_text("".join(out), encoding="utf-8")

if __name__ == "__main__":
    for arg in sys.argv[1:]:
        p = Path(arg)
        if p.is_file():
            normalize_comment_prefix(p)
