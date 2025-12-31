#!/usr/bin/env python3
"""
Stage comment-only changes from working tree into the Git index, leaving code changes unstaged.

Notes:
- Targets C/C++-like files (.c/.cc/.cpp/.h/.hpp/.cxx/.hxx etc.).
- "Comment-only line" means: line is blank/whitespace OR begins with //, /*, *, */, or is inside a /* */ block.
- Inline comments on code lines (e.g. 'int x; // comment') are treated as CODE lines and will NOT be staged.
- This is heuristic, not a full parser.
"""

from __future__ import annotations

import argparse
import difflib
import os
import subprocess
import sys
from typing import List, Optional, Tuple


def run_git(args: List[str], input_bytes: Optional[bytes] = None, check: bool = True) -> Tuple[int, bytes, bytes]:
    p = subprocess.run(
        ["git"] + args,
        input=input_bytes,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if check and p.returncode != 0:
        sys.stderr.write(p.stderr.decode("utf-8", errors="replace"))
    return p.returncode, p.stdout, p.stderr


def read_worktree(path: str) -> bytes:
    with open(path, "rb") as f:
        return f.read()


def read_head_blob(path: str) -> bytes:
    rc, out, _ = run_git(["show", f"HEAD:{path}"], check=False)
    if rc != 0:
        raise SystemExit(f"ERROR: Cannot read HEAD:{path} (is the file tracked?)")
    return out


def decode_utf8(b: bytes, path: str) -> str:
    try:
        return b.decode("utf-8")
    except UnicodeDecodeError:
        raise SystemExit(f"ERROR: Non-UTF8 file refused: {path}")


def is_c_like(path: str) -> bool:
    ext = os.path.splitext(path)[1].lower()
    return ext in {
        ".c", ".h", ".cc", ".cpp", ".cxx", ".hh", ".hpp", ".hxx",
        ".ipp", ".tpp", ".inc", ".java", ".js", ".ts", ".tsx", ".jsx",
    }


def classify_comment_only_cpp(lines: List[str]) -> List[bool]:
    """
    Return a boolean list: True if the line is comment-only (or blank), False otherwise.
    Handles /* */ blocks when they start at the beginning of a line (after whitespace).
    """
    flags: List[bool] = []
    in_block = False

    for line in lines:
        stripped = line.lstrip()

        # Treat blank/whitespace-only lines as "non-code" (allowed to stage with comment changes).
        if stripped.strip() == "":
            flags.append(True)
            continue

        if in_block:
            flags.append(True)
            if "*/" in stripped:
                in_block = False
            continue

        if stripped.startswith("//"):
            flags.append(True)
            continue

        if stripped.startswith("/*"):
            flags.append(True)
            # If it doesn't end on the same line, enter block mode.
            if "*/" not in stripped[2:]:
                in_block = True
            continue

        if stripped.startswith("*") or stripped.startswith("*/"):
            flags.append(True)
            continue

        # Otherwise, treat as code (even if it contains inline comments).
        flags.append(False)

    return flags


def build_comment_only_version(base_lines: List[str], work_lines: List[str], base_is_comment: List[bool], work_is_comment: List[bool]) -> List[str]:
    """
    Create a synthetic file version:
    - Keep BASE content for code lines.
    - For hunks that touch only comment-only lines, take WORK content.
    - For replace blocks with equal lengths, selectively take WORK line-by-line for comment-only lines.
    """
    sm = difflib.SequenceMatcher(a=base_lines, b=work_lines)
    out: List[str] = []

    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        a = base_lines[i1:i2]
        b = work_lines[j1:j2]
        a_flags = base_is_comment[i1:i2]
        b_flags = work_is_comment[j1:j2]

        if tag == "equal":
            out.extend(a)
            continue

        if tag == "delete":
            # Only delete if all deleted lines are comment-only.
            if all(a_flags):
                continue
            out.extend(a)
            continue

        if tag == "insert":
            # Only insert if all inserted lines are comment-only.
            if all(b_flags):
                out.extend(b)
            # else ignore insertion
            continue

        if tag == "replace":
            # If lengths match, do per-line selection to salvage comment-only edits even when mixed.
            if len(a) == len(b):
                for k in range(len(a)):
                    if a_flags[k] and b_flags[k]:
                        out.append(b[k])
                    else:
                        out.append(a[k])
                continue

            # If entire block is comment-only on both sides, take WORK, else keep BASE.
            if all(a_flags) and all(b_flags):
                out.extend(b)
            else:
                out.extend(a)
            continue

    return out


def make_patch(path: str, base_lines: List[str], new_lines: List[str]) -> bytes:
    if base_lines == new_lines:
        return b""

    # difflib wants lines with line endings preserved
    diff_lines = list(difflib.unified_diff(
        base_lines,
        new_lines,
        fromfile=f"a/{path}",
        tofile=f"b/{path}",
        n=3,
        lineterm="\n",
    ))

    # Prepend a minimal git-style header to keep git apply happy.
    header = [f"diff --git a/{path} b/{path}\n"]
    return ("".join(header + diff_lines)).encode("utf-8")


def stage_patch(patch: bytes) -> None:
    if not patch:
        return
    rc, _, err = run_git(["apply", "--cached", "--whitespace=nowarn", "--recount", "-"], input_bytes=patch, check=False)
    if rc != 0:
        sys.stderr.write(err.decode("utf-8", errors="replace"))
        raise SystemExit("ERROR: git apply --cached failed (patch did not apply).")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="+", help="Files to stage comment-only changes for (working tree -> index).")
    ap.add_argument("--dry-run", action="store_true", help="Do not stage; just report.")
    args = ap.parse_args()

    # Refuse if there are already staged changes (avoid accidental mixing).
    rc, _, _ = run_git(["diff", "--cached", "--quiet"], check=False)
    if rc != 0:
        sys.stderr.write("ERROR: You already have staged changes. Commit or unstage them first.\n")
        return 2

    staged_any = False

    for path in args.paths:
        if not os.path.exists(path):
            sys.stderr.write(f"WARN: file not found in worktree, skipping: {path}\n")
            continue
        if not is_c_like(path):
            sys.stderr.write(f"WARN: not a C-like file, skipping: {path}\n")
            continue

        base_txt = decode_utf8(read_head_blob(path), path)
        work_txt = decode_utf8(read_worktree(path), path)

        base_lines = base_txt.splitlines(keepends=True)
        work_lines = work_txt.splitlines(keepends=True)

        base_flags = classify_comment_only_cpp(base_lines)
        work_flags = classify_comment_only_cpp(work_lines)

        comment_only_lines = build_comment_only_version(base_lines, work_lines, base_flags, work_flags)
        patch = make_patch(path, base_lines, comment_only_lines)

        if not patch:
            print(f"[SKIP] {path}: no comment-only changes to stage")
            continue

        if args.dry_run:
            print(f"[DRY]  {path}: would stage comment-only changes")
            continue

        stage_patch(patch)
        staged_any = True
        print(f"[OK]   {path}: staged comment-only changes")

    if staged_any and not args.dry_run:
        print("\nStaged diff (cached):")
        subprocess.run(["git", "diff", "--cached", "--name-only"])

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
