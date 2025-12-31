#!/usr/bin/env python3
"""
Stage comment-only changes from working tree into the Git index, leaving code changes unstaged.

Key behavior:
- Compares HEAD vs working tree for each given file.
- Builds a synthetic "comment-only" version:
  - Code lines are kept from HEAD.
  - Comment-only lines (or blank lines) are taken from the working tree.
- Applies the resulting patch to the index only (git apply --cached).
- Code changes remain unstaged in the working tree.

Notes / limitations:
- Targets C/C++-like files by extension.
- "Comment-only line" means:
    * blank/whitespace-only line, OR
    * line that begins (after whitespace) with //, /*, *, or */, OR
    * line inside a /* ... */ block (when the block started at line-begin).
- Inline comments on code lines (e.g. 'int x; // comment') are treated as CODE lines and will NOT be staged.
- This is a heuristic, not a full parser.
- Requires files to be tracked in HEAD (reads HEAD:<path>).
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


def git_text(args: List[str]) -> str:
    rc, out, err = run_git(args, check=False)
    if rc != 0:
        sys.stderr.write(err.decode("utf-8", errors="replace"))
        raise SystemExit(f"ERROR: git {' '.join(args)} failed")
    return out.decode("utf-8", errors="replace").strip()


def resolve_input_path(arg: str, repo_root: str) -> Tuple[str, str]:
    """
    Return (repo_rel_path, fs_abs_path).
    Accepts:
      - path relative to current directory (e.g. Clagsy.cpp from matgen dir)
      - path relative to repo root (e.g. mplapack/test/matgen/Clagsy.cpp)
      - absolute path
    """
    repo_root_abs = os.path.abspath(repo_root)

    # Candidate A: treat as path from current working directory
    abs_a = os.path.abspath(arg)
    if os.path.exists(abs_a):
        if os.path.commonpath([repo_root_abs, abs_a]) != repo_root_abs:
            raise SystemExit(f"ERROR: Path is outside repo: {arg}")
        repo_rel = os.path.relpath(abs_a, repo_root_abs).replace(os.sep, "/")
        return repo_rel, abs_a

    # Candidate B: treat as repo-root-relative path
    abs_b = os.path.abspath(os.path.join(repo_root_abs, arg))
    if os.path.exists(abs_b):
        repo_rel = os.path.relpath(abs_b, repo_root_abs).replace(os.sep, "/")
        return repo_rel, abs_b

    raise SystemExit(f"ERROR: File not found: {arg}")


def read_worktree(abs_path: str) -> bytes:
    with open(abs_path, "rb") as f:
        return f.read()


def read_head_blob(repo_rel_path: str) -> bytes:
    rc, out, _ = run_git(["show", f"HEAD:{repo_rel_path}"], check=False)
    if rc != 0:
        raise SystemExit(f"ERROR: Cannot read HEAD:{repo_rel_path} (is the file tracked?)")
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
        ".ipp", ".tpp", ".inc",
        ".java", ".js", ".jsx", ".ts", ".tsx",
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

        # Blank/whitespace-only lines are considered "non-code" and are safe to stage with comments.
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
            # If it does not end on the same line, enter block mode.
            if "*/" not in stripped[2:]:
                in_block = True
            continue

        if stripped.startswith("*") or stripped.startswith("*/"):
            flags.append(True)
            continue

        # Otherwise treat as code (even if it contains inline comments).
        flags.append(False)

    return flags


def build_comment_only_version(
    base_lines: List[str],
    work_lines: List[str],
    base_is_comment: List[bool],
    work_is_comment: List[bool],
) -> List[str]:
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
            # else ignore insertion (keep base)
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


def make_patch(repo_rel_path: str, base_lines: List[str], new_lines: List[str]) -> bytes:
    if base_lines == new_lines:
        return b""

    # difflib wants lines with line endings preserved
    diff_lines = list(
        difflib.unified_diff(
            base_lines,
            new_lines,
            fromfile=f"a/{repo_rel_path}",
            tofile=f"b/{repo_rel_path}",
            n=3,
            lineterm="\n",
        )
    )

    # Prepend a minimal git-style header to keep git apply happy.
    header = [f"diff --git a/{repo_rel_path} b/{repo_rel_path}\n"]
    return ("".join(header + diff_lines)).encode("utf-8")


def stage_patch(patch: bytes) -> None:
    if not patch:
        return
    rc, _, err = run_git(
        ["apply", "--cached", "--whitespace=nowarn", "--recount", "-"],
        input_bytes=patch,
        check=False,
    )
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

    repo_root = git_text(["rev-parse", "--show-toplevel"])
    staged_any = False

    for arg in args.paths:
        try:
            repo_path, fs_path = resolve_input_path(arg, repo_root)
        except SystemExit as e:
            sys.stderr.write(str(e) + "\n")
            continue

        if not is_c_like(repo_path):
            sys.stderr.write(f"WARN: not a C-like file, skipping: {repo_path}\n")
            continue

        base_txt = decode_utf8(read_head_blob(repo_path), repo_path)
        work_txt = decode_utf8(read_worktree(fs_path), fs_path)

        base_lines = base_txt.splitlines(keepends=True)
        work_lines = work_txt.splitlines(keepends=True)

        base_flags = classify_comment_only_cpp(base_lines)
        work_flags = classify_comment_only_cpp(work_lines)

        comment_only_lines = build_comment_only_version(base_lines, work_lines, base_flags, work_flags)
        patch = make_patch(repo_path, base_lines, comment_only_lines)

        if not patch:
            print(f"[SKIP] {repo_path}: no comment-only changes to stage")
            continue

        if args.dry_run:
            print(f"[DRY]  {repo_path}: would stage comment-only changes")
            continue

        stage_patch(patch)
        staged_any = True
        print(f"[OK]   {repo_path}: staged comment-only changes")

    if staged_any and not args.dry_run:
        print("\nStaged files (cached):")
        subprocess.run(["git", "diff", "--cached", "--name-only"], check=False)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
