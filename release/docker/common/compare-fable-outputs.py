#!/usr/bin/env python3
"""Compare generated library/test sources, public headers and source lists."""

import difflib
from pathlib import Path
import sys


def outputs(root):
    paths = set()
    for directory in (
        "mpblas/reference", "mplapack/reference",
        "mplapack/test/eig/common", "mplapack/test/lin/common",
        "mplapack/test/matgen",
    ):
        paths.update(p.relative_to(root) for p in (root / directory).glob("*.cpp"))
    for directory in (
        "mpblas/reference", "mplapack/reference", "mpblas/optimized",
        "mplapack/test/eig", "mplapack/test/lin", "mplapack/test/matgen",
    ):
        paths.update(p.relative_to(root) for p in (root / directory).rglob("*.am"))
    for flavor in ("gmp", "mpfr", "qd", "dd", "double", "binary80", "binary128"):
        for prefix in ("mpblas", "mplapack", "mplapack_eig", "mplapack_lin", "mplapack_matgen"):
            path = Path("include") / f"{prefix}_{flavor}.h"
            if (root / path).is_file():
                paths.add(path)
    for name in (
        "include/mplapack_generic.h", "mpblas/reference/mpblas_generic.h",
        "mplapack/reference/mplapack_generic.h",
    ):
        if (root / name).is_file():
            paths.add(Path(name))
    for directory in ("mplapack/test/eig", "mplapack/test/lin"):
        for pattern in ("*.in", "*.rfp"):
            paths.update(
                p.relative_to(root) for p in (root / directory).glob(pattern)
                if p.name != "Makefile.in"
            )
    return paths


def main():
    baseline, generated = (Path(arg).resolve() for arg in sys.argv[1:])
    expected, actual = outputs(baseline), outputs(generated)
    if not expected or not actual:
        raise SystemExit("FAIL: no Fable outputs found")
    failures = 0
    for path in sorted(expected | actual):
        before, after = baseline / path, generated / path
        if path not in expected or path not in actual:
            print(f"FAIL: {'added' if path not in expected else 'removed'} {path}")
            failures += 1
        elif before.read_bytes() != after.read_bytes():
            print(f"FAIL: changed {path}")
            print("".join(difflib.unified_diff(
                before.read_text().splitlines(keepends=True),
                after.read_text().splitlines(keepends=True),
                fromfile=f"baseline/{path}", tofile=f"generated/{path}",
            )), end="")
            failures += 1
    print(f"Compared {len(expected | actual)} files; differences: {failures}")
    return bool(failures)


if __name__ == "__main__":
    sys.exit(main())
