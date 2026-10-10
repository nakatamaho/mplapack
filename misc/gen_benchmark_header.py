#!/usr/bin/env python3
"""Write include/mplapack_benchmark_<name>.h from a built backend library.

The header maps SYMBOL_GCC_<ROUTINE> to the mangled name of each routine in
libmplapack_<name>.  Mangled names cannot be derived by renaming another
backend's header: libQD3 complex types are aliases of qd3_complex<T>, which
changes both the type encoding and the substitution indices.  This script
takes the routine list from an existing backend header and looks up every
routine in the library's symbol table.

    misc/gen_benchmark_header.py <name> <libmplapack_name.a> [--like dd]

Routines that the library does not define are dropped and reported.
"""

import argparse
import pathlib
import re
import subprocess
import sys

TOP = pathlib.Path(__file__).resolve().parent.parent


def demangled_name(symbol):
    m = re.match(r"_Z(\d+)", symbol)
    if not m:
        return None
    n = int(m.group(1))
    start = len(m.group(0))
    return symbol[start:start + n]


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("name", help="backend name, e.g. td")
    parser.add_argument("library", help="path to libmplapack_<name>.a or .so")
    parser.add_argument("--like", default="dd",
                        help="backend whose header supplies the routine list (default: dd)")
    args = parser.parse_args()

    nm = subprocess.run(["nm", "--defined-only", args.library],
                        capture_output=True, text=True, check=True).stdout
    by_name = {}
    for line in nm.splitlines():
        fields = line.split()
        if len(fields) == 3 and fields[1] in ("T", "W") and fields[2].startswith("_Z"):
            by_name.setdefault(demangled_name(fields[2]), set()).add(fields[2])

    like = args.like
    source = (TOP / "include" / f"mplapack_benchmark_{like}.h").read_text().split("\n")
    suffix_re = re.compile(rf"_({re.escape(like)}|qd)$")
    macro_re = re.compile(rf"_({re.escape(like.upper())}|QD)$")
    out, missing = [], []
    for line in source:
        m = re.match(r'#define (SYMBOL_GCC_\w+) "(_Z\d+\w+?)"$', line)
        if not m:
            out.append(line)
            continue
        routine = suffix_re.sub(f"_{args.name}", demangled_name(m.group(2)))
        macro = macro_re.sub(f"_{args.name.upper()}", m.group(1))
        candidates = sorted(by_name.get(routine, ()))
        if len(candidates) != 1:
            missing.append((routine, candidates))
            continue
        out.append(f'#define {macro} "{candidates[0]}"')

    path = TOP / "include" / f"mplapack_benchmark_{args.name}.h"
    path.write_text("\n".join(out))
    print(f"wrote {path.relative_to(TOP)}")
    for routine, candidates in missing:
        what = "not defined" if not candidates else f"ambiguous: {candidates}"
        print(f"  dropped {routine} ({what})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
