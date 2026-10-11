#!/usr/bin/env python3
"""Rename one backend to another in a file, as the first step of a copy.

    misc/backend_rename.py <model> <new> <file>      # writes to stdout

<model> and <new> are backend names from backends.txt (for example dd and
td).  The model's name is replaced where it is a token of an identifier,
macro or path: dd_real, Rlamch_dd, _dd_opt, MPLAPACK_BUILD_WITH_DD, DD_LABEL,
Rgemm.dd, mpblas/optimized/dd.  It is not replaced inside other words, so
LDADD, top_builddir and "add" stay intact.

Names that all libQD3 backends share are never renamed (KEEP): the libQD3
install variables and library, the qd/ header directory, and DD_CXXFLAGS,
which holds the no-FMA flags used by every backend with the nofma trait.

See docs/adding-a-backend.md, sections 1 and 2.2.
"""

import re
import sys

KEEP = [
    r"QD_LIBDIR", r"QD_INCLUDEDIR", r"QD_RPATH_LDFLAGS", r"QD_EXTRA_CXXFLAGS",
    r"qd_extra_cxxflags", r"DD_CXXFLAGS", r"-lqd\b", r"\blibqd\b",
    r"(?<=[<\"])qd/",   # #include <qd/...>
    r"/include/qd\b",   # -I$(prefix)/include/qd
]


def rename(text, model, new):
    saved = []

    def protect(m):
        saved.append(m.group(0))
        return f"\0{len(saved) - 1}\0"

    for pattern in KEEP:
        text = re.sub(pattern, protect, text)
    for old, rep in ((model, new), (model.upper(), new.upper())):
        # A token boundary is anything but a letter or digit; "_" and "."
        # separate tokens (dd_real, Rlamch_dd, Rgemm.dd).
        text = re.sub(rf"(?<![A-Za-z0-9]){old}(?![A-Za-z0-9])", rep, text)
    return re.sub(r"\0(\d+)\0", lambda m: saved[int(m.group(1))], text)


def main():
    if len(sys.argv) != 4:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    model, new, path = sys.argv[1:]
    with open(path) as f:
        sys.stdout.write(rename(f.read(), model, new))
    return 0


if __name__ == "__main__":
    sys.exit(main())
