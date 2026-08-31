#!/usr/bin/env python3
"""Load an MPFR backend shared library in a clean process."""

import ctypes
import os
import sys


if len(sys.argv) != 2:
    print(f"usage: {sys.argv[0]} <MPFR shared library>", file=sys.stderr)
    sys.exit(2)

library = os.path.abspath(sys.argv[1])
try:
    ctypes.CDLL(library, mode=getattr(os, "RTLD_LOCAL", 0))
except OSError as error:
    print(f"FAIL: unable to load {library}: {error}", file=sys.stderr)
    sys.exit(1)

print(f"PASS: clean-process load of {library}")
