#!/bin/bash
set -eu

# The per-backend Makefile.am files are generated from backends.txt.
if command -v python3 >/dev/null 2>&1; then
    python3 misc/gen_backend_makefiles.py --check
else
    echo "gen_configure.sh: python3 not found; not checking generated Makefile.am files" >&2
fi

aclocal
autoheader
automake -a -v --add-missing
autoconf
