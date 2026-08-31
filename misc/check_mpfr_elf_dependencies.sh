#!/bin/sh
# Check the ELF dependency closure of an Autotools MPFR shared library.
# usage: check_mpfr_elf_dependencies.sh <MPFR shared library>
set -eu

if test "$#" -ne 1; then
    echo "usage: $0 <MPFR shared library>" >&2
    exit 2
fi

library=$1
if ! command -v readelf >/dev/null 2>&1 ||
   ! command -v ldd >/dev/null 2>&1; then
    echo "SKIP: readelf and ldd are required for ELF dependency QA"
    exit 77
fi

needed=$(readelf -d "$library")
status=0
for dependency in libmpc.so libmpfr.so libgmp.so; do
    if ! printf '%s\n' "$needed" | grep -Eq "Shared library: \[$dependency"; then
        echo "FAIL: $library has no DT_NEEDED entry for $dependency" >&2
        status=1
    fi
done

relocations=$(ldd -r "$library" 2>&1 || true)
if printf '%s\n' "$relocations" | grep -Eq 'undefined symbol|not found'; then
    echo "FAIL: unresolved runtime relocations for $library" >&2
    printf '%s\n' "$relocations" >&2
    status=1
fi

if test "$status" -eq 0; then
    echo "PASS: DT_NEEDED and runtime relocations for $library"
fi
exit "$status"
