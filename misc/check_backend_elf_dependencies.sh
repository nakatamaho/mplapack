#!/bin/sh
# Check the ELF dependency closure of an MPLAPACK backend shared library.
# usage: check_backend_elf_dependencies.sh <library> <dependency-stem> ...
set -eu

if test "$#" -lt 2; then
    echo "usage: $0 <library> <dependency-stem> ..." >&2
    exit 2
fi

library=$1
shift
if ! command -v readelf >/dev/null 2>&1 ||
   ! command -v ldd >/dev/null 2>&1; then
    echo "SKIP: readelf and ldd are required for ELF dependency QA"
    exit 77
fi

needed=$(readelf -d "$library")
status=0
for dependency in "$@"; do
    if ! printf '%s\n' "$needed" | grep -Eq "Shared library: \[lib${dependency}\.so"; then
        echo "FAIL: $library has no DT_NEEDED entry for lib${dependency}.so" >&2
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
