#!/bin/sh
# Check backend/component relocations without executing numerical tests.
# FABLE common-block globals are supplied by the final test driver.
set -eu
if test "$#" -lt 1; then
    echo "usage: $0 <top-builddir> [dependency-library-dir ...]" >&2
    exit 2
fi
case $(uname -s) in
    Linux) ;;
    *) echo "SKIP: support relocation QA requires Linux"; exit 77 ;;
esac
command -v ldd >/dev/null 2>&1 || exit 77
builddir=$(cd "$1" && pwd)
shift
runtime_dirs="$builddir/mplapack/reference/.libs"
for backend in gmp mpfr qd dd; do
    runtime_dirs="$runtime_dirs:$builddir/mpblas/optimized/$backend/.libs"
done
for family in matgen eig lin; do
    runtime_dirs="$runtime_dirs:$builddir/mplapack/test/$family/.libs"
done
for directory in "$@"; do
    runtime_dirs="$runtime_dirs:$directory"
done
checked=0
for family in matgen eig lin; do
    for backend in gmp mpfr qd dd; do
        for suffix in '' _opt; do
            library="$builddir/mplapack/test/$family/.libs/lib${family}_${backend}${suffix}.so"
            test -f "$library" || continue
            checked=$((checked + 1))
            relocations=$(LD_LIBRARY_PATH="$runtime_dirs${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
                ldd -r "$library" 2>&1 || true)
            if printf '%s\n' "$relocations" | grep -q 'not found'; then
                printf '%s\n' "$relocations" >&2
                exit 1
            fi
            unresolved=$(printf '%s\n' "$relocations" | \
                sed -n 's/.*undefined symbol: //p' | sed 's/[[:space:]].*//')
            unexpected=$(printf '%s\n' "$unresolved" | \
                grep -Ev '^(|selval|m|lerr|mplusn|n|fs|nout|nunit|infot|selwr|selopt|seldim|srnamt|ok|selwi|i|iparms|k)$' || true)
            if test -n "$unexpected"; then
                echo "FAIL: $library has unresolved backend/component symbols" >&2
                printf '%s\n' "$unexpected" >&2
                exit 1
            fi
            echo "PASS: $library (test-driver globals excluded)"
        done
    done
done
if test "$checked" -eq 0; then
    echo "FAIL: no backend support libraries found" >&2
    exit 1
fi
