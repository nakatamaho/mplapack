#!/bin/sh
# Verify test-library ownership and the MinGW object-before-archive link order.
# Usage: sh misc/check_test_link_order.sh <top-srcdir> [--mingw]
set -eu
test "$#" -ge 1 && test "$#" -le 2 || exit 2
srcdir=$(cd "$1" && pwd)
if test -f "$srcdir/backends.txt"; then
    backends=$(awk '!/^#/ && NF { print $1 }' "$srcdir/backends.txt")
else
    backends="gmp mpfr qd dd double binary80 binary128"
fi
for family in lin eig; do
    for backend in $backends; do
        makefile="$srcdir/mplapack/test/$family/$backend/Makefile.am"
        if grep -Eq '_LDFLAGS.*\$\((.*libdepends|mplibs)\)|whole-archive' "$makefile"; then
            echo "FAIL: libraries precede executable objects in $makefile" >&2
            exit 1
        fi
        awk '
            /^x.*_SOURCES[[:space:]]*=/ {
                target=$1; sub(/_SOURCES$/, "", target); sources[target]=1
            }
            /^x.*_LDADD[[:space:]]*=/ {
                target=$1; sub(/_LDADD$/, "", target); deps[target]=1
            }
            END {
                for (target in sources)
                    if (!deps[target]) {
                        print "FAIL: missing LDADD for " target; failed=1
                    }
                exit failed
            }
        ' "$makefile"
        generated="${makefile%.am}.in"
        if test -f "$generated"; then
            awk '
                /AM_V_CXXLD/ && /_OBJECTS/ {
                    if (!index($0, "_LDADD") ||
                        index($0, "_OBJECTS") > index($0, "_LDADD")) {
                        print "FAIL: generated recipe puts libraries first: " $0
                        failed=1
                    }
                    count++
                }
                END { exit failed || count == 0 }
            ' "$generated"
        fi
    done
done
echo "PASS: LIN/EIG reference/optimized libraries follow executable objects"
test "${2:-}" = --mingw || exit 0
cxx=${MINGW_CXX:-x86_64-w64-mingw32-g++-posix}
ar=${MINGW_AR:-x86_64-w64-mingw32-ar}
runner=${WINE:-wine}
command -v "$cxx" >/dev/null
command -v "$ar" >/dev/null
command -v "$runner" >/dev/null
tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/mplapack-mingw-overrides.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT HUP INT TERM
# The two release lines use different backend macro conventions.
macro=___MPLAPACK_BUILD_WITH_DOUBLE___
if grep -q 'defined MPLAPACK_BUILD_WITH_DOUBLE' "$srcdir/include/mpblas.h"; then
    macro=MPLAPACK_BUILD_WITH_DOUBLE
fi
"$cxx" -std=gnu++17 -D"$macro" -I"$srcdir/include" -I"$srcdir/fable" \
    -c "$srcdir/misc/test_mingw_override.cpp" -o "$tmpdir/main.o"
"$cxx" -std=gnu++17 -D"$macro" -I"$srcdir/include" -I"$srcdir/fable" \
    -c "$srcdir/mpblas/reference/Mxerbla.cpp" -o "$tmpdir/production.o"
"$ar" rcs "$tmpdir/libproduction.a" "$tmpdir/production.o"
for family in lin eig; do
    "$cxx" -std=gnu++17 -D"$macro" -I"$srcdir/include" -I"$srcdir/fable" \
        -c "$srcdir/mplapack/test/$family/common/Mxerbla.cpp" -o "$tmpdir/test.o"
    "$cxx" -Wl,--allow-multiple-definition "$tmpdir/main.o" "$tmpdir/test.o" \
        "$tmpdir/libproduction.a" -o "$tmpdir/objects-first.exe"
    "$runner" "$tmpdir/objects-first.exe"
    # Negative control: reproduce the previous release's archive-first link.
    "$cxx" -Wl,--allow-multiple-definition \
        -Wl,--whole-archive,"$tmpdir/libproduction.a",--no-whole-archive \
        "$tmpdir/main.o" "$tmpdir/test.o" -o "$tmpdir/archive-first.exe"
    rc=0
    "$runner" "$tmpdir/archive-first.exe" || rc=$?
    test "$rc" -eq 1 || {
        echo "FAIL: archive-first negative control returned $rc, expected 1" >&2
        exit 1
    }
    echo "PASS: MinGW/Wine $family test Mxerbla wins; old order exits 1"
done
