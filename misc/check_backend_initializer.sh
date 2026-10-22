#!/bin/sh
# Test actual libQD3 backend constructor extraction from a static archive.
# Usage: CPPFLAGS='-I<qd/include>' LIBS='-L<qd/lib> -lqd' \
#   sh misc/check_backend_initializer.sh <top-srcdir> <backend>
# For MinGW set CXX, AR, NM, LOG_COMPILER=wine and the required WINEPATH.
set -eu
test "$#" -eq 2 || exit 2
srcdir=$(cd "$1" && pwd)
backend=$2
# Validate the name without maintaining a second backend table.
case "$backend" in ''|*[!a-z0-9]*) exit 2 ;; esac
initializer="mplapack_initialize_$backend"
grep -Fq "void $initializer(void) { fpu_fix_start(" \
    "$srcdir/mpblas/reference/mplapackinit.cpp" || {
    echo "SKIP: $backend has no libQD3 FPU constructor"; exit 77;
}
cxx=${CXX:-c++}
ar=${AR:-ar}
nm=${NM:-nm}
runner=${LOG_COMPILER:-}
tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/mplapack-backend-init.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT HUP INT TERM
macro=$(printf '%s' "$backend" | tr '[:lower:]' '[:upper:]')
$cxx ${CPPFLAGS:-} ${CXXFLAGS:-} -std=gnu++17 -ffp-contract=off \
    -DMPLAPACK_BUILD_WITH_${macro} -I"$srcdir/include" \
    -c "$srcdir/mpblas/reference/mplapackinit.cpp" -o "$tmpdir/init.o"
$ar rcs "$tmpdir/libinit.a" "$tmpdir/init.o"
NM="$nm" sh "$srcdir/misc/check_unique_symbols.sh" "$tmpdir/libinit.a"
# Read the object symbol rather than assuming an ELF/Mach-O/PE prefix.
symbol=$($nm "$tmpdir/init.o" | awk -v name="$initializer" \
    '$NF ~ name && $(NF-1) == "T" {print $NF}')
test -n "$symbol"
for mode in ordinary retained; do
    retain=
    test "$mode" != retained || retain="-Wl,-u,$symbol"
    $cxx ${CXXFLAGS:-} -std=gnu++17 -ffp-contract=off \
        "$srcdir/misc/test_backend_initializer.cpp" $retain \
        "$tmpdir/libinit.a" ${LDFLAGS:-} ${LIBS:--lqd} \
        -o "$tmpdir/$mode.exe"
    rc=0
    $runner "$tmpdir/$mode.exe" || rc=$?
    test "$rc" -ne 77 || { echo "SKIP: x87 probe requires x86"; exit 77; }
    expected=1
    test "$mode" != retained || expected=0
    test "$rc" -eq "$expected" || {
        echo "FAIL: $backend $mode returned $rc, expected $expected" >&2
        exit 1
    }
done
echo "PASS: $backend constructor retained; unretained negative control fails"
