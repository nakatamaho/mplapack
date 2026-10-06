#!/bin/sh
# Check LIN/EIG test override ownership against built shared libraries.
# Usage: check_test_override_link.sh <top-builddir> <top-srcdir>
set -eu
test "$#" -eq 2 || { echo "usage: $0 <top-builddir> <top-srcdir>" >&2; exit 2; }
builddir=$(cd "$1" && pwd)
srcdir=$(cd "$2" && pwd)
case $(uname -s) in
    Linux) ;;
    *) echo "SKIP: override as-needed regression requires Linux"; exit 77 ;;
esac
cxx=${CXX:-c++}
command -v readelf >/dev/null 2>&1 || exit 77
pkg_config=${PKG_CONFIG:-pkg-config}
tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/mplapack-test-overrides.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT HUP INT TERM
export PKG_CONFIG_PATH="$builddir${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
checked=0
for backend in gmp mpfr qd dd; do
    test -f "$builddir/mplapack_$backend.pc" || continue
    checked=$((checked + 1))
    macro=$(printf '%s' "$backend" | tr '[:lower:]' '[:upper:]')
    cflags=$($pkg_config --cflags "mplapack_$backend")
    for family in lin eig; do
        objects=
        for source in "$srcdir/misc/test_override_link.cpp" \
            "$srcdir/mplapack/test/$family/common/Mxerbla.cpp" \
            "$srcdir/mplapack/test/$family/common/iMlaenv.cpp" \
            "$srcdir/mplapack/test/$family/common/Mxlaenv.cpp"; do
            object="$tmpdir/${backend}_${family}_$(basename "$source" .cpp).o"
            $cxx ${CXXFLAGS:-} -std=gnu++17 -I"$builddir/include" \
                -I"$srcdir/include" -I"$srcdir/fable"  \
                -DMPLAPACK_BUILD_WITH_${macro} -DMPLAPACK_INTERNAL \
                $cflags -c "$source" -o "$object"
            objects="$objects $object"
        done
        if test "$family" = eig; then
            $cxx ${CXXFLAGS:-} -std=gnu++17 -I"$builddir/include" \
                -I"$srcdir/include" -I"$srcdir/fable"  \
                -DMPLAPACK_BUILD_WITH_${macro} -DMPLAPACK_INTERNAL \
                -DTEST_DMD $cflags -c "$srcdir/misc/test_override_link.cpp" \
                -o "$tmpdir/${backend}_dmd.o"
        fi
        for suffix in '' _opt; do
            primary="$builddir/mplapack/reference/.libs"
            if test -n "$suffix"; then
                primary="$builddir/mpblas/optimized/$backend/.libs"
            fi
            support="$builddir/mplapack/test/$family/.libs"
            matgen="$builddir/mplapack/test/matgen/.libs"
            libs=$($pkg_config --libs "mplapack_${backend}${suffix}")
            executable="$tmpdir/consumer"
            $cxx $objects -Wl,--as-needed -L"$support" \
                -l${family}_${backend}${suffix} -l${family}_override_${backend}${suffix} \
                -L"$matgen" -lmatgen_${backend}${suffix} -L"$primary" $libs -o "$executable"
            LD_LIBRARY_PATH="$support:$matgen:$primary${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
                "$executable"
            echo "PASS: $family $backend$suffix executable-owned overrides"
            if test "$family" = eig; then
                # Mirror the DMD graph: no override archive or test iMlaenv object.
                $cxx "$tmpdir/${backend}_dmd.o" \
                    "$tmpdir/${backend}_eig_Mxerbla.o" \
                    "$tmpdir/${backend}_eig_Mxlaenv.o" -Wl,--as-needed \
                    -L"$support" -leig_${backend}${suffix} \
                    -L"$matgen" -lmatgen_${backend}${suffix} \
                    -L"$primary" $libs -o "$executable"
                if readelf -d "$support/libeig_${backend}${suffix}.so" | grep -q 'NEEDED.*override'; then
                    echo "FAIL: DMD support graph pulls in a test override" >&2
                    exit 1
                fi
                LD_LIBRARY_PATH="$support:$matgen:$primary${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
                    "$executable"
                echo "PASS: dmd $backend$suffix production iMlaenv (ISPEC=9: 25)"
            fi
        done
    done
done

test "$checked" -gt 0 || { echo "FAIL: no configured backend found" >&2; exit 1; }

# Exercise actual LIN drivers when they have been built. Zero-size matrices
# retain error-exit checks while avoiding the release numerical workload.
for suffix in '' _opt; do
    driver="$builddir/mplapack/test/lin/mpfr/xlintstR_mpfr${suffix}"
    test -x "$driver" || continue
    output="$tmpdir/lin-error-exits.out"
    MPFRXX_DEFAULT_EMIN=-32765 MPFRXX_DEFAULT_EMAX=32768 \
        "$driver" < "$srcdir/misc/test_override_error.in" > "$output" 2>&1
    if ! grep -q 'DGE routines passed the tests of the error exits' "$output"; then
        cat "$output" >&2
        exit 1
    fi
    echo "PASS: actual LIN mpfr$suffix DGE error exits (zero-size matrices)"
done
