#!/bin/sh
# Check the shared dependency closure of backend test/support libraries.
# usage: check_backend_support_libraries.sh <top-builddir> <top-srcdir> [dependency-dir ...]
set -eu

if test "$#" -lt 2; then
    echo "usage: $0 <top-builddir> <top-srcdir> [dependency-dir ...]" >&2
    exit 2
fi

top_builddir=$1
top_srcdir=$2
shift 2
runtime_dirs="$top_builddir/mplapack/reference/.libs"
runtime_dirs="$runtime_dirs:$top_builddir/mpblas/optimized/gmp/.libs"
runtime_dirs="$runtime_dirs:$top_builddir/mpblas/optimized/qd/.libs"
runtime_dirs="$runtime_dirs:$top_builddir/mpblas/optimized/dd/.libs"
runtime_dirs="$runtime_dirs:$top_builddir/mplapack/test/matgen/.libs"
runtime_dirs="$runtime_dirs:$top_builddir/mplapack/test/lin/.libs"
runtime_dirs="$runtime_dirs:$top_builddir/mplapack/test/eig/.libs"
for dependency_dir in "$@"; do
    if test -n "$dependency_dir"; then
        runtime_dirs="$runtime_dirs:$dependency_dir"
    fi
done

status=0
for suite in matgen lin eig; do
    library_dir="$top_builddir/mplapack/test/$suite/.libs"
    for backend in gmp gmp_opt qd qd_opt dd dd_opt; do
        case "$backend" in
            gmp*) dependency=gmp ;;
            qd*|dd*) dependency=qd ;;
        esac
        for library in \
            "$library_dir"/lib${suite}_${backend}.so \
            "$library_dir"/lib${suite}_${backend}.dylib \
            "$library_dir"/lib${suite}_${backend}.*.dylib \
            "$library_dir"/lib${suite}_${backend}.dll; do
            if test ! -f "$library"; then
                continue
            fi
            echo "Checking $library"
            if test "$(uname -s)" = Linux; then
                if ! command -v readelf >/dev/null 2>&1 ||
                   ! command -v ldd >/dev/null 2>&1; then
                    echo "SKIP: ELF tools are required for $library"
                    continue
                fi
                needed=$(readelf -d "$library")
                if ! printf '%s\n' "$needed" |
                    grep -Eq "Shared library: \[lib${dependency}\.so"; then
                    echo "FAIL: $library has no DT_NEEDED entry for lib${dependency}.so" >&2
                    status=1
                fi
                relocations=$(LD_LIBRARY_PATH="$runtime_dirs${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
                    ldd -r "$library" 2>&1 || true)
                if printf '%s\n' "$relocations" | grep -Eq 'not found'; then
                    echo "FAIL: missing runtime dependency for $library" >&2
                    printf '%s\n' "$relocations" >&2
                    status=1
                fi
                # These support libraries are test harness components.  They
                # intentionally refer to harness globals and override
                # routines supplied by the final test executable.  Reject
                # only unresolved symbols that the declared backend library
                # itself should provide; those indicate underlinking.
                unresolved=$(printf '%s\n' "$relocations" |
                    sed -n 's/.*undefined symbol: //p')
                dependency_paths=$(printf '%s\n' "$relocations" |
                    awk -v soname="^lib${dependency}\\.so" '$1 ~ soname {print $3}')
                for dependency_path in $dependency_paths; do
                    if ! command -v nm >/dev/null 2>&1; then
                        continue
                    fi
                    exported=$(nm -D --defined-only -P "$dependency_path" 2>/dev/null |
                        awk '{print $1}')
                    for symbol in $unresolved; do
                        if printf '%s\n' "$exported" | grep -Fxq "$symbol"; then
                            echo "FAIL: $library has unresolved backend symbol $symbol" >&2
                            status=1
                        fi
                    done
                done
                echo "PASS: backend dependency and relocation check for $library"
            fi
        done
    done
done
exit "$status"
