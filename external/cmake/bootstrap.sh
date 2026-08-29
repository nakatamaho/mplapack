#!/bin/sh
set -eu

version=3.20.6
expected_sha256=a0bd485e1a38dd13c0baec89d5f4adbf61c7fd32fddb38eabc69a75bc0b65d72

if test "$#" -ne 4; then
    echo "usage: $0 ARCHIVE WORKDIR PREFIX MAKE" >&2
    exit 2
fi

archive=$1
workdir=$2
prefix=$3
make_command=$4

jobs=${CMAKE_BUILD_PARALLEL_LEVEL-}
if test -z "$jobs"; then
    jobs=`getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1`
fi
case "$jobs" in
    ''|*[!0-9]*) jobs=1 ;;
esac
if test "$jobs" -lt 1; then
    jobs=1
fi

case "$workdir" in
    ""|/)
        echo "error: refusing to use unsafe CMake work directory '$workdir'" >&2
        exit 2
        ;;
esac

if test ! -f "$archive"; then
    echo "error: bundled CMake archive not found: $archive" >&2
    exit 1
fi

if command -v sha256sum >/dev/null 2>&1; then
    printf '%s  %s\n' "$expected_sha256" "$archive" | sha256sum -c -
elif command -v shasum >/dev/null 2>&1; then
    printf '%s  %s\n' "$expected_sha256" "$archive" | shasum -a 256 -c -
else
    echo "error: no SHA256 checksum tool found (sha256sum or shasum)" >&2
    exit 1
fi

rm -rf "$workdir"
mkdir -p "$workdir/src"
tar -xzf "$archive" -C "$workdir/src"
source_dir="$workdir/src/cmake-$version"

if test ! -d "$source_dir"; then
    echo "error: bundled CMake archive has no cmake-$version source directory" >&2
    exit 1
fi

rm -rf "$prefix"
mkdir -p "$prefix"
cd "$source_dir"

./bootstrap --prefix="$prefix" --parallel="$jobs" -- \
    -DBUILD_TESTING=OFF \
    -DCMAKE_USE_OPENSSL=OFF

make_args=${MAKEFLAGS-}
if test -n "$make_args"; then
    # shellcheck disable=SC2086
    "$make_command" -C "$source_dir" -j"$jobs" $make_args
    # shellcheck disable=SC2086
    "$make_command" -C "$source_dir" -j"$jobs" $make_args install
else
    "$make_command" -C "$source_dir" -j"$jobs"
    "$make_command" -C "$source_dir" -j"$jobs" install
fi

test -x "$prefix/bin/cmake"
test -x "$prefix/bin/ctest"
