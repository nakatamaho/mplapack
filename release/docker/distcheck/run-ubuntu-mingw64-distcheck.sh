#!/bin/bash
set -euxo pipefail

echo '=== toolchain ==='
command -v x86_64-w64-mingw32-gcc-posix
command -v x86_64-w64-mingw32-g++-posix
command -v x86_64-w64-mingw32-gfortran-posix
x86_64-w64-mingw32-gcc-posix -dumpmachine
x86_64-w64-mingw32-g++-posix -dumpversion

get_make_jobs() {
    if [ -n "${MPLAPACK_MAKE_JOBS:-}" ]; then
        printf '%s\n' "${MPLAPACK_MAKE_JOBS}"
        return
    fi
    if command -v lscpu >/dev/null 2>&1; then
        cores="$(lscpu | awk -F: '/^Core\(s\) per socket:/ {gsub(/^[ \t]+/, "", $2); print $2; exit}')"
        sockets="$(lscpu | awk -F: '/^Socket\(s\):/ {gsub(/^[ \t]+/, "", $2); print $2; exit}')"
        if [ -n "${cores}" ] && [ -n "${sockets}" ] && [ "${cores}" -gt 0 ] 2>/dev/null && [ "${sockets}" -gt 0 ] 2>/dev/null; then
            printf '%s\n' "$((cores * sockets))"
            return
        fi
    fi
    nproc
}

MAKE_JOBS="$(get_make_jobs)"
echo "MAKE_JOBS=${MAKE_JOBS}"

echo '=== ccache stats (before) ==='
ccache -s || true

rm -rf /work/mplapack
if [ -f "$MPLAPACK_REPO" ]; then
    git clone --no-checkout "$MPLAPACK_REPO" /work/mplapack
    cd /work/mplapack
    git checkout "$MPLAPACK_REF"
else
    git clone --depth 1 --branch "$MPLAPACK_REF" "$MPLAPACK_REPO" /work/mplapack || {
        git clone "$MPLAPACK_REPO" /work/mplapack
        cd /work/mplapack
        git checkout "$MPLAPACK_REF"
    }
    cd /work/mplapack
fi
git log -1

cd mplapack/test/compare
bash gen.Makefile.am.sh
cd /work/mplapack
autoreconf --force --install

ARCH=$(dpkg --print-architecture)
BUILD_TRIPLE=$(gcc -dumpmachine)
echo "Detected architecture: $ARCH"
echo "Detected build triple: $BUILD_TRIPLE"

echo '=== Preparing Wine runtime ==='
GCC_VER=$(x86_64-w64-mingw32-gcc-posix -dumpversion | cut -d- -f1)
export WINEPATH="/usr/x86_64-w64-mingw32/lib/;/usr/lib/gcc/x86_64-w64-mingw32/${GCC_VER}-posix/;/usr/local/bin"
echo "WINEPATH=$WINEPATH"
echo '=== Enabling Wine dot files ==='
before=$(stat -c '%Y' "$HOME/.wine/user.reg" 2>/dev/null || echo 0)
wine reg add 'HKEY_CURRENT_USER\Software\Wine' /v ShowDotFiles /d Y /f >/dev/null
if [ -f "$HOME/.wine/user.reg" ]; then
    while [ "$(stat -c '%Y' "$HOME/.wine/user.reg")" = "$before" ]; do sleep 1; done
fi
echo '=== Checking winepath ==='
command -v winepath
WP_TEST=$(winepath -w /work 2>/dev/null || true)
test -n "$WP_TEST"
echo "winepath(/work)=$WP_TEST"

COMMON_OPTS="--host=x86_64-w64-mingw32 --build=$BUILD_TRIPLE --enable-gmp=yes --enable-mpfr=yes --enable-binary128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-test=yes"
if [ "$ARCH" = "amd64" ] || [ "$ARCH" = "i386" ]; then
    CONFIGURE_OPTS="$COMMON_OPTS --enable-benchmark=yes --enable-binary80=yes"
else
    CONFIGURE_OPTS="$COMMON_OPTS"
fi
echo "Configure options: $CONFIGURE_OPTS"

RESULTS_VERSION=$(sed -n 's/^m4_define(\[MPLAPACK_VER_MAJOR\], \[\([^]]*\)\])/\1/p' configure.ac | tr -d '[:space:]').$(sed -n 's/^m4_define(\[MPLAPACK_VER_MINOR\], \[\([^]]*\)\])/\1/p' configure.ac | tr -d '[:space:]').$(sed -n 's/^m4_define(\[MPLAPACK_VER_PATCH\], \[\([^]]*\)\])/\1/p' configure.ac | tr -d '[:space:]')
: "${MPLAPACK_TEST_RESULTS_BASE:=/results}"
export MPLAPACK_TEST_RESULTS_STAGING="$MPLAPACK_TEST_RESULTS_BASE/$RESULTS_VERSION"
rm -rf "$MPLAPACK_TEST_RESULTS_STAGING"
mkdir -p "$MPLAPACK_TEST_RESULTS_STAGING"
echo "MPLAPACK_TEST_RESULTS_STAGING=$MPLAPACK_TEST_RESULTS_STAGING"

./configure $CONFIGURE_OPTS
make -j"${MAKE_JOBS}"
make install
INSTALL_PREFIX="$(sed -n 's/^prefix = //p' Makefile | head -n 1)"
bash release/check-installed-examples.sh "${INSTALL_PREFIX}" Makefile.mingw "${MAKE_JOBS}" CXX=x86_64-w64-mingw32-g++-posix
bash release/check-installed-benchmarks.sh "${INSTALL_PREFIX}"

echo '=== Running make distcheck ==='
env CC="ccache x86_64-w64-mingw32-gcc-posix" CXX="ccache x86_64-w64-mingw32-g++-posix" FC="ccache x86_64-w64-mingw32-gfortran-posix" \
    make distcheck LOG_COMPILER=wine MAKEFLAGS="-j${MAKE_JOBS}" DISTCHECK_CONFIGURE_FLAGS="$CONFIGURE_OPTS"

echo '=== ccache stats (after) ==='
ccache -s || true
