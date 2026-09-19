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
/usr/local/bin/checkout-source.sh "${MPLAPACK_REPO}" "${MPLAPACK_REF}" /work/mplapack
cd /work/mplapack
git log -1

autoreconf --force --install

ARCH=$(dpkg --print-architecture)
BUILD_TRIPLE=$(gcc -dumpmachine)
echo "Detected architecture: ${ARCH}"
echo "Detected build triple: ${BUILD_TRIPLE}"

COMMON_OPTS="--host=x86_64-w64-mingw32 --build=${BUILD_TRIPLE} --enable-gmp=yes --enable-mpfr=yes --enable-binary128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-test=yes"
if [ "${ARCH}" = "amd64" ] || [ "${ARCH}" = "i386" ]; then
    CONFIGURE_OPTS="${COMMON_OPTS} --enable-benchmark=yes --enable-binary80=yes"
else
    CONFIGURE_OPTS="${COMMON_OPTS}"
fi
echo "Configure options: ${CONFIGURE_OPTS}"

./configure ${CONFIGURE_OPTS}
make -j"${MAKE_JOBS}"
make install

echo '=== ccache stats (after) ==='
ccache -s || true
