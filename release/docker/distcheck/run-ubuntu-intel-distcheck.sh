#!/bin/bash
set -euxo pipefail

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

rm -rf /work/mplapack /work/mplapack-src
SOURCE_KIND=git
if [ -n "${MPLAPACK_SOURCE_TARBALL:-}" ]; then
    test -f "${MPLAPACK_SOURCE_TARBALL}"
    mkdir -p /work/mplapack-src
    tar xf "${MPLAPACK_SOURCE_TARBALL}" -C /work/mplapack-src
    src_dir="$(find /work/mplapack-src -mindepth 1 -maxdepth 1 -type d | head -n 1)"
    test -n "${src_dir}"
    mv "${src_dir}" /work/mplapack
    cd /work/mplapack
    SOURCE_KIND=tarball
    echo "Using source tarball: ${MPLAPACK_SOURCE_TARBALL}"
else
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
fi

set +u
source /opt/intel/oneapi/setvars.sh
set -u

if [ "$SOURCE_KIND" = "git" ]; then
    cd mplapack/test/compare
    bash gen.Makefile.am.sh
    cd /work/mplapack
    autoreconf --force --install
else
    echo "Using distributed configure files from source tarball; skipping autoreconf."
fi

ARCH=$(dpkg --print-architecture)
echo "Detected architecture: $ARCH"
COMMON_OPTS="--enable-gmp=yes --enable-mpfr=yes --enable-binary128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-test=yes"
if [ "$ARCH" = "amd64" ] || [ "$ARCH" = "i386" ]; then
    CONFIGURE_OPTS="$COMMON_OPTS --enable-benchmark=yes --enable-binary80=yes --with-openblas=-qmkl"
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

OPENBLAS_CONFIGURE_OPT=--with-openblas=-qmkl bash misc/reconfig.ubuntu24.04.intel.sh
make -j"${MAKE_JOBS}"
make install
INSTALL_PREFIX="$(sed -n 's/^prefix = //p' Makefile | head -n 1)"
bash release/check-installed-examples.sh "${INSTALL_PREFIX}" Makefile.linux.inteloneAPI "${MAKE_JOBS}"
bash release/check-installed-benchmarks.sh "${INSTALL_PREFIX}"
make distcheck MAKEFLAGS="-j${MAKE_JOBS}" DISTCHECK_CONFIGURE_FLAGS="$CONFIGURE_OPTS"

echo '=== ccache stats (after) ==='
ccache -s || true
