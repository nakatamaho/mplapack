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

get_distro_version() {
    if [ -n "${MPLAPACK_DISTRO_VERSION:-}" ]; then
        printf '%s\n' "${MPLAPACK_DISTRO_VERSION}"
        return
    fi
    if [ -r /etc/os-release ]; then
        . /etc/os-release
        printf '%s\n' "${VERSION_ID:-}"
        return
    fi
    printf '%s\n' "unknown"
}

get_m4_value() {
    sed -n "s/^m4_define(\[$1\], \[\([^]]*\)\])/\1/p" configure.ac | tr -d '[:space:]'
}

MAKE_JOBS="$(get_make_jobs)"
DISTRO_VERSION="$(get_distro_version)"
RECONFIG_SCRIPT="misc/reconfig.ubuntu${DISTRO_VERSION}.nvidia.sh"
NVIDIA_CONFIGURE_OPTS="--enable-gmp=yes --enable-mpfr=yes --enable-binary128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-binary80=yes --enable-cuda=yes --enable-test=yes --enable-benchmark=yes"

echo "MAKE_JOBS=${MAKE_JOBS}"
echo "MPLAPACK_DISTRO_VERSION=${DISTRO_VERSION}"
echo "RECONFIG_SCRIPT=${RECONFIG_SCRIPT}"
echo '=== nvcc version ==='
nvcc --version

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

if [ ! -f "$RECONFIG_SCRIPT" ]; then
    echo "ERROR: NVIDIA reconfig script not found: $RECONFIG_SCRIPT" >&2
    exit 1
fi

RESULTS_VERSION="$(get_m4_value MPLAPACK_VER_MAJOR).$(get_m4_value MPLAPACK_VER_MINOR).$(get_m4_value MPLAPACK_VER_PATCH)"
: "${MPLAPACK_TEST_RESULTS_BASE:=/results}"
export MPLAPACK_TEST_RESULTS_STAGING="$MPLAPACK_TEST_RESULTS_BASE/$RESULTS_VERSION"
rm -rf "$MPLAPACK_TEST_RESULTS_STAGING"
mkdir -p "$MPLAPACK_TEST_RESULTS_STAGING"
echo "MPLAPACK_TEST_RESULTS_STAGING=$MPLAPACK_TEST_RESULTS_STAGING"

if [ "$SOURCE_KIND" = "git" ]; then
    bash "$RECONFIG_SCRIPT"
else
    echo "Using distributed configure files from source tarball; skipping autoreconf."
    env CC="ccache gcc" CXX="ccache g++" FC="gfortran" F77="gfortran" \
        ./configure --prefix="$HOME/MPLAPACK_NVIDIA" $NVIDIA_CONFIGURE_OPTS --with-opencl=/usr/local/cuda
fi
INSTALL_PREFIX="$(sed -n 's/^prefix = //p' Makefile | head -n 1)"
make -j"${MAKE_JOBS}"
make install
bash release/check-installed-examples.sh "${INSTALL_PREFIX}" Makefile.linux_cuda,Makefile.linux "${MAKE_JOBS}"
bash release/check-installed-benchmarks.sh "${INSTALL_PREFIX}"

echo "DISTCHECK_CONFIGURE_FLAGS=$NVIDIA_CONFIGURE_OPTS"
echo '=== Running make distcheck with NVIDIA reconfig options ==='
make distcheck MAKEFLAGS="-j${MAKE_JOBS}" DISTCHECK_CONFIGURE_FLAGS="$NVIDIA_CONFIGURE_OPTS"

echo '=== ccache stats (after) ==='
ccache -s || true
