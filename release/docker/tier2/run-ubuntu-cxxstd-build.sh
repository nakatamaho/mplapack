#!/bin/bash
set -euxo pipefail

get_make_jobs() {
    if [ -n "${MPLAPACK_MAKE_JOBS:-}" ]; then
        printf '%s\n' "${MPLAPACK_MAKE_JOBS}"
        return
    fi
    nproc
}

prepare_source() {
    rm -rf /work/mplapack /work/mplapack-src /work/mplapack-source
    SOURCE_KIND=git
    if [ -n "${MPLAPACK_SOURCE_TARBALL:-}" ]; then
        test -f "${MPLAPACK_SOURCE_TARBALL}"
        mkdir -p /work/mplapack-src
        tar xf "${MPLAPACK_SOURCE_TARBALL}" -C /work/mplapack-src
        src_dir="$(find /work/mplapack-src -mindepth 1 -maxdepth 1 -type d | head -n 1)"
        test -n "${src_dir}"
        mv "${src_dir}" /work/mplapack-source
        SOURCE_KIND=tarball
        echo "Using source tarball: ${MPLAPACK_SOURCE_TARBALL}"
    elif [ -f "${MPLAPACK_REPO}" ]; then
        /usr/local/bin/checkout-source.sh "${MPLAPACK_REPO}" "${MPLAPACK_REF}" /work/mplapack-source
    else
        /usr/local/bin/checkout-source.sh "${MPLAPACK_REPO}" "${MPLAPACK_REF}" /work/mplapack-source
    fi
    cd /work/mplapack-source
    if [ "${SOURCE_KIND}" = "git" ]; then
        git log -1
    fi
}

MAKE_JOBS="$(get_make_jobs)"
echo "MAKE_JOBS=${MAKE_JOBS}"

echo '=== ccache stats (before) ==='
ccache -s || true

prepare_source

RECONFIG_SCRIPTS="${MPLAPACK_TIER2_RECONFIG_SCRIPTS:-misc/reconfig.ubuntu26.04.c++17.sh misc/reconfig.ubuntu26.04.c++20.sh misc/reconfig.ubuntu26.04.c++23.sh misc/reconfig.ubuntu26.04.c++26.sh misc/reconfig.ubuntu26.04.gnuc++20.sh misc/reconfig.ubuntu26.04.gnuc++23.sh misc/reconfig.ubuntu26.04.gnuc++26.sh}"

for reconfig_script in ${RECONFIG_SCRIPTS}; do
    echo "=== TIER2 BUILD START: ${reconfig_script} ==="
    rm -rf /work/mplapack
    mkdir -p /work/mplapack
    cp -a /work/mplapack-source/. /work/mplapack/
    cd /work/mplapack
    test -f "${reconfig_script}"
    bash -x "${reconfig_script}"
    make -j"${MAKE_JOBS}"
    echo "=== TIER2 BUILD OK: ${reconfig_script} ==="
done

echo '=== ccache stats (after) ==='
ccache -s || true
echo "=== ALL TIER2 C++ STANDARD BUILDS COMPLETED SUCCESSFULLY ==="
