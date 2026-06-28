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

rm -rf /work/mplapack
if [ -f "${MPLAPACK_REPO}" ]; then
    git clone --no-checkout "${MPLAPACK_REPO}" /work/mplapack
    cd /work/mplapack
    git checkout "${MPLAPACK_REF}"
else
    git clone --depth 1 --branch "${MPLAPACK_REF}" "${MPLAPACK_REPO}" /work/mplapack || {
        git clone "${MPLAPACK_REPO}" /work/mplapack
        cd /work/mplapack
        git checkout "${MPLAPACK_REF}"
    }
    cd /work/mplapack
fi
git log -1

set +u
source /opt/intel/oneapi/setvars.sh
set -u

cd mplapack/test/compare
bash gen.Makefile.am.sh
cd /work/mplapack
autoreconf --force --install

OPENBLAS_CONFIGURE_OPT=--with-openblas=-qmkl bash misc/reconfig.ubuntu24.04.intel.sh
make -j"${MAKE_JOBS}"
make install

echo '=== ccache stats (after) ==='
ccache -s || true
