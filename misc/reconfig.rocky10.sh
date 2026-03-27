#!/bin/sh

USE_CCACHE=yes
CCACHE_VERSION="${CCACHE_VERSION:-4.13.2}"
CCACHE_PREFIX="${HOME}/opt/ccache-${CCACHE_VERSION}"
CCACHE_URL="https://github.com/ccache/ccache/releases/download/v${CCACHE_VERSION}/ccache-${CCACHE_VERSION}.tar.gz"

ensure_ccache() {
    if [ -x "${CCACHE_PREFIX}/bin/ccache" ]; then
        PATH="${CCACHE_PREFIX}/bin:${PATH}"
        export PATH
        return 0
    fi

    if command -v ccache >/dev/null 2>&1; then
        return 0
    fi

    if command -v dnf >/dev/null 2>&1; then
        if [ "$(id -u)" -eq 0 ]; then
            dnf -y install ccache || true
        elif command -v sudo >/dev/null 2>&1; then
            sudo dnf -y install ccache || true
        fi
    fi

    if command -v ccache >/dev/null 2>&1; then
        return 0
    fi

    command -v cmake >/dev/null 2>&1 || {
        echo "Error: cmake is required to build ccache from source." 1>&2
        return 1
    }

    CCACHE_BUILD_ROOT="${HOME}/tmp/ccache-build-${CCACHE_VERSION}"
    CCACHE_TARBALL="${CCACHE_BUILD_ROOT}/ccache-${CCACHE_VERSION}.tar.gz"
    CCACHE_SRC_ROOT="${CCACHE_BUILD_ROOT}/src"
    CCACHE_BUILD_DIR="${CCACHE_BUILD_ROOT}/build"

    rm -rf "${CCACHE_BUILD_ROOT}"
    mkdir -p "${CCACHE_SRC_ROOT}" "${CCACHE_BUILD_DIR}" "${CCACHE_PREFIX}"

    python3 - "${CCACHE_URL}" "${CCACHE_TARBALL}" <<'PY'
import sys
import urllib.request

url = sys.argv[1]
path = sys.argv[2]
with urllib.request.urlopen(url) as response, open(path, "wb") as output:
    output.write(response.read())
PY

    tar -xzf "${CCACHE_TARBALL}" -C "${CCACHE_SRC_ROOT}"
    (
        cd "${CCACHE_BUILD_DIR}" && \
        cmake -D CMAKE_BUILD_TYPE=Release \
              -D CMAKE_INSTALL_PREFIX="${CCACHE_PREFIX}" \
              "${CCACHE_SRC_ROOT}/ccache-${CCACHE_VERSION}" && \
        make -j"$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)" && \
        make install
    ) || return 1

    PATH="${CCACHE_PREFIX}/bin:${PATH}"
    export PATH
    command -v ccache >/dev/null 2>&1
}

if [ x$USE_CCACHE = x"yes" ] ; then
ensure_ccache || exit 1
CXX="ccache g++" ; export CXX
CC="ccache gcc" ; export CC
FC="ccache gfortran"; export FC
F77="ccache gfortran"; export F77
ccache -M 60G
else
CXX="g++" ; export CXX
CC="gcc" ; export CC
FC="gfortran"; export FC
F77="gfortran"; export F77
fi

pushd mplapack/test/compare ; bash gen.Makefile.am.sh ; popd

autoreconf --force --install

arch="$(uname -m)"
case "$arch" in
    x86_64|i386|i486|i586|i686)
        ./configure --prefix="$HOME/MPLAPACK" --enable-gmp=yes --enable-mpfr=yes --enable-binary128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-binary80=yes --enable-test=yes --enable-benchmark=yes
        ;;
    *)
        ./configure --prefix="$HOME/MPLAPACK" --enable-gmp=yes --enable-mpfr=yes --enable-binary128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-test=yes --enable-benchmark=yes
        ;;
esac
