#!/bin/sh

USE_CCACHE=yes
if [ "x$USE_CCACHE" = x"yes" ] ; then
command -v ccache >/dev/null 2>&1 || {
    echo "Error: ccache is not installed in the image." 1>&2
    exit 1
}
CXX="ccache g++" ; export CXX
CC="ccache gcc" ; export CC
FC="ccache gfortran"; export FC
F77="ccache gfortran"; export F77
"$(cd "$(dirname "$0")" && pwd)/ccache_set_min_maxsize.sh" 80G
else
CXX="g++" ; export CXX
CC="gcc" ; export CC
FC="gfortran"; export FC
F77="gfortran"; export F77
fi

pushd mplapack/test/compare ; bash gen.Makefile.am.sh ; popd

autoreconf --force --install
arch="$(uname -m)"
enable_benchmark="${ENABLE_BENCHMARK:-yes}"
configure_opts="--prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-binary128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-test=yes"

case "$arch" in
    x86_64|i386|i486|i586|i686)
        configure_opts="$configure_opts --enable-binary80=yes"
        ;;
    *)
        ;;
esac

case "$enable_benchmark" in
    yes)
        configure_opts="$configure_opts --enable-benchmark=yes --with-openblas=system"
        ;;
    no)
        ;;
    *)
        echo "Error: ENABLE_BENCHMARK must be yes or no, got: $enable_benchmark" >&2
        exit 1
        ;;
esac

./configure $configure_opts
