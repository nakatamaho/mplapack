#!/bin/bash

USE_CCACHE=yes
source /opt/intel/oneapi/setvars.sh
OPENBLAS_CONFIGURE_OPT="${OPENBLAS_CONFIGURE_OPT:---with-openblas=-qmkl}"

#https://gmplib.org/list-archives/gmp-bugs/2014-September/003526.html
if [ x$USE_CCACHE = x"yes" ] ; then
CXX="ccache icpx" ; export CXX
CC="ccache icx" ; export CC
FC="ccache ifx"; export FC
F77="ccache ifx"; export F77
"$(cd "$(dirname "$0")" && pwd)/ccache_set_min_maxsize.sh" 80G
else
CXX="icpx" ; export CXX
CC="icx" ; export CC
FC="ifx"; export FC
F77="ifx"; export F77
fi

pushd mplapack/test/compare ; bash gen.Makefile.am.sh ; popd

autoreconf --force --install
./configure --prefix=$HOME/MPLAPACK_INTELONEAPI --enable-gmp=yes --enable-mpfr=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-test=yes --enable-benchmark=yes --enable-binary128=yes --enable-binary80=yes "$OPENBLAS_CONFIGURE_OPT"
