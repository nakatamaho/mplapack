#!/bin/sh

USE_CCACHE=yes

if [ x$USE_CCACHE = x"yes" ] ; then
CXX="ccache /opt/rh/devtoolset-9/root/usr/bin/g++" ; export CXX
CC="ccache /opt/rh/devtoolset-9/root/usr/bin/gcc" ; export CC
FC="/opt/rh/devtoolset-9/root/usr/bin/gfortran"; export FC
F77="/opt/rh/devtoolset-9/root/usr/bin/gfortran"; export F77
ccache -M 20G
else
CXX="/opt/rh/devtoolset-9/root/usr/bin/g++" ; export CXX
CC="/opt/rh/devtoolset-9/root/usr/bin/gcc" ; export CC
FC="/opt/rh/devtoolset-9/root/usr/bin/gfortran"; export FC
F77="/opt/rh/devtoolset-9/root/usr/bin/gfortran"; export F77
fi

pushd mplapack/test/compare ; bash gen.Makefile.am.sh ; popd

aclocal ; autoconf ; automake --add-missing
autoreconf --force --install
if [ `uname -m` = "x86_64" ]; then
./configure --prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-_Float128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-_Float64x=yes --enable-test=yes
else
./configure --prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-_Float128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-test=yes
fi
