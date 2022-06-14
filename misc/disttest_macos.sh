TOPDIR=/Volumes/Users/maho/tmp/mplapack
VERSION=2.0.0
####################
cd $TOPDIR
pushd mplapack/test/compare ; bash gen.Makefile.am.sh ; popd
aclocal ; autoconf ; automake --add-missing
autoreconf --force --install

CXX="ccache g++-mp-10" ; export CXX
CC="ccache gcc-mp-10" ; export CC
FC="ccache gfortran-mp-10"; export FC

./configure --prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-_Float128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-_Float64x=yes --enable-test=yes
make dist-xz
####################
rm -rf $TOPDIR/tmp
rm -rf $HOME/MPLAPACK
mkdir $TOPDIR/tmp
cd $TOPDIR/tmp
rm mplapack-${VERSION}.tar.xz
cp $TOPDIR/mplapack-${VERSION}.tar.xz .
tar xvfJ mplapack-${VERSION}.tar.xz
cd mplapack-${VERSION}

./configure --prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-_Float128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-_Float64x=yes --enable-test=yes
make -j4
make install
