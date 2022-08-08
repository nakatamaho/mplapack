TOPDIR=/home/docker/mplapack
VERSION=2.0.1
####################
cd $TOPDIR
pushd mplapack/test/compare ; bash gen.Makefile.am.sh ; popd
aclocal ; autoconf ; automake --add-missing
autoreconf --force --install
./configure --prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-_Float128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-_Float64x=yes --enable-test=yes --enable-cuda
make dist-xz
####################
rm -rf $TOPDIR/tmp
rm -rf $HOME/MPLAPACK
mkdir $TOPDIR/tmp
cd $TOPDIR/tmp
rm mplapack-${VERSION}.tar.xz
cp $TOPDIR/mplapack-${VERSION}.tar.xz .
tar xvf mplapack-${VERSION}.tar.xz
cd mplapack-${VERSION}
export CC='ccache gcc'
export CXX='ccache g++'
export FC='ccache gfortran'
./configure --prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-_Float128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-_Float64x=yes --enable-test=yes --enable-cuda
make -j`getconf _NPROCESSORS_ONLN`
make install
