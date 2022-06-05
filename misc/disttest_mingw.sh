TOPDIR=/home/docker/mplapack
VERSION=2.0.0
####################
cd $TOPDIR
pushd mplapack/test/compare ; bash gen.Makefile.am.sh ; popd
aclocal ; autoconf ; automake --add-missing
autoreconf --force --install
./configure --prefix=$HOME/MPLAPACK --host=x86_64-w64-mingw32 --target=x86_64-w64-mingw32 --enable-gmp=yes --enable-mpfr=yes --enable-_Float128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-_Float64x=yes --enable-test=yes
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
#
CXX="ccache x86_64-w64-mingw32-g++" ; export CXX
CC="ccache x86_64-w64-mingw32-gcc" ; export CC
FC="x86_64-w64-mingw32-gfortran"; export FC
F77="x86_64-w64-mingw32-gfortran"; export F77
NM="x86_64-w64-mingw32-nm" ; export NM
RANLIB="x86_64-w64-mingw32-ranlib" ; export RANLIB
AR="x86_64-w64-mingw32-ar" ; export AR
export WINEPATH="/usr/x86_64-w64-mingw32/lib/;/usr/lib/gcc/x86_64-w64-mingw32/9.3-win32/;/usr/lib/gcc/x86_64-w64-mingw32/9.3-posix"
export WINEDEBUG=-all
#
rm -rf ~/.wine
mkdir -p ~/.wine

#
./configure --prefix=$HOME/MPLAPACK --host=x86_64-w64-mingw32 --target=x86_64-w64-mingw32 --enable-gmp=yes --enable-mpfr=yes --enable-_Float128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-_Float64x=yes --enable-test=yes
make -j`getconf _NPROCESSORS_ONLN`
make install
