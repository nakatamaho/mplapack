TOPDIR=/home/docker/mplapack
VERSION=2.0.0

cd $TOPDIR
pushd mplapack/test/compare ; bash gen.Makefile.am.sh ; popd
aclocal ; autoconf ; automake --add-missing
autoreconf --force --install
make dist
rm -rf $TOPTDIR/tmp
cd $TOPTDIR/tmp
cp $TOPDIR/mplapack-${VERSION}.tar.gz .
tar xvfz mplapack-${VERSION}.tar.gz
./configure --prefix=$HOME/MPLAPACK --enable-gmp=yes --enable-mpfr=yes --enable-_Float128=yes --enable-qd=yes --enable-dd=yes --enable-double=yes --enable-_Float64x=yes --enable-test=yes

