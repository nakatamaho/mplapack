
MPLAPACKVER=0.8.0
rm -rf build
mkdir build
cd build

mkdir gcc
cd gcc
md5sum ~/mplapack-$MPLAPACKVER.tar.gz
tar xvfz ~/mplapack-$MPLAPACKVER.tar.gz
cd mplapack-$MPLAPACKVER
./configure --enable-optimization=yes --enable-debug=yes --enable-mpfr=yes --enable-qd=yes --enable-double=yes --enable-__float128=yes --prefix=/home/docker/MPLAPACK.GCC
/usr/bin/time make -j20
/usr/bin/time make install
cd ..
cd ..

mkdir icc
cd icc
md5sum ~/mplapack-$MPLAPACKVER.tar.gz
tar xvfz ~/mplapack-$MPLAPACKVER.tar.gz
cd mplapack-$MPLAPACKVER
source /opt/intel/composer_xe_2013/bin/compilervars.sh intel64
export CC="icc"
export CXX="icpc"
export FC="ifort"
./configure --enable-optimization=yes --enable-debug=yes --enable-mpfr=yes --enable-qd=yes --enable-double=yes --enable-__float128=yes --prefix=/home/docker/MPLAPACK.ICC --with-external-blas="-mkl" --with-external-lapack="-mkl"
/usr/bin/time make -j20
/usr/bin/time make install
cd ..
cd ..
