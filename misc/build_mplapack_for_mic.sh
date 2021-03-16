#!/bin/bash
# a sample build script for building mplapack
# on Intel Xeon Phi
# how to use
# after doing cd mplapack/misc ; bash prepare_mplapacklibs_for_mic.sh,
# then 
# $ bash build_mplapack_for_mic.sh

source /opt/intel/composer_xe_2013/bin/compilervars.sh intel64

MICLIBS_HOME="`pwd`/MIC"
HOSTDUMMYLIBS_HOME="`pwd`/HOST"
MPLAPACKINSTALL_HOME="`pwd`/MPLAPACK_MIC"

QDVER=2.3.7
GMPVER=5.0.5
MPFRVER=3.1.1
MPCVER=1.0.1
MPLAPACKVER=0.8.0

export CC="icc"
export CXX="icpc"
export FC="ifort"
export CFLAGS="-DMMIC -I$MICLIBS_HOME/include"
export CXXFLAGS="-DMMIC -I$MICLIBS_HOME/include"
export FCFLAGS="-DMMIC -I$MICLIBS_HOME/include"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MICLIBS_HOME/lib:$HOSTDUMMYLIBS_HOME/lib

rm -rf mplapack-$MPLAPACKVER
tar xfz ../../mplapack-$MPLAPACKVER.tar.gz
cd mplapack-$MPLAPACKVER

perl -p -i -e 's/\@IS_MIC_TRUE\@//g' external/lapack/Makefile.in

LDFLAGS="-L$HOSTDUMMYLIBS_HOME/lib" ./configure --with-system-gmp --with-system-mpfr --with-system-mpc --with-system-qd --enable-debug=yes --enable-optimization=yes --enable-mpfr=yes --enable-_Float128=yes --with-gccname=/usr/linux-k1om-4.7/bin/x86_64-k1om-linux-gcc --enable-qd=yes --enable-double=yes --prefix=$MPLAPACKINSTALL_HOME --with-external-blas="-mkl" --with-external-lapack="-mkl"

files=$(find ./* -name Makefile)
perl -p -i -e 's/-DMMIC/-mmic/g' $files
perl -p -i -e 's/HOST/MIC/g' $files
perl -p -i -e 's/intel64/mic/g' $files
perl -p -i -e 's/\/lib64/\/usr\/linux-k1om-4.7\/linux-k1om\/lib64\//g' $files
perl -p -i -e 's/FORTRAN = \$\(FC\)/FORTRAN = \$\(FC\) -mmic/g' $files
perl -p -i -e 's/-xHost//g' $files
files2=$(find ./* -name libtool)
perl -p -i -e 's/intel64/mic/g' $files2

/usr/bin/time make -j20

