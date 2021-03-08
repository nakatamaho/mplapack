#!/bin/bash
# a sample build script for preparing libraires required
# on Intel Xeon Phi
# how to use
# cd mplapack/misc ; bash prepare_mplapacklibs_for_mic.sh

source /opt/intel/composer_xe_2013/bin/compilervars.sh intel64
MICLIBS_HOME="`pwd`/MIC"
HOSTDUMMYLIBS_HOME="`pwd`/HOST"

QDVER=2.3.7
GMPVER=5.0.5
MPFRVER=3.1.1
MPCVER=1.0.1

export CC="icc"
export CXX="icpc"
export FC="ifort"
export CFLAGS="-DMMIC -I$MICLIBS_HOME/include"
export CXXFLAGS="-DMMIC -I$MICLIBS_HOME/include"
export FFLAGS="-DMMIC -I$MICLIBS_HOME/include"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MICLIBS_HOME/lib:$HOSTDUMMYLIBS_HOME/lib

rm -rf $MICLIBS_HOME
rm -rf $HOSTDUMMYLIBS_HOME
rm -rf tmp
mkdir tmp
cd tmp

mkdir mic
mkdir host

cd mic
rm -rf qd-$QDVER
tar xfz ../../../external/qd/download/qd-$QDVER.tar.gz
cd qd-$QDVER
LDFLAGS="-L$MICLIBS_HOME/lib" ./configure --prefix=$MICLIBS_HOME/ # --enable-fma
files=$(find ./* -name Makefile)
perl -p -i -e 's/-DMMIC/-mmic/g' $files
make ; make install
cd ..
cd ..

cd host
rm -rf qd-$QDVER
tar xfz ../../../external/qd/download/qd-$QDVER.tar.gz
cd qd-$QDVER
LDFLAGS="-L$HOSTDUMMYLIBS_HOME/lib" ./configure --prefix=$HOSTDUMMYLIBS_HOME/
make ; make install
cd ..
cd ..

cd mic
rm -rf gmp-$GMPVER
tar xfj ../../../external/gmp/download/gmp-$GMPVER.tar.bz2
cd gmp-$GMPVER
rm mpn/x86_64/divrem_1.asm
rm mpn/x86_64/mod_1_4.asm
rm mpn/x86_64/sqr_basecase.asm
rm mpn/x86_64/*/popcount.asm
rm mpn/x86_64/gcd_1.asm
LDFLAGS="-L$MICLIBS_HOME/lib" ./configure --prefix=$MICLIBS_HOME/ --enable-cxx
files=$(find ./* -name Makefile)
perl -p -i -e 's/-DMMIC/-mmic/g' $files
make ; make install
cd ..
cd ..

cd host
rm -rf gmp-$GMPVER
tar xfj ../../../external/gmp/download/gmp-$GMPVER.tar.bz2
cd gmp-$GMPVER
LDFLAGS="-L$HOSTDUMMYLIBS_HOME/lib" ./configure --prefix=$HOSTDUMMYLIBS_HOME/ --enable-cxx
make ; make install
cd ..
cd ..

cd mic
rm -rf mpfr-$MPFRVER
tar xfj ../../../external/mpfr/download/mpfr-$MPFRVER.tar.bz2
cd mpfr-$MPFRVER
LDFLAGS="-L$MICLIBS_HOME/lib" ./configure --prefix=$MICLIBS_HOME/ --with-gmp=$HOSTDUMMYLIBS_HOME/
files=$(find ./* -name Makefile)
perl -p -i -e 's/-DMMIC/-mmic/g' $files
Qmake; make install
cd ..
cd ..

cd host
rm -rf mpfr-$MPFRVER
tar xfj ../../../external/mpfr/download/mpfr-$MPFRVER.tar.bz2
cd mpfr-$MPFRVER
LDFLAGS="-L$HOSTDUMMYLIBS_HOME/lib" ./configure --prefix=$HOSTDUMMYLIBS_HOME/ --with-gmp=$HOSTDUMMYLIBS_HOME/
make; make install
cd ..
cd ..

cd mic
rm -rf mpc-$MPCVER
tar xfz ../../../external/mpc/download/mpc-$MPCVER.tar.gz
cd mpc-$MPCVER
LDFLAGS="-L$MICLIBS_HOME/lib" ./configure --prefix=$MICLIBS_HOME/ --with-gmp=$HOSTDUMMYLIBS_HOME/ --with-mpfr=$HOSTDUMMYLIBS_HOME/
files=$(find ./* -name Makefile)
perl -p -i -e 's/-DMMIC/-mmic/g' $files
perl -p -i -e 's/HOST/MIC/g' $files
make ; make install
cd ..
cd ..

cd host
rm -rf mpc-$MPCVER
tar xfz ../../../external/mpc/download/mpc-$MPCVER.tar.gz
cd mpc-$MPCVER
LDFLAGS="-L$HOSTDUMMYLIBS_HOME/lib" ./configure --prefix=$HOSTDUMMYLIBS_HOME/ --with-gmp=$HOSTDUMMYLIBS_HOME/ --with-mpfr=$HOSTDUMMYLIBS_HOME/
make ; make install
cd ..
cd ..



