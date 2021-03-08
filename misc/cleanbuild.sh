git clone https://github.com/nakatamaho/mplapack.git
cd mplapack
libtoolize -c
aclocal
autoheader
automake --add-missing
autoreconf
autoconf
automake
./configure --enable-optimization=yes --enable-debug=yes --enable-mpfr=yes --enable-qd=yes --enable-double=yes --enable-__float128=yes --prefix=/home/docker/MPLAPACK.GCC
/usr/bin/time make -j20
/usr/bin/time make install
