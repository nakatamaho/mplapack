cp Makefile.am ../qd
cp Makefile.am ../dd
sed -i.bak "s/gmp/qd/g" ../qd/Makefile.am
sed -i.bak "s/gmp/dd/g" ../dd/Makefile.am
sed -i.bak "s/-lqdxx //g" ../qd/Makefile.am
sed -i.bak "s/-lddxx //g" ../dd/Makefile.am
sed -i.bak "s/-ldd/-lqd/g" ../dd/Makefile.am
sed -i.bak "s/___MPLAPACK_BUILD_WITH_GMP___/___MPLAPACK_BUILD_WITH_QD___/g" ../qd/Makefile.am
sed -i.bak "s/___MPLAPACK_BUILD_WITH_GMP___/___MPLAPACK_BUILD_WITH_DD___/g" ../dd/Makefile.am

