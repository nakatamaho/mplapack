#!/bin/bash

NAME="Rpotrf"
cp "$NAME"_mpfr.cpp "$NAME"_gmp.cpp
cp "$NAME"_mpfr.cpp "$NAME"_dd.cpp
cp "$NAME"_mpfr.cpp "$NAME"_qd.cpp
cp "$NAME"_mpfr.cpp "$NAME"_double.cpp
cp "$NAME"_mpfr.cpp "$NAME"__Float128.cpp

sed -ibak 's/mpreal/dd_real/g' "$NAME"_dd.cpp
sed -ibak 's/MPFR/DD/g' "$NAME"_dd.cpp
sed -ibak 's/mpfr/dd/g' "$NAME"_dd.cpp
sed -ibak 's/___MPLAPACK_BUILD_WITH_MPFR___/___MPLAPACK_BUILD_WITH_DD___/g' "$NAME"_dd.cpp

sed -ibak 's/MPFR/GMP/g' "$NAME"_gmp.cpp
sed -ibak 's/mpfr/gmp/g' "$NAME"_gmp.cpp
sed -ibak 's/mpreal/mpf_class/g' "$NAME"_gmp.cpp
sed -ibak 's/___MPLAPACK_BUILD_WITH_MPFR___/___MPLAPACK_BUILD_WITH_GMP___/g' "$NAME"_gmp.cpp

sed -ibak 's/mpreal/qd_real/g' "$NAME"_qd.cpp
sed -ibak 's/MPFR/QD/g' "$NAME"_qd.cpp
sed -ibak 's/mpfr/qd/g' "$NAME"_qd.cpp
sed -ibak 's/___MPLAPACK_BUILD_WITH_MPFR___/___MPLAPACK_BUILD_WITH_QD___/g' "$NAME"_qd.cpp

sed -ibak 's/mpreal/_Float128/g' "$NAME"__Float128.cpp
sed -ibak 's/MPFR/_FLOAT128/g' "$NAME"__Float128.cpp
sed -ibak 's/mpfr/_Float128/g' "$NAME"__Float128.cpp
sed -ibak 's/___MPLAPACK_BUILD_WITH_MPFR___/___MPLAPACK_BUILD_WITH__FLOAT128___/g' "$NAME"__Float128.cpp

sed -ibak 's/mpreal/double/g' "$NAME"_double.cpp
sed -ibak 's/MPFR/DOUBLE/g' "$NAME"_double.cpp
sed -ibak 's/mpfr/double/g' "$NAME"_double.cpp
sed -ibak 's/___MPLAPACK_BUILD_WITH_MPFR___/___MPLAPACK_BUILD_WITH_DOUBLE___/g' "$NAME"_double.cpp
