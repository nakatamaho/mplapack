rm mlapack_gmp.h
rm mlapack_qd.h
rm mlapack_dd.h

cp mlapack_generic.h mlapack_mpfr.h 
sed -i.bak "s/MLAPACK_GENERIC_H/MLAPACK_MPFR_H/g" mlapack_mpfr.h
sed -i.bak "s/INTEGER/mpackint/g" mlapack_mpfr.h
sed -i.bak "s/LOGICAL/mpacklogical/g" mlapack_mpfr.h
sed -i.bak "s/REAL/mpreal/g" mlapack_mpfr.h
sed -i.bak "s/COMPLEX/mpcomplex/g" mlapack_mpfr.h
sed -i.bak "/_gmp/d" mlapack_mpfr.h
sed -i.bak "/_qd/d" mlapack_mpfr.h
sed -i.bak "/_dd/d" mlapack_mpfr.h
sed -i.bak "/_double/d" mlapack_mpfr.h

cp mlapack_generic.h mlapack_gmp.h 
sed -i.bak "s/MLAPACK_GENERIC_H/MLAPACK_GMP_H/g" mlapack_gmp.h
sed -i.bak "s/INTEGER/mpackint/g" mlapack_gmp.h
sed -i.bak "s/LOGICAL/mpacklogical/g" mlapack_gmp.h
sed -i.bak "s/REAL/mpf_class/g" mlapack_gmp.h
sed -i.bak "s/COMPLEX/mpc_class/g" mlapack_gmp.h
sed -i.bak "/_mpfr/d" mlapack_gmp.h
sed -i.bak "/_qd/d" mlapack_gmp.h
sed -i.bak "/_dd/d" mlapack_gmp.h
sed -i.bak "/_double/d" mlapack_gmp.h

cp mlapack_generic.h mlapack_qd.h 
sed -i.bak "s/MLAPACK_GENERIC_H/MLAPACK_QD_H/g" mlapack_qd.h
sed -i.bak "s/INTEGER/mpackint/g" mlapack_qd.h
sed -i.bak "s/LOGICAL/mpacklogical/g" mlapack_qd.h
sed -i.bak "s/REAL/qd_real/g" mlapack_qd.h
sed -i.bak "s/COMPLEX/qd_complex/g" mlapack_qd.h
sed -i.bak "/_gmp/d" mlapack_qd.h
sed -i.bak "/_mpfr/d" mlapack_qd.h
sed -i.bak "/_dd/d" mlapack_qd.h
sed -i.bak "/_double/d" mlapack_qd.h

cp mlapack_generic.h mlapack_dd.h 
sed -i.bak "s/MLAPACK_GENERIC_H/MLAPACK_DD_H/g" mlapack_dd.h
sed -i.bak "s/INTEGER/mpackint/g" mlapack_dd.h
sed -i.bak "s/LOGICAL/mpacklogical/g" mlapack_dd.h
sed -i.bak "s/REAL/dd_real/g" mlapack_dd.h
sed -i.bak "s/COMPLEX/dd_complex/g" mlapack_dd.h
sed -i.bak "/_gmp/d" mlapack_dd.h
sed -i.bak "/_mpfr/d" mlapack_dd.h
sed -i.bak "/_qd/d" mlapack_dd.h
sed -i.bak "/_double/d" mlapack_dd.h

cp mlapack_generic.h mlapack_double.h 
sed -i.bak "s/MLAPACK_GENERIC_H/MLAPACK_DOUBLE_H/g" mlapack_double.h
sed -i.bak "s/INTEGER/mpackint/g" mlapack_double.h
sed -i.bak "s/LOGICAL/mpacklogical/g" mlapack_double.h
sed -i.bak "s/REAL/double/g" mlapack_double.h
sed -i.bak "s/COMPLEX/std::complex<double>/g" mlapack_double.h
sed -i.bak "/_gmp/d" mlapack_double.h
sed -i.bak "/_mpfr/d" mlapack_double.h
sed -i.bak "/_qd/d" mlapack_double.h
sed -i.bak "/_dd/d" mlapack_double.h
