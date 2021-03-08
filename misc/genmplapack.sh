rm mplapack_gmp.h
rm mplapack_qd.h
rm mplapack_dd.h

cp mplapack_generic.h mplapack_mpfr.h 
sed -i.bak "s/MPLAPACK_GENERIC_H/MPLAPACK_MPFR_H/g" mplapack_mpfr.h
sed -i.bak "s/INTEGER/mplapackint/g" mplapack_mpfr.h
sed -i.bak "s/LOGICAL/mplapacklogical/g" mplapack_mpfr.h
sed -i.bak "s/REAL/mpreal/g" mplapack_mpfr.h
sed -i.bak "s/COMPLEX/mpcomplex/g" mplapack_mpfr.h
sed -i.bak "/_gmp/d" mplapack_mpfr.h
sed -i.bak "/_qd/d" mplapack_mpfr.h
sed -i.bak "/_dd/d" mplapack_mpfr.h
sed -i.bak "/_double/d" mplapack_mpfr.h

cp mplapack_generic.h mplapack_gmp.h 
sed -i.bak "s/MPLAPACK_GENERIC_H/MPLAPACK_GMP_H/g" mplapack_gmp.h
sed -i.bak "s/INTEGER/mplapackint/g" mplapack_gmp.h
sed -i.bak "s/LOGICAL/mplapacklogical/g" mplapack_gmp.h
sed -i.bak "s/REAL/mpf_class/g" mplapack_gmp.h
sed -i.bak "s/COMPLEX/mpc_class/g" mplapack_gmp.h
sed -i.bak "/_mpfr/d" mplapack_gmp.h
sed -i.bak "/_qd/d" mplapack_gmp.h
sed -i.bak "/_dd/d" mplapack_gmp.h
sed -i.bak "/_double/d" mplapack_gmp.h

cp mplapack_generic.h mplapack_qd.h 
sed -i.bak "s/MPLAPACK_GENERIC_H/MPLAPACK_QD_H/g" mplapack_qd.h
sed -i.bak "s/INTEGER/mplapackint/g" mplapack_qd.h
sed -i.bak "s/LOGICAL/mplapacklogical/g" mplapack_qd.h
sed -i.bak "s/REAL/qd_real/g" mplapack_qd.h
sed -i.bak "s/COMPLEX/qd_complex/g" mplapack_qd.h
sed -i.bak "/_gmp/d" mplapack_qd.h
sed -i.bak "/_mpfr/d" mplapack_qd.h
sed -i.bak "/_dd/d" mplapack_qd.h
sed -i.bak "/_double/d" mplapack_qd.h

cp mplapack_generic.h mplapack_dd.h 
sed -i.bak "s/MPLAPACK_GENERIC_H/MPLAPACK_DD_H/g" mplapack_dd.h
sed -i.bak "s/INTEGER/mplapackint/g" mplapack_dd.h
sed -i.bak "s/LOGICAL/mplapacklogical/g" mplapack_dd.h
sed -i.bak "s/REAL/dd_real/g" mplapack_dd.h
sed -i.bak "s/COMPLEX/dd_complex/g" mplapack_dd.h
sed -i.bak "/_gmp/d" mplapack_dd.h
sed -i.bak "/_mpfr/d" mplapack_dd.h
sed -i.bak "/_qd/d" mplapack_dd.h
sed -i.bak "/_double/d" mplapack_dd.h

cp mplapack_generic.h mplapack_double.h 
sed -i.bak "s/MPLAPACK_GENERIC_H/MPLAPACK_DOUBLE_H/g" mplapack_double.h
sed -i.bak "s/INTEGER/mplapackint/g" mplapack_double.h
sed -i.bak "s/LOGICAL/mplapacklogical/g" mplapack_double.h
sed -i.bak "s/REAL/double/g" mplapack_double.h
sed -i.bak "s/COMPLEX/std::complex<double>/g" mplapack_double.h
sed -i.bak "/_gmp/d" mplapack_double.h
sed -i.bak "/_mpfr/d" mplapack_double.h
sed -i.bak "/_qd/d" mplapack_double.h
sed -i.bak "/_dd/d" mplapack_double.h
